/*
    LensSerious — build the database from the upstream XML.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file import_lensfun_xml.c
 *
 * @brief The periodic rebuild: upstream XML in, one SQLite file out.
 *
 * @details This tool reads the XML **through liblensfun**, and that is a deliberate
 * choice worth defending, because at first glance it looks like the one dependency this
 * project exists to remove.
 *
 * It is not the same dependency. What LensSerious replaces is lensfun's RUNTIME: the
 * ~90 ms XML parse in every session, the process-wide database object, the mutex every
 * lookup serialises on, the callback machinery that can only produce a displacement map
 * on one CPU thread. None of that is here. This is an offline tool, run when upstream
 * publishes new calibrations, on a maintainer's machine or in CI. Its dependencies cost
 * a rendered frame nothing.
 *
 * What using liblensfun buys is correctness that would otherwise have to be re-earned.
 * The XML is not a straightforward document: a lens's focal and aperture RANGES are
 * usually absent and inferred from its model string by lfLens::GuessParameters(), and
 * those ranges feed the vignetting interpolation's distance metric -- get them wrong and
 * every vignetting lookup on that lens is silently off. Reimplementing that heuristic
 * blind, against 1562 lenses, is a large surface of quiet errors. Reading through the
 * library that already does it means the database contains, by construction, exactly what
 * lensfun itself would have computed. tests/parity_db.c then asserts that end to end.
 *
 * The exit is open: nothing in the schema or in lensserious_db.c knows liblensfun exists.
 * When a native XML reader is written it replaces this file alone, and the parity test is
 * already there to hold it to the same standard.
 *
 * @section rebuild The rebuild contract
 *
 * The database is written to a temporary file beside the target and renamed into place.
 * Readers open with `immutable=1` and take no locks (see lensserious_db.h); that promise
 * only holds if a live file is never edited. rename(2) is atomic within a filesystem, so
 * a reader gets one whole version or the other, and a reader that already had the old
 * file open keeps reading the inode it opened.
 */

#include "lensserious.h"
#include "lensserious_db.h"

#include <lensfun.h>
#include <sqlite3.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static const char *SCHEMA_PATH_HINT = "src/schema.sql";

/* ------------------------------------------------------------------------- */

static void die(const char *what, sqlite3 *db)
{
  fprintf(stderr, "import: %s%s%s\n", what,
          db ? ": " : "", db ? sqlite3_errmsg(db) : "");
  exit(1);
}

static void exec_or_die(sqlite3 *db, const char *sql)
{
  char *err = NULL;
  if(sqlite3_exec(db, sql, NULL, NULL, &err) != SQLITE_OK)
  {
    fprintf(stderr, "import: %s\n", err ? err : "(no message)");
    sqlite3_free(err);
    exit(1);
  }
}

static sqlite3_stmt *prep(sqlite3 *db, const char *sql)
{
  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db, sql, -1, &st, NULL) != SQLITE_OK) die(sql, db);
  return st;
}

static void step_reset(sqlite3 *db, sqlite3_stmt *st)
{
  if(sqlite3_step(st) != SQLITE_DONE) die("insert", db);
  sqlite3_reset(st);
  sqlite3_clear_bindings(st);
}

/* ------------------------------------------------------------------------- */

/** @brief Walk an lfMLstr: the default string, then (lang, translation) pairs. */
static void for_each_mlstr(const lfMLstr s,
                           void (*fn)(const char *lang, const char *value, void *ud),
                           void *ud)
{
  if(!s) return;
  const char *p = s;
  fn(NULL, p, ud);                    /* the default, untranslated spelling */
  p += strlen(p) + 1;
  while(*p)
  {
    const char *lang = p;
    p += strlen(p) + 1;
    if(!*p && !*(p)) { /* malformed tail: stop rather than run off the end */ }
    const char *value = p;
    fn(lang, value, ud);
    p += strlen(p) + 1;
  }
}

typedef struct
{
  sqlite3 *db;
  sqlite3_stmt *st;
  sqlite3_stmt *tok;      /* NULL for cameras: only lenses are fuzzy-matched */
  long long owner_id;
  const char *kind;
} name_ctx_t;

static void insert_name(const char *lang, const char *value, void *ud)
{
  name_ctx_t *c = (name_ctx_t *)ud;
  if(!value || !*value) return;

  char norm[512];
  ls_db_normalize(value, norm, sizeof(norm));

  sqlite3_bind_int64(c->st, 1, c->owner_id);
  sqlite3_bind_text(c->st, 2, c->kind, -1, SQLITE_STATIC);
  if(lang) sqlite3_bind_text(c->st, 3, lang, -1, SQLITE_TRANSIENT);
  else sqlite3_bind_null(c->st, 3);
  sqlite3_bind_text(c->st, 4, value, -1, SQLITE_TRANSIENT);
  sqlite3_bind_text(c->st, 5, norm, -1, SQLITE_TRANSIENT);
  /* Only lens names carry a token digest -- camera lookups are exact and never score -- and
   * camera_name has no such column, so this statement has five parameters, not six. */
  if(c->tok)
  {
    unsigned char digest[2 + 32 * 5];
    const size_t dn = ls_db_token_digest(norm, digest, sizeof(digest));
    sqlite3_bind_blob(c->st, 6, digest, (int)dn, SQLITE_TRANSIENT);
  }
  step_reset(c->db, c->st);

  if(!c->tok) return;
  /* Index this name's tokens, so a lookup can ask which of ITS tokens is rarest instead of
   * scoring the whole catalogue. Tokenised here, offline, by the same function the matcher
   * uses on the query. */
  char tokens[32][48];
  const int n = ls_db_tokenize(norm, &tokens[0][0], 32, 48);
  for(int i = 0; i < n; i++)
  {
    sqlite3_bind_int64(c->tok, 1, c->owner_id);
    sqlite3_bind_text(c->tok, 2, c->kind, -1, SQLITE_STATIC);
    sqlite3_bind_text(c->tok, 3, tokens[i], -1, SQLITE_TRANSIENT);
    step_reset(c->db, c->tok);
  }
}

/* ------------------------------------------------------------------------- */
/* Model vocabularies. The XML's names become the evaluator's integers here,   */
/* once, offline -- never at render time.                                     */
/* ------------------------------------------------------------------------- */

static int dist_model(enum lfDistortionModel m)
{
  switch(m)
  {
    case LF_DIST_MODEL_POLY3:  return LS_DIST_POLY3;
    case LF_DIST_MODEL_POLY5:  return LS_DIST_POLY5;
    case LF_DIST_MODEL_PTLENS: return LS_DIST_PTLENS;
    default: return LS_DIST_NONE;
  }
}

static int tca_model(enum lfTCAModel m)
{
  switch(m)
  {
    case LF_TCA_MODEL_LINEAR: return LS_TCA_LINEAR;
    case LF_TCA_MODEL_POLY3:  return LS_TCA_POLY3;
    default: return LS_TCA_NONE;
  }
}

static int lens_type(enum lfLensType t)
{
  switch(t)
  {
    case LF_RECTILINEAR:             return LS_LENS_RECTILINEAR;
    case LF_FISHEYE:                 return LS_LENS_FISHEYE;
    case LF_PANORAMIC:               return LS_LENS_PANORAMIC;
    case LF_EQUIRECTANGULAR:         return LS_LENS_EQUIRECTANGULAR;
    case LF_FISHEYE_ORTHOGRAPHIC:    return LS_LENS_FISHEYE_ORTHOGRAPHIC;
    case LF_FISHEYE_STEREOGRAPHIC:   return LS_LENS_FISHEYE_STEREOGRAPHIC;
    case LF_FISHEYE_EQUISOLID:       return LS_LENS_FISHEYE_EQUISOLID;
    case LF_FISHEYE_THOBY:           return LS_LENS_FISHEYE_THOBY;
    default:                         return LS_LENS_UNKNOWN;
  }
}

/* ------------------------------------------------------------------------- */

static long long mount_id_for(sqlite3 *db, sqlite3_stmt *sel, sqlite3_stmt *ins,
                              const char *name)
{
  if(!name || !*name) return 0;

  sqlite3_bind_text(sel, 1, name, -1, SQLITE_TRANSIENT);
  long long id = 0;
  if(sqlite3_step(sel) == SQLITE_ROW) id = sqlite3_column_int64(sel, 0);
  sqlite3_reset(sel);
  sqlite3_clear_bindings(sel);
  if(id) return id;

  sqlite3_bind_text(ins, 1, name, -1, SQLITE_TRANSIENT);
  step_reset(db, ins);
  return sqlite3_last_insert_rowid(db);
}

static char *read_file(const char *path)
{
  FILE *f = fopen(path, "rb");
  if(!f) return NULL;
  fseek(f, 0, SEEK_END);
  const long n = ftell(f);
  fseek(f, 0, SEEK_SET);
  char *buf = (char *)malloc((size_t)n + 1);
  if(!buf || fread(buf, 1, (size_t)n, f) != (size_t)n)
  {
    free(buf);
    fclose(f);
    return NULL;
  }
  buf[n] = '\0';
  fclose(f);
  return buf;
}

int main(int argc, char **argv)
{
  if(argc < 3)
  {
    fprintf(stderr,
            "usage: %s <schema.sql> <out.db> [lensfun-xml-dir]\n"
            "\n"
            "Reads the lensfun database (the system one, or the directory given) and\n"
            "writes it as a LensSerious SQLite database. The output is written to a\n"
            "temporary file and renamed into place, so readers using immutable=1 never\n"
            "see a partial file -- see lensserious_db.h.\n",
            argv[0]);
    return 2;
  }
  const char *schema_path = argv[1];
  const char *out_path = argv[2];
  const char *xml_dir = (argc > 3) ? argv[3] : NULL;

  char *schema = read_file(schema_path);
  if(!schema)
  {
    fprintf(stderr, "import: cannot read schema `%s' (expected %s)\n",
            schema_path, SCHEMA_PATH_HINT);
    return 1;
  }

  /* --- read upstream ---------------------------------------------------- */

  lfDatabase *ldb = lf_db_new();
  if(!ldb) die("lf_db_new", NULL);

  if(xml_dir)
  {
    if(lf_db_load_directory(ldb, xml_dir) != 0)
    {
      fprintf(stderr, "import: no database loaded from `%s'\n", xml_dir);
      return 1;
    }
  }
  else if(lf_db_load(ldb) != LF_NO_ERROR)
  {
    fprintf(stderr, "import: no system lensfun database found\n");
    return 1;
  }

  const lfCamera *const *cameras = lf_db_get_cameras(ldb);
  const lfLens *const *lenses = lf_db_get_lenses(ldb);
  const lfMount *const *mounts = lf_db_get_mounts(ldb);

  /* --- write ------------------------------------------------------------ */

  char tmp_path[4096];
  snprintf(tmp_path, sizeof(tmp_path), "%s.tmp-%ld", out_path, (long)getpid());
  unlink(tmp_path);

  sqlite3 *db = NULL;
  if(sqlite3_open_v2(tmp_path, &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL) != SQLITE_OK)
    die("open output", db);

  /* Nothing here needs crash-safety: on failure the temporary file is discarded and the
   * live database is untouched, so durability guarantees are pure cost. */
  exec_or_die(db, "PRAGMA journal_mode = OFF");
  exec_or_die(db, "PRAGMA synchronous = OFF");
  exec_or_die(db, schema);
  free(schema);
  exec_or_die(db, "BEGIN");

  sqlite3_stmt *ins_meta = prep(db, "INSERT INTO meta(key, value) VALUES (?1, ?2)");
  sqlite3_stmt *sel_mount = prep(db, "SELECT id FROM mount WHERE name = ?1");
  sqlite3_stmt *ins_mount = prep(db, "INSERT INTO mount(name) VALUES (?1)");
  sqlite3_stmt *ins_compat = prep(db, "INSERT OR IGNORE INTO mount_compat(mount_id, compat_id) VALUES (?1, ?2)");
  sqlite3_stmt *ins_cam = prep(db, "INSERT INTO camera(maker, model, variant, mount_id, crop_factor)"
                                   " VALUES (?1, ?2, ?3, ?4, ?5)");
  sqlite3_stmt *ins_cam_name = prep(db, "INSERT INTO camera_name(camera_id, kind, lang, value, norm)"
                                        " VALUES (?1, ?2, ?3, ?4, ?5)");
  sqlite3_stmt *ins_lens = prep(db, "INSERT INTO lens(maker, model, type, crop_factor, aspect_ratio,"
                                    " center_x, center_y, min_focal, max_focal, min_aperture, max_aperture)"
                                    " VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10,?11)");
  sqlite3_stmt *ins_lens_name = prep(db, "INSERT INTO lens_name(lens_id, kind, lang, value, norm, tokens)"
                                         " VALUES (?1, ?2, ?3, ?4, ?5, ?6)");
  sqlite3_stmt *ins_token = prep(db, "INSERT INTO lens_token(lens_id, kind, token)"
                                     " VALUES (?1, ?2, ?3)");
  sqlite3_stmt *ins_lens_mount = prep(db, "INSERT OR IGNORE INTO lens_mount(lens_id, mount_id) VALUES (?1, ?2)");
  sqlite3_stmt *ins_dist = prep(db, "INSERT INTO calib_distortion(lens_id, ord, model, focal, t0, t1, t2)"
                                    " VALUES (?1,?2,?3,?4,?5,?6,?7)");
  sqlite3_stmt *ins_tca = prep(db, "INSERT INTO calib_tca(lens_id, ord, model, focal, t0,t1,t2,t3,t4,t5)"
                                   " VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10)");
  sqlite3_stmt *ins_vig = prep(db, "INSERT INTO calib_vignetting(lens_id, ord, model, focal, aperture,"
                                   " distance, t0, t1, t2) VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9)");
  sqlite3_stmt *ins_real_focal = prep(db, "INSERT INTO lens_real_focal(lens_id, focal, real_focal)"
                                          " VALUES (?1, ?2, ?3)");

  /* mounts, then their compatibility lists (which reference mounts by name) */
  int n_mounts = 0;
  for(const lfMount *const *m = mounts; m && *m; m++, n_mounts++)
    mount_id_for(db, sel_mount, ins_mount, lf_mlstr_get((*m)->Name));

  for(const lfMount *const *m = mounts; m && *m; m++)
  {
    const long long id = mount_id_for(db, sel_mount, ins_mount, lf_mlstr_get((*m)->Name));
    if(!id || !(*m)->Compat) continue;
    for(char **c = (*m)->Compat; *c; c++)
    {
      const long long cid = mount_id_for(db, sel_mount, ins_mount, *c);
      if(!cid) continue;
      sqlite3_bind_int64(ins_compat, 1, id);
      sqlite3_bind_int64(ins_compat, 2, cid);
      step_reset(db, ins_compat);
    }
  }

  int n_cameras = 0;
  for(const lfCamera *const *c = cameras; c && *c; c++, n_cameras++)
  {
    const lfCamera *cam = *c;
    const long long mid = mount_id_for(db, sel_mount, ins_mount, cam->Mount);

    sqlite3_bind_text(ins_cam, 1, lf_mlstr_get(cam->Maker), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(ins_cam, 2, lf_mlstr_get(cam->Model), -1, SQLITE_TRANSIENT);
    if(cam->Variant) sqlite3_bind_text(ins_cam, 3, lf_mlstr_get(cam->Variant), -1, SQLITE_TRANSIENT);
    else sqlite3_bind_null(ins_cam, 3);
    if(mid) sqlite3_bind_int64(ins_cam, 4, mid); else sqlite3_bind_null(ins_cam, 4);
    sqlite3_bind_double(ins_cam, 5, cam->CropFactor);
    step_reset(db, ins_cam);

    const long long id = sqlite3_last_insert_rowid(db);
    name_ctx_t ctx = { db, ins_cam_name, NULL, id, "maker" };
    for_each_mlstr(cam->Maker, insert_name, &ctx);
    ctx.kind = "model";
    for_each_mlstr(cam->Model, insert_name, &ctx);
  }

  int n_lenses = 0, n_dist = 0, n_tca = 0, n_vig = 0;
  for(const lfLens *const *l = lenses; l && *l; l++, n_lenses++)
  {
    const lfLens *lens = *l;

    sqlite3_bind_text(ins_lens, 1, lf_mlstr_get(lens->Maker), -1, SQLITE_TRANSIENT);
    sqlite3_bind_text(ins_lens, 2, lf_mlstr_get(lens->Model), -1, SQLITE_TRANSIENT);
    sqlite3_bind_int(ins_lens, 3, lens_type(lens->Type));
    sqlite3_bind_double(ins_lens, 4, lens->CropFactor);
    sqlite3_bind_double(ins_lens, 5, lens->AspectRatio);
    sqlite3_bind_double(ins_lens, 6, lens->CenterX);
    sqlite3_bind_double(ins_lens, 7, lens->CenterY);
    sqlite3_bind_double(ins_lens, 8, lens->MinFocal);
    sqlite3_bind_double(ins_lens, 9, lens->MaxFocal);
    sqlite3_bind_double(ins_lens, 10, lens->MinAperture);
    sqlite3_bind_double(ins_lens, 11, lens->MaxAperture);
    step_reset(db, ins_lens);

    const long long id = sqlite3_last_insert_rowid(db);

    name_ctx_t ctx = { db, ins_lens_name, ins_token, id, "maker" };
    for_each_mlstr(lens->Maker, insert_name, &ctx);
    ctx.kind = "model";
    for_each_mlstr(lens->Model, insert_name, &ctx);

    if(lens->Mounts)
      for(char **m = lens->Mounts; *m; m++)
      {
        const long long mid = mount_id_for(db, sel_mount, ins_mount, *m);
        if(!mid) continue;
        sqlite3_bind_int64(ins_lens_mount, 1, id);
        sqlite3_bind_int64(ins_lens_mount, 2, mid);
        step_reset(db, ins_lens_mount);
      }

    if(lens->CalibDistortion)
      for(lfLensCalibDistortion **c = lens->CalibDistortion; *c; c++, n_dist++)
      {
        sqlite3_bind_int64(ins_dist, 1, id);
        sqlite3_bind_int(ins_dist, 2, (int)(c - lens->CalibDistortion));
        sqlite3_bind_int(ins_dist, 3, dist_model((*c)->Model));
        sqlite3_bind_double(ins_dist, 4, (*c)->Focal);
        for(int t = 0; t < 3; t++) sqlite3_bind_double(ins_dist, 5 + t, (*c)->Terms[t]);
        step_reset(db, ins_dist);
      }

    if(lens->CalibTCA)
      for(lfLensCalibTCA **c = lens->CalibTCA; *c; c++, n_tca++)
      {
        sqlite3_bind_int64(ins_tca, 1, id);
        sqlite3_bind_int(ins_tca, 2, (int)(c - lens->CalibTCA));
        sqlite3_bind_int(ins_tca, 3, tca_model((*c)->Model));
        sqlite3_bind_double(ins_tca, 4, (*c)->Focal);
        for(int t = 0; t < 6; t++) sqlite3_bind_double(ins_tca, 5 + t, (*c)->Terms[t]);
        step_reset(db, ins_tca);
      }

    if(lens->CalibVignetting)
      for(lfLensCalibVignetting **c = lens->CalibVignetting; *c; c++, n_vig++)
      {
        sqlite3_bind_int64(ins_vig, 1, id);
        sqlite3_bind_int(ins_vig, 2, (int)(c - lens->CalibVignetting));
        sqlite3_bind_int(ins_vig, 3, LS_VIG_PA);
        sqlite3_bind_double(ins_vig, 4, (*c)->Focal);
        sqlite3_bind_double(ins_vig, 5, (*c)->Aperture);
        sqlite3_bind_double(ins_vig, 6, (*c)->Distance);
        for(int t = 0; t < 3; t++) sqlite3_bind_double(ins_vig, 7 + t, (*c)->Terms[t]);
        step_reset(db, ins_vig);
      }

    /* <real-focal-length>. Rare -- a few dozen lenses -- and load-bearing where it exists:
     * it is the focal the PROJECTION runs on, and on the Sigma 4.5mm circular fisheye it is
     * 0.47x the nominal one. Zero entries are dropped here rather than at read time because
     * upstream skips them when interpolating, so a stored zero is a row nothing can use. */
    if(lens->CalibRealFocal)
      for(lfLensCalibRealFocal **c = lens->CalibRealFocal; *c; c++)
      {
        if((*c)->RealFocal == 0.f) continue;
        sqlite3_bind_int64(ins_real_focal, 1, id);
        sqlite3_bind_double(ins_real_focal, 2, (*c)->Focal);
        sqlite3_bind_double(ins_real_focal, 3, (*c)->RealFocal);
        step_reset(db, ins_real_focal);
      }
  }

  /* provenance, so a database can always say where it came from */
  {
    char buf[256];
    const time_t now = time(NULL);
    /* gmtime_r is POSIX and absent on MinGW, where this tool now builds because a consumer
     * generates its database at build time. gmtime() returns a pointer into static storage,
     * which matters not at all here: this is a single-threaded offline importer writing one
     * timestamp. The three-way #ifdef the "thread-safe" variants would need is more risk
     * than the race it would prevent. */
    const struct tm *tm_utc = gmtime(&now);
    if(tm_utc)
      strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", tm_utc);
    else
      snprintf(buf, sizeof(buf), "unknown");

    char lf_version[64];
    snprintf(lf_version, sizeof(lf_version), "lensfun %d.%d.%d.%d",
             LF_VERSION_MAJOR, LF_VERSION_MINOR, LF_VERSION_MICRO, LF_VERSION_BUGFIX);

    const struct { const char *k; const char *v; } rows[] = {
      { "built_utc", buf },
      { "source", xml_dir ? xml_dir : "system lensfun database" },
      { "importer", lf_version },
    };
    for(size_t i = 0; i < sizeof(rows) / sizeof(rows[0]); i++)
    {
      sqlite3_bind_text(ins_meta, 1, rows[i].k, -1, SQLITE_STATIC);
      sqlite3_bind_text(ins_meta, 2, rows[i].v, -1, SQLITE_TRANSIENT);
      step_reset(db, ins_meta);
    }
  }

  /* Derive the token frequencies from what was just inserted, rather than counting in C:
   * one statement, and it cannot disagree with the table it summarises. */
  exec_or_die(db, "INSERT INTO token_df(kind, token, df)"
                  " SELECT kind, token, COUNT(DISTINCT lens_id) FROM lens_token"
                  " GROUP BY kind, token");

  exec_or_die(db, "COMMIT");
  exec_or_die(db, "ANALYZE");
  /* One contiguous file with no free pages: it is never written again, and every reader
   * mmaps it. */
  exec_or_die(db, "VACUUM");

  sqlite3_finalize(ins_meta);
  sqlite3_finalize(sel_mount);
  sqlite3_finalize(ins_mount);
  sqlite3_finalize(ins_compat);
  sqlite3_finalize(ins_cam);
  sqlite3_finalize(ins_cam_name);
  sqlite3_finalize(ins_lens);
  sqlite3_finalize(ins_lens_name);
  sqlite3_finalize(ins_token);
  sqlite3_finalize(ins_lens_mount);
  sqlite3_finalize(ins_dist);
  sqlite3_finalize(ins_real_focal);
  sqlite3_finalize(ins_tca);
  sqlite3_finalize(ins_vig);
  if(sqlite3_close(db) != SQLITE_OK) die("close output", db);

  lf_db_destroy(ldb);

  /* Atomic swap: see @ref rebuild. */
  /* POSIX rename() replaces an existing destination atomically; the Windows CRT's does
   * NOT -- it fails outright when the target exists, so a rebuild over yesterday's database
   * would leave the .tmp behind and report failure. Removing the target first gives up
   * atomicity on Windows only, which is the correct trade for a build-time tool: the file
   * is generated, not user data, and a failed build is more visible than a torn write.
   * Everywhere else the rename stays atomic, which is what lets a reader hold the database
   * open with immutable=1 while it is being replaced. */
#ifdef _WIN32
  unlink(out_path);
#endif
  if(rename(tmp_path, out_path) != 0)
  {
    fprintf(stderr, "import: cannot rename `%s' to `%s': %s\n",
            tmp_path, out_path, strerror(errno));
    unlink(tmp_path);
    return 1;
  }

  printf("import: %d mounts, %d cameras, %d lenses "
         "(%d distortion, %d tca, %d vignetting calibrations) -> %s\n",
         n_mounts, n_cameras, n_lenses, n_dist, n_tca, n_vig, out_path);
  return 0;
}
