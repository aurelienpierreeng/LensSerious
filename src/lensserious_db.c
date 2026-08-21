/*
    LensSerious — reading the lens database.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.

    See include/lensserious_db.h for the contract. The short version: the file is opened
    read-only and immutable, so SQLite takes no file locks and keeps no shared state; the
    handle is opened NOMUTEX, so SQLite takes no mutex either. Neither does this file.
    What makes that safe is that a database is replaced by rename and never edited.
*/

#include "lensserious_db.h"

#include <sqlite3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LS_DB_SCHEMA_VERSION 1

struct ls_db_t
{
  sqlite3 *sql;
  int schema_version;
  char error[256];
};

static void _db_err(ls_db_t *db, const char *what)
{
  if(!db) return;
  const char *detail = db->sql ? sqlite3_errmsg(db->sql) : "no connection";
  snprintf(db->error, sizeof(db->error), "%s: %s", what, detail);
}

/**
 * @brief Build the URI this library insists on, whatever the caller passed.
 *
 * @details `mode=ro` and `immutable=1` are not the caller's to choose: every promise in
 * the header rests on them. A caller's own parameters are kept, so a path on a network
 * share can still carry e.g. `nolock=1`.
 *
 * A bare path has to be percent-encoded before it can become a URI -- '?' and '#' in a
 * directory name would otherwise be read as the start of the query, and SQLite would open
 * some other file, or none.
 */
static char *_db_build_uri(const char *path)
{
  if(!path || !*path) return NULL;

  const int is_uri = strncmp(path, "file:", 5) == 0;
  const size_t len = strlen(path);
  /* Worst case: every byte percent-encoded, plus the scheme and our parameters. */
  char *uri = (char *)malloc(len * 3 + 64);
  if(!uri) return NULL;

  size_t w = 0;
  if(is_uri)
  {
    memcpy(uri, path, len);
    w = len;
  }
  else
  {
    memcpy(uri, "file:", 5);
    w = 5;
    for(size_t i = 0; i < len; i++)
    {
      const unsigned char c = (unsigned char)path[i];
      if(c == '?' || c == '#' || c == '%')
      {
        static const char hex[] = "0123456789ABCDEF";
        uri[w++] = '%';
        uri[w++] = hex[c >> 4];
        uri[w++] = hex[c & 0xF];
      }
      else
        uri[w++] = (char)c;
    }
  }

  const char *sep = strchr(uri, '?') ? "&" : "?";
  w += (size_t)snprintf(uri + w, 64, "%smode=ro&immutable=1", sep);
  uri[w] = '\0';
  return uri;
}

ls_db_t *ls_db_open(const char *path)
{
  char *uri = _db_build_uri(path);
  if(!uri) return NULL;

  ls_db_t *db = (ls_db_t *)calloc(1, sizeof(ls_db_t));
  if(!db)
  {
    free(uri);
    return NULL;
  }
  db->schema_version = -1;

  /* NOMUTEX: this handle belongs to one thread and SQLite need not defend it. READONLY
   * and the immutable URI parameter between them mean no locking at all. */
  const int rc = sqlite3_open_v2(uri, &db->sql,
                                 SQLITE_OPEN_READONLY | SQLITE_OPEN_NOMUTEX | SQLITE_OPEN_URI,
                                 NULL);
  free(uri);
  if(rc != SQLITE_OK)
  {
    ls_db_close(db);
    return NULL;
  }

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql, "PRAGMA user_version", -1, &st, NULL) == SQLITE_OK
     && sqlite3_step(st) == SQLITE_ROW)
    db->schema_version = sqlite3_column_int(st, 0);
  sqlite3_finalize(st);

  if(db->schema_version != LS_DB_SCHEMA_VERSION)
  {
    ls_db_close(db);
    return NULL;
  }
  return db;
}

void ls_db_close(ls_db_t *db)
{
  if(!db) return;
  if(db->sql) sqlite3_close(db->sql);
  free(db);
}

const char *ls_db_error(const ls_db_t *db)
{
  return (db && db->error[0]) ? db->error : NULL;
}

int ls_db_schema_version(const ls_db_t *db)
{
  return db ? db->schema_version : -1;
}

/* ------------------------------------------------------------------------- */
/* Name normalisation. Must match tools/import_lensfun_xml.c byte for byte,   */
/* since that is what filled the `norm` columns being compared against.       */
/* ------------------------------------------------------------------------- */

/**
 * @brief Case-fold ASCII, drop punctuation, collapse whitespace.
 *
 * @details Deliberately ASCII-only and byte-wise: it must produce the same bytes in the
 * importer and here, on every platform, without a locale or an ICU version in the middle
 * of it. Non-ASCII bytes pass through untouched, so a UTF-8 maker name still matches
 * itself; it just does not case-fold, which no upstream name needs.
 */
size_t ls_db_normalize(const char *in, char *out, size_t out_size)
{
  if(!out || out_size == 0) return 0;
  size_t w = 0;
  int pending_space = 0;

  for(const unsigned char *p = (const unsigned char *)(in ? in : ""); *p; p++)
  {
    unsigned char c = *p;
    if(c <= ' ' || c == '_' || c == '-' || c == '/' || c == '\\' || c == ',' || c == '.'
       || c == '(' || c == ')' || c == '[' || c == ']' || c == '\'' || c == '"' || c == ':'
       || c == ';' || c == '+' || c == '*')
    {
      if(w > 0) pending_space = 1;
      continue;
    }
    if(c >= 'A' && c <= 'Z') c = (unsigned char)(c - 'A' + 'a');

    if(pending_space)
    {
      if(w + 1 >= out_size) break;
      out[w++] = ' ';
      pending_space = 0;
    }
    if(w + 1 >= out_size) break;
    out[w++] = (char)c;
  }
  out[w] = '\0';
  return w;
}

/* ------------------------------------------------------------------------- */

static int _fill_calibrations(ls_db_t *db, long long lens_id, ls_lens_t *out)
{
  static const struct
  {
    const char *sql;
    int terms;
  } q[3] = {
    { "SELECT model, focal, t0, t1, t2 FROM calib_distortion WHERE lens_id = ?1", 3 },
    { "SELECT model, focal, t0, t1, t2, t3, t4, t5 FROM calib_tca WHERE lens_id = ?1", 6 },
    { "SELECT model, focal, aperture, distance, t0, t1, t2 FROM calib_vignetting WHERE lens_id = ?1", 3 },
  };

  for(int kind = 0; kind < 3; kind++)
  {
    sqlite3_stmt *st = NULL;
    if(sqlite3_prepare_v2(db->sql, q[kind].sql, -1, &st, NULL) != SQLITE_OK)
    {
      _db_err(db, "prepare calibration");
      return -1;
    }
    sqlite3_bind_int64(st, 1, lens_id);

    int n = 0;
    int overflow = 0;
    while(sqlite3_step(st) == SQLITE_ROW)
    {
      if(n >= LS_MAX_CALIB)
      {
        overflow = 1;
        break;
      }
      const int model = sqlite3_column_int(st, 0);
      const float focal = (float)sqlite3_column_double(st, 1);

      if(kind == 0)
      {
        ls_calib_dist_t *c = &out->dist[n];
        c->model = (ls_dist_model_t)model;
        c->focal = focal;
        for(int t = 0; t < 3; t++) c->terms[t] = (float)sqlite3_column_double(st, 2 + t);
      }
      else if(kind == 1)
      {
        ls_calib_tca_t *c = &out->tca[n];
        c->model = (ls_tca_model_t)model;
        c->focal = focal;
        for(int t = 0; t < 6; t++) c->terms[t] = (float)sqlite3_column_double(st, 2 + t);
      }
      else
      {
        ls_calib_vig_t *c = &out->vig[n];
        c->model = (ls_vig_model_t)model;
        c->focal = focal;
        c->aperture = (float)sqlite3_column_double(st, 2);
        c->distance = (float)sqlite3_column_double(st, 3);
        for(int t = 0; t < 3; t++) c->terms[t] = (float)sqlite3_column_double(st, 4 + t);
      }
      n++;
    }
    sqlite3_finalize(st);

    if(overflow)
    {
      /* Truncating would corrupt both the IDW weighting and the spline neighbours, so this
       * is an error rather than a partial answer. LS_MAX_CALIB is 512 against a densest
       * real lens of 440 vignetting points; a database that trips this needs the constant
       * raised, not the rows dropped. */
      snprintf(db->error, sizeof(db->error),
               "lens %lld has more than %d calibration entries of one kind",
               lens_id, (int)LS_MAX_CALIB);
      return -1;
    }

    if(kind == 0) out->n_dist = n;
    else if(kind == 1) out->n_tca = n;
    else out->n_vig = n;
  }
  return 0;
}

int ls_db_lens_by_id(ls_db_t *db, long long lens_id, ls_lens_t *out)
{
  if(!db || !db->sql || !out) return -1;
  memset(out, 0, sizeof(*out));

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql,
                        "SELECT type, crop_factor, aspect_ratio, center_x, center_y,"
                        "       min_focal, max_focal FROM lens WHERE id = ?1",
                        -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare lens");
    return -1;
  }
  sqlite3_bind_int64(st, 1, lens_id);

  if(sqlite3_step(st) != SQLITE_ROW)
  {
    sqlite3_finalize(st);
    return 0;
  }
  out->type = (ls_lens_type_t)sqlite3_column_int(st, 0);
  out->crop_factor = (float)sqlite3_column_double(st, 1);
  out->aspect_ratio = (float)sqlite3_column_double(st, 2);
  out->center_x = (float)sqlite3_column_double(st, 3);
  out->center_y = (float)sqlite3_column_double(st, 4);
  out->min_focal = (float)sqlite3_column_double(st, 5);
  out->max_focal = (float)sqlite3_column_double(st, 6);
  sqlite3_finalize(st);

  return (_fill_calibrations(db, lens_id, out) == 0) ? 1 : -1;
}

/**
 * @brief The lens id whose names match, preferring the calibration crop nearest @p crop.
 * @return the id, 0 for no match, -1 on error.
 */
static long long _find_lens_id(ls_db_t *db, const char *maker, const char *model, float crop)
{
  char nmodel[512], nmaker[512];
  ls_db_normalize(model, nmodel, sizeof(nmodel));
  ls_db_normalize(maker, nmaker, sizeof(nmaker));

  /* One statement whether or not a maker was given: ?3 < 0 disables the maker test, so
   * there is a single query plan to reason about rather than two that can drift. */
  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql,
                        "SELECT l.id FROM lens l"
                        " JOIN lens_name m ON m.lens_id = l.id AND m.kind = 'model' AND m.norm = ?1"
                        " WHERE (?3 = 0 OR EXISTS (SELECT 1 FROM lens_name k"
                        "        WHERE k.lens_id = l.id AND k.kind = 'maker' AND k.norm = ?2))"
                        " ORDER BY abs(l.crop_factor - ?4) ASC, l.id ASC LIMIT 1",
                        -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare lens lookup");
    return -1;
  }
  sqlite3_bind_text(st, 1, nmodel, -1, SQLITE_STATIC);
  sqlite3_bind_text(st, 2, nmaker, -1, SQLITE_STATIC);
  sqlite3_bind_int(st, 3, (maker && *nmaker) ? 1 : 0);
  sqlite3_bind_double(st, 4, (crop > 0.f) ? (double)crop : 0.0);

  long long id = 0;
  if(sqlite3_step(st) == SQLITE_ROW) id = sqlite3_column_int64(st, 0);
  sqlite3_finalize(st);
  return id;
}

int ls_db_find_lens(ls_db_t *db, const char *maker, const char *model, float crop,
                    ls_lens_t *out)
{
  if(!db || !db->sql || !model || !out) return -1;

  const long long id = _find_lens_id(db, maker, model, crop);
  if(id < 0) return -1;
  if(id == 0) return 0;
  return ls_db_lens_by_id(db, id, out);
}

int ls_db_find_camera(ls_db_t *db, const char *maker, const char *model, ls_camera_t *out)
{
  if(!db || !db->sql || !model || !out) return -1;
  memset(out, 0, sizeof(*out));

  char nmodel[512], nmaker[512];
  ls_db_normalize(model, nmodel, sizeof(nmodel));
  ls_db_normalize(maker, nmaker, sizeof(nmaker));

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql,
                        "SELECT c.crop_factor, ifnull(c.mount_id, 0) FROM camera c"
                        " JOIN camera_name m ON m.camera_id = c.id AND m.kind = 'model' AND m.norm = ?1"
                        " WHERE (?3 = 0 OR EXISTS (SELECT 1 FROM camera_name k"
                        "        WHERE k.camera_id = c.id AND k.kind = 'maker' AND k.norm = ?2))"
                        " ORDER BY c.id ASC LIMIT 1",
                        -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare camera lookup");
    return -1;
  }
  sqlite3_bind_text(st, 1, nmodel, -1, SQLITE_STATIC);
  sqlite3_bind_text(st, 2, nmaker, -1, SQLITE_STATIC);
  sqlite3_bind_int(st, 3, (maker && *nmaker) ? 1 : 0);

  int found = 0;
  if(sqlite3_step(st) == SQLITE_ROW)
  {
    out->crop_factor = (float)sqlite3_column_double(st, 0);
    out->mount_id = sqlite3_column_int64(st, 1);
    found = 1;
  }
  sqlite3_finalize(st);
  return found;
}

int ls_db_lens_fits_mount(ls_db_t *db, long long lens_id, long long mount_id)
{
  if(!db || !db->sql) return -1;
  if(mount_id <= 0) return 0;

  sqlite3_stmt *st = NULL;
  /* Direct fit, or the camera's mount declares the lens's mount compatible. */
  if(sqlite3_prepare_v2(db->sql,
                        "SELECT 1 FROM lens_mount lm WHERE lm.lens_id = ?1 AND ("
                        "  lm.mount_id = ?2"
                        "  OR EXISTS (SELECT 1 FROM mount_compat mc"
                        "             WHERE mc.mount_id = ?2 AND mc.compat_id = lm.mount_id))"
                        " LIMIT 1",
                        -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare mount test");
    return -1;
  }
  sqlite3_bind_int64(st, 1, lens_id);
  sqlite3_bind_int64(st, 2, mount_id);
  const int fits = (sqlite3_step(st) == SQLITE_ROW) ? 1 : 0;
  sqlite3_finalize(st);
  return fits;
}

int ls_db_list_lenses(ls_db_t *db, long long *out_ids, int max)
{
  if(!db || !db->sql) return -1;

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql, "SELECT id FROM lens ORDER BY id ASC", -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare lens list");
    return -1;
  }
  int n = 0;
  while(sqlite3_step(st) == SQLITE_ROW)
  {
    if(out_ids)
    {
      if(n >= max) break;
      out_ids[n] = sqlite3_column_int64(st, 0);
    }
    n++;
  }
  sqlite3_finalize(st);
  return n;
}

int ls_db_lens_name(ls_db_t *db, long long lens_id, char *maker, size_t maker_size,
                    char *model, size_t model_size)
{
  if(!db || !db->sql) return -1;
  if(maker && maker_size) maker[0] = '\0';
  if(model && model_size) model[0] = '\0';

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql, "SELECT maker, model FROM lens WHERE id = ?1", -1, &st, NULL)
     != SQLITE_OK)
  {
    _db_err(db, "prepare lens name");
    return -1;
  }
  sqlite3_bind_int64(st, 1, lens_id);

  int written = 0;
  if(sqlite3_step(st) == SQLITE_ROW)
  {
    const char *a = (const char *)sqlite3_column_text(st, 0);
    const char *b = (const char *)sqlite3_column_text(st, 1);
    if(maker && maker_size && a)
    {
      snprintf(maker, maker_size, "%s", a);
      written += (int)strlen(maker);
    }
    if(model && model_size && b)
    {
      snprintf(model, model_size, "%s", b);
      written += (int)strlen(model);
    }
  }
  sqlite3_finalize(st);
  return written;
}

int ls_db_meta(ls_db_t *db, const char *key, char *out, size_t out_size)
{
  if(!db || !db->sql || !key || !out || out_size == 0) return -1;
  out[0] = '\0';

  sqlite3_stmt *st = NULL;
  if(sqlite3_prepare_v2(db->sql, "SELECT value FROM meta WHERE key = ?1", -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare meta");
    return -1;
  }
  sqlite3_bind_text(st, 1, key, -1, SQLITE_STATIC);

  int n = -1;
  if(sqlite3_step(st) == SQLITE_ROW)
  {
    const char *v = (const char *)sqlite3_column_text(st, 0);
    if(v)
    {
      snprintf(out, out_size, "%s", v);
      n = (int)strlen(out);
    }
  }
  sqlite3_finalize(st);
  return n;
}
