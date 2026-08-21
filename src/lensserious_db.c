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

#define LS_DB_SCHEMA_VERSION 2

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
    { "SELECT model, focal, t0, t1, t2 FROM calib_distortion WHERE lens_id = ?1 ORDER BY ord", 3 },
    { "SELECT model, focal, t0, t1, t2, t3, t4, t5 FROM calib_tca WHERE lens_id = ?1 ORDER BY ord", 6 },
    { "SELECT model, focal, aperture, distance, t0, t1, t2 FROM calib_vignetting WHERE lens_id = ?1 ORDER BY ord", 3 },
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

/* ------------------------------------------------------------------------- */
/* Fuzzy matching.                                                            */
/*                                                                            */
/* A raw file names a lens the way its vendor abbreviates it, and upstream     */
/* names it the way upstream chose. "16-35mm f/4G ED VR" has to reach "Nikon   */
/* AF-S Nikkor 16-35mm f/4G ED VR" with most of the tokens missing.            */
/*                                                                            */
/* The weights below are not derived from anything: they were calibrated       */
/* against liblensfun's own decisions, which tests/match_lensfun.c re-checks    */
/* over the whole database. Change one and run that test.                      */
/* ------------------------------------------------------------------------- */

enum { LS_MAX_TOKENS = 32, LS_TOKEN_LEN = 48 };

/**
 * @brief A tokenised name, with everything the comparison needs precomputed.
 *
 * @details The hash and the length are not an optimisation detail, they are most of the
 * function's cost. Comparing two names is O(tokens x tokens), and doing that with strcmp()
 * -- calling strlen() inside the inner loop, no less -- measured 610 ns per catalogue name,
 * against 82 ns for SQLite to hand the row over. Scoring was 88% of a lookup.
 *
 * `bloom` is the OR of (1 << (hash % 64)) over the tokens: if two names share no bit they
 * share no token, so the whole quadratic comparison can be skipped on one AND.
 */
typedef struct
{
  char t[LS_MAX_TOKENS][LS_TOKEN_LEN];
  unsigned h[LS_MAX_TOKENS];       /**< FNV-1a of the token */
  unsigned char len[LS_MAX_TOKENS];
  unsigned long long bloom;
  int n;
} ls_tokens_t;

/** @brief Split a normalised name on spaces, and split letter/digit runs apart.
 *
 * @details "16-35mm" normalises to "1635mm" -- no. Normalisation drops '-', so the two
 * halves of a focal range would fuse into one meaningless token and stop matching the
 * catalogue's "16 35mm". Splitting where a digit meets a letter, and where a letter meets
 * a digit, keeps "16", "35", "mm", "f", "4g" as separate comparable units. */
int ls_db_tokenize(const char *norm, char *out_tokens, int max, int stride)
{
  int n = 0;
  const char *p = norm ? norm : "";
  while(*p && n < max)
  {
    while(*p == ' ') p++;
    if(!*p) break;

    char *w = out_tokens + (size_t)n * stride;
    int len = 0;
    int prev_digit = -1;
    while(*p && *p != ' ' && len < stride - 1)
    {
      const int is_digit = (*p >= '0' && *p <= '9');
      if(prev_digit >= 0 && is_digit != prev_digit) break;   /* letter<->digit boundary */
      w[len++] = *p++;
      prev_digit = is_digit;
    }
    w[len] = '\0';
    if(len) n++;
  }
  return n;
}

static void _tokenize(const char *norm, ls_tokens_t *out)
{
  out->n = ls_db_tokenize(norm, &out->t[0][0], LS_MAX_TOKENS, LS_TOKEN_LEN);
  out->bloom = 0;
  for(int i = 0; i < out->n; i++)
  {
    unsigned hash = 2166136261u;      /* FNV-1a */
    const char *c = out->t[i];
    int len = 0;
    for(; *c; c++, len++)
    {
      hash ^= (unsigned char)*c;
      hash *= 16777619u;
    }
    out->h[i] = hash;
    out->len[i] = (unsigned char)len;
    out->bloom |= 1ULL << (hash & 63u);
  }
}

/**
 * @brief How well @p pat's tokens are covered by @p cand's, 0..100.
 *
 * @details Asymmetric on purpose. The pattern is what the camera wrote and the candidate
 * is the catalogue entry, so every pattern token ought to appear in the candidate --
 * missing one is evidence against the match. The reverse is not true: the catalogue is
 * more verbose than any EXIF field, and penalising it for that would favour the shortest
 * name in the database for every query.
 *
 * Unmatched candidate tokens are still worth a small penalty, or "35mm" would score the
 * same against every 35mm lens ever made; it is just much smaller than the forward one.
 */
static float _score_tokens(const ls_tokens_t *pat, const ls_tokens_t *cand)
{
  if(pat->n == 0 || cand->n == 0) return 0.f;

  unsigned char used[LS_MAX_TOKENS] = { 0 };
  float got = 0.f;

  for(int i = 0; i < pat->n; i++)
  {
    const unsigned hi = pat->h[i];
    const int li = pat->len[i];
    float best = 0.f;
    int best_j = -1;
    for(int j = 0; j < cand->n; j++)
    {
      if(used[j]) continue;
      float s = 0.f;
      const int lj = cand->len[j];

      /* Hash and length first. Equal tokens must agree on both, so an inequality here
       * rejects the pair in two compares instead of a call into strcmp(). */
      if(hi == cand->h[j] && li == lj && memcmp(pat->t[i], cand->t[j], (size_t)li) == 0)
        s = 1.f;
      else
      {
        /* A prefix is weak evidence: "nikkor" vs "nikkor" is a match, "af" vs "afs" is a
         * hint. Anything shorter than three characters is noise, not a hint. The lengths
         * are already known, so this costs no strlen(). */
        const int lmin = (li < lj) ? li : lj;
        /* First byte before the call. Without it memcmp() runs for EVERY pair of tokens
         * that are not equal -- ~120 calls per catalogue name -- and prefix matches are
         * rare, so almost all of them return on their first byte anyway. Doing that compare
         * inline is the difference between ~700 ns and ~200 ns per name. */
        if(lmin >= 3 && pat->t[i][0] == cand->t[j][0]
           && memcmp(pat->t[i] + 1, cand->t[j] + 1, (size_t)(lmin - 1)) == 0)
          s = 0.5f * (float)lmin / (float)((li > lj) ? li : lj);
      }
      if(s > best)
      {
        best = s;
        best_j = j;
      }
    }
    if(best_j >= 0 && best > 0.f)
    {
      used[best_j] = 1;
      got += best;
    }
  }

  const float forward = got / (float)pat->n;                    /* pattern covered */
  int unmatched = 0;
  for(int j = 0; j < cand->n; j++) if(!used[j]) unmatched++;
  const float verbosity = (float)unmatched / (float)cand->n;    /* candidate's extra words */

  float score = 100.f * (forward - 0.15f * verbosity * forward);
  if(score < 0.f) score = 0.f;
  return score;
}

int ls_db_match_lens(ls_db_t *db, const char *maker, const char *model, long long mount_id,
                     ls_db_match_t *out, int max)
{
  if(!db || !db->sql || !model || !out || max <= 0) return -1;

  char nmodel[512], nmaker[512];
  ls_db_normalize(model, nmodel, sizeof(nmodel));
  ls_db_normalize(maker, nmaker, sizeof(nmaker));

  ls_tokens_t pat, pat_maker;
  _tokenize(nmodel, &pat);
  _tokenize(nmaker, &pat_maker);
  if(pat.n == 0) return 0;

  /* Two phases, and the first one is what makes this cheap.
   *
   * Scoring every name in the catalogue costs 610 ns each -- 88% of a lookup, measured --
   * and there are ~4700 of them. The fix is not to score faster (hashing the tokens first
   * was tried and moved nothing: the per-comparison cost was already ~5 ns, there were
   * simply half a million comparisons). The fix is to score far fewer names.
   *
   * So: ask lens_token which of the QUERY's tokens is rarest, and gather only the lenses
   * carrying it. "16" or "nikkor" reaches tens of lenses where "mm" or "f" reaches
   * thousands, which is exactly why gathering on ALL the query's tokens -- also tried, also
   * measured -- is no better than the full scan it replaces.
   *
   * The frequencies come from token_df, precomputed at import. Deriving them here with a
   * GROUP BY over lens_token was the first attempt and cost 0.6 ms on its own: counting how
   * often "mm" occurs means walking every one of its index rows.
   *
   * If that yields nothing (a query whose rarest token no catalogue name shares), the scan
   * runs after all, so the answer never depends on the pruning. */
  char in_list[LS_MAX_TOKENS * (LS_TOKEN_LEN + 4) + 8];
  {
    size_t w = 0;
    in_list[w++] = '(';
    for(int i = 0; i < pat.n; i++)
    {
      if(i) in_list[w++] = ',';
      in_list[w++] = '\'';
      for(const char *c = pat.t[i]; *c; c++)
        if(*c != '\'') in_list[w++] = *c;   /* tokens are [a-z0-9] after normalisation */
      in_list[w++] = '\'';
    }
    in_list[w++] = ')';
    in_list[w] = '\0';
  }

  char rarest[LS_TOKEN_LEN] = { 0 };
  {
    char sqlbuf[sizeof(in_list) + 256];
    snprintf(sqlbuf, sizeof(sqlbuf),
             "SELECT token FROM token_df WHERE kind = 'model' AND token IN %s"
             " ORDER BY df ASC LIMIT 1", in_list);
    sqlite3_stmt *rq = NULL;
    if(sqlite3_prepare_v2(db->sql, sqlbuf, -1, &rq, NULL) == SQLITE_OK)
    {
      if(sqlite3_step(rq) == SQLITE_ROW)
      {
        const char *t = (const char *)sqlite3_column_text(rq, 0);
        if(t) snprintf(rarest, sizeof(rarest), "%s", t);
      }
      sqlite3_finalize(rq);
    }
  }

  sqlite3_stmt *st = NULL;
  char sqlbuf[1024];
  const char *mount_clause =
      (mount_id > 0)
          ? " AND EXISTS (SELECT 1 FROM lens_mount lm WHERE lm.lens_id = n.lens_id AND ("
            "   lm.mount_id = ?1 OR EXISTS (SELECT 1 FROM mount_compat mc"
            "     WHERE mc.mount_id = ?1 AND mc.compat_id = lm.mount_id)))"
          : "";
  if(rarest[0])
    snprintf(sqlbuf, sizeof(sqlbuf),
             "SELECT n.lens_id, n.norm, n.kind FROM lens_name n"
             " WHERE n.lens_id IN (SELECT lens_id FROM lens_token"
             "                     WHERE kind = 'model' AND token = ?2)%s", mount_clause);
  else
    snprintf(sqlbuf, sizeof(sqlbuf),
             "SELECT n.lens_id, n.norm, n.kind FROM lens_name n WHERE 1%s", mount_clause);

  if(sqlite3_prepare_v2(db->sql, sqlbuf, -1, &st, NULL) != SQLITE_OK)
  {
    _db_err(db, "prepare match");
    return -1;
  }
  if(rarest[0]) sqlite3_bind_text(st, 2, rarest, -1, SQLITE_STATIC);
  if(mount_id > 0) sqlite3_bind_int64(st, 1, mount_id);

  /* Best score per lens, kept in a small open-addressed table: a lens has several names
   * and several of them may score, but only its best counts. */
  enum { SLOTS = 4096 };
  long long *ids = (long long *)calloc(SLOTS, sizeof(long long));
  float *best = (float *)calloc(SLOTS, sizeof(float));
  float *maker_bonus = (float *)calloc(SLOTS, sizeof(float));
  if(!ids || !best || !maker_bonus)
  {
    free(ids); free(best); free(maker_bonus);
    sqlite3_finalize(st);
    return -1;
  }

  while(sqlite3_step(st) == SQLITE_ROW)
  {
    const long long id = sqlite3_column_int64(st, 0);
    const char *norm = (const char *)sqlite3_column_text(st, 1);
    const char *kind = (const char *)sqlite3_column_text(st, 2);
    if(!norm || !kind) continue;

    ls_tokens_t cand;
    _tokenize(norm, &cand);

    /* No shared token means no score, and the quadratic comparison below would only
     * discover that the expensive way. One AND settles it for most of the catalogue. */
    if(!(cand.bloom & ((kind[0] == 'm' && kind[1] == 'o') ? pat.bloom : pat_maker.bloom)))
      continue;

    size_t slot = (size_t)id % SLOTS;
    while(ids[slot] && ids[slot] != id) slot = (slot + 1) % SLOTS;
    ids[slot] = id;

    if(kind[0] == 'm' && kind[1] == 'o')          /* "model" */
    {
      const float s = _score_tokens(&pat, &cand);
      if(s > best[slot]) best[slot] = s;
    }
    else if(pat_maker.n)                          /* "maker" */
    {
      /* The maker is corroboration, not a filter: vendors and upstream disagree about
       * their own names often enough ("Nikon Corporation" vs "Nikon") that requiring it
       * loses more matches than the false positives it prevents. */
      const float s = _score_tokens(&pat_maker, &cand);
      if(s > maker_bonus[slot]) maker_bonus[slot] = s;
    }
  }
  sqlite3_finalize(st);

  int n = 0;
  for(size_t slot = 0; slot < SLOTS; slot++)
  {
    if(!ids[slot] || best[slot] <= 0.f) continue;
    const float score = best[slot] + 0.10f * maker_bonus[slot];

    /* Insertion sort into the caller's top-N: max is small (a GUI shows a handful). */
    int at = n;
    while(at > 0 && out[at - 1].score < score)
    {
      if(at < max) out[at] = out[at - 1];
      at--;
    }
    if(at < max)
    {
      out[at].lens_id = ids[slot];
      out[at].score = score;
      if(n < max) n++;
    }
  }

  free(ids);
  free(best);
  free(maker_bonus);
  return n;
}
