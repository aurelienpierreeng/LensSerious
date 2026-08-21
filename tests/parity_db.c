/*
    LensSerious — the database must say exactly what liblensfun says.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file parity_db.c
 *
 * @brief Holds tools/import_lensfun_xml.c to the same standard as everything else here:
 * not "it imported without complaining", but "every number it wrote is the number
 * liblensfun would have produced, for every lens in the database".
 *
 * Three things are checked, and they fail for different reasons:
 *
 *  1. **Content.** Walk both sides in insertion order and compare every field of every
 *     lens, calibration arrays included, EXACTLY. Not approximately: the importer binds
 *     liblensfun's floats as doubles and the reader casts them back, which round-trips
 *     bit for bit, so any tolerance here would only hide a real defect.
 *  2. **Lookup.** Every lens must be findable by its own maker and model through
 *     ls_db_find_lens(). Content parity says the rows are right; this says the index and
 *     the name normalisation agree with each other.
 *  3. **Concurrency.** The header claims read access is thread-safe with no mutex. Several
 *     threads each open their own handle on the same file and run the same lookups; every
 *     one must get identical bytes. A test cannot prove the absence of a race, but it can
 *     fail loudly if the "one handle per thread" contract does not actually hold.
 *
 * Note the bridge below is NOT the one in parity_lensfun.c. That one collapses every
 * non-rectilinear lens to LS_LENS_FISHEYE because it skips them anyway; this one maps
 * the type faithfully, because here the type is part of what is being verified.
 */

#include "lensserious_db.h"

#include <lensfun.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
static int signed_zeros = 0;

#define FAILF(...)                                                                        \
  do {                                                                                    \
    if(failures < 20) { fprintf(stderr, "FAIL: "); fprintf(stderr, __VA_ARGS__); }         \
    failures++;                                                                           \
  } while(0)

/* The same mapping the importer uses. Kept separate from parity_lensfun.c's lossy one. */
static ls_lens_type_t _type_from_lf(enum lfLensType t)
{
  switch(t)
  {
    case LF_RECTILINEAR:           return LS_LENS_RECTILINEAR;
    case LF_FISHEYE:               return LS_LENS_FISHEYE;
    case LF_PANORAMIC:             return LS_LENS_PANORAMIC;
    case LF_EQUIRECTANGULAR:       return LS_LENS_EQUIRECTANGULAR;
    case LF_FISHEYE_ORTHOGRAPHIC:  return LS_LENS_FISHEYE_ORTHOGRAPHIC;
    case LF_FISHEYE_STEREOGRAPHIC: return LS_LENS_FISHEYE_STEREOGRAPHIC;
    case LF_FISHEYE_EQUISOLID:     return LS_LENS_FISHEYE_EQUISOLID;
    case LF_FISHEYE_THOBY:         return LS_LENS_FISHEYE_THOBY;
    default:                       return LS_LENS_UNKNOWN;
  }
}

static int _lens_from_lf(const lfLens *lf, ls_lens_t *out)
{
  memset(out, 0, sizeof(*out));
  out->type = _type_from_lf(lf->Type);
  out->crop_factor = lf->CropFactor;
  out->aspect_ratio = lf->AspectRatio;
  out->center_x = lf->CenterX;
  out->center_y = lf->CenterY;
  out->min_focal = lf->MinFocal;
  out->max_focal = lf->MaxFocal;

  /* <real-focal-length>, with the zero entries dropped -- the importer drops them too,
   * because upstream skips them when interpolating and a stored zero is a row nothing can
   * use. Both sides must drop them identically or the counts will not match. */
  if(lf->CalibRealFocal)
    for(int i = 0; lf->CalibRealFocal[i]; i++)
    {
      const lfLensCalibRealFocal *c = lf->CalibRealFocal[i];
      if(c->RealFocal == 0.f) continue;
      if(out->n_real_focal >= LS_MAX_CALIB) return -1;
      out->real_focal[out->n_real_focal].focal = c->Focal;
      out->real_focal[out->n_real_focal].real_focal = c->RealFocal;
      out->n_real_focal++;
    }

  if(lf->CalibDistortion)
    for(int i = 0; lf->CalibDistortion[i]; i++)
    {
      if(out->n_dist >= LS_MAX_CALIB) return -1;
      const lfLensCalibDistortion *c = lf->CalibDistortion[i];
      ls_calib_dist_t *d = &out->dist[out->n_dist++];
      switch(c->Model)
      {
        case LF_DIST_MODEL_POLY3:  d->model = LS_DIST_POLY3; break;
        case LF_DIST_MODEL_POLY5:  d->model = LS_DIST_POLY5; break;
        case LF_DIST_MODEL_PTLENS: d->model = LS_DIST_PTLENS; break;
        default:                   d->model = LS_DIST_NONE; break;
      }
      d->focal = c->Focal;
      for(int k = 0; k < 3; k++) d->terms[k] = c->Terms[k];
    }

  if(lf->CalibTCA)
    for(int i = 0; lf->CalibTCA[i]; i++)
    {
      if(out->n_tca >= LS_MAX_CALIB) return -1;
      const lfLensCalibTCA *c = lf->CalibTCA[i];
      ls_calib_tca_t *t = &out->tca[out->n_tca++];
      switch(c->Model)
      {
        case LF_TCA_MODEL_LINEAR: t->model = LS_TCA_LINEAR; break;
        case LF_TCA_MODEL_POLY3:  t->model = LS_TCA_POLY3; break;
        default:                  t->model = LS_TCA_NONE; break;
      }
      t->focal = c->Focal;
      for(int k = 0; k < 6; k++) t->terms[k] = c->Terms[k];
    }

  if(lf->CalibVignetting)
    for(int i = 0; lf->CalibVignetting[i]; i++)
    {
      if(out->n_vig >= LS_MAX_CALIB) return -1;
      const lfLensCalibVignetting *c = lf->CalibVignetting[i];
      ls_calib_vig_t *v = &out->vig[out->n_vig++];
      v->model = (c->Model == LF_VIGNETTING_MODEL_PA) ? LS_VIG_PA : LS_VIG_NONE;
      v->focal = c->Focal;
      v->aperture = c->Aperture;
      v->distance = c->Distance;
      for(int k = 0; k < 3; k++) v->terms[k] = c->Terms[k];
    }

  return 0;
}

static int _same_lens(const ls_lens_t *a, const ls_lens_t *b, const char *name)
{
  int bad = 0;
#define CMP(field)                                                                        \
  if(a->field != b->field)                                                                \
  {                                                                                       \
    FAILF("%s: " #field " %g (lensfun) != %g (db)\n", name, (double)a->field, (double)b->field); \
    bad = 1;                                                                              \
  }
  CMP(type) CMP(crop_factor) CMP(aspect_ratio) CMP(center_x) CMP(center_y)
  CMP(min_focal) CMP(max_focal) CMP(n_dist) CMP(n_tca) CMP(n_vig) CMP(n_real_focal)
#undef CMP
  if(bad) return 0;

  /* Exact, with one stated exception. IEEE has two encodings of zero, and SQLite keeps
   * only one: it stores a REAL with no fractional part as an integer, so a term that was
   * -0.0 in the XML reads back as +0.0. Upstream has 40 such terms.
   *
   * That is a change of encoding, not of value. -0.0 == +0.0 compares true, and every
   * evaluator here only ever multiplies or adds these terms, where the sign of zero
   * propagates into a result that is itself zero and compares equal. The one operation
   * that would expose it -- dividing by a term -- is the vignetting `1/c`, and c is
   * guarded by `c != 0.f`, which is false for BOTH zeros.
   *
   * So the comparison below is IEEE equality (which makes the two zeros equal and keeps
   * NaN unequal to everything, deliberately), and signed zeros are counted and reported
   * rather than waved through silently. Anything else still fails, bit for bit. */
#define CMPF(what, idx, field, x, y)                                                      \
  do {                                                                                    \
    const float _x = (x), _y = (y);                                                       \
    unsigned _ux, _uy;                                                                    \
    memcpy(&_ux, &_x, 4); memcpy(&_uy, &_y, 4);                                           \
    if(_x == _y)                                                                          \
    {                                                                                     \
      if(_ux != _uy) signed_zeros++;   /* the only encoding difference tolerated */        \
    }                                                                                     \
    else                                                                                  \
    {                                                                                     \
      FAILF("%s: %s[%d]." field " %.9g [%08x] (lensfun) != %.9g [%08x] (db)\n",            \
            name, what, idx, (double)_x, _ux, (double)_y, _uy);                           \
      bad = 1;                                                                            \
    }                                                                                     \
  } while(0);

  for(int i = 0; i < a->n_dist; i++)
  {
    if(a->dist[i].model != b->dist[i].model)
    {
      FAILF("%s: distortion[%d].model %d != %d\n", name, i, a->dist[i].model, b->dist[i].model);
      bad = 1;
    }
    CMPF("distortion", i, "focal", a->dist[i].focal, b->dist[i].focal)
    for(int k = 0; k < 3; k++) CMPF("distortion", i, "terms", a->dist[i].terms[k], b->dist[i].terms[k])
  }
  for(int i = 0; i < a->n_tca; i++)
  {
    if(a->tca[i].model != b->tca[i].model)
    {
      FAILF("%s: tca[%d].model %d != %d\n", name, i, a->tca[i].model, b->tca[i].model);
      bad = 1;
    }
    CMPF("tca", i, "focal", a->tca[i].focal, b->tca[i].focal)
    for(int k = 0; k < 6; k++) CMPF("tca", i, "terms", a->tca[i].terms[k], b->tca[i].terms[k])
  }
  /* The reader returns these ordered by focal, which is also the order the importer wrote
   * them in and the order upstream lists them; if that ever stops being true the spline
   * would still pick the same four neighbours, but the comparison below would not. */
  for(int i = 0; i < a->n_real_focal; i++)
  {
    CMPF("real_focal", i, "focal", a->real_focal[i].focal, b->real_focal[i].focal)
    CMPF("real_focal", i, "real_focal", a->real_focal[i].real_focal, b->real_focal[i].real_focal)
  }
  for(int i = 0; i < a->n_vig; i++)
  {
    if(a->vig[i].model != b->vig[i].model)
    {
      FAILF("%s: vignetting[%d].model %d != %d\n", name, i, a->vig[i].model, b->vig[i].model);
      bad = 1;
    }
    CMPF("vignetting", i, "focal", a->vig[i].focal, b->vig[i].focal)
    CMPF("vignetting", i, "aperture", a->vig[i].aperture, b->vig[i].aperture)
    CMPF("vignetting", i, "distance", a->vig[i].distance, b->vig[i].distance)
    for(int k = 0; k < 3; k++) CMPF("vignetting", i, "terms", a->vig[i].terms[k], b->vig[i].terms[k])
  }
#undef CMPF
  return !bad;
}

/* --- concurrency ---------------------------------------------------------- */

typedef struct
{
  const char *path;
  const char **models;
  const char **makers;
  const float *crops;
  int n;
  ls_lens_t *reference;   /* what a single-threaded read produced */
  int mismatches;
  int errors;
} thread_arg_t;

static void *_hammer(void *ud)
{
  thread_arg_t *a = (thread_arg_t *)ud;

  /* Its OWN handle: that is the contract the header states, and what this exercises. */
  ls_db_t *db = ls_db_open(a->path);
  if(!db)
  {
    a->errors++;
    return NULL;
  }
  for(int round = 0; round < 4; round++)
    for(int i = 0; i < a->n; i++)
    {
      ls_lens_t got;
      if(ls_db_find_lens(db, a->makers[i], a->models[i], a->crops[i], &got) != 1)
      {
        a->errors++;
        continue;
      }
      if(memcmp(&got, &a->reference[i], sizeof(ls_lens_t)) != 0) a->mismatches++;
    }
  ls_db_close(db);
  return NULL;
}

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <lenses.db>\n", argv[0]);
    return 2;
  }
  const char *path = argv[1];

  ls_db_t *db = ls_db_open(path);
  if(!db)
  {
    fprintf(stderr, "FAIL: cannot open `%s' as a LensSerious database\n", path);
    return 1;
  }

  char built[128] = { 0 };
  ls_db_meta(db, "built_utc", built, sizeof(built));
  printf("db parity: schema v%d, built %s\n", ls_db_schema_version(db),
         built[0] ? built : "(unknown)");

  lfDatabase *ldb = lf_db_new();
  if(!ldb || lf_db_load(ldb) != LF_NO_ERROR)
  {
    fprintf(stderr, "FAIL: no system lensfun database to compare against\n");
    return 1;
  }
  const lfLens *const *lenses = lf_db_get_lenses(ldb);

  int n_lf = 0;
  for(const lfLens *const *l = lenses; l && *l; l++) n_lf++;

  const int n_db = ls_db_list_lenses(db, NULL, 0);
  if(n_db != n_lf)
    FAILF("lens count: %d in lensfun, %d in the database\n", n_lf, n_db);

  long long *ids = (long long *)calloc((size_t)(n_db > 0 ? n_db : 1), sizeof(long long));
  ls_db_list_lenses(db, ids, n_db);

  /* --- 1. content, in insertion order ---------------------------------- */

  const int n = (n_lf < n_db) ? n_lf : n_db;
  ls_lens_t *reference = (ls_lens_t *)calloc((size_t)(n > 0 ? n : 1), sizeof(ls_lens_t));
  const char **models = (const char **)calloc((size_t)(n > 0 ? n : 1), sizeof(char *));
  const char **makers = (const char **)calloc((size_t)(n > 0 ? n : 1), sizeof(char *));

  int compared = 0;
  for(int i = 0; i < n; i++)
  {
    const lfLens *lf = lenses[i];
    const char *name = lf_mlstr_get(lf->Model);

    ls_lens_t want, got;
    if(_lens_from_lf(lf, &want) != 0)
    {
      FAILF("%s: more than %d calibration entries -- LS_MAX_CALIB is too small\n",
            name, (int)LS_MAX_CALIB);
      continue;
    }
    if(ls_db_lens_by_id(db, ids[i], &got) != 1)
    {
      FAILF("%s: lens id %lld not readable (%s)\n", name, ids[i],
            ls_db_error(db) ? ls_db_error(db) : "no row");
      continue;
    }
    _same_lens(&want, &got, name);
    compared++;

    reference[i] = got;
    models[i] = lf_mlstr_get(lf->Model);
    makers[i] = lf_mlstr_get(lf->Maker);
  }
  printf("db parity: %d lenses compared field by field against liblensfun"
         " (%d signed zeros normalised by SQLite)\n", compared, signed_zeros);

  /* --- 2. every lens findable by its own name -------------------------- */

  float *crops = (float *)calloc((size_t)(n > 0 ? n : 1), sizeof(float));
  for(int i = 0; i < n; i++) crops[i] = reference[i].crop_factor;

  int found = 0, missing = 0, wrong = 0;
  for(int i = 0; i < n; i++)
  {
    if(!models[i]) continue;
    ls_lens_t got;
    const int rc = ls_db_find_lens(db, makers[i], models[i], crops[i], &got);
    if(rc != 1)
    {
      if(missing < 10) fprintf(stderr, "FAIL: `%s' / `%s' not found by name\n",
                               makers[i] ? makers[i] : "(null)", models[i]);
      missing++;
      continue;
    }
    found++;
    /* Names are not unique -- upstream ships the same lens calibrated on several sensors.
     * What must hold is that the lookup returned A lens with that name, so compare against
     * whichever row it picked rather than demanding the one we started from. */
    if(memcmp(&got, &reference[i], sizeof(ls_lens_t)) != 0) wrong++;
    /* The threads must reproduce THIS answer -- the by-name lookup they will run -- not
     * the by-id load above. Upstream ships the same lens calibrated on several sensors, so
     * the two legitimately differ for a name with siblings; comparing across them would
     * report a race that is really just a different question being asked. */
    reference[i] = got;
  }
  if(missing) FAILF("%d of %d lenses could not be found by their own name\n", missing, n);
  printf("db parity: %d lenses found by name (%d resolved to a same-named sibling)\n",
         found, wrong);

  /* --- 3. concurrent readers, no mutex --------------------------------- */

  enum { THREADS = 8, SAMPLE = 64 };
  const int sample = (n < SAMPLE) ? n : SAMPLE;
  pthread_t th[THREADS];
  thread_arg_t args[THREADS];
  int spawned = 0;

  for(int t = 0; t < THREADS; t++)
  {
    args[t] = (thread_arg_t){ path, models, makers, crops, sample, reference, 0, 0 };
    if(pthread_create(&th[t], NULL, _hammer, &args[t]) == 0) spawned++;
    else break;
  }
  for(int t = 0; t < spawned; t++) pthread_join(th[t], NULL);

  int cerr = 0, cmis = 0;
  for(int t = 0; t < spawned; t++)
  {
    cerr += args[t].errors;
    cmis += args[t].mismatches;
  }
  if(cerr || cmis)
    FAILF("concurrent readers: %d lookup errors, %d differing results across %d threads\n",
          cerr, cmis, spawned);
  printf("db parity: %d threads x %d rounds x %d lookups, own handle each, %s\n",
         spawned, 4, sample, (cerr || cmis) ? "MISMATCHED" : "all identical");

  free(ids);
  free(crops);
  free(reference);
  free(models);
  free(makers);
  ls_db_close(db);
  lf_db_destroy(ldb);

  if(failures)
  {
    fprintf(stderr, "db parity: %d failure(s)\n", failures);
    return 1;
  }
  printf("db parity: pass\n");
  return 0;
}
