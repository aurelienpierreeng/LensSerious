/*
    LensSerious — does the matcher decide what lensfun decides?

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file match_lensfun.c
 *
 * @brief lensfun's scorer could not be ported: its sources are not part of an installed
 * lensfun, only the headers and the shared object. So this does not check that
 * ls_db_match_lens() is the same ALGORITHM -- it checks that it reaches the same ANSWER,
 * which is the only part a user can observe.
 *
 * liblensfun is the oracle. For every lens in the database, both are asked the same
 * question and their top pick compared. The queries are the shapes a raw file actually
 * produces, generated mechanically from each catalogue name:
 *
 *   full      the catalogue name, verbatim -- the easy case, and a floor: any failure here
 *             is a defect in normalisation, not in scoring
 *   nomaker   the name with its leading maker word removed, which is what most vendors
 *             write in EXIF
 *   tail      the last four tokens, standing in for the terse "16-35mm f/4G ED VR" form
 *   lower     lower-cased and stripped of punctuation
 *
 * Agreement is reported per shape rather than as one number, because the shapes fail for
 * different reasons and averaging them hides which.
 *
 * A disagreement is not automatically a defect on this side: upstream ships the same lens
 * calibrated on several sensors, so two ids can be equally right. Those are counted
 * separately, by comparing the NAMES the two sides landed on rather than the ids.
 */

#include "lensserious_db.h"

#include <lensfun.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
  const char *label;
  int asked, agree, same_name, miss, wrong;
} shape_t;

static void _strip_leading_word(const char *in, char *out, size_t n)
{
  const char *sp = strchr(in, ' ');
  snprintf(out, n, "%s", sp ? sp + 1 : in);
}

static void _last_tokens(const char *in, int want, char *out, size_t n)
{
  const char *starts[64];
  int count = 0;
  for(const char *p = in; *p && count < 64;)
  {
    while(*p == ' ') p++;
    if(!*p) break;
    starts[count++] = p;
    while(*p && *p != ' ') p++;
  }
  const int from = (count > want) ? count - want : 0;
  snprintf(out, n, "%s", count ? starts[from] : in);
}

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <lenses.db>\n", argv[0]);
    return 2;
  }

  ls_db_t *db = ls_db_open(argv[1]);
  if(!db)
  {
    fprintf(stderr, "FAIL: cannot open `%s'\n", argv[1]);
    return 1;
  }
  lfDatabase *ldb = lf_db_new();
  if(!ldb || lf_db_load(ldb) != LF_NO_ERROR)
  {
    fprintf(stderr, "FAIL: no system lensfun database\n");
    return 1;
  }
  const lfLens *const *lenses = lf_db_get_lenses(ldb);

  shape_t shapes[4] = {
    { "full   ", 0, 0, 0, 0, 0 },
    { "nomaker", 0, 0, 0, 0, 0 },
    { "tail   ", 0, 0, 0, 0, 0 },
    { "lower  ", 0, 0, 0, 0, 0 },
  };

  int shown = 0;
  for(const lfLens *const *l = lenses; l && *l; l++)
  {
    const lfLens *lf = *l;
    const char *maker = lf_mlstr_get(lf->Maker);
    const char *model = lf_mlstr_get(lf->Model);
    if(!model || !*model) continue;

    char q[4][512];
    snprintf(q[0], sizeof(q[0]), "%s", model);
    _strip_leading_word(model, q[1], sizeof(q[1]));
    _last_tokens(model, 4, q[2], sizeof(q[2]));
    ls_db_normalize(model, q[3], sizeof(q[3]));

    for(int k = 0; k < 4; k++)
    {
      shape_t *sh = &shapes[k];
      sh->asked++;

      /* the oracle */
      const lfLens **lf_hits = lf_db_find_lenses_hd(ldb, NULL, maker, q[k], 0);
      const char *lf_pick = (lf_hits && lf_hits[0]) ? lf_mlstr_get(lf_hits[0]->Model) : NULL;

      /* us */
      ls_db_match_t got[4];
      const int n = ls_db_match_lens(db, maker, q[k], 0, 0.f, got, 4);
      char ls_maker[256] = { 0 }, ls_model[256] = { 0 };
      if(n > 0) ls_db_lens_name(db, got[0].lens_id, ls_maker, sizeof(ls_maker),
                                ls_model, sizeof(ls_model));

      if(!lf_pick)
      {
        /* lensfun found nothing; nothing to agree or disagree with */
        if(lf_hits) lf_free(lf_hits);
        sh->asked--;
        continue;
      }
      if(n <= 0)
      {
        sh->miss++;
        if(shown < 12)
        {
          fprintf(stderr, "MISS  [%s] `%s' -> lensfun `%s', us (nothing)\n",
                  sh->label, q[k], lf_pick);
          shown++;
        }
      }
      else if(strcmp(ls_model, lf_pick) == 0)
      {
        sh->agree++;
      }
      else
      {
        /* Different name. Upstream duplicates a lens across sensors, so check whether the
         * two names are the same string under normalisation before calling it wrong. */
        char a[512], b[512];
        ls_db_normalize(ls_model, a, sizeof(a));
        ls_db_normalize(lf_pick, b, sizeof(b));
        if(strcmp(a, b) == 0)
          sh->same_name++;
        else
        {
          sh->wrong++;
          if(shown < 12)
          {
            fprintf(stderr, "DIFF  [%s] `%s'\n         lensfun: `%s'\n         us     : `%s' (%.1f)\n",
                    sh->label, q[k], lf_pick, ls_model, (double)got[0].score);
            shown++;
          }
        }
      }
      if(lf_hits) lf_free(lf_hits);
    }
  }

  /* --- phase two: the same question, with a CAMERA ---------------------------------
   *
   * Everything above passes NULL for lensfun's camera and 0 for our crop factor, so it
   * never exercises the rule that picks BETWEEN lenses of the same name. Upstream rejects a
   * calibration measured on a sensor more than 4% larger than the camera's and grades the
   * rest by how closely the two match, and several lenses in the database differ ONLY in
   * that. Without this phase the matcher can return a name-identical lens carrying the
   * wrong calibration, which is what happened on a real Nikon D5300: a full-frame row
   * instead of the 1.528 one, and 0.19 px of geometry error that no name comparison could
   * see. */
  const lfCamera *const *cameras = lf_db_get_cameras(ldb);
  int crop_asked = 0, crop_agree = 0, crop_shown = 0;
  for(int i = 0; lenses && lenses[i]; i++)
  {
    const lfLens *lf = lenses[i];
    const char *maker = lf_mlstr_get(lf->Maker);
    const char *model = lf_mlstr_get(lf->Model);
    if(!model || !*model || lf->CropFactor <= 0.f) continue;

    /* A REAL camera, the one whose sensor is closest to what this lens was calibrated on.
     * Synthesising an lfCamera would need its Mount set, and the setter is C++-only; using
     * a real one also makes this the path Ansel actually takes -- resolve the camera from
     * EXIF, then resolve the lens against that camera. */
    const lfCamera *cam = NULL;
    float best_d = 1e9f;
    for(int c = 0; cameras && cameras[c]; c++)
    {
      const float d = fabsf(cameras[c]->CropFactor - lf->CropFactor);
      if(d < best_d) { best_d = d; cam = cameras[c]; }
    }
    if(!cam) continue;

    /* Our side resolves the same camera by name first, exactly as a consumer does. */
    ls_camera_t ls_cam;
    if(ls_db_find_camera(db, lf_mlstr_get(cam->Maker), lf_mlstr_get(cam->Model), &ls_cam) != 1)
      continue;

    const lfLens **lf_hits = lf_db_find_lenses_hd(ldb, cam, maker, model, 0);
    ls_db_match_t got[1];
    const int n = ls_db_match_lens(db, maker, model, ls_cam.mount_id, ls_cam.crop_factor,
                                   got, 1);

    if(lf_hits && lf_hits[0] && n > 0)
    {
      crop_asked++;
      ls_lens_t ours;
      const float theirs = lf_hits[0]->CropFactor;
      const float mine = (ls_db_lens_by_id(db, got[0].lens_id, &ours) == 1)
                             ? ours.crop_factor : 0.f;
      /* The CALIBRATION crop of the pick, not its name: two rows can share a name and this
       * is the field that distinguishes them, so it is the field to compare. */
      if(theirs == mine) crop_agree++;
      else if(crop_shown < 8)
      {
        fprintf(stderr, "CROP  `%s' on a %.3f body: lensfun picked crop %.3f, us %.3f\n",
                model, (double)lf->CropFactor, (double)theirs, (double)mine);
        crop_shown++;
      }
    }
    if(lf_hits) lf_free(lf_hits);
  }
  printf("\nmatch: with a camera, the calibration crop of the pick agrees %d/%d (%.1f%%)\n",
         crop_agree, crop_asked, crop_asked ? 100.0 * crop_agree / crop_asked : 0.0);

  printf("\nmatch: agreement with liblensfun's top pick, whole database\n");
  int total_asked = 0, total_ok = 0;
  for(int k = 0; k < 4; k++)
  {
    const shape_t *s = &shapes[k];
    if(!s->asked) continue;
    const int ok = s->agree + s->same_name;
    total_asked += s->asked;
    total_ok += ok;
    printf("  %s  %5d asked  %5d agree (%5.1f%%)  %4d same-name  %4d differ  %4d missed\n",
           s->label, s->asked, ok, 100.0 * ok / s->asked, s->same_name, s->wrong, s->miss);
  }
  printf("  overall  %5d asked  %5d agree (%.1f%%)\n",
         total_asked, total_ok, total_asked ? 100.0 * total_ok / total_asked : 0.0);

  ls_db_close(db);
  lf_db_destroy(ldb);

  /* A floor, not a target. The verbatim shape must be essentially perfect -- anything less
   * is a normalisation bug -- and the abbreviated shapes are where a scorer earns its
   * keep. These thresholds exist so a regression fails the build; raising them is the
   * point of any further work on the scorer. */
  const double overall = total_asked ? 100.0 * total_ok / total_asked : 0.0;
  const double full_rate = shapes[0].asked
      ? 100.0 * (shapes[0].agree + shapes[0].same_name) / shapes[0].asked : 0.0;
  int bad = 0;
  /* The crop rule is what picks BETWEEN lenses of the same name, and getting it wrong
   * returns a name-identical lens carrying the wrong calibration -- invisible to every
   * name-based check above. It was found by running a real raw through Ansel, not here,
   * which is why the floor is high: this phase exists to keep it found. */
  /* Not 100%: two entries in the current database are cases where liblensfun picks a
   * DIFFERENTLY NAMED lens and this picks the exact match -- "Sigma 70-200mm f/2.8 EX DG
   * HSM" against upstream's choice of the EX DG *OS* HSM, and the same shape for a Tokina
   * 11-20. Disagreeing with the oracle there is the right answer, so the floor sits below
   * them rather than the code being bent to reproduce them. */
  const double crop_rate = crop_asked ? 100.0 * crop_agree / crop_asked : 0.0;
  if(crop_rate < 99.0)
  {
    fprintf(stderr, "FAIL: with a camera, the pick's calibration crop agrees only %.1f%%\n",
            crop_rate);
    bad = 1;
  }
  if(full_rate < 99.0)
  {
    fprintf(stderr, "FAIL: verbatim names agree only %.1f%% of the time\n", full_rate);
    bad = 1;
  }
  if(overall < 90.0)
  {
    fprintf(stderr, "FAIL: overall agreement %.1f%% is below the 90%% floor\n", overall);
    bad = 1;
  }
  return bad;
}
