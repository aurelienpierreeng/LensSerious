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
      const int n = ls_db_match_lens(db, maker, q[k], 0, got, 4);
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
