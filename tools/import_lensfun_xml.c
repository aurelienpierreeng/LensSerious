/*
    LensSerious — the periodic rebuild, as a command.

    Copyright (C) 2026 Aurélien PIERRE.
    License: LGPL-3.0-or-later (see LICENSE).
*/

/** @file import_lensfun_xml.c
 *
 * @brief Command-line front end to ls_import_run().
 *
 * @details The conversion itself is in src/lensserious_import.c, because this is not its
 * only caller: a consumer shipping the database also has to be able to rebuild it on a
 * user's machine, against that machine's own lensfun installation and their own profiles.
 * Both go through the same function.
 */

#include "lensserious_import.h"

#include <stdio.h>

int main(int argc, char **argv)
{
  if(argc < 3)
  {
    fprintf(stderr,
            "usage: %s <schema.sql> <out.db> [lensfun-xml-dir]\n"
            "\n"
            "Reads the lensfun database and writes it as a LensSerious SQLite database.\n"
            "The output is written to a temporary file and renamed into place, so readers\n"
            "using immutable=1 never see a partial file -- see lensserious_db.h.\n"
            "\n"
            "With no directory, lensfun's own lf_db_load() is used, which searches all\n"
            "four standard locations and applies lensfun's own override rules -- including\n"
            "the user's hand-written profiles. Give a directory to read exactly that one,\n"
            "which is what a reproducible package build wants.\n",
            argv[0]);
    return 2;
  }

  return ls_import_run(argv[1], argv[2], NULL, (argc > 3) ? argv[3] : NULL);
}
