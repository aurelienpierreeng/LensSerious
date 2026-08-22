/*
    LensSerious — the offline importer, as a function rather than a program.

    Copyright (C) 2026 Aurélien PIERRE.
    License: LGPL-3.0-or-later (see LICENSE).
*/

/** @file lensserious_import.h
 *
 * @brief Turn a lensfun database into a LensSerious one.
 *
 * @details This is the only part of the project that links liblensfun, and it lives in its
 * own target (`lensserious_import`) for exactly that reason: a consumer that renders pixels
 * links `lensserious` and `lensserious_db` and gets no lensfun dependency with them. Nothing
 * in a pixel pipeline should reach this header.
 *
 * It exists as a function because there are two callers. One is the offline tool in this
 * repository, run when upstream publishes calibrations. The other is a consumer's own
 * updater -- Ansel ships `ansel-lens-db-update` -- which needs the same conversion with
 * different default paths. A second copy of three hundred lines of schema-filling is a
 * second copy to drift, and the drift would be silent: both would produce a database, and
 * only one of them would be right.
 *
 * @see tools/import_lensfun_xml.c, which is now a five-line front end to this.
 */

#ifndef LENSSERIOUS_IMPORT_H
#define LENSSERIOUS_IMPORT_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Read a lensfun database and write it as a LensSerious SQLite file.
 *
 * @param schema_path the SQL schema to create the file from (`src/schema.sql`).
 * @param out_path where to write. The file is written to a temporary and RENAMED into
 * place, so a reader opening it with `immutable=1` never sees a partial database and a
 * killed run cannot leave one that looks complete. See lensserious_db.h.
 * @param base_xml_dir a directory of lensfun XML to load FIRST, or NULL. Everything loaded
 * afterwards can override a lens defined here, and a lens nobody else defines survives --
 * lensfun's own "later definitions override earlier" rule, applied by lensfun.
 *
 * This is what makes an update additive rather than a replacement. A consumer that ships a
 * database passes its own calibrations here and the machine's profiles are merged over the
 * top; without it, a machine with no system-wide lensfun yields a database holding only the
 * few profiles its user wrote by hand, which is a catastrophic thing to install in place of
 * fifteen hundred lenses.
 * @param xml_dir a directory of lensfun XML to read, or NULL.
 *
 * NULL is the interesting value, and is what a consumer's updater wants: it calls
 * lensfun's own `lf_db_load()`, which searches all four standard locations -- the system
 * database, the system updates directory, the user's updates directory, and the user's own
 * `~/.local/share/lensfun`. That last one is where hand-written profiles live, and the
 * aggregation and override rules between the four are lensfun's, not this project's
 * approximation of them. Passing a directory instead reads exactly that directory, which is
 * what a reproducible package build wants, since a build machine has no user profiles.
 *
 * @return 0 on success, non-zero on failure (already reported on stderr).
 *
 * @warning Not thread-safe and not re-entrant: this is a batch conversion that exits the
 * process on a database error, which is the right contract for the two command-line tools
 * that call it and the wrong one for anything else.
 */
int ls_import_run(const char *schema_path, const char *out_path,
                  const char *base_xml_dir, const char *xml_dir);

#ifdef __cplusplus
}
#endif

#endif /* LENSSERIOUS_IMPORT_H */
