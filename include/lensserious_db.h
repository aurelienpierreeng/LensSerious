/*
    LensSerious — reading the lens database.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file lensserious_db.h
 *
 * @brief A read-only database, and an API with nothing behind it.
 *
 * @details What this deliberately is not: lensfun's `lfDatabase`, a process-wide object
 * that must be constructed before anything can be asked, that parses ~8 MB of XML into
 * 1051 cameras and 1562 lenses on the way up (measured 89-102 ms and +4 MB resident), and
 * that every caller must then serialise against. In Ansel that object became one global
 * plus one mutex held across every lookup, which is precisely the shape this library
 * exists to remove.
 *
 * Here the database is a file, and this is how you read it:
 *
 *  - **Stateless.** No initialisation, no shutdown, no global. Every function is a pure
 *    function of the file and its arguments, and every result is a plain value the caller
 *    owns outright. Nothing is cached between calls, so nothing can go stale.
 *  - **Read-only, and lock-free.** The file is opened `mode=ro&immutable=1`, which tells
 *    SQLite the bytes cannot change underneath it: no WAL, no shared-memory segment, no
 *    file locks, no rollback journal. Combined with `SQLITE_OPEN_NOMUTEX` there is no
 *    mutex anywhere on the read path -- not one this library takes, and not one SQLite
 *    takes for it.
 *  - **Any path.** A filesystem path or a full `file:` URI, wherever the caller keeps it.
 *    A process can read several databases at once with no interaction between them.
 *
 * @section threading The one rule
 *
 * `immutable=1` is a promise about the FILE, and `SQLITE_OPEN_NOMUTEX` is a promise about
 * the HANDLE. So:
 *
 *  - **One handle per thread.** A handle is cheap (@ref ls_db_open does no parsing; it
 *    opens a file). Threads must not share one -- there is no mutex to make that safe,
 *    which is the point. Handles to the same file from different threads are fine and do
 *    not interact.
 *  - **A database file is replaced, never edited.** The rebuild writes a new file and
 *    renames it into place, so an open handle keeps reading the inode it opened and a
 *    handle opened during the rename gets one whole version or the other. Editing a live
 *    database in place breaks the immutable promise and the reader will return wrong
 *    answers rather than fail -- do not do it. tools/import_lensfun_xml.c does the
 *    write-then-rename itself.
 */

#ifndef LENSSERIOUS_DB_H
#define LENSSERIOUS_DB_H

#include "lensserious.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief The comparison form of a name: ASCII case-folded, punctuation dropped, runs of
 * whitespace collapsed to one space.
 *
 * @details Exposed because the importer and the reader must produce identical bytes -- the
 * importer stores it, the reader compares against it -- so there is exactly one
 * implementation and no locale, ICU version or platform in the middle of it. Non-ASCII
 * bytes pass through untouched.
 *
 * @return the length written, excluding the terminator.
 */
size_t ls_db_normalize(const char *in, char *out, size_t out_size);

/** @brief The tokens ls_db_match_lens() compares, from an already-normalised name.
 *
 * @details Split on spaces, and again wherever a letter meets a digit -- normalisation
 * drops '-', so "16-35mm" would otherwise fuse into one meaningless token instead of
 * matching the catalogue's "16 35mm".
 *
 * Exposed for the same reason as ls_db_normalize(): the importer writes these into the
 * token index and the matcher looks them up, so there must be exactly one implementation.
 *
 * @param norm a name already put through ls_db_normalize().
 * @param out_tokens caller's array of @p max buffers, each at least @p stride bytes.
 * @param max how many tokens to write at most.
 * @param stride the size of one buffer in @p out_tokens, in bytes.
 * @return how many tokens were written.
 */
int ls_db_tokenize(const char *norm, char *out_tokens, int max, int stride);

/** @brief FNV-1a of a token, the hash the stored digest and the matcher both use. */
unsigned ls_db_token_hash(const char *token);

/**
 * @brief Pack a normalised name's tokens into the digest stored in `lens_name.tokens`.
 *
 * @details Layout: `uint16 n`, then n x `uint32` hash, then n x `uint8` length. The matcher
 * scores on those alone -- it never sees the token text -- which is what lets the candidate
 * query read a blob out of a covering index instead of decoding and re-tokenising text on
 * every lookup.
 *
 * @return bytes written, or 0 if @p out is too small.
 */
size_t ls_db_token_digest(const char *norm, unsigned char *out, size_t out_size);

/** An open database. Not shared between threads; see @ref threading. */
typedef struct ls_db_t ls_db_t;

/** A camera, as far as a correction is concerned. */
typedef struct ls_camera_t
{
  float crop_factor;
  /** The camera's mount, as an id into the database it was read from. Meaningless against
   *  any other database, and not stable across rebuilds. Pass it to ls_db_lens_fits_mount(). */
  long long mount_id;
} ls_camera_t;

/**
 * @brief Open a database for reading.
 *
 * @param path a filesystem path, or a `file:` URI. Either way it is opened read-only and
 * immutable; any query parameters the caller supplies in a URI are preserved, but `mode`
 * and `immutable` are forced.
 * @return NULL if the file cannot be opened or is not a LensSerious database of a schema
 * version this build understands. Never fails for concurrency reasons: there are none.
 */
ls_db_t *ls_db_open(const char *path);

/** @brief Release a handle. Safe on NULL. */
void ls_db_close(ls_db_t *db);

/** @brief The last error on @p db, as a string owned by @p db, or NULL. */
const char *ls_db_error(const ls_db_t *db);

/** @brief Schema version of the open file, or -1. Bumped when the layout changes. */
int ls_db_schema_version(const ls_db_t *db);

/** @brief A `meta` value by key (`built_utc`, `source`, `lensfun_db_version`, ...).
 *  @return bytes written excluding the terminator, or -1 if absent. */
int ls_db_meta(ls_db_t *db, const char *key, char *out, size_t out_size);

/**
 * @brief Find a camera by maker and model.
 *
 * @details Matching is exact on the normalised form of each name -- case-folded, with
 * punctuation and runs of whitespace collapsed -- against every spelling the database
 * holds for that camera, translations included. @p maker may be NULL to search on the
 * model alone.
 *
 * This is NOT lensfun's fuzzy scorer, which also tolerates missing words, reordered
 * tokens and vendor prefixes. Porting it is separate work; until then a caller that needs
 * that behaviour should treat a miss here as "not found" rather than as "no such camera".
 *
 * @return 1 on a match, 0 if nothing matched, -1 on error.
 */
int ls_db_find_camera(ls_db_t *db, const char *maker, const char *model, ls_camera_t *out);

/**
 * @brief Find a lens by maker and model, and fill @p out with its coefficients.
 *
 * @details Same matching as ls_db_find_camera(). When several lenses share a name --
 * upstream distinguishes them by the sensor they were calibrated on -- the one whose
 * crop factor is closest to @p crop wins, matching how a caller would pick by hand;
 * pass 0 to take the first.
 *
 * @p out is filled completely, calibration arrays included, and is thereafter independent
 * of @p db: the handle may be closed and the value stays valid. That is the whole shape of
 * this library -- a lens is data, not a handle into a database.
 *
 * @return 1 on a match, 0 if nothing matched, -1 on error (including a lens carrying more
 * than LS_MAX_CALIB entries of any kind, which is a database this build cannot represent).
 */
int ls_db_find_lens(ls_db_t *db, const char *maker, const char *model, float crop,
                    ls_lens_t *out);

/** One candidate from ls_db_match_lens(), best score first. */
typedef struct ls_db_match_t
{
  long long lens_id;
  float score;        /**< 0..100; 100 is an exact match of every token, both ways. */
} ls_db_match_t;

/**
 * @brief Find the lenses a free-text name most likely refers to.
 *
 * @details ls_db_find_lens() answers "is there a lens called exactly this". This answers
 * the question a raw file actually asks, where the EXIF string is a vendor's abbreviation
 * of the name upstream chose -- "16-35mm f/4G ED VR" against "Nikon AF-S Nikkor 16-35mm
 * f/4G ED VR", with tokens missing, reordered, or spelled differently.
 *
 * Scoring is token-based and deliberately simple; see the implementation for what each
 * weight is and why. It is calibrated against liblensfun's own decisions rather than
 * against a specification: tests/match_lensfun.c asks both this and lf_db_find_lenses_hd()
 * the same questions over the whole database and reports where they disagree.
 *
 * @param maker may be NULL. When given it is scored, not required -- vendors disagree with
 * upstream about their own name often enough that requiring it loses more than it saves.
 * @param mount_id when > 0, only lenses that fit this mount are considered
 * (ls_db_find_camera() supplies it). 0 considers every lens.
 * @param db an open database.
 * @param model the free-text name to resolve. Required.
 * @param out caller's array, @p max entries, filled best-first.
 * @param max the length of @p out.
 * @return how many candidates were written, or -1 on error.
 */
int ls_db_match_lens(ls_db_t *db, const char *maker, const char *model, long long mount_id,
                     ls_db_match_t *out, int max);

/** @brief Load a lens by its database id, for a caller that already resolved one. */
int ls_db_lens_by_id(ls_db_t *db, long long lens_id, ls_lens_t *out);

/**
 * @brief Does a lens fit a camera's mount, upstream's compatibility table included?
 * @return 1 yes, 0 no, -1 on error.
 */
int ls_db_lens_fits_mount(ls_db_t *db, long long lens_id, long long mount_id);

/**
 * @brief Enumerate lens ids, oldest-inserted first, for tests and for a GUI's lens picker.
 *
 * @param db an open database.
 * @param out_ids caller's array, or NULL to count only.
 * @param max the length of @p out_ids.
 * @return how many ids were written (or how many exist, if @p out_ids is NULL), or -1.
 */
int ls_db_list_lenses(ls_db_t *db, long long *out_ids, int max);

/** @brief The lens's maker/model, as stored. @return bytes written, or -1. */
int ls_db_lens_name(ls_db_t *db, long long lens_id, char *maker, size_t maker_size,
                    char *model, size_t model_size);

#ifdef __cplusplus
}
#endif

#endif /* LENSSERIOUS_DB_H */
