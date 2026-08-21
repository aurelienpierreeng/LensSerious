-- LensSerious — the lens database as a table of numbers.
--
-- Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
--
-- Written once by tools/import_lensfun_xml.c from the upstream XML, then only ever read.
-- The reader opens it with immutable=1, which is what lets it skip file locking entirely
-- (see include/lensserious_db.h): a file this schema describes is never updated in place.
-- A rebuild writes a new file and renames it over the old one.
--
-- Calibration rows store the model as the LS_DIST_*/LS_TCA_*/LS_VIG_* integer and the
-- terms in the order ls_calib_*_t declares them, so a lookup is a copy rather than a
-- translation. That is deliberate: every conversion between the XML's vocabulary and the
-- evaluator's happens once, offline, where it can be checked, rather than per render.

PRAGMA user_version = 2;

CREATE TABLE meta (
  key   TEXT PRIMARY KEY,
  value TEXT NOT NULL
) WITHOUT ROWID;

CREATE TABLE mount (
  id   INTEGER PRIMARY KEY,
  name TEXT NOT NULL UNIQUE
);

-- Directed: upstream lists compatibility per mount and does not assume symmetry.
CREATE TABLE mount_compat (
  mount_id  INTEGER NOT NULL REFERENCES mount(id),
  compat_id INTEGER NOT NULL REFERENCES mount(id),
  PRIMARY KEY (mount_id, compat_id)
) WITHOUT ROWID;

CREATE TABLE camera (
  id          INTEGER PRIMARY KEY,
  maker       TEXT NOT NULL,          -- the untranslated <maker>, as matched against EXIF
  model       TEXT NOT NULL,
  variant     TEXT,
  mount_id    INTEGER REFERENCES mount(id),
  crop_factor REAL NOT NULL
);

CREATE TABLE lens (
  id           INTEGER PRIMARY KEY,
  maker        TEXT NOT NULL,
  model        TEXT NOT NULL,
  type         INTEGER NOT NULL,      -- ls_lens_type_t
  crop_factor  REAL NOT NULL,         -- of the CALIBRATION sensor
  aspect_ratio REAL NOT NULL,
  center_x     REAL NOT NULL,
  center_y     REAL NOT NULL,
  min_focal    REAL NOT NULL,         -- the vignetting IDW metric needs both, and most
  max_focal    REAL NOT NULL,         -- lenses only carry them in their model string
  min_aperture REAL,
  max_aperture REAL
);

CREATE TABLE lens_mount (
  lens_id  INTEGER NOT NULL REFERENCES lens(id),
  mount_id INTEGER NOT NULL REFERENCES mount(id),
  PRIMARY KEY (lens_id, mount_id)
) WITHOUT ROWID;

-- Every spelling a lens or camera answers to: the untranslated name plus each <... lang="">
-- variant. The matcher needs all of them, and `norm` is the case-folded, punctuation-free
-- form it compares on, computed once here rather than per lookup.
CREATE TABLE lens_name (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  kind    TEXT NOT NULL,              -- 'maker' | 'model'
  lang    TEXT,                       -- NULL for the untranslated name
  value   TEXT NOT NULL,
  norm    TEXT NOT NULL
);

CREATE TABLE camera_name (
  camera_id INTEGER NOT NULL REFERENCES camera(id),
  kind      TEXT NOT NULL,
  lang      TEXT,
  value     TEXT NOT NULL,
  norm      TEXT NOT NULL
);

-- Terms are packed exactly as ls_calib_dist_t/ls_calib_tca_t/ls_calib_vig_t declare them.
--
-- WITHOUT ROWID, keyed on (lens_id, ord). The table IS its primary-key b-tree, so one
-- lens's calibration is physically contiguous and reading it never leaves that b-tree --
-- where a plain table plus an index on lens_id costs a separate table seek for every row
-- the index finds.
--
-- `ord` is the position within the lens as upstream listed it, and it turns a convention
-- into a guarantee: the reader fills ls_lens_t's arrays in row order, which under a rowid
-- table was insertion order only by habit, and is now the key order the storage engine
-- promises.
--
-- Measured on a 70-row fetch, all three against the same surrounding schema:
--   WITHOUT ROWID (lens_id, ord)     0.0385 ms   3.76 MB   <- this
--   plain table + index(lens_id)     0.0436 ms   4.19 MB
--   plain table + index(lens_id,ord) 0.0474 ms   4.35 MB   -- the column AND a wider index
-- Faster and smaller: such a table repeats its key in every row, but dropping the three
-- lens_id indexes more than pays for that. An earlier round rejected this form on size,
-- having compared it against a baseline taken before lens_token and token_df existed --
-- which is what a comparison between two different schemas is worth.
CREATE TABLE calib_distortion (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  ord     INTEGER NOT NULL,           -- position within this lens, as upstream listed it
  model   INTEGER NOT NULL,           -- ls_dist_model_t
  focal   REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL
,
  PRIMARY KEY (lens_id, ord)
) WITHOUT ROWID;

CREATE TABLE calib_tca (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  ord     INTEGER NOT NULL,           -- position within this lens, as upstream listed it
  model   INTEGER NOT NULL,           -- ls_tca_model_t
  focal   REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL,
  t3 REAL NOT NULL, t4 REAL NOT NULL, t5 REAL NOT NULL
,
  PRIMARY KEY (lens_id, ord)
) WITHOUT ROWID;

CREATE TABLE calib_vignetting (
  lens_id  INTEGER NOT NULL REFERENCES lens(id),
  ord      INTEGER NOT NULL,          -- position within this lens, as upstream listed it
  model    INTEGER NOT NULL,          -- ls_vig_model_t
  focal    REAL NOT NULL,
  aperture REAL NOT NULL,
  distance REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL
,
  PRIMARY KEY (lens_id, ord)
) WITHOUT ROWID;

-- The inverted index the fuzzy matcher prunes with. One row per (name, token).
--
-- It is NOT used to gather every candidate that shares any token -- that was tried and
-- measured, and it is no better than a full scan, because the common tokens ("mm", "f",
-- "ed", "vr", the focal digits) each appear in a large fraction of the catalogue. What the
-- matcher does instead is ask this table which of the QUERY's tokens is rarest, and gather
-- only the lenses carrying that one. See ls_db_match_lens().
CREATE TABLE lens_token (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  kind    TEXT NOT NULL,              -- 'maker' | 'model'
  token   TEXT NOT NULL
);

-- How many lenses each token appears in, so the matcher can find its query's rarest token
-- with one point lookup per token instead of counting index rows. Counting is what made the
-- first version of that lookup slow: "mm" appears in thousands of names, and COUNT(*) has
-- to walk every one of them.
CREATE TABLE token_df (
  kind  TEXT NOT NULL,                -- 'maker' | 'model'
  token TEXT NOT NULL,
  df    INTEGER NOT NULL,
  PRIMARY KEY (kind, token)
) WITHOUT ROWID;

-- A focal-length remap a few compacts carry, kept so nothing from upstream is lost even
-- though the evaluators do not consult it yet.
CREATE TABLE lens_real_focal (
  lens_id    INTEGER NOT NULL REFERENCES lens(id),
  focal      REAL NOT NULL,
  real_focal REAL NOT NULL
);

-- Lookups are by name, and the calibration fetch is by lens_id. Nothing else is indexed:
-- an index the reader never uses is bytes every reader pays to mmap.
CREATE INDEX idx_lens_name_norm ON lens_name(norm);
-- Every name of one lens, for the fuzzy matcher's second phase: without it the query that
-- gathers the pruned candidates' names has no way in and scans the whole table, which is
-- the scan the pruning existed to avoid.
--
-- COVERING: kind and norm are in the index, not just lens_id, so the search never touches
-- the table. An index that only carries the search key still costs one table seek per row
-- found, which is what `SEARCH ... USING INDEX` hides -- the row it lands on is an index
-- entry, and every column outside the index is another b-tree descent.
CREATE INDEX idx_lens_name_lens ON lens_name(lens_id, kind, norm);
-- COVERING for the same reason: the pruning subquery wants lens_id, so carrying it in the
-- index means the token lookup never descends into lens_token itself.
CREATE INDEX idx_lens_token ON lens_token(kind, token, lens_id);
CREATE INDEX idx_camera_name_norm ON camera_name(norm);
CREATE INDEX idx_lens_mount_mount ON lens_mount(mount_id);
