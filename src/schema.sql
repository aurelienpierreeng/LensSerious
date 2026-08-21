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

PRAGMA user_version = 1;

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
CREATE TABLE calib_distortion (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  model   INTEGER NOT NULL,           -- ls_dist_model_t
  focal   REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL
);

CREATE TABLE calib_tca (
  lens_id INTEGER NOT NULL REFERENCES lens(id),
  model   INTEGER NOT NULL,           -- ls_tca_model_t
  focal   REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL,
  t3 REAL NOT NULL, t4 REAL NOT NULL, t5 REAL NOT NULL
);

CREATE TABLE calib_vignetting (
  lens_id  INTEGER NOT NULL REFERENCES lens(id),
  model    INTEGER NOT NULL,          -- ls_vig_model_t
  focal    REAL NOT NULL,
  aperture REAL NOT NULL,
  distance REAL NOT NULL,
  t0 REAL NOT NULL, t1 REAL NOT NULL, t2 REAL NOT NULL
);

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
CREATE INDEX idx_camera_name_norm ON camera_name(norm);
CREATE INDEX idx_calib_distortion_lens ON calib_distortion(lens_id);
CREATE INDEX idx_calib_tca_lens ON calib_tca(lens_id);
CREATE INDEX idx_calib_vignetting_lens ON calib_vignetting(lens_id);
CREATE INDEX idx_lens_mount_mount ON lens_mount(mount_id);
