"""
io.py
=====
Load and validate the seqlet annotation table, and load pattern motif matrices
(PPM / CWM) from a MoDISco HDF5.

The whole toolkit trusts this table, so validation is strict and **fails loudly**
rather than continuing silently on bad input (the user's stated preference).
"""

import logging

import numpy as np
import pandas as pd

from . import REQUIRED_COLUMNS, pattern_to_h5

# hdf5plugin MUST be imported before h5py, or many production MoDISco .h5 files
# (which use compressed HDF5 filters) will not read. See TOOLKIT_SPEC.md §4.5.
try:
    import hdf5plugin  # noqa: F401
    _HDF5PLUGIN = True
except ModuleNotFoundError:
    _HDF5PLUGIN = False

log = logging.getLogger("motifdist.io")

# Columns that must never contain NaN. `subpattern` is required to be present but
# is allowed to be empty, because most downstream steps do not use it.
_NONNULL_COLUMNS = ("chrom", "start", "end", "strand", "pattern")


def load_annotation(path, add_mid=True):
    """Read and validate a seqlet annotation CSV; return a pandas DataFrame.

    Required columns: chrom, start, end (forward-strand BED coords, start < end),
    strand (+/-), pattern, subpattern. Validation is strict: any violation raises
    ValueError with a message naming the problem. On success, logs the per-pattern
    seqlet counts and (by default) adds a `mid` column = (start + end) // 2, which
    every spatial step uses as the seqlet's position.
    """
    df = pd.read_csv(path)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: missing required column(s): {missing}. "
            f"Expected all of {list(REQUIRED_COLUMNS)}."
        )

    for col in _NONNULL_COLUMNS:
        n_na = df[col].isna().sum()
        if n_na:
            raise ValueError(f"{path}: column {col!r} has {n_na} NaN/empty value(s).")

    for col in ("start", "end"):
        if not np.issubdtype(df[col].dtype, np.integer):
            # allow integer-valued floats, but reject genuinely fractional coords
            as_float = pd.to_numeric(df[col], errors="coerce")
            if as_float.isna().any() or not np.all(as_float == as_float.astype("int64")):
                raise ValueError(f"{path}: column {col!r} must be integer coordinates.")
            df[col] = as_float.astype("int64")

    bad_span = df["start"] >= df["end"]
    if bad_span.any():
        i = int(bad_span.idxmax())
        raise ValueError(
            f"{path}: {int(bad_span.sum())} row(s) have start >= end "
            f"(first at row {i}: start={df.at[i, 'start']}, end={df.at[i, 'end']})."
        )

    bad_strand = ~df["strand"].isin(["+", "-"])
    if bad_strand.any():
        vals = sorted(df.loc[bad_strand, "strand"].unique())
        raise ValueError(f"{path}: strand column has values outside {{+,-}}: {vals}.")

    if add_mid:
        df["mid"] = (df["start"] + df["end"]) // 2

    counts = df["pattern"].value_counts().sort_index()
    log.info("loaded %d seqlets across %d patterns from %s", len(df), len(counts), path)
    for pat, n in counts.items():
        log.debug("  %-20s %d seqlets", pat, n)

    return df


def open_h5(path):
    """Open a MoDISco HDF5 read-only, with a clear message if a compression
    filter is missing. Returns an h5py.File the caller must close (use as a
    context manager)."""
    import h5py

    try:
        return h5py.File(path, "r")
    except OSError as exc:
        msg = str(exc).lower()
        if ("filter" in msg or "plugin" in msg or "blosc" in msg) and not _HDF5PLUGIN:
            raise OSError(
                f"Could not open {path!r}. It likely uses HDF5 compression filters "
                "that need the optional package 'hdf5plugin'. Install it in this "
                "environment and rerun."
            ) from exc
        raise


def pattern_ppm(modisco_h5, pattern):
    """Position-probability matrix (base frequencies, shape (W, 4)) for a pattern.

    `pattern` is a table label like 'pos/pattern_0'. Stored under the MoDISco key
    'sequence'.
    """
    arm, pat = pattern_to_h5(pattern)
    with open_h5(modisco_h5) as f:
        return f[arm][pat]["sequence"][:]


def pattern_cwm(modisco_h5, pattern):
    """Contribution-weighted matrix (the pattern's average motif, shape (W, 4)).

    Stored under the MoDISco key 'contrib_scores'. This is signed attribution,
    not a probability matrix.
    """
    arm, pat = pattern_to_h5(pattern)
    with open_h5(modisco_h5) as f:
        return f[arm][pat]["contrib_scores"][:]


def list_patterns(modisco_h5):
    """All pattern table-labels present in the HDF5, e.g. ['pos/pattern_0', ...]."""
    out = []
    with open_h5(modisco_h5) as f:
        for arm in ("pos_patterns", "neg_patterns"):
            if arm not in f:
                continue
            short_arm = arm.split("_")[0]  # 'pos' / 'neg'
            for pat in f[arm]:
                out.append(f"{short_arm}/{pat}")
    return out
