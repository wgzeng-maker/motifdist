"""
nulls.py  — null model #2: is a spacing peak real, or just lumpiness?
====================================================================
This is a *different* null from the co-occurrence null in cooccurrence.py. Do not
conflate them.

A tall histogram bin is cheap. Throw ~600 distances into ~100 bins and the tallest
bin will be several times the typical bin *by luck alone* — randomness is lumpy.
So a "peak" in a spacing histogram is not evidence of preferred spacing until you
show it beats what chance produces.

Method: `peak_height_ratio` measures how peaky an observed histogram is (tallest
bin / median non-empty bin). The null shuffles seqlet positions within each
chromosome, rebuilds the spacing histogram, and records how tall the tallest bin
gets by chance. The empirical p-value is the fraction of shuffles whose peak is at
least as tall as the observed one.

Reference cautionary tale: the NFI->Zic "peak" (4.6x) was exactly what chance
produces (~4.7x). Without this null you would report a spacing preference that
isn't there.
"""

import logging

import numpy as np
import pandas as pd

from .spacing import spacing_analysis
# The two nulls shuffle positions the same way (permute midpoints within each
# chromosome, rebuild start/end); they differ only in the *statistic* they then
# compute — co-occurrence count (null #1) vs spacing peak height (null #2). Share
# the one shuffle so they cannot drift.
from .cooccurrence import shuffled_copy as _shuffle_positions

log = logging.getLogger("motifdist.nulls")

DEFAULT_BINS = np.arange(-250, 251, 5)


def peak_height_ratio(distances, bins=DEFAULT_BINS):
    """Tallest bin / median non-empty bin height. NaN if < 20 distances.

    This is the single "peakiness" statistic used for the observed data and every
    shuffle, so observed and null are directly comparable.
    """
    if len(distances) < 20:
        return float("nan")
    h, _ = np.histogram(distances, bins=bins)
    nz = h[h > 0]
    bg = np.median(nz) if len(nz) else 0
    return float(h.max() / bg) if bg > 0 else float("nan")


def spacing_null_test(ann, pairs, ceiling=250, n_shuffles=200, seed=0,
                      bins=DEFAULT_BINS, min_distances=20, exclude_overlap=True):
    """For each (anchor, partner) pair, test the same- and opposite-strand spacing
    peaks against the shuffle null. Returns a tidy DataFrame.

    Columns: anchor, partner, strand, n, obs_peak, null_mean, null_max,
    n_valid_nulls, pvalue.

    `n_shuffles` defaults to 200 (the recommended minimum for a reported result);
    drop it to ~30 only for a quick exploratory pass. Strands with fewer than
    `min_distances` observed distances are skipped (reported with n but no
    p-value).
    """
    rng = np.random.default_rng(seed)
    results = []

    for anchor, partner in pairs:
        same_obs, opp_obs = spacing_analysis(ann, anchor, partner, ceiling=ceiling,
                                             exclude_overlap=exclude_overlap)
        for strand_label, obs in (("same", same_obs), ("opposite", opp_obs)):
            if len(obs) < min_distances:
                results.append({
                    "anchor": anchor, "partner": partner, "strand": strand_label,
                    "n": int(len(obs)), "obs_peak": np.nan, "null_mean": np.nan,
                    "null_max": np.nan, "n_valid_nulls": 0, "pvalue": np.nan,
                })
                log.info("%s -> %s [%s]: n=%d, too few to test",
                         anchor, partner, strand_label, len(obs))
                continue

            obs_ratio = peak_height_ratio(obs, bins=bins)
            null_ratios = []
            for _ in range(n_shuffles):
                shuf = _shuffle_positions(ann, rng)
                sm, op = spacing_analysis(shuf, anchor, partner, ceiling=ceiling,
                                          exclude_overlap=exclude_overlap)
                d = sm if strand_label == "same" else op
                r = peak_height_ratio(d, bins=bins)
                if not np.isnan(r):
                    null_ratios.append(r)

            null_ratios = np.array(null_ratios)
            if len(null_ratios) == 0:
                # every shuffle produced too few distances to score
                results.append({
                    "anchor": anchor, "partner": partner, "strand": strand_label,
                    "n": int(len(obs)), "obs_peak": round(obs_ratio, 2),
                    "null_mean": np.nan, "null_max": np.nan,
                    "n_valid_nulls": 0, "pvalue": np.nan,
                })
                log.info("%s -> %s [%s]: null empty after shuffling, skipped",
                         anchor, partner, strand_label)
                continue

            pval = (int(np.sum(null_ratios >= obs_ratio)) + 1) / (len(null_ratios) + 1)
            results.append({
                "anchor": anchor, "partner": partner, "strand": strand_label,
                "n": int(len(obs)),
                "obs_peak": round(obs_ratio, 2),
                "null_mean": round(float(null_ratios.mean()), 2),
                "null_max": round(float(null_ratios.max()), 2),
                "n_valid_nulls": int(len(null_ratios)),
                "pvalue": round(pval, 4),
            })
            log.info("%s -> %s [%s]: obs=%.1fx null=%.1fx (max %.1f) p=%.3f [%d nulls]",
                     anchor, partner, strand_label, obs_ratio, null_ratios.mean(),
                     null_ratios.max(), pval, len(null_ratios))

    return pd.DataFrame(results)
