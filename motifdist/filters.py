"""
filters.py  — Stage B: the four artifact filters (the heart of the toolkit)
==========================================================================
Applied to Stage-A co-occurrence hits to drop artifacts. Each filter catches
something the others miss; keep all four.

  1. Small-number noise   (`filter_small_number`)
        The top raw hits were 40x "enriched" but built on 1-2 observations.
        Require both patterns to be reasonably large and the pair to have enough
        observed co-occurrences.

  2. Sibling-motif filter (`sibling_qvals` + `drop_siblings`)
        Two patterns can be the *same TF* captured twice; their co-occurrence is
        then trivial. Pattern-vs-pattern TOMTOM; drop pairs whose motifs match
        (q < 0.05). This correctly caught pattern_11 <-> pattern_16 (both ZNF143,
        q=2.7e-06). Uses real TOMTOM, not the notebook's unreliable homemade CWM
        cosine score.

  3. Promiscuity filter   (`filter_promiscuity`)
        A pattern that pairs with *everything* (Zic=11, YY2=10 partners in the
        reference run) marks a busy regulatory neighborhood, not a specific
        composite. Keep only pairs where BOTH partners are selective.

  4. Overlap-exclusion    (`overlap_spike_qc`; exclusion itself lives in spacing/
        cooccurrence via the *_no_overlap functions). NOVEL contribution.
        Cross-pattern seqlets can physically overlap (share up to ~49% of bases)
        and produce a phantom ~15 bp spacing spike that survives the sibling test.
        Exclude physically-overlapping seqlet pairs before any spacing analysis.
"""

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("motifdist.filters")


# --- Filter 1: small-number noise ------------------------------------------

def filter_small_number(cooc, min_seqlets=200, min_pairs=30, min_enrichment=1.5):
    """Keep pairs where both patterns have >= `min_seqlets` seqlets, the pair has
    >= `min_pairs` observed co-occurrences, and enrichment > `min_enrichment`.

    Expects the tidy DataFrame from cooccurrence_enrichment. Returns the surviving
    rows (a copy) and logs how many were dropped.
    """
    before = len(cooc)
    keep = cooc[
        (cooc["anchor_size"] >= min_seqlets)
        & (cooc["partner_size"] >= min_seqlets)
        & (cooc["obs_pairs"] >= min_pairs)
        & (cooc["enrichment"] > min_enrichment)
    ].copy()
    log.info("small-number filter: kept %d / %d pairs "
             "(min_seqlets=%d, min_pairs=%d, enrichment>%.1f)",
             len(keep), before, min_seqlets, min_pairs, min_enrichment)
    return keep


# --- Filter 2: sibling motifs (real TOMTOM) --------------------------------

def sibling_qvals(pairs, tomtom_qvals):
    """Attach a `motif_qval` column: the TOMTOM q-value between each pair's two
    patterns (low q = similar motifs = likely siblings).

    `tomtom_qvals` is the DataFrame from tomtom.pattern_vs_pattern_qvals (columns
    Query_ID, Target_ID, q-value). Missing pairs get NaN (treated as "different").
    """
    lut = {}
    for _, r in tomtom_qvals.iterrows():
        lut[(r["Query_ID"], r["Target_ID"])] = r["q-value"]

    def q(a, b):
        # TOMTOM is asymmetric in reporting; take the smaller of both directions.
        vals = [lut.get((a, b)), lut.get((b, a))]
        vals = [v for v in vals if v is not None and not pd.isna(v)]
        return min(vals) if vals else np.nan

    out = pairs.copy()
    out["motif_qval"] = [q(r["anchor"], r["partner"]) for _, r in out.iterrows()]
    return out


def drop_siblings(pairs_with_qval, q_thresh=0.05):
    """Drop pairs whose motifs match each other (motif_qval < `q_thresh`).

    Pairs with NaN q-value (not reported by TOMTOM, i.e. clearly different) are
    kept. Requires the `motif_qval` column from `sibling_qvals`.
    """
    before = len(pairs_with_qval)
    is_sibling = pairs_with_qval["motif_qval"] < q_thresh  # NaN -> False -> kept
    keep = pairs_with_qval[~is_sibling].copy()
    dropped = pairs_with_qval[is_sibling]
    log.info("sibling filter: kept %d / %d pairs (dropped %d as same-motif, q<%.3g)",
             len(keep), before, len(dropped), q_thresh)
    for _, r in dropped.iterrows():
        log.debug("  dropped sibling %s <-> %s (q=%.2e)",
                  r["anchor"], r["partner"], r["motif_qval"])
    return keep


# --- Filter 3: promiscuity --------------------------------------------------

def partner_counts(pairs):
    """How many distinct partners each pattern appears with among `pairs`.

    Counts unordered involvement: a pattern is counted once per distinct other
    pattern it pairs with (in either anchor or partner role). Returns a dict
    {pattern: n_partners}.
    """
    partners = {}
    for _, r in pairs.iterrows():
        a, b = r["anchor"], r["partner"]
        partners.setdefault(a, set()).add(b)
        partners.setdefault(b, set()).add(a)
    return {p: len(s) for p, s in partners.items()}


def filter_promiscuity(pairs, max_partners=2):
    """Keep only pairs where BOTH patterns have <= `max_partners` distinct
    partners among the surviving pairs. Adds anchor/partner promiscuity columns
    and logs the partner-count distribution.
    """
    counts = partner_counts(pairs)
    dist = {}
    for c in counts.values():
        dist[c] = dist.get(c, 0) + 1
    log.info("promiscuity: partner-count distribution %s",
             {k: dist[k] for k in sorted(dist)})

    out = pairs.copy()
    out["anchor_partners"] = out["anchor"].map(counts)
    out["partner_partners"] = out["partner"].map(counts)
    keep = out[(out["anchor_partners"] <= max_partners)
               & (out["partner_partners"] <= max_partners)].copy()
    log.info("promiscuity filter: kept %d / %d pairs (both partners <= %d partners)",
             len(keep), len(pairs), max_partners)
    return keep


# --- Filter 4: overlap-exclusion QC diagnostic ------------------------------
# The exclusion itself is applied by the *_no_overlap functions in spacing.py and
# cooccurrence.py. This is the QC report that shows the phantom spike is an
# overlap artifact.

def overlap_spike_qc(ann, patternA, patternB, spike_window=5, ceiling=100):
    """Diagnose the phantom tight-spacing spike: what fraction of anchor seqlets
    *in the spike* (|distance| <= `spike_window` bp) physically overlap their
    partner, versus those outside the spike (the background)?

    In the reference run the spike was ~100% overlapping pairs vs ~20% in the
    background — a clean tell that the spike is an overlap artifact, not spacing.
    Returns a dict with the two fractions and counts.
    """
    A = ann[ann["pattern"] == patternA]
    B = ann[ann["pattern"] == patternB]
    spike, background = [], []
    for chrom in set(A["chrom"]) & set(B["chrom"]):
        a = A[A["chrom"] == chrom]
        b = B[B["chrom"] == chrom]
        bm, bs, be = b["mid"].values, b["start"].values, b["end"].values
        if len(bm) == 0:
            continue
        for _, ar in a.iterrows():
            d_all = bm - ar["mid"]
            k = int(np.argmin(np.abs(d_all)))
            d = d_all[k]
            if abs(d) > ceiling:
                continue
            if ar["strand"] == "-":
                d = -d
            overlaps = bool((bs[k] < ar["end"]) and (be[k] > ar["start"]))
            (spike if abs(d) <= spike_window else background).append(overlaps)

    spike = np.array(spike, dtype=bool)
    background = np.array(background, dtype=bool)
    report = {
        "anchor": patternA, "partner": patternB,
        "n_spike": int(len(spike)),
        "spike_overlap_frac": round(float(spike.mean()), 3) if len(spike) else float("nan"),
        "n_background": int(len(background)),
        "background_overlap_frac": round(float(background.mean()), 3) if len(background) else float("nan"),
    }
    log.info("overlap QC %s -> %s: spike %.0f%% overlap (n=%d) vs background %.0f%% (n=%d)",
             patternA, patternB, 100 * report["spike_overlap_frac"], report["n_spike"],
             100 * report["background_overlap_frac"], report["n_background"])
    return report
