"""
cooccurrence.py  — Stage A: the co-occurrence screen + null model #1
====================================================================
For every ordered pattern pair (A, B), how often does an A-seqlet have a
B-seqlet within `ceiling` bp on the same chromosome? Compare that observed count
against a null that shuffles *which pattern label* sits at each existing seqlet
position, within each chromosome.

Why that null (null model #1, the "co-occurrence null"): every seqlet already
lives inside an accessible peak, so any two patterns are automatically somewhat
close just from peak geometry. Shuffling labels among the existing positions
preserves that geometry and destroys only pattern-*specific* placement. So
enrichment = observed / null-mean isolates real pattern affinity from "they both
sit in peaks."

Expected, documented behaviors (don't treat as bugs):
  * The median enrichment across all pairs should sit near 1 (~0.81 in the
    reference run). We report it as a sanity check — far from 1 means the null is
    biased.
  * Most pairs come out *depleted*, not enriched: patterns compete for limited
    space inside peaks.
"""

import logging

import numpy as np
import pandas as pd

log = logging.getLogger("motifdist.cooccurrence")


def count_near_pairs(dfA, dfB, ceiling):
    """Count how many dfA seqlets ("anchors") have >=1 dfB seqlet within `ceiling`
    bp on the same chromosome. Directional: A is the anchor.

    Both frames must carry a `mid` column (seqlet midpoint). Vectorized per
    chromosome via binary search: for each anchor we only check its nearest
    neighbor to the left and right in the sorted partner positions.
    """
    total = 0
    for chrom in set(dfA["chrom"]) & set(dfB["chrom"]):
        a = np.sort(dfA.loc[dfA["chrom"] == chrom, "mid"].values)
        b = np.sort(dfB.loc[dfB["chrom"] == chrom, "mid"].values)
        if len(a) == 0 or len(b) == 0:
            continue
        idx = np.searchsorted(b, a)
        right = np.where(idx < len(b), np.abs(b[np.clip(idx, 0, len(b) - 1)] - a), np.inf)
        left = np.where(idx > 0, np.abs(a - b[np.clip(idx - 1, 0, len(b) - 1)]), np.inf)
        total += int((np.minimum(left, right) <= ceiling).sum())
    return total


def count_near_pairs_no_overlap(dfA, dfB, ceiling):
    """Like `count_near_pairs`, but a partner that physically *overlaps* the
    anchor on the genome does not count as "near".

    Two seqlets in different patterns can share up to ~49% of their bases (MoDISco
    only suppresses >50% overlap within a single extraction pass, not across
    patterns). Counting those as co-occurrences inflates pairs that merely detect
    the same DNA. See the overlap-exclusion filter and TOOLKIT_SPEC.md §3.B.4.
    """
    total = 0
    for chrom in set(dfA["chrom"]) & set(dfB["chrom"]):
        a = dfA[dfA["chrom"] == chrom]
        b = dfB[dfB["chrom"] == chrom]
        bm, bs, be = b["mid"].values, b["start"].values, b["end"].values
        if len(bm) == 0:
            continue
        for _, ar in a.iterrows():
            ok = ~((bs < ar["end"]) & (be > ar["start"]))  # non-overlapping partners
            if not ok.any():
                continue
            if np.min(np.abs(bm[ok] - ar["mid"])) <= ceiling:
                total += 1
    return total


def shuffled_copy(ann_df, rng):
    """Null model #1: return a copy of `ann_df` with the `mid` of each seqlet
    randomly reassigned *within its chromosome*, drawn from the pool of all seqlet
    midpoints on that chromosome.

    This preserves where seqlets can be (peak geometry) and destroys only which
    pattern sits where. `rng` is a numpy Generator, so results are reproducible
    for a fixed seed.

    start/end are rebuilt around the shuffled midpoint (keeping each seqlet's
    width) so the copy stays internally consistent — this matters when the
    overlap-aware counter (`count_near_pairs_no_overlap`) reads them. The default
    counter uses only `mid`, so its results are unaffected.
    """
    out = ann_df.copy()
    for _, g in ann_df.groupby("chrom"):
        out.loc[g.index, "mid"] = rng.permutation(g["mid"].values)
    half = (out["end"] - out["start"]) // 2
    out["start"] = out["mid"] - half
    out["end"] = out["mid"] + half
    return out


def cooccurrence_enrichment(ann, ceiling=1000, n_shuffle=20, seed=0, counter=count_near_pairs):
    """Run the full Stage-A screen and return a tidy DataFrame, one row per
    ordered pair (anchor, partner).

    Columns: anchor, partner, anchor_size, partner_size, obs_pairs, null_mean,
    enrichment (= obs / null_mean), log2_enrich, pvalue (two-sided empirical).

    The empirical p-value uses the +1 correction, so with `n_shuffle` shuffles the
    smallest achievable p is 1/(n_shuffle+1). Use enough shuffles that this floor
    is below the significance you care about (>=200 for reported results).

    `counter` lets you swap in `count_near_pairs_no_overlap` to run the whole
    screen with overlapping seqlet pairs excluded.
    """
    patterns = sorted(ann["pattern"].unique())
    n = len(patterns)
    pat_dfs = {p: ann[ann["pattern"] == p] for p in patterns}
    sizes = np.array([len(pat_dfs[p]) for p in patterns])

    obs = np.zeros((n, n))
    for i, pA in enumerate(patterns):
        for j, pB in enumerate(patterns):
            if i != j:
                obs[i, j] = counter(pat_dfs[pA], pat_dfs[pB], ceiling)

    rng = np.random.default_rng(seed)
    null_stack = np.zeros((n_shuffle, n, n))
    for s in range(n_shuffle):
        shuf = shuffled_copy(ann, rng)
        shuf_dfs = {p: shuf[shuf["pattern"] == p] for p in patterns}
        for i, pA in enumerate(patterns):
            for j, pB in enumerate(patterns):
                if i != j:
                    null_stack[s, i, j] = counter(shuf_dfs[pA], shuf_dfs[pB], ceiling)
        log.info("co-occurrence null shuffle %d/%d", s + 1, n_shuffle)

    null_mean = null_stack.mean(axis=0)
    enrichment = np.divide(obs, null_mean, out=np.zeros_like(obs), where=null_mean > 0)

    rows = []
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            null_vals = null_stack[:, i, j]
            mu = null_vals.mean()
            # two-sided empirical p: how often is the shuffle as far from its mean
            # as the observed value is?
            as_extreme = int(np.sum(np.abs(null_vals - mu) >= abs(obs[i, j] - mu)))
            pval = (as_extreme + 1) / (n_shuffle + 1)
            e = enrichment[i, j]
            rows.append({
                "anchor": patterns[i], "partner": patterns[j],
                "anchor_size": int(sizes[i]), "partner_size": int(sizes[j]),
                "obs_pairs": int(obs[i, j]),
                "null_mean": round(float(null_mean[i, j]), 2),
                "enrichment": round(float(e), 3),
                "log2_enrich": round(float(np.log2(e)), 3) if e > 0 else np.nan,
                "pvalue": round(pval, 4),
            })
    cooc = pd.DataFrame(rows)

    # Sanity check to surface, not hide (TOOLKIT_SPEC.md §3.A).
    off = cooc.loc[cooc["enrichment"] > 0, "enrichment"]
    med = float(off.median()) if len(off) else float("nan")
    frac_depleted = float((off < 0.67).mean()) if len(off) else float("nan")
    log.info("median enrichment across pairs = %.2f (should sit near 1)", med)
    log.info("fraction of pairs depleted (<0.67x) = %.0f%% "
             "(most pairs depleted is expected: patterns compete for peak space)",
             100 * frac_depleted)

    return cooc
