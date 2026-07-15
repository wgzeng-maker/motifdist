"""
spacing.py  — Stage C: spacing / distance geometry
==================================================
For a pattern pair that survived the filters, measure the geometry: for each
anchor seqlet, the signed distance to its nearest partner, split by relative
strand.

Signed distance is measured in the **anchor's reading frame**:

    signed_dist = partner_mid - anchor_mid            # + = partner downstream on genome
    if anchor.strand == '-':  signed_dist *= -1       # read the anchor "forward"

Only the anchor is flipped. Relative strand ("same" vs "opposite") is decided by
strand *equality*, which is frame-independent, so the partner is never flipped.
Positive = partner downstream in the anchor's reading direction; negative =
upstream.

Overlap exclusion is ON by default: a partner seqlet that physically overlaps the
anchor is skipped before choosing the nearest one. Overlapping cross-pattern
seqlets (they can share up to ~49% of their bases) create a phantom ~15 bp
spacing spike that survives the sibling test; excluding them removes it. Set
`exclude_overlap=False` to reproduce the raw (inflated) view for comparison.
"""

import numpy as np


def spacing_analysis(ann, patternA, patternB, ceiling=1000, exclude_overlap=True):
    """Signed anchor-frame distances from each patternA seqlet to its nearest
    patternB seqlet within `ceiling` bp, on the same chromosome.

    Returns (same_strand, opposite_strand): two numpy arrays of signed distances.
    With `exclude_overlap=True` (default), partner seqlets that overlap the anchor
    on the genome are excluded before the nearest one is chosen.
    """
    A = ann[ann["pattern"] == patternA]
    B = ann[ann["pattern"] == patternB]
    same, opp = [], []

    for chrom in set(A["chrom"]) & set(B["chrom"]):
        a = A[A["chrom"] == chrom]
        b = B[B["chrom"] == chrom]
        bm, bs, be = b["mid"].values, b["start"].values, b["end"].values
        bstr = b["strand"].values
        if len(bm) == 0:
            continue
        for _, ar in a.iterrows():
            if exclude_overlap:
                keep = ~((bs < ar["end"]) & (be > ar["start"]))
            else:
                keep = np.ones(len(bm), dtype=bool)
            if not keep.any():
                continue
            d_all = bm[keep] - ar["mid"]          # + = partner downstream on genome
            k = int(np.argmin(np.abs(d_all)))
            d = d_all[k]
            if abs(d) > ceiling:
                continue
            if ar["strand"] == "-":
                d = -d                            # put into the anchor's reading frame
            partner_strand = bstr[keep][k]
            (same if ar["strand"] == partner_strand else opp).append(int(d))

    return np.array(same), np.array(opp)


def spacing_to_edges(ann, patternA, patternB, ceiling=100, exclude_overlap=False):
    """Diagnostic: distance from the anchor midpoint to the nearest *edge* (start
    or end) of the partner seqlet, in the anchor's reading frame.

    A bipartite motif split into two patterns has a gap in its middle, so a
    midpoint-to-midpoint distance can land in that empty gap and read as a fake
    offset. Measuring to the nearest edge is a cross-check on whether a spike is
    real spacing or a midpoint artifact. Returns one signed-distance array (not
    split by strand).
    """
    A = ann[ann["pattern"] == patternA]
    B = ann[ann["pattern"] == patternB]
    dists = []
    for chrom in set(A["chrom"]) & set(B["chrom"]):
        a = A[A["chrom"] == chrom]
        b = B[B["chrom"] == chrom]
        bs, be = b["start"].values, b["end"].values
        for _, ar in a.iterrows():
            am = ar["mid"]
            keep = np.ones(len(bs), dtype=bool)
            if exclude_overlap:
                keep = ~((bs < ar["end"]) & (be > ar["start"]))
            if not keep.any():
                continue
            d_start = bs[keep] - am
            d_end = be[keep] - am
            closest_edge = np.where(np.abs(d_start) < np.abs(d_end), d_start, d_end)
            k = int(np.argmin(np.abs(closest_edge)))
            d = closest_edge[k]
            if abs(d) > ceiling:
                continue
            if ar["strand"] == "-":
                d = -d
            dists.append(int(d))
    return np.array(dists)


def peak_stats(distances, window=25, bins=None):
    """Summarize a signed-distance array: (concentration, std, peak_location_bp).

    concentration = fraction of distances within +/-`window` bp of the modal bin
    (0 = flat/billboard, high = tight spike). peak_location = center of the most
    populated bin. Returns (nan, nan, nan) when there are fewer than 10 distances.
    """
    if len(distances) < 10:
        return float("nan"), float("nan"), float("nan")
    if bins is None:
        bins = np.arange(-1000, 1001, 50)
    hist, edges = np.histogram(distances, bins=bins)
    peak_bin = int(np.argmax(hist))
    peak_loc = (edges[peak_bin] + edges[peak_bin + 1]) / 2
    within = int(np.sum(np.abs(distances - peak_loc) <= window))
    return round(within / len(distances), 3), round(float(np.std(distances)), 1), round(float(peak_loc), 0)
