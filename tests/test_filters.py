import numpy as np
import pandas as pd

from conftest import make_ann
from motifdist.filters import (
    filter_small_number, sibling_qvals, drop_siblings,
    partner_counts, filter_promiscuity, overlap_spike_qc,
)


def _cooc(rows):
    cols = ["anchor", "partner", "anchor_size", "partner_size", "obs_pairs", "enrichment"]
    return pd.DataFrame(rows, columns=cols)


def test_small_number_filter():
    cooc = _cooc([
        ("A", "B", 300, 300, 50, 3.0),   # keep
        ("C", "D", 300, 300, 2, 40.0),   # drop: too few obs pairs (noise)
        ("E", "F", 10, 300, 50, 3.0),    # drop: E too small
        ("G", "H", 300, 300, 50, 1.1),   # drop: not enriched
    ])
    keep = filter_small_number(cooc, min_seqlets=200, min_pairs=30, min_enrichment=1.5)
    assert set(zip(keep["anchor"], keep["partner"])) == {("A", "B")}


def test_sibling_filter_drops_matching_motifs():
    pairs = _cooc([
        ("A", "B", 300, 300, 50, 3.0),
        ("C", "D", 300, 300, 50, 3.0),
    ])
    # TOMTOM says A<->B are the same motif (q tiny); C/D not reported (different).
    qvals = pd.DataFrame({
        "Query_ID": ["A"], "Target_ID": ["B"], "q-value": [2.7e-06],
    })
    with_q = sibling_qvals(pairs, qvals)
    assert with_q.loc[with_q["anchor"] == "A", "motif_qval"].iloc[0] == 2.7e-06
    assert np.isnan(with_q.loc[with_q["anchor"] == "C", "motif_qval"].iloc[0])
    kept = drop_siblings(with_q, q_thresh=0.05)
    assert set(zip(kept["anchor"], kept["partner"])) == {("C", "D")}


def test_sibling_qval_symmetric_lookup():
    # TOMTOM reported only the B->A direction; lookup must still find it for A->B.
    pairs = _cooc([("A", "B", 300, 300, 50, 3.0)])
    qvals = pd.DataFrame({"Query_ID": ["B"], "Target_ID": ["A"], "q-value": [1e-8]})
    with_q = sibling_qvals(pairs, qvals)
    assert with_q["motif_qval"].iloc[0] == 1e-8


def test_promiscuity_filter():
    # A pairs with B, C, D (3 partners = promiscuous). B pairs only with A.
    pairs = _cooc([
        ("A", "B", 300, 300, 50, 3.0),
        ("A", "C", 300, 300, 50, 3.0),
        ("A", "D", 300, 300, 50, 3.0),
        ("E", "F", 300, 300, 50, 3.0),   # both selective -> survive
    ])
    counts = partner_counts(pairs)
    assert counts["A"] == 3 and counts["B"] == 1 and counts["E"] == 1
    keep = filter_promiscuity(pairs, max_partners=2)
    # every pair touching A is dropped (A has 3 partners); only E-F survives
    assert set(zip(keep["anchor"], keep["partner"])) == {("E", "F")}


def test_overlap_spike_qc_flags_artifact():
    # Two separated regions. In the "spike" region every anchor's nearest partner
    # overlaps it (mids 10bp apart). In the "background" region the nearest partner
    # is 90bp away and does not overlap. The QC should show the spike is ~100%
    # overlapping pairs while the background is ~0%.
    rows = []
    for i in range(15):
        s = 1000 * i
        rows.append(("chr1", s, s + 40, "+", "A"))        # mid s+20
        rows.append(("chr1", s + 10, s + 50, "+", "B"))   # mid s+30, overlaps A, d=10 -> spike
    for i in range(15):
        s = 500000 + 1000 * i
        rows.append(("chr1", s, s + 40, "+", "A"))        # mid s+20
        rows.append(("chr1", s + 90, s + 130, "+", "B"))  # mid s+110, d=90, no overlap -> background
    ann = make_ann(rows)
    report = overlap_spike_qc(ann, "A", "B", spike_window=10, ceiling=200)
    assert report["n_spike"] > 0 and report["n_background"] > 0
    assert report["spike_overlap_frac"] >= 0.9         # spike is an overlap artifact
    assert report["background_overlap_frac"] <= 0.1    # background is not
