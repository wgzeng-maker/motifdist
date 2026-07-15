import numpy as np

from conftest import make_ann
from motifdist.spacing import spacing_analysis, peak_stats


def test_minus_strand_frame_flip():
    # Anchor on '-' at mid 100; partner on '+' at mid 130 (30bp downstream on genome).
    # Genome distance = +30. Anchor is '-', so it flips to -30 (upstream in anchor
    # reading frame). Relative strand is opposite (anchor '-', partner '+').
    ann = make_ann([
        ("chr1", 90, 110, "-", "A"),   # mid 100, anchor
        ("chr1", 120, 140, "+", "B"),  # mid 130, partner
    ])
    same, opp = spacing_analysis(ann, "A", "B", ceiling=1000)
    assert same.size == 0
    assert list(opp) == [-30]


def test_plus_strand_no_flip_same_strand():
    # Anchor '+' at 100, partner '+' at 130 -> +30, same strand.
    ann = make_ann([
        ("chr1", 90, 110, "+", "A"),
        ("chr1", 120, 140, "+", "B"),
    ])
    same, opp = spacing_analysis(ann, "A", "B", ceiling=1000)
    assert list(same) == [30]
    assert opp.size == 0


def test_overlap_excluded_by_default():
    # Only partner is an overlapping one -> excluded by default -> no distance.
    ann = make_ann([
        ("chr1", 100, 140, "+", "A"),   # 100-140
        ("chr1", 130, 170, "+", "B"),   # overlaps A
    ])
    same, opp = spacing_analysis(ann, "A", "B", ceiling=1000, exclude_overlap=True)
    assert same.size == 0 and opp.size == 0
    # with exclusion off, the (midpoint) distance is measured
    same2, _ = spacing_analysis(ann, "A", "B", ceiling=1000, exclude_overlap=False)
    assert same2.size == 1


def test_nearest_nonoverlapping_partner_chosen():
    # Anchor at 100-140. A closer partner overlaps (excluded); a farther one does not.
    ann = make_ann([
        ("chr1", 100, 140, "+", "A"),   # mid 120
        ("chr1", 130, 170, "+", "B"),   # overlaps -> excluded
        ("chr1", 300, 340, "+", "B"),   # mid 320 -> +200, kept
    ])
    same, _ = spacing_analysis(ann, "A", "B", ceiling=1000)
    assert list(same) == [200]


def test_peak_stats_flat_vs_tight():
    tight = np.zeros(50, dtype=int)              # all at 0
    conc, std, loc = peak_stats(tight)
    # fully concentrated, zero spread; peak lands in the bin straddling 0
    # (50bp bins, so the 0-value bin is [0,50) -> center 25, within one binwidth)
    assert conc == 1.0 and std == 0.0 and abs(loc) <= 50
    assert peak_stats(np.arange(3))[0] != peak_stats(np.arange(3))[0]  # nan for <10
