import numpy as np

from conftest import make_ann
from motifdist.cooccurrence import (
    count_near_pairs, count_near_pairs_no_overlap, shuffled_copy, cooccurrence_enrichment,
)


def test_count_near_pairs_hand_example():
    # A anchors at mid 100 (chr1), 1000 (chr1), 100 (chr2)
    A = make_ann([
        ("chr1", 90, 110, "+", "A"),    # mid 100
        ("chr1", 990, 1010, "+", "A"),  # mid 1000
        ("chr2", 90, 110, "+", "A"),    # mid 100
    ])
    # B seqlets: near the first A (150), far from the second, none on chr2
    B = make_ann([
        ("chr1", 140, 160, "+", "B"),   # mid 150 -> 50bp from A@100
        ("chr1", 5000, 5020, "+", "B"),  # mid 5010 -> far
    ])
    # ceiling 100: only A@100 on chr1 has a B within 100bp
    assert count_near_pairs(A, B, ceiling=100) == 1
    # ceiling 5000: A@100 (50) and A@1000 (min(150? no) ... nearest is 150=850? )
    # A@1000 nearest B is 150 (850bp) or 5010 (4010bp) -> 850 <= 5000 -> counts
    assert count_near_pairs(A, B, ceiling=5000) == 2


def test_count_near_pairs_no_overlap_excludes_overlap():
    # One A and one B that physically overlap on the genome.
    A = make_ann([("chr1", 100, 140, "+", "A")])   # 100-140
    B = make_ann([("chr1", 130, 170, "+", "B")])   # 130-170 overlaps A
    # midpoints are 120 and 150 -> 30bp apart, so the plain counter counts it
    assert count_near_pairs(A, B, ceiling=100) == 1
    # but they overlap, so the no-overlap counter excludes it
    assert count_near_pairs_no_overlap(A, B, ceiling=100) == 0


def test_shuffled_copy_is_deterministic_and_within_chrom():
    ann = make_ann([
        ("chr1", 0, 10, "+", "A"), ("chr1", 100, 110, "+", "B"),
        ("chr1", 200, 210, "+", "A"), ("chr2", 0, 10, "+", "A"),
    ])
    s1 = shuffled_copy(ann, np.random.default_rng(0))
    s2 = shuffled_copy(ann, np.random.default_rng(0))
    assert list(s1["mid"]) == list(s2["mid"])            # seeded -> deterministic
    # chr2 has a single seqlet, so its mid cannot move
    chr2 = s1[s1["chrom"] == "chr2"]
    assert list(chr2["mid"]) == [5]
    # the multiset of chr1 midpoints is preserved (only permuted)
    assert sorted(s1[s1["chrom"] == "chr1"]["mid"]) == [5, 105, 205]


def test_enrichment_deterministic_and_reports_pairs():
    ann = make_ann([
        ("chr1", 0, 10, "+", "A"), ("chr1", 20, 30, "+", "B"),
        ("chr1", 500, 510, "+", "A"), ("chr1", 520, 530, "+", "B"),
        ("chr2", 0, 10, "+", "A"), ("chr2", 900, 910, "+", "B"),
    ])
    d1 = cooccurrence_enrichment(ann, ceiling=100, n_shuffle=10, seed=1)
    d2 = cooccurrence_enrichment(ann, ceiling=100, n_shuffle=10, seed=1)
    assert d1.equals(d2)
    # both ordered pairs (A->B and B->A) present
    assert set(zip(d1["anchor"], d1["partner"])) == {("A", "B"), ("B", "A")}
