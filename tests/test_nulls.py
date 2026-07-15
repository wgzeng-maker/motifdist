import numpy as np

from conftest import make_ann
from motifdist.nulls import peak_height_ratio, spacing_null_test


def test_peak_height_ratio_hand_value():
    # 20 distances stacked in the 0-bin (one tall bin of height 20), plus 20 more
    # spread one-per-bin across distinct 5bp bins (each height 1). Median nonzero
    # bin height = 1, tallest = 20 -> ratio 20.
    bins = np.arange(-250, 251, 5)
    spread = [(-245 + 5 * k) for k in range(20)]   # 20 distinct bins, one each
    distances = np.array([0] * 20 + spread)
    assert peak_height_ratio(distances, bins=bins) == 20.0


def test_peak_height_ratio_too_few():
    assert np.isnan(peak_height_ratio(np.arange(5)))


def test_spacing_null_test_seeded_deterministic():
    # A small dataset with a couple of patterns; just needs to run and be stable.
    rng = np.random.default_rng(0)
    rows = []
    for i in range(60):
        pos = int(rng.integers(0, 10000))
        rows.append(("chr1", pos, pos + 20, "+", "A"))
        rows.append(("chr1", pos + 40, pos + 60, "+", "B"))
    ann = make_ann(rows)
    d1 = spacing_null_test(ann, [("A", "B")], ceiling=250, n_shuffles=10, seed=42)
    d2 = spacing_null_test(ann, [("A", "B")], ceiling=250, n_shuffles=10, seed=42)
    assert d1.equals(d2)
    assert set(d1["strand"]) <= {"same", "opposite"}
    assert {"obs_peak", "null_mean", "null_max", "pvalue"}.issubset(d1.columns)
