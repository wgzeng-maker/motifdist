import pandas as pd
import pytest

from motifdist.io import load_annotation

VALID = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr2"],
    "start": [100, 500, 100],
    "end": [120, 520, 130],
    "strand": ["+", "-", "+"],
    "pattern": ["pos/pattern_0", "pos/pattern_1", "pos/pattern_0"],
    "subpattern": ["", "", ""],
})


def _write(df, tmp_path, name="ann.csv"):
    p = tmp_path / name
    df.to_csv(p, index=False)
    return str(p)


def test_valid_loads_and_adds_mid(tmp_path):
    df = load_annotation(_write(VALID, tmp_path))
    assert list(df["mid"]) == [110, 510, 115]
    assert len(df) == 3


def test_missing_column_raises(tmp_path):
    bad = VALID.drop(columns=["strand"])
    with pytest.raises(ValueError, match="missing required column"):
        load_annotation(_write(bad, tmp_path))


def test_start_ge_end_raises(tmp_path):
    bad = VALID.copy()
    bad.loc[1, "end"] = bad.loc[1, "start"]  # start == end
    with pytest.raises(ValueError, match="start >= end"):
        load_annotation(_write(bad, tmp_path))


def test_bad_strand_raises(tmp_path):
    bad = VALID.copy()
    bad.loc[0, "strand"] = "*"
    with pytest.raises(ValueError, match="strand"):
        load_annotation(_write(bad, tmp_path))


def test_nan_in_required_raises(tmp_path):
    bad = VALID.copy()
    bad.loc[0, "chrom"] = None
    with pytest.raises(ValueError, match="NaN"):
        load_annotation(_write(bad, tmp_path))


def test_fractional_coordinate_raises(tmp_path):
    bad = VALID.copy().astype({"start": float})
    bad.loc[0, "start"] = 100.5
    with pytest.raises(ValueError, match="integer"):
        load_annotation(_write(bad, tmp_path))
