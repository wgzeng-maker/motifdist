import sys
from pathlib import Path

import pandas as pd
import pytest

# Make the package importable without an install step.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def make_ann(rows):
    """Build a validated-shape annotation DataFrame from (chrom, start, end,
    strand, pattern) tuples, filling subpattern and the derived `mid` column."""
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "pattern"])
    df["subpattern"] = ""
    df["mid"] = (df["start"] + df["end"]) // 2
    return df


@pytest.fixture
def make_ann_fixture():
    return make_ann
