from pathlib import Path

import numpy as np

from motifdist.tomtom import parse_tomtom

FIX = Path(__file__).parent / "fixtures"


def test_parse_headerless_txt():
    # The environment this ran in wrote a headerless `tomtom.txt` (TOOLKIT_SPEC §4.2).
    df = parse_tomtom(FIX / "tomtom_headerless.txt")
    assert len(df) == 3
    # dot-form IDs are normalized to slash form
    assert df.iloc[0]["Query_ID"] == "pos/pattern_11"
    assert df.iloc[0]["Target_ID"] == "pos/pattern_16"
    # the ZNF143 sibling pair has a tiny q-value, parsed as a float
    assert np.isclose(df.iloc[0]["q-value"], 2.7e-06)
    assert df["q-value"].dtype.kind == "f"


def test_parse_with_header_tsv():
    # Newer MEME writes `tomtom.tsv` with a header (TOOLKIT_SPEC §4.1).
    df = parse_tomtom(FIX / "tomtom_with_header.tsv")
    assert len(df) == 2
    assert df.iloc[0]["Query_ID"] == "pos/pattern_11"
    assert np.isclose(df.iloc[1]["q-value"], 0.2)


def test_both_formats_agree_on_qvalues():
    a = parse_tomtom(FIX / "tomtom_headerless.txt")
    b = parse_tomtom(FIX / "tomtom_with_header.tsv")
    key = ("pos/pattern_11", "pos/pattern_16")
    qa = a[(a.Query_ID == key[0]) & (a.Target_ID == key[1])]["q-value"].iloc[0]
    qb = b[(b.Query_ID == key[0]) & (b.Target_ID == key[1])]["q-value"].iloc[0]
    assert np.isclose(qa, qb)
