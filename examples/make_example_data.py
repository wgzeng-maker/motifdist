"""
make_example_data.py
====================
Generate a small, synthetic seqlet annotation table that reproduces the two
outcomes the toolkit is built to tell apart end-to-end:

  * pos/pattern_0  (anchor)
  * pos/pattern_1  (COMPOSITE) -- in ~60% of peaks it sits exactly +40 bp
                                  downstream of the anchor, same strand; in the
                                  rest it sits at a random in-peak distance. The
                                  fixed-spacing majority makes a real spike that
                                  should SURVIVE the spacing null (#2); the random
                                  minority provides realistic background.
  * pos/pattern_2  (BILLBOARD) -- always a random in-peak distance from the
                                  anchor. It co-occurs with the anchor (so it
                                  passes the co-occurrence screen) but has no
                                  preferred spacing, so it should FAIL the null.

Peaks are sparse (a few seqlets each) and far apart, so (a) co-occurrence is
local to a peak and (b) two patterns that reliably share a peak come out enriched
above the label-shuffle null. This is what lets both real and billboard pairs
survive to Stage C, where the spacing null separates them.

The output CSV is checked in as example_annotation.csv; rerun only to regenerate.
"""

import numpy as np
import pandas as pd

SEED = 7
N_PEAKS = 300
PEAK_SPACING = 6000
COMPOSITE_OFFSET = 40          # fixed anchor->composite spacing (bp), same strand
COMPOSITE_FRACTION = 0.6       # fraction of peaks where pattern_1 takes that spacing
WIDTH = 20
MIN_OFFSET, MAX_OFFSET = 60, 240   # random in-peak offsets (avoid overlapping the anchor)


def _rand_offset(rng):
    mag = int(rng.integers(MIN_OFFSET, MAX_OFFSET))
    return mag if rng.random() < 0.5 else -mag


def main(out_path="example_annotation.csv"):
    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(N_PEAKS):
        center = 100_000 + i * PEAK_SPACING
        chrom = f"chr{i % 3 + 1}"

        anchor_strand = "+" if rng.random() < 0.5 else "-"
        direction = 1 if anchor_strand == "+" else -1
        rows.append((chrom, center, center + WIDTH, anchor_strand, "pos/pattern_0"))

        # composite partner (pattern_1): fixed +40 same-strand in most peaks,
        # random distance in the rest.
        if rng.random() < COMPOSITE_FRACTION:
            b_start = center + direction * COMPOSITE_OFFSET
            b_strand = anchor_strand
        else:
            b_start = center + _rand_offset(rng)
            b_strand = "+" if rng.random() < 0.5 else "-"
        rows.append((chrom, b_start, b_start + WIDTH, b_strand, "pos/pattern_1"))

        # billboard partner (pattern_2): always a random in-peak distance/strand.
        c_start = center + _rand_offset(rng)
        c_strand = "+" if rng.random() < 0.5 else "-"
        rows.append((chrom, c_start, c_start + WIDTH, c_strand, "pos/pattern_2"))

    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "pattern"])
    df["subpattern"] = ""
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)
    df.to_csv(out_path, index=False)
    print(f"wrote {out_path}: {len(df)} seqlets, {df['pattern'].nunique()} patterns")
    print(df["pattern"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "example_annotation.csv")
