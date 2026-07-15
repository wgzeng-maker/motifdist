# Runnable example

A tiny synthetic dataset that reproduces the two outcomes the toolkit exists to
tell apart: one **real composite** (significant preferred spacing) and one
**billboard** (co-occurs, but no preferred spacing — a null result).

`example_annotation.csv` is checked in (900 seqlets, 3 patterns). Regenerate it
with `python make_example_data.py` if you like; it is fully seeded.

The three patterns:

| pattern | role | placement |
|---|---|---|
| `pos/pattern_0` | anchor | one per peak |
| `pos/pattern_1` | **composite** | +40 bp downstream of the anchor, same strand, in ~60% of peaks (random distance in the rest) |
| `pos/pattern_2` | **billboard** | random in-peak distance from the anchor, always |

## Run the whole pipeline

```bash
# from this examples/ directory, with the package importable
# (pip install -e .. , or set PYTHONPATH=.. )

# Stage A — co-occurrence screen + null #1
python -m motifdist cooccur \
  --ann example_annotation.csv --ceiling 1000 --n-shuffles 200 --seed 0 \
  --out cooc.csv

# Stage B — the four filters (sibling filter needs a MoDISco .h5 + Docker;
# skipped here because this synthetic set has no motif matrices)
python -m motifdist filter \
  --cooc cooc.csv --min-seqlets 200 --min-pairs 30 --max-partners 2 \
  --skip-sibling --out candidates.csv

# Stage C — spacing geometry + null #2
python -m motifdist spacing \
  --ann example_annotation.csv --pairs candidates.csv \
  --ceiling 250 --n-shuffles 200 --seed 0 --out spacing_out

# assemble the PDF
python -m motifdist report --spacing spacing_out --out report.pdf
```

## What you should see

Both partner pairs pass the co-occurrence screen (each is ~1.8× enriched — both
genuinely share peaks with the anchor). The spacing null (`spacing_out/
spacing_null_test.csv`) is where they separate:

```
anchor          partner         strand    n    obs_peak  null_mean  pvalue
pos/pattern_0   pos/pattern_1   same      240  184.0     ~12        0.005   <- COMPOSITE (real +40bp spacing)
pos/pattern_0   pos/pattern_1   opposite  60   4.0       ~12        1.0
pos/pattern_0   pos/pattern_2   same      153  3.5       ~12        1.0     <- BILLBOARD (co-occurs, no spacing)
pos/pattern_0   pos/pattern_2   opposite  147  2.5       ~12        1.0
```

The composite's same-strand peak (184×) towers over the null (~12×), p = 0.005.
The billboard's tallest bin (2.5–3.5×) is *below* what chance produces, p = 1.0.
That contrast — not the co-occurrence enrichment, which is nearly identical for
both — is the whole point of the spacing null.
