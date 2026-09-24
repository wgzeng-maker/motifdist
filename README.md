# motifdist

Python tools for analyzing co-occurrence, spacing, and orientation of transcription-factor motifs mapped to genome coordinates.

Starting from a ChromBPNet / TF-MoDISco seqlet annotation table, the pipeline screens motif pairs for co-occurrence, filters potential artifacts, and tests spacing against a separate shuffle-based null model.

| At a glance | |
|---|---|
| Input | A genome-coordinate seqlet annotation CSV |
| Workflow | Co-occurrence screen → artifact filters → spacing analysis → report |
| Outputs | Co-occurrence and candidate tables, spacing results, figures, and a PDF report |
| Example | A bundled, seeded synthetic dataset with planted spacing and null cases |
| Upstream tool | [Modisco_Annotation](https://github.com/wgzeng-maker/Modisco_Annotation) produces the annotation table |

## Start here

- **Try the workflow:** [Runnable example](examples/README.md), using the bundled dataset.
- **Use your own data:** [Install](#install), [input contract](#input-contract), and [CLI](#cli).
- **Inspect the method:** [Scientific idea](#the-scientific-idea), [artifact filters](#the-four-filters-stage-b), and [signed distances](#the-signed-distance-convention-stage-c).
- **Inspect the implementation:** [Module layout](#module-layout) and [tests](tests/).

The tool starts from existing annotations; it does not call or back-annotate seqlets.
Spacing results are leads for validation, as discussed under [Interpretation and limitations](#interpretation-and-limitations).

## Install

```bash
pip install -e .            # core: numpy + pandas
pip install -e ".[viz,h5]"  # add figures (matplotlib, logomaker) + MoDISco I/O (h5py, hdf5plugin)
```

`hdf5plugin` is imported before `h5py` because many production MoDISco `.h5`
files use compressed HDF5 filters. TOMTOM runs inside the
`kundajelab/chrombpnet:latest` Docker image, so the sibling-motif filter needs
Docker available.

## Input contract

A seqlet annotation CSV with these columns (validated strictly on load; the
loader fails loudly rather than continuing on bad input):

| column | meaning |
|---|---|
| `chrom` | chromosome, e.g. `chr7` |
| `start`, `end` | **forward-strand** BED coordinates, always, `start < end` |
| `strand` | `+` / `-`; reading direction only |
| `pattern` | pattern label, e.g. `pos/pattern_2` |
| `subpattern` | subcluster label (present but may be empty) |

Optional secondary inputs (only for the sibling filter and logos): a TF-MoDISco
`.h5` and a JASPAR `.meme` database. **The JASPAR path is always an argument —
never a hardcoded version.**

## CLI

```bash
# Stage A — co-occurrence screen + null #1
motifdist cooccur  --ann TABLE --ceiling 1000 --n-shuffles 200 --out cooc.csv

# Stage B — the four filters
motifdist filter   --cooc cooc.csv --modisco MODISCO.h5 \
                   --host-mount /abs/path/mounted/at/work \
                   --min-seqlets 200 --min-pairs 30 --max-partners 2 \
                   --out candidates.csv
#   (omit --modisco or pass --skip-sibling to skip the TOMTOM sibling filter)

# Stage C — spacing geometry + null #2
motifdist spacing  --ann TABLE --pairs candidates.csv \
                   --n-shuffles 200 --out spacing/

# assemble the figures + tables into one PDF
motifdist report   --spacing spacing/ --out report.pdf
```

Every subcommand takes explicit inputs, a deterministic `--seed`, stop-fails on
missing/invalid files, and logs what was dropped at each step and why. Use
**≥200 shuffles for any reported result** (30 is only for a quick look; with 30
shuffles the smallest achievable p-value is 1/31 ≈ 0.032).

See [`examples/`](examples/) for a complete, seeded run on a small bundled dataset
that reproduces one significant (composite) and one null (billboard) result.

## Library

```python
from motifdist.io import load_annotation
from motifdist.cooccurrence import cooccurrence_enrichment
from motifdist.filters import filter_small_number, filter_promiscuity
from motifdist.spacing import spacing_analysis
from motifdist.nulls import spacing_null_test

ann  = load_annotation("annotation.csv")
cooc = cooccurrence_enrichment(ann, ceiling=1000, n_shuffle=200, seed=0)
cand = filter_promiscuity(filter_small_number(cooc), max_partners=2)
same, opp = spacing_analysis(ann, "pos/pattern_0", "pos/pattern_1")   # overlap-excluded
nulls = spacing_null_test(ann, [("pos/pattern_0", "pos/pattern_1")], n_shuffles=200)
```

## Module layout

```
motifdist/
  io.py            load + validate the annotation table; PPM/CWM loader from a MoDISco .h5
  cooccurrence.py  count_near_pairs, shuffled_copy (null #1), enrichment table
  spacing.py       signed-distance spacing_analysis, edge distances, peak stats
  nulls.py         peak_height_ratio + the spacing shuffle-null p-value (null #2)
  filters.py       small-number, sibling (TOMTOM), promiscuity, overlap-exclusion QC
  tomtom.py        Dockerized TOMTOM wrapper + version-robust output parsing
  viz.py           spacing histograms, motif logos, co-occurrence heatmap
  cli.py           the `motifdist` subcommands
```

## The scientific idea

### Billboard vs composite

Two motifs that appear near each other can mean two very different things:

- **Billboard** — both TFs bind the same open region because it's accessible, but
  in no particular arrangement. Distances between them are effectively random.
- **Composite element** — the two motifs bind at a *preferred spacing and
  orientation*, because the proteins physically cooperate. Distances pile up at a
  specific value.

Co-occurrence alone cannot tell these apart: both produce "these two patterns are
often near each other." Only the *geometry* of the spacing does.

### Why peaks make co-occurrence deceptive (null model #1)

Every seqlet already sits inside an accessible peak, so any two patterns are
automatically somewhat close just from peak geometry. **Null model #1** shuffles
*which pattern label* sits at each existing seqlet position, within each
chromosome. That preserves peak geometry and destroys only pattern-*specific*
placement, so enrichment = observed ÷ null isolates real affinity.

Two behaviors this surfaces (features, not bugs):
- The **median enrichment across all pairs sits near 1** (~0.8 in our reference
  run). If it drifts far from 1, the null is biased — `motifdist` reports it so
  you can see that.
- **Most pairs are depleted, not enriched**: patterns compete for limited space
  inside peaks.

### Why a tall histogram bin is cheap (null model #2)

Throw ~600 distances into ~100 bins and the tallest bin will be several times the
typical bin *by luck alone* — randomness is lumpy. So a "peak" in a spacing
histogram is not evidence of preferred spacing until it beats chance. **Null model
#2** shuffles seqlet positions, rebuilds the spacing histogram, and records how
tall the tallest bin gets by chance. The reported p-value is the fraction of
shuffles whose peak is at least as tall as the observed one.

> In our reference run the NFI→Zic "peak" (4.6×) was *exactly* what chance
> produces (4.7×). Without this null you would report a spacing preference that
> does not exist. **The two nulls are different and must not be conflated:** #1 is
> about co-occurrence, #2 is about spacing.

## The four filters (Stage B)

Applied to co-occurrence hits to drop artifacts. Each catches something the
others miss — keep all four.

1. **Small-number noise.** The top raw hits in our reference run were 40×
   "enriched" but built on **1–2 observations**. Require both patterns to be
   reasonably large (`--min-seqlets`, default 200) and the pair to have enough
   observed co-occurrences (`--min-pairs`, default 30).

2. **Sibling-motif filter.** Two patterns can be the *same TF* captured twice;
   their co-occurrence is then trivial. `motifdist` runs **pattern-vs-pattern
   TOMTOM** and drops pairs whose motifs match (q < 0.05). This correctly caught
   pattern_11 ↔ pattern_16 (both ZNF143, q = 2.7e-06). It uses **real TOMTOM**,
   not a homemade similarity score (an earlier cosine-style CWM score ranked that
   very ZNF143 sibling pair as *less* similar than unrelated pairs — backwards —
   so it is not used for filtering).

3. **Promiscuity filter.** A pattern that pairs with *everything* (Zic = 11, YY2 =
   10 partners in our reference run) marks a busy regulatory neighborhood, not a
   specific composite. Keep only pairs where **both** partners are selective
   (`--max-partners`).

4. **Overlap-exclusion filter — novel.** Two seqlets in *different* patterns can
   physically overlap on the genome (sharing up to ~49% of their bases), because
   MoDISco's overlap suppression only prevents >50% overlap *within a single
   extraction pass*, not across patterns. Overlapping pairs produce a **phantom
   tight-spacing spike (~15 bp)** that survives the sibling test. `motifdist`
   **excludes physically-overlapping seqlet pairs before any spacing analysis**
   (on by default). The QC diagnostic is a clean tell: in our reference run the
   phantom spike was **100% overlapping pairs vs ~20% in the background**.

## The signed-distance convention (Stage C)

Distances are measured in the **anchor's reading frame**:

```
signed_dist = partner_mid - anchor_mid          # + = partner downstream on genome
if anchor.strand == '-':  signed_dist *= -1      # read the anchor "forward"
```

Only the anchor is flipped. Relative strand ("same" vs "opposite") is decided by
strand *equality*, which is frame-independent, so the partner is never flipped.
Positive = downstream in the anchor's reading direction; negative = upstream.
Same-strand and opposite-strand distances are reported as separate histograms.

## Interpretation and limitations

When you test several pairs, a single significant hit is roughly what multiple
testing predicts by chance. Treat one significant spacing result among many as a
lead to validate, not a finding. Report results with **≥200 shuffles**, and read
the co-occurrence sanity check (median enrichment near 1) before trusting the
enrichment numbers. `motifdist` is built to help you *not* oversell a peak.

## Development

```bash
pip install -e ".[dev]"
pytest
```

Tests cover the load/validate contract, both null models (seeded → deterministic),
the signed-distance frame math (a hand-checked minus-strand example), each
filter's drop logic, and the TOMTOM output parser (against saved `tomtom.txt` /
`tomtom.tsv` fixtures — the live-Docker path is not exercised in unit tests).

## License

MIT — see [LICENSE](LICENSE).
