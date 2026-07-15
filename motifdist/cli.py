"""
cli.py  — the `motifdist` command-line interface
================================================
Four subcommands mirror the three pipeline stages plus report assembly:

    motifdist cooccur  --ann TABLE --ceiling 1000 --n-shuffles 200 --out cooc.csv
    motifdist filter   --cooc cooc.csv [--modisco H5 --host-mount DIR] \
                       --min-seqlets 200 --min-pairs 30 --max-partners 2 \
                       --out candidates.csv
    motifdist spacing  --ann TABLE --pairs candidates.csv --n-shuffles 200 \
                       --out spacing/
    motifdist report   --spacing spacing/ --out report.pdf

Every subcommand: explicit inputs, a deterministic --seed, stop-fail on missing
or invalid files, and structured logging of what was dropped at each step.
"""

import argparse
import logging
import os
import sys

from . import __version__


def _setup_logging(verbose):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )


def _require_file(path, what):
    if not os.path.isfile(path):
        raise SystemExit(f"error: {what} not found: {path}")


# --- cooccur ----------------------------------------------------------------

def cmd_cooccur(args):
    from .io import load_annotation
    from .cooccurrence import (
        cooccurrence_enrichment, count_near_pairs, count_near_pairs_no_overlap,
    )

    _require_file(args.ann, "annotation table (--ann)")
    ann = load_annotation(args.ann)
    counter = count_near_pairs_no_overlap if args.no_overlap else count_near_pairs
    cooc = cooccurrence_enrichment(
        ann, ceiling=args.ceiling, n_shuffle=args.n_shuffles, seed=args.seed, counter=counter,
    )
    cooc.sort_values("enrichment", ascending=False).to_csv(args.out, index=False)
    logging.getLogger("motifdist").info("wrote %d ordered pairs -> %s", len(cooc), args.out)


# --- filter -----------------------------------------------------------------

def cmd_filter(args):
    import pandas as pd
    from .filters import (
        filter_small_number, sibling_qvals, drop_siblings, filter_promiscuity,
    )

    log = logging.getLogger("motifdist")
    _require_file(args.cooc, "co-occurrence table (--cooc)")
    cooc = pd.read_csv(args.cooc)

    # Filter 1: small-number noise
    pairs = filter_small_number(
        cooc, min_seqlets=args.min_seqlets, min_pairs=args.min_pairs,
        min_enrichment=args.min_enrichment,
    )

    # Filter 2: sibling motifs (needs MoDISco h5 + Docker TOMTOM)
    if args.modisco and not args.skip_sibling:
        _require_file(args.modisco, "MoDISco HDF5 (--modisco)")
        from .tomtom import pattern_vs_pattern_qvals
        host_mount = args.host_mount or os.path.dirname(os.path.abspath(args.modisco))
        work_dir = args.tomtom_dir or os.path.join(host_mount, "pattern_vs_pattern_tomtom")
        log.info("running pattern-vs-pattern TOMTOM (host mount %s)", host_mount)
        qvals = pattern_vs_pattern_qvals(
            args.modisco, work_dir, host_mount, docker_image=args.docker_image,
        )
        pairs = drop_siblings(sibling_qvals(pairs, qvals), q_thresh=args.sibling_q)
    else:
        log.warning("sibling-motif filter SKIPPED (no --modisco given or --skip-sibling set); "
                    "co-occurring copies of the same TF will not be removed")

    # Filter 3: promiscuity
    pairs = filter_promiscuity(pairs, max_partners=args.max_partners)

    pairs.sort_values("enrichment", ascending=False).to_csv(args.out, index=False)
    log.info("wrote %d candidate composite pairs -> %s", len(pairs), args.out)


# --- spacing ----------------------------------------------------------------

def cmd_spacing(args):
    import pandas as pd
    from .io import load_annotation
    from .spacing import spacing_analysis, peak_stats
    from .nulls import spacing_null_test
    from .filters import overlap_spike_qc

    log = logging.getLogger("motifdist")
    _require_file(args.ann, "annotation table (--ann)")
    _require_file(args.pairs, "candidate pairs (--pairs)")
    ann = load_annotation(args.ann)
    cand = pd.read_csv(args.pairs)
    if not {"anchor", "partner"}.issubset(cand.columns):
        raise SystemExit("error: --pairs CSV must have 'anchor' and 'partner' columns")
    pair_list = list(zip(cand["anchor"], cand["partner"]))
    if not pair_list:
        raise SystemExit("error: --pairs CSV has no rows; nothing to analyze")

    os.makedirs(args.out, exist_ok=True)
    tf_of = {}  # optional TF labels; annotation may carry them
    if "pattern_tf_top3" in ann.columns:
        tf_of = ann.drop_duplicates("pattern").set_index("pattern")["pattern_tf_top3"].to_dict()

    # Per-pair signed distances (overlap-excluded), long form, + peak stats
    dist_rows, stat_rows, qc_rows = [], [], []
    for anchor, partner in pair_list:
        same, opp = spacing_analysis(ann, anchor, partner, ceiling=args.ceiling)
        for strand, arr in (("same", same), ("opposite", opp)):
            for d in arr:
                dist_rows.append({"anchor": anchor, "partner": partner,
                                  "strand": strand, "distance": int(d)})
            conc, std, loc = peak_stats(arr)
            stat_rows.append({"anchor": anchor, "partner": partner, "strand": strand,
                              "n": int(len(arr)), "peak_concentration": conc,
                              "spread_std": std, "peak_location_bp": loc})
        qc_rows.append(overlap_spike_qc(ann, anchor, partner))

    pd.DataFrame(dist_rows).to_csv(os.path.join(args.out, "spacing_distances.csv"), index=False)
    pd.DataFrame(stat_rows).to_csv(os.path.join(args.out, "spacing_stats.csv"), index=False)
    pd.DataFrame(qc_rows).to_csv(os.path.join(args.out, "overlap_qc.csv"), index=False)

    # Null model #2: is each peak real, or lumpiness?
    null_df = spacing_null_test(
        ann, pair_list, ceiling=args.ceiling, n_shuffles=args.n_shuffles, seed=args.seed,
    )
    null_df.to_csv(os.path.join(args.out, "spacing_null_test.csv"), index=False)

    # Figures: one histogram PNG per pair + a combined PDF
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.backends.backend_pdf import PdfPages
        from .viz import plot_spacing_fine

        pdf_path = os.path.join(args.out, "spacing_histograms.pdf")
        with PdfPages(pdf_path) as pdf:
            for anchor, partner in pair_list:
                from . import short_label
                png = os.path.join(args.out, f"spacing_{short_label(anchor)}_{short_label(partner)}.png")
                plot_spacing_fine(ann, anchor, partner, tf_of=tf_of,
                                  ceiling=args.ceiling, out_path=png, pdf=pdf)
        log.info("saved figures -> %s", pdf_path)
    except ImportError:
        log.warning("matplotlib not available; skipped figures (CSVs still written)")

    log.info("spacing analysis complete -> %s/", args.out)


# --- report -----------------------------------------------------------------

def cmd_report(args):
    import glob
    import pandas as pd

    log = logging.getLogger("motifdist")
    if not os.path.isdir(args.spacing):
        raise SystemExit(f"error: spacing directory not found: {args.spacing}")
    null_csv = os.path.join(args.spacing, "spacing_null_test.csv")
    _require_file(null_csv, "spacing_null_test.csv (run `motifdist spacing` first)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    null_df = pd.read_csv(null_csv)
    qc_path = os.path.join(args.spacing, "overlap_qc.csv")
    qc_df = pd.read_csv(qc_path) if os.path.isfile(qc_path) else None
    pngs = sorted(glob.glob(os.path.join(args.spacing, "spacing_*.png")))

    with PdfPages(args.out) as pdf:
        # Page 1: null-test summary table + interpretation
        fig, ax = plt.subplots(figsize=(12, max(3, len(null_df) * 0.4 + 3)))
        ax.axis("off")
        if len(null_df):
            tbl = ax.table(cellText=null_df.round(3).values, colLabels=null_df.columns,
                           loc="center", cellLoc="center")
            tbl.auto_set_font_size(False); tbl.set_fontsize(7); tbl.scale(1, 1.4)
        ax.set_title(
            "motifdist spacing report — null model #2 (spacing shuffle null)\n"
            "obs_peak = tallest bin / median bin. A pair is a real composite only "
            "when obs_peak beats the null (low p-value).\n"
            "Reminder: with several pairs tested, one significant hit is roughly "
            "what multiple testing predicts by chance. Use >=200 shuffles.",
            fontsize=10, pad=16)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Page 2: overlap-exclusion QC
        if qc_df is not None and len(qc_df):
            fig, ax = plt.subplots(figsize=(12, max(3, len(qc_df) * 0.4 + 2)))
            ax.axis("off")
            tbl = ax.table(cellText=qc_df.round(3).values, colLabels=qc_df.columns,
                           loc="center", cellLoc="center")
            tbl.auto_set_font_size(False); tbl.set_fontsize(7); tbl.scale(1, 1.4)
            ax.set_title("Overlap-exclusion QC — spike vs background overlap fraction\n"
                         "A spike that is ~100% overlapping pairs (vs ~20% background) is a "
                         "physical-overlap artifact, not real spacing.", fontsize=10, pad=16)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Remaining pages: the per-pair histogram PNGs
        for png in pngs:
            img = plt.imread(png)
            fig, ax = plt.subplots(figsize=(12, 5))
            ax.imshow(img); ax.axis("off")
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    log.info("wrote report -> %s (%d pairs, %d figures)", args.out, len(null_df), len(pngs))


# --- parser -----------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(prog="motifdist", description=__doc__.split("\n")[2])
    p.add_argument("--version", action="version", version=f"motifdist {__version__}")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("cooccur", help="Stage A: co-occurrence screen + null #1")
    c.add_argument("--ann", required=True, help="seqlet annotation CSV")
    c.add_argument("--ceiling", type=int, default=1000, help="max bp between seqlets (default 1000)")
    c.add_argument("--n-shuffles", type=int, default=200, help="null shuffles (default 200)")
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("--no-overlap", action="store_true",
                   help="exclude physically-overlapping seqlet pairs from the count")
    c.add_argument("--out", required=True, help="output co-occurrence CSV")
    c.set_defaults(func=cmd_cooccur)

    f = sub.add_parser("filter", help="Stage B: the four artifact filters")
    f.add_argument("--cooc", required=True, help="co-occurrence CSV from `cooccur`")
    f.add_argument("--modisco", help="MoDISco HDF5 (enables the sibling-motif filter)")
    f.add_argument("--jaspar", help="JASPAR .meme path (optional; never version-hardcoded)")
    f.add_argument("--min-seqlets", type=int, default=200)
    f.add_argument("--min-pairs", type=int, default=30)
    f.add_argument("--min-enrichment", type=float, default=1.5)
    f.add_argument("--max-partners", type=int, default=2)
    f.add_argument("--sibling-q", type=float, default=0.05)
    f.add_argument("--skip-sibling", action="store_true", help="skip the TOMTOM sibling filter")
    f.add_argument("--host-mount", help="absolute host dir to mount at /work for Docker TOMTOM")
    f.add_argument("--tomtom-dir", help="working dir for the TOMTOM MEME + output")
    f.add_argument("--docker-image", default="kundajelab/chrombpnet:latest")
    f.add_argument("--out", required=True, help="output candidate-pairs CSV")
    f.set_defaults(func=cmd_filter)

    s = sub.add_parser("spacing", help="Stage C: spacing geometry + null #2")
    s.add_argument("--ann", required=True, help="seqlet annotation CSV")
    s.add_argument("--pairs", required=True, help="candidate pairs CSV (anchor,partner)")
    s.add_argument("--ceiling", type=int, default=250, help="spacing window bp (default 250)")
    s.add_argument("--n-shuffles", type=int, default=200, help="spacing null shuffles (default 200)")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--out", required=True, help="output directory")
    s.set_defaults(func=cmd_spacing)

    r = sub.add_parser("report", help="assemble spacing outputs into a PDF")
    r.add_argument("--spacing", required=True, help="directory produced by `motifdist spacing`")
    r.add_argument("--out", required=True, help="output report.pdf")
    r.set_defaults(func=cmd_report)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    _setup_logging(args.verbose)
    args.func(args)


if __name__ == "__main__":
    main()
