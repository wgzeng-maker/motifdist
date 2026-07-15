"""
viz.py  — figures: co-occurrence heatmap, spacing histograms, motif logos
=========================================================================
All plotting lives here so the analysis modules stay import-light. matplotlib is
imported lazily inside each function, so you can use the rest of motifdist without
it installed.
"""

import logging

import numpy as np

from . import short_label
from .spacing import spacing_analysis, peak_stats

log = logging.getLogger("motifdist.viz")

_SAME_COLOR = "#4C72B0"
_OPP_COLOR = "#C44E52"


def plot_cooccurrence_heatmap(cooc, out_path=None, min_size=0, ax=None):
    """Heatmap of log2 co-occurrence enrichment (red = attract, blue = avoid,
    white = chance). `cooc` is the tidy DataFrame from cooccurrence_enrichment.
    `min_size` restricts to patterns with at least that many seqlets.
    """
    import matplotlib.pyplot as plt

    patterns = sorted(set(cooc["anchor"]) | set(cooc["partner"]))
    size_of = dict(zip(cooc["anchor"], cooc["anchor_size"]))
    size_of.update(dict(zip(cooc["partner"], cooc["partner_size"])))
    if min_size:
        patterns = [p for p in patterns if size_of.get(p, 0) >= min_size]
    patterns = sorted(patterns, key=lambda p: -size_of.get(p, 0))
    idx = {p: i for i, p in enumerate(patterns)}
    n = len(patterns)

    mat = np.full((n, n), np.nan)
    for _, r in cooc.iterrows():
        if r["anchor"] in idx and r["partner"] in idx and r["enrichment"] > 0:
            mat[idx[r["anchor"]], idx[r["partner"]]] = np.log2(r["enrichment"])

    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(max(8, n * 0.4), max(7, n * 0.35)))
    vmax = np.nanpercentile(np.abs(mat), 98) if np.isfinite(mat).any() else 1.0
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    labels = [f"{short_label(p)} ({size_of.get(p, 0)})" for p in patterns]
    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("partner pattern (B)"); ax.set_ylabel("anchor pattern (A)")
    ax.set_title("Co-occurrence enrichment: log2(observed / null)\n"
                 "red = near more than chance, blue = avoid, white = chance")
    if created:
        fig.colorbar(im, ax=ax, label="log2 enrichment", shrink=0.6)
        fig.tight_layout()
        if out_path:
            fig.savefig(out_path, dpi=150)
            log.info("saved %s", out_path)
        return fig
    return im


def _two_panel_spacing(ann, anchor, partner, tf_of, ceiling, bins, exclude_overlap, title_suffix):
    import matplotlib.pyplot as plt

    same, opp = spacing_analysis(ann, anchor, partner, ceiling=ceiling,
                                 exclude_overlap=exclude_overlap)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, data, title, color in (
        (axes[0], same, "SAME strand  (anchor →, partner →)", _SAME_COLOR),
        (axes[1], opp, "OPPOSITE strand  (anchor →, partner ←)", _OPP_COLOR),
    ):
        if len(data):
            ax.hist(data, bins=bins, color=color, edgecolor="white")
            conc, _, loc = peak_stats(data, bins=bins)
            ax.axvline(loc, color="black", ls="--", lw=1, alpha=0.6)
            ax.set_title(f"{title}\nn={len(data)}  conc={conc}  peak@{loc:.0f}bp", fontsize=10)
        else:
            ax.set_title(f"{title}\nn=0", fontsize=10)
        ax.axvline(0, color="black", lw=1)
        ax.set_xlabel("signed distance to partner (bp)\nneg=upstream, pos=downstream (anchor frame)")
    axes[0].set_ylabel("count of anchor seqlets")
    tfA = (tf_of or {}).get(anchor) or anchor
    tfB = (tf_of or {}).get(partner) or partner
    fig.suptitle(f"Spacing: {anchor} ({tfA})  →  {partner} ({tfB}){title_suffix}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


def plot_spacing(ann, anchor, partner, tf_of=None, ceiling=1000, binwidth=40,
                 exclude_overlap=True, out_path=None, pdf=None):
    """Coarse two-panel spacing histogram (same-strand | opposite-strand)."""
    bins = np.arange(-ceiling, ceiling + binwidth, binwidth)
    suffix = "  [overlap-excluded]" if exclude_overlap else ""
    fig = _two_panel_spacing(ann, anchor, partner, tf_of, ceiling, bins, exclude_overlap, suffix)
    _finish(fig, out_path, pdf)
    return fig


def plot_spacing_fine(ann, anchor, partner, tf_of=None, ceiling=250, binwidth=5,
                      exclude_overlap=True, out_path=None, pdf=None):
    """Fine-resolution spacing histogram. A tight composite shows a sharp spike;
    a billboard smears into a broad hump."""
    bins = np.arange(-ceiling, ceiling + binwidth, binwidth)
    suffix = f"  [fine ±{ceiling}bp, overlap-excluded]" if exclude_overlap else f"  [fine ±{ceiling}bp]"
    fig = _two_panel_spacing(ann, anchor, partner, tf_of, ceiling, bins, exclude_overlap, suffix)
    _finish(fig, out_path, pdf)
    return fig


def plot_pattern_logos(modisco_h5, patterns, tf_of=None, out_path=None, pdf=None):
    """Plot each pattern's CWM as a sequence logo (one row per pattern). Useful
    for eyeballing whether two 'different' patterns are halves of one motif."""
    import matplotlib.pyplot as plt
    import logomaker
    import pandas as pd
    from .io import pattern_cwm

    fig, axes = plt.subplots(len(patterns), 1, figsize=(10, 2.2 * len(patterns)))
    if len(patterns) == 1:
        axes = [axes]
    for ax, p in zip(axes, patterns):
        cwm = pattern_cwm(modisco_h5, p)
        df = pd.DataFrame(cwm, columns=["A", "C", "G", "T"])
        logomaker.Logo(df, ax=ax)
        name = (tf_of or {}).get(p) or "unnamed"
        ax.set_title(f"{p}  ({name})", fontsize=10)
    fig.tight_layout()
    _finish(fig, out_path, pdf)
    return fig


def _finish(fig, out_path, pdf):
    import matplotlib.pyplot as plt

    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    if out_path:
        fig.savefig(out_path, dpi=140, bbox_inches="tight")
        log.info("saved %s", out_path)
    if pdf is not None or out_path:
        plt.close(fig)
