"""
motifdist
=========
Analyze the spatial relationships between transcription-factor motifs that have
already been mapped to genome coordinates (a ChromBPNet / TF-MoDISco seqlet
annotation table). Given *where* motifs land, motifdist asks:

  which motif pairs co-occur more than chance, and do any of them have a
  preferred spacing/orientation (a real **composite element**) rather than
  merely sharing the same accessible region (a **billboard**)?

The scientific payload is four artifact filters plus two distinct null models
(see the README). This package starts from an annotation table; producing that
table lives in the upstream Modisco_Annotation repo and is out of scope here.

Public API (the pieces most people import):

    from motifdist.io import load_annotation, pattern_cwm, pattern_ppm
    from motifdist.cooccurrence import count_near_pairs, cooccurrence_enrichment
    from motifdist.spacing import spacing_analysis, spacing_to_edges
    from motifdist.nulls import peak_height_ratio, spacing_null_test
    from motifdist.filters import (
        filter_small_number, filter_promiscuity, sibling_qvals, overlap_spike_qc,
    )
"""

__version__ = "0.1.0"

# Canonical pattern-label helpers, used across modules and tests.
REQUIRED_COLUMNS = ("chrom", "start", "end", "strand", "pattern", "subpattern")


def short_label(pattern):
    """'pos/pattern_2' -> 'p2', 'neg/pattern_10' -> 'n10'. Stable for odd names."""
    if "/" not in pattern:
        return pattern
    arm, pat = pattern.split("/", 1)
    return f"{arm[0]}{pat.replace('pattern_', '')}"


def pattern_to_h5(pattern):
    """'pos/pattern_2' -> ('pos_patterns', 'pattern_2').

    The annotation table labels patterns 'pos/pattern_N'; inside a MoDISco HDF5
    they live under the group 'pos_patterns'/'neg_patterns'. This is the single
    place that conversion happens.
    """
    arm, pat = pattern.split("/", 1)
    arm = arm.replace("pos", "pos_patterns").replace("neg", "neg_patterns")
    return arm, pat
