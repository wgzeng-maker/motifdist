"""
tomtom.py  — motif comparison via TOMTOM (Dockerized MEME)
=========================================================
Two jobs:
  1. Write MoDISco patterns to a MEME file (PPMs), so they can be compared.
  2. Run TOMTOM (in the ChromBPNet Docker image) and parse its output robustly.

This wrapper exists because the *sibling-motif filter* needs real motif
comparison. The notebook also contained a homemade cosine-style CWM similarity
score; it is unreliable (it ranked a known ZNF143 sibling pair as *less* similar
than unrelated pairs — backwards), so it is NOT used for filtering. Use TOMTOM.

Bugs this module fixes (they bit the real analysis; see TOOLKIT_SPEC.md §4):
  * Output filename is version-dependent: newer MEME writes `tomtom.tsv`, the box
    this ran on writes `tomtom.txt`. We glob for both.
  * The `.txt` output can be headerless: we parse with explicit column names and
    detect whether a `#`/`Query_ID` header line is present.
  * JASPAR version drift: the `.meme` database path is always a caller argument,
    never a hardcoded version.
  * Docker hygiene: mount paths are absolute (never `~`); TOMTOM results are read
    from `-text` stdout so no root-owned files are created on the host.
"""

import glob
import logging
import os
import subprocess

import numpy as np
import pandas as pd

log = logging.getLogger("motifdist.tomtom")

# TOMTOM's standard column order, used when the output has no header line.
TOMTOM_COLUMNS = [
    "Query_ID", "Target_ID", "Optimal_offset", "p-value", "E-value", "q-value",
    "Overlap", "Query_consensus", "Target_consensus", "Orientation",
]

DEFAULT_DOCKER_IMAGE = "kundajelab/chrombpnet:latest"


def _trim_pfm(pfm, frac=0.3):
    """Trim low-information flanks so short motifs match cleanly."""
    ic = 2 + (pfm * np.log2(pfm + 1e-9)).sum(axis=1)
    keep = np.where(ic >= frac * ic.max())[0]
    return pfm[keep.min(): keep.max() + 1] if len(keep) else pfm


def write_patterns_meme(modisco_h5, out_path, trim_frac=0.3):
    """Write every MoDISco pattern's PPM to a MEME-format file at `out_path`.

    Motif names are 'pos.pattern_0' / 'neg.pattern_3' (dot form), which
    `parse_tomtom` converts back to the table's 'pos/pattern_0' labels.
    """
    from .io import open_h5

    with open_h5(modisco_h5) as f, open(out_path, "w") as out:
        out.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n")
        out.write("Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        for arm in ("pos_patterns", "neg_patterns"):
            if arm not in f:
                continue
            for pat in f[arm]:
                pfm = _trim_pfm(f[arm][pat]["sequence"][:], frac=trim_frac)
                name = f"{arm.split('_')[0]}.{pat}"  # 'pos.pattern_0'
                out.write(f"MOTIF {name}\n")
                out.write(f"letter-probability matrix: alength= 4 w= {len(pfm)}\n")
                for row in pfm:
                    r = row / row.sum()
                    out.write(f" {r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {r[3]:.6f}\n")
                out.write("\n")
    log.info("wrote patterns MEME: %s", out_path)
    return out_path


def _dot_to_slash(x):
    """'pos.pattern_0' -> 'pos/pattern_0'; leave JASPAR/other IDs unchanged."""
    if isinstance(x, str) and x.count(".") == 1 and x.split(".")[0] in ("pos", "neg"):
        arm, pat = x.split(".")
        return f"{arm}/{pat}"
    return x


def parse_tomtom(path):
    """Parse a TOMTOM text output into a DataFrame, robust to a missing header.

    Works for both `-text` stdout captured to a file and TOMTOM's own
    `tomtom.tsv`/`tomtom.txt`. Pattern IDs in dot form are normalized to slash
    form. Numeric columns (offset, p/E/q-value, overlap) are coerced to numbers.
    """
    # Real TOMTOM output carries a trailing footer of '#'-comment lines (version
    # and format notes), so we always skip '#' lines. The header row, when present,
    # is the first non-comment line (`tomtom.tsv`); a headerless `tomtom.txt` has a
    # data row there instead.
    with open(path) as fh:
        first_data = next((ln for ln in fh if ln.strip() and not ln.startswith("#")), "")
    has_header = first_data.startswith("Query_ID") or "Query_ID\t" in first_data

    if has_header:
        df = pd.read_csv(path, sep="\t", comment="#")
        df = df.rename(columns={c: c.strip() for c in df.columns})
    else:
        df = pd.read_csv(path, sep="\t", header=None, names=TOMTOM_COLUMNS, comment="#")

    # tolerate either 'q-value' or 'qval'-style names
    rename = {}
    for c in df.columns:
        lc = c.lower().replace("-", "").replace("_", "")
        if lc in ("qvalue", "qval"):
            rename[c] = "q-value"
        elif lc in ("pvalue", "pval"):
            rename[c] = "p-value"
    df = df.rename(columns=rename)

    for col in ("Query_ID", "Target_ID"):
        if col in df.columns:
            df[col] = df[col].map(_dot_to_slash)
    for col in ("Optimal_offset", "p-value", "E-value", "q-value", "Overlap"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _find_tomtom_output(outdir):
    """Return the TOMTOM output file in `outdir`, tolerating .tsv vs .txt."""
    for name in ("tomtom.tsv", "tomtom.txt"):
        hits = glob.glob(os.path.join(outdir, name))
        if hits:
            return hits[0]
    raise FileNotFoundError(
        f"No tomtom.tsv or tomtom.txt found in {outdir!r}. TOMTOM may have failed; "
        "check that the Docker image ran and the input MEME file is valid."
    )


def run_tomtom(query_meme, target_meme, out_tsv, host_mount, docker_image=DEFAULT_DOCKER_IMAGE,
               min_overlap=5, dist="pearson", thresh=1.0):
    """Run TOMTOM inside Docker and capture results to `out_tsv` on the host.

    query_meme / target_meme are paths *inside* the container (under /work).
    `host_mount` is the absolute host directory mounted at /work — it must be
    absolute (never '~'); this is enforced. Results are read from TOMTOM's `-text`
    stdout, so no root-owned output files are written to the host.

    Returns `out_tsv`. Requires Docker; raises if the run fails.
    """
    host_mount = os.path.abspath(os.path.expanduser(host_mount))
    if "~" in host_mount:
        raise ValueError("host_mount must be an absolute path, not contain '~'.")

    cmd = [
        "docker", "run", "--rm", "-v", f"{host_mount}:/work", docker_image,
        "tomtom", "-no-ssc", "-oc", ".", "--verbosity", "1", "-text",
        "-min-overlap", str(min_overlap), "-dist", dist, "-thresh", str(thresh),
        query_meme, target_meme,
    ]
    log.info("running TOMTOM: %s", " ".join(cmd))
    with open(out_tsv, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"TOMTOM failed (rc={proc.returncode}):\n{proc.stderr}")
    return out_tsv


def pattern_vs_pattern_qvals(modisco_h5, work_dir, host_mount, docker_image=DEFAULT_DOCKER_IMAGE):
    """Run pattern-vs-pattern TOMTOM (every pattern against every other) and
    return a DataFrame with columns [Query_ID, Target_ID, q-value].

    `work_dir` is a host directory (under `host_mount`) for the MEME file and the
    captured output. Convenience wrapper around write_patterns_meme + run_tomtom +
    parse_tomtom; the sibling filter calls this.
    """
    os.makedirs(work_dir, exist_ok=True)
    meme_host = os.path.join(work_dir, "all_patterns.meme")
    write_patterns_meme(modisco_h5, meme_host)

    host_mount = os.path.abspath(os.path.expanduser(host_mount))
    rel = os.path.relpath(meme_host, host_mount)
    meme_container = f"/work/{rel}"
    out_tsv = os.path.join(work_dir, "pattern_vs_pattern.tsv")

    run_tomtom(meme_container, meme_container, out_tsv, host_mount, docker_image=docker_image)
    df = parse_tomtom(out_tsv)
    return df[["Query_ID", "Target_ID", "q-value"]]
