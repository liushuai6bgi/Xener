#!/usr/bin/env python3
"""Post-run quality gate for xener pipelines.

This script is invoked automatically at the end of run_pipeline.py (in-process,
via ``run_gate()`` — no subprocess, no h5ad re-read), and can also be run
standalone:

    python scripts/check_output.py --outdir output/edf/

It reads celltype_weight.csv and the run log in --outdir, computes five
quality signals, and exits non-zero if any threshold is breached. The
thresholds are intentionally tight: a passing run is a run whose output
is likely to be biologically useful, not just a run that produced files.

Cluster sizes for the weak-cluster check are read from the lightweight
``{dataset}_annotation.csv`` (auto-detected in --outdir, or pass
``--annotation-csv``), NOT from the multi-GB h5ad — see mandatory-rules.md
sec.10 (single inspection / I/O discipline). ``--h5ad`` remains as a legacy
fallback only.

Skill context: this is the mandatory Step 5.5 quality gate defined in
references/workflows/self-tuning-protocol.md. Do not declare an xener
run done without running this and seeing it pass.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


# Failure thresholds. Tightened to catch the common "ran clean, output
# garbage" failure mode where model_species is too narrow for the organ.
THRESH_MEAN_KG_MISS = 0.30
THRESH_TAIL_KG_MISS = 0.80
THRESH_TAIL_KG_MISS_FRAC = 0.05
THRESH_MIN_UNIQUE_CELLTYPES = 5
THRESH_MIN_CLUSTERS_FOR_DIVERSITY = 10
THRESH_WEAK_INIT_WEIGHT = 50.0
THRESH_ALLOW_WEAK_CLUSTERS = 1


def find_run_log(outdir: Path) -> Path | None:
    """Locate the most recent xener run log in outdir.

    run_pipeline.py mirrors xener's stdout to outdir/xener.log, which is
    the primary location. Older runs may have the log elsewhere; we look
    for the most plausible candidate in outdir first, then fall back to
    common places.
    """
    candidates = [
        outdir / "xener.log",
        outdir / "run.log",
        outdir.parent / "xener.log",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def parse_kg_miss_per_cluster(log_path: Path) -> dict[int, float]:
    """Parse `total X% homolos of organ[...] not in kg` lines from the log.

    Returns dict mapping cluster_id -> miss fraction (0.0 to 1.0).
    Cluster IDs come from the most recent "processing cluster N" header.
    """
    text = log_path.read_text(errors="ignore")
    cov: dict[int, float] = {}
    current_cluster: int | None = None
    in_annotation = False
    for line in text.splitlines():
        if "cell annotation organ" in line:
            in_annotation = True
            continue
        if not in_annotation:
            continue
        m = re.search(r"processing cluster (\d+)", line)
        if m:
            current_cluster = int(m.group(1))
        m = re.search(r"total ([\d.]+) % homolos of organ\[[^\]]+\] not in kg", line)
        if m and current_cluster is not None:
            cov[current_cluster] = float(m.group(1)) / 100.0
    return cov


def check_kg_miss(cov: dict[int, float]) -> tuple[bool, list[str]]:
    msgs = []
    ok = True
    if not cov:
        msgs.append("[WARN] check 1+2: no KG miss data in log; "
                    "could not parse 'total X% homolos ... not in kg' lines.")
        return True, msgs
    n = len(cov)
    mean_miss = sum(cov.values()) / n
    tail = sum(1 for v in cov.values() if v > THRESH_TAIL_KG_MISS)
    tail_frac = tail / n

    msgs.append(f"[INFO] check 1+2: mean KG miss = {mean_miss:.1%}, "
                f"clusters with >{THRESH_TAIL_KG_MISS:.0%} miss = {tail}/{n} "
                f"({tail_frac:.1%})")
    if mean_miss > THRESH_MEAN_KG_MISS:
        ok = False
        msgs.append(f"[FAIL] check 1: mean KG miss {mean_miss:.1%} "
                    f"> {THRESH_MEAN_KG_MISS:.0%}. "
                    "model_species likely too narrow for the chosen organ. "
                    "See references/workflows/species-selection.md worked example.")
    if tail_frac > THRESH_TAIL_KG_MISS_FRAC:
        ok = False
        msgs.append(f"[FAIL] check 2: {tail}/{n} clusters "
                    f"({tail_frac:.1%}) have KG miss > "
                    f"{THRESH_TAIL_KG_MISS:.0%} (threshold "
                    f"{THRESH_TAIL_KG_MISS_FRAC:.0%}).")
    return ok, msgs


def check_celltype_diversity(ct_path: Path) -> tuple[bool, list[str]]:
    msgs = []
    if not ct_path.exists():
        return False, [f"[FAIL] celltype_weight.csv not found at {ct_path}"]
    df = pd.read_csv(ct_path)
    if "cluster" not in df.columns or "celltype" not in df.columns \
            or "init_weight" not in df.columns:
        return False, ["[FAIL] celltype_weight.csv missing required columns "
                       "(cluster, celltype, init_weight)"]
    top1 = df.loc[df.groupby("cluster")["init_weight"].idxmax()]
    n_clusters = top1["cluster"].nunique()
    n_unique_types = top1["celltype"].nunique()
    msgs.append(f"[INFO] check 3: {n_unique_types} unique top-1 cell types "
                f"across {n_clusters} clusters")
    if n_clusters > THRESH_MIN_CLUSTERS_FOR_DIVERSITY \
            and n_unique_types < THRESH_MIN_UNIQUE_CELLTYPES:
        return False, [f"[FAIL] check 3: only {n_unique_types} unique top-1 "
                       f"cell types for {n_clusters} clusters "
                       f"(expected >= {THRESH_MIN_UNIQUE_CELLTYPES}). "
                       "model_species is too narrow or organ filter is wrong."]
    return True, msgs


def cluster_sizes_from_annotation(
    annotation_csv: Path | None, cluster_key: str | None
) -> tuple[dict, list[str]]:
    """Derive per-cluster cell counts from the lightweight annotation CSV.

    This is the I/O-cheap replacement for re-reading the multi-GB h5ad: the
    `{dataset}_annotation.csv` written by run_pipeline.py carries one row per
    cell with its cluster label, so a value_counts() of the cluster column
    yields the exact same sizes the h5ad would (a few MB vs. ~1 GB, and it is
    already on disk next to celltype_weight.csv). See mandatory-rules.md sec.10.

    Returns (sizes_dict, messages). Keys are normalized to int where possible
    so they match the integer cluster ids in celltype_weight.csv.
    """
    msgs: list[str] = []
    if not annotation_csv or not Path(annotation_csv).exists():
        return {}, msgs
    try:
        ann = pd.read_csv(annotation_csv, index_col=0)
    except Exception as e:
        return {}, [f"[WARN] check 4: could not read annotation CSV "
                    f"{annotation_csv}: {e}"]
    # Pick the cluster column: explicit cluster_key, else first known default.
    key = None
    if cluster_key and cluster_key in ann.columns:
        key = cluster_key
    else:
        for cand in ("leiden", "louvain", "cluster", "Cluster"):
            if cand in ann.columns:
                key = cand
                break
    if key is None:
        return {}, [f"[WARN] check 4: no cluster column in {annotation_csv} "
                    f"(looked for {cluster_key!r}, leiden, louvain, cluster)."]
    counts = ann[key].astype(str).value_counts().to_dict()
    sizes = {int(k) if str(k).isdigit() else k: int(v) for k, v in counts.items()}
    return sizes, msgs


def check_weak_clusters(
    ct_path: Path,
    annotation_csv: Path | None = None,
    cluster_key: str | None = None,
    h5ad_path: Path | None = None,
) -> tuple[bool, list[str]]:
    msgs = []
    df = pd.read_csv(ct_path)
    top1 = df.loc[df.groupby("cluster")["init_weight"].idxmax()].copy()
    top1["init_weight"] = top1["init_weight"].astype(float)

    # Cluster sizes come from the lightweight annotation CSV (no h5ad re-read).
    cluster_sizes: dict = {}
    if annotation_csv:
        cluster_sizes, size_msgs = cluster_sizes_from_annotation(
            annotation_csv, cluster_key
        )
        msgs.extend(size_msgs)
    # Legacy fallback: only if no annotation CSV was available and an h5ad was
    # explicitly passed (kept for standalone/back-compat use; not used by the
    # in-process gate path).
    if not cluster_sizes and h5ad_path and Path(h5ad_path).exists():
        try:
            import scanpy as sc
            adata = sc.read(h5ad_path)
            for key in ("leiden", "louvain", "cluster", "Cluster"):
                if key in adata.obs.columns:
                    cluster_sizes = (
                        adata.obs[key].astype(str).value_counts().to_dict()
                    )
                    cluster_sizes = {
                        int(k) if k.isdigit() else k: v
                        for k, v in cluster_sizes.items()
                    }
                    break
        except Exception as e:
            msgs.append(f"[WARN] check 4: could not read h5ad for cluster "
                        f"sizes: {e}")

    weak = []
    for _, row in top1.iterrows():
        cid = row["cluster"]
        try:
            cid_int = int(cid)
        except (ValueError, TypeError):
            cid_int = cid
        size = cluster_sizes.get(cid_int, 0)
        if row["init_weight"] < THRESH_WEAK_INIT_WEIGHT and size > 200:
            weak.append((cid, row["init_weight"], size))
    msgs.append(f"[INFO] check 4: clusters with weak top-1 init_weight "
                f"(< {THRESH_WEAK_INIT_WEIGHT}) and >200 cells: {len(weak)}")
    if len(weak) > THRESH_ALLOW_WEAK_CLUSTERS:
        msgs.append(f"[FAIL] check 4: {len(weak)} clusters have near-zero "
                    f"confidence annotations (init_weight < "
                    f"{THRESH_WEAK_INIT_WEIGHT}, n_cells > 200).")
        return False, msgs
    return True, msgs


def check_empty_annotations(ct_path: Path) -> tuple[bool, list[str]]:
    df = pd.read_csv(ct_path)
    n_clusters_with_rows = df["cluster"].nunique()
    # Heuristic: if config declared N clusters but we got far fewer, warn.
    # Without access to config, just report.
    return True, [f"[INFO] check 5: {n_clusters_with_rows} clusters have "
                  "annotation rows in celltype_weight.csv"]


def run_gate(
    outdir,
    annotation_csv=None,
    cluster_key=None,
    log_path=None,
    h5ad_path=None,
    as_json=False,
):
    """Run all five quality checks for a finished run; return (ok, results).

    This is the importable core of the gate. ``run_pipeline.py`` calls it
    in-process right after writing the annotation artifact, so the gate adds
    ZERO extra h5ad reads (cluster sizes come from ``annotation_csv``, a few-MB
    CSV). The standalone ``main()`` below is a thin CLI wrapper over it.

    Parameters
    ----------
    outdir : str | Path
        Run output directory (holds celltype_weight.csv and xener.log).
    annotation_csv : str | Path | None
        The ``{dataset}_annotation.csv`` for cheap cluster-size lookup. If
        omitted, auto-detected as the first ``*_annotation.csv`` in ``outdir``.
    cluster_key : str | None
        Cluster column name in the annotation CSV (e.g. ``leiden``).
    log_path : str | Path | None
        Run log; auto-detected in ``outdir`` if omitted.
    h5ad_path : str | Path | None
        Legacy fallback for cluster sizes only when no annotation CSV exists.

    Returns
    -------
    (overall_ok: bool, results: dict)
    """
    outdir = Path(outdir)
    log_path = Path(log_path) if log_path else find_run_log(outdir)
    ct_path = outdir / "celltype_weight.csv"

    # Auto-detect the annotation CSV in outdir when not supplied.
    if annotation_csv is None:
        cands = sorted(outdir.glob("*_annotation.csv"))
        annotation_csv = cands[0] if cands else None
    annotation_csv = Path(annotation_csv) if annotation_csv else None
    h5ad_path = Path(h5ad_path) if h5ad_path else None

    results = {}
    overall_ok = True

    # Check 1+2: KG miss (parsed from the run log)
    if log_path and Path(log_path).exists():
        cov = parse_kg_miss_per_cluster(Path(log_path))
        ok, msgs = check_kg_miss(cov)
        results["kg_miss"] = {
            "per_cluster": {str(k): v for k, v in cov.items()},
            "messages": msgs,
            "ok": ok,
        }
        overall_ok = overall_ok and ok
    else:
        results["kg_miss"] = {
            "ok": False,
            "messages": [f"[FAIL] could not locate run log in {outdir}. "
                         "Pass --log explicitly if the log is elsewhere."],
        }
        overall_ok = False

    # Check 3: cell-type diversity
    ok, msgs = check_celltype_diversity(ct_path)
    results["celltype_diversity"] = {"messages": msgs, "ok": ok}
    overall_ok = overall_ok and ok

    # Check 4: weak clusters (sizes from annotation CSV, NOT the h5ad)
    ok, msgs = check_weak_clusters(ct_path, annotation_csv, cluster_key, h5ad_path)
    results["weak_clusters"] = {"messages": msgs, "ok": ok}
    overall_ok = overall_ok and ok

    # Check 5: empty annotations
    ok, msgs = check_empty_annotations(ct_path)
    results["empty_annotations"] = {"messages": msgs, "ok": ok}
    overall_ok = overall_ok and ok

    # Emit results (shared by the in-process and CLI callers).
    if as_json:
        print(json.dumps(results, indent=2))
    else:
        for section in ("kg_miss", "celltype_diversity", "weak_clusters",
                        "empty_annotations"):
            for m in results[section]["messages"]:
                print(m)
        print()
        if overall_ok:
            print("Quality gate PASSED.")
        else:
            print("Quality gate FAILED. Re-run with adjusted config "
                  "(typically widen model_species).")

    return overall_ok, results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="xener output directory")
    ap.add_argument("--log", default=None,
                    help="Path to run log (auto-detected if omitted)")
    ap.add_argument("--annotation-csv", default=None,
                    help="Lightweight {dataset}_annotation.csv for cluster "
                         "sizes (auto-detected in --outdir if omitted). "
                         "Preferred over --h5ad: avoids re-reading the h5ad.")
    ap.add_argument("--cluster-key", default=None,
                    help="Cluster column name in the annotation CSV (e.g. leiden)")
    ap.add_argument("--h5ad", default=None,
                    help="LEGACY fallback for cluster sizes, used only when no "
                         "annotation CSV is available. Prefer --annotation-csv.")
    ap.add_argument("--json", action="store_true",
                    help="Emit results as JSON")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        print(f"[FAIL] outdir does not exist: {outdir}", file=sys.stderr)
        sys.exit(2)

    overall_ok, _ = run_gate(
        outdir=outdir,
        annotation_csv=args.annotation_csv,
        cluster_key=args.cluster_key,
        log_path=args.log,
        h5ad_path=args.h5ad,
        as_json=args.json,
    )
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
