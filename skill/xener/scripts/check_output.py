#!/usr/bin/env python3
"""Post-run quality gate for xener pipelines.

This script is invoked automatically at the end of run_pipeline.py, and
can also be run standalone:

    python scripts/check_output.py --outdir output/edf/

It reads celltype_weight.csv and the run log in --outdir, computes five
quality signals, and exits non-zero if any threshold is breached. The
thresholds are intentionally tight: a passing run is a run whose output
is likely to be biologically useful, not just a run that produced files.

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


def check_weak_clusters(ct_path: Path, h5ad_path: Path | None) -> tuple[bool, list[str]]:
    msgs = []
    df = pd.read_csv(ct_path)
    top1 = df.loc[df.groupby("cluster")["init_weight"].idxmax()].copy()
    top1["init_weight"] = top1["init_weight"].astype(float)

    # If h5ad provided, weight by cluster size
    cluster_sizes: dict[int, int] = {}
    if h5ad_path and h5ad_path.exists():
        try:
            import scanpy as sc
            adata = sc.read(h5ad_path)
            # cluster_key is the first column name whose values look like
            # the cluster ids in celltype_weight.csv. Best effort: try
            # the most common default 'leiden' first.
            for key in ("leiden", "louvain", "cluster", "Cluster"):
                if key in adata.obs.columns:
                    cluster_sizes = (
                        adata.obs[key].astype(str).value_counts().to_dict()
                    )
                    # normalize keys to int if possible
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="xener output directory")
    ap.add_argument("--log", default=None,
                    help="Path to run log (auto-detected if omitted)")
    ap.add_argument("--h5ad", default=None,
                    help="Original h5ad (used to weight weak-cluster check)")
    ap.add_argument("--json", action="store_true",
                    help="Emit results as JSON")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        print(f"[FAIL] outdir does not exist: {outdir}", file=sys.stderr)
        sys.exit(2)

    log_path = Path(args.log) if args.log else find_run_log(outdir)
    ct_path = outdir / "celltype_weight.csv"
    h5ad_path = Path(args.h5ad) if args.h5ad else None

    results = {}
    overall_ok = True

    # Check 1+2: KG miss
    if log_path and log_path.exists():
        cov = parse_kg_miss_per_cluster(log_path)
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

    # Check 4: weak clusters
    ok, msgs = check_weak_clusters(ct_path, h5ad_path)
    results["weak_clusters"] = {"messages": msgs, "ok": ok}
    overall_ok = overall_ok and ok

    # Check 5: empty annotations
    ok, msgs = check_empty_annotations(ct_path)
    results["empty_annotations"] = {"messages": msgs, "ok": ok}
    overall_ok = overall_ok and ok

    if args.json:
        print(json.dumps(results, indent=2))
    else:
        for section in ("kg_miss", "celltype_diversity", "weak_clusters",
                        "empty_annotations"):
            for m in results[section]["messages"]:
                print(m)
        print()
        if overall_ok:
            print("Quality gate PASSED.")
            sys.exit(0)
        else:
            print("Quality gate FAILED. Re-run with adjusted config "
                  "(typically widen model_species).")
            sys.exit(1)


if __name__ == "__main__":
    main()
