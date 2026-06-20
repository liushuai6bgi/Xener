#!/usr/bin/env python3
"""Post-run quality gate for xener pipelines (BASELINE-RELATIVE).

This script is invoked automatically at the end of run_pipeline.py (in-process,
via ``run_gate()`` -- no subprocess, no h5ad re-read), and can also be run
standalone:

    python scripts/check_output.py --outdir output/edf/
    python scripts/check_output.py --outdir output/edf_v2/ --baseline output/edf/

It reads celltype_weight.csv and the run log in --outdir and computes a small
set of quality SIGNALS, then decides PASS/FAIL with this rule:

    PASS = structurally sound
           AND ( absolute targets met OR not worsened vs --baseline )

In words ("improvement is enough"):

  * Absolute thresholds (mean KG miss <= 30%, etc.) are **advisory targets**,
    not hard failures. Being above a target emits a [WARN], never a [FAIL].
  * The only HARD failures are
      (a) STRUCTURAL breakage -- celltype_weight.csv missing/empty, no cluster
          has annotation rows, or the run log is absent/unparseable; and
      (b) a measured REGRESSION of a key signal versus a supplied --baseline
          (you changed the config and made the signals worse).
  * With no --baseline, a structurally-sound run PASSES even if it is below an
    advisory target -- there is nothing to improve against yet. Run it again
    with a soft-lever change and pass the first run as --baseline to prove the
    change improved the signals.

Each run writes its metrics to ``<outdir>/gate_metrics.json`` so a later run
can be compared against it via ``--baseline <that outdir or json>``.

This honours mandatory-rules.md sec.11: the intended response to weak signals is
to tune SOFT levers (model_species -- add a KG-rich relative, even a more
distant one; BLAST thresholds; top_num; ...) and keep the config that IMPROVES
the signals. HARD constraints -- the sample's real ``organ``, the cluster_key,
the input data, anything the user fixed -- are NEVER changed to move a signal.
The gate never reads or touches ``organ``; it cannot tempt you into the organ
trap.

Cluster sizes are no longer needed (the old absolute weak-cluster check on
``init_weight`` has been removed -- it was scale-dependent and brittle at its
boundary). ``--annotation-csv`` / ``--h5ad`` are accepted for backward
compatibility but unused.

Skill context: this is the mandatory Step 5.5 quality gate defined in
references/workflows/self-tuning-protocol.md. Do not declare an xener run done
without running it; a structurally-broken run, or a soft-lever change that
regresses the signals, is still a failure.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


# Advisory targets. Being above these emits a [WARN], not a [FAIL]; they are
# goals the agent should tune SOFT levers toward, not gate-blocking thresholds.
TARGET_MEAN_KG_MISS = 0.30
TARGET_TAIL_KG_MISS = 0.80
TARGET_TAIL_KG_MISS_FRAC = 0.05
TARGET_MIN_UNIQUE_CELLTYPES = 5
MIN_CLUSTERS_FOR_DIVERSITY = 10

# A baseline comparison must move a fractional signal by more than this margin
# to count as improvement/regression (suppresses floating-point and tie noise).
REGRESSION_EPS_FRAC = 0.01  # 1 percentage point on miss fractions

METRICS_FILENAME = "gate_metrics.json"


def find_run_log(outdir: Path) -> Path | None:
    """Locate the most recent xener run log in outdir."""
    candidates = [
        outdir / "xener.log",
        outdir / "run.log",
        outdir.parent / "xener.log",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def parse_kg_miss_per_cluster(log_path: Path, n_expected_clusters: int | None = None) -> dict[int, float]:
    """Parse `total X% homolos of organ[...] not in kg` lines from the log.

    Returns dict mapping cluster_id -> miss fraction (0.0 to 1.0).
    Cluster IDs come from the most recent "processing cluster N" header.

    Parameters
    ----------
    log_path:
        Path to the xener run log.
    n_expected_clusters:
        Number of clusters in the dataset (from celltype_weight.csv).
        Pass this to enable a sanity check when KG miss count mismatches.
    """
    text = log_path.read_text(errors="replace")
    cov: dict[int, float] = {}
    current_cluster: int | None = None
    in_annotation = False
    unmatched_miss_lines = 0
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
        if m:
            if current_cluster is not None:
                cov[current_cluster] = float(m.group(1)) / 100.0
            else:
                unmatched_miss_lines += 1

    if unmatched_miss_lines:
        print(f"[WARN] {unmatched_miss_lines} KG-miss line(s) could not be "
              "assigned to a cluster (logged before any 'processing cluster' "
              "header).", file=sys.stderr)

    if n_expected_clusters is not None and cov and len(cov) != n_expected_clusters:
        print(f"[WARN] KG-miss entries parsed for {len(cov)} cluster(s), but "
              f"celltype_weight.csv has {n_expected_clusters} cluster(s). The "
              "log format may have changed, or some clusters have no KG-miss "
              "output.", file=sys.stderr)

    return cov


def compute_metrics(cov: dict[int, float], ct_path: Path) -> dict:
    """Compute the quality signals. Values are None when not computable.

    Returns a dict with:
      mean_kg_miss, tail_kg_miss_frac           -- coverage (from the log)
      n_unique_top1_celltypes, n_clusters       -- diversity (from celltype_weight.csv)
      n_clusters_with_rows                       -- structure
    """
    m: dict = {
        "mean_kg_miss": None,
        "tail_kg_miss_frac": None,
        "n_unique_top1_celltypes": None,
        "n_clusters": None,
        "n_clusters_with_rows": None,
    }
    if cov:
        n = len(cov)
        m["mean_kg_miss"] = sum(cov.values()) / n
        m["tail_kg_miss_frac"] = sum(1 for v in cov.values() if v > TARGET_TAIL_KG_MISS) / n
    if ct_path.exists():
        df = pd.read_csv(ct_path)
        if {"cluster", "celltype", "init_weight"}.issubset(df.columns):
            top1 = df.loc[df.groupby("cluster")["init_weight"].idxmax()]
            m["n_unique_top1_celltypes"] = int(top1["celltype"].nunique())
            m["n_clusters"] = int(top1["cluster"].nunique())
        if "cluster" in df.columns:
            m["n_clusters_with_rows"] = int(df["cluster"].nunique())
    return m


def load_baseline(baseline: str | Path | None) -> tuple[dict | None, list[str]]:
    """Resolve --baseline (a gate_metrics.json file or an outdir holding one)."""
    if not baseline:
        return None, []
    p = Path(baseline)
    if p.is_dir():
        p = p / METRICS_FILENAME
    if not p.exists():
        return None, [f"[WARN] baseline: {p} not found; comparison skipped."]
    try:
        return json.loads(p.read_text()), [f"[INFO] baseline: comparing against {p}"]
    except Exception as e:  # noqa: BLE001
        return None, [f"[WARN] baseline: could not read {p}: {e}; comparison skipped."]


def check_structure(ct_path: Path, log_path: Path | None,
                    cov: dict[int, float], metrics: dict) -> tuple[bool, list[str]]:
    """Hard structural floor: the run must have produced parseable annotations."""
    msgs: list[str] = []
    ok = True
    if not ct_path.exists():
        return False, [f"[FAIL] structure: celltype_weight.csv not found at {ct_path}"]
    if metrics["n_clusters_with_rows"] in (None, 0):
        ok = False
        msgs.append("[FAIL] structure: no cluster has annotation rows in "
                    "celltype_weight.csv (empty annotation).")
    else:
        msgs.append(f"[INFO] structure: {metrics['n_clusters_with_rows']} clusters "
                    "have annotation rows.")
    if not (log_path and Path(log_path).exists()):
        ok = False
        msgs.append("[FAIL] structure: could not locate run log; KG-miss signal "
                    "is unobservable. Pass --log if it lives elsewhere.")
    elif not cov:
        msgs.append("[WARN] structure: run log found but no 'total X% homolos ... "
                    "not in kg' lines parsed; KG-miss signal unavailable.")
    return ok, msgs


def report_targets(metrics: dict) -> tuple[bool, list[str]]:
    """Advisory comparison against absolute targets. Never FAILs; sets met flag."""
    msgs: list[str] = []
    met = True

    mm = metrics["mean_kg_miss"]
    tf = metrics["tail_kg_miss_frac"]
    if mm is not None:
        line = f"[INFO] signal: mean KG miss = {mm:.1%} (target <= {TARGET_MEAN_KG_MISS:.0%})"
        if mm > TARGET_MEAN_KG_MISS:
            met = False
            line = (f"[WARN] signal: mean KG miss {mm:.1%} ABOVE target "
                    f"{TARGET_MEAN_KG_MISS:.0%}. Soft fix: add a KG-rich relative to "
                    "model_species (closest first, or a more distant well-covered "
                    "species) and re-run with this run as --baseline. Never change a "
                    "confirmed organ to lower it (mandatory-rules.md sec.11).")
        msgs.append(line)
    if tf is not None:
        line = (f"[INFO] signal: clusters with KG miss > {TARGET_TAIL_KG_MISS:.0%} "
                f"= {tf:.1%} (target <= {TARGET_TAIL_KG_MISS_FRAC:.0%})")
        if tf > TARGET_TAIL_KG_MISS_FRAC:
            met = False
            line = (f"[WARN] signal: {tf:.1%} of clusters have severe KG miss "
                    f"(> {TARGET_TAIL_KG_MISS:.0%}); target <= {TARGET_TAIL_KG_MISS_FRAC:.0%}.")
        msgs.append(line)

    nu = metrics["n_unique_top1_celltypes"]
    nc = metrics["n_clusters"]
    if nu is not None and nc is not None:
        line = f"[INFO] signal: {nu} unique top-1 cell types across {nc} clusters"
        if nc > MIN_CLUSTERS_FOR_DIVERSITY and nu < TARGET_MIN_UNIQUE_CELLTYPES:
            met = False
            line = (f"[WARN] signal: only {nu} unique top-1 cell types for {nc} "
                    f"clusters (target >= {TARGET_MIN_UNIQUE_CELLTYPES}); model_species "
                    "may be too narrow for the organ.")
        msgs.append(line)
    return met, msgs


def compare_baseline(metrics: dict, base: dict) -> tuple[bool, bool, list[str]]:
    """Compare current metrics to a baseline.

    Returns (improved, worsened, messages). A signal counts only when present in
    both. Fractional signals use REGRESSION_EPS_FRAC; the diversity count uses a
    strict integer change.
    """
    msgs: list[str] = []
    improved = False
    worsened = False

    def cmp_frac(key, label, lower_is_better=True):
        nonlocal improved, worsened
        cur, old = metrics.get(key), base.get(key)
        if cur is None or old is None:
            return
        delta = cur - old
        if abs(delta) <= REGRESSION_EPS_FRAC:
            msgs.append(f"[INFO] baseline: {label} {old:.1%} -> {cur:.1%} (flat)")
            return
        better = (delta < 0) if lower_is_better else (delta > 0)
        arrow = "improved" if better else "REGRESSED"
        if better:
            improved = True
        else:
            worsened = True
        msgs.append(f"[INFO] baseline: {label} {old:.1%} -> {cur:.1%} ({arrow})")

    cmp_frac("mean_kg_miss", "mean KG miss", lower_is_better=True)
    cmp_frac("tail_kg_miss_frac", "severe-miss cluster frac", lower_is_better=True)

    cur, old = metrics.get("n_unique_top1_celltypes"), base.get("n_unique_top1_celltypes")
    if cur is not None and old is not None:
        if cur > old:
            improved = True
            msgs.append(f"[INFO] baseline: unique top-1 types {old} -> {cur} (improved)")
        elif cur < old:
            worsened = True
            msgs.append(f"[INFO] baseline: unique top-1 types {old} -> {cur} (REGRESSED)")
        else:
            msgs.append(f"[INFO] baseline: unique top-1 types {old} -> {cur} (flat)")
    return improved, worsened, msgs


def run_gate(
    outdir,
    annotation_csv=None,   # accepted for backward-compat; unused (check 4 removed)
    cluster_key=None,      # accepted for backward-compat; unused
    log_path=None,
    h5ad_path=None,        # accepted for backward-compat; unused
    as_json=False,
    baseline=None,
):
    """Run the baseline-relative quality gate; return (ok, results).

    PASS = structurally sound AND (absolute targets met OR not worsened vs
    --baseline). Absolute thresholds are advisory. See the module docstring and
    references/workflows/self-tuning-protocol.md.

    Always writes ``<outdir>/gate_metrics.json`` so a later run can pass this
    one as ``--baseline``.
    """
    outdir = Path(outdir)
    log_path = Path(log_path) if log_path else find_run_log(outdir)
    ct_path = outdir / "celltype_weight.csv"

    # Probe celltype_weight.csv for expected cluster count, so
    # parse_kg_miss_per_cluster can sanity-check its results.
    n_clusters_expected: int | None = None
    if ct_path.exists():
        try:
            tmp = pd.read_csv(ct_path)
            if "cluster" in tmp.columns:
                n_clusters_expected = int(tmp["cluster"].nunique())
        except Exception:
            pass

    cov = (parse_kg_miss_per_cluster(Path(log_path), n_clusters_expected)
           if log_path and Path(log_path).exists() else {})
    metrics = compute_metrics(cov, ct_path)

    # Persist metrics for use as a future baseline (best-effort).
    try:
        (outdir / METRICS_FILENAME).write_text(json.dumps(metrics, indent=2))
    except Exception:  # noqa: BLE001
        pass

    structure_ok, struct_msgs = check_structure(ct_path, log_path, cov, metrics)
    met_targets, target_msgs = report_targets(metrics)

    base, base_load_msgs = load_baseline(baseline)
    improved = worsened = False
    cmp_msgs: list[str] = []
    if base is not None:
        improved, worsened, cmp_msgs = compare_baseline(metrics, base)

    # Decision: structural floor first, then improvement-is-enough.
    if not structure_ok:
        overall_ok = False
        verdict = "FAILED (structural breakage)"
    elif met_targets:
        overall_ok = True
        verdict = "PASSED (absolute targets met)"
    elif base is not None and worsened:
        overall_ok = False
        verdict = "FAILED (regression vs baseline)"
    elif base is not None:
        overall_ok = True
        verdict = ("PASSED (improved vs baseline)" if improved
                   else "PASSED (held vs baseline; below target but not worsened)")
    else:
        overall_ok = True
        verdict = ("PASSED (below target; no baseline to judge improvement -- "
                   "re-run with a soft-lever change and --baseline to confirm a gain)")

    results = {
        "metrics": metrics,
        "per_cluster_kg_miss": {str(k): v for k, v in cov.items()},
        "met_targets": met_targets,
        "improved_vs_baseline": improved,
        "worsened_vs_baseline": worsened,
        "ok": overall_ok,
        "verdict": verdict,
        "messages": struct_msgs + target_msgs + base_load_msgs + cmp_msgs,
    }

    if as_json:
        print(json.dumps(results, indent=2))
    else:
        for m in results["messages"]:
            print(m)
        print()
        print(f"Quality gate {'PASSED' if overall_ok else 'FAILED'}: {verdict}.")
        if not overall_ok:
            print("Fix: if structural, inspect the pipeline output/log. If a "
                  "regression, revert the soft-lever change. Never change a HARD "
                  "constraint (organ, cluster_key, inputs) to move a signal "
                  "(mandatory-rules.md sec.11).")

    return overall_ok, results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="xener output directory")
    ap.add_argument("--baseline", default=None,
                    help="A prior run's gate_metrics.json (or its outdir) to "
                         "compare against. With a baseline, a structurally-sound "
                         "run that improved or held its signals PASSES even below "
                         "target; a run that regressed FAILS.")
    ap.add_argument("--log", default=None,
                    help="Path to run log (auto-detected if omitted)")
    ap.add_argument("--annotation-csv", default=None,
                    help="Accepted for backward compatibility; unused.")
    ap.add_argument("--cluster-key", default=None,
                    help="Accepted for backward compatibility; unused.")
    ap.add_argument("--h5ad", default=None,
                    help="Accepted for backward compatibility; unused.")
    ap.add_argument("--json", action="store_true", help="Emit results as JSON")
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
        baseline=args.baseline,
    )
    sys.exit(0 if overall_ok else 1)


if __name__ == "__main__":
    main()
