# Post-Run Quality Gate (REQUIRED) — baseline-relative

**This is a mandatory step, not optional polish.** Every xener run, whether
invoked by `run_pipeline.py` or step-by-step, MUST pass through the quality
gate before the run is declared "done".

The gate is **baseline-relative** ("improvement is enough"). It does **not**
hard-fail a run merely for being above an absolute threshold. Its rule is:

```
PASS = structurally sound
       AND ( absolute targets met OR not worsened vs --baseline )
```

- **Absolute targets** (mean KG miss ≤ 30%, severe-miss tail ≤ 5%, ≥ 5 unique
  top-1 cell types) are **advisory**. Being above a target prints a `[WARN]`,
  never a `[FAIL]`. They are goals you tune SOFT levers toward, not blockers.
- The only **HARD failures** are:
  1. **Structural breakage** — `celltype_weight.csv` missing/empty, no cluster
     has annotation rows, or the run log is absent/unparseable.
  2. **Regression vs `--baseline`** — you changed the config and a key signal
     got measurably worse (mean KG miss or severe-miss tail up by > 1pp, or
     unique top-1 types down).
- With **no `--baseline`**, a structurally-sound run PASSES even below target —
  there is nothing yet to improve against. To *prove* a soft-lever change
  helped, re-run with the previous run's outdir as `--baseline`.

> **The old absolute weak-cluster check on `init_weight` (`< 50` for n_cells >
> 200) has been REMOVED.** It was scale-dependent (`init_weight` is an
> unnormalized sum that ranged ~14 to ~1280 across clusters in one real
> dataset) and brittle at its boundary (a cluster at 50.5 flipping to 49.8 from
> a trivial perturbation). Confidence is now judged *relatively* (did a
> soft-lever change improve the signals?), not by an absolute cutoff. Note:
> `init_weight` is **still** used for *refinement* eligibility (top-2 ratio);
> that is a different, scale-invariant ratio and is unaffected.

## Why a gate at all

A run can complete with exit code 0 and still produce biologically useless
annotations — most often because the KG had few edges for the chosen
`model_species` + `organ`, so homologs map to nothing and weights are thin.
The gate surfaces that as a `[WARN]` and, crucially, lets you **measure whether
a fix helped** by comparing against a baseline.

## Signals the gate computes

Run `scripts/check_output.py --outdir <outdir> [--baseline <prev_outdir>]`
(or let `run_pipeline.py` call it in-process).

| Signal | Source | Advisory target | Role |
|--------|--------|-----------------|------|
| Mean KG miss (per cluster) | `total X% homolos of organ[...] not in kg` in the log | ≤ 0.30 | primary, comparable |
| Severe-miss tail | same | ≤ 5% of clusters > 0.80 | comparable |
| Unique top-1 cell types | `celltype_weight.csv` | ≥ 5 (when > 10 clusters) | comparable |
| Annotation rows present | `celltype_weight.csv` | every cluster | **structural floor** |

Each run writes these to `<outdir>/gate_metrics.json` so a later run can be
compared against it.

## The intended response to a `[WARN]` (mandatory-rules.md §11)

A high mean KG miss is almost always `model_species` being too narrow / shallow
for the organ. The fix is a **SOFT lever**, and the gate now *rewards* it:

1. **Add a KG-present relative to `model_species`** — closest first; or, when
   close relatives are absent or KG-shallow for the organ, a **more distant but
   well-covered species** (e.g. Arabidopsis for a plant stem). `model_species`
   has no hard phylogenetic boundary. Keep the close relatives in the list.
2. Re-run from Step 3 (reuse `marker_gene.csv` / `marker_weight.csv`), passing
   the **previous run's outdir as `--baseline`**.
3. If the signal **improved (or held)** → the gate PASSES; adopt the new
   config. If it **regressed** → the gate FAILS; revert.

Other soft levers: relax BLAST (`pident`/`evalue`/`bitscore`), raise `top_num`,
adjust `threshold`/`decay_factor`/`mode`, restrict `candidate_annotation`.

> **The HARD constraints are never touched to move a signal.** The gate does
> not read `organ` and cannot be passed by changing it. Switching a stem sample
> to Root, or to `None`/`Unknown` just to drop the filter, is the **organ
> trap** — forbidden by `mandatory-rules.md` §11 even though it would lower KG
> miss. If, after adding the best available species (close *and* distant
> well-covered), the mean KG miss is still above target but did not regress, the
> run **PASSES** under the correct organ; the residual is shallow KG coverage
> for that organ, documented in `autonomous_log.md`, not a config bug to chase.

## What the gate looks like in practice

```bash
# First run (no baseline): structurally sound, below target -> PASS (advisory).
$ python scripts/check_output.py --outdir output/edf_v1/
[INFO] structure: 36 clusters have annotation rows.
[WARN] signal: mean KG miss 47.4% ABOVE target 30%. Soft fix: add a KG-rich
       relative to model_species ... Never change a confirmed organ (§11).
[INFO] signal: 3 unique top-1 cell types across 36 clusters
Quality gate PASSED: PASSED (below target; no baseline to judge improvement).

# Re-run after adding a KG-rich species, compared to the first run.
$ python scripts/check_output.py --outdir output/edf_v2/ --baseline output/edf_v1/
[INFO] baseline: mean KG miss 47.4% -> 31.2% (improved)
[INFO] baseline: unique top-1 types 3 -> 9 (improved)
Quality gate PASSED: PASSED (improved vs baseline).
```

## Automatic invocation

`run_pipeline.py` calls the gate at the end, in-process. Pass `--baseline
<prev_outdir>` to `run_pipeline.py` to forward it to the gate. The pipeline
exits non-zero only on a **structural** failure or a **regression** vs the
baseline — not for being above an advisory target.

## Manual re-run shortcut (improve, then compare)

```bash
# v1 (baseline): closest relatives only.
python scripts/run_pipeline.py --config config_v1.yaml --init-config init.yaml

# v2: add a KG-rich (possibly distant) species, compare to v1.
python scripts/run_pipeline.py --config config_v2.yaml --init-config init.yaml \
    --baseline output/edf_v1
# Gate PASSES iff v2 improved-or-held vs v1 (and is structurally sound).
```

Or run the steps individually (reuse `marker_weight.csv`) and call
`check_output.py --outdir output/edf_v2 --baseline output/edf_v1` at the end.

## Legacy "self-tuning" guidance

The original heuristics (e.g. "if median init_weight is low, re-run with
`--threshold null --decay-factor 0.5`") still apply as **second-line soft-lever
adjustments**. Apply one, re-run with the prior run as `--baseline`, and keep
the change only if the gate confirms the signals improved or held.
