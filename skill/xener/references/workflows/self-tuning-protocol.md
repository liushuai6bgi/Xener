# Post-Run Quality Gate (REQUIRED)

**This is a mandatory step, not optional polish.** Every xener run,
whether invoked by `run_pipeline.py` or step-by-step, MUST pass through
the quality gate below before the run is declared "done". A run that
produces output files but fails these checks is a failed run.

## Why this is mandatory

A run can complete with exit code 0 and still produce biologically
useless annotations:

- The pipeline writes `celltype_weight.csv`, but the KG may have had
  no edges for the chosen model_species + organ combination, so most
  cell types are "unknown" or the same label is repeated.
- The pipeline writes `gene_homolo_weight.csv`, but BLAST was too
  stringent (or pointed at a species with poor annotation depth), so
  few homologs survived.
- `init_weight` may look "high" for some clusters simply because a
  few genes voted overwhelmingly, while 70%+ of marker genes had no
  KG entry at all.

The signals for all of these are already in the run log and the
output zips. The agent must read them, not skip them.

## Required checks (must run after Step 5)

Run `scripts/check_output.py --outdir <outdir>` (or read the
diagnostic inline). It performs these five checks; the run is
**failed** if any threshold is breached.

| # | Signal | Source | Failure threshold |
|---|--------|--------|-------------------|
| 1 | Mean KG miss rate (per cluster) | `total X% homolos of organ[...] not in kg` in the run log | `mean > 0.30` |
| 2 | Tail of clusters with severe KG miss | same | `>5% of clusters with miss > 0.80` |
| 3 | Unique top-1 cell types across clusters | `celltype_weight.csv`, top-1 per cluster | `<5 unique types for >10 clusters` |
| 4 | Weak-confidence clusters | `celltype_weight.csv`, max per cluster | `>1 cluster with top-1 init_weight < 50` for n_cells > 200 |
| 5 | Empty annotations | same | `>0 clusters with no celltype_weight rows` |

The first three are the most important. A failure of (1) or (3) almost
always means `model_species` is too narrow for the chosen organ -- see
`workflows/species-selection.md` for the worked example.

## What the gate looks like in practice

```bash
$ python scripts/check_output.py --outdir output/edf/

[FAIL] check 1: mean KG miss 47.4% > 30%
        -> model_species too narrow for organ[Root].
        -> add target species itself (if it's a model organism)
           or a more well-annotated close relative.
[FAIL] check 3: only 3 unique top-1 cell types for 36 clusters
        -> expected 8-15+ for a root tip atlas.
        -> see species-selection.md worked example.

Quality gate FAILED. Re-run with adjusted config.
exit 1
```

## Automatic invocation

`run_pipeline.py` calls `check_output.py` automatically at the end.
If the gate fails, the pipeline exits non-zero and the agent must
adjust the config (most commonly `model_species`) and re-run from
Step 3 (re-using `marker_gene.csv` and `marker_weight.csv` from
the first run, since those don't depend on the BLAST database).

## Manual re-run shortcut

```bash
# Re-run only step 3-5 with new model_species, reusing earlier outputs:
python scripts/step3_mapping.py \
    --input output/edf/marker_weight.csv \
    --fasta abc.fasta \
    --species Arabidopsis_thaliana Brassica_rapa \
    --pident 60 --evalue 0.05 --bitscore 200 \
    --multihomolo \
    --outdir output/edf/

python scripts/step4_topk.py --input output/edf/gene_homolo_weight.csv \
    --top-num 20 --multihomolo --outdir output/edf/

python scripts/step5_annotate.py --input output/edf/topk_markers.csv \
    --organ Root --threshold null --mode path --decay-factor 0.7 \
    --outdir output/edf/

python scripts/check_output.py --outdir output/edf/   # must pass
```

## Legacy "self-tuning" guidance

The original self-tuning heuristics (e.g., "if median init_weight
< 0.3, re-run with --threshold null --decay-factor 0.5") still
apply **as second-line adjustments** after the post-run gate has
identified a real problem. They are no longer the primary
diagnostic -- the gate above is.
