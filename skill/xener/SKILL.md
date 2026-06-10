---
name: xener
version: 0.1.10
description: |
  Use this skill whenever the user wants to annotate single-cell RNA-seq data
  with cell type labels, especially in cross-species scenarios where the
  target species lacks a well-annotated reference. Xener maps marker genes
  from a non-model species to one or more model species via BLAST homology,
  then propagates weights through a curated knowledge graph to predict cell
  types per cluster. Supports sub-cluster refinement with Moran's I gene
  filtering for splitting "mixed" clusters into subtypes.

  Trigger on any of: single-cell annotation, cell type prediction, cell type
  labeling, cross-species gene mapping, h5ad cell type analysis, marker gene
  annotation for clusters, scRNA-seq cluster annotation, "annotate this
  dataset", "what cell types are in my data", "label clusters in this
  single-cell data", or any mention of BLAST + knowledge-graph based
  annotation. Make sure to invoke this skill even when the user just says
  "annotate" in a single-cell context — do not fall back to generic Python
  scripts that import xener directly.

allowed-tools: Bash(python scripts/*), Read, Write, Edit, Glob, Grep
metadata:
  category: bioinformatics
  language: python
  dependencies:
    - xener>=0.1.10
    - scanpy
    - pandas
    - pyyaml
  domain: single-cell RNA-seq, cross-species cell type annotation
---

# Xener — Cross-species Cell Type Annotation

## What it does
Xener annotates cell types in single-cell RNA-seq data from non-model species.
Given an `.h5ad` file and a FASTA file, it:

1. Identifies marker genes per cluster (Step 1)
2. Calculates gene weights from differential expression (Step 2)
3. Maps them to model species via BLAST homology (Step 3)
4. Selects top-k genes per cluster (Step 4)
5. Propagates gene→celltype weights through a curated knowledge graph (Step 5)
6. Optionally refines "mixed" clusters into subtypes (Refinement)

## When to use it
- Target organism has no high-quality reference atlas
- Need cross-species transfer (e.g., annotate bamboo using Arabidopsis + Oryza)
- Want to refine a "mixed" cluster into subtypes
- Have an h5ad + fasta pair and want a YAML-driven annotation pipeline

## Autonomous agent mode (preferred)

The agent should make **most parameter decisions on its own** based on
dataset inspection and result quality. Read
`references/autonomous-decision-making.md` for the full decision matrix
and heuristics. The high-level flow is:

1. **Inspect** the h5ad (`workflows/inspection.md`) to detect species,
   cluster_key, organ hints
2. **Decide** model_species, organ, top_num, BLAST thresholds, mode,
   decay_factor — see `workflows/species-selection.md` for species
   heuristics and `autonomous-decision-making.md` for the rest
3. **Compose a config.yaml** with reasoning
4. **Confirm biological choices** (species, organ) with the user —
   propose 1-3 candidates with reasoning, not open-ended questions
5. **Run** the pipeline (`workflows/full-pipeline.md`)
6. **Run the mandatory post-run quality gate** (`scripts/check_output.py`,
   invoked automatically by `run_pipeline.py`). If it fails, the
   pipeline exits non-zero. The agent must diagnose the failure
   (typically: widen `model_species` to include the target species
   itself when it is a model organism) and re-run from Step 3.
   See `workflows/self-tuning-protocol.md` for the gate logic and
   `workflows/species-selection.md` for the most common fix. When
   the gate fails, **also consult `references/log-interpretation.md`**
   for the stage-by-stage grep-driven diagnosis — the gate identifies
   the class of failure, the log reference identifies the exact
   failure inside the class.
7. **Refine** mixed clusters (`workflows/refinement.md`) — in autonomous
   mode (complete-annotation), refine **every** cluster where the
   top-2 distinct cell types have an `init_weight` ratio > 0.5, in
   order of descending ratio. The previous "up to 3" guidance is
   withdrawn; see `mandatory-rules.md` §8 and
   `autonomous-decision-making.md` §"Completeness vs. demonstration"
   for the full reasoning. In manual mode, present the eligible
   list and ask the user which to refine before running.
8. **Visualize** (`workflows/refinement.md` Step D) — run `plot_umap.py
   --mode overview` once, then `--mode refine` for **every** cluster
   that split into >1 subtype. Visualization coverage must match
   refinement coverage; plotting only a "representative" few is the
   same demo-mode shortcut withdrawn in step 7. Clusters that refined
   to a single uniform label have nothing to plot and are skipped.

## Quick start (manual mode)

If the user wants full control, the manual flow is:

1. **Install**: `python scripts/install.py` (or `pip install xener`)
2. **Write a config.yaml** (see `references/config-schema.md`)
3. **Validate** with `references/workflows/config-validation.md`
4. **Run the full pipeline**:
   ```bash
   python scripts/run_pipeline.py --config config.yaml
   ```
5. **Refine mixed clusters** (optional): see `references/workflows/refinement.md`
6. **Visualize** the annotation and/or refinement as UMAP PNGs. The input can be
   the lightweight `{dataset}_annotation.csv` (UMAP coords + labels) — no need to
   reload the multi-GB h5ad:
   ```bash
   # Overview: cluster | xener | xener_max | xener_refine in one figure
   python scripts/plot_umap.py --input output/edf/edf_annotation.csv \
       --mode overview --cluster-key leiden --outdir output/
   # Whole-dataset annotation (2 panels)
   python scripts/plot_umap.py --input output/edf/edf_annotation.csv \
       --mode annotation --cluster-key leiden --outdir output/
   # One refined cluster (repeat for EVERY cluster that split; see below)
   python scripts/plot_umap.py --input output/edf/edf_annotation.csv \
       --mode refine --cluster-key leiden --cluster-id 4 \
       --refine-key xener_refine --outdir output/
   ```
   In complete-annotation mode, plot a `--mode refine` figure for **every**
   cluster that split into >1 subtype — not a representative subset.
   Visualization coverage must match refinement coverage; see
   `workflows/refinement.md` → Step D.

## Custom initialization (own KG / BLAST database)

By default Xener uses the public cloud Knowledge Graph and a bundled BLAST
database — no setup needed. If the user wants their **own** infrastructure
(on-prem Neo4j KG, a locally-built BLAST DB, or a BLASTP cache), they can
supply an **init-config** that says *where* Xener gets its data — distinct from
the run `config.yaml`, which says *what* to annotate.

Validate it first, then pass the same `--init-config` to the pipeline (every
Xener-building script accepts the flag):

```bash
python scripts/init_xener.py --init-config xener-init.yaml      # validate
python scripts/run_pipeline.py --config config.yaml --init-config xener-init.yaml
```

All init keys are optional and may instead be inlined into `config.yaml`. See
`references/workflows/initialization.md` and
`examples/xener-init.example.yaml`.

## Mandatory rules (read `references/mandatory-rules.md` first)

The full list of do's and don'ts is in `references/mandatory-rules.md`. The
three most important constraints:

- **Use ONLY scripts in `scripts/`. Never `import xener` directly.**
- **Never pass config via stdin.** Always use a config file path.
- **Always confirm `model_species` and `organ` with the user before running.**

## Reference docs (load on demand)

| File | When to read |
|------|--------------|
| `references/mandatory-rules.md` | Before any operation — full list of constraints |
| `references/autonomous-decision-making.md` | **Start here in autonomous mode** — decision matrix |
| `references/config-schema.md` | When writing or validating a YAML config |
| `references/parameters.md` | When tuning CLI parameters (e.g., pident threshold, top-k) |
| `references/workflows/inspection.md` | **First step in autonomous mode** — inspect h5ad |
| `references/workflows/initialization.md` | When the user wants a custom KG / BLAST database (`--init-config`) instead of the cloud defaults |
| `references/workflows/species-selection.md` | When picking model_species for a target |
| `references/workflows/full-pipeline.md` | When running end-to-end with `run_pipeline.py` |
| `references/workflows/step-by-step.md` | When iterating on individual steps |
| `references/workflows/config-validation.md` | Before running the pipeline, to verify config |
| `references/workflows/refinement.md` | When splitting a mixed cluster into subtypes |
| `references/workflows/self-tuning-protocol.md` | **After initial run** — auto-adjust parameters |
| `references/log-interpretation.md` | **When a run fails or output looks wrong** — grep-driven diagnosis from the run log |
| `references/output-files.md` | When interpreting `.csv` outputs |
| `references/troubleshooting.md` | When a step fails or produces unexpected results |

## Visualization

After Step 5 (and optionally after refinement), plot the result with
`scripts/plot_umap.py`:

| Mode | Output | Description |
|------|--------|-------------|
| `annotation` | `umap_annotation.png` | Side-by-side: cluster labels + Xener annotation |
| `overview` | `umap_overview.png` | One figure, 4 panels: cluster, `xener`, `xener_max`, `xener_refine` (unrefined cells gray) |
| `refine --cluster-id N` | `umap_refine_cluster_N.png` | Cluster N highlighted + refined subtypes |

In complete-annotation mode, run `overview` once, then `refine` for **every**
cluster that split into >1 subtype (not a representative subset). The column
meanings (`xener` vs `xener_max` vs `xener_refine`, and why `xener_refine` is
empty for unrefined cells) are in `references/output-files.md`.

## Examples

See `examples/config.example.yaml` for a complete run config (what to
annotate), and `examples/xener-init.example.yaml` for an init config (custom
KG / BLAST database — only needed for non-default infrastructure).
