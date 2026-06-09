# Output Files

The pipeline writes the following files to `outdir/`.

| File | Source step | Description |
|------|-------------|-------------|
| `marker_gene.csv` | Step 1 | Marker genes per cluster |
| `marker_weight.csv` | Step 2 | Weighted marker genes |
| `blastp_{species}.csv` | Step 3 | BLASTP alignment results (one per model species) |
| `gene_homolo_weight.csv` | Step 3 | Homology-mapped genes with weights |
| `topk_markers.csv` | Step 4 | Top-k genes per cluster |
| `celltype_weight.csv` | Step 5 | All predicted cell types per cluster |
| `debug_params.yaml` | Pipeline | Actual parameters used in each step (for reproducibility) |
| `config.yaml` | Pipeline | Copy of the input config (from `run_from_yaml` only) |
| `annotation/cluster_{id}_gene2celltype.xml` | Step 5 | Per-cluster annotation graph |
| `refine_suggestions.json` | `suggest_refine.py` | Suggested clusters for refinement |
| `xener.log` | Pipeline | Mirror of xener's stdout (created by `run_pipeline.py`); consumed by `scripts/check_output.py` for the mandatory post-run quality gate |
| `umap_annotation.png` | `plot_umap.py --mode annotation` | Side-by-side UMAP: cluster + Xener |
| `umap_refine_cluster_{N}.png` | `plot_umap.py --mode refine` | UMAP highlighting cluster N + refine |
| `refine_output/refined_{cluster_id}.csv` | `refine_cluster.py` | Refined annotation DataFrame (CSV) |
| `refine_output/refined_{cluster_id}_gene2celltype.gexf` | `refine_cluster.py` | Gene→homolo→celltype graph (Gephi/Cytoscape/networkx) |

## `debug_params.yaml` structure

This file records the actual parameter values used in each key step.
**Always inspect it after a run** to confirm the pipeline used the values
you intended (vs. the package defaults):

```yaml
cell_annotation:
  decay_factor: 0.7
  mode: path
  organ: leaf
  threshold: null
get_gene_weight:
  marker_weight_method: prod
get_topk_gene:
  multihomolo: true
  top_num: 30
mapping:
  bitscore: 200
  evalue: 0.05
  homolo_weight_key: pident
  model_species:
  - Oryza_sativa
  pident: 60
```

## Intermediates and checkpointing

Step scripts save intermediate files (`.csv`) so a step can be re-run
independently with new parameters without re-running earlier steps.
For example, to re-run step 3 with a different `--pident` threshold
without re-running step 1 and step 2, simply re-invoke
`step3_mapping.py` with the same `marker_weight.csv` as input.

## What "good" output looks like (sanity-check after a run)

The pipeline can complete cleanly and still produce biologically
useless output. Before declaring a run done, verify these signals
against `celltype_weight.csv` and the run log. The post-run gate
(`scripts/check_output.py`) does this automatically; the table below
is for manual inspection.

| Signal | Healthy range | Red flag | What to do |
|--------|---------------|----------|------------|
| Mean `total X% homolos of organ[...] not in kg` across clusters | < 30% | > 30% | Add a more well-annotated `model_species` (often the target species itself) |
| Tail: clusters with > 80% KG miss | < 5% of clusters | > 5% | Same as above |
| Number of unique top-1 cell types across clusters | 8-15+ for typical plant/animal atlases | < 5 unique types for > 10 clusters | Same as above |
| Per-cluster top-1 `init_weight` | > 50 for n_cells > 200 | One or more clusters near zero | Re-run step 5 with `--mode node` for those clusters, or trigger refinement |
| Clusters with no `celltype_weight` rows | 0 | > 0 | Investigate: organ filter, BLAST stringency, or KG coverage |
