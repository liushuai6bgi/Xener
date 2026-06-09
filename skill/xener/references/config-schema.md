# Config Schema (YAML)

A complete reference for the YAML configuration file passed to
`scripts/run_pipeline.py --config <path>`.

## Full example

```yaml
# Required fields
cluster_key: leiden
model_species:
  - Brassica_rapa
non_model_fasta: Arabidopsis_thaliana.fasta
non_model_h5ad: ERP132245.h5ad
organ: leaf
outdir: output/ERP132245

# Optional fields (defaults shown)
marker_weight_method: prod          # prod | sum
top_num: 30                         # Top N genes per cluster
homolo_weight_key: pident           # pident | evalue | bitscore
multihomolo: true                   # Keep multiple homologs per gene
decay_factor: 0.7                   # Weight decay for graph propagation
mode: path                          # node | path
threshold: null                     # Z-score threshold for filtering
candidate_annotation:               # Restrict cell types (optional)
  - type1
  - type2
```

## Field reference

| Field | Type | Required | Default | Description |
|-------|------|----------|---------|-------------|
| `cluster_key` | string | yes | — | Column in `adata.obs` containing cluster labels |
| `model_species` | list[string] | yes | — | One or more model species for BLAST. Must appear in `list_species.py` output |
| `non_model_fasta` | string | yes | — | Path to FASTA file of the non-model species |
| `non_model_h5ad` | string | yes | — | Path to h5ad file of the non-model species |
| `organ` | string | yes | — | Organ filter for the knowledge graph. Must appear in `list_organs.py` output |
| `outdir` | string | yes | — | Output directory (created if missing) |
| `marker_weight_method` | string | no | `prod` | How to combine marker scores: `prod` (product) or `sum` (sum) |
| `top_num` | int | no | `30` | Top N genes retained per cluster after step 4 |
| `homolo_weight_key` | string | no | `pident` | Which BLAST column to use as the homology weight |
| `multihomolo` | bool | no | `true` | If `true`, keep multiple homologs per non-model gene |
| `decay_factor` | float | no | `0.7` | Graph weight decay at each propagation hop |
| `mode` | string | no | `path` | `node` (single type) or `path` (trajectory-aware) |
| `threshold` | float \| null | no | `null` | Z-score threshold for filtering low-weight cell types |
| `candidate_annotation` | list[string] | no | `null` | Restrict the output cell-type set |

## Minimal example

If you do not need to override any defaults:

```yaml
cluster_key: leiden
model_species:
  - Brassica_rapa
non_model_fasta: Arabidopsis_thaliana.fasta
non_model_h5ad: ERP132245.h5ad
organ: leaf
outdir: output/ERP132245
```

## Naming convention for `outdir`

Recommended: include the dataset label in the path (e.g.,
`output/PRJNA662627/`) so multiple runs do not collide.
