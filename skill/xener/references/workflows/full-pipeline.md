# Full Pipeline Workflow

Run the entire Xener pipeline end-to-end from a YAML config.

## Prerequisites

- `xener` is installed (`pip install xener` or `python scripts/install.py`)
- A valid `config.yaml` exists (see `config-schema.md`)
- `model_species` and `organ` have been confirmed with the user
  (see `config-validation.md`)

## Command

```bash
python scripts/run_pipeline.py --config config.yaml
```

## What happens

1. **Step 1** — Marker genes are computed per cluster
2. **Step 2** — Gene weights are calculated (`prod` by default)
3. **Step 3** — BLAST maps each gene to the model species
4. **Step 4** — Top-k genes (default 30) are selected per cluster
5. **Step 5** — Knowledge-graph propagation predicts cell types

All intermediates are saved as `.csv` files in `outdir/`. The final
`debug_params.yaml` records every parameter actually used.

## Outputs

See `output-files.md` for the full list. Key files:
- `celltype_weight.csv` — all predicted cell types per cluster
- `annotation/cluster_{id}_gene2celltype.xml` — per-cluster graph
- `debug_params.yaml` — for reproducibility

## Next step

After the full pipeline, consider running refinement on clusters
with mixed cell types. See `refinement.md`.
