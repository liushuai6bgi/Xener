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

> **Waiting for the run (mandatory-rules.md §9):** this step takes minutes.
> Either run it foreground with a generous Bash `timeout` (up to the 10-min
> max), or background it and **end your turn** — the completion notification
> re-invokes you. **Never** poll with `sleep N && grep`; to peek at progress
> use a bare `grep outdir/xener.log` (no `sleep`), to wait synchronously use
> `TaskOutput block:true`.

## What happens

1. **Step 1** — Marker genes are computed per cluster
2. **Step 2** — Gene weights are calculated (`prod` by default)
3. **Step 3** — BLAST maps each gene to the model species
4. **Step 4** — Top-k genes (default 30) are selected per cluster
5. **Step 5** — Knowledge-graph propagation predicts cell types
6. **Step 5.5** — The mandatory **quality gate runs in-process** (no
   subprocess, no h5ad re-read). It is **baseline-relative**: it hard-fails
   (pipeline exits non-zero) only on *structural breakage* or a *regression vs
   `--baseline`*. Being above an advisory target (e.g. mean KG miss > 30%) is a
   `[WARN]`, not a failure. To improve a WARN, add a KG-rich relative to
   `model_species` and re-run passing the previous outdir as `--baseline`
   (`run_pipeline.py --baseline <prev_outdir>`); keep the change if signals
   improved or held. Never change the confirmed `organ` to move a signal. See
   `self-tuning-protocol.md`.

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
