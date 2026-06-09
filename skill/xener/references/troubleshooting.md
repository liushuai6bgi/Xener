# Troubleshooting

Common errors and how to resolve them.

## "ModuleNotFoundError: No module named 'xener'"

You have not installed the package yet. Run:

```bash
pip install xener
# or
python scripts/install.py
```

Then verify with `python -c "import xener"` (the only allowed inline import).

## "Invalid species" / "Invalid organ"

The `model_species` or `organ` value is not in the BLAST / KG database.
Run the validators first:

```bash
python scripts/list_species.py
python scripts/list_organs.py
```

The output values must be used **exactly** (case-sensitive, no
abbreviation). Cross-reference with `references/config-schema.md`.

## Step 3 fails with "BLAST database not found"

The BLAST database for the chosen model species is not installed
locally. Check the `xener` package documentation for BLAST database
installation. This is a package-level issue, not a skill issue.

## Step 5 returns an empty annotation

Likely causes (in order of likelihood):
1. The chosen `organ` does not match the cluster's biology
2. The BLAST homologs in step 3 are too stringent (raise `--pident`,
   lower `--evalue`)
3. `top_num` is too small — try increasing from 30 to 50
4. The cluster has too few differentially expressed genes

Re-run individual steps with relaxed parameters — see
`workflows/step-by-step.md`.

## Refinement fails with "celltype not found"

`--celltype` values must be exact matches in
`celltype_weight.csv` for the target cluster. Run:

```bash
python scripts/suggest_refine.py \
    --input output/celltype_weight.csv \
    --topk 5 \
    --outdir output/
```

Then use the exact cell type names from the output.

## Pipeline runs but UMAP looks wrong

The h5ad may be missing `X_pca`, `connectivities`, or `X_umap`.
The pipeline will compute them automatically. If clustering was
done with a non-default resolution, the resulting `leiden` clusters
may merge or split biologically meaningful groups.

## Out of memory

Try:
- Reducing `top_num`
- Running steps one at a time (each step's output is checkpointed)
- Subsampling the h5ad before running

## "Permission denied" writing to `outdir`

Ensure the directory exists and is writable. The skill does not
auto-create parent directories for `--outdir` in all scripts —
create it manually first if nested paths are used.

## `plot_umap.py` says "Embedding 'X_umap' not in adata.obsm"

The h5ad has not been UMAP-projected. Run `sc.tl.umap(adata)` first
(or use the `--embedding-key` flag if a different key is present,
e.g., `X_tsne`).

## `plot_umap.py` says "Column 'xener' not found"

The h5ad was saved before annotation. Re-run the pipeline (it
writes the annotation into `adata.obs['xener']` when using
`Xener.run_from_yaml`). For step-by-step usage, manually add the
column:

```python
adata.obs['xener'] = adata.obs[cluster_key].map(cluster2celltype)
```

## Refinement panel is missing in the refine plot

The `--refine-key` column (default `xener_refine`) is not in
`adata.obs`. Run `refine_cluster.py` first; it adds this column.

## Need to compare multiple model species?

Run `step3_mapping.py` once per species, then concatenate the
`gene_homolo_weight.csv` outputs (or pass multiple `--species`
arguments in a single step-3 invocation).
