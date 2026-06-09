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

## Self-debug: when the pipeline "completes" but the result is wrong

The pipeline can exit 0 and still produce biologically useless
output. In that case, **do not re-run blindly** — diagnose first.

1. **Read `references/log-interpretation.md`** for the full
   log-driven debug protocol. The xener logger emits a small fixed
   set of strings; you can `grep` them.
2. **Run the post-run quality gate** explicitly:
   `python scripts/check_output.py --outdir <outdir>`. It catches
   the five most common silent failures (mean KG miss, tail of
   severe miss, too-few unique top-1 cell types, weak-confidence
   clusters, empty annotations).
3. **Grep pipeline landmarks**:
   ```bash
   grep -E '>>>|ERROR|WARNING|Traceback' outdir/xener.log
   ```
4. **Quick decision table** — match the symptom to the cause:

| Symptom in log / output | Likely cause | Fix |
|---|---|---|
| `total X% homolos of organ[...] not in kg` with mean `X > 0.30` | KG coverage gap | Add target species (if model organism) or a closer relative — see `workflows/species-selection.md` |
| `no candidate type after threshold, set to "unknown"` repeated | `threshold` too strict | Re-run step 5 with `--threshold null` |
| `multiple mapping detected!` | `mapping_strict=1` collapsing ties | Set `mapping_strict=0` |
| `>>>Xener pipeline finished` but all clusters have weight ≈ 0 | Empty BLAST result (wrong fasta / wrong model_species) | Verify `non_model_fasta` contains marker genes; verify `model_species` are in `list_species.py` |
| `BLASTP failed with returncode=N` | BLASTP exit error | Re-read stderr in error message; check disk space and FASTA format |
| `>>>cell annotation` line shows `0 celltypes before threshold` | Empty `topk_markers.csv` | Re-run step 4 with higher `--top-num` (try 50) |
| `refine summary: 0 clusters queued for refinement` | Top celltypes have no clear weight drop | Increase `topk` in step 4 or accept dominant types |

For the full stage-by-stage table (Stages 1–6) and the common
patterns cheat sheet, see `references/log-interpretation.md`.
