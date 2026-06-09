# Log Interpretation and Debugging

**Read this whenever a xener run produces output that looks wrong, or when
the post-run gate (`scripts/check_output.py`) fails.** The xener logger
emits a small, fixed set of strings. You can `grep` them to reconstruct
what happened without reading the entire log. This reference tells you
which string indicates which failure mode, and what to change.

## How to capture the log

`run_pipeline.py` mirrors xener's stdout to `outdir/xener.log`. For
fine-grained control, the user can call `setup_logger` directly:

```python
from xener import Xener
from xener.utils.logger import setup_logger, add_file_handler

logger = setup_logger(level=20)            # 20 = logging.INFO
add_file_handler(logger, log_file='xener.log', unbuffered=True)
logger.setLevel(10)                        # 10 = logging.DEBUG
```

Log format is `[YYYY-MM-DD HH:MM:SS] [xener] [LEVEL] message`.
After a run, the first thing to do is:

```bash
# Show pipeline landmarks and any error/warning
grep -E '>>>|ERROR|WARNING|Traceback' outdir/xener.log | head -100
```

## Top-level pipeline landmarks (always grep first)

These strings mark the boundaries of the pipeline phases. Their
**presence/absence** is the fastest health check.

| Grep target | Meaning | What it tells you |
|---|---|---|
| `>>>Xener pipeline started` | Entry of `__call__` | Tail lists resolved config (outdir, model_species, organ, top_num, mapping_strict, ann_strict) |
| `All 4 checkpoints present, skipping marker/weight/mapping/topk steps` | Cache hit, Step 5 only | Reused earlier `.csv` files. Safe but not a fresh end-to-end run. |
| `Resuming from step N` | Partial cache | Steps before N came from `.csv`, N onwards is fresh. |
| `No usable checkpoint, running the full pipeline from scratch` | Full fresh run | Every step re-ran; safe to inspect all intermediates. |
| `>>>cell annotation` | Start of Step 5 | Per-cluster logs are indented to its context. |
| `>>>Xener pipeline finished in NNs` | End of `__call__` | Wall-clock duration; tail lists top-5 celltype frequencies. |
| `KG get_gene2celltype_kg done in` / `KG get_celltype2celltype_kg done in` | One KG round trip | Slow network call. Count occurrences and sum elapsed times. |
| `BLASTP done in` / `BLASTP cache hit` | One `blastp` invocation | Real run vs. cache. |
| `>>>refine started` / `>>>refine finished` | Sub-cluster refinement span | Both present ⇒ refinement ran. |
| `>>>refine_single_cluster` | Entry to per-cluster refinement | Cross-check `cluster_id`, `candidate_celltype` against your request. |

The wall-clock in the `finished` line is the whole `__call__` excluding
Python import and setup.

## Stage 1: `Xener.__init__` — environment problems

| Log line | What it means | What to do |
|---|---|---|
| `KG backend: bolt (url=...)` / `KG backend: http (url=...)` | Confirms which backend was picked from the URL scheme | If wrong, fix `KG_url` in config |
| `Bundled blastdb missing, downloading blastdb.zip...` | First-time setup, fetching bundled DB | One-time. If it hangs, check `xenor.dcs.cloud` reachability. |
| `Provided blastdb_path ... does not exist, falling back to bundled database.` | Your `blastdb_path` was wrong, but xener kept going | Verify path; run still succeeds with bundled DB. |
| `KG species_organ_cell loaded: N unique organs` | KG round trip completed | If `N=0`, the KG is empty or unreachable — look above for a stack trace. |

## Stage 2: `get_markers` — input data sanity

| Log line | What it means | What to do |
|---|---|---|
| `raw_available[True/False]` | Whether `adata.raw` is set | Required for `use_raw=True` downstream. |
| `use_raw[True/False]` | Whether `rank_genes_groups` ran on `.raw` | If `False` and you expected `True`, the h5ad was saved without raw. |
| `Unavailable data! cann't find available counts.` | Both `.X` and `.raw.X` are non-sparse | Re-save the h5ad with raw counts in `.raw` (e.g. `adata.raw = adata` before writing). |
| `highly_variable_genes[N]` | Triggered by `adata.shape[1] > 4000` (or `force_HVG=True`) | Informational; gene count reduced to `N`. |
| `skip highly_variable_genes.` | Skipped because gene count is small | Informational. |

## Stage 3: `get_gene_weight` — numeric edge cases

| Log line | What it means | What to do |
|---|---|---|
| `replace_inf logfc: X` | `X` was the max-abs `logfc` × 10; real `+inf` was clipped to it | If `X` is huge, your data has extreme fold-changes — review upstream DE. |
| `replace_zero pvals_adj: Y` | Real zero `pv` was replaced with the smallest non-zero value × 0.5 | Sanity check `Y` is not a weird number. |
| `logfc has inf values` | Residual inf slipped through after clipping | Almost always bad DE input (zero-variance clusters, all-zero pct1, etc.). |
| `pts_delta has inf values` | Same — inf in pct delta | Same as above. |
| `-log10_pvals_adj has inf values` | Same — inf in `-log10(pv)` | Same as above. |
| `weight has inf values` | Same — final weight is inf | Same as above. |

If you see any of the residual-inf warnings, **stop and fix the input**
(marker table) — downstream steps will be unreliable.

## Stage 4: `mapping` — BLASTP

| Log line | What it means | What to do |
|---|---|---|
| `BLASTP starting: query=..., db=..., threads=N` | New BLASTP run started | Wait for the next `BLASTP done in` line. |
| `BLASTP failed with returncode=N` | BLASTP exited non-zero | xener raises `RuntimeError` with stderr in the message. Re-run after fixing the cause. |
| `BLASTP cache hit: .../blastp_<species>.csv (N rows, skipping alignment)` | Reused a previous result | If you expected a fresh run, delete `blastp_<species>.csv` or pass a different `outdir`. |
| `BLASTP best-hits per (qseqid,sseqid): N rows` | Unique query-subject pairs after `idxmax` on bitscore | If dramatically lower than your marker count, `pident`/`evalue`/`bitscore` are too strict. |
| `N groups merged to M groups in mapping` | Inner join on BLAST result dropped some marker rows | Check that `non_model_fasta` actually contains the marker genes. |
| `mapping_strict[-N] is too loose!` | `mapping_strict < 0` is active | See "Strict modes" in main README. |
| `multiple mapping detected!` | `mapping_strict=1` collapsed multi-copy families within a group | See "Strict modes" in main README. |

## Stage 5: `cell_annotation` — per-cluster decisions

Each cluster logs in this fixed order:

1. **KG query** — `KG get_gene2celltype_kg done in Ts: G genes, C celltypes, matrix=..., nnz=N`
2. **Decay** — `total X% homolos of organ[O] not in kg` (**the most common cause of bad annotations**)
3. **Branch** — one of:
   - `single candidate "X" after threshold, returning directly`
   - `top1 z-score>3 and top2<3, returning top1 "X" without aggregation`
   - `ambiguous top types (top1 z=A, top2 z=B), running celltype2celltype aggregation with N candidates`
   - `no candidate type after threshold, set to "unknown"`
4. **Path mode only** — `path step N, current="X", M parent candidates, kept K (softmax>0.15)`
5. **Save** — `gexf graph saved to .../cluster_X_gene2celltype.xml (nodes=N, edges=M)`

If many clusters end up `unknown`:

1. **Check the decay line first** — `total X% homolos of organ[O] not in kg`. If `X > 30%`, the KG has no edges for your `model_species` + `organ`. Try a different `organ` (or `organ=None`) or expand `model_species`.
2. Check the threshold — `celltype numbers after threshold` in the line `X: Y celltypes before threshold, Z after threshold=...` tells you how aggressive the filter is. If `Z=0`, threshold is too strict.
3. Inspect the saved `cluster_X_gene2celltype.xml` in Gephi / Cytoscape / `networkx` to see which homolos connected to which cell types.

## Stage 6: `refine` / `refine_single_cluster` — sub-cluster debugging

| Log line | What it means | What to do |
|---|---|---|
| `>>>refine started: cluster_key=..., topk=..., organ=..., strict=..., key_added=...` | Effective refine parameters | Cross-check against the values you passed. |
| `refine summary: N clusters queued for refinement, M skipped (single candidate).` | How many entered `refine_single_cluster` | If `N=0`, the top-`topk` types have no clear weight drop, or all are singletons. |
| `refine: skip cluster X (only 1 candidate celltype, need >=2 to refine).` | Single candidate, nothing to disambiguate | Informational. |
| `refine: skip cluster X (no clear weight drop in top-5 celltypes, dominant type is unambiguous).` | Top celltype has no competitor | Increase `topk` or accept the dominant type. |
| `>>>refine_single_cluster: cluster_id[X], cluster_key[Y], candidate_celltype[Z]` | Entry to single-cluster refinement | If `cluster_id` is wrong, fix your `celltypes_weight` upstream. |
| `refine_single_cluster effective params: moranI_threshold=..., strict=..., split_method=..., markergene_method=...` | Resolved values of all four switches | Double-check — these drive the branch taken below. |
| `strict=N: keep_max per column on KG matrix (shape ...)` | `strict>0` column-wise filtering active | See "Strict modes" in main README. |
| `markergene_method=...: N celltypes in queue, M unique markers` | Per-cluster queue size | Empty queue ⇒ no annotation. |
| `note: bindiv mode will add temporary columns to adata.obs: ...` | `adata.obs` will be mutated | Use `split_method='argmax'` if you need a clean adata. |
| `argmax result distribution: [...]` | argmax branch — cell counts per type | Concentrated ⇒ confident; spread ⇒ ambiguous. |
| `refine_single_cluster X total: [...]` | Final unique labels in this cluster | If all `waitting` or `unknown`, no cell passed the score threshold. |
| `gene2celltype_g built: N nodes, M edges` | Returned NetworkX graph | `M=0` is suspicious — refinement produced no structure. |
| `ValueError: refine_single_cluster: cluster_id[X] not in group_gene_homolo_weight groups. Available groups: [...]` | Cluster has no `gene_homolo_weight` rows | Pass a `cluster_id` that exists in `group_gene_homolo_weight['group']`. |

## Common failure patterns and quick fixes

This is the agent's first-stop decision tree. `grep` the run log for the
pattern in column 2; if you see it, apply the fix in column 4.

| Symptom | Grep this | Likely cause | Fix |
|---|---|---|---|
| Pipeline hangs after `BLASTP starting` | `BLASTP done` not seen | BLASTP is genuinely slow OR crashed silently | Wait several minutes; if `BLASTP failed` appears, you have a real error. Otherwise `kill -9` and reduce `top_num` for a smoke test. |
| Many `KG get_gene2celltype_kg returned an empty matrix` | `returned an empty matrix` | KG has no homolo→celltype edges for your `organ` filter | Drop `organ` (set to `None`) or use a different organ; verify gene names match species conventions. |
| All `celltype=unknown` | `set to "unknown"` | Threshold too strict, or no candidates after graph propagation | Inspect the gexf file; lower `threshold`; try `mode='node'` (less aggressive than `mode='path'`). |
| `multiple mapping detected!` warnings | `multiple mapping` | `mapping_strict=1` keeping ties across multi-copy families | Set `mapping_strict=0` or check for duplicated gene entries upstream. |
| `top1 z-score>3 and top2<3, returning top1` for every cluster | `returning top1` | KG graph propagation is not differentiating between clusters | Verify species-homolog overlap is non-trivial; check `gene_homolo_weight.shape` is not tiny. |
| Refinement modifies `adata.obs` with `_EXP` columns | `will add temporary columns to adata.obs` | Expected behavior of `split_method='bindiv'` | Use `split_method='argmax'` if you need a clean adata. |
| `Checkpoint invalid (expected N gene-group combos, got M)` | `Checkpoint invalid` | `top_num` changed since last run, so cached `topk_markers.csv` no longer matches | Delete `topk_markers.csv` (and downstream) before re-running with new `top_num`. |
| `KG HTTP <METHOD> <PATH> returned <CODE>` | `KG HTTP ... returned` | KG server is unhealthy or URL is wrong | Retry; verify `KG_url`; check upstream KG health. |
| Pipeline takes much longer than usual | `BLASTP done in` / `KG get_gene2celltype_kg done in` | One stage has slowed down | Compare per-stage elapsed times — KG vs. BLASTP vs. `cell_annotation` — to pinpoint the bottleneck. |
| Mean `total ... not in kg` > 30% | `not in kg` | KG coverage gap for this organ | Add the target species itself if it's a model organism, or pick a closer relative (see `workflows/species-selection.md`). |

## Debug workflow (apply after a failed run)

1. **Check the gate** — did `scripts/check_output.py` fail? If yes, jump
   straight to the table in `workflows/self-tuning-protocol.md`. The
   gate identifies the *class* of failure; this log-interpretation
   reference identifies the *exact* failure inside the class.
2. **Grep landmarks** — `grep -E '>>>|ERROR|WARNING|Traceback' outdir/xener.log`.
3. **Stage-by-stage** — find the first `WARNING`/`ERROR` and consult the
   matching section above.
4. **Verify input** — if the warning is in Stages 2–3, the issue is
   data quality. Re-save the h5ad (with raw) and re-run from Step 1.
5. **Adjust and re-run from the failing step** — do NOT re-run the
   whole pipeline. Use the `.csv` checkpointing to skip earlier steps.
   See `workflows/step-by-step.md` for per-step commands.

## Cross-references

- `workflows/self-tuning-protocol.md` — the **post-run quality gate**
  (run this first when a run "completes but is wrong")
- `references/troubleshooting.md` — installation / setup errors (not
  log-related; read this if xener won't even start)
- `references/output-files.md` — what each `.csv` in `outdir/` contains
- `references/autonomous-decision-making.md` — when to re-run with new
  parameters vs. accept the result
