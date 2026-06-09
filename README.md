# Xener
> This is the public version, containing only the necessary code.

A cross-species single-cell cell type annotation tool using knowledge graph.

## Installation

```bash
pip install .
# or
pip install xener
```

## Quick Start

### With a YAML config file

```python
from xener import Xener

annor = Xener()
cluster2celltype, _, debug_params = annor.run_from_yaml('config.yaml')
```

`config.yaml` example:

```yaml
cluster_key: leiden
model_species:
- Brassica_rapa
non_model_fasta: Arabidopsis_thaliana.fasta
non_model_h5ad: ERP132245.h5ad
organ: leaf
outdir: output/ERP132245
mapping_strict: 0       # <0=ignore BLAST quality, 0=balanced, 1=suppress multi-copy families
ann_strict: 0           # <0=exploratory (more types), 0=balanced, 1=cleaner (1 type/marker), 2=strictest (1 type/cluster)
```

### Programmatic API

```python
from xener import Xener

annor = Xener()
cluster2celltype, _, debug_params = annor(
    non_model_h5ad='ERP132245.h5ad',
    cluster_key='leiden',
    outdir='output/ERP132245',
    non_model_fasta='Arabidopsis_thaliana.fasta',
    model_species=['Brassica_rapa'],
    organ='leaf',
    mapping_strict=0,   # <0=ignore BLAST quality, 0=balanced, 1=suppress multi-copy families
    ann_strict=0,       # <0=exploratory (more types), 0=balanced, 1=cleaner (1 type/marker), 2=strictest (1 type/cluster)
)
```

Defaults for `marker_weight_method`, `mode`, `decay_factor`, `multihomolo`, `top_num`, etc. come from the default config and can be overridden as keyword arguments.

The third return value `debug_params` is a dict recording the actual parameters used in each key step, saved as `debug_params.yaml` in the output directory. It helps with reproducibility.

## Step-by-step

The `__call__` API above is the simplest way to run the full pipeline. If you need fine-grained control, you can call each step individually:

```python
from xener import Xener
import scanpy as sc

annor = Xener()
adata = sc.read('ERP132245.h5ad')
cluster_key = 'leiden'
non_model_fasta = 'Arabidopsis_thaliana.fasta'
model_species = ['Brassica_rapa']
organ = 'leaf'
outdir = 'output/ERP132245'

marker_gene = annor.get_markers(adata, cluster_key)

marker_weight, debug_gw = annor.get_gene_weight(marker_gene)

gene_homolo_weight, debug_map = annor.mapping(marker_weight, non_model_fasta, model_species, outdir)

topk_markers, debug_topk = annor.get_topk_gene(gene_homolo_weight, top_num=30)
# Only the top 30 genes will be retained for the subsequent steps.

cluster2celltype, _, celltype_weight, debug_ann = annor.cell_annotation(
    topk_markers, outdir / 'annotation', organ)

# Collect and save debug_params for reproducibility
import yaml
debug_params = {}
debug_params['get_gene_weight'] = debug_gw
debug_params['mapping'] = debug_map
debug_params['get_topk_gene'] = debug_topk
debug_params['cell_annotation'] = debug_ann
with open(outdir / 'debug_params.yaml', 'w') as f:
    yaml.dump(debug_params, f, default_flow_style=False)
```

Each step function (except `get_markers`) returns `(result, debug_params)` — a dict of the actual parameters used internally. Unlike `__call__`, step-by-step mode requires you to collect and save them explicitly.

## Output directory

```
outdir/
├── marker_gene.csv
├── marker_weight.csv
├── blastp_{species}.csv            # one per model species
├── gene_homolo_weight.csv
├── topk_markers.csv
├── celltype_weight.csv
├── debug_params.yaml               # actual parameters used in each step
├── config.yaml                     # from run_from_yaml only
└── annotation/
    ├── cluster_{id}_gene2celltype.xml  # The annotation path of this cluster
    └── ...
```

## debug_params.yaml

This file records the actual parameter values used in each key step of the pipeline, making results reproducible:

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

## Strict modes

The pipeline supports three levels of strictness (`mapping_strict`, `ann_strict`) for filtering ambiguous homology and annotation results:

| Value | mapping_strict (BLAST) | ann_strict (KG query) |
|-------|----------------------|----------------------|
| **< 0** | All `homolo_weight` = 1. Every homolog carries equal weight regardless of identity. **Result: more genes retained, but BLAST quality is ignored — use when BLAST scores are unreliable.** | All non-zero weights → 1. Every marker–celltype link gets equal vote. **Result: more cell types predicted per cluster, but with lower precision — use for broad exploratory annotation.** |
| **0** (default) | Raw BLAST weights as-is. **Result: balanced — strong homologs contribute more, weak ones less.** | Raw KG weights as-is. **Result: balanced default.** |
| **1** | Keep only the top weight per (gene group, gene). **Result: prevents multi-copy gene families from dominating within a species group.** | Keep only the max-confidence cell type per marker gene (row-wise). **Result: each marker votes for one cell type — fewer, cleaner predictions per cluster.** |
| **2** | Same as 0. | Keep only the single global max-confidence cell type. **Result: only one cell type survives — use when you need a single definitive answer per cluster.** |

In the `refine_single_cluster` step, `strict>0` keeps only the max-confidence cell type per gene (column-wise) and removes zero rows/columns. **Result: cleaner sub-cluster annotation with fewer ambiguous assignments.**

```python
# Example: strict annotation
cluster2celltype, _, _, _ = annor.cell_annotation(
    topk_markers, outdir / 'annotation', organ,
    ann_strict=1)  # keep only the best cell type per marker gene
```

## Sub-cluster refinement

```python
cluster_id = 0
candidate_celltype = ['type1', 'type2']
# Only support the values that appear in celltype_weight[celltype_weight['cluster'] == cluster_id]['celltype'].unique()
key_added = 'xener_refine'
moranI_threshold = 0.5
# moranI_threshold: Moran's I filtering threshold in [-1, 1]. The closer to 1, the stricter it is.
# Values outside [-1, 1] skip the screening step.
strict = 0                 # strict>0: keep only the max-confidence cell type per gene (column-wise)
split_method = 'bindiv'    # 'bindiv' = binary division, 'argmax' = argmax over cell-type scores
markergene_method = 'diff' # 'diff' = use only differential markers, 'all' = use all markers
# celltype_geneCount_gene: optional precomputed (celltype, gene_count, gene_list) list.
# If None, it is computed from the knowledge graph.

geneCount, diffgeneCount, annotation, gene2celltype_g = annor.refine_single_cluster(
    adata, topk_markers, cluster_key, cluster_id, candidate_celltype,
    key_added, organ,
    moranI_threshold=moranI_threshold, strict=strict,
    split_method=split_method, markergene_method=markergene_method)
# geneCount: list of (celltype, gene_count, gene_list) — markers per cell type.
# diffgeneCount: list of (celltype, diff_gene_count, diff_gene_list) — unique markers per cell type.
# annotation: DataFrame of adata_sub.obs[[key_added]]; the refined sub-cluster labels.
# gene2celltype_g: NetworkX graph linking gene -> homolo -> celltype. use `nx.write_gexf` to save.
# The results can be found in the returned annotation[key_added] DataFrame.
```

## Log interpretation and debugging

By default, xener logs at `INFO` level to the console with the format
`[YYYY-MM-DD HH:MM:SS] [xener] [LEVEL] message`. For debugging or post-mortem
analysis, redirect logs to a file and/or raise verbosity:

```python
from xener import Xener
from xener.utils.logger import setup_logger, add_file_handler

# File handler (unbuffered, appends by default, see xener/utils/logger.py).
logger = setup_logger(level=20)  # 20 = logging.INFO
add_file_handler(logger, log_file='xener.log', unbuffered=True)

# Switch to DEBUG to also see raw KG HTTP status / bytes, BLASTP command,
# build_graph_from_adjust_matrix internals, etc.
logger.setLevel(10)  # 10 = logging.DEBUG
```

### Pipeline landmarks (top-level grep targets)

The pipeline emits a small set of unambiguous markers you can `grep` to
reconstruct what happened without reading the full log:

| Grep target | Meaning |
|---|---|
| `>>>Xener pipeline started` | Entry of `__call__`. Tail is the resolved config (outdir, model_species, organ, top_num, mapping_strict, ann_strict). |
| `Resuming from step` / `All 4 checkpoints present` / `No usable checkpoint` | Tells you whether this run re-used cached zip files or started from scratch. |
| `>>>cell annotation` | Start of step 5; subsequent per-cluster logs are indented to its context. |
| `>>>Xener pipeline finished in NNs` | End of `__call__`. Tail lists top-5 celltype frequencies. |
| `KG get_gene2celltype_kg done in` / `KG get_celltype2celltype_kg done in` | Each KG round trip — slow network calls. |
| `BLASTP done in` / `BLASTP cache hit` | Each `blastp` invocation (real run vs. cache). |
| `>>>refine started` / `>>>refine finished` | Sub-cluster refinement span. |

The elapsed time printed in the `finished` line is the wall-clock duration
of the whole `__call__` (excluding Python import / setup).

### Stage-by-stage debug recipes

**1. `Xener.__init__` — environment problems**

| Log line | What it means | What to do |
|---|---|---|
| `KG backend: bolt (url=...)` / `KG backend: http (url=...)` | Confirms which backend the client picked based on the URL scheme. | If you intended the other backend, fix `KG_url` in your config. |
| `Bundled blastdb missing, downloading blastdb.zip...` | First-time setup is fetching the bundled database. | One-time; wait. If it hangs, check network / `xenor.dcs.cloud` reachability. |
| `Provided blastdb_path ... does not exist, falling back to bundled database.` | Your `blastdb_path` was wrong, but xener kept going. | Verify the path; the run will still succeed but use the bundled DB. |
| `KG species_organ_cell loaded: N unique organs` | KG round trip completed. | If 0, the KG is empty or unreachable — earlier lines will show a stack trace. |

**2. `get_markers` — input data sanity**

| Log line | What it means | What to do |
|---|---|---|
| `raw_available[True/False]` | Whether `adata.raw` is set. | Required for `use_raw=True` downstream. |
| `use_raw[True/False]` | Whether `rank_genes_groups` ran on `.raw`. | If `False` and you expected `True`, your h5ad may have been saved without raw. |
| `Unavailable data! cann't find available counts.` | Both `.X` and `.raw.X` are non-sparse. | Re-save the h5ad with raw counts (`sc.pp.log1p` after `normalize_total` won't work — you need raw counts in `.raw`). |
| `highly_variable_genes[N]` | Triggered by `adata.shape[1] > 4000` (or `force_HVG=True`). | Informational; reduces gene count to `N`. |

**3. `get_gene_weight` — numeric edge cases**

| Log line | What it means | What to do |
|---|---|---|
| `replace_inf logfc: X` / `replace_zero pvals_adj: Y` | Real inf/zero values were clipped. | The chosen X/Y shows the magnitude; if X is huge your data has extreme fold-changes — review upstream DE. |
| `logfc has inf values` / `pts_delta has inf values` / `-log10_pvals_adj has inf values` / `weight has inf values` | A residual inf slipped through after clipping. | Almost always a symptom of bad DE input (zero-variance clusters, all-zero pct1, etc.). |

**4. `mapping` — BLASTP**

| Log line | What it means | What to do |
|---|---|---|
| `BLASTP starting: query=..., db=..., threads=N` | New BLASTP run started. | Wait; the next `BLASTP done in` line will report wall time. |
| `BLASTP failed with returncode=N` | BLASTP exited non-zero. | xener now raises `RuntimeError`; the message includes stderr. |
| `BLASTP cache hit: ...csv (N rows, skipping alignment)` | Reused a previous BLASTP result. | If you expected a fresh run, delete `blastp_<species>.csv` or pass a different `outdir`. |
| `BLASTP best-hits per (qseqid,sseqid): N rows` | Number of unique query-subject pairs after `idxmax` on bitscore. | Sanity check vs. your marker count — if dramatically lower, your pident/evalue/bitscore thresholds are too strict. |
| `N groups merged to M groups in mapping` | Inner join on BLAST result dropped some marker rows → clusters have no BLAST hits. | Check your `non_model_fasta` actually contains the marker genes. |
| `mapping_strict[N] is too loose!` / `multiple mapping detected!` | `mapping_strict` collapsed multi-copy families. | See "Strict modes" section. |

**5. `cell_annotation` — per-cluster decisions**

Each `cluster_X` logs:

1. **KG query** — `KG get_gene2celltype_kg done in Ts: G genes, C celltypes, matrix=..., nnz=N`
2. **Decay** — `total X% homolos of organ[O] not in kg` (likely the most common reason for low-quality annotations)
3. **Branch** — one of:
   - `single candidate "X" after threshold, returning directly`
   - `top1 z-score>3 and top2<3, returning top1 "X" without aggregation`
   - `ambiguous top types (top1 z=A, top2 z=B), running celltype2celltype aggregation with N candidates`
   - `no candidate type after threshold, set to "unknown"`
4. **Path mode only** — `path step N, current="X", M parent candidates, kept K (softmax>0.15)`
5. **Save** — `gexf graph saved to .../cluster_X_gene2celltype.xml (nodes=N, edges=M)`

If many clusters end up `unknown`:

1. Check `total X% homolos of organ[O] not in kg` — if X is high (>30%), the KG doesn't have your species' homologs at all. Try a different `organ` or a different `model_species`.
2. Check the threshold — `celltype` numbers after threshold tell you how aggressive the filter is.
3. Inspect the saved `cluster_X_gene2celltype.xml` in Gephi to see which homolos connected to which celltypes.

**6. `refine` / `refine_single_cluster` — sub-cluster debugging**

| Log line | What it means | What to do |
|---|---|---|
| `>>>refine started: cluster_key=..., topk=..., organ=..., strict=..., key_added=...` | Effective refine parameters. | Cross-check against the values you passed in. |
| `refine summary: N clusters queued for refinement, M skipped (single candidate).` | How many clusters actually entered `refine_single_cluster`. | If 0, your `celltypes_weight` is dominated by singletons or the top-`topk` types have no clear weight drop. |
| `refine: skip cluster X (only 1 candidate celltype, need >=2 to refine).` | Cluster has only one candidate; nothing to disambiguate. | Informational. |
| `refine: skip cluster X (no clear weight drop in top-5 celltypes, dominant type is unambiguous).` | Top celltype has no competitor within `topk`. | Increase `topk` or accept the dominant type. |
| `>>>refine_single_cluster: cluster_id[X], cluster_key[Y], candidate_celltype[Z]` | Entry to single-cluster refinement. | If cluster_id is wrong, fix your `celltypes_weight` upstream. |
| `refine_single_cluster effective params: moranI_threshold=..., strict=..., split_method=..., markergene_method=...` | Resolved values of all four refinement switches. | Double-check these — they drive the branch taken below. |
| `strict=N: keep_max per column on KG matrix ...` | `strict>0` column-wise filtering active. | See "Strict modes" section. |
| `markergene_method=...: N celltypes in queue, M unique markers` | Per-cluster queue size. | Empty queue ⇒ no annotation will happen. |
| `note: bindiv mode will add temporary columns to adata.obs: f"{idx}_{type}_EXP", ...` | Heads-up that `adata.obs` is mutated. | If you don't want these columns, use `split_method='argmax'`. |
| `argmax result distribution: [...]` | argmax branch — count of cells assigned to each celltype. | Concentrated in one type → confident. Spread → ambiguous. |
| `refine_single_cluster X total: [...]` | Final unique labels in this cluster after refinement. | If all `waitting` or `unknown`, no cell passed the score threshold. |
| `gene2celltype_g built: N nodes, M edges` | The returned NetworkX graph. | M=0 is suspicious — refinement produced no structure. |
| `ValueError: refine_single_cluster: cluster_id[X] not in group_gene_homolo_weight groups. Available groups: [...]` | The cluster you're trying to refine has no `gene_homolo_weight` rows. | Pass the correct `cluster_id` (must exist in `group_gene_homolo_weight['group']`). |

### Common patterns and quick fixes

| Symptom | Grep this | Likely cause | Fix |
|---|---|---|---|
| Pipeline hangs after `BLASTP starting` | `BLASTP done` not seen | BLASTP is genuinely slow; or crashed silently. | Wait several minutes; if you see `BLASTP failed` you have a real error. Otherwise `kill -9` and reduce `top_num` for a smoke test. |
| Many `KG get_gene2celltype_kg returned an empty matrix (organ=..., N homolos)` | as written | KG has no homolo→celltype edges for your `organ` filter. | Drop `organ` (set to `None`) or use a different organ; verify the gene names match species conventions. |
| All `celltype=unknown` | `set to "unknown"` | Threshold too strict, or no candidates after graph propagation. | Inspect the gexf file; lower `threshold`; try `mode='node'` (less aggressive than `mode='path'`). |
| `multiple mapping detected!` warnings | `multiple mapping` | Your `mapping_strict=1` is keeping ties across multi-copy families. | Set `mapping_strict=0` or check for duplicated gene entries upstream. |
| `top1 z-score>3 and top2<3, returning top1` for every cluster | `returning top1` | KG graph propagation is not differentiating between clusters. | Verify the species-homolog overlap is non-trivial; check `gene_homolo_weight.shape` is not tiny. |
| Refinement modifies `adata.obs` with `_EXP` columns | `will add temporary columns to adata.obs` | Expected behavior of `split_method='bindiv'`. | Use `split_method='argmax'` if you need a clean adata. |
| `checkpoint invalid (expected N gene-group combos, got M)` | `Checkpoint invalid` | Your `top_num` changed since the last run, so the cached `topk_markers.csv` no longer matches. | Delete `topk_markers.csv` (and downstream) before re-running with new `top_num`. |
| `KG HTTP <METHOD> <PATH> returned <CODE>` | `KG HTTP ... returned` | KG server is unhealthy or URL is wrong. | Retry; verify `KG_url`; check upstream KG health. |
| Pipeline takes much longer than usual | `BLASTP done in` / `KG get_gene2celltype_kg done in` | One stage has slowed down. | Check the per-stage elapsed times — KG vs. BLASTP vs. `cell_annotation` pinpoint the bottleneck. |

## Claude Code / AI Agent Skill

An agent skill definition is provided at [skill/xener/SKILL.md](skill/xener/SKILL.md) for use with Claude Code and compatible AI coding assistants. The skill guides the agent through the full xener workflow: installation → config validation → species/organ selection → pipeline execution → refinement suggestion.

### Scripts

All available scripts are in [skill/xener/scripts/](skill/xener/scripts/):

| Script | Purpose |
|--------|---------|
| `install.py` | Install xener package |
| `list_species.py` | Query available reference species from BLAST database |
| `list_organs.py` | Query available organs from knowledge graph |
| `run_pipeline.py` | Run the full annotation pipeline from a config file |
| `step1_markers.py` | Get marker genes per cluster |
| `step2_weight.py` | Calculate gene weights |
| `step3_mapping.py` | BLAST homology mapping |
| `step4_topk.py` | Get top-k genes |
| `step5_annotate.py` | Cell type annotation |
| `refine_cluster.py` | Sub-cluster refinement |
| `suggest_refine.py` | Suggest refinement candidates from annotation results |

For agent usage details, refer to [skill/xener/SKILL.md](skill/xener/SKILL.md).

## Links

[Homepage](https://xenor.dcs.cloud/): https://xenor.dcs.cloud/

[PyPI](https://pypi.org/project/xener/): https://pypi.org/project/xener

[Github](https://github.com/liushuai6bgi/Xener): https://github.com/liushuai6bgi/Xener
