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
├── marker_gene.zip
├── marker_weight.zip
├── blastp_{species}.zip            # one per model species
├── gene_homolo_weight.zip
├── topk_markers.zip
├── celltype_weight.zip
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
