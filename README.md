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

## Sub-cluster refinement

```python
cluster_id = 0
candidate_celltype = ['type1', 'type2']
# Only support the values that appear in celltype_weight[celltype_weight['cluster'] == cluster_id]['celltype'].unique()
key_added = 'xener_refine'
moranI_threshold = 0.5
# moranI_threshold used for gene screening, the effective value ranges from [-1, 1].
# The closer to 1, the stricter it is. If an invalid value is input, the screening step will be skipped.

geneCount, diffgeneCount, annotation = annor.refine_single_cluster(
    adata, topk_markers, cluster_key, cluster_id, candidate_celltype,
    key_added, organ, moranI_threshold)
# The results can be found in the returned annotation[key_added] DataFrame.
```

## Links

[Homepage](https://xenor.dcs.cloud/): https://xenor.dcs.cloud/

[PyPI](https://pypi.org/project/xener/): https://pypi.org/project/xener

[Github](https://github.com/liushuai6bgi/Xener): https://github.com/liushuai6bgi/Xener
