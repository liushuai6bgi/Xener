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

```python
from xener import Xener

# Initialize
annor = Xener()

# Run full pipeline
cluster2celltype, _ = annor.run_from_yaml('config.yaml')
```

`config.yaml` example.

```yaml
cluster_key: leiden
model_species:
- Brassica_rapa
non_model_fasta: Arabidopsis_thaliana.fasta
non_model_h5ad: ERP132245.h5ad
organ: leaf
outdir: output/ERP132245
```

## Step-by-step

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

marker_weight = annor.get_gene_weight(marker_gene)

gene_homolo_weight = annor.mapping(marker_weight, non_model_fasta, model_species, outdir)

topk_markers = annor.get_topk_gene(gene_homolo_weight, k=30)
# Only the top 30 genes will be retained for the subsequent steps.

cluster2celltype, _, _, _, _ = self.cell_annotation(
            topk_markers, annotation_info_path, organ)
```

## Sub-cluster refinement

```python
cluster_id = 0
candidate_celltype = ['type1', 'type2']
key_added = 'xener_refine'
moranI_threshold = 0.5
# moranI_threshold used for gene screening, the effective value ranges from [-1, 1]. The closer to 1, the stricter it is. If an invalid value is input, the screening step will be skipped.

annor.refine_single_cluster(adata, topk_markers, 
            cluster_key, cluster_id, candidate_celltype, 
            key_added, organ, moranI_threshold)
# The results can be found in  adata.obs[key_added]
```

## Links

Homepage: https://xenor.dcs.cloud/
