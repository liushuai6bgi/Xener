# Parameters Reference

Every CLI flag exposed by the scripts. Use this file when tuning
individual steps.

## Step 1: Get marker genes (`step1_markers.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to h5ad file |
| `--cluster-key` | required | Column in `adata.obs` with cluster labels |
| `--outdir` | required | Output directory |

## Step 2: Calculate gene weights (`step2_weight.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | `marker_gene.csv` from step 1 |
| `--method` | `prod` | `prod` (product) or `sum` (sum of log scores) |
| `--outdir` | required | Output directory |

## Step 3: BLAST homology mapping (`step3_mapping.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | `marker_weight.csv` from step 2 |
| `--fasta` | required | Path to non-model species FASTA |
| `--species` | required | Model species name (repeat for multiple) |
| `--weight-key` | `pident` | Column to use: `pident`, `evalue`, or `bitscore` |
| `--pident` | `60` | Minimum percent identity filter |
| `--evalue` | `0.05` | Maximum e-value filter |
| `--bitscore` | `200` | Minimum bitscore filter |
| `--outdir` | required | Output directory |

## Step 4: Top-k genes (`step4_topk.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | `gene_homolo_weight.csv` from step 3 |
| `--k` | `30` | Top N genes per cluster |
| `--multihomolo` | `true` | Keep multiple homologs per gene |
| `--outdir` | required | Output directory |

## Step 5: Cell type annotation (`step5_annotate.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | `topk_markers.csv` from step 4 |
| `--outdir` | required | Output directory (recommended: `output/annotation`) |
| `--organ` | `None` | Organ filter for the knowledge graph |
| `--mode` | `path` | `node` (single type) or `path` (trajectory) |
| `--decay-factor` | `0.7` | Graph weight decay |
| `--threshold` | `None` | Z-score threshold for cell-type filtering |
| `--candidate-annotation` | `None` | Restrict to specific cell types (space-separated) |

## Refinement (`refine_cluster.py`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Path to h5ad file |
| `--markers` | required | `gene_homolo_weight.csv` from step 3 |
| `--cluster-key` | required | Column in `adata.obs` |
| `--cluster-id` | required | ID of the cluster to refine |
| `--celltype` | required | Comma-separated candidate cell types (must come from `celltype_weight.csv`) |
| `--organ` | required | Organ filter |
| `--moran-i` | `0.5` | Moran's I threshold (range [-1, 1]; >0.5 = stricter) |
| `--split-method` | `argmax` | `argmax` (assign by max score) or `bindiv` (binary division) |
| `--markergene-method` | `all` | `all` (all top-k genes) or `diff` (differential only) |
| `--key-added` | `xener_refine` | Column name in `adata.obs` for refined annotation |
| `--outdir` | required | Output directory |

> **Defaults reflect the empirically best combination:**
> `--markergene-method all --split-method argmax` has the **highest
> refinement success rate** (it actually splits mixed clusters into
> sub-populations rather than collapsing them to a single label). The old
> defaults (`diff` / `bindiv`) frequently failed to split: `diff` returns 0
> markers whenever the two candidate cell types share homologs in the KG, and
> `bindiv` then assigns every cell to one candidate (a `[N, 0]` result).
> Override to `diff` / `bindiv` only when you specifically want the stricter
> differential/binary behavior.

## Full pipeline (`run_pipeline.py`)

Accepts a YAML config and runs all five steps in order. The CLI flags
above map to YAML fields; see `config-schema.md` for the mapping.
