# Step-by-Step Workflow

Run individual steps to fine-tune parameters. Each step reads its
input `.csv` from the previous step and writes its own `.csv`
output. You can re-run a step with new parameters without
re-running earlier steps.

## Step 1: Get marker genes

```bash
python scripts/step1_markers.py --input data.h5ad --cluster-key leiden --outdir output/
```

→ writes `marker_gene.csv`

## Step 2: Calculate gene weights

```bash
python scripts/step2_weight.py --input output/marker_gene.csv --method prod --outdir output/
```

→ writes `marker_weight.csv`

**Tuning tip**: switch `--method` between `prod` and `sum` to see
which gives more balanced weights across clusters.

## Step 3: BLAST homology mapping

```bash
python scripts/step3_mapping.py --input output/marker_weight.csv \
    --fasta Arabidopsis_thaliana.fasta \
    --species Brassica_rapa \
    --weight-key pident \
    --pident 60 --evalue 0.05 --bitscore 200 \
    --outdir output/
```

→ writes `blastp_{species}.csv` and `gene_homolo_weight.csv`

**Tuning tips**:
- Lower `--pident` (e.g., 40) for distantly related species
- Raise `--bitscore` (e.g., 300) for higher confidence homologs
- Switch `--weight-key` to `bitscore` if `pident` produces
  uniform weights

## Step 4: Top-k genes

```bash
python scripts/step4_topk.py --input output/gene_homolo_weight.csv --k 30 \
    --multihomolo --outdir output/
```

→ writes `topk_markers.csv`

**Tuning tip**: try `--k 50` if step 5 returns empty annotations.

## Step 5: Cell type annotation

```bash
python scripts/step5_annotate.py --input output/topk_markers.csv \
    --outdir output/annotation \
    --organ leaf \
    --mode path --decay-factor 0.7
```

→ writes `celltype_weight.csv` and per-cluster XMLs in `annotation/`

**Tuning tips**:
- Switch `--mode` to `node` if `path` produces trajectory-like
  cell types you don't want
- Lower `--decay-factor` (e.g., 0.5) to weight local cell types
  more strongly
- Add `--candidate-annotation type1 type2` to restrict the output
  set (useful for debugging)

## Re-running a step

To re-run step 3 with different BLAST thresholds, simply invoke
`step3_mapping.py` again with the same `marker_weight.csv` and
new flags. The output `gene_homolo_weight.csv` is overwritten.
