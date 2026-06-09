# Species Selection Heuristics

The agent should select `model_species` automatically based on the
target (non-model) species. This document provides the decision
heuristic.

## Quick algorithm

```text
1. Detect target species from:
   - h5ad var_names prefix (e.g., "AT" → Arabidopsis, "Zm" → maize)
   - h5ad uns['organism'] or similar metadata
   - filename hints
   - user prompt keywords

2. If detected, pick 1-3 species from the SAME taxonomic family:
   - Phylogenetic distance: same genus > same family > same order
   - Trade-off: more species = more homolog evidence, but more noise

3. If not detected, list available species to the user and ask
   for the target species first.
```

## Phylogenetic distance table (plants)

Use this table to map target species → recommended model species.

| Target family / species | Recommended model_species |
|-------------------------|---------------------------|
| Brassicaceae (Brassica, Arabidopsis) | Arabidopsis_thaliana, Brassica_rapa |
| Poaceae / grasses (rice, wheat, maize, bamboo) | Oryza_sativa, Zea_mays, Triticum_aestivum, Sorghum_bicolor |
| Fabaceae / legumes (soybean, Medicago, pea) | Glycine_max, Medicago_truncatula, Lotus_japonicus |
| Solanaceae (tomato, tobacco, potato) | Solanum_lycopersicum, Nicotiana_tabacum |
| Asteraceae (sunflower, lettuce) | Helianthus_annuus, Lactuca_sativa |
| Rosaceae (apple, strawberry, peach) | Fragaria_vesca, Prunus_persica, Malus_domestica |
| Cucurbitaceae (cucumber, melon) | Cucumis_sativus, Cucumis_melo |
| Malvaceae (cotton) | Gossypium_hirsutum |
| Vitaceae (grape) | Vitis_vinifera |
| Salicaceae (poplar) | Populus_trichocarpa |
| Others / unknown | Ask the user, or use Arabidopsis_thaliana + Oryza_sativa as a safe cross-kingdom pair |

## Phylogenetic distance table (animals)

| Target | Recommended model_species |
|--------|---------------------------|
| Primates | Homo_sapiens, Macaca_mulatta |
| Rodents | Mus_musculus, Rattus_norvegicus |
| Other mammals | Homo_sapiens, Mus_musculus |
| Zebrafish / fish | Danio_rerio |
| Chicken | Gallus_gallus |
| Drosophila | Drosophila_melanogaster |
| C. elegans | Caenorhabditis_elegans |

## BLAST threshold guidance by distance

Adjust `--pident` and `--evalue` based on phylogenetic distance:

| Distance | --pident | --evalue | --bitscore |
|----------|----------|----------|------------|
| Same genus | 80 | 0.01 | 300 |
| Same family | 60 | 0.05 | 200 (default) |
| Same order | 45 | 0.1 | 150 |
| Same kingdom (distant) | 30 | 1.0 | 100 |

## Multiple model_species strategy

For non-model species without close relatives, using 2-3 distantly
related model species can rescue weak annotations:

```bash
python scripts/step3_mapping.py \
    --input output/marker_weight.csv \
    --fasta target.fasta \
    --species Oryza_sativa Arabidopsis_thaliana Zea_mays \
    --pident 40 --evalue 0.1 --bitscore 150 \
    --multihomolo \
    --outdir output/
```

The pipeline automatically merges homologs from multiple species.
This is **strongly recommended** for novel species with no close
reference genome.

## Always-validate step

After picking candidates, **always run**:

```bash
python scripts/list_species.py
```

and verify each candidate appears in the output before using it in
the config. Cross-reference with `references/config-schema.md` for
the exact field format.

## Worked example: target IS a model species

**Common pitfall:** the target species is already a well-annotated
model organism (e.g., the user has an Arabidopsis h5ad and treats it as
`non_model_h5ad`). The agent may wrongly exclude the target species
from `model_species` thinking "self-mapping defeats the purpose".

**Correct decision:** include the target species in `model_species`
together with 1-2 close relatives. The xener KG has the deepest
gene->celltype edges for model organisms; self-mapping gives the highest
KG hit rate. Cross-species independence is served by adding the
relatives, not by excluding self.

Concrete config for an Arabidopsis root dataset:

```yaml
cluster_key: leiden
model_species:
  - Arabidopsis_thaliana   # target species itself - highest KG coverage
  - Brassica_rapa          # close relative - independent corroboration
non_model_fasta: abc.fasta
non_model_h5ad: edf.h5ad
organ: Root
outdir: output/edf
```

Excluding `Arabidopsis_thaliana` here typically produces 30-75% of
homologs "not in KG" and reduces the observed cell-type diversity by
half or more. The `scripts/check_output.py` post-run gate will catch
this.
