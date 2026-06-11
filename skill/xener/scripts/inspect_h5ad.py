#!/usr/bin/env python3
"""One-shot inspection of a single-cell .h5ad for the Xener agent skill.

CLI wrapper used during references/workflows/inspection.md. Loads the h5ad
**exactly once** and prints everything the agent needs to compose a config:
cluster_key, species hint, organ hint, cluster sizes, a recommended top_num,
and an .X / raw sanity check.

Why a dedicated script: the target h5ad is frequently multi-GB, and reading it
is the single dominant I/O cost of the whole workflow. The old guidance was a
loose inline `python -c` snippet, which in practice led the agent to read the
file 2-3 separate times (once for columns, again for organ value-counts, again
for cluster sizes / .X). This script folds all of that into ONE load so no
follow-up read is ever required. See mandatory-rules.md sec.10 (single
inspection / I/O discipline).

Mandatory-rules note: this is the inspection analogue of the "only inline
scanpy import allowed" exception (mandatory-rules.md sec.1). It imports scanpy,
never `import xener`, and needs no --init-config (pure h5ad inspection, no KG /
BLAST).

Usage:
    python scripts/inspect_h5ad.py edf.h5ad
    python scripts/inspect_h5ad.py edf.h5ad --json     # machine-readable
"""

import argparse
import json
import sys


# Species-hint detection by var_names prefix (mirrors the table in
# references/workflows/inspection.md). First match wins.
SPECIES_HINTS = [
    ("AT", "Arabidopsis_thaliana"),       # AT1G01010 etc. (TAIR locus IDs)
    ("Zm", "Zea_mays"),                   # Zm00001 / Zm...
    ("LOC_Os", "Oryza_sativa"),
    ("Os", "Oryza_sativa"),
    ("Bra", "Brassica_rapa"),
    ("Brara", "Brassica_rapa"),
    ("Glyma.", "Glycine_max"),
    ("Medtr", "Medicago_truncatula"),
    ("MTR_", "Medicago_truncatula"),
]


def detect_species(var_names):
    """Return (detected_species_or_None, fraction_matched, matched_prefix)."""
    names = [str(v) for v in var_names]
    n = len(names)
    if n == 0:
        return None, 0.0, None
    for prefix, species in SPECIES_HINTS:
        hits = sum(1 for g in names if g.startswith(prefix))
        if hits / n > 0.5:
            return species, hits / n, prefix
    return None, 0.0, None


def recommend_cluster_key(obs):
    """Pick a clustering column: leiden > louvain > *cluster*, skipping
    existing annotation columns. Returns (key_or_None, cardinality)."""
    cols = list(obs.columns)
    lower = {c: c.lower() for c in cols}

    def card(c):
        try:
            return int(obs[c].nunique())
        except Exception:
            return -1

    # Skip columns that are clearly an existing annotation (don't re-annotate).
    def is_annotation(c):
        lc = lower[c]
        return ("celltype" in lc) or ("annotation" in lc) or (lc == "cell_type")

    for want in ("leiden", "louvain"):
        for c in cols:
            if lower[c] == want and 2 <= card(c) <= 100:
                return c, card(c)
    # Any column containing "cluster" (e.g. Seurat_clusters) with sane cardinality.
    candidates = [
        c for c in cols
        if "cluster" in lower[c] and not is_annotation(c) and 2 <= card(c) <= 100
    ]
    if candidates:
        # Prefer the lowest-cardinality clustering (more aggregated).
        candidates.sort(key=card)
        return candidates[0], card(candidates[0])
    return None, -1


def recommend_top_num(median_size):
    """Documented heuristic (inspection.md): smaller clusters need more
    candidate genes; larger clusters can afford fewer."""
    if median_size < 50:
        return 50
    if median_size < 200:
        return 30
    if median_size < 1000:
        return 20
    return 15


def detect_organ(obs, uns):
    """Best-effort organ hint from obs Organ/Tissue value_counts or uns."""
    hints = {}
    for col in ("Organ", "organ", "Tissue", "tissue"):
        if col in obs.columns:
            vc = obs[col].value_counts()
            hints[col] = {str(k): int(v) for k, v in vc.head(5).items()}
    for key in ("organism", "species", "organ", "tissue"):
        if key in uns:
            try:
                hints[f"uns[{key}]"] = str(uns[key])
            except Exception:
                pass
    return hints


def inspect(path):
    import numpy as np
    import scanpy as sc
    from scipy.sparse import issparse

    adata = sc.read(path)  # the ONE load
    obs = adata.obs
    uns = adata.uns

    # obs cardinality (likely-clustering range)
    cardinality = {}
    for col in obs.columns:
        try:
            ncat = int(obs[col].nunique())
        except Exception:
            continue
        if 2 <= ncat <= 100:
            cardinality[col] = ncat

    cluster_key, ck_card = recommend_cluster_key(obs)

    # cluster sizes for the recommended key
    sizes_summary = {}
    if cluster_key is not None:
        vc = obs[cluster_key].astype(str).value_counts()
        med = float(np.median(vc.values))
        sizes_summary = {
            "n_clusters": int(vc.shape[0]),
            "min": int(vc.min()),
            "median": int(med),
            "mean": round(float(vc.mean()), 1),
            "max": int(vc.max()),
        }
    top_num = recommend_top_num(sizes_summary.get("median", 200)) if sizes_summary else 30

    # species hint
    species, frac, prefix = detect_species(adata.var_names)

    # .X sanity (small slice only)
    X = adata.X
    sub = X[:200].toarray() if issparse(X) else np.asarray(X[:200])
    x_min, x_max = float(sub.min()), float(sub.max())
    looks_lognorm = bool(x_max < 50 and x_min >= 0)

    report = {
        "path": str(path),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_columns": list(obs.columns),
        "obs_cardinality_2_100": cardinality,
        "recommended_cluster_key": cluster_key,
        "recommended_cluster_key_cardinality": ck_card,
        "uns_keys": list(uns.keys()),
        "obsm_keys": list(adata.obsm.keys()),
        "has_X_umap": "X_umap" in adata.obsm,
        "raw_set": adata.raw is not None,
        "var_names_sample": [str(v) for v in adata.var_names[:25]],
        "species_hint": species,
        "species_hint_fraction": round(frac, 3),
        "species_hint_prefix": prefix,
        "cluster_sizes": sizes_summary,
        "recommended_top_num": top_num,
        "organ_hints": detect_organ(obs, uns),
        "X_check": {
            "sparse": bool(issparse(X)),
            "dtype": str(X.dtype),
            "min": round(x_min, 4),
            "max": round(x_max, 4),
            "looks_lognorm": looks_lognorm,
        },
    }
    return report


def print_human(r):
    print(f"=== {r['path']} ===")
    print(f"shape: {r['n_obs']} cells x {r['n_vars']} genes")
    print()
    print("=== clustering ===")
    print(f"recommended cluster_key: {r['recommended_cluster_key']} "
          f"({r['recommended_cluster_key_cardinality']} clusters)")
    print(f"obs columns (cardinality 2..100): {r['obs_cardinality_2_100']}")
    if r["cluster_sizes"]:
        s = r["cluster_sizes"]
        print(f"cluster sizes: n={s['n_clusters']} min={s['min']} "
              f"median={s['median']} mean={s['mean']} max={s['max']}")
    print(f"=> recommended top_num: {r['recommended_top_num']}")
    print()
    print("=== species hint ===")
    if r["species_hint"]:
        print(f"{r['species_hint']}  (prefix '{r['species_hint_prefix']}', "
              f"{r['species_hint_fraction']:.0%} of var_names)")
    else:
        print("none detected from var_names prefix; decide from metadata / prompt")
    print(f"var_names sample: {r['var_names_sample'][:10]}")
    print()
    print("=== organ hints ===")
    print(r["organ_hints"] if r["organ_hints"] else "none detected")
    print()
    print("=== embedding / raw / .X ===")
    print(f"obsm keys: {r['obsm_keys']}  | X_umap present: {r['has_X_umap']}")
    print(f"raw set: {r['raw_set']}")
    xc = r["X_check"]
    print(f".X: sparse={xc['sparse']} dtype={xc['dtype']} "
          f"range=[{xc['min']}, {xc['max']}] looks_lognorm={xc['looks_lognorm']}")
    print()
    print("Next: compose config.yaml from these values (see "
          "references/workflows/inspection.md). Do NOT read the h5ad again - "
          "every field needed is above (mandatory-rules.md sec.10).")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("h5ad", help="Path to the .h5ad file to inspect")
    ap.add_argument("--json", action="store_true",
                    help="Emit the report as JSON instead of human-readable text")
    args = ap.parse_args()

    try:
        report = inspect(args.h5ad)
    except Exception as e:
        print(f"[ERROR] failed to inspect {args.h5ad}: {e}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print_human(report)


if __name__ == "__main__":
    main()
