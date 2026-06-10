#!/usr/bin/env python3
"""Run the full Xener pipeline end-to-end from a YAML config.

CLI wrapper used by the Xener agent skill. Reads config.yaml, runs all five
pipeline steps in order (markers → weights → BLAST mapping → top-k →
annotation), and writes a debug_params.yaml recording the actual parameters
used in each step.

Skill context: this is the main entry point invoked when the user wants a
one-shot annotation run. The YAML schema is documented in
skill/xener/references/config-schema.md. Validation workflow is in
references/workflows/config-validation.md — always run that workflow before
invoking this script.
"""

import argparse
import yaml
import os
import sys
from pathlib import Path

# Sibling import: scripts/ is sys.path[0] when invoked as `python scripts/...`.
# build_xener honors an optional init-config (KG endpoint / BLAST DB); with
# none it is exactly Xener(), so default runs are unchanged.
from _xener_init import build_xener, load_init_config, INIT_CONFIG_KEYS


def main():
    parser = argparse.ArgumentParser(description="Run xener full pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--init-config", default=None, metavar="PATH",
        help="Path to a SEPARATE init-config YAML describing WHERE Xener gets "
             "its data (KG_url / KG_usr / KG_pwd, blastdb_path, "
             "blastp_result_path). Optional: you may instead inline those same "
             "keys directly into --config. Omit both to use the public cloud "
             "KG + bundled BLAST database. See "
             "references/workflows/initialization.md.",
    )
    args = parser.parse_args()

    # Read the config as UTF-8 explicitly. Without an explicit encoding,
    # Python uses the platform default (e.g. GBK on zh-CN Windows), which
    # raises UnicodeDecodeError the moment a YAML comment contains a non-ASCII
    # byte such as an em-dash. Configs are authored as UTF-8, so decode as UTF-8.
    with open(args.config, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    os.makedirs(config["outdir"], exist_ok=True)

    # Resolve the init-config (WHERE Xener gets its data) from one of two
    # places, in priority order:
    #   1. a separate --init-config file, or
    #   2. KG_url / blastdb_path / ... keys inlined directly into config.yaml.
    # Either is optional; with neither, build_xener() == Xener() (cloud KG +
    # bundled BLAST DB), so existing configs run exactly as before.
    if args.init_config:
        init_config = load_init_config(args.init_config)
        print(f"Using init-config from {args.init_config}: "
              f"{ {k: ('***' if k == 'KG_pwd' and v else v) for k, v in init_config.items()} }")
    else:
        init_config = {k: config[k] for k in INIT_CONFIG_KEYS if k in config}
        if init_config:
            print(f"Using init-config inlined in {args.config}: "
                  f"{ {k: ('***' if k == 'KG_pwd' and v else v) for k, v in init_config.items()} }")
        else:
            init_config = None

    # Mirror all stdout (including xener's own logger output, which prints
    # the per-cluster "total X% homolos of organ[...] not in kg" lines the
    # quality gate depends on) to outdir/xener.log so the gate can always
    # find it. We do this by replacing sys.stdout with a tee.
    log_path = Path(config["outdir"]) / "xener.log"
    _real_stdout = sys.stdout
    log_fh = open(log_path, "w", encoding="utf-8")

    class _Tee:
        def __init__(self, *streams):
            self._streams = streams
        def write(self, s):
            for st in self._streams:
                st.write(s)
        def flush(self):
            for st in self._streams:
                try:
                    st.flush()
                except Exception:
                    pass

    sys.stdout = _Tee(_real_stdout, log_fh)

    try:
        annor = build_xener(init_config)
        cluster2celltype, cluster2max, debug_params = annor.run_from_yaml(args.config)

        print("Pipeline complete.")
        print("Cluster annotations:", cluster2celltype)
        print("Cluster max-init-weight cell types:", cluster2max)
    finally:
        sys.stdout = _real_stdout
        log_fh.close()

    # Persist a LIGHTWEIGHT annotation artifact (per-cell UMAP coordinates +
    # label columns) instead of a full h5ad copy, to avoid duplicating a
    # multi-GB expression matrix on disk. scripts/plot_umap.py reconstructs a
    # minimal AnnData from this CSV. Neither run_from_yaml nor refine_cluster.py
    # persists these labels, so this is the single consolidation point for
    # 'xener' / 'xener_refine'.
    try:
        import scanpy as sc
        import pandas as pd
        from glob import glob

        h5ad_in = config["non_model_h5ad"]
        cluster_key = config["cluster_key"]
        adata = sc.read(h5ad_in)

        c2c = {str(k): v for k, v in cluster2celltype.items()}
        c2m = {str(k): v for k, v in cluster2max.items()}
        clusters = adata.obs[cluster_key].astype(str)

        out = pd.DataFrame(index=adata.obs_names)
        out[cluster_key] = clusters.values
        out["xener"] = clusters.map(c2c).astype(str).values
        out["xener_max"] = clusters.map(c2m).astype(str).values

        # UMAP coordinates, so plotting needs no access to the source h5ad.
        if "X_umap" in adata.obsm:
            umap = adata.obsm["X_umap"]
            out["UMAP_1"] = umap[:, 0]
            out["UMAP_2"] = umap[:, 1]
        else:
            print("[WARN] X_umap not in adata.obsm; annotation CSV will have "
                  "no UMAP coordinates.", file=sys.stderr)

        # Merge refinement: per-cell subtype label. This column carries ONLY
        # the refinement result - a cell is labeled iff its cluster was
        # refined; cells in unrefined clusters stay EMPTY (NaN), not copied
        # from `xener`. Conflating the two would silently overstate how much of
        # the dataset was actually split. Coalesce downstream if you need a
        # fully-populated column: out["xener_refine"].fillna(out["xener"]).
        refine_dir = Path(config["outdir"]) / "refine_output"
        refine_csvs = sorted(glob(str(refine_dir / "refined_*.csv")))
        if refine_csvs:
            refine_col = pd.Series(
                pd.NA, index=out.index, dtype=object
            )
            n_refined_cells = 0
            for csv in refine_csvs:
                df = pd.read_csv(csv, index_col=0)
                if "xener_refine" not in df.columns:
                    continue
                common = out.index.intersection(df.index)
                refine_col.loc[common] = df.loc[common, "xener_refine"].astype(str)
                n_refined_cells += len(common)
            out["xener_refine"] = refine_col
            n_unrefined = int(out["xener_refine"].isna().sum())
            print(f"Merged {len(refine_csvs)} refinement CSV(s) into "
                  f"xener_refine ({n_refined_cells} cells labeled, "
                  f"{n_unrefined} left empty = not refined).")

        basename = Path(h5ad_in).stem
        annot_path = Path(config["outdir"]) / f"{basename}_annotation.csv"
        out.to_csv(annot_path)
        print(f"Annotation artifact written to {annot_path}")
    except Exception as e:
        print(f"[WARN] Failed to persist annotation artifact: {e}", file=sys.stderr)

    # Mandatory post-run quality gate. If this fails, the pipeline exits
    # non-zero and the run is not considered successful. The agent must
    # inspect the gate output, adjust the config (most commonly
    # model_species), and re-run from step 3.
    import subprocess
    h5ad_path = config.get("non_model_h5ad")
    gate = subprocess.run(
        [sys.executable, str(Path(__file__).parent / "check_output.py"),
         "--outdir", config["outdir"], *(["--h5ad", h5ad_path] if h5ad_path else [])],
        capture_output=True, text=True,
    )
    print(gate.stdout, end="")
    if gate.returncode != 0:
        print(gate.stderr, end="", file=sys.stderr)
        sys.exit(gate.returncode)


if __name__ == "__main__":
    main()