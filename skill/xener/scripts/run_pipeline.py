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

from xener import Xener


def main():
    parser = argparse.ArgumentParser(description="Run xener full pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    os.makedirs(config["outdir"], exist_ok=True)

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
                    pass)

    sys.stdout = _Tee(_real_stdout, log_fh)

    try:
        annor = Xener()
        cluster2celltype, cluster2max, debug_params = annor.run_from_yaml(args.config)

        print("Pipeline complete.")
        print("Cluster annotations:", cluster2celltype)
        print("Cluster max-init-weight cell types:", cluster2max)
    finally:
        sys.stdout = _real_stdout
        log_fh.close()

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