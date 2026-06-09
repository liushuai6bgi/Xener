#!/usr/bin/env python3
"""Install xener from PyPI.

CLI wrapper used by the Xener agent skill. Equivalent to `pip install xener`
but invoked as a script so the agent does not need to call pip directly.

Skill context: first step in the workflow. Run this (or `pip install xener`)
before any other script. Verify with `python -c "import xener"` (the only
inline import allowed by skill/xener/references/mandatory-rules.md).
"""

import subprocess
import sys


def main():
    print("Installing xener from PyPI...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "xener"])
    print("xener installed successfully.")


if __name__ == "__main__":
    main()