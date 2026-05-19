#!/usr/bin/env python3
"""Install xener from PyPI."""

import subprocess
import sys


def main():
    print("Installing xener from PyPI...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "xener"])
    print("xener installed successfully.")


if __name__ == "__main__":
    main()