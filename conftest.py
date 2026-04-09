"""
conftest.py — Ensures the project root is on sys.path for all pytest runs
and direct script executions from any working directory.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
