"""
results_store.py — Shared results database for the Gerizim analysis pipeline.

OVERVIEW
--------
Each analysis script writes its computed p-values, counts, and statistics
into a single JSON file (data/store/results.json) by calling:

    from lib.results_store import ResultsStore
    store = ResultsStore()
    store.write("clusterApBinom", 0.0103)      # key matches \newcommand name
    store.write("pCircAp", 0.0009)

Downstream scripts (bonferroni_correction.py, fdr_multiple_comparisons.py)
then read from the store instead of hardcoding values:

    store = ResultsStore()
    p = store.read("clusterApBinom")

The store is a flat JSON dict mapping string keys to scalar values.
Keys match the \newcommand macro names used in the manuscript, so there
is a single source of truth from computation → store → LaTeX.

STORE LOCATION
--------------
    data/store/results.json

This file is NOT committed to git (add it to .gitignore).  It is rebuilt
from scratch each time `bash manuscript/reproduce_all_macros.sh` is run.

THREAD SAFETY
-------------
Uses a file lock (via a .lock file) so parallel script runs don't corrupt
the store.  In practice scripts run sequentially, so this is just a safety net.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Location
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).parent.parent
STORE_PATH = _REPO_ROOT / "data" / "store" / "results.json"
_LOCK_PATH  = _REPO_ROOT / "data" / "store" / "results.json.lock"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

class ResultsStore:
    """
    Thin wrapper around the results JSON file.

    Usage
    -----
        store = ResultsStore()

        # Write individual values (idempotent — safe to call multiple times)
        store.write("pCircAp", 0.0009)
        store.write("NcircTotal", 90)

        # Write a dict of values at once
        store.write_many({"pCircAp": 0.0009, "NcircTotal": 90})

        # Read a value (KeyError if missing)
        p = store.read("pCircAp")

        # Read with a default
        p = store.read("pCircAp", default=None)

        # Check existence
        if store.has("pCircAp"):
            ...

        # Dump entire store
        d = store.all()
    """

    def __init__(self, path: Path | str = STORE_PATH):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    # ── Internal helpers ────────────────────────────────────────────────────

    def _load(self) -> dict:
        if not self.path.exists():
            return {}
        with open(self.path) as f:
            try:
                return json.load(f)
            except json.JSONDecodeError:
                return {}

    def _save(self, data: dict) -> None:
        tmp = self.path.with_suffix(".json.tmp")
        with open(tmp, "w") as f:
            json.dump(data, f, indent=2, sort_keys=True)
            f.write("\n")
        tmp.replace(self.path)

    def _acquire_lock(self, timeout: float = 10.0) -> None:
        deadline = time.monotonic() + timeout
        while True:
            try:
                fd = os.open(str(_LOCK_PATH), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.close(fd)
                return
            except FileExistsError:
                if time.monotonic() > deadline:
                    # Stale lock — remove and retry once
                    try:
                        _LOCK_PATH.unlink()
                    except FileNotFoundError:
                        pass
                    continue
                time.sleep(0.05)

    def _release_lock(self) -> None:
        try:
            _LOCK_PATH.unlink()
        except FileNotFoundError:
            pass

    # ── Public methods ──────────────────────────────────────────────────────

    def write(self, key: str, value: Any) -> None:
        """Write a single key → value to the store."""
        self._acquire_lock()
        try:
            data = self._load()
            data[key] = value
            self._save(data)
        finally:
            self._release_lock()

    def write_many(self, mapping: dict) -> None:
        """Write multiple key → value pairs atomically."""
        self._acquire_lock()
        try:
            data = self._load()
            data.update(mapping)
            self._save(data)
        finally:
            self._release_lock()

    _SENTINEL = object()

    def read(self, key: str, default: Any = _SENTINEL) -> Any:
        """Read a value.  Raises KeyError if missing and no default given."""
        data = self._load()
        if key not in data:
            if default is ResultsStore._SENTINEL:
                raise KeyError(
                    f"Key '{key}' not found in results store ({self.path}).\n"
                    f"Run the analysis script that produces this value first."
                )
            return default
        return data[key]

    def read_many(self, keys: list[str]) -> dict:
        """Read multiple keys at once.  Raises KeyError on any missing key."""
        data = self._load()
        missing = [k for k in keys if k not in data]
        if missing:
            raise KeyError(
                f"Keys missing from results store: {missing}\n"
                f"Run the producing analysis scripts first."
            )
        return {k: data[k] for k in keys}

    def has(self, key: str) -> bool:
        return key in self._load()

    def all(self) -> dict:
        """Return a copy of the entire store."""
        return dict(self._load())

    def clear(self) -> None:
        """Wipe the store.  Used at the start of a full pipeline run."""
        self._acquire_lock()
        try:
            self._save({})
        finally:
            self._release_lock()


# ---------------------------------------------------------------------------
# Module-level convenience functions (use the default store path)
# ---------------------------------------------------------------------------

def write(key: str, value: Any) -> None:
    ResultsStore().write(key, value)

def write_many(mapping: dict) -> None:
    ResultsStore().write_many(mapping)

def read(key: str, default: Any = None) -> Any:
    return ResultsStore().read(key, default=default)

def has(key: str) -> bool:
    return ResultsStore().has(key)
