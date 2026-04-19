"""
generate_audit_fine_sweep.py
============================
Produce supplementary/audit/fine_sweep_audit.txt

Runs the fine unit sweep (±1% of canonical 0.10 beru) for both the
dome/spherical population and the full UNESCO corpus, and writes the
results table to the supplementary audit directory.

Run from repo root:
    python3 tools/generate_audit_fine_sweep.py
"""

import io
import sys
import contextlib
from pathlib import Path
from datetime import datetime, timezone

_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

OUT = _ROOT / "supplementary" / "audit" / "fine_sweep_audit.txt"


def run():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "fine_sweep_audit",
            _ROOT / "analysis" / "unesco" / "fine_sweep_audit.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.main()

    raw = buf.getvalue()

    lines_out = raw.splitlines()

    while lines_out and not lines_out[-1].strip():
        lines_out.pop()

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines_out += [
        "",
        "─" * 80,
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_fine_sweep.py",
        "Source    : analysis/unesco/fine_sweep_audit.py",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines_out), encoding="utf-8")
    print(f"Written -> {OUT}")


if __name__ == "__main__":
    run()
