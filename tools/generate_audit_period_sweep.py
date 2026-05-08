"""
generate_audit_period_sweep.py
==============================
Produce supplementary/audit/period_sweep_audit.txt

Runs the full period sweep dissociation analysis (1.5° to 6.0° in 0.1°
steps, 46 periods) and writes the complete table — all periods, both
Rayleigh and binomial proximity results — to the supplementary audit
directory.

This is the full sweep underlying Table 1 of the paper. The four periods
shown in Table 1 (2.4°, 3.0°, 3.5°, 3.6°) are flagged with a marker.

Run from repo root:
    python3 tools/generate_audit_period_sweep.py
"""

import io
import sys
import contextlib
from pathlib import Path
from datetime import datetime, timezone

_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

OUT = _ROOT / "supplementary" / "audit" / "period_sweep_audit.txt"


def run():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "period_sweep_dissociation",
            _ROOT / "analysis" / "unesco" / "period_sweep_dissociation.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)

    raw = buf.getvalue()

    lines_out = raw.splitlines()

    # Strip trailing blank lines
    while lines_out and not lines_out[-1].strip():
        lines_out.pop()

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines_out += [
        "",
        "─" * 80,
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_period_sweep.py",
        "Source    : analysis/unesco/period_sweep_dissociation.py",
        "Note      : Four periods marked '<-- Table 1' are the rows shown in",
        "            Table 1 of the manuscript (reframed_paper.tex).",
        "            Tier-A threshold = ±0.10 × T; geometric null p0 = 0.20.",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines_out), encoding="utf-8")
    print(f"Written -> {OUT}")


if __name__ == "__main__":
    run()
