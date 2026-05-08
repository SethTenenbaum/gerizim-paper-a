"""
generate_audit_null_b.py
========================
Produce supplementary/audit/null_b_uniform_draw_audit.txt

Runs the Null B uniform random draw test (analysis/unesco/null_b_uniform_draw.py)
and writes the full output to the supplementary audit directory.

Null B asks: if N longitudes were drawn uniformly at random from [−180°, +180°],
how often would the resulting sample produce at least as many harmonic-proximate
sites as the observed dome/stupa corpus?  This is the complement to Null A
(within-dome bootstrap, which conditions on observed clustering).

Run from repo root:
    python3 tools/generate_audit_null_b.py
"""

import io
import sys
import contextlib
from pathlib import Path
from datetime import datetime, timezone

_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

OUT = _ROOT / "supplementary" / "audit" / "null_b_uniform_draw_audit.txt"


def run():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "null_b_uniform_draw",
            _ROOT / "analysis" / "unesco" / "null_b_uniform_draw.py",
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
        f"Script    : tools/generate_audit_null_b.py",
        "Source    : analysis/unesco/null_b_uniform_draw.py",
        "Note      : Null B complements Null A (within-dome bootstrap).",
        "            Null A (ns) — hit count reproducible from dome geography.",
        "            Null B (**) — hit count not reproducible from uniform random.",
        "            The two nulls address orthogonal questions; see manuscript",
        "            Section 5.3 (Negative Controls) and Section 6.2 (Discussion).",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines_out), encoding="utf-8")
    print(f"Written -> {OUT}")


if __name__ == "__main__":
    run()
