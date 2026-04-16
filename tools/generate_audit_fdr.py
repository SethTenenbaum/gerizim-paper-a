"""
generate_audit_fdr.py
=====================
Produce supplementary/audit/fdr.txt

Wrapper around analysis/unesco/fdr_multiple_comparisons.py that:
  1. Captures all stdout from the FDR script
  2. Strips the LaTeX \\newcommand lines at the end (those belong in
     generated_macros.tex, not in the supplementary audit file)
  3. Writes the clean audit text to supplementary/audit/fdr.txt

Run from repo root:
    python3 tools/generate_audit_fdr.py
"""

import io
import sys
import contextlib
from pathlib import Path
from datetime import datetime, timezone

# Ensure repo root is importable
_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_ROOT))

OUT = _ROOT / "supplementary" / "audit" / "fdr.txt"


def run():
    # Capture stdout from the FDR script
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        # Import and call main() directly so we share the same process
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "fdr_multiple_comparisons",
            _ROOT / "analysis" / "unesco" / "fdr_multiple_comparisons.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.main()

    raw = buf.getvalue()

    # Strip LaTeX macro lines — everything from the first \newcommand onwards
    # that appears after the FDR table content.
    # Strategy: drop any line that starts with optional whitespace + \newcommand
    # or is a % LaTeX macros comment, but only after the interpretation section.
    lines_in   = raw.splitlines()
    lines_out  = []
    in_macros  = False

    for line in lines_in:
        stripped = line.strip()
        # Detect the start of the macro block
        if stripped.startswith("% LaTeX macros") or stripped.startswith("% Per-test BH"):
            in_macros = True
        if in_macros:
            # Keep blank lines only if they close a previous content block,
            # but once we are in macro territory, stop.
            continue
        lines_out.append(line)

    # Trim trailing blank lines, then add a clean footer
    while lines_out and not lines_out[-1].strip():
        lines_out.pop()

    ts = datetime.now(timezone.utc).strftime("%a %b %d %H:%M:%S UTC %Y")
    lines_out += [
        "",
        "─" * 110,
        f"Generated : {ts}",
        f"Script    : tools/generate_audit_fdr.py",
        "LaTeX macros for this analysis are in manuscript/generated_macros.tex",
        "",
    ]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines_out), encoding="utf-8")
    print(f"Written -> {OUT}")


if __name__ == "__main__":
    run()
