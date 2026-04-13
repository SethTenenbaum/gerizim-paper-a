#!/usr/bin/env python3
"""
expand_and_convert.py
Expands all \newcommand macros from generated_macros.tex into the manuscript,
then converts the result to a Word .docx via pandoc.

Usage (from repo root):
    python tools/expand_and_convert.py
"""

import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).parent.parent
MANUSCRIPT = REPO / "manuscript" / "paper_a_primary_unesco.tex"
MACROS_TEX = REPO / "manuscript" / "generated_macros.tex"
LOCAL_MACROS = REPO / "manuscript" / "paper_a_primary_unesco.tex"  # for \newcommand{\degr}
OUTPUT_DOCX = REPO / "manuscript" / "paper_a_primary_unesco.docx"


def parse_macros(tex_path: Path) -> dict[str, str]:
    """Return {macro_name: replacement_text} from a file of \\newcommand definitions.
    Handles values that contain nested braces (e.g. 62{,}785 thousands separators)."""
    macros = {}
    text = tex_path.read_text(encoding="utf-8")
    # Find \newcommand{\Name} then extract the balanced-brace value
    header = re.compile(r'\\newcommand\s*\{\\([A-Za-z@]+)\}\s*')
    pos = 0
    while pos < len(text):
        m = header.search(text, pos)
        if not m:
            break
        name = m.group(1)
        # Find the opening { of the value
        i = m.end()
        while i < len(text) and text[i] != '{':
            if text[i] == '\n':
                break  # no value on this line
            i += 1
        if i >= len(text) or text[i] != '{':
            pos = m.end()
            continue
        # Extract balanced braces
        depth = 0
        start = i
        for j in range(i, len(text)):
            if text[j] == '{':
                depth += 1
            elif text[j] == '}':
                depth -= 1
                if depth == 0:
                    value = text[start+1:j]
                    macros[name] = value
                    break
        pos = m.end()
    return macros


def expand_macros(source: str, macros: dict[str, str], max_passes: int = 6) -> str:
    """
    Iteratively substitute \\MacroName with its value.
    Multiple passes handle macros whose values reference other macros.
    Longest-name-first ordering prevents partial matches.
    """
    for _ in range(max_passes):
        prev = source
        for name in sorted(macros, key=len, reverse=True):
            repl = macros[name]
            source = re.sub(r'\\' + re.escape(name) + r'(?![A-Za-z@])',
                            lambda m, r=repl: r, source)
        if source == prev:
            break
    return source


def main() -> None:
    print("Parsing macros …")
    macros = parse_macros(MACROS_TEX)
    # Also grab any \newcommand lines in the manuscript header (e.g. \degr)
    macros.update(parse_macros(MANUSCRIPT))
    print(f"  {len(macros)} macros loaded")

    print("Expanding macros in manuscript …")
    source = MANUSCRIPT.read_text(encoding="utf-8")

    # Remove the \input{generated_macros} line — values are inlined below
    source = re.sub(r'\\input\{generated_macros\}\s*\n?', '', source)

    # Expand macros (skip 'degr' — its value \textdegree{} causes a regex
    # collision with \textdegree already in the source)
    expand_macros_map = {k: v for k, v in macros.items() if k != 'degr'}
    expanded = expand_macros(source, expand_macros_map)

    # Leave \newcommand lines in place — pandoc tolerates them and removing
    # them risks eating adjacent braces (e.g. {,} thousands separators).

    # Replace LaTeX thousands-separator {,} with plain comma for pandoc/Word
    expanded = expanded.replace('{,}', ',')

    # Write to a temp file next to the manuscript so relative \includegraphics paths work
    tmp = MANUSCRIPT.parent / "_expanded_tmp.tex"
    tmp.write_text(expanded, encoding="utf-8")
    print(f"  Written expanded source → {tmp.name}")

    print("Running pandoc …")
    result = subprocess.run(
        [
            "pandoc", str(tmp),
            "--from=latex",
            "--to=docx",
            f"--output={OUTPUT_DOCX}",
            "--resource-path=.",
            "--wrap=none",
        ],
        cwd=MANUSCRIPT.parent,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print("pandoc stderr:\n", result.stderr)
        print(f"  (expanded source preserved at {tmp} for inspection)")
        sys.exit(result.returncode)

    tmp.unlink()  # clean up temp file only on success

    if result.stderr.strip():
        print("pandoc warnings:\n", result.stderr.strip())

    size_kb = OUTPUT_DOCX.stat().st_size // 1024
    print(f"\n✅  {OUTPUT_DOCX.relative_to(REPO)}  ({size_kb} KB)")


if __name__ == "__main__":
    main()
