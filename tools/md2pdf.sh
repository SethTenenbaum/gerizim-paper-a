#!/usr/bin/env bash
set -euo pipefail

# md2pdf.sh
# Simple markdown to PDF converter using pandoc and a TeX engine.
# Usage: ./tools/md2pdf.sh [IN.md] [OUT.pdf]

INFILE="${1:-guide/statistical_methods_guide.md}"
OUTFILE="${2:-${INFILE%.md}.pdf}"

if [ ! -f "$INFILE" ]; then
  echo "Input file not found: $INFILE" >&2
  echo "Usage: $0 [IN.md] [OUT.pdf]" >&2
  exit 2
fi

# Check for pandoc
if ! command -v pandoc >/dev/null 2>&1; then
  cat <<EOF >&2
Pandoc is not installed. Install it with Homebrew:
  brew install pandoc
or get it from https://pandoc.org/installing.html
If you do not have Homebrew, see https://brew.sh
EOF
  exit 3
fi

# Prefer xelatex for good font handling; fall back to pdflatex
PDF_ENGINE=""
if command -v xelatex >/dev/null 2>&1; then
  PDF_ENGINE="xelatex"
elif command -v pdflatex >/dev/null 2>&1; then
  PDF_ENGINE="pdflatex"
else
  cat <<EOF >&2
No LaTeX engine found. Install BasicTeX or MacTeX. With Homebrew you can install BasicTeX via:
  brew install --cask basictex
After installing BasicTeX you may need to install pandoc dependencies. See https://tug.org/mactex/
EOF
  exit 4
fi

# Run pandoc
pandoc "$INFILE" \
  --from markdown+yaml_metadata_block+footnotes+link_attributes+smart \
  --pdf-engine="$PDF_ENGINE" \
  --toc --toc-depth=2 --number-sections \
  -V geometry:margin=1in \
  -V lang:en-US \
  --highlight-style=tango \
  -o "$OUTFILE"

echo "Wrote $OUTFILE"
