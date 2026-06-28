#!/usr/bin/env bash
# Build assignment-style formal report (Word + PDF) with embedded figures.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPORTS="$ROOT/reports"
SRC="$REPORTS/assignment_report.md"
DOCX="$REPORTS/assignment_report.docx"
PDF="$REPORTS/assignment_report.pdf"
TEX="$REPORTS/assignment_report.tex"

if [[ ! -f "$SRC" ]]; then
  echo "Missing $SRC"
  exit 1
fi

if ! command -v pandoc >/dev/null 2>&1; then
  echo "pandoc not found"
  exit 1
fi

echo "Building assignment report Word..."
pandoc "$SRC" \
  -o "$DOCX" \
  --from markdown-smart \
  --resource-path="$ROOT:$ROOT/figures" \
  --metadata title="Final Project Report: TCP Modeling in Glioblastoma (CFB-GBM)" \
  --metadata author="TCP Modeling GBM Project Group"

echo "Building assignment report LaTeX..."
pandoc "$SRC" \
  -o "$TEX" \
  --standalone \
  --from markdown-smart \
  --resource-path="$ROOT:$ROOT/figures" \
  --metadata title="Final Project Report: TCP Modeling in Glioblastoma (CFB-GBM)" \
  --metadata author="TCP Modeling GBM Project Group" \
  --variable geometry:margin=2.5cm \
  --variable fontsize=11pt

if command -v pdflatex >/dev/null 2>&1; then
  echo "Building assignment report PDF..."
  (cd "$REPORTS" && pdflatex -interaction=nonstopmode assignment_report.tex >/dev/null && \
   pdflatex -interaction=nonstopmode assignment_report.tex >/dev/null) || true
  [[ -f "$PDF" ]] && echo "Wrote $PDF" || echo "PDF build had warnings — check $REPORTS/assignment_report.log"
else
  echo "pdflatex not found — DOCX and TEX written."
fi

echo "Wrote $DOCX"
echo "Wrote $TEX"
