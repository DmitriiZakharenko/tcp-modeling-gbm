#!/usr/bin/env bash
# Build assignment-style formal report (Word + PDF) with embedded figures.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPORTS="$ROOT/reports"
FIG_STAGE="$REPORTS/figures"
SRC="$REPORTS/assignment_report.md"
DOCX="$REPORTS/assignment_report.docx"
PDF="$REPORTS/assignment_report.pdf"
TEX="$REPORTS/assignment_report.tex"
HEADER="$REPORTS/assignment_report_header.tex"

if [[ ! -f "$SRC" ]]; then
  echo "Missing $SRC"
  exit 1
fi

if ! command -v pandoc >/dev/null 2>&1; then
  echo "pandoc not found"
  exit 1
fi

# Stage figures next to LaTeX output (pdflatex/xelatex look in reports/figures/)
mkdir -p "$FIG_STAGE"
FIG_LIST=(
  03_tcp_curves_os_proxy.png
  04_clinical_prognosis.png
  05_rano_vs_os_tcp_auc.png
  06_within_arm_rano_tcp.png
  07_rano_logistic_roc_40gy.png
  07_rano_volume_validation_40gy.png
  08_pooled_rano_roc.png
  08_pyradiomics_vs_volume_auc.png
  08_pyradiomics_nested_cv_auc.png
)
for f in "${FIG_LIST[@]}"; do
  cp -f "$ROOT/figures/$f" "$FIG_STAGE/$f"
done
echo "Staged ${#FIG_LIST[@]} figures in $FIG_STAGE"

PANDOC_COMMON=(
  "$SRC"
  --from markdown-smart
  --standalone
  --toc
  --number-sections
  --resource-path="$REPORTS:$FIG_STAGE:$ROOT/figures"
  --metadata title="Final Project Report: TCP Modeling in Glioblastoma (CFB-GBM)"
  --metadata author="TCP Modeling GBM Project Group"
  --metadata date="2026-06-28"
)

echo "Building assignment report Word..."
pandoc "${PANDOC_COMMON[@]}" -o "$DOCX"

echo "Building assignment report LaTeX..."
pandoc "${PANDOC_COMMON[@]}" \
  -o "$TEX" \
  --include-in-header="$HEADER" \
  -V geometry:margin=2.5cm \
  -V fontsize=11pt \
  -V linestretch=1.12 \
  -V documentclass=article

PDF_ENGINE="xelatex"
if ! command -v xelatex >/dev/null 2>&1; then
  PDF_ENGINE="pdflatex"
fi

echo "Building assignment report PDF via pandoc ($PDF_ENGINE)..."
pandoc "${PANDOC_COMMON[@]}" \
  -o "$PDF" \
  --pdf-engine="$PDF_ENGINE" \
  --include-in-header="$HEADER" \
  -V geometry:margin=2.5cm \
  -V fontsize=11pt \
  -V linestretch=1.12 \
  -V documentclass=article

echo "Wrote $DOCX"
echo "Wrote $TEX"
echo "Wrote $PDF"
