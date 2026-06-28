#!/usr/bin/env bash
# Build Word + LaTeX/PDF manuscript with embedded numbered figures.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPORTS="$ROOT/reports"
DRAFT="$REPORTS/manuscript_with_figures.md"
EQ_FRAGMENT="$REPORTS/manuscript_equations_fragment.tex"
TEX_OUT="$REPORTS/manuscript.tex"
DOCX_OUT="$REPORTS/manuscript.docx"
PDF_OUT="$REPORTS/manuscript.pdf"

python3 "$ROOT/scripts/build_manuscript_with_figures.py"

if [[ ! -f "$DRAFT" ]]; then
  echo "Missing $DRAFT"
  exit 1
fi

if command -v pandoc >/dev/null 2>&1; then
  echo "Building LaTeX via pandoc..."
  pandoc "$DRAFT" \
    -o "$TEX_OUT" \
    --standalone \
    --metadata title="Dosimetric-Response Modeling in Glioblastoma (CFB-GBM)" \
    --metadata author="TCP Modeling GBM Project" \
    --variable geometry:margin=2.5cm \
    --variable fontsize=11pt \
    --from markdown-smart \
    --resource-path="$ROOT:$ROOT/figures"

  if [[ -f "$EQ_FRAGMENT" ]]; then
    python3 - <<'PY'
from pathlib import Path
tex = Path("reports/manuscript.tex")
frag = Path("reports/manuscript_equations_fragment.tex")
if tex.exists() and frag.exists():
    text = tex.read_text()
    marker = "\\section{2. Materials and Methods}"
    insert = (
        "\n\\subsection{Mathematical models (equations)}\n"
        "\\input{manuscript_equations_fragment.tex}\n"
    )
    if marker in text and "manuscript_equations_fragment" not in text:
        text = text.replace(marker, marker + insert, 1)
        tex.write_text(text)
PY
  fi

  echo "Building Word via pandoc..."
  pandoc "$DRAFT" \
    -o "$DOCX_OUT" \
    --from markdown-smart \
    --resource-path="$ROOT:$ROOT/figures"

  if command -v pdflatex >/dev/null 2>&1; then
    echo "Building PDF via pdflatex..."
    (cd "$REPORTS" && pdflatex -interaction=nonstopmode manuscript.tex >/dev/null && \
     pdflatex -interaction=nonstopmode manuscript.tex >/dev/null) || true
    [[ -f "$PDF_OUT" ]] && echo "Wrote $PDF_OUT" || echo "PDF build skipped (pdflatex warnings)"
  else
    echo "pdflatex not found — LaTeX and DOCX written; compile PDF manually."
  fi
  echo "Wrote $TEX_OUT"
  echo "Wrote $DOCX_OUT"
else
  echo "pandoc not found. Install: brew install pandoc"
  exit 1
fi
