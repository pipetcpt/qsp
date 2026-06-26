#!/usr/bin/env bash
# render.sh — build assets, then render both decks (RevealJS HTML + Beamer PDF).
# Requires: quarto on PATH. RevealJS needs no LaTeX; Beamer needs a LaTeX engine
# (install once with:  quarto install tinytex).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

echo ">> building assets"
bash scripts/prepare-assets.sh

echo ">> rendering RevealJS (HTML)"
quarto render slides.qmd --to revealjs

echo ">> rendering Beamer (PDF)"
quarto render slides.qmd --to beamer || {
  echo "!! Beamer render failed (LaTeX missing?). Try: quarto install tinytex"
  exit 1
}

echo ">> done. Outputs in _output/"
ls -la _output/ 2>/dev/null || true
