#!/usr/bin/env bash
# render.sh — build assets, then render the RevealJS deck (self-contained HTML).
# Requires: quarto on PATH (RevealJS needs no LaTeX).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

echo ">> building assets"
bash scripts/prepare-assets.sh

echo ">> rendering RevealJS (HTML)"
quarto render slides.qmd --to revealjs

echo ">> done. Output: _output/slides.html"
ls -la _output/slides.html 2>/dev/null || true
