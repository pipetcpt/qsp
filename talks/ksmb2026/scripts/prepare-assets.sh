#!/usr/bin/env bash
# prepare-assets.sh — build every figure the slides reference.
#
# Primary path (used here): a pure-Python builder (cairosvg + Pillow + matplotlib
# + numpy) that needs NO R / Graphviz / ImageMagick:
#     pip install cairosvg pillow matplotlib numpy
#     bash prepare-assets.sh
#
# It (1) rasterizes hero maps from each disease's small SVG, (2) tiles a 24-map
# montage, (3) draws the pipeline & CCR-loop diagrams, (4) integrates reduced
# illustrative ODE plots, (5) plots library-wide distributions, and (6) writes
# labeled placeholders for Shiny screenshots and institutional logos.
#
# Alternatives if you have the classic toolchain:
#   maps:    rsvg-convert -w 1600 in.svg -o out.png   (or:  dot -Tpng in.dot)
#   montage: montage <pngs> -tile 6x4 -geometry +6+6 montage.png   (ImageMagick)
#   plots:   Rscript make-figs.R                       (real mrgsolve simulations)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# Refresh the library statistics CSV + _variables.yml first.
bash "$HERE/library-stats.sh" > "$HERE/../_variables.yml" 2>/dev/null || true
cat >> "$HERE/../_variables.yml" <<'YAML'
n_categories: 15
max_nodes: 685
talk_date: "2026"
YAML

# (Re)generate the per-model statistics CSV consumed by the meta-analysis figure.
ROOT="$(cd "$HERE/../../.." && pwd)"
CSV="$HERE/../assets/plots/library_stats.csv"
mkdir -p "$(dirname "$CSV")"
( cd "$ROOT"
  set +e   # grep -c exits 1 on zero matches; that is not an error here
  echo "model,clusters,refs,nodes" > "$CSV"
  for d in */; do
    d=${d%/}; [ "$d" = "talks" ] && continue
    dot=$(ls "$d"/*_qsp*.dot 2>/dev/null | head -1); [ -z "$dot" ] && continue
    ref=$(ls "$d"/*references.md 2>/dev/null | head -1)
    cl=$(grep -c 'subgraph cluster' "$dot"); cl=${cl:-0}
    nodes=$(grep -cE '\[label=|shape=' "$dot"); nodes=${nodes:-0}
    if [ -n "$ref" ]; then rf=$(grep -coE 'https://pubmed' "$ref"); else rf=0; fi
    printf '%s,%s,%s,%s\n' "$d" "$cl" "$rf" "$nodes" >> "$CSV"
  done
)

python3 "$HERE/build_assets.py"
python3 "$HERE/build_merigolix.py"   # faithful preview of the live Merigolix dashboard
echo "prepare-assets: done."
