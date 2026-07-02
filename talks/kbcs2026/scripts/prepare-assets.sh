#!/usr/bin/env bash
# prepare-assets.sh — refresh live counts + regenerate the generated figures
# (montage, bridge diagram, webapp mock, QR). The 5 poster figures under
# assets/poster/ are committed source images — this script does NOT touch them.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

bash "$HERE/library-stats.sh" > "$HERE/../_variables.yml" 2>/dev/null || true
cat >> "$HERE/../_variables.yml" <<'YAML'
n_categories: 18
talk_date: "2026"
YAML

python3 "$HERE/build_assets.py"
echo "prepare-assets: done."
