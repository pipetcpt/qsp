# KSMB 2026 — *An LLM-Augmented R/Shiny Platform for QSP*

Conference talk (~25 min) for the **2026 Annual Conference of the Korean Society
for Mathematical Biology**, Special Session B1 — *Mathematical Modeling in New
Drug Development and Clinical Pharmacotherapy*. Speaker: **Sungpil Han, MD, PhD**.

A single Quarto source (`slides.qmd`) renders to a self-contained **RevealJS**
HTML deck. Theme: clean academic light, **Inter** typeface (embedded, offline).
Figures are built from the repo's mechanistic-map SVGs plus generated diagrams,
plots, and a repository QR code.

## Contents

| Path | What |
|------|------|
| `slides.qmd` | the deck (RevealJS) |
| `_variables.yml` | centralized live counts (models, refs, clusters…) — refresh anytime |
| `_quarto.yml` | project config (output dir) |
| `theme/custom.scss` | clean-academic-light theme |
| `theme/fonts/` | embedded Inter `.woff2` files |
| `scripts/library-stats.sh` | recompute counts from the repo → `_variables.yml` |
| `scripts/prepare-assets.sh` → `build_assets.py` | build all figures (Python only) |
| `scripts/make-figs.R` | *optional* real-`mrgsolve` simulation figures |
| `assets/` | generated `figs/` (maps, montage, QR), `plots/` (diagrams, curves), plus `shiny/` & `logos/` placeholders |
| `render.sh`, `Makefile` | one-command build |

## Build

```bash
# tooling
#   quarto:  https://quarto.org/docs/get-started/   (RevealJS needs no LaTeX)
#   figures: pip install cairosvg pillow matplotlib numpy "qrcode[pil]"

make all            # assets + RevealJS  ->  _output/slides.html (self-contained)
# or:
make assets         # regenerate figures + counts only
make revealjs       # render the HTML only
```

Output: `_output/slides.html` (open in any browser; fonts and images embedded).

## Updating the numbers (the library grows daily)

```bash
bash scripts/library-stats.sh            # prints live YAML
make assets                              # refresh _variables.yml + figures
```
The deck reads counts via Quarto `{{< var ... >}}` shortcodes, so a single
re-run updates every slide.

## Assets to drop in (optional — replace placeholders, then re-render)

- `assets/logos/logos_row.png` — institutional logos (Catholic Univ. / PIPET)
- `assets/shiny/igan_dashboard.png` — a real Shiny dashboard screenshot

## Note on the simulation plots

`build_assets.py` draws **illustrative reduced-model** dynamics (Euler
integration of a simplified subsystem) so the deck is self-contained without R.
For figures from the **actual** `mrgsolve` models, use `scripts/make-figs.R`.
These models are educational/research QSP models — **not** for clinical use.
