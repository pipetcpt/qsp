# KBCS 2026 — *Advancing Drug Discovery and Development through Clinical Pharmacological Approaches on Organs-on-Chips*

Invited lecture for the **2026 Annual Conference of the Korean BioChip Society (KBCS)**,
session *Medical Unmet Needs for BioChips*. Speaker: **Sungpil Han, MD, PhD**
(clinical pharmacologist, The Catholic University of Korea / Seoul St. Mary's Hospital).

A single Quarto source (`slides.qmd`) renders to a self-contained **RevealJS** HTML
deck. Theme: clean academic light, **Inter** typeface (embedded, offline), native
MathML equations (no CDN dependency). A small QR code linking to
`github.com/pipetcpt/qsp` persists in the top-right corner of every slide.

## Structure (2026-07 revision)

The deck opens with the NAM translation loop (Cao, *Clin Pharmacol Ther* 2025)
to frame PBPK/QSP as the mechanistic layer that qualifies NAM readouts for
clinical use, then moves through an expanded clinical-pharmacology / QSP core
(adapted and substantially grown from `talks/ksmb2026`): what QSP/PBPK are,
the recurring math, MIDD, the Claude Code Routine engine, the 258-model open
QSP library, and many drug–disease case studies (a drugs/targets/endpoints
table, twelve one-insight vignettes, six full deep dives, and the Merigolix
live-app flagship case). The original Gut–Liver–Kidney MPS + PBPK case study
is condensed to ~1/3 its original length and moved to the back, keeping only
the essential mechanism and a single "predicted vs. observed" figure as
evidence the model fits well — everything else in that section is a table or
bullet list. There is no Acknowledgements slide and no "why a clinician"
slide. Every slide uses one uniform body font size — `{.smaller}` is reserved
for the handful of slides with a genuinely large data table (the 12-row
drugs/targets table, the 6-drug validation table, the references list); no
`stat-grid` boxes anywhere, only tables.

## Contents

| Path | What |
|------|------|
| `slides.qmd` | the deck (RevealJS) |
| `_variables.yml` | centralized live counts (QSP library size) — refresh anytime |
| `theme/custom.scss`, `theme/fonts.css`, `theme/fonts/` | shared clean-academic-light theme (same as `talks/ksmb2026`) |
| `theme/qr-overlay.html` | fixed-position QR code injected via `include-after-body` |
| `assets/poster/fig1_fit_only.png` | the ONE retained poster figure (cropped to just the predicted-vs-observed AUC/Cmax scatter panels) |
| `assets/figs/{igan,scd,mm,ra,dmd,ad}_map.png` | real Graphviz mechanistic maps for the six deep-dive diseases |
| `assets/figs/montage.png` | real mechanistic-map gallery montage |
| `assets/shiny/merigolix_dashboard.webp` | real Shiny dashboard screenshot (Merigolix flagship case) |
| `assets/figs/repo_qr.png` | QR code to the QSP library repo |
| `scripts/prepare-assets.sh` → `build_assets.py` | build the montage + deep-dive maps + QR (Python only: cairosvg, Pillow, qrcode) |
| `scripts/library-stats.sh` | recompute QSP-library counts from the repo → `_variables.yml` |
| `render.sh`, `Makefile` | one-command build |

## Content provenance

The condensed case-study section at the back is drawn from the original
standalone poster `poster_MPS2026_standalone.html` (Sungpil Han's Gut–Liver–
Kidney MPS + PBPK work) — mechanism, the falsifiable RAF scaling rule,
validation numbers, and disease-state results, all now presented as bullets
and tables rather than pasted figures, except the one retained fit-quality
scatter plot (cropped from the poster's own figure, not redrawn).

## Build

```bash
# tooling: quarto (no LaTeX needed for RevealJS);
#          pip install cairosvg pillow qrcode
make all            # assets + RevealJS -> _output/slides.html (self-contained)
```

## Updating the QSP-library numbers

```bash
bash scripts/library-stats.sh
make assets
```

## Publishing to GitHub Pages

The rendered deck is also published at `docs/kbcs2026/index.html` in the repo
root, served by GitHub Pages alongside `docs/index.html` (the KSMB 2026 talk).
After re-rendering, copy `_output/slides.html` to `../../docs/kbcs2026/index.html`
and commit.
