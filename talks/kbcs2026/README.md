# KBCS 2026 — *Advancing Drug Discovery and Development through Clinical Pharmacological Approaches on Organs-on-Chips*

Invited lecture for the **2026 Annual Conference of the Korean BioChip Society (KBCS)**,
session *Medical Unmet Needs for BioChips*. Speaker: **Sungpil Han, MD, PhD**
(clinical pharmacologist, The Catholic University of Korea / Seoul St. Mary's Hospital).

A single Quarto source (`slides.qmd`) renders to a self-contained **RevealJS** HTML
deck. Theme: clean academic light, **Inter** typeface (embedded, offline), native
MathML equations (no CDN dependency). A small QR code linking to
`github.com/pipetcpt/qsp` persists in the bottom-left corner of every slide.

## Contents

| Path | What |
|------|------|
| `slides.qmd` | the deck (RevealJS) |
| `_variables.yml` | centralized live counts (QSP library size) — refresh anytime |
| `theme/custom.scss`, `theme/fonts.css`, `theme/fonts/` | shared clean-academic-light theme (same as `talks/ksmb2026`) |
| `theme/qr-overlay.html` | fixed-position QR code injected via `include-after-body` |
| `assets/poster/` | the 5 figures from the original MPS Gut–Liver–Kidney poster (source images, committed as-is, downscaled for web) |
| `assets/figs/` | generated figures: QSP-library montage, chip→PBPK→QSP bridge diagram, future-webapp mockup, repo QR |
| `scripts/prepare-assets.sh` → `build_assets.py` | build all generated figures (Python only: cairosvg, Pillow, matplotlib, numpy, qrcode) |
| `scripts/library-stats.sh` | recompute QSP-library counts from the repo → `_variables.yml` |
| `render.sh`, `Makefile` | one-command build |

## Content provenance

The case-study section (platform, methods ①–⑤, Results 1–2, limitations table,
conclusions, references) is merged in full from the original standalone poster
`poster_MPS2026_standalone.html` (Sungpil Han's Gut–Liver–Kidney MPS + PBPK work).
Figures are the poster's own rendered PNGs, downscaled for a lighter self-contained
HTML file — not redrawn.

## Build

```bash
# tooling: quarto (no LaTeX needed for RevealJS);
#          pip install cairosvg pillow matplotlib numpy qrcode
make all            # assets + RevealJS -> _output/slides.html (self-contained)
```

## Updating the QSP-library numbers

```bash
bash scripts/library-stats.sh
make assets
```
