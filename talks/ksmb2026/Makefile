# KSMB 2026 talk — build targets (RevealJS HTML)
# Requires: quarto (PATH); python3 with cairosvg, pillow, matplotlib, numpy, qrcode.

QUARTO ?= quarto

.PHONY: all assets revealjs publish clean

all: revealjs

assets:
	bash scripts/prepare-assets.sh

revealjs: assets
	$(QUARTO) render slides.qmd --to revealjs

# Copy the rendered deck to the repo's GitHub Pages folder (served at
# https://pipetcpt.github.io/qsp/). Run after `make revealjs`, then commit docs/.
publish: revealjs
	cp _output/slides.html ../../docs/index.html
	@echo "Published to docs/index.html — commit & push, Pages serves it at /qsp/"

clean:
	rm -rf _output .quarto slides_files
