# KSMB 2026 talk — build targets
# Requires: quarto (PATH); python3 with cairosvg, pillow, matplotlib, numpy.
# Beamer additionally needs a LaTeX engine: `quarto install tinytex`.

QUARTO ?= quarto

.PHONY: all assets revealjs beamer clean

all: assets revealjs beamer

assets:
	bash scripts/prepare-assets.sh

revealjs: assets
	$(QUARTO) render slides.qmd --to revealjs

beamer: assets
	$(QUARTO) render slides.qmd --to beamer

clean:
	rm -rf _output .quarto slides_files
