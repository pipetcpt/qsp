# KBCS 2026 talk — build targets (RevealJS HTML)
QUARTO ?= quarto

.PHONY: all assets revealjs clean

all: revealjs

assets:
	bash scripts/prepare-assets.sh

revealjs: assets
	$(QUARTO) render slides.qmd --to revealjs

clean:
	rm -rf _output .quarto slides_files
