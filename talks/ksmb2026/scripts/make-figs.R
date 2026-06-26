#!/usr/bin/env Rscript
# make-figs.R — OPTIONAL: regenerate the three deep-dive simulation figures from
# the REAL mrgsolve models (instead of the illustrative reduced-model curves that
# build_assets.py draws). Requires R with: mrgsolve, dplyr, ggplot2.
#
#   Rscript make-figs.R
#
# This is provided for users who have the full R toolchain. The slide deck renders
# fine without it — build_assets.py already produced labeled illustrative plots.

suppressMessages({
  ok <- requireNamespace("mrgsolve", quietly = TRUE) &&
        requireNamespace("ggplot2",  quietly = TRUE)
})
if (!ok) stop("Install mrgsolve + ggplot2 to use make-figs.R")

library(mrgsolve); library(ggplot2)
root  <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..", ".."))
plots <- file.path(dirname(sys.frame(1)$ofile), "..", "assets", "plots")
dir.create(plots, showWarnings = FALSE, recursive = TRUE)

run_one <- function(model_path, end = 365, ...) {
  mod <- mread(model_path)
  as.data.frame(mrgsim(mod, end = end, ...))
}

# Example: IgA nephropathy UPCR / eGFR trajectories.
# (Scenario columns/parameters depend on the model's own definitions; adapt the
#  `param()`/`ev()` calls to the scenario encoding in igan_mrgsolve_model.R.)
igan <- file.path(root, "iga-nephropathy", "igan_mrgsolve_model.R")
if (file.exists(igan)) {
  out <- run_one(igan, end = 720)
  # ... select UPCR/eGFR columns and ggsave to plots/plot_igan.png ...
  message("igan simulated: ", paste(head(names(out)), collapse = ", "))
}
message("make-figs.R: adapt scenario calls per model, then ggsave to assets/plots/.")
