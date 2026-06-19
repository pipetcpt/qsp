## ============================================================
## validate.R — reproducibility & validation checks for the NAFLD/NASH QSP model
## Run from repo root:   Rscript validation/validate.R
## Requires: R (4.5/4.6) + Rtools + {mrgsolve, dplyr, ggplot2, tidyr, patchwork}
## Exits 0 if all checks pass, 1 otherwise.
## ============================================================
suppressWarnings({
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  ROOT <- if (length(f)) normalizePath(file.path(dirname(f), "..")) else normalizePath(".")
})
MODEL <- file.path(ROOT, "nafld_mrgsolve_model.R")
FIG   <- file.path(ROOT, "validation", "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

owd <- setwd(FIG)                          # model auto-saves scenario PNGs here
cat("Sourcing model:", MODEL, "\n\n")
suppressWarnings(suppressMessages(source(MODEL, chdir = FALSE, echo = FALSE)))
setwd(owd)

g <- function(df, sc, col, wk) df[[col]][df$Scenario == sc & abs(df$Week - wk) < 0.1][1]
pass <- TRUE
chk  <- function(name, ok, detail) {
  cat(sprintf("[%s] %-32s %s\n", ifelse(ok, "PASS", "FAIL"), name, detail)); if (!ok) pass <<- FALSE
}

cat("── Validation checks (model vs published, placebo-corrected) ──\n")
## 1) placebo is an exact flat steady state
pl <- results[results$Scenario == "Placebo", ]
drift <- max(abs(pl$LF_PCT - pl$LF_PCT[which.min(pl$time)]))
chk("placebo steady-state (flat)", drift < 1e-2, sprintf("max LF drift = %.3g", drift))

## 2) no divergence; states finite & physiological
mx <- max(results$LF_PCT, na.rm = TRUE); fin <- all(is.finite(results$LF_PCT))
chk("no divergence / finite", fin && mx < 60, sprintf("max LF = %.1f%%, finite = %s", mx, fin))

## 3) resmetirom 100 mg liver fat ≈ MAESTRO-NAFLD-1 wk52 -33.9%
rb <- g(results, "Resmetirom 100 mg QD", "PDFF", 0); r52 <- g(results, "Resmetirom 100 mg QD", "PDFF", 52)
rsm <- 100 * (r52 - rb) / rb
chk("resmetirom PDFF @wk52", abs(rsm - (-33.9)) < 4, sprintf("%.1f%%  (MAESTRO-NAFLD-1 -33.9%%)", rsm))

## 4) empagliflozin liver fat ≈ E-LIFT placebo-corrected -24.7%
ev_emp <- build_regimen(dose_emp = 10, duration_wk = 72)
emp <- as.data.frame(mrgsim(param(mod, DOSE_EMP = 1), ev = ev_emp, tgrid = tgrid, obsonly = TRUE))
emp$Week <- emp$time / 168
empa <- 100 * (emp$PDFF[abs(emp$Week - 72) < 0.1][1] / emp$PDFF[abs(emp$Week - 0) < 0.1][1] - 1)
chk("empagliflozin liver fat @wk72", abs(empa - (-24.7)) < 4, sprintf("%.1f%%  (E-LIFT placebo-corr -24.7%%)", empa))

## 5) semaglutide 2.4 mg weight ≈ ESSENCE -10.5% / STEP -14.9%
sb <- g(results, "Semaglutide 2.4 mg QW", "BODY_WT", 0); s72 <- g(results, "Semaglutide 2.4 mg QW", "BODY_WT", 72)
sema <- 100 * (s72 - sb) / sb
chk("semaglutide weight @wk72", sema <= -9 && sema >= -16, sprintf("%.1f%%  (ESSENCE -10.5%% / STEP -14.9%%)", sema))

cat("\n", ifelse(pass, "==> ALL CHECKS PASSED", "==> SOME CHECKS FAILED"), "\n", sep = "")
cat("Scenario figures written to: ", FIG, "\n", sep = "")
quit(status = ifelse(pass, 0L, 1L))
