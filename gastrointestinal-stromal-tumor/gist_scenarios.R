## =============================================================================
##  GIST QSP model — scenario driver
## =============================================================================
##  Runs the 26 scenarios of gist_mrgsolve_model.R and computes the endpoints
##  quoted in README.md.  The dependency-free reference implementation of the
##  same 49-state system, with the same parameters and the same scenarios, is
##  gist_python_twin.py — it prints every number and asserts 41 of them, so run
##  that first if you want to check the model without R.
##
##    source("gist_scenarios.R")
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

mod <- mread_cache("gist", "gist_mrgsolve_model.R")

## -----------------------------------------------------------------------------
## dosing helpers
## -----------------------------------------------------------------------------
qd <- function(cmt, amt, start = 0, days = 2600) {
  if (amt <= 0 || days <= 0) return(NULL)
  ev(time = start, amt = amt, ii = 1, addl = days - 1, cmt = cmt)
}
## on/off cycling: sunitinib 4 weeks on 2 off, regorafenib 3 on 1 off
cyc <- function(cmt, amt, on, off, start = 0, days = 2600) {
  if (amt <= 0 || days <= 0) return(NULL)
  ncyc <- ceiling(days / (on + off))
  Reduce(`+`, lapply(seq_len(ncyc) - 1, function(k)
    ev(time = start + k * (on + off), amt = amt, ii = 1, addl = on - 1, cmt = cmt)))
}

## -----------------------------------------------------------------------------
## endpoints.  RECIST 1.1 progression is relative to the nadir observed on the
## CURRENT line, assessed every 7 days — the twin uses the same interval so the
## two implementations are comparable.  Real trials assess every 8-12 weeks,
## which quantises PFS more coarsely.
## -----------------------------------------------------------------------------
recist_pd <- function(d, t0 = 0, assess = 7) {
  d <- dplyr::filter(as.data.frame(d), time >= t0)
  if (nrow(d) < 3) return(NA_real_)
  tt  <- seq(min(d$time) + assess, max(d$time), by = assess)
  s   <- stats::approx(d$time, d$SLD, tt)$y
  nad <- cummin(s)
  hit <- which(s >= 1.2 * nad & s >= nad + 5)
  if (!length(hit)) return(NA_real_)
  tt[hit[1]] - t0
}
best_resp <- function(d, t0 = 0) {
  d <- dplyr::filter(as.data.frame(d), time >= t0)
  100 * (min(d$SLD) - d$SLD[1]) / d$SLD[1]
}
mo <- function(d) d / 30.4375

## -----------------------------------------------------------------------------
## Run a sequence of lines of therapy, switching at RECIST progression.
## `lines` is a list of lines; each line is a list of agents, and an empty line
## is best supportive care.
## -----------------------------------------------------------------------------
mk_ev <- function(agent, start, days) {
  cmt <- paste0("A_", agent$drug)
  switch(agent$sched,
         "qd"  = qd(cmt, agent$dose, start, days),
         "4/2" = cyc(cmt, agent$dose, 28, 14, start, days),
         "3/1" = cyc(cmt, agent$dose, 21, 7, start, days),
         stop("unknown schedule"))
}
run_sequence <- function(pars = list(), lines, end = 4000) {
  m <- do.call(param, c(list(mod), pars))
  t0 <- 0; evs <- NULL; res <- list()
  for (i in seq_along(lines)) {
    e <- NULL
    for (a in lines[[i]]) {
      e2 <- mk_ev(a, t0, end - t0)
      e  <- if (is.null(e)) e2 else e + e2
    }
    try_ev <- if (is.null(evs)) e else if (is.null(e)) evs else evs + e
    out <- mrgsim(m, if (is.null(try_ev)) ev(amt = 0, cmt = 1) else try_ev,
                  end = end, delta = 1)
    pd <- recist_pd(out, t0)
    res[[i]] <- list(t0 = t0, pfs = pd, best = best_resp(out, t0))
    if (is.na(pd)) { evs <- try_ev; break }
    e_trunc <- NULL
    for (a in lines[[i]]) {
      e2 <- mk_ev(a, t0, ceiling(pd))
      e_trunc <- if (is.null(e_trunc)) e2 else e_trunc + e2
    }
    evs <- if (is.null(evs)) e_trunc else if (is.null(e_trunc)) evs else evs + e_trunc
    t0  <- t0 + ceiling(pd)
  }
  sim <- mrgsim(m, if (is.null(evs)) ev(amt = 0, cmt = 1) else evs,
                end = end, delta = 1)
  list(sim = as.data.frame(sim), lines = res)
}
total_time <- function(r) mo(sum(vapply(r$lines,
  function(x) ifelse(is.na(x$pfs), 0, x$pfs), numeric(1))))

IM4  <- list(drug = "IM", dose = 400,  sched = "qd")
IM8  <- list(drug = "IM", dose = 800,  sched = "qd")
SU5  <- list(drug = "SU", dose = 50,   sched = "4/2")
SU37 <- list(drug = "SU", dose = 37.5, sched = "qd")
RE1  <- list(drug = "RE", dose = 160,  sched = "3/1")
RI1  <- list(drug = "RI", dose = 150,  sched = "qd")
RI0  <- list(drug = "RI", dose = 100,  sched = "qd")
AV3  <- list(drug = "AV", dose = 300,  sched = "qd")

## =============================================================================
## BLOCK 1 — natural history and the three time scales                     [C2]
## =============================================================================
s1 <- run_sequence(list(), list(list()),    end = 900)     # untreated
s2 <- run_sequence(list(), list(list(IM4)), end = 2600)    # KIT exon 11
with(s2$sim, {
  i48 <- which.min(abs(time - 2)); i56 <- which.min(abs(time - 56))
  cat(sprintf("SUV 48h %+.1f%%  Ki67 d7 %+.1f%%  SLD 48h %+.1f%%  SLD 8wk %+.1f%%\n",
              100 * (SUV[i48] / SUV[1] - 1),
              100 * (KI67[which.min(abs(time - 7))] / KI67[1] - 1),
              100 * (SLD[i48] / SLD[1] - 1), 100 * (SLD[i56] / SLD[1] - 1)))
})
## Three variables, three clocks: signalling in hours, cycling fraction in days,
## imaged mass in months.

## =============================================================================
## BLOCK 2 — first line by genotype, and the dose threshold                [C3]
## =============================================================================
## Exon 11 is flat above ~300 mg: escalation buys depth of response, not time.
## Exon 9 crosses the proliferation threshold between 400 and 800 mg, which is
## the MetaGIST result as a mechanism rather than as a subgroup finding.
grid <- expand.grid(GENO = 1:4, DOSE = c(200, 300, 400, 600, 800, 1200))
first_line <- lapply(seq_len(nrow(grid)), function(i)
  run_sequence(list(GENO = grid$GENO[i]),
               list(list(list(drug = "IM", dose = grid$DOSE[i], sched = "qd"))),
               end = 2600))
grid$pfs_mo  <- vapply(first_line, function(r) mo(r$lines[[1]]$pfs), numeric(1))
grid$best_pc <- vapply(first_line, function(r) r$lines[[1]]$best, numeric(1))
print(pivot_wider(grid, names_from = DOSE, values_from = pfs_mo, id_cols = GENO))

## PDGFRA D842V needs avapritinib, not more imatinib
d842_im <- run_sequence(list(GENO = 3), list(list(IM4)), end = 900)
d842_av <- run_sequence(list(GENO = 3), list(list(AV3)), end = 2600)

## =============================================================================
## BLOCK 3 — exposure: autoinduction, AGP, drug interaction, adherence     [C3]
## =============================================================================
low_exp <- run_sequence(list(CLF = 3.0, AGPF = 1.6), list(list(IM4)), end = 2600)
rifamp  <- run_sequence(list(INDF = 2.6),            list(list(IM4)), end = 2600)
## Prediction, and a MISS worth reporting: for exon 11 the model makes first-line
## PFS almost exposure-independent above a trough of ~250 ng/mL — only the DEPTH
## of response tracks exposure — because progression is driven by a clone
## imatinib never covers at any dose.  That predicts MetaGIST and under-predicts
## the Cmin/TTP association Demetri reported.

## =============================================================================
## BLOCK 4 — the INTRIGUE crossover                                        [C1]
## =============================================================================
## ctDNA subgroups are represented by which secondary-mutation class the tumour
## is able to generate, holding the TOTAL mutation flux constant.
mu_t    <- 2.0e-7 + 3.0e-7
atponly <- list(MU_ATP = mu_t,  MU_AL = 1e-11)
alonly  <- list(MU_ATP = 1e-11, MU_AL = mu_t)
crossover <- bind_rows(lapply(
  list(list(sub = "exon 11 + 13/14", p = atponly, d = SU5),
       list(sub = "exon 11 + 13/14", p = atponly, d = RI1),
       list(sub = "exon 11 + 17/18", p = alonly,  d = SU5),
       list(sub = "exon 11 + 17/18", p = alonly,  d = RI1)),
  function(x) {
    r <- run_sequence(x$p, list(list(IM4), list(x$d)), end = 3400)
    data.frame(subgroup = x$sub, drug = x$d$drug,
               pfs1 = mo(r$lines[[1]]$pfs), pfs2 = mo(r$lines[[2]]$pfs))
  }))
print(crossover)
## Observed (INTRIGUE ctDNA, Heinrich 2024 Nat Med):
##   exon 11+13/14  sunitinib 15.0  ripretinib  4.0 months
##   exon 11+17/18  sunitinib  1.5  ripretinib 14.2 months
##   ITT            sunitinib  8.3  ripretinib  8.0 months
## The two subgroups were mutually exclusive, so the ITT median is the median of
## a mixture of two opposite-signed effects — a property of the enrolled
## population, not of any patient.

## =============================================================================
## BLOCK 5 — the full sequence, and imatinib rechallenge                   [C1]
## =============================================================================
seq_full <- run_sequence(list(), list(list(IM4), list(SU5), list(RE1), list(RI1)))
seq_bsc  <- run_sequence(list(), list(list(IM4), list(SU5), list(RE1), list()))
seq_re   <- run_sequence(list(), list(list(IM4), list(SU5), list(RE1), list(IM4)))
## The fourth-line imatinib arm does anything at all because ~16% of viable cells
## still carry only the primary mutation (RIGHT: 1.8 vs 0.9 months).  Clone
## composition at that point:
with(seq_bsc$sim, {
  i <- which.min(abs(time - seq_bsc$lines[[4]]$t0))
  cat(sprintf("line 4 composition: primary %.1f%%  ATP %.1f%%  loop %.1f%%  bypass %.1f%%\n",
              100 * N_PRIM[i] / N_TOT[i], 100 * N_ATP[i] / N_TOT[i],
              100 * N_LOOP[i] / N_TOT[i], 100 * N_BYP[i] / N_TOT[i]))
})

## =============================================================================
## BLOCK 6 — ctDNA lead time
## =============================================================================
with(seq_bsc$sim, {
  t1 <- time[which(VAF_R >= 0.01)[1]]
  cat(sprintf("VAF 1%% on day %.0f; first-line PD on day %.0f; lead %.1f months\n",
              t1, seq_bsc$lines[[1]]$pfs, mo(seq_bsc$lines[[1]]$pfs - t1)))
})

## =============================================================================
## BLOCK 7 — stopping the drug: BFR14                                      [C2]
## =============================================================================
## Patients randomised after 1, 3 or 5 years were still non-progressing, which
## selects progressively harder.  The phenotype of each cohort is set from its
## OWN continuation arm, so the interruption arm is a prediction.
bfr_pheno <- list(`1` = list(TURNOVER = 85),
                  `3` = list(TURNOVER = 0),
                  `5` = list(TURNOVER = 0, KPMAX = 0.0225))
bfr <- bind_rows(lapply(names(bfr_pheno), function(y) {
  d <- as.numeric(y) * 365
  m <- do.call(param, c(list(mod), bfr_pheno[[y]]))
  s_stop <- mrgsim(m, qd("A_IM", 400, 0, d),        end = d + 1500, delta = 1)
  s_go   <- mrgsim(m, qd("A_IM", 400, 0, d + 3000), end = d + 3000, delta = 1)
  data.frame(years = y, interruption = mo(recist_pd(s_stop, d)),
             continuation = mo(recist_pd(s_go, d)))
}))
print(bfr)
## Observed (BFR14): 6.1 vs 27.8, 7.0 vs 67.0, 12.0 vs not reached months.

## =============================================================================
## BLOCK 8 — schedule and the toxicity/dose-intensity loop
## =============================================================================
su_42   <- run_sequence(list(), list(list(IM4), list(SU5)),  end = 3000)
su_cont <- run_sequence(list(), list(list(IM4), list(SU37)), end = 3000)
## DIFIX = 1 turns the dose-reduction feedback off, which separates how much of
## regorafenib's short third-line PFS is coverage and how much is tolerability.
reg_fixed <- run_sequence(list(DIFIX = 1), list(list(IM4), list(SU5), list(RE1)))

## =============================================================================
## BLOCK 9 — adjuvant therapy over a virtual population
## =============================================================================
## RFS at 5 and 10 years for 12 / 36 / 60 months of adjuvant imatinib.  Observed
## (SSGXVIII, Joensuu 2020): 5-year RFS 53.0% (12 mo) vs 71.4% (36 mo); 10-year
## 41.8% vs 52.5%.  The model has no cure term, so the curves separate while the
## drug is given and converge afterwards.
set.seed(20260730)
n   <- 400
pop <- data.frame(ID = seq_len(n),
                  N_MICRO  = 10^rnorm(n, log10(2.2e6), 0.55),
                  TURNOVER = pmax(0, rnorm(n, 85, 45)),
                  KPMAX    = rlnorm(n, log(0.029), 0.25),
                  CLF      = rlnorm(n, 0, 0.30))
adj <- bind_rows(lapply(c(365, 1095, 1826), function(dur) {
  out <- as.data.frame(mrgsim(param(mod, SETTING = 2), idata = pop,
                              events = qd("A_IM", 400, 0, dur),
                              end = 3652, delta = 28))
  out |> group_by(time) |>
    summarise(RFS = mean(DETECT == 0), .groups = "drop") |>
    mutate(duration = sprintf("%.0f months", dur / 30.4375))
}))
ggplot(adj, aes(time / 365, 100 * RFS, colour = duration)) +
  geom_step(linewidth = 0.8) +
  geom_vline(xintercept = c(1, 3, 5), linetype = 3, colour = "grey60") +
  labs(x = "years since resection", y = "recurrence-free (%)", colour = NULL) +
  theme_minimal()

## =============================================================================
## BLOCK 10 — cover the clones early?  The model says no, and says why     [C1]
## =============================================================================
strat_seq  <- run_sequence(list(), list(list(IM4), list(RI1), list(SU5)))
strat_comb <- run_sequence(list(), list(list(IM4, RI0), list(SU5)))
cat(sprintf("sequential total %.1f mo   upfront combination total %.1f mo\n",
            total_time(strat_seq), total_time(strat_comb)))
## The combination buys first-line PFS and loses total time, because covering the
## activation-loop clone before it is selected removes the heterogeneity the
## later lines exploit.  Not the intuitive answer, and a falsifiable statement
## about any future combination trial.

## =============================================================================
## BLOCK 11 — falsification: one clone, resistance as potency drift
## =============================================================================
## The fair single-compartment competitor is not "no resistance" but resistance
## as a continuous loss of sensitivity under drug pressure, REFITTED to the same
## first-line PFS anchor.  In gist_python_twin.py this is FALSIFY = 1 with
## KDRIFT refitted by bisection to 2.1e-3/day.  It then loses two things
## decisively — the INTRIGUE crossover (with one clone the subgroup parameters
## are inert, so the crossover is unreachable at ANY parameter value) and
## fourth-line ripretinib — while the imatinib-rechallenge benefit does NOT
## discriminate between the two structures.  See README.md section 9.
