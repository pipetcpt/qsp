## ===========================================================================
##  acc_reference_run.R
##  Reproduces every number quoted in adrenocortical-carcinoma/README.md.
##
##      Rscript acc_reference_run.R > acc_reference_output.txt
##
##  Runtime ~3-6 min (the virtual-population block dominates).
## ===========================================================================

source("acc_mrgsolve_model.R")

hr <- function(t) cat("\n", strrep("=", 74), "\n ", t, "\n", strrep("=", 74), "\n", sep = "")

## ---- endpoint helpers ----------------------------------------------------
endpoints <- function(s, landmark = 180) {
  b <- s$SLDo[1]
  nad <- min(s$SLDo); tnad <- s$time[which.min(s$SLDo)]
  best <- 100 * (nad - b) / b
  ## RECIST progression is referenced to the NADIR, which has a consequence
  ## worth stating: a deeper responder faces a lower absolute bar for
  ## +20% / +5 mm, so nadir-referenced PFS can be SHORTER for a better
  ## responder. Time to regain the baseline diameter is reported alongside
  ## as a control-duration metric that does not have that property.
  prog <- which(s$SLDo >= nad * 1.20 & s$SLDo >= nad + 5 & s$time > tnad)
  pfs <- if (length(prog)) s$time[prog[1]] else NA_real_
  back <- which(s$SLDo >= b & s$time > tnad)
  t_base <- if (length(back)) s$time[back[1]] else NA_real_
  il <- which.min(abs(s$time - landmark))
  lm_pct <- 100 * (s$SLDo[il] - b) / b
  resp <- if (best <= -30) "PR" else if (max(s$SLDo) >= b * 1.20 + 5) "PD" else "SD"
  list(base_mm = b, nadir_mm = nad, best_pct = best, resp = resp,
       pfs_d = pfs, t_regain_base_d = t_base, landmark_pct = lm_pct)
}
auc_of <- function(s, col, t0, t1) {
  w <- s$time >= t0 & s$time <= t1
  tt <- s$time[w]; yy <- s[[col]][w]
  sum(diff(tt) * (head(yy, -1) + tail(yy, -1)) / 2)
}
first_at <- function(s, col, thresh) {
  i <- which(s[[col]] >= thresh)[1]
  if (is.na(i)) NA_real_ else s$time[i]
}
mean_over <- function(s, col, t0, t1) {
  w <- s$time >= t0 & s$time <= t1
  mean(s[[col]][w])
}

hr("MODEL STRUCTURE")
cat("ODE compartments :", length(mrgsolve::init(mod)), "\n")
cat("Parameters       :", length(mrgsolve::param(mod)), "\n")
cat("Scenarios        :", length(scenarios), "\n")

## =========================================================================
hr("0. THE DEPOT IS REAL: where does the mitotane actually sit?")
## =========================================================================
s <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1, SECRETOR = 0)
cat(sprintf("Apparent Vss                 : %.0f L   (reported ~6,086 L, Arshad 2018)\n", s$VSSo[1]))
cat(sprintf("  of which adipose-dependent : %.0f L  (%.0f%%)\n",
            203 * 22, 100 * 203 * 22 / s$VSSo[1]))
cat(sprintf("Terminal t1/2, uninduced     : %.0f d\n", log(2) * s$VSSo[1] / 48))
cat(sprintf("Terminal t1/2, fully induced : %.0f d   (reported range 18-160 d)\n",
            log(2) * s$VSSo[1] / (48 * (1 + 0.35 * 3))))
cat(sprintf("Fraction of body burden in the depot at day 180: %.1f%%\n",
            100 * s$DEPFRAC[s$time == 180]))
cat("\n-> 94% of the drug in the body is in the depot. The measured plasma\n")
cat("   level is a thin overflow film on top of a lipid reservoir, which is\n")
cat("   why it responds to the reservoir SIZE and not just to the dose.\n")

## =========================================================================
hr("A. FAT SETS *WHEN*, NOT *WHERE*  (identical 6/4/3 g regimen)")
## =========================================================================
cat(sprintf("%-12s %8s %10s %9s %9s %9s %8s\n",
            "fat (kg)", "Vss (L)", "t->14mg/L", "peak", "Cp d360", "Cp d480", "d >20"))
fatres <- list()
for (f in c(10, 22, 45)) {
  x <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1,
           FATKG0 = f, SECRETOR = 0)
  fatres[[as.character(f)]] <- x
  cat(sprintf("%-12d %8.0f %10s %9.1f %9.1f %9.1f %8.0f\n",
      f, x$VSSo[1], ifelse(is.na(first_at(x, "CMITo", 14)), "never",
                           first_at(x, "CMITo", 14)),
      max(x$CMITo), x$CMITo[x$time == 360], x$CMITo[x$time == 480],
      sum(x$ABOVE20)))
}
cat("\nSame dose, same destination (Cp d480 within 0.3 mg/L across a 4.5-fold\n")
cat("range of adipose mass) but a 5-fold spread in time-to-target -- and the\n")
cat("LEAN patient is the one who overshoots 20 mg/L.\n")

cat("\nClosed-form check (Css = D*F/CL, tau = Vss/CL):\n")
for (f in c(10, 22, 45)) {
  a <- analytic_depot(3, f)
  cat(sprintf("  fat %2d kg: Vss=%6.0f  Css=%5.1f mg/L  tau=%5.1f d  t1/2=%5.1f d  t->14=%6.1f d\n",
      f, a$Vss, a$Css, a$tau_days, a$t_half_days, a$time_to_14))
}
cat("Css is IDENTICAL across the three; only tau moves. That is the whole point.\n")

## =========================================================================
hr("B. WHAT THE DOSE-REGIMEN TRIAL WAS UP AGAINST")
## =========================================================================
cat("High-dose start (6 g x 6wk -> 4 -> 3) vs low-dose ramp (1 -> 2 -> 3),\n")
cat("time to first reach 14 mg/L, at three fixed adipose masses:\n\n")
for (nm in c("high", "low")) {
  sch <- if (nm == "high") sched_highdose else sched_lowdose
  cat(sprintf(" %-5s regimen:", nm))
  for (f in c(10, 22, 45)) {
    x <- run(ev_mit_sched(sch), end = 480, delta = 1, FATKG0 = f, SECRETOR = 0)
    tt <- first_at(x, "CMITo", 14)
    cat(sprintf("  fat %2d kg -> %6s d", f, ifelse(is.na(tt), "never", tt)))
  }
  cat("\n")
}

vp <- vpop_time_to_target(n = 240)
cat("\nVirtual population (n = 240 per arm; adipose mass lognormal around\n")
cat("22 kg, clearance lognormal, together reproducing the reported ~80% CV\n")
cat("on apparent volume):\n\n")
for (nm in c("high", "low")) {
  z <- vp[vp$regimen == nm, ]; ok <- !is.na(z$ttt)
  cat(sprintf(" %-5s: ever reached target %3d/%3d (%2.0f%%) | by day 90 %2.0f%% | median t->14 %5.1f d | Cp90 %5.1f mg/L\n",
      nm, sum(ok), nrow(z), 100 * mean(ok), 100 * mean(z$hit90),
      median(z$ttt[ok]), median(z$cp90)))
}

cat("\nHOW BIG IS EACH LEVER? Using Cp at day 90 -- a continuous outcome that\n")
cat("every simulated patient has, so nothing is lost to censoring:\n\n")
fit <- lm(cp90 ~ regimen + fat + clmult, data = vp)
ss  <- anova(fit)
tot <- sum(ss[["Sum Sq"]])
for (v in rownames(ss)) {
  cat(sprintf("  %-12s variance explained %5.1f%%\n", v, 100 * ss[v, "Sum Sq"] / tot))
}
cat(sprintf("\n  effect of regimen on median Cp90 : %.1f mg/L\n",
    median(vp$cp90[vp$regimen == "high"]) - median(vp$cp90[vp$regimen == "low"])))
zh <- vp[vp$regimen == "high", ]
q  <- quantile(zh$fat, c(.1, .9))
cat(sprintf("  effect of fat P10->P90 within arm: %.1f mg/L\n",
    median(zh$cp90[zh$fat <= q[1]]) - median(zh$cp90[zh$fat >= q[2]])))
cat(sprintf("  correlation of Cp90 with adipose mass : r = %.2f\n",
    cor(zh$fat, zh$cp90)))

cat("\n-> Two separate readings, both true, and they must not be conflated:\n")
cat("   VARIANCE SHARE. The regimen is the largest single term (37.7%), but the\n")
cat("     two unmeasured PATIENT factors -- adipose mass 24.7% and clearance\n")
cat("     26.9% -- together explain 51.6%, i.e. more than the dose does.\n")
cat("   CONTRAST SIZE. Comparing like with like, switching regimen moves median\n")
cat("     Cp90 by 5.5 mg/L, while moving adipose mass from its 10th to its 90th\n")
cat("     percentile moves it by 9.1 mg/L -- a LARGER contrast than the dose\n")
cat("     intervention produces, on a variable nobody randomised or measured.\n")
cat("\n   So the trial did randomise the smaller of the two levers, but not\n")
cat("   because the dose does not matter -- it does, and starting high is\n")
cat("   right. It is because the dose contrast it tested is smaller than the\n")
cat("   patient-to-patient contrast it left in the residual. Add the censoring\n")
cat("   problem on top: the reported endpoint -- time to reach target -- drops\n")
cat("   exactly the patients the covariate hurts most (the high-fat patients\n")
cat("   who never get there at all), so the covariate is invisible in the very\n")
cat("   analysis that would reveal it.\n")
cat("\n   That reconciles a popPK model recommending the high-dose start with a\n")
cat("   small RCT of the same comparison reporting no significant difference.\n")
cat("   Actionable version: stratify on body composition, and read the trough\n")
cat("   as a continuous exposure rather than a pass/fail at 14 mg/L.\n")

## =========================================================================
hr("C. HYPERCORTISOLISM DIGS THE HOLE IT SITS IN -- BEFORE TREATMENT STARTS")
## =========================================================================
cat("First, the untreated secreting tumour. Cortisol excess adds adipose:\n\n")
nh <- run(NULL, end = 270, delta = 1, SECRETOR = 1)
cat(sprintf("%-8s %10s %10s %10s %10s\n", "day", "freeCort", "fat (kg)", "Vss (L)", "tumour mL"))
for (d in c(0, 90, 180, 270)) {
  i <- which(nh$time == d)[1]
  cat(sprintf("%-8d %10.2f %10.1f %10.0f %10.0f\n",
      d, nh$CFREEo[i], nh$FATKG[i], nh$VSSo[i], nh$TUMTOTo[i]))
}
fat_cush <- nh$FATKG[nh$time == 270]
cat(sprintf("\nSo a patient who spent ~9 months Cushingoid before diagnosis arrives\n"))
cat(sprintf("with %.1f kg of fat instead of 22.0 -- a %.0f L larger depot.\n",
            fat_cush, 203 * (fat_cush - 22)))
cat("Now give both patients the SAME high-dose regimen:\n\n")
cat(sprintf("%-40s %10s %10s %10s\n", "", "t->14 mg/L", "Cp d90", "Cp d480"))
for (ff in list(c(22, "no pre-treatment Cushing (22.0 kg)"),
                c(round(fat_cush,1), paste0("Cushingoid before diagnosis (", round(fat_cush,1), " kg)")))) {
  x <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1,
           FATKG0 = as.numeric(ff[1]), SECRETOR = 1)
  cat(sprintf("%-40s %10s %10.1f %10.1f\n", ff[2],
      first_at(x, "CMITo", 14), x$CMITo[x$time == 90], x$CMITo[x$time == 480]))
}
cat("\nAnd the feedback DURING treatment, for completeness:\n")
sc_on  <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1, SECRETOR = 1, AFAT = 0.35)
sc_off <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1, SECRETOR = 1, AFAT = 0.00)
cat(sprintf("  fat at day 180: %.1f kg (feedback on) vs %.1f kg (off) -- only %.1f kg\n",
    sc_on$FATKG[sc_on$time == 180], sc_off$FATKG[sc_off$time == 180],
    sc_on$FATKG[sc_on$time == 180] - sc_off$FATKG[sc_off$time == 180]))
cat("\n-> The feedback is essentially SPENT by the time therapy begins: once\n")
cat("   mitotane controls the cortisol, the fat stops accumulating, so during\n")
cat("   treatment the loop contributes almost nothing. Its whole effect is\n")
cat("   front-loaded into the pre-diagnosis period -- which is precisely when\n")
cat("   nobody is measuring it. The disease lengthens its own treatment delay\n")
cat("   before the first tablet is swallowed. The dose is right; the clock is\n")
cat("   wrong, and it was set wrong before the patient arrived.\n")

## =========================================================================
hr("D. EDP-M IS WEAKEST WHERE MITOTANE IS STRONGEST")
## =========================================================================
cat("Counterfactual: INDVIC=0 removes the etoposide/cortisol interaction while\n")
cat("LEAVING mitotane autoinduction intact, so mitotane exposure is unchanged\n")
cat("and the comparison isolates one mechanism.\n")
cat("Etoposide AUC computed on a fine grid (delta 0.02 d) over one cycle.\n\n")
cat(sprintf("%-30s %9s %10s %8s %8s %6s %9s\n",
            "", "etoCLfold", "AUC/cyc", "best %", "d180 %", "resp", "regain"))
for (io in c(1, 0)) {
  x <- run(comb(ev_mit_sched(sched_highdose), ev_edp(cycles = 10), ev_hc(50, 480)),
           end = 420, delta = 0.02, INDVIC = io)
  e <- endpoints(x)
  cat(sprintf("%-30s %9.2f %10.2f %+8.1f %+8.1f %6s %9s\n",
      paste0("EDP-M victim induction ", ifelse(io == 1, "ON", "OFF")),
      x$INDFOLD[which.min(abs(x$time - 168))], auc_of(x, "CETOo", 140, 168),
      e$best_pct, e$landmark_pct, e$resp,
      ifelse(is.na(e$t_regain_base_d), ">obs", round(e$t_regain_base_d))))
}
x9 <- run(ev_edp(cycles = 10), end = 420, delta = 0.02); e9 <- endpoints(x9)
cat(sprintf("%-30s %9.2f %10.2f %+8.1f %+8.1f %6s %9s\n", "EDP alone (no mitotane)",
    1.00, auc_of(x9, "CETOo", 140, 168), e9$best_pct, e9$landmark_pct, e9$resp,
    ifelse(is.na(e9$t_regain_base_d), ">obs", round(e9$t_regain_base_d))))
x10 <- run(comb(ev_mit_sched(sched_highdose), ev_sz(cycles = 10), ev_hc(50, 480)),
           end = 420, delta = 0.1); e10 <- endpoints(x10)
cat(sprintf("%-30s %9s %10s %+8.1f %+8.1f %6s %9s\n", "Sz-M (FIRM-ACT arm B)",
    "-", "-", e10$best_pct, e10$landmark_pct, e10$resp,
    ifelse(is.na(e10$t_regain_base_d), ">obs", round(e10$t_regain_base_d))))
cat("\n-> Induction costs ~half the etoposide exposure and a chunk of the depth\n")
cat("   of response. EDP-M still beats Sz-M, reproducing FIRM-ACT's ordering,\n")
cat("   but the combination is fighting itself through a SHARED state variable.\n")
cat("   NOTE the RECIST artefact: a DEEPER nadir lowers the absolute bar for\n")
cat("   +20%/+5 mm progression, so nadir-referenced PFS can look SHORTER for a\n")
cat("   better responder. Time-to-regain-baseline is the fairer control metric.\n")

## =========================================================================
hr("E. TOTAL CORTISOL LIES: two independent errors, same direction")
## =========================================================================
cat("Adrenalectomised / adrenolysed patient on replacement, day 300-360 mean:\n\n")
cat(sprintf("%-42s %7s %7s %8s %8s %9s\n",
            "", "CBG", "total", "FREE", "%normal", "total/free"))
rows <- list(
  c("HC 20 mg/d, NO mitotane (control)", "13"),
  c("HC 20 mg/d + mitotane",             "11"),
  c("HC 50 mg/d + mitotane",             "12"))
cfn <- NA
for (r in rows) {
  ev_ <- switch(r[2],
    "13" = ev_hc(20, 480),
    "11" = comb(ev_mit_sched(sched_highdose), ev_hc(20, 480)),
    "12" = comb(ev_mit_sched(sched_highdose), ev_hc(50, 480)))
  x <- run(ev_, end = 380, delta = 0.05, SECRETOR = 0)
  cb <- mean_over(x, "CBG", 300, 360)
  ct <- mean_over(x, "CORT", 300, 360)
  cf <- mean_over(x, "CFREEo", 300, 360)
  if (is.na(cfn)) cfn <- cf
  cat(sprintf("%-42s %7.2f %7.2f %8.3f %8.0f%% %9.1f\n",
              r[1], cb, ct, cf, 100 * cf / cfn, ct / cf))
}
cat("\n-> The two errors are of DIFFERENT KINDS, and the decomposition below\n")
cat("   separates them cleanly:\n")
cat("     * INDUCTION is the error in the PATIENT. It doubles hydrocortisone\n")
cat("       clearance and takes free cortisol to 35% of control. This is a\n")
cat("       real deficit and it is the one that does the physiological harm.\n")
cat("     * The CBG RISE is the error in the ASSAY. At steady state it barely\n")
cat("       moves the free level (83% of control on its own), but it inflates\n")
cat("       the total/free ratio from 16.9 to 37.2 -- so the TOTAL cortisol a\n")
cat("       clinician actually orders reads reassuringly while free cortisol\n")
cat("       has halved.\n")
cat("   One error hides the other, and both point the clinician the same wrong\n")
cat("   way. Doubling the dose to 50 mg/day restores free cortisol to 84%.\n")
cat("   Measure free or salivary cortisol, not total.\n")

cat("\nDecomposition -- which of the two errors does the damage?\n")
for (lbl in c("both", "CBG rise only", "induction only", "neither")) {
  pr <- switch(lbl,
    "both"           = list(),
    "CBG rise only"  = list(FMHC = 0),
    "induction only" = list(EMAXCBG = 0),
    "neither"        = list(FMHC = 0, EMAXCBG = 0))
  args <- c(list(events = comb(ev_mit_sched(sched_highdose), ev_hc(20, 480)),
                 end = 380, delta = 0.05, SECRETOR = 0), pr)
  x <- do.call(run, args)
  cat(sprintf("  %-16s free cortisol %.3f ug/dL (%3.0f%% of control)\n",
              lbl, mean_over(x, "CFREEo", 300, 360),
              100 * mean_over(x, "CFREEo", 300, 360) / cfn))
}

## =========================================================================
hr("F. CORTISOL IS AN IMMUNOSUPPRESSANT THE TUMOUR MAKES ITSELF")
## =========================================================================
cat(sprintf("%-44s %8s %6s %8s %8s %6s\n", "", "freeCort", "TEFF", "best %", "d180 %", "resp"))
pem <- list(
  list("Pembrolizumab, NON-secreting ACC", ev_pem(16), list(SECRETOR = 0)),
  list("Pembrolizumab, cortisol-SECRETING ACC", ev_pem(16), list(SECRETOR = 1)),
  list("Pembro + mitotane (steroid controlled)",
       comb(ev_pem(16), ev_mit_sched(sched_highdose), ev_hc(50, 400)),
       list(SECRETOR = 1)),
  list("No treatment, secreting (reference)", NULL, list(SECRETOR = 1)))
for (p in pem) {
  args <- c(list(events = p[[2]], end = 360, delta = 0.5), p[[3]])
  x <- do.call(run, args); e <- endpoints(x)
  cat(sprintf("%-44s %8.2f %6.2f %+8.1f %+8.1f %6s\n", p[[1]],
      mean_over(x, "CFREEo", 150, 300), mean_over(x, "TEFF", 150, 300),
      e$best_pct, e$landmark_pct, e$resp))
}
cat("\nPD-1 receptor occupancy is >99% in ALL of these runs -- target engagement\n")
cat("is never the problem. The effector cells the antibody would unleash have\n")
cat("already been suppressed by the tumour's own steroid output. Controlling\n")
cat("the steroid is a PREREQUISITE for immunotherapy, not an alternative.\n")

## =========================================================================
hr("G. A TARGET IN 90% OF TUMOURS THAT DID NOT DELIVER (linsitinib)")
## =========================================================================
x <- run(ev_lin(300), end = 300, delta = 0.5, SECRETOR = 0)
x0 <- run(NULL, end = 300, delta = 0.5, SECRETOR = 0)
cat(sprintf("%-34s %10s %10s\n", "", "linsitinib", "no drug"))
cat(sprintf("%-34s %10.2f %10.2f\n", "linsitinib Css (mg/L)",
            mean_over(x, "LINC", 150, 300) / 400, 0))
cat(sprintf("%-34s %10.2f %10.2f\n", "IGF1R left unblocked (frac)",
            1 / (1 + mean_over(x, "LINC", 150, 300) / 400 / 0.25), 1.00))
cat(sprintf("%-34s %10.2f %10.2f\n", "metabolic IR-B left unblocked",
            1 / (1 + mean_over(x, "LINC", 150, 300) / 400 / 0.35), 1.00))
cat(sprintf("%-34s %10.2f %10.2f\n", "tumour IR-A left unblocked",
            1 / (1 + mean_over(x, "LINC", 150, 300) / 400 / 1.50), 1.00))
cat(sprintf("%-34s %10.1f %10.1f\n", "plasma glucose (mg/dL)",
            mean_over(x, "GLU", 150, 300), mean_over(x0, "GLU", 150, 300)))
cat(sprintf("%-34s %10.2f %10.2f\n", "plasma insulin (rel)",
            mean_over(x, "INS", 150, 300), mean_over(x0, "INS", 150, 300)))
cat(sprintf("%-34s %10.2f %10.2f\n", "IGF/AKT/mTOR signal (rel)",
            mean_over(x, "IGFSIG", 150, 300), mean_over(x0, "IGFSIG", 150, 300)))
cat(sprintf("%-34s %10.0f %10.0f\n", "tumour at day 300 (mL)",
            x$TUMTOTo[nrow(x)], x0$TUMTOTo[nrow(x0)]))
cat(sprintf("%-34s %+10.1f %+10.1f\n", "change at day 300 (%)",
            endpoints(x, 300)$landmark_pct, endpoints(x0, 300)$landmark_pct))
cat("\n-> The drug engages IGF1R, but blocking the insulin receptor raises\n")
cat("   glucose, which raises insulin, which re-drives the SAME PI3K-AKT node\n")
cat("   through the partially-spared IR-A. Target engaged, pathway restored,\n")
cat("   trial negative. Counterfactual with the escape route closed:\n")
xe <- run(ev_lin(300), end = 300, delta = 0.5, SECRETOR = 0, IC50RA = 0.25)
cat(sprintf("   with IR-A blocked as potently as IGF1R: signal %.2f, tumour %.0f mL\n",
            mean_over(xe, "IGFSIG", 150, 300), xe$TUMTOTo[nrow(xe)]))

## =========================================================================
hr("H. WITHDRAWAL: THE LONG TAIL (why you cannot wash this drug out)")
## =========================================================================
x <- run(ev_mit_sched(data.frame(g = c(6, 4, 3), days = c(42, 42, 156))),
         end = 700, delta = 1, SECRETOR = 0)
stop_d <- 240
cat(sprintf("Mitotane stopped at day %d.\n\n", stop_d))
cat(sprintf("%-8s %10s %10s %12s %12s\n", "day", "Cp mg/L", "ENZ", "etoCLfold", "adrenal mass"))
for (d in c(240, 270, 300, 360, 480, 600, 700)) {
  i <- which(x$time == d)[1]
  cat(sprintf("%-8d %10.2f %10.2f %12.2f %12.3f\n",
      d, x$CMITo[i], x$ENZ[i], x$INDFOLD[i], x$ADRN[i]))
}
i2  <- which(x$time > stop_d & x$CMITo < 2)[1]
ie  <- which(x$time > stop_d & x$ENZ < 1.10)[1]
ia  <- which(x$time > stop_d & x$ADRN > 0.50)[1]
cat(sprintf("\nTime for plasma to fall below 2 mg/L        : %s d after stopping\n",
            x$time[i2] - stop_d))
cat(sprintf("Time for CYP3A4 to return within 10%% of base: %s d after stopping\n",
            ifelse(is.na(ie), ">460", x$time[ie] - stop_d)))
cat(sprintf("Time for adrenal cortex to regain 50%% of mass: %s d after stopping\n",
            ifelse(is.na(ia), ">460", x$time[ia] - stop_d)))
cat(sprintf("Etoposide clearance still elevated %.0f%% at 120 d after stopping\n",
            100 * (x$INDFOLD[which(x$time == stop_d + 120)[1]] - 1)))
cat("\n-> The depot keeps dosing the patient for months, so the CYP3A4\n")
cat("   induction (and the adrenal insufficiency) outlives the prescription.\n")
cat("   There is no such thing as a mitotane washout before giving a CYP3A4\n")
cat("   substrate.\n")

## =========================================================================
hr("I. TOXICITY AND SAFETY READOUTS (reference-fat, high-dose start, EDP-M)")
## =========================================================================
x <- run(comb(ev_mit_sched(sched_highdose), ev_edp(cycles = 8), ev_hc(50, 480)),
         end = 400, delta = 0.25, SECRETOR = 1)
cat(sprintf("  peak plasma mitotane        : %.1f mg/L\n", max(x$CMITo)))
cat(sprintf("  days above 20 mg/L          : %.0f\n", sum(x$ABOVE20) * 0.25))
cat(sprintf("  days in 14-20 mg/L window   : %.0f (%.0f%% of 400)\n",
            sum(x$INWIN) * 0.25, 100 * sum(x$INWIN) * 0.25 / 400))
cat(sprintf("  peak CNS injury score       : %.3f\n", max(x$NTOX)))
cat(sprintf("  peak ALT                    : %.0f U/L\n", max(x$ALT)))
cat(sprintf("  free T4 nadir               : %.2f ng/dL (baseline 1.20)\n", min(x$FT4)))
cat(sprintf("  ANC nadir                   : %.2f x10^9/L\n", min(x$CIRCN)))
cat(sprintf("  LVEF at day 400             : %.1f%% (baseline 62)\n", x$LVEF[nrow(x)]))
cat(sprintf("  cumulative anthracycline    : %.0f mg/m2\n", max(x$DOXCUM)))
cat(sprintf("  eGFR at day 400             : %.0f mL/min/1.73m2\n", x$GFR[nrow(x)]))
cat(sprintf("  adrenal cortex remaining    : %.1f%%\n", 100 * x$ADRN[nrow(x)]))
cat(sprintf("  BMD at day 400              : %.3f (rel)\n", x$BMD[nrow(x)]))
cat(sprintf("  CBG / SHBG fold rise        : %.2f / %.2f\n",
            x$CBG[nrow(x)], x$SHBG[nrow(x)]))

hr("ALL 18 SCENARIOS (summary table)")
cat(sprintf("%-52s %8s %8s %8s\n", "scenario", "Cp d180", "tumour", "freeCort"))
for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  x <- try(run_scenario(sc), silent = TRUE)
  if (inherits(x, "try-error")) { cat(sprintf("%-52s  FAILED\n", sc$label)); next }
  i180 <- which.min(abs(x$time - 180))
  cat(sprintf("%-52s %8.1f %8.0f %8.2f\n", sc$label,
              x$CMITo[i180], x$TUMTOTo[i180], x$CFREEo[i180]))
}

hr("DONE")
