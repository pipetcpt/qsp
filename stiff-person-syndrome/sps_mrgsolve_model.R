## =============================================================
## Stiff Person Syndrome (SPS) — Quantitative Systems Pharmacology
## mrgsolve ODE model (26 compartments)
##
## Disease: anti-GAD65 autoimmunity → ↓ GAD65 enzyme → ↓ CNS GABA
##          → α-MN hyper-excitability → axial rigidity + spasms
##
## Drugs implemented (PK + PD coupling):
##   • Diazepam (oral) + active metabolite desmethyldiazepam
##   • Baclofen (oral or intrathecal pump)
##   • Gabapentin (oral, α2δ binding)
##   • IVIG (IV, anti-idiotype + FcRn-mediated catabolism)
##   • Rituximab (IV, CD20-mediated B-cell depletion)
##   • Plasmapheresis (PLEX, instantaneous removal events)
##   • Prednisolone (oral)
##
## Scenarios in `sps_run()`:
##   1. Newly-diagnosed, BZD only
##   2. BZD + baclofen oral combo
##   3. BZD + IVIG q4w cycle
##   4. BZD + rituximab induction (375 mg/m2 x4)
##   5. PLEX rescue (refractory crisis)
##   6. Intrathecal baclofen pump
##
## All parameter values are illustrative, anchored to published
## ranges; see `sps_references.md`.
## =============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

sps_code <- '
$PROB
# Stiff Person Syndrome QSP (SPS-QSP v1.0)
# Reference parameter values approximate adult, 70 kg, GAD-Ab+ SPS

$PARAM @annotated
// ---- Diazepam PK (2-cmt + active metabolite) ----
KA_DIAZ   :  1.2     : Diazepam absorption rate (1/h)
CL_DIAZ   :  1.6     : Diazepam clearance (L/h)
V1_DIAZ   : 70       : Diazepam Vc (L)
Q_DIAZ    :  3.5     : Diazepam Q (L/h)
V2_DIAZ   : 120      : Diazepam Vp (L)
FM_DMD    :  0.55    : Fraction metabolized to desmethyldiazepam
CL_DMD    :  0.45    : DMD clearance (L/h)
V_DMD     : 90       : DMD volume (L)
FBP_DIAZ  :  0.015   : Diazepam free fraction (1.5%)
FBP_DMD   :  0.020   : DMD free fraction
BBB_DIAZ  :  0.90    : Diazepam brain partition

// ---- Baclofen PK (oral + intrathecal) ----
KA_BAC    :  0.8     : Baclofen oral Ka (1/h)
F_BAC     :  0.80    : Baclofen oral bioavailability
CL_BAC    : 11.5     : Baclofen CL (L/h, mostly renal)
V_BAC     : 50       : Baclofen V (L)
CSF_KIN   :  0.025   : Plasma->CSF rate (1/h)
CSF_KOUT  :  0.06    : CSF clearance rate (1/h)
V_CSF     :  0.140   : CSF volume (L)
IT_F      :  1.0     : Intrathecal bioavailability

// ---- Gabapentin PK ----
KA_GAB    :  0.7     : Gabapentin Ka
F_GAB     :  0.45    : Gabapentin saturable F (lumped at clinical doses)
CL_GAB    :  9       : Gabapentin CL (L/h)
V_GAB     :  65      : Gabapentin V (L)

// ---- IVIG / IgG PK (anti-idiotype) ----
CL_IGG    :  0.32    : Total IgG clearance (L/d)
V_IGG     :  3.5     : IgG distribution volume (L)
KEL_IGG   :  0.0046  : IgG kel (1/h) ~ t1/2 ~21d

// ---- Rituximab PK (TMDD-light, 2-cmt) ----
CL_RTX    :  0.27    : Rituximab CL (L/d)
V1_RTX    :  3.0     : Rituximab Vc (L)
V2_RTX    :  2.7     : Rituximab Vp (L)
Q_RTX     :  0.40    : Rituximab Q (L/d)

// ---- Prednisolone PK ----
KA_PRED   :  1.5
CL_PRED   :  17.3
V_PRED    :  35

// ---- B-cell / plasma-cell dynamics ----
B0        : 250      : Baseline CD20+ B cells (cells/uL)
P0        :  15      : Baseline plasmablasts (cells/uL)
LL0       :  10      : Long-lived plasma cells (relative units)
KPROD_B   :  2.5     : B-cell production (cells/uL/d)
KDEG_B    :  0.010   : B-cell turnover (1/d)
KRTX_KILL :  0.45    : Rituximab Emax kill (1/d at saturation)
KC50_RTX  :  0.5     : Rituximab Cp at half-Emax kill (ug/mL)
KMAT_PB   :  0.05    : Naive B -> plasmablast (1/d)
KMAT_LL   :  0.008   : Plasmablast -> LLPC (1/d)
KDEG_PB   :  0.07    : Plasmablast death (1/d)
KDEG_LL   :  0.002   : LLPC death (1/d) (long-lived)

// ---- Anti-GAD65 antibody dynamics ----
ABAB0     : 5000     : Baseline anti-GAD65 (relative U/mL)
KPROD_AB_PB : 30     : Plasmablast Ab production rate (U/mL per cell/uL per d)
KPROD_AB_LL : 12     : LLPC Ab production rate
KDEG_AB   :  0.04    : Ab elimination (1/d, ~17 d half-life)
IVIG_NEUT :  0.0010  : IVIG anti-idiotype neutralization potency
KCSF_AB   :  0.001   : Serum->CSF Ab transfer (1/d)
KCSF_AB_OUT : 0.05   : CSF Ab clearance (1/d)

// ---- GAD65 enzyme / GABA pool ----
GAD0      : 100      : GAD65 activity (% baseline) target
KIN_GAD   :  0.5     : GAD synthesis (1/d)
KOUT_GAD  :  0.005   : GAD turnover (1/d) at baseline = KIN_GAD/GAD0
AB_KILL   :  0.00006 : Ab inhibitory potency on GAD (per U/mL per d)
KIN_GABA  :  1.0     : GABA synthesis const (1/d) baseline 100
KOUT_GABA :  0.01    : GABA degradation (1/d)
GABA0     : 100      : Steady-state GABA (%)

// ---- Excitability cascade ----
MN0       : 100      : alpha-MN excitability baseline (%)
KIN_MN    :  4.0     : MN synthesis (1/d)
KOUT_MN   :  0.04    : MN excitability decay (1/d)
EMAX_GABA :  0.85    : Max GABA-mediated MN suppression
EC50_GABA : 60       : GABA at half max inhibition (% baseline)
EMAX_BZD  :  0.40    : BZD allosteric Emax
EC50_BZD  :  0.20    : Free brain diazepam-equiv (mg/L) for half Emax
EMAX_BAC  :  0.35    : Baclofen Emax
EC50_BAC  :  0.08    : Baclofen CSF EC50 (mg/L)
EMAX_GAB  :  0.20    : Gabapentin Emax (α2δ)
EC50_GAB  :   4.0    : Gabapentin plasma EC50 (mg/L)

// ---- Symptom transducers ----
STIFF0    : 0        : Stiffness score baseline (HSI scale 0-100)
SPASM0    : 0        : Spasm frequency baseline (events/d)
KIN_STIFF :  6       : Stiffness build-up (1/d) when MN > threshold
KOUT_STIFF:  0.10    : Stiffness resolution (1/d)
KIN_SPASM :  4       : Spasm forcing (1/d)
KOUT_SPASM:  0.30    : Spasm resolution (1/d)
MN_THR    : 100      : MN threshold for symptom generation

// ---- Steroid effect ----
EMAX_PRED :  0.30    : Prednisolone effect on B-cell production
EC50_PRED :  0.05    : Prednisolone (mg/L) for half effect

// ---- Bone / safety surrogate ----
BMD0      : 1.0      : Baseline BMD (g/cm^2)
KBMD_PRED :  0.00008 : Steroid-driven BMD loss per ng/mL/d

// ---- Body weight (covariate) ----
WT        : 70       : Body weight (kg)

$CMT @annotated
DIAZ_GUT  : Diazepam gut (mg)
DIAZ_C    : Diazepam central (mg)
DIAZ_P    : Diazepam peripheral (mg)
DMD_C     : Desmethyldiazepam (mg)
BAC_GUT   : Baclofen gut (mg)
BAC_C     : Baclofen central (mg)
BAC_CSF   : Baclofen CSF (mg)
GAB_GUT   : Gabapentin gut (mg)
GAB_C     : Gabapentin central (mg)
PRED_GUT  : Prednisolone gut (mg)
PRED_C    : Prednisolone central (mg)
IVIG_C    : IVIG plasma (g)
RTX_C     : Rituximab central (mg)
RTX_P     : Rituximab peripheral (mg)
BNAIVE    : Naive B cells (cells/uL)
PBLAST    : Plasmablasts (cells/uL)
LLPC      : Long-lived plasma cells (rel U)
ABAB_S    : Anti-GAD65 serum (U/mL)
ABAB_CSF  : Anti-GAD65 CSF (U/mL)
GAD       : GAD65 enzyme activity (%)
GABA      : CNS GABA pool (%)
MN_EXC    : alpha-MN excitability (%)
STIFF     : Stiffness score (0-100)
SPASM     : Spasm freq (events/d)
BMD       : Lumbar BMD (g/cm^2)
HSI       : Heightened Sensitivity Index

$GLOBAL
#define DIAZ_FREE_BRAIN ( (DIAZ_C/V1_DIAZ) * FBP_DIAZ * BBB_DIAZ + (DMD_C/V_DMD) * FBP_DMD * BBB_DIAZ )
#define BAC_CSF_C       ( BAC_CSF / V_CSF )
#define GAB_C_C         ( GAB_C / V_GAB )
#define PRED_C_C        ( PRED_C / V_PRED )
#define RTX_C_C         ( RTX_C / V1_RTX )

$MAIN
F_BAC_GUT = F_BAC;
F_GAB_GUT = F_GAB;
ABAB_S_0   = ABAB0;
ABAB_CSF_0 = 0.02 * ABAB0;
GAD_0      = GAD0;
GABA_0     = GABA0;
MN_EXC_0   = MN0;
STIFF_0    = STIFF0;
SPASM_0    = SPASM0;
BMD_0      = BMD0;
BNAIVE_0   = B0;
PBLAST_0   = P0;
LLPC_0     = LL0;
HSI_0      = 0;

$ODE
// ---- Diazepam PK + active metabolite ----
dxdt_DIAZ_GUT = -KA_DIAZ * DIAZ_GUT;
dxdt_DIAZ_C   =  KA_DIAZ * DIAZ_GUT
                - (CL_DIAZ/V1_DIAZ)*DIAZ_C
                - (Q_DIAZ /V1_DIAZ)*DIAZ_C
                + (Q_DIAZ /V2_DIAZ)*DIAZ_P;
dxdt_DIAZ_P   =  (Q_DIAZ /V1_DIAZ)*DIAZ_C - (Q_DIAZ /V2_DIAZ)*DIAZ_P;
dxdt_DMD_C    =  FM_DMD * (CL_DIAZ/V1_DIAZ)*DIAZ_C - (CL_DMD/V_DMD)*DMD_C;

// ---- Baclofen PK (oral + IT) ----
dxdt_BAC_GUT = -KA_BAC * BAC_GUT;
dxdt_BAC_C   =  KA_BAC * BAC_GUT - (CL_BAC/V_BAC)*BAC_C - CSF_KIN*BAC_C + CSF_KOUT*BAC_CSF*(V_CSF/V_BAC);
dxdt_BAC_CSF =  CSF_KIN*BAC_C - CSF_KOUT*BAC_CSF;

// ---- Gabapentin PK ----
dxdt_GAB_GUT = -KA_GAB * GAB_GUT;
dxdt_GAB_C   =  KA_GAB * GAB_GUT - (CL_GAB/V_GAB)*GAB_C;

// ---- Prednisolone ----
dxdt_PRED_GUT = -KA_PRED * PRED_GUT;
dxdt_PRED_C   =  KA_PRED * PRED_GUT - (CL_PRED/V_PRED)*PRED_C;

// ---- IVIG PK (g; treat as 1-cmt) ----
dxdt_IVIG_C = -KEL_IGG * IVIG_C * 24;     // convert per d -> per h scale

// ---- Rituximab PK ----
dxdt_RTX_C  =  (Q_RTX/V2_RTX)*RTX_P - (Q_RTX/V1_RTX + CL_RTX/V1_RTX)*RTX_C;
dxdt_RTX_P  =  (Q_RTX/V1_RTX)*RTX_C - (Q_RTX/V2_RTX)*RTX_P;

// ---- B-cell / Plasma-cell PD ----
// Rituximab kill (saturable in Cp)
double E_RTX = KRTX_KILL * RTX_C_C / (KC50_RTX + RTX_C_C);
// Prednisone effect on B-cell production
double E_PRED = EMAX_PRED * PRED_C_C / (EC50_PRED + PRED_C_C);
// MMF/steroid placeholder reduction in production (here uses pred)
double PROD_B = KPROD_B * (1 - E_PRED);

dxdt_BNAIVE = PROD_B - KDEG_B*BNAIVE - E_RTX*BNAIVE - KMAT_PB*BNAIVE/(B0/100);
dxdt_PBLAST = KMAT_PB*BNAIVE/(B0/100) - KDEG_PB*PBLAST - 0.5*E_RTX*PBLAST;
dxdt_LLPC   = KMAT_LL*PBLAST - KDEG_LL*LLPC;

// ---- Anti-GAD65 antibody (per day terms) ----
// IVIG anti-idiotype neutralization
double IVIG_CONC = IVIG_C / V_IGG;        // g/L
dxdt_ABAB_S   = KPROD_AB_PB*PBLAST + KPROD_AB_LL*LLPC
                 - KDEG_AB*ABAB_S
                 - IVIG_NEUT*IVIG_CONC*ABAB_S
                 - KCSF_AB*ABAB_S;
dxdt_ABAB_CSF = KCSF_AB*ABAB_S - KCSF_AB_OUT*ABAB_CSF;

// ---- GAD65 enzyme activity (CSF Ab-driven) ----
dxdt_GAD  = KIN_GAD*GAD0/100 - KOUT_GAD*GAD - AB_KILL*ABAB_CSF*GAD/100;

// ---- GABA pool driven by GAD activity ----
dxdt_GABA = KIN_GABA * (GAD/GAD0) * 100 - KOUT_GABA*GABA;

// ---- alpha-MN excitability ----
double GABA_EFF  = EMAX_GABA * GABA / (EC50_GABA + GABA);
double BZD_EFF   = EMAX_BZD  * DIAZ_FREE_BRAIN / (EC50_BZD + DIAZ_FREE_BRAIN);
double BAC_EFF   = EMAX_BAC  * BAC_CSF_C / (EC50_BAC + BAC_CSF_C);
double GAB_EFF   = EMAX_GAB  * GAB_C_C   / (EC50_GAB + GAB_C_C);
double INH       = GABA_EFF + BZD_EFF + BAC_EFF + GAB_EFF;
if(INH > 0.95) INH = 0.95;

dxdt_MN_EXC = KIN_MN*(1 - INH) - KOUT_MN*MN_EXC;

// ---- Stiffness ----
double EXC_DRIVE = (MN_EXC > MN_THR) ? (MN_EXC - MN_THR) : 0;
dxdt_STIFF = KIN_STIFF*EXC_DRIVE - KOUT_STIFF*STIFF;
if(STIFF < 0) STIFF = 0;
if(STIFF > 100) STIFF = 100;

// ---- Spasms ----
dxdt_SPASM = KIN_SPASM*EXC_DRIVE/30 - KOUT_SPASM*SPASM;
if(SPASM < 0) SPASM = 0;

// ---- Bone (steroid effect) ----
dxdt_BMD = -KBMD_PRED * PRED_C_C * 1000;   // mg/L scale

// ---- HSI surrogate ----
dxdt_HSI = 0.5*STIFF - 0.1*HSI;

$TABLE
double diaz_total_ngml = (DIAZ_C/V1_DIAZ)*1000;
double dmd_ngml        = (DMD_C/V_DMD)*1000;
double bac_csf_ngml    = BAC_CSF_C*1000;
double gab_mgL         = GAB_C_C;
double rtx_ugml        = RTX_C_C;
double abab_pct        = 100*ABAB_S/ABAB0;
double gad_pct         = GAD;
double gaba_pct        = GABA;
double mn_pct          = MN_EXC;
double inh_total       = INH;

$CAPTURE
diaz_total_ngml dmd_ngml bac_csf_ngml gab_mgL rtx_ugml
abab_pct gad_pct gaba_pct mn_pct inh_total
STIFF SPASM BMD HSI BNAIVE PBLAST LLPC ABAB_S ABAB_CSF
'

sps_mod <- mread("sps", code = sps_code)

## --------------------------------------------------------------
## Steady-state seeding: run model from disease-OFF to disease-ON
## (illustrative): for the public release we seed disease at t=0
## with elevated anti-GAD65 (ABAB0) and partially depressed GAD.
## --------------------------------------------------------------

## ---------- Helper to set initial conditions for "active SPS" --
sps_active_init <- function(mod, severity = "moderate") {
  init_list <- switch(severity,
    mild = list(ABAB_S = 4000,  GAD = 80, GABA = 80, MN_EXC = 120, STIFF = 25, SPASM = 0.5),
    moderate = list(ABAB_S = 10000, GAD = 55, GABA = 55, MN_EXC = 150, STIFF = 55, SPASM = 1.5),
    severe = list(ABAB_S = 25000, GAD = 35, GABA = 35, MN_EXC = 180, STIFF = 80, SPASM = 4),
    crisis = list(ABAB_S = 40000, GAD = 25, GABA = 25, MN_EXC = 200, STIFF = 92, SPASM = 7)
  )
  init(mod, init_list)
}

## ---------- Build event objects for each scenario --------------
ev_diaz <- function(daily_mg, days = 28) {
  ev(amt = daily_mg/3, ii = 8, addl = ceiling(days*3) - 1, cmt = "DIAZ_GUT", time = 0)
}
ev_baclo_oral <- function(daily_mg, days = 28) {
  ev(amt = daily_mg/3, ii = 8, addl = ceiling(days*3) - 1, cmt = "BAC_GUT", time = 0)
}
ev_baclo_IT <- function(daily_mcg, days = 28) {
  ev(amt = daily_mcg/1000, ii = 0.001, addl = 0, cmt = "BAC_CSF",   # continuous infusion approx
     time = 0, rate = (daily_mcg/1000)/24)
}
ev_gabapentin <- function(daily_mg, days = 28) {
  ev(amt = daily_mg/3, ii = 8, addl = ceiling(days*3) - 1, cmt = "GAB_GUT", time = 0)
}
ev_ivig_cycle <- function(total_g_per_cycle = 140, n_cycles = 6, ii_d = 28) {
  ev(amt = total_g_per_cycle/2, ii = 24*ii_d, addl = n_cycles - 1,
     cmt = "IVIG_C", time = 24)   # 2 g/kg over 2 d -> simplified daily-bolus
}
ev_rituximab <- function(mg = 1000, ii_d = 14, addl = 1) {
  ev(amt = mg, ii = 24*ii_d, addl = addl, cmt = "RTX_C", time = 0, rate = mg/4)
}
ev_prednisone <- function(daily_mg, days = 14) {
  ev(amt = daily_mg, ii = 24, addl = days - 1, cmt = "PRED_GUT", time = 0)
}
ev_plex_event <- function(fractions = c(0.6, 0.4, 0.3, 0.25, 0.2),
                          start_d = 0, interval_d = 2) {
  # Use observation events with bioavail modifier? Simpler: dose negative? Not allowed.
  # Implementation: simulate with rate of removal via large transient KEL bumps.
  # Here we return a data frame of "events" that the simulator post-processes.
  data.frame(
    time = start_d*24 + (seq_along(fractions)-1)*interval_d*24,
    fraction_remaining = fractions
  )
}

## --------------------------------------------------------------
## Helper to apply PLEX (manual override of ABAB_S between events)
## --------------------------------------------------------------
sps_apply_plex <- function(mod, plex_tbl, end_d = 90, ev_other = NULL,
                           severity = "severe") {
  mod_i <- sps_active_init(mod, severity)
  out_full <- NULL
  prev_t <- 0
  for(i in seq_len(nrow(plex_tbl))){
    t_i <- plex_tbl$time[i]
    if(t_i > prev_t){
      out_seg <- mrgsim(mod_i, ev = ev_other, end = t_i, delta = 0.5)
      out_full <- rbind(out_full, as.data.frame(out_seg))
      tail_state <- tail(as.data.frame(out_seg), 1)
      # Apply removal: multiply ABAB_S and ABAB_CSF
      f <- plex_tbl$fraction_remaining[i]
      mod_i <- init(mod_i,
                    ABAB_S   = tail_state$ABAB_S * f,
                    ABAB_CSF = tail_state$ABAB_CSF * f)
    }
    prev_t <- t_i
  }
  out_seg <- mrgsim(mod_i, ev = ev_other, end = end_d*24, delta = 0.5)
  out_full <- rbind(out_full, as.data.frame(out_seg))
  out_full
}

## --------------------------------------------------------------
## Scenario factory
## --------------------------------------------------------------
sps_run <- function(scenario = c("dx_bzd_only",
                                 "bzd_baclofen_combo",
                                 "bzd_ivig_q4w",
                                 "rtx_induction",
                                 "plex_rescue",
                                 "intrathecal_baclofen"),
                    severity = "moderate",
                    horizon_d = 180) {
  scenario <- match.arg(scenario)
  mod      <- sps_active_init(sps_mod, severity)

  if(scenario == "dx_bzd_only"){
    ev_in <- ev_diaz(daily_mg = 30, days = horizon_d)
    mrgsim(mod, ev = ev_in, end = horizon_d*24, delta = 0.25) %>% as.data.frame()

  } else if(scenario == "bzd_baclofen_combo"){
    ev_in <- c(ev_diaz(20, horizon_d), ev_baclo_oral(60, horizon_d))
    mrgsim(mod, ev = ev_in, end = horizon_d*24, delta = 0.25) %>% as.data.frame()

  } else if(scenario == "bzd_ivig_q4w"){
    ev_in <- c(ev_diaz(15, horizon_d),
               ev_ivig_cycle(total_g_per_cycle = 140,
                             n_cycles = ceiling(horizon_d/28),
                             ii_d = 28))
    mrgsim(mod, ev = ev_in, end = horizon_d*24, delta = 0.25) %>% as.data.frame()

  } else if(scenario == "rtx_induction"){
    ev_in <- c(ev_diaz(15, horizon_d),
               ev_rituximab(mg = 1000, ii_d = 14, addl = 1),
               ev_prednisone(100, 14))
    mrgsim(mod, ev = ev_in, end = horizon_d*24, delta = 0.25) %>% as.data.frame()

  } else if(scenario == "plex_rescue"){
    ev_in <- ev_diaz(20, horizon_d)
    plex_tbl <- ev_plex_event(c(0.55,0.45,0.40,0.35,0.30), start_d = 0, interval_d = 2)
    sps_apply_plex(sps_mod, plex_tbl, end_d = horizon_d,
                   ev_other = ev_in, severity = severity)

  } else if(scenario == "intrathecal_baclofen"){
    ev_in <- c(ev_diaz(10, horizon_d),
               ev_baclo_IT(daily_mcg = 400, days = horizon_d))
    mrgsim(mod, ev = ev_in, end = horizon_d*24, delta = 0.25) %>% as.data.frame()
  }
}

## --------------------------------------------------------------
## QC / smoke test (only runs if invoked directly)
## --------------------------------------------------------------
if(sys.nframe() == 0L && interactive()){
  res <- sps_run("bzd_ivig_q4w", severity = "moderate", horizon_d = 84)
  ggplot(res, aes(time/24)) +
    geom_line(aes(y = STIFF, color = "Stiffness")) +
    geom_line(aes(y = SPASM*10, color = "Spasm x10")) +
    geom_line(aes(y = abab_pct, color = "Anti-GAD65 %")) +
    geom_line(aes(y = gaba_pct,  color = "GABA %")) +
    labs(x = "Days", y = "Magnitude", color = "Variable",
         title = "SPS QSP — IVIG q4w on moderate disease") +
    theme_minimal()
}
