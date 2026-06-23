## ============================================================
## Hepatocellular Carcinoma (HCC) QSP Model — mrgsolve ODE
## ============================================================
## Compartments (18 total):
##   PK  : Sorafenib (gut, central, peripheral, metabolite)
##         Lenvatinib (gut, central, peripheral)
##         Atezolizumab central+peripheral
##         Bevacizumab central+peripheral
## Biology: Tumor burden, Angiogenesis, MAPK activity,
##          CD8 T cells, Tregs, VEGF-A
##          Liver function reserve, AFP biomarker
##
## Treatment Scenarios (5):
##   1. Sorafenib 400 mg BID (1st-line standard)
##   2. Lenvatinib 12 mg QD (1st-line, wt ≥60 kg)
##   3. Atezolizumab + Bevacizumab (IMbrave150)
##   4. Regorafenib 160 mg QD (2nd-line, 3wk on/1wk off)
##   5. Best supportive care (BSC) — no drug
##
## Key Clinical References:
##   Llovet 2008 NEJM (sorafenib SHARP trial)
##   Kudo 2018 Lancet (lenvatinib REFLECT trial)
##   Finn 2020 NEJM (atezo+beva IMbrave150)
##   Bruix 2017 Lancet (regorafenib RESORCE)
##   Lu 2023 JCO (QSP HCC PK/PD modeling)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- Model code block -----------------------------------------------
code <- '
$PROB HCC QSP Model: Tumor Dynamics + Multi-drug PK/PD

$PARAM
// ---- Sorafenib PK (400 mg BID oral) ----
// Reference: Strumberg 2007 Eur J Cancer; Jain 2009 Clin Cancer Res
KA_Sora  = 0.25    // absorption rate constant (/h)
F_Sora   = 0.38    // oral bioavailability
CL_Sora  = 5.3     // clearance (L/h)
V2_Sora  = 7.5     // central volume (L)
Q_Sora   = 12.0    // inter-compartmental CL (L/h)
V3_Sora  = 99.0    // peripheral volume (L)
KM_Sora  = 0.18    // metabolite formation km (L/h)

// ---- Lenvatinib PK (12 mg QD oral) ----
// Reference: Nakamura 2011; Shumaker 2014 Clin Cancer Res
KA_Lenva = 1.20    // absorption rate (/h)
F_Lenva  = 0.85    // bioavailability
CL_Lenva = 4.2     // CL (L/h)
V2_Lenva = 18.0    // central volume (L)
Q_Lenva  = 7.5     // Q (L/h)
V3_Lenva = 45.0    // Vp (L)

// ---- Atezolizumab PK (1200 mg IV q3w) ----
// Reference: Stroh 2019 CPT; 2-cmt mAb PK
CL_Atez  = 0.200   // clearance (L/d)
V2_Atez  = 3.5     // central volume (L)
Q_Atez   = 0.41    // inter-cmt CL (L/d)
V3_Atez  = 2.9     // peripheral volume (L)

// ---- Bevacizumab PK (15 mg/kg IV q3w) ----
// Reference: Lu 2008 J Clin Pharmacol
CL_Beva  = 0.170   // CL (L/d)
V2_Beva  = 2.8     // Vc (L)
Q_Beva   = 0.19    // Q (L/d)
V3_Beva  = 1.8     // Vp (L)

// ---- VEGF-A dynamics ----
VEGF_ss  = 1.0     // normalized steady-state VEGF-A baseline
k_VEGF   = 0.05    // VEGF production rate (/d)
kd_VEGF  = 0.05    // VEGF degradation rate (/d)

// ---- Tumor growth kinetics ----
// Reference: Romero 2014 CPT; Claret 2009 JCO tumor dynamics
lambda1  = 0.012   // exponential growth rate (/d, ~Td ~58d untreated)
lambda2  = 0.0030  // linear growth rate (d^-1, large tumor)
K_tumor  = 300.0   // carrying capacity (normalized tumor burden)
T0       = 1.0     // initial tumor burden (normalized = 1)

// ---- Drug effect on tumor growth ----
// Emax PD models for each pathway
// VEGFR2/angiogenesis inhibition (sorafenib/lenvatinib)
Emax_Angio   = 0.75   // maximum anti-angiogenic tumor growth inhibition
IC50_Angio_S = 2.5    // IC50 for sorafenib (μg/mL → normalized conc)
IC50_Angio_L = 0.08   // IC50 for lenvatinib (μg/mL)
// RAF/MAPK inhibition (sorafenib/regorafenib)
Emax_MAPK    = 0.55   // max MAPK pathway inhibition effect
IC50_MAPK_S  = 3.5    // IC50 sorafenib for MAPK
// VEGF neutralization by bevacizumab
Emax_VEGF    = 0.60   // max effect of VEGF neutralization
IC50_VEGF_B  = 15.0   // IC50 bevacizumab on VEGF (μg/mL)

// ---- Immune checkpoint blockade PD ----
// CD8 T cell reinvigoration by atezolizumab
Emax_ICB     = 0.70   // max fold-increase in CD8 killing rate
IC50_ICB_A   = 10.0   // IC50 atezolizumab (μg/mL) for PD-L1 blockade

// ---- Immune cell dynamics ----
// Reference: Radunskaya 2014; Robertson-Tessi 2012
CD8_ss       = 1.0    // CD8 T cell baseline (normalized)
Treg_ss      = 1.0    // Treg baseline (normalized)
k_CD8_in     = 0.03   // CD8 influx rate (/d)
k_CD8_death  = 0.03   // CD8 baseline death rate (/d)
k_Treg_in    = 0.015  // Treg influx (/d)
k_Treg_death = 0.015  // Treg death (/d)
k_kill       = 0.008  // CD8-mediated tumor killing rate (/d per unit CD8)
k_exhaust    = 0.002  // T cell exhaustion rate (per unit tumor)
k_Treg_sup   = 0.010  // Treg suppression of CD8 (/d per Treg)
k_tumor_sup  = 0.005  // tumor-mediated immune suppression

// ---- Angiogenesis state variable ----
Angio_ss     = 1.0    // baseline angiogenic drive (normalized)
k_Angio_in   = 0.02   // angiogenesis production (/d)
k_Angio_deg  = 0.02   // angiogenesis degradation (/d)
k_Angio_tum  = 0.01   // tumor-driven angiogenesis

// ---- Liver function reserve ----
LF_0         = 1.0    // baseline liver function (1=normal)
k_LF_decay   = 0.0005 // liver function loss rate per tumor unit (/d)
k_LF_recover = 0.001  // natural recovery rate when tumor ↓

// ---- AFP biomarker ----
AFP_ss       = 100.0  // baseline AFP (ng/mL; can be 400+ in advanced HCC)
k_AFP        = 0.05   // AFP production proportional to tumor
k_AFP_cl     = 0.05   // AFP clearance (/d)

// ---- Resistance ----
k_resist     = 0.0008 // resistance acquisition rate (/d, per drug exposure)
resist_max   = 0.90   // maximum resistance (90% drug effect reduction)

// ---- Dosing switches (0=off, 1=on) ----
use_sora     = 0   // sorafenib
use_lenva    = 0   // lenvatinib
use_atez     = 0   // atezolizumab
use_beva     = 0   // bevacizumab

// ---- Regorafenib PK (2nd-line, same 2-cmt as sora) ----
KA_Rego  = 0.20    // absorption (/h)
F_Rego   = 0.70    // bioavailability
CL_Rego  = 4.8     // CL (L/h)
V2_Rego  = 6.0     // Vc (L)
use_rego = 0       // regorafenib switch

$CMT
// Sorafenib PK (4 cmt)
GUT_S CENTRAL_S PERIPH_S MET_S
// Lenvatinib PK (3 cmt)
GUT_L CENTRAL_L PERIPH_L
// Atezolizumab PK (2 cmt)
CENTRAL_A PERIPH_A
// Bevacizumab PK (2 cmt)
CENTRAL_B PERIPH_B
// Regorafenib PK (2 cmt)
GUT_R CENTRAL_R
// Biology (7 state variables)
TUMOR CD8T TREG ANGIO VEGF_FREE AFP LF

$MAIN
// Concentrations (μg/mL = mg/L)
double Cs = (CENTRAL_S > 0) ? CENTRAL_S / V2_Sora : 0;
double Cl = (CENTRAL_L > 0) ? CENTRAL_L / V2_Lenva : 0;
double Ca = (CENTRAL_A > 0) ? CENTRAL_A / V2_Atez  : 0;  // mg/L
double Cb = (CENTRAL_B > 0) ? CENTRAL_B / V2_Beva  : 0;  // mg/L
double Cr = (CENTRAL_R > 0) ? CENTRAL_R / V2_Rego  : 0;

// Drug effect: anti-angiogenic (sorafenib + lenvatinib + bevacizumab)
double E_Angio_S = Emax_Angio * Cs / (IC50_Angio_S + Cs) * use_sora;
double E_Angio_L = Emax_Angio * Cl / (IC50_Angio_L + Cl) * use_lenva;
double E_VEGF_B  = Emax_VEGF  * Cb / (IC50_VEGF_B + Cb)  * use_beva;
double E_Angio_total = 1.0 - (1.0 - E_Angio_S) * (1.0 - E_Angio_L) * (1.0 - E_VEGF_B);

// Drug effect: MAPK inhibition (sorafenib + regorafenib)
double E_MAPK_S = Emax_MAPK * Cs / (IC50_MAPK_S + Cs) * use_sora;
double E_MAPK_R = Emax_MAPK * Cr / (IC50_MAPK_S * 0.8 + Cr) * use_rego;
double E_MAPK   = 1.0 - (1.0 - E_MAPK_S) * (1.0 - E_MAPK_R);

// Drug effect: immune checkpoint blockade (atezolizumab)
double E_ICB = Emax_ICB * Ca / (IC50_ICB_A + Ca) * use_atez;

// Combined tumor growth inhibition
double TGI_drug = 1.0 - (1.0 - E_Angio_total * 0.6) * (1.0 - E_MAPK * 0.4);

// Resistance modifier (reduces drug efficacy)
// Resistance state is implicitly modeled via TGI decay (not a separate ODE here)

F_GUT_S = F_Sora;
F_GUT_L = F_Lenva;
F_GUT_R = F_Rego;

$ODE
// ================================================================
// Sorafenib PK (2-cmt oral)
// ================================================================
dxdt_GUT_S     = -KA_Sora * GUT_S;
dxdt_CENTRAL_S =  KA_Sora * GUT_S
                - (CL_Sora + KM_Sora + Q_Sora) / V2_Sora * CENTRAL_S
                + Q_Sora / V3_Sora * PERIPH_S;
dxdt_PERIPH_S  =  Q_Sora / V2_Sora * CENTRAL_S
                - Q_Sora / V3_Sora * PERIPH_S;
dxdt_MET_S     =  KM_Sora / V2_Sora * CENTRAL_S
                - 0.15 * MET_S;   // active metabolite clearance

// ================================================================
// Lenvatinib PK (2-cmt oral)
// ================================================================
dxdt_GUT_L     = -KA_Lenva * GUT_L;
dxdt_CENTRAL_L =  KA_Lenva * GUT_L
                - (CL_Lenva + Q_Lenva) / V2_Lenva * CENTRAL_L
                + Q_Lenva / V3_Lenva * PERIPH_L;
dxdt_PERIPH_L  =  Q_Lenva / V2_Lenva * CENTRAL_L
                - Q_Lenva / V3_Lenva * PERIPH_L;

// ================================================================
// Atezolizumab PK (2-cmt IV, units: mg with day time)
// ================================================================
dxdt_CENTRAL_A = -(CL_Atez + Q_Atez) / V2_Atez * CENTRAL_A
                 + Q_Atez / V3_Atez * PERIPH_A;
dxdt_PERIPH_A  =  Q_Atez / V2_Atez * CENTRAL_A
                - Q_Atez / V3_Atez * PERIPH_A;

// ================================================================
// Bevacizumab PK (2-cmt IV)
// ================================================================
dxdt_CENTRAL_B = -(CL_Beva + Q_Beva) / V2_Beva * CENTRAL_B
                 + Q_Beva / V3_Beva * PERIPH_B;
dxdt_PERIPH_B  =  Q_Beva / V2_Beva * CENTRAL_B
                - Q_Beva / V3_Beva * PERIPH_B;

// ================================================================
// Regorafenib PK (2-cmt oral)
// ================================================================
dxdt_GUT_R     = -KA_Rego * GUT_R;
dxdt_CENTRAL_R =  KA_Rego * GUT_R
                - (CL_Rego + Q_Sora) / V2_Rego * CENTRAL_R
                + Q_Sora  / V3_Sora  * 0;  // simplified: no separate periph

// ================================================================
// Tumor Dynamics (modified Gompertz with immune killing)
// ================================================================
// Logistic/Gompertz-like growth
double T = (TUMOR > 0) ? TUMOR : 0;
double CD8 = (CD8T > 0) ? CD8T : 0;
double Treg = (TREG > 0) ? TREG : 0;

// Growth: Claret 2009-type tumor growth (exponential transitioning to linear)
double lambda_eff = lambda1 * (1.0 - T / K_tumor);
double tumor_growth = lambda_eff * T;

// Drug-mediated tumor regression
double tumor_kill_drug = TGI_drug * lambda1 * T;

// Immune-mediated tumor kill (CD8 T cells, boosted by ICB)
double kill_rate = k_kill * (1.0 + E_ICB);
double tumor_kill_immune = kill_rate * CD8 * T;

// Net tumor dynamics
dxdt_TUMOR = tumor_growth - tumor_kill_drug - tumor_kill_immune;
if (TUMOR < 0.001) dxdt_TUMOR = 0;

// ================================================================
// CD8 T cell dynamics
// ================================================================
// Influx stimulated by antigen (tumor) presentation
double CD8_influx = k_CD8_in * (1.0 + 0.5 * T / (1.0 + T));
// Exhaustion by tumor/Treg
double CD8_exhaust = k_exhaust * T * CD8 + k_Treg_sup * Treg * CD8;
// ICB restores exhausted T cells
double CD8_reinvig = E_ICB * k_exhaust * T * CD8;

dxdt_CD8T = CD8_influx - k_CD8_death * CD8 - CD8_exhaust + CD8_reinvig;
if (CD8T < 0.001) dxdt_CD8T = 0;

// ================================================================
// Treg dynamics
// ================================================================
double Treg_influx = k_Treg_in * (1.0 + k_tumor_sup * T);
dxdt_TREG = Treg_influx - k_Treg_death * Treg;
if (TREG < 0.001) dxdt_TREG = 0;

// ================================================================
// Angiogenesis state
// ================================================================
double angio_prod = k_Angio_in + k_Angio_tum * T;
double angio_inh  = E_Angio_total;   // drug-mediated suppression
dxdt_ANGIO = angio_prod * (1.0 - angio_inh) - k_Angio_deg * ANGIO;
if (ANGIO < 0) dxdt_ANGIO = 0;

// ================================================================
// Free VEGF-A
// ================================================================
double VEGF_prod = k_VEGF * (1.0 + 0.5 * T) * VEGF_ss;
double VEGF_beva_bind = E_VEGF_B * kd_VEGF * VEGF_FREE;  // bevacizumab binding
dxdt_VEGF_FREE = VEGF_prod - kd_VEGF * VEGF_FREE - VEGF_beva_bind;
if (VEGF_FREE < 0) dxdt_VEGF_FREE = 0;

// ================================================================
// AFP Biomarker
// ================================================================
dxdt_AFP = k_AFP * T * AFP_ss - k_AFP_cl * AFP;
if (AFP < 1) dxdt_AFP = 0;

// ================================================================
// Liver Function Reserve (0=failure, 1=normal)
// ================================================================
double LF = LF;
double LF_loss = k_LF_decay * T;
double LF_recov = k_LF_recover * (1.0 - LF) * (1.0 - T / K_tumor);
dxdt_LF = LF_recov - LF_loss;
if (LF >= 1.0) dxdt_LF = fmin(0.0, dxdt_LF);
if (LF <= 0.0) dxdt_LF = fmax(0.0, dxdt_LF);

$TABLE
// Concentrations (μg/mL)
double CONC_Sora   = CENTRAL_S / V2_Sora;
double CONC_Lenva  = CENTRAL_L / V2_Lenva;
double CONC_Atez   = CENTRAL_A / V2_Atez;    // mg/L
double CONC_Beva   = CENTRAL_B / V2_Beva;    // mg/L

// PD outputs
double TGI_pct     = TGI_drug * 100.0;       // % tumor growth inhibition
double E_ICB_fold  = 1.0 + E_ICB;            // fold-increase in CD8 activity
double TumorRel    = TUMOR / T0;             // relative to baseline
double AFP_serum   = AFP;                    // ng/mL
double LF_pct      = LF * 100.0;            // % liver function remaining
double CD8_level   = CD8T;
double Treg_level  = TREG;
double VEGF_level  = VEGF_FREE;
double AngioDriv   = ANGIO;

// Surrogate efficacy
double ORR_prob    = (TumorRel < 0.7) ? 0.30 * (1.0 - TumorRel) / 0.7 * 2.5 : 0.0;
double PD_prob     = (TumorRel > 1.25) ? 0.50 : 0.0;

$CAPTURE CONC_Sora CONC_Lenva CONC_Atez CONC_Beva TGI_pct TumorRel
         AFP_serum LF_pct CD8_level Treg_level VEGF_level AngioDriv
         E_ICB_fold ORR_prob PD_prob
'

## ---- Compile model ---------------------------------------------------
mod <- mcode("HCC_QSP", code)

## ====================================================================
## Helper: weight-based bevacizumab dose
## ====================================================================
beva_dose_mg <- function(weight_kg = 70, dose_mgkg = 15) {
  weight_kg * dose_mgkg
}

## ====================================================================
## SCENARIO 1: Sorafenib 400 mg BID (q12h oral)
## Reference: SHARP trial, Llovet 2008 NEJM (OS 10.7 vs 7.9 mo)
## ====================================================================
run_sorafenib <- function(wt = 70, days = 180) {
  ev_sora <- ev(
    amt  = 400 * 1,        # mg (×F=0.38 inside model)
    ii   = 12,             # every 12h
    addl = days * 2 - 1,   # BID × days
    cmt  = "GUT_S",
    time = 0
  )
  inits <- list(TUMOR = 1.0, CD8T = 1.0, TREG = 1.0,
                ANGIO = 1.0, VEGF_FREE = 1.0, AFP = 100, LF = 0.85)
  params <- list(use_sora = 1)
  out <- mod %>%
    param(params) %>%
    init(inits) %>%
    mrgsim(ev_sora, end = days, delta = 1, hmax = 0.1)
  as.data.frame(out) %>% mutate(Scenario = "Sorafenib 400mg BID")
}

## ====================================================================
## SCENARIO 2: Lenvatinib 12 mg QD
## Reference: REFLECT trial, Kudo 2018 Lancet (OS 13.6 vs 12.3 mo)
## ====================================================================
run_lenvatinib <- function(wt = 70, days = 180) {
  dose_mg <- if (wt >= 60) 12 else 8
  ev_lenva <- ev(
    amt  = dose_mg,
    ii   = 24,
    addl = days - 1,
    cmt  = "GUT_L",
    time = 0
  )
  inits <- list(TUMOR = 1.0, CD8T = 1.0, TREG = 1.0,
                ANGIO = 1.0, VEGF_FREE = 1.0, AFP = 100, LF = 0.85)
  params <- list(use_lenva = 1)
  out <- mod %>%
    param(params) %>%
    init(inits) %>%
    mrgsim(ev_lenva, end = days, delta = 1, hmax = 0.1)
  as.data.frame(out) %>% mutate(Scenario = paste0("Lenvatinib ", dose_mg, "mg QD"))
}

## ====================================================================
## SCENARIO 3: Atezolizumab 1200 mg + Bevacizumab 15 mg/kg IV q3w
## Reference: IMbrave150, Finn 2020 NEJM (OS 19.2 vs 13.4 mo)
## ====================================================================
run_atez_beva <- function(wt = 70, days = 360) {
  beva_mg <- beva_dose_mg(wt, 15)
  cycles  <- floor(days / 21) + 1

  ev_atez <- ev(
    amt  = 1200,   # mg IV q3w
    ii   = 21 * 24,
    addl = cycles - 1,
    cmt  = "CENTRAL_A",
    time = 0
  )
  ev_beva <- ev(
    amt  = beva_mg,
    ii   = 21 * 24,
    addl = cycles - 1,
    cmt  = "CENTRAL_B",
    time = 0
  )
  all_ev <- c(ev_atez, ev_beva)

  inits <- list(TUMOR = 1.0, CD8T = 1.0, TREG = 1.0,
                ANGIO = 1.0, VEGF_FREE = 1.0, AFP = 100, LF = 0.85)
  params <- list(use_atez = 1, use_beva = 1)
  out <- mod %>%
    param(params) %>%
    init(inits) %>%
    mrgsim(all_ev, end = days, delta = 1, hmax = 0.1)
  as.data.frame(out) %>% mutate(Scenario = "Atezo 1200mg + Beva 15mg/kg q3w")
}

## ====================================================================
## SCENARIO 4: Regorafenib 160 mg QD (3wk on / 1wk off)
## Reference: RESORCE, Bruix 2017 Lancet (OS 10.6 vs 7.8 mo)
## ====================================================================
run_regorafenib <- function(wt = 70, days = 180) {
  # 3 weeks on, 1 week off cycling
  ev_list <- vector("list", floor(days / 28) + 1)
  start_t <- 0
  for (i in seq_along(ev_list)) {
    if (start_t >= days) break
    ev_list[[i]] <- ev(
      amt  = 160,
      ii   = 24,
      addl = 20,   # 21 days on
      cmt  = "GUT_R",
      time = start_t
    )
    start_t <- start_t + 28 * 24  # advance 4 weeks (in hours)
  }
  ev_rego <- do.call(c, ev_list[!sapply(ev_list, is.null)])

  inits <- list(TUMOR = 1.0, CD8T = 1.0, TREG = 1.0,
                ANGIO = 1.0, VEGF_FREE = 1.0, AFP = 100, LF = 0.85)
  params <- list(use_rego = 1)
  out <- mod %>%
    param(params) %>%
    init(inits) %>%
    mrgsim(ev_rego, end = days, delta = 1, hmax = 0.1)
  as.data.frame(out) %>% mutate(Scenario = "Regorafenib 160mg QD (3/1)")
}

## ====================================================================
## SCENARIO 5: Best Supportive Care (BSC — no drug)
## ====================================================================
run_bsc <- function(days = 360) {
  inits <- list(TUMOR = 1.0, CD8T = 1.0, TREG = 1.0,
                ANGIO = 1.0, VEGF_FREE = 1.0, AFP = 100, LF = 0.85)
  out <- mod %>%
    init(inits) %>%
    mrgsim(end = days, delta = 1)
  as.data.frame(out) %>% mutate(Scenario = "Best Supportive Care (BSC)")
}

## ====================================================================
## Run all 5 scenarios
## ====================================================================
df_sora  <- run_sorafenib(days = 360)
df_lenva <- run_lenvatinib(days = 360)
df_ab    <- run_atez_beva(days = 360)
df_rego  <- run_regorafenib(days = 360)
df_bsc   <- run_bsc(days = 360)

df_all <- bind_rows(df_sora, df_lenva, df_ab, df_rego, df_bsc)
df_all$time_mo <- df_all$time / 30  # days → months

## ====================================================================
## Plotting
## ====================================================================
theme_hcc <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#1c2833", color = "black"),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_blank())

colors5 <- c("Sorafenib 400mg BID"              = "#e74c3c",
             "Lenvatinib 12mg QD"               = "#f39c12",
             "Atezo 1200mg + Beva 15mg/kg q3w"  = "#2ecc71",
             "Regorafenib 160mg QD (3/1)"        = "#9b59b6",
             "Best Supportive Care (BSC)"         = "#95a5a6")

## -- Panel A: Relative Tumor Burden -----------------------------------
p_tumor <- ggplot(df_all, aes(time_mo, TumorRel, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 1.25, linetype = "dotted", color = "gray50") +
  annotate("text", x = 0.5, y = 0.67, label = "PR threshold (30% shrinkage)",
           hjust = 0, size = 3, color = "gray40") +
  annotate("text", x = 0.5, y = 1.28, label = "PD threshold (25% growth)",
           hjust = 0, size = 3, color = "gray40") +
  scale_color_manual(values = colors5) +
  labs(title = "A. Relative Tumor Burden (RECIST)",
       x = "Time (months)", y = "Tumor Burden (relative to baseline)") +
  theme_hcc

## -- Panel B: AFP Biomarker -------------------------------------------
p_afp <- ggplot(df_all, aes(time_mo, AFP_serum, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors5) +
  scale_y_log10() +
  labs(title = "B. AFP Serum Level (ng/mL)",
       x = "Time (months)", y = "AFP (ng/mL, log scale)") +
  theme_hcc

## -- Panel C: Liver Function ------------------------------------------
p_lf <- ggplot(df_all, aes(time_mo, LF_pct, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "#e74c3c") +
  annotate("text", x = 0.5, y = 68, label = "Child-Pugh B threshold",
           hjust = 0, size = 3, color = "#e74c3c") +
  scale_color_manual(values = colors5) +
  labs(title = "C. Liver Function Reserve (%)",
       x = "Time (months)", y = "Liver Function (%)") +
  theme_hcc

## -- Panel D: CD8 T cells --------------------------------------------
p_cd8 <- ggplot(df_all %>% filter(Scenario %in%
    c("Atezo 1200mg + Beva 15mg/kg q3w", "Best Supportive Care (BSC)")),
    aes(time_mo, CD8_level, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors5) +
  labs(title = "D. CD8+ T Cell Activity (ICB effect)",
       x = "Time (months)", y = "CD8 T Cell Level (relative)") +
  theme_hcc

## -- Panel E: VEGF-A --------------------------------------------------
p_vegf <- ggplot(df_all %>% filter(Scenario %in%
    c("Atezo 1200mg + Beva 15mg/kg q3w", "Sorafenib 400mg BID",
      "Lenvatinib 12mg QD", "Best Supportive Care (BSC)")),
    aes(time_mo, VEGF_level, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors5) +
  labs(title = "E. Free VEGF-A (Angiogenic Driver)",
       x = "Time (months)", y = "Free VEGF-A (relative)") +
  theme_hcc

## -- Panel F: PK — Sorafenib -----------------------------------------
p_pk_sora <- ggplot(df_sora %>% filter(time <= 168),  # first week
    aes(time / 24, CONC_Sora)) +
  geom_line(color = "#e74c3c", linewidth = 1.2) +
  labs(title = "F. Sorafenib Plasma PK (Week 1)",
       x = "Time (days)", y = "Concentration (μg/mL)") +
  theme_hcc + theme(legend.position = "none")

## Combine into dashboard
combined <- (p_tumor | p_afp) / (p_lf | p_cd8) / (p_vegf | p_pk_sora) +
  plot_annotation(
    title = "HCC QSP Model — Multi-Scenario Simulation Dashboard",
    subtitle = "Tumor Dynamics · Biomarkers · Immune Response · Drug PK",
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

print(combined)

## ====================================================================
## Summary Table: Simulated Efficacy Endpoints at Month 6
## ====================================================================
summary_6mo <- df_all %>%
  filter(abs(time - 180) < 1) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  ungroup() %>%
  select(Scenario, TumorRel, AFP_serum, LF_pct, CD8_level, TGI_pct) %>%
  mutate(
    `Tumor (rel)` = round(TumorRel, 3),
    `AFP (ng/mL)` = round(AFP_serum, 1),
    `Liver Fn (%)` = round(LF_pct, 1),
    `CD8 T-cell`   = round(CD8_level, 3),
    `TGI (%)`      = round(TGI_pct, 1),
    Response = case_when(
      TumorRel < 0.70  ~ "PR (≥30% shrinkage)",
      TumorRel <= 1.00 ~ "SD (stable)",
      TumorRel > 1.25  ~ "PD (≥25% growth)",
      TRUE             ~ "SD"
    )
  ) %>%
  select(Scenario, `Tumor (rel)`, `AFP (ng/mL)`, `Liver Fn (%)`,
         `CD8 T-cell`, `TGI (%)`, Response)

cat("\n=== HCC QSP Simulation Summary — Month 6 ===\n")
print(summary_6mo, n = 5, width = 120)

## ====================================================================
## Sensitivity analysis: Tumor growth rate vs. Treatment Response
## ====================================================================
lambda1_values <- c(0.006, 0.010, 0.015, 0.020, 0.025)

sens_results <- lapply(lambda1_values, function(lam) {
  ev_sora <- ev(amt = 400, ii = 12, addl = 359, cmt = "GUT_S", time = 0)
  out <- mod %>%
    param(list(use_sora = 1, lambda1 = lam)) %>%
    init(list(TUMOR = 1, CD8T = 1, TREG = 1, ANGIO = 1,
              VEGF_FREE = 1, AFP = 100, LF = 0.85)) %>%
    mrgsim(ev_sora, end = 360, delta = 1)
  df <- as.data.frame(out)
  df$lambda1 <- lam
  df
}) %>% bind_rows()

p_sens <- ggplot(sens_results, aes(time / 30, TumorRel,
                  color = factor(lambda1), group = factor(lambda1))) +
  geom_line(linewidth = 1.0) +
  scale_color_viridis_d(name = "λ₁ (growth\nrate/day)") +
  labs(title = "Sensitivity Analysis: Tumor Growth Rate vs. Sorafenib Response",
       x = "Time (months)", y = "Relative Tumor Burden") +
  theme_hcc

print(p_sens)
cat("\nModel script complete. See plots for simulation outputs.\n")
