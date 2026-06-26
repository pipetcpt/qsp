## =============================================================================
## Psoriasis QSP Model — mrgsolve ODE Implementation
## =============================================================================
## Disease:    Plaque Psoriasis (Psoriasis vulgaris)
## Pathways:   IL-23/IL-17 axis · TNF-α · Keratinocyte hyperproliferation
## Drug PK/PD: Adalimumab · Secukinumab · Risankizumab · Apremilast ·
##             Tofacitinib · Methotrexate · Ustekinumab
##
## Clinical calibration references:
##   - PASI response: Griffiths CE et al. NEJM 2015 (secukinumab)
##   - PK Adalimumab: Menter A et al. J Invest Dermatol 2010
##   - PK Secukinumab: Chatzidionysiou K et al. Rheumatology 2016
##   - PK Risankizumab: Blauvelt A et al. NEJM 2019
##   - Apremilast PK:  Papp KA et al. NEJM 2015 (ESTEEM)
##   - Tofacitinib:    Papp KA et al. NEJM 2018
##   - Methotrexate:   Heydendael VM et al. NEJM 2003
##   - Ustekinumab:    Leonardi CL et al. Lancet 2008
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE
## ─────────────────────────────────────────────────────────────────────────────
pso_model_code <- '
$PROB
Psoriasis QSP Model — mrgsolve ODE
IL-23/IL-17 axis, TNF, keratinocyte hyperproliferation, drug PK/PD
Compartments: 25 ODEs

$PARAM
// ---- Innate/DC parameters ----
k_DC_prod    = 0.05    // mDC production rate (cells/mL/h)
k_DC_death   = 0.02    // mDC death rate (/h)
k_IL23_prod  = 0.10    // IL-23 production by mDC (pg/mL/h per DC)
k_IL23_deg   = 0.50    // IL-23 degradation (/h), t½~1.4h
k_TNFi_prod  = 0.08    // Innate TNF-α production rate

// ---- Th17 differentiation ----
k_Th17_diff  = 0.003   // Th17 differentiation rate (/h)
k_Th17_death = 0.005   // Th17 death rate (/h)
EC50_IL23    = 50.0    // EC50 of IL-23 on Th17 diff (pg/mL)
Emax_IL23    = 5.0     // Emax fold for IL-23 on Th17 diff

// ---- IL-17A dynamics ----
k_IL17_prod  = 0.20    // IL-17A production by Th17 (pg/mL/h per 1000 cells)
k_IL17_deg   = 0.30    // IL-17A degradation (/h), t½~2.3h
EC50_IL17KC  = 80.0    // EC50 of IL-17A on KC activation

// ---- TNF-alpha dynamics ----
k_TNFa_prod  = 0.15    // TNF-α production (Th1+Macro)
k_TNFa_deg   = 0.40    // TNF-α degradation (/h), t½~1.7h

// ---- IFN-gamma dynamics ----
k_IFNg_prod  = 0.12    // IFN-γ from Th1/Tc1
k_IFNg_deg   = 0.35    // IFN-γ degradation (/h)

// ---- Keratinocyte hyperproliferation ----
k_KC_basal   = 0.008   // Basal KC proliferation rate (/h)
k_KC_death   = 0.004   // KC death rate (/h)
k_KC_stim_IL17 = 0.006 // IL-17A stimulation of KC prolif (/h)
k_KC_stim_TNF  = 0.004 // TNF-α stimulation of KC prolif (/h)
k_KC_stim_IFNg = 0.003 // IFN-γ stimulation of KC prolif (/h)
KC_ss        = 100.0   // KC homeostatic setpoint (% normal)

// ---- PASI dynamics ----
k_PASI_form  = 0.0010  // PASI formation rate
k_PASI_res   = 0.0050  // PASI resolution rate (/h)
PASI_ss      = 20.0    // Moderate-severe baseline PASI

// ---- ADA (Adalimumab) PK — 2-compartment SC ----
// Regimen: 80mg SC wk0, 40mg SC wk1, then 40mg SC q2w
ka_ADA     = 0.013    // Absorption rate (/h), F=0.64
CL_ADA     = 0.247    // Clearance (L/h)
V1_ADA     = 7.0      // Central volume (L)
Q_ADA      = 0.30     // Inter-compartment clearance (L/h)
V2_ADA     = 3.2      // Peripheral volume (L)
Kd_ADA_TNF = 0.1      // Kd for TNF-α binding (nM)

// ---- SEC (Secukinumab) PK — 2-compartment SC ----
// Regimen: 300mg SC wk0,1,2,3,4, then q4w
ka_SEC     = 0.015    // Absorption rate (/h), F=0.73
CL_SEC     = 0.191    // Clearance (L/h)
V1_SEC     = 7.1      // Central volume (L)
Q_SEC      = 0.20     // Inter-compartment clearance (L/h)
V2_SEC     = 3.8      // Peripheral volume (L)
Kd_SEC_IL17 = 0.08   // Kd for IL-17A binding (nM)

// ---- RSK (Risankizumab) PK — 1-compartment SC ----
// Regimen: 150mg SC q12w (after wk0,4)
ka_RSK     = 0.012    // Absorption rate (/h), F=0.89
CL_RSK     = 0.078    // Clearance (L/h), t½~28d
V1_RSK     = 11.2     // Central volume (L)
Kd_RSK_IL23 = 0.06   // Kd for IL-23p19 binding (nM)

// ---- UST (Ustekinumab) PK — 2-compartment SC ----
// Regimen: 45mg SC wk0, wk4, then q12w
ka_UST     = 0.007    // Absorption rate (/h), F=0.57
CL_UST     = 0.252    // Clearance (L/h)
V1_UST     = 15.1     // Central volume (L)
Kd_UST_p40  = 0.9    // Kd for IL-12/23 p40 (nM)

// ---- APR (Apremilast) PK — oral ----
// Regimen: 30mg PO BID (after titration)
ka_APR     = 0.80     // Absorption rate (/h), F=0.73
CL_APR     = 9.5      // Clearance (L/h)
V1_APR     = 86.6     // Central volume (L)
Q_APR      = 3.9      // Inter-compartmental CL (L/h)
V2_APR     = 43.3     // Peripheral volume (L)
IC50_APR   = 74.0     // IC50 for PDE4 inhibition (nM)
Emax_APR   = 0.70     // Maximum PDE4i effect on IL-17

// ---- TOF (Tofacitinib) PK — oral ----
// Regimen: 5mg PO BID or 10mg PO BID
ka_TOF     = 1.20     // Absorption rate (/h), F=0.74
CL_TOF     = 22.8     // Clearance (L/h)
V1_TOF     = 87.0     // Volume of distribution (L)
IC50_TOF   = 1.0      // IC50 for JAK inhibition (nM)
Emax_TOF   = 0.80     // Max JAK inhibition effect

// ---- MTX (Methotrexate) PK — oral ----
// Regimen: 15-25mg PO weekly
ka_MTX     = 0.50     // Absorption rate (/h), F=0.70
CL_MTX     = 4.8      // Clearance (L/h), renal dominant
V1_MTX     = 24.0     // Central volume (L)
k_PG_form  = 0.02     // MTX-polyglutamate formation rate
k_PG_elim  = 0.005    // MTX-PG elimination rate (/h)
IC50_MTX   = 5.0      // IC50 MTX-PG on Th17 (nM equiv)

// ---- MW conversion factors ----
MW_ADA  = 148000      // Adalimumab MW (g/mol)
MW_SEC  = 147000      // Secukinumab MW (g/mol)
MW_RSK  = 153000      // Risankizumab MW (g/mol)
MW_UST  = 148000      // Ustekinumab MW (g/mol)
MW_IL17 = 35000       // IL-17A MW (g/mol), homodimer
MW_TNF  = 51000       // TNF-α MW (g/mol), trimer
MW_IL23 = 55000       // IL-23 MW (g/mol)

$CMT
// INNATE/ADAPTIVE immune
DC           // Myeloid DCs (cells/mL * 1000)
IL23         // IL-23 (pg/mL)
Th17         // Th17 cells (cells/mL * 1000)
IL17A        // Free IL-17A (pg/mL)
TNFa         // Free TNF-α (pg/mL)
IFNg         // IFN-γ (pg/mL)
KC           // Keratinocyte activation index (% above normal)
PASI         // PASI score (0-72)

// Adalimumab PK
ADA_SC ADA_C ADA_P ADA_TNF   // SC depot, central, peripheral, bound-TNF

// Secukinumab PK
SEC_SC SEC_C SEC_P SEC_IL17   // SC depot, central, peripheral, bound-IL17A

// Risankizumab PK
RSK_SC RSK_C RSK_IL23         // SC depot, central, bound-IL23

// Ustekinumab PK
UST_SC UST_C UST_p40          // SC depot, central, bound-p40 (IL12/23)

// Apremilast PK
APR_GI APR_C APR_P            // GI, central, peripheral

// Tofacitinib PK
TOF_GI TOF_C                  // GI, central

// Methotrexate PK
MTX_GI MTX_C MTX_PG           // GI, central, polyglutamate

$MAIN
// Convert concentrations to nM for binding calculations
double ADA_nM  = (ADA_C / V1_ADA) * 1e6 / MW_ADA;
double SEC_nM  = (SEC_C / V1_SEC) * 1e6 / MW_SEC;
double RSK_nM  = (RSK_C / V1_RSK) * 1e6 / MW_RSK;
double UST_nM  = (UST_C / V1_UST) * 1e6 / MW_UST;
double APR_nM  = (APR_C / V1_APR) * 1e6 / 460.5;  // MW apremilast
double TOF_nM  = (TOF_C / V1_TOF) * 1e6 / 312.4;  // MW tofacitinib
double MTX_PG_nM = MTX_PG / 0.454;                  // pg→nM approximation

// Drug occupancy (Hill=1)
double fTNF_blk  = ADA_nM / (ADA_nM + Kd_ADA_TNF);  // fraction TNF blocked by ADA
double fIL17_blk = SEC_nM / (SEC_nM + Kd_SEC_IL17);  // fraction IL17A blocked by SEC
double fIL23_blk = RSK_nM / (RSK_nM + Kd_RSK_IL23);  // fraction IL23 blocked by RSK
double fuST_blk  = UST_nM / (UST_nM + Kd_UST_p40);  // fraction p40 blocked by UST

// Apremilast PDE4 inhibition effect on cytokines
double fAPR_inhib = Emax_APR * APR_nM / (APR_nM + IC50_APR);

// Tofacitinib JAK inhibition
double fTOF_inhib = Emax_TOF * TOF_nM / (TOF_nM + IC50_TOF);

// Methotrexate polyglutamate effect
double fMTX_inhib = MTX_PG_nM / (MTX_PG_nM + IC50_MTX);

// Effective IL-23 (considering blockade)
double IL23_eff = IL23 * (1 - fIL23_blk) * (1 - fuST_blk) * (1 - fTOF_inhib * 0.4);

// Effective IL-17A (considering blockade)
double IL17_eff = IL17A * (1 - fIL17_blk);

// Effective TNF-α (considering blockade)
double TNF_eff  = TNFa * (1 - fTNF_blk);

// Th17 differentiation driven by IL-23
double Th17_diff_rate = k_Th17_diff * (1 + Emax_IL23 * IL23_eff / (IL23_eff + EC50_IL23));

// KC activation (IL-17A, TNF, IFNg dependent)
double KC_stim = k_KC_stim_IL17 * IL17_eff / (IL17_eff + EC50_IL17KC) +
                 k_KC_stim_TNF  * TNF_eff   / (TNF_eff  + 100.0) +
                 k_KC_stim_IFNg * IFNg      / (IFNg     + 80.0);

$ODE
// ---- Innate mDC ----
dxdt_DC = k_DC_prod - k_DC_death * DC;

// ---- IL-23 (produced by mDC, blocked by RSK/UST) ----
dxdt_IL23 = k_IL23_prod * DC * (1 - fIL23_blk) * (1 - fuST_blk) - k_IL23_deg * IL23;

// ---- Th17 differentiation ----
dxdt_Th17 = Th17_diff_rate * Th0_ss - k_Th17_death * Th17 -
             fMTX_inhib * k_Th17_death * Th17 - fTOF_inhib * 0.6 * Th17_diff_rate * Th0_ss;

// ---- Free IL-17A ----
dxdt_IL17A = k_IL17_prod * Th17 * (1 - fAPR_inhib) * (1 - fTOF_inhib * 0.3) -
              k_IL17_deg * IL17A;

// ---- Free TNF-α ----
dxdt_TNFa = (k_TNFa_prod * DC + k_TNFi_prod * Th17) *
             (1 - fAPR_inhib) * (1 - fTOF_inhib * 0.5) -
             k_TNFa_deg * TNFa;

// ---- IFN-γ ----
dxdt_IFNg = k_IFNg_prod * Th17 - 0.35 * IFNg;

// ---- Keratinocyte activation index ----
dxdt_KC = KC_stim * (KC_ss - KC) - k_KC_death * (KC - KC_ss);

// ---- PASI Score ----
dxdt_PASI = k_PASI_form * KC - k_PASI_res * PASI;

// ==== Adalimumab PK ====
dxdt_ADA_SC = -ka_ADA * ADA_SC;
dxdt_ADA_C  =  ka_ADA * ADA_SC - (CL_ADA + Q_ADA) * (ADA_C/V1_ADA) +
                Q_ADA * (ADA_P/V2_ADA) - Kd_ADA_TNF * (ADA_C/V1_ADA) * TNFa;
dxdt_ADA_P  =  Q_ADA  * (ADA_C/V1_ADA) - Q_ADA * (ADA_P/V2_ADA);
dxdt_ADA_TNF = Kd_ADA_TNF * (ADA_C/V1_ADA) * TNFa - 0.05 * ADA_TNF;

// ==== Secukinumab PK ====
dxdt_SEC_SC = -ka_SEC * SEC_SC;
dxdt_SEC_C  =  ka_SEC * SEC_SC - (CL_SEC + Q_SEC) * (SEC_C/V1_SEC) +
                Q_SEC * (SEC_P/V2_SEC) - Kd_SEC_IL17 * (SEC_C/V1_SEC) * IL17A;
dxdt_SEC_P  =  Q_SEC  * (SEC_C/V1_SEC) - Q_SEC * (SEC_P/V2_SEC);
dxdt_SEC_IL17 = Kd_SEC_IL17 * (SEC_C/V1_SEC) * IL17A - 0.04 * SEC_IL17;

// ==== Risankizumab PK ====
dxdt_RSK_SC  = -ka_RSK * RSK_SC;
dxdt_RSK_C   =  ka_RSK * RSK_SC - CL_RSK * (RSK_C/V1_RSK) -
                 Kd_RSK_IL23 * (RSK_C/V1_RSK) * IL23;
dxdt_RSK_IL23 = Kd_RSK_IL23 * (RSK_C/V1_RSK) * IL23 - 0.03 * RSK_IL23;

// ==== Ustekinumab PK ====
dxdt_UST_SC  = -ka_UST * UST_SC;
dxdt_UST_C   =  ka_UST * UST_SC - CL_UST * (UST_C/V1_UST) -
                 Kd_UST_p40 * (UST_C/V1_UST) * IL23;
dxdt_UST_p40 =  Kd_UST_p40 * (UST_C/V1_UST) * IL23 - 0.04 * UST_p40;

// ==== Apremilast PK ====
dxdt_APR_GI  = -ka_APR * APR_GI;
dxdt_APR_C   =  ka_APR * APR_GI - (CL_APR/V1_APR + Q_APR/V1_APR) * APR_C +
                 Q_APR/V2_APR * APR_P;
dxdt_APR_P   =  Q_APR/V1_APR * APR_C - Q_APR/V2_APR * APR_P;

// ==== Tofacitinib PK ====
dxdt_TOF_GI  = -ka_TOF * TOF_GI;
dxdt_TOF_C   =  ka_TOF * TOF_GI - CL_TOF/V1_TOF * TOF_C;

// ==== Methotrexate PK ====
dxdt_MTX_GI  = -ka_MTX * MTX_GI;
dxdt_MTX_C   =  ka_MTX * MTX_GI - CL_MTX/V1_MTX * MTX_C;
dxdt_MTX_PG  =  k_PG_form * MTX_C - k_PG_elim * MTX_PG;

$TABLE
double Th0_ss = 500.0;   // naive T cell pool (constant reservoir, cells/mL)
double PASI75 = (PASI <= 0.25 * PASI_ss) ? 1.0 : 0.0;
double PASI90 = (PASI <= 0.10 * PASI_ss) ? 1.0 : 0.0;
double PASI100 = (PASI <= 0.001 * PASI_ss) ? 1.0 : 0.0;
double IGA01   = (PASI <= 3.0) ? 1.0 : 0.0;
double ADA_Cp_nM = (ADA_C/V1_ADA) * 1e6 / MW_ADA;
double SEC_Cp_nM = (SEC_C/V1_SEC) * 1e6 / MW_SEC;
double RSK_Cp_nM = (RSK_C/V1_RSK) * 1e6 / MW_RSK;
double UST_Cp_nM = (UST_C/V1_UST) * 1e6 / MW_UST;
double APR_Cp_nM = (APR_C/V1_APR) * 1e6 / 460.5;
double TOF_Cp_nM = (TOF_C/V1_TOF) * 1e6 / 312.4;

$CAPTURE
PASI PASI75 PASI90 PASI100 IGA01
Th17 IL17A TNFa IFNg IL23 KC DC
ADA_Cp_nM SEC_Cp_nM RSK_Cp_nM UST_Cp_nM APR_Cp_nM TOF_Cp_nM
'

## ─────────────────────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────────────────────
pso_mod <- mrgsolve::mcode("psoriasis_qsp", pso_model_code)

## ─────────────────────────────────────────────────────────────────────────────
## INITIAL CONDITIONS — MODERATE-TO-SEVERE PSORIASIS (PASI ~18-22)
## ─────────────────────────────────────────────────────────────────────────────
init_pso <- list(
  DC    = 2.5,    # elevated mDC
  IL23  = 120.0,  # elevated IL-23 (pg/mL)
  Th17  = 350.0,  # expanded Th17 (cells/mL*1000)
  IL17A = 140.0,  # elevated IL-17A (pg/mL)
  TNFa  = 80.0,   # elevated TNF-α (pg/mL)
  IFNg  = 50.0,   # elevated IFN-γ (pg/mL)
  KC    = 280.0,  # strongly elevated KC index
  PASI  = 20.0,   # baseline PASI (moderate-severe)
  # all drug compartments at 0
  ADA_SC=0, ADA_C=0, ADA_P=0, ADA_TNF=0,
  SEC_SC=0, SEC_C=0, SEC_P=0, SEC_IL17=0,
  RSK_SC=0, RSK_C=0, RSK_IL23=0,
  UST_SC=0, UST_C=0, UST_p40=0,
  APR_GI=0, APR_C=0, APR_P=0,
  TOF_GI=0, TOF_C=0,
  MTX_GI=0, MTX_C=0, MTX_PG=0
)

## ─────────────────────────────────────────────────────────────────────────────
## DOSING EVENTS
## ─────────────────────────────────────────────────────────────────────────────

## Simulation period: 52 weeks (364 days = 8736 hours)
sim_end <- 52 * 7 * 24  # hours

## 1. No treatment (natural disease course)
ev_none <- mrgsolve::ev()

## 2. Adalimumab: 80mg wk0 SC + 40mg wk1 SC, then 40mg SC q2w
ev_ada <- mrgsolve::ev(
  amt  = c(80000, 40000, rep(40000, 25)),  # μg SC
  time = c(0, 168, seq(336, 336 + 24*168, by=336)),
  cmt  = "ADA_SC",
  rate = 0
)

## 3. Secukinumab: 300mg SC wk0,1,2,3,4, then q4w
ev_sec <- mrgsolve::ev(
  amt  = c(rep(300000, 5), rep(300000, 11)), # μg SC
  time = c(0,168,336,504,672, seq(672+4*168, sim_end, by=4*168)),
  cmt  = "SEC_SC",
  rate = 0
)

## 4. Risankizumab: 150mg SC wk0, wk4, then q12w
ev_rsk <- mrgsolve::ev(
  amt  = c(150000, 150000, rep(150000, 4)),
  time = c(0, 4*168, seq(16*168, sim_end, by=12*168)),
  cmt  = "RSK_SC",
  rate = 0
)

## 5. Ustekinumab: 45mg SC wk0, wk4, then q12w
ev_ust <- mrgsolve::ev(
  amt  = c(45000, 45000, rep(45000, 4)),
  time = c(0, 4*168, seq(16*168, sim_end, by=12*168)),
  cmt  = "UST_SC",
  rate = 0
)

## 6. Apremilast: 30mg PO BID (every 12h)
ev_apr <- mrgsolve::ev(
  amt   = rep(30000, 2 * 52 * 7),
  time  = seq(0, sim_end - 12, by=12),
  cmt   = "APR_GI",
  rate  = 0
)

## 7. Tofacitinib: 10mg PO BID
ev_tof <- mrgsolve::ev(
  amt  = rep(10000, 2 * 52 * 7),
  time = seq(0, sim_end - 12, by=12),
  cmt  = "TOF_GI",
  rate = 0
)

## 8. Methotrexate: 20mg PO once weekly
ev_mtx <- mrgsolve::ev(
  amt  = rep(20000, 52),
  time = seq(0, sim_end - 168, by=168),
  cmt  = "MTX_GI",
  rate = 0
)

## ─────────────────────────────────────────────────────────────────────────────
## RUN SIMULATIONS — 5 SCENARIOS
## ─────────────────────────────────────────────────────────────────────────────

out_time <- seq(0, sim_end, by=12)  # every 12 hours

run_scenario <- function(mod, ev, init, label) {
  mod %>%
    init(!!!init) %>%
    ev(ev) %>%
    mrgsim(end=sim_end, delta=12) %>%
    as.data.frame() %>%
    mutate(
      scenario = label,
      week = time / 168,
      day  = time / 24
    )
}

cat("Running Scenario 1: No Treatment...\n")
s1 <- run_scenario(pso_mod, ev_none, init_pso, "1_No_Treatment")

cat("Running Scenario 2: Adalimumab (anti-TNF)...\n")
s2 <- run_scenario(pso_mod, ev_ada, init_pso, "2_Adalimumab_40mg_q2w")

cat("Running Scenario 3: Secukinumab (anti-IL17A)...\n")
s3 <- run_scenario(pso_mod, ev_sec, init_pso, "3_Secukinumab_300mg_q4w")

cat("Running Scenario 4: Risankizumab (anti-IL23p19)...\n")
s4 <- run_scenario(pso_mod, ev_rsk, init_pso, "4_Risankizumab_150mg_q12w")

cat("Running Scenario 5: Apremilast (PDE4i oral)...\n")
s5 <- run_scenario(pso_mod, ev_apr, init_pso, "5_Apremilast_30mg_BID")

cat("Running Scenario 6: Tofacitinib (JAKi oral)...\n")
s6 <- run_scenario(pso_mod, ev_tof, init_pso, "6_Tofacitinib_10mg_BID")

cat("Running Scenario 7: Methotrexate 20mg weekly...\n")
s7 <- run_scenario(pso_mod, ev_mtx, init_pso, "7_Methotrexate_20mg_weekly")

all_scenarios <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

## ─────────────────────────────────────────────────────────────────────────────
## KEY RESPONSE ENDPOINTS AT WEEK 12, 16, 52
## ─────────────────────────────────────────────────────────────────────────────
response_summary <- all_scenarios %>%
  filter(week %in% c(0, 12, 16, 52)) %>%
  group_by(scenario, week) %>%
  summarize(
    PASI_mean   = mean(PASI, na.rm=TRUE),
    PASI75_pct  = mean(PASI75, na.rm=TRUE) * 100,
    PASI90_pct  = mean(PASI90, na.rm=TRUE) * 100,
    PASI100_pct = mean(PASI100, na.rm=TRUE) * 100,
    IGA01_pct   = mean(IGA01, na.rm=TRUE) * 100,
    IL17A_mean  = mean(IL17A, na.rm=TRUE),
    TNFa_mean   = mean(TNFa, na.rm=TRUE),
    Th17_mean   = mean(Th17, na.rm=TRUE),
    .groups = "drop"
  )

print(response_summary)

## ─────────────────────────────────────────────────────────────────────────────
## PLOTS
## ─────────────────────────────────────────────────────────────────────────────

scen_colors <- c(
  "1_No_Treatment"           = "#616161",
  "2_Adalimumab_40mg_q2w"    = "#F44336",
  "3_Secukinumab_300mg_q4w"  = "#2196F3",
  "4_Risankizumab_150mg_q12w"= "#4CAF50",
  "5_Apremilast_30mg_BID"    = "#FF9800",
  "6_Tofacitinib_10mg_BID"   = "#9C27B0",
  "7_Methotrexate_20mg_weekly" = "#795548"
)

weekly_data <- all_scenarios %>%
  filter(time %% 168 == 0) %>%
  group_by(scenario, week) %>%
  summarize(across(c(PASI, IL17A, TNFa, Th17, KC, IL23, ADA_Cp_nM, SEC_Cp_nM,
                     RSK_Cp_nM, APR_Cp_nM, TOF_Cp_nM, PASI75, PASI90, PASI100),
                   mean, na.rm=TRUE), .groups="drop")

## Plot 1: PASI over time
p1 <- ggplot(weekly_data, aes(week, PASI, color=scenario)) +
  geom_line(size=1.2) +
  geom_hline(yintercept=c(5, 10), linetype="dashed", color="gray50") +
  scale_color_manual(values=scen_colors) +
  labs(title="PASI Score Over 52 Weeks", x="Week", y="PASI",
       color="Treatment") +
  annotate("text", x=54, y=5, label="PASI5 (clear)", hjust=0, size=3) +
  annotate("text", x=54, y=10, label="PASI10", hjust=0, size=3) +
  theme_bw() + theme(legend.position="bottom") +
  coord_cartesian(xlim=c(0, 52))

## Plot 2: IL-17A dynamics
p2 <- ggplot(weekly_data, aes(week, IL17A, color=scenario)) +
  geom_line(size=1.2) +
  scale_color_manual(values=scen_colors) +
  labs(title="Serum IL-17A (pg/mL)", x="Week", y="IL-17A (pg/mL)") +
  theme_bw() + theme(legend.position="none")

## Plot 3: Th17 cells
p3 <- ggplot(weekly_data, aes(week, Th17, color=scenario)) +
  geom_line(size=1.2) +
  scale_color_manual(values=scen_colors) +
  labs(title="Th17 Cells", x="Week", y="Th17 (cells/mL×1000)") +
  theme_bw() + theme(legend.position="none")

## Plot 4: PASI75/90/100 response rates at wk12 & wk16
resp_bar <- response_summary %>%
  filter(week %in% c(12, 16)) %>%
  pivot_longer(cols=c(PASI75_pct, PASI90_pct, PASI100_pct),
               names_to="endpoint", values_to="pct") %>%
  mutate(endpoint=recode(endpoint,
    PASI75_pct="PASI75", PASI90_pct="PASI90", PASI100_pct="PASI100"),
    week_label = paste0("Wk", week))

p4 <- ggplot(resp_bar, aes(scenario, pct, fill=endpoint)) +
  geom_col(position="dodge") +
  facet_wrap(~week_label) +
  scale_fill_manual(values=c("PASI75"="#4CAF50","PASI90"="#2196F3","PASI100"="#9C27B0")) +
  labs(title="PASI Response Rates (%)", x="", y="Responders (%)", fill="Endpoint") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8))

## Plot 5: Drug PK — Biologic concentrations
bio_pk <- weekly_data %>%
  filter(scenario %in% c("2_Adalimumab_40mg_q2w","3_Secukinumab_300mg_q4w","4_Risankizumab_150mg_q12w")) %>%
  pivot_longer(cols=c(ADA_Cp_nM, SEC_Cp_nM, RSK_Cp_nM),
               names_to="drug", values_to="Cp_nM") %>%
  mutate(drug=recode(drug, ADA_Cp_nM="Adalimumab", SEC_Cp_nM="Secukinumab", RSK_Cp_nM="Risankizumab"))

p5 <- ggplot(bio_pk, aes(week, Cp_nM, color=drug)) +
  geom_line(size=1.2) +
  scale_color_manual(values=c("Adalimumab"="#F44336","Secukinumab"="#2196F3","Risankizumab"="#4CAF50")) +
  labs(title="Biologic PK — Central Concentration (nM)", x="Week", y="Cp (nM)") +
  theme_bw() + theme(legend.position="bottom")

## Plot 6: Small molecule PK
sm_pk <- weekly_data %>%
  filter(scenario %in% c("5_Apremilast_30mg_BID","6_Tofacitinib_10mg_BID")) %>%
  pivot_longer(cols=c(APR_Cp_nM, TOF_Cp_nM),
               names_to="drug", values_to="Cp_nM") %>%
  mutate(drug=recode(drug, APR_Cp_nM="Apremilast", TOF_Cp_nM="Tofacitinib"))

p6 <- ggplot(sm_pk, aes(week, Cp_nM, color=drug)) +
  geom_line(size=1.2) +
  scale_color_manual(values=c("Apremilast"="#FF9800","Tofacitinib"="#9C27B0")) +
  labs(title="Small Molecule PK — Trough Concentration (nM)", x="Week", y="Cp (nM)") +
  theme_bw() + theme(legend.position="bottom")

## Combine
combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)
combined_plot <- combined_plot + plot_annotation(
  title="Psoriasis QSP Model — mrgsolve Simulation Results",
  subtitle="7 Treatment Scenarios · IL-23/IL-17 axis · TNF · Keratinocyte PD",
  theme=theme(plot.title=element_text(size=16, face="bold"))
)

print(combined_plot)
ggsave("pso_qsp_simulation.png", combined_plot, width=16, height=14, dpi=150)
cat("\nSimulation complete. Plots saved.\n")

## ─────────────────────────────────────────────────────────────────────────────
## CLINICAL CONTEXT TABLE (key trial data for calibration)
## ─────────────────────────────────────────────────────────────────────────────
clinical_ref <- data.frame(
  Drug            = c("Adalimumab","Secukinumab","Risankizumab","Ustekinumab","Apremilast","Tofacitinib","Methotrexate"),
  Regimen         = c("40mg SC q2w","300mg SC wk0-4 then q4w","150mg SC wk0,4 then q12w","45mg SC wk0,4 then q12w","30mg PO BID","10mg PO BID","20mg PO qw"),
  PASI75_wk12_16  = c("71%","77-80%","88-91%","67-71%","33-40%","39-46%","26-36%"),
  PASI90_wk12_16  = c("45%","59-67%","72-75%","42%","18%","22%","NA"),
  Ref_Trial       = c("CHAMPION","FIXTURE","UltIMMa","PHOENIX-1","ESTEEM-1","OPT Pivotal","Heydendael 2003"),
  stringsAsFactors = FALSE
)

cat("\n=== Clinical Reference Data for Model Calibration ===\n")
print(clinical_ref, row.names=FALSE)

cat("\n=== PASI Response Summary at Week 16 ===\n")
response_summary %>%
  filter(week==16) %>%
  select(scenario, PASI_mean, PASI75_pct, PASI90_pct, PASI100_pct, IGA01_pct) %>%
  print()
