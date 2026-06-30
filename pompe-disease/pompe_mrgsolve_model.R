## =================================================================
## Pompe Disease (GSDII) — Quantitative Systems Pharmacology
## mrgsolve ODE model (22 compartments)
##
## Disease: GAA enzyme deficiency → lysosomal glycogen accumulation →
##          autophagic buildup → cardiomyopathy (IOPD) + skeletal myopathy
##          + diaphragmatic / respiratory failure (LOPD)
##
## Drugs implemented (PK + PD coupling):
##   • Alglucosidase alfa (Lumizyme/Myozyme) 20 mg/kg IV q2w
##   • Avalglucosidase alfa (Nexviazyme, COMET trial) 20 mg/kg IV q2w
##   • Cipaglucosidase alfa + Miglustat (Pombiliti+Opfolda, PROPEL)
##   • AAV9-hGAA gene-therapy bolus (research stage)
##   • Rituximab-based ITI (CRIM- IOPD prophylaxis)
##
## Endpoints:
##   • LV mass index (g/m^2)  – IOPD primary
##   • FVC upright (% pred)   – LOPD primary (COMET)
##   • 6-minute walk distance (m)
##   • GMFM-88, ventilator-free survival, anti-GAA ADA titre, Hex4
##
## Parameter values approximate adult (70 kg) LOPD unless noted; IOPD
## variant flagged by IOPD_FLAG = 1. All values illustrative; see
## `pompe_references.md` for source ranges.
## =================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

pompe_code <- '
$PROB
# Pompe Disease QSP v1.0
# 22-compartment ODE; adult LOPD baseline, IOPD switchable.

$PARAM @annotated
// ---- Patient & disease ----
WT        :  70    : Body weight (kg)
BSA       :  1.80  : Body surface area (m^2)
IOPD_FLAG :   0    : 1 = infantile-onset, 0 = LOPD adult
CRIM_NEG  :   0    : 1 = CRIM-negative (higher ADA risk)
GAA_BASE  :   0.10 : Residual GAA activity (fraction of normal)
SEX_M     :   1    : 1 = male, 0 = female (covariate placeholder)

// ---- Alglucosidase alfa PK (2-cmt) ----
KA_ALGLU  :   0.0  : Not used (IV infusion)
CL_ALGLU  :  21    : Clearance L/d (~0.27 mL/min/kg * 70 * 1440)
V1_ALGLU  :   3.5  : Central volume (L) ~50 mL/kg
Q_ALGLU   :   2.0  : Intercompartmental clearance (L/d)
V2_ALGLU  :   4.5  : Peripheral volume (L)
M6P_ALGLU :   2.5  : Mannose-6-P per mole alglucosidase (mol/mol)

// ---- Avalglucosidase alfa PK ----
CL_AVAL   :  17    : Clearance (L/d) lower than alglu
V1_AVAL   :   3.3  : Central V (L)
Q_AVAL    :   2.0
V2_AVAL   :   4.5
M6P_AVAL  :  37.5  : ~15-fold higher M6P content

// ---- Cipaglucosidase alfa PK ----
CL_CIPA   :  22
V1_CIPA   :   3.6
Q_CIPA    :   1.8
V2_CIPA   :   4.5
M6P_CIPA  :  25

// ---- Miglustat (stabiliser) ----
KA_MIG    :   3.5  : Miglustat absorption (1/d)
CL_MIG    :  84    : Clearance (L/d)
V_MIG     : 100    : Volume (L)
F_MIG     :   0.83 : Bioavailability
MIG_STAB  :   0.35 : Plasma stabilisation factor (Imax)
MIG_EC50  :   1.5  : Plasma stabilisation EC50 (mg/L)

// ---- M6P-receptor uptake (Michaelis–Menten) ----
VMAX_UPT  :   8.0  : Maximal tissue uptake rate (mg/d)
KM_UPT    :   0.7  : Km (mg/L) for CI-MPR
RHO_AVAL  :   3.0  : Relative uptake potency Aval vs Alglu (M6P x15 → uptake gain x3)
RHO_CIPA  :   1.8  : Relative uptake potency Cipa vs Alglu
ADA_KI    :  10    : ADA neutralising IC50 (titre units)
ADA_BLOCK :   0.95 : Maximal ADA-mediated block fraction

// ---- Tissue lysosomal pool ----
KOUT_TISS :   0.07 : Lysosomal enzyme degradation (1/d)
GAA_TURN  :   0.05 : Endogenous GAA recovery rate (1/d)
GLYC_SS_M :  10    : Steady-state lysosomal glycogen muscle (a.u.)
GLYC_SS_C :  15    : Steady-state cardiac glycogen (a.u., IOPD ref)
KIN_GLYC  :   0.20 : Glycogen accumulation rate (1/d)
KOUT_GLYC :   0.04 : Baseline glycogen clearance (1/d)
KGAA_GLYC :   1.20 : GAA-driven glycogen hydrolysis (1/d / nmol)
HX4_GAIN  :   0.30 : Hex4 generation gain
KOUT_HX4  :   0.6  : Urinary Hex4 clearance (1/d)

// ---- ADA dynamics ----
KP_AB     :   0.012 : Plasmablast generation per dose (1/d)
KP_AB_CRIM:   0.045 : CRIM- plasmablast surge (1/d)
KOUT_AB   :   0.05  : ADA decay (1/d, ~14 d t1/2)
ADA_AMP   :   1.0   : Amplification factor on antigen exposure

// ---- Disease physiology ----
MM_BASE   :  20.0   : Muscle mass index baseline (kg)
KM_LOSS   :   0.0012 : Muscle loss rate (1/d) per a.u. glycogen
KM_GAIN   :   0.0006 : Muscle regeneration rate per ERT delivery (1/d)
DIAPH_BASE:   1.0    : Diaphragm function (1 = normal)
DIAPH_LOSS:   0.0010
DIAPH_GAIN:   0.0006
FVC_MAX   :  95      : FVC upright max (% predicted)
FVC_MIN   :  20      : Floor
LVMI_BASE :  60      : LV mass index baseline (g/m^2, adult)
LVMI_IOPD : 250      : IOPD presenting LVMI
KLV_GAIN  :   0.05   : LV hypertrophy rate per cardiac glycogen
KLV_REVERSE:  0.04   : Reverse remodel rate per delivered enzyme
SMWT_MAX  : 600
SMWT_MIN  : 100

// ---- Endpoint thresholds ----
VENT_FVC_THR : 35    : FVC < 35% triggers vent_failure hazard
VENT_HAZARD  :  0.0008 : Daily hazard for vent_failure if FVC < threshold

// ---- AAV gene therapy ----
AAV_DOSE  :   0      : Single bolus (vg/kg, set externally)
AAV_kexp  :   0.0035 : Expression rise (1/d)
AAV_DECAY :   0.00010: Vector dilution / immune loss (1/d)
AAV_GAIN  :   0.40   : Max contribution to tissue GAA (fraction of normal)
AAV_NAB   :   1.0    : Pre-existing AAV NAb modifier (0–1, 0 blocks)

// ---- Rituximab ITI (for CRIM-) ----
RTX_KIN   :   0.0    : Set externally per dose
RTX_KOUT  :   0.04   : Rituximab decay (1/d, t1/2 ~ 20 d)
RTX_ADA_K :   0.20   : Maximal inhibition of plasmablast generation
RTX_EC50  :  10      : Rituximab EC50 (mg/L)

$CMT @annotated
ALGLU_C   : alglucosidase central (mg)
ALGLU_P   : alglucosidase peripheral (mg)
AVAL_C    : avalglucosidase central (mg)
AVAL_P    : avalglucosidase peripheral (mg)
CIPA_C    : cipaglucosidase central (mg)
CIPA_P    : cipaglucosidase peripheral (mg)
MIG_A     : miglustat absorption depot (mg)
MIG_C     : miglustat central (mg)
GAA_M     : muscle lysosomal GAA pool (a.u.)
GAA_C     : cardiac lysosomal GAA pool (a.u.)
GAA_D     : diaphragm GAA pool (a.u.)
GLYC_M    : muscle glycogen (a.u.)
GLYC_C    : cardiac glycogen (a.u.)
GLYC_D    : diaphragm glycogen (a.u.)
HEX4      : plasma Hex4 biomarker (a.u.)
ADA_T     : anti-GAA ADA titre (units)
LVMI      : LV mass index (g/m^2)
MM_IDX    : muscle mass index (kg)
DIAPH_F   : diaphragm function (0-1)
FVC_UP    : FVC upright (% predicted)
AAV_X     : AAV-driven tissue GAA expression (a.u.)
RTX_C     : Rituximab central (mg)

$MAIN
// initial conditions tuned for IOPD vs LOPD presentation
double GLYC_init_M  = GLYC_SS_M * (1 - GAA_BASE) + 1e-6;
double GLYC_init_C  = (IOPD_FLAG > 0.5 ? GLYC_SS_C : 2.0) * (1 - GAA_BASE);
double GLYC_init_D  = GLYC_SS_M * (1 - GAA_BASE);
double LVMI_init    = (IOPD_FLAG > 0.5 ? LVMI_IOPD : LVMI_BASE);

GAA_M_0   = GAA_BASE * 10.0;
GAA_C_0   = GAA_BASE * 10.0;
GAA_D_0   = GAA_BASE * 10.0;
GLYC_M_0  = GLYC_init_M;
GLYC_C_0  = GLYC_init_C;
GLYC_D_0  = GLYC_init_D;
HEX4_0    = 8.0 * (1 - GAA_BASE);
ADA_T_0   = 0.0;
LVMI_0    = LVMI_init;
MM_IDX_0  = MM_BASE * (IOPD_FLAG > 0.5 ? 0.4 : 0.85);
DIAPH_F_0 = DIAPH_BASE * (IOPD_FLAG > 0.5 ? 0.30 : 0.80);
FVC_UP_0  = FVC_MAX * (IOPD_FLAG > 0.5 ? 0.35 : 0.65);
AAV_X_0   = 0.0;
RTX_C_0   = 0.0;

$ODE
// ---- Plasma concentrations ----
double Cp_alglu = ALGLU_C / V1_ALGLU;
double Cp_aval  = AVAL_C  / V1_AVAL;
double Cp_cipa  = CIPA_C  / V1_CIPA;
double Cp_mig   = MIG_C   / V_MIG;
double Cp_rtx   = RTX_C   / V1_ALGLU;   // approx

// Miglustat stabilisation (Imax on Cipa CL)
double mig_eff  = MIG_STAB * Cp_mig / (MIG_EC50 + Cp_mig);
double Cp_cipa_stab = Cp_cipa * (1 + mig_eff);

// ADA-mediated neutralisation
double ada_block = ADA_BLOCK * pow(ADA_T,2) / (pow(ADA_KI,2) + pow(ADA_T,2));

// CI-MPR uptake (Michaelis–Menten, sum across drugs)
double uptake_alglu = VMAX_UPT * Cp_alglu / (KM_UPT + Cp_alglu) * (1 - ada_block);
double uptake_aval  = VMAX_UPT * RHO_AVAL * Cp_aval / (KM_UPT/RHO_AVAL + Cp_aval) * (1 - 0.5*ada_block);
double uptake_cipa  = VMAX_UPT * RHO_CIPA * Cp_cipa_stab / (KM_UPT + Cp_cipa_stab) * (1 - 0.5*ada_block);

double tissue_supply = uptake_alglu + uptake_aval + uptake_cipa;
double aav_supply    = AAV_GAIN * AAV_X * AAV_NAB;

// Drug PK ODEs (IV infusions/bolus drive in vivo dosing events)
dxdt_ALGLU_C = - (CL_ALGLU / V1_ALGLU) * ALGLU_C - Q_ALGLU/V1_ALGLU * ALGLU_C + Q_ALGLU/V2_ALGLU * ALGLU_P;
dxdt_ALGLU_P =   Q_ALGLU/V1_ALGLU * ALGLU_C - Q_ALGLU/V2_ALGLU * ALGLU_P;
dxdt_AVAL_C  = - (CL_AVAL / V1_AVAL) * AVAL_C - Q_AVAL/V1_AVAL * AVAL_C + Q_AVAL/V2_AVAL * AVAL_P;
dxdt_AVAL_P  =   Q_AVAL/V1_AVAL * AVAL_C - Q_AVAL/V2_AVAL * AVAL_P;
dxdt_CIPA_C  = - (CL_CIPA / V1_CIPA) * CIPA_C - Q_CIPA/V1_CIPA * CIPA_C + Q_CIPA/V2_CIPA * CIPA_P;
dxdt_CIPA_P  =   Q_CIPA/V1_CIPA * CIPA_C - Q_CIPA/V2_CIPA * CIPA_P;
dxdt_MIG_A   = - KA_MIG * MIG_A;
dxdt_MIG_C   =   KA_MIG * MIG_A * F_MIG - (CL_MIG / V_MIG) * MIG_C;
dxdt_RTX_C   = - RTX_KOUT * RTX_C;

// Tissue lysosomal enzyme pools (muscle, cardiac, diaphragm)
double rec_GAA   = GAA_TURN * (1 - GAA_BASE);  // small endogenous recovery only if GAA_BASE>0
dxdt_GAA_M   = 0.55 * tissue_supply + 0.6 * aav_supply + rec_GAA - KOUT_TISS * GAA_M;
dxdt_GAA_C   = 0.25 * tissue_supply + 0.3 * aav_supply + rec_GAA - KOUT_TISS * GAA_C;
dxdt_GAA_D   = 0.20 * tissue_supply + 0.1 * aav_supply + rec_GAA - KOUT_TISS * GAA_D;

// Lysosomal glycogen dynamics
dxdt_GLYC_M = KIN_GLYC - (KOUT_GLYC + KGAA_GLYC * GAA_M / 10.0) * GLYC_M;
dxdt_GLYC_C = KIN_GLYC * (IOPD_FLAG > 0.5 ? 1.5 : 0.4)
              - (KOUT_GLYC + KGAA_GLYC * GAA_C / 10.0) * GLYC_C;
dxdt_GLYC_D = KIN_GLYC - (KOUT_GLYC + KGAA_GLYC * GAA_D / 10.0) * GLYC_D;

// Hex4 biomarker (sum proxy)
double hex4_gen = HX4_GAIN * (GLYC_M + 0.5*GLYC_D);
dxdt_HEX4 = hex4_gen - KOUT_HX4 * HEX4;

// ADA dynamics (per-dose Ag triggers via mtime/event)
double ag_drive = (Cp_alglu + Cp_aval + Cp_cipa) * ADA_AMP;
double rtx_eff  = RTX_ADA_K * Cp_rtx / (RTX_EC50 + Cp_rtx);
double KP_AB_eff = (CRIM_NEG > 0.5 ? KP_AB_CRIM : KP_AB) * (1 - rtx_eff);
dxdt_ADA_T   = KP_AB_eff * ag_drive - KOUT_AB * ADA_T;

// AAV expression rise then dilution
dxdt_AAV_X   = AAV_kexp * (AAV_DOSE > 0 ? 1.0 - AAV_X : -AAV_X)
               - AAV_DECAY * AAV_X;

// LVMI dynamics (IOPD primary)
dxdt_LVMI = KLV_GAIN * GLYC_C - KLV_REVERSE * (GAA_C + aav_supply);

// Skeletal muscle mass index
dxdt_MM_IDX = -KM_LOSS * GLYC_M * MM_IDX + KM_GAIN * (GAA_M + aav_supply) * (MM_BASE - MM_IDX);

// Diaphragm function (0–1)
dxdt_DIAPH_F = -DIAPH_LOSS * GLYC_D * DIAPH_F + DIAPH_GAIN * (GAA_D + aav_supply) * (1.0 - DIAPH_F);

// FVC upright tracks diaphragm function and muscle mass
double FVC_target = FVC_MIN + (FVC_MAX - FVC_MIN) * (0.7 * DIAPH_F + 0.3 * MM_IDX / MM_BASE);
dxdt_FVC_UP = 0.05 * (FVC_target - FVC_UP);

$TABLE
// 6MWT depends on muscle + diaphragm + FVC
double SMWT = SMWT_MIN + (SMWT_MAX - SMWT_MIN) * (0.5*MM_IDX/MM_BASE + 0.3*DIAPH_F + 0.2*FVC_UP/FVC_MAX);
if (SMWT < SMWT_MIN) SMWT = SMWT_MIN;
if (SMWT > SMWT_MAX) SMWT = SMWT_MAX;

// Composite endpoints
double VENT_RISK = (FVC_UP < VENT_FVC_THR ? VENT_HAZARD * (VENT_FVC_THR - FVC_UP) : 0.0);
double SF36_PCS  = 30 + 25 * (MM_IDX/MM_BASE) + 15 * (FVC_UP/FVC_MAX);
double NTproBNP  = 100 + 25 * LVMI;
double EF_LV     = 65 - 0.05 * fmax(0, LVMI - LVMI_BASE);
double CK        = 200 + 80 * GLYC_M;

capture Cp_alglu;
capture Cp_aval;
capture Cp_cipa;
capture Cp_mig;
capture Cp_rtx;
capture SMWT;
capture VENT_RISK;
capture SF36_PCS;
capture NTproBNP;
capture EF_LV;
capture CK;
capture tissue_supply;
capture ada_block;
'

pompe_mod <- mcode("pompe_qsp", pompe_code)

## ----------------------------------------------------------------
## Convenience scenario runner
## ----------------------------------------------------------------

#' Run a Pompe-disease QSP scenario
#' @param scenario  one of c("no_tx","alglu","aval","cipa_mig","aav_gt","alglu_iti")
#' @param iopd      logical; 1 = infantile, 0 = LOPD
#' @param years     simulation duration (years)
#' @return long-format data.frame with simulated outputs
pompe_run <- function(scenario = "alglu", iopd = FALSE, years = 3) {
  if (!requireNamespace("mrgsolve", quietly = TRUE))
    stop("mrgsolve package required")

  param <- list(IOPD_FLAG = as.integer(iopd),
                CRIM_NEG  = as.integer(iopd) * 0.2,  # ~20% IOPD CRIM- assumption
                AAV_DOSE  = 0)

  ev <- mrgsolve::ev()

  q2w <- function(amt_per_kg, days, drug_cmt) {
    n <- floor(days / 14)
    times <- (0:n) * 14
    mrgsolve::ev(time = times, amt = amt_per_kg * 70, rate = -1, cmt = drug_cmt, evid = 1)
  }

  if (scenario == "alglu")     ev <- q2w(20, years*365, "ALGLU_C")
  if (scenario == "aval")      ev <- q2w(20, years*365, "AVAL_C")
  if (scenario == "cipa_mig") {
    ev_cipa <- q2w(20, years*365, "CIPA_C")
    n   <- years*365
    ev_mig <- mrgsolve::ev(time = seq(0, n, 1), amt = 195, cmt = "MIG_A", evid = 1)
    ev <- c(ev_cipa, ev_mig)
  }
  if (scenario == "aav_gt") {
    param$AAV_DOSE <- 1
    ev <- mrgsolve::ev(time = 0, amt = 0, cmt = "AAV_X", evid = 0)
  }
  if (scenario == "alglu_iti") {
    ev_a  <- q2w(20, years*365, "ALGLU_C")
    ev_rt <- mrgsolve::ev(time = c(-14, -7, 0, 7), amt = 700, cmt = "RTX_C", evid = 1)
    ev    <- c(ev_a, ev_rt)
    param$CRIM_NEG <- 1
  }
  if (scenario == "no_tx") {
    ev <- mrgsolve::ev(time = 0, amt = 0, cmt = "ALGLU_C", evid = 0)
  }

  out <- pompe_mod |>
    mrgsolve::param(param) |>
    mrgsolve::ev(ev) |>
    mrgsolve::mrgsim(end = years * 365, delta = 1) |>
    as.data.frame()

  out$scenario <- scenario
  out
}

## ----------------------------------------------------------------
## Calibration notes (anchors)
## ----------------------------------------------------------------
## • Adult LOPD baseline 6MWT ~ 300 m (Kishnani 2019 NEJM COMET)
## • Alglucosidase alfa CL ~ 0.27 mL/min/kg, V_ss ~ 100 mL/kg (Hahn 2008)
## • Avalglucosidase alfa ~3-fold higher M6P-tagged binding & uptake (Pena 2019)
## • Cipaglucosidase+Miglustat (PROPEL): +14 m 6MWT vs alglucosidase (Schoser Lancet Neurol 2021)
## • IOPD untreated mortality <12 mo, ERT brings 1-yr ventilator-free survival to ~88% (Kishnani 2007)
## • CRIM- IOPD HSAT ≥51,200 → loss of clinical response (Banugaria 2011, Messinger 2012 ITI)
## • Hex4 declines on ERT (Young 2009, An 2005)

## Example: out <- pompe_run("alglu", iopd = FALSE, years = 3)
##          ggplot(out, aes(time, FVC_UP)) + geom_line()
