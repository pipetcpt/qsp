## ============================================================
## NAFLD/NASH QSP Model — mrgsolve Implementation
## Non-Alcoholic Fatty Liver Disease / Steatohepatitis
##
## Model features:
##  - Drug PK: Resmetirom (THRβ agonist), OCA (FXR agonist),
##             Semaglutide (GLP-1 RA), Empagliflozin (SGLT2i)
##  - Disease PD: Hepatic lipid metabolism, inflammation,
##                fibrosis, insulin resistance
##  - ≥18 ODE compartments
##  - ≥5 treatment scenarios
##
## Parameters INFORMED BY (effect sizes qualitatively consistent with; NOT formally fitted):
##  - Resmetirom : MAESTRO-NASH Ph3 (Harrison SA et al., NEJM 2024;390:497-509);
##                 MAESTRO-NAFLD-1 (Harrison SA et al., Nat Med 2023;29:2919-2928)
##  - Obeticholic acid : REGENERATE Ph3 (Younossi ZM et al., Lancet 2019;394:2184-2196;
##                 final analysis Sanyal AJ et al., J Hepatol 2023;79:1110-1120)
##  - Semaglutide : Ph2 (Newsome PN et al., NEJM 2021;384:1113-1124);
##                 ESSENCE Ph3 (Sanyal AJ et al., NEJM 2025;392:2089-2099)
##  - Empagliflozin : E-LIFT (Kuchay MS et al., Diabetes Care 2018;41:1801-1808);
##                 EMPA-REG OUTCOME (Zinman B et al., NEJM 2015) for CV/safety context only
## NOTE: disease-side parameters are phenomenological — chosen for steady-state stability
##       (kin = KOUT*baseline) and loop gain < 1, NOT individually fitted to trial data.
##       See nafld_model_design_brief.md.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ── Model code block ─────────────────────────────────────────
code <- '
$PROB NAFLD/NASH QSP Model

$PARAM
// ── Resmetirom PK (100 mg QD oral) ──────────────────────
KA_RSM    = 0.8      // absorption rate constant (1/h)
CL_RSM    = 4.5      // clearance (L/h)
V1_RSM    = 25       // central Vd (L)
Q_RSM     = 8.0      // intercompartmental clearance (L/h)
V2_RSM    = 60       // peripheral Vd (L)
FTUP_RSM  = 0.70     // hepatic uptake fraction (OATP1B)
EC50_RSM  = 0.08     // EC50 for THRβ activation (µg/L)
EMAX_RSM  = 0.85     // max THRβ effect (fraction)

// ── OCA PK (25 mg QD oral) ───────────────────────────────
KA_OCA    = 1.2      // (1/h)
CL_OCA    = 15       // (L/h) — high hepatic clearance
V1_OCA    = 20       // (L)
EC50_OCA  = 0.05     // EC50 for FXR activation (µg/L)
EMAX_OCA  = 0.75     // max FXR effect

// ── Semaglutide PK (2.4 mg SC weekly) ───────────────────
KA_SEM    = 0.0045   // SC absorption (1/h); slow
CL_SEM    = 0.034    // (L/h)
V1_SEM    = 3.2      // (L)
EC50_SEM  = 0.025    // EC50 for GLP-1R (µg/L)
EMAX_SEM  = 0.80     // max weight/IR reduction

// ── Empagliflozin PK (10 mg QD oral) ────────────────────
KA_EMP    = 1.5      // (1/h)
CL_EMP    = 10       // (L/h)
V1_EMP    = 70       // (L)
EC50_EMP  = 0.015    // EC50 SGLT2 (same conc. units as Cp_EMP) — recal: was 0.15 (10 mg gave E_EMP~0.15, near-inert)
EMAX_EMP  = 0.90     // max glycosuric effect
WEMP_LF   = 0.47     // empagliflozin hepatic-fat efflux gain (NEW; phenomenological surrogate for SGLT2i-
                     // induced hepatic FAO). Calibrated to E-LIFT PLACEBO-CORRECTED (between-group) PDFF
                     // effect ~ -24.7% (empa 16.2->11.3 minus control 16.4->15.5 = -4.0 abs), NOT the raw
                     // within-arm -30%, because the QSP placebo arm is flat by construction. Single 20-wk
                     // n=50 trial anchor -> treat as uncertain; sensitivity-analyse WEMP_LF.

// ── Disease biology parameters ───────────────────────────
// Every disease pool uses a steady-state-consistent turnover (indirect-response)
// form: dxdt_X = KOUT_X*X0*(driver, =1 at baseline) - KOUT_X*X, so the initialized
// baseline IS a steady state (kin = KOUT*baseline) and placebo stays flat
// (Dayneka 1993; Jusko & Ko 1994; Woo 2009). Cross-talk uses dimensionless
// fold-change gains sized so the steatosis→Kupffer→TNFα→IR→DNL→fat loop gain < 1
// (Angeli-Ferrell-Sontag 2004). See nafld_model_design_brief.md.

// Hepatic fat / steatosis (influx split: DNL + adipose-NEFA uptake)
KOUT_LF   = 0.003    // liver fat turnover (1/h; t½~10 d). kin = KOUT_LF*LF0
LF0       = 0.20     // baseline hepatic fat fraction (20%, NASH)
WDNL      = 0.4      // fraction of fat influx from ENDOGENOUS+DIETARY (DNL 26% + diet 15% per Donnelly 2005;
                     // lumped because both are small/IR-linked). NOTE: drug DNL-suppression also acts on the
                     // ~15% dietary part -> modest over-credit; future: split WDNL=0.25 (DNL) + WDIET=0.15 const.
WUPT      = 0.6      // fraction from systemic NEFA uptake (∝ body weight) — Donnelly 2005 NEFA 59%

// De novo lipogenesis
DNL0      = 1.0      // baseline DNL (relative, =1)
KDNL_IR   = 0.4      // IR → DNL fold-change sensitivity

// Insulin resistance (set-point form; already SS-consistent)
IR0       = 2.5      // baseline HOMA-IR
KOUT_IR   = 0.02     // IR resolution rate (1/h)
KFFA_IR   = 0.15     // liver-fat → IR sensitivity
KTNF_IR   = 0.25     // TNF-α → IR sensitivity

// Kupffer cell activation (saturable lipotoxicity drive)
KUP0      = 0.5      // baseline Kupffer activation (0-1)
KOUT_KUP  = 0.05     // (1/h)
GLIP_KUP  = 0.5      // lipotoxicity → Kupffer sensitivity
GSAT_KUP  = 0.5      // saturation of Kupffer drive (ceiling = 1+GLIP_KUP/GSAT_KUP)

// TNF-alpha (Kupffer-driven, adiponectin-protected)
KOUT_TNF  = 0.12     // (1/h)
GKUP_TNF  = 0.5      // Kupffer → TNF-α sensitivity
TNF0      = 0.5      // baseline

// IL-6 (Kupffer-driven)
KOUT_IL6  = 0.15     // (1/h)
GKUP_IL6  = 0.5      // Kupffer → IL-6 sensitivity
IL60      = 0.33     // baseline

// TGF-β1 (Kupffer + lipotoxicity drive)
KOUT_TGF  = 0.08     // (1/h)
GKUP_TGF  = 0.4      // Kupffer → TGF-β sensitivity
GLIP_TGF  = 0.3      // lipotoxicity → TGF-β sensitivity
TGF0      = 0.15     // baseline

// Hepatic Stellate Cell (HSC) activation (TGF-β-driven)
KOUT_HSC  = 0.005    // quiescence rate (1/h) — slow (days)
GTGF_HSC  = 0.6      // TGF-β → HSC sensitivity
HSC0      = 0.10     // baseline (low in healthy)

// Collagen / Fibrosis (slowest pool — months)
KOUT_COL  = 0.0008   // collagen turnover (1/h; t½~36 d)
GHSC_COL  = 0.6      // HSC → collagen sensitivity
COL0      = 0.15     // baseline collagen (rel)
KFIBREG   = 0.8      // fibrosis-collagen conversion (stage per unit)

// Liver enzymes (ALT injury = sum of normalized TNF + lipotoxicity)
KOUT_ALT  = 0.030    // ALT turnover (1/h; t½~23h). kin = KOUT_ALT*ALT0
GTNF_ALT  = 0.5      // TNF-α → ALT injury sensitivity
GLIP_ALT  = 0.5      // lipotoxicity → ALT injury sensitivity
ALT0      = 45       // baseline ALT (U/L) — elevated NASH

// Adiponectin (protective; inversely with obesity)
ADIPON0   = 6.0      // baseline (µg/mL; reduced in NASH)
KOUT_ADI  = 0.05     // (1/h)

// Reference values (steatosis/lipotoxicity scaling)
FXR0      = 0.5      // baseline FXR activation (0-1)
LIPOTOX0  = 0.3      // baseline lipotoxicity index

// Body weight
WT0       = 95       // baseline weight (kg)
KOUT_WT   = 0.0015   // (1/h)

$PARAM
// Dosing flags (0=off, 1=on)
DOSE_RSM  = 0        // Resmetirom 100 mg QD
DOSE_OCA  = 0        // OCA 25 mg QD
DOSE_SEM  = 0        // Semaglutide 2.4 mg SC weekly
DOSE_EMP  = 0        // Empagliflozin 10 mg QD

$CMT
// PK compartments
RSM_GUT RSM_CENT RSM_PERI RSM_LIVER
OCA_GUT OCA_CENT
SEM_SC  SEM_CENT
EMP_GUT EMP_CENT

// Disease state compartments
LIVER_FAT  // hepatic fat fraction
INS_RES    // insulin resistance (HOMA-IR proxy)
KUPFFER    // Kupffer cell activation (0-1)
TNFA       // TNF-α (rel units)
IL6C       // IL-6 (rel units)
TGFB       // TGF-β1 (rel units)
HSC        // hepatic stellate cell activation (0-1)
COLLAGEN   // collagen / ECM content (rel)
ALT_CMT    // serum ALT (U/L)
ADIPONECTIN // adiponectin (µg/mL)
BODY_WT    // body weight (kg)

$MAIN
// ── Steady-state initialization ──────────────────────────
LIVER_FAT_0   = LF0;
INS_RES_0     = IR0;
KUPFFER_0     = KUP0;
TNFA_0        = TNF0;
IL6C_0        = IL60;
TGFB_0        = TGF0;
HSC_0         = HSC0;
COLLAGEN_0    = COL0;
ALT_CMT_0     = ALT0;
ADIPONECTIN_0 = ADIPON0;
BODY_WT_0     = WT0;

$ODE
// ─────────────────────────────────────────────────────────
// Drug effects (Emax models)
// ─────────────────────────────────────────────────────────
double Cp_RSM = RSM_CENT / V1_RSM;  // µg/L central
double Cl_RSM_h = FTUP_RSM * Cp_RSM; // liver conc approx
double E_RSM = EMAX_RSM * Cl_RSM_h / (EC50_RSM + Cl_RSM_h);   // THRβ activation [0,1]

double Cp_OCA = OCA_CENT / V1_OCA;
double E_OCA  = EMAX_OCA * Cp_OCA / (EC50_OCA + Cp_OCA);       // FXR activation

double Cp_SEM = SEM_CENT / V1_SEM;
double E_SEM  = EMAX_SEM * Cp_SEM / (EC50_SEM + Cp_SEM);       // GLP-1R effect

double Cp_EMP = EMP_CENT / V1_EMP;
double E_EMP  = EMAX_EMP * Cp_EMP / (EC50_EMP + Cp_EMP);       // SGLT2 inhibition

// ─────────────────────────────────────────────────────────
// Resmetirom PK ODEs
// ─────────────────────────────────────────────────────────
dxdt_RSM_GUT   = -KA_RSM * RSM_GUT;
dxdt_RSM_CENT  = KA_RSM * RSM_GUT
               - (CL_RSM + Q_RSM) / V1_RSM * RSM_CENT
               + Q_RSM / V2_RSM * RSM_PERI;
dxdt_RSM_PERI  = Q_RSM / V1_RSM * RSM_CENT
               - Q_RSM / V2_RSM * RSM_PERI;
dxdt_RSM_LIVER = FTUP_RSM * CL_RSM / V1_RSM * RSM_CENT
               - 0.5 * RSM_LIVER;  // hepatic metabolism

// ─────────────────────────────────────────────────────────
// OCA PK ODEs
// ─────────────────────────────────────────────────────────
dxdt_OCA_GUT  = -KA_OCA * OCA_GUT;
dxdt_OCA_CENT = KA_OCA * OCA_GUT - CL_OCA / V1_OCA * OCA_CENT;

// ─────────────────────────────────────────────────────────
// Semaglutide PK ODEs
// ─────────────────────────────────────────────────────────
dxdt_SEM_SC   = -KA_SEM * SEM_SC;
dxdt_SEM_CENT = KA_SEM * SEM_SC - CL_SEM / V1_SEM * SEM_CENT;

// ─────────────────────────────────────────────────────────
// Empagliflozin PK ODEs
// ─────────────────────────────────────────────────────────
dxdt_EMP_GUT  = -KA_EMP * EMP_GUT;
dxdt_EMP_CENT = KA_EMP * EMP_GUT - CL_EMP / V1_EMP * EMP_CENT;

// ─────────────────────────────────────────────────────────
// Body weight — reduced by semaglutide and empagliflozin
// ─────────────────────────────────────────────────────────
double WT_target = WT0 * (1 - 0.15 * E_SEM - 0.05 * E_EMP);
dxdt_BODY_WT = KOUT_WT * (WT_target - BODY_WT);

// ─────────────────────────────────────────────────────────
// Adiponectin — increases with weight loss
// ─────────────────────────────────────────────────────────
double ADIPON_ss = ADIPON0 * (WT0 / BODY_WT) * (1 + 0.3 * E_SEM);
dxdt_ADIPONECTIN = KOUT_ADI * (ADIPON_ss - ADIPONECTIN);

// ─────────────────────────────────────────────────────────
// Insulin resistance (set-point form; IR_ss = IR0 at baseline)
// ─────────────────────────────────────────────────────────
double IR_from_FFA  = KFFA_IR  * (LIVER_FAT / LF0 - 1);
double IR_from_TNF  = KTNF_IR  * (TNFA / TNF0 - 1);
double IR_adipon_inh= 0.2 * (ADIPONECTIN / ADIPON0 - 1); // adiponectin protects
double IR_drug_inh  = 0.3 * E_SEM + 0.20 * E_EMP + 0.10 * E_RSM;
double IR_ss = IR0 * (1 + IR_from_FFA + IR_from_TNF - IR_adipon_inh) * (1 - IR_drug_inh);
IR_ss = (IR_ss < 0.5) ? 0.5 : IR_ss;
dxdt_INS_RES = KOUT_IR * (IR_ss - INS_RES);

// ─────────────────────────────────────────────────────────
// De novo lipogenesis (normalized; DNL_n = 1 at baseline)
// ─────────────────────────────────────────────────────────
double DNL_n = DNL0 * (1 + KDNL_IR * (INS_RES / IR0 - 1));
DNL_n = DNL_n * (1 - (0.30 * E_RSM + 0.25 * E_OCA));  // THRβ/FXR inhibit DNL (RSM 0.4->0.30: recalibrate PDFF toward MAESTRO -34/-39%)

// ─────────────────────────────────────────────────────────
// Liver fat (turnover; kin = KOUT_LF*LF0, influx = DNL + adipose-NEFA)
// ─────────────────────────────────────────────────────────
double WT_n      = BODY_WT / WT0;                       // adipose NEFA flux (dominant)
double LF_kin    = KOUT_LF * LF0 * (WDNL * DNL_n + WUPT * WT_n);
double LF_efflux = KOUT_LF
                 * (1 + 0.48 * E_RSM)   // THRβ → FAO ↑ (recal 0.6->0.48: model PDFF -40% -> ~-36%)
                 * (1 + 0.15 * E_SEM)   // GLP-1 → hepatic fat ↓
                 * (1 + WEMP_LF * E_EMP); // SGLT2i → hepatic fat ↓ (NEW; E-LIFT MRI-PDFF 16.2->11.3%)
dxdt_LIVER_FAT = LF_kin - LF_efflux * LIVER_FAT;

// ─────────────────────────────────────────────────────────
// Lipotoxicity (fold-change proxy from liver fat; = 1 at baseline)
// ─────────────────────────────────────────────────────────
double LIP_n = LIVER_FAT / LF0;

// ─────────────────────────────────────────────────────────
// Kupffer cell activation (saturable lipotoxicity drive)
// ─────────────────────────────────────────────────────────
double S_KUP  = 1 + GLIP_KUP * (LIP_n - 1) / (1 + GSAT_KUP * (LIP_n - 1));
double KUP_inh = 1 + 0.3 * E_OCA + 0.2 * E_SEM; // drug suppression (FXR/GLP-1)
dxdt_KUPFFER = KOUT_KUP * KUP0 * S_KUP / KUP_inh - KOUT_KUP * KUPFFER;

// ─────────────────────────────────────────────────────────
// TNF-α (Kupffer-driven, adiponectin-protected)
// ─────────────────────────────────────────────────────────
double S_TNF   = 1 + GKUP_TNF * (KUPFFER / KUP0 - 1);
double TNF_inh = 1 + 0.2 * E_OCA + 0.2 * E_SEM;
dxdt_TNFA = KOUT_TNF * TNF0 * S_TNF * (ADIPON0 / ADIPONECTIN) / TNF_inh
          - KOUT_TNF * TNFA;

// ─────────────────────────────────────────────────────────
// IL-6 (Kupffer-driven)
// ─────────────────────────────────────────────────────────
double S_IL6 = 1 + GKUP_IL6 * (KUPFFER / KUP0 - 1);
dxdt_IL6C = KOUT_IL6 * IL60 * S_IL6 * (1 - 0.2 * E_OCA) - KOUT_IL6 * IL6C;

// ─────────────────────────────────────────────────────────
// TGF-β1 (Kupffer + lipotoxicity drive; FXR/THRβ anti-fibrotic)
// ─────────────────────────────────────────────────────────
double S_TGF   = 1 + GKUP_TGF * (KUPFFER / KUP0 - 1) + GLIP_TGF * (LIP_n - 1);
double TGF_inh = 1 + 0.4 * E_OCA + 0.2 * E_RSM;
dxdt_TGFB = KOUT_TGF * TGF0 * S_TGF / TGF_inh - KOUT_TGF * TGFB;

// ─────────────────────────────────────────────────────────
// Hepatic Stellate Cell activation (TGF-β-driven)
// ─────────────────────────────────────────────────────────
double S_HSC   = 1 + GTGF_HSC * (TGFB / TGF0 - 1);
double HSC_inh = 1 + 0.4 * E_OCA + 0.15 * E_RSM + 0.1 * E_SEM;
dxdt_HSC = KOUT_HSC * HSC0 * S_HSC / HSC_inh - KOUT_HSC * HSC;

// ─────────────────────────────────────────────────────────
// Collagen / ECM (HSC-driven; slowest pool)
// ─────────────────────────────────────────────────────────
double S_COL   = 1 + GHSC_COL * (HSC / HSC0 - 1);
double COL_inh = 1 + 0.35 * E_OCA + 0.15 * E_RSM;
dxdt_COLLAGEN = KOUT_COL * COL0 * S_COL / COL_inh - KOUT_COL * COLLAGEN;

// ─────────────────────────────────────────────────────────
// ALT (hepatocellular injury = SUM of normalized TNF + lipotoxicity)
// ─────────────────────────────────────────────────────────
double S_ALT = 1 + GTNF_ALT * (TNFA / TNF0 - 1) + GLIP_ALT * (LIP_n - 1);
dxdt_ALT_CMT = KOUT_ALT * ALT0 * S_ALT - KOUT_ALT * ALT_CMT;

$TABLE
// Derived PK metrics
double Cp_RSM_out = RSM_CENT / V1_RSM;
double Cp_OCA_out = OCA_CENT / V1_OCA;
double Cp_SEM_out = SEM_CENT / V1_SEM;
double Cp_EMP_out = EMP_CENT / V1_EMP;

// Liver fat %
double LF_PCT = LIVER_FAT * 100;    // convert to percentage

// PDFF (MRI-PDFF proxy: correlated with histological fat %)
double PDFF = LF_PCT * 0.85;        // slightly lower than histological

// Fibrosis stage (0-4, continuous)
double FIB_SCORE = KFIBREG * COLLAGEN;
FIB_SCORE = (FIB_SCORE > 4) ? 4 : FIB_SCORE;

// NAS score components (simplified)
double NAS = 0;
NAS += (LF_PCT >= 5) ? ((LF_PCT < 33) ? 1 : (LF_PCT < 66) ? 2 : 3) : 0;  // steatosis
NAS += (ALT_CMT > 40) ? 1 : 0;   // lobular inflammation proxy
NAS += (TNFA > TNF0 * 1.5) ? 1 : 0; // ballooning proxy
NAS = (NAS > 8) ? 8 : NAS;

// HOMA-IR estimate
double HOMA_IR = INS_RES;

// Serum TG proxy (VLDL output ∝ liver fat + IR)
double TG_SERUM = 150 * (LIVER_FAT / LF0) * (INS_RES / IR0) * (1 - 0.5 * E_RSM);

// LDL-C proxy
double LDL_C = 120 * (1 + 0.1 * (LIVER_FAT / LF0 - 1)) * (1 - 0.25 * E_RSM);

// FIB-4 index (proxy; simplified)
double FIB4 = 1.8 * (COLLAGEN / COL0) * (ALT_CMT / ALT0);

// ELF score proxy
double ELF = 9 + 0.5 * FIB_SCORE;

// Drug effects visible
double EFFECT_RSM = E_RSM;
double EFFECT_OCA = E_OCA;
double EFFECT_SEM = E_SEM;
double EFFECT_EMP = E_EMP;

$CAPTURE
// NOTE: compartments (ALT_CMT, TNFA, IL6C, TGFB, COLLAGEN, HSC, KUPFFER,
// ADIPONECTIN, BODY_WT, INS_RES, LIVER_FAT) are returned automatically and
// must NOT be listed here — mrgsolve >=1.0 rejects compartments in $CAPTURE.
Cp_RSM_out Cp_OCA_out Cp_SEM_out Cp_EMP_out
LF_PCT PDFF FIB_SCORE NAS HOMA_IR
TG_SERUM LDL_C FIB4 ELF
EFFECT_RSM EFFECT_OCA EFFECT_SEM EFFECT_EMP
'

## ── Compile the model ──────────────────────────────────────
mod <- mcode("NAFLD_QSP", code)

## ── Helper: build event table ──────────────────────────────
build_regimen <- function(dose_rsm = 0, dose_oca = 0,
                          dose_sem = 0, dose_emp = 0,
                          duration_wk = 72) {
  hours <- duration_wk * 7 * 24
  ev_list <- list()

  if (dose_rsm > 0) {
    ev_list[["rsm"]] <- ev(
      amt  = dose_rsm, cmt = "RSM_GUT",
      ii   = 24, addl = ceiling(hours / 24) - 1
    )
  }
  if (dose_oca > 0) {
    ev_list[["oca"]] <- ev(
      amt  = dose_oca, cmt = "OCA_GUT",
      ii   = 24, addl = ceiling(hours / 24) - 1
    )
  }
  if (dose_sem > 0) {
    ev_list[["sem"]] <- ev(
      amt  = dose_sem, cmt = "SEM_SC",
      ii   = 168, addl = ceiling(hours / 168) - 1  # weekly
    )
  }
  if (dose_emp > 0) {
    ev_list[["emp"]] <- ev(
      amt  = dose_emp, cmt = "EMP_GUT",
      ii   = 24, addl = ceiling(hours / 24) - 1
    )
  }
  if (length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = "RSM_GUT"))
  Reduce("+", ev_list)
}

## ── Simulation time grid ───────────────────────────────────
tgrid <- tgrid(0, 72 * 7 * 24, delta = 12)  # 72 weeks, every 12h

## ── Scenario definitions ───────────────────────────────────
scenarios <- list(
  list(
    label     = "Placebo",
    params    = list(DOSE_RSM = 0, DOSE_OCA = 0, DOSE_SEM = 0, DOSE_EMP = 0),
    regimen   = build_regimen(),
    color     = "#999999",
    linetype  = "dashed"
  ),
  list(
    label     = "Resmetirom 100 mg QD",
    params    = list(DOSE_RSM = 1, DOSE_OCA = 0, DOSE_SEM = 0, DOSE_EMP = 0),
    regimen   = build_regimen(dose_rsm = 100),
    color     = "#E91E8C",
    linetype  = "solid"
  ),
  list(
    label     = "OCA 25 mg QD",
    params    = list(DOSE_RSM = 0, DOSE_OCA = 1, DOSE_SEM = 0, DOSE_EMP = 0),
    regimen   = build_regimen(dose_oca = 25),
    color     = "#1565C0",
    linetype  = "solid"
  ),
  list(
    label     = "Semaglutide 2.4 mg QW",
    params    = list(DOSE_RSM = 0, DOSE_OCA = 0, DOSE_SEM = 1, DOSE_EMP = 0),
    regimen   = build_regimen(dose_sem = 2.4),
    color     = "#388E3C",
    linetype  = "solid"
  ),
  list(
    label     = "Resmetirom + Semaglutide",
    params    = list(DOSE_RSM = 1, DOSE_OCA = 0, DOSE_SEM = 1, DOSE_EMP = 0),
    regimen   = build_regimen(dose_rsm = 100, dose_sem = 2.4),
    color     = "#F57F17",
    linetype  = "solid"
  ),
  list(
    label     = "Triple (RSM + OCA + Sema)",
    params    = list(DOSE_RSM = 1, DOSE_OCA = 1, DOSE_SEM = 1, DOSE_EMP = 0),
    regimen   = build_regimen(dose_rsm = 100, dose_oca = 25, dose_sem = 2.4),
    color     = "#B71C1C",
    linetype  = "solid"
  )
)

## ── Run simulations ────────────────────────────────────────
run_scenario <- function(sc) {
  param_update <- do.call(param, c(list(mod), sc$params))
  out <- mrgsim(param_update, ev = sc$regimen, tgrid = tgrid, obsonly = TRUE)
  df  <- as.data.frame(out)
  df$Scenario  <- sc$label
  df$Color     <- sc$color
  df$Linetype  <- sc$linetype
  df$Week      <- df$time / (7 * 24)
  df
}

cat("Running simulations...\n")
results <- bind_rows(lapply(scenarios, run_scenario))

## ── Plotting ───────────────────────────────────────────────
theme_qsp <- theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.background= element_rect(fill = "#ECEFF1"),
    panel.grid.minor= element_blank()
  )

scen_colors   <- setNames(sapply(scenarios, `[[`, "color"),
                           sapply(scenarios, `[[`, "label"))
scen_linetypes<- setNames(sapply(scenarios, `[[`, "linetype"),
                           sapply(scenarios, `[[`, "label"))

## Filter to weekly snapshots for cleaner plots
res_wk <- results %>% filter(abs(Week - round(Week)) < 0.1)

p_lf <- ggplot(res_wk, aes(Week, LF_PCT, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey50") +
  labs(title = "Hepatic Fat Fraction", y = "Liver Fat (%)", x = "Week") +
  theme_qsp

p_fib <- ggplot(res_wk, aes(Week, FIB_SCORE, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(title = "Fibrosis Score (0–4)", y = "Fibrosis Stage", x = "Week") +
  theme_qsp

p_alt <- ggplot(res_wk, aes(Week, ALT_CMT, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  geom_hline(yintercept = 35, linetype = "dotted", color = "grey50") +
  labs(title = "Serum ALT", y = "ALT (U/L)", x = "Week") +
  theme_qsp

p_nas <- ggplot(res_wk, aes(Week, NAS, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  scale_y_continuous(breaks = 0:8, limits = c(0, 8)) +
  labs(title = "NAFLD Activity Score (NAS)", y = "NAS (0–8)", x = "Week") +
  theme_qsp

p_homa <- ggplot(res_wk, aes(Week, HOMA_IR, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  labs(title = "HOMA-IR (Insulin Resistance)", y = "HOMA-IR", x = "Week") +
  theme_qsp

p_tg <- ggplot(res_wk, aes(Week, TG_SERUM, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  geom_hline(yintercept = 150, linetype = "dotted", color = "grey50") +
  labs(title = "Serum Triglycerides", y = "TG (mg/dL)", x = "Week") +
  theme_qsp

p_tnf <- ggplot(res_wk, aes(Week, TNFA, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  labs(title = "TNF-α (Inflammation)", y = "TNF-α (rel. units)", x = "Week") +
  theme_qsp

p_wt <- ggplot(res_wk, aes(Week, BODY_WT, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_linetype_manual(values = scen_linetypes) +
  labs(title = "Body Weight", y = "Weight (kg)", x = "Week") +
  theme_qsp

## Combined dashboard plot
combined <- (p_lf | p_fib) / (p_alt | p_nas) / (p_homa | p_tg) / (p_tnf | p_wt) +
  plot_annotation(
    title    = "NAFLD/NASH QSP Model — Treatment Scenario Comparison",
    subtitle = "72-week simulation · Resmetirom · OCA · Semaglutide · Combination",
    caption  = "Effect sizes qualitatively consistent with MAESTRO-NASH, REGENERATE, ESSENCE; disease parameters phenomenological (not formally fitted)",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 12))
  )

cat("Saving combined plot...\n")
ggsave("nafld_simulation_results.png", combined, width = 14, height = 18, dpi = 150)
cat("Plot saved: nafld_simulation_results.png\n")

## ── PK Profile (Resmetirom steady state) ──────────────────
pk_ev <- ev(amt = 100, cmt = "RSM_GUT", ii = 24, addl = 13)
pk_grid <- tgrid(300, 336, delta = 0.5)  # last 36h of 14 days SS
pk_sim <- mrgsim(param(mod, DOSE_RSM = 1), ev = pk_ev,
                 tgrid = pk_grid, obsonly = TRUE)
pk_df <- as.data.frame(pk_sim)
pk_df$Week <- pk_df$time / 24

p_pk <- ggplot(pk_df %>% filter(time > 300), aes(time - 312, Cp_RSM_out)) +
  geom_line(color = "#E91E8C", linewidth = 1.2) +
  labs(
    title    = "Resmetirom Steady-State PK (Day 14)",
    subtitle = "100 mg QD oral; 2-compartment model",
    x        = "Hours post-dose",
    y        = "Plasma Concentration (µg/L)"
  ) +
  theme_qsp

ggsave("nafld_resmetirom_pk.png", p_pk, width = 7, height = 4, dpi = 150)
cat("PK plot saved: nafld_resmetirom_pk.png\n")

## ── Endpoint summary table at Week 72 ──────────────────────
summary_tab <- res_wk %>%
  filter(abs(Week - 72) < 0.1) %>%
  select(Scenario, LF_PCT, PDFF, FIB_SCORE, NAS, ALT_CMT,
         HOMA_IR, TG_SERUM, LDL_C, BODY_WT) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

cat("\n── Week-72 Endpoint Summary ──────────────────────────────\n")
print(summary_tab)

## ── Dose-response for Resmetirom ──────────────────────────
doses_rsm <- c(0, 25, 50, 80, 100, 150)
dose_response <- lapply(doses_rsm, function(d) {
  ev_d  <- build_regimen(dose_rsm = d, duration_wk = 52)
  param_d <- param(mod, DOSE_RSM = if (d > 0) 1 else 0)
  out_d <- mrgsim(param_d, ev = ev_d,
                  tgrid = tgrid(0, 52 * 168, delta = 168),
                  obsonly = TRUE)
  df_d <- as.data.frame(out_d) %>%
    filter(abs(time - 52 * 168) < 5) %>%
    mutate(Dose = d)
  df_d
})
dr_df <- bind_rows(dose_response)

p_dr <- ggplot(dr_df, aes(Dose, LF_PCT)) +
  geom_line(color = "#E91E8C", linewidth = 1.2) +
  geom_point(size = 3, color = "#880E4F") +
  labs(
    title    = "Resmetirom Dose–Response (Week 52)",
    subtitle = "Hepatic fat fraction as primary endpoint",
    x        = "Daily Dose (mg)",
    y        = "Hepatic Fat (%)"
  ) +
  theme_qsp

ggsave("nafld_dose_response.png", p_dr, width = 7, height = 4, dpi = 150)
cat("Dose-response plot saved: nafld_dose_response.png\n")
cat("NAFLD/NASH mrgsolve simulation complete.\n")
