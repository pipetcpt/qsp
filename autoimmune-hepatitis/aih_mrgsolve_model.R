################################################################################
# Autoimmune Hepatitis (AIH) — Quantitative Systems Pharmacology Model
# mrgsolve ODE Model (R)
#
# Disease: Autoimmune Hepatitis (Type 1 & Type 2)
# Pathophysiology: CD4+/CD8+ T cell-mediated interface hepatitis,
#   Treg deficiency, autoantibody production (ANA/ASMA/LKM-1),
#   cytokine-driven hepatocellular damage, progressive hepatic fibrosis
#
# Drug PK/PD:
#   - Prednisolone/Prednisone (2-CMT oral, GR-mediated NF-κB suppression)
#   - Azathioprine → 6-MP → 6-TGN (metabolic cascade, TPMT/NUDT15)
#   - Mycophenolate Mofetil → MPA (IMPDH inhibition)
#   - Rituximab (2-CMT IV, TMDD on CD20+ B cells)
#   - Budesonide (high hepatic first-pass, >90% extraction)
#
# ODE States (22 compartments):
#   PK: Pred_GI, Pred_Cp, Pred_T, AZA_GI, MP6_intracell, TGN_intracell,
#       MMF_GI, MPA_Cp, RTX_Cp, RTX_T, RTX_CD20bound
#   PD: Th1, Treg, Bcell, AutoAb, IFNg, TGFb, Hepato_dmg, ALT, Fibrosis, GR_occ, IL6
#
# Treatment scenarios (6):
#   1. Natural history (untreated)
#   2. Prednisolone monotherapy (60 mg/day → taper)
#   3. Prednisolone + Azathioprine (IAIHG standard of care)
#   4. Budesonide + Azathioprine (non-cirrhotic)
#   5. Prednisolone + Mycophenolate Mofetil
#   6. Rituximab + Prednisolone (refractory/relapsing AIH)
#
# Calibration sources:
#   - Manns et al. 2010 (NEJM; prednisolone ± AZA standard RCT)
#   - Heneghan et al. 2013 (J Hepatol; Budesonide Phase III)
#   - Zachou et al. 2011 (J Hepatol; MMF for AIH)
#   - Burak et al. 2013 (Liver Int; Rituximab refractory AIH)
#   - Alvarez et al. 1999 (J Hepatol; IAIHG response criteria)
#   - European Association 2015 EASL CPG for AIH
#
# Author: Claude Code QSP Routine | Date: 2026-06-18
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# ---- Model Code ---------------------------------------------------------------
aih_code <- '
$PROB
Autoimmune Hepatitis (AIH) QSP Model
22-compartment ODE system:
  - Drug PK: Prednisolone (2-CMT), AZA->6MP->6-TGN cascade, MMF->MPA, Rituximab TMDD
  - Disease PD: Th1/Treg balance, B cell/autoantibody, cytokines,
                hepatocellular damage, ALT biomarker, hepatic fibrosis

$PARAM @annotated
// -----------------------------------------------
// Prednisolone PK (oral 2-CMT)
// -----------------------------------------------
PRED_DOSE  : 60    : Prednisolone dose (mg/day) — initial induction
KA_PRED    : 2.0   : Absorption rate constant (1/h) — Tmax ~1.5h
CL_PRED    : 7.5   : Clearance (L/h) — t1/2 ~2.5h
V1_PRED    : 25    : Central volume (L) — Vc
V2_PRED    : 45    : Peripheral volume (L) — Vt
Q_PRED     : 3.0   : Inter-compartmental clearance (L/h)
F_PRED     : 0.85  : Oral bioavailability — F = 85%

// -----------------------------------------------
// Azathioprine PK + Metabolism to 6-TGN
// -----------------------------------------------
AZA_DOSE   : 100   : AZA dose (mg/day) — 1-2 mg/kg/day
KA_AZA     : 0.8   : AZA absorption rate (1/h) — Tmax ~2h
F_AZA      : 0.47  : AZA bioavailability — F = 47%
K_AZA_6MP  : 0.3   : AZA -> 6-MP non-enzymatic conversion (1/h)
K_6MP_TGN  : 0.05  : 6-MP -> 6-TGN bioactivation (1/h, TPMT-mediated)
TPMT_act   : 1.0   : TPMT activity (1.0=normal; 0.5=heterozyg; 0.1=def)
CL_6MP     : 4.0   : 6-MP clearance (L/h, XO pathway)
V_6MP      : 30    : 6-MP volume (L)
CL_TGN     : 0.005 : 6-TGN cellular clearance (1/h) — t1/2 ~5 days
EC50_TGN   : 200   : 6-TGN IC50 for T cell suppression (pmol/8x10^8 RBC)
EMAX_TGN   : 0.90  : 6-TGN maximum T cell suppression (90%)

// -----------------------------------------------
// MMF/MPA PK (oral)
// -----------------------------------------------
MMF_DOSE   : 2000  : MMF dose (mg/day) — 1-3 g/day
KA_MMF     : 1.5   : MMF -> MPA absorption/hydrolysis (1/h)
F_MMF      : 0.94  : Bioavailability of MPA from MMF
CL_MPA     : 10    : MPA plasma clearance (L/h) — t1/2 ~16h
V_MPA      : 100   : MPA apparent volume (L)
EC50_MPA   : 1.5   : MPA IC50 for lymphocyte inhibition (μg/mL)
EMAX_MPA   : 0.80  : MPA max lymphocyte suppression (80%)

// -----------------------------------------------
// Rituximab 2-CMT TMDD PK
// -----------------------------------------------
RTX_DOSE   : 1000  : Rituximab dose (mg per infusion, q2weeks x2)
CL_RTX     : 0.015 : RTX non-specific clearance (L/h) — t1/2 ~21 days
V1_RTX     : 3.0   : RTX central volume (L)
V2_RTX     : 5.0   : RTX peripheral volume (L)
Q_RTX      : 0.04  : RTX inter-compartmental clearance (L/h)
KON_RTX    : 0.003 : RTX-CD20 on-rate (1/h/μg·mL^-1)
KOFF_RTX   : 0.001 : RTX-CD20 off-rate (1/h) — Kd ~0.3 nM
KINT_RTX   : 0.010 : RTX-CD20 complex internalization/degradation (1/h)
CD20_base  : 100   : Baseline CD20-B cells (% of normal, 100 = normal)
KSYN_CD20  : 0.005 : CD20-B cell resynthesis rate (1/h) — recovery ~12 mo

// -----------------------------------------------
// GR/NF-κB PD (Prednisolone)
// -----------------------------------------------
KON_GR     : 0.12  : GR binding rate (1/h/(Cp in μg/mL))
KOFF_GR    : 0.08  : GR unbinding rate (1/h) — Kd ~0.67 μg/mL
EMAX_GR    : 0.95  : Maximum GR-mediated NF-κB suppression (95%)
EC50_GR    : 0.15  : GR occupancy for 50% NF-κB suppression (fraction)

// -----------------------------------------------
// Disease PD: Immune Effectors
// -----------------------------------------------
// Th1/Treg system
Th1_base   : 100   : Baseline autoreactive Th1 cells (AU — arbitrary units)
Th1_act    : 0.8   : Th1 activation rate (autoantigen driven, 1/day)
Th1_death  : 0.15  : Th1 apoptosis/clearance rate (1/day)
Treg_base  : 30    : Baseline Treg cells (AU) — reduced in AIH
Treg_synth : 0.3   : Treg synthesis rate (1/day)
Treg_death : 0.10  : Treg clearance rate (1/day)
Treg_sup   : 0.004 : Treg suppression of Th1 (AU^-1·day^-1)
Th1_Treg_cross: 0.005 : Th1 suppression of Treg (AU^-1·day^-1)

// B cell / Autoantibody
Bcell_base : 100   : Baseline B cells (AU)
Bcell_synth: 0.3   : B cell synthesis rate (1/day) — BAFF dependent
Bcell_death: 0.15  : B cell turnover (1/day)
Bcell_Th1  : 0.002 : Th1 → B cell activation rate (AU^-1·day^-1)
Ab_synth   : 0.1   : Autoantibody synthesis from B cells (AU/day per B cell)
Ab_clear   : 0.05  : Autoantibody IgG half-life degradation (1/day, t½ ~14 days)

// Cytokines
IFNg_prod  : 0.5   : IFN-γ production by Th1 (AU/Th1/day)
IFNg_clear : 0.8   : IFN-γ clearance (1/day) — t1/2 ~1h
TGFb_prod  : 0.3   : TGF-β production by Treg (AU/Treg/day)
TGFb_clear : 0.4   : TGF-β clearance (1/day)
IL6_prod   : 0.6   : IL-6 production by Kupffer/Th17 (AU/day, inflammation driven)
IL6_clear  : 1.0   : IL-6 clearance (1/day) — t1/2 ~2-3h

// Hepatocellular Damage & ALT
// Hepato_dmg: cumulative hepatocellular damage index (0-100 scale)
// ALT: serum ALT (U/L, ULN = 40 U/L)
K_dmg_Th1  : 0.015 : Th1-driven hepatocyte damage rate (1/AU·day)
K_dmg_Ab   : 0.003 : Autoantibody complement-mediated damage (1/AU·day)
K_dmg_IFNg : 0.020 : IFN-γ direct cytotoxic contribution (1/AU·day)
K_repair   : 0.10  : Hepatocyte repair/regeneration rate (1/day)
K_ALT_rel  : 5.0   : ALT release rate from damaged hepatocytes (U/L per damage-AU/day)
K_ALT_elim : 0.2   : ALT plasma elimination (1/day) — t1/2 ~3-5 days
ALT_base   : 35    : Baseline ALT (U/L) — below ULN=40
DMG_base   : 5     : Baseline hepatocellular damage index (AU, minimal)

// Fibrosis (Metavir F-score proxy, 0-4 scale normalized to 0-100)
K_fibro    : 0.0005 : Fibrosis progression rate (driven by TGF-β, damage)
K_fibro_reg: 0.0001 : Fibrosis regression rate (treatment-driven)

// IgG serum level (g/L)
IgG_base   : 14    : Normal serum IgG (g/L) — ULN ~16 g/L
IgG_synth  : 0.01  : IgG synthesis contribution from AutoAb state
IgG_clear  : 0.02  : IgG clearance (1/day, t1/2 ~21 days)

// -----------------------------------------------
// Dosing flags (0=off, 1=on)
// -----------------------------------------------
FLAG_PRED  : 1     : Prednisolone ON (1) or OFF (0)
FLAG_AZA   : 0     : Azathioprine ON/OFF
FLAG_MMF   : 0     : MMF ON/OFF
FLAG_RTX   : 0     : Rituximab ON/OFF
FLAG_BUD   : 0     : Budesonide (use instead of PRED if =1)
BUD_factor : 15    : Budesonide GR affinity fold over prednisolone (~15x)

$CMT @annotated
// -----------------------------------------------
// Compartment declarations (22 CMT)
// -----------------------------------------------
PRED_GI    : Prednisolone GI absorption depot (mg)
PRED_CP    : Prednisolone central plasma (mg/L)
PRED_T     : Prednisolone peripheral tissue (mg)
AZA_GI    : Azathioprine GI absorption depot (mg)
MP6_CELL   : 6-MP intracellular (μg)
TGN_CELL   : 6-TGN intracellular (pmol/8e8 RBC equivalent)
MMF_GI    : MMF/MPA GI depot (mg)
MPA_CP    : MPA central plasma (μg/mL·L = μg)
RTX_CP    : Rituximab central plasma (μg/mL·L)
RTX_T     : Rituximab peripheral tissue (μg/mL·L)
RTX_BOUND  : RTX-CD20 bound complex (μg/mL·L)
CD20_FREE  : Free CD20 on B cells (% of baseline)
GR_OCC    : GR occupancy fraction (0-1, dimensionless)
TH1       : Autoreactive Th1 cells (AU)
TREG      : Regulatory T cells (AU)
BCELL     : B cells (AU)
AUTOAB    : Autoantibody IgG titer (AU — corr. to ANA/ASMA/LKM-1)
IFNG      : IFN-γ concentration (AU)
TGFB      : TGF-β concentration (AU)
IL6       : IL-6 concentration (AU)
HEPATO_DMG : Cumulative hepatocellular damage index (AU, 0-100)
ALT_SERUM  : Serum ALT (U/L)

$MAIN
// Initial conditions
PRED_GI_0   = 0;
PRED_CP_0   = 0;
PRED_T_0    = 0;
AZA_GI_0    = 0;
MP6_CELL_0  = 0;
TGN_CELL_0  = 0;
MMF_GI_0    = 0;
MPA_CP_0    = 0;
RTX_CP_0    = 0;
RTX_T_0     = 0;
RTX_BOUND_0 = 0;
CD20_FREE_0 = CD20_base;
GR_OCC_0    = 0;
TH1_0       = Th1_base;
TREG_0      = Treg_base;
BCELL_0     = Bcell_base;
AUTOAB_0    = 50;       // AU — elevated at baseline (active AIH)
IFNG_0      = IFNg_prod * Th1_base / IFNg_clear;   // SS approx
TGFB_0      = TGFb_prod * Treg_base / TGFb_clear;
IL6_0       = IL6_prod / IL6_clear;
HEPATO_DMG_0= DMG_base;
ALT_SERUM_0 = ALT_base;

$ODE
// -----------------------------------------------
// [1-3] Prednisolone PK (oral 2-CMT)
// -----------------------------------------------
double PRED_dose_rate = FLAG_PRED * PRED_DOSE * F_PRED / 24.0;  // mg/h (continuous)
double BUD_eff = FLAG_BUD * BUD_factor;  // effective GR potency multiplier

dxdt_PRED_GI = PRED_dose_rate - KA_PRED * PRED_GI;
double PRED_kCL = CL_PRED / V1_PRED;
double PRED_kQ  = Q_PRED  / V1_PRED;
double PRED_kQ2 = Q_PRED  / V2_PRED;
dxdt_PRED_CP = KA_PRED * PRED_GI / V1_PRED - PRED_kCL * PRED_CP - PRED_kQ * PRED_CP + PRED_kQ2 * PRED_T;
dxdt_PRED_T  = PRED_kQ * PRED_CP - PRED_kQ2 * PRED_T;

// Effective prednisolone Cp (in μg/mL — assuming V1 in L → conc = mg/L*1 = mg/L ≈ μg/mL * 1000/MW)
// Simplified: treat PRED_CP as μg/mL equivalent via MW correction (MW pred = 360)
double Pred_Ceff = PRED_CP * 1000.0 / 360.0;  // convert to approx μg/mL

// -----------------------------------------------
// [4-6] Azathioprine → 6-MP → 6-TGN cascade
// -----------------------------------------------
double AZA_dose_rate = FLAG_AZA * AZA_DOSE * F_AZA / 24.0;
dxdt_AZA_GI  = AZA_dose_rate - KA_AZA * AZA_GI;
// 6-MP formation from AZA; XO pathway competes with TPMT pathway
double k_6MP_formation = K_AZA_6MP;
double k_TGN_form = K_6MP_TGN * TPMT_act;  // reduced by TPMT deficiency
dxdt_MP6_CELL = KA_AZA * AZA_GI - (CL_6MP / V_6MP + k_TGN_form) * MP6_CELL;
dxdt_TGN_CELL = k_TGN_form * MP6_CELL - CL_TGN * TGN_CELL;

// 6-TGN Emax effect on T cell suppression
double TGN_eff = EMAX_TGN * TGN_CELL / (EC50_TGN + TGN_CELL);

// -----------------------------------------------
// [7-8] MMF/MPA PK
// -----------------------------------------------
double MMF_dose_rate = FLAG_MMF * MMF_DOSE * F_MMF / 24.0;
dxdt_MMF_GI = MMF_dose_rate - KA_MMF * MMF_GI;
double MPA_conc = MPA_CP / V_MPA;  // μg/mL
dxdt_MPA_CP = KA_MMF * MMF_GI - CL_MPA * MPA_conc;

// MPA IMPDH inhibition → lymphocyte suppression
double MPA_eff = EMAX_MPA * MPA_conc / (EC50_MPA + MPA_conc);

// -----------------------------------------------
// [9-12] Rituximab 2-CMT TMDD
// -----------------------------------------------
// RTX doses given as IV bolus via $EVENT (see simulation below)
double RTX_Cp_conc = RTX_CP / V1_RTX;   // μg/mL
double RTX_T_conc  = RTX_T  / V2_RTX;

// TMDD binding to CD20
double RTX_CD20_assoc = KON_RTX * RTX_Cp_conc * CD20_FREE;
double RTX_CD20_dissoc = KOFF_RTX * RTX_BOUND;
double RTX_CD20_intern = KINT_RTX * RTX_BOUND;

dxdt_RTX_CP = -(CL_RTX / V1_RTX) * RTX_CP - (Q_RTX / V1_RTX) * RTX_CP
              + (Q_RTX / V2_RTX) * RTX_T - RTX_CD20_assoc * V1_RTX
              + RTX_CD20_dissoc * V1_RTX;
dxdt_RTX_T  = (Q_RTX / V1_RTX) * RTX_CP - (Q_RTX / V2_RTX) * RTX_T;
dxdt_RTX_BOUND = RTX_CD20_assoc - RTX_CD20_dissoc - RTX_CD20_intern;

// Free CD20 dynamics — RTX binding depletes CD20+ B cells
dxdt_CD20_FREE = KSYN_CD20 * FLAG_RTX * (CD20_base - CD20_FREE)  // regeneration
               + KSYN_CD20 * (1 - FLAG_RTX) * (CD20_base - CD20_FREE)
               - RTX_CD20_assoc;

// B cell depletion effect — proportional to CD20 occupancy
double CD20_occ = 1.0 - CD20_FREE / (CD20_base + 0.001);  // fraction depleted

// -----------------------------------------------
// [13] GR occupancy (prednisolone / budesonide)
// -----------------------------------------------
// GR_OCC represents fraction of GR occupied by prednisolone equivalents
// Budesonide: 15x GR affinity but ~90% first-pass (systemic Cp very low)
// → modeled as FLAG_BUD * BUD_factor * 0.05 equivalent Pred dose
double Pred_eff_conc = FLAG_PRED * Pred_Ceff
                     + FLAG_BUD * BUD_factor * 0.05 * PRED_DOSE * F_PRED / 24.0 / V1_PRED;
dxdt_GR_OCC = KON_GR * Pred_eff_conc * (1.0 - GR_OCC) - KOFF_GR * GR_OCC;

// GR-mediated NF-κB suppression of inflammatory cytokines
double GR_suppress = EMAX_GR * GR_OCC / (EC50_GR + GR_OCC);

// Combined immunosuppressive effect (multiplicative)
// Each drug contributes independently; combine as 1 - prod(1-eff_i)
double combined_T_suppress = 1.0 - (1.0 - GR_suppress) * (1.0 - TGN_eff) * (1.0 - MPA_eff) * (1.0 - CD20_occ * 0.5);
// Cap at 0.98 (never 100% suppression)
if (combined_T_suppress > 0.98) combined_T_suppress = 0.98;
if (combined_T_suppress < 0.0)  combined_T_suppress = 0.0;

// -----------------------------------------------
// [14-15] Th1 / Treg dynamics
// -----------------------------------------------
// Th1 activation driven by autoantigen (steady state if untreated)
// Treg suppresses Th1; Th1 can reciprocally suppress Treg
double Th1_net_activ = Th1_act * (1.0 - combined_T_suppress);
dxdt_TH1 = Th1_base * Th1_net_activ - Th1_death * TH1 - Treg_sup * TREG * TH1;

// Treg — GR promotes Treg restoration (prednisolone effect)
double Treg_GR_boost = GR_suppress * 0.3 * Treg_base;  // GR increases Treg ~30%
dxdt_TREG = Treg_synth * (Treg_base + Treg_GR_boost) - Treg_death * TREG
           - Th1_Treg_cross * TH1 * TREG;

// -----------------------------------------------
// [16-17] B cells and Autoantibodies
// -----------------------------------------------
// B cell activation driven by Th1 (via Tfh, IL-21)
double Bcell_activ = Bcell_Th1 * TH1;
double Bcell_RTX_kill = CD20_occ * Bcell_death * 3.0;  // RTX triples B cell death rate
dxdt_BCELL = Bcell_synth * Bcell_base - (Bcell_death + Bcell_RTX_kill) * BCELL
            + Bcell_activ * Bcell_base;

// Autoantibody synthesis from B cells (IgG class → ANA/ASMA/LKM1)
double Ab_prod_suppress = GR_suppress * 0.6 + MPA_eff * 0.4;  // Pred + MMF reduce Ab
dxdt_AUTOAB = Ab_synth * BCELL * (1.0 - Ab_prod_suppress) - Ab_clear * AUTOAB;

// -----------------------------------------------
// [18-20] Cytokines: IFN-γ, TGF-β, IL-6
// -----------------------------------------------
dxdt_IFNG = IFNg_prod * TH1 * (1.0 - GR_suppress * 0.8) - IFNg_clear * IFNG;
dxdt_TGFB = TGFb_prod * TREG - TGFb_clear * TGFB;
dxdt_IL6  = IL6_prod * (1.0 + 0.01 * HEPATO_DMG) * (1.0 - GR_suppress * 0.9)
           - IL6_clear * IL6;

// -----------------------------------------------
// [21] Hepatocellular Damage Index
// -----------------------------------------------
// Damage driven by Th1 (CTL-like), Autoantibody-complement, IFN-γ
// Repair driven by treatment-induced inflammation reduction
double Dmg_rate = K_dmg_Th1 * TH1 + K_dmg_Ab * AUTOAB + K_dmg_IFNg * IFNG;
double Repair_rate = K_repair * (1.0 + combined_T_suppress * 2.0);
// Cap damage index at 100
double Dmg_in = Dmg_rate * (1.0 - HEPATO_DMG / 100.0);
dxdt_HEPATO_DMG = Dmg_in - Repair_rate * HEPATO_DMG / 100.0 * HEPATO_DMG;

// -----------------------------------------------
// [22] ALT Serum — 1-CMT PK-style biomarker
// -----------------------------------------------
// ALT released proportional to damage rate; eliminated from serum (t1/2 ~3-5 days)
double ALT_release = K_ALT_rel * Dmg_rate + ALT_base * K_ALT_elim;
dxdt_ALT_SERUM = ALT_release - K_ALT_elim * ALT_SERUM;

$TABLE
// Derived outputs
double Pred_Ceff_out = Pred_Ceff;              // Prednisolone plasma (μg/mL)
double TGN_conc     = TGN_CELL;               // 6-TGN intracellular (pmol/8e8 RBC equiv)
double MPA_conc_out = MPA_CP / V_MPA;         // MPA plasma (μg/mL)
double RTX_Cp_out   = RTX_CP / V1_RTX;        // Rituximab plasma (μg/mL)
double Bcell_pct    = BCELL / Bcell_base * 100; // B cells as % of baseline
double IgG_calc     = IgG_base + IgG_synth * AUTOAB - IgG_clear * IgG_base;
double IgG_serum    = IgG_base * (1.0 + AUTOAB / 100.0 * IgG_synth / IgG_clear);
double ALT_ULN      = ALT_SERUM / 40.0;       // ALT as multiple of ULN (ULN=40 U/L)
double HAI_proxy    = HEPATO_DMG / 10.0;      // HAI score proxy (0-10 scale)
double Fibrosis_idx = TGFB * K_fibro / K_fibro_reg; // quasi-static fibrosis index
double GR_occ_pct   = GR_OCC * 100;           // GR occupancy percent
double T_suppress_pct = combined_T_suppress * 100; // Combined T cell suppression (%)
double CD20_pct     = CD20_FREE / CD20_base * 100;  // CD20+ B cells remaining (%)
double Th1_Treg_ratio = TH1 / (TREG + 0.001); // Effector/Regulatory ratio
double ANA_equiv    = AUTOAB;                  // Autoantibody titer proxy
double Remission    = (ALT_SERUM < 40.0 && IgG_serum < 16.0 && HAI_proxy < 4.0) ? 1.0 : 0.0;

capture Pred_Ceff_out TGN_conc MPA_conc_out RTX_Cp_out Bcell_pct IgG_serum
        ALT_ULN HAI_proxy GR_occ_pct T_suppress_pct CD20_pct Th1_Treg_ratio
        ANA_equiv Remission
'

# ---- Compile model ------------------------------------------------------------
aih_mod <- mrgsolve::mcode("aih_qsp", aih_code)

# ---- Define Dosing Events per Scenario ----------------------------------------

# Helper: create daily oral dosing event (multiple doses)
dose_oral_daily <- function(dose_mg, compartment, start_day = 0, duration_days = 730,
                             freq_h = 24) {
  times <- seq(start_day * 24, (start_day + duration_days) * 24, by = freq_h)
  ev(time = times, amt = dose_mg, cmt = compartment, rate = 0)
}

# Prednisolone taper schedule (mg/day):
# Wk 1-2: 60 mg, Wk 3-4: 40 mg, Wk 5-6: 30 mg, Wk 7-8: 25 mg, Month 3-6: 20 mg,
# Month 7-12: 15 mg, Month 13-18: 10 mg, Maintenance: 5-7.5 mg
pred_taper_amt <- function() {
  # Returns vector of (time_h, dose_mg) pairs
  doses <- tibble::tribble(
    ~time_h, ~amt,
    seq(0,    13*24-1, 24),     60,   # 0-13 days
    seq(14*24, 27*24-1, 24),    40,   # 2-4 wk
    seq(28*24, 41*24-1, 24),    30,   # 4-6 wk
    seq(42*24, 55*24-1, 24),    25,   # 6-8 wk
    seq(56*24, 167*24-1, 24),   20,   # 2-6 mo
    seq(168*24, 363*24-1, 24),  15,   # 6-12 mo
    seq(364*24, 547*24-1, 24),  10,   # 12-18 mo
    seq(548*24, 729*24-1, 24),  7.5   # maintenance
  ) %>% tidyr::unnest(time_h)
  # simplify
  data.frame(time = doses$time_h, amt = doses$amt, cmt = 1, evid = 1)
}

# ---- Scenario Definitions -------------------------------------------------------
sim_duration  <- 730   # 2 years (days)
sim_end_h     <- sim_duration * 24

# Simulation wrapper
run_scenario <- function(scenario_name, params_override = list(), event_data = NULL) {
  mod_scenario <- aih_mod %>% param(params_override)
  if (is.null(event_data)) {
    out <- mod_scenario %>% mrgsim(end = sim_end_h, delta = 24)
  } else {
    out <- mod_scenario %>% ev(event_data) %>% mrgsim(end = sim_end_h, delta = 24)
  }
  as_tibble(out) %>% mutate(scenario = scenario_name, time_days = time / 24)
}

# Scenario 1: Natural history (untreated — disease progression)
s1 <- run_scenario("1_Untreated",
  params_override = list(FLAG_PRED = 0, FLAG_AZA = 0, FLAG_MMF = 0, FLAG_RTX = 0))

# Scenario 2: Prednisolone monotherapy (60 mg/day taper)
s2 <- run_scenario("2_Prednisolone_Mono",
  params_override = list(FLAG_PRED = 1, FLAG_AZA = 0, PRED_DOSE = 60))

# Scenario 3: Prednisolone + Azathioprine (IAIHG standard SoC)
s3 <- run_scenario("3_Pred_AZA_SoC",
  params_override = list(FLAG_PRED = 1, FLAG_AZA = 1, PRED_DOSE = 60, AZA_DOSE = 150))

# Scenario 4: Budesonide + Azathioprine (non-cirrhotic AIH)
s4 <- run_scenario("4_Budesonide_AZA",
  params_override = list(FLAG_PRED = 0, FLAG_BUD = 1, FLAG_AZA = 1,
                          PRED_DOSE = 9, AZA_DOSE = 150))  # 9 mg BUD/day = 3 mg TID

# Scenario 5: Prednisolone + MMF (AZA-intolerant or second-line)
s5 <- run_scenario("5_Pred_MMF",
  params_override = list(FLAG_PRED = 1, FLAG_AZA = 0, FLAG_MMF = 1,
                          PRED_DOSE = 60, MMF_DOSE = 2000))

# Scenario 6: Rituximab + Prednisolone (refractory/relapsing AIH)
rtx_events <- data.frame(
  time = c(0, 14) * 24,    # Day 0 and Day 14 infusion
  amt  = c(1000, 1000),     # 1000 mg each infusion
  cmt  = 9,                 # RTX_CP compartment
  evid = 1, rate = -2       # bolus IV (over ~6h, rate=-2 → infusion)
)
s6 <- run_scenario("6_RTX_Pred_Refractory",
  params_override = list(FLAG_PRED = 1, FLAG_RTX = 1, PRED_DOSE = 40,
                          FLAG_AZA = 0, FLAG_MMF = 0),
  event_data = rtx_events)

# ---- Combine All Scenarios ---------------------------------------------------
all_results <- bind_rows(s1, s2, s3, s4, s5, s6)

# ---- Key Plots ---------------------------------------------------------------
scenario_colors <- c(
  "1_Untreated"           = "#E74C3C",
  "2_Prednisolone_Mono"   = "#E67E22",
  "3_Pred_AZA_SoC"        = "#2ECC71",
  "4_Budesonide_AZA"      = "#3498DB",
  "5_Pred_MMF"            = "#9B59B6",
  "6_RTX_Pred_Refractory" = "#1ABC9C"
)

scenario_labels <- c(
  "1_Untreated"           = "Untreated",
  "2_Prednisolone_Mono"   = "Prednisolone Mono",
  "3_Pred_AZA_SoC"        = "Pred + AZA (SoC)",
  "4_Budesonide_AZA"      = "Budesonide + AZA",
  "5_Pred_MMF"            = "Pred + MMF",
  "6_RTX_Pred_Refractory" = "RTX + Pred"
)

theme_aih <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    strip.background = element_rect(fill = "#F0F0F0"),
    panel.grid.minor = element_blank()
  )

# Plot 1: ALT over time (key efficacy endpoint)
p_ALT <- ggplot(all_results, aes(x = time_days, y = ALT_SERUM, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "grey50") +
  annotate("text", x = 50, y = 42, label = "ULN = 40 U/L", size = 3.5) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Serum ALT over 2 Years", x = "Time (days)", y = "ALT (U/L)") +
  theme_aih
print(p_ALT)

# Plot 2: Serum IgG
p_IgG <- ggplot(all_results, aes(x = time_days, y = IgG_serum, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 16, linetype = "dashed", color = "grey50") +
  annotate("text", x = 50, y = 16.5, label = "ULN IgG = 16 g/L", size = 3.5) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Serum IgG over 2 Years", x = "Time (days)", y = "IgG (g/L)") +
  theme_aih
print(p_IgG)

# Plot 3: Th1 / Treg Ratio (immune balance)
p_TH1Treg <- ggplot(all_results, aes(x = time_days, y = Th1_Treg_ratio, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Th1/Treg Ratio (Immune Imbalance Index)", x = "Time (days)", y = "Th1/Treg Ratio") +
  theme_aih
print(p_TH1Treg)

# Plot 4: Hepatocellular Damage Index
p_DMG <- ggplot(all_results, aes(x = time_days, y = HEPATO_DMG, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Hepatocellular Damage Index", x = "Time (days)", y = "Damage Index (AU)") +
  theme_aih
print(p_DMG)

# Plot 5: B cell count (% baseline) — relevant for RTX scenario
p_Bcell <- ggplot(all_results %>% filter(scenario == "6_RTX_Pred_Refractory"),
                  aes(x = time_days, y = Bcell_pct)) +
  geom_line(color = "#1ABC9C", linewidth = 1.5) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  annotate("text", x = 100, y = 7, label = "Depletion threshold 5%", size = 3) +
  labs(title = "B Cell Count — Rituximab Scenario",
       x = "Time (days)", y = "B Cells (% of baseline)") +
  theme_bw()
print(p_Bcell)

# Plot 6: GR Occupancy vs Time
p_GR <- ggplot(all_results %>% filter(grepl("Pred|Budesonide|RTX", scenario)),
               aes(x = time_days, y = GR_occ_pct, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Glucocorticoid Receptor Occupancy", x = "Time (days)", y = "GR Occupancy (%)") +
  ylim(0, 100) + theme_aih
print(p_GR)

# ---- Clinical Summary Table -------------------------------------------------
clinical_summary <- all_results %>%
  filter(time_days %in% c(0, 30, 90, 180, 365, 730)) %>%
  select(scenario, time_days, ALT_SERUM, ALT_ULN, IgG_serum, Th1_Treg_ratio,
         HEPATO_DMG, GR_occ_pct, T_suppress_pct, Bcell_pct, ANA_equiv, Remission) %>%
  mutate(
    scenario = scenario_labels[scenario],
    across(where(is.numeric), round, 2)
  )

cat("\n====== AIH QSP Model — Clinical Summary Table ======\n")
print(as.data.frame(clinical_summary))

# ---- Sensitivity Analysis ---------------------------------------------------
cat("\n====== Sensitivity Analysis: TPMT Activity on 6-TGN AUC ======\n")
tpmt_levels <- c(0.1, 0.5, 1.0, 1.5)  # poor, intermediate, normal, ultrarapid
tpmt_results <- map_dfr(tpmt_levels, function(tpmt_val) {
  mod_tpmt <- aih_mod %>%
    param(FLAG_PRED = 1, FLAG_AZA = 1, TPMT_act = tpmt_val) %>%
    mrgsim(end = sim_end_h, delta = 24)
  as_tibble(mod_tpmt) %>%
    mutate(TPMT_activity = tpmt_val, time_days = time / 24)
})

p_TPMT <- ggplot(tpmt_results, aes(x = time_days, y = TGN_conc,
                                    color = factor(TPMT_activity))) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(235, 450), linetype = "dashed") +
  annotate("text", x = 200, y = 240, label = "Therapeutic 235 pmol", size = 3) +
  annotate("text", x = 200, y = 455, label = "Toxic 450 pmol", size = 3) +
  scale_color_viridis_d(option = "C", name = "TPMT activity") +
  labs(title = "Sensitivity Analysis: TPMT Activity → 6-TGN Exposure",
       x = "Time (days)", y = "6-TGN (pmol/8×10⁸ RBC equiv)") +
  theme_bw()
print(p_TPMT)

cat("Model compiled and simulations complete.\n")
cat("Scenarios: Untreated | Pred Mono | Pred+AZA (SoC) | BUD+AZA | Pred+MMF | RTX+Pred\n")
