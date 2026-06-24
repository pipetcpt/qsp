## =============================================================================
## ATTR Amyloidosis QSP Model — mrgsolve ODE Implementation
## Transthyretin Amyloidosis: Protein Misfolding, Organ Damage & Drug PK/PD
## =============================================================================
## Disease:  Transthyretin Amyloidosis (ATTR)
##           — ATTRwt (wild-type, cardiac predominant)
##           — ATTRv  (hereditary, neuropathic/cardiac)
## Drugs:    Tafamidis (TTR stabilizer, oral, QD)
##           Vutrisiran (GalNAc-siRNA, SC, Q3M)
##           Inotersen  (2'-MOE ASO, SC, QW)
##           Patisiran  (LNP-siRNA, IV, Q3W)
## ODEs:     25 compartments
## Scenarios: 7 treatment scenarios
## References: ATTR-ACT (Maurer 2018 NEJM), APOLLO (Adams 2018 NEJM),
##             HELIOS-A (Gillmore 2021 NEJM), NEURO-TTR (Benson 2018 Lancet)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE BLOCK
## ─────────────────────────────────────────────────────────────────────────────
code <- '
$PROB
ATTR Amyloidosis QSP Model
Protein misfolding cascade + organ damage + PK of 4 drug classes
25-compartment ODE system

$PARAM @annotated
// ── Patient characteristics ──────────────────────────────────────
ATTR_type : 1   : Disease type (1=ATTRwt-CM, 2=ATTRv neuropathic)
WT         : 75 : Body weight (kg)
AGE        : 72 : Age at diagnosis (years)

// ── TTR synthesis / mRNA dynamics ─────────────────────────────────
kin_mRNA   : 0.01440 : TTR mRNA synthesis rate constant (h^-1, t1/2~48h)
kout_mRNA  : 0.01440 : TTR mRNA degradation rate (h^-1)

// ── Tafamidis PK (oral, 61mg QD) ─────────────────────────────────
ka_TAF     : 0.42   : Absorption rate constant TAF (h^-1)
CL_TAF     : 0.96   : Tafamidis clearance (L/h)
Vc_TAF     : 16.0   : Central volume TAF (L)
Vp_TAF     : 32.0   : Peripheral volume TAF (L)
Q_TAF      : 1.20   : Inter-compartment CL TAF (L/h)
F_TAF      : 0.99   : Tafamidis oral bioavailability
// Tafamidis PD
Emax_stab  : 0.80   : Maximum TTR stabilization fraction (Emax)
EC50_stab  : 0.80   : EC50 for TTR stabilization (ug/mL)
n_stab     : 1.50   : Hill coefficient stabilization

// ── Vutrisiran PK (SC, 25mg Q3M) ─────────────────────────────────
ka_VUT     : 0.080  : SC absorption rate vutrisiran (h^-1)
CL_VUT     : 0.120  : Vutrisiran clearance (L/h)
Vc_VUT     : 5.80   : Central volume vutrisiran (L)
F_VUT      : 0.820  : SC bioavailability vutrisiran
// Vutrisiran PD (RISC-mediated TTR mRNA knockdown)
Imax_VUT   : 0.830  : Maximum mRNA inhibition by vutrisiran
IC50_VUT   : 0.450  : IC50 vutrisiran on TTR mRNA (ng/mL)

// ── Inotersen PK (SC, 300mg QW) ──────────────────────────────────
ka_INO     : 0.100  : SC absorption rate inotersen (h^-1)
CL_INO     : 0.040  : Inotersen clearance (L/h)
Vc_INO     : 78.0   : Central volume inotersen (L) [extensive tissue dist]
F_INO      : 0.700  : SC bioavailability inotersen
// Inotersen PD (RNase H1 mediated)
Imax_INO   : 0.720  : Maximum mRNA inhibition by inotersen
IC50_INO   : 28.0   : IC50 inotersen on TTR mRNA (ng/mL)

// ── Patisiran PK (IV, 0.3 mg/kg Q3W) ─────────────────────────────
CL_PAT     : 0.180  : Patisiran clearance (L/h)
Vc_PAT     : 3.30   : Central volume patisiran (L)
Vp_PAT     : 5.60   : Peripheral volume patisiran (L)
Q_PAT      : 0.960  : Inter-compartment CL patisiran (L/h)
// Patisiran PD (RISC/Ago2 mediated)
Imax_PAT   : 0.800  : Maximum mRNA inhibition by patisiran
IC50_PAT   : 0.600  : IC50 patisiran on TTR mRNA (ng/mL)

// ── TTR misfolding / aggregation ──────────────────────────────────
ksyn_TET   : 0.180  : TTR tetramer synthesis rate (AU/h, normalized)
kout_TET   : 0.0098 : TTR tetramer clearance (h^-1, t1/2~72h plasma)
kdis_TET   : 0.020  : Tetramer dissociation rate (h^-1)
kagg_MONO  : 0.100  : Monomer aggregation → oligomers (h^-1)
kdeg_MONO  : 0.150  : Misfolded monomer proteolytic degradation (h^-1)
kfib_OLIGO : 0.050  : Oligomer → fibril elongation rate (h^-1)
kdeg_OLIGO : 0.030  : Oligomer degradation (h^-1)
kdeg_FIB   : 0.0001 : Fibril degradation (h^-1, very slow!)

// ── Tissue deposition fractions ───────────────────────────────────
f_heart_wt : 0.400  : Fraction of fibrils depositing in heart (ATTRwt)
f_nerve_wt : 0.050  : Fraction depositing in PNS (ATTRwt)
f_heart_v  : 0.200  : Fraction depositing in heart (ATTRv)
f_nerve_v  : 0.350  : Fraction depositing in PNS (ATTRv)

// ── Cardiac PD ────────────────────────────────────────────────────
kdet_EF    : 0.00015 : LVEF deterioration rate per fibril·inflam unit (h^-1)
krec_EF    : 0.00050 : LVEF recovery rate constant (h^-1)
LVEF_base  : 62.0    : Baseline LVEF (%) [ATTRwt ~60-65%]
LVEF_min   : 20.0    : Minimum LVEF (%)
kin_BNP    : 0.500   : NT-proBNP production rate (pg/mL/h per unit)
kout_BNP   : 0.0800  : NT-proBNP elimination rate (h^-1, t1/2~8h)
BNP_base   : 200.0   : Baseline NT-proBNP (pg/mL) [normal <300 pg/mL]
kin_inflam : 0.0100  : Cardiac inflammation production rate (h^-1)
kout_inflam: 0.0500  : Cardiac inflammation resolution rate (h^-1)
FIB50_inf  : 1.50    : FIB_HRT at half-max inflammation induction (AU)

// ── Neurological PD ───────────────────────────────────────────────
kin_NIS    : 0.00200 : NIS progression rate (score/h per fibril unit)
kout_NIS   : 0.000050: NIS regression rate (h^-1, very slow spontaneous)
NIS_base   : 5.0     : Baseline NIS score (mild neuropathy start)
kdet_mBMI  : 0.000040: mBMI decline rate per NIS unit (h^-1)
NIS50_mBMI : 50.0    : NIS at half-max mBMI decline
mBMI_base  : 1000.0  : Baseline mBMI (g/m^2; BMI 25 × albumin 40 g/L)

// ── Renal PD ──────────────────────────────────────────────────────
kdet_eGFR  : 0.000020: eGFR decline rate per systemic fibril (h^-1)
eGFR_base  : 72.0    : Baseline eGFR (mL/min/1.73m^2)

$OMEGA @labels ETACL_TAF ETAVC_TAF ETACL_VUT ETACL_INO
0.09 0.04 0.09 0.09

$SIGMA
0.04

$MAIN
// ── Active fraction ──────────────────────────────────────────────
double f_heart = (ATTR_type == 1) ? f_heart_wt : f_heart_v ;
double f_nerve = (ATTR_type == 1) ? f_nerve_wt : f_nerve_v ;
double f_sys   = 1.0 - f_heart - f_nerve ;
if(f_sys < 0.05) f_sys = 0.05 ;

// ── Baseline tetramer ─────────────────────────────────────────────
// At SS: ksyn_TET * mRNA_ss = (kdis_TET + kout_TET) * TET_ss
double TET_ss  = ksyn_TET / (kdis_TET + kout_TET) ;
double MONO_ss = (2.0 * kdis_TET * TET_ss) / (kagg_MONO + kdeg_MONO) ;

// ── IIV ───────────────────────────────────────────────────────────
double CL_TAF_i  = CL_TAF  * exp(ETACL_TAF) ;
double Vc_TAF_i  = Vc_TAF  * exp(ETAVC_TAF) ;
double CL_VUT_i  = CL_VUT  * exp(ETACL_VUT) ;
double CL_INO_i  = CL_INO  * exp(ETACL_INO) ;

// ── Initial conditions ────────────────────────────────────────────
double A_TAF_GUT_0 = 0 ;
double A_TAF_C_0   = 0 ;
double A_TAF_P_0   = 0 ;
double A_VUT_SC_0  = 0 ;
double A_VUT_C_0   = 0 ;
double A_INO_SC_0  = 0 ;
double A_INO_C_0   = 0 ;
double A_PAT_C_0   = 0 ;
double A_PAT_P_0   = 0 ;
double TTR_MRNA_0  = 1.0 ;
double TTR_TET_0   = TET_ss ;
double TTR_MONO_0  = MONO_ss ;
double TTR_OLIGO_0 = 0.01 ;
double FIB_HRT_0   = 0.10 ;  // small initial cardiac load (pre-clinical)
double FIB_NRV_0   = 0.05 ;
double FIB_SYS_0   = 0.05 ;
double INFLAM_0    = kin_inflam * (FIB_HRT_0/(FIB50_inf+FIB_HRT_0)) / kout_inflam ;
double LVEF_0      = LVEF_base ;
double NT_proBNP_0 = BNP_base ;
double NIS_0       = NIS_base ;
double mBMI_0      = mBMI_base ;
double eGFR_0      = eGFR_base ;

$INIT
A_TAF_GUT  = A_TAF_GUT_0,
A_TAF_C    = A_TAF_C_0,
A_TAF_P    = A_TAF_P_0,
A_VUT_SC   = A_VUT_SC_0,
A_VUT_C    = A_VUT_C_0,
A_INO_SC   = A_INO_SC_0,
A_INO_C    = A_INO_C_0,
A_PAT_C    = A_PAT_C_0,
A_PAT_P    = A_PAT_P_0,
TTR_MRNA   = TTR_MRNA_0,
TTR_TET    = TET_ss,
TTR_MONO   = MONO_ss,
TTR_OLIGO  = TTR_OLIGO_0,
FIB_HRT    = FIB_HRT_0,
FIB_NRV    = FIB_NRV_0,
FIB_SYS    = FIB_SYS_0,
INFLAM     = INFLAM_0,
LVEF       = LVEF_0,
NT_proBNP  = NT_proBNP_0,
NIS        = NIS_0,
mBMI       = mBMI_0,
eGFR       = eGFR_0,
SYMP_CARD  = 0,
SYMP_NEURO = 0

$ODE
// ═══════════════════════════════════════════════════════════════
// DRUG PK
// ═══════════════════════════════════════════════════════════════

// ── Tafamidis PK (oral 2-cpt) ────────────────────────────────
dxdt_A_TAF_GUT = -ka_TAF * A_TAF_GUT ;
double C_TAF   = A_TAF_C / Vc_TAF_i ;
dxdt_A_TAF_C   = F_TAF * ka_TAF * A_TAF_GUT
                 + Q_TAF * (A_TAF_P/Vp_TAF - C_TAF)
                 - (CL_TAF_i/Vc_TAF_i) * A_TAF_C ;
dxdt_A_TAF_P   = Q_TAF * (C_TAF - A_TAF_P/Vp_TAF) ;

// ── Vutrisiran PK (SC 1-cpt) ─────────────────────────────────
dxdt_A_VUT_SC  = -ka_VUT * A_VUT_SC ;
double C_VUT   = A_VUT_C / Vc_VUT ;
dxdt_A_VUT_C   = F_VUT * ka_VUT * A_VUT_SC
                 - (CL_VUT_i/Vc_VUT) * A_VUT_C ;

// ── Inotersen PK (SC 1-cpt) ──────────────────────────────────
dxdt_A_INO_SC  = -ka_INO * A_INO_SC ;
double C_INO   = A_INO_C / Vc_INO ;
dxdt_A_INO_C   = F_INO * ka_INO * A_INO_SC
                 - (CL_INO_i/Vc_INO) * A_INO_C ;

// ── Patisiran PK (IV 2-cpt; input via R_PAT dosing) ──────────
double C_PAT   = A_PAT_C / Vc_PAT ;
dxdt_A_PAT_C   = R_PAT                             // zero-order IV infusion
                 + Q_PAT * (A_PAT_P/Vp_PAT - C_PAT)
                 - (CL_PAT/Vc_PAT) * A_PAT_C ;
dxdt_A_PAT_P   = Q_PAT * (C_PAT - A_PAT_P/Vp_PAT) ;

// ═══════════════════════════════════════════════════════════════
// DRUG EFFECTS
// ═══════════════════════════════════════════════════════════════

// TTR stabilization (tafamidis): reduces tetramer dissociation
double E_stab = Emax_stab * pow(C_TAF, n_stab) /
                (pow(EC50_stab, n_stab) + pow(C_TAF, n_stab)) ;
if(E_stab > 0.95) E_stab = 0.95 ;

// RNA-based: TTR mRNA knockdown
double E_VUT  = Imax_VUT * C_VUT / (IC50_VUT + C_VUT) ;
double E_PAT  = Imax_PAT * C_PAT / (IC50_PAT + C_PAT) ;
double E_INO  = Imax_INO * C_INO / (IC50_INO + C_INO) ;
double E_RNA  = E_VUT + E_PAT + E_INO ;
if(E_RNA > 0.95) E_RNA = 0.95 ;

// ═══════════════════════════════════════════════════════════════
// DISEASE PD — TTR misfolding cascade
// ═══════════════════════════════════════════════════════════════

// TTR mRNA (normalized, 1 = normal)
dxdt_TTR_MRNA  = kin_mRNA * (1.0 - E_RNA) - kout_mRNA * TTR_MRNA ;

// TTR tetramer (plasma level, AU)
// Synthesis ∝ mRNA; dissociation inhibited by stabilizer
dxdt_TTR_TET   = ksyn_TET * TTR_MRNA
                 - kdis_TET * (1.0 - E_stab) * TTR_TET
                 - kout_TET * TTR_TET ;

// Misfolded monomer (from dissociation; 2 monomers per tetramer split)
dxdt_TTR_MONO  = 2.0 * kdis_TET * (1.0 - E_stab) * TTR_TET
                 - kagg_MONO  * TTR_MONO
                 - kdeg_MONO  * TTR_MONO ;

// Oligomers (toxic intermediates)
dxdt_TTR_OLIGO = kagg_MONO  * TTR_MONO
                 - kfib_OLIGO * TTR_OLIGO
                 - kdeg_OLIGO * TTR_OLIGO ;

// Mature fibrils in heart (very slow clearance)
dxdt_FIB_HRT   = kfib_OLIGO * TTR_OLIGO * f_heart
                 - kdeg_FIB * FIB_HRT ;

// Mature fibrils in peripheral nerves
dxdt_FIB_NRV   = kfib_OLIGO * TTR_OLIGO * f_nerve
                 - kdeg_FIB * FIB_NRV ;

// Systemic fibrils (kidneys, spleen, GI, etc.)
dxdt_FIB_SYS   = kfib_OLIGO * TTR_OLIGO * f_sys
                 - kdeg_FIB * FIB_SYS ;

// ═══════════════════════════════════════════════════════════════
// DISEASE PD — Cardiac
// ═══════════════════════════════════════════════════════════════

// Cardiac inflammation (macrophage, IL-1β, NLRP3-driven)
dxdt_INFLAM    = kin_inflam * FIB_HRT / (FIB50_inf + FIB_HRT)
                 - kout_inflam * INFLAM ;

// LVEF (%) — deterioration by fibrils×inflammation; recovery if load reduces
double EF_curr = LVEF ;
double EF_det  = kdet_EF * FIB_HRT * (1.0 + INFLAM) * EF_curr ;
double EF_rec  = krec_EF * (LVEF_base - EF_curr) ;
dxdt_LVEF      = EF_rec - EF_det ;
if(LVEF < LVEF_min) dxdt_LVEF = 0 ;

// NT-proBNP (pg/mL) — wall stress + inflammation + inverse LVEF
double NT_in   = kin_BNP * (INFLAM + 1.0/(EF_curr/100.0 + 0.01)) ;
dxdt_NT_proBNP = NT_in - kout_BNP * NT_proBNP ;

// ═══════════════════════════════════════════════════════════════
// DISEASE PD — Neurological
// ═══════════════════════════════════════════════════════════════

// NIS (0–244): progressive nerve damage from fibril deposition
dxdt_NIS       = kin_NIS * FIB_NRV
                 - kout_NIS * NIS ;

// mBMI (g/m^2): nutrition decline driven by NIS (GI dysmotility)
dxdt_mBMI      = -kdet_mBMI * (NIS/(NIS + NIS50_mBMI)) * mBMI ;

// ═══════════════════════════════════════════════════════════════
// DISEASE PD — Renal
// ═══════════════════════════════════════════════════════════════
dxdt_eGFR      = -kdet_eGFR * FIB_SYS * eGFR ;

// ═══════════════════════════════════════════════════════════════
// COMPOSITE SYMPTOM SCORES (cumulative)
// ═══════════════════════════════════════════════════════════════
// Cardiac symptom accumulator (integrates NYHA class analogue)
double NYHA_analogue = 1.0 + (1.0 - LVEF/LVEF_base) * 3.0 ;
dxdt_SYMP_CARD  = NYHA_analogue - SYMP_CARD * 0.001 ;

// Neurological disability accumulator (FAP stage analogue)
double FAP_analogue  = 0.5 + NIS/60.0 ;
dxdt_SYMP_NEURO = FAP_analogue - SYMP_NEURO * 0.001 ;

$TABLE
capture C_TAF_ugmL  = A_TAF_C / Vc_TAF ;
capture C_VUT_ngmL  = (A_VUT_C / Vc_VUT) * 1000 ;
capture C_INO_ngmL  = (A_INO_C / Vc_INO) * 1000 ;
capture C_PAT_ugmL  = A_PAT_C / Vc_PAT ;
capture E_stab_frac = Emax_stab * pow(C_TAF_ugmL, n_stab) /
                      (pow(EC50_stab, n_stab) + pow(C_TAF_ugmL, n_stab)) ;
capture E_RNA_frac  = E_VUT + E_PAT + E_INO ;
capture TTR_MRNA_rel = TTR_MRNA ;
capture TTR_TET_AU   = TTR_TET ;
capture TTR_OLIGO_AU = TTR_OLIGO ;
capture FIB_HRT_AU   = FIB_HRT ;
capture FIB_NRV_AU   = FIB_NRV ;
capture LVEF_pct     = LVEF ;
capture NTproBNP     = NT_proBNP ;
capture NIS_score    = NIS ;
capture mBMI_val     = mBMI ;
capture eGFR_val     = eGFR ;
capture INFLAM_lvl   = INFLAM ;

$CAPTURE C_TAF_ugmL C_VUT_ngmL C_INO_ngmL C_PAT_ugmL
         E_stab_frac E_RNA_frac
         TTR_MRNA_rel TTR_TET_AU TTR_OLIGO_AU
         FIB_HRT_AU FIB_NRV_AU
         LVEF_pct NTproBNP NIS_score mBMI_val eGFR_val INFLAM_lvl
'

## ─────────────────────────────────────────────────────────────────────────────
## Compile model
## ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("ATTR_QSP", code)

## ─────────────────────────────────────────────────────────────────────────────
## Treatment Scenario Definitions
## ─────────────────────────────────────────────────────────────────────────────
SIM_DURATION <- 365 * 3   # 3-year simulation (hours = 26280)
SIM_DELTA    <- 6         # Output every 6 hours

## S1: Natural history — ATTRwt (cardiac predominant)
ev_S1_wt <- ev(time = 0, amt = 0, cmt = 1, ii = 0, addl = 0)

## S2: Natural history — ATTRv (neuropathic)
ev_S2_v  <- ev(time = 0, amt = 0, cmt = 1, ii = 0, addl = 0)

## S3: Tafamidis 61 mg QD oral (ATTR-ACT: ATTRwt-CM)
# cmt=1 = A_TAF_GUT; amt in mg → μg for concentration in L
# C(μg/mL) = A(mg)/V(L) × (1000 μg/mg) / (1000 mL/L) = A/V [same]
ev_S3_taf <- ev(time = 0, amt = 61, cmt = 1, ii = 24, addl = 1094)  # 3yr daily

## S4: Patisiran 0.3 mg/kg IV Q3W (APOLLO: ATTRv polyneuropathy)
# Infusion: 60-min (0.017 h^-1 fraction? use R_PAT as cmt infusion)
# Patisiran central cmt = A_PAT_C (cmt=8), IV bolus approximated
# Use tinf=1h: amt / tinf → rate
patisiran_mg <- 0.3 * 75      # dose = 22.5 mg for 75 kg
ev_S4_pat <- ev(time = 0, amt = patisiran_mg, cmt = 8, rate = -2,
                ii = 504, addl = 155)  # Q3W = 504h, 3yr

## S5: Vutrisiran 25 mg SC Q3M (HELIOS-A: ATTRv polyneuropathy)
# 1 injection every ~91 days = 2184h
ev_S5_vut <- ev(time = 0, amt = 25, cmt = 4,
                ii = 2184, addl = 11)  # Q3M × 12 = 3yr

## S6: Inotersen 300 mg SC QW (NEURO-TTR: ATTRv polyneuropathy)
ev_S6_ino <- ev(time = 0, amt = 300, cmt = 6,
                ii = 168, addl = 155)  # QW × 156 weeks = 3yr

## S7: Tafamidis + Vutrisiran combination (hypothetical combination therapy)
ev_S7_comb <- ev_S3_taf + ev_S5_vut

## ─────────────────────────────────────────────────────────────────────────────
## Run Simulations
## ─────────────────────────────────────────────────────────────────────────────

run_scenario <- function(mod, ev_obj, attr_type = 1,
                         label = "Scenario", n_subj = 1) {
  out <- mod %>%
    param(ATTR_type = attr_type) %>%
    ev(ev_obj) %>%
    mrgsim(end = SIM_DURATION, delta = SIM_DELTA, nid = n_subj) %>%
    as_tibble() %>%
    mutate(time_yr = time / 8760,
           Scenario = label)
  return(out)
}

message("Running 7 treatment scenarios (3-year simulation each)...")
out_S1 <- run_scenario(mod, ev_S1_wt, attr_type=1, label="S1: ATTRwt Natural History")
out_S2 <- run_scenario(mod, ev_S2_v,  attr_type=2, label="S2: ATTRv Natural History")
out_S3 <- run_scenario(mod, ev_S3_taf, attr_type=1, label="S3: Tafamidis 61mg QD")
out_S4 <- run_scenario(mod, ev_S4_pat, attr_type=2, label="S4: Patisiran 0.3mg/kg Q3W")
out_S5 <- run_scenario(mod, ev_S5_vut, attr_type=2, label="S5: Vutrisiran 25mg Q3M")
out_S6 <- run_scenario(mod, ev_S6_ino, attr_type=2, label="S6: Inotersen 300mg QW")
out_S7 <- run_scenario(mod, ev_S7_comb, attr_type=1, label="S7: Tafamidis+Vutrisiran (combo)")

all_out <- bind_rows(out_S1, out_S2, out_S3, out_S4, out_S5, out_S6, out_S7)

## ─────────────────────────────────────────────────────────────────────────────
## Visualization
## ─────────────────────────────────────────────────────────────────────────────
scen_colors <- c(
  "S1: ATTRwt Natural History"       = "#e74c3c",
  "S2: ATTRv Natural History"        = "#c0392b",
  "S3: Tafamidis 61mg QD"            = "#2980b9",
  "S4: Patisiran 0.3mg/kg Q3W"       = "#27ae60",
  "S5: Vutrisiran 25mg Q3M"          = "#9b59b6",
  "S6: Inotersen 300mg QW"           = "#f39c12",
  "S7: Tafamidis+Vutrisiran (combo)" = "#1abc9c"
)

# Plot 1: Drug PK profiles (Year 1, first 90 days)
pk_data <- all_out %>% filter(time <= 90*24, Scenario %in%
  c("S3: Tafamidis 61mg QD","S4: Patisiran 0.3mg/kg Q3W",
    "S5: Vutrisiran 25mg Q3M","S6: Inotersen 300mg QW"))

p_taf <- pk_data %>% filter(Scenario == "S3: Tafamidis 61mg QD") %>%
  ggplot(aes(time/24, C_TAF_ugmL)) +
  geom_line(color="#2980b9", size=1.2) +
  geom_hline(yintercept=0.8, linetype="dashed", color="gray50") +
  annotate("text", x=60, y=1.0, label="EC50=0.8 μg/mL", size=3, color="gray40") +
  labs(x="Days", y="Conc (μg/mL)", title="Tafamidis PK (61mg QD)") +
  theme_bw()

p_vut <- pk_data %>% filter(Scenario == "S5: Vutrisiran 25mg Q3M") %>%
  ggplot(aes(time/24, C_VUT_ngmL)) +
  geom_line(color="#9b59b6", size=1.2) +
  labs(x="Days", y="Conc (ng/mL)", title="Vutrisiran PK (25mg Q3M)") +
  theme_bw()

p_ino <- pk_data %>% filter(Scenario == "S6: Inotersen 300mg QW") %>%
  ggplot(aes(time/24, C_INO_ngmL)) +
  geom_line(color="#f39c12", size=1.2) +
  geom_hline(yintercept=28, linetype="dashed", color="gray50") +
  annotate("text", x=60, y=35, label="IC50=28 ng/mL", size=3, color="gray40") +
  labs(x="Days", y="Conc (ng/mL)", title="Inotersen PK (300mg QW)") +
  theme_bw()

p_pat <- pk_data %>% filter(Scenario == "S4: Patisiran 0.3mg/kg Q3W") %>%
  ggplot(aes(time/24, C_PAT_ugmL)) +
  geom_line(color="#27ae60", size=1.2) +
  labs(x="Days", y="Conc (μg/mL)", title="Patisiran PK (0.3mg/kg Q3W)") +
  theme_bw()

p_pk <- (p_taf + p_vut) / (p_ino + p_pat) +
  plot_annotation(title="Drug PK Profiles — First 90 Days",
                  theme=theme(plot.title=element_text(size=14, face="bold")))

# Plot 2: TTR mRNA knockdown
p_mrna <- all_out %>%
  filter(Scenario %in% c("S1: ATTRwt Natural History",
                          "S4: Patisiran 0.3mg/kg Q3W",
                          "S5: Vutrisiran 25mg Q3M",
                          "S6: Inotersen 300mg QW")) %>%
  ggplot(aes(time_yr, TTR_MRNA_rel * 100, color=Scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=scen_colors) +
  geom_hline(yintercept=c(20, 30), linetype="dashed", color="gray60") +
  annotate("text", x=2.5, y=25, label="~80% knockdown target (siRNA)", size=3.5) +
  labs(x="Time (years)", y="TTR mRNA (% baseline)",
       title="TTR mRNA Expression — RNA-based Therapies",
       color="Scenario") +
  scale_x_continuous(breaks=0:3) +
  ylim(0, 105) +
  theme_bw() + theme(legend.position="bottom")

# Plot 3: Cardiac outcomes
p_cardiac <- all_out %>%
  filter(Scenario %in% c("S1: ATTRwt Natural History",
                          "S3: Tafamidis 61mg QD",
                          "S7: Tafamidis+Vutrisiran (combo)")) %>%
  select(time_yr, Scenario, LVEF_pct, NTproBNP, FIB_HRT_AU) %>%
  pivot_longer(c(LVEF_pct, NTproBNP, FIB_HRT_AU), names_to="endpoint") %>%
  mutate(endpoint_lab = case_when(
    endpoint == "LVEF_pct"   ~ "LVEF (%)",
    endpoint == "NTproBNP"   ~ "NT-proBNP (pg/mL)",
    endpoint == "FIB_HRT_AU" ~ "Cardiac Fibril Load (AU)"
  )) %>%
  ggplot(aes(time_yr, value, color=Scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=scen_colors) +
  facet_wrap(~endpoint_lab, scales="free_y", ncol=1) +
  labs(x="Time (years)", y="", title="Cardiac Outcomes — ATTRwt",
       color="Scenario") +
  scale_x_continuous(breaks=0:3) +
  theme_bw() + theme(legend.position="bottom")

# Plot 4: Neurological outcomes (ATTRv)
p_neuro <- all_out %>%
  filter(Scenario %in% c("S2: ATTRv Natural History",
                          "S4: Patisiran 0.3mg/kg Q3W",
                          "S5: Vutrisiran 25mg Q3M",
                          "S6: Inotersen 300mg QW")) %>%
  select(time_yr, Scenario, NIS_score, mBMI_val, FIB_NRV_AU) %>%
  pivot_longer(c(NIS_score, mBMI_val, FIB_NRV_AU), names_to="endpoint") %>%
  mutate(endpoint_lab = case_when(
    endpoint == "NIS_score"  ~ "NIS Score (0-244)",
    endpoint == "mBMI_val"   ~ "Modified BMI (g/m²)",
    endpoint == "FIB_NRV_AU" ~ "Nerve Fibril Load (AU)"
  )) %>%
  ggplot(aes(time_yr, value, color=Scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=scen_colors) +
  facet_wrap(~endpoint_lab, scales="free_y", ncol=1) +
  labs(x="Time (years)", y="", title="Neurological Outcomes — ATTRv",
       color="Scenario") +
  scale_x_continuous(breaks=0:3) +
  theme_bw() + theme(legend.position="bottom")

# Plot 5: Scenario comparison at 18 months
compare_18mo <- all_out %>%
  filter(abs(time - 18*730.5) == min(abs(time - 18*730.5))) %>%
  group_by(Scenario) %>%
  summarise(
    LVEF_pct     = mean(LVEF_pct),
    NTproBNP     = mean(NTproBNP),
    NIS_score    = mean(NIS_score),
    mBMI_val     = mean(mBMI_val),
    TTR_MRNA_pct = mean(TTR_MRNA_rel * 100),
    FIB_HRT_AU   = mean(FIB_HRT_AU),
    FIB_NRV_AU   = mean(FIB_NRV_AU),
    eGFR_val     = mean(eGFR_val),
    .groups = "drop"
  )

message("\n── 18-Month Summary (Simulation Output) ────────────────────────")
print(compare_18mo, n=7)
message("────────────────────────────────────────────────────────────────")

# Print all plots
print(p_pk)
print(p_mrna)
print(p_cardiac)
print(p_neuro)

message("\n✓ ATTR QSP simulation complete.")
message("  7 scenarios × 3 years × 6h delta = ", nrow(all_out)/7, " timepoints/scenario")
message("  Key results at 18-months printed above.")
