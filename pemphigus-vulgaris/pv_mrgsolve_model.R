## =============================================================================
## Pemphigus Vulgaris (PV) — Quantitative Systems Pharmacology Model
## mrgsolve ODE-based PK/PD simulation
##
## Disease: Autoimmune blistering disease driven by anti-Dsg3/Dsg1 IgG
## Key drugs: Prednisone, Rituximab (RTX), Mycophenolate Mofetil (MMF),
##            Azathioprine (AZA), IVIg, Efgartigimod (FcRn blocker)
##
## Compartments (18 ODE states):
##   Drug PK  : Cp1_pred, Cp2_pred, CR1, CR2, CMPA
##   Immune   : BN, BGC, BM, SLPC, LLPC, Tfh, Treg
##   Antibody : Ab3, Ab1
##   Disease  : Dsg3_loss, PDAI_state, Comp_act, Cort_bone
##
## References:
##   - Joly et al. NEJM 2017 (RITUX3 trial)
##   - Meijer et al. JACI 2016 (PK of rituximab in PV)
##   - Murrell et al. JAMA Derm 2018 (MMF trial)
##   - Amber et al. BJD 2019 (PV disease models)
##   - Colliou et al. Sci Transl Med 2013 (Dsg3 tolerance)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─────────────────────────────────────────────────────────────────────────────
## 1. MODEL CODE
## ─────────────────────────────────────────────────────────────────────────────

pv_model_code <- '
$PROB Pemphigus Vulgaris QSP model — mrgsolve

$PARAM @annotated
// ── Prednisolone PK (2-compartment, oral) ────────────────────────────────
Cl_pred  : 5.10 : Prednisolone clearance (L/h) [Pickup 1977, Czock 2005]
Vc_pred  : 35.0 : Central volume of distribution (L)
Q_pred   : 2.50 : Inter-compartmental clearance (L/h)
Vp_pred  : 60.0 : Peripheral volume (L)
ka_pred  : 1.20 : Absorption rate constant (h-1, F~0.82)
F_pred   : 0.82 : Oral bioavailability

// ── Rituximab PK (2-compartment, IV, target-mediated) ────────────────────
Cl_RTX   : 0.014 : RTX linear clearance (L/h) → t1/2 ~22 days [Meijer 2016]
Vc_RTX   : 3.10  : RTX central volume (L)
Q_RTX    : 0.18  : RTX inter-compartmental CL (L/h)
Vp_RTX   : 4.40  : RTX peripheral volume (L)
kTMDD    : 0.002 : TMDD binding rate (L/h per B-cell unit) [Reff 1994]

// ── MPA (from MMF) PK (1-compartment) ────────────────────────────────────
Cl_MPA   : 15.0  : MPA apparent clearance (L/h) [Bullingham 1998]
Vc_MPA   : 50.0  : MPA volume of distribution (L)
ka_MPA   : 0.80  : MMF absorption rate (h-1)
F_MPA    : 0.94  : MMF de-esterification bioavailability

// ── IVIg PK ──────────────────────────────────────────────────────────────
Cl_IVIg  : 0.005 : IVIg clearance baseline (L/h) [Bleeker 2001]
Vc_IVIg  : 3.50  : IVIg volume (L)

// ── B cell dynamics ───────────────────────────────────────────────────────
kBN_prod : 0.015  : Naive B cell production (cells/mL/h × 10^-6)
kBN_die  : 0.004  : Naive B cell natural death rate (h-1)
kBN_act  : 0.0003 : Naive B→GC activation rate (h-1 per Tfh unit)
kBGC_die : 0.010  : GC B cell apoptosis rate (h-1)
kBM_form : 0.008  : Memory B cell formation from GC (h-1)
kBM_die  : 0.0004 : Memory B cell death rate (h-1; t1/2 ~70 days)
kSLPC_form: 0.012 : Short-lived PC formation rate (h-1)
kSLPC_die: 0.040  : SLPC death rate (h-1; t1/2 ~17 h)
kLLPC_form: 0.002 : Long-lived PC formation from memory B (h-1)
kLLPC_die: 0.0003 : LLPC death rate (h-1; t1/2 ~months)

// ── Tfh/Treg dynamics ────────────────────────────────────────────────────
kTfh_prod: 0.002  : Tfh production rate (h-1)
kTfh_die : 0.008  : Tfh natural death (h-1)
kTreg_prod:0.001  : Treg production rate (h-1)
kTreg_die: 0.005  : Treg death rate (h-1)

// ── Anti-Dsg antibody kinetics ────────────────────────────────────────────
kAb3_prod: 0.0006 : Anti-Dsg3 IgG secretion per PC (U/mL/h per cell unit)
kAb3_die : 0.0014 : Anti-Dsg3 IgG degradation rate (h-1; t1/2 ~20 days)
kAb1_prod: 0.0003 : Anti-Dsg1 IgG secretion rate
kAb1_die : 0.0014 : Anti-Dsg1 IgG degradation rate

// ── Disease dynamics ──────────────────────────────────────────────────────
Dsg3_base : 100.0 : Baseline Dsg3 protein expression (arb units)
kDsg3_int : 0.0004: Dsg3 internalization rate per unit Ab3 (h-1)
kDsg3_rest: 0.0008: Dsg3 re-expression rate (h-1)
kPDAI_up  : 0.018 : PDAI increase rate (per unit Dsg3 loss fraction)
kPDAI_dn  : 0.004 : PDAI natural resolution rate (h-1)

// ── Complement activation ─────────────────────────────────────────────────
kComp_on  : 0.0003: Complement activation rate per unit IgG1-Ab
kComp_off : 0.010  : Complement resolution rate (h-1)

// ── Corticosteroid bone effect ────────────────────────────────────────────
kBone_loss: 0.0001 : Bone mineral density loss rate per pred concentration
kBone_rest: 0.00005: Natural bone restoration rate (h-1)

// ── Drug PD parameters ───────────────────────────────────────────────────
EC50_pred_T  : 5.0   : Pred EC50 for T cell suppression (mg/L)
Emax_pred_T  : 0.90  : Pred Emax for T cell suppression
EC50_pred_B  : 8.0   : Pred EC50 for B cell/GC suppression (mg/L)
Emax_pred_B  : 0.85  : Pred Emax for B cell suppression
EC50_RTX_B   : 0.05  : RTX EC50 for B cell depletion (mg/L)
Emax_RTX_B   : 0.99  : RTX Emax for B cell depletion
Hill_RTX     : 2.0   : RTX Hill coefficient
EC50_MPA_GC  : 2.0   : MPA EC50 for GC B cell proliferation block (mg/L)
Emax_MPA_GC  : 0.80  : MPA Emax GC block
EC50_IVIg_Ab : 1.0   : IVIg EC50 for FcRn saturation (mg/mL)
Emax_IVIg_Ab : 0.75  : IVIg Emax for Ab catabolism acceleration
EC50_EFG_Ab  : 0.10  : Efgartigimod EC50 for IgG reduction (mg/L)
Emax_EFG_Ab  : 0.90  : Efgartigimod Emax

// ── Baseline steady-state values ─────────────────────────────────────────
BN0    : 50.0  : Baseline naive B cells (x10^6/L)
BGC0   : 5.0   : Baseline GC B cells
BM0    : 20.0  : Baseline memory B cells
SLPC0  : 3.0   : Baseline short-lived plasma cells
LLPC0  : 2.0   : Baseline long-lived plasma cells
Tfh0   : 10.0  : Baseline Tfh cells (arb)
Treg0  : 8.0   : Baseline Treg cells (arb)
Ab3_0  : 180.0 : Initial anti-Dsg3 IgG (U/mL) — active disease
Ab1_0  : 80.0  : Initial anti-Dsg1 IgG (U/mL)

$CMT @annotated
// ── Drug PK (5 compartments) ─────────────────────────────────────────────
DEPOT_pred  : Prednisolone gut depot (mg)
Cp1_pred    : Prednisolone central plasma (mg/L)
Cp2_pred    : Prednisolone peripheral (mg/L)
CR1         : Rituximab central plasma (mg/L)
CR2         : Rituximab peripheral (mg/L)
DEPOT_MPA   : MPA gut depot (mg)
CMPA        : MPA central plasma (mg/L)
DEPOT_IVIg  : IVIg infusion rate depot
CIVG        : IVIg plasma concentration (mg/mL)
CEFG        : Efgartigimod plasma concentration (mg/L)

// ── Immune compartments (7) ───────────────────────────────────────────────
BN    : Naive B cells (x10^6/L)
BGC   : Germinal center B cells (x10^6/L)
BM    : Memory B cells (x10^6/L)
SLPC  : Short-lived plasma cells (x10^6/L)
LLPC  : Long-lived plasma cells (x10^6/L)
Tfh   : Tfh cells (arb units)
Treg  : Treg cells (arb units)

// ── Antibody biomarkers (2) ───────────────────────────────────────────────
Ab3   : Anti-Dsg3 IgG (U/mL) — disease-driving
Ab1   : Anti-Dsg1 IgG (U/mL) — superficial split

// ── Disease state variables (4) ──────────────────────────────────────────
Dsg3_loss  : Dsg3 protein loss (0=normal, ↑=disease)
PDAI_state : Disease Activity Index (0-250)
Comp_act   : Complement activation level (0-100)
Cort_bone  : Cumulative bone effect (steroid-induced)

$MAIN
// ── PD effect functions ───────────────────────────────────────────────────
// Prednisolone effects
double Epred_T = Emax_pred_T * Cp1_pred / (EC50_pred_T + Cp1_pred + 1e-10);
double Epred_B = Emax_pred_B * Cp1_pred / (EC50_pred_B + Cp1_pred + 1e-10);

// Rituximab Hill depletion
double CR1_Hill = pow(CR1, Hill_RTX);
double EC50_Hill = pow(EC50_RTX_B, Hill_RTX);
double ERTX   = Emax_RTX_B * CR1_Hill / (EC50_Hill + CR1_Hill + 1e-20);

// MPA effect on GC proliferation
double EMPA   = Emax_MPA_GC * CMPA / (EC50_MPA_GC + CMPA + 1e-10);

// IVIg effect on IgG catabolism (FcRn saturation)
double EIVIg  = Emax_IVIg_Ab * CIVG / (EC50_IVIg_Ab + CIVG + 1e-10);

// Efgartigimod effect
double EEFG   = Emax_EFG_Ab * CEFG / (EC50_EFG_Ab + CEFG + 1e-10);

// Combined IgG catabolism enhancement (IVIg + efgartigimod)
double E_Ab_cat = fmax(EIVIg, EEFG);

// B cell TMDD term (target-mediated disposition of RTX)
double B_total = BN + BGC + BM;

// Available Dsg3 for internalization
double Dsg3_avail = fmax(0.0, Dsg3_base - Dsg3_loss);

// Treg relative inhibition on B cell activation (Treg0 = homeostasis)
double Treg_suppress = (Treg0 > 0) ? Treg / (Treg + Treg0) : 0.5;

// Tfh stimulation factor
double Tfh_stimul = Tfh / (Tfh0 + 1e-10);

$ODE
// ═══════════════════════════════════════════════════════════════════════════
// DRUG PK
// ═══════════════════════════════════════════════════════════════════════════

// Prednisolone (2-compartment, oral)
dxdt_DEPOT_pred = -ka_pred * DEPOT_pred;
dxdt_Cp1_pred   = (ka_pred * DEPOT_pred * F_pred) / Vc_pred
                  - (Cl_pred / Vc_pred) * Cp1_pred
                  - (Q_pred  / Vc_pred) * Cp1_pred
                  + (Q_pred  / Vp_pred) * Cp2_pred;
dxdt_Cp2_pred   =  (Q_pred  / Vc_pred) * Cp1_pred
                  - (Q_pred  / Vp_pred) * Cp2_pred;

// Rituximab (2-compartment + target binding TMDD)
dxdt_CR1 = -(Cl_RTX / Vc_RTX) * CR1
           - (Q_RTX  / Vc_RTX) * CR1
           + (Q_RTX  / Vp_RTX) * CR2
           - kTMDD * CR1 * B_total;          // target-mediated elimination
dxdt_CR2 =  (Q_RTX  / Vc_RTX) * CR1
           - (Q_RTX  / Vp_RTX) * CR2;

// MPA (1-compartment effective, oral MMF)
dxdt_DEPOT_MPA = -ka_MPA * DEPOT_MPA;
dxdt_CMPA      = (ka_MPA * DEPOT_MPA * F_MPA) / Vc_MPA
                 - (Cl_MPA / Vc_MPA) * CMPA;

// IVIg (IV bolus/infusion, 1-compartment)
dxdt_DEPOT_IVIg = -0.5 * DEPOT_IVIg;        // slow IVIg infusion approximation
dxdt_CIVG       = 0.5 * DEPOT_IVIg / Vc_IVIg
                  - (Cl_IVIg / Vc_IVIg) * CIVG * (1.0 + 3.0 * EIVIg);

// Efgartigimod (IV, simple 1-comp)
dxdt_CEFG = -0.04 * CEFG;                   // t1/2 ~17 h (SC 10mg/kg)

// ═══════════════════════════════════════════════════════════════════════════
// IMMUNE COMPARTMENTS
// ═══════════════════════════════════════════════════════════════════════════

// Naive B cells
dxdt_BN = kBN_prod
         - kBN_die  * BN * (1.0 + Epred_B + ERTX)
         - kBN_act  * BN * Tfh_stimul * (1.0 - Treg_suppress)
                        * (1.0 - Epred_T);

// Germinal Center B cells
dxdt_BGC = kBN_act * BN * Tfh_stimul * (1.0 - Treg_suppress)
                        * (1.0 - Epred_T)
          - kBGC_die  * BGC
          - kBM_form  * BGC * (1.0 - EMPA) * (1.0 - Epred_B)
          - kSLPC_form * BGC * (1.0 - EMPA)
          - ERTX * kBGC_die * BGC;           // RTX-mediated GC B depletion

// Memory B cells
dxdt_BM = kBM_form  * BGC * (1.0 - ERTX)
         - kBM_die   * BM * (1.0 + ERTX)
         - kLLPC_form * BM;                  // differentiation to LLPC

// Short-lived plasma cells (CD20-, relatively RTX-resistant)
dxdt_SLPC = kSLPC_form * BGC * (1.0 - EMPA)
           - kSLPC_die  * SLPC * (1.0 + Epred_B);

// Long-lived plasma cells (bone marrow niche, CD20-, RTX-resistant)
dxdt_LLPC = kLLPC_form * BM
           - kLLPC_die  * LLPC;              // RTX-resistant (no CD20)

// Tfh cells
dxdt_Tfh  = kTfh_prod * (1.0 - Epred_T)
           - kTfh_die  * Tfh;

// Regulatory T cells
dxdt_Treg = kTreg_prod + 0.3 * Epred_T      // corticosteroids partially expand Treg
           - kTreg_die  * Treg;

// ═══════════════════════════════════════════════════════════════════════════
// ANTIBODY BIOMARKERS
// ═══════════════════════════════════════════════════════════════════════════

// Anti-Dsg3 IgG (primary pathogenic antibody)
dxdt_Ab3 = kAb3_prod * (SLPC + LLPC * 2.0)              // LLPC more productive
           - kAb3_die  * Ab3 * (1.0 + E_Ab_cat);        // FcRn blockade accelerates

// Anti-Dsg1 IgG (cutaneous pemphigus component)
dxdt_Ab1 = kAb1_prod * (SLPC + LLPC)
           - kAb1_die  * Ab1 * (1.0 + E_Ab_cat);

// ═══════════════════════════════════════════════════════════════════════════
// DISEASE STATE
// ═══════════════════════════════════════════════════════════════════════════

// Dsg3 protein loss (driven by Ab3-mediated internalization)
dxdt_Dsg3_loss = kDsg3_int  * Ab3 * Dsg3_avail
                - kDsg3_rest * Dsg3_loss;

// PDAI disease activity (driven by fraction of Dsg3 lost)
dxdt_PDAI_state = kPDAI_up * (Dsg3_loss / Dsg3_base)    // blister formation
                 - kPDAI_dn * PDAI_state;

// Complement activation (IgG1 subclass Ab3 drives complement)
double Ab3_IgG1_frac = 0.25;  // IgG1 ~25% of total Ab3
dxdt_Comp_act = kComp_on  * Ab3 * Ab3_IgG1_frac
               - kComp_off * Comp_act;

// Bone mineral density loss (cumulative corticosteroid effect)
dxdt_Cort_bone = kBone_loss * Cp1_pred
                - kBone_rest * Cort_bone;

$TABLE
// Derived outputs for clinical readout
double Anti_Dsg3   = Ab3;
double Anti_Dsg1   = Ab1;
double PDAI        = PDAI_state;
double B_total_obs = BN + BGC + BM;
double PC_total    = SLPC + LLPC;
double Pred_Cp     = Cp1_pred;
double RTX_Cp      = CR1;
double MPA_Cp      = CMPA;

// Blister area approximation (0-100% BSA)
double BSA_blisters = fmin(100.0, 0.3 * PDAI_state);

// Clinical response categories
double CR_off = (PDAI_state < 2.0 && Cp1_pred < 0.1) ? 1.0 : 0.0;
double PR_min = (PDAI_state < 8.0 && Cp1_pred < 0.2) ? 1.0 : 0.0;

// Bone density Z-score estimated loss
double BMD_loss_pct = fmin(15.0, Cort_bone * 0.5);

$INIT
DEPOT_pred = 0
Cp1_pred   = 0
Cp2_pred   = 0
CR1        = 0
CR2        = 0
DEPOT_MPA  = 0
CMPA       = 0
DEPOT_IVIg = 0
CIVG       = 0
CEFG       = 0
BN         = 50.0
BGC        = 5.0
BM         = 20.0
SLPC       = 3.0
LLPC       = 2.0
Tfh        = 10.0
Treg       = 8.0
Ab3        = 180.0
Ab1        = 80.0
Dsg3_loss  = 30.0
PDAI_state = 20.0
Comp_act   = 5.0
Cort_bone  = 0.0
'

## ─────────────────────────────────────────────────────────────────────────────
## 2. COMPILE MODEL
## ─────────────────────────────────────────────────────────────────────────────

mod <- mcode("PV_QSP", pv_model_code)

## ─────────────────────────────────────────────────────────────────────────────
## 3. TREATMENT SCENARIOS
## ─────────────────────────────────────────────────────────────────────────────

# Time grid: 2 years (8760 hours), observations every 24h
sim_time <- seq(0, 8760, by = 24)

# ── Scenario 1: High-dose corticosteroids alone (historical standard) ──────
# Prednisone 1.5 mg/kg/day (assume 70 kg → 105 mg/day oral)
ev_cs_high <- ev(time = 0, amt = 105, cmt = "DEPOT_pred", ii = 24,
                 addl = 364, rate = 0)   # daily for 1 year

# Taper schedule: halve every 8 weeks after week 8
ev_cs_taper1 <- ev(time = 8*7*24,  amt = 52.5, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_taper2 <- ev(time = 16*7*24, amt = 25.0, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_taper3 <- ev(time = 24*7*24, amt = 12.5, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_taper4 <- ev(time = 32*7*24, amt = 5.0,  cmt = "DEPOT_pred", ii = 24, addl = 363 - 32*7)

ev_sc1 <- c(ev_cs_high, ev_cs_taper1, ev_cs_taper2, ev_cs_taper3, ev_cs_taper4)

# ── Scenario 2: Rituximab + Low-dose prednisone (RITUX3 protocol) ──────────
# RTX 1000mg IV at weeks 0 and 2, + prednisone 0.5 mg/kg → 35 mg/day
ev_rtx1 <- ev(time = 0,       amt = 1000, cmt = "CR1",       rate = 1000/6)
ev_rtx2 <- ev(time = 2*7*24,  amt = 1000, cmt = "CR1",       rate = 1000/6)
ev_rtx_re1 <- ev(time = 26*7*24, amt = 500, cmt = "CR1",     rate = 500/6) # re-treat at 6 mo

ev_cs_low <- ev(time = 0, amt = 35, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_t1  <- ev(time = 8*7*24, amt = 17.5, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_t2  <- ev(time = 16*7*24, amt = 5.0, cmt = "DEPOT_pred", ii = 24, addl = 363 - 16*7)

ev_sc2 <- c(ev_rtx1, ev_rtx2, ev_rtx_re1, ev_cs_low, ev_cs_t1, ev_cs_t2)

# ── Scenario 3: MMF + Moderate prednisone (Murrell protocol) ───────────────
# MMF 2g/day oral + prednisone 1 mg/kg (70 mg/day) tapering
ev_mmf <- ev(time = 0, amt = 1000, cmt = "DEPOT_MPA", ii = 12, addl = 728,
             rate = 0)   # 1g BID
ev_cs_mod <- ev(time = 0, amt = 70, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_mt1 <- ev(time = 8*7*24, amt = 35, cmt = "DEPOT_pred", ii = 24, addl = 55)
ev_cs_mt2 <- ev(time = 16*7*24, amt = 15, cmt = "DEPOT_pred", ii = 24, addl = 363 - 16*7)

ev_sc3 <- c(ev_mmf, ev_cs_mod, ev_cs_mt1, ev_cs_mt2)

# ── Scenario 4: Rituximab + MMF combination ────────────────────────────────
ev_sc4 <- c(ev_rtx1, ev_rtx2, ev_rtx_re1, ev_mmf, ev_cs_low, ev_cs_t1, ev_cs_t2)

# ── Scenario 5: IV Pulse steroid + rituximab (severe/refractory PV) ────────
# IV methylprednisolone 1g × 3 days at induction, then rituximab
ev_ivmp1 <- ev(time = 0, amt = 1000 * 1.25, cmt = "DEPOT_pred",   # 1g IVMP in pred equivalents
               rate = 1000/6, addl = 2, ii = 24)
ev_sc5 <- c(ev_ivmp1, ev_rtx1, ev_rtx2, ev_cs_low, ev_cs_t1, ev_cs_t2)

# ── Scenario 6: Efgartigimod + low-dose pred (PEMPHIX trial concept) ───────
# Efgartigimod 10 mg/kg IV q4w; FcRn blocker → rapid IgG reduction
ev_efg_wk0  <- ev(time = 0,       amt = 700, cmt = "CEFG", rate = 700/2)  # 70 kg × 10 mg/kg
ev_efg_wk4  <- ev(time = 4*7*24,  amt = 700, cmt = "CEFG", rate = 700/2)
ev_efg_wk8  <- ev(time = 8*7*24,  amt = 700, cmt = "CEFG", rate = 700/2)
ev_efg_wk12 <- ev(time = 12*7*24, amt = 700, cmt = "CEFG", rate = 700/2)
ev_efg_wk16 <- ev(time = 16*7*24, amt = 700, cmt = "CEFG", rate = 700/2)
ev_efg_wk20 <- ev(time = 20*7*24, amt = 700, cmt = "CEFG", rate = 700/2)

ev_sc6 <- c(ev_efg_wk0, ev_efg_wk4, ev_efg_wk8, ev_efg_wk12,
            ev_efg_wk16, ev_efg_wk20, ev_cs_low, ev_cs_t1, ev_cs_t2)

## ─────────────────────────────────────────────────────────────────────────────
## 4. SIMULATION
## ─────────────────────────────────────────────────────────────────────────────

run_scenario <- function(ev_obj, label) {
  out <- mod %>%
    ev(ev_obj) %>%
    mrgsim(end = 8760, delta = 24) %>%
    as.data.frame()
  out$scenario <- label
  out
}

cat("Running 6 treatment scenarios...\n")
sc1 <- run_scenario(ev_sc1, "High-dose CS (historical)")
sc2 <- run_scenario(ev_sc2, "RTX + low CS (RITUX3)")
sc3 <- run_scenario(ev_sc3, "MMF + moderate CS")
sc4 <- run_scenario(ev_sc4, "RTX + MMF combination")
sc5 <- run_scenario(ev_sc5, "IV pulse + RTX (severe)")
sc6 <- run_scenario(ev_sc6, "Efgartigimod + low CS")

results <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6)
results$week <- results$time / 168   # hours to weeks

## ─────────────────────────────────────────────────────────────────────────────
## 5. PLOTS
## ─────────────────────────────────────────────────────────────────────────────

theme_pv <- theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E8F4FF"),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 10))

palette_sc <- c(
  "High-dose CS (historical)"  = "#E74C3C",
  "RTX + low CS (RITUX3)"      = "#2980B9",
  "MMF + moderate CS"          = "#27AE60",
  "RTX + MMF combination"      = "#8E44AD",
  "IV pulse + RTX (severe)"    = "#F39C12",
  "Efgartigimod + low CS"      = "#16A085"
)

## ── Figure 1: Disease Activity (PDAI) ──────────────────────────────────────
p1 <- ggplot(results, aes(week, PDAI, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 2,  linetype = "dashed", color = "black", size = 0.8) +
  geom_hline(yintercept = 8,  linetype = "dotted", color = "grey50", size = 0.8) +
  annotate("text", x = 100, y = 2.5, label = "CR-off threshold", size = 3.2) +
  annotate("text", x = 100, y = 8.5, label = "PR-min threshold", size = 3.2) +
  scale_color_manual(values = palette_sc) +
  scale_x_continuous(breaks = seq(0, 52, 13), labels = paste0("W", seq(0, 52, 13))) +
  coord_cartesian(xlim = c(0, 52)) +
  labs(title = "Pemphigus Vulgaris — Disease Activity (PDAI)",
       subtitle = "First year; dashed = CR-off target (PDAI < 2)",
       x = "Week", y = "PDAI Score", color = "Treatment") +
  theme_pv
print(p1)

## ── Figure 2: Anti-Dsg3 IgG dynamics ──────────────────────────────────────
p2 <- ggplot(results, aes(week, Anti_Dsg3, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#CC0000") +
  annotate("text", x = 80, y = 22, label = "Pathogenic threshold (20 U/mL)", size = 3.2) +
  scale_color_manual(values = palette_sc) +
  coord_cartesian(xlim = c(0, 104)) +
  labs(title = "Anti-Dsg3 IgG Kinetics",
       subtitle = "2-year follow-up; dashed = pathogenic cutoff",
       x = "Week", y = "Anti-Dsg3 IgG (U/mL)", color = "Treatment") +
  theme_pv
print(p2)

## ── Figure 3: B cell dynamics under RTX ────────────────────────────────────
p3 <- results %>%
  filter(scenario %in% c("RTX + low CS (RITUX3)", "RTX + MMF combination",
                         "High-dose CS (historical)")) %>%
  pivot_longer(cols = c("BN", "BGC", "BM"), names_to = "Bcell_type",
               values_to = "count") %>%
  ggplot(aes(week, count, color = scenario, linetype = Bcell_type)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = palette_sc) +
  scale_linetype_manual(values = c("BN" = "solid", "BGC" = "dashed", "BM" = "dotted")) +
  coord_cartesian(xlim = c(0, 104)) +
  labs(title = "B Cell Dynamics (Naive, GC, Memory)",
       x = "Week", y = "B Cell Count (×10⁶/L)",
       color = "Treatment", linetype = "B Cell Type") +
  theme_pv
print(p3)

## ── Figure 4: Prednisolone PK ──────────────────────────────────────────────
p4 <- results %>%
  filter(scenario %in% c("High-dose CS (historical)", "RTX + low CS (RITUX3)",
                         "MMF + moderate CS")) %>%
  ggplot(aes(week, Pred_Cp, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = palette_sc) +
  coord_cartesian(xlim = c(0, 52)) +
  labs(title = "Prednisolone Plasma Concentration",
       x = "Week", y = "Prednisolone Cp (mg/L)", color = "Treatment") +
  theme_pv
print(p4)

## ── Figure 5: Bone loss risk (cumulative steroid effect) ───────────────────
p5 <- ggplot(results, aes(week, BMD_loss_pct, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#8B4513") +
  annotate("text", x = 60, y = 10.5, label = "Osteoporosis risk threshold (~10%)", size = 3.2) +
  scale_color_manual(values = palette_sc) +
  labs(title = "Cumulative Bone Mineral Density Loss Risk",
       subtitle = "Steroid-induced osteoporosis surrogate endpoint",
       x = "Week", y = "Estimated BMD Loss (%)", color = "Treatment") +
  theme_pv
print(p5)

## ── Figure 6: Complement activation ───────────────────────────────────────
p6 <- ggplot(results, aes(week, Comp_act, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = palette_sc) +
  coord_cartesian(xlim = c(0, 52)) +
  labs(title = "Complement Activation Dynamics",
       subtitle = "Driven by IgG1 anti-Dsg3 subclass",
       x = "Week", y = "Complement Activity (arb units)", color = "Treatment") +
  theme_pv
print(p6)

## ─────────────────────────────────────────────────────────────────────────────
## 6. SUMMARY TABLE
## ─────────────────────────────────────────────────────────────────────────────

summary_tbl <- results %>%
  group_by(scenario) %>%
  summarize(
    PDAI_wk4     = round(PDAI[week >= 4  & week < 5][1], 1),
    PDAI_wk12    = round(PDAI[week >= 12 & week < 13][1], 1),
    PDAI_wk24    = round(PDAI[week >= 24 & week < 25][1], 1),
    PDAI_wk52    = round(PDAI[week >= 52 & week < 53][1], 1),
    Ab3_wk12     = round(Anti_Dsg3[week >= 12 & week < 13][1], 0),
    Ab3_wk52     = round(Anti_Dsg3[week >= 52 & week < 53][1], 0),
    CR_off_wk52  = CR_off[week >= 52 & week < 53][1],
    BMD_loss_wk52= round(BMD_loss_pct[week >= 52 & week < 53][1], 1),
    .groups = "drop"
  )

cat("\n=== PV QSP Model — Treatment Comparison Summary ===\n")
print(as.data.frame(summary_tbl))

## ─────────────────────────────────────────────────────────────────────────────
## 7. DOSE-RESPONSE ANALYSIS (RTX dose)
## ─────────────────────────────────────────────────────────────────────────────

cat("\nRunning RTX dose-response analysis...\n")

rtx_doses <- c(100, 250, 500, 750, 1000, 1500)
dr_results <- lapply(rtx_doses, function(dose) {
  ev_rtx_d1 <- ev(time = 0,      amt = dose, cmt = "CR1", rate = dose/6)
  ev_rtx_d2 <- ev(time = 2*7*24, amt = dose, cmt = "CR1", rate = dose/6)
  ev_d <- c(ev_rtx_d1, ev_rtx_d2, ev_cs_low, ev_cs_t1, ev_cs_t2)

  out <- mod %>% ev(ev_d) %>% mrgsim(end = 52*7*24, delta = 24) %>% as.data.frame()
  out$RTX_dose_mg <- dose
  out$week <- out$time / 168
  out
})
dr_df <- bind_rows(dr_results)

p_dr <- dr_df %>%
  ggplot(aes(week, PDAI, group = RTX_dose_mg, color = factor(RTX_dose_mg))) +
  geom_line(size = 1.0) +
  scale_color_brewer(palette = "YlOrRd", direction = -1) +
  labs(title = "RTX Dose-Response (PDAI over 1 year)",
       subtitle = "Each pair dose × 2 infusions (week 0, 2) + low-dose prednisone",
       x = "Week", y = "PDAI Score", color = "RTX Dose (mg)") +
  theme_pv
print(p_dr)

cat("\nPemphigus Vulgaris QSP Model simulation complete.\n")
cat("Key outputs: PDAI trajectory, anti-Dsg3 IgG kinetics,\n")
cat("             B cell dynamics, bone loss risk, dose-response\n")
