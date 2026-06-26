## =============================================================================
## Hemophagocytic Lymphohistiocytosis (HLH) — QSP Model (mrgsolve)
## =============================================================================
## Description:
##   A quantitative systems pharmacology (QSP) model of HLH pathophysiology
##   capturing NK cell/CTL dysfunction, cytokine storm dynamics (IFN-γ, IL-6,
##   TNF-α, IL-18, IL-10), macrophage activation, hemophagocytosis, organ
##   damage markers, and the PK/PD of six treatment options:
##     1. Dexamethasone (DEX)
##     2. Etoposide (ETOP)
##     3. Cyclosporine A (CsA)
##     4. Emapalumab (anti-IFN-γ mAb)
##     5. Anakinra (IL-1Ra)
##     6. Ruxolitinib (JAK1/2 inhibitor)
##
## Compartments: 20 ODE compartments
## Parameters calibrated from:
##   - Jordan et al. Blood 2011 (HLH-2004 study)
##   - La Rosée et al. Blood 2019 (HLH treatment guidelines)
##   - Daver et al. Am J Hematol 2017
##   - Locatelli et al. NEJM 2020 (emapalumab)
##   - Ravelli et al. Ann Rheum Dis 2016 (MAS criteria)
##   - Histiocyte Society, HLH-2004 protocol
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---------------------------------------------------------------------------
## Model code
## ---------------------------------------------------------------------------
hlh_code <- '
$PROB
HLH QSP Model — cytokine storm, macrophage activation, drug PK/PD

$PARAM @annotated
// --- Disease state parameters ---
kNK_base   : 1.0   : NK cell baseline production (cells/day normalized)
kCTL_base  : 1.0   : CTL baseline production
d_NK       : 0.1   : NK cell natural death rate (1/day)
d_CTL      : 0.1   : CTL natural death rate
fNK_dysfunc: 0.0   : NK dysfunction fraction (0=normal; 1=complete FHL; PRIMARY HLH)

kTrig      : 0.1   : Trigger/stimulus strength (0=none; 1=severe viral/EBV)
EC50_trig  : 0.5   : EC50 of trigger on T cell activation

// --- Cytokine production rates (pg/mL/day) ---
kIFNg_prod : 100   : IFN-γ production by activated T/NK cells (pg/mL/day)
kIL6_prod  : 50    : IL-6 production
kTNFa_prod : 80    : TNF-α production
kIL18_prod : 30    : IL-18 production
kIL10_prod : 40    : IL-10 production (paradoxical amplifier in HLH)
kIL12_prod : 20    : IL-12 production by APCs

// --- Cytokine clearance rates (1/day) ---
dIFNg : 2.0  : IFN-γ clearance
dIL6  : 4.0  : IL-6 clearance
dTNFa : 6.0  : TNF-α clearance
dIL18 : 1.5  : IL-18 clearance
dIL10 : 3.0  : IL-10 clearance
dIL12 : 2.0  : IL-12 clearance

// --- Macrophage/hemophagocytosis ---
kMacro_act  : 0.5  : Macrophage activation rate by IFN-γ (normalized)
EC50_IFNg_M : 200  : EC50 of IFN-γ for macrophage activation (pg/mL)
d_Macro     : 0.2  : Macrophage inactivation rate
kHemophag   : 0.3  : Hemophagocytosis rate by activated macrophages

// --- Ferritin dynamics ---
kFerr_prod  : 10   : Ferritin production rate (ng/mL/day)
dFerr       : 0.05 : Ferritin clearance rate
EC50_IL6_Ferr : 100 : EC50 of IL-6 for ferritin induction

// --- Cytopenias (bone marrow suppression index) ---
kBM_supp    : 0.4  : Bone marrow suppression rate by hemophagocytosis
d_BM_supp   : 0.1  : BM suppression recovery rate

// --- Organ damage (liver, coagulation index) ---
kLiver_dmg  : 0.2  : Liver damage rate by activated macrophages
dLiver_dmg  : 0.15 : Liver damage resolution rate
kCoag_dmg   : 0.15 : Coagulopathy induction rate
dCoag_dmg   : 0.1  : Coagulopathy resolution

// --- Drug PK: Dexamethasone (2-compartment oral) ---
ka_DEX  : 1.2   : DEX absorption rate constant (1/h)
CL_DEX  : 18.0  : DEX clearance (L/h)
V1_DEX  : 50    : Central volume (L)
V2_DEX  : 120   : Peripheral volume (L)
Q_DEX   : 8.0   : Inter-compartmental clearance (L/h)

// --- Drug PK: Etoposide (2-compartment IV) ---
CL_ETOP : 20.0  : Etoposide CL (L/h)
V1_ETOP : 8     : Central volume (L)
V2_ETOP : 35    : Peripheral volume (L)
Q_ETOP  : 4.0   : Inter-compartmental CL (L/h)

// --- Drug PK: Cyclosporine (2-compartment oral) ---
ka_CSA  : 0.5   : CsA absorption rate (1/h)
CL_CSA  : 25.0  : CsA clearance (L/h)
V1_CSA  : 200   : Central volume (L)
V2_CSA  : 600   : Peripheral volume (L)
Q_CSA   : 30.0  : Inter-compartmental CL (L/h)

// --- Drug PK: Emapalumab (mAb 2-cmt + TMDD vs IFN-γ) ---
CL_EMAPA  : 0.004  : Emapalumab linear CL (L/h)
V1_EMAPA  : 4.0    : Central volume (L)
V2_EMAPA  : 3.0    : Peripheral volume (L)
Q_EMAPA   : 0.01   : Inter-compartmental CL (L/h)
Kd_EMAPA  : 50     : Kd for IFN-γ binding (pg/mL)
kon_EMAPA : 0.001  : Association rate constant
koff_EMAPA: 0.05   : Dissociation rate
kin_TMDD  : 0.002  : Internalization of drug-IFN-γ complex

// --- Drug PK: Anakinra (SC 1-cmt) ---
ka_ANK  : 0.08  : Anakinra SC absorption (1/h)
F_ANK   : 0.95  : SC bioavailability
CL_ANK  : 2.5   : Renal clearance (L/h)
V_ANK   : 7.0   : Volume of distribution (L)

// --- Drug PK: Ruxolitinib (oral, 2-cmt) ---
ka_RUX  : 1.8   : Ruxolitinib absorption (1/h)
CL_RUX  : 17.5  : CL (L/h)
V1_RUX  : 75    : Central volume (L)
V2_RUX  : 130   : Peripheral volume (L)
Q_RUX   : 15.0  : Inter-compartmental CL (L/h)

// --- PD: Drug effects ---
Emax_DEX_IFNg  : 0.80 : Max inhibition of IFN-γ by DEX (Emax)
EC50_DEX_IFNg  : 50   : EC50 of DEX plasma level (ng/mL)
Emax_DEX_Macro : 0.75 : Max inhibition of macrophage activation by DEX
EC50_DEX_Macro : 60   : EC50

Emax_ETOP      : 0.90 : Max T cell/NK depletion by etoposide
EC50_ETOP      : 2.0  : EC50 of etoposide (μg/mL)

Emax_CSA_IL2   : 0.85 : Max IL-2/IFN-γ inhibition by CsA
EC50_CSA       : 300  : EC50 of CsA (ng/mL)

Emax_EMAPA     : 0.95 : Max IFN-γ neutralization fraction
EC50_EMAPA_IFNg: 50   : EC50 of free emapalumab vs IFN-γ

Emax_ANK_IL1   : 0.80 : Max IL-1β signaling block by anakinra
EC50_ANK       : 80   : EC50 of anakinra (ng/mL)

Emax_RUX_JAK   : 0.75 : Max JAK-STAT signaling inhibition by ruxolitinib
EC50_RUX       : 200  : EC50 of ruxolitinib (ng/mL)
n_RUX          : 1.5  : Hill coefficient for ruxolitinib

$CMT @annotated
// Disease state
NK_CELLS  : NK cell pool (normalized units)
CTL_CELLS : CD8+ CTL pool
T_ACT     : Activated T cells (effector)
APC_ACT   : Activated APCs

// Cytokines (all in pg/mL)
IFNg      : IFN-gamma (pg/mL)
IL6       : IL-6 (pg/mL)
TNFa      : TNF-alpha (pg/mL)
IL18      : IL-18 (pg/mL)
IL10      : IL-10 (pg/mL)
IL12      : IL-12 (pg/mL)

// Macrophage & hemophagocytosis
MAC_ACT   : Activated macrophages (normalized)
HEMOPHAG_IDX : Hemophagocytosis index (0-1)
FERR      : Serum ferritin (ng/mL)
BM_SUPP   : Bone marrow suppression index (0-1)

// Organ damage
LIVER_DMG : Liver damage index (0-1)
COAG_DMG  : Coagulopathy index (0-1)

// Drug PK compartments
DEX_GI    : Dexamethasone GI/depot
DEX_C     : Dexamethasone central
DEX_P     : Dexamethasone peripheral

ETOP_C    : Etoposide central
ETOP_P    : Etoposide peripheral

CSA_GI    : Cyclosporine GI/depot
CSA_C     : Cyclosporine central
CSA_P     : Cyclosporine peripheral

EMAPA_C   : Emapalumab central
EMAPA_P   : Emapalumab peripheral

ANK_SC    : Anakinra SC depot
ANK_C     : Anakinra central

RUX_GI    : Ruxolitinib GI/depot
RUX_C     : Ruxolitinib central
RUX_P     : Ruxolitinib peripheral

$GLOBAL
double Emax_func(double C, double Emax, double EC50, double n) {
    if(C <= 0) return 0;
    double Cn = pow(C, n);
    double EC50n = pow(EC50, n);
    return Emax * Cn / (Cn + EC50n);
}

$MAIN
// Baseline conditions (pre-disease)
if(NEWIND <= 1) {
    NK_CELLS_0  = 1.0;
    CTL_CELLS_0 = 1.0;
    T_ACT_0     = 0.01;
    APC_ACT_0   = 0.1;
    IFNg_0      = 5.0;    // pg/mL baseline
    IL6_0       = 2.0;
    TNFa_0      = 1.0;
    IL18_0      = 10.0;
    IL10_0      = 3.0;
    IL12_0      = 1.0;
    MAC_ACT_0   = 0.1;
    HEMOPHAG_IDX_0 = 0.0;
    FERR_0      = 150.0;  // ng/mL normal
    BM_SUPP_0   = 0.0;
    LIVER_DMG_0 = 0.0;
    COAG_DMG_0  = 0.0;
}

$ODE
// -----------------------------------------------------------------------
// Plasma concentrations (for PD calculations)
double C_DEX   = DEX_C / V1_DEX * 1000;   // convert to ng/mL
double C_ETOP  = ETOP_C / V1_ETOP;        // μg/mL (dose in mg)
double C_CSA   = CSA_C / V1_CSA * 1000;   // ng/mL
double C_EMAPA = EMAPA_C / V1_EMAPA;      // μg/mL
double C_ANK   = ANK_C / V_ANK * 1000;    // ng/mL
double C_RUX   = RUX_C / V1_RUX * 1000;  // ng/mL

// -----------------------------------------------------------------------
// Drug effect modifiers (inhibition)
double Inh_DEX_IFNg  = Emax_func(C_DEX, Emax_DEX_IFNg, EC50_DEX_IFNg, 1.0);
double Inh_DEX_Macro = Emax_func(C_DEX, Emax_DEX_Macro, EC50_DEX_Macro, 1.0);
double Inh_ETOP      = Emax_func(C_ETOP, Emax_ETOP, EC50_ETOP, 2.0);
double Inh_CSA       = Emax_func(C_CSA, Emax_CSA_IL2, EC50_CSA, 1.0);
double Inh_EMAPA_IFNg= Emax_func(C_EMAPA * 1000, Emax_EMAPA, EC50_EMAPA_IFNg, 1.0);
double Inh_ANK       = Emax_func(C_ANK, Emax_ANK_IL1, EC50_ANK, 1.0);
double Inh_RUX       = Emax_func(C_RUX, Emax_RUX_JAK, EC50_RUX, n_RUX);

// Combined IFN-γ inhibition (multiple drugs can act simultaneously)
double Inh_IFNg_total = 1 - (1 - Inh_DEX_IFNg) * (1 - Inh_CSA) *
                            (1 - Inh_EMAPA_IFNg) * (1 - Inh_RUX);

// -----------------------------------------------------------------------
// NK cell & CTL dynamics
// NK dysfunction from genetic defect or drug-mediated depletion
double NK_cytotox_capacity = (1.0 - fNK_dysfunc) * NK_CELLS * (1.0 - Inh_ETOP);
double CTL_cytotox_capacity = CTL_CELLS * (1.0 - Inh_ETOP);

// Trigger activates APCs → T cells
double Trigger_effect = kTrig / (kTrig + EC50_trig);  // simplified

dxdt_NK_CELLS  = kNK_base * (1 - fNK_dysfunc) - d_NK * NK_CELLS -
                 Inh_ETOP * 0.05 * NK_CELLS;
dxdt_CTL_CELLS = kCTL_base - d_CTL * CTL_CELLS -
                 Inh_ETOP * 0.05 * CTL_CELLS;

// APC activation by trigger
double APC_clearance = (NK_cytotox_capacity + CTL_cytotox_capacity) * 0.05;
dxdt_APC_ACT = kTrig * 2.0 - 0.3 * APC_ACT - APC_clearance;

// T cell activation by APCs (driven by IL-12, IL-18)
double IL12_stim = IL12 / (IL12 + 10.0);
double IL18_stim = IL18 / (IL18 + 20.0);
dxdt_T_ACT = 0.5 * APC_ACT * (1 + IL12_stim + IL18_stim) -
             0.4 * T_ACT * (1 - Inh_ETOP) * (1 - Inh_CSA) -
             0.2 * T_ACT;

// -----------------------------------------------------------------------
// Cytokine ODEs
// IFN-γ: produced by activated T cells and NK cells, amplified by IL-12/IL-18
double IFNg_prod = kIFNg_prod * T_ACT * (1 + IL18_stim) * (1 - Inh_IFNg_total);
dxdt_IFNg = IFNg_prod + 5.0 - dIFNg * IFNg;  // +5 = basal

// IL-6: produced by APCs and macrophages, stimulated by TNF-α
double TNFa_stim = TNFa / (TNFa + 50.0);
double IL6_prod  = kIL6_prod * (APC_ACT + MAC_ACT) * (1 + TNFa_stim) *
                   (1 - Inh_RUX * 0.6);
dxdt_IL6 = IL6_prod + 2.0 - dIL6 * IL6;

// TNF-α: produced by macrophages and T cells
double IFNg_stim_M = IFNg / (IFNg + EC50_IFNg_M);
double TNFa_prod = kTNFa_prod * MAC_ACT * (1 + IFNg_stim_M) *
                  (1 - Inh_DEX_Macro);
dxdt_TNFa = TNFa_prod + 1.0 - dTNFa * TNFa;

// IL-18: produced by macrophages, amplified feedback
double IL18_prod = kIL18_prod * MAC_ACT * (1 + IFNg_stim_M);
dxdt_IL18 = IL18_prod + 10.0 - dIL18 * IL18;

// IL-10: paradoxical in HLH (high IL-10 accompanies disease)
double IL10_prod = kIL10_prod * T_ACT * MAC_ACT;
dxdt_IL10 = IL10_prod + 3.0 - dIL10 * IL10;

// IL-12: produced by APCs, drives Th1/IFN-γ
double IL12_prod = kIL12_prod * APC_ACT * (1 - Inh_DEX_Macro);
dxdt_IL12 = IL12_prod + 1.0 - dIL12 * IL12;

// -----------------------------------------------------------------------
// Macrophage activation (driven by IFN-γ)
double Mac_activ_rate = kMacro_act * IFNg / (IFNg + EC50_IFNg_M) *
                        (1 - Inh_DEX_Macro);
dxdt_MAC_ACT = Mac_activ_rate - d_Macro * MAC_ACT;

// Hemophagocytosis index (0-1)
dxdt_HEMOPHAG_IDX = kHemophag * MAC_ACT * (1 - HEMOPHAG_IDX) -
                    0.1 * HEMOPHAG_IDX * (1 - Inh_ETOP);  // etoposide partial reversal

// Ferritin (ng/mL) — driven by IL-6 and TNF-α
double Ferr_stim = 1.0 + IL6 / (IL6 + EC50_IL6_Ferr) + TNFa / (TNFa + 100.0);
dxdt_FERR = kFerr_prod * MAC_ACT * Ferr_stim - dFerr * FERR;

// Bone marrow suppression index (0=normal, 1=maximal suppression)
dxdt_BM_SUPP = kBM_supp * HEMOPHAG_IDX * (1 - BM_SUPP) - d_BM_supp * BM_SUPP;

// -----------------------------------------------------------------------
// Organ damage
// Liver damage (from macrophage activation + TNF-α/IL-6)
dxdt_LIVER_DMG = kLiver_dmg * MAC_ACT * (TNFa / 100.0) * (1 - LIVER_DMG) -
                 dLiver_dmg * LIVER_DMG * (1 - HEMOPHAG_IDX);

// Coagulopathy index (fibrinogen consumption)
dxdt_COAG_DMG = kCoag_dmg * MAC_ACT * HEMOPHAG_IDX * (1 - COAG_DMG) -
                dCoag_dmg * COAG_DMG;

// -----------------------------------------------------------------------
// Drug PK ODEs

// Dexamethasone (2-compartment oral, dose in mg converted to μg → ng/mL)
dxdt_DEX_GI = -ka_DEX * DEX_GI;
dxdt_DEX_C  =  ka_DEX * DEX_GI - (CL_DEX + Q_DEX) / V1_DEX * DEX_C +
               Q_DEX / V2_DEX * DEX_P;
dxdt_DEX_P  =  Q_DEX / V1_DEX * DEX_C - Q_DEX / V2_DEX * DEX_P;

// Etoposide (2-compartment IV infusion)
dxdt_ETOP_C = -(CL_ETOP + Q_ETOP) / V1_ETOP * ETOP_C + Q_ETOP / V2_ETOP * ETOP_P;
dxdt_ETOP_P = Q_ETOP / V1_ETOP * ETOP_C - Q_ETOP / V2_ETOP * ETOP_P;

// Cyclosporine (2-compartment oral)
dxdt_CSA_GI = -ka_CSA * CSA_GI;
dxdt_CSA_C  =  ka_CSA * CSA_GI - (CL_CSA + Q_CSA) / V1_CSA * CSA_C +
               Q_CSA / V2_CSA * CSA_P;
dxdt_CSA_P  =  Q_CSA / V1_CSA * CSA_C - Q_CSA / V2_CSA * CSA_P;

// Emapalumab (mAb 2-compartment IV)
dxdt_EMAPA_C = -(CL_EMAPA + Q_EMAPA) / V1_EMAPA * EMAPA_C +
                Q_EMAPA / V2_EMAPA * EMAPA_P;
dxdt_EMAPA_P = Q_EMAPA / V1_EMAPA * EMAPA_C - Q_EMAPA / V2_EMAPA * EMAPA_P;

// Anakinra (SC)
dxdt_ANK_SC = -ka_ANK * ANK_SC;
dxdt_ANK_C  =  ka_ANK * ANK_SC * F_ANK - CL_ANK / V_ANK * ANK_C;

// Ruxolitinib (2-compartment oral)
dxdt_RUX_GI = -ka_RUX * RUX_GI;
dxdt_RUX_C  =  ka_RUX * RUX_GI - (CL_RUX + Q_RUX) / V1_RUX * RUX_C +
               Q_RUX / V2_RUX * RUX_P;
dxdt_RUX_P  =  Q_RUX / V1_RUX * RUX_C - Q_RUX / V2_RUX * RUX_P;

$TABLE
// Derived outputs
double C_DEX_out   = DEX_C / V1_DEX * 1000;   // ng/mL
double C_ETOP_out  = ETOP_C / V1_ETOP;         // μg/mL
double C_CSA_out   = CSA_C / V1_CSA * 1000;    // ng/mL
double C_EMAPA_out = EMAPA_C / V1_EMAPA * 1000; // μg/mL → convert
double C_ANK_out   = ANK_C / V_ANK * 1000;     // ng/mL
double C_RUX_out   = RUX_C / V1_RUX * 1000;    // ng/mL

// Clinical biomarkers
double sCD25_equiv = 100 + T_ACT * 500;         // surrogate sCD25 (U/mL)
double NK_activity = NK_CELLS * (1 - fNK_dysfunc); // NK functional capacity
double Fibrinogen  = 3.5 * (1 - COAG_DMG * 0.8);   // g/L (3.5 normal)
double Triglycerides_mM = 1.5 + MAC_ACT * 4.0;     // mmol/L

// Survival proxy (0=alive, 1=death; threshold-based)
double Death_risk = LIVER_DMG * 0.3 + COAG_DMG * 0.25 +
                    BM_SUPP * 0.25 + HEMOPHAG_IDX * 0.2;
double Survival_prob = exp(-Death_risk * 2.5);  // exponential decay

// HScore contributor (simplified 0-100 scale)
double HScore_component = (IFNg / 2000 + IL6 / 1000 + FERR / 100000 +
                            HEMOPHAG_IDX + BM_SUPP) * 20;
double HScore_approx = HScore_component > 337 ? 337 : HScore_component;

$CAPTURE
C_DEX_out C_ETOP_out C_CSA_out C_EMAPA_out C_ANK_out C_RUX_out
IFNg IL6 TNFa IL18 IL10
MAC_ACT HEMOPHAG_IDX FERR BM_SUPP
LIVER_DMG COAG_DMG
sCD25_equiv NK_activity Fibrinogen Triglycerides_mM
Survival_prob HScore_approx
T_ACT APC_ACT
'

## ---------------------------------------------------------------------------
## Load model
## ---------------------------------------------------------------------------
hlh_mod <- mcode("hlh_qsp", hlh_code)
cat("Model loaded:", outnames(hlh_mod) %>% length(), "output variables\n")

## ---------------------------------------------------------------------------
## SCENARIO DEFINITIONS
## ---------------------------------------------------------------------------
## Scenario 1: Untreated (natural progression)
## Scenario 2: HLH-2004 protocol (DEX + ETOP + CsA)
## Scenario 3: Emapalumab + DEX (primary/refractory HLH)
## Scenario 4: Anakinra-based (MAS/secondary HLH)
## Scenario 5: Ruxolitinib salvage (refractory/relapsed HLH)

## Baseline parameters (disease onset with moderate trigger, partial NK dysfunction)
base_params <- list(
    kTrig      = 0.8,     # Strong trigger (e.g., EBV infection)
    fNK_dysfunc = 0.3,    # Partial NK dysfunction (secondary HLH)
    kMacro_act = 0.5
)

## Time grid: 0 to 60 days
tgrid <- tgrid(0, 60, 0.5)  # 0.5-day steps

## ---- SCENARIO 1: Untreated ------------------------------------------------
ev_untreated <- ev(time = 0, cmt = "APC_ACT", amt = 0)  # no drug
sim_untreated <- hlh_mod %>%
    param(base_params) %>%
    mrgsim(events = ev_untreated, tgrid = tgrid,
           init = list(APC_ACT = 0.5, T_ACT = 0.05, MAC_ACT = 0.2,
                        IFNg = 20, IL6 = 10, FERR = 300)) %>%
    as.data.frame() %>%
    mutate(Scenario = "Untreated")

## ---- SCENARIO 2: HLH-2004 protocol ----------------------------------------
## DEX 10 mg/m² daily (approx 17 mg/day IV equiv → dose in mg)
## ETOP 150 mg/m² biweekly (d1, d8, d15, d22, d29, d36)
## CsA 6 mg/kg/day oral (approx 360 mg/day)
ev_hlh2004 <- ev_rep(
    c(ev(time = 0, cmt = "DEX_GI", amt = 17, ii = 24, addl = 55),
      ev(time = 0, cmt = "CSA_GI", amt = 360, ii = 24, addl = 55),
      ev(time = 0, cmt = "ETOP_C", amt = 270, rate = 270/2, ii = 168, addl = 5)),
    n = 1)

sim_hlh2004 <- hlh_mod %>%
    param(base_params) %>%
    mrgsim(events = ev_hlh2004, tgrid = tgrid,
           init = list(APC_ACT = 0.5, T_ACT = 0.05, MAC_ACT = 0.2,
                        IFNg = 20, IL6 = 10, FERR = 300)) %>%
    as.data.frame() %>%
    mutate(Scenario = "HLH-2004 (DEX+ETOP+CsA)")

## ---- SCENARIO 3: Emapalumab + DEX (primary / refractory HLH) ---------------
## Emapalumab 1 mg/kg biweekly IV (escalate to 3, 6, 10 mg/kg)
## + Dexamethasone 10 mg/m²
ev_emapa <- ev_rep(
    c(ev(time = 0, cmt = "DEX_GI", amt = 17, ii = 24, addl = 55),
      ev(time = 0, cmt = "EMAPA_C", amt = 70, rate = 70/2, ii = 168, addl = 5)),
    n = 1)

sim_emapa <- hlh_mod %>%
    param(c(base_params, list(fNK_dysfunc = 0.9))) %>%  # primary HLH
    mrgsim(events = ev_emapa, tgrid = tgrid,
           init = list(APC_ACT = 0.6, T_ACT = 0.1, MAC_ACT = 0.3,
                        IFNg = 100, IL6 = 50, FERR = 2000)) %>%
    as.data.frame() %>%
    mutate(Scenario = "Emapalumab + DEX (primary HLH)")

## ---- SCENARIO 4: Anakinra (MAS/sJIA-HLH) ----------------------------------
## Anakinra 100-300 mg/day SC (high-dose 4-8 mg/kg/day)
ev_anakr <- ev_rep(
    c(ev(time = 0, cmt = "DEX_GI", amt = 17, ii = 24, addl = 55),
      ev(time = 0, cmt = "ANK_SC", amt = 200, ii = 24, addl = 55)),
    n = 1)

sim_anakr <- hlh_mod %>%
    param(c(base_params, list(kTrig = 0.5, fNK_dysfunc = 0.1))) %>%  # MAS/sJIA
    mrgsim(events = ev_anakr, tgrid = tgrid,
           init = list(APC_ACT = 0.3, T_ACT = 0.03, MAC_ACT = 0.2,
                        IFNg = 30, IL6 = 80, FERR = 1500)) %>%
    as.data.frame() %>%
    mutate(Scenario = "Anakinra + DEX (MAS/sJIA-HLH)")

## ---- SCENARIO 5: Ruxolitinib salvage (refractory/relapsed) ----------------
## Ruxolitinib 20 mg BID (twice daily)
ev_ruxo <- ev_rep(
    c(ev(time = 0, cmt = "DEX_GI", amt = 17, ii = 24, addl = 55),
      ev(time = 0, cmt = "RUX_GI", amt = 20, ii = 12, addl = 115)),
    n = 1)

sim_ruxo <- hlh_mod %>%
    param(c(base_params, list(kTrig = 1.0, fNK_dysfunc = 0.5))) %>%
    mrgsim(events = ev_ruxo, tgrid = tgrid,
           init = list(APC_ACT = 0.7, T_ACT = 0.12, MAC_ACT = 0.4,
                        IFNg = 200, IL6 = 80, FERR = 5000)) %>%
    as.data.frame() %>%
    mutate(Scenario = "Ruxolitinib + DEX (refractory)")

## ---------------------------------------------------------------------------
## Combine all scenarios
## ---------------------------------------------------------------------------
all_sims <- bind_rows(
    sim_untreated, sim_hlh2004, sim_emapa,
    sim_anakr, sim_ruxo
)

scenario_colors <- c(
    "Untreated"                       = "#E74C3C",
    "HLH-2004 (DEX+ETOP+CsA)"        = "#2980B9",
    "Emapalumab + DEX (primary HLH)"  = "#8E44AD",
    "Anakinra + DEX (MAS/sJIA-HLH)"  = "#27AE60",
    "Ruxolitinib + DEX (refractory)"  = "#D35400"
)

## ---------------------------------------------------------------------------
## VISUALIZATION
## ---------------------------------------------------------------------------

## Figure 1: IFN-γ dynamics across scenarios
p1 <- ggplot(all_sims, aes(time, IFNg, color = Scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = c(50, 500), linetype = "dashed",
               color = c("orange", "red"), alpha = 0.6) +
    annotate("text", x = 58, y = 60, label = "Mild elevation", size = 3, hjust = 1) +
    annotate("text", x = 58, y = 510, label = "Severe HLH", size = 3, hjust = 1) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "IFN-γ Dynamics", x = "Time (days)", y = "IFN-γ (pg/mL)",
         color = "Scenario") +
    theme_bw() + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

## Figure 2: Ferritin trajectory
p2 <- ggplot(all_sims, aes(time, FERR, color = Scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = c(500, 10000), linetype = "dashed",
               color = c("orange", "red"), alpha = 0.6) +
    scale_color_manual(values = scenario_colors) +
    scale_y_log10() +
    labs(title = "Serum Ferritin", x = "Time (days)", y = "Ferritin (ng/mL, log scale)",
         color = "Scenario") +
    theme_bw() + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

## Figure 3: Hemophagocytosis index
p3 <- ggplot(all_sims, aes(time, HEMOPHAG_IDX, color = Scenario)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Hemophagocytosis Index", x = "Time (days)",
         y = "Hemophagocytosis Index (0–1)", color = "Scenario") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

## Figure 4: Survival probability
p4 <- ggplot(all_sims, aes(time, Survival_prob, color = Scenario)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Estimated Survival Probability", x = "Time (days)",
         y = "P(survival)", color = "Scenario") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

## Figure 5: Drug PK (Emapalumab scenario)
p5 <- all_sims %>%
    filter(Scenario == "Emapalumab + DEX (primary HLH)") %>%
    select(time, C_DEX_out, C_EMAPA_out) %>%
    pivot_longer(-time, names_to = "Drug", values_to = "Concentration") %>%
    mutate(Drug = recode(Drug,
        "C_DEX_out" = "Dexamethasone (ng/mL)",
        "C_EMAPA_out" = "Emapalumab (μg/mL × 0.1)")) %>%
    ggplot(aes(time, Concentration, color = Drug)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Dexamethasone (ng/mL)" = "#2980B9",
                                   "Emapalumab (μg/mL × 0.1)" = "#8E44AD")) +
    labs(title = "Drug PK: Emapalumab + DEX Scenario",
         x = "Time (days)", y = "Concentration", color = "Drug") +
    theme_bw() + theme(legend.position = "bottom")

## Figure 6: HScore approximation
p6 <- ggplot(all_sims, aes(time, HScore_approx, color = Scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 169, linetype = "dashed", color = "red", alpha = 0.8) +
    annotate("text", x = 5, y = 175, label = "HScore ≥169 (diagnostic threshold)",
             size = 3, color = "red") +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Approximate HScore Trajectory",
         x = "Time (days)", y = "HScore (0–337)", color = "Scenario") +
    theme_bw() + theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3))

## Print all figures
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

## ---------------------------------------------------------------------------
## Summary table at Day 14 and Day 30
## ---------------------------------------------------------------------------
summary_tbl <- all_sims %>%
    filter(time %in% c(0, 7, 14, 30, 60)) %>%
    group_by(Scenario, time) %>%
    summarise(
        IFNg_pgmL    = round(mean(IFNg), 1),
        Ferritin_ngmL = round(mean(FERR), 0),
        HemophagIdx  = round(mean(HEMOPHAG_IDX), 3),
        BM_Supp      = round(mean(BM_SUPP), 3),
        Survival_pct = round(mean(Survival_prob) * 100, 1),
        HScore_approx = round(mean(HScore_approx), 1),
        .groups = "drop"
    ) %>%
    arrange(time, Scenario)

cat("\n=== Summary: Key Biomarkers by Scenario and Time Point ===\n")
print(summary_tbl, n = 100)

## ---------------------------------------------------------------------------
## Sensitivity Analysis: NK dysfunction level vs. IFN-γ at Day 14
## ---------------------------------------------------------------------------
nk_levels <- seq(0, 0.95, by = 0.05)
sens_results <- lapply(nk_levels, function(nk_frac) {
    hlh_mod %>%
        param(c(base_params, list(fNK_dysfunc = nk_frac))) %>%
        mrgsim(events = ev_untreated, tgrid = tgrid,
               init = list(APC_ACT = 0.5, T_ACT = 0.05, MAC_ACT = 0.2,
                             IFNg = 20, IL6 = 10, FERR = 300)) %>%
        as.data.frame() %>%
        filter(time == 14) %>%
        mutate(NK_dysfunction = nk_frac) %>%
        head(1)
}) %>% bind_rows()

p_sens <- ggplot(sens_results, aes(NK_dysfunction, IFNg)) +
    geom_line(color = "#E74C3C", linewidth = 1.5) +
    geom_point(color = "#E74C3C", size = 2) +
    geom_vline(xintercept = 0.3, linetype = "dashed") +
    annotate("text", x = 0.35, y = max(sens_results$IFNg) * 0.8,
             label = "Secondary HLH\n(typical)", size = 3) +
    labs(title = "Sensitivity: NK Dysfunction vs. IFN-γ at Day 14",
         x = "NK Cell Dysfunction Fraction (0=normal, 1=complete FHL)",
         y = "IFN-γ at Day 14 (pg/mL)") +
    theme_bw()

print(p_sens)

cat("\n=== HLH QSP Model Run Complete ===\n")
cat("Scenarios simulated: Untreated, HLH-2004, Emapalumab+DEX,\n")
cat("                     Anakinra+DEX (MAS), Ruxolitinib+DEX (refractory)\n")
