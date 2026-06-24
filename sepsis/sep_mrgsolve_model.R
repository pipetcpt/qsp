## =============================================================================
## Sepsis & Septic Shock — QSP mrgsolve Model
## =============================================================================
## Disease    : Sepsis / Septic Shock (Gram-negative bacteremia prototype)
## Framework  : mrgsolve ODE (R)
## Compartments: 24 ODE compartments
## Scenarios  : 7 treatment scenarios
## References : Singer 2016 (Sepsis-3), Rivers 2001 (EGDT), ADRENAL 2018,
##              REMAP-CAP 2021, ATTACC/ACTIV 2021, Kumar 2006 (antibiotic timing)
## Author     : Claude Code Routine (CCR) — 2026-06-24
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL DEFINITION
## ─────────────────────────────────────────────────────────────────────────────
code <- '
$PROB
Sepsis & Septic Shock QSP Model
24-compartment ODE model integrating:
  - Bacterial PK dynamics (pathogen burden)
  - 2-cmpt antibiotic PK (meropenem prototype)
  - Cytokine cascade: TNFalpha, IL6, IL10, IL1beta
  - Innate immune: Neutrophils (blood/tissue), Macrophages
  - Complement effector (C5a proxy)
  - Coagulation: Thrombin, Fibrin, PAI1
  - Endothelial damage index
  - 6-domain SOFA (Lung, Renal, Liver, CNS, Cardio, Coag)
  - Vasopressor (norepinephrine) PK/PD
  - Corticosteroid (hydrocortisone) PK/PD
  - IL-6R inhibitor (tocilizumab) PK/PD
  - Lactate & Metabolic stress

$PARAM
// ── Bacterial dynamics ────────────────────────────────────────────────────────
kgrow    = 1.2    // bacterial growth rate (1/h), ~doubling 35 min
kdeath0  = 0.05   // spontaneous bacterial death (1/h)
Bmax     = 1e8    // maximum bacterial load (CFU/mL proxy)
B0       = 1e5    // initial inoculum (CFU/mL)

// ── Antibiotic PK (meropenem prototype, 1g q8h) ────────────────────────────
CL_abx   = 10.0   // clearance (L/h)
V1_abx   = 15.0   // central volume (L)
V2_abx   = 20.0   // peripheral volume (L)
Q_abx    = 8.0    // intercompartmental CL (L/h)
MIC      = 0.5    // minimum inhibitory concentration (mg/L)
Emax_abx = 0.95   // maximum kill rate fraction
EC50_abx = 2.0    // EC50 for kill effect (mg/L above MIC)
kmax_abx = 5.0    // maximum antibiotic kill rate (1/h)

// ── Cytokine parameters ────────────────────────────────────────────────────
// TNFalpha
kprod_TNF = 0.8   // production (ng/mL/h stimulated)
kdeg_TNF  = 0.6   // degradation (1/h, t1/2 ~70 min)
EC50_TNF  = 1e4   // bacteria CFU for half-max TNF induction
kAmp_TNF  = 0.3   // TNF autocrine amplification
TNF0      = 0.02  // baseline TNF (ng/mL)

// IL-6
kprod_IL6 = 1.2
kdeg_IL6  = 0.15  // t1/2 ~4.6 h
EC50_IL6  = 50.0  // TNF for IL-6 induction
IL6_0     = 0.01

// IL-10 (anti-inflammatory)
kprod_IL10 = 0.5
kdeg_IL10  = 0.2
EC50_IL10  = 5.0  // TNF for IL-10 induction
IL10_0     = 0.02

// IL-1beta
kprod_IL1b = 0.6
kdeg_IL1b  = 0.8  // t1/2 ~50 min
EC50_IL1b  = 2e4
IL1b_0     = 0.01

// ── Neutrophils ────────────────────────────────────────────────────────────
Neut_blood0 = 5000  // baseline circulating neutrophils (cells/uL)
krecruit_N  = 0.4   // recruitment rate to tissue
kmarg_N     = 0.2   // margination from blood to tissue
krestore_N  = 0.3   // restoration toward baseline
kprod_N     = 200   // bone marrow output (cells/uL/h)
kdeath_N    = 0.08  // tissue neutrophil death

// Macrophages
Mac0        = 100   // baseline activated macrophages (arbitrary)
kact_Mac    = 0.5   // LPS/cytokine activation rate
kdact_Mac   = 0.1   // deactivation rate
kprod_Mac   = 10    // monocyte-derived production

// ── Complement (C5a effector) ──────────────────────────────────────────────
kprod_C5a  = 0.3
kdeg_C5a   = 1.2   // rapid half-life (~10 min half-life)
C5a_0      = 0.1

// ── Coagulation ────────────────────────────────────────────────────────────
kprod_Thr  = 0.4   // thrombin generation rate
kdeg_Thr   = 0.5
kprod_Fib  = 0.2   // fibrin deposition
kdeg_Fib   = 0.1
kprod_PAI  = 0.3   // PAI-1 production
kdeg_PAI   = 0.2
Thrombin0  = 1.0
Fibrin0    = 1.0
PAI1_0     = 1.0

// ── Endothelial damage ─────────────────────────────────────────────────────
kdam_End  = 0.15   // endothelial damage rate (cytokine-driven)
krep_End  = 0.05   // endothelial repair rate
Endot0    = 0.0    // initial damage (0 = intact)

// ── SOFA sub-scores dynamics ───────────────────────────────────────────────
// Lung (PaO2/FiO2 proxy): baseline 400 mmHg
PF_base   = 400
kdam_lung = 0.02
krep_lung = 0.01

// Renal (Creatinine): baseline 0.9 mg/dL
Cr_base   = 0.9
kprod_Cr  = 0.05
kclr_Cr   = 0.04

// Bilirubin: baseline 0.5 mg/dL
Bil_base  = 0.5
kprod_Bil = 0.03
kclr_Bil  = 0.04

// Lactate (mmol/L): baseline 1.0
Lac_base  = 1.0
kprod_Lac = 0.08
kclr_Lac  = 0.12

// MAP (mmHg): baseline 90
MAP_base  = 90
kdam_MAP  = 0.25   // cytokine-driven vasodilation
krep_MAP  = 0.15   // autoregulatory recovery

// Platelet (×10^9/L): baseline 250
Plt_base  = 250
kcons_Plt = 0.04   // consumptive drop (DIC)
kregen_Plt= 0.02   // megakaryocyte regeneration

// ── Vasopressor PK (norepinephrine, mcg/kg/min → simplified PK) ───────────
CL_NE    = 120.0   // clearance (L/h) — very rapid
V1_NE    = 18.0    // volume (L)
Emax_NE  = 30.0    // max MAP increase (mmHg)
EC50_NE  = 0.05    // plasma conc for half-max (mcg/mL equiv)

// ── Corticosteroid PK/PD (hydrocortisone 200 mg/day = 8.33 mg/h) ──────────
CL_HC    = 25.0    // clearance (L/h)
V1_HC    = 40.0    // volume (L)
Emax_HC  = 0.65    // max cytokine suppression (fraction)
EC50_HC  = 200.0   // conc for half-max (ng/mL)

// ── IL-6R inhibitor PK (tocilizumab 8 mg/kg IV) ───────────────────────────
CL_Toci  = 0.28    // clearance (L/h)
V1_Toci  = 4.2     // volume (L)
Emax_Toci = 0.90
EC50_Toci = 1.5    // mcg/mL

// ── Treatment switches (0 = off, 1 = on) ──────────────────────────────────
useAbx    = 0
useNE     = 0
useHC     = 0
useToci   = 0
useFluid  = 0

// ── Fluid resuscitation (MAP boost, L in first 3h) ─────────────────────────
FluidBoost = 0.0  // L given (activates early MAP restoration)

$CMT
// Bacterial burden
BACT

// Antibiotic PK
ABX1      // central compartment (mg)
ABX2      // peripheral compartment (mg)

// Cytokines (ng/mL)
TNF
IL6
IL10
IL1B

// Immune cells
NEUT_B    // blood neutrophils (cells/uL)
NEUT_T    // tissue neutrophils
MACS      // activated macrophages

// Complement
C5A

// Coagulation
THROMBIN
FIBRIN
PAI1

// Endothelial damage index (0–1)
ENDOT

// Organ function endpoints
PF_RATIO   // PaO2/FiO2 (mmHg)
CREATININE // mg/dL
BILIRUBIN  // mg/dL
LACTATE    // mmol/L
MAP_val    // mmHg
PLT_COUNT  // ×10^9/L

// Vasopressor plasma (mcg/mL equiv)
NE_C

// Corticosteroid plasma (ng/mL)
HC_C

// Tocilizumab plasma (mcg/mL)
TOCI_C

$INIT
BACT      = 1e5
ABX1      = 0
ABX2      = 0
TNF       = 0.02
IL6       = 0.01
IL10      = 0.02
IL1B      = 0.01
NEUT_B    = 5000
NEUT_T    = 200
MACS      = 10
C5A       = 0.1
THROMBIN  = 1.0
FIBRIN    = 1.0
PAI1      = 1.0
ENDOT     = 0.0
PF_RATIO  = 400
CREATININE= 0.9
BILIRUBIN = 0.5
LACTATE   = 1.0
MAP_val   = 90.0
PLT_COUNT = 250.0
NE_C      = 0.0
HC_C      = 0.0
TOCI_C    = 0.0

$MAIN
// ── Antibiotic concentration (mg/L) ──────────────────────────────────────
double Cp_abx = ABX1 / V1_abx;

// Antibiotic kill: sigmoidal Emax above MIC
double dC_abx = (Cp_abx > MIC) ? (Cp_abx - MIC) : 0.0;
double Ekill   = kmax_abx * pow(dC_abx, 2.0) / (pow(EC50_abx, 2.0) + pow(dC_abx, 2.0));

// ── IL-10 inhibition of pro-inflammatory cytokines ─────────────────────────
double Inh_IL10 = 1.0 / (1.0 + IL10 / 5.0);  // Hill inhibition

// ── Hydrocortisone effect ──────────────────────────────────────────────────
double E_HC = Emax_HC * HC_C / (EC50_HC + HC_C);

// ── Tocilizumab effect on IL-6 signaling ──────────────────────────────────
double E_Toci = Emax_Toci * TOCI_C / (EC50_Toci + TOCI_C);

// ── Vasopressor MAP effect ──────────────────────────────────────────────────
double E_NE = Emax_NE * NE_C / (EC50_NE + NE_C);

// ── Combined cytokine suppression (HC + IL-10) ────────────────────────────
double Inh_combined = (1.0 - E_HC) * Inh_IL10;

// ── Bacterial driven NF-kB signal ────────────────────────────────────────
double NFkB_signal = BACT / (BACT + 1e5);

// ── Cytokine inflammation index (drives organ damage) ─────────────────────
double CytokineIndex = (TNF / 10.0 + IL6 / 200.0 + IL1B / 5.0) / 3.0;
double CI_norm = CytokineIndex / (1.0 + CytokineIndex);

$ODE
// ─── BACTERIAL DYNAMICS ─────────────────────────────────────────────────────
double BactLogistic = kgrow * BACT * (1.0 - BACT / Bmax);
double BactPhago    = 0.002 * NEUT_T * BACT / (BACT + 1e4);  // neutrophil kill
double BactAbxKill  = useAbx * Ekill * BACT;
dxdt_BACT = BactLogistic - kdeath0 * BACT - BactPhago - BactAbxKill;
if (BACT < 1.0) dxdt_BACT = -BACT;  // floor

// ─── ANTIBIOTIC PK ──────────────────────────────────────────────────────────
// Dosing handled via ADDL/II in event table
dxdt_ABX1 = -CL_abx * Cp_abx - Q_abx * (ABX1/V1_abx - ABX2/V2_abx);
dxdt_ABX2 =  Q_abx * (ABX1/V1_abx - ABX2/V2_abx);

// ─── CYTOKINES ───────────────────────────────────────────────────────────────
// TNFalpha
double Stim_TNF = kprod_TNF * NFkB_signal * Inh_combined
                  + kAmp_TNF * TNF * Inh_combined;
dxdt_TNF = Stim_TNF - kdeg_TNF * (TNF - TNF0);
if (TNF < 0) dxdt_TNF = 0;

// IL-6
double Stim_IL6 = kprod_IL6 * TNF / (EC50_IL6 + TNF) * Inh_combined * (1.0 - E_Toci);
dxdt_IL6 = Stim_IL6 - kdeg_IL6 * (IL6 - IL6_0);
if (IL6 < 0) dxdt_IL6 = 0;

// IL-10 (anti-inflammatory, induced by TNF and macrophages)
double Stim_IL10 = kprod_IL10 * TNF / (EC50_IL10 + TNF);
dxdt_IL10 = Stim_IL10 - kdeg_IL10 * (IL10 - IL10_0);
if (IL10 < 0) dxdt_IL10 = 0;

// IL-1beta
double Stim_IL1B = kprod_IL1b * NFkB_signal * Inh_combined;
dxdt_IL1B = Stim_IL1B - kdeg_IL1b * (IL1B - IL1b_0);
if (IL1B < 0) dxdt_IL1B = 0;

// ─── INNATE IMMUNE CELLS ────────────────────────────────────────────────────
// Blood neutrophils (recruitment to tissue driven by IL-8 proxy = IL6)
double NeutRecruit = kmarg_N * NEUT_B * IL6 / (50.0 + IL6);
double NeutProd    = kprod_N * (1.0 + 2.0 * NFkB_signal); // emergency granulopoiesis
dxdt_NEUT_B = NeutProd - kdeath_N * NEUT_B - NeutRecruit
              + krestore_N * (Neut_blood0 - NEUT_B) * (BACT < 100 ? 1.0 : 0.0);

// Tissue neutrophils
dxdt_NEUT_T = NeutRecruit - kdeath_N * NEUT_T * (1.0 + CI_norm);

// Activated macrophages
double MacAct = kact_Mac * NFkB_signal * (1.0 - E_HC);
dxdt_MACS = MacAct * Mac0 + kprod_Mac - kdact_Mac * MACS;

// ─── COMPLEMENT (C5a proxy) ──────────────────────────────────────────────────
dxdt_C5A = kprod_C5a * NFkB_signal - kdeg_C5a * (C5A - C5a_0);

// ─── COAGULATION ─────────────────────────────────────────────────────────────
// Thrombin (driven by TF expression from TNF + endothelial damage)
double ThrombStim = kprod_Thr * CI_norm * (1.0 + ENDOT);
dxdt_THROMBIN = ThrombStim - kdeg_Thr * (THROMBIN - Thrombin0);

// Fibrin (driven by thrombin)
double FibrinForm = kprod_Fib * THROMBIN;
double FibrinLysis = kdeg_Fib * (1.0 / PAI1); // tPA activity (inhibited by PAI-1)
dxdt_FIBRIN = FibrinForm - FibrinLysis * FIBRIN;

// PAI-1 (driven by TNF)
double PAI_stim = kprod_PAI * TNF / (5.0 + TNF);
dxdt_PAI1 = PAI_stim - kdeg_PAI * (PAI1 - PAI1_0);

// ─── ENDOTHELIAL DAMAGE (0 = intact, 1 = full failure) ───────────────────────
double EndDam  = kdam_End * CI_norm * (1.0 - ENDOT);
double EndRep  = krep_End * ENDOT * (1.0 - CI_norm);
dxdt_ENDOT = EndDam - EndRep;
if (ENDOT > 1.0) dxdt_ENDOT = 0;
if (ENDOT < 0.0) dxdt_ENDOT = 0;

// ─── ORGAN ENDPOINTS ─────────────────────────────────────────────────────────
// PaO2/FiO2 ratio (mmHg) — decreases with ARDS
double PF_dam = kdam_lung * CI_norm * ENDOT * PF_RATIO;
double PF_rep = krep_lung * (PF_base - PF_RATIO) * (CI_norm < 0.2 ? 1.0 : 0.0);
dxdt_PF_RATIO = PF_rep - PF_dam;
if (PF_RATIO < 50) dxdt_PF_RATIO = 0;

// Creatinine (mg/dL)
double CrProd = kprod_Cr * (1.0 + 3.0 * CI_norm * ENDOT);
double CrClr  = kclr_Cr * CREATININE;
dxdt_CREATININE = CrProd - CrClr;

// Bilirubin (mg/dL)
double BilProd = kprod_Bil * (1.0 + 2.0 * CI_norm);
double BilClr  = kclr_Bil * BILIRUBIN;
dxdt_BILIRUBIN = BilProd - BilClr;

// Lactate (mmol/L)
double LacProd = kprod_Lac * CI_norm * (1.0 + ENDOT);
double LacClr  = kclr_Lac * LACTATE;
dxdt_LACTATE = LacProd - LacClr;

// MAP (mmHg) — vasodilatory shock from NO/cytokines
double MAP_dam = kdam_MAP * CI_norm * ENDOT;
double MAP_NE_eff = useNE * E_NE;
double MAP_HC_eff = useHC * 0.3 * E_HC;           // steroid vasopressor sensitization
double MAP_fluid  = useFluid * FluidBoost * 2.0;   // simplified fluid effect
double MAP_rep    = krep_MAP * (MAP_base - MAP_val) * 0.1;
dxdt_MAP_val = MAP_NE_eff + MAP_fluid + MAP_rep - MAP_dam;
if (MAP_val < 30) dxdt_MAP_val = 0;

// Platelet count (×10^9/L) — consumptive with DIC
double PltCons = kcons_Plt * FIBRIN * PLT_COUNT / 200.0;
double PltRegen= kregen_Plt * (Plt_base - PLT_COUNT) * 0.1;
dxdt_PLT_COUNT = PltRegen - PltCons;
if (PLT_COUNT < 10) dxdt_PLT_COUNT = 0;

// ─── VASOPRESSOR PK (1-cmpt, infusion via event table) ──────────────────────
dxdt_NE_C = -CL_NE / V1_NE * NE_C;

// ─── CORTICOSTEROID PK (1-cmpt) ──────────────────────────────────────────────
dxdt_HC_C = -CL_HC / V1_HC * HC_C;

// ─── TOCILIZUMAB PK (1-cmpt, single or repeat dose) ─────────────────────────
dxdt_TOCI_C = -CL_Toci / V1_Toci * TOCI_C;

$TABLE
// ── SOFA score calculation (Sepsis-3, range 0–24) ──────────────────────────

// Lung SOFA (PaO2/FiO2)
int sofa_lung;
if (PF_RATIO >= 400)       sofa_lung = 0;
else if (PF_RATIO >= 300)  sofa_lung = 1;
else if (PF_RATIO >= 200)  sofa_lung = 2;
else if (PF_RATIO >= 100)  sofa_lung = 3;
else                        sofa_lung = 4;

// Renal SOFA (Creatinine)
int sofa_renal;
if (CREATININE < 1.2)      sofa_renal = 0;
else if (CREATININE < 2.0) sofa_renal = 1;
else if (CREATININE < 3.5) sofa_renal = 2;
else if (CREATININE < 5.0) sofa_renal = 3;
else                        sofa_renal = 4;

// Liver SOFA (Bilirubin)
int sofa_liver;
if (BILIRUBIN < 1.2)       sofa_liver = 0;
else if (BILIRUBIN < 2.0)  sofa_liver = 1;
else if (BILIRUBIN < 6.0)  sofa_liver = 2;
else if (BILIRUBIN < 12.0) sofa_liver = 3;
else                        sofa_liver = 4;

// Cardiovascular SOFA (MAP + vasopressor)
int sofa_cardio;
if (MAP_val >= 70)          sofa_cardio = 0;
else if (MAP_val >= 65)     sofa_cardio = 1;
else if (useNE == 1.0 && NE_C > 0.005) sofa_cardio = 3;
else                        sofa_cardio = 2;

// Coagulation SOFA (Platelets)
int sofa_coag;
if (PLT_COUNT >= 150)       sofa_coag = 0;
else if (PLT_COUNT >= 100)  sofa_coag = 1;
else if (PLT_COUNT >= 50)   sofa_coag = 2;
else if (PLT_COUNT >= 20)   sofa_coag = 3;
else                        sofa_coag = 4;

// CNS SOFA (proxy from lactate/encephalopathy marker)
int sofa_cns;
double enc_marker = LACTATE + CI_norm * 3.0;
if (enc_marker < 1.5)       sofa_cns = 0;
else if (enc_marker < 2.5)  sofa_cns = 1;
else if (enc_marker < 4.0)  sofa_cns = 2;
else if (enc_marker < 6.0)  sofa_cns = 3;
else                        sofa_cns = 4;

// Total SOFA
int SOFA = sofa_lung + sofa_renal + sofa_liver + sofa_cardio + sofa_coag + sofa_cns;

// 28-day mortality probability (log-odds model: Seymour 2017 JAMA)
double SOFA_logit = -6.5 + 0.45 * SOFA;
double Mortality_Prob = 1.0 / (1.0 + exp(-SOFA_logit));

// Shock criteria
double shock_flag = ((MAP_val < 65.0) && (LACTATE > 2.0)) ? 1.0 : 0.0;

// Bacteria exposure (log CFU)
double log_BACT = (BACT > 0) ? log10(BACT) : 0.0;

double TNF_ng   = TNF;
double IL6_pg   = IL6 * 1000.0;  // convert ng/mL → pg/mL for reporting
double Cp_abx_report = ABX1 / V1_abx;
double CytokineIdx = (TNF/10.0 + IL6/200.0 + IL1B/5.0)/3.0;

$CAPTURE
BACT ABX1 ABX2 TNF IL6 IL10 IL1B
NEUT_B NEUT_T MACS C5A
THROMBIN FIBRIN PAI1 ENDOT
PF_RATIO CREATININE BILIRUBIN LACTATE MAP_val PLT_COUNT
NE_C HC_C TOCI_C
SOFA sofa_lung sofa_renal sofa_liver sofa_cardio sofa_coag sofa_cns
Mortality_Prob shock_flag log_BACT
TNF_ng IL6_pg Cp_abx_report CytokineIdx
'

## ─────────────────────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("sepsis_qsp", code)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: BUILD EVENT TABLE
## ─────────────────────────────────────────────────────────────────────────────
make_events <- function(
    use_abx   = FALSE, abx_dose = 1000, abx_ii = 8, abx_addl = 83,  # mg, q8h × 28d
    use_ne    = FALSE, ne_rate  = 0.002,                              # mcg/mL/h equiv
    use_hc    = FALSE, hc_dose  = 50,   hc_ii  = 6,  hc_addl = 111, # mg q6h × 7d
    use_toci  = FALSE, toci_dose = 672,                               # mg = 8 mg/kg × 84 kg
    use_fluid = FALSE, fluid_dose = 30                                # mL/kg = 2.1L equiv
) {
  ev <- ev()

  if (use_abx)
    ev <- ev + ev(amt = abx_dose, cmt = "ABX1", ii = abx_ii, addl = abx_addl, time = 0)

  if (use_ne)
    ev <- ev + ev(amt = ne_rate * 18, cmt = "NE_C", rate = ne_rate, time = 2, to = 72)

  if (use_hc)
    ev <- ev + ev(amt = hc_dose, cmt = "HC_C", ii = hc_ii, addl = hc_addl, time = 1)

  if (use_toci)
    ev <- ev + ev(amt = toci_dose, cmt = "TOCI_C", time = 6)

  if (use_fluid)
    ev <- ev + ev(amt = fluid_dose * 1000, cmt = "MAP_val", time = 0, evid = 2)  # bolus effect

  ev
}

## ─────────────────────────────────────────────────────────────────────────────
## 7 TREATMENT SCENARIOS
## ─────────────────────────────────────────────────────────────────────────────
scenarios <- list(
  S1_NoTx = list(
    label  = "S1: No Treatment (Natural History)",
    color  = "#C0392B",
    params = list(useAbx=0, useNE=0, useHC=0, useToci=0, useFluid=0),
    events = make_events()
  ),
  S2_AbxOnly = list(
    label  = "S2: Antibiotics Only",
    color  = "#E67E22",
    params = list(useAbx=1, useNE=0, useHC=0, useToci=0, useFluid=0),
    events = make_events(use_abx = TRUE)
  ),
  S3_AbxNE = list(
    label  = "S3: Antibiotics + Norepinephrine",
    color  = "#F1C40F",
    params = list(useAbx=1, useNE=1, useHC=0, useToci=0, useFluid=0),
    events = make_events(use_abx = TRUE, use_ne = TRUE)
  ),
  S4_AbxNEFluid = list(
    label  = "S4: Antibiotics + NE + Fluid (Bundle)",
    color  = "#27AE60",
    params = list(useAbx=1, useNE=1, useHC=0, useToci=0, useFluid=1, FluidBoost=3),
    events = make_events(use_abx = TRUE, use_ne = TRUE, use_fluid = TRUE)
  ),
  S5_AbxNEFluidHC = list(
    label  = "S5: Bundle + Hydrocortisone",
    color  = "#2980B9",
    params = list(useAbx=1, useNE=1, useHC=1, useToci=0, useFluid=1, FluidBoost=3),
    events = make_events(use_abx = TRUE, use_ne = TRUE, use_hc = TRUE, use_fluid = TRUE)
  ),
  S6_BundleToci = list(
    label  = "S6: Bundle + HC + Tocilizumab",
    color  = "#8E44AD",
    params = list(useAbx=1, useNE=1, useHC=1, useToci=1, useFluid=1, FluidBoost=3),
    events = make_events(use_abx = TRUE, use_ne = TRUE, use_hc = TRUE,
                         use_toci = TRUE, use_fluid = TRUE)
  ),
  S7_ImmunoSupp = list(
    label  = "S7: Immunocompromised (Abx + NE only, high inoculum)",
    color  = "#7F8C8D",
    params = list(useAbx=1, useNE=1, useHC=0, useToci=0, useFluid=1,
                  B0=1e6, FluidBoost=2),  # higher inoculum, blunted immune
    events = make_events(use_abx = TRUE, use_ne = TRUE, use_fluid = TRUE)
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## RUN SIMULATIONS
## ─────────────────────────────────────────────────────────────────────────────
sim_time <- seq(0, 168, by = 0.5)  # 0–168 h (7 days)

results_list <- lapply(names(scenarios), function(sname) {
  sc <- scenarios[[sname]]

  # Build parameter list (override model defaults)
  pmod <- mod %>% param(sc$params)

  # Run simulation
  out <- mrgsim(pmod, events = sc$events, tgrid = sim_time, obsonly = TRUE)

  as.data.frame(out) %>%
    mutate(scenario = sc$label, color = sc$color)
})

results_all <- bind_rows(results_list)

## ─────────────────────────────────────────────────────────────────────────────
## PLOTTING FUNCTIONS
## ─────────────────────────────────────────────────────────────────────────────
scenario_colors <- setNames(
  sapply(scenarios, function(s) s$color),
  sapply(scenarios, function(s) s$label)
)

#' Plot a single endpoint across all scenarios
plot_endpoint <- function(data, yvar, ylab, title, ytrans = "identity",
                          hline = NULL, hline_label = NULL) {
  p <- ggplot(data, aes(x = time, y = .data[[yvar]],
                         color = scenario, group = scenario)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_color_manual(values = scenario_colors, name = "Treatment Scenario") +
    scale_x_continuous(
      breaks = seq(0, 168, 24),
      labels = paste0("Day ", seq(0, 7)),
      name   = "Time (Hours from Sepsis Onset)"
    ) +
    scale_y_continuous(trans = ytrans, name = ylab) +
    labs(title = title,
         subtitle = "Sepsis QSP Model — 7-day simulation",
         caption = "mrgsolve ODE model | CCR 2026-06-24") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))

  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                         color = "firebrick", linewidth = 0.8) +
      annotate("text", x = 5, y = hline * 1.05,
               label = hline_label, color = "firebrick", size = 3.5)
  }
  p
}

## Generate key plots
p_bact   <- plot_endpoint(results_all, "log_BACT", "Log₁₀ Bacterial Load (CFU/mL)",
                          "Bacteremia — Pathogen Clearance",
                          hline = 2, hline_label = "Clearance threshold")

p_tnf    <- plot_endpoint(results_all, "TNF_ng", "TNFα (ng/mL)",
                          "TNFα — Cytokine Storm Kinetics")

p_il6    <- plot_endpoint(results_all, "IL6_pg", "IL-6 (pg/mL)",
                          "IL-6 — Acute Phase Response",
                          hline = 100, hline_label = "Elevated (>100 pg/mL)")

p_sofa   <- plot_endpoint(results_all, "SOFA", "SOFA Score (0–24)",
                          "SOFA Score — Organ Failure Index",
                          hline = 11, hline_label = "High mortality risk (SOFA ≥11)")

p_map    <- plot_endpoint(results_all, "MAP_val", "MAP (mmHg)",
                          "Mean Arterial Pressure",
                          hline = 65, hline_label = "Shock threshold (<65 mmHg)")

p_lactate<- plot_endpoint(results_all, "LACTATE", "Lactate (mmol/L)",
                          "Blood Lactate — Tissue Hypoperfusion",
                          hline = 2, hline_label = "Septic shock threshold (>2 mmol/L)")

p_pf     <- plot_endpoint(results_all, "PF_RATIO", "PaO2/FiO2 (mmHg)",
                          "Respiratory Function (ARDS Risk)",
                          hline = 300, hline_label = "Mild ARDS (<300)")

p_mort   <- plot_endpoint(results_all, "Mortality_Prob", "28-day Mortality Probability",
                          "Estimated 28-day Mortality Risk")

## ─────────────────────────────────────────────────────────────────────────────
## SUMMARY TABLE AT 72h
## ─────────────────────────────────────────────────────────────────────────────
summary_72h <- results_all %>%
  filter(abs(time - 72) < 0.6) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, log_BACT, TNF_ng, IL6_pg, SOFA, MAP_val, LACTATE,
         CREATININE, PF_RATIO, Mortality_Prob, shock_flag) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

print(summary_72h)

## ─────────────────────────────────────────────────────────────────────────────
## ANTIBIOTIC PK PROFILE
## ─────────────────────────────────────────────────────────────────────────────
abx_data <- results_all %>%
  filter(scenario == "S2: Antibiotics Only" & time <= 48) %>%
  mutate(Cp_mg_L = Cp_abx_report)

p_abx_pk <- ggplot(abx_data, aes(x = time, y = Cp_mg_L)) +
  geom_line(color = "#E67E22", linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  annotate("text", x = 5, y = 0.7, label = "MIC = 0.5 mg/L", color = "red", size = 3.5) +
  scale_x_continuous(breaks = seq(0, 48, 8), name = "Time (h)") +
  labs(title = "Meropenem PK Profile (1g q8h IV)",
       y = "Plasma Concentration (mg/L)",
       caption = "fT>MIC drives bacterial kill (Emax PD model)") +
  theme_bw(base_size = 12)

cat("\nSepsis QSP model compiled and 7 scenarios simulated successfully.\n")
cat("Use p_bact, p_tnf, p_il6, p_sofa, p_map, p_lactate, p_pf, p_mort to view plots.\n")
cat("summary_72h contains 72h endpoint comparison across scenarios.\n")
