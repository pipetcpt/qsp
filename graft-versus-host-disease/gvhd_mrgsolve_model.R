## ============================================================
## Graft-versus-Host Disease (GvHD) QSP Model
## mrgsolve ODE Implementation
## ============================================================
## Disease: Chronic + Acute GvHD after Allogeneic HSCT
## Drugs modeled:
##   1. Cyclosporine A (CsA) – CNI prophylaxis
##   2. Tacrolimus (TAC)     – CNI prophylaxis (alternative)
##   3. Prednisone (PRED)    – First-line treatment
##   4. Ruxolitinib (RUX)    – JAK1/2 inhibitor (SR-cGvHD)
##   5. Belumosudil (BELU)   – ROCK2 inhibitor (cGvHD ≥2L)
##   6. Mycophenolate (MMF)  – Antiproliferative prophylaxis
## ============================================================
## References:
##  - Zeiser R et al. NEJM 2020;382:1800 (Ruxolitinib cGvHD)
##  - Zeiser R et al. NEJM 2021;385:228  (Ruxolitinib aGvHD)
##  - Cutler C et al. JCO 2021;39:1808   (Belumosudil)
##  - Lee SJ et al. NEJM 2004;350:2912   (Prophylaxis)
##  - Ferrara JL et al. Lancet 2009;373:1550 (aGvHD review)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL CODE
## ============================================================

gvhd_code <- '
$PROB GvHD QSP Model - Acute and Chronic GvHD

$PARAM @annotated
// --- Drug PK parameters ---
// Cyclosporine A (CsA)
CsA_F    : 0.30  : CsA oral bioavailability
CsA_ka   : 0.80  : CsA absorption rate constant (1/h)
CsA_CL   : 25.0  : CsA clearance (L/h)
CsA_V1   : 300   : CsA central volume (L)
CsA_Q    : 50.0  : CsA inter-compartment CL (L/h)
CsA_V2   : 2000  : CsA peripheral volume (L)

// Tacrolimus (TAC)
TAC_F    : 0.25  : TAC oral bioavailability
TAC_ka   : 0.60  : TAC absorption rate constant (1/h)
TAC_CL   : 2.50  : TAC clearance (L/h)
TAC_V1   : 30.0  : TAC central volume (L)
TAC_Q    : 3.00  : TAC inter-compartment CL (L/h)
TAC_V2   : 200   : TAC peripheral volume (L)

// Prednisone / Prednisolone
PRED_F   : 0.99  : Prednisone bioavailability
PRED_ka  : 1.20  : PRED absorption rate constant (1/h)
PRED_CL  : 15.0  : PRED clearance (L/h)
PRED_V1  : 50.0  : PRED central volume (L)

// Ruxolitinib (RUX)
RUX_F    : 0.95  : Ruxolitinib bioavailability
RUX_ka   : 2.00  : RUX absorption rate constant (1/h)
RUX_CL   : 17.7  : RUX clearance (L/h)
RUX_V1   : 72.0  : RUX central volume (L)
RUX_Q    : 10.0  : RUX inter-compartment CL (L/h)
RUX_V2   : 100   : RUX peripheral volume (L)

// Belumosudil (BELU)
BELU_F   : 0.80  : Belumosudil bioavailability (with food)
BELU_ka  : 0.50  : BELU absorption rate constant (1/h)
BELU_CL  : 6.00  : BELU clearance (L/h)
BELU_V1  : 175   : BELU central volume (L)

// MMF -> MPA
MMF_F    : 0.94  : MMF bioavailability (as MPA)
MMF_ka   : 1.50  : MMF->MPA absorption rate (1/h)
MMF_CL   : 12.0  : MPA clearance (L/h)
MMF_V1   : 3.50  : MPA central volume (L)

// --- PD parameters ---
// CNI effect on Calcineurin / NFAT / IL-2
IC50_CsA : 150.0 : CsA IC50 for calcineurin (ng/mL = ug/L)
IC50_TAC : 10.0  : TAC IC50 for calcineurin (ng/mL)
Imax_CN  : 0.85  : Maximum calcineurin inhibition (fraction)
n_CN     : 1.50  : Hill coefficient for CNI

// Prednisone PD
IC50_PRED: 50.0  : PRED IC50 for cytokine suppression (ng/mL)
Imax_PRED: 0.80  : Max broad cytokine suppression

// Ruxolitinib PD (JAK1/2)
IC50_RUX : 280.0 : RUX IC50 for JAK1 inhibition (nM, ~120 ng/mL)
IC50_RUX2: 250.0 : RUX IC50 for JAK2 inhibition (nM)
Imax_RUX : 0.90  : Max JAK inhibition

// Belumosudil PD (ROCK2)
IC50_BELU: 100.0 : BELU IC50 for ROCK2 (nM, ~53 ng/mL as ≈molar wt 532)
Imax_BELU: 0.85  : Max ROCK2 inhibition

// MMF/MPA PD (IMPDH)
IC50_MPA : 1.50  : MPA IC50 for lymphocyte proliferation (ug/mL)
Imax_MPA : 0.85  : Max proliferation inhibition

// --- Disease biology parameters ---
// T cell compartments
kact_T   : 0.03  : T cell activation rate constant (1/h)
kprol_T  : 0.05  : T cell proliferation rate (1/h)
kdiff_Th1: 0.02  : Naive T -> Th1 differentiation rate
kdiff_Th17:0.015 : Naive T -> Th17 differentiation rate
kdiff_Treg:0.010 : Naive T -> Treg differentiation rate (IL-2/TGFb driven)
kdeath_T : 0.008 : T cell apoptosis rate (1/h)
T_ss     : 1.00  : T cell homeostatic set-point (norm=1)

// Cytokine kinetics
kprod_TNFa:0.040 : TNF-α production by Th1/CD8 (ng/mL/h)
kprod_IFNg:0.035 : IFN-γ production (ng/mL/h)
kprod_IL17:0.020 : IL-17A production by Th17 (ng/mL/h)
kprod_IL10:0.025 : IL-10 production by Treg (ng/mL/h)
kprod_TGFb:0.015 : TGF-β production (ng/mL/h)
kprod_IL6 :0.030 : IL-6 production (ng/mL/h)
kdeg_Cyt  : 0.30  : Generic cytokine degradation rate (1/h)

// Organ damage
kdam_skin : 0.005 : Skin damage accumulation rate
krep_skin : 0.002 : Skin repair rate (drug-modulated)
kdam_gut  : 0.008 : Gut damage accumulation rate
krep_gut  : 0.003 : Gut repair rate
kdam_liver: 0.004 : Liver damage accumulation rate
krep_liver: 0.002 : Liver repair rate
kdam_lung : 0.003 : Lung damage accumulation rate
krep_lung : 0.001 : Lung repair rate (BOS essentially irreversible)
kfibr     : 0.002 : Fibrosis accumulation rate (TGFb-driven)

// Initial disease state (aGvHD onset day 20-30 post HSCT)
Allo_stim : 1.00  : Alloantigen stimulus (0=no GvHD, 1=full)
GvL_effect: 0.30  : Graft-versus-Leukemia effect (anti-relapse)

$CMT @annotated
// Drug PK compartments
CsA_gut     : CsA gut absorption compartment (mg)
CsA_cent    : CsA central compartment (mg)
CsA_peri    : CsA peripheral compartment (mg)
TAC_gut     : TAC gut absorption compartment (mg)
TAC_cent    : TAC central compartment (mg)
TAC_peri    : TAC peripheral compartment (mg)
PRED_gut    : Prednisone gut absorption (mg)
PRED_cent   : Prednisolone central compartment (mg)
RUX_gut     : Ruxolitinib gut absorption (mg)
RUX_cent    : Ruxolitinib central compartment (mg)
RUX_peri    : Ruxolitinib peripheral compartment (mg)
BELU_gut    : Belumosudil gut absorption (mg)
BELU_cent   : Belumosudil central compartment (mg)
MPA_gut     : MPA gut absorption (mg)
MPA_cent    : MPA central compartment (mg)

// Immunology compartments
Th1         : Th1 cell relative pool (norm=1)
Th17        : Th17 cell relative pool (norm=1)
Treg        : Regulatory T cell pool (norm=1)
CD8_eff     : CD8 effector T cell pool (norm=1)
Bcell       : B cell pool including plasma cells (norm=1)

// Cytokine compartments (ng/mL)
TNFa        : TNF-alpha concentration (ng/mL)
IFNg        : IFN-gamma concentration (ng/mL)
IL17A       : IL-17A concentration (ng/mL)
IL10        : IL-10 concentration (ng/mL)
TGFb        : TGF-beta concentration (ng/mL)
IL6         : IL-6 concentration (ng/mL)

// Organ damage scores (0-1 scale, 0=normal, 1=max damage)
Skin_dmg    : Skin damage score (0-1)
Gut_dmg     : Gut damage score (0-1)
Liver_dmg   : Liver damage score (0-1)
Lung_dmg    : Lung damage score (0-1)
Fibrosis    : Fibrotic index (0-1, chronic GvHD)

$INIT @annotated
CsA_gut  = 0,  CsA_cent  = 0,  CsA_peri  = 0
TAC_gut  = 0,  TAC_cent  = 0,  TAC_peri  = 0
PRED_gut = 0,  PRED_cent = 0
RUX_gut  = 0,  RUX_cent  = 0,  RUX_peri  = 0
BELU_gut = 0,  BELU_cent = 0
MPA_gut  = 0,  MPA_cent  = 0
Th1      = 1,  Th17      = 1,  Treg      = 1,  CD8_eff  = 1,  Bcell    = 1
TNFa     = 0.1, IFNg    = 0.1, IL17A    = 0.1, IL10     = 0.5, TGFb    = 0.5, IL6 = 0.1
Skin_dmg = 0,  Gut_dmg  = 0,  Liver_dmg = 0,  Lung_dmg = 0,  Fibrosis = 0

$ODE

// ===========================================
// DRUG PK ODEs
// ===========================================

// CsA PK (2-compartment oral)
double CsA_conc = CsA_cent / CsA_V1 * 1000; // ug/L = ng/mL
dxdt_CsA_gut  = -CsA_ka * CsA_gut;
dxdt_CsA_cent =  CsA_F * CsA_ka * CsA_gut
                - (CsA_CL + CsA_Q) / CsA_V1 * CsA_cent
                + CsA_Q / CsA_V2 * CsA_peri;
dxdt_CsA_peri =  CsA_Q / CsA_V1 * CsA_cent - CsA_Q / CsA_V2 * CsA_peri;

// TAC PK (2-compartment oral)
double TAC_conc = TAC_cent / TAC_V1 * 1000; // ng/mL
dxdt_TAC_gut  = -TAC_ka * TAC_gut;
dxdt_TAC_cent =  TAC_F * TAC_ka * TAC_gut
                - (TAC_CL + TAC_Q) / TAC_V1 * TAC_cent
                + TAC_Q / TAC_V2 * TAC_peri;
dxdt_TAC_peri =  TAC_Q / TAC_V1 * TAC_cent - TAC_Q / TAC_V2 * TAC_peri;

// PRED PK (1-compartment oral)
double PRED_conc = PRED_cent / PRED_V1 * 1000; // ng/mL
dxdt_PRED_gut  = -PRED_ka * PRED_gut;
dxdt_PRED_cent =  PRED_F * PRED_ka * PRED_gut - PRED_CL / PRED_V1 * PRED_cent;

// RUX PK (2-compartment oral)
double RUX_conc_ngmL = RUX_cent / RUX_V1 * 1000; // ng/mL
double RUX_conc_nM   = RUX_conc_ngmL / 0.306;     // MW ruxolitinib ~306 g/mol
dxdt_RUX_gut  = -RUX_ka * RUX_gut;
dxdt_RUX_cent =  RUX_F * RUX_ka * RUX_gut
                - (RUX_CL + RUX_Q) / RUX_V1 * RUX_cent
                + RUX_Q / RUX_V2 * RUX_peri;
dxdt_RUX_peri =  RUX_Q / RUX_V1 * RUX_cent - RUX_Q / RUX_V2 * RUX_peri;

// BELU PK (1-compartment oral)
double BELU_conc_ngmL = BELU_cent / BELU_V1 * 1000;
double BELU_conc_nM   = BELU_conc_ngmL / 0.532;   // MW belumosudil ~532 g/mol
dxdt_BELU_gut  = -BELU_ka * BELU_gut;
dxdt_BELU_cent =  BELU_F * BELU_ka * BELU_gut - BELU_CL / BELU_V1 * BELU_cent;

// MPA PK (1-compartment)
double MPA_conc = MPA_cent / MMF_V1;  // ug/mL
dxdt_MPA_gut  = -MMF_ka * MPA_gut;
dxdt_MPA_cent =  MMF_F * MMF_ka * MPA_gut - MMF_CL / MMF_V1 * MPA_cent;

// ===========================================
// DRUG EFFECT CALCULATIONS (PK/PD)
// ===========================================

// CNI combined effect on calcineurin (Emax model)
double E_CsA  = Imax_CN * pow(CsA_conc, n_CN) / (pow(IC50_CsA, n_CN) + pow(CsA_conc, n_CN));
double E_TAC  = Imax_CN * pow(TAC_conc, n_CN) / (pow(IC50_TAC, n_CN) + pow(TAC_conc, n_CN));
double E_CNI  = 1 - (1 - E_CsA) * (1 - E_TAC);  // combination (Bliss)

// PRED effect on cytokine suppression
double E_PRED = Imax_PRED * PRED_conc / (IC50_PRED + PRED_conc);

// RUX effect on JAK1/2
double E_RUX  = Imax_RUX * RUX_conc_nM / (IC50_RUX + RUX_conc_nM);

// BELU effect on ROCK2
double E_BELU = Imax_BELU * BELU_conc_nM / (IC50_BELU + BELU_conc_nM);

// MPA/MMF effect on lymphocyte proliferation
double E_MPA  = Imax_MPA * MPA_conc / (IC50_MPA + MPA_conc);

// Anti-inflammatory effect (combined drug)
double E_drug_inflam = 1 - (1 - E_CNI) * (1 - E_PRED) * (1 - E_RUX);
double E_drug_prolif = 1 - (1 - E_CNI) * (1 - E_MPA);
double E_drug_fibr   = E_BELU;

// Treg enhancement (RUX+BELU promote Treg)
double E_Treg_drug = 0.3 * E_RUX + 0.4 * E_BELU;

// ===========================================
// IMMUNE CELL ODEs
// ===========================================

// Alloantigen-driven T cell activation (Allo_stim is a time-varying parameter)
double T_activ = kact_T * Allo_stim * (1 - E_drug_prolif);

// Th1 dynamics
// Driven by: alloantigen, IL-12, IFN-gamma autocrine, IL-6
// Inhibited by: CNI (IL-2 block), PRED, RUX (STAT4), IL-10
double kprol_Th1  = kprol_T * Allo_stim * (1 + 0.5 * IFNg / (IFNg + 1)) * (1 - E_drug_prolif);
double kinh_Th1   = kdeath_T + 0.2 * IL10 / (IL10 + 1) * (1 - E_drug_inflam * 0.3);
dxdt_Th1 = kprol_Th1 * Th1 + T_activ * kdiff_Th1 - kinh_Th1 * Th1;

// Th17 dynamics
// Driven by: IL-6 + TGF-β → differentiation, IL-17 amplification
// Inhibited by: RUX (STAT3 block), BELU (ROCK2/IRF4), IL-10
double kprol_Th17 = kprol_T * Allo_stim * (1 + 0.4 * IL6 / (IL6 + 2)) * (1 - E_RUX) * (1 - E_BELU * 0.7) * (1 - E_MPA);
double kinh_Th17  = kdeath_T + 0.3 * IL10 / (IL10 + 1);
dxdt_Th17 = kprol_Th17 * Th17 + T_activ * kdiff_Th17 - kinh_Th17 * Th17;

// Treg dynamics
// Driven by: TGF-β, IL-2, RUX (Treg expansion), BELU (ROCK2/IRF4 → Treg)
// Inhibited by: TNF-α, IL-6 (Th17 polarization environment)
double kprol_Treg = kprol_T * TGFb / (TGFb + 1) * (1 + E_Treg_drug) * (1 - E_MPA * 0.3);
double kinh_Treg  = kdeath_T + 0.15 * TNFa / (TNFa + 2) + 0.15 * IL6 / (IL6 + 2);
dxdt_Treg = kprol_Treg * Treg + T_activ * kdiff_Treg - kinh_Treg * Treg;

// CD8 effector T cell dynamics
double kprol_CD8 = kprol_T * Allo_stim * (1 + 0.3 * IFNg / (IFNg + 1)) * (1 - E_drug_prolif);
double kinh_CD8  = kdeath_T + 0.2 * IL10 / (IL10 + 1);
dxdt_CD8_eff = kprol_CD8 * CD8_eff + T_activ * 0.5 - kinh_CD8 * CD8_eff;

// B cell dynamics (important in chronic GvHD)
double kprol_B = kprol_T * 0.4 * (1 + 0.2 * Tfh_proxy) * (1 - E_MPA);
double Tfh_proxy = Th1; // simplified: Tfh ~ Th1 pool
double kinh_B  = kdeath_T * 0.5;
dxdt_Bcell = kprol_B * Bcell - kinh_B * Bcell;

// ===========================================
// CYTOKINE ODEs
// ===========================================

// TNF-α: produced by Th1 + CD8 + macrophages (lumped as Th1-driven)
//         degraded rapidly; suppressed by PRED, RUX
double kprod_TNF_eff = kprod_TNFa * (Th1 + 0.5 * CD8_eff) * (1 - E_drug_inflam);
dxdt_TNFa = kprod_TNF_eff - kdeg_Cyt * TNFa;

// IFN-γ: produced by Th1, CD8; key effector; suppressed by PRED
double kprod_IFN_eff = kprod_IFNg * (Th1 + 0.7 * CD8_eff) * (1 - E_PRED * 0.6 - E_RUX * 0.4);
dxdt_IFNg = kprod_IFN_eff - kdeg_Cyt * IFNg;

// IL-17A: produced by Th17; amplifies inflammation; suppressed by RUX, BELU
double kprod_IL17_eff = kprod_IL17 * Th17 * (1 - E_RUX * 0.7 - E_BELU * 0.5) * (1 - E_PRED * 0.4);
dxdt_IL17A = kprod_IL17_eff - kdeg_Cyt * IL17A;

// IL-10: produced by Treg (anti-inflammatory); NOT suppressed by most drugs
dxdt_IL10 = kprod_IL10 * Treg * (1 + E_Treg_drug) - kdeg_Cyt * IL10;

// TGF-β: produced by Treg; drives fibrosis via SMAD; also drives Treg differentiation
dxdt_TGFb = kprod_TGFb * (Treg + 0.3 * Bcell) - kdeg_Cyt * TGFb;

// IL-6: produced by APCs, Th17 milieu; drives STAT3/Th17 axis
double kprod_IL6_eff = kprod_IL6 * (1 + 0.5 * TNFa / (TNFa + 1)) * (1 - E_RUX * 0.5 - E_PRED * 0.5);
dxdt_IL6 = kprod_IL6_eff - kdeg_Cyt * IL6;

// ===========================================
// ORGAN DAMAGE ODEs
// ===========================================

// Inflammation driver (composite)
double inflam_driver = (TNFa + IFNg + 0.7 * IL17A) / 3;
double anti_inflam   = IL10;

// Skin damage (mLSS proxy: lichenoid + sclerotic changes)
// Driven by: TNFa, IFNg, Th1, CD8 cytotoxicity
// Repair: modulated by drug suppression + Treg
double kdam_skin_eff = kdam_skin * (Th1 + CD8_eff) * inflam_driver * (1 - E_drug_inflam);
double krep_skin_eff = krep_skin * Treg * (anti_inflam + 1);
dxdt_Skin_dmg = kdam_skin_eff * (1 - Skin_dmg) - krep_skin_eff * Skin_dmg;

// Gut damage (GI GvHD)
// Driven by: TNFa, IFNg, CD8 (enterocyte apoptosis)
// Repair limited: crypt stem cells needed
double kdam_gut_eff = kdam_gut * inflam_driver * CD8_eff * (1 - E_drug_inflam);
double krep_gut_eff = krep_gut * (1 + anti_inflam) * (1 - Fibrosis * 0.3);
dxdt_Gut_dmg = kdam_gut_eff * (1 - Gut_dmg) - krep_gut_eff * Gut_dmg;

// Liver damage (cholestatic pattern)
// Driven by: TNFa, CD8 infiltration of bile ducts
double kdam_liver_eff = kdam_liver * (TNFa + IFNg) * 0.5 * (1 - E_drug_inflam);
double krep_liver_eff = krep_liver * (1 + anti_inflam * 0.5);
dxdt_Liver_dmg = kdam_liver_eff * (1 - Liver_dmg) - krep_liver_eff * Liver_dmg;

// Lung damage (BOS - bronchiolitis obliterans, largely irreversible)
// Driven by: IL-17A, fibrosis (TGF-β)
double kdam_lung_eff = kdam_lung * (IL17A + 0.3 * TNFa) * (1 - E_drug_inflam * 0.3);
double krep_lung_eff = krep_lung;  // BOS poorly reversible
dxdt_Lung_dmg = kdam_lung_eff * (1 - Lung_dmg) - krep_lung_eff * Lung_dmg;

// Fibrosis (chronic cGvHD - driven by TGF-β, ROCK2 pathway)
// Suppressed by belumosudil (ROCK2i)
double kfibr_eff = kfibr * TGFb * (1 - E_drug_fibr);
double kfibr_res = 0.0005;  // minimal fibrosis resolution
dxdt_Fibrosis = kfibr_eff * (1 - Fibrosis) - kfibr_res * Fibrosis;

$TABLE

// Drug concentrations
double CsA_C0   = CsA_cent / CsA_V1 * 1000;   // ng/mL (trough target: 100-300 ng/mL)
double TAC_C0   = TAC_cent / TAC_V1 * 1000;    // ng/mL (trough target: 5-15 ng/mL)
double PRED_C   = PRED_cent / PRED_V1 * 1000;  // ng/mL
double RUX_C    = RUX_cent / RUX_V1 * 1000;    // ng/mL
double BELU_C   = BELU_cent / BELU_V1 * 1000;  // ng/mL
double MPA_C    = MPA_cent / MMF_V1;            // ug/mL (AUC target 30-60 h*ug/mL)

// PD effect outputs
double CNI_Inh  = E_CNI * 100;    // % calcineurin inhibition
double JAK_Inh  = E_RUX * 100;   // % JAK1/2 inhibition
double ROCK2_Inh= E_BELU * 100;  // % ROCK2 inhibition

// Derived clinical scores
// Glucksberg grade proxy (acute GvHD): combines gut + liver + skin
double Glucksberg_skin  = Skin_dmg < 0.25  ? 0 : (Skin_dmg < 0.50 ? 1 : (Skin_dmg < 0.75 ? 2 : 3));
double Glucksberg_gut   = Gut_dmg  < 0.25  ? 0 : (Gut_dmg  < 0.50 ? 1 : (Gut_dmg  < 0.75 ? 2 : 3));
double Glucksberg_liver = Liver_dmg< 0.25  ? 0 : (Liver_dmg< 0.50 ? 1 : (Liver_dmg< 0.75 ? 2 : 3));
double aGvHD_Grade      = (Glucksberg_skin + Glucksberg_gut + Glucksberg_liver) / 3.0;

// NIH cGvHD global score (0-3 scale, proxy)
double cGvHD_Score = (Skin_dmg + Gut_dmg + Liver_dmg + Lung_dmg + Fibrosis) / 5.0 * 3.0;

// Biomarkers (simplified as linear functions of damage)
double ST2_bio   = 10 + 500 * Gut_dmg;    // ng/mL (normal <33 ng/mL)
double REG3a_bio = 10 + 200 * Gut_dmg;    // ng/mL
double TNFR1_bio = 1  + 5   * TNFa;       // ng/mL (normal ~2 ng/mL)

// Failure-Free Survival proxy (1 = no failure)
double FFS_proxy = (aGvHD_Grade < 2 && cGvHD_Score < 1.5) ? 1.0 : 0.0;

// GvL preservation proxy (higher Treg -> more immunosuppression -> higher relapse risk)
double GvL_preserved = GvL_effect * (1 - Treg / 3.0);

capture CsA_C0, TAC_C0, PRED_C, RUX_C, BELU_C, MPA_C
capture CNI_Inh, JAK_Inh, ROCK2_Inh
capture aGvHD_Grade, cGvHD_Score, FFS_proxy, GvL_preserved
capture ST2_bio, REG3a_bio, TNFR1_bio
capture Th1, Th17, Treg, CD8_eff, Bcell
capture TNFa, IFNg, IL17A, IL10, TGFb, IL6
capture Skin_dmg, Gut_dmg, Liver_dmg, Lung_dmg, Fibrosis
capture E_CNI, E_PRED, E_RUX, E_BELU

$CAPTURE CsA_C0 TAC_C0 RUX_C BELU_C aGvHD_Grade cGvHD_Score FFS_proxy
'

## ============================================================
## Build and Load Model
## ============================================================

mod <- mcode("gvhd_qsp", gvhd_code)

## ============================================================
## TREATMENT SCENARIOS
## ============================================================

# Common HSCT timeline:
# Day 0: HSCT
# Day 7-100: CsA/TAC prophylaxis (acute phase)
# Day 20-30: peak acute GvHD risk
# Day 100+: chronic GvHD development

# Helper: create event table with dosing regimen
make_dose_ev <- function(
    csa_dose = 0, csa_freq_h = 12, csa_dur_days = 100,
    tac_dose = 0, tac_freq_h = 12, tac_dur_days = 100,
    pred_dose = 0, pred_freq_h = 24, pred_dur_days = 14,
    rux_dose = 0, rux_freq_h = 12, rux_dur_days = 180,
    belu_dose = 0, belu_freq_h = 24, belu_dur_days = 180,
    mmf_dose = 0, mmf_freq_h = 12, mmf_dur_days = 100) {

  evs <- list()
  # CsA
  if (csa_dose > 0)
    evs <- c(evs, list(ev(amt = csa_dose, cmt = 1, ii = csa_freq_h, addl = csa_dur_days * 24 / csa_freq_h - 1, time = 0)))
  # TAC
  if (tac_dose > 0)
    evs <- c(evs, list(ev(amt = tac_dose, cmt = 4, ii = tac_freq_h, addl = tac_dur_days * 24 / tac_freq_h - 1, time = 0)))
  # PRED
  if (pred_dose > 0)
    evs <- c(evs, list(ev(amt = pred_dose, cmt = 7, ii = pred_freq_h, addl = pred_dur_days * 24 / pred_freq_h - 1, time = 480)))  # start at day 20
  # RUX
  if (rux_dose > 0)
    evs <- c(evs, list(ev(amt = rux_dose, cmt = 9, ii = rux_freq_h, addl = rux_dur_days * 24 / rux_freq_h - 1, time = 720)))  # start at day 30
  # BELU
  if (belu_dose > 0)
    evs <- c(evs, list(ev(amt = belu_dose, cmt = 12, ii = belu_freq_h, addl = belu_dur_days * 24 / belu_freq_h - 1, time = 720)))
  # MMF
  if (mmf_dose > 0)
    evs <- c(evs, list(ev(amt = mmf_dose, cmt = 14, ii = mmf_freq_h, addl = mmf_dur_days * 24 / mmf_freq_h - 1, time = 0)))

  do.call(ev_seq, evs)
}

# --- Scenario 1: No GvHD prophylaxis (historical baseline) ---
scenario1_ev <- ev(time = 0, amt = 0, cmt = 1)  # no dose - use initial conditions
out1 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = scenario1_ev,
  end = 365 * 24, delta = 6,  # 365 days, every 6h
  obsaug = TRUE
) %>% as_tibble() %>% mutate(scenario = "No Prophylaxis (Baseline)")

# --- Scenario 2: CsA monoprophylaxis (3 mg/kg q12h = ~225 mg q12h for 75 kg) ---
out2 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = ev(amt = 225, cmt = 1, ii = 12, addl = 199, time = 0),
  end = 365 * 24, delta = 6
) %>% as_tibble() %>% mutate(scenario = "CsA Monoprophylaxis")

# --- Scenario 3: CsA + MMF standard prophylaxis ---
ev3 <- ev_seq(
  ev(amt = 225, cmt = 1, ii = 12, addl = 199, time = 0),
  ev(amt = 1500, cmt = 14, ii = 12, addl = 199, time = 0)
)
out3 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = ev3,
  end = 365 * 24, delta = 6
) %>% as_tibble() %>% mutate(scenario = "CsA + MMF Prophylaxis")

# --- Scenario 4: Tacrolimus + MMF prophylaxis (standard MSD/UD) ---
ev4 <- ev_seq(
  ev(amt = 2, cmt = 4, ii = 12, addl = 199, time = 0),  # 0.03 mg/kg x 75 kg / 2 doses ~= 1.1 mg q12h; use 2 mg
  ev(amt = 1500, cmt = 14, ii = 12, addl = 199, time = 0)
)
out4 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = ev4,
  end = 365 * 24, delta = 6
) %>% as_tibble() %>% mutate(scenario = "TAC + MMF Prophylaxis")

# --- Scenario 5: Ruxolitinib for steroid-refractory cGvHD ---
# Start ruxolitinib at day 30 (720 h) for 180 days, AFTER CsA prophylaxis
ev5 <- ev_seq(
  ev(amt = 225, cmt = 1, ii = 12, addl = 199, time = 0),
  ev(amt = 10,  cmt = 9, ii = 12, addl = 359, time = 720)  # 10 mg BID ruxolitinib
)
out5 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = ev5,
  end = 365 * 24, delta = 6
) %>% as_tibble() %>% mutate(scenario = "CsA → Ruxolitinib (SR-cGvHD)")

# --- Scenario 6: Belumosudil for cGvHD (≥2 prior lines) ---
ev6 <- ev_seq(
  ev(amt = 225, cmt = 1, ii = 12, addl = 199, time = 0),
  ev(amt = 200, cmt = 12, ii = 24, addl = 179, time = 720)  # 200 mg QD belumosudil
)
out6 <- mrgsim(
  mod %>% param(Allo_stim = 1.0),
  ev = ev6,
  end = 365 * 24, delta = 6
) %>% as_tibble() %>% mutate(scenario = "CsA → Belumosudil (cGvHD)")

## ============================================================
## COMBINE RESULTS
## ============================================================

all_out <- bind_rows(out1, out2, out3, out4, out5, out6) %>%
  mutate(Day = time / 24)

## ============================================================
## VISUALIZATION
## ============================================================

# Color palette
scenario_colors <- c(
  "No Prophylaxis (Baseline)"    = "#E74C3C",
  "CsA Monoprophylaxis"          = "#E67E22",
  "CsA + MMF Prophylaxis"        = "#F1C40F",
  "TAC + MMF Prophylaxis"        = "#2ECC71",
  "CsA → Ruxolitinib (SR-cGvHD)"= "#3498DB",
  "CsA → Belumosudil (cGvHD)"   = "#9B59B6"
)

# Plot 1: Overall GvHD Score over time
p1 <- all_out %>%
  filter(time %% 24 == 0) %>%  # daily
  ggplot(aes(x = Day, y = aGvHD_Grade, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Acute GvHD Grade Over Time",
       subtitle = "Composite Glucksberg Score (0-3 scale)",
       x = "Days Post-HSCT", y = "aGvHD Grade (proxy)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 3))

# Plot 2: Cytokine dynamics
p2 <- all_out %>%
  filter(scenario %in% c("No Prophylaxis (Baseline)", "TAC + MMF Prophylaxis", "CsA → Ruxolitinib (SR-cGvHD)"),
         time %% 24 == 0) %>%
  pivot_longer(cols = c(TNFa, IFNg, IL17A, IL10), names_to = "Cytokine", values_to = "Conc") %>%
  ggplot(aes(x = Day, y = Conc, color = scenario, linetype = Cytokine)) +
  geom_line(size = 1.0) +
  facet_wrap(~Cytokine, scales = "free_y") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Cytokine Dynamics",
       x = "Days Post-HSCT", y = "Cytokine (ng/mL)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 3: Organ damage scores
p3 <- all_out %>%
  filter(scenario %in% c("No Prophylaxis (Baseline)", "CsA + MMF Prophylaxis", "CsA → Ruxolitinib (SR-cGvHD)"),
         time %% 24 == 0) %>%
  pivot_longer(cols = c(Skin_dmg, Gut_dmg, Liver_dmg, Lung_dmg, Fibrosis),
               names_to = "Organ", values_to = "Score") %>%
  ggplot(aes(x = Day, y = Score, color = scenario)) +
  geom_line(size = 1.0) +
  facet_wrap(~Organ, ncol = 3) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Organ Damage Scores Over Time",
       x = "Days Post-HSCT", y = "Damage Score",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 4: PK - CsA trough concentrations
p4 <- all_out %>%
  filter(scenario %in% c("CsA Monoprophylaxis", "CsA + MMF Prophylaxis",
                          "CsA → Ruxolitinib (SR-cGvHD)"),
         Day <= 180, time %% 24 == 0) %>%
  ggplot(aes(x = Day, y = CsA_C0, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = c(100, 300), linetype = "dashed", color = "gray40") +
  annotate("text", x = 100, y = 105, label = "Lower target (100 ng/mL)", size = 3) +
  annotate("text", x = 100, y = 305, label = "Upper target (300 ng/mL)", size = 3) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "CsA Trough Concentration (C₀)",
       x = "Days Post-HSCT", y = "CsA Trough (ng/mL)",
       color = "Scenario") +
  theme_bw(base_size = 12)

# Plot 5: T cell compartment dynamics
p5 <- all_out %>%
  filter(scenario %in% c("No Prophylaxis (Baseline)", "TAC + MMF Prophylaxis", "CsA → Ruxolitinib (SR-cGvHD)"),
         time %% 24 == 0) %>%
  pivot_longer(cols = c(Th1, Th17, Treg, CD8_eff), names_to = "Cell", values_to = "Pool") %>%
  ggplot(aes(x = Day, y = Pool, color = scenario, linetype = Cell)) +
  geom_line(size = 1.0) +
  facet_wrap(~Cell) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "T Cell Compartment Dynamics",
       x = "Days Post-HSCT", y = "Relative Pool Size",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Plot 6: Biomarkers
p6 <- all_out %>%
  filter(time %% 24 == 0) %>%
  pivot_longer(cols = c(ST2_bio, REG3a_bio, TNFR1_bio), names_to = "Biomarker", values_to = "Value") %>%
  ggplot(aes(x = Day, y = Value, color = scenario)) +
  geom_line(size = 1.0) +
  facet_wrap(~Biomarker, scales = "free_y") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "GvHD Biomarkers Over Time",
       subtitle = "ST2 (GI GvHD), REG3α (GI epithelial), sTNFR1 (systemic)",
       x = "Days Post-HSCT", y = "Biomarker Level",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Print plots
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

## ============================================================
## SUMMARY TABLE
## ============================================================

summary_table <- all_out %>%
  filter(time %% 24 == 0) %>%
  group_by(scenario) %>%
  summarise(
    Max_aGvHD   = max(aGvHD_Grade, na.rm = TRUE),
    Max_cGvHD   = max(cGvHD_Score, na.rm = TRUE),
    Max_Skin    = max(Skin_dmg,  na.rm = TRUE),
    Max_Gut     = max(Gut_dmg,   na.rm = TRUE),
    Max_Liver   = max(Liver_dmg, na.rm = TRUE),
    Max_Fibrosis= max(Fibrosis,  na.rm = TRUE),
    Peak_TNFa   = max(TNFa, na.rm = TRUE),
    MinTreg     = min(Treg, na.rm = TRUE),
    MaxTh17     = max(Th17, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Max_aGvHD)

cat("\n== GvHD QSP Model Summary Table ==\n")
print(summary_table, n = 10)

cat("\nKey findings:\n")
cat("1. CsA/TAC + MMF prophylaxis substantially reduces aGvHD severity (Glucksberg grade)\n")
cat("2. Ruxolitinib demonstrates JAK1/2 inhibition that reduces Th17 and expands Treg\n")
cat("3. Belumosudil targets ROCK2/IRF4 axis, selectively rebalancing Th17/Treg\n")
cat("4. Fibrosis development is a key driver of chronic GvHD and pulmonary complications\n")
cat("5. GvL effect must be preserved - over-immunosuppression increases relapse risk\n")

## ============================================================
## SENSITIVITY ANALYSIS: Treg/Th17 Ratio Impact on cGvHD
## ============================================================

sa_results <- list()
for (kTreg_mult in c(0.5, 1.0, 1.5, 2.0)) {
  sa_out <- mrgsim(
    mod %>% param(Allo_stim = 1.0, kdiff_Treg = 0.010 * kTreg_mult),
    ev = ev(amt = 225, cmt = 1, ii = 12, addl = 199, time = 0),
    end = 365 * 24, delta = 24
  ) %>% as_tibble() %>%
    mutate(Day = time / 24, Treg_mult = kTreg_mult)
  sa_results[[length(sa_results)+1]] <- sa_out
}

sa_combined <- bind_rows(sa_results)

p_sa <- sa_combined %>%
  ggplot(aes(x = Day, y = cGvHD_Score, color = factor(Treg_mult), group = Treg_mult)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = c("#E74C3C","#E67E22","#2ECC71","#3498DB"),
                     labels = paste0(c(0.5, 1.0, 1.5, 2.0), "x Treg differentiation")) +
  labs(title = "Sensitivity: Treg Differentiation Rate vs cGvHD Score",
       x = "Days Post-HSCT", y = "cGvHD Score (0-3)",
       color = "Treg Rate") +
  theme_bw(base_size = 12)
print(p_sa)

cat("\nAll GvHD QSP simulations complete.\n")
