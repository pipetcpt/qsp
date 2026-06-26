## ============================================================================
## Gastric Cancer (위선암) — QSP mrgsolve Model
##
## Description:
##   A multi-compartment ODE-based QSP model for gastric cancer integrating
##   drug PK (trastuzumab, ramucirumab, nivolumab, T-DXd, zolbetuximab,
##   capecitabine/FOLFOX) with disease PD (tumor growth inhibition, HER2
##   signaling, VEGF angiogenesis, immune dynamics, CEA biomarker).
##
## Clinical Trial Parameter Sources:
##   - ToGA trial (Bang et al. Lancet 2010; PMID: 20728210): Trastuzumab + FOLFOX/XELOX HER2+
##   - RAINBOW trial (Wilke et al. Lancet Oncol 2014; PMID: 25043997): Ramucirumab + paclitaxel
##   - CheckMate 649 (Janjigian et al. Nature 2021; PMID: 33811764): Nivolumab + chemo CPS≥5
##   - SPOTLIGHT (Shitara et al. Lancet 2023; PMID: 36764316): Zolbetuximab + mFOLFOX6
##   - DESTINY-Gastric01 (Shitara et al. NEJM 2020; PMID: 32469182): T-DXd 2L HER2+
##   - FLOT4 (Al-Batran et al. Lancet 2019; PMID: 30982686): FLOT perioperative
##   - KEYNOTE-811 (Janjigian et al. Nature 2024; PMID: 38811132): Pembro + Tras + chemo HER2+
##
## Compartments (18 total):
##   Drug PK (12): Trastuzumab C/P, Ramucirumab C/P, Nivolumab C/P,
##                 Capecitabine (gut/plasma/active), T-DXd plasma, Zolbetuximab
##   Disease PD (6): Tumor burden, HER2 signaling, VEGF-free,
##                   CD8 T effector, Treg, TAM_M2, CEA, Cancer Stem Cell
##
## Author: QSP Library (CCR auto-generated)
## Date: 2026-06-23
## ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================================
## MODEL CODE
## ============================================================================

gc_model_code <- '
$PROB
Gastric Cancer QSP Model — Integrated PK/PD with Immune Dynamics
Clinical Subtypes: HER2+, MSI-H, CLDN18.2+, Standard

$PARAM @annotated
// ---- Trastuzumab PK (2-compartment, IV) ----
// Source: Hayashi et al. CPT 2015; Gibiansky et al. J Pharmacokinet 2012
CL_tras   : 0.218  : Trastuzumab clearance (L/day)
V1_tras   : 2.91   : Trastuzumab central volume (L)
V2_tras   : 2.45   : Trastuzumab peripheral volume (L)
Q_tras    : 0.742  : Trastuzumab inter-compartmental CL (L/day)
F_tras    : 1.0    : Trastuzumab bioavailability (IV=1)
kon_tras  : 0.36   : Trastuzumab-HER2 on-rate (L/nmol/day)
koff_tras : 0.0024 : Trastuzumab-HER2 off-rate (1/day)
HER2_total: 2.5    : Total HER2 receptor density (nmol/L, IHC3+)

// ---- Ramucirumab PK (2-compartment, IV) ----
// Source: Tabernero et al. Ann Oncol 2015; RAINBOW PK analysis
CL_ramu   : 0.282  : Ramucirumab clearance (L/day)
V1_ramu   : 3.14   : Ramucirumab central volume (L)
V2_ramu   : 1.89   : Ramucirumab peripheral volume (L)
Q_ramu    : 0.562  : Ramucirumab inter-compartmental CL (L/day)
kon_ramu  : 0.48   : Ramucirumab-VEGFR2 on-rate (L/nmol/day)
koff_ramu : 0.0018 : Ramucirumab-VEGFR2 off-rate (1/day)
VEGFR2_tot: 1.8    : Total VEGFR2 concentration (nmol/L)

// ---- Nivolumab PK (2-compartment, IV) ----
// Source: Bajaj et al. CPT:PSP 2017; CheckMate 649 PK
CL_nivo   : 0.312  : Nivolumab clearance (L/day)
V1_nivo   : 3.75   : Nivolumab central volume (L)
V2_nivo   : 2.94   : Nivolumab peripheral volume (L)
Q_nivo    : 0.724  : Nivolumab inter-compartmental CL (L/day)
kon_nivo  : 0.55   : Nivolumab-PD1 on-rate (L/nmol/day)
koff_nivo : 0.0035 : Nivolumab-PD1 off-rate (1/day)
PD1_tot   : 3.2    : Total PD-1 on T cells (nmol/L)

// ---- Capecitabine/5-FU PK ----
// Source: Miwa et al. Eur J Cancer 1998; REAL2 trial PK
ka_cape   : 1.92   : Capecitabine absorption rate (1/day)
F_cape    : 0.70   : Capecitabine bioavailability (oral)
CL_cape   : 38.4   : Capecitabine clearance (L/day)
V_cape    : 25.2   : Capecitabine volume (L)
k_conv    : 2.16   : Conversion to 5-FU (1/day, 3-step via CES, CDA, TP)
CL_FU     : 86.4   : 5-FU (F-dUMP) clearance (L/day)
V_FU      : 18.0   : 5-FU volume (L)

// ---- T-DXd PK (ADC) ----
// Source: Ogitani et al. Clin Cancer Res 2016; DESTINY-Gastric01 PK
CL_tdxd   : 0.198  : T-DXd clearance (L/day)
V1_tdxd   : 3.28   : T-DXd central volume (L)
kon_tdxd  : 0.42   : T-DXd-HER2 on-rate (L/nmol/day)
koff_tdxd : 0.0031 : T-DXd-HER2 off-rate (1/day)
k_payload : 0.144  : DXd payload release rate (1/day)
CL_dxd    : 46.8   : Free DXd clearance (L/day)
V_dxd     : 22.5   : Free DXd volume (L)

// ---- Zolbetuximab PK ----
// Source: Türeci et al. Ann Oncol 2019; SPOTLIGHT trial PK
CL_zolbe  : 0.245  : Zolbetuximab clearance (L/day)
V1_zolbe  : 3.05   : Zolbetuximab central volume (L)
kon_zolbe : 0.39   : Zolbetuximab-CLDN18.2 on-rate (L/nmol/day)
koff_zolbe: 0.0028 : Zolbetuximab-CLDN18.2 off-rate (1/day)
CLDN_tot  : 4.2    : Total CLDN18.2 (nmol/L, CLDN18.2+ tumor)

// ---- Tumor Growth Inhibition (Simeoni + Logistic) ----
// Source: Simeoni et al. Cancer Res 2004; adapted for GC
lambda0   : 0.0215 : Tumor exponential growth rate (1/day)
lambda1   : 0.385  : Tumor linear growth rate (mm3/day)
k1_tgi    : 0.034  : Drug-induced tumor cell death rate (1/day/conc)
k2_tgi    : 0.018  : Damage-to-death transit rate (1/day)
kDE       : 0.0048 : Natural tumor cell death rate (1/day)
TV0       : 500    : Initial tumor volume (mm3)
TV_max    : 15000  : Maximum tumor volume (mm3)

// ---- HER2 Signaling ----
ksynHER2  : 0.024  : HER2 signaling synthesis rate (1/day)
kdegHER2  : 0.096  : HER2 signaling degradation rate (1/day)
EC50_tras : 0.85   : Trastuzumab EC50 for HER2 inhibition (nmol/L)
n_tras    : 1.4    : Trastuzumab Hill coefficient

// ---- VEGF / Angiogenesis ----
ksynVEGF  : 0.185  : VEGF production rate (ng/mL/day)
kdegVEGF  : 0.432  : VEGF degradation rate (1/day)
EC50_ramu : 1.12   : Ramucirumab EC50 for VEGF-angiogenesis inh. (nmol/L)

// ---- Immune Dynamics ----
// CD8+ T effector cells
ksynCD8   : 0.048  : CD8 T-cell source rate (cells/μL/day)
kprolCD8  : 0.24   : CD8 T-cell proliferation rate (1/day)
kdegCD8   : 0.072  : CD8 T-cell death rate (1/day)
kexhCD8   : 0.018  : CD8 T-cell exhaustion rate (1/day, PD1-mediated)
EC50_nivo_CD8: 0.65: Nivolumab EC50 for CD8 restoration (nmol/L)

// Treg dynamics
ksynTreg  : 0.024  : Treg source rate (cells/μL/day)
kprolTreg : 0.12   : Treg proliferation rate (1/day)
kdegTreg  : 0.048  : Treg death rate (1/day)
kTregSup  : 0.085  : Treg suppression of CD8 rate (1/day per Treg cell)

// TAM M2 dynamics
ksynTAM   : 0.036  : TAM M2 source (cells/μL/day)
kdegTAM   : 0.072  : TAM M2 death rate (1/day)
kTAMVEGF  : 0.045  : TAM→VEGF production coupling (1/day)

// ---- CEA Biomarker ----
ksynCEA   : 0.0054 : CEA synthesis rate (ng/mL/day, tumor-derived)
kdegCEA   : 0.0288 : CEA clearance rate (1/day)
CEA0      : 25.0   : Baseline CEA (ng/mL)
corr_TV_CEA: 0.015 : Tumor volume to CEA coupling (per mm3)

// ---- Cancer Stem Cell ----
kCSC_self  : 0.0048: CSC self-renewal rate (1/day)
kCSC_diff  : 0.018 : CSC differentiation to progenitor rate (1/day)
kdegCSC    : 0.0072: CSC death rate (1/day)
CSC0       : 50    : Initial CSC count (cells/μL)

// ---- Killing Efficacy Parameters ----
Emax_tras  : 0.72  : Max tumor kill by trastuzumab
Emax_ramu  : 0.55  : Max tumor kill by ramucirumab (antiangiogenic)
Emax_nivo  : 0.65  : Max CD8 restoration by nivolumab (immune-mediated)
Emax_chemo : 0.88  : Max tumor kill by chemotherapy (5-FU/oxali/taxane)
Emax_tdxd  : 0.85  : Max tumor kill by T-DXd
Emax_zolbe : 0.68  : Max tumor kill by zolbetuximab (ADCC)
EC50_FU    : 0.45  : 5-FU EC50 for tumor kill (μmol/L)
EC50_dxd   : 0.28  : DXd EC50 for tumor kill (μmol/L)
EC50_zolbe_kill: 0.92: Zolbetuximab EC50 for ADCC killing (nmol/L)

$CMT @annotated
// ---- Drug PK compartments ----
TRAS_C    : Trastuzumab central (nmol/L)
TRAS_P    : Trastuzumab peripheral (nmol/L)
TRAS_B    : Trastuzumab-HER2 bound complex (nmol/L)
RAMU_C    : Ramucirumab central (nmol/L)
RAMU_P    : Ramucirumab peripheral (nmol/L)
RAMU_B    : Ramucirumab-VEGFR2 bound (nmol/L)
NIVO_C    : Nivolumab central (nmol/L)
NIVO_P    : Nivolumab peripheral (nmol/L)
NIVO_B    : Nivolumab-PD1 bound (nmol/L)
CAPE_GUT  : Capecitabine gut (mg)
CAPE_C    : Capecitabine central plasma (mg/L)
FU5_C     : 5-FU active (F-dUMP) plasma (μmol/L)
TDXD_C    : T-DXd central (nmol/L)
TDXD_B    : T-DXd-HER2 bound (nmol/L)
DXD_FREE  : Free DXd payload (μmol/L)
ZOLBE_C   : Zolbetuximab central (nmol/L)
ZOLBE_B   : Zolbetuximab-CLDN18.2 bound (nmol/L)
PEMBRO_C  : Pembrolizumab central (nmol/L)

// ---- Disease PD compartments ----
TV        : Tumor volume (mm3)
TV_DMG1   : Damaged tumor cells stage 1 (mm3)
TV_DMG2   : Damaged tumor cells stage 2 (mm3)
HER2_SIG  : HER2 downstream signaling activity (AU)
VEGF_FREE : Free VEGF-A (ng/mL)
CD8_TEFF  : CD8+ effector T cells (cells/μL)
TREG      : Regulatory T cells (cells/μL)
TAM_M2    : M2 tumor-associated macrophages (cells/μL)
CEA_BM    : CEA biomarker (ng/mL)
CSC       : Cancer stem cells (cells/μL)

$MAIN
// ---- Initial Conditions ----
TV_0      = TV0;
TV_DMG1_0 = 0;
TV_DMG2_0 = 0;
HER2_SIG_0= 1.0;     // normalized to 1 (baseline active)
VEGF_FREE_0 = ksynVEGF / kdegVEGF;  // steady state
CD8_TEFF_0 = ksynCD8 / kdegCD8;     // steady state ~0.67 cells/μL (relative)
TREG_0    = ksynTreg / kdegTreg;     // steady state Treg
TAM_M2_0  = ksynTAM / kdegTAM;      // steady state TAM M2
CEA_BM_0  = CEA0;
CSC_0     = CSC0;

$ODE
// ================================================================
// DRUG PK — Trastuzumab (2-cmt + target binding)
// ================================================================
double HER2_free = HER2_total - TRAS_B - TDXD_B;  // free HER2
if(HER2_free < 0) HER2_free = 0;

double R_tras_on  = kon_tras * TRAS_C * HER2_free;
double R_tras_off = koff_tras * TRAS_B;
double CL_tras_tot = CL_tras + koff_tras;  // net clearance with TMDD

dxdt_TRAS_C  = -CL_tras/V1_tras * TRAS_C
               - Q_tras/V1_tras * TRAS_C
               + Q_tras/V2_tras * TRAS_P
               - R_tras_on + R_tras_off;
dxdt_TRAS_P  =  Q_tras/V1_tras * TRAS_C - Q_tras/V2_tras * TRAS_P;
dxdt_TRAS_B  =  R_tras_on - R_tras_off - koff_tras * TRAS_B;

// ================================================================
// DRUG PK — Ramucirumab (2-cmt + VEGFR2 binding)
// ================================================================
double VEGFR2_free = VEGFR2_tot - RAMU_B;
if(VEGFR2_free < 0) VEGFR2_free = 0;

double R_ramu_on  = kon_ramu * RAMU_C * VEGFR2_free;
double R_ramu_off = koff_ramu * RAMU_B;

dxdt_RAMU_C  = -CL_ramu/V1_ramu * RAMU_C
               - Q_ramu/V1_ramu * RAMU_C
               + Q_ramu/V2_ramu * RAMU_P
               - R_ramu_on + R_ramu_off;
dxdt_RAMU_P  =  Q_ramu/V1_ramu * RAMU_C - Q_ramu/V2_ramu * RAMU_P;
dxdt_RAMU_B  =  R_ramu_on - R_ramu_off - koff_ramu * RAMU_B;

// ================================================================
// DRUG PK — Nivolumab (2-cmt + PD-1 binding)
// ================================================================
double PD1_free = PD1_tot - NIVO_B;
if(PD1_free < 0) PD1_free = 0;

double R_nivo_on  = kon_nivo * NIVO_C * PD1_free;
double R_nivo_off = koff_nivo * NIVO_B;

dxdt_NIVO_C  = -CL_nivo/V1_nivo * NIVO_C
               - Q_nivo/V1_nivo * NIVO_C
               + Q_nivo/V2_nivo * NIVO_P
               - R_nivo_on + R_nivo_off;
dxdt_NIVO_P  =  Q_nivo/V1_nivo * NIVO_C - Q_nivo/V2_nivo * NIVO_P;
dxdt_NIVO_B  =  R_nivo_on - R_nivo_off - koff_nivo * NIVO_B;

// ================================================================
// DRUG PK — Capecitabine / 5-FU oral
// ================================================================
dxdt_CAPE_GUT = -ka_cape * CAPE_GUT;
dxdt_CAPE_C   =  F_cape * ka_cape * CAPE_GUT - CL_cape/V_cape * CAPE_C;
dxdt_FU5_C    =  k_conv * CAPE_C - CL_FU/V_FU * FU5_C;

// ================================================================
// DRUG PK — T-DXd (ADC with payload release)
// ================================================================
double HER2_free_tdxd = HER2_total - TRAS_B - TDXD_B;
if(HER2_free_tdxd < 0) HER2_free_tdxd = 0;

double R_tdxd_on  = kon_tdxd * TDXD_C * HER2_free_tdxd;
double R_tdxd_off = koff_tdxd * TDXD_B;

dxdt_TDXD_C  = -CL_tdxd/V1_tdxd * TDXD_C - R_tdxd_on + R_tdxd_off;
dxdt_TDXD_B  =  R_tdxd_on - R_tdxd_off - k_payload * TDXD_B;
dxdt_DXD_FREE=  k_payload * TDXD_B - CL_dxd/V_dxd * DXD_FREE;

// ================================================================
// DRUG PK — Zolbetuximab (CLDN18.2 binding)
// ================================================================
double CLDN_free = CLDN_tot - ZOLBE_B;
if(CLDN_free < 0) CLDN_free = 0;

double R_zolbe_on  = kon_zolbe * ZOLBE_C * CLDN_free;
double R_zolbe_off = koff_zolbe * ZOLBE_B;

dxdt_ZOLBE_C = -CL_zolbe/V1_zolbe * ZOLBE_C - R_zolbe_on + R_zolbe_off;
dxdt_ZOLBE_B =  R_zolbe_on - R_zolbe_off - koff_zolbe * ZOLBE_B;

// ================================================================
// DRUG PK — Pembrolizumab (simplified 1-cmt)
// ================================================================
dxdt_PEMBRO_C = -CL_nivo/V1_nivo * PEMBRO_C;  // similar PK to nivolumab

// ================================================================
// DISEASE PD — Tumor Growth Inhibition (Simeoni transit model)
// ================================================================
// Drug effect functions (Emax model)
double E_tras  = Emax_tras  * pow(TRAS_B, n_tras)   / (pow(EC50_tras, n_tras) + pow(TRAS_B, n_tras));
double E_ramu  = Emax_ramu  * RAMU_B   / (EC50_ramu  + RAMU_B);
double E_FU    = Emax_chemo * FU5_C    / (EC50_FU   + FU5_C);
double E_dxd   = Emax_tdxd  * DXD_FREE / (EC50_dxd  + DXD_FREE);
double E_zolbe = Emax_zolbe * ZOLBE_B  / (EC50_zolbe_kill + ZOLBE_B);

// CD8 T cell killing effect on tumor
double E_CD8   = 0.12 * CD8_TEFF / (0.5 + CD8_TEFF);  // immune-mediated kill

// Combined drug effect (additive-independent Bliss)
double E_total = 1.0 - (1.0 - E_tras) * (1.0 - E_ramu) * (1.0 - E_FU)
                       * (1.0 - E_dxd) * (1.0 - E_zolbe) * (1.0 - E_CD8);

// Simeoni tumor growth model with logistic correction
double TV_total = TV + TV_DMG1 + TV_DMG2;
double growth_rate;
if(TV_total < TV_max) {
  double logistic_factor = 1.0 - TV_total / TV_max;
  growth_rate = 2.0 * lambda0 * lambda1 / (lambda1 + 2.0 * lambda0 * TV) * TV * logistic_factor;
} else {
  growth_rate = 0;
}

double k1_eff  = k1_tgi * E_total;
dxdt_TV      =  growth_rate - k1_eff * TV - kDE * TV;
dxdt_TV_DMG1 =  k1_eff * TV  - k2_tgi * TV_DMG1;
dxdt_TV_DMG2 =  k2_tgi * TV_DMG1 - k2_tgi * TV_DMG2;

// ================================================================
// HER2 Signaling Dynamics
// ================================================================
double InhHER2 = pow(TRAS_B, n_tras) / (pow(EC50_tras, n_tras) + pow(TRAS_B, n_tras));
dxdt_HER2_SIG = ksynHER2 - kdegHER2 * HER2_SIG - kdegHER2 * InhHER2 * HER2_SIG;

// ================================================================
// VEGF-A Free Concentration
// ================================================================
// VEGF produced by tumor and TAM M2, consumed by VEGFR2 binding
double VEGF_prod = ksynVEGF * (1.0 + 0.1 * TV_total / TV0) + kTAMVEGF * TAM_M2;
double InhVEGF   = E_ramu;  // ramucirumab blocks VEGFR2
dxdt_VEGF_FREE = VEGF_prod - kdegVEGF * VEGF_FREE * (1.0 + InhVEGF);

// ================================================================
// Immune Dynamics — CD8+ Effector T cells
// ================================================================
// PD-1 occupancy by nivolumab or pembrolizumab
double PD1_occ_nivo  = NIVO_B   / PD1_tot;
double PD1_occ_pembro= PEMBRO_C / (PEMBRO_C + EC50_nivo_CD8);
double PD1_total_occ = 1.0 - (1.0 - PD1_occ_nivo) * (1.0 - PD1_occ_pembro);

// CD8 exhaustion relieved by anti-PD1
double kexh_eff  = kexhCD8 * (1.0 - Emax_nivo * PD1_total_occ);
// Treg suppression of CD8
double kTreg_sup = kTregSup * TREG;
// Tumor antigen-driven CD8 proliferation
double CD8_prolif = kprolCD8 * TV_total / (TV_total + 2000.0) * CD8_TEFF;

dxdt_CD8_TEFF = ksynCD8 + CD8_prolif - kdegCD8 * CD8_TEFF
                - kexh_eff * CD8_TEFF - kTreg_sup * CD8_TEFF;

// ================================================================
// Immune Dynamics — Regulatory T cells
// ================================================================
// Treg expanded by TGF-β (proxy: tumor volume & TAM M2)
double Treg_expand = kprolTreg * TV_total / (TV_total + 5000.0) + 0.02 * TAM_M2;
dxdt_TREG = ksynTreg + Treg_expand * TREG - kdegTreg * TREG;

// ================================================================
// TME — TAM M2 Macrophages
// ================================================================
// TAM M2 recruited by tumor-derived CCL2, expanded by IL-10
double TAM_recruit = 0.005 * TV_total;
dxdt_TAM_M2 = ksynTAM + TAM_recruit - kdegTAM * TAM_M2;

// ================================================================
// CEA Biomarker
// ================================================================
// CEA proportional to viable tumor volume
dxdt_CEA_BM = ksynCEA * (corr_TV_CEA * TV_total + 1.0) - kdegCEA * CEA_BM;

// ================================================================
// Cancer Stem Cells
// ================================================================
// CSC self-renewal, inhibited by chemotherapy
double E_CSC_chemo = 0.4 * FU5_C / (0.8 + FU5_C) + 0.3 * DXD_FREE / (0.5 + DXD_FREE);
dxdt_CSC = kCSC_self * (1.0 - E_CSC_chemo) * CSC
           - kCSC_diff * CSC - kdegCSC * CSC;

$TABLE
// ---- Derived outputs for plotting ----
double TumorVol    = TV + TV_DMG1 + TV_DMG2;
double Conc_tras   = TRAS_C;
double Conc_ramu   = RAMU_C;
double Conc_nivo   = NIVO_C;
double Conc_FU5    = FU5_C;
double Conc_dxd    = DXD_FREE;
double Conc_zolbe  = ZOLBE_C;
double PD1_occupancy = NIVO_B / PD1_tot * 100.0;  // % PD-1 occupied
double HER2_free_pct = (HER2_total - TRAS_B - TDXD_B) / HER2_total * 100.0;
double VEGFR2_blocked= RAMU_B / VEGFR2_tot * 100.0;  // % VEGFR2 blocked

$CAPTURE TumorVol Conc_tras Conc_ramu Conc_nivo Conc_FU5 Conc_dxd
         Conc_zolbe PD1_occupancy HER2_free_pct VEGFR2_blocked
         CEA_BM CD8_TEFF TREG TAM_M2 HER2_SIG VEGF_FREE CSC
'

## ============================================================================
## COMPILE MODEL
## ============================================================================
gc_mod <- mcode("gastric_cancer_qsp", gc_model_code)

## ============================================================================
## TREATMENT SCENARIOS
## ============================================================================

## Helper: create dosing events
make_doses <- function(drug_cmt, dose, interval, n_doses, tinf = 0.5) {
  # dose in nmol for mAbs, mg for chemo
  ev(cmt = drug_cmt, amt = dose, ii = interval, addl = n_doses - 1,
     rate = dose / tinf)
}

## Cycle lengths (days)
# FLOT: 14-day cycle (EOW) x 8 cycles perioperative = 112 days total
# FOLFOX: 14-day cycle x 12 cycles = 168 days
# Ramucirumab: 14-day cycle
# Nivolumab: 14-day (Q2W) or 28-day (Q4W)

## ============================================================================
## Scenario 1: FLOT perioperative (5-FU + Leucovorin + Oxaliplatin + Docetaxel)
## Ref: FLOT4 (Al-Batran et al. Lancet 2019)
## 4 pre-op cycles + 4 post-op cycles (every 2 weeks)
## Primary endpoint: OS median 50 vs 35 months (FLOT vs ECF/ECX)
## ============================================================================
scenario1_FLOT <- function(mod) {
  # Capecitabine proxy for 5-FU (oral, 1250 mg/m2 BID d1-14 for CAPE arm)
  # For FLOT IV 5-FU 2600 mg/m2 (24h infusion) + leucovorin 200 mg/m2
  # Simplified as capecitabine oral equivalent
  cape_dose <- 3500  # mg (approximate dose)
  cape_ev   <- ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13)  # 14 days

  # 8 cycles, 14-day intervals
  ev_list <- lapply(0:7, function(cycle) {
    start_day <- cycle * 14
    ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13, time = start_day)
  })
  ev_flot <- do.call(c, ev_list)

  sim <- mod %>%
    ev(ev_flot) %>%
    mrgsim(end = 168, delta = 0.5)

  as_tibble(sim) %>% mutate(scenario = "Scenario 1: FLOT Perioperative")
}

## ============================================================================
## Scenario 2: Trastuzumab + FOLFOX/XELOX (HER2+ first-line)
## Ref: ToGA trial (Bang et al. Lancet 2010) — mOS 13.8 vs 11.1 mo
## Ref: KEYNOTE-811 update (2024) — Pembro + Tras + FOLFOX, ORR 74.4%
## Dosing: Trastuzumab 8mg/kg loading → 6mg/kg Q3W IV
##         FOLFOX6: oxaliplatin 100mg/m2 + 5-FU 400mg/m2 bolus + 2400mg/m2 46h Q2W
## ============================================================================
scenario2_Tras_FOLFOX <- function(mod) {
  # Convert mg/kg to nmol (MW trastuzumab = 148,000 g/mol, BSA 1.7m2, weight 65kg)
  # Loading dose 8mg/kg x 65kg = 520mg = 520,000 μg / 148,000 g/mol = 3.51 nmol
  # Adjusted to plasma concentration: ~3.51/2.91 nmol/L
  tras_loading <- 3.51  # nmol
  tras_maint   <- 2.63  # nmol (6 mg/kg)

  ev_tras <- ev(cmt = "TRAS_C", amt = tras_loading, time = 0, rate = tras_loading/0.5) +
             ev(cmt = "TRAS_C", amt = tras_maint, ii = 21, addl = 5, time = 21, rate = tras_maint/0.5)

  cape_dose <- 3500  # mg XELOX cape equivalent
  ev_list   <- lapply(0:5, function(cycle) {
    ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13, time = cycle * 21)
  })
  ev_cape   <- do.call(c, ev_list)

  ev_combo  <- ev_tras + ev_cape

  sim <- mod %>%
    param(HER2_total = 4.0) %>%  # HER2+ high expression
    ev(ev_combo) %>%
    mrgsim(end = 168, delta = 0.5)

  as_tibble(sim) %>% mutate(scenario = "Scenario 2: Trastuzumab + FOLFOX (HER2+ 1L)")
}

## ============================================================================
## Scenario 3: Ramucirumab + Paclitaxel (second-line)
## Ref: RAINBOW trial (Wilke et al. Lancet Oncol 2014)
## mOS 9.6 vs 7.4 mo (ramu+pac vs pac); ORR 28% vs 16%
## Dosing: Ramucirumab 8mg/kg D1,D15 Q28D; Paclitaxel 80mg/m2 D1,8,15 Q28D
## ============================================================================
scenario3_Ramu_Pac <- function(mod) {
  # Ramucirumab 8mg/kg x 65kg = 520mg / MW 147,000 = 3.54 nmol
  ramu_dose <- 3.54  # nmol
  # Q2W dosing within 28-day cycle
  ev_ramu <- ev(cmt = "RAMU_C", amt = ramu_dose, ii = 14, addl = 7, time = 0,
                rate = ramu_dose/0.5)

  # Paclitaxel as chemo proxy (simplified via CAPE for TGI)
  pac_dose  <- 2800  # mg equivalent
  ev_pac_list <- lapply(0:3, function(cycle) {
    base <- cycle * 28
    ev(cmt = "CAPE_GUT", amt = pac_dose/3, time = base) +
    ev(cmt = "CAPE_GUT", amt = pac_dose/3, time = base + 7) +
    ev(cmt = "CAPE_GUT", amt = pac_dose/3, time = base + 14)
  })
  ev_pac <- do.call(c, ev_pac_list)
  ev_combo <- ev_ramu + ev_pac

  sim <- mod %>%
    param(TV0 = 1500) %>%  # 2L = larger tumor burden
    ev(ev_combo) %>%
    mrgsim(end = 168, delta = 0.5)

  as_tibble(sim) %>% mutate(scenario = "Scenario 3: Ramucirumab + Paclitaxel (2L)")
}

## ============================================================================
## Scenario 4: Nivolumab + Chemotherapy (CPS≥5 first-line)
## Ref: CheckMate 649 (Janjigian et al. Nature 2021)
## CPS≥5: mOS 14.4 vs 11.1 mo; mPFS 7.7 vs 6.0 mo; ORR 60% vs 45%
## Dosing: Nivolumab 360mg Q3W + FOLFOX6/XELOX
## ============================================================================
scenario4_Nivo_Chemo <- function(mod) {
  # Nivolumab 360mg / MW 146,000 = 2.47 nmol → per V1=3.75L = 0.66 nmol/L
  nivo_dose <- 2.47  # nmol
  ev_nivo   <- ev(cmt = "NIVO_C", amt = nivo_dose, ii = 21, addl = 7, time = 0,
                  rate = nivo_dose/0.5)

  cape_dose <- 3500  # mg
  ev_list   <- lapply(0:7, function(cycle) {
    ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13, time = cycle * 21)
  })
  ev_cape   <- do.call(c, ev_list)
  ev_combo  <- ev_nivo + ev_cape

  sim <- mod %>%
    param(PD1_tot = 4.0,   # CPS≥5 = higher PD-L1
          kexhCD8 = 0.025) %>%
    ev(ev_combo) %>%
    mrgsim(end = 252, delta = 0.5)  # extended follow-up

  as_tibble(sim) %>% mutate(scenario = "Scenario 4: Nivolumab + Chemo (CPS≥5, 1L)")
}

## ============================================================================
## Scenario 5: T-DXd (Trastuzumab Deruxtecan) — HER2+ second-line
## Ref: DESTINY-Gastric01 (Shitara et al. NEJM 2020)
## ORR 51.3% vs 14.3% (T-DXd vs physician's choice); mOS 12.5 vs 8.4 mo
## Dosing: T-DXd 6.4mg/kg Q3W IV
## ============================================================================
scenario5_TDXd <- function(mod) {
  # T-DXd 6.4mg/kg x 65kg = 416mg / MW 184,000 = 2.26 nmol
  tdxd_dose <- 2.26  # nmol
  ev_tdxd   <- ev(cmt = "TDXD_C", amt = tdxd_dose, ii = 21, addl = 7, time = 0,
                  rate = tdxd_dose/0.5)

  sim <- mod %>%
    param(HER2_total = 2.5,  # HER2 low-mod (IHC 1+/2+) also eligible
          TV0 = 2000) %>%    # 2L larger tumor
    ev(ev_tdxd) %>%
    mrgsim(end = 252, delta = 0.5)

  as_tibble(sim) %>% mutate(scenario = "Scenario 5: T-DXd (HER2+ 2L)")
}

## ============================================================================
## Scenario 6: Zolbetuximab + mFOLFOX6 (CLDN18.2+ first-line)
## Ref: SPOTLIGHT (Shitara et al. Lancet 2023) — mOS 18.2 vs 15.5 mo
## Ref: GLOW (Shah et al. ESMO 2023) — mOS 14.4 vs 12.2 mo (CAPOX)
## Dosing: Zolbetuximab 800mg/m2 loading → 600mg/m2 Q3W + mFOLFOX6 Q2W
## ============================================================================
scenario6_Zolbe_mFOLFOX <- function(mod) {
  # Zolbetuximab 800mg/m2 x 1.7m2 = 1360mg loading / MW ~148,000 = 9.19 nmol
  zolbe_load  <- 9.19  # nmol loading
  zolbe_maint <- 6.89  # nmol (600mg/m2)

  ev_zolbe_load  <- ev(cmt = "ZOLBE_C", amt = zolbe_load, time = 0, rate = zolbe_load/1.5)
  ev_zolbe_maint <- ev(cmt = "ZOLBE_C", amt = zolbe_maint, ii = 21, addl = 7, time = 21,
                       rate = zolbe_maint/1.5)
  ev_zolbe <- ev_zolbe_load + ev_zolbe_maint

  cape_dose <- 3500  # mg (mFOLFOX6 proxy)
  ev_list   <- lapply(0:7, function(cycle) {
    ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13, time = cycle * 21)
  })
  ev_cape  <- do.call(c, ev_list)
  ev_combo <- ev_zolbe + ev_cape

  sim <- mod %>%
    param(CLDN_tot = 5.5,   # CLDN18.2+ high expression
          HER2_total = 0.5) %>%  # HER2- assumed for CLDN18.2+ selection
    ev(ev_combo) %>%
    mrgsim(end = 252, delta = 0.5)

  as_tibble(sim) %>% mutate(scenario = "Scenario 6: Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)")
}

## ============================================================================
## SIMULATE ALL SCENARIOS
## ============================================================================
simulate_all_scenarios <- function() {
  cat("Compiling gastric cancer QSP model...\n")
  mod <- gc_mod

  cat("Running Scenario 1: FLOT perioperative...\n")
  s1 <- tryCatch(scenario1_FLOT(mod), error = function(e) {
    message("Scenario 1 error: ", e$message); NULL
  })

  cat("Running Scenario 2: Trastuzumab + FOLFOX (HER2+ 1L)...\n")
  s2 <- tryCatch(scenario2_Tras_FOLFOX(mod), error = function(e) {
    message("Scenario 2 error: ", e$message); NULL
  })

  cat("Running Scenario 3: Ramucirumab + Paclitaxel (2L)...\n")
  s3 <- tryCatch(scenario3_Ramu_Pac(mod), error = function(e) {
    message("Scenario 3 error: ", e$message); NULL
  })

  cat("Running Scenario 4: Nivolumab + Chemo (CPS>=5, 1L)...\n")
  s4 <- tryCatch(scenario4_Nivo_Chemo(mod), error = function(e) {
    message("Scenario 4 error: ", e$message); NULL
  })

  cat("Running Scenario 5: T-DXd (HER2+ 2L)...\n")
  s5 <- tryCatch(scenario5_TDXd(mod), error = function(e) {
    message("Scenario 5 error: ", e$message); NULL
  })

  cat("Running Scenario 6: Zolbetuximab + mFOLFOX6 (CLDN18.2+ 1L)...\n")
  s6 <- tryCatch(scenario6_Zolbe_mFOLFOX(mod), error = function(e) {
    message("Scenario 6 error: ", e$message); NULL
  })

  results <- bind_rows(Filter(Negate(is.null), list(s1, s2, s3, s4, s5, s6)))
  return(results)
}

## ============================================================================
## PLOTTING FUNCTIONS
## ============================================================================

plot_tumor_kinetics <- function(results) {
  p <- results %>%
    ggplot(aes(x = time, y = TumorVol, color = scenario)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 500, linetype = "dashed", color = "gray50") +
    scale_y_continuous(trans = "log10",
                       labels = scales::comma,
                       breaks = c(100, 500, 1000, 2500, 5000, 10000)) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Gastric Cancer — Tumor Volume Kinetics Across Treatment Scenarios",
      subtitle = "Simeoni TGI model with Bliss independence drug combination",
      x = "Time (days)",
      y = "Tumor Volume (mm³, log scale)",
      color = "Treatment Scenario",
      caption = "Reference lines: TV0 = 500 mm³. Parameters from ToGA, RAINBOW, CheckMate 649, SPOTLIGHT, DESTINY-Gastric01"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8))
  return(p)
}

plot_pk_curves <- function(results) {
  pk_data <- results %>%
    select(time, scenario, Conc_tras, Conc_ramu, Conc_nivo, Conc_FU5, Conc_dxd, Conc_zolbe) %>%
    pivot_longer(cols = c(Conc_tras, Conc_ramu, Conc_nivo, Conc_FU5, Conc_dxd, Conc_zolbe),
                 names_to = "Drug", values_to = "Concentration") %>%
    mutate(Drug = recode(Drug,
      "Conc_tras"  = "Trastuzumab (nmol/L)",
      "Conc_ramu"  = "Ramucirumab (nmol/L)",
      "Conc_nivo"  = "Nivolumab (nmol/L)",
      "Conc_FU5"   = "5-FU active (μmol/L)",
      "Conc_dxd"   = "DXd payload (μmol/L)",
      "Conc_zolbe" = "Zolbetuximab (nmol/L)"
    ))

  p <- pk_data %>%
    filter(Concentration > 0.001) %>%
    ggplot(aes(x = time, y = Concentration, color = scenario)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~Drug, scales = "free_y", ncol = 3) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Drug PK Profiles — Gastric Cancer Treatment Scenarios",
      x = "Time (days)", y = "Concentration",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
  return(p)
}

plot_immune_dynamics <- function(results) {
  imm_data <- results %>%
    select(time, scenario, CD8_TEFF, TREG, TAM_M2, PD1_occupancy) %>%
    pivot_longer(cols = c(CD8_TEFF, TREG, TAM_M2, PD1_occupancy),
                 names_to = "Immune", values_to = "Value") %>%
    mutate(Immune = recode(Immune,
      "CD8_TEFF"      = "CD8+ T Effector (cells/μL)",
      "TREG"          = "Treg (cells/μL)",
      "TAM_M2"        = "TAM M2 (cells/μL)",
      "PD1_occupancy" = "PD-1 Occupancy (%)"
    ))

  p <- imm_data %>%
    ggplot(aes(x = time, y = Value, color = scenario)) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~Immune, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Immune Cell Dynamics in Gastric Cancer TME",
      x = "Time (days)", y = "Value",
      color = "Scenario",
      caption = "PD-1 occupancy reflects nivolumab/pembrolizumab binding to PD-1 receptor"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  return(p)
}

plot_biomarkers <- function(results) {
  bm_data <- results %>%
    select(time, scenario, CEA_BM, HER2_SIG, VEGF_FREE) %>%
    pivot_longer(cols = c(CEA_BM, HER2_SIG, VEGF_FREE),
                 names_to = "Biomarker", values_to = "Value") %>%
    mutate(Biomarker = recode(Biomarker,
      "CEA_BM"    = "CEA (ng/mL)",
      "HER2_SIG"  = "HER2 Signaling (AU)",
      "VEGF_FREE" = "Free VEGF-A (ng/mL)"
    ))

  p <- bm_data %>%
    ggplot(aes(x = time, y = Value, color = scenario)) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~Biomarker, scales = "free_y") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Biomarker Dynamics — Gastric Cancer QSP Model",
      x = "Time (days)", y = "Value",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  return(p)
}

## ============================================================================
## PARAMETER SENSITIVITY ANALYSIS
## ============================================================================

sensitivity_analysis <- function(mod) {
  cat("Running parameter sensitivity analysis...\n")

  # Vary HER2 expression level for scenario 2
  her2_levels <- c(0.5, 1.0, 2.5, 4.0, 6.0)

  sens_results <- lapply(her2_levels, function(her2) {
    cape_dose <- 3500
    tras_loading <- 3.51
    tras_maint   <- 2.63

    ev_tras <- ev(cmt = "TRAS_C", amt = tras_loading, time = 0, rate = tras_loading/0.5) +
               ev(cmt = "TRAS_C", amt = tras_maint, ii = 21, addl = 5, time = 21, rate = tras_maint/0.5)

    ev_list <- lapply(0:5, function(cycle) {
      ev(cmt = "CAPE_GUT", amt = cape_dose, ii = 1, addl = 13, time = cycle * 21)
    })
    ev_cape <- do.call(c, ev_list)

    sim <- mod %>%
      param(HER2_total = her2) %>%
      ev(ev_tras + ev_cape) %>%
      mrgsim(end = 168, delta = 1.0)

    as_tibble(sim) %>% mutate(HER2_level = her2)
  })

  do.call(rbind, sens_results)
}

## ============================================================================
## SUMMARY: Response Rates at Day 84 (12 weeks)
## ============================================================================

compute_response_summary <- function(results) {
  results %>%
    group_by(scenario) %>%
    summarize(
      TV_baseline = first(TumorVol),
      TV_week12   = approx(time, TumorVol, xout = 84)$y,
      TV_week24   = approx(time, TumorVol, xout = 168)$y,
      CEA_baseline = first(CEA_BM),
      CEA_week12   = approx(time, CEA_BM, xout = 84)$y,
      .groups = "drop"
    ) %>%
    mutate(
      pct_change_12w = (TV_week12 - TV_baseline) / TV_baseline * 100,
      pct_change_24w = (TV_week24 - TV_baseline) / TV_baseline * 100,
      RECIST_response = case_when(
        pct_change_12w <= -30 ~ "Partial Response (PR)",
        pct_change_12w <= 20  ~ "Stable Disease (SD)",
        TRUE                  ~ "Progressive Disease (PD)"
      ),
      CEA_response = case_when(
        (CEA_week12 - CEA_baseline) / CEA_baseline * 100 <= -30 ~ "CEA responder",
        TRUE ~ "CEA non-responder"
      )
    )
}

## ============================================================================
## MAIN EXECUTION
## ============================================================================

if(interactive()) {
  cat("=== Gastric Cancer QSP Model ===\n")
  cat("Compiling model and running all 6 treatment scenarios...\n\n")

  results <- simulate_all_scenarios()

  cat("\n=== Response Summary at 12 and 24 weeks ===\n")
  summary_tbl <- compute_response_summary(results)
  print(summary_tbl, n = Inf)

  cat("\n=== Generating plots ===\n")
  p1 <- plot_tumor_kinetics(results)
  p2 <- plot_pk_curves(results)
  p3 <- plot_immune_dynamics(results)
  p4 <- plot_biomarkers(results)

  print(p1)
  print(p2)
  print(p3)
  print(p4)

  cat("\nDone! Gastric Cancer QSP simulation complete.\n")
}
