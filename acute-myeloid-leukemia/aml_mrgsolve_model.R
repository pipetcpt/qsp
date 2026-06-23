# =============================================================================
# Acute Myeloid Leukemia (AML) QSP mrgsolve Model
# =============================================================================
# File   : aml_mrgsolve_model.R
# Version: 2.0
# Date   : 2026-06-23
# Author : Claude Code Routine (CCR) — Auto-generated QSP model
#
# Description:
#   A multi-drug, multi-compartment QSP model for AML capturing:
#     - Leukemia stem cell (LSC) → progenitor (LPC) → blast (LBC) hierarchy
#     - Drug PK: venetoclax (3-cmt), azacitidine (2-cmt SC), gilteritinib
#       (3-cmt), enasidenib (2-cmt), cytarabine (1-cmt IV),
#       ATRA (1-cmt, APL), ATO (1-cmt, APL)
#     - Drug PD: BCL-2 occupancy, FLT3 inhibition, DNMT inhibition/
#       differentiation induction, S-phase cytotoxicity,
#       ATRA/ATO-driven APL differentiation and apoptosis
#     - Myelosuppression (Friberg semi-mechanistic): ANC, platelets, Hgb
#     - MRD (log10 scale), bone marrow blast %, CR/MRD-neg status flags
#
# Key Clinical Trials Parameterized Against:
#   - VIALE-A (DiNardo 2020 NEJM): VEN + AZA, newly diagnosed unfit AML
#   - ADMIRAL (Perl 2019 NEJM): gilteritinib vs salvage chemo FLT3+ R/R AML
#   - IDHENTIFY (Stein 2019 Lancet Oncol): enasidenib IDH2+ R/R AML
#   - 7+3 standard induction (Dohner 2017 NEJM guidelines)
#   - VIALE-C (Wei 2020 JCO): VEN + LDAC elderly AML
#   - APL (Lo-Coco 2013 NEJM): ATRA + ATO
#
# References:
#   - Friberg LE et al. J Clin Oncol 2002; 20(24):4713-4721
#   - Gibiansky E et al. CPT Pharmacometrics Syst Pharmacol 2022
#   - Hamed SS et al. CPT Pharmacometrics Syst Pharmacol 2020
#   - Salem AH et al. Leukemia 2017; 31(9):2059-2066 (VEN PK)
#   - Loke J et al. Br J Haematol 2021 (gilteritinib PK/PD)
#   - Laille E et al. J Clin Pharmacol 2014 (AZA PK)
#   - Yen K et al. Leukemia 2017 (enasidenib PK)
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# =============================================================================
# MODEL CODE — passed to mcode()
# =============================================================================

AML_model_code <- '
[PROB]
// ==========================================================================
// Acute Myeloid Leukemia (AML) QSP Model
// ==========================================================================
// Multi-drug QSP model for AML covering:
//   * Leukemia cell hierarchy: LSC -> LPC -> LBC (blasts)
//   * Drug PK: venetoclax (3-cmt), azacitidine (2-cmt SC),
//              gilteritinib (3-cmt), enasidenib (2-cmt), cytarabine (1-cmt IV)
//              ATRA (1-cmt oral), ATO (1-cmt IV/oral)
//   * Drug PD: BCL-2 occupancy, FLT3 inhibition, DNMT inhibition,
//              S-phase cytotoxicity, ATRA/ATO APL differentiation
//   * Myelosuppression: Friberg semi-mechanistic ANC/PLT/Hgb
//   * Efficacy outputs: BM blast %, CR status, MRD (log10)
// ==========================================================================

[PARAM]
// --------------------------------------------------------------------------
// VENETOCLAX PK (3-compartment, oral, food-effect)
// Salem AH et al. Leukemia 2017; CPT:PSP 2019
// --------------------------------------------------------------------------
VEN_CL       = 15.6     // L/h  apparent clearance (fed state)
VEN_V1       = 98.4     // L    central volume
VEN_V2       = 289.0    // L    peripheral volume
VEN_Q        = 11.2     // L/h  inter-compartmental CL
VEN_ka       = 0.54     // 1/h  first-order absorption rate
VEN_F        = 0.57     // --   bioavailability (fed, low-fat meal)
VEN_F_fasted = 0.22     // --   bioavailability (fasted)

// --------------------------------------------------------------------------
// AZACITIDINE PK (2-compartment SC)
// Laille E et al. J Clin Pharmacol 2014
// --------------------------------------------------------------------------
AZA_CL = 147.0   // L/h apparent clearance
AZA_V1 = 76.0    // L   central volume
AZA_V2 = 118.0   // L   peripheral volume
AZA_Q  = 68.0    // L/h inter-compartmental CL
AZA_ka = 2.1     // 1/h SC absorption rate constant
AZA_F  = 0.89    // --  SC bioavailability

// --------------------------------------------------------------------------
// GILTERITINIB PK (3-compartment oral)
// Loke J et al. Br J Haematol 2021; pop-PK analysis
// --------------------------------------------------------------------------
GILT_CL = 14.5    // L/h apparent clearance
GILT_V1 = 1093.0  // L   central volume
GILT_V2 = 814.0   // L   peripheral volume
GILT_Q  = 7.2     // L/h
GILT_ka = 0.32    // 1/h

// --------------------------------------------------------------------------
// ENASIDENIB PK (2-compartment oral)
// Yen K et al. Leukemia 2017; pop-PK model
// --------------------------------------------------------------------------
ENASID_CL = 1.86  // L/h
ENASID_V  = 55.8  // L
ENASID_ka = 0.41  // 1/h

// --------------------------------------------------------------------------
// CYTARABINE PK (IV, 1-compartment)
// Capizzi RL et al. Semin Hematol 1991
// --------------------------------------------------------------------------
CYTARAB_CL = 47.0  // L/h
CYTARAB_V  = 32.0  // L

// --------------------------------------------------------------------------
// VENETOCLAX PD -- BCL-2 inhibition / apoptosis induction
// Souers AJ et al. Nat Med 2013; DiNardo CD et al. NEJM 2020
// --------------------------------------------------------------------------
VEN_EC50 = 0.15   // ug/mL BCL-2 occupancy EC50
VEN_Emax = 0.92   // --    max apoptosis induction fraction
VEN_hill = 1.5    // --    Hill coefficient

// --------------------------------------------------------------------------
// GILTERITINIB PD -- FLT3 kinase inhibition
// Lee LY et al. Nat Med 2017
// --------------------------------------------------------------------------
GILT_IC50 = 0.05  // ug/mL FLT3 inhibition IC50 (pFLT3)
GILT_Imax = 0.88  // --    max FLT3 pathway inhibition

// --------------------------------------------------------------------------
// AZACITIDINE PD -- DNMT inhibition -> differentiation induction
// Silverman LR et al. J Clin Oncol 2006
// --------------------------------------------------------------------------
AZA_IC50 = 0.8    // ug/mL differentiation induction EC50
AZA_Imax = 0.75   // --    max differentiation induction

// --------------------------------------------------------------------------
// LEUKEMIC STEM CELL (LSC) DYNAMICS
// Leder K et al. PLoS Comput Biol 2010; Craig M et al. J Pharmacokinet 2016
// --------------------------------------------------------------------------
k_LSC_prolif = 0.18   // 1/day net self-renewal rate constant
K_LSC_carry  = 1e9    // cells carrying capacity
k_LSC_diff   = 0.08   // 1/day differentiation rate to LPC
k_LSC_apop   = 0.02   // 1/day baseline apoptosis rate

// --------------------------------------------------------------------------
// LEUKEMIC PROGENITOR CELL (LPC) DYNAMICS
// --------------------------------------------------------------------------
k_LPC_prolif = 0.35   // 1/day proliferation rate
k_LPC_diff   = 0.12   // 1/day differentiation to blasts
k_LPC_apop   = 0.05   // 1/day baseline apoptosis

// --------------------------------------------------------------------------
// LEUKEMIC BLAST (LBC) DYNAMICS -- peripheral/BM blasts
// --------------------------------------------------------------------------
k_LBC_clear  = 0.15   // 1/day clearance / senescent death
k_LBC_apop   = 0.04   // 1/day baseline apoptosis

// --------------------------------------------------------------------------
// NORMAL HEMATOPOIESIS -- Friberg semi-mechanistic myelosuppression
// Friberg LE et al. J Clin Oncol 2002; 20:4713-4721
// --------------------------------------------------------------------------
k_ANC_base  = 3.5     // x10^9/L baseline ANC
k_PLT_base  = 200.0   // x10^9/L baseline platelets
k_Hgb_base  = 13.5    // g/dL   baseline hemoglobin
ANC_MTT     = 5.9     // days   ANC mean transit time
PLT_MTT     = 7.1     // days   platelet mean transit time
k_circ_0    = 5.3     // 1/day  baseline proliferating ANC

// --------------------------------------------------------------------------
// MRD PARAMETERS
// --------------------------------------------------------------------------
MRD_sensitivity = -4.5   // log10 PCR assay sensitivity floor

// --------------------------------------------------------------------------
// TUMOR LYSIS SYNDROME / TOXICITY THRESHOLDS
// --------------------------------------------------------------------------
TLS_threshold = 25.0    // % blast elevated TLS risk above this level

// --------------------------------------------------------------------------
// MUTATION STATUS FLAGS (1 = present)
// --------------------------------------------------------------------------
FLT3_status  = 1   // FLT3-ITD or TKD mutation
IDH2_status  = 0   // IDH2 mutation
NPM1_status  = 1   // NPM1 mutation
APL_status   = 0   // PML-RARA (APL subtype)

// --------------------------------------------------------------------------
// TREATMENT REGIMEN FLAGS (1 = active)
// --------------------------------------------------------------------------
use_VEN     = 1
use_AZA     = 1
use_GILT    = 0
use_ENASID  = 0
use_CYTARAB = 0

// ATRA/ATO parameters (APL)
ATRA_Emax  = 0.90
ATRA_EC50  = 0.12   // ug/mL
ATO_Emax   = 0.88
ATO_EC50   = 0.08   // ug/mL

// Body surface area for dose normalization
BSA = 1.8

[CMT]
// Drug PK compartments
VEN_gut      // venetoclax GI depot (mg)
VEN_cent     // venetoclax central (mg)
VEN_peri     // venetoclax peripheral (mg)
AZA_depot    // azacitidine SC depot (mg)
AZA_cent     // azacitidine central (mg)
AZA_peri     // azacitidine peripheral (mg)
GILT_gut     // gilteritinib GI depot (mg)
GILT_cent    // gilteritinib central (mg)
GILT_peri    // gilteritinib peripheral (mg)
ENASID_gut   // enasidenib GI depot (mg)
ENASID_cent  // enasidenib central (mg)
CYTARAB_cent // cytarabine central (mg)

// ATRA / ATO (APL arms) -- simple 1-cmt
ATRA_cent    // all-trans retinoic acid central (mg)
ATO_cent     // arsenic trioxide central (mg)

// Disease -- leukemia cell hierarchy (x10^6 cells)
LSC          // leukemic stem cells
LPC          // leukemic progenitor cells
LBC          // leukemic blasts (circulating/BM)

// PD intermediates
BCL2_occ     // BCL-2 occupancy fraction (0-1)

// Normal hematopoiesis -- Friberg ANC chain
ANC_prol     // proliferating cells
ANC_trans1   // transit compartment 1
ANC_trans2   // transit compartment 2
ANC_circ     // circulating ANC (x10^9/L)

// Platelets and hemoglobin
PLT_circ     // circulating platelets (x10^9/L)
Hgb_circ     // hemoglobin (g/dL)

// MRD and total tumor
MRD_log      // log10(total tumor cells)
tumor_vol    // total tumor burden (x10^6 cells)

[GLOBAL]
double VEN_C1, AZA_C1, GILT_C1, ENASID_C1, CYTARAB_C1, ATRA_C1, ATO_C1;
double VEN_eff, FLT3_inh, AZA_diff_eff, CYTARAB_kill;
double ATRA_eff, ATO_eff, APL_apop;
double apop_drug, prolif_inh;
double Edrug_ANC, FB_ANC, FB_PLT;

[MAIN]
// Initial conditions at diagnosis
VEN_gut_0      = 0.0;
VEN_cent_0     = 0.0;
VEN_peri_0     = 0.0;
AZA_depot_0    = 0.0;
AZA_cent_0     = 0.0;
AZA_peri_0     = 0.0;
GILT_gut_0     = 0.0;
GILT_cent_0    = 0.0;
GILT_peri_0    = 0.0;
ENASID_gut_0   = 0.0;
ENASID_cent_0  = 0.0;
CYTARAB_cent_0 = 0.0;
ATRA_cent_0    = 0.0;
ATO_cent_0     = 0.0;

LSC_0          = 500.0;   // x10^6 cells (~60-70% BM blasts at dx)
LPC_0          = 2000.0;
LBC_0          = 5000.0;

BCL2_occ_0     = 0.0;

ANC_prol_0     = k_circ_0;
ANC_trans1_0   = k_circ_0;
ANC_trans2_0   = k_circ_0;
ANC_circ_0     = k_ANC_base;
PLT_circ_0     = k_PLT_base;
Hgb_circ_0     = k_Hgb_base;

MRD_log_0      = log10(LSC_0 + LPC_0 + LBC_0);
tumor_vol_0    = LSC_0 + LPC_0 + LBC_0;

[ODE]
// ==========================================================================
// DRUG PK
// ==========================================================================

// --- Venetoclax (3-cmt, oral) ---
VEN_C1 = VEN_cent / VEN_V1;
double VEN_C2 = VEN_peri / VEN_V2;
dxdt_VEN_gut  = -VEN_ka * VEN_gut;
dxdt_VEN_cent =  VEN_F * VEN_ka * VEN_gut
                - (VEN_CL / VEN_V1) * VEN_cent
                - (VEN_Q  / VEN_V1) * VEN_cent
                + (VEN_Q  / VEN_V2) * VEN_peri;
dxdt_VEN_peri =  (VEN_Q  / VEN_V1) * VEN_cent
                - (VEN_Q  / VEN_V2) * VEN_peri;

// --- Azacitidine (2-cmt SC) ---
AZA_C1 = AZA_cent / AZA_V1;
double AZA_C2 = AZA_peri / AZA_V2;
dxdt_AZA_depot =  -AZA_ka * AZA_depot;
dxdt_AZA_cent  =   AZA_F * AZA_ka * AZA_depot
                  - (AZA_CL / AZA_V1) * AZA_cent
                  - (AZA_Q  / AZA_V1) * AZA_cent
                  + (AZA_Q  / AZA_V2) * AZA_peri;
dxdt_AZA_peri  =   (AZA_Q  / AZA_V1) * AZA_cent
                  - (AZA_Q  / AZA_V2) * AZA_peri;

// --- Gilteritinib (3-cmt, oral) ---
GILT_C1 = GILT_cent / GILT_V1;
double GILT_C2 = GILT_peri / GILT_V2;
dxdt_GILT_gut  = -GILT_ka * GILT_gut;
dxdt_GILT_cent =  GILT_ka * GILT_gut
                 - (GILT_CL / GILT_V1) * GILT_cent
                 - (GILT_Q  / GILT_V1) * GILT_cent
                 + (GILT_Q  / GILT_V2) * GILT_peri;
dxdt_GILT_peri =  (GILT_Q  / GILT_V1) * GILT_cent
                 - (GILT_Q  / GILT_V2) * GILT_peri;

// --- Enasidenib (2-cmt, oral) ---
ENASID_C1 = ENASID_cent / ENASID_V;
dxdt_ENASID_gut  = -ENASID_ka * ENASID_gut;
dxdt_ENASID_cent =  ENASID_ka * ENASID_gut
                   - (ENASID_CL / ENASID_V) * ENASID_cent;

// --- Cytarabine (IV, 1-cmt) ---
CYTARAB_C1 = CYTARAB_cent / CYTARAB_V;
dxdt_CYTARAB_cent = -(CYTARAB_CL / CYTARAB_V) * CYTARAB_cent;

// --- ATRA (APL; simplified 1-cmt oral) ---
// CL ~250 L/h, V ~700 L, ka ~0.85 1/h, F ~0.18 (auto-induction at day 14)
ATRA_C1 = ATRA_cent / 700.0;
dxdt_ATRA_cent = -(250.0 / 700.0) * ATRA_cent;

// --- ATO (arsenic trioxide; simplified 1-cmt) ---
// CL ~66 L/h, V ~400 L; plasma half-life ~10h
ATO_C1 = ATO_cent / 400.0;
dxdt_ATO_cent = -(66.0 / 400.0) * ATO_cent;

// ==========================================================================
// PD EFFECT CALCULATIONS
// ==========================================================================

// BCL-2 occupancy by venetoclax (E-max Hill, biophase equilibration tau~0.5h)
double VEN_eff_inst = VEN_Emax * pow(VEN_C1, VEN_hill)
                      / (pow(VEN_EC50, VEN_hill) + pow(VEN_C1, VEN_hill) + 1e-12);
// tau = 0.5h = 0.021 day; biophase ODE
dxdt_BCL2_occ = (VEN_eff_inst - BCL2_occ) / 0.021;
VEN_eff = BCL2_occ;

// FLT3 inhibition by gilteritinib
FLT3_inh = GILT_Imax * GILT_C1 / (GILT_IC50 + GILT_C1 + 1e-12);
FLT3_inh = FLT3_inh * FLT3_status * use_GILT;

// AZA: DNMT inhibition -> differentiation induction
AZA_diff_eff = AZA_Imax * AZA_C1 / (AZA_IC50 + AZA_C1 + 1e-12);
AZA_diff_eff = AZA_diff_eff * use_AZA;

// Cytarabine: S-phase-specific kill (Michaelis-Menten)
CYTARAB_kill = 0.55 * CYTARAB_C1 / (0.8 + CYTARAB_C1 + 1e-12);
CYTARAB_kill = CYTARAB_kill * use_CYTARAB;

// ATRA: differentiation induction in APL (PML-RARa degradation)
ATRA_eff = ATRA_Emax * ATRA_C1 / (ATRA_EC50 + ATRA_C1 + 1e-12);
ATRA_eff = ATRA_eff * APL_status;

// ATO: pro-apoptotic + PML-RARa degradation synergy with ATRA
ATO_eff = ATO_Emax * ATO_C1 / (ATO_EC50 + ATO_C1 + 1e-12);
ATO_eff = ATO_eff * APL_status;

// Combined APL-specific effect (synergistic)
APL_apop = 1.0 - (1.0 - ATRA_eff) * (1.0 - ATO_eff);
if(APL_apop > 0.98) APL_apop = 0.98;

// Total drug-induced apoptosis on leukemia cells
apop_drug = VEN_eff * use_VEN + 0.4 * CYTARAB_kill + APL_apop;
if(apop_drug > 0.97) apop_drug = 0.97;

// FLT3 inhibition reduces proliferation (max 70% reduction)
prolif_inh = 1.0 - 0.70 * FLT3_inh;

// ==========================================================================
// LEUKEMIA CELL HIERARCHY DYNAMICS (time in days)
// ==========================================================================

double LSC_net = (LSC > 0.0) ? LSC : 0.0;
double LPC_net = (LPC > 0.0) ? LPC : 0.0;
double LBC_net = (LBC > 0.0) ? LBC : 0.0;

// Leukemic Stem Cells
// BCL-2 inhibition targets LSC effectively (Salem 2019 Nat Med)
dxdt_LSC = (k_LSC_prolif * prolif_inh
            - k_LSC_diff
            - k_LSC_apop
            - apop_drug * 0.45
            - AZA_diff_eff * 0.08
            ) * LSC_net;

// Leukemic Progenitor Cells
dxdt_LPC = k_LSC_diff * LSC_net
           + (k_LPC_prolif * prolif_inh
              - k_LPC_diff
              - k_LPC_apop
              - apop_drug * 0.65
              ) * LPC_net
           - AZA_diff_eff * LPC_net * 0.20;

// Leukemic Blasts (peripheral + BM)
dxdt_LBC = k_LPC_diff * LPC_net
           - (k_LBC_clear
              + k_LBC_apop
              + apop_drug * 0.90
              + CYTARAB_kill * 0.9
              ) * LBC_net
           - AZA_diff_eff * LBC_net * 0.25;

// MRD (log10 scale; tau = 0.2 day for equilibration)
double total_tumor = LSC_net + LPC_net + LBC_net;
double total_safe  = (total_tumor < 1e-3) ? 1e-3 : total_tumor;
dxdt_MRD_log = (log10(total_safe) - MRD_log) / 0.2;

// Total tumor burden
dxdt_tumor_vol = dxdt_LSC + dxdt_LPC + dxdt_LBC;

// ==========================================================================
// NORMAL HEMATOPOIESIS -- Friberg semi-mechanistic myelosuppression
// Friberg LE et al. J Clin Oncol 2002
// ==========================================================================

double ANC_circ_safe = (ANC_circ > 0.01) ? ANC_circ : 0.01;

// Feedback factor (exponent 0.17 from Friberg 2002)
FB_ANC = pow(k_ANC_base / ANC_circ_safe, 0.17);

// Combined myelosuppressive drug effect
Edrug_ANC = CYTARAB_kill * 0.50
            + (AZA_C1 / (AZA_IC50 + AZA_C1 + 1e-12)) * 0.25 * use_AZA;
if(Edrug_ANC > 0.90) Edrug_ANC = 0.90;

double ktr_ANC = 4.0 / ANC_MTT;

dxdt_ANC_prol   = k_circ_0 * FB_ANC * (1.0 - Edrug_ANC) - ktr_ANC * ANC_prol;
dxdt_ANC_trans1 = ktr_ANC * (ANC_prol   - ANC_trans1);
dxdt_ANC_trans2 = ktr_ANC * (ANC_trans1 - ANC_trans2);
dxdt_ANC_circ   = ktr_ANC * ANC_trans2  - ktr_ANC * ANC_circ;

// Platelets (simplified Friberg 2-transit)
double PLT_safe = (PLT_circ > 1.0) ? PLT_circ : 1.0;
FB_PLT = pow(k_PLT_base / PLT_safe, 0.19);
double ktr_PLT = 4.0 / PLT_MTT;
dxdt_PLT_circ = k_PLT_base * FB_PLT * (1.0 - Edrug_ANC * 0.75)
                - ktr_PLT * PLT_circ;

// Hemoglobin (slow dynamics, ~30-day kinetics)
dxdt_Hgb_circ = (k_Hgb_base * (1.0 - Edrug_ANC * 0.30) - Hgb_circ) / 30.0;

[TABLE]
// Bone marrow blast percentage
double LBC_safe     = (LBC > 0.0) ? LBC : 0.0;
double BM_blast_pct = 100.0 * LBC_safe / (LBC_safe + 8000.0);
if(BM_blast_pct > 100.0) BM_blast_pct = 100.0;
if(BM_blast_pct <   0.0) BM_blast_pct = 0.0;

// Complete Remission (CR): blasts < 5%, ANC >= 1.0, PLT >= 100
double CR_status  = 0.0;
if(BM_blast_pct < 5.0 && ANC_circ >= 1.0 && PLT_circ >= 100.0) CR_status = 1.0;

// CRi (CR with incomplete count recovery)
double CRi_status = 0.0;
if(BM_blast_pct < 5.0 && (ANC_circ < 1.0 || PLT_circ < 100.0)) CRi_status = 1.0;

// MRD negativity (<=10^-3 detection threshold)
double MRD_neg = (MRD_log <= MRD_sensitivity) ? 1.0 : 0.0;

// Drug concentrations (ug/mL)
double VEN_conc     = VEN_cent  / VEN_V1;
double AZA_conc     = AZA_cent  / AZA_V1;
double GILT_conc    = GILT_cent / GILT_V1;
double ENASID_conc  = ENASID_cent / ENASID_V;
double CYTARAB_conc = CYTARAB_cent / CYTARAB_V;
double ATRA_conc    = ATRA_cent / 700.0;
double ATO_conc     = ATO_cent  / 400.0;

// Clinical toxicity flags
double neutropenia_g3   = (ANC_circ < 1.0)  ? 1.0 : 0.0;
double neutropenia_g4   = (ANC_circ < 0.5)  ? 1.0 : 0.0;
double thrombo_g3       = (PLT_circ < 50.0) ? 1.0 : 0.0;
double anemia_transfuse = (Hgb_circ < 8.0)  ? 1.0 : 0.0;
double TLS_risk         = (BM_blast_pct > TLS_threshold) ? 1.0 : 0.0;
double diff_syndrome    = (APL_status > 0.5 && ANC_circ > 12.0) ? 1.0 : 0.0;

[CAPTURE]
BM_blast_pct CR_status CRi_status MRD_log MRD_neg
ANC_circ PLT_circ Hgb_circ
VEN_conc AZA_conc GILT_conc ENASID_conc CYTARAB_conc ATRA_conc ATO_conc
LSC LPC LBC tumor_vol
TLS_risk neutropenia_g3 neutropenia_g4 thrombo_g3 anemia_transfuse diff_syndrome
BCL2_occ FLT3_inh AZA_diff_eff CYTARAB_kill
ATRA_eff ATO_eff APL_apop
'

# =============================================================================
# BUILD MODEL
# =============================================================================

mod <- mcode("AML_QSP", AML_model_code, compile = TRUE)

cat("Model compiled successfully.\n")
cat("Compartments:", length(mod@cmtL), "\n")
cat("Parameters:", nrow(param(mod)), "\n")

# =============================================================================
# ============ R SIMULATION CODE ============
# =============================================================================

# Color palette for consistent plotting across scenarios
drug_cols <- c(
  "7+3 Induction"         = "#E41A1C",
  "VEN + AZA (VIALE-A)"   = "#377EB8",
  "Gilteritinib 120mg"    = "#4DAF4A",
  "Enasidenib 100mg"      = "#984EA3",
  "VEN+AZA+GILT (Triple)" = "#FF7F00",
  "LDAC + VEN"            = "#A65628",
  "ATRA + ATO (APL)"      = "#F781BF"
)

# Helper: pull a single row at time point closest to t
get_time_point <- function(df, t) {
  df %>%
    filter(abs(time - t) == min(abs(time - t))) %>%
    slice(1)
}

# Helper: summarise key endpoints
summarise_response <- function(sim_df, scenario_name) {
  d28  <- get_time_point(sim_df, 28)
  d168 <- get_time_point(sim_df, 168)
  tibble(
    Scenario     = scenario_name,
    CR_d28       = d28$CR_status,
    CRi_d28      = if ("CRi_status" %in% names(d28)) d28$CRi_status else 0,
    Blast_d28    = round(d28$BM_blast_pct, 1),
    MRD_d168     = round(d168$MRD_log, 2),
    MRD_neg_d168 = d168$MRD_neg,
    ANC_nadir    = round(min(sim_df$ANC_circ), 2),
    PLT_nadir    = round(min(sim_df$PLT_circ), 1)
  )
}

# =============================================================================
# SCENARIO 1: Standard 7+3 Induction
#   Cytarabine 200 mg/m2 CI d1-7 + Idarubicin 12 mg/m2 d1-3
#   Simulation: 28 days (1 induction cycle)
#   Reference: Dohner H et al. NEJM 2017; Stone RM et al. NEJM 2017
# =============================================================================
cat("\n=== Scenario 1: Standard 7+3 Induction ===\n")

# Cytarabine CI: 200 mg/m2 x 1.8 m2 = 360 mg/day as continuous infusion
# Modelled as half-hourly micro-boluses into central compartment
cytar_daily_mg  <- 200 * 1.8          # 360 mg/day
cytar_times_h   <- seq(0, 7*24 - 0.5, by = 0.5)  # hours over 7 days
cytar_times_d   <- cytar_times_h / 24              # convert to days
cytar_per_bolus <- cytar_daily_mg * 7 / length(cytar_times_d)

e1_cytar <- ev(
  time = cytar_times_d,
  amt  = cytar_per_bolus,
  cmt  = "CYTARAB_cent"
)

# Idarubicin 12 mg/m2 d1-3: modelled as cytarabine-equivalent kill boost
# IDR 12 mg/m2 x 1.8 m2 x potency multiplier ~1.5
idar_equiv   <- 1.8 * 12 * 1.5
idar_d13_d   <- c(0, 1, 2)

e1_idar <- ev(
  time = idar_d13_d,
  amt  = idar_equiv,
  cmt  = "CYTARAB_cent"
)

e1 <- e1_cytar + e1_idar

p1 <- param(mod,
  use_VEN     = 0,
  use_AZA     = 0,
  use_GILT    = 0,
  use_ENASID  = 0,
  use_CYTARAB = 1,
  FLT3_status = 1,
  IDH2_status = 0,
  APL_status  = 0
)

sim1 <- mrgsim(p1, events = e1, end = 28, delta = 0.25, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "7+3 Induction", CRi_status = 0)

cat("Day 28 BM blast%:", get_time_point(sim1, 28)$BM_blast_pct %>% round(1), "\n")
cat("CR at day 28:", get_time_point(sim1, 28)$CR_status, "\n")
cat("ANC nadir:", min(sim1$ANC_circ) %>% round(2), "x10^9/L\n")

# --- Plot Scenario 1 ---
p1_blast <- ggplot(sim1, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#E41A1C", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  annotate("text", x = 18, y = 7.5, label = "CR threshold (5%)", size = 3.2, color = "darkgreen") +
  labs(title = "BM Blast % — 7+3 Induction", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12) +
  coord_cartesian(ylim = c(0, 80))

p1_anc <- ggplot(sim1, aes(x = time)) +
  geom_line(aes(y = ANC_circ), color = "#E41A1C", linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "darkgreen") +
  annotate("text", x = 20, y = 0.7, label = "G3 (<1.0)", size = 3, color = "darkgreen") +
  annotate("text", x = 20, y = 0.3, label = "G4 (<0.5)", size = 3, color = "orange") +
  labs(title = "ANC — Friberg Model", x = "Time (days)", y = "ANC (x10^9/L)") +
  theme_bw(base_size = 12)

p1_plt <- ggplot(sim1, aes(x = time)) +
  geom_line(aes(y = PLT_circ), color = "#E41A1C", linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  annotate("text", x = 20, y = 60, label = "G3 (<50)", size = 3, color = "red") +
  labs(title = "Platelets", x = "Time (days)", y = "PLT (x10^9/L)") +
  theme_bw(base_size = 12)

fig1 <- (p1_blast | p1_anc | p1_plt) +
  plot_annotation(
    title    = "Scenario 1: Standard 7+3 Induction",
    subtitle = "Cytarabine 200 mg/m2 CI d1-7 + Idarubicin 12 mg/m2 d1-3",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig1)
ggsave("aml_scenario1_73induction.png", fig1, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 2: VEN + AZA (VIALE-A Regimen)
#   Venetoclax 400mg PO QD d1-28 each cycle (ramp-up in cycle 1)
#   Azacitidine 75 mg/m2 SC d1-7 each cycle
#   6 cycles (168 days)
#   Reference: DiNardo CD et al. NEJM 2020 (VIALE-A)
# =============================================================================
cat("\n=== Scenario 2: VEN + AZA (VIALE-A) ===\n")

n_cycles     <- 6
cycle_length <- 28

# Venetoclax: ramp-up cycle 1, then 400mg
ven_times <- numeric(0)
ven_amts  <- numeric(0)

for (cyc in 0:(n_cycles - 1)) {
  offset <- cyc * cycle_length
  if (cyc == 0) {
    ven_times <- c(ven_times, offset + 0, offset + 1, offset + 2:27)
    ven_amts  <- c(ven_amts, 100, 200, rep(400, 26))
  } else {
    ven_times <- c(ven_times, offset + 0:27)
    ven_amts  <- c(ven_amts, rep(400, 28))
  }
}

e2_ven <- ev(time = ven_times, amt = ven_amts, cmt = "VEN_gut")

# AZA: 75 mg/m2 x 1.8 m2 = 135 mg SC d1-7 per cycle
aza_dose_per_inj <- 75 * 1.8  # 135 mg

aza_times <- numeric(0)
for (cyc in 0:(n_cycles - 1)) {
  offset <- cyc * cycle_length
  aza_times <- c(aza_times, offset + 0:6)
}

e2_aza <- ev(time = aza_times, amt = aza_dose_per_inj, cmt = "AZA_depot")

e2 <- e2_ven + e2_aza

p2 <- param(mod,
  use_VEN     = 1,
  use_AZA     = 1,
  use_GILT    = 0,
  use_ENASID  = 0,
  use_CYTARAB = 0,
  FLT3_status = 1,
  IDH2_status = 0,
  APL_status  = 0
)

sim2 <- mrgsim(p2, events = e2, end = 168, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "VEN + AZA (VIALE-A)")

cat("Day 28 CR:", get_time_point(sim2, 28)$CR_status, "\n")
cat("Day 168 MRD log10:", get_time_point(sim2, 168)$MRD_log %>% round(2), "\n")
cat("Day 168 MRD negative:", get_time_point(sim2, 168)$MRD_neg, "\n")

# --- Plot Scenario 2 ---
p2_cr <- ggplot(sim2, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#377EB8", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = seq(28, 168, by = 28), linetype = "dotted", color = "gray60") +
  labs(title = "BM Blast % — VEN+AZA", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p2_mrd <- ggplot(sim2, aes(x = time)) +
  geom_line(aes(y = MRD_log), color = "#377EB8", linewidth = 1.1) +
  geom_hline(yintercept = -3,   linetype = "dashed", color = "purple") +
  geom_hline(yintercept = MRD_sensitivity, linetype = "dotted", color = "gray50") +
  annotate("text", x = 100, y = -2.6, label = "MRD- threshold", size = 3, color = "purple") +
  labs(title = "MRD Trajectory", x = "Time (days)", y = "log10(Tumor Cells)") +
  theme_bw(base_size = 12)

p2_anc <- ggplot(sim2, aes(x = time)) +
  geom_ribbon(aes(ymin = 0, ymax = ANC_circ), fill = "#377EB8", alpha = 0.25) +
  geom_line(aes(y = ANC_circ), color = "#377EB8", linewidth = 0.9) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
  labs(title = "ANC over 6 Cycles", x = "Time (days)", y = "ANC (x10^9/L)") +
  theme_bw(base_size = 12)

fig2 <- (p2_cr | p2_mrd | p2_anc) +
  plot_annotation(
    title    = "Scenario 2: VEN + AZA (VIALE-A) — 6 Cycles",
    subtitle = "Venetoclax 400mg QD (ramp-up) + Azacitidine 75 mg/m2 SC d1-7",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig2)
ggsave("aml_scenario2_venaza.png", fig2, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 3: Gilteritinib 120mg QD (FLT3+ Relapsed/Refractory)
#   Continuous oral dosing for 6 months (180 days)
#   R/R patient: higher baseline tumor burden
#   Reference: Perl AE et al. NEJM 2019 (ADMIRAL trial)
# =============================================================================
cat("\n=== Scenario 3: Gilteritinib 120mg QD (FLT3+ R/R) ===\n")

gilt_times <- seq(0, 179, by = 1)  # once daily for 180 days

e3 <- ev(time = gilt_times, amt = 120, cmt = "GILT_gut")

# R/R patient: more aggressive disease (higher LSC, slower intrinsic clearance)
p3 <- param(mod,
  use_VEN      = 0,
  use_AZA      = 0,
  use_GILT     = 1,
  use_ENASID   = 0,
  use_CYTARAB  = 0,
  FLT3_status  = 1,
  IDH2_status  = 0,
  APL_status   = 0,
  k_LSC_prolif = 0.22,  # more aggressive in R/R
  k_LBC_clear  = 0.10   # slower clearance in R/R
)

mod3 <- mod %>%
  param(p3) %>%
  init(
    LSC       = 800,
    LPC       = 3500,
    LBC       = 8000,
    MRD_log   = log10(800 + 3500 + 8000),
    tumor_vol = 800 + 3500 + 8000
  )

sim3 <- mrgsim(mod3, events = e3, end = 180, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "Gilteritinib 120mg", CRi_status = 0)

cat("Day 56 BM blast%:", get_time_point(sim3, 56)$BM_blast_pct %>% round(1), "\n")
cat("Day 180 CR:", get_time_point(sim3, 180)$CR_status, "\n")
cat("GILT Css day 14:", get_time_point(sim3, 14)$GILT_conc %>% round(3), "ug/mL\n")

# --- Plot Scenario 3 ---
p3_blast <- ggplot(sim3, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#4DAF4A", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  labs(title = "BM Blast % — Gilteritinib", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p3_flt3 <- ggplot(sim3, aes(x = time)) +
  geom_line(aes(y = FLT3_inh * 100), color = "#4DAF4A", linewidth = 1) +
  geom_hline(yintercept = 85, linetype = "dashed", color = "navy", alpha = 0.7) +
  annotate("text", x = 130, y = 87, label = "Target >85%", size = 3, color = "navy") +
  labs(title = "FLT3 Pathway Inhibition (%)", x = "Time (days)", y = "FLT3 Inhibition (%)") +
  theme_bw(base_size = 12)

p3_conc <- ggplot(sim3, aes(x = time)) +
  geom_line(aes(y = GILT_conc), color = "#4DAF4A", linewidth = 1) +
  labs(title = "Gilteritinib Concentration", x = "Time (days)", y = "Gilteritinib (ug/mL)") +
  theme_bw(base_size = 12)

fig3 <- (p3_blast | p3_flt3 | p3_conc) +
  plot_annotation(
    title    = "Scenario 3: Gilteritinib 120mg QD — FLT3+ R/R AML (ADMIRAL)",
    subtitle = "Continuous daily dosing, 6-month simulation, R/R patient",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig3)
ggsave("aml_scenario3_gilteritinib.png", fig3, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 4: Enasidenib 100mg QD (IDH2+ Relapsed/Refractory)
#   Continuous oral dosing, 6 months
#   IDH2 inhibition restores normal differentiation (2-HG pathway)
#   Monitor for differentiation syndrome (early WBC-like rise)
#   Reference: Stein EM et al. Blood 2017; IDHENTIFY (Stein 2019 Lancet Oncol)
# =============================================================================
cat("\n=== Scenario 4: Enasidenib 100mg QD (IDH2+ R/R) ===\n")

enasid_times <- seq(0, 179, by = 1)

# Enasidenib PK/PD: differentiation induction channelled through AZA_diff_eff
# Override AZA channel with enasidenib PK parameters and adjusted IC50/Imax
p4 <- param(mod,
  use_VEN     = 0,
  use_AZA     = 1,   # AZA channel used as enasidenib differentiation proxy
  use_GILT    = 0,
  use_ENASID  = 1,
  use_CYTARAB = 0,
  FLT3_status = 0,
  IDH2_status = 1,
  APL_status  = 0,
  AZA_IC50    = 1.5,   # enasidenib: requires higher conc for full effect
  AZA_Imax    = 0.68,  # enasidenib differentiation induction
  AZA_CL      = 1.86,  # use enasidenib CL
  AZA_V1      = 55.8,
  AZA_V2      = 0.0,
  AZA_Q       = 0.0,
  AZA_ka      = 0.41,
  AZA_F       = 1.0
)

# Dose into AZA_depot as enasidenib surrogate (absorbed dose)
e4_enasid <- ev(time = enasid_times, amt = 100, cmt = "AZA_depot")

mod4 <- mod %>%
  param(p4) %>%
  init(
    LSC       = 700,
    LPC       = 3000,
    LBC       = 7000,
    MRD_log   = log10(700 + 3000 + 7000),
    tumor_vol = 700 + 3000 + 7000
  )

sim4 <- mrgsim(mod4, events = e4_enasid, end = 180, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(
    Scenario   = "Enasidenib 100mg",
    CRi_status = 0
  )

cat("Day 56 BM blast%:", get_time_point(sim4, 56)$BM_blast_pct %>% round(1), "\n")
cat("Day 180 CR:", get_time_point(sim4, 180)$CR_status, "\n")
cat("Peak ANC (diff syndrome proxy):", max(sim4$ANC_circ) %>% round(2), "x10^9/L\n")

# --- Plot Scenario 4 ---
p4_blast <- ggplot(sim4, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#984EA3", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  labs(title = "BM Blast % — Enasidenib", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p4_diff <- ggplot(sim4, aes(x = time)) +
  geom_line(aes(y = AZA_diff_eff * 100), color = "#984EA3", linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(title = "Differentiation Induction\n(IDH2 pathway restoration)",
       x = "Time (days)", y = "Differentiation Effect (%)") +
  theme_bw(base_size = 12)

p4_anc <- ggplot(sim4, aes(x = time)) +
  geom_line(aes(y = ANC_circ), color = "#984EA3", linewidth = 1) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "darkorange") +
  annotate("text", x = 90, y = 11, label = "Diff. syndrome threshold", size = 3, color = "darkorange") +
  labs(title = "ANC — Differentiation Syndrome Monitor",
       x = "Time (days)", y = "ANC (x10^9/L)") +
  theme_bw(base_size = 12)

fig4 <- (p4_blast | p4_diff | p4_anc) +
  plot_annotation(
    title    = "Scenario 4: Enasidenib 100mg QD — IDH2+ R/R AML (IDHENTIFY)",
    subtitle = "IDH2 inhibition restores differentiation; differentiation syndrome monitoring",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig4)
ggsave("aml_scenario4_enasidenib.png", fig4, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 5: VEN + AZA + Gilteritinib Triple (Investigational)
#   FLT3+ newly diagnosed AML, unfit for intensive chemotherapy
#   VEN 400mg QD d1-28 + AZA 75mg/m2 SC d1-7 + GILT 120mg QD continuous
#   Reference: Pratz KW et al. JCO 2022 (phase 1b/2 LACEWING-like)
# =============================================================================
cat("\n=== Scenario 5: VEN + AZA + Gilteritinib Triple (Investigational) ===\n")

gilt_times_168 <- seq(0, 167, by = 1)
e5_gilt <- ev(time = gilt_times_168, amt = 120, cmt = "GILT_gut")
e5 <- e2_ven + e2_aza + e5_gilt  # reuse VEN+AZA from Scenario 2

p5 <- param(mod,
  use_VEN     = 1,
  use_AZA     = 1,
  use_GILT    = 1,
  use_ENASID  = 0,
  use_CYTARAB = 0,
  FLT3_status = 1,
  IDH2_status = 0,
  APL_status  = 0
)

sim5 <- mrgsim(p5, events = e5, end = 168, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "VEN+AZA+GILT (Triple)")

cat("Day 28 CR:", get_time_point(sim5, 28)$CR_status, "\n")
cat("Day 168 MRD log10:", get_time_point(sim5, 168)$MRD_log %>% round(2), "\n")
cat("ANC nadir:", min(sim5$ANC_circ) %>% round(2), "\n")

# --- Plot Scenario 5 ---
p5_blast <- ggplot(sim5, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#FF7F00", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = seq(28, 168, by = 28), linetype = "dotted", color = "gray60") +
  labs(title = "BM Blast % — Triple Combo", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p5_mrd <- ggplot(sim5, aes(x = time)) +
  geom_line(aes(y = MRD_log), color = "#FF7F00", linewidth = 1.1) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "purple") +
  geom_hline(yintercept = MRD_sensitivity, linetype = "dotted", color = "gray50") +
  labs(title = "MRD Trajectory", x = "Time (days)", y = "log10(Tumor Cells)") +
  theme_bw(base_size = 12)

p5_tox <- ggplot(sim5, aes(x = time)) +
  geom_line(aes(y = ANC_circ, color = "ANC"), linewidth = 1) +
  geom_line(aes(y = PLT_circ / 50, color = "PLT/50"), linewidth = 1) +
  scale_color_manual(values = c("ANC" = "#FF7F00", "PLT/50" = "#A65628")) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "darkgreen") +
  labs(title = "Hematologic Toxicity", x = "Time (days)",
       y = "ANC (x10^9/L) / PLT/50", color = "Parameter") +
  theme_bw(base_size = 12)

fig5 <- (p5_blast | p5_mrd | p5_tox) +
  plot_annotation(
    title    = "Scenario 5: VEN + AZA + Gilteritinib Triple — FLT3+ ND AML",
    subtitle = "Investigational 3-drug combination, 6 cycles",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig5)
ggsave("aml_scenario5_triple.png", fig5, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 6: LDAC + VEN (Non-intensive, Elderly)
#   Low-dose cytarabine 20mg SC BID d1-10 + Venetoclax 600mg QD d1-28
#   (Pre-2020 dosing, VIALE-C trial)
#   Reference: Wei AH et al. JCO 2020 (VIALE-C); Pollyea DA Nat Med 2018
#   Target: elderly (>=75) unfit for intensive chemo
# =============================================================================
cat("\n=== Scenario 6: LDAC + VEN (Non-intensive, Elderly) ===\n")

# LDAC 20mg SC BID d1-10 per cycle, 6 cycles
ldac_times <- numeric(0)
for (cyc in 0:(n_cycles - 1)) {
  offset <- cyc * cycle_length
  for (dd in 0:9) {
    ldac_times <- c(ldac_times, offset + dd, offset + dd + 0.5)
  }
}
e6_ldac <- ev(time = ldac_times, amt = 20, cmt = "CYTARAB_cent")

# VEN 600mg QD d1-28 with ramp-up (100/200/400/600mg d1-4)
ven_600_times <- numeric(0)
ven_600_amts  <- numeric(0)
for (cyc in 0:(n_cycles - 1)) {
  offset <- cyc * cycle_length
  if (cyc == 0) {
    ven_600_times <- c(ven_600_times, offset + 0, offset + 1, offset + 2,
                       offset + 3:27)
    ven_600_amts  <- c(ven_600_amts, 100, 200, 400, rep(600, 25))
  } else {
    ven_600_times <- c(ven_600_times, offset + 0:27)
    ven_600_amts  <- c(ven_600_amts, rep(600, 28))
  }
}
e6_ven <- ev(time = ven_600_times, amt = ven_600_amts, cmt = "VEN_gut")

e6 <- e6_ldac + e6_ven

# Elderly patient: reduced CL (renal function decline)
p6 <- param(mod,
  use_VEN     = 1,
  use_AZA     = 0,
  use_GILT    = 0,
  use_ENASID  = 0,
  use_CYTARAB = 1,
  FLT3_status = 0,
  IDH2_status = 0,
  NPM1_status = 1,
  APL_status  = 0,
  VEN_CL      = 12.0,   # reduced CL in elderly
  CYTARAB_CL  = 35.0    # reduced renal function
)

sim6 <- mrgsim(p6, events = e6, end = 168, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "LDAC + VEN")

cat("Day 28 CR+CRi:", {
  d28 <- get_time_point(sim6, 28)
  d28$CR_status + d28$CRi_status
}, "\n")
cat("Day 168 MRD log10:", get_time_point(sim6, 168)$MRD_log %>% round(2), "\n")

# --- Plot Scenario 6 ---
p6_blast <- ggplot(sim6, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#A65628", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = seq(28, 168, by = 28), linetype = "dotted", color = "gray60") +
  labs(title = "BM Blast % — LDAC+VEN", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p6_ven <- ggplot(sim6, aes(x = time)) +
  geom_line(aes(y = VEN_conc), color = "#A65628", linewidth = 1) +
  labs(title = "Venetoclax Concentration (600mg, Elderly PK)",
       x = "Time (days)", y = "VEN (ug/mL)") +
  theme_bw(base_size = 12)

p6_anc <- ggplot(sim6, aes(x = time)) +
  geom_line(aes(y = ANC_circ), color = "#A65628", linewidth = 1) +
  geom_ribbon(aes(ymin = 0, ymax = pmin(ANC_circ, 1.0)),
              fill = "#A65628", alpha = 0.20) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "darkgreen") +
  labs(title = "ANC — Elderly Patient", x = "Time (days)", y = "ANC (x10^9/L)") +
  theme_bw(base_size = 12)

fig6 <- (p6_blast | p6_ven | p6_anc) +
  plot_annotation(
    title    = "Scenario 6: LDAC + Venetoclax 600mg QD — Elderly Non-intensive (VIALE-C)",
    subtitle = "LDAC 20mg SC BID d1-10 + VEN 600mg QD with ramp-up, 6 cycles",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig6)
ggsave("aml_scenario6_ldacven.png", fig6, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# SCENARIO 7: ATRA + ATO (APL — Acute Promyelocytic Leukemia)
#   ATRA 45 mg/m2/day PO divided q12h + ATO 0.15 mg/kg IV d1-5 weekly
#   Induction: 8 weeks; Consolidation: 4 x 28-day cycles
#   Reference: Lo-Coco F et al. NEJM 2013; Platzbecker U et al. JCO 2017
#   APL subtype: PML-RARa drives promyelocyte differentiation block
# =============================================================================
cat("\n=== Scenario 7: ATRA + ATO (APL) ===\n")

# ATRA: 45 mg/m2 x 1.8 m2 = 81 mg/day divided BID
# Pre-systemic bioavailability ~18%; absorbed dose per dose ~7.3 mg
atra_abs_per_dose <- 81 * 0.18 / 2   # ~7.3 mg absorbed per BID dose

# ATO: 0.15 mg/kg x 70 kg = 10.5 mg IV d1-5 weekly
ato_dose_per_day <- 0.15 * 70   # 10.5 mg

# Induction: 8 weeks (56 days)
atra_ind_times <- seq(0, 55, by = 0.5)           # BID for 56 days
ato_ind_weeks  <- 0:7                             # 8 weekly 5-day courses
ato_ind_times  <- unlist(lapply(ato_ind_weeks, function(w) w * 7 + 0:4))

e7_atra_ind <- ev(time = atra_ind_times, amt = atra_abs_per_dose, cmt = "ATRA_cent")
e7_ato_ind  <- ev(time = ato_ind_times,  amt = ato_dose_per_day,  cmt = "ATO_cent")

# Consolidation: 4 cycles x 28 days
consol_start <- 56
e7_atra_con_times <- numeric(0)
e7_ato_con_times  <- numeric(0)

for (cyc in 0:3) {
  offset_con <- consol_start + cyc * 28
  e7_atra_con_times <- c(e7_atra_con_times, seq(offset_con, offset_con + 13.5, by = 0.5))
  for (wk in 0:3) {
    e7_ato_con_times <- c(e7_ato_con_times, offset_con + wk * 7 + 0:4)
  }
}

e7_atra_con <- ev(time = e7_atra_con_times, amt = atra_abs_per_dose, cmt = "ATRA_cent")
e7_ato_con  <- ev(time = e7_ato_con_times,  amt = ato_dose_per_day,  cmt = "ATO_cent")

e7 <- e7_atra_ind + e7_ato_ind + e7_atra_con + e7_ato_con

# APL-specific parameterization
p7 <- param(mod,
  use_VEN     = 0,
  use_AZA     = 0,
  use_GILT    = 0,
  use_ENASID  = 0,
  use_CYTARAB = 0,
  FLT3_status = 0,
  IDH2_status = 0,
  NPM1_status = 0,
  APL_status  = 1,     # activate ATRA/ATO pathway
  ATRA_Emax   = 0.92,
  ATRA_EC50   = 0.12,
  ATO_Emax    = 0.90,
  ATO_EC50    = 0.06,
  k_LBC_apop  = 0.08,  # APL promyelocytes more susceptible
  k_LSC_apop  = 0.04
)

mod7 <- mod %>%
  param(p7) %>%
  init(
    LSC       = 300,
    LPC       = 1500,
    LBC       = 3500,  # APL: ~45% BM promyelocytes/blasts
    MRD_log   = log10(300 + 1500 + 3500),
    tumor_vol = 300 + 1500 + 3500
  )

sim7 <- mrgsim(mod7, events = e7, end = 168, delta = 0.5, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(Scenario = "ATRA + ATO (APL)", CRi_status = 0)

cat("Day 28 CR:", get_time_point(sim7, 28)$CR_status, "\n")
cat("Day 56 MRD log10:", get_time_point(sim7, 56)$MRD_log %>% round(2), "\n")
cat("Day 168 MRD negative:", get_time_point(sim7, 168)$MRD_neg, "\n")
cat("Note: ATRA auto-induction reduces plasma levels ~50% by day 14 in vivo\n")

# --- Plot Scenario 7 ---
p7_blast <- ggplot(sim7, aes(x = time)) +
  geom_line(aes(y = BM_blast_pct), color = "#F781BF", linewidth = 1.1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = 56, linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = 62, y = 55, label = "Consol start", size = 3, color = "blue") +
  labs(title = "BM Blast % — APL (ATRA+ATO)", x = "Time (days)", y = "BM Blast (%)") +
  theme_bw(base_size = 12)

p7_apop <- ggplot(sim7, aes(x = time)) +
  geom_line(aes(y = APL_apop * 100), color = "#F781BF", linewidth = 1) +
  labs(title = "ATRA+ATO Combined Apoptotic\nEffect (%)",
       x = "Time (days)", y = "Apoptosis Induction (%)") +
  theme_bw(base_size = 12)

p7_mrd <- ggplot(sim7, aes(x = time)) +
  geom_line(aes(y = MRD_log), color = "#F781BF", linewidth = 1.1) +
  geom_hline(yintercept = -3,           linetype = "dashed", color = "purple") +
  geom_hline(yintercept = MRD_sensitivity, linetype = "dotted", color = "gray50") +
  annotate("text", x = 115, y = -2.6, label = "MRD- threshold", size = 3, color = "purple") +
  labs(title = "MRD Trajectory (PML-RARa proxy)",
       x = "Time (days)", y = "log10(Tumor Cells)") +
  theme_bw(base_size = 12)

fig7 <- (p7_blast | p7_apop | p7_mrd) +
  plot_annotation(
    title    = "Scenario 7: ATRA + ATO — APL (PML-RARa+ AML M3 Subtype)",
    subtitle = "Induction 8 wk + Consolidation 4 cycles; ATRA 45mg/m2/d + ATO 0.15mg/kg d1-5/wk",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig7)
ggsave("aml_scenario7_aplatraato.png", fig7, width = 12, height = 4.5, dpi = 150)

# =============================================================================
# FINAL COMPARISON FIGURE — All 7 Scenarios
# =============================================================================
cat("\n=== Generating Final Comparison Figure ===\n")

# Build summary table
summary_table <- bind_rows(
  summarise_response(sim1, "7+3 Induction"),
  summarise_response(sim2, "VEN + AZA (VIALE-A)"),
  summarise_response(sim3, "Gilteritinib 120mg"),
  summarise_response(sim4, "Enasidenib 100mg"),
  summarise_response(sim5, "VEN+AZA+GILT (Triple)"),
  summarise_response(sim6, "LDAC + VEN"),
  summarise_response(sim7, "ATRA + ATO (APL)")
)

cat("\n--- Summary Table ---\n")
print(summary_table, n = Inf)

# Combine all simulations to 168 days for cross-scenario plots
all_sims <- bind_rows(
  sim1 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim2 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim3 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim4 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim5 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim6 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario),
  sim7 %>% select(time, BM_blast_pct, ANC_circ, PLT_circ, MRD_log, MRD_neg,
                  CR_status, CRi_status, Scenario)
) %>%
  mutate(
    Scenario_label = case_when(
      grepl("7\\+3",      Scenario) ~ "7+3 Induction",
      grepl("VIALE",      Scenario) ~ "VEN + AZA (VIALE-A)",
      grepl("Gilter",     Scenario) ~ "Gilteritinib 120mg",
      grepl("Enasidenib", Scenario) ~ "Enasidenib 100mg",
      grepl("Triple",     Scenario) ~ "VEN+AZA+GILT (Triple)",
      grepl("LDAC",       Scenario) ~ "LDAC + VEN",
      grepl("APL|ATRA",   Scenario) ~ "ATRA + ATO (APL)",
      TRUE ~ Scenario
    )
  ) %>%
  filter(time <= 168)

# ---- Panel A: BM Blast % trajectories ----
pA <- ggplot(all_sims, aes(x = time, y = BM_blast_pct, color = Scenario_label)) +
  geom_line(linewidth = 1.0, alpha = 0.85) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray30", linewidth = 0.6) +
  scale_color_manual(values = drug_cols) +
  labs(title = "A. Bone Marrow Blast % Over Time",
       x = "Time (days)", y = "BM Blast (%)", color = "Regimen") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.key.size  = unit(0.5, "cm"),
        legend.text      = element_text(size = 9))

# ---- Panel B: MRD trajectories ----
pB <- ggplot(all_sims, aes(x = time, y = MRD_log, color = Scenario_label)) +
  geom_line(linewidth = 1.0, alpha = 0.85) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "purple", linewidth = 0.6) +
  geom_hline(yintercept = MRD_sensitivity, linetype = "dotted", color = "gray50", linewidth = 0.5) +
  scale_color_manual(values = drug_cols) +
  annotate("text", x = 140, y = -2.55, label = "MRD- threshold", size = 3, color = "purple") +
  labs(title = "B. MRD (log10 Tumor Cells) Over Time",
       x = "Time (days)", y = "log10(Tumor Cells)", color = "Regimen") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ---- Panel C: ANC nadir profiles ----
pC <- ggplot(all_sims, aes(x = time, y = ANC_circ, color = Scenario_label)) +
  geom_line(linewidth = 0.9, alpha = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange", linewidth = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "darkgreen", linewidth = 0.5) +
  scale_color_manual(values = drug_cols) +
  labs(title = "C. ANC Profiles (Myelosuppression)",
       x = "Time (days)", y = "ANC (x10^9/L)", color = "Regimen") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ---- Panel D: Platelet trajectories ----
pD <- ggplot(all_sims, aes(x = time, y = PLT_circ, color = Scenario_label)) +
  geom_line(linewidth = 0.9, alpha = 0.75) +
  geom_hline(yintercept = 50,  linetype = "dashed", color = "red",       linewidth = 0.6) +
  geom_hline(yintercept = 100, linetype = "dotted", color = "darkgreen", linewidth = 0.5) +
  scale_color_manual(values = drug_cols) +
  labs(title = "D. Platelet Profiles",
       x = "Time (days)", y = "Platelets (x10^9/L)", color = "Regimen") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ---- Panel E: MRD negativity step-function ----
pE <- ggplot(all_sims, aes(x = time, y = MRD_neg, color = Scenario_label)) +
  geom_step(linewidth = 1.0, alpha = 0.85) +
  scale_color_manual(values = drug_cols) +
  scale_y_continuous(breaks = c(0, 1), labels = c("MRD+", "MRD-")) +
  labs(title = "E. MRD Negativity Over Time",
       x = "Time (days)", y = "MRD Status", color = "Regimen") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ---- Panel F: Summary bar chart at day 168 ----
final_summary <- all_sims %>%
  group_by(Scenario_label) %>%
  filter(abs(time - 168) == min(abs(time - 168))) %>%
  slice(1) %>%
  ungroup() %>%
  select(Scenario_label, CR_status, MRD_neg) %>%
  pivot_longer(cols = c(CR_status, MRD_neg),
               names_to = "Endpoint", values_to = "Value") %>%
  mutate(
    Endpoint = recode(Endpoint,
                      "CR_status" = "CR at Day 168",
                      "MRD_neg"   = "MRD- at Day 168")
  )

pF <- ggplot(final_summary, aes(x = Scenario_label, y = Value, fill = Scenario_label)) +
  geom_col(alpha = 0.85) +
  facet_wrap(~Endpoint, ncol = 2) +
  scale_fill_manual(values = drug_cols) +
  scale_y_continuous(breaks = c(0, 1), labels = c("No", "Yes")) +
  labs(title = "F. CR and MRD Negativity at Day 168",
       x = NULL, y = "Achieved (0/1)", fill = "Regimen") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 35, hjust = 1, size = 8)
  )

# ---- Assemble final comparison figure ----
fig_compare <- (pA + pB) / (pC + pD) / (pE + pF) +
  plot_annotation(
    title    = "AML QSP Model — Cross-Scenario Comparison (7 Treatment Regimens)",
    subtitle = paste0(
      "Newly diagnosed: 7+3, VEN+AZA, VEN+AZA+GILT, LDAC+VEN, ATRA+ATO(APL)\n",
      "Relapsed/Refractory: Gilteritinib (FLT3+), Enasidenib (IDH2+) | ",
      "Simulation period: 168 days"
    ),
    theme = theme(
      plot.title    = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray30")
    )
  )

print(fig_compare)
ggsave("aml_all_scenarios_comparison.png", fig_compare,
       width = 14, height = 16, dpi = 150)

cat("\n=== Simulation Complete ===\n")
cat("Output figures written:\n")
cat("  aml_scenario1_73induction.png\n")
cat("  aml_scenario2_venaza.png\n")
cat("  aml_scenario3_gilteritinib.png\n")
cat("  aml_scenario4_enasidenib.png\n")
cat("  aml_scenario5_triple.png\n")
cat("  aml_scenario6_ldacven.png\n")
cat("  aml_scenario7_aplatraato.png\n")
cat("  aml_all_scenarios_comparison.png\n")

# =============================================================================
# BONUS: BCL-2 EC50 Sensitivity Analysis
#   Vary VEN_EC50 over 8-fold range; examine BM blast% and MRD at day 56
#   for VEN+AZA regimen (Scenario 2 conditions)
# =============================================================================
cat("\n=== Bonus: BCL-2 EC50 Sensitivity (VEN+AZA) ===\n")

ec50_vals <- c(0.05, 0.10, 0.15, 0.25, 0.40)  # ug/mL

sens_results <- lapply(ec50_vals, function(ec50_val) {
  p_sens <- param(mod,
    use_VEN     = 1,
    use_AZA     = 1,
    use_GILT    = 0,
    use_ENASID  = 0,
    use_CYTARAB = 0,
    FLT3_status = 1,
    IDH2_status = 0,
    APL_status  = 0,
    VEN_EC50    = ec50_val
  )
  sim_s <- mrgsim(p_sens, events = e2, end = 56, delta = 0.5, obsonly = TRUE) %>%
    as.data.frame()
  d56 <- get_time_point(sim_s, 56)
  tibble(
    VEN_EC50    = ec50_val,
    Blast_d56   = round(d56$BM_blast_pct, 2),
    CR_d56      = d56$CR_status,
    MRD_log_d56 = round(d56$MRD_log, 2)
  )
})

sens_df <- bind_rows(sens_results)
cat("\nBCL-2 EC50 Sensitivity Table (VEN+AZA at day 56):\n")
print(sens_df)

p_sens <- ggplot(sens_df, aes(x = VEN_EC50, y = Blast_d56)) +
  geom_line(color = "#377EB8", linewidth = 1.2) +
  geom_point(color = "#377EB8", size = 3) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  scale_x_log10(breaks = ec50_vals, labels = ec50_vals) +
  labs(
    title    = "BCL-2 EC50 Sensitivity — BM Blast % at Day 56 (VEN+AZA)",
    subtitle = "Reference EC50 = 0.15 ug/mL (VIALE-A calibration)",
    x        = "VEN_EC50 (ug/mL, log scale)",
    y        = "BM Blast % at Day 56"
  ) +
  theme_bw(base_size = 12)

print(p_sens)
ggsave("aml_sensitivity_ec50.png", p_sens, width = 8, height = 5, dpi = 150)

cat("\nAll AML QSP simulations and sensitivity analysis completed successfully.\n")
