# =============================================================================
# Crohn's Disease (CD) — Quantitative Systems Pharmacology (QSP) Model
# mrgsolve R implementation
#
# 22-state ODE model
#   Drug PK   (10 states): IFX (2-cpt), ADA (SC+central), UST (2-cpt),
#                            VDZ (1-cpt), 6-TGN (AZA metabolite),
#                            Prednisone (1-cpt), Upadacitinib (1-cpt)
#   Disease PD (12 states): TNF, IL12/23, IL-17, Th17, Th1, Treg,
#                            Neutrophil, MucosalInflam, CRP, FecalCalprotectin,
#                            BMD, Hemoglobin
#
# Calibration references
#   IFX PK  : Fasanmade AA et al. Int J Clin Pharmacol Ther. 2009
#             Ordas I et al. Clin Gastroenterol Hepatol. 2012
#   ADA PK  : Ternant D et al. Br J Clin Pharmacol. 2015
#   UST PK  : Colombel JF et al. Gut. 2017; UNIFI pharmacokinetics
#   VDZ PK  : Rosario M et al. Inflamm Bowel Dis. 2015
#   Disease : Vodovotz Y & An G (eds.) Complex Systems and Computational
#             Biology Approaches to Acute Inflammation; ACCENT I, CHARM,
#             UNIFI, VARSITY, U-ACHIEVE trials
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ---------------------------------------------------------------------------
# MODEL CODE
# ---------------------------------------------------------------------------
cd_code <- '
$PROB Crohn Disease QSP Model — 22 compartment PK/PD

$PARAM
  // --- IFX 2-compartment IV PK ---
  Vc_IFX  = 4.0,    // L, central volume (population mean ~70 kg patient)
  Vp_IFX  = 3.5,    // L, peripheral volume
  CL_IFX  = 0.55,   // L/day, linear clearance
  Q_IFX   = 1.40,   // L/day, inter-compartmental CL

  // --- ADA SC 1-compartment PK ---
  V_ADA   = 7.8,    // L, distribution volume
  CL_ADA  = 0.35,   // L/day, clearance (t½ ~ 15 days)
  ka_ADA  = 0.38,   // day⁻¹, SC absorption rate constant
  F_ADA   = 0.64,   // SC bioavailability

  // --- UST 2-compartment PK (IV induction, SC maintenance) ---
  Vc_UST  = 3.5,    // L
  Vp_UST  = 7.3,    // L
  CL_UST  = 0.19,   // L/day  (t½ IV ~ 25 days)
  Q_UST   = 0.91,   // L/day
  ka_UST  = 0.17,   // day⁻¹ for SC maintenance
  F_UST   = 0.57,   // SC bioavailability

  // --- VDZ 1-compartment IV PK ---
  Vc_VDZ  = 3.3,    // L
  CL_VDZ  = 0.49,   // L/day (t½ ~ 26 days)

  // --- 6-TGN (AZA metabolite) PK ---
  V_TGN   = 1500,   // L, apparent volume (RBC-distributed)
  CL_TGN  = 10.5,   // L/day, effective clearance
  ka_TGN  = 0.040,  // day⁻¹, conversion from AZA to 6-TGN
  F_AZA   = 0.47,   // AZA oral bioavailability

  // --- Prednisone PK ---
  V_PRED  = 45.0,   // L
  CL_PRED = 7.2,    // L/day  (t½ ~ 4.4 h)
  ka_PRED = 12.0,   // day⁻¹ (rapid oral absorption)
  F_PRED  = 0.80,

  // --- Upadacitinib (JAK1i) PK ---
  V_UPA   = 85.0,   // L
  CL_UPA  = 24.0,   // L/day (t½ ~ 2.5 h, 45 mg QD)
  ka_UPA  = 2.5,    // day⁻¹
  F_UPA   = 0.79,

  // --- Anti-drug antibody (ATI for IFX) ---
  kATI    = 0.00012, // day⁻¹, ATI formation rate (immunogenicity)
  kATI_cl = 0.015,   // day⁻¹, ATI clearance
  EC50_ATI = 1.0,    // mcg/mL, IFX conc for 50% ATI formation suppression

  // ----------------------------------------------------------------
  // DISEASE PD PARAMETERS
  // ----------------------------------------------------------------

  // TNF-α (pg/mL)  — baseline 100 in active CD, ~15 in remission
  TNF0        = 100,
  kprod_TNF   = 5.00, // pg/mL/day, baseline prod (= kdeg*TNF0)
  kdeg_TNF    = 0.050, // day⁻¹
  // Anti-TNF Emax/EC50
  Emax_IFX_TNF = 0.93,  EC50_IFX_TNF = 1.0,  // mcg/mL
  Emax_ADA_TNF = 0.90,  EC50_ADA_TNF = 0.8,

  // IL-12/23 composite (pg/mL)  — baseline 80
  IL12230     = 80,
  kprod_IL1223 = 4.00,
  kdeg_IL1223  = 0.050,
  Emax_UST    = 0.88,  EC50_UST = 0.8,   // mcg/mL

  // IL-17A (pg/mL)  — baseline 60
  IL170       = 60,
  kprod_IL17  = 3.00,
  kdeg_IL17   = 0.050,
  kIL23_IL17  = 0.020, // IL-23 amplification of IL-17 production (normalized)

  // Th17 cells (normalized, SS = 1.0)
  Th170       = 1.0,
  kprol_Th17  = 0.060, // day⁻¹ (= kdeath at SS to maintain baseline)
  kdeg_Th17   = 0.060, // day⁻¹
  kIL23_Th17  = 0.030, // normalized IL-23 drive on Th17 proliferation
  kTreg_Th17  = 0.020, // Treg-mediated suppression

  // Th1 cells (normalized, SS = 1.0)
  Th10        = 1.0,
  kprol_Th1   = 0.060,
  kdeg_Th1    = 0.060,
  kIL12_Th1   = 0.025,

  // Treg cells (normalized, SS = 1.0)
  Treg0       = 1.0,
  kprol_Treg  = 0.050,
  kdeg_Treg   = 0.050,

  // Mucosal neutrophil (normalized, SS = 1.0)
  Neut0       = 1.0,
  kprod_Neut  = 0.050, // day⁻¹  (= kdeg at SS)
  kdeg_Neut   = 0.050,
  kIL17_Neut  = 0.030, // IL-17-driven neutrophil production
  kTNF_Neut   = 0.020,

  // Mucosal Inflammation Score 0–100 (SS = 65 for moderate active CD)
  MI0         = 65,
  kprod_MI    = 3.25,  // (= kdeg_MI * MI0 at SS)
  kdeg_MI     = 0.050,
  kTNF_MI     = 0.015, // TNF contribution (normalized to TNF0=100)
  kIL17_MI    = 0.010,
  kTreg_MI    = 0.010, // Treg anti-inflammatory effect
  kNeut_MI    = 0.008,

  // CRP (mg/L)  — SS = 25 in active CD
  CRP0        = 25,
  kprod_CRP   = 0.096, // day⁻¹  (kprod_CRP * MI0 = kdeg_CRP * CRP0)
  kdeg_CRP    = 0.250,

  // Fecal Calprotectin (μg/g)  — SS = 800
  FC0         = 800,
  kprod_FC    = 2.0,   // base production
  kdeg_FC     = 0.060,
  kNeut_FC    = 25.0,  // Neut-driven FC release
  kMI_FC      = 0.55,  // MI-driven FC

  // BMD (normalized T-score = 0 → 1.0; active CD → slow loss)
  BMD0        = 1.00,
  kform_BMD   = 0.000195, // day⁻¹ bone formation
  kresorb_BMD = 0.000200, // day⁻¹ base resorption
  kTNF_BMD    = 0.000080, // TNF-mediated extra resorption (normalized)
  kPRED_BMD   = 0.000150, // steroid-induced BMD loss (per mg/L PRED)

  // Hemoglobin (g/dL)  — SS = 10.5 (mild anemia of chronic disease)
  Hgb0        = 10.5,
  kprod_Hgb   = 0.0913, // day⁻¹ erythropoiesis
  kdeg_Hgb    = 0.0070,
  kMI_Hgb     = 0.00195, // MI-mediated erythropoiesis suppression

  // Combined drug PD effects on MI (Emax, EC50)
  Emax_VDZ_MI = 0.60,  EC50_VDZ_MI = 0.4,  // mcg/mL
  Emax_TGN_MI = 0.45,  EC50_TGN_MI = 235,  // pmol/8×10^8 RBC
  Emax_PRED   = 0.72,  EC50_PRED   = 0.05, // mg/L
  Emax_UPA    = 0.74,  EC50_UPA    = 0.08  // mg/L

$CMT
  // Drug PK compartments (amounts in mg)
  A1_IFX     // IFX central  (mg)
  A2_IFX     // IFX peripheral
  ADA_dep    // ADA SC depot
  A_ADA      // ADA central
  A1_UST     // UST central
  A2_UST     // UST peripheral
  A_VDZ      // VDZ central
  A_TGN      // 6-TGN (pmol-equivalent scaled)
  A_PRED     // Prednisone central
  A_UPA      // Upadacitinib central
  // Disease biology compartments
  TNF        // free TNF-α (pg/mL)
  IL1223     // IL-12/23 composite (pg/mL)
  IL17       // IL-17A (pg/mL)
  Th17       // Th17 cells (normalized)
  Th1        // Th1 cells  (normalized)
  Treg       // Regulatory T cells (normalized)
  Neut       // Mucosal neutrophils (normalized)
  MI         // Mucosal Inflammation Score (0–100)
  CRP        // C-reactive protein (mg/L)
  FC         // Fecal calprotectin (μg/g)
  BMD        // Bone mineral density (normalized)
  Hgb        // Hemoglobin (g/dL)

$INIT
  A1_IFX = 0, A2_IFX = 0,
  ADA_dep = 0, A_ADA = 0,
  A1_UST = 0, A2_UST = 0,
  A_VDZ = 0,
  A_TGN = 0,
  A_PRED = 0,
  A_UPA = 0,
  TNF    = 100,
  IL1223 = 80,
  IL17   = 60,
  Th17   = 1.0,
  Th1    = 1.0,
  Treg   = 1.0,
  Neut   = 1.0,
  MI     = 65,
  CRP    = 25,
  FC     = 800,
  BMD    = 1.00,
  Hgb    = 10.5

$ODE
  // ============================================================
  // DRUG CONCENTRATIONS (from amounts)
  // ============================================================
  double C1_IFX  = A1_IFX  / Vc_IFX;      // mcg/mL
  double C2_IFX  = A2_IFX  / Vp_IFX;
  double C_ADA   = A_ADA   / V_ADA;
  double C1_UST  = A1_UST  / Vc_UST;
  double C_VDZ   = A_VDZ   / Vc_VDZ;
  double C_TGN_c = A_TGN   / V_TGN;       // pmol-eq/L (scale to RBC units)
  double C_PRED  = A_PRED  / V_PRED;       // mg/L
  double C_UPA   = A_UPA   / V_UPA;        // mg/L

  // ============================================================
  // IFX 2-COMPARTMENT PK (amounts)
  // ============================================================
  double k10_IFX  = CL_IFX / Vc_IFX;
  double k12_IFX  = Q_IFX  / Vc_IFX;
  double k21_IFX  = Q_IFX  / Vp_IFX;
  // ATI increases clearance (maximum 3-fold at high ATI)
  double ATI_eff  = 0;  // placeholder; ATI tracked separately below as scalar
  dxdt_A1_IFX = -k10_IFX * A1_IFX - k12_IFX * A1_IFX + k21_IFX * A2_IFX;
  dxdt_A2_IFX =  k12_IFX * A1_IFX - k21_IFX * A2_IFX;

  // ============================================================
  // ADA SC DEPOT → CENTRAL
  // ============================================================
  dxdt_ADA_dep = -ka_ADA * ADA_dep;
  dxdt_A_ADA   =  ka_ADA * F_ADA * ADA_dep - (CL_ADA / V_ADA) * A_ADA;

  // ============================================================
  // UST 2-COMPARTMENT
  // ============================================================
  double k10_UST = CL_UST / Vc_UST;
  double k12_UST = Q_UST  / Vc_UST;
  double k21_UST = Q_UST  / Vp_UST;
  dxdt_A1_UST = -k10_UST * A1_UST - k12_UST * A1_UST + k21_UST * A2_UST;
  dxdt_A2_UST =  k12_UST * A1_UST - k21_UST * A2_UST;

  // ============================================================
  // VDZ 1-COMPARTMENT
  // ============================================================
  dxdt_A_VDZ = -(CL_VDZ / Vc_VDZ) * A_VDZ;

  // ============================================================
  // 6-TGN FROM AZA (single depot → RBC pool)
  // ============================================================
  dxdt_A_TGN = ka_TGN * F_AZA * (A_TGN > 0 ? A_TGN : 0) - (CL_TGN / V_TGN) * A_TGN;
  // NOTE: AZA dosing goes into A_TGN depot via events (cmt="A_TGN")

  // ============================================================
  // PREDNISONE ORAL
  // ============================================================
  dxdt_A_PRED = -ka_PRED * A_PRED - (CL_PRED / V_PRED) * A_PRED;
  // Oral dose added to depot → combined absorption+central (simplified 1-cpt)

  // ============================================================
  // UPADACITINIB ORAL
  // ============================================================
  dxdt_A_UPA = -ka_UPA * A_UPA - (CL_UPA / V_UPA) * A_UPA;

  // ============================================================
  // DRUG EFFECTS (Emax models, inhibitory)
  // ============================================================
  double E_IFX   = Emax_IFX_TNF * C1_IFX / (EC50_IFX_TNF + C1_IFX);
  double E_ADA   = Emax_ADA_TNF * C_ADA   / (EC50_ADA_TNF + C_ADA);
  // Combined anti-TNF (independent inhibition on log scale → additive on linear)
  double E_antiTNF = 1.0 - (1.0 - E_IFX) * (1.0 - E_ADA);
  double E_antiTNF_clamp = E_antiTNF > 0.97 ? 0.97 : E_antiTNF;

  double E_UST   = Emax_UST   * C1_UST / (EC50_UST   + C1_UST);
  double E_VDZ   = Emax_VDZ_MI * C_VDZ  / (EC50_VDZ_MI + C_VDZ);
  double C_TGN_rbc = C_TGN_c * 1000;  // convert to pmol/8×10^8 RBC units
  double E_TGN   = Emax_TGN_MI * C_TGN_rbc / (EC50_TGN_MI + C_TGN_rbc);
  double E_PRED  = Emax_PRED   * C_PRED     / (EC50_PRED   + C_PRED);
  double E_UPA   = Emax_UPA    * C_UPA      / (EC50_UPA    + C_UPA);

  // Combined anti-inflammatory effect on MI (independent inhibition)
  double E_combined = 1.0 - (1.0 - E_antiTNF_clamp) *
                             (1.0 - E_UST * 0.55) *
                             (1.0 - E_VDZ * 0.40) *
                             (1.0 - E_TGN * 0.30) *
                             (1.0 - E_PRED) *
                             (1.0 - E_UPA  * 0.85);

  // ============================================================
  // DISEASE BIOLOGY ODEs
  // ============================================================

  // --- TNF-α ---
  double stim_TNF  = 1.0 + 0.40 * (MI / MI0 - 1.0);  // MI amplifies TNF production
  if(stim_TNF < 0.1) stim_TNF = 0.1;
  dxdt_TNF = kprod_TNF * stim_TNF * (1.0 - E_antiTNF_clamp) - kdeg_TNF * TNF;

  // --- IL-12/23 composite ---
  double stim_IL1223 = 1.0 + 0.30 * (MI / MI0 - 1.0);
  if(stim_IL1223 < 0.1) stim_IL1223 = 0.1;
  dxdt_IL1223 = kprod_IL1223 * stim_IL1223 * (1.0 - E_UST) * (1.0 - E_UPA * 0.30)
                - kdeg_IL1223 * IL1223;

  // --- IL-17A ---
  double drive_IL17 = (1.0 + kIL23_IL17 * (IL1223 / IL12230 - 1.0));
  if(drive_IL17 < 0.1) drive_IL17 = 0.1;
  dxdt_IL17 = kprod_IL17 * drive_IL17 * Th17 / Th170
              * (1.0 - E_UPA * 0.55)
              - kdeg_IL17 * IL17;

  // --- Th17 cells ---
  double Th17_stim  = 1.0 + kIL23_Th17 * (IL1223 / IL12230 - 1.0);
  double Th17_supp  = kTreg_Th17 * (Treg / Treg0 - 1.0);
  if(Th17_stim < 0.05) Th17_stim = 0.05;
  dxdt_Th17 = kprol_Th17 * Th17 * Th17_stim * (1.0 - Th17_supp)
              * (1.0 - E_VDZ * 0.75)
              * (1.0 - E_TGN * 0.45)
              * (1.0 - E_UPA * 0.65)
              - kdeg_Th17 * Th17;

  // --- Th1 cells ---
  double Th1_stim  = 1.0 + kIL12_Th1 * (IL1223 / IL12230 - 1.0);
  if(Th1_stim < 0.05) Th1_stim = 0.05;
  dxdt_Th1 = kprol_Th1 * Th1 * Th1_stim
             * (1.0 - E_VDZ * 0.75)
             * (1.0 - E_TGN * 0.40)
             * (1.0 - E_UPA * 0.50)
             - kdeg_Th1 * Th1;

  // --- Regulatory T cells ---
  // Treg mildly expands with effective treatment (IL-10 loop)
  double Treg_stim = 1.0 + E_combined * 0.15;
  dxdt_Treg = kprol_Treg * Treg * Treg_stim - kdeg_Treg * Treg;

  // --- Mucosal neutrophils ---
  double Neut_drive = 1.0 + kIL17_Neut * (IL17 / IL170 - 1.0)
                           + kTNF_Neut  * (TNF  / TNF0  - 1.0);
  if(Neut_drive < 0.1) Neut_drive = 0.1;
  dxdt_Neut = kprod_Neut * Neut * Neut_drive * (1.0 - E_combined * 0.80)
              - kdeg_Neut * Neut;

  // --- Mucosal Inflammation Score (0–100) ---
  double MI_stim = kTNF_MI  * (TNF    / TNF0    - 1.0)
                 + kIL17_MI * (IL17   / IL170   - 1.0)
                 + kNeut_MI * (Neut   / Neut0   - 1.0);
  double MI_supp = kTreg_MI * (Treg   / Treg0   - 1.0);
  double MI_net  = 1.0 + MI_stim - MI_supp;
  if(MI_net < 0.05) MI_net = 0.05;
  dxdt_MI = kprod_MI * MI_net * (1.0 - E_combined) - kdeg_MI * MI;
  if(MI < 0)   MI = 0;
  if(MI > 100) MI = 100;

  // --- CRP (acute-phase response, reflects IL-6 which reflects MI) ---
  dxdt_CRP = kprod_CRP * MI - kdeg_CRP * CRP;
  if(CRP < 0) CRP = 0;

  // --- Fecal Calprotectin ---
  dxdt_FC = kprod_FC + kNeut_FC * Neut + kMI_FC * MI * (1.0 - E_combined * 0.90)
            - kdeg_FC * FC;
  if(FC < 0) FC = 0;

  // --- Bone Mineral Density ---
  double BMD_resorb = (kresorb_BMD + kTNF_BMD * TNF / TNF0
                       + kPRED_BMD * C_PRED) * BMD;
  dxdt_BMD = kform_BMD * BMD - BMD_resorb;

  // --- Hemoglobin ---
  double Hgb_suppress = kMI_Hgb * MI;
  double Hgb_prod     = kprod_Hgb * (1.0 - Hgb_suppress);
  if(Hgb_prod < 0) Hgb_prod = 0;
  dxdt_Hgb = Hgb_prod - kdeg_Hgb * Hgb;
  if(Hgb < 0) Hgb = 0;

$TABLE
  // Drug concentrations (mcg/mL)
  double C1_IFX_c  = A1_IFX / Vc_IFX;
  double C_ADA_c   = A_ADA  / V_ADA;
  double C1_UST_c  = A1_UST / Vc_UST;
  double C_VDZ_c   = A_VDZ  / Vc_VDZ;
  double C_TGN_r   = (A_TGN / V_TGN) * 1000; // pmol/8×10^8 RBC equivalent
  double C_PRED_c  = A_PRED / V_PRED;
  double C_UPA_c   = A_UPA  / V_UPA;

  // Clinical indices
  double CDAI      = 150 + (MI - MI0) * 2.2;
  if(CDAI < 0) CDAI = 0;

  double HBI       = 5.0 + (MI - MI0) * 0.13;
  if(HBI < 0) HBI = 0;

  double SES_CD    = (MI / 100.0) * 24.0;
  if(SES_CD < 0) SES_CD = 0;

  double Remission    = (CDAI < 150) ? 1.0 : 0.0;
  double Response100  = ((150 + (MI0 - MI0) * 2.2) - CDAI >= 100) ? 1.0 : 0.0;
  double MucHeal      = (SES_CD <= 2.0) ? 1.0 : 0.0;

  // α4β7 receptor occupancy (VDZ)
  double VDZ_RO    = C_VDZ_c / (0.7 + C_VDZ_c);  // Kd ~ 0.7 mcg/mL

  capture C_IFX    = C1_IFX_c;
  capture C_ADA    = C_ADA_c;
  capture C_UST    = C1_UST_c;
  capture C_VDZ    = C_VDZ_c;
  capture VDZ_RO   = VDZ_RO;
  capture C_TGN    = C_TGN_r;
  capture C_PRED   = C_PRED_c;
  capture C_UPA    = C_UPA_c;
  capture TNF_c    = TNF;
  capture IL1223_c = IL1223;
  capture IL17_c   = IL17;
  capture Th17_c   = Th17;
  capture Treg_c   = Treg;
  capture MI_c     = MI;
  capture CRP_c    = CRP;
  capture FC_c     = FC;
  capture BMD_c    = BMD;
  capture Hgb_c    = Hgb;
  capture CDAI_c   = CDAI;
  capture HBI_c    = HBI;
  capture SES_CD_c = SES_CD;
  capture Remission_c = Remission;
  capture MucHeal_c   = MucHeal;

$CAPTURE C_IFX C_ADA C_UST C_VDZ VDZ_RO C_TGN C_PRED C_UPA
         TNF_c IL1223_c IL17_c Th17_c Treg_c
         MI_c CRP_c FC_c BMD_c Hgb_c
         CDAI_c HBI_c SES_CD_c Remission_c MucHeal_c
'

# ---------------------------------------------------------------------------
# COMPILE MODEL
# ---------------------------------------------------------------------------
mod <- mcode("crohn_qsp", cd_code)

# ===========================================================================
# DOSING HELPERS
# ===========================================================================

# IFX induction + maintenance: 5 mg/kg IV at w0, w2, w6, then q8w
# amount in mg = dose_mgkg * weight (70 kg default)
make_IFX_events <- function(weeks = 52, dose_mgkg = 5, wt = 70) {
  amt <- dose_mgkg * wt  # mg
  # Induction: w0, w2, w6
  ind  <- c(0, 14, 42)
  # Maintenance q8w from w14 to end
  maint <- seq(98, weeks * 7, by = 56)
  times <- c(ind, maint)
  ev(amt = amt, cmt = "A1_IFX", time = times, rate = -2)  # rate=-2 => bolus
}

# ADA: 160/80/40 mg → 40 mg q2w
make_ADA_events <- function(weeks = 52) {
  times <- c(0, 14, 28, seq(42, weeks * 7, by = 14))
  amts  <- c(160, 80, 40, rep(40, length(seq(42, weeks * 7, by = 14))))
  ev(amt = amts, cmt = "ADA_dep", time = times)
}

# UST: IV induction (6 mg/kg → ~420 mg for 70 kg) at w0, then SC 90 mg q8w
make_UST_events <- function(weeks = 52, wt = 70) {
  iv_dose  <- round(6 * wt / 10) * 10  # rounded to nearest 10
  sc_times <- seq(56, weeks * 7, by = 56)
  c(ev(amt = iv_dose, cmt = "A1_UST", time = 0, rate = -2),
    ev(amt = 90,      cmt = "A1_UST", time = sc_times))
}

# VDZ: 300 mg IV at w0, w2, w6, then q8w
make_VDZ_events <- function(weeks = 52) {
  ind   <- c(0, 14, 42)
  maint <- seq(98, weeks * 7, by = 56)
  ev(amt = 300, cmt = "A_VDZ", time = c(ind, maint), rate = -2)
}

# AZA 150 mg/day (continuous; modeled as daily bolus into TGN depot)
make_AZA_events <- function(weeks = 52, dose_mg = 150) {
  ev(amt = dose_mg, cmt = "A_TGN", time = 0:(weeks * 7), ii = 1, addl = 0)
}

# Prednisone 40 mg QD for 4 weeks → 30→20→10→off (4-week taper)
make_PRED_events <- function() {
  ev(amt = c(rep(40, 28), rep(30, 14), rep(20, 14), rep(10, 7)),
     cmt = "A_PRED",
     time = 0:62)
}

# Upadacitinib 45 mg QD (induction) → 30 mg QD (maintenance)
make_UPA_events <- function(weeks = 52) {
  ind_times  <- 0:55
  maint_times <- 56:(weeks * 7)
  c(ev(amt = 45, cmt = "A_UPA", time = ind_times),
    ev(amt = 30, cmt = "A_UPA", time = maint_times))
}

# ===========================================================================
# SIMULATION SCENARIOS
# ===========================================================================
sim_years   <- 1.5
sim_end_day <- sim_years * 365
obs_times   <- c(seq(0, 90, by = 7), seq(91, sim_end_day, by = 14))

run_scenario <- function(events, label, mod = mod, end = sim_end_day) {
  mrgsim(mod, events, end = end, delta = 1, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(Scenario = label)
}

# Scenario 1: No Treatment (disease course, control)
ev_none <- ev(time = 0, amt = 0, cmt = "A1_IFX")  # null event

# Scenario 2: Infliximab 5 mg/kg induction + maintenance
ev_ifx <- make_IFX_events(weeks = sim_years * 52)

# Scenario 3: Adalimumab 160/80/40 mg
ev_ada <- make_ADA_events(weeks = sim_years * 52)

# Scenario 4: Ustekinumab
ev_ust <- make_UST_events(weeks = sim_years * 52)

# Scenario 5: Vedolizumab
ev_vdz <- make_VDZ_events(weeks = sim_years * 52)

# Scenario 6: Upadacitinib (JAK1 inhibitor)
ev_upa <- make_UPA_events(weeks = sim_years * 52)

# Scenario 7: IFX + AZA combination (combination therapy)
ev_ifx_aza <- as.ev(bind_rows(as.data.frame(ev_ifx),
                               as.data.frame(make_AZA_events(weeks = sim_years * 52))))

cat("Running simulations...\n")
scenarios <- list(
  "No Treatment"      = ev_none,
  "Infliximab"        = ev_ifx,
  "Adalimumab"        = ev_ada,
  "Ustekinumab"       = ev_ust,
  "Vedolizumab"       = ev_vdz,
  "Upadacitinib"      = ev_upa,
  "IFX + Azathioprine" = ev_ifx_aza
)

results <- lapply(names(scenarios), function(nm) {
  cat("  Simulating:", nm, "\n")
  mrgsim(mod, scenarios[[nm]], end = sim_end_day, delta = 1, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(Scenario = nm)
}) %>% bind_rows()

cat("Simulations complete.\n")

# ===========================================================================
# PLOT RESULTS
# ===========================================================================
scenario_colors <- c(
  "No Treatment"       = "#757575",
  "Infliximab"         = "#E53935",
  "Adalimumab"         = "#8E24AA",
  "Ustekinumab"        = "#1E88E5",
  "Vedolizumab"        = "#00897B",
  "Upadacitinib"       = "#F4511E",
  "IFX + Azathioprine" = "#3949AB"
)

results$time_weeks <- results$time / 7

# Helper plot function
plot_var <- function(data, yvar, ylab, title, pct = FALSE) {
  ggplot(data, aes(x = time_weeks, y = .data[[yvar]], color = Scenario)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_color_manual(values = scenario_colors) +
    labs(x = "Time (weeks)", y = ylab, title = title) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(size = 11, face = "bold"))
}

p1  <- plot_var(results, "CDAI_c", "CDAI", "Crohn's Disease Activity Index") +
       geom_hline(yintercept = 150, linetype = "dashed", color = "gray40") +
       annotate("text", x = max(results$time_weeks)*0.85, y = 155, label = "Remission", size = 3)

p2  <- plot_var(results, "MI_c",   "MI Score (0–100)", "Mucosal Inflammation Score")
p3  <- plot_var(results, "CRP_c",  "CRP (mg/L)",       "C-Reactive Protein") +
       geom_hline(yintercept = 5, linetype = "dashed", color = "gray40")
p4  <- plot_var(results, "FC_c",   "FC (μg/g)",        "Fecal Calprotectin") +
       geom_hline(yintercept = 150, linetype = "dashed", color = "gray40")
p5  <- plot_var(results, "TNF_c",  "TNF-α (pg/mL)",    "Mucosal TNF-α")
p6  <- plot_var(results, "IL17_c", "IL-17A (pg/mL)",   "Mucosal IL-17A")
p7  <- plot_var(results, "Hgb_c",  "Hemoglobin (g/dL)","Hemoglobin")
p8  <- plot_var(results, "BMD_c",  "BMD (normalized)", "Bone Mineral Density") +
       geom_hline(yintercept = 0.97, linetype = "dashed", color = "gray40")

# Drug PK comparison (anti-TNF subset)
drug_pk <- results %>% filter(Scenario %in% c("Infliximab", "Adalimumab", "IFX + Azathioprine"))
p9  <- ggplot(drug_pk, aes(x = time_weeks, y = C_IFX, color = Scenario)) +
       geom_line(linewidth = 0.8, alpha = 0.9) +
       scale_color_manual(values = scenario_colors) +
       geom_hline(yintercept = 3, linetype = "dashed", color = "gray40") +
       annotate("text", x = 5, y = 3.4, label = "IFX trough target (≥3 μg/mL)", size = 3) +
       labs(x = "Time (weeks)", y = "IFX conc. (mcg/mL)", title = "Infliximab PK") +
       theme_bw(base_size = 11) + theme(legend.position = "bottom", legend.title = element_blank())

ada_pk <- results %>% filter(Scenario %in% c("Adalimumab"))
p10 <- ggplot(ada_pk, aes(x = time_weeks, y = C_ADA, color = Scenario)) +
       geom_line(linewidth = 0.8, color = "#8E24AA") +
       geom_hline(yintercept = 5, linetype = "dashed", color = "gray40") +
       labs(x = "Time (weeks)", y = "ADA conc. (mcg/mL)", title = "Adalimumab PK") +
       theme_bw(base_size = 11) + theme(legend.position = "none")

# Composite figure
fig_main <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8) +
  plot_annotation(
    title   = "Crohn's Disease QSP Model — Treatment Scenarios",
    subtitle = "Infliximab · Adalimumab · Ustekinumab · Vedolizumab · Upadacitinib · IFX+AZA",
    theme   = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(fig_main)

# ===========================================================================
# SUMMARY TABLE AT WEEK 52
# ===========================================================================
week52 <- results %>%
  filter(abs(time - 364) <= 1) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  select(Scenario, CDAI_c, MI_c, CRP_c, FC_c, TNF_c, Hgb_c, BMD_c,
         Remission_c, MucHeal_c) %>%
  rename(
    CDAI     = CDAI_c,
    MI       = MI_c,
    CRP      = CRP_c,
    FC       = FC_c,
    TNF      = TNF_c,
    Hgb      = Hgb_c,
    BMD      = BMD_c,
    Remission = Remission_c,
    Mucosal_Healing = MucHeal_c
  ) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

cat("\n========== Week 52 Summary ==========\n")
print(as.data.frame(week52))

# ===========================================================================
# POPULATION VARIABILITY (IFX scenario, 200 patients)
# ===========================================================================
# Introduce PK variability: CL_IFX CV=35%, Vc_IFX CV=25%
set.seed(42)
n_pat <- 200

idata_pop <- tibble(
  ID      = 1:n_pat,
  CL_IFX  = rlnorm(n_pat, log(0.55), 0.35),
  Vc_IFX  = rlnorm(n_pat, log(4.0),  0.25),
  Vp_IFX  = rlnorm(n_pat, log(3.5),  0.25),
  CL_ADA  = rlnorm(n_pat, log(0.35), 0.35),
  V_ADA   = rlnorm(n_pat, log(7.8),  0.25)
)

cat("\nRunning population simulation (N=200, IFX)...\n")
pop_sim <- mrgsim(mod, ev_ifx,
                  idata = idata_pop,
                  end   = 365, delta = 7,
                  obsonly = TRUE) %>%
  as.data.frame()

# 5th, 50th, 95th percentile CRP over time
pop_crp <- pop_sim %>%
  group_by(time) %>%
  summarise(
    p05 = quantile(CRP_c,  0.05),
    p50 = quantile(CRP_c,  0.50),
    p95 = quantile(CRP_c,  0.95),
    .groups = "drop"
  ) %>%
  mutate(time_weeks = time / 7)

p_pop <- ggplot(pop_crp, aes(x = time_weeks)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "#E53935", alpha = 0.25) +
  geom_line(aes(y = p50), color = "#E53935", linewidth = 1.2) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  labs(x = "Time (weeks)", y = "CRP (mg/L)",
       title = "IFX Population PK/PD Simulation (N=200)",
       subtitle = "Median (solid), 5th–95th percentile (ribbon)") +
  theme_bw(base_size = 12)

print(p_pop)

cat("\n=== Crohn's Disease QSP Model simulation complete ===\n")
