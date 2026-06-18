################################################################################
# Chronic Kidney Disease (CKD) — Quantitative Systems Pharmacology Model
# mrgsolve ODE-based PK/PD Model
#
# Disease Subsystems:
#   1. Kidney pathophysiology (nephron loss, eGFR, proteinuria)
#   2. RAAS cascade (Ang II, Aldosterone, MR activation)
#   3. Inflammation (NF-κB, IL-6, TNF-α, macrophage infiltration)
#   4. Fibrosis (TGF-β1, Smad2/3, collagen deposition)
#   5. CKD-Mineral Bone Disease (FGF-23, Klotho, PTH, Vitamin D)
#   6. Anemia (EPO, HIF-PHI, hepcidin, hemoglobin)
#   7. Cardiovascular (endothelial dysfunction, LVH)
#
# Drug PK: ACE inhibitors, ARBs, finerenone, SGLT2 inhibitors,
#          ESAs (epoetin/darbepoetin), HIF-PHI (roxadustat),
#          phosphate binders, calcimimetics
#
# Key References:
#   - Heerspink HJL et al. DAPA-CKD. NEJM 2020;383:1436-1446.
#   - Bakris GL et al. FIDELIO-DKD. NEJM 2020;383:2219-2229.
#   - Perkovic V et al. CREDENCE. NEJM 2019;380:2295-2306.
#   - Provenzano R et al. DOLOMITES. JASN 2021;32:2021-2032.
#   - Remuzzi G et al. Lancet 2002;359:1309-1315.
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────
# MODEL DEFINITION
# ─────────────────────────────────────────────────────────────
ckd_model_code <- '
$PROB
CKD QSP Model — RAAS, Inflammation, Fibrosis, CKD-MBD, Anemia, CV

$PARAM @annotated
// ── Patient Characteristics ──────────────────────────────
eGFR0       : 45     : Baseline eGFR (mL/min/1.73m2)
UACR0       : 300    : Baseline UACR (mg/g Cr)
Hb0         : 11.0   : Baseline hemoglobin (g/dL)
SBP0        : 150    : Baseline systolic BP (mmHg)
Age         : 65     : Patient age (yr)
T2DM_flag   : 1      : Type 2 diabetes (0=no, 1=yes)

// ── Disease Progression Parameters ──────────────────────────
kGS         : 0.0003 : Rate of glomerulosclerosis progression (/day)
kTA         : 0.0002 : Rate of tubular atrophy (/day)
kNephron    : 0.0004 : Nephron loss rate constant (/day)
kFib        : 0.0008 : Fibrosis rate constant (/day)
kEGFR_loss  : 0.0006 : eGFR decline rate constant (/day)
kPU_eGFR    : 0.01   : Proteinuria → eGFR injury coupling

// ── RAAS Parameters ───────────────────────────────────────
AngII_ss    : 1.0    : Angiotensin II at steady-state (normalized)
Aldo_ss     : 1.0    : Aldosterone at steady-state (normalized)
kAngII_prod : 0.5    : Ang II production rate constant (/h)
kAngII_deg  : 0.5    : Ang II degradation rate constant (/h)
kAldo_prod  : 0.3    : Aldosterone production rate from Ang II
kAldo_deg   : 0.3    : Aldosterone degradation rate (/h)

// ── Inflammation Parameters ───────────────────────────────
IL6_ss      : 1.0    : Baseline IL-6 (normalized)
TNFa_ss     : 1.0    : Baseline TNF-α (normalized)
kIL6_prod   : 0.1    : IL-6 production rate constant (/h)
kIL6_deg    : 0.15   : IL-6 degradation rate constant (/h)
kTNFa_prod  : 0.08   : TNF-α production rate constant (/h)
kTNFa_deg   : 0.12   : TNF-α degradation rate constant (/h)
kMacro      : 0.05   : Macrophage infiltration rate (/h)
kMacro_out  : 0.08   : Macrophage clearance rate (/h)
NF_kB_base  : 1.0    : Baseline NF-κB activity

// ── TGF-β/Fibrosis Parameters ────────────────────────────
TGFb1_ss    : 1.0    : Baseline TGF-β1 (normalized)
kTGF_prod   : 0.06   : TGF-β1 production rate (/h)
kTGF_deg    : 0.08   : TGF-β1 degradation rate (/h)
kCol_prod   : 0.04   : Collagen deposition rate
kCol_deg    : 0.01   : Collagen degradation rate (/h)
Smad7_inh   : 0.3    : Smad7 inhibitory coefficient
kEMT        : 0.02   : EMT rate constant

// ── CKD-MBD Parameters ───────────────────────────────────
FGF23_ss    : 100    : Baseline FGF-23 (pg/mL)
PTH_ss      : 65     : Baseline PTH (pg/mL)
Klotho_ss   : 1.0    : Baseline Klotho (normalized)
Phos_ss     : 3.8    : Baseline serum phosphate (mg/dL)
Ca_ss       : 9.2    : Baseline serum calcium (mg/dL)
kFGF23_prod : 0.08   : FGF-23 production rate
kFGF23_deg  : 0.08   : FGF-23 degradation rate (/h)
kPTH_prod   : 0.12   : PTH production rate
kPTH_deg    : 0.15   : PTH degradation rate (/h)
kKlotho_deg : 0.05   : Klotho degradation rate
kKlotho_syn : 0.05   : Klotho synthesis rate
kPhos_retent: 0.02   : Phosphate retention rate with eGFR loss
kVitD_prod  : 0.04   : Active Vitamin D production rate
kVitD_deg   : 0.06   : Active Vitamin D degradation (/h)

// ── Anemia Parameters ────────────────────────────────────
EPO_ss      : 12     : Baseline EPO (mIU/mL)
Hep_ss      : 25     : Baseline Hepcidin (ng/mL)
kEPO_prod   : 0.10   : EPO production rate
kEPO_deg    : 0.30   : EPO degradation (/h)
kHb_prod    : 0.15   : Hemoglobin production rate (g/dL/day)
kHb_deg     : 0.008  : Hemoglobin degradation rate (/day)
kHep_prod   : 0.08   : Hepcidin production rate (IL-6 driven)
kHep_deg    : 0.05   : Hepcidin degradation (/h)

// ── Cardiovascular Parameters ─────────────────────────────
LVH_ss      : 1.0    : Baseline LVH index (normalized)
kLVH_prod   : 0.001  : LVH progression rate (/day)
kLVH_reg    : 0.0005 : LVH regression rate (/day)
kVC         : 0.002  : Vascular calcification rate (/day)

// ── Drug PK Parameters (ACE inhibitor — Ramipril) ────────
F_ACEi      : 0.28   : Bioavailability ramipril
ka_ACEi     : 0.8    : Absorption rate constant (/h)
CL_ACEi     : 2.5    : Clearance ramiprilat (L/h)
V_ACEi      : 12     : Volume of distribution (L)
IC50_ACEi   : 0.5    : IC50 for ACE inhibition (normalized)
Imax_ACEi   : 0.90   : Maximum ACE inhibition

// ── Drug PK Parameters (ARB — Losartan/EXP-3174) ─────────
F_ARB       : 0.33   : Bioavailability losartan
ka_ARB      : 0.9    : Absorption rate (/h)
CL_ARB      : 7.0    : Clearance EXP-3174 (L/h)
V_ARB       : 34     : Volume of distribution (L)
IC50_ARB    : 0.3    : IC50 AT1R blockade (normalized)
Imax_ARB    : 0.85   : Maximum AT1R blockade

// ── Drug PK Parameters (Finerenone — nsMRA) ──────────────
F_Fin       : 0.44   : Bioavailability finerenone
ka_Fin      : 1.2    : Absorption rate (/h)
CL_Fin      : 5.0    : Total clearance (L/h)
V1_Fin      : 22     : Central volume (L)
V2_Fin      : 18     : Peripheral volume (L)
Q_Fin       : 3.0    : Intercompartmental clearance (L/h)
IC50_Fin    : 0.8    : IC50 MR blockade (normalized)
Imax_Fin    : 0.90   : Maximum MR blockade

// ── Drug PK Parameters (SGLT2i — Dapagliflozin) ──────────
F_SGLT2     : 0.78   : Bioavailability dapagliflozin
ka_SGLT2    : 1.5    : Absorption rate (/h)
CL_SGLT2    : 7.2    : Total clearance (L/h)
V_SGLT2     : 22     : Volume of distribution (L)
IC50_SGLT2  : 0.6    : IC50 SGLT2 inhibition (normalized)
Imax_SGLT2  : 0.88   : Maximum SGLT2 inhibition

// ── Drug PK Parameters (ESA — Epoetin alfa) ──────────────
F_ESA       : 0.62   : Bioavailability epoetin SC
ka_ESA      : 0.016  : Absorption rate SC (/h)
CL_ESA      : 0.12   : Total clearance (L/h)
V_ESA       : 8.0    : Volume of distribution (L)
Emax_ESA    : 3.0    : Maximum EPO stimulation factor
EC50_ESA    : 10     : EC50 epoetin (mIU/mL)

// ── Drug PK Parameters (HIF-PHI — Roxadustat) ────────────
F_PHI       : 0.88   : Bioavailability roxadustat
ka_PHI      : 1.1    : Absorption rate (/h)
CL_PHI      : 3.5    : Total clearance (L/h)
V_PHI       : 16     : Volume of distribution (L)
Imax_PHI    : 0.80   : Maximum PHD inhibition
IC50_PHI    : 1.5    : IC50 PHD inhibition (normalized)

// ── Drug Doses (default = 0, set via event table) ────────
DOSE_ACEi   : 0      : Ramipril dose (mg)
DOSE_ARB    : 0      : Losartan dose (mg)
DOSE_Fin    : 0      : Finerenone dose (mg)
DOSE_SGLT2  : 0      : Dapagliflozin dose (mg)
DOSE_ESA    : 0      : Epoetin alfa dose (IU)
DOSE_PHI    : 0      : Roxadustat dose (mg)

$CMT @annotated
// Drug PK compartments
ACEi_gut   : Ramipril gut compartment (mg)
ACEi_c     : Ramiprilat central (mg/L equiv)
ARB_gut    : Losartan gut compartment (mg)
ARB_c      : EXP-3174 central (mg/L)
Fin_gut    : Finerenone gut (mg)
Fin_c      : Finerenone central (mg/L)
Fin_p      : Finerenone peripheral (mg/L)
SGLT2_gut  : Dapagliflozin gut (mg)
SGLT2_c    : Dapagliflozin central (mg/L)
ESA_sc     : Epoetin SC depot (IU)
ESA_c      : Epoetin plasma (mIU/mL equiv)
PHI_gut    : Roxadustat gut (mg)
PHI_c      : Roxadustat central (mg/L)

// Disease state compartments
Nephron    : Remaining nephrons (fraction 0–1)
eGFR       : Estimated GFR (mL/min/1.73m2)
UACR_st    : UACR state (mg/g)
AngII      : Angiotensin II (normalized)
Aldo       : Aldosterone (normalized)
IL6        : IL-6 (normalized)
TNFa       : TNF-alpha (normalized)
Macro      : Infiltrating macrophages (normalized)
TGFb       : TGF-β1 (normalized)
Collagen   : Collagen/fibrosis index (0–1)
FGF23      : FGF-23 (pg/mL)
Klotho     : Alpha-Klotho (normalized)
Phos       : Serum phosphate (mg/dL)
PTH        : Intact PTH (pg/mL)
VitD       : Active Vitamin D / Calcitriol (pg/mL)
EPO        : Plasma EPO (mIU/mL)
Hepcidin   : Hepcidin (ng/mL)
Hemoglobin : Hemoglobin (g/dL)
LVH_idx    : LVH index (normalized)
VC_idx     : Vascular calcification index (0–1)
BP         : Systolic blood pressure (mmHg)

$MAIN
// Normalize Ang II and Aldo at baseline
double AngII_norm = AngII / AngII_ss;
double Aldo_norm  = Aldo  / Aldo_ss;
double IL6_norm   = IL6   / IL6_ss;
double TGFb_norm  = TGFb  / TGFb1_ss;
double FGF23_norm = FGF23 / FGF23_ss;
double PTH_norm   = PTH   / PTH_ss;

// ── Drug concentrations ──────────────────────────────────────
double C_ACEi  = ACEi_c  / V_ACEi;
double C_ARB   = ARB_c   / V_ARB;
double C_Fin   = Fin_c   / V1_Fin;
double C_SGLT2 = SGLT2_c / V_SGLT2;
double C_ESA   = ESA_c;
double C_PHI   = PHI_c   / V_PHI;

// ── Drug effects (Emax/Imax Hill functions) ──────────────────
double E_ACEi  = Imax_ACEi  * C_ACEi  / (IC50_ACEi  + C_ACEi);
double E_ARB   = Imax_ARB   * C_ARB   / (IC50_ARB   + C_ARB);
double E_Fin   = Imax_Fin   * C_Fin   / (IC50_Fin   + C_Fin);
double E_SGLT2 = Imax_SGLT2 * C_SGLT2 / (IC50_SGLT2 + C_SGLT2);
double E_ESA   = Emax_ESA   * C_ESA   / (EC50_ESA   + C_ESA);
double E_PHI   = Imax_PHI   * C_PHI   / (IC50_PHI   + C_PHI);

// Combined RAS blockade (ACEi OR ARB)
double E_RAS   = 1.0 - (1.0 - E_ACEi) * (1.0 - E_ARB);

// ── Glomerular filtration pressure effect ────────────────────
double GFP_effect = 1.0 + 0.5 * AngII_norm * (1.0 - E_RAS);

// SGLT2i reduces glomerular hypertension via TGF
double SGLT2_GH_red = 0.3 * E_SGLT2;

// ── NF-κB activation (from AngII, Aldo, uremia) ─────────────
double Uremia_tox  = (1.0 - eGFR / eGFR0);
double NF_kB_act   = NF_kB_base * (1.0 + 0.5 * AngII_norm + 0.3 * Aldo_norm
                                   + 0.4 * Uremia_tox)
                                * (1.0 - 0.5 * E_RAS) * (1.0 - 0.3 * E_Fin);

// ── HIF-PHI effect on EPO ────────────────────────────────────
double HIF_stim = 1.0 + 2.5 * E_PHI;

// ── Smad2/3 activity (inhibited by Smad7 ~30%) ───────────────
double Smad23_act = TGFb_norm * (1.0 - Smad7_inh) * (1.0 - 0.4 * E_Fin);

$ODE
// ── DRUG PK ODEs ─────────────────────────────────────────────
// ACE inhibitor (Ramipril → Ramiprilat)
dxdt_ACEi_gut = -ka_ACEi * ACEi_gut;
dxdt_ACEi_c   = F_ACEi * ka_ACEi * ACEi_gut - CL_ACEi * C_ACEi;

// ARB (Losartan → EXP-3174)
dxdt_ARB_gut  = -ka_ARB * ARB_gut;
dxdt_ARB_c    = F_ARB * ka_ARB * ARB_gut - CL_ARB * C_ARB;

// Finerenone (2-compartment)
dxdt_Fin_gut  = -ka_Fin * Fin_gut;
dxdt_Fin_c    = F_Fin * ka_Fin * Fin_gut - (CL_Fin + Q_Fin) * C_Fin + Q_Fin * (Fin_p / V2_Fin);
dxdt_Fin_p    = Q_Fin * C_Fin * V1_Fin   - Q_Fin * (Fin_p / V2_Fin);

// SGLT2 inhibitor (Dapagliflozin)
dxdt_SGLT2_gut = -ka_SGLT2 * SGLT2_gut;
dxdt_SGLT2_c   = F_SGLT2 * ka_SGLT2 * SGLT2_gut - CL_SGLT2 * C_SGLT2;

// ESA (Epoetin alfa SC)
dxdt_ESA_sc   = -ka_ESA * ESA_sc;
dxdt_ESA_c    = F_ESA * ka_ESA * ESA_sc - CL_ESA * C_ESA;

// HIF-PHI (Roxadustat)
dxdt_PHI_gut  = -ka_PHI * PHI_gut;
dxdt_PHI_c    = F_PHI * ka_PHI * PHI_gut - CL_PHI * C_PHI;

// ── NEPHRON MASS & eGFR ───────────────────────────────────────
// Nephron loss driven by glomerulosclerosis + fibrosis
double nephron_loss = kNephron * Nephron * (GFP_effect - SGLT2_GH_red)
                    * (1.0 + 0.5 * Collagen);
dxdt_Nephron  = -nephron_loss;

// eGFR proportional to nephron mass, modified by drug effects
double eGFR_target = eGFR0 * Nephron * (1.0 + 0.12 * E_SGLT2 + 0.08 * E_RAS);
dxdt_eGFR     = 0.05 * (eGFR_target - eGFR);

// UACR — proteinuria driven by GFP and reduced by RAS blockade
double UACR_ss_cur = UACR0 * (GFP_effect - SGLT2_GH_red)
                   * (1.0 - 0.45 * E_RAS) * (1.0 - 0.30 * E_SGLT2)
                   * (1.0 - 0.25 * E_Fin);
dxdt_UACR_st  = 0.02 * (UACR_ss_cur - UACR_st);

// ── RAAS CASCADE ─────────────────────────────────────────────
// Ang II: produced from renin/ACE activity, cleared; ACEi/ARB reduce
double AngII_ss_new = AngII_ss * (1.0 - E_RAS);
dxdt_AngII = kAngII_prod * AngII_ss_new - kAngII_deg * AngII;

// Aldosterone: driven by Ang II; MRA/nsMRA block MR not Aldo production
double Aldo_drive = Aldo_ss * AngII_norm * (1.0 - 0.5 * E_ACEi);
dxdt_Aldo  = kAldo_prod * Aldo_drive - kAldo_deg * Aldo;

// Blood pressure
double BP_ss_target = SBP0 - 18.0 * E_RAS - 5.0 * E_SGLT2 - 3.0 * E_Fin
                    + 20.0 * (AngII_norm - 1.0) + 5.0 * Uremia_tox;
dxdt_BP    = 0.005 * (BP_ss_target - BP);

// ── INFLAMMATION ─────────────────────────────────────────────
// Macrophage infiltration driven by MCP-1 (AngII, uremia)
double Mac_drive = kMacro * (1.0 + 0.5 * AngII_norm + 0.5 * Uremia_tox);
dxdt_Macro = Mac_drive - kMacro_out * Macro;

// IL-6: driven by NF-κB + macrophages
double IL6_prod = kIL6_prod * NF_kB_act * (1.0 + 0.3 * Macro);
dxdt_IL6   = IL6_prod - kIL6_deg * IL6;

// TNF-α: driven by macrophages, NF-κB
double TNFa_prod = kTNFa_prod * NF_kB_act * (1.0 + 0.4 * Macro);
dxdt_TNFa  = TNFa_prod - kTNFa_deg * TNFa;

// ── TGF-β1 / FIBROSIS ────────────────────────────────────────
// TGF-β1: driven by Ang II, M2 macrophages, IL-1β (IL6 proxy); MRA reduces
double TGFb_prod = kTGF_prod * (AngII_norm + 0.3 * IL6_norm + 0.2 * Macro)
                 * (1.0 - 0.35 * E_Fin) * (1.0 - 0.15 * E_RAS);
dxdt_TGFb  = TGFb_prod - kTGF_deg * TGFb;

// Collagen/fibrosis — driven by Smad2/3 activity, cleared slowly
double Col_prod = kCol_prod * Smad23_act;
double Col_deg  = kCol_deg  * Collagen;
dxdt_Collagen = Col_prod - Col_deg;

// ── CKD-MBD ──────────────────────────────────────────────────
// Phosphate rises as eGFR falls
double Phos_target = Phos_ss * (1.0 + kPhos_retent * (eGFR0 / (eGFR + 1.0) - 1.0));
dxdt_Phos  = 0.01 * (Phos_target - Phos);

// Klotho falls with CKD progression; FGF23 co-receptor
double Klotho_target = Klotho_ss * (eGFR / eGFR0);
dxdt_Klotho = kKlotho_syn * Klotho_target - kKlotho_deg * Klotho;

// FGF-23: rises with phosphate and CKD; Klotho required for signaling
double FGF23_drive = kFGF23_prod * (Phos / Phos_ss) * (eGFR0 / (eGFR + 0.1));
dxdt_FGF23 = FGF23_drive - kFGF23_deg * FGF23;

// Active Vitamin D: reduced by FGF-23 and low eGFR
double VitD_ss_cur = 30.0 * (eGFR / eGFR0) / (1.0 + 0.02 * (FGF23 - FGF23_ss));
dxdt_VitD  = kVitD_prod * VitD_ss_cur - kVitD_deg * VitD;

// PTH: driven by low Ca (simplified), high Phos, low Vit D; Klotho modulates
double VitD_norm   = VitD / 30.0;
double PTH_stim    = kPTH_prod * (1.0 + 0.5 * (Phos / Phos_ss - 1.0))
                               / (1.0 + 0.5 * VitD_norm)
                               / (1.0 + 0.3 * Klotho);
dxdt_PTH   = PTH_stim - kPTH_deg * PTH;

// ── ANEMIA ───────────────────────────────────────────────────
// EPO: production falls with eGFR loss; HIF-PHI restores
double EPO_prod = kEPO_prod * EPO_ss * (eGFR / eGFR0) * HIF_stim
                + E_ESA * kEPO_prod * EPO_ss;
dxdt_EPO   = EPO_prod - kEPO_deg * EPO;

// Hepcidin: driven by IL-6 (inflammation); sequesters iron
double Hep_prod = kHep_prod * (1.0 + 1.5 * IL6_norm);
dxdt_Hepcidin = Hep_prod - kHep_deg * Hepcidin;

// Hemoglobin: driven by EPO (erythropoiesis); lost via uremia
double EPO_stim = (EPO / EPO_ss) + E_ESA + 0.8 * E_PHI;
double Hb_prod  = kHb_prod * EPO_stim / (1.0 + 0.2 * (Hepcidin / Hep_ss));
double Hb_loss  = kHb_deg  * Hemoglobin * (1.0 + 0.3 * Uremia_tox);
dxdt_Hemoglobin = Hb_prod - Hb_loss;

// ── CARDIOVASCULAR ────────────────────────────────────────────
// LVH: driven by BP and Aldo; reversed partially by treatment
double LVH_drive = kLVH_prod * (BP / SBP0) * Aldo_norm * (1.0 - 0.4 * E_Fin);
dxdt_LVH_idx = LVH_drive - kLVH_reg * LVH_idx * (E_RAS + 0.5 * E_Fin + 0.3 * E_SGLT2);

// Vascular calcification: driven by Phos, PTH, uremia
double VC_drive = kVC * (Phos / Phos_ss) * PTH_norm * (1.0 + 0.5 * Uremia_tox);
dxdt_VC_idx = VC_drive * (1.0 - VC_idx);

$TABLE
// Output variables
double UACR_out    = UACR_st;
double eGFR_out    = eGFR;
double Hb_out      = Hemoglobin;
double PTH_out     = PTH;
double SBP_out     = BP;
double FGF23_out   = FGF23;
double Phos_out    = Phos;
double Collagen_out = Collagen;
double LVH_out     = LVH_idx;
double Nephron_out = Nephron;

// CKD Stage classification
double CKD_stage;
if      (eGFR >= 90)  CKD_stage = 1;
else if (eGFR >= 60)  CKD_stage = 2;
else if (eGFR >= 45)  CKD_stage = 3;
else if (eGFR >= 30)  CKD_stage = 3;
else if (eGFR >= 15)  CKD_stage = 4;
else                  CKD_stage = 5;

// Anemia classification
double Anemia_grade;
if      (Hemoglobin >= 12) Anemia_grade = 0;
else if (Hemoglobin >= 10) Anemia_grade = 1;
else if (Hemoglobin >= 8)  Anemia_grade = 2;
else                        Anemia_grade = 3;

// CV risk composite
double CV_risk = 0.3 * (BP / SBP0 - 0.5) + 0.3 * (VC_idx) + 0.2 * (LVH_idx - 0.5)
               + 0.2 * (UACR_out / UACR0);

$CAPTURE
eGFR_out UACR_out Hb_out PTH_out SBP_out FGF23_out Phos_out
Collagen_out LVH_out VC_idx Nephron_out CKD_stage Anemia_grade CV_risk
AngII Aldo IL6 TNFa TGFb Klotho VitD EPO Hepcidin
'

# ─────────────────────────────────────────────────────────────
# COMPILE MODEL
# ─────────────────────────────────────────────────────────────
ckd_mod <- mrgsolve::mcode("CKD_QSP", ckd_model_code)

# ─────────────────────────────────────────────────────────────
# INITIAL CONDITIONS FUNCTION
# ─────────────────────────────────────────────────────────────
ckd_init <- function(egfr0 = 45, uacr0 = 300, hb0 = 11.0, sbp0 = 150) {
  list(
    Nephron   = 1.0,
    eGFR      = egfr0,
    UACR_st   = uacr0,
    AngII     = 1.0,
    Aldo      = 1.0,
    IL6       = 1.0,
    TNFa      = 1.0,
    Macro     = 1.0,
    TGFb      = 1.0,
    Collagen  = 0.2,
    FGF23     = 100,
    Klotho    = 1.0,
    Phos      = 3.8,
    PTH       = 65,
    VitD      = 30.0,
    EPO       = 12,
    Hepcidin  = 25,
    Hemoglobin = hb0,
    LVH_idx   = 1.0,
    VC_idx    = 0.05,
    BP        = sbp0
  )
}

# ─────────────────────────────────────────────────────────────
# EVENT TABLES (dosing regimens)
# ─────────────────────────────────────────────────────────────
sim_duration <- 365 * 3  # 3 years in days
sim_step     <- 1        # daily output

# Helper: create dosing event table
make_events <- function(drug, dose_mg, freq_h, start_day = 0, duration_days = sim_duration) {
  cmt_map <- list(
    ACEi  = "ACEi_gut",
    ARB   = "ARB_gut",
    Fin   = "Fin_gut",
    SGLT2 = "SGLT2_gut",
    PHI   = "PHI_gut"
  )
  ev_start  <- start_day * 24
  ev_end    <- (start_day + duration_days) * 24
  mrgsolve::ev(amt = dose_mg,
               cmt = cmt_map[[drug]],
               ii  = freq_h,
               addl = floor((ev_end - ev_start) / freq_h),
               time = ev_start)
}

# ─────────────────────────────────────────────────────────────
# SCENARIO DEFINITIONS (5 clinical scenarios)
# ─────────────────────────────────────────────────────────────
# Scenario parameters: DKD patient, CKD G3a, T2DM
base_params <- list(
  eGFR0 = 45, UACR0 = 300, Hb0 = 11.0, SBP0 = 150, T2DM_flag = 1
)

scenarios <- list(

  # 1. Untreated (natural history)
  Natural_History = list(
    label  = "Natural History (No Treatment)",
    events = mrgsolve::ev(amt = 0, cmt = "ACEi_gut", time = 0),
    color  = "#E74C3C"
  ),

  # 2. Standard of care: ACEi + ARB (not recommended together long-term, but historical)
  ACEi_Monotherapy = list(
    label  = "ACEi Monotherapy (Ramipril 10mg/day)",
    events = make_events("ACEi", dose_mg = 10, freq_h = 24),
    color  = "#3498DB"
  ),

  # 3. Standard of care: ACEi + Finerenone
  ACEi_Finerenone = list(
    label  = "ACEi + Finerenone (Ramipril + Finerenone 20mg)",
    events = mrgsolve::ev_c(
      make_events("ACEi", dose_mg = 10, freq_h = 24),
      make_events("Fin",  dose_mg = 20, freq_h = 24)
    ),
    color  = "#9B59B6"
  ),

  # 4. ACEi + SGLT2i (DAPA-CKD regimen)
  ACEi_SGLT2i = list(
    label  = "ACEi + Dapagliflozin 10mg (DAPA-CKD)",
    events = mrgsolve::ev_c(
      make_events("ACEi",  dose_mg = 10, freq_h = 24),
      make_events("SGLT2", dose_mg = 10, freq_h = 24)
    ),
    color  = "#27AE60"
  ),

  # 5. Triple therapy: ACEi + SGLT2i + Finerenone
  Triple_Therapy = list(
    label  = "Triple Therapy (ACEi + Dapa + Finerenone)",
    events = mrgsolve::ev_c(
      make_events("ACEi",  dose_mg = 10, freq_h = 24),
      make_events("SGLT2", dose_mg = 10, freq_h = 24),
      make_events("Fin",   dose_mg = 20, freq_h = 24)
    ),
    color  = "#F39C12"
  )
)

# ─────────────────────────────────────────────────────────────
# RUN SIMULATIONS
# ─────────────────────────────────────────────────────────────
run_scenario <- function(scen_name, scen) {
  init_vals <- ckd_init(egfr0 = base_params$eGFR0,
                        uacr0 = base_params$UACR0,
                        hb0   = base_params$Hb0,
                        sbp0  = base_params$SBP0)

  tobs <- seq(0, sim_duration * 24, by = 24)  # hourly time, daily output

  out <- ckd_mod %>%
    mrgsolve::init(init_vals) %>%
    mrgsolve::ev(scen$events) %>%
    mrgsolve::mrgsim(end = sim_duration * 24, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(
      Day      = time / 24,
      Year     = Day / 365,
      Scenario = scen$label,
      Color    = scen$color
    )
  out
}

cat("Running CKD QSP simulations...\n")
results_list <- lapply(names(scenarios), function(sn) {
  cat("  Scenario:", sn, "\n")
  run_scenario(sn, scenarios[[sn]])
})
results <- bind_rows(results_list)

# ─────────────────────────────────────────────────────────────
# ANEMIA SCENARIO: ESA vs HIF-PHI vs Natural
# ─────────────────────────────────────────────────────────────
anemia_scenarios <- list(
  No_Anemia_Tx = list(
    label  = "No Anemia Treatment",
    events = mrgsolve::ev(amt = 0, cmt = "ACEi_gut", time = 0),
    color  = "#E74C3C"
  ),
  ESA_TIW = list(
    label  = "Epoetin alfa 4000 IU TIW (SC)",
    events = {
      doses_per_week <- 3
      ev_sc <- mrgsolve::ev(
        amt  = 4000,
        cmt  = "ESA_sc",
        ii   = 24 * 7 / doses_per_week,
        addl = floor(sim_duration / (7 / doses_per_week)) - 1,
        time = 0
      )
      ev_sc
    },
    color  = "#3498DB"
  ),
  HIF_PHI_TIW = list(
    label  = "Roxadustat 100mg TIW (oral)",
    events = mrgsolve::ev(
      amt  = 100,
      cmt  = "PHI_gut",
      ii   = 24 * 7 / 3,
      addl = floor(sim_duration / (7/3)) - 1,
      time = 0
    ),
    color  = "#27AE60"
  )
)

anemia_results <- bind_rows(lapply(names(anemia_scenarios), function(sn) {
  cat("  Anemia Scenario:", sn, "\n")
  init_vals <- ckd_init()
  out <- ckd_mod %>%
    mrgsolve::init(init_vals) %>%
    mrgsolve::ev(anemia_scenarios[[sn]]$events) %>%
    mrgsolve::mrgsim(end = sim_duration * 24, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(
      Day      = time / 24,
      Scenario = anemia_scenarios[[sn]]$label,
      Color    = anemia_scenarios[[sn]]$color
    )
  out
}))

# ─────────────────────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────────────────────
theme_ckd <- function() {
  theme_bw(base_size = 11) +
    theme(
      legend.position   = "bottom",
      legend.text       = element_text(size = 8),
      legend.title      = element_blank(),
      strip.background  = element_rect(fill = "#2C3E50"),
      strip.text        = element_text(color = "white", face = "bold"),
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 12),
      plot.subtitle     = element_text(size = 9, color = "gray40")
    )
}

scenario_colors <- setNames(
  sapply(scenarios, `[[`, "color"),
  sapply(scenarios, `[[`, "label")
)

# Helper function for scenario plots
plot_var <- function(data, yvar, ylabel, title, hline = NULL) {
  p <- ggplot(data, aes(x = Day, y = .data[[yvar]],
                        color = Scenario, group = Scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    labs(title = title, x = "Day", y = ylabel) +
    theme_ckd()
  if (!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline$val, linetype = hline$lty, color = hline$col)
  }
  p
}

# Main outcome plots
p_eGFR <- plot_var(results, "eGFR_out", "eGFR (mL/min/1.73m²)",
                   "A. eGFR Trajectory", hline = list(val = 15, lty = "dashed", col = "#E74C3C"))

p_UACR <- plot_var(results, "UACR_out", "UACR (mg/g Cr)",
                   "B. Proteinuria (UACR)") +
  geom_hline(yintercept = 300, linetype = "dashed", color = "gray50") +
  scale_y_log10()

p_BP   <- plot_var(results, "SBP_out", "SBP (mmHg)", "C. Systolic Blood Pressure") +
  geom_hline(yintercept = 130, linetype = "dashed", color = "gray50")

p_PTH  <- plot_var(results, "PTH_out", "Intact PTH (pg/mL)", "D. PTH (CKD-MBD)")

p_Phos <- plot_var(results, "Phos_out", "Phosphate (mg/dL)", "E. Serum Phosphate") +
  geom_hline(yintercept = 5.5, linetype = "dashed", color = "#E74C3C")

p_LVH  <- plot_var(results, "LVH_out", "LVH Index (normalized)", "F. Left Ventricular Hypertrophy")

p_Fib  <- plot_var(results, "Collagen_out", "Fibrosis Index (0–1)", "G. Renal Fibrosis (Collagen)")

p_CV   <- plot_var(results, "CV_risk", "CV Risk Score (composite)", "H. Cardiovascular Risk")

# Anemia plots
anemia_colors <- setNames(
  sapply(anemia_scenarios, `[[`, "color"),
  sapply(anemia_scenarios, `[[`, "label")
)

p_Hb <- ggplot(anemia_results, aes(x = Day, y = Hb_out, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = anemia_colors) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 12, linetype = "dashed", color = "#27AE60") +
  labs(title = "I. Hemoglobin Response: ESA vs. HIF-PHI",
       x = "Day", y = "Hemoglobin (g/dL)") +
  theme_ckd()

p_EPO <- ggplot(anemia_results, aes(x = Day, y = EPO, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = anemia_colors) +
  labs(title = "J. Plasma EPO Dynamics", x = "Day", y = "EPO (mIU/mL)") +
  theme_ckd()

# CKD-MBD plots
p_FGF23 <- plot_var(results, "FGF23_out", "FGF-23 (pg/mL)", "K. FGF-23 Progression")
p_Klotho <- plot_var(results, "Klotho", "Klotho (normalized)", "L. Klotho Decline")

# Combine main figure
main_fig <- (p_eGFR | p_UACR | p_BP | p_PTH) /
            (p_Phos | p_LVH | p_Fib | p_CV) /
            (p_Hb   | p_EPO | p_FGF23 | p_Klotho) +
  plot_annotation(
    title    = "Chronic Kidney Disease QSP Model — Treatment Scenario Simulation",
    subtitle = paste0("Patient: eGFR₀=45 mL/min/1.73m², UACR₀=300 mg/g, T2DM, SBP₀=150 mmHg ",
                      "| Duration: 3 years"),
    theme = theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10))
  )

# ─────────────────────────────────────────────────────────────
# DOSE-RESPONSE: finerenone dose vs UACR reduction at 12 months
# ─────────────────────────────────────────────────────────────
fin_doses <- c(5, 10, 20, 40)
fin_dr <- lapply(fin_doses, function(d) {
  ev_d <- mrgsolve::ev_c(
    make_events("ACEi", dose_mg = 10, freq_h = 24),
    mrgsolve::ev(amt = d, cmt = "Fin_gut", ii = 24,
                 addl = 365 - 1, time = 0)
  )
  out <- ckd_mod %>%
    mrgsolve::init(ckd_init()) %>%
    mrgsolve::ev(ev_d) %>%
    mrgsolve::mrgsim(end = 365 * 24, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    filter(abs(time - 365 * 24) < 1) %>%
    summarise(
      Dose     = d,
      UACR_pct = 100 * (UACR0 - UACR_out) / UACR0,
      eGFR_ch  = eGFR_out - eGFR0,
      PTH_ch   = PTH_out  - PTH_ss
    )
  out
})
fin_dr_df <- bind_rows(fin_dr)

p_dr_fin <- ggplot(fin_dr_df, aes(x = Dose, y = UACR_pct)) +
  geom_line(color = "#9B59B6", linewidth = 1.2) +
  geom_point(color = "#9B59B6", size = 3) +
  labs(title = "Finerenone Dose–Response: % UACR Reduction at 12 Months",
       x = "Finerenone Dose (mg/day)", y = "% UACR Reduction vs Baseline") +
  theme_ckd()

# ─────────────────────────────────────────────────────────────
# SUMMARY TABLE: 2-year outcomes
# ─────────────────────────────────────────────────────────────
summary_2yr <- results %>%
  filter(abs(Day - 730) < 2) %>%
  group_by(Scenario) %>%
  summarise(
    eGFR_2yr    = round(mean(eGFR_out), 1),
    UACR_2yr    = round(mean(UACR_out), 0),
    SBP_2yr     = round(mean(SBP_out),  1),
    Hb_2yr      = round(mean(Hb_out),   1),
    PTH_2yr     = round(mean(PTH_out),  0),
    CKD_Stage   = round(mean(CKD_stage), 1),
    CV_risk_2yr = round(mean(CV_risk),  3),
    .groups     = "drop"
  )

cat("\n═══════════════════════════════════════════════════════\n")
cat("CKD QSP Model — 2-Year Outcome Summary\n")
cat("═══════════════════════════════════════════════════════\n")
print(summary_2yr, n = Inf)

cat("\n─────────────────────────────────────────────────────\n")
cat("Finerenone Dose–Response (vs ACEi background, 12 months)\n")
cat("─────────────────────────────────────────────────────\n")
print(fin_dr_df)

# ─────────────────────────────────────────────────────────────
# SAVE OUTPUTS
# ─────────────────────────────────────────────────────────────
if (!dir.exists("output")) dir.create("output", recursive = TRUE)

ggsave("output/ckd_qsp_main_panel.pdf", plot = main_fig,
       width = 18, height = 14, dpi = 150)
ggsave("output/ckd_dose_response_finerenone.pdf", plot = p_dr_fin,
       width = 7, height = 5)

cat("\nSimulation complete. Plots saved to output/.\n")
