## =============================================================================
## Mixed Connective Tissue Disease (MCTD) — mrgsolve QSP Model
## =============================================================================
## Disease: Mixed Connective Tissue Disease (혼합결합조직병)
## Hallmark: Anti-U1-RNP autoantibodies; overlap of SLE, SSc, PM features
## Key organs: Vasculature, Lung (PAH, ILD), Muscle (myositis), Joints
## Treatments modeled:
##   1. Hydroxychloroquine (HCQ) 400 mg/day
##   2. Prednisone (Pred) 5–60 mg/day
##   3. Mycophenolate mofetil (MMF) 1–3 g/day
##   4. Rituximab (RTX) 1000 mg x2 doses
##   5. Bosentan (ERA) 125 mg BID — for PAH component
##
## Key References:
##   - Alarcon-Segovia et al. 1987 (MCTD definition)
##   - Kim P et al. Am J Med 2021 (MCTD natural history)
##   - Hajas A et al. Arthritis Res Ther 2013 (MCTD outcomes)
##   - Sharp GC et al. Am J Med 1972 (original MCTD description)
##   - Ruaro B et al. J Immunol Res 2021 (MCTD PAH)
##   - Yaniv G et al. Semin Arthritis Rheum 2016 (overlap disease)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─── Model Definition ────────────────────────────────────────────────────────

mctd_model <- '
$PROB MCTD QSP Model: Anti-U1-RNP Autoimmunity + Multi-organ PK/PD

$PARAM
// ─── Disease Biology Parameters ────────────────────────────────
// Anti-U1-RNP autoantibody dynamics
kAb_prod = 0.005    // Ab production rate (AU/mL/day) — baseline
kAb_elim = 0.0693   // Ab elimination (t1/2 = 10 days, long-lived plasma cells)
kAb_LLPC  = 0.002   // Long-lived plasma cell contribution rate
Ab0       = 300     // Initial anti-U1-RNP titer (AU/mL; MCTD range: 100-1000)

// T cell compartments (cells/μL)
kTh1_act  = 0.15    // Th1 activation rate
kTh17_act = 0.12    // Th17 activation rate
kTreg_act = 0.08    // Treg activation rate
kT_death  = 0.10    // T cell death rate (all effectors)
Th1_0     = 150     // Baseline Th1 cells/μL
Th17_0    = 60      // Baseline Th17 cells/μL
Treg_0    = 80      // Baseline Treg cells/μL

// B cell compartments (cells/μL)
kB_act    = 0.20    // Naive B cell activation rate
kGC       = 0.15    // GC reaction rate
kPC       = 0.10    // Plasma cell differentiation rate
kB_death  = 0.08    // B cell death rate
Bnai0     = 200     // Naive B cells/μL
GCB0      = 30      // GC B cells/μL
PC0       = 20      // Plasmablasts/μL

// Cytokine dynamics (pg/mL)
kTNF_prod = 2.0     // TNF-α production (pg/mL/day)
kIL6_prod = 3.0     // IL-6 production (pg/mL/day)
kIL17_prod= 1.5     // IL-17A production (pg/mL/day)
kIFNg_prod= 1.0     // IFN-γ production (pg/mL/day)
kTGFb_prod= 0.5     // TGF-β production (pg/mL/day)
kCyto_elim= 0.50    // Cytokine elimination rate (t1/2 ~ 1.4 days)
TNF0      = 25      // Baseline TNF-α (pg/mL)
IL6_0     = 15      // Baseline IL-6 (pg/mL)
IL17_0    = 18      // Baseline IL-17A (pg/mL)
IFNg_0    = 8       // Baseline IFN-γ (pg/mL)
TGFb_0    = 5       // Baseline TGF-β (pg/mL)
IFNa_0    = 20      // Baseline IFN-α (pg/mL) — elevated in MCTD

// ─── Vascular Pathology (Raynaud's, PAH) ───────────────────────
ET1_0     = 2.5     // Baseline endothelin-1 (pg/mL; normal <2, MCTD: 2-8)
kET1_prod = 0.05    // ET-1 production rate
kET1_elim = 0.30    // ET-1 elimination (t1/2 ~ 2.3 min, but plasma levels slow)
RVSP_0    = 45      // Baseline RVSP (mmHg; elevated if PAH)
kPVR_ET1  = 0.08    // ET-1 effect on PVR
PVR_0     = 4.0     // Baseline PVR (Wood units; normal <3, PAH: >3)

// ─── Lung compartment (ILD) ────────────────────────────────────
FVC_0     = 75      // Baseline FVC %predicted (MCTD-ILD: often 60-80%)
DLCO_0    = 65      // Baseline DLCO %predicted
kFVC_decl = 0.0003  // Natural FVC decline rate per day
kCollagen = 0.002   // Collagen deposition rate (TGF-β driven)
Collagen0 = 0.3     // Baseline collagen index (0-1 scale)

// ─── Muscle compartment (PM-like myositis) ─────────────────────
CK_0      = 800     // Baseline serum CK (U/L; MCTD myositis: 200-5000)
kCK_prod  = 5.0     // CK release from inflamed muscle (U/L/day)
kCK_elim  = 0.15    // CK elimination (t1/2 ~ 4.6 days)
MMT_0     = 70      // Baseline MMT-8 score (max 80; MCTD: 60-75)

// ─── Joint inflammation (Arthritis) ────────────────────────────
SJC_0     = 6       // Swollen joint count (MCTD: 4-12)
TJC_0     = 8       // Tender joint count
DAS28_0   = 4.5     // Baseline DAS28 score (moderate-high)
kJoint_inflam = 0.03 // Joint inflammation rate

// ─── Complement ────────────────────────────────────────────────
C3_0      = 80      // Baseline C3 (mg/dL; normal 90-180; MCTD: 60-100)
C4_0      = 15      // Baseline C4 (mg/dL; normal 16-47; MCTD: 10-20)
kComp_cons = 0.002  // Complement consumption rate (IC-mediated)
kComp_prod = 0.20   // Complement synthesis rate

// ─── Inflammatory markers ──────────────────────────────────────
ESR_0     = 55      // Baseline ESR (mm/hr; MCTD: 40-80)
CRP_0     = 18      // Baseline CRP (mg/L; MCTD: 10-30)

// ─── PK: Hydroxychloroquine (HCQ) ──────────────────────────────
// Blood compartment (WBC concentration more relevant)
HCQ_F    = 0.79     // Oral bioavailability (79%)
HCQ_Ka   = 0.35     // Absorption rate constant (h^-1)
HCQ_Vd   = 800      // Volume of distribution (L/kg, very high)
HCQ_CL   = 19.5     // Total clearance (L/h)
HCQ_Kp   = 2000     // WBC:plasma partitioning ratio
HCQ_t12  = 1344     // Terminal t1/2 = 40-56 days (hours)

// HCQ PD (TLR inhibition)
HCQ_IC50 = 0.8      // IC50 for TLR7/8 inhibition (μg/mL WBC)
HCQ_Emax = 0.75     // Maximum effect (75% suppression)
HCQ_hill = 1.5      // Hill coefficient

// ─── PK: Mycophenolate (MPA — active metabolite) ───────────────
MPA_F    = 0.94     // Bioavailability of MMF→MPA (~94%)
MPA_Ka   = 1.1      // Absorption rate (h^-1)
MPA_Vd   = 3.6      // Volume of distribution (L/kg)
MPA_CL   = 11.6     // Clearance (L/h)
MPA_t12  = 17.9     // Effective t1/2 (hours)
MPA_EHC  = 0.30     // Enterohepatic circulation fraction

// MPA PD (IMPDH inhibition)
MPA_IC50 = 1.5      // IC50 for IMPDH (mg/L)
MPA_Emax = 0.85     // Maximum B/T cell proliferation inhibition
MPA_hill = 1.2

// ─── PK: Prednisone/Prednisolone ───────────────────────────────
PRED_F   = 0.92     // Oral bioavailability of prednisolone
PRED_Ka  = 0.85     // Absorption rate (h^-1)
PRED_Vd  = 0.55     // Volume of distribution (L/kg)
PRED_CL  = 12.0     // Clearance (L/h)
PRED_t12 = 3.5      // t1/2 = 2-4 h

// Prednisolone PD (NF-κB inhibition)
PRED_IC50 = 0.05    // IC50 for GR-mediated NF-κB suppression (mg/L)
PRED_Emax = 0.90    // Maximum cytokine suppression
PRED_hill = 1.8

// ─── PK: Rituximab ─────────────────────────────────────────────
RTX_Vd1  = 3.5      // Central volume (L)
RTX_Vd2  = 2.1      // Peripheral volume (L)
RTX_CL   = 0.24     // Linear clearance (L/day)
RTX_Q    = 0.51     // Inter-compartmental clearance (L/day)
RTX_CL2  = 0.69     // Target-mediated clearance (L/day; CD20-bound)
RTX_t12  = 22       // Terminal t1/2 (days)

// RTX PD (B cell depletion)
RTX_EC50 = 5.0      // EC50 for B cell depletion (μg/mL)
RTX_Emax = 0.99     // 99% B cell depletion at saturation
RTX_hill = 2.0

// ─── PK: Bosentan (ERA for PAH) ────────────────────────────────
BOS_F    = 0.50     // Oral bioavailability (50%)
BOS_Ka   = 0.80     // Absorption (h^-1)
BOS_Vd   = 0.18     // Vd (L/kg)
BOS_CL   = 2.0      // Clearance (L/h)
BOS_t12  = 5.4      // t1/2 = 5 h

// Bosentan PD (ET-A receptor blockade)
BOS_IC50 = 0.15     // IC50 for ET-1 effects (mg/L)
BOS_Emax = 0.80     // Maximum PVR reduction
BOS_hill = 1.5

// ─── Body weight for dose calculation ──────────────────────────
WT       = 60       // Body weight (kg)

// ─── Dose (mg; switched on by external dose table) ─────────────
DOSE_HCQ  = 0       // mg/day HCQ
DOSE_PRED = 0       // mg/day prednisone
DOSE_MMF  = 0       // mg/day MMF
DOSE_RTX  = 0       // mg rituximab (IV bolus)
DOSE_BOS  = 0       // mg/day bosentan (total)

$CMT
// HCQ compartments
HCQ_GUT HCQ_CENT HCQ_WBC

// MPA compartments
MPA_GUT MPA_CENT MPA_EHC_CMT

// Prednisolone compartments
PRED_GUT PRED_CENT

// Rituximab compartments (2-cpt + target)
RTX_C1 RTX_C2

// Bosentan compartments
BOS_GUT BOS_CENT

// Immune cell compartments
Th1_CMT Th17_CMT Treg_CMT Bnai_CMT GCB_CMT PC_CMT

// Autoantibody
AbU1RNP

// Cytokines (pg/mL equiv.)
TNF_CMT IL6_CMT IL17_CMT IFNg_CMT TGFb_CMT IFNa_CMT

// Vascular
ET1_CMT PVR_CMT

// Organ function
FVC_CMT DLCO_CMT Collagen_CMT CK_CMT MMT_CMT SJC_CMT

// Complement
C3_CMT C4_CMT

$MAIN
// Drug concentrations
double HCQ_Cp = HCQ_CENT / (HCQ_Vd * WT);    // μg/mL plasma
double HCQ_Cw = HCQ_WBC / (HCQ_Vd * WT * 0.001); // μg/mL WBC (much higher)
double MPA_Cp = MPA_CENT / (MPA_Vd * WT);     // mg/L
double PRED_Cp = PRED_CENT / (PRED_Vd * WT);  // mg/L
double RTX_Cp = RTX_C1 / RTX_Vd1;             // μg/mL
double BOS_Cp = BOS_CENT / (BOS_Vd * WT);     // mg/L

// Drug effect functions (Emax model)
double HCQ_eff = HCQ_Emax * pow(HCQ_Cw, HCQ_hill) /
                 (pow(HCQ_IC50, HCQ_hill) + pow(HCQ_Cw, HCQ_hill));
double MPA_eff = MPA_Emax * pow(MPA_Cp, MPA_hill) /
                 (pow(MPA_IC50, MPA_hill) + pow(MPA_Cp, MPA_hill));
double PRED_eff = PRED_Emax * pow(PRED_Cp, PRED_hill) /
                  (pow(PRED_IC50, PRED_hill) + pow(PRED_Cp, PRED_hill));
double RTX_eff = RTX_Emax * pow(RTX_Cp, RTX_hill) /
                 (pow(RTX_EC50, RTX_hill) + pow(RTX_Cp, RTX_hill));
double BOS_eff = BOS_Emax * pow(BOS_Cp, BOS_hill) /
                 (pow(BOS_IC50, BOS_hill) + pow(BOS_Cp, BOS_hill));

// Feedback from disease (immune cells driving inflammation)
double Th_inflam = (Th1_CMT + Th17_CMT) / (Th1_0 + Th17_0);
double Ab_drive  = AbU1RNP / Ab0;

// Initial conditions (set at start)
_F(HCQ_GUT)  = HCQ_F;
_F(MPA_GUT)  = MPA_F;
_F(PRED_GUT) = PRED_F;
_F(BOS_GUT)  = BOS_F;

$ODE
// ─── HCQ PK ─────────────────────────────────────────────
dxdt_HCQ_GUT  = -HCQ_Ka * HCQ_GUT;
dxdt_HCQ_CENT = HCQ_Ka * HCQ_GUT - (HCQ_CL / (HCQ_Vd * WT)) * HCQ_CENT;
dxdt_HCQ_WBC  = HCQ_Kp * (HCQ_Ka * HCQ_GUT) -
                (HCQ_CL / (HCQ_Vd * WT)) * HCQ_WBC;

// ─── MPA PK ─────────────────────────────────────────────
double MPA_abs = MPA_Ka * MPA_GUT;
double MPA_EHC_rel = MPA_EHC * (MPA_CL / (MPA_Vd * WT)) * MPA_CENT;
dxdt_MPA_GUT    = -MPA_abs;
dxdt_MPA_CENT   = MPA_abs - (MPA_CL / (MPA_Vd * WT)) * MPA_CENT + MPA_EHC_rel;
dxdt_MPA_EHC_CMT= (1 - MPA_EHC) * (MPA_CL / (MPA_Vd * WT)) * MPA_CENT -
                   MPA_EHC_rel;

// ─── Prednisolone PK ─────────────────────────────────────
dxdt_PRED_GUT  = -PRED_Ka * PRED_GUT;
dxdt_PRED_CENT = PRED_Ka * PRED_GUT - (PRED_CL / (PRED_Vd * WT)) * PRED_CENT;

// ─── Rituximab PK (2-cpt + TMDD) ─────────────────────────
dxdt_RTX_C1 = -(RTX_CL / RTX_Vd1) * RTX_C1 - RTX_Q * RTX_C1 / RTX_Vd1 +
               RTX_Q * RTX_C2 / RTX_Vd2 -
               RTX_CL2 * RTX_eff * RTX_C1 / RTX_Vd1;
dxdt_RTX_C2 = RTX_Q * RTX_C1 / RTX_Vd1 - RTX_Q * RTX_C2 / RTX_Vd2;

// ─── Bosentan PK ─────────────────────────────────────────
dxdt_BOS_GUT  = -BOS_Ka * BOS_GUT;
dxdt_BOS_CENT = BOS_Ka * BOS_GUT - (BOS_CL / (BOS_Vd * WT)) * BOS_CENT;

// ─── T Cell Dynamics ─────────────────────────────────────
// Th1: driven by IFN-α, IL-12; suppressed by HCQ, Pred, Treg
double Th1_prod = kTh1_act * Th1_0 * (1 + 0.5 * IFNa_CMT / IFNa_0) *
                  (1 - PRED_eff * 0.7) * (1 - HCQ_eff * 0.4);
double Th1_death= kT_death * Th1_CMT * (1 + 0.3 * Treg_CMT / Treg_0);
dxdt_Th1_CMT  = Th1_prod - Th1_death;

// Th17: driven by IL-6, TGF-β; suppressed by Pred, Treg
double Th17_prod = kTh17_act * Th17_0 * (1 + 0.3 * IL6_CMT / IL6_0) *
                   (1 - PRED_eff * 0.8) * (1 - MPA_eff * 0.5);
double Th17_death= kT_death * Th17_CMT;
dxdt_Th17_CMT = Th17_prod - Th17_death;

// Treg: driven by TGF-β, IL-2; boosted by HCQ (restore tolerance)
double Treg_prod = kTreg_act * Treg_0 * (1 + 0.2 * TGFb_CMT / TGFb_0) *
                   (1 + HCQ_eff * 0.2);
double Treg_death= kT_death * 0.5 * Treg_CMT;
dxdt_Treg_CMT = Treg_prod - Treg_death;

// ─── B Cell Dynamics ─────────────────────────────────────
// Naive B: activated by BCR + IL-21; suppressed by RTX, MPA
double Bnai_prod  = kB_act * Bnai0 * (1 - RTX_eff) * (1 - MPA_eff * 0.8);
double Bnai_death = kB_death * Bnai_CMT;
double Bnai_to_GC = kGC * Bnai_CMT * (1 - MPA_eff * 0.7);
dxdt_Bnai_CMT = Bnai_prod - Bnai_death - Bnai_to_GC;

// GC B cells → plasma cells
double GCB_prod  = Bnai_to_GC;
double GCB_death = kB_death * GCB_CMT;
double GCB_to_PC = kPC * GCB_CMT;
dxdt_GCB_CMT = GCB_prod - GCB_death - GCB_to_PC;

// Plasmablasts/Long-lived plasma cells (anti-U1-RNP producers)
double PC_prod  = GCB_to_PC + kAb_LLPC * Ab0;
double PC_death = kB_death * PC_CMT;
dxdt_PC_CMT = PC_prod - PC_death;

// ─── Anti-U1-RNP Autoantibody ─────────────────────────────
// Produced by PC; suppressed by RTX (B depletion), HCQ
double Ab_prod = kAb_prod * PC_CMT * (1 - HCQ_eff * 0.4);
double Ab_elim = kAb_elim * AbU1RNP;
dxdt_AbU1RNP = Ab_prod - Ab_elim;

// ─── Cytokine Dynamics ────────────────────────────────────
// TNF-α
dxdt_TNF_CMT = kTNF_prod * Th_inflam * (1 - PRED_eff * 0.85) -
               kCyto_elim * TNF_CMT;

// IL-6 (drives Th17, acute phase)
dxdt_IL6_CMT = kIL6_prod * Th_inflam * (1 + 0.2 * TNF_CMT / TNF0) *
               (1 - PRED_eff * 0.80) -
               kCyto_elim * IL6_CMT;

// IL-17A
dxdt_IL17_CMT= kIL17_prod * Th17_CMT / Th17_0 *
               (1 - PRED_eff * 0.75) -
               kCyto_elim * IL17_CMT;

// IFN-γ (Th1-derived, myositis driver)
dxdt_IFNg_CMT= kIFNg_prod * Th1_CMT / Th1_0 *
               (1 - PRED_eff * 0.70) -
               kCyto_elim * IFNg_CMT;

// TGF-β (fibrosis driver)
dxdt_TGFb_CMT= kTGFb_prod * (1 + 0.3 * IL6_CMT / IL6_0) *
               (1 - PRED_eff * 0.50) -
               kCyto_elim * TGFb_CMT;

// IFN-α (type I interferon, pDC derived, suppressed by HCQ)
dxdt_IFNa_CMT= 0.15 * Ab_drive * (1 - HCQ_eff * 0.70) - 0.4 * IFNa_CMT;

// ─── Vascular Pathology ───────────────────────────────────
// ET-1: increased by IL-6, IFN-α; blocked by bosentan (ERA)
dxdt_ET1_CMT = kET1_prod * (1 + 0.3 * IL6_CMT / IL6_0 +
                              0.2 * IFNa_CMT / IFNa_0) *
               (1 - BOS_eff) - kET1_elim * ET1_CMT;

// PVR (Wood units): ET-1 driven; bosentan reduces
dxdt_PVR_CMT = kPVR_ET1 * ET1_CMT * (1 - BOS_eff * 0.6) -
               0.05 * (PVR_CMT - 1.5); // remodeling plateau

// ─── Lung Function ────────────────────────────────────────
// Collagen deposition (ILD fibrosis): TGF-β, IL-13 driven
dxdt_Collagen_CMT = kCollagen * TGFb_CMT / TGFb_0 *
                    (1 - PRED_eff * 0.40) * (1 - MPA_eff * 0.50) -
                    0.001 * Collagen_CMT;

// FVC decline: TGF-β driven fibrosis; partially reversible with treatment
dxdt_FVC_CMT = -kFVC_decl * (1 + 5 * Collagen_CMT) * FVC_CMT +
               0.01 * PRED_eff * (FVC_0 - FVC_CMT);

// DLCO decline: similar, correlates with vascular + parenchymal disease
dxdt_DLCO_CMT= -kFVC_decl * 0.8 * (1 + 3 * PVR_CMT / PVR_0) * DLCO_CMT +
               0.008 * PRED_eff * (DLCO_0 - DLCO_CMT);

// ─── Myositis / CK ────────────────────────────────────────
// CK: IFN-γ, CD8+ T cells damage muscle
dxdt_CK_CMT = kCK_prod * IFNg_CMT / IFNg_0 *
              (1 - PRED_eff * 0.70) * (1 - MPA_eff * 0.50) -
              kCK_elim * CK_CMT;

// MMT score: improves with CK reduction
dxdt_MMT_CMT = 0.005 * (80 - MMT_CMT) * PRED_eff -
               0.002 * CK_CMT / CK_0 * (80 - MMT_0);

// ─── Joint Inflammation ────────────────────────────────────
dxdt_SJC_CMT = kJoint_inflam * Th_inflam * (1 - PRED_eff * 0.75) *
               (1 - HCQ_eff * 0.50) -
               0.05 * SJC_CMT;

// ─── Complement (consumption vs. synthesis) ───────────────
dxdt_C3_CMT = kComp_prod * C3_0 - kComp_cons * AbU1RNP / Ab0 *
              C3_CMT - 0.10 * C3_CMT;
dxdt_C4_CMT = kComp_prod * 0.5 * C4_0 - kComp_cons * 1.5 *
              AbU1RNP / Ab0 * C4_CMT - 0.12 * C4_CMT;

$TABLE
// ─── Drug concentrations ───────────────────────────────────
double HCQ_plasma = HCQ_CENT / (HCQ_Vd * WT);   // μg/mL
double MPA_plasma = MPA_CENT / (MPA_Vd * WT);    // mg/L
double PRED_plasma= PRED_CENT / (PRED_Vd * WT);  // mg/L
double RTX_plasma = RTX_C1 / RTX_Vd1;            // μg/mL
double BOS_plasma = BOS_CENT / (BOS_Vd * WT);    // mg/L

// ─── Derived PD Endpoints ──────────────────────────────────
// ESR (mm/hr): driven by fibrinogen (IL-6 effect)
double ESR   = ESR_0 * (0.5 + 0.5 * IL6_CMT / IL6_0) *
               (1 - PRED_eff * 0.50);

// CRP (mg/L): acute phase, IL-6 driven
double CRP   = CRP_0 * (0.4 + 0.6 * IL6_CMT / IL6_0) *
               (1 - PRED_eff * 0.60);

// DAS28 (joints + ESR)
double DAS28 = 0.56 * sqrt(SJC_CMT) + 0.28 * sqrt(SJC_CMT * 1.3) +
               0.70 * log(ESR + 1) + 0.014 * 50;

// RVSP estimated (mmHg): from PVR and cardiac output (CO=5 L/min)
double RVSP  = PVR_CMT * 5 + 8; // RVSP ≈ PVR × CO + RVEDP

// 6MWT estimate (m): function of WHO FC (PVR/RVSP proxy)
double SixMWT= 450 - 80 * (PVR_CMT - 2.0);

// WHO FC estimate: rough conversion from RVSP
double WHO_FC_est = (RVSP < 40) ? 1.0 : (RVSP < 55) ? 2.0 :
                   (RVSP < 70) ? 3.0 : 4.0;

// MCTD Activity Index (simplified Yamanaka criteria)
double MCTD_AI = (AbU1RNP > Ab0 * 1.5 ? 1 : 0) +
                 (CK_CMT > 1000 ? 2 : CK_CMT > 500 ? 1 : 0) +
                 (RVSP > 40 ? 2 : RVSP > 35 ? 1 : 0) +
                 (FVC_CMT < 70 ? 2 : FVC_CMT < 80 ? 1 : 0) +
                 (SJC_CMT > 8 ? 1 : 0) +
                 (ESR > 60 ? 1 : 0);

capture HCQ_plasma MPA_plasma PRED_plasma RTX_plasma BOS_plasma
capture Th1_CMT Th17_CMT Treg_CMT Bnai_CMT GCB_CMT PC_CMT
capture AbU1RNP TNF_CMT IL6_CMT IL17_CMT IFNg_CMT TGFb_CMT IFNa_CMT
capture ET1_CMT PVR_CMT RVSP SixMWT WHO_FC_est
capture FVC_CMT DLCO_CMT Collagen_CMT
capture CK_CMT MMT_CMT SJC_CMT DAS28
capture C3_CMT C4_CMT ESR CRP MCTD_AI

$INIT
HCQ_GUT = 0
HCQ_CENT = 0
HCQ_WBC = 0
MPA_GUT = 0
MPA_CENT = 0
MPA_EHC_CMT = 0
PRED_GUT = 0
PRED_CENT = 0
RTX_C1 = 0
RTX_C2 = 0
BOS_GUT = 0
BOS_CENT = 0
Th1_CMT = 150
Th17_CMT = 60
Treg_CMT = 80
Bnai_CMT = 200
GCB_CMT = 30
PC_CMT = 20
AbU1RNP = 300
TNF_CMT = 25
IL6_CMT = 15
IL17_CMT = 18
IFNg_CMT = 8
TGFb_CMT = 5
IFNa_CMT = 20
ET1_CMT = 2.5
PVR_CMT = 4.0
FVC_CMT = 75
DLCO_CMT = 65
Collagen_CMT = 0.3
CK_CMT = 800
MMT_CMT = 70
SJC_CMT = 6
C3_CMT = 80
C4_CMT = 15
'

## ─── Compile the model ───────────────────────────────────────────────────────
mod <- mcode("MCTD_QSP", mctd_model)

## ─── Simulation Helper ───────────────────────────────────────────────────────
run_sim <- function(mod, dose_hcq = 0, dose_pred = 0, dose_mmf = 0,
                    dose_rtx = 0, dose_bos = 0, duration_days = 365) {

  # Build dosing events
  ev_list <- list()

  if (dose_hcq > 0) {
    ev_list[["HCQ"]] <- ev(amt = dose_hcq, cmt = "HCQ_GUT",
                           ii = 24, addl = duration_days - 1, time = 0)
  }
  if (dose_pred > 0) {
    ev_list[["PRED"]] <- ev(amt = dose_pred, cmt = "PRED_GUT",
                            ii = 24, addl = duration_days - 1, time = 0)
  }
  if (dose_mmf > 0) {
    # BID dosing
    ev_list[["MMF"]] <- ev(amt = dose_mmf / 2, cmt = "MPA_GUT",
                           ii = 12, addl = 2 * duration_days - 1, time = 0)
  }
  if (dose_rtx > 0) {
    # Two doses 2 weeks apart
    ev_list[["RTX"]] <- ev(amt = dose_rtx, cmt = "RTX_C1", time = 0) +
                         ev(amt = dose_rtx, cmt = "RTX_C1", time = 14 * 24)
  }
  if (dose_bos > 0) {
    # BID dosing
    ev_list[["BOS"]] <- ev(amt = dose_bos / 2, cmt = "BOS_GUT",
                           ii = 12, addl = 2 * duration_days - 1, time = 0)
  }

  if (length(ev_list) == 0) {
    ev_total <- ev(time = 0, amt = 0, cmt = 1)  # dummy
  } else {
    ev_total <- Reduce("+", ev_list)
  }

  out <- mod %>%
    ev(ev_total) %>%
    mrgsim(end = duration_days * 24, delta = 24) %>%
    as_tibble() %>%
    mutate(time_days = time / 24)

  return(out)
}

## ─── 5 Treatment Scenarios ───────────────────────────────────────────────────

cat("\n=== MCTD QSP Model: 5 Treatment Scenarios ===\n")

# Scenario 1: No treatment (disease natural history)
cat("Running Scenario 1: No treatment...\n")
sc1 <- run_sim(mod, duration_days = 365)
sc1$Scenario <- "1. No Treatment (Natural History)"

# Scenario 2: HCQ monotherapy (standard of care, mild MCTD)
cat("Running Scenario 2: HCQ 400 mg/day...\n")
sc2 <- run_sim(mod, dose_hcq = 400, duration_days = 365)
sc2$Scenario <- "2. HCQ 400 mg/day"

# Scenario 3: Moderate disease — HCQ + Prednisone + MMF
cat("Running Scenario 3: HCQ + Pred 20 mg/day + MMF 2 g/day...\n")
sc3 <- run_sim(mod, dose_hcq = 400, dose_pred = 20, dose_mmf = 2000,
               duration_days = 365)
sc3$Scenario <- "3. HCQ + Pred (20mg) + MMF (2g/day)"

# Scenario 4: Severe/refractory — add Rituximab (B cell depletion)
cat("Running Scenario 4: HCQ + Pred + MMF + Rituximab 1000 mg x2...\n")
sc4 <- run_sim(mod, dose_hcq = 400, dose_pred = 20, dose_mmf = 2000,
               dose_rtx = 1000, duration_days = 365)
sc4$Scenario <- "4. HCQ + Pred + MMF + RTX 1000mg×2"

# Scenario 5: PAH dominant — Add Bosentan ERA therapy
cat("Running Scenario 5: HCQ + Pred + MMF + Bosentan 250 mg/day (PAH)...\n")
sc5 <- run_sim(mod, dose_hcq = 400, dose_pred = 20, dose_mmf = 2000,
               dose_bos = 250, duration_days = 365)
sc5$Scenario <- "5. HCQ + Pred + MMF + Bosentan (PAH)"

## ─── Combine results ──────────────────────────────────────────────────────────
all_sc <- bind_rows(sc1, sc2, sc3, sc4, sc5)

## ─── Summary Table ───────────────────────────────────────────────────────────
summary_tbl <- all_sc %>%
  filter(time_days %in% c(0, 30, 90, 180, 365)) %>%
  select(Scenario, time_days, AbU1RNP, CK_CMT, FVC_CMT, DLCO_CMT,
         RVSP, PVR_CMT, SJC_CMT, DAS28, ESR, CRP, MCTD_AI,
         Bnai_CMT, Th17_CMT) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

cat("\n=== Clinical Summary Table ===\n")
print(summary_tbl, n = Inf)

## ─── Plots ───────────────────────────────────────────────────────────────────
theme_mctd <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "#2E4057", color = NA),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid.minor = element_blank())

colors <- c("#E74C3C", "#3498DB", "#2ECC71", "#9B59B6", "#E67E22")

p1 <- ggplot(all_sc, aes(x = time_days, y = AbU1RNP, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = c(100, 300, 600), linetype = "dashed",
             alpha = 0.4, color = "gray50") +
  annotate("text", x = 5, y = 300, label = "Baseline", size = 3, hjust = 0) +
  labs(title = "Anti-U1-RNP Antibody Titer (MCTD Biomarker)",
       x = "Time (days)", y = "Anti-U1-RNP (AU/mL)",
       color = NULL) +
  theme_mctd

p2 <- ggplot(all_sc, aes(x = time_days, y = CK_CMT, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = c(200, 1000), linetype = "dashed",
             alpha = 0.4, color = "red") +
  labs(title = "Serum CK — Myositis Marker",
       x = "Time (days)", y = "CK (U/L)", color = NULL) +
  theme_mctd

p3 <- ggplot(all_sc, aes(x = time_days, y = FVC_CMT, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = c(80, 70), linetype = "dashed",
             alpha = 0.5, color = "navy") +
  annotate("text", x = 5, y = 80, label = "80% predicted",
           size = 3, hjust = 0, color = "navy") +
  labs(title = "FVC % Predicted — ILD Lung Function",
       x = "Time (days)", y = "FVC (%pred)", color = NULL) +
  theme_mctd

p4 <- ggplot(all_sc, aes(x = time_days, y = RVSP, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = c(35, 50, 70), linetype = "dashed",
             alpha = 0.4, color = c("gold", "orange", "red")) +
  annotate("text", x = 5, y = 35, label = "PAH threshold",
           size = 3, hjust = 0) +
  labs(title = "RVSP — Pulmonary Arterial Hypertension",
       x = "Time (days)", y = "RVSP (mmHg)", color = NULL) +
  theme_mctd

p5 <- all_sc %>%
  select(time_days, Scenario, Th1_CMT, Th17_CMT, Treg_CMT) %>%
  pivot_longer(cols = c(Th1_CMT, Th17_CMT, Treg_CMT),
               names_to = "Cell", values_to = "Count") %>%
  mutate(Cell = recode(Cell,
    "Th1_CMT" = "Th1", "Th17_CMT" = "Th17", "Treg_CMT" = "Treg")) %>%
  filter(Scenario %in% c("1. No Treatment (Natural History)",
                          "4. HCQ + Pred + MMF + RTX 1000mg×2")) %>%
  ggplot(aes(x = time_days, y = Count, color = Cell, linetype = Scenario)) +
  geom_line(size = 1.2) +
  labs(title = "T Cell Dynamics: No Treatment vs. Intensive Therapy",
       x = "Time (days)", y = "Cells/μL", color = "Cell Type") +
  theme_mctd

p6 <- ggplot(all_sc, aes(x = time_days, y = MCTD_AI, color = Scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = c(2, 5), linetype = "dashed",
             alpha = 0.5) +
  annotate("text", x = 5, y = 2, label = "Low activity",
           size = 3, hjust = 0) +
  annotate("text", x = 5, y = 5, label = "High activity",
           size = 3, hjust = 0) +
  labs(title = "MCTD Activity Index (Composite Score)",
       x = "Time (days)", y = "MCTD Activity Index", color = NULL) +
  theme_mctd

## Print all plots
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

## ─── Dose-Response Analysis ───────────────────────────────────────────────────
cat("\n=== Dose-Response: Prednisone Dose vs. MCTD Endpoints at Day 90 ===\n")
pred_doses <- c(5, 10, 20, 40, 60)
dr_results <- lapply(pred_doses, function(d) {
  out <- run_sim(mod, dose_hcq = 400, dose_pred = d, dose_mmf = 2000,
                 duration_days = 90)
  out[nrow(out), c("MCTD_AI", "CK_CMT", "FVC_CMT", "RVSP", "ESR",
                   "AbU1RNP", "SJC_CMT")] %>%
    mutate(Prednisone_mg = d)
})
dr_tbl <- bind_rows(dr_results)
cat("Day 90 Endpoints by Prednisone Dose:\n")
print(dr_tbl)

## ─── Bosentan PAH Dose-Response ───────────────────────────────────────────────
cat("\n=== ERA Dose-Response: Bosentan vs. PAH Endpoints at 6 months ===\n")
bos_doses <- c(0, 125, 250, 500)
bos_dr <- lapply(bos_doses, function(d) {
  out <- run_sim(mod, dose_hcq = 400, dose_pred = 20, dose_mmf = 2000,
                 dose_bos = d, duration_days = 180)
  out[nrow(out), c("RVSP", "PVR_CMT", "SixMWT", "ET1_CMT", "FVC_CMT")] %>%
    mutate(Bosentan_mg_day = d)
})
bos_tbl <- bind_rows(bos_dr)
cat("6-Month PAH Endpoints by Bosentan Dose:\n")
print(bos_tbl)

cat("\nMCTD QSP Model simulation complete.\n")
cat("All 5 treatment scenarios and dose-response analyses computed.\n")
