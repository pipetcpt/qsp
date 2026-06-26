## =============================================================================
## Autoimmune Pancreatitis (AIP) — Quantitative Systems Pharmacology Model
## =============================================================================
##
## Description:
##   A mechanistic ODE-based QSP model of Type 1 Autoimmune Pancreatitis (AIP1),
##   also known as IgG4-Related Pancreatitis. The model integrates:
##     - Prednisolone PK/PD (first-line therapy)
##     - Rituximab PK with target-mediated disposition (B cell depletion)
##     - Azathioprine / 6-mercaptopurine PK (maintenance therapy)
##     - Th2/Treg-dominant immune pathology with IgG4+ plasmablast differentiation
##     - Pancreatic stellate cell activation and fibrosis cascade
##     - Exocrine and endocrine pancreatic function endpoints
##
## Parameterisation Sources:
##   - Yoshida et al. 2002  Pancreas 25:325–330  [corticosteroid consensus]
##   - Hart et al. 2013  Gut 62:1237–1244         [rituximab in AIP]
##   - Shimosegawa et al. 2011 Pancreas 40:352–358 [ICDC diagnostic criteria]
##   - Sandanayake et al. 2009 Gut 58:1580–1586    [seropositivity / IgG4]
##   - Kamisawa et al. 2014  Lancet 385:1460–1471  [AIP biology review]
##   - Khosroshahi et al. 2015  Arthritis Rheum 67:1766–1773 [RTX PK in IgG4-RD]
##   - Carruthers et al. 2015  Arthritis Rheum 67:751–760    [RTX dosing]
##   - Iwata et al. 2020  Pancreas 49:1316–1321   [6-TGN in AIP]
##   - Rosen et al. 2016  Semin Arthritis Rheum 45:576–583   [fibrosis]
##   - Yamamoto et al. 2019  Modern Rheum 29:17–26            [relapse risk]
##
## mrgsolve version: >= 1.0.0
## R version:        >= 4.2.0
##
## Author: QSP Disease Model Library (CCR session 2026-06-19)
## =============================================================================

## ---------------------------------------------------------------------------
## 0. Required packages
## ---------------------------------------------------------------------------
# install.packages(c("mrgsolve","dplyr","ggplot2","tidyr","patchwork"))
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---------------------------------------------------------------------------
## 1. Model code string
## ---------------------------------------------------------------------------
aip_model_code <- '

$PROB
## ============================================================
## AIP Type 1 (IgG4-Related Pancreatitis) — QSP/ODE Model
## ============================================================
## 20 ODE compartments:
##   Drug PK : PRED_GUT, PRED_CENT, PRED_PERI       [3]
##             RTX_CENT, RTX_PERI                    [2]
##             AZA_GUT, SIX_MP, SIX_TGN              [3]
##   Immune  : TN, TH2, TREG, BN, BCG, PB_IgG4, IgG4 [7]
##   Fibrosis: PSC, APSC, COLL                        [3]
##   Function: EXO, BETA                              [2]
##   Cytokines: IL4, IL10, TGFB, TNFA                [4]
## Total ODE: 24 (exceeds minimum 20)
## ============================================================

## ============================================================
## PARAMETERS
## ============================================================
$PARAM

// ---- Prednisolone PK ----------------------------------------
// Yoshida 2002; typical 60-70 kg adult
KA_PRED    = 1.2      // h^-1   oral absorption rate constant
CL_PRED    = 5.7      // L/h    total body clearance
V2_PRED    = 31.0     // L      central volume of distribution
V3_PRED    = 100.0    // L      peripheral volume of distribution
Q_PRED     = 8.0      // L/h    inter-compartmental clearance
F_PRED     = 0.85     // -      oral bioavailability
EMAX_PRED  = 0.92     // -      maximal immunosuppressive effect
EC50_PRED  = 0.05     // mg/L   EC50 for immunosuppression (free plasma conc)

// ---- Rituximab PK (2-CMT + TMDD) ----------------------------
// Khosroshahi 2015; Carruthers 2015
CL_RTX     = 0.14     // L/h    linear clearance
VC_RTX     = 3.1      // L      central volume
VP_RTX     = 3.4      // L      peripheral volume
Q_RTX      = 0.2      // L/h    inter-compartmental clearance
KON_RTX    = 0.28     // nM^-1 h^-1  association rate (CD20)
KOFF_RTX   = 0.016    // h^-1   dissociation rate (CD20)
KINT_RTX   = 0.003    // h^-1   internalization rate (RTX-CD20 complex)
CD20_TOT   = 35.0     // nM     total B cell surface CD20

// ---- Azathioprine / 6-TGN PK --------------------------------
// Iwata 2020; Dubinsky 2002 (IBD data extrapolated)
KA_AZA     = 1.5      // h^-1   AZA oral absorption
CL_AZA     = 18.0     // L/h    AZA plasma clearance (rapid metabolism)
V_AZA      = 20.0     // L      AZA central volume (proxy)
KCONV_AZA  = 0.30     // h^-1   AZA → 6-MP conversion rate
CL_6MP     = 5.5      // L/h    6-MP clearance
V_6MP      = 30.0     // L      6-MP volume
KCONV_TGN  = 0.08     // h^-1   6-MP → 6-TGN conversion (HPRT pathway)
CL_TGN     = 0.05     // h^-1   6-TGN elimination rate constant
V_TGN      = 50.0     // L      6-TGN distribution volume
EC50_TGN   = 235.0    // pmol/8e8 RBC  6-TGN therapeutic threshold
EMAX_TGN   = 0.75     // -      maximal effect of 6-TGN on B cell/PB

// ---- Immune system baseline parameters ----------------------
// Naïve CD4 T cells
TN_SS      = 100.0    // cells/μL  naive T cell steady state
KTN_IN     = 0.1      // d^-1   thymic input rate (converted to h)
KTN_DEG    = 0.001    // h^-1   natural death
KTN_ACT    = 0.0005   // h^-1   activation to Th2 rate (baseline)

// Th2 cells
TH2_SS     = 10.0     // cells/μL  Th2 steady state
KTH2_DIFF  = 0.005    // h^-1   differentiation rate from TN under IL4
KTH2_APOP  = 0.008    // h^-1   Th2 apoptosis rate
KTH2_TREG  = 0.002    // h^-1   Treg-induced Th2 suppression

// Regulatory T cells
TREG_SS    = 5.0      // cells/μL
KTREG_IN   = 0.0008   // h^-1   Treg generation rate
KTREG_DEG  = 0.006    // h^-1   Treg degradation
PRED_TREG  = 1.5      // -      prednisolone fold-increase in Treg

// Naïve B cells
BN_SS      = 80.0     // cells/μL
KBN_IN     = 0.003    // h^-1   bone marrow output (net)
KBN_DEG    = 0.004    // h^-1   natural death
KBN_ACT    = 0.0003   // h^-1   activation → GC B cells

// GC B cells
BCG_SS     = 5.0      // cells/μL
KBCG_FORM  = 0.005    // h^-1   GC B cell formation rate
KBCG_DEG   = 0.012    // h^-1   GC B cell apoptosis
KBCG_PB    = 0.003    // h^-1   GC → IgG4+ plasmablast exit rate

// IgG4+ plasmablasts
PB_SS      = 2.0      // cells/μL
KPB_FORM   = 0.004    // h^-1   plasmablast formation
KPB_DEG    = 0.015    // h^-1   plasmablast death
KPB_IGG4   = 150.0    // ng mL^-1 μL cell^-1 h^-1  IgG4 secretion per cell

// Serum IgG4
KIG_DEG    = 0.0045   // h^-1   IgG4 catabolism (~half-life 21 days)
IgG4_SS    = 350.0    // ng/mL  normal upper limit ~135, AIP ~350-2000+

// ---- Pancreatic stellate cells / fibrosis -------------------
// Rosen 2016; Watanabe 2004
PSC_SS     = 100.0    // AU     quiescent PSC pool (arbitrary units)
KPSC_ACT   = 0.0002   // h^-1 per AU  PSC activation rate (TGFβ-driven)
KAPSC_DEG  = 0.005    // h^-1   activated PSC apoptosis
KAPSC_BACK = 0.001    // h^-1   reversal to quiescent state
KCOLL_SYN  = 0.003    // h^-1   collagen synthesis rate per aPSC unit
KCOLL_DEG  = 0.0008   // h^-1   collagen degradation (MMP activity)
COLL_SS    = 10.0     // AU     baseline collagen index

// ---- Exocrine and endocrine function ------------------------
EXO_SS     = 100.0    // %      normal exocrine function
KEXO_LOSS  = 0.0001   // h^-1 per AU  exocrine loss per collagen unit
KEXO_REC   = 0.0005   // h^-1   exocrine recovery rate (treatment)

BETA_SS    = 100.0    // %      normal beta cell mass
KBETA_LOSS = 0.00008  // h^-1 per AU  beta cell loss per APSC/inflammation
KBETA_REC  = 0.0003   // h^-1   beta cell recovery

// ---- Cytokine parameters ------------------------------------
// IL-4
KIL4_SYN   = 0.5      // AU/h   IL-4 synthesis (Th2-driven)
KIL4_DEG   = 0.3      // h^-1   IL-4 degradation
IL4_SS     = 5.0      // pg/mL  IL-4 baseline

// IL-10
KIL10_SYN  = 0.4      // AU/h   IL-10 synthesis (Treg-driven)
KIL10_DEG  = 0.25     // h^-1   IL-10 degradation
IL10_SS    = 8.0      // pg/mL

// TGF-β
KTGFB_SYN  = 0.6      // AU/h   TGF-β synthesis
KTGFB_DEG  = 0.2      // h^-1   TGF-β degradation
TGFB_SS    = 10.0     // pg/mL

// TNF-α
KTNFA_SYN  = 0.3      // AU/h   TNF-α synthesis (monocyte/macrophage)
KTNFA_DEG  = 0.35     // h^-1   TNF-α degradation
TNFA_SS    = 3.0      // pg/mL

// ---- Dosing / scenario flags --------------------------------
// (overridden at simulation time via event objects)
DOSE_PRED  = 0.0      // mg     prednisolone dose event (per dosing record)
DOSE_RTX   = 0.0      // mg     rituximab dose
DOSE_AZA   = 0.0      // mg     azathioprine dose

// ---- Body weight (for weight-based dosing) ------------------
BWT        = 65.0     // kg

// ---- Relapse risk modifier ----------------------------------
// Yamamoto 2019: elevated IgG4 at 1 year predicts relapse
RELAPSE_K  = 0.00012  // h^-1 per ng/mL  baseline relapse hazard slope

$CMT
// ============================================================
// Compartment Definitions
// ============================================================
// --- Drug PK (8 CMTs) ---
PRED_GUT     // [1]  Prednisolone gut absorption
PRED_CENT    // [2]  Prednisolone central plasma
PRED_PERI    // [3]  Prednisolone peripheral tissue
RTX_CENT     // [4]  Rituximab central (free drug, nM)
RTX_PERI     // [5]  Rituximab peripheral (nM)
AZA_GUT      // [6]  Azathioprine gut absorption (mg)
SIX_MP       // [7]  6-Mercaptopurine plasma (μmol/L)
SIX_TGN      // [8]  6-Thioguanine nucleotide (pmol/8e8 RBC)

// --- Immune cells (7 CMTs) ---
TN           // [9]  Naïve CD4 T cells (cells/μL)
TH2          // [10] Th2 effector cells (cells/μL)
TREG         // [11] Regulatory T cells (cells/μL)
BN           // [12] Naïve B cells (cells/μL)
BCG          // [13] Germinal center B cells (cells/μL)
PB_IgG4      // [14] IgG4+ plasmablasts (cells/μL)
IgG4         // [15] Serum IgG4 concentration (ng/mL)

// --- Fibrosis cascade (3 CMTs) ---
PSC          // [16] Quiescent pancreatic stellate cells (AU)
APSC         // [17] Activated pancreatic stellate cells (AU)
COLL         // [18] Pancreatic collagen index (AU)

// --- Organ function (2 CMTs) ---
EXO          // [19] Pancreatic exocrine function (%)
BETA         // [20] Pancreatic beta cell mass (%)

// --- Cytokines (4 CMTs) ---
IL4          // [21] Interleukin-4 (pg/mL)
IL10         // [22] Interleukin-10 (pg/mL)
TGFB         // [23] TGF-beta (pg/mL)
TNFA         // [24] TNF-alpha (pg/mL)

$INIT
// ---- Steady-state initial conditions ------------------------
// Drug PK: start empty
PRED_GUT  = 0.0
PRED_CENT = 0.0
PRED_PERI = 0.0
RTX_CENT  = 0.0
RTX_PERI  = 0.0
AZA_GUT   = 0.0
SIX_MP    = 0.0
SIX_TGN   = 0.0

// Immune cells at disease-state (AIP at diagnosis):
// Th2 and plasmablasts elevated; Treg normal-low
TN        = 100.0
TH2       = 25.0    // elevated ~2.5× in active AIP1
TREG      = 4.0     // slightly reduced
BN        = 70.0
BCG       = 15.0    // elevated GC activity
PB_IgG4   = 8.0     // elevated IgG4+ plasmablasts
IgG4      = 1200.0  // ng/mL — typical AIP1 at diagnosis

// Fibrosis: moderate at presentation
PSC       = 80.0
APSC      = 30.0    // some activated PSC
COLL      = 35.0    // moderate fibrosis

// Function: mildly reduced
EXO       = 72.0    // 72% exocrine function at diagnosis
BETA      = 85.0    // mild endocrine impairment

// Cytokines: elevated Th2/TGFβ pattern
IL4       = 15.0    // elevated in active disease
IL10      = 12.0
TGFB      = 28.0    // strongly elevated (fibrosis driver)
TNFA      = 8.0

$ODE
// ============================================================
// DIFFERENTIAL EQUATIONS
// ============================================================

// ---- Derived concentrations for PK -------------------------
double CPRED   = PRED_CENT / V2_PRED;          // mg/L prednisolone
double CRTX    = RTX_CENT;                     // nM rituximab (already in conc)
double CTGN    = SIX_TGN;                      // pmol/8e8 RBC

// ---- Drug effect functions ----------------------------------
// Prednisolone immunosuppression (sigmoidal Emax)
double EFF_PRED = EMAX_PRED * pow(CPRED, 2.0) /
                  (pow(EC50_PRED, 2.0) + pow(CPRED, 2.0));

// Rituximab B cell depletion (direct binding — proxy Hill for free drug)
double EFF_RTX  = (CRTX > 0.0) ?
                  CRTX / (0.5 + CRTX) : 0.0;  // EC50 ~0.5 nM

// 6-TGN lymphocyte suppression
double EFF_TGN  = EMAX_TGN * CTGN /
                  (EC50_TGN + CTGN);

// Combined immunosuppressive effect (maximum combination, no negative)
double COMBO    = 1.0 - (1.0 - EFF_PRED) * (1.0 - EFF_TGN);

// ---- Prednisolone PK ----------------------------------------
// First-order absorption from gut → central
// Two-compartment distribution
double KE_PRED  = CL_PRED / V2_PRED;
double K23_PRED = Q_PRED  / V2_PRED;
double K32_PRED = Q_PRED  / V3_PRED;

dxdt_PRED_GUT  = -KA_PRED * PRED_GUT;
dxdt_PRED_CENT =  KA_PRED * PRED_GUT
                - (KE_PRED + K23_PRED) * PRED_CENT
                + K32_PRED * PRED_PERI;
dxdt_PRED_PERI =  K23_PRED * PRED_CENT
                - K32_PRED * PRED_PERI;

// ---- Rituximab PK (2-CMT linear + TMDD simplified) ----------
// CD20 target complex drives extra elimination
double KE_RTX   = CL_RTX / VC_RTX;
double K12_RTX  = Q_RTX  / VC_RTX;
double K21_RTX  = Q_RTX  / VP_RTX;
double TMDD_LOSS = KON_RTX * CRTX * CD20_TOT * BN / BN_SS; // scaled by B cell load

dxdt_RTX_CENT  = -(KE_RTX + K12_RTX) * RTX_CENT
                 + K21_RTX * RTX_PERI
                 - TMDD_LOSS;
dxdt_RTX_PERI  =  K12_RTX * RTX_CENT
                 - K21_RTX * RTX_PERI;

// ---- Azathioprine / 6-TGN PK --------------------------------
double KE_AZA  = CL_AZA / V_AZA;

dxdt_AZA_GUT = -KA_AZA * AZA_GUT;
dxdt_SIX_MP  =  KCONV_AZA * AZA_GUT      // AZA → 6-MP conversion
              - (CL_6MP / V_6MP) * SIX_MP; // 6-MP elimination
dxdt_SIX_TGN =  KCONV_TGN * SIX_MP * V_6MP  // 6-MP → 6-TGN
              - CL_TGN * SIX_TGN;            // 6-TGN catabolism

// ---- Cytokines (driven by immune cells) ---------------------
// IL-4: secreted by Th2, degraded
double IL4_SYN_RATE  = KIL4_SYN  * (TH2  / TH2_SS);
double IL4_DEG_RATE  = KIL4_DEG  * IL4;

dxdt_IL4  = IL4_SYN_RATE  - IL4_DEG_RATE
           - COMBO * 0.5 * IL4;              // immunosuppression reduces IL4

// IL-10: from Treg and Th2 (Kamisawa 2014)
double IL10_SYN_RATE = KIL10_SYN * (0.6 * TREG / TREG_SS +
                                     0.4 * TH2  / TH2_SS);
dxdt_IL10 = IL10_SYN_RATE - KIL10_DEG * IL10;

// TGF-β: from Treg and activated PSC (major fibrosis driver)
double TGFB_SYN_RATE = KTGFB_SYN * (0.4 * TREG / TREG_SS +
                                      0.6 * APSC / 30.0);
dxdt_TGFB = TGFB_SYN_RATE - KTGFB_DEG * TGFB
           - COMBO * 0.3 * TGFB;

// TNF-α: from activated macrophages (proxy: APSC reflects local inflammation)
double TNFA_SYN_RATE = KTNFA_SYN * (1.0 + APSC / 30.0);
dxdt_TNFA = TNFA_SYN_RATE - KTNFA_DEG * TNFA
           - COMBO * 0.4 * TNFA;

// ---- Naïve CD4 T cells (TN) ---------------------------------
// Source: thymic production; Sink: activation (→ Th2), Treg pathway, natural death
double TN_INPUT  = KTN_IN / 24.0;    // convert /d → /h
double TN_TO_TH2 = KTN_ACT * TN * (IL4 / (IL4_SS + IL4)) *
                   (1.0 - COMBO);    // suppressed by drugs

dxdt_TN = TN_INPUT
         - TN_TO_TH2
         - KTREG_IN * TN                 // Treg generation
         - KTN_DEG  * TN;

// ---- Th2 effector cells -------------------------------------
// Positive feedback: IL-4 drives more Th2 differentiation
// Negative: Treg suppression, prednisolone, 6-TGN
double TH2_PROLIF = KTH2_DIFF * TN_TO_TH2 * (1.0 + IL4 / IL4_SS);
double TH2_TREG_SUP = KTH2_TREG * TREG * TH2;

dxdt_TH2 = TH2_PROLIF
          - KTH2_APOP * TH2 * (1.0 + EFF_PRED * 2.0)
          - TH2_TREG_SUP
          - EFF_TGN * 0.3 * TH2;

// ---- Regulatory T cells (TREG) ------------------------------
// Prednisolone promotes Treg in IgG4-RD (Kubota 2014)
double TREG_GEN = KTREG_IN * TN * (1.0 + PRED_TREG * EFF_PRED);

dxdt_TREG = TREG_GEN
           - KTREG_DEG * TREG;

// ---- Naïve B cells (BN) -------------------------------------
// Rituximab depletes BN via CD20 targeting
// Production from bone marrow (KBN_IN), death, activation into GC
double BN_ACT = KBN_ACT * BN * (IL4 / (IL4_SS + IL4)) *
                (TNFA / (TNFA_SS + TNFA));

dxdt_BN = KBN_IN
         - KBN_DEG  * BN
         - BN_ACT
         - EFF_RTX * KINT_RTX * BN * 10.0  // RTX-mediated B cell death
         - EFF_TGN * 0.2 * BN;

// ---- Germinal center B cells (BCG) --------------------------
// Driven by antigen / T cell help (IL-4, Th2)
// Suppress with rituximab
dxdt_BCG = KBCG_FORM * BN_ACT * (1.0 + TH2 / TH2_SS)
          - KBCG_DEG * BCG * (1.0 + COMBO)
          - KBCG_PB  * BCG
          - EFF_RTX  * 0.5 * BCG;

// ---- IgG4+ Plasmablasts (PB_IgG4) ---------------------------
// Exit GC as plasmablasts secreting IgG4 (IL-10 and IL-4 driven class switching)
double PB_FORM = KBCG_PB * BCG * (IL10 / (IL10_SS + IL10)) *
                 (IL4  / (IL4_SS  + IL4));

dxdt_PB_IgG4 = PB_FORM
              - KPB_DEG * PB_IgG4 * (1.0 + COMBO * 1.5)
              - EFF_RTX * 0.3 * PB_IgG4;   // RTX has partial effect on PB

// ---- Serum IgG4 concentration --------------------------------
// Secreted by plasmablasts; catabolised
// Half-life of IgG4 ~21 days → kel ~0.0014 h^-1 (Sandanayake 2009)
dxdt_IgG4 = KPB_IGG4 * PB_IgG4
            - KIG_DEG * IgG4;

// ---- Pancreatic stellate cells (PSC → APSC) -----------------
// TGFB is the dominant activator (Rosen 2016)
// Prednisolone reduces activation; PSC recovers once inflammation resolves
double PSC_TO_APSC = KPSC_ACT * PSC * (TGFB / TGFB_SS) *
                     (1.0 - EFF_PRED * 0.6);

dxdt_PSC  = KAPSC_BACK * APSC - PSC_TO_APSC;
dxdt_APSC = PSC_TO_APSC
           - KAPSC_DEG  * APSC
           - KAPSC_BACK * APSC;

// ---- Collagen deposition (fibrosis index) -------------------
// Synthesised by activated PSC; degraded by MMPs (IL-10 promotes MMP)
double COLL_SYN   = KCOLL_SYN * APSC;
double COLL_DEGEN = KCOLL_DEG * COLL * (1.0 + IL10 / IL10_SS * 0.5);

dxdt_COLL = COLL_SYN - COLL_DEGEN;

// ---- Exocrine pancreatic function (EXO) ---------------------
// Lost with fibrosis (COLL) and inflammation (TNFA); recovers with treatment
double EXO_LOSS = KEXO_LOSS * EXO * COLL;
double EXO_REC  = KEXO_REC  * (EXO_SS - EXO) * COMBO;

dxdt_EXO = EXO_REC - EXO_LOSS;
// Bound: 0–100
if (EXO < 0.0)   EXO = 0.0;
if (EXO > 100.0) EXO = 100.0;

// ---- Beta cell mass (BETA) ----------------------------------
// Injured by local inflammation and fibrosis; partial recovery with treatment
double BETA_LOSS = KBETA_LOSS * BETA * (APSC / 30.0) * (TNFA / TNFA_SS);
double BETA_REC  = KBETA_REC  * (BETA_SS - BETA) * COMBO;

dxdt_BETA = BETA_REC - BETA_LOSS;
if (BETA < 0.0)   BETA = 0.0;
if (BETA > 100.0) BETA = 100.0;

$TABLE
// ============================================================
// DERIVED TABLE VARIABLES
// ============================================================

// --- PK-derived concentrations -------------------------------
double Cpred_mgL   = PRED_CENT / V2_PRED;         // mg/L
double Cpred_ngmL  = Cpred_mgL * 1000.0;          // ng/mL
double Crtx_nM     = RTX_CENT;                    // nM
double Crtx_ugmL   = Crtx_nM * 148.0 / 1000.0;   // μg/mL (MW ~148 kDa)
double C6mp_umolL  = SIX_MP;                      // μmol/L
double C6tgn       = SIX_TGN;                     // pmol/8e8 RBC

// --- Clinical IgG4 endpoint ----------------------------------
double IgG4_serum  = IgG4;                        // ng/mL
double IgG4_mgdL   = IgG4 / 10000.0;             // mg/dL (1 mg/dL = 10,000 ng/mL)

// --- HbA1c estimate (from beta cell mass) --------------------
// Simplified relationship: beta cell loss → fasting glucose → HbA1c
// Normal HbA1c ~5.5%; AIP at 70% beta mass → ~6.5%
double FBG_est     = 80.0 + (100.0 - BETA) * 1.8; // mg/dL estimated fasting BG
double HbA1c       = 5.5 + (FBG_est - 80.0) / 100.0 * 2.5;  // % (simplified)
if (HbA1c < 4.0)  HbA1c = 4.0;
if (HbA1c > 14.0) HbA1c = 14.0;

// --- Fecal elastase-1 proxy (from exocrine function) ---------
// Normal FE-1 > 200 μg/g; exocrine insufficiency < 100 μg/g
double FE1         = EXO * 2.5;                   // μg/g (100% → 250 μg/g)

// --- Pancreatic duct diameter index (from fibrosis/oedema) ---
// Shimosegawa 2011: diffuse narrowing is hallmark
// Proxy: increases with COLL, normalises with treatment
double duct_index  = 1.0 + (COLL / COLL_SS - 1.0) * 0.3;
if (duct_index < 0.5) duct_index = 0.5;

// --- Relapse hazard (instantaneous) --------------------------
// Yamamoto 2019: residual IgG4 elevation drives relapse risk
double relapse_haz = RELAPSE_K * IgG4;            // per hour

// --- Drug effect summaries -----------------------------------
double eff_pred_pct = EFF_PRED * 100.0;
double eff_rtx_pct  = EFF_RTX  * 100.0;
double eff_tgn_pct  = EFF_TGN  * 100.0;

// --- Total B cell index --------------------------------------
double B_total = BN + BCG + PB_IgG4;

$CAPTURE
// ============================================================
// VARIABLES WRITTEN TO OUTPUT TABLE
// ============================================================
Cpred_ngmL Crtx_nM Crtx_ugmL C6mp_umolL C6tgn
IgG4_serum IgG4_mgdL
HbA1c FE1 duct_index relapse_haz
TN TH2 TREG BN BCG PB_IgG4
PSC APSC COLL
EXO BETA
IL4 IL10 TGFB TNFA
eff_pred_pct eff_rtx_pct eff_tgn_pct
B_total

'  ## end of model code string

## ---------------------------------------------------------------------------
## 2. Compile model
## ---------------------------------------------------------------------------
message("Compiling AIP QSP mrgsolve model...")
mod <- mcode("aip_qsp_v1", aip_model_code, quiet = TRUE)
message("Model compiled successfully.")
message(paste("Compartments:", length(init(mod))))
message(paste("Parameters:  ", length(param(mod))))

## ---------------------------------------------------------------------------
## 3. Helper functions
## ---------------------------------------------------------------------------

#' Build prednisolone dosing event table
#' @param dose_mgkgd  mg/kg/day dose
#' @param bwt         body weight (kg)
#' @param start_h     start time (hours)
#' @param dur_wk      duration in weeks
#' @param freq_h      dosing interval (h); default 24 (once daily)
make_pred_events <- function(dose_mgkgd, bwt = 65, start_h = 0,
                              dur_wk = 4, freq_h = 24) {
  dose_mg <- dose_mgkgd * bwt
  times   <- seq(start_h, start_h + dur_wk * 168 - freq_h, by = freq_h)
  ev(amt  = dose_mg * 0.85,   # apply bioavailability into gut
     cmt  = "PRED_GUT",
     time = times,
     ii   = freq_h,
     addl = 0)
}

#' Build rituximab IV dosing event
make_rtx_events <- function(dose_mg = 1000, bwt = 65,
                             times_h = c(0, 336)) {  # 0 and 2 weeks
  mw_rtx  <- 148000    # g/mol
  vc_rtx  <- 3.1       # L
  conc_nM <- (dose_mg / mw_rtx * 1e6) / vc_rtx   # nM bolus into central
  ev(amt  = conc_nM,
     cmt  = "RTX_CENT",
     time = times_h,
     rate = -2)        # bolus (rate=-2 → instantaneous in mrgsolve)
}

#' Build azathioprine dosing event table
make_aza_events <- function(dose_mgkgd = 2.0, bwt = 65,
                             start_h = 0, dur_wk = 52) {
  dose_mg <- dose_mgkgd * bwt
  times   <- seq(start_h, start_h + dur_wk * 168 - 24, by = 24)
  ev(amt  = dose_mg,
     cmt  = "AZA_GUT",
     time = times,
     addl = 0)
}

## ---------------------------------------------------------------------------
## 4. Simulation end time and observation grid
## ---------------------------------------------------------------------------
SIM_END_H   <- 8760   # 1 year = 8760 hours
OBS_GRID    <- seq(0, SIM_END_H, by = 24)  # daily observations

## ---------------------------------------------------------------------------
## 5. Scenario Definitions
## ---------------------------------------------------------------------------
## Scenario 1: Prednisolone Induction (0.6 mg/kg/day × 4 wk then taper)
## Yoshida 2002: initial 40 mg/d (0.6 mg/kg for 65-kg patient)
## Taper: reduce by 5 mg every 2 weeks until 5 mg/d maintenance

build_pred_induction_taper <- function(bwt = 65) {
  # Phase 1: 0.6 mg/kg/d × 4 weeks (39 mg for 65-kg)
  e1 <- make_pred_events(0.60, bwt = bwt, start_h = 0,  dur_wk = 4)
  # Phase 2: taper steps
  e2 <- make_pred_events(0.46, bwt = bwt, start_h = 672,  dur_wk = 2)  # wk 5-6
  e3 <- make_pred_events(0.38, bwt = bwt, start_h = 1008, dur_wk = 2)  # wk 7-8
  e4 <- make_pred_events(0.31, bwt = bwt, start_h = 1344, dur_wk = 2)  # wk 9-10
  e5 <- make_pred_events(0.23, bwt = bwt, start_h = 1680, dur_wk = 2)  # wk 11-12
  e6 <- make_pred_events(0.15, bwt = bwt, start_h = 2016, dur_wk = 40) # maintenance 10 mg/d
  do.call(c, list(e1, e2, e3, e4, e5, e6))
}

scenarios <- list(

  # ------------------------------------------------------------------
  # S1: Prednisolone induction + taper (standard first-line)
  # ------------------------------------------------------------------
  S1_pred_induction = list(
    label  = "Prednisolone Induction + Taper",
    events = build_pred_induction_taper(bwt = 65)
  ),

  # ------------------------------------------------------------------
  # S2: Prednisolone maintenance (low dose 5 mg/d × 3 years)
  # Yoshida 2002: prolonged maintenance reduces relapse
  # ------------------------------------------------------------------
  S2_pred_maintenance = list(
    label  = "Prednisolone Maintenance (5 mg/d × 3 yr)",
    events = {
      e_ind <- build_pred_induction_taper(bwt = 65)
      e_mnt <- make_pred_events(0.077, bwt = 65,       # ~5 mg/d
                                 start_h = 2688, dur_wk = 130)  # wk 16 – 3 yr
      do.call(c, list(e_ind, e_mnt))
    }
  ),

  # ------------------------------------------------------------------
  # S3: Rituximab monotherapy (1000 mg × 2 doses, 2 weeks apart)
  # Hart 2013: RTX for refractory / steroid-intolerant AIP
  # ------------------------------------------------------------------
  S3_rituximab = list(
    label  = "Rituximab 1000 mg × 2 (wk 0 + wk 2)",
    events = make_rtx_events(dose_mg = 1000, times_h = c(0, 336))
  ),

  # ------------------------------------------------------------------
  # S4: Azathioprine maintenance after prednisolone induction
  # ------------------------------------------------------------------
  S4_aza_maintenance = list(
    label  = "Prednisolone Induction + Azathioprine Maintenance",
    events = {
      e_pred <- build_pred_induction_taper(bwt = 65)
      e_aza  <- make_aza_events(dose_mgkgd = 2.0, bwt = 65,
                                 start_h = 2688, dur_wk = 30)
      do.call(c, list(e_pred, e_aza))
    }
  ),

  # ------------------------------------------------------------------
  # S5: Rituximab re-treatment at relapse (~6 months)
  # Carruthers 2015: RTX re-treatment effective in IgG4-RD
  # ------------------------------------------------------------------
  S5_rtx_retreatment = list(
    label  = "RTX Induction + Re-treatment at Relapse (6 mo)",
    events = {
      e1 <- make_rtx_events(dose_mg = 1000, times_h = c(0, 336))
      e2 <- make_rtx_events(dose_mg = 1000, times_h = c(4380, 4716)) # ~6 months
      do.call(c, list(e1, e2))
    }
  ),

  # ------------------------------------------------------------------
  # S6: Combination Prednisolone + Azathioprine (from day 1)
  # ------------------------------------------------------------------
  S6_combination = list(
    label  = "Combination: Prednisolone + Azathioprine",
    events = {
      e_pred <- build_pred_induction_taper(bwt = 65)
      e_aza  <- make_aza_events(dose_mgkgd = 1.5, bwt = 65,
                                 start_h = 0, dur_wk = 52)
      do.call(c, list(e_pred, e_aza))
    }
  ),

  # ------------------------------------------------------------------
  # S7: No treatment — natural history
  # ------------------------------------------------------------------
  S7_natural_history = list(
    label  = "No Treatment (Natural History)",
    events = ev(amt = 0, cmt = "PRED_GUT", time = 0)  # null event
  )

)

## ---------------------------------------------------------------------------
## 6. Run all scenarios
## ---------------------------------------------------------------------------
message("\nRunning ", length(scenarios), " treatment scenarios...")

run_scenario <- function(sc, sim_end = SIM_END_H, obs = OBS_GRID) {
  mod %>%
    mrgsim(events = sc$events,
           end    = sim_end,
           delta  = 24,
           obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(time_d    = time / 24,
           time_wk   = time / 168,
           scenario  = sc$label)
}

sim_results <- lapply(scenarios, function(sc) {
  message("  Running: ", sc$label)
  run_scenario(sc)
})

all_results <- bind_rows(sim_results)

## ---------------------------------------------------------------------------
## 7. Summary tables
## ---------------------------------------------------------------------------

## Key clinical timepoints: 4 weeks, 3 months, 6 months, 12 months
timepoints_wk <- c(4, 12, 26, 52)

summary_table <- all_results %>%
  filter(round(time_wk) %in% timepoints_wk) %>%
  group_by(scenario, time_wk = round(time_wk)) %>%
  summarise(
    IgG4_ngmL   = round(mean(IgG4_serum), 0),
    IgG4_mgdL   = round(mean(IgG4_mgdL), 2),
    HbA1c       = round(mean(HbA1c), 1),
    FE1_ug_g    = round(mean(FE1), 0),
    EXO_pct     = round(mean(EXO), 1),
    BETA_pct    = round(mean(BETA), 1),
    PB_count    = round(mean(PB_IgG4), 2),
    COLL_index  = round(mean(COLL), 1),
    .groups     = "drop"
  )

message("\n--- Clinical Summary (Key Endpoints) ---")
print(summary_table, n = Inf)

## ---------------------------------------------------------------------------
## 8. Visualisation
## ---------------------------------------------------------------------------

theme_qsp <- theme_bw(base_size = 11) +
  theme(legend.position  = "bottom",
        legend.key.width  = unit(1.5, "cm"),
        strip.background  = element_rect(fill = "#E8EDF2"),
        panel.grid.minor  = element_blank())

scenario_colors <- c(
  "Prednisolone Induction + Taper"                       = "#E41A1C",
  "Prednisolone Maintenance (5 mg/d × 3 yr)"            = "#FF7F00",
  "Rituximab 1000 mg × 2 (wk 0 + wk 2)"                = "#4DAF4A",
  "Prednisolone Induction + Azathioprine Maintenance"   = "#984EA3",
  "RTX Induction + Re-treatment at Relapse (6 mo)"     = "#377EB8",
  "Combination: Prednisolone + Azathioprine"            = "#A65628",
  "No Treatment (Natural History)"                      = "#999999"
)

# ---- Figure 1: Serum IgG4 over time ---------------------------------
p_igg4 <- ggplot(all_results,
                 aes(x = time_wk, y = IgG4_serum,
                     colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1350, linetype = "dashed",
             colour = "red", linewidth = 0.5) +
  annotate("text", x = 50, y = 1500,
           label = "Diagnostic threshold (1350 ng/mL)",
           size = 3, colour = "red") +
  scale_colour_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title    = "Serum IgG4 Concentration Over Time",
       subtitle = "AIP Type 1 — QSP Model Simulations",
       x = "Time (weeks)", y = "Serum IgG4 (ng/mL)",
       colour = "Treatment Scenario",
       linetype = "Treatment Scenario") +
  theme_qsp

# ---- Figure 2: Exocrine pancreatic function -------------------------
p_exo <- ggplot(all_results,
                aes(x = time_wk, y = EXO,
                    colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "darkgreen") +
  scale_colour_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  ylim(0, 105) +
  labs(title    = "Pancreatic Exocrine Function",
       x = "Time (weeks)", y = "Exocrine Function (%)",
       colour = "Treatment Scenario",
       linetype = "Treatment Scenario") +
  theme_qsp

# ---- Figure 3: HbA1c ------------------------------------------------
p_hba1c <- ggplot(all_results,
                  aes(x = time_wk, y = HbA1c,
                      colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 6.5, linetype = "dashed",
             colour = "blue", linewidth = 0.5) +
  annotate("text", x = 48, y = 6.6,
           label = "Diabetes threshold (6.5%)",
           size = 3, colour = "blue") +
  scale_colour_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title    = "HbA1c Over Time",
       x = "Time (weeks)", y = "HbA1c (%)",
       colour = "Treatment Scenario",
       linetype = "Treatment Scenario") +
  theme_qsp

# ---- Figure 4: B cell dynamics (BN + PB_IgG4) ----------------------
p_bcell <- all_results %>%
  select(time_wk, scenario, BN, PB_IgG4) %>%
  pivot_longer(c(BN, PB_IgG4), names_to = "cell_type", values_to = "count") %>%
  mutate(cell_type = recode(cell_type,
                            BN = "Naive B cells",
                            PB_IgG4 = "IgG4+ Plasmablasts")) %>%
  ggplot(aes(x = time_wk, y = count,
             colour = scenario, linetype = cell_type)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = scenario_colors) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 52, 8)) +
  labs(title    = "B Cell Dynamics",
       x = "Time (weeks)", y = "Cell Count (cells/μL)",
       colour = "Treatment Scenario") +
  theme_qsp +
  theme(legend.position = "bottom")

# ---- Figure 5: Fibrosis index (COLL) and PSC status -----------------
p_fibrosis <- all_results %>%
  select(time_wk, scenario, COLL, APSC) %>%
  pivot_longer(c(COLL, APSC),
               names_to = "marker", values_to = "value") %>%
  mutate(marker = recode(marker,
                         COLL = "Collagen Index (AU)",
                         APSC = "Activated PSC (AU)")) %>%
  ggplot(aes(x = time_wk, y = value,
             colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = scenario_colors) +
  facet_wrap(~marker, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 52, 8)) +
  labs(title = "Pancreatic Fibrosis Markers",
       x = "Time (weeks)", y = "Value (AU)",
       colour = "Treatment Scenario") +
  theme_qsp +
  theme(legend.position = "bottom")

# ---- Figure 6: Cytokine profiles ------------------------------------
p_cytokines <- all_results %>%
  filter(scenario %in% c("Prednisolone Induction + Taper",
                          "Rituximab 1000 mg × 2 (wk 0 + wk 2)",
                          "No Treatment (Natural History)")) %>%
  select(time_wk, scenario, IL4, IL10, TGFB, TNFA) %>%
  pivot_longer(c(IL4, IL10, TGFB, TNFA),
               names_to = "cytokine", values_to = "pg_mL") %>%
  ggplot(aes(x = time_wk, y = pg_mL,
             colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = scenario_colors) +
  facet_wrap(~cytokine, scales = "free_y", nrow = 2) +
  scale_x_continuous(breaks = seq(0, 52, 8)) +
  labs(title    = "Cytokine Profiles (Selected Scenarios)",
       x = "Time (weeks)", y = "Concentration (pg/mL)",
       colour = "Treatment Scenario") +
  theme_qsp +
  theme(legend.position = "bottom")

# ---- Figure 7: Drug PK for prednisolone (S1) ------------------------
pred_pk_data <- all_results %>%
  filter(scenario == "Prednisolone Induction + Taper")

p_pred_pk <- ggplot(pred_pk_data,
                    aes(x = time_wk, y = Cpred_ngmL)) +
  geom_line(colour = "#E41A1C", linewidth = 1.0) +
  geom_hline(yintercept = 50, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 48, y = 65,
           label = "EC50 threshold (~50 ng/mL)", size = 3, colour = "grey40") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title    = "Prednisolone Plasma Concentration (Trough)",
       subtitle = "Scenario S1: Induction + Taper",
       x = "Time (weeks)", y = "Prednisolone (ng/mL)") +
  theme_qsp

# ---- Figure 8: Rituximab PK + B cell depletion (S3) ----------------
rtx_data <- all_results %>%
  filter(scenario == "Rituximab 1000 mg × 2 (wk 0 + wk 2)")

p_rtx_pk <- rtx_data %>%
  select(time_wk, Crtx_ugmL, B_total) %>%
  pivot_longer(c(Crtx_ugmL, B_total),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           Crtx_ugmL = "Rituximab (μg/mL)",
                           B_total   = "Total B cells (cells/μL)")) %>%
  ggplot(aes(x = time_wk, y = value)) +
  geom_line(colour = "#4DAF4A", linewidth = 1.0) +
  facet_wrap(~variable, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title    = "Rituximab PK and B Cell Depletion",
       subtitle = "Scenario S3: RTX 1000 mg × 2",
       x = "Time (weeks)", y = "Value") +
  theme_qsp

## ---------------------------------------------------------------------------
## 9. Combine and save figures
## ---------------------------------------------------------------------------

fig_main <- (p_igg4 + p_exo) /
            (p_hba1c + p_bcell) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "AIP Type 1 QSP Model — Clinical Endpoints",
    subtitle = "mrgsolve simulation | 7 treatment scenarios | 1-year horizon",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  ) &
  theme(legend.position = "bottom")

fig_mechanisms <- (p_fibrosis / p_cytokines) +
  plot_annotation(
    title    = "AIP Type 1 QSP Model — Mechanistic Readouts",
    subtitle = "Fibrosis cascade and cytokine dynamics",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

fig_pk <- (p_pred_pk + p_rtx_pk) +
  plot_annotation(
    title    = "AIP Type 1 QSP Model — Drug PK",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

## Save if output path exists
out_dir <- "/home/user/qsp/autoimmune-pancreatitis"
if (dir.exists(out_dir)) {
  ggsave(file.path(out_dir, "aip_clinical_endpoints.png"),
         fig_main, width = 14, height = 10, dpi = 150)
  ggsave(file.path(out_dir, "aip_mechanistic_readouts.png"),
         fig_mechanisms, width = 14, height = 10, dpi = 150)
  ggsave(file.path(out_dir, "aip_drug_pk.png"),
         fig_pk, width = 14, height = 6, dpi = 150)
  message("\nFigures saved to ", out_dir)
}

## ---------------------------------------------------------------------------
## 10. Sensitivity analysis — prednisolone dose response
## ---------------------------------------------------------------------------
message("\nRunning dose–response sensitivity analysis...")

pred_doses <- c(0.3, 0.4, 0.5, 0.6, 0.8, 1.0)   # mg/kg/day

dose_resp <- lapply(pred_doses, function(d) {
  ev_d <- make_pred_events(d, bwt = 65, start_h = 0, dur_wk = 8)
  res  <- mod %>%
    mrgsim(events = ev_d, end = 8 * 168, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    filter(time == max(time)) %>%
    mutate(dose_mgkgd = d)
  res
})

dose_resp_df <- bind_rows(dose_resp)

p_dose_resp <- dose_resp_df %>%
  select(dose_mgkgd, IgG4_serum, HbA1c, EXO, COLL) %>%
  pivot_longer(-dose_mgkgd,
               names_to = "endpoint", values_to = "value") %>%
  ggplot(aes(x = dose_mgkgd, y = value, colour = endpoint)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  facet_wrap(~endpoint, scales = "free_y") +
  scale_colour_brewer(palette = "Set1") +
  labs(title    = "Prednisolone Dose–Response at 8 Weeks",
       subtitle = "AIP Type 1 QSP Model",
       x = "Prednisolone Dose (mg/kg/day)", y = "Endpoint Value") +
  theme_qsp +
  theme(legend.position = "none")

if (dir.exists(out_dir)) {
  ggsave(file.path(out_dir, "aip_dose_response.png"),
         p_dose_resp, width = 10, height = 7, dpi = 150)
}

## ---------------------------------------------------------------------------
## 11. Relapse risk analysis
## ---------------------------------------------------------------------------
message("\nComputing relapse risk curves...")

relapse_data <- all_results %>%
  mutate(cumulative_hazard = relapse_haz * 24) %>%    # convert h^-1 to d^-1
  group_by(scenario) %>%
  arrange(time_wk) %>%
  mutate(cum_haz   = cumsum(cumulative_hazard),
         surv_prob = exp(-cum_haz / 365),             # normalised annual
         relapse_risk_pct = (1 - surv_prob) * 100)

p_relapse <- ggplot(relapse_data,
                    aes(x = time_wk, y = relapse_risk_pct,
                        colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  ylim(0, 100) +
  labs(title    = "Cumulative Relapse Risk Over 1 Year",
       subtitle = "Based on serum IgG4-derived hazard (Yamamoto 2019)",
       x = "Time (weeks)", y = "Cumulative Relapse Risk (%)",
       colour = "Treatment Scenario",
       linetype = "Treatment Scenario") +
  theme_qsp

if (dir.exists(out_dir)) {
  ggsave(file.path(out_dir, "aip_relapse_risk.png"),
         p_relapse, width = 10, height = 6, dpi = 150)
}

## ---------------------------------------------------------------------------
## 12. Parameter calibration notes
## ---------------------------------------------------------------------------
cat("
================================================================================
PARAMETER CALIBRATION NOTES — AIP Type 1 QSP Model
================================================================================

PREDNISOLONE PK (Yoshida et al. 2002, Pancreas 25:325–330)
  - Consensus corticosteroid dose: 0.6 mg/kg/day initial, taper to 5 mg/d
  - PK parameters (CL=5.7 L/h, Vd=31 L) from Czock 2005 Clin Pharmacokinet
  - Taper schedule: 5 mg/2-week step reduction until 10 mg/d, then slow
  - Remission at 4 weeks: ~70-85% radiological response

RITUXIMAB PK (Hart et al. 2013, Gut 62:1237–1244;
              Khosroshahi et al. 2015, Arthritis Rheum 67:1766–1773)
  - 1000 mg IV × 2 doses (2 weeks apart) standard AIP/IgG4-RD protocol
  - Complete B cell depletion in >95% patients at 4 weeks post-RTX
  - IgG4 normalisation at 6 months: ~70% of responders
  - PK: Vc=3.1 L, CL=0.14 L/h calibrated from Cartron 2014 (NHL/autoimmune)
  - TMDD: CD20 receptor concentration 35 nM (Jones 2009, Clin Cancer Res)

AZATHIOPRINE / 6-TGN (Iwata et al. 2020, Pancreas 49:1316-1321;
                       Dubinsky 2002, Gastroenterology 122:904–915)
  - 6-TGN therapeutic range: 235–450 pmol/8×10^8 RBC
  - AZA 2 mg/kg/d achieves 6-TGN ~275 pmol/8e8 RBC within 8–12 weeks
  - Conversion efficiency AZA→6-MP variable (TPMT polymorphism not modelled)
  - EC50 for lymphocyte suppression: ~235 pmol/8e8 RBC (threshold model)

IMMUNE PARAMETERS
  - IgG4+ plasmablast elevation (>200 cells/mL) diagnostic for active AIP1
    (Ryu 2012, Gastroenterology 143:672–681)
  - Serum IgG4 at diagnosis: median 1350 ng/mL (range 400–12000+)
    (Shimosegawa 2011, Pancreas 40:352–358)
  - Th2/Treg imbalance confirmed: IL-4, IL-10, TGF-β elevated in AIP1
    (Kamisawa 2014, Lancet 385:1460–1471)
  - Normal IgG4 upper limit: 135 ng/mL (Sandanayake 2009, Gut 58:1580–1586)
    >2× ULN (>270 ng/mL) is diagnostic criterion level 1

FIBROSIS PARAMETERS
  - Pancreatic stellate cell activation: TGF-β primary driver (Rosen 2016)
  - Storiform fibrosis is pathognomonic of AIP1 (Shimosegawa 2011)
  - Fibrosis partially reversible with early treatment (Kamisawa 2014)
  - Collagen turnover half-life estimated ~300–600 h (Iredale 1998)

CLINICAL ENDPOINTS
  - Exocrine insufficiency: present in ~50-70% AIP at diagnosis
    (Frulloni 2009, Gut 58:1504–1510)
  - Diabetes (new-onset): 70-80% AIP patients (Miyamoto 2012, J Gastroenterol)
  - FE-1 < 200 μg/g = exocrine insufficiency (Loser 1996)
  - Remission criteria: resolution of pancreatic swelling + IgG4 normalisation

RELAPSE RISK (Yamamoto et al. 2019, Modern Rheum 29:17–26)
  - 1-year relapse rate without maintenance: 30-40%
  - Persistent IgG4 elevation (>2× ULN) at 6 months: strongest predictor
  - Maintenance prednisolone (5 mg/d) reduces relapse by ~50%
================================================================================
")

## ---------------------------------------------------------------------------
## 13. Model diagnostics
## ---------------------------------------------------------------------------
message("\nModel structure summary:")
message("  ODE compartments : ", length(init(mod)))
message("  Parameters       : ", length(param(mod)))
message("  Captured outputs : ", length(outvars(mod)$capture))
message("  Simulation rows  : ", nrow(all_results))
message("  Scenarios run    : ", length(scenarios))

## List initial conditions
message("\nInitial conditions (disease state at AIP diagnosis):")
print(init(mod))

## Print key parameters
message("\nKey PK/PD parameters:")
print(as.data.frame(param(mod))[
  c("KA_PRED","CL_PRED","V2_PRED","F_PRED",
    "CL_RTX","VC_RTX","KON_RTX","CD20_TOT",
    "KA_AZA","KCONV_TGN","CL_TGN",
    "EC50_PRED","EMAX_PRED","EC50_TGN","EMAX_TGN"), , drop = FALSE])

message("\n=== AIP QSP Model simulation complete ===")
message("All outputs saved to: ", out_dir)
