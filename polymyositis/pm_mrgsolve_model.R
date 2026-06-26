## =============================================================================
## Polymyositis (PM) QSP Model — mrgsolve ODE
## =============================================================================
## Disease: Polymyositis (IIM subtype, CD8+ T-cell mediated inflammatory myopathy)
## Author : QSP Library (CCR)
## Date   : 2026-06-19
##
## Model Structure (26 compartments):
##   PK: Prednisone (3-CMT), Methotrexate (2-CMT + polyglutamate),
##       Azathioprine/6-TGN (2-CMT), Rituximab (2-CMT TMDD),
##       IVIG (1-CMT FcRn), Baricitinib/JAKi (2-CMT)
##   PD: CD8+ T cells, CD4+/Th1 cells, B cells, Cytokine network,
##       Muscle Inflammation State, CK dynamics, MMT-8 score
##
## Key References:
##   Lundberg et al. Nat Rev Dis Primers 2021; 7:9
##   Zong et al. Ann Rheum Dis 2021; 80:1293
##   Aggarwal et al. Semin Arthritis Rheum 2020; 50:1031
##   Oddis & Aggarwal, Clin Immunol 2018; 192:64
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─── Model Definition ────────────────────────────────────────────────────────
pm_model_code <- '
$PROB Polymyositis QSP Model — CD8+ T-cell Mediated Myopathy

$PARAM
// ── Prednisone PK ──────────────────────────────────────────────────────
  PRED_KA   = 1.5,     // Absorption rate constant (h-1)
  PRED_F    = 0.82,    // Oral bioavailability
  PRED_CL   = 15.0,    // Clearance (L/h)
  PRED_V1   = 45.0,    // Central volume (L)
  PRED_Q    = 8.0,     // Inter-compartmental clearance (L/h)
  PRED_V2   = 90.0,    // Peripheral volume (L)
  KHB       = 0.7,     // Prednisolone conversion (hepatic activation fraction)

// ── Methotrexate PK ─────────────────────────────────────────────────────
  MTX_KA    = 0.9,     // Absorption (h-1) — weekly SC/PO
  MTX_F     = 0.72,    // Bioavailability
  MTX_CL    = 5.2,     // Clearance (L/h)
  MTX_V1    = 28.0,    // Central volume (L)
  MTX_KPG   = 0.05,    // Polyglutamation rate (h-1)
  MTX_KPGOUT= 0.003,   // Polyglutamate elimination (h-1; long t½ ~7d)

// ── AZA/6-TGN PK ─────────────────────────────────────────────────────
  AZA_KA    = 0.8,     // AZA absorption (h-1)
  AZA_F     = 0.60,    // AZA bioavailability
  AZA_K12   = 0.25,    // AZA → 6-MP conversion (h-1)
  AZA_6MP_CL= 2.5,     // 6-MP clearance (L/h)
  AZA_V2    = 40.0,    // 6-MP/6-TGN volume (L)
  TPMT_act  = 0.6,     // TPMT activity (normalized; 0=poor, 1=normal)

// ── Rituximab PK (TMDD) ──────────────────────────────────────────────
  RTX_CL    = 0.008,   // Linear clearance (L/h)
  RTX_V1    = 3.5,     // Central volume (L) — ~3.1-4 L reported
  RTX_Q     = 0.025,   // Inter-compartmental clearance (L/h)
  RTX_V2    = 5.0,     // Peripheral volume (L)
  RTX_KON   = 0.60,    // CD20 binding (nM-1 h-1)
  RTX_KOFF  = 0.01,    // CD20 unbinding (h-1)
  RTX_KINT  = 0.03,    // Internalization rate (h-1)
  RTX_CD20_0= 150.0,   // Baseline CD20 (nM)

// ── IVIG PK (1-CMT FcRn) ─────────────────────────────────────────────
  IVIG_CL0  = 0.006,   // Baseline IgG clearance (L/h)
  IVIG_V    = 4.5,     // Distribution volume (L)
  FCRN_CAP  = 500.0,   // FcRn capacity (arbitrary units)

// ── JAK inhibitor PK (Baricitinib) ──────────────────────────────────
  JAKI_KA   = 1.2,     // Absorption (h-1)
  JAKI_F    = 0.79,    // Bioavailability
  JAKI_CL   = 5.8,     // Clearance (L/h)
  JAKI_V1   = 76.0,    // Central volume (L)
  JAKI_Q    = 3.5,     // Inter-compt CL (L/h)
  JAKI_V2   = 40.0,    // Peripheral volume (L)

// ── CD8+ T cell dynamics ─────────────────────────────────────────────
  CD8_PROD  = 5.0,     // Naïve CD8+ production (cells/μL/d)
  CD8_KDEATH= 0.02,    // Naïve CD8+ death rate (d-1)
  CD8_KACT  = 0.15,    // CD8 activation rate (d-1, MHC-I dependent)
  CD8_KEXP  = 0.8,     // Effector expansion factor
  CD8_KMIG  = 0.20,    // Migration to muscle (d-1)
  CD8_KDEATH_E= 0.08,  // Effector CD8 death (d-1)

// ── CD4+ / Th1 dynamics ──────────────────────────────────────────────
  CD4_PROD  = 4.0,     // Th1 production (cells/μL/d)
  CD4_KACT  = 0.10,    // Activation rate (d-1)
  CD4_KDEATH= 0.05,    // Death rate (d-1)
  TH1_IFNg_PROD = 0.8, // IFN-γ production per Th1 cell (pg/mL/cell/d)

// ── B cell & Autoantibody dynamics ───────────────────────────────────
  B_PROD    = 2.0,     // Naïve B cell production (cells/μL/d)
  B_KACT    = 0.05,    // B activation rate (d-1)
  B_KDEATH  = 0.02,    // Naïve B death (d-1)
  PC_PROD_RATE= 0.10,  // Plasma cell production from B cells (d-1)
  PC_KDEATH = 0.005,   // Long-lived plasma cell death (d-1) — t½ ~6mo
  AB_PROD   = 0.5,     // Autoantibody production per PC (U/mL/cell/d)
  AB_CL     = 0.015,   // Autoantibody clearance (d-1) — t½ ~21d

// ── Cytokine dynamics ────────────────────────────────────────────────
  IFNG_PROD_0= 0.08,   // Basal IFN-γ production (pg/mL/d)
  IFNG_KD   = 0.30,    // IFN-γ clearance (d-1)
  TNF_PROD_0 = 0.06,   // Basal TNF-α production (pg/mL/d)
  TNF_KD    = 0.50,    // TNF-α clearance (d-1)
  IL6_PROD_0 = 0.10,   // Basal IL-6 production (pg/mL/d)
  IL6_KD    = 0.60,    // IL-6 clearance (d-1)
  IL1_PROD_0 = 0.04,   // Basal IL-1β production (pg/mL/d)
  IL1_KD    = 0.40,    // IL-1β clearance (d-1)

// ── Muscle Inflammation State ─────────────────────────────────────────
  MUS_INF_KMAX = 0.20, // Max muscle inflammation induction rate (d-1)
  MUS_INF_KD  = 0.05,  // Spontaneous resolution rate (d-1)
  MUS_INF_EC50 = 15.0, // EC50 for effector CD8 in muscle (cells/μL)
  MHC_IFNG_EC50= 2.0,  // EC50 for IFN-γ on MHC-I expression (pg/mL)
  MHCI_KMAX = 5.0,     // Max MHC-I upregulation (fold)

// ── CK dynamics ──────────────────────────────────────────────────────
  CK_BASE   = 120.0,   // Normal CK baseline (U/L)
  CK_KPROD  = 0.10,    // CK release rate driven by muscle damage (d-1)
  CK_KD     = 0.15,    // CK elimination from serum (d-1; t½~1d)
  CK_MAX_FACTOR = 400.0, // Max CK = 400× baseline in severe PM

// ── MMT-8 score dynamics ─────────────────────────────────────────────
  MMT_MAX   = 80.0,    // Maximum MMT-8 (full strength)
  MMT_KDECLINE=0.008,  // Rate of decline per unit muscle inflammation (d-1)
  MMT_KREC  = 0.015,   // Rate of recovery with treatment response (d-1)

// ── Steroid PD parameters ────────────────────────────────────────────
  PRED_IC50 = 2.0,     // Prednisolone IC50 for inflammation (ng/mL)
  PRED_IMAX = 0.90,    // Max inhibition by prednisolone
  PRED_HILL = 1.5,     // Hill coefficient

// ── MTX PD parameters ────────────────────────────────────────────────
  MTX_IC50  = 0.05,    // MTX polyglutamate IC50 for CD8 prolif (μg/L)
  MTX_IMAX  = 0.75,    // Max inhibition

// ── AZA/6-TGN PD parameters ──────────────────────────────────────────
  TGN_IC50  = 180.0,   // 6-TGN IC50 for lymphocyte proliferation (fmol/8×10^8 RBC)
  TGN_IMAX  = 0.70,    // Max inhibition

// ── RTX PD parameters ────────────────────────────────────────────────
  RTX_DEPL_EC50 = 0.5, // RTX-CD20 complex EC50 for B depletion (nM)
  RTX_DEPL_IMAX = 0.97,// Max B cell depletion

// ── IVIG PD parameters ───────────────────────────────────────────────
  IVIG_IC50 = 8.0,     // IVIG IC50 for macrophage FcR (g/L)
  IVIG_IMAX = 0.60,    // Max IVIG anti-inflammatory effect

// ── JAKi PD parameters ───────────────────────────────────────────────
  JAKI_IC50  = 15.0,   // JAKi IC50 for JAK1/2 (nM)
  JAKI_IMAX  = 0.85,   // Max JAKi effect on IFN-γ / IL-6 signaling

// ── Disease modifier ─────────────────────────────────────────────────
  PM_SEVERITY = 1.0    // 1=moderate, 2=severe; scales initial conditions

$CMT
// PK compartments (12)
  PRED_GUT    // Prednisone gut compartment
  PRED_C1     // Prednisone central (active prednisolone)
  PRED_C2     // Prednisone peripheral
  MTX_GUT     // MTX gut
  MTX_C1      // MTX central
  MTX_PG      // MTX polyglutamates (intracellular)
  AZA_GUT     // AZA gut
  SixTGN      // 6-TGN blood (surrogate for 6-TGN in RBC)
  RTX_C1      // Rituximab central
  RTX_C2      // Rituximab peripheral
  IVIG_C      // IVIG (IgG) central
  JAKI_C1     // JAKi central
  JAKI_C2     // JAKi peripheral

// Immunology compartments (8)
  CD8N        // Naïve CD8+ T cells (cells/μL)
  CD8E        // Effector CD8+ T cells in muscle (cells/μL equiv.)
  CD4TH1      // Th1 CD4+ T cells (cells/μL)
  BCELL       // Naïve + activated B cells (cells/μL)
  PLASMA_C    // Long-lived plasma cells (cells/μL)
  AUTOAB      // Autoantibody level (U/mL)
  RTX_CD20    // Free CD20 receptor (nM)
  RTX_BOUND   // RTX-CD20 complex (nM)

// Disease/PD compartments (5)
  IFNG        // IFN-γ (pg/mL)
  TNFa        // TNF-α (pg/mL)
  IL6C        // IL-6 (pg/mL)
  MHCI        // MHC-I expression on muscle (fold-change)
  MUS_INF     // Muscle inflammation index (0–10 scale)
  CK_S        // Serum CK (U/L)
  MMT8        // Manual Muscle Test score (0–80)

$INIT
  PRED_GUT  = 0,    PRED_C1   = 0,    PRED_C2   = 0,
  MTX_GUT   = 0,    MTX_C1    = 0,    MTX_PG    = 0,
  AZA_GUT   = 0,    SixTGN    = 0,
  RTX_C1    = 0,    RTX_C2    = 0,
  IVIG_C    = 8.0,  // Normal endogenous IgG ~8 g/L
  JAKI_C1   = 0,    JAKI_C2   = 0,
  CD8N      = 250,  // Normal naïve CD8 ~250 cells/μL
  CD8E      = 5,    // Small resting effector pool
  CD4TH1    = 150,  // Normal Th1 ~150 cells/μL
  BCELL     = 90,   // Normal B cells ~90/μL
  PLASMA_C  = 20,   // Baseline plasma cells
  AUTOAB    = 5,    // Low baseline autoAb
  RTX_CD20  = 150,  // Free CD20 (nM)
  RTX_BOUND = 0,    // RTX-CD20 complex
  IFNG      = 0.5,  // Basal IFN-γ ~0.5 pg/mL
  TNFa      = 0.3,  // Basal TNF-α ~0.3 pg/mL
  IL6C      = 1.0,  // Basal IL-6 ~1 pg/mL
  MHCI      = 1.0,  // Normal MHC-I expression (fold-change = 1)
  MUS_INF   = 0.5,  // Low baseline inflammation
  CK_S      = 120,  // Normal CK ~120 U/L
  MMT8      = 78    // Near-normal strength at baseline

$MAIN
  // Active prednisolone concentration driving PD
  double PREDNIS = KHB * PRED_C1 / PRED_V1 * 1000; // ng/mL

  // MTX polyglutamate for PD
  double MTX_PG_conc = MTX_PG;  // μg/L equivalent

  // 6-TGN for PD (fmol/8×10^8 RBC equivalent from plasma surrogate)
  double TGN_conc = SixTGN * 10.0; // scaling to physiological units

  // RTX free concentration
  double RTX_free = RTX_C1 / RTX_V1 * 1e3; // nM

  // JAKi concentration (nM)
  double JAKI_Cp = JAKI_C1 / JAKI_V1 * 1e6; // convert mg/L → nM (MW~371 g/mol)

  // Drug PD effects (inhibition; 0=no effect, 1=full inhibition)
  double E_PRED = PRED_IMAX * pow(PREDNIS, PRED_HILL) /
                  (pow(PRED_IC50, PRED_HILL) + pow(PREDNIS, PRED_HILL));
  double E_MTX  = MTX_IMAX * MTX_PG_conc / (MTX_IC50 + MTX_PG_conc);
  double E_TGN  = TGN_IMAX * TGN_conc / (TGN_IC50 + TGN_conc);
  double E_JAKI = JAKI_IMAX * JAKI_Cp / (JAKI_IC50 + JAKI_Cp);

  // RTX-mediated B cell depletion effect
  double RTX_RC = RTX_BOUND; // nM of complex
  double E_RTX_B = RTX_DEPL_IMAX * RTX_RC / (RTX_DEPL_EC50 + RTX_RC);

  // IVIG anti-inflammatory effect (on macrophage/Fc signaling)
  double IVIG_Cp = IVIG_C / IVIG_V;  // g/L
  double E_IVIG  = IVIG_IMAX * IVIG_Cp / (IVIG_IC50 + IVIG_Cp);

  // Combined immunosuppression on CD8 activation/expansion
  double E_COMBO_CD8 = 1.0 - (1.0 - E_PRED) * (1.0 - E_MTX) * (1.0 - E_TGN);
  // Combined on cytokine production
  double E_COMBO_CYT = 1.0 - (1.0 - E_PRED) * (1.0 - E_JAKI) * (1.0 - E_IVIG);
  // Combined on B cell/autoAb
  double E_COMBO_B   = 1.0 - (1.0 - E_RTX_B) * (1.0 - E_PRED) * (1.0 - E_MTX);

  // MHC-I upregulation driven by IFN-γ
  double MHCI_eq = 1.0 + (MHCI_KMAX - 1.0) * IFNG /
                   (MHC_IFNG_EC50 + IFNG);

  // Muscle inflammation drives CK
  double CK_eq = CK_BASE * (1.0 + CK_MAX_FACTOR * MUS_INF / (1.0 + MUS_INF));

$ODE
  // ── Prednisone PK ──────────────────────────────────────────────
  dxdt_PRED_GUT = -PRED_KA * PRED_GUT;
  dxdt_PRED_C1  =  PRED_KA * PRED_GUT
                  - (PRED_CL + PRED_Q) / PRED_V1 * PRED_C1
                  + PRED_Q / PRED_V2 * PRED_C2;
  dxdt_PRED_C2  =  PRED_Q / PRED_V1 * PRED_C1
                  - PRED_Q / PRED_V2 * PRED_C2;

  // ── MTX PK ──────────────────────────────────────────────────────
  dxdt_MTX_GUT = -MTX_KA * MTX_GUT;
  dxdt_MTX_C1  =  MTX_KA * MTX_GUT
                  - (MTX_CL / 28.0 + MTX_KPG) * MTX_C1; // CL/V first order
  dxdt_MTX_PG  =  MTX_KPG * MTX_C1 - MTX_KPGOUT * MTX_PG;

  // ── AZA/6-TGN PK ────────────────────────────────────────────────
  dxdt_AZA_GUT  = -AZA_KA * AZA_GUT;
  double AZA_to_6TGN = AZA_K12 * AZA_GUT * (1.0 - 0.3 * TPMT_act);
  dxdt_SixTGN   =  AZA_to_6TGN - (AZA_6MP_CL / AZA_V2) * SixTGN;

  // ── Rituximab TMDD PK ─────────────────────────────────────────────
  dxdt_RTX_C1   = -(RTX_CL + RTX_Q) / RTX_V1 * RTX_C1
                  + RTX_Q / RTX_V2 * RTX_C2
                  - RTX_KON * RTX_C1 / RTX_V1 * RTX_CD20 * RTX_V1
                  + RTX_KOFF * RTX_BOUND;
  dxdt_RTX_C2   =  RTX_Q / RTX_V1 * RTX_C1 - RTX_Q / RTX_V2 * RTX_C2;
  dxdt_RTX_CD20 = -RTX_KON * (RTX_C1 / RTX_V1) * RTX_CD20
                  + RTX_KOFF * RTX_BOUND
                  + RTX_KINT * (RTX_CD20_0 - RTX_CD20); // CD20 regeneration
  dxdt_RTX_BOUND=  RTX_KON * (RTX_C1 / RTX_V1) * RTX_CD20
                  - (RTX_KOFF + RTX_KINT) * RTX_BOUND;

  // ── IVIG PK (FcRn model) ─────────────────────────────────────────
  double IVIG_FcRn_CL = IVIG_CL0 * (1.0 + IVIG_C / FCRN_CAP); // Nonlinear CL
  dxdt_IVIG_C   = -IVIG_FcRn_CL * IVIG_C;

  // ── JAKi PK ──────────────────────────────────────────────────────
  dxdt_JAKI_C1  = -(JAKI_CL + JAKI_Q) / JAKI_V1 * JAKI_C1
                  + JAKI_Q / JAKI_V2 * JAKI_C2;
  dxdt_JAKI_C2  =  JAKI_Q / JAKI_V1 * JAKI_C1 - JAKI_Q / JAKI_V2 * JAKI_C2;

  // ── CD8+ T cell dynamics ─────────────────────────────────────────
  // Activation amplified by MHC-I and IFN-γ; inhibited by drugs
  double MHCI_stim = MHCI / (1.0 + MHCI);
  dxdt_CD8N  =  CD8_PROD - CD8_KDEATH * CD8N
               - CD8_KACT * CD8N * MHCI_stim * (1.0 - E_COMBO_CD8);
  dxdt_CD8E  =  CD8_KACT * CD8N * MHCI_stim * (1.0 - E_COMBO_CD8) * CD8_KEXP
               - CD8_KMIG * CD8E
               - CD8_KDEATH_E * CD8E;

  // ── CD4+/Th1 dynamics ────────────────────────────────────────────
  dxdt_CD4TH1=  CD4_PROD * (1.0 + IFNG / (2.0 + IFNG))  // IL-12 drives Th1
               - CD4_KDEATH * CD4TH1
               - CD4_KACT * E_COMBO_CD8 * CD4TH1;  // Drugs suppress

  // ── B cell & Plasma cell dynamics ────────────────────────────────
  dxdt_BCELL =  B_PROD - B_KDEATH * BCELL
               - B_KACT * (1.0 - E_COMBO_B) * BCELL
               - E_RTX_B * B_KACT * BCELL;  // RTX depletion
  dxdt_PLASMA_C = PC_PROD_RATE * BCELL * (1.0 - E_COMBO_B)
                  - PC_KDEATH * PLASMA_C
                  - E_RTX_B * 0.3 * PLASMA_C;  // RTX partial effect on PC

  // ── Autoantibody dynamics ─────────────────────────────────────────
  dxdt_AUTOAB = AB_PROD * PLASMA_C * (1.0 - E_COMBO_B)
               - AB_CL * AUTOAB;

  // ── Cytokine dynamics ─────────────────────────────────────────────
  // IFN-γ: produced by Th1 and effector CD8; inhibited by steroids/JAKi
  dxdt_IFNG  =  IFNG_PROD_0
               + TH1_IFNg_PROD * CD4TH1 / 100.0
               + 0.05 * CD8E
               - IFNG_KD * IFNG
               - E_COMBO_CYT * IFNG_PROD_0 * 5.0;

  // TNF-α: from macrophages and Th1; suppressed by drugs
  dxdt_TNFa  =  TNF_PROD_0 * (1.0 + MUS_INF)
               + 0.01 * CD4TH1
               - TNF_KD * TNFa
               - E_COMBO_CYT * TNF_PROD_0 * 4.0;

  // IL-6: from macrophages; drives CRP; suppressed by drugs
  dxdt_IL6C  =  IL6_PROD_0 * (1.0 + 0.5 * MUS_INF)
               + 0.008 * CD4TH1
               - IL6_KD * IL6C
               - E_COMBO_CYT * IL6_PROD_0 * 5.0;

  // ── MHC-I expression on muscle ────────────────────────────────────
  dxdt_MHCI  =  0.5 * (MHCI_eq - MHCI);  // Slow equilibration (t½ ~1.4d)

  // ── Muscle Inflammation Index ─────────────────────────────────────
  // Driven by effector CD8 and autoAb; suppressed by drugs
  double CD8E_stim = MUS_INF_KMAX * CD8E / (MUS_INF_EC50 + CD8E);
  double AB_stim   = 0.02 * AUTOAB / (20.0 + AUTOAB);
  double Drug_suppr = (E_PRED * 0.5 + E_MTX * 0.15 + E_TGN * 0.10 +
                       E_RTX_B * 0.15 + E_IVIG * 0.05 + E_JAKI * 0.20);
  dxdt_MUS_INF = CD8E_stim + AB_stim
                - MUS_INF_KD * MUS_INF * (1.0 + Drug_suppr);

  // ── Serum CK dynamics ─────────────────────────────────────────────
  dxdt_CK_S  =  CK_KPROD * MUS_INF * CK_BASE
               - CK_KD * (CK_S - CK_BASE);  // Mean-reversion + inflammation

  // ── MMT-8 Score dynamics ──────────────────────────────────────────
  dxdt_MMT8  = -MMT_KDECLINE * MUS_INF * MMT8
               + MMT_KREC * (MMT_MAX - MMT8) * Drug_suppr;

$TABLE
  // PK derived quantities
  double Cpred   = PRED_C1 / PRED_V1 * KHB * 1000; // ng/mL prednisolone
  double Cmtx    = MTX_C1 / 28.0 * 1000;            // ng/mL MTX
  double C6TGN   = SixTGN * 10.0;                   // fmol equivalent
  double Cjaki   = JAKI_C1 / JAKI_V1 * 1e6;        // nM

  // PD derived
  double E_PRED_out = PRED_IMAX * pow(Cpred, PRED_HILL) /
                      (pow(PRED_IC50, PRED_HILL) + pow(Cpred, PRED_HILL));
  double E_JAKI_out = JAKI_IMAX * Cjaki / (JAKI_IC50 + Cjaki);

  // Clinical biomarkers
  double CK_fold   = CK_S / CK_BASE;         // CK fold-change over ULN
  double MMT8_pct  = MMT8 / MMT_MAX * 100;   // MMT-8 % of max
  double TIS_approx = (1.0 - MUS_INF / 10.0) * 100; // approx TIS (0-100)
  double IFNsig = IFNG / 0.5;                // IFN signature (normalized)

  // Disease Activity Categories
  double CK_ULN = CK_S / 200.0;              // ULN = 200 U/L
  int Active_Dis = (CK_ULN > 5 || MMT8 < 60) ? 1 : 0;
  int Remission  = (CK_ULN < 1.5 && MMT8 > 72) ? 1 : 0;

$CAPTURE
  Cpred MTX_PG C6TGN Cjaki
  E_PRED_out E_JAKI_out
  CD8E CD4TH1 BCELL PLASMA_C AUTOAB
  IFNG TNFa IL6C MHCI MUS_INF
  CK_S CK_fold MMT8 MMT8_pct TIS_approx
  CK_ULN Active_Dis Remission IFNsig
'

pm_mod <- mcode("polymyositis_qsp", pm_model_code)

## =============================================================================
## Treatment Scenarios
## =============================================================================

## Helper: build event table
make_events <- function(scenario) {
  ev <- switch(scenario,

    # 1. Untreated PM — natural disease course
    "untreated" = ev(amt=0, time=0, cmt=1, evid=2),  # no drug

    # 2. Standard-of-care: High-dose prednisone monotherapy
    # Pred 60 mg/day PO for 4 weeks, then taper
    "pred_mono" = {
      pred_initial <- ev(amt=60, ii=24, addl=27, cmt=1, evid=1) # 60 mg/d × 28d
      pred_taper1  <- ev(amt=40, ii=24, addl=27, cmt=1, evid=1) # 40 mg/d × 28d
      pred_taper2  <- ev(amt=20, ii=24, addl=55, cmt=1, evid=1) # 20 mg/d × 56d
      pred_maint   <- ev(amt=10, ii=24, addl=179,cmt=1, evid=1) # 10 mg/d maint
      pred_initial$time <- 0
      pred_taper1$time  <- 28*24
      pred_taper2$time  <- 56*24
      pred_maint$time   <- 112*24
      bind(pred_initial, pred_taper1, pred_taper2, pred_maint)
    },

    # 3. Pred + MTX combination (standard): Pred + MTX 15 mg/wk SC
    "pred_mtx" = {
      pred_ev <- ev(amt=60, ii=24, addl=27, cmt=1, evid=1)
      pred_t1  <- ev(amt=40, ii=24, addl=27, cmt=1, evid=1)
      pred_t2  <- ev(amt=20, ii=24, addl=111,cmt=1, evid=1)
      pred_t3  <- ev(amt=10, ii=24, addl=179,cmt=1, evid=1)
      pred_ev$time <- 0; pred_t1$time <- 28*24
      pred_t2$time <- 56*24; pred_t3$time <- 112*24
      mtx_ev <- ev(amt=15, ii=168, addl=51, cmt=4, evid=1) # wkly × 52wk
      mtx_ev$time <- 0
      bind(pred_ev, pred_t1, pred_t2, pred_t3, mtx_ev)
    },

    # 4. Pred + AZA combination: AZA 2 mg/kg/day (~150 mg for 75 kg patient)
    "pred_aza" = {
      pred_ev  <- ev(amt=60, ii=24, addl=27, cmt=1, evid=1)
      pred_t1  <- ev(amt=20, ii=24, addl=167,cmt=1, evid=1)
      pred_ev$time <- 0; pred_t1$time <- 28*24
      aza_ev   <- ev(amt=150, ii=24, addl=363, cmt=7, evid=1) # daily × 52wk
      aza_ev$time <- 0
      bind(pred_ev, pred_t1, aza_ev)
    },

    # 5. Rituximab for refractory PM: 1000 mg IV × 2 doses (2-wk apart)
    "rtx_combo" = {
      pred_ev  <- ev(amt=40, ii=24, addl=363, cmt=1, evid=1) # background pred
      pred_ev$time <- 0
      rtx_1    <- ev(amt=1000, cmt=9, evid=1, rate=-2) # infusion over ~4h
      rtx_2    <- ev(amt=1000, cmt=9, evid=1, rate=-2)
      rtx_1$time <- 0; rtx_2$time <- 14*24
      bind(pred_ev, rtx_1, rtx_2)
    },

    # 6. JAK inhibitor (Baricitinib 4 mg/day PO): refractory/ILD-associated PM
    "jaki_combo" = {
      pred_ev  <- ev(amt=20, ii=24, addl=363, cmt=1, evid=1)
      pred_ev$time <- 0
      jaki_ev  <- ev(amt=4, ii=24, addl=363, cmt=12, evid=1)
      jaki_ev$time <- 0
      bind(pred_ev, jaki_ev)
    }
  )
  ev
}

## Simulate all 6 scenarios
scen_names <- c("untreated", "pred_mono", "pred_mtx",
                "pred_aza",  "rtx_combo", "jaki_combo")
scen_labels <- c("Untreated", "Prednisone Mono",
                 "Pred + MTX", "Pred + AZA",
                 "Rituximab + Pred", "JAKi + Pred")

# Disease onset: simulate with elevated initial conditions to represent PM disease state
pm_disease_init <- param(pm_mod,
  PM_SEVERITY = 2.0   # Active PM
)

# Use modified initial conditions for disease state
pm_start <- init(pm_disease_init,
  CD8E    = 80,    # Active effector CD8 in muscle
  CD8N    = 300,
  CD4TH1  = 250,   # Elevated Th1
  BCELL   = 120,   # Elevated B cells
  PLASMA_C= 50,
  AUTOAB  = 45,    # Elevated autoantibodies (anti-Jo-1)
  IFNG    = 18,    # Elevated IFN-γ (PM signature)
  TNFa    = 8,
  IL6C    = 12,
  MHCI    = 6.5,   // MHC-I upregulated on muscle
  MUS_INF = 7.0,   // Active muscle inflammation
  CK_S    = 4800,  // CK ~4800 U/L (40× ULN, typical active PM)
  MMT8    = 38     // Severe weakness (MMT8 38/80)
)

sim_results <- lapply(seq_along(scen_names), function(i) {
  e <- make_events(scen_names[i])
  if (scen_names[i] == "untreated") {
    # Untreated: just run forward with no dosing
    out <- pm_start %>%
      mrgsim(end = 365*24, delta = 12) %>%
      as.data.frame()
  } else {
    out <- pm_start %>%
      ev(e) %>%
      mrgsim(end = 365*24, delta = 12) %>%
      as.data.frame()
  }
  out$Scenario <- scen_labels[i]
  out$time_d   <- out$time / 24
  out
})

sim_all <- bind_rows(sim_results)
sim_all$Scenario <- factor(sim_all$Scenario, levels = scen_labels)

## =============================================================================
## Visualization
## =============================================================================

theme_qsp <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E3F2FD"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

cols6 <- c("#E53935","#1E88E5","#43A047","#FB8C00","#8E24AA","#00ACC1")
names(cols6) <- scen_labels

# 1. Serum CK over time
p_ck <- ggplot(sim_all, aes(time_d, CK_S, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = cols6) +
  scale_y_log10(labels = scales::comma) +
  labs(title = "Serum CK Dynamics", x = "Time (days)", y = "CK (U/L, log scale)",
       color = NULL) +
  annotate("text", x = 10, y = 250, label = "ULN (200 U/L)", size = 3,
           color = "gray50") +
  theme_qsp

# 2. MMT-8 score over time
p_mmt <- ggplot(sim_all, aes(time_d, MMT8, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 72, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = cols6) +
  ylim(0, 80) +
  labs(title = "Manual Muscle Test (MMT-8)", x = "Time (days)", y = "MMT-8 Score (0–80)",
       color = NULL) +
  annotate("text", x = 10, y = 74, label = "Remission threshold (72)", size = 3,
           color = "gray50") +
  theme_qsp

# 3. IFN-γ (disease biomarker)
p_ifng <- ggplot(sim_all, aes(time_d, IFNG, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols6) +
  labs(title = "IFN-γ (Interferon Signature)", x = "Time (days)",
       y = "IFN-γ (pg/mL)", color = NULL) +
  theme_qsp

# 4. Muscle inflammation index
p_inf <- ggplot(sim_all, aes(time_d, MUS_INF, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols6) +
  labs(title = "Muscle Inflammation Index", x = "Time (days)",
       y = "Inflammation Index (0–10)", color = NULL) +
  theme_qsp

# 5. Effector CD8+ T cells in muscle
p_cd8 <- ggplot(sim_all, aes(time_d, CD8E, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols6) +
  labs(title = "Effector CD8+ T Cells (Muscle)", x = "Time (days)",
       y = "CD8E (cells/μL equiv.)", color = NULL) +
  theme_qsp

# 6. B cells (rituximab effect)
p_bcell <- ggplot(sim_all, aes(time_d, BCELL, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = cols6) +
  labs(title = "B Cells (CD19+)", x = "Time (days)",
       y = "B cells (cells/μL)", color = NULL) +
  theme_qsp

# Combine plots
p_combined <- (p_ck + p_mmt) / (p_ifng + p_inf) / (p_cd8 + p_bcell) +
  plot_annotation(
    title = "Polymyositis QSP Model — Treatment Scenario Comparison",
    subtitle = "6-scenario simulation: Untreated vs. Standard & Advanced Therapies",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(p_combined)

## =============================================================================
## PK Profile: Prednisone single dose
## =============================================================================
pk_pred <- pm_mod %>%
  init(IVIG_C = 8) %>%
  ev(amt = 60, cmt = 1, evid = 1) %>%
  mrgsim(end = 72, delta = 0.5) %>%
  as.data.frame()

p_pk_pred <- ggplot(pk_pred, aes(time, Cpred)) +
  geom_line(color = "#1E88E5", linewidth = 1.2) +
  labs(title = "Prednisolone PK after 60 mg Oral Prednisone",
       x = "Time (h)", y = "Prednisolone Concentration (ng/mL)") +
  theme_qsp

print(p_pk_pred)

## =============================================================================
## Parameter Sensitivity — CK at 90 days
## =============================================================================
param_sens <- expand.grid(
  PRED_IC50 = c(1, 2, 4),
  CD8_KACT  = c(0.10, 0.15, 0.25)
)

sens_results <- lapply(seq_len(nrow(param_sens)), function(i) {
  mod_i <- param(pm_start,
    PRED_IC50 = param_sens$PRED_IC50[i],
    CD8_KACT  = param_sens$CD8_KACT[i]
  )
  e <- make_events("pred_mtx")
  out <- mod_i %>% ev(e) %>% mrgsim(end = 90*24, delta = 24) %>% as.data.frame()
  last <- tail(out, 1)
  data.frame(
    PRED_IC50 = param_sens$PRED_IC50[i],
    CD8_KACT  = param_sens$CD8_KACT[i],
    CK_day90  = last$CK_S,
    MMT8_day90 = last$MMT8
  )
})

sens_df <- bind_rows(sens_results)
sens_df$CD8_KACT_label <- paste0("CD8 k_act = ", sens_df$CD8_KACT)

p_sens <- ggplot(sens_df, aes(factor(PRED_IC50), CK_day90, fill = CD8_KACT_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Sensitivity: CK at Day 90 (Pred + MTX)",
       x = "Pred IC50 (ng/mL)", y = "CK at Day 90 (U/L)", fill = NULL) +
  theme_qsp

print(p_sens)

cat("\n=== Polymyositis QSP Simulation Summary ===\n")
cat("Scenarios simulated:", length(scen_labels), "\n")
cat("Compartments: 26 (13 PK + 8 Immunology + 5 Disease/PD)\n")
cat("Simulation horizon: 365 days\n")
for (s in scen_labels) {
  d <- sim_all[sim_all$Scenario == s & sim_all$time_d > 360, ]
  if (nrow(d) > 0) {
    d <- d[1,]
    cat(sprintf("  %-25s | CK: %6.0f U/L | MMT8: %4.1f | Remission: %s\n",
                s, d$CK_S, d$MMT8, ifelse(d$Remission == 1, "YES", "NO")))
  }
}
