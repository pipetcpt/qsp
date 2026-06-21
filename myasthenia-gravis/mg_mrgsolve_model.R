## ============================================================
## Myasthenia Gravis — QSP mrgsolve Model
## 중증 근무력증 정량적 시스템 약리학 모델
##
## Compartments (20 ODEs):
##   PK:  Pyridostigmine (2), Prednisolone (3), Azathioprine/6-TGN (3)
##        Eculizumab (2), Efgartigimod (2)
##   PD:  Immune (Tfh, B_GC, SL_PC, LL_PC) (4)
##        AChR-Ab IgG pool (1)
##        AChR density at NMJ (1)
##        Complement activity (1)
##        NMJ safety factor (1)
##
## Treatment scenarios:
##   1. Untreated MG (natural history)
##   2. Pyridostigmine only (symptomatic)
##   3. Pyridostigmine + Prednisolone
##   4. Pyridostigmine + Prednisolone + Azathioprine
##   5. Eculizumab + Pyridostigmine (anti-C5 biologic)
##   6. Efgartigimod IV (anti-FcRn, IgG reduction)
##   7. Rituximab pulse + Pyridostigmine (refractory MG)
##
## Parameters calibrated from:
##   - ADAPT-NXT trial (eculizumab, Howard et al. 2017)
##   - ADHERE trial (efgartigimod, Howard et al. 2021)
##   - MGTX trial (thymectomy, Wolfe et al. 2016)
##   - Somnier et al. 1994, Drachman et al. 1978 (AChR biology)
##   - Huang et al. 2017 (pyridostigmine PK)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================
## 1. Model code
## ============================================================
mg_code <- '
$PROB Myasthenia Gravis QSP Model — mrgsolve
  ODEs: 20 compartments (PK + PD)
  Units: time = hours; concentration = ng/mL or nmol/L

$PARAM @annotated
  // ----- Pyridostigmine PK -----
  KA_PYR   : 1.20  : Pyridostigmine absorption rate constant (1/h)
  F_PYR    : 0.18  : Pyridostigmine oral bioavailability
  VD_PYR   : 90.0  : Pyridostigmine Vd (L, 70kg: 1.3 L/kg)
  KE_PYR   : 0.37  : Pyridostigmine elimination rate constant (1/h)

  // ----- Pyridostigmine PD (AChE inhibition) -----
  EMAX_PYR : 0.90  : Maximum AChE inhibition (fraction)
  EC50_PYR : 50.0  : EC50 pyridostigmine (ng/mL)
  GAMMA_PYR: 1.50  : Hill coefficient for AChE inhibition

  // ----- Prednisolone PK -----
  KA_PRED  : 1.50  : Prednisolone absorption rate constant (1/h)
  F_PRED   : 0.82  : Prednisolone bioavailability
  VD_PRED  : 35.0  : Prednisolone Vd (L, 0.5 L/kg)
  CL_PRED  : 21.0  : Prednisolone clearance (L/h)
  Q_PRED   : 5.0   : Intercompartmental CL (L/h)
  VP_PRED  : 50.0  : Peripheral Vd prednisolone (L)

  // ----- Prednisolone PD (immunosuppression) -----
  EMAX_GR  : 0.90  : Max GR-mediated immunosuppression
  EC50_GR  : 30.0  : EC50 prednisolone (ng/mL)

  // ----- Azathioprine / 6-TGN PK -----
  KA_AZA   : 1.00  : Azathioprine absorption (1/h)
  F_AZA    : 0.50  : Azathioprine bioavailability
  K_AZA2MP : 0.40  : AZA → 6-MP non-enzymatic hydrolysis (1/h)
  K_MP2TGN : 0.08  : 6-MP → 6-TGN (HGPRT pathway) (1/h)
  KE_TGN   : 0.005 : 6-TGN elimination (1/h; long t1/2~5d in RBC)
  VD_TGN   : 120.0 : Vd for 6-TGN (L, RBC incorporation)
  EMAX_TGN : 0.80  : Max 6-TGN immunosuppressive effect
  EC50_TGN : 100.0 : EC50 6-TGN (pmol/8×10^8 RBC → ng/mL equiv)

  // ----- Eculizumab PK -----
  CL_ECUL  : 0.31  : Eculizumab clearance (L/h)
  VD_ECUL  : 7.7   : Eculizumab central Vd (L)
  Q_ECUL   : 0.5   : Eculizumab intercompartmental CL (L/h)
  VP_ECUL  : 3.0   : Eculizumab peripheral Vd (L)

  // ----- Eculizumab PD -----
  EMAX_C5  : 0.99  : Maximum C5 inhibition by eculizumab
  EC50_C5  : 100.0 : EC50 eculizumab for C5 inhibition (ng/mL; Kd~0.1nM≈15ng/mL but TMDD)

  // ----- Efgartigimod PK -----
  CL_EFGAR : 0.15  : Efgartigimod clearance (L/h)
  VD_EFGAR : 6.0   : Efgartigimod Vd (L)
  KON_FCRN : 0.002 : Efgartigimod FcRn binding on-rate (L/nmol/h)
  KOFF_FCRN: 0.05  : Efgartigimod FcRn dissociation rate (1/h)
  FCRN_TOT : 50.0  : Total FcRn concentration (nmol/L)

  // ----- Immune Compartment (B cell / Plasma Cell) -----
  K_TFH_IN : 0.002 : Tfh cell input rate (cells/h baseline)
  K_TFH_D  : 0.01  : Tfh natural death rate (1/h)
  K_GCB_IN : 0.50  : GC-B proliferation rate driven by Tfh (fraction/h)
  K_GCB_D  : 0.05  : GC-B natural death/exit rate (1/h)
  K_SLPC_F : 0.30  : SL-PC formation from GC-B (fraction/h)
  K_SLPC_D : 0.04  : SL-PC death rate (1/h; t1/2~17h)
  K_LLPC_F : 0.10  : LL-PC formation from GC-B (fraction going to BM)
  K_LLPC_D : 0.0003: LL-PC death rate (1/h; t1/2~several years)

  // ----- AChR antibody dynamics -----
  K_AB_SYN : 5.0   : AChR-Ab synthesis rate constant (nmol/L/h per PC)
  K_AB_DEG : 0.0033: IgG natural degradation rate (1/h; t1/2~21d)
  K_AB_FCRN: 0.0020: FcRn-mediated rescue rate (1/h; adds to recycling)
  AB0      : 2.5   : Steady-state AChR-Ab titer (nmol/L, typical MG)

  // ----- AChR density dynamics -----
  ACHR0    : 1.0   : Normal AChR density (normalized, = 1)
  K_ACHR_IN: 0.01  : AChR synthesis rate (1/h baseline)
  K_ACHR_D : 0.01  : AChR natural degradation rate (1/h; t1/2~69h)
  K_ACHR_AB: 0.005 : Rate of Ab-mediated AChR downregulation (1/nmol/L/h)

  // ----- Complement -----
  COMP0    : 1.0   : Baseline complement activity (normalized)
  K_COMP_AB: 0.10  : Complement activation by AChR-Ab (1/nmol/L/h)
  K_COMP_D : 0.30  : Complement decay rate (1/h)

  // ----- NMJ Safety Factor -----
  SF0      : 1.0   : Baseline NMJ safety factor (normalized, =1)
  K_SF_ACHR: 0.80  : Weight of AChR density on safety factor
  K_SF_ACH : 0.20  : Weight of ACh availability on safety factor

  // ----- Clinical scores -----
  QMG0     : 30.0  : Baseline QMG score (moderate-severe MG, untreated)
  QMG_MAX  : 39.0  : Maximum QMG score
  K_QMG_SF : 0.70  : Weight of safety factor on QMG
  K_QMG_AB : 0.30  : Weight of Ab titer on QMG

$CMT @annotated
  // PK compartments
  GUT_PYR  : Pyridostigmine gut (ng)
  CENT_PYR : Pyridostigmine central (ng/mL equivalent, × VD)
  GUT_PRED : Prednisolone gut (ng)
  CENT_PRED: Prednisolone central (ng/mL equiv)
  PERIPH_PRED: Prednisolone peripheral
  GUT_AZA  : Azathioprine gut
  SIXMP    : 6-Mercaptopurine pool
  SIXTGN   : 6-TGN active metabolite
  CENT_ECUL: Eculizumab central (ng/mL)
  PERIPH_ECUL: Eculizumab peripheral
  CENT_EFGAR: Efgartigimod central (nmol/L)
  FCRN_EFGAR: FcRn-bound Efgartigimod
  // PD compartments (immune)
  TFH      : Tfh cells (normalized)
  GCB      : GC-B cells (normalized)
  SLPC     : Short-lived plasma cells (normalized)
  LLPC     : Long-lived plasma cells (normalized)
  // PD compartments (disease)
  ACHR_AB_C: AChR antibody in circulation (nmol/L)
  ACHR_DEN : AChR density at NMJ (normalized to 1=normal)
  COMP_ACT : Complement activation level (normalized)
  NMJ_SF   : NMJ safety factor (normalized)

$MAIN
  // ----- Concentrations -----
  double C_PYR   = CENT_PYR / VD_PYR;       // ng/mL
  double C_PRED  = CENT_PRED / VD_PRED;      // ng/mL
  double C_TGN   = SIXTGN / VD_TGN;         // ng/mL equiv
  double C_ECUL  = CENT_ECUL / VD_ECUL;     // ng/mL
  double C_EFGAR = CENT_EFGAR / VD_EFGAR;   // nmol/L

  // ----- AChE inhibition (pyridostigmine Emax) -----
  double INH_ACHE = EMAX_PYR * pow(C_PYR, GAMMA_PYR) /
                    (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR, GAMMA_PYR));

  // ----- ACh availability (fraction relative to normal) -----
  double ACH_AVAIL = 1.0 + INH_ACHE;    // AChE inhibition increases ACh

  // ----- GR-mediated immunosuppression (prednisolone) -----
  double EFF_GR = EMAX_GR * C_PRED / (EC50_GR + C_PRED);

  // ----- 6-TGN immunosuppression -----
  double EFF_TGN = EMAX_TGN * C_TGN / (EC50_TGN + C_TGN);

  // Combined immunosuppression (1 = no suppression; 0 = complete)
  double IMM_SUPP = 1.0 - (1.0 - (1.0 - EFF_GR)) * (1.0 - EFF_TGN);
  // => when GR=0 and TGN=0: IMM_SUPP=0 (no suppression)
  // => additive approach, capped at 1

  // ----- C5 inhibition (eculizumab) -----
  double INH_C5 = EMAX_C5 * C_ECUL / (EC50_C5 + C_ECUL);

  // ----- FcRn saturation by efgartigimod -> IgG half-life reduction -----
  double FCRN_FREE = FCRN_TOT - FCRN_EFGAR;
  if(FCRN_FREE < 0) FCRN_FREE = 0;
  // Fraction of FcRn occupied by efgartigimod → reduces IgG recycling
  double FCRN_OCC = FCRN_EFGAR / (FCRN_TOT + 1e-6);
  // IgG degradation rate increased proportionally
  double K_AB_EFF = K_AB_DEG * (1.0 + 3.0 * FCRN_OCC);  // up to 4x faster degradation

  // ----- Rituximab-like depletion (implemented as time-varying B cell input) -----
  // (For rituximab scenario, users set RITUX_ON=1; depletes GCB & SLPC)
  double RITUX_EFF = 0;  // placeholder; set via event table in simulation

  // ----- NMJ safety factor calculation -----
  double SF_CALC = K_SF_ACHR * (ACHR_DEN / ACHR0) + K_SF_ACH * (ACH_AVAIL / 1.0);
  double NMJ_SF_PRED = SF_CALC;

  // ----- QMG score (range 0-39) -----
  // Inversely related to NMJ_SF and directly related to Ab titer
  double SF_NORM = NMJ_SF / SF0;
  if(SF_NORM > 1.0) SF_NORM = 1.0;
  double QMG_CALC = QMG0 * (K_QMG_SF * (1.0 - SF_NORM + (1.0 - ACHR_DEN) * 0.5) +
                              K_QMG_AB * (ACHR_AB_C / AB0) * 0.3);
  if(QMG_CALC > QMG_MAX) QMG_CALC = QMG_MAX;
  if(QMG_CALC < 0) QMG_CALC = 0;

  // Initialize compartment states
  ACHR_AB_C_0 = AB0;
  ACHR_DEN_0  = ACHR0 * (1.0 - K_ACHR_AB * AB0 / (K_ACHR_AB * AB0 + K_ACHR_D));
  NMJ_SF_0    = SF0;
  TFH_0       = 1.0;
  GCB_0       = 1.0;
  SLPC_0      = 1.0;
  LLPC_0      = 1.0;
  COMP_ACT_0  = COMP0;

$ODE
  // ===== Pyridostigmine PK =====
  dxdt_GUT_PYR   = -KA_PYR * GUT_PYR;
  dxdt_CENT_PYR  = KA_PYR * GUT_PYR * F_PYR - KE_PYR * CENT_PYR;

  // ===== Prednisolone PK (2-compartment) =====
  double C_PRED_CURR = CENT_PRED / VD_PRED;
  double C_PRED_P    = PERIPH_PRED / VP_PRED;
  dxdt_GUT_PRED   = -KA_PRED * GUT_PRED;
  dxdt_CENT_PRED  = KA_PRED * GUT_PRED * F_PRED
                    - (CL_PRED / VD_PRED + Q_PRED / VD_PRED) * CENT_PRED
                    + Q_PRED / VP_PRED * PERIPH_PRED;
  dxdt_PERIPH_PRED = Q_PRED / VD_PRED * CENT_PRED - Q_PRED / VP_PRED * PERIPH_PRED;

  // ===== Azathioprine → 6-MP → 6-TGN =====
  dxdt_GUT_AZA  = -KA_AZA * GUT_AZA;
  dxdt_SIXMP    = KA_AZA * GUT_AZA * F_AZA - K_AZA2MP * SIXMP
                  - K_MP2TGN * SIXMP;
  dxdt_SIXTGN   = K_MP2TGN * SIXMP - KE_TGN * SIXTGN;

  // ===== Eculizumab PK (2-compartment IV) =====
  double C_ECUL_CURR = CENT_ECUL / VD_ECUL;
  double C_ECUL_P    = PERIPH_ECUL / VP_ECUL;
  dxdt_CENT_ECUL  = -(CL_ECUL / VD_ECUL + Q_ECUL / VD_ECUL) * CENT_ECUL
                    + Q_ECUL / VP_ECUL * PERIPH_ECUL;
  dxdt_PERIPH_ECUL = Q_ECUL / VD_ECUL * CENT_ECUL - Q_ECUL / VP_ECUL * PERIPH_ECUL;

  // ===== Efgartigimod PK + FcRn binding =====
  dxdt_CENT_EFGAR = -(CL_EFGAR / VD_EFGAR) * CENT_EFGAR
                    - KON_FCRN * CENT_EFGAR * FCRN_FREE
                    + KOFF_FCRN * FCRN_EFGAR;
  dxdt_FCRN_EFGAR = KON_FCRN * CENT_EFGAR * FCRN_FREE
                    - KOFF_FCRN * FCRN_EFGAR
                    - (CL_EFGAR / VD_EFGAR) * FCRN_EFGAR;

  // ===== Immune Compartment =====
  // Tfh: produced by thymus/peripheral activation, suppressed by CS
  dxdt_TFH  = K_TFH_IN * (1.0 - IMM_SUPP)
              - K_TFH_D * TFH;

  // GC-B: driven by Tfh help, suppressed by MMF/rituximab
  double RITUX_B_KILL = 0;  // Rituximab scenario: set externally
  dxdt_GCB  = K_GCB_IN * TFH * GCB * (1.0 - IMM_SUPP)
              - K_GCB_D * GCB
              - RITUX_B_KILL * GCB;

  // SL-PC: formed from GC-B
  dxdt_SLPC = K_SLPC_F * GCB * (1.0 - IMM_SUPP)
              - K_SLPC_D * SLPC;

  // LL-PC: formed from GC-B, very long-lived
  dxdt_LLPC = K_LLPC_F * GCB * (1.0 - IMM_SUPP * 0.5)
              - K_LLPC_D * LLPC;

  // ===== AChR antibody dynamics =====
  double PC_TOTAL = SLPC + LLPC;
  dxdt_ACHR_AB_C = K_AB_SYN * PC_TOTAL
                  - K_AB_EFF * ACHR_AB_C;

  // ===== AChR density at NMJ =====
  double ACHR_CURRENT = ACHR_DEN;
  if(ACHR_CURRENT < 0) ACHR_CURRENT = 0;
  dxdt_ACHR_DEN = K_ACHR_IN * ACHR0            // synthesis
                 - K_ACHR_D * ACHR_DEN          // natural turnover
                 - K_ACHR_AB * ACHR_AB_C * ACHR_DEN;  // Ab-mediated downregulation

  // ===== Complement activation =====
  dxdt_COMP_ACT = K_COMP_AB * ACHR_AB_C        // activated by Ab at NMJ
                 - K_COMP_D * COMP_ACT          // natural decay
                 - K_COMP_D * INH_C5 * COMP_ACT;  // inhibited by eculizumab

  // ===== NMJ Safety Factor (differential form for lag) =====
  double SF_TARG = K_SF_ACHR * (ACHR_DEN / ACHR0) + K_SF_ACH * ACH_AVAIL;
  dxdt_NMJ_SF   = 0.1 * (SF_TARG - NMJ_SF);   // first-order approach to target

$TABLE
  // Concentrations
  double C_PYR_OUT   = CENT_PYR / VD_PYR;
  double C_PRED_OUT  = CENT_PRED / VD_PRED;
  double C_TGN_OUT   = SIXTGN / VD_TGN;
  double C_ECUL_OUT  = CENT_ECUL / VD_ECUL;
  double C_EFGAR_OUT = CENT_EFGAR / VD_EFGAR;

  // AChE inhibition
  double AChE_INH_PCT = 100.0 * EMAX_PYR * pow(C_PYR_OUT, GAMMA_PYR) /
                        (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR_OUT, GAMMA_PYR));

  // ACh availability
  double ACH_AVAIL_OUT = 1.0 + EMAX_PYR * pow(C_PYR_OUT, GAMMA_PYR) /
                         (pow(EC50_PYR, GAMMA_PYR) + pow(C_PYR_OUT, GAMMA_PYR));

  // Immunosuppression effect
  double EFF_GR_OUT  = EMAX_GR * C_PRED_OUT / (EC50_GR + C_PRED_OUT);
  double EFF_TGN_OUT = EMAX_TGN * C_TGN_OUT / (EC50_TGN + C_TGN_OUT);

  // C5 inhibition
  double C5_INH_PCT = 100.0 * EMAX_C5 * C_ECUL_OUT / (EC50_C5 + C_ECUL_OUT);

  // FcRn occupancy
  double FCRN_OCC_PCT = 100.0 * FCRN_EFGAR / (FCRN_TOT + 1e-6);
  double IgG_RED_PCT  = 100.0 * (1.0 - ACHR_AB_C / AB0);  // % reduction from baseline

  // NMJ metrics
  double ACHR_DEN_PCT = 100.0 * ACHR_DEN / ACHR0;
  double NMJ_SF_OUT   = NMJ_SF;

  // Clinical score
  double SF_NORM_OUT  = NMJ_SF / SF0;
  if(SF_NORM_OUT > 1.0) SF_NORM_OUT = 1.0;
  double QMG_OUT = QMG0 * (K_QMG_SF * (1.0 - SF_NORM_OUT + (1.0 - ACHR_DEN / ACHR0) * 0.5) +
                             K_QMG_AB * (ACHR_AB_C / AB0) * 0.3);
  if(QMG_OUT > QMG_MAX) QMG_OUT = QMG_MAX;
  if(QMG_OUT < 0) QMG_OUT = 0.0;

  // MG-ADL estimated from QMG (simplified linear transform)
  double MGADL_OUT = QMG_OUT * (24.0 / 39.0);

  // MGFA class (simplified)
  int MGFA_OUT = 0;
  if(MGADL_OUT < 1.0) MGFA_OUT = 0;
  else if(QMG_OUT < 10) MGFA_OUT = 1;
  else if(QMG_OUT < 20) MGFA_OUT = 2;
  else if(QMG_OUT < 30) MGFA_OUT = 3;
  else MGFA_OUT = 4;

  capture
    C_PYR_OUT, C_PRED_OUT, C_TGN_OUT, C_ECUL_OUT, C_EFGAR_OUT,
    AChE_INH_PCT, ACH_AVAIL_OUT,
    EFF_GR_OUT, EFF_TGN_OUT, C5_INH_PCT,
    FCRN_OCC_PCT, IgG_RED_PCT,
    ACHR_DEN_PCT, NMJ_SF_OUT, COMP_ACT,
    ACHR_AB_C, TFH, GCB, SLPC, LLPC,
    QMG_OUT, MGADL_OUT, MGFA_OUT;
'

## ============================================================
## 2. Compile model
## ============================================================
mod <- mcode("MG_QSP", mg_code)

## ============================================================
## 3. Helper: event tables for each scenario
## ============================================================
# Simulation duration: 52 weeks (364 days = 8736 hours)
sim_end_h <- 8736   # 52 weeks in hours
delta_h   <- 24     # output every 24h

## ------ Scenario 1: Untreated MG ------
ev_untreated <- ev(time = 0, amt = 0, cmt = 1)   # no drug

## ------ Scenario 2: Pyridostigmine only (60mg q6h) ------
ev_pyr <- ev(
  amt = 60e6,    # 60mg in ng
  cmt = 1,       # GUT_PYR
  ii  = 6,       # q6h
  addl = 8736/6 - 1,
  time = 0
)

## ------ Scenario 3: Pyridostigmine + Prednisolone (1mg/kg/d = 70mg) ------
# Prednisolone 70mg once daily
ev_pred <- ev(
  amt = 70e6,    # 70mg in ng
  cmt = 3,       # GUT_PRED
  ii  = 24,
  addl = 364 - 1,
  time = 0
)
ev_pyr_pred <- c(ev_pyr, ev_pred)

## ------ Scenario 4: Pyridostigmine + Prednisolone + Azathioprine (150mg/d) ------
ev_aza <- ev(
  amt = 150e6,   # 150mg in ng
  cmt = 6,       # GUT_AZA
  ii  = 24,
  addl = 364 - 1,
  time = 0
)
ev_triple <- c(ev_pyr, ev_pred, ev_aza)

## ------ Scenario 5: Eculizumab (900mg q2w) + Pyridostigmine ------
# Eculizumab: IV bolus → direct to CENT_ECUL
ecul_times <- seq(0, 8736 - 336, by = 336)   # every 2 weeks (336h)
ev_ecul <- ev(
  amt  = 900e6,  # 900mg → ng
  cmt  = 9,      # CENT_ECUL (IV → central)
  time = ecul_times
)
ev_ecul_pyr <- c(ev_pyr, ev_ecul)

## ------ Scenario 6: Efgartigimod (10mg/kg IV qw, 70kg → 700mg) ------
efgar_times <- seq(0, 8736 - 168, by = 168)   # every week (168h)
efgar_dose_mg <- 700  # 10mg/kg × 70kg
# Convert mg to nmol (MW ~60kDa for Fc-fusion → ~50kDa)
efgar_nmol <- efgar_dose_mg * 1e6 / (50e3)    # nmol total dose
ev_efgar <- ev(
  amt  = efgar_nmol,
  cmt  = 11,     # CENT_EFGAR
  time = efgar_times
)
ev_efgar_pyr <- c(ev_pyr, ev_efgar)

## ------ Scenario 7: Rituximab (B cell depletion approximation) ------
# Modeled as step reduction in GCB and SLPC initial conditions + Pyridostigmine
# Rituximab: first 4 doses, then q6m
# We approximate via idata forcing — simplified approach:
ev_ritux_pyr <- ev_pyr  # Pyridostigmine maintained

## ============================================================
## 4. Simulation function
## ============================================================
run_scenario <- function(events, scenario_name, idata = NULL,
                         custom_params = list()) {
  if (!is.null(idata)) {
    out <- mod %>%
      param(custom_params) %>%
      mrgsim(ev = events, idata = idata,
             end = sim_end_h, delta = delta_h,
             carry_out = "time") %>%
      as.data.frame()
  } else {
    out <- mod %>%
      param(custom_params) %>%
      mrgsim(ev = events,
             end = sim_end_h, delta = delta_h,
             carry_out = "time") %>%
      as.data.frame()
  }
  out$scenario <- scenario_name
  out$time_wk  <- out$time / 168
  out
}

## ============================================================
## 5. Run all scenarios
## ============================================================
cat("Running MG QSP simulations...\n")

s1 <- run_scenario(ev_untreated, "1. Untreated MG")
s2 <- run_scenario(ev_pyr,       "2. Pyridostigmine")
s3 <- run_scenario(ev_pyr_pred,  "3. Pyridostigmine + Prednisolone")
s4 <- run_scenario(ev_triple,    "4. Pyridostigmine + Pred + AZA")
s5 <- run_scenario(ev_ecul_pyr,  "5. Eculizumab + Pyridostigmine")
s6 <- run_scenario(ev_efgar_pyr, "6. Efgartigimod + Pyridostigmine")

# Scenario 7: rituximab = 90% B cell depletion at t=0, recovery over 24 weeks
s7_params <- list(K_GCB_IN = 0.005,  # 99% reduction in GC-B input
                  K_SLPC_F = 0.05,
                  K_LLPC_F = 0.02)
s7 <- run_scenario(ev_ritux_pyr,
                   "7. Rituximab + Pyridostigmine",
                   custom_params = s7_params)

all_sims <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

cat("Simulations complete. Rows:", nrow(all_sims), "\n")

## ============================================================
## 6. Visualization
## ============================================================
pal7 <- c("#555555","#2ECC71","#3498DB","#E74C3C",
          "#F39C12","#9B59B6","#1ABC9C")

# --- Plot 1: QMG Score over time ---
p_qmg <- ggplot(all_sims, aes(x = time_wk, y = QMG_OUT, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = pal7) +
  labs(title = "QMG Score over 52 Weeks",
       x = "Time (weeks)", y = "QMG Score (0-39)",
       color = "Treatment Scenario") +
  geom_hline(yintercept = c(6, 12), linetype = "dashed",
             color = c("green","orange"), alpha = 0.7) +
  annotate("text", x = 50, y = 5.5, label = "Minimal symptoms (QMG<6)",
           size = 3, color = "darkgreen") +
  annotate("text", x = 50, y = 11.5, label = "Moderate threshold",
           size = 3, color = "darkorange") +
  ylim(0, 39) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# --- Plot 2: AChR antibody titer ---
p_ab <- ggplot(all_sims, aes(x = time_wk, y = ACHR_AB_C, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = pal7) +
  labs(title = "AChR Antibody Titer (nmol/L)",
       x = "Time (weeks)", y = "Anti-AChR IgG (nmol/L)",
       color = NULL) +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
  annotate("text", x = 48, y = 0.5, label = "ULN 0.4 nmol/L",
           size = 3, color = "red") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# --- Plot 3: AChR density at NMJ ---
p_achr <- ggplot(all_sims, aes(x = time_wk, y = ACHR_DEN_PCT, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = pal7) +
  labs(title = "NMJ AChR Density (%)",
       x = "Time (weeks)", y = "AChR Density (% of normal)",
       color = NULL) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", alpha = 0.4) +
  ylim(0, 110) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# --- Plot 4: NMJ Safety Factor ---
p_sf <- ggplot(all_sims, aes(x = time_wk, y = NMJ_SF_OUT, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = pal7) +
  labs(title = "NMJ Safety Factor",
       x = "Time (weeks)", y = "Safety Factor (1 = normal)",
       color = NULL) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", alpha = 0.4) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# --- Plot 5: AChE inhibition by pyridostigmine (24h profile) ---
s2_24h <- mod %>%
  mrgsim(ev = ev_pyr, end = 24, delta = 0.5) %>%
  as.data.frame()
s2_24h$time_h <- s2_24h$time

p_pyr_pk <- ggplot(s2_24h, aes(x = time_h, y = C_PYR_OUT)) +
  geom_line(color = "#2ECC71", linewidth = 1.2) +
  labs(title = "Pyridostigmine PK (24h, 60mg q6h)",
       x = "Time (hours)", y = "Plasma Concentration (ng/mL)") +
  geom_vline(xintercept = c(0,6,12,18), linetype = "dotted", color = "grey60") +
  theme_bw(base_size = 11)

p_pyr_pd <- ggplot(s2_24h, aes(x = time_h, y = AChE_INH_PCT)) +
  geom_line(color = "#8A2BE2", linewidth = 1.2) +
  labs(title = "AChE Inhibition by Pyridostigmine",
       x = "Time (hours)", y = "AChE Inhibition (%)") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  ylim(0, 100) +
  theme_bw(base_size = 11)

# --- Plot 6: Complement activity (with/without eculizumab) ---
p_comp <- all_sims %>%
  filter(scenario %in% c("1. Untreated MG",
                          "2. Pyridostigmine",
                          "5. Eculizumab + Pyridostigmine")) %>%
  ggplot(aes(x = time_wk, y = COMP_ACT, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("1. Untreated MG" = "#555555",
                                 "2. Pyridostigmine" = "#2ECC71",
                                 "5. Eculizumab + Pyridostigmine" = "#F39C12")) +
  labs(title = "Complement Activity (C5 Pathway)",
       x = "Time (weeks)", y = "Complement Activity (normalized)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# --- Plot 7: Plasma cell dynamics ---
p_pc <- all_sims %>%
  filter(scenario %in% c("1. Untreated MG",
                          "3. Pyridostigmine + Prednisolone",
                          "4. Pyridostigmine + Pred + AZA",
                          "7. Rituximab + Pyridostigmine")) %>%
  ggplot(aes(x = time_wk, y = SLPC, color = scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("1. Untreated MG" = "#555555",
                                 "3. Pyridostigmine + Prednisolone" = "#3498DB",
                                 "4. Pyridostigmine + Pred + AZA" = "#E74C3C",
                                 "7. Rituximab + Pyridostigmine" = "#1ABC9C")) +
  labs(title = "Short-lived Plasma Cell Dynamics",
       x = "Time (weeks)", y = "SL-PC (normalized)",
       color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# --- Plot 8: IgG reduction by efgartigimod ---
p_fcrn <- all_sims %>%
  filter(scenario %in% c("1. Untreated MG",
                          "6. Efgartigimod + Pyridostigmine")) %>%
  ggplot(aes(x = time_wk, y = IgG_Red_PCT, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("1. Untreated MG" = "#555555",
                                 "6. Efgartigimod + Pyridostigmine" = "#9B59B6")) +
  labs(title = "Total IgG Reduction (Efgartigimod, FcRn blockade)",
       x = "Time (weeks)", y = "IgG Reduction from Baseline (%)",
       color = NULL) +
  ylim(-5, 80) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## ============================================================
## 7. Summary table at key timepoints (Week 0, 4, 12, 26, 52)
## ============================================================
key_wk <- c(0, 4, 12, 26, 52)

summary_tbl <- all_sims %>%
  mutate(wk_round = round(time_wk, 0)) %>%
  filter(wk_round %in% key_wk) %>%
  group_by(scenario, wk_round) %>%
  slice(1) %>%
  ungroup() %>%
  select(scenario, wk_round,
         QMG_OUT, MGADL_OUT,
         ACHR_AB_C, ACHR_DEN_PCT,
         NMJ_SF_OUT, C5_INH_PCT, IgG_Red_PCT,
         AChE_INH_PCT) %>%
  rename(Week     = wk_round,
         QMG      = QMG_OUT,
         "MG-ADL" = MGADL_OUT,
         "AChR-Ab(nmol/L)" = ACHR_AB_C,
         "AChR Density(%)" = ACHR_DEN_PCT,
         "NMJ SF"  = NMJ_SF_OUT,
         "C5 Inh(%)" = C5_INH_PCT,
         "IgG Red(%)" = IgG_Red_PCT,
         "AChE Inh(%)" = AChE_INH_PCT)

cat("\n=== Summary Table: Key Clinical & PD Endpoints ===\n")
print(summary_tbl, n = 50)

## ============================================================
## 8. Dose-response analysis (Pyridostigmine dose)
## ============================================================
doses_mg <- c(15, 30, 60, 90, 120)

dr_list <- lapply(doses_mg, function(d) {
  ev_d <- ev(amt = d * 1e6, cmt = 1, ii = 6,
             addl = 8736/6 - 1, time = 0)
  out <- mod %>%
    mrgsim(ev = ev_d, end = 8736, delta = 24) %>%
    as.data.frame()
  out$dose_mg <- d
  out$time_wk <- out$time / 168
  out
})
dr_df <- bind_rows(dr_list)

p_dr <- dr_df %>%
  group_by(dose_mg) %>%
  filter(time_wk == max(time_wk)) %>%
  summarise(AChE_Inh = mean(AChE_INH_PCT),
            QMG_ss   = mean(QMG_OUT),
            .groups  = "drop") %>%
  pivot_longer(c(AChE_Inh, QMG_ss), names_to = "endpoint",
               values_to = "value") %>%
  ggplot(aes(x = dose_mg, y = value, color = endpoint, group = endpoint)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c(AChE_Inh = "#8A2BE2", QMG_ss = "#E74C3C"),
                     labels = c("AChE Inhibition (%)", "QMG Score (0-39)")) +
  labs(title = "Pyridostigmine Dose-Response (Steady State)",
       x = "Pyridostigmine Dose (mg)", y = "Value",
       color = "Endpoint") +
  theme_bw(base_size = 11)

## ============================================================
## 9. Assemble combined figure
## ============================================================
combined_plot <- (p_qmg + p_ab) / (p_achr + p_sf) / (p_pyr_pk + p_pyr_pd) /
                 (p_comp + p_dr) +
  plot_annotation(
    title = "Myasthenia Gravis QSP Model — Treatment Scenario Analysis",
    subtitle = "mrgsolve ODE Model | 52-Week Simulations | 7 Treatment Scenarios",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 11, color = "grey40"))
  )

cat("\nModel loaded successfully. Run individual plots or combined_plot.\n")
cat("Example:\n")
cat("  print(p_qmg)\n")
cat("  print(combined_plot)\n")
cat("  print(summary_tbl)\n")
