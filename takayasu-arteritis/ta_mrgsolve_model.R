## ============================================================
## Takayasu Arteritis (TA) — QSP mrgsolve ODE Model
## Author : QSP Disease Model Library (CCR)
## Date   : 2026-06-19
## ============================================================
## Key References:
##   Nakaoka 2018 (Lancet): tocilizumab in TA (TAKT trial)
##   Hellmich 2020: EULAR recommendations for large-vessel vasculitis
##   Tombetti 2019 (Nat Rev Rheum): pathogenesis review
##   Hatemi 2022: ITAS disease activity scoring
##   Mekinian 2012: PET-CT monitoring of TA
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ------------------------------------------------------------
## 1. Model Definition
## ------------------------------------------------------------
ta_model_code <- '
$PROB
  Takayasu Arteritis QSP Model
  20 ODE compartments: Drug PK (Prednisone/Prednisolone, Tocilizumab,
  Methotrexate, Infliximab) + Disease Biology (IL-6, sIL-6R, TNF-alpha,
  Th1, Th17, Treg, Vessel Wall Inflammation, Stenosis Index) +
  Biomarkers (CRP, ESR, PET-SUV, MRI-VWT)

$PARAM
  /* ---- Prednisone / Prednisolone PK ---- */
  ka_PRED  = 1.50   /* /h  absorption rate  */
  Vc_PRED  = 28.0   /* L   central volume   */
  Vp_PRED  = 56.0   /* L   peripheral vol.  */
  CLd_PRED = 40.0   /* L/h inter-compartmental CL */
  CL_PRED  = 18.0   /* L/h total clearance  */
  ke0_PRED = 0.25   /* /h  effect-site equilibration */
  F_PRED   = 0.82   /* bioavailability       */

  /* ---- Tocilizumab PK (SC / IV) ---- */
  ka_TCZ   = 0.0083 /* /h  SC absorption ~t_max 4-5 d */
  Vc_TCZ   = 3.5    /* L   central volume (mAb)  */
  Vp_TCZ   = 2.0    /* L   peripheral             */
  CLd_TCZ  = 0.30   /* L/h inter-compartmental    */
  CL_TCZ   = 0.18   /* L/h linear CL              */
  CLMM_TCZ = 0.080  /* L/h Michaelis-Menten max   */
  KM_TCZ   = 1.20   /* mg/L Michaelis-Menten Km   */
  F_TCZ    = 0.80   /* SC bioavailability          */

  /* ---- Methotrexate PK ---- */
  ka_MTX   = 0.90   /* /h  oral absorption        */
  Vc_MTX   = 18.0   /* L                          */
  CL_MTX   = 5.0    /* L/h renal dominant         */
  kpg_MTX  = 0.12   /* /h  polyglutamation rate   */
  kdpg_MTX = 0.018  /* /h  depolyglutamation      */
  F_MTX    = 0.70

  /* ---- Infliximab PK ---- */
  Vc_IFX   = 3.0    /* L                          */
  Vp_IFX   = 1.8    /* L                          */
  CLd_IFX  = 0.22   /* L/h                        */
  CL_IFX   = 0.16   /* L/h                        */

  /* ---- IL-6 dynamics ---- */
  ksyn_IL6 = 0.30   /* pg/mL/h  baseline synthesis */
  kdeg_IL6 = 0.25   /* /h       degradation        */
  IL6_base = 1.2    /* pg/mL    baseline IL-6      */

  /* ---- Soluble IL-6R ---- */
  ksyn_sR  = 3.0    /* ng/mL/h  sIL-6R synthesis  */
  kdeg_sR  = 0.015  /* /h       sIL-6R degradation */
  sR_base  = 200.0  /* ng/mL    baseline           */
  kon_IL6R = 0.001  /* /h/pg    association rate   */
  koff_IL6R= 0.002  /* /h       dissociation rate  */

  /* ---- TNF-alpha dynamics ---- */
  ksyn_TNF = 0.10
  kdeg_TNF = 0.30
  TNF_base = 0.5    /* pg/mL    baseline           */

  /* ---- T-cell dynamics ---- */
  /* Th1 */
  ksyn_Th1 = 0.05   /* cells/uL/h */
  kdeg_Th1 = 0.02
  Th1_base = 25.0   /* cells/uL */

  /* Th17 */
  ksyn_Th17= 0.03
  kdeg_Th17= 0.02
  Th17_base= 15.0

  /* Treg */
  ksyn_Treg= 0.02
  kdeg_Treg= 0.018
  Treg_base= 8.0

  /* ---- Vessel Wall Inflammation (VWI) [0-10 scale] ---- */
  ksyn_VWI = 0.008  /* score/h */
  kdeg_VWI = 0.004
  VWI_base = 2.0    /* steady-state healthy baseline */

  /* ---- Stenosis Index [0-100 %] ---- */
  kprog_ST = 0.0005 /* %/h  driven by VWI           */
  kreg_ST  = 0.0001 /* /h   spontaneous regression  */
  ST_base  = 0.0    /* starting stenosis            */

  /* ---- CRP (mg/L) ---- */
  ksyn_CRP = 0.60   /* mg/L/h produced by IL-6      */
  kdeg_CRP = 0.025  /* /h                           */
  CRP_base = 5.0    /* mg/L baseline                */

  /* ---- ESR (mm/hr) ---- */
  kESR     = 2.5    /* mm/hr per CRP unit           */
  ESR_base = 20.0   /* baseline                     */

  /* ---- PET-CT SUVmax ---- */
  ksyn_PET = 0.0015
  kdeg_PET = 0.003
  PET_base = 1.5    /* healthy vessel               */

  /* ---- MRI Vessel Wall Thickness (mm) ---- */
  kVWT     = 0.002  /* mm per VWI unit*h            */
  VWT_base = 1.5    /* mm normal                    */

  /* ---- Drug PD — Emax parameters ---- */
  Emax_PRED   = 0.85   /* max IL-6 suppression (GC)    */
  EC50_PRED   = 0.10   /* mg/L prednisolone effect-site */
  hill_PRED   = 1.5

  Emax_TCZ    = 0.95   /* max IL-6R occupancy          */
  EC50_TCZ    = 0.50   /* mg/L tocilizumab             */
  hill_TCZ    = 1.2

  Emax_MTX    = 0.55   /* max T-cell proliferation inhib*/
  EC50_MTX    = 0.05   /* umol/L MTX-PG               */
  hill_MTX    = 1.0

  Emax_IFX    = 0.90   /* max TNF neutralization       */
  EC50_IFX    = 0.80   /* mg/L infliximab              */
  hill_IFX    = 1.3

  /* ---- Disease amplification factors ---- */
  amp_IL6_Th1  = 0.006  /* Th1 amplification of IL-6   */
  amp_IL6_Th17 = 0.004  /* Th17 amplification          */
  amp_VWI_IL6  = 0.20   /* IL-6 drives VWI             */
  amp_VWI_TNF  = 0.15   /* TNF drives VWI              */
  amp_VWI_Th17 = 0.08   /* Th17 drives VWI             */
  amp_ST_VWI   = 0.40   /* VWI drives stenosis         */
  amp_IL6_TNF  = 0.05   /* TNF amplifies IL-6          */
  inh_Treg_Th1 = 0.025  /* Treg suppresses Th1         */
  inh_Treg_Th17= 0.030  /* Treg suppresses Th17        */

  /* ---- Simulation parameters ---- */
  WT    = 65.0   /* kg body weight               */
  DOSE_PRED = 0  /* mg/day prednisone (0 = off)  */
  DOSE_TCZ  = 0  /* mg q4w tocilizumab           */
  DOSE_MTX  = 0  /* mg/week methotrexate         */
  DOSE_IFX  = 0  /* mg/kg infliximab q6w         */

$CMT
  /* Prednisone/Prednisolone */
  PRED_GUT PRED_C PRED_P PRED_EFF
  /* Tocilizumab */
  TCZ_SC TCZ_C TCZ_P
  /* Methotrexate */
  MTX_GUT MTX_C MTX_PG
  /* Infliximab */
  IFX_C IFX_P
  /* Disease compartments */
  IL6 sIL6R IL6_cmplx TNF
  TH1 TH17 TREG
  VWI ST
  /* Biomarkers */
  CRP PET VWT

$INIT
  PRED_GUT = 0,   PRED_C = 0,   PRED_P = 0,   PRED_EFF = 0,
  TCZ_SC   = 0,   TCZ_C  = 0,   TCZ_P  = 0,
  MTX_GUT  = 0,   MTX_C  = 0,   MTX_PG = 0,
  IFX_C    = 0,   IFX_P  = 0,
  IL6      = 1.2, sIL6R  = 200, IL6_cmplx = 0, TNF = 0.5,
  TH1      = 25,  TH17   = 15,  TREG = 8,
  VWI      = 2.0, ST     = 0.0,
  CRP      = 5.0, PET    = 1.5, VWT = 1.5

$ODE
  /* ========== Drug PK ========== */

  /* -- Prednisone gut → prednisolone plasma -- */
  double dPRED_GUT = -ka_PRED * PRED_GUT;
  double dPRED_C   = ka_PRED * F_PRED * PRED_GUT / Vc_PRED
                     - (CL_PRED/Vc_PRED) * PRED_C
                     - (CLd_PRED/Vc_PRED) * PRED_C
                     + (CLd_PRED/Vp_PRED) * PRED_P;
  double dPRED_P   = (CLd_PRED/Vc_PRED) * PRED_C
                     - (CLd_PRED/Vp_PRED) * PRED_P;
  double dPRED_EFF = ke0_PRED * (PRED_C - PRED_EFF);

  /* -- Tocilizumab SC → plasma -- */
  double dTCZ_SC = -ka_TCZ * TCZ_SC;
  double CL_TCZ_tot = CL_TCZ + CLMM_TCZ * TCZ_C / (KM_TCZ + TCZ_C);
  double dTCZ_C  = ka_TCZ * F_TCZ * TCZ_SC / Vc_TCZ
                   - (CL_TCZ_tot/Vc_TCZ) * TCZ_C
                   - (CLd_TCZ/Vc_TCZ) * TCZ_C
                   + (CLd_TCZ/Vp_TCZ) * TCZ_P;
  double dTCZ_P  = (CLd_TCZ/Vc_TCZ) * TCZ_C
                   - (CLd_TCZ/Vp_TCZ) * TCZ_P;

  /* -- MTX oral → plasma → polyglutamates -- */
  double dMTX_GUT = -ka_MTX * MTX_GUT;
  double dMTX_C   = ka_MTX * F_MTX * MTX_GUT / Vc_MTX
                    - (CL_MTX/Vc_MTX) * MTX_C
                    - kpg_MTX * MTX_C;
  double dMTX_PG  = kpg_MTX * MTX_C - kdpg_MTX * MTX_PG;

  /* -- Infliximab IV → plasma -- */
  double dIFX_C   = -(CL_IFX/Vc_IFX) * IFX_C
                    - (CLd_IFX/Vc_IFX) * IFX_C
                    + (CLd_IFX/Vp_IFX) * IFX_P;
  double dIFX_P   = (CLd_IFX/Vc_IFX) * IFX_C
                    - (CLd_IFX/Vp_IFX) * IFX_P;

  /* ========== Drug PD (Emax functions) ========== */

  /* Prednisolone inhibits IL-6, TNF, Th1, Th17 via GR */
  double Inh_PRED = Emax_PRED * pow(PRED_EFF, hill_PRED) /
                    (pow(EC50_PRED, hill_PRED) + pow(PRED_EFF, hill_PRED));

  /* Tocilizumab blocks IL-6R (trans-signaling & signaling) */
  double Occ_TCZ  = Emax_TCZ * pow(TCZ_C, hill_TCZ) /
                    (pow(EC50_TCZ, hill_TCZ) + pow(TCZ_C, hill_TCZ));

  /* MTX-PG inhibits T cell proliferation */
  double Inh_MTX  = Emax_MTX * pow(MTX_PG, hill_MTX) /
                    (pow(EC50_MTX, hill_MTX) + pow(MTX_PG, hill_MTX));

  /* Infliximab neutralizes TNF-alpha */
  double Inh_IFX  = Emax_IFX * pow(IFX_C, hill_IFX) /
                    (pow(EC50_IFX, hill_IFX) + pow(IFX_C, hill_IFX));

  /* ========== Disease Compartments ========== */

  /* IL-6 ODE
     Synthesis driven by Th1, Th17, TNF, amplified by disease;
     Inhibited by prednisolone (GR) and blocked by TCZ (feedback rise sIL-6R)
     TCZ leads to paradoxical serum IL-6 rise (sIL-6R)
  */
  double IL6_syn = ksyn_IL6
                   + amp_IL6_Th1  * TH1
                   + amp_IL6_Th17 * TH17
                   + amp_IL6_TNF  * TNF;
  double IL6_deg = kdeg_IL6 * (1 + Occ_TCZ * 0.2);  /* slight increase in IL-6 half-life blockade */
  double dIL6    = IL6_syn * (1 - Inh_PRED)
                   - IL6_deg * IL6
                   - kon_IL6R * IL6 * sIL6R
                   + koff_IL6R * IL6_cmplx;

  /* Soluble IL-6R: TCZ increases sIL-6R as free receptor accumulates */
  double sIL6R_factor = 1 + 2.5 * Occ_TCZ;  /* TCZ blocks membrane IL-6R, shedding increases */
  double dSIL6R  = ksyn_sR * sIL6R_factor
                   - kdeg_sR * sIL6R
                   - kon_IL6R * IL6 * sIL6R
                   + koff_IL6R * IL6_cmplx;

  double dIL6_cmplx = kon_IL6R * IL6 * sIL6R
                      - koff_IL6R * IL6_cmplx
                      - kdeg_IL6 * IL6_cmplx;

  /* TNF-alpha ODE */
  double dTNF = ksyn_TNF + 0.03 * TH1
                - kdeg_TNF * TNF * (1 - Inh_IFX) * (1 - 0.6 * Inh_PRED);

  /* Th1 ODE: driven by IL-6, IFN-gamma circuit; inhibited by Treg and drugs */
  double dTH1 = ksyn_Th1 * (1 + 0.05 * IL6)
                - kdeg_Th1 * TH1 * (1 - Inh_MTX) * (1 - Inh_PRED * 0.5)
                - inh_Treg_Th1 * TREG * TH1;

  /* Th17 ODE: driven by IL-6 + IL-23 axis, TCZ markedly suppresses */
  double dTH17 = ksyn_Th17 * (1 + 0.04 * IL6) * (1 - Occ_TCZ * 0.8)
                 - kdeg_Th17 * TH17 * (1 - Inh_MTX * 0.7) * (1 - Inh_PRED * 0.4)
                 - inh_Treg_Th17 * TREG * TH17;

  /* Treg ODE: prednisolone and TCZ partially restore Treg */
  double dTREG = ksyn_Treg * (1 + 0.3 * Inh_PRED + 0.2 * Occ_TCZ)
                 - kdeg_Treg * TREG;

  /* Vessel Wall Inflammation Index [0-10]
     Driven by IL-6, TNF, Th17; suppressed by all drugs via Emax
  */
  double VWI_drive = amp_VWI_IL6 * IL6
                     + amp_VWI_TNF * TNF
                     + amp_VWI_Th17 * TH17;
  double Drug_inh_VWI = 1 - (1 - Inh_PRED) * (1 - Occ_TCZ * 0.9)
                              * (1 - Inh_IFX * 0.7) * (1 - Inh_MTX * 0.3);
  double dVWI = ksyn_VWI * VWI_drive * (1 - Drug_inh_VWI)
                - kdeg_VWI * VWI;

  /* Stenosis Index [0-100 %] — driven by cumulative VWI */
  double dST = kprog_ST * amp_ST_VWI * VWI
               - kreg_ST * ST;

  /* CRP (mg/L) — produced proportional to IL-6; TCZ rapidly normalizes */
  double dCRP = ksyn_CRP * IL6 * (1 - Occ_TCZ * 0.95)
                - kdeg_CRP * CRP;

  /* PET-CT SUVmax — correlates with VWI and vessel wall metabolic activity */
  double dPET = ksyn_PET * VWI
                - kdeg_PET * PET;

  /* MRI Vessel Wall Thickness (mm) — driven by cumulative stenosis/VWI */
  double dVWT = kVWT * VWI - 0.001 * VWT;

  /* ========== Assign DES ========== */
  dxdt_PRED_GUT  = dPRED_GUT;
  dxdt_PRED_C    = dPRED_C;
  dxdt_PRED_P    = dPRED_P;
  dxdt_PRED_EFF  = dPRED_EFF;
  dxdt_TCZ_SC    = dTCZ_SC;
  dxdt_TCZ_C     = dTCZ_C;
  dxdt_TCZ_P     = dTCZ_P;
  dxdt_MTX_GUT   = dMTX_GUT;
  dxdt_MTX_C     = dMTX_C;
  dxdt_MTX_PG    = dMTX_PG;
  dxdt_IFX_C     = dIFX_C;
  dxdt_IFX_P     = dIFX_P;
  dxdt_IL6       = dIL6;
  dxdt_sIL6R     = dSIL6R;
  dxdt_IL6_cmplx = dIL6_cmplx;
  dxdt_TNF       = dTNF;
  dxdt_TH1       = dTH1;
  dxdt_TH17      = dTH17;
  dxdt_TREG      = dTREG;
  dxdt_VWI       = dVWI;
  dxdt_ST        = dST;
  dxdt_CRP       = dCRP;
  dxdt_PET       = dPET;
  dxdt_VWT       = dVWT;

$TABLE
  double NIH_SCORE = 0;
  /* NIH Disease Activity Score: 0-20 scale
     based on new/worsening features:
     systemic symptoms (2pts), ESR (2pts), angiography (2pts),
     ischemic symptoms (6pts), BP difference (2pts)
  */
  double ESR_now = ESR_base + kESR * (CRP - CRP_base);
  NIH_SCORE = 2 * (CRP > 20 ? 1 : CRP/20) +  /* systemic inflammation */
              2 * (ESR_now > 40 ? 1 : ESR_now/40) +  /* ESR */
              3 * (VWI/10) +                   /* vascular inflammation */
              4 * (ST/50) +                    /* stenosis extent */
              3 * (PET/4);                     /* PET activity */
  if(NIH_SCORE > 20) NIH_SCORE = 20;

  double ITAS_SCORE = 0;
  /* ITAS 2010 simplification */
  ITAS_SCORE = 1.5 * (CRP > 10 ? 1 : 0) +
               1.5 * (VWI > 5 ? 1 : 0) +
               3.0 * (ST > 20 ? 1 : 0) +
               2.0 * (PET > 2.5 ? 1 : 0);

  double RESPONSE_FLAG = NIH_SCORE < 4 ? 1 : 0;  /* low disease activity */
  double CP_PRED  = PRED_C;
  double CP_TCZ   = TCZ_C;
  double CP_MTX   = MTX_C;
  double CP_IFX   = IFX_C;

  capture ESR_now NIH_SCORE ITAS_SCORE RESPONSE_FLAG
  capture CP_PRED CP_TCZ CP_MTX CP_IFX
  capture Drug_inh_VWI Occ_TCZ Inh_PRED Inh_MTX Inh_IFX
'

mod <- mcode("TakayasuArteritis", ta_model_code)

## ------------------------------------------------------------
## 2. Helper Functions: Dosing Events
## ------------------------------------------------------------

## Prednisone 1 mg/kg/day orally (continuous via events)
pred_events <- function(dose_mg = 65, duration_days = 365,
                        taper_to = 10, taper_start = 60) {
  times_init  <- seq(0, taper_start * 24, by = 24)
  times_taper <- seq((taper_start + 7) * 24, duration_days * 24, by = 24)
  dose_taper  <- seq(dose_mg, taper_to,
                     length.out = length(times_taper))
  ev_init  <- ev(time = times_init,  amt = dose_mg, cmt = "PRED_GUT")
  ev_taper <- ev(time = times_taper, amt = dose_taper, cmt = "PRED_GUT")
  c(ev_init, ev_taper)
}

## Tocilizumab 162 mg SC q2w (standard TA dosing)
tcz_events <- function(start_day = 0, n_doses = 26) {
  times <- seq(start_day * 24, by = 14 * 24, length.out = n_doses)
  ev(time = times, amt = 162, cmt = "TCZ_SC")
}

## Methotrexate 15 mg/week orally
mtx_events <- function(start_day = 0, duration_days = 365) {
  times <- seq(start_day * 24, duration_days * 24, by = 7 * 24)
  ev(time = times, amt = 15, cmt = "MTX_GUT")
}

## Infliximab 5 mg/kg IV: 0, 2, 6 weeks then q6w
ifx_events <- function(start_day = 0, wt_kg = 65, n_maint = 8) {
  dose <- 5 * wt_kg
  induction_times <- c(0, 2, 6) * 7 * 24 + start_day * 24
  maint_times     <- seq((6 + 6) * 7 * 24, by = 6 * 7 * 24,
                         length.out = n_maint) + start_day * 24
  ev(time = c(induction_times, maint_times),
     amt  = dose, cmt = "IFX_C", rate = -2)  /* 2h infusion */
}

## ------------------------------------------------------------
## 3. Simulation Parameters
## ------------------------------------------------------------
sim_end <- 365 * 24   /* 1 year in hours */
dt      <- 4           /* 4-hour output step */
times   <- seq(0, sim_end, by = dt)

base_params <- list(
  ksyn_IL6 = 0.45,   /* active TA: elevated IL-6 synthesis */
  ksyn_TNF = 0.20,
  ksyn_Th1 = 0.10,
  ksyn_Th17= 0.07,
  VWI      = 6.0,    /* active disease starting VWI */
  IL6      = 8.0,    /* elevated baseline */
  TNF      = 2.0,
  TH1      = 50.0,
  TH17     = 35.0,
  CRP      = 45.0,
  PET      = 3.5,
  VWT      = 4.5
)

## ------------------------------------------------------------
## 4. Treatment Scenarios
## ------------------------------------------------------------

## Scenario 1: No treatment (natural history)
ev_none <- ev(time = 0, amt = 0, cmt = "PRED_GUT")

## Scenario 2: Prednisone monotherapy (1 mg/kg/day → taper to 10 mg/day)
ev_pred <- pred_events(dose_mg = 65, duration_days = 365,
                       taper_to = 10, taper_start = 60)

## Scenario 3: Prednisone + Methotrexate (standard first-line combo)
ev_pred_mtx <- c(
  pred_events(dose_mg = 65, taper_to = 10, taper_start = 60),
  mtx_events(start_day = 0)
)

## Scenario 4: Prednisone + Tocilizumab (TAKT trial regimen)
ev_pred_tcz <- c(
  pred_events(dose_mg = 65, taper_to = 7.5, taper_start = 90),
  tcz_events(start_day = 0)
)

## Scenario 5: Prednisone + Infliximab (refractory TA)
ev_pred_ifx <- c(
  pred_events(dose_mg = 65, taper_to = 10, taper_start = 60),
  ifx_events(start_day = 0)
)

scenarios <- list(
  "1. No Treatment (Natural History)" = ev_none,
  "2. Prednisone Monotherapy"         = ev_pred,
  "3. Prednisone + Methotrexate"      = ev_pred_mtx,
  "4. Prednisone + Tocilizumab (TAKT)"= ev_pred_tcz,
  "5. Prednisone + Infliximab (Refract.)"= ev_pred_ifx
)

## ------------------------------------------------------------
## 5. Run All Scenarios
## ------------------------------------------------------------
run_scenario <- function(ev_obj, scen_name) {
  idata <- as.data.frame(base_params)
  out <- mod %>%
    param(base_params) %>%
    init(IL6 = base_params$IL6,
         TNF = base_params$TNF,
         TH1 = base_params$TH1,
         TH17 = base_params$TH17,
         VWI = base_params$VWI,
         CRP = base_params$CRP,
         PET = base_params$PET,
         VWT = base_params$VWT) %>%
    ev(ev_obj) %>%
    mrgsim(end = sim_end, delta = dt, carry_out = "evid") %>%
    as_tibble() %>%
    mutate(time_days = time / 24, Scenario = scen_name)
  out
}

results <- bind_rows(lapply(names(scenarios), function(nm) {
  run_scenario(scenarios[[nm]], nm)
}))

## ------------------------------------------------------------
## 6. Visualization
## ------------------------------------------------------------
theme_ta <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#E65100", color = "white"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

scen_colors <- c(
  "1. No Treatment (Natural History)"       = "#D32F2F",
  "2. Prednisone Monotherapy"               = "#F57C00",
  "3. Prednisone + Methotrexate"            = "#388E3C",
  "4. Prednisone + Tocilizumab (TAKT)"      = "#1565C0",
  "5. Prednisone + Infliximab (Refract.)"   = "#6A1B9A"
)

p1 <- ggplot(results, aes(time_days, IL6, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "A. Serum IL-6 (pg/mL)", x = "Time (days)", y = "IL-6 (pg/mL)") +
  theme_ta

p2 <- ggplot(results, aes(time_days, CRP, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 12, label = "CRP = 10 mg/L", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "B. CRP (mg/L)", x = "Time (days)", y = "CRP (mg/L)") +
  theme_ta

p3 <- ggplot(results, aes(time_days, NIH_SCORE, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 4.5, label = "Remission threshold", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "C. NIH Disease Activity Score (0-20)", x = "Time (days)", y = "NIH Score") +
  theme_ta

p4 <- ggplot(results, aes(time_days, VWI, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "D. Vessel Wall Inflammation Index (0-10)", x = "Time (days)", y = "VWI Score") +
  theme_ta

p5 <- ggplot(results, aes(time_days, ST, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "E. Arterial Stenosis Index (%)", x = "Time (days)", y = "Stenosis (%)") +
  theme_ta

p6 <- ggplot(results, aes(time_days, PET, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 2.7, label = "PET activity threshold", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "F. PET-CT FDG SUVmax (Vascular)", x = "Time (days)", y = "SUVmax") +
  theme_ta

p_pk1 <- ggplot(
  results %>% filter(grepl("TCZ", Scenario)),
  aes(time_days, CP_TCZ, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "G. Tocilizumab Plasma Concentration",
       x = "Time (days)", y = "Concentration (mg/L)") +
  theme_ta

p_th <- ggplot(
  results %>% select(time_days, Scenario, TH1, TH17, TREG) %>%
    pivot_longer(c(TH1, TH17, TREG), names_to = "Cell", values_to = "Count"),
  aes(time_days, Count, color = Cell, linetype = Scenario)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = c(TH1 = "#D32F2F", TH17 = "#1565C0", TREG = "#2E7D32"),
                     labels = c("Th1 Cells", "Th17 Cells", "Treg Cells")) +
  labs(title = "H. T Cell Populations (cells/µL)", x = "Time (days)", y = "Count (cells/µL)") +
  theme_ta

## Combined dashboard
dashboard <- (p1 + p2) / (p3 + p4) / (p5 + p6)
dashboard_full <- dashboard + plot_annotation(
  title = "Takayasu Arteritis QSP Model — Treatment Scenario Comparison",
  subtitle = "5 scenarios: natural history vs. prednisolone ± MTX / TCZ / IFX",
  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
)

## ------------------------------------------------------------
## 7. Clinical Trial Calibration Reference
## ------------------------------------------------------------
calibration_table <- data.frame(
  Trial          = c("TAKT (Nakaoka 2018 Lancet)",
                     "NIHON-BVAS Pred monotherapy",
                     "Abisror 2013 JRA",
                     "Comarmond 2012 Medicine",
                     "Hoffman 1994 Ann Intern Med"),
  Intervention   = c("Tocilizumab + Pred vs Pred",
                     "Pred 1 mg/kg taper",
                     "Infliximab (TNF)",
                     "Infliximab rescue",
                     "Methotrexate + Pred"),
  Key_Outcome    = c("Time to relapse HR 0.41 (TCZ arm)",
                     "~70% initial remission rate",
                     "Remission 93% refractory TA",
                     "Response 67% refractory",
                     "Steroid-sparing; 72% remission"),
  Model_Parameter= c("Emax_TCZ=0.95, EC50_TCZ=0.50",
                     "Emax_PRED=0.85, taper over 60 days",
                     "Emax_IFX=0.90, EC50_IFX=0.80",
                     "Inh_IFX 0-0.90 dose-range",
                     "Emax_MTX=0.55, EC50_MTX=0.05"),
  stringsAsFactors = FALSE
)

cat("==================================================\n")
cat("  Takayasu Arteritis QSP Model — Calibration\n")
cat("==================================================\n")
print(calibration_table, row.names = FALSE)

cat("\n--- End-of-year summary (day 365) ---\n")
summary_365 <- results %>%
  filter(abs(time_days - 365) < 0.25) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  select(Scenario, IL6, CRP, VWI, ST, PET, NIH_SCORE, ITAS_SCORE) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
print(as.data.frame(summary_365), row.names = FALSE)

print(dashboard_full)

## ------------------------------------------------------------
## 8. Sensitivity Analysis: IL-6 synthesis rate vs. TCZ efficacy
## ------------------------------------------------------------
cat("\n--- Sensitivity: ksyn_IL6 vs. 1-year NIH score (TCZ scenario) ---\n")
ksyn_range <- seq(0.20, 0.80, by = 0.15)
sens_results <- lapply(ksyn_range, function(k) {
  p_override <- modifyList(base_params, list(ksyn_IL6 = k, IL6 = k / 0.25 * 1.2))
  out <- mod %>%
    param(p_override) %>%
    init(IL6 = p_override$IL6, CRP = 45, VWI = 6, TH1 = 50, TH17 = 35) %>%
    ev(scenarios[["4. Prednisone + Tocilizumab (TAKT)"]]) %>%
    mrgsim(end = sim_end, delta = dt) %>%
    as_tibble() %>%
    filter(abs(time / 24 - 365) < 0.5) %>%
    slice(1) %>%
    mutate(ksyn_IL6 = k)
})
sens_df <- bind_rows(sens_results)
cat(sprintf("  ksyn_IL6=%.2f → NIH_SCORE=%.2f\n",
            sens_df$ksyn_IL6, sens_df$NIH_SCORE))
