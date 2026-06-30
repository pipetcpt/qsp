# =============================================================================
# Cluster Headache (CH) — Quantitative Systems Pharmacology (QSP) ODE model
# mrgsolve format
#
# Scope
#   - Hypothalamic generator (circadian gate) driving trigeminovascular attacks
#   - Attack-rate state controlled by CGRP & PACAP tone; modulated by treatments
#   - PK/PD of 7 drugs : sumatriptan SC, zolmitriptan IN, verapamil PR,
#                        lithium, topiramate, galcanezumab (CGRP mAb), prednisone
#   - Devices/non-PK : O2 (effect compartment), GON block (decaying effect)
#   - Clinical endpoints : attacks/week, mean weekly attack rate, response
#
# Citations used to anchor parameters (see ch_references.md)
#   - Goadsby et al. NEJM 2019 (galcanezumab eCH) -- mean attack reduction 8.7 vs 5.2
#   - Cohen et al. JAMA 2009 (high-flow O2) -- 78% pain-free @15 min
#   - Ekbom et al. NEJM 1991 (SC sumatriptan) -- 74% relief @15 min
#   - Leone et al. Neurology 2000 (verapamil RCT) -- 240 mg/d efficacy
#   - Steiner et al. Cephalalgia 1997 (lithium vs verapamil)
#   - Wei DY, Goadsby PJ. Curr Pain Headache Rep 2021 (pharmacology review)
#
# Compartments (23 main ODE + auxiliary):
#   1  SUMA_DEPOT     sumatriptan SC depot                        (mg)
#   2  SUMA_CENT      sumatriptan central                         (mg)
#   3  ZOL_DEPOT      zolmitriptan IN depot                        (mg)
#   4  ZOL_CENT       zolmitriptan central                         (mg)
#   5  ZOL_PER        zolmitriptan peripheral                      (mg)
#   6  VERA_DEPOT     verapamil PR depot                           (mg)
#   7  VERA_CENT      verapamil central                            (mg)
#   8  VERA_PER       verapamil peripheral                         (mg)
#   9  LI_DEPOT       lithium depot                                (mg)
#   10 LI_CENT        lithium central                              (mg)
#   11 TOPI_DEPOT     topiramate depot                             (mg)
#   12 TOPI_CENT      topiramate central                           (mg)
#   13 GALCA_SC       galcanezumab SC depot                        (mg)
#   14 GALCA_CENT     galcanezumab central                         (mg)
#   15 PRED_DEPOT     prednisone PO depot                          (mg)
#   16 PRED_CENT      prednisolone central                         (mg)
#   17 HYPO_DRIVE     hypothalamic excitability (0=quiet 1=active)
#   18 CGRP           trigeminal CGRP tone (a.u., reference 1.0)
#   19 PACAP          PACAP / VIP tone (a.u.)
#   20 PIAL           pial vasodilation effect site (a.u.)
#   21 ATTACK_HZ      instantaneous attack hazard /h
#   22 CUM_ATTACKS    cumulative attacks count
#   23 BOUT_TIMER     within-bout time (d) for episodic CH gate
#   24 O2_EFFECT      bolus O2 effect compartment (transient)
#   25 GON_EFF        GON block effect compartment (decays over weeks)
#
# Five+ scenarios (see "events" examples at bottom):
#   S0 No treatment
#   S1 O2 + sumatriptan abortive
#   S2 Verapamil 240->480 mg/d preventive
#   S3 Verapamil + lithium (chronic CH)
#   S4 Galcanezumab 300 mg SC monthly
#   S5 Prednisone bridge + verapamil (induction)
#   S6 GON block + verapamil
# =============================================================================

library(mrgsolve)
library(tidyverse)

ch_code <- '
$PROB Cluster Headache QSP — hypothalamic-trigeminovascular axis with 7 drugs

$PARAM @annotated
// ---- sumatriptan SC PK
KA_SUMA   :  2.5   : SC absorption rate /h
CL_SUMA   :  18    : Sumatriptan CL (L/h)
V_SUMA    :  9     : Central V (L)
EC50_SUMA :  20    : Conc giving 50% acute abortive effect (ng/mL)
EMAX_SUMA :  0.92  : Max fraction of attack hazard reduction

// ---- zolmitriptan IN PK
KA_ZOL    :  1.5   : IN absorption /h
CL_ZOL    :  9     : CL (L/h)
V_ZOL     :  90    : Central V (L)
Q_ZOL     :  6     : intercompartmental (L/h)
VP_ZOL    :  60    : peripheral V (L)
EMAX_ZOL  :  0.80  : abortive Emax
EC50_ZOL  :  6     : ng/mL

// ---- verapamil PR PK (extended release, racemic)
KA_VERA   :  0.20  : absorption /h
CL_VERA   :  60    : apparent CL (L/h)
V_VERA    :  300   : central V (L)
Q_VERA    :  35    : intercompartmental (L/h)
VP_VERA   :  500   : peripheral V (L)
EC50_VERA :  90    : ng/mL effective for prevention
EMAX_VERA :  0.55  : max fractional attack-rate reduction

// ---- lithium PK (single compartment, renal)
KA_LI     :  1.5   : absorption /h
CL_LI     :  1.5   : CL (L/h, ≈25 mL/min)
V_LI      :  45    : V (L)
EC50_LI   :  0.6   : mEq/L
EMAX_LI   :  0.50  : preventive

// ---- topiramate PK
KA_TOPI   :  0.6   : /h
CL_TOPI   :  1.2   : L/h
V_TOPI    :  60    : L
EC50_TOPI :  5     : µg/mL
EMAX_TOPI :  0.35

// ---- galcanezumab PK (SC mAb, FcRn recycled, linear)
KA_GALCA  :  0.012 : SC absorption /h (≈ F 0.75 t_abs days)
CL_GALCA  :  0.008 : L/h ≈ 0.19 L/d
V_GALCA   :  7.5   : L
EC50_GALCA:  1.5   : µg/mL drug at site, IC50 of CGRP binding
EMAX_GALCA:  0.65  : max attack-rate reduction (eCH ECH trial scale)

// ---- prednisolone PK (active form from prednisone)
KA_PRED   :  2.5   : /h (rapid)
CL_PRED   :  6     : L/h
V_PRED    :  40    : L
EC50_PRED :  20    : ng/mL
EMAX_PRED :  0.70  : transitional bridge

// ---- disease parameters (hypothalamic generator)
BASE_HYPO :  0.05  : tonic hypothalamic drive baseline (quiet remission)
BOUT_AMP  :  0.95  : in-bout drive amplitude
KIN_HYPO  :  0.020 : drive build-up /h
KOUT_HYPO :  0.010 : drive decay /h
KIN_CGRP  :  0.15  : CGRP production rate driven by hypothalamus
KOUT_CGRP :  0.30  : CGRP elimination /h
KIN_PACAP :  0.08  : PACAP production /h
KOUT_PACAP:  0.20  : PACAP elimination /h
PIAL_HALF :  3     : pial response half-time (h)
ATTACK_K0 :  0.001 : baseline attack hazard /h (~ 0.024 /d ~ no attacks)
ATTACK_KSAT:  0.18 : max hazard /h ( ≈ 4 attacks/day cap )
CGRP_SET  :  1.0   : reference CGRP tone (mAb at this point ≈50% effect)
PACAP_SET :  1.0   : reference PACAP tone
KO2_ON    :  3.0   : O2 effect time-constant on (/h)
KO2_OFF   :  6.0   : O2 effect decay /h
GON_HL    :  336   : GON-block effect half-life (h ~ 2 weeks)

// ---- circadian / bout gate
CIRC_AMP  :  0.35  : 24-h amplitude on hazard
CIRC_PHASE:  3     : night peak (~03:00 hr-of-day)
BOUT_LEN  :  42    : days in bout
REMISSION :  240   : days in remission
BOUT_ON   :  1     : start in active bout (1) or not (0)

// ---- patient covariates
CrCL      :  100   : creatinine clearance (mL/min)
WT        :  78    : body weight (kg)
SEX       :  1     : 1=M 0=F
SMOKER    :  1     : 1 active 0 no
CHRONIC   :  0     : 1 chronic CH (no remission)

$CMT @annotated
SUMA_DEPOT  : sumatriptan SC depot
SUMA_CENT   : sumatriptan central
ZOL_DEPOT   : zolmitriptan IN depot
ZOL_CENT    : zolmitriptan central
ZOL_PER     : zolmitriptan peripheral
VERA_DEPOT  : verapamil PR depot
VERA_CENT   : verapamil central
VERA_PER    : verapamil peripheral
LI_DEPOT    : lithium depot
LI_CENT     : lithium central
TOPI_DEPOT  : topiramate depot
TOPI_CENT   : topiramate central
GALCA_SC    : galcanezumab SC depot
GALCA_CENT  : galcanezumab central
PRED_DEPOT  : prednisone depot
PRED_CENT   : prednisolone central
HYPO_DRIVE  : hypothalamic excitability
CGRP        : trigeminal CGRP tone (a.u.)
PACAP       : PACAP tone (a.u.)
PIAL        : pial vasodilation effect site
ATTACK_HZ   : instantaneous attack hazard /h
CUM_ATTACKS : cumulative attacks
BOUT_TIMER  : days within current bout
O2_EFFECT   : O2 effect site
GON_EFF     : GON block effect site

$MAIN
double cl_li = CL_LI * (CrCL/100);   // renal scaling for lithium
double bout_drive = BOUT_ON * BOUT_AMP;
// chronic CH never enters remission gate
if(CHRONIC == 1) bout_drive = BOUT_AMP;

$ODE
// concentrations
double CP_SUMA  = SUMA_CENT * 1000 / V_SUMA;          // ng/mL (mg->µg & ng adj.)
double CP_ZOL   = ZOL_CENT  * 1000 / V_ZOL;           // ng/mL
double CP_VERA  = VERA_CENT * 1000 / V_VERA;          // ng/mL
double CP_LI    = LI_CENT   / V_LI;                   // mEq/L (mg -> mEq approx)
double CP_TOPI  = TOPI_CENT / V_TOPI;                 // µg/mL
double CP_GALCA = GALCA_CENT / V_GALCA;               // µg/mL
double CP_PRED  = PRED_CENT * 1000 / V_PRED;          // ng/mL

// ---- PK ODEs
dxdt_SUMA_DEPOT = -KA_SUMA  * SUMA_DEPOT;
dxdt_SUMA_CENT  =  KA_SUMA  * SUMA_DEPOT - (CL_SUMA / V_SUMA) * SUMA_CENT;

dxdt_ZOL_DEPOT  = -KA_ZOL  * ZOL_DEPOT;
dxdt_ZOL_CENT   =  KA_ZOL  * ZOL_DEPOT
                  - (CL_ZOL / V_ZOL) * ZOL_CENT
                  - (Q_ZOL  / V_ZOL) * ZOL_CENT
                  + (Q_ZOL  / VP_ZOL) * ZOL_PER;
dxdt_ZOL_PER    =  (Q_ZOL / V_ZOL) * ZOL_CENT - (Q_ZOL / VP_ZOL) * ZOL_PER;

dxdt_VERA_DEPOT = -KA_VERA * VERA_DEPOT;
dxdt_VERA_CENT  =  KA_VERA * VERA_DEPOT
                  - (CL_VERA / V_VERA) * VERA_CENT
                  - (Q_VERA  / V_VERA) * VERA_CENT
                  + (Q_VERA  / VP_VERA) * VERA_PER;
dxdt_VERA_PER   =  (Q_VERA / V_VERA) * VERA_CENT - (Q_VERA / VP_VERA) * VERA_PER;

dxdt_LI_DEPOT   = -KA_LI * LI_DEPOT;
dxdt_LI_CENT    =  KA_LI * LI_DEPOT - (cl_li / V_LI) * LI_CENT;

dxdt_TOPI_DEPOT = -KA_TOPI * TOPI_DEPOT;
dxdt_TOPI_CENT  =  KA_TOPI * TOPI_DEPOT - (CL_TOPI / V_TOPI) * TOPI_CENT;

dxdt_GALCA_SC   = -KA_GALCA * GALCA_SC;
dxdt_GALCA_CENT =  KA_GALCA * GALCA_SC - (CL_GALCA / V_GALCA) * GALCA_CENT;

dxdt_PRED_DEPOT = -KA_PRED * PRED_DEPOT;
dxdt_PRED_CENT  =  KA_PRED * PRED_DEPOT - (CL_PRED / V_PRED) * PRED_CENT;

// ---- circadian modulation (24-h cycle, fits peak in early hours)
double hr     = fmod(SOLVERTIME, 24.0);
double circ   = 1.0 + CIRC_AMP * cos(2*M_PI*(hr - CIRC_PHASE) / 24.0);

// ---- bout timer (episodic vs chronic)
dxdt_BOUT_TIMER = 1.0/24.0;   // accrue days; reset via $TABLE if needed (here simple)

// ---- hypothalamic drive (target circ * bout_drive)
double hypo_target = BASE_HYPO + circ * bout_drive;
dxdt_HYPO_DRIVE = KIN_HYPO * (hypo_target - HYPO_DRIVE);

// ---- drug effects on hypothalamic drive (preventive)
double E_vera   = EMAX_VERA  * CP_VERA  / (CP_VERA  + EC50_VERA);
double E_li     = EMAX_LI    * CP_LI    / (CP_LI    + EC50_LI);
double E_topi   = EMAX_TOPI  * CP_TOPI  / (CP_TOPI  + EC50_TOPI);
double E_galca  = EMAX_GALCA * CP_GALCA / (CP_GALCA + EC50_GALCA);
double E_pred   = EMAX_PRED  * CP_PRED  / (CP_PRED  + EC50_PRED);
double E_gon    = GON_EFF;                                // 0..1 fraction
double preventive = 1 - (1-E_vera)*(1-E_li)*(1-E_topi)*(1-E_galca)*(1-E_pred)*(1-E_gon);
if(preventive > 0.95) preventive = 0.95;

// ---- trigeminovascular tone driven by hypothalamus, dampened by preventive
double drive_eff = HYPO_DRIVE * (1 - preventive);
dxdt_CGRP  = KIN_CGRP  * drive_eff - KOUT_CGRP  * CGRP;
dxdt_PACAP = KIN_PACAP * drive_eff - KOUT_PACAP * PACAP;

// ---- pial response (effect compartment of CGRP+PACAP)
double pial_drive = 0.6*CGRP + 0.4*PACAP;
dxdt_PIAL  = (log(2)/PIAL_HALF) * (pial_drive - PIAL);

// ---- O2 effect decay
dxdt_O2_EFFECT = -KO2_OFF * O2_EFFECT;

// ---- GON block decay (single exponential, k = ln2/HL)
dxdt_GON_EFF = -(log(2)/GON_HL) * GON_EFF;

// ---- attack hazard (CGRP-driven, capped, acute abortive abates)
double E_suma_acute = EMAX_SUMA * CP_SUMA / (CP_SUMA + EC50_SUMA);
double E_zol_acute  = EMAX_ZOL  * CP_ZOL  / (CP_ZOL  + EC50_ZOL);
double E_O2_acute   = O2_EFFECT;                          // 0..1
double acute_abort  = 1 - (1 - E_suma_acute)*(1 - E_zol_acute)*(1 - E_O2_acute);
if(acute_abort > 0.99) acute_abort = 0.99;

double hz_raw = ATTACK_KSAT * (PIAL / (PIAL + 1.0));
double hz_now = (ATTACK_K0 + hz_raw) * (1 - acute_abort);
dxdt_ATTACK_HZ = (hz_now - ATTACK_HZ) * 5.0;   // smoothing

dxdt_CUM_ATTACKS = ATTACK_HZ;                  // integrate hazard

$TABLE
capture ConcSuma    = CP_SUMA;
capture ConcZol     = CP_ZOL;
capture ConcVera    = CP_VERA;
capture ConcLi      = CP_LI;
capture ConcTopi    = CP_TOPI;
capture ConcGalca   = CP_GALCA;
capture ConcPred    = CP_PRED;
capture HypoDrive   = HYPO_DRIVE;
capture CGRPtone    = CGRP;
capture PACAPtone   = PACAP;
capture Pialtone    = PIAL;
capture HazardPerH  = ATTACK_HZ;
capture AttacksWeek = ATTACK_HZ * 168.0;
capture Preventive  = preventive;
capture AcuteAbort  = acute_abort;
'

# ---- compile ----------------------------------------------------------------
ch_mod <- mcode("cluster_headache_qsp", ch_code)

# ============================================================================
# Helper builders for the 6 standard scenarios
# ============================================================================

build_scenario <- function(scenario = "S0") {
  ev <- ev(amt = 0, cmt = 1, time = 0)   # placeholder

  if (scenario == "S0") {
    # no drug
    return(ev)
  }
  if (scenario == "S1") {
    # 1 attack at day 1 morning: 6 mg sumatriptan SC + O2 bolus (set effect to 0.78)
    e1 <- ev(time = 24 + 4,  amt = 6,    cmt = "SUMA_DEPOT")
    e2 <- ev(time = 24 + 4,  amt = 0.78, cmt = "O2_EFFECT", evid = 1)  # bolus
    return(c(e1, e2))
  }
  if (scenario == "S2") {
    # verapamil 240 mg PR BID for 14 d, then 480 mg/d (i.e. 240 BID)
    e1 <- ev(amt = 240, cmt = "VERA_DEPOT", ii = 12, addl = 27, time = 0)
    return(e1)
  }
  if (scenario == "S3") {
    # verapamil + lithium 600 mg/d split BID
    e1 <- ev(amt = 240, cmt = "VERA_DEPOT", ii = 12, addl = 55, time = 0)
    e2 <- ev(amt = 300, cmt = "LI_DEPOT",   ii = 12, addl = 55, time = 0)
    return(c(e1, e2))
  }
  if (scenario == "S4") {
    # galcanezumab 300 mg SC monthly x 3
    e1 <- ev(amt = 300, cmt = "GALCA_SC", ii = 24*28, addl = 2, time = 0)
    return(e1)
  }
  if (scenario == "S5") {
    # prednisone 60 mg/d x 5 d -> taper -> verapamil 240 mg BID
    e1 <- ev(amt = 60, cmt = "PRED_DEPOT", ii = 24, addl = 4, time = 0)
    e2 <- ev(amt = 40, cmt = "PRED_DEPOT", ii = 24, addl = 2, time = 24*5)
    e3 <- ev(amt = 20, cmt = "PRED_DEPOT", ii = 24, addl = 2, time = 24*8)
    e4 <- ev(amt = 240, cmt = "VERA_DEPOT", ii = 12, addl = 50, time = 0)
    return(c(e1, e2, e3, e4))
  }
  if (scenario == "S6") {
    # GON block (one-shot effect 0.65) + verapamil
    e1 <- ev(time = 0, amt = 0.65, cmt = "GON_EFF", evid = 1)
    e2 <- ev(amt = 240, cmt = "VERA_DEPOT", ii = 12, addl = 27, time = 0)
    return(c(e1, e2))
  }
  ev
}

# ============================================================================
# Example simulation: run all six scenarios (12-week horizon)
# ============================================================================

simulate_all <- function(weeks = 12) {
  scenarios <- c("S0","S1","S2","S3","S4","S5","S6")
  out_list <- list()
  for (sc in scenarios) {
    ev_sc <- build_scenario(sc)
    out <- ch_mod %>%
      ev(ev_sc) %>%
      mrgsim(end = 24*7*weeks, delta = 1) %>%
      as_tibble() %>%
      mutate(scenario = sc)
    out_list[[sc]] <- out
  }
  bind_rows(out_list)
}

# ============================================================================
# Quick sanity helpers ---------------------------------------------------------
# (run with: source("ch_mrgsolve_model.R"); simulate_all() %>% ...)
# ============================================================================
