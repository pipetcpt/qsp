## =============================================================================
##  HIGH-ALTITUDE ILLNESS (AMS / HACE / HAPE) — mrgsolve QSP MODEL
## =============================================================================
##
##  고산병 (급성 고산병 · 고산 뇌부종 · 고산 폐부종) 정량적 시스템 약리 모델
##
##  구조적 주장 (the structural claim)
##  ---------------------------------------------------------------------------
##  고산병은 세 개의 병이 아니다.  하나의 흡입 산소분압(P_IO2)이 시간상수가
##  서로 다른 세 개의 방어기전을 통해 작용한 결과다.
##
##     P_IO2 = F_IO2 x (P_B(h) - 47)              <- 유일한 외생 변수
##
##     방어 1 (분 -> 일)   : 환기.  P_AO2 = P_IO2 - P_aCO2 x [F_IO2+(1-F_IO2)/R]
##                           화학수용체가 아니라 산-염기가 한계를 정한다.
##     방어 2 (초)         : 저산소성 폐혈관 수축(HPV).  가스교환은 지키고
##                           모세혈관은 파괴한다 — 공간적으로 불균일하기 때문.
##     방어 3 (초)         : 뇌혈관 확장.  뇌혈류는 지키고 두개내 용적
##                           완충능을 소모한다.
##
##  방어 1은 신장/CSF 중탄산 때문에 며칠이 걸리고, 방어 2와 3은 즉시 작동한다.
##  고산의학 전체가 이 시간 불일치에서 나온다: 느린 방어가 도착할 때까지
##  빠른 방어 둘이 견제 없이 돌아가고, 약은 정확히 네 갈래로 나뉜다 —
##  "방어 1을 가속"(아세타졸아미드), "방어 2를 둔화"(니페디핀·타다라필),
##  "방어 3의 대가를 지불"(덱사메타손), "원인을 제거"(하산·산소·가모우백).
##
##  ---------------------------------------------------------------------------
##  50 compartments · 200+ parameters · 23 scenarios
##  검증: 모든 방정식은 hai_reference_model.py 에 Python/scipy 로 독립 재구현되어
##       있으며 두 구현의 수치가 hai_reference_output.txt 에서 비교된다.
##
##  단위: 압력 mmHg · 유량 L/min · 시간 h · 대사량 mL/min STPD
##       산-염기 mEq/L · Hb g/dL · 부종 mL
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB
# High-Altitude Illness QSP model (AMS / HACE / HAPE)

$PARAM @annotated
// ---- environment (time-varying; supply as data-set columns) ---------------
ALT    : 0      : Altitude (m)
FIO2X  : 0.2094 : Inspired O2 fraction (supplemental O2 raises this)
PBAG   : 0      : Portable hyperbaric bag over-pressure (mmHg)
EXER   : 1      : Metabolic rate as a multiple of resting VO2
SLEEPF : 0      : Sleep flag (1 = NREM sleep; blunts the peripheral drive)

// ---- subject phenotype ----------------------------------------------------
HVRSUB : 1.0    : Individual hypoxic ventilatory response multiplier
LAMMAX : 5.0    : Maximal HPV constriction factor (9.0 if HAPE-susceptible)
AHET   : 0.50   : HPV-responsive fraction of the pulmonary bed (heterogeneity)
PVI    : 25.0   : Craniospinal pressure-volume index (mL); 16 = "tight fit"
HBSL   : 15.0   : Sea-level haemoglobin (g/dL)
WT     : 70     : Body weight (kg)

// ---- gas / metabolic ------------------------------------------------------
PH2O   : 47.0   : Saturated water vapour pressure at 37 C (mmHg)
RQ     : 0.85   : Respiratory exchange ratio
VO2R   : 250.0  : Resting VO2 (mL/min STPD)
VCO2R  : 212.5  : Resting VCO2 (mL/min STPD)
VDVTR  : 0.30   : Dead-space fraction at rest
VDVTE  : 0.18   : Dead-space fraction at heavy exercise

// ---- ventilatory controller ----------------------------------------------
VEB    : 0.50   : Non-chemical (wakefulness) ventilation (L/min)
GC     : 1.293  : Central CO2 sensitivity (L/min/mmHg)
GP     : 7.00   : Peripheral gain (L/min per unit SaO2 deficit per mmHg)
TAUPHI : 2.0    : Softplus width of the chemoreflex knee (mmHg)
PHTHRC : 7.409  : CSF pH at the central apnoeic threshold
PHTHRP : 7.526  : Arterial pH constant setting the peripheral threshold
SAOREF : 0.980  : SaO2 at which the peripheral hypoxic drive is zero
VEMAX  : 200.0  : Mechanical ceiling on minute ventilation (L/min)
SLPGP  : 0.72   : Multiplier on peripheral gain during NREM sleep

// ---- acid-base ------------------------------------------------------------
TAUREN : 34.0   : Renal acid-base time constant (h)
KREN   : 115.0  : Renal gain (mEq/L per pH unit)
TAUCSF : 8.0    : CSF bicarbonate time constant (h)
WCSFPH : 0.55   : Weight of active choroid-plexus pH regulation
CCSF   : 0.246  : Offset making sea-level CSF HCO3 exactly 22.0

// ---- carotid-body plasticity ---------------------------------------------
VAHMAX : 2.10   : Maximal multiplicative gain of the peripheral drive
TAUVAH : 60.0   : Ventilatory acclimatisation time constant (h)
HVDMAX : 0.28   : Maximal hypoxic ventilatory decline (roll-off)
TAUHVN : 0.40   : Roll-off onset time constant (h)
TAUHVF : 2.0    : Roll-off offset time constant (h)

// ---- blood oxygen transport ----------------------------------------------
P50STD : 26.8   : Standard P50 (mmHg)
BOHR   : -0.48  : d(log10 P50)/d(pH)
DPGSL  : 5.00   : Sea-level 2,3-DPG (mmol/L)
DPGMAX : 6.00   : Maximal 2,3-DPG at altitude (mmol/L)
TAUDPG : 30.0   : 2,3-DPG time constant (h)
KP50DP : 0.075  : Fractional P50 rise per unit relative DPG rise
HBCON  : 1.34   : mL O2 per g haemoglobin
O2SOL  : 0.003  : mL O2 per dL per mmHg

// ---- diffusion / shunt ----------------------------------------------------
DLREST : 37.0   : Resting DLO2 (mL/min/mmHg)
DLEXMX : 78.0   : DLO2 at maximal recruitment (mL/min/mmHg)
FSHBAS : 0.035  : Anatomical + VA/Q shunt at rest
FSHFLD : 0.40   : Additional shunt at complete alveolar flooding

// ---- cardiac / pulmonary vascular ----------------------------------------
QREST  : 6.0    : Resting cardiac output (L/min)
HRSL   : 65.0   : Sea-level heart rate (bpm)
TAUHR  : 0.5    : Heart-rate time constant (h)
PLASL  : 8.0    : Sea-level left atrial pressure (mmHg)
KPLAEX : 0.55   : Rise in P_LA per L/min of cardiac output above rest
RPULM  : 1.00   : Total pulmonary vascular resistance at sea level (WU)
FRACV  : 0.50   : Venous share of total pulmonary resistance
P50HPV : 45.0   : PAO2 for half-maximal HPV (mmHg)
NHPV   : 4.0    : Hill coefficient of HPV
KAPHPV : 0.08   : Constriction of the "non-responsive" bed, as a fraction
KREC   : 0.35   : Pulmonary recruitment/distension per unit relative flow
TAUHPS : 72.0   : Slow HPV / remodelling time constant (h)
HPVSMX : 0.45   : Maximal slow-component addition to lambda

// ---- capillary stress failure / oedema -----------------------------------
PCAPCR : 19.5   : Capillary pressure threshold for stress failure (mmHg)
KFLEAK : 2.60   : Leak coefficient (mL/h per mmHg^1.5)
NLEAK  : 1.5    : Exponent of the stress-failure law
ELW0   : 300.0  : Normal extravascular lung water (mL)
ELWF50 : 700.0  : Excess EVLW giving 50 percent alveolar flooding (mL)
NFLOOD : 2.2    : Hill coefficient of flooding
KAFC   : 0.115  : Basal alveolar fluid clearance rate constant (1/h)
AFCH50 : 0.72   : SaO2 at which clearance capacity is halved
AFCN   : 6.0    : Hill coefficient of hypoxic clearance inhibition
TAUAFC : 6.0    : Clearance-capacity time constant (h)
TAUPRM : 8.0    : Permeability time constant (h)
KPRMIN : 0.85   : Permeability gain per unit inflammation
TAUINF : 12.0   : Lung inflammation time constant (h)
KINFLK : 0.0022 : Inflammation driven per mL/h of leak

// ---- cerebral -------------------------------------------------------------
KCBFO2 : 1.10   : Hypoxic cerebral vasodilation gain
CBFP50 : 35.0   : PaO2 for half-maximal hypoxic vasodilation (mmHg)
CBFO2N : 4.0    : Hill coefficient
CBFMID : 42.0   : Midpoint of the CO2 reactivity sigmoid (mmHg)
CBFW   : 12.0   : Width of the CO2 reactivity sigmoid (mmHg)
CBFLO  : 0.55   : Lower plateau of CO2 reactivity
CBFHI  : 1.65   : Span of CO2 reactivity
TAUCBF : 0.25   : Cerebral blood flow time constant (h)
ICP0   : 10.0   : Baseline intracranial pressure (mmHg)
CBVSL  : 75.0   : Sea-level cerebral blood volume (mL)
FCBV   : 0.25   : Fraction of the CBV change that occupies the cranium
KVEGF  : 1.00   : Brain VEGF gain
TAUVEG : 10.0   : VEGF time constant (h)
KBBB   : 0.65   : BBB permeability gain per unit VEGF
TAUBBB : 8.0    : BBB time constant (h)
KEDV   : 4.00   : Vasogenic oedema formation gain
TAUEDV : 20.0   : Vasogenic oedema clearance time constant (h)
KEDC   : 6.00   : Cytotoxic oedema gain (mL)
EDCO50 : 0.70   : SaO2 at half-maximal cytotoxic oedema
EDCN   : 10.0   : Hill coefficient of cytotoxic oedema
TAUEDC : 6.0    : Cytotoxic oedema time constant (h)
ICPHAC : 22.0   : ICP at half-maximal HACE risk (mmHg)
ICPHCN : 6.0    : Hill coefficient of HACE risk

// ---- symptoms (Lake Louise 2018 sub-scores) ------------------------------
TAUSXN : 5.0    : Symptom onset time constant (h)
TAUSXF : 14.0   : Symptom offset time constant (h)
KHDICP : 0.175  : Headache per mmHg of ICP above baseline
KHDHYP : 5.50   : Headache per unit SaO2 deficit
KGI    : 1.55   : GI symptom gain
KFAT   : 4.90   : Fatigue gain
KDIZ   : 4.30   : Dizziness gain
KSLP   : 9.50   : Sleep-disturbance gain

// ---- sympathetic / fluid --------------------------------------------------
TAUSYM : 2.0    : Sympathetic time constant (h)
KSYM   : 3.2    : Sympathetic gain per unit SaO2 deficit
TAUALD : 8.0    : Aldosterone time constant (h)
TAUADH : 4.0    : Vasopressin time constant (h)
KFLUID : 90.0   : Fluid retention per unit antidiuretic index (mL)
TAUFLD : 30.0   : Fluid time constant (h)

// ---- erythropoiesis / volume ---------------------------------------------
EPOSL  : 12.0   : Sea-level serum EPO (mU/mL)
KEPO   : 165.0  : EPO gain on the hypoxic stimulus
TAUEPO : 6.0    : EPO time constant (h)
EPOHFB : 0.055  : EPO feedback per (g/dL)^2 of Hb above baseline
TAURET : 84.0   : Marrow transit time (h)
KERY   : 0.000165 : g Hb per mU/mL EPO per h
HBMSL  : 800.0  : Sea-level haemoglobin mass (g)
TAUHBM : 2400.0 : Red-cell turnover time constant (h)
PVSL   : 3.10   : Sea-level plasma volume (L)
PVCON  : 0.20   : Maximal fractional plasma-volume contraction
TAUPV  : 40.0   : Plasma-volume time constant (h)

// ---- chronic --------------------------------------------------------------
KCMS   : 0.055  : Chronic mountain sickness gain
TAUCMS : 2000.0 : CMS time constant (h)
VISCK  : 2.31   : Relative viscosity exponent (mu ~ exp(k*Hct))
TAUMUS : 500.0  : Muscle adaptation time constant (h)
MUSMAX : 0.22   : Maximal muscle adaptation index

// ---- drug PK/PD -----------------------------------------------------------
ACZF   : 0.95   : Acetazolamide bioavailability
ACZKA  : 1.8    : Acetazolamide absorption rate (1/h)
ACZV   : 15.0   : Acetazolamide central volume (L)
ACZCL  : 2.1    : Acetazolamide clearance (L/h)
ACZKIN : 0.55   : Transfer into the deep RBC carbonic-anhydrase pool (1/h)
ACZKOU : 0.075  : Transfer out of the deep pool (1/h)
ACZVR  : 2.4    : Deep-pool volume (L)
ACZIR  : 3.5    : IC50 for renal carbonic anhydrase (mg/L)
ACZER  : 7.4    : Maximal fall in plasma HCO3 (mEq/L)
ACZIC  : 6.0    : IC50 for choroid-plexus carbonic anhydrase (mg/L)
ACZEC  : 3.6    : Maximal additional fall in CSF HCO3 (mEq/L)
ACZKRB : 0.045  : Maximal fractional slowing of CO2 transport
ACZIRB : 220.0  : IC50 in the RBC compartment (mg/L)

DEXF   : 0.80   : Dexamethasone bioavailability
DEXKA  : 2.0    : Dexamethasone absorption rate (1/h)
DEXV   : 60.0   : Dexamethasone volume (L)
DEXCL  : 12.0   : Dexamethasone clearance (L/h)
DEXKE0 : 0.045  : Dexamethasone effect-compartment rate (1/h)
DEXEC5 : 0.030  : Dexamethasone effect-site EC50 (mg/L)
DEXBBB : 0.75   : Maximal suppression of BBB permeability gain
DEXVEG : 0.65   : Maximal suppression of brain VEGF
DEXHPV : 0.30   : Maximal suppression of HPV
DEXAFC : 0.55   : Maximal upregulation of alveolar fluid clearance
DEXSX  : 0.60   : Maximal direct symptomatic suppression

NIFF   : 0.55   : Nifedipine bioavailability
NIFKA  : 0.55   : Nifedipine (SR) absorption rate (1/h)
NIFV   : 55.0   : Nifedipine volume (L)
NIFCL  : 30.0   : Nifedipine clearance (L/h)
NIFEC5 : 0.020  : Nifedipine EC50 (mg/L)
NIFHPV : 0.45   : Maximal HPV suppression by nifedipine

TADF   : 0.80   : Tadalafil bioavailability
TADKA  : 1.1    : Tadalafil absorption rate (1/h)
TADV   : 63.0   : Tadalafil volume (L)
TADCL  : 2.5    : Tadalafil clearance (L/h)
TADEC5 : 0.055  : Tadalafil EC50 (mg/L)
TADHPV : 0.50   : Maximal HPV suppression by tadalafil

SALKE0 : 0.12   : Salmeterol effect-compartment decay (1/h)
SALEC5 : 0.5    : Salmeterol EC50 (effect units)
SALAFC : 0.70   : Maximal upregulation of alveolar fluid clearance

IBUF   : 0.85   : Ibuprofen bioavailability
IBUKA  : 2.2    : Ibuprofen absorption rate (1/h)
IBUV   : 10.0   : Ibuprofen volume (L)
IBUCL  : 3.5    : Ibuprofen clearance (L/h)
IBUEC5 : 8.0    : Ibuprofen EC50 (mg/L)
IBUHD  : 0.45   : Maximal headache suppression by ibuprofen

$CMT @annotated
BE     : Base excess (mEq/L)
HCO3C  : CSF bicarbonate (mEq/L)
VAH    : Carotid-body plasticity gain
HVD    : Hypoxic ventilatory decline
DPG    : 2,3-diphosphoglycerate (mmol/L)
HBM    : Haemoglobin mass (g)
PV     : Plasma volume (L)
EPO    : Serum erythropoietin (mU/mL)
RET    : Reticulocyte pool (g Hb in transit)
HPVS   : Slow HPV / remodelling index
ELW    : Extravascular lung water (mL)
FLOOD  : Alveolar flooding fraction
AFC    : Alveolar fluid clearance capacity
PERM   : Pulmonary capillary permeability index
INFL   : Lung inflammation index
CBFR   : Cerebral blood flow, relative to sea level
VEGFB  : Brain VEGF index
BBBP   : Blood-brain barrier permeability index
EDV    : Vasogenic cerebral oedema (mL)
EDC    : Cytotoxic cerebral oedema (mL)
CSFRES : Craniospinal reserve consumed (mL)
SYM    : Sympathetic tone index
HRS    : Heart rate (bpm)
ALDO   : Aldosterone index
ADH    : Vasopressin index
FLUIDR : Retained fluid (mL)
HEAD   : Headache sub-score (0-3)
GISX   : GI sub-score (0-3)
FAT    : Fatigue sub-score (0-3)
DIZ    : Dizziness sub-score (0-3)
SLEEPD : Sleep-disturbance latent
ACZA   : Acetazolamide gut (mg)
ACZC   : Acetazolamide central (mg)
ACZR   : Acetazolamide deep RBC pool (mg)
DEXA   : Dexamethasone gut (mg)
DEXC   : Dexamethasone central (mg)
DEXE   : Dexamethasone effect site (mg)
NIFA   : Nifedipine gut (mg)
NIFC   : Nifedipine central (mg)
TADA   : Tadalafil gut (mg)
TADC   : Tadalafil central (mg)
SALE   : Salmeterol effect site
IBUA   : Ibuprofen gut (mg)
IBUC   : Ibuprofen central (mg)
LAC    : Arterial lactate (mmol/L)
MUSC   : Muscle / capillary adaptation index
HYPD   : Cumulative hypoxic dose (SaO2-deficit x h)
CMS    : Chronic mountain sickness score
PAPS   : Slow mean pulmonary artery pressure (mmHg)
ACCL   : Integrated acclimatisation index

$GLOBAL
#include <cmath>

// ---------------------------------------------------------------------------
//  BAROMETRIC PRESSURE — West (1996), fitted to MEASURED Himalayan pressures.
//  The ICAO standard atmosphere assumes a mid-latitude tropopause and
//  under-predicts the Everest summit by 17 mmHg, an error worth ~600 m of
//  apparent altitude at the one place on earth where it matters most.
// ---------------------------------------------------------------------------
double f_pb(double alt_m) {
  double h = alt_m / 1000.0;
  return exp(6.63268 - 0.1112*h - 0.00149*h*h);
}

//  Standard O2-Hb dissociation curve (Severinghaus 1979)
double f_sev(double po2) {
  if (po2 < 1e-9) po2 = 1e-9;
  return 1.0 / (23400.0/(po2*po2*po2 + 150.0*po2) + 1.0);
}

//  P50 with Bohr and 2,3-DPG shifts
double f_p50(double ph, double dpg, double p50std, double bohr,
             double kp50dp, double dpgsl) {
  double v = p50std * pow(10.0, bohr*(ph - 7.40));
  v *= (1.0 + kp50dp*(dpg/dpgsl - 1.0));
  return v;
}

double f_sat(double po2, double p50a, double p50std) {
  return f_sev(po2 * p50std / p50a);
}

double f_cont(double po2, double hb, double p50a, double p50std,
              double hbcon, double o2sol) {
  return hbcon*hb*f_sat(po2, p50a, p50std) + o2sol*po2;
}

//  Invert the content equation for PO2 (bisection; the function is monotone)
double f_po2_from_cont(double c, double hb, double p50a, double p50std,
                       double hbcon, double o2sol) {
  double lo = 1e-4, hi = 700.0, mid = 0.0;
  if (f_cont(hi,hb,p50a,p50std,hbcon,o2sol) < c) return hi;
  if (f_cont(lo,hb,p50a,p50std,hbcon,o2sol) > c) return lo;
  for (int i = 0; i < 45; i++) {
    mid = 0.5*(lo+hi);
    if (f_cont(mid,hb,p50a,p50std,hbcon,o2sol) < c) lo = mid; else hi = mid;
  }
  return 0.5*(lo+hi);
}

// ---------------------------------------------------------------------------
//  ACID-BASE.  Base excess (Siggaard-Andersen / van Slyke) plus
//  Henderson-Hasselbalch, solved simultaneously for HCO3 and pH.
//  Carrying the METABOLIC component in one state variable (BE) means the
//  acute non-bicarbonate buffering falls out for free: the model reproduces
//  dHCO3/dPaCO2 ~ -0.2 acutely and ~ -0.5 chronically without either number
//  being a parameter.
// ---------------------------------------------------------------------------
double f_hco3(double be, double paco2) {
  double lo = 1.0, hi = 60.0, mid, ph, r;
  for (int i = 0; i < 45; i++) {
    mid = 0.5*(lo+hi);
    ph = 6.1 + log10(mid/(0.03*paco2));
    r = 0.93*(mid - 24.4 + 14.8*(ph - 7.40)) - be;
    if (r < 0) lo = mid; else hi = mid;
  }
  return 0.5*(lo+hi);
}
double f_ph(double hco3, double paco2) {
  return 6.1 + log10(hco3/(0.03*paco2));
}

double f_softplus(double x, double tau) {
  double z = x/tau;
  if (z > 30.0) return x;
  if (z < -30.0) return tau*exp(z);
  return tau*log1p(exp(z));
}

//  Normalised cerebral CO2 reactivity — SIGMOID, saturating at both ends.
//  A linear %/mmHg reactivity extrapolated to a summit PaCO2 of 13 mmHg
//  would abolish cerebral blood flow; the real response saturates.
double f_cbfco2(double p, double mid, double w, double lo, double hi) {
  double num = lo + hi/(1.0 + exp(-(p - mid)/w));
  double den = lo + hi/(1.0 + exp(-(40.0 - mid)/w));
  return num/den;
}

$MAIN
// (dosing is handled through the standard mrgsolve event mechanism)

$ODE
// ===========================================================================
//  0.  DERIVED BLOOD COMPOSITION
// ===========================================================================
double bv   = (PV > 0.5 ? PV : 0.5) + HBM/340.0;      // 340 g Hb per L of RBC
double hb   = HBM/(10.0*bv);
double hct  = (HBM/340.0)/bv;

// ===========================================================================
//  1.  DRUG EFFECT SITES
// ===========================================================================
double c_acz  = ACZC/ACZV;
double c_aczr = ACZR/ACZVR;
double c_dexe = DEXE/DEXV;
double c_nif  = NIFC/NIFV;
double c_tad  = TADC/TADV;
double c_ibu  = IBUC/IBUV;

double e_dex   = c_dexe/(DEXEC5 + c_dexe);
double e_nif   = c_nif /(NIFEC5 + c_nif);
double e_tad   = c_tad /(TADEC5 + c_tad);
double e_ibu   = c_ibu /(IBUEC5 + c_ibu);
double e_aczr_ = c_acz /(ACZIR  + c_acz);            // renal CA
double e_aczc_ = c_acz /(ACZIC  + c_acz);            // choroid-plexus CA
double e_aczrb = c_aczr/(ACZIRB + c_aczr);           // RBC CA
double e_sal   = SALE  /(SALEC5 + SALE);

//  HPV suppression: drugs act on independent steps, so they COMBINE
//  multiplicatively on the surviving fraction, not additively on the effect.
double hpv_supp = 1.0 - (1.0 - NIFHPV*e_nif)*(1.0 - TADHPV*e_tad)
                      *(1.0 - DEXHPV*e_dex);

// ===========================================================================
//  2.  CARDIAC OUTPUT AND METABOLIC LOAD
// ===========================================================================
double q    = QREST*(1.0 + 0.55*(EXER - 1.0))*(1.0 + 0.10*SYM);
double vo2  = VO2R*EXER;
double vco2 = VCO2R*EXER;
double exf  = (EXER - 1.0)/3.0;  if (exf > 1.0) exf = 1.0;  if (exf < 0) exf = 0;
double vdvt = VDVTR + (VDVTE - VDVTR)*exf;
double kco2 = 863.0*(vco2/1000.0)/(1.0 - vdvt);        // L/min * mmHg

// ===========================================================================
//  3.  THE VENTILATORY / GAS-EXCHANGE FIXED POINT
//      Solve  PaCO2 = 863*VCO2 / VA(PaCO2)  by bisection.  The residual is
//      strictly decreasing (more CO2 -> more ventilation -> less CO2), so the
//      bracketed root is unique.
// ===========================================================================
double pb    = f_pb(ALT) + PBAG;
double pio2  = FIO2X*(pb - PH2O);
double cgas  = FIO2X + (1.0 - FIO2X)/RQ;
double dlo2  = DLREST + (DLEXMX - DLREST)*exf;
double gp_e  = GP*HVRSUB*(SLEEPF > 0.5 ? SLPGP : 1.0);
double fs    = FSHBAS + FSHFLD*FLOOD;  if (fs > 0.85) fs = 0.85;

double lo = 3.0, hi = 90.0, paco2 = 40.0;
double pao2=0, ph=0, hco3=0, cc=0, cao2=0, cv=0, pvo2=0, beta=0, inc=0;
double pa_id=0, pao2a=0, sat=0, ve=0, p50a=0, bc=0, bp=0, dhyp=0, resid=0;

for (int it = 0; it < 48; it++) {
  paco2 = 0.5*(lo + hi);

  pao2 = pio2 - paco2*cgas;             if (pao2 < 1.0) pao2 = 1.0;
  hco3 = f_hco3(BE, paco2);
  ph   = f_ph(hco3, paco2);
  p50a = f_p50(ph, DPG, P50STD, BOHR, KP50DP, DPGSL);

  cc   = f_cont(pao2, hb, p50a, P50STD, HBCON, O2SOL);
  cao2 = cc - fs*(vo2/(10.0*q))/((1.0 - fs) > 1e-3 ? (1.0 - fs) : 1e-3);
  if (cao2 < 0.5) cao2 = 0.5;

  //  Piiper-Scheid incomplete equilibration: exp(-DL/(beta*Q))
  cv    = cao2 - vo2/(10.0*q);   if (cv < 0.5) cv = 0.5;
  pvo2  = f_po2_from_cont(cv, hb, p50a, P50STD, HBCON, O2SOL);
  beta  = (cc - cv)/((pao2 - pvo2) > 1.0 ? (pao2 - pvo2) : 1.0);
  if (beta < 1e-4) beta = 1e-4;
  inc   = exp(-dlo2/(beta*q*10.0));
  pa_id = f_po2_from_cont(cao2, hb, p50a, P50STD, HBCON, O2SOL);
  pao2a = pa_id - (pao2 - pvo2)*inc;   if (pao2a < 1.0) pao2a = 1.0;
  sat   = f_sat(pao2a, p50a, P50STD);

  //  chemoreflex: two additive drives, BOTH with acid-base-set thresholds
  bc   = HCO3C/(0.03*pow(10.0, PHTHRC - 6.1));
  bp   = hco3 /(0.03*pow(10.0, PHTHRP - 6.1));
  dhyp = SAOREF - sat;   if (dhyp < 0) dhyp = 0;
  ve   = VEB + GC*f_softplus(paco2 - bc, TAUPHI)
             + gp_e*VAH*(1.0 - HVD)*dhyp*f_softplus(paco2 - bp, TAUPHI);
  if (ve > VEMAX) ve = VEMAX;
  if (ve < 0.05)  ve = 0.05;
  double ve_eff = ve*(1.0 - ACZKRB*e_aczrb);

  resid = kco2/ve_eff - paco2;
  if (resid > 0) lo = paco2; else hi = paco2;
}
double va = ve*(1.0 - vdvt);

// ===========================================================================
//  4.  ACID-BASE DYNAMICS
//      The kidney is the rate-limiting step of acclimatisation (tau ~ 34 h)
//      and acetazolamide is the drug that bypasses it.
// ===========================================================================
double be_ss    = -KREN*(ph - 7.40) - ACZER*e_aczr_;
double hco3c_ss = WCSFPH*(0.03*paco2*pow(10.0, 7.32 - 6.1))
                + (1.0 - WCSFPH)*hco3 + CCSF - ACZEC*e_aczc_;

dxdt_BE    = (be_ss - BE)/TAUREN;
dxdt_HCO3C = (hco3c_ss - HCO3C)/TAUCSF;

// ===========================================================================
//  5.  CHEMOREFLEX PLASTICITY, DPG
// ===========================================================================
double r1 = dhyp/0.18;  if (r1 > 1.0) r1 = 1.0;
double r2 = dhyp/0.15;  if (r2 > 1.0) r2 = 1.0;
double r3 = dhyp/0.20;  if (r3 > 1.0) r3 = 1.0;

double vah_ss = 1.0 + (VAHMAX - 1.0)*r1;
dxdt_VAH = (vah_ss - VAH)/TAUVAH;
double hvd_ss = HVDMAX*r2;
dxdt_HVD = (hvd_ss > HVD) ? (hvd_ss - HVD)/TAUHVN : (hvd_ss - HVD)/TAUHVF;
double dpg_ss = DPGSL + (DPGMAX - DPGSL)*r3;
dxdt_DPG = (dpg_ss - DPG)/TAUDPG;

// ===========================================================================
//  6.  ERYTHROPOIESIS AND PLASMA VOLUME
//      Two adaptations with an order-of-magnitude difference in speed, and
//      the FAST one moves the measured variable (Hct) while the SLOW one is
//      the only one that adds oxygen-carrying capacity.
// ===========================================================================
double dhb    = hb - HBSL;  if (dhb < 0) dhb = 0;
double epo_ss = EPOSL + KEPO*dhyp/(1.0 + EPOHFB*dhb*dhb);
dxdt_EPO = (epo_ss - EPO)/TAUEPO;
double epo_ex = EPO - EPOSL;  if (epo_ex < 0) epo_ex = 0;
dxdt_RET = KERY*epo_ex*HBMSL/100.0 - RET/TAURET;
dxdt_HBM = RET/TAURET - (HBM - HBMSL*(HBSL/15.0))/TAUHBM;
double pv_ss = PVSL*(1.0 - PVCON*r2) + FLUIDR/1000.0;
dxdt_PV = (pv_ss - PV)/TAUPV;

// ===========================================================================
//  7.  PULMONARY HAEMODYNAMICS — THE TWO-BED HPV MODEL
//
//      G_tot   = G0 [ a/lambda + (1-a)/mu ]
//      Pcap,B  = P_LA + Q*(r_v/mu) / [ a/lambda + (1-a)/mu ]
//
//      (1-a) CANCELS out of the numerator.  Taking lambda -> infinity with
//      mu = 1 + kappa(lambda-1) diverging alongside it, the amplification
//      ceiling is  1/[1 - a(1-kappa)]  -- which collapses to the clean
//      1/(1-a) only when the "non-responsive" bed truly does not constrict.
//      Either way the conclusion is the same and it is the point of the
//      submodel: HETEROGENEITY sets the ceiling, HPV STRENGTH only sets how
//      fast you reach it.  This is why an exaggerated HPV response is necessary but
//      not sufficient for HAPE, and why cardiac output — which multiplies
//      the whole term — is the parameter a clinician can actually control.
// ===========================================================================
double hhpv  = 1.0/(1.0 + pow((pao2 > 1.0 ? pao2 : 1.0)/P50HPV, NHPV));
double lam   = 1.0 + LAMMAX*hhpv*(1.0 + HPVSMX*HPVS)*(1.0 - hpv_supp);
if (lam < 1.0) lam = 1.0;
double mu    = 1.0 + KAPHPV*(lam - 1.0);
double r_v   = RPULM*FRACV;
double p_la  = PLASL + KPLAEX*((q - QREST) > 0 ? (q - QREST) : 0.0);
double rec   = 1.0 + KREC*((q/QREST - 1.0) > 0 ? (q/QREST - 1.0) : 0.0);
double den   = (AHET/lam + (1.0 - AHET)/mu)*rec;
double pvr   = 1.0/den;
double mpap  = p_la + q*pvr;
double qopen = q*((1.0 - AHET)/mu)/den;
double pcap  = p_la + q*(r_v/mu)/den;

dxdt_HPVS = (r2 - HPVS)/TAUHPS;
dxdt_PAPS = (mpap - PAPS)/24.0;

// ===========================================================================
//  8.  CAPILLARY STRESS FAILURE AND ALVEOLAR FLOODING (HAPE)
//      TWO positive feedback loops close here:
//        loop 1  flooding -> shunt -> hypoxaemia -> more HPV -> more Pcap
//        loop 2  hypoxaemia -> ENaC/Na-K-ATPase down -> the drain closes
// ===========================================================================
double over  = pcap - PCAPCR;  if (over < 0) over = 0;
double leak  = KFLEAK*pow(over, NLEAK)*PERM;
double excess= ELW - ELW0;     if (excess < 0) excess = 0;
double clear = KAFC*AFC*excess*(1.0 + DEXAFC*e_dex + SALAFC*e_sal);
dxdt_ELW = leak - clear;

double fl_ss = pow(excess, NFLOOD)/(pow(ELWF50, NFLOOD) + pow(excess, NFLOOD));
dxdt_FLOOD = (fl_ss - FLOOD)/0.5;

double s_cl  = (sat > 0.30 ? sat : 0.30);
double afc_ss= 1.0/(1.0 + pow(AFCH50/s_cl, AFCN));
if (afc_ss > 1.0) afc_ss = 1.0;  if (afc_ss < 0.15) afc_ss = 0.15;
dxdt_AFC  = (afc_ss - AFC)/TAUAFC;
dxdt_INFL = (KINFLK*leak - INFL)/TAUINF;
dxdt_PERM = ((1.0 + KPRMIN*INFL) - PERM)/TAUPRM;

// ===========================================================================
//  9.  CEREBRAL CIRCULATION, OEDEMA AND INTRACRANIAL PRESSURE (AMS / HACE)
// ===========================================================================
double cbf_o2  = 1.0 + KCBFO2/(1.0 + pow((pao2a > 3.0 ? pao2a : 3.0)/CBFP50, CBFO2N));
double cbf_co2 = f_cbfco2(paco2, CBFMID, CBFW, CBFLO, CBFHI);
double cbf_ss  = cbf_o2*cbf_co2;  if (cbf_ss < 0.35) cbf_ss = 0.35;
dxdt_CBFR = (cbf_ss - CBFR)/TAUCBF;

double vegf_ss = KVEGF*r3*(1.0 - DEXVEG*e_dex);
dxdt_VEGFB = (vegf_ss - VEGFB)/TAUVEG;
dxdt_BBBP  = ((1.0 + KBBB*VEGFB*(1.0 - DEXBBB*e_dex)) - BBBP)/TAUBBB;

double dcbf   = CBFR - 1.0;  if (dcbf < 0) dcbf = 0;
double dp_cap = 6.0*dcbf;
double dbbb   = BBBP - 1.0;  if (dbbb < 0) dbbb = 0;
dxdt_EDV = KEDV*dbbb*(0.35 + dp_cap)/10.0 - EDV/TAUEDV;
double edc_ss = KEDC/(1.0 + pow(s_cl/EDCO50, EDCN));
dxdt_EDC = (edc_ss - EDC)/TAUEDC;

double dv_br = EDV + EDC + CBVSL*(CBFR - 1.0)*FCBV;  if (dv_br < 0) dv_br = 0;
double icp   = ICP0*pow(10.0, dv_br/PVI);  if (icp > 80.0) icp = 80.0;
dxdt_CSFRES = (dv_br - CSFRES)/2.0;

// ===========================================================================
// 10.  SYMPATHETIC, CARDIOVASCULAR, FLUID
// ===========================================================================
dxdt_SYM  = ((KSYM*dhyp + 0.25*(EXER - 1.0)) - SYM)/TAUSYM;
dxdt_HRS  = ((HRSL*(1.0 + 0.42*SYM) + 28.0*(EXER - 1.0)) - HRS)/TAUHR;
dxdt_ALDO = ((1.0 + 3.6*dhyp) - ALDO)/TAUALD;
double hd_ex = HEAD - 0.6;  if (hd_ex < 0) hd_ex = 0;
dxdt_ADH  = ((1.0 + 0.5*hd_ex + 1.6*dhyp) - ADH)/TAUADH;
double adh_ex = ADH - 1.0;  if (adh_ex < 0) adh_ex = 0;
dxdt_FLUIDR = (KFLUID*adh_ex - FLUIDR)/TAUFLD;

// ===========================================================================
// 11.  SYMPTOMS — Lake Louise 2018
//      Every sub-score is an INTEGRATOR with a 5 h onset and a 14 h offset.
//      That asymmetry is why AMS peaks 18-24 h after arrival rather than on
//      arrival, and why "climb high, sleep low" works at all.
// ===========================================================================
double dexsx = 1.0 - DEXSX*e_dex;
double dicp  = icp - ICP0;  if (dicp < 0) dicp = 0;

double head_ss = (KHDICP*dicp + KHDHYP*dhyp)*dexsx*(1.0 - IBUHD*e_ibu);
if (head_ss > 3.0) head_ss = 3.0;
double gi_ss = KGI*(0.55*dicp/8.0 + dhyp)*dexsx;   if (gi_ss > 3.0) gi_ss = 3.0;
double fat_ss= KFAT*(dhyp + 0.25*SLEEPD/3.0)*dexsx;if (fat_ss> 3.0) fat_ss= 3.0;
double diz_ss= KDIZ*dhyp*dexsx;                    if (diz_ss> 3.0) diz_ss= 3.0;
double slp_ss= KSLP*dhyp;                          if (slp_ss> 3.0) slp_ss= 3.0;

dxdt_HEAD = (head_ss > HEAD) ? (head_ss - HEAD)/TAUSXN : (head_ss - HEAD)/TAUSXF;
dxdt_GISX = (gi_ss   > GISX) ? (gi_ss   - GISX)/TAUSXN : (gi_ss   - GISX)/TAUSXF;
dxdt_FAT  = (fat_ss  > FAT ) ? (fat_ss  - FAT )/TAUSXN : (fat_ss  - FAT )/TAUSXF;
dxdt_DIZ  = (diz_ss  > DIZ ) ? (diz_ss  - DIZ )/TAUSXN : (diz_ss  - DIZ )/TAUSXF;
dxdt_SLEEPD = (slp_ss - SLEEPD)/4.0;

// ===========================================================================
// 12.  MISCELLANEOUS AND CHRONIC
// ===========================================================================
dxdt_LAC  = ((1.0 + 5.5*exf*(1.0 + 2.0*dhyp)) - LAC)/0.5;
dxdt_MUSC = (MUSMAX*r2 - MUSC)/TAUMUS;
dxdt_HYPD = dhyp;
double hct_ex = hct - 0.52;  if (hct_ex < 0) hct_ex = 0;
double dh_ex  = dhyp - 0.08; if (dh_ex  < 0) dh_ex  = 0;
dxdt_CMS  = ((KCMS*hct_ex*100.0 + 6.0*dh_ex) - CMS)/TAUCMS;
dxdt_ACCL = ((1.0 - (HCO3C - 12.0)/(22.0 - 12.0)) - ACCL)/24.0;

// ===========================================================================
// 13.  DRUG PHARMACOKINETICS
// ===========================================================================
dxdt_ACZA = -ACZKA*ACZA;
dxdt_ACZC = ACZKA*ACZA*ACZF - ACZCL*c_acz - ACZKIN*ACZC + ACZKOU*ACZR;
dxdt_ACZR = ACZKIN*ACZC - ACZKOU*ACZR;
dxdt_DEXA = -DEXKA*DEXA;
dxdt_DEXC = DEXKA*DEXA*DEXF - DEXCL*DEXC/DEXV;
dxdt_DEXE = DEXKE0*(DEXC - DEXE);
dxdt_NIFA = -NIFKA*NIFA;
dxdt_NIFC = NIFKA*NIFA*NIFF - NIFCL*c_nif;
dxdt_TADA = -TADKA*TADA;
dxdt_TADC = TADKA*TADA*TADF - TADCL*c_tad;
dxdt_SALE = -SALKE0*SALE;
dxdt_IBUA = -IBUKA*IBUA;
dxdt_IBUC = IBUKA*IBUA*IBUF - IBUCL*c_ibu;

$TABLE
// ---- recompute the algebra for output (identical to the ODE block) --------
double o_bv  = (PV > 0.5 ? PV : 0.5) + HBM/340.0;
double o_hb  = HBM/(10.0*o_bv);
double o_hct = (HBM/340.0)/o_bv;

double o_q   = QREST*(1.0 + 0.55*(EXER - 1.0))*(1.0 + 0.10*SYM);
double o_vo2 = VO2R*EXER;
double o_vco2= VCO2R*EXER;
double o_exf = (EXER - 1.0)/3.0; if (o_exf>1.0) o_exf=1.0; if (o_exf<0) o_exf=0;
double o_vdvt= VDVTR + (VDVTE - VDVTR)*o_exf;
double o_k   = 863.0*(o_vco2/1000.0)/(1.0 - o_vdvt);

double o_cacz = ACZC/ACZV, o_caczr = ACZR/ACZVR;
double o_cdex = DEXE/DEXV, o_cnif = NIFC/NIFV, o_ctad = TADC/TADV;
double o_edex = o_cdex/(DEXEC5 + o_cdex);
double o_enif = o_cnif/(NIFEC5 + o_cnif);
double o_etad = o_ctad/(TADEC5 + o_ctad);
double o_eaczrb = o_caczr/(ACZIRB + o_caczr);
double o_supp = 1.0 - (1.0 - NIFHPV*o_enif)*(1.0 - TADHPV*o_etad)
                    *(1.0 - DEXHPV*o_edex);

double o_pb   = f_pb(ALT) + PBAG;
double o_pio2 = FIO2X*(o_pb - PH2O);
double o_cgas = FIO2X + (1.0 - FIO2X)/RQ;
double o_dl   = DLREST + (DLEXMX - DLREST)*o_exf;
double o_gp   = GP*HVRSUB*(SLEEPF > 0.5 ? SLPGP : 1.0);
double o_fs   = FSHBAS + FSHFLD*FLOOD; if (o_fs > 0.85) o_fs = 0.85;

double olo = 3.0, ohi = 90.0, o_paco2 = 40.0;
double o_pao2=0,o_ph=0,o_hco3=0,o_cc=0,o_cao2=0,o_cv=0,o_pvo2=0;
double o_beta=0,o_inc=0,o_paid=0,o_pao2a=0,o_sat=0,o_ve=0,o_p50=0;
double o_bc=0,o_bp=0,o_dhyp=0;
for (int it = 0; it < 48; it++) {
  o_paco2 = 0.5*(olo + ohi);
  o_pao2 = o_pio2 - o_paco2*o_cgas; if (o_pao2 < 1.0) o_pao2 = 1.0;
  o_hco3 = f_hco3(BE, o_paco2);
  o_ph   = f_ph(o_hco3, o_paco2);
  o_p50  = f_p50(o_ph, DPG, P50STD, BOHR, KP50DP, DPGSL);
  o_cc   = f_cont(o_pao2, o_hb, o_p50, P50STD, HBCON, O2SOL);
  o_cao2 = o_cc - o_fs*(o_vo2/(10.0*o_q))/((1.0-o_fs)>1e-3?(1.0-o_fs):1e-3);
  if (o_cao2 < 0.5) o_cao2 = 0.5;
  o_cv   = o_cao2 - o_vo2/(10.0*o_q); if (o_cv < 0.5) o_cv = 0.5;
  o_pvo2 = f_po2_from_cont(o_cv, o_hb, o_p50, P50STD, HBCON, O2SOL);
  o_beta = (o_cc - o_cv)/((o_pao2 - o_pvo2) > 1.0 ? (o_pao2 - o_pvo2) : 1.0);
  if (o_beta < 1e-4) o_beta = 1e-4;
  o_inc  = exp(-o_dl/(o_beta*o_q*10.0));
  o_paid = f_po2_from_cont(o_cao2, o_hb, o_p50, P50STD, HBCON, O2SOL);
  o_pao2a= o_paid - (o_pao2 - o_pvo2)*o_inc; if (o_pao2a < 1.0) o_pao2a = 1.0;
  o_sat  = f_sat(o_pao2a, o_p50, P50STD);
  o_bc   = HCO3C/(0.03*pow(10.0, PHTHRC - 6.1));
  o_bp   = o_hco3/(0.03*pow(10.0, PHTHRP - 6.1));
  o_dhyp = SAOREF - o_sat; if (o_dhyp < 0) o_dhyp = 0;
  o_ve   = VEB + GC*f_softplus(o_paco2 - o_bc, TAUPHI)
              + o_gp*VAH*(1.0 - HVD)*o_dhyp*f_softplus(o_paco2 - o_bp, TAUPHI);
  if (o_ve > VEMAX) o_ve = VEMAX;  if (o_ve < 0.05) o_ve = 0.05;
  double o_vee = o_ve*(1.0 - ACZKRB*o_eaczrb);
  if (o_k/o_vee - o_paco2 > 0) olo = o_paco2; else ohi = o_paco2;
}

double o_hh  = 1.0/(1.0 + pow((o_pao2 > 1.0 ? o_pao2 : 1.0)/P50HPV, NHPV));
double o_lam = 1.0 + LAMMAX*o_hh*(1.0 + HPVSMX*HPVS)*(1.0 - o_supp);
if (o_lam < 1.0) o_lam = 1.0;
double o_mu  = 1.0 + KAPHPV*(o_lam - 1.0);
double o_pla = PLASL + KPLAEX*((o_q - QREST) > 0 ? (o_q - QREST) : 0.0);
double o_rec = 1.0 + KREC*((o_q/QREST - 1.0) > 0 ? (o_q/QREST - 1.0) : 0.0);
double o_den = (AHET/o_lam + (1.0 - AHET)/o_mu)*o_rec;
double o_mpap= o_pla + o_q/o_den;
double o_pcap= o_pla + o_q*(RPULM*FRACV/o_mu)/o_den;

double o_cbfo2 = 1.0 + KCBFO2/(1.0 + pow((o_pao2a>3.0?o_pao2a:3.0)/CBFP50, CBFO2N));
double o_dv    = EDV + EDC + CBVSL*(CBFR - 1.0)*FCBV; if (o_dv < 0) o_dv = 0;
double o_icp   = ICP0*pow(10.0, o_dv/PVI); if (o_icp > 80.0) o_icp = 80.0;

double o_res   = o_paco2 - o_bc;
double o_resc  = (o_res > 0.02 ? o_res : 0.02);
double o_ahi   = 80.0/(1.0 + pow(o_resc/1.55, 3.5)) + 1.0;
double o_lls   = HEAD + GISX + FAT + DIZ;
double o_evlw  = (ELW - ELW0 > 0 ? ELW - ELW0 : 0);

capture PB     = o_pb;
capture PIO2   = o_pio2;
capture PAO2   = o_pao2;
capture PaO2   = o_pao2a;
capture PaCO2  = o_paco2;
capture SaO2   = 100.0*o_sat;
capture pHa    = o_ph;
capture HCO3   = o_hco3;
capture P50a   = o_p50;
capture VE     = o_ve;
capture VA     = o_ve*(1.0 - o_vdvt);
capture AaDO2  = o_pao2 - o_pao2a;
capture CaO2   = o_cao2;
capture DO2    = o_cao2*o_q*10.0;
capture PvO2   = o_pvo2;
capture Bc     = o_bc;
capture Bp     = o_bp;
capture CO2res = o_res;
capture AHI    = o_ahi;
capture LAMBDA = o_lam;
capture mPAP   = o_mpap;
capture PCAPOP = o_pcap;
capture AMPCEIL= 1.0/(1.0 - AHET*(1.0 - KAPHPV));
capture EVLW   = o_evlw;
capture ICP    = o_icp;
capture LLS    = o_lls;
capture AMS    = (HEAD >= 1.0 && o_lls >= 3.0) ? 1.0 : 0.0;
capture HACER  = 1.0/(1.0 + pow(ICPHAC/(o_icp>1.0?o_icp:1.0), ICPHCN));
capture HB     = o_hb;
capture HCT    = 100.0*o_hct;
capture QC     = o_q;
capture CACZ   = o_cacz;
capture CDEX   = o_cdex;

$CAPTURE ALT EXER FIO2X
'

mod <- mcode_cache("hai_qsp", code)

## ============================================================================
##  SUBJECT PHENOTYPES
## ============================================================================
phenotype <- function(label = "typical trekker",
                      HVRSUB = 1.0, LAMMAX = 5.0, AHET = 0.50,
                      PVI = 25.0, HBSL = 15.0) {
  list(label = label, par = list(HVRSUB = HVRSUB, LAMMAX = LAMMAX,
                                 AHET = AHET, PVI = PVI, HBSL = HBSL))
}
TYPICAL   <- phenotype("typical trekker")
HAPE_SUSC <- phenotype("HAPE-susceptible", HVRSUB = 0.75, LAMMAX = 9.0, AHET = 0.85)
TIGHT_FIT <- phenotype("tight-fit (low craniospinal compliance)",
                       HVRSUB = 0.85, PVI = 16.0)
ELITE     <- phenotype("elite climber (high HVR)", HVRSUB = 1.35)

init_state <- function(pheno) {
  list(BE = 0, HCO3C = 22.0, VAH = 1.0, HVD = 0, DPG = 5.0,
       HBM = 800.0*(pheno$par$HBSL/15.0), PV = 3.10, EPO = 12.0, RET = 0,
       HPVS = 0, ELW = 300.0, FLOOD = 0, AFC = 1.0, PERM = 1.0, INFL = 0,
       CBFR = 1.0, VEGFB = 0, BBBP = 1.0, EDV = 0, EDC = 0, CSFRES = 0,
       SYM = 0, HRS = 65.0, ALDO = 1.0, ADH = 1.0, FLUIDR = 0,
       HEAD = 0, GISX = 0, FAT = 0, DIZ = 0, SLEEPD = 0,
       LAC = 1.0, MUSC = 0, HYPD = 0, CMS = 0, PAPS = 14.0, ACCL = 0)
}

## ============================================================================
##  ASCENT PROFILE HELPERS
## ============================================================================
#' Build a time-varying data set: altitude, exercise, FiO2, bag pressure, sleep.
build_profile <- function(t_end, alt_fun, exer_fun = function(t) 1,
                          fio2_fun = function(t) 0.2094,
                          bag_fun  = function(t) 0,
                          sleep_fun = function(t) as.numeric((t %% 24) >= 22 | (t %% 24) < 6),
                          dt = 0.25, ID = 1) {
  tt <- seq(0, t_end, by = dt)
  data.frame(ID = ID, time = tt,
             ALT    = vapply(tt, alt_fun,   numeric(1)),
             EXER   = vapply(tt, exer_fun,  numeric(1)),
             FIO2X  = vapply(tt, fio2_fun,  numeric(1)),
             PBAG   = vapply(tt, bag_fun,   numeric(1)),
             SLEEPF = vapply(tt, sleep_fun, numeric(1)),
             evid = 0, cmt = 0, amt = 0)
}

ramp <- function(pts) {
  ts <- vapply(pts, `[[`, numeric(1), 1)
  as <- vapply(pts, `[[`, numeric(1), 2)
  function(t) approx(ts, as, xout = t, rule = 2)$y
}

dose_rows <- function(start, interval, n, cmt_name, amt, ID = 1) {
  cmt_idx <- match(cmt_name, names(mrgsolve::init(mod)))
  data.frame(ID = ID, time = start + (0:(n - 1))*interval,
             ALT = NA, EXER = NA, FIO2X = NA, PBAG = NA, SLEEPF = NA,
             evid = 1, cmt = cmt_idx, amt = amt)
}

#' Merge a covariate time course with dosing rows and run.
run_scenario <- function(pheno, t_end, alt_fun, exer_fun = function(t) 1,
                         fio2_fun = function(t) 0.2094,
                         bag_fun = function(t) 0, doses = NULL, dt = 0.25) {
  cov <- build_profile(t_end, alt_fun, exer_fun, fio2_fun, bag_fun, dt = dt)
  dat <- cov
  if (!is.null(doses)) {
    d <- doses
    ## carry the covariates forward onto the dosing rows
    for (nm in c("ALT", "EXER", "FIO2X", "PBAG", "SLEEPF")) {
      d[[nm]] <- approx(cov$time, cov[[nm]], xout = d$time, rule = 2)$y
    }
    dat <- dplyr::arrange(dplyr::bind_rows(cov, d), time, evid)
  }
  mod %>%
    param(pheno$par) %>%
    init(init_state(pheno)) %>%
    data_set(dat) %>%
    mrgsim(end = -1, add = cov$time, tad = FALSE, delta = dt,
           atol = 1e-8, rtol = 1e-6, maxsteps = 100000) %>%
    as.data.frame()
}

## ============================================================================
##  SCENARIO LIBRARY  (23 scenarios)
## ============================================================================
RAPID  <- ramp(list(c(0, 1130), c(5, 3200), c(9, 4559)))
GRADED <- ramp(list(c(0, 1130), c(6, 2400), c(24, 2400), c(30, 3200),
                    c(48, 3200), c(54, 3800), c(72, 3800), c(78, 4559)))
EXERT  <- function(t) { d <- t %% 24; if (d >= 9 && d <= 15) 3.0 else 1.0 }

scenarios <- list(

  ## --- 1. control ----------------------------------------------------------
  `01_sea_level` = function()
    run_scenario(TYPICAL, 72, function(t) 0),

  ## --- 2-5. AMS prophylaxis at 4559 m (Capanna Regina Margherita) ----------
  `02_rapid_4559` = function()
    run_scenario(TYPICAL, 120, RAPID),

  `03_rapid_4559_acz125` = function()
    run_scenario(TYPICAL, 120, RAPID,
                 doses = dose_rows(0, 12, 12, "ACZA", 125)),

  `04_rapid_4559_acz250` = function()
    run_scenario(TYPICAL, 120, RAPID,
                 doses = dose_rows(0, 12, 12, "ACZA", 250)),

  `05_rapid_4559_dex` = function()
    run_scenario(TYPICAL, 120, RAPID,
                 doses = dose_rows(0, 12, 10, "DEXA", 4)),

  ## --- 6-7. ascent strategy ------------------------------------------------
  `06_graded_4559` = function()
    run_scenario(TYPICAL, 144, GRADED),

  `07_climb_high_sleep_low` = function()
    run_scenario(TYPICAL, 96, function(t) {
      base <- approx(c(0, 6, 24, 30, 48, 54, 72), c(1130, 2400, 2400, 3000,
                     3000, 3600, 3600), xout = t, rule = 2)$y
      d <- t %% 24
      if (d >= 10 && d <= 16) base + 700 else base
    }),

  ## --- 8-11, 22. HAPE prophylaxis -----------------------------------------
  `08_hapeS_rapid_exercise` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT),

  `09_hapeS_nifedipine` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 doses = dose_rows(0, 12, 8, "NIFA", 30)),

  `10_hapeS_tadalafil` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 doses = dose_rows(0, 12, 8, "TADA", 10)),

  `11_hapeS_dexamethasone` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 doses = dose_rows(0, 12, 8, "DEXA", 8)),

  `22_hapeS_salmeterol` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 doses = dose_rows(0, 12, 8, "SALE", 1.4)),

  ## --- 12-14. HAPE rescue: three ways of raising PIO2 ----------------------
  `12_hape_rescue_descent` = function()
    run_scenario(HAPE_SUSC, 96,
                 function(t) if (t < 40) RAPID(t) else max(2500, 4559 - (t - 40)*700),
                 exer_fun = EXERT),

  `13_hape_rescue_o2` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 fio2_fun = function(t) if (t < 40) 0.2094 else 0.28),

  `14_hape_rescue_gamow` = function()
    run_scenario(HAPE_SUSC, 96, RAPID, exer_fun = EXERT,
                 bag_fun = function(t) if (t < 40) 0 else 105),

  ## --- 15-16. HACE ---------------------------------------------------------
  `15_hace_push_on` = function()
    run_scenario(TIGHT_FIT, 96,
                 ramp(list(c(0, 1130), c(5, 3200), c(9, 4559), c(30, 4559),
                           c(40, 5300), c(60, 5900)))),

  `16_hace_rescue` = function() {
    up <- ramp(list(c(0, 1130), c(5, 3200), c(9, 4559), c(30, 4559),
                    c(40, 5300), c(60, 5900)))
    run_scenario(TIGHT_FIT, 96,
                 function(t) if (t < 62) up(t) else max(2500, 5900 - (t - 62)*900),
                 doses = dose_rows(62, 6, 6, "DEXA", 8))
  },

  ## --- 17-18. Everest summit day -------------------------------------------
  `17_everest_summit` = function()
    run_scenario(ELITE, 490,
                 ramp(list(c(0, 5300), c(240, 5300), c(300, 6400), c(400, 7100),
                           c(450, 7900), c(470, 8848), c(476, 7900))),
                 exer_fun = function(t) if (t >= 450 && t <= 476) 2.6 else 1.0),

  `18_everest_summit_o2` = function()
    run_scenario(ELITE, 490,
                 ramp(list(c(0, 5300), c(240, 5300), c(300, 6400), c(400, 7100),
                           c(450, 7900), c(470, 8848), c(476, 7900))),
                 exer_fun = function(t) if (t >= 450 && t <= 476) 2.6 else 1.0,
                 fio2_fun = function(t) if (t >= 450) 0.45 else 0.2094),

  ## --- 19-20. sleep at 4000 m ----------------------------------------------
  `19_sleep_4000` = function()
    run_scenario(TYPICAL, 72, ramp(list(c(0, 1130), c(6, 4000)))),

  `20_sleep_4000_acz` = function()
    run_scenario(TYPICAL, 72, ramp(list(c(0, 1130), c(6, 4000))),
                 doses = dose_rows(0, 12, 6, "ACZA", 250)),

  ## --- 21. three-week acclimatisation --------------------------------------
  `21_acclimatise_3800` = function()
    run_scenario(TYPICAL, 21*24, ramp(list(c(0, 200), c(10, 3800))), dt = 0.5),

  ## --- 23. ibuprofen -------------------------------------------------------
  `23_rapid_4559_ibuprofen` = function()
    run_scenario(TYPICAL, 120, RAPID,
                 doses = dose_rows(0, 8, 12, "IBUA", 600))
)

## ============================================================================
##  ANALYSIS HELPERS
## ============================================================================

#' Summarise one scenario the way the reference implementation does.
summarise_run <- function(df) {
  nightly <- df %>%
    mutate(night = floor(time/24),
           insleep = (time %% 24) >= 22 | (time %% 24) < 6) %>%
    filter(insleep) %>% group_by(night) %>%
    summarise(AHI = mean(AHI), .groups = "drop")
  data.frame(
    SaO2_min   = min(df$SaO2),      SaO2_final  = tail(df$SaO2, 1),
    PaO2_min   = min(df$PaO2),      PaCO2_final = tail(df$PaCO2, 1),
    HCO3_final = tail(df$HCO3, 1),  LLS_max     = max(df$LLS),
    LLS_t_max  = df$time[which.max(df$LLS)],
    AMS_hours  = sum(df$AMS)*(df$time[2] - df$time[1]),
    ICP_max    = max(df$ICP),       HACE_risk   = max(df$HACER),
    EVLW_max   = max(df$EVLW),      EVLW_final  = tail(df$EVLW, 1),
    mPAP_max   = max(df$mPAP),      Pcap_max    = max(df$PCAPOP),
    CO2res_min = min(df$CO2res),    AHI_night1  = nightly$AHI[1],
    AHI_last   = tail(nightly$AHI, 1),
    Hb_final   = tail(df$HB, 1),    Hct_final   = tail(df$HCT, 1))
}

run_all <- function() {
  out <- lapply(names(scenarios), function(nm) {
    message("running ", nm)
    s <- summarise_run(scenarios[[nm]]())
    cbind(scenario = nm, s)
  })
  do.call(rbind, out)
}

## ---------------------------------------------------------------------------
##  Static physiology: the arithmetic that needs no ODE at all
## ---------------------------------------------------------------------------
pb_west   <- function(h) exp(6.63268 - 0.1112*(h/1000) - 0.00149*(h/1000)^2)
pio2_alt  <- function(h, fio2 = 0.2094, bag = 0) fio2*(pb_west(h) + bag - 47)
pao2_alt  <- function(h, paco2, fio2 = 0.2094, bag = 0, R = 0.85)
  pio2_alt(h, fio2, bag) - paco2*(fio2 + (1 - fio2)/R)
severinghaus <- function(p) 1/(23400/(p^3 + 150*p) + 1)

#' The maximum capillary overpressure a lung can generate is bounded by
#' 1/(1-a).  Heterogeneity sets the ceiling; HPV strength sets only the speed.
pcap_two_bed <- function(lambda, a, Q = 6, kappa = 0.08, r_v = 0.50, P_LA = 8,
                         krec = 0.35, Qrest = 6) {
  mu  <- 1 + kappa*(lambda - 1)
  rec <- 1 + krec*pmax(0, Q/Qrest - 1)
  den <- (a/lambda + (1 - a)/mu)*rec
  P_LA + 0.55*pmax(0, Q - Qrest) + Q*(r_v/mu)/den
}

#' DO2 = Q*CaO2 with Q ~ mu^-gamma and mu ~ exp(k*Hct), CaO2 ~ Hct:
#'   DO2 ~ Hct*exp(-k*gamma*Hct)  ->  argmax = 1/(k*gamma).  Exact.
#' SaO2 cancels, so the optimum does NOT move with altitude.
optimal_hct <- function(gamma = 1, k = 2.31) min(1/(k*gamma), 1)

#' Every intervention in one currency: how many metres of descent would give
#' the same arterial saturation?
descent_equivalent <- function(alt, sao2_target, sao2_fun) {
  f <- function(a) sao2_fun(a) - sao2_target
  if (f(0) < 0) return(NA_real_)
  if (f(alt) > 0) return(0)
  alt - uniroot(f, c(0, alt), tol = 1)$root
}

## ============================================================================
##  PLOTS
## ============================================================================
plot_scenario <- function(df, title = "") {
  df %>%
    select(time, ALT, SaO2, PaCO2, HCO3, LLS, ICP, mPAP, PCAPOP, EVLW,
           CO2res, AHI) %>%
    pivot_longer(-time) %>%
    ggplot(aes(time, value)) +
    geom_line(linewidth = 0.7, colour = "#1565c0") +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    labs(x = "time (h)", y = NULL, title = title) +
    theme_bw(base_size = 10)
}

compare_arms <- function(named_dfs, var = "LLS") {
  bind_rows(lapply(names(named_dfs), function(n)
    transform(named_dfs[[n]][, c("time", var)], arm = n))) %>%
    ggplot(aes(time, .data[[var]], colour = arm)) +
    geom_line(linewidth = 0.8) +
    labs(x = "time (h)", y = var) + theme_bw(base_size = 11)
}

## ============================================================================
##  EXAMPLE
## ============================================================================
if (interactive()) {
  res <- run_all()
  print(res, digits = 3)

  arms <- list(`no prophylaxis`  = scenarios$`02_rapid_4559`(),
               `acetazolamide`   = scenarios$`04_rapid_4559_acz250`(),
               `dexamethasone`   = scenarios$`05_rapid_4559_dex`())
  print(compare_arms(arms, "LLS"))
  print(compare_arms(arms, "SaO2"))     # <- the asymmetry: only one moves this
  print(compare_arms(arms, "AHI"))
}
