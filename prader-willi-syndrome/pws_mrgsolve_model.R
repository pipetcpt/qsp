## =====================================================================
##  pws_mrgsolve_model.R
##  Prader-Willi Syndrome (PWS) — QSP model for mrgsolve
##  프래더-윌리 증후군 — mrgsolve QSP 모델
##
##  64 ODEs:  5 prohormone pools + 7 neuropeptide/gut signals
##          + 7 energy & growth + 3 food-seeking/behaviour
##          + 5 somatotropic + 4 glucose/insulin + 4 gonadal/bone
##          + 4 airway + 3 musculoskeletal + 20 drug PK + 2 diagnostics
##
##  THE ONE STRUCTURAL IDEA
##  -----------------------
##  The paternal 15q11-q13 lesion is written as ONE SCALAR — prohormone
##  convertase PC1/3 activity, PC13 — feeding five prohormone branches.
##  Each branch carries a single dimensionless ESCAPE RATIO
##
##      eps_i = d_i * Km_i / kcat_i     (precursor escape / processing)
##
##  and eps_i alone fixes both how much product that branch loses and
##  how much precursor it accumulates, in OPPOSITE directions:
##
##      L_i = 1 - (1+eps_i)/(1+eps_i/PC13)      product loss     [Eq. A]
##      R_i = (1+eps_i)/(PC13+eps_i)            precursor gain   [Eq. B]
##
##  At PC13 = 0.40 that gives, with no per-branch fitting at all:
##      pro-oxytocin (eps 6.0)   product -56%   precursor 1.09x
##      POMC         (eps 1.5)   product -47%   precursor 1.32x
##      pro-GHRH     (eps 1.0)   product -43%   precursor 1.43x
##      proinsulin   (eps 0.15)  product -16%   precursor 2.09x
##      proghrelin   (eps 0.10)  product -12%   precursor 2.20x
##
##  Two exact consequences follow, and they are a measurement theory:
##      precursor + product = S        -> blind to PC13 (conservation)
##      precursor / product = 1/PC13   -> identical for ALL branches,
##                                        blind to eps
##  So the two easy assays are each blind to exactly what the other
##  measures, and neither identifies the failing branch.
##
##  FOUR MORE STRUCTURAL CHOICES
##  ----------------------------
##  (1) alpha-MSH is NOT a parallel satiety arm.  MC4R satiety signalling
##      is relayed through PVN oxytocin neurons, so alpha-MSH enters as a
##      SATURATING INPUT GAIN on the oxytocin arm.  PWS blocks pro-oxytocin
##      processing, i.e. BELOW MC4R, so MC4R agonism is bounded above by
##      the very quantity PWS lacks.  Oxytocin-receptor agonism acts on
##      BOTH sides of the synapse (post-synaptic satiety arm AND
##      pre-synaptic AgRP inhibition) and is not so bounded.  That
##      asymmetry alone separates the carbetocin and setmelanotide
##      results, with no drug-specific parameter.
##  (2) Satiety is a HARMONIC MEAN, SAT = 1/sum(w_i/x_i).  The weakest arm
##      rate-limits: normalizing all four non-oxytocin arms moves SAT from
##      0.46 to only 0.49 of a possible 1.00.
##  (3) Food-seeking SEEK is BISTABLE (Hill-4 self-sensitization,
##      reinforced by food actually obtained).  At food-secure inputs it
##      has three fixed points; at free access the LOW one is annihilated
##      (saddle-node at CUE* = 0.58).  Food security does not lower the
##      drive — it restores the EXISTENCE of a low state.
##  (4) Growth hormone reaches the upper airway TWICE on different clocks:
##      IGF-1 -> lymphoid hypertrophy (tau ~20 d, with 100-d involution)
##      and lean mass -> respiratory muscle (tau 75 d) plus falling fat.
##      The transient AHI worsening at weeks 4-8 is the interference
##      between them, not a fitted curve.
##
##  TIME SCALE
##  ----------
##  DAY scale throughout.  Every species with a sub-hour plasma half-life
##  (ghrelin, insulin, GH, LH) is a day-scale pool whose fast kinetics are
##  pre-averaged; no rate constant exceeds ~14 /d.  Drug concentrations
##  are therefore APPARENT day-scale exposures calibrated to reproduce the
##  correct 24-h AUC, not the true Cmax.  FPOTGH < 1 carries the
##  pulsatility correction for somatropin: pulsatile GH is more
##  IGF-1-efficient per unit AUC than a flat exposure.
##
##  VALIDATION
##  ----------
##  No R runtime was available where this was written, so every equation
##  below was first implemented and executed in dependency-free Python
##  (pws_reference_model.py, fixed-step RK4) and this file is a
##  line-by-line port of what that produced.  The Python run exposed
##  several real defects and each is flagged [DEFECT n] at the point of
##  the fix:
##    [DEFECT 1] neuropeptide tones were not normalized to the control, so
##               the CONTROL had SAT = 0.05 and latched into hyperphagia.
##    [DEFECT 2] the Cunningham resting-energy intercept (370 kcal) is an
##               adult regression and starves a neonate; replaced by a
##               Schofield reference scaled by body composition.
##    [DEFECT 3] GH -> IGF-1 as a power law gave IGF-1 SDS +20 and a
##               height of 8 metres; replaced by a saturating hyperbola.
##    [DEFECT 4] lean mass was accrued BOTH from the energy partition and
##               from a height-growth term, double-counting it; lean is
##               now a developmental target and fat is the energy buffer.
##    [DEFECT 5] mean gastric content cannot distinguish meal-related
##               distension from chronic overfill on a day scale, which
##               turned PWS's delayed gastric emptying into a permanent
##               satiety bonus; the distension signal is now meal load.
##    [DEFECT 6] activity cost was charged on TOTAL fat, so a control
##               losing fat gained expenditure and drifted to 2% body fat;
##               only EXCESS fat over the age reference costs mobility.
##
##  Numerical convergence: halving dt from 0.125 to 0.0625 d over a
##  12-year run moves every reported endpoint by < 1.1e-3 relative.
##
##  ==  USAGE  ==========================================================
##      library(mrgsolve); library(dplyr); library(tidyr); library(ggplot2)
##      mod <- mread("pws_mrgsolve_model.R")
##      # then see the scenario driver at the foot of this file
##  =====================================================================

$PROB
# Prader-Willi Syndrome QSP model
- One convertase, five branches, a harmonic-mean satiety integrator,
  a bistable food-seeking state, and two clocks on the upper airway.
- 64 ODEs. Time in DAYS from birth. Run `pws_reference_model.py` for the
  verified reference output that every number in README.md quotes.

$GLOBAL
#define NAGE 25

// ---- reference (non-PWS) developmental scaffold -----------------------
static const double AGE_T[NAGE] = {0,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                                   15,16,17,18,20,25,30,40,60};
static const double LFM_T[NAGE] = {2.7,5.4,7.2,9.6,11.4,13.0,14.5,16.2,18.0,
                                   20.0,22.2,24.8,27.5,30.5,34.5,39.0,43.5,
                                   47.0,49.5,51.0,53.0,55.0,56.0,56.0,54.0};
static const double FM_T[NAGE]  = {0.5,2.3,2.9,3.2,3.3,3.4,3.6,4.0,4.5,5.2,
                                   6.0,6.8,7.6,8.2,8.5,8.8,9.2,9.8,10.5,11.2,
                                   12.5,14.0,15.5,18.0,20.0};
static const double HT_T[NAGE]  = {50,68,76,87,96,103,110,116,122,128,133,138,
                                   144,150,157,164,170,174,176,177,177,177,
                                   177,177,176};
static const double IGF_T[NAGE] = {50,55,60,68,76,85,95,110,124,140,160,180,
                                   215,260,310,350,348,340,310,280,250,210,
                                   190,160,120};
static const double GHN_T[NAGE] = {3.6,3.5,3.4,3.2,3.1,3.0,3.0,2.9,2.9,2.9,
                                   3.0,3.1,3.4,3.6,4.0,4.2,4.0,3.8,3.4,3.0,
                                   2.4,1.8,1.5,1.1,0.8};
static const double LY_T[NAGE]  = {0.55,0.70,0.85,1.05,1.20,1.28,1.30,1.26,
                                   1.18,1.10,1.02,0.96,0.90,0.86,0.82,0.79,
                                   0.77,0.75,0.74,0.73,0.72,0.70,0.70,0.70,
                                   0.70};
static const double CAP_T[NAGE] = {1.00,0.85,0.62,0.42,0.30,0.23,0.18,0.15,
                                   0.12,0.10,0.09,0.08,0.07,0.07,0.06,0.06,
                                   0.06,0.05,0.05,0.05,0.05,0.05,0.05,0.06,
                                   0.08};
static const double SUCK_T[NAGE]= {0.42,0.50,0.66,0.86,0.96,1.00,1.00,1.00,
                                   1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,
                                   1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,
                                   1.00};
static const double GVB_T[NAGE] = {36,20,13,10,8.0,7.2,6.5,6.2,6.0,5.6,5.0,
                                   4.6,4.0,3.4,2.6,1.8,1.0,0.45,0.15,0.05,
                                   0,0,0,0,0};
static const double GVS_T[NAGE] = {0,0,0,0,0,0,0,0,0,0.2,0.8,1.8,2.6,3.8,4.6,
                                   4.4,3.2,1.6,0.7,0.2,0,0,0,0,0};
static const double SXR_T[NAGE] = {0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,
                                   0.03,0.03,0.05,0.10,0.24,0.46,0.70,0.86,
                                   0.94,0.98,1.00,1.00,1.00,1.00,1.00,0.95,
                                   0.80};

static double itp(const double* tab, double a) {
  if (a <= AGE_T[0]) return tab[0];
  if (a >= AGE_T[NAGE-1]) return tab[NAGE-1];
  int lo = 0, hi = NAGE - 1;
  while (hi - lo > 1) { int m = (lo + hi) / 2;
    if (AGE_T[m] <= a) lo = m; else hi = m; }
  double f = (a - AGE_T[lo]) / (AGE_T[hi] - AGE_T[lo]);
  return tab[lo] + f * (tab[hi] - tab[lo]);
}

// Schofield resting energy expenditure by age band and body weight
static double schofield(double a, double bw) {
  if (a <  3.0) { double v = 60.9*bw - 54.0; return v < 90.0 ? 90.0 : v; }
  if (a < 10.0) return 22.7*bw + 495.0;
  if (a < 18.0) return 17.5*bw + 651.0;
  if (a < 30.0) return 15.3*bw + 679.0;
  return 11.6*bw + 879.0;
}

// height-age, for the weight-FOR-HEIGHT target a PWS clinic titrates to
static double height_age(double ht) {
  if (ht <= HT_T[0]) return 0.0;
  for (int i = 0; i < NAGE - 1; ++i) {
    if (HT_T[i] <= ht && ht <= HT_T[i+1]) {
      if (HT_T[i+1] == HT_T[i]) return AGE_T[i];
      double f = (ht - HT_T[i]) / (HT_T[i+1] - HT_T[i]);
      return AGE_T[i] + f * (AGE_T[i+1] - AGE_T[i]);
    }
  }
  return 18.0;
}

// quantities shared between $ODE, $TABLE and $CAPTURE
double AGEY, PC13, BW, FMx, LFMx;
double CGH, CLG, CCB, CDZ, CSM, CSG, COC, CLV, CTS, CMF;
double ECB, EDZx, EDZI, ESMx, ESGG, ESGS, ESGA, EOCg, EOCgh, EOCin, EMFx;
double GHREL, GHX, IGFT, IGFdrv, EXC, GHSEC;
double REE, ACTt, ATh, EIREQ, EIDRV, EICAP, EIx, TEEx, ESUFF, DITh;
double ASKx, FRUSTx, REINFx, SELFGx, TITRx, BWTGT;
double MLOAD, FASTs, FASTg, RI, RELAYx, MCIN;
double X1,X2,X3,X4,X5, SATx, DRVI, DRVE, GHRARM, AGEFF, GATEHP;
double SEXTOT, SEXREL, SEXFAC, EPIg, GVx, GVF, FMODx, LFMT;
double SPOMCx, SPOXTx, SPINSx, SPGHRx, PHIG;
double fPOMC, fPOXT, fPGHRH, fPINS, fPGHR, FACx;
double SIx, INSSEC, AHIT, AHIOBS, RMST, TONET, TONE0, LEPSAT, LEPF;

$PARAM @annotated
// ---- 1. genetic lesion layer ----------------------------------------
LES    :  1.0  : 1 = PWS (paternal 15q11-q13 silent), 0 = control (-)
DPC13  :  0.60 : fractional loss of hypothalamic/islet PC1/3 activity (-)
FOXTN  :  0.75 : surviving PVN oxytocin neuron fraction (-)
FPYYS  :  0.60 : blunted post-prandial PYY response (-)
FVAG   :  0.70 : reduced vagal afferent / distension gain (-)
FGHRC  :  2.20 : increased gastric X/A cell secretory drive (-)
FACYLP :  1.35 : raised acyl:unacyl ghrelin ratio in PWS (-)
FADPN  :  1.45 : relative hyperadiponectinaemia for a given fat mass (-)
FHYPO  :  0.22 : hypogonadotropic GnRH amplitude factor (-)
FSOM   :  0.75 : somatotroph secretory reserve factor (-)
DMKRN3 :  1.20 : MKRN3 loss advances the pubertal gate (y)
FTONE0 :  0.60 : neonatal hypotonia factor (-)
FTONEP :  0.78 : persistent (adult) muscle-tone deficit in PWS (-)
FSUB   :  1.00 : subtype modifier on behaviour (DEL 1.15 / mUPD 0.90) (-)

// ---- 2. prohormone processing: kc normalized to 1/d, so eps == d ----
KCPOMC :  1.0  : POMC processing rate constant (1/d)
EPSPOMC:  1.50 : POMC escape ratio (-)
KCPOXT :  1.0  : pro-oxytocin processing rate constant (1/d)
EPSPOXT:  6.00 : pro-oxytocin escape ratio (HIGH) (-)
KCPGHRH:  1.0  : pro-GHRH processing rate constant (1/d)
EPSPGHRH: 1.00 : pro-GHRH escape ratio (-)
KCPINS :  1.0  : proinsulin processing rate constant (1/d)
EPSPINS:  0.15 : proinsulin escape ratio (LOW) (-)
KCPGHR :  1.0  : proghrelin processing rate constant (1/d)
EPSPGHR:  0.10 : proghrelin escape ratio (LOW) (-)
SPOMC0 :  1.0  : POMC synthesis rate (-/d)
SPOXT0 :  1.0  : pro-oxytocin synthesis rate (-/d)
SPGHRH0:  1.0  : pro-GHRH synthesis rate (-/d)
SPINS0 :  1.0  : proinsulin synthesis rate (-/d)
SPGHR0 :  1.0  : proghrelin synthesis rate (-/d)
KLEPP  :  0.45 : leptin drive on POMC transcription (-)
KLEPR0 : 12.0  : CNS leptin signal half-saturation (ng/mL)

// ---- 3. neuropeptide tones (normalized so control = 1.0) ------------
KSAMSH :  6.25 : alpha-MSH formation gain = KDAMSH*(1+EPSPOMC) (-)
KDAMSH :  2.50 : alpha-MSH tone turnover (1/d)
KSOXT  : 49.0  : oxytocin formation gain = KDOXT*(1+EPSPOXT) (-)
KDOXT  :  7.00 : oxytocin tone turnover (1/d)
SAGRP  :  3.0  : AgRP/NPY synthesis (-/d)
KDAGRP :  1.0  : AgRP/NPY turnover (1/d)
KINHM  :  0.00 : alpha-MSH inhibition of AgRP - ZERO, action is post-synaptic (-)
KINHO  :  2.00 : oxytocin pre-synaptic inhibition of AgRP (-)
KLEPA  : 25.0  : leptin suppression of AgRP (ng/mL)
KFASTA :  0.40 : fasting amplification of AgRP (-)
X1CTL  :  1.000: control oxytocin-arm value, for relay integrity (-)
X1REF  :  0.285: reference PWS oxytocin-arm value, for behaviour term (-)

// ---- 4. ghrelin ------------------------------------------------------
KSAG   : 12400.0 : acyl-ghrelin formation gain (pg/mL/d per unit flux)
FACYL  :  0.22   : GOAT acylation fraction (-)
KCLAG  :  4.0    : acyl-ghrelin day-scale pool turnover (1/d)
KDEACYL:  0.55   : plasma de-acylation AG -> UAG (1/d)
KSUAG  : 10500.0 : unacylated ghrelin formation gain (pg/mL/d)
KCLUAG :  4.0    : unacylated ghrelin turnover (1/d)
KIINS  : 18.0    : insulin suppression of ghrelin secretion (uU/mL)
KFASTG :  0.45   : fasting amplification of ghrelin secretion (-)
KAG50  : 300.0   : GHS-R1a occupancy half-saturation - SATURATION (pg/mL)
KGHRD  :  0.10   : weight of the ghrelin arm in orexigenic drive (-)

// ---- 5. gut peptides / gastric handling -----------------------------
SPYY   : 18.0  : PYY formation gain (pmol/L/d)
KDPYY  :  2.0  : PYY turnover (1/d)
PYY0   :  9.0  : reference PYY (pmol/L)
SGLP1  : 26.0  : GLP-1 formation gain (pmol/L/d)
KDGLP1 :  2.6  : GLP-1 turnover (1/d)
GLP10  : 10.0  : reference GLP-1 (pmol/L)
KEMPT  :  5.0  : day-scale gastric emptying (1/d)
FGE    :  0.70 : PWS delayed gastric emptying (-)

// ---- 6. leptin -------------------------------------------------------
KSLEP  :  1.37 : leptin formation per kg fat (ng/mL/kg/d)
KDLEP  :  1.6  : leptin turnover (1/d)

// ---- 7. energy balance and body composition ------------------------
WREEL  :  0.88 : lean-mass weight in the REE composition term (-)
WREEF  :  0.12 : fat-mass weight in the REE composition term (-)
PALX   :  0.62 : activity thermogenesis as fraction of REE at ACT = 1 (-)
ACTREF :  1.000: reference activity level, defines EIREQ (-)
FDIT   :  0.10 : diet-induced thermogenesis (-)
KGHRMR :  0.16 : GH signal to resting energy expenditure (-)
ACT0   :  1.00 : activity scale (-)
TAUACT : 30.0  : activity adaptation time constant (d)
KFATACT:  0.45 : EXCESS fat penalty on activity (-)
KFA    : 22.0  : half-constant of the excess-fat activity penalty (kg)
RHOF   : 9440.0: energy density of fat mass (kcal/kg)
RHOL   : 1800.0: energy density of fat-free mass (kcal/kg)
FMMIN  :  1.2  : lower guard on fat mass (kg)
TAULFM : 180.0 : lean mass tracks its developmental target (d)
WLIGF0 :  0.78 : floor of the IGF-1 term in the lean target (-)
LIGFCAP:  1.25 : cap on the IGF-1 drive in the lean target (-)
WLTON0 :  0.80 : floor of the muscle-tone term in the lean target (-)
FMODCAP:  1.05 : ceiling on the lean-target modifier (-)
KSEXPART: 0.35 : gonadal-steroid contribution to the lean target (-)
EIBF   :  0.933: appetite set-point as fraction of the requirement (-)
KEIS   :  0.55 : realized food-seeking drive on intake (-)
KEIH   :  0.35 : hyperphagia score drive on intake (-)
KLEPEI :  0.25 : adiposity negative feedback on intake (-)
KTITR  :  1.20 : caregiver titration gain toward the weight target (-)
EICAPF :  1.40 : absolute physiological ceiling on EI/EIREQ (-)
AVAIL  :  1.00 : portions OFFERED, as a multiple of the requirement (-)
CUE    :  0.28 : external food-cue and access exposure (-)
TITR   :  1.00 : 0 = portions ignore the growth chart, 1 = titrated (-)
BWTGTR :  1.15 : weight-for-height target the family is holding (-)
FOSA   :  1.00 : upper-airway (anatomical) vulnerability multiplier (-)
GVFL   :  0.72 : GH-independent floor of growth velocity (-)
GVCAP  :  1.45 : ceiling on the IGF-1 drive to growth velocity (-)
KSEXH  : 380.0 : sex-steroid half-constant for the growth spurt (ng/dL)
BACLOSE: 17.2  : bone age at epiphyseal fusion (y)
KBASEX :  0.30 : sex-steroid acceleration of bone age (-)

// ---- 8. food-seeking bistability -----------------------------------
KSON   :  0.120: cue-driven activation of food-seeking (1/d)
KSELF  :  3.10 : self-sensitization amplitude (1/d)
KSHALF :  0.450: Hill-4 midpoint of self-sensitization (-)
NSELF  :  4.0  : Hill coefficient of self-sensitization (-)
KSOFF  :  2.20 : satiety-driven termination (1/d)
RFL    :  0.55 : floor of the reinforcement factor (-)
W1     :  0.44 : weight of the PVN-oxytocin satiety arm (-)
W2     :  0.16 : weight of the vagal/distension arm (-)
W3     :  0.13 : weight of the PYY arm (-)
W4     :  0.14 : weight of the GLP-1 arm (-)
W5     :  0.13 : weight of the leptin/insulin CNS arm (-)
KREL   :  0.20 : MC4R to PVN relay half-saturation - THE SERIES GAIN (-)
NAGRP  :  1.00 : AgRP exponent in the drive (-)
HQMAX  : 36.0  : HQ-CT maximum score (points)
WHQS   :  0.35 : weight of realized seeking in HQ-CT (-)
WHQD   :  0.65 : weight of internal drive in HQ-CT (-)
KHQ    :  2.75 : HQ-CT half-constant in drive units (-)
NHQ    :  1.60 : HQ-CT Hill coefficient (-)
TAUHQ  : 14.0  : HQ-CT adaptation time constant (d)
AGEHP  :  6.50 : midpoint of the Miller phase-3 hyperphagia gate (y)
WAGEHP :  1.40 : width of the hyperphagia gate (y)
AGEHP4 : 32.0  : midpoint of the partial phase-4 decline (y)
WAGEHP4:  6.0  : width of the phase-4 decline (y)
FHP4   :  0.22 : depth of the phase-4 decline (-)

// ---- 9. behaviour ----------------------------------------------------
BEH0   : 34.0  : behaviour-score scale (points)
TAUBEH : 45.0  : behaviour adaptation time constant (d)
KBFR   :  0.78 : frustration drive on behaviour (-)
KBOXT  :  0.16 : oxytocin-arm relief of behaviour (-)

// ---- 10. somatotropic axis ------------------------------------------
KSGHRH :  4.0  : GHRH tone formation gain (-)
KDGHRH :  2.0  : GHRH tone turnover (1/d)
TAUSST :  5.0  : somatostatin adaptation (d)
KSSTI  :  0.55 : IGF-1 drive on somatostatin (-)
KSSTG  :  0.60 : somatostatin inhibition of GH secretion (-)
KCLGH  :  5.5  : GH pool turnover (1/d)
KIGFFB : 260.0 : IGF-1 negative feedback on GH secretion (ng/mL)
KFMGH  : 26.0  : fat-mass suppression of GH secretion (kg)
KGIGF  :  3.00 : saturating GH to IGF-1 half-constant (GHREL units)
KDIGF  :  1.15 : IGF-1 turnover (1/d)
KBP3   :  2.9  : IGFBP-3 formation gain (mg/L/d)
KDBP3  :  1.05 : IGFBP-3 turnover (1/d)
FPOTGH :  0.63 : potency of a SMOOTH day-scale somatropin exposure (-)

// ---- 11. glucose / insulin ------------------------------------------
KEGP   : 2100.0: hepatic glucose production scale (mg/dL*L/d)
KIEGP  : 14.0  : insulin suppression of hepatic glucose output (uU/mL)
KGABS  :  0.055: glucose appearance per kcal emptied (mg/dL*L/kcal)
KGU0   :  6.5  : insulin-independent glucose disposal (1/d)
KGUI   :  0.93 : insulin-dependent glucose disposal (1/d per uU/mL)
VG     :  2.6  : glucose distribution volume (L/dL-scaled)
KGGLU  : 105.0 : glucose half-saturation of insulin secretion (mg/dL)
KSINS  : 88.0  : insulin secretion gain (uU/mL/d per unit flux)
VINS   :  1.0  : insulin distribution scale (-)
KCLINS :  6.0  : insulin day-scale pool turnover (1/d)
SI0    :  1.00 : reference insulin sensitivity (-)
KADPNSI:  0.42 : adiponectin drive on insulin sensitivity (-)
KFMSI  :  0.85 : fat-mass suppression of insulin sensitivity (-)
KFS    : 24.0  : half-constant of fat-mass insulin resistance (kg)
KGHSI  :  0.30 : GH-induced insulin resistance (-)
TAUA1C : 35.0  : HbA1c adaptation time constant (d)
KSADPN : 17.8  : adiponectin formation gain (ug/mL/d)
KDADPN :  1.1  : adiponectin turnover (1/d)
KFADP  : 30.0  : fat-mass suppression of adiponectin (kg)

// ---- 12. gonadal axis / bone ----------------------------------------
BAPUB  : 11.5  : bone age at the pubertal gate (y)
WPUB   :  1.10 : width of the pubertal gate (y)
SGNRH  :  1.0  : GnRH tone scale (-)
TAUGN  :  8.0  : GnRH adaptation (d)
KLH    :  4.6  : LH formation gain (IU/L/d)
KDLH   :  2.4  : LH turnover (1/d)
KSEX   : 360.0 : sex-steroid formation gain (ng/dL/d per IU/L)
KDSEX  :  1.15 : sex-steroid turnover (1/d)
KBMDU  :  0.00135 : bone formation gain (g/cm2/d)
KBMDR  :  0.00120 : bone resorption rate (1/d)
KBMDL  :  0.35 : mechanical-loading drive on bone formation (-)
KBMDI  :  0.22 : IGF-1 drive on bone formation (-)
KBMDH  :  0.55 : hypogonadal amplification of resorption (-)

// ---- 13. upper airway -----------------------------------------------
KLYU   :  0.075: IGF-1 drive on lymphoid hypertrophy (1/d)
KLYD   :  0.052: lymphoid involution rate (1/d)
LMAX   :  2.30 : maximum lymphoid volume (-)
KADL   :  3.4  : adaptation amplification of involution (-)
TAUAD  : 100.0 : adaptation time constant (d)
RMS0   :  1.00 : respiratory muscle strength scale (-)
TAURMS : 75.0  : respiratory muscle adaptation - THE SLOW CLOCK (d)
KRMSIGF:  0.28 : IGF-1 drive on respiratory muscle strength (-)
NRMS   :  0.60 : lean-mass exponent on respiratory muscle strength (-)
AHI0   :  1.6  : baseline apnoea-hypopnoea index (events/h)
KAHIL  :  7.0  : lymphoid drive on obstructive events (events/h)
KAHIF  : 11.0  : fat/lean-ratio drive on obstructive events (events/h)
FLR0   :  0.35 : fat:lean ratio below which fat costs nothing (-)
KAHIR  :  6.5  : respiratory-muscle protection (events/h)
KAHIC  :  6.0  : central (hypothalamic) apnoea drive (events/h)
TAUAHI :  7.0  : AHI adaptation time constant (d)

// ---- 14. muscle tone / scoliosis / thermoregulation ---------------
TAUTONE: 60.0  : muscle-tone adaptation (d)
KTIGF  :  0.15 : IGF-1 drive on muscle tone (-)
TONEA  :  0.78 : reference (control) muscle tone (-)
KCOB   :  0.058: growth-velocity drive on Cobb angle (deg/cm)
NCOB   :  1.30 : hypotonia exponent on Cobb progression (-)
TAUTEMP: 420.0 : thermoregulatory adaptation (d)

// ---- 15. drug PK (apparent, day-scale, AUC-calibrated) ------------
KAGH   :  5.5  : somatropin absorption (1/d)
KEGH   :  1.4  : somatropin elimination (1/d)
VGH    :  2.60 : somatropin apparent volume (L/kg)
KALG   :  0.42 : long-acting GH absorption (1/d)
KELG   :  0.55 : long-acting GH elimination (1/d)
VLG    :  0.95 : long-acting GH apparent volume (L/kg)
KACB   : 14.0  : carbetocin intranasal absorption (1/d)
KECB   :  8.0  : carbetocin elimination (1/d)
VCB    : 32.0  : carbetocin apparent volume (L)
FCB    :  0.055: carbetocin intranasal bioavailability (-)
KADZ   :  1.5  : diazoxide choline ER absorption (1/d)
KEDZ   :  0.594: diazoxide elimination, t1/2 ~ 28 h (1/d)
VDZ    :  0.62 : diazoxide apparent volume (L/kg)
KASM   :  3.0  : setmelanotide absorption (1/d)
KESM   :  1.45 : setmelanotide elimination (1/d)
VSM    : 52.0  : setmelanotide apparent volume (L)
KASG   :  0.30 : semaglutide absorption (1/d)
KESG   :  0.099: semaglutide elimination, t1/2 ~ 7 d (1/d)
VSG    : 12.5  : semaglutide apparent volume (L)
KAOC   :  0.075: octreotide LAR release (1/d)
KEOC   :  0.34 : octreotide elimination (1/d)
VOC    : 22.0  : octreotide apparent volume (L)
KALV   :  6.0  : livoletide absorption (1/d)
KELV   :  2.8  : livoletide elimination (1/d)
VLV    :  0.20 : livoletide apparent volume (L/kg)
KATS   :  0.35 : testosterone enanthate release (1/d)
KETS   :  0.90 : testosterone elimination (1/d)
VTS    : 1400.0: testosterone apparent volume (L)
KAMF   :  6.0  : metformin absorption (1/d)
KEMF   :  2.9  : metformin elimination (1/d)
VMF    : 63.0  : metformin apparent volume (L)

// ---- 16. drug PD -----------------------------------------------------
ECBMAX :  1.35 : carbetocin OXTR Emax (-)
ECB50  :  0.55 : carbetocin OXTR EC50 (ng/mL)
EV1AMX :  1.30 : carbetocin V1a counter-effect Emax (-)
EV1A50 :  2.40 : carbetocin V1a EC50 (ng/mL)
EDZMAX :  0.40 : diazoxide AgRP K-ATP Emax (-)
EDZ50  : 15.0  : diazoxide AgRP K-ATP EC50 (ug/mL)
EDZIMX :  0.70 : diazoxide beta-cell K-ATP Emax - THE AE (-)
EDZI50 : 22.0  : diazoxide beta-cell K-ATP EC50 (ug/mL)
ESMMAX :  3.20 : setmelanotide MC4R Emax (-)
ESM50  :  9.0  : setmelanotide MC4R EC50 (ng/mL)
ESGGE  :  1.20 : semaglutide gastric-emptying delay Emax (-)
ESG50  : 18.0  : semaglutide EC50 (ng/mL)
ESGSAT :  1.60 : semaglutide contribution to the GLP-1 arm (-)
ESGEI  :  0.22 : semaglutide direct intake suppression (-)
ESGSI  :  0.30 : semaglutide insulin-sensitivity gain (-)
ESGAM  :  0.45 : GLP-1R agonist inhibition of AgRP neurons (-)
EOCG   :  3.10 : octreotide ghrelin suppression Emax (-)
EOC50  :  1.40 : octreotide EC50 (ng/mL)
EOCGH  :  2.60 : octreotide GH-secretion suppression Emax (-)
EOCINS :  0.80 : octreotide insulin suppression Emax (-)
KLVAG  :  0.008: livoletide acyl-ghrelin antagonism (per ng/mL)
EMFSI  :  0.22 : metformin insulin-sensitivity Emax (-)
EMF50  :  1.10 : metformin EC50 (ug/mL)
EMFEI  :  0.035: metformin intake effect (-)

// ---- 17. simulation control ----------------------------------------
AGE0   :  0.0  : chronological age at t = 0 (y)

$CMT @annotated
// ---- prohormone pools (PC1/3 substrates) ---------------------------
PPOMC  : POMC precursor pool (-)
PPOXT  : pro-oxytocin precursor pool (-)
PPGHRH : pro-GHRH precursor pool (-)
PPINS  : proinsulin pool (-)
PPGHR  : proghrelin pool (-)
// ---- neuropeptide products / gut-brain signals ---------------------
AMSH   : alpha-MSH tone (-)
OXT    : PVN oxytocin tone (-)
AGRP   : AgRP/NPY orexigenic tone (-)
AG     : acylated ghrelin (pg/mL)
UAG    : unacylated ghrelin (pg/mL)
PYY    : PYY 3-36 (pmol/L)
GLP1E  : endogenous GLP-1 (pmol/L)
// ---- adiposity, energy, growth -------------------------------------
LEP    : leptin (ng/mL)
FM     : fat mass - THE ENERGY BUFFER (kg)
LFM    : fat-free mass (kg)
GUT    : upper-gut nutrient pool (kcal)
HT     : height (cm)
ACT    : physical activity level (-)
BA     : bone age (y)
// ---- food-seeking / behaviour --------------------------------------
SEEK   : food-seeking state - BISTABLE (-)
HQ     : HQ-CT hyperphagia score (points)
BEH    : maladaptive behaviour score (points)
// ---- somatotropic axis ---------------------------------------------
GHRH   : GHRH tone (-)
SST    : somatostatin tone (-)
GH     : growth hormone, 24-h mean (ug/L)
IGF1   : IGF-1 (ng/mL)
IGFBP3 : IGFBP-3 (mg/L)
// ---- glucose / insulin ---------------------------------------------
GLU    : plasma glucose (mg/dL)
INS    : plasma insulin (uU/mL)
HBA1C  : HbA1c (%)
ADPN   : adiponectin (ug/mL)
// ---- gonadal axis / bone -------------------------------------------
GNRH   : GnRH tone (-)
LH     : LH (IU/L)
SEX    : endogenous sex steroid (ng/dL)
BMD    : lumbar bone mineral density (g/cm2)
// ---- upper airway ---------------------------------------------------
LYMPH  : adenotonsillar/lymphoid volume - FAST CLOCK (-)
ADAPT  : lymphoid adaptation/involution signal (-)
RMS    : respiratory muscle strength - SLOW CLOCK (-)
AHI    : apnoea-hypopnoea index (events/h)
// ---- musculoskeletal / thermoregulation ----------------------------
TONE   : muscle tone (-)
COBB   : scoliosis Cobb angle (deg)
TEMP   : thermoregulatory instability index (-)
// ---- drug PK --------------------------------------------------------
AGHD   : somatropin SC depot (mg)
AGHC   : somatropin central (mg)
ALGD   : long-acting GH depot (mg)
ALGC   : long-acting GH central (mg)
ACBD   : carbetocin intranasal depot (mg)
ACBC   : carbetocin central (mg)
ADZD   : diazoxide choline depot (mg)
ADZC   : diazoxide central (mg)
ASMD   : setmelanotide depot (mg)
ASMC   : setmelanotide central (mg)
ASGD   : semaglutide depot (mg)
ASGC   : semaglutide central (mg)
AOCD   : octreotide LAR depot (mg)
AOCC   : octreotide central (mg)
ALVD   : livoletide depot (mg)
ALVC   : livoletide central (mg)
ATSD   : testosterone depot (mg)
ATSC   : testosterone central (mg)
AMFD   : metformin depot (mg)
AMFC   : metformin central (mg)
// ---- diagnostics ----------------------------------------------------
CUMEI  : cumulative energy intake (kcal)
CUMTEE : cumulative energy expenditure (kcal)

$MAIN
// Initial conditions are FUNCTIONS OF THE LESION, so they live here and
// not in $CMT: a control and a PWS patient start from different
// steady states of the same equations.
double pc  = 1.0 - DPC13 * LES;
double fox = (LES > 0.5) ? FOXTN : 1.0;

PPOMC_0  = 1.0 / (pc + EPSPOMC);
PPOXT_0  = fox / (pc + EPSPOXT);
PPGHRH_0 = 1.0 / (pc + EPSPGHRH);
PPINS_0  = 1.0 / (pc + EPSPINS);
PPGHR_0  = 1.0 / (pc + EPSPGHR);
AMSH_0   = KSAMSH / KDAMSH * pc * PPOMC_0;
OXT_0    = KSOXT  / KDOXT  * pc * PPOXT_0;
AGRP_0   = 1.0;
AG_0     = 350.0;
UAG_0    = 1250.0;
PYY_0    = PYY0;
GLP1E_0  = GLP10;
LEP_0    = KSLEP * itp(FM_T, AGE0) / KDLEP;
FM_0     = 0.5;
LFM_0    = 2.7;
GUT_0    = (schofield(AGE0, itp(LFM_T,AGE0) + itp(FM_T,AGE0))
            * (1.0 + PALX*ACTREF) / (1.0 - FDIT)) / KEMPT;
HT_0     = 50.0;
ACT_0    = 0.60;
BA_0     = 0.0;
SEEK_0   = 0.05;
HQ_0     = 1.0;
BEH_0    = 8.0;
GHRH_0   = KSGHRH / KDGHRH * pc * PPGHRH_0;
SST_0    = 1.0;
GH_0     = itp(GHN_T, AGE0);
IGF1_0   = itp(IGF_T, AGE0);
IGFBP3_0 = 2.0;
GLU_0    = 82.0;
INS_0    = 8.0;
HBA1C_0  = 4.6;
ADPN_0   = 12.0;
GNRH_0   = 0.0;
LH_0     = 0.0;
SEX_0    = 8.0;
BMD_0    = 0.42;
LYMPH_0  = itp(LY_T, AGE0);
ADAPT_0  = 0.0;
RMS_0    = 1.0;
AHI_0    = 2.0;
TONE_0   = (LES > 0.5) ? TONEA * FTONE0 : TONEA;
COBB_0   = 2.0;
TEMP_0   = (LES > 0.5) ? 1.0 : 0.0;

$ODE
// =====================================================================
//  ORDER MATTERS AND IS NOT ARBITRARY.
//  The oxytocin relay arm (X1) must be known before the energy block,
//  because the adiposity negative feedback on intake is TRANSDUCED by
//  that relay; and the realized intake must be known before the
//  remaining satiety arms, because meal load drives distension, PYY and
//  GLP-1.  Getting this order wrong silently changes the model.
// =====================================================================
AGEY = AGE0 + SOLVERTIME / 365.25;
PC13 = 1.0 - DPC13 * LES;
FMx  = (FM  > FMMIN) ? FM  : FMMIN;
LFMx = (LFM > 1.0)   ? LFM : 1.0;
BW   = FMx + LFMx;
double pws = (LES > 0.5) ? 1.0 : 0.0;

double LFMN = itp(LFM_T, AGEY);
double FMN  = itp(FM_T,  AGEY);
double IGFN = itp(IGF_T, AGEY);
double GHN  = itp(GHN_T, AGEY);
double LYN  = itp(LY_T,  AGEY);
double LEPN = KSLEP * FMN / KDLEP;
double REEN = schofield(AGEY, LFMN + FMN);

// ---- (1) drug concentrations, apparent day-scale --------------------
CGH = AGHC / (VGH * BW) * 1000.0;                 // ug/L
CLG = ALGC / (VLG * BW) * 1000.0;                 // ug/L
CCB = ACBC / VCB * 1000.0;                        // ng/mL
CDZ = ADZC / (VDZ * BW);                          // ug/mL
CSM = ASMC / VSM * 1000.0;                        // ng/mL
CSG = ASGC / VSG * 1000.0;                        // ng/mL
COC = AOCC / VOC * 1000.0;                        // ng/mL
CLV = ALVC / (VLV * BW) * 1000.0;                 // ng/mL
CTS = ATSC / VTS * 100000.0;                      // ng/dL
CMF = AMFC / VMF;                                 // ug/mL

// ---- (2) receptor-level effects.  Note ECB is BIPHASIC: agonism minus
//      V1a cross-activation.  Its optimum is analytic (see the header
//      of pws_calibration.py) and is what inverts the dose-response.
ECB   = ECBMAX*CCB/(ECB50+CCB) - EV1AMX*CCB/(EV1A50+CCB);
EDZx  = EDZMAX*CDZ/(EDZ50+CDZ);
EDZI  = EDZIMX*CDZ/(EDZI50+CDZ);
ESMx  = ESMMAX*CSM/(ESM50+CSM);
ESGG  = ESGGE *CSG/(ESG50+CSG);
ESGS  = ESGSAT*CSG/(ESG50+CSG);
ESGA  = ESGAM *CSG/(ESG50+CSG);
EOCg  = EOCG  *COC/(EOC50+COC);
EOCgh = EOCGH *COC/(EOC50+COC);
EOCin = EOCINS*COC/(EOC50+COC);
EMFx  = EMFSI *CMF/(EMF50+CMF);

// ---- (3) somatotropic axis, reference-normalized --------------------
double fbI = (1.0 + IGFN/KIGFFB) / (1.0 + IGF1/KIGFFB);
double fbF = (1.0 + FMN /KFMGH)  / (1.0 + FMx /KFMGH);
double sstx = (SST > 0.05) ? SST : 0.05;
double fbS = (1.0 + KSSTG) / (1.0 + KSSTG*sstx);
GHSEC = KCLGH * GHN * GHRH * ((pws>0.5)?FSOM:1.0) * fbI*fbF*fbS/(1.0+EOCgh);
GHREL = (GH + FPOTGH*(CGH + CLG)) / GHN;
double gex = (GHREL > 1.0) ? GHREL - 1.0 : 0.0;
GHX = gex / (1.0 + gex);
double gg = (GHREL > 1e-4) ? GHREL : 1e-4;
// [DEFECT 3] a power law here gave IGF-1 SDS +20; GH -> IGF-1 saturates.
IGFT   = IGFN * (gg/(KGIGF+gg)) * (KGIGF+1.0);
IGFdrv = IGF1 / ((IGFN > 1.0) ? IGFN : 1.0);
EXC = IGFdrv - 1.0; if (EXC < 0.0) EXC = 0.0; if (EXC > 1.5) EXC = 1.5;

// ---- (4) the melanocortin -> PVN-oxytocin relay ---------------------
// alpha-MSH is an INPUT GAIN here, not a parallel arm.  relay(1) = 1 and
// relay(inf) = KREL+1, so MC4R agonism has a hard CEILING of 1.2x and
// starts from PWS's own alpha-MSH tone.  This is the whole reason MC4R
// agonism works in POMC/LEPR deficiency (lesion ABOVE MC4R) and not in
// PWS (lesion BELOW it).
MCIN   = AMSH + ESMx;
RELAYx = (MCIN/(KREL+MCIN)) * (KREL+1.0);
X1 = OXT * RELAYx * (1.0 + ECB);
if (X1 < 1e-4) X1 = 1e-4;
// RELAY INTEGRITY: leptin's anorexigenic signal is transduced by this
// same relay, so PWS is FUNCTIONALLY LEPTIN-RESISTANT with normal leptin
// and normal receptors, at no extra parameter.
RI = X1 / X1CTL; if (RI > 1.20) RI = 1.20;

// ---- (5) energy balance and what is actually eaten ------------------
// [DEFECT 2] Schofield reference scaled by composition, not Cunningham.
REE = REEN * (WREEL*LFMx/LFMN + WREEF*FMx/FMN) * (1.0 + KGHRMR*GHX);
// [DEFECT 6] only EXCESS fat costs mobility, else a control losing fat
// gains expenditure and drifts to 2% body fat.
double fx = FMx - FMN; if (fx < 0.0) fx = 0.0;
ACTt = ACT0 * TONE/TONEA * (1.0 - KFATACT*fx/(fx+KFA));
ATh  = REE * PALX * ACT;
EIREQ = REEN * (1.0 + PALX*ACTREF) / (1.0 - FDIT);
double offs = EIBF + KEIS*SEEK + KEIH*HQ/HQMAX - KLEPEI*RI*(LEP/LEPN - 1.0);
offs *= (1.0 - ESGEI*ESGS/ESGSAT - EMFEI*EMFx/EMFSI);
EIDRV = EIREQ * offs; if (EIDRV < 0.0) EIDRV = 0.0;
// What is OFFERED.  TITR = 0: age-normative portions that take no notice
// of the child's weight - the phase-2a situation, and the reason a NORMAL
// intake is a surplus for a low-expenditure child.  TITR = 1: caregivers
// titrate toward a weight-FOR-HEIGHT target, the way a PWS clinic does;
// then WEIGHT is pinned and COMPOSITION is set by the lean target, which
// is how growth hormone lowers %fat while barely moving BMI.
double ha = height_age(HT);
BWTGT = BWTGTR * (itp(LFM_T,ha) + itp(FM_T,ha));
TITRx = 1.0 + TITR * KTITR * (BWTGT/BW - 1.0);
if (TITRx > 1.50) TITRx = 1.50; if (TITRx < 0.45) TITRx = 0.45;
EICAP = AVAIL * EIREQ * TITRx;
double suck = (pws > 0.5) ? itp(SUCK_T, AGEY) : 1.0;
EIx = EIDRV; if (EICAP < EIx) EIx = EICAP;
double hard = EICAPF * EIREQ; if (hard < EIx) EIx = hard;
EIx *= suck; if (EIx < 0.0) EIx = 0.0;
ASKx = EIREQ * (EIBF + KEIS*SEEK + KEIH*HQ/HQMAX);
FRUSTx = (ASKx - EIx) / ((ASKx > 1.0) ? ASKx : 1.0);
if (FRUSTx < 0.0) FRUSTx = 0.0; if (FRUSTx > 1.0) FRUSTx = 1.0;
REINFx = 1.0 - FRUSTx;
DITh = FDIT * EIx;
TEEx = REE + ATh + DITh;
ESUFF = EIx / ((TEEx > 1.0) ? TEEx : 1.0);
if (ESUFF > 1.0) ESUFF = 1.0; if (ESUFF < 0.45) ESUFF = 0.45;

// ---- (6) meal load.  [DEFECT 5] the mean gastric pool cannot separate
//      meal-related distension from chronic overfill on a day scale.
MLOAD = EIx / EIREQ;
if (MLOAD > 1.30) MLOAD = 1.30; if (MLOAD < 0.50) MLOAD = 0.50;
MLOAD *= (1.0 + 0.25*ESGG/ESGGE);
FASTs = 1.0 + KFASTA*(1.0 - MLOAD); if (FASTs < 0.35) FASTs = 0.35;
FASTg = 1.0 + KFASTG*(1.0 - MLOAD); if (FASTg < 0.35) FASTg = 0.35;

// ---- (7) prohormone synthesis and processing fluxes ----------------
LEPSAT = (LEP/(KLEPR0+LEP)) / (LEPN/(KLEPR0+LEPN));
SPOMCx = SPOMC0 * (1.0 + KLEPP*(LEPSAT - 1.0));
SPOXTx = SPOXT0 * ((pws>0.5) ? FOXTN : 1.0);
PHIG   = GLU*GLU/(KGGLU*KGGLU + GLU*GLU);
SPINSx = SPINS0 * PHIG / 0.3959;
SPGHRx = SPGHR0 * ((pws>0.5)?FGHRC:1.0) * FASTg
         / (1.0 + INS/KIINS) / (1.0 + EOCg);
fPOMC  = KCPOMC  * PC13 * PPOMC;
fPOXT  = KCPOXT  * PC13 * PPOXT;
fPGHRH = KCPGHRH * PC13 * PPGHRH;
fPINS  = KCPINS  * PC13 * PPINS;
fPGHR  = KCPGHR  * PC13 * PPGHR;

dxdt_PPOMC  = SPOMCx  - (KCPOMC  * PC13 + EPSPOMC ) * PPOMC;
dxdt_PPOXT  = SPOXTx  - (KCPOXT  * PC13 + EPSPOXT ) * PPOXT;
dxdt_PPGHRH = SPGHRH0 - (KCPGHRH * PC13 + EPSPGHRH) * PPGHRH;
dxdt_PPINS  = SPINSx  - (KCPINS  * PC13 + EPSPINS ) * PPINS;
dxdt_PPGHR  = SPGHRx  - (KCPGHR  * PC13 + EPSPGHR ) * PPGHR;

// ---- (8) neuropeptide tones.  [DEFECT 1] KSAMSH/KSOXT normalize these
//      to 1.0 in a control; without that the CONTROL latches.
dxdt_AMSH = KSAMSH * fPOMC - KDAMSH * AMSH;
dxdt_OXT  = KSOXT  * fPOXT - KDOXT  * OXT;
// OXTR agonism acts on BOTH sides of the synapse: post-synaptically on
// the satiety arm (X1 above) and PRE-synaptically on AgRP (here).  This
// is the asymmetry that separates carbetocin from setmelanotide.
double oxt_eff = OXT * (1.0 + ECB);
LEPF = (1.0 + LEPN/KLEPA) / (1.0 + RI*LEP/KLEPA);
dxdt_AGRP = SAGRP*FASTs*LEPF/(1.0 + EDZx + ESGA)
            - KDAGRP*AGRP*(1.0 + KINHM*AMSH + KINHO*oxt_eff);

// ---- (9) ghrelin ----------------------------------------------------
FACx = FACYL * ((pws>0.5)?FACYLP:1.0); if (FACx > 0.60) FACx = 0.60;
dxdt_AG  = KSAG*fPGHR*FACx - (KCLAG + KDEACYL)*AG;
dxdt_UAG = KSUAG*fPGHR*(1.0-FACx) + KDEACYL*AG - KCLUAG*UAG;

// ---- (10) gut peptides ---------------------------------------------
dxdt_PYY   = SPYY  * MLOAD - KDPYY  * PYY;
dxdt_GLP1E = SGLP1 * MLOAD - KDGLP1 * GLP1E;

// ---- (11) the satiety integrator: a HARMONIC MEAN -------------------
X2 = ((pws>0.5)?FVAG:1.0) * MLOAD;              if (X2 < 1e-4) X2 = 1e-4;
X3 = ((pws>0.5)?FPYYS:1.0) * PYY/PYY0;          if (X3 < 1e-4) X3 = 1e-4;
X4 = GLP1E/GLP10 + ESGS;                        if (X4 < 1e-4) X4 = 1e-4;
double lepq = (LEP > 0.15) ? LEP : 0.15;
double insq = (INS > 1.0)  ? INS : 1.0;
X5 = pow(lepq/LEPN, 0.28) * pow(insq/11.0, 0.22);
if (X5 < 1e-4) X5 = 1e-4;
SATx = 1.0 / (W1/X1 + W2/X2 + W3/X3 + W4/X4 + W5/X5);

// ---- (12) orexigenic drive; the ghrelin arm SATURATES ---------------
AGEFF = AG / (1.0 + KLVAG*CLV);
double ghr  = AGEFF/(KAG50+AGEFF);
double ghr0 = 350.0/(KAG50+350.0);
GHRARM = ghr/ghr0;
GATEHP = (1.0/(1.0+exp(-(AGEY-AGEHP)/WAGEHP)))
         * (1.0 - FHP4/(1.0+exp(-(AGEY-AGEHP4)/WAGEHP4)));
double agq = (AGRP > 0.02) ? AGRP : 0.02;
DRVI = pow(agq, NAGRP) * (1.0 + KGHRD*(GHRARM - 1.0)) * GATEHP;
if (DRVI < 0.0) DRVI = 0.0;
DRVE = DRVI * CUE;

// ---- (13) food-seeking: BISTABLE -----------------------------------
// Self-sensitization is REINFORCED by food actually obtained, so removing
// ACCESS (not merely cues) is what can destroy the latched upper state.
SELFGx = KSELF * (RFL + (1.0-RFL)*REINFx);
double sh = pow(SEEK, NSELF);
dxdt_SEEK = KSON*DRVE*(1.0-SEEK)
            + SELFGx * sh/(pow(KSHALF,NSELF)+sh) * (1.0-SEEK)
            - KSOFF*SEEK*SATx;

// ---- (14) HQ-CT -----------------------------------------------------
double dh = pow(DRVI, NHQ);
double HQt = HQMAX * (WHQS*SEEK + WHQD*dh/(pow(KHQ,NHQ)+dh));
dxdt_HQ = (HQt - HQ) / TAUHQ;

// ---- (15) behaviour: the unrealized drive IS the frustration term ---
double BEHt = BEH0 * FSUB * (1.0 + KBFR*FRUSTx)
              * (1.0 - KBOXT*(X1/X1REF - 1.0)) * (0.55 + 0.45*GATEHP);
if (BEHt < 0.0) BEHt = 0.0;
dxdt_BEH = (BEHt - BEH) / TAUBEH;

// ---- (16) leptin, gut pool, activity, bone age ----------------------
dxdt_LEP = KSLEP*FMx - KDLEP*LEP;
dxdt_GUT = EIx - KEMPT*((pws>0.5)?FGE:1.0)/(1.0+ESGG)*GUT;
dxdt_ACT = (ACTt - ACT) / TAUACT;

// ---- (17) growth and the lean-mass target --------------------------
SEXTOT = SEX + CTS;
SEXREL = (SEXTOT/(SEXTOT+KSEXH)) / (620.0/(620.0+KSEXH));
double sxr = itp(SXR_T, AGEY); if (sxr < 0.03) sxr = 0.03;
SEXFAC = SEXREL / sxr; if (SEXFAC > 1.30) SEXFAC = 1.30;
EPIg = 1.0/(1.0 + exp((BA - BACLOSE)/0.35));
double igq = (IGFdrv > 0.15) ? IGFdrv : 0.15;
double igc = (igq < GVCAP) ? igq : GVCAP;
GVF = GVFL + (1.0-GVFL)*igc;
GVx = (itp(GVB_T,AGEY)*GVF + itp(GVS_T,AGEY)*SEXFAC*sqrt(GVF)) * ESUFF * EPIg;
dxdt_HT = GVx / 365.25;
dxdt_BA = (1.0/365.25) * (1.0 + KBASEX*((SEXREL<1.2)?SEXREL:1.2));
// Every lean-target factor equals 1.0 for the non-PWS reference at every
// age, so the control tracks LFMN by construction and every deviation is
// a NAMED mechanism.  GH and testosterone reach lean mass THROUGH these
// factors; the model gives them no private channel.
double igl = (IGFdrv < LIGFCAP) ? IGFdrv : LIGFCAP;
FMODx = (WLIGF0 + (1.0-WLIGF0)*igl)
        * (WLTON0 + (1.0-WLTON0)*TONE/TONEA)
        * sqrt(ESUFF)
        * (1.0 - 0.5*KSEXPART + 0.5*KSEXPART*SEXFAC);
if (FMODx > FMODCAP) FMODx = FMODCAP;
LFMT = LFMN * FMODx;
// [DEFECT 4] lean mass used to be accrued BOTH from the energy partition
// and from a height-growth term, double-counting it.  Lean now follows a
// developmental target and FAT IS THE ENERGY BUFFER, so energy is exactly
// conserved and PWS's lean deficit is DERIVED, not imposed.
dxdt_LFM = (LFMT - LFM) / TAULFM;
dxdt_FM  = (EIx - TEEx - RHOL*dxdt_LFM) / RHOF;
dxdt_CUMEI  = EIx;
dxdt_CUMTEE = TEEx;

// ---- (18) somatotropic ODEs ----------------------------------------
dxdt_GHRH = KSGHRH*fPGHRH - KDGHRH*GHRH;
dxdt_SST  = (1.0 + KSSTI*(IGFdrv - 1.0) - SST) / TAUSST;
dxdt_GH   = GHSEC - KCLGH*GH;
double liver = 0.55 + 0.45*ESUFF;
dxdt_IGF1   = KDIGF * (IGFT*liver - IGF1);
dxdt_IGFBP3 = KBP3 * pow(gg, 0.65) - KDBP3*IGFBP3;

// ---- (19) glucose / insulin ----------------------------------------
dxdt_ADPN = KSADPN*((pws>0.5)?FADPN:1.0)/(1.0+FMx/KFADP) - KDADPN*ADPN;
SIx = SI0 * (1.0 + KADPNSI*(ADPN/11.0 - 1.0))
      / (1.0 + KFMSI*FMx/(FMx+KFS))
      / (1.0 + KGHSI*GHX)
      * (1.0 + EMFx + ESGSI*ESGS/ESGSAT);
if (SIx < 0.15) SIx = 0.15;
INSSEC = KSINS*fPINS/(1.0+EDZI)/(1.0+EOCin);
dxdt_INS = INSSEC/VINS - KCLINS*INS;
double EGP = KEGP/(1.0 + INS/KIEGP);
double Ra  = KGABS*KEMPT*((pws>0.5)?FGE:1.0)*GUT;
double Rd  = (KGU0 + KGUI*INS*SIx)*GLU;
dxdt_GLU   = (EGP + Ra - Rd) / VG;
dxdt_HBA1C = ((GLU*1.14 + 46.7)/28.7 - HBA1C) / TAUA1C;

// ---- (20) gonadal axis / bone --------------------------------------
// MKRN3 is a pubertal BRAKE, so losing it ADVANCES the gate, while the
// hypothalamic lesion LOWERS the amplitude: early onset, low amplitude,
// non-progressive puberty, both signs from the same deletion.
double bapub = BAPUB - ((pws>0.5)?DMKRN3:0.0);
double gateP = 1.0/(1.0 + exp(-(BA - bapub)/WPUB));
dxdt_GNRH = (SGNRH*gateP*((pws>0.5)?FHYPO:1.0) - GNRH) / TAUGN;
dxdt_LH   = KLH*GNRH - KDLH*LH;
dxdt_SEX  = KSEX*LH - KDSEX*SEX;
double sxq = (SEXTOT > 0.02) ? SEXTOT : 0.02;
double hypo = 1.0 - sxq/620.0; if (hypo < 0.0) hypo = 0.0;
dxdt_BMD = KBMDU*sqrt(sxq/620.0)*(1.0 + KBMDL*LFMx/LFMN)*(1.0 + KBMDI*EXC)
           - KBMDR*BMD*(1.0 + KBMDH*hypo);

// ---- (21) upper airway: TWO CLOCKS ---------------------------------
dxdt_ADAPT = (EXC - ADAPT) / TAUAD;
dxdt_LYMPH = KLYU*EXC*(LMAX - LYMPH) - KLYD*(LYMPH - LYN)*(1.0 + KADL*ADAPT);
RMST = pow(LFMx/LFMN, NRMS) * (1.0 + KRMSIGF*EXC);
dxdt_RMS = (RMS0*RMST - RMS) / TAURMS;
// FOSA is an ANATOMICAL susceptibility, so it multiplies only the
// OBSTRUCTIVE terms.  That is what makes it amplify the GH WINDOW and
// not the plateau.
double dly = LYMPH - LYN; if (dly < 0.0) dly = 0.0;
double flr = FMx/LFMx - FLR0; if (flr < 0.0) flr = 0.0;
AHIOBS = KAHIL*dly + KAHIF*flr;
AHIT = AHI0 + FOSA*AHIOBS - KAHIR*(RMS - 1.0)
       + KAHIC*itp(CAP_T,AGEY)*LES*(1.0 - 0.45*GHX);
if (AHIT < 0.3) AHIT = 0.3;
dxdt_AHI = (AHIT - AHI) / TAUAHI;

// ---- (22) tone, scoliosis, thermoregulation ------------------------
// hypotonia improves through infancy but never resolves in PWS
double amin = (AGEY/5.0 < 1.0) ? AGEY/5.0 : 1.0;
TONE0 = (pws>0.5) ? TONEA*(FTONE0 + (FTONEP-FTONE0)*amin) : TONEA;
TONET = TONE0 * (1.0 + KTIGF*EXC) * pow(LFMx/LFMN, 0.40);
dxdt_TONE = (TONET - TONE) / TAUTONE;
double toneq = (TONE > 0.25) ? TONE : 0.25;
// GH raises growth velocity (worse) AND muscle tone (better); the two
// nearly cancel, which is what the GH trials found.
dxdt_COBB = KCOB*(GVx/365.25)*pow(TONEA/toneq, NCOB);
double TEMPT = LES * (0.35 + 0.65*exp(-AGEY/6.0));
dxdt_TEMP = (TEMPT - TEMP) / TAUTEMP;

// ---- (23) drug PK ---------------------------------------------------
dxdt_AGHD = -KAGH*AGHD;   dxdt_AGHC = KAGH*AGHD - KEGH*AGHC;
dxdt_ALGD = -KALG*ALGD;   dxdt_ALGC = KALG*ALGD - KELG*ALGC;
dxdt_ACBD = -KACB*ACBD;   dxdt_ACBC = KACB*ACBD*FCB - KECB*ACBC;
dxdt_ADZD = -KADZ*ADZD;   dxdt_ADZC = KADZ*ADZD - KEDZ*ADZC;
dxdt_ASMD = -KASM*ASMD;   dxdt_ASMC = KASM*ASMD - KESM*ASMC;
dxdt_ASGD = -KASG*ASGD;   dxdt_ASGC = KASG*ASGD - KESG*ASGC;
dxdt_AOCD = -KAOC*AOCD;   dxdt_AOCC = KAOC*AOCD - KEOC*AOCC;
dxdt_ALVD = -KALV*ALVD;   dxdt_ALVC = KALV*ALVD - KELV*ALVC;
dxdt_ATSD = -KATS*ATSD;   dxdt_ATSC = KATS*ATSD - KETS*ATSC;
dxdt_AMFD = -KAMF*AMFD;   dxdt_AMFC = KAMF*AMFD - KEMF*AMFC;

$TABLE
double AGEOUT = AGE0 + SOLVERTIME/365.25;
double LFMNO  = itp(LFM_T, AGEOUT);
double FMNO   = itp(FM_T,  AGEOUT);
double HTNO   = itp(HT_T,  AGEOUT);
double IGFNO  = itp(IGF_T, AGEOUT);
double FMo  = (FM  > FMMIN) ? FM  : FMMIN;
double LFMo = (LFM > 1.0)   ? LFM : 1.0;
double BWo  = FMo + LFMo;

capture AGEY_OUT = AGEOUT;
capture BODYWT   = BWo;
capture PBF      = 100.0 * FMo / BWo;
capture BMI      = BWo / pow(HT/100.0, 2.0);
capture HTSDS    = (HT - HTNO) / (1.9 + 0.031*HTNO);
capture IGFSDS   = (IGF1 - IGFNO) / (0.32*IGFNO);
capture WTSDS    = (BWo - (LFMNO+FMNO)) / (0.13*(LFMNO+FMNO));
capture AGUAG    = AG / ((UAG > 1e-6) ? UAG : 1e-6);
capture AGUAGEFF = AGEFF / ((UAG > 1e-6) ? UAG : 1e-6);
capture HOMAIR   = GLU * INS / 405.0;
capture LEANDEF  = 100.0 * (LFMo/LFMNO - 1.0);
capture TEEFRAC  = TEEx / EIREQ;
capture EIFRAC   = EIx  / EIREQ;

// the structural quantities, exposed so a Shiny app can plot the
// mechanism and not merely the endpoints
capture SAT_OUT   = SATx;
capture ARM_OXT   = X1;
capture ARM_VAGAL = X2;
capture ARM_PYY   = X3;
capture ARM_GLP1  = X4;
capture ARM_LEPINS= X5;
capture RELAY_INT = RI;
capture DRIVE_INT = DRVI;
capture DRIVE_EFF = DRVE;
capture GHREL_OUT = GHREL;
capture GHRELIN_ARM = GHRARM;
capture SELF_GAIN = SELFGx;
capture REINFORCE = REINFx;
capture FRUSTRATION = FRUSTx;
capture MEALLOAD  = MLOAD;
capture LEAN_TGT  = LFMT;
capture WT_TGT    = BWTGT;
capture EI_OUT    = EIx;
capture TEE_OUT   = TEEx;
capture REE_OUT   = REE;
capture GV_OUT    = GVx;
capture PC13_OUT  = PC13;

// the five branch losses and precursor ratios, in closed form, so the
// map in cluster 3 of pws_qsp_model.dot can be checked at run time
capture L_POXT = 1.0 - (1.0+EPSPOXT )/(1.0+EPSPOXT /PC13);
capture L_POMC = 1.0 - (1.0+EPSPOMC )/(1.0+EPSPOMC /PC13);
capture L_GHRH = 1.0 - (1.0+EPSPGHRH)/(1.0+EPSPGHRH/PC13);
capture L_PINS = 1.0 - (1.0+EPSPINS )/(1.0+EPSPINS /PC13);
capture L_PGHR = 1.0 - (1.0+EPSPGHR )/(1.0+EPSPGHR /PC13);
capture PREC_PROD = 1.0/PC13;   // identical for EVERY branch, exactly

// drug exposures
capture C_GH = CGH;  capture C_CB = CCB;  capture C_DZ = CDZ;
capture C_SM = CSM;  capture C_SG = CSG;  capture C_OC = COC;
capture C_LV = CLV;  capture C_TS = CTS;  capture C_MF = CMF;
capture E_OXTR = ECB;  capture E_AGRP_KATP = EDZx;
capture E_BCELL_KATP = EDZI;  capture E_MC4R = ESMx;

$CAPTURE @annotated
AG      : acylated ghrelin (pg/mL)
HQ      : HQ-CT hyperphagia score (points)
SEEK    : food-seeking state (-)
AHI     : apnoea-hypopnoea index (events/h)
IGF1    : IGF-1 (ng/mL)
HT      : height (cm)
FM      : fat mass (kg)
LFM     : fat-free mass (kg)
HBA1C   : HbA1c (%)
INS     : insulin (uU/mL)
ADPN    : adiponectin (ug/mL)
BMD     : bone mineral density (g/cm2)
COBB    : Cobb angle (deg)
TONE    : muscle tone (-)
BEH     : behaviour score (points)
LYMPH   : lymphoid volume (-)
RMS     : respiratory muscle strength (-)
SEX     : sex steroid (ng/dL)
LH      : LH (IU/L)
OXT     : oxytocin tone (-)
AMSH    : alpha-MSH tone (-)
AGRP    : AgRP tone (-)

## =====================================================================
##  SCENARIO DRIVER
##  ---------------
##  Every scenario below is reproduced numerically in
##  pws_reference_output.txt by pws_reference_model.py.  Compartment
##  numbers for dosing:
##      AGHD 43  ALGD 45  ACBD 47  ADZD 49  ASMD 51
##      ASGD 53  AOCD 55  ALVD 57  ATSD 59  AMFD 61
##  (the depot of each pair; use `cmt="AGHD"` if you prefer names)
##
##  library(mrgsolve); library(dplyr); library(tidyr); library(ggplot2)
##  mod <- mread("pws_mrgsolve_model.R")
##
##  ## ---- the four food environments -------------------------------
##  ENV_CTRL  <- list(AVAIL = 1.40, CUE = 1.00, TITR = 0)   # non-PWS
##  ENV_FREE  <- list(AVAIL = 1.40, CUE = 1.00, TITR = 0)   # unmanaged PWS
##  ENV_MGMT  <- list(AVAIL = 1.00, CUE = 0.28, TITR = 1)   # standard care
##  ENV_TIGHT <- list(AVAIL = 0.90, CUE = 0.10, TITR = 1)   # strict
##
##  YEARS <- function(y) y * 365.25
##  gh_daily <- function(from_y, to_y, mgkg = 0.035)
##    ev(time = YEARS(from_y), amt = 0, cmt = "AGHD", ii = 1,
##       addl = floor(YEARS(to_y - from_y)), rate = 0) |>
##    mutate(amt = mgkg * 30)   # or set amt per body weight in a loop
##
##  ## ---- S1  non-PWS control, birth -> 25 y ------------------------
##  s1 <- mod |> param(c(LES = 0, ENV_CTRL)) |>
##          mrgsim(end = YEARS(25), delta = 7)
##
##  ## ---- S2  PWS, unmanaged (free access) --------------------------
##  s2 <- mod |> param(c(LES = 1, ENV_FREE)) |>
##          mrgsim(end = YEARS(25), delta = 7)
##
##  ## ---- S3  PWS, standard management, no growth hormone -----------
##  s3 <- mod |> param(c(LES = 1, ENV_MGMT)) |>
##          mrgsim(end = YEARS(25), delta = 7)
##
##  ## ---- S4  PWS, management + somatropin from 1.0 y ---------------
##  ##   expect: %fat 16.5 vs 29.0, height 169.8 vs 150.9 cm at 25 y,
##  ##           HQ-CT essentially unchanged (11.5 vs 12.0)
##  gh <- ev(time = YEARS(1), amt = 0.7, cmt = "AGHD", ii = 1, addl = 8760)
##  s4 <- mod |> param(c(LES = 1, ENV_MGMT)) |> ev(gh) |>
##          mrgsim(end = YEARS(25), delta = 7)
##
##  ## ---- S5  the Miller nutritional phases -------------------------
##  ##   normal portions, NOT titrated: weight-for-age SDS turns positive
##  ##   years before HQ-CT rises.  This is phase 2a, and it needs no
##  ##   appetite parameter to change.
##  s5 <- mod |> param(LES = 1, AVAIL = 1.00, CUE = 1.00, TITR = 0) |>
##          mrgsim(end = YEARS(20), delta = 7)
##
##  ## ---- S6  THE AIRWAY WINDOW -------------------------------------
##  ##   start GH at 6 y and sample weekly; AHI peaks at weeks 4-8 and
##  ##   then falls below baseline.  Re-run with FOSA = 1.9 for a
##  ##   vulnerable airway: 11.2 -> 22.1 -> 4.2 events/h.
##  ghw <- ev(time = YEARS(6), amt = 0.7, cmt = "AGHD", ii = 1, addl = 7300)
##  s6a <- mod |> param(c(LES = 1, ENV_MGMT, FOSA = 1.0)) |> ev(ghw) |>
##           mrgsim(end = YEARS(8.5), delta = 3.5)
##  s6b <- mod |> param(c(LES = 1, ENV_MGMT, FOSA = 1.9)) |> ev(ghw) |>
##           mrgsim(end = YEARS(8.5), delta = 3.5)
##
##  ## ---- S7  the drug panel at 12 y (all on background GH) ---------
##  ##   expected placebo-corrected dHQ-CT after 8-13 weeks:
##  ##     carbetocin 3.2 mg TID   -1.85    (AG   +1.5%)
##  ##     carbetocin 9.6 mg TID   -1.38    (AG   +1.1%)
##  ##     DCCR 5.1 mg/kg          -1.54    (AG   +5.3%)
##  ##     setmelanotide 3 mg      -0.30    (AG   +0.5%)
##  ##     livoletide 60 ug/kg     -0.21    (AG:UAG -21.9%)
##  ##     octreotide LAR 30 mg    -0.50    (AG   -72.5%)
##  ##     semaglutide 2.4 mg QW   -2.78    (AG   +4.4%)
##  ##   The ghrelin columns and the HQ-CT column are nearly ORTHOGONAL.
##  arm <- function(dose_ev, wk)
##    mod |> param(c(LES = 1, ENV_MGMT)) |>
##      ev(c(gh, dose_ev)) |> mrgsim(end = YEARS(12) + wk*7, delta = 7)
##  cb32 <- ev(time = YEARS(12), amt = 3.2,  cmt = "ACBD", ii = 1/3, addl = 167)
##  cb96 <- ev(time = YEARS(12), amt = 9.6,  cmt = "ACBD", ii = 1/3, addl = 167)
##  dccr <- ev(time = YEARS(12), amt = 200,  cmt = "ADZD", ii = 1,   addl = 90)
##  setm <- ev(time = YEARS(12), amt = 3.0,  cmt = "ASMD", ii = 1,   addl = 83)
##  livo <- ev(time = YEARS(12), amt = 2.4,  cmt = "ALVD", ii = 1,   addl = 83)
##  octr <- ev(time = YEARS(12), amt = 30,   cmt = "AOCD", ii = 28,  addl = 2)
##  sema <- ev(time = YEARS(12), amt = 2.4,  cmt = "ASGD", ii = 7,   addl = 12)
##  metf <- ev(time = YEARS(12), amt = 1500, cmt = "AMFD", ii = 1,   addl = 90)
##
##  ## ---- S8  carbetocin dose-ranging: the INVERTED response --------
##  ##   the optimum is analytic, C* = 1.21 ng/mL; 3.2 mg TID gives
##  ##   Cavg 2.06 ng/mL (98% of E(C*)) and 9.6 mg gives 6.19 (62%).
##  ##   Scan 0.4 -> 19.2 mg TID and the response peaks and then FALLS.
##  ##
##  ## ---- S9  DCCR dose-ranging: a FIXED therapeutic index ----------
##  ##   efficacy and hyperglycaemia are the same K-ATP channel in two
##  ##   tissues, so dHQ-CT and dHbA1c rise together at a ratio set once
##  ##   by EDZ50/EDZI50 = 0.68.  No dose separates them.
##  ##
##  ## ---- S10 subtype: deletion vs mUPD -----------------------------
##  ##   param(FSUB = 1.15) for deletion (worse behaviour)
##  ##   param(FSUB = 0.90) for mUPD
##  ##
##  ## ---- S11 sex-steroid replacement from 14 y ---------------------
##  test <- ev(time = YEARS(14), amt = 100, cmt = "ATSD", ii = 14, addl = 120)
##  ##
##  ## ---- S12 long-acting weekly growth hormone --------------------
##  lagh <- ev(time = YEARS(1), amt = 4.9, cmt = "ALGD", ii = 7, addl = 1250)
##  ##
##  ## ---- S13 the LATCH: environment gates the pharmacology ---------
##  ##   Let a patient live 10-12 y with free access (ENV_FREE), then
##  ##   restore ENV_MGMT.  The fat comes off; the STATE does not
##  ##   (SEEK 0.71 -> 0.69 over 24 months).  Nothing in the panel
##  ##   annihilates the upper fixed point: that needs an oxytocin-arm
##  ##   gain of 2.4x and carbetocin's analytic ceiling is 1.49x.
##  ##   => the model's drug-development target is EFFICACY, not potency,
##  ##      and its trial-design prediction is a BIMODAL HQ-CT response.
##  ##
##  ## ---- S14 free access + everything -----------------------------
##  ##   In a free-access environment the low state does not exist, so
##  ##   carbetocin + DCCR together still leave SEEK at 0.63.  This is
##  ##   why DESTINY-PWS and CARE-PWS both required a stable food-secure
##  ##   environment to enrol, and why their effect sizes must not be
##  ##   extrapolated to families who do not have one.
##  ##
##  ## ---- S15 dt / solver check ------------------------------------
##  ##   mrgsim(..., hmax = 0.125) reproduces the Python RK4 run to
##  ##   < 1.1e-3 relative on every reported endpoint.
## =====================================================================
