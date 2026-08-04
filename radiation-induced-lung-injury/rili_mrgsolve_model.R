## =====================================================================
##  rili_mrgsolve_model.R
##  Radiation-Induced Lung Injury (RILI) — QSP model for mrgsolve
##  방사선 유발 폐손상 — mrgsolve QSP 모델
##
##  79 ODEs = 6 DVH dose bins x 10 states + 19 global
##
##  THE ONE STRUCTURAL IDEA
##  -----------------------
##  A radiotherapy plan is not a dose, it is a DISTRIBUTION.  The lung is
##  discretised into 6 dose-volume-histogram bins; bin b holds fractional
##  lung volume V1..V6 and receives total physical dose D_b = DRX * f_b.
##  All ten biological states exist per bin and every organ-level
##  read-out (DLCO, FVC, CTCAE grade, fibrosis score) is a volume-
##  weighted integral over bins.  Two plans with the same mean lung dose
##  can therefore differ in every other moment of the same histogram —
##  and they do (see scenario B below).
##
##  TWO LOOPS, TWO TIME CONSTANTS, ONE INSULT
##  -----------------------------------------
##  FAST  DAMP -> NF-kB -> cytokine -> oedema.  Weeks.  Sub-critical, so
##        it resolves on its own.  Steroid-sensitive.
##  SLOW  TGF-b1 -> myofibroblast -> collagen -> stiffness -> TGF-b1.
##        Months.  BISTABLE, so it latches.  Steroid-INSENSITIVE.
##  A drug that reaches one cannot be detected with an endpoint that
##  reads the other.  That is the model's central clinical claim.
##
##  VALIDATION
##  ----------
##  No R runtime was available where this was written, so every equation
##  below was first implemented and executed in dependency-free Python
##  (rili_reference_model.py) with a fixed-step RK4 integrator, and this
##  file is a line-by-line port of what that produced.  The Python run
##  exposed TWELVE real defects; each is documented at the point of
##  the fix in that file, and the ones that changed model structure (not
##  just parameter values) are flagged [DEFECT n] here too.
##
##  Fixed points of the slow subsystem, computed analytically and
##  reproduced by simulation:
##      healthy     TGFB 0.0401  MFB 0.00057  COL 1.0000   stable
##      separatrix  TGFB 0.2115  MFB 0.3563   COL 1.5864   UNSTABLE
##      fibrotic    TGFB 0.7767  MFB 0.9692   COL 2.9400   stable
##
##  USAGE
##      library(mrgsolve); library(dplyr)
##      mod <- mread("rili_mrgsolve_model", "path/to/dir")
##      out <- mod %>% mrgsim(end = 730, delta = 1) %>% as_tibble()
##  Scenarios are at the bottom of this file.
## =====================================================================

$PROB
# Radiation-Induced Lung Injury (RILI)
# 6 DVH bins x 10 states + 19 global = 79 ODEs

$PARAM @annotated
// ---- plan geometry --------------------------------------------------
DRX    : 60.0  : Prescription dose (Gy)
NFX    : 30    : Number of fractions
FXWK   :  5.0  : Fractions per week
V1     : 0.40  : Fractional lung volume in dose bin 1
V2     : 0.20  : Fractional lung volume in dose bin 2
V3     : 0.12  : Fractional lung volume in dose bin 3
V4     : 0.11  : Fractional lung volume in dose bin 4
V5     : 0.09  : Fractional lung volume in dose bin 5
V6     : 0.08  : Fractional lung volume in dose bin 6
// bin representative dose as a fraction of the prescription; these are
// the mid-points of bin edges 0, .08, .24, .44, .66, .88, 1.05
F1     : 0.040 : Bin 1 dose fraction of DRX
F2     : 0.160 : Bin 2 dose fraction of DRX
F3     : 0.340 : Bin 3 dose fraction of DRX
F4     : 0.550 : Bin 4 dose fraction of DRX
F5     : 0.770 : Bin 5 dose fraction of DRX
F6     : 0.965 : Bin 6 dose fraction of DRX
// regional perfusion weight; lower it where emphysema has destroyed
// tissue, so dose delivered there costs less gas exchange
PF1    : 1.0   : Perfusion weight bin 1
PF2    : 1.0   : Perfusion weight bin 2
PF3    : 1.0   : Perfusion weight bin 3
PF4    : 1.0   : Perfusion weight bin 4
PF5    : 1.0   : Perfusion weight bin 5
PF6    : 1.0   : Perfusion weight bin 6

// ---- radiobiology ---------------------------------------------------
ABL    :  3.0  : alpha/beta of late-responding lung (Gy)
ABT    : 10.0  : alpha/beta of NSCLC (Gy)
ALPHAT :  0.31 : Tumour alpha (1/Gy)

// ---- AT2 (type II pneumocyte) ---------------------------------------
KDAM   : 0.0250 : AT2 lethal hits per Gy3 of lung BED
KMIT   : 0.0250 : Mitosis-linked death of lethally hit AT2 (1/d), tau 40 d
KREP   : 0.0300 : Logistic AT2 repopulation rate (1/d)
KDCYT  : 0.0150 : Cytokine-mediated bystander AT2 loss (1/d)
KMCYT  : 1.0000 : Half-max cytokine for bystander loss

// ---- microvascular endothelium --------------------------------------
KECDAM : 0.0300 : Endothelial radiation kill per Gy3/d
KECREP : 0.0120 : Endothelial repair rate (1/d) - slow in irradiated lung
KECCYT : 0.0300 : Cytokine-mediated endothelial loss
KECIRR : 0.0043 : Permanent loss of repair ceiling per Gy3 cumulative BED

// ---- surfactant ------------------------------------------------------
KSYN   : 0.5000 : Surfactant synthesis per unit (AT2 + DOOM)
KSDEG  : 0.5000 : Surfactant turnover

// ---- alveolar-capillary leak / oedema -------------------------------
KPRMEC : 0.0250 : Leak from endothelial loss
KPRMEP : 0.0900 : Leak from surfactant/epithelial loss (dominant arm)
KPCYT  : 0.1000 : Leak from cytokine tone
KPCL   : 0.0450 : Lymphatic clearance of alveolar oedema (1/d), tau 22 d

// ---- FAST loop: DAMP -> NF-kB -> cytokine ---------------------------
KDAMP  : 2.0000 : DAMP yield per unit dying-cell flux
KDIR   : 0.0050 : Direct NF-kB activation per Gy3/d during delivery
KAMP   : 0.0250 : Autocrine cytokine amplification
CMAX   : 2.0000 : Saturation of the amplification term
KCE    : 0.1500 : Cytokine elimination (1/d)
KSYSIN : 0.1000 : Systemic -> local cytokine spill-in
KSOUT  : 0.0500 : Local -> systemic cytokine spill-out
KCLS   : 0.3000 : Systemic cytokine clearance

// ---- SLOW loop: TGF-beta1 -------------------------------------------
KTGF0  : 0.0100 : Constitutive latent-TGFb activation
KTGFR  : 0.0120 : Radiation/ROS-driven activation per Gy3/d
KTGFI  : 0.0500 : Inflammation-driven activation (consequential effect)
KTGFM  : 0.1600 : Mechanotransduction-driven activation
KMCOL  : 1.0000 : Collagen excess at half-max force transmission
HTGFM  : 3.0000 : Hill coefficient of the stiffness threshold
KTGFA  : 0.0450 : Myofibroblast autocrine TGFb
KTGFEL : 0.2500 : Active TGFb elimination

// ---- myofibroblast ---------------------------------------------------
KMFB   : 1.6000 : Myofibroblast differentiation rate
KMTGF  : 0.5500 : Half-max TGFb for differentiation
HMFB   : 4.0000 : Hill coefficient; > 2 is what buys the bistability
KMFBS  : 0.5000 : Stiffness boost to differentiation
KMFBD  : 0.0800 : Myofibroblast apoptosis (1/d)

// ---- collagen and cross-linking -------------------------------------
KCOL0  : 0.0477625 : Constitutive synthesis, set so COL* = 1.0000 exactly
KCOL   : 0.0500 : Myofibroblast-driven collagen synthesis
KCDEG  : 0.0500 : MMP-mediated degradation
KTIMP  : 1.2000 : TGFb-driven TIMP inhibition of degradation
KXL    : 0.0080 : LOX maturation of labile -> cross-linked collagen
XMAX   : 1.0000 : Maximum cross-linkable collagen density
KXDEG  : 0.0015 : Cross-linked collagen turnover (1/d), tau ~ 1.8 y

// ---- function mapping ------------------------------------------------
PEC    : 0.7000 : Exponent of endothelial integrity in diffusing capacity
PSURF  : 0.3000 : Exponent of surfactant in diffusing capacity
KPF    : 0.5000 : Oedema half-effect on diffusing capacity
KCF    : 1.5000 : Collagen half-effect on diffusing capacity
KCV    : 2.5000 : Collagen half-effect on vital capacity
KPV    : 1.2000 : Oedema half-effect on vital capacity
DLCO0  : 85.0   : Baseline DLCO (% predicted)
FVC0   : 95.0   : Baseline FVC (% predicted)

// ---- clinical scoring ------------------------------------------------
A1     : 1.0000 : Weight of oedema in the pneumonitis index
A2     : 0.6000 : Weight of cytokine tone in the pneumonitis index
DLREF  : 85.0   : Reference baseline DLCO for reserve scaling
PNI50  : 1.2900 : Index giving 50% risk of grade >=2 RP (fitted, QUANTEC)
PNISL  : 0.4100 : Logistic slope on ln(index)

// ---- drug PK / PD ----------------------------------------------------
KADEX  : 6.0000 : Prednisolone absorption (1/d)
CLDEX  : 8.0000 : Prednisolone clearance (L/d)
VDEX   : 50.000 : Prednisolone volume (L)
IMXDEX : 0.8500 : Max NF-kB suppression by steroid
IC50DEX: 0.0300 : Steroid IC50 (mg/L)
EMXCLR : 1.2000 : Max fold-increase in alveolar fluid clearance by steroid
KAPIR  : 4.0000 : Pirfenidone absorption (1/d)
CLPIR  : 60.000 : Pirfenidone clearance (L/d)
VPIR   : 70.000 : Pirfenidone volume (L)
IMXPIR : 0.6000 : Max TGFb-signalling / collagen-synthesis inhibition
IC50PIR: 3.0000 : Pirfenidone IC50 (mg/L)
KANIN  : 1.5000 : Nintedanib absorption (1/d)
CLNIN  : 1390.0 : Nintedanib clearance (L/d)
VNIN   : 1050.0 : Nintedanib volume (L)
IMXNIN : 0.6500 : Max fibroblast-proliferation inhibition
IC50NIN: 0.0150 : Nintedanib IC50 (mg/L)
CLDUR  : 0.2320 : Durvalumab clearance (L/d)
VDUR   : 3.6000 : Durvalumab central volume (L)
QDUR   : 0.4800 : Durvalumab intercompartmental clearance (L/d)
VPDUR  : 3.0000 : Durvalumab peripheral volume (L)
EMXDUR : 0.3500 : Max fast-loop gain increase by checkpoint blockade
EC50DUR: 8.0000 : Durvalumab EC50 (mg/L)
KTUMIMM: 0.0060 : Immune-mediated tumour kill (1/d)
KAACE  : 1.0000 : Lisinopril absorption (1/d)
CLACE  : 6.0000 : Lisinopril clearance (L/d)
VACE   : 100.00 : Lisinopril volume (L)
IMXACE : 0.4500 : Max inhibition of latent-TGFb activation
IC50ACE: 0.0200 : Lisinopril IC50 (mg/L)

// ---- radioprotector: handled algebraically, not as an ODE ------------
// A thiol radioprotector only protects tissue irradiated WHILE the drug
// is present.  WR-1065 has a plasma half-life of ~8 min, so protection
// decays exponentially with the drug-to-beam interval.  That one
// parameter separates the amifostine trials that found benefit from
// those that did not, so it is exposed directly.
AMION  : 0     : Amifostine given (0/1)
AMIDEL : 20.0  : Minutes between amifostine and the beam
PFAMI  : 0.6500: Maximal amifostine radioprotection factor
THAMI  : 8.0   : WR-1065 plasma half-life (min)
AVAON  : 0     : Avasopasem given (0/1)
AVADEL : 30.0  : Minutes between avasopasem and the beam
PFAVA  : 0.4500: Maximal avasopasem radioprotection factor
THAVA  : 27.0  : Avasopasem plasma half-life (min)

// ---- oral dosing rates (mg/day) and windows (day) --------------------
// Oral drugs enter as a mg/day rate.  This is exact for the daily
// averaged exposure and adequate because every PD time constant
// downstream is >= days.  Durvalumab is a true q28d IV bolus and is
// given through $PKMODEL-style dosing records into DURC instead.
RDEX   : 0.0   : Prednisolone rate (mg/d)
TDEX0  : 1e9   : Prednisolone start day
TDEXD  : 14.0  : Days at full steroid dose before taper
TDEXT  : 42.0  : Days of linear steroid taper
RPIR   : 0.0   : Pirfenidone rate (mg/d)
TPIR0  : 1e9   : Pirfenidone start day
TPIR1  : 1e9   : Pirfenidone stop day
RNIN   : 0.0   : Nintedanib rate (mg/d)
TNIN0  : 1e9   : Nintedanib start day
TNIN1  : 1e9   : Nintedanib stop day
RACE   : 0.0   : Lisinopril rate (mg/d)
TACE0  : 1e9   : Lisinopril start day
TACE1  : 1e9   : Lisinopril stop day

// ---- tumour ----------------------------------------------------------
KGROW  : 0.0347 : Tumour clonogen repopulation rate (1/d)
TUMLN0 : 20.72  : ln(clonogens) at start, ln(1e9)

$CMT @annotated
// bin 1
AT2_1 : Viable type II pneumocytes, bin 1
DOM_1 : Lethally hit AT2 awaiting mitosis, bin 1
EC_1  : Microvascular endothelial integrity, bin 1
SUR_1 : Surfactant pool, bin 1
PRM_1 : Alveolar-capillary leak / oedema, bin 1
CYT_1 : Local pro-inflammatory cytokine pool, bin 1
TGF_1 : Active TGF-beta1, bin 1
MFB_1 : Myofibroblast density, bin 1
COL_1 : Labile collagen, bin 1
XCL_1 : Cross-linked collagen, bin 1
// bin 2
AT2_2 : Viable type II pneumocytes, bin 2
DOM_2 : Lethally hit AT2 awaiting mitosis, bin 2
EC_2  : Microvascular endothelial integrity, bin 2
SUR_2 : Surfactant pool, bin 2
PRM_2 : Alveolar-capillary leak / oedema, bin 2
CYT_2 : Local pro-inflammatory cytokine pool, bin 2
TGF_2 : Active TGF-beta1, bin 2
MFB_2 : Myofibroblast density, bin 2
COL_2 : Labile collagen, bin 2
XCL_2 : Cross-linked collagen, bin 2
// bin 3
AT2_3 : Viable type II pneumocytes, bin 3
DOM_3 : Lethally hit AT2 awaiting mitosis, bin 3
EC_3  : Microvascular endothelial integrity, bin 3
SUR_3 : Surfactant pool, bin 3
PRM_3 : Alveolar-capillary leak / oedema, bin 3
CYT_3 : Local pro-inflammatory cytokine pool, bin 3
TGF_3 : Active TGF-beta1, bin 3
MFB_3 : Myofibroblast density, bin 3
COL_3 : Labile collagen, bin 3
XCL_3 : Cross-linked collagen, bin 3
// bin 4
AT2_4 : Viable type II pneumocytes, bin 4
DOM_4 : Lethally hit AT2 awaiting mitosis, bin 4
EC_4  : Microvascular endothelial integrity, bin 4
SUR_4 : Surfactant pool, bin 4
PRM_4 : Alveolar-capillary leak / oedema, bin 4
CYT_4 : Local pro-inflammatory cytokine pool, bin 4
TGF_4 : Active TGF-beta1, bin 4
MFB_4 : Myofibroblast density, bin 4
COL_4 : Labile collagen, bin 4
XCL_4 : Cross-linked collagen, bin 4
// bin 5
AT2_5 : Viable type II pneumocytes, bin 5
DOM_5 : Lethally hit AT2 awaiting mitosis, bin 5
EC_5  : Microvascular endothelial integrity, bin 5
SUR_5 : Surfactant pool, bin 5
PRM_5 : Alveolar-capillary leak / oedema, bin 5
CYT_5 : Local pro-inflammatory cytokine pool, bin 5
TGF_5 : Active TGF-beta1, bin 5
MFB_5 : Myofibroblast density, bin 5
COL_5 : Labile collagen, bin 5
XCL_5 : Cross-linked collagen, bin 5
// bin 6
AT2_6 : Viable type II pneumocytes, bin 6
DOM_6 : Lethally hit AT2 awaiting mitosis, bin 6
EC_6  : Microvascular endothelial integrity, bin 6
SUR_6 : Surfactant pool, bin 6
PRM_6 : Alveolar-capillary leak / oedema, bin 6
CYT_6 : Local pro-inflammatory cytokine pool, bin 6
TGF_6 : Active TGF-beta1, bin 6
MFB_6 : Myofibroblast density, bin 6
COL_6 : Labile collagen, bin 6
XCL_6 : Cross-linked collagen, bin 6
// global
CYTS  : Systemic cytokine pool
BEL_1 : Cumulative lung BED3, bin 1
BEL_2 : Cumulative lung BED3, bin 2
BEL_3 : Cumulative lung BED3, bin 3
BEL_4 : Cumulative lung BED3, bin 4
BEL_5 : Cumulative lung BED3, bin 5
BEL_6 : Cumulative lung BED3, bin 6
BEDT  : Cumulative tumour BED10
TUMLN : ln(surviving tumour clonogens)
DEXD  : Prednisolone depot (mg)
DEXC  : Prednisolone central (mg)
PIRD  : Pirfenidone depot (mg)
PIRC  : Pirfenidone central (mg)
NIND  : Nintedanib depot (mg)
NINC  : Nintedanib central (mg)
DURC  : Durvalumab central (mg)
DURP  : Durvalumab peripheral (mg)
ACED  : Lisinopril depot (mg)
ACEC  : Lisinopril central (mg)

$MAIN
// the healthy lung is a genuine fixed point of this system; the
// analytic values are used as initial conditions so that a beam-off
// control run returns exactly baseline (verified: COL drifts by 1.2e-5
// over 730 days)
AT2_1_0 = 1.0; DOM_1_0 = 0.0; EC_1_0 = 1.0; SUR_1_0 = 1.0; PRM_1_0 = 0.0;
CYT_1_0 = 0.0; TGF_1_0 = 0.040102; MFB_1_0 = 0.000565; COL_1_0 = 1.0;
XCL_1_0 = 0.0;
AT2_2_0 = 1.0; DOM_2_0 = 0.0; EC_2_0 = 1.0; SUR_2_0 = 1.0; PRM_2_0 = 0.0;
CYT_2_0 = 0.0; TGF_2_0 = 0.040102; MFB_2_0 = 0.000565; COL_2_0 = 1.0;
XCL_2_0 = 0.0;
AT2_3_0 = 1.0; DOM_3_0 = 0.0; EC_3_0 = 1.0; SUR_3_0 = 1.0; PRM_3_0 = 0.0;
CYT_3_0 = 0.0; TGF_3_0 = 0.040102; MFB_3_0 = 0.000565; COL_3_0 = 1.0;
XCL_3_0 = 0.0;
AT2_4_0 = 1.0; DOM_4_0 = 0.0; EC_4_0 = 1.0; SUR_4_0 = 1.0; PRM_4_0 = 0.0;
CYT_4_0 = 0.0; TGF_4_0 = 0.040102; MFB_4_0 = 0.000565; COL_4_0 = 1.0;
XCL_4_0 = 0.0;
AT2_5_0 = 1.0; DOM_5_0 = 0.0; EC_5_0 = 1.0; SUR_5_0 = 1.0; PRM_5_0 = 0.0;
CYT_5_0 = 0.0; TGF_5_0 = 0.040102; MFB_5_0 = 0.000565; COL_5_0 = 1.0;
XCL_5_0 = 0.0;
AT2_6_0 = 1.0; DOM_6_0 = 0.0; EC_6_0 = 1.0; SUR_6_0 = 1.0; PRM_6_0 = 0.0;
CYT_6_0 = 0.0; TGF_6_0 = 0.040102; MFB_6_0 = 0.000565; COL_6_0 = 1.0;
XCL_6_0 = 0.0;
TUMLN_0 = TUMLN0;

$GLOBAL
#define NBIN 6

// per-bin geometry, filled in $PREAMBLE-equivalent code each call
double gV[NBIN], gD[NBIN], gRB[NBIN], gW[NBIN];
double gTCOURSE, gRBT, gBEDT;

// Vx: fractional lung volume receiving >= x Gy, interpolated inside the
// bin that contains x.  Bin edges as fractions of DRX.
double edgef(int i) {
  static const double e[7] = {0.0, 0.08, 0.24, 0.44, 0.66, 0.88, 1.05};
  return e[i];
}

$PREAMBLE
// nothing time-varying here; geometry is rebuilt in $ODE from parameters
// so that DRX/NFX/V1..V6 can be swept from an idata set

$ODE
// ------------------------------------------------------------------ //
// 0.  plan geometry                                                  //
// ------------------------------------------------------------------ //
double VV[NBIN]; VV[0]=V1; VV[1]=V2; VV[2]=V3; VV[3]=V4; VV[4]=V5; VV[5]=V6;
double FF[NBIN]; FF[0]=F1; FF[1]=F2; FF[2]=F3; FF[3]=F4; FF[4]=F5; FF[5]=F6;
double PP[NBIN]; PP[0]=PF1; PP[1]=PF2; PP[2]=PF3; PP[3]=PF4; PP[4]=PF5;
PP[5]=PF6;

double vsum = 0.0, wsum = 0.0;
int b;
for (b = 0; b < NBIN; ++b) vsum += VV[b];
if (vsum <= 0) vsum = 1.0;

double TCOURSE = NFX / FXWK * 7.0;
double RTON = (SOLVERTIME >= 0.0 && SOLVERTIME <= TCOURSE) ? 1.0 : 0.0;

// radioprotection factor: TEMPORAL COINCIDENCE, not concentration
double PROT = 0.0;
if (AMION > 0.5) PROT += PFAMI * exp(-log(2.0) * AMIDEL / THAMI);
if (AVAON > 0.5) PROT += PFAVA * exp(-log(2.0) * AVADEL / THAVA);
if (PROT > 0.95) PROT = 0.95;

double DB[NBIN], DFX[NBIN], BEDL[NBIN], RB[NBIN], WT[NBIN];
for (b = 0; b < NBIN; ++b) {
  DB[b]   = DRX * FF[b];
  DFX[b]  = DB[b] / NFX;
  // linear-quadratic: BED = D (1 + d/(alpha/beta)).  This is the only
  // place fractionation enters the normal-tissue side, and it enters
  // through the LUNG alpha/beta of 3 Gy while the tumour sees 10 Gy.
  BEDL[b] = DB[b] * (1.0 + DFX[b] / ABL);
  RB[b]   = BEDL[b] / TCOURSE * RTON * (1.0 - PROT);
  WT[b]   = VV[b] / vsum * PP[b];
  wsum   += WT[b];
}
if (wsum <= 0) wsum = 1.0;

double dT = DRX / NFX;
double BEDTt = DRX * (1.0 + dT / ABT);
double RBT = BEDTt / TCOURSE * RTON;

// ------------------------------------------------------------------ //
// 1.  drug PK                                                        //
// ------------------------------------------------------------------ //
// oral inputs as mg/day rates within their windows; the steroid has a
// full-dose period followed by a linear taper
double rdex = 0.0;
if (RDEX > 0.0 && SOLVERTIME >= TDEX0) {
  double el = SOLVERTIME - TDEX0;
  if (el <= TDEXD)                 rdex = RDEX;
  else if (el <= TDEXD + TDEXT)    rdex = RDEX * (1.0 - (el - TDEXD) / TDEXT);
}
double rpir = (SOLVERTIME >= TPIR0 && SOLVERTIME <= TPIR1) ? RPIR : 0.0;
double rnin = (SOLVERTIME >= TNIN0 && SOLVERTIME <= TNIN1) ? RNIN : 0.0;
double race = (SOLVERTIME >= TACE0 && SOLVERTIME <= TACE1) ? RACE : 0.0;

dxdt_DEXD = rdex - KADEX * DEXD;
dxdt_DEXC = KADEX * DEXD - CLDEX / VDEX * DEXC;
dxdt_PIRD = rpir - KAPIR * PIRD;
dxdt_PIRC = KAPIR * PIRD - CLPIR / VPIR * PIRC;
dxdt_NIND = rnin - KANIN * NIND;
dxdt_NINC = KANIN * NIND - CLNIN / VNIN * NINC;
dxdt_ACED = race - KAACE * ACED;
dxdt_ACEC = KAACE * ACED - CLACE / VACE * ACEC;
dxdt_DURC = -CLDUR / VDUR * DURC - QDUR / VDUR * DURC + QDUR / VPDUR * DURP;
dxdt_DURP =  QDUR / VDUR * DURC - QDUR / VPDUR * DURP;

double CDEX = DEXC / VDEX;
double CPIR = PIRC / VPIR;
double CNIN = NINC / VNIN;
double CACE = ACEC / VACE;
double CDUR = DURC / VDUR;

// ------------------------------------------------------------------ //
// 2.  drug effects, each labelled by WHICH LOOP it reaches           //
// ------------------------------------------------------------------ //
// steroid, route 1: transcriptional suppression of the FAST loop
double SDEX = 1.0 - IMXDEX * CDEX / (IC50DEX + CDEX);
// steroid, route 2: ENaC / Na,K-ATPase driven alveolar fluid clearance.
// [DEFECT 12] without this second route the steroid reached the oedema
// state only through the cytokine arm, and since oedema is carried
// mostly by the epithelial arm the entire benefit was a 16% fall in the
// index -- far less than the observed clinical response.
double ECLR = EMXCLR * CDEX / (IC50DEX + CDEX);
// pirfenidone: TGFb signalling + collagen synthesis.  SLOW loop only.
double IPIR = IMXPIR * CPIR / (IC50PIR + CPIR);
// nintedanib: PDGFR/FGFR/VEGFR, fibroblast proliferation.  SLOW only.
double ININ = IMXNIN * CNIN / (IC50NIN + CNIN);
// ACE inhibitor: latent TGFb activation.  SLOW loop only.
double IACE = IMXACE * CACE / (IC50ACE + CACE);
// durvalumab: raises the GAIN of the fast loop.  Not an added term --
// a multiplier -- which is why the excess risk grows with mean lung dose
// instead of being a constant offset.
double GIMM = 1.0 + EMXDUR * CDUR / (EC50DUR + CDUR);

// ------------------------------------------------------------------ //
// 3.  per-bin biology                                                //
// ------------------------------------------------------------------ //
double AT2v[NBIN], DOMv[NBIN], ECv[NBIN], SURv[NBIN], PRMv[NBIN];
double CYTv[NBIN], TGFv[NBIN], MFBv[NBIN], COLv[NBIN], XCLv[NBIN];
double BELv[NBIN];
AT2v[0]=AT2_1; DOMv[0]=DOM_1; ECv[0]=EC_1; SURv[0]=SUR_1; PRMv[0]=PRM_1;
CYTv[0]=CYT_1; TGFv[0]=TGF_1; MFBv[0]=MFB_1; COLv[0]=COL_1; XCLv[0]=XCL_1;
BELv[0]=BEL_1;
AT2v[1]=AT2_2; DOMv[1]=DOM_2; ECv[1]=EC_2; SURv[1]=SUR_2; PRMv[1]=PRM_2;
CYTv[1]=CYT_2; TGFv[1]=TGF_2; MFBv[1]=MFB_2; COLv[1]=COL_2; XCLv[1]=XCL_2;
BELv[1]=BEL_2;
AT2v[2]=AT2_3; DOMv[2]=DOM_3; ECv[2]=EC_3; SURv[2]=SUR_3; PRMv[2]=PRM_3;
CYTv[2]=CYT_3; TGFv[2]=TGF_3; MFBv[2]=MFB_3; COLv[2]=COL_3; XCLv[2]=XCL_3;
BELv[2]=BEL_3;
AT2v[3]=AT2_4; DOMv[3]=DOM_4; ECv[3]=EC_4; SURv[3]=SUR_4; PRMv[3]=PRM_4;
CYTv[3]=CYT_4; TGFv[3]=TGF_4; MFBv[3]=MFB_4; COLv[3]=COL_4; XCLv[3]=XCL_4;
BELv[3]=BEL_4;
AT2v[4]=AT2_5; DOMv[4]=DOM_5; ECv[4]=EC_5; SURv[4]=SUR_5; PRMv[4]=PRM_5;
CYTv[4]=CYT_5; TGFv[4]=TGF_5; MFBv[4]=MFB_5; COLv[4]=COL_5; XCLv[4]=XCL_5;
BELv[4]=BEL_5;
AT2v[5]=AT2_6; DOMv[5]=DOM_6; ECv[5]=EC_6; SURv[5]=SUR_6; PRMv[5]=PRM_6;
CYTv[5]=CYT_6; TGFv[5]=TGF_6; MFBv[5]=MFB_6; COLv[5]=COL_6; XCLv[5]=XCL_6;
BELv[5]=BEL_6;

double dAT2[NBIN], dDOM[NBIN], dEC[NBIN], dSUR[NBIN], dPRM[NBIN];
double dCYT[NBIN], dTGF[NBIN], dMFB[NBIN], dCOL[NBIN], dXCL[NBIN];
double spill = 0.0;

for (b = 0; b < NBIN; ++b) {
  double AT2 = AT2v[b], DOM = DOMv[b], EC = ECv[b], SUR = SURv[b];
  double PRM = PRMv[b], CYT = CYTv[b], TGF = TGFv[b], MFB = MFBv[b];
  double COL = COLv[b], XCL = XCLv[b];
  double COLX = (COL - 1.0 > 0.0) ? (COL - 1.0) : 0.0;
  // stiffness, mechanotransduction and every functional read-out see
  // labile PLUS cross-linked collagen
  double CXT = COLX + XCL;

  // ---- AT2.  Radiation moves cells into a DOOMED pool that only dies
  // when it next attempts mitosis.  This delay is the reason
  // late-responding tissue injury is delayed at all.
  double JRAD = KDAM * RB[b] * AT2;
  double JCYT = KDCYT * (CYT / (KMCYT + CYT)) * AT2;
  double cap  = 1.0 - (AT2 + DOM);
  if (cap < 0.0) cap = 0.0;
  double grow = KREP * AT2 * cap;
  dAT2[b] = grow - JRAD - JCYT;
  dDOM[b] = JRAD + JCYT - KMIT * DOM;
  double DTH = KMIT * DOM;          // dying-cell flux -> the DAMP source

  // ---- microvascular endothelium.  Repair is toward a CEILING that
  // falls with cumulative dose: capillary dropout and loss of
  // endothelial progenitors are not repaired, so an irradiated bin can
  // never return to EC = 1.  This is the permanent, dose-graded
  // component of DLCO loss, and it exists in bins that never fibrose.
  double ECMAX = 1.0 / (1.0 + KECIRR * BELv[b]);
  dEC[b] = KECREP * (ECMAX - EC) - KECDAM * RB[b] * EC - KECCYT * CYT * EC;

  // ---- surfactant.  [DEFECT 9] the producing pool is AT2 + DOOM, not
  // viable AT2 alone: a lethally hit pneumocyte keeps making surfactant
  // until it dies, so the alveolar lining keeps thinning for weeks after
  // the last fraction.  With viable AT2 only, the whole injury peaked at
  // the end of RT instead of 4-12 weeks later.
  dSUR[b] = KSYN * (AT2 + DOM) - KSDEG * SUR;

  // ---- alveolar-capillary leak.  [DEFECT 10] the endothelial and
  // epithelial sources are SEPARATE coefficients: endothelium recovers
  // on KECREP while the epithelial pool keeps falling to day ~84, and
  // weighting them alike made the two arms cancel.
  dPRM[b] = KPRMEC * (1.0 - EC) + KPRMEP * (1.0 - SUR) + KPCYT * CYT
            - KPCL * (1.0 + ECLR) * PRM;

  // ---- FAST loop.  Sub-critical by construction:
  //   GIMM*(KAMP + KDAMP*KDCYT/KMCYT) + KSYSIN*KSOUT/KCLS  <  KCE
  //   no drug     0.0717 < 0.150      durvalumab  0.0910 < 0.150
  // so pneumonitis resolves on its own, and checkpoint blockade
  // lengthens the effective time constant from 12.8 d to 17.0 d rather
  // than adding a source term.
  double amp = KAMP * CYT / (1.0 + CYT / CMAX);
  dCYT[b] = SDEX * GIMM * (KDAMP * DTH + KDIR * RB[b] + amp)
            + KSYSIN * CYTS - KCE * CYT;
  spill += VV[b] / vsum * CYT;

  // ---- SLOW loop.  Four activators of latent TGF-beta1; only the
  // mechanotransduction and autocrine limbs close a positive feedback.
  // [DEFECT 2] the stiffness limb is a Hill-3 THRESHOLD, not a ratio:
  // latent TGF-b1 is unfolded by integrin-transmitted force only once
  // the matrix is stiff enough to resist myofibroblast contraction.
  // Written as km*CX/(1+CX) it was near-linear at low collagen, which
  // put enough gain around the loop to destroy the healthy fixed point
  // entirely -- the unirradiated lung had nowhere to sit.
  double mech = KTGFM * pow(CXT, HTGFM) / (pow(KMCOL, HTGFM)
                                           + pow(CXT, HTGFM));
  double act  = KTGF0 + KTGFR * RB[b] + KTGFI * CYT + mech + KTGFA * MFB;
  dTGF[b] = act * (1.0 - IACE) - KTGFEL * TGF;

  double drive = pow(TGF, HMFB) / (pow(KMTGF, HMFB) + pow(TGF, HMFB));
  double room  = 1.0 - MFB;
  if (room < 0.0) room = 0.0;
  dMFB[b] = KMFB * drive * (1.0 - IPIR) * (1.0 + KMFBS * CXT)
            * (1.0 - ININ) * room - KMFBD * MFB;

  double deg = KCDEG * COL / (1.0 + KTIMP * TGF / (1.0 + TGF));
  // [DEFECT 11] LOX/LOXL2 cross-linked collagen is not an MMP substrate.
  // Without this pool an antifibrotic started a YEAR after the switch
  // flipped abolished the fibrosis completely, i.e. the model claimed
  // pirfenidone reverses established radiation fibrosis.  Mature matrix
  // is mechanical MEMORY, and cross-linking must also be SLOW (0.008/d)
  // so the labile loop alone decides whether the switch flips while the
  // cross-linked pool accumulates afterwards and locks it in.
  double xroom = XMAX - XCL;
  if (xroom < 0.0) xroom = 0.0;
  double xflux = KXL * COLX * xroom;
  dCOL[b] = KCOL0 + KCOL * MFB * (1.0 - IPIR) - deg - xflux;
  dXCL[b] = xflux - KXDEG * XCL;
}

dxdt_AT2_1=dAT2[0]; dxdt_DOM_1=dDOM[0]; dxdt_EC_1=dEC[0];
dxdt_SUR_1=dSUR[0]; dxdt_PRM_1=dPRM[0]; dxdt_CYT_1=dCYT[0];
dxdt_TGF_1=dTGF[0]; dxdt_MFB_1=dMFB[0]; dxdt_COL_1=dCOL[0];
dxdt_XCL_1=dXCL[0];
dxdt_AT2_2=dAT2[1]; dxdt_DOM_2=dDOM[1]; dxdt_EC_2=dEC[1];
dxdt_SUR_2=dSUR[1]; dxdt_PRM_2=dPRM[1]; dxdt_CYT_2=dCYT[1];
dxdt_TGF_2=dTGF[1]; dxdt_MFB_2=dMFB[1]; dxdt_COL_2=dCOL[1];
dxdt_XCL_2=dXCL[1];
dxdt_AT2_3=dAT2[2]; dxdt_DOM_3=dDOM[2]; dxdt_EC_3=dEC[2];
dxdt_SUR_3=dSUR[2]; dxdt_PRM_3=dPRM[2]; dxdt_CYT_3=dCYT[2];
dxdt_TGF_3=dTGF[2]; dxdt_MFB_3=dMFB[2]; dxdt_COL_3=dCOL[2];
dxdt_XCL_3=dXCL[2];
dxdt_AT2_4=dAT2[3]; dxdt_DOM_4=dDOM[3]; dxdt_EC_4=dEC[3];
dxdt_SUR_4=dSUR[3]; dxdt_PRM_4=dPRM[3]; dxdt_CYT_4=dCYT[3];
dxdt_TGF_4=dTGF[3]; dxdt_MFB_4=dMFB[3]; dxdt_COL_4=dCOL[3];
dxdt_XCL_4=dXCL[3];
dxdt_AT2_5=dAT2[4]; dxdt_DOM_5=dDOM[4]; dxdt_EC_5=dEC[4];
dxdt_SUR_5=dSUR[4]; dxdt_PRM_5=dPRM[4]; dxdt_CYT_5=dCYT[4];
dxdt_TGF_5=dTGF[4]; dxdt_MFB_5=dMFB[4]; dxdt_COL_5=dCOL[4];
dxdt_XCL_5=dXCL[4];
dxdt_AT2_6=dAT2[5]; dxdt_DOM_6=dDOM[5]; dxdt_EC_6=dEC[5];
dxdt_SUR_6=dSUR[5]; dxdt_PRM_6=dPRM[5]; dxdt_CYT_6=dCYT[5];
dxdt_TGF_6=dTGF[5]; dxdt_MFB_6=dMFB[5]; dxdt_COL_6=dCOL[5];
dxdt_XCL_6=dXCL[5];

// systemic cytokine pool: this is the only route by which an
// unirradiated bin becomes inflamed, i.e. out-of-field pneumonitis
dxdt_CYTS = KSOUT * spill - KCLS * CYTS;

dxdt_BEL_1 = BEDL[0] / TCOURSE * RTON;
dxdt_BEL_2 = BEDL[1] / TCOURSE * RTON;
dxdt_BEL_3 = BEDL[2] / TCOURSE * RTON;
dxdt_BEL_4 = BEDL[3] / TCOURSE * RTON;
dxdt_BEL_5 = BEDL[4] / TCOURSE * RTON;
dxdt_BEL_6 = BEDL[5] / TCOURSE * RTON;
dxdt_BEDT  = RBT;

// ------------------------------------------------------------------ //
// 4.  tumour                                                         //
// ------------------------------------------------------------------ //
// [DEFECT 8] repopulation is gated off below one surviving clonogen.
// Applied unconditionally, a plan that drove the clonogen number to
// 0.85 cells still "regrew" to e^23 by day 730 and every TCP read 0%.
double alive = (TUMLN > 0.0) ? 1.0 : 0.0;
dxdt_TUMLN = KGROW * alive - ALPHAT * RBT
             - KTUMIMM * CDUR / (EC50DUR + CDUR);

$TABLE
double VVt[NBIN]; VVt[0]=V1; VVt[1]=V2; VVt[2]=V3; VVt[3]=V4; VVt[4]=V5;
VVt[5]=V6;
double FFt[NBIN]; FFt[0]=F1; FFt[1]=F2; FFt[2]=F3; FFt[3]=F4; FFt[4]=F5;
FFt[5]=F6;
double PPt[NBIN]; PPt[0]=PF1; PPt[1]=PF2; PPt[2]=PF3; PPt[3]=PF4;
PPt[4]=PF5; PPt[5]=PF6;
double ECt[NBIN]; ECt[0]=EC_1; ECt[1]=EC_2; ECt[2]=EC_3; ECt[3]=EC_4;
ECt[4]=EC_5; ECt[5]=EC_6;
double SUt[NBIN]; SUt[0]=SUR_1; SUt[1]=SUR_2; SUt[2]=SUR_3; SUt[3]=SUR_4;
SUt[4]=SUR_5; SUt[5]=SUR_6;
double PRt[NBIN]; PRt[0]=PRM_1; PRt[1]=PRM_2; PRt[2]=PRM_3; PRt[3]=PRM_4;
PRt[4]=PRM_5; PRt[5]=PRM_6;
double CYt[NBIN]; CYt[0]=CYT_1; CYt[1]=CYT_2; CYt[2]=CYT_3; CYt[3]=CYT_4;
CYt[4]=CYT_5; CYt[5]=CYT_6;
double COt[NBIN]; COt[0]=COL_1; COt[1]=COL_2; COt[2]=COL_3; COt[3]=COL_4;
COt[4]=COL_5; COt[5]=COL_6;
double XCt[NBIN]; XCt[0]=XCL_1; XCt[1]=XCL_2; XCt[2]=XCL_3; XCt[3]=XCL_4;
XCt[4]=XCL_5; XCt[5]=XCL_6;

int j;
double vs = 0.0, ws = 0.0;
for (j = 0; j < NBIN; ++j) vs += VVt[j];
if (vs <= 0) vs = 1.0;
double Wt[NBIN];
for (j = 0; j < NBIN; ++j) { Wt[j] = VVt[j] / vs * PPt[j]; ws += Wt[j]; }
if (ws <= 0) ws = 1.0;

// ---- plan descriptors ----------------------------------------------
double MLD = 0.0;
for (j = 0; j < NBIN; ++j) MLD += VVt[j] / vs * DRX * FFt[j];

// Vx by linear interpolation inside the bin containing x
double V5G = 0.0, V20G = 0.0, V40G = 0.0;
for (j = 0; j < NBIN; ++j) {
  double lo = DRX * edgef(j), hi = DRX * edgef(j + 1);
  double frac = VVt[j] / vs;
  if (hi > 5.0)  V5G  += (lo >= 5.0)  ? frac : frac * (hi - 5.0)  / (hi - lo);
  if (hi > 20.0) V20G += (lo >= 20.0) ? frac : frac * (hi - 20.0) / (hi - lo);
  if (hi > 40.0) V40G += (lo >= 40.0) ? frac : frac * (hi - 40.0) / (hi - lo);
}

// ---- function ------------------------------------------------------
double DL = 0.0, FV = 0.0, PNI = 0.0, FIB = 0.0, VFIB = 0.0;
for (j = 0; j < NBIN; ++j) {
  double ec = (ECt[j] > 1e-9) ? ECt[j] : 1e-9;
  double su = (SUt[j] > 1e-9) ? SUt[j] : 1e-9;
  double pr = (PRt[j] > 0.0) ? PRt[j] : 0.0;
  double cy = (CYt[j] > 0.0) ? CYt[j] : 0.0;
  double cx = ((COt[j] - 1.0 > 0.0) ? (COt[j] - 1.0) : 0.0) + XCt[j];
  double fd = pow(ec, PEC) * pow(su, PSURF)
              / (1.0 + pr / KPF) / (1.0 + cx / KCF);
  DL  += Wt[j] / ws * fd;
  FV  += Wt[j] / ws * (1.0 / (1.0 + cx / KCV) / (1.0 + pr / KPV));
  // the index is LINEAR in the states.  Using PRM/(1+PRM) would impose
  // its own saturation on top of the states' own, which is an
  // assumption about the volume effect smuggled into the read-out.
  PNI += VVt[j] / vs * (A1 * pr + A2 * cy);
  FIB += VVt[j] / vs * cx;
  if (COt[j] + XCt[j] > 1.586) VFIB += VVt[j] / vs;   // past the separatrix
}

double DLCO = DLCO0 * DL;
double FVC  = FVC0 * FV;

// low baseline reserve amplifies the SAME absolute injury into a higher
// grade; this is the only place the patient's starting lung enters
double RF   = sqrt(DLREF / ((DLCO0 > 20.0) ? DLCO0 : 20.0));
double PNIE = PNI * RF;

// NTCP for CTCAE grade >= 2 radiation pneumonitis.  Logistic in
// ln(index), so a sigmoid dose-response emerges rather than being
// assumed.  PNI50/PNISL are the model's only two fitted constants.
double NTCP = (PNIE > 1e-9)
  ? 1.0 / (1.0 + exp(-(log(PNIE) - log(PNI50)) / PNISL)) : 0.0;

// CTCAE grade as an OUTPUT of the continuous state
double DROP = 100.0 * (DLCO0 - DLCO) / DLCO0;
double GRADE = 0.0;
if      (PNIE >= 2.30 && DROP >= 45.0) GRADE = 4.0;
else if (PNIE >= 1.55)                 GRADE = 3.0;
else if (PNIE >= 0.95)                 GRADE = 2.0;
else if (PNIE >= 0.45)                 GRADE = 1.0;

double TCP = exp(-exp((TUMLN < 50.0) ? TUMLN : 50.0)) * 100.0;
double UCP = TCP / 100.0 * (1.0 - NTCP) * 100.0;

double CDEXo = DEXC / VDEX;
double CPIRo = PIRC / VPIR;
double CNINo = NINC / VNIN;
double CDURo = DURC / VDUR;
double CACEo = ACEC / VACE;

$CAPTURE @annotated
MLD   : Mean lung dose (Gy)
V5G   : Fractional lung volume >= 5 Gy
V20G  : Fractional lung volume >= 20 Gy
V40G  : Fractional lung volume >= 40 Gy
DLCO  : Diffusing capacity (% predicted)
FVC   : Forced vital capacity (% predicted)
PNI   : Pneumonitis index
PNIE  : Reserve-adjusted pneumonitis index
NTCP  : Probability of CTCAE grade >= 2 pneumonitis
GRADE : CTCAE radiation-pneumonitis grade
FIB   : Volume-weighted fibrosis score (collagen excess)
VFIB  : Fraction of lung past the fibrotic separatrix
TCP   : Tumour control probability (%)
UCP   : Uncomplicated cure probability (%)
CDEXo : Prednisolone concentration (mg/L)
CPIRo : Pirfenidone concentration (mg/L)
CNINo : Nintedanib concentration (mg/L)
CDURo : Durvalumab concentration (mg/L)
CACEo : Lisinopril concentration (mg/L)

## =====================================================================
##  SCENARIOS
##  Every number quoted below is what rili_reference_model.py produced;
##  see rili_reference_output.txt for the full run.
## =====================================================================
##
## $SETUP -------------------------------------------------------------
##   library(mrgsolve); library(dplyr); library(tidyr)
##   mod <- mread("rili_mrgsolve_model", ".")
##   run <- function(mod, ..., end = 730) {
##     mod %>% param(...) %>% mrgsim(end = end, delta = 1) %>% as_tibble()
##   }
##   # plan library: V1..V6 are the DVH, DRX/NFX the prescription
##   plans <- list(
##     imrt60   = list(DRX=60, NFX=30, V1=.40, V2=.20, V3=.12, V4=.11,
##                     V5=.09, V6=.08),
##     crt60    = list(DRX=60, NFX=30, V1=.52, V2=.13, V3=.08, V4=.08,
##                     V5=.09, V6=.10),
##     proton60 = list(DRX=60, NFX=30, V1=.58, V2=.16, V3=.09, V4=.07,
##                     V5=.055, V6=.045),
##     sbrt54   = list(DRX=54, NFX=3,  V1=.860, V2=.075, V3=.032,
##                     V4=.018, V5=.010, V6=.005),
##     hypo60   = list(DRX=60, NFX=8,  V1=.700, V2=.130, V3=.070,
##                     V4=.045, V5=.033, V6=.022),
##     imrt74   = list(DRX=74, NFX=37, V1=.34, V2=.20, V3=.13, V4=.12,
##                     V5=.11, V6=.10),
##     bath     = list(DRX=60, NFX=30, V1=.06, V2=.55, V3=.29, V4=.05,
##                     V5=.03, V6=.02),
##     hot      = list(DRX=60, NFX=30, V1=.600, V2=.110, V3=.050,
##                     V4=.050, V5=.070, V6=.120))
##
## [A] NATURAL HISTORY — the plan is the only variable -----------------
##   do.call(run, c(list(mod), plans$sbrt54))     # MLD  4.30  NTCP  2.1%
##   do.call(run, c(list(mod), plans$hypo60))     # MLD  8.64  NTCP  6.3%
##   do.call(run, c(list(mod), plans$proton60))   # MLD 12.22  NTCP  8.0%
##   do.call(run, c(list(mod), plans$imrt60))     # MLD 17.75  NTCP 16.9%
##   do.call(run, c(list(mod), plans$crt60))      # MLD 16.72  NTCP 14.4%
##   do.call(run, c(list(mod), plans$imrt74))     # MLD 24.94  NTCP 26.7%
##   # Calibration: QUANTEC MLD 13/20/24 Gy -> 10/20/30%; the model gives
##   # 8.9/21.2/29.3%.  Peak pneumonitis at 36-59 d AFTER the last
##   # fraction, i.e. inside the observed 4-12 week window.
##
## [B] SAME MEAN LUNG DOSE, DIFFERENT HISTOGRAM ------------------------
##   do.call(run, c(list(mod), plans$bath))
##   do.call(run, c(list(mod), plans$hot))
##   # MLD 15.53 vs 15.35 Gy — matched.  V40 4.9% vs 18.8% — 3.8x apart.
##   #   pneumonitis  NTCP 15.1%  vs 11.7%   (the BATH plan is worse)
##   #   fibrosis     Vfib  2.0%  vs 12.0%   (the HOT plan is 6x worse)
##   # The two endpoints rank the two plans in OPPOSITE directions at
##   # matched mean dose, because pneumonitis is a volume-weighted SUM
##   # over a parallel organ while fibrosis is a local THRESHOLD.
##
## [C] CORTICOSTEROID — fast loop only --------------------------------
##   run(mod, RDEX = 60, TDEX0 = 56)   # started before the peak
##   run(mod, RDEX = 60, TDEX0 = 120)  # started after the peak
##   # Index falls 0.6715 -> 0.6043 when started at day 56 and not at all
##   # when started at day 120.  Vfib is 8.0% in EVERY steroid arm: the
##   # drug cannot touch fibrosis, at any dose, at any time.
##   # 30 mg and 60 mg are indistinguishable (both ~125x the IC50).
##
## [D] ANTIFIBROTICS — slow loop only ---------------------------------
##   run(mod, RPIR = 2403, TPIR0 = 0,   TPIR1 = 365)   # pirfenidone
##   run(mod, RNIN = 300,  TNIN0 = 0,   TNIN1 = 365)   # nintedanib
##   run(mod, RACE = 20,   TACE0 = 0,   TACE1 = 730)   # lisinopril
##   # All three take Vfib from 8.0% to 0.0% and leave the peak
##   # pneumonitis index at 0.6715 — UNCHANGED to four decimals.  A trial
##   # with a pneumonitis endpoint cannot detect any of them.
##
## [E] RADIOPROTECTOR — temporal coincidence ---------------------------
##   run(mod, AMION = 1, AMIDEL = 0)    # NTCP  3.3%
##   run(mod, AMION = 1, AMIDEL = 15)   # NTCP 13.1%
##   run(mod, AMION = 1, AMIDEL = 30)   # NTCP 15.9%
##   run(mod, AMION = 1, AMIDEL = 60)   # NTCP 16.8%  (control 16.9%)
##   run(mod, AVAON = 1, AVADEL = 30)   # NTCP 12.5%
##   # Amifostine's entire effect is gone by 30-60 min because WR-1065
##   # has an 8-minute half-life.  The benefit is not a property of the
##   # drug, it is a property of the SCHEDULE — which is a mechanistic
##   # reading of why the amifostine trials disagree with each other.
##
## [F] CONSOLIDATION DURVALUMAB x MEAN LUNG DOSE (PACIFIC) -------------
##   for (p in c("proton60","imrt60","imrt74"))
##     do.call(run, c(list(mod), plans[[p]], list(DURDOSE = 1500)))
##   # NTCP without -> with durvalumab:
##   #   MLD 12.2   8.0% -> 9.2%    (+1.2 points)
##   #   MLD 17.8  16.9% -> 19.2%   (+2.3 points)
##   #   MLD 24.9  26.7% -> 30.2%   (+3.5 points)
##   # Because checkpoint blockade is a GAIN on the fast loop and not an
##   # added source, the excess risk GROWS with mean lung dose.  PACIFIC
##   # reports any-grade 24.8% -> 33.9% and grade 3-4 2.6% -> 3.4%; the
##   # model's grade >=2 excess sits between them.  The prediction that
##   # goes beyond the trial is that the interaction is concentrated in
##   # the high-MLD subgroup rather than spread evenly.
##   # Durvalumab dosing: ev(amt = 1500, cmt = "DURC", ii = 28, addl = 12,
##   #                       time = 56)
##
## [G] BASELINE RESERVE — identical plan, different lung ---------------
##   for (d in c(100, 85, 60, 45)) run(mod, DLCO0 = d)
##   # NTCP 14.3 / 16.9 / 23.7 / 30.6%.  Vfib is 8.0% in all four: the
##   # DAMAGE is identical and only the GRADE moves, because grade is a
##   # threshold on injury divided by reserve.
##
## [H] THERAPEUTIC RATIO ----------------------------------------------
##   # tumour BED10 and lung BED3 both rise with fraction size, but not
##   # by the same factor, and the irradiated VOLUME falls at the same
##   # time.  UCP = TCP x (1 - NTCP) is the only endpoint that sees both.
##   #   IMRT 60/30    d/fx 2.0  BEDT10  72.0  TCP  42.5  NTCP 16.9  UCP 35.4
##   #   IMRT 74/37    d/fx 2.0  BEDT10  88.8  TCP  99.5  NTCP 26.7  UCP 72.9
##   #   hypofx 60/8   d/fx 7.5  BEDT10 105.0  TCP 100.0  NTCP  6.3  UCP 93.7
##   #   SBRT 54/3     d/fx 18.0 BEDT10 151.2  TCP 100.0  NTCP  2.1  UCP 97.9
##   # SBRT wins on BOTH axes at once, which cannot be a radiobiological
##   # effect: alpha/beta 3 vs 10 moves the ratio AGAINST the lung at
##   # large fraction size.  The gain is GEOMETRIC — the mean lung dose
##   # falls from 17.8 to 4.3 Gy — and the model separates the two
##   # contributions because they enter in different places.
##
## [I] THE ANTIFIBROTIC TIME WINDOW -----------------------------------
##   for (t in c(0, 7, 14, 21, 28, 35, 42, 49, 56, 70))
##     run(mod, RPIR = 2403, TPIR0 = t, TPIR1 = t + 365, end = 900)
##   # start day  0..42 -> Vfib 0.0%, FIB900 0.007-0.013  (prevented)
##   # start day  49+   -> Vfib 8.0%, FIB900 0.235        (no effect)
##   # RT ends on day 42.  The window shuts BETWEEN day 42 and day 49 —
##   # within one week of the last fraction — and it shuts completely,
##   # not gradually: day 42 gives FIB900 0.0125 and day 49 gives 0.2347.
##   # The slow loop is bistable and cross-linked collagen is not an MMP
##   # substrate, so an antifibrotic is not a dose-response drug: it is a
##   # race against the separatrix.  Given concurrently it prevents
##   # fibrosis completely; given two months later it does nothing.  This
##   # is a falsifiable prediction, and it is the opposite of how such
##   # drugs are usually trialled (started once fibrosis is diagnosed).
##
## [J] STEROID SYMPTOM RESPONSE ---------------------------------------
##   # Index at day 56 / +14 / +28 / +56 after starting at day 56:
##   #   no steroid            0.6004  0.6586  0.6711  0.6123
##   #   prednisolone 60 mg/d  0.6004  0.5991  0.5805  0.5148
##   # Vfib 8.0% in both.  The steroid bends the acute curve and leaves
##   # the late one exactly where it was.
## =====================================================================
