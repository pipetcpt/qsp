## =====================================================================
##  ida_mrgsolve_model.R
##  Iron Deficiency Anaemia (IDA) — QSP / PK-PD model
##  37 ODE compartments · 122 parameters · 15 therapeutic scenarios
##
##  철결핍성 빈혈 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  THESIS 1 — oral iron absorption is not a fraction, it is a PRODUCT of
##  three factors, and the dose closes the third one behind itself:
##
##      absorbed flux  =  A_LUM  x  DMT1 capacity  x  FPN_ENT
##                          |            |              |
##       luminal soluble ---+            |              |
##       iron (saturable, Km 11 mg) -----+              |
##       enterocyte export capacity, set by the -------+
##       hepcidin the PREVIOUS dose generated
##
##  A 60 mg elemental dose raises serum iron 27.9 -> 224 ug/dL by 5.8 h
##  and hepcidin 0.32 -> 1.24 ng/mL (3.9x) by 10.3 h.  Hepcidin then
##  degrades enterocyte export capacity, which recovers on a ~1-2 day
##  clock (KSYN_FPE = 0.004/h), not the ~2.5 h clock of hepcidin itself.
##  So an identical probe dose absorbs only:
##
##      +4 h   75.6 %      +24 h   91.5 %      +48 h   94.2 %
##      +8 h   88.4 %      +36 h   93.1 %      +72 h   95.8 %
##
##  of what it would absorb on a rested gut.  This is the classical
##  "mucosal block" (Hahn 1943) written as a rate equation, and it makes
##  the DOSING INTERVAL, not the dose, the decisive design variable.
##
##  THESIS 2 — fractional absorption and total absorbed iron are ordered
##  in OPPOSITE directions by dosing frequency, so "best regimen" depends
##  entirely on which currency you are paying in (note C, 14-day totals,
##  60 mg elemental per dose, perfect adherence):
##
##      q72h   FIA 20.9 %   4.48 mg/day absorbed
##      q48h   FIA 20.1 %   6.03 mg/day        <- efficiency optimum
##      q24h   FIA 17.6 %  10.55 mg/day
##      q12h   FIA 13.9 %  16.72 mg/day        <- delivery optimum
##
##  Alternate-day dosing is therefore an EFFICIENCY and TOLERABILITY
##  optimum, not a delivery optimum.  Every milligram swallowed on an
##  alternate-day schedule buys 1.14x more absorbed iron than on a daily
##  schedule, but the daily schedule still delivers 1.75x more iron per
##  day.  A model that only reports fractional absorption would recommend
##  the wrong regimen for a symptomatic patient, and a model that only
##  reports total absorption would recommend the wrong regimen for a
##  patient who cannot tolerate the tablets.
##
##  THESIS 3 — haemoglobin is limited by a DIFFERENT constraint than
##  absorption, so the two come apart.  Marrow iron consumption is
##  13.3 mg/day at baseline and ~85 % of it is supplied by macrophage
##  recycling of senescent red cells, not by the gut.  Oral iron adds
##  only ~10 mg/day even when it works.  Consequently intravenous dose
##  escalation buys STORES, not SPEED (note F):
##
##      IV dose    500     1000     2000     3000  mg
##      dHb wk 1  0.35     0.37     0.38     0.39  g/dL
##      dHb wk 3  1.55     1.99     2.23     2.30  g/dL
##      max slope 0.73     0.89     0.96     0.98  g/dL/week
##      stores      72      314     1021     1720  mg
##
##  A six-fold dose range (500 -> 3000 mg) changes the maximum weekly
##  haemoglobin slope by 36 % and week-1 response by 11 %, while changing
##  iron stores 24-fold.  The marrow, not the iron supply, sets the clock.
##
##  THESIS 4 — the same 1000 mg of intravenous iron buys IDENTICAL
##  haematology and very different harm depending on the carbohydrate
##  shell (note H).  Ferric carboxymaltose and ferric derisomaltose both
##  give dHb +3.20 g/dL at 12 weeks, but carboxymaltose inhibits FGF23
##  cleavage:
##
##                        iFGF23 peak   PO4 nadir   days < 2.0 mg/dL
##      carboxymaltose      231 pg/mL    1.98        4.2
##      derisomaltose        40 pg/mL    3.02        0.0
##
##  Equal efficacy, unequal cost — a comparison only a mechanistic model
##  can make, because the efficacy endpoint cannot distinguish them.
##
##  THESIS 5 — inflammation does not merely blunt absorption, it moves the
##  iron to the wrong place.  At IL-6 40 pg/mL the gut still absorbs
##  4.5 mg/day, yet haemoglobin does not rise at all (dHb -0.07 g/dL over
##  12 weeks) while ferritin climbs 5.8 -> 55 ng/mL and TSAT stays ~12 %.
##  The absorbed iron is real; it is simply retained behind macrophage
##  ferroportin.  This is the ferritin/TSAT dissociation of functional
##  iron deficiency, produced here by mechanism rather than assumed.
## ---------------------------------------------------------------------
##  CALIBRATION SUMMARY  (60 kg woman, 3.9 L blood, 2.4 L plasma)
##
##  Replete reference state (VBLEED = 0.6 mL/day):
##      Hb 13.86 g/dL · MCH 30.3 pg · RBC 4.57 e12/L · serum iron 110 ug/dL
##      TIBC 316 · TSAT 34.8 % · hepcidin 6.43 ng/mL · ferritin 66.9 ng/mL
##      storage iron 463 mg · tissue iron 386 mg · absorption 1.24 mg/day
##      reticulocytes 38 e9/L · CHr 30.3 pg · sTfR 1.36 mg/L
##      marrow iron consumption 16.3 mg/day · total body iron 2769 mg
##
##  IDA reference state (VBLEED = 6.375 mL/day ~ 178 mL/cycle):
##      Hb 9.02 g/dL · MCH 21.2 pg · RBC 4.26 e12/L · serum iron 28 ug/dL
##      TIBC 427 · TSAT 6.5 % · hepcidin 0.32 ng/mL · ferritin 5.8 ng/mL
##      storage iron 15 mg · tissue iron 316 mg · absorption 2.53 mg/day
##      reticulocytes 42 e9/L · CHr 20.5 pg · sTfR 3.63 mg/L · EPO 94 mIU/mL
##      marrow iron consumption 13.3 mg/day · total body iron 1585 mg
##
##  Both states are the model's own equilibria, not imposed values: the
##  IDA state is what VBLEED = 6.375 mL/day produces when the model is
##  integrated to steady state from the replete state.
##
##  Literature anchors reproduced (see ida_references.md):
##    - fractional absorption of 60 mg elemental in iron-depleted women
##      ~20-22 % (model 22.0 %)                     [Stoffel 2017/2020]
##    - hepcidin still elevated 24 h after a 60 mg dose (model 1.36x)
##                                                   [Moretti 2015]
##    - a second dose given 4 h later is poorly absorbed (model 75.6 %)
##                                                   [Moretti 2015]
##    - daily : alternate-day fractional absorption ratio 0.75 observed,
##      0.87 in this model — the model is CONSERVATIVE about the size of
##      the refractory penalty (documented limitation)
##                                                   [Stoffel 2017]
##    - serum iron peak 2-5 h after an oral dose (model 5.8 h, late)
##    - reticulocyte peak day 7-14 (model day 12-14)
##    - FCM hypophosphataemia: nadir ~1.9-2.0 mg/dL near day 7-14,
##      iFGF23 peak day 1-2, absent with derisomaltose
##                                                   [Wolf 2020 PHOSPHARE]
##    - mass balance verified exact (0.0000 % error over 12 weeks)
##
##  IMPORTANT — this is an educational/research model. Parameters are
##  literature-informed approximations, not a validated patient model.
## =====================================================================

library(mrgsolve)
library(dplyr)

ida_code <- '
$PROB
# Iron Deficiency Anaemia QSP model
- 37 compartments: gut lumen, enterocyte labile + ferritin pools, plasma
  transferrin iron, NTBI, IV colloid, macrophage and hepatocyte stores,
  non-erythroid tissue iron, transferrin, hepcidin, two ferroportin
  pools, erythroferrone, IL-6, a 6-stage erythron carrying cells AND
  their haemoglobin mass separately, biomarkers, and the FGF23-phosphate
  axis.
- Time unit: HOURS. Iron amounts: mg. Cells: 1e12 cells. Hb mass: g.

$GLOBAL
#define PV_dL (PV_L * 10.0)
#define BV_dL (BV_L * 10.0)
// serum iron (ug/dL) from the plasma transferrin-iron amount (mg)
#define SI_    (A_TF * 1000.0 / PV_dL)
#define TSAT_  (100.0 * SI_ / (TIBC > 1.0 ? TIBC : 1.0))
#define HBMASS (HRBC + HRET)
#define HB_    (HBMASS / BV_dL)
#define NCELL  (NRBC + NRET)

$PARAM @annotated
// ---- body size -----------------------------------------------------
BV_L      : 3.9      : Blood volume (L)
PV_L      : 2.4      : Plasma volume (L)

// ---- diet and oral absorption --------------------------------------
DIET      : 16.0     : Dietary iron intake (mg/day)
F_DIET    : 0.30     : Bioavailable fraction of dietary iron
F_ORAL    : 0.90     : Fraction of an oral salt dose reaching the window
KTR       : 0.35     : Luminal iron swept past the duodenum (1/h)
VMAX_D    : 5.57031  : Apical DMT1 maximum uptake velocity (mg/h)
KM_D      : 11.0     : DMT1 Michaelis constant (mg in window)
K_APIC    : 45.0     : Hepcidin IC50 on apical uptake (ng/mL)
K_IRP     : 30.25    : Enterocyte iron for IRP down-regulation (mg)
KEXP      : 2.035    : Basolateral export rate constant at FPN_ENT=1 (1/h)
KFT       : 0.120    : Labile iron into enterocyte ferritin (1/h)
KMOB      : 0.136    : Enterocyte ferritin back to labile at FPN_ENT=1 (1/h)
KSHED     : 0.0140   : Enterocyte ferritin iron lost by shedding (1/h)

// ---- plasma iron and storage kinetics ------------------------------
K_LIV0    : 0.020    : Basal plasma to hepatocyte uptake (1/h)
K_LIV1    : 0.289258 : Holo-Tf driven hepatocyte uptake (1/h)
KM_LIVT   : 45.0     : TSAT at half-maximal hepatocyte uptake (percent)
K_LIV_EXP : 0.01111  : Hepatocyte iron export at FPN_RES=1 (1/h)
K_RES_EXP : 0.0724   : Macrophage storage mobilisation at FPN_RES=1 (1/h)
KM_FPN_REL: 0.0969   : FPN_RES at which half of recycled iron goes direct
K_RES_LIV : 0.00013  : Macrophage to hepatocyte transfer (1/h)
F_COL_RES : 0.98     : Fraction of IV colloid handled by the RES
K_COL     : 0.0693   : IV colloid clearance from plasma (1/h)
K_TISS    : 0.030    : Plasma to non-erythroid tissue (1/h)
TISS_CAP  : 450.0    : Non-erythroid tissue iron capacity (mg)
K_TISS_OUT: 2.95e-05 : Tissue iron back to plasma (1/h)
TSAT_NTBI : 75.0     : TSAT above which NTBI appears (percent)
K_NTBI_ON : 0.60     : NTBI formation above threshold (1/h)
K_NTBI_CL : 1.20     : NTBI clearance to liver and tissue (1/h)
LOSS_BASAL: 0.0417   : Obligatory epithelial iron loss (mg/h)

// ---- transferrin ---------------------------------------------------
TIBC_MIN  : 265.0    : TIBC when replete (ug/dL)
TIBC_MAX  : 470.0    : TIBC when fully depleted (ug/dL)
KI_TIBC   : 22.0     : Ferritin for half-maximal TIBC up-regulation
IL6_TIBC  : 0.30     : Maximal fractional TIBC suppression by IL-6
KI6_TIBC  : 25.0     : IL-6 for half-maximal TIBC suppression
KT_TIBC   : 0.0060   : Transferrin turnover (1/h)

// ---- hepcidin ------------------------------------------------------
KOUT_HEP  : 0.277    : Hepcidin elimination (1/h), t1/2 2.5 h
KSYN_HEP  : 0.725    : Hepcidin synthesis scale (ng/mL/h)
HEP_STORE0: 0.30     : Store-independent floor of the BMP6 signal
HEP_STORE_E:2.50     : Maximal store (BMP6/SMAD) drive
KM_HEP_LIV: 260.0    : Hepatocyte iron at half-maximal BMP6 signal (mg)
HEP_TSAT0 : 0.30     : Floor of the holo-Tf (TfR2/HFE) signal
HEP_TSAT_E: 1.70     : Maximal holo-Tf drive
KM_HEP_TSAT:25.0     : TSAT at half-maximal holo-Tf drive (percent)
HILL_TSAT : 2.0      : Hill coefficient of holo-Tf sensing
EMAX_IL6  : 8.0      : Maximal IL-6/STAT3 induction of hepcidin
KM_IL6    : 18.0     : IL-6 at half-maximal STAT3 induction (pg/mL)
KI_ERFE   : 9.0      : Erythroferrone for half-maximal suppression
TMPRSS6   : 1.0      : Matriptase-2 activity (1 normal, <1 IRIDA)

// ---- ferroportin: two clocks ---------------------------------------
KSYN_FPE  : 0.0040   : Recovery of enterocyte export capacity (1/h)
KDEG_FPE  : 0.0300   : Hepcidin-driven loss of enterocyte capacity
KSYN_FPR  : 0.0350   : Macrophage/hepatocyte ferroportin resynthesis (1/h)
KDEG_FPR  : 0.1200   : Hepcidin-driven ferroportin degradation
HILL_FPN  : 0.60     : Power-law exponent of the hepcidin-FPN relationship

// ---- erythroferrone and inflammation -------------------------------
KSYN_ERFE : 0.30     : ERFE synthesis scale (ng/mL/h)
KOUT_ERFE : 0.150    : ERFE elimination (1/h)
KOUT_IL6  : 0.050    : IL-6 elimination (1/h)
IL6_IN    : 0.0      : IL-6 zero-order input (pg/mL/h)
IMAX_IL6_EPO:0.45    : Maximal IL-6 blunting of EPO production
KM_IL6_EPO: 25.0     : IL-6 for half-maximal EPO blunting
IMAX_IL6_ERY:0.35    : Maximal cytokine suppression of progenitors
KM_IL6_ERY: 25.0     : IL-6 for half-maximal progenitor suppression
E_IL6_RBC : 0.50     : Maximal inflammatory shortening of RBC survival
KM_IL6_RBC: 25.0     : IL-6 for half-maximal RBC survival shortening

// ---- erythropoiesis ------------------------------------------------
KIN_PROG  : 0.00684375: Baseline progenitor influx (1e12 cells/h)
KPROG     : 0.02083  : Progenitor to erythroblast transit (1/h)
KEB       : 0.025    : Erythroblast stage transit (1/h), 3 x 40 h
KRET      : 0.04167  : Blood reticulocyte maturation (1/h), 24 h
KRBC      : 0.000347 : RBC senescence (1/h), 120-day lifespan
EPO0      : 10.0     : Baseline EPO (mIU/mL)
HB_REF    : 13.5     : Hb at which EPO equals EPO0 (g/dL)
KEPO      : 0.50     : Log-linear EPO sensitivity to Hb (dL/g)
EPO_MAX   : 1200.0   : EPO cap (mIU/mL)
EMAX_EPO  : 2.50     : Maximal fold expansion of progenitor influx
EC50_EPO  : 100.0    : EPO above baseline for half-maximal expansion
K_EXP_FE  : 12.0     : TSAT at which half the EPO expansion is realised
VHB       : 0.267539 : Maximal Hb synthesis per cell (pg/cell/h)
K_TSAT_FE : 3.50     : TSAT at half-maximal haemoglobinisation (percent)
W1        : 0.50     : Hb synthesis weight, proerythroblast
W2        : 1.00     : Hb synthesis weight, basophilic-polychromatic
W3        : 1.50     : Hb synthesis weight, orthochromatic
WRET      : 0.30     : Hb synthesis weight, reticulocyte
KAPO_MAX  : 0.0060   : Maximal iron-restriction erythroblast apoptosis (1/h)
FE_PER_HB : 3.47     : mg iron per g haemoglobin

// ---- blood loss ----------------------------------------------------
VBLEED    : 6.375    : Ongoing blood loss (mL/day)

// ---- biomarkers ----------------------------------------------------
KT_FERR   : 0.030    : Serum ferritin turnover (1/h)
FERR_FLOOR: 3.0      : Ferritin floor (ng/mL)
K_FERR_RES: 3.50     : mg macrophage iron per ng/mL ferritin
K_FERR_LIV: 8.50     : mg hepatocyte iron per ng/mL ferritin
FERR_IL6  : 32.0     : Maximal acute-phase ferritin rise (ng/mL)
KM_FERR_IL6:20.0     : IL-6 for half-maximal acute-phase rise
KT_STFR   : 0.010    : Soluble transferrin receptor turnover (1/h)
STFR0     : 0.90     : sTfR with no iron restriction (mg/L)
STFR_E    : 5.00     : Maximal iron-restriction driven sTfR rise

// ---- FGF23 / phosphate axis ----------------------------------------
KSYN_FGF  : 8.0      : iFGF23 synthesis (pg/mL/h)
KOUT_FGF  : 0.20     : iFGF23 elimination (1/h)
E_CLV     : 6.00     : Maximal FGF23 cleavage-inhibition effect
KM_CLV    : 210.0    : Cleavage-inhibition signal at half-maximal effect
K_CLV_IN  : 1.0      : Colloid to cleavage-inhibition signal (1/h)
KOUT_CLV  : 0.0050   : Decay of the cleavage-inhibition signal (1/h)
FCM_FGF   : 1.0      : 1 = carboxymaltose, 0 = derisomaltose/sucrose
PHOS0     : 3.60     : Baseline serum phosphate (mg/dL)
KOUT_PHOS : 0.025    : Phosphate turnover (1/h)
IMAX_FGF_P: 0.72     : Maximal FGF23-driven fall in phosphate input
IC50_FGF_P: 190.0    : iFGF23 for half-maximal phosphaturia (pg/mL)
CTRIOL0   : 42.0     : Baseline 1,25(OH)2D (pg/mL)
KOUT_CTRIOL:0.060    : Calcitriol turnover (1/h)
IMAX_FGF_D: 0.60     : Maximal FGF23 suppression of 1-alpha-hydroxylase
IC50_FGF_D: 160.0    : iFGF23 for half-maximal calcitriol suppression
F_CTRIOL_P: 0.35     : Calcitriol-dependent share of phosphate input
PTH0      : 42.0     : Baseline PTH (pg/mL)
KOUT_PTH  : 0.120    : PTH turnover (1/h)
E_PTH_D   : 1.50     : PTH rise as calcitriol falls

// ---- GI tolerability -----------------------------------------------
KT_GI     : 0.080    : GI symptom score turnover (1/h)
K_GI      : 0.55     : GI score per mg/h of unabsorbed luminal iron
EMAX_ADH  : 0.45     : Maximal loss of adherence
K50_ADH   : 3.0      : GI score at half-maximal loss of adherence

// ---- symptom endpoints ---------------------------------------------
FACIT_MAX : 52.0     : FACIT-Fatigue score ceiling
KM_FACIT_HB:2.0      : Hb for half-maximal fatigue benefit (g/dL)
W_FACIT_TISS:0.35    : Weight of tissue iron on the fatigue score
TISS_REF  : 380.0    : Reference non-erythroid tissue iron (mg)
IRLS_MAX  : 32.0     : IRLS restless-legs severity ceiling

$CMT @annotated
A_LUM   : Luminal soluble iron in the absorptive window (mg)
A_ENT   : Enterocyte labile iron (mg)
A_EFT   : Enterocyte ferritin iron, the mucosal-block memory (mg)
A_TF    : Plasma transferrin-bound iron (mg)
A_NTBI  : Non-transferrin-bound iron (mg)
A_COL   : Plasma iron-carbohydrate colloid (mg)
A_RES   : Macrophage storage iron (mg)
A_LIV   : Hepatocyte storage iron (mg)
A_TISS  : Non-erythroid tissue iron (mg)
TIBC    : Total iron-binding capacity (ug/dL)
HEP     : Serum hepcidin-25 (ng/mL)
FPN_ENT : Enterocyte iron-export capacity (fraction of maximum)
FPN_RES : Macrophage and hepatocyte ferroportin (fraction of maximum)
ERFE    : Erythroferrone (ng/mL)
IL6     : Interleukin-6 (pg/mL)
PROG    : BFU-E and CFU-E progenitors (1e12 cells)
N1      : Proerythroblasts (1e12 cells)
N2      : Basophilic-polychromatic erythroblasts (1e12 cells)
N3      : Orthochromatic erythroblasts (1e12 cells)
H1      : Haemoglobin mass in stage-1 erythroblasts (g)
H2      : Haemoglobin mass in stage-2 erythroblasts (g)
H3      : Haemoglobin mass in stage-3 erythroblasts (g)
NRET    : Blood reticulocytes (1e12 cells)
HRET    : Haemoglobin mass in reticulocytes (g)
NRBC    : Circulating mature red cells (1e12 cells)
HRBC    : Haemoglobin mass in mature red cells (g)
FERR    : Serum ferritin (ng/mL)
STFR    : Soluble transferrin receptor (mg/L)
CLV     : FGF23 cleavage-inhibition signal (mg-equivalents)
FGF23   : Intact FGF23 (pg/mL)
PHOS    : Serum phosphate (mg/dL)
CTRIOL  : 1,25(OH)2 vitamin D (pg/mL)
PTH     : Parathyroid hormone (pg/mL)
GI      : GI symptom score (arbitrary units)
CUM_ABS : Cumulative iron absorbed across the basolateral membrane (mg)
CUM_IV  : Cumulative intravenous iron delivered (mg)
CUM_LOSS: Cumulative systemic iron loss (mg)

$MAIN
// ---- default initial condition: the calibrated IDA reference state ----
// (re-equilibrate by simulating with no dose for ~2-3 years to obtain the
//  steady state belonging to any other parameter set)
A_LUM_0   = 0.241813;  A_ENT_0  = 0.247647;  A_EFT_0 = 0.699696;
A_TF_0    = 0.668953;  A_NTBI_0 = 0.0;       A_COL_0 = 0.0;
A_RES_0   = 5.78149;   A_LIV_0  = 9.59032;   A_TISS_0 = 315.948;
TIBC_0    = 427.236;   HEP_0    = 0.317682;
FPN_ENT_0 = 0.209516;  FPN_RES_0 = 0.367185;
ERFE_0    = 2.77898;   IL6_0    = 0.0;
PROG_0    = 0.418298;  N1_0 = 0.321552;  N2_0 = 0.29668;  N3_0 = 0.273745;
H1_0      = 1.03353;   H2_0 = 2.8612;    H3_0 = 5.28086;
NRET_0    = 0.163967;  HRET_0 = 3.36915;
NRBC_0    = 16.4533;   HRBC_0 = 348.531;
FERR_0    = 5.78322;   STFR_0 = 3.63274;
CLV_0     = 0.0;       FGF23_0 = 40.0;   PHOS_0 = 3.01695;
CTRIOL_0  = 36.96;     PTH_0  = 49.56;   GI_0   = 0.649231;

// ---- oral bioavailability, reduced by GI-driven non-adherence -------
F_A_LUM = F_ORAL * (1.0 - EMAX_ADH * GI / (K50_ADH + GI));

$ODE
// =====================================================================
//  derived quantities
// =====================================================================
double si   = SI_;
double tsat = TSAT_;
double hb   = HB_;
double f_fe = tsat / (tsat + K_TSAT_FE);          // haemoglobinisation gate
double vhb  = VHB * f_fe;                         // pg/cell/h

// EPO: log-linear in the haemoglobin deficit, blunted by inflammation
double epo = EPO0 * exp(KEPO * (HB_REF - hb))
             * (1.0 - IMAX_IL6_EPO * IL6 / (KM_IL6_EPO + IL6));
if (epo > EPO_MAX) epo = EPO_MAX;
double dep = epo - EPO0; if (dep < 0.0) dep = 0.0;
double epo_eff = 1.0 + EMAX_EPO * dep / (EC50_EPO + dep);

// iron-restricted erythropoiesis: EPO cannot expand a marrow it cannot
// supply, which is why the reticulocyte count in IDA is inappropriately
// low for the degree of anaemia
double f_exp = tsat / (tsat + K_EXP_FE);
double epo_real = 1.0 + (epo_eff - 1.0) * f_exp;

// marrow haemoglobin synthesis (g/h) and its iron cost (mg/h)
double wsum   = W1 * N1 + W2 * N2 + W3 * N3 + WRET * NRET;
double hbsyn  = wsum * vhb;
double fe_marrow = FE_PER_HB * hbsyn;

// =====================================================================
//  gut: lumen, enterocyte labile pool, enterocyte ferritin pool
// =====================================================================
double f_apic = (1.0 / (1.0 + HEP / K_APIC))
                * (1.0 / (1.0 + (A_ENT + A_EFT) / K_IRP));
double v_dmt1 = VMAX_D * A_LUM / (KM_D + A_LUM) * f_apic;
double v_abs  = KEXP * FPN_ENT * A_ENT;      // THE absorbed flux
double v_ft   = KFT * A_ENT;                 // labile -> ferritin
double v_mob  = KMOB * FPN_ENT * A_EFT;      // ferritin -> labile
double v_shed = KSHED * A_EFT;               // lost with the shed cell

dxdt_A_LUM = DIET * F_DIET / 24.0 - KTR * A_LUM - v_dmt1;
dxdt_A_ENT = v_dmt1 - v_abs - v_ft + v_mob;
dxdt_A_EFT = v_ft - v_mob - v_shed;

// =====================================================================
//  systemic iron transport
// =====================================================================
double liv_up  = (K_LIV0 + K_LIV1 * tsat / (KM_LIVT + tsat)) * A_TF;
double liv_exp = K_LIV_EXP * FPN_RES * A_LIV;
double res_exp = K_RES_EXP * FPN_RES * A_RES;
double res_liv = K_RES_LIV * A_RES;
double col_out = K_COL * A_COL;
double exc     = tsat - TSAT_NTBI;
double ntbi_on = (exc > 0.0) ? K_NTBI_ON * (exc / 100.0) * A_TF : 0.0;
double ntbi_cl = K_NTBI_CL * A_NTBI;
double tiss_up = K_TISS * A_TF * (1.0 - A_TISS / TISS_CAP);
if (tiss_up < 0.0) tiss_up = 0.0;
double tiss_out = K_TISS_OUT * A_TISS;

// erythrophagocytosis of senescent red cells: the dominant iron flux
double krbc_eff = KRBC * (1.0 + E_IL6_RBC * IL6 / (KM_IL6_RBC + IL6));
double ephago   = FE_PER_HB * krbc_eff * HRBC;
double f_rel    = FPN_RES / (FPN_RES + KM_FPN_REL);
double rec_rel  = ephago * f_rel;      // released straight back to plasma
double rec_sto  = ephago - rec_rel;    // stored as macrophage ferritin

// iron recovered from iron-restricted erythroblast apoptosis
double kapo    = KAPO_MAX * (1.0 - f_fe);
double apo_fe  = FE_PER_HB * kapo * (H1 + H2 + H3);

double fbleed = VBLEED / (1000.0 * 24.0 * BV_L);   // fraction of blood/h

dxdt_A_TF = v_abs + rec_rel + res_exp + liv_exp + tiss_out
            + (1.0 - F_COL_RES) * col_out
            - fe_marrow - liv_up - tiss_up - ntbi_on - LOSS_BASAL;
dxdt_A_NTBI = ntbi_on - ntbi_cl;
dxdt_A_COL  = -col_out;
dxdt_A_RES  = F_COL_RES * col_out + rec_sto + apo_fe - res_exp - res_liv;
dxdt_A_LIV  = liv_up + res_liv + 0.7 * ntbi_cl - liv_exp;
dxdt_A_TISS = tiss_up + 0.3 * ntbi_cl - tiss_out;

// =====================================================================
//  transferrin
// =====================================================================
double tibc_t = (TIBC_MIN + (TIBC_MAX - TIBC_MIN) * KI_TIBC / (KI_TIBC + FERR))
                * (1.0 - IL6_TIBC * IL6 / (KI6_TIBC + IL6));
dxdt_TIBC = KT_TIBC * (tibc_t - TIBC);

// =====================================================================
//  hepcidin: four inputs on one promoter
// =====================================================================
double f_store = HEP_STORE0 + HEP_STORE_E * A_LIV / (KM_HEP_LIV + A_LIV);
double tsn     = pow(tsat, HILL_TSAT);
double f_tsat  = HEP_TSAT0 + HEP_TSAT_E * tsn / (pow(KM_HEP_TSAT, HILL_TSAT) + tsn);
double f_il6   = 1.0 + EMAX_IL6 * IL6 / (KM_IL6 + IL6);
double f_erfe  = 1.0 / (1.0 + ERFE / KI_ERFE);
dxdt_HEP = KSYN_HEP * f_store * f_tsat * f_il6 * f_erfe / TMPRSS6
           - KOUT_HEP * HEP;

// ferroportin, two clocks: the enterocyte recovers on a 1-2 day clock
// (mucosal block), macrophages within hours (acute hypoferraemia)
double hpow = (HEP > 0.0) ? pow(HEP, HILL_FPN) : 0.0;
dxdt_FPN_ENT = KSYN_FPE * (1.0 - FPN_ENT) - KDEG_FPE * hpow * FPN_ENT;
dxdt_FPN_RES = KSYN_FPR * (1.0 - FPN_RES) - KDEG_FPR * hpow * FPN_RES;

dxdt_ERFE = KSYN_ERFE * (epo_eff - 0.75) - KOUT_ERFE * ERFE;
dxdt_IL6  = IL6_IN - KOUT_IL6 * IL6;

// =====================================================================
//  erythron: cells and their haemoglobin tracked in parallel chains
// =====================================================================
double f_il6_ery = 1.0 - IMAX_IL6_ERY * IL6 / (KM_IL6_ERY + IL6);

dxdt_PROG = KIN_PROG * f_il6_ery * epo_real - KPROG * PROG - kapo * PROG;
dxdt_N1   = KPROG * PROG - KEB * N1 - kapo * N1;
dxdt_N2   = KEB * N1 - KEB * N2 - kapo * N2;
dxdt_N3   = KEB * N2 - KEB * N3 - kapo * N3;
dxdt_H1   = W1 * N1 * vhb - KEB * H1 - kapo * H1;
dxdt_H2   = KEB * H1 + W2 * N2 * vhb - KEB * H2 - kapo * H2;
dxdt_H3   = KEB * H2 + W3 * N3 * vhb - KEB * H3 - kapo * H3;
dxdt_NRET = KEB * N3 - KRET * NRET - fbleed * NRET;
dxdt_HRET = KEB * H3 + WRET * NRET * vhb - KRET * HRET - fbleed * HRET;
dxdt_NRBC = KRET * NRET - krbc_eff * NRBC - fbleed * NRBC;
dxdt_HRBC = KRET * HRET - krbc_eff * HRBC - fbleed * HRBC;

// =====================================================================
//  biomarkers
// =====================================================================
double ferr_t = FERR_FLOOR + A_RES / K_FERR_RES + A_LIV / K_FERR_LIV
                + FERR_IL6 * IL6 / (KM_FERR_IL6 + IL6);
dxdt_FERR = KT_FERR * (ferr_t - FERR);
double stfr_t = STFR0 + STFR_E * (1.0 - f_fe) * (0.5 + 0.5 * epo_eff);
dxdt_STFR = KT_STFR * (stfr_t - STFR);

// =====================================================================
//  FGF23 - phosphate - vitamin D axis (carboxymaltose-specific)
// =====================================================================
dxdt_CLV = K_CLV_IN * FCM_FGF * col_out - KOUT_CLV * CLV;
double e_clv = 1.0 + E_CLV * CLV / (KM_CLV + CLV);
dxdt_FGF23 = KSYN_FGF * e_clv - KOUT_FGF * FGF23;
double f_phos_d = (1.0 - F_CTRIOL_P) + F_CTRIOL_P * CTRIOL / CTRIOL0;
double inh_p = 1.0 - IMAX_FGF_P * FGF23 / (IC50_FGF_P + FGF23);
dxdt_PHOS = KOUT_PHOS * (PHOS0 * f_phos_d * inh_p - PHOS);
double inh_d = 1.0 - IMAX_FGF_D * FGF23 / (IC50_FGF_D + FGF23);
dxdt_CTRIOL = KOUT_CTRIOL * (CTRIOL0 * inh_d - CTRIOL);
double pth_t = PTH0 * (1.0 + E_PTH_D * (1.0 - CTRIOL / CTRIOL0));
dxdt_PTH = KOUT_PTH * (pth_t - PTH);

// =====================================================================
//  tolerability and bookkeeping
// =====================================================================
double unabs = KTR * A_LUM + v_shed;
dxdt_GI = K_GI * unabs - KT_GI * GI;

dxdt_CUM_ABS  = v_abs;
dxdt_CUM_IV   = 0.0;
dxdt_CUM_LOSS = LOSS_BASAL + FE_PER_HB * fbleed * HBMASS;

$TABLE
si   = SI_;
tsat = TSAT_;
hb   = HB_;
double ncell = (NCELL > 1e-9) ? NCELL : 1e-9;
f_fe = tsat / (tsat + K_TSAT_FE);
epo = EPO0 * exp(KEPO * (HB_REF - hb))
             * (1.0 - IMAX_IL6_EPO * IL6 / (KM_IL6_EPO + IL6));
if (epo > EPO_MAX) epo = EPO_MAX;
dep = epo - EPO0; if (dep < 0.0) dep = 0.0;
epo_eff  = 1.0 + EMAX_EPO * dep / (EC50_EPO + dep);
f_exp = tsat / (tsat + K_EXP_FE);
epo_real = 1.0 + (epo_eff - 1.0) * f_exp;
wsum  = W1 * N1 + W2 * N2 + W3 * N3 + WRET * NRET;
hbsyn = wsum * VHB * f_fe;

capture SI      = si;
capture TSAT    = tsat;
capture HB      = hb;
capture RBCC    = ncell / BV_L;                       // 1e12/L
capture MCH     = HBMASS / ncell;                     // pg
capture RETPC   = 100.0 * NRET / ncell;               // percent
capture RET_ABS = NRET / BV_L * 1000.0;               // 1e9/L
capture CHR     = HRET / ((NRET > 1e-9) ? NRET : 1e-9);
capture F_FE    = f_fe;
capture F_EXP   = f_exp;
capture EPO     = epo;
capture EPO_EFF = epo_eff;
capture EPO_REAL= epo_real;
capture HBSYN   = hbsyn;
capture FE_MARROW = FE_PER_HB * hbsyn;
capture V_ABS   = KEXP * FPN_ENT * A_ENT;
capture ABS_DAY = KEXP * FPN_ENT * A_ENT * 24.0;
capture ADH     = 1.0 - EMAX_ADH * GI / (K50_ADH + GI);
capture STORES  = A_RES + A_LIV;
capture BODY_FE = A_TF + A_NTBI + A_COL + A_RES + A_LIV + A_TISS
                  + FE_PER_HB * (HBMASS + H1 + H2 + H3);
capture FACIT   = FACIT_MAX * ((1.0 - W_FACIT_TISS) * (hb / (hb + KM_FACIT_HB))
                  + W_FACIT_TISS * ((A_TISS / TISS_REF < 1.0) ? A_TISS / TISS_REF : 1.0));
capture IRLS    = IRLS_MAX * (1.0 - ((A_TISS / TISS_REF < 1.0) ? A_TISS / TISS_REF : 1.0));
'

ida <- mcode("ida", ida_code)

## =====================================================================
##  dosing helpers
## =====================================================================
##  Oral iron enters A_LUM (compartment 1) with bioavailability F_ORAL
##  further reduced by GI-driven non-adherence.  IV iron enters A_COL
##  (compartment 6) and is additionally recorded in CUM_IV.
## ---------------------------------------------------------------------
oral_ev <- function(dose_mg, interval_h, n_dose, start_h = 0) {
  ev(time = start_h, amt = dose_mg, cmt = "A_LUM",
     ii = interval_h, addl = n_dose - 1)
}
iv_ev <- function(dose_mg, times_h = 0) {
  Reduce(c, lapply(times_h, function(t)
    ev(time = t, amt = dose_mg, cmt = "A_COL") +
    ev(time = t, amt = dose_mg, cmt = "CUM_IV")))
}

run_sim <- function(mod = ida, events = NULL, end_h = 12 * 7 * 24,
                    delta = 6, ...) {
  pars <- list(...)
  m <- if (length(pars)) do.call(param, c(list(mod), pars)) else mod
  if (is.null(events)) mrgsim(m, end = end_h, delta = delta, hmax = 0.5)
  else                 mrgsim(m, events = events, end = end_h,
                              delta = delta, hmax = 0.5)
}

## =====================================================================
##  SCENARIO 0 — reference states
## =====================================================================
##  The shipped initial condition IS the calibrated IDA steady state, so
##  an untreated run should be nearly flat apart from the slow decline
##  driven by ongoing blood loss.
## ---------------------------------------------------------------------
sc0 <- run_sim(end_h = 12 * 7 * 24)

## Replete reference state: switch off the bleeding and equilibrate.
## (2 years is enough for stores, tissue iron and transferrin to settle.)
replete <- ida %>% param(VBLEED = 0.6) %>%
  mrgsim(end = 730 * 24, delta = 24, hmax = 1) %>% as_tibble() %>% tail(1)

## =====================================================================
##  SCENARIOS 1-15 — twelve weeks of therapy
## =====================================================================
WK <- 7 * 24
scenarios <- list(
  "1  no treatment"             = list(ev = NULL,                          p = list()),
  "2  oral 65 mg daily"         = list(ev = oral_ev(65,  24, 84),          p = list()),
  "3  oral 65 mg BID"           = list(ev = oral_ev(65,  12, 168),         p = list()),
  "4  oral 130 mg daily"        = list(ev = oral_ev(130, 24, 84),          p = list()),
  "5  oral 130 mg alternate"    = list(ev = oral_ev(130, 48, 42),          p = list()),
  "6  oral 65 mg alternate"     = list(ev = oral_ev(65,  48, 42),          p = list()),
  "7  oral 195 mg alternate"    = list(ev = oral_ev(195, 48, 42),          p = list()),
  "8  IV FCM 1000 mg"           = list(ev = iv_ev(1000, 0),                p = list()),
  "9  IV FCM 750 mg x2"         = list(ev = iv_ev(750, c(0, WK)),          p = list()),
  "10 IV derisomaltose 1000 mg" = list(ev = iv_ev(1000, 0),                p = list(FCM_FGF = 0)),
  "11 IV sucrose 200 mg x5"     = list(ev = iv_ev(200, WK * (0:4)),        p = list(FCM_FGF = 0.10)),
  "12 oral daily + IL-6 20"     = list(ev = oral_ev(65, 24, 84),           p = list(IL6_IN = 1.0)),
  "13 IV FCM 1000 + IL-6 20"    = list(ev = iv_ev(1000, 0),               p = list(IL6_IN = 1.0)),
  "14 oral daily + haemostasis" = list(ev = oral_ev(65, 24, 84),           p = list(VBLEED = 0.6)),
  "15 haemostasis alone"        = list(ev = NULL,                          p = list(VBLEED = 0.6))
)

sim_scenario <- function(s) {
  m <- ida
  if (length(s$p)) m <- do.call(param, c(list(m), s$p))
  out <- if (is.null(s$ev)) mrgsim(m, end = 12 * WK, delta = 12, hmax = 0.5)
         else               mrgsim(m, events = s$ev, end = 12 * WK,
                                   delta = 12, hmax = 0.5)
  as_tibble(out)
}

at_h <- function(d, h) d[which.min(abs(d$time - h)), ]

scenario_table <- lapply(names(scenarios), function(nm) {
  d <- sim_scenario(scenarios[[nm]])
  hb0 <- d$HB[1]
  tibble(
    scenario = nm,
    Hb_w2  = at_h(d, 2  * WK)$HB,
    Hb_w4  = at_h(d, 4  * WK)$HB,
    Hb_w8  = at_h(d, 8  * WK)$HB,
    Hb_w12 = at_h(d, 12 * WK)$HB,
    dHb    = at_h(d, 12 * WK)$HB - hb0,
    TSAT   = at_h(d, 12 * WK)$TSAT,
    ferritin = at_h(d, 12 * WK)$FERR,
    MCH    = at_h(d, 12 * WK)$MCH,
    stores = at_h(d, 12 * WK)$STORES,
    ret_peak = max(d$RET_ABS),
    PO4_nadir = min(d$PHOS),
    absorbed  = at_h(d, 12 * WK)$CUM_ABS - d$CUM_ABS[1]
  )
}) %>% bind_rows()

## =====================================================================
##  NOTE A — a single 60 mg elemental dose: the gate closing
## =====================================================================
noteA <- ida %>% mrgsim(events = oral_ev(60, 24, 1), end = 72,
                        delta = 0.25, hmax = 0.05) %>% as_tibble()
##  expected: serum iron 27.9 -> 224 ug/dL at ~5.8 h; hepcidin
##  0.32 -> 1.24 ng/mL (3.9x) at ~10.3 h and still 1.36x at 24 h;
##  FPN_ENT 0.210 -> 0.170 (-19 %), only 81 % recovered at 24 h and
##  85 % at 48 h; fractional absorption 22.0 %.

## =====================================================================
##  NOTE B — the refractory window, measured with a probe dose
## =====================================================================
##  Each probe is scored against a paired run WITHOUT the probe, so the
##  number is the iron the probe itself delivered.
##  EMAX_ADH is switched off here so that the number reflects absorptive
##  PHYSIOLOGY alone; the tolerability penalty is quantified separately.
probe_absorbed <- function(gap_h) {
  m <- ida %>% param(EMAX_ADH = 0)
  both <- m %>% mrgsim(events = ev(time = 0,     amt = 60, cmt = "A_LUM") +
                                ev(time = gap_h, amt = 60, cmt = "A_LUM"),
                       end = gap_h + 72, delta = 24, hmax = 0.05) %>% as_tibble()
  cond <- m %>% mrgsim(events = ev(time = 0, amt = 60, cmt = "A_LUM"),
                       end = gap_h + 72, delta = 24, hmax = 0.05) %>% as_tibble()
  tail(both$CUM_ABS, 1) - tail(cond$CUM_ABS, 1)
}
rested_probe <- probe_absorbed(1e5)   # a gap long enough to be a rested gut
noteB <- tibble(gap_h = c(4, 8, 12, 24, 36, 48, 72)) %>%
  mutate(absorbed_mg = vapply(gap_h, probe_absorbed, numeric(1)),
         pct_of_rested = 100 * absorbed_mg / 13.181)
##  expected, as a percentage of a rested-gut probe (13.18 mg):
##  +4 h 75.6 · +8 h 88.4 · +12 h 89.5 · +24 h 91.5 · +36 h 93.1 ·
##  +48 h 94.2 · +72 h 95.8

## =====================================================================
##  NOTE C — 14-day regimens: fractional and total absorption diverge
## =====================================================================
regimens <- tribble(
  ~label,          ~dose, ~interval, ~n,
  "60 mg q12h",       60,        12, 28,
  "60 mg q24h",       60,        24, 14,
  "60 mg q48h",       60,        48,  7,
  "60 mg q72h",       60,        72,  5,
  "120 mg q24h",     120,        24, 14,
  "120 mg q48h",     120,        48,  7,
  "180 mg q48h",     180,        48,  7,
  "30 mg q24h",       30,        24, 14
)
regimen_absorbed <- function(dose, interval, n, adhere = TRUE) {
  m <- if (adhere) ida else ida %>% param(EMAX_ADH = 0)
  d <- m %>% mrgsim(events = oral_ev(dose, interval, n), end = 15 * 24,
                    delta = 24, hmax = 0.05) %>% as_tibble()
  b <- m %>% mrgsim(end = 15 * 24, delta = 24, hmax = 0.05) %>% as_tibble()
  tail(d$CUM_ABS, 1) - tail(b$CUM_ABS, 1)
}
noteC <- regimens %>%
  mutate(total_mg  = dose * n,
         absorbed  = mapply(regimen_absorbed, dose, interval, n,
                            MoreArgs = list(adhere = FALSE)),
         FIA_pct   = 100 * absorbed / total_mg,
         mg_per_day = absorbed / 14,
         per_dose   = absorbed / n)
##  expected (perfect adherence): FIA falls monotonically as the interval
##  shortens (20.9 / 20.1 / 17.6 / 13.9 % for q72/q48/q24/q12h) while
##  absorbed mg/day rises monotonically (4.48 / 6.03 / 10.55 / 16.72).

## =====================================================================
##  NOTE D — IV dose sweep: stores scale, the Hb clock does not
## =====================================================================
noteD <- lapply(c(0, 200, 500, 750, 1000, 1500, 2000, 3000), function(dose) {
  d <- if (dose == 0) mrgsim(ida, end = 12 * WK, delta = 12, hmax = 0.5)
       else mrgsim(ida, events = iv_ev(dose, 0), end = 12 * WK,
                   delta = 12, hmax = 0.5)
  d <- as_tibble(d); hb0 <- d$HB[1]
  slopes <- vapply(c(0, WK, 2 * WK, 3 * WK),
                   function(t) at_h(d, t + WK)$HB - at_h(d, t)$HB, numeric(1))
  tibble(dose_mg = dose,
         dHb_w1 = at_h(d, WK)$HB - hb0, dHb_w3 = at_h(d, 3 * WK)$HB - hb0,
         dHb_w12 = at_h(d, 12 * WK)$HB - hb0, max_weekly_slope = max(slopes),
         stores_w12 = at_h(d, 12 * WK)$STORES,
         ferritin_w12 = at_h(d, 12 * WK)$FERR)
}) %>% bind_rows()

## =====================================================================
##  NOTE E — inflammation: absorbed iron that never becomes haemoglobin
## =====================================================================
noteE <- lapply(c(0, 0.25, 0.5, 1, 2, 4), function(il) {
  o <- ida %>% param(IL6_IN = il) %>%
       mrgsim(events = oral_ev(65, 24, 84), end = 12 * WK, delta = 24, hmax = 0.5) %>%
       as_tibble()
  i <- ida %>% param(IL6_IN = il) %>%
       mrgsim(events = iv_ev(1000, 0), end = 12 * WK, delta = 24, hmax = 0.5) %>%
       as_tibble()
  w8 <- at_h(o, 8 * WK); w7 <- at_h(o, 8 * WK - WK)
  tibble(IL6_IN = il, IL6 = w8$IL6, hepcidin = w8$HEP, FPN_ENT = w8$FPN_ENT,
         absorbed_per_day = (w8$CUM_ABS - w7$CUM_ABS) / 7,
         loss_per_day     = (w8$CUM_LOSS - w7$CUM_LOSS) / 7,
         ferritin = w8$FERR, TSAT = w8$TSAT,
         dHb_oral = at_h(o, 12 * WK)$HB - o$HB[1],
         dHb_IV   = at_h(i, 12 * WK)$HB - i$HB[1])
}) %>% bind_rows()

## =====================================================================
##  NOTE F — IRIDA: TMPRSS6 loss of function, its own baseline
## =====================================================================
irida_base <- ida %>% param(TMPRSS6 = 0.30, VBLEED = 0.6) %>%
  mrgsim(end = 730 * 24, delta = 24, hmax = 1) %>% as_tibble() %>% tail(1)
##  Oral iron is nearly futile here and IV iron works, because the lesion
##  is inappropriately HIGH hepcidin rather than a depleted store.

## =====================================================================
##  NOTE G — the phosphate penalty of carboxymaltose
## =====================================================================
noteG <- lapply(list(
    list(l = "FCM 1000 mg",        ev = iv_ev(1000, 0),          p = list()),
    list(l = "FCM 750 mg x2",      ev = iv_ev(750, c(0, WK)),    p = list()),
    list(l = "FCM 500 mg",         ev = iv_ev(500, 0),           p = list()),
    list(l = "derisomaltose 1000", ev = iv_ev(1000, 0),          p = list(FCM_FGF = 0)),
    list(l = "sucrose 200 x5",     ev = iv_ev(200, WK * (0:4)),  p = list(FCM_FGF = 0.10))
  ), function(s) {
  m <- ida; if (length(s$p)) m <- do.call(param, c(list(m), s$p))
  d <- m %>% mrgsim(events = s$ev, end = 8 * WK, delta = 6, hmax = 0.5) %>% as_tibble()
  tibble(regimen = s$l, iFGF23_peak = max(d$FGF23),
         t_peak_d = d$time[which.max(d$FGF23)] / 24,
         PO4_nadir = min(d$PHOS), t_nadir_d = d$time[which.min(d$PHOS)] / 24,
         days_below_2.0 = sum(d$PHOS < 2.0) * 6 / 24,
         days_below_2.5 = sum(d$PHOS < 2.5) * 6 / 24,
         calcitriol_nadir = min(d$CTRIOL), PTH_peak = max(d$PTH))
}) %>% bind_rows()

## =====================================================================
##  NOTE H — tissue iron lags haemoglobin (fatigue and restless legs)
## =====================================================================
noteH <- lapply(c("1  no treatment", "2  oral 65 mg daily",
                  "5  oral 130 mg alternate", "8  IV FCM 1000 mg"), function(nm) {
  d <- sim_scenario(scenarios[[nm]])
  z <- d[1, ]; a <- at_h(d, 4 * WK); b <- at_h(d, 12 * WK)
  tibble(scenario = nm, Hb_w4 = a$HB, tissue_w4 = a$A_TISS, IRLS_w4 = a$IRLS,
         Hb_w12 = b$HB, tissue_w12 = b$A_TISS, IRLS_w12 = b$IRLS,
         FACIT_w12 = b$FACIT,
         tissue_repaired_pct = 100 * (b$A_TISS - z$A_TISS) / (386 - z$A_TISS),
         Hb_repaired_pct     = 100 * (b$HB - z$HB) / (13.86 - z$HB))
}) %>% bind_rows()

## =====================================================================
##  NOTE I — mass balance audit
## =====================================================================
##  Body iron must satisfy:
##     BODY_FE(t) = BODY_FE(0) + CUM_IV + CUM_ABS - CUM_LOSS
##  The enterocyte pools are deliberately EXCLUDED from BODY_FE: iron
##  taken up apically but shed with the enterocyte never entered the
##  body, so counting it would inflate both sides.
mass_balance <- function(events = NULL, end_h = 12 * WK) {
  d <- (if (is.null(events)) mrgsim(ida, end = end_h, delta = end_h, hmax = 0.5)
        else mrgsim(ida, events = events, end = end_h, delta = end_h, hmax = 0.5)) %>%
       as_tibble()
  z <- d[1, ]; e <- d[nrow(d), ]
  expected <- z$BODY_FE + (e$CUM_IV - z$CUM_IV) +
              (e$CUM_ABS - z$CUM_ABS) - (e$CUM_LOSS - z$CUM_LOSS)
  tibble(body_t0 = z$BODY_FE, body_end = e$BODY_FE, expected = expected,
         error_pct = 100 * (e$BODY_FE - expected) / expected)
}
noteI <- bind_rows(
  mutate(mass_balance(NULL),           run = "no treatment"),
  mutate(mass_balance(iv_ev(1000, 0)), run = "IV FCM 1000 mg"),
  mutate(mass_balance(oral_ev(65, 24, 84)), run = "oral 65 mg daily")
)

## =====================================================================
##  print everything
## =====================================================================
if (interactive() || !is.null(getOption("ida.report"))) {
  cat("\n=== 12-week scenarios ===\n");        print(as.data.frame(scenario_table), digits = 4)
  cat("\n=== NOTE B: refractory window ===\n");print(as.data.frame(noteB), digits = 4)
  cat("\n=== NOTE C: 14-day regimens ===\n");  print(as.data.frame(noteC), digits = 4)
  cat("\n=== NOTE D: IV dose sweep ===\n");    print(as.data.frame(noteD), digits = 4)
  cat("\n=== NOTE E: inflammation ===\n");     print(as.data.frame(noteE), digits = 4)
  cat("\n=== NOTE G: phosphate ===\n");        print(as.data.frame(noteG), digits = 4)
  cat("\n=== NOTE H: tissue iron lag ===\n");  print(as.data.frame(noteH), digits = 4)
  cat("\n=== NOTE I: mass balance ===\n");     print(as.data.frame(noteI), digits = 6)
}
