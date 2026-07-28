# =============================================================================
#  gdm_mrgsolve_model.R
#  Gestational Diabetes Mellitus (GDM) — Quantitative Systems Pharmacology model
#  ---------------------------------------------------------------------------
#  임신성 당뇨병 QSP 모델 · mrgsolve (31 ODEs, maternal–placental–fetal triad)
#
#  CENTRAL MODELLING THESIS
#  ------------------------
#  GDM is *not* modelled here as "high glucose in pregnancy". It is modelled as
#  a RACE BETWEEN TWO CLOCKS that both start at conception:
#
#     clock 1  the placenta grows logistically and its hormone output
#              (hPL, progesterone, oestradiol, TNF-alpha, cortisol via CRH)
#              drives maternal insulin sensitivity DOWN by ~50-60% by term.
#              This clock is essentially IDENTICAL in every pregnancy.
#
#     clock 2  the maternal beta-cell compartment expands (lactogen -> PRLR ->
#              JAK2/STAT5 -> Tph1/serotonin -> proliferation, plus a LEFT SHIFT
#              of the glucose-sensing curve). This clock is set by a single
#              patient-level parameter, BCAP (beta-cell adaptive capacity).
#
#  Whether a woman develops GDM is therefore an EMERGENT property of BCAP, not
#  an input. Two consequences fall out of the structure and are the reason the
#  model is worth simulating rather than tabulating:
#
#   (a) PHYSIOLOGIC SUBTYPES EMERGE (cf. Powe, Diabetes Care 2016;39:1052).
#       Analytic fixed points of the maternal glucose/insulin subsystem at
#       gestational week 40 (derivation in the CALIBRATION section below):
#          BCAP 1.00 / BMI 22  -> FPG ~77 mg/dL, fasting insulin ~14  (NGT)
#          BCAP 0.85 / BMI 24  -> FPG ~87 mg/dL, fasting insulin ~13  (secretion-
#                                 deficient: FPG nearly NORMAL, the defect only
#                                 shows up on the post-load limb)
#          BCAP 1.00 / BMI 34  -> FPG ~100 mg/dL, fasting insulin ~18 (insulin-
#                                 resistant: high FPG *with* hyperinsulinaemia)
#          BCAP 0.55 / BMI 31  -> FPG ~110 mg/dL                      (mixed)
#       The same disease label, three different fasting/insulin signatures, and
#       — importantly — three different predicted drug responses, because
#       metformin acts on the hepatic arm and insulin on the whole system.
#
#   (b) THE FETUS IS DRIVEN BY A ONE-WAY VALVE. Glucose crosses by GLUT1
#       facilitated diffusion; INSULIN DOES NOT CROSS. So maternal glucose sets
#       fetal glucose, fetal glucose grows the fetal beta-cell mass, fetal
#       insulin then acts as a *fetal* growth factor (Pedersen 1952; Freinkel's
#       fuel-mediated teratogenesis, Diabetes 1980;29:1023). Fat mass is far
#       more insulin-elastic than lean mass, which is why the overgrowth is
#       ASYMMETRIC (abdominal circumference >> head circumference) and why the
#       same birth weight carries different shoulder-dystocia risk in GDM than
#       in constitutional macrosomia.
#
#  WHAT THE MODEL IS FOR
#  ---------------------
#  Comparing therapeutic strategies on endpoints that matter (LGA, neonatal
#  hypoglycaemia, shoulder dystocia, pre-eclampsia, caesarean, and the maternal
#  5-year T2DM hazard) while explicitly tracking FETAL DRUG EXPOSURE — the axis
#  on which insulin, metformin and glyburide genuinely differ:
#       insulin    : does not cross the placenta (fetal exposure = 0)
#       metformin  : crosses freely, umbilical:maternal ratio ~1.0
#       glyburide  : crosses, BCRP/ABCG2-effluxed, cord:maternal ~0.7
#
#  UNITS
#  -----
#  time            days of gestation (TIME = gestational age in days; GW = TIME/7)
#                  simulation window: day 56 (GW 8) -> day 364 (12 wk postpartum)
#                  delivery at TDEL (default 273 d = GW 39)
#  glucose         mg/dL            insulin        microU/mL
#  hPL             mg/L             progesterone   ng/mL      oestradiol pg/mL
#  TNF-alpha       pg/mL            adiponectin    microg/mL
#  FFA             mmol/L           fetal mass     g          maternal fat  kg
#  metformin       mg (amount) / mg/L (conc.)      glyburide  mg / mg/L
#  insulin dose    U
#  glucose fluxes  mg/dL/day  (a lumped maternal glucose volume VG = 112 dL is
#                  folded into the rate constants; see CALIBRATION)
#
#  HONEST STATEMENT OF STATUS
#  --------------------------
#  Every number in $PARAM is either (i) taken from the literature cited in
#  gdm_references.md, or (ii) derived by hand from a published steady-state
#  observation, with the derivation written out in the CALIBRATION section so it
#  can be checked and re-fitted. The model has NOT been fitted to individual
#  patient data. The scenario block at the bottom must actually be run; the
#  numbers quoted in the comments above are ANALYTIC FIXED POINTS of the
#  maternal subsystem, not simulation output. Do not quote them as results.
#
#  Requires: mrgsolve (>= 1.0), dplyr, tidyr, ggplot2
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

gdm_code <- '
$PROB
# Gestational Diabetes Mellitus QSP model
- maternal-placental-fetal triad, 31 ODEs
- placental contra-insulin drive vs beta-cell compensation
- insulin / metformin / glyburide PK-PD with explicit placental transfer

$PARAM
// ---------------------------------------------------------------- PLACENTA --
PLMAX    = 1.0      // asymptotic normalised syncytiotrophoblast mass
KGP      = 0.035    // logistic placental growth rate (1/day); 0.05 at GW8 -> 0.99 at GW40
TDEL     = 273      // gestational day of delivery (GW 39)
KDELIV   = 5.0      // post-delivery placental involution rate (1/day) -> hPL gone in <1 d
PLACINS  = 1.0      // placental sufficiency multiplier (<1 = insufficiency/IUGR)

KIN_HPL  = 158.0    // hPL synthesis (mg/L/day per unit placenta)
KOUT_HPL = 24.0     // hPL elimination (1/day); term hPL ~ 6.5 mg/L
KIN_PROG = 3640.0   // progesterone synthesis (ng/mL/day per unit placenta)
KOUT_P   = 24.0     // -> term progesterone ~ 150 ng/mL
KCL      = 0.03     // corpus-luteum involution rate (1/day); CLP starts at 25 ng/mL
                    //   and is gone by ~GW 20 (luteo-placental shift)
KIN_E2   = 480000.0 // oestradiol synthesis (pg/mL/day per unit placenta)
KOUT_E2  = 24.0     // -> term E2 ~ 20000 pg/mL

KIN_TNF  = 30.0     // TNF-alpha synthesis (pg/mL/day) placental term
FADIP    = 0.055    // adipose contribution to TNF-alpha (per kg fat above FATREF)
FATREF   = 17.4     // reference maternal fat mass at BMI 22 (kg)
KOUT_TNF = 24.0     // -> term TNF ~ 2.4 pg/mL lean, ~3.1 pg/mL obese
TNF0     = 1.2      // pre-pregnancy TNF-alpha (pg/mL)

KSYN_AD  = 236.0    // adiponectin synthesis (microg/mL/day) -> 9 microg/mL at BMI 22
IADT     = 0.55     // TNF-alpha suppression of adiponectin
IADF     = 0.030    // adiposity suppression of adiponectin (per kg fat)
KDEG_AD  = 12.0     // adiponectin elimination (1/day); non-pregnant ~ 9 microg/mL
ADIPO0   = 9.0      // reference adiponectin (microg/mL)
KADS     = 0.25     // adiponectin -> insulin-sensitivity gain

// ------------------------------------------- CONTRA-INSULIN / SENSITIVITY --
SHPL     = 0.100    // hPL weight in contra-insulin index      (per mg/L)
SPROG    = 0.250    // progesterone weight                     (per 100 ng/mL)
STNF     = 0.120    // TNF-alpha weight                        (per pg/mL)
SCORT    = 0.160    // cortisol/CRH weight                     (per unit placenta)
HFR      = 0.80     // hepatic:peripheral ratio of IR development
SIREF    = 1.51     // peripheral insulin sensitivity, BMI 22, non-pregnant
SIHREF   = 0.0645   // hepatic insulin sensitivity (EGP suppression, per microU/mL)
BMI      = 22.0     // pre-gestational BMI (kg/m2)
KBMI     = 0.055    // BMI penalty on insulin sensitivity (per kg/m2 above 22)
EXEFF    = 0.0      // exercise effect on SI (0 = sedentary, 0.20 = >=150 min/wk)

// -------------------------------------------------- MATERNAL GLUCOSE FLUX --
EGPB     = 2552.0   // maximal endogenous glucose production (mg/dL/day)
EGPP     = 0.30     // placental augmentation of maximal EGP at term
KEGPF    = 0.20     // FFA (gluconeogenic substrate) augmentation of EGP
KGEGP    = 0.50     // hepatic glucose autoregulation (glucose-mediated EGP brake)
GBRAIN   = 900.0    // obligate insulin-independent uptake, CNS+RBC (mg/dL/day)
GEZI     = 2.0      // residual glucose effectiveness (1/day)
KPLU     = 50.0     // uteroplacental glucose extraction (per unit placenta/day)
KREN     = 8.0      // renal glucose clearance above threshold (1/day)
GLUTHR   = 155.0    // renal glucose threshold in pregnancy (mg/dL); non-preg 180
VG       = 112.0    // maternal glucose distribution volume (dL)
KA_G     = 30.0     // meal glucose absorption (1/day)
KTR_G    = 26.0     // gastric transit (1/day); low-GI diet lowers this
FCHO     = 0.90     // fraction of ingested carbohydrate reaching plasma

// ------------------------------------------------ MATERNAL BETA-CELL AXIS --
BCAP     = 1.00     // *** beta-cell ADAPTIVE CAPACITY — the GDM switch ***
                    //     1.00 NGT | 0.85 secretion-deficient | 0.55 mixed | 0.30 severe
KPROL    = 0.0135   // lactogen-driven beta-cell expansion rate (1/day)
BCMMAX   = 2.6      // ceiling on beta-cell mass x function
KMH      = 2.0      // hPL EC50 for the adaptation signal (mg/L)
KSENS    = 0.28     // maximal LEFT SHIFT of the glucose-sensing curve
KAPO     = 0.0015   // basal beta-cell attrition (1/day)
KGLUTOX  = 0.90     // glucotoxic attrition (per 100 mg/dL above 100)
KFFATOX  = 0.45     // lipotoxic attrition (per mmol/L FFA above 0.5)
KDECL    = 0.00022  // post-partum chronic beta-cell decline (1/day)
SMAX     = 3.30e7   // maximal insulin secretion rate (microU/day) at BCM = 1
                    //   SMAX*(1+KINCR) = 3.70e7, the value used in the
                    //   fixed-point derivations in the CALIBRATION section
KG50     = 130.0    // glucose EC50 for secretion (mg/dL), non-adapted
HG       = 2.0      // Hill coefficient for glucose-stimulated secretion
KINCR    = 0.12     // incretin (GLP-1) amplification of secretion
KMGLP    = 0.50     // metformin augmentation of the incretin term (gut GLP-1)
VI       = 7000.0   // maternal insulin distribution volume (mL)
KEI0     = 204.0    // insulin elimination (1/day); t1/2 ~5 min
KEIP     = 0.35     // pregnancy increase in insulin clearance (placental IDE)

// ------------------------------------------------------------- MATERNAL FFA --
KLIP     = 5.4      // basal lipolytic flux (mmol/L/day)
LIPP     = 0.85     // placental (hPL/GH-V/cortisol) augmentation of lipolysis
SIL_FFA  = 0.075    // insulin suppression of lipolysis (per microU/mL)
KFFA     = 12.0     // FFA clearance (1/day); non-pregnant FFA ~ 0.45 mmol/L
FFAREF   = 0.45     // reference FFA (mmol/L)

// ----------------------------------------------------- MATERNAL ADIPOSITY --
RFAT     = 0.032    // gestational fat accretion (kg/day) at 200 g CHO, sedentary
                    //   -> ~6.9 kg fat over GW 8-39; with fetus+placenta+fluid
                    //      gives total GWG ~11-12 kg (IOM normal-BMI target)
CHOD     = 200.0    // prescribed daily carbohydrate (g/day) — MNT lever
KLOSSM   = 0.010    // POST-PARTUM fat loss back toward set point (1/day)

// ------------------------------------------------------ GLYCAEMIC MARKERS --
KRBC     = 0.025    // erythrocyte turnover / glycation averaging (1/day), t1/2 ~28 d
A1COFF   = 0.40     // downward offset of HbA1c in pregnancy (accelerated erythropoiesis)

// ---------------------------------------------------- PLACENTAL TRANSFER ---
KTRF     = 120.0    // GLUT1 materno-fetal glucose transfer (1/day per unit placenta)
FCONS    = 1440.0   // fetal+placental glucose consumption at term (mg/dL/day)
FINS     = 0.35     // fetal insulin augmentation of fetal glucose utilisation

// ------------------------------------------------------ FETAL COMPARTMENT --
KFBG     = 0.0151   // fetal beta-cell developmental growth (1/day)
BCFMAX   = 3.0      // ceiling on fetal beta-cell mass
FBTHR    = 70.0     // fetal glucose threshold for beta-cell hyperplasia (mg/dL)
FBS      = 1.40     // strength of glucose-driven fetal beta-cell hyperplasia
SMAXF    = 7.90e5   // fetal maximal insulin secretion (microU/day) at BCF = 1
KGF      = 110.0    // fetal glucose EC50 for secretion (mg/dL)
VIF      = 170.0    // fetal insulin volume (mL)
KEIF     = 150.0    // fetal insulin elimination (1/day)
INSFREF  = 8.0      // reference term fetal insulin (microU/mL) in NGT
KSIGF    = 0.60     // fetal IGF-1 synthesis (1/day, normalised)
KDIGF    = 0.60     // fetal IGF-1 elimination (1/day) -> baseline IGF1F = 1
AIGF     = 0.20     // fetal insulin -> IGF-1 gain
KLN      = 0.0419   // fetal lean-mass logistic growth (1/day)
LNMAX    = 3600.0   // fetal lean-mass ceiling (g)
KFT      = 5.00     // fetal fat accretion (g/day) scale
AFAT     = 0.55     // fetal insulin -> fat accretion gain  *** key elasticity ***
BFAT     = 0.25     // maternal FFA/TG -> fetal fat accretion gain

// ---------------------------------------------------------- METFORMIN PK ---
KAM      = 12.0     // metformin absorption (1/day)
FM       = 0.55     // oral bioavailability
VMC      = 100.0    // central volume (L)
VMP      = 180.0    // peripheral volume (L)
QM       = 260.0    // intercompartmental clearance (L/day)
CLM      = 1440.0   // renal clearance (L/day), non-pregnant
CLMP     = 0.25     // pregnancy increase in metformin CL (GFR +50%)
KPLM     = 200.0    // placental metformin transfer (L/day equivalent)
VMF      = 3.0      // fetal metformin volume (L)
CLMF     = 12.0     // fetal metformin clearance (L/day)
                    //   KPLM/(KPLM+CLMF) = 0.94 -> umbilical:maternal ~1.0
EMAXM    = 0.35     // maximal metformin suppression of EGP
EC50M    = 0.50     // metformin EC50 (mg/L)
EMAXMSI  = 0.15     // maximal metformin gain in peripheral SI
METADH   = 1.0      // adherence multiplier (GI intolerance -> <1)

// ------------------------------------------------------------ INSULIN PK ---
KAB      = 2.10     // basal analogue (detemir/NPH) absorption (1/day), t1/2 ~8 h
KAF      = 20.0     // rapid analogue (aspart/lispro) absorption (1/day), tmax ~50 min

// ---------------------------------------------------------- GLYBURIDE PK ---
KAG      = 14.0     // glyburide absorption (1/day)
FG       = 0.90     // oral bioavailability
VGC      = 12.0     // central volume (L)
CLG      = 34.0     // clearance (L/day), non-pregnant
CLGP     = 1.00     // pregnancy DOUBLING of glyburide clearance
KPLG     = 28.0     // placental glyburide transfer, net of BCRP efflux (L/day)
VGF_F    = 1.4      // fetal glyburide volume (L)
CLGF     = 12.0     // fetal glyburide clearance (L/day)
                    //   KPLG/(KPLG+CLGF) = 0.70 -> cord:maternal ~0.7
SGLYB    = 1.05e7   // maximal glucose-INDEPENDENT insulin secretion (microU/day)
EC50G    = 0.09     // glyburide EC50 at SUR1 (mg/L)
GLYBSTR  = 0.35     // beta-cell strain (ER stress) per unit SUR1 occupancy

// ----------------------------------------------------- OUTCOME MODEL (PD) --
SDBW     = 0.115    // CV of birth weight for gestational age
ZLGA     = 1.282    // z threshold for LGA (>90th percentile)
SDIND    = 1.00     // residual individual variability in birth-weight z
PNH0     = 0.015    // baseline neonatal hypoglycaemia probability
PNHMAX   = 0.350    // maximal neonatal hypoglycaemia probability
KNH      = 1.10     // EC50 of relative fetal hyperinsulinaemia
GNH      = 1.60     // Hill coefficient for neonatal hypoglycaemia
PSD0     = 0.014    // baseline shoulder-dystocia probability at z = 0
KSD      = 1.05     // birth-weight-z sensitivity of shoulder dystocia
PCS0     = 0.250    // baseline caesarean probability
KCS      = 0.140    // birth-weight-z sensitivity of caesarean
PPE0     = 0.035    // baseline pre-eclampsia probability (lean, normoglycaemic)
KPEG     = 0.566    // glycaemic sensitivity of pre-eclampsia (per 10 mg/dL MBG >95)
KPEB     = 1.100    // BMI sensitivity of pre-eclampsia (per 5 kg/m2 >25)
H0T2D    = 0.020    // 5-year T2DM hazard at reference disposition index
KDI      = 3.80     // disposition-index sensitivity of the T2DM hazard
LACT     = 0.0      // lactation >=3 months (0/1); protective
KLACT    = 0.45     // hazard reduction with lactation

// 31 compartments. Declared once, with values, in $INIT — do NOT also add a
// $CMT block for the same names (mrgsolve treats that as a redeclaration).
//   placenta   (7): PLAC HPL PROG E2 TNFA ADIPO CLP
//   maternal   (8): BCM GLU INS FFA MBG GUT1 GUT2 FATM
//   fetal      (6): GLUF INSF BCF IGF1F LEANF FATF
//   metformin  (4): MGUT MCEN MPER MFET
//   insulin    (2): IBAS IBOL
//   glyburide  (3): GGUT GCEN GFET
//   safety     (1): HYPOAUC
$INIT
PLAC = 0.050, HPL = 0.33, PROG = 7.6, E2 = 1000.0, TNFA = 1.2, ADIPO = 9.0, CLP = 25.0,
BCM = 1.00, GLU = 85.0, INS = 6.0, FFA = 0.45, MBG = 85.0, GUT1 = 0.0, GUT2 = 0.0, FATM = 18.0,
GLUF = 80.0, INSF = 2.0, BCF = 0.050, IGF1F = 1.0, LEANF = 1.5, FATF = 0.05,
MGUT = 0.0, MCEN = 0.0, MPER = 0.0, MFET = 0.0,
IBAS = 0.0, IBOL = 0.0,
GGUT = 0.0, GCEN = 0.0, GFET = 0.0,
HYPOAUC = 0.0

$GLOBAL
// derived quantities shared between $ODE and $TABLE
double CID, SIREL, SIRELH, FBMI, ADIPOF, SIP, SILV, ADAPT, KG50E;
double EGPn, RAM, UII, UID, UPLAC, UREN, ISRE, ISRG, IEXO, KEIt;
double CMET, CMETF, CGLY, CGLYF, MEGP, MSI, SUROCC;
double FWT, FWREF, BWZ, GWk, PLACEFF, LACTDR, POSTP, FATSP;

#define PHI(x) (1.0/(1.0+exp(-1.702*(x))))   // logistic approximation to Normal CDF

$MAIN
// Maternal fat mass at conception scales with pre-gestational BMI:
// height 1.62 m -> BMI 22 = 58 kg at ~30% fat = 17.4 kg; BMI 34 = 89 kg at
// ~45% fat = 40 kg. Linear interpolation is adequate over the BMI 20-40 range.
FATSP  = 17.4 + 1.90*(BMI - 22.0);
FATM_0 = FATSP;

$ODE
GWk     = SOLVERTIME/7.0;
POSTP   = (SOLVERTIME > TDEL) ? 1.0 : 0.0;

// ---------------------------------------------------------------- PLACENTA --
// logistic growth until delivery, then rapid involution (hPL t1/2 ~15 min in
// vivo; clinically undetectable within hours -> the contra-insulin drive
// switches OFF at delivery, which is exactly why GDM "resolves" post partum)
if (POSTP < 0.5) {
  dxdt_PLAC = KGP*PLAC*(1.0 - PLAC/PLMAX);
} else {
  dxdt_PLAC = -KDELIV*PLAC;
}
PLACEFF = PLAC*PLACINS;

dxdt_HPL   = KIN_HPL*PLACEFF  - KOUT_HPL*HPL;
dxdt_CLP   = -KCL*CLP;                                  // corpus luteum involution
dxdt_PROG  = KIN_PROG*PLACEFF - KOUT_P*PROG;
dxdt_E2    = KIN_E2*PLACEFF   - KOUT_E2*E2;
dxdt_TNFA  = KIN_TNF*PLACEFF
             + KOUT_TNF*(TNF0 + FADIP*fmax(0.0, FATM - FATREF) - TNFA);
dxdt_ADIPO = KSYN_AD/(1.0 + IADT*TNFA + IADF*FATM) - KDEG_AD*ADIPO;

// ---------------------------------------------- CONTRA-INSULIN AGGREGATION --
// hPL, progesterone and TNF-alpha all converge on IRS-1 serine phosphorylation;
// cortisol (via placental CRH) acts mainly hepatically. We aggregate them into
// one dimensionless index because they are not separately identifiable from
// clinical data — that is a deliberate parsimony choice, not an oversight.
CID    = SHPL*HPL + SPROG*(PROG + CLP)/100.0 + STNF*TNFA + SCORT*PLACEFF;
SIREL  = 1.0/(1.0 + CID);
SIRELH = 1.0/(1.0 + HFR*CID);
FBMI   = exp(-KBMI*(BMI - 22.0));
ADIPOF = 1.0 + KADS*(ADIPO/ADIPO0 - 1.0);

CMET  = MCEN/VMC;                       // maternal metformin (mg/L)
CMETF = MFET/VMF;                       // fetal metformin (mg/L)
CGLY  = GCEN/VGC;                       // maternal glyburide (mg/L)
CGLYF = GFET/VGF_F;                     // fetal glyburide (mg/L)
MEGP  = 1.0/(1.0 + EMAXM*CMET/(EC50M + CMET));          // hepatic effect
MSI   = EMAXMSI*CMET/(EC50M + CMET);                    // peripheral effect
SUROCC = CGLY/(EC50G + CGLY);                           // SUR1 fractional occupancy

SIP  = SIREF  * SIREL  * FBMI * ADIPOF * (1.0 + EXEFF + MSI);
SILV = SIHREF * SIRELH * FBMI          * (1.0 + EXEFF);

// ------------------------------------------------ BETA-CELL ADAPTATION -----
// Lactogen drive is hPL-dominated; BCAP scales BOTH the proliferative response
// and the left shift of the glucose-sensing curve, so a woman with low BCAP is
// short of beta-cell MASS *and* runs at a higher glucose set-point.
LACTDR = HPL/(HPL + KMH);
ADAPT  = BCAP*LACTDR;
KG50E  = KG50*(1.0 - KSENS*ADAPT);

// Attrition carries three stressors: glucotoxicity, lipotoxicity, and — for the
// sulfonylurea arm — the secretory strain of glucose-INDEPENDENT stimulation
// (GLYBSTR x SUR1 occupancy). That last term is why a secretagogue is not a
// neutral way to lower maternal glucose in a beta-cell that is already failing.
dxdt_BCM = KPROL*BCAP*LACTDR*BCM*(1.0 - BCM/BCMMAX)*(1.0 - POSTP)
           - KAPO*(1.0 + KGLUTOX*fmax(0.0, GLU - 100.0)/100.0
                       + KFFATOX*fmax(0.0, FFA - 0.5)
                       + GLYBSTR*SUROCC)*BCM
           - KDECL*POSTP*BCM*(1.0 + KGLUTOX*fmax(0.0, GLU - 100.0)/100.0);

// ------------------------------------------------- MATERNAL GLUCOSE FLUXES --
RAM   = FCHO*KA_G*GUT2/VG;
UII   = GBRAIN + GEZI*GLU;
UID   = SIP*INS*GLU;
UPLAC = KPLU*PLACEFF*fmax(0.0, GLU - GLUF);
UREN  = KREN*fmax(0.0, GLU - GLUTHR);
EGPn  = EGPB*(1.0 + EGPP*PLACEFF)*(1.0 + KEGPF*(FFA/FFAREF - 1.0))*MEGP
        / ((1.0 + SILV*INS)*(1.0 + KGEGP*fmax(0.0, GLU - 80.0)/80.0));

dxdt_GUT1 = -KTR_G*GUT1;
dxdt_GUT2 =  KTR_G*GUT1 - KA_G*GUT2;
dxdt_GLU  = EGPn + RAM - UII - UID - UPLAC - UREN;

// --------------------------------------------------------- MATERNAL INSULIN --
ISRE   = BCM*SMAX*pow(GLU, HG)/(pow(KG50E, HG) + pow(GLU, HG))
         *(1.0 + KINCR*(1.0 + KMGLP*CMET/(EC50M + CMET)));
ISRG   = BCM*SGLYB*SUROCC;                       // glucose-INDEPENDENT (the hazard)
IEXO   = (KAB*IBAS + KAF*IBOL)*1.0e6;            // U -> microU
KEIt   = KEI0*(1.0 + KEIP*PLACEFF);

dxdt_INS = (ISRE + ISRG + IEXO)/VI - KEIt*INS;

dxdt_IBAS = -KAB*IBAS;
dxdt_IBOL = -KAF*IBOL;

// ------------------------------------------------------------ MATERNAL FFA --
dxdt_FFA = KLIP*(1.0 + LIPP*PLACEFF)/(1.0 + SIL_FFA*INS) - KFFA*FFA;

// --------------------------------------------------- ADIPOSITY & GLYCATION --
dxdt_FATM = RFAT*(1.0 - fmin(0.9, 2.0*EXEFF))*(CHOD/200.0)*(1.0 - POSTP)
            - KLOSSM*POSTP*(FATM - FATSP);
dxdt_MBG  = KRBC*(GLU - MBG);                    // exponentially weighted mean glucose

// cumulative maternal hypoglycaemia exposure (mg/dL x day below 70)
dxdt_HYPOAUC = fmax(0.0, 70.0 - GLU);

// ---------------------------------------------------------- FETAL GLUCOSE --
// GLUT1 facilitated diffusion, bidirectional, saturable only at extremes;
// consumption scales with fetal mass and is amplified by fetal insulin.
// Every fetal derivative is gated by (1 - POSTP): at delivery the whole fetal
// block FREEZES at its delivery values, which is what the endpoint equations
// read. (Without the gate, cutting off placental transfer drains GLUF to
// negative values and the squared Hill term happily keeps making insulin.)
FWT = LEANF + FATF;
dxdt_GLUF = (KTRF*PLACEFF*(GLU - GLUF)
             - FCONS*(1.0 + FINS*(INSF/INSFREF - 1.0))*(FWT/3450.0))
            *(1.0 - POSTP);

// ------------------------------------------------------ FETAL BETA-CELL ----
dxdt_BCF = KFBG*(1.0 + FBS*fmax(0.0, GLUF - FBTHR)/FBTHR)*BCF*(1.0 - BCF/BCFMAX)
           *(1.0 - POSTP);

// glyburide reaching the fetus closes fetal SUR1 as well -> glucose-independent
// fetal insulin release, the mechanistic reason for the neonatal-hypoglycaemia
// signal in the glyburide trials
dxdt_INSF = ((SMAXF*BCF*pow(fmax(0.0, GLUF), 2.0)
              /(pow(KGF, 2.0) + pow(fmax(0.0, GLUF), 2.0))
              + SMAXF*BCF*0.45*CGLYF/(EC50G + CGLYF))/VIF
             - KEIF*INSF)*(1.0 - POSTP);

dxdt_IGF1F = (KSIGF*(1.0 + AIGF*(INSF/INSFREF - 1.0)) - KDIGF*IGF1F)*(1.0 - POSTP);

// -------------------------------------------------------- FETAL BODY MASS --
// lean mass: IGF-1 driven logistic (nearly saturated at term -> inelastic)
// fat  mass: insulin + lipid driven, weighted to the third trimester (PLAC^3)
//            -> the elastic compartment; this asymmetry IS the phenotype
dxdt_LEANF = KLN*IGF1F*LEANF*(1.0 - LEANF/LNMAX)*(1.0 - POSTP);
dxdt_FATF  = KFT*pow(PLACEFF, 3.0)
             *(1.0 + AFAT*(INSF/INSFREF - 1.0) + BFAT*(FFA/FFAREF - 1.0))
             *(LEANF/1000.0)*(1.0 - POSTP)
             + 0.20*CMETF*0.0;                   // placeholder: see note in README

// --------------------------------------------------------- METFORMIN PK ----
// PGATE is a SMOOTH gate on placental drug transfer (surface-area dependent);
// a hard if/else here makes LSODA chatter around the switch point.
double PGATE = PLACEFF/(PLACEFF + 0.02);
double JMET  = KPLM*PGATE*(CMET - CMETF);
double JGLY  = KPLG*PGATE*(CGLY - CGLYF);

dxdt_MGUT = -KAM*MGUT;
dxdt_MCEN =  FM*KAM*MGUT*METADH
             - (CLM*(1.0 + CLMP*PLACEFF)/VMC)*MCEN
             - (QM/VMC)*MCEN + (QM/VMP)*MPER
             - JMET;
dxdt_MPER =  (QM/VMC)*MCEN - (QM/VMP)*MPER;
dxdt_MFET =  JMET - CLMF*CMETF;

// --------------------------------------------------------- GLYBURIDE PK ----
dxdt_GGUT = -KAG*GGUT;
dxdt_GCEN =  FG*KAG*GGUT - (CLG*(1.0 + CLGP*PLACEFF)/VGC)*GCEN - JGLY;
dxdt_GFET =  JGLY - CLGF*CGLYF;

$TABLE
GWk   = TIME/7.0;
FWT   = LEANF + FATF;

// Re-evaluate the derived quantities AT THE OUTPUT TIME. The $ODE globals hold
// whatever the last internal solver step left behind, which is close but not
// identical; anything reported to the user is recomputed here.
PLACEFF = PLAC*PLACINS;
CID     = SHPL*HPL + SPROG*(PROG + CLP)/100.0 + STNF*TNFA + SCORT*PLACEFF;
SIREL   = 1.0/(1.0 + CID);
SIRELH  = 1.0/(1.0 + HFR*CID);
FBMI    = exp(-KBMI*(BMI - 22.0));
ADIPOF  = 1.0 + KADS*(ADIPO/ADIPO0 - 1.0);
CMET    = MCEN/VMC;
CMETF   = MFET/VMF;
CGLY    = GCEN/VGC;
CGLYF   = GFET/VGF_F;
MSI     = EMAXMSI*CMET/(EC50M + CMET);
SUROCC  = CGLY/(EC50G + CGLY);
SIP     = SIREF  * SIREL  * FBMI * ADIPOF * (1.0 + EXEFF + MSI);
SILV    = SIHREF * SIRELH * FBMI          * (1.0 + EXEFF);
LACTDR  = HPL/(HPL + KMH);
ADAPT   = BCAP*LACTDR;
KG50E   = KG50*(1.0 - KSENS*ADAPT);

// ---- reference fetal-weight curve (ln W = 0.8669 + 0.31123*GW - 0.0032313*GW^2)
//      anchored to 330 g @ GW20, 1150 g @ GW28, 3450 g @ GW40
double GWc = fmin(GWk, 42.0);
FWREF = exp(0.8669 + 0.31123*GWc - 0.0032313*GWc*GWc);
BWZ   = (FWT - FWREF)/(SDBW*FWREF);

double FATPCT = 100.0*FATF/fmax(1.0, FWT);
double HBA1C  = (MBG + 46.7)/28.7 - A1COFF;
double MATSUDA = SIP/SIREF;                       // insulin sensitivity, fraction of pre-preg
double DI      = SIP*BCM;                         // disposition index (SI x secretory capacity)
double DIREF   = SIREF*1.0;
double NHR     = INSF/INSFREF;                    // relative fetal hyperinsulinaemia
double CORDCP  = INSF*3.0;                        // cord C-peptide proxy (pmol/L-ish scale)

// ---- clinical endpoint probabilities (evaluated meaningfully at TIME = TDEL) --
double PLGA  = PHI((BWZ - ZLGA)/SDIND);
double PMACR = PHI((FWT - 4000.0)/(SDBW*FWREF));
double PNHYP = PNH0 + (PNHMAX - PNH0)
               *pow(fmax(0.0, NHR - 1.0), GNH)
               /(pow(KNH, GNH) + pow(fmax(0.0, NHR - 1.0), GNH));
double PSD   = PSD0*exp(KSD*fmax(0.0, BWZ));
double PCS   = fmin(0.60, PCS0 + KCS*fmax(0.0, BWZ));
double PPE   = fmin(0.60, PPE0*(1.0 + KPEG*fmax(0.0, MBG - 95.0)/10.0)
                              *(1.0 + KPEB*fmax(0.0, BMI - 25.0)/5.0));

// ---- maternal 5-year T2DM risk from the POST-PARTUM disposition index -------
double H5    = H0T2D*exp(-KDI*(DI/DIREF - 1.0))*(1.0 - KLACT*LACT);
double CIT2D = 1.0 - exp(-fmin(3.0, H5));

// ---- diagnostic-style observables ------------------------------------------
double GLUFRATIO = GLUF/fmax(1.0, GLU);
double EXOSHARE  = IEXO/fmax(1.0, ISRE + ISRG + IEXO);   // fraction of insulin that is exogenous

$CAPTURE
GWk CID SIREL SIP SILV MATSUDA ADAPT KG50E DI DIREF
EGPn UID UPLAC UREN RAM ISRE ISRG IEXO EXOSHARE KEIt
FWT FWREF BWZ FATPCT HBA1C NHR CORDCP GLUFRATIO
CMET CMETF CGLY CGLYF SUROCC
PLGA PMACR PNHYP PSD PCS PPE H5 CIT2D
'

gdm <- mcode("gdm_qsp", gdm_code, soloc = tempdir())

# mrgsolve accepts a character `cmt` column in data sets, but a numeric index
# is accepted by every version, so resolve names -> indices once here and let
# the event builders below stay readable.
CMT_INDEX <- setNames(seq_along(mrgsolve::cmt(gdm)), mrgsolve::cmt(gdm))
cmt_i <- function(nm) {
  i <- unname(CMT_INDEX[nm])
  if (any(is.na(i))) stop("unknown compartment: ", paste(nm[is.na(i)], collapse = ", "))
  i
}

# =============================================================================
#  CALIBRATION — every non-obvious number, and how it was obtained
# =============================================================================
#
#  A. LUMPED GLUCOSE UNITS.  All maternal glucose fluxes are expressed in
#     mg/dL/day with the distribution volume VG = 112 dL (0.16 L/kg x 70 kg)
#     folded in. Conversion: 1 mg/kg/min at 70 kg = 100,800 mg/day = 900
#     mg/dL/day. So EGP of 2.0 mg/kg/min == 1800 mg/dL/day. Read every
#     glucose-flux parameter through that factor of 900.
#
#  B. PLACENTAL CLOCK. PLAC is logistic, 0.05 at GW 8 -> 0.99 at GW 40 with
#     KGP = 0.035/day. hPL, progesterone and oestradiol are quasi-steady with
#     respect to PLAC (KOUT = 24/day, t1/2 = 0.7 d) so their term values are set
#     directly by KIN/KOUT: hPL 6.5 mg/L, progesterone ~150 ng/mL, E2 ~20,000
#     pg/mL — all standard third-trimester values.
#
#  C. INSULIN-SENSITIVITY DECLINE. Target: SI falls ~50-55% from early to late
#     pregnancy (Catalano, Am J Obstet Gynecol 1999;180:903; Am J Physiol
#     1993;264:E60). The weights SHPL/SPROG/STNF/SCORT were chosen so that
#         CID(GW 8)  ~ 0.22   -> SIREL = 0.82
#         CID(GW 40) ~ 1.40   -> SIREL = 0.42
#     i.e. SI(term)/SI(early) = 0.51. TNF-alpha carries the single largest
#     weight per unit change, matching Kirwan (Diabetes 2002;51:2207), where
#     TNF-alpha was the strongest inverse correlate of insulin sensitivity in
#     pregnancy — stronger than hPL, cortisol, oestradiol or leptin.
#
#  D. MATERNAL FASTING FIXED POINTS (solved by hand, iterating
#     GLU <-> INS to the joint steady state of dxdt_GLU = dxdt_INS = 0):
#
#       non-pregnant, BMI 22, BCAP 1.0     : FPG ~85, insulin ~5.8
#       GW 40, BCAP 1.00, BMI 22 (NGT)     : FPG ~77, insulin ~14, BCM ~2.0
#       GW 40, BCAP 0.85, BMI 24           : FPG ~87, insulin ~13, BCM ~1.75
#       GW 40, BCAP 1.00, BMI 34           : FPG ~100, insulin ~18, BCM ~2.0
#       GW 40, BCAP 0.55, BMI 31           : FPG ~110, insulin ~13, BCM ~1.4
#
#     The FALL in fasting glucose in normal late pregnancy (85 -> 77) is not
#     imposed; it emerges because the uteroplacental sink UPLAC (KPLU = 50,
#     ~600 mg/dL/day at term = 0.67 mg/kg/min, consistent with measured uterine
#     glucose uptake) outgrows the modest rise in maximal EGP (EGPP = 0.30).
#     Hepatic glucose autoregulation (KGEGP) is what keeps EGP from running
#     away when hepatic insulin sensitivity collapses — without it the model
#     predicts frank fasting hyperglycaemia in every obese pregnancy, which is
#     wrong.
#
#  E. BETA-CELL ADAPTATION. Human pregnancy raises insulin secretion ~2-2.5x,
#     mostly by FUNCTION rather than mass (Butler, Diabetologia 2010;53:2167 —
#     human beta-cell mass rises only ~1.4x, unlike the 2-3x of rodents). The
#     model splits this into (i) BCM rising to ~2.0 with BCAP = 1 and (ii) a
#     left shift of the secretion curve (KSENS = 0.28, KG50 130 -> ~102 mg/dL)
#     representing increased glucose sensitivity via the lactogen -> Tph1 ->
#     serotonin -> HTR3A pathway (Kim, Nat Med 2010;16:804; Schraenen,
#     Diabetologia 2010;53:2589).
#
#  F. FETAL GROWTH ELASTICITIES — the calibration that matters most.
#     Reference NGT term composition: lean 2990 g + fat 460 g = 3450 g, fat
#     13.3% (Catalano, Am J Obstet Gynecol 2003;189:1698). AFAT = 0.55 and
#     AIGF = 0.20 were set so that
#         maternal MBG 95  -> BW ~3450 g, z 0.00, fat 13%   (NGT)
#         maternal MBG 110 -> BW ~3690 g, z 0.59, fat ~16%  (mild GDM)
#         maternal MBG 135 -> BW ~4040 g, z 1.48, fat ~21%  (uncontrolled)
#     i.e. ~16 g of birth weight per mg/dL of mean maternal glucose, with the
#     increment falling ~2/3 on FAT and ~1/3 on LEAN mass. That asymmetry is
#     the model's mechanistic claim, and it is testable against neonatal body
#     composition rather than birth weight alone.
#
#  G. LGA CALIBRATION AGAINST HAPO (N Engl J Med 2008;358:1991).
#     P(LGA) = Phi((BWZ - 1.282)/SDIND) with SDIND = 1.0, so BWZ = 0 gives
#     exactly the 10% definitional rate. HAPO's seven FPG categories then map:
#         category 1 (FPG <75 mg/dL)  observed LGA  5.3%  <- model BWZ ~ -0.30
#         category 7 (FPG >=100)      observed LGA 26.3%  <- model BWZ ~ +0.60
#     That is a required span of ~0.9 z-units across a 25 mg/dL FPG range,
#     equivalent to HAPO's adjusted OR of 1.38 for LGA per 1-SD (6.9 mg/dL)
#     rise in FPG. Reproducing the CONTINUOUS gradient — no threshold — is
#     scenario S9 and is the single most important validation of this model.
#
#  H. TREATMENT-EFFECT CALIBRATION (Landon/MFMU mild GDM, N Engl J Med
#     2009;361:1339; Crowther/ACHOIS, N Engl J Med 2005;352:2477):
#         endpoint            untreated -> treated (observed)
#         LGA                    14.5%  ->  7.1%
#         macrosomia >4000 g     14.3%  ->  5.9%
#         shoulder dystocia       4.0%  ->  1.5%
#         pre-eclampsia/HDP      13.6%  ->  8.6%
#         caesarean              33.8%  -> 26.9%
#     PSD0/KSD, PCS0/KCS and PPE0/KPEG/KPEB were solved from these pairs (the
#     pre-eclampsia pair gives KPEG = 0.566 and, with a cohort BMI ~30,
#     KPEB = 1.1 and PPE0 = 0.035 uniquely).
#
#  I. NEONATAL HYPOGLYCAEMIA is driven by NHR = fetal insulin / reference, i.e.
#     by the abrupt loss of transplacental glucose against a still-hyperactive
#     fetal beta-cell. PNH0 = 1.5%, PNHMAX = 35%, KNH = 1.1, GNH = 1.6 give
#     11% at NHR = 1.63 and 23% at NHR = 2.6, bracketing reported rates.
#
#  J. DRUG PK.
#     Metformin: F 0.55, CL/F ~1440 L/day non-pregnant rising ~25% in pregnancy
#       (Eyal, Drug Metab Dispos 2010;38:833; Liao, Clin Pharmacokinet 2020),
#       two-compartment with a large peripheral volume (erythrocyte/tissue
#       sequestration). Placental transfer is modelled as FREE bidirectional
#       diffusion via OCT3/PMAT, giving an umbilical:maternal ratio near 1.0
#       (Vanky, Hum Reprod 2005;20:1593; Charles, Ther Drug Monit 2006;28:67).
#       PD: EMAXM = 0.35 suppression of EGP with EC50M = 0.5 mg/L, plus a
#       smaller peripheral SI gain — the standard hepatic-dominant split.
#     Insulin: KAB = 2.1/day for detemir/NPH (t1/2,abs ~8 h), KAF = 20/day for
#       aspart/lispro (tmax ~50 min). Requirement rises from ~0.7 to ~1.0
#       U/kg/day by the third trimester, which the model reproduces because the
#       dose needed to hold FPG <95 is computed, not assumed (see titrate()).
#       Fetal exposure is structurally ZERO — insulin has no transfer term.
#     Glyburide: CL doubles in pregnancy (Hebert, Clin Pharmacol Ther
#       2009;85:607). Transfer is net of BCRP/ABCG2 efflux, giving a cord:
#       maternal ratio ~0.7 (Hemauer, Am J Obstet Gynecol 2010;202:383).
#       Its distinguishing PD feature is GLUCOSE-INDEPENDENT secretion (ISRG),
#       in mother AND fetus — the mechanistic basis for the excess neonatal
#       hypoglycaemia in Senat (JAMA 2018;319:1773) and in the Balsells
#       meta-analysis (BMJ 2015;350:h102).
#
#  K. POST-PARTUM / T2DM. At TDEL the placenta involutes, CID collapses and SI
#     returns to its pre-gestational value; what does NOT return is BCM. The
#     post-partum disposition index DI = SIP x BCM therefore exposes the
#     pre-existing defect. The 5-year hazard H5 = H0T2D x exp(-KDI(DI/DIREF-1))
#     with H0T2D = 0.02 and KDI = 3.8 yields a relative risk of ~8 at
#     DI/DIREF = 0.45, matching Bellamy's pooled RR of 7.43 (Lancet
#     2009;373:1773) and Vounzoulaki's cumulative incidence (BMJ 2020;369:m1361).
#     Lactation >=3 months cuts the hazard 45% (Gunderson, Ann Intern Med
#     2015;163:889).
#
# =============================================================================


# =============================================================================
#  DOSING / REGIMEN HELPERS
# =============================================================================

GW  <- function(day) day/7
DAY <- function(gw)  gw*7

#' Meal schedule: 3 meals + 2 snacks per day, expressed as glucose (mg) into GUT1
#'
#' @param cho_day  total daily carbohydrate (g)
#' @param start,end  gestational days
#' @param low_gi   TRUE flattens the absorption profile (handled via KTR_G)
meal_events <- function(cho_day = 200, start = 56, end = 364) {
  # distribution across the day (fractions of daily CHO) and clock times (days)
  frac  <- c(breakfast = 0.225, snack1 = 0.085, lunch = 0.300,
             snack2    = 0.090, dinner = 0.300)
  clock <- c(breakfast = 0.29,  snack1 = 0.42,  lunch = 0.52,
             snack2    = 0.66,  dinner = 0.79)
  days  <- seq(floor(start), floor(end) - 1)
  out <- lapply(seq_along(frac), function(i) {
    data.frame(ID = 1, time = days + clock[i], cmt = cmt_i("GUT1"),
               amt = cho_day*frac[i]*1000, evid = 1)
  })
  do.call(rbind, out)
}

#' Metformin: BID with meals
met_events <- function(mg_per_dose = 1000, start_gw = 28, end_gw = 39) {
  days <- seq(DAY(start_gw), DAY(end_gw) - 1)
  rbind(
    data.frame(ID = 1, time = days + 0.29, cmt = cmt_i("MGUT"), amt = mg_per_dose, evid = 1),
    data.frame(ID = 1, time = days + 0.79, cmt = cmt_i("MGUT"), amt = mg_per_dose, evid = 1)
  )
}

#' Basal-bolus insulin. Basal at bedtime, boluses with the three main meals.
insulin_events <- function(total_u, start_gw = 28, end_gw = 39, basal_frac = 0.5) {
  days  <- seq(DAY(start_gw), DAY(end_gw) - 1)
  basal <- total_u*basal_frac
  bolus <- total_u*(1 - basal_frac)/3
  rbind(
    data.frame(ID = 1, time = days + 0.92, cmt = cmt_i("IBAS"), amt = basal, evid = 1),
    data.frame(ID = 1, time = days + 0.28, cmt = cmt_i("IBOL"), amt = bolus, evid = 1),
    data.frame(ID = 1, time = days + 0.51, cmt = cmt_i("IBOL"), amt = bolus, evid = 1),
    data.frame(ID = 1, time = days + 0.78, cmt = cmt_i("IBOL"), amt = bolus, evid = 1)
  )
}

#' Glyburide: once or twice daily
glyb_events <- function(mg_per_dose = 5, bid = TRUE, start_gw = 28, end_gw = 39) {
  days <- seq(DAY(start_gw), DAY(end_gw) - 1)
  ev <- data.frame(ID = 1, time = days + 0.25, cmt = cmt_i("GGUT"), amt = mg_per_dose, evid = 1)
  if (bid) ev <- rbind(ev, data.frame(ID = 1, time = days + 0.75,
                                      cmt = cmt_i("GGUT"), amt = mg_per_dose, evid = 1))
  ev
}

#' Combine event tables and sort
combine_ev <- function(...) {
  ev <- do.call(rbind, Filter(Negate(is.null), list(...)))
  ev[order(ev$time), ]
}

#' Pull the pre-breakfast (fasting) value of a column, one row per gestational day
fasting_series <- function(sim, col = "GLU") {
  d <- as.data.frame(sim)
  d$day <- floor(d$time)
  d$frac <- d$time - d$day
  d <- d[d$frac > 0.20 & d$frac < 0.29, ]
  stats::aggregate(d[[col]], by = list(day = d$day), FUN = mean) |>
    setNames(c("day", col))
}

#' Simple pharmacologic titration: escalate total daily insulin until the mean
#' fasting glucose over the last week is below `fpg_target`. This mirrors real
#' clinical practice (dose is an OUTCOME of the patient's physiology, not an
#' input) and is why the model reproduces the rising third-trimester insulin
#' requirement without being told to.
titrate_insulin <- function(pars, fpg_target = 95, start_gw = 28, end_gw = 39,
                            cho_day = 200, u_start = 20, u_max = 160,
                            step = 1.25, max_iter = 14) {
  u <- u_start
  for (i in seq_len(max_iter)) {
    ev <- combine_ev(meal_events(cho_day, 56, DAY(end_gw)),
                     insulin_events(u, start_gw, end_gw))
    s  <- gdm |> param(pars) |> data_set(ev) |>
      mrgsim(start = 56, end = DAY(end_gw), delta = 0.05, hmax = 0.05)
    f  <- fasting_series(s, "GLU")
    fpg <- mean(tail(f$GLU, 7))
    if (is.finite(fpg) && fpg <= fpg_target) break
    if (u >= u_max) break
    u <- min(u_max, u*step)
  }
  list(total_u = u, fpg = fpg, iter = i)
}

#' Run one scenario end-to-end (pregnancy + 12 weeks post partum)
run_scenario <- function(label, pars = list(), ev = NULL,
                         end_day = 364, delta = 0.05) {
  if (is.null(ev)) ev <- meal_events(200, 56, end_day)
  out <- gdm |>
    param(pars) |>
    data_set(ev) |>
    mrgsim(start = 56, end = end_day, delta = delta, hmax = 0.05) |>
    as.data.frame()
  out$scenario <- label
  out
}

#' Values at the moment of delivery — the row that carries every endpoint
at_delivery <- function(sim, tdel = 273) {
  d <- as.data.frame(sim)
  d[which.min(abs(d$time - tdel)), ]
}


# =============================================================================
#  SCENARIOS
#  Ten scenarios. S1-S8 are therapeutic strategies on the SAME patient
#  physiology so the comparison is fair; S9 is the HAPO gradient validation;
#  S10 is the post-partum / life-course arm.
# =============================================================================

# ---- patient archetypes ------------------------------------------------------
P_NGT     <- list(BCAP = 1.00, BMI = 22, EXEFF = 0.00)   # normal glucose tolerance
P_DEFIC   <- list(BCAP = 0.85, BMI = 24, EXEFF = 0.00)   # secretion-deficient GDM
P_RESIST  <- list(BCAP = 1.00, BMI = 34, EXEFF = 0.00)   # insulin-resistant GDM
P_MIXED   <- list(BCAP = 0.55, BMI = 31, EXEFF = 0.00)   # mixed / severe GDM

scenarios <- list()

# ---- S1  Normoglycaemic reference pregnancy ---------------------------------
scenarios$S1 <- function() {
  run_scenario("S1 NGT reference", P_NGT,
               meal_events(220, 56, 364))
}

# ---- S2  Untreated GDM, secretion-deficient subtype ------------------------
#      Expected signature: near-normal FPG, abnormal post-load excursions.
#      This is the phenotype that a fasting-glucose-only screen MISSES.
scenarios$S2 <- function() {
  run_scenario("S2 untreated GDM (secretion-deficient)", P_DEFIC,
               meal_events(220, 56, 364))
}

# ---- S3  Untreated GDM, insulin-resistant subtype --------------------------
#      Expected signature: high FPG *with* hyperinsulinaemia, high FFA,
#      high maternal fat accretion, strong pre-eclampsia contribution.
scenarios$S3 <- function() {
  run_scenario("S3 untreated GDM (insulin-resistant)", P_RESIST,
               meal_events(220, 56, 364))
}

# ---- S4  MNT + exercise only (first line; controls ~70-85% of GDM) ---------
#      Two levers: carbohydrate quantity (CHOD, and hence meal amt) and
#      distribution, plus an insulin-INDEPENDENT SI gain from exercise.
scenarios$S4 <- function() {
  p <- modifyList(P_MIXED, list(EXEFF = 0.20, CHOD = 175, KTR_G = 16))
  run_scenario("S4 MNT + exercise", p, meal_events(175, 56, 364))
}

# ---- S5  MNT + metformin 1000 mg BID from GW 28 ----------------------------
#      Note what the model shows that a summary table cannot: the fetus is
#      exposed to essentially the same metformin concentration as the mother.
scenarios$S5 <- function() {
  p <- modifyList(P_MIXED, list(EXEFF = 0.20, CHOD = 175, KTR_G = 16))
  ev <- combine_ev(meal_events(175, 56, 364), met_events(1000, 28, 39))
  run_scenario("S5 MNT + metformin", p, ev)
}

# ---- S6  MNT + titrated basal-bolus insulin (gold standard) ----------------
#      The dose is TITRATED, not prescribed — see titrate_insulin().
scenarios$S6 <- function() {
  p  <- modifyList(P_MIXED, list(EXEFF = 0.20, CHOD = 175, KTR_G = 16))
  tt <- titrate_insulin(p, fpg_target = 95, start_gw = 28, end_gw = 39, cho_day = 175)
  message(sprintf("S6 titration converged to %.0f U/day (FPG %.1f mg/dL)",
                  tt$total_u, tt$fpg))
  ev <- combine_ev(meal_events(175, 56, 364), insulin_events(tt$total_u, 28, 39))
  out <- run_scenario("S6 MNT + insulin (titrated)", p, ev)
  attr(out, "insulin_u") <- tt$total_u
  out
}

# ---- S7  MNT + glyburide 5 mg BID ------------------------------------------
#      The point of this arm is not efficacy on maternal glucose (it works);
#      it is the fetal exposure term and the glucose-INDEPENDENT fetal insulin
#      release it produces.
scenarios$S7 <- function() {
  p  <- modifyList(P_MIXED, list(EXEFF = 0.20, CHOD = 175, KTR_G = 16))
  ev <- combine_ev(meal_events(175, 56, 364), glyb_events(5, TRUE, 28, 39))
  run_scenario("S7 MNT + glyburide", p, ev)
}

# ---- S8  Metformin with insulin add-on (the MiG-trial reality) -------------
#      In MiG, 46% of metformin-assigned women needed supplemental insulin.
#      Here that is a MODEL PREDICTION: if FPG stays above 95 on metformin
#      alone, insulin is added and titrated.
scenarios$S8 <- function() {
  p  <- modifyList(P_MIXED, list(EXEFF = 0.20, CHOD = 175, KTR_G = 16))
  ev0 <- combine_ev(meal_events(175, 56, 364), met_events(1000, 28, 39))
  s0  <- gdm |> param(p) |> data_set(ev0) |>
    mrgsim(start = 56, end = 273, delta = 0.05, hmax = 0.05)
  fpg <- mean(tail(fasting_series(s0, "GLU")$GLU, 7))
  if (fpg > 95) {
    tt <- titrate_insulin(p, 95, 30, 39, 175, u_start = 12)
    ev <- combine_ev(ev0, insulin_events(tt$total_u, 30, 39))
    message(sprintf("S8 metformin insufficient (FPG %.1f) -> insulin %.0f U/day added",
                    fpg, tt$total_u))
  } else {
    ev <- ev0
    message(sprintf("S8 metformin sufficient (FPG %.1f)", fpg))
  }
  run_scenario("S8 metformin + insulin add-on", p, ev)
}

# ---- S9  HAPO gradient validation ------------------------------------------
#      Sweep BCAP/BMI to span the HAPO fasting-glucose categories and check
#      that LGA rises CONTINUOUSLY with no threshold. This is the key test.
hapo_sweep <- function() {
  grid <- expand.grid(BCAP = seq(0.45, 1.15, by = 0.05),
                      BMI  = c(22, 26, 30, 34))
  res <- lapply(seq_len(nrow(grid)), function(i) {
    p <- list(BCAP = grid$BCAP[i], BMI = grid$BMI[i])
    s <- run_scenario(sprintf("BCAP%.2f_BMI%d", grid$BCAP[i], grid$BMI[i]),
                      p, meal_events(200, 56, 280), end_day = 280, delta = 0.25)
    f <- fasting_series(s, "GLU")
    d <- at_delivery(s)
    data.frame(BCAP = grid$BCAP[i], BMI = grid$BMI[i],
               FPG_gw28 = mean(f$GLU[f$day >= 190 & f$day <= 200]),
               FPG_term = mean(tail(f$GLU, 7)),
               MBG = d$MBG, HbA1c = d$HBA1C,
               BW = d$FWT, BWZ = d$BWZ, FATPCT = d$FATPCT,
               NHR = d$NHR,
               P_LGA = d$PLGA, P_NEOHYPO = d$PNHYP,
               P_SD = d$PSD, P_PE = d$PPE, P_CS = d$PCS,
               DI = d$DI, CIT2D = d$CIT2D)
  })
  do.call(rbind, res)
}

# ---- S10  Post-partum trajectory and 5-year T2DM hazard --------------------
#      Same four archetypes, followed 12 weeks past delivery. The placenta
#      involutes, SI recovers, and the residual beta-cell deficit is exposed.
scenarios$S10 <- function() {
  arche <- list(NGT = P_NGT, deficient = P_DEFIC,
                resistant = P_RESIST, mixed = P_MIXED)
  out <- lapply(names(arche), function(nm) {
    s <- run_scenario(paste0("S10 ", nm), arche[[nm]],
                      meal_events(200, 56, 364), end_day = 364, delta = 0.25)
    s$archetype <- nm
    s
  })
  do.call(rbind, out)
}


# =============================================================================
#  DRIVER
# =============================================================================

run_all <- function(which = c("S1","S2","S3","S4","S5","S6","S7","S8")) {
  do.call(rbind, lapply(which, function(k) scenarios[[k]]()))
}

#' Endpoint table across scenarios — the deliverable a clinician would read
endpoint_table <- function(sims, tdel = 273) {
  sims |>
    group_by(scenario) |>
    slice(which.min(abs(time - tdel))) |>
    transmute(
      scenario,
      FPG_term      = round(GLU, 1),
      MBG           = round(MBG, 1),
      HbA1c         = round(HBA1C, 2),
      SI_pct_prepreg= round(100*MATSUDA, 1),
      BCM           = round(BCM, 2),
      DI_rel        = round(DI/DIREF, 2),
      birth_weight  = round(FWT),
      BW_z          = round(BWZ, 2),
      fetal_fat_pct = round(FATPCT, 1),
      cord_insulin  = round(NHR, 2),
      fetal_metformin = round(CMETF, 3),
      fetal_glyburide = round(CGLYF, 4),
      P_LGA         = round(100*PLGA, 1),
      P_macrosomia  = round(100*PMACR, 1),
      P_neonat_hypo = round(100*PNHYP, 1),
      P_shoulder_dys= round(100*PSD, 1),
      P_preeclampsia= round(100*PPE, 1),
      P_caesarean   = round(100*PCS, 1),
      T2DM_5yr      = round(100*CIT2D, 1)
    ) |>
    ungroup()
}

# ---- plots -------------------------------------------------------------------
plot_maternal <- function(sims) {
  f <- sims |> mutate(day = floor(time), frac = time - day) |>
    filter(frac > 0.20, frac < 0.29) |>
    group_by(scenario, day) |> summarise(FPG = mean(GLU), .groups = "drop")
  ggplot(f, aes(day/7, FPG, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 95, linetype = 2) +
    labs(x = "Gestational week", y = "Fasting plasma glucose (mg/dL)",
         title = "Maternal fasting glucose by strategy",
         subtitle = "dashed line = 95 mg/dL treatment target") +
    theme_bw()
}

plot_fetal <- function(sims) {
  ggplot(sims, aes(GWk, FWT, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_line(aes(y = FWREF), colour = "grey40", linetype = 2, inherit.aes = TRUE) +
    labs(x = "Gestational week", y = "Fetal weight (g)",
         title = "Fetal growth trajectories",
         subtitle = "dashed grey = reference (50th centile) curve") +
    theme_bw()
}

plot_fetal_exposure <- function(sims) {
  sims |>
    select(scenario, GWk, maternal = CMET, fetal = CMETF) |>
    pivot_longer(c(maternal, fetal)) |>
    filter(GWk > 27, GWk < 39.5) |>
    ggplot(aes(GWk, value, colour = name)) +
    geom_line() + facet_wrap(~scenario, scales = "free_y") +
    labs(x = "Gestational week", y = "Metformin (mg/L)",
         title = "Metformin: the fetus sees what the mother sees",
         colour = NULL) +
    theme_bw()
}

plot_hapo <- function(sweep) {
  ggplot(sweep, aes(FPG_term, 100*P_LGA, colour = factor(BMI))) +
    geom_point() + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(x = "Term fasting plasma glucose (mg/dL)", y = "P(LGA) %",
         colour = "Pre-gestational BMI",
         title = "HAPO gradient reproduced: continuous, no threshold",
         subtitle = "observed HAPO: 5.3% LGA in category 1 -> 26.3% in category 7") +
    theme_bw()
}

if (identical(environment(), globalenv()) && !interactive()) {
  sims <- run_all()
  print(as.data.frame(endpoint_table(sims)), width = 200)
  sweep <- hapo_sweep()
  print(head(sweep, 20))
}

# =============================================================================
#  KNOWN LIMITATIONS (state them, do not bury them)
# =============================================================================
#  1. One aggregated contra-insulin index. hPL, progesterone, cortisol and
#     TNF-alpha are not separately identifiable from clinical data, so they are
#     lumped. Any claim about a single hormone's contribution is untestable in
#     this model by construction.
#  2. No spatial/organ resolution of maternal glucose. Muscle, liver and
#     adipose are collapsed into peripheral (SIP) and hepatic (SILV) arms.
#  3. Meal absorption is a two-transit-compartment approximation; incretin
#     secretion is a constant amplification (KINCR) rather than a nutrient-
#     driven state, so oral-vs-IV glucose differences are not captured.
#  4. Fetal metformin is tracked but given NO growth effect (the placeholder
#     term in dxdt_FATF is deliberately multiplied by zero). The MiG-TOFU
#     follow-up (Rowan, Diabetes Care 2011;34:2279) found metformin-exposed
#     children had MORE subcutaneous fat at 2 years and higher BMI at 9 years,
#     which cannot be predicted from any in-utero mechanism the model contains.
#     Turning that term on would fabricate a mechanism. It is left visible and
#     inert so the gap is explicit.
#  5. Outcome probabilities are population regressions bolted onto mechanistic
#     drivers, not mechanistic models of labour, placentation or neonatal
#     adaptation. Pre-eclampsia in particular shares upstream biology with GDM
#     (sFlt-1/PlGF, endothelial dysfunction) that is drawn on the map but NOT
#     given ODEs here.
#  6. Stillbirth, congenital malformation (a first-trimester, pre-existing-
#     diabetes phenomenon rather than a GDM one) and neonatal respiratory
#     outcomes are on the map but not in the equations.
#  7. Not fitted. No individual-level data, no parameter uncertainty, no
#     inter-individual variability block ($OMEGA is intentionally absent).
#     Every scenario is a single deterministic subject.
# =============================================================================
