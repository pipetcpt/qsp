## =====================================================================
##  PRIMARY POSTPARTUM HAEMORRHAGE — QSP model (mrgsolve)
##  41 ODE compartments · 10 drug PK · 31 disease/physiology · 12 scenario blocks
##  70 kg term parturient, immediately after delivery of the placenta.
##  Time unit: MINUTES.  Volumes mL.  Masses g (haemoglobin, fibrinogen).
## =====================================================================
##
##  THE THESIS THIS MODEL IS BUILT AROUND
##  -------------------------------------
##  Postpartum haemorrhage is usually taught as a coagulation problem with a
##  mechanical preamble.  This model asserts the reverse, and makes the
##  reversal arithmetic.
##
##  (1) PPH IS A FLUX PROBLEM AGAINST A STORE THAT ONLY LOOKS LARGE.
##      Term uterine blood flow is UBF0 = 750 mL/min — 10-15% of a 7 L/min
##      cardiac output — delivered to ~120 spiral-artery stumps that
##      trophoblast remodelling has stripped of their own smooth muscle, so
##      they CANNOT constrict.  Pregnancy's celebrated volume expansion is
##      2100 mL (70 -> 100 mL/kg):
##
##          2100 mL / 750 mL/min  =  2.8 MINUTES of unopposed flow.
##
##      The model's totally atonic uterus (ATN = 1) reaches 516 mL/min and:
##
##          500 mL at  2.0 min        2100 mL (ATLS class III) at 6.0 min
##         1000 mL at  4.0 min        exsanguination           at  10  min
##
##      Every therapeutic claim below is a claim about how much of that clock
##      an intervention buys.
##
##  (2) THE VALVE IS POISEUILLE, NOT LINEAR.
##      Haemostasis at the placental bed is mechanical: interlacing
##      myometrial bundles ("the living ligature") compress each vessel.
##      Flow through a compressed vessel goes as the fourth power of radius,
##      so the model writes patency = (1 - TONE)^4:
##
##          TONE 0.91 -> 6e-5   (0.05 mL/min)   a normal, firm uterus
##          TONE 0.50 -> 0.063  (  47 mL/min)
##          TONE 0.20 -> 0.410  ( 307 mL/min)   a "boggy" uterus
##
##      The clinical distinction between firm and boggy is a few percent of
##      contractile capacity, and the model is correspondingly BISTABLE:
##
##          ATN 0.60 -> 2218 mL, survives, MAP nadir 78 mmHg
##          ATN 0.63 -> DEAD at 70 min (3544 mL)
##
##  (3) COAGULATION IS LOCKED OUT UNTIL THE MECHANICS WORK.
##      Clot formation is shear-gated, fSH = 1/(1 + (Q/80)^2), and a clot
##      over an uncompressed sinus can never take more than 35% of the
##      patency (CSEAL0) however good the coagulation profile is.  At 10 min:
##
##          normal            Q   0 mL/min  fSH 1.00  CLT -> 0.98, MAT -> 1.00
##          moderate atony    Q  39         fSH 0.81  CLT    0.83
##          severe UNTREATED  Q 168         fSH 0.19  CLT plateaus at 0.22
##          severe TREATED    Q  22         fSH 0.93  CLT    0.87
##
##      So "uterotonics and tamponade before blood products" is not a
##      tradition, it is the precondition for the products to do anything.
##      Durable haemostasis additionally requires the clot to survive long
##      enough at low shear to ORGANISE (state MAT) — which is what every
##      mechanical manoeuvre is really buying.
##
##  (4) OXYTOCIN DELETES ITS OWN TARGET.
##      OTR internalises on agonist exposure (KDES 0.08/min, occupancy-
##      driven) and recovers with t1/2 2.6 h.  Consequences, both derived:
##
##      · 8 h of augmentation at 12 mU/min -> steady-state OTR 0.39, matching
##        the 30-60% loss of myometrial binding sites reported after
##        oxytocin-augmented labour;
##      · treatment itself desensitises — from OTR 1.00, a standard infusion
##        reaches OTR 0.09 by 120 min.
##
##      In an augmented patient with severe atony (2605 mL on oxytocin +
##      massage), the two available moves are NOT equivalent:
##
##          DOUBLE the oxytocin      2462 mL   (-143)
##          QUADRUPLE the oxytocin   2328 mL   (-277)
##          add carboprost (FP)      1680 mL   (-925)
##          add ergometrine (a1)     1781 mL   (-824)
##          add misoprostol (EP2/3)  1953 mL   (-652)
##
##      Crossing to a different receptor is worth 3-6x escalating the dose of
##      the drug whose receptor has already gone.
##
##  (5) FIBRINOGEN CROSSES ITS CRITICAL LINE FIRST — BY RATIO ARITHMETIC.
##      What decides the order in which a diluted, bleeding patient becomes
##      deficient is critical value / starting value, nothing else:
##
##          fibrinogen 2.0 / 4.5  = 0.44   <- largest ratio, crosses first
##          pooled factors 0.30 / 1.00 = 0.30
##          platelets 50 / 250 = 0.20
##
##      In the crystalloid-first arm the model confirms the ordering:
##      fibrinogen crosses 2.0 g/L at 3476 mL of cumulative loss, Hb crosses
##      7 g/dL at 3555 mL, and platelets (nadir 71) and pooled factors
##      (nadir 0.33) never cross at all.
##
##      Corollary, also arithmetic: FFP CANNOT RAISE A LOW FIBRINOGEN.  FFP
##      contains 2 g/L; raising 7 L of plasma by 1 g/L needs 7 g = 14 units =
##      3.5 L of FFP.  FFP's concentration IS the transfusion target, so it
##      drives fibrinogen toward 2 g/L from either side.
##
##  (6) THE RESUSCITATION THAT SAVES THE CIRCULATION DEGRADES THE LIGATURE.
##      Calcium is a coagulation cofactor AND the contractile ion of smooth
##      muscle, so citrate-driven hypocalcaemia is a double hit unique to
##      obstetric haemorrhage; cold products cool both enzymes and muscle;
##      and crystalloid raises MAP, which — because the uterine circulation
##      is pressure-passive — raises the leak.  Isolated in the model:
##
##          goal-directed 1:1 + fibrinogen + calcium + warmer   3489 mL
##          crystalloid-first, RBC only                         3953 mL
##            (fibrinogen 1.37 g/L, 34.6 C, pH 7.31, 6.2 units)
##          MTP without calcium replacement            +922 mL, iCa 0.45
##          MTP without a fluid warmer     +2873 mL, 31.8 C, 16 vs 11 units
##
##  (7) TXA'S WINDOW IS NOT THREE HOURS. IT IS THE DURATION OF BLEEDING.
##      TXA makes no clot; it protects the clot the mechanics allowed to
##      form.  In a refractory patient whose bleeding lasts ~40 min:
##
##          no TXA 4099 mL · 5 min 3400 (-699) · 15 min 3645 (-454) ·
##          30 min 4030 (-69) · 45/90/180 min 4099 (0)
##
##      The WOMAN trial's 3-hour boundary is the population's bleeding
##      duration, not a property of the molecule.
##
##  (8) TIMING BEATS DEVICE.  Catastrophic atony (ATN 0.94) on the full
##      pharmacological bundle dies at 80 min having lost 10 L.  Adding:
##
##          balloon at 20 min    4528 mL, survives, O2 debt   152 mL
##          balloon at 45 min    7617 mL, survives, O2 debt  1195 mL
##          balloon at 75 min    DEAD at 150 min, O2 debt    8260 mL
##          aortic compression 8-35 min then balloon at 35:
##                               2782 mL, NO transfusion, O2 debt 0
##          uterine artery ligation at 50 min  9624 mL, 10 units, iCa 0.55
##          hysterectomy at 60 min             8451 mL, survives
##          hysterectomy at 120 min            a post-mortem procedure
##
##  (9) HAEMOGLOBIN IS THE LAST THING TO MOVE.  Whole blood leaves at
##      constant concentration, so the untreated arms reach death with Hb
##      11.2-11.4 g/dL.  A haemoglobin-triggered transfusion policy is
##      structurally too late; the model's MTP triggers on ONGOING LOSS
##      (>1500 mL with Q > 25 mL/min), which is what protocols actually do.
##
##  WHAT IS FITTED AND WHAT IS DERIVED
##  ----------------------------------
##  FITTED (6 things):
##    1. UBF0 and the tone->patency exponent NP, set so that a normal third
##       stage loses ~270 mL (measured average 300-500 mL) and closes in
##       ~5 min.
##    2. The ATN -> (DRVI, CAP) map, set so that the treatment ladder
##       reproduces clinical experience: single agents suffice in moderate
##       atony, oxytocin alone fails in severe atony, and a fraction of cases
##       are refractory to the whole pharmacological ladder.
##    3. KFORM, QS50, CSEAL0, KMAT, set so that a low-flow bed thromboses
##       durably over ~15-30 min while a high-flow bed does not thrombose at
##       all.
##    4. KDES/KREC, set so that 8 h of augmentation leaves OTR ~0.39
##       (Phaneuf/Robinson-type binding-site data, 30-60% loss).
##    5. KFDEG/KLYS/KTPA*, set so that severe PPH reaches fibrinogen <2 g/L
##       and 3-4x baseline plasmin, and so that TXA's benefit is confined to
##       the bleeding period.
##    6. Citrate/calcium and thermal constants, set so that an obstetric MTP
##       produces ionised calcium 0.5-0.9 mmol/L and ~2.8 C/h of cooling per
##       100 mL/min of unwarmed 4 C product.
##
##  DERIVED, i.e. left as falsifiable output:
##    - the 2.8-minute reserve and the 6-minute class-III crossing
##    - the ATN 0.60 / 0.63 knife edge
##    - the whole who-needs-what escalation ladder (scenario block 3)
##    - dose-escalation futility vs receptor-crossing benefit
##    - fibrinogen-crosses-first, and the OBS2 empiric-fibrinogen null result
##    - TXA's window as the duration of bleeding
##    - timing dominating device for mechanical control
##    - haemoglobin as the last variable to move
##
##  PROVENANCE OF THE NUMBERS QUOTED ABOVE
##  --------------------------------------
##  R and mrgsolve were NOT available in the environment where this file was
##  written, so this .R file has not been executed here.  Every number quoted
##  in these comments and in README.md was produced by an independent
##  reference implementation of the IDENTICAL equations and parameters
##  (pure-python RK4, dt = 0.01 min).  Discrete doses are given as
##  instantaneous boluses there; the event tables below reproduce that, with
##  infusions where the label specifies one.  If you re-run this file and a
##  number disagrees, trust the R output and please report the difference.
##
##  Calibration literature is listed in pph_references.md; the clinical
##  anchors used most heavily are the WOMAN trial and its timing IPD
##  analysis, OBS2, FIB-PPH, CHAMPION, Charbit/Cortet on fibrinogen,
##  Hiippala on dilutional order, Suarez on balloon tamponade, Gillissen and
##  Henriquez (TeMpOH) on fluid management, and Pacheco on obstetric MTP.
##
##  NOT FOR CLINICAL USE.  Educational / research QSP model.
## =====================================================================

library(mrgsolve)
library(dplyr)

pph_code <- '
$PROB
# Primary postpartum haemorrhage — mechanical haemostasis, shear-gated
# coagulation, receptor desensitisation and the cost of resuscitation.

$PARAM @annotated
// ---------------- subject ------------------------------------------------
WT      :  70.0  : body weight (kg)
V0      : 7000.0 : blood volume at term (mL, 100 mL/kg)
CO0     :   7.0  : cardiac output at term (L/min)
MAP0    :  85.0  : baseline mean arterial pressure (mmHg)
HB0     :  11.5  : baseline haemoglobin (g/dL)
FIB0    :   4.5  : baseline fibrinogen (g/L)
PLT0    : 250.0  : baseline platelet count (10^9/L)
FCT0    :   1.0  : baseline pooled clotting-factor activity (fraction)
LOSS0   : 250.0  : blood lost during placental separation (mL)

// ---------------- uterus: severity and mechanics -------------------------
UBF0    : 750.0  : term uterine blood flow (mL/min)
ATN     :   0.0  : atony severity 0-1 (the single disease knob)
DRVMAX  :   4.5  : intrinsic postpartum contractile drive when ATN = 0
KCAP    :   0.45 : loss of contractile CEILING per unit ATN
KDRV    :   1.4  : drive giving half-maximal tone
HDRV    :   2.0  : Hill coefficient on drive
KON     :   0.5  : rate constant of tone change (/min)
TONE0   :   0.55 : tone at the moment of placental separation
NP      :   4.0  : tone -> patency exponent (Poiseuille)
RET0    :   0.0  : retained tissue fraction (caps maximum tone)
KTR     :   0.0  : trauma (laceration) leak at MAP0 (mL/min)

// ---------------- oxytocin -----------------------------------------------
V_OXY   :  12.0  : oxytocin distribution volume (L)
CL_OXY  :   2.08 : oxytocin clearance (L/min), t1/2 ~4 min
EMAX_OXY:   2.6  : maximum oxytocin contractile drive
EC50_OXY:   0.06 : oxytocin EC50 at OTR = 1 (IU/L)

// ---------------- carbetocin ---------------------------------------------
V_CBT   :  22.0  : carbetocin distribution volume (L)
KE_CBT  :   0.01733 : carbetocin elimination (/min), t1/2 40 min
EMAX_CBT:   2.4  : maximum carbetocin contractile drive
EC50_CBT:   0.9  : carbetocin EC50 at OTR = 1 (ug/L)

// ---------------- ergometrine (alpha1 / 5-HT2A) --------------------------
KA_ERG  :   0.30 : IM absorption rate constant (/min)
F_ERG   :   0.85 : IM bioavailability
V_ERG   :  60.0  : distribution volume (L)
KE_ERG  :   0.0077 : elimination (/min), t1/2 90 min
EMAX_ERG:   1.6  : maximum ergometrine contractile drive
EC50_ERG:   2.0  : ergometrine EC50 (ug/L)

// ---------------- carboprost (15-methyl-PGF2alpha, FP receptor) ----------
KA_PGF  :   0.08 : IM absorption rate constant (/min)
F_PGF   :   0.90 : IM bioavailability
V_PGF   :  20.0  : distribution volume (L)
KE_PGF  :   0.0231 : elimination (/min), t1/2 30 min
EMAX_PGF:   2.2  : maximum carboprost contractile drive
EC50_PGF:   1.5  : carboprost EC50 (ug/L)

// ---------------- misoprostol (PGE1, EP2/EP3) ---------------------------
KA_MSO  :   0.05 : sublingual absorption rate constant (/min)
F_MSO   :   0.70 : sublingual bioavailability
V_MSO   :  60.0  : distribution volume (L)
KE_MSO  :   0.023 : elimination of misoprostol acid (/min), t1/2 30 min
EMAX_MSO:   1.3  : maximum misoprostol contractile drive
EC50_MSO:   2.5  : misoprostol EC50 (ug/L)

// ---------------- oxytocin receptor dynamics -----------------------------
OTR0    :   1.0  : available OTR fraction at delivery (1 = agonist-naive)
KREC    :   0.0045 : receptor recovery (/min), t1/2 2.6 h
KDES    :   0.08 : occupancy-driven desensitisation (/min)

// ---------------- clot formation and the shear gate ----------------------
KFORM   :   0.18 : maximal placental-bed clot formation rate (/min)
KFIB    :   1.8  : fibrinogen giving half-maximal clot rate (g/L), Hill 3
KPLT    :  60.0  : platelet count giving half-maximal clot rate (10^9/L)
KFCT    :   0.35 : factor activity giving half-maximal clot rate
QS50    :  80.0  : leak flow at which the shear gate is half closed (mL/min)
CSEAL0  :   0.35 : fraction of patency a clot can seal with NO tone
KMAT    :   0.05 : thrombus organisation rate at zero shear (/min)
KLYS    :   0.06 : clot lysis per unit excess plasmin (/min)

// ---------------- fibrinolysis -------------------------------------------
KTPAB   :   0.14 : basal tPA appearance (units/min)
KTPACL  :   0.14 : tPA clearance (/min), t1/2 5 min
KTPAINJ :   1.4  : flow-driven tPA release
KTPASHK :   1.6  : shock-driven tPA release
KACT    :   0.35 : plasminogen activation rate
KPLSCL  :   0.35 : plasmin inactivation (/min), t1/2 2 min
KFDEG   :   0.0045 : fibrinogenolysis per unit excess plasmin (/min)

// ---------------- tranexamic acid ---------------------------------------
V1_TX   :  10.0  : TXA central volume (L)
V2_TX   :   8.0  : TXA peripheral volume (L)
K12_TX  :   0.12 : TXA central -> peripheral (/min)
K21_TX  :   0.15 : TXA peripheral -> central (/min)
KE_TX   :   0.011 : TXA elimination (/min)
ITX50   :   8.0  : TXA concentration for half-maximal effect (mg/L)
ETXACT  :   0.85 : maximum inhibition of plasminogen activation
ETXPROT :   0.95 : maximum protection of formed clot

// ---------------- blood products (per unit) ------------------------------
RBC_VOL : 280.0  : red-cell unit volume (mL)
RBC_HB  :  55.0  : haemoglobin per red-cell unit (g)
RBC_CIT :   3.0  : citrate per red-cell unit (mmol)
FFP_VOL : 250.0  : FFP unit volume (mL)
FFP_FIB :   0.5  : fibrinogen per FFP unit (g) — i.e. 2 g/L
FFP_FCT :   0.8  : factor activity of FFP
FFP_CIT :  12.0  : citrate per FFP unit (mmol)
PLT_VOL : 250.0  : platelet unit volume (mL)
PLT_CNT : 300.0  : platelets per unit (10^9)
PLT_FIB :   0.2  : fibrinogen per platelet unit (g)
PLT_CIT :   8.0  : citrate per platelet unit (mmol)
FGC_VOL : 100.0  : volume per 2 g of fibrinogen concentrate (mL)
FGC_FIB :   2.0  : fibrinogen per 100 mL of concentrate (g)

// ---------------- citrate / calcium -------------------------------------
CAI0    :   1.15 : baseline ionised calcium (mmol/L)
KCITCA  :  70.0  : citrate load halving ionised calcium (mmol)
VMAXC   :   3.0  : maximal citrate metabolism (mmol/min)
KMC     :   3.0  : Km of citrate metabolism (mmol)
KCAEQ   :   0.25 : rate of approach to the citrate-set calcium (/min)
VCA     :  25.0  : effective distribution volume of ionised calcium (L)

// ---------------- temperature -------------------------------------------
CBODY   : 245.0  : body heat capacity (kJ/degC)
QMET    :   5.4  : metabolic heat production (kJ/min)
QENV    :   5.4  : baseline environmental loss (kJ/min)
EXPO    :   0.0  : additional loss when exposed in theatre (kJ/min)
TFCOLD  :   4.0  : temperature of unwarmed blood products (degC)
TFROOM  :  20.0  : temperature of unwarmed crystalloid (degC)
TWARM   :  38.0  : temperature delivered by a fluid warmer (degC)
WARMER  :   0.0  : 1 = fluid warmer in use

// ---------------- perfusion, oxygen debt, acid-base ---------------------
KCO     :   1.5  : fall in cardiac output per unit volume deficit
VO2REQ  : 300.0  : oxygen requirement at term (mL O2/min)
EXTMAX  :   0.70 : maximal oxygen extraction ratio
LDDEBT  : 8400.0 : lethal accumulated oxygen debt (mL O2, 120 mL/kg)
KLACP   :   0.0035 : lactate produced per unit oxygen-debt flux
KLACCL  :   0.023 : lactate clearance (/min) at normal liver flow
LAC0    :   1.2  : baseline lactate (mmol/L)
HCO0    :  20.0  : baseline bicarbonate at term (mmol/L)
KBUF    :   0.85 : bicarbonate consumed per unit lactate produced
CRYSTBIC:   0.0  : buffer content of the crystalloid (0 saline, 29 balanced)
UO0     :   1.5  : urine output at normal perfusion (mL/min)
KREF    :   0.006 : transcapillary refill per mL of deficit (/min)
REFMAX  :  12.0  : maximum transcapillary refill (mL/min)

// ---------------- non-drug interventions (times in min; 1e9 = never) ----
TOXYI   : 1e9    : start of the oxytocin infusion (min)
ROXYI   :   0.1667: oxytocin infusion rate (IU/min, = 40 IU / 4 h)
TOXYE   : 1e9    : end of the oxytocin infusion (min)
TMASS   : 1e9    : start of uterine massage / bimanual compression (min)
TMASSE  : 1e9    : massage stopped (min)
EMECH   :   1.2  : contractile drive contributed by massage
KMECHON :   0.6  : onset rate of the massage effect (/min)
KMECHOF :   0.25 : decay rate once massage stops (/min)
TBAL    : 1e9    : intrauterine balloon tamponade inserted (min)
EBAL    :   0.80 : fraction of patency occluded by the balloon
TAOC    : 1e9    : aortic compression started (min)
TAOCE   : 1e9    : aortic compression released (min)
EAOC    :   0.85 : fraction of uterine perfusion pressure removed
TREPAIR : 1e9    : laceration sutured (min)
TMANUAL : 1e9    : manual removal of retained tissue (min)
TUAL    : 1e9    : uterine artery ligation (min)
EUAL    :   0.45 : fraction of UBF0 removed by ligation
THYST   : 1e9    : hysterectomy (min)

// ---------------- resuscitation protocol --------------------------------
TRESUS  : 1e9    : time the resuscitation protocol starts (min)
RCRYST  :   0.0  : crystalloid infusion rate (mL/min)
TCRYSTE : 1e9    : crystalloid stopped (min)
HBTRIG  :   7.0  : transfuse if haemoglobin below this (g/dL)
MTPLOSS : 1e9    : activate MTP above this cumulative loss (mL)
QSTOP   :  25.0  : ... provided the leak still exceeds this (mL/min)
RBCRATE :   0.2  : red cells transfused (units/min) when triggered
MAXRBC  :  10.0  : red-cell units available
FFPRATIO:   0.0  : FFP units per red-cell unit
PLTTRIG :   0.0  : transfuse platelets below this count (10^9/L)
PLTRATE :   0.0  : platelet units/min when triggered
MAXPLT  :   2.0  : platelet units available
FIBTRIG :   0.0  : give fibrinogen below this concentration (g/L)
RFGC    :   0.0  : fibrinogen concentrate rate (g/min) when triggered
MAXFGC  :   6.0  : grams of fibrinogen concentrate available
TFGC    : 1e9    : time of an EMPIRIC fibrinogen dose (min)
DFGC    :   0.0  : size of that empiric dose (g, given over 10 min)
RCA     :   0.0  : calcium replacement (mmol/min) when triggered
CATRIG  :   0.0  : replace calcium below this ionised calcium (mmol/L)
MAXCA   :  20.0  : mmol of calcium available

$CMT @annotated
// ---- drug PK (10) ----
OXY   : oxytocin in the central compartment (IU)
CBT   : carbetocin in the central compartment (ug)
MSD   : misoprostol at the sublingual site (ug)
MSC   : misoprostol acid, central (ug)
ERD   : ergometrine at the intramuscular site (mg)
ERC   : ergometrine, central (mg)
PGD   : carboprost at the intramuscular site (ug)
PGC   : carboprost, central (ug)
TX1   : tranexamic acid, central (mg)
TX2   : tranexamic acid, peripheral (mg)
// ---- circulating volume and its contents (5) ----
V     : circulating blood volume (mL)
CIV   : crystalloid still intravascular (mL)
HBM   : circulating haemoglobin mass (g)
FBM   : circulating fibrinogen mass (g)
PLTM  : circulating platelets (10^9)
// ---- pooled coagulation activity (1) ----
FCT   : pooled clotting-factor activity times volume (activity.L)
// ---- uterus (4) ----
TONE  : myometrial tone (0-1)
OTR   : available oxytocin receptor fraction (0-1)
RET   : retained tissue fraction (0-1)
MECH  : mechanical (massage) contractile drive
// ---- haemostasis at the two bleeding sites (5) ----
CLT   : placental-bed clot occlusion (0-1)
CLTR  : trauma-site clot occlusion (0-1)
MAT   : thrombus organisation / durability of the seal (0-1)
TPA   : endothelial and decidual plasminogen activator (normalised)
PLS   : active plasmin (normalised, baseline 1)
// ---- internal milieu (5) ----
CIT   : citrate load (mmol)
CAI   : ionised calcium (mmol/L)
TMP   : core temperature (degC)
LAC   : lactate (mmol/L)
HCO   : bicarbonate (mmol/L)
// ---- integrated outcomes and counters (9) ----
LOSS  : cumulative blood loss (mL)
UO    : cumulative urine output (mL)
ODEBT : accumulated oxygen debt (mL O2)
TSEV  : cumulative minutes with MAP < 50 mmHg (Sheehan / organ risk)
TAKI  : cumulative minutes with MAP < 65 mmHg (AKI risk)
XSEV  : cumulative minutes with fibrinogen < 2 g/L
RBCU  : red-cell units transfused
FFPU  : FFP units transfused
PLTU  : platelet units transfused
FGCU  : grams of fibrinogen concentrate given
CAU   : calcium given (mmol)

$MAIN
// individual severity and receptor availability: the two ETAs act on the
// logit of the two knobs that matter, so both stay inside (0, 1).
ATNi = ATN;
if (ATN > 1e-6 && ATN < 1.0 - 1e-6) {
  double la = log(ATN / (1.0 - ATN)) + ETA(1);
  ATNi = 1.0 / (1.0 + exp(-la));
}
OTRi = OTR0;
if (OTR0 > 1e-6 && OTR0 < 1.0 - 1e-6) {
  double lo = log(OTR0 / (1.0 - OTR0)) + ETA(2);
  OTRi = 1.0 / (1.0 + exp(-lo));
}

TONE_0  = TONE0;
OTR_0   = OTRi;
RET_0   = RET0;
TPA_0   = 1.0;
PLS_0   = 1.0;
CAI_0   = CAI0;
TMP_0   = 37.0;
LAC_0   = LAC0;
HCO_0   = HCO0;
// the separation-phase loss is already banked at t = 0
V_0     = V0 - LOSS0;
LOSS_0  = LOSS0;
double f0 = 1.0 - LOSS0 / V0;
HBM_0   = HB0  * V0 / 100.0  * f0;
FBM_0   = FIB0 * V0 / 1000.0 * f0;
PLTM_0  = PLT0 * V0 / 1000.0 * f0;
FCT_0   = FCT0 * V0 / 1000.0 * f0;

$GLOBAL
#define HILL(x, k, n) (pow(x, n) / (pow(k, n) + pow(x, n)))
#define CLAMP(x, lo, hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))
// derived quantities shared between $ODE and $TABLE
double Hb_, FIB_, PLT_, FACT_, DEF_, pH_, PaCO2_, MAPn, COn, HRn, SIn;
double DO2n, DEBTF, M_PERF, M_LIV, GFRn, M_CA, M_PHC, M_TMPC, M_TMPT, M_PHT;
double Coxy, Ccbt, Cerg, Cpgf, Cmso, Ctxa, DRV_, TONETGT, PATn, SEALn;
double UBFn, Qut_, Qtr_, Qleak, fSH, fFIB, fPLT, fFCT, KFORME, PLSEFF;
double TXAp, TXAa, Rrbc, Rffp, Rplt, Rfgc, Rca, Rcry, Vprod;
double ATNi, OTRi;

$ODE
double Vsafe = (V < 1500.0 ? 1500.0 : V);
double VL    = Vsafe / 1000.0;
Hb_   = HBM  / (Vsafe / 100.0);
FIB_  = FBM  / VL;
PLT_  = PLTM / VL;
FACT_ = FCT  / VL;
DEF_  = CLAMP((V0 - Vsafe) / V0, 0.0, 0.95);

// ---------------- acid-base ---------------------------------------------
double HCOs = (HCO < 2.0 ? 2.0 : HCO);
PaCO2_ = CLAMP(10.0 + HCOs, 16.0, 40.0);
pH_    = 6.1 + log10(HCOs / (0.03 * PaCO2_));

// ---------------- modifiers on tone and on clot -------------------------
double CAIs = (CAI < 0.2 ? 0.2 : CAI);
M_CA   = HILL(CAIs, 0.55, 2.0) / HILL(1.15, 0.55, 2.0);
M_PHC  = (1.0 / (1.0 + exp(-(pH_ - 7.05) / 0.07))) /
         (1.0 / (1.0 + exp(-(7.44 - 7.05) / 0.07)));
M_TMPC = (TMP < 37.0 ? exp(-0.12 * (37.0 - TMP)) : 1.0);
M_TMPT = CLAMP(1.0 - 0.07 * (37.0 - TMP), 0.30, 1.0);
M_PHT  = CLAMP(1.0 - 0.55 * (7.35 - pH_) / 0.35, 0.25, 1.0);

// ---------------- haemodynamics -----------------------------------------
double fphv = (1.0 / (1.0 + exp(-(pH_ - 7.05) / 0.08))) /
              (1.0 / (1.0 + exp(-(7.44 - 7.05) / 0.08)));
double M_VASO = (0.55 + 0.45 * fphv) * (0.85 + 0.15 * CAIs / CAI0);
double D50 = 0.36 * M_VASO;
MAPn = MAP0 / (1.0 + pow(DEF_ / (D50 < 0.05 ? 0.05 : D50), 6.0));
COn  = CO0 * CLAMP(1.0 - KCO * DEF_, 0.10, 1.05);
HRn  = CLAMP(85.0 + 190.0 * pow(DEF_, 1.2), 85.0, 175.0);
SIn  = HRn / (MAPn * 1.35 < 20.0 ? 20.0 : MAPn * 1.35);
DO2n = COn * Hb_ * 1.34 * 0.98 * 10.0;
DEBTF = VO2REQ - EXTMAX * DO2n;
if (DEBTF < 0.0) DEBTF = 0.0;

M_PERF = HILL(MAPn, 42.0, 3.0) / HILL(MAP0, 42.0, 3.0);
M_LIV  = CLAMP(HILL(MAPn, 45.0, 2.0) / HILL(MAP0, 45.0, 2.0), 0.15, 1.0);
GFRn   = CLAMP(HILL(MAPn, 55.0, 4.0) / HILL(MAP0, 55.0, 4.0), 0.0, 1.2);

// ---------------- drug concentrations and contractile drive -------------
Coxy = OXY / V_OXY;
Ccbt = CBT / V_CBT;
Cerg = ERC * 1000.0 / V_ERG;
Cpgf = PGC / V_PGF;
Cmso = MSC / V_MSO;
Ctxa = TX1 / V1_TX;

double OTRs = CLAMP(OTR, 0.02, 1.0);
double EMO = EMAX_OXY * (OTRs / 0.4 < 1.0 ? OTRs / 0.4 : 1.0);
double EMC = EMAX_CBT * (OTRs / 0.4 < 1.0 ? OTRs / 0.4 : 1.0);
double D_oxy = EMO * Coxy / (EC50_OXY / OTRs + Coxy);
double D_cbt = EMC * Ccbt / (EC50_CBT / OTRs + Ccbt);
double D_erg = EMAX_ERG * Cerg / (EC50_ERG + Cerg);
double D_pgf = EMAX_PGF * Cpgf / (EC50_PGF + Cpgf);
double D_mso = EMAX_MSO * Cmso / (EC50_MSO + Cmso);

double DRVI_ = DRVMAX * (1.0 - ATNi);
double CAP_  = 1.0 - KCAP * ATNi;
DRV_ = DRVI_ + CAP_ * (D_oxy + D_cbt + D_erg + D_pgf + D_mso + MECH);
double TMAXn = CLAMP(1.0 - RET, 0.0, 1.0);
TONETGT = TMAXn * CAP_ * HILL(DRV_, KDRV, HDRV) * M_CA * M_TMPT * M_PHT * M_PERF;

// ---------------- the leak ----------------------------------------------
double hyst = (SOLVERTIME >= THYST ? 1.0 : 0.0);
double aoc  = ((SOLVERTIME >= TAOC && SOLVERTIME < TAOCE) ? EAOC : 0.0);
double ual  = (SOLVERTIME >= TUAL ? EUAL : 0.0);
double bal  = (SOLVERTIME >= TBAL ? EBAL * CLAMP((SOLVERTIME - TBAL) / 3.0, 0.0, 1.0) : 0.0);
UBFn = UBF0 * (MAPn / MAP0) * (1.0 - aoc) * (1.0 - ual);

double TONEs = CLAMP(TONE, 0.0, 1.0);
double CLTs  = CLAMP(CLT, 0.0, 1.0);
double MATs  = CLAMP(MAT, 0.0, 1.0);
PATn  = pow(1.0 - TONEs, NP);
// a clot only seals a sinus the myometrium (or an organised thrombus) holds shut
SEALn = (CSEAL0 + (1.0 - CSEAL0) * (TONEs > MATs ? TONEs : MATs)) * CLTs;
Qut_  = UBFn * PATn * (1.0 - SEALn) * (1.0 - bal);
double rep = (SOLVERTIME >= TREPAIR ? 1.0 : 0.0);
Qtr_  = KTR * (MAPn / MAP0) * (1.0 - CLAMP(CLTR, 0.0, 1.0)) * (1.0 - rep);
Qleak = (Qut_ + Qtr_) * (1.0 - hyst);

// ---------------- clot formation, the shear gate, lysis -----------------
fFIB = HILL((FIB_ > 0.0 ? FIB_ : 0.0), KFIB, 3.0);
fPLT = HILL((PLT_ > 0.0 ? PLT_ : 0.0), KPLT, 2.0);
fFCT = HILL((FACT_ > 0.0 ? FACT_ : 0.0), KFCT, 2.0);
fSH  = 1.0 / (1.0 + pow(Qleak / QS50, 2.0));
KFORME = KFORM * fFIB * fPLT * fFCT * M_CA * M_PHC * M_TMPC;
PLSEFF = (PLS > 1.0 ? PLS - 1.0 : 0.0);
TXAp = ETXPROT * Ctxa / (Ctxa + ITX50);
TXAa = ETXACT  * Ctxa / (Ctxa + ITX50);

// ---------------- product and fluid delivery ----------------------------
int resus = (SOLVERTIME >= TRESUS ? 1 : 0);
int need  = ((Hb_ < HBTRIG) || (LOSS > MTPLOSS && Qleak > QSTOP)) ? 1 : 0;
Rrbc = (resus && need && RBCU < MAXRBC)                ? RBCRATE : 0.0;
Rffp = Rrbc * FFPRATIO;
Rplt = (resus && PLT_ < PLTTRIG && PLTU < MAXPLT)       ? PLTRATE : 0.0;
Rfgc = (resus && FIB_ < FIBTRIG && FGCU < MAXFGC)      ? RFGC : 0.0;
if (SOLVERTIME >= TFGC && SOLVERTIME < TFGC + 10.0) Rfgc = Rfgc + DFGC / 10.0;
Rca  = (resus && CAI < CATRIG && CAU < MAXCA)            ? RCA : 0.0;
Rcry = ((SOLVERTIME >= TRESUS && SOLVERTIME < TCRYSTE)              ? RCRYST : 0.0);

Vprod = Rrbc * RBC_VOL + Rffp * FFP_VOL + Rplt * PLT_VOL +
        Rfgc / FGC_FIB * FGC_VOL;
double Vin    = Vprod + Rcry;
double refill = KREF * (V0 - Vsafe > 0.0 ? V0 - Vsafe : 0.0);
if (refill > REFMAX) refill = REFMAX;
double uo = UO0 * GFRn;

// ================= DIFFERENTIAL EQUATIONS ==============================
// ---- drug PK ----
dxdt_OXY = -CL_OXY / V_OXY * OXY +
           ((SOLVERTIME >= TOXYI && SOLVERTIME < TOXYE) ? ROXYI : 0.0);
dxdt_CBT = -KE_CBT * CBT;
dxdt_MSD = -KA_MSO * MSD;
dxdt_MSC =  KA_MSO * MSD * F_MSO - KE_MSO * MSC;
dxdt_ERD = -KA_ERG * ERD;
dxdt_ERC =  KA_ERG * ERD * F_ERG - KE_ERG * ERC;
dxdt_PGD = -KA_PGF * PGD;
dxdt_PGC =  KA_PGF * PGD * F_PGF - KE_PGF * PGC;
dxdt_TX1 = -(KE_TX + K12_TX) * TX1 + K21_TX * TX2;
dxdt_TX2 =  K12_TX * TX1 - K21_TX * TX2;

// ---- volume and its contents (masses, so dilution is automatic) ----
dxdt_V   = -Qleak + Vin + refill - uo;
dxdt_CIV =  Rcry - 0.035 * CIV - (Qleak / Vsafe) * CIV;
dxdt_HBM = -(Qleak / Vsafe) * HBM + Rrbc * RBC_HB;
dxdt_FBM = -(Qleak / Vsafe) * FBM
           - KFDEG * PLSEFF * FBM * (1.0 - TXAp)
           - 2.0 * KFORME * fSH * (1.0 - CLTs)
           + 0.0014
           + Rffp * FFP_FIB + Rplt * PLT_FIB + Rfgc;
dxdt_PLTM = -(Qleak / Vsafe) * PLTM
            - 200.0 * KFORME * fSH * (1.0 - CLTs)
            + Rplt * PLT_CNT;
dxdt_FCT  = -(Qleak / Vsafe) * FCT
            - 0.30 * KFORME * fSH * (1.0 - CLTs)
            + Rffp * FFP_FCT * FFP_VOL / 1000.0;

// ---- uterus ----
dxdt_TONE = KON * (TONETGT - TONE);
double occ = Coxy / (EC50_OXY + Coxy) + Ccbt / (EC50_CBT + Ccbt);
if (occ > 1.0) occ = 1.0;
dxdt_OTR  = KREC * (1.0 - OTR) - KDES * occ * OTR;
dxdt_RET  = (SOLVERTIME >= TMANUAL ? -2.0 * RET : 0.0);
dxdt_MECH = ((SOLVERTIME >= TMASS && SOLVERTIME < TMASSE)
             ? KMECHON * (EMECH - MECH) : -KMECHOF * MECH);

// ---- haemostasis ----
dxdt_CLT  = KFORME * fSH * (1.0 - CLTs) - KLYS * PLSEFF * CLTs * (1.0 - TXAp);
dxdt_MAT  = KMAT * fSH * CLTs * (1.0 - MATs) - 0.03 * PLSEFF * MATs * (1.0 - TXAp);
dxdt_CLTR = 0.6 * KFORME * (1.0 - CLAMP(CLTR, 0.0, 1.0))
            - KLYS * PLSEFF * CLAMP(CLTR, 0.0, 1.0) * (1.0 - TXAp);
dxdt_TPA  = KTPAB + KTPAINJ * (Qleak / UBF0)
            + KTPASHK * pow((MAPn < 65.0 ? (65.0 - MAPn) / 65.0 : 0.0), 1.5)
            - KTPACL * TPA;
dxdt_PLS  = KACT * TPA * (FACT_ > 0.0 ? FACT_ : 0.0) * (1.0 - TXAa)
            - KPLSCL * PLS;

// ---- citrate, calcium, temperature, perfusion markers ----
double Rcit = Rrbc * RBC_CIT + Rffp * FFP_CIT + Rplt * PLT_CIT;
dxdt_CIT = Rcit - VMAXC * (CIT / (KMC + CIT)) * M_LIV;
double CAIeq = CAI0 / (1.0 + (CIT > 0.0 ? CIT : 0.0) / KCITCA);
dxdt_CAI = KCAEQ * (CAIeq - CAI) + Rca / VCA;

double Tfp = (WARMER > 0.5 ? TWARM : TFCOLD);
double Tfc = (WARMER > 0.5 ? TWARM : TFROOM);
double Qfl = 0.00418 * (Vprod * (TMP - Tfp) + Rcry * (TMP - Tfc));
dxdt_TMP = (QMET * (0.7 + 0.3 * M_PERF) - QENV - EXPO - Qfl) / CBODY;

double lprod = KLACP * DEBTF;
dxdt_LAC = lprod + KLACCL * LAC0 - KLACCL * LAC * M_LIV;
dxdt_HCO = -KBUF * lprod + Rcry * (CRYSTBIC - HCO) / Vsafe;

// ---- outcomes and counters ----
dxdt_LOSS  = Qleak;
dxdt_UO    = uo;
dxdt_ODEBT = DEBTF;
dxdt_TSEV  = 1.0 / (1.0 + exp((MAPn - 50.0) / 1.5));
dxdt_TAKI  = 1.0 / (1.0 + exp((MAPn - 65.0) / 1.5));
dxdt_XSEV  = 1.0 / (1.0 + exp((FIB_ - 2.0) / 0.05));
dxdt_RBCU  = Rrbc;
dxdt_FFPU  = Rffp;
dxdt_PLTU  = Rplt;
dxdt_FGCU  = Rfgc;
dxdt_CAU   = Rca;

$TABLE
capture AtonySeverity = ATNi;
capture OTRstart = OTRi;
capture Hb      = Hb_;
capture Fib     = FIB_;
capture Plt     = PLT_;
capture Factors = FACT_;
capture MAP     = MAPn;
capture CO      = COn;
capture HR      = HRn;
capture ShockIx = SIn;
capture pH      = pH_;
capture Deficit = DEF_;
capture Qbleed  = Qleak;
capture Quterus = Qut_;
capture Qtrauma = Qtr_;
capture Patency = PATn;
capture Seal    = SEALn;
capture ShearGate = fSH;
capture ToneTgt = TONETGT;
capture Drive   = DRV_;
capture DO2     = DO2n;
capture DebtFlux= DEBTF;
capture Oxytocin= Coxy;
capture Carbopr = Cpgf;
capture Ergomet = Cerg;
capture Misopr  = Cmso;
capture TXA     = Ctxa;
capture TXAprot = TXAp;
capture Plasmin = PLS;
capture UrineRate = UO0 * GFRn;
capture Uflow   = UBFn;
capture PPH1000 = (LOSS >= 1000.0 ? 1.0 : 0.0);
capture PPH1500 = (LOSS >= 1500.0 ? 1.0 : 0.0);
capture Exsang  = (DEF_ > 0.45 ? 1.0 : 0.0);
capture Lethal  = ((DEF_ > 0.45 || ODEBT > LDDEBT) ? 1.0 : 0.0);

$CAPTURE @annotated
ATN     : population atony severity for this subject
OTR0    : population receptor availability for this subject

$OMEGA @annotated @block
// between-subject variability, added to the logit/log of the two knobs that
// matter most: contractile capacity and receptor availability
eATN  : 0.09 : variability in atony severity (logit scale)
eOTR  : 0.04 0.16 : variability in receptor availability (logit scale)

$SIGMA 0
'

pph <- mcode("pph", pph_code, atol = 1e-8, rtol = 1e-8, maxsteps = 100000)

## =====================================================================
##  DOSING / INTERVENTION BUILDERS
## =====================================================================

## uterotonic events. Doses are the ones in the WHO / ACOG / RCOG bundles.
ev_oxy_proph <- function(t = 1)   ev(time = t,  amt = 10,   cmt = "OXY")   # 10 IU IM
ev_oxy_treat <- function(t = 5)   ev(time = t,  amt = 5,    cmt = "OXY")   # 5 IU slow IV
ev_carbetocin<- function(t = 5)   ev(time = t,  amt = 100,  cmt = "CBT")   # 100 ug IV
ev_ergometrine<- function(t = 15) ev(time = t,  amt = 0.5,  cmt = "ERD")   # 0.5 mg IM
ev_carboprost<- function(t = 20)  ev(time = t,  amt = 250,  cmt = "PGD")   # 250 ug IM
ev_misoprostol<- function(t = 5, amt = 800) ev(time = t, amt = amt, cmt = "MSD")
ev_txa       <- function(t = 10)  ev(time = t,  amt = 1000, cmt = "TX1")   # 1 g IV

## the parameter blocks that describe non-drug care
p_oxy_infusion <- list(TOXYI = 5, TOXYE = 245, ROXYI = 0.1667)   # 40 IU / 4 h
p_massage      <- list(TMASS = 4, TMASSE = 90)
p_goal_directed <- list(TRESUS = 12, RCRYST = 12, TCRYSTE = 40,
                        HBTRIG = 7, MTPLOSS = 1500, QSTOP = 25,
                        RBCRATE = 0.2, FFPRATIO = 1.0,
                        PLTTRIG = 75, PLTRATE = 0.05,
                        FIBTRIG = 2.0, RFGC = 0.4,
                        RCA = 0.6, CATRIG = 1.0, WARMER = 1)
p_crystalloid_first <- list(TRESUS = 12, RCRYST = 70, TCRYSTE = 90,
                            HBTRIG = 7, MTPLOSS = 1500, QSTOP = 25,
                            RBCRATE = 0.2, FFPRATIO = 0, FIBTRIG = 0,
                            RCA = 0, WARMER = 0)

run <- function(..., par = list(), events = NULL, end = 240, delta = 1) {
  extra <- list(...)
  p <- modifyList(par, extra)
  m <- param(pph, p)
  if (is.null(events)) events <- ev(time = 0, amt = 0, cmt = "OXY")
  mrgsim(m, events = events, end = end, delta = delta) %>% as_tibble()
}

summarise_run <- function(d) {
  dead <- d %>% filter(Lethal > 0) %>% slice(1)
  tibble(
    loss_mL    = max(d$LOSS),
    controlled = nrow(dead) == 0 && tail(d$Qbleed, 1) < 25,
    death_min  = if (nrow(dead)) dead$time[1] else NA_real_,
    Hb_min     = min(d$Hb),
    Fib_min    = min(d$Fib),
    MAP_min    = min(d$MAP),
    tone_end   = tail(d$TONE, 1),
    RBC_units  = max(d$RBCU),
    FFP_units  = max(d$FFPU),
    Fg_g       = max(d$FGCU),
    iCa_min    = min(d$CAI),
    temp_min   = min(d$TMP),
    pH_min     = min(d$pH),
    O2debt     = max(d$ODEBT),
    min_MAP_lt50 = max(d$TSEV),
    min_Fib_lt2  = max(d$XSEV)
  )
}

## =====================================================================
##  THE 20 SCENARIOS
##  Losses quoted in the comments are from the reference implementation
##  described in the header (see PROVENANCE).
## =====================================================================

## --- block 1: natural history — one parameter spans the whole disease ---
## ATN 0.00 ->  270 mL   ATN 0.63 -> DEAD 70 min   ATN 0.88 -> DEAD 12 min
## ATN 0.50 ->  827 mL   ATN 0.74 -> DEAD 22 min   ATN 0.94 -> DEAD 10 min
## ATN 0.60 -> 2218 mL                             ATN 1.00 -> DEAD 10 min
S01_natural_history <- function() {
  lapply(c(0, 0.5, 0.6, 0.63, 0.74, 0.88, 0.94, 1.0), function(a)
    summarise_run(run(ATN = a)) %>% mutate(ATN = a)) %>% bind_rows()
}

## --- block 2: prophylaxis (AMTSL) on a uterus destined for atony --------
## nothing 2218 · oxytocin 10 IU 1212 · carbetocin 946 · misoprostol 1118
S02_prophylaxis <- function() bind_rows(
  summarise_run(run(ATN = 0.6)) %>% mutate(arm = "none"),
  summarise_run(run(ATN = 0.6, events = ev_oxy_proph(1))) %>% mutate(arm = "oxytocin 10 IU IM"),
  summarise_run(run(ATN = 0.6, events = ev_carbetocin(1))) %>% mutate(arm = "carbetocin 100 ug"),
  summarise_run(run(ATN = 0.6, events = ev_misoprostol(1, 600))) %>% mutate(arm = "misoprostol 600 ug")
)

## --- block 3: who needs what (the escalation ladder) -------------------
## ATN 0.60: any single agent suffices (massage 940, oxytocin 1078)
## ATN 0.70: oxytocin ALONE fails (D@77); + massage 1446
## ATN 0.74: massage alone D@58, oxytocin alone D@51, TOGETHER 2054 (survives)
## ATN 0.78: only the four-agent ladder holds (1877)
## ATN 0.84+: pharmacology cannot reach the threshold at all
S03_ladder <- function() {
  arms <- list(
    untreated   = function(a) run(ATN = a),
    massage     = function(a) run(ATN = a, par = p_massage),
    oxytocin    = function(a) run(ATN = a, par = p_oxy_infusion,
                                  events = ev_oxy_treat()),
    oxy_massage = function(a) run(ATN = a, par = c(p_oxy_infusion, p_massage),
                                  events = ev_oxy_treat()),
    full_ladder = function(a) run(ATN = a, par = c(p_oxy_infusion, p_massage),
                                  events = ev_oxy_treat() + ev_ergometrine(15) +
                                           ev_carboprost(20) + ev_carboprost(35))
  )
  expand.grid(ATN = c(0.60, 0.63, 0.66, 0.70, 0.74, 0.78, 0.84, 0.90),
              arm = names(arms), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(res = list(summarise_run(arms[[arm]](ATN)))) %>%
    tidyr::unnest(res) %>% ungroup()
}

## --- block 4: receptor desensitisation --------------------------------
## severe atony (ATN 0.74) on oxytocin + massage:
##   OTR 1.00 -> 2054 · 0.60 -> 2316 · 0.38 -> 2605 · 0.25 -> 2923 mL
## from OTR 0.38: double dose -143, quadruple -277,
##   ergometrine -824, carboprost -925, misoprostol -652 mL
S04_desensitisation <- function() {
  base <- function(otr, ...) run(ATN = 0.74, OTR0 = otr,
                                 par = c(p_oxy_infusion, p_massage), ...)
  bind_rows(
    lapply(c(1.0, 0.60, 0.38, 0.25), function(o)
      summarise_run(base(o, events = ev_oxy_treat())) %>%
        mutate(arm = sprintf("OTR0 %.2f", o))) %>% bind_rows(),
    summarise_run(run(ATN = 0.74, OTR0 = 0.38,
                      par = c(p_massage, list(TOXYI = 5, TOXYE = 245, ROXYI = 0.3334)),
                      events = ev(time = 5, amt = 10, cmt = "OXY"))) %>%
      mutate(arm = "OTR0 0.38 + double oxytocin"),
    summarise_run(run(ATN = 0.74, OTR0 = 0.38,
                      par = c(p_massage, list(TOXYI = 5, TOXYE = 245, ROXYI = 0.6668)),
                      events = ev(time = 5, amt = 20, cmt = "OXY"))) %>%
      mutate(arm = "OTR0 0.38 + quadruple oxytocin"),
    summarise_run(base(0.38, events = ev_oxy_treat() + ev_ergometrine(15))) %>%
      mutate(arm = "OTR0 0.38 + ergometrine at 15"),
    summarise_run(base(0.38, events = ev_oxy_treat() + ev_carboprost(15))) %>%
      mutate(arm = "OTR0 0.38 + carboprost at 15"),
    summarise_run(base(0.38, events = ev_oxy_treat() + ev_misoprostol(15))) %>%
      mutate(arm = "OTR0 0.38 + misoprostol at 15")
  )
}

## --- block 5: TXA timing ----------------------------------------------
## no TXA 4099 · 5 min 3400 · 15 min 3645 · 30 min 4030 · >=45 min 4099 mL
S05_txa_timing <- function() {
  lapply(c(NA, 5, 15, 30, 45, 90, 180), function(tt) {
    e <- ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20) + ev_carboprost(35)
    if (!is.na(tt)) e <- e + ev_txa(tt) + ev_txa(tt + 30)
    summarise_run(run(ATN = 0.88,
                      par = c(p_oxy_infusion, p_massage, p_goal_directed),
                      events = e)) %>%
      mutate(txa_min = tt)
  }) %>% bind_rows()
}

## --- block 6: resuscitation strategy ----------------------------------
## none D@29 · crystalloid-first 3953 (Fib 1.37, 34.6 C, pH 7.31)
## goal-directed 3489 · no warmer 4021 · no FFP 3460
## empiric 4 g fibrinogen 3473 (vs 3489 — the OBS2 null result)
S06_resuscitation <- function() {
  e <- ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20) +
       ev_carboprost(35) + ev_txa(10) + ev_txa(40)
  base <- c(p_oxy_infusion, p_massage)
  bind_rows(
    summarise_run(run(ATN = 0.88, par = base, events = e)) %>% mutate(arm = "no resuscitation"),
    summarise_run(run(ATN = 0.88, par = c(base, p_crystalloid_first), events = e)) %>%
      mutate(arm = "crystalloid-first, RBC only"),
    summarise_run(run(ATN = 0.88, par = c(base, p_crystalloid_first, list(CRYSTBIC = 29)),
                      events = e)) %>% mutate(arm = "crystalloid-first, balanced"),
    summarise_run(run(ATN = 0.88, par = c(base, p_crystalloid_first,
                                          list(FIBTRIG = 2.0, RFGC = 0.4)), events = e)) %>%
      mutate(arm = "crystalloid-first + fibrinogen"),
    summarise_run(run(ATN = 0.88, par = c(base, p_goal_directed), events = e)) %>%
      mutate(arm = "goal-directed"),
    summarise_run(run(ATN = 0.88, par = c(base, p_goal_directed, list(WARMER = 0)),
                      events = e)) %>% mutate(arm = "goal-directed, no warmer"),
    summarise_run(run(ATN = 0.88, par = c(base, p_goal_directed, list(FFPRATIO = 0)),
                      events = e)) %>% mutate(arm = "goal-directed, no FFP"),
    summarise_run(run(ATN = 0.88, par = c(base, p_goal_directed,
                                          list(TFGC = 15, DFGC = 4)), events = e)) %>%
      mutate(arm = "goal-directed + empiric 4 g Fg (OBS2)")
  )
}

## --- block 7: mechanical / surgical escalation (ATN 0.94) -------------
## bundle only DEAD@80 (10044 mL) · balloon 20 min 4528 · 45 min 7617
## balloon 75 min DEAD@150 · aortic bridge + balloon 2782 (no transfusion)
## uterine artery ligation at 50 -> 9624 · hysterectomy at 60 -> 8451
S07_mechanical <- function() {
  e <- ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20) +
       ev_carboprost(35) + ev_txa(10) + ev_txa(40)
  base <- c(p_oxy_infusion, p_massage, p_goal_directed)
  bind_rows(
    summarise_run(run(ATN = 0.94, par = base, events = e)) %>% mutate(arm = "bundle only"),
    lapply(c(20, 45, 75), function(t)
      summarise_run(run(ATN = 0.94, par = c(base, list(TBAL = t)), events = e)) %>%
        mutate(arm = sprintf("balloon at %d min", t))) %>% bind_rows(),
    summarise_run(run(ATN = 0.94, par = c(base, list(TAOC = 8, TAOCE = 35, TBAL = 35)),
                      events = e)) %>% mutate(arm = "aortic compression 8-35 + balloon 35"),
    summarise_run(run(ATN = 0.94, par = c(base, list(TUAL = 50)), events = e)) %>%
      mutate(arm = "uterine artery ligation at 50"),
    summarise_run(run(ATN = 0.94, par = c(base, list(THYST = 60)), events = e)) %>%
      mutate(arm = "hysterectomy at 60"),
    summarise_run(run(ATN = 0.94, par = c(base, list(THYST = 120)), events = e)) %>%
      mutate(arm = "hysterectomy at 120")
  )
}

## --- block 8: what the products themselves cost ------------------------
## with Ca + warmer 9243 mL, 11 units, iCa 0.53, 36.9 C
## no calcium      10165 mL, 13 units, iCa 0.45
## no warmer       12116 mL, 16 units, 31.8 C
S08_transfusion_cost <- function() {
  e <- ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20) +
       ev_carboprost(35) + ev_txa(10) + ev_txa(40)
  base <- c(p_oxy_infusion, p_massage, p_goal_directed,
            list(TBAL = 60, MAXRBC = 16, MAXFGC = 8))
  bind_rows(
    summarise_run(run(ATN = 0.94, par = base, events = e)) %>% mutate(arm = "Ca + warmer"),
    summarise_run(run(ATN = 0.94, par = c(base, list(RCA = 0)), events = e)) %>%
      mutate(arm = "no calcium"),
    summarise_run(run(ATN = 0.94, par = c(base, list(WARMER = 0)), events = e)) %>%
      mutate(arm = "no warmer"),
    summarise_run(run(ATN = 0.94, par = c(base, list(RCA = 0, WARMER = 0)), events = e)) %>%
      mutate(arm = "neither")
  )
}

## --- block 9: the other three Ts --------------------------------------
## laceration 260 mL/min: uterotonics only DEAD@79; full escalation DEAD@80;
##   repair at 25 min 3794 mL; repair at 60 min 8179 mL
## retained tissue 35%: full ladder DEAD@106; manual removal at 30 -> 4407
## fibrinogen 1.6 g/L at presentation: 2193 mL, or 1677 with replacement
S09_other_Ts <- function() {
  e0 <- ev_oxy_treat()
  base <- c(p_oxy_infusion, p_massage, p_goal_directed)
  bind_rows(
    summarise_run(run(ATN = 0.30, KTR = 260, par = base, events = e0)) %>%
      mutate(arm = "laceration, uterotonics only"),
    summarise_run(run(ATN = 0.30, KTR = 260, par = base,
                      events = e0 + ev_ergometrine(15) + ev_carboprost(20))) %>%
      mutate(arm = "laceration + full escalation"),
    summarise_run(run(ATN = 0.30, KTR = 260, par = c(base, list(TREPAIR = 25)),
                      events = e0)) %>% mutate(arm = "laceration, repair at 25"),
    summarise_run(run(ATN = 0.30, KTR = 260, par = c(base, list(TREPAIR = 60)),
                      events = e0)) %>% mutate(arm = "laceration, repair at 60"),
    summarise_run(run(ATN = 0.74, RET0 = 0.35, par = base,
                      events = e0 + ev_ergometrine(15) + ev_carboprost(20) +
                               ev_carboprost(35))) %>%
      mutate(arm = "retained tissue 35%"),
    summarise_run(run(ATN = 0.74, RET0 = 0.35, par = c(base, list(TMANUAL = 30)),
                      events = e0 + ev_ergometrine(15) + ev_carboprost(20) +
                               ev_carboprost(35))) %>%
      mutate(arm = "retained tissue, removal at 30"),
    summarise_run(run(ATN = 0.74, FIB0 = 1.6,
                      par = c(base, list(FIBTRIG = 0, RFGC = 0)),
                      events = e0 + ev_ergometrine(15) + ev_carboprost(20))) %>%
      mutate(arm = "fibrinogen 1.6, no replacement"),
    summarise_run(run(ATN = 0.74, FIB0 = 1.6,
                      par = c(base, list(FIBTRIG = 2.5, RFGC = 0.4)),
                      events = e0 + ev_ergometrine(15) + ev_carboprost(20))) %>%
      mutate(arm = "fibrinogen 1.6, replaced to 2.5")
  )
}

## --- block 10: which factor crosses its critical value first ----------
## crystalloid-first arm: fibrinogen <2.0 g/L at 3476 mL, Hb <7 at 3555 mL,
## platelets never (<50), pooled factors never (<0.30).
## Ratio critical/initial: fibrinogen 0.44 > factors 0.30 > platelets 0.20
S10_dilution_order <- function() {
  d <- run(ATN = 0.88,
           par = c(p_oxy_infusion, p_massage, p_crystalloid_first),
           events = ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20) +
                    ev_carboprost(35) + ev_txa(10) + ev_txa(40))
  crossing <- function(col, crit) {
    i <- which(d[[col]] < crit)
    if (!length(i)) return(tibble(variable = col, crossed_at_mL = NA, nadir = min(d[[col]])))
    tibble(variable = col, crossed_at_mL = d$LOSS[i[1]], nadir = min(d[[col]]))
  }
  bind_rows(crossing("Fib", 2.0), crossing("Plt", 50), crossing("Factors", 0.30),
            crossing("Hb", 7.0))
}

## --- block 11: the reserve arithmetic ---------------------------------
## total atony: peak 516 mL/min; 500 mL at 2.0 min, 1000 at 4.0,
## 2100 (class III) at 6.0, exsanguination at 10 min.
## 2100 mL of pregnancy volume expansion / 750 mL/min = 2.8 min.
S11_reserve <- function() {
  d <- run(ATN = 1.0, end = 30, delta = 0.5)
  tibble(peak_leak = max(d$Qbleed),
         t500 = d$time[which(d$LOSS >= 500)[1]],
         t1000 = d$time[which(d$LOSS >= 1000)[1]],
         t2100 = d$time[which(d$LOSS >= 2100)[1]],
         t_lethal = d$time[which(d$Lethal > 0)[1]],
         reserve_min = 2100 / 750)
}

## --- block 12: population simulation ----------------------------------
## 500 subjects with variability in contractile capacity and receptor
## availability, all on the WHO first-line bundle: the distribution of blood
## loss, the PPH and severe-PPH rates, and the refractory fraction.
S12_population <- function(n = 500, seed = 20260804) {
  set.seed(seed)
  idata <- tibble(
    ID   = seq_len(n),
    ATN  = plogis(qlogis(0.55) + rnorm(n, 0, 0.9)),
    OTR0 = plogis(qlogis(0.75) + rnorm(n, 0, 0.8))
  )
  m <- zero_re(param(pph, c(p_oxy_infusion, p_massage)))
  out <- mrgsim(m, idata = idata,
                events = ev_oxy_treat() + ev_ergometrine(15) + ev_carboprost(20),
                end = 240, delta = 2) %>% as_tibble()
  out %>% group_by(ID) %>%
    summarise(loss = max(LOSS), lethal = max(Lethal), .groups = "drop") %>%
    summarise(median_loss = median(loss),
              p_PPH_1000 = mean(loss >= 1000),
              p_severe_1500 = mean(loss >= 1500),
              p_massive_2500 = mean(loss >= 2500),
              p_lethal = mean(lethal > 0))
}

## =====================================================================
##  Run everything (uncomment to execute)
## =====================================================================
# print(S01_natural_history())
# print(S02_prophylaxis())
# print(S03_ladder())
# print(S04_desensitisation())
# print(S05_txa_timing())
# print(S06_resuscitation())
# print(S07_mechanical())
# print(S08_transfusion_cost())
# print(S09_other_Ts())
# print(S10_dilution_order())
# print(S11_reserve())
# print(S12_population())

## =====================================================================
##  SENSITIVITY / FALSIFICATION NOTES
##  ---------------------------------
##  The three structural claims that most deserve attack, and how to attack
##  them with this file:
##
##  1. NP = 4 (the Poiseuille valve).  Set NP = 2 (area-proportional) and the
##     knife edge disappears: tone becomes a graded, forgiving variable and
##     the model loses its explanation for why "firm vs boggy" is a usable
##     bedside distinction.  Set NP = 1 and a normal third stage bleeds
##     3 L/h.  If real placental-bed flow is closer to linear in tone, the
##     model is wrong about the whole shape of the disease.
##
##  2. CSEAL0 = 0.35 (a clot over an uncompressed sinus).  Raise it to 1.0
##     and severe atony becomes self-limiting: the bed thromboses shut
##     regardless of tone, uterotonics stop mattering, and refractory PPH
##     ceases to exist.  Lower it to 0 and no amount of coagulation support
##     ever helps.  This single parameter is the model's claim about why PPH
##     is a mechanical disease, so it is the one to measure.
##
##  3. KDES/KREC (receptor turnover).  If OTR recovery is much faster than
##     t1/2 2.6 h, dose escalation should work and the model's advice to
##     cross receptor classes is wrong.  The prediction is testable in vitro
##     on myometrial strips from augmented and non-augmented labours.
##
##  Two further quantitative predictions worth checking against registry
##  data: that time-to-tamponade should dominate choice-of-device in
##  refractory atony, and that ionised calcium during obstetric MTP should
##  correlate with subsequent uterine tone, not only with clot firmness.
## =====================================================================
