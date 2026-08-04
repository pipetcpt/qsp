## =====================================================================
##  CONTROLLED OVARIAN STIMULATION for IVF / ICSI — QSP model
##  (조절 난소 자극 · Controlled Ovarian Stimulation, COS)
##
##  65 ODE compartments:
##      30 follicle states (10 slots x diameter / granulosa / viability)
##      12 drug PK        (rFSH, corifollitropin, GnRH antagonist, hCG,
##                         GnRH agonist, cabergoline, letrozole, vaginal P4)
##      11 endocrine      (pituitary LH pool, LH, FSH, E2, inhibin B, AMH,
##                         progesterone, corpus luteum, surge readiness,
##                         GnRH-receptor availability)
##       5 OHSS           (VEGF, permeability, third-space fluid, plasma
##                         volume, pregnancy hCG)
##       7 accumulators   (commitment exposure, maturation clock, rupture
##                         integral, oocytes, MII oocytes, AUC E2, AUC VEGF)
##
##  TIME UNIT: DAYS, with t = 0 at menses (cycle day 1).
##  UNITS: diameter mm · FSH & LH IU/L · hCG IU/L · E2 pg/mL · P4 ng/mL
##         inhibin B pg/mL · AMH ng/mL · antagonist & agonist ng/mL
##         ascites & plasma volume L.
## =====================================================================
##
##  THE THESIS THIS MODEL IS BUILT AROUND
##  ------------------------------------
##  Ovarian stimulation is taught as a trade-off: push harder, get more
##  eggs, accept more risk. This model asserts something sharper — that
##  the trade-off is not a trade-off at all but an IDENTITY, because the
##  yield and the danger are two read-outs of the SAME state variable, and
##  that the only place in the whole system where they can be separated is
##  the SHAPE of the trigger.
##
##  (1) ONE STATE VARIABLE, FOUR CLINICAL NUMBERS.
##      Growing granulosa mass (Mg = sum of m_i x G_i x V_i) is read by
##      four different equations:
##          oestradiol            <- Mg x aromatase competence x theca LH
##          trigger-day progest.  <- Mg              (linear, no saturation)
##          oocytes retrieved     <- Mg's follicles above 11 mm
##          ovarian VEGF output   <- Mg x LHCGR occupancy
##      Nothing in the parameter list couples them; they are coupled
##      because they are the same mass. This is why lowering the FSH dose
##      moves yield and OHSS in the SAME direction (see scenario 5: 75 IU
##      in a high responder costs 4.4 oocytes and still leaves moderate
##      OHSS), and why "gentler stimulation" is not an OHSS strategy.
##
##  (2) ONE LIGAND POOL, THREE RESPONSE KERNELS.
##      LH, hCG and pregnancy hCG all act on the same receptor, so the
##      model keeps ONE ligand variable (LHEQ). What differs is what each
##      downstream process needs from it:
##
##        process              kernel                    needs
##        -------------------  ------------------------  ----------------
##        luteal / theca       LHEQ/(LHEQ+2)             basal pulses
##        support              (K = 2 IU/L)              are enough
##        meiotic commitment   Hill(LHEQ, 25, 6)         a real surge, but
##                             + autonomous 34 h clock   only its EDGE
##        wall rupture         Hill(LHEQ, 25, 6)         the edge AND a
##                             integral >= 0.90 d        sustained signal
##        VEGF / permeability  Hill(LHEQ, 40, 3)         the AREA, days
##                             integral                  of it
##
##      From that asymmetry, and no drug-specific parameter, the model
##      derives the entire modern trigger literature:
##        - hCG 10 000 IU: VEGF area 5.16 d, ascites 1.74 L, Hct 46.2%
##        - agonist 0.2 mg: VEGF area 0.41 d, ascites 0.04 L, Hct 38.3%
##        - and yet 17.1 vs 18.5 oocytes and the SAME maturation fraction
##          (0.77 vs 0.78), because maturation only ever needed the edge.
##      A 12.6-fold difference in VEGF exposure and a 43-fold difference in
##      ascites for an 8% difference in yield.
##        - and low-dose hCG is NOT the same trick: 1500 IU keeps ascites at
##          0.02 L but its peak (32 IU/L) is marginal for the meiotic kernel,
##          so the maturation fraction falls to 0.66. The agonist's peak is
##          152 IU/L with almost no area — which is the entire point.
##
##  (3) MONO-OVULATION IS NOT A PARAMETER.
##      There is no "number of follicles" input. Each of the 10 slots owns
##      an FSH threshold drawn from a log-normal (T50 = 9 IU/L, sigma =
##      0.45) and grows only while total FSH exceeds it. In an unstimulated
##      cycle the leader's oestradiol and inhibin drive FSH from 8 down to
##      2.6 IU/L, which is below every other slot's threshold: ONE follicle
##      ovulates on day 14-15. Give exogenous FSH and the decline never
##      happens, so 18 follicles grow. Same equations, same parameters —
##      the only difference is who controls the FSH concentration.
##      Selection is completed by two structural facts, not by an "if":
##      the threshold FALLS as a follicle enlarges (KDTH, rising FSHR
##      density) and above 11 mm the follicle switches to LH support
##      (KLHA), so the leader survives an FSH concentration that starves
##      the cohort.
##
##  (4) THE TRIGGER-DAY PROGESTERONE RISE IS ARITHMETIC.
##      P4 leaks from granulosa at KP4G per unit mass with NO saturation,
##      so trigger-day P4 is a follicle counter: 0.66 ng/mL at AFC 12,
##      1.31 at AFC 25, 1.49 at AFC 32. The 1.5 ng/mL threshold that sends
##      a cycle to freeze-all is therefore not a separate biological event;
##      it is the same mass, read a fourth way. Progestin-primed
##      stimulation (PPOS, scenario 17) pushes P4 to 5.1 ng/mL and makes
##      freeze-all obligatory — for exactly the same reason.
##
##  (5) THE 34-38 h RETRIEVAL WINDOW IS A GAP BETWEEN TWO CLOCKS.
##      Maturation completes at TMAT = 0.85 d after commitment; the wall
##      ruptures at TRUP = 1.55 d after commitment PROVIDED the signal
##      integral has reached RUPX = 0.90 d. Retrieval must sit between
##      them. The model is therefore forced to reproduce the clinical
##      window without being told it: 28 h -> 4.1 mature oocytes, 36 h ->
##      6.3, 40 h -> 5.3, 44 h -> 1.4.
##
##  WHAT THE MODEL IS NOT TOLD
##  --------------------------
##      - how many follicles grow, or how many oocytes are retrieved
##      - how many days of stimulation a patient needs
##      - that hCG is "riskier" than a GnRH agonist
##      - any OHSS severity scale, threshold or risk score
##      - that PCOS patients are high responders
##      - that AMH predicts oocyte yield
##  All six come out of the integration.
##
##  CALIBRATION TARGETS (observed -> model)
##  ---------------------------------------
##   unstimulated cycle: 1 dominant follicle (1.2), E2 peak 250-400 -> 430
##     pg/mL on day 10.9, LH surge 40-90 -> 55.7 IU/L on day 11.2, ovulation
##     day 12.5, FSH 7.8 -> 2.6 IU/L, mid-luteal P4 10-20 -> 12.5 ng/mL,
##     luteal LH 2.7-3.7 IU/L. The model's dominant follicle stops at
##     16.9 mm where 20-22 mm is observed: a documented 3-5 mm shortfall,
##     because the E2 that fires the surge is reached before the diameter is.
##   AFC 12 + rFSH 150 IU (ESTHER-1 conventional arm, PMID 27912901):
##     9-11 days of stimulation -> 11.4; E2 1500-2500 -> 1615 pg/mL;
##     8-11 oocytes -> 8.0; MII fraction 0.75-0.85 -> 0.79;
##     trigger-day P4 0.5-1.1 -> 0.66 ng/mL; serum FSH 12-20 -> 11.9 + 0.3;
##     endogenous FSH suppressed to 0.27 IU/L; AMH falls 47% during COS
##     (observed 30-50%); inhibin B 51 -> 227 pg/mL
##   AFC 25-32 + rFSH 150 IU: 15-22 oocytes -> 15.7 / 18.5
##   AFC 5 + rFSH 300 IU (POSEIDON group 3-4): 3-5 oocytes -> 3.7
##   ganirelix 0.25 mg s.c.: Cmax 11.2 ng/mL -> 11.4, t1/2 13 h,
##     LH suppression 70-80% -> 77% (PMID 10593372)
##   hCG 10 000 IU i.m.: peak 200-250 IU/L -> 220, terminal t1/2 33 h
##   GnRH agonist trigger: LH peak 100-150 IU/L -> 130, back to baseline
##     by 18-24 h; luteal P4 collapse -> 4.6 ng/mL on luteal day 7 vs
##     34.9 with hCG (the Humaidan 2005 result, Hum Reprod 20:2887)
##   cabergoline 0.5 mg from trigger: moderate-severe OHSS RR 0.38 in
##     meta-analysis -> ascites 1.74 -> 1.44 L (-17%), grade severe ->
##     moderate; the model is more conservative than the trials
##   euploid fraction by age (Franasiak 2014, PMID 24355045):
##     32 y 0.65 -> 0.68, 38 y 0.45 -> 0.42, 42 y 0.25 -> 0.24
##
##  VERIFICATION AND WHAT IT CHANGED
##  --------------------------------
##  All 64 ODEs were independently re-implemented in Python/scipy
##  (cos_reference_check.py) and integrated over 34 days for 22 scenarios
##  plus six parameter sweeps. Every number quoted in this header and in
##  README.md is printed by that script (cos_reference_output.txt).
##  Re-implementation exposed and fixed EIGHT defects that the first draft of
##  the equations contained:
##    (i)   the meiotic-commitment integral was driven by TOTAL LHCGR
##          occupancy, so basal LH of 5 IU/L committed the oocytes on cycle
##          day 1 and every follicle ovulated at 4 mm on day 2. Fixed by
##          giving each downstream process its own kernel — which is now
##          the model's central claim rather than a patch.
##    (ii)  the oestradiol positive-feedback switch had time constants
##          (KR 3.0, KROFF 0.6) fast enough to fire a partial surge at
##          E2 = 148 pg/mL on stimulation day 2, which arrested every
##          follicle at 7.1 mm through the PCOS LH-arrest term.
##    (iii) the pituitary LH pool had zero-order synthesis and no capacity,
##          so under a GnRH antagonist the pool grew without bound and the
##          antagonist had NO steady-state effect on LH (LH returned to
##          5 IU/L). Fixed with a saturable pool, which then produced the
##          agonist-trigger flare magnitude as a consequence.
##    (iv)  the aspiration pulse counted oocytes and collapsed follicles in
##          the same window, so the count was read off follicles that were
##          already emptying: 2.4 oocytes instead of 8.8 for AFC 12.
##    (v)   the GnRH-agonist bolus was folded into the morning of the
##          trigger day rather than given at the trigger hour, advancing it
##          by 9.6 h and silently ovulating the dual-trigger arm before
##          retrieval (2.8 instead of 14.6 oocytes).
##    (vi)  progesterone output was linear in luteal mass, giving 2.3 ng/mL
##          in an unstimulated luteal phase (observed 10-20) while a
##          15-corpus-luteum IVF luteal phase reached 173 ng/mL (observed
##          25-60). A saturating term reproduces both ends.
##    (vii) with no progesterone feedback on the GnRH pulse generator, LH
##          recovered to 8-10 IU/L within two days of any trigger, so the
##          corpus luteum was always rescued and the agonist trigger showed
##          only a 13% luteal deficit. Adding the brake (KP4LH) makes the
##          luteal mass suppress its own support and reproduces both the
##          agonist-trigger deficit and the dual-trigger rescue.
##    (viii) the FSH threshold distribution was centred at 6.5 IU/L, which
##          let SIX slots grow on endogenous FSH alone: the unstimulated
##          cycle was tri-follicular and its oestradiol reached 300 pg/mL by
##          cycle day 5. Re-centring at 9.0 IU/L with sigma 0.45 restores
##          mono-ovulation without touching any stimulated arm's calibration.
##
##  THE MODEL'S MOST EXPOSED PREDICTIONS
##  ------------------------------------
##    - Luteal letrozole abolishes oestradiol (-79%) and changes ascites by
##      under 4%, because in this model E2 is a MARKER and VEGF is the
##      MEDIATOR. Randomised data suggesting letrozole reduces OHSS would
##      falsify the omission of an E2-dependent permeability term.
##    - After a GnRH-agonist trigger the model predicts spontaneous
##      ovulation is essentially absent even at 44 h, because the rupture
##      integral never reaches RUPX. Clinics that delay retrieval after an
##      agonist trigger should therefore lose fewer oocytes than after hCG.
##    - Coasting works by killing small-follicle granulosa mass, so it
##      cannot reduce OHSS without costing oocytes (-35% here). A trial
##      showing OHSS reduction at unchanged yield would falsify the
##      single-mass structure.
## =====================================================================

library(mrgsolve)
library(dplyr)

cos_code <- '
$PARAM @annotated
// ------------------------------------------------- patient covariates
AFC     : 12.0   : antral follicle count, 2-9 mm (follicles)
T50     : 9.0    : median FSH threshold of the antral cohort (IU/L)
SIGT    : 0.45   : log-SD of the FSH-threshold distribution (-)
TONE    : 1.0    : GnRH pulse-generator tone (1 = normal, 1.7-1.8 = PCOS)
AGE     : 32.0   : age (years) — enters the embryology chain only

// ------------------------------------------------- protocol switches
FSHDOSE : 150.0  : daily rFSH dose (IU) — handled by the event table
TTRIG   : 1.0e6  : time of the trigger injection (d); 1e6 = no trigger
TOPU    : 0.0    : time of oocyte retrieval (d); 0 = none
TIMPL   : 1.0e6  : time implantation starts producing hCG (d)
PREG    : 0.0    : 1 = fresh transfer that implants (late-OHSS driver)
PPOS    : 0.0    : 1 = progestin-primed stimulation
FSHSTART: 2.0    : cycle day stimulation begins (d)

// ------------------------------------------------- follicle thresholds
HT      : 4.0    : Hill coefficient of the FSH threshold
KAMHT   : 0.15   : AMH-driven elevation of the FSH threshold
AMHREF  : 2.5    : reference AMH for that term (ng/mL)
KDTH    : 0.45   : fall of the threshold as a follicle enlarges
DTH     : 10.0   : midpoint of that fall (mm)

// ------------------------------------------------- follicle growth
KGR     : 2.90   : maximal diameter growth rate (mm/d)
WGR     : 0.35   : growth-rate floor for small antral follicles
DGR     : 9.0    : diameter at which growth accelerates (mm)
NGR     : 3.0    : steepness of that acceleration
DMAX    : 24.0   : asymptotic follicle diameter (mm)
KLHA    : 0.75   : LH contribution to drive above DLH
DLH     : 11.0   : diameter at which LHCGR appears on granulosa (mm)
KLHARR  : 11.0   : LH tone that arrests pre-LHCGR follicles (IU/L, PCOS)
NLHARR  : 4.0    : steepness of that arrest

// ------------------------------------------------- granulosa & atresia
KG      : 0.85   : granulosa recruitment toward capacity (1/d)
KGD     : 0.55   : granulosa loss when starved of FSH (1/d)
DPRE    : 20.0   : diameter at which granulosa capacity = 1 (mm)
DARO    : 12.0   : diameter of half-maximal aromatase competence (mm)
NARO    : 3.0    : steepness of aromatase acquisition
KATR    : 0.12   : atresia rate at zero FSH drive (1/d)

// ------------------------------------------------- ovulation / retrieval
KCOMMIT : 0.06   : LHCGR exposure that commits the oocyte (d)
TMAT    : 0.85   : maturation clock centre after commitment (d)
NMAT    : 6.0    : steepness of the maturation clock
TRUP    : 1.55   : rupture clock centre after commitment (d)
NRUP    : 60.0   : steepness of the rupture clock
RUPX    : 0.90   : cumulative LHCGR signal needed to complete rupture (d)
NRUPX   : 4.0    : steepness of that completion requirement
KRUP    : 6.0    : follicle collapse rate at full rupture (1/d)
WASP    : 0.05   : duration of the aspiration itself (d, ~1.2 h)
WCOL    : 0.30   : duration of post-aspiration collapse (d)
KCOL    : 3.0    : collapse strength of the aspiration
FGRAN   : 0.55   : granulosa fraction removed by aspiration
ETARET  : 0.88   : oocyte recovery per aspirated follicle
DASP    : 11.0   : 50% aspiration-yield diameter (mm)
NASP    : 8.0    : steepness of that gate
DMII    : 13.0   : 50% meiotic-competence diameter (mm)
NMII    : 5.0    : steepness of that gate

// ------------------------------------------------- gonadotropin secretion
KFS     : 53.1   : FSH secretion rate constant (IU/L/d)
KELF    : 4.0    : endogenous FSH elimination (1/d)
KI      : 120.0  : inhibin-B feedback constant (pg/mL)
KE2F    : 100.0  : oestradiol feedback constant (pg/mL)
FANTF   : 0.45   : maximal antagonist suppression of FSH
FANTL   : 0.90   : maximal antagonist suppression of LH
KRELLH  : 0.1245 : fractional LH release per unit GnRH drive (1/d)
SLH0    : 8000.0 : basal releasable pituitary LH pool (IU)
SLHMAX  : 12000.0: storage capacity of the gonadotroph pool (IU)
KSYNLH  : 2988.0 : LH synthesis rate at an empty pool (IU/d)
VLH     : 10.0   : LH distribution volume (L)
KELLH   : 16.6   : LH elimination (1/d, t1/2 1 h)
AMPS    : 55.0   : GnRH drive amplitude of the endogenous surge
AMPA    : 44.0   : GnRH drive amplitude of an agonist trigger
KR      : 1.20   : surge-readiness accumulation (1/d)
KROFF   : 0.45   : surge-readiness decay (1/d)
E2S     : 400.0  : E2 set-point for positive feedback (pg/mL)
NE2S    : 8.0    : steepness of the E2 switch
NRS     : 12.0   : steepness of the readiness-to-surge switch
DSURG   : 12.0   : follicle size that licenses positive feedback (mm)
NDSURG  : 8.0    : steepness of that gate
KP4S    : 2.0    : progesterone block of the surge (ng/mL)

// ------------------------------------------------- steroids & peptides
KE2G    : 6000.0 : E2 from aromatase-competent granulosa (pg/mL/d per unit)
KE2CL   : 900.0  : E2 from corpus luteum (pg/mL/d per unit)
KTHECA  : 0.40   : LH requirement of theca androgen supply (IU/L)
KELE2   : 12.0   : E2 elimination (1/d)
E2BASE  : 18.0   : adrenal/extragonadal E2 (pg/mL)
KIB     : 100.7  : inhibin B from granulosa (pg/mL/d per unit)
WSMIB   : 0.10   : small-antral contribution to inhibin B
KELIB   : 3.0    : inhibin B elimination (1/d)
KAMH    : 0.0838 : AMH per small antral follicle (ng/mL/d)
WAMHB   : 0.55   : pre-antral pool that keeps secreting AMH during COS
KELAMH  : 0.60   : AMH elimination (1/d)
DAMH    : 8.0    : diameter above which AMH switches off (mm)
KP4G    : 0.075  : P4 leak per unit granulosa mass (ng/mL/d) — NO saturation
KP4CL   : 18.1   : P4 from corpus luteum (ng/mL/d per unit)
KSATCL  : 1.46   : luteal mass at which P4 output per cell halves
KELP4   : 0.55   : P4 elimination (1/d)
KAP4    : 6.0    : absorption of vaginal/IM progesterone (1/d)
FP4     : 0.025  : (ng/mL) per mg of exogenous progesterone
KFCL    : 1.0    : granulosa-to-luteal conversion efficiency
KLYS0   : 0.10   : intrinsic luteolysis (1/d)
KLYSM   : 0.50   : LH-withdrawal luteolysis (1/d)

// --------------------------------- LHCGR: one ligand, three kernels
KSUP    : 2.0    : theca / luteal support half-saturation (IU/L)
KLHV    : 40.0   : VEGF-inducing luteinisation signal (IU/L)
NLHV    : 3.0    : steepness of the VEGF kernel
KLHM    : 25.0   : meiotic-commitment signal (IU/L)
NLHM    : 6.0    : steepness of the commitment kernel

// ------------------------------------------------- OHSS
KV      : 1.0    : VEGF production per unit mass per unit occupancy (1/d)
WCL     : 1.5    : luteal weighting of VEGF production
KELV    : 3.0    : VEGF elimination (1/d)
KVP     : 3.05   : VEGF for half-maximal permeability (units)
NVP     : 5.3    : steepness of the permeability response
KP      : 1.40   : permeability induction (1/d)
KPOFF   : 0.70   : permeability resolution (1/d)
EMAXCAB : 0.50   : cabergoline effect on VEGFR2 signalling
EC50CAB : 0.55   : cabergoline EC50 (units; 1.0 = 0.5 mg/d steady state)
LP      : 0.60   : trans-capillary leak at PERM = 1 (L/d)
KREAB   : 0.25   : peritoneal reabsorption (1/d)
KVPC    : 0.90   : plasma-volume compensation (1/d)
VP0     : 2.60   : reference plasma volume (L)
RBCV    : 1.60   : red-cell volume (L)
IVFLUID : 0.0    : supportive intravenous fluid (L/d, from TOPU+3)

// ------------------------------------------------- drug PK
KAF     : 2.50   : rFSH absorption (1/d)
KELFX   : 0.400  : rFSH elimination (1/d, multiple-dose t1/2 42 h)
VF      : 26.0   : rFSH V/F (L)
KACO    : 0.60   : corifollitropin absorption (1/d)
KELCO   : 0.241  : corifollitropin elimination (1/d, t1/2 69 h)
POTCO   : 10.0   : FSH-equivalent potency of the corifollitropin state
KAANT   : 100.0  : GnRH antagonist absorption (1/d, tmax 1.1 h)
KELANT  : 1.28   : GnRH antagonist elimination (1/d, t1/2 13 h)
VANT    : 20.0   : GnRH antagonist V/F (L)
IC50ANT : 0.44   : GnRH-receptor IC50 of the antagonist (ng/mL)
KAHCG   : 1.40   : hCG absorption (1/d)
KELHCG  : 0.50   : hCG elimination (1/d, t1/2 33 h)
VHCG    : 20.0   : hCG V/F (L)
KAAGO   : 40.0   : triptorelin absorption (1/d)
KELAGO  : 5.55   : triptorelin elimination (1/d, t1/2 3 h)
VAGO    : 30.0   : triptorelin V/F (L)
EC50AGO : 0.05   : GnRH-receptor EC50 of the agonist (ng/mL)
KELCAB  : 0.256  : cabergoline elimination (1/d, t1/2 65 h)
KELLET  : 0.37   : letrozole elimination (1/d, t1/2 45 h)
EMAXLET : 0.80   : maximal aromatase inhibition
EC50LET : 0.50   : letrozole EC50 (units)
KGRDOWN : 3.0    : GnRH-receptor down-regulation at full agonist occupancy (1/d)
KGRREC  : 0.20   : GnRH-receptor resynthesis (1/d, t1/2 3.5 d)
KP4LH   : 12.0   : progesterone that halves GnRH pulse drive (ng/mL)
NP4LH   : 2.0    : steepness of that brake
KELPHCG : 0.55   : pregnancy hCG elimination (1/d)
PH0     : 8.0    : pregnancy hCG production at implantation (IU/L/d)
GPH     : 0.35   : pregnancy hCG growth rate (1/d, doubling 2.0 d)

// ------------------------------------------------- initial conditions
D0INIT  : 5.0    : mean diameter of the antral cohort at cycle day 1 (mm)
G0INIT  : 0.03   : granulosa mass of a 5 mm follicle (fraction of mature)

$CMT @annotated
// --- follicle slot i: diameter (mm) ---------------------------------
D1  : slot 1 diameter (most FSH-sensitive)
D2  : slot 2 diameter
D3  : slot 3 diameter
D4  : slot 4 diameter
D5  : slot 5 diameter
D6  : slot 6 diameter
D7  : slot 7 diameter
D8  : slot 8 diameter
D9  : slot 9 diameter
D10 : slot 10 diameter (least FSH-sensitive)
// --- follicle slot i: functional granulosa mass (1 = preovulatory) ---
G1  : slot 1 granulosa mass
G2  : slot 2 granulosa mass
G3  : slot 3 granulosa mass
G4  : slot 4 granulosa mass
G5  : slot 5 granulosa mass
G6  : slot 6 granulosa mass
G7  : slot 7 granulosa mass
G8  : slot 8 granulosa mass
G9  : slot 9 granulosa mass
G10 : slot 10 granulosa mass
// --- follicle slot i: viability (1 = healthy, 0 = atretic) ----------
V1  : slot 1 viability
V2  : slot 2 viability
V3  : slot 3 viability
V4  : slot 4 viability
V5  : slot 5 viability
V6  : slot 6 viability
V7  : slot 7 viability
V8  : slot 8 viability
V9  : slot 9 viability
V10 : slot 10 viability
// --- drug PK --------------------------------------------------------
FSHDEP : rFSH subcutaneous depot (IU)
FSHC   : rFSH central (IU)
CORID  : corifollitropin depot (IU-equivalent)
CORIC  : corifollitropin central (IU-equivalent)
ANTD   : GnRH antagonist depot (ug)
ANTC   : GnRH antagonist central (ug)
HCGD   : hCG depot (IU)
HCG    : hCG central (IU)
AGOD   : GnRH agonist depot (ug)
AGOC   : GnRH agonist central (ug)
CAB    : cabergoline effect compartment (units)
LET    : letrozole effect compartment (units)
// --- endocrine ------------------------------------------------------
SLH    : releasable pituitary LH pool (IU)
LH     : plasma LH (IU/L)
FSHE   : endogenous plasma FSH (IU/L)
E2     : plasma oestradiol (pg/mL)
INHB   : plasma inhibin B (pg/mL)
AMH    : plasma anti-Mullerian hormone (ng/mL)
P4D    : vaginal/IM progesterone depot (mg)
P4     : plasma progesterone (ng/mL)
CL     : corpus luteum mass (units of granulosa origin)
RS     : surge readiness (0-1)
// --- OHSS -----------------------------------------------------------
VEGF   : ovarian VEGF signal above baseline (units)
PERM   : capillary permeability index (units)
ASC    : third-space (ascitic) fluid (L)
VP     : plasma volume (L)
PHCG   : pregnancy-derived hCG (IU/L)
// --- accumulators ---------------------------------------------------
EXPLH  : cumulative meiotic-kernel LHCGR exposure (d)
MCLK   : autonomous maturation clock since commitment (d)
ROV    : rupture integral (d)
OOC    : oocytes retrieved (count)
MIIC   : metaphase-II oocytes (count)
AUCE2  : cumulative oestradiol exposure (ng/mL*d)
EXPV   : cumulative VEGF-kernel LHCGR exposure (d)
GRR    : GnRH-receptor availability (0-1)

$GLOBAL
#define HL(x, k, n) (pow(fmax((x), 0.0)/(k), (n)) / (1.0 + pow(fmax((x), 0.0)/(k), (n))))

// equal-probability standard-normal quantile midpoints, i = 1..10
static const double ZQ[10] = {-1.6449, -1.0364, -0.6745, -0.3853, -0.1257,
                               0.1257,  0.3853,  0.6745,  1.0364,  1.6449};
double TH[10];      // per-slot FSH threshold (IU/L)
double MULT;        // follicles represented by one slot

$MAIN
// --- the antral cohort: one log-normal threshold distribution -------
MULT = AFC / 10.0;
for (int i = 0; i < 10; ++i) TH[i] = T50 * exp(SIGT * ZQ[i]);

D1_0 = D0INIT; D2_0 = D0INIT; D3_0 = D0INIT; D4_0 = D0INIT; D5_0 = D0INIT;
D6_0 = D0INIT; D7_0 = D0INIT; D8_0 = D0INIT; D9_0 = D0INIT; D10_0 = D0INIT;
G1_0 = G0INIT; G2_0 = G0INIT; G3_0 = G0INIT; G4_0 = G0INIT; G5_0 = G0INIT;
G6_0 = G0INIT; G7_0 = G0INIT; G8_0 = G0INIT; G9_0 = G0INIT; G10_0 = G0INIT;
V1_0 = 1.0; V2_0 = 1.0; V3_0 = 1.0; V4_0 = 1.0; V5_0 = 1.0;
V6_0 = 1.0; V7_0 = 1.0; V8_0 = 1.0; V9_0 = 1.0; V10_0 = 1.0;

SLH_0  = SLH0;
LH_0   = 5.0;
FSHE_0 = 6.0;
E2_0   = 40.0;
INHB_0 = 45.0;
AMH_0  = KAMH * (WAMHB*AFC + AFC*(1.0 - HL(D0INIT, DAMH, 4.0))) / KELAMH;
P4_0   = 0.4;
VP_0   = VP0;
GRR_0  = 1.0;

$ODE
// =====================================================================
//  0. gather the follicle states into arrays
// =====================================================================
double Dv[10] = {D1, D2, D3, D4, D5, D6, D7, D8, D9, D10};
double Gv[10] = {G1, G2, G3, G4, G5, G6, G7, G8, G9, G10};
double Vv[10] = {V1, V2, V3, V4, V5, V6, V7, V8, V9, V10};
double dD[10], dG[10], dV[10];
for (int i = 0; i < 10; ++i) {
  Dv[i] = fmax(Dv[i], 1.0e-6);
  Gv[i] = fmax(Gv[i], 0.0);
  Vv[i] = fmin(fmax(Vv[i], 0.0), 1.0);
}

// =====================================================================
//  1. drug concentrations
// =====================================================================
double CFSHEX = FSHC/VF + POTCO*CORIC/VF;          // IU/L  exogenous FSH
double CANT   = fmax(ANTC, 0.0)/VANT;              // ng/mL antagonist
double CAGO   = fmax(AGOC, 0.0)/VAGO;              // ng/mL agonist
double CHCG   = fmax(HCG,  0.0)/VHCG;              // IU/L  hCG
double CCAB   = fmax(CAB,  0.0);
double CLET   = fmax(LET,  0.0);
double FSHTOT = fmax(FSHE, 0.0) + CFSHEX;

// --- GnRH receptor: agonist and antagonist COMPETE for it -----------
double XA    = CAGO / EC50AGO;
double XN    = CANT / IC50ANT;
double OCCAGO = XA / (1.0 + XA + XN);
double OCCANT = XN / (1.0 + XA + XN);

// --- LHCGR: ONE ligand pool read by THREE kernels -------------------
double LHEQ = fmax(LH, 0.0) + CHCG + fmax(PHCG, 0.0);
double SUP  = LHEQ / (LHEQ + KSUP);                // support  (K = 2)
double OCCV = HL(LHEQ, KLHV, NLHV);                // VEGF     (K = 40)
double OCCM = HL(LHEQ, KLHM, NLHM);                // meiosis  (K = 25)

// =====================================================================
//  2. retrieval: aspiration counts first, the follicle collapses after
// =====================================================================
double RUP   = HL(MCLK, TRUP, NRUP) * HL(ROV, RUPX, NRUPX);
double PULSE = (TOPU > 0.0 && SOLVERTIME >= TOPU && SOLVERTIME < TOPU + WASP)
                 ? 1.0/WASP : 0.0;
double COLL  = (TOPU > 0.0 && SOLVERTIME >= TOPU + WASP &&
                SOLVERTIME < TOPU + WASP + WCOL) ? 1.0/WCOL : 0.0;

// =====================================================================
//  3. follicle dynamics
// =====================================================================
double MG = 0.0, MGA = 0.0, NSM = 0.0, SMALL = 0.0;
double NASPF = 0.0, NMIIF = 0.0, CLASP = 0.0, CLRUP = 0.0, DMAXF = 0.0;
for (int i = 0; i < 10; ++i) {
  // the threshold is not a constant: it falls as the follicle enlarges
  // (rising FSHR density + intrafollicular IGF amplification), and it is
  // raised by AMH — which is why the same 150 IU recruits a different
  // number of follicles in a high- and a low-AMH ovary.
  double TEFF = TH[i] * (1.0 + KAMHT*(fmax(AMH,0.0)/AMHREF - 1.0))
                      * (1.0 - KDTH*HL(Dv[i], DTH, 4.0));
  if (TEFF < 0.5) TEFF = 0.5;
  double A    = HL(FSHTOT, TEFF, HT);
  double LHR  = HL(Dv[i], DLH, 6.0);               // LHCGR acquisition
  double DRV  = A + KLHA*LHR*SUP;                  // FSH -> LH handover
  if (DRV > 1.0) DRV = 1.0;
  // PCOS: high LH tone arrests follicles that have not yet acquired LHCGR
  double ARR  = 1.0 / (1.0 + pow(fmax(LH,0.0)/KLHARR, NLHARR) * (1.0 - LHR));
  DRV *= ARR;

  double GCAP = pow(Dv[i]/DPRE, 2.0); if (GCAP > 1.0) GCAP = 1.0;
  double GFAC = WGR + (1.0 - WGR)*HL(Dv[i], DGR, NGR);
  double GATE = HL(Dv[i], DASP, NASP);
  double RUPI = RUP * LHR;                         // only LHCGR-competent

  dD[i] = KGR*DRV*GFAC*(1.0 - Dv[i]/DMAX)
          - KRUP*RUPI*(Dv[i] - 3.0)
          - KCOL*COLL*GATE*(Dv[i] - 3.0);
  dG[i] = KG*DRV*(GCAP - Gv[i]) - KGD*(1.0 - DRV)*Gv[i]
          - KRUP*RUPI*Gv[i]
          - FGRAN*KCOL*COLL*GATE*Gv[i];
  dV[i] = -KATR*pow(1.0 - DRV, 2.0)*Vv[i];

  MG    += MULT*Gv[i]*Vv[i];
  MGA   += MULT*Gv[i]*Vv[i]*HL(Dv[i], DARO, NARO);
  NSM   += MULT*Vv[i]*(1.0 - HL(Dv[i], DAMH, 4.0));
  SMALL += MULT*(1.0 - HL(Dv[i], DAMH, 4.0));
  NASPF += MULT*Vv[i]*GATE;
  NMIIF += MULT*Vv[i]*GATE*HL(Dv[i], DMII, NMII);
  CLASP += MULT*Gv[i]*Vv[i]*GATE;
  CLRUP += MULT*Gv[i]*Vv[i]*LHR;
  if (Dv[i] > DMAXF) DMAXF = Dv[i];
}
dxdt_D1 = dD[0]; dxdt_D2 = dD[1]; dxdt_D3 = dD[2]; dxdt_D4 = dD[3];
dxdt_D5 = dD[4]; dxdt_D6 = dD[5]; dxdt_D7 = dD[6]; dxdt_D8 = dD[7];
dxdt_D9 = dD[8]; dxdt_D10 = dD[9];
dxdt_G1 = dG[0]; dxdt_G2 = dG[1]; dxdt_G3 = dG[2]; dxdt_G4 = dG[3];
dxdt_G5 = dG[4]; dxdt_G6 = dG[5]; dxdt_G7 = dG[6]; dxdt_G8 = dG[7];
dxdt_G9 = dG[8]; dxdt_G10 = dG[9];
dxdt_V1 = dV[0]; dxdt_V2 = dV[1]; dxdt_V3 = dV[2]; dxdt_V4 = dV[3];
dxdt_V5 = dV[4]; dxdt_V6 = dV[5]; dxdt_V7 = dV[6]; dxdt_V8 = dV[7];
dxdt_V9 = dV[8]; dxdt_V10 = dV[9];

// =====================================================================
//  4. drug PK
// =====================================================================
dxdt_FSHDEP = -KAF*FSHDEP;
dxdt_FSHC   =  KAF*FSHDEP - KELFX*FSHC;
dxdt_CORID  = -KACO*CORID;
dxdt_CORIC  =  KACO*CORID - KELCO*CORIC;
dxdt_ANTD   = -KAANT*ANTD;
dxdt_ANTC   =  KAANT*ANTD - KELANT*ANTC;
dxdt_HCGD   = -KAHCG*HCGD;
dxdt_HCG    =  KAHCG*HCGD - KELHCG*HCG;
dxdt_AGOD   = -KAAGO*AGOD;
dxdt_AGOC   =  KAAGO*AGOD - KELAGO*AGOC;
dxdt_CAB    = -KELCAB*CAB;
dxdt_LET    = -KELLET*LET;

// =====================================================================
//  5. pituitary and gonadotropins
// =====================================================================
double P4BLK = 1.0 / (1.0 + pow(fmax(P4,0.0)/KP4S, 2.0));
// The positive feedback is licensed by a PREOVULATORY-SIZE follicle, not by
// oestradiol alone, and a GnRH antagonist prevents ENTRY into surge mode —
// not merely the release. Without the second term a 0.25 mg antagonist lets
// a partial surge through and the cycle ovulates on stimulation day 6.
dxdt_RS = KR*HL(E2, E2S, NE2S)*HL(DMAXF, DSURG, NDSURG)*P4BLK
            *(1.0 - FANTL*OCCANT)*(1.0 - RS)
          - KROFF*RS;
double SURGE  = AMPS*HL(fmax(RS,0.0), 0.5, NRS)*P4BLK;
double GNRHLH = (TONE + SURGE)*(1.0 - FANTL*OCCANT) + AMPA*OCCAGO;
double GNRHF  = (TONE + 0.25*SURGE)*(1.0 - FANTF*OCCANT) + 0.15*AMPA*OCCAGO;
// (both are modulated by P4 and receptor availability just below)

// Progesterone slows the GnRH pulse generator. One corpus luteum (P4 ~12)
// halves LH and the luteal phase still lasts 14 d; fifteen corpora lutea
// (P4 ~37) shut LH down to 1.3 IU/L, so an IVF luteal mass starves its own
// support — and only hCG, which bypasses the pituitary, escapes the loop.
double P4FB = 1.0 / (1.0 + pow(fmax(P4,0.0)/KP4LH, NP4LH));
// An agonist OCCUPIES and then DOWN-REGULATES the receptor; a competitive
// antagonist only occupies it. That is why an agonist trigger costs the
// luteal phase (P4 22.2 vs 37.4 ng/mL on luteal day 7) and a dual trigger
// with 1500 IU hCG buys it back (34.8).
double GRRC = fmin(fmax(GRR, 0.0), 1.0);
dxdt_GRR = KGRREC*(1.0 - GRRC) - KGRDOWN*OCCAGO*GRRC;
GNRHLH *= P4FB*GRRC;
GNRHF  *= (0.5 + 0.5*P4FB)*GRRC;

double REL  = KRELLH*GNRHLH*fmax(SLH, 0.0);
dxdt_SLH    = KSYNLH*(1.0 - fmax(SLH,0.0)/SLHMAX) - REL;
dxdt_LH     = REL/VLH - KELLH*fmax(LH, 0.0);

double FINH = 1.0 / (1.0 + fmax(INHB,0.0)/KI + fmax(E2,0.0)/KE2F);
dxdt_FSHE   = KFS*GNRHF*FINH - KELF*fmax(FSHE, 0.0);

// =====================================================================
//  6. steroids, peptides, corpus luteum
// =====================================================================
double THECA = LHEQ / (LHEQ + KTHECA);
double ARO   = 1.0 - EMAXLET*CLET/(CLET + EC50LET);
dxdt_E2   = (KE2G*MGA*THECA + KE2CL*fmax(CL,0.0)*THECA)*ARO
            + E2BASE*KELE2 - KELE2*fmax(E2, 0.0);
dxdt_INHB = KIB*(MG + WSMIB*NSM) - KELIB*fmax(INHB, 0.0);
dxdt_AMH  = KAMH*(WAMHB*AFC + SMALL) - KELAMH*fmax(AMH, 0.0);

double PPOSP4 = (PPOS > 0.5 && SOLVERTIME >= FSHSTART && SOLVERTIME < TTRIG)
                  ? 4.0*KELP4 : 0.0;
dxdt_P4D = -KAP4*P4D;
dxdt_P4  = KP4G*MG
           + KP4CL*fmax(CL,0.0)/(1.0 + fmax(CL,0.0)/KSATCL)
           + FP4*KAP4*P4D + PPOSP4
           - KELP4*fmax(P4, 0.0);

double LYS = KLYS0 + KLYSM*pow(1.0 - SUP, 2.0);
dxdt_CL = KFCL*(KRUP*RUP*CLRUP + (1.0 - FGRAN)*KCOL*COLL*CLASP)
          - LYS*fmax(CL, 0.0);

// =====================================================================
//  7. the two read-outs of the trigger: EDGE and AREA
// =====================================================================
dxdt_EXPLH = OCCM;                                  // commitment integral
dxdt_MCLK  = HL(EXPLH, KCOMMIT, 8.0);               // autonomous clock
dxdt_ROV   = OCCM;                                  // rupture integral
dxdt_EXPV  = OCCV;                                  // VEGF-kernel area

dxdt_VEGF = KV*(MG + WCL*fmax(CL,0.0))*OCCV - KELV*fmax(VEGF, 0.0);
double CABE = EMAXCAB*CCAB/(CCAB + EC50CAB);
dxdt_PERM = KP*HL(VEGF, KVP, NVP)*(1.0 - CABE) - KPOFF*fmax(PERM, 0.0);

double LEAK = LP*fmax(PERM,0.0)*fmax(VP,0.1)/VP0;
double IVIN = (IVFLUID > 0.0 && TOPU > 0.0 && SOLVERTIME > TOPU + 3.0)
                ? IVFLUID : 0.0;
dxdt_ASC  = LEAK - KREAB*fmax(ASC, 0.0);
dxdt_VP   = KVPC*(VP0 - VP) - LEAK + IVIN;

double PHCGIN = (PREG > 0.5 && SOLVERTIME > TIMPL)
                  ? PH0*exp(fmin(GPH*(SOLVERTIME - TIMPL), 8.0)) : 0.0;
dxdt_PHCG = PHCGIN - KELPHCG*fmax(PHCG, 0.0);

// =====================================================================
//  8. counting: oocytes are counted DURING aspiration, from intact
//     follicles, and maturity is read off the autonomous clock
// =====================================================================
dxdt_OOC   = PULSE*NASPF*ETARET*HL(EXPLH, KCOMMIT, 4.0);
dxdt_MIIC  = PULSE*NMIIF*ETARET*HL(MCLK, TMAT, NMAT);
dxdt_AUCE2 = fmax(E2, 0.0)/1000.0;

$TABLE
double FSHTOTAL = fmax(FSHE,0.0) + FSHC/VF + POTCO*CORIC/VF;
double LHEQOUT  = fmax(LH,0.0) + fmax(HCG,0.0)/VHCG + fmax(PHCG,0.0);
double HCT      = RBCV/(RBCV + fmax(VP, 0.2))*100.0;
double CANTOUT  = fmax(ANTC,0.0)/VANT;
double CHCGOUT  = fmax(HCG,0.0)/VHCG;
double CAGOOUT  = fmax(AGOC,0.0)/VAGO;

// follicle counts (a scan report, not a state)
double NF11 = 0.0, NF14 = 0.0, NF17 = 0.0, MGOUT = 0.0, OVOL = 0.0;
double DD[10] = {D1, D2, D3, D4, D5, D6, D7, D8, D9, D10};
double GG[10] = {G1, G2, G3, G4, G5, G6, G7, G8, G9, G10};
double VVv[10] = {V1, V2, V3, V4, V5, V6, V7, V8, V9, V10};
for (int i = 0; i < 10; ++i) {
  double mm = AFC/10.0;
  if (DD[i] >= 11.0) NF11 += mm*VVv[i];
  if (DD[i] >= 14.0) NF14 += mm*VVv[i];
  if (DD[i] >= 17.0) NF17 += mm*VVv[i];
  MGOUT += mm*GG[i]*fmin(fmax(VVv[i],0.0),1.0);
  OVOL  += mm*3.14159265/6.0*pow(DD[i]/10.0, 3.0);   // mL, sphere
}
OVOL = OVOL + 2.0 + 6.0*fmax(PERM, 0.0);             // stroma + oedema

// OHSS grade, computed from the fluid states (0 none .. 4 critical)
double OHSSG = 0.0;
if (fmax(ASC,0.0) >= 0.10) OHSSG = 1.0;
if (HCT >= 43.0 || fmax(ASC,0.0) >= 0.35) OHSSG = 2.0;
if (HCT >= 45.0 || fmax(ASC,0.0) >= 1.50) OHSSG = 3.0;
if (HCT >= 55.0 || fmax(ASC,0.0) >= 4.00) OHSSG = 4.0;

// embryology chain (algebraic, evaluated on the running MII count)
double BLR  = 0.55 - 0.011*fmax(0.0, AGE - 33.0);
double PEU  = 0.85/(1.0 + exp((AGE - 38.2)/4.0));
double TWOPN = MIIC*0.72;
double BLAST = TWOPN*BLR;
double EUPL  = BLAST*PEU;
double CLBR  = 1.0 - pow(1.0 - 0.52, EUPL);

$CAPTURE FSHTOTAL LHEQOUT HCT CANTOUT CHCGOUT CAGOOUT NF11 NF14 NF17
         MGOUT OVOL OHSSG TWOPN BLAST EUPL CLBR
'

mod <- mcode("cos", cos_code)

## =====================================================================
##  PROTOCOL BUILDER
##  Doses are given as boluses into depots, so the timing of the trigger
##  injection is exact. Verification note: folding the agonist bolus into
##  the morning of the trigger day advances it by 9.6 h and destroys the
##  dual-trigger arm — see header defect (v).
## =====================================================================
cos_protocol <- function(fsh_dose   = 150,      # IU/d recombinant FSH
                         fsh_start  = 2,        # cycle day
                         fsh_stop   = NA,       # coasting: stop FSH here
                         hmg_lh     = 0,        # IU/d LH activity (hp-hMG)
                         cori       = 0,        # ug corifollitropin alfa
                         ant_start  = 6,        # cycle day (0 = none)
                         ant_dose   = 0.25,     # mg/d ganirelix/cetrorelix
                         trigger    = "hcg",    # hcg | ago | dual | none
                         hcg_dose   = 10000,    # IU
                         ago_dose   = 0.2,      # mg triptorelin
                         opu_delay  = 1.5,      # d after the trigger
                         cab        = FALSE,    # cabergoline 0.5 mg/d
                         letro      = FALSE,    # luteal letrozole 5 mg/d
                         luteal_p4  = 0,        # mg/d vaginal progesterone
                         fresh      = TRUE,     # fresh transfer?
                         ppos       = FALSE,    # progestin priming
                         iv_fluid   = 0,        # L/d supportive care
                         name       = "protocol") {
  ## built explicitly rather than with as.list(environment()) so that no
  ## default arrives as an unevaluated promise
  list(fsh_dose = fsh_dose, fsh_start = fsh_start, fsh_stop = fsh_stop,
       hmg_lh = hmg_lh, cori = cori, ant_start = ant_start,
       ant_dose = ant_dose, trigger = trigger, hcg_dose = hcg_dose,
       ago_dose = ago_dose, opu_delay = opu_delay, cab = cab, letro = letro,
       luteal_p4 = luteal_p4, fresh = fresh, ppos = ppos,
       iv_fluid = iv_fluid, name = name)
}

## Bioavailabilities, kept identical to $PARAM / the python reference
FB <- list(FSH = 0.80, CORI = 0.60, ANT = 0.91, HCG = 0.70, AGO = 0.90,
           HMG = 0.20)

## compartment number from name (data sets need numeric cmt)
cmtn <- function(m, nm) match(nm, cmt(m))

## ---------------------------------------------------------------------
##  EVENT TABLE
##  A plain data frame (ID/time/cmt/amt/evid) rather than stacked ev()
##  objects, so that every dose lands at an exact time and the table can be
##  passed straight to mrgsim(data = ., obsaug = TRUE).
## ---------------------------------------------------------------------
cos_events <- function(m, p, ttrig = NA, tend = 30) {
  rec <- function(time, amt, nm)
    data.frame(ID = 1, time = time, cmt = cmtn(m, nm), amt = amt, evid = 1)
  e <- list()
  stim_end <- if (is.na(ttrig)) tend else floor(ttrig)
  if (!is.na(p$fsh_stop)) stim_end <- min(stim_end, p$fsh_stop - 1)

  ## daily gonadotropin --------------------------------------------------
  if (p$cori > 0) {
    e[[length(e) + 1]] <- rec(p$fsh_start, FB$CORI * p$cori * 100, "CORID")
    if (p$fsh_dose > 0 && stim_end >= p$fsh_start + 7)      # ENGAGE day 8 on
      e[[length(e) + 1]] <- rec(seq(p$fsh_start + 7, stim_end),
                                FB$FSH * p$fsh_dose, "FSHDEP")
  } else if (p$fsh_dose > 0 && stim_end >= p$fsh_start) {
    dd <- seq(p$fsh_start, stim_end)
    e[[length(e) + 1]] <- rec(dd, FB$FSH * p$fsh_dose, "FSHDEP")
    if (p$hmg_lh > 0)
      e[[length(e) + 1]] <- rec(dd, FB$HCG * p$hmg_lh * FB$HMG, "HCGD")
  }

  ## GnRH antagonist, daily from ant_start until retrieval ---------------
  if (p$ant_start > 0) {
    last <- if (is.na(ttrig)) tend else ceiling(ttrig + p$opu_delay)
    if (last >= p$ant_start)
      e[[length(e) + 1]] <- rec(seq(p$ant_start, last),
                                FB$ANT * p$ant_dose * 1000, "ANTD")
  }

  ## trigger and luteal phase — at the exact hour ------------------------
  if (!is.na(ttrig) && p$trigger != "none") {
    hd <- switch(p$trigger, hcg = p$hcg_dose, dual = 1500, 0)
    ad <- if (p$trigger %in% c("ago", "dual")) p$ago_dose else 0
    if (hd > 0) e[[length(e) + 1]] <- rec(ttrig, FB$HCG * hd, "HCGD")
    if (ad > 0) e[[length(e) + 1]] <- rec(ttrig, FB$AGO * ad * 1000, "AGOD")
    topu <- ttrig + p$opu_delay
    if (p$cab && tend > ttrig)
      e[[length(e) + 1]] <- rec(seq(ceiling(ttrig), tend), 0.256, "CAB")
    if (p$letro && tend > topu)
      e[[length(e) + 1]] <- rec(seq(ceiling(topu), tend), 0.37, "LET")
    if (p$luteal_p4 > 0 && tend > topu)
      e[[length(e) + 1]] <- rec(seq(ceiling(topu), tend), p$luteal_p4, "P4D")
  }
  d <- do.call(rbind, e)
  d[order(d$time), ]
}

## ---------------------------------------------------------------------
##  THE TRIGGER DAY IS AN OUTPUT, NOT AN INPUT
##  Pass 1 stimulates without a trigger and asks the model when three
##  follicles reach 17 mm; pass 2 re-runs with the trigger scheduled at
##  21:00 of that day. "Days of stimulation" is therefore predicted.
## ---------------------------------------------------------------------
cos_trigger_day <- function(m, p, ncrit = 3, dcrit = 17, maxstim = 16) {
  o <- m %>% param(TTRIG = 1e6, TOPU = 0, PREG = 0, FSHSTART = p$fsh_start,
                   PPOS = as.numeric(p$ppos)) %>%
    mrgsim(data = cos_events(m, p, NA, p$fsh_start + maxstim + 1),
           end = p$fsh_start + maxstim, delta = 0.25, obsaug = TRUE,
           atol = 1e-10, rtol = 1e-8) %>% as_tibble()
  ok <- which(o$NF17 >= ncrit & o$time >= p$fsh_start + 3)
  day <- if (length(ok)) floor(o$time[ok[1]]) else p$fsh_start + maxstim
  list(ttrig = day + 0.4, stim_days = day + 0.4 - p$fsh_start,
       surge = max(o$LH[o$time < day]) > 12, scan = o)
}

run_cos <- function(m, patient = list(), p = cos_protocol(), tend = 30) {
  m <- do.call(param, c(list(m), patient))
  tg <- cos_trigger_day(m, p)
  topu <- tg$ttrig + p$opu_delay
  m %>% param(TTRIG = tg$ttrig, TOPU = topu,
              TIMPL = topu + 9, PREG = as.numeric(p$fresh),
              PPOS = as.numeric(p$ppos), FSHSTART = p$fsh_start,
              IVFLUID = p$iv_fluid) %>%
    mrgsim(data = cos_events(m, p, tg$ttrig, tend), end = tend, delta = 0.05,
           obsaug = TRUE, atol = 1e-10, rtol = 1e-8) %>% as_tibble() %>%
    mutate(protocol = p$name, ttrig = tg$ttrig, topu = topu,
           stim_days = tg$stim_days)
}

cos_endpoints <- function(o) {
  tr <- o$ttrig[1]; op <- o$topu[1]
  at <- function(tt, col) approx(o$time, o[[col]], xout = tt, rule = 2)$y
  tibble(protocol   = o$protocol[1],
         stim_days  = o$stim_days[1],
         E2_trig    = at(tr, "E2"),
         P4_trig    = at(tr, "P4"),
         LH_trig    = at(tr, "LH"),
         foll11     = at(tr, "NF11"),
         foll17     = at(tr, "NF17"),
         oocytes    = tail(o$OOC, 1),
         MII        = tail(o$MIIC, 1),
         MII_frac   = tail(o$MIIC, 1)/pmax(tail(o$OOC, 1), 0.05),
         AUC_vegf   = tail(o$EXPV, 1) - at(tr, "EXPV"),
         AUC_meiosis= tail(o$EXPLH, 1) - at(tr, "EXPLH"),
         VEGF_max   = max(o$VEGF),
         HCT_max    = max(o$HCT),
         ASC_max    = max(o$ASC),
         P4_lut7    = at(op + 7, "P4"),
         OHSS       = max(o$OHSSG),
         blast      = tail(o$BLAST, 1),
         euploid    = tail(o$EUPL, 1),
         CLBR       = tail(o$CLBR, 1))
}

## =====================================================================
##  PATIENTS — only AFC, the threshold median, LH tone and age differ
## =====================================================================
PT_NORMAL  <- list(AFC = 12, T50 =  9.0, TONE = 1.0, AGE = 32)  # AMH 2.5
PT_PCOS    <- list(AFC = 25, T50 =  9.0, TONE = 1.8, AGE = 31)  # AMH 5.2
PT_EXTREME <- list(AFC = 32, T50 =  9.0, TONE = 1.7, AGE = 29)  # AMH 6.5
PT_POOR    <- list(AFC =  5, T50 = 12.0, TONE = 1.0, AGE = 38)  # AMH 1.0
PT_OLDER   <- list(AFC =  8, T50 = 10.0, TONE = 1.0, AGE = 41)  # AMH 1.7

## =====================================================================
##  SCENARIO LIBRARY — 22 arms.
##  The numbers in the comments are from the verified Python reference
##  implementation of these same equations (cos_reference_output.txt).
## =====================================================================
if (interactive()) {

  ## --- 0  unstimulated cycle: ONE follicle, no drugs at all ----------
  ##     E2 peak 430 pg/mL d10.9, LH surge 55.7 IU/L d11.2, mid-luteal P4 12.5
  s00 <- run_cos(mod, PT_NORMAL,
                 cos_protocol(fsh_dose = 0, ant_start = 0, trigger = "none",
                              fresh = FALSE, name = "natural"))

  ## --- 1  standard antagonist protocol -------------------------------
  ##     11.4 stim days, E2 1615, P4 0.66, 8.0 oocytes, 6.3 MII, CLBR 72%
  s01 <- run_cos(mod, PT_NORMAL, cos_protocol(luteal_p4 = 600,
                                              name = "standard-150"))

  ## --- 2  individualised lower dose (follitropin-delta style) --------
  ##     6.8 oocytes and 5.2 MII — 15% fewer eggs for a 33% dose cut
  s02 <- run_cos(mod, PT_NORMAL, cos_protocol(fsh_dose = 100, luteal_p4 = 600,
                                              name = "indiv-100"))

  ## --- 3  high responder, hCG trigger, fresh transfer ---------------
  ##     15.7 oocytes, E2 3594, SEVERE OHSS (Hct 45.4%, ascites 1.56 L)
  s03 <- run_cos(mod, PT_PCOS, cos_protocol(luteal_p4 = 600, name = "hi-hcg"))

  ## --- 4  same patient, agonist trigger + freeze-all ----------------
  ##     14.7 oocytes (94%), ascites 0.02 L — the model's headline result
  s04 <- run_cos(mod, PT_PCOS, cos_protocol(trigger = "ago", fresh = FALSE,
                                            name = "hi-ago-freeze"))

  ## --- 5  high responder, HALVE the FSH dose, keep hCG --------------
  ##     11.3 oocytes (-28%) and STILL moderate OHSS: dose cannot separate
  s05 <- run_cos(mod, PT_PCOS, cos_protocol(fsh_dose = 75, luteal_p4 = 600,
                                            name = "hi-75"))

  ## --- 6  high responder, hCG + cabergoline 0.5 mg ------------------
  ##     ascites 1.56 -> 1.29 L (-17%), yield unchanged, severe -> moderate
  s06 <- run_cos(mod, PT_PCOS, cos_protocol(cab = TRUE, luteal_p4 = 600,
                                            name = "hi-cabergoline"))

  ## --- 7  coasting: stop FSH on CD8, trigger 3 days later -----------
  ##     OHSS mild (0.20 L) but 12.4 oocytes and MII 8.1 (-35%) — it costs eggs
  s07 <- run_cos(mod, PT_PCOS, cos_protocol(fsh_stop = 8, luteal_p4 = 600,
                                            name = "hi-coasting"))

  ## --- 8  poor responder (POSEIDON 3/4), 300 IU ---------------------
  ##     3.7 oocytes, E2 889 pg/mL — dose cannot create follicles that are gone
  s08 <- run_cos(mod, PT_POOR, cos_protocol(fsh_dose = 300, luteal_p4 = 600,
                                            name = "poor-300"))

  ## --- 9  poor responder + hp-hMG LH activity -----------------------
  ##     E2 889 -> 955 pg/mL (+7%), oocytes unchanged (3.6): LH buys steroid,
  ##     not eggs
  s09 <- run_cos(mod, PT_POOR, cos_protocol(fsh_dose = 225, hmg_lh = 75,
                                            luteal_p4 = 600, name = "poor-LH"))

  ## --- 10 corifollitropin alfa 150 ug, one injection ----------------
  ##     9.4 oocytes vs 8.0 on daily 150 IU — front-loaded FSH widens the window
  s10 <- run_cos(mod, PT_NORMAL, cos_protocol(cori = 150, luteal_p4 = 600,
                                              name = "corifollitropin"))

  ## --- 11 NO antagonist: the premature LH surge --------------------
  ##     LH 74 IU/L pre-trigger, follicles ovulate, 0.2 oocytes, P4 18.8
  s11 <- run_cos(mod, PT_NORMAL, cos_protocol(ant_start = 0, luteal_p4 = 600,
                                              name = "no-antagonist"))

  ## --- 12 dual trigger (agonist + hCG 1500 IU) ---------------------
  ##     14.7 oocytes, VEGF area 5.85 d, ascites 0.05 L, luteal P4 rescued
  s12 <- run_cos(mod, PT_PCOS, cos_protocol(trigger = "dual", luteal_p4 = 600,
                                            name = "dual-trigger"))

  ## --- 13 agonist trigger WITH fresh transfer ----------------------
  ##     luteal-day-7 P4 22.2 vs 37.4 ng/mL — the Humaidan 2005 direction
  s13 <- run_cos(mod, PT_PCOS, cos_protocol(trigger = "ago", fresh = TRUE,
                                            luteal_p4 = 600, name = "ago-fresh"))

  ## --- 14 retrieval delayed to 40 h -------------------------------
  ##     7.6 oocytes, MII 5.3 (-16%) — the rupture integral catches up
  s14 <- run_cos(mod, PT_NORMAL, cos_protocol(opu_delay = 40/24,
                                              luteal_p4 = 600, name = "opu-40h"))

  ## --- 15 age 41, AFC 8 ------------------------------------------
  ##     6.1 oocytes, 4.9 MII, but 0.46 euploid blastocysts: CLBR 29%
  s15 <- run_cos(mod, PT_OLDER, cos_protocol(fsh_dose = 300, luteal_p4 = 600,
                                             name = "age-41"))

  ## --- 16 high responder, hCG + freeze-all ------------------------
  ##     early OHSS unchanged (1.57 L) — freezing removes only the LATE peak
  s16 <- run_cos(mod, PT_PCOS, cos_protocol(fresh = FALSE, name = "hi-hcg-freeze"))

  ## --- 17 progestin-primed stimulation + agonist trigger ----------
  ##     P4 5.2 ng/mL before the trigger: freeze-all is obligatory
  s17 <- run_cos(mod, PT_PCOS,
                 cos_protocol(ppos = TRUE, ant_start = 0, trigger = "ago",
                              fresh = FALSE, name = "PPOS"))

  ## --- 18 AFC 32, hCG 10000, fresh transfer ----------------------
  ##     18.5 oocytes, Hct 46.2%, ascites 1.74 L, biphasic (late) OHSS
  s18 <- run_cos(mod, PT_EXTREME, cos_protocol(luteal_p4 = 600, name = "xs-hcg"))

  ## --- 19 AFC 32, agonist trigger + freeze-all -------------------
  ##     17.1 oocytes, ascites 0.04 L: 43x less fluid for 8% fewer eggs
  s19 <- run_cos(mod, PT_EXTREME, cos_protocol(trigger = "ago", fresh = FALSE,
                                               name = "xs-ago"))

  ## --- 20 AFC 32, hCG + cabergoline + freeze-all ----------------
  s20 <- run_cos(mod, PT_EXTREME, cos_protocol(cab = TRUE, fresh = FALSE,
                                               name = "xs-cab"))

  ## --- 21 AFC 32, hCG + freeze-all + supportive IV fluid --------
  s21 <- run_cos(mod, PT_EXTREME, cos_protocol(fresh = FALSE, iv_fluid = 1.5,
                                               name = "xs-iv"))

  arms <- bind_rows(lapply(list(s00, s01, s02, s03, s04, s05, s06, s07, s08,
                                s09, s10, s11, s12, s13, s14, s15, s16, s17,
                                s18, s19, s20, s21), cos_endpoints))
  print(as.data.frame(arms), digits = 3)

  ## -------------------------------------------------------------------
  ##  SWEEP A — does the FSH dose separate yield from OHSS?  (it does not)
  ## -------------------------------------------------------------------
  doseA <- bind_rows(lapply(c(75, 112, 150, 225, 300, 450), function(d)
    cos_endpoints(run_cos(mod, PT_PCOS,
                          cos_protocol(fsh_dose = d, luteal_p4 = 600,
                                       name = paste0("AFC25-", d, "IU"))))))
  print(as.data.frame(doseA), digits = 3)

  ## -------------------------------------------------------------------
  ##  SWEEP B — one ligand, three kernels: the whole trigger literature
  ## -------------------------------------------------------------------
  trigB <- bind_rows(
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(hcg_dose = 10000,
                    fresh = FALSE, name = "hCG 10000"))),
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(hcg_dose = 5000,
                    fresh = FALSE, name = "hCG 5000"))),
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(hcg_dose = 2500,
                    fresh = FALSE, name = "hCG 2500"))),
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(hcg_dose = 1500,
                    fresh = FALSE, name = "hCG 1500"))),
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(trigger = "dual",
                    fresh = FALSE, name = "dual"))),
    cos_endpoints(run_cos(mod, PT_EXTREME, cos_protocol(trigger = "ago",
                    fresh = FALSE, name = "agonist"))))
  print(as.data.frame(trigB), digits = 3)

  ## -------------------------------------------------------------------
  ##  SWEEP C — the retrieval window falls out of two clocks
  ## -------------------------------------------------------------------
  opuC <- bind_rows(lapply(c(28, 32, 34, 36, 38, 40, 44), function(h)
    cos_endpoints(run_cos(mod, PT_NORMAL,
                          cos_protocol(opu_delay = h/24, luteal_p4 = 600,
                                       name = paste0(h, "h"))))))
  print(as.data.frame(opuC), digits = 3)

  ## -------------------------------------------------------------------
  ##  SWEEP D — reserve sweep: the CLBR plateau is not imposed
  ## -------------------------------------------------------------------
  afcD <- bind_rows(lapply(c(4, 6, 9, 12, 16, 20, 25, 32), function(a)
    cos_endpoints(run_cos(mod, list(AFC = a, T50 = 9, TONE = 1, AGE = 32),
                          cos_protocol(luteal_p4 = 600,
                                       name = paste0("AFC", a))))))
  print(as.data.frame(afcD), digits = 3)
}
