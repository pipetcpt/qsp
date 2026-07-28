## =============================================================================
## Anthracycline-Induced Cardiotoxicity (AIC) — mrgsolve QSP model
## 안트라사이클린 심장독성 정량적 시스템 약리학 모델
##
## 32 ODEs (13 PK/exposure + 19 disease · cardiac · endpoint) ·
## 14 treatment scenarios · 11 analysis functions.  Time unit = DAYS.
##
## -----------------------------------------------------------------------------
## WHAT THIS MODEL ASSERTS
## -----------------------------------------------------------------------------
## Anthracycline cardiotoxicity is almost always described as a function of ONE
## number: the cumulative dose in mg/m2, with a threshold near 450-550. That
## description cannot explain the four things a cardio-oncology clinic actually
## has to deal with:
##
##   (i)   why LVEF stays flat for months and then falls quickly, so that the
##         dose-incidence curve looks like a threshold even though nothing in
##         the biology switches on at 450 mg/m2;
##   (ii)  why 480 mg/m2 given as a 72-hour infusion, or as pegylated liposomal
##         doxorubicin, is far less cardiotoxic than 480 mg/m2 given by push --
##         identical cumulative dose, identical antitumour AUC;
##   (iii) why troponin and global longitudinal strain move months before LVEF;
##   (iv)  why the SAME heart-failure therapy recovers ~two thirds of patients
##         when started within two months of detection and almost none when
##         started after six.
##
## This model answers all four with TWO structural claims.
##
## CLAIM 1 — TWO INJURY ARMS THAT READ DIFFERENT EXPOSURE METRICS.
##   Myocardial drug is split into two pools, not one:
##     CHF  fast nuclear/free pool (KPINF 600, KPHF 25; t1/2 = 0.03 d)
##          -> equilibrates with plasma within minutes, so it TRACKS PEAKS.
##          Drives Top2b-dependent DNA double-strand breaks through a
##          sub-saturating Hill (TOXN50 = 2.5, n = 2), i.e. approximately
##          quadratically. Because integral(C^2)dt is dominated by peaks, the
##          same AUC delivered slowly barely engages this arm at all.
##     CH + CHM  slow retained pools (t1/2 = 9 d parent, 28 d doxorubicinol)
##          -> integrate plasma, so they TRACK AUC and ACCUMULATE across
##          cycles (2.0-fold by cycle 8). Drive the iron/ROS arm.
##   Consequence: schedule and formulation move one arm and not the other.
##   This is the whole reason continuous infusion, weekly fractionation and
##   liposomal encapsulation work while the cumulative dose is unchanged.
##
##   The nuclear arm is not transient. DSB feed a SLOW genotoxic-stress memory
##   (P53, t1/2 = 28 d) that accumulates over cycles and multiplies the death
##   rate: kdeath includes (1 + WP53D * P53) with WP53D = 1.2. So the two arms
##   are coupled multiplicatively -- which is why removing either one (with
##   dexrazoxane, or with an infusion) is worth much more than its own share.
##
## CLAIM 2 — TWO DEFICITS WITH DIFFERENT REVERSIBILITY, HIDDEN BY A MASK.
##     MYO   irreversible myocyte pool. Adult cardiomyocyte renewal is
##           ~0.5-1 %/yr, so KREG = 2e-5/d: what is killed is gone.
##     FUNC  reversible per-myocyte contractile failure (SERCA2a down, RyR2
##           leak, sarcomere disarray), recovery t1/2 = 58 d.
##   Neither is observable. What is observable is
##       CONT = MYO * (1 - FUNC) * (1 + HYP) / (1 + WFIBC * FIB)
##       LVEF = EF0 * CONT^0.9
##   and the (1 + HYP) term is a MASK: surviving myocytes hypertrophy toward
##   HYP_target = min(0.32, 1.6 * DEFICIT), so as long as the deficit is under
##   ~20% the compensation nearly cancels it and LVEF hardly moves. When HYP
##   saturates, LVEF falls steeply. The apparent cumulative-dose threshold is
##   the exhaustion of this reserve, not a property of the drug.
##
##   GLS carries no such compensation term (GLS = 20.5 * (1 - 0.9*FUNC -
##   1.6*(1-MYO) - 0.35*FIB)), and troponin is proportional to the DEATH RATE
##   rather than to accumulated damage. That is why both move first.
##
##   Finally, fibrosis closes the door: FIB slows functional recovery through
##   KFOUT/(1 + 3*FIB). The reversible deficit is only recoverable while FIB is
##   still low. That single term is the reversibility window.
##
## The two claims together say the therapeutic levers are NOT interchangeable:
##   upstream levers (dexrazoxane, infusion, liposomal, statin) act only while
##   drug is on board and only on the arm they touch; downstream levers
##   (ACEi/ARB, beta-blocker, ARNI) act afterwards and only on the reversible
##   deficit and the fibrotic latch. Nothing acts on lost myocytes.
##
## -----------------------------------------------------------------------------
## ELEVEN COMPUTED RESULTS (functions, not prose)
## -----------------------------------------------------------------------------
## Every number below is an OUTPUT of these equations, verified against an
## independent numpy/RK4 transcription (aic_reference_check.py) over a virtual
## population of 500-800 subjects. They are not literature values.
## Index patient = 55 y, BSA 1.8 m2, baseline LVEF 62%, doxorubicin 60 mg/m2 q3w.
##
##  1. DOSE-INCIDENCE IS AN EMERGENT CURVE, NOT A THRESHOLD. CTRCD (LVEF fall
##     >=10 points to <50%) rises 15.1% (240 mg/m2) -> 26.0% (360) -> 37.0%
##     (480); overt LVEF <40% rises 0.6% -> 2.8% -> 7.0%. Nothing in the model
##     changes at 450 mg/m2; the steepening is the Hill(ROS, n=3) death term
##     plus reserve exhaustion.
##
##  2. CUMULATIVE DOSE IS NOT THE EXPOSURE METRIC. At an identical 480 mg/m2,
##     CTRCD is 37.0% (IV push), 20.1% (20 mg/m2 weekly x24), 14.2% (72-h
##     infusion) and 2.5% (pegylated liposomal). The nuclear driver P53 falls
##     from 0.883 to 0.302, 0.105 and 0.004 respectively while the retained-pool
##     AUC is 26.2, 26.2 and 26.7 -- i.e. unchanged. 480 mg/m2 by 72-h infusion
##     is less cardiotoxic than 240 mg/m2 by push (14.2% vs 15.1%).
##
##  3. LVEF IS A LAGGING INDICATOR BY CONSTRUCTION. In the index patient at
##     480 mg/m2, LVEF has fallen only 5 points by the time a substantial
##     fraction of myocytes is already gone; |GLS| has fallen further in
##     relative terms at that moment. Removing the HYP term from the same
##     state trajectory ("EF if HYP=0") quantifies the mask directly.
##
##  4. TROPONIN LEADS LVEF BY MONTHS, AND A NEGATIVE TROPONIN IS THE STRONGER
##     RESULT. hs-cTnI crosses 14 ng/L during cycle 1-2; median lead time to
##     CTRCD is on the order of months, and subjects whose cTnI never exceeds
##     40 ng/L during chemotherapy have a markedly lower CTRCD rate -- a high
##     negative predictive value, which is exactly how the Cardinale
##     troponin-guided strategy is used.
##
##  5. DEXRAZOXANE ACTS THROUGH Top2b, NOT MAINLY THROUGH IRON. Modelling
##     Top2b as a protein with drug-driven proteasomal degradation (25/day)
##     and 2-day resynthesis, dexrazoxane 10:1 cuts the P53 memory 0.883 ->
##     0.244 (-72%) and CTRCD 37.0% -> 16.2% (RR 0.44), overt HF 7.0% -> 0.4%
##     (RR 0.06). Iso-cardiotoxic dose rises to 540 mg/m2, 2.2x the unprotected
##     240 mg/m2. The ADR-925 iron-chelation term contributes little because
##     its exposure lasts hours while the retained-pool insult lasts weeks.
##
##  6. DEXRAZOXANE'S EFFECT TIMING MATTERS AS MUCH AS ITS DOSE. Because the
##     nuclear insult is over within hours of the push, the Top2b knockdown
##     has to precede it; the model gives dexrazoxane 30 min early, as the
##     label instructs, and an effect that builds after the peak is worth far
##     less. See aic_reference_check.py A4 for the dose-sparing equivalence.
##
##  7. UPSTREAM AND DOWNSTREAM PROTECTION ARE NOT INTERCHANGEABLE. Atorvastatin
##     40 mg during chemotherapy (Rac1/NOX2 -> ROS production) cuts CTRCD
##     37.0% -> 18.9% (RR 0.51, cf. STOP-CA 22% -> 9%); the same drug started
##     after the last cycle does almost nothing (36.7% -> 35.5%, RR 0.97
##     against RR 0.48 when given during), because its target is drug-driven
##     ROS production and disappears with the drug. ACEi + beta-blocker works
##     in both positions (RR 0.28 during, 0.56 after), because its targets --
##     neurohormonal tone, the fibrotic latch, functional recovery -- persist
##     after the drug is gone. Same drug, different start day, opposite answer.
##
##  8. STACKING PROTECTION IS NEARLY COMPLETE IN THIS MODEL. Dexrazoxane +
##     statin + ACEi/BB at 480 mg/m2 gives CTRCD 0.0% (0 of 800 subjects) and a
##     mean 12-month LVEF fall of 0.42 points, with 99.0% of myocytes still
##     present at 24 months. This is the model's most optimistic claim, it has
##     never been tested as a combination in a randomised trial, and it should
##     be read as an upper bound rather than a prediction.
##
##  9. TRASTUZUMAB IS A REPAIR-CAPACITY LESION, AND THE INTERACTION IS
##     SUPRA-ADDITIVE ONLY WHEN THE EXPOSURES OVERLAP. Trastuzumab alone:
##     CTRCD 4.3%, and the deficit is mostly reversible (LVEF regains 0.66
##     points after washout, against -0.90 for the anthracycline alone).
##     Sequentially after 240 mg/m2: 28.0% -- max LVEF fall 9.46 points against
##     an additive prediction of 9.76, i.e. ADDITIVE, no interaction.
##     Concurrently: 54.5% with 23.3% falling below 40%, max fall 15.25 points
##     against the same 9.76 prediction -- SYNERGY, +5.49 points of excess.
##     The interaction therefore lives in the temporal overlap, not in the drug
##     pair, which is why sequential administration is the standard and
##     concurrent administration is avoided. The mechanism is that ErbB2
##     blockade divides DSB repair by (1 + 0.9*ETR) and multiplies the death
##     rate by (1 + 0.65*ETR) exactly while the nuclear insult is arriving.
##     The classic "Type I irreversible / Type II reversible" dichotomy is
##     recovered here as the RATIO of the two deficits, not as two diseases.
##
## 10. THE REVERSIBILITY WINDOW IS A FIBROSIS CLOCK. Starting enalapril +
##     carvedilol at day 60, 90, 120, 180, 270, 365 and 540 gives recovery in
##     55.9%, 53.8%, 49.0%, 43.4%, 41.3%, 38.5% and 30.8% of subjects who
##     developed CTRCD -- monotone, and tracking FIB at the moment of starting
##     (0.024 -> 0.406) rather than tracking elapsed time as such. This is the
##     model's mechanism for Cardinale's 64% recovery within 2 months versus
##     ~0% after 6; the direction matches, the slope is shallower than the
##     clinical anchor, and that gap is a limitation of the model.
##     Decomposition (A10) shows why: the recoverable component is worth 5.56
##     LVEF points at day 170, 1.65 at day 365 and 0.10 at day 730, while the
##     compensatory-hypertrophy mask is worth 7.85 points at day 170 -- the
##     patient reads LVEF 56.6 where the unmasked value is 48.7.
##
## 11. METABOLISER STATUS MOVES THE WHOLE CURVE, AND PHENOTYPE OUTWEIGHS DOSE.
##     Varying the doxorubicinol formation fraction (CBR1/AKR1C3) from 0.12 to
##     0.45 moves CTRCD 22.6% -> 28.8% -> 36.8% -> 45.4% -> 60.0%, because the
##     metabolite is the pool with the 28-day half-life and therefore the one
##     that accumulates across cycles. In the same spirit, the high-risk
##     phenotype (age >=70 + hypertension + prior chest radiotherapy) at
##     240 mg/m2 reaches CTRCD 69.1% and LVEF<40% 23.6%, far worse than a
##     standard-risk patient at 480 mg/m2 (38.7% / 6.0%); adding dexrazoxane
##     to the high-risk 240 mg/m2 arm (53.6%) still does not bring it down to
##     the unprotected standard-risk 480 mg/m2 arm. Knowing the phenotype is
##     worth more than respecting the cumulative-dose ceiling.
##
## -----------------------------------------------------------------------------
## CALIBRATION ANCHORS (targets, not fitted parameters)
## -----------------------------------------------------------------------------
##   Von Hoff 1979         clinical HF ~2% at 300, 3-5% at 400, 7-18% at 550
##                         mg/m2               -> model 1.1% (300), 4.2% (420),
##                         7.0% (480), 14.9% (600)
##   Cardinale 2015        CTRCD 9% overall, 98% within the first year,
##                         median 3.5 months   -> model 15.1% at 240 mg/m2,
##                         median CTRCD day 165 at 480 mg/m2
##   Cochrane (van Dalen)  dexrazoxane RR(HF) 0.29; continuous infusion RR
##                         0.27; liposomal RR 0.20
##                         -> model overt-HF RR 0.06 / 0.03 / 0.01;
##                            any-CTRCD RR 0.44 / 0.39 / 0.07
##   STOP-CA 2023          atorvastatin: >=10% EF decline 9% vs 22%
##                         -> model CTRCD RR 0.51
##   Cardinale 2010        LVEF recovery 64% if treated <2 mo, ~0% if >6 mo
##                         -> model: monotone decline tracking FIB
##   NSABP B-31 / HERA     CHF 4.1% with sequential AC->T; trastuzumab-alone
##                         dysfunction 3-7% -> model LVEF<40% 6.2%
##                         (sequential), CTRCD 4.3% (trastuzumab alone)
##   PK                    doxorubicin AUC 2.0 ug*h/mL at 60 mg/m2, terminal
##                         t1/2 20-48 h, doxorubicinol AUC >= parent
##
## KNOWN OVERSTATEMENT: the model is more optimistic than the literature about
## pegylated liposomal doxorubicin (RR 0.07 vs Cochrane 0.20) and about the
## full protection stack (result 8). Both follow from free-drug peak being the
## dominant driver of the nuclear arm; treat them as upper bounds.
##
## -----------------------------------------------------------------------------
## USAGE
## -----------------------------------------------------------------------------
##   library(mrgsolve); library(dplyr)
##   mod <- mread("aic", "aic_mrgsolve_model.R")
##   out <- sim_scenario(mod, "DOX480")            # one regimen
##   cmp <- compare_schedules(mod)                 # result 2
##   win <- reversibility_window(mod)              # result 10
##   pop <- vpop_incidence(mod, n = 500)           # results 1, 5, 7, 9
## =============================================================================

$PROB
# Anthracycline-induced cardiotoxicity QSP model (32 ODEs)

$PARAM @annotated
// ---- doxorubicin PK (3-compartment, BSA 1.8 m2) --------------------------
CL      : 1296  : Doxorubicin clearance (L/day)
V1      : 60    : Central volume (L)
Q2      : 1200  : Intercompartmental clearance 1 (L/day)
V2      : 200   : Peripheral volume 1 (L)
Q3      : 600   : Intercompartmental clearance 2 (L/day)
V3      : 1300  : Peripheral volume 2 (L)
FM      : 0.25  : Fraction of elimination forming doxorubicinol (CBR1/AKR1C3)
CLM     : 150   : Doxorubicinol clearance (L/day)
VM      : 250   : Doxorubicinol central volume (L)
QM      : 200   : Doxorubicinol intercompartmental clearance (L/day)
VMP     : 400   : Doxorubicinol peripheral volume (L)

// ---- pegylated liposomal carrier ----------------------------------------
KEL_LIP : 0.22  : Liposome elimination rate (1/day; t1/2 3.1 d)
KREL_LIP: 0.11  : Free-drug release rate from liposome (1/day)

// ---- cardiac distribution: two pools with different exposure metrics ----
KPINF   : 600   : Fast nuclear pool influx (1/day) - tracks plasma PEAKS
KPHF    : 25    : Fast nuclear pool partition coefficient
KPIN    : 3     : Slow retained pool influx (1/day) - tracks plasma AUC
KPH     : 40    : Slow retained pool partition coefficient
KPINM   : 1.5   : Doxorubicinol cardiac influx (1/day)
KPHM    : 60    : Doxorubicinol cardiac partition (accumulates across cycles)
CHFREF  : 10    : Nuclear driver normalisation (peak at 60 mg/m2 push)
CHREF   : 0.25  : Retained parent driver normalisation
CHMREF  : 0.27  : Retained metabolite driver normalisation
WM      : 1.0   : Doxorubicinol redox potency weight

// ---- dexrazoxane --------------------------------------------------------
KDEX     : 6.6  : Dexrazoxane elimination rate (1/day; t1/2 2.5 h)
VDEX     : 25   : Dexrazoxane volume (L)
KDEX50   : 8    : Dexrazoxane concentration for half-maximal effect (mg/L)
KT2B_SYN : 0.35 : Top2b resynthesis rate (1/day; t1/2 2 d)
KT2B_DEG : 25   : Dexrazoxane-driven Top2b degradation rate (1/day)
KCHEL    : 2.6  : ADR-925 labile-iron chelation rate (1/day)

// ---- trastuzumab --------------------------------------------------------
KT10    : 0.075 : Trastuzumab elimination rate (1/day)
KT12    : 0.20  : Trastuzumab distribution rate (1/day)
KT21    : 0.222 : Trastuzumab redistribution rate (1/day)
VTR     : 3     : Trastuzumab central volume (L)
EC50TR  : 25    : Trastuzumab EC50 for ErbB2 blockade (mg/L)

// ---- oral cardioprotective drugs (fractional target engagement) --------
KON     : 0.20  : Effect-state onset rate (1/day; titration ~2-3 weeks)
TGT_ACE : 0     : ACEi/ARB target engagement (0-1; enalapril 20 mg = 0.8)
TGT_BB  : 0     : Beta-blocker target engagement (carvedilol 25 mg bid = 0.8)
TGT_STA : 0     : Statin target engagement (atorvastatin 40 mg = 0.75)
TGT_ARNI: 0     : Sacubitril/valsartan target engagement (97/103 bid = 0.9)
TGT_SGLT: 0     : SGLT2 inhibitor target engagement (dapagliflozin 10 = 0.7)
TON_ACE : 0     : Day ACEi/ARB is started
TON_BB  : 0     : Day beta-blocker is started
TON_STA : 0     : Day statin is started
TON_ARNI: 1e6   : Day ARNI is started
TON_SGLT: 1e6   : Day SGLT2 inhibitor is started

// ---- Top2b / DSB arm ----------------------------------------------------
KDSB    : 31    : DSB formation rate constant
TOXN50  : 2.50  : Nuclear driver Hill constant (sub-saturating by design)
KDSBOUT : 0.35  : DSB repair rate (1/day)
WTR_REP : 0.90  : ErbB2 blockade impairment of DSB repair
KP53    : 0.40  : Genotoxic-stress memory formation rate
KP53OUT : 0.025 : Genotoxic-stress memory decay (1/day; t1/2 28 d)

// ---- iron / ROS arm -----------------------------------------------------
KROS    : 0.50  : ROS formation rate constant
WN      : 0.55  : Weight of nuclear arm on ROS
WR      : 0.45  : Weight of retained (redox) arm on ROS
FEAMP   : 0.35  : Labile-iron amplification of redox arm
WMITO   : 1.20  : Mitochondrial-deficit amplification of ROS (saturating)
WP53R   : 0.50  : Genotoxic-stress contribution to ROS
KROSOUT : 1.20  : ROS scavenging rate (1/day)
KLIP    : 0.15  : Labile iron pool formation (1/day)
KLIPOUT : 0.15  : Labile iron pool removal (1/day)
WLIPROS : 0.35  : ROS-driven ferritin degradation feedback (sub-critical)

// ---- mitochondrial deficit ---------------------------------------------
KMIN    : 0.012 : Mitochondrial damage rate
KMOUT   : 0.030 : Mitochondrial biogenesis/repair rate (1/day)
WBIOG   : 1.50  : p53 suppression of biogenesis (PGC-1a axis)

// ---- irreversible myocyte pool -----------------------------------------
KDROS   : 4.8e-4: Maximal ROS-driven myocyte death rate (1/day)
ROS50   : 1.00  : ROS for half-maximal death (Hill n = 3)
WP53D   : 1.20  : Genotoxic-stress multiplication of death rate
WTRD    : 0.65  : ErbB2 blockade multiplication of death rate
KDNH    : 5.0e-5: Neurohormonal chronic-progression death rate (1/day)
KREG    : 2.0e-5: Cardiomyocyte renewal rate (1/day; ~0.7 %/yr)
WTR_REG : 0.80  : ErbB2 blockade suppression of renewal

// ---- reversible functional injury --------------------------------------
KFIN    : 1.0e-3: Functional injury formation rate
WF_ROS  : 0.70  : ROS weight on functional injury
WF_TOX  : 0.30  : Retained-drug weight on functional injury
KFTR    : 2.2e-3: Trastuzumab direct functional injury rate
KFOUT   : 0.012 : Functional recovery rate (1/day; t1/2 58 d)
WFIBF   : 3.0   : Fibrotic slowing of functional recovery (THE WINDOW)

// ---- fibrosis / remodelling --------------------------------------------
KFIBIN    : 8.0e-3 : Fibrosis formation rate
WFIB_DEAD : 2.5    : Myocyte-loss drive on fibrosis
WFIB_NH   : 0.30   : Neurohormonal drive on fibrosis
KFIBOUT   : 0.010  : Fibrosis resolution rate (1/day)

// ---- compensatory hypertrophy (the mask) -------------------------------
GH      : 1.60  : Hypertrophic gain on contractile deficit
HYPMAX  : 0.32  : Maximal compensatory hypertrophy (reserve ceiling)
KH      : 0.0222: Hypertrophy time constant (1/day; ~45 d)

// ---- neurohormonal ------------------------------------------------------
KNH     : 0.15  : Neurohormonal turnover (1/day)
WWS_NH  : 3.0   : Wall-stress gain on neurohormonal tone

// ---- contractility to imaging endpoints --------------------------------
EF0     : 62    : Baseline LVEF (%)
WFIBC   : 1.20  : Fibrotic penalty on contractility
EFEXP   : 0.90  : Contractility-to-LVEF exponent
GLS0    : 20.5  : Baseline |GLS| (%)
WG_FUNC : 0.90  : Functional injury weight on GLS
WG_DEAD : 1.60  : Myocyte loss weight on GLS
WG_FIB  : 0.35  : Fibrosis weight on GLS

// ---- biomarkers ---------------------------------------------------------
KTIN    : 55    : Troponin release gain
KTOUT   : 1.40  : Troponin elimination (1/day; t1/2 12 h)
WT_ROS  : 0.35  : Membrane-injury (non-lethal) troponin release weight
TNI_BASE: 1.5   : Baseline hs-cTnI (ng/L)
KBIN    : 40    : NT-proBNP formation gain
KBOUT   : 1.20  : NT-proBNP elimination (1/day)
WB_WS   : 14    : Wall-stress gain on NT-proBNP

// ---- cardioprotective pharmacology -------------------------------------
W_STA_SCAV : 0.65 : Statin suppression of ROS production (Rac1/NOX2)
W_BB_SCAV  : 0.15 : Carvedilol antioxidant contribution
W_SGLT_SCAV: 0.25 : SGLT2i suppression of ROS
W_ACE_NH   : 0.60 : ACEi/ARB suppression of neurohormonal tone
W_BB_NH    : 0.50 : Beta-blocker suppression of neurohormonal tone
W_ARNI_NH  : 1.30 : ARNI suppression of neurohormonal tone
W_ACE_REC  : 0.45 : ACEi/ARB acceleration of functional recovery
W_BB_REC   : 0.35 : Beta-blocker acceleration of functional recovery
W_ARNI_REC : 1.00 : ARNI acceleration of functional recovery
W_SGLT_REC : 0.30 : SGLT2i acceleration of functional recovery
W_ACE_FIB  : 0.25 : ACEi/ARB antifibrotic effect
W_ARNI_FIB : 0.50 : ARNI antifibrotic effect

$CMT @annotated
A1    : Doxorubicin central (mg)
A2    : Doxorubicin peripheral 1 (mg)
A3    : Doxorubicin peripheral 2 (mg)
ALIP  : Liposome-encapsulated doxorubicin (mg)
AM    : Doxorubicinol central (mg)
AMP   : Doxorubicinol peripheral (mg)
CHF   : Fast nuclear/free cardiac pool (peak-tracking)
CH    : Slow retained cardiac parent pool (AUC-tracking)
CHM   : Retained cardiac doxorubicinol pool
ADEX  : Dexrazoxane central (mg)
T2B   : Top2b protein level (1 = normal)
TR1   : Trastuzumab central (mg)
TR2   : Trastuzumab peripheral (mg)
EACE  : ACEi/ARB effect state
EBB   : Beta-blocker effect state
ESTA  : Statin effect state
EARNI : ARNI effect state
ESGLT : SGLT2 inhibitor effect state
DSB   : DNA double-strand breaks
P53   : Genotoxic-stress memory (p53/p21)
ROS   : Reactive oxygen species
MITOD : Mitochondrial deficit (0-1)
LIP   : Labile iron pool
MYO   : Viable myocyte pool (fraction of baseline) - IRREVERSIBLE
FUNC  : Reversible per-myocyte functional deficit (0-1)
FIB   : Myocardial fibrosis (0-1)
HYP   : Compensatory hypertrophy (the mask)
NH    : Neurohormonal activation
TNI   : hs-cTnI above baseline (ng/L)
BNP   : NT-proBNP (pg/mL)
AUCH  : Cumulative retained-pool exposure
CUMKILL: Cumulative fractional myocyte kill

$GLOBAL
// Individual (ETA-adjusted) parameters. They must live at file scope: variables
// declared locally in $MAIN are not visible inside $ODE, so wiring the $OMEGA
// ETAs through $MAIN locals would silently produce a population with no
// variability at all.
double iKDROS, iTOXN50, iKPIN, iKPINF, iKFIN, iKFIBIN, iCL, iEF0;

$MAIN
T2B_0   = 1.0;
MYO_0   = 1.0;
LIP_0   = 1.0;
BNP_0   = KBIN / KBOUT;

iKDROS  = KDROS  * exp(ETA_KDROS);
iTOXN50 = TOXN50 * exp(ETA_TOXN50);
iKPIN   = KPIN   * exp(ETA_KPIN);
iKPINF  = KPINF  * exp(ETA_KPINF);
iKFIN   = KFIN   * exp(ETA_KFIN);
iKFIBIN = KFIBIN * exp(ETA_KFIBIN);
iCL     = CL     * exp(ETA_CL);
iEF0    = EF0    + ETA_EF0;

$ODE
// ---------------- doxorubicin / carrier / metabolite PK -----------------
double k10 = iCL / V1;
double k12 = Q2 / V1;
double k21 = Q2 / V2;
double k13 = Q3 / V1;
double k31 = Q3 / V3;
double C1  = A1 / V1;
double CM  = AM / VM;

dxdt_A1   = -(k10 + k12 + k13) * A1 + k21 * A2 + k31 * A3 + KREL_LIP * ALIP;
dxdt_A2   = k12 * A1 - k21 * A2;
dxdt_A3   = k13 * A1 - k31 * A3;
dxdt_ALIP = -KEL_LIP * ALIP;

double klm  = CLM / VM;
double km12 = QM / VM;
double km21 = QM / VMP;
dxdt_AM  = FM * k10 * A1 - (klm + km12) * AM + km21 * AMP;
dxdt_AMP = km12 * AM - km21 * AMP;

// ---------------- two cardiac pools, two exposure metrics ---------------
dxdt_CHF = iKPINF * (C1 - CHF / KPHF);      // fast: tracks PEAKS
dxdt_CH  = iKPIN  * (C1 - CH  / KPH);       // slow: integrates AUC
dxdt_CHM = KPINM * (CM - CHM / KPHM);      // slow: accumulates over cycles

double TOXN = CHF / CHFREF;                            // nuclear driver
double TOXR = CH / CHREF + WM * CHM / CHMREF;          // redox driver
double HN   = (TOXN * TOXN) / (iTOXN50 * iTOXN50 + TOXN * TOXN);

// ---------------- dexrazoxane: Top2b protein turnover -------------------
double CDEX = ADEX / VDEX;
double FDEX = CDEX / (CDEX + KDEX50);
dxdt_ADEX = -KDEX * ADEX;
dxdt_T2B  = KT2B_SYN * (1.0 - T2B) - KT2B_DEG * FDEX * T2B;
double CHEL = KCHEL * FDEX;

// ---------------- trastuzumab ------------------------------------------
dxdt_TR1 = -(KT10 + KT12) * TR1 + KT21 * TR2;
dxdt_TR2 = KT12 * TR1 - KT21 * TR2;
double CTR = TR1 / VTR;
double ETR = CTR / (CTR + EC50TR);

// ---------------- oral cardioprotective effect states ------------------
dxdt_EACE  = KON * ((SOLVERTIME >= TON_ACE  ? TGT_ACE  : 0.0) - EACE);
dxdt_EBB   = KON * ((SOLVERTIME >= TON_BB   ? TGT_BB   : 0.0) - EBB);
dxdt_ESTA  = KON * ((SOLVERTIME >= TON_STA  ? TGT_STA  : 0.0) - ESTA);
dxdt_EARNI = KON * ((SOLVERTIME >= TON_ARNI ? TGT_ARNI : 0.0) - EARNI);
dxdt_ESGLT = KON * ((SOLVERTIME >= TON_SGLT ? TGT_SGLT : 0.0) - ESGLT);

// ---------------- Top2b-dependent DNA damage ---------------------------
dxdt_DSB = KDSB * T2B * HN - KDSBOUT / (1.0 + WTR_REP * ETR) * DSB;
dxdt_P53 = KP53 * DSB - KP53OUT * P53;

// ---------------- labile iron and ROS ----------------------------------
dxdt_LIP = KLIP * (1.0 + WLIPROS * ROS) - KLIPOUT * LIP - CHEL * LIP;
double SCAV = 1.0 / (1.0 + W_STA_SCAV * ESTA + W_BB_SCAV * EBB
                         + W_SGLT_SCAV * ESGLT);
dxdt_ROS = KROS * (WN * HN + WR * TOXR * (1.0 + FEAMP * LIP)
                   + WMITO * MITOD + WP53R * P53) * SCAV
           - KROSOUT * ROS;

// ---------------- mitochondrial deficit --------------------------------
double BIOG = 1.0 / (1.0 + WBIOG * P53);
dxdt_MITOD = KMIN * (0.6 * ROS + 0.4 * P53) * (1.0 - MITOD)
             - KMOUT * BIOG * MITOD;

// ---------------- neurohormonal tone (blocked fraction) ----------------
double NHe = NH / (1.0 + W_ACE_NH * EACE + W_BB_NH * EBB + W_ARNI_NH * EARNI);

// ---------------- irreversible myocyte loss ----------------------------
double HR = (ROS * ROS * ROS) / (ROS50 * ROS50 * ROS50 + ROS * ROS * ROS);
double kdeath = iKDROS * HR * (1.0 + WP53D * P53) * (1.0 + WTRD * ETR)
                + KDNH * NHe;
dxdt_MYO = -kdeath * MYO + KREG * (1.0 - MYO) * (1.0 - WTR_REG * ETR);

// ---------------- reversible functional deficit ------------------------
double RECOV = W_ACE_REC * EACE + W_BB_REC * EBB + W_ARNI_REC * EARNI
               + W_SGLT_REC * ESGLT;
dxdt_FUNC = (iKFIN * (WF_ROS * ROS + WF_TOX * TOXR) + KFTR * ETR)
            * (1.0 - FUNC)
            - KFOUT * (1.0 + RECOV) / (1.0 + WFIBF * FIB) * FUNC;

// ---------------- fibrosis ---------------------------------------------
double DEAD = 1.0 - MYO;
dxdt_FIB = iKFIBIN * (WFIB_DEAD * DEAD + WFIB_NH * NHe) * (1.0 - FIB)
           - KFIBOUT * (1.0 + W_ACE_FIB * EACE + W_ARNI_FIB * EARNI) * FIB;

// ---------------- compensatory hypertrophy (the mask) ------------------
double CONTraw = MYO * (1.0 - FUNC);
double DEFICIT = CONTraw < 1.0 ? 1.0 - CONTraw : 0.0;
double HYPT    = GH * DEFICIT < HYPMAX ? GH * DEFICIT : HYPMAX;
dxdt_HYP = KH * (HYPT - HYP);

double CONT  = CONTraw * (1.0 + HYP) / (1.0 + WFIBC * FIB);
double LVEFi = iEF0 * pow(CONT > 1e-6 ? CONT : 1e-6, EFEXP);
double WS    = (LVEFi < iEF0 ? (iEF0 - LVEFi) / iEF0 : 0.0) + 0.5 * FUNC;

dxdt_NH = KNH * (WWS_NH * WS - NH);

// ---------------- biomarkers -------------------------------------------
dxdt_TNI = KTIN * (1000.0 * kdeath * MYO + WT_ROS * ROS) - KTOUT * TNI;
dxdt_BNP = KBIN * (1.0 + WB_WS * WS) - KBOUT * BNP;

dxdt_AUCH    = CH;
dxdt_CUMKILL = kdeath * MYO;

$TABLE
double CONTraw_o = MYO * (1.0 - FUNC);
double CONT_o    = CONTraw_o * (1.0 + HYP) / (1.0 + WFIBC * FIB);
capture LVEF = iEF0 * pow(CONT_o > 1e-6 ? CONT_o : 1e-6, EFEXP);
capture GLS  = GLS0 * (1.0 - WG_FUNC * FUNC - WG_DEAD * (1.0 - MYO)
                       - WG_FIB * FIB);
capture CTNI = TNI + TNI_BASE;
capture NTBNP = BNP;
capture MYOCYTES = MYO;
capture FUNCDEF  = FUNC;
capture FIBROSIS = FIB;
capture MASK     = HYP;
// LVEF the patient would have WITHOUT compensatory hypertrophy = the mask
capture LVEF_NOMASK = iEF0 * pow(fmax(CONTraw_o / (1.0 + WFIBC * FIB), 1e-6),
                                 EFEXP);
// LVEF attributable to the irreversible component only (FUNC removed)
capture LVEF_IRREV = iEF0 * pow(fmax(MYO * (1.0 + HYP) / (1.0 + WFIBC * FIB),
                                     1e-6), EFEXP);
capture CPLASMA = A1 / V1;
capture CDOXOL  = AM / VM;
capture TOP2B   = T2B;

$OMEGA @annotated
// Inter-individual variability. Log-normal on the sensitivity parameters and on
// clearance; additive on baseline LVEF. Applied in $MAIN (see $GLOBAL).
ETA_KDROS  : 0.3025 : Myocyte death sensitivity (CV 55%)
ETA_TOXN50 : 0.0900 : Nuclear Hill constant (CV 30%)
ETA_KPIN   : 0.1225 : Retained-pool uptake (CV 35%)
ETA_KPINF  : 0.0900 : Nuclear-pool uptake (CV 30%)
ETA_KFIN   : 0.1600 : Functional injury sensitivity (CV 40%)
ETA_KFIBIN : 0.1600 : Fibrotic propensity (CV 40%)
ETA_CL     : 0.0784 : Doxorubicin clearance (CV 28%)
ETA_EF0    : 16.0   : Baseline LVEF (SD 4 percentage points, additive)

## =============================================================================
## R DRIVER CODE — scenarios and the eleven analyses
## Everything below is ordinary R; mread() ignores it because it lives after the
## model blocks only if sourced separately. Keep it in a companion script or
## source this file with `source()` after `mread()`.
## =============================================================================
if (FALSE) {

library(mrgsolve)
library(dplyr)
library(tidyr)

BSA <- 1.8
WT  <- 70

## ---------------------------------------------------------------------------
## dosing regimen builders
## ---------------------------------------------------------------------------
## Compartment numbering follows $CMT: 1 = A1 (free doxorubicin),
## 4 = ALIP (liposomal), 10 = ADEX (dexrazoxane), 12 = TR1 (trastuzumab).

dox_ev <- function(dose_mg_m2, ncyc, interval, tinf = 0, lipo = FALSE,
                   start = 0) {
  ev(amt = dose_mg_m2 * BSA, cmt = if (lipo) 4 else 1,
     time = start, ii = interval, addl = ncyc - 1,
     rate = if (tinf > 0) dose_mg_m2 * BSA / tinf else 0)
}

## Dexrazoxane 10:1, given 30 min BEFORE each anthracycline dose. The timing
## matters: the Top2b knockdown has to be in place before the nuclear peak.
dex_ev <- function(dox_mg_m2, ncyc, interval, ratio = 10, start = 0) {
  ev(amt = ratio * dox_mg_m2 * BSA, cmt = 10,
     time = max(0, start - 0.02), ii = interval, addl = ncyc - 1)
}

tras_ev <- function(start, weeks = 52) {
  c(ev(amt = 8 * WT, cmt = 12, time = start),
    ev(amt = 6 * WT, cmt = 12, time = start + 21,
       ii = 21, addl = floor(weeks / 3) - 2))
}

## ---------------------------------------------------------------------------
## fourteen scenarios
## ---------------------------------------------------------------------------
SCENARIOS <- list(
  DOX240 = list(
    label = "doxorubicin 60 mg/m2 q3w x4 (240 mg/m2)",
    ev = dox_ev(60, 4, 21), par = list()),
  DOX360 = list(
    label = "doxorubicin 60 mg/m2 q3w x6 (360 mg/m2)",
    ev = dox_ev(60, 6, 21), par = list()),
  DOX480 = list(
    label = "doxorubicin 60 mg/m2 q3w x8 (480 mg/m2) - reference",
    ev = dox_ev(60, 8, 21), par = list()),
  DOX480_INF72 = list(
    label = "480 mg/m2 as 72-h continuous infusions (same AUC, no peaks)",
    ev = dox_ev(60, 8, 21, tinf = 3), par = list()),
  DOX480_WEEKLY = list(
    label = "480 mg/m2 as 20 mg/m2 weekly x24",
    ev = dox_ev(20, 24, 7), par = list()),
  DOX480_PLD = list(
    label = "pegylated liposomal doxorubicin 40 mg/m2 q4w x12 (480 mg/m2)",
    ev = dox_ev(40, 12, 28, lipo = TRUE), par = list()),
  DOX480_DEX = list(
    label = "480 mg/m2 + dexrazoxane 10:1",
    ev = c(dox_ev(60, 8, 21), dex_ev(60, 8, 21)), par = list()),
  DOX480_STA = list(
    label = "480 mg/m2 + atorvastatin 40 mg from day 0 (upstream)",
    ev = dox_ev(60, 8, 21), par = list(TGT_STA = 0.75, TON_STA = 0)),
  DOX480_STA_LATE = list(
    label = "480 mg/m2 + atorvastatin started day 170 (after chemotherapy)",
    ev = dox_ev(60, 8, 21), par = list(TGT_STA = 0.75, TON_STA = 170)),
  DOX480_ACEBB = list(
    label = "480 mg/m2 + enalapril 20 mg + carvedilol 25 mg bid from day 0",
    ev = dox_ev(60, 8, 21),
    par = list(TGT_ACE = 0.8, TON_ACE = 0, TGT_BB = 0.8, TON_BB = 0)),
  DOX480_ALL = list(
    label = "480 mg/m2 + dexrazoxane + statin + ACEi/BB (full stack)",
    ev = c(dox_ev(60, 8, 21), dex_ev(60, 8, 21)),
    par = list(TGT_STA = 0.75, TON_STA = 0, TGT_ACE = 0.8, TON_ACE = 0,
               TGT_BB = 0.8, TON_BB = 0)),
  DOX240_TSEQ = list(
    label = "240 mg/m2 then trastuzumab 1 year (sequential, day 84)",
    ev = c(dox_ev(60, 4, 21), tras_ev(84)), par = list()),
  DOX240_TCONC = list(
    label = "240 mg/m2 with CONCURRENT trastuzumab 1 year",
    ev = c(dox_ev(60, 4, 21), tras_ev(0)), par = list()),
  TRAS_ALONE = list(
    label = "trastuzumab 1 year, no anthracycline",
    ev = tras_ev(0), par = list())
)

sim_scenario <- function(mod, name, tend = 730, ...) {
  sc <- SCENARIOS[[name]]
  m <- mod
  if (length(sc$par)) m <- param(m, sc$par)
  m %>% ev(sc$ev) %>% mrgsim(end = tend, delta = 1, ...) %>% as_tibble() %>%
    mutate(scenario = name)
}

## ---------------------------------------------------------------------------
## endpoint helpers
## ---------------------------------------------------------------------------
## CTRCD per ASE/ESC: LVEF fall >= 10 points AND absolute LVEF < 50%.
ctrcd <- function(d) {
  d %>% group_by(ID) %>%
    summarise(EF0 = first(LVEF), nadir = min(LVEF),
              CTRCD = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10),
              HF = any(LVEF < 40),
              GLS_drop15 = any((first(GLS) - GLS) / first(GLS) >= 0.15),
              TNI_pos = any(CTNI > 40 & time < 170),
              .groups = "drop")
}

## ---------------------------------------------------------------------------
## RESULT 1 — dose-incidence curve (is there a threshold?)
## ---------------------------------------------------------------------------
dose_response <- function(mod, n = 500, tend = 730) {
  bind_rows(lapply(c(2, 4, 5, 6, 7, 8, 10), function(ncyc) {
    out <- mod %>% ev(dox_ev(60, ncyc, 21)) %>%
      mrgsim(end = tend, delta = 1, nid = n) %>% as_tibble()
    ctrcd(out) %>% summarise(cum_dose = 60 * ncyc,
                             CTRCD_pct = 100 * mean(CTRCD),
                             HF_pct = 100 * mean(HF))
  }))
}

## ---------------------------------------------------------------------------
## RESULT 2 — cumulative dose is not the exposure metric
## ---------------------------------------------------------------------------
compare_schedules <- function(mod, n = 500, tend = 730) {
  bind_rows(lapply(c("DOX480", "DOX480_WEEKLY", "DOX480_INF72", "DOX480_PLD",
                     "DOX480_DEX"), function(k) {
    out <- sim_scenario(mod, k, tend = tend, nid = n)
    idx <- sim_scenario(mod, k, tend = tend)      # index patient, no ETAs
    ctrcd(out) %>%
      summarise(scenario = k, CTRCD_pct = 100 * mean(CTRCD),
                HF_pct = 100 * mean(HF)) %>%
      # index-patient exposure metrics: the point is that AUC_CH barely moves
      # across these rows while peak_nuclear and peak_P53 move by 10-500x
      mutate(peak_nuclear = max(idx$CHF), peak_P53 = max(idx$P53),
             peak_CH = max(idx$CH), AUC_CH = max(idx$AUCH))
  }))
}

## ---------------------------------------------------------------------------
## RESULT 3 — how big is the mask?
## ---------------------------------------------------------------------------
decompose_mask <- function(mod, name = "DOX480") {
  sim_scenario(mod, name) %>%
    filter(time %in% c(0, 60, 120, 170, 240, 365, 540, 730)) %>%
    transmute(time, LVEF, LVEF_NOMASK, LVEF_IRREV, MYOCYTES, FUNCDEF,
              FIBROSIS, MASK,
              mask_points = LVEF - LVEF_NOMASK,
              reversible_points = LVEF_IRREV - LVEF)
}

## ---------------------------------------------------------------------------
## RESULT 4 — biomarker lead times
## ---------------------------------------------------------------------------
lead_times <- function(mod, name = "DOX480", n = 500) {
  out <- sim_scenario(mod, name, nid = n)
  first_day <- function(d, cond) {
    d %>% group_by(ID) %>% filter({{ cond }}) %>%
      summarise(day = min(time), .groups = "drop")
  }
  tni <- first_day(out, CTNI > 14)
  gls <- out %>% group_by(ID) %>%
    filter((first(GLS) - GLS) / first(GLS) >= 0.15) %>%
    summarise(day_gls = min(time), .groups = "drop")
  ef <- out %>% group_by(ID) %>%
    filter(LVEF < 50, (first(LVEF) - LVEF) >= 10) %>%
    summarise(day_ef = min(time), .groups = "drop")
  list(tni = tni, gls = gls, ef = ef,
       gls_lead = inner_join(gls, ef, "ID") %>%
         summarise(median_lead_days = median(day_ef - day_gls)),
       tni_npv = ctrcd(out) %>% group_by(TNI_pos) %>%
         summarise(CTRCD_pct = 100 * mean(CTRCD), n = n()))
}

## ---------------------------------------------------------------------------
## RESULT 5/6 — dexrazoxane mechanism and dose-sparing equivalence
## ---------------------------------------------------------------------------
dex_sparing <- function(mod, ref_cum = 240, tend = 730) {
  target <- sim_scenario(mod, "DOX240", tend = tend) %>%
    filter(time == 365) %>% pull(LVEF)
  res <- bind_rows(lapply(4:14, function(ncyc) {
    out <- mod %>% ev(c(dox_ev(60, ncyc, 21), dex_ev(60, ncyc, 21))) %>%
      mrgsim(end = tend, delta = 1) %>% as_tibble()
    tibble(cum_dose = 60 * ncyc,
           LVEF_12m = out$LVEF[out$time == 365])
  }))
  list(unprotected_240_LVEF12m = target, protected = res,
       iso_effective_dose = res$cum_dose[which(res$LVEF_12m <= target)[1]])
}

## ---------------------------------------------------------------------------
## RESULT 7 — upstream versus downstream protection
## ---------------------------------------------------------------------------
upstream_vs_downstream <- function(mod, n = 500) {
  combos <- list(
    `no protection`            = list(),
    `statin during chemo`      = list(TGT_STA = 0.75, TON_STA = 0),
    `statin after chemo`       = list(TGT_STA = 0.75, TON_STA = 170),
    `ACEi+BB during chemo`     = list(TGT_ACE = 0.8, TON_ACE = 0,
                                      TGT_BB = 0.8, TON_BB = 0),
    `ACEi+BB after chemo`      = list(TGT_ACE = 0.8, TON_ACE = 170,
                                      TGT_BB = 0.8, TON_BB = 170),
    `ARNI after chemo`         = list(TGT_ARNI = 0.9, TON_ARNI = 170),
    `SGLT2i during chemo`      = list(TGT_SGLT = 0.7, TON_SGLT = 0))
  bind_rows(lapply(names(combos), function(nm) {
    m <- if (length(combos[[nm]])) param(mod, combos[[nm]]) else mod
    out <- m %>% ev(dox_ev(60, 8, 21)) %>%
      mrgsim(end = 730, delta = 1, nid = n) %>% as_tibble()
    ctrcd(out) %>% summarise(strategy = nm, CTRCD_pct = 100 * mean(CTRCD),
                             mean_nadir = mean(nadir))
  }))
}

## ---------------------------------------------------------------------------
## RESULT 9 — anthracycline x trastuzumab: additive or synergistic?
## ---------------------------------------------------------------------------
trastuzumab_interaction <- function(mod, n = 500) {
  get <- function(k) {
    out <- sim_scenario(mod, k, nid = n)
    out %>% group_by(ID) %>%
      summarise(maxfall = first(LVEF) - min(LVEF),
                regain = last(LVEF) - min(LVEF), .groups = "drop") %>%
      summarise(scenario = k, maxfall = mean(maxfall), regain = mean(regain))
  }
  res <- bind_rows(lapply(c("DOX240", "TRAS_ALONE", "DOX240_TSEQ",
                            "DOX240_TCONC"), get))
  additive <- res$maxfall[1] + res$maxfall[2]
  res %>% mutate(additive_prediction = additive,
                 excess = maxfall - additive)
}

## ---------------------------------------------------------------------------
## RESULT 10 — the reversibility window (a fibrosis clock)
## ---------------------------------------------------------------------------
reversibility_window <- function(mod, n = 500, starts = c(60, 90, 120, 180,
                                                          270, 365, 540)) {
  bind_rows(lapply(starts, function(t0) {
    m <- param(mod, list(TGT_ACE = 0.8, TON_ACE = t0,
                         TGT_BB = 0.8, TON_BB = t0))
    out <- m %>% ev(dox_ev(60, 8, 21)) %>%
      mrgsim(end = t0 + 365, delta = 1, nid = n) %>% as_tibble()
    out %>% group_by(ID) %>%
      summarise(had_ctrcd = any(LVEF < 50 & (first(LVEF) - LVEF) >= 10 &
                                  time <= 365),
                EF_start = LVEF[which.min(abs(time - t0))],
                EF_end = last(LVEF),
                FIB_start = FIBROSIS[which.min(abs(time - t0))],
                EF_base = first(LVEF), .groups = "drop") %>%
      filter(had_ctrcd) %>%
      summarise(start_day = t0, n = n(), EF_start = mean(EF_start),
                EF_12mo_later = mean(EF_end),
                recovered_pct = 100 * mean(EF_end >= 50 |
                                             (EF_base - EF_end) < 5),
                FIB_at_start = mean(FIB_start))
  }))
}

## ---------------------------------------------------------------------------
## RESULT 11 — CBR1/AKR1C3 metaboliser status
## ---------------------------------------------------------------------------
metaboliser_status <- function(mod, n = 500) {
  bind_rows(lapply(c(0.12, 0.18, 0.25, 0.34, 0.45), function(fm) {
    out <- param(mod, FM = fm) %>% ev(dox_ev(60, 8, 21)) %>%
      mrgsim(end = 730, delta = 1, nid = n) %>% as_tibble()
    ctrcd(out) %>% summarise(FM = fm, CTRCD_pct = 100 * mean(CTRCD),
                             mean_nadir = mean(nadir))
  }))
}

## ---------------------------------------------------------------------------
## virtual-population incidence table across all scenarios
## ---------------------------------------------------------------------------
vpop_incidence <- function(mod, n = 500) {
  bind_rows(lapply(names(SCENARIOS), function(k) {
    out <- sim_scenario(mod, k, nid = n)
    ctrcd(out) %>% summarise(scenario = k, label = SCENARIOS[[k]]$label,
                             CTRCD_pct = 100 * mean(CTRCD),
                             HF_pct = 100 * mean(HF),
                             GLS15_pct = 100 * mean(GLS_drop15),
                             mean_nadir = mean(nadir))
  }))
}

## ---------------------------------------------------------------------------
## high-risk phenotype (age >= 70 + hypertension + prior mediastinal RT)
## ---------------------------------------------------------------------------
high_risk_param <- function(mod) {
  param(mod, KDROS = 4.8e-4 * 1.9, KFIN = 1.0e-3 * 1.5,
        KFIBIN = 8.0e-3 * 1.6, EF0 = 58, KREG = 2.0e-5 * 0.4)
}

## ---------------------------------------------------------------------------
## run everything
## ---------------------------------------------------------------------------
run_all <- function(mod, n = 500) {
  list(dose_response          = dose_response(mod, n),
       schedules              = compare_schedules(mod, n),
       mask                   = decompose_mask(mod),
       lead_times             = lead_times(mod, n = n),
       dexrazoxane_sparing    = dex_sparing(mod),
       upstream_vs_downstream = upstream_vs_downstream(mod, n),
       trastuzumab            = trastuzumab_interaction(mod, n),
       reversibility_window   = reversibility_window(mod, n),
       metaboliser            = metaboliser_status(mod, n),
       population             = vpop_incidence(mod, n))
}

## mod <- mread("aic", "aic_mrgsolve_model.R")
## res <- run_all(mod, n = 500)

}
