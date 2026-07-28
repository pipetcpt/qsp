## =============================================================================
## Tardive Dyskinesia (TD) — mrgsolve QSP model
## 지연성 운동이상증 정량적 시스템 약리학 모델
##
## 40 ODEs (19 PK/exposure + 21 disease · basal-ganglia · endpoint) ·
## 10 treatment scenarios · 11 analysis functions.  Time unit = DAYS.
##
## -----------------------------------------------------------------------------
## WHAT THIS MODEL ASSERTS
## -----------------------------------------------------------------------------
## Tardive dyskinesia is usually described as a *consequence*: chronic D2
## blockade upregulates postsynaptic D2 receptors, the receptors become
## supersensitive, and involuntary movements follow. That description cannot
## explain the three things clinicians actually have to deal with:
##
##   (i)   why stopping the offending drug makes the movements WORSE first,
##   (ii)  why stopping works in some patients and in others never works again,
##   (iii) why a VMAT2 inhibitor reliably lowers AIMS and yet the disease is
##         exactly where it was when the drug is stopped.
##
## This model answers all three with ONE structural claim: TD is carried by
## TWO memories with different time constants and different reversibility, and
## the drug levers act on different ones.
##
##   RUP  — postsynaptic D2 supersensitivity (receptor density, D2High
##          fraction, Gi/beta-arrestin gain).  Driven ONLY above an occupancy
##          threshold: KIN_R * Hill(OCC; OCC50_R = 0.70, n = 4).
##          Fully reversible with tau = 300 d (slowed by damage:
##          KOUT_R / (1 + 0.8 * SDAM)).
##
##   SDAM — structural striatal damage (cholinergic/parvalbumin interneuron
##          loss, spine loss, GAD67 down, gliosis).  Driven by
##          0.55*ROS + 0.45*RUP, and — the essential term — it feeds its own
##          driver back through mitochondrial ROS:
##              ros_drive += W_S_ROS * Hill(SDAM; 0.50, n = 4),  W_S_ROS = 2.2
##          Because that feedback is SATURATING (negligible below SDAM ~ 0.3,
##          near-maximal above ~0.7), the SDAM equation is BISTABLE:
##              lower stable state  SDAM = 0        (repair wins)
##              unstable threshold  SDAM ~ 0.50-0.55
##              upper stable state  SDAM ~ 0.85-0.90  (self-sustaining)
##
## The observed movements are not read off either state directly. They are the
## output of a basal-ganglia loop driven by an effective over-stimulation term
##
##      D2STIM = DA_SYN * (1 + 2.2*RUP) * (1 - 0.55*OCC) * (1+0.5*E_ACH)
##                                                        / (1+0.35*E_AMA)
##      EXC    = max(D2STIM - 1, 0) + 1.2*SDAM
##      IND -> GPe -> STN -> GPi -> THAL -> AIMS = 26*Hill(THAL-1; 0.55, 2)
##
## Two consequences fall out of the algebra and are the reason the model is
## written this way:
##
##   * the (1 - 0.55*OCC) factor means ONGOING BLOCKADE MASKS THE SYMPTOM it
##     is creating — remove the blockade and AIMS rises before it falls;
##   * the +1.2*SDAM term is OCCUPANCY-INDEPENDENT — once the latch is set,
##     there is an AIMS floor that no manipulation of the offending drug can
##     reach.
##
## And the two available levers sit on opposite sides of the same synapse:
## lowering D2 blockade (dose reduction / clozapine switch) removes the drive
## to RUP but unmasks the symptom and costs psychosis control; lowering
## dopamine supply (VMAT2 inhibition) removes the symptom but touches neither
## RUP nor SDAM and costs parkinsonism / depression / sedation.
##
## -----------------------------------------------------------------------------
## ELEVEN COMPUTED RESULTS (functions, not prose)
## -----------------------------------------------------------------------------
## All numbers below are OUTPUTS of these equations, verified against an
## independent numpy/scipy transcription (td_reference_check.py), not
## literature values. Index patient = 58 y, risperidone-equivalent 8 mg/day,
## FGA-level oxidative burden (RISK_FGA = 1.6), OCC = 0.86, AIMS(day 730) =
## 10.59.
##
##  1. TD_scenarios() — 10 regimens over 5 years. Withdrawal at y2 takes AIMS
##     10.59 -> 6.16 at y5 but PSYCH 0.30 -> 1.00; valbenazine 80 mg takes
##     AIMS to 8.37 with PSYCH unchanged (0.31) but PARK 0.27 -> 0.61.
##
##  2. TD_exposure_threshold() — CUMULATIVE DOSE IS NOT THE EXPOSURE METRIC.
##     16 mg for 300 d (280 occupancy-days at OCC 0.93) leaves persistent TD
##     (AIMS 4.87 six years after stopping); 2 mg for 1500 d (813
##     occupancy-days at OCC 0.54) resolves completely (AIMS 0.18). Nearly
##     3x the occupancy-time is harmless if it is delivered below the
##     plasticity threshold.
##
##  3. TD_reversibility_window() — withdrawal succeeds if it happens by day
##     285 (AIMS at +6 y = 0.90) and fails from day 290 (4.82). The point of
##     no return is a 5-day-wide step, and it occurs while SDAM is only 0.34
##     — well BELOW the 0.50-0.55 self-sustaining threshold — because RUP
##     keeps feeding damage for months after the drug is gone. The window
##     cannot be read off the patient's current state.
##
##  4. TD_withdrawal_crossover() — the paradox, quantified. Stopping raises
##     AIMS from 10.59 to a peak of 12.02 at day +30 and only crosses below
##     the continue-treatment curve at day +135. Halving the dose crosses at
##     +237 d, a 25% reduction at +318 d. Slower reduction = smaller peak,
##     later crossover, smaller ultimate benefit.
##
##  5. TD_vmat2_dose_response() — valbenazine 80 mg/day gives AIMS -3.51
##     (-33%) at week 6; deutetrabenazine 36 mg/day gives -3.09 (-29%);
##     both with RUP and SDAM unchanged at year 5 (0.885).
##
##  6. TD_cyp2d6_panel() — the same 80 mg gives -22% in a CYP2D6 ultrarapid
##     metabolizer and -54% in a poor metabolizer, with QTc 1.20 vs 3.84 ms
##     and PARK 0.48 vs 0.74: genotype moves both efficacy and toxicity along
##     one exposure axis, which is why the label caps the dose in PMs.
##
##  7. TD_suppression_vs_modification() — 2 y of valbenazine then washout:
##     AIMS 12.36 -> 8.88 on drug, back to 12.33 eight weeks after stopping,
##     with RUP 1.41 vs 1.45 untreated and SDAM identical. Clozapine switch:
##     AIMS 7.70 at y4 and STILL 6.00 at y6, RUP 1.45 -> 0.30. Suppression and
##     modification are different verbs.
##
##  8. TD_combination_interaction() — clozapine x valbenazine. At week 6 the
##     switch alone makes AIMS WORSE (+1.32) and the combination makes it
##     better (-2.35): the VMAT2 inhibitor is a BRIDGE that pays for the
##     withdrawal-emergent worsening of the switch. By year 7 the interaction
##     term is +2.54 (sub-additive) because both are pushing against the same
##     structural floor.
##
##  9. TD_opposed_levers() — the frontier. At day 1095, stopping the
##     antipsychotic buys AIMS 11.82 -> 9.44 at the price of PSYCH 1.00;
##     valbenazine 80 buys 11.82 -> 8.26 at PSYCH 0.27 but PARK 0.52,
##     DEPR 0.28, ADHER 0.69. Clozapine is the only strategy that improves
##     AIMS, RUP and psychosis control at once.
##
## 10. TD_risk_scan() — identical regimen, different host: latch day 585 (25 y,
##     SGA) vs 366 (75 y, FGA) vs 357 (58 y + diabetes). Benztropine takes
##     AIMS at 1 y from 7.10 to 10.69. Ginkgo (EGb761) is the only adjunct
##     that moves the latch itself: latch day 400 -> 780, AIMS at y1 7.10 ->
##     4.64 — and by y5 the difference has shrunk to 12.61 vs 12.43, so the
##     model says antioxidants DELAY the latch rather than prevent it.
##
## 11. TD_latch_bistability() — the nullcline, printed: three roots, with the
##     unstable one between 0.50 and 0.55 and the upper stable one between
##     0.85 and 0.90.
##
## -----------------------------------------------------------------------------
## CALIBRATION ANCHORS (literature -> parameter)
## -----------------------------------------------------------------------------
##  * Risperidone active-moiety PK: CL/F 120 L/day, V2 100 L, t1/2 ~ 20 h;
##    4 mg/day -> Cave ~ 33 ng/mL.                        (KA_AP, CL_AP, V2_AP)
##  * Striatal D2 occupancy vs active-moiety concentration, EC50 ~ 12 ng/mL;
##    4 mg -> ~73%, 8 mg -> ~86%.                            (EC50_D2, HILL_D2)
##  * Therapeutic occupancy threshold ~65%, EPS threshold ~78-80%.
##                                                     (OCC50_P, OCC50_PK)
##  * Clozapine 350 mg/day -> ~390 ng/mL with D2 occupancy ~28%.
##                                        (CL_CLZ, V_CLZ, EC50_D2_CLZ)
##  * Valbenazine 80 mg/day -> NBI-98782 Cave ~ 24 ng/mL; CYP2D6 PM roughly
##    doubles metabolite exposure.                (FM_VAL, CL_NBI, CYP2D6)
##  * Deutetrabenazine 36 mg/day -> total (alpha+beta)-HTBZ ~ 15 ng/mL;
##    QTc effect small (~4 ms).                     (FM_DTB, CL_HTB, QT_HTB)
##  * KINECT-3: valbenazine 80 mg, AIMS -3.2 at week 6 from baseline ~10.
##    Model: -3.51 from 10.59.                               (EC50_VMAT_NBI)
##  * ARM-TD / AIM-TD: deutetrabenazine 36 mg, AIMS ~ -3.0 at week 12.
##    Model: -3.09 at week 6.                                (EC50_VMAT_HTB)
##  * Annual TD incidence ~5-6%/y (FGA) and ~3-4%/y (SGA), rising to
##    25-30%/y over age 55.                    (RISK_FGA, AGE in RISKMOD)
##  * Withdrawal-emergent dyskinesia peaks within ~4-6 weeks of dose
##    reduction.                                       (MASK, K_AIMS, KBG)
##  * Persistence: roughly a third of cases remit after withdrawal, the rest
##    persist for years.                            (KOUT_S, S50, W_S_ROS)
##  * Ginkgo biloba EGb761 240 mg/day: AIMS improvement in RCT.
##                                             (ANTIOX_MAX, KOUT_ROS)
##
## -----------------------------------------------------------------------------
## DOSING CONVENTION
## -----------------------------------------------------------------------------
## Every drug is given as a ZERO-ORDER DAILY INPUT (DOSE mg/day entering the
## absorption compartment continuously) rather than as discrete events. For
## chronic daily therapy over years this is equivalent to bolus dosing at the
## resolution of every read-out in this file, and it makes the R model and the
## Python reference implementation numerically identical. Regimen changes are
## smooth unit steps, step(t; t0, tau) = 0.5*(1 + tanh((t - t0)/tau)) with
## TAPER = 1 day, so the solver never sees a discontinuity.
## To use event-based dosing instead, set DOSE_* = 0 and pass
##   ev(amt = 4, ii = 1, addl = 1824, cmt = "AP_gut")
## =============================================================================

library(mrgsolve)
suppressMessages(library(dplyr))

td_code <- '
$PARAM @annotated
// ---- antipsychotic PK (risperidone active-moiety equivalents) -------------
KA_AP    : 20.0  : absorption rate constant (1/day)
CL_AP    : 120.0 : apparent clearance (L/day)
V2_AP    : 100.0 : central volume (L)
V3_AP    : 80.0  : peripheral volume (L)
Q_AP     : 30.0  : intercompartmental clearance (L/day)
KA_LAI   : 0.16  : long-acting-injectable release rate (1/day)
F_LAI    : 0.0   : fraction of antipsychotic dose given as LAI (-)
KE0_AP   : 12.0  : striatal effect-site equilibration (1/day)
EC50_D2  : 12.0  : concentration for 50% striatal D2 occupancy (ng/mL)
HILL_D2  : 1.25  : Hill coefficient for D2 occupancy (-)
// ---- clozapine ------------------------------------------------------------
KA_CLZ   : 15.0  : clozapine absorption (1/day)
CL_CLZ   : 900.0 : clozapine clearance (L/day)
V_CLZ    : 500.0 : clozapine volume (L)
EC50_D2_CLZ  : 900.0 : clozapine concentration for 50% D2 occupancy (ng/mL)
EC50_CLZ_EFF : 250.0 : clozapine concentration for 50% non-D2 efficacy (ng/mL)
// ---- valbenazine -> NBI-98782 --------------------------------------------
KA_VAL   : 10.0  : valbenazine absorption (1/day)
CL_VAL   : 500.0 : valbenazine clearance (L/day)
V_VAL    : 500.0 : valbenazine volume (L)
FM_VAL   : 0.30  : fraction converted to NBI-98782 (-)
CL_NBI   : 1000.0 : NBI-98782 clearance (L/day)
V_NBI    : 900.0 : NBI-98782 central volume (L)
Q_NBI    : 200.0 : NBI-98782 intercompartmental clearance (L/day)
V3_NBI   : 600.0 : NBI-98782 peripheral volume (L)
EC50_VMAT_NBI : 26.0 : NBI-98782 concentration for 50% VMAT2 occupancy (ng/mL)
// ---- deutetrabenazine -> (alpha+beta)-HTBZ -------------------------------
KA_DTB   : 12.0  : deutetrabenazine absorption (1/day)
FM_DTB   : 0.50  : fraction converted to active HTBZ metabolites (-)
CL_HTB   : 1200.0 : HTBZ clearance (L/day)
V_HTB    : 800.0 : HTBZ central volume (L)
Q_HTB    : 150.0 : HTBZ intercompartmental clearance (L/day)
V3_HTB   : 500.0 : HTBZ peripheral volume (L)
EC50_VMAT_HTB : 19.0 : HTBZ concentration for 50% VMAT2 occupancy (ng/mL)
CYP2D6   : 1.0   : CYP2D6 clearance multiplier for NBI/HTBZ (-)
// ---- adjuncts -------------------------------------------------------------
KA_AMA   : 8.0   : amantadine absorption (1/day)
CL_AMA   : 400.0 : amantadine clearance (L/day)
V_AMA    : 350.0 : amantadine volume (L)
EC50_AMA : 400.0 : amantadine EC50 (ng/mL)
KA_ACH   : 10.0  : benztropine absorption (1/day)
CL_ACH   : 1000.0 : benztropine clearance (L/day)
V_ACH    : 800.0 : benztropine volume (L)
EC50_ACH : 1.5   : benztropine EC50 (ng/mL)
KIN_GKB  : 0.143 : ginkgo effect-compartment on-rate (1/day)
KOUT_GKB : 0.143 : ginkgo effect-compartment off-rate (1/day)
ANTIOX_MAX : 0.85 : maximal fractional increase in ROS elimination (-)
KOUT_BONT : 0.0111 : botulinum effect decay (1/day)
BONT_MAX : 0.35  : maximal local AIMS attenuation by botulinum (-)
BONT_INT : 90.0  : botulinum injection interval (day)
// ---- nigrostriatal dopamine handling -------------------------------------
SYN0     : 1.0   : basal dopamine synthesis rate (unit/day)
A_AUTO   : 0.60  : autoreceptor-blockade gain on synthesis (-)
KVES     : 5.0   : VMAT2-mediated vesicular uptake (1/day)
KMAO     : 3.0   : cytosolic MAO/COMT degradation (1/day)
KREL     : 2.0   : vesicular release rate constant (1/day)
KDAT     : 8.0   : DAT reuptake (1/day)
KDIFF    : 1.0   : synaptic diffusion loss (1/day)
KBACK    : 7.2   : recycling of reuptaken dopamine to cytosol (1/day)
A_FIRE   : 0.50  : autoreceptor-blockade gain on firing (-)
KIN_DEP  : 0.0167 : depolarization-block on-rate (1/day)
KOUT_DEP : 0.0167 : depolarization-block off-rate (1/day)
DEPOL_MAX : 0.45 : maximal depolarization block (-)
DA_SYN0  : 0.13889 : analytic baseline synaptic dopamine (unit)
DA_CYT0  : 0.25  : analytic baseline cytosolic dopamine (unit)
// ---- postsynaptic D2 supersensitivity ------------------------------------
KIN_R    : 0.0035 : D2 upregulation zero-order rate (1/day)
KOUT_R   : 0.003333 : D2 upregulation loss (1/day, tau = 300 d)
OCC50_R  : 0.70  : occupancy giving half-maximal plasticity drive (-)
HILL_R   : 4.0   : steepness of the plasticity threshold (-)
W_ACH_R  : 0.50  : anticholinergic amplification of upregulation (-)
W_SDAM_R : 0.80  : slowing of recovery by structural damage (-)
// ---- oxidative stress ----------------------------------------------------
KIN_ROS  : 0.03333 : ROS formation rate constant (1/day)
KOUT_ROS : 0.03333 : ROS elimination rate constant (1/day)
W_CYT_ROS : 0.80 : weight of cytosolic dopamine on ROS drive (-)
W_OCC_ROS : 0.50 : weight of D2 blockade on ROS drive (-)
W_RUP_ROS : 0.25 : weight of upregulation on ROS drive (-)
W_S_ROS  : 2.20  : damage-to-ROS positive feedback gain (-)
S_ROS50  : 0.50  : damage at which the ROS feedback is half-maximal (-)
HILL_S_ROS : 4.0 : steepness of the damage-to-ROS feedback (-)
ROS_CONST : 0.0  : additional constitutive ROS drive (-)
RISK_FGA : 1.0   : first-generation-antipsychotic oxidative multiplier (-)
W_AMA_ROS : 0.30 : amantadine contribution to ROS elimination (-)
// ---- structural damage latch ---------------------------------------------
KIN_S    : 0.0025 : damage accrual rate (1/day)
KOUT_S   : 0.00030 : damage repair rate (1/day, tau ~ 9 y)
W_ROS_S  : 0.55  : weight of ROS on damage drive (-)
W_RUP_S  : 0.45  : weight of upregulation on damage drive (-)
S50      : 0.85  : damage drive giving half-maximal accrual (-)
HILL_S   : 8.0   : steepness of the damage threshold (-)
// ---- postsynaptic stimulation and basal ganglia --------------------------
G_RUP    : 2.20  : gain of upregulation on effective D2 stimulation (-)
MASK     : 0.55  : symptomatic masking by ongoing D2 blockade (-)
G_ACH    : 0.50  : anticholinergic amplification of stimulation (-)
G_AMA    : 0.35  : amantadine attenuation of stimulation (-)
G_S      : 1.20  : occupancy-independent contribution of damage (-)
K_IND    : 0.80  : D2 inhibition of indirect-pathway MSNs (-)
K_DIR    : 0.50  : D1 facilitation of direct-pathway output (-)
KBG      : 0.50  : basal-ganglia loop relaxation rate (1/day)
K_THAL_GPI : 0.80 : pallidal weight on thalamic gating (-)
AIMS_MAX : 26.0  : maximal AIMS dyskinesia subscore (points)
TH50     : 0.55  : thalamic excess giving half-maximal AIMS (-)
HILL_TH  : 2.0   : steepness of the AIMS transfer function (-)
K_AIMS   : 0.1429 : AIMS observation/rating lag (1/day)
// ---- psychiatric, motor, mood, cardiac endpoints -------------------------
OCC50_P  : 0.58  : occupancy for half-maximal psychosis control (-)
HILL_P   : 3.0   : steepness of psychosis control (-)
E_MAX_P  : 1.25  : maximal attainable psychosis control (-)
W_SUPER_P : 0.20 : erosion of control by supersensitivity (-)
W_CLZ_P  : 0.75  : clozapine non-D2 contribution to control (-)
K_PSY    : 0.0333 : psychosis state relaxation (1/day)
OCC50_PK : 0.78  : occupancy threshold for parkinsonism (-)
HILL_PK  : 4.0   : steepness of the parkinsonism threshold (-)
VOCC50_PK : 0.70 : VMAT2 occupancy threshold for parkinsonism (-)
HILL_VPK : 3.0   : steepness of VMAT2 parkinsonism (-)
W_RUP_PK : 0.30  : protection from parkinsonism by supersensitivity (-)
K_PARK   : 0.0714 : parkinsonism relaxation (1/day)
VOCC50_DEP : 0.75 : VMAT2 occupancy threshold for depression (-)
HILL_VDEP : 2.0  : steepness of depletion-related depression (-)
W_VMAT_DEP : 0.60 : maximal depletion-related mood burden (-)
W_OCC_DEP : 0.20 : D2-blockade contribution to mood burden (-)
K_DEPR   : 0.0476 : mood state relaxation (1/day)
QT_NBI   : 0.080 : QTc slope for NBI-98782 (ms per ng/mL)
QT_HTB   : 0.250 : QTc slope for HTBZ (ms per ng/mL)
K_QTC    : 1.0   : QTc equilibration (1/day)
W_PARK_AD : 0.50 : parkinsonism effect on adherence (-)
W_DEPR_AD : 0.40 : mood effect on adherence (-)
W_PSY_AD : 0.30  : psychosis effect on adherence (-)
K_ADHER  : 0.0333 : adherence relaxation (1/day)
// ---- clinician dose-escalation policy ------------------------------------
ESC_ON   : 0.0   : enable the escalation policy (0/1)
K_ESC    : 0.004 : escalation gain (1/day)
PSY_TARGET : 0.20 : psychosis burden that triggers escalation (-)
DOSEMULT_MAX : 2.0 : maximal dose multiplier (-)
// ---- host risk modifiers -------------------------------------------------
AGE      : 58.0  : age (year)
ESTROGEN : 0.0   : estrogen-replete protection flag (0/1)
DM_RISK  : 0.0   : diabetes/metabolic risk flag (0/1)
GEN_RISK : 1.0   : genetic risk multiplier (-)
// ---- regimen -------------------------------------------------------------
DOSE_AP  : 8.0   : antipsychotic dose before switch (mg/day risperidone-eq)
DOSE_AP2 : 8.0   : antipsychotic dose after switch (mg/day risperidone-eq)
TSW_AP   : 1000000 : day of antipsychotic dose change (day)
DOSE_CLZ : 0.0   : clozapine dose (mg/day)
TSTART_CLZ : 1000000 : clozapine start day (day)
DOSE_VAL : 0.0   : valbenazine dose (mg/day)
TSTART_VAL : 1000000 : valbenazine start day (day)
TSTOP_VAL : 1000000 : valbenazine stop day (day)
DOSE_DTB : 0.0   : deutetrabenazine dose (mg/day)
TSTART_DTB : 1000000 : deutetrabenazine start day (day)
TSTOP_DTB : 1000000 : deutetrabenazine stop day (day)
DOSE_AMA : 0.0   : amantadine dose (mg/day)
TSTART_AMA : 1000000 : amantadine start day (day)
DOSE_ACH : 0.0   : benztropine dose (mg/day)
TSTART_ACH : 1000000 : benztropine start day (day)
GKB_ON   : 0.0   : ginkgo EGb761 flag (0/1)
TSTART_GKB : 1000000 : ginkgo start day (day)
BONT_ON  : 0.0   : botulinum toxin flag (0/1)
TSTART_BONT : 1000000 : botulinum start day (day)
TAPER    : 1.0   : smoothing time for regimen changes (day)

$CMT @annotated
AP_gut   : oral antipsychotic absorption compartment (mg)
AP_cen   : antipsychotic central compartment (mg)
AP_per   : antipsychotic peripheral compartment (mg)
AP_lai   : long-acting-injectable depot (mg)
VAL_gut  : valbenazine absorption compartment (mg)
VAL_cen  : valbenazine central compartment (mg)
NBI_cen  : NBI-98782 central compartment (mg)
NBI_per  : NBI-98782 peripheral compartment (mg)
DTB_gut  : deutetrabenazine absorption compartment (mg)
HTB_cen  : HTBZ central compartment (mg)
HTB_per  : HTBZ peripheral compartment (mg)
CLZ_gut  : clozapine absorption compartment (mg)
CLZ_cen  : clozapine central compartment (mg)
AMA_gut  : amantadine absorption compartment (mg)
AMA_cen  : amantadine central compartment (mg)
ACH_gut  : benztropine absorption compartment (mg)
ACH_cen  : benztropine central compartment (mg)
GKB      : ginkgo antioxidant effect compartment (-)
BONT     : botulinum local effect compartment (-)
CE_AP    : striatal effect-site concentration (ng/mL)
DA_CYT   : cytosolic dopamine (unit)
DA_VES   : vesicular dopamine (unit)
DA_SYN   : synaptic dopamine (unit)
DEPOL    : depolarization block of dopamine neurons (-)
RUP      : postsynaptic D2 supersensitivity (-)
ROS      : striatal oxidative stress index (-)
SDAM     : structural striatal damage (-)
IND      : indirect-pathway (striatopallidal) activity (-)
GPE      : external pallidal activity (-)
STN      : subthalamic activity (-)
GPI      : internal pallidal / SNr output (-)
THAL     : thalamic gating (-)
AIMS     : AIMS dyskinesia subscore (points)
PSYCH    : psychosis burden (-)
PARK     : drug-induced parkinsonism burden (-)
DEPR     : depressive/mood burden (-)
QTC      : QTc change from baseline (ms)
ADHER    : adherence (-)
DOSEMULT : clinician dose multiplier (-)
CUMOCC   : cumulative D2 occupancy-time (occupancy-day)

$GLOBAL
double hillf(double x, double x50, double n) {
  if (x < 0.0) x = 0.0;
  double xn = pow(x, n);
  double d  = xn + pow(x50, n);
  return (d > 0.0) ? xn / d : 0.0;
}
double stepf(double t, double t0, double tau) {
  return 0.5 * (1.0 + tanh((t - t0) / tau));
}

$MAIN
DA_CYT_0   = DA_CYT0;
DA_VES_0   = 0.625;
DA_SYN_0   = DA_SYN0;
IND_0      = 1.0;
GPE_0      = 1.0;
STN_0      = 1.0;
GPI_0      = 1.0;
THAL_0     = 1.0;
ADHER_0    = 1.0;
DOSEMULT_0 = 1.0;

$ODE
// ---------------- regimen drivers ----------------------------------------
double ap_dose = (DOSE_AP + (DOSE_AP2 - DOSE_AP) * stepf(SOLVERTIME, TSW_AP, TAPER))
                 * DOSEMULT;
double adh     = ADHER; if (adh < 0.0) adh = 0.0; if (adh > 1.0) adh = 1.0;
double ap_oral = ap_dose * (1.0 - F_LAI) * adh;
double ap_dep  = ap_dose * F_LAI;
double clz_dose = DOSE_CLZ * stepf(SOLVERTIME, TSTART_CLZ, TAPER) * adh;
double val_dose = DOSE_VAL * (stepf(SOLVERTIME, TSTART_VAL, TAPER)
                              - stepf(SOLVERTIME, TSTOP_VAL, TAPER));
double dtb_dose = DOSE_DTB * (stepf(SOLVERTIME, TSTART_DTB, TAPER)
                              - stepf(SOLVERTIME, TSTOP_DTB, TAPER));
double ama_dose = DOSE_AMA * stepf(SOLVERTIME, TSTART_AMA, TAPER);
double ach_dose = DOSE_ACH * stepf(SOLVERTIME, TSTART_ACH, TAPER);
double gkb_in   = GKB_ON  * stepf(SOLVERTIME, TSTART_GKB, TAPER);
double bont_in  = BONT_ON * stepf(SOLVERTIME, TSTART_BONT, TAPER);

// ---------------- pharmacokinetics ---------------------------------------
double C_AP  = AP_cen  / V2_AP * 1000.0;      // mg/L -> ng/mL
double C_APp = AP_per  / V3_AP * 1000.0;
dxdt_AP_gut = ap_oral - KA_AP * AP_gut;
dxdt_AP_lai = ap_dep  - KA_LAI * AP_lai;
dxdt_AP_cen = KA_AP * AP_gut + KA_LAI * AP_lai
              - CL_AP * C_AP / 1000.0 - Q_AP * (C_AP - C_APp) / 1000.0;
dxdt_AP_per = Q_AP * (C_AP - C_APp) / 1000.0;

double C_CLZ = CLZ_cen / V_CLZ * 1000.0;
dxdt_CLZ_gut = clz_dose - KA_CLZ * CLZ_gut;
dxdt_CLZ_cen = KA_CLZ * CLZ_gut - CL_CLZ * C_CLZ / 1000.0;

double C_VAL  = VAL_cen / V_VAL * 1000.0;
double C_NBI  = NBI_cen / V_NBI * 1000.0;
double C_NBIp = NBI_per / V3_NBI * 1000.0;
double cl_nbi = CL_NBI * CYP2D6;
dxdt_VAL_gut = val_dose - KA_VAL * VAL_gut;
dxdt_VAL_cen = KA_VAL * VAL_gut - CL_VAL * C_VAL / 1000.0;
dxdt_NBI_cen = FM_VAL * CL_VAL * C_VAL / 1000.0 - cl_nbi * C_NBI / 1000.0
               - Q_NBI * (C_NBI - C_NBIp) / 1000.0;
dxdt_NBI_per = Q_NBI * (C_NBI - C_NBIp) / 1000.0;

double C_HTB  = HTB_cen / V_HTB * 1000.0;
double C_HTBp = HTB_per / V3_HTB * 1000.0;
double cl_htb = CL_HTB * CYP2D6;
dxdt_DTB_gut = dtb_dose - KA_DTB * DTB_gut;
dxdt_HTB_cen = FM_DTB * KA_DTB * DTB_gut - cl_htb * C_HTB / 1000.0
               - Q_HTB * (C_HTB - C_HTBp) / 1000.0;
dxdt_HTB_per = Q_HTB * (C_HTB - C_HTBp) / 1000.0;

double C_AMA = AMA_cen / V_AMA * 1000.0;
dxdt_AMA_gut = ama_dose - KA_AMA * AMA_gut;
dxdt_AMA_cen = KA_AMA * AMA_gut - CL_AMA * C_AMA / 1000.0;

double C_ACH = ACH_cen / V_ACH * 1000.0;
dxdt_ACH_gut = ach_dose - KA_ACH * ACH_gut;
dxdt_ACH_cen = KA_ACH * ACH_gut - CL_ACH * C_ACH / 1000.0;

dxdt_GKB  = KIN_GKB * gkb_in - KOUT_GKB * GKB;
dxdt_BONT = bont_in / BONT_INT - KOUT_BONT * BONT;

// ---------------- receptor occupancies -----------------------------------
dxdt_CE_AP = KE0_AP * (C_AP - CE_AP);
double occ_ap  = hillf(CE_AP, EC50_D2, HILL_D2);
double occ_clz = hillf(C_CLZ, EC50_D2_CLZ, 1.0);
double OCC     = occ_ap + occ_clz * (1.0 - occ_ap);
if (OCC > 0.995) OCC = 0.995;
double EFF_CLZ = hillf(C_CLZ, EC50_CLZ_EFF, 2.0);
double onbi    = hillf(C_NBI, EC50_VMAT_NBI, 1.0);
double ohtb    = hillf(C_HTB, EC50_VMAT_HTB, 1.0);
double OCCV    = onbi + ohtb * (1.0 - onbi);
if (OCCV > 0.98) OCCV = 0.98;
double E_AMA   = hillf(C_AMA, EC50_AMA, 1.0);
double E_ACH   = hillf(C_ACH, EC50_ACH, 1.0);
double ANTIOX  = ANTIOX_MAX * ((GKB  < 1.0) ? GKB  : 1.0);
double E_BONT  = BONT_MAX   * ((BONT < 1.0) ? BONT : 1.0);

// ---------------- nigrostriatal dopamine pools ---------------------------
double syn  = SYN0 * (1.0 + A_AUTO * OCC);
double fire = 1.0 + A_FIRE * OCC - DEPOL; if (fire < 0.05) fire = 0.05;
double kves = KVES * (1.0 - OCCV);
dxdt_DA_CYT = syn + KBACK * DA_SYN - (kves + KMAO) * DA_CYT;
dxdt_DA_VES = kves * DA_CYT - KREL * fire * DA_VES;
dxdt_DA_SYN = KREL * fire * DA_VES - (KDAT + KDIFF) * DA_SYN;
dxdt_DEPOL  = KIN_DEP * DEPOL_MAX * hillf(OCC, 0.75, 4.0) - KOUT_DEP * DEPOL;

double DA_N  = DA_SYN / DA_SYN0;
double CYT_N = DA_CYT / DA_CYT0;

// ---------------- plasticity: supersensitivity, ROS, structural latch ----
double risk_age = 1.0 + 0.015 * ((AGE > 40.0) ? (AGE - 40.0) : 0.0);
double RISKMOD  = risk_age * GEN_RISK * (1.0 + 0.25 * DM_RISK)
                  * (1.0 - 0.25 * ESTROGEN);
dxdt_RUP = KIN_R * hillf(OCC, OCC50_R, HILL_R) * RISKMOD * (1.0 + W_ACH_R * E_ACH)
           - KOUT_R * RUP / (1.0 + W_SDAM_R * SDAM);

double cyt_term = pow(CYT_N, 1.5) - 1.0; if (cyt_term < 0.0) cyt_term = 0.0;
double ros_drive = W_CYT_ROS * cyt_term
                   + W_OCC_ROS * hillf(OCC, 0.70, 4.0) * RISK_FGA
                   + W_RUP_ROS * RUP
                   + W_S_ROS * hillf(SDAM, S_ROS50, HILL_S_ROS)
                   + 0.25 * ((AGE > 40.0) ? (AGE - 40.0) : 0.0) / 30.0
                   + 0.20 * DM_RISK + ROS_CONST;
dxdt_ROS = KIN_ROS * ros_drive
           - KOUT_ROS * (1.0 + ANTIOX + W_AMA_ROS * E_AMA) * ROS;

double drive_s = W_ROS_S * ROS + W_RUP_S * RUP;
dxdt_SDAM = KIN_S * hillf(drive_s, S50, HILL_S) * (1.0 - SDAM) - KOUT_S * SDAM;

// ---------------- effective postsynaptic D2 stimulation ------------------
double D2STIM = DA_N * (1.0 + G_RUP * RUP) * (1.0 - MASK * OCC)
                * (1.0 + G_ACH * E_ACH) / (1.0 + G_AMA * E_AMA);
double EXC = ((D2STIM > 1.0) ? (D2STIM - 1.0) : 0.0) + G_S * SDAM;

// ---------------- basal-ganglia loop ------------------------------------
double ind_t  = 1.0 / (1.0 + K_IND * EXC);
dxdt_IND  = KBG * (ind_t - IND);
double gpe_t  = 2.0 - IND;
dxdt_GPE  = KBG * (gpe_t - GPE);
double stn_t  = 1.0 / (0.5 + 0.5 * GPE);
dxdt_STN  = KBG * (stn_t - STN);
double gpi_t  = 0.55 * STN + 0.45 / (1.0 + K_DIR * EXC);
dxdt_GPI  = KBG * (gpi_t - GPI);
double thal_t = 1.0 / (1.0 - K_THAL_GPI + K_THAL_GPI * GPI);
dxdt_THAL = KBG * (thal_t - THAL);

double aims_t = AIMS_MAX * hillf(THAL - 1.0, TH50, HILL_TH) * (1.0 - E_BONT);
dxdt_AIMS = K_AIMS * (aims_t - AIMS);

// ---------------- psychiatric, motor, mood, cardiac ---------------------
double control = E_MAX_P * hillf(OCC, OCC50_P, HILL_P)
                 / (1.0 + W_SUPER_P * pow(RUP, 1.5)) + W_CLZ_P * EFF_CLZ;
if (control > 1.0) control = 1.0;
double psy_t = 1.0 - control; if (psy_t < 0.0) psy_t = 0.0;
dxdt_PSYCH = K_PSY * (psy_t - PSYCH);

double park_t = 1.2 * hillf(OCC, OCC50_PK, HILL_PK)
                + 0.9 * hillf(OCCV, VOCC50_PK, HILL_VPK) - W_RUP_PK * RUP;
if (park_t < 0.0) park_t = 0.0;
dxdt_PARK = K_PARK * (park_t - PARK);

double depr_t = W_VMAT_DEP * hillf(OCCV, VOCC50_DEP, HILL_VDEP)
                + W_OCC_DEP * hillf(OCC, 0.80, 4.0);
dxdt_DEPR = K_DEPR * (depr_t - DEPR);

double qtc_t = QT_NBI * C_NBI + QT_HTB * C_HTB;
dxdt_QTC = K_QTC * (qtc_t - QTC);

double adh_t = 1.0 / (1.0 + W_PARK_AD * PARK + W_DEPR_AD * DEPR
                      + W_PSY_AD * PSYCH);
dxdt_ADHER = K_ADHER * (adh_t - ADHER);

double esc = ESC_ON * K_ESC * ((PSYCH > PSY_TARGET) ? (PSYCH - PSY_TARGET) : 0.0);
dxdt_DOSEMULT = esc * (DOSEMULT_MAX - DOSEMULT);

dxdt_CUMOCC = OCC;

$TABLE
double C_APo   = AP_cen  / V2_AP * 1000.0;
double C_NBIo  = NBI_cen / V_NBI * 1000.0;
double C_HTBo  = HTB_cen / V_HTB * 1000.0;
double C_CLZo  = CLZ_cen / V_CLZ * 1000.0;
double occ_apo = hillf(CE_AP, EC50_D2, HILL_D2);
double occ_clzo = hillf(C_CLZo, EC50_D2_CLZ, 1.0);
double OCCo = occ_apo + occ_clzo * (1.0 - occ_apo);
if (OCCo > 0.995) OCCo = 0.995;
double onbio = hillf(C_NBIo, EC50_VMAT_NBI, 1.0);
double ohtbo = hillf(C_HTBo, EC50_VMAT_HTB, 1.0);
double OCCVo = onbio + ohtbo * (1.0 - onbio);
if (OCCVo > 0.98) OCCVo = 0.98;
double DA_No = DA_SYN / DA_SYN0;
double PLAST_DRIVE = hillf(OCCo, OCC50_R, HILL_R);

$CAPTURE @annotated
C_APo : antipsychotic active-moiety concentration (ng/mL)
C_NBIo : NBI-98782 concentration (ng/mL)
C_HTBo : HTBZ concentration (ng/mL)
C_CLZo : clozapine concentration (ng/mL)
OCCo : striatal D2 occupancy (-)
OCCVo : VMAT2 occupancy (-)
DA_No : synaptic dopamine relative to baseline (-)
PLAST_DRIVE : plasticity drive Hill(OCC) (-)
'

mod <- mcode("td_qsp", td_code)

## =============================================================================
## helpers
## =============================================================================
sim <- function(days = 1825, ...) {
  p <- list(...)
  m <- if (length(p)) param(mod, p) else mod
  as.data.frame(mrgsim(m, end = days, delta = 1, hmax = 8,
                       rtol = 1e-6, atol = 1e-9))
}

at_day <- function(d, col, day) {
  approx(d$time, d[[col]], xout = day, rule = 2)$y
}

fmt <- function(x, n = 2) formatC(x, format = "f", digits = n, width = n + 6)

## index patient used throughout: 58 y, risperidone-equivalent 8 mg/day,
## first-generation-level oxidative burden
INDEX <- list(DOSE_AP = 8, DOSE_AP2 = 8, RISK_FGA = 1.6)

## =============================================================================
## 1. treatment scenarios
## =============================================================================
TD_SCENARIOS <- list(
  "S1 natural history (risp 4 mg, 5 y)" =
    list(DOSE_AP = 4, DOSE_AP2 = 4, RISK_FGA = 1.0),
  "S2 FGA-equivalent high occupancy (8 mg)" =
    list(),
  "S3 50% dose reduction at y2" =
    list(DOSE_AP2 = 4, TSW_AP = 730),
  "S4 full withdrawal at y2" =
    list(DOSE_AP2 = 0, TSW_AP = 730),
  "S5 switch to clozapine at y2" =
    list(DOSE_AP2 = 0, TSW_AP = 730, DOSE_CLZ = 350, TSTART_CLZ = 730),
  "S6 valbenazine 80 mg from y2" =
    list(DOSE_VAL = 80, TSTART_VAL = 730),
  "S7 deutetrabenazine 36 mg from y2" =
    list(DOSE_DTB = 36, TSTART_DTB = 730),
  "S8 benztropine add-on from y2" =
    list(DOSE_ACH = 2, TSTART_ACH = 730),
  "S9 ginkgo + amantadine from y2" =
    list(GKB_ON = 1, TSTART_GKB = 730, DOSE_AMA = 200, TSTART_AMA = 730),
  "S10 clinician escalation policy" =
    list(ESC_ON = 1)
)

TD_scenarios <- function(days = 1825) {
  cat("\n", strrep("=", 104), "\n",
      "1. TREATMENT SCENARIOS (index patient; AIMS = dyskinesia subscore 0-26)\n",
      strrep("=", 104), "\n", sep = "")
  cat(sprintf("%-42s%7s%7s%7s%7s%7s%7s%7s%7s%7s\n", "scenario",
              "AIMSy2", "AIMSy5", "RUP", "SDAM", "OCC", "PSY", "PARK",
              "DEPR", "QTc"))
  out <- list()
  for (nm in names(TD_SCENARIOS)) {
    p <- modifyList(INDEX, TD_SCENARIOS[[nm]])
    d <- do.call(sim, c(list(days = days), p))
    out[[nm]] <- d
    cat(sprintf("%-42s%7.2f%7.2f%7.3f%7.3f%7.3f%7.3f%7.3f%7.3f%7.2f\n", nm,
                at_day(d, "AIMS", 730), at_day(d, "AIMS", days),
                at_day(d, "RUP", days), at_day(d, "SDAM", days),
                at_day(d, "OCCo", days), at_day(d, "PSYCH", days),
                at_day(d, "PARK", days), at_day(d, "DEPR", days),
                at_day(d, "QTC", days)))
  }
  invisible(out)
}

## =============================================================================
## 2. is it dose or is it time?  matched cumulative occupancy, split differently
## =============================================================================
TD_exposure_threshold <- function() {
  cat("\n", strrep("=", 104), "\n",
      "2. IS IT DOSE OR IS IT TIME?  same drug, different dose x duration ",
      "splits\n", strrep("=", 104), "\n", sep = "")
  cat(sprintf("%8s%10s%8s%10s%12s%10s%10s%13s%12s\n", "dose", "duration",
              "OCC", "occ-days", "drive-days", "SDAMend", "AIMSend",
              "AIMSstop+6y", "outcome"))
  grid <- list(c(16, 300), c(8, 460), c(6, 560), c(4, 700), c(3, 900),
               c(2, 1500), c(1.5, 1825))
  res <- data.frame()
  for (g in grid) {
    dose <- g[1]; dur <- g[2]
    d <- sim(days = dur + 2190, DOSE_AP = dose, DOSE_AP2 = 0,
             TSW_AP = dur, RISK_FGA = 1.6)
    m <- d[d$time <= dur, ]
    drv <- sum(diff(m$time) * (head(m$PLAST_DRIVE, -1) +
                                 tail(m$PLAST_DRIVE, -1)) / 2)
    a6 <- at_day(d, "AIMS", dur + 2190)
    res <- rbind(res, data.frame(dose = dose, dur = dur,
                                 occ = at_day(d, "OCCo", dur - 5),
                                 occd = at_day(d, "CUMOCC", dur),
                                 drive = drv,
                                 sdam = at_day(d, "SDAM", dur),
                                 aims = at_day(d, "AIMS", dur), aims6 = a6))
    cat(sprintf("%8.1f%10.0f%8.3f%10.1f%12.1f%10.3f%10.2f%13.2f%12s\n",
                dose, dur, tail(res$occ, 1), tail(res$occd, 1), drv,
                tail(res$sdam, 1), tail(res$aims, 1), a6,
                ifelse(a6 > 2, "PERSISTENT", "resolved")))
  }
  invisible(res)
}

## =============================================================================
## 3. reversibility window — when does withdrawal stop working?
## =============================================================================
TD_reversibility_window <- function() {
  cat("\n", strrep("=", 104), "\n",
      "3. REVERSIBILITY WINDOW — withdrawal at increasing exposure duration",
      "\n", strrep("=", 104), "\n", sep = "")
  cat(sprintf("%9s%11s%10s%11s%11s%10s%10s%10s%13s\n", "stop day",
              "SDAM@stop", "RUP@stop", "AIMS@stop", "AIMSpeak", "peak d",
              "AIMS+2y", "AIMS+6y", "outcome"))
  stops <- c(60, 90, 120, 150, 180, 210, 240, 260, 270, 275, 280, 285,
             290, 300, 330, 365, 460, 550, 730, 1095, 1460)
  res <- data.frame()
  for (s in stops) {
    d <- sim(days = s + 2190, DOSE_AP = 8, DOSE_AP2 = 0, TSW_AP = s,
             RISK_FGA = 1.6)
    post <- d[d$time >= s, ]
    ip <- which.max(post$AIMS)
    a6 <- at_day(d, "AIMS", s + 2190)
    res <- rbind(res, data.frame(stop = s, sdam = at_day(d, "SDAM", s),
                                 rup = at_day(d, "RUP", s),
                                 aims = at_day(d, "AIMS", s),
                                 peak = post$AIMS[ip],
                                 peakd = post$time[ip] - s,
                                 a2 = at_day(d, "AIMS", s + 730), a6 = a6))
    r <- tail(res, 1)
    cat(sprintf("%9.0f%11.3f%10.3f%11.2f%11.2f%10.0f%10.2f%10.2f%13s\n",
                r$stop, r$sdam, r$rup, r$aims, r$peak, r$peakd, r$a2, r$a6,
                ifelse(a6 > 2, "PERSISTENT", "resolved")))
  }
  invisible(res)
}

## =============================================================================
## 4. the withdrawal paradox and its crossover time
## =============================================================================
TD_withdrawal_crossover <- function() {
  cat("\n", strrep("=", 104), "\n",
      "4. WITHDRAWAL PARADOX — AIMS after a dose change on day 730\n",
      strrep("=", 104), "\n", sep = "")
  ref <- do.call(sim, c(list(days = 2555), INDEX))
  arms <- list("continue 8 mg" = list(DOSE_AP2 = 8),
               "reduce to 6 mg" = list(DOSE_AP2 = 6),
               "reduce to 4 mg" = list(DOSE_AP2 = 4),
               "reduce to 2 mg" = list(DOSE_AP2 = 2),
               "stop" = list(DOSE_AP2 = 0),
               "stop + clozapine 350" = list(DOSE_AP2 = 0, DOSE_CLZ = 350,
                                             TSTART_CLZ = 730))
  cat(sprintf("%-24s%9s%8s%8s%8s%8s%8s%8s%8s%12s\n", "strategy", "d730",
              "d744", "d760", "d820", "d1095", "d1825", "d2555", "peak",
              "crossover"))
  res <- data.frame()
  for (nm in names(arms)) {
    p <- modifyList(modifyList(INDEX, list(TSW_AP = 730)), arms[[nm]])
    d <- do.call(sim, c(list(days = 2555), p))
    post <- d$time >= 730
    dif <- d$AIMS[post] - ref$AIMS[post]
    cross <- if (any(dif < 0)) d$time[post][which(dif < 0)[1]] - 730 else NA
    vals <- sapply(c(730, 744, 760, 820, 1095, 1825, 2555),
                   function(x) at_day(d, "AIMS", x))
    cat(sprintf("%-24s%9.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%12.0f\n",
                nm, vals[1], vals[2], vals[3], vals[4], vals[5], vals[6],
                vals[7], max(d$AIMS[post]), cross))
    res <- rbind(res, data.frame(arm = nm, t(vals),
                                 peak = max(d$AIMS[post]), crossover = cross))
  }
  invisible(res)
}

## =============================================================================
## 5. VMAT2 inhibitor dose-response, short and long run
## =============================================================================
TD_vmat2_dose_response <- function() {
  cat("\n", strrep("=", 104), "\n",
      "5. VMAT2 INHIBITOR DOSE-RESPONSE (added at y2; week-6 read-out = ",
      "day 772)\n", strrep("=", 104), "\n", sep = "")
  ref <- do.call(sim, c(list(days = 1825), INDEX))
  base <- at_day(ref, "AIMS", 730)
  cat(sprintf("%-26s%8s%10s%8s%8s%8s%7s%7s%7s%7s%7s%9s%9s\n",
              "drug / dose", "Cave", "VMAT2occ", "DA_syn", "AIMS", "dAIMS",
              "%", "PARK", "DEPR", "QTc", "TI", "SDAMy5", "AIMSy5"))
  grid <- c(lapply(c(20, 40, 60, 80, 120),
                   function(x) list(drug = "valbenazine", dose = x,
                                    p = list(DOSE_VAL = x, TSTART_VAL = 730))),
            lapply(c(12, 24, 36, 48, 72),
                   function(x) list(drug = "deutetrabenazine", dose = x,
                                    p = list(DOSE_DTB = x, TSTART_DTB = 730))))
  res <- data.frame()
  for (g in grid) {
    d <- do.call(sim, c(list(days = 1825), modifyList(INDEX, g$p)))
    conc <- at_day(d, if (g$drug == "valbenazine") "C_NBIo" else "C_HTBo", 772)
    a <- at_day(d, "AIMS", 772); da <- a - base
    pk <- at_day(d, "PARK", 772); dp <- at_day(d, "DEPR", 772)
    cat(sprintf("%-26s%8.1f%10.3f%8.3f%8.2f%8.2f%7.1f%7.3f%7.3f%7.2f%7.2f%9.3f%9.2f\n",
                paste(g$drug, g$dose, "mg"), conc, at_day(d, "OCCVo", 772),
                at_day(d, "DA_No", 772), a, da, 100 * da / base, pk, dp,
                at_day(d, "QTC", 772), -da / max(pk + dp, 1e-6),
                at_day(d, "SDAM", 1825), at_day(d, "AIMS", 1825)))
    res <- rbind(res, data.frame(drug = g$drug, dose = g$dose, conc = conc,
                                 aims = a, daims = da, park = pk, depr = dp))
  }
  cat(sprintf("  reference (no VMAT2 inhibitor): AIMS d730 = %.2f, ",
              base),
      sprintf("d772 = %.2f, SDAM y5 = %.3f, AIMS y5 = %.2f\n",
              at_day(ref, "AIMS", 772), at_day(ref, "SDAM", 1825),
              at_day(ref, "AIMS", 1825)), sep = "")
  invisible(res)
}

## =============================================================================
## 6. CYP2D6 phenotype panel
## =============================================================================
TD_cyp2d6_panel <- function() {
  cat("\n", strrep("=", 104), "\n",
      "6. CYP2D6 PHENOTYPE x VMAT2 INHIBITOR (label doses, day 772)\n",
      strrep("=", 104), "\n", sep = "")
  ref <- do.call(sim, c(list(days = 800), INDEX))
  base <- at_day(ref, "AIMS", 730)
  phen <- list("UM (1.6x)" = 1.6, "normal (1.0x)" = 1.0, "IM (0.7x)" = 0.7,
               "PM (0.5x)" = 0.5, "PM + CYP3A4 inhibitor (0.35x)" = 0.35)
  drugs <- list("valbenazine 80" = list(DOSE_VAL = 80, TSTART_VAL = 730),
                "deutetrabenazine 36" = list(DOSE_DTB = 36, TSTART_DTB = 730))
  cat(sprintf("%-30s%-22s%8s%10s%8s%9s%7s%7s%7s\n", "phenotype", "drug",
              "Cave", "VMAT2occ", "AIMS", "dAIMS%", "PARK", "DEPR", "QTc"))
  res <- data.frame()
  for (ph in names(phen)) for (dr in names(drugs)) {
    d <- do.call(sim, c(list(days = 800),
                        modifyList(INDEX, c(drugs[[dr]],
                                            list(CYP2D6 = phen[[ph]])))))
    conc <- at_day(d, if (grepl("valb", dr)) "C_NBIo" else "C_HTBo", 772)
    a <- at_day(d, "AIMS", 772)
    cat(sprintf("%-30s%-22s%8.1f%10.3f%8.2f%9.1f%7.3f%7.3f%7.2f\n", ph, dr,
                conc, at_day(d, "OCCVo", 772), a, 100 * (a - base) / base,
                at_day(d, "PARK", 772), at_day(d, "DEPR", 772),
                at_day(d, "QTC", 772)))
    res <- rbind(res, data.frame(phenotype = ph, drug = dr, conc = conc,
                                 aims = a))
  }
  invisible(res)
}

## =============================================================================
## 7. suppression vs. disease modification (washout test)
## =============================================================================
TD_suppression_vs_modification <- function() {
  cat("\n", strrep("=", 104), "\n",
      "7. SUPPRESSION vs DISEASE MODIFICATION — 2 y of therapy from y2, ",
      "washout at y4\n", strrep("=", 104), "\n", sep = "")
  arms <- list(
    "no treatment" = list(),
    "valbenazine 80 (stop y4)" = list(DOSE_VAL = 80, TSTART_VAL = 730,
                                      TSTOP_VAL = 1460),
    "deutetrabenazine 36 (stop y4)" = list(DOSE_DTB = 36, TSTART_DTB = 730,
                                           TSTOP_DTB = 1460),
    "clozapine switch y2" = list(DOSE_AP2 = 0, TSW_AP = 730, DOSE_CLZ = 350,
                                 TSTART_CLZ = 730),
    "clozapine + valbenazine" = list(DOSE_AP2 = 0, TSW_AP = 730,
                                     DOSE_CLZ = 350, TSTART_CLZ = 730,
                                     DOSE_VAL = 80, TSTART_VAL = 730,
                                     TSTOP_VAL = 1460))
  cat(sprintf("%-32s%9s%9s%8s%9s%12s%9s%14s\n", "arm", "AIMSy2", "AIMSy4",
              "RUPy4", "SDAMy4", "AIMSy4+8w", "AIMSy6", "rebound vs y2"))
  res <- data.frame()
  for (nm in names(arms)) {
    d <- do.call(sim, c(list(days = 2190), modifyList(INDEX, arms[[nm]])))
    v <- c(at_day(d, "AIMS", 730), at_day(d, "AIMS", 1460),
           at_day(d, "RUP", 1460), at_day(d, "SDAM", 1460),
           at_day(d, "AIMS", 1516), at_day(d, "AIMS", 2190))
    cat(sprintf("%-32s%9.2f%9.2f%8.3f%9.3f%12.2f%9.2f%14.2f\n", nm, v[1],
                v[2], v[3], v[4], v[5], v[6], v[5] - v[1]))
    res <- rbind(res, data.frame(arm = nm, t(v)))
  }
  invisible(res)
}

## =============================================================================
## 8. clozapine x valbenazine interaction (2 x 2)
## =============================================================================
TD_combination_interaction <- function() {
  cat("\n", strrep("=", 104), "\n",
      "8. INTERACTION — clozapine switch x valbenazine (effect = AIMS ",
      "change from day 730)\n", strrep("=", 104), "\n", sep = "")
  get <- function(clz, val) {
    p <- INDEX
    if (clz) p <- modifyList(p, list(DOSE_AP2 = 0, TSW_AP = 730,
                                     DOSE_CLZ = 350, TSTART_CLZ = 730))
    if (val) p <- modifyList(p, list(DOSE_VAL = 80, TSTART_VAL = 730))
    d <- do.call(sim, c(list(days = 2555), p))
    sapply(c(730, 772, 1095, 1825, 2555), function(x) at_day(d, "AIMS", x))
  }
  a00 <- get(0, 0); a01 <- get(0, 1); a10 <- get(1, 0); a11 <- get(1, 1)
  cat(sprintf("%-12s%10s%10s%10s%10s%9s%9s%13s\n", "endpoint", "neither",
              "valben", "clozap", "both", "E_val", "E_clz", "interaction"))
  res <- data.frame()
  for (i in 2:5) {
    b <- a00[1]
    e00 <- a00[i] - b; e01 <- a01[i] - b; e10 <- a10[i] - b; e11 <- a11[i] - b
    inter <- e11 - (e01 + e10 - e00)
    cat(sprintf("%-12s%10.2f%10.2f%10.2f%10.2f%9.2f%9.2f%13.2f\n",
                paste("day", c(730, 772, 1095, 1825, 2555)[i]), e00, e01,
                e10, e11, e01 - e00, e10 - e00, inter))
    res <- rbind(res, data.frame(day = c(730, 772, 1095, 1825, 2555)[i],
                                 e00, e01, e10, e11, interaction = inter))
  }
  invisible(res)
}

## =============================================================================
## 9. the opposed-levers frontier
## =============================================================================
TD_opposed_levers <- function() {
  cat("\n", strrep("=", 104), "\n",
      "9. OPPOSED LEVERS — what each way of lowering AIMS costs (day 1095)\n",
      strrep("=", 104), "\n", sep = "")
  strat <- list(
    "continue 8 mg" = list(),
    "lower D2 block: 4 mg" = list(DOSE_AP2 = 4, TSW_AP = 730),
    "lower D2 block: 2 mg" = list(DOSE_AP2 = 2, TSW_AP = 730),
    "lower D2 block: stop" = list(DOSE_AP2 = 0, TSW_AP = 730),
    "clozapine 350 switch" = list(DOSE_AP2 = 0, TSW_AP = 730,
                                  DOSE_CLZ = 350, TSTART_CLZ = 730),
    "lower DA supply: VBZ 40" = list(DOSE_VAL = 40, TSTART_VAL = 730),
    "lower DA supply: VBZ 80" = list(DOSE_VAL = 80, TSTART_VAL = 730),
    "lower DA supply: VBZ 120" = list(DOSE_VAL = 120, TSTART_VAL = 730),
    "both: 4 mg + VBZ 80" = list(DOSE_AP2 = 4, TSW_AP = 730, DOSE_VAL = 80,
                                 TSTART_VAL = 730))
  cat(sprintf("%-30s%7s%7s%8s%8s%7s%7s%7s%7s%7s\n", "strategy from y2",
              "OCC", "VMAT2", "AIMS", "PSYCH", "PARK", "DEPR", "ADHER",
              "RUP", "SDAM"))
  res <- data.frame()
  for (nm in names(strat)) {
    d <- do.call(sim, c(list(days = 1200), modifyList(INDEX, strat[[nm]])))
    v <- sapply(c("OCCo", "OCCVo", "AIMS", "PSYCH", "PARK", "DEPR", "ADHER",
                  "RUP", "SDAM"), function(cc) at_day(d, cc, 1095))
    cat(sprintf("%-30s%7.3f%7.3f%8.2f%8.3f%7.3f%7.3f%7.3f%7.3f%7.3f\n", nm,
                v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]))
    res <- rbind(res, data.frame(strategy = nm, t(v)))
  }
  invisible(res)
}

## =============================================================================
## 10. host risk-factor scan
## =============================================================================
TD_risk_scan <- function() {
  cat("\n", strrep("=", 104), "\n",
      "10. HOST RISK MODIFIERS — identical regimen (8 mg risp-eq, 5 y)\n",
      strrep("=", 104), "\n", sep = "")
  hosts <- list(
    "25 y, SGA" = list(AGE = 25, RISK_FGA = 1.0),
    "45 y, SGA" = list(AGE = 45, RISK_FGA = 1.0),
    "58 y, SGA" = list(AGE = 58, RISK_FGA = 1.0),
    "58 y, FGA (index patient)" = list(AGE = 58, RISK_FGA = 1.6),
    "25 y, FGA" = list(AGE = 25, RISK_FGA = 1.6),
    "75 y, FGA" = list(AGE = 75, RISK_FGA = 1.6),
    "58 y, FGA, diabetes" = list(AGE = 58, RISK_FGA = 1.6, DM_RISK = 1),
    "58 y, FGA, estrogen-replete" = list(AGE = 58, RISK_FGA = 1.6,
                                         ESTROGEN = 1),
    "58 y, FGA, high-risk genotype" = list(AGE = 58, RISK_FGA = 1.6,
                                           GEN_RISK = 1.3),
    "58 y, FGA + benztropine 2 mg" = list(AGE = 58, RISK_FGA = 1.6,
                                          DOSE_ACH = 2, TSTART_ACH = 0),
    "58 y, FGA + ginkgo EGb761" = list(AGE = 58, RISK_FGA = 1.6, GKB_ON = 1,
                                       TSTART_GKB = 0))
  cat(sprintf("%-32s%9s%9s%9s%8s%9s%11s\n", "host", "AIMSy1", "AIMSy3",
              "AIMSy5", "RUPy5", "SDAMy5", "latch day"))
  res <- data.frame()
  for (nm in names(hosts)) {
    d <- do.call(sim, c(list(days = 1825),
                        modifyList(list(DOSE_AP = 8, DOSE_AP2 = 8),
                                   hosts[[nm]])))
    lat <- d$time[which(d$SDAM > 0.49)[1]]
    cat(sprintf("%-32s%9.2f%9.2f%9.2f%8.3f%9.3f%11.0f\n", nm,
                at_day(d, "AIMS", 365), at_day(d, "AIMS", 1095),
                at_day(d, "AIMS", 1825), at_day(d, "RUP", 1825),
                at_day(d, "SDAM", 1825), ifelse(is.na(lat), NA, lat)))
    res <- rbind(res, data.frame(host = nm, latch = lat))
  }
  invisible(res)
}

## =============================================================================
## 11. bistability of the structural latch (drug-free nullcline)
## =============================================================================
TD_latch_bistability <- function() {
  cat("\n", strrep("=", 104), "\n",
      "11. STRUCTURAL LATCH — dSDAM/dt with no external drive (drug-free)\n",
      strrep("=", 104), "\n", sep = "")
  p <- as.list(param(mod))
  hl <- function(x, x50, n) {
    x <- pmax(x, 0); x^n / (x^n + x50^n)
  }
  s <- seq(0, 1, by = 0.05)
  ros <- p$W_S_ROS * hl(s, p$S_ROS50, p$HILL_S_ROS)
  drive <- p$W_ROS_S * ros
  ds <- p$KIN_S * hl(drive, p$S50, p$HILL_S) * (1 - s) - p$KOUT_S * s
  cat(sprintf("%8s%9s%9s%17s%12s\n", "SDAM", "drive", "Hill", "dSDAM/dt",
              "direction"))
  for (i in seq_along(s))
    cat(sprintf("%8.2f%9.3f%9.4f%17.3e%12s\n", s[i], drive[i],
                hl(drive[i], p$S50, p$HILL_S), ds[i],
                ifelse(ds[i] > 0, "up", "down")))
  up <- which(diff(sign(ds)) > 0)
  dn <- which(diff(sign(ds)) < 0)
  if (length(up)) {
    cat(sprintf("\n  lower stable state : SDAM = 0 (repair wins)\n"))
    cat(sprintf("  unstable threshold : SDAM between %.2f and %.2f\n",
                s[up[1]], s[up[1] + 1]))
    dn2 <- dn[dn > up[1]]
    if (length(dn2))
      cat(sprintf("  upper stable state : SDAM between %.2f and %.2f",
                  s[dn2[1]], s[dn2[1] + 1]),
          "-> self-sustaining damage (irreversible TD)\n")
  }
  invisible(data.frame(SDAM = s, dSDAM = ds))
}

## =============================================================================
## run everything
## =============================================================================
TD_run_all <- function() {
  TD_scenarios()
  TD_exposure_threshold()
  TD_reversibility_window()
  TD_withdrawal_crossover()
  TD_vmat2_dose_response()
  TD_cyp2d6_panel()
  TD_suppression_vs_modification()
  TD_combination_interaction()
  TD_opposed_levers()
  TD_risk_scan()
  TD_latch_bistability()
  cat("\n", strrep("=", 104), "\n",
      "All values are model outputs, not literature values. ",
      "See td_references.md for the calibration sources\n",
      "and td_reference_check.py for an independent numerical ",
      "reproduction of every number.\n",
      strrep("=", 104), "\n", sep = "")
  invisible(NULL)
}

if (!interactive() && !exists("TD_SKIP_RUN")) {
  ## Rscript td_mrgsolve_model.R
  ## (set TD_SKIP_RUN <- TRUE before source() to load the model only,
  ##  which is what td_shiny_app.R does)
  TD_run_all()
}
