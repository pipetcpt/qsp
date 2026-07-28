# =============================================================================
#  ppd_mrgsolve_model.R
#  Postpartum depression (PPD) — Quantitative Systems Pharmacology model
#  ---------------------------------------------------------------------------
#  산후우울증 QSP 모델 · mrgsolve · 38 ODEs · 9 treatment scenarios
#
#  ===========================================================================
#  CENTRAL MODELLING THESIS
#  ===========================================================================
#  This model is NOT organised around "allopregnanolone falls after delivery
#  and brexanolone puts it back".  It is organised around a single structural
#  fact about a PRODUCT of two factors with different time constants:
#
#      TONIC INHIBITION IS A PRODUCT OF A LIGAND THAT COLLAPSES IN HOURS AND
#      A RECEPTOR POOL THAT RE-EXPANDS OVER WEEKS.
#
#      G_tonic(t) = R_delta(t) * [ 1 + Emax * PAM(t)^h / (EC50^h + PAM(t)^h) ]
#                   \__________/   \_______________________________________/
#                    receptor arm            ligand arm
#                    t1/2 ~ 120 HOURS        t1/2 ~ 12 HOURS
#
#  During pregnancy the receptor arm is homeostatically down-regulated so that
#  the product stays at its non-pregnant set-point:
#
#      R_delta,ss(PAM) = G_TARGET / (1 + potentiation(PAM))
#
#  This is why EVERY steady state in the model has the same tonic inhibition
#  (Table 1 of the reference output: 1.000 non-pregnant, 1.000 at term, 1.000
#  five days postpartum).  Postpartum depression is therefore NOT a steady
#  state of this system and no steady-state analysis can produce it.  It is
#  the TRANSIENT between two identical steady states.
#
#  ---------------------------------------------------------------------------
#  THESIS 1.  THE ONSET WINDOW IS DERIVED FROM TWO TIME CONSTANTS AND NOTHING
#             ELSE.
#
#  Plasma allopregnanolone falls ~21-fold in 72 h (81.5 -> 3.9 nM in the
#  simulation).  The delta-subunit pool, sitting at its term value of 0.400,
#  has moved only ~18 % of the way back to its non-pregnant value of 0.964 by
#  then.  The product is therefore 0.411 against a set-point of 1.000:
#
#      A 59 % LOSS OF TONIC INHIBITION THAT NEITHER FACTOR SHOWS ALONE.
#
#  Solving for when the deficit recovers gives a window, not a moment:
#      deficit > 25 % .......... for  6.2 days
#      deficit > 10 % .......... for 12.8 days
#      deficit >  5 % .......... for 17.8 days
#  i.e. a 2-3 week window, which is where the epidemiology puts PPD onset
#  (PMID 23487258).  Nothing in the code was told about postpartum weeks.
#  The WIDTH is set by KR, which has no human measurement (A8), so read it
#  as an order of magnitude and not as a prediction.
#
#  ---------------------------------------------------------------------------
#  THESIS 2.  THE SAME HORMONE TRAJECTORY PRODUCES BLUES OR AN EPISODE
#             DEPENDING ON ONE PARAMETER.
#
#  Every simulated woman gets the identical steroid collapse.  Vulnerability V
#  scales symptom accrual, slows delta-pool recovery, and narrows the E/I
#  reserve.  The result is a graded population outcome from ONE parameter:
#      V = 0.80  peak EPDS  2.3   day-42 HAM-D  4.0   (nothing happens)
#      V = 1.00  peak EPDS  7.1   day-42 HAM-D  7.0   (self-limited dip)
#      V = 1.40  peak EPDS 16.1   day-42 HAM-D 23.9   (sustained episode)
#      V = 1.70  peak EPDS 19.7   day-42 HAM-D 27.9   (severe, sustained)
#  The average-V woman stays under the 12/13 EPDS screening cut-off and
#  recovers on her own — the model's account of "baby blues" — while the
#  high-V woman latches into a sustained episode.  This is the
#  formal version of Bloch's 2000 experiment (PMID 10831472), in which
#  identical steroid withdrawal triggered symptoms only in women with a
#  history of postpartum depression.
#
#  ---------------------------------------------------------------------------
#  THESIS 3.  THE FLAT 60 -> 90 ug/kg/h DOSE-RESPONSE IS A PREDICTION OF THE
#             HILL FUNCTION, NOT A TRIAL ANOMALY.
#
#  With EC50 = 120 nM brain allopregnanolone-equivalents and h = 1.4:
#      60 ug/kg/h -> Css 188 nM plasma -> 377 nM brain-eq -> potentiation 2.08
#      90 ug/kg/h -> Css 283 nM plasma -> 565 nM brain-eq -> potentiation 2.24
#  A 50 % increase in exposure buys +7.8 % potentiation.  Site-2 direct
#  channel activation (sedation, PMID 17108970) has no such ceiling in this
#  range.  So the model predicts a flat efficacy comparison with more
#  sedation at the higher rate, which is what phase 3 study 1 reported
#  (-19.5 for BRX60 vs -17.7 for BRX90, PMID 30177236).
#
#  The same arithmetic says the infusion does NOT "restore pregnancy levels":
#  at 90 ug/kg/h it OVERSHOOTS third-trimester allopregnanolone ~3.5-fold.
#  That is why sedation, not efficacy, is the dose-limiting toxicity.
#
#  ---------------------------------------------------------------------------
#  THESIS 4.  DURABILITY BELONGS TO THE RECEPTOR, NOT TO THE DRUG.
#
#  Brexanolone stops at 60 h and zuranolone at day 14, yet both keep working.
#  Freezing receptor plasticity (KR = 0) leaves the ON-DRUG response intact
#  and costs 3-4 HAM-D points at day 45 in BOTH arms, with no parameter of the
#  symptom equations touched.  So part of what looks like drug durability is
#  receptor recovery that would have happened anyway.
#
#  BUT THIS IS ALSO WHERE THE MODEL FAILS, AND THE FAILURE IS INFORMATIVE.
#  The model does NOT reproduce the MAINTENANCE of drug-placebo separation at
#  days 28-45: every arm converges to about HAM-D 14, because that partially-
#  recovered state is the system's late-time attractor once the trough has
#  closed, and no drug parameter moves where the attractor sits.  Two obvious
#  repairs were implemented as switches (KR_BOOST, W_PAM_BDNF) and both make
#  the fit WORSE, for one reason that is worth stating:
#
#      BREXANOLONE *IS* ALLOPREGNANOLONE, SO THE DRUG AND THE ENDOGENOUS
#      LIGAND ENTER THE MODEL THROUGH THE SAME TERM.  ANY PARAMETER THAT
#      GIVES THE DRUG EXTRA CREDIT ALSO GIVES PREGNANCY EXTRA CREDIT.
#
#  Raising W_PAM_BDNF drops the simulated enrolment HAM-D from 28.8 to 23.9,
#  shrinking the room a drug has to work in.  KR_BOOST fails even more
#  directly: in a homeostatic set-point model a PAM cannot accelerate recovery
#  of the pool, because it LOWERS the set-point the pool is chasing.
#  The missing durability mechanism therefore cannot live on the
#  allopregnanolone-potentiation axis at all — it has to be something the drug
#  does that pregnancy does not.
#
#  Consequently the course-length prediction does NOT come out either: the
#  3-day zuranolone course is worse at day 15 (it has stopped) but
#  indistinguishable at day 45.  This model supports only the weak claim that
#  a course must still be running to hold an acute advantage; it cannot say
#  that shorter courses relapse.
#
#  ---------------------------------------------------------------------------
#  THESIS 5.  THE ENORMOUS PLACEBO RESPONSE IN THESE TRIALS IS PARTLY
#             MECHANISTIC, AND THE MODEL SAYS WHICH PART.
#
#  The placebo arms improved by 8.8 to 14.0 HAM-D points, which is larger than
#  the drug-placebo difference in every one of these trials.  Three of the
#  contributors are state variables in this model, not noise:
#    (i)   an enrolled woman is sitting on a delta-pool that is still
#          re-expanding, so she would have improved anyway;
#    (ii)  the brexanolone trials were INPATIENT — 60 h of continuous nursing
#          care protects maternal sleep, which is an active input here;
#    (iii) infant night-waking load itself decays over the first months.
#  Only the residual (expectancy, structured contact) is fitted as a
#  non-mechanistic term.  ONE half of this works and the other does not:
#    SUPPORTED — an inpatient placebo arm beats an outpatient one at every
#    enrolment day, which is the observed ordering between the brexanolone
#    (-14.0 at 60 h, inpatient) and zuranolone (-11.6 at day 15, outpatient)
#    programmes, and the model was not told about it.
#    NOT SUPPORTED — after calibration the predicted placebo response is
#    nearly FLAT in enrolment day (~1 HAM-D point from day 7 to day 180) even
#    though the delta-pool moves from 0.66 to 0.97, because the fitted care
#    term dominates the mechanistic one.  The attractive hypothesis that
#    receptor recovery explains much of the placebo response is NOT borne out
#    at these weights.
#
#  ---------------------------------------------------------------------------
#  THESIS 6.  ONSET KINETICS SEPARATE MECHANISMS THAT WEEK-6 EFFICACY CANNOT.
#
#  Sertraline's SERT occupancy is immediate and its final effect is
#  competitive, yet it cannot move HAM-D on day 3, because its path to the
#  circuit runs through two slow states in series: 5-HT1A autoreceptor
#  desensitisation (t1/2 10 d) and structural plasticity (t1/2 10 d).  A
#  GABA_A PAM bypasses both.  The model reproduces a day-3 difference of more
#  than 10 HAM-D points from time constants alone.  One caveat to name: in
#  this calibration sertraline keeps improving past day 28 and ends up the
#  best arm at day 45, which is an artefact of SYN entering the excitatory
#  load as a divisor with no ceiling — not a claim about comparative efficacy.
#
#  ---------------------------------------------------------------------------
#  THESIS 7.  MATERNAL AND INFANT EXPOSURE ARE DIFFERENT QUESTIONS.
#
#  The model carries a milk compartment and an infant PK compartment, so the
#  relative infant dose is computed rather than asserted, and the infant
#  plasma concentration is compared with the model's OWN potentiation EC50
#  instead of with the maternal dose.  Both neurosteroids land in the low
#  single-digit-percent RID region, and an orally dosed neurosteroid is
#  additionally crippled by first-pass metabolism.
#
#  ===========================================================================
#  WHERE THIS MODEL IS WEAK — READ BEFORE QUOTING ANY NUMBER
#  ===========================================================================
#  A1  Neurosteroids are handled as TOTAL plasma/brain concentrations and the
#      potentiation EC50 is expressed in the same currency, so >99 % protein
#      binding and brain partitioning are LUMPED INTO EC50.  The EC50 here is
#      not comparable to a patch-clamp EC50.
#  A2  Brexanolone IS allopregnanolone, so its partition coefficient and
#      intrinsic potency are set equal to the endogenous values BY IDENTITY.
#      This is the one place the model gets a parameter for free.
#  A3  ZUR_EQ (zuranolone -> allopregnanolone equivalents) is the ONLY
#      parameter calibrated against a zuranolone endpoint.
#  A4  Third-trimester plasma allopregnanolone is taken as ~80 nM and the
#      postpartum floor as ~2 nM.  Assay methods disagree severalfold
#      (PMID 11238543); the RATIO (~40x) is what was matched, and a 2025 IPD
#      meta-analysis (PMID 39511449) finds absolute concentrations do NOT
#      separate depressed from non-depressed women — which is precisely this
#      model's premise.
#  A5  K_CARE / K_NSP are CALIBRATED NON-MECHANISTIC terms for expectancy and
#      structured clinical contact, fitted once to the placebo arms and then
#      held fixed in every active arm.  The model does not explain expectancy.
#  A6  Reference enrolment is day 21 postpartum; real trials enrolled up to
#      6-12 months postpartum, and the predicted placebo response depends
#      strongly on this choice.
#  A7  THE MODEL OPERATES NEAR A BIFURCATION.  The depressed and recovered
#      states are separated by a threshold, and a trajectory passing close to
#      it moves slowly (critical slowing down).  This is deliberate — it is
#      how the model explains why modest interventions can flip outcomes and
#      why placebo trajectories are so gradual and so variable — but it also
#      means late endpoints (day 45) are intrinsically sensitive to THR0 and
#      KR.  The sensitivity table quantifies this rather than hiding it.
#  A8  KR (delta-pool recovery t1/2 ~ 7 d) has NO human measurement.  It is
#      inferred from rodent work (PMID 9789080) and it is the single most
#      influential parameter in the model.  Everything downstream of it should
#      be read as a structural argument, not a prediction with error bars.
#
#  ===========================================================================
#  All numbers in this header were produced by ppd_reference_check.py, a
#  dependency-free Python re-implementation of the SAME equations, because no
#  R runtime was available in the build environment.  Its verbatim output is
#  ppd_reference_output.txt.  If the two disagree, that file was the one run.
#  ===========================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

ppd_code <- '
$PROB
# Postpartum depression QSP model
# 38 ODEs · neurosteroid withdrawal x GABA_A receptor plasticity
# Time is in HOURS.  Delivery occurs at t = TDEL.

$PARAM @annotated
// ---- body and conversions -------------------------------------------------
WT      :  70.0  : maternal weight (kg)
MW_ALLO : 318.5  : allopregnanolone / brexanolone molar mass (g/mol)
MW_ZUR  : 376.5  : zuranolone molar mass (g/mol)
MW_SER  : 306.2  : sertraline molar mass (g/mol)
TDEL    : 1000.0 : simulation time of placental expulsion (h)

// ---- vulnerability --------------------------------------------------------
V       :  1.70  : vulnerability gain (1.0 = population average)

// ---- placental / neurosteroid synthesis ----------------------------------
KDEL    :  4.6    : placental involution rate after expulsion (1/h)
KEL_P4  :  0.1733  : progesterone elimination (1/h, t1/2 4 h)
KSYN_P4 : 68.98    : placental progesterone output at term (nM/h)
IN_P4   :  0.3466  : non-placental progesterone input (nM/h)
K5A     :  0.30    : SRD5A1 flux P4 to 5a-DHP (1/h)
K3A     :  0.90    : AKR1C2 flux 5a-DHP to ALLO (1/h)
KEL_DHP :  0.30    : other 5a-DHP loss (1/h)
F_ALLO  :  0.0514  : fraction of 3a-HSD flux reaching plasma ALLO
KEL_ALLO:  0.0578  : plasma ALLO elimination (1/h, t1/2 12 h)
IN_ALLO :  0.0867  : adrenal/CNS ALLO input (nM/h)
KOUT_B  :  0.35    : brain ALLO equilibration (1/h, t1/2 2 h)
KP_ALLO :  2.0    : brain:plasma partition (identical for brexanolone, A2)
IN_ALLOB:  0.35    : local astrocytic brain ALLO synthesis (nM/h)
KEL_E2  :  0.1155  : estradiol elimination (1/h, t1/2 6 h)
KSYN_E2 :  6.92    : placental estradiol output at term (nM/h)
IN_E2   :  0.00924 : ovarian estradiol floor (nM/h)
KEL_PCRH:  0.693   : placental CRH clearance (1/h)

// ---- HPA axis -------------------------------------------------------------
KH      :  0.000990 : hypothalamic CRH recovery (1/h, t1/2 700 h)
PC50    :  0.55     : placental CRH suppressing hypothalamic CRH 50 percent
KA_ACTH :  1.0     : ACTH turnover (1/h)
W_PCRH  :  1.60     : placental CRH contribution to ACTH drive
CORT_FB : 12.0      : cortisol for half-maximal GR feedback (ug/dL)
N_FB    :  1.5     : Hill coefficient of cortisol feedback
KC      :  6.73     : cortisol production per unit ACTH ((ug/dL)/h)
KEL_CORT:  0.462    : cortisol elimination (1/h, t1/2 1.5 h)
KG      :  0.00963  : GR sensitivity adaptation (1/h, t1/2 72 h)
GR50    : 25.0      : cortisol for half GR down-regulation (ug/dL)
CORT_NP : 10.0      : non-pregnant reference cortisol (ug/dL)

// ---- GABA_A receptor plasticity (the slow arm) ---------------------------
G_TARGET:  1.0     : homeostatic tonic-conductance set-point
EMAX_PAM:  2.5     : maximal fractional potentiation (extrasynaptic)
EC50_PAM: 120.0     : nM brain ALLO-equivalents (lumps binding, see A1)
H_PAM   :  1.4     : Hill coefficient of potentiation
EMAX_SYN:  0.60     : maximal potentiation of synaptic gamma2 receptors
EC50_SYN: 400.0     : nM, synaptic receptors are less steroid-sensitive
KR      :  0.0057566 : delta-pool plasticity (1/h, t1/2 168 h = 7 d)
V_KR    :  0.2808     : vulnerability slowing of receptor recovery
W_TONIC :  0.90     : weight of tonic inhibition in effective inhibition
W_PHASIC:  0.10     : weight of phasic inhibition
RG_GAIN :  0.70     : reciprocal up-regulation of the synaptic pool
KK      :  0.0289   : KCC2 turnover (1/h, t1/2 24 h)
KCC2_CORT: 0.60     : cortisol sensitivity of KCC2 loss
PAM_NP  :  6.0     : nM brain ALLO-equivalents, non-pregnant reference

// ---- monoamine / plasticity ---------------------------------------------
KM_MAOA :  0.01444  : MAO-A adaptation (1/h, t1/2 48 h)
F_MAOA  :  0.21     : maximal MAO-A rise on complete E2 withdrawal
E250    :  0.50     : estradiol restraining MAO-A half-maximally (nM)
K5HT    :  0.0578   : serotonergic tone turnover (1/h)
EMAX_SSRI: 1.20     : maximal 5-HT gain from full SERT blockade
EC50_SER:  9.0     : sertraline for 50 percent SERT occupancy (ng/mL)
W_KYN_5HT: 0.30     : tryptophan-diversion penalty on 5-HT synthesis
KAUTO   :  0.00289  : 5-HT1A autoreceptor desensitisation (1/h, t1/2 240 h)
AUTO_BRAKE:0.85     : fraction of SSRI effect gated by the autoreceptor
OCC50_AUTO:0.50     : SERT occupancy desensitising the autoreceptor 50 percent
KB_BDNF :  0.00578  : BDNF turnover (1/h, t1/2 120 h)
KSYN_PL :  0.00289  : structural plasticity turnover (1/h, t1/2 240 h)
W_AMPA  :  0.80     : AMPA surge to BDNF gain (ketamine arm)
W_PAM_BDNF: 0.0     : does the PAM have a plasticity arm? (tested, rejected)

// ---- inflammation / kynurenine ------------------------------------------
KI_INFL :  0.0289   : IL-6 turnover (1/h, t1/2 24 h)
INFL_BASE: 2.0     : baseline IL-6-equivalent (pg/mL)
W_SD_INFL: 0.50     : sleep debt to inflammation
KK_KYN  :  0.01444  : KYN/TRP turnover (1/h, t1/2 48 h)
W_INFL_KYN:1.00     : inflammation to IDO1 to KYN/TRP
W_CORT_KYN:0.40     : cortisol to TDO to KYN/TRP

// ---- sleep ---------------------------------------------------------------
SLP_MAX :  8.0     : achievable sleep (h/night)
SLP_NEED:  7.5     : sleep requirement (h/night)
A_WAKE  :  0.18     : infant-waking penalty on sleep
A_SYMP  :  1.482    : insomnia symptom penalty on sleep (closes the loop)
KS_SLP  :  0.10     : sleep-state adaptation (1/h)
KDEC_SD :  0.010    : sleep-debt repayment (1/h, t1/2 69 h)
DEL_INFL: 22.0      : delivery inflammatory surge, IL-6-eq bolus (pg/mL)
WAKE_LP :  0.30     : third-trimester nocturnal disturbance
WAKE_PP0:  1.80     : amplitude of newborn night-waking load
WAKE_FL :  1.00     : asymptotic night-waking load
TAU_WAKE: 1080.0    : consolidation of infant sleep (h)
WAKE_PROT: 1.00     : night-waking load under protected sleep

// ---- excitation / symptom transfer --------------------------------------
KE_EXC  :  0.0578   : E/I index adaptation (1/h, t1/2 12 h)
W_SD    :  2.0     : sleep debt to excitatory load
W_INFL  :  0.50     : inflammation to excitatory load
W_KYN   :  0.30     : QUIN/KYNA shift to excitatory load
W_5HT   :  0.60     : low serotonergic tone to excitatory load
W_HPA   :  0.30     : blunted-HPA contribution
W_SELF  :  0.22651    : symptom self-reinforcement (rumination, DMN, glutamate)
W_SYN   :  0.80     : structural plasticity buffering excitatory load
THR0    :  2.634     : E/I threshold for symptom accrual at V = 1
V_THR   :  1.00     : exponent by which vulnerability narrows the reserve
KON     :  0.0036432  : symptom accrual rate (1/h)
KOFF    :  0.0017154 : intrinsic symptom resolution (1/h, scaled by 1/V)
KFAST   : 84.83      : acceleration of resolution below threshold
SMAX    :  1.0     : maximum latent symptom load
SMIN    :  0.155    : residual symptom floor (HAM-D ~7, remission bound)
KBOND   :  0.00963  : bonding impairment turnover (1/h)
W_BOND  :  0.85     : symptom to bonding impairment

// ---- endpoint scales -----------------------------------------------------
HAMD_0  :  2.0     : HAM-D17 intercept
HAMD_SC : 32.0      : HAM-D17 span
EPDS_0  :  1.0     : EPDS intercept
EPDS_SC : 22.0      : EPDS span

// ---- nonspecific care (A5) ----------------------------------------------
K_CARE  :  0.019967  : inpatient continuous-care effect (1/h)
K_NSP   :  0.0058718 : outpatient trial-contact effect (1/h)
K_CBT   :  0.0030   : psychotherapy effect at full engagement (1/h)
KCBT_ON :  0.00248  : psychotherapy engagement build-up (1/h)

// ---- drug PK -------------------------------------------------------------
BRX_CL  : 70.0   : brexanolone clearance (L/h, 1.0 L/h/kg)
BRX_V1  : 105.0  : brexanolone central volume (L)
BRX_V2  : 300.0  : brexanolone peripheral volume (L)
BRX_Q   : 40.0   : brexanolone intercompartmental clearance (L/h)
ZUR_KA  :  0.45  : zuranolone absorption (1/h)
ZUR_CL  : 12.0   : zuranolone CL/F (L/h)
ZUR_V1  : 200.0  : zuranolone V1/F (L)
ZUR_V2  : 100.0  : zuranolone V2/F (L)
ZUR_Q   : 15.0   : zuranolone Q/F (L/h)
ZUR_EQ  :  0.10464  : brain ALLO-eq per nM plasma zuranolone (A3)
SER_KA  :  0.60  : sertraline absorption (1/h)
SER_CL  : 96.0   : sertraline CL/F (L/h)
SER_V   : 3000.0 : sertraline V/F (L)
ESK_CL  : 90.0   : esketamine clearance (L/h)
ESK_V   : 250.0  : esketamine volume (L)
K_AMPA_IN : 0.020   : AMPA surge input per 100 ng/mL esketamine (1/h)
K_AMPA_OUT: 0.01925 : AMPA surge decay (1/h, t1/2 36 h)

// ---- lactation transfer --------------------------------------------------
MP_BRX  :  1.50  : milk:plasma ratio, brexanolone
MP_ZUR  :  1.50  : milk:plasma ratio, zuranolone
MP_SER  :  1.80  : milk:plasma ratio, sertraline
MILK_INTAKE: 150.0 : infant milk intake (mL/kg/day)
F_ORAL_NS: 0.05  : infant oral bioavailability of a neurosteroid
V_INF   :  2.0  : infant volume of distribution (L/kg)
CL_INF  :  0.25  : infant clearance (L/h/kg)

// ---- treatment-window switches (times are HOURS of simulation time) -----
CARE_T0 : 1e9 : start of inpatient continuous care
CARE_T1 : 1e9 : end of inpatient continuous care
PROT_T0 : 1e9 : start of protected-sleep conditions
PROT_T1 : 1e9 : end of protected-sleep conditions
NSP_T0  : 1e9 : start of outpatient trial contact
CBT_T0  : 1e9 : start of psychotherapy

$CMT @annotated
PLAC   : placental functional mass (rel, 1 = term)
P4     : plasma progesterone (nM)
DHP    : plasma 5a-dihydroprogesterone (nM)
ALLOP  : plasma allopregnanolone (nM)
ALLOB  : brain allopregnanolone, effector species (nM)
E2     : plasma estradiol (nM)
PCRH   : placental CRH (rel)
HCRH   : hypothalamic CRH drive (rel)
ACTH   : ACTH (rel)
CORT   : cortisol (ug/dL)
GRFN   : GR feedback sensitivity (rel)
RD     : extrasynaptic delta-GABA_A receptor pool (rel)
RG     : synaptic gamma2-GABA_A receptor pool (rel)
KCC2   : KCC2 surface expression (rel)
MAOA   : MAO-A binding density (rel)
FIVEHT : serotonergic tone (rel)
AUTO   : 5-HT1A autoreceptor sensitivity (rel)
BDNF   : BDNF signalling (rel)
SYN    : structural plasticity index (rel)
INFL   : inflammation, IL-6-equivalent (pg/mL)
KYNR   : kynurenine/tryptophan ratio (rel)
SLP    : sleep obtained (h/night)
SDEBT  : cumulative sleep debt (h)
EXC    : excitation/inhibition imbalance index (rel)
SYMP   : latent depressive symptom load (0-1)
BONDI  : mother-infant bonding impairment (0-1)
BRX1   : brexanolone central (ug)
BRX2   : brexanolone peripheral (ug)
ZURA   : zuranolone absorption depot (mg)
ZUR1   : zuranolone central (mg)
ZUR2   : zuranolone peripheral (mg)
SERA   : sertraline absorption depot (mg)
SERC   : sertraline central (mg)
ESKC   : esketamine central (mg)
AMPAS  : AMPA-mediated plasticity surge (rel)
CBTP   : psychotherapy progress (0-1)
MILKD  : cumulative drug delivered to infant via milk (ug/kg)
INFP   : infant plasma concentration, ALLO-equivalents (nM)

$GLOBAL
#define POT_T(pam)  (EMAX_PAM * pow(fmax(pam,0.0), H_PAM) \\
                     / (pow(EC50_PAM, H_PAM) + pow(fmax(pam,0.0), H_PAM)))
#define POT_S(pam)  (EMAX_SYN * pow(fmax(pam,0.0), H_PAM) \\
                     / (pow(EC50_SYN, H_PAM) + pow(fmax(pam,0.0), H_PAM)))
#define RDSS(pam)   (fmin(1.0, G_TARGET / (1.0 + POT_T(pam))))

$MAIN
// initial conditions = the PREGNANCY STEADY STATE.  Running the model from
// t = 0 with TDEL = 1000 h reproduces these values to 3 decimals, which is
// the check that the analytic and simulated pregnancy states agree.
PLAC_0   = 1.0;
P4_0     = 400.0;
DHP_0    = 100.0;
ALLOP_0  = 81.5;
ALLOB_0  = 164.0;
E2_0     = 60.0;
PCRH_0   = 1.0;
HCRH_0   = 0.348;
ACTH_0   = 1.10;
CORT_0   = 16.3;
GRFN_0   = 0.61;
RD_0     = 0.397;
RG_0     = 1.42;
KCC2_0   = 0.90;
MAOA_0   = 1.00;
FIVEHT_0 = 0.93;
AUTO_0   = 1.00;
BDNF_0   = 0.73;
SYN_0    = 0.73;
INFL_0   = 2.00;
KYNR_0   = 1.25;
SLP_0    = 7.57;
SDEBT_0  = 0.0;
EXC_0    = 1.72;
SYMP_0   = 0.0;
BONDI_0  = 0.0;
KCC2_0   = 0.90;

// non-pregnant reference effective inhibition (normalisation constant)
double rd_np = RDSS(PAM_NP);
double rg_np = 1.0 + RG_GAIN * (1.0 - rd_np);
GREF = W_TONIC * rd_np * (1.0 + POT_T(PAM_NP))
     + W_PHASIC * rg_np * (1.0 + POT_S(PAM_NP));

$ODE
double tad = SOLVERTIME - TDEL;          // hours since delivery
double post = (tad >= 0.0) ? 1.0 : 0.0;

// ---- 1. placenta and steroid trajectories -------------------------------
double plac = fmax(PLAC, 0.0);
dxdt_PLAC  = -KDEL * plac * post;
dxdt_P4    = KSYN_P4 * plac + IN_P4 - KEL_P4 * P4;
dxdt_DHP   = K5A * P4 - (K3A + KEL_DHP) * DHP;
dxdt_ALLOP = F_ALLO * K3A * DHP + IN_ALLO - KEL_ALLO * ALLOP;
dxdt_ALLOB = KOUT_B * KP_ALLO * ALLOP + IN_ALLOB - KOUT_B * ALLOB;
dxdt_E2    = KSYN_E2 * plac + IN_E2 - KEL_E2 * E2;
dxdt_PCRH  = KEL_PCRH * (plac - PCRH);

// ---- 2. drug PK ---------------------------------------------------------
dxdt_BRX1 = -(BRX_CL / BRX_V1) * BRX1 - (BRX_Q / BRX_V1) * BRX1
            + (BRX_Q / BRX_V2) * BRX2;
dxdt_BRX2 =  (BRX_Q / BRX_V1) * BRX1 - (BRX_Q / BRX_V2) * BRX2;
dxdt_ZURA = -ZUR_KA * ZURA;
dxdt_ZUR1 =  ZUR_KA * ZURA - (ZUR_CL / ZUR_V1) * ZUR1
            - (ZUR_Q / ZUR_V1) * ZUR1 + (ZUR_Q / ZUR_V2) * ZUR2;
dxdt_ZUR2 =  (ZUR_Q / ZUR_V1) * ZUR1 - (ZUR_Q / ZUR_V2) * ZUR2;
dxdt_SERA = -SER_KA * SERA;
dxdt_SERC =  SER_KA * SERA - (SER_CL / SER_V) * SERC;
dxdt_ESKC = -(ESK_CL / ESK_V) * ESKC;

double c_brx_ng = BRX1 / BRX_V1;                    // ng/mL
double c_brx_nM = c_brx_ng * 1000.0 / MW_ALLO;
double c_zur_ng = ZUR1 * 1000.0 / ZUR_V1;           // mg/L -> ng/mL
double c_zur_nM = c_zur_ng * 1000.0 / MW_ZUR;
double c_ser_ng = SERC * 1000.0 / SER_V;            // ng/mL
double c_esk_ng = ESKC * 1000.0 / ESK_V;            // ng/mL

// ---- 3. total PAM in brain allopregnanolone equivalents ----------------
// brexanolone IS allopregnanolone: same Kp, potency exactly 1 (assumption A2)
double pam = ALLOB + KP_ALLO * c_brx_nM + ZUR_EQ * KP_ALLO * c_zur_nM;

// ---- 4. HPA axis -------------------------------------------------------
dxdt_HCRH = KH * (1.0 / (1.0 + PCRH / PC50) - HCRH);
double fb = 1.0 / (1.0 + pow(GRFN * CORT / CORT_FB, N_FB));
dxdt_ACTH = KA_ACTH * ((HCRH + W_PCRH * PCRH) * fb - ACTH);
dxdt_CORT = KC * ACTH - KEL_CORT * CORT;
double grss = (1.0 / (1.0 + CORT / GR50)) / (1.0 + 0.40 * (V - 1.0));
dxdt_GRFN = KG * (grss - GRFN);
double cort_rel = CORT / CORT_NP;

// ---- 5. receptor plasticity (THE SLOW ARM) ----------------------------
double rdss  = RDSS(pam);
double kr_eff = KR / (1.0 + V_KR * (V - 1.0))
              * (1.0 + KR_BOOST * POT_T(pam));   // KR_BOOST = 0 by default
dxdt_RD = kr_eff * (rdss - RD);
double rgss = 1.0 + RG_GAIN * (1.0 - rdss);
dxdt_RG = kr_eff * (rgss - RG);
double kcc2ss = (1.0 / (1.0 + KCC2_CORT * pow(fmax(0.0, cort_rel - 1.0), 1.2)))
              * (0.70 + 0.30 * BDNF);
dxdt_KCC2 = KK * (kcc2ss - KCC2);

double g_tonic  = RD * (1.0 + POT_T(pam)) * sqrt(fmax(KCC2, 1e-6));
double g_phasic = RG * (1.0 + POT_S(pam));
double g_eff    = W_TONIC * g_tonic + W_PHASIC * g_phasic;

// ---- 6. monoamine / plasticity ----------------------------------------
double maoass = 1.0 + F_MAOA * (1.0 - E2 / (E2 + E250));
dxdt_MAOA = KM_MAOA * (maoass - MAOA);
double occ = c_ser_ng / (c_ser_ng + EC50_SER);
double ssri_gain = 1.0 + EMAX_SSRI * occ * (1.0 - AUTO_BRAKE * AUTO);
double fivess = (1.0 / fmax(MAOA, 1e-6))
              / (1.0 + W_KYN_5HT * fmax(0.0, KYNR - 1.0)) * ssri_gain;
dxdt_FIVEHT = K5HT * (fivess - FIVEHT);
double autoss = 1.0 / (1.0 + pow(occ / OCC50_AUTO, 2.0));
dxdt_AUTO = KAUTO * (autoss - AUTO);
dxdt_AMPAS = K_AMPA_IN * c_esk_ng / 100.0 - K_AMPA_OUT * AMPAS;
double bdnfss = pow(fmax(FIVEHT, 1e-6), 0.60) * (1.0 + W_AMPA * AMPAS)
              * (1.0 + W_PAM_BDNF * POT_T(pam))   // W_PAM_BDNF = 0 by default
              / (1.0 + 0.50 * fmax(0.0, cort_rel - 1.0))
              / (1.0 + 0.30 * fmax(0.0, INFL - INFL_BASE) / 10.0);
dxdt_BDNF = KB_BDNF * (bdnfss - BDNF);
dxdt_SYN  = KSYN_PL * (BDNF - SYN);

// ---- 7. inflammation / kynurenine -------------------------------------
double inflss = INFL_BASE * (1.0 + W_SD_INFL * SDEBT / 15.0);
dxdt_INFL = KI_INFL * (inflss - INFL);
double kynss = 1.0 + W_INFL_KYN * fmax(0.0, INFL - INFL_BASE) / 10.0
             + W_CORT_KYN * fmax(0.0, cort_rel - 1.0);
dxdt_KYNR = KK_KYN * (kynss - KYNR);

// ---- 8. sleep ---------------------------------------------------------
double protect = (SOLVERTIME >= PROT_T0 && SOLVERTIME < PROT_T1) ? 1.0 : 0.0;
double wk;
if (tad < 0.0)        wk = WAKE_LP;
else if (protect > 0) wk = WAKE_PROT;
else                  wk = WAKE_FL + WAKE_PP0 * exp(-tad / TAU_WAKE);
double slpss = SLP_MAX / (1.0 + A_WAKE * wk + A_SYMP * SYMP);
dxdt_SLP   = KS_SLP * (slpss - SLP);
dxdt_SDEBT = (SLP_NEED - SLP) / 24.0 - KDEC_SD * SDEBT;

// ---- 9. excitatory load and symptoms ---------------------------------
double drive = 1.0
             + W_SD   * SDEBT / 48.0
             + W_INFL * fmax(0.0, INFL - INFL_BASE) / 10.0
             + W_KYN  * fmax(0.0, KYNR - 1.0)
             + W_5HT  * fmax(0.0, 1.0 / fmax(FIVEHT, 1e-6) - 1.0)
             + W_HPA  * fmax(0.0, 1.0 - cort_rel)
             + W_SELF * SYMP;
drive /= pow(fmax(SYN, 1e-6), W_SYN);
dxdt_EXC = KE_EXC * (drive * GREF / fmax(g_eff, 1e-6) - EXC);

double thr  = THR0 / pow(V, V_THR);
double over = fmax(0.0, EXC - thr);
double fsat = over / (1.0 + over);
double koff = (KOFF / V) * (1.0 + KFAST * fmax(0.0, thr - EXC));
double k_care = (SOLVERTIME >= CARE_T0 && SOLVERTIME < CARE_T1) ? K_CARE : 0.0;
double k_nsp  = (SOLVERTIME >= NSP_T0) ? K_NSP : 0.0;
dxdt_CBTP = (SOLVERTIME >= CBT_T0) ? KCBT_ON * (1.0 - CBTP) : 0.0;
double k_cbt = K_CBT * CBTP;
dxdt_SYMP = KON * V * fsat * (SMAX - SYMP)
          - (koff + k_care + k_nsp + k_cbt) * fmax(0.0, SYMP - SMIN);
dxdt_BONDI = KBOND * (W_BOND * SYMP - BONDI);

// ---- 10. lactation transfer ------------------------------------------
double c_milk = MP_BRX * c_brx_ng + MP_ZUR * c_zur_ng + MP_SER * c_ser_ng;
double inf_rate = c_milk * MILK_INTAKE / 24.0 / 1000.0;      // ug/kg/h
dxdt_MILKD = inf_rate;
dxdt_INFP  = F_ORAL_NS * inf_rate / V_INF * 1000.0 / MW_ALLO
           - (CL_INF / V_INF) * INFP;

$TABLE
double TAD   = SOLVERTIME - TDEL;
double DAYPP = TAD / 24.0;
double PAM   = ALLOB + KP_ALLO * (BRX1 / BRX_V1) * 1000.0 / MW_ALLO
             + ZUR_EQ * KP_ALLO * (ZUR1 * 1000.0 / ZUR_V1) * 1000.0 / MW_ZUR;
double POTENT  = POT_T(PAM);
double GTONIC  = RD * (1.0 + POTENT) * sqrt(fmax(KCC2, 1e-6));
double GPHASIC = RG * (1.0 + POT_S(PAM));
double GEFF    = W_TONIC * GTONIC + W_PHASIC * GPHASIC;
double RDSET   = RDSS(PAM);
double DEFICIT = 100.0 * (1.0 - GTONIC / G_TARGET);   // percent tonic deficit
double THRESH  = THR0 / pow(V, V_THR);
double RESERVE = THRESH - EXC;                        // > 0 means sub-threshold
double HAMD    = HAMD_0 + HAMD_SC * SYMP;
double EPDS    = EPDS_0 + EPDS_SC * SYMP;
double REMIT   = (HAMD <= 7.0) ? 1.0 : 0.0;
double CBRX    = BRX1 / BRX_V1;                       // ng/mL
double CZUR    = ZUR1 * 1000.0 / ZUR_V1;              // ng/mL
double CSER    = SERC * 1000.0 / SER_V;               // ng/mL
double SEROCC  = 100.0 * CSER / (CSER + EC50_SER);    // percent SERT occupancy
double CMILK   = MP_BRX * CBRX + MP_ZUR * CZUR + MP_SER * CSER;
// site-2 (direct gating) sedation index: no saturation in the clinical range
double SEDIX   = (PAM > 0.0) ? PAM / 1500.0 : 0.0;

$CAPTURE @annotated
TAD     : hours since delivery
DAYPP   : days postpartum
PAM     : total PAM, brain ALLO-equivalents (nM)
POTENT  : fractional potentiation of tonic conductance
GTONIC  : tonic inhibitory conductance (rel to set-point)
GPHASIC : phasic inhibition (rel)
GEFF    : effective inhibition (weighted)
RDSET   : delta-pool homeostatic set-point at current PAM
DEFICIT : tonic-inhibition deficit (percent of set-point)
THRESH  : E/I threshold for symptom accrual
RESERVE : threshold minus EXC (positive = sub-threshold)
HAMD    : HAM-D17 total score
EPDS    : Edinburgh Postnatal Depression Scale
REMIT   : remission indicator (HAM-D <= 7)
CBRX    : brexanolone plasma concentration (ng/mL)
CZUR    : zuranolone plasma concentration (ng/mL)
CSER    : sertraline plasma concentration (ng/mL)
SEROCC  : SERT occupancy (percent)
CMILK   : total drug concentration in milk (ng/mL)
SEDIX   : site-2 sedation index (rel)
GREF    : non-pregnant reference effective inhibition
'

mod <- mcode("ppd", ppd_code)

# =============================================================================
#  SIMULATION HELPERS
# =============================================================================
TDEL   <- 1000                       # simulation time of delivery (h)
ENROL  <- TDEL + 21 * 24             # reference enrolment: day 21 postpartum
TEND   <- ENROL + 46 * 24

#' Delivery event: the one discrete consequence that is not a rate — the
#' sterile inflammatory surge of parturition (IL-6-equivalent bolus).
ev_delivery <- function() ev(time = TDEL, amt = 22, cmt = "INFL")

#' Brexanolone titration exactly as labelled: 30 -> 60 -> (60 or 90) -> 30
#' ug/kg/h over 4 + 20 + 28 + 8 hours, 60 hours total.
ev_brexanolone <- function(top = 90, wt = 70, start = ENROL) {
  rates <- c(30, min(60, top), top, 30)
  durs  <- c(4, 20, 28, 8)
  t0 <- start
  out <- NULL
  for (i in seq_along(rates)) {
    e <- ev(time = t0, amt = rates[i] * wt * durs[i],
            rate = rates[i] * wt, cmt = "BRX1")
    out <- if (is.null(out)) e else c(out, e)
    t0 <- t0 + durs[i]
  }
  out
}

ev_zuranolone <- function(dose = 50, days = 14, start = ENROL)
  ev(time = start, amt = dose, cmt = "ZURA", ii = 24, addl = days - 1)

ev_sertraline <- function(start = ENROL, days = 46)
  c(ev(time = start, amt = 50, cmt = "SERA", ii = 24, addl = 6),
    ev(time = start + 7 * 24, amt = 100, cmt = "SERA", ii = 24,
       addl = days - 8))

ev_esketamine <- function(dose_mgkg = 0.25, wt = 70, start = ENROL)
  ev(time = start, amt = dose_mgkg * wt, rate = dose_mgkg * wt / (2/3),
     cmt = "ESKC")

#' Window switches for the non-pharmacological / trial-setting inputs.
win <- function(care = NULL, protect = NULL, nsp = NULL, cbt = NULL) {
  p <- list(CARE_T0 = 1e9, CARE_T1 = 1e9, PROT_T0 = 1e9, PROT_T1 = 1e9,
            NSP_T0 = 1e9, CBT_T0 = 1e9)
  if (!is.null(care))    { p$CARE_T0 <- care[1];    p$CARE_T1 <- care[2] }
  if (!is.null(protect)) { p$PROT_T0 <- protect[1]; p$PROT_T1 <- protect[2] }
  if (!is.null(nsp))       p$NSP_T0  <- nsp
  if (!is.null(cbt))       p$CBT_T0  <- cbt
  p
}

run_arm <- function(events = NULL, params = list(), V = 1.70, end = TEND,
                    delta = 1) {
  e <- ev_delivery()
  if (!is.null(events)) e <- c(e, events)
  mod %>%
    param(c(list(V = V), params)) %>%
    ev(e) %>%
    mrgsim(end = end, delta = delta) %>%
    as_tibble()
}

# =============================================================================
#  SCENARIO 0 — NATURAL HISTORY: the trough, and what V does to it
#  This is the scenario that carries the argument.  No treatment at all.
# =============================================================================
sc0 <- bind_rows(lapply(
  c(0.80, 1.00, 1.40, 1.70, 2.00),
  function(v) run_arm(V = v) %>% mutate(V = v)
))

# The whole thesis in one figure: tonic inhibition troughs because a
# fast-falling ligand multiplies a slow-recovering receptor pool.
p_trough <- sc0 %>%
  filter(V == 1.70, DAYPP >= -2, DAYPP <= 45) %>%
  select(DAYPP, `plasma ALLO (nM)` = ALLOP, `delta pool R` = RD,
         `tonic G` = GTONIC, `HAM-D17` = HAMD) %>%
  pivot_longer(-DAYPP) %>%
  ggplot(aes(DAYPP, value)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~name, scales = "free_y") +
  labs(x = "days postpartum", y = NULL,
       title = "Fast ligand x slow receptor = a transient trough",
       subtitle = "no treatment; V = 1.70") +
  theme_bw()

# =============================================================================
#  SCENARIO 1 — PLACEBO, INPATIENT (the brexanolone-trial setting)
#  60 h of continuous nursing care.  Protected sleep is an ACTIVE input here,
#  which is the model's explanation of the -14.0 HAM-D placebo response.
# =============================================================================
sc1 <- run_arm(params = win(care = c(ENROL, ENROL + 60),
                            protect = c(ENROL, ENROL + 60),
                            nsp = ENROL))

# =============================================================================
#  SCENARIO 2 — BREXANOLONE 60 and 90 ug/kg/h (60-hour infusion)
#  Observed: 60 h HAM-D -19.5 (BRX60) · -17.7 (BRX90) · -14.0 (placebo)
#  (Meltzer-Brody 2018 Lancet, PMID 30177236, study 1)
# =============================================================================
sc2_60 <- run_arm(ev_brexanolone(top = 60),
                  win(care = c(ENROL, ENROL + 60),
                      protect = c(ENROL, ENROL + 60), nsp = ENROL))
sc2_90 <- run_arm(ev_brexanolone(top = 90),
                  win(care = c(ENROL, ENROL + 60),
                      protect = c(ENROL, ENROL + 60), nsp = ENROL))

# =============================================================================
#  SCENARIO 3 — ZURANOLONE 50 mg x 14 d (SKYLARK) and 30 mg x 14 d (ROBIN)
#  Observed day 15: -15.6 vs -11.6 (SKYLARK, PMID 37491938)
#                   -17.8 vs -13.6 (ROBIN,  PMID 34190962)
#  The model predicts 30 and 50 mg to be nearly indistinguishable because
#  potentiation is saturated — which is what the two trials' drug-placebo
#  differences (-4.2 and -4.0) actually show.
# =============================================================================
sc3_50 <- run_arm(ev_zuranolone(50, 14), win(nsp = ENROL))
sc3_30 <- run_arm(ev_zuranolone(30, 14), win(nsp = ENROL))

# =============================================================================
#  SCENARIO 4 — COURSE LENGTH: 3-day versus 14-day zuranolone
#  Same daily dose.  RESULT: the short course is worse at day 15 and
#  indistinguishable at day 45 — the model does NOT support a relapse claim.
# =============================================================================
sc4_3d <- run_arm(ev_zuranolone(50, 3), win(nsp = ENROL))

# =============================================================================
#  SCENARIO 5 — THE BRIDGE TEST (the model's falsifiable structural claim)
#  Freeze receptor plasticity (KR = 0) and re-run the active arms.  The
#  ON-DRUG response should survive and the POST-WASHOUT benefit should not.
# =============================================================================
sc5_brx_norec <- run_arm(ev_brexanolone(top = 60),
                         c(win(care = c(ENROL, ENROL + 60),
                               protect = c(ENROL, ENROL + 60), nsp = ENROL),
                           list(KR = 0)))
sc5_zur_norec <- run_arm(ev_zuranolone(50, 14),
                         c(win(nsp = ENROL), list(KR = 0)))

# =============================================================================
#  SCENARIO 6 — SERTRALINE 50 -> 100 mg: the same endpoint, a different clock
#  SERT occupancy is immediate; the clinical effect is not, because it runs
#  through autoreceptor desensitisation and structural plasticity in series.
# =============================================================================
sc6 <- run_arm(ev_sertraline(), win(nsp = ENROL))

# =============================================================================
#  SCENARIO 7 — ESKETAMINE 0.25 mg/kg single dose (glutamatergic arm)
# =============================================================================
sc7 <- run_arm(ev_esketamine(), win(nsp = ENROL))

# =============================================================================
#  SCENARIO 8 — NON-PHARMACOLOGICAL: protected sleep block, and CBT/IPT
#  Protected sleep is not a control in this model — it acts on a real state.
# =============================================================================
sc8_sleep <- run_arm(params = win(protect = c(ENROL, TEND), nsp = ENROL))
sc8_cbt   <- run_arm(params = win(nsp = ENROL, cbt = ENROL))

# =============================================================================
#  SCENARIO 9 — ENROLMENT TIMING: why placebo response should shrink with
#  later enrolment (a testable prediction about trial design, not biology)
# =============================================================================
sc9 <- bind_rows(lapply(c(7, 14, 21, 42, 90, 180), function(day) {
  en <- TDEL + day * 24
  run_arm(params = win(care = c(en, en + 60), protect = c(en, en + 60),
                       nsp = en),
          end = en + 46 * 24) %>% mutate(enrol_day = day)
}))

# =============================================================================
#  SCENARIO 10 — TWO REJECTED REPAIRS (kept in the file because a failed
#  hypothesis that was actually run is worth more than one that was argued)
#  KR_BOOST : let PAM exposure accelerate delta-pool plasticity
#  W_PAM_BDNF : give the PAM a plasticity arm through BDNF
#  Both make the fit worse, because brexanolone IS allopregnanolone: any
#  parameter that rewards the drug also rewards pregnancy and lowers the
#  baseline the patient starts from.
# =============================================================================
sc10 <- bind_rows(
  run_arm(ev_brexanolone(60), c(win(care = c(ENROL, ENROL + 60),
                                    protect = c(ENROL, ENROL + 60),
                                    nsp = ENROL), list(KR_BOOST = 1.0))) %>%
    mutate(repair = "KR_BOOST = 1.0"),
  run_arm(ev_brexanolone(60), c(win(care = c(ENROL, ENROL + 60),
                                    protect = c(ENROL, ENROL + 60),
                                    nsp = ENROL), list(W_PAM_BDNF = 0.15))) %>%
    mutate(repair = "W_PAM_BDNF = 0.15"),
  sc2_60 %>% mutate(repair = "reference (both off)")
)

# =============================================================================
#  ENDPOINT TABLE — change from enrolment, all arms
# =============================================================================
delta_hamd <- function(d, at_days, enrol = ENROL) {
  base <- d$HAMD[which.min(abs(d$time - enrol))]
  sapply(at_days, function(dd)
    d$HAMD[which.min(abs(d$time - (enrol + dd * 24)))] - base)
}
days <- c(2.5, 3, 15, 28, 45)
endpoints <- tibble(
  arm = c("placebo (inpatient)", "brexanolone 60", "brexanolone 90",
          "zuranolone 50 x 14 d", "zuranolone 30 x 14 d",
          "zuranolone 50 x 3 d", "sertraline", "esketamine",
          "protected sleep", "CBT/IPT",
          "brexanolone 60 (KR = 0)", "zuranolone 50 (KR = 0)")
) %>% bind_cols(
  as_tibble(t(sapply(list(sc1, sc2_60, sc2_90, sc3_50, sc3_30, sc4_3d, sc6,
                          sc7, sc8_sleep, sc8_cbt, sc5_brx_norec,
                          sc5_zur_norec),
                     delta_hamd, at_days = days)),
            .name_repair = ~ paste0("d", days))
)
print(endpoints, n = Inf)

# =============================================================================
#  SENSITIVITY — one parameter at a time on the day-15 zuranolone endpoint
#  Expected ordering: KR and THR0 dominate, ZUR_EQ does not, because
#  potentiation is saturated at therapeutic exposure.
# =============================================================================
sens_keys <- c("KR", "EC50_PAM", "EMAX_PAM", "ZUR_EQ", "THR0", "KON", "KOFF",
               "KFAST", "W_SD", "A_SYMP", "KDEC_SD", "K_NSP", "V_KR")
sens <- bind_rows(lapply(sens_keys, function(k) {
  base <- as.numeric(param(mod)[[k]])
  vals <- sapply(c(0.7, 1.3), function(mult) {
    pl <- win(nsp = ENROL); pl[[k]] <- base * mult
    delta_hamd(run_arm(ev_zuranolone(50, 14), pl, end = ENROL + 16 * 24), 15)
  })
  tibble(parameter = k, lo = vals[1], hi = vals[2],
         span = abs(vals[2] - vals[1]))
})) %>% arrange(desc(span))
print(sens, n = Inf)

# =============================================================================
#  FIGURES
# =============================================================================
p_arms <- bind_rows(
  sc1    %>% mutate(arm = "placebo (inpatient)"),
  sc2_60 %>% mutate(arm = "brexanolone 60"),
  sc3_50 %>% mutate(arm = "zuranolone 50 x 14 d"),
  sc4_3d %>% mutate(arm = "zuranolone 50 x 3 d"),
  sc6    %>% mutate(arm = "sertraline")
) %>%
  filter(time >= ENROL - 24, time <= TEND) %>%
  mutate(day_from_enrol = (time - ENROL) / 24) %>%
  ggplot(aes(day_from_enrol, HAMD, colour = arm)) +
  geom_hline(yintercept = 7, linetype = 3) +
  geom_line(linewidth = 0.8) +
  labs(x = "days from enrolment", y = "HAM-D17",
       title = "Treatment arms; dotted line = remission threshold",
       subtitle = "V = 1.70, enrolment day 21 postpartum") +
  theme_bw()

p_bridge <- bind_rows(
  sc2_60        %>% mutate(arm = "brexanolone, KR normal"),
  sc5_brx_norec %>% mutate(arm = "brexanolone, KR = 0"),
  sc3_50        %>% mutate(arm = "zuranolone, KR normal"),
  sc5_zur_norec %>% mutate(arm = "zuranolone, KR = 0")
) %>%
  filter(time >= ENROL - 24) %>%
  mutate(day_from_enrol = (time - ENROL) / 24) %>%
  ggplot(aes(day_from_enrol, HAMD, colour = arm)) +
  geom_line(linewidth = 0.8) +
  labs(x = "days from enrolment", y = "HAM-D17",
       title = "The bridge test: durability belongs to the receptor",
       subtitle = "freezing delta-pool plasticity spares the on-drug response") +
  theme_bw()

# print(p_trough); print(p_arms); print(p_bridge)

# =============================================================================
#  NOTE ON RUNNING THIS FILE
#  -------------------------------------------------------------------------
#  Requires R with mrgsolve, dplyr, ggplot2, tidyr.  The build environment for
#  this repository has no R runtime, so the numbers quoted in the header and in
#  README.md were generated by ppd_reference_check.py — an independent
#  re-implementation of these identical equations in pure-stdlib Python with a
#  fixed-step RK4 integrator.  Discrepancies between this file and that one are
#  bugs in this file; ppd_reference_output.txt is the executed artefact.
# =============================================================================
