## =============================================================================
##  ros_mrgsolve_model.R
##  QSP model of ROSACEA
##  (erythematotelangiectatic · papulopustular · phymatous · ocular)
##
##  44 ODEs · KLK5/LL-37 innate amplifier · Demodex follicular ecology ·
##  TRPV1-CGRP neurovascular tone · VEGF angiogenesis · Th17/neutrophil
##  infiltrate · lymphatic stasis and phyma fibrosis · ocular surface ·
##  topical + systemic drug PK/PD · device modalities · 18 scenarios.
##
## -----------------------------------------------------------------------------
##  THE MODELLING IDEA
## -----------------------------------------------------------------------------
##  Rosacea is built here as ONE upstream amplifier feeding FOUR effector
##  states that differ only in TIME CONSTANT and in whether they remember:
##
##      STATE 1  TONE    (tau ~ 1 h,      no memory)  -> flushing
##      STATE 2  VDEN    (tau ~ 3 months, memory)     -> persistent erythema
##      STATE 3  PAP     (tau ~ 3 weeks,  memory)     -> papules/pustules
##      STATE 4  FIB+GLND(tau ~ years,    HYSTERESIS) -> phyma
##
##  There is no subtype switch anywhere in the code. Only four susceptibility
##  parameters ever move:
##
##      SPROT  stratum-corneum protease / KLK5 set-point   (1 = normal)
##      SNEUR  neurovascular (TRPV1-CGRP) gain             (1 = normal)
##      SMITE  follicular carrying capacity for Demodex    (1 = normal)
##      SFIBR  fibrogenic / glandular propensity           (1 = normal)
##
##  ETR is high SNEUR with low SMITE; PPR is high SMITE with high SPROT;
##  phyma is high SFIBR plus TIME. Everything else — including which drug
##  works on which endpoint — is an OUTPUT.
##
##  Consequences the model produces rather than assumes:
##
##    * brimonidine reads STATE 1 only. It gives the largest same-day drop in
##      CEA of any agent in the file and changes VDEN, PAP, FIB by nothing.
##      Two adaptation states sit under it — A2AR (receptor availability,
##      internalised by chronic occupancy) and VDILC (compensatory vasodilator
##      drive) — so after 8 weeks of daily use the TROUGH erythema, and the
##      erythema after withdrawal, OVERSHOOT baseline. Rebound is nowhere
##      coded as an adverse event; it is what those two ODEs do.
##
##    * doxycycline is anti-protease, not anti-parasitic. Across a TEN-FOLD dose
##      range (20 -> 200 mg) the lesion count falls -46% -> -79% while Demodex
##      density moves only -16% -> -23% -- and that residual movement is
##      INDIRECT (less KLK5 -> better barrier -> smaller follicular carrying
##      capacity), not acaricidal. At 40 mg the antibacterial term (DOXBIC50 =
##      2 mg/L) is engaged only ~23%, so the lesion response is carried by the
##      MMP-9/KLK5/IL-1beta arm. NOTE this is also the model prediction most
##      exposed to refutation: it still gives 100 mg a further lesion benefit
##      over 40 mg (-70% vs -57%), whereas the clinical claim for the
##      sub-antimicrobial dose is near-equivalence. If that equivalence is the
##      truth, IMAXK/IMAXM/IMAXI are too high, not the IC50s.
##
##    * ivermectin carries two actions on two clocks — a fast GluCl-mediated
##      mite kill and a slower LL-37/TLR2 damping. Because DEMO has a small
##      re-immigration term (IMMIG), clearing the reservoir buys TIME TO
##      RELAPSE rather than speed of response: the ivermectin-vs-metronidazole
##      separation in the model is larger after withdrawal than during
##      treatment, which is the shape of the ATTRACT extension data.
##
##    * pulsed-dye laser is the only intervention that removes STATE 2, so it
##      is the only one that lowers CEA without touching lesion count; and FIB
##      has KFL ~ 0, so no drug in the file reverses phyma once built.
##
##    * flush FREQUENCY has a floor set by TRIGEF x TRPV. No licensed agent in
##      the model touches TRPV, so every pharmacological scenario leaves flush
##      frequency near that floor — only trigger avoidance moves it.
##
## -----------------------------------------------------------------------------
##  UNITS
## -----------------------------------------------------------------------------
##  time            days
##  DEMO            mites / cm2 (standardised skin surface biopsy)
##  PAP             inflammatory lesion count (papules + pustules, full face)
##  TONE            dimensionless 0-1 (0 = maximally constricted)
##  BARR            dimensionless 0-1 (barrier integrity)
##  FIB, GLND, MGDX dimensionless 0-1 (fraction of maximal remodelling)
##  all mediators   dimensionless, 1.0 = healthy reference level
##  drug amounts    mg (systemic) or "applications" (topical, 1 = one dose)
##  CDOX, CISO      mg/L
##
## -----------------------------------------------------------------------------
##  CALIBRATION TARGETS (approximate literature anchors; see ros_references.md)
## -----------------------------------------------------------------------------
##  Demodex density (mites/cm2, SSSB)      healthy ~0.7  ·  PPR ~10-15
##    Forton 1993 Br J Dermatol; Casas 2012 Exp Dermatol
##  Cathelicidin LL-37 / KLK5 in lesional skin   ~10x and ~2-3x normal
##    Yamasaki 2007 Nat Med
##  Ivermectin 1% od, week 12 IGA success (0/1)  ~38-40% vs vehicle ~12-19%
##    inflammatory lesion reduction ~76% vs ~50%   Stein 2014 JAAD (2 trials)
##  Ivermectin 1% vs metronidazole 0.75%, week 16 lesion reduction ~83% vs ~74%
##    ATTRACT, Taieb 2015 Br J Dermatol
##  Median time to relapse after withdrawal   ~115 d (IVM) vs ~85 d (MTZ)
##    Taieb 2016 JEADV (ATTRACT extension)
##  Azelaic acid 15% bid, week 12 lesion reduction ~55-60% vs ~40% vehicle
##    Thiboutot 2003 JAAD
##  Doxycycline 40 mg MR, week 16 lesion reduction ~45-60% vs ~20-48% placebo
##    Del Rosso 2007 JAAD (2 trials); Cmax ~0.6 mg/L, sub-MIC
##  Brimonidine 0.33% gel: onset ~30 min, peak effect h 3-6, duration ~12 h;
##    composite (CEA+PSA >= 2 grades) at h 3 ~ 20-30% vs ~ 8-10% vehicle
##    Fowler 2012/2013 Br J Dermatol; worsening erythema ~9% over 1 year
##    Moore 2014 J Drugs Dermatol
##  Oxymetazoline 1% cream: composite at h 3-6 ~12-15% vs ~6% vehicle
##    REVEAL, Baumann 2018 JAAD
##  PDL 595 nm: ~30-50% erythema/telangiectasia clearance per session
##    Alam 2003; Tan 2017 systematic reviews
##  Ocular involvement in 50-75% of rosacea, poorly correlated with skin score
##    Vieira 2012 Clin Dermatol
##
##  These are ANCHORS, not a fit. The model is semi-quantitative: it is meant
##  to reproduce the ORDERING and the TIME COURSE of these observations, and
##  to be re-fitted before any quantitative claim is made.
##
## -----------------------------------------------------------------------------
##  WHAT THE MODEL ACTUALLY PRODUCES (verified with mrgsolve 1.5.2, 16 weeks,
##  PPR-moderate unless stated; reproduce with the calls at the end of the file)
## -----------------------------------------------------------------------------
##  phenotype (10-y burn-in)   DEMO   CEA   ILC   PHYGR  FLFREQ  OSDI
##    healthy                   0.8   0.16   0.0   0.00    0.24    6.1
##    ETR-moderate              1.9   2.92   5.6   0.17    4.64   25.3
##    PPR-moderate             11.5   2.63  14.1   0.47    3.64   39.1
##    phyma-prone male         19.4   2.87  17.9   1.66    4.05   41.3
##    ocular-dominant          10.1   2.32  12.5   0.20    3.03   48.5
##
##  16-week lesion-count change        model     published anchor
##    ivermectin 1% od                 -84%      -76 to -83%
##    metronidazole 0.75% bid          -66%      -74%
##    azelaic acid 15% bid             -61%      -55 to -60%
##    doxycycline 40 mg MR             -57%      -46 to -61%
##    minocycline 1.5% foam            -53%      -50 to -60%
##    isotretinoin 20 mg/day           -71%      -70 to -90%
##
##  emergent, not coded:
##    * laser CEA -18% with lesion count UNCHANGED; ivermectin lesion count
##      -83% with CEA only -9%  (the endpoint dissociation)
##    * brimonidine day-1 CEA 2.92 -> 2.18, and at DESENS = 2 the week-8
##      TROUGH (3.13) sits ABOVE the untreated baseline (2.92) with a
##      post-withdrawal peak of 3.25
##    * time to lose half the treatment gain: 59 d (ivermectin) vs 18 d
##      (metronidazole) -- the reservoir, not the potency, sets the clock
##    * flush frequency 4.64/day is unmoved by brimonidine, ivermectin,
##      doxycycline and laser alike; only trigger avoidance (1.75) and the
##      investigational TRPV1 antagonist (1.30) lower it
##    * phyma grade rises monotonically 1.66 -> 2.39 over 10 years; starting
##      isotretinoin at year 2 bends the slope (2.39 -> 1.41) but only
##      debulking steps it down
##
## -----------------------------------------------------------------------------
##  IMPORTANT — the model is initialised at the HEALTHY steady state. Setting
##  SPROT/SNEUR/SMITE/SFIBR away from 1 starts a transient that takes months
##  (STATE 2) to years (STATE 4) to settle. ALWAYS burn in with ros_steady()
##  or ros_init_at() before reading a phenotype or starting a drug.
## =============================================================================

library(mrgsolve)
library(dplyr)

ros_code <- '

$PROB
# Rosacea QSP model (ros)
- 44 ODEs. One innate amplifier (KLK5 -> LL-37) feeding four effector states
  with four different time constants. Subtypes are outputs of four
  susceptibility parameters, not switches.

$SET end = 112, delta = 1, maxsteps = 200000, rtol = 1e-6, atol = 1e-8

$GLOBAL
// ---- helpers -------------------------------------------------------------
//  The healthy steady state sits EXACTLY on the kink of every "excess over
//  reference" term (MMP9 - 1, ROS - 1, ...), so a hard max(0,x) puts a
//  derivative discontinuity right at the operating point and lsoda chatters
//  its step size down to machine epsilon. POS() is therefore a smooth
//  rectifier (softplus-like, agrees with max(0,x) to ~1e-4) and the 0-1
//  clamps used inside $ODE are smooth minima. The $TABLE clamps can stay
//  hard: they are not integrated.
#define SMOOTH_EPS 1e-6
#define POS(x)  (0.5 * ((x) + sqrt((x) * (x) + SMOOTH_EPS)))
#define SMIN(a,b) (0.5 * ((a) + (b) - sqrt(((a) - (b)) * ((a) - (b)) + SMOOTH_EPS)))
#define SCLAMP01(x) SMIN(POS(x), 1.0)
#define CAP1(x) ((x) > 1.0 ? 1.0 : ((x) < 0.0 ? 0.0 : (x)))

// assigned in $MAIN, used in $ODE and $TABLE
double VDOXL;      // doxycycline volume of distribution (L)
double VISOL;      // isotretinoin apparent volume (L)
double TONEREF;    // healthy tone reference
double VDREF;      // healthy vessel-density reference

$PARAM @annotated
// ---------------- susceptibility (the ONLY phenotype knobs) ---------------
SPROT   :  1.0  : stratum-corneum protease / KLK5 set-point multiplier ()
SNEUR   :  1.0  : neurovascular TRPV1-CGRP gain multiplier ()
SMITE   :  1.0  : follicular carrying capacity for Demodex multiplier ()
SFIBR   :  1.0  : fibrogenic and glandular propensity multiplier ()
ANDROG  :  1.0  : androgen / DHT drive on sebaceous gland ()
DESENS  :  1.0  : individual alpha-2A desensitisation rate multiplier ()

// ---------------- exposome -----------------------------------------------
//  Every exposome variable is read as an EXCESS over its healthy reference
//  (UV0, TRIG0, STRESS0). That is what makes the healthy state an exact fixed
//  point of the whole 44-ODE system instead of an approximate one: with all
//  susceptibilities at 1 and all exposures at reference, every mediator drive
//  evaluates to 1.0 and every derivative to 0.
TRIGB   :  0.30 : baseline trigger load (0-1 exposome score)
UVLOAD  :  0.20 : chronic UV exposure (0-1)
AVOID   :  0.0  : trigger-avoidance fraction achieved (0-1)
STRESSL :  0.20 : psychological stress load (0-1)
SKINCARE:  0.0  : gentle cleanser + ceramide moisturiser adherence (0-1)
PENB    :  0.50 : barrier-loss amplification of trigger penetration ()
UV0     :  0.20 : reference UV load (healthy skin) ()
TRIG0   :  0.40 : reference effective trigger load (healthy skin) ()
STRESS0 :  0.20 : reference stress load (healthy skin) ()

// ---------------- barrier -------------------------------------------------
BMAX    :  0.95 : maximal barrier integrity ()
KBARR   :  0.15 : barrier turnover rate (1/day)
BSENSM  :  0.12 : MMP-driven barrier loss sensitivity ()
BSENSK  :  0.10 : KLK-driven barrier loss sensitivity ()
ESKIN   :  0.06 : skincare effect on barrier target ()

// ---------------- sebaceous habitat --------------------------------------
SEB0    :  1.0  : reference sebum output ()
KSEB    :  0.04 : sebaceous gland adaptation rate (1/day)
ASEBAND :  0.60 : androgen effect on sebum ()

// ---------------- Demodex ecology ----------------------------------------
DCAP0   :  0.80 : healthy Demodex carrying capacity (mites/cm2)
KDG     :  0.10 : Demodex intrinsic growth rate (1/day)
KDCL    :  0.005: immune clearance of Demodex (1/day)
IMMIG   :  0.004: re-immigration / untreated-reservoir influx (mites/cm2/day)
ACAPPLUG:  1.20 : follicular plugging expansion of carrying capacity ()
KBOLIN  :  1.0  : Bacillus oleronius antigen production per mite (1/day)
KBOLOUT :  1.0  : antigen clearance (1/day)
ABURST  :  6.0  : dead-mite antigen burst amplification ()
BOLN    :  0.80 : healthy antigen reference level ()

// ---------------- innate amplifier ---------------------------------------
//  Every mediator target has the SATURATING form  1 + Emax * D/(D + K),
//  where D is a sum of reference-corrected drives. Two reasons:
//   (i) a chain of linear sums multiplies gains and the cascade runs away —
//       the first version of this model diverged at ~day 45;
//  (ii) real cytokine responses saturate, so the bound is also the physiology.
//  The loop gain around KLK5 (barrier + TLR2 + MMP-9 + ROS + mast cell) is
//  deliberately kept near 0.5, giving a closed-loop amplification of ~2x:
//  strong enough to be an amplifier, small enough to have a fixed point.
KKLK    :  0.50 : KLK5 activity turnover (1/day)
AKMAX   :  2.00 : maximal feedback contribution to KLK5 activity ()
KAKQ    :  3.00 : half-effect of the KLK5 feedback drive sum ()
AKPROT  :  0.90 : barrier-loss drive on KLK5 ()
AKROS   :  0.30 : ROS drive on KLK5 ()
AKMMP   :  0.35 : MMP-9 drive on KLK5 (feedback loop 2) ()
AKTLR   :  0.45 : TLR2 drive on KLK5 (feedback loop 1) ()
AKMC    :  0.20 : mast-cell tryptase drive on KLK5 ()
KLL     :  0.60 : LL-37 turnover (1/day)
LLMAX   : 12.0  : maximal LL-37 induction over healthy ()
KLLQ    :  6.00 : half-effect of KLK5 excess on LL-37 ()
AVITD   :  0.70 : UV/vitamin-D induction of cathelicidin precursor ()
KTL     :  0.40 : TLR2 expression turnover (1/day)
TLMAX   :  3.00 : maximal TLR2 induction ()
KTLQ    :  4.00 : half-effect of the TLR2 drive sum ()
ATLL    :  0.50 : LL-37 induction of TLR2 (loop 1) ()
ATLB    :  0.60 : Bacillus antigen induction of TLR2 ()
KROS    :  1.00 : ROS turnover (1/day)
RSMAX   :  3.00 : maximal ROS induction ()
KRSQ    :  3.00 : half-effect of the ROS drive sum ()
AROSUV  :  1.50 : UV drive on ROS ()
AROSN   :  0.50 : neutrophil/NET drive on ROS ()
AROSM   :  0.40 : dead-mite drive on ROS ()
KI1     :  0.60 : IL-1beta turnover (1/day)
I1MAX   :  6.00 : maximal IL-1beta induction ()
KI1Q    :  8.00 : half-effect of the IL-1beta drive sum ()
AI1T    :  0.45 : TLR2 x LL-37 (NLRP3) drive on IL-1beta ()
AI1R    :  0.30 : ROS drive on IL-1beta ()
KMM     :  0.50 : MMP-9 turnover (1/day)
MMMAX   :  5.00 : maximal MMP-9 induction ()
KMMQ    :  5.00 : half-effect of the MMP-9 drive sum ()
AMI1    :  0.45 : IL-1beta drive on MMP-9 ()
AMIL17  :  0.35 : IL-17 drive on MMP-9 ()
AMROS   :  0.30 : ROS drive on MMP-9 ()

// ---------------- cellular infiltrate ------------------------------------
KMC     :  0.08 : mast-cell density turnover (1/day)
MCMAX   :  3.00 : maximal mast-cell expansion ()
KMCQ    :  4.00 : half-effect of the mast-cell drive sum ()
AMCLL   :  0.70 : LL-37 / MRGPRX2 drive on mast cells ()
AMCSP   :  1.50 : stress / substance-P drive on mast cells ()
AMCTRIG :  1.00 : trigger drive on mast cells ()
KNE     :  0.30 : neutrophil turnover (1/day)
NEMAX   :  5.00 : maximal neutrophil influx ()
KNEQ    :  5.00 : half-effect of the neutrophil drive sum ()
ANI1    :  0.60 : IL-1beta / CXCL8 drive on neutrophils ()
ANIL17  :  0.55 : IL-17 drive on neutrophils (loop 3) ()
ANBOL   :  0.12 : Bacillus antigen drive on neutrophils ()
KTH     :  0.06 : Th17 expansion rate (1/day)
THMAX   :  3.00 : maximal Th17 expansion ()
KTHQ    :  4.00 : half-effect of the Th17 drive sum ()
ATHI1   :  0.55 : IL-1beta/IL-6 drive on Th17 ()
ATHMC   :  0.35 : mast-cell drive on Th17 ()
KI17    :  0.25 : IL-17A turnover (1/day)

// ---------------- angiogenesis / vessel structure ------------------------
KVE     :  0.30 : VEGF-A turnover (1/day)
VEMAX   :  3.00 : maximal VEGF-A induction ()
KVEQ    :  3.00 : half-effect of the VEGF drive sum ()
AVIL17  :  0.55 : IL-17 drive on VEGF ()
AVMC    :  0.50 : mast-cell drive on VEGF ()
AVROS   :  0.25 : ROS/HIF drive on VEGF ()
KVG     :  0.0035: vessel growth rate constant (1/day)
HVE     :  1.40 : Hill power of VEGF on vessel growth ()
KVCG    :  0.0022: neurotrophic (CGRP/NGF) vessel growth (1/day)
KVL     :  0.005: spontaneous vessel regression (1/day)
VDMAX   :  3.00 : maximal vessel density (x healthy)
KLAS    :  0.80 : laser ablation efficiency per unit fluence signal (1/day)
KLASEL  :  1.00 : laser signal elimination (1/day)

// ---------------- neurovascular ------------------------------------------
KTR     :  0.05 : TRPV1 expression turnover (1/day)
ATRSEN  :  0.90 : flush-memory sensitisation of TRPV1 (loop 4) ()
KSEN    : 60.0  : flush-memory half-effect (episode-days)
FLTHR   :  2.00 : flush rate above which neuroplastic sensitisation starts (/day)
ATRLL   :  0.40 : LL-37 drive on TRPV1 (burning in PPR) ()
KTRLQ   :  3.00 : half-effect of LL-37 excess on TRPV1 ()
KFM     :  0.008: flush-memory decay (1/day)
FRQMX   :  6.00 : maximal flush frequency (episodes/day)
KFRQ    :  1.20 : trigger x TRPV half-effect for flush frequency ()
FLOFF   :  0.35 : trigger x TRPV product below which flushing is incidental ()
KCG     :  2.00 : CGRP turnover (1/day)
CGMAX   :  4.00 : maximal CGRP induction ()
KCGQ    :  2.50 : half-effect of the CGRP drive ()
ACGT    :  1.50 : trigger x TRPV1 drive on CGRP ()
KNO     :  4.00 : nitric-oxide turnover (1/day)
NOMAX   :  2.00 : maximal NO induction ()
KNOQ    :  1.50 : half-effect of the NO drive sum ()
ANOMC   :  0.35 : mast-cell/histamine drive on NO ()
ANOI    :  0.20 : IL-1beta/iNOS drive on NO ()
ANOT    :  1.20 : direct trigger drive on NO ()
KTONE   : 24.0  : vascular tone equilibration rate (1/day, tau ~1 h)
TONE0   :  0.30 : healthy resting tone (0-1)
ATCG    :  0.55 : CGRP effect on tone ()
ATNO    :  0.35 : NO effect on tone ()
ATVD    :  0.55 : compensatory vasodilator drive effect on tone ()

// ---------------- alpha-adrenergic pharmacodynamics + adaptation ---------
BRMEC50 :  0.55 : brimonidine effect-site EC50 (applications)
EMXA2   :  0.70 : maximal alpha-2A vasoconstriction ()
KA2DES  :  0.09 : alpha-2A desensitisation rate (1/day)
KA2RES  :  0.05 : alpha-2A resensitisation rate (1/day)
KVCUP   :  0.06 : build-up of compensatory vasodilator drive (1/day)
KVCDN   :  0.05 : decay of compensatory vasodilator drive (1/day)
OXYEC50 :  0.60 : oxymetazoline effect-site EC50 (applications)
EMXA1   :  0.45 : maximal alpha-1A vasoconstriction ()
CARV    :  0.0  : carvedilol flag (0/1)
EMXBB   :  0.20 : beta-blocker effect on tone ()
CLONF   :  0.0  : clonidine flag (0/1)
EMXCL   :  0.15 : clonidine effect on tone ()

// ---------------- lesions ------------------------------------------------
KPG     :  0.080: lesion formation rate constant (lesions/day)
HNE     :  1.00 : Hill power of neutrophils on lesions ()
HIL     :  1.00 : Hill power of IL-17 on lesions ()
APDEM   :  1.50 : perifollicular (Demodex) amplification of lesions ()
KPDEM   :  6.00 : Demodex half-effect for lesion amplification (mites/cm2)
KPL     :  0.06 : lesion resolution rate (1/day)

// ---------------- oedema / lymphatics / phyma ----------------------------
KOE     :  0.010: oedema accumulation rate (1/day)
AOEV    :  0.60 : VEGF permeability drive on oedema ()
AOEM    :  0.50 : mast-cell permeability drive on oedema ()
ASTAS   :  0.80 : lymphatic-stasis self-amplification ()
KSTAS   :  0.40 : stasis half-effect ()
KOEL    :  0.020: oedema resolution (1/day)
KFG     :  0.00008: fibrosis accumulation rate (1/day)
ATGMC   :  0.60 : mast-cell TGF-beta drive ()
ATGI17  :  0.35 : IL-17 TGF-beta drive ()
ATGOE   :  0.70 : lymphstasis TGF-beta drive ()
FIBMAX  :  1.00 : maximal fibrosis ()
TGTHR   :  1.10 : TGF-beta drive below which fibrosis does not accumulate ()
KFL     :  0.00005: spontaneous fibrosis regression (1/day, ~0 = hysteresis)
KGG     :  0.00005: sebaceous hyperplasia accumulation rate (1/day)
AGLE    :  0.30 : LL-37/EGFR drive on glandular hyperplasia ()
GLMAX   :  1.00 : maximal glandular hyperplasia ()
GLTHR   :  1.00 : glandular drive below which hyperplasia does not accumulate ()
KGL     :  0.0002: glandular regression (1/day)
KDBLK   :  0.60 : debulking efficiency per unit signal (1/day)
KDBEL   :  1.00 : debulking signal elimination (1/day)

// ---------------- ocular -------------------------------------------------
KMGG    :  0.010: meibomian dysfunction accumulation (1/day)
ADBREV  :  0.80 : D. brevis drive on meibomian dysfunction (on Demodex EXCESS) ()
KMGD    :  8.00 : Demodex-excess half-effect at the lid (mites/cm2)
AMGK    :  0.10 : KLK5 drive on meibomian dysfunction ()
KMGL    :  0.010: meibomian recovery (1/day)
LIDHYG  :  0.0  : lid hygiene adherence (0-1)
ELID    :  2.00 : lid-hygiene amplification of recovery ()
IPLMG   :  0.0  : IPL-to-lid-margin flag (0/1)
EIPLM   :  1.50 : lid IPL amplification of recovery ()
KOC     :  0.15 : ocular surface inflammation turnover (1/day)
AOSM    :  1.20 : meibomian dysfunction drive on ocular inflammation ()
AOCI17  :  0.12 : IL-17 drive on ocular inflammation ()
AOCM    :  0.10 : MMP-9 drive on ocular inflammation ()

// ---------------- drug PK ------------------------------------------------
KIVMA   :  0.80 : ivermectin cream release into skin (1/day)
KIVMF   :  0.25 : ivermectin follicular elimination (1/day)
IVMEC50 :  1.50 : ivermectin GluCl EC50 for mite kill (applications)
IVMKMX  :  0.16 : maximal ivermectin mite kill rate (1/day)
IVMAIC50:  2.20 : ivermectin anti-inflammatory IC50 (applications)
KMTZ    :  1.20 : metronidazole skin elimination (1/day)
MTZRIC50:  0.35 : metronidazole ROS-scavenging IC50 (applications)
MTZNIC50:  0.22 : metronidazole neutrophil-suppression IC50 (applications)
KAZA    :  1.50 : azelaic acid skin elimination (1/day)
AZAKIC50:  0.55 : azelaic acid KLK5-inhibition IC50 (applications)
AZARIC50:  1.10 : azelaic acid antioxidant IC50 (applications)
KBRA    :  6.00 : brimonidine release into dermis (1/day)
KBRE    :  2.80 : brimonidine effect-site elimination (1/day, tau ~8.5 h)
KOXA    :  4.00 : oxymetazoline release into dermis (1/day)
KOXE    :  2.00 : oxymetazoline effect-site elimination (1/day)
KMIN    :  1.00 : minocycline foam follicular elimination (1/day)
MINEC50 :  1.20 : minocycline mite-kill IC50 (applications)
MINK    :  0.10 : maximal minocycline mite kill (1/day)
MINNIC50:  1.00 : minocycline neutrophil IC50 (applications)
MINBIC50:  1.00 : minocycline antibacterial IC50 (applications)
KADOX   :  4.00 : doxycycline absorption (1/day)
FDOXB   :  0.90 : doxycycline bioavailability ()
VDOX    : 60.0  : doxycycline volume of distribution (L)
CLDOX   : 48.0  : doxycycline clearance (L/day)
DOXMIC50:  0.10 : doxycycline MMP-9 inhibition IC50 (mg/L)
IMAXM   :  0.80 : maximal MMP-9 inhibition by doxycycline ()
DOXKIC50:  0.12 : doxycycline KLK5 inhibition IC50 (mg/L)
IMAXK   :  0.60 : maximal KLK5 inhibition by doxycycline ()
DOXIIC50:  0.15 : doxycycline IL-1beta inhibition IC50 (mg/L)
IMAXI   :  0.55 : maximal IL-1beta inhibition by doxycycline ()
DOXOIC50:  0.30 : doxycycline ocular MMP-9 IC50 (mg/L)
DOXBIC50:  2.00 : doxycycline ANTIBACTERIAL IC50 (mg/L) -- sub-MIC at 40 mg
KAISO   :  2.00 : isotretinoin absorption (1/day)
VISO    : 40.0  : isotretinoin apparent volume (L)
CLISO   : 33.0  : isotretinoin apparent clearance (L/day)
ISOEC50 :  0.30 : isotretinoin sebosuppression EC50 (mg/L)
ISOEMAX :  0.85 : maximal sebosuppression ()
ISOAIC50:  0.30 : isotretinoin innate (TLR2/KLK5) damping IC50 (mg/L)
PERMD   :  0.0  : permethrin/crotamiton acaricide flag (0-1)
PERMK   :  0.25 : permethrin mite kill rate (1/day)
TRPANT  :  0.0  : investigational TRPV1 antagonist exposure ()
TRPIC50 :  1.00 : TRPV1 antagonist IC50 ()
SECU    :  0.0  : investigational IL-17 blockade exposure ()
SECIC50 :  1.00 : IL-17 blockade IC50 ()
HCQ     :  0.0  : hydroxychloroquine exposure ()
HCQIC50 :  1.00 : hydroxychloroquine mast-cell IC50 ()

// ---------------- endpoint mapping weights -------------------------------
WS1     :  0.45 : weight of reversible tone on CEA ()
WS2     :  0.45 : weight of vessel structure on CEA ()
WS3     :  0.10 : weight of inflammatory redness on CEA ()
ERYT1   :  0.70 : tone excursion giving a full CEA unit-scale ()
ERYT2   :  1.50 : vessel-density excursion giving full CEA unit-scale ()
CEAOFF  :  0.15 : healthy-skin CEA offset ()
OSDI0   :  6.00 : normal-population OSDI floor ()
KQOL    : 12.0  : half-effect of the composite burden on DLQI ()

$INIT @annotated
//  HEALTHY initial conditions live HERE (as $INIT values), NOT in $MAIN:
//  ICs assigned in $MAIN are re-imposed on every simulation and would silently
//  override init(), so every "treated" run would start from healthy skin
//  instead of from the burned-in chronic state. ros_init_at() depends on this.
// ---- topical / systemic drug PK -----------------------------------------
IVMSK : 0    : ivermectin cream skin depot (applications)
IVMFO : 0    : ivermectin follicular effect site (applications)
MTZSK : 0    : metronidazole skin depot (applications)
AZASK : 0    : azelaic acid skin depot (applications)
BRMSK : 0    : brimonidine skin depot (applications)
BRMEF : 0    : brimonidine dermal effect site (applications)
OXYSK : 0    : oxymetazoline skin depot (applications)
OXYEF : 0    : oxymetazoline dermal effect site (applications)
MINSK : 0    : minocycline foam follicular depot (applications)
DOXG  : 0    : doxycycline gut depot (mg)
DOXP  : 0    : doxycycline central amount (mg)
ISOG  : 0    : isotretinoin gut depot (mg)
ISOP  : 0    : isotretinoin central amount (mg)
// ---- device signals -----------------------------------------------------
LASX  : 0    : vascular laser/IPL ablation signal (units)
DBLK  : 0    : phyma debulking signal (units)
// ---- barrier and habitat -----------------------------------------------
BARR  : 0.95 : barrier integrity (0-1)
SEB   : 1    : sebum output (x healthy)
DEMO  : 0.8  : Demodex density (mites/cm2)
BOL   : 0.8  : Bacillus oleronius antigen load ()
// ---- innate amplifier ---------------------------------------------------
KLK   : 1    : active KLK5 activity (x healthy)
LL37  : 1    : LL-37 cathelicidin (x healthy)
TLR2  : 1    : TLR2 expression (x healthy)
IL1B  : 1    : IL-1beta (x healthy)
MMP9  : 1    : MMP-9 activity (x healthy)
ROS   : 1    : oxidative stress (x healthy)
// ---- infiltrate ---------------------------------------------------------
MC    : 1    : mast cells (x healthy)
NEU   : 1    : neutrophils (x healthy)
TH17  : 1    : Th17 cells (x healthy)
IL17  : 1    : IL-17A (x healthy)
// ---- vasculature -------------------------------------------------------
VEGF  : 1    : VEGF-A (x healthy)
VDEN  : 1    : vessel density and calibre STATE 2 (x healthy)
// ---- neurovascular -----------------------------------------------------
TRPV  : 1    : TRPV1 expression / sensitisation (x healthy)
CGRP  : 1    : CGRP neuropeptide tone (x healthy)
NOX   : 1    : nitric oxide tone (x healthy)
TONE  : 0.3  : vascular tone STATE 1 (0-1)
A2AR  : 1    : alpha-2A receptor available fraction (0-1)
VDILC : 0    : compensatory vasodilator drive (rebound state)
FLMEM : 0    : cumulative flush memory (episode-days)
// ---- lesions and remodelling -------------------------------------------
PAP   : 0    : inflammatory lesion count STATE 3 (lesions)
OEDE  : 0    : dermal oedema / lymphatic stasis (0-1)
FIB   : 0    : dermal fibrosis STATE 4 (0-1)
GLND  : 0    : sebaceous hyperplasia STATE 4 (0-1)
// ---- ocular ------------------------------------------------------------
MGDX  : 0    : meibomian gland dysfunction (0-1)
OCUL  : 0    : ocular surface inflammation (0-1)
$MAIN
VDOXL   = VDOX;
VISOL   = VISO;
TONEREF = TONE0;
VDREF   = 1.0;

$ODE
// ===========================================================================
//  0. drug exposure
// ===========================================================================
double CDOX = DOXP / VDOXL;                 // mg/L
double CISO = ISOP / VISOL;                 // mg/L

// inhibition / engagement fractions ----------------------------------------
double FAZAK = 1.0 / (1.0 + AZASK / AZAKIC50);        // azelaic -> KLK5
double FAZAR = 1.0 / (1.0 + AZASK / AZARIC50);        // azelaic -> ROS
double FDOXK = 1.0 - IMAXK * CDOX / (CDOX + DOXKIC50);// doxy   -> KLK5
double FDOXM = 1.0 - IMAXM * CDOX / (CDOX + DOXMIC50);// doxy   -> MMP-9
double FDOXI = 1.0 - IMAXI * CDOX / (CDOX + DOXIIC50);// doxy   -> IL-1beta
double FDOXO = 1.0 / (1.0 + CDOX  / DOXOIC50);        // doxy   -> ocular
double FABX  = 1.0 / (1.0 + CDOX / DOXBIC50
                          + MINSK / MINBIC50);        // antibacterial term
double FIVMA = 1.0 / (1.0 + IVMFO / IVMAIC50);        // ivermectin -> LL-37
double FIVMT = 1.0 / (1.0 + IVMFO / (2.0 * IVMAIC50));// ivermectin -> TLR2
double FISOA = 1.0 / (1.0 + CISO / ISOAIC50);         // retinoid -> innate arm
double FMTZR = 1.0 / (1.0 + MTZSK / MTZRIC50);        // metronidazole -> ROS
double FMTZN = 1.0 / (1.0 + MTZSK / MTZNIC50);        // metronidazole -> PMN
double FMINN = 1.0 / (1.0 + MINSK / MINNIC50);        // minocycline -> PMN
double FTRPB = 1.0 / (1.0 + TRPANT / TRPIC50);        // TRPV1 antagonist
double FSECU = 1.0 / (1.0 + SECU / SECIC50);          // IL-17 blockade
double FHCQ  = 1.0 / (1.0 + HCQ / HCQIC50);           // HCQ -> mast cell

// ===========================================================================
//  1. barrier and effective trigger load
// ===========================================================================
double BTGT = BMAX * (1.0 + ESKIN * SKINCARE)
              / (1.0 + BSENSM * POS(MMP9 - 1.0) + BSENSK * POS(KLK - 1.0));
dxdt_BARR = KBARR * (BTGT - BARR);

double BLOSS = POS(BMAX - BARR) / BMAX;               // 0 = intact
double TRIGEF = SCLAMP01((TRIGB + 0.5 * UVLOAD) * (1.0 - AVOID)
                     * (1.0 + PENB * BLOSS));

// ===========================================================================
//  2. sebaceous habitat and Demodex ecology
// ===========================================================================
double EISO = ISOEMAX * CISO / (ISOEC50 + CISO);
double SEBT = SEB0 * (1.0 + ASEBAND * (ANDROG - 1.0));
dxdt_SEB = KSEB * (SEBT * (1.0 - EISO) - SEB);

double CAPD = DCAP0 * SMITE * SEB * (1.0 + ACAPPLUG * BLOSS);
double KILLI = IVMKMX * IVMFO / (IVMEC50 + IVMFO);    // ivermectin GluCl
double KILLM = MINK   * MINSK / (MINEC50 + MINSK);    // minocycline (weak)
double KILLP = PERMK  * PERMD;                        // permethrin
double KILLT = KILLI + KILLM + KILLP;
dxdt_DEMO = KDG * DEMO * (1.0 - DEMO / CAPD)
            - KILLT * DEMO - KDCL * DEMO + IMMIG;

double KFLUX = KILLT * DEMO;                          // dead-mite antigen burst
dxdt_BOL = KBOLIN * (DEMO + ABURST * KFLUX) * FABX - KBOLOUT * BOL;
double BOLR = BOL / BOLN;

// ===========================================================================
//  3. innate amplifier: KLK5 -> LL-37 -> TLR2 (loop 1) and MMP-9 (loop 2)
// ===========================================================================
//  KLK5 set-point = susceptibility PLUS a saturating feedback contribution.
//  Note SPROT is ADDITIVE, not multiplicative, on the feedback: if it scaled
//  the loop as well, loop gain would grow with susceptibility and severe
//  patients would have no fixed point at all.
double SKLK = AKPROT * BLOSS
            + AKROS * POS(ROS  - 1.0)
            + AKMMP * POS(MMP9 - 1.0)
            + AKTLR * POS(TLR2 - 1.0)
            + AKMC  * POS(MC   - 1.0);
double KLKD = SPROT + AKMAX * SKLK / (SKLK + KAKQ);
dxdt_KLK = KKLK * (KLKD * FAZAK * FDOXK * FISOA - KLK);

double QLL = POS(KLK - 1.0) * (1.0 + AVITD * (UVLOAD - UV0)) * FIVMA;
dxdt_LL37 = KLL * (1.0 + LLMAX * QLL / (QLL + KLLQ) - LL37);

double QTL = (ATLL * POS(LL37 - 1.0) + ATLB * POS(BOLR - 1.0)) * FIVMT * FISOA;
dxdt_TLR2 = KTL * (1.0 + TLMAX * QTL / (QTL + KTLQ) - TLR2);

double QRS = (AROSUV * (UVLOAD - UV0) + AROSN * POS(NEU - 1.0)
              + AROSM * KFLUX) * FMTZR * FAZAR;
dxdt_ROS = KROS * (1.0 + RSMAX * POS(QRS) / (POS(QRS) + KRSQ) - ROS);

double QI1 = (AI1T * POS(TLR2 * LL37 - 1.0) + AI1R * POS(ROS - 1.0)) * FDOXI;
dxdt_IL1B = KI1 * (1.0 + I1MAX * QI1 / (QI1 + KI1Q) - IL1B);

double QMM = (AMI1 * POS(IL1B - 1.0) + AMIL17 * POS(IL17 - 1.0)
              + AMROS * POS(ROS - 1.0)) * FDOXM;
dxdt_MMP9 = KMM * (1.0 + MMMAX * QMM / (QMM + KMMQ) - MMP9);

// ===========================================================================
//  4. cellular infiltrate (STATE 3 upstream)
// ===========================================================================
double QMC = (AMCLL * POS(LL37 - 1.0) + AMCSP * (STRESSL - STRESS0)
              + AMCTRIG * (TRIGEF - TRIG0)) * FHCQ;
dxdt_MC = KMC * (1.0 + MCMAX * POS(QMC) / (POS(QMC) + KMCQ) - MC);

double QNE = (ANI1 * POS(IL1B - 1.0) + ANIL17 * POS(IL17 - 1.0)
              + ANBOL * POS(BOLR - 1.0)) * FMTZN * FMINN;
dxdt_NEU = KNE * (1.0 + NEMAX * QNE / (QNE + KNEQ) - NEU);

double QTH = ATHI1 * POS(IL1B - 1.0) + ATHMC * POS(MC - 1.0);
dxdt_TH17 = KTH * (1.0 + THMAX * QTH / (QTH + KTHQ) - TH17);
dxdt_IL17 = KI17 * (1.0 + (TH17 - 1.0) * FSECU - IL17);

// ===========================================================================
//  5. neurovascular axis — STATE 1 plus its sensitisation memory (loop 4)
// ===========================================================================
//  Sensitisation has a THRESHOLD: only a flush rate above FLTHR lays down
//  neuroplastic memory, so a healthy person with a hot climate flushes without
//  becoming a rosacea patient — and above the threshold, loop 4 closes.
double TRD = 1.0 + ATRSEN * (FLMEM / (FLMEM + KSEN))
             + ATRLL * POS(LL37 - 1.0) / (POS(LL37 - 1.0) + KTRLQ);
dxdt_TRPV = KTR * (SNEUR * TRD * FTRPB - TRPV);

double QCG = ACGT * POS(TRIGEF * TRPV - TRIG0);
dxdt_CGRP = KCG * (1.0 + CGMAX * QCG / (QCG + KCGQ) - CGRP);

double QNO = ANOMC * POS(MC - 1.0) + ANOI * POS(IL1B - 1.0)
             + ANOT * POS(TRIGEF - TRIG0);
dxdt_NOX = KNO * (1.0 + NOMAX * QNO / (QNO + KNOQ) - NOX);

// alpha-adrenergic engagement and its two adaptation states -----------------
double OCCA2 = BRMEF / (BRMEF + BRMEC50);
double OCCA1 = OXYEF / (OXYEF + OXYEC50);
double EA2 = EMXA2 * A2AR * OCCA2;
double EA1 = EMXA1 * OCCA1;
double EBB = EMXBB * CARV + EMXCL * CLONF;

dxdt_A2AR  = KA2RES * (1.0 - A2AR) - KA2DES * DESENS * A2AR * OCCA2;
dxdt_VDILC = KVCUP * DESENS * OCCA2 - KVCDN * VDILC;

double ADIL = 1.0 + ATCG * POS(CGRP - 1.0) + ATNO * POS(NOX - 1.0)
              + ATVD * VDILC;
double TTGT = SCLAMP01(TONE0 * ADIL * (1.0 - EA2) * (1.0 - EA1) * (1.0 - EBB));
dxdt_TONE = KTONE * (TTGT - TONE);

double QFL = POS(TRIGEF * TRPV - FLOFF);
double FLRATE = FRQMX * QFL / (KFRQ + QFL);
dxdt_FLMEM = POS(FLRATE - FLTHR) - KFM * FLMEM;

// ===========================================================================
//  6. angiogenesis — STATE 2 (memory) and its only eraser, the laser
// ===========================================================================
double QVE = AVIL17 * POS(IL17 - 1.0) + AVMC * POS(MC - 1.0)
             + AVROS * POS(ROS - 1.0);
dxdt_VEGF = KVE * (1.0 + VEMAX * QVE / (QVE + KVEQ) - VEGF);

double HEADROOM = POS(1.0 - (VDEN - 1.0) / (VDMAX - 1.0));
dxdt_VDEN = (KVG * POS(pow(VEGF, HVE) - 1.0) + KVCG * POS(CGRP - 1.0))
            * HEADROOM
            - KVL * POS(VDEN - 1.0)
            - KLAS * LASX * POS(VDEN - 1.0);

// ===========================================================================
//  7. lesions — STATE 3
// ===========================================================================
double LESD = POS(pow(NEU, HNE) * pow(IL17, HIL) - 1.0)
              * (1.0 + APDEM * DEMO / (DEMO + KPDEM));
dxdt_PAP = KPG * LESD - KPL * PAP;

// ===========================================================================
//  8. oedema, lymphatic stasis, phyma — STATE 4 (hysteresis)
// ===========================================================================
double PERMD2 = AOEV * POS(VEGF - 1.0) + AOEM * POS(MC - 1.0);
dxdt_OEDE = KOE * PERMD2 * (1.0 + ASTAS * OEDE / (OEDE + KSTAS))
            * POS(1.0 - OEDE) - KOEL * OEDE;

double TGFD = ATGMC * POS(MC - 1.0) + ATGI17 * POS(IL17 - 1.0)
              + ATGOE * OEDE;
dxdt_FIB = KFG * SFIBR * POS(TGFD - TGTHR) * POS(1.0 - FIB / FIBMAX)
           - KFL * FIB - KDBLK * DBLK * FIB;

double GLD = POS(SEB - 1.0) + AGLE * POS(LL37 - 1.0);
dxdt_GLND = KGG * SFIBR * POS(GLD - GLTHR) * POS(1.0 - GLND / GLMAX)
            - KGL * GLND - KDBLK * DBLK * GLND;

// ===========================================================================
//  9. ocular surface
// ===========================================================================
double DEXC = POS(DEMO - DCAP0);                      // Demodex EXCESS
double MGD_IN = ADBREV * DEXC / (DEXC + KMGD) + AMGK * POS(KLK - 1.0);
dxdt_MGDX = KMGG * MGD_IN * POS(1.0 - MGDX)
            - KMGL * MGDX * (1.0 + ELID * LIDHYG + EIPLM * IPLMG);

double OCD = AOSM * MGDX + AOCI17 * POS(IL17 - 1.0) + AOCM * POS(MMP9 - 1.0);
dxdt_OCUL = KOC * (OCD * FDOXO - OCUL);

// ===========================================================================
// 10. drug PK
// ===========================================================================
dxdt_IVMSK = -KIVMA * IVMSK;
dxdt_IVMFO =  KIVMA * IVMSK - KIVMF * IVMFO;
dxdt_MTZSK = -KMTZ * MTZSK;
dxdt_AZASK = -KAZA * AZASK;
dxdt_BRMSK = -KBRA * BRMSK;
dxdt_BRMEF =  KBRA * BRMSK - KBRE * BRMEF;
dxdt_OXYSK = -KOXA * OXYSK;
dxdt_OXYEF =  KOXA * OXYSK - KOXE * OXYEF;
dxdt_MINSK = -KMIN * MINSK;
dxdt_DOXG  = -KADOX * DOXG;
dxdt_DOXP  =  KADOX * FDOXB * DOXG - (CLDOX / VDOXL) * DOXP;
dxdt_ISOG  = -KAISO * ISOG;
dxdt_ISOP  =  KAISO * ISOG - (CLISO / VISOL) * ISOP;
dxdt_LASX  = -KLASEL * LASX;
dxdt_DBLK  = -KDBEL * DBLK;

$TABLE
// --------------------------------------------------------------------------
//  clinical endpoints. Recomputed here from state values so that they are
//  exact at the output time.
// --------------------------------------------------------------------------
double CDOXO  = DOXP / VDOX;
double CISOO  = ISOP / VISO;
double BLOSSO = (BMAX - BARR) > 0.0 ? (BMAX - BARR) / BMAX : 0.0;
double TRIGO  = CAP1((TRIGB + 0.5 * UVLOAD) * (1.0 - AVOID)
                     * (1.0 + PENB * BLOSSO));

// erythema decomposition ---------------------------------------------------
double ERYS1 = CAP1((TONE - TONE0) / ERYT1);        // reversible / blanchable
double ERYS2 = CAP1((VDEN - 1.0) / ERYT2);          // structural
double ERYS3 = CAP1(PAP / 40.0);                    // inflammatory redness
double CEA = CEAOFF + 4.0 * (WS1 * ERYS1 + WS2 * ERYS2 + WS3 * ERYS3);
if (CEA > 4.0) CEA = 4.0;

// flushing ----------------------------------------------------------------
double QFLT = POS(TRIGO * TRPV - FLOFF);
double FLFREQ = FRQMX * QFLT / (KFRQ + QFLT);

// sensory / stinging ------------------------------------------------------
//  saturating rather than hard-capped, so severity keeps ordering patients
//  instead of piling them all on the ceiling
double SDRV = 3.0 * POS(TRPV - 1.0) + 5.0 * BLOSSO;
double STING = 10.0 * SDRV / (SDRV + 9.0);

double PSA = 0.85 * CEA + 0.12 * STING;
if (PSA > 4.0) PSA = 4.0;

// lesions and IGA ---------------------------------------------------------
double ILC = PAP;
double IGA = 0.22 * pow(POS(PAP), 0.70) + 0.25 * CEA;
if (IGA > 4.0) IGA = 4.0;

// telangiectasia and phyma ------------------------------------------------
double TELSC = 3.0 * CAP1((VDEN - 1.0) / 1.60);
double PHYGR = 3.0 * CAP1((0.6 * FIB + 0.4 * GLND) / 0.60);

// ocular ------------------------------------------------------------------
//  OSDI0 is the normal-population floor (a healthy lid is not OSDI 0)
double OSDI = OSDI0 + (100.0 - OSDI0)
              * CAP1(0.5 * OCUL / (OCUL + 2.0) + 0.5 * MGDX);

// erythema index (reflectance-spectrometry surrogate, a.u.) ---------------
double ERYIDX = 12.0 + 18.0 * (0.45 * ERYS1 + 0.55 * ERYS2);

// quality of life ---------------------------------------------------------
double QOLD = 0.80 * CEA + 0.06 * PAP + 1.20 * PHYGR + 0.25 * STING
              + 0.50 * FLFREQ + 0.03 * OSDI;
double DLQI = 30.0 * QOLD / (QOLD + KQOL);

// target engagement readouts ----------------------------------------------
double OCC_A2 = BRMEF / (BRMEF + BRMEC50);
double OCC_A1 = OXYEF / (OXYEF + OXYEC50);
double TE_A2  = EMXA2 * A2AR * OCC_A2;              // realised constriction
double TE_MMP = IMAXM * CDOXO / (CDOXO + DOXMIC50);
double TE_ABX = 1.0 - 1.0 / (1.0 + CDOXO / DOXBIC50);
double TE_KILL = IVMKMX * IVMFO / (IVMEC50 + IVMFO);

$CAPTURE @annotated
CEA    : clinician erythema assessment (0-4)
PSA    : patient self-assessment of erythema (0-4)
IGA    : investigator global assessment (0-4, success = 0/1)
ILC    : inflammatory lesion count (papules + pustules)
TELSC  : telangiectasia grade (0-3)
PHYGR  : phyma grade (0-3)
FLFREQ : flush frequency (episodes/day)
STING  : stinging / burning VAS (0-10)
DLQI   : dermatology life quality index (0-30)
OSDI   : ocular surface disease index (0-100)
ERYIDX : erythema index (a.u.)
ERYS1  : reversible (blanchable) erythema fraction ()
ERYS2  : structural (vessel) erythema fraction ()
TRIGO  : effective trigger load ()
CDOXO  : doxycycline plasma concentration (mg/L)
CISOO  : isotretinoin plasma concentration (mg/L)
OCC_A2 : alpha-2A receptor occupancy ()
TE_A2  : realised alpha-2A vasoconstriction ()
TE_MMP : MMP-9 inhibition fraction ()
TE_ABX : antibacterial target engagement fraction ()
TE_KILL: ivermectin mite-kill rate (1/day)
'

ros <- mread_cache("ros", code = ros_code, soloc = tempdir())

## =============================================================================
##  PHENOTYPE LIBRARY
##  Four susceptibility parameters. No subtype switch — the clinical subtype is
##  whatever those four numbers plus TIME produce.
## =============================================================================
ros_phenotypes <- function() {
  data.frame(
    phenotype = c("healthy", "ETR-mild", "ETR-moderate", "ETR-severe",
                  "PPR-mild", "PPR-moderate", "PPR-severe",
                  "mixed ETR+PPR", "phyma-prone male", "ocular-dominant"),
    SPROT  = c(1.0, 1.6, 2.2, 2.8, 2.4, 3.0, 3.6, 3.0, 3.0, 2.4),
    SNEUR  = c(1.0, 2.0, 2.8, 3.6, 1.4, 1.8, 2.0, 3.0, 2.0, 1.6),
    SMITE  = c(1.0, 1.4, 1.9, 2.4, 6.0, 11.0, 16.0, 11.0, 13.0, 10.0),
    SFIBR  = c(1.0, 1.0, 1.2, 1.5, 1.2, 1.5, 2.0, 1.8, 4.0, 1.2),
    ANDROG = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.1, 1.7, 1.0),
    ## lid-margin susceptibility: the ocular arm has its own slow state, which
    ## is why ocular severity tracks the cutaneous score so poorly
    ADBREV = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 2.2),
    TRIGB  = c(0.30, 0.45, 0.55, 0.70, 0.35, 0.40, 0.45, 0.60, 0.40, 0.35),
    UVLOAD = c(0.20, 0.30, 0.35, 0.45, 0.25, 0.30, 0.35, 0.35, 0.45, 0.25),
    stringsAsFactors = FALSE
  )
}

ros_pheno <- function(name, ...) {
  p <- ros_phenotypes()
  i <- match(name, p$phenotype)
  if (is.na(i)) stop("unknown phenotype: ", name,
                     "\navailable: ", paste(p$phenotype, collapse = ", "))
  out <- as.list(p[i, setdiff(names(p), "phenotype")])
  utils::modifyList(out, list(...))
}

## =============================================================================
##  BURN-IN HELPERS
##  Susceptibility parameters describe a chronic state; the four effector states
##  need months (STATE 2) to years (STATE 4) to express it. Never read a
##  phenotype without burning in first.
## =============================================================================
## NOTE on the integrator: TONE has tau ~ 1 h while FIB has tau ~ years, so a
## multi-year burn-in is a genuinely stiff problem. Keep the output grid coarse
## but NOT a single interval (lsoda counts steps per output interval), and leave
## maxsteps at the elevated $SET value.
ros_steady <- function(mod = ros, pars = list(), days = 3650, grid = 30) {
  m <- mod
  if (length(pars)) m <- param(m, pars)
  out <- as.data.frame(mrgsim(m, end = days, delta = grid))
  out[nrow(out), ]
}

ros_init_at <- function(mod = ros, pars = list(), days = 3650, grid = 30) {
  m <- mod
  if (length(pars)) m <- param(m, pars)
  last <- as.data.frame(mrgsim(m, end = days, delta = grid))
  last <- last[nrow(last), ]
  cmts <- names(init(m))
  iv <- as.numeric(last[1, cmts])
  names(iv) <- cmts
  init(m, as.list(iv))
}

## =============================================================================
##  REGIMEN BUILDERS
##  Topical doses are in "applications": 1 unit = one full-face application.
## =============================================================================
dose_ivermectin <- function(start = 0, weeks = 16, per_day = 1) {
  ev(amt = 1, cmt = "IVMSK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_metronidazole <- function(start = 0, weeks = 16, per_day = 2) {
  ev(amt = 1, cmt = "MTZSK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_azelaic <- function(start = 0, weeks = 16, per_day = 2) {
  ev(amt = 1, cmt = "AZASK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_brimonidine <- function(start = 0, weeks = 8, per_day = 1) {
  ev(amt = 1, cmt = "BRMSK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_oxymetazoline <- function(start = 0, weeks = 8, per_day = 1) {
  ev(amt = 1, cmt = "OXYSK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_minocycline_foam <- function(start = 0, weeks = 12, per_day = 1) {
  ev(amt = 1, cmt = "MINSK", time = start, ii = 1 / per_day,
     addl = round(weeks * 7 * per_day) - 1)
}
dose_doxycycline <- function(mg = 40, start = 0, weeks = 16) {
  ev(amt = mg, cmt = "DOXG", time = start, ii = 1,
     addl = round(weeks * 7) - 1)
}
dose_isotretinoin <- function(mg = 20, start = 0, weeks = 20) {
  ev(amt = mg, cmt = "ISOG", time = start, ii = 1,
     addl = round(weeks * 7) - 1)
}
## Laser / IPL: each session deposits one unit of ablation signal, which the
## model turns into a step reduction of VDEN (STATE 2). Nothing else moves.
dose_laser <- function(sessions = 3, interval = 28, start = 0, fluence = 1) {
  ev(amt = fluence, cmt = "LASX", time = start, ii = interval,
     addl = sessions - 1)
}
## Debulking (CO2 laser / electrosurgery) for phyma: the only exit from STATE 4.
dose_debulk <- function(sessions = 1, interval = 180, start = 0, extent = 1) {
  ev(amt = extent, cmt = "DBLK", time = start, ii = interval,
     addl = sessions - 1)
}

combine_rx <- function(...) {
  rs <- Filter(Negate(is.null), list(...))
  if (!length(rs)) return(NULL)
  Reduce(function(a, b) c(a, b), rs)
}

## =============================================================================
##  CORE RUNNER
## =============================================================================
##  pars     = PHENOTYPE parameters. Applied during the burn-in, so they
##             describe who the patient is before treatment starts.
##  rx_pars  = INTERVENTION parameters (trigger avoidance, lid hygiene, an
##             investigational antagonist). Applied ONLY AFTER the burn-in.
##  Mixing the two up is easy and silently wrong: an intervention left in the
##  burn-in is already at steady state at t = 0 and shows no effect at all.
ros_run <- function(phenotype = "PPR-moderate", regimen = NULL,
                    days = 365, delta = 1, pars = list(), rx_pars = list(),
                    burnin = 3650, mod = ros) {
  ph <- if (is.null(phenotype)) list() else ros_pheno(phenotype)
  ph <- utils::modifyList(ph, pars)
  m  <- ros_init_at(mod, ph, days = burnin)          # chronic pre-treatment state
  m  <- param(m, utils::modifyList(ph, rx_pars))     # then switch treatment on
  if (is.null(regimen)) {
    out <- mrgsim(m, end = days, delta = delta)
  } else {
    out <- mrgsim(m, events = regimen, end = days, delta = delta)
  }
  as.data.frame(out)
}

## =============================================================================
##  SCENARIOS (18)
##  S1-S3   untreated phenotypes
##  S4-S12  monotherapies
##  S13-S15 combinations and devices
##  S16-S18 withdrawal / relapse, phyma trajectory, ocular arm
## =============================================================================
ros_scenarios <- function(weeks = 16) {
  W <- weeks
  list(
    S1  = list(label = "S1 healthy skin",
               ph = "healthy", rx = NULL, days = 7 * W),
    S2  = list(label = "S2 ETR moderate, untreated",
               ph = "ETR-moderate", rx = NULL, days = 7 * W),
    S3  = list(label = "S3 PPR moderate, untreated",
               ph = "PPR-moderate", rx = NULL, days = 7 * W),
    S4  = list(label = "S4 ivermectin 1% od",
               ph = "PPR-moderate", rx = dose_ivermectin(weeks = W),
               days = 7 * W),
    S5  = list(label = "S5 metronidazole 0.75% bid",
               ph = "PPR-moderate", rx = dose_metronidazole(weeks = W),
               days = 7 * W),
    S6  = list(label = "S6 azelaic acid 15% bid",
               ph = "PPR-moderate", rx = dose_azelaic(weeks = W),
               days = 7 * W),
    S7  = list(label = "S7 doxycycline 40 mg MR (sub-antimicrobial)",
               ph = "PPR-moderate", rx = dose_doxycycline(40, weeks = W),
               days = 7 * W),
    S8  = list(label = "S8 doxycycline 100 mg (antimicrobial)",
               ph = "PPR-moderate", rx = dose_doxycycline(100, weeks = W),
               days = 7 * W),
    S9  = list(label = "S9 minocycline 1.5% foam od",
               ph = "PPR-moderate", rx = dose_minocycline_foam(weeks = W),
               days = 7 * W),
    S10 = list(label = "S10 brimonidine 0.33% gel od (ETR)",
               ph = "ETR-moderate", rx = dose_brimonidine(weeks = W),
               days = 7 * W),
    S11 = list(label = "S11 oxymetazoline 1% cream od (ETR)",
               ph = "ETR-moderate", rx = dose_oxymetazoline(weeks = W),
               days = 7 * W),
    S12 = list(label = "S12 isotretinoin 20 mg/day",
               ph = "PPR-severe", rx = dose_isotretinoin(20, weeks = W),
               days = 7 * W),
    S13 = list(label = "S13 ivermectin + brimonidine (different states)",
               ph = "mixed ETR+PPR",
               rx = combine_rx(dose_ivermectin(weeks = W),
                               dose_brimonidine(weeks = W)),
               days = 7 * W),
    S14 = list(label = "S14 PDL x3 (STATE 2 deletion)",
               ph = "ETR-moderate", rx = dose_laser(3, 28), days = 7 * W),
    S15 = list(label = "S15 doxycycline 40 mg + PDL x3",
               ph = "mixed ETR+PPR",
               rx = combine_rx(dose_doxycycline(40, weeks = W),
                               dose_laser(3, 28)),
               days = 7 * W),
    S16 = list(label = "S16 trigger avoidance + skincare only",
               ph = "ETR-moderate", rx = NULL,
               rx_pars = list(AVOID = 0.6, SKINCARE = 1, UVLOAD = 0.05),
               days = 7 * W),
    S17 = list(label = "S17 ocular rosacea: doxy 40 mg + lid hygiene",
               ph = "ocular-dominant",
               rx = dose_doxycycline(40, weeks = W),
               rx_pars = list(LIDHYG = 1), days = 7 * W),
    S18 = list(label = "S18 investigational: TRPV1 antagonist",
               ph = "ETR-moderate", rx = NULL,
               rx_pars = list(TRPANT = 3), days = 7 * W)
  )
}

ros_run_scenario <- function(key, weeks = 16, delta = 1) {
  sc <- ros_scenarios(weeks)[[key]]
  if (is.null(sc)) stop("unknown scenario: ", key)
  df <- ros_run(sc$ph, sc$rx, days = sc$days, delta = delta,
                pars = if (is.null(sc$pars)) list() else sc$pars,
                rx_pars = if (is.null(sc$rx_pars)) list() else sc$rx_pars)
  df$scenario <- sc$label
  df$key <- key
  df
}

ros_run_all <- function(keys = names(ros_scenarios()), weeks = 16, delta = 1) {
  do.call(rbind, lapply(keys, ros_run_scenario, weeks = weeks, delta = delta))
}

## =============================================================================
##  SPECIAL-PURPOSE EXPERIMENTS
##  Each one exists to test a claim the map makes, not to illustrate a drug.
## =============================================================================

## E1 — brimonidine within-day profile plus the withdrawal rebound.
##      Look for: a large trough/rebound CEA above the untreated baseline at
##      day (weeks*7 + 1..21), largest when DESENS is high. Nothing in the
##      code says "rebound"; A2AR and VDILC produce it.
ros_rebound <- function(weeks = 8, follow = 28, desens = c(0.5, 1, 2),
                        phenotype = "ETR-moderate") {
  base <- ros_run(phenotype, NULL, days = weeks * 7 + follow, delta = 0.05)
  base$arm <- "untreated"
  arms <- lapply(desens, function(d) {
    df <- ros_run(phenotype, dose_brimonidine(weeks = weeks),
                  days = weeks * 7 + follow, delta = 0.05,
                  pars = list(DESENS = d))
    df$arm <- paste0("brimonidine, DESENS = ", d)
    df
  })
  do.call(rbind, c(list(base), arms))
}

## E2 — the ATTRACT question: ivermectin vs metronidazole DURING treatment and
##      AFTER withdrawal. The mite reservoir (IMMIG) is what separates the
##      relapse curves; on-treatment separation is much smaller.
ros_relapse <- function(weeks = 16, follow = 180,
                        phenotype = "PPR-moderate") {
  arms <- list(
    ivermectin    = dose_ivermectin(weeks = weeks),
    metronidazole = dose_metronidazole(weeks = weeks)
  )
  out <- lapply(names(arms), function(nm) {
    df <- ros_run(phenotype, arms[[nm]], days = weeks * 7 + follow, delta = 1)
    df$arm <- nm
    df
  })
  df <- do.call(rbind, out)
  df$phase <- ifelse(df$time <= weeks * 7, "treatment", "withdrawal")
  df
}

## Time to lose half of the gain achieved at the end of treatment. Reported
## this way rather than "half of baseline" because arms that achieved different
## depths of response are otherwise not comparable.
ros_relapse_times <- function(df, weeks = 16, var = "ILC") {
  tstop <- weeks * 7
  do.call(rbind, lapply(split(df, df$arm), function(s) {
    b <- s[[var]][1]
    e <- s[[var]][which.min(abs(s$time - tstop))]
    thr <- e + 0.5 * (b - e)
    post <- s[s$time > tstop, ]
    hit <- post$time[post[[var]] >= thr][1]
    data.frame(arm = s$arm[1], baseline = round(b, 1),
               end_of_treatment = round(e, 1),
               pct_reduction = round(100 * (e - b) / b, 1),
               days_to_lose_half_the_gain =
                 if (is.na(hit)) NA_real_ else hit - tstop)
  }))
}

## E3 — the falsifiable doxycycline prediction: lesion count falls, Demodex
##      density does not move at 40 mg; 100 mg engages the antibacterial term
##      and still adds little to the lesion endpoint.
ros_doxy_dissociation <- function(weeks = 16, phenotype = "PPR-moderate") {
  do.call(rbind, lapply(c(20, 40, 100, 200), function(mg) {
    df <- ros_run(phenotype, dose_doxycycline(mg, weeks = weeks),
                  days = weeks * 7, delta = 1)
    df$dose_mg <- mg
    df
  }))
}

## E4 — endpoint dissociation: laser moves CEA and not ILC; ivermectin the
##      reverse. Same patient, same duration.
ros_endpoint_cross <- function(weeks = 16, phenotype = "mixed ETR+PPR") {
  arms <- list(
    none        = NULL,
    ivermectin  = dose_ivermectin(weeks = weeks),
    laser       = dose_laser(3, 28),
    both        = combine_rx(dose_ivermectin(weeks = weeks), dose_laser(3, 28))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    df <- ros_run(phenotype, arms[[nm]], days = weeks * 7, delta = 1)
    df$arm <- nm
    df
  }))
}

## E5 — the phyma clock: 10 years of STATE 4 with and without early control,
##      then a debulking procedure. KFL is ~0, so the fibrotic state does not
##      run backwards on its own.
ros_phyma <- function(years = 10, phenotype = "phyma-prone male") {
  days <- round(years * 365)
  arms <- list(
    untreated = NULL,
    isotretinoin_from_year2 = dose_isotretinoin(20, start = 730,
                                                weeks = round((days - 730) / 7)),
    debulk_at_year8 = dose_debulk(1, start = round(8 * 365))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    df <- ros_run(phenotype, arms[[nm]], days = days, delta = 7)
    df$arm <- nm
    df
  }))
}

## E6 — flush frequency floor: no pharmacological arm in the model lowers
##      flush frequency much, because none of them touch TRPV. Only AVOID and
##      the investigational TRPV1 antagonist do.
ros_flush_floor <- function(weeks = 16, phenotype = "ETR-moderate") {
  arms <- list(
    untreated    = list(rx = NULL, rp = list()),
    brimonidine  = list(rx = dose_brimonidine(weeks = weeks), rp = list()),
    ivermectin   = list(rx = dose_ivermectin(weeks = weeks), rp = list()),
    doxycycline  = list(rx = dose_doxycycline(40, weeks = weeks), rp = list()),
    laser        = list(rx = dose_laser(3, 28), rp = list()),
    avoidance    = list(rx = NULL, rp = list(AVOID = 0.6, UVLOAD = 0.05)),
    TRPV1_antag  = list(rx = NULL, rp = list(TRPANT = 3))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    a <- arms[[nm]]
    df <- ros_run(phenotype, a$rx, days = weeks * 7, delta = 1, rx_pars = a$rp)
    df$arm <- nm
    df
  }))
}

## =============================================================================
##  VIRTUAL POPULATION — sweep the four susceptibility parameters and let the
##  clinical subtype fall out. Nothing in the model assigns a subtype label.
## =============================================================================
ros_vpop <- function(n = 200, seed = 20260728, weeks = 0) {
  set.seed(seed)
  d <- data.frame(
    SPROT  = exp(rnorm(n, log(2.2), 0.35)),
    SNEUR  = exp(rnorm(n, log(2.0), 0.40)),
    SMITE  = exp(rnorm(n, log(6.0), 0.90)),
    SFIBR  = exp(rnorm(n, log(1.5), 0.50)),
    TRIGB  = pmin(1, pmax(0.1, rnorm(n, 0.45, 0.15))),
    UVLOAD = pmin(1, pmax(0.0, rnorm(n, 0.30, 0.12)))
  )
  res <- do.call(rbind, lapply(seq_len(n), function(i) {
    ss <- ros_steady(ros, as.list(d[i, ]), days = 3650)
    data.frame(id = i, d[i, ],
               DEMO = round(ss$DEMO, 2), CEA = round(ss$CEA, 2),
               ILC = round(ss$ILC, 1), TELSC = round(ss$TELSC, 2),
               FLFREQ = round(ss$FLFREQ, 2), PHYGR = round(ss$PHYGR, 2),
               OSDI = round(ss$OSDI, 1), DLQI = round(ss$DLQI, 1),
               stringsAsFactors = FALSE)
  }))
  ## subtype is DERIVED, never assigned
  res$subtype <- with(res, ifelse(ILC >= 10 & CEA >= 2, "mixed",
                           ifelse(ILC >= 10, "PPR",
                             ifelse(CEA >= 2, "ETR", "subclinical"))))
  res
}

## =============================================================================
##  QUICK-LOOK HELPERS
## =============================================================================
ros_summary <- function(df) {
  last <- df[df$time == max(df$time), ]
  cols <- c("scenario", "CEA", "PSA", "IGA", "ILC", "TELSC", "PHYGR",
            "FLFREQ", "STING", "DLQI", "OSDI", "DEMO", "VDEN", "TONE",
            "KLK", "LL37", "IL17")
  cols <- intersect(cols, names(last))
  last[, cols]
}

## Percent change from the first output row, per arm — the way phase-3 papers
## report inflammatory lesion counts.
ros_pct_change <- function(df, var = "ILC", group = "arm") {
  g <- if (group %in% names(df)) df[[group]] else "all"
  do.call(rbind, lapply(split(df, g), function(s) {
    b <- s[[var]][which.min(s$time)]
    data.frame(group = s[[group]][1], time = s$time,
               value = s[[var]],
               pct = 100 * (s[[var]] - b) / ifelse(b == 0, NA, b))
  }))
}

## =============================================================================
##  WHAT TO LOOK FOR (none of it is coded as a rule)
## -----------------------------------------------------------------------------
##   ros_summary(ros_run_all(c("S3","S4","S5","S7")))
##     -> ivermectin beats metronidazole on ILC, and doxycycline 40 mg gets
##        most of the way there with DEMO untouched.
##
##   d <- ros_doxy_dissociation(); tapply(d$DEMO, d$dose_mg, function(x) tail(x,1))
##     -> Demodex density essentially flat across 20-200 mg, while
##        tapply(d$ILC, ...) falls: anti-protease, not anti-parasitic.
##
##   r <- ros_rebound(); subset(r, time > 56)   # after withdrawal
##     -> CEA in the brimonidine arms crosses ABOVE the untreated arm; the
##        crossing is larger and later for DESENS = 2.
##
##   x <- ros_endpoint_cross(); ros_summary(x)
##     -> laser: CEA down, ILC unchanged. ivermectin: ILC down, CEA barely.
##
##   p <- ros_phyma(); subset(p, arm == "untreated" & time %% 365 == 0)$PHYGR
##     -> monotone rise; the isotretinoin arm bends the slope but never
##        returns; only the debulking arm steps down.
##
##   f <- ros_flush_floor(); tapply(f$FLFREQ, f$arm, function(x) tail(x, 1))
##     -> every drug arm sits at the same floor; only avoidance and the
##        investigational TRPV1 antagonist move it.
##
##   v <- ros_vpop(); table(v$subtype)
##     -> ETR / PPR / mixed emerge from four continuous parameters.
## =============================================================================
