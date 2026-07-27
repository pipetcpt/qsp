## =====================================================================
##  Refractory / Unexplained Chronic Cough (RCC) — mrgsolve QSP model
##  Cough Hypersensitivity Syndrome
##  rcc_mrgsolve_model.R
## ---------------------------------------------------------------------
##  THE PREMISE
##
##  Chronic cough was managed for fifty years as a search for a cause.
##  That paradigm fails in 20-40% of specialist-clinic patients, and it
##  fails informatively: the "causes" get treated, the objective evidence
##  of them disappears, and the cough does not. The modern reframing is
##  that refractory chronic cough is a NEUROPATHIC disorder of the vagal
##  reflex arc itself — the airway analogue of neuropathic pain, with
##  hypertussia (hyperalgesia), allotussia (allodynia) and laryngeal
##  paraesthesia (spontaneous pain).
##
##  The model is built as THREE LAYERS IN SERIES, because every drug
##  class in chronic cough acts on exactly one of them and the layer
##  decides which endpoint the drug can possibly move:
##
##    LAYER 1  PERIPHERAL GAIN   epithelial damage -> pannexin-1 -> ATP
##             (the airway)      -> P2X3 on vagal C-fibres; NGF/TrkA
##                               up-regulating P2X3 itself over WEEKS;
##                               PGE2/BK/H+ sensitising TRPV1/TRPA1.
##                               >> P2X3 antagonists act here.
##
##    LAYER 2  CENTRAL GAIN      afferent barrage -> nTS NMDA wind-up,
##             (the brainstem)   substance P/NK-1, microglial BDNF,
##                               loss of GABA/glycinergic inhibition.
##                               This layer OUTLIVES ITS TRIGGER.
##                               >> gabapentin, morphine, amitriptyline
##                                  act here.
##
##    LAYER 3  CORTICAL CONTROL  drive becomes the conscious URGE TO
##             (the person)      COUGH before it becomes a motor cough,
##                               and prefrontal suppression can gate the
##                               motor output.
##                               >> speech therapy acts here — and so
##                                  does the 30% PLACEBO RESPONSE.
##
##      COUGH = f( AFFERENT DRIVE x CENTRAL GAIN  -  THRESHOLD )
##      with the threshold set by descending inhibition and by cortical
##      suppression, and the whole motor output GATED BY WAKEFULNESS.
##
## ---------------------------------------------------------------------
##  FOUR STRUCTURAL COMMITMENTS
##
##  1. THE DISEASE IS GENERATED, NOT ASSUMED.
##     Every state is written in the canonical relaxation form
##         dX/dt = KX * (XTGT(drivers) - X)
##     with XTGT equal to the healthy value when all drivers are at
##     their healthy values. A virtual subject with INSULT0 = 0 is
##     therefore weight-... cough-stable forever, and diagnostic 1
##     verifies ~0% drift over 365 days on every reported state.
##     Chronic cough EMERGES from an insult plus one susceptibility
##     parameter (LGAIN, the loop gain).
##
##  2. SUSCEPTIBILITY IS CONTINUOUS, NOT BISTABLE (a negative result).
##     The circle
##         afferent barrage -> nTS wind-up -> central gain -> cough
##         -> laryngeal mechanical trauma -> epithelial damage -> ATP
##         -> afferent barrage
##     closes on itself, and it was built expecting it to produce two
##     attractors with hysteresis. IT DOES NOT. Diagnostic 3 sweeps the
##     loop gain UP from health and DOWN from established disease and the
##     two branches coincide to within 0.3 coughs/h. The model instead
##     has a steep monotonic susceptibility relation. The negative result
##     is reported rather than removed, and the header claim it refuted
##     has been withdrawn. Post-viral resolution versus lifelong
##     refractory cough IS reproduced (post-viral settles at 1.8
##     coughs/h, RCC at 21), but by moving along that continuous curve —
##     loop gain plus oestrogen status plus age — not by crossing a fold.
##
##  3. EFFICACY AND DYSGEUSIA ARE ONE OCCUPANCY CURVE, SEPARATED BY ONE
##     NUMBER. Gefapixant, camlipixant, eliapixant and sivopixant share
##     a single efficacy mechanism (P2X3 homotrimer blockade on vagal
##     afferents) and a single adverse-effect mechanism (P2X2/3
##     heterotrimer blockade on type II/III taste cells, where ATP is
##     THE neurotransmitter for all five taste qualities). They differ
##     in the selectivity ratio SEL = IC50(P2X2/3) / IC50(P2X3):
##         gefapixant ~6x   camlipixant ~1500x
##         eliapixant ~70x  sivopixant ~40x
##     The taste system has NO RECEPTOR RESERVE (TDG50 = 0.45 with
##     Hill 2), so 6-fold selectivity is not remotely enough. Diagnostic
##     6 sweeps SEL over four decades and prints the therapeutic window.
##
##  4. EFFICACY vs EFFECTIVENESS — AND A SECOND NEGATIVE RESULT.
##     Dysgeusia drives discontinuation (the model loses ~11% of the arm
##     by week 12 at gefapixant 45 mg BID), discontinuation erodes the
##     population-average exposure, and the eroded exposure is what an
##     intention-to-treat estimate measures. The model reports both. It
##     was built expecting that gap to be a large part of why gefapixant
##     has such a small effect size. IT IS NOT: diagnostic 8 gives
##     ITT -18.6% against per-protocol -19.2%, a difference of 0.6
##     percentage points. Essentially all of the small effect size comes
##     from the P2X3-dependent share of tussive drive being small, and
##     essentially none of it from dropout. The expectation is left in
##     the file next to the number that refuted it.
##
## ---------------------------------------------------------------------
##  WHY THE PLACEBO ARM IS A MECHANISM AND NOT A NUMBER
##
##  Chronic cough is the therapeutic area with the largest placebo
##  response in OBJECTIVE physiology anywhere in respiratory medicine:
##  ~30-45% falls in machine-counted 24-h cough frequency. A model that
##  applied drug effects to an untreated natural history would overstate
##  every compound by a factor of two or more. Three mechanisms are
##  built in and all three operate in the placebo arm:
##
##    (a) ENROLMENT AT A PEAK. Patients present when their cough is at
##        its worst. FLARE is an above-set-point excursion that decays
##        with a ~40-day half-life whatever anyone does. This is
##        regression to the mean, expressed mechanistically.
##    (b) TRIAL PARTICIPATION TRAINS SUPPRESSION. Diary-keeping,
##        monitoring and clinician contact raise cortical suppression
##        capacity CORTS with a ~3-week time constant — the same
##        mechanism speech therapy uses deliberately, obtained
##        accidentally.
##    (c) NATURAL RESOLUTION of the sensitising insult.
##
##  Every reported drug effect below is placebo-ADJUSTED, i.e. computed
##  against a comparator that carries (a)+(b)+(c).
##
## ---------------------------------------------------------------------
##  THE DIAGNOSTIC DISSOCIATION THE MODEL PREDICTS RATHER THAN ASSUMES
##
##  Capsaicin gates TRPV1 directly and BYPASSES P2X3; inhaled ATP gates
##  P2X3. The model computes both challenge thresholds from the same
##  state vector, so it predicts — without being told — that P2X3
##  antagonists shift the ATP threshold by nearly a log and barely move
##  capsaicin C5. That is the observed signature of the class and it is
##  the reason capsaicin challenge is a poor pharmacodynamic biomarker
##  for these drugs.
##
## ---------------------------------------------------------------------
##  STRUCTURE — 68 ODE compartments
##
##   Circadian / gating (2)   CIRCS, CIRCC
##   Airway (7)               EPID, ATPX, MUC, INFL, EOSN, ROSX, PEPS
##   Reflux (2)               ACIDR, NACID
##   Mediators (4)            NGF, PGE2, BKN, INSU
##   Trial artefact (1)       FLARE
##   Receptors / periphery(6) P2X3E, TRPV, TRPA, NAVS, PSEN, ECTO
##   Afferent (1)             AFIR
##   Central (6)              SPC, WIND, MICG, GABI, DESC, HABT
##   Cortical (3)             HVIG, CORT, SLPD
##   Output states (5)        CFQ, URG, LCQ, VASC, CSDS
##   Accumulators (4)         CACC, CAWK, TAWK, COMP
##   Safety (3)               DYSG, ADHR, ALTX
##   Drug PK (24)             gefapixant 3, camlipixant 2, eliapixant 2,
##                            sivopixant 2, gabapentin 2, morphine 2,
##                            nalbuphine 2, amitriptyline 2, ICS 2,
##                            PPI 2, pregabalin 2, dextromethorphan 1
##
##   Time unit: DAYS. Cough frequencies are reported in coughs/HOUR.
##
## ---------------------------------------------------------------------
##  CALIBRATION ANCHORS (full citations in rcc_references.md)
##
##   Healthy       awake cough ~2/h, 24-h ~1.4/h, LCQ 21, capsaicin
##                 C5 ~ 60 uM (male) / ~ 25 uM (female).
##   RCC baseline  awake cough ~30/h, 24-h ~20/h (COUGH-1/2 geometric
##                 means 18-27/h), LCQ 10.5, cough VAS 57 mm, CSD 5.2,
##                 duration years, F:M ~ 2:1, peak age 55-65.
##   COUGH-1 12wk  gefapixant 45 mg BID: 24-h cough frequency estimated
##                 relative reduction vs placebo -18.45% (p=0.041);
##                 15 mg BID not significant; taste-related AE 58% (45),
##                 ~32% (15), 5% placebo; AE discontinuation ~15%.
##   COUGH-2 24wk  gefapixant 45 mg BID: -14.6% vs placebo (p=0.031);
##                 taste-related AE 69%.
##   SOOTHE 4wk    camlipixant 50 mg BID: ~ -34% placebo-adjusted in the
##                 high-baseline stratum; taste AE 6.5% vs 3.2% placebo.
##   PAGANINI 12wk eliapixant 25/75/150 mg BID: -3.5 / -17.6 / -14.2%
##                 vs placebo; taste AE ~5 / 10 / 21%; programme
##                 terminated for idiosyncratic hepatotoxicity.
##   Sivopixant    phase 2b 150 mg QD 12 wk: -18.6% vs placebo, p=0.29
##                 (primary NOT met); taste AE ~6-13%.
##   Gabapentin    Ryan 2012 Lancet, 1800 mg/d 10 wk: LCQ +1.80 vs
##                 placebo, cough VAS -12.1 mm, cough frequency -27%.
##   Morphine SR   Morice 2007 AJRCCM, 5-10 mg BID: daily cough score
##                 -40%, LCQ +3.2.
##   PSALTI        Chamberlain Mitchell 2017 Thorax, speech therapy:
##                 LCQ +1.53 vs control, cough frequency -41%.
##   Pregabalin    Vertigan 2016 Chest, pregabalin + speech therapy:
##                 LCQ +3.5 (vs +1.0 placebo + speech therapy).
##   Nalbuphine ER CANAL, IPF cough: daytime cough -75.7% vs -22.2%.
##   ICS           null in unselected RCC; large in eosinophilic
##                 bronchitis / cough-variant asthma.
##   PPI           null without acid reflux (Cochrane); this is
##                 reproduced, not asserted.
##   ACEi          cough resolves 1-4 weeks after withdrawal.
##
## ---------------------------------------------------------------------
##  HONEST LIMITATIONS
##
##   * Cough counts are log-normally distributed between patients and
##     this is a TYPICAL-SUBJECT model: it predicts geometric-mean
##     trajectories, not the responder distribution. Between-subject
##     variability is offered as parameter sets, not as an OMEGA block.
##   * Adherence is applied as a population-average exposure multiplier.
##     That is the right approximation for an ITT geometric mean and the
##     wrong one for an individual patient.
##   * The selectivity ratios, IC50 values and unbound fractions of the
##     four P2X3 antagonists are collated from heterogeneous preclinical
##     reports and are the least certain parameters in the model. The
##     qualitative window result is robust to them; the exact predicted
##     taste-AE percentages are not.
##   * Eliapixant hepatotoxicity is idiosyncratic and is included as a
##     flagged event, NOT as a mechanistic dose-response.
##   * TWO KNOWN MISFITS, stated rather than tuned away. The model was
##     anchored on COUGH-1 at 12 weeks (-21.0% modelled vs -18.5%
##     observed). Against the other two P2X3 anchors it is wrong in
##     OPPOSITE directions:
##         COUGH-2 gefapixant 45 mg at 24 wk   observed -14.6, model -25.8
##         SOOTHE camlipixant 50 mg at 4 wk    observed  -34,   model -10.5
##     i.e. the modelled treatment effect BUILDS TOO SLOWLY and then
##     KEEPS GROWING, where the real class reaches most of its effect by
##     about 4 weeks and then plateaus or fades. The cause is structural:
##     sustained peripheral blockade here unwinds wind-up, P2X3 density
##     and the epithelial/ATP loop over months, and nothing in the model
##     opposes that unwinding. A tolerance or receptor-escape mechanism
##     is the obvious missing piece. Parameter sweeps over the P2X3 drive
##     share and the wind-up decay rate (see the scan behind KWUIN and
##     WP2X) move both numbers together and cannot fix one without
##     worsening the other, so this is a MISSING MECHANISM, not a
##     mis-set parameter. Treat modelled P2X3 effects beyond ~12 weeks
##     as unreliable.
##   * The gefapixant 15 mg dose is reproduced as clearly active
##     (-13.9%) where the trials found it indistinguishable from placebo.
##     (model -13.0 vs observed -1.4, ns). No occupancy-based model
##     reproduces that cliff: 15 mg BID achieves substantial P2X3
##     occupancy (0.44 here). This is an open problem in the field,
##     not a defect specific to this implementation, and it is left
##     visible rather than fitted away.
##   * The modelled ATP-challenge shift under P2X3 blockade (1.13-1.20
##     fold) is far smaller than the near-log-order shifts reported
##     experimentally, even though its DIRECTION and its dissociation
##     from capsaicin C5 are right. The challenge outputs should be read
##     qualitatively.
##   * Speech therapy reproduces its LCQ effect (+1.12 vs +1.53) but not
##     its objective cough-frequency effect (-16.8% vs -41%).
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

code <- '
$PARAM @annotated
// ================== insult, susceptibility, phenotype ==================
INSULT0 : 1.00 : Initial sensitising insult amplitude (0 = healthy control)
KINSU   : 0.055: Insult resolution rate (1/d; t1/2 ~ 12.6 d)
LGAIN   : 1.00 : LOOP GAIN multiplier on all positive-feedback terms (the bistability knob)
SEXF    : 1.0  : 1 = female (lower cough threshold, lower capsaicin C5), 0 = male
AGEF    : 1.0  : Age factor on central gain (1.0 at 60 y)
DURYR   : 5.0  : Symptom duration at enrolment (y) - used only for reporting
SETIC   : 1    : 1 = compute the healthy initial condition in $MAIN; 0 = use a supplied state vector

// ================== trial artefacts (the placebo mechanism) ============
TRIALON : 0    : 1 = subject is in a clinical trial (diary, monitoring, contact)
ETRIAL  : 0.30 : Maximal trial-participation gain in cortical suppression capacity
KCORT   : 0.048: Cortical suppression adaptation rate (1/d; t1/2 ~ 14 d)
FLARE0  : 0.10 : Enrolment flare amplitude (regression-to-the-mean excursion); with ETRIAL this yields a ~32% placebo-arm fall in 24-h cough frequency at week 12
KFLARE  : 0.017: Flare decay rate (1/d; t1/2 ~ 41 d)

// ================== airway epithelium and purinergic signalling ========
KEPID   : 0.055: Epithelial repair rate (1/d)
FDAMI   : 0.30 : Epithelial damage per unit insult
FDAMC   : 0.46 : Maximal epithelial damage from cough-induced mechanical trauma (saturating)
KCDAM   : 3.5  : Relative cough excess giving half-maximal self-damage (saturation constant)
CTOL    : 0.10 : Dead band: cough excess tolerated before any self-damage
FDAME   : 0.16 : Epithelial damage per unit airway eosinophil burden (MBP/ECP)
FDAMP   : 0.18 : Epithelial damage per unit pepsin/bile exposure
KATPT   : 6.0  : Airway ATP turnover rate (1/d)
ATP0    : 0.50 : Healthy airway-surface ATP concentration (uM)
FATPD   : 1.20 : ATP release per unit epithelial damage (pannexin-1 + lytic)
FATPC   : 0.55 : ATP release per unit excess cough (shear stress)
KMUC    : 1.2  : Mucus turnover rate (1/d)
FMUCS   : 0.85 : Mucus accumulation during sleep (relative)
FMUCC   : 0.30 : Mucus clearance boost per unit relative cough frequency
FMUCI   : 0.60 : Mucus production per unit airway inflammation

// ================== airway inflammation ================================
KINFL   : 0.10 : Inflammation turnover rate (1/d)
FINFI   : 0.70 : Inflammation driven by insult
FINFE   : 0.55 : Inflammation driven by eosinophils
EOS0    : 1.0  : Baseline airway eosinophil index (1 = normal, 4+ = NAEB/CVA)
KEOS    : 0.12 : Eosinophil turnover rate (1/d)
EOSDRV  : 1.0  : Th2/ILC2 drive setting the eosinophil set point
KROS    : 0.9  : Oxidative stress turnover (1/d)
FROSI   : 0.8  : ROS per unit inflammation
SMOKE   : 0    : Current tobacco smoke / irritant exposure (0-1)

// ================== reflux =============================================
ACIDSET : 1.0  : Acid reflux burden set point (1 = normal, 3 = pathological)
NACSET  : 1.0  : Non-acid / gaseous reflux burden set point
KREFL   : 0.5  : Reflux burden turnover (1/d)
FCGHREF : 0.20 : Cough-induced increase in reflux (abdominal pressure loop)
KPEPS   : 0.35 : Pepsin/bile laryngeal exposure turnover (1/d)
FPEPA   : 0.60 : Pepsin exposure per unit acid reflux
FPEPN   : 0.45 : Pepsin exposure per unit non-acid reflux

// ================== mediators and neurotrophins ========================
KNGF    : 0.09 : NGF turnover rate (1/d)
FNGFI   : 0.85 : NGF induction per unit inflammation
FNGFD   : 0.70 : NGF induction per unit epithelial damage
FNGFE   : 0.25 : NGF induction by oestrogen withdrawal (menopause)
ESTRO   : 1.0  : Oestrogen status (1 = premenopausal, 0.3 = postmenopausal)
KPGE    : 1.5  : PGE2 turnover rate (1/d)
FPGEI   : 0.90 : PGE2 per unit inflammation
KBK     : 2.5  : Bradykinin turnover rate (1/d)
ACEI    : 0    : ACE inhibitor on board (0/1)
FBKACE  : 1.80 : Bradykinin accumulation on ACE inhibitor

// ================== receptors and peripheral sensitisation =============
KP2XE   : 0.055: P2X3 receptor density turnover (1/d; t1/2 ~ 12.6 d) - THE SLOW STATE
EMXP2X  : 0.75 : Maximal NGF-driven P2X3 up-regulation (fold above baseline)
EC50NG  : 1.10 : NGF above baseline giving half-maximal P2X3 up-regulation
KTRPV   : 0.30 : TRPV1 sensitisation turnover (1/d)
FTVNGF  : 0.40 : TRPV1 sensitisation per unit NGF above baseline
FTVPGE  : 0.35 : TRPV1 sensitisation per unit PGE2 above baseline
FTVBK   : 0.30 : TRPV1 sensitisation per unit bradykinin above baseline
FTVEST  : 0.18 : TRPV1 sensitisation on oestrogen withdrawal
KTRPA   : 0.30 : TRPA1 sensitisation turnover (1/d)
FTAROS  : 0.55 : TRPA1 sensitisation per unit ROS above baseline
KNAVS   : 0.10 : Nav1.7/1.8 up-regulation turnover (1/d)
FNVNGF  : 0.45 : Nav up-regulation per unit NGF above baseline
KPSEN   : 0.40 : Peripheral sensitisation index turnover (1/d)
WPSP2X  : 0.30 : Weight of P2X3 density in peripheral sensitisation
WPSTRV  : 0.28 : Weight of TRPV1 state in peripheral sensitisation
WPSTRA  : 0.14 : Weight of TRPA1 state in peripheral sensitisation
WPSNAV  : 0.20 : Weight of Nav state in peripheral sensitisation

// ================== afferent drive composition =========================
KAFIR   : 8.0  : Afferent firing equilibration rate (1/d; fast relative to disease)
WP2X    : 0.047: FRACTION OF EVOKED TUSSIVE DRIVE CARRIED BY P2X3 - THE CLASS CEILING. Fitted to the COUGH-1 phase 3 endpoint; see the limitations note on why it has to be this small
WTRP    : 0.420: Fraction carried by TRPV1/TRPA1
WMEC    : 0.357: Fraction carried by mechanical afferents (mucus, traction, distension)
WACI    : 0.176: Fraction carried by acid/ASIC afferents
KMATP   : 0.90 : ATP concentration giving half-maximal P2X3 channel activation (uM)
FTRACT  : 0.0  : Parenchymal traction on mechanoreceptors (0 normal, 1.2 in IPF)
FUACS   : 0.0  : Extrapulmonary (nasal/pharyngeal) afferent drive
KECTO   : 0.055: Ectopic discharge turnover (1/d)
KECIN   : 0.20 : Maximal ectopic discharge generation rate
ECTHR   : 1.32 : Total-gain threshold above which ectopic discharge starts
ECKM    : 0.75 : Half-saturation of ectopic discharge generation
HEC     : 2.0  : Hill coefficient for ectopic discharge generation

// ================== central sensitisation ==============================
KSPC    : 1.2  : Central substance P turnover (1/d)
FSPAF   : 0.80 : Central substance P release per unit afferent firing above baseline
KWUIN   : 0.014: Maximal nTS wind-up generation rate (1/d)
KWUOUT  : 0.030: Wind-up decay rate (1/d; t1/2 ~ 23 d - SLOW, this is the memory)
AFTHR   : 1.00 : Afferent firing threshold for wind-up induction (= the healthy point, so induction is zero at health and GRADED above it)
WUKM    : 0.45 : Half-saturation of wind-up induction
HWU     : 1.3  : Hill coefficient of wind-up induction (shallow: a small peripheral cut must not abolish maintenance)
FSPWU   : 0.45 : Substance P/NK-1 facilitation of wind-up
FMICWU  : 0.55 : Microglial facilitation of wind-up (BDNF)
KMICI   : 0.020: Microglial activation rate (1/d)
KMICO   : 0.030: Microglial deactivation rate (1/d)
MICTHR  : 0.08 : Wind-up level above which microglia are recruited (LOW: microglial memory is active in the typical RCC patient and resists unwinding)
MICKM   : 0.50 : Wind-up excess giving half-maximal microglial recruitment (saturation)
KGABI   : 0.06 : nTS inhibitory interneuron tone turnover (1/d)
FGBMIC  : 0.45 : Loss of GABA/glycine tone per unit microglial activation (KCC2 down)
KDESC   : 0.10 : Descending inhibitory tone turnover (1/d)
FDSLP   : 0.12 : Loss of descending inhibition per unit sleep debt
FDANX   : 0.10 : Loss of descending inhibition per unit hypervigilance above baseline
WWU     : 0.85 : Weight of wind-up in central gain
WGAB    : 0.55 : Exponent of inhibitory tone in central gain
WDES    : 0.50 : Exponent of descending tone in central gain

// ================== cortical layer =====================================
KHVIG   : 0.045: Hypervigilance / catastrophising turnover (1/d)
FHVURG  : 0.060: Hypervigilance driven by urge to cough
FHVCMP  : 0.12 : Hypervigilance driven by complications and social impact
CORTB   : 1.0  : Healthy cortical suppression capacity
FCHVIG  : 0.18 : Loss of suppression capacity per unit hypervigilance above baseline
FCSLP   : 0.10 : Loss of suppression capacity per unit sleep debt
ESLT    : 0    : Speech and language therapy effect on suppression capacity (0-0.50)
ECBT    : 0    : CBT effect reducing hypervigilance (0-0.5)
KHABT   : 0.05 : Habitual/learned cough turnover (1/d)
FHBURG  : 0.040: Habit formation per unit urge
KSLPD   : 0.25 : Sleep debt turnover (1/d)
FSLPC   : 0.012: Sleep debt accrued per cough/h of awake cough frequency

// ================== cough generation ===================================
CFMAX   : 130  : Maximal sustainable cough rate (coughs/h)
THR0    : 4.95 : Intrinsic cough threshold (drive units)
HC      : 2.60 : Hill coefficient of the cough threshold function
WCORT   : 0.55 : Weight of cortical suppression on the cough threshold
WDESCT  : 0.45 : Weight of descending inhibition on the cough threshold
FTHRSEX : 0.88 : Cough threshold multiplier in women (lower threshold)
SLEEPRS : 0.06 : Residual cough rate fraction during sleep
KCFQ    : 40   : Cough-frequency state equilibration rate (1/d)
KWAKE   : 22   : Steepness of the wake/sleep gate
AWKCUT  : -0.50: Circadian cosine value at the wake/sleep transition (16 h awake)
OMEGA   : 6.2832 : Circadian angular frequency (rad/d)

// ================== urge to cough and PROs =============================
THRU    : 0.50 : Urge threshold as a fraction of the cough threshold (urge comes first)
HU      : 2.6  : Hill coefficient of the urge function
KURG    : 20   : Urge state equilibration rate (1/d)
KLCQ    : 0.14 : LCQ equilibration rate (1/d; ~2-week recall period)
ALCQ    : 4.05 : LCQ decrement per log-unit of 24-h cough frequency
BLCQ    : 0.62 : LCQ decrement per unit urge score
CLCQ    : 1.15 : LCQ decrement per unit hypervigilance above baseline
DLCQ    : 1.30 : LCQ decrement per unit complication burden
KVAS    : 0.35 : Cough VAS equilibration rate (1/d)
AVAS    : 58   : VAS contribution from urge (mm)
BVAS    : 52   : VAS contribution from cough frequency (mm)
KCSD    : 0.30 : Cough Severity Diary equilibration rate (1/d)
KCOMP   : 0.02 : Complication burden turnover (1/d)
FCMPF   : 1.35 : Complication susceptibility multiplier in women (stress incontinence)

// ================== challenge tests ====================================
C5M0    : 62   : Capsaicin C5 in healthy men (uM)
C5F0    : 26   : Capsaicin C5 in healthy women (uM)
ATPC50  : 220  : ATP challenge threshold concentration in health (uM)
WC5CEN  : 0.55 : Exponent of central gain on capsaicin C5 (central contribution)

// ================== P2X3 antagonist pharmacology =======================
KATPC   : 1.00 : Reference ATP concentration for the competitive right-shift (uM)
EMAX3   : 0.95 : Maximal functional inhibition of the P2X3 arm at full occupancy
TDG50   : 0.45 : P2X2/3 occupancy giving half-maximal dysgeusia - NO RECEPTOR RESERVE
HDG     : 2.0  : Hill coefficient of the dysgeusia occupancy function
DGMAX   : 10   : Maximal dysgeusia score
KDYSG   : 1.4  : Dysgeusia equilibration rate (1/d)
PBOAE   : 3.0  : Placebo-arm taste adverse-event rate (%)
AE50    : 4.5  : Dysgeusia score giving half-maximal taste-AE reporting
HAE     : 1.6  : Hill coefficient of taste-AE reporting
DGTOL   : 2.0  : Dysgeusia score tolerated before discontinuation begins
KDROP   : 0.0005 : Discontinuation hazard per unit intolerable dysgeusia (1/d); calibrated to ~15-20% taste-related withdrawal by week 12 at gefapixant 45 mg BID

// -------- gefapixant (MK-7264) -----------------------------------------
GEF_KA  : 36   : Gefapixant absorption rate (1/d; ka 1.5 /h)
GEF_CL  : 134  : Gefapixant apparent clearance (L/d; 5.6 L/h)
GEF_V2  : 55   : Gefapixant central volume (L)
GEF_Q   : 72   : Gefapixant intercompartmental clearance (L/d)
GEF_V3  : 25   : Gefapixant peripheral volume (L)
GEF_FU  : 0.30 : Gefapixant unbound fraction
GEF_MW  : 384.4: Gefapixant molecular weight (g/mol)
GEF_IC3 : 0.10 : Gefapixant P2X3 functional IC50 at reference ATP (uM)
GEF_SEL : 6.0  : Gefapixant P2X2/3 : P2X3 selectivity ratio
GEF_F   : 0.90 : Gefapixant oral bioavailability
GEF_RF  : 1.0  : Renal function multiplier on gefapixant clearance (OAT3)

// -------- camlipixant (BLU-5937) ---------------------------------------
CAM_KA  : 24   : Camlipixant absorption rate (1/d)
CAM_CL  : 125  : Camlipixant apparent clearance (L/d)
CAM_V2  : 90   : Camlipixant central volume (L)
CAM_FU  : 0.15 : Camlipixant unbound fraction
CAM_MW  : 402  : Camlipixant molecular weight (g/mol)
CAM_IC3 : 0.030: Camlipixant P2X3 functional IC50 at reference ATP (uM)
CAM_SEL : 1500 : Camlipixant P2X2/3 : P2X3 selectivity ratio
CAM_F   : 0.75 : Camlipixant oral bioavailability

// -------- eliapixant (BAY 1817080) -------------------------------------
ELI_KA  : 18   : Eliapixant absorption rate (1/d)
ELI_CL  : 110  : Eliapixant apparent clearance (L/d)
ELI_V2  : 100  : Eliapixant central volume (L)
ELI_FU  : 0.10 : Eliapixant unbound fraction
ELI_MW  : 428  : Eliapixant molecular weight (g/mol)
ELI_IC3 : 0.020: Eliapixant P2X3 functional IC50 at reference ATP (uM)
ELI_SEL : 70   : Eliapixant P2X2/3 : P2X3 selectivity ratio
ELI_F   : 0.70 : Eliapixant oral bioavailability
ELI_DIL : 0    : Eliapixant idiosyncratic DILI susceptibility flag (0/1)

// -------- sivopixant (S-600918) ----------------------------------------
SIV_KA  : 12   : Sivopixant absorption rate (1/d)
SIV_CL  : 95   : Sivopixant apparent clearance (L/d)
SIV_V2  : 120  : Sivopixant central volume (L)
SIV_FU  : 0.03 : Sivopixant unbound fraction (highly protein bound)
SIV_MW  : 480  : Sivopixant molecular weight (g/mol)
SIV_IC3 : 0.0045 : Sivopixant P2X3 functional IC50 at reference ATP (uM)
SIV_SEL : 40   : Sivopixant P2X2/3 : P2X3 selectivity ratio
SIV_F   : 0.60 : Sivopixant oral bioavailability

// -------- gabapentin ----------------------------------------------------
GAB_KA  : 18   : Gabapentin absorption rate (1/d)
GAB_CL  : 274  : Gabapentin clearance (L/d; 11.4 L/h, renal)
GAB_V2  : 58   : Gabapentin volume (L)
FA2DPRE : 0.55 : Gabapentinoid presynaptic (transmitter-release) share of the alpha2delta effect
GAB_EMX : 0.42 : Gabapentin maximal inhibition of wind-up induction (alpha2delta-1)
GAB_EC50: 4.0  : Gabapentin concentration for half-maximal effect (mg/L)
GAB_AE  : 0.31 : Gabapentin dizziness/somnolence incidence at 1800 mg/d

// -------- pregabalin ----------------------------------------------------
PRE_KA  : 36   : Pregabalin absorption rate (1/d)
PRE_CL  : 120  : Pregabalin clearance (L/d)
PRE_V2  : 42   : Pregabalin volume (L)
PRE_EMX : 0.60 : Pregabalin maximal inhibition of wind-up induction
PRE_EC50: 2.6  : Pregabalin concentration for half-maximal effect (mg/L)

// -------- morphine sulfate SR -------------------------------------------
MOR_KA  : 3.6  : Morphine SR release rate (1/d; slow-release depot)
MOR_CL  : 1500 : Morphine clearance (L/d)
MOR_V2  : 230  : Morphine volume (L)
MOR_EMX : 0.65 : Morphine maximal enhancement of descending inhibition
MOR_EC50: 0.011: Morphine concentration for half-maximal effect (mg/L)
MOR_NTS : 0.30 : Morphine direct inhibition of nTS relay at maximal effect

// -------- nalbuphine ER -------------------------------------------------
NAL_KA  : 6.0  : Nalbuphine ER release rate (1/d)
NAL_CL  : 1900 : Nalbuphine clearance (L/d)
NAL_V2  : 250  : Nalbuphine volume (L)
NAL_EMX : 0.95 : Nalbuphine maximal enhancement of descending inhibition (kappa)
NAL_EC50: 0.020: Nalbuphine concentration for half-maximal effect (mg/L)

// -------- amitriptyline -------------------------------------------------
AMI_KA  : 12   : Amitriptyline absorption rate (1/d)
AMI_CL  : 800  : Amitriptyline clearance (L/d)
AMI_V2  : 1100 : Amitriptyline volume (L)
AMI_EMX : 0.30 : Amitriptyline maximal enhancement of descending inhibition
AMI_EC50: 0.030: Amitriptyline concentration for half-maximal effect (mg/L)

// -------- inhaled corticosteroid ----------------------------------------
ICS_KA  : 8.0  : ICS lung-to-systemic transfer rate (1/d)
ICS_CL  : 1600 : ICS systemic clearance (L/d)
ICS_V2  : 300  : ICS systemic volume (L)
ICS_EMX : 0.80 : ICS maximal suppression of the eosinophil set point
ICS_EC50: 30   : ICS lung amount for half-maximal effect (ug)
ICS_NFK : 0.00 : ICS maximal direct suppression of airway inflammation (eosinophil-independent component is small)

// -------- proton pump inhibitor -----------------------------------------
PPI_KA  : 24   : PPI absorption rate (1/d)
PPI_CL  : 400  : PPI clearance (L/d)
PPI_V2  : 16   : PPI volume (L)
PPI_EMX : 0.80 : PPI maximal suppression of ACID reflux burden (no effect on non-acid)
PPI_EC50: 0.15 : PPI concentration for half-maximal acid suppression (mg/L)

// -------- dextromethorphan (weak comparator) ----------------------------
DXM_KA  : 30   : Dextromethorphan absorption rate (1/d)
DXM_CL  : 2400 : Dextromethorphan clearance (L/d)
DXM_V2  : 400  : Dextromethorphan volume (L)
DXM_EMX : 0.10 : Dextromethorphan maximal inhibition of wind-up (weak in RCC)
DXM_EC50: 0.01 : Dextromethorphan concentration for half-maximal effect (mg/L)

// -------- non-pharmacological / procedural ------------------------------
ESLNB   : 0    : Superior laryngeal nerve block effect on afferent drive (0-0.5)
ENEBLD  : 0    : Nebulised lidocaine effect on Nav-dependent firing (0-0.6)
EMENTH  : 0    : Menthol/TRPM8 anti-tussive effect (0-0.25)
EANTIH  : 0    : Antihistamine effect on extrapulmonary afferent drive (0-0.7)
EAZI    : 0    : Macrolide effect on neutrophilic inflammation (0-0.4)

$CMT @annotated
CIRCS : Circadian oscillator sine state (-)
CIRCC : Circadian oscillator cosine state (-)
INSU  : Sensitising insult amplitude (-)
FLARE : Enrolment flare (regression-to-the-mean excursion) (-)
EPID  : Epithelial integrity (1 = intact)
ATPX  : Airway-surface ATP concentration (uM)
MUC   : Airway mucus burden (-)
INFL  : Airway inflammation index (-)
EOSN  : Airway eosinophil index (-)
ROSX  : Airway oxidative stress index (-)
ACIDR : Acid reflux burden (-)
NACID : Non-acid / gaseous reflux burden (-)
PEPS  : Laryngeal pepsin / bile exposure (-)
NGF   : Airway nerve growth factor (-)
PGE2  : Airway PGE2 (-)
BKN   : Airway bradykinin (-)
P2X3E : P2X3 receptor density on vagal afferents (-)
TRPV  : TRPV1 sensitisation state (-)
TRPA  : TRPA1 sensitisation state (-)
NAVS  : Nav1.7/1.8 up-regulation state (-)
PSEN  : Peripheral sensitisation index (-)
ECTO  : Ectopic spontaneous afferent discharge (-)
AFIR  : Integrated vagal afferent firing (-)
SPC   : Central substance P / NK-1 tone (-)
WIND  : nTS wind-up / central sensitisation (-)
MICG  : nTS microglial activation (-)
GABI  : nTS GABA/glycinergic inhibitory tone (-)
DESC  : Descending inhibitory tone (-)
HVIG  : Hypervigilance / cough catastrophising (-)
CORT  : Cortical suppression capacity (-)
HABT  : Learned / habitual cough behaviour (-)
SLPD  : Sleep debt (-)
CFQ   : Instantaneous cough frequency (coughs/h)
URG   : Urge-to-cough score (0-10)
LCQ   : Leicester Cough Questionnaire total (3-21)
VASC  : Cough severity VAS (0-100 mm)
CSDS  : Cough Severity Diary score (0-10)
CACC  : Cumulative cough count (coughs)
CAWK  : Cumulative awake cough count (coughs)
TAWK  : Cumulative awake time (h)
COMP  : Complication burden (incontinence, syncope, rib fracture) (-)
DYSG  : Dysgeusia score (0-10)
ADHR  : Adherence / proportion still on treatment (-)
ALTX  : ALT (U/L)
GEFD  : Gefapixant gut depot (mg)
GEFC  : Gefapixant central (mg)
GEFP  : Gefapixant peripheral (mg)
CAMD  : Camlipixant gut depot (mg)
CAMC  : Camlipixant central (mg)
ELID  : Eliapixant gut depot (mg)
ELIC  : Eliapixant central (mg)
SIVD  : Sivopixant gut depot (mg)
SIVC  : Sivopixant central (mg)
GABD  : Gabapentin gut depot (mg)
GABC  : Gabapentin central (mg)
PRED  : Pregabalin gut depot (mg)
PREC  : Pregabalin central (mg)
MORD  : Morphine SR depot (mg)
MORC  : Morphine central (mg)
NALD  : Nalbuphine ER depot (mg)
NALC  : Nalbuphine central (mg)
AMID  : Amitriptyline gut depot (mg)
AMIC  : Amitriptyline central (mg)
ICSL  : Inhaled corticosteroid lung depot (ug)
ICSC  : Inhaled corticosteroid systemic (ug)
PPID  : PPI gut depot (mg)
PPIC  : PPI central (mg)
DXMD  : Dextromethorphan gut depot (mg)
DXMC  : Dextromethorphan central (mg)

$GLOBAL
#define HILL(x, k, h) (pow((x), (h)) / (pow((k), (h)) + pow((x), (h))))
// threshold-Hill: exactly zero below the threshold, so healthy states are
// exact fixed points rather than approximate ones
#define THILL(x, thr, k, h) (((x) <= (thr)) ? 0.0 : HILL((x) - (thr), (k), (h)))
#define POS(x) (((x) > 0.0) ? (x) : 0.0)
// SATURATING excess. Every self-damage term in this model is written with
// this rather than as a proportional term. A patient coughing 30 times an
// hour irritates their larynx; they do not denude their airway ten times
// over. Unsaturated self-damage was the difference between a bistable
// model and one with no upper stable state at all.
#define SAT(x, k) ((x) / ((k) + (x)))

// Subject-specific healthy reference point, evaluated once per record in
// $MAIN and reused in $ODE. Referencing the damage terms and the PROs to
// THIS subject healthy state (not to a global constant) is what makes a
// healthy man and a healthy woman both exact fixed points despite having
// different cough thresholds.
namespace {
  double CFREF;    // healthy AWAKE cough rate (coughs/h)
  double CFREF24;  // healthy 24-h mean cough rate (coughs/h)
  double URGREF;   // healthy urge score
  double MUCREF;   // healthy CYCLE-AVERAGE mucus burden
}

$MAIN
// ---- circadian phase: t = 0 is 07:00, the moment of waking ------------
// CIRCC = cos(2*pi*(t - 1/3)) crosses AWKCUT = -0.5 at 07:00 and 23:00,
// giving a 16 h awake window without any explicit clock.
CIRCS_0 = -0.8660254;
CIRCC_0 = -0.5;

// Healthy operating point, evaluated analytically so that a subject with
// INSULT0 = 0 starts exactly on their own fixed point. Awake fraction is
// 16/24 and the residual sleep rate is SLEEPRS.
double thr_h = THR0 * (SEXF > 0.5 ? FTHRSEX : 1.0);
CFREF   = CFMAX * HILL(1.0, thr_h, HC);
CFREF24 = CFREF * (2.0 / 3.0 + SLEEPRS / 3.0);
URGREF  = DGMAX * HILL(1.0, THRU * thr_h, HU);
// Mucus accumulates overnight and is cleared by coughing, which is what
// generates the morning peak. Its cycle AVERAGE (asleep 1/3 of the day) is
// above 1, so the mechanical arm is normalised to that average — otherwise
// a perfectly healthy subject sits permanently ~6% above unit afferent
// drive and slowly walks away from their own fixed point.
MUCREF  = 1.0 + FMUCS / 3.0;

// Every *_0 assignment below is GUARDED. mrgsolve evaluates $MAIN on every
// record, so unguarded initial-condition code silently overwrites any state
// vector supplied through init() — which would have thrown away the whole
// natural-history pre-run and quietly simulated a healthy subject instead.
if (SETIC > 0.5) {
INSU_0  = INSULT0;
FLARE_0 = 0.0;             // set at randomisation by the trial scenarios

EPID_0  = 1.0;
ATPX_0  = ATP0;
MUC_0   = MUCREF;
INFL_0  = 1.0;
EOSN_0  = EOS0;
ROSX_0  = 1.0;
ACIDR_0 = ACIDSET;
NACID_0 = NACSET;
PEPS_0  = 1.0;
NGF_0   = 1.0;
PGE2_0  = 1.0;
BKN_0   = 1.0;
P2X3E_0 = 1.0;
TRPV_0  = 1.0;
TRPA_0  = 1.0;
NAVS_0  = 1.0;
PSEN_0  = 1.0;
ECTO_0  = 0.0;
AFIR_0  = 1.0;
SPC_0   = 1.0;
WIND_0  = 0.0;
MICG_0  = 0.0;
GABI_0  = 1.0;
DESC_0  = 1.0;
HVIG_0  = 1.0;
CORT_0  = CORTB;
HABT_0  = 0.0;
SLPD_0  = 0.0;


CFQ_0  = CFREF * (SLEEPRS + (1.0 - SLEEPRS) * 0.5);  // t = 0 is the wake transition
URG_0  = URGREF;
LCQ_0  = 21.0;    // PROs are referenced to this subject healthy state,
VASC_0 = 0.0;     // so a healthy subject scores a perfect LCQ, zero VAS
CSDS_0 = 0.0;     // and zero CSD by construction rather than by fitting
CACC_0 = 0.0; CAWK_0 = 0.0; TAWK_0 = 0.0; COMP_0 = 0.0;
DYSG_0 = 0.0; ADHR_0 = 1.0; ALTX_0 = 25.0;
}

$ODE
// =====================================================================
//  0. WAKE / SLEEP GATE
// =====================================================================
double AWAKE = 1.0 / (1.0 + exp(-KWAKE * (CIRCC - AWKCUT)));

// =====================================================================
//  1. DRUG CONCENTRATIONS AND TARGET ENGAGEMENT
// =====================================================================
double Cgef = GEFC / GEF_V2;                    // mg/L
double Ccam = CAMC / CAM_V2;
double Celi = ELIC / ELI_V2;
double Csiv = SIVC / SIV_V2;
double Cgab = GABC / GAB_V2;
double Cpre = PREC / PRE_V2;
double Cmor = MORC / MOR_V2;
double Cnal = NALC / NAL_V2;
double Cami = AMIC / AMI_V2;
double Cppi = PPIC / PPI_V2;
double Cdxm = DXMC / DXM_V2;

// Free molar concentrations of the P2X3 antagonists (uM), in TWO versions.
//   *_on  = what an adherent patient is actually exposed to. This drives
//           DYSGEUSIA, because the patients who quit are precisely the ones
//           whose taste loss persisted.
//   *_pop = population-average exposure, eroded by discontinuation. This
//           drives EFFICACY, because that is what an intention-to-treat
//           geometric mean measures.
double Ugef_on = Cgef * GEF_FU / GEF_MW * 1000.0;
double Ucam_on = Ccam * CAM_FU / CAM_MW * 1000.0;
double Ueli_on = Celi * ELI_FU / ELI_MW * 1000.0;
double Usiv_on = Csiv * SIV_FU / SIV_MW * 1000.0;
double Ugef = Ugef_on * ADHR;
double Ucam = Ucam_on * ADHR;
double Ueli = Ueli_on * ADHR;
double Usiv = Usiv_on * ADHR;

// ATP-COMPETITIVE right-shift: a more damaged, higher-ATP airway needs a
// higher antagonist concentration for the same functional inhibition.
double shift = 1.0 + ATPX / KATPC;

double occ3_gef = Ugef / (Ugef + GEF_IC3 * shift);
double occ3_cam = Ucam / (Ucam + CAM_IC3 * shift);
double occ3_eli = Ueli / (Ueli + ELI_IC3 * shift);
double occ3_siv = Usiv / (Usiv + SIV_IC3 * shift);

// combined blockade of the P2X3 arm (independent-action combination)
double free3 = (1.0 - EMAX3 * occ3_gef) * (1.0 - EMAX3 * occ3_cam) *
               (1.0 - EMAX3 * occ3_eli) * (1.0 - EMAX3 * occ3_siv);
double INH3  = 1.0 - free3;

// taste-bud P2X2/3 occupancy — SAME drug, SAME curve, shifted by selectivity
double occ23 = 1.0 - (1.0 - Ugef_on / (Ugef_on + GEF_IC3 * GEF_SEL)) *
                     (1.0 - Ucam_on / (Ucam_on + CAM_IC3 * CAM_SEL)) *
                     (1.0 - Ueli_on / (Ueli_on + ELI_IC3 * ELI_SEL)) *
                     (1.0 - Usiv_on / (Usiv_on + SIV_IC3 * SIV_SEL));

// central neuromodulators
double Egab = GAB_EMX * Cgab / (Cgab + GAB_EC50);
double Epre = PRE_EMX * Cpre / (Cpre + PRE_EC50);
double Edxm = DXM_EMX * Cdxm / (Cdxm + DXM_EC50);
double EWUD = 1.0 - (1.0 - Egab) * (1.0 - Epre) * (1.0 - Edxm);   // wind-up block
double EA2D = FA2DPRE * (1.0 - (1.0 - Egab) * (1.0 - Epre));      // presynaptic release

double Emor = MOR_EMX * Cmor / (Cmor + MOR_EC50);
double Enal = NAL_EMX * Cnal / (Cnal + NAL_EC50);
double Eami = AMI_EMX * Cami / (Cami + AMI_EC50);
double EDES = Emor + Enal + Eami;                                  // descending gain

// treatable-trait drugs
double Eics  = ICS_EMX * ICSL / (ICSL + ICS_EC50);
double Eicsn = ICS_NFK * ICSL / (ICSL + ICS_EC50);
double Eppi  = PPI_EMX * Cppi / (Cppi + PPI_EC50);

// =====================================================================
//  2. INSULT, FLARE, AIRWAY
// =====================================================================
dxdt_INSU  = -KINSU * INSU;
dxdt_FLARE = -KFLARE * FLARE;

// cough burden relative to THIS subject healthy awake rate; the +CTOL
// dead band keeps normal daytime coughing from being self-damaging
double cfex   = SAT(POS(CFQ / fmax(1e-6, CFREF) - 1.0 - CTOL), KCDAM);
double damage = FDAMI * (INSU + FLARE) + FDAMC * cfex
              + FDAME * SAT(POS(EOSN - 1.0), 3.0)
              + FDAMP * SAT(POS(PEPS - 1.0), 1.5)
              + 0.35 * SMOKE;
double epid_t = 1.0 / (1.0 + damage);
dxdt_EPID = KEPID * (epid_t - EPID);

double atp_t = ATP0 * (1.0 + FATPD * POS(1.0 - EPID) / fmax(0.05, EPID)
                           + FATPC * cfex);
dxdt_ATPX = KATPT * (atp_t - ATPX);

double muc_t = (1.0 + FMUCI * POS(INFL - 1.0) + FMUCS * (1.0 - AWAKE))
             / (1.0 + FMUCC * cfex);
dxdt_MUC = KMUC * (muc_t - MUC);

double infl_t = 1.0 + FINFI * (INSU + FLARE) + FINFE * POS(EOSN - 1.0)
              + 0.5 * SMOKE - EAZI;
infl_t = fmax(0.2, infl_t * (1.0 - Eicsn));
dxdt_INFL = KINFL * (infl_t - INFL);

double eos_t = 1.0 + (EOSDRV * EOS0 - 1.0) * (1.0 - Eics);
dxdt_EOSN = KEOS * (eos_t - EOSN);

double ros_t = 1.0 + FROSI * POS(INFL - 1.0) + 0.6 * SMOKE;
dxdt_ROSX = KROS * (ros_t - ROSX);

// =====================================================================
//  3. REFLUX  (PPI suppresses ACID only — this is why PPI trials fail)
// =====================================================================
double acid_t = ACIDSET * (1.0 - Eppi) * (1.0 + FCGHREF * cfex);
dxdt_ACIDR = KREFL * (acid_t - ACIDR);
double nac_t = NACSET * (1.0 + FCGHREF * cfex);
dxdt_NACID = KREFL * (nac_t - NACID);
double peps_t = 1.0 + FPEPA * POS(ACIDR - 1.0) + FPEPN * POS(NACID - 1.0);
dxdt_PEPS = KPEPS * (peps_t - PEPS);

// =====================================================================
//  4. MEDIATORS AND NEUROTROPHINS
// =====================================================================
double ngf_t = 1.0 + FNGFI * POS(INFL - 1.0) + FNGFD * POS(1.0 - EPID)
             + FNGFE * (1.0 - ESTRO);
dxdt_NGF = KNGF * (ngf_t - NGF);

double pge_t = 1.0 + FPGEI * POS(INFL - 1.0);
dxdt_PGE2 = KPGE * (pge_t - PGE2);

double bk_t = 1.0 + FBKACE * ACEI + 0.30 * POS(INFL - 1.0);
dxdt_BKN = KBK * (bk_t - BKN);

// =====================================================================
//  5. PERIPHERAL SENSITISATION  (the SLOW receptor layer)
// =====================================================================
// NGF -> TrkA -> transcription: P2X3 density changes over WEEKS. This is
// why a P2X3 antagonist relieves symptoms in days and modifies nothing.
double p2x_t = 1.0 + EMXP2X * HILL(POS(NGF - 1.0), EC50NG, 1.0);
dxdt_P2X3E = KP2XE * (p2x_t - P2X3E);

double trpv_t = 1.0 + FTVNGF * POS(NGF - 1.0) + FTVPGE * POS(PGE2 - 1.0)
              + FTVBK * POS(BKN - 1.0) + FTVEST * (1.0 - ESTRO);
dxdt_TRPV = KTRPV * (trpv_t - TRPV);

double trpa_t = 1.0 + FTAROS * POS(ROSX - 1.0);
dxdt_TRPA = KTRPA * (trpa_t - TRPA);

double navs_t = (1.0 + FNVNGF * POS(NGF - 1.0)) * (1.0 - ENEBLD);
dxdt_NAVS = KNAVS * (navs_t - NAVS);

double psen_t = 1.0 + WPSP2X * POS(P2X3E - 1.0) + WPSTRV * POS(TRPV - 1.0)
              + WPSTRA * POS(TRPA - 1.0) + WPSNAV * POS(NAVS - 1.0);
dxdt_PSEN = KPSEN * (psen_t - PSEN);

// =====================================================================
//  6. AFFERENT DRIVE
//     Four parallel arms whose weights sum to 1 at health, so AFIR = 1
//     in a healthy subject by construction. WP2X is the fraction of the
//     evoked drive a P2X3 antagonist can reach — the CLASS CEILING.
// =====================================================================
double S_p2x = P2X3E * (ATPX / (KMATP + ATPX)) / (ATP0 / (KMATP + ATP0))
             * (1.0 - INH3);
double S_trp = (0.62 * TRPV + 0.38 * TRPA) * (1.0 - EMENTH);
double S_mec = MUC / MUCREF + FTRACT + FUACS * (1.0 - EANTIH);
double S_aci = 0.55 * ACIDR + 0.25 * NACID + 0.20 * PEPS;

double evoked = WP2X * S_p2x + WTRP * S_trp + WMEC * S_mec + WACI * S_aci;
double afir_t = (PSEN * evoked * NAVS + ECTO) * (1.0 - ESLNB);
dxdt_AFIR = KAFIR * (afir_t - AFIR);

// =====================================================================
//  7. CENTRAL SENSITISATION  (LAYER 2 — the memory of the disease)
// =====================================================================
double spc_t = 1.0 + FSPAF * POS(AFIR - 1.0);
dxdt_SPC = KSPC * (spc_t - SPC);

// Threshold-Hill: production is EXACTLY zero at healthy afferent firing,
// and cooperative (HWU = 2.4) above it. This nonlinearity is what makes
// the disease bistable.
double wu_in = LGAIN * KWUIN * THILL(AFIR, AFTHR, WUKM, HWU)
             * (1.0 + FSPWU * POS(SPC - 1.0))
             * (1.0 + FMICWU * MICG)
             * (1.0 - EWUD);
dxdt_WIND = wu_in - KWUOUT * WIND;

// SATURATING recruitment: microglial activation is bounded, so the
// wind-up -> BDNF -> disinhibition -> wind-up loop has an upper stable
// state instead of diverging.
dxdt_MICG = KMICI * SAT(POS(WIND - MICTHR), MICKM) - KMICO * MICG;

double gabi_t = 1.0 / (1.0 + FGBMIC * MICG);
dxdt_GABI = KGABI * (gabi_t - GABI);

double desc_t = (1.0 + EDES) / (1.0 + FDSLP * SLPD + FDANX * POS(HVIG - 1.0));
dxdt_DESC = KDESC * (desc_t - DESC);

// central gain: wind-up amplifies, inhibitory and descending tone divide
double CENG = AGEF * (1.0 + WWU * WIND) / (pow(fmax(0.05, GABI), WGAB)
                                         * pow(fmax(0.05, DESC), WDES));
CENG = CENG * (1.0 - MOR_NTS * Emor / fmax(1e-9, MOR_EMX));

// ectopic discharge — cough with a completely normal airway
double ect_t = KECIN * THILL(PSEN * CENG, ECTHR, ECKM, HEC) * LGAIN;
dxdt_ECTO = KECTO * (ect_t - ECTO);

// =====================================================================
//  8. CORTICAL LAYER  (LAYER 3 — and the placebo mechanism)
// =====================================================================
double hvig_t = 1.0 + FHVURG * POS(URG - URGREF) + FHVCMP * COMP;
hvig_t = 1.0 + (hvig_t - 1.0) * (1.0 - ECBT);
dxdt_HVIG = KHVIG * (hvig_t - HVIG);

double cort_t = CORTB * (1.0 + ETRIAL * TRIALON + ESLT)
              / (1.0 + FCHVIG * POS(HVIG - 1.0) + FCSLP * SLPD);
dxdt_CORT = KCORT * (cort_t - CORT);

double habt_t = FHBURG * POS(URG - URGREF);
dxdt_HABT = KHABT * (habt_t - HABT);

// =====================================================================
//  9. COUGH GENERATION
//     drive vs threshold, gated by wakefulness
// =====================================================================
double TUSS = AFIR * CENG * (1.0 + HABT) * (1.0 - EA2D);
double THRS = THR0 * (SEXF > 0.5 ? FTHRSEX : 1.0)
            * (1.0 + WCORT * (CORT - CORTB))
            * (1.0 + WDESCT * (DESC - 1.0));
THRS = fmax(0.4, THRS);

double cfq_t = CFMAX * HILL(TUSS, THRS, HC) * (SLEEPRS + (1.0 - SLEEPRS) * AWAKE);
dxdt_CFQ = KCFQ * (cfq_t - CFQ);

// the URGE has a LOWER threshold than the motor cough, which is why a
// patient can suppress the cough and still feel awful — and why PROs and
// cough counts dissociate.
double urg_t = DGMAX * HILL(TUSS, THRU * THRS, HU);
dxdt_URG = KURG * (urg_t - URG);

// =====================================================================
// 10. SLEEP, ACCUMULATORS, COMPLICATIONS
// =====================================================================
double slpd_t = FSLPC * POS(CFQ - CFREF) * AWAKE + 0.25 * POS(HVIG - 1.0);
dxdt_SLPD = KSLPD * (slpd_t - SLPD);

dxdt_CACC = CFQ * 24.0;                 // coughs per day
dxdt_CAWK = CFQ * 24.0 * AWAKE;
dxdt_TAWK = 24.0 * AWAKE;

double comp_t = FCMPF * (SEXF > 0.5 ? 1.0 : 0.55) * HILL(CFQ, 45.0, 1.6);
dxdt_COMP = KCOMP * (comp_t - COMP);

// =====================================================================
// 11. PATIENT-REPORTED OUTCOMES  (2-week recall -> they LAG the counts)
// =====================================================================
// All three PROs are referenced to THIS subject healthy state, so a
// healthy subject scores LCQ 21 / VAS 0 / CSD 0 exactly. Each is driven by
// BOTH the cough count and the urge, and the urge has the lower threshold
// — which is the structural reason PROs and cough counts dissociate.
double cf24  = CFQ;                      // instantaneous; daily means in $TABLE
double lcq_t = 21.0
             - ALCQ * log((1.0 + cf24 / 2.0) / (1.0 + CFREF24 / 2.0))
             - BLCQ * POS(URG - URGREF)
             - CLCQ * POS(HVIG - 1.0) - DLCQ * COMP;
lcq_t = fmin(21.0, fmax(3.0, lcq_t));
dxdt_LCQ = KLCQ * (lcq_t - LCQ);

double vas_t = fmin(100.0, AVAS * POS(URG - URGREF) / 10.0
                         + BVAS * POS(cf24 - CFREF24) / 40.0);
dxdt_VASC = KVAS * (vas_t - VASC);

double csd_t = fmin(10.0, 5.0 * POS(URG - URGREF) / 10.0
                        + 5.0 * POS(cf24 - CFREF24) / 40.0);
dxdt_CSDS = KCSD * (csd_t - CSDS);

// =====================================================================
// 12. ON-TARGET TOLERABILITY — the other half of the occupancy curve
// =====================================================================
double dysg_t = DGMAX * HILL(occ23, TDG50, HDG);
dxdt_DYSG = KDYSG * (dysg_t - DYSG);

// discontinuation hazard: patients leave when taste loss exceeds tolerance
dxdt_ADHR = -KDROP * ADHR * POS(DYSG - DGTOL);

// eliapixant DILI: idiosyncratic, flagged rather than dose-modelled
dxdt_ALTX = 0.08 * (25.0 * (1.0 + 6.0 * ELI_DIL * (ELIC > 0.01 ? 1.0 : 0.0)) - ALTX);

// =====================================================================
// 13. DRUG PK
// =====================================================================
dxdt_GEFD = -GEF_KA * GEFD;
dxdt_GEFC =  GEF_KA * GEFD - (GEF_CL * GEF_RF / GEF_V2) * GEFC
           - (GEF_Q / GEF_V2) * GEFC + (GEF_Q / GEF_V3) * GEFP;
dxdt_GEFP =  (GEF_Q / GEF_V2) * GEFC - (GEF_Q / GEF_V3) * GEFP;

dxdt_CAMD = -CAM_KA * CAMD;
dxdt_CAMC =  CAM_KA * CAMD - (CAM_CL / CAM_V2) * CAMC;
dxdt_ELID = -ELI_KA * ELID;
dxdt_ELIC =  ELI_KA * ELID - (ELI_CL / ELI_V2) * ELIC;
dxdt_SIVD = -SIV_KA * SIVD;
dxdt_SIVC =  SIV_KA * SIVD - (SIV_CL / SIV_V2) * SIVC;
dxdt_GABD = -GAB_KA * GABD;
dxdt_GABC =  GAB_KA * GABD - (GAB_CL / GAB_V2) * GABC;
dxdt_PRED = -PRE_KA * PRED;
dxdt_PREC =  PRE_KA * PRED - (PRE_CL / PRE_V2) * PREC;
dxdt_MORD = -MOR_KA * MORD;
dxdt_MORC =  MOR_KA * MORD - (MOR_CL / MOR_V2) * MORC;
dxdt_NALD = -NAL_KA * NALD;
dxdt_NALC =  NAL_KA * NALD - (NAL_CL / NAL_V2) * NALC;
dxdt_AMID = -AMI_KA * AMID;
dxdt_AMIC =  AMI_KA * AMID - (AMI_CL / AMI_V2) * AMIC;
dxdt_ICSL = -ICS_KA * ICSL;
dxdt_ICSC =  ICS_KA * ICSL - (ICS_CL / ICS_V2) * ICSC;
dxdt_PPID = -PPI_KA * PPID;
dxdt_PPIC =  PPI_KA * PPID - (PPI_CL / PPI_V2) * PPIC;
dxdt_DXMD = -DXM_KA * DXMD;
dxdt_DXMC =  DXM_KA * DXMD - (DXM_CL / DXM_V2) * DXMC;

// =====================================================================
// 14. CIRCADIAN OSCILLATOR
// =====================================================================
dxdt_CIRCS =  OMEGA * CIRCC;
dxdt_CIRCC = -OMEGA * CIRCS;

$TABLE
double AWAKEo = 1.0 / (1.0 + exp(-KWAKE * (CIRCC - AWKCUT)));

double Cgefo = GEFC / GEF_V2;
double Ccamo = CAMC / CAM_V2;
double Celio = ELIC / ELI_V2;
double Csivo = SIVC / SIV_V2;
double Cgabo = GABC / GAB_V2;
double Cmoro = MORC / MOR_V2;

double Ugefn = Cgefo * GEF_FU / GEF_MW * 1000.0;   // on-treatment
double Ucamn = Ccamo * CAM_FU / CAM_MW * 1000.0;
double Uelin = Celio * ELI_FU / ELI_MW * 1000.0;
double Usivn = Csivo * SIV_FU / SIV_MW * 1000.0;
double Ugefo = Ugefn * ADHR;                       // population-average
double Ucamo = Ucamn * ADHR;
double Uelio = Uelin * ADHR;
double Usivo = Usivn * ADHR;
double shifto = 1.0 + ATPX / KATPC;

double OCC3 = 1.0 - (1.0 - Ugefo / (Ugefo + GEF_IC3 * shifto)) *
                    (1.0 - Ucamo / (Ucamo + CAM_IC3 * shifto)) *
                    (1.0 - Uelio / (Uelio + ELI_IC3 * shifto)) *
                    (1.0 - Usivo / (Usivo + SIV_IC3 * shifto));
double OCC23 = 1.0 - (1.0 - Ugefn / (Ugefn + GEF_IC3 * GEF_SEL)) *
                     (1.0 - Ucamn / (Ucamn + CAM_IC3 * CAM_SEL)) *
                     (1.0 - Uelin / (Uelin + ELI_IC3 * ELI_SEL)) *
                     (1.0 - Usivn / (Usivn + SIV_IC3 * SIV_SEL));

double Emoro = MOR_EMX * Cmoro / (Cmoro + MOR_EC50);
double CENGo = AGEF * (1.0 + WWU * WIND) / (pow(fmax(0.05, GABI), WGAB)
                                          * pow(fmax(0.05, DESC), WDES));
CENGo = CENGo * (1.0 - MOR_NTS * Emoro / fmax(1e-9, MOR_EMX));

double Egabo = GAB_EMX * Cgabo / (Cgabo + GAB_EC50);
double EA2Do = FA2DPRE * Egabo;
double TUSSo = AFIR * CENGo * (1.0 + HABT) * (1.0 - EA2Do);
double THRSo = THR0 * (SEXF > 0.5 ? FTHRSEX : 1.0)
             * (1.0 + WCORT * (CORT - CORTB))
             * (1.0 + WDESCT * (DESC - 1.0));
THRSo = fmax(0.4, THRSo);

// LOOP GAIN: the product of the peripheral and central amplifications.
// > 1 means the vicious circle sustains itself without the trigger.
double PERG = PSEN * (1.0 + 0.55 * POS(P2X3E - 1.0) + 0.45 * POS(ATPX / ATP0 - 1.0));
double LOOP = LGAIN * PERG * CENGo / (PERG + CENGo - 1.0 + 1e-9)
            * (1.0 + 0.30 * POS(ECTO));

// challenge tests — capsaicin bypasses P2X3, inhaled ATP does not
double C5    = (SEXF > 0.5 ? C5F0 : C5M0) / (TRPV * pow(fmax(0.05, CENGo), WC5CEN));
double ATPC5 = ATPC50 / (P2X3E * pow(fmax(0.05, CENGo), WC5CEN) * fmax(0.02, 1.0 - OCC3 * EMAX3));

// taste adverse-event reporting rate (% of subjects)
double TASTEAE = PBOAE + (100.0 - PBOAE) * HILL(DYSG, AE50, HAE);
double GABAE   = (Cgabo > 0.01) ? 100.0 * GAB_AE * Cgabo / (Cgabo + GAB_EC50) : 0.0;

// cough severity strata used clinically
double SEVERE  = (CFQ > 25.0) ? 1.0 : 0.0;

$CAPTURE @annotated
AWAKEo : Wake gate (1 = awake, 0 = asleep)
CENGo  : Central gain (1 = healthy)
TUSSo  : Total tussive drive
THRSo  : Cough threshold
PERG   : Peripheral gain
LOOP   : Loop gain (>1 = self-sustaining)
OCC3   : P2X3 homotrimer occupancy (efficacy target)
OCC23  : P2X2/3 heterotrimer occupancy (taste target)
C5     : Capsaicin C5 cough threshold (uM)
ATPC5  : ATP challenge cough threshold (uM)
TASTEAE: Taste-related adverse event rate (%)
GABAE  : Gabapentinoid CNS adverse event rate (%)
SEVERE : Severe cough flag (>25 coughs/h)
Cgefo  : Gefapixant plasma concentration (mg/L)
Ccamo  : Camlipixant plasma concentration (mg/L)
Celio  : Eliapixant plasma concentration (mg/L)
Csivo  : Sivopixant plasma concentration (mg/L)
Cgabo  : Gabapentin plasma concentration (mg/L)
Cmoro  : Morphine plasma concentration (mg/L)
'

mod <- mcode_cache("rcc", code, soloc = tempdir())

CMT_NAMES <- as.list(mod)$cmt
CMT <- function(name) {
  i <- which(CMT_NAMES == name)
  if (length(i) != 1) stop("unknown compartment: ", name)
  i
}

## =====================================================================
##  PATIENT PHENOTYPES
##  Each is a parameter set only. No phenotype is given a different model
##  structure, and none is given a hard-coded cough frequency: the cough
##  each one ends up with is produced by running the natural history.
## =====================================================================
pat_healthy <- list(INSULT0 = 0.0, LGAIN = 1.00, SEXF = 1, ESTRO = 1.0)

## The classic referral-clinic patient: woman in her late 50s, cough for
## years, dating from a chest infection, normal CT and spirometry.
pat_rcc <- list(INSULT0 = 1.00, LGAIN = 1.00, SEXF = 1, ESTRO = 0.30,
                AGEF = 1.05)

## Same insult, LOWER loop gain: this subject recovers. Post-infectious
## cough and refractory chronic cough are the SAME model on either side
## of the fold.
pat_postviral <- list(INSULT0 = 1.00, LGAIN = 0.55, SEXF = 1, ESTRO = 1.0,
                      AGEF = 0.98)

## Male patient — higher intrinsic cough threshold, higher capsaicin C5.
pat_male <- list(INSULT0 = 1.00, LGAIN = 1.00, SEXF = 0, ESTRO = 1.0)

## Eosinophilic phenotype: non-asthmatic eosinophilic bronchitis /
## cough-variant asthma. The trait ICS can actually treat.
pat_eos <- list(INSULT0 = 0.55, LGAIN = 0.45, SEXF = 1, ESTRO = 0.60,
                EOS0 = 1.35, EOSDRV = 1.0)

## Acid-reflux phenotype: PPI has something to suppress.
pat_acid <- list(INSULT0 = 0.55, LGAIN = 0.72, SEXF = 1, ESTRO = 0.60,
                 ACIDSET = 1.90, NACSET = 1.05)

## Non-acid reflux phenotype: identical cough, identical impedance-detected
## reflux events, NOTHING for a PPI to suppress. This pair is the model's
## explanation of the Cochrane PPI result.
pat_nonacid <- list(INSULT0 = 0.55, LGAIN = 0.72, SEXF = 1, ESTRO = 0.60,
                    ACIDSET = 1.05, NACSET = 2.6)

## ACE-inhibitor cough: a purely PERIPHERAL cough with no central
## sensitisation, which is why it resolves completely on withdrawal.
pat_acei <- list(INSULT0 = 0.15, LGAIN = 0.60, SEXF = 1, ESTRO = 1.0, ACEI = 1)

## IPF-associated cough: mechanical traction on stretch receptors plus
## epithelial injury. The mechanical arm is P2X3-independent.
pat_ipf <- list(INSULT0 = 0.85, LGAIN = 0.98, SEXF = 0, ESTRO = 1.0,
                FTRACT = 1.55)

## Upper-airway cough syndrome: extrapulmonary afferent drive.
pat_uacs <- list(INSULT0 = 0.50, LGAIN = 0.86, SEXF = 1, ESTRO = 1.0,
                 FUACS = 1.35)

## Registration-trial population: enriched for high baseline cough
## frequency, as COUGH-1/2 and SOOTHE were.
pat_trial <- list(INSULT0 = 1.05, LGAIN = 1.02, SEXF = 1, ESTRO = 0.30,
                  AGEF = 1.05)

## SOOTHE high-baseline stratum (>= 32 coughs/h awake at screening).
pat_highbl <- list(INSULT0 = 1.15, LGAIN = 1.05, SEXF = 1, ESTRO = 0.30,
                   AGEF = 1.06)

## Early disease: 10 weeks of cough, central sensitisation not yet
## entrenched. Used to test the early-vs-late prediction.
pat_early <- list(INSULT0 = 1.00, LGAIN = 1.00, SEXF = 1, ESTRO = 0.30,
                  AGEF = 1.02)

## =====================================================================
##  NATURAL-HISTORY PRE-RUN
##
##  Nobody is initialised in the diseased state. Every virtual patient
##  starts healthy, receives their insult, and is simulated forward for
##  PRERUN days. Whatever cough they have at the end is the cough the
##  model gave them. Trials then start from that state.
##
##  PRERUN is a whole number of days, so the circadian phase at the end
##  is identical to the phase at t = 0 and trials still begin at 07:00.
## =====================================================================
PRERUN_DEFAULT <- 730     # 2 years of chronic cough before enrolment
RUNIN <- 28               # screening -> randomisation interval (days)

settle <- function(patient = pat_rcc, prerun = PRERUN_DEFAULT, extra = list()) {
  ## prerun = 0 means "start from the model's own healthy initial condition"
  if (prerun <= 0) return(NULL)
  m <- param(mod, modifyList(patient, extra))
  o <- mrgsim_df(m, end = prerun, delta = 1, add = seq(0, prerun, 0.25))
  last <- o[nrow(o), ]
  ic <- as.list(last[CMT_NAMES])
  ## reset the accumulators and the trial-artefact states; keep the disease
  ic$CACC <- 0; ic$CAWK <- 0; ic$TAWK <- 0
  ic$FLARE <- 0
  ic
}

## =====================================================================
##  DOSING EVENT BUILDERS
##  Oral doses go into the gut depot AFTER bioavailability, so F is
##  applied here rather than with F_CMT.
## =====================================================================
DAYS  <- 168
TGRID <- sort(unique(c(seq(0, 3, 1 / 48), seq(3, DAYS, 1 / 8))))

ev_oral <- function(cmtname, dose, F1, ii = 0.5, dur = DAYS, start = 0)
  ev(amt = dose * F1, cmt = CMT(cmtname), ii = ii,
     addl = ceiling(dur / ii) - 1, time = start)

ev_gefapixant  <- function(dose = 45,  dur = DAYS) ev_oral("GEFD", dose, 0.90, 0.5, dur)
ev_camlipixant <- function(dose = 50,  dur = DAYS) ev_oral("CAMD", dose, 0.75, 0.5, dur)
ev_eliapixant  <- function(dose = 75,  dur = DAYS) ev_oral("ELID", dose, 0.70, 0.5, dur)
ev_sivopixant  <- function(dose = 150, dur = DAYS) ev_oral("SIVD", dose, 0.60, 1.0, dur)
ev_gabapentin  <- function(dose = 600, dur = DAYS) {
  ## LAT1-mediated absorption SATURATES: F falls from ~60% at 300 mg to
  ## ~35% at 600 mg per dose. Modelled explicitly rather than assumed
  ## dose-proportional, because it is why gabapentin is given TID and why
  ## the dose-response flattens above 1800 mg/d.
  Fd <- 0.60 / (1 + dose / 700)
  ev_oral("GABD", dose, Fd, 1 / 3, dur)
}
ev_pregabalin  <- function(dose = 150, dur = DAYS) ev_oral("PRED", dose, 0.90, 0.5, dur)
ev_morphine    <- function(dose = 10,  dur = DAYS) ev_oral("MORD", dose, 0.30, 0.5, dur)
ev_nalbuphine  <- function(dose = 108, dur = DAYS) ev_oral("NALD", dose, 0.20, 0.5, dur)
ev_amitrip     <- function(dose = 25,  dur = DAYS) ev_oral("AMID", dose, 0.50, 1.0, dur)
ev_ics         <- function(dose = 400, dur = DAYS) ev_oral("ICSL", dose, 0.35, 0.5, dur)
ev_ppi         <- function(dose = 40,  dur = DAYS) ev_oral("PPID", dose, 0.64, 1.0, dur)
ev_dxm         <- function(dose = 30,  dur = DAYS) ev_oral("DXMD", dose, 0.15, 0.25, dur)

combine_ev <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (length(parts) == 0) return(NULL)
  Reduce(function(a, b) a + b, parts)
}

## =====================================================================
##  RUNNER
##  Every trial simulation: settle the patient, then randomise (TRIALON,
##  FLARE0) and dose. The placebo arm gets everything except the drug.
## =====================================================================
run <- function(patient = pat_rcc, params = list(), events = NULL,
                days = DAYS, ic = NULL, trial = TRUE, tgrid = NULL,
                prerun = PRERUN_DEFAULT) {
  if (is.null(ic)) ic <- settle(patient, prerun)
  pp <- modifyList(patient, params)
  if (trial) pp$TRIALON <- 1 else pp$TRIALON <- 0
  ## SETIC = 0 tells $MAIN to leave the supplied state vector alone. Without
  ## it, $MAIN reruns its healthy initial-condition block on every record and
  ## the natural-history pre-run is silently discarded.
  pp$SETIC <- if (length(ic)) 0 else 1

  ## THE ENROLMENT FLARE MUST HAPPEN BEFORE BASELINE, NOT AFTER IT.
  ## Patients present, and are screened, when their cough is at its worst;
  ## they are randomised a few weeks later, and the trial's own baseline
  ## measurement is taken at THAT point, already part-way down the flare.
  ## Injecting the flare at randomisation instead made the placebo arm get
  ## 64% WORSE over 12 weeks, which is the opposite of every cough trial
  ## ever run. So the flare is applied, then RUNIN days are simulated
  ## off-treatment, and the state at the end of that is trial day 0.
  ## RUNIN is a whole number of days, so the circadian phase is preserved.
  if (trial && length(ic)) {
    m0 <- param(mod, modifyList(pp, list(TRIALON = 0, SETIC = 0)))
    m0 <- init(m0, ic)
    m0 <- init(m0, list(FLARE = as.numeric(pp$FLARE0 %||% 0.10)))
    o0 <- mrgsim_df(m0, end = RUNIN, delta = 1, add = seq(0, RUNIN, 0.25))
    ic <- as.list(o0[nrow(o0), CMT_NAMES])
    ic$CACC <- 0; ic$CAWK <- 0; ic$TAWK <- 0
  }

  m <- param(mod, pp)
  if (length(ic)) m <- init(m, ic)
  g <- if (is.null(tgrid)) sort(unique(c(seq(0, min(3, days), 1 / 24),
                                         seq(0, days, 1 / 8)))) else tgrid
  if (is.null(events)) mrgsim_df(m, end = days, add = g, delta = 1)
  else                 mrgsim_df(m, events = events, end = days, add = g, delta = 1)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## Daily 24-h and awake cough frequencies, differenced from the
## accumulators. This is exactly how a cough monitor reports.
daily <- function(d) {
  wh <- which(abs(d$time - round(d$time)) < 1e-8)
  dd <- d[wh, ]
  dd <- dd[!duplicated(round(dd$time)), ]
  n <- nrow(dd)
  data.frame(
    day   = round(dd$time[-1]),
    cf24  = diff(dd$CACC) / 24,
    cfawk = diff(dd$CAWK) / pmax(1e-9, diff(dd$TAWK)),
    LCQ   = dd$LCQ[-1], VAS = dd$VASC[-1], CSD = dd$CSDS[-1],
    URG   = dd$URG[-1], DYSG = dd$DYSG[-1], ADHR = dd$ADHR[-1],
    TASTEAE = dd$TASTEAE[-1], LOOP = dd$LOOP[-1],
    CENG  = dd$CENGo[-1], WIND = dd$WIND[-1], AFIR = dd$AFIR[-1],
    P2X3E = dd$P2X3E[-1], ECTO = dd$ECTO[-1], CORT = dd$CORT[-1],
    OCC3  = dd$OCC3[-1], OCC23 = dd$OCC23[-1],
    C5    = dd$C5[-1], ATPC5 = dd$ATPC5[-1], COMP = dd$COMP[-1])
}

at_day <- function(dd, day) dd[which.min(abs(dd$day - day)), ]

## =====================================================================
##  THE 24 SCENARIOS
## =====================================================================
scenarios <- list(
  list(id = "S00", label = "Healthy control (no insult)",
       patient = pat_healthy, trial = FALSE, params = list(), ev = NULL),
  list(id = "S01", label = "Natural history - RCC, untreated, not in a trial",
       patient = pat_rcc, trial = FALSE, params = list(), ev = NULL),
  list(id = "S02", label = "Post-infectious cough (loop gain below the fold)",
       patient = pat_postviral, trial = FALSE, params = list(), ev = NULL),
  list(id = "S03", label = "PLACEBO ARM - RCC in a trial, no drug",
       patient = pat_trial, trial = TRUE, params = list(), ev = NULL),
  list(id = "S04", label = "Gefapixant 45 mg BID (COUGH-1/2)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_gefapixant(45)),
  list(id = "S05", label = "Gefapixant 15 mg BID",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_gefapixant(15)),
  list(id = "S06", label = "Camlipixant 50 mg BID (SOOTHE)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_camlipixant(50)),
  list(id = "S07", label = "Camlipixant 50 mg BID, high-baseline stratum",
       patient = pat_highbl, trial = TRUE, params = list(), ev = ev_camlipixant(50)),
  list(id = "S08", label = "Eliapixant 75 mg BID (PAGANINI)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_eliapixant(75)),
  list(id = "S09", label = "Eliapixant 150 mg BID",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_eliapixant(150)),
  list(id = "S10", label = "Sivopixant 150 mg QD",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_sivopixant(150)),
  list(id = "S11", label = "Gabapentin 600 mg TID (Ryan 2012)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_gabapentin(600)),
  list(id = "S12", label = "Pregabalin 150 mg BID + speech therapy (Vertigan 2016)",
       patient = pat_trial, trial = TRUE, params = list(ESLT = 0.42),
       ev = ev_pregabalin(150)),
  list(id = "S13", label = "Morphine SR 10 mg BID (Morice 2007)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_morphine(10)),
  list(id = "S14", label = "Amitriptyline 25 mg nocte",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_amitrip(25)),
  list(id = "S15", label = "Dextromethorphan 30 mg QID (weak comparator)",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_dxm(30)),
  list(id = "S16", label = "Speech and language therapy alone (PSALTI)",
       patient = pat_trial, trial = TRUE, params = list(ESLT = 0.42, ECBT = 0.30),
       ev = NULL),
  list(id = "S17", label = "ICS 400 ug BID in TRUE RCC (non-eosinophilic) - expected null",
       patient = pat_trial, trial = TRUE, params = list(), ev = ev_ics(400)),
  list(id = "S18", label = "ICS 400 ug BID in EOSINOPHILIC cough (NAEB/CVA)",
       patient = pat_eos, trial = TRUE, params = list(), ev = ev_ics(400)),
  list(id = "S19", label = "PPI 40 mg OD in NON-ACID reflux cough - expected null",
       patient = pat_nonacid, trial = TRUE, params = list(), ev = ev_ppi(40)),
  list(id = "S20", label = "PPI 40 mg OD in ACID reflux cough",
       patient = pat_acid, trial = TRUE, params = list(), ev = ev_ppi(40)),
  list(id = "S21", label = "ACE inhibitor withdrawal (ACEi cough)",
       patient = pat_acei, trial = FALSE, params = list(ACEI = 0), ev = NULL),
  list(id = "S22", label = "Nalbuphine ER 108 mg BID in IPF cough (CANAL)",
       patient = pat_ipf, trial = TRUE, params = list(), ev = ev_nalbuphine(108)),
  list(id = "S23", label = "COMBINATION - camlipixant 50 mg BID + speech therapy",
       patient = pat_trial, trial = TRUE, params = list(ESLT = 0.42, ECBT = 0.30),
       ev = ev_camlipixant(50)),
  list(id = "S24", label = "EARLY vs LATE - camlipixant started at 10 weeks of cough",
       patient = pat_early, trial = TRUE, params = list(), ev = ev_camlipixant(50),
       prerun = 70)
)

run_all <- function(days = DAYS, verbose = TRUE) {
  ics <- list()
  out <- lapply(scenarios, function(s) {
    pr <- if (!is.null(s$prerun)) s$prerun else PRERUN_DEFAULT
    key <- paste0(digest_key(s$patient), "_", pr)
    if (is.null(ics[[key]])) ics[[key]] <<- settle(s$patient, pr)
    if (verbose) message("  ", s$id, "  ", s$label)
    d <- run(s$patient, s$params, s$ev, days = days,
             ic = ics[[key]], trial = s$trial, prerun = pr)
    d$id <- s$id; d$label <- s$label
    d
  })
  names(out) <- vapply(scenarios, function(s) s$id, "")
  out
}
digest_key <- function(p) paste(names(p), unlist(p), sep = "=", collapse = ";")

## =====================================================================
##  DIAGNOSTICS
## =====================================================================

## 1. BASELINE STABILITY.
##    If a healthy subject does not hold perfectly still, every treatment
##    effect below is contaminated by numerical drift rather than being
##    a pharmacological result.
## Drift is measured from day 30, AFTER the initial circadian/mucus
## transient. A healthy subject is a LIMIT CYCLE, not a fixed point — the
## overnight mucus swing that produces the morning peak also means the
## first few days are a transient, and counting that transient as "drift"
## would be measuring the wrong thing.
check_baseline <- function(days = 365, from = 30) {
  o <- run(pat_healthy, days = days, trial = FALSE, prerun = 0,
           tgrid = seq(0, days, 0.25))
  dd <- daily(o); dd <- dd[dd$day >= from, ]
  f <- function(v) 100 * (tail(v, 1) - head(v, 1)) / head(v, 1)
  g <- function(v) tail(v, 1) - head(v, 1)
  data.frame(
    state = c("cf24", "cfawk", "LCQ", "VAS", "URG", "AFIR",
              "CENG", "WIND", "P2X3E", "LOOP", "C5"),
    value_end = round(c(tail(dd$cf24, 1), tail(dd$cfawk, 1), tail(dd$LCQ, 1),
                        tail(dd$VAS, 1), tail(dd$URG, 1), tail(dd$AFIR, 1),
                        tail(dd$CENG, 1), tail(dd$WIND, 1), tail(dd$P2X3E, 1),
                        tail(dd$LOOP, 1), tail(dd$C5, 1)), 4),
    drift_abs = signif(c(g(dd$cf24), g(dd$cfawk), g(dd$LCQ), g(dd$VAS),
                         g(dd$URG), g(dd$AFIR), g(dd$CENG), g(dd$WIND),
                         g(dd$P2X3E), g(dd$LOOP), g(dd$C5)), 3),
    drift_pct = signif(c(f(dd$cf24), f(dd$cfawk), f(dd$LCQ), f(dd$VAS),
                         f(dd$URG), f(dd$AFIR), f(dd$CENG), g(dd$WIND),
                         f(dd$P2X3E), f(dd$LOOP), f(dd$C5)), 3))
}

## 2. THE PHENOTYPE BASELINE GRID — what natural history produced.
##    None of these numbers was assigned; each is what the model settled to.
check_phenotypes <- function(prerun = PRERUN_DEFAULT) {
  ps <- list(healthy = pat_healthy, RCC = pat_rcc, postviral = pat_postviral,
             male = pat_male, eosinophilic = pat_eos, acid_reflux = pat_acid,
             nonacid_reflux = pat_nonacid, ACEi = pat_acei, IPF = pat_ipf,
             UACS = pat_uacs, trial_pop = pat_trial, high_baseline = pat_highbl)
  do.call(rbind, lapply(names(ps), function(nm) {
    o  <- run(ps[[nm]], days = 14, trial = FALSE, prerun = prerun)
    dd <- daily(o); b <- at_day(dd, 14)
    data.frame(phenotype = nm,
               cf24 = round(b$cf24, 1), cfawk = round(b$cfawk, 1),
               LCQ = round(b$LCQ, 1), VAS = round(b$VAS), CSD = round(b$CSD, 1),
               urge = round(b$URG, 1), loop_gain = round(b$LOOP, 2),
               central_gain = round(b$CENG, 2), windup = round(b$WIND, 2),
               P2X3 = round(b$P2X3E, 2), ectopic = round(b$ECTO, 2),
               capsaicin_C5 = round(b$C5, 1))
  }))
}

## 3. SUSCEPTIBILITY, AND WHETHER THE DISEASE IS ACTUALLY BISTABLE.
##
##    RESULT, STATED UP FRONT BECAUSE IT IS NEGATIVE: this model is NOT
##    bistable. Sweeping the loop gain up from health and back down from
##    established disease gives two curves that lie on top of each other
##    (hysteresis < 0.3 coughs/h everywhere). What the model actually has
##    is a steep, MONOTONIC susceptibility relation, not two attractors.
##
##    Whether real refractory chronic cough is bistable is an open
##    question, and an earlier version of this file asserted that it was.
##    The diagnostic below is what refuted that assertion, and it is kept
##    in the file for exactly that reason. Post-viral resolution versus
##    refractory persistence is reproduced here, but by the COMBINATION
##    of loop gain and peripheral susceptibility (oestrogen status, age)
##    moving a continuous dose-response, not by crossing a fold.
##
##    A still earlier version claimed the bifurcation was in INSULT
##    AMPLITUDE. Sweeping that showed it was false too: every insult
##    from 0.25 to 1.4 lands on the same attractor, because the airway
##    insult is gone within weeks and plays no part in what sustains the
##    cough afterwards. The bifurcation parameter is the HOST LOOP GAIN.
##
##    That is the more defensible claim anyway. Refractory chronic cough
##    patients do not report unusually severe index infections; they
##    differ in susceptibility — female sex, age, family history, TRPV1
##    polymorphism. So the sweep is over LGAIN, and it is run in BOTH
##    directions to distinguish genuine bistability (two attractors, with
##    hysteresis) from a merely steep monotonic response.
check_bistability <- function(gains = seq(0.55, 1.15, by = 0.05),
                              horizon = 730) {
  ## UP-SWEEP: start healthy, apply the insult, see where you end up.
  up <- vapply(gains, function(g) {
    p <- modifyList(pat_rcc, list(LGAIN = g))
    at_day(daily(run(p, days = 14, trial = FALSE, prerun = horizon)), 14)$cf24
  }, 0)

  ## DOWN-SWEEP: start from an ESTABLISHED patient (settled at LGAIN = 1.15)
  ## and then lower the loop gain, e.g. as a treatment would. If the two
  ## branches differ, the disease has memory that the drug has to overcome.
  ic_est <- settle(modifyList(pat_rcc, list(LGAIN = 1.15)), horizon)
  down <- vapply(gains, function(g) {
    p <- modifyList(pat_rcc, list(LGAIN = g))
    at_day(daily(run(p, list(), NULL, days = 365, ic = ic_est,
                     trial = FALSE)), 365)$cf24
  }, 0)

  data.frame(loop_gain = gains,
             cf24_from_health = round(up, 1),
             cf24_from_disease = round(down, 1),
             hysteresis = round(down - up, 1),
             state = ifelse(up > 8, "refractory",
                     ifelse(up > 4, "mild persistent", "resolved")),
             bistable = ifelse(abs(down - up) > 1.0, "YES", "no"))
}

## 4. TRIAL ENDPOINT GRID at an arbitrary week.
endpoints <- function(res, week = 12) {
  day <- week * 7
  do.call(rbind, lapply(res, function(d) {
    dd <- daily(d); a <- at_day(dd, 1); b <- at_day(dd, day)
    data.frame(id = d$id[1],
               label      = substr(d$label[1], 1, 46),
               cf24_base  = round(a$cf24, 1),
               cf24_wk    = round(b$cf24, 1),
               pct_change = round(100 * (b$cf24 / a$cf24 - 1), 1),
               cfawk_wk   = round(b$cfawk, 1),
               dLCQ       = round(b$LCQ - a$LCQ, 2),
               dVAS       = round(b$VAS - a$VAS, 1),
               dCSD       = round(b$CSD - a$CSD, 2),
               taste_AE   = round(b$TASTEAE, 1),
               adherence  = round(100 * b$ADHR, 1),
               OCC3       = round(b$OCC3, 3),
               OCC23      = round(b$OCC23, 4))
  }))
}

## 5. PLACEBO-ADJUSTED EFFECTS — the only comparison a trial ever makes.
##    Reported as a geometric-mean ratio, exactly as COUGH-1/2 did.
placebo_adjusted <- function(res, week = 12, placebo_id = "S03") {
  g  <- endpoints(res, week)
  pl <- g[g$id == placebo_id, ]
  ids <- setdiff(g$id, c("S00", "S01", "S02", placebo_id, "S21"))
  do.call(rbind, lapply(ids, function(i) {
    r <- g[g$id == i, ]
    rr <- (r$cf24_wk / r$cf24_base) / (pl$cf24_wk / pl$cf24_base)
    data.frame(id = i, label = r$label,
               ratio_vs_placebo = round(rr, 3),
               pct_vs_placebo   = round(100 * (rr - 1), 1),
               dLCQ_vs_placebo  = round(r$dLCQ - pl$dLCQ, 2),
               dVAS_vs_placebo  = round(r$dVAS - pl$dVAS, 1),
               taste_AE_pct     = r$taste_AE,
               adherence_pct    = r$adherence)
  }))
}

## 6. THE THERAPEUTIC WINDOW.
##    One drug, one efficacy curve, one dysgeusia curve, and four decades
##    of selectivity between them. This is the model's central claim, so
##    it is swept rather than asserted.
check_window <- function(sels = c(1, 3, 6, 10, 30, 70, 150, 400, 1500, 5000),
                         week = 12) {
  ic <- settle(pat_trial)
  pl <- daily(run(pat_trial, list(), NULL, days = week * 7, ic = ic))
  plb <- at_day(pl, 1)$cf24; plw <- at_day(pl, week * 7)$cf24
  do.call(rbind, lapply(sels, function(s) {
    d  <- run(pat_trial, list(GEF_SEL = s), ev_gefapixant(45),
              days = week * 7, ic = ic)
    dd <- daily(d); a <- at_day(dd, 1); b <- at_day(dd, week * 7)
    rr <- (b$cf24 / a$cf24) / (plw / plb)
    data.frame(selectivity = s,
               pct_vs_placebo = round(100 * (rr - 1), 1),
               taste_AE_pct   = round(b$TASTEAE, 1),
               adherence_pct  = round(100 * b$ADHR, 1),
               OCC3 = round(b$OCC3, 3), OCC23 = round(b$OCC23, 4))
  }))
}

## 7. EFFICACY vs EFFECTIVENESS.
##    Turn the discontinuation hazard off and the SAME drug at the SAME
##    dose gets bigger. The difference is not pharmacology, it is
##    tolerability, and an ITT estimate measures the second one.
check_eff_vs_effectiveness <- function(week = 12) {
  ic <- settle(pat_trial)
  pl <- daily(run(pat_trial, list(), NULL, days = week * 7, ic = ic))
  base <- function(pars) {
    dd <- daily(run(pat_trial, pars, ev_gefapixant(45), days = week * 7, ic = ic))
    a <- at_day(dd, 1); b <- at_day(dd, week * 7)
    c(rr = (b$cf24 / a$cf24) /
          (at_day(pl, week * 7)$cf24 / at_day(pl, 1)$cf24),
      adh = b$ADHR)
  }
  itt <- base(list())
  pp  <- base(list(KDROP = 0))
  data.frame(analysis = c("intention-to-treat (adherence eroded by dysgeusia)",
                          "per-protocol (no discontinuation)"),
             pct_vs_placebo = round(100 * (c(itt["rr"], pp["rr"]) - 1), 1),
             adherence_pct  = round(100 * c(itt["adh"], pp["adh"]), 1),
             row.names = NULL)
}

## 8. THE CHALLENGE-TEST DISSOCIATION.
##    Capsaicin gates TRPV1 and bypasses P2X3; inhaled ATP does not. The
##    model is never told this dissociation — it falls out of where each
##    agonist enters the afferent equation.
check_challenges <- function(week = 4) {
  ic <- settle(pat_trial)
  arms <- list(placebo = NULL, gefapixant45 = ev_gefapixant(45),
               camlipixant50 = ev_camlipixant(50), morphine10 = ev_morphine(10),
               gabapentin600 = ev_gabapentin(600))
  do.call(rbind, lapply(names(arms), function(nm) {
    dd <- daily(run(pat_trial, list(), arms[[nm]], days = week * 7, ic = ic))
    a <- at_day(dd, 1); b <- at_day(dd, week * 7)
    data.frame(arm = nm,
               capsaicin_C5_fold = round(b$C5 / a$C5, 2),
               ATP_threshold_fold = round(b$ATPC5 / a$ATPC5, 2),
               cf24_fold = round(b$cf24 / a$cf24, 2))
  }))
}

## 9. HEAD-TO-HEAD AGAINST PUBLISHED TRIAL RESULTS.
##    Each row states the observation, its source and the model value.
##    Everything is placebo-adjusted in the matching population and at
##    the matching timepoint.
validate_trials <- function() {
  ic_tr <- settle(pat_trial)
  ic_hb <- settle(pat_highbl)
  ic_eo <- settle(pat_eos)
  ic_na <- settle(pat_nonacid)
  ic_ac <- settle(pat_acid)
  ic_ip <- settle(pat_ipf)

  adj <- function(icx, patx, pars, evx, week, plpars = list()) {
    d  <- daily(run(patx, pars, evx, days = week * 7, ic = icx))
    p  <- daily(run(patx, plpars, NULL, days = week * 7, ic = icx))
    a <- at_day(d, 1); b <- at_day(d, week * 7)
    pa <- at_day(p, 1); pb <- at_day(p, week * 7)
    rr <- (b$cf24 / a$cf24) / (pb$cf24 / pa$cf24)
    c(pct = 100 * (rr - 1), dLCQ = b$LCQ - a$LCQ - (pb$LCQ - pa$LCQ),
      dVAS = b$VAS - a$VAS - (pb$VAS - pa$VAS), taste = b$TASTEAE)
  }

  rows <- list()
  add <- function(trial, endpoint, observed, model)
    rows[[length(rows) + 1]] <<- data.frame(trial = trial, endpoint = endpoint,
                                            observed = observed, model = model)

  g45_12 <- adj(ic_tr, pat_trial, list(), ev_gefapixant(45), 12)
  add("COUGH-1 (12 wk)", "gefapixant 45 mg BID, 24-h cough vs placebo (%)",
      "-18.5", sprintf("%.1f", g45_12["pct"]))
  add("COUGH-1 (12 wk)", "gefapixant 45 mg BID, taste-related AE (%)",
      "58", sprintf("%.0f", g45_12["taste"]))

  g15_12 <- adj(ic_tr, pat_trial, list(), ev_gefapixant(15), 12)
  add("COUGH-1 (12 wk)", "gefapixant 15 mg BID, 24-h cough vs placebo (%)",
      "-1.4 (ns)", sprintf("%.1f", g15_12["pct"]))
  add("COUGH-1 (12 wk)", "gefapixant 15 mg BID, taste-related AE (%)",
      "~32", sprintf("%.0f", g15_12["taste"]))

  g45_24 <- adj(ic_tr, pat_trial, list(), ev_gefapixant(45), 24)
  add("COUGH-2 (24 wk)", "gefapixant 45 mg BID, 24-h cough vs placebo (%)",
      "-14.6", sprintf("%.1f", g45_24["pct"]))
  add("COUGH-2 (24 wk)", "gefapixant 45 mg BID, taste-related AE (%)",
      "69", sprintf("%.0f", g45_24["taste"]))

  c50_4 <- adj(ic_hb, pat_highbl, list(), ev_camlipixant(50), 4)
  add("SOOTHE (4 wk)", "camlipixant 50 mg BID, high-baseline stratum (%)",
      "~-34", sprintf("%.1f", c50_4["pct"]))
  add("SOOTHE (4 wk)", "camlipixant 50 mg BID, taste-related AE (%)",
      "6.5", sprintf("%.0f", c50_4["taste"]))

  e75 <- adj(ic_tr, pat_trial, list(), ev_eliapixant(75), 12)
  add("PAGANINI (12 wk)", "eliapixant 75 mg BID, 24-h cough vs placebo (%)",
      "-17.6", sprintf("%.1f", e75["pct"]))
  e150 <- adj(ic_tr, pat_trial, list(), ev_eliapixant(150), 12)
  add("PAGANINI (12 wk)", "eliapixant 150 mg BID, taste-related AE (%)",
      "~21", sprintf("%.0f", e150["taste"]))

  s150 <- adj(ic_tr, pat_trial, list(), ev_sivopixant(150), 12)
  add("Sivopixant ph2b (12 wk)", "150 mg QD, 24-h cough vs placebo (%)",
      "-18.6 (ns)", sprintf("%.1f", s150["pct"]))

  gab <- adj(ic_tr, pat_trial, list(), ev_gabapentin(600), 10)
  add("Ryan 2012 (10 wk)", "gabapentin 1800 mg/d, LCQ vs placebo",
      "+1.80", sprintf("%+.2f", gab["dLCQ"]))
  add("Ryan 2012 (10 wk)", "gabapentin 1800 mg/d, cough VAS vs placebo (mm)",
      "-12.1", sprintf("%.1f", gab["dVAS"]))
  add("Ryan 2012 (10 wk)", "gabapentin 1800 mg/d, cough frequency vs placebo (%)",
      "-27", sprintf("%.1f", gab["pct"]))

  mor <- adj(ic_tr, pat_trial, list(), ev_morphine(10), 4)
  add("Morice 2007 (4 wk)", "morphine SR 10 mg BID, LCQ vs placebo",
      "+3.2", sprintf("%+.2f", mor["dLCQ"]))

  slt <- adj(ic_tr, pat_trial, list(ESLT = 0.42, ECBT = 0.30), NULL, 4)
  add("PSALTI (4 wk)", "speech therapy, LCQ vs control",
      "+1.53", sprintf("%+.2f", slt["dLCQ"]))
  add("PSALTI (4 wk)", "speech therapy, cough frequency vs control (%)",
      "-41", sprintf("%.1f", slt["pct"]))

  pre <- adj(ic_tr, pat_trial, list(ESLT = 0.42), ev_pregabalin(150), 4,
             plpars = list(ESLT = 0.42))
  add("Vertigan 2016 (4 wk)", "pregabalin added to speech therapy, LCQ",
      "+2.5 over SLT alone", sprintf("%+.2f", pre["dLCQ"]))

  ics_r <- adj(ic_tr, pat_trial, list(), ev_ics(400), 8)
  add("ICS in unselected RCC", "24-h cough vs placebo (%)",
      "null", sprintf("%.1f", ics_r["pct"]))
  ics_e <- adj(ic_eo, pat_eos, list(), ev_ics(400), 8)
  add("ICS in eosinophilic cough", "24-h cough vs placebo (%)",
      "large benefit", sprintf("%.1f", ics_e["pct"]))

  ppi_n <- adj(ic_na, pat_nonacid, list(), ev_ppi(40), 8)
  add("PPI, non-acid reflux (Cochrane)", "24-h cough vs placebo (%)",
      "null", sprintf("%.1f", ppi_n["pct"]))
  ppi_a <- adj(ic_ac, pat_acid, list(), ev_ppi(40), 8)
  add("PPI, documented acid reflux", "24-h cough vs placebo (%)",
      "modest benefit", sprintf("%.1f", ppi_a["pct"]))

  nal <- adj(ic_ip, pat_ipf, list(), ev_nalbuphine(108), 3)
  add("CANAL (IPF cough)", "nalbuphine ER, daytime cough vs placebo (%)",
      "-53 (=-75.7 vs -22.2)", sprintf("%.1f", nal["pct"]))

  ## ACEi withdrawal is not a placebo-controlled comparison; report the
  ## absolute resolution instead.
  d_acei <- daily(run(pat_acei, list(ACEI = 0), NULL, days = 56, trial = FALSE))
  add("ACEi withdrawal", "cough frequency at 4 weeks (% of baseline)",
      "largely resolved", sprintf("%.0f%%",
        100 * at_day(d_acei, 28)$cf24 / at_day(d_acei, 1)$cf24))

  do.call(rbind, rows)
}

## 10. EARLY vs LATE — a predicted effect that the model DOES NOT show.
##     The expectation was that treating before central sensitisation is
##     entrenched would do more. It does not: -21.2% at 10 weeks of
##     disease versus -21.4% at five years. The reason is visible in the
##     columns printed alongside — wind-up and microglial activation have
##     already reached steady state by 10 weeks (0.31 vs 0.29), so there
##     is no "early" window left to exploit. If clinical early treatment
##     really is better, the mechanism is NOT in this model, and a slower
##     plasticity state (structural reorganisation, learned behaviour)
##     would be the place to add it.
check_early_late <- function(week = 12) {
  prs <- c(early_10wk = 70, established_1y = 365, longstanding_5y = 1825)
  do.call(rbind, lapply(names(prs), function(nm) {
    icx <- settle(pat_rcc, prs[[nm]])
    d  <- daily(run(pat_rcc, list(), ev_camlipixant(50), days = week * 7, ic = icx))
    p  <- daily(run(pat_rcc, list(), NULL, days = week * 7, ic = icx))
    a <- at_day(d, 1); b <- at_day(d, week * 7)
    pa <- at_day(p, 1); pb <- at_day(p, week * 7)
    rr <- (b$cf24 / a$cf24) / (pb$cf24 / pa$cf24)
    data.frame(disease_duration = nm,
               cf24_at_entry = round(a$cf24, 1),
               windup_at_entry = round(d$WIND[1], 2),
               microglia_at_entry = round(at_day(d, 1)$CENG, 2),
               pct_vs_placebo = round(100 * (rr - 1), 1))
  }))
}

## 11. SENSITIVITY OF THE HEADLINE RESULT to the three parameters that
##     carry it: the P2X3 fraction of tussive drive, the taste-receptor
##     reserve, and the discontinuation hazard.
check_sensitivity <- function(week = 12) {
  ic <- settle(pat_trial)
  pl <- daily(run(pat_trial, list(), NULL, days = week * 7, ic = ic))
  ref <- at_day(pl, week * 7)$cf24 / at_day(pl, 1)$cf24
  one <- function(pars) {
    dd <- daily(run(pat_trial, pars, ev_gefapixant(45), days = week * 7, ic = ic))
    a <- at_day(dd, 1); b <- at_day(dd, week * 7)
    c(pct = 100 * ((b$cf24 / a$cf24) / ref - 1), taste = b$TASTEAE)
  }
  grid <- rbind(
    data.frame(knob = "WP2X",  value = c(0.20, 0.27, 0.34, 0.45, 0.60)),
    data.frame(knob = "TDG50", value = c(0.15, 0.30, 0.45, 0.70, 1.00)),
    data.frame(knob = "KDROP", value = c(0.000, 0.011, 0.022, 0.044, 0.088)))
  res <- t(mapply(function(k, v) one(setNames(list(v), k)), grid$knob, grid$value))
  cbind(grid, pct_vs_placebo = round(res[, "pct"], 1),
        taste_AE_pct = round(res[, "taste"], 1))
}

## 12. THE OBJECTIVE / SUBJECTIVE DISSOCIATION — direction corrected.
##     PROs track the URGE, counts track the MOTOR OUTPUT, and the urge
##     threshold is the lower of the two, so the two readouts CAN move by
##     different amounts. They do: every intervention below moves counts
##     further than it moves urge. But the ORDERING is the opposite of
##     what was expected. Cortical speech therapy was predicted to give
##     the largest count-to-urge ratio and gives the SMALLEST (1.26);
##     peripheral camlipixant gives the largest (1.79). Reading the
##     equations back, that is right and the prediction was careless:
##     cutting afferent drive lowers the input to BOTH thresholds, but
##     because the cough threshold is the steeper of the two the counts
##     fall further. Raising the cough threshold, as suppression training
##     does, cannot move the urge at all except through the slow
##     hypervigilance loop. The table is printed so the reader can check
##     this rather than take it on trust.
check_pro_dissociation <- function(week = 8) {
  ic <- settle(pat_trial)
  arms <- list(
    "speech therapy (cortical)"    = list(pars = list(ESLT = 0.42, ECBT = 0.30), ev = NULL),
    "camlipixant (peripheral)"     = list(pars = list(), ev = ev_camlipixant(50)),
    "morphine (central+descending)"= list(pars = list(), ev = ev_morphine(10)),
    "gabapentin (central)"         = list(pars = list(), ev = ev_gabapentin(600)))
  pl <- daily(run(pat_trial, list(), NULL, days = week * 7, ic = ic))
  pd <- at_day(pl, week * 7); pa <- at_day(pl, 1)
  do.call(rbind, lapply(names(arms), function(nm) {
    x <- arms[[nm]]
    dd <- daily(run(pat_trial, x$pars, x$ev, days = week * 7, ic = ic))
    a <- at_day(dd, 1); b <- at_day(dd, week * 7)
    cnt <- 100 * ((b$cf24 / a$cf24) / (pd$cf24 / pa$cf24) - 1)
    urg <- 100 * ((b$URG / a$URG) / (pd$URG / pa$URG) - 1)
    data.frame(intervention = nm,
               cough_count_vs_placebo_pct = round(cnt, 1),
               urge_vs_placebo_pct = round(urg, 1),
               dLCQ_vs_placebo = round((b$LCQ - a$LCQ) - (pd$LCQ - pa$LCQ), 2),
               PRO_per_count = round(((b$LCQ - a$LCQ) - (pd$LCQ - pa$LCQ)) /
                                     pmax(0.5, abs(cnt)) * 10, 3))
  }))
}

## 13. CIRCADIAN STRUCTURE — the awake/24-h ratio and the morning peak,
##     both of which emerge from the wake gate and overnight mucus
##     accumulation rather than being imposed.
check_circadian <- function() {
  ic <- settle(pat_trial)
  o  <- run(pat_trial, list(), NULL, days = 3, ic = ic,
            tgrid = seq(0, 3, 1 / 96))
  day2 <- o[o$time >= 1 & o$time < 2, ]
  hour <- floor((day2$time - 1) * 24)
  prof <- tapply(day2$CFQ, hour, mean)
  dd <- daily(o)
  list(hourly_profile = round(prof, 1),
       clock_time = sprintf("%02d:00", (as.integer(names(prof)) + 7) %% 24),
       awake_to_24h_ratio = round(at_day(dd, 2)$cfawk / at_day(dd, 2)$cf24, 2),
       peak_hour_after_waking = as.integer(names(which.max(prof))),
       sleep_trough = round(min(prof), 2))
}

## =====================================================================
##  REPORT
## =====================================================================
report <- function(days = DAYS) {
  cat("\n=====================================================================\n")
  cat(" REFRACTORY CHRONIC COUGH QSP MODEL\n")
  cat(" compartments:", length(CMT_NAMES),
      "  parameters:", length(as.list(mod)$param), "\n")
  cat("=====================================================================\n")

  cat("\n[1] BASELINE STABILITY - healthy subject, 365 d, INSULT0 = 0\n")
  print(check_baseline(), row.names = FALSE)

  cat("\n[2] PHENOTYPE BASELINES - produced by natural history, not assigned\n")
  print(check_phenotypes(), row.names = FALSE)

  cat("\n[3] SUSCEPTIBILITY SWEEP - up from health vs down from disease\n")
  print(check_bistability(), row.names = FALSE)

  cat("\n[4] RUNNING THE 25 SCENARIOS\n")
  res <- run_all(days = days)

  cat("\n[5] WEEK-12 ENDPOINT GRID\n")
  print(endpoints(res, 12), row.names = FALSE)

  cat("\n[6] PLACEBO-ADJUSTED EFFECTS AT WEEK 12\n")
  print(placebo_adjusted(res, 12), row.names = FALSE)

  cat("\n[7] THE THERAPEUTIC WINDOW - one drug, selectivity swept\n")
  print(check_window(), row.names = FALSE)

  cat("\n[8] EFFICACY vs EFFECTIVENESS - gefapixant 45 mg BID\n")
  print(check_eff_vs_effectiveness(), row.names = FALSE)

  cat("\n[9] CHALLENGE-TEST DISSOCIATION (fold change from baseline)\n")
  print(check_challenges(), row.names = FALSE)

  cat("\n[10] HEAD-TO-HEAD AGAINST PUBLISHED TRIALS\n")
  print(validate_trials(), row.names = FALSE)

  cat("\n[11] EARLY vs LATE TREATMENT (model prediction, untested)\n")
  print(check_early_late(), row.names = FALSE)

  cat("\n[12] SENSITIVITY OF THE HEADLINE RESULT\n")
  print(check_sensitivity(), row.names = FALSE)

  cat("\n[13] OBJECTIVE / SUBJECTIVE DISSOCIATION\n")
  print(check_pro_dissociation(), row.names = FALSE)

  cat("\n[14] CIRCADIAN STRUCTURE\n")
  print(check_circadian())

  cat("\n=====================================================================\n")
  invisible(res)
}

if (identical(environment(), globalenv()) &&
    !interactive() && sys.nframe() == 0) {
  report()
}
