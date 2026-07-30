# =====================================================================
# Retinopathy of Prematurity (ROP) — mrgsolve QSP model
#   Author : Claude Code Routine (2026-07-30)
#
#   Scope  : the two-phase natural history of ROP in the extremely
#            preterm retina — phase 1 vaso-obliteration / arrested
#            vascularisation, the metabolic switch, phase 2 hypoxia-
#            driven neovascularisation — together with every lever that
#            is actually pulled in a NICU: the SpO2 target range,
#            supplemental oxygen at prethreshold, packed-red-cell
#            transfusion, nutrition/IGF-1, AA:DHA supplementation,
#            diode laser ablation, and intravitreal bevacizumab /
#            ranibizumab / aflibercept with their systemic escape.
#
# ---------------------------------------------------------------------
#   WHY THIS MODEL IS BUILT THE WAY IT IS
#
#   ROP is routinely described as a disease of "too much then too little
#   oxygen".  That phrasing hides the only two quantities that matter,
#   and it hides the fact that they are set by different people:
#
#     SUPPLY = PCHOR, the choroidal oxygen tension.  Set by the
#              ventilator.  Not autoregulated — the choroid tracks PaO2
#              almost linearly.
#     DEMAND = PCRIT = MO2 · h² / 2k, the choroidal oxygen tension
#              REQUIRED to oxygenate an avascular retina of thickness h
#              consuming MO2.  Set by maturation.  No drug moves it.
#
#   Write those two separately and the whole disease falls out of their
#   collision.  PINNER = max(0, PCHOR − PCRIT) is the oxygen tension in
#   the inner avascular retina; the date at which PCHOR = PCRIT is the
#   phase 1 → phase 2 switch; and the sign of every oxygen intervention
#   is the sign of (t − t_switch).
#
#   The second structural decision is that there are TWO VEGF POOLS, not
#   one, and they read opposite sides of that collision:
#
#     VFRONT — VEGF at the growing vascular front.  Reads the LOCAL pO2
#              there (choroid minus the front region's own consumption).
#              Sets vessel growth (gate KG) and endothelial survival
#              (gate KSURV).  This is the pool hyperoxia collapses.
#     VEGFR  — VEGF in the vitreous.  Reads hypoxic intensity times the
#              AREA of avascular retina.  Sets neovascularisation (gate
#              KNV, Hill n=4) and plus disease.  This is the pool that
#              explodes after the switch.
#
#   A single well-mixed VEGF variable cannot do this: it would have to be
#   simultaneously suppressed by hyperoxia (to arrest growth in phase 1)
#   and elevated by it (to produce more NV in phase 2).  Splitting the
#   pool is what makes the sign flip representable at all.
#
#   Every number quoted below is produced by running this file
#   (run_all_rop()), or is stated arithmetic from published trial data
#   with the source named.
#
#   ---------------------------------------------------------------
#   AXIS 1 — THE PULSE OXIMETER IS NOT AN OXYGEN SENSOR, AND THE ERROR
#            IS BIGGER THAN THE TRIALS' EXPERIMENTAL CONTRAST.
#     With Hill n = 2.7 and P50 shifting from 19.0 mmHg (HbF-dominant)
#     to 26.6 mmHg (post-transfusion HbA), the SAME SpO2 maps to PaO2
#     values 1.40-fold apart at every point in the range.  Put the other
#     way round: a PaO2 of 50 mmHg reads as
#             SpO2 93.2%  in an untransfused preterm infant
#             SpO2 84.6%  after replacement with adult blood
#     — an 8.6-point offset.  The entire prescribed contrast of the
#     NeOProM trials (85–89% vs 91–95%) is 6 points.  The uncontrolled
#     drift of the haemoglobin switch across a NICU stay therefore moves
#     retinal oxygen by MORE than the intervention those trials tested.
#     Worse, dPaO2/dSpO2 is 1.3 mmHg per point at SpO2 87% and 45 mmHg
#     per point at 99%: the instrument loses all resolution exactly where
#     the therapeutic threshold of AXIS 2 lies.
#
#   AXIS 2 — A SQUARE-ROOT LAW PUTS A CEILING ON EVERY OXYGEN-BASED
#            PHASE 2 THERAPY, AND THE CEILING IS ANATOMICAL.
#     The avascular retina is a slab fed from one face.  Its oxygenated
#     depth is L = sqrt(2k·PCHOR/MO2): to double the reach you need FOUR
#     times the choroidal pO2.  Because MO2 and h both rise with
#     maturation, PCRIT climbs from 6.7 mmHg at PMA 28 wk to 93.7 mmHg
#     at PMA 36 wk.  At PMA 35 wk the model needs PaO2 = 94.6 mmHg to
#     re-oxygenate the full thickness, i.e. SpO2 98.7% in an untransfused
#     infant and 96.9% after transfusion.  STOP-ROP prescribed 96–99%,
#     which in HbF-rich blood spans PaO2 61.7–104.2 mmHg — the band
#     STRADDLES the requirement.  That is the model's account of why
#     STOP-ROP was a near miss (48% → 41%, adjusted OR 0.72,
#     95% CI 0.52–1.01) rather than either a success or a null: only the
#     upper part of the prescribed range clears the diffusion threshold,
#     and the oximeter cannot tell you which infants are in it.
#
#   AXIS 3 — REPRODUCING NeOProM REQUIRES AN ACHIEVED SEPARATION OF
#            ~2 SpO2 POINTS, NOT THE 6 THAT WERE PRESCRIBED.
#     One parameter (KNV) was calibrated so that the higher-target arm
#     shows 15.6% treatment-requiring ROP against NeOProM's observed
#     14.9%.  Nothing else about ROP was fitted.  Sweeping the achieved
#     separation then gives
#        separation   1.0    2.0    3.0    4.0    6.5 points
#        ROP RR      0.872  0.743  0.674  0.604  0.471
#     NeOProM observed RR 0.74.  The model lands on it at 2.0 points,
#     and at that separation ALSO predicts the low-arm absolute
#     incidence as 11.6% against 10.9% observed — a number nothing was
#     fitted to.  At the prescribed 6.5-point separation the model says
#     the effect should have been RR 0.47.  The mechanistic effect of
#     the POLICIES is roughly three times the effect the trials
#     measured, and the difference is exposure, not biology: an achieved
#     separation of ~2 points is what ~50–60% time-in-target-range with
#     overlapping distributions produces.  This is a falsifiable claim
#     about the achieved oximetry, not about ROP.
#
#   AXIS 4 — THE TWO STANDARD ANTI-VEGF DOSES ARE EQUIMOLAR, SO NOTHING
#            THAT DIFFERS BETWEEN THEM IS THE DOSE.
#     bevacizumab 0.625 mg / 149 kDa = 4.195 nmol
#     ranibizumab 0.200 mg /  48 kDa = 4.167 nmol   (ratio 1.007)
#     What differs is valency (2 sites vs 1 → capacity ratio 2.01),
#     vitreous half-life (9.82 d vs 7.19 d in adult human eyes), and —
#     the only large difference — the Fc domain.  Capacity × half-life
#     favours bevacizumab 2.75-fold ocularly; FcRn recycling favours it
#     41.7-fold systemically in this model, against the 35-fold serum
#     AUC ratio Avery et al. measured in adults.  One molecule is a
#     modestly better ocular agent and a far worse systemic one, and
#     both facts come from the same Fc.
#
#   AXIS 5 — THE DOSE IS ADULT-SCALED BUT ONLY ONE OF THE TWO RELEVANT
#            VOLUMES IS.
#     Ocular effect scales with dose / vitreous volume; systemic risk
#     scales with dose / body weight.  Preterm vitreous 1.1 mL vs adult
#     4.0 mL; preterm body 1 kg vs adult 70 kg.  Therefore
#        ocular concentration : 1.8× the adult exposure
#        systemic dose per kg : 35× the adult exposure
#     from the same "half the adult dose".  And in the vitreous the drug
#     is in ~229,000-fold molar excess over its target (3813 nM of
#     binding sites against 0.033 nM of VEGF at 1500 pg/mL).  Efficacy
#     is therefore logarithmic in dose — each halving costs one
#     half-life of coverage — while systemic exposure is strictly
#     linear.  That asymmetry, not caution, is the argument for
#     de-escalation.
#
#   AXIS 6 — THE PEDIG DOSE FLOOR IS A STOICHIOMETRIC LIMIT, NOT AN
#            AFFINITY LIMIT, AND THE MODEL PUTS IT WHERE IT WAS FOUND.
#     PEDIG de-escalated bevacizumab over a 312-fold range and found the
#     first failure only at the bottom: 0.031 mg 9/9, 0.016 mg 13/13,
#     0.008 mg 9/9, 0.004 mg 9/10, 0.002 mg 17/23.  Three constraints
#     could set that floor.  Computed for a 1.1 mL eye:
#       affinity  — 4-week trough falls to 0.227 nM (0.004 mg) and
#                   0.114 nM (0.002 mg) with the volume-scaled preterm
#                   half-life, i.e. 2–4× the reported Kd of 0.058 nM.
#       capacity  — total binding sites are 53.7 pmol (0.004 mg) and
#                   26.9 pmol (0.002 mg) against a cumulative 4-week
#                   vitreous VEGF production of 17.1 pmol if vitreous
#                   VEGF turns over with a 1 h half-life (2.9 pmol at
#                   6 h).  The margin falls to 1.6× exactly at the dose
#                   that failed.
#       duration  — coverage is logarithmic in dose, so the ladder buys
#                   only ~1 half-life per halving.
#     Both affinity and capacity bind within a factor of ~2 at the
#     observed floor, so this model CANNOT discriminate them, and says
#     so.  The measurement that would: the vitreous clearance rate of
#     VEGF in the preterm eye.  The capacity reading makes a
#     falsifiable prediction the affinity reading does not — the floor
#     must move UP in aggressive posterior ROP, because capacity is
#     consumed by VEGF production and APROP produces more.
#
#   AXIS 7 — MEASURED SERUM VEGF SUPPRESSION AFTER IVB IS TENFOLD LESS
#            THAN 1:1 BINDING PREDICTS, AND THE GAP IS THE COMPLEX.
#     Sato et al. measured, after 0.5 mg total IVB in infants, serum
#     bevacizumab 1214 ng/mL and serum free VEGF 269 pg/mL against a
#     1628 pg/mL baseline (16.5% of baseline).  Naive 1:1 binding at the
#     measured drug concentration and the published Kd predicts 1.4% of
#     baseline — twelvefold too much suppression.  The discrepancy is
#     not affinity: it is that drug-bound VEGF is protected from
#     clearance, so TOTAL serum VEGF rises while the bound FRACTION
#     stays high.  Setting complex clearance 17-fold slower than free
#     VEGF reproduces 252 pg/mL at day 14 against 269 observed.  This
#     makes a directly testable prediction that no published study
#     reports: TOTAL (bound + free) serum VEGF after IVB should be
#     roughly an order of magnitude ABOVE baseline while free VEGF is
#     below it.
#
#   AXIS 8 — LOW DOSE IS NOT A SAFETY COMPROMISE; IT IS ALSO THE BETTER
#            OCULAR STRATEGY, BECAUSE OF THE VEGF WINDOW.
#     Physiologic vascularisation needs VEGF above KG; neovascularisation
#     needs it above KNV.  A dose that drives the shared vitreous pool
#     below KG stops the normal vessels too, leaving persistent
#     avascular retina — the substrate for late reactivation.  In the
#     severe reference patient:
#        dose      weeks VEGF<KG   reactivation    residual avascular
#        0.625 mg      6.0           +8.1 wk            0.414
#        0.125 mg      4.4           +6.4 wk            0.342
#        0.031 mg      3.0           +5.0 wk            0.281
#        0.004 mg      1.4           +3.1 wk            0.203
#     Each fourfold dose reduction costs ~1.5 weeks of quiescence and
#     buys ~0.06 of retina that goes on to vascularise.  The high dose
#     does not merely carry systemic risk; it freezes the retina in the
#     avascular state that makes it need the drug again.
#
#   AXIS 9 — THE COMPETING RISK IS ARITHMETIC, AND IT IS DECISIVE.
#     From NeOProM's own absolute rates, lowering the target range per
#     1000 infants: 40 fewer ROP treatments, 28 more deaths, 23 more
#     cases of severe necrotising enterocolitis.  Carrying the 40
#     avoided treatments through ETROP's outcome rates for treated eyes
#     (9.1% unfavourable structural, 14.5% unfavourable acuity) gives
#     3.6 structural and 5.8 acuity outcomes avoided — that is
#     7.7 deaths per unfavourable structural outcome prevented, or
#     4.8 per unfavourable acuity outcome.  No parameter of this model
#     enters that calculation; it is the trial data divided.
#
#   ---------------------------------------------------------------
#   WHAT THIS MODEL GETS WRONG (reported, not repaired)
#
#   F1  It CANNOT reproduce BEAT-ROP's laser recurrence rate.  Observed:
#       bevacizumab 4% vs laser 22% before 54 wk PMA.  The model gives
#       laser 0% at every ablation completeness from 95% down to 45%
#       (residual avascular fraction 0.008 → 0.059, never enough to
#       re-cross KNV) and anti-VEGF 100%.  The ordering is inverted.
#       The diagnosis is structural and specific: laser failure is
#       SPATIALLY FOCAL — a skipped sector stays locally hypoxic — and a
#       single well-mixed avascular compartment cannot represent a
#       sector.  Reproducing 22% requires a sectorised avascular
#       compartment, which this model does not have.  Note also that the
#       model's near-universal late anti-VEGF reactivation is closer to
#       the long-follow-up literature than to BEAT-ROP's 54-week
#       primary endpoint.
#
#   F2  It CANNOT reproduce Mega Donna Mega (AA:DHA supplementation
#       halved severe ROP, 33.3% → 15.8%, aRR 0.50).  The model gives
#       RR 1.00.  The reason is quantitative and worth stating: the only
#       channel the model offers AA/DHA is oxidative-stress-driven
#       vaso-obliteration, and at calibrated ROS levels that channel is
#       at most 5.8% of the obliteration rate.  Even abolishing it
#       entirely could not produce the observed effect.  The AA/DHA
#       effect must therefore act through something absent here —
#       membrane phospholipid composition at the growth front,
#       resolvin/lipoxin mediators, or a direct effect on angiogenic
#       tone — and the model is evidence that it is NOT simple
#       antioxidant protection.
#
#   F3  It over-predicts the incidence of ANY ROP: 99% of the simulated
#       cohort reaches at least stage 1, against ~60–70% reported.  The
#       model has no infants who complete vascularisation before the
#       metabolic switch, because vasculogenesis is too slow relative to
#       maturation at the mild end of the cohort.  The
#       treatment-requiring endpoint — the one calibrated, and the one
#       every drug axis rests on — matches (15.6% vs 14.9%).
#
#   F4  The STOP-ROP subgroup finding is neither confirmed nor refuted.
#       The trial found benefit only without plus disease (46% → 32%)
#       and none with it (52% → 57%).  At the model's enrolment trigger
#       almost no infant yet has plus disease (25/621), and its
#       "threshold" surrogate is essentially unreachable in the
#       calibrated cohort, so the subgroup contrast cannot be evaluated.
#       With a Type-1 endpoint the model gives 23.2% → 13.6%
#       (RR 0.59), a LARGER phase-2 oxygen benefit than STOP-ROP's
#       OR 0.72 (CI up to 1.01).
#
#   F5  Serum free VEGF after ranibizumab transiently exceeds baseline
#       in the model (120% at day 14) because complexed VEGF accumulates
#       faster at low drug levels than binding removes free VEGF.  Avery
#       et al. measured a slight DECREASE (trough 14.4 vs 17 pg/mL
#       baseline).  The direction of the bevacizumab/ranibizumab
#       contrast is right; the ranibizumab absolute level overshoots.
#
#   F6  The achieved-separation inference (AXIS 3) rests on the model's
#       phase-1 oxygen sensitivity, which is calibrated to ROP incidence
#       levels rather than to any direct measurement of front-pool VEGF.
#       An independent measurement of retinal VEGF versus PaO2 in phase
#       1 would test it; none exists in human preterm eyes.
#
#   F7  Vitreous volume in the preterm eye (1.1 mL) and the exponent
#       scaling intravitreal half-life with eye volume (2/3) are
#       assumptions, not measurements, and AXIS 6's floor moves with
#       both.  Under the unscaled adult half-life the 4-week troughs are
#       15-fold higher and the affinity constraint no longer binds.
#
# ---------------------------------------------------------------------
#   CALIBRATION LEDGER — what was fitted, and to what
#     KNV, KPLV      : treatment-requiring ROP = 14.9% in the higher-
#                      target arm (NeOProM, PMID 29872859)
#     achieved sep.  : 2.0 SpO2 points, to NeOProM ROP RR 0.74
#     H0, AHYPOX     : death 17.1% and RR 1.17 (NeOProM)
#     FSYS, CLSYS    : serum bevacizumab 1214 ng/mL at day 14 after
#                      0.5 mg (Sato, PMID 21930258); serum t1/2 21 d and
#                      tmax 14 d (Kong, PMID 25613938)
#     KDVSB          : serum free VEGF 269 pg/mL at day 14 (Sato)
#     CLSYS(ranib)   : bev/ranibizumab serum AUC ratio 35 (Avery,
#                      PMID 25001321)
#     KMYOP          : zone I laser spherical equivalent −8.44 D
#                      (Geloneck, PMID 25103848)
#     KDIFF, MO2,    : physiological literature values, NOT fitted
#     P50A/F, NHILL,   (Krogh constant, retinal O2 consumption,
#     THOC, KDOC       oxyhaemoglobin dissociation, human intravitreal
#                      half-lives, VEGF binding affinities)
#   PREDICTIONS (nothing fitted to them)
#     low-arm ROP incidence 11.6% vs 10.9% observed
#     serum free VEGF 12.7% of baseline at 0.625 mg vs 80.8% at
#       0.031 mg — the systemic safety case for de-escalation
#     bev myopia −1.78 D vs −1.51 D observed (Geloneck)
#     PaO2 requirement at PMA 35 straddled by STOP-ROP's 96–99% band
# =====================================================================

library(mrgsolve)
suppressMessages(library(dplyr))

rop_code <- '
$PROB
# Retinopathy of Prematurity QSP model (36 ODEs)
# supply/demand collision in the avascular retina; two VEGF pools;
# intravitreal anti-VEGF with target-mediated binding and systemic escape.

$PLUGIN autodec

$PARAM @annotated
// ---------------- patient ------------------------------------------
GA      : 25.0  : Gestational age at birth (weeks)
BW      : 700   : Birth weight (g)

// ---------------- oxygen carriage ---------------------------------
P50A    : 26.6  : P50 of adult haemoglobin (mmHg)
P50F    : 19.0  : P50 of fetal haemoglobin (mmHg)
NHILL   : 2.7   : Hill coefficient of the oxyhaemoglobin dissociation curve
FHBF0   : 0.90  : Fetal haemoglobin fraction at birth (-)
KHBF    : 0.011 : Endogenous gamma-to-beta globin switch rate (1/day)
CHOR    : 0.85  : Choroidal pO2 as a fraction of arterial pO2 (-)
SPO2TGT : 93.5  : Achieved SpO2 midpoint while oxygen-dependent (%)
SPO2RA  : 98.0  : SpO2 once in room air (%)
PMAWEAN : 36.0  : PMA at which supplemental oxygen is weaned (weeks)
TAUWEAN : 1.5   : Weaning transition width (weeks)

// ---------------- retinal oxygen diffusion ------------------------
KDIFF   : 4.7e-10 : Krogh diffusion constant (mL O2 per cm per s per mmHg)
MO2MIN  : 0.30e-4 : Peripheral retinal O2 consumption floor (mL O2/mL/s)
MO2MAX  : 3.20e-4 : Peripheral retinal O2 consumption span (mL O2/mL/s)
PMA50M  : 32.5    : PMA at half-maximal retinal O2 consumption (weeks)
TAUM    : 1.6     : Width of the metabolic maturation sigmoid (weeks)
HRET0   : 0.0105  : Avascular retinal thickness at baseline (cm)
HRETMX  : 0.0185  : Mature avascular retinal thickness (cm)
PMA50H  : 33.0    : PMA at half-maximal retinal thickness (weeks)
TAUH    : 2.5     : Width of the thickness maturation sigmoid (weeks)
FTHK    : 0.60    : Front-region thickness as a fraction of h (-)
PFMIN   : 2.0     : Floor on front-region pO2 (mmHg)

// ---------------- hypoxia signalling / VEGF pools -----------------
KO2     : 25.0  : Oxygen tension for half-maximal HIF stabilisation (mmHg)
KO2F    : 25.0  : Same, at the vascular front (mmHg)
KINH    : 2.0   : HIF-1alpha formation rate (1/day)
KOUTH   : 2.0   : HIF-1alpha degradation rate (1/day)
KSV     : 224.6 : Constitutive vitreous VEGF synthesis (pg/mL/day)
AHIF    : 1.703 : HIF amplification of vitreous VEGF synthesis (-)
NHIFV   : 4.0   : Cooperativity of the HIF-to-VEGF map (-)
AVREF   : 0.25  : Reference avascular fraction for area scaling (-)
AINFL   : 0.6   : Inflammatory amplification of VEGF synthesis (-)
KDV     : 1.0   : Free vitreous VEGF elimination (1/day)
KDVB    : 0.06  : Drug-bound vitreous VEGF elimination (1/day)
KSVF    : 48.9  : Constitutive front-pool VEGF set point (pg/mL)
AHIFF   : 239.0 : Hypoxic amplification of the front pool (-)
NHIFF   : 4.0   : Cooperativity of the front-pool oxygen response (-)
TAUVF   : 0.5   : Front-pool equilibration time constant (day)
KSE     : 1.0   : Retinal erythropoietin synthesis (1/day)
KDE     : 1.0   : Retinal erythropoietin elimination (1/day)
AHIFE   : 2.5   : HIF amplification of retinal erythropoietin (-)
KSA     : 1.0   : Angiopoietin-2 synthesis (1/day)
KDA     : 1.0   : Angiopoietin-2 elimination (1/day)
AHIFA   : 1.5   : HIF amplification of angiopoietin-2 (-)

// ---------------- vascularisation / obliteration ------------------
KGROW   : 0.070 : Maximal retinal vascularisation rate (1/day)
KG      : 300.0 : Front-pool VEGF gate for vessel growth (pg/mL)
KIGF    : 25.0  : Serum IGF-1 for half-maximal permissiveness (ng/mL)
KOBLV   : 0.075 : Vaso-obliteration rate via VEGF withdrawal (1/day)
KOBLR   : 0.035 : Vaso-obliteration rate via oxidative stress (1/day)
KSURV   : 180.0 : Front-pool VEGF gate for endothelial survival (pg/mL)
PMAOBL  : 30.5  : PMA at half-loss of obliteration susceptibility (weeks)
TAUOBL  : 1.5   : Width of the susceptibility window (weeks)
PMAGSTOP: 45.0  : PMA at which vasculogenesis intrinsically stops (weeks)
TAUGSTOP: 2.0   : Width of that stop (weeks)

// ---------------- ridge / NV / plus disease -----------------------
KSH     : 0.30  : Ridge formation rate (1/day)
KDSH    : 0.10  : Ridge regression rate (1/day)
KSHV    : 880.0 : Vitreous VEGF threshold for ridge formation (pg/mL)
NSH     : 4.0   : Cooperativity of ridge formation (-)
KNVR    : 0.85  : Neovascularisation formation rate (1/day)
KNV     : 1170  : Vitreous VEGF threshold for neovascularisation (pg/mL)
NNV     : 4.0   : Cooperativity of the neovascularisation threshold (-)
NVMAX   : 1.0   : Neovascular capacity per unit avascular retina (-)
KNVREG  : 0.30  : Neovascular regression rate (1/day)
KPL     : 0.30  : Plus-disease formation rate (1/day)
KPLV    : 1300  : Vitreous VEGF threshold for plus disease (pg/mL)
NPL     : 4.0   : Cooperativity of plus disease (-)
KPLO    : 0.16  : Plus-disease resolution rate (1/day)

// ---------------- fibrosis / detachment ---------------------------
TAUVS   : 3.0   : Time constant of the slow VEGF trace (day)
KTGF    : 0.55  : TGF-beta formation rate (1/day)
KDTGF   : 0.35  : TGF-beta elimination rate (1/day)
KCR     : 0.0025: Crunch sensitivity to abrupt VEGF withdrawal (per pg/mL)
KFIB    : 0.32  : Fibrovascular membrane formation rate (1/day)
KDFIB   : 0.05  : Fibrovascular membrane resolution rate (1/day)
KDET    : 0.12  : Tractional detachment formation rate (1/day)
KDDET   : 0.02  : Tractional detachment resolution rate (1/day)

// ---------------- growth / IGF-1 / inflammation -------------------
KSI     : 3.2   : Serum IGF-1 synthesis (ng/mL/day)
KDI     : 0.10  : Serum IGF-1 elimination (1/day)
KIINF   : 0.45  : Inflammatory suppression of IGF-1 (-)
NUTR    : 1.0   : Nutritional adequacy index (-)
KWT     : 0.016 : Logistic weight gain rate (1/day)
WTMAX   : 4200  : Asymptotic weight (g)
KDINF   : 0.25  : Resolution of systemic inflammation (1/day)

// ---------------- lung / oxidative stress ------------------------
KBPD    : 0.055 : Hyperoxic lung injury accrual (1/day)
KDBPD   : 0.020 : Lung injury resolution (1/day)
KROS    : 0.030 : Oxidative stress formation (1/day)
KDROS   : 0.45  : Oxidative stress clearance (1/day)
PROS    : 65.0  : PaO2 above which oxidative stress accrues (mmHg)
KAADHA  : 0.0   : AA:DHA enteral supplementation switch (0/1)

// ---------------- drug (defaults: bevacizumab) --------------------
MWD     : 149.0 : Drug molecular weight (kDa)
VALD    : 2.0   : VEGF binding sites per drug molecule (-)
KDOC    : 0.30  : Ocular drug-VEGF dissociation constant (nM)
KDSY    : 0.30  : Systemic drug-VEGF dissociation constant (nM)
THOC    : 9.82  : Intravitreal half-life in the adult human eye (day)
OCEXP   : 0.667 : Exponent scaling ocular half-life with vitreous volume (-)
VVITAD  : 4.0   : Adult vitreous volume (mL)
VVIT0   : 1.10  : Preterm vitreous volume (mL)
VVITMX  : 2.20  : Vitreous volume at the end of follow-up (mL)
PMA50V  : 40.0  : PMA at half-maximal vitreous volume growth (weeks)
TAUV    : 8.0   : Width of vitreous volume growth (weeks)
KAQ     : 6.0   : Aqueous transit rate constant (1/day)
FSYS    : 0.44  : Systemic bioavailability of ocular efflux (-)
CLSYS   : 3.60  : Systemic clearance (mL/day)
VC      : 50.0  : Central volume of distribution (mL)
KPT     : 0.35  : Central-to-peripheral rate constant (1/day)
KTP     : 0.30  : Peripheral-to-central rate constant (1/day)
VS0     : 1628  : Baseline serum VEGF (pg/mL)
KDVS    : 1.0   : Free serum VEGF elimination (1/day)
KDVSB   : 0.06  : Drug-bound serum VEGF elimination (1/day)

// ---------------- outcomes ---------------------------------------
KMYOP   : 16.0  : Dioptres of myopia per unit ablated retinal fraction (D)
KMYOPB  : -1.2  : Baseline spherical equivalent (D)
KMYOPR  : 0.010 : ROP-associated myopic drift (D/day per unit NV)
KVFL    : 0.85  : Visual field loss per unit ablated retinal fraction (-)
H0      : 0.001122 : Baseline daily mortality hazard (1/day)
AHYPOX  : 0.150 : Hazard amplification by hypoxaemia (-)
ABPD    : 0.0022: Hazard amplification by lung injury (-)
PHYPOX  : 55.0  : PaO2 below which mortality hazard rises (mmHg)

// ---------------- interventions ----------------------------------
LASERT  : 1e6   : Time of laser photocoagulation (day)
LASDUR  : 0.25  : Duration of the laser event (day)
LFRAC   : 0.85  : Fraction of avascular retina ablated (-)
O2SUPP  : 0     : Supplemental oxygen switch (0/1)
O2SUPPT : 1e6   : Time supplemental oxygen starts (day)
SPO2SUP : 97.5  : SpO2 target during supplemental oxygen (%)
TXT1    : 1e6   : Time of first red-cell transfusion (day)
TXT2    : 1e6   : Time of second red-cell transfusion (day)
TXDUR   : 0.25  : Duration of a transfusion event (day)
FTX     : 0.20  : Fraction of fetal haemoglobin replaced per transfusion (-)
SEPT    : 1e6   : Time of a sepsis episode (day)
SEPDUR  : 0.50  : Duration of the sepsis stimulus (day)
SEPAMT  : 1.5   : Magnitude of the inflammatory stimulus (-)

$CMT @annotated
FHBF    : Fetal haemoglobin fraction (-)
HIF     : HIF-1alpha activity in the avascular retina (-)
VEGFR   : Total vitreous VEGF (pg/mL)
EPOR    : Retinal erythropoietin (-)
ANG2    : Angiopoietin-2 (-)
VASC    : Vascularised fraction of the retina (-)
ABLA    : Laser-ablated fraction of the retina (-)
SHUNT   : Ridge / mesenchymal shunt index (-)
NV      : Extraretinal neovascular mass (-)
PLUS    : Plus-disease index (-)
VSLOW   : Slow trace of free vitreous VEGF (pg/mL)
TGFB    : TGF-beta / fibrogenic drive (-)
FIBRO   : Fibrovascular membrane (-)
DETACH  : Tractional detachment index (-)
IGF1    : Serum IGF-1 (ng/mL)
WT      : Body weight (g)
INFL    : Systemic inflammation index (-)
BPD     : Lung injury index (-)
DVIT    : Free drug in the vitreous (nmol)
DAQ     : Drug in the aqueous transit compartment (nmol)
DCEN    : Drug in the central compartment (nmol)
DPER    : Drug in the peripheral compartment (nmol)
VEGFS   : Total serum VEGF (pg/mL)
AUCDS   : Cumulative systemic drug exposure (nM*day)
AUCVS   : Cumulative serum VEGF suppression (pg/mL*day)
HAZD    : Cumulative mortality hazard (-)
HYPBUR  : Cumulative hypoxaemia burden (day)
HYPRBUR : Cumulative hyperoxia burden (day)
MYOP    : Spherical equivalent refraction (D)
VFLOSS  : Peripheral visual field loss (-)
AVBUR   : Cumulative avascular hypoxia burden (day)
NVPEAK  : Running peak of neovascular mass (-)
VADEF   : Visual acuity deficit index (logMAR)
ROS     : Oxidative stress index (-)
RETPR   : Reactivation pressure (-)
VFRONT  : VEGF at the growing vascular front (pg/mL)

$MAIN
if (NEWIND <= 1) {
  FHBF_0   = FHBF0;
  HIF_0    = 0.30;
  VEGFR_0  = 400.0;
  EPOR_0   = 0.50;
  ANG2_0   = 0.50;
  // vascularised retinal AREA at birth: radius linear in PMA from 16 wk,
  // reaching the temporal periphery at 40 wk; area goes as radius squared
  rbirth   = (GA - 16.0) / 24.0;
  if (rbirth < 0.02) rbirth = 0.02;
  if (rbirth > 1.00) rbirth = 1.00;
  VASC_0   = rbirth * rbirth;
  IGF1_0   = 22.0;
  WT_0     = BW;
  VEGFS_0  = VS0;
  MYOP_0   = KMYOPB;
  VFRONT_0 = 300.0;
}

$ODE
// ================= algebraic block ================================
PMA   = GA + SOLVERTIME / 7.0;

// ---- oxygen policy: target while oxygen-dependent, then room air
o2dep = 1.0 - 1.0 / (1.0 + exp(-(PMA - PMAWEAN) / TAUWEAN));
spo2  = SPO2TGT * o2dep + SPO2RA * (1.0 - o2dep);
if (O2SUPP > 0.5 && SOLVERTIME >= O2SUPPT) spo2 = SPO2SUP;
if (spo2 > 99.4) spo2 = 99.4;

// ---- oxyhaemoglobin dissociation: SpO2 -> PaO2 depends on P50
P50h  = P50A + (P50F - P50A) * FHBF;
Ssat  = spo2 / 100.0;
PAO2  = P50h * pow(Ssat / (1.0 - Ssat), 1.0 / NHILL);
PCHOR = CHOR * PAO2;

// ---- maturation of retinal metabolic demand and thickness
MO2   = MO2MIN + MO2MAX / (1.0 + exp(-(PMA - PMA50M) / TAUM));
HRET  = HRET0 + (HRETMX - HRET0) / (1.0 + exp(-(PMA - PMA50H) / TAUH));

// ---- Krogh slab: choroid-fed avascular retina, no-flux inner boundary
PCRIT  = MO2 * HRET * HRET / (2.0 * KDIFF);
LPEN   = sqrt(2.0 * KDIFF * PCHOR / MO2) * 1e4;          // um
PINNER = PCHOR - PCRIT;
if (PINNER < 0.0) PINNER = 0.0;
FHYP0  = KO2 / (KO2 + PINNER);

// ---- local pO2 at the growing vascular front: the front region is
//      thinner, so its own consumption subtracts less -- and it is the
//      RISE in that consumption with maturation that later shields the
//      front from systemic hyperoxia
HF     = FTHK * HRET;
PFRONT = 0.90 * PCHOR - MO2 * HF * HF / (2.0 * KDIFF);
if (PFRONT < PFMIN) PFRONT = PFMIN;
FHYPF  = KO2F / (KO2F + PFRONT);
VFSS   = KSVF * (1.0 + AHIFF * pow(FHYPF, NHIFF));

// ---- eye growth and intravitreal drug
VVIT   = VVIT0 + (VVITMX - VVIT0) / (1.0 + exp(-(PMA - PMA50V) / TAUV));
THOCE  = THOC * pow(VVIT / VVITAD, OCEXP);
CVIT   = 1000.0 * DVIT / VVIT;                            // nM
SITE   = VALD * CVIT;
VEGFF  = VEGFR * KDOC / (KDOC + SITE);                    // free vitreous VEGF
VFR    = VFRONT * KDOC / (KDOC + SITE);                   // free front VEGF
CSYS   = 1000.0 * DCEN / VC;                              // nM
VEGFSF = VEGFS * KDSY / (KDSY + VALD * CSYS);             // free serum VEGF

// ---- retinal geometry: the avascular area is the VEGF source term
AVASC  = 1.0 - VASC - ABLA;
if (AVASC < 0.0) AVASC = 0.0;
VSRC   = (AVASC / AVREF) * pow(FHYP0, NHIFV);

// ---- event pulses (finite duration, so everything stays in the ODE)
klas = 0.0;
if (SOLVERTIME >= LASERT && SOLVERTIME < LASERT + LASDUR)
  klas = -log(1.0 - LFRAC) / LASDUR;
ktx = 0.0;
if (SOLVERTIME >= TXT1 && SOLVERTIME < TXT1 + TXDUR)
  ktx = -log(1.0 - FTX) / TXDUR;
if (SOLVERTIME >= TXT2 && SOLVERTIME < TXT2 + TXDUR)
  ktx = -log(1.0 - FTX) / TXDUR;
ksep = 0.0;
if (SOLVERTIME >= SEPT && SOLVERTIME < SEPT + SEPDUR)
  ksep = SEPAMT / SEPDUR;

// ================= differential equations =========================
// 1 fetal haemoglobin: endogenous switch plus transfusion dilution
dxdt_FHBF = -KHBF * FHBF - ktx * FHBF;

// 2 HIF-1alpha tracks hypoxia INTENSITY in the avascular retina
dxdt_HIF = KINH * FHYP0 - KOUTH * HIF;

// 3 vitreous VEGF: source is intensity x AREA; bound VEGF clears slowly
dxdt_VEGFR = KSV * (1.0 + AHIF * VSRC) * (1.0 + AINFL * INFL)
             - KDV * VEGFF - KDVB * (VEGFR - VEGFF);

// 4-5 retinal erythropoietin and angiopoietin-2
dxdt_EPOR = KSE * (1.0 + AHIFE * HIF) - KDE * EPOR;
dxdt_ANG2 = KSA * (1.0 + AHIFA * HIF) - KDA * ANG2;

// 6 vascularised fraction: VEGF-gated, IGF-1-permissive growth at the
//   front, minus obliteration of vessels that lose their survival signal
gIGF  = IGF1 / (KIGF + IGF1);
fGROW = VFR / (KG + VFR);
gstop = 1.0 - 1.0 / (1.0 + exp(-(PMA - PMAGSTOP) / TAUGSTOP));
grow  = KGROW * gIGF * fGROW * gstop * AVASC;
sobl  = 1.0 - 1.0 / (1.0 + exp(-(PMA - PMAOBL) / TAUOBL));
fwd   = KSURV / (KSURV + VFR);
obl   = sobl * VASC * (KOBLV * fwd + KOBLR * ROS);
dxdt_VASC = grow - obl;

// 7 laser-ablated fraction: deletes the VEGF source outright
dxdt_ABLA = klas * AVASC;

// 8 ridge / mesenchymal shunt (stage 1-2)
fSH = pow(VEGFF, NSH) / (pow(KSHV, NSH) + pow(VEGFF, NSH));
dxdt_SHUNT = KSH * AVASC * fSH * 4.0 * (1.0 - SHUNT) - KDSH * SHUNT;

// 9 extraretinal neovascularisation: steep VEGF threshold, capacity set
//   by the amount of avascular retina still present
fNV = pow(VEGFF, NNV) / (pow(KNV, NNV) + pow(VEGFF, NNV));
rnv = KNVR * fNV * (AVASC * NVMAX - NV) - KNVREG * NV * (1.0 - fNV);
dxdt_NV = rnv;

// 10 plus disease
fPL = pow(VEGFF, NPL) / (pow(KPLV, NPL) + pow(VEGFF, NPL));
dxdt_PLUS = KPL * fPL * (1.0 + 0.4 * ANG2) * (1.0 - PLUS) - KPLO * PLUS;

// 11-12 slow VEGF trace and the crunch-driven fibrogenic signal
dxdt_VSLOW = (VEGFF - VSLOW) / TAUVS;
crunch = VSLOW - VEGFF;
if (crunch < 0.0) crunch = 0.0;
dxdt_TGFB = KTGF * (0.2 + KCR * crunch) - KDTGF * TGFB;

// 13-14 fibrovascular contraction and tractional detachment
dxdt_FIBRO  = KFIB * NV * TGFB * (1.0 - FIBRO) - KDFIB * FIBRO;
dxdt_DETACH = KDET * FIBRO * NV * (1.0 - DETACH) - KDDET * DETACH;

// 15-17 serum IGF-1, weight, systemic inflammation
dxdt_IGF1 = KSI * NUTR * (1.0 + 0.020 * (PMA - 25.0)) / (1.0 + KIINF * INFL)
            - KDI * IGF1;
dxdt_WT   = KWT * WT * (1.0 - WT / WTMAX) * NUTR;
dxdt_INFL = ksep - KDINF * INFL;

// 18 hyperoxic lung injury
fio2eq = (PAO2 - 40.0) / 260.0 + 0.21;
if (fio2eq < 0.21) fio2eq = 0.21;
if (fio2eq > 1.00) fio2eq = 1.00;
exc = fio2eq - 0.30;
if (exc < 0.0) exc = 0.0;
dxdt_BPD = KBPD * exc * 10.0 - KDBPD * BPD;

// 19-22 drug: vitreous depot -> aqueous -> plasma (two compartments)
kel = log(2.0) / THOCE;
dxdt_DVIT = -kel * DVIT;
dxdt_DAQ  = kel * DVIT - KAQ * DAQ;
dxdt_DCEN = FSYS * KAQ * DAQ - CLSYS / VC * DCEN - KPT * DCEN + KTP * DPER;
dxdt_DPER = KPT * DCEN - KTP * DPER;

// 23 serum VEGF: complexed VEGF is protected from clearance, which is
//    why measured free VEGF falls far less than 1:1 binding predicts
dxdt_VEGFS = KDVS * VS0 - KDVS * VEGFSF - KDVSB * (VEGFS - VEGFSF);

// 24-25 exposure metrics
dxdt_AUCDS = CSYS;
supp = VS0 - VEGFSF;
if (supp < 0.0) supp = 0.0;
dxdt_AUCVS = supp;

// 26 cumulative mortality hazard (competing risk)
hyp = (PHYPOX - PAO2) / PHYPOX;
if (hyp < 0.0) hyp = 0.0;
dxdt_HAZD = H0 * (1.0 + AHYPOX * hyp * 100.0 + ABPD * BPD * 100.0);

// 27-28 oxygen burden trackers
dxdt_HYPBUR = hyp;
hyper = (PAO2 - 80.0) / 80.0;
if (hyper < 0.0) hyper = 0.0;
dxdt_HYPRBUR = hyper;

// 29-30 refraction and visual field: laser ablation drives both
dxdt_MYOP   = -KMYOP * klas * AVASC - KMYOPR * NV;
dxdt_VFLOSS = KVFL * klas * AVASC;

// 31 cumulative avascular hypoxia burden
dxdt_AVBUR = AVASC * FHYP0;

// 32 running peak of neovascular mass
pk = rnv;
if (pk < 0.0) pk = 0.0;
dxdt_NVPEAK = pk;

// 33 visual acuity deficit accrual
dxdt_VADEF = 0.004 * DETACH + 0.0010 * FIBRO;

// 34 oxidative stress
rosdrv = (PAO2 - PROS) / PROS;
if (rosdrv < 0.0) rosdrv = 0.0;
dxdt_ROS = KROS * rosdrv * (1.0 - 0.45 * KAADHA) - KDROS * ROS;

// 35 reactivation pressure: avascular retina x approach to the NV gate
dxdt_RETPR = AVASC * fNV - 0.25 * RETPR;

// 36 front-pool VEGF equilibrates fast with local pO2
dxdt_VFRONT = (VFSS - VFRONT) / TAUVF;

$TABLE
PMAo   = GA + TIME / 7.0;
o2depo = 1.0 - 1.0 / (1.0 + exp(-(PMAo - PMAWEAN) / TAUWEAN));
spo2o  = SPO2TGT * o2depo + SPO2RA * (1.0 - o2depo);
if (O2SUPP > 0.5 && TIME >= O2SUPPT) spo2o = SPO2SUP;
if (spo2o > 99.4) spo2o = 99.4;
P50o   = P50A + (P50F - P50A) * FHBF;
So     = spo2o / 100.0;
PAO2o  = P50o * pow(So / (1.0 - So), 1.0 / NHILL);
PCHORo = CHOR * PAO2o;
MO2o   = MO2MIN + MO2MAX / (1.0 + exp(-(PMAo - PMA50M) / TAUM));
HRETo  = HRET0 + (HRETMX - HRET0) / (1.0 + exp(-(PMAo - PMA50H) / TAUH));
PCRITo = MO2o * HRETo * HRETo / (2.0 * KDIFF);
LPENo  = sqrt(2.0 * KDIFF * PCHORo / MO2o) * 1e4;
PINNo  = PCHORo - PCRITo;
if (PINNo < 0.0) PINNo = 0.0;
HFo    = FTHK * HRETo;
PFRONTo = 0.90 * PCHORo - MO2o * HFo * HFo / (2.0 * KDIFF);
if (PFRONTo < PFMIN) PFRONTo = PFMIN;
VVITo  = VVIT0 + (VVITMX - VVIT0) / (1.0 + exp(-(PMAo - PMA50V) / TAUV));
CVITo  = 1000.0 * DVIT / VVITo;
VEGFFo = VEGFR * KDOC / (KDOC + VALD * CVITo);
VFRo   = VFRONT * KDOC / (KDOC + VALD * CVITo);
CSYSo  = 1000.0 * DCEN / VC;
CSYSNG = CSYSo * MWD;
VEGFSFo = VEGFS * KDSY / (KDSY + VALD * CSYSo);
AVASCo = 1.0 - VASC - ABLA;
if (AVASCo < 0.0) AVASCo = 0.0;

// ICROP zone from the vascularised retinal area
ZONE = 3.0;
if (VASC < 0.75) ZONE = 2.0;
if (VASC < 0.30) ZONE = 1.0;

// ICROP stage
STAGE = 0.0;
if (SHUNT  > 0.15) STAGE = 1.0;
if (SHUNT  > 0.50) STAGE = 2.0;
if (NV     > 0.25) STAGE = 3.0;
if (DETACH > 0.15) STAGE = 4.0;
if (DETACH > 0.55) STAGE = 5.0;

PLUSD = (PLUS > 0.50) ? 1.0 : 0.0;

// ETROP type 1: zone I any stage with plus, zone I stage 3, or
// zone II stage 2-3 with plus
TYPE1 = 0.0;
if (ZONE < 1.5 && PLUSD > 0.5) TYPE1 = 1.0;
if (ZONE < 1.5 && STAGE > 2.5) TYPE1 = 1.0;
if (ZONE > 1.5 && ZONE < 2.5 && STAGE > 1.5 && PLUSD > 0.5) TYPE1 = 1.0;

// CRYO-ROP threshold surrogate
THRESH = (NV > 0.35 && PLUSD > 0.5) ? 1.0 : 0.0;

PDEATH = 1.0 - exp(-HAZD);

$CAPTURE @annotated
PMAo    : Postmenstrual age (weeks)
spo2o   : SpO2 (%)
PAO2o   : Arterial pO2 (mmHg)
PCHORo  : Choroidal pO2 (mmHg)
PCRITo  : Choroidal pO2 required for full-thickness oxygenation (mmHg)
PINNo   : Inner avascular retinal pO2 (mmHg)
LPENo   : Oxygen penetration depth into the avascular retina (um)
PFRONTo : pO2 at the growing vascular front (mmHg)
P50o    : Haemoglobin P50 (mmHg)
MO2o    : Retinal oxygen consumption (mL O2/mL/s)
VEGFFo  : Free vitreous VEGF (pg/mL)
VFRo    : Free front-pool VEGF (pg/mL)
AVASCo  : Avascular retinal fraction (-)
VVITo   : Vitreous volume (mL)
CVITo   : Free intravitreal drug (nM)
CSYSo   : Free plasma drug (nM)
CSYSNG  : Free plasma drug (ng/mL)
VEGFSFo : Free serum VEGF (pg/mL)
ZONE    : ICROP zone (1-3)
STAGE   : ICROP stage (0-5)
PLUSD   : Plus disease present (0/1)
TYPE1   : ETROP type 1 ROP (0/1)
THRESH  : CRYO-ROP threshold surrogate (0/1)
PDEATH  : Cumulative probability of death (-)
'

# ---------------------------------------------------------------------
#  build
# ---------------------------------------------------------------------
mod_rop <- mcode("rop", rop_code)

# ---------------------------------------------------------------------
#  drug parameter sets
# ---------------------------------------------------------------------
DRUG <- list(
  bevacizumab = list(MWD = 149.0, VALD = 2.0, KDOC = 0.30, KDSY = 0.30,
                     THOC = 9.82, CLSYS = 3.60, FSYS = 0.44),
  # ranibizumab: 48 kDa Fab, no Fc, so no FcRn recycling. CLSYS set from
  # the 35-fold bevacizumab/ranibizumab serum AUC ratio (Avery 2014).
  ranibizumab = list(MWD = 48.0, VALD = 1.0, KDOC = 0.10, KDSY = 0.10,
                     THOC = 7.19, CLSYS = 156.0, FSYS = 0.44),
  # aflibercept: 115 kDa Fc-fusion decoy receptor, sub-picomolar affinity
  aflibercept = list(MWD = 115.0, VALD = 2.0, KDOC = 0.0005, KDSY = 0.0005,
                     THOC = 9.00, CLSYS = 9.00, FSYS = 0.44)
)

mg_to_nmol <- function(dose_mg, mw_kda) dose_mg * 1e-3 / (mw_kda * 1e3) * 1e9

#' Simulate one infant.
#' @param dose_mg intravitreal dose in mg (0 = none)
#' @param drug one of names(DRUG)
#' @param tdose time of injection (days after birth)
sim_rop <- function(..., dose_mg = 0, drug = "bevacizumab", tdose = NA,
                    end = 175, delta = 0.5) {
  p <- list(...)
  if (dose_mg > 0) p <- modifyList(p, DRUG[[drug]])
  m <- if (length(p)) update(mod_rop, param = p) else mod_rop
  ev <- if (dose_mg > 0 && !is.na(tdose)) {
    ev(time = tdose, amt = mg_to_nmol(dose_mg, DRUG[[drug]]$MWD), cmt = "DVIT")
  } else {
    ev(time = 0, amt = 0, cmt = "DVIT")
  }
  mrgsim_df(m, events = ev, end = end, delta = delta, atol = 1e-10, rtol = 1e-8)
}

# ---------------------------------------------------------------------
#  SCENARIOS
# ---------------------------------------------------------------------

# S1 — natural history and the oxygen sign flip -----------------------
scenario_oxygen_sweep <- function(targets = c(85, 87, 89, 91, 93, 95, 97, 99)) {
  bind_rows(lapply(targets, function(tg) {
    d <- sim_rop(SPO2TGT = tg)
    j <- which.max(d$VEGFFo)
    tibble(SpO2 = tg,
           PFRONT_30wk   = approx(d$PMAo, d$PFRONTo, 30)$y,
           VFRONT_30wk   = approx(d$PMAo, d$VFRo, 30)$y,
           VASC_33wk     = approx(d$PMAo, d$VASC, 33)$y,
           switch_PMA    = d$PMAo[j],
           AVASC_switch  = d$AVASCo[j],
           peak_VEGF     = max(d$VEGFFo),
           peak_NV       = max(d$NV),
           peak_PLUS     = max(d$PLUS),
           type1         = max(d$TYPE1),
           max_stage     = max(d$STAGE))
  }))
}

# S2 — supplemental oxygen at prethreshold (STOP-ROP design) ----------
scenario_supplemental_o2 <- function(spo2sup = c(94, 96, 97.5, 99),
                                     enrol_pma = 34.5, GA = 24, NUTR = 0.65) {
  base <- sim_rop(GA = GA, NUTR = NUTR, SPO2TGT = 95)
  tenr <- (enrol_pma - GA) * 7
  bind_rows(lapply(c(NA, spo2sup), function(s) {
    d <- if (is.na(s)) base else
      sim_rop(GA = GA, NUTR = NUTR, SPO2TGT = 95,
              O2SUPP = 1, O2SUPPT = tenr, SPO2SUP = s)
    post <- d$time >= tenr
    tibble(arm = if (is.na(s)) "conventional" else sprintf("supplemental %.1f%%", s),
           PaO2_after   = approx(d$PMAo, d$PAO2o, enrol_pma + 1)$y,
           PCRIT_after  = approx(d$PMAo, d$PCRITo, enrol_pma + 1)$y,
           reoxygenated = as.numeric(approx(d$PMAo, d$PINNo, enrol_pma + 1)$y > 0),
           peak_NV_post = max(d$NV[post]),
           type1_post   = max(d$TYPE1[post]),
           BPD          = max(d$BPD))
  }))
}

# S3 — the anti-VEGF dose ladder -------------------------------------
scenario_dose_ladder <- function(doses = c(0.625, 0.25, 0.125, 0.063, 0.031,
                                           0.016, 0.008, 0.004, 0.002),
                                 drug = "bevacizumab", tdose = 60) {
  bind_rows(lapply(doses, function(dm) {
    d <- sim_rop(dose_mg = dm, drug = drug, tdose = tdose, end = 200, delta = 0.25)
    f <- function(day, col) approx(d$time, d[[col]], tdose + day)$y
    tibble(drug = drug, dose_mg = dm,
           nmol           = mg_to_nmol(dm, DRUG[[drug]]$MWD),
           sites_pmol     = mg_to_nmol(dm, DRUG[[drug]]$MWD) * DRUG[[drug]]$VALD * 1000,
           vitreous_nM_d1 = f(1, "CVITo"),
           vitreous_nM_d28 = f(28, "CVITo"),
           serum_ngmL_d14 = f(14, "CSYSNG"),
           serumVEGF_pct_d14 = 100 * f(14, "VEGFSFo") / 1628,
           serumVEGF_pct_d56 = 100 * f(56, "VEGFSFo") / 1628,
           days_VEGF_below_half = sum(d$VEGFSFo < 0.5 * 1628 & d$time > tdose) * 0.25,
           serum_AUC_nM_d = max(d$AUCDS))
  }))
}

# S4 — laser versus anti-VEGF in a treatment-requiring eye -----------
scenario_treatment <- function(GA = 24, NUTR = 0.65, SPO2TGT = 95) {
  sev <- list(GA = GA, NUTR = NUTR, SPO2TGT = SPO2TGT,
              KSV = 224.6 * 1.25, KOBLV = 0.075 * 1.4)
  base <- do.call(sim_rop, c(sev, list(end = 200)))
  i <- which(base$TYPE1 > 0.5)
  if (!length(i)) return(tibble(arm = "no type 1 ROP reached"))
  tt <- base$time[i[1]]
  arms <- list(
    untreated  = c(sev, list(end = 300)),
    laser      = c(sev, list(LASERT = tt, end = 300)),
    bev_0.625  = c(sev, list(dose_mg = 0.625, drug = "bevacizumab", tdose = tt, end = 300)),
    bev_0.031  = c(sev, list(dose_mg = 0.031, drug = "bevacizumab", tdose = tt, end = 300)),
    ranib_0.2  = c(sev, list(dose_mg = 0.200, drug = "ranibizumab", tdose = tt, end = 300)),
    aflib_0.4  = c(sev, list(dose_mg = 0.400, drug = "aflibercept", tdose = tt, end = 300))
  )
  bind_rows(lapply(names(arms), function(nm) {
    d <- do.call(sim_rop, arms[[nm]])
    post <- d$time > tt + 7
    ri <- which(d$TYPE1 > 0.5 & post)
    tibble(arm = nm, treat_PMA = base$PMAo[i[1]],
           reactivation      = as.numeric(length(ri) > 0),
           react_PMA         = if (length(ri)) d$PMAo[ri[1]] else NA_real_,
           react_wk_after    = if (length(ri)) (d$PMAo[ri[1]] - base$PMAo[i[1]]) else NA_real_,
           wk_VEGF_below_KG  = sum(d$VEGFFo < 300 & d$time > tt) * 0.5 / 7,
           residual_avascular = tail(d$AVASCo, 1),
           VASC_final        = tail(d$VASC, 1),
           myopia_D          = tail(d$MYOP, 1),
           field_loss        = tail(d$VFLOSS, 1),
           max_detach        = max(d$DETACH))
  }))
}

# S5 — transfusion shifts P50 at a fixed SpO2 target -----------------
scenario_transfusion <- function() {
  arms <- list(none = list(),
               one_d21 = list(TXT1 = 21),
               two = list(TXT1 = 18, TXT2 = 32),
               three_heavy = list(TXT1 = 18, TXT2 = 32, FTX = 0.33))
  bind_rows(lapply(names(arms), function(nm) {
    d <- do.call(sim_rop, c(arms[[nm]], list(SPO2TGT = 93.5)))
    tibble(arm = nm,
           P50_35wk   = approx(d$PMAo, d$P50o, 35)$y,
           FHbF_35wk  = approx(d$PMAo, d$FHBF, 35)$y,
           PaO2_30wk  = approx(d$PMAo, d$PAO2o, 30)$y,
           PaO2_35wk  = approx(d$PMAo, d$PAO2o, 35)$y,
           AVASC_switch = d$AVASCo[which.max(d$VEGFFo)],
           peak_VEGF  = max(d$VEGFFo),
           peak_NV    = max(d$NV))
  }))
}

# S6 — nutrition / IGF-1 (and the rhIGF-1 replacement question) ------
scenario_igf1 <- function() {
  arms <- list(adequate = list(NUTR = 1.0),
               growth_failure = list(NUTR = 0.55),
               rhIGF1_replacement = list(KSI = 3.2 * 2.4),
               replacement_in_growth_failure = list(NUTR = 0.55, KSI = 3.2 * 2.4))
  bind_rows(lapply(names(arms), function(nm) {
    d <- do.call(sim_rop, c(arms[[nm]], list(SPO2TGT = 95)))
    tibble(arm = nm, IGF1_mean = mean(d$IGF1),
           VASC_33wk = approx(d$PMAo, d$VASC, 33)$y,
           AVASC_switch = d$AVASCo[which.max(d$VEGFFo)],
           peak_NV = max(d$NV), type1 = max(d$TYPE1))
  }))
}

# S7 — the Krogh diffusion ceiling -----------------------------------
scenario_diffusion_ceiling <- function(pma = c(28, 30, 32, 33, 34, 35, 36, 38, 40)) {
  p <- as.list(param(mod_rop))
  sg <- function(x) 1 / (1 + exp(-x))
  bind_rows(lapply(pma, function(a) {
    mo2 <- p$MO2MIN + p$MO2MAX * sg((a - p$PMA50M) / p$TAUM)
    h   <- p$HRET0 + (p$HRETMX - p$HRET0) * sg((a - p$PMA50H) / p$TAUH)
    pc  <- mo2 * h^2 / (2 * p$KDIFF)
    s4p <- function(pao2, p50) { r <- (pao2 / p50)^p$NHILL; 100 * r / (1 + r) }
    tibble(PMA = a, MO2 = mo2, h_um = h * 1e4, PCRIT_mmHg = pc,
           PaO2_required = pc / p$CHOR,
           SpO2_required_HbF = s4p(pc / p$CHOR, p$P50F),
           SpO2_required_HbA = s4p(pc / p$CHOR, p$P50A),
           Lpen_at_PaO2_60_um  = sqrt(2 * p$KDIFF * 0.85 * 60 / mo2) * 1e4,
           Lpen_at_PaO2_110_um = sqrt(2 * p$KDIFF * 0.85 * 110 / mo2) * 1e4)
  }))
}

# S8 — the oximetry offset (arithmetic, no ODE needed) ---------------
scenario_oximetry <- function(spo2 = c(85, 87, 89, 91, 93, 95, 96, 97, 98, 99)) {
  p <- as.list(param(mod_rop))
  f <- function(s, p50) p50 * (s / 100 / (1 - s / 100))^(1 / p$NHILL)
  tibble(SpO2 = spo2, PaO2_HbF = f(spo2, p$P50F), PaO2_HbA = f(spo2, p$P50A),
         ratio = f(spo2, p$P50A) / f(spo2, p$P50F),
         dPaO2_per_point = f(spo2 + 0.5, p$P50F) - f(spo2 - 0.5, p$P50F))
}

# S9 — virtual cohort under two SpO2 targets -------------------------
#      Reproduces NeOProM at an ACHIEVED separation of 2.0 points.
cohort_rop <- function(n = 200, seed = 7, targets = c(91.5, 93.5), end = 140) {
  set.seed(seed)
  ga <- pmin(pmax(rnorm(n, 26.0, 1.1), 23), 27.9)
  idata <- tibble(
    ID = seq_len(n), GA = ga,
    BW = pmin(pmax(832 * exp(0.16 * (ga - 26)) * rlnorm(n, 0, 0.14), 380), 1500),
    NUTR   = rlnorm(n, 0, 0.22),
    KGROW  = 0.070 * rlnorm(n, 0, 0.20),
    KOBLV  = 0.075 * rlnorm(n, 0, 0.30),
    KNVR   = 0.850 * rlnorm(n, 0, 0.25),
    KSV    = 224.6 * rlnorm(n, 0, 0.20),
    FHBF0  = pmin(pmax(rnorm(n, 0.90, 0.04), 0.6), 0.97),
    MO2MAX = 3.20e-4 * rlnorm(n, 0, 0.12),
    PMA50M = rnorm(n, 32.5, 0.8),
    TXT1   = ifelse(runif(n) < 0.60, runif(n, 12, 30), 1e6),
    TXT2   = ifelse(runif(n) < 0.35, runif(n, 32, 55), 1e6),
    SEPT   = ifelse(runif(n) < 0.30, runif(n, 10, 45), 1e6))
  bind_rows(lapply(targets, function(tg) {
    out <- mrgsim_df(update(mod_rop, param = list(SPO2TGT = tg)), idata = idata,
                     end = end, delta = 1, atol = 1e-9, rtol = 1e-7)
    out %>% group_by(ID) %>%
      summarise(type1 = max(TYPE1), plus = max(PLUSD), stage = max(STAGE),
                pdeath = max(PDEATH), .groups = "drop") %>%
      summarise(SpO2_target = tg, n = n(),
                treated_ROP_pct = 100 * mean(type1),
                plus_pct = 100 * mean(plus),
                stage3plus_pct = 100 * mean(stage >= 3),
                death_pct = 100 * mean(pdeath))
  }))
}

# ---------------------------------------------------------------------
#  run everything
# ---------------------------------------------------------------------
run_all_rop <- function(cohort_n = 200) {
  cat("\n=== S8  oximetry: the same SpO2 is two different PaO2 ===\n")
  print(as.data.frame(scenario_oximetry()), digits = 4)
  cat("\n=== S7  Krogh ceiling: choroidal pO2 needed to reach the inner retina ===\n")
  print(as.data.frame(scenario_diffusion_ceiling()), digits = 4)
  cat("\n=== S1  oxygen target sweep: the phase-1 dose-response ===\n")
  print(as.data.frame(scenario_oxygen_sweep()), digits = 4)
  cat("\n=== S2  supplemental oxygen at prethreshold (STOP-ROP design) ===\n")
  print(as.data.frame(scenario_supplemental_o2()), digits = 4)
  cat("\n=== S3  bevacizumab dose ladder (PEDIG de-escalation) ===\n")
  print(as.data.frame(scenario_dose_ladder()), digits = 4)
  cat("\n=== S3b ranibizumab 0.2 mg for comparison ===\n")
  print(as.data.frame(scenario_dose_ladder(doses = c(0.2, 0.1), drug = "ranibizumab")),
        digits = 4)
  cat("\n=== S4  laser versus anti-VEGF ===\n")
  print(as.data.frame(scenario_treatment()), digits = 4)
  cat("\n=== S5  red-cell transfusion at a fixed SpO2 target ===\n")
  print(as.data.frame(scenario_transfusion()), digits = 4)
  cat("\n=== S6  nutrition, IGF-1 and rhIGF-1 replacement ===\n")
  print(as.data.frame(scenario_igf1()), digits = 4)
  cat(sprintf("\n=== S9  virtual cohort, n=%d, achieved separation 2.0 points ===\n", cohort_n))
  print(as.data.frame(cohort_rop(n = cohort_n)), digits = 4)
  cat("\nNeOProM observed: treated ROP 14.9%% (higher target) vs 10.9%% (lower),",
      "\n                  RR 0.74; death 17.1%% vs 19.9%%, RR 1.17.\n")
  invisible(NULL)
}

if (identical(environment(), globalenv()) && !interactive()) {
  # run_all_rop()
}
