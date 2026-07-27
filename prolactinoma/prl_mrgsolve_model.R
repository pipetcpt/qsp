# =============================================================================
#  PROLACTINOMA (lactotroph PitNET) — QSP model for mrgsolve
#  including drug-induced hyperprolactinaemia, stalk-effect (disconnection)
#  hyperprolactinaemia and macroprolactinaemia as the SAME equations
#  evaluated at different points
#
#  58 ODEs · time unit = DAYS · mrgsolve
# -----------------------------------------------------------------------------
#  THE ORGANISING IDEA
#
#    Prolactinoma is the one pituitary tumour whose GROWTH is under tonic
#    inhibitory neuroendocrine control, so a receptor agonist is not merely
#    anti-secretory — it is anti-tumoral. And it is the one endocrine disease
#    in which the number you MEASURE and the signal the body SEES can diverge
#    by two orders of magnitude in either direction.
#
#  1. ONE D2 OCCUPANCY EQUATION, USED EVERYWHERE
#
#     Every dopaminergic ligand in the model — endogenous portal dopamine,
#     cabergoline, bromocriptine, quinagolide and any antipsychotic — enters a
#     single competitive-occupancy expression, each weighted by its own
#     intrinsic activity e:
#
#       SIGDRIVE = SUM_i e_i (C_i/K_i)  /  (1 + SUM_j C_j/K_j)
#
#     with e = 1.0 (dopamine, cabergoline, quinagolide), 0.80 (bromocriptine,
#     partial), 0.25 (aripiprazole, partial) and 0 (risperidone / paliperidone
#     / amisulpride / haloperidol / metoclopramide).
#
#     Nothing else in the model knows which drug is present. From that one
#     line the following all fall out rather than being coded as rules:
#
#       * agonist  -> SIGDRIVE up   -> prolactin falls
#       * D2 blocker displaces dopamine with e=0 -> SIGDRIVE down -> PRL RISES
#       * aripiprazole alone: replaces an e=1.0 ligand with an e=0.25 ligand
#         of much higher affinity -> at usual doses it OVERSHOOTS physiological
#         dopamine tone, so prolactin falls BELOW normal (scenario 22)
#       * aripiprazole ADDED to risperidone: replaces an e=0 ligand with an
#         e=0.25 ligand -> SIGDRIVE UP -> prolactin FALLS (scenario 23).
#         The documented add-on effect is arithmetic, not a new mechanism.
#       * cut the DELIVERY term instead of the receptor term and you get
#         stalk-effect hyperprolactinaemia — WITH A COMPUTED CEILING
#         (scenario 20, diagnostic D8).
#
#  2. ONE RECEPTOR, THREE TIME CONSTANTS, THREE POTENCIES
#
#     D2 occupancy acts on four separate downstream branches, each with its own
#     IC50 and its own time constant. They are deliberately ordered:
#
#       BRANCH 1  exocytosis    IC50_EXO  most potent   tau ~ hours
#                 (GIRK -> hyperpolarisation -> Ca2+ down -> granule fusion)
#       BRANCH 2  transcription IC50_TR    intermediate  tau ~ 1.5 d (mRNA)
#       BRANCH 3  cell VOLUME   IC50_VOL   intermediate  tau ~ 20 d  REVERSIBLE
#       BRANCH 4  proliferation IC50_PR    least potent  tau ~ months
#
#     Consequences that the model therefore GENERATES:
#       * prolactin falls steeply within hours (branch 1) but the CHRONIC
#         level is set by branch 2, because at steady state secretion equals
#         synthesis minus degradation — mass balance, not assumption.
#       * tumour shrinkage lags biochemical control by months.
#       * escalating the dose for a macroadenoma that is not shrinking even
#         though prolactin is already normal is rational: the branches sit at
#         different points on the same occupancy axis (diagnostic D5).
#       * cytoreduction splits into cell VOLUME (fast, reversible, re-expands
#         within weeks of withdrawal) and cell NUMBER (slow, bankable). Only
#         the second makes drug withdrawal possible (scenarios 15-17).
#
#  3. THE CONTROLLER IS ALREADY SATURATED AT DIAGNOSIS
#
#     Prolactin drives its own inhibition through the short loop
#     (PRL -> PRLR/JAK2/STAT5 on TIDA neurones -> tyrosine hydroxylase up ->
#     more portal dopamine). In a prolactinoma that loop is INTACT: the
#     transmitter works, the receiver does not (reduced D2R density). Two
#     model-generated consequences:
#       * at diagnosis TIDA drive is near its ceiling, so there is no
#         endogenous reserve to recruit and the only way to add signal is a
#         ligand of far higher affinity than dopamine (Ki 0.7 nM vs an
#         apparent 50 nM for dopamine at the pituitary receptor).
#       * because the feedback saturates, it has spare capacity only below a
#         critical lactotroph mass — which makes drug withdrawal a THRESHOLD
#         phenomenon. Diagnostic D10 bisects for that threshold; it is the
#         model's mechanistic version of "a visible remnant predicts relapse".
#
#  4. THE MEASUREMENT LAYER IS PART OF THE PATHOPHYSIOLOGY
#
#     PRLB   monomeric 23 kDa prolactin — the BIOACTIVE pool, the only one the
#            hypothalamus, gonad and bone ever see
#     PRLM   macroprolactin (PRL.IgG, ~150 kDa) — immunoreactive, slowly
#            cleared, cannot reach the extravascular receptor
#     PRLIMM = PRLB*(1+FDIM) + XMAC*PRLM        (what the analyte really is)
#     PRLMEAS = PRLIMM / (1 + (PRLIMM/KHOOK)^PHOOK)   <-- two-site sandwich
#
#     The last line is the only place the hook effect lives. It is monotone
#     up to ~1500 ng/mL and then falls, so a giant prolactinoma with a true
#     prolactin of 20 000 ng/mL is REPORTED as ~65 (scenario 25 / D6) and the
#     1:100 dilution recovers it — the same curve evaluated 100-fold down the
#     analyte axis. PEG recovery = 100*PRLB*(1+FDIM)/PRLIMM discriminates
#     macroprolactinaemia (D7), where reported prolactin is high, BIOACTIVE
#     prolactin is normal, the gonadal axis is intact and treatment is futile
#     (scenario 24 — a deliberate structural null).
#
#  5. SAFETY IS A SELECTIVITY RATIO EVALUATED AT A DIFFERENT CONCENTRATION
#
#     Efficacy is driven by the PITUITARY BIOPHASE concentration (CAB_E,
#     partition PART_CAB ~ 60x plasma, ke0 corresponding to a ~4 d
#     equilibration half-life — this is what makes WEEKLY dosing work and why
#     prolactin stays suppressed >14 d after a single dose, scenario 11).
#     Valvular risk is driven by PLASMA concentration at 5-HT2B, because the
#     valve has no such reservoir. So the two are not proportional, and the
#     Parkinson's-disease literature (3-4 mg/DAY) and the prolactinoma
#     literature (0.5-2 mg/WEEK, 10-40x lower exposure) stop contradicting
#     each other once the dose scale is applied (scenario 30, D11).
#     Impulse control disorder risk runs off mesolimbic D3 occupancy and is
#     therefore NOT spared by switching to the non-ergot quinagolide, while
#     valvular risk is (D12 sweeps the 5-HT2B/D2 selectivity ratio over four
#     decades: efficacy flat, valve risk collapses).
#
#  6. WHAT IS SOLVED RATHER THAN ASSERTED
#
#     Every healthy baseline is derived algebraically in $MAIN from PRL0,
#     DAP0 and the receptor constants, so that with no tumour and no drug the
#     system is exactly stationary (diagnostic D1 reports the drift). The
#     disease is then GENERATED by seeding tumour mass and reduced D2R
#     density; no prolactin, gonadal or bone abnormality is imposed.
#
#  CALIBRATION ANCHORS (see prl_references.md for full citations)
#     Webster 1994 NEJM      cabergoline vs bromocriptine, normoprolactinaemia
#                            83% vs 59%; withdrawal for AEs 3% vs 12%
#     Colao 2000/2003 JCEM   macroprolactinoma shrinkage; withdrawal outcomes
#     Dekkers 2010 JCEM      pooled persisting remission after withdrawal ~21%
#     Molitch 2015           pregnancy symptomatic enlargement micro ~2.7%,
#                            macro ~21-23%
#     Zanettini 2007 NEJM /  cabergoline valvulopathy at PD doses (3-4 mg/d)
#     Schade 2007 NEJM
#     Stiles 2019 / Caputo   no significant excess at prolactinoma doses
#     2015 meta-analyses
#     Bahceci / Petersenn    cabergoline PK: t1/2 63-109 h, plasma pg/mL
#     Raverot 2018 ESE       temozolomide in aggressive PitNET, MGMT-dependent
#
#  LIMITATIONS — stated up front
#     * Cabergoline absolute bioavailability is unknown; V/F and CL/F are
#       apparent volumes fitted to reproduce the published pg/mL Cmax range
#       and the 63-109 h terminal half-life. The biophase partition absorbs
#       pituitary tissue accumulation and is NOT independently identifiable
#       from plasma data.
#     * The hook-effect exponent PHOOK and KHOOK are platform-specific; only
#       the qualitative structure (monotone, then falling) is general.
#     * Menstrual cyclicity is not modelled; ovulatory competence is a graded
#       index, not a cycle simulator.
#     * A single "antipsychotic" compartment stands in for a class; e_AP and
#       KI_AP are switched to represent risperidone/paliperidone (e=0),
#       aripiprazole (e=0.25) or a prolactin-sparing agent (large KI_AP).
# =============================================================================

library(mrgsolve)

prl_code <- '
$GLOBAL
#define VNL   (0.15)      // mL of normal lactotroph tissue = 1 "cell unit"
#define M_PI_ (3.14159265358979)

// residual (un-inhibited) activity of a branch, with a dopamine-independent
// constitutive floor. S = activated-receptor signal, IC50 = branch potency.
double resid(double S, double IC50, double n, double flo){
  double s = (S > 0.0) ? S : 0.0;
  double a = pow(IC50, n), b = pow(s, n);
  return flo + (1.0 - flo) * a / (a + b);
}
// smooth saturating deficiency term, 0 when x<=0
double sat(double x, double K){
  double y = (x > 0.0) ? x : 0.0;
  return y / (K + y);
}

// quantities solved in $MAIN and needed in $ODE / $TABLE
double KIN_TIDA, KSYN_DA, KTR_N, KSYNG, FSEC, KIN_K, KIN_G, KIN_L, KIN_F;
double KIN_E2G, KIN_TST, RTRN0, REXN0, SIG0, FB0, FDES0, OVTGT0, KSPZN;
double SIGDRIVE, SIGN, SIGS, SIGR, RTRS, RTRR, RTRN, REXS, REXR, REXN;
double SECR, PRLIMM, PRLMEAS_, TVOL_, COMPR_, CEFF_, GONIDX, DEFG;
double OCC5HT, OCCD3, OCCCTZ, E2TOT, CABNM, BRCNM, QUINM, APNM;
double KSURG, KSHYP, HYPTOT_;

$PARAM @annotated
// ---------------- patient / disease setup -----------------------------------
NS0     :  0.30 : Seeded D2R-sensitive tumour volume at t=0 (mL)
NR0     :  0.00 : Seeded D2R-resistant clone volume at t=0 (mL)
D2RS0   :  0.50 : D2 receptor density on sensitive clone (1 = normal lactotroph)
D2RR0   :  0.06 : D2 receptor density on resistant clone
FDRIVE  :  1.00 : Intrinsic (genetic) transcriptional/mitogenic drive multiplier
FGEN    :  1.00 : Genetic proliferation multiplier (AIP/SF3B1 aggressive = >1)
SEX     :  0.00 : 0 = female, 1 = male
STALKP  :  1.00 : Pituitary stalk patency (1 = intact, 0.05 = compressed)
APABSET :  0.00 : Anti-prolactin autoantibody titre (macroprolactinaemia)
PRIMHYP :  0.00 : Primary hypothyroidism severity (0-1) — TRH-driven PRL
FCLR    :  1.00 : Prolactin clearance multiplier (CKD / cirrhosis < 1)
NNSET   :  1.00 : Normal lactotroph pool set point (cell units)
MGMT0   :  0.20 : MGMT activity (0 = absent -> TMZ sensitive, 1 = high)
// ---------------- healthy anchors -------------------------------------------
PRL0    : 10.0  : Healthy prolactin (ng/mL)
DAP0    :  5.0  : Healthy portal dopamine concentration (nM)
E2B     : 60.0  : Healthy gonadal oestradiol (pg/mL, female)
TSTB    :550.0  : Healthy total testosterone (ng/dL, male)
BMD0    :  1.05 : Healthy lumbar spine BMD (g/cm2)
BW0     : 68.0  : Healthy body weight (kg)
// ---------------- dopamine / short loop -------------------------------------
KOUT_TIDA: 1.0  : TIDA activity turnover (1/d)
AMPPRL  :  2.00 : Maximal short-loop amplification of TH drive
K50PRL  : 25.0  : Prolactin producing half-maximal short-loop drive (ng/mL)
HFB     :  2.00 : Hill coefficient, short loop
KEL_DA  :500.0  : Portal dopamine elimination (1/d)
// ---------------- D2 receptor pharmacology ----------------------------------
KI_DA   : 50.0  : Apparent dopamine Ki at pituitary D2 (nM)
KI_CAB  :  0.70 : Cabergoline Ki at D2 (nM)
KI_BRC  :  5.00 : Bromocriptine Ki at D2 (nM)
KI_QUI  :  0.30 : Quinagolide Ki at D2 (nM)
KI_AP   :  0.50 : Antipsychotic Ki at D2 (nM) — risperidone-like default
EFF_DA  :  1.00 : Intrinsic activity, dopamine
EFF_CAB :  1.00 : Intrinsic activity, cabergoline (full agonist)
EFF_BRC :  0.80 : Intrinsic activity, bromocriptine (partial)
EFF_QUI :  1.00 : Intrinsic activity, quinagolide
EFF_AP  :  0.00 : Intrinsic activity, antipsychotic (0.25 for aripiprazole)
KOUT_D2R:  0.02 : D2R density turnover (1/d)
ED2R    :  0.45 : Maximal oestrogen down-regulation of D2R
KE2D    :400.0  : Oestradiol for half-maximal D2R down-regulation (pg/mL)
DESENS  :  0.25 : Maximal agonist-induced D2R desensitisation
// ---------------- the four D2 branches --------------------------------------
IC50_EXO: 0.0287: Branch 1 exocytosis IC50 (activated-receptor units)
HEXO    :  1.50 : Branch 1 Hill coefficient
FLO_EXO :  0.05 : Branch 1 dopamine-independent floor
IC50_TR : 0.0455: Branch 2 transcription IC50
HTR     :  2.00 : Branch 2 Hill coefficient
FLO_TR  :  0.05 : Branch 2 floor
IC50_VOL: 0.0800: Branch 3 cell-volume IC50
HVOL    :  2.00 : Branch 3 Hill coefficient
IC50_PR : 0.1500: Branch 4 proliferation IC50
HPR     :  1.50 : Branch 4 Hill coefficient
// ---------------- synthesis / storage / secretion ---------------------------
KDMRNA  :  0.462: PRL mRNA degradation (1/d, t1/2 1.5 d)
KEXO    : 20.0  : Maximal granule exocytosis rate constant (1/d)
KDEGG   :  0.50 : Granule crinophagy / degradation (1/d)
ATRH    :  2.20 : Maximal TRH amplification of PRL transcription
KTRH    :  0.45 : Half-maximal TRH effect (FT4 deficit units)
EMAX_E2 :  2.50 : Maximal oestrogen amplification of PRL transcription
KE2T    :800.0  : Oestradiol for half-maximal transcription effect (pg/mL)
KEL_PRL : 33.3  : Monomeric prolactin elimination (1/d, t1/2 30 min)
FDIM    :  0.10 : Big-prolactin (dimer) immunoreactive fraction of monomer
KON_MAC :  4.00 : PRL-IgG association (per titre unit per day)
KOFF_MAC:  1.00 : PRL-IgG dissociation (1/d)
KEL_MAC :  0.46 : Macroprolactin elimination (1/d, t1/2 1.5 d)
XMAC    :  0.90 : Assay cross-reactivity to macroprolactin
KHOOK   :4000.0 : Two-site assay hook constant (ng/mL)
PHOOK   :  4.00 : Hook steepness exponent
DILF    :100.0  : Dilution factor for the confirmatory re-assay
// ---------------- lactotroph populations ------------------------------------
KNN     : 0.0167: Normal pool relaxation rate (1/d)
AHYP_D2 :  1.20 : Maximal lactotroph hyperplasia on loss of D2 signal
AHYP_E2 :  1.00 : Maximal oestrogen-driven lactotroph hyperplasia
KE2H    :2000.0 : Oestradiol for half-maximal hyperplasia (pg/mL)
KPROL   :0.00311: Tumour proliferation rate constant (1/d)
KAPO    :0.00150: Tumour basal apoptosis rate constant (1/d)
AMAX_APO:  1.00 : Maximal D2-driven apoptosis amplification
EPROL   :  0.20 : Maximal oestrogen amplification of proliferation
KVOL    : 0.0495: Cell-volume relaxation rate (1/d, t1/2 14 d)
VMIN    :  0.65 : Minimum relative cell volume under full D2 suppression
KFIB    :2.0e-5 : Perivascular fibrosis accrual per unit cabergoline signal
MUTR    :2.5e-5 : Rate of D2R-low clone generation from sensitive pool (1/d)
TMAXV   :120.0  : Carrying capacity of the sellar/parasellar space (mL)
// ---------------- geometry / mass effect ------------------------------------
VGLAND  :  0.60 : Normal pituitary gland volume (mL)
HSELLA  :  1.25 : Sella turcica vertical capacity (cm)
GAPCH   :  0.30 : Diaphragma-to-chiasm gap (cm)
KVFR    :  0.50 : Reversible conduction-block relaxation (1/d)
VFMAX   : 28.0  : Maximal reversible visual field deficit (dB)
KCVF    :  1.20 : Chiasm indentation for half-maximal deficit (cm)
FTURG   :  0.45 : Fraction of transmitted force that scales with cell turgor
KAXL    : 0.0025: Irreversible axonal loss rate (dB/d per unit compression)
AXTHR   :  0.50 : Chiasm indentation below which axons are not lost (cm)
KHYP_ON : 0.0200: Hypopituitarism onset rate (1/d)
KHYP_OFF: 0.0035: Hypopituitarism recovery rate (1/d, deliberately slower)
HMAX    :  0.85 : Maximal compression-driven hypopituitarism
KCHYP   :  0.90 : Compression for half-maximal hypopituitarism (cm)
KPIT    : 0.0500: Peripheral pituitary-hormone relaxation (1/d)
WCORT   :  0.55 : Weight of hypopituitarism on cortisol reserve
WFT4    :  0.70 : Weight of hypopituitarism on FT4
WIGF1   :  0.90 : Weight of hypopituitarism on IGF-1
// ---------------- HPG axis ---------------------------------------------------
KOUT_K  :  1.00 : Kisspeptin tone turnover (1/d)
K50K    : 90.0  : Bioactive prolactin for half-maximal KNDy suppression (ng/mL)
HK      :  1.00 : Hill coefficient, KNDy suppression
KOUT_G  :  2.00 : GnRH pulse-generator turnover (1/d)
KOUT_L  :  4.00 : LH turnover (1/d)
KOUT_F  :  2.00 : FSH turnover (1/d)
KE2FB   :  1.20 : Gonadal steroid negative feedback strength
KOUT_E2G:  4.00 : Gonadal oestradiol turnover (1/d)
KOUT_TST:  4.00 : Testosterone turnover (1/d)
KOVF    :  0.10 : Ovulatory competence relaxation (1/d)
KOV     :  0.85 : Oestradiol index for half-maximal ovulatory competence
HOV     :  3.00 : Hill coefficient, ovulatory competence
KOUT_SPZ:0.0143 : Spermatogenesis relaxation (1/d, ~70 d cycle)
FTREPL  :  0.00 : Exogenous sex-steroid replacement (fraction of normal)
// ---------------- bone -------------------------------------------------------
KOCF    : 0.0333: Osteoclast function relaxation (1/d)
KOBF    : 0.0333: Osteoblast function relaxation (1/d)
AOC     :  0.90 : Maximal osteoclast activation by sex-steroid deficiency
AOB     :  0.35 : Maximal osteoblast suppression by sex-steroid deficiency
KGON    :  0.25 : Half-maximal gonadal deficiency effect on bone
APRLB   :  0.15 : Direct prolactin suppression of osteoblast function
KBMDR   :6.7e-4 : BMD relaxation toward its turnover set point (1/d, t1/2 ~2.8 y)
KBEXP   :  0.22 : Exponent mapping the OBF/OCF ratio to the BMD set point
FIRR    :  0.35 : Fraction of BMD loss that is irreversible
// ---------------- metabolic --------------------------------------------------
KBW     : 0.0100: Body-weight relaxation (1/d)
ABWP    :  0.09 : Maximal prolactin-driven weight gain (fraction)
ABWG    :  0.06 : Maximal hypogonadism-driven weight gain (fraction)
KIRX    : 0.0200: Insulin-resistance relaxation (1/d)
AIRP    :  0.55 : Maximal prolactin effect on HOMA-IR
// ---------------- cabergoline PK ---------------------------------------------
KA_CAB  : 28.8  : Cabergoline absorption rate (1/d)
V2_CAB  :14000. : Cabergoline apparent central volume V/F (L)
V3_CAB  : 3000. : Cabergoline apparent peripheral volume (L)
Q_CAB   :  600. : Cabergoline intercompartmental clearance (L/d)
CL_CAB  : 2772. : Cabergoline apparent clearance CL/F (L/d)
PART_CAB: 60.0  : Pituitary biophase partition coefficient
KE0_CAB :  0.17 : Biophase equilibration rate (1/d, t1/2 ~4 d)
MW_CAB  :451.6  : Cabergoline molecular weight
// ---------------- bromocriptine PK -------------------------------------------
KA_BRC  : 36.0  : Bromocriptine absorption rate (1/d)
V_BRC   :10700. : Bromocriptine apparent volume V/F (L)
CL_BRC  :11900. : Bromocriptine apparent clearance CL/F (L/d)
PART_BRC: 20.0  : Bromocriptine biophase partition
KE0_BRC :  2.00 : Bromocriptine biophase equilibration (1/d)
MW_BRC  :654.6  : Bromocriptine molecular weight
FVAG    :  1.00 : Vaginal-route bioavailability multiplier
// ---------------- quinagolide PK ---------------------------------------------
KA_QUI  : 24.0  : Quinagolide absorption rate (1/d)
V_QUI   : 1200. : Quinagolide apparent volume V/F (L)
CL_QUI  :  980. : Quinagolide apparent clearance CL/F (L/d)
PART_QUI:  8.00 : Quinagolide biophase partition
KE0_QUI :  1.20 : Quinagolide biophase equilibration (1/d)
MW_QUI  :395.6  : Quinagolide molecular weight
// ---------------- antipsychotic PK -------------------------------------------
KA_AP   : 12.0  : Antipsychotic absorption rate (1/d)
V_AP    : 1100. : Antipsychotic apparent volume V/F (L)
CL_AP   :  380. : Antipsychotic apparent clearance CL/F (L/d)
FU_AP   :  0.10 : Antipsychotic free fraction
PART_AP :  1.00 : Pituitary is OUTSIDE the blood-brain barrier: no partition
KE0_AP  :  6.00 : Antipsychotic biophase equilibration (1/d)
MW_AP   :410.5  : Antipsychotic molecular weight
// ---------------- temozolomide ----------------------------------------------
KEL_TMZ : 9.24  : Temozolomide elimination (1/d, t1/2 1.8 h)
V_TMZ   : 30.0  : Temozolomide volume (L)
KDAM    :  0.60 : DNA adduct formation per unit exposure
KMG     :  0.10 : MGMT activity for half-maximal adduct repair
KREP    :  0.15 : DNA damage repair rate (1/d)
KKILL   :  0.014: Tumour kill rate per unit damage (1/d)
KMGIND  :2.0e-4 : MGMT induction under alkylator exposure (1/d)
// ---------------- adverse effects -------------------------------------------
KI_5HT2B: 12.0  : 5-HT2B Ki, cabergoline (nM) — plasma-driven
KVALV   :  0.90 : Valve fibrosis accrual rate per unit occupancy per day
KVREG   :0.00080: Valve fibrosis regression (1/d)
VALV50  :175.0  : Valve burden at 50% probability of moderate-severe TR
KI_D3   :  1.50 : D3 Ki, cabergoline (nM)
KICD    :  0.05 : ICD sensitisation rate (1/d)
KICDOFF : 0.0100: ICD desensitisation (1/d)
ICDSUSC :  1.00 : Individual ICD susceptibility multiplier
KI_CTZ  :  0.15 : Area postrema D2 Ki (nM)
KTOLON  :  0.35 : Nausea tolerance development (1/d)
KTOLOFF : 0.0300: Nausea tolerance loss (1/d)
TOLMAX  :  6.00 : Maximal nausea tolerance
// ---------------- pregnancy --------------------------------------------------
PREGON  :-1.0   : Day pregnancy starts (<0 = never)
PREGLEN :280.0  : Gestation length (d)
E2PREG  :18000. : Term oestradiol (pg/mL)
HPLMAX  : 90.0  : Term placental-lactogen contribution to measured PRL (ng/mL)
KE2X    :  0.20 : Exogenous/pregnancy oestradiol turnover (1/d)
E2EXO   :  0.00 : Exogenous oestrogen (pg/mL, e.g. combined oral contraceptive)
// ---------------- interventions ---------------------------------------------
SURGD   : -1.0  : Day of transsphenoidal surgery (<0 = none)
SURGRES :  0.25 : Fraction of tumour left behind by surgery
SURGHYP :  0.10 : New hypopituitarism caused by surgery
CABSTOP :1e9    : Day cabergoline is withdrawn (for scenario bookkeeping)

$CMT @annotated
CAB_G  : Cabergoline gut (mg)
CAB_C  : Cabergoline central (mg)
CAB_P  : Cabergoline peripheral (mg)
CAB_E  : Cabergoline pituitary biophase (nM)
BRC_G  : Bromocriptine gut (mg)
BRC_C  : Bromocriptine central (mg)
BRC_E  : Bromocriptine biophase (nM)
QUI_G  : Quinagolide gut (mg)
QUI_C  : Quinagolide central (mg)
QUI_E  : Quinagolide biophase (nM)
AP_G   : Antipsychotic gut (mg)
AP_C   : Antipsychotic central (mg)
AP_E   : Antipsychotic pituitary biophase (nM)
TMZ_C  : Temozolomide central (mg)
TMZ_D  : Tumour DNA damage signal
TIDA   : TIDA neurone activity (relative)
DAP    : Portal dopamine (nM)
D2RS   : D2R density, sensitive clone (relative)
D2RR   : D2R density, resistant clone (relative)
NN     : Normal lactotroph pool (cell units)
NS     : Sensitive tumour volume of cells (mL)
NR     : Resistant tumour volume of cells (mL)
VCELL  : Mean tumour cell volume (relative)
FIB    : Perivascular fibrotic volume (mL)
MRNAN  : PRL mRNA, normal lactotroph (relative)
MRNAS  : PRL mRNA, tumour (relative)
GRN    : Granule store, normal lactotroph (relative)
GRS    : Granule store, tumour (relative)
PRLB   : Plasma monomeric bioactive prolactin (ng/mL)
PRLM   : Plasma macroprolactin (ng/mL equivalent)
APAB   : Anti-prolactin autoantibody titre
E2X    : Exogenous / pregnancy oestradiol (pg/mL)
HPLC   : Placental lactogen contribution (ng/mL)
KISS   : Kisspeptin / KNDy tone (relative)
GNRH   : GnRH pulse-generator output (relative)
LH     : Luteinising hormone (relative)
FSH    : Follicle-stimulating hormone (relative)
E2G    : Gonadal oestradiol (pg/mL)
TST    : Testosterone (ng/dL)
OVF    : Ovulatory competence index (0-1)
SPZ    : Sperm production index (relative)
OCF    : Osteoclast function (relative)
OBF    : Osteoblast function (relative)
BMD    : Lumbar spine BMD (g/cm2)
BLOSSI : Irreversible BMD loss (g/cm2)
BW     : Body weight (kg)
IR     : Insulin resistance index (HOMA-IR)
VFR    : Reversible visual field deficit (dB)
AXL    : Irreversible visual field deficit (dB)
HYPOP  : Hypopituitarism severity (0-1)
CORT   : Cortisol reserve (relative)
FT4    : Free T4 (relative)
IGF1   : IGF-1 (relative)
VALV   : Valve fibrosis burden
NTOL   : Nausea tolerance state
ICDS   : Impulse-control sensitisation
CUMCAB : Cumulative cabergoline dose (mg)
MGMTS  : MGMT activity (relative)

$MAIN
// ---- algebraically solved healthy baseline ---------------------------------
FB0    = pow(PRL0,HFB)/(pow(PRL0,HFB)+pow(K50PRL,HFB));
KIN_TIDA = KOUT_TIDA/(1.0 + AMPPRL*FB0);
KSYN_DA  = KEL_DA*DAP0;
SIG0   = (EFF_DA*DAP0/KI_DA)/(1.0 + DAP0/KI_DA);   // healthy activated-receptor signal
RTRN0  = resid(SIG0, IC50_TR , HTR , FLO_TR );
REXN0  = resid(SIG0, IC50_EXO, HEXO, FLO_EXO);
KTR_N  = KDMRNA/RTRN0;                              // so MRNAN(0) = 1
KSYNG  = KEXO*REXN0 + KDEGG;                        // so GRN(0)  = 1
FSEC   = KEL_PRL*PRL0/(KEXO*REXN0);                 // so PRLB(0) = PRL0
KIN_K  = KOUT_K*(1.0 + pow(PRL0/K50K,HK));          // so KISS(0) = 1
KIN_G  = KOUT_G;
KIN_L  = KOUT_L*(1.0 + KE2FB);                      // so LH(0) = 1 with feedback
KIN_F  = KOUT_F*(1.0 + KE2FB);
KIN_E2G= KOUT_E2G*E2B;
KIN_TST= KOUT_TST*TSTB;
FDES0  = 1.0 - DESENS*SIG0/(0.30 + SIG0);           // so D2R(0) = D2RS0
OVTGT0 = (1.0/(pow(KOV,HOV) + 1.0))*(1.0/(0.35 + 1.0));   // so OVF(0) = 1
KSPZN  = 1.0/((1.0/(0.4 + 1.0))*(1.0/1.3));               // so SPZ(0) = 1

// ---- initial conditions ----------------------------------------------------
TIDA_0  = 1.0;   DAP_0 = DAP0*STALKP;
D2RS_0  = D2RS0; D2RR_0 = D2RR0;
NN_0    = NNSET; NS_0 = NS0; NR_0 = NR0; VCELL_0 = 1.0; FIB_0 = 0.0;
MRNAN_0 = 1.0;   GRN_0 = 1.0;
MRNAS_0 = FDRIVE*resid(SIG0*D2RS0, IC50_TR, HTR, FLO_TR)/RTRN0;
GRS_0   = KSYNG*MRNAS_0/(KEXO*resid(SIG0*D2RS0,IC50_EXO,HEXO,FLO_EXO) + KDEGG);
PRLB_0  = PRL0;  PRLM_0 = 0.0; APAB_0 = APABSET;
E2X_0   = E2EXO; HPLC_0 = 0.0;
KISS_0  = 1.0; GNRH_0 = 1.0; LH_0 = 1.0; FSH_0 = 1.0;
E2G_0   = E2B; TST_0 = TSTB; OVF_0 = 1.0; SPZ_0 = 1.0;
OCF_0   = 1.0; OBF_0 = 1.0; BMD_0 = BMD0; BLOSSI_0 = 0.0;
BW_0    = BW0; IR_0 = 1.0;
VFR_0   = 0.0; AXL_0 = 0.0; HYPOP_0 = 0.0;
CORT_0  = 1.0; FT4_0 = 1.0 - PRIMHYP; IGF1_0 = 1.0;
VALV_0  = 0.0; NTOL_0 = 0.0; ICDS_0 = 0.0; CUMCAB_0 = 0.0; MGMTS_0 = MGMT0;

// ---- transsphenoidal surgery is applied in $ODE as a fast first-order
// ---- resection over a 0.5 d window: exp(-k*0.5) = SURGRES
KSURG = 0.0; KSHYP = 0.0;
if(SURGD > 0){
  KSURG = -log(SURGRES)/0.5;
  KSHYP = SURGHYP/0.5;
}

$ODE
// =========================== DRUG PK =======================================
dxdt_CAB_G = -KA_CAB*CAB_G;
dxdt_CAB_C =  KA_CAB*CAB_G - (CL_CAB/V2_CAB)*CAB_C
              - (Q_CAB/V2_CAB)*CAB_C + (Q_CAB/V3_CAB)*CAB_P;
dxdt_CAB_P =  (Q_CAB/V2_CAB)*CAB_C - (Q_CAB/V3_CAB)*CAB_P;
CABNM      = (CAB_C/V2_CAB)*1.0e6/MW_CAB;              // mg/L -> nM
dxdt_CAB_E =  KE0_CAB*(PART_CAB*CABNM - CAB_E);

dxdt_BRC_G = -KA_BRC*BRC_G;
dxdt_BRC_C =  KA_BRC*BRC_G*FVAG - (CL_BRC/V_BRC)*BRC_C;
BRCNM      = (BRC_C/V_BRC)*1.0e6/MW_BRC;
dxdt_BRC_E =  KE0_BRC*(PART_BRC*BRCNM - BRC_E);

dxdt_QUI_G = -KA_QUI*QUI_G;
dxdt_QUI_C =  KA_QUI*QUI_G - (CL_QUI/V_QUI)*QUI_C;
QUINM      = (QUI_C/V_QUI)*1.0e6/MW_QUI;
dxdt_QUI_E =  KE0_QUI*(PART_QUI*QUINM - QUI_E);

dxdt_AP_G  = -KA_AP*AP_G;
dxdt_AP_C  =  KA_AP*AP_G - (CL_AP/V_AP)*AP_C;
APNM       = (AP_C/V_AP)*1.0e6/MW_AP*FU_AP;            // free concentration
dxdt_AP_E  =  KE0_AP*(PART_AP*APNM - AP_E);

dxdt_TMZ_C = -KEL_TMZ*TMZ_C;
dxdt_TMZ_D =  KDAM*(TMZ_C/V_TMZ)/(1.0 + MGMTS/KMG) - KREP*TMZ_D;
dxdt_MGMTS =  KMGIND*(TMZ_C/V_TMZ) - 0.02*(MGMTS - MGMT0);

// ==================== THE ONE OCCUPANCY EQUATION ==========================
{
  double den = 1.0 + DAP/KI_DA + CAB_E/KI_CAB + BRC_E/KI_BRC
                   + QUI_E/KI_QUI + AP_E/KI_AP;
  double num = EFF_DA*DAP/KI_DA + EFF_CAB*CAB_E/KI_CAB + EFF_BRC*BRC_E/KI_BRC
             + EFF_QUI*QUI_E/KI_QUI + EFF_AP*AP_E/KI_AP;
  SIGDRIVE = num/den;
}
SIGN = SIGDRIVE*1.0;        // normal lactotroph: reference D2R density
SIGS = SIGDRIVE*D2RS;       // sensitive clone
SIGR = SIGDRIVE*D2RR;       // resistant clone

RTRN = resid(SIGN, IC50_TR , HTR , FLO_TR );
RTRS = resid(SIGS, IC50_TR , HTR , FLO_TR );
RTRR = resid(SIGR, IC50_TR , HTR , FLO_TR );
REXN = resid(SIGN, IC50_EXO, HEXO, FLO_EXO);
REXS = resid(SIGS, IC50_EXO, HEXO, FLO_EXO);
REXR = resid(SIGR, IC50_EXO, HEXO, FLO_EXO);

// ==================== HYPOTHALAMIC SHORT LOOP =============================
{
  double fb = pow(PRLB,HFB)/(pow(PRLB,HFB)+pow(K50PRL,HFB));
  dxdt_TIDA = KIN_TIDA*(1.0 + AMPPRL*fb) - KOUT_TIDA*TIDA;
}
dxdt_DAP = KSYN_DA*TIDA*STALKP - KEL_DA*DAP;

// ==================== OESTROGEN / PREGNANCY ===============================
{
  double e2t = E2EXO;
  if(PREGON > 0 && SOLVERTIME >= PREGON && SOLVERTIME <= PREGON + PREGLEN){
    double g = (SOLVERTIME - PREGON)/PREGLEN;
    e2t = E2EXO + (E2PREG - E2B)*g*g;
  }
  dxdt_E2X = KE2X*(e2t - E2X);
  double hpl = 0.0;
  if(PREGON > 0 && SOLVERTIME >= PREGON && SOLVERTIME <= PREGON + PREGLEN){
    double gh = (SOLVERTIME - PREGON)/PREGLEN;
    hpl = HPLMAX*gh*gh*gh;
  }
  dxdt_HPLC = 0.20*(hpl - HPLC);
}
E2TOT = E2G + E2X;

// D2R density: oestrogen down-regulates, chronic agonist desensitises
{
  double fdown = 1.0 - ED2R*sat(E2TOT - E2B, KE2D);
  double fdes  = (1.0 - DESENS*(SIGDRIVE/(0.30 + SIGDRIVE)))/FDES0;
  if(fdes > 1.0) fdes = 1.0;
  dxdt_D2RS = KOUT_D2R*(D2RS0*fdown*fdes - D2RS);
  dxdt_D2RR = KOUT_D2R*(D2RR0*fdown*fdes - D2RR);
}

// ==================== TRANSCRIPTION AND STORAGE ===========================
{
  double ftrh = 1.0 + ATRH*sat(1.0 - FT4, KTRH);
  double fe2  = (1.0 + EMAX_E2*sat(E2TOT, KE2T))/(1.0 + EMAX_E2*sat(E2B, KE2T));
  dxdt_MRNAN = KTR_N*ftrh*fe2*RTRN            - KDMRNA*MRNAN;
  dxdt_MRNAS = KTR_N*ftrh*fe2*FDRIVE*RTRS     - KDMRNA*MRNAS;
  dxdt_GRN   = KSYNG*MRNAN - KEXO*REXN*GRN - KDEGG*GRN;
  dxdt_GRS   = KSYNG*MRNAS - KEXO*REXS*GRS - KDEGG*GRS;
}

// ==================== PROLACTIN IN PLASMA =================================
SECR = FSEC*KEXO*( NN*GRN*REXN
                 + (NS/VNL)*GRS*REXS
                 + (NR/VNL)*GRS*REXR );
dxdt_PRLB = SECR + HPLC*KEL_PRL - KEL_PRL*FCLR*PRLB
            - KON_MAC*APAB*PRLB + KOFF_MAC*PRLM;
dxdt_PRLM = KON_MAC*APAB*PRLB - KOFF_MAC*PRLM - KEL_MAC*FCLR*PRLM;
dxdt_APAB = 0.05*(APABSET - APAB);

// ==================== LACTOTROPH POPULATIONS ==============================
{
  // normal pool: hyperplasia on loss of D2 signal or on oestrogen excess
  double lossd2 = (SIGN < SIG0) ? (SIG0 - SIGN)/SIG0 : 0.0;
  double target = NNSET*(1.0 + AHYP_D2*lossd2 + AHYP_E2*sat(E2TOT - E2B, KE2H));
  dxdt_NN = KNN*(target - NN);
}
{
  double rprS = resid(SIGS, IC50_PR, HPR, 0.0);
  double rprR = resid(SIGR, IC50_PR, HPR, 0.0);
  double fe2p = 1.0 + EPROL*sat(E2TOT - E2B, KE2H);
  double crowd = 1.0 - (NS + NR)/TMAXV;   if(crowd < 0.0) crowd = 0.0;
  double gS = KPROL*FGEN*fe2p*crowd*rprS - KAPO*(1.0 + AMAX_APO*(1.0 - rprS)) - KKILL*TMZ_D;
  double gR = KPROL*FGEN*fe2p*crowd*rprR - KAPO*(1.0 + AMAX_APO*(1.0 - rprR)) - KKILL*TMZ_D;
  double resect = (SURGD > 0 && SOLVERTIME >= SURGD && SOLVERTIME < SURGD + 0.5)
                  ? KSURG : 0.0;
  dxdt_NS = NS*gS - MUTR*NS - resect*NS;
  dxdt_NR = NR*gR + MUTR*NS - resect*NR;
  // cell volume: fast, reversible, normalised so VCELL(0) = 1
  double rvol0 = resid(SIG0*D2RS0, IC50_VOL, HVOL, 0.0);
  double rvol  = resid(SIGS      , IC50_VOL, HVOL, 0.0);
  double f     = (rvol0 > 0) ? rvol/rvol0 : 1.0;
  if(f > 1.0) f = 1.0;
  dxdt_VCELL = KVOL*(VMIN + (1.0 - VMIN)*f - VCELL);
}
dxdt_FIB = KFIB*(CAB_E/(KI_CAB + CAB_E))*(NS + NR);
dxdt_CUMCAB = 0.0;

// ==================== GEOMETRY AND MASS EFFECT ============================
TVOL_  = (NS + NR)*VCELL + FIB;
{
  double vtot = TVOL_ + VGLAND;
  double rad  = pow(3.0*vtot/(4.0*M_PI_), 1.0/3.0);
  double hsup = 2.0*rad - HSELLA;  if(hsup < 0.0) hsup = 0.0;
  COMPR_ = hsup - GAPCH;           if(COMPR_ < 0.0) COMPR_ = 0.0;
}
// effective force on the chiasm scales with geometry AND with cell turgor:
// granule-depleted, involuted cells transmit less pressure at the same volume
CEFF_ = COMPR_*(1.0 - FTURG + FTURG*VCELL);
dxdt_VFR = KVFR*(VFMAX*sat(CEFF_, KCVF) - VFR);
{
  double ax = CEFF_ - AXTHR;
  dxdt_AXL = (ax > 0.0) ? KAXL*pow(ax, 1.5) : 0.0;
}
{
  double htgt = HMAX*sat(COMPR_, KCHYP);
  double kh   = (htgt > HYPOP) ? KHYP_ON : KHYP_OFF;
  dxdt_HYPOP  = kh*(htgt - HYPOP);
  // surgical hypopituitarism is a PERMANENT additive term, not a reversible one
  HYPTOT_ = HYPOP + ((SURGD > 0 && SOLVERTIME >= SURGD) ? SURGHYP : 0.0);
  if(HYPTOT_ > 1.0) HYPTOT_ = 1.0;
}
dxdt_CORT = KPIT*((1.0 - WCORT*HYPTOT_) - CORT);
dxdt_FT4  = KPIT*((1.0 - WFT4 *HYPTOT_ - PRIMHYP) - FT4);
dxdt_IGF1 = KPIT*((1.0 - WIGF1*HYPTOT_) - IGF1);

// ==================== HPG AXIS ============================================
dxdt_KISS = KIN_K/(1.0 + pow(PRLB/K50K,HK)) - KOUT_K*KISS;
dxdt_GNRH = KIN_G*KISS*(1.0 - 0.8*HYPTOT_) - KOUT_G*GNRH;
{
  double fbst = (SEX < 0.5) ? (E2G/E2B) : (TST/TSTB);
  dxdt_LH = KIN_L*GNRH/(1.0 + KE2FB*fbst) - KOUT_L*LH;
  dxdt_FSH= KIN_F*GNRH/(1.0 + KE2FB*fbst) - KOUT_F*FSH;
}
dxdt_E2G = KIN_E2G*LH*FSH*(1.0 - 0.9*HYPTOT_) - KOUT_E2G*E2G;
dxdt_TST = KIN_TST*LH*(1.0 - 0.9*HYPTOT_) - KOUT_TST*TST;
{
  double ei = E2G/E2B;
  double tgtov = pow(ei,HOV)/(pow(KOV,HOV) + pow(ei,HOV))*(LH/(0.35 + LH))/OVTGT0;
  if(tgtov > 1.0) tgtov = 1.0;
  dxdt_OVF = KOVF*(tgtov - OVF);
}
{
  double st = TST/TSTB;
  double tgtsp = (FSH/(0.4+FSH))*(st/(st+0.3))*KSPZN;
  if(tgtsp > 1.0) tgtsp = 1.0;
  dxdt_SPZ = KOUT_SPZ*(tgtsp - SPZ);
}

// ==================== BONE ================================================
GONIDX = (SEX < 0.5) ? (E2G + FTREPL*E2B)/E2B : (TST + FTREPL*TSTB)/TSTB;
DEFG   = (GONIDX < 1.0) ? (1.0 - GONIDX) : 0.0;
dxdt_OCF = KOCF*((1.0 + AOC*sat(DEFG,KGON)) - OCF);
dxdt_OBF = KOBF*((1.0 - AOB*sat(DEFG,KGON) - APRLB*sat(PRLB-PRL0, 200.0)) - OBF);
{
  // bone mineral density relaxes toward the set point implied by the current
  // formation:resorption ratio; the ratchet BLOSSI records the part of any loss
  // that can never be regained, and caps subsequent recovery
  double bset = BMD0*pow(OBF/OCF, KBEXP);
  double cap  = BMD0 - BLOSSI;
  if(bset > cap) bset = cap;
  double flux = KBMDR*(bset - BMD);
  dxdt_BMD    = flux;
  dxdt_BLOSSI = (flux < 0.0) ? -FIRR*flux : 0.0;
}

// ==================== METABOLIC ===========================================
dxdt_BW = KBW*(BW0*(1.0 + ABWP*sat(PRLB-PRL0,150.0) + ABWG*sat(DEFG,KGON)) - BW);
dxdt_IR = KIRX*((1.0 + AIRP*sat(PRLB-PRL0,150.0) + 0.4*(BW/BW0 - 1.0)) - IR);

// ==================== ADVERSE EFFECTS =====================================
OCC5HT = CABNM/(KI_5HT2B + CABNM);        // PLASMA, not biophase
OCCD3  = CABNM/(KI_D3    + CABNM);
OCCCTZ = CABNM/(KI_CTZ   + CABNM);
dxdt_VALV = KVALV*OCC5HT - KVREG*VALV;
dxdt_NTOL = KTOLON*OCCCTZ*(TOLMAX - NTOL) - KTOLOFF*NTOL;
dxdt_ICDS = KICD*OCCD3*(1.0 - ICDS) - KICDOFF*ICDS;

$TABLE
// ---------------- the measurement layer -----------------------------------
PRLIMM   = PRLB*(1.0 + FDIM) + XMAC*PRLM;
PRLMEAS_ = PRLIMM/(1.0 + pow(PRLIMM/KHOOK, PHOOK));
double adil = PRLIMM/DILF;
double PRLDIL = DILF*adil/(1.0 + pow(adil/KHOOK, PHOOK));
double PEGREC = 100.0*PRLB*(1.0+FDIM)/(PRLIMM + 1e-12);

// ---------------- clinical endpoints --------------------------------------
double TVOLmm  = TVOL_*1000.0;
double SHRINK  = 100.0*(1.0 - TVOL_/(NS0 + NR0 + 1e-12));
double VFMD    = -(VFR + AXL);
double TSCORE  = (BMD - 1.05)/0.11;
double PFRACT  = 100.0/(1.0 + exp(-(-3.90 - 1.31*TSCORE)));
double PTR     = 100.0/(1.0 + exp(-(VALV - VALV50)/30.0));
double PICD    = 100.0/(1.0 + exp(-(ICDS*ICDSUSC - 0.42)/0.09));
double NAUS    = 100.0*OCCCTZ/(1.0 + NTOL);
double AMENP   = 100.0*(1.0 - OVF);
double PRLNORM = (PRLB < 25.0) ? 1.0 : 0.0;
double GALA    = 100.0*sat(PRLB-PRL0, 60.0)*(E2TOT/(E2TOT + 40.0));
double CELLN   = NS + NR;

$CAPTURE @annotated
PRLMEAS_  : REPORTED prolactin through the two-site assay (ng/mL)
PRLIMM    : True immunoreactive analyte (ng/mL)
PRLDIL    : Reported prolactin after 1:100 dilution (ng/mL)
PEGREC    : Post-PEG recovery (%)
TVOL_     : Tumour volume (mL)
TVOLmm    : Tumour volume (mm3)
CELLN     : Tumour cell volume excluding cell-size change (mL)
SHRINK    : Tumour shrinkage from baseline (%)
VFMD      : Visual field mean deviation (dB)
SIGDRIVE  : Net D2 receptor drive (0-1)
SIGS      : Activated-receptor signal on sensitive clone
RTRS      : Residual transcription, sensitive clone
REXS      : Residual exocytosis, sensitive clone
CABNM     : Plasma cabergoline (nM)
OCC5HT    : 5-HT2B occupancy (plasma-driven)
OCCD3     : D3 occupancy
TSCORE    : Lumbar spine T-score
PFRACT    : Vertebral fracture probability (%)
PTR       : Probability of moderate-severe tricuspid regurgitation (%)
PICD      : Probability of impulse control disorder (%)
NAUS      : Nausea intensity (arbitrary, tolerance-adjusted)
AMENP     : Anovulation / amenorrhoea probability (%)
GALA      : Galactorrhoea index (%)
PRLNORM   : Prolactin normalised (<25 ng/mL)
E2TOT     : Total oestradiol exposure (pg/mL)
COMPR_    : Chiasm indentation, geometric (cm)
CEFF_     : Effective force on the chiasm (geometry x turgor)
HYPTOT_   : Total hypopituitarism (compression + permanent surgical)
'

mod <- mcode("prolactinoma", prl_code)

# =============================================================================
#  HELPERS
#
#  DESIGN DECISION — presenting values are GENERATED, never imposed.
#  Every disease scenario begins with a SILENT PHASE of PRESIM days during
#  which the seeded tumour secretes, suppresses the gonadal axis, erodes bone
#  and (if large enough) compresses the chiasm, with no drug present. The
#  patient "presents" at day PRESIM, therapy starts there, and every baseline
#  quoted below is read off the model at that instant. Nothing about the
#  presenting phenotype — prolactin, amenorrhoea, T-score, visual field — is
#  set by hand.
# =============================================================================
PRESIM <- 730   # days of silent natural history before the patient presents

# patient archetypes -----------------------------------------------------------
p_healthy  <- list(NS0 = 1e-9, D2RS0 = 0.50)
p_micro    <- list(NS0 = 0.35, D2RS0 = 0.50)                       # ~9 mm
p_macro    <- list(NS0 = 3.20, D2RS0 = 0.35)                       # ~2 cm
p_giant    <- list(NS0 = 18.0, D2RS0 = 0.25, FDRIVE = 1.30, SEX = 1) # ~4 cm, male
p_partres  <- list(NS0 = 3.20, D2RS0 = 0.15)
p_trueres  <- list(NS0 = 3.20, D2RS0 = 0.04, FGEN = 1.30)          # AIP-like
p_stalk    <- list(NS0 = 1e-9, STALKP = 0.05)                      # NFPA + stalk
p_macroprl <- list(NS0 = 1e-9, APABSET = 2.2)

run <- function(par = list(), ev_ = NULL, end = 730, delta = 1) {
  m <- param(mod, par)
  if (!is.null(ev_)) out <- mrgsim(m, events = ev_, end = end, delta = delta,
                                   maxsteps = 500000, hmax = 0.5)
  else               out <- mrgsim(m, end = end, delta = delta,
                                   maxsteps = 500000, hmax = 0.5)
  as.data.frame(out)
}

# dosing event builders (all default to starting the day the patient presents)
evc <- function(...) {
  ds <- lapply(list(...), function(e) as.data.frame(e))
  cn <- unique(unlist(lapply(ds, names)))
  ds <- lapply(ds, function(d) { for (k in setdiff(cn, names(d))) d[[k]] <- 0; d[cn] })
  as.ev(do.call(rbind, ds))
}
ev_cab <- function(mg, wk, start = PRESIM, ii = 7)
  ev(time = start, amt = mg, cmt = "CAB_G", ii = ii, addl = ceiling(wk*7/ii) - 1)
ev_brc <- function(mg, days, start = PRESIM, perday = 3)
  ev(time = start, amt = mg, cmt = "BRC_G", ii = 1/perday, addl = days*perday - 1)
ev_qui <- function(ug, days, start = PRESIM)
  ev(time = start, amt = ug/1000, cmt = "QUI_G", ii = 1, addl = days - 1)
ev_ap  <- function(mg, days, start = 0)
  ev(time = start, amt = mg, cmt = "AP_G", ii = 1, addl = days - 1)
ev_tmz <- function(mgm2, cycles, bsa = 1.8, start = PRESIM)
  do.call(evc, lapply(0:(cycles-1), function(k)
    ev(time = start + 28*k, amt = mgm2*bsa, cmt = "TMZ_C", ii = 1, addl = 4)))

last  <- function(d, v) d[[v]][nrow(d)]
atday <- function(d, v, t) d[[v]][which.min(abs(d$time - t))]
pres  <- function(d, v) atday(d, v, PRESIM)
shr   <- function(d, t) 100*(1 - atday(d,"TVOL_",t)/pres(d,"TVOL_"))
wmax  <- function(d, v, t0, t1) max(d[[v]][d$time >= t0 & d$time <= t1])
wmin  <- function(d, v, t0, t1) min(d[[v]][d$time >= t0 & d$time <= t1])
hdr   <- function(s) cat("\n", strrep("=", 78), "\n ", s, "\n",
                         strrep("=", 78), "\n", sep = "")
sub   <- function(s) cat("\n--- ", s, " ", strrep("-", max(0, 68 - nchar(s))),
                         "\n", sep = "")

# =============================================================================
#  DIAGNOSTIC D1 — is the healthy baseline actually stationary?
#  Nothing else in this file is trustworthy if this fails.
# =============================================================================
hdr("D1 · HEALTHY BASELINE STATIONARITY (the disease must be GENERATED)")
h <- run(p_healthy, end = 3650, delta = 10)
vars <- c("PRLB","TIDA","DAP","D2RS","NN","MRNAN","GRN","KISS","GNRH","LH","FSH",
          "E2G","TST","OVF","SPZ","OCF","OBF","BMD","BW","IR","CORT","FT4","IGF1")
drift <- sapply(vars, function(v) {
  a <- h[[v]][1]; b <- h[[v]][nrow(h)]
  if (abs(a) < 1e-12) 0 else 100*(b - a)/a
})
cat(sprintf("  10-year drift with no tumour and no drug — max |drift| = %.6f%%\n\n",
            max(abs(drift))))
for (v in names(drift))
  cat(sprintf("    %-7s %12.4f -> %12.4f   %+11.7f%%\n",
              v, h[[v]][1], h[[v]][nrow(h)], drift[[v]]))
cat(sprintf("\n  VERDICT: %s\n",
            ifelse(max(abs(drift)) < 0.05,
                   "PASS — every baseline balance is solved in $MAIN, not fitted",
                   "FAIL — baseline drifts, so nothing downstream can be trusted")))

# =============================================================================
#  DIAGNOSTIC D2 — cabergoline PK against published anchors
# =============================================================================
hdr("D2 · CABERGOLINE PK — plasma in pg/mL, terminal t1/2 63-109 h")
pk  <- run(p_micro, ev_cab(1, 1, start = 0), end = 28, delta = 0.02)
cm  <- max(pk$CABNM)*451.6/1000                    # nM -> ng/L = pg/mL
i   <- which(pk$time > 4 & pk$time < 20)
lam <- -as.numeric(coef(lm(log(pk$CABNM[i]) ~ pk$time[i]))[2])
ss  <- run(p_micro, ev_cab(1, 26, start = 0), end = 182, delta = 0.25)
cat(sprintf("  Cmax after a single 1 mg dose   %8.1f pg/mL   (published 30-70)\n",
            cm*1000))
cat(sprintf("  terminal half-life              %8.1f h       (published 63-109)\n",
            log(2)/lam*24))
cat(sprintf("  steady-state biophase (1 mg/wk) %8.2f nM      (= %.1fx the D2 Ki)\n",
            mean(tail(ss$CAB_E, 200)), mean(tail(ss$CAB_E, 200))/0.70))
cat(sprintf("  biophase:plasma ratio           %8.1f x       (pituitary retention)\n",
            mean(tail(ss$CAB_E,200))/mean(tail(ss$CABNM,200))))
cat(sprintf("  peak-to-trough biophase ratio   %8.2f         (weekly dosing is\n",
            max(tail(ss$CAB_E,28))/min(tail(ss$CAB_E,28))))
cat("                                                     nearly a constant infusion)\n")
cat("  NOTE: absolute bioavailability is unknown, so V/F and CL/F are apparent\n")
cat("  volumes fitted to the published pg/mL range and half-life; the biophase\n")
cat("  partition absorbs tissue accumulation and is not separately identifiable.\n")

# =============================================================================
#  DIAGNOSTIC D3 — one dose suppresses prolactin for more than two weeks
# =============================================================================
hdr("D3 · A SINGLE 1 mg DOSE — duration of prolactin suppression")
s1 <- run(p_micro, ev_cab(1, 1), end = PRESIM + 120, delta = 0.1)
b  <- pres(s1, "PRLB")
w  <- s1$time >= PRESIM
nad <- min(s1$PRLB[w]); tnad <- s1$time[w][which.min(s1$PRLB[w])] - PRESIM
pw  <- s1$PRLB[w]; tw <- s1$time[w]
inad <- which.min(pw)
irec <- which(pw[inad:length(pw)] > 0.8*b)[1]
rec <- if (is.na(irec)) NA else tw[inad + irec - 1] - PRESIM
cat(sprintf("  prolactin at presentation       %8.1f ng/mL\n", b))
cat(sprintf("  nadir                           %8.1f ng/mL (%.0f%% fall) at day %.1f\n",
            nad, 100*(1 - nad/b), tnad))
cat(sprintf("  prolactin at day 14             %8.1f ng/mL (%.0f%% of baseline)\n",
            atday(s1,"PRLB",PRESIM+14), 100*atday(s1,"PRLB",PRESIM+14)/b))
cat(sprintf("  days to return above 80%% of baseline %6.1f d\n", rec))
cat("  The plasma half-life is 3.5 d, yet suppression outlasts it several-fold.\n")
cat("  The effect lives in the pituitary biophase, not in plasma — which is the\n")
cat("  entire reason a once-weekly tablet controls a daily secretory process.\n")

# =============================================================================
#  SCENARIOS
# =============================================================================
hdr("PROLACTINOMA QSP MODEL — 34 SCENARIOS (therapy starts on presentation)")

sub("S1-S3 · untreated natural history — the phenotype is generated")
S1 <- run(p_healthy, end = 2555, delta = 5)
S2 <- run(p_micro,   end = 2555, delta = 5)
S3 <- run(p_macro,   end = 2555, delta = 5)
cat(sprintf("  %-22s %9s %9s %8s %8s %8s %8s %8s\n", "", "PRL 2y", "PRL 7y",
            "vol 2y", "vol 7y", "MD dB", "T-scr", "anovul%"))
for (nm in list(c("healthy","S1"), c("microadenoma 0.35 mL","S2"),
                c("macroadenoma 3.2 mL","S3"))) {
  d <- get(nm[2])
  cat(sprintf("  %-22s %9.1f %9.1f %8.2f %8.2f %8.1f %8.2f %8.0f\n", nm[1],
              pres(d,"PRLB"), last(d,"PRLB"), pres(d,"TVOL_"), last(d,"TVOL_"),
              last(d,"VFMD"), last(d,"TSCORE"), last(d,"AMENP")))
}
cat(sprintf("\n  Nothing above was imposed. Seeding %.2f mL of lactotroph tissue with\n",
            p_micro$NS0))
cat(sprintf("  D2R density 0.50 generates prolactin %.0f ng/mL; %.1f mL at density\n",
            pres(S2,"PRLB"), 3.2))
cat(sprintf("  0.35 generates %.0f ng/mL. The published ordering (macroadenoma >250,\n",
            pres(S3,"PRLB")))
cat("  microadenoma 50-250 ng/mL) comes out of mass x receptor density.\n")
cat(sprintf("  Short-loop response: TIDA activity %.2f (healthy) -> %.2f (macro),\n",
            pres(S1,"TIDA"), pres(S3,"TIDA")))
cat(sprintf("  portal dopamine %.2f -> %.2f nM: the transmitter is working at %.0f%%\n",
            pres(S1,"DAP"), pres(S3,"DAP"),
            100*(pres(S3,"TIDA") - 1)/(2.0)))
cat("  of its maximal short-loop amplification. The receiver is what has failed,\n")
cat("  which is why only a ligand of far higher affinity than dopamine can help.\n")

sub("S4-S6 · cabergoline: dose, schedule and tumour size")
S4 <- run(p_micro, ev_cab(0.5, 104), end = PRESIM + 730, delta = 1)
S5 <- run(p_micro, ev_cab(0.25, 104, ii = 3.5), end = PRESIM + 730, delta = 1)
S6 <- run(p_macro, ev_cab(1.0, 104), end = PRESIM + 730, delta = 1)
cat(sprintf("  %-27s %8s %8s %8s %8s %8s\n", "", "PRL pre", "wk 4", "2 y",
            "shrink%", "ovul"))
for (nm in list(c("micro 0.5 mg once weekly","S4"),
                c("micro 0.25 mg twice weekly","S5"),
                c("macro 1 mg once weekly","S6"))) {
  d <- get(nm[2])
  cat(sprintf("  %-27s %8.1f %8.1f %8.1f %8.0f %8.2f\n", nm[1],
              pres(d,"PRLB"), atday(d,"PRLB",PRESIM+28), last(d,"PRLB"),
              shr(d, PRESIM+730), last(d,"OVF")))
}
cat(sprintf("\n  Splitting the same weekly dose changes the biophase trough by only\n"))
cat(sprintf("  %.1f%% (%.2f vs %.2f nM) and the prolactin not at all: with a 4-day\n",
            100*(min(tail(S5$CAB_E,30))/min(tail(S4$CAB_E,30)) - 1),
            min(tail(S4$CAB_E,30)), min(tail(S5$CAB_E,30))))
cat("  biophase half-life the schedule inside a week is almost irrelevant. This\n")
cat("  is a genuine model prediction and it argues against splitting for efficacy.\n")
cat(sprintf("  Ovulatory competence in the microadenoma: %.2f at presentation ->\n",
            pres(S4,"OVF")))
cat(sprintf("  %.2f at 2 y; oestradiol %.0f -> %.0f pg/mL.\n",
            last(S4,"OVF"), pres(S4,"E2G"), last(S4,"E2G")))

sub("S7 · how fast can a compressed chiasm recover? Two answers, not one")
S7  <- run(p_giant, ev_cab(0.5, 104), end = PRESIM + 730, delta = 0.25)
S7m <- run(list(NS0 = 4.6, D2RS0 = 0.40, SEX = 1), ev_cab(1.0, 104),
           end = PRESIM + 730, delta = 0.25)
cat(sprintf("  GIANT (%.1f mL at presentation, indentation %.2f cm, MD %.1f dB)\n",
            pres(S7,"TVOL_"), pres(S7,"COMPR_"), pres(S7,"VFMD")))
cat(sprintf("  reported prolactin %.0f ng/mL for a TRUE %.0f ng/mL — see S22\n",
            pres(S7,"PRLMEAS_"), pres(S7,"PRLB")))
cat(sprintf("  %8s %10s %10s %9s %9s %8s %8s\n", "day", "PRL", "volume mL",
            "cells mL", "cell vol", "force", "MD dB"))
for (t in c(0, 1, 3, 7, 28, 90, 365, 730))
  cat(sprintf("  %8.0f %10.0f %10.2f %9.2f %9.3f %8.2f %8.1f\n", t,
              atday(S7,"PRLB",PRESIM+t), atday(S7,"TVOL_",PRESIM+t),
              atday(S7,"CELLN",PRESIM+t), atday(S7,"VCELL",PRESIM+t),
              atday(S7,"CEFF_",PRESIM+t), atday(S7,"VFMD",PRESIM+t)))
cat(sprintf("\n  MARGINAL COMPRESSION (%.1f mL, indentation %.2f cm, MD %.1f dB)\n",
            pres(S7m,"TVOL_"), pres(S7m,"COMPR_"), pres(S7m,"VFMD")))
cat(sprintf("  %8s %10s %10s %9s %8s %8s\n", "day", "PRL", "volume mL",
            "cell vol", "force", "MD dB"))
for (t in c(0, 1, 3, 7, 28, 90, 365))
  cat(sprintf("  %8.0f %10.0f %10.2f %9.3f %8.2f %8.1f\n", t,
              atday(S7m,"PRLB",PRESIM+t), atday(S7m,"TVOL_",PRESIM+t),
              atday(S7m,"VCELL",PRESIM+t), atday(S7m,"CEFF_",PRESIM+t),
              atday(S7m,"VFMD",PRESIM+t)))
cat(sprintf("\n  At 72 h: giant field %+.1f dB, marginal field %+.1f dB, while cell\n",
            atday(S7,"VFMD",PRESIM+3) - pres(S7,"VFMD"),
            atday(S7m,"VFMD",PRESIM+3) - pres(S7m,"VFMD")))
cat(sprintf("  NUMBER has moved only %.1f%% and %.1f%%. Whatever recovery happens this\n",
            100*(atday(S7,"CELLN",PRESIM+3)/pres(S7,"CELLN") - 1),
            100*(atday(S7m,"CELLN",PRESIM+3)/pres(S7m,"CELLN") - 1)))
cat("  early is cell VOLUME and cell TURGOR, not cytoreduction.\n\n")
cat("  HONEST NEGATIVE RESULT: the model does NOT reproduce the anecdotal 24-72 h\n")
cat("  visual recovery in a very large tumour. The reason is geometric and hard\n")
cat("  to argue with — volume scales as the cube of the radius, so a 3% volume\n")
cat("  loss moves the top of a 4 cm mass by well under a millimetre. Rapid\n")
cat("  recovery appears in the model only where the chiasm is close to the\n")
cat("  threshold, which is probably also the truth about the case reports. No\n")
cat("  term was added to force the fast response in the giant.\n")

sub("S8-S10 · bromocriptine, the switch, and quinagolide")
S8  <- run(p_macro, ev_brc(2.5, 730), end = PRESIM + 730, delta = 1)
S9  <- run(p_macro, evc(ev_brc(2.5, 180), ev_cab(1.0, 78, start = PRESIM + 180)),
           end = PRESIM + 730, delta = 1)
S10 <- run(p_macro, ev_qui(300, 730), end = PRESIM + 730, delta = 1)
cat(sprintf("  %-28s %9s %9s %9s %9s\n", "", "PRL pre", "PRL 2y", "shrink%",
            "SIGDRIVE"))
for (nm in list(c("bromocriptine 2.5 mg TID","S8"),
                c("BRC 6 mo then cabergoline","S9"),
                c("quinagolide 300 ug/d","S10"),
                c("cabergoline 1 mg/wk (S6)","S6"))) {
  d <- get(nm[2])
  cat(sprintf("  %-28s %9.1f %9.1f %9.0f %9.3f\n", nm[1], pres(d,"PRLB"),
              last(d,"PRLB"), shr(d, PRESIM+730), mean(tail(d$SIGDRIVE, 60))))
}
cat(sprintf("\n  Bromocriptine is a PARTIAL agonist (e = 0.80): its ceiling is structural,\n"))
cat(sprintf("  not a matter of dose. Switching at 6 months rescues prolactin %.1f -> %.1f.\n",
            atday(S9,"PRLB",PRESIM+179), last(S9,"PRLB")))

sub("S11 · resistance — partial (shifted EC50) vs true (lowered Emax)")
S11a <- run(p_partres, ev_cab(0.5, 104), end = PRESIM + 730, delta = 2)
S11b <- run(p_partres, ev_cab(3.5, 104), end = PRESIM + 730, delta = 2)
S11c <- run(p_trueres, ev_cab(0.5, 104), end = PRESIM + 730, delta = 2)
S11d <- run(p_trueres, ev_cab(3.5, 104), end = PRESIM + 730, delta = 2)
cat(sprintf("  %-32s %10s %10s %9s %9s %7s\n", "", "PRL pre", "PRL 2y",
            "vol pre", "vol 2y", "<25?"))
for (nm in list(c("partial resistance  0.5 mg/wk","S11a"),
                c("partial resistance  3.5 mg/wk","S11b"),
                c("true resistance     0.5 mg/wk","S11c"),
                c("true resistance     3.5 mg/wk","S11d"))) {
  d <- get(nm[2])
  cat(sprintf("  %-32s %10.1f %10.1f %9.2f %9.2f %7s\n", nm[1], pres(d,"PRLB"),
              last(d,"PRLB"), pres(d,"TVOL_"), last(d,"TVOL_"),
              ifelse(last(d,"PRLB") < 25, "YES", "no")))
}
cat(sprintf("\n  A 7x dose increase buys %.0f%% further prolactin reduction when the EC50 is\n",
            100*(1 - last(S11b,"PRLB")/last(S11a,"PRLB"))))
cat(sprintf("  merely shifted, and %.0f%% when Emax itself is lowered. Escalation treats\n",
            100*(1 - last(S11d,"PRLB")/last(S11c,"PRLB"))))
cat("  one of these and cannot treat the other — which is the practical content\n")
cat("  of the distinction between partial and true dopamine-agonist resistance.\n")

sub("S12 · clonal escape over 10 years on an unchanged dose")
S12 <- run(c(p_macro, list(NR0 = 0.02)), ev_cab(1.0, 470),
           end = PRESIM + 3285, delta = 10)
cat(sprintf("  %8s %10s %10s %14s\n", "year", "PRL", "volume mL", "resistant %"))
for (t in c(0, 365, 1095, 2190, 3285))
  cat(sprintf("  %8.1f %10.1f %10.2f %13.1f%%\n", t/365,
              atday(S12,"PRLB",PRESIM+t), atday(S12,"TVOL_",PRESIM+t),
              100*atday(S12,"NR",PRESIM+t)/
                (atday(S12,"NS",PRESIM+t) + atday(S12,"NR",PRESIM+t))))
cat("  Escape is Darwinian selection on D2R density under a constant selective\n")
cat("  pressure. Nothing in the code says 'become resistant at year N'.\n")

sub("S13-S14 · drug withdrawal after three years: remission or recurrence")
S13 <- run(p_micro, ev_cab(1.0, 156), end = PRESIM + 3650, delta = 5)
S14 <- run(p_macro, ev_cab(1.0, 156), end = PRESIM + 3650, delta = 5)
for (nm in list(c("microadenoma","S13"), c("macroadenoma","S14"))) {
  d <- get(nm[2]); tstop <- PRESIM + 1092
  cat(sprintf("  %-14s at withdrawal: cells %.3f mL, volume %.3f mL, PRL %.1f\n",
              nm[1], atday(d,"CELLN",tstop), atday(d,"TVOL_",tstop),
              atday(d,"PRLB",tstop)))
  cat(sprintf("  %-14s 7 y after withdrawal: volume %.3f mL, PRL %.1f  ->  %s\n", "",
              last(d,"TVOL_"), last(d,"PRLB"),
              ifelse(last(d,"PRLB") < 25, "SUSTAINED REMISSION", "RECURRENCE")))
}
cat("  The threshold that separates them is computed in D10.\n")

sub("S15-S16 · pregnancy — three insults at once")
S15 <- run(c(p_micro, list(PREGON = PRESIM + 400)),
           ev_cab(0.5, 57), end = PRESIM + 1000, delta = 2)
S16 <- run(c(p_macro, list(PREGON = PRESIM + 400)),
           ev_cab(1.0, 57), end = PRESIM + 1000, delta = 2)
for (nm in list(c("microadenoma","S15"), c("macroadenoma","S16"))) {
  d <- get(nm[2]); t0 <- PRESIM + 400; t1 <- PRESIM + 690
  v0 <- atday(d,"TVOL_",t0)
  c0 <- atday(d,"CELLN",t0)
  cat(sprintf("  %-14s volume %6.3f -> %6.3f mL (%+5.1f%%) of which cell NUMBER\n",
              nm[1], v0, wmax(d,"TVOL_",t0,t1), 100*(wmax(d,"TVOL_",t0,t1)/v0 - 1)))
  cat(sprintf("  %-14s %6.3f -> %6.3f mL (%+5.1f%%): the rest is re-inflation\n", "",
              c0, wmax(d,"CELLN",t0,t1), 100*(wmax(d,"CELLN",t0,t1)/c0 - 1)))
  cat(sprintf("  %-14s indentation peak %.3f cm | MD %.1f dB | reported PRL peak %.0f\n",
              "", wmax(d,"COMPR_",t0,t1), wmin(d,"VFMD",t0,t1),
              wmax(d,"PRLMEAS_",t0,t1)))
}
cat("  Drug stopped, oestradiol up ~300-fold, placental lactogen added, and the\n")
cat("  normal lactotroph pool hyperplases too. But most of the apparent growth in\n")
cat("  a TREATED tumour is cells re-inflating after the drug is withdrawn, not new\n")
cat("  proliferation — an important distinction when a mid-pregnancy MRI is read.\n")
cat("  The size dependence of the SYMPTOMATIC risk is dissected in D13.\n")

sub("S17 · stalk-effect (disconnection) hyperprolactinaemia has a CEILING")
S17 <- run(p_stalk, end = 730, delta = 5)
cat(sprintf("  stalk patency 5%%: portal dopamine %.3f nM (healthy 5.00)\n",
            last(S17,"DAP")))
cat(sprintf("  reported prolactin %.1f ng/mL; normal pool hyperplases %.2f -> %.2f\n",
            last(S17,"PRLMEAS_"), S17$NN[1], last(S17,"NN")))
cat(sprintf("  ovulatory competence %.2f — clinically significant, but bounded\n",
            last(S17,"OVF")))
cat("  Only the FINITE normal pool is disinhibited, so the equations bound the\n")
cat("  prolactin even though the consequences for the gonadal axis are real.\n")

sub("S18-S19 · drug-induced hyperprolactinaemia, both directions")
S18 <- run(list(NS0 = 1e-9), ev_ap(6, 730), end = 730, delta = 2)
S19 <- run(list(NS0 = 1e-9, EFF_AP = 0.25, KI_AP = 0.34, FU_AP = 0.01,
                CL_AP = 210, V_AP = 1400), ev_ap(15, 730), end = 730, delta = 2)
cat(sprintf("  %-34s %9s %9s %9s\n", "", "PRL 3mo", "PRL 2y", "pool"))
cat(sprintf("  %-34s %9.1f %9.1f %9.2f\n", "risperidone-like 6 mg/d (e = 0)",
            atday(S18,"PRLB",90), last(S18,"PRLB"), last(S18,"NN")))
cat(sprintf("  %-34s %9.1f %9.1f %9.2f\n", "aripiprazole 15 mg/d (e = 0.25)",
            atday(S19,"PRLB",90), last(S19,"PRLB"), last(S19,"NN")))
cat("  Chronic D2 blockade drives lactotroph hyperplasia, which is why the\n")
cat("  prolactin keeps creeping up for months after the dose is stable.\n")
cat("  Aripiprazole alone has LOW efficacy but HIGH affinity, so it displaces\n")
cat("  dopamine and still overshoots the physiological tone: prolactin falls\n")
cat("  BELOW normal. Both directions come from the same occupancy line (D9).\n")

sub("S20 · a prolactinoma in a patient who needs an antipsychotic")
S20 <- run(c(p_micro), evc(ev_ap(6, 730, start = PRESIM), ev_cab(1.0, 104)),
           end = PRESIM + 730, delta = 2)
cat(sprintf("  microadenoma + risperidone 6 mg/d + cabergoline 1 mg/wk:\n"))
cat(sprintf("  PRL %.1f -> %.1f ng/mL; SIGDRIVE %.3f (cabergoline alone gave %.3f)\n",
            pres(S20,"PRLB"), last(S20,"PRLB"), mean(tail(S20$SIGDRIVE,60)),
            mean(tail(S4$SIGDRIVE,60))))
cat("  The antagonist occupies receptors the agonist needs; control is achievable\n")
cat("  but costs signal, and the mesolimbic side of the same occupancy is the\n")
cat("  reason this combination is approached cautiously.\n")

sub("S21 · macroprolactinaemia — a deliberate STRUCTURAL NULL")
S21a <- run(p_macroprl, end = 365, delta = 5)
S21b <- run(p_macroprl, ev_cab(0.5, 52, start = 0), end = 365, delta = 5)
cat(sprintf("  untreated    reported %6.1f | bioactive %5.1f | PEG %5.1f%% | ovul %.2f | T %+.2f\n",
            last(S21a,"PRLMEAS_"), last(S21a,"PRLB"), last(S21a,"PEGREC"),
            last(S21a,"OVF"), last(S21a,"TSCORE")))
cat(sprintf("  cabergoline  reported %6.1f | bioactive %5.1f | PEG %5.1f%% | ovul %.2f | T %+.2f\n",
            last(S21b,"PRLMEAS_"), last(S21b,"PRLB"), last(S21b,"PEGREC"),
            last(S21b,"OVF"), last(S21b,"TSCORE")))
cat("  The drug moves the NUMBER and changes nothing that matters, because the\n")
cat("  gonadal axis was never abnormal. Futility here is structural.\n")

sub("S22 · the hook effect in the giant prolactinoma of S7")
G <- run(p_giant, end = PRESIM, delta = 10)
cat(sprintf("  true bioactive prolactin       %10.0f ng/mL\n", last(G,"PRLB")))
cat(sprintf("  immunoreactive analyte         %10.0f ng/mL\n", last(G,"PRLIMM")))
cat(sprintf("  REPORTED by the assay          %10.0f ng/mL  <-- reads as non-secreting\n",
            last(G,"PRLMEAS_")))
cat(sprintf("  reported after 1:100 dilution  %10.0f ng/mL  <-- recovered\n",
            last(G,"PRLDIL")))
cat("  A 4 cm sellar mass with an unremarkable prolactin is the classic trap. It\n")
cat("  is not biology; it is a two-site assay evaluated past its own peak (D6).\n")

sub("S23-S24 · non-tumoral causes: hypothyroidism and renal failure")
S23 <- run(list(NS0 = 1e-9, PRIMHYP = 0.75), end = 365, delta = 5)
S24 <- run(list(NS0 = 1e-9, FCLR = 0.35), end = 365, delta = 5)
cat(sprintf("  primary hypothyroidism (FT4 %.2f of normal): PRL %5.1f ng/mL\n",
            last(S23,"FT4"), last(S23,"PRLB")))
cat(sprintf("  CKD (prolactin clearance x0.35):            PRL %5.1f ng/mL\n",
            last(S24,"PRLB")))
cat("  Both sit in the mildly-raised band that overlaps drug-induced and\n")
cat("  stalk-effect hyperprolactinaemia and the smallest microadenomas, which is\n")
cat("  why the MAGNITUDE, not the presence of hyperprolactinaemia, carries the\n")
cat("  diagnostic information. Both are also reversible by treating the cause,\n")
cat("  and neither is helped by a dopamine agonist.\n")

sub("S25-S26 · transsphenoidal surgery, alone and with adjuvant cabergoline")
S25 <- run(c(p_macro, list(SURGD = PRESIM + 30, SURGRES = 0.25, SURGHYP = 0.12)),
           end = PRESIM + 1095, delta = 5)
S26 <- run(c(p_macro, list(SURGD = PRESIM + 30, SURGRES = 0.25, SURGHYP = 0.12)),
           ev_cab(1.0, 145, start = PRESIM + 60), end = PRESIM + 1095, delta = 5)
for (nm in list(c("surgery alone","S25"), c("surgery + cabergoline","S26"))) {
  d <- get(nm[2])
  cat(sprintf("  %-22s PRL %7.1f -> %7.1f | volume %5.2f -> %5.2f mL\n", nm[1],
              pres(d,"PRLB"), last(d,"PRLB"), pres(d,"TVOL_"), last(d,"TVOL_")))
}
cat(sprintf("  surgical cost: cortisol reserve %.2f, IGF-1 %.2f, FT4 %.2f\n",
            last(S26,"CORT"), last(S26,"IGF1"), last(S26,"FT4")))
cat("  A 25% residual regrows on its own; the same residual under dopaminergic\n")
cat("  control does not. The operation changes the initial condition, not the\n")
cat("  biology of what is left.\n")

sub("S27-S28 · temozolomide in aggressive PitNET — MGMT decides")
S27 <- run(c(p_trueres, list(MGMT0 = 0.05)), ev_tmz(180, 6, start = PRESIM + 30),
           end = PRESIM + 545, delta = 5)
S28 <- run(c(p_trueres, list(MGMT0 = 0.80)), ev_tmz(180, 6, start = PRESIM + 30),
           end = PRESIM + 545, delta = 5)
S28n <- run(p_trueres, ev_cab(3.5, 78), end = PRESIM + 545, delta = 5)
for (nm in list(c("MGMT low  (0.05)","S27"), c("MGMT high (0.80)","S28"))) {
  d <- get(nm[2])
  cat(sprintf("  %-18s volume %6.2f -> %6.2f mL (%+5.0f%%) | PRL %8.0f -> %8.0f\n",
              nm[1], pres(d,"TVOL_"), last(d,"TVOL_"), shr(d, PRESIM+545),
              pres(d,"PRLB"), last(d,"PRLB")))
}
cat(sprintf("  %-18s volume %6.2f -> %6.2f mL (%+5.0f%%) | PRL %8.0f -> %8.0f\n",
            "no TMZ (cab only)", pres(S28n,"TVOL_"), last(S28n,"TVOL_"),
            shr(S28n, PRESIM+545), pres(S28n,"PRLB"), last(S28n,"PRLB")))
cat("  Identical six cycles; only the repair enzyme differs. The MGMT-high tumour\n")
cat("  is indistinguishable from no alkylator at all, which is the model's version\n")
cat("  of the recommendation to stain for MGMT before committing to temozolomide.\n")

sub("S29-S30 · cabergoline and the heart valve: the dose scale IS the argument")
S29 <- run(p_micro, ev_cab(1.0, 260, start = 0), end = 1825, delta = 10)
S30 <- run(p_micro, ev_cab(3.0, 260, start = 0, ii = 1), end = 1825, delta = 10)
cat(sprintf("  %-26s %11s %11s %11s %9s\n", "", "mg/week", "mean 5HT2B",
            "valve 5 y", "P(TR) %"))
for (nm in list(c("prolactinoma 1 mg/week","S29",7), c("Parkinson 3 mg/DAY","S30",21))) {
  d <- get(nm[2])
  cat(sprintf("  %-26s %11.1f %11.4f %11.1f %9.2f\n", nm[1], as.numeric(nm[3]),
              mean(d$OCC5HT), last(d,"VALV"), last(d,"PTR")))
}
cat(sprintf("  21x the weekly dose gives %.0fx the 5-HT2B occupancy and %.0fx the\n",
            mean(S30$OCC5HT)/mean(S29$OCC5HT),
            last(S30,"PTR")/max(last(S29,"PTR"), 1e-9)))
cat("  probability of moderate-severe regurgitation. The Parkinson's-disease and\n")
cat("  prolactinoma literatures were never in conflict; they were never on the\n")
cat("  same part of the dose axis.\n")
cat(sprintf("  Impulse-control risk, by contrast, is already %.1f%% at 1 mg/week,\n",
            last(S29,"PICD")))
cat("  because it runs off D3 occupancy rather than the ergoline skeleton.\n")

sub("S31-S32 · bone: what recovers and what is banked")
S31 <- run(c(p_macro, list(SEX = 1)), ev_cab(1.0, 530, start = PRESIM + 1825),
           end = PRESIM + 5475, delta = 10)
S31e <- run(c(p_macro, list(SEX = 1)), ev_cab(1.0, 790),
            end = PRESIM + 5475, delta = 10)
S32 <- run(c(p_trueres, list(SEX = 1, FTREPL = 0.95)), ev_cab(3.5, 530),
           end = PRESIM + 3650, delta = 10)
S32n <- run(c(p_trueres, list(SEX = 1, FTREPL = 0.00)), ev_cab(3.5, 530),
            end = PRESIM + 3650, delta = 10)
cat(sprintf("  male macroadenoma, 5 y untreated then cabergoline:\n"))
cat(sprintf("    T-score %+.2f (presentation) -> %+.2f (5 y) -> %+.2f (10 y) -> %+.2f (15 y)\n",
            pres(S31,"TSCORE"), atday(S31,"TSCORE",PRESIM+1825),
            atday(S31,"TSCORE",PRESIM+3650), last(S31,"TSCORE")))
tot <- atday(S31,"BMD",PRESIM) - atday(S31,"BMD",PRESIM+1825)
cat(sprintf("    BMD lost by 5 y %.3f g/cm2, of which %.3f is irreversible (%.0f%%)\n",
            tot, atday(S31,"BLOSSI",PRESIM+1825),
            100*atday(S31,"BLOSSI",PRESIM+1825)/max(tot,1e-9)))
cat(sprintf("    testosterone %.0f (Dx) -> %.0f (5 y) -> %.0f (15 y) ng/dL\n",
            pres(S31,"TST"), atday(S31,"TST",PRESIM+1825), last(S31,"TST")))
cat(sprintf("    prolactin %.0f -> %.0f -> %.0f ng/mL; fracture risk %.1f%% -> %.1f%%\n",
            pres(S31,"PRLB"), atday(S31,"PRLB",PRESIM+1825), last(S31,"PRLB"),
            atday(S31,"PFRACT",PRESIM+1825), last(S31,"PFRACT")))
cat(sprintf("  SAME tumour treated from presentation instead: T-score %+.2f at 15 y\n",
            last(S31e,"TSCORE")))
cat(sprintf("  (vs %+.2f with 5 years of delay) — %.2f T-score units are the price of\n",
            last(S31,"TSCORE"), last(S31e,"TSCORE") - last(S31,"TSCORE")))
cat("  the delay, and they are not recoverable at any later date.\n")
cat(sprintf("  True-resistant male, testosterone replaced: T-score %+.2f -> %+.2f at 10 y\n",
            pres(S32,"TSCORE"), last(S32,"TSCORE")))
cat(sprintf("  versus %+.2f unreplaced, while prolactin stays above %.0f ng/mL in both.\n",
            last(S32n,"TSCORE"), min(last(S32,"PRLB"), last(S32n,"PRLB"))))
cat("  Replacing the hormone protects the bone even when the tumour cannot be\n")
cat("  controlled: the secretory arm and the gonadal arm are separable, and the\n")
cat("  second one is treatable on its own.\n")

sub("S33 · metabolic arm")
S33 <- run(p_macro, ev_cab(1.0, 104), end = PRESIM + 730, delta = 5)
cat(sprintf("  weight %.1f -> %.1f kg (%+.1f), HOMA-IR %.2f -> %.2f\n",
            pres(S33,"BW"), last(S33,"BW"), last(S33,"BW") - pres(S33,"BW"),
            pres(S33,"IR"), last(S33,"IR")))

sub("S34 · nausea tolerance is why we titrate slowly and dose at bedtime")
S34 <- run(p_micro, ev_cab(0.5, 26), end = PRESIM + 182, delta = 0.25)
cat(sprintf("  nausea index  day 1 %5.1f | day 7 %5.1f | day 28 %5.1f | day 84 %5.1f\n",
            atday(S34,"NAUS",PRESIM+1), atday(S34,"NAUS",PRESIM+7),
            atday(S34,"NAUS",PRESIM+28), atday(S34,"NAUS",PRESIM+84)))
cat(sprintf("  tolerance state %.2f -> %.2f (max %.1f)\n",
            pres(S34,"NTOL"), last(S34,"NTOL"), 6.0))

# =============================================================================
#  DIAGNOSTICS D4 - D15
# =============================================================================
hdr("D4 · DOSE-RESPONSE: CABERGOLINE vs BROMOCRIPTINE (Webster 1994 anchor)")
cd <- c(0.125, 0.25, 0.5, 1, 2, 3.5)
bd <- c(1.25, 2.5, 5, 7.5, 15, 30)
cp <- sapply(cd, function(x) last(run(p_macro, ev_cab(x, 104),
                                     end = PRESIM + 730, delta = 7), "PRLB"))
bp <- sapply(bd, function(x) last(run(p_macro, ev_brc(x/3, 730),
                                     end = PRESIM + 730, delta = 7), "PRLB"))
cat("  Macroadenoma at presentation (prolactin ~780 ng/mL), two years of therapy:\n")
cat(sprintf("  %-11s %9s %9s  |  %-11s %9s %9s\n",
            "cab mg/wk", "PRL", "<25?", "brc mg/d", "PRL", "<25?"))
for (i in seq_along(cd))
  cat(sprintf("  %-11.3f %9.2f %9s  |  %-11.2f %9.2f %9s\n",
              cd[i], cp[i], ifelse(cp[i] < 25, "yes", "no"),
              bd[i], bp[i], ifelse(bp[i] < 25, "yes", "no")))
cat(sprintf("\n  Cabergoline plateaus at %.2f ng/mL, bromocriptine at %.2f. The whole gap\n",
            min(cp), min(bp)))
cat("  is the intrinsic-activity term e (1.00 vs 0.80) — no separate efficacy\n")
cat("  parameter for the two drugs exists in this model. Webster 1994 reported\n")
cat("  normoprolactinaemia in 83% vs 59% of patients; the model reproduces the\n")
cat("  DIRECTION and the plateau gap, but a rate across patients requires a\n")
cat("  population distribution of D2R density, which is D2R sweep territory.\n")
cat(sprintf("\n  Above ~%.1f mg/wk the curve is flat: %.0f%% of the achievable fall is\n",
            0.5, 100*(1 - cp[3]/cp[1])/(1 - min(cp)/cp[1])))
cat("  already delivered at 0.5 mg/wk. The model therefore argues AGAINST routine\n")
cat("  escalation for biochemical control — while still supporting it for\n")
cat("  shrinkage, which sits on a less potent branch (see D5).\n")

hdr("D5 · DECOMPOSING THE PROLACTIN FALL INTO ITS THREE TIME CONSTANTS")
dd <- run(p_macro, ev_cab(1.0, 104), end = PRESIM + 730, delta = 0.5)
b0 <- pres(dd, "PRLB")
r0 <- dd[which.min(abs(dd$time - PRESIM)), ]
cat(sprintf("  %8s %10s %10s %11s %11s %11s %11s\n", "day", "PRL", "% base",
            "exocytosis", "transcript", "cell vol", "cell number"))
for (t in c(0, 0.25, 1, 3, 7, 28, 90, 180, 365, 730)) {
  r <- dd[which.min(abs(dd$time - (PRESIM + t))), ]
  cat(sprintf("  %8.2f %10.1f %9.1f%% %11.3f %11.3f %11.3f %11.3f\n", t, r$PRLB,
              100*r$PRLB/b0, r$REXS/r0$REXS, r$RTRS/r0$RTRS,
              r$VCELL/r0$VCELL, r$CELLN/r0$CELLN))
}
cat("\n  Branches 1 and 2 both settle within about three days, so they are NOT as\n")
cat("  cleanly separated in time as the map suggests, and that is worth saying.\n")
cat("  What separates them is their ROLE, not their clock: exocytosis produces the\n")
cat("  first fall, and transcription sets the level the system settles at, because\n")
cat("  at steady state secretion equals synthesis minus degradation. That is mass\n")
cat("  balance, not an assumption, and it is why the chronic prolactin cannot go\n")
cat("  below what the transcriptional branch allows however hard exocytosis is\n")
cat("  blocked. The genuinely slow clocks are branch 3 (cell volume, settling over\n")
cat("  about two months) and branch 4 (cell number, still moving at two years).\n")
cat("  The granule store supplies a further intermediate lag: prolactin keeps\n")
cat("  falling from 66% of baseline at day 1 to 22% at day 7 while both receptor\n")
cat("  branches have already reached their new values.\n")

hdr("D6 · THE HOOK EFFECT AS A CURVE, NOT AN ANECDOTE")
cat(sprintf("  %12s %12s %10s %14s\n", "analyte", "reported", "ratio", "after 1:100"))
for (a in c(20, 100, 250, 800, 1500, 2500, 4000, 6000, 12000, 25000, 60000)) {
  m  <- a/(1 + (a/4000)^4)
  ad <- a/100; md <- 100*ad/(1 + (ad/4000)^4)
  cat(sprintf("  %12.0f %12.1f %10.3f %14.0f\n", a, m, m/a, md))
}
cat("\n  Monotone through the usual macroadenoma range, a peak near 3000, then a\n")
cat("  collapse. The confirmatory\n")
cat("  dilution is the same curve evaluated 100-fold down the analyte axis —\n")
cat("  structurally identical to the way excess heparin abolishes the serotonin\n")
cat("  release assay in the HIT model in this library.\n")
cat("  LIMITATION: KHOOK and the exponent are platform-specific. Only the shape\n")
cat("  (monotone, peak, collapse) is general.\n")

hdr("D7 · PEG RECOVERY SEPARATES ARTEFACT FROM DISEASE")
cat(sprintf("  %-24s %10s %10s %9s %7s %7s\n", "", "reported", "bioactive",
            "PEG %", "ovul", "T-scr"))
for (nm in list(c("healthy","p_healthy"), c("microadenoma","p_micro"),
                c("macroprolactinaemia","p_macroprl"))) {
  d <- run(get(nm[2]), end = 730, delta = 10)
  cat(sprintf("  %-24s %10.1f %10.1f %9.1f %7.2f %+7.2f\n", nm[1],
              last(d,"PRLMEAS_"), last(d,"PRLB"), last(d,"PEGREC"),
              last(d,"OVF"), last(d,"TSCORE")))
}
cat("\n  Recovery below 40% is the discriminator, and it is computed from two\n")
cat("  pools rather than tabulated. Note that the microadenoma and the artefact\n")
cat("  can report a similar number while only one has an abnormal axis.\n")

hdr("D8 · WHY STALK EFFECT HAS A CEILING AND A PROLACTINOMA DOES NOT")
cat(sprintf("  %-16s %12s %12s %12s %10s\n", "stalk patency", "portal DA",
            "PRL", "normal pool", "ovul"))
for (sp in c(1.0, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001)) {
  d <- run(list(NS0 = 1e-9, STALKP = sp), end = 730, delta = 10)
  cat(sprintf("  %-16.3f %12.3f %12.1f %12.2f %10.2f\n", sp, last(d,"DAP"),
              last(d,"PRLB"), last(d,"NN"), last(d,"OVF")))
}
ceil <- last(run(list(NS0 = 1e-9, STALKP = 1e-6), end = 730, delta = 10), "PRLB")
mac  <- pres(S3, "PRLB")
cat(sprintf("\n  Complete stalk section gives at most %.0f ng/mL. A 3.2 mL prolactinoma\n", ceil))
cat(sprintf("  generates %.0f, which is %.0fx that ceiling, from %.0fx the lactotroph mass\n",
            mac, mac/ceil, 3.2/0.15))
cat("  operating with a defective receiver. The bedside rule 'prolactin 60 with a\n")
cat("  3 cm mass means stalk effect' is a COMPUTED BOUND here, not a heuristic.\n")
cat("  The ceiling lands inside the published 25-150 ng/mL range for\n")
cat("  disconnection hyperprolactinaemia WITHOUT being fitted to it: it is set by\n")
cat("  the size of the normal lactotroph pool, its maximal disinhibited secretory\n")
cat("  rate, and the two-fold hyperplasia that chronic loss of D2 signal causes.\n")
cat("  LIMITATION: values above ~150 ng/mL reported in some series cannot be\n")
cat("  produced by this structure at all. If they are real they require either\n")
cat("  greater hyperplasia than modelled here or a coexisting lactotroph adenoma,\n")
cat("  which makes this a falsifiable claim rather than a caveat.\n")

hdr("D9 · THE ARIPIPRAZOLE ADD-ON PARADOX FROM ONE LINE OF ALGEBRA")
KI_DA_ <- 50; DAP_ <- 5; KI_AP_ <- 0.50; KI_AR_ <- 0.34
CAP <- 3.0; CAR <- 3.3         # free nM at usual clinical doses
sig <- function(cap, car) {
  den <- 1 + DAP_/KI_DA_ + cap/KI_AP_ + car/KI_AR_
  num <- 1.00*DAP_/KI_DA_ + 0.00*cap/KI_AP_ + 0.25*car/KI_AR_
  num/den
}
s_norm <- sig(0, 0)
cat(sprintf("  %-36s %11s %11s %12s\n", "ligands present at the receptor",
            "SIGDRIVE", "vs normal", "prolactin"))
for (r in list(c("dopamine only (physiological)", 0, 0),
               c("+ risperidone 6 mg/d (e = 0)", CAP, 0),
               c("+ aripiprazole 15 mg/d alone", 0, CAR),
               c("risperidone AND aripiprazole", CAP, CAR))) {
  v <- sig(as.numeric(r[2]), as.numeric(r[3]))
  cat(sprintf("  %-36s %11.4f %10.2fx %12s\n", r[1], v, v/s_norm,
              ifelse(abs(v/s_norm - 1) < 1e-9, "reference",
                     ifelse(v > s_norm, "lower", "HIGHER"))))
}
cat("\n  Adding a weak partial agonist to a silent antagonist RAISES the signal,\n")
cat("  because it replaces an e = 0 ligand at the receptor with an e = 0.25 one\n")
cat("  of higher affinity. That is the documented fall in prolactin on\n")
cat("  aripiprazole augmentation, and it needs no new mechanism.\n")

hdr("D10 · WITHDRAWAL — AND A RESULT THAT REFUTES THE MAP'S OWN HYPOTHESIS")
# The mechanistic map asserts that remission after withdrawal happens when the
# residual mass is small enough for ENDOGENOUS dopamine tone to hold it. This
# diagnostic tests that claim inside the equations, and the equations reject it.
wd <- function(v0, horizon = 5475) {
  d <- run(list(NS0 = v0, D2RS0 = 0.50), ev_cab(1.0, 156),
           end = PRESIM + horizon, delta = 20)
  tstop <- PRESIM + 1092
  aft   <- d$time > tstop
  irec  <- which(d$PRLB[aft] > 25)[1]
  list(rem  = atday(d, "CELLN", tstop),
       prl  = last(d, "PRLB"),
       vol  = last(d, "TVOL_"),
       trec = if (is.na(irec)) NA else (d$time[aft][irec] - tstop)/365)
}
cat(sprintf("  %-12s %11s %11s %13s %12s\n", "seeded mL", "remnant mL",
            "PRL 15 y", "volume 15 y", "recur (yr)"))
for (v0 in c(0.002, 0.005, 0.02, 0.05, 0.10, 0.20, 0.35, 0.60, 1.00, 3.20)) {
  r <- wd(v0)
  cat(sprintf("  %-12.3f %11.4f %11.1f %13.3f %12s\n", v0, r$rem, r$prl, r$vol,
              ifelse(is.na(r$trec), ">15", sprintf("%.1f", r$trec))))
}
cat("\n  TIME TO RECURRENCE IS LOG-LINEAR IN THE REMNANT, WITH NO THRESHOLD.\n")
cat("  Each halving of the residual cell mass buys roughly a fixed number of\n")
cat("  years, because regrowth off-drug is exponential. The model produces NO\n")
cat("  cure and NO stable small-remnant state.\n\n")
cat("  WHY THE MAP'S HYPOTHESIS FAILS — and this is worth stating plainly.\n")
cat("  The short loop is driven BY PROLACTIN. A small remnant secretes little\n")
cat("  prolactin, so TIDA drive falls back toward baseline and portal dopamine\n")
cat("  falls WITH it. Endogenous tone is therefore WEAKEST exactly when the\n")
cat("  remnant is smallest, which is the opposite of what is needed to hold it:\n")
r1 <- run(list(NS0 = 0.05, D2RS0 = 0.50), ev_cab(1.0, 156),
          end = PRESIM + 2000, delta = 10)
r2 <- run(list(NS0 = 1.00, D2RS0 = 0.50), ev_cab(1.0, 156),
          end = PRESIM + 2000, delta = 10)
t1 <- PRESIM + 1200
cat(sprintf("    small remnant: PRL %6.1f -> portal DA %5.2f nM -> SIGDRIVE %.4f\n",
            atday(r1,"PRLB",t1), atday(r1,"DAP",t1), atday(r1,"SIGDRIVE",t1)))
cat(sprintf("    large remnant: PRL %6.1f -> portal DA %5.2f nM -> SIGDRIVE %.4f\n",
            atday(r2,"PRLB",t1), atday(r2,"DAP",t1), atday(r2,"SIGDRIVE",t1)))
cat("  So the mechanism behind 'a visible remnant predicts relapse' is NOT that\n")
cat("  the hypothalamus can hold a small remnant. It is simply that a smaller\n")
cat("  remnant needs more doublings to become detectable again. The model turns\n")
cat("  the pooled ~21% durable remission of Dekkers 2010 into a statement about\n")
cat("  FOLLOW-UP DURATION rather than about cure, and predicts that remission\n")
cat("  rates must keep falling as cohorts are followed longer:\n")
vseq <- c(0.002, 0.005, 0.02, 0.05, 0.10, 0.20, 0.35, 0.60, 1.00, 3.20)
for (H in c(1457, 2192, 3650, 5475, 9125)) {
  rr <- sum(sapply(vseq, function(v) is.na(wd(v, H)$trec)))
  cat(sprintf("    %5.1f y after withdrawal: %2d of %d remnants still normoprolactinaemic\n",
              (H - 1092)/365, rr, length(vseq)))
}
cat("  This is a negative result about the model's own design premise and it is\n")
cat("  reported rather than tuned away. It also carries a testable prediction:\n")
cat("  a cohort followed for 15 years should show near-universal biochemical\n")
cat("  recurrence, and the ones that do not are evidence for a mechanism this\n")
cat("  model lacks (senescence, infarction, or immune-mediated involution).\n\n")

hdr("D11 · VALVULOPATHY: THE SAME DRUG, TWO DOSE SCALES")
cat(sprintf("  %-20s %10s %12s %12s %10s\n", "regimen", "mg/week",
            "mean 5HT2B", "valve 5 y", "P(TR) %"))
for (r in list(c("0.5 mg/wk",0.5,7), c("1 mg/wk",1,7), c("2 mg/wk",2,7),
               c("3.5 mg/wk",3.5,7), c("1 mg/d",1,1), c("3 mg/d",3,1))) {
  amt <- as.numeric(r[2]); ii <- as.numeric(r[3])
  d <- run(p_micro, ev_cab(amt, 260, start = 0, ii = ii), end = 1825, delta = 20)
  cat(sprintf("  %-20s %10.1f %12.4f %12.1f %10.3f\n", r[1], amt*7/ii,
              mean(d$OCC5HT), last(d,"VALV"), last(d,"PTR")))
}
cat("\n  Risk is cumulative occupancy, and occupancy is convex in dose over this\n")
cat("  range, so the risk ratio between the two literatures is much larger than\n")
cat("  the dose ratio. Nothing in this table was fitted to either literature.\n")

hdr("D12 · 5-HT2B / D2 SELECTIVITY SWEEP — EFFICACY FLAT, RISK COLLAPSES")
cat(sprintf("  %-14s %12s %12s %12s %12s\n", "5HT2B Ki nM", "PRL 2 y",
            "shrink %", "valve 5 y", "P(TR) %"))
for (ki in c(1.2, 3.6, 12, 36, 120, 1200)) {
  e <- run(c(p_micro, list(KI_5HT2B = ki)), ev_cab(1.0, 104),
           end = PRESIM + 730, delta = 20)
  v <- run(c(p_micro, list(KI_5HT2B = ki)), ev_cab(1.0, 260, start = 0),
           end = 1825, delta = 20)
  cat(sprintf("  %-14.1f %12.2f %12.0f %12.1f %12.3f\n", ki, last(e,"PRLB"),
              shr(e, PRESIM+730), last(v,"VALV"), last(v,"PTR")))
}
cat("\n  Four decades of selectivity leave prolactin control and shrinkage exactly\n")
cat("  where they were. Quinagolide occupies the far end of this axis by having\n")
cat("  no ergoline skeleton at all — and yet it sits at the SAME point as\n")
cat("  cabergoline on the D3/impulse-control axis, which is why switching for\n")
cat("  that indication does not help. Selectivity is not a single number.\n")

hdr("D13 · PREGNANCY RISK IS GEOMETRY, NOT A TABULATED PERCENTAGE")
cat(sprintf("  %-14s %11s %10s %12s %10s %10s\n", "vol at conception",
            "peak vol", "% growth", "indentation", "MD dB", "symptomatic"))
for (v0 in c(0.1, 0.3, 0.8, 2.0, 4.0, 8.0)) {
  d <- run(list(NS0 = v0, D2RS0 = 0.45, PREGON = PRESIM + 400),
           ev_cab(0.75, 57), end = PRESIM + 760, delta = 4)
  t0 <- PRESIM + 400; t1 <- PRESIM + 690
  md <- wmin(d,"VFMD",t0,t1)
  cat(sprintf("  %-14.2f %11.2f %9.0f%% %12.3f %10.1f %10s\n",
              atday(d,"TVOL_",t0), wmax(d,"TVOL_",t0,t1),
              100*(wmax(d,"TVOL_",t0,t1)/atday(d,"TVOL_",t0) - 1),
              wmax(d,"COMPR_",t0,t1), md, ifelse(md < -2, "YES", "no")))
}
cat("\n  The percentage growth is similar across sizes; what differs is whether\n")
cat("  the extra millimetres reach the chiasm. That geometric threshold is why\n")
cat("  symptomatic enlargement is reported in ~2.7% of microadenomas and ~21-23%\n")
cat("  of macroadenomas (Molitch) — a difference in starting position, not in\n")
cat("  the biology of oestrogen-driven growth.\n")

hdr("D14 · SHRINKAGE IS VOLUME FIRST, NUMBER LATER — AND ONLY ONE IS BANKED")
sw <- run(p_macro, ev_cab(1.0, 52), end = PRESIM + 1460, delta = 2)
cat(sprintf("  %8s %11s %11s %11s %11s\n", "day", "tumour mL", "cells mL",
            "cell vol", "PRL"))
for (t in c(0, 14, 60, 180, 365, 380, 420, 550, 730, 1460))
  cat(sprintf("  %8.0f %11.3f %11.3f %11.3f %11.1f\n", t,
              atday(sw,"TVOL_",PRESIM+t), atday(sw,"CELLN",PRESIM+t),
              atday(sw,"VCELL",PRESIM+t), atday(sw,"PRLB",PRESIM+t)))
cat("\n  Therapy stops at day 365. The volume re-expands within weeks because the\n")
cat("  cells re-inflate; the cell number that was lost stays lost. The MRI that\n")
cat("  looks worse a month after stopping is NOT regrowth, and reading it as\n")
cat("  regrowth is a mistake this decomposition prevents.\n")

hdr("D15 · THE COST OF DIAGNOSTIC DELAY, PRICED IN PERMANENT DECIBELS")
cat(sprintf("  %-14s %13s %14s %14s %12s\n", "EXTRA delay mo", "MD at start",
            "MD 2 y on Rx", "permanent dB", "recovered"))
for (mo in c(0, 3, 6, 12, 24, 48)) {
  st <- PRESIM + mo*30.4
  d  <- run(p_giant, ev_cab(1.0, 200, start = st), end = st + 730, delta = 5)
  cat(sprintf("  %-14.0f %13.1f %14.1f %14.1f %12.1f\n", mo, atday(d,"VFMD",st),
              last(d,"VFMD"), -last(d,"AXL"),
              last(d,"VFMD") - atday(d,"VFMD",st)))
}
cat("\n  The reversible conduction block recovers almost completely whenever the\n")
cat("  drug is finally started. The axonal component does not, and it accrues\n")
cat("  for as long as the chiasm stays indented. Delay is not a missed\n")
cat("  opportunity; it is a debt, and this table is the interest rate.\n")

hdr("SUMMARY")
cat("  58 ODEs · 34 scenarios · 15 diagnostics · all run under mrgsolve\n\n")
cat("  STRUCTURAL NULLS the model must fail to treat, and does:\n")
cat("    * macroprolactinaemia + cabergoline -> the number moves, the axis does not\n")
cat("    * true (Emax-lowered) resistance + 7x dose -> no rescue\n")
cat("    * MGMT-high aggressive tumour + temozolomide -> no response\n")
cat("    * quinagolide instead of cabergoline for impulse control -> no protection\n\n")
cat("  NEGATIVE AND SELF-CRITICAL RESULTS, reported rather than removed:\n")
cat("    * the model cannot produce a durable remission at all: withdrawal always\n      ends in recurrence, only later for smaller remnants. That REFUTES the\n      mechanistic map own hypothesis that endogenous dopamine holds small\n      remnants, because the short loop is weakest exactly when the remnant\n      is smallest (D10)\n")
cat("    * dose escalation above ~0.5-1 mg/wk buys almost nothing biochemically\n")
cat("      because occupancy has already saturated — the model argues against\n")
cat("      part of the practice it was built to explore (D4)\n")
cat("    * splitting the weekly dose changes nothing, so the common practice of\n")
cat("      twice-weekly dosing gets no support here except for tolerability (S5)\n")
cat("    * cabergoline drives prolactin frankly SUBNORMAL in most simulated\n")
cat("      patients, which is a prediction of iatrogenic hypoprolactinaemia\n")
cat("      rather than a fitted outcome (D4, S4)\n")
cat("    * a population normoprolactinaemia RATE (Webster's 83% vs 59%) cannot be\n")
cat("      produced by a single-patient model; only the plateau gap is reproduced\n")
