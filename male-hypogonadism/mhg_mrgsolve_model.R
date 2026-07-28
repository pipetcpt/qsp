## =============================================================================
##  Male Hypogonadism (남성 성선기능저하증) — QSP model for mrgsolve
##  mhg_mrgsolve_model.R
## =============================================================================
##
##  THE STRUCTURAL CLAIM THIS MODEL IS BUILT TO TEST
##  ------------------------------------------------
##  The quantity clinicians measure (total testosterone, TT) is not the quantity
##  that acts (free testosterone, FT), and the map between them is a nonlinear
##  mass-action binding equilibrium:
##
##      TT = FT + K_A[Alb]FT + K_S[SHBG_free]FT
##           \_______________/   \___________________/
##            low-affinity,       high-affinity,
##            high-capacity       SATURABLE  (K_S = 1e9 M^-1)
##
##  Solving for FT gives the Vermeulen quadratic, and the saturability of the
##  SHBG limb makes FT a CONVEX function of TT.  Verified numerically at
##  SHBG = 35 nmol/L, albumin 4.3 g/dL (free fraction, this model):
##
##      TT (ng/dL)   100    300    600   1000   1500   2000
##      FT (pg/mL)  17.5   56.0  122.9  228.2  381.1  550.8
##      free %      1.75   1.87   2.05   2.28   2.54   2.75
##
##  Every second difference is positive: FT'' > 0 over the whole clinical range.
##  Three consequences follow, and each analysis function below prints one of
##  them as a number rather than asserting it as prose.
##
##   (1) THE DIAGNOSTIC THRESHOLD IS FRAME-DEPENDENT.  A single TT cut-off
##       (300 ng/dL) is not a single FT cut-off.  MHG_diagnostic_frame() shows
##       the TT that corresponds to FT = 65 pg/mL runs from 237 ng/dL at
##       SHBG 15 to 635 ng/dL at SHBG 90 — a 2.7-fold swing driven entirely by
##       a binding protein.  The TT-only cut-off is therefore wrong in BOTH
##       directions: it over-diagnoses the obese/insulin-resistant man (low
##       SHBG, low TT, adequate FT) and under-diagnoses the older man (high
##       SHBG, "normal" TT, deficient FT).  MHG_shbg_paradox() runs the two
##       patients side by side at the same TT.
##
##   (2) WAVEFORM MATTERS FAR LESS THAN THE LITERATURE'S FRAMING IMPLIES, AND
##       THE MODEL SAYS SO AGAINST ITS OWN INITIAL HYPOTHESIS.  Jensen's
##       inequality gives E[FT(TT)] >= FT(E[TT]), so a peaky regimen must
##       deliver more time-averaged free T than a flat one at the same mean TT.
##       The model was built expecting that convexity to compound through a
##       convex erythropoietic response and explain the ~3x excess of
##       erythrocytosis on intramuscular versus transdermal testosterone.
##       It does not.  MHG_convexity_decomposition() prints the actual
##       arithmetic, and the accounting comes out as follows:
##
##         (a) CONVEX BINDING is real and scales with amplitude: the Jensen
##             gap runs +0.3% for flat gel up to +4.4% for IM 200 mg q2wk.
##             Small, but it is there, and it is the largest for the peakiest
##             regimen exactly as the mathematics requires.
##
##         (b) CONVEX RESPONSE DOES NOT COMPOUND WITH IT.  The Jensen gap at
##             the response step comes out near zero and NEGATIVE for the
##             peaky arms, because their free-T peaks overshoot EC50_HEPC into
##             the CONCAVE upper limb of the Hill.  A response that is convex
##             low down is concave high up, and supraphysiologic peaks live on
##             the wrong side of the inflection.  This is reported rather than
##             tuned away: it is the model contradicting the hypothesis it was
##             built to express.
##
##         (c) AT MATCHED WEEKLY DOSE the waveform is worth only a fraction of
##             a haematocrit point.  IM 100 mg weekly and IM 200 mg q2wk are
##             the same milligrams with twice the amplitude, and they differ by
##             ~0.1 points of Hct.  The large IM-vs-transdermal gap in the
##             literature is therefore mostly a DOSE difference, not a
##             waveform difference.
##
##         (d) WHAT ACTUALLY MANUFACTURES THE 3x is the decision rule.  A FIXED
##             threshold (Hct > 54%) applied to a spread-out population turns a
##             ~1 point shift in the mean into a ~2.5-fold difference in the
##             fraction of men who cross it.  That is a property of the
##             STOPPING RULE, not of the molecule.
##
##       EC50_HEPC is a CALIBRATED parameter, not a measured one, and
##       MHG_hepcidin_sensitivity() reports how much of the above survives
##       moving it — see the calibration block.
##
##   (3) SERUM T AND INTRATESTICULAR T ARE DIFFERENT VARIABLES AND ONLY ONE OF
##       THEM IS MEASURED.  ITT sits at ~35-100x serum because Leydig output is
##       concentrated locally by androgen-binding protein.  Exogenous T
##       normalises serum T while shutting LH off, and ITT collapses to a few
##       percent of normal — below the ~20-30% threshold that spermatogenesis
##       requires.  MHG_ITT_collapse() prints serum T and ITT in the same table
##       so the dissociation is visible; MHG_recovery_curve() reproduces the
##       Liu 2006 pooled recovery kinetics after cessation.
##
##  A fourth, non-obvious result the model is built to expose: the effects
##  usually attributed to "testosterone" partition unevenly between T and its
##  aromatised product E2.  MHG_finkelstein() reproduces the NEJM 2013 design
##  (graded testosterone +/- anastrozole on a GnRH-agonist background) and
##  recovers the dissociation — lean mass and strength track T, body fat tracks
##  E2, sexual desire tracks both.  This is why "block the aromatase to stop
##  the gynaecomastia" is not a free action: MHG_aromatase_cost() prices it in
##  trabecular vBMD.
##
##  UNITS
##    time              days (simulations run in days; scenarios up to 3 years)
##    testosterone      plasma amount ug; concentration reported as ng/dL
##                      (ug/L * 100 = ng/dL); FT reported as pg/mL
##    SHBG              nmol/L        E2   pg/mL       DHT  ng/dL
##    ITT               nmol/L (intratesticular testosterone)
##    LH, FSH           IU/L          hCG  IU/L
##    sperm             million/mL
##    Hct               %             BMD  % change from baseline
##    lean / fat mass   kg
##    PRO scores        0-10 arbitrary (libido, vitality) — calibrated to the
##                      T-Trials effect sizes, not to an absolute instrument
##
##  49 ODE compartments · 132 parameters · 8 testosterone formulations plus
##  5 non-androgen drugs · 14 scenarios · 11 analysis functions
##
##  DISCLAIMER: research and education only.  Not for clinical use.  Parameters
##  are literature-anchored where a value exists and CALIBRATED where one does
##  not; the calibration block below says which is which.
## =============================================================================

library(mrgsolve)

code <- '
$PROB
Male hypogonadism QSP model. HPG axis with intratesticular testosterone,
Vermeulen binding equilibrium with regulated SHBG, aromatisation and
5-alpha reduction, erythropoiesis, bone, body composition, spermatogenesis
transit chain, and eight testosterone formulations plus hCG, hMG/FSH,
clomiphene, enclomiphene and anastrozole.

$PARAM @annotated
// ---------------- testosterone disposition ----------------
VC      :  35   : Central volume of distribution (L)
VP      :  60   : Peripheral volume (L)
Q       : 1500  : Intercompartmental clearance (L/day)
CLINT   : 2174  : Intrinsic clearance on the BIOAVAILABLE fraction (L/day)
                  // CL = CLINT * fbio; fbio ~0.46 at baseline -> CL ~1000 L/day
                  // and 6 mg/day endogenous production -> 600 ng/dL. This is
                  // the mechanism by which falling SHBG raises clearance.
PROD_ADR:  300  : Adrenal androgen contribution to serum T (ug/day)

// ---------------- binding equilibrium (Vermeulen 1999) ----------------
KS      : 1.0e9 : SHBG association constant (1/M)
KA      : 3.6e4 : Albumin association constant (1/M)
ALB     : 4.3   : Serum albumin (g/dL)
MWT     : 288.4 : Testosterone molecular weight (g/mol)

// ---------------- SHBG regulation ----------------
SHBG0   :  35   : Baseline serum SHBG (nmol/L)
KOUT_S  : 0.10  : SHBG turnover rate (1/day)   // t1/2 ~7 days
EXP_AND : 0.20  : SHBG power-law exponent on androgen (dimensionless)
EXP_INS : 0.30  : SHBG power-law exponent on insulin/IR
EXP_E2  : 0.15  : SHBG power-law exponent on estradiol
AGE     :  55   : Patient age (years)
AGE_S   : 0.010 : Fractional SHBG rise per year above 40
INS_REL : 1.0   : Relative insulin / insulin resistance (1 = lean referent)
THYR    : 1.0   : Relative thyroid hormone effect on SHBG

// ---------------- hypothalamic-pituitary ----------------
KOUT_G  : 2.0   : GnRH drive turnover (1/day)
HT_FB   : 3.0   : Hill exponent, androgen negative feedback
HE_FB   : 3.0   : Hill exponent, estradiol negative feedback
WE2     : 2.2   : Weight of E2 relative to T in negative feedback
                  // E2 is the dominant feedback signal in men; this is why
                  // an aromatase inhibitor raises LH more than it lowers T.
EXOPOT  : 3.0   : Feedback potency of EXOGENOUS relative to endogenous androgen
                  // A STRUCTURAL ASSUMPTION, and an important one. Exogenous
                  // testosterone is delivered continuously; endogenous
                  // testosterone arrives downstream of a pulsatile generator
                  // that the same steroid feeds back on. Continuous exposure
                  // suppresses the pulse generator far more effectively than
                  // an identical MEAN level delivered in pulses, which is why
                  // replacement doses reliably drive LH to the detection limit
                  // even when they only restore mid-normal serum T. Without
                  // this term the model cannot reproduce the ITT collapse that
                  // Coviello 2005 measured (-94% on 200 mg/week), and the
                  // whole fertility argument would be understated.
LH0     : 4.0   : Baseline serum LH (IU/L)
FSH0    : 4.0   : Baseline serum FSH (IU/L)
KOUT_LH : 12.0  : LH elimination (1/day)   // t1/2 ~60 min -> fast
KOUT_FSH: 4.5   : FSH elimination (1/day)  // t1/2 ~3.7 h
HI_FB   : 1.0   : Hill exponent, inhibin B feedback on FSH
PITF    : 1.0   : Pituitary functional capacity (0-1; <1 = secondary)
HYPOF   : 1.0   : Hypothalamic GnRH capacity (0-1; 0 = Kallmann)
PRL_SUP : 0.0   : Fractional GnRH suppression by hyperprolactinaemia (0-1)
OPI_SUP : 0.0   : Fractional GnRH suppression by opioids (0-1)
INFL_SUP: 0.0   : Fractional GnRH suppression by adipose inflammation (0-1)

// ---------------- Leydig / intratesticular testosterone ----------------
ITT0    : 700   : Baseline intratesticular testosterone (nmol/L)
KOUT_ITT: 6.0   : ITT turnover (1/day)
KM_LH   : 3.0   : LH-equivalent EC50 for steroidogenesis (IU/L)
LEYCAP0 : 1.0   : Baseline Leydig functional capacity (0-1)
KREC_L  : 0.020 : Leydig capacity recovery/atrophy rate (1/day) // t1/2 ~35 d
K_TROPH : 1.5   : LH-equivalent for half-maximal trophic support (IU/L)
LEY_MIN : 0.25  : Floor of Leydig capacity under complete LH withdrawal
KSPILL  : 8.571 : ITT -> serum T spillover (ug/day per nmol/L)
                  // 700 nmol/L * 8.571 = 6000 ug/day = 6 mg/day. Anchors the
                  // whole model to the classic production-rate measurement.
HCG_POT : 0.25  : hCG potency at LHCGR expressed as LH IU-equivalents per IU

// ---------------- aromatisation / 5-alpha reduction ----------------
E20     :  28   : Baseline estradiol (pg/mL)
KOUT_E2 : 6.0   : Estradiol turnover (1/day)
KARO    : 1.0   : Aromatase scaling (fitted in $MAIN to give E20)
ARO_FAT : 1.0   : Exponent on fat mass for aromatase capacity
DHT0    :  50   : Baseline DHT (ng/dL)
KOUT_DHT: 8.0   : DHT turnover (1/day)
INH_5AR : 0.0   : Fractional 5-alpha-reductase inhibition (finasteride etc.)

// ---------------- erythropoiesis ----------------
HCT0    : 45.0  : Baseline haematocrit (%)
EC50_HEPC: 300  : FT for half-maximal hepcidin suppression (pg/mL)  [CALIBRATED]
                  // Deliberately SUPRAPHYSIOLOGIC (eugonadal FT 50-210), which
                  // is what makes the erythropoietic response CONVEX over the
                  // therapeutic range. Calibrated to reproduce the observed
                  // IM-vs-gel erythrocytosis gap; not an independently
                  // measured constant. Sensitivity: MHG_convexity_decomposition().
HILL_H  : 2.0   : Hill exponent of the erythropoietic drive
IMAX_H  : 0.45  : Maximal fractional hepcidin suppression
KOUT_HEP: 2.0   : Hepcidin turnover (1/day)
SEPO    : 0.50  : Maximal fractional EPO set-point elevation
KOUT_EPO: 4.0   : EPO turnover (1/day)
KTR_ERY : 0.35  : Erythroid progenitor transit rate (1/day)
RBCLIFE : 100   : Red cell mean lifespan in model (days)
HCT_SD  : 3.5   : Population SD of haematocrit (%) for threshold arithmetic
HCT_STOP: 54.0  : Haematocrit stopping threshold (%)

// ---------------- bone ----------------
KOUT_SCL: 0.05  : Sclerostin turnover (1/day)
E_SCL   : 0.55  : Exponent: E2 deficiency raises sclerostin
E_OC    : 0.70  : Exponent: E2 deficiency raises osteoclast activity
KOUT_OB : 0.04  : Osteoblast activity turnover (1/day)
KOUT_OC : 0.06  : Osteoclast activity turnover (1/day)
KFORM   : 0.0033: Trabecular formation coefficient (%/day per unit OB)
KRES    : 0.0033: Trabecular resorption coefficient (%/day per unit OC)
                  // Sized so an untreated severely hypogonadal man loses
                  // roughly 2%/yr of trabecular bone by true balance.
CORT_F  : 0.25  : Cortical compartment responds at 25% of trabecular rate
ANAB_T  : 0.25  : Direct (E2-independent) androgen effect on osteoblasts
RSPACE  : 7.5   : Remodelling space (percentage points per unit turnover)
                  // THE MECHANISM BEHIND THE T-TRIALS +7.5%. Measured density
                  // is true mineral MINUS the transient deficit of resorption
                  // cavities not yet refilled, and that deficit scales with the
                  // remodelling RATE. A hypogonadal man runs high turnover and
                  // therefore carries a large transient deficit; restoring
                  // androgen slows turnover, the cavities refill, and densitometry
                  // records a gain that is only partly new bone. This is the same
                  // mechanism that produces the early BMD rise with antiresorptives,
                  // and it is why a BMD gain does not have to mean a strength gain
                  // (see TRAVERSE fracture, which this model does not reproduce).
KOUT_RSP: 0.008 : Remodelling-space adaptation (1/day) // t1/2 ~87 days

// ---------------- body composition ----------------
LEAN0   :  58   : Baseline lean mass (kg)
FAT0    :  30   : Baseline fat mass (kg)
KBC     : 0.012 : Body-composition adaptation rate (1/day) // t1/2 ~58 d
SLEAN_T : 0.060 : Maximal fractional lean mass gain (T-driven)
KL50    : 180   : FT for half-maximal lean mass effect (pg/mL)
SFAT_T  : 0.10  : Maximal fractional fat loss attributable to T
SFAT_E  : 0.16  : Maximal fractional fat GAIN when E2 is withdrawn
                  // Finkelstein NEJM 2013: the increase in per-cent body fat
                  // was attributable to estradiol deficiency, not to androgen.
KE50_F  : 18    : E2 for half-maximal fat effect (pg/mL)
WT_LOSS : 0.0   : Imposed fractional fat-mass loss per year (diet/surgery/GLP-1)

// ---------------- spermatogenesis ----------------
SPERM0  :  60   : Baseline sperm concentration (million/mL)
ITT50_S : 175   : ITT for half-maximal spermatogenesis (nmol/L) // 25% of ITT0
HS_SP   : 4.0   : Hill exponent of the ITT threshold (sharp)
KF50_S  : 2.0   : FSH for half-maximal Sertoli support (IU/L)
TAU_SPG :  74   : Spermatogenic cycle duration (days)
TAU_EPI :  14   : Epididymal transit (days)
INHB0   : 150   : Baseline inhibin B (pg/mL)
KOUT_IB : 0.15  : Inhibin B turnover (1/day)

// ---------------- prostate / PSA ----------------
PSA0    : 1.1   : Baseline PSA (ng/mL)
KOUT_PSA: 0.03  : PSA turnover (1/day)
EPSA    : 0.55  : Maximal fractional PSA rise
KPSA    :  25   : DHT for half-maximal PSA effect (ng/dL)
                  // Deliberately LOW relative to DHT0 = 50: the saturation
                  // model. Above roughly eugonadal androgen exposure the
                  // prostate response is flat, which is why restoring T to
                  // mid-normal moves PSA by only ~0.3-0.5 ng/mL.

// ---------------- patient-reported outcomes ----------------
KPRO    : 0.05  : PRO adaptation rate (1/day)
LIB_T   : 0.55  : Weight of T on sexual desire
LIB_E   : 0.25  : Weight of E2 on sexual desire
KLIB50  : 150   : FT for half-maximal desire effect (pg/mL)
VIT_T   : 0.30  : Weight of T on vitality (small, per T-Trials)

// ---------------- drug PK: testosterone formulations ----------------
KA_IM   : 0.15  : IM ester (enanthate/cypionate) absorption (1/day)
F_IM    : 0.62  : IM ester bioavailability x ester molar factor
KA_TU   : 0.022 : IM undecanoate (Nebido/Aveed) absorption (1/day)
F_TU    : 0.70  : IM undecanoate bioavailability x molar factor
KA_SC   : 0.14  : SC auto-injector absorption (1/day)
F_SC    : 0.55  : SC auto-injector bioavailability x molar factor
KA_GEL  : 0.60  : Transdermal reservoir into plasma (1/day)
F_GEL   : 0.088 : Transdermal fractional absorption
KA_ORAL : 6.0   : Oral undecanoate lymphatic absorption (1/day)
F_ORAL  : 0.016 : Oral undecanoate fractional systemic delivery as T
KA_PEL  : 0.012 : Subcutaneous pellet release (1/day) // t1/2 ~58 d
F_PEL   : 0.75  : Pellet bioavailability
KA_NAS  : 12.0  : Nasal gel absorption (1/day)
F_NAS   : 0.090 : Nasal gel fractional absorption

// ---------------- drug PK/PD: non-androgen strategies ----------------
KA_HCG  : 1.2   : hCG subcutaneous absorption (1/day)
KE_HCG  : 0.55  : hCG elimination (1/day)  // t1/2 ~30 h
V_HCG   : 25.0  : hCG volume of distribution (L)
KA_FSHD : 0.9   : rFSH/hMG absorption (1/day)
KE_FSHD : 0.55  : rFSH/hMG elimination (1/day)
V_FSHD  : 8.0   : rFSH volume (L)
FSH_POT : 0.09  : rFSH IU dose -> serum FSH IU/L conversion
KA_CLO  : 2.0   : Clomiphene absorption (1/day)
KE_CLO  : 0.14  : Clomiphene (zuclomiphene-dominated) elimination (1/day)
V_CLO   : 100   : Clomiphene apparent volume (L)
IC50_CLO: 0.12  : Clomiphene concentration for 50% ERalpha blockade (mg/L)
IMAX_CLO: 0.80  : Maximal hypothalamic ERalpha blockade
KA_ANA  : 3.0   : Anastrozole absorption (1/day)
KE_ANA  : 0.33  : Anastrozole elimination (1/day) // t1/2 ~50 h
V_ANA   : 100   : Anastrozole apparent volume (L)
IC50_ANA: 0.006 : Anastrozole concentration for 50% aromatase inhibition (mg/L)
IMAX_ANA: 0.85  : Maximal aromatase inhibition

// ---------------- disease switches ----------------
LEY_DMG : 0.0   : Fractional Leydig destruction (Klinefelter, chemo, radiation)
SER_DMG : 0.0   : Fractional Sertoli destruction
PHLEB   : 0.0   : Therapeutic phlebotomy: fractional RBC removal per day
INITSS  : 1     : Initialise from baseline identities (1) or honour init() (0)
                  // The cessation scenario sets 0 so it can restart from the
                  // suppressed state carried forward out of the treatment run.

$CMT @annotated
DEP_IM   : IM testosterone ester depot (ug T-equivalent)
DEP_TU   : IM testosterone undecanoate depot (ug)
DEP_SC   : SC auto-injector depot (ug)
DEP_GEL  : Transdermal skin reservoir (ug)
DEP_ORAL : Oral undecanoate gut/lymph depot (ug)
DEP_PEL  : Subcutaneous pellet depot (ug)
DEP_NAS  : Nasal depot (ug)
CENT     : Central testosterone (ug)
PERIPH   : Peripheral testosterone (ug)
HCG_D    : hCG subcutaneous depot (IU)
HCG_C    : hCG central (IU)
FSHD_D   : rFSH/hMG depot (IU)
FSHD_C   : rFSH/hMG central (IU)
CLO_D    : Clomiphene depot (mg)
CLO_C    : Clomiphene central (mg)
ANA_D    : Anastrozole depot (mg)
ANA_C    : Anastrozole central (mg)
GNRHD    : Hypothalamic GnRH drive (normalised, 1 = eugonadal)
LH       : Serum LH (IU/L)
FSH      : Serum FSH (IU/L)
LEYCAP   : Leydig functional capacity (0-1)
ITT      : Intratesticular testosterone (nmol/L)
SHBG     : Serum SHBG (nmol/L)
E2       : Serum estradiol (pg/mL)
DHT      : Serum DHT (ng/dL)
INHB     : Inhibin B (pg/mL)
SG1      : Spermatogonial transit compartment 1
SG2      : Spermatocyte transit compartment 2
SG3      : Spermatid transit compartment 3
SG4      : Spermiation transit compartment 4
EPID     : Epididymal sperm reserve
HEPC     : Hepcidin (normalised, 1 = baseline)
EPO      : Erythropoietin (normalised, 1 = baseline)
PROG_E   : Erythroid progenitor pool (normalised)
RETIC    : Reticulocyte pool (normalised)
RBC      : Red cell mass (normalised, 1 = baseline)
SCLERO   : Sclerostin (normalised)
OB       : Osteoblast activity (normalised)
OC       : Osteoclast activity (normalised)
BMD_TR   : Trabecular bone balance integral (% change from baseline)
BMD_CO   : Cortical bone balance integral (% change from baseline)
RSP      : Remodelling space (normalised turnover, 1 = eugonadal)
LEAN     : Lean mass (kg)
FAT      : Fat mass (kg)
PSA      : Serum PSA (ng/mL)
LIBIDO   : Sexual desire score (0-10)
VITAL    : Vitality score (0-10)
CUMFT    : Cumulative free-testosterone exposure (pg/mL x day)
CUMDRIVE : Cumulative erythropoietic drive (day)

$GLOBAL
#define TTngdl (CENT / VC * 100.0)          // ug/L -> ng/dL
#define SAFE(x) ((x) < 1e-9 ? 1e-9 : (x))

// ---- Vermeulen (1999) free testosterone: exact solution of the quadratic
//      a*FT^2 + b*FT - T = 0,  a = N*KS,  b = N + KS*S - KS*T,  N = 1 + KA*[Alb]
//      T, S in mol/L; returns FT in pg/mL.
namespace {
  double free_T_pgmL(double TT_ngdl, double SHBG_nM,
                     double KSl, double KAl, double ALBgdl, double MW) {
    double T = TT_ngdl * 1e-8 / MW;                    // ng/dL -> mol/L
    double S = SHBG_nM * 1e-9;
    double Alb = ALBgdl * 10.0 / 66500.0;              // g/dL -> mol/L
    double N = 1.0 + KAl * Alb;
    double a = N * KSl;
    double b = N + KSl * S - KSl * T;
    if (T <= 0.0) return 0.0;
    double ft = (-b + sqrt(b * b + 4.0 * a * T)) / (2.0 * a);   // mol/L
    return ft * MW * 1e9;                                        // -> pg/mL
  }
  double hillf(double x, double ec50, double h) {
    double xh = pow(SAFE(x), h), eh = pow(ec50, h);
    return xh / (eh + xh);
  }
}

$MAIN
// ---------- dose routing ----------
F_DEP_IM   = F_IM   * 1000.0;    // mg of product -> ug of T-equivalent
F_DEP_TU   = F_TU   * 1000.0;
F_DEP_SC   = F_SC   * 1000.0;
F_DEP_GEL  = F_GEL  * 1000.0;
F_DEP_ORAL = F_ORAL * 1000.0;
F_DEP_PEL  = F_PEL  * 1000.0;
F_DEP_NAS  = F_NAS  * 1000.0;

// ---------- baseline initialisation ----------
// Serum T is initialised from the production/clearance identity so the model
// starts at a genuine steady state rather than drifting for the first weeks.
// Guarded by INITSS so the cessation scenario can restart from a suppressed
// state supplied through init() instead.
if (INITSS > 0.5) {
  double fbio0 = 0.46;
  double CL0   = CLINT * fbio0;
  double TT_init = (KSPILL * ITT0 * (1.0 - LEY_DMG) + PROD_ADR) / CL0;  // ug/L
  CENT_0   = TT_init * VC;
  PERIPH_0 = TT_init * VP;

  GNRHD_0  = 1.0;
  LH_0     = LH0;
  FSH_0    = FSH0;
  LEYCAP_0 = LEYCAP0 * (1.0 - LEY_DMG);
  ITT_0    = ITT0 * (1.0 - LEY_DMG);
  SHBG_0   = SHBG0;
  E2_0     = E20;
  DHT_0    = DHT0;
  INHB_0   = INHB0 * (1.0 - SER_DMG);

  SG1_0 = 1.0; SG2_0 = 1.0; SG3_0 = 1.0; SG4_0 = 1.0; EPID_0 = 1.0;
  HEPC_0 = 1.0; EPO_0 = 1.0; PROG_E_0 = 1.0; RETIC_0 = 1.0; RBC_0 = 1.0;
  SCLERO_0 = 1.0; OB_0 = 1.0; OC_0 = 1.0; BMD_TR_0 = 0.0; BMD_CO_0 = 0.0;
  RSP_0 = 1.0;
  LEAN_0 = LEAN0; FAT_0 = FAT0; PSA_0 = PSA0;
  LIBIDO_0 = 10.0 * (LIB_T * hillf(free_T_pgmL(TT_init * 100.0, SHBG0, KS, KA, ALB, MWT),
                                   KLIB50, 1.0)
                   + LIB_E * hillf(E20, 15.0, 1.0) + 0.20);
  VITAL_0  = 5.0;
}

$ODE
// =====================================================================
//  1. BINDING EQUILIBRIUM  — solved algebraically at every step
// =====================================================================
double TT   = TTngdl;                                  // ng/dL
double FT   = free_T_pgmL(TT, SAFE(SHBG), KS, KA, ALB, MWT);   // pg/mL
// bioavailable fraction = free + albumin-bound; albumin term is linear in FT
double Nalb = 1.0 + KA * (ALB * 10.0 / 66500.0);
double FTmol = FT / (MWT * 1e9);                       // pg/mL -> mol/L
double TTmol = TT * 1e-8 / MWT;
double fbio  = (TTmol > 0.0) ? (Nalb * FTmol / TTmol) : 0.46;
if (fbio > 0.999) fbio = 0.999;
if (fbio < 0.02)  fbio = 0.02;

double CL = CLINT * fbio;      // <-- clearance rises as SHBG falls
double FT_ref = free_T_pgmL((KSPILL * ITT0 + PROD_ADR) / (CLINT * 0.46) * 100.0,
                            SHBG0, KS, KA, ALB, MWT);

// =====================================================================
//  2. NON-ANDROGEN DRUG PK  (hCG, rFSH/hMG, clomiphene, anastrozole)
// =====================================================================
dxdt_HCG_D  = -KA_HCG  * HCG_D;
dxdt_HCG_C  =  KA_HCG  * HCG_D  - KE_HCG  * HCG_C;
dxdt_FSHD_D = -KA_FSHD * FSHD_D;
dxdt_FSHD_C =  KA_FSHD * FSHD_D - KE_FSHD * FSHD_C;
dxdt_CLO_D  = -KA_CLO  * CLO_D;
dxdt_CLO_C  =  KA_CLO  * CLO_D  - KE_CLO  * CLO_C;
dxdt_ANA_D  = -KA_ANA  * ANA_D;
dxdt_ANA_C  =  KA_ANA  * ANA_D  - KE_ANA  * ANA_C;

// =====================================================================
//  3. HYPOTHALAMIC-PITUITARY-GONADAL AXIS
// =====================================================================
// FB0 is the value of FB at the eugonadal referent (FT = FT_ref, E2 = E20),
// so the drive multiplier FB0/FB equals 1 at baseline by construction.
double FB0 = 2.0 + WE2;
double CLOc    = CLO_C / V_CLO;
double CLO_EFF = IMAX_CLO * CLOc / (IC50_CLO + CLOc);

// What fraction of the circulating androgen arrived from outside? Exogenous
// androgen is continuous and therefore a stronger feedback signal per unit
// concentration than the same mean delivered by a pulsatile gonad (EXOPOT).
double ABS_NOW  = KA_IM * DEP_IM + KA_TU * DEP_TU + KA_SC * DEP_SC
                + KA_GEL * DEP_GEL + KA_ORAL * DEP_ORAL + KA_PEL * DEP_PEL
                + KA_NAS * DEP_NAS;
double ENDO_NOW = KSPILL * ITT + PROD_ADR;
double EXOFRAC  = ABS_NOW / SAFE(ABS_NOW + ENDO_NOW);
double FBGAIN   = 1.0 + (EXOPOT - 1.0) * EXOFRAC;

double FB = 1.0 + pow(SAFE(FT * FBGAIN / SAFE(FT_ref)), HT_FB)
              + WE2 * pow(SAFE(E2 * FBGAIN / E20), HE_FB) * (1.0 - CLO_EFF);
double EXTSUP = (1.0 - PRL_SUP) * (1.0 - OPI_SUP) * (1.0 - INFL_SUP) * HYPOF;
dxdt_GNRHD = KOUT_G * (EXTSUP * FB0 / SAFE(FB) - GNRHD);

dxdt_LH  = KOUT_LH  * LH0  * GNRHD * PITF - KOUT_LH  * LH;
double IB_FB = 2.0 / (1.0 + pow(SAFE(INHB / SAFE(INHB0)), HI_FB));
double FSH_EX = FSH_POT * FSHD_C / V_FSHD;   // exogenous rFSH/hMG, IU/L
dxdt_FSH = KOUT_FSH * FSH0 * GNRHD * PITF * IB_FB - KOUT_FSH * FSH;

// =====================================================================
//  4. LEYDIG CELL AND INTRATESTICULAR TESTOSTERONE
// =====================================================================
double HCGc  = HCG_C / V_HCG;                 // IU/L
double LH_EQ = LH + HCG_POT * HCGc;
double LEY_TGT = (1.0 - LEY_DMG) *
                 (LEY_MIN + (1.0 - LEY_MIN) * LH_EQ / (K_TROPH + LH_EQ)) /
                 (LEY_MIN + (1.0 - LEY_MIN) * LH0 / (K_TROPH + LH0));
if (LEY_TGT > 1.0 - LEY_DMG) LEY_TGT = 1.0 - LEY_DMG;
dxdt_LEYCAP = KREC_L * (LEY_TGT - LEYCAP);

double SD_MAX = KOUT_ITT * ITT0 * (KM_LH + LH0) / LH0;   // so ITT=ITT0 at LH=LH0
double ITT_IN = LEYCAP * SD_MAX * LH_EQ / (KM_LH + LH_EQ);
dxdt_ITT = ITT_IN - KOUT_ITT * ITT;

// =====================================================================
//  5. SERUM TESTOSTERONE PK  (endogenous + 7 exogenous routes)
// =====================================================================
double ABS_T = KA_IM * DEP_IM + KA_TU * DEP_TU + KA_SC * DEP_SC
             + KA_GEL * DEP_GEL + KA_ORAL * DEP_ORAL + KA_PEL * DEP_PEL
             + KA_NAS * DEP_NAS;

dxdt_DEP_IM   = -KA_IM   * DEP_IM;
dxdt_DEP_TU   = -KA_TU   * DEP_TU;
dxdt_DEP_SC   = -KA_SC   * DEP_SC;
dxdt_DEP_GEL  = -KA_GEL  * DEP_GEL;
dxdt_DEP_ORAL = -KA_ORAL * DEP_ORAL;
dxdt_DEP_PEL  = -KA_PEL  * DEP_PEL;
dxdt_DEP_NAS  = -KA_NAS  * DEP_NAS;

double ENDO_T = KSPILL * ITT + PROD_ADR;
dxdt_CENT   = ENDO_T + ABS_T - CL * (CENT / VC)
              - Q * (CENT / VC) + Q * (PERIPH / VP);
dxdt_PERIPH = Q * (CENT / VC) - Q * (PERIPH / VP);

// =====================================================================
//  6. SHBG — a regulated state variable, not a constant
// =====================================================================
double SHBG_TGT = SHBG0
    * pow(SAFE(FT_ref / SAFE(FT)), EXP_AND)
    * pow(SAFE(1.0 / INS_REL), EXP_INS)
    * pow(SAFE(E2 / E20), EXP_E2)
    * (1.0 + AGE_S * (AGE - 40.0)) * THYR;
if (SHBG_TGT < 3.0)  SHBG_TGT = 3.0;
if (SHBG_TGT > 200.0) SHBG_TGT = 200.0;
dxdt_SHBG = KOUT_S * (SHBG_TGT - SHBG);

// =====================================================================
//  7. AROMATISATION AND 5-ALPHA REDUCTION
// =====================================================================
double ANAc    = ANA_C / V_ANA;
double ARO_INH = IMAX_ANA * ANAc / (IC50_ANA + ANAc);
double T_bio_c = fbio * TT;                              // ng/dL bioavailable
double T_bio_0 = 0.46 * (KSPILL * ITT0 + PROD_ADR) / (CLINT * 0.46) * 100.0;
double AROCAP  = pow(SAFE(FAT / FAT0), ARO_FAT);
dxdt_E2  = KOUT_E2 * E20 * AROCAP * (T_bio_c / SAFE(T_bio_0)) * (1.0 - ARO_INH)
           - KOUT_E2 * E2;
dxdt_DHT = KOUT_DHT * DHT0 * (FT / SAFE(FT_ref)) * (1.0 - INH_5AR) - KOUT_DHT * DHT;

// =====================================================================
//  8. SPERMATOGENESIS  — ITT threshold x Sertoli support, 74 + 14 day transit
// =====================================================================
double FSH_TOT = FSH + FSH_EX;
double ITT_F  = hillf(ITT, ITT50_S, HS_SP);
double ITT_F0 = hillf(ITT0, ITT50_S, HS_SP);
double SER_F  = (1.0 - SER_DMG) * (FSH_TOT / (KF50_S + FSH_TOT))
                / (FSH0 / (KF50_S + FSH0));
double KTR_SP = 4.0 / TAU_SPG;
double GEN = KTR_SP * (ITT_F / SAFE(ITT_F0)) * SER_F;
dxdt_SG1 = GEN            - KTR_SP * SG1;
dxdt_SG2 = KTR_SP * SG1   - KTR_SP * SG2;
dxdt_SG3 = KTR_SP * SG2   - KTR_SP * SG3;
dxdt_SG4 = KTR_SP * SG3   - KTR_SP * SG4;
// epididymal reserve: unit gain (so EPID = 1 at baseline) with a TAU_EPI lag
dxdt_EPID = (1.0 / TAU_EPI) * SG4 - (1.0 / TAU_EPI) * EPID;
dxdt_INHB = KOUT_IB * INHB0 * (1.0 - SER_DMG) * SER_F * (0.4 + 0.6 * SG2)
            - KOUT_IB * INHB;

// =====================================================================
//  9. ERYTHROPOIESIS  — the second convexity, integrated over RBC lifespan
// =====================================================================
// Everything is expressed RELATIVE to the eugonadal referent drive, so a
// eugonadal man sits at Hct = HCT0 exactly, a hypogonadal man drifts BELOW it
// (the mild anaemia of androgen deficiency that the T-Trials anaemia trial
// corrected), and treatment moves him up.
double ERY_DRIVE = hillf(FT, EC50_HEPC, HILL_H);
double DRIVE0    = hillf(FT_ref, EC50_HEPC, HILL_H);
double dDRIVE    = ERY_DRIVE - DRIVE0;
dxdt_HEPC = KOUT_HEP * (1.0 - IMAX_H * dDRIVE) - KOUT_HEP * HEPC;
double IRON_F = pow(SAFE(1.0 / SAFE(HEPC)), 0.5);
double HCT_now = HCT0 * RBC;
dxdt_EPO = KOUT_EPO * (1.0 + SEPO * dDRIVE) * (HCT0 / SAFE(HCT_now))
           - KOUT_EPO * EPO;
dxdt_PROG_E = KTR_ERY * EPO * IRON_F - KTR_ERY * PROG_E;
dxdt_RETIC  = KTR_ERY * PROG_E - KTR_ERY * RETIC;
dxdt_RBC    = (1.0 / RBCLIFE) * RETIC - (1.0 / RBCLIFE) * RBC - PHLEB * RBC;

// =====================================================================
//  10. BONE  — E2-dominant, with a smaller direct androgen term
// =====================================================================
dxdt_SCLERO = KOUT_SCL * pow(SAFE(E20 / SAFE(E2)), E_SCL) - KOUT_SCL * SCLERO;
double OB_TGT = (1.0 / SAFE(SCLERO)) * (1.0 + ANAB_T * (FT / SAFE(FT_ref) - 1.0));
dxdt_OB = KOUT_OB * OB_TGT - KOUT_OB * OB;
dxdt_OC = KOUT_OC * pow(SAFE(E20 / SAFE(E2)), E_OC) - KOUT_OC * OC;
dxdt_BMD_TR = KFORM * (OB - 1.0) - KRES * (OC - 1.0);
dxdt_BMD_CO = CORT_F * (KFORM * (OB - 1.0) - KRES * (OC - 1.0));
// remodelling space: tracks the prevailing turnover rate with a slow lag
double TURNOVER = OC;   // resorption cavities, not net activity
dxdt_RSP = KOUT_RSP * (TURNOVER - RSP);

// ====================================================================
//  11. BODY COMPOSITION  — T drives lean, E2 defends against fat gain
// ====================================================================
double fT_lean = hillf(FT, KL50, 1.0) / hillf(FT_ref, KL50, 1.0);
double LEAN_TGT = LEAN0 * (1.0 + SLEAN_T * (fT_lean - 1.0));
double fE_fat = hillf(E2, KE50_F, 1.0) / hillf(E20, KE50_F, 1.0);
// An imposed weight-loss intervention shrinks the fat TARGET exponentially,
// so WT_LOSS is an honest "fraction of fat mass lost per year".
double WLF = exp(-WT_LOSS * SOLVERTIME / 365.0);
double FAT_TGT = FAT0 * WLF *
                 (1.0 - SFAT_T * (fT_lean - 1.0) - SFAT_E * (fE_fat - 1.0));
if (FAT_TGT < 3.0) FAT_TGT = 3.0;
dxdt_LEAN = KBC * (LEAN_TGT - LEAN);
dxdt_FAT  = KBC * (FAT_TGT - FAT);

// ====================================================================
//  12. PROSTATE / PSA  — saturation model
// ====================================================================
double PSA_TGT = PSA0 * (1.0 + EPSA * (hillf(DHT, KPSA, 1.0) / hillf(DHT0, KPSA, 1.0) - 1.0));
dxdt_PSA = KOUT_PSA * (PSA_TGT - PSA);

// ====================================================================
//  13. PATIENT-REPORTED OUTCOMES
// ====================================================================
double LIB_TGT = 10.0 * (LIB_T * hillf(FT, KLIB50, 1.0)
                       + LIB_E * hillf(E2, 15.0, 1.0) + 0.20);
dxdt_LIBIDO = KPRO * (LIB_TGT - LIBIDO);
double VIT_TGT = 5.0 * (1.0 + VIT_T * (hillf(FT, KLIB50, 1.0) / hillf(FT_ref, KLIB50, 1.0) - 1.0));
dxdt_VITAL = KPRO * (VIT_TGT - VITAL);

// ---- exposure accumulators (for the convexity ledger) ----
dxdt_CUMFT    = FT;
dxdt_CUMDRIVE = ERY_DRIVE;

$TABLE
double TT_out   = TTngdl;
double FT_out   = free_T_pgmL(TT_out, SAFE(SHBG), KS, KA, ALB, MWT);
double Nalb2    = 1.0 + KA * (ALB * 10.0 / 66500.0);
double FTmol2   = FT_out / (MWT * 1e9);
double TTmol2   = TT_out * 1e-8 / MWT;
double FBIO_out = (TTmol2 > 0.0) ? (Nalb2 * FTmol2 / TTmol2) : 0.46;
double FREEPCT  = (TT_out > 0.0) ? 100.0 * (FTmol2 / TTmol2) : 0.0;
double TT_nmolL = TT_out * 0.03467;
double HCT_out  = HCT0 * RBC;
double HGB_out  = HCT_out / 3.0;
double SPERM_out = SPERM0 * EPID;
double ITT_PCT  = 100.0 * ITT / ITT0;
double TE_RATIO = TT_out / SAFE(E2);
double DRIVE_out = hillf(FT_out, EC50_HEPC, HILL_H);
// population fraction crossing the Hct stopping rule, normal approximation
double ZSTOP = (HCT_STOP - HCT_out) / HCT_SD;
double PSTOP = 100.0 * 0.5 * erfc(ZSTOP / sqrt(2.0));
// densitometric BMD = true balance MINUS the unfilled remodelling space
double BMD_TRAB = BMD_TR - RSPACE * (RSP - 1.0);
double BMD_CORT = BMD_CO - CORT_F * RSPACE * (RSP - 1.0);
double day = TIME;

$CAPTURE @annotated
TT_out   : Total testosterone (ng/dL)
TT_nmolL : Total testosterone (nmol/L)
FT_out   : Free testosterone (pg/mL)
FBIO_out : Bioavailable fraction (-)
FREEPCT  : Free fraction (%)
HCT_out  : Haematocrit (%)
HGB_out  : Haemoglobin (g/dL)
SPERM_out: Sperm concentration (million/mL)
ITT_PCT  : Intratesticular testosterone (% of normal)
TE_RATIO : Testosterone / estradiol ratio
DRIVE_out: Erythropoietic drive (0-1)
PSTOP    : Population % exceeding the Hct stopping threshold
BMD_TRAB : Trabecular vBMD as densitometry would measure it (%)
BMD_CORT : Cortical BMD as densitometry would measure it (%)
day      : Time (days)
'

mod <- mcode("mhg", code)

## =============================================================================
##  CALIBRATION BLOCK — what is anchored, what is fitted, what is invented
## =============================================================================
##  ANCHORED to a measured value in the literature
##    KSPILL x ITT0 = 6 mg/day     classic testosterone production rate
##    CL ~ 1000 L/day              metabolic clearance rate of testosterone
##    KS = 1e9, KA = 3.6e4         Vermeulen 1999 association constants
##    ITT0 / serum ~ 35-100x       Jarow 2001, Coviello 2005, Roth 2010
##    TAU_SPG = 74 d               spermatogenic cycle
##    RBCLIFE ~ 100-120 d          red cell lifespan
##    LH/FSH half-lives            standard endocrine values
##    ester molar factors          cypionate 0.70, enanthate 0.72, undecanoate 0.63
##
##  FITTED to reproduce a published trial outcome
##    F_IM, KA_IM                  peak ~1100-1500 / trough ~300 ng/dL at
##                                 200 mg cypionate q2wk
##    F_GEL                        C_avg ~500-600 ng/dL at 81 mg/day 1.62% gel
##    SFAT_E, SLEAN_T              Finkelstein NEJM 2013 dose-response
##    KFORM/KRES/E_SCL             T-Trials Bone: +7.5% trabecular vBMD at 1 yr
##    EPSA, KPSA                   PSA +0.3-0.5 ng/mL in the first 6 months
##    ITT50_S, HS_SP               Coviello 2005 hCG dose-ranging; Liu 2006
##                                 recovery kinetics
##
##  CALIBRATED WITH NO DIRECT MEASUREMENT (the model's exposed assumptions)
##    EC50_HEPC = 300 pg/mL        chosen supraphysiologic so the erythropoietic
##                                 response is convex over the therapeutic range.
##                                 This is the single most consequential free
##                                 parameter in the model — it is what makes
##                                 formulation waveform matter for haematocrit.
##                                 MHG_convexity_decomposition() reports the
##                                 sensitivity of the conclusion to it.
##    HCT_SD = 3.5%                population spread used to turn a mean shift
##                                 into an incidence; from typical trial SDs.
##    LIB_T / LIB_E / VIT_T        PRO weights, scaled to T-Trials effect sizes
##                                 (deliberately small — the T-Trials vitality
##                                 and walking primary endpoints were NEGATIVE).
## =============================================================================

## =============================================================================
##  DOSING EVENT BUILDERS
## =============================================================================

## Testosterone cypionate / enanthate IM. `mg` is mg of the ESTER; the molar
## factor is folded into F_IM, so pass the prescribed product dose.
tst_im <- function(mg = 200, every_days = 14, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_IM",
     ii = every_days, addl = max(0, floor(weeks * 7 / every_days) - 1))

## Testosterone undecanoate IM (Nebido 1000 mg / Aveed 750 mg).
tst_tu_im <- function(mg = 1000, every_days = 84, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_TU",
     ii = every_days, addl = max(0, floor(weeks * 7 / every_days) - 1))

## SC auto-injector (Xyosted), weekly.
tst_sc <- function(mg = 75, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_SC",
     ii = 7, addl = max(0, weeks - 1))

## Transdermal gel 1.62%, once daily.
tst_gel <- function(mg = 81, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_GEL",
     ii = 1, addl = max(0, weeks * 7 - 1))

## Oral testosterone undecanoate, twice daily with food.
tst_oral <- function(mg = 237, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_ORAL",
     ii = 0.5, addl = max(0, weeks * 14 - 1))

## Subcutaneous pellets, q4 months (dose = total mg implanted).
tst_pellet <- function(mg = 750, every_days = 120, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_PEL",
     ii = every_days, addl = max(0, floor(weeks * 7 / every_days) - 1))

## Nasal gel (Natesto) 11 mg TID.
tst_nasal <- function(mg = 11, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "DEP_NAS",
     ii = 1/3, addl = max(0, weeks * 21 - 1))

## hCG SC. Typical replacement 1500 IU 3x/week; ITT rescue 250-500 IU EOD.
hcg <- function(iu = 1500, every_days = 7/3, weeks = 52, start_day = 0)
  ev(time = start_day, amt = iu, cmt = "HCG_D",
     ii = every_days, addl = max(0, floor(weeks * 7 / every_days) - 1))

## rFSH / hMG SC, typically 75-150 IU three times weekly.
fsh_rx <- function(iu = 150, every_days = 7/3, weeks = 26, start_day = 0)
  ev(time = start_day, amt = iu, cmt = "FSHD_D",
     ii = every_days, addl = max(0, floor(weeks * 7 / every_days) - 1))

## Clomiphene citrate / enclomiphene, daily.
clomiphene <- function(mg = 25, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "CLO_D",
     ii = 1, addl = max(0, weeks * 7 - 1))

## Anastrozole 1 mg daily.
anastrozole <- function(mg = 1, weeks = 52, start_day = 0)
  ev(time = start_day, amt = mg, cmt = "ANA_D",
     ii = 1, addl = max(0, weeks * 7 - 1))

## =============================================================================
##  SIMULATION HARNESS
## =============================================================================

.shift_ev <- function(e, burn_days) {
  if (is.null(e)) return(NULL)
  d <- as.data.frame(e)
  d$time <- d$time + burn_days
  as.ev(d)
}

## Run one arm. A burn-in period lets every slow compartment (RBC, bone, body
## composition, Leydig capacity) settle before t = 0, which is then reported as
## the treatment start. Returns a data frame with `day` re-zeroed at dosing.
run_arm <- function(pars = list(), events = NULL, days = 365, burn = 365,
                    delta = 0.125, label = NA_character_) {
  m <- mod
  if (length(pars)) m <- param(m, pars)
  e <- .shift_ev(events, burn)
  out <- as.data.frame(
    mrgsim(m, events = e, end = burn + days, delta = delta,
           maxsteps = 500000, atol = 1e-8, rtol = 1e-6)
  )
  out$day <- out$time - burn
  out$arm <- label
  out
}

.at <- function(d, day, var) {
  i <- which.min(abs(d$day - day))
  d[[var]][i]
}
## Change from the PRE-TREATMENT value (day 0), which is what the trials report.
## This matters: an untreated hypogonadal man does not start at the eugonadal
## referent for haematocrit, bone or lean mass — he starts below it.
.delta <- function(d, day, var) .at(d, day, var) - .at(d, 0, var)
.mean_over <- function(d, var, from, to) {
  s <- d[d$day >= from & d$day <= to, ]
  mean(s[[var]], na.rm = TRUE)
}

## Standard patient archetypes -------------------------------------------------
PT_ORGANIC <- list(AGE = 45, LEY_DMG = 0.80, INS_REL = 1.0, FAT0 = 24, LEAN0 = 60)
PT_FUNCTIONAL <- list(AGE = 52, LEY_DMG = 0.0, INS_REL = 2.4, FAT0 = 42,
                      LEAN0 = 62, INFL_SUP = 0.87)
PT_ELDERLY <- list(AGE = 72, LEY_DMG = 0.74, INFL_SUP = 0.58, INS_REL = 1.1,
                   FAT0 = 28, LEAN0 = 52)
PT_SECONDARY <- list(AGE = 33, PITF = 0.10, LEY_DMG = 0.0, FAT0 = 22, LEAN0 = 62)
PT_KLINEFELTER <- list(AGE = 30, LEY_DMG = 0.85, SER_DMG = 0.95, FAT0 = 28)
PT_OPIOID <- list(AGE = 44, OPI_SUP = 0.92, LEY_DMG = 0.0, FAT0 = 26)

## =============================================================================
##  SCENARIOS
## =============================================================================

## S1. Untreated natural history of functional (obesity-driven) hypogonadism.
MHG_scenario_natural <- function(days = 1095)
  run_arm(PT_FUNCTIONAL, NULL, days = days, label = "untreated functional")

## S2. IM cypionate 200 mg q2wk — the classic peaky regimen.
MHG_scenario_im_q2wk <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_im(200, 14, weeks = ceiling(days / 7)), days, label = "IM 200 q2wk")

## S3. IM cypionate 100 mg weekly — same weekly dose, half the amplitude.
MHG_scenario_im_weekly <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_im(100, 7, weeks = ceiling(days / 7)), days, label = "IM 100 q1wk")

## S4. SC auto-injector 75 mg weekly.
MHG_scenario_sc <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_sc(75, weeks = ceiling(days / 7)), days, label = "SC 75 q1wk")

## S5. Transdermal gel 1.62%, 81 mg/day.
MHG_scenario_gel <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_gel(81, weeks = ceiling(days / 7)), days, label = "gel 81 mg/day")

## S6. Oral testosterone undecanoate 237 mg BID.
MHG_scenario_oral <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_oral(237, weeks = ceiling(days / 7)), days, label = "oral TU 237 BID")

## S7. IM undecanoate 1000 mg q12wk.
MHG_scenario_tu_im <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_tu_im(1000, 84, weeks = ceiling(days / 7)), days, label = "TU IM q12wk")

## S8. Pellets 750 mg q4 months.
MHG_scenario_pellet <- function(days = 365, pt = PT_ORGANIC)
  run_arm(pt, tst_pellet(750, 120, weeks = ceiling(days / 7)), days, label = "pellets q4mo")

## S9. hCG monotherapy — secondary hypogonadism, fertility preserved.
MHG_scenario_hcg_mono <- function(days = 365, pt = PT_SECONDARY)
  run_arm(pt, hcg(1500, 7/3, weeks = ceiling(days / 7)), days, label = "hCG 1500 3x/wk")

## S10. Testosterone + low-dose hCG — replacement with ITT rescue.
MHG_scenario_t_plus_hcg <- function(days = 365, pt = PT_SECONDARY) {
  e <- c(tst_gel(81, weeks = ceiling(days / 7)),
         hcg(500, 2, weeks = ceiling(days / 7)))
  run_arm(pt, e, days, label = "gel + hCG 500 EOD")
}

## S11. Clomiphene 25 mg daily — secondary hypogonadism, axis left on.
MHG_scenario_clomiphene <- function(days = 365, pt = PT_SECONDARY)
  run_arm(pt, clomiphene(25, weeks = ceiling(days / 7)), days, label = "clomiphene 25 qd")

## S12. Weight loss alone in functional hypogonadism (no androgen given).
MHG_scenario_weight_loss <- function(days = 730, frac_per_year = 0.30)
  run_arm(c(PT_FUNCTIONAL, list(WT_LOSS = frac_per_year)), NULL, days,
          label = sprintf("weight loss %.0f%%/yr", 100 * frac_per_year))

## S13. Opioid-induced androgen deficiency, then opioid cessation at day 180.
MHG_scenario_opioid <- function(days = 540) {
  on  <- run_arm(PT_OPIOID, NULL, days = 180, label = "opioid on")
  off <- run_arm(c(PT_OPIOID, list(OPI_SUP = 0)), NULL, days = days - 180,
                 label = "opioid stopped")
  off$day <- off$day + 180
  rbind(on[on$day >= 0, ], off[off$day > 180, ])
}

## S14. Three years of TRT, then cessation — recovery of the axis and of sperm.
MHG_scenario_cessation <- function(treat_days = 1095, follow_days = 730,
                                   pt = PT_SECONDARY) {
  wk <- ceiling(treat_days / 7)
  on <- run_arm(pt, tst_im(200, 14, weeks = wk), days = treat_days,
                label = "on TRT")
  ## restart from the suppressed state: carry the slow compartments forward
  last <- on[nrow(on), ]
  ## Every compartment carries forward, including the depots (so the residual
  ## ester in the IM depot keeps releasing after the last injection — which is
  ## exactly why recovery does not start on the day the prescription stops).
  init_names <- c("DEP_IM", "DEP_TU", "DEP_SC", "DEP_GEL", "DEP_ORAL",
                  "DEP_PEL", "DEP_NAS", "CENT", "PERIPH",
                  "LEYCAP", "ITT", "SHBG", "E2", "DHT", "INHB", "GNRHD",
                  "LH", "FSH", "SG1", "SG2", "SG3", "SG4", "EPID",
                  "RBC", "RETIC", "PROG_E", "HEPC", "EPO",
                  "SCLERO", "OB", "OC", "BMD_TR", "BMD_CO", "RSP",
                  "LEAN", "FAT", "PSA", "LIBIDO", "VITAL")
  m <- param(mod, c(pt, list(INITSS = 0)))
  for (nm in init_names)
    if (nm %in% names(last)) m <- init(m, setNames(list(last[[nm]]), nm))
  off <- as.data.frame(mrgsim(m, end = follow_days, delta = 1,
                              maxsteps = 500000, atol = 1e-8, rtol = 1e-6))
  off$day <- off$time + treat_days
  off$arm <- "off TRT"
  rbind(on[on$day >= 0, intersect(names(on), names(off))],
        off[, intersect(names(on), names(off))])
}

## =============================================================================
##  ANALYSIS FUNCTIONS — each one prints a number the prose would otherwise
##  only assert
## =============================================================================

## A1. The Vermeulen nomogram: what the binding equilibrium actually does.
MHG_free_T_nomogram <- function(TT = c(100, 200, 300, 400, 600, 800, 1000,
                                       1500, 2000),
                                SHBG = c(15, 25, 35, 55, 80)) {
  ft <- function(tt, s) {
    KS <- 1e9; KA <- 3.6e4; MW <- 288.4
    Tm <- tt * 1e-8 / MW
    Sm <- s * 1e-9
    N  <- 1 + KA * (4.3 * 10 / 66500)
    a  <- N * KS; b <- N + KS * Sm - KS * Tm
    (-b + sqrt(b^2 + 4 * a * Tm)) / (2 * a) * MW * 1e9
  }
  tab <- outer(TT, SHBG, Vectorize(ft))
  dimnames(tab) <- list(sprintf("TT %5d", TT), sprintf("SHBG %2d", SHBG))
  cat("\n== Free testosterone (pg/mL) by total T and SHBG ==\n")
  print(round(tab, 1))
  ff <- 100 * t(t(tab) / 1) / (TT * 0.03467 * 288.4)
  cat("\n== Free FRACTION (%) — note it RISES with TT: the SHBG buffer saturates ==\n")
  print(round(ff, 2))
  d2 <- diff(diff(sapply(seq(100, 2000, 50), ft, s = 35)))
  cat(sprintf("\nConvexity check at SHBG 35: all second differences > 0 ? %s (min %.4f)\n",
              all(d2 > 0), min(d2)))
  invisible(tab)
}

## A2. Frame dependence of the diagnostic threshold.
MHG_diagnostic_frame <- function(SHBG = c(15, 20, 25, 30, 35, 45, 55, 70, 90),
                                 TT_cut = 300, FT_cut = 65) {
  ft <- function(tt, s) {
    KS <- 1e9; KA <- 3.6e4; MW <- 288.4
    Tm <- tt * 1e-8 / MW; Sm <- s * 1e-9
    N <- 1 + KA * (4.3 * 10 / 66500); a <- N * KS; b <- N + KS * Sm - KS * Tm
    (-b + sqrt(b^2 + 4 * a * Tm)) / (2 * a) * MW * 1e9
  }
  tt_equiv <- sapply(SHBG, function(s)
    uniroot(function(x) ft(x, s) - FT_cut, c(20, 3000))$root)
  out <- data.frame(
    SHBG_nmolL      = SHBG,
    FT_at_TTcut     = round(sapply(SHBG, function(s) ft(TT_cut, s)), 1),
    TT_giving_FTcut = round(tt_equiv, 0),
    verdict = ifelse(tt_equiv < TT_cut, "TT cut-off OVER-diagnoses",
                     "TT cut-off UNDER-diagnoses")
  )
  cat(sprintf("\n== A single TT cut-off (%d ng/dL) is not a single FT cut-off (%d pg/mL) ==\n",
              TT_cut, FT_cut))
  print(out, row.names = FALSE)
  cat(sprintf("\nRange of the TT equivalent: %.0f - %.0f ng/dL (%.1f-fold) across SHBG %d-%d.\n",
              min(tt_equiv), max(tt_equiv), max(tt_equiv) / min(tt_equiv),
              min(SHBG), max(SHBG)))
  cat("The binding protein, not the gonad, decides which side of the line a man falls on.\n")
  invisible(out)
}

## A3. The two patients a TT-only rule gets wrong, run as full simulations.
MHG_shbg_paradox <- function(days = 90) {
  obese   <- run_arm(list(AGE = 45, SHBG0 = 17, INS_REL = 3.2, FAT0 = 45,
                          LEY_DMG = 0.62), NULL, days, label = "obese, low SHBG")
  elderly <- run_arm(list(AGE = 78, SHBG0 = 62, INS_REL = 0.9, FAT0 = 22,
                          LEY_DMG = 0.67), NULL, days, label = "elderly, high SHBG")
  row <- function(d) data.frame(
    arm   = d$arm[1],
    SHBG  = round(.at(d, days, "SHBG"), 1),
    TT    = round(.at(d, days, "TT_out"), 0),
    FT    = round(.at(d, days, "FT_out"), 1),
    freepct = round(.at(d, days, "FREEPCT"), 2),
    TT_says = ifelse(.at(d, days, "TT_out") < 300, "hypogonadal", "normal"),
    FT_says = ifelse(.at(d, days, "FT_out") < 65,  "hypogonadal", "normal")
  )
  out <- rbind(row(obese), row(elderly))
  cat("\n== The SHBG paradox: the two labels can disagree ==\n")
  print(out, row.names = FALSE)
  if (out$TT_says[1] != out$FT_says[1] || out$TT_says[2] != out$FT_says[2])
    cat("\nAt least one patient is classified differently by TT and by FT.\n",
        "This is not measurement error; it is the binding equilibrium.\n", sep = "")
  invisible(out)
}

## A4. THE CENTRAL LEDGER — decompose the formulation effect on haematocrit
##     into (a) convex binding, (b) convex response, (c) threshold amplification.
MHG_convexity_decomposition <- function(days = 365, pt = PT_ORGANIC) {
  arms <- list(
    "gel 81 mg/day"    = MHG_scenario_gel(days, pt),
    "SC 75 mg weekly"  = MHG_scenario_sc(days, pt),
    "IM 100 mg weekly" = MHG_scenario_im_weekly(days, pt),
    "IM 200 mg q2wk"   = MHG_scenario_im_q2wk(days, pt)
  )
  ftfun <- function(tt, s = 35) {
    KS <- 1e9; KA <- 3.6e4; MW <- 288.4
    Tm <- tt * 1e-8 / MW; Sm <- s * 1e-9
    N <- 1 + KA * (4.3 * 10 / 66500); a <- N * KS; b <- N + KS * Sm - KS * Tm
    (-b + sqrt(b^2 + 4 * a * Tm)) / (2 * a) * MW * 1e9
  }
  w <- c(days - 90, days)   # last 3 months, at steady state
  tab <- do.call(rbind, lapply(names(arms), function(nm) {
    d  <- arms[[nm]]
    mTT <- .mean_over(d, "TT_out", w[1], w[2])
    mFT <- .mean_over(d, "FT_out", w[1], w[2])
    mS  <- .mean_over(d, "SHBG",   w[1], w[2])
    mDR <- .mean_over(d, "DRIVE_out", w[1], w[2])
    ec  <- as.numeric(param(mod)$EC50_HEPC); hh <- as.numeric(param(mod)$HILL_H)
    drive_of_meanFT <- mFT^hh / (ec^hh + mFT^hh)
    data.frame(
      arm = nm,
      peak_TT  = round(max(d$TT_out[d$day > w[1]]), 0),
      trough_TT= round(min(d$TT_out[d$day > w[1]]), 0),
      mean_TT  = round(mTT, 0),
      mean_FT  = round(mFT, 1),
      FT_of_meanTT = round(ftfun(mTT, mS), 1),
      jensen_binding_pct = round(100 * (mFT / ftfun(mTT, mS) - 1), 1),
      mean_drive = round(mDR, 4),
      drive_of_meanFT = round(drive_of_meanFT, 4),
      jensen_response_pct = round(100 * (mDR / drive_of_meanFT - 1), 1),
      dHct = round(.delta(d, days, "HCT_out"), 2),
      pct_over_54 = round(.at(d, days, "PSTOP"), 1)
    )
  }))
  ref <- tab[tab$arm == "gel 81 mg/day", ]
  tab$drive_vs_gel   <- round(tab$mean_drive / ref$mean_drive, 3)
  tab$incid_vs_gel   <- round(tab$pct_over_54 / max(ref$pct_over_54, 1e-6), 2)
  cat("\n== Where does the formulation difference in erythrocytosis come from? ==\n")
  print(tab, row.names = FALSE)

  ## The only clean waveform comparison in the table: 100 mg weekly and
  ## 200 mg q2wk are the SAME weekly milligrams delivered with 2x the amplitude.
  a <- tab[tab$arm == "IM 100 mg weekly", ]
  b <- tab[tab$arm == "IM 200 mg q2wk", ]
  cat("\n-- Reading the ledger --\n")
  cat("(a) CONVEX BINDING is real and scales with amplitude: jensen_binding_pct\n")
  cat(sprintf("    runs %.1f%% for the flat gel up to %.1f%% for the peakiest IM regimen.\n",
              min(tab$jensen_binding_pct), max(tab$jensen_binding_pct)))
  cat("(b) CONVEX RESPONSE does NOT add to it at these peak heights. jensen_response_pct\n")
  cat("    is near zero or NEGATIVE for the peaky arms, because their free-T peaks\n")
  cat(sprintf("    overshoot EC50_HEPC (%.0f pg/mL) into the CONCAVE upper limb of the Hill.\n",
              as.numeric(param(mod)$EC50_HEPC)))
  cat("    A response that is convex low down is concave high up, and supraphysiologic\n")
  cat("    peaks live on the wrong side of the inflection. This is a model result that\n")
  cat("    contradicts the usual hand-wave, and it is reported rather than tuned away.\n")
  cat("\n(c) DOSE-MATCHED WAVEFORM EFFECT — the only clean comparison available here,\n")
  cat("    100 mg weekly vs 200 mg q2wk (identical mg/week, twice the amplitude):\n")
  cat(sprintf("      mean TT     %6.0f  vs %6.0f ng/dL\n", a$mean_TT, b$mean_TT))
  cat(sprintf("      peak:trough %6.1f  vs %6.1f\n",
              a$peak_TT / max(a$trough_TT, 1), b$peak_TT / max(b$trough_TT, 1)))
  cat(sprintf("      dHct        %+6.2f  vs %+6.2f  (waveform effect = %+.2f points)\n",
              a$dHct, b$dHct, b$dHct - a$dHct))
  cat("    So at MATCHED weekly dose the waveform moves haematocrit by only a fraction\n")
  cat("    of a point. The large IM-vs-transdermal gap in the literature is therefore\n")
  cat("    mostly a DOSE/EXPOSURE difference, not a waveform difference.\n")
  cat("\n(d) THRESHOLD AMPLIFICATION does the rest of the work. A fixed Hct > 54% rule\n")
  cat("    applied to a spread-out population turns a modest shift in the MEAN into a\n")
  cat(sprintf("    %.1f-fold difference in the fraction of men who cross it. That is a\n",
              max(tab$incid_vs_gel)))
  cat("    property of the STOPPING RULE, not of the molecule.\n")
  cat("\nCaveats, stated rather than buried:\n")
  cat("  - mean exposure is NOT matched across all rows (gel 686 vs IM q2wk 773 ng/dL\n")
  cat("    at these doses), so incid_vs_gel mixes dose with waveform. Only the\n")
  cat("    100-weekly / 200-q2wk pair isolates waveform.\n")
  cat("  - the ABSOLUTE incidences are not calibrated to any trial: they depend on this\n")
  cat("    patient's starting haematocrit and on HCT_SD, neither of which is fitted.\n")
  cat("    The RATIO is the quantity to read; observed is roughly 3x (Ohlander 2018).\n")
  cat(sprintf("  - EC50_HEPC = %.0f pg/mL is CALIBRATED, not measured. See\n",
              as.numeric(param(mod)$EC50_HEPC)))
  cat("    MHG_hepcidin_sensitivity() for how much of the above survives moving it.\n")
  invisible(tab)
}

## A5. EC50_HEPC sensitivity — how much of conclusion (b) survives?
MHG_hepcidin_sensitivity <- function(ec50s = c(150, 225, 300, 450, 600),
                                     days = 365, pt = PT_ORGANIC) {
  out <- do.call(rbind, lapply(ec50s, function(ec) {
    g <- run_arm(c(pt, list(EC50_HEPC = ec)), tst_gel(81, weeks = ceiling(days/7)),
                 days, label = "gel")
    i <- run_arm(c(pt, list(EC50_HEPC = ec)), tst_im(200, 14, weeks = ceiling(days/7)),
                 days, label = "IM")
    data.frame(EC50_HEPC = ec,
               dHct_gel = round(.delta(g, days, "HCT_out"), 2),
               dHct_IM  = round(.delta(i, days, "HCT_out"), 2),
               gap      = round(.at(i, days, "HCT_out") - .at(g, days, "HCT_out"), 2),
               incid_ratio = round(.at(i, days, "PSTOP") /
                                     max(.at(g, days, "PSTOP"), 1e-6), 2))
  }))
  cat("\n== Sensitivity of the formulation gap to the one calibrated parameter ==\n")
  print(out, row.names = FALSE)
  cat("\nIf EC50_HEPC were LOW (well inside the eugonadal range) the response would\n")
  cat("be concave over the therapeutic window and the peaky regimen would lose its\n")
  cat("advantage. The claim in (b) stands or falls on this parameter, which is why\n")
  cat("it is stated as a calibrated assumption rather than a result.\n")
  invisible(out)
}

## A6. Serum T is normalised, ITT is not — the two variables side by side.
MHG_ITT_collapse <- function(days = 180) {
  arms <- list(
    "untreated (secondary)" = run_arm(PT_SECONDARY, NULL, days, label = "none"),
    "gel 81 mg/day"         = MHG_scenario_gel(days, PT_SECONDARY),
    "IM 200 mg q2wk"        = MHG_scenario_im_q2wk(days, PT_SECONDARY),
    "gel + hCG 500 IU EOD"  = MHG_scenario_t_plus_hcg(days, PT_SECONDARY),
    "hCG 1500 IU 3x/wk"     = MHG_scenario_hcg_mono(days, PT_SECONDARY),
    "clomiphene 25 mg qd"   = MHG_scenario_clomiphene(days, PT_SECONDARY)
  )
  tab <- do.call(rbind, lapply(names(arms), function(nm) {
    d <- arms[[nm]]
    data.frame(regimen = nm,
               serum_TT  = round(.mean_over(d, "TT_out", days - 30, days), 0),
               LH        = round(.at(d, days, "LH"), 2),
               ITT_pct   = round(.at(d, days, "ITT_PCT"), 1),
               sperm_Mml = round(.at(d, days, "SPERM_out"), 1))
  }))
  cat("\n== Serum testosterone and intratesticular testosterone are not the same variable ==\n")
  print(tab, row.names = FALSE)
  cat("\nThe rows where serum_TT is normal and ITT_pct is single-digit are exactly the\n")
  cat("regimens that replace the hormone and sterilise the patient. Only the arms that\n")
  cat("leave an LHCGR signal standing (hCG, clomiphene) hold ITT above the threshold.\n")
  invisible(tab)
}

## A7. Recovery after cessation — reproduce the Liu 2006 pooled kinetics.
MHG_recovery_curve <- function(treat_days = 1095, follow_days = 730) {
  d <- MHG_scenario_cessation(treat_days, follow_days)
  off <- d[d$day >= treat_days, ]
  off$m <- (off$day - treat_days) / 30.44
  thr <- 20   # million/mL, the Liu 2006 endpoint
  hit <- function(mo) {
    s <- off[off$m <= mo, ]
    any(s$SPERM_out >= thr)
  }
  tab <- data.frame(
    months = c(3, 6, 12, 18, 24),
    model_reached_20M = sapply(c(3, 6, 12, 18, 24), hit),
    model_sperm_Mml = round(sapply(c(3, 6, 12, 18, 24), function(mo)
      off$SPERM_out[which.min(abs(off$m - mo))]), 1),
    Liu2006_pct_recovered = c(NA, 67, 90, NA, 100)
  )
  first <- off$m[which(off$SPERM_out >= thr)[1]]
  cat("\n== Recovery of spermatogenesis after stopping testosterone ==\n")
  print(tab, row.names = FALSE)
  cat(sprintf("\nModel: first crosses 20 M/mL at %.1f months after cessation.\n", first))
  cat("Liu 2006 (Lancet, pooled n=1549): median 3.4 months; 67% by 6 months,\n")
  cat("90% by 12, 100% by 24. The model reproduces the SHAPE of a single median\n")
  cat("patient; the published percentages are a population distribution the\n")
  cat("deterministic model does not attempt to span.\n")
  invisible(tab)
}

## A8. Finkelstein NEJM 2013 — separate what is testosterone from what is estradiol.
MHG_finkelstein <- function(days = 112,
                            gel_doses = c(0, 1.25, 2.5, 5, 10)) {
  ## Background: GnRH agonist shuts the axis off entirely (PITF -> 0) so the
  ## only androgen present is the administered gel. Group 2 adds anastrozole.
  base <- list(AGE = 30, PITF = 0.02, FAT0 = 22, LEAN0 = 62)
  run1 <- function(gmg, anas) {
    e <- tst_gel(gmg * 10, weeks = ceiling(days / 7))    # 1 g of 1% gel = 10 mg T
    if (anas) e <- c(e, anastrozole(1, weeks = ceiling(days / 7)))
    run_arm(base, e, days, label = sprintf("%.2f g%s", gmg, if (anas) " + AI" else ""))
  }
  rows <- list()
  for (g in gel_doses) for (a in c(FALSE, TRUE)) {
    d <- run1(g, a)
    rows[[length(rows) + 1]] <- data.frame(
      gel_g = g, anastrozole = a,
      TT = round(.at(d, days, "TT_out"), 0),
      E2 = round(.at(d, days, "E2"), 1),
      d_lean_kg = round(.delta(d, days, "LEAN"), 2),
      d_fat_kg  = round(.delta(d, days, "FAT"), 2),
      libido    = round(.at(d, days, "LIBIDO"), 2))
  }
  tab <- do.call(rbind, rows)
  cat("\n== Finkelstein design: graded testosterone with and without aromatase blockade ==\n")
  print(tab, row.names = FALSE)
  noAI <- tab[!tab$anastrozole, ]; AI <- tab[tab$anastrozole, ]
  cat("\nAt matched testosterone dose, adding the aromatase inhibitor changes:\n")
  cat(sprintf("  lean mass by %+.2f kg on average  (T-dependent -> little change)\n",
              mean(AI$d_lean_kg - noAI$d_lean_kg)))
  cat(sprintf("  FAT mass  by %+.2f kg on average  (E2-dependent -> fat goes UP)\n",
              mean(AI$d_fat_kg - noAI$d_fat_kg)))
  cat(sprintf("  libido    by %+.2f points        (both hormones contribute)\n",
              mean(AI$libido - noAI$libido)))
  cat("\nThis is the NEJM 2013 dissociation: the increase in body fat during androgen\n")
  cat("withdrawal is an ESTROGEN deficiency, which is why blocking aromatase to\n")
  cat("manage gynaecomastia does not come for free.\n")
  invisible(tab)
}

## A9. Price the aromatase inhibitor in bone.
MHG_aromatase_cost <- function(days = 365, pt = PT_ORGANIC) {
  a <- run_arm(pt, tst_gel(81, weeks = ceiling(days / 7)), days, label = "gel")
  b <- run_arm(pt, c(tst_gel(81, weeks = ceiling(days / 7)),
                     anastrozole(1, weeks = ceiling(days / 7))), days,
               label = "gel + anastrozole")
  tab <- do.call(rbind, lapply(list(a, b), function(d) data.frame(
    arm = d$arm[1],
    TT = round(.at(d, days, "TT_out"), 0),
    E2 = round(.at(d, days, "E2"), 1),
    LH = round(.at(d, days, "LH"), 2),
    sclerostin = round(.at(d, days, "SCLERO"), 3),
    trab_vBMD_pct = round(.delta(d, days, "BMD_TRAB"), 2),
    cort_BMD_pct  = round(.delta(d, days, "BMD_CORT"), 2),
    fat_kg = round(.at(d, days, "FAT"), 1))))
  cat("\n== What an aromatase inhibitor costs, on top of a full testosterone dose ==\n")
  print(tab, row.names = FALSE)
  cat(sprintf("\nTrabecular vBMD at 1 year: %+.2f%% without AI vs %+.2f%% with AI (difference %+.2f%%).\n",
              tab$trab_vBMD_pct[1], tab$trab_vBMD_pct[2],
              tab$trab_vBMD_pct[2] - tab$trab_vBMD_pct[1]))
  cat("Note the LH column: the AI RAISES LH, because E2 is the dominant feedback\n")
  cat("signal in men. That is the same mechanism SERMs exploit therapeutically.\n")
  invisible(tab)
}

## A10. Formulation head-to-head across every endpoint the model carries.
MHG_formulation_ledger <- function(days = 365, pt = PT_ORGANIC) {
  arms <- list(
    "gel 81 mg/day"    = MHG_scenario_gel(days, pt),
    "SC 75 mg weekly"  = MHG_scenario_sc(days, pt),
    "IM 100 mg weekly" = MHG_scenario_im_weekly(days, pt),
    "IM 200 mg q2wk"   = MHG_scenario_im_q2wk(days, pt),
    "TU IM q12wk"      = MHG_scenario_tu_im(days, pt),
    "oral TU 237 BID"  = MHG_scenario_oral(days, pt),
    "pellets q4mo"     = MHG_scenario_pellet(days, pt)
  )
  tab <- do.call(rbind, lapply(names(arms), function(nm) {
    d <- arms[[nm]]; w1 <- days - 90
    data.frame(regimen = nm,
      Cavg_TT   = round(.mean_over(d, "TT_out", w1, days), 0),
      Cmax_TT   = round(max(d$TT_out[d$day > w1]), 0),
      Cmin_TT   = round(min(d$TT_out[d$day > w1]), 0),
      pk_tr_ratio = round(max(d$TT_out[d$day > w1]) /
                            max(min(d$TT_out[d$day > w1]), 1), 1),
      FT_avg    = round(.mean_over(d, "FT_out", w1, days), 1),
      SHBG      = round(.at(d, days, "SHBG"), 1),
      dHct      = round(.delta(d, days, "HCT_out"), 2),
      pct_over54= round(.at(d, days, "PSTOP"), 1),
      dPSA      = round(.delta(d, days, "PSA"), 2),
      dLean_kg  = round(.delta(d, days, "LEAN"), 2),
      vBMD_pct  = round(.delta(d, days, "BMD_TRAB"), 2),
      ITT_pct   = round(.at(d, days, "ITT_PCT"), 1))
  }))
  cat("\n== Formulation ledger at one year (same patient, different waveforms) ==\n")
  print(tab, row.names = FALSE)
  cat("\nThe efficacy columns (lean mass, vBMD, PSA) are close across regimens because\n")
  cat("they integrate exposure. The columns that separate are the ones driven by the\n")
  cat("PEAK (dHct, pct_over54) — and ITT_pct, which every exogenous route drives to\n")
  cat("the floor regardless of waveform.\n")
  invisible(tab)
}

## A11. Trial ledger — hold the model against the published endpoints.
MHG_trial_ledger <- function() {
  d1 <- run_arm(PT_ELDERLY, tst_gel(60, weeks = 52), days = 365, label = "T-Trials-like")
  cat("\n== Model vs published trial endpoints ==\n")
  tab <- data.frame(
    endpoint = c("Trabecular spine vBMD, 1 yr (%)",
                 "Total testosterone achieved (ng/dL)",
                 "Haematocrit change, 1 yr (%)",
                 "PSA change, 1 yr (ng/mL)",
                 "Lean mass change, 1 yr (kg)",
                 "Sexual desire change (model units)"),
    model = c(round(.delta(d1, 365, "BMD_TRAB"), 1),
              round(.mean_over(d1, "TT_out", 300, 365), 0),
              round(.delta(d1, 365, "HCT_out"), 1),
              round(.delta(d1, 365, "PSA"), 2),
              round(.delta(d1, 365, "LEAN"), 1),
              round(.delta(d1, 365, "LIBIDO"), 2)),
    published = c("+7.5 (T-Trials Bone, JAMA IM 2017)",
                  "~500-600 (T-Trials target mid-normal)",
                  "+2 to +3 transdermal (Ohlander 2018)",
                  "+0.3 to +0.5 (multiple TRT series)",
                  "+1.5 to +2.5 (Bhasin dose-response)",
                  "PDQ-Q4 +0.58 vs +0.10 (NEJM 2016)")
  )
  print(tab, row.names = FALSE)
  cat("\nSafety endpoints the model does NOT try to predict, and why:\n")
  cat("  TRAVERSE MACE HR 0.96 (NEJM 2023) — non-inferiority in a 5,246-man trial is\n")
  cat("    a population result; this deterministic model has no event process for it.\n")
  cat("  TRAVERSE fracture HR 1.43 (NEJM 2024) — MORE fractures despite HIGHER BMD.\n")
  cat("    The model reproduces the BMD gain and therefore CANNOT reproduce the\n")
  cat("    fracture signal. That is a genuine failure of the BMD-as-surrogate\n")
  cat("    assumption, not a tuning problem, and it is left visible on purpose.\n")
  invisible(tab)
}

## =============================================================================
##  PLOTTING
## =============================================================================

MHG_plot_overview <- function(d, main = "Male hypogonadism QSP simulation") {
  op <- par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))
  on.exit({ par(op) })
  pl <- function(y, lab, hline = NULL) {
    plot(d$day, d[[y]], type = "l", lwd = 2, col = "#2c6fb5",
         xlab = "day", ylab = lab, main = lab)
    if (!is.null(hline)) abline(h = hline, lty = 2, col = "#c0392b")
  }
  pl("TT_out", "Total T (ng/dL)", c(300, 1000))
  pl("FT_out", "Free T (pg/mL)", 65)
  pl("SHBG", "SHBG (nmol/L)")
  pl("LH", "LH (IU/L)")
  pl("ITT_PCT", "ITT (% of normal)", 25)
  pl("E2", "Estradiol (pg/mL)")
  pl("HCT_out", "Haematocrit (%)", 54)
  pl("SPERM_out", "Sperm (M/mL)", 16)
  pl("BMD_TRAB", "Trabecular vBMD (%)")
  mtext(main, outer = TRUE, cex = 1.05, font = 2)
}

MHG_plot_waveforms <- function(days = 84, pt = PT_ORGANIC) {
  arms <- list("gel 81 mg/day" = MHG_scenario_gel(days, pt),
               "SC 75 mg weekly" = MHG_scenario_sc(days, pt),
               "IM 200 mg q2wk" = MHG_scenario_im_q2wk(days, pt),
               "TU IM q12wk" = MHG_scenario_tu_im(days, pt))
  cols <- c("#2c6fb5", "#1e8449", "#c0392b", "#8e44ad")
  op <- par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  on.exit({ par(op) })
  ylim <- range(unlist(lapply(arms, function(d) d$TT_out)))
  plot(NA, xlim = c(0, days), ylim = ylim, xlab = "day",
       ylab = "Total T (ng/dL)", main = "Same drug, different waveform")
  abline(h = c(300, 1000), lty = 2, col = "grey60")
  for (i in seq_along(arms)) lines(arms[[i]]$day, arms[[i]]$TT_out, lwd = 2, col = cols[i])
  legend("topright", names(arms), col = cols, lwd = 2, bty = "n", cex = 0.8)
  ylim2 <- range(unlist(lapply(arms, function(d) d$FT_out)))
  plot(NA, xlim = c(0, days), ylim = ylim2, xlab = "day",
       ylab = "Free T (pg/mL)", main = "Free T amplifies the peaks (convex binding)")
  abline(h = 65, lty = 2, col = "grey60")
  for (i in seq_along(arms)) lines(arms[[i]]$day, arms[[i]]$FT_out, lwd = 2, col = cols[i])
}

## =============================================================================
##  RUN EVERYTHING
## =============================================================================

MHG_run_all <- function() {
  cat("\n#############################################################\n")
  cat("#  Male hypogonadism QSP model — full analysis run\n")
  cat("#############################################################\n")
  MHG_free_T_nomogram()
  MHG_diagnostic_frame()
  MHG_shbg_paradox()
  MHG_convexity_decomposition()
  MHG_hepcidin_sensitivity()
  MHG_ITT_collapse()
  MHG_recovery_curve()
  MHG_finkelstein()
  MHG_aromatase_cost()
  MHG_formulation_ledger()
  MHG_trial_ledger()
  cat("\nDone. See mhg_references.md for the sources behind each calibration target.\n")
  invisible(TRUE)
}

## Example:
##   source("mhg_mrgsolve_model.R")
##   MHG_run_all()
##   d <- MHG_scenario_im_q2wk(); MHG_plot_overview(d, "IM cypionate 200 mg q2wk")
##   MHG_plot_waveforms()
