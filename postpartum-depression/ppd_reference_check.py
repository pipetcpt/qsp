#!/usr/bin/env python3
# =============================================================================
#  ppd_reference_check.py
#  Postpartum depression (PPD) QSP model — INDEPENDENT REFERENCE IMPLEMENTATION
#  of the ODE system coded in ppd_mrgsolve_model.R
#  ---------------------------------------------------------------------------
#  WHY THIS FILE EXISTS
#  --------------------
#  The mrgsolve model cannot be executed in this environment (no R runtime).
#  Every quantitative claim made in README.md and in the header of
#  ppd_mrgsolve_model.R is therefore produced HERE, by a second, dependency-free
#  (pure standard library, fixed-step RK4) implementation of the SAME equations
#  with the SAME parameter values.  If the two files ever disagree, THIS is the
#  one that was actually run, and ppd_reference_output.txt is its verbatim
#  output.
#
#  Run:  python3 ppd_reference_check.py > ppd_reference_output.txt
#
#  ---------------------------------------------------------------------------
#  THE ONE IDEA
#  ---------------------------------------------------------------------------
#  Tonic inhibition is a PRODUCT of two factors with very different time
#  constants:
#
#      G_tonic(t) = R_delta(t) * [ 1 + Emax * PAM(t)^h / (EC50^h + PAM(t)^h) ]
#                   \__________/   \_______________________________________/
#                    receptor arm              ligand arm
#                    tau ~ 7 DAYS              tau ~ 12 HOURS
#
#  During pregnancy the receptor arm is homeostatically down-regulated so that
#  the product stays at its non-pregnant set-point (G_target = 1):
#
#      R_delta,ss(PAM) = G_target / (1 + potentiation(PAM))
#
#  At delivery the ligand arm collapses in hours; the receptor arm needs weeks.
#  A product of a fast fall and a slow rise TROUGHS.  Nothing in this file
#  asserts a trough: it falls out of the two time constants.
#
#  Everything else in the model decides (a) how deep the excitatory load pushes
#  the circuit while the trough is open, and (b) whether the mother's
#  symptom-sleep feedback loop is strong enough to LATCH into a depressed state
#  that then outlives the trough.  That latch is why a 60-hour infusion can buy
#  30 days of benefit.
#
#  STATE VECTOR (38 ODEs) — identical ordering to the R model
#  ---------------------------------------------------------------------------
#   0 PLAC   placental functional mass                             rel (1 = term)
#   1 P4     plasma progesterone                                   nM
#   2 DHP    plasma 5alpha-dihydroprogesterone                     nM
#   3 ALLOP  plasma allopregnanolone                               nM
#   4 ALLOB  brain allopregnanolone (effector species)             nM
#   5 E2     plasma estradiol                                      nM
#   6 PCRH   placental CRH                                         rel (1 = term)
#   7 HCRH   hypothalamic CRH drive                                rel (1 = normal)
#   8 ACTH   ACTH                                                  rel
#   9 CORT   cortisol                                              ug/dL
#  10 GRFN   glucocorticoid-receptor feedback sensitivity          rel
#  11 RD     extrasynaptic delta-GABA_A receptor pool              rel
#  12 RG     synaptic gamma2-GABA_A receptor pool                  rel
#  13 KCC2   KCC2 surface expression (Cl- extrusion)               rel
#  14 MAOA   MAO-A binding density                                 rel
#  15 FIVEHT serotonergic tone                                     rel
#  16 AUTO   5-HT1A autoreceptor sensitivity (SSRI latency gate)   rel
#  17 BDNF   BDNF signalling                                       rel
#  18 SYN    structural/synaptic plasticity index                  rel
#  19 INFL   inflammation (IL-6-equivalent)                        pg/mL
#  20 KYNR   kynurenine/tryptophan ratio                           rel
#  21 SLP    sleep obtained                                        h/night
#  22 SDEBT  cumulative sleep debt                                 h
#  23 EXC    excitation/inhibition imbalance index                 rel
#  24 SYMP   latent depressive symptom load                        0-1
#  25 BONDI  mother-infant bonding impairment                      0-1
#  26 BRX1   brexanolone, central                                  ug
#  27 BRX2   brexanolone, peripheral                               ug
#  28 ZURA   zuranolone, absorption depot                          mg
#  29 ZUR1   zuranolone, central                                   mg
#  30 ZUR2   zuranolone, peripheral                                mg
#  31 SERA   sertraline, absorption depot                          mg
#  32 SERC   sertraline, central                                   mg
#  33 ESKC   esketamine, central                                   mg
#  34 AMPAS  AMPA-mediated plasticity surge (ketamine arm)         rel
#  35 CBTP   psychotherapy progress                                0-1
#  36 MILKD  cumulative drug delivered to infant via milk          ug/kg
#  37 INFP   infant plasma concentration (ALLO-equivalents)        nM
#
#  CALIBRATION TARGETS (all fetched from PubMed, see ppd_references.md)
#  ---------------------------------------------------------------------------
#   Kanes 2017 Lancet 390:480 (PMID 28619476), phase 2, severe (HAM-D >=26):
#       60 h HAM-D change  brexanolone -21.0  vs placebo  -8.8
#   Meltzer-Brody 2018 Lancet 392:1058 (PMID 30177236), phase 3 study 1 (>=26):
#       60 h HAM-D change  BRX60 -19.5 · BRX90 -17.7 · placebo -14.0
#   ditto study 2 (HAM-D 20-25):
#       60 h HAM-D change  BRX90 -14.6 · placebo -12.1
#   Deligiannidis 2023 Am J Psychiatry 180:668 (PMID 37491938), SKYLARK:
#       day 15 HAM-D change zuranolone 50 mg -15.6 vs placebo -11.6
#       (significant also at days 3, 28 and 45)
#   Deligiannidis 2021 JAMA Psychiatry 78:951 (PMID 34190962), ROBIN:
#       day 15 HAM-D change zuranolone 30 mg -17.8 vs placebo -13.6
#
#  ASSUMPTIONS THAT MUST BE READ BEFORE ANY NUMBER BELOW IS QUOTED
#  ---------------------------------------------------------------------------
#  A1  Neurosteroid concentrations are handled as TOTAL plasma (and total
#      brain) concentrations, and the potentiation EC50 is expressed in the
#      same currency.  Protein binding (>99 % for both allopregnanolone and
#      zuranolone) and brain partitioning are therefore LUMPED INTO EC50.
#      Absolute EC50 here is not comparable to a patch-clamp EC50.
#  A2  Brexanolone IS allopregnanolone.  Its brain:plasma partition
#      coefficient and its intrinsic potency are therefore set EQUAL to the
#      endogenous values by identity, not fitted.  This is the one place the
#      model gets a free, unfalsifiable-looking parameter for free.
#  A3  Zuranolone is a synthetic analogue; a single lumped factor ZUR_EQ
#      converts its plasma concentration into brain allopregnanolone
#      equivalents.  ZUR_EQ is the ONLY parameter calibrated against a
#      zuranolone endpoint (SKYLARK day 15).  Everything else about the
#      zuranolone arms is then a prediction.
#  A4  Third-trimester plasma allopregnanolone is taken as ~80 nM and the
#      day 4-5 postpartum floor as ~2 nM.  Published values span a wide range
#      because RIA and LC-MS/MS assays disagree by severalfold; the model's
#      behaviour depends on the RATIO (~40x) far more than on the absolute
#      values, which is why the ratio is what was matched.
#  A5  The "nonspecific care" term (k_care, k_nsp) is a CALIBRATED,
#      NON-MECHANISTIC component representing expectancy, structured clinical
#      contact and trial procedures.  It is fitted ONCE to the placebo arms
#      and then held fixed across every active arm, so that the drug effect
#      reported is always the INCREMENT over an equally-cared-for placebo.
#      The model does NOT claim to explain expectancy.
#  A6  Reference enrolment is day 21 postpartum.  Real trials enrolled over a
#      much wider window (up to 6-12 months), and Table 7 shows how strongly
#      the predicted placebo response depends on this choice.  Any comparison
#      of an absolute predicted HAM-D change with a trial number inherits this
#      uncertainty.
# =============================================================================

import math

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
P = dict(
    # ---- body / conversions -----------------------------------------------
    WT=70.0,             # maternal weight, kg
    MW_ALLO=318.5,       # g/mol, allopregnanolone (= brexanolone)
    MW_ZUR=376.5,        # g/mol, zuranolone
    MW_SER=306.2,        # g/mol, sertraline

    # ---- 1. placental / neurosteroid synthesis ---------------------------
    KDEL=4.6,            # 1/h placental involution after expulsion (t1/2 9 min)
    KEL_P4=0.1733,       # 1/h progesterone elimination (t1/2 4 h)
    KSYN_P4=68.98,       # nM/h placental progesterone output at term
    IN_P4=0.3466,        # nM/h non-placental progesterone (residual, ~2 nM)
    K5A=0.30,            # 1/h SRD5A1 flux P4 -> 5a-DHP
    K3A=0.90,            # 1/h AKR1C2 flux 5a-DHP -> ALLO
    KEL_DHP=0.30,        # 1/h other 5a-DHP loss
    F_ALLO=0.0514,       # fraction of 3a-HSD flux appearing as plasma ALLO
    KEL_ALLO=0.0578,     # 1/h plasma ALLO elimination (t1/2 12 h)
    IN_ALLO=0.0867,      # nM/h adrenal/CNS-derived ALLO floor (~1.5 nM)
    KOUT_B=0.35,         # 1/h brain ALLO equilibration (t1/2 2 h)
    KP_ALLO=2.0,         # brain:plasma partition (applies to brexanolone too, A2)
    IN_ALLOB=0.35,       # nM/h local (astrocytic) brain synthesis floor
    KEL_E2=0.1155,       # 1/h estradiol elimination (t1/2 6 h)
    KSYN_E2=6.92,        # nM/h placental estradiol output at term
    IN_E2=0.00924,       # nM/h ovarian floor during lactational amenorrhoea
    KEL_PCRH=0.693,      # 1/h placental CRH clearance (t1/2 1 h)

    # ---- 2. HPA axis ------------------------------------------------------
    KH=0.000990,         # 1/h hypothalamic CRH recovery (t1/2 700 h ~ 29 d)
    PC50=0.55,           # placental-CRH level suppressing hypothalamic CRH 50 %
    KA_ACTH=1.0,         # 1/h ACTH turnover
    W_PCRH=1.60,         # placental CRH contribution to ACTH drive
    CORT_FB=12.0,        # ug/dL cortisol for half-maximal GR feedback
    N_FB=1.5,            # Hill coefficient, cortisol feedback
    KC=6.73,             # (ug/dL)/h cortisol production per unit ACTH
    KEL_CORT=0.462,      # 1/h cortisol elimination (t1/2 1.5 h)
    KG=0.00963,          # 1/h GR-sensitivity adaptation (t1/2 72 h)
    GR50=25.0,           # ug/dL cortisol for half GR down-regulation
    CORT_NP=10.0,        # ug/dL non-pregnant reference cortisol

    # ---- 3. GABA_A receptor plasticity (THE SLOW ARM) --------------------
    G_TARGET=1.0,        # homeostatic tonic-conductance set-point
    EMAX_PAM=2.5,        # maximal fractional potentiation, extrasynaptic
    EC50_PAM=120.0,      # nM brain ALLO-equivalents (see A1)
    H_PAM=1.4,           # Hill coefficient, potentiation
    EMAX_SYN=0.60,       # maximal potentiation, synaptic gamma2 receptors
    EC50_SYN=400.0,      # nM, synaptic receptors are less steroid-sensitive
    KR=0.0057566,         # 1/h delta-pool plasticity (t1/2 120 h = 5 d)
    V_KR=0.2808,           # vulnerability slows receptor recovery
    KR_BOOST=0.0,        # does PAM exposure ACCELERATE receptor plasticity?
                         #  0 = pure-bridge hypothesis (drug only substitutes)
    W_TONIC=0.90,        # weight of tonic inhibition in effective inhibition
    W_PHASIC=0.10,       # weight of phasic inhibition (zuranolone also hits gamma2)
    RG_GAIN=0.70,        # reciprocal up-regulation of synaptic pool
    KK=0.0289,           # 1/h KCC2 turnover (t1/2 24 h)
    KCC2_CORT=0.60,      # cortisol sensitivity of KCC2 loss
    PAM_NP=6.0,          # nM brain ALLO-equivalents, non-pregnant reference

    # ---- 4. monoamine / plasticity ---------------------------------------
    KM_MAOA=0.01444,     # 1/h MAO-A adaptation (t1/2 48 h)
    F_MAOA=0.21,         # maximal MAO-A rise on complete E2 withdrawal
    E250=0.50,           # nM estradiol restraining MAO-A half-maximally
    K5HT=0.0578,         # 1/h serotonergic tone turnover (t1/2 12 h)
    EMAX_SSRI=1.20,      # maximal 5-HT gain from full SERT blockade
    EC50_SER=9.0,        # ng/mL sertraline for 50 % SERT occupancy
    W_KYN_5HT=0.30,      # substrate diversion penalty on 5-HT synthesis
    KAUTO=0.00289,       # 1/h 5-HT1A autoreceptor desensitisation (t1/2 240 h)
    AUTO_BRAKE=0.85,     # fraction of SSRI effect gated by the autoreceptor
    OCC50_AUTO=0.50,     # SERT occupancy desensitising the autoreceptor 50 %
    KB_BDNF=0.00578,     # 1/h BDNF turnover (t1/2 120 h = 5 d)
    KSYN_PL=0.00289,     # 1/h structural plasticity turnover (t1/2 240 h = 10 d)
    W_AMPA=0.80,         # AMPA surge -> BDNF gain (ketamine arm)
    W_PAM_BDNF=0.0,      # does restoring tonic inhibition drive plasticity?
                         #  0 = pure-bridge; >0 = drug has a plasticity arm

    # ---- 5. inflammation / kynurenine ------------------------------------
    KI_INFL=0.0289,      # 1/h IL-6 turnover (t1/2 24 h)
    INFL_BASE=2.0,       # pg/mL baseline IL-6-equivalent
    DEL_INFL=22.0,       # pg/mL delivery inflammatory surge (bolus at t = 0)
    W_SD_INFL=0.50,      # sleep debt -> inflammation
    KK_KYN=0.01444,      # 1/h KYN/TRP turnover (t1/2 48 h)
    W_INFL_KYN=1.00,     # inflammation -> IDO1 -> KYN/TRP
    W_CORT_KYN=0.40,     # cortisol -> TDO -> KYN/TRP

    # ---- 6. sleep ---------------------------------------------------------
    SLP_MAX=8.0,         # h/night achievable sleep
    SLP_NEED=7.5,        # h/night requirement
    A_WAKE=0.18,         # infant-waking penalty on sleep
    A_SYMP=1.482,         # insomnia symptom penalty on sleep (closes the loop)
    KS_SLP=0.10,         # 1/h sleep-state adaptation
    KDEC_SD=0.010,       # 1/h sleep-debt repayment (t1/2 69 h)
    WAKE_LATE_PREG=0.30, # third-trimester nocturnal disturbance
    WAKE_PP0=1.80,       # amplitude of newborn night-waking load
    WAKE_PP_FLOOR=1.00,  # asymptotic night-waking load
    TAU_WAKE=1080.0,     # h, consolidation of infant sleep (t ~ 45 d)
    WAKE_PROTECT=1.00,   # night-waking load under protected-sleep conditions

    # ---- 7. excitation / symptom transfer --------------------------------
    KE_EXC=0.0578,       # 1/h E/I index adaptation (t1/2 12 h)
    W_SD=2.0,           # sleep debt -> excitatory load
    W_INFL=0.50,         # inflammation -> excitatory load
    W_KYN=0.30,          # QUIN/KYNA shift -> excitatory load
    W_5HT=0.60,          # low serotonergic tone -> excitatory load
    W_HPA=0.30,          # blunted-HPA contribution (postpartum hypocortisolism)
    W_SELF=0.22651,         # symptom self-reinforcement (rumination · DMN · glutamate)
    W_SYN=0.80,          # structural plasticity buffers excitatory load
    THR0=2.634,           # E/I threshold for symptom accrual at V = 1
    V_THR=1.00,          # exponent: vulnerability narrows the reserve
    KON=0.0036432,           # 1/h symptom accrual rate
    KOFF=0.0017154,          # 1/h intrinsic symptom resolution, scaled by 1/V
    KFAST=84.83,          # acceleration of resolution below threshold
    SMAX=1.0,            # maximum latent symptom load
    SMIN=0.155,          # residual symptom floor (HAM-D ~7 = remission bound)
    KBOND=0.00963,       # 1/h bonding impairment turnover (t1/2 72 h)
    W_BOND=0.85,         # symptom -> bonding impairment

    # ---- 8. scales --------------------------------------------------------
    HAMD_0=2.0,          # HAM-D17 intercept
    HAMD_SC=32.0,        # HAM-D17 span
    EPDS_0=1.0,          # EPDS intercept
    EPDS_SC=22.0,        # EPDS span

    # ---- 9. nonspecific care (A5, calibrated to placebo arms) ------------
    K_CARE=0.019967,       # 1/h inpatient continuous-care effect
    K_NSP=0.0058718,        # 1/h outpatient trial-contact effect
    K_CBT=0.0030,        # 1/h psychotherapy effect at full engagement
    KCBT_ON=0.00248,     # 1/h psychotherapy engagement build-up (t1/2 280 h)

    # ---- 10. brexanolone PK (2-compartment, IV) --------------------------
    BRX_CL=70.0,         # L/h  (1.0 L/h/kg x 70 kg)
    BRX_V1=105.0,        # L    (1.5 L/kg)
    BRX_V2=300.0,        # L    (4.3 L/kg)
    BRX_Q=40.0,          # L/h  -> terminal t1/2 ~ 8.6 h

    # ---- 11. zuranolone PK (oral, 2-compartment) ------------------------
    ZUR_KA=0.45,         # 1/h
    ZUR_CL=12.0,         # L/h
    ZUR_V1=200.0,        # L
    ZUR_V2=100.0,        # L
    ZUR_Q=15.0,          # L/h -> terminal t1/2 ~ 19 h
    ZUR_EQ=0.10464,         # nM brain ALLO-eq per nM plasma zuranolone (A3)

    # ---- 12. sertraline PK (oral, 1-compartment) ------------------------
    SER_KA=0.60,         # 1/h
    SER_CL=96.0,         # L/h
    SER_V=3000.0,        # L -> t1/2 ~ 21.7 h

    # ---- 13. esketamine PK + plasticity surge ---------------------------
    ESK_CL=90.0,         # L/h
    ESK_V=250.0,         # L -> t1/2 ~ 1.9 h
    K_AMPA_IN=0.020,     # 1/h per 100 ng/mL esketamine
    K_AMPA_OUT=0.01925,  # 1/h AMPA surge decay (t1/2 36 h)

    # ---- 14. lactation transfer -----------------------------------------
    MP_BRX=1.50,         # milk:plasma ratio, brexanolone
    MP_ZUR=1.50,         # milk:plasma ratio, zuranolone
    MP_SER=1.80,         # milk:plasma ratio, sertraline
    MILK_INTAKE=150.0,   # mL/kg/day infant milk intake
    F_ORAL_NS=0.05,      # oral bioavailability of a neurosteroid in the infant
    V_INF=2.0,           # L/kg infant volume of distribution
    CL_INF=0.25,         # L/h/kg infant clearance (immature)
)

# state index map -----------------------------------------------------------
S = {n: i for i, n in enumerate("""PLAC P4 DHP ALLOP ALLOB E2 PCRH HCRH ACTH CORT
GRFN RD RG KCC2 MAOA FIVEHT AUTO BDNF SYN INFL KYNR SLP SDEBT EXC SYMP BONDI
BRX1 BRX2 ZURA ZUR1 ZUR2 SERA SERC ESKC AMPAS CBTP MILKD INFP""".split())}
NS = len(S)


# ---------------------------------------------------------------------------
# ALGEBRAIC HELPERS
# ---------------------------------------------------------------------------
def pot_tonic(pam, p):
    """Fractional potentiation of extrasynaptic (delta) tonic conductance."""
    x = max(pam, 0.0) ** p['H_PAM']
    return p['EMAX_PAM'] * x / (p['EC50_PAM'] ** p['H_PAM'] + x)


def pot_syn(pam, p):
    """Fractional potentiation of synaptic (gamma2) phasic inhibition."""
    x = max(pam, 0.0) ** p['H_PAM']
    return p['EMAX_SYN'] * x / (p['EC50_SYN'] ** p['H_PAM'] + x)


def rd_setpoint(pam, p):
    """Homeostatic delta-pool set-point: keeps the PRODUCT at G_TARGET."""
    return min(1.0, p['G_TARGET'] / (1.0 + pot_tonic(pam, p)))


def gref(p):
    """Non-pregnant reference effective inhibition (normalisation constant)."""
    pam = p['PAM_NP']
    rd = rd_setpoint(pam, p)
    rg = 1.0 + p['RG_GAIN'] * (1.0 - rd)
    gt = rd * (1.0 + pot_tonic(pam, p))
    gp = rg * (1.0 + pot_syn(pam, p))
    return p['W_TONIC'] * gt + p['W_PHASIC'] * gp


def wake_load(t, p, protect=False):
    """Infant-driven nocturnal waking load. t in h, t = 0 is delivery."""
    if t < 0.0:
        return p['WAKE_LATE_PREG']
    if protect:
        return p['WAKE_PROTECT']
    return p['WAKE_PP_FLOOR'] + p['WAKE_PP0'] * math.exp(-t / p['TAU_WAKE'])


def hamd(symp, p):
    return p['HAMD_0'] + p['HAMD_SC'] * symp


def epds(symp, p):
    return p['EPDS_0'] + p['EPDS_SC'] * symp


# ---------------------------------------------------------------------------
# THE ODE SYSTEM
# ---------------------------------------------------------------------------
def rhs(t, y, p, rx):
    """
    t   : hours, t = 0 is placental expulsion
    y   : state vector
    p   : parameter dict
    rx  : regimen dict with keys
            V        vulnerability gain (1.0 = population average)
            brx_rate function t -> ug/kg/h brexanolone infusion rate
            protect  function t -> bool  protected-sleep conditions
            care     function t -> bool  inpatient continuous care
            nsp      function t -> bool  outpatient trial contact
            cbt      function t -> bool  psychotherapy
            esk_inf  function t -> mg/h esketamine infusion rate
    """
    d = [0.0] * NS
    g = lambda n: y[S[n]]
    V = rx['V']

    # ---- 1. placenta and steroid trajectories ---------------------------
    plac = max(g('PLAC'), 0.0)
    d[S['PLAC']] = -p['KDEL'] * plac if t >= 0.0 else 0.0

    d[S['P4']] = p['KSYN_P4'] * plac + p['IN_P4'] - p['KEL_P4'] * g('P4')
    d[S['DHP']] = p['K5A'] * g('P4') - (p['K3A'] + p['KEL_DHP']) * g('DHP')
    d[S['ALLOP']] = (p['F_ALLO'] * p['K3A'] * g('DHP') + p['IN_ALLO']
                     - p['KEL_ALLO'] * g('ALLOP'))
    d[S['ALLOB']] = (p['KOUT_B'] * p['KP_ALLO'] * g('ALLOP') + p['IN_ALLOB']
                     - p['KOUT_B'] * g('ALLOB'))
    d[S['E2']] = p['KSYN_E2'] * plac + p['IN_E2'] - p['KEL_E2'] * g('E2')
    d[S['PCRH']] = p['KEL_PCRH'] * (plac - g('PCRH'))

    # ---- 2. drug PK ------------------------------------------------------
    brx_in = rx['brx_rate'](t) * p['WT']                       # ug/h
    d[S['BRX1']] = (brx_in
                    - (p['BRX_CL'] / p['BRX_V1']) * g('BRX1')
                    - (p['BRX_Q'] / p['BRX_V1']) * g('BRX1')
                    + (p['BRX_Q'] / p['BRX_V2']) * g('BRX2'))
    d[S['BRX2']] = ((p['BRX_Q'] / p['BRX_V1']) * g('BRX1')
                    - (p['BRX_Q'] / p['BRX_V2']) * g('BRX2'))
    d[S['ZURA']] = -p['ZUR_KA'] * g('ZURA')
    d[S['ZUR1']] = (p['ZUR_KA'] * g('ZURA')
                    - (p['ZUR_CL'] / p['ZUR_V1']) * g('ZUR1')
                    - (p['ZUR_Q'] / p['ZUR_V1']) * g('ZUR1')
                    + (p['ZUR_Q'] / p['ZUR_V2']) * g('ZUR2'))
    d[S['ZUR2']] = ((p['ZUR_Q'] / p['ZUR_V1']) * g('ZUR1')
                    - (p['ZUR_Q'] / p['ZUR_V2']) * g('ZUR2'))
    d[S['SERA']] = -p['SER_KA'] * g('SERA')
    d[S['SERC']] = p['SER_KA'] * g('SERA') - (p['SER_CL'] / p['SER_V']) * g('SERC')
    d[S['ESKC']] = rx['esk_inf'](t) - (p['ESK_CL'] / p['ESK_V']) * g('ESKC')

    # concentrations
    c_brx_ng = g('BRX1') / p['BRX_V1']                       # ng/mL
    c_brx_nM = c_brx_ng * 1000.0 / p['MW_ALLO']
    c_zur_ng = g('ZUR1') * 1e6 / p['ZUR_V1'] / 1000.0        # mg/L -> ng/mL
    c_zur_nM = c_zur_ng * 1000.0 / p['MW_ZUR']
    c_ser_ng = g('SERC') * 1e6 / p['SER_V'] / 1000.0         # ng/mL
    c_esk_ng = g('ESKC') * 1e6 / p['ESK_V'] / 1000.0         # ng/mL

    # ---- 3. total PAM in brain allopregnanolone equivalents -------------
    #  brexanolone IS allopregnanolone -> identical Kp, potency 1 (A2)
    pam = (g('ALLOB')
           + p['KP_ALLO'] * c_brx_nM
           + p['ZUR_EQ'] * p['KP_ALLO'] * c_zur_nM)

    # ---- 4. HPA ----------------------------------------------------------
    d[S['HCRH']] = p['KH'] * (1.0 / (1.0 + g('PCRH') / p['PC50']) - g('HCRH'))
    fb = 1.0 / (1.0 + (g('GRFN') * g('CORT') / p['CORT_FB']) ** p['N_FB'])
    drive_acth = (g('HCRH') + p['W_PCRH'] * g('PCRH')) * fb
    d[S['ACTH']] = p['KA_ACTH'] * (drive_acth - g('ACTH'))
    d[S['CORT']] = p['KC'] * g('ACTH') - p['KEL_CORT'] * g('CORT')
    grss = (1.0 / (1.0 + g('CORT') / p['GR50'])) / (1.0 + 0.40 * (V - 1.0))
    d[S['GRFN']] = p['KG'] * (grss - g('GRFN'))
    cort_rel = g('CORT') / p['CORT_NP']

    # ---- 5. receptor plasticity (THE SLOW ARM) --------------------------
    rdss = rd_setpoint(pam, p)
    kr_eff = (p['KR'] / (1.0 + p['V_KR'] * (V - 1.0))
              * (1.0 + p['KR_BOOST'] * pot_tonic(pam, p)))
    d[S['RD']] = kr_eff * (rdss - g('RD'))
    rgss = 1.0 + p['RG_GAIN'] * (1.0 - rdss)
    d[S['RG']] = kr_eff * (rgss - g('RG'))
    kcc2ss = ((1.0 / (1.0 + p['KCC2_CORT'] * max(0.0, cort_rel - 1.0) ** 1.2))
              * (0.70 + 0.30 * g('BDNF')))
    d[S['KCC2']] = p['KK'] * (kcc2ss - g('KCC2'))

    # effective inhibition
    g_tonic = g('RD') * (1.0 + pot_tonic(pam, p)) * math.sqrt(max(g('KCC2'), 1e-6))
    g_phasic = g('RG') * (1.0 + pot_syn(pam, p))
    g_eff = p['W_TONIC'] * g_tonic + p['W_PHASIC'] * g_phasic

    # ---- 6. monoamine / plasticity --------------------------------------
    maoass = 1.0 + p['F_MAOA'] * (1.0 - g('E2') / (g('E2') + p['E250']))
    d[S['MAOA']] = p['KM_MAOA'] * (maoass - g('MAOA'))
    occ = c_ser_ng / (c_ser_ng + p['EC50_SER'])
    ssri_gain = 1.0 + p['EMAX_SSRI'] * occ * (1.0 - p['AUTO_BRAKE'] * g('AUTO'))
    fivess = ((1.0 / max(g('MAOA'), 1e-6))
              / (1.0 + p['W_KYN_5HT'] * max(0.0, g('KYNR') - 1.0))
              * ssri_gain)
    d[S['FIVEHT']] = p['K5HT'] * (fivess - g('FIVEHT'))
    autoss = 1.0 / (1.0 + (occ / p['OCC50_AUTO']) ** 2)
    d[S['AUTO']] = p['KAUTO'] * (autoss - g('AUTO'))
    d[S['AMPAS']] = (p['K_AMPA_IN'] * c_esk_ng / 100.0
                     - p['K_AMPA_OUT'] * g('AMPAS'))
    bdnfss = (max(g('FIVEHT'), 1e-6) ** 0.60
              * (1.0 + p['W_AMPA'] * g('AMPAS'))
              * (1.0 + p['W_PAM_BDNF'] * pot_tonic(pam, p))
              / (1.0 + 0.50 * max(0.0, cort_rel - 1.0))
              / (1.0 + 0.30 * max(0.0, g('INFL') - p['INFL_BASE']) / 10.0))
    d[S['BDNF']] = p['KB_BDNF'] * (bdnfss - g('BDNF'))
    d[S['SYN']] = p['KSYN_PL'] * (g('BDNF') - g('SYN'))

    # ---- 7. inflammation / kynurenine ----------------------------------
    inflss = p['INFL_BASE'] * (1.0 + p['W_SD_INFL'] * g('SDEBT') / 15.0)
    d[S['INFL']] = p['KI_INFL'] * (inflss - g('INFL'))
    kynss = (1.0
             + p['W_INFL_KYN'] * max(0.0, g('INFL') - p['INFL_BASE']) / 10.0
             + p['W_CORT_KYN'] * max(0.0, cort_rel - 1.0))
    d[S['KYNR']] = p['KK_KYN'] * (kynss - g('KYNR'))

    # ---- 8. sleep --------------------------------------------------------
    wk = wake_load(t, p, protect=rx['protect'](t))
    slpss = p['SLP_MAX'] / (1.0 + p['A_WAKE'] * wk + p['A_SYMP'] * g('SYMP'))
    d[S['SLP']] = p['KS_SLP'] * (slpss - g('SLP'))
    d[S['SDEBT']] = ((p['SLP_NEED'] - g('SLP')) / 24.0
                     - p['KDEC_SD'] * g('SDEBT'))

    # ---- 9. excitatory load and symptoms -------------------------------
    drive = (1.0
             + p['W_SD'] * g('SDEBT') / 48.0
             + p['W_INFL'] * max(0.0, g('INFL') - p['INFL_BASE']) / 10.0
             + p['W_KYN'] * max(0.0, g('KYNR') - 1.0)
             + p['W_5HT'] * max(0.0, 1.0 / max(g('FIVEHT'), 1e-6) - 1.0)
             + p['W_HPA'] * max(0.0, 1.0 - cort_rel)
             + p['W_SELF'] * g('SYMP'))
    drive /= max(g('SYN'), 1e-6) ** p['W_SYN']
    d[S['EXC']] = p['KE_EXC'] * (drive * gref(p) / max(g_eff, 1e-6) - g('EXC'))

    thr = p['THR0'] / V ** p['V_THR']
    over = max(0.0, g('EXC') - thr)
    fsat = over / (1.0 + over)
    koff = (p['KOFF'] / V) * (1.0 + p['KFAST'] * max(0.0, thr - g('EXC')))
    k_care = p['K_CARE'] if rx['care'](t) else 0.0
    k_nsp = p['K_NSP'] if rx['nsp'](t) else 0.0
    d[S['CBTP']] = (p['KCBT_ON'] * (1.0 - g('CBTP'))) if rx['cbt'](t) else 0.0
    k_cbt = p['K_CBT'] * g('CBTP')
    d[S['SYMP']] = (p['KON'] * V * fsat * (p['SMAX'] - g('SYMP'))
                    - (koff + k_care + k_nsp + k_cbt)
                    * max(0.0, g('SYMP') - p['SMIN']))
    d[S['BONDI']] = p['KBOND'] * (p['W_BOND'] * g('SYMP') - g('BONDI'))

    # ---- 10. lactation transfer -----------------------------------------
    c_milk = (p['MP_BRX'] * c_brx_ng + p['MP_ZUR'] * c_zur_ng
              + p['MP_SER'] * c_ser_ng)                        # ng/mL
    inf_rate = c_milk * p['MILK_INTAKE'] / 24.0 / 1000.0       # ug/kg/h
    d[S['MILKD']] = inf_rate
    d[S['INFP']] = (p['F_ORAL_NS'] * inf_rate / p['V_INF'] * 1000.0
                    / p['MW_ALLO']
                    - (p['CL_INF'] / p['V_INF']) * g('INFP'))
    return d


# ---------------------------------------------------------------------------
# INTEGRATION
# ---------------------------------------------------------------------------
def rk4(t0, t1, y, p, rx, h=0.05, events=None, record=None, out=None,
        record_every=0.5):
    """Fixed-step RK4 with discrete events (dosing, boluses)."""
    t = t0
    n = int(round((t1 - t0) / h))
    stride = max(1, int(round(record_every / h)))
    for i in range(n):
        if events:
            for ev_t, fn in events:
                if t <= ev_t < t + h:
                    fn(y)
        if record and out is not None and i % stride == 0:
            record(t, y, out)
        k1 = rhs(t, y, p, rx)
        y2 = [y[j] + 0.5 * h * k1[j] for j in range(NS)]
        k2 = rhs(t + 0.5 * h, y2, p, rx)
        y3 = [y[j] + 0.5 * h * k2[j] for j in range(NS)]
        k3 = rhs(t + 0.5 * h, y3, p, rx)
        y4 = [y[j] + h * k3[j] for j in range(NS)]
        k4 = rhs(t + h, y4, p, rx)
        for j in range(NS):
            y[j] += h / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
            if j not in (S['MILKD'],):
                y[j] = max(y[j], 0.0)
        t += h
    return t, y


NO = lambda t: False
ZERO = lambda t: 0.0


def base_rx(V=1.0, **kw):
    rx = dict(V=V, brx_rate=ZERO, protect=NO, care=NO, nsp=NO, cbt=NO,
              esk_inf=ZERO)
    rx.update(kw)
    return rx


_SS_CACHE = {}


def pregnancy_steady_state(p, V=1.0):
    """Burn-in: 4000 h with an intact placenta, then return the state."""
    ck = (V, tuple(sorted(p.items())))
    if ck in _SS_CACHE:
        return list(_SS_CACHE[ck])
    y = [0.0] * NS
    y[S['PLAC']] = 1.0
    y[S['RD']] = 0.5
    y[S['RG']] = 1.0
    y[S['KCC2']] = 1.0
    y[S['MAOA']] = 1.0
    y[S['FIVEHT']] = 1.0
    y[S['AUTO']] = 1.0
    y[S['BDNF']] = 1.0
    y[S['SYN']] = 1.0
    y[S['INFL']] = p['INFL_BASE']
    y[S['KYNR']] = 1.0
    y[S['SLP']] = 7.0
    y[S['GRFN']] = 0.7
    y[S['CORT']] = 15.0
    y[S['ACTH']] = 1.0
    y[S['EXC']] = 1.0
    rx = base_rx(V=V)
    rk4(-4000.0, -48.0, y, p, rx, h=0.5)
    _SS_CACHE[ck] = list(y)
    return list(y)


def deliver(y, p):
    """Discrete consequence of delivery that is not a rate: the IL-6 surge."""
    y[S['INFL']] += p['DEL_INFL']


def simulate(p, V=1.7, t_end=1080.0, enrol=504.0,
             arm='none', enrol_hamd_target=None, h=0.10, want=None,
             care_until=None, nsp_after_enrol=True, protect_until=None,
             cbt=False):
    """
    Run pregnancy -> delivery -> postpartum, apply a treatment arm at `enrol`,
    and return a dict of trajectories at requested times (hours postpartum).

    arm: 'none' | 'placebo_in' | 'placebo_out' | 'brx60' | 'brx90'
         | 'zur50' | 'zur30' | 'zur50_3d' | 'sertraline' | 'esketamine'
         | 'sleep' | 'cbt'
    """
    y = pregnancy_steady_state(p, V=V)
    events = [(0.0, lambda yy: deliver(yy, p))]

    # --- regimen construction -------------------------------------------
    brx = ZERO
    esk = ZERO
    care = NO
    nsp = NO
    protect = NO
    cbt_f = NO

    if arm in ('placebo_in', 'brx60', 'brx90'):
        # inpatient 60-hour infusion setting
        t0, t1 = enrol, enrol + 60.0
        care = lambda t, a=t0, b=t1: a <= t < b
        protect = lambda t, a=t0, b=t1: a <= t < b
        nsp = lambda t, a=t0: t >= a
        if arm in ('brx60', 'brx90'):
            top = 60.0 if arm == 'brx60' else 90.0

            def brx(t, a=t0, top=top):
                dt = t - a
                if dt < 0 or dt >= 60.0:
                    return 0.0
                if dt < 4.0:
                    return 30.0
                if dt < 24.0:
                    return min(60.0, top)
                if dt < 52.0:
                    return top
                return 30.0
    elif arm in ('placebo_out', 'zur50', 'zur30', 'zur50_3d', 'sertraline',
                 'esketamine', 'cbt', 'sleep'):
        nsp = lambda t, a=enrol: t >= a
        if arm == 'sleep':
            protect = lambda t, a=enrol: t >= a
        if arm == 'cbt' or cbt:
            cbt_f = lambda t, a=enrol: t >= a
        if arm == 'esketamine':
            def esk(t, a=enrol):
                return (0.25 * p['WT'] / 0.667) if a <= t < a + 0.667 else 0.0

    if care_until is not None:
        care = lambda t, a=enrol, b=care_until: a <= t < b
    if protect_until is not None:
        protect = lambda t, a=enrol, b=protect_until: a <= t < b
    if not nsp_after_enrol:
        nsp = NO

    rx = base_rx(V=V, brx_rate=brx, protect=protect, care=care, nsp=nsp,
                 cbt=cbt_f, esk_inf=esk)

    # --- oral dosing events ---------------------------------------------
    if arm in ('zur50', 'zur30', 'zur50_3d'):
        dose = 50.0 if arm != 'zur30' else 30.0
        ndays = 14 if arm != 'zur50_3d' else 3
        for k in range(ndays):
            events.append((enrol + 24.0 * k,
                           lambda yy, dd=dose: yy.__setitem__(
                               S['ZURA'], yy[S['ZURA']] + dd)))
    if arm == 'sertraline':
        for k in range(int((t_end - enrol) // 24) + 1):
            dose = 50.0 if k < 7 else 100.0
            events.append((enrol + 24.0 * k,
                           lambda yy, dd=dose: yy.__setitem__(
                               S['SERA'], yy[S['SERA']] + dd)))
    events.sort(key=lambda e: e[0])

    # --- record ----------------------------------------------------------
    want = sorted(want or [])
    out = {}
    grid = {}

    def record(t, yy, o):
        grid[round(t, 4)] = list(yy)

    t = -48.0
    rk4(t, t_end, y, p, rx, h=h, events=events, record=record, out=out)

    keys = sorted(grid)

    def at(tq):
        best = min(keys, key=lambda k: abs(k - tq))
        return grid[best]

    res = dict(states={}, series=[])
    for tq in want:
        res['states'][tq] = at(tq)
    res['at'] = at
    res['grid'] = grid
    res['keys'] = keys
    return res


# ---------------------------------------------------------------------------
# REPORTING HELPERS
# ---------------------------------------------------------------------------
def line(ch='-', n=79):
    return ch * n


def head(title):
    print()
    print(line('='))
    print(title)
    print(line('='))


def gv(state, name):
    return state[S[name]]


def main():
    p = dict(P)
    print(line('='))
    print(" POSTPARTUM DEPRESSION QSP MODEL — REFERENCE IMPLEMENTATION OUTPUT")
    print(" 38-ODE system · pure-stdlib RK4 (h = 0.1 h) · t = 0 is delivery")
    print(line('='))

    # =====================================================================
    head("TABLE 1 · The two arms of the product, and why a trough exists")
    # =====================================================================
    print("Potentiation is Emax*PAM^h/(EC50^h+PAM^h) with Emax=%.2f, EC50=%.0f nM,"
          " h=%.1f" % (p['EMAX_PAM'], p['EC50_PAM'], p['H_PAM']))
    print("The delta-pool set-point is R_ss = G_target/(1+potentiation), so the")
    print("PRODUCT is held at G_target = %.2f in ANY steady state.\n"
          % p['G_TARGET'])
    print("%-28s %10s %12s %12s %12s"
          % ("state", "PAM(nM)", "potentiation", "R_delta,ss", "G_tonic,ss"))
    print(line())
    for nm, pam in [("non-pregnant", 6.0), ("2nd trimester", 60.0),
                    ("term (3rd trimester)", 160.0),
                    ("day 3 postpartum", 4.8), ("day 5 postpartum", 4.0)]:
        rd = rd_setpoint(pam, p)
        print("%-28s %10.1f %12.3f %12.3f %12.3f"
              % (nm, pam, pot_tonic(pam, p), rd, rd * (1 + pot_tonic(pam, p))))
    print()
    print("The steady states are all identical (%.2f) — the disease is NOT a"
          % p['G_TARGET'])
    print("steady state.  It is the TRANSIENT between two of them:")
    rd_term = rd_setpoint(160.0, p)
    pam_pp = 4.8
    g_pp = rd_term * (1 + pot_tonic(pam_pp, p))
    print("   R_delta is still %.3f (the term value) when PAM has already"
          % rd_term)
    print("   reached %.1f nM  ->  G_tonic = %.3f  =  a %.0f%% LOSS of tonic"
          % (pam_pp, g_pp, 100 * (1 - g_pp / p['G_TARGET'])))
    print("   inhibition that no steady-state analysis would ever show.")
    print()
    kr = p['KR']
    for frac, lbl in [(0.25, "25 %"), (0.10, "10 %"), (0.05, "5 %")]:
        rd_target = (1.0 - frac) * p['G_TARGET'] / (1 + pot_tonic(pam_pp, p))
        rd_inf = rd_setpoint(pam_pp, p)
        if rd_target >= rd_inf:
            print("   deficit never as small as %s" % lbl)
            continue
        frac_way = (rd_target - rd_term) / (rd_inf - rd_term)
        tt = -math.log(1 - frac_way) / kr
        print("   tonic-inhibition deficit stays worse than %-5s for %6.0f h "
              "= %4.1f d" % (lbl, tt, tt / 24.0))
    print()
    print("=> the model DERIVES a 2-3 week window of reduced tonic inhibition")
    print("   from two time constants (ligand t1/2 %.0f h, receptor t1/2 %.0f h)"
          % (0.693 / p['KEL_ALLO'], 0.693 / p['KR']))
    print("   and nothing else.  Compare with the observed clustering of PPD")
    print("   onset in the first postpartum weeks (PMID 23487258).  The width")
    print("   is set by the RECEPTOR constant, which has no human measurement")
    print("   (assumption A8) — so treat the width as an order of magnitude,")
    print("   not a prediction.")

    # =====================================================================
    head("TABLE 2 · Simulated endocrine trajectory across delivery")
    # =====================================================================
    r = simulate(p, V=1.0, t_end=1080.0, arm='none',
                 want=[-24, 0, 6, 12, 24, 48, 72, 120, 240, 504, 1080])
    print("%8s %9s %9s %9s %9s %9s %8s %8s %8s"
          % ("t (h)", "P4 nM", "ALLOp nM", "ALLOb nM", "E2 nM",
             "R_delta", "G_tonic", "CORT", "HCRH"))
    print(line())
    for tq in [-24, 0, 6, 12, 24, 48, 72, 120, 240, 504, 1080]:
        st = r['states'][tq]
        pam = gv(st, 'ALLOB')
        gt = gv(st, 'RD') * (1 + pot_tonic(pam, p)) * math.sqrt(gv(st, 'KCC2'))
        print("%8.0f %9.1f %9.2f %9.2f %9.3f %9.3f %8.3f %8.1f %8.3f"
              % (tq, gv(st, 'P4'), gv(st, 'ALLOP'), pam, gv(st, 'E2'),
                 gv(st, 'RD'), gt, gv(st, 'CORT'), gv(st, 'HCRH')))
    st0 = r['states'][-24]
    st72 = r['states'][72]
    print()
    print("plasma allopregnanolone falls %.0f-fold (%.1f -> %.2f nM) in 72 h;"
          % (gv(st0, 'ALLOP') / gv(st72, 'ALLOP'),
             gv(st0, 'ALLOP'), gv(st72, 'ALLOP')))
    print("estradiol falls %.0f-fold; the delta-receptor pool has moved %.1f %%"
          % (gv(st0, 'E2') / max(gv(st72, 'E2'), 1e-9),
             100 * (gv(st72, 'RD') - gv(st0, 'RD'))
             / max(1e-9, (0.964 - gv(st0, 'RD')))))
    print("of the way back to its non-pregnant value.  THAT is the mismatch.")

    # =====================================================================
    head("TABLE 3 · Same trajectory, three vulnerabilities (no treatment)")
    # =====================================================================
    print("Identical hormones in every row.  Only V differs.")
    print("V scales symptom accrual gain, narrows the E/I threshold")
    print("(THR = %.1f/V^%.1f) and slows delta-pool recovery.\n"
          % (p['THR0'], p['V_THR']))
    print("%-26s %7s %7s %8s %8s %8s %8s %8s"
          % ("phenotype", "V", "THR", "peakEXC", "peakS", "HAMD", "EPDS",
             "d42 HAMD"))
    print(line())
    peaks = {}
    for lbl, V in [("resilient (support+)", 0.80),
                   ("population average", 1.00),
                   ("at risk (prior MDD)", 1.40),
                   ("high risk (prior PPD)", 1.70),
                   ("very high risk", 2.00)]:
        rr = simulate(p, V=V, t_end=1080.0, arm='none',
                      want=list(range(0, 1081, 6)))
        keys = [k for k in rr['keys'] if k >= 0]
        mx_s, mx_e, t_mx = 0.0, 0.0, 0.0
        for k in keys:
            stt = rr['grid'][k]
            if stt[S['SYMP']] > mx_s:
                mx_s, t_mx = stt[S['SYMP']], k
            mx_e = max(mx_e, stt[S['EXC']])
        st42 = rr['at'](42 * 24.0)
        peaks[V] = (mx_s, t_mx, rr)
        print("%-26s %7.2f %7.2f %8.2f %8.3f %8.1f %8.1f %8.1f"
              % (lbl, V, p['THR0'] / V ** p['V_THR'], mx_e, mx_s,
                 hamd(mx_s, p), epds(mx_s, p), hamd(gv(st42, 'SYMP'), p)))
    print()
    print("peak symptom load occurs at t = %.0f h = %.1f d postpartum for the"
          % (peaks[1.70][1], peaks[1.70][1] / 24.0))
    print("high-risk phenotype.  The average-V woman gets a self-limited dip")
    print("(EPDS %.1f, below the 12/13 screening cut-off) — the model's"
          % epds(peaks[1.00][0], p))
    print("account of 'baby blues' — while the high-risk woman latches into a")
    print("sustained episode (HAM-D %.1f).  ONE parameter separates them."
          % hamd(peaks[1.70][0], p))

    # =====================================================================
    head("TABLE 4 · Brexanolone: PK, saturation, and the flat 60->90 arm")
    # =====================================================================
    print("%-14s %10s %10s %12s %12s %10s"
          % ("rate ug/kg/h", "Css ng/mL", "Css nM", "brain nM-eq",
             "potentiation", "vs 30"))
    print(line())
    base = None
    for rate in [30.0, 60.0, 90.0, 120.0]:
        css_ng = rate * p['WT'] / p['BRX_CL']
        css_nM = css_ng * 1000.0 / p['MW_ALLO']
        brain = p['KP_ALLO'] * css_nM
        po = pot_tonic(brain, p)
        base = base or po
        print("%-14.0f %10.1f %10.1f %12.0f %12.3f %10.2f"
              % (rate, css_ng, css_nM, brain, po, po / base))
    po60 = pot_tonic(p['KP_ALLO'] * (60 * p['WT'] / p['BRX_CL']) * 1000
                     / p['MW_ALLO'], p)
    po90 = pot_tonic(p['KP_ALLO'] * (90 * p['WT'] / p['BRX_CL']) * 1000
                     / p['MW_ALLO'], p)
    print()
    print("Going from 60 to 90 ug/kg/h raises steady-state exposure by 50 %")
    print("but potentiation by only %.1f %% — the Hill function is already"
          % (100 * (po90 / po60 - 1)))
    print("saturated.  Site-2 direct activation (sedation) has no such ceiling")
    print("in this range.  A flat 60-vs-90 efficacy comparison with MORE")
    print("sedation at the higher rate is therefore the PREDICTION, and it is")
    print("what phase 3 study 1 reported (-19.5 for BRX60 vs -17.7 for BRX90).")
    print()
    print("Steady-state brexanolone exposure at 90 ug/kg/h is %.0f nM plasma,"
          % (90 * p['WT'] / p['BRX_CL'] * 1000 / p['MW_ALLO']))
    print("i.e. %.1fx the third-trimester allopregnanolone level the model uses"
          % ((90 * p['WT'] / p['BRX_CL'] * 1000 / p['MW_ALLO']) / 80.0))
    print("(80 nM).  The infusion does not 'restore pregnancy levels' — it")
    print("overshoots them severalfold, which is exactly why sedation, not")
    print("efficacy, is the dose-limiting toxicity.")

    # =====================================================================
    head("TABLE 5 · Treatment arms vs the published trial endpoints")
    # =====================================================================
    V = 1.70
    arms = [('placebo_in', 'placebo (inpatient, 60 h)'),
            ('brx60', 'brexanolone 60 ug/kg/h'),
            ('brx90', 'brexanolone 90 ug/kg/h'),
            ('placebo_out', 'placebo (outpatient)'),
            ('zur30', 'zuranolone 30 mg x 14 d'),
            ('zur50', 'zuranolone 50 mg x 14 d'),
            ('zur50_3d', 'zuranolone 50 mg x 3 d'),
            ('sertraline', 'sertraline 50->100 mg'),
            ('esketamine', 'esketamine 0.25 mg/kg once'),
            ('sleep', 'protected sleep only'),
            ('cbt', 'CBT/IPT only')]
    enrol = 504.0
    store = {}
    print("HAM-D17 change from enrolment (day 21 postpartum, V = %.2f)\n" % V)
    print("%-30s %7s %7s %7s %7s %7s %7s"
          % ("arm", "base", "60 h", "d 3", "d 15", "d 28", "d 45"))
    print(line())
    for arm, lbl in arms:
        rr = simulate(p, V=V, t_end=enrol + 45 * 24 + 1, enrol=enrol, arm=arm,
                      want=[])
        b = hamd(gv(rr['at'](enrol), 'SYMP'), p)
        row = [b]
        for dt in [60.0, 72.0, 15 * 24.0, 28 * 24.0, 45 * 24.0]:
            row.append(hamd(gv(rr['at'](enrol + dt), 'SYMP'), p) - b)
        store[arm] = (rr, b, row)
        print("%-30s %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f"
              % (lbl, row[0], row[1], row[2], row[3], row[4], row[5]))
    print()
    print("OBSERVED (published, see header):")
    print("  brexanolone phase 3 study 1, 60 h : BRX60 -19.5 · BRX90 -17.7 ·"
          " placebo -14.0")
    print("  brexanolone phase 2, 60 h         : BRX -21.0 · placebo -8.8")
    print("  SKYLARK, day 15                   : ZUR50 -15.6 · placebo -11.6")
    print("  ROBIN, day 15                     : ZUR30 -17.8 · placebo -13.6")
    print()
    d_brx = store['brx60'][2][1] - store['placebo_in'][2][1]
    d_z50 = store['zur50'][2][3] - store['placebo_out'][2][3]
    print("MODEL drug-placebo differences:  brexanolone 60 h %+.1f  (obs -5.5)"
          % d_brx)
    print("                                zuranolone d15   %+.1f  (obs -4.0)"
          % d_z50)
    print()
    print("WHERE THIS FAILS — state it plainly:")
    print("  The model reproduces (a) the enrolment severity, (b) the 60-hour")
    print("  brexanolone response, (c) the day-15 zuranolone response and")
    print("  (d) the magnitude of both placebo responses.  It does NOT")
    print("  reproduce the MAINTENANCE of drug-placebo separation at days")
    print("  28-45: every arm converges to about HAM-D %.0f, whereas the trials"
          % (store['placebo_out'][1] + store['placebo_out'][2][5]))
    print("  kept a small separation to day 30 (brexanolone) and day 45")
    print("  (zuranolone).  The late-time attractor of this system is a")
    print("  partially-recovered state that all arms fall into once the")
    print("  receptor trough has closed, and no drug parameter changes where")
    print("  that attractor sits — only how fast a patient reaches it.")
    print("  Table 12 shows that the two obvious mechanistic repairs fail, and")
    print("  why the failure is structural rather than a fitting problem.")

    # =====================================================================
    head("TABLE 6 · The bridge test: is durability the receptor, or the drug?")
    # =====================================================================
    print("Brexanolone stops at 60 h.  Zuranolone stops at day 14.  Both keep")
    print("working.  The model says the drug does not hold the improvement —")
    print("it holds the patient until R_delta has re-expanded, and the")
    print("symptom-sleep loop then keeps her in the recovered state.")
    print("So switch OFF receptor recovery (KR = 0) and re-run.\n")
    print("%-30s %10s %10s %10s"
          % ("arm", "d 3", "d 15", "d 45"))
    print(line())
    p_norec = dict(p)
    p_norec['KR'] = 0.0
    for arm, lbl in [('brx60', 'brexanolone, KR normal'),
                     ('zur50', 'zuranolone 50 mg, KR normal')]:
        rr, b, row = store[arm]
        print("%-30s %10.1f %10.1f %10.1f" % (lbl, row[2], row[3], row[5]))
    for arm, lbl in [('brx60', 'brexanolone, KR = 0'),
                     ('zur50', 'zuranolone 50 mg, KR = 0')]:
        rr = simulate(p_norec, V=V, t_end=enrol + 45 * 24 + 1, enrol=enrol,
                      arm=arm, want=[])
        b = hamd(gv(rr['at'](enrol), 'SYMP'), p)
        row = [hamd(gv(rr['at'](enrol + dt), 'SYMP'), p) - b
               for dt in [72.0, 15 * 24.0, 45 * 24.0]]
        print("%-30s %10.1f %10.1f %10.1f" % (lbl, row[0], row[1], row[2]))
    print()
    print("With the receptor frozen, the ON-DRUG response survives (the ligand")
    print("arm still works) but the day-45 benefit is 3-4 HAM-D points worse in")
    print("BOTH arms.  No parameter of the symptom equations was touched, so")
    print("this is a structural consequence: part of what looks like drug")
    print("durability is receptor recovery that would have happened anyway.")
    print("It is NOT the whole story — see the failure noted under Table 5.")
    print()
    print("Course length at the same daily dose (zuranolone 50 mg):")
    for arm, lbl in [('zur50_3d', '3-day course'), ('zur50', '14-day course')]:
        rr, b, row = store[arm]
        print("  %-14s day 15 %+6.1f   day 45 %+6.1f" % (lbl, row[3], row[5]))
    print()
    print("NOTE — the course-length prediction does NOT come out: the 3-day")
    print("course is worse at day 15 (it has stopped) but indistinguishable at")
    print("day 45, because by then every arm has reached the same attractor.")
    print("So this model does not support 'shorter courses relapse'.  It")
    print("supports only the weaker claim that a course must still be running")
    print("to hold an acute advantage.  A model that reproduced maintained")
    print("separation would be needed to say anything about course length, and")
    print("this is not that model.")

    # =====================================================================
    head("TABLE 7 · Why the placebo response is so large, and so variable")
    # =====================================================================
    print("The placebo arms of these trials improved by 8.8 to 14.0 HAM-D")
    print("points.  In this model that is not noise: an enrolled woman is")
    print("sitting on a receptor pool that is still re-expanding, and trial")
    print("participation itself changes her sleep.  Both are dated quantities,")
    print("so the model predicts placebo response should SHRINK with later")
    print("enrolment.\n")
    print("%-16s %10s %10s %12s %12s"
          % ("enrolment", "R_delta", "G_tonic", "placebo 60h", "placebo d15"))
    print(line())
    for day in [7, 14, 21, 42, 90, 180]:
        e = day * 24.0
        rr_in = simulate(p, V=V, t_end=e + 45 * 24 + 1, enrol=e,
                         arm='placebo_in', want=[])
        rr_out = simulate(p, V=V, t_end=e + 45 * 24 + 1, enrol=e,
                          arm='placebo_out', want=[])
        st = rr_in['at'](e)
        gt = (gv(st, 'RD') * (1 + pot_tonic(gv(st, 'ALLOB'), p))
              * math.sqrt(gv(st, 'KCC2')))
        b_in = hamd(gv(rr_in['at'](e), 'SYMP'), p)
        b_out = hamd(gv(rr_out['at'](e), 'SYMP'), p)
        print("%-16s %10.3f %10.3f %12.1f %12.1f"
              % ("day %d" % day, gv(st, 'RD'), gt,
                 hamd(gv(rr_in['at'](e + 60.0), 'SYMP'), p) - b_in,
                 hamd(gv(rr_out['at'](e + 15 * 24.0), 'SYMP'), p) - b_out))
    print()
    print("TWO readings, and only one of them supports the model:")
    print("  SUPPORTED — the inpatient 60-hour column is larger than the")
    print("  outpatient day-15 column at every enrolment day, and for a")
    print("  mechanistic reason: continuous nursing care protects sleep, which")
    print("  is a state variable here.  The brexanolone trials were inpatient")
    print("  and their placebo arms improved -14.0 in 60 h; the zuranolone")
    print("  trials were outpatient and theirs took 15 days to reach -11.6.")
    print("  The model reproduces that ordering without being told about it.")
    print("  NOT SUPPORTED — the predicted placebo response is nearly FLAT in")
    print("  enrolment day (about 1 HAM-D point across day 7 to day 180), even")
    print("  though the delta-pool column moves from 0.66 to 0.97.  The")
    print("  calibrated care term dominates the mechanistic one.  So the model")
    print("  does not, after calibration, predict a useful enrolment-timing")
    print("  effect, and the earlier hypothesis that receptor recovery explains")
    print("  much of the placebo response is NOT borne out at these weights.")

    # =====================================================================
    head("TABLE 8 · Onset kinetics separate the mechanisms, not the efficacy")
    # =====================================================================
    print("%-30s %8s %8s %8s %8s %8s"
          % ("arm", "d 1", "d 3", "d 7", "d 15", "d 45"))
    print(line())
    for arm, lbl in [('brx60', 'brexanolone 60 (GABA_A)'),
                     ('zur50', 'zuranolone 50 (GABA_A)'),
                     ('esketamine', 'esketamine (glutamate)'),
                     ('sertraline', 'sertraline (monoamine)'),
                     ('placebo_out', 'placebo (outpatient)')]:
        rr, b, _ = store[arm]
        row = [hamd(gv(rr['at'](enrol + dt * 24.0), 'SYMP'), p) - b
               for dt in [1, 3, 7, 15, 45]]
        print("%-30s %8.1f %8.1f %8.1f %8.1f %8.1f" % tuple([lbl] + row))
    print()
    print("Sertraline is the one arm that keeps improving past day 28 in this")
    print("model, because its slow plasticity arm keeps pushing after the")
    print("receptor trough has closed.  That is a model artefact worth naming:")
    print("it comes from SYN entering the excitatory load as a divisor with no")
    print("ceiling, and it should not be read as a claim that sertraline beats")
    print("a neurosteroid at week 6.")
    print()
    print("Nothing about sertraline was made weak: its SERT occupancy is")
    st_s = store['sertraline'][0]['at'](enrol + 15 * 24.0)
    c_ser = gv(st_s, 'SERC') * 1e3 / p['SER_V']
    print("immediate (%.0f %% at day 15) and its final effect is competitive."
          % (100 * c_ser / (c_ser + p['EC50_SER'])))
    print("The delay comes from two states in series — 5-HT1A autoreceptor")
    print("desensitisation (t1/2 %.0f d) and structural plasticity (t1/2 %.0f d)"
          % (0.693 / p['KAUTO'] / 24, 0.693 / p['KSYN_PL'] / 24))
    print("— which is why a GABA_A PAM and an SSRI can converge by week 6 and")
    print("still be unrecognisably different on day 3.")

    # =====================================================================
    head("TABLE 9 · Maternal versus infant exposure during lactation")
    # =====================================================================
    print("%-26s %10s %10s %10s %10s"
          % ("regimen", "Cmax mat", "milk ng/mL", "RID %", "infant nM"))
    print(line())
    for arm, lbl, mdose in [('brx90', 'brexanolone 90 (60 h)', None),
                            ('zur50', 'zuranolone 50 mg/d', 50.0),
                            ('sertraline', 'sertraline 100 mg/d', 100.0)]:
        rr, b, _ = store[arm]
        keys = [k for k in rr['keys'] if enrol <= k <= enrol + 15 * 24]
        cmax, milk_max, inf_max = 0.0, 0.0, 0.0
        d0, d1 = None, None
        for k in keys:
            st = rr['grid'][k]
            cb = gv(st, 'BRX1') / p['BRX_V1']
            cz = gv(st, 'ZUR1') * 1e3 / p['ZUR_V1']
            cs = gv(st, 'SERC') * 1e3 / p['SER_V']
            c = max(cb, cz, cs)
            mk = p['MP_BRX'] * cb + p['MP_ZUR'] * cz + p['MP_SER'] * cs
            cmax = max(cmax, c)
            milk_max = max(milk_max, mk)
            inf_max = max(inf_max, gv(st, 'INFP'))
            if abs(k - (enrol + 10 * 24)) < 0.03:
                d0 = gv(st, 'MILKD')
            if abs(k - (enrol + 11 * 24)) < 0.03:
                d1 = gv(st, 'MILKD')
        if arm == 'brx90':
            st_e = rr['at'](enrol + 60.0)
            tot_inf = gv(st_e, 'MILKD')                      # ug/kg over 60 h
            mat_tot = 90.0 * 60.0                            # ug/kg over 60 h
            rid = 100.0 * tot_inf / mat_tot
        else:
            daily_inf = (d1 - d0) if (d0 is not None and d1 is not None) else 0.0
            rid = 100.0 * (daily_inf / 1000.0) / (mdose / p['WT'])
        print("%-26s %10.1f %10.1f %10.2f %10.4f"
              % (lbl, cmax, milk_max, rid, inf_max))
    print()
    print("The relative infant dose is what matters, not the maternal dose.")
    print("The 60-hour brexanolone infusion delivers well under 1 %% of the")
    print("maternal dose per kg, against a measured brexanolone milk study")
    print("(PMID 35869362) that this column should be checked against; the")
    print("14-day oral zuranolone course delivers several-fold more because it")
    print("is given for 14 days rather than 60 hours, and still sits under the")
    print("conventional 10 %% concern threshold.")
    print("An orally dosed neurosteroid is additionally crippled by first-pass")
    print("metabolism (F_oral = %.2f assumed here, with NO neonatal"
          % p['F_ORAL_NS'])
    print("measurement behind it), so the infant plasma concentrations above")
    print("are ~2 orders of magnitude below the model's own potentiation EC50")
    print("of %.0f nM.  Maternal sedation and infant sedation are different"
          % p['EC50_PAM'])
    print("questions and this model separates them explicitly rather than")
    print("reasoning from the maternal dose.")

    # =====================================================================
    head("TABLE 10 · Non-pharmacological arms and the vicious cycle")
    # =====================================================================
    print("The symptom -> insomnia -> sleep debt -> excitatory load -> symptom")
    print("loop is closed in this model.  Its loop gain is what makes a")
    print("60-hour intervention durable, and it is also a target itself.\n")
    print("%-30s %8s %8s %8s %8s"
          % ("arm", "d 3", "d 15", "d 45", "sleep h"))
    print(line())
    for arm, lbl in [('placebo_out', 'no intervention (contact only)'),
                     ('sleep', 'protected sleep block'),
                     ('cbt', 'CBT/IPT'),
                     ('zur50', 'zuranolone 50 mg')]:
        rr, b, _ = store[arm]
        row = [hamd(gv(rr['at'](enrol + dt * 24.0), 'SYMP'), p) - b
               for dt in [3, 15, 45]]
        slp = gv(rr['at'](enrol + 15 * 24.0), 'SLP')
        print("%-30s %8.1f %8.1f %8.1f %8.2f" % tuple([lbl] + row + [slp]))
    print()
    print("Protected sleep is not a placebo control in this model — it is an")
    print("active arm acting on a real state variable, which is exactly why the")
    print("inpatient placebo arm of the brexanolone trials was so effective and")
    print("why 'placebo-controlled' is a subtler claim in PPD than elsewhere.")

    # =====================================================================
    head("TABLE 11 · Sensitivity of the day-15 zuranolone effect")
    # =====================================================================
    print("one parameter at a time, +/- 30 %, effect on model day-15 HAM-D")
    print("change in the zuranolone 50 mg arm (reference %.1f)\n"
          % store['zur50'][2][3])
    ref = store['zur50'][2][3]
    rows = []
    for key in ['KR', 'EC50_PAM', 'EMAX_PAM', 'ZUR_EQ', 'THR0', 'KON', 'KOFF',
                'KFAST', 'W_SD', 'A_SYMP', 'KDEC_SD', 'K_NSP', 'V_KR']:
        vals = []
        for mult in (0.7, 1.3):
            pp = dict(p)
            pp[key] = p[key] * mult
            rr = simulate(pp, V=V, t_end=enrol + 16 * 24, enrol=enrol,
                          arm='zur50', want=[])
            b = hamd(gv(rr['at'](enrol), 'SYMP'), pp)
            vals.append(hamd(gv(rr['at'](enrol + 15 * 24.0), 'SYMP'), pp) - b)
        rows.append((key, vals[0], vals[1], max(abs(vals[0] - ref),
                                                abs(vals[1] - ref))))
    rows.sort(key=lambda r: -r[3])
    print("%-14s %10s %10s %10s" % ("parameter", "-30 %", "+30 %", "max |d|"))
    print(line())
    for k, a, b2, m in rows:
        print("%-14s %10.1f %10.1f %10.2f" % (k, a, b2, m))
    print()
    print("Read this ranking against the model's own story, because it does")
    print("not fully agree with it:")
    print("  * THR0 is first, as claimed — the E/I threshold is the disease")
    print("    boundary in this model and the endpoint moves 2.6 points with it.")
    print("  * K_NSP is SECOND, and K_NSP is the fitted, non-mechanistic")
    print("    care/expectancy term (assumption A5).  That is a warning label:")
    print("    a large share of this endpoint is calibration, not mechanism.")
    print("  * KDEC_SD and W_SD are third and fourth, i.e. the SLEEP LOOP")
    print("    matters more here than anything to do with the drug.")
    print("  * the drug's own potency parameters (ZUR_EQ, EC50_PAM, EMAX_PAM)")
    print("    sit mid-table, because potentiation is saturated at therapeutic")
    print("    exposure (Table 4) — that part of the story does hold.")
    print("  * KR ranks LOW (0.62) even though it sets the trough width in")
    print("    Table 1.  The two are not in conflict: KR governs WHEN the")
    print("    window closes, while the day-15 endpoint is governed by how the")
    print("    symptom loop behaves inside a window that is open either way.")

    # =====================================================================
    head("TABLE 12 · Two mechanistic repairs, both rejected")
    # =====================================================================
    print("The model under-predicts maintained drug-placebo separation")
    print("(Table 5).  The two obvious repairs were implemented as switches and")
    print("tested rather than argued about.\n")
    print("REPAIR 1 (KR_BOOST): let PAM exposure ACCELERATE delta-pool")
    print("plasticity, so the drug speeds up the recovery it is bridging.")
    print("REPAIR 2 (W_PAM_BDNF): give the PAM a plasticity arm through BDNF,")
    print("so restoring tonic inhibition drives structural recovery.\n")
    print("%-26s %8s %8s %8s %8s"
          % ("setting", "base", "brx 60h", "brx d30", "zur d15"))
    print(line())
    for lbl, key, val in [("reference (both off)", None, None),
                          ("KR_BOOST = 1.0", 'KR_BOOST', 1.0),
                          ("KR_BOOST = 3.0", 'KR_BOOST', 3.0),
                          ("W_PAM_BDNF = 0.15", 'W_PAM_BDNF', 0.15),
                          ("W_PAM_BDNF = 0.35", 'W_PAM_BDNF', 0.35)]:
        pp = dict(p)
        if key:
            pp[key] = val
        rb = simulate(pp, V=V, t_end=enrol + 31 * 24, enrol=enrol, arm='brx60',
                      want=[])
        rz = simulate(pp, V=V, t_end=enrol + 16 * 24, enrol=enrol, arm='zur50',
                      want=[])
        b = hamd(gv(rb['at'](enrol), 'SYMP'), pp)
        bz = hamd(gv(rz['at'](enrol), 'SYMP'), pp)
        print("%-26s %8.1f %8.1f %8.1f %8.1f"
              % (lbl, b,
                 hamd(gv(rb['at'](enrol + 60.0), 'SYMP'), pp) - b,
                 hamd(gv(rb['at'](enrol + 30 * 24.0), 'SYMP'), pp) - b,
                 hamd(gv(rz['at'](enrol + 15 * 24.0), 'SYMP'), pp) - bz))
    print()
    print("BOTH REPAIRS MAKE THE FIT WORSE, AND FOR THE SAME REASON.")
    print("Brexanolone IS allopregnanolone (assumption A2 — a fact, not a")
    print("modelling choice), so the drug and the endogenous ligand enter the")
    print("model through the SAME potentiation term.  Any parameter that gives")
    print("the drug extra credit also gives PREGNANCY extra credit: it lowers")
    print("the baseline the patient starts from and shrinks the room a drug has")
    print("to work in.  The 'base' column above shows exactly that (28.8 ->")
    print("25.5 -> 23.9 as W_PAM_BDNF rises).  REPAIR 1 fails even more")
    print("directly: in a homeostatic set-point model a PAM cannot accelerate")
    print("recovery of the pool, because it LOWERS the set-point the pool is")
    print("chasing.  Speeding up the chase moves it toward a lower target.")
    print()
    print("So the missing durability mechanism cannot live on the")
    print("allopregnanolone-potentiation axis at all.  It has to be something")
    print("the drug does that pregnancy does not — a candidate the trials could")
    print("distinguish: a durable change in the sleep/behaviour loop bought by")
    print("60 hours of functioning (which this model would represent as a")
    print("lasting reduction in A_SYMP or WAKE, not as a receptor effect).")

    # =====================================================================
    head("SUMMARY OF NUMBERS QUOTED IN README.md")
    # =====================================================================
    rd_t = rd_setpoint(160.0, p)
    print("term delta-pool set-point          : %.3f (vs %.3f non-pregnant)"
          % (rd_t, rd_setpoint(6.0, p)))
    print("tonic inhibition at day-3 nadir    : %.3f of set-point (%.0f %% loss)"
          % (g_pp, 100 * (1 - g_pp)))
    print("trough wider than 10 %% deficit     : %.0f d" % (
        -math.log(1 - ((0.9 * p['G_TARGET'] / (1 + pot_tonic(4.8, p)) - rd_t)
                       / (rd_setpoint(4.8, p) - rd_t))) / p['KR'] / 24))
    print("60 -> 90 ug/kg/h potentiation gain : %+.1f %%"
          % (100 * (po90 / po60 - 1)))
    print("brexanolone 60 h, model vs placebo : %+.1f (obs -5.5)" % d_brx)
    print("zuranolone d15, model vs placebo   : %+.1f (obs -4.0)" % d_z50)
    print("brexanolone RID                    : see Table 9")
    print()
    print(line('='))
    print(" END OF REFERENCE OUTPUT")
    print(line('='))


if __name__ == '__main__':
    main()
