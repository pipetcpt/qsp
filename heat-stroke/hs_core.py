"""
hs_core.py -- Heat Stroke QSP model, executable core.

This is the reference implementation of the ODE system that is also written in
mrgsolve form in `hs_mrgsolve_model.R`.  Every number quoted in README.md and in
hs_references.md is produced by running this file (or hs_analysis.py, which
imports it).  R/mrgsolve is not installed in the build environment, so this
Python/scipy implementation is the executable source of truth and the R file is
the transcription; the two were written from the same equation sheet and the
symbol names match one-to-one.

Time unit: MINUTES throughout.
Temperature: degrees Celsius.
Heat: watts (J/s) internally, converted to J/min at the integrator boundary.

Structure
---------
 A. Psychrometrics + biophysical heat transfer   (Gagge/ASHRAE two-node, extended
                                                  to three nodes: core, muscle, skin)
 B. Thermoregulatory effectors                   (sweating, skin blood flow, shivering)
 C. Thermal dose                                 (Sapareto-Dewey CEM43)
 D. Cardiovascular + splanchnic steal
 E. Gut barrier -> endotoxin
 F. Cytokine cascade + HMGB1 commitment switch
 G. Endothelium + coagulation (DIC)
 H. Organ injury (liver, muscle, kidney, CNS)
 I. Drug PK/PD (paracetamol/NAPQI, ibuprofen, dantrolene, hydrocortisone,
                thrombomodulin alfa) and cooling as a conductance
"""

import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# state vector index map
# ----------------------------------------------------------------------------
SNAMES = [
    "TCR", "TMU", "TSK",          # 0-2   thermal nodes
    "CEM43", "CEM43E",            # 3-4   raw and HSP-protected thermal dose
    "HSP",                        # 5     HSP70 (fold of unstressed)
    "WDEF", "VP", "SWFAT",        # 6-8   water deficit (L), plasma volume (L), gland fatigue
    "GUT", "LPS",                 # 9-10  enterocyte injury (0-1), endotoxin (EU/mL)
    "TNF", "IL6", "IL1B", "IL10", # 11-14 pg/mL
    "HMGB1", "NETS",              # 15-16 ng/mL, arbitrary
    "SDC1", "TF", "THR",          # 17-19 syndecan-1 ng/mL, tissue factor, thrombin (nM)
    "FIB", "PLT", "PC", "DDIM",   # 20-23 mg/dL, 10^9/L, % activity, ug/mL FEU
    "ALT", "AST", "CK", "MB",     # 24-27 U/L, U/L, U/L, ug/L
    "SCR", "GFRF",                # 28-29 mg/dL, GFR fraction of normal
    "NSE", "CNSD",                # 30-31 ug/L, encephalopathy 0-1
    "LIVF",                       # 32    functional liver mass fraction
    "MUINJ",                      #       injured myocyte pool (CK/Mb source)
    "CAST",                       # 33    myoglobin cast burden
    "KPOT", "LAC",                # 34-35 mmol/L
    "PARA_A", "PARA_C", "PARA_P", # 36-38 paracetamol (mg)
    "NAPQI", "GSH",               # 39-40 NAPQI burden, hepatic GSH (fraction)
    "IBU_A", "IBU_C",             # 41-42 ibuprofen (mg)
    "DAN_C", "DAN_P",             # 43-44 dantrolene (mg)
    "HC_C",                       # 45    hydrocortisone (mg)
    "RTM_C",                      # 46    thrombomodulin alfa (U/kg-equivalent conc)
    "CORT",                       # 47    endogenous cortisol (ug/dL)
    "FLUID",                      # 48    cumulative IV crystalloid (L)
    "DOSEPRE",                    # 49    CEM43 accumulated before cooling started
]
IX = {n: i for i, n in enumerate(SNAMES)}
NSTATE = len(SNAMES)

# ----------------------------------------------------------------------------
# default parameters
# ----------------------------------------------------------------------------
P0 = dict(
    # --- anthropometry -------------------------------------------------------
    BW=70.0,            # kg
    HT=1.75,            # m
    CP_BODY=3470.0,     # J/(kg.K) whole-body specific heat
    FMU=0.40,           # muscle mass fraction
    FSK=0.10,           # skin/shell mass fraction
    # core fraction = 1 - FMU - FSK

    # --- metabolism ----------------------------------------------------------
    MREST=1.20,         # W/kg basal heat production
    Q10=2.30,           # van't Hoff coefficient for metabolic heat
    ETA_EX=0.20,        # mechanical efficiency of exercise
    FMU_EX=0.85,        # fraction of exercise heat deposited in muscle
    FMU_REST=0.22,      # fraction of resting heat deposited in muscle
    SHIVMAX=350.0,      # W maximal shivering heat production
    TSK_SHIV=31.0,      # skin temperature at which shivering engages
    TCR_SHIV=38.6,      # core temperature below which shivering is permitted

    # --- passive conductances (W/K) -----------------------------------------
    KCM0=70.0,          # core<->muscle tissue conduction
    KCS0=10.0,          # core<->skin tissue conduction
    KBLOOD=1.13,        # W/K per (L/h) of blood flow  (rho*c = 4081 J/(L.K))

    # --- skin blood flow -----------------------------------------------------
    SKBF0=0.30,         # L/min baseline
    SKBFMAX=7.50,       # L/min maximum
    TCVD0=36.80,        # core threshold for active vasodilation
    KVD=1.20,           # deg C for full vasodilation
    WSK_VD=0.15,        # weight of skin temperature in vasodilation drive
    FVD_AGE=1.00,       # multiplier: elderly ~0.60
    FVD_DRUG=1.00,      # multiplier: beta-blocker / anticholinergic
    TSK_VC=28.0,        # skin temperature at which cold vasoconstriction is complete
    # Cold-induced cutaneous vasoconstriction is INCOMPLETE while the central
    # hyperthermic drive is on.  This is the quantitative form of the standard
    # answer to the "ice water causes vasoconstriction and traps heat" objection
    # (Proulx 2003; Casa 2007): in a hyperthermic subject the central drive
    # overrides the local cold signal and immersion cooling is not impaired.
    FCOLD_FLOOR0=0.25,  # residual skin flow fraction with no central drive
    FCOLD_FLOOR1=0.45,  # extra residual fraction at full central drive
    KVC_COLD=0.45,      # tissue-conductance fraction retained when vasoconstricted

    # --- sweating ------------------------------------------------------------
    SWMAX=20.0,         # g/min maximal sweat output (=1.2 L/h, unacclimatised)
    SWGAIN=26.0,        # g/min per deg C of integrated drive
    TCSW0=37.00,        # core sweat onset threshold
    TSKSW0=33.0,
    WSK_SW=0.10,
    FSW_AGE=1.00,       # elderly ~0.65
    FSW_DRUG=1.00,      # anticholinergic ~0.35
    KSWFAT=0.0016,      # /min gland fatigue accrual at full wettedness
    KSWREC=0.0007,      # /min recovery
    SWFAT_MAX=0.30,     # maximal hidromeiosis (30% loss of output)
    KDEHY_SW=0.055,     # fractional sweat loss per % body-mass deficit

    # --- environment (may be overridden by a phase schedule) -----------------
    TA=25.0,            # air temperature deg C
    TR=None,            # mean radiant temperature; None -> TA
    RH=0.50,            # relative humidity fraction
    VAIR=1.0,           # m/s air velocity
    ICL=0.30,           # clo
    QSOL=0.0,           # W direct solar/radiant gain
    WMAX=0.85,          # maximum skin wettedness

    # --- cooling intervention -----------------------------------------------
    UA_COOL=0.0,        # W/K conductance of the cooling device (0 = none)
    T_COOL=2.0,         # deg C temperature of the cooling medium
    COOL_ON=1e9,        # min: time cooling starts
    COOL_STOP_TC=38.6,  # deg C: stop rule
    IVF_RATE=0.0,       # L/min cold-fluid infusion rate
    IVF_TEMP=4.0,       # deg C infusate temperature
    IVF_ON=1e9,
    IVF_DUR=0.0,        # min

    # --- thermal dose --------------------------------------------------------
    CEM_R_LO=0.25,      # Sapareto-Dewey R below 43 deg C
    CEM_R_HI=0.50,      # R at/above 43 deg C
    HSP50=2.5,          # HSP70 fold-induction giving 2x dose protection

    # --- HSP70 ---------------------------------------------------------------
    HSP_BASE=1.0,
    KHSP_IN=0.0025,     # /min
    KHSP_OUT=0.0009,    # /min  (t1/2 ~ 13 h)
    THSF=39.6,          # deg C midpoint of HSF1 activation
    NHSF=8.0,

    # --- water / plasma ------------------------------------------------------
    TBW_FRAC=0.60,      # total body water fraction
    VP0=3.0,            # L plasma volume
    KPV=1.30,           # fraction of water deficit borne by plasma (relative)
    INSENS=0.45,        # g/min insensible loss
    ORAL_RATE=0.0,      # L/min oral intake

    # --- cardiovascular ------------------------------------------------------
    CO0=5.0,            # L/min at rest
    KCO_EX=0.0110,      # L/min per W of exercise metabolic rate
    KCO_T=0.90,         # L/min per deg C of core hyperthermia
    COMAX=24.0,
    MAP0=90.0,          # mmHg
    KMAP_VP=1.60,       # sensitivity of MAP to fractional plasma volume loss
    KMAP_SKBF=0.045,    # mmHg lost per (L/min) of skin vasodilation
    QSPL0=1.50,         # L/min splanchnic flow at rest
    KSPL_EX=0.55,       # maximal fractional reduction by exercise
    KSPL_SK=0.45,       # maximal fractional reduction by cutaneous steal
    MEX_REF=900.0,      # W exercise metabolic rate giving maximal splanchnic steal
    QSPL_CRIT=0.90,     # L/min below which the mucosa is ischaemic

    # --- gut -----------------------------------------------------------------
    KG_ISCH=0.0060,     # /min gain of ischaemic enterocyte injury
    KG_HEAT=0.0130,     # per CEM43-minute of direct thermal injury
    KG_REP=0.00090,     # /min repair
    KPERM=0.35,         # GUT giving half-maximal permeability
    FLUX_LPS=0.08,      # EU/mL/min at full permeability
    CL_LPS=0.075,       # /min hepatic clearance at full liver function

    # --- cytokines (pg/mL, /min) --------------------------------------------
    KTNF=5.0, KDTNF=0.030, KM_LPS=0.55, WD_TNF=0.55,
    KIL6=11.0, KDIL6=0.0075, WT_IL6=0.0020, KI10=380.0,
    KIL1=0.66, KDIL1=0.016,
    KIL10=0.24, KDIL10=0.014, WCORT_IL10=0.020,
    KNET=0.020, KDNET=0.010,

    # --- HMGB1 commitment switch --------------------------------------------
    # Saddle-node bistability.  Necrotic cells release HMGB1; extracellular
    # HMGB1 drives RAGE/TLR4 signalling that produces more necrosis, and the
    # autocatalytic production SATURATES while clearance stays linear.  That is
    # the standard two-stable-state form:
    #     dH/dt = drive + A*H^n/(H^n+K^n) - c_eff*H,   c_eff = KH_OUT - KH_NEC*KNEC_H
    # The drive term is calibrated so that the switch latches at a thermal
    # dose of roughly 10-15 CEM43, which is where the epidemiology puts it:
    # exertional heat stroke cooled within ~30 min survives essentially
    # always, 60+ min of delay does not, and classic heat stroke arrives
    # having paid many times that.
    # With the values below the saddle-node sits at c_eff = 0.00498/min -- so an
    # agent that adds >0.00235/min of HMGB1 clearance abolishes the ON state
    # entirely rather than merely shifting it.  That is the structural claim the
    # model makes about thrombomodulin alfa (see KH_RTM).
    KH_NEC=14.5,        # ng/mL of HMGB1 per unit necrosis
    KH_AUTO=0.16,       # ng/mL/min saturating autocatalytic production
    KH_HALF=17.0,       # ng/mL midpoint of the autocatalytic term
    NH=3.0,
    KH_OUT=0.0035,      # /min clearance (t1/2 ~ 3.3 h)
    KH_RTM=0.0020,      # /min extra clearance per unit of rTM concentration
    KNEC_H=3.0e-5,      # necrosis driven per (ng/mL) HMGB1 per min
    KNEC_TH=0.080,      # necrosis driven per CEM43-minute

    # --- endothelium / coagulation ------------------------------------------
    SDC1_0=18.0, KSDC=0.85, KDSDC=0.0055,
    KTF=0.0055, KDTF=0.0060, KTF_LPS=0.75,
    KTHR=1.30, KDTHR=0.22, KTHR_PC=0.85,
    FIB0=300.0, KFIB_SYN=0.055, KFIB_APF=1.60, KFIB_CONS=0.055,
    PLT0=250.0, KPLT_SYN=0.00080, KPLT_CONS=0.019,
    PC0=100.0, KPC_SYN=0.00090, KPC_CONS=0.030,
    KDD=0.85, KDDD=0.0028,

    # --- organ injury --------------------------------------------------------
    ALT0=25.0, KALT=1350.0, KDALT=0.00040,       # t1/2 ~ 29 h
    AST0=25.0, KAST=1950.0, KDAST=0.00096,       # t1/2 ~ 12 h
    KAST_MU=0.0130,                              # muscle contribution to AST
    CK0=120.0, KCK=48000.0, KDCK=0.00032,        # t1/2 ~ 36 h
    MB0=30.0, KMB=8200.0, KDMB=0.0077,           # t1/2 ~ 90 min
    KMU_REL=0.0012,     # /min release from injured myocytes (t1/2 ~9.6 h)
    KRHAB_TH=0.145,                              # rhabdo per CEM43-min in muscle
    KRHAB_ISCH=0.0022,
    SCR0=0.90, KSCR_OUT=0.0125,   # KSCR_IN is derived: KSCR_OUT*SCR0 (see rhs)
    KGFR_VP=1.10, KGFR_CAST=0.55, KGFR_DIC=0.30, KGFR_REC=0.00035,
    KCAST=1.0e-7, KDCAST=0.00050,
    KLIV_TH=0.00185, KLIV_ISCH=0.00062, KLIV_NAP=0.00450, KLIV_REC=0.00030,
    KLIV_HMGB=0.00160,  # sustained hepatic injury while the switch is ON
    KGFR_HMGB=0.35,     # HMGB1 contribution to the GFR deficit
    KCNS_HMGB=0.00220,  # persistent encephalopathy while the switch is ON
    NSE0=8.0, KNSE=140.0, KDNSE=0.00048,
    KCNS_TH=0.0180, KCNS_INFL=0.00030, KCNS_REC=0.0035,
    KPOT0=4.0, KK_RHAB=0.000030, KK_REN=0.030,
    LAC0=1.0, KLAC=0.0130, KDLAC=0.0130,

    # --- endogenous cortisol -------------------------------------------------
    CORT0=12.0, KCORT_IN=0.030, KCORT_OUT=0.0075, KCORT_STIM=2.6,

    # --- paracetamol ---------------------------------------------------------
    PARA_KA=0.030, PARA_CL=21.0/60.0*1000.0,   # mL/min  (21 L/h)
    PARA_V1=32.0*1000.0, PARA_V2=18.0*1000.0, PARA_Q=8.0/60.0*1000.0,
    FNAPQI=0.055,        # fraction of clearance to NAPQI (raised by CYP2E1 induction)
    FNAPQI_HEAT=1.9,     # multiplier at full heat/fasting CYP2E1 induction
    KGSH_REC=0.00090, KGSH_USE=0.0130,
    KTSET_PARA=0.40,     # deg C set-point depression at 20 mg/L
    PARA_EC50=12.0,

    # --- ibuprofen -----------------------------------------------------------
    IBU_KA=0.030, IBU_CL=3.5/60.0*1000.0, IBU_V=10.0*1000.0,
    KTSET_IBU=0.60, IBU_EC50=15.0,
    KGFR_IBU=0.16,       # fractional GFR loss from renal prostaglandin blockade
    KPLTFN_IBU=0.55,     # platelet function loss (bleeding, not count)

    # --- dantrolene ----------------------------------------------------------
    DAN_CL=0.60/60.0*1000.0, DAN_V1=36.0*1000.0, DAN_V2=40.0*1000.0,
    DAN_Q=1.2/60.0*1000.0,
    DAN_EMAX=0.12,       # maximal fractional cut in muscle heat production
    DAN_EC50=2.5,        # mg/L

    # --- hydrocortisone ------------------------------------------------------
    HC_CL=18.0/60.0*1000.0, HC_V=35.0*1000.0,
    HC_EMAX=0.35, HC_EC50=0.15,   # fractional cytokine suppression

    # --- thrombomodulin alfa (ART-123) ---------------------------------------
    # Concentration is normalised so that one standard daily dose (380 U/kg)
    # is entered as an amount of 3.5, giving C_rtm = 1.0 unit/L at peak.
    RTM_CL=2.02, RTM_V=3.5*1000.0,     # CL in mL/min -> t1/2 ~ 20 h
    RTM_EAPC=1.60,       # fold increase in protein C activation
    RTM_EHMGB=1.00,      # scaling of the lectin-domain HMGB1 term

    # --- misc ----------------------------------------------------------------
    ACCLIM=0.0,          # 0 = unacclimatised, 1 = fully heat-acclimatised
)


def make_params(**kw):
    p = dict(P0)
    p.update(kw)
    if p["TR"] is None:
        p["TR"] = p["TA"]
    # heat acclimatisation is a coordinated set of changes, applied here so a
    # single knob moves all of them together (Periard 2021, Physiol Rev)
    a = p["ACCLIM"]
    if a:
        p["SWMAX"] = p["SWMAX"] * (1.0 + 0.85 * a)     # 1.2 -> 2.2 L/h
        p["TCSW0"] = p["TCSW0"] - 0.35 * a             # earlier onset
        p["TCVD0"] = p["TCVD0"] - 0.25 * a
        p["VP0"] = p["VP0"] * (1.0 + 0.115 * a)        # +11.5% plasma volume
        p["HSP_BASE"] = p["HSP_BASE"] * (1.0 + 1.30 * a)
        p["SWGAIN"] = p["SWGAIN"] * (1.0 + 0.30 * a)
    return p


# ----------------------------------------------------------------------------
# A. psychrometrics and biophysics
# ----------------------------------------------------------------------------
def psat(T):
    """Saturation vapour pressure in kPa (Buck 1981)."""
    return 0.61121 * np.exp((18.678 - T / 234.5) * (T / (257.14 + T)))


def dubois(bw, ht):
    return 0.202 * bw ** 0.425 * ht ** 0.725


def wet_bulb(Ta, RH):
    """Stull (2011) approximation to the psychrometric wet-bulb temperature."""
    rh = 100.0 * RH
    return (Ta * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
            + np.arctan(Ta + rh) - np.arctan(rh - 1.676331)
            + 0.00391838 * rh ** 1.5 * np.arctan(0.023101 * rh)
            - 4.686035)


def heat_exchange(TSK, p, env):
    """Dry (convective+radiative) loss in W and maximal evaporative capacity in W.

    Positive Q_dry means the body LOSES heat.  When ambient exceeds skin
    temperature Q_dry is negative: the environment is a heat source.
    """
    Ta, RH, v, Icl = env["TA"], env["RH"], env["VAIR"], env["ICL"]
    Tr = env.get("TR", Ta)
    AD = p["_AD"]
    hc = max(3.1, 8.3 * v ** 0.6)
    hr = 4.7
    h = hc + hr
    fcl = 1.0 + 0.31 * Icl
    Rcl = 0.155 * Icl
    To = (hc * Ta + hr * Tr) / h                      # operative temperature
    Q_dry = AD * (TSK - To) / (Rcl + 1.0 / (fcl * h))
    he = 16.5 * hc                                     # Lewis relation, W/(m2.kPa)
    Recl = 0.0276 * Icl
    Pa = RH * psat(Ta)
    E_max = AD * (psat(TSK) - Pa) / (Recl + 1.0 / (fcl * he))
    return Q_dry, max(E_max, 0.0), Pa


# ----------------------------------------------------------------------------
# environment / activity schedule
# ----------------------------------------------------------------------------
class Schedule:
    """Piecewise-constant environment and workload.

    phases: list of dicts, each with 'dur' (min) plus any of
    TA, TR, RH, VAIR, ICL, QSOL, MEX, UA_COOL, T_COOL, IVF_RATE, ORAL_RATE.
    """

    def __init__(self, phases, base):
        self.phases = phases
        self.base = base
        self.edges = np.cumsum([0.0] + [ph["dur"] for ph in phases])

    def at(self, t):
        env = dict(self.base)
        i = int(np.searchsorted(self.edges, t, side="right") - 1)
        i = min(max(i, 0), len(self.phases) - 1)
        env.update({k: v for k, v in self.phases[i].items() if k != "dur"})
        if env.get("TR") is None:
            env["TR"] = env["TA"]
        return env

    @property
    def tend(self):
        return float(self.edges[-1])


def hill(x, k, n):
    x = max(x, 0.0)
    return x ** n / (x ** n + k ** n)


# ----------------------------------------------------------------------------
# the right-hand side
# ----------------------------------------------------------------------------
def rhs(t, y, p, sched, aux=None):
    s = {n: y[i] for i, n in enumerate(SNAMES)}
    env = sched.at(t)
    BW = p["BW"]

    TCR, TMU, TSK = s["TCR"], s["TMU"], s["TSK"]

    # ---- drug concentrations (mg/L) ---------------------------------------
    C_para = s["PARA_C"] / p["PARA_V1"] * 1000.0
    C_ibu = s["IBU_C"] / p["IBU_V"] * 1000.0
    C_dan = s["DAN_C"] / p["DAN_V1"] * 1000.0
    C_hc = s["HC_C"] / p["HC_V"] * 1000.0
    C_rtm = s["RTM_C"] / p["RTM_V"] * 1000.0

    # ---- B. metabolic heat production -------------------------------------
    MEX = env.get("MEX", 0.0)                       # total metabolic rate of exercise, W
    # van't Hoff acceleration of resting metabolism.  This is the autocatalytic
    # term: hotter tissue makes more heat.  Above ~42 C it reverses as enzymes
    # denature, so the factor is rolled off rather than extrapolated.
    q10f = p["Q10"] ** ((min(TCR, 42.0) - 37.0) / 10.0)
    q10f *= 1.0 / (1.0 + max(0.0, TCR - 42.0) ** 2)
    M_rest = p["MREST"] * BW * q10f
    dan_eff = 1.0 - p["DAN_EMAX"] * C_dan / (C_dan + p["DAN_EC50"])
    H_ex = MEX * (1.0 - p["ETA_EX"]) * dan_eff
    # shivering: only when the shell is cold and the core is not hot
    shiv_drive = max(0.0, p["TSK_SHIV"] - TSK) / 4.0 * max(0.0, p["TCR_SHIV"] - TCR) / 1.5
    H_shiv = p["SHIVMAX"] * min(1.0, shiv_drive)
    H_prod = M_rest + H_ex + H_shiv

    f_mu = p["FMU_EX"] if MEX > 50.0 else p["FMU_REST"]
    H_mu = H_ex * f_mu + (M_rest + H_shiv) * p["FMU_REST"]
    H_cr = H_prod - H_mu

    # ---- effector: skin blood flow ----------------------------------------
    # antipyretic set-point depression enters HERE and in the sweat drive, i.e.
    # on the *thresholds* -- which is exactly why it cannot help once the
    # effectors are saturated.
    dTset = (p["KTSET_PARA"] * C_para / (C_para + p["PARA_EC50"])
             + p["KTSET_IBU"] * C_ibu / (C_ibu + p["IBU_EC50"]))
    pv_frac = s["VP"] / p["VP0"]
    vd_drive = ((TCR - (p["TCVD0"] - dTset)) + p["WSK_VD"] * (TSK - 33.0)) / p["KVD"]
    vd_drive = min(max(vd_drive, 0.0), 1.0)
    # cold-induced vasoconstriction, overridden in proportion to central drive
    fc_floor = p["FCOLD_FLOOR0"] + p["FCOLD_FLOOR1"] * vd_drive
    f_cold = fc_floor + (1.0 - fc_floor) * min(1.0, max(0.0, (TSK - p["TSK_VC"]) / 5.0))
    SKBF_raw = p["SKBF0"] + (p["SKBFMAX"] - p["SKBF0"]) * vd_drive * \
        p["FVD_AGE"] * p["FVD_DRUG"] * f_cold
    # one-pass (non-iterated) baroreflex: cutaneous vasodilation and hypovolaemia
    # both lower mean arterial pressure, and low pressure withdraws skin flow.
    map_now = (p["MAP0"] * (1.0 - p["KMAP_VP"] * max(0.0, 1.0 - pv_frac))
               - p["KMAP_SKBF"] * p["MAP0"] * (SKBF_raw - p["SKBF0"]))
    f_bp = min(1.0, max(0.25, map_now / p["MAP0"]))
    SKBF = p["SKBF0"] + (SKBF_raw - p["SKBF0"]) * f_bp
    # tissue conduction is itself reduced by vasoconstriction of the shell
    kcs_tissue = p["KCS0"] * (p["KVC_COLD"] + (1.0 - p["KVC_COLD"]) * vd_drive)
    KCS = kcs_tissue + p["KBLOOD"] * SKBF * 60.0
    Q_mu_flow = 0.75 + 0.0085 * MEX + 0.25 * max(0.0, TCR - 37.0)
    KCM = p["KCM0"] + p["KBLOOD"] * Q_mu_flow * 60.0

    # ---- effector: sweating -----------------------------------------------
    pct_def = 100.0 * s["WDEF"] / BW
    f_dehy = max(0.35, 1.0 - p["KDEHY_SW"] * pct_def)
    sw_drive = (TCR - (p["TCSW0"] - dTset)) + p["WSK_SW"] * (TSK - p["TSKSW0"])
    m_sw = p["SWGAIN"] * max(0.0, sw_drive)
    m_sw = min(m_sw, p["SWMAX"])
    m_sw *= p["FSW_AGE"] * p["FSW_DRUG"] * f_dehy * (1.0 - s["SWFAT"])

    # ---- A. environmental heat exchange ------------------------------------
    Q_dry, E_max_env, Pa = heat_exchange(TSK, p, env)
    Q_sol = env.get("QSOL", 0.0)
    LATENT = 2426.0                                  # J/g
    E_sweat_cap = m_sw * LATENT / 60.0               # W available if fully evaporated
    E_env_ceiling = E_max_env * p["WMAX"]
    E_req = max(0.0, H_prod - Q_dry + Q_sol)         # evaporation needed for balance
    # Actual evaporative cooling is the smaller of what the glands deliver and
    # what the air can absorb, and is never more than balance requires.  When
    # E_req exceeds capacity the subject still sweats maximally and the surplus
    # DRIPS: cooling is capped, water loss is not.  That asymmetry is why
    # uncompensable heat dehydrates faster than it cools.
    E_actual = min(E_sweat_cap, E_env_ceiling, E_req)
    Q_res = (0.0014 * (34.0 - env["TA"]) + 0.0173 * (5.87 - Pa)) * H_prod

    # ---- cooling device ----------------------------------------------------
    # Every modality is represented as one empirical lumped conductance UA
    # (W/K) to a medium at T_COOL.  UA is fitted per modality to published
    # whole-body core cooling rates (see hs_calibrate.py); it is not a
    # first-principles surface coefficient.  Making modality a single potency
    # number is the point: it puts "which cooler" and "how long until cooling"
    # on the same axis.
    UA = env.get("UA_COOL", 0.0)
    TCOOL = env.get("T_COOL", p["T_COOL"])
    # the stop rule is a property of the protocol, applied on the core
    if UA > 0.0 and TCR <= env.get("COOL_STOP_TC", p["COOL_STOP_TC"]):
        UA = 0.0
    Q_dev = UA * (TSK - TCOOL)
    if env.get("IMMERSE", 0.0) > 0.5 and UA > 0.0:
        # a submerged skin exchanges with the water, not with the air
        Q_dry = 0.0
        E_actual = 0.0
        Q_sol = 0.0

    # Cold intravenous fluid: an enthalpy sink proportional to infusion rate.
    # The infusate mixes with blood and is therefore distributed to the three
    # nodes in proportion to their perfusion, NOT dumped entirely on the core.
    # Loading it all on the core over-predicts the measured core fall by ~50%.
    ivr = env.get("IVF_RATE", 0.0)                   # L/min
    Q_ivf_tot = ivr * 1000.0 * 4.18 * (TCR - env.get("IVF_TEMP", p["IVF_TEMP"])) / 60.0
    CO = min(p["COMAX"], p["CO0"] + p["KCO_EX"] * MEX
             + p["KCO_T"] * max(0.0, TCR - 37.0))
    CO = max(CO, SKBF + Q_mu_flow + 0.5)
    f_sk_perf = SKBF / CO
    f_mu_perf = Q_mu_flow / CO
    f_cr_perf = 1.0 - f_sk_perf - f_mu_perf
    Q_ivf = Q_ivf_tot * f_cr_perf

    # ---- thermal node balances ---------------------------------------------
    m_cr = BW * (1.0 - p["FMU"] - p["FSK"])
    m_mu = BW * p["FMU"]
    m_sk = BW * p["FSK"]
    C_cr = m_cr * p["CP_BODY"]
    C_mu = m_mu * p["CP_BODY"]
    C_sk = m_sk * p["CP_BODY"]

    Q_cr_mu = KCM * (TCR - TMU)
    Q_cr_sk = KCS * (TCR - TSK)

    dTCR = 60.0 * (H_cr - Q_cr_mu - Q_cr_sk - Q_res - Q_ivf) / C_cr
    dTMU = 60.0 * (H_mu + Q_cr_mu - Q_ivf_tot * f_mu_perf) / C_mu
    dTSK = 60.0 * (Q_cr_sk - Q_dry - E_actual + Q_sol - Q_dev
                   - Q_ivf_tot * f_sk_perf) / C_sk

    # ---- C. thermal dose ----------------------------------------------------
    # Sapareto-Dewey cumulative equivalent minutes at 43 C.  The exponent is
    # clipped so a runaway core cannot overflow the integrator before the
    # terminal event fires.
    R = p["CEM_R_LO"] if TCR < 43.0 else p["CEM_R_HI"]
    dose_rate = min(R ** float(np.clip(43.0 - TCR, -6.0, 30.0)), 64.0)
    prot = 1.0 + (s["HSP"] - 1.0) / p["HSP50"]
    dose_eff = dose_rate / max(prot, 1.0)
    R_mu = p["CEM_R_LO"] if TMU < 43.0 else p["CEM_R_HI"]
    dose_mu = min(R_mu ** float(np.clip(43.0 - TMU, -6.0, 30.0)), 64.0) / max(prot, 1.0)

    dCEM43 = dose_rate
    dCEM43E = dose_eff
    dDOSEPRE = dose_rate if UA <= 0.0 and ivr <= 0.0 else 0.0

    # ---- HSP70 --------------------------------------------------------------
    hsf = hill(TCR - 36.0, p["THSF"] - 36.0, p["NHSF"])
    dHSP = p["KHSP_IN"] * hsf - p["KHSP_OUT"] * (s["HSP"] - p["HSP_BASE"])

    # ---- water balance ------------------------------------------------------
    # ORAL_MATCH replaces losses exactly (ad-libitum drinking to balance), which
    # isolates the THERMAL boundary from the HYDRATION boundary.  Without it, a
    # multi-hour exposure becomes uncompensable through dehydration at every
    # ambient temperature and the critical-environment question has no answer.
    if env.get("ORAL_MATCH", 0.0) > 0.5:
        intake = (m_sw + p["INSENS"]) / 1000.0 + ivr
    else:
        intake = env.get("ORAL_RATE", 0.0) + ivr      # L/min
    dWDEF = (m_sw + p["INSENS"]) / 1000.0 - intake
    tbw = p["TBW_FRAC"] * BW
    dVP = -p["KPV"] * (p["VP0"] / tbw) * ((m_sw + p["INSENS"]) / 1000.0 - intake)
    wet = 1.0 if E_env_ceiling <= 1e-9 else min(1.0, E_sweat_cap / max(E_env_ceiling, 1e-9))
    dSWFAT = p["KSWFAT"] * wet * (p["SWFAT_MAX"] - s["SWFAT"]) / p["SWFAT_MAX"] \
        - p["KSWREC"] * s["SWFAT"]
    dFLUID = ivr

    # ---- D. cardiovascular / splanchnic steal -------------------------------
    steal_ex = p["KSPL_EX"] * min(1.0, MEX / p["MEX_REF"])
    steal_sk = p["KSPL_SK"] * (SKBF - p["SKBF0"]) / (p["SKBFMAX"] - p["SKBF0"])
    QSPL = p["QSPL0"] * max(0.10, (1.0 - steal_ex) * (1.0 - steal_sk)) * f_bp
    ISCH = min(1.0, max(0.0, 1.0 - QSPL / p["QSPL_CRIT"]))

    # ---- E. gut barrier -----------------------------------------------------
    dGUT = (p["KG_ISCH"] * ISCH ** 2 + p["KG_HEAT"] * dose_eff) * (1.0 - s["GUT"]) \
        - p["KG_REP"] * s["GUT"]
    PERM = s["GUT"] / (s["GUT"] + p["KPERM"])
    livf = max(s["LIVF"], 0.05)
    dLPS = p["FLUX_LPS"] * PERM - p["CL_LPS"] * s["LPS"] * livf * min(1.0, QSPL / p["QSPL0"])

    # ---- F. cytokines and the commitment switch -----------------------------
    hc_sup = 1.0 - p["HC_EMAX"] * C_hc / (C_hc + p["HC_EC50"])
    cort_sup = 1.0 / (1.0 + 0.010 * max(0.0, s["CORT"] - p["CORT0"]))
    lps_sig = s["LPS"] / (s["LPS"] + p["KM_LPS"])
    hdamp = s["HMGB1"] / (s["HMGB1"] + 25.0)
    damp = hdamp + 0.5 * (s["HSP"] - 1.0) / 6.0
    dTNF = p["KTNF"] * (lps_sig + p["WD_TNF"] * damp) * hc_sup * cort_sup \
        - p["KDTNF"] * s["TNF"]
    il10_brake = 1.0 / (1.0 + s["IL10"] / p["KI10"])
    dIL6 = p["KIL6"] * (lps_sig + p["WT_IL6"] * s["TNF"] + 0.35 * damp) * il10_brake \
        * hc_sup - p["KDIL6"] * s["IL6"]
    dIL1B = p["KIL1"] * (lps_sig + 0.45 * damp) * hc_sup * cort_sup - p["KDIL1"] * s["IL1B"]
    dIL10 = p["KIL10"] * s["IL6"] / 100.0 + p["WCORT_IL10"] * s["CORT"] \
        - p["KDIL10"] * s["IL10"]
    dNETS = p["KNET"] * (s["IL1B"] / 100.0 + lps_sig) - p["KDNET"] * s["NETS"]

    # necrosis: direct thermal + HMGB1-driven microvascular injury
    NEC = p["KNEC_TH"] * dose_eff + p["KNEC_H"] * s["HMGB1"]
    rtm_h = p["KH_RTM"] * C_rtm * p["RTM_EHMGB"]
    # NOTE: the autocatalytic term is NOT multiplied by H.  Writing it as
    # A*H*hill(H) makes production super-linear at every H, which destroys the
    # OFF state and makes the switch latch unconditionally.  Saturating
    # production against linear clearance is what creates two stable states.
    dHMGB1 = (p["KH_NEC"] * NEC
              + p["KH_AUTO"] * hill(s["HMGB1"], p["KH_HALF"], p["NH"])
              - (p["KH_OUT"] + rtm_h) * s["HMGB1"])

    # ---- G. endothelium and coagulation -------------------------------------
    dSDC1 = p["KSDC"] * (dose_eff + 0.010 * s["TNF"] + 0.0035 * s["IL6"]) \
        - p["KDSDC"] * (s["SDC1"] - p["SDC1_0"])
    dTF = p["KTF"] * (p["KTF_LPS"] * lps_sig + 0.0075 * s["TNF"] + 0.0030 * s["IL1B"]) \
        - p["KDTF"] * s["TF"]
    pc_act = s["PC"] / p["PC0"] * (1.0 + (p["RTM_EAPC"] - 1.0) * C_rtm / (C_rtm + 0.5))
    dTHR = p["KTHR"] * s["TF"] - p["KDTHR"] * s["THR"] * (1.0 + p["KTHR_PC"] * pc_act)
    apf = s["IL6"] / (s["IL6"] + 120.0)               # acute-phase drive
    dFIB = (p["KFIB_SYN"] * (p["FIB0"] - s["FIB"]) * livf * (1.0 + p["KFIB_APF"] * apf)
            - p["KFIB_CONS"] * s["THR"] * s["FIB"] / 100.0)
    dPLT = (p["KPLT_SYN"] * (p["PLT0"] - s["PLT"])
            - p["KPLT_CONS"] * s["THR"] * s["PLT"] / 100.0)
    dPC = (p["KPC_SYN"] * (p["PC0"] - s["PC"]) * livf
           - p["KPC_CONS"] * s["THR"] * s["PC"] / 100.0
           * (1.0 + (p["RTM_EAPC"] - 1.0) * C_rtm / (C_rtm + 0.5)))
    dDDIM = p["KDD"] * s["THR"] * s["FIB"] / 300.0 - p["KDDD"] * s["DDIM"]

    # ---- H. organ injury ----------------------------------------------------
    nap_inj = p["KLIV_NAP"] * max(0.0, 1.0 - s["GSH"]) * s["NAPQI"] / 50.0
    liv_hit = (p["KLIV_TH"] * dose_eff + p["KLIV_ISCH"] * ISCH ** 2 * (1 + 2.0 * damp)
               + p["KLIV_HMGB"] * hdamp + nap_inj)
    dLIVF = -liv_hit * s["LIVF"] + p["KLIV_REC"] * (1.0 - s["LIVF"])
    hep_nec = liv_hit * s["LIVF"]
    dALT = p["KALT"] * hep_nec - p["KDALT"] * (s["ALT"] - p["ALT0"])
    rhab = p["KRHAB_TH"] * dose_mu + p["KRHAB_ISCH"] * ISCH ** 2 * min(1.0, MEX / 400.0)
    mu_rel = p["KMU_REL"] * s["MUINJ"]
    dMUINJ = rhab - mu_rel
    dAST = p["KAST"] * hep_nec + p["KAST_MU"] * p["KCK"] * mu_rel / 60.0 \
        - p["KDAST"] * (s["AST"] - p["AST0"])
    dCK = p["KCK"] * mu_rel - p["KDCK"] * (s["CK"] - p["CK0"])
    dMB = p["KMB"] * mu_rel - p["KDMB"] * (s["MB"] - p["MB0"])
    dCAST = p["KCAST"] * s["MB"] * (2.0 - min(1.0, pv_frac)) - p["KDCAST"] * s["CAST"]
    dic_burden = max(0.0, 1.0 - s["PLT"] / p["PLT0"])
    ibu_gfr = p["KGFR_IBU"] * C_ibu / (C_ibu + p["IBU_EC50"])
    gfr_target = max(0.05, 1.0 - p["KGFR_VP"] * max(0.0, 1.0 - pv_frac)
                     - p["KGFR_CAST"] * min(1.0, s["CAST"])
                     - p["KGFR_DIC"] * dic_burden - p["KGFR_HMGB"] * hdamp
                     - ibu_gfr)
    dGFRF = p["KGFR_REC"] * 160.0 * (gfr_target - s["GFRF"])
    # creatinine: constant generation, clearance proportional to GFR.  The
    # generation term is pinned to KSCR_OUT*SCR0 so that a healthy subject with
    # GFRF = 1 sits exactly at SCR0 instead of drifting upward.
    dSCR = p["KSCR_OUT"] * p["SCR0"] - p["KSCR_OUT"] * s["GFRF"] * s["SCR"]
    dNSE = p["KNSE"] * dose_eff * (1.0 + 1.5 * damp) - p["KDNSE"] * (s["NSE"] - p["NSE0"])
    dCNSD = (p["KCNS_TH"] * dose_eff + p["KCNS_INFL"] * (s["IL6"] / 100.0 + damp)
             + p["KCNS_HMGB"] * hdamp) * (1.0 - s["CNSD"]) - p["KCNS_REC"] * s["CNSD"]
    dKPOT = p["KK_RHAB"] * p["KCK"] * rhab - p["KK_REN"] * s["GFRF"] * (s["KPOT"] - p["KPOT0"])
    dLAC = p["KLAC"] * (ISCH + 0.4 * min(1.0, MEX / 600.0)) * 8.0 \
        - p["KDLAC"] * (s["LAC"] - p["LAC0"]) * livf

    # ---- endogenous cortisol -------------------------------------------------
    dCORT = p["KCORT_IN"] * p["KCORT_STIM"] * (dose_eff + s["IL6"] / 200.0) \
        - p["KCORT_OUT"] * (s["CORT"] - p["CORT0"])

    # ---- I. drug PK ----------------------------------------------------------
    dPARA_A = -p["PARA_KA"] * s["PARA_A"]
    cyp = 1.0 + (p["FNAPQI_HEAT"] - 1.0) * min(1.0, dose_eff / 0.25)
    cl_para = p["PARA_CL"] * livf
    dPARA_C = (p["PARA_KA"] * s["PARA_A"] - cl_para * C_para / 1000.0
               - p["PARA_Q"] * (C_para - s["PARA_P"] / p["PARA_V2"] * 1000.0) / 1000.0)
    dPARA_P = p["PARA_Q"] * (C_para - s["PARA_P"] / p["PARA_V2"] * 1000.0) / 1000.0
    dNAPQI = p["FNAPQI"] * cyp * cl_para * C_para / 1000.0 - 0.10 * s["NAPQI"]
    dGSH = p["KGSH_REC"] * (1.0 - s["GSH"]) - p["KGSH_USE"] * s["NAPQI"] / 50.0 * s["GSH"]

    dIBU_A = -p["IBU_KA"] * s["IBU_A"]
    dIBU_C = p["IBU_KA"] * s["IBU_A"] - p["IBU_CL"] * C_ibu / 1000.0 * livf

    dDAN_C = (-p["DAN_CL"] * C_dan / 1000.0 * livf
              - p["DAN_Q"] * (C_dan - s["DAN_P"] / p["DAN_V2"] * 1000.0) / 1000.0)
    dDAN_P = p["DAN_Q"] * (C_dan - s["DAN_P"] / p["DAN_V2"] * 1000.0) / 1000.0
    dHC_C = -p["HC_CL"] * C_hc / 1000.0
    dRTM_C = -p["RTM_CL"] * C_rtm / 1000.0

    d = np.zeros(NSTATE)
    d[IX["TCR"]] = dTCR
    d[IX["TMU"]] = dTMU
    d[IX["TSK"]] = dTSK
    d[IX["CEM43"]] = dCEM43
    d[IX["CEM43E"]] = dCEM43E
    d[IX["HSP"]] = dHSP
    d[IX["WDEF"]] = dWDEF
    d[IX["VP"]] = dVP
    d[IX["SWFAT"]] = dSWFAT
    d[IX["GUT"]] = dGUT
    d[IX["LPS"]] = dLPS
    d[IX["TNF"]] = dTNF
    d[IX["IL6"]] = dIL6
    d[IX["IL1B"]] = dIL1B
    d[IX["IL10"]] = dIL10
    d[IX["HMGB1"]] = dHMGB1
    d[IX["NETS"]] = dNETS
    d[IX["SDC1"]] = dSDC1
    d[IX["TF"]] = dTF
    d[IX["THR"]] = dTHR
    d[IX["FIB"]] = dFIB
    d[IX["PLT"]] = dPLT
    d[IX["PC"]] = dPC
    d[IX["DDIM"]] = dDDIM
    d[IX["ALT"]] = dALT
    d[IX["AST"]] = dAST
    d[IX["CK"]] = dCK
    d[IX["MB"]] = dMB
    d[IX["SCR"]] = dSCR
    d[IX["GFRF"]] = dGFRF
    d[IX["NSE"]] = dNSE
    d[IX["CNSD"]] = dCNSD
    d[IX["LIVF"]] = dLIVF
    d[IX["MUINJ"]] = dMUINJ
    d[IX["CAST"]] = dCAST
    d[IX["KPOT"]] = dKPOT
    d[IX["LAC"]] = dLAC
    d[IX["PARA_A"]] = dPARA_A
    d[IX["PARA_C"]] = dPARA_C
    d[IX["PARA_P"]] = dPARA_P
    d[IX["NAPQI"]] = dNAPQI
    d[IX["GSH"]] = dGSH
    d[IX["IBU_A"]] = dIBU_A
    d[IX["IBU_C"]] = dIBU_C
    d[IX["DAN_C"]] = dDAN_C
    d[IX["DAN_P"]] = dDAN_P
    d[IX["HC_C"]] = dHC_C
    d[IX["RTM_C"]] = dRTM_C
    d[IX["CORT"]] = dCORT
    d[IX["FLUID"]] = dFLUID
    d[IX["DOSEPRE"]] = dDOSEPRE

    if aux is not None:
        aux.update(dict(t=t, SKBF=SKBF, QSPL=QSPL, ISCH=ISCH, PERM=PERM,
                        m_sw=m_sw, E_max=E_env_ceiling, E_act=E_actual,
                        Q_dry=Q_dry, H_prod=H_prod, dose_rate=dose_rate,
                        MAP=map_now, KCS=KCS, E_req=E_req,
                        HSI=100.0 * E_req / max(E_env_ceiling, 1e-9),
                        Q_dev=Q_dev, Q_ivf=Q_ivf, wet=wet, C_para=C_para,
                        C_ibu=C_ibu, C_rtm=C_rtm, dTset=dTset))
    return d


# ----------------------------------------------------------------------------
# initial condition
# ----------------------------------------------------------------------------
def y0(p, TCR=37.0, TSK=33.0, TMU=37.2):
    y = np.zeros(NSTATE)
    y[IX["TCR"]] = TCR
    y[IX["TMU"]] = TMU
    y[IX["TSK"]] = TSK
    y[IX["HSP"]] = p["HSP_BASE"]
    y[IX["VP"]] = p["VP0"]
    y[IX["LIVF"]] = 1.0
    y[IX["GSH"]] = 1.0
    y[IX["FIB"]] = p["FIB0"]
    y[IX["PLT"]] = p["PLT0"]
    y[IX["PC"]] = p["PC0"]
    y[IX["SDC1"]] = p["SDC1_0"]
    y[IX["ALT"]] = p["ALT0"]
    y[IX["AST"]] = p["AST0"]
    y[IX["CK"]] = p["CK0"]
    y[IX["MB"]] = p["MB0"]
    y[IX["SCR"]] = p["SCR0"]
    y[IX["GFRF"]] = 1.0
    y[IX["NSE"]] = p["NSE0"]
    y[IX["KPOT"]] = p["KPOT0"]
    y[IX["LAC"]] = p["LAC0"]
    y[IX["CORT"]] = p["CORT0"]
    return y


# ----------------------------------------------------------------------------
# bolus dosing events
# ----------------------------------------------------------------------------
TCR_LETHAL = 44.0     # deg C: no survivor is reported above this; integration stops


def _lethal_event(t, y, p, sched):
    return y[IX["TCR"]] - TCR_LETHAL


_lethal_event.terminal = True
_lethal_event.direction = 1


def simulate(p, sched, doses=(), t_eval=None, y_init=None, max_step=2.0,
             stop_tc=None):
    """Integrate the system.

    doses: iterable of (time_min, state_name, amount) bolus additions.

    Integration terminates if the core reaches TCR_LETHAL.  Without this an
    uncompensable exposure has no fixed point and the solver runs the core to
    physically meaningless temperatures before overflowing -- the absence of an
    equilibrium is the disease, so it has to be handled explicitly rather than
    integrated through.
    """
    p = dict(p)
    p["_AD"] = dubois(p["BW"], p["HT"])
    y = y0(p) if y_init is None else np.array(y_init, dtype=float)
    tend = sched.tend
    if t_eval is None:
        t_eval = np.arange(0.0, tend + 0.5, 1.0)

    breaks = sorted(set([0.0, tend] + list(sched.edges)
                        + [d[0] for d in doses if 0.0 < d[0] < tend]))
    out_t, out_y = [], []
    t_death = None
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        if t1 <= t0:
            continue
        for (td, nm, amt) in doses:
            if abs(td - t0) < 1e-9:
                y[IX[nm]] += amt
        te = t_eval[(t_eval >= t0) & (t_eval <= t1)]
        if len(te) == 0 or te[0] > t0:
            te = np.concatenate(([t0], te))
        if te[-1] < t1:
            te = np.concatenate((te, [t1]))
        events = [_lethal_event]
        if stop_tc is not None:
            def _tc_event(t, yy, p_, s_, _tgt=stop_tc):
                return yy[IX["TCR"]] - _tgt
            _tc_event.terminal = True
            _tc_event.direction = 1
            events.append(_tc_event)
        sol = solve_ivp(rhs, (t0, t1), y, args=(p, sched), t_eval=te,
                        method="LSODA", rtol=1e-6, atol=1e-8, max_step=max_step,
                        events=events)
        if not sol.success:
            raise RuntimeError(f"integration failed on [{t0},{t1}]: {sol.message}")
        last = (i == len(breaks) - 2) or (sol.status == 1)
        out_t.append(sol.t if last else sol.t[:-1])
        out_y.append(sol.y if last else sol.y[:, :-1])
        if sol.status == 1:                      # a terminal event fired
            # events[0] is the lethal-temperature guard, events[1] (if present)
            # is the caller's stop_tc probe.  Only the former means death.
            if len(sol.t_events[0]):
                t_death = float(sol.t_events[0][0])
            break
        y = sol.y[:, -1].copy()

    T = np.concatenate(out_t)
    Y = np.concatenate(out_y, axis=1)
    p["_t_lethal"] = t_death
    return T, Y, p


def get(Y, name):
    return Y[IX[name], :]


def derived(T, Y, p, sched):
    """Recompute algebraic quantities along a solution."""
    keys = ["SKBF", "QSPL", "ISCH", "PERM", "m_sw", "E_max", "E_act", "Q_dry",
            "H_prod", "dose_rate", "MAP", "HSI", "Q_dev", "Q_ivf", "E_req",
            "C_para", "C_ibu", "C_rtm", "dTset"]
    out = {k: np.zeros(len(T)) for k in keys}
    aux = {}
    for i, t in enumerate(T):
        rhs(t, Y[:, i], p, sched, aux)
        for k in keys:
            out[k][i] = aux[k]
    return out


# ----------------------------------------------------------------------------
# clinical scores computed as OUTPUTS, never as states
# ----------------------------------------------------------------------------
def isth_dic_score(plt_, ddim, fib, ptratio):
    sc = 0
    sc += 2 if plt_ < 50 else (1 if plt_ < 100 else 0)
    sc += 3 if ddim > 5.0 else (2 if ddim > 1.0 else 0)
    sc += 1 if fib < 100 else 0
    sc += 2 if ptratio > 1.5 else (1 if ptratio > 1.25 else 0)
    return sc


def gcs_from_cnsd(c):
    """Map the encephalopathy state onto a GCS-like 3-15 scale."""
    return float(np.clip(15.0 - 12.0 * c, 3.0, 15.0))
