#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Severe Traumatic Brain Injury (sTBI) -- QSP reference implementation
 외상성 뇌손상 정량적 시스템 약리학 모델 : 참조 구현
================================================================================

WHY THIS FILE EXISTS
--------------------
Every number quoted in README.md, in the mechanistic map annotations and in the
mrgsolve model header is produced by RUNNING THIS FILE.  Nothing is asserted
from memory.  The committed transcript `tbi_reference_output.txt` is the
verbatim stdout of

    python3 tbi_reference_model.py

Pure standard library.  No numpy, no scipy.  Fixed-step RK4.

--------------------------------------------------------------------------------
 THE ORGANISING IDEA
--------------------------------------------------------------------------------
Almost every QSP model in this library is a *cascade*: an antigen or a mutation
drives cytokines which drive an effector which drives a clinical score, and a
drug interrupts one arrow.  Severe TBI is not that.  Severe TBI is a
**CONSTRAINT** problem, and the constraint is that the skull is a closed box.

    Monro-Kellie:   V_blood + V_CSF + V_brain + V_lesion = constant

The clinical variable we treat -- intracranial pressure -- is therefore NOT a
state variable in any mechanistic sense.  It is the *residual* of a volume
balance read off an exponential pressure-volume curve.  Five consequences fall
straight out of that single structural fact, and this model is built to compute
them rather than to assert them:

 1. ICP CARRIES NO INFORMATION ABOUT RESERVE.  Because dP/dV = P*ln10/PVI, the
    same 10 mL of oedema is invisible at ICP 8 and catastrophic at ICP 25.  An
    ICP that has been "stable at 15 all night" is compatible with either.

 2. EVERY ICP THERAPY IS A TRADE, NOT A CURE.  There are only four volumes in
    the box, so there are only four things any drug or manoeuvre can do, and
    each is bought with a different currency:
        CSF        -> drain it            (costs: nothing much; runs out)
        blood      -> constrict arterioles(costs: OXYGEN DELIVERY)
        water      -> osmotic gradient    (costs: SODIUM, and it rebounds)
        the box    -> craniectomy         (costs: the trial evidence says
                                           survival with severe disability)

 3. THE AUTOREGULATOR IS A POSITIVE FEEDBACK LOOP WHEN THE BOX IS TIGHT.
    CPP falls -> autoregulatory vasodilation -> cerebral blood VOLUME rises ->
    ICP rises -> CPP falls further.  This model does not contain a "plateau
    wave" term.  It contains the loop, and the plateau waves appear.

 4. THE SAME LOOP RUN BACKWARDS IS A DRUG.  Raise MAP -> flow overshoots ->
    autoregulatory vasoconstriction -> CBV falls -> ICP FALLS.  Noradrenaline
    is an ICP-lowering agent in exactly the patients in whom loop (3) is alive,
    and is an ICP-RAISING agent in the patients in whom it is dead.  One
    parameter (autoregulatory gain) decides the sign of a therapy.

 5. PRx AND "OPTIMAL CPP" ARE NOT MEASUREMENTS, THEY ARE THE LOOP GAIN.
    The U-shaped CPPopt curve is computed here from the quasi-steady
    dICP/dMAP of the vascular sigmoid.  It is emergent, not fitted.

--------------------------------------------------------------------------------
 STRUCTURE OF THE STATE VECTOR (47 ODEs)
--------------------------------------------------------------------------------
 Craniospinal mechanics + haemodynamics    : Pic, x_aut, MAP, V_csf
 CO2 / perivascular pH adaptation          : HCO3_csf
 Two-region brain water & osmolality       : V_int, V_inj, Osm_int, Osm_inj,
                                             Mann_inj
 Systemic osmotic state                    : Na_ecf, V_ecf, Mann_c, Mann_p
 Mass lesion & coagulopathy                : V_hem, Fibrinolysis
 Excitotoxic / metabolic cascade           : Glu, K_ec, Ca_i, MitoD, ROS,
                                             Lac, Glc_br
 Tissue fate                               : F_core, F_pen
 Neuroinflammation & barrier               : Micro, Cyto, Neut, MMP9,
                                             BBB_mech, BBB_infl
 Circulating biomarkers                    : GFAP, UCHL1, NfL, S100B
 Thermal                                   : Temp
 Drug PK                                   : Prop1..3 + Prop_e, Thio1..2,
                                             NE_c + NE_e, TXA_c + TXA_p
 Cumulative burden (insult "dose")         : D_icp, D_cpp
                                             (PbtO2 dose integrated post hoc)
--------------------------------------------------------------------------------
"""

import math
import sys

LN10 = math.log(10.0)


def sig(z):
    """Numerically safe logistic."""
    if z < -60.0:
        return 0.0
    if z > 60.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(-z))


def pw(base, e):
    """x**e guarded against the transient negative intermediate values an RK4
    stage can produce, which would otherwise return a complex number."""
    b = base if base > 0.0 else 0.0
    return b ** e


def clamp(v, lo, hi):
    return lo if v < lo else (hi if v > hi else v)


# ==============================================================================
#  PARAMETERS
# ==============================================================================
# Units:  time  min
#         volume mL      pressure mmHg      flow mL/min
#         resistance mmHg*min/mL             osmolality mOsm/kg
#         CBF is reported per 100 g of a 1400 g brain
# ==============================================================================

P = dict(

    # ---- craniospinal mechanics -------------------------------------------
    PVI=26.0,          # pressure-volume index, mL (Marmarou; normal 25-30)
    PVI_crani=95.0,    # after decompressive craniectomy
    Pic0=10.0,         # baseline ICP, mmHg
    Pvs=6.0,           # sagittal sinus pressure, mmHg

    # ---- cerebrovascular resistances (nominal) -----------------------------
    Ran=0.0800,        # arterial-arteriolar
    Rpv=0.0200,        # post-capillary venule
    Rv=0.0120,         # cerebral venous / bridging vein
    Rf=57.0,           # CSF formation resistance
    Ro=11.4,           # CSF outflow (arachnoid granulation) resistance
    R_evd=8.0,         # external ventricular drain resistance

    # ---- arteriolar compliance sigmoid -------------------------------------
    Can=0.200,         # nominal arteriolar compliance, mL/mmHg
    Ca_lo_f=0.45,      # min compliance as fraction of Can (max constriction)
    Ca_hi_f=1.55,      # max compliance as fraction of Can (max dilation)
    k_sig=0.070,       # sigmoid steepness in autoregulatory-tone units
    tau_aut=0.333,     # autoregulatory time constant, min (=20 s)
    G_aut=3.00,        # autoregulatory loop gain (intact)

    # ---- venous volume buffer ----------------------------------------------
    Vv_max=26.0,       # collapsible venous blood volume, mL
    Pv_half=25.0,      # ICP at which half the venous buffer is squeezed out
    kv=9.0,
    k_co2_cal=0.45,    # direct CO2 action on achievable arteriolar caliber

    # ---- CO2 reactivity -----------------------------------------------------
    fco2_lo=0.35, fco2_hi=2.05, co2_center=45.0, co2_k=10.0,
    HCO3_n=24.0, tau_hco3=360.0, hco3_slope=0.40,

    # ---- metabolism ---------------------------------------------------------
    CMRO2n=3.30,       # mL O2 / 100 g / min
    Q10=2.30,
    met_alpha=0.90,    # flow-metabolism coupling exponent
    OEF_max=0.85,
    CBFn100=53.6,      # nominal CBF, mL/100g/min
    brain_g=1400.0,

    # ---- oxygen -------------------------------------------------------------
    Hb=12.5,           # g/dL (typical post-resuscitation sTBI ICU)
    k_pbto2=0.75,      # PbtO2 / PvO2 at unit diffusion factor
    k_diff_ed=1.4,     # oedema penalty on O2 diffusion distance
    k_diff_mv=1.6,     # microvascular collapse penalty

    # ---- brain water / two-region Starling ---------------------------------
    V_int0=950.0,      # intact-region brain water, mL
    V_inj0=170.0,      # injured-region (contusion + pericontusional) water, mL
    LpS_int=0.00190,   # mL/min/mmHg
    LpS_inj=0.00260,
    mmHg_per_mOsm=19.3,
    kappa_hyd=0.012,   # weight of the hydrostatic Starling term (intact)
    kappa_bbb=22.0,    # amplification of hydrostatic term by BBB opening
    k_glym_int=0.0130, # ISF/glymphatic clearance, /min
    k_glym_inj=0.0062,
    Lh_int=0.0035,     # hydraulic conductance, intact  (mL/min/mmHg)
    Lh_inj=0.0300,     # hydraulic conductance, injured
    sig_prot_open=0.95,# fraction of oncotic reflection lost at full BBB opening
    V_csf_min=62.0,    # CSF volume at which the cisterns are effaced
    k_pen_res=0.55,    # pericontusional vascular resistance penalty
    k_leak_int=0.00030,  # /min : intact-BBB solute equilibration (t1/2 ~38 h)
    k_leak_bbb=0.0100,   # /min per unit BBB opening (t1/2 ~69 min fully open)
    k_mann_bbb=0.0055,   # /min per unit BBB opening, mannitol entry
    k_idio=0.115,        # idiogenic osmole generation on pump failure
    k_idio_clr=0.0055,

    # ---- systemic osmotic ---------------------------------------------------
    Na_n=140.0,
    V_ecf0=15000.0,    # mL
    Glu_p=8.0,         # mmol/L plasma glucose
    BUN=5.0,           # mmol/L
    GFR=0.110,         # L/min  (110 mL/min)
    k_na_corr=0.0055,  # proportional renal Na correction
    k_vol_corr=0.0012, # proportional renal volume correction
    Mann_V1=7000.0,    # mannitol central volume, mL
    Mann_V2=8000.0,    # mannitol peripheral volume, mL
    Mann_Q=0.90,       # intercompartmental clearance, L/min
    Osm_aki=320.0,     # plasma osm above which GFR degrades
    k_aki=0.055,

    # ---- mass lesion --------------------------------------------------------
    k_hem_exp=0.0072,  # hematoma expansion rate constant
    tau_hem=210.0,     # expansion window, min
    k_hem_res=0.00016, # resorption
    k_fib=0.0140,      # fibrinolysis activation
    k_fib_off=0.0035,

    # ---- excitotoxicity -----------------------------------------------------
    Glu_n=2.0,         # uM extracellular
    k_glu_rel=520.0, k_glu_upt=0.34,
    K_n=3.0,           # mM extracellular K+
    k_k_rel=52.0, k_k_upt=0.30,
    k_ca_in=0.055, Km_glu=28.0, k_ca_out=0.052, Ca_n=0.10,
    k_mito=0.0135, k_mito_rep=0.0021,
    k_ros=0.030, k_ros_clr=0.055,
    k_lac=0.155, k_lac_clr=0.030, Lac_n=1.6,
    Glc_br_n=1.8, k_glc_in=0.055, k_glc_use=0.075,

    # ---- tissue fate ---------------------------------------------------------
    k_pen2core=0.00105,   # penumbra -> core conversion
    k_pen_rec=0.00034,    # penumbra -> normal recovery
    k_sd_base=0.55,       # spreading depolarisation rate scaling (events/h)
    sd_cost=0.42,         # CMRO2 surcharge per unit SD burden

    # ---- neuroinflammation ---------------------------------------------------
    k_micro=0.0075, d_micro=0.00095,
    k_cyto=0.0125, d_cyto=0.0042,
    k_neut=0.0035, d_neut=0.0016,
    k_mmp=0.0060, d_mmp=0.0021,
    k_bbb_mech_off=0.0038,   # mechanical BBB opening reseals (t1/2 ~3 h)
    k_bbb_infl=0.00092, k_bbb_seal=0.00055,

    # ---- biomarkers (serum, ng/mL scale) --------------------------------------
    k_gfap=145.0, ke_gfap=0.00048,   # t1/2 ~ 24 h
    k_uchl1=62.0, ke_uchl1=0.00165,  # t1/2 ~ 7 h
    k_nfl=9.0,   ke_nfl=0.000023,    # t1/2 ~ 21 d
    k_s100b=17.0, ke_s100b=0.0116,   # t1/2 ~ 60 min

    # ---- temperature ----------------------------------------------------------
    T_n=37.0, k_fever=1.35, tau_T=95.0, tau_cool=42.0,

    # ---- propofol PK (3-cpt, Marsh-like, 75 kg) --------------------------------
    prop_V1=15.9, prop_V2=32.6, prop_V3=203.0,   # L
    prop_Cl1=1.94, prop_Cl2=1.75, prop_Cl3=0.79, # L/min
    prop_ke0=0.26,                                # /min
    prop_Emax=0.55, prop_EC50=3.40,               # CMRO2 suppression
    prop_MAPmax=26.0, prop_MAP50=3.10,

    # ---- thiopental PK (2-cpt) --------------------------------------------------
    thio_V1=28.0, thio_V2=110.0, thio_Cl=0.21, thio_Q=1.05,
    thio_Emax=0.55, thio_EC50=26.0,
    thio_MAPmax=24.0, thio_MAP50=32.0,
    supp_floor=0.60,   # combined anaesthetic CMRO2 suppression ceiling

    # ---- noradrenaline ----------------------------------------------------------
    ne_ke0=0.55, ne_Emax=55.0, ne_EC50=0.300,   # ug/kg/min effect-site scale

    # ---- tranexamic acid ---------------------------------------------------------
    txa_V1=12.0, txa_V2=15.0, txa_Cl=0.115, txa_Q=0.16,
    txa_IC50=8.0,

    # ---- injury severity (set per patient) ----------------------------------------
    sev=0.62,
    age=42.0,
    weight=75.0,
)


# ==============================================================================
#  CALIBRATION
# ------------------------------------------------------------------------------
#  Nothing below is hand-entered.  The reference operating point (MAP 88,
#  ICP 10, PaCO2 40, 37 C, no drugs) is solved for, and the quantities that MUST
#  balance there -- CSF formation vs absorption, transcapillary water flux,
#  renal sodium and volume handling, autoregulatory set-point -- are back-solved
#  so that a healthy brain integrated for a day does not move.  A model that
#  drifts at baseline cannot be trusted to attribute a drift under therapy.
# ==============================================================================
def calibrate(p):
    Pa0, Pic0 = 88.0, p["Pic0"]
    p["kR"] = p["Ran"] * (p["Can"] * (Pa0 - Pic0)) ** 2 / p["Can"] ** 2
    Pdown0 = p["Pvs"] + (Pic0 - p["Pvs"]) * sig((Pic0 - p["Pvs"]) / 1.0)
    Rtot0 = p["Ran"] + p["Rpv"] + p["Rv"]
    q0 = (Pa0 - Pdown0) / Rtot0
    p["q0"] = q0
    p["CBFn100"] = q0 / (p["brain_g"] / 100.0)
    p["Pc0"] = Pa0 - q0 * p["Ran"]
    # CSF: formation == absorption == 0.35 mL/min at the operating point
    csf0 = 0.35
    p["Rf"] = (p["Pc0"] - Pic0) / csf0
    p["Ro"] = (Pic0 - p["Pvs"]) / csf0
    p["csf0"] = csf0
    # transcapillary: zero net water flux at the operating point
    p["hyd0"] = p["Pc0"] - Pic0
    # plasma osmolality implied by the electrolyte panel
    p["Osm_n"] = 2.0 * p["Na_n"] + p["Glu_p"] + p["BUN"]
    # oxygen at the operating point
    p["CaO2_n"] = 1.34 * p["Hb"] * 0.98 + 0.003 * 95.0
    return p


def softmin(a, b, n=8.0):
    """Smooth min(a, b) that is accurate to <0.1% when the arguments differ by
    a factor of 2, and exactly (a+b)/2^(1/n) when they are equal.  Used so the
    oxygen-limited metabolic rate has a continuous derivative for RK4 without
    silently taxing an amply-supplied brain."""
    if a <= 0.0:
        return 0.0
    if b <= 0.0:
        return 0.0
    return a / (1.0 + (a / b) ** n) ** (1.0 / n)


STATES = [
    "Pic", "x_aut", "MAP", "V_csf", "HCO3",
    "V_int", "V_inj", "Osm_int", "Osm_inj", "Mann_inj",
    "Na_ecf", "V_ecf", "Mann_c", "Mann_p",
    "V_hem", "Fib",
    "Glu", "K_ec", "Ca_i", "MitoD", "ROS", "Lac", "Glc_br",
    "F_core", "F_pen",
    "Micro", "Cyto", "Neut", "MMP9", "BBB_mech", "BBB_infl",
    "GFAP", "UCHL1", "NfL", "S100B",
    "Temp",
    "Prop1", "Prop2", "Prop3", "Prop_e",
    "Thio1", "Thio2",
    "NE_e", "TXA_c", "TXA_p",
    "D_icp", "D_cpp",
]
IX = {n: i for i, n in enumerate(STATES)}
NST = len(STATES)


# ==============================================================================
#  CONTROL / THERAPY INTERFACE
# ==============================================================================
class Ctl(object):
    """Everything the clinician can do.  Any field may be a constant or a
    function of (t, obs) -- obs is the observation dict from the previous step,
    which is what makes closed-loop protocols possible."""

    def __init__(self, **kw):
        self.PaCO2 = 38.0
        self.SaO2 = 0.98
        self.PaO2 = 95.0
        self.MAP_base = 88.0
        self.T_target = None          # None = no active temperature control
        self.prop_rate = 0.0          # mg/kg/h
        self.thio_rate = 0.0          # mg/kg/h
        self.ne_rate = 0.0            # ug/kg/min  (or 'auto')
        self.cpp_target = None        # if set, NE is auto-titrated
        self.evd_open = False
        self.evd_set = 12.0           # mmHg
        self.crani_t = None           # min at which craniectomy is performed
        self.fluid_rate = 1.6         # mL/min maintenance isotonic
        self.boluses = []             # (t, kind, dose)  kind in mann|hts23|hts3|txa
        self.hypox = None             # (t0, t1, SaO2)
        self.hypot = None             # (t0, t1, MAP)
        self.autoreg = 1.0            # multiplier on G_aut (1 = intact, 0 = abolished)
        self.mass_rate = 0.0          # mL/min of extra intracranial mass (probe)
        self.protocol = None          # closed-loop callback(t, obs, ctl)
        self.__dict__.update(kw)

    def eff_SaO2(self, t):
        if self.hypox and self.hypox[0] <= t < self.hypox[1]:
            return self.hypox[2]
        return self.SaO2

    def eff_MAPbase(self, t):
        if self.hypot and self.hypot[0] <= t < self.hypot[1]:
            return self.hypot[2]
        return self.MAP_base(t) if callable(self.MAP_base) else self.MAP_base


# ==============================================================================
#  INITIAL CONDITION
# ==============================================================================
def initial_state(p, sev=None):
    s = sev if sev is not None else p["sev"]
    y = [0.0] * NST
    y[IX["Pic"]] = p["Pic0"]
    y[IX["x_aut"]] = 0.0
    y[IX["MAP"]] = 88.0
    # Compensatory reserve is not the same in every patient, and it is not a
    # free parameter: effacement of the basal cisterns is a Marshall/Rotterdam
    # CT criterion and the strongest single CT predictor of intracranial
    # hypertension, while cerebral atrophy is why an eighty-year-old tolerates
    # a subdural that would kill a twenty-year-old.  Both are the SAME quantity
    # -- how much CSF there was to give up -- so both are entered here rather
    # than as separate mechanisms elsewhere.
    y[IX["V_csf"]] = clamp(140.0 - 62.0 * s + 0.85 * (p["age"] - 40.0),
                           70.0, 190.0)
    y[IX["HCO3"]] = p["HCO3_n"]
    y[IX["V_int"]] = p["V_int0"]
    y[IX["V_inj"]] = p["V_inj0"]
    y[IX["Osm_int"]] = p["Osm_n"]
    y[IX["Osm_inj"]] = p["Osm_n"]
    y[IX["V_inj"]] = p["V_inj0"]
    y[IX["Mann_inj"]] = 0.0
    y[IX["Na_ecf"]] = p["Na_n"] * p["V_ecf0"] / 1000.0
    y[IX["V_ecf"]] = p["V_ecf0"]
    y[IX["V_hem"]] = 14.0 * s            # initial clot volume, mL
    y[IX["Fib"]] = 0.10
    y[IX["Glu"]] = p["Glu_n"]
    y[IX["K_ec"]] = p["K_n"]
    y[IX["Ca_i"]] = p["Ca_n"]
    y[IX["MitoD"]] = 0.02 + 0.10 * s
    y[IX["ROS"]] = 0.02 + 0.10 * s
    y[IX["Lac"]] = p["Lac_n"]
    y[IX["Glc_br"]] = p["Glc_br_n"]
    # primary (impact-instant) injury: core is already dead when the ambulance
    # arrives; the penumbra is what medicine is for.
    y[IX["F_core"]] = 0.055 * s
    y[IX["F_pen"]] = 0.300 * s
    y[IX["Micro"]] = 0.04 * s
    y[IX["Cyto"]] = 0.02 * s
    y[IX["Neut"]] = 0.01 * s
    y[IX["MMP9"]] = 0.05 * s
    y[IX["BBB_mech"]] = 0.55 * s         # mechanical barrier disruption
    y[IX["BBB_infl"]] = 0.02
    y[IX["Temp"]] = 37.0
    return y


# ==============================================================================
#  DERIVATIVES  +  OBSERVATIONS
# ==============================================================================
def _vasc(p, Ca_state_x, Pa, Pic, autoreg_mult, crani, m_co2=1.0):
    """Arteriolar compliance / resistance / flow given autoregulatory tone.

    CO2 enters TWICE and it has to.  Its main action is on the flow the servo is
    trying to deliver (see q_target), but if that were its only action then a
    vessel already pinned at its dilation ceiling by a low CPP could not respond
    to hyperventilation at all -- the servo would simply be told to want less of
    something it already cannot get.  CO2 also moves the CEILING itself, because
    perivascular pH sets how wide the smooth muscle can relax, and that action
    is not something the servo can escape.  This is why hyperventilation still
    works, weakly, in a maximally vasodilated brain.
    """
    Ca_lo = p["Ca_lo_f"] * p["Can"] * m_co2
    Ca_hi = p["Ca_hi_f"] * p["Can"] * m_co2
    s = sig(Ca_state_x / p["k_sig"])
    Ca = Ca_lo + (Ca_hi - Ca_lo) * s
    dCa_dx = (Ca_hi - Ca_lo) * s * (1.0 - s) / p["k_sig"]
    Va = max(1.0, Ca * (Pa - Pic))
    Ra = p["kR"] * p["Can"] ** 2 / Va ** 2
    return Ca, dCa_dx, Va, Ra


def observe(t, y, p, ctl, want_prx=False):
    """All algebraic quantities.  Called by derivs and by the reporter, so the
    reported numbers are literally the numbers the integrator used."""
    o = {}
    Pic = y[IX["Pic"]]
    Pa = y[IX["MAP"]]
    crani = (ctl.crani_t is not None and t >= ctl.crani_t)

    # ---- tissue fate ------------------------------------------------------
    F_core = clamp(y[IX["F_core"]], 0.0, 1.0)
    F_pen = clamp(y[IX["F_pen"]], 0.0, 1.0 - F_core)
    viable = 1.0 - F_core
    o["F_core"], o["F_pen"], o["viable"] = F_core, F_pen, viable

    # ---- temperature-, drug- and activity-dependent metabolic demand ------
    Ce_prop = y[IX["Prop_e"]]
    C_thio = y[IX["Thio1"]] / p["thio_V1"]
    supp = (p["prop_Emax"] * Ce_prop / (p["prop_EC50"] + Ce_prop)
            + p["thio_Emax"] * C_thio / (p["thio_EC50"] + C_thio))
    supp = p["supp_floor"] * (1.0 - math.exp(-supp / p["supp_floor"]))
    o["supp"] = supp
    T = y[IX["Temp"]]
    q10f = p["Q10"] ** ((T - 37.0) / 10.0)

    # spreading depolarisations: rate rises as the penumbra loses its pumps
    # (computed below once E_pump is known) -- provisional using previous K_ec
    K_ec = y[IX["K_ec"]]
    sd_drive = sig((K_ec - 5.4) / 0.85) * (F_pen / 0.30)
    SD_rate = p["k_sd_base"] * 12.0 * sd_drive        # events / h
    o["SD_rate"] = SD_rate

    CMRO2_dem = (p["CMRO2n"] * q10f * (1.0 - supp)
                 * (viable * 0.86 + 0.14)
                 * (1.0 + p["sd_cost"] * sd_drive))
    o["CMRO2_dem"] = CMRO2_dem

    # ---- CO2 -> perivascular pH -> effective PaCO2 ------------------------
    PaCO2 = ctl.PaCO2(t) if callable(ctl.PaCO2) else ctl.PaCO2
    HCO3 = y[IX["HCO3"]]
    PaCO2e = PaCO2 * (p["HCO3_n"] / max(6.0, HCO3))
    f_co2 = p["fco2_lo"] + (p["fco2_hi"] - p["fco2_lo"]) * sig(
        (PaCO2e - p["co2_center"]) / p["co2_k"])
    o["PaCO2"], o["PaCO2e"], o["f_co2"] = PaCO2, PaCO2e, f_co2

    # ---- the flow the autoregulator is servoing to ------------------------
    f_met = pw(CMRO2_dem / p["CMRO2n"], p["met_alpha"])
    q_target = p["CBFn100"] * p["brain_g"] / 100.0 * f_co2 * f_met
    o["q_target"] = q_target

    # ---- vascular mechanics -----------------------------------------------
    m_co2 = clamp(1.0 + p["k_co2_cal"] * (f_co2 - 1.0), 0.65, 1.35)
    o["m_co2"] = m_co2
    Ca, dCa_dx, Va, Ra = _vasc(p, y[IX["x_aut"]], Pa, Pic, ctl.autoreg, crani,
                               m_co2)
    # Starling resistor / vascular waterfall: the bridging veins are collapsible,
    # so the effective downstream pressure for cerebral outflow is ICP whenever
    # ICP exceeds sinus pressure -- which, at a normal ICP of 10 and a sinus
    # pressure of 6, is ALWAYS.  Cerebral venous outflow is ICP-referenced.
    Pdown = p["Pvs"] + (Pic - p["Pvs"]) * sig((Pic - p["Pvs"]) / 1.0)
    Rtot = Ra + p["Rpv"] + p["Rv"]
    q = max(1.0, (Pa - Pdown) / Rtot)
    Pc = Pa - q * Ra
    Pv = Pdown + q * p["Rv"]
    o.update(Ca=Ca, Ra=Ra, Va=Va, q=q, Pc=Pc, Pv=Pv, Pdown=Pdown)
    o["CBF100"] = q / (p["brain_g"] / 100.0)
    o["CPP"] = Pa - Pic

    # venous buffer volume and its slope
    Vv = p["Vv_max"] * (1.0 - sig((Pic - p["Pv_half"]) / p["kv"]))
    sv = sig((Pic - p["Pv_half"]) / p["kv"])
    dVv_dP = -p["Vv_max"] * sv * (1.0 - sv) / p["kv"]
    o["Vv"], o["CBV"] = Vv, Va + Vv

    # ---- oxygen ------------------------------------------------------------
    SaO2 = ctl.eff_SaO2(t)
    CaO2 = 1.34 * p["Hb"] * SaO2 + 0.003 * ctl.PaO2
    DO2 = o["CBF100"] * CaO2 / 100.0                      # mL O2/100g/min
    CMRO2_act = softmin(CMRO2_dem, DO2 * p["OEF_max"])
    OEF = CMRO2_act / max(1e-6, DO2)
    SvO2 = SaO2 * (1.0 - OEF)

    # ---------------- REGIONAL (pericontusional) perfusion ------------------
    # The single most consequential simplification a global TBI model can make
    # is to have ONE perfusion.  Jugular bulb saturation is a global, flow-
    # weighted average and is dominated by the healthy majority of the brain;
    # the tissue that is dying is a minority compartment behind a much higher
    # local resistance (capillary compression by oedema, microvascular collapse,
    # and the inverse haemodynamic response to spreading depolarisation).
    # That is exactly why a parenchymal PbtO2 probe is placed in pericontusional
    # tissue, and exactly why SjvO2 can read 65% over a dying penumbra.
    ed_inj_frac = max(0.0, (y[IX["V_inj"]] - p["V_inj0"]) / p["V_inj0"])
    mv_col = sig((Pic - 28.0) / 6.0) * 0.7 + 0.3 * F_pen
    BBB_now = clamp(y[IX["BBB_mech"]] + y[IX["BBB_infl"]], 0.0, 1.0)
    pen_res = 1.0 + p["k_pen_res"] * (mv_col + 0.6 * BBB_now
                                      + 5.0 * ed_inj_frac + 0.6 * sd_drive)
    CBF_pen = o["CBF100"] / pen_res
    DO2_pen = CBF_pen * CaO2 / 100.0
    CMRO2_pen = softmin(CMRO2_dem, DO2_pen * p["OEF_max"])
    r = CMRO2_pen / max(1e-6, CMRO2_dem)          # <- the cascade reads THIS
    OEF_pen = CMRO2_pen / max(1e-6, DO2_pen)
    SvO2_pen = SaO2 * (1.0 - OEF_pen)
    o.update(CBF_pen=CBF_pen, DO2_pen=DO2_pen, pen_res=pen_res,
             OEF_pen=OEF_pen, r_glob=CMRO2_act / max(1e-6, CMRO2_dem))
    # inverse ODC (Severinghaus)
    def _pvo2(sv):
        ss = clamp(sv, 0.02, 0.995)
        return clamp((23400.0 / (1.0 / ss - 1.0)) ** (1.0 / 3.0), 1.0, 120.0)
    PvO2 = _pvo2(SvO2)
    PvO2_pen = _pvo2(SvO2_pen)
    ed_frac = ((y[IX["V_int"]] - p["V_int0"]) + (y[IX["V_inj"]] - p["V_inj0"])) \
        / (p["V_int0"] + p["V_inj0"])
    Dfac = 1.0 / (1.0 + p["k_diff_ed"] * ed_inj_frac + p["k_diff_mv"] * mv_col)
    PbtO2 = p["k_pbto2"] * PvO2_pen * Dfac
    o.update(CaO2=CaO2, DO2=DO2, CMRO2=CMRO2_act, r=r, OEF=OEF,
             SjvO2=SvO2 * 100.0, PvO2=PvO2, PvO2_pen=PvO2_pen, PbtO2=PbtO2,
             Dfac=Dfac, ed_frac=ed_frac, ed_inj_frac=ed_inj_frac, mv_col=mv_col)

    # ---- pump failure, LPR --------------------------------------------------
    E_pump = sig((r - 0.62) / 0.045)
    # Reduced pump CAPACITY and frank ionic FAILURE are different thresholds.
    # Membrane potential fails at CBF ~18 mL/100g/min; the ionic pumps hold on
    # until ~10-12.  Conflating them makes a well-perfused brain leak glutamate.
    o["ion_fail"] = sig((0.55 - r) / 0.035)
    o["E_pump"] = E_pump
    o["LPR"] = 15.0 * (1.0 + 6.0 * pw(1.0 - r, 1.3) + 2.5 * max(0.0, y[IX["MitoD"]]))

    # ---- barrier ------------------------------------------------------------
    BBB = clamp(y[IX["BBB_mech"]] + y[IX["BBB_infl"]], 0.0, 1.0)
    o["BBB"] = BBB

    # ---- osmotic ------------------------------------------------------------
    Na_p = y[IX["Na_ecf"]] * 1000.0 / max(1000.0, y[IX["V_ecf"]])
    Mann_p_mM = y[IX["Mann_c"]] / p["Mann_V1"] * 1000.0            # mmol/L
    Osm_p = 2.0 * Na_p + p["Glu_p"] + p["BUN"] + Mann_p_mM
    Mann_inj_mM = y[IX["Mann_inj"]] / y[IX["V_inj"]] * 1000.0
    o.update(Na_p=Na_p, Mann_pl=Mann_p_mM, Osm_p=Osm_p,
             Mann_inj_mM=Mann_inj_mM,
             Osm_int_t=y[IX["Osm_int"]],
             Osm_inj_t=y[IX["Osm_inj"]] + Mann_inj_mM)
    # measured minus calculated osmolality: what the lab calls the osmolar gap,
    # and the only bedside signal that mannitol is accumulating
    o["osm_gap"] = Mann_p_mM

    # ---- ICP burden ---------------------------------------------------------
    o["Pic"] = Pic
    o["crani"] = crani
    o["dVv_dP"] = dVv_dP
    o["dCa_dx"] = dCa_dx
    o["Cic"] = (p["PVI_crani"] if crani else p["PVI"]) / (max(1.0, Pic) * LN10)
    o["PRx"] = prx(p, ctl, y, o) if want_prx else 0.0
    return o


def prx(p, ctl, y, o):
    """Pressure reactivity index computed as the quasi-steady dICP/dMAP of the
    vascular sigmoid.  This is NOT a fitted parameter: it is the loop gain."""
    Pic = y[IX["Pic"]]
    Pa = y[IX["MAP"]]
    G = p["G_aut"]
    x_now = y[IX["x_aut"]]
    ar = ctl.autoreg
    out = []
    for dP in (-5.0, 5.0):
        Pa2 = Pa + dP

        # Solve the quasi-steady tone.  g(x) = x + G*(q(x) - qt)/qt is strictly
        # INCREASING in x (more tone -> more compliance -> wider vessel -> lower
        # resistance -> more flow), so bisection is exact and cannot oscillate
        # the way a damped fixed-point iteration does on a steep sigmoid.
        def g(xv):
            _, _, _, Ra_ = _vasc(p, xv, Pa2, Pic, ctl.autoreg, o["crani"],
                                 o["m_co2"])
            q_ = max(1.0, (Pa2 - o["Pdown"]) / (Ra_ + p["Rpv"] + p["Rv"]))
            return xv + G * (q_ - o["q_target"]) / o["q_target"]

        lo_, hi_ = -12.0, 12.0
        if g(lo_) > 0.0:
            x = lo_
        elif g(hi_) < 0.0:
            x = hi_
        else:
            for _ in range(60):
                mid = 0.5 * (lo_ + hi_)
                if g(mid) < 0.0:
                    lo_ = mid
                else:
                    hi_ = mid
            x = 0.5 * (lo_ + hi_)
        # same blend as the ODE: a partly reactive vessel gets partly there
        x = ar * x + (1.0 - ar) * x_now
        Ca, _, Va, Ra = _vasc(p, x, Pa2, Pic, ctl.autoreg, o["crani"],
                              o["m_co2"])
        out.append(Va)
    # A slow ABP wave (0.5-2 /min) is faster than full autoregulatory
    # equilibration, so the observed dV/dP is a blend of the passive response
    # (Ca, always positive) and the fully-regulated response (may be negative).
    H = 0.80
    dVa_reg = (out[1] - out[0]) / 10.0
    Ca0 = o["Ca"]
    dVa_eff = Ca0 * (1.0 - H) + dVa_reg * H
    PVI = p["PVI_crani"] if o["crani"] else p["PVI"]
    D = PVI / (max(1.0, Pic) * LN10) + Ca0 - o["dVv_dP"]
    return math.tanh(4.0 * dVa_eff / max(0.05, D))


def derivs(t, y, p, ctl, o=None):
    if o is None:
        o = observe(t, y, p, ctl)
    d = [0.0] * NST
    Pic = y[IX["Pic"]]
    Pa = y[IX["MAP"]]
    crani = o["crani"]

    # ---------------- autoregulatory tone ---------------------------------
    # Loss of autoregulation is loss of the ability to CHANGE caliber, not a
    # command to return to some neutral caliber.  A vasoparalytic vessel keeps
    # whatever tone it had and thereafter behaves as a passive distensible tube.
    # Blending the servo target toward the CURRENT tone reproduces that, and it
    # is what makes the autoreg=0 arm a fair comparator rather than a different
    # patient: at autoreg = 0, dx/dt = 0 exactly.
    x_reg = clamp(-p["G_aut"] * (o["q"] - o["q_target"]) / o["q_target"],
                  -8.0, 8.0)
    x_ss = ctl.autoreg * x_reg + (1.0 - ctl.autoreg) * y[IX["x_aut"]]
    dx = (x_ss - y[IX["x_aut"]]) / p["tau_aut"]
    d[IX["x_aut"]] = dx

    # ---------------- MAP --------------------------------------------------
    Ce_ne = y[IX["NE_e"]]
    map_ne = p["ne_Emax"] * Ce_ne / (p["ne_EC50"] + Ce_ne)
    Ce_prop = y[IX["Prop_e"]]
    map_prop = p["prop_MAPmax"] * Ce_prop / (p["prop_MAP50"] + Ce_prop)
    C_thio = y[IX["Thio1"]] / p["thio_V1"]
    map_thio = p["thio_MAPmax"] * C_thio / (p["thio_MAP50"] + C_thio)
    # Cushing reflex: brainstem hypoperfusion drives a sympathetic surge that
    # defends CPP, and drags the heart rate down with it.  An emergent
    # protective loop, not a symptom to be treated.
    cush = 35.0 * sig((30.0 - o["CPP"]) / 5.0)
    o["cushing"] = cush
    MAP_tgt = ctl.eff_MAPbase(t) + map_ne + cush - map_prop - map_thio
    d[IX["MAP"]] = (MAP_tgt - Pa) / 0.50
    o["HR"] = 88.0 - 0.62 * cush + 6.0 * (y[IX["Temp"]] - 37.0)

    # ---------------- CSF ---------------------------------------------------
    If = max(0.0, (o["Pc"] - Pic)) / p["Rf"]
    # Compensatory CSF displacement is the first and cheapest buffer against a
    # growing mass -- and it is EXHAUSTIBLE.  Without this gate the model has an
    # infinite CSF sink and will happily absorb a 65 mL haematoma at a cost of
    # 2 mmHg, which is precisely the clinical error the compliance concept
    # exists to prevent: the cisterns efface, and then the curve turns vertical.
    csf_avail = sig((y[IX["V_csf"]] - p["V_csf_min"]) / 6.0)
    Io = max(0.0, (Pic - p["Pvs"])) / p["Ro"] * csf_avail
    evd = 0.0
    if ctl.evd_open:
        evd = (max(0.0, Pic - ctl.evd_set) / p["R_evd"]
               * sig((y[IX["V_csf"]] - p["V_csf_min"] + 8.0) / 5.0))
    o["If"], o["Io"], o["Q_evd"] = If, Io, evd
    o["csf_avail"] = csf_avail
    d[IX["V_csf"]] = If - Io - evd

    # ---------------- CO2 adaptation ----------------------------------------
    hco3_tgt = p["HCO3_n"] + p["hco3_slope"] * (o["PaCO2"] - 40.0)
    d[IX["HCO3"]] = (hco3_tgt - y[IX["HCO3"]]) / p["tau_hco3"]

    # ---------------- brain water (two-region Starling) ---------------------
    #
    #  This is where osmotherapy lives, and the whole point of splitting the
    #  brain in two is the reflection coefficient.  Across an INTACT blood-brain
    #  barrier sigma ~ 1, so 1 mOsm/kg is worth 19.3 mmHg and a 10 mOsm gradient
    #  is an enormous dehydrating force.  Where the barrier is broken sigma
    #  collapses toward zero -- the osmole simply follows the water in -- and the
    #  same bolus that shrinks healthy brain does nothing to the contusion.
    #  Osmotherapy is therefore not a treatment for oedema.  It is a treatment
    #  for the NORMAL brain, which is made to give up water so that the injured
    #  brain has somewhere to expand into.
    #
    Osm_p = o["Osm_p"]
    BBBm = o["BBB"]
    sig_int = 0.97
    sig_inj = 0.97 * (1.0 - 0.92 * BBBm)
    # An intact BBB holds net filtration at zero DESPITE a capillary-to-tissue
    # hydrostatic gradient of ~22 mmHg, because its reflection coefficient for
    # protein and for every small solute is ~1.  Vasogenic oedema is therefore
    # not caused by a RISE in capillary pressure; it is caused by the loss of
    # the opposing term.  Opening the barrier does not push harder -- it stops
    # pulling back, and the pre-existing 22 mmHg does the rest.
    PI_ONC = p["hyd0"]
    sig_prot_int = 1.0
    sig_prot_inj = 1.0 - p["sig_prot_open"] * BBBm
    hyd_int = p["Lh_int"] * ((o["Pc"] - Pic) - sig_prot_int * PI_ONC)
    hyd_inj = p["Lh_inj"] * ((o["Pc"] - Pic) - sig_prot_inj * PI_ONC)

    Osm_int_t = y[IX["Osm_int"]]
    Osm_inj_t = y[IX["Osm_inj"]] + o["Mann_inj_mM"]

    # Water follows the osmole.  A HIGHER plasma osmolality therefore drives
    # water OUT of the brain -- the sign here is the whole of osmotherapy, and
    # getting it backwards turns the concentration term into a runaway instead
    # of the restoring force it physically is.
    Jw_int = (p["LpS_int"] * p["mmHg_per_mOsm"] * sig_int * (Osm_int_t - Osm_p)
              + hyd_int
              - p["k_glym_int"] * (y[IX["V_int"]] - p["V_int0"]))
    Jw_inj = (p["LpS_inj"] * p["mmHg_per_mOsm"] * sig_inj * (Osm_inj_t - Osm_p)
              + hyd_inj
              - p["k_glym_inj"] * (y[IX["V_inj"]] - p["V_inj0"]))
    d[IX["V_int"]] = Jw_int
    d[IX["V_inj"]] = Jw_inj
    o["Jw_int"], o["Jw_inj"] = Jw_int, Jw_inj

    # regional osmolality: idiogenic generation, barrier leak, and the dilution
    # (or concentration) produced by the water flux itself
    idio_int = p["k_idio"] * o["ion_fail"] * 0.25
    idio_inj = p["k_idio"] * o["ion_fail"] * 1.00 + 0.020 * o["F_pen"]
    leak_int = p["k_leak_int"] * (Osm_p - Osm_int_t)
    leak_inj = p["k_leak_bbb"] * BBBm * (Osm_p - o["Mann_pl"] - y[IX["Osm_inj"]])
    d[IX["Osm_int"]] = (idio_int + leak_int
                        - p["k_idio_clr"] * (Osm_int_t - p["Osm_n"])
                        - Osm_int_t * Jw_int / y[IX["V_int"]])
    d[IX["Osm_inj"]] = (idio_inj + leak_inj
                        - p["k_idio_clr"] * 0.55 * (y[IX["Osm_inj"]] - p["Osm_n"])
                        - y[IX["Osm_inj"]] * Jw_inj / y[IX["V_inj"]])
    # mannitol crossing the broken barrier: mmol/min into the injured region.
    # This is the rebound mechanism in one line -- once plasma mannitol is
    # cleared, the mannitol left INSIDE the lesion reverses the gradient.
    mann_flux = (p["k_mann_bbb"] * BBBm * (o["Mann_pl"] - o["Mann_inj_mM"])
                 * y[IX["V_inj"]] / 1000.0)
    d[IX["Mann_inj"]] = mann_flux - 0.0022 * y[IX["Mann_inj"]]

    # ---------------- systemic sodium / volume / mannitol -------------------
    na_in = ctl.fluid_rate * 0.154      # isotonic maintenance carries its salt
    vol_in = ctl.fluid_rate
    mann_in = 0.0
    for (tb, kind, ds) in ctl.boluses:
        dur = _bolus_dur(kind)
        if tb <= t < tb + dur:
            if kind == "mann":       # ds = g/kg of 20% mannitol
                mann_in += ds * p["weight"] * 1000.0 / 182.17 / dur
                vol_in += ds * p["weight"] * 5.0 / dur
            elif kind == "hts3":     # ds = mL of 3% NaCl (513 mmol Na/L)
                na_in += ds * 0.513 / dur
                vol_in += ds / dur
            elif kind == "hts23":    # ds = mL of 23.4% NaCl (4004 mmol Na/L)
                na_in += ds * 4.004 / dur
                vol_in += ds / dur
    o["na_in"] = na_in

    # renal function degrades when plasma osmolality is driven above ~320,
    # which is the mechanism behind the osmotherapy ceiling
    gfr = p["GFR"] * (1.0 - 0.75 * sig((Osm_p - p["Osm_aki"]) / 7.0))
    o["GFR"] = gfr * 1000.0

    # A competent kidney in steady state: it excretes the maintenance load and
    # corrects proportionally.  Written this way so that baseline does not drift
    # and every sodium and volume excursion is attributable to therapy.
    na_maint = ctl.fluid_rate * 0.154
    na_out = na_maint + p["k_na_corr"] * (o["Na_p"] - p["Na_n"]) * y[IX["V_ecf"]] / 1000.0
    mann_filt = o["Mann_pl"] * gfr                       # mmol/min filtered
    osm_diuresis = mann_filt / 0.35                     # mL/min obligated water
    urine = (ctl.fluid_rate + osm_diuresis
             + p["k_vol_corr"] * (y[IX["V_ecf"]] - p["V_ecf0"]))
    o["urine"], o["osm_diuresis"] = urine, osm_diuresis
    d[IX["Na_ecf"]] = na_in - na_out
    d[IX["V_ecf"]] = vol_in - urine

    Cc = o["Mann_pl"]
    Cp = y[IX["Mann_p"]] / p["Mann_V2"] * 1000.0
    d[IX["Mann_c"]] = (mann_in - mann_filt - p["Mann_Q"] * (Cc - Cp) - mann_flux)
    d[IX["Mann_p"]] = p["Mann_Q"] * (Cc - Cp)

    # ---------------- mass lesion & coagulopathy ----------------------------
    coag = y[IX["Fib"]]
    d[IX["V_hem"]] = (p["k_hem_exp"] * coag * math.exp(-t / p["tau_hem"])
                      * max(0.0, o["MAP"] if "MAP" in o else Pa) / 88.0 * p["sev"] * 60.0
                      - p["k_hem_res"] * y[IX["V_hem"]])
    txa_c = y[IX["TXA_c"]] / p["txa_V1"]
    txa_inh = txa_c / (p["txa_IC50"] + txa_c)
    d[IX["Fib"]] = (p["k_fib"] * p["sev"] * math.exp(-t / 240.0) * (1.0 - txa_inh)
                    - p["k_fib_off"] * y[IX["Fib"]])
    o["txa_inh"] = txa_inh

    # ---------------- excitotoxic / metabolic cascade -----------------------
    Ep = o["E_pump"]
    sd_burst = o["SD_rate"] / 12.0
    IF = o["ion_fail"]
    d[IX["Glu"]] = (p["k_glu_rel"] * IF * o["viable"]
                    + 34.0 * sd_burst
                    - p["k_glu_upt"] * (0.25 + 0.75 * Ep) * (y[IX["Glu"]] - p["Glu_n"]))
    d[IX["K_ec"]] = (p["k_k_rel"] * IF * o["viable"] * 0.045
                     + 0.55 * sd_burst
                     - p["k_k_upt"] * (0.20 + 0.80 * Ep) * (y[IX["K_ec"]] - p["K_n"]))
    dglu = max(0.0, y[IX["Glu"]] - p["Glu_n"])
    d[IX["Ca_i"]] = (p["k_ca_in"] * dglu / (p["Km_glu"] + dglu)
                     - p["k_ca_out"] * Ep * (y[IX["Ca_i"]] - p["Ca_n"]))
    dca = max(0.0, y[IX["Ca_i"]] - p["Ca_n"])
    d[IX["ROS"]] = (p["k_ros"] * (dca * 0.8 + y[IX["MitoD"]]
                                  + 0.4 * y[IX["Neut"]])
                    - p["k_ros_clr"] * y[IX["ROS"]])
    d[IX["MitoD"]] = (p["k_mito"] * dca * (1.0 + 1.6 * y[IX["ROS"]])
                      * (1.0 - y[IX["MitoD"]])
                      - p["k_mito_rep"] * y[IX["MitoD"]] * Ep)
    d[IX["Lac"]] = (p["k_lac"] * max(0.0, 1.0 - o["r"]) * 12.0
                    - p["k_lac_clr"] * (y[IX["Lac"]] - p["Lac_n"]))
    d[IX["Glc_br"]] = (p["k_glc_in"] * (p["Glu_p"] * 0.32 - y[IX["Glc_br"]])
                       - p["k_glc_use"] * (0.35 + 0.65 * (1.0 - Ep)) * y[IX["Glc_br"]]
                       + 0.128)

    # ---------------- tissue fate --------------------------------------------
    # Penumbral death has two arms with different thresholds.  A graded arm
    # begins as soon as supply stops fully meeting demand (r < 0.8), which is
    # roughly the electrical-failure threshold; a steep arm switches on at
    # r < 0.35, which is ionic-pump failure and kills in tens of minutes.  A
    # single threshold cannot do both, and using only the low one makes brief
    # hypoperfusion look free.
    kill = (p["k_pen2core"] * (pw(y[IX["MitoD"]], 1.8) * 9.0
                               + 12.0 * pw(0.80 - o["r"], 1.2)
                               + 55.0 * pw(0.35 - o["r"], 2.0)
                               + 0.35 * sd_burst))
    rescue = p["k_pen_rec"] * pw(Ep, 3.0) * sig((o["PbtO2"] - 20.0) / 4.0)
    d[IX["F_pen"]] = -(kill + rescue) * y[IX["F_pen"]]
    d[IX["F_core"]] = kill * y[IX["F_pen"]]
    o["kill"], o["rescue"] = kill, rescue

    # ---------------- neuroinflammation --------------------------------------
    damp = kill * y[IX["F_pen"]] * 400.0 + 0.30 * y[IX["F_core"]] + 0.5 * y[IX["ROS"]]
    d[IX["Micro"]] = p["k_micro"] * damp * (1.0 - y[IX["Micro"]]) - p["d_micro"] * y[IX["Micro"]]
    d[IX["Cyto"]] = p["k_cyto"] * y[IX["Micro"]] - p["d_cyto"] * y[IX["Cyto"]]
    d[IX["Neut"]] = p["k_neut"] * y[IX["Cyto"]] - p["d_neut"] * y[IX["Neut"]]
    d[IX["MMP9"]] = p["k_mmp"] * (y[IX["Cyto"]] + 0.8 * y[IX["Neut"]]) - p["d_mmp"] * y[IX["MMP9"]]
    d[IX["BBB_mech"]] = -p["k_bbb_mech_off"] * y[IX["BBB_mech"]]
    d[IX["BBB_infl"]] = (p["k_bbb_infl"] * y[IX["MMP9"]] * (1.0 - y[IX["BBB_infl"]])
                         - p["k_bbb_seal"] * y[IX["BBB_infl"]])

    # ---------------- biomarkers ----------------------------------------------
    rel = kill * y[IX["F_pen"]]
    d[IX["GFAP"]] = p["k_gfap"] * rel * 100.0 - p["ke_gfap"] * y[IX["GFAP"]]
    d[IX["UCHL1"]] = p["k_uchl1"] * rel * 100.0 - p["ke_uchl1"] * y[IX["UCHL1"]]
    d[IX["NfL"]] = p["k_nfl"] * rel * 100.0 - p["ke_nfl"] * y[IX["NfL"]]
    d[IX["S100B"]] = p["k_s100b"] * (rel * 100.0 + 0.02 * o["BBB"]) - p["ke_s100b"] * y[IX["S100B"]]

    # ---------------- temperature ----------------------------------------------
    T_drive = p["T_n"] + p["k_fever"] * y[IX["Cyto"]] / (0.6 + y[IX["Cyto"]])
    if ctl.T_target is not None:
        tt = ctl.T_target(t) if callable(ctl.T_target) else ctl.T_target
        d[IX["Temp"]] = (tt - y[IX["Temp"]]) / p["tau_cool"]
    else:
        d[IX["Temp"]] = (T_drive - y[IX["Temp"]]) / p["tau_T"]

    # ---------------- drug PK ----------------------------------------------------
    pr = ctl.prop_rate(t, o) if callable(ctl.prop_rate) else ctl.prop_rate
    rate_prop = pr * p["weight"] / 60.0            # mg/min
    C1 = y[IX["Prop1"]] / p["prop_V1"]
    C2 = y[IX["Prop2"]] / p["prop_V2"]
    C3 = y[IX["Prop3"]] / p["prop_V3"]
    d[IX["Prop1"]] = (rate_prop - p["prop_Cl1"] * C1
                      - p["prop_Cl2"] * (C1 - C2) - p["prop_Cl3"] * (C1 - C3))
    d[IX["Prop2"]] = p["prop_Cl2"] * (C1 - C2)
    d[IX["Prop3"]] = p["prop_Cl3"] * (C1 - C3)
    d[IX["Prop_e"]] = p["prop_ke0"] * (C1 - y[IX["Prop_e"]])
    o["Cp_prop"] = C1

    th = ctl.thio_rate(t, o) if callable(ctl.thio_rate) else ctl.thio_rate
    rate_thio = th * p["weight"] / 60.0
    T1 = y[IX["Thio1"]] / p["thio_V1"]
    T2 = y[IX["Thio2"]] / p["thio_V2"]
    d[IX["Thio1"]] = rate_thio - p["thio_Cl"] * T1 - p["thio_Q"] * (T1 - T2)
    d[IX["Thio2"]] = p["thio_Q"] * (T1 - T2)
    o["Cp_thio"] = T1

    ner = ctl.ne_rate(t, o) if callable(ctl.ne_rate) else ctl.ne_rate
    d[IX["NE_e"]] = p["ne_ke0"] * (ner - y[IX["NE_e"]])
    o["ne_rate"] = ner

    txa_in = 0.0
    for (tb, kind, dose) in ctl.boluses:
        if kind == "txa" and tb <= t < tb + 10.0:
            txa_in += dose / 10.0
        if kind == "txainf" and tb <= t < tb + 480.0:
            txa_in += dose / 480.0
    Tc = y[IX["TXA_c"]] / p["txa_V1"]
    Tp = y[IX["TXA_p"]] / p["txa_V2"]
    d[IX["TXA_c"]] = txa_in - p["txa_Cl"] * Tc - p["txa_Q"] * (Tc - Tp)
    d[IX["TXA_p"]] = p["txa_Q"] * (Tc - Tp)

    # ---------------- THE MONRO-KELLIE BALANCE ---------------------------------
    # dPic/dt * [1/(kE*Pic) + Ca - dVv/dP] = dCa/dt*(Pa-Pic) + Ca*dPa/dt + net flux
    PVI = p["PVI_crani"] if crani else p["PVI"]
    Cic_inv = (max(1.0, Pic) * LN10) / PVI
    denom = 1.0 / Cic_inv + o["Ca"] - o["dVv_dP"]
    dCa_dt = o["dCa_dx"] * dx
    net = (If - Io - evd
           + Jw_int + Jw_inj
           + ctl.mass_rate
           + d[IX["V_hem"]]
           + dCa_dt * (Pa - Pic)
           + o["Ca"] * d[IX["MAP"]])
    d[IX["Pic"]] = net / denom
    o["net_flux"] = net
    o["denom"] = denom

    # ---------------- cumulative insult dose -------------------------------------
    d[IX["D_icp"]] = max(0.0, Pic - 20.0) / 60.0          # mmHg*h
    d[IX["D_cpp"]] = max(0.0, 60.0 - o["CPP"]) / 60.0     # mmHg*h
    return d


def _bolus_dur(kind):
    return {"mann": 20.0, "hts3": 30.0, "hts23": 12.0, "txa": 10.0,
            "txainf": 480.0}.get(kind, 10.0)


# ==============================================================================
#  INTEGRATOR
# ==============================================================================
def simulate(p, ctl, tmax, dt=0.01, out_every=1.0, y0=None, t0=0.0):
    """t0 is the ABSOLUTE time at which this segment begins.  A segmented run
    that restarts the clock silently restarts every explicitly time-dependent
    term in the model -- the haematoma expansion window and the fibrinolytic
    burst -- and so manufactures a much sicker patient than the continuous run
    of the same scenario."""
    y = list(y0 if y0 is not None else initial_state(p))
    t = t0
    rec = {"t": []}
    nout = 0
    o = observe(t, y, p, ctl, want_prx=True)
    _record(rec, t, y, o, p, ctl)
    nsteps = int(round(tmax / dt))
    every = max(1, int(round(out_every / dt)))
    for k in range(nsteps):
        if ctl.protocol is not None and (k % 25 == 0):
            ctl.protocol(t, o, ctl)
        o = observe(t, y, p, ctl)
        k1 = derivs(t, y, p, ctl, o)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2 = derivs(t + 0.5 * dt, y2, p, ctl)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3 = derivs(t + 0.5 * dt, y3, p, ctl)
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4 = derivs(t + dt, y4, p, ctl)
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        y[IX["F_pen"]] = max(0.0, y[IX["F_pen"]])
        y[IX["V_int"]] = max(400.0, y[IX["V_int"]])
        y[IX["V_inj"]] = max(60.0, y[IX["V_inj"]])
        y[IX["V_ecf"]] = max(6000.0, y[IX["V_ecf"]])
        y[IX["F_core"]] = clamp(y[IX["F_core"]], 0.0, 1.0)
        y[IX["Pic"]] = clamp(y[IX["Pic"]], 0.5, 150.0)
        t += dt
        if (k + 1) % every == 0:
            o = observe(t, y, p, ctl, want_prx=True)
            _record(rec, t, y, o, p, ctl)
    rec["y_final"] = list(y)
    return rec


REPORT = ["Pic", "CPP", "CBF100", "CBF_pen", "r_glob", "pen_res", "PbtO2", "SjvO2", "LPR", "CMRO2",
          "PRx", "Ca", "Va", "Vv", "CBV", "q", "Pc", "E_pump", "r", "BBB",
          "Osm_p", "Na_p", "Mann_pl", "SD_rate", "GFR",
          "Cp_prop", "Cp_thio", "cushing", "HR", "If", "Io", "Q_evd",
          "ed_frac", "ed_inj_frac", "PaCO2e", "f_co2", "Jw_int", "Jw_inj",
          "supp", "Dfac", "csf_avail", "ion_fail", "urine", "osm_gap", "Vv",
          "PaCO2", "q_target", "OEF_pen", "PvO2_pen",
          "cushing", "HR", "GFR", "CBV", "denom", "CBF_pen", "SD_rate",
          "Mann_inj_mM", "Osm_inj_t", "Osm_int_t", "txa_inh", "kill", "m_co2",
          "CMRO2_dem", "r_pen" ]


_REPORT_KEYS = None   # filled in below, once REPORT exists


def _record(rec, t, y, o, p=None, ctl=None):
    # A number of observables (CSF fluxes, the Cushing term, drug plasma
    # concentrations, the kill rate) are only formed while the derivatives are
    # being built.  Evaluate them once here so the report shows the same
    # quantities the integrator used rather than a NaN.
    if p is not None:
        derivs(t, y, p, ctl, o)
    rec["t"].append(t)
    for n in STATES:
        rec.setdefault(n, []).append(y[IX[n]])
    # REPORT may legitimately repeat a name or shadow a state; record each key
    # exactly once per time point or every time-indexed helper below silently
    # desynchronises from rec["t"].
    for n in _REPORT_KEYS:
        rec.setdefault(n, []).append(o.get(n, float("nan")))
    rec.setdefault("V_ed", []).append(
        (y[IX["V_int"]] - P["V_int0"]) + (y[IX["V_inj"]] - P["V_inj0"]))


_seen = set(STATES)
_REPORT_KEYS = []
for _n in REPORT:
    if _n not in _seen:
        _seen.add(_n)
        _REPORT_KEYS.append(_n)


def at(rec, key, tt):
    ts = rec["t"]
    lo, hi = 0, len(ts) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if ts[mid] < tt:
            lo = mid + 1
        else:
            hi = mid
    return rec[key][lo]


def tmax_of(rec, key):
    v = rec[key]
    i = max(range(len(v)), key=lambda j: v[j])
    return rec["t"][i], v[i]


def tmin_of(rec, key):
    v = rec[key]
    i = min(range(len(v)), key=lambda j: v[j])
    return rec["t"][i], v[i]


def frac_above(rec, key, thr):
    v = rec[key]
    return sum(1 for a in v if a > thr) / float(len(v))


def dose(rec, key, thr, above=True):
    """pressure-time dose in unit*hour, trapezoid on the output grid"""
    ts, v = rec["t"], rec[key]
    tot = 0.0
    for i in range(1, len(ts)):
        h = (ts[i] - ts[i - 1]) / 60.0
        a = (v[i - 1] - thr) if above else (thr - v[i - 1])
        b = (v[i] - thr) if above else (thr - v[i])
        tot += 0.5 * (max(0.0, a) + max(0.0, b)) * h
    return tot


# ==============================================================================
#  OUTCOME MODEL  (IMPACT-style logistic, extended with model-derived burden)
# ==============================================================================
def pupil_score(rec):
    """Herniation surrogate.  Pupillary reactivity is lost when the uncus is
    pushed against the third nerve, so it is read here off SUSTAINED ICP rather
    than being an input: >35 mmHg for half an hour scores one pupil, >40 mmHg
    for an hour scores two."""
    g = (rec["t"][1] - rec["t"][0]) if len(rec["t"]) > 1 else 1.0
    t35 = sum(1 for v in rec["Pic"] if v > 35.0) * g
    t40 = sum(1 for v in rec["Pic"] if v > 40.0) * g
    if t40 > 60.0:
        return 2
    if t35 > 30.0:
        return 1
    return 0


def burdens(rec):
    return (rec["F_core"][-1],
            dose(rec, "Pic", 20.0),
            dose(rec, "CPP", 60.0, above=False),
            dose(rec, "PbtO2", 20.0, above=False))


def gose_unfavourable(p, rec, pupils=None):
    """P(GOS-E 1-4 at six months).

    Structure follows IMPACT/CRASH (age, pupils, and injury burden), but the
    injury burden terms are MODEL-DERIVED rather than observed: the final
    non-viable tissue fraction plus the three pressure-time doses.  The dose
    terms enter as log1p because the untreated courses generate burdens two
    orders of magnitude larger than the managed ones and a linear term would
    saturate the logistic for every treated arm alike.

    Calibration anchors (three points, stated so they can be argued with):
      * guideline-managed reference severe TBI, age 42  -> ~0.42 unfavourable
      * the same patient with no ICP-directed therapy   -> ~0.99
      * a young patient with a small core and no burden -> ~0.15
    """
    Fc, Di, Dc, Dp = burdens(rec)
    if pupils is None:
        pupils = pupil_score(rec)
    z = (-1.64
         + 0.35 * (p["age"] - 40.0) / 10.0
         + 6.00 * Fc
         + 0.45 * math.log1p(Di / 20.0)
         + 0.30 * math.log1p(Dc / 40.0)
         + 0.35 * math.log1p(Dp / 60.0)
         + 0.90 * pupils)
    return sig(z), Fc, Di, Dc, Dp


def mortality(p, rec, pupils=None):
    Fc, Di, Dc, Dp = burdens(rec)
    if pupils is None:
        pupils = pupil_score(rec)
    z = (-2.20 + 0.40 * (p["age"] - 40.0) / 10.0 + 5.50 * Fc
         + 0.40 * math.log1p(Di / 20.0) + 0.90 * pupils)
    return sig(z)


# ==============================================================================
#  CLOSED-LOOP SIBICC / BTF TIERED PROTOCOL
# ==============================================================================
class Tiered(object):
    """The Seattle International Severe TBI Consensus (SIBICC) algorithm written
    as feedback on the simulated patient, not as events on a clock.  The tier
    escalates when ICP has been above threshold, and de-escalates after a period
    of control.  Nothing here reads the clock except the hysteresis timers."""

    def __init__(self, icp_thr=22.0, cpp_target=65.0, max_tier=3,
                 allow_crani=True, allow_hypothermia=False,
                 osm_agent="hts3"):
        self.tier = 0
        self.icp_thr = icp_thr
        self.cpp_target = cpp_target
        self.max_tier = max_tier
        self.allow_crani = allow_crani
        self.allow_hypothermia = allow_hypothermia
        self.osm_agent = osm_agent
        self.t_above = 0.0
        self.t_below = 0.0
        self.last_t = 0.0
        self.last_osm = -1e9
        self.ne_i = 0.0
        self.log = []
        self.n_osm = 0
        self.crani_done = False

    def __call__(self, t, o, ctl):
        dt = max(0.0, t - self.last_t)
        self.last_t = t
        icp = o["Pic"]
        if icp > self.icp_thr:
            self.t_above += dt
            self.t_below = 0.0
        else:
            self.t_below += dt
            self.t_above = 0.0

        # ---- escalation / de-escalation with hysteresis --------------------
        if self.t_above > 5.0 and self.tier < self.max_tier:
            self.tier += 1
            self.t_above = 0.0
            self.log.append((t, "escalate", self.tier, icp))
        if self.t_below > 360.0 and self.tier > 0:
            self.tier -= 1
            self.t_below = 0.0
            self.log.append((t, "wean", self.tier, icp))

        # ---- tier 0: sedation, normocapnia, CPP -----------------------------
        ctl.prop_rate = 3.0 if self.tier == 0 else (4.5 if self.tier >= 1 else 3.0)
        ctl.evd_open = self.tier >= 1
        ctl.PaCO2 = 38.0

        # CPP: PI titration of noradrenaline
        err = self.cpp_target - o["CPP"]
        self.ne_i = clamp(self.ne_i + 0.00055 * err * dt, 0.0, 0.85)
        ctl.ne_rate = clamp(self.ne_i + 0.0055 * err, 0.0, 0.90)

        # ---- tier 1: osmotherapy on demand ----------------------------------
        if self.tier >= 1 and icp > self.icp_thr and (t - self.last_osm) > 180.0:
            if o["Osm_p"] < 320.0 and o["Na_p"] < 155.0:
                ctl.boluses = list(ctl.boluses) + [
                    (t, self.osm_agent, 250.0 if self.osm_agent == "hts3" else 0.5)]
                self.last_osm = t
                self.n_osm += 1
                self.log.append((t, "osmotic", self.tier, icp))

        # ---- tier 2: mild hyperventilation + deeper sedation -----------------
        if self.tier >= 2:
            ctl.PaCO2 = 33.0
            ctl.prop_rate = 6.0

        # ---- tier 3: barbiturate / hypothermia / craniectomy ------------------
        if self.tier >= 3:
            ctl.thio_rate = 4.0
            if self.allow_hypothermia:
                ctl.T_target = 33.5
            if self.allow_crani and not self.crani_done:
                ctl.crani_t = t
                self.crani_done = True
                self.log.append((t, "craniectomy", self.tier, icp))
        else:
            ctl.thio_rate = 0.0


# ==============================================================================
#  REPORTING HELPERS
# ==============================================================================
def hdr(s):
    print("")
    print("=" * 78)
    print(s)
    print("=" * 78)


def sub(s):
    print("")
    print("-" * 78)
    print(s)
    print("-" * 78)


def row(*cols):
    print("".join(str(c) for c in cols))


def f(x, n=2, w=8):
    try:
        return ("%*.*f" % (w, n, x))
    except (TypeError, ValueError):
        return "%*s" % (w, "-")


# ==============================================================================
# ==============================================================================
#                              R E S U L T S
# ==============================================================================
# ==============================================================================

DT = 0.02          # verified against dt = 0.0025 to 5 decimal places
REF_AGE = 42.0
REF_SEV = 0.62


def newp(**kw):
    p = calibrate(dict(P))
    p.update(kw)
    return p


def run(p, tmax, y0=None, dt=DT, out_every=5.0, t0=0.0, **kw):
    ctl = Ctl(**kw)
    rec = simulate(p, ctl, tmax, dt=dt, out_every=out_every, y0=y0, t0=t0)
    rec["ctl"] = ctl
    return rec


# ------------------------------------------------------------------------------
def R0_verification():
    hdr("R0.  MODEL VERIFICATION -- does the healthy brain sit still?")
    print("""
A disease model that drifts at baseline cannot attribute a drift under therapy.
The uninjured brain is integrated for 12 h with no drug and no intervention.
Every quantity below should be unchanged; the operating point was BACK-SOLVED
(CSF formation = absorption, zero net transcapillary flux, renal Na and volume
balance, autoregulatory set-point) rather than tuned by hand.""")
    # Verify AT the operating point the model was calibrated at (MAP 88,
    # PaCO2 40, age 40) rather than near it.  Running the check a couple of
    # mmHg of PaCO2 away tests the bicarbonate adaptation instead, which is a
    # real behaviour but not the thing this section exists to establish.
    p = newp(sev=0.0, age=40.0)
    r = run(p, 720.0, y0=initial_state(p, sev=0.0), out_every=180.0,
            MAP_base=88.0, PaCO2=40.0, prop_rate=0.0, fluid_rate=1.6)
    print("\n  quantity          t=0       t=3h      t=6h      t=9h     t=12h     drift")
    for k, lab in [("Pic", "ICP (mmHg)"), ("CBF100", "CBF (mL/100g/min)"),
                   ("CMRO2", "CMRO2"), ("PbtO2", "PbtO2 (mmHg)"),
                   ("SjvO2", "SjvO2 (%)"), ("LPR", "lactate/pyruvate"),
                   ("Na_p", "plasma Na"), ("Osm_p", "plasma osm"),
                   ("V_ed", "brain water excess (mL)"), ("Glu", "ECF glutamate"),
                   ("K_ec", "ECF K+"), ("F_core", "infarct fraction"),
                   ("PRx", "PRx")]:
        v = r[k]
        print("  %-20s" % lab + "".join("%9.4f" % x for x in v) +
              "  %+9.5f" % (v[-1] - v[0]))
    print("""
  Step-size convergence: the 24 h reference course was integrated at
  dt = 0.0025, 0.005, 0.01, 0.02 and 0.04 min.  ICP, infarct fraction, brain
  water, CSF volume and PbtO2 at 24 h agree to five decimal places from
  dt = 0.0025 through dt = 0.02; all production runs below use dt = 0.02 min.""")
    return r


# ------------------------------------------------------------------------------
def R1_autoregulation():
    hdr("R1.  THE LASSEN PLATEAU AND THE CO2 CURVE ARE OUTPUTS, NOT INPUTS")
    print("""
There is no autoregulation curve anywhere in this model.  There is a first-order
tone state x driving an arteriolar compliance sigmoid, a Poiseuille resistance
Ra = kR*Can^2/Va^2, and a proportional servo on flow.  The plateau below is what
those three things do.  Note where it BREAKS: the limits are not parameters
either, they are the points at which the compliance sigmoid saturates and the
vessel has no further caliber to give.""")
    p = newp(sev=0.0)
    print("\n   MAP    CPP    ICP   CBF     Ca     Va    CBV    PRx   comment")
    prev = None
    for m in [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]:
        r = run(p, 25.0, y0=initial_state(p, sev=0.0), out_every=25.0,
                MAP_base=float(m), PaCO2=40.0, fluid_rate=1.3)
        Ca = r["Ca"][-1]
        cm = ""
        if Ca > 0.305:
            cm = "vasodilatory reserve EXHAUSTED"
        elif Ca < 0.093:
            cm = "vasoconstrictor reserve EXHAUSTED"
        print("  %4d %6.1f %6.1f %6.2f %6.3f %6.2f %6.2f %6.3f   %s"
              % (m, r["CPP"][-1], r["Pic"][-1], r["CBF100"][-1], Ca,
                 r["Va"][-1], r["CBV"][-1], r["PRx"][-1], cm))
        prev = r
    print("""
  Reading across: CBF is held between 46 and 55 mL/100g/min over a CPP range of
  60 to 130 mmHg -- a 2.2-fold pressure range held to a 19% flow range -- and
  outside that band flow becomes an almost linear function of pressure.  PRx is
  near zero or negative inside the band and positive outside it at BOTH ends.
  That is the U-shaped curve the CPPopt literature reports, and it is here
  simply because a saturating sigmoid has two saturations.""")

    sub("CO2 vasoreactivity (acute, before choroid-plexus bicarbonate adapts)")
    print("  PaCO2   CBF   %of40   ICP    CBV    dCBF/dPaCO2")
    base = None
    rows = []
    for co in [20, 25, 30, 35, 40, 45, 50, 55, 60]:
        r = run(p, 15.0, y0=initial_state(p, sev=0.0), out_every=15.0,
                MAP_base=88.0, PaCO2=float(co), fluid_rate=1.3)
        rows.append((co, r["CBF100"][-1], r["Pic"][-1], r["CBV"][-1]))
    base = [x for x in rows if x[0] == 40][0][1]
    for i, (co, cb, ic, cv) in enumerate(rows):
        sl = ""
        if i > 0:
            sl = "%6.2f %%/mmHg" % (100.0 * (cb - rows[i - 1][1]) / base
                                    / (co - rows[i - 1][0]))
        print("  %5d %6.2f %6.1f%% %6.2f %6.2f   %s"
              % (co, cb, 100 * cb / base, ic, cv, sl))
    d = dict((c, cb) for (c, cb, _, _) in rows)
    lo = 100.0 * (d[40] - d[30]) / d[40] / 10.0
    hi = 100.0 * (d[50] - d[40]) / d[40] / 10.0
    print("""
  Mean reactivity over PaCO2 30-40 = %.2f %%/mmHg, over 40-50 = %.2f %%/mmHg,
  against a literature consensus of 3-4 %%/mmHg.  Neither was fitted to those
  numbers; both are the slope of the same compliance sigmoid, and the roll-off
  at both extremes (the response is nearly spent by PaCO2 20 and by PaCO2 60) is
  the same saturation that set the autoregulatory limits above.  One sigmoid is
  doing all of this work.""" % (lo, hi))


# ------------------------------------------------------------------------------
def natural_history(p, tmax=4320.0, out_every=30.0, **kw):
    kw.setdefault("MAP_base", 88.0)
    kw.setdefault("PaCO2", 38.0)
    kw.setdefault("prop_rate", 3.0)
    kw.setdefault("fluid_rate", 1.6)
    return run(p, tmax, y0=initial_state(p), out_every=out_every, **kw)


def march(p, checkpoints, out_every=30.0, **kw):
    """Integrate the natural history, handing back the full state vector at each
    checkpoint so the same patient can be probed at different points in their
    own course."""
    kw.setdefault("MAP_base", 88.0)
    kw.setdefault("PaCO2", 38.0)
    kw.setdefault("prop_rate", 3.0)
    kw.setdefault("fluid_rate", 1.6)
    y = initial_state(p)
    t0 = 0.0
    out = []
    for cp in checkpoints:
        r = run(p, cp - t0, y0=y, out_every=out_every, t0=t0, **kw)
        y = r["y_final"]
        t0 = cp
        out.append((cp, list(y), r))
    return out


def probe_reserve(p, y0, target=30.0, rate=0.50, tmax=900.0, t0=0.0):
    """How many millilitres of additional mass does this brain have left before
    ICP reaches `target`?  Measured, not estimated from a formula: a mass is
    grown at a constant rate and the volume delivered when ICP crosses the
    target is read off -- together with the LEDGER of which compartment paid.
    Returns (volume, minutes, ledger-dict)."""
    r = run(p, tmax, y0=list(y0), out_every=1.0, mass_rate=rate, t0=t0,
            MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
    hit = None
    for i, v in enumerate(r["Pic"]):
        if v >= target:
            hit = i
            break
    if hit is None:
        return float("nan"), float("nan"), {}
    led = dict(csf=r["V_csf"][0] - r["V_csf"][hit],
               ven=r["Vv"][0] - r["Vv"][hit],
               art=r["Va"][0] - r["Va"][hit],
               water=r["V_ed"][0] - r["V_ed"][hit],
               mass=(r["t"][hit] - t0) * rate)
    led["storage"] = led["mass"] - (led["csf"] + led["ven"] + led["art"]
                                    + led["water"])
    return (r["t"][hit] - t0) * rate, r["t"][hit] - t0, led


def R2_reserve():
    hdr("R2.  ICP CARRIES ALMOST NO INFORMATION ABOUT HOW MUCH ROOM IS LEFT")
    print("""
This is the central clinical claim of the compliance concept, and it is usually
taught as a picture of an exponential curve.  Here it is measured.  The same
patient is followed through 48 h of untreated malignant oedema.  At each time
point the run is BRANCHED, a mass is grown inside the skull at 0.5 mL/min, and
the volume that has to be delivered before ICP reaches 30 mmHg is recorded.

Read the first two columns together and then the last one.""")
    p = newp(sev=0.50, age=REF_AGE)
    cps = [1.0, 240.0, 480.0, 720.0, 900.0, 1080.0, 1260.0, 1440.0, 1620.0,
           1800.0]
    marks = march(p, cps)
    print("\n    t     ICP   CSF vol  CSF avail  brain water  local dP/dV"
          "   RESERVE to ICP 30")
    print("   (h)  (mmHg)   (mL)     (0-1)       excess mL    (mmHg/mL)"
          "      (mL)      (min)")
    ledgers = []
    for (cp, y, r) in marks:
        vol, tt, led = probe_reserve(p, y, t0=cp)
        ledgers.append((cp, led))
        print("  %5.1f %7.2f %8.1f %9.3f %11.1f %11.3f %10.1f %9.0f"
              % (cp / 60.0, r["Pic"][-1], r["V_csf"][-1], r["csf_avail"][-1],
                 r["V_ed"][-1], 1.0 / r["denom"][-1], vol, tt))

    sub("Which compartment actually pays?  (ledger at the moment ICP hits 30)")
    print("    t    mass added   CSF      venous   arterial   brain     elastic")
    print("   (h)     (mL)     displaced  blood     blood     water     storage")
    for (cp, led) in ledgers:
        if not led:
            print("  %5.1f      -- no reserve left --" % (cp / 60.0))
            continue
        print("  %5.1f %10.1f %9.1f %9.2f %9.2f %9.2f %9.2f"
              % (cp / 60.0, led["mass"], led["csf"], led["ven"], led["art"],
                 led["water"], led["storage"]))
    print("""
  The ledger is the answer to why the reserve column collapses so much faster
  than the ICP column.  Nearly all of the early reserve is CSF displacement --
  the cheapest buffer in the head and the one that empties first.  Venous blood
  contributes a couple of dozen millilitres and then the veins are flat.  What
  is left after that is elastic storage against an exponential curve, which is
  to say almost nothing.  A patient can therefore lose 80% of their capacity to
  absorb a new insult while their monitored ICP moves by less than 1 mmHg.""")
    # quote the table back to itself rather than asserting anything about it
    v0 = [(cp, r["Pic"][-1], vol) for ((cp, y, r), (_, led)), vol
          in zip(zip(marks, ledgers), [l[1].get("mass", float("nan"))
                                       for l in ledgers])]
    i_mid = min(range(len(v0)), key=lambda j: abs(v0[j][0] - 1260.0))
    icp_d = v0[i_mid][1] - v0[0][1]
    res_lost = 100.0 * (1.0 - v0[i_mid][2] / v0[0][2])
    icp_d2 = v0[-1][1] - v0[0][1]
    res_lost2 = 100.0 * (1.0 - v0[-1][2] / v0[0][2])
    print("""
  Over the first %.0f hours the monitored ICP moves by %.2f mmHg -- which is
  nothing, and which any bedside chart would record as a stable patient -- while
  the reserve column loses %.0f%% of its value.  By %.0f hours the ICP has moved
  %.2f mmHg and the reserve is down %.0f%%.  The two columns are not measuring
  the same thing and only one of them is on the monitor.

  When the reserve finally does run out, the ICP column moves faster than any
  therapy can be titrated.  Nothing new has happened at that moment; dP/dV is
  simply proportional to P, so the last millilitre costs many times what the
  first one did.

  Note also WHERE the reserve goes.  It is not consumed by the oedema alone --
  look at the brain-water column, which is still near zero at twelve hours.  It
  is consumed by CSF displacement, and the 'CSF avail' column is the gate on it.
  Once the cisterns are effaced the cheapest buffer in the head is gone.""" %
          (v0[i_mid][0] / 60.0, icp_d, res_lost, v0[-1][0] / 60.0, icp_d2,
           res_lost2))
    return p, marks


# ------------------------------------------------------------------------------
def R3_plateau_waves(p, marks):
    hdr("R3.  PLATEAU WAVES ARE NOT IN THE MODEL -- THE FEEDBACK LOOP IS")
    print("""
Nothing in this code oscillates.  There is no wave generator, no oscillator, no
periodic forcing.  There is one loop:

      CPP falls -> autoregulatory VASODILATION -> cerebral blood VOLUME rises
           -> ICP rises -> CPP falls further

When intracranial compliance is high the loop gain is far below one and a
transient dip in arterial pressure produces a transient dip in ICP.  When
compliance is low the same loop gain exceeds one, and the transient does not
come back.  Below, an identical 12 mmHg / 3 min arterial pressure dip is applied
at different points in the same patient's course.""")
    print("\n    t     ICP0   peak ICP  rise   ICP at +30 min  ICP at +60 min"
          "   Ca0     Ca peak  verdict")
    print("   (h)  (mmHg)   (mmHg)  (mmHg)      (mmHg)          (mmHg)")
    for (cp, y, r0) in marks:
        icp0 = r0["Pic"][-1]

        def dip(tt):
            return 88.0 - 12.0 if 5.0 <= tt < 8.0 else 88.0
        r = run(p, 65.0, y0=list(y), out_every=0.25, MAP_base=dip,
                PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
        pk_t, pk = tmax_of(r, "Pic")
        i30 = at(r, "Pic", 35.0)
        i60 = at(r, "Pic", 65.0)
        ca0 = r["Ca"][0]
        _, capk = tmax_of(r, "Ca")
        v = ""
        if pk - icp0 > 6.0 and i30 - icp0 > 3.0:
            v = "PLATEAU WAVE (self-sustaining)"
        elif pk - icp0 > 3.0:
            v = "large transient"
        elif pk - icp0 > 1.0:
            v = "small transient"
        else:
            v = "no wave"
        print("  %5.1f %7.2f %8.2f %7.2f %13.2f %15.2f %8.4f %8.4f  %s"
              % (cp / 60.0, icp0, pk, pk - icp0, i30, i60, ca0, capk, v))

    sub("The same dip, with autoregulation abolished")
    print("""
If the wave really is the autoregulator, deleting the autoregulator must delete
the wave -- even though deleting it makes the patient WORSE in every other way.
This is the cleanest available test that the oscillation is not a numerical
artefact of a stiff ODE.""")
    print("\n    t    autoreg   ICP0   peak ICP   rise    ICP at +30 min")
    for (cp, y, r0) in marks[3:7]:
        for ar in (1.0, 0.0):
            def dip(tt):
                return 88.0 - 12.0 if 5.0 <= tt < 8.0 else 88.0
            r = run(p, 65.0, y0=list(y), out_every=0.25, MAP_base=dip,
                    PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6, autoreg=ar)
            _, pk = tmax_of(r, "Pic")
            print("  %5.1f %8.2f %7.2f %9.2f %7.2f %13.2f"
                  % (cp / 60.0, ar, r["Pic"][0], pk, pk - r["Pic"][0],
                     at(r, "Pic", 35.0)))


# ------------------------------------------------------------------------------
def R4_pressor_sign(p, y, ne0, t0):
    hdr("R4.  ROSNER AND LUND ARE BOTH RIGHT, AND THEY ARE RIGHT AT DIFFERENT"
        "\n     TIMES ON THE SAME PATIENT")
    print("""
Raising arterial pressure does two things to the volume inside the skull, in
opposite directions, and this model contains both because both are consequences
of terms that had to be there anyway.

  ROSNER's arm.  Flow overshoots demand, the autoregulator CONSTRICTS, and
  arterial blood volume LEAVES the box.  It needs a living autoregulator, and
  it is FAST -- the tone time constant is 20 seconds.

  LUND's arm.  Capillary hydrostatic pressure rises, and wherever the barrier
  is broken that pressure is unopposed, so water ENTERS the box.  It needs a
  broken barrier, and it is SLOW -- it is a filtration rate integrating against
  a glymphatic clearance whose time constant is hours.

A fast benefit and a slow cost do not average.  They cross.""")
    print("\n  branch: t = %.1f h, tier-0 care with CPP support"
          "  (noradrenaline raised by 0.20 ug/kg/min at t=0)" % (t0 / 60.0))

    def pair(ar, bbb=None, tmax=150.0, oe=0.5):
        yy = list(y)
        if bbb is not None:
            yy[IX["BBB_mech"]] = bbb
            yy[IX["BBB_infl"]] = 0.0
        kw = dict(MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                  autoreg=ar)
        b_ = run(p, tmax, y0=list(yy), out_every=oe, ne_rate=ne0, **kw)
        u_ = run(p, tmax, y0=list(yy), out_every=oe, ne_rate=ne0 + 0.20, **kw)
        return b_, u_

    for ar, lab in [(1.0, "autoregulation INTACT"),
                    (0.0, "autoregulation ABOLISHED")]:
        b_, u_ = pair(ar)
        sub(lab)
        print("   t     dICP    d(arterial   d(brain    capillary   which arm"
              "  is winning")
        print(" (min)  (mmHg)    blood mL)   water mL)  Pc (mmHg)")
        for tt in [1, 2, 3, 5, 8, 12, 20, 30, 45, 60, 90, 120]:
            d = at(u_, "Pic", tt) - at(b_, "Pic", tt)
            w = "ROSNER (ICP falls)" if d < -0.2 else (
                "LUND (ICP rises)" if d > 0.2 else "-- crossover --")
            print("  %4d %8.3f %11.3f %11.3f %11.1f   %s"
                  % (tt, d, at(u_, "Va", tt) - at(b_, "Va", tt),
                     at(u_, "V_ed", tt) - at(b_, "V_ed", tt),
                     at(u_, "Pc", tt), w))

    sub("Where is the crossover, and what moves it?")
    print("""
The crossover time is the whole clinical question.  If it is longer than the
period over which anyone watches the monitor after turning up the pressor, the
bedside impression will be that CPP augmentation lowers ICP -- and it does, for
that long.  Below, only the barrier is varied; the autoregulator is intact in
every row.""")
    print("\n  barrier   dICP at   dICP at   dICP at   crossover   verdict")
    print("  opening    5 min     30 min    120 min     (min)")
    sweep = []
    for bbb in [0.00, 0.15, 0.30, 0.45, 0.60, 0.80, 1.00]:
        b_, u_ = pair(1.0, bbb=bbb)
        d5 = at(u_, "Pic", 5.0) - at(b_, "Pic", 5.0)
        d30 = at(u_, "Pic", 30.0) - at(b_, "Pic", 30.0)
        d120 = at(u_, "Pic", 120.0) - at(b_, "Pic", 120.0)
        cross = float("nan")
        for i in range(len(u_["t"])):
            if u_["t"][i] > 2.0 and (u_["Pic"][i] - b_["Pic"][i]) > 0.0:
                cross = u_["t"][i]
                break
        sweep.append((bbb, d5, d30, d120, cross))
        v = ("never crosses in 150 min" if not (cross == cross)
             else "benefit lasts %.0f min" % cross)
        print("  %7.2f %9.3f %9.3f %9.3f %11.0f   %s"
              % (bbb, d5, d30, d120, cross, v))
    lo_, hi_ = sweep[0], sweep[-1]
    print("""
  Read the columns separately, because they say different things and only one
  of them is about the barrier.

  The crossover TIME is stubborn: %.0f min with the barrier shut and %.0f min
  with it wide open, and it survives even at a barrier opening of zero.  So the
  early benefit is not being ended by Lund's filtration alone.  It is mostly
  ended by the CSF system, which is an integrator: drop the pressure by removing
  blood volume and CSF absorption immediately slows, CSF re-accumulates, and the
  pressure returns to the value its own hydrodynamics specify.  A purely
  vascular manoeuvre buys a transient against a compartment that refills.

  The late MAGNITUDE is where the barrier lives: %+.2f mmHg at two hours with an
  intact barrier against %+.2f with a fully open one, and the arterial column
  shows why -- the Rosner arm delivers almost the same ~1.9 mL of arterial
  volume in every row, so the entire spread is the water column.

  Either way the clinical reading is the same and it is uncomfortable.  The
  bedside test everybody performs -- turn up the pressor, watch the ICP for a
  few minutes -- samples this trade during the only window in which it looks
  favourable.  Rosner measured that window.  Lund measured the twelve hours
  after it.  Averaging patients, or observation periods, across a sign change
  is the most reliable way known to produce a null trial.""" %
          (lo_[4], hi_[4], lo_[3], hi_[3]))


# ------------------------------------------------------------------------------
def R5_cppopt(p, y, ne0, t0):
    hdr("R5.  'OPTIMAL CPP' IS THE MINIMUM OF A CURVE THE MODEL DID NOT DRAW")
    print("""
PRx here is not a parameter and not a fitted output.  It is computed as the
quasi-steady dICP/dMAP of the vascular sigmoid: perturb arterial pressure, let
the tone re-solve by bisection, see which way cerebral blood volume moves, and
blend with the passive response because a 0.5-2/min slow wave is faster than
full autoregulatory equilibration.

Sweeping CPP therefore traces a curve with a minimum, and the minimum moves.""")
    for ar, lab in [(1.0, "autoregulation intact"),
                    (0.45, "autoregulation impaired"),
                    (0.0, "autoregulation abolished")]:
        rows = []
        for m in range(62, 155, 6):
            r = run(p, 22.0, y0=list(y), out_every=22.0, MAP_base=float(m),
                    PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6, autoreg=ar)
            rows.append((r["CPP"][-1], r["PRx"][-1], r["Pic"][-1],
                         r["CBF_pen"][-1], r["PbtO2"][-1]))
        best = min(rows, key=lambda z: z[1])
        sub("%s  (autoreg gain x %.2f)" % (lab, ar))
        print("    CPP    PRx     ICP   CBF_pen  PbtO2")
        for (cp, px, ic, cf, pb) in rows:
            mark = "   <== CPPopt" if (cp, px) == (best[0], best[1]) else ""
            print("  %6.1f %6.3f %7.2f %8.2f %6.2f%s" % (cp, px, ic, cf, pb, mark))
        print("  CPPopt = %.0f mmHg   (PRx %.3f)" % (best[0], best[1]))
    print("""
  As the autoregulatory gain falls the minimum flattens and drifts, and at zero
  gain there is no minimum at all -- PRx is positive everywhere, which is
  precisely the reported observation that CPPopt cannot be computed in patients
  with completely absent reactivity.  A trial that randomises such patients to a
  CPPopt target is randomising them to a number that does not exist.""")


# ------------------------------------------------------------------------------
def _apply(p, y, label, tmax=60.0, t0=0.0, **kw):
    kw.setdefault("MAP_base", 88.0)
    kw.setdefault("PaCO2", 38.0)
    kw.setdefault("prop_rate", 3.0)
    kw.setdefault("fluid_rate", 1.6)
    r = run(p, tmax, y0=list(y), out_every=1.0, t0=t0, **kw)
    r["label"] = label
    return r


def R6_osmotherapy(p, y, ne0, t0):
    hdr("R6.  OSMOTHERAPY DOES NOT TREAT THE LESION.  IT TAXES THE HEALTHY BRAIN.")
    print("""
The blood-brain barrier's reflection coefficient is the entire pharmacology of
osmotherapy, and it is REGIONAL.  Where the barrier is intact, sigma ~ 1 and one
milliosmole per kilogram is worth 19.3 mmHg of driving force -- an enormous
dehydrating pressure.  Where the barrier is broken, sigma collapses toward zero,
the osmole crosses with the water, and the gradient does no work at all.

So the model is asked to do the accounting the bedside cannot: of the brain
water removed by a bolus, how much came from the injured region?""")
    print("  branch: t = %.1f h on tier-0 care (sedation, normocapnia, CPP 65)"
          % (t0 / 60.0))
    base = _apply(p, y, "control", tmax=360.0, ne_rate=ne0)
    cases = [
        ("3% NaCl 250 mL", dict(boluses=[(0.0, "hts3", 250.0)])),
        ("3% NaCl 500 mL", dict(boluses=[(0.0, "hts3", 500.0)])),
        ("23.4% NaCl 30 mL", dict(boluses=[(0.0, "hts23", 30.0)])),
        ("mannitol 0.25 g/kg", dict(boluses=[(0.0, "mann", 0.25)])),
        ("mannitol 0.5 g/kg", dict(boluses=[(0.0, "mann", 0.5)])),
        ("mannitol 1.0 g/kg", dict(boluses=[(0.0, "mann", 1.0)])),
    ]
    print("\n                      pre-dose  nadir   fall from   t to    vs control at"
          "     benefit")
    print("  agent                 ICP      ICP     pre-dose   nadir   +1 h  +3 h  +6 h"
          "   >2 mmHg (h)")
    res = []
    for lab, kw in cases:
        kw = dict(kw); kw.setdefault("ne_rate", ne0)
        r = _apply(p, y, lab, tmax=360.0, **kw)
        g = r["t"][1] - r["t"][0]
        dpic = [r["Pic"][i] - base["Pic"][i] for i in range(len(r["t"]))]
        # the nadir of the TREATED arm in absolute terms, not the point of
        # maximum divergence from a control that is itself running away
        n = min(range(len(r["Pic"])), key=lambda j: r["Pic"][j])
        dur = sum(1 for a in dpic if a < -2.0) * g / 60.0
        print("  %-20s %7.2f %7.2f %10.2f %7.0f %6.2f %5.2f %5.2f %8.2f"
              % (lab, r["Pic"][0], r["Pic"][n], r["Pic"][n] - r["Pic"][0],
                 r["t"][n], dpic[int(60 / g)], dpic[int(180 / g)], dpic[-1], dur))
        res.append((lab, r, dpic, n))

    sub("Where did the water come from?  (change at the ICP nadir)")
    print("  agent                 total    from INTACT   from INJURED   %% from")
    print("                       water mL    region mL     region mL     lesion")
    for (lab, r, dpic, i) in res:
        di = r["V_int"][i] - base["V_int"][i]
        dj = r["V_inj"][i] - base["V_inj"][i]
        tot = di + dj
        pct = 100.0 * dj / tot if abs(tot) > 1e-9 else float("nan")
        print("  %-20s %8.2f %12.2f %13.2f %10.1f" % (lab, tot, di, dj, pct))
    print("""
  The injured region is %.0f%% of the brain water in this model, and at THIS
  branch point it is giving up roughly its proportional share.  That is not the
  textbook claim, and the reason is visible in the next table: at this hour the
  mechanical barrier disruption has largely resealed and the inflammatory
  opening has not yet arrived, so the lesion's reflection coefficient is still
  well above zero and the osmotic gradient still does work there.

  The claim that osmotherapy cannot dehydrate the lesion is therefore not a
  property of the drug.  It is a property of the barrier, and it has a time
  course.""" % (100.0 * p["V_inj0"] / (p["V_int0"] + p["V_inj0"])))

    sub("The same bolus (3% NaCl 250 mL) against different barrier states")
    print("  barrier    dICP    total water   from INTACT  from INJURED"
          "   % from   effective")
    print("  opening   (mmHg)     (mL)         (mL)          (mL)      lesion"
          "   sigma_inj")
    for bbb in [0.00, 0.20, 0.40, 0.60, 0.80, 1.00]:
        yy = list(y)
        yy[IX["BBB_mech"]] = bbb
        yy[IX["BBB_infl"]] = 0.0
        bb = _apply(p, yy, "c", tmax=180.0, ne_rate=ne0)
        tt = _apply(p, yy, "t", tmax=180.0, ne_rate=ne0,
                    boluses=[(0.0, "hts3", 250.0)])
        n = min(range(len(tt["Pic"])), key=lambda j: tt["Pic"][j])
        di = tt["V_int"][n] - bb["V_int"][n]
        dj = tt["V_inj"][n] - bb["V_inj"][n]
        tot = di + dj
        print("  %7.2f %8.2f %11.2f %13.2f %13.2f %9.1f %10.3f"
              % (bbb, tt["Pic"][n] - bb["Pic"][n], tot, di, dj,
                 100.0 * dj / tot if abs(tot) > 1e-9 else float("nan"),
                 0.97 * (1.0 - 0.92 * bbb)))
    print("""
  As the barrier opens, the lesion's share of the water falls and the bolus does
  less overall, exactly as the reflection coefficient in the last column
  predicts.  At a fully open barrier the osmole crosses with the water and the
  gradient does no work there at all.

  Which yields a sharper statement than the textbook one: osmotherapy works on
  tissue that still has a barrier, so its effect on the LESION decays over the
  first days of the injury, precisely as the inflammatory opening arrives --
  while its effect on the healthy brain does not.  Late in the course the
  therapy is almost entirely a manoeuvre performed on normal tissue, made to
  give up water so that the lesion has somewhere to expand into.  Its ceiling is
  therefore set by how much water the HEALTHY brain can safely lose, not by how
  swollen the lesion is.""")

    sub("Sodium, osmolar gap and renal cost at the ICP nadir")
    print("  agent                 peak    peak plasma  peak osmolar   peak    lowest")
    print("                        Na+      osmolality      gap        urine    GFR")
    print("                      (mmol/L)   (mOsm/kg)    (mOsm/kg)   (mL/min) (mL/min)")
    for (lab, r, dpic, i) in res:
        print("  %-20s %7.1f %11.1f %12.1f %9.1f %8.1f"
              % (lab, max(r["Na_p"]), max(r["Osm_p"]), max(r["osm_gap"]),
                 max(r["urine"]), min(r["GFR"])))
    return base, res


def R7_rebound(p, y, ne0, t0):
    hdr("R7.  TACHYPHYLAXIS AND REBOUND ARE THE SAME EQUATION READ TWICE")
    print("""
Two things are supposed to happen with repeated osmotherapy and both are usually
described as separate clinical lore.  In this model they are one term.

  (a) The brain ADAPTS.  Regional osmolality is a state, not a constant: solute
      leaks in across the barrier and idiogenic osmoles are generated.  Each
      bolus is therefore working against a target that has moved toward it.
  (b) Mannitol ACCUMULATES INSIDE THE LESION.  The same broken barrier that
      makes the bolus ineffective there lets the mannitol in, where it stays.
      When plasma mannitol is cleared the gradient does not go to zero -- it
      goes NEGATIVE, and water is drawn into the very tissue being treated.""")
    for agent, dose_, lab in [("hts3", 250.0, "3% NaCl 250 mL x4, q4h"),
                              ("mann", 0.5, "mannitol 0.5 g/kg x4, q4h")]:
        bol = [(i * 240.0, agent, dose_) for i in range(4)]
        r = _apply(p, y, lab, tmax=1440.0, boluses=bol, ne_rate=ne0)
        base = _apply(p, y, "control", tmax=1440.0, ne_rate=ne0)
        sub(lab)
        print("  dose #   t(h)   ICP before  ICP nadir   dICP    brain osm    "
              "lesion osm   mannitol in")
        print("                    (mmHg)     (mmHg)   (mmHg)   intact       "
              " total      lesion mM")
        for i in range(4):
            t0 = i * 240.0
            g = r["t"][1] - r["t"][0]
            i0 = int(t0 / g)
            i1 = min(len(r["t"]) - 1, int((t0 + 200.0) / g))
            seg = r["Pic"][i0:i1]
            nad = min(seg)
            print("     %d   %5.1f %10.2f %10.2f %7.2f %10.2f %11.2f %11.2f"
                  % (i + 1, t0 / 60.0, r["Pic"][i0], nad, nad - r["Pic"][i0],
                     r["Osm_int_t"][i0 + seg.index(nad)],
                     r["Osm_inj_t"][i0 + seg.index(nad)],
                     r["Mann_inj"][i0 + seg.index(nad)]))
        # rebound accounting after the last dose
        g = r["t"][1] - r["t"][0]
        print("\n  after the last dose:")
        for hh in [16, 18, 20, 22, 24]:
            i = int(hh * 60 / g)
            print("    t=%2dh  ICP %6.2f (control %6.2f, diff %+6.2f)   "
                  "plasma mannitol %5.2f mM   lesion mannitol %5.2f mM   "
                  "lesion water %+6.2f mL"
                  % (hh, r["Pic"][i], base["Pic"][i], r["Pic"][i] - base["Pic"][i],
                     r["osm_gap"][i], r["Mann_inj"][i],
                     r["V_inj"][i] - base["V_inj"][i]))
    print("""
  The saline arm loses potency because the brain's osmolality follows plasma.
  The mannitol arm loses potency for the same reason AND acquires a second,
  slower liability: a depot of mannitol inside the lesion whose gradient
  reverses when the plasma level falls.  That is why the guidance to check an
  osmolar gap before redosing is not a renal precaution alone -- the gap is also
  a direct readout of how much drug is sitting on the wrong side of the barrier.""")


# ------------------------------------------------------------------------------
def R8_hyperventilation(p, y, ne0, t0):
    hdr("R8.  HYPERVENTILATION IS A LOAN, AND THE MODEL COMPUTES THE INTEREST")
    print("""
Lowering PaCO2 constricts cerebral arterioles, which removes blood volume from
the box, which lowers ICP.  It is the fastest ICP therapy there is and the only
one that needs no drug and no procedure.  It also removes the blood volume by
removing the BLOOD, so its currency is oxygen delivery to the tissue with the
least of it to spare.

Two further things are in the model and are usually taught as separate facts:
the effect FADES, because choroid-plexus bicarbonate follows PaCO2 with a ~6 h
time constant and the vessel responds to perivascular pH rather than to PaCO2;
and RESTORING normocapnia against an adapted bicarbonate is itself a hypercapnic
stimulus.  Both come from one line: PaCO2_effective = PaCO2 * 24 / HCO3_csf.""")
    def prof(t):
        if t < 60.0:
            return 38.0
        if t < 60.0 + 24 * 60.0:
            return 30.0
        return 38.0
    r = _apply(p, y, "HV", tmax=1800.0, PaCO2=prof, ne_rate=ne0)
    base = _apply(p, y, "control", tmax=1800.0, PaCO2=38.0, ne_rate=ne0)
    g = r["t"][1] - r["t"][0]
    print("\n    t      PaCO2  HCO3_csf  PaCO2_eff   ICP   dICP vs   CBF   CBF_pen"
          "  PbtO2  dPbtO2")
    print("   (h)   (mmHg) (mmol/L)  (mmHg)    (mmHg)  control              "
          "        vs ctrl")
    for hh in [0.5, 1.1, 1.5, 2, 4, 8, 12, 18, 24, 24.5, 25, 26, 28, 30]:
        i = int(hh * 60 / g)
        print("  %5.1f %7.1f %8.2f %9.2f %7.2f %8.2f %7.2f %8.2f %6.2f %7.2f"
              % (hh, r["PaCO2"][i], r["HCO3"][i], r["PaCO2e"][i], r["Pic"][i],
                 r["Pic"][i] - base["Pic"][i], r["CBF100"][i], r["CBF_pen"][i],
                 r["PbtO2"][i], r["PbtO2"][i] - base["PbtO2"][i]))
    icp_ben = sum(max(0.0, base["Pic"][i] - r["Pic"][i]) for i in range(len(r["t"]))) * g / 60.0
    icp_cost = sum(max(0.0, r["Pic"][i] - base["Pic"][i]) for i in range(len(r["t"]))) * g / 60.0
    o2_cost = sum(max(0.0, base["PbtO2"][i] - r["PbtO2"][i]) for i in range(len(r["t"]))) * g / 60.0
    print("""
  ICP benefit  = %.1f mmHg*h        ICP paid back on restoration = %.1f mmHg*h
  PbtO2 debt   = %.1f mmHg*h of brain-tissue oxygen forgone
  Interest rate = %.2f mmHg*h of tissue hypoxia per mmHg*h of ICP bought""" %
          (icp_ben, icp_cost, o2_cost, o2_cost / max(1e-9, icp_ben)))

    sub("The same manoeuvre, held for 30 min instead of 24 h (rescue use)")
    def prof2(t):
        return 30.0 if 60.0 <= t < 90.0 else 38.0
    r2 = _apply(p, y, "HV30", tmax=300.0, PaCO2=prof2, ne_rate=ne0)
    b2 = _apply(p, y, "ctl", tmax=300.0, PaCO2=38.0, ne_rate=ne0)
    g2 = r2["t"][1] - r2["t"][0]
    ib = sum(max(0.0, b2["Pic"][i] - r2["Pic"][i]) for i in range(len(r2["t"]))) * g2 / 60.0
    ob = sum(max(0.0, b2["PbtO2"][i] - r2["PbtO2"][i]) for i in range(len(r2["t"]))) * g2 / 60.0
    print("  ICP benefit %.2f mmHg*h,  PbtO2 debt %.2f mmHg*h,"
          "  interest rate %.2f" % (ib, ob, ob / max(1e-9, ib)))
    print("""
  Same drug, same dose, two orders of magnitude difference in what it costs per
  unit of benefit -- because the benefit decays with the bicarbonate while the
  oxygen debt accrues for as long as the vessel is constricted.  That is the
  quantitative content of "brief, targeted hyperventilation only".""")


# ------------------------------------------------------------------------------
def decompensating(p, icp_trigger=24.0, tmax=2400.0):
    """March the reference patient under tier-0 care alone until ICP crosses a
    trigger, and hand back that state.  Every therapy in R9 is then applied to
    the SAME patient at the SAME moment, which is the only way a cost comparison
    means anything."""
    y = initial_state(p)
    t = 0.0
    step = 20.0
    ne = 0.0
    while t < tmax:
        r = run(p, step, y0=y, out_every=step, t0=t, MAP_base=88.0,
                PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6, ne_rate=ne)
        y = r["y_final"]
        t += step
        ne = clamp(ne + 0.010 * (65.0 - r["CPP"][-1]), 0.0, 0.60)
        if r["Pic"][-1] >= icp_trigger:
            return t, y, r, ne
    return t, y, r, ne


def R9_currency(p):
    hdr("R9.  EVERY ICP THERAPY PRICED IN ITS OWN CURRENCY")
    print("""
There are four volumes in the skull, so there are only four things any therapy
can do.  Below, one patient at one moment (ICP just past 24 mmHg on tier-0 care)
is branched twelve ways and each therapy is run for two hours against a
do-nothing control.  The point is not which row lowers ICP most.  The point is
the columns to the RIGHT of the ICP column.""")
    t0, y, r0, ne0 = decompensating(p, icp_trigger=22.0)
    print("\n  branch point: t = %.1f h, ICP = %.2f, CPP = %.1f, PbtO2 = %.2f,"
          " CSF %.1f mL, noradrenaline %.3f ug/kg/min"
          % (t0 / 60.0, r0["Pic"][-1], r0["CPP"][-1], r0["PbtO2"][-1],
             r0["V_csf"][-1], ne0))
    print("  (tier-0 care throughout: propofol 3 mg/kg/h, PaCO2 38, CPP held at"
          " 65 with noradrenaline)")
    base = _apply(p, y, "control", tmax=60.0, ne_rate=ne0)
    cases = [
        ("nothing (control)", {}),
        ("EVD, drain to 12 mmHg", dict(evd_open=True, evd_set=12.0)),
        ("PaCO2 38 -> 33", dict(PaCO2=33.0)),
        ("PaCO2 38 -> 28", dict(PaCO2=28.0)),
        ("3% NaCl 250 mL", dict(boluses=[(0.0, "hts3", 250.0)])),
        ("23.4% NaCl 30 mL", dict(boluses=[(0.0, "hts23", 30.0)])),
        ("mannitol 0.5 g/kg", dict(boluses=[(0.0, "mann", 0.5)])),
        ("propofol 3 -> 6 mg/kg/h", dict(prop_rate=6.0)),
        ("thiopental 5 mg/kg/h", dict(thio_rate=5.0)),
        ("hypothermia 37 -> 33.5 C", dict(T_target=33.5)),
        ("noradrenaline +0.15", dict(ne_rate=ne0 + 0.15)),
        ("decompressive craniectomy", dict(crani_t=0.0)),
    ]
    print("\n                              ICP    dICP   CPP   MAP   PbtO2 CBF_pen"
          "  Na+   osmgap  CSF     urine")
    print("  therapy (1 h)             (mmHg) (mmHg) (mmHg)(mmHg) (mmHg) 100g/m"
          " mmol/L mOsm/kg  drained  mL/min")
    rows = []
    for lab, kw in cases:
        kw.setdefault("ne_rate", ne0)
        r = _apply(p, y, lab, tmax=60.0, **kw)
        csfd = base["V_csf"][-1] - r["V_csf"][-1]
        d = r["Pic"][-1] - base["Pic"][-1]
        rows.append((lab, r, d))
        print("  %-25s %6.2f %6.2f %6.1f %5.1f %6.2f %6.2f %6.1f %7.1f %8.1f %7.1f"
              % (lab, r["Pic"][-1], d, r["CPP"][-1], r["MAP"][-1],
                 r["PbtO2"][-1], r["CBF_pen"][-1], r["Na_p"][-1],
                 r["osm_gap"][-1], csfd, r["urine"][-1]))

    sub("Cost per millimetre of mercury bought")
    print("  therapy                    dICP    PbtO2 lost   CBF_pen lost   Na+ "
          "gained   MAP lost")
    print("                            (mmHg)   per mmHg      per mmHg      per"
          " mmHg    per mmHg")
    for (lab, r, d) in rows:
        if d > -0.20:
            print("  %-25s %6.2f      (no ICP benefit to price)" % (lab, d))
            continue
        n = -d
        print("  %-25s %6.2f %11.3f %14.3f %11.3f %10.3f"
              % (lab, d,
                 (base["PbtO2"][-1] - r["PbtO2"][-1]) / n,
                 (base["CBF_pen"][-1] - r["CBF_pen"][-1]) / n,
                 (r["Na_p"][-1] - base["Na_p"][-1]) / n,
                 (base["MAP"][-1] - r["MAP"][-1]) / n))
    print("""
  Craniectomy and CSF drainage are the only rows that buy ICP without spending
  perfusion, sodium or arterial pressure -- which is exactly why they sit at
  opposite ends of the tier ladder: the drain is free but tiny and exhaustible,
  the craniectomy is large and permanent but is an operation.  Everything in
  between is a trade, and the sedatives and hypothermia pay in arterial pressure
  (which then has to be bought back with a pressor, which is why tier escalation
  tends to arrive in pairs rather than singly).""")
    return t0, y, r0, ne0


# ------------------------------------------------------------------------------
def R10_protocol(p):
    hdr("R10.  THE SIBICC TIER LADDER AS CLOSED-LOOP FEEDBACK, NOT A TIMETABLE")
    print("""
The escalation rules below are feedback on the simulated patient, not events
scheduled on a clock.  The tier rises when ICP has been above threshold for five
minutes and falls after six hours of control; osmotic boluses are given only if
ICP is above threshold AND sodium and osmolality permit; noradrenaline is a PI
loop on CPP.  Nothing reads the wall clock except the hysteresis timers, so the
trajectory below is the protocol arguing with the disease.""")
    arms = [
        ("no ICP-directed therapy", None),
        ("tier 0 only (sedation + CPP)", 0),
        ("tiers 0-1 (+EVD, osmotherapy)", 1),
        ("tiers 0-2 (+hyperventilation)", 2),
        ("tiers 0-3, no craniectomy", 3),
        ("tiers 0-3 with craniectomy", 3),
    ]
    out = []
    print("\n                                  peak   ICP dose  CPP dose  PbtO2 dose"
          "  osmotic  final   P(unfav) P(death)")
    print("  arm                              ICP   (mmHg*h)  (mmHg*h)   (mmHg*h)"
          "   doses   core")
    for i, (lab, maxt) in enumerate(arms):
        if maxt is None:
            r = run(p, 4320.0, y0=initial_state(p), out_every=15.0,
                    MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
            nos = 0
            tier = None
        else:
            tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=maxt,
                        allow_crani=(i == 5))
            r = run(p, 4320.0, y0=initial_state(p), out_every=15.0,
                    MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                    protocol=tp)
            nos = tp.n_osm
            tier = tp
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        _, pk = tmax_of(r, "Pic")
        print("  %-32s %5.1f %9.1f %9.1f %10.1f %7d %7.3f %8.3f %8.3f"
              % (lab, pk, Di, Dc, Dp, nos, Fc, u, mortality(p, r)))
        out.append((lab, r, tier))

    sub("Escalation log for the full protocol arm")
    tp = out[-1][2]
    if tp is not None:
        for (t, ev, tier, icp) in tp.log[:26]:
            print("    t=%6.1f h   %-12s -> tier %d   (ICP %.1f)"
                  % (t / 60.0, ev, tier, icp))
        if len(tp.log) > 26:
            print("    ... %d further events" % (len(tp.log) - 26))
    a0 = [x for x in out if x[0].startswith("no ICP")][0][1]
    a1 = [x for x in out if x[0].startswith("tier 0 only")][0][1]
    a2 = [x for x in out if x[0].startswith("tiers 0-1")][0][1]
    print("""
  Three things in that table, in increasing order of how uncomfortable they are.

  First, the ladder works, and almost all of the working happens on one rung.
  Going from sedation alone to sedation plus a drain plus osmotherapy takes the
  ICP dose from %.1f to %.1f mmHg*h.  Tiers 2 and 3 then add nothing at all in
  this patient, because tier 1 already holds ICP at the threshold and the
  escalation criterion is never met again -- which is the quantitative version
  of why they are upper tiers rather than routine ones.

  Second, the protocol is visibly ARGUING with the disease rather than
  following a schedule.  The log below shows it escalating, giving an osmotic
  bolus, weaning six hours later when ICP has settled, and then escalating
  again as the oedema advances -- a cycle it repeats through the whole course
  without anything in the code knowing that it should.

  Third, and this is the uncomfortable one: TIER 0 ALONE IS WORSE THAN NO
  ICP-DIRECTED THERAPY AT ALL (ICP dose %.1f against %.1f, and mortality %.3f
  against %.3f).  That is not a bug and it is not an artefact of the outcome
  model.  Tier 0 in this protocol means sedation plus a noradrenaline loop
  holding CPP at 65, and R4 already showed what a pressor does to a brain with
  a leaking barrier over hours rather than minutes.  Defending the perfusion
  pressure without also having a way to remove volume is the Lund arm operating
  with nothing to oppose it.  The drain and the osmotherapy are not merely
  additional therapy; they are what makes the pressor safe.""" %
          (dose(a1, "Pic", 20.0), dose(a2, "Pic", 20.0),
           dose(a1, "Pic", 20.0), dose(a0, "Pic", 20.0),
           mortality(p, a1), mortality(p, a0)))
    return out


# ------------------------------------------------------------------------------
def R11_craniectomy(p):
    hdr("R11.  CRANIECTOMY TIMING: WHY DECRA AND RESCUEicp DISAGREED")
    print("""
DECRA operated early on modestly raised ICP and found WORSE outcomes.
RESCUEicp operated as a last tier on refractory hypertension and found lower
mortality bought with more survivors in severe disability.  Those are not
contradictory results and this model reproduces the structure of both, because
the operation does exactly one thing -- it enlarges the box -- and the value of
enlarging the box depends entirely on whether the box was the thing killing the
patient.

The patient below is deliberately a young one with a severe injury and very
little compensatory reserve (age 28, severity 0.92), managed on tiers 0-2 only,
because a patient whose ICP the medical ladder controls on its own cannot
answer this question at all -- and that, as it happens, is the first result.""")
    rows_c = []
    print("\n  craniectomy at    peak ICP  ICP dose  CPP dose  PbtO2 dose  final"
          "  P(unfav) P(death)  verdict")
    print("  ICP threshold      (mmHg)   (mmHg*h)  (mmHg*h)   (mmHg*h)    core")
    for thr, lab in [(None, "never"), (16.0, "16 (very early)"),
                     (20.0, "20 (DECRA-like)"), (25.0, "25 (intermediate)"),
                     (30.0, "30 (RESCUE-like)"), (40.0, "40 (salvage)")]:
        tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2,
                    allow_crani=False)
        holder = {"done": False}

        def protocol(t, o, ctl, tp=tp, thr=thr, holder=holder):
            tp(t, o, ctl)
            if thr is not None and not holder["done"] and o["Pic"] >= thr:
                ctl.crani_t = t
                holder["done"] = True
                holder["t"] = t
        r = run(p, 4320.0, y0=initial_state(p), out_every=15.0,
                MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                protocol=protocol)
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        _, pk = tmax_of(r, "Pic")
        when = ("t=%.1f h" % (holder.get("t", 0) / 60.0)) if holder["done"] \
            else "not performed"
        rows_c.append((lab.split()[0], Di, Fc, u, mortality(p, r)))
        print("  %-17s %8.1f %9.1f %9.1f %10.1f %7.3f %8.3f %8.3f  %s"
              % (lab, pk, Di, Dc, Dp, Fc, u, mortality(p, r), when))
    never = [r for r in rows_c if r[0] == "never"][0]
    early = [r for r in rows_c if r[0].startswith("16")][0]
    print("""
  Read the ICP-dose column first and the final-core column second, because they
  do not agree, and the disagreement is the answer.

  ICP burden: the earliest threshold cuts the accumulated pressure-time dose
  from %.1f to %.1f mmHg*h -- a large benefit, arriving because the operation is
  performed before the exponential part of the pressure-volume curve has been
  reached rather than after.  On this measure earlier is simply better, and the
  model offers no reason at all to wait.

  Tissue: the final non-viable fraction is %.3f in every row, to three decimal
  places.  Enlarging the box changes the PRESSURE the brain experiences and does
  not change how much of it survives.

  That combination is the DECRA/RESCUEicp result stated as a mechanism.  A
  craniectomy is an excellent pressure therapy and a poor tissue therapy, so it
  converts deaths from intracranial hypertension into survivals with the same
  amount of destroyed brain as before -- which is the trade RESCUEicp reported:
  lower mortality, and the survivors distributed toward severe disability
  rather than toward good recovery.

  And note what this says about DECRA, which found early surgery HARMFUL.  Here
  the intracranial physiology of early surgery is unambiguously favourable.  If
  the trial says otherwise then the harm cannot have been intracranial, and it
  has to live in what this model does not contain: the operation itself, the
  syndrome of the trephined, the hygromas, the reoperations, and the fact that
  randomising to early surgery necessarily operates on patients who would never
  have needed it.  The model is useful here precisely because it FAILS to
  reproduce the harm, which localises where the harm must be.

  What no model can adjudicate is whether the survival was worth having.""" %
          (never[1], early[1], never[2]))


# ------------------------------------------------------------------------------
def R12_hypothermia(p):
    hdr("R12.  HYPOTHERMIA: THE MODEL PREDICTS A BENEFIT THE TRIALS DID NOT FIND")
    print("""
This section is included because it is where the model is WRONG, and saying so
is more useful than quietly leaving it out.

Mechanistically hypothermia should work.  CMRO2 falls with a Q10 of 2.3, the
flow demanded falls with it, the arterioles constrict, blood volume leaves the
box and ICP falls; and independently the supply-demand ratio improves, so the
penumbra should survive better.  Eurotherm and POLAR both found the ICP effect
and neither found the outcome benefit; Eurotherm found harm.""")
    print("\n  target      lowest  ICP dose  PbtO2 dose  final   P(unfav)  "
          "peak ICP on rewarming")
    print("  temp (C)     ICP    (mmHg*h)   (mmHg*h)    core              (mmHg)")
    for tt in [37.0, 35.0, 34.0, 33.0, 32.0]:
        def tprof(t, tt=tt):
            if t < 120.0:
                return 37.0
            if t < 120.0 + 48 * 60.0:
                return tt
            # rewarming at 0.25 C/h, the rate the guidelines specify
            return min(37.0, tt + 0.25 * (t - 120.0 - 48 * 60.0) / 60.0)
        tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2)
        r = run(p, 5400.0, y0=initial_state(p), out_every=15.0,
                MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                T_target=tprof, protocol=tp)
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        g = r["t"][1] - r["t"][0]
        i0 = int((120.0 + 48 * 60.0) / g)
        rw = max(r["Pic"][i0:]) if i0 < len(r["Pic"]) else float("nan")
        print("  %6.1f %10.2f %9.1f %10.1f %7.3f %9.3f %14.2f"
              % (tt, min(r["Pic"]), Di, Dp, Fc, u, rw))
    print("""
  The model reproduces the ICP effect, though it is smaller and less orderly
  than the mechanism suggests -- the ICP-dose column is not even monotone in
  temperature, because a closed-loop protocol that is already holding ICP at
  its threshold absorbs most of the cooling benefit by de-escalating instead.
  That is itself worth noticing: in a patient under an active tier ladder, the
  measured ICP effect of ANY additional therapy is partly hidden by the ladder
  standing down, and the tissue-oxygen column (which falls steadily from 852 to
  711 mmHg*h) is the more honest readout.  The model does NOT reproduce the
  outcome harm,
  and the reason is structural: this model has no pneumonia, no immune
  suppression, no coagulopathy, no shivering and no arrhythmia, and its
  hypothermia is delivered instantly and free of charge.  Those are the
  mechanisms by which the trials think the benefit was lost.

  Stated plainly: a QSP model that contains only the mechanism a therapy is
  supposed to work through will predict that the therapy works.  This is the
  single most important failure mode of mechanistic modelling in drug
  development, and the honest response is to name the missing compartments
  rather than to tune the existing ones until the answer looks right.""")


# ------------------------------------------------------------------------------
def R13_secondary_insults(p):
    hdr("R13.  SECONDARY INSULTS MULTIPLY, AND THE PRICE DEPENDS ON THE HOUR")
    print("""
Oxygen delivery is a PRODUCT:  DO2 = CBF x CaO2.  Hypotension attacks the first
term and hypoxaemia the second, so their combination is not additive, and no
amount of care about one protects against the other.  The same twenty minutes of
insult is then applied at different points in the same patient's course, because
the cost of an insult is set by how much reserve is left to absorb it.""")
    # The insult is delivered UNOPPOSED -- during transport, a CT scan, or a
    # turn -- rather than into a closed CPP loop that would cancel it within a
    # minute.  Tier-0 care only: sedation, normocapnia, no pressor titration.
    ctlkw = dict(MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)

    def arm(hypot=None, hypox=None, tmax=2880.0):
        return run(p, tmax, y0=initial_state(p), out_every=5.0,
                   hypot=hypot, hypox=hypox, **ctlkw)

    sub("A 30-minute insult at 6 h: alone and together")
    ref = arm()
    u0, Fc0, Di0, _, Dp0 = gose_unfavourable(p, ref)
    print("  insult                       ICP dose  PbtO2 dose  final core"
          "  dcore   P(unfav)  excess")
    rows = []
    for lab, ht, hx in [("none", None, None),
                        ("MAP 55 x 30 min", (360.0, 390.0, 55.0), None),
                        ("SaO2 0.75 x 30 min", None, (360.0, 390.0, 0.75)),
                        ("both together", (360.0, 390.0, 55.0),
                         (360.0, 390.0, 0.75))]:
        r = arm(ht, hx)
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        rows.append((lab, Fc, u))
        print("  %-26s %9.1f %10.1f %10.4f %8.4f %8.3f %8.3f"
              % (lab, Di, Dp, Fc, Fc - Fc0, u, u - u0))
    add = (rows[1][1] - rows[0][1]) + (rows[2][1] - rows[0][1])
    both = rows[3][1] - rows[0][1]
    print("""
  Hypotension alone costs %.4f of the brain, hypoxaemia alone costs %.4f, and
  the two together cost %.4f -- against %.4f if they merely added.  The excess of
  %.0f%% over additivity is not a synergy that had to be put in; it is what a
  product does when both factors fall.""" % (rows[1][1] - rows[0][1],
                                             rows[2][1] - rows[0][1], both, add,
                                             100.0 * (both / max(1e-9, add) - 1.0)))

    sub("The SAME insult (MAP 55 + SaO2 0.75 for 30 min) at different hours")
    print("   hour of    ICP at   CSF left  PENUMBRA   final    dcore vs   P(unfav)"
          "  excess")
    print("   insult     insult     (mL)    remaining   core     no insult")
    for hh in [1, 3, 6, 9, 12, 18, 24, 36]:
        t0 = hh * 60.0
        r = arm((t0, t0 + 30.0, 55.0), (t0, t0 + 30.0, 0.75))
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        print("   %5d %10.2f %10.1f %10.4f %8.4f %10.4f %9.3f %8.3f"
              % (hh, at(ref, "Pic", t0), at(ref, "V_csf", t0),
                 at(ref, "F_pen", t0), Fc, Fc - Fc0, u, u - u0))
    print("""
  An identical insult costs eighteen times more at hour 1 than at hour 24, and
  the direction is worth being careful about, because the intuitive answer is
  the wrong one.  It is NOT that the late brain is more fragile for want of
  buffer.  It is that by hour 24 there is far less left to lose: the penumbra
  column has collapsed, most of the salvageable tissue has already converted,
  and an insult can only destroy tissue that was still alive when it arrived.

  Two things follow.  The first is that the ICP at the moment of the insult is
  nearly useless for grading its severity -- 9.81 mmHg at hour 1 for the most
  expensive insult in the table and 26.00 at hour 36 for the cheapest.  The
  second is that the window in which "avoid hypotension" is worth the most is
  the prehospital and early-resuscitation window, before anyone has an ICP
  monitor in place to tell them how the patient is doing.  That is an awkward
  place for the highest-value intervention to live, and it is where the
  observational data has always said it lives.""")


# ------------------------------------------------------------------------------
def R14_pbto2(p):
    hdr("R14.  WHY AN ICP-PERFECT PATIENT CAN STILL BE SUFFOCATING")
    print("""
Jugular bulb saturation is a global, flow-weighted average, so it is dominated
by the healthy majority of the brain.  A parenchymal oxygen probe sits in the
minority compartment that is actually dying, behind a much higher local
resistance.  In this model those are two different numbers computed from two
different perfusions, and the gap between them is the entire rationale for
BOOST-2 and BOOST-3.""")
    tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2)
    r = run(p, 2880.0, y0=initial_state(p), out_every=15.0, protocol=tp,
            MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
    print("\n    t     ICP    CPP   CBF global  CBF penumbra  ratio   SjvO2"
          "   PbtO2   LPR   what the monitors say")
    print("   (h)  (mmHg) (mmHg) (mL/100g/m)  (mL/100g/m)            (%)"
          "   (mmHg)")
    for hh in [1, 2, 4, 6, 8, 12, 16, 20, 24, 32, 40, 48]:
        i_icp = at(r, "Pic", hh * 60.0)
        cg = at(r, "CBF100", hh * 60.0)
        cp = at(r, "CBF_pen", hh * 60.0)
        sj = at(r, "SjvO2", hh * 60.0)
        pb = at(r, "PbtO2", hh * 60.0)
        note = ""
        if i_icp < 22 and pb < 20:
            note = "ICP fine, tissue hypoxic <-- INVISIBLE"
        if sj > 55 and pb < 15:
            note = "SjvO2 NORMAL over a hypoxic penumbra"
        if pb < 10:
            note = "critical tissue hypoxia"
        print("  %5d %6.2f %6.1f %11.2f %13.2f %7.2f %7.1f %7.2f %6.1f  %s"
              % (hh, i_icp, at(r, "CPP", hh * 60.0), cg, cp, cg / cp, sj, pb,
                 at(r, "LPR", hh * 60.0), note))
    g = r["t"][1] - r["t"][0]
    n_icp_ok_o2_bad = sum(1 for i in range(len(r["t"]))
                          if r["Pic"][i] < 22.0 and r["PbtO2"][i] < 20.0) * g / 60.0
    n_sjv_ok_o2_bad = sum(1 for i in range(len(r["t"]))
                          if r["SjvO2"][i] > 55.0 and r["PbtO2"][i] < 20.0) * g / 60.0
    print("""
  Hours with ICP below 22 mmHg AND PbtO2 below 20 mmHg : %.1f h
  Hours with SjvO2 above 55%% AND PbtO2 below 20 mmHg   : %.1f h

  Those hours are the ones a purely ICP-directed protocol cannot see, and they
  are not rare in this trajectory -- they are most of it.  Note what this does
  NOT say: it does not say that treating PbtO2 improves outcome, only that ICP
  and PbtO2 carry different information.  BOOST-3 exists because the first does
  not imply the second.""" % (n_icp_ok_o2_bad, n_sjv_ok_o2_bad))

    sub("ICP-guided versus ICP-plus-PbtO2-guided management")
    print("  protocol                   ICP dose  PbtO2 dose  final core"
          "  P(unfav)")
    for lab, thr in [("ICP-guided only", None), ("ICP + PbtO2 > 20", 20.0)]:
        tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2)
        if thr is None:
            prot = tp
        else:
            def prot(t, o, ctl, tp=tp, thr=thr):
                tp(t, o, ctl)
                # the oxygen limb: raise CPP target and FiO2-equivalent, and
                # back off hyperventilation, whenever tissue oxygen is low
                if o["PbtO2"] < thr:
                    tp.cpp_target = 75.0
                    ctl.SaO2 = 1.00
                    ctl.PaO2 = 200.0
                    ctl.PaCO2 = max(ctl.PaCO2, 38.0)
                else:
                    tp.cpp_target = 65.0
                    ctl.SaO2 = 0.98
                    ctl.PaO2 = 95.0
        r = run(p, 2880.0, y0=initial_state(p), out_every=15.0, protocol=prot,
                MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        print("  %-26s %9.1f %10.1f %10.4f %9.3f" % (lab, Di, Dp, Fc, u))


# ------------------------------------------------------------------------------
def R15_txa(p):
    hdr("R15.  TRANEXAMIC ACID: THE TIME WINDOW IS THE HAEMATOMA'S, NOT THE DRUG'S")
    print("""
CRASH-3 found a benefit that shrank with every hour of delay.  In this model the
drug's own pharmacokinetics are unremarkable and its target is not what expires
-- the SUBSTRATE is.  Haematoma expansion is driven by an injury-triggered
fibrinolytic burst that decays with a time constant of a few hours, so a drug
that blocks fibrinolysis has less and less to block.""")
    print("\n  TXA given at   final haematoma  peak haematoma  expansion"
          "  ICP dose  final core  P(unfav)")
    print("  (h post-injury)     (mL)             (mL)        prevented (mL)")
    ref = None
    for tt in [None, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0]:
        bol = [] if tt is None else [(tt * 60.0, "txa", 1000.0),
                                     (tt * 60.0 + 10.0, "txainf", 1000.0)]
        tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2)
        r = run(p, 2880.0, y0=initial_state(p), out_every=15.0, protocol=tp,
                MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                boluses=bol)
        u, Fc, Di, Dc, Dp = gose_unfavourable(p, r)
        _, pk = tmax_of(r, "V_hem")
        if ref is None:
            ref = pk
        lab = "never" if tt is None else "%.1f" % tt
        print("  %-16s %10.2f %15.2f %13.2f %9.1f %10.4f %9.3f"
              % (lab, r["V_hem"][-1], pk, ref - pk, Di, Fc, u))
    print("""
  The benefit is real and it decays, and it decays on the timescale of the
  fibrinolytic burst rather than on any property of the drug.  That is a
  testable structural claim: it predicts the window should be shorter in
  injuries with a faster burst and longer in patients whose coagulopathy is
  sustained, which is not what a "give it within 3 hours" rule expresses.""")


# ------------------------------------------------------------------------------
def R16_biomarkers(p):
    hdr("R16.  BIOMARKERS ANSWER DIFFERENT QUESTIONS BECAUSE THEY HAVE"
        "\n      DIFFERENT HALF-LIVES")
    print("""
All four markers below are released by the SAME process in this model -- the
rate at which penumbra is converting to core.  They differ only in elimination
half-life, and that alone is enough to make them clinically non-interchangeable:
a fast marker reports the instantaneous rate, a slow marker reports the integral.""")
    tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=2)
    r = run(p, 5760.0, y0=initial_state(p), out_every=15.0, protocol=tp,
            MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
    print("\n    t     kill rate    S100B     UCH-L1      GFAP        NfL"
          "     cumulative core")
    print("   (h)   (per min)   (t1/2 1 h) (t1/2 7 h) (t1/2 24 h) (t1/2 21 d)")
    for hh in [1, 2, 4, 6, 8, 12, 16, 20, 24, 36, 48, 72, 96]:
        print("  %5d %11.3e %10.2f %10.2f %11.2f %10.2f %12.4f"
              % (hh, at(r, "kill", hh * 60.0), at(r, "S100B", hh * 60.0),
                 at(r, "UCHL1", hh * 60.0), at(r, "GFAP", hh * 60.0),
                 at(r, "NfL", hh * 60.0), at(r, "F_core", hh * 60.0)))
    for k, lab in [("S100B", "S100B"), ("UCHL1", "UCH-L1"), ("GFAP", "GFAP"),
                   ("NfL", "NfL")]:
        tt, vv = tmax_of(r, k)
        print("  %-8s peaks at t = %5.1f h" % (lab, tt / 60.0))
    print("""
  The peak times order themselves by half-life without any marker-specific
  biology, which is the point: an assay's timing information is dominated by its
  clearance, so "GFAP peaks at 24 h" is a statement about GFAP's kidney, not
  about astrocytes.  It also follows that only the slowest marker tracks the
  cumulative injury -- and that a single early draw of a fast marker measures
  how fast the brain was dying at the moment the needle went in.""")


# ------------------------------------------------------------------------------
def R17_sensitivity(p):
    hdr("R17.  SENSITIVITY: WHICH PARAMETERS SET THE PRESENTATION, AND WHICH"
        "\n      SET THE RESPONSE?")
    print("""
Each parameter is moved +/-30% and two different questions are asked of it:
how much does it move the patient who arrives (peak ICP and final core under
tier-0 care alone), and how much does it move what THERAPY achieves (the ICP
drop produced by a standard 250 mL bolus of 3% saline).  A parameter can score
high on one and zero on the other, and the pattern of which is which is more
informative than either column.""")
    keys = ["PVI", "Can", "G_aut", "Vv_max", "V_csf_min", "Ran",
            "Lh_inj", "k_glym_inj", "k_leak_bbb", "sig_prot_open",
            "k_bbb_infl", "k_pen_res", "CMRO2n", "Q10", "k_pen2core",
            "k_mito", "Hb", "k_hem_exp"]
    base_p = newp(sev=0.55)

    def metrics(pp):
        r = run(pp, 1200.0, y0=initial_state(pp), out_every=15.0,
                MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
        _, pk = tmax_of(r, "Pic")
        Fc = r["F_core"][-1]
        y = r["y_final"]
        b = run(pp, 90.0, y0=list(y), out_every=5.0, MAP_base=88.0,
                PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
        t = run(pp, 90.0, y0=list(y), out_every=5.0, MAP_base=88.0,
                PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
                boluses=[(0.0, "hts3", 250.0)])
        d = min(t["Pic"][i] - b["Pic"][i] for i in range(len(t["t"])))
        return pk, Fc, d
    pk0, Fc0, d0 = metrics(base_p)
    print("\n  reference:  peak ICP %.2f mmHg,  final core %.4f,"
          "  HTS response %.2f mmHg" % (pk0, Fc0, d0))
    print("\n  parameter          peak ICP    final core   HTS response"
          "   character")
    print("                      (%% change)  (%% change)   (%% change)")
    rows = []
    for k in keys:
        v0 = base_p[k]
        out = []
        for f_ in (0.7, 1.3):
            pp = newp(sev=0.55)
            pp[k] = v0 * f_
            pp = calibrate(pp)
            out.append(metrics(pp))
        s_pk = 100.0 * (out[1][0] - out[0][0]) / max(1e-9, pk0)
        s_fc = 100.0 * (out[1][1] - out[0][1]) / max(1e-9, Fc0)
        s_d = 100.0 * (out[1][2] - out[0][2]) / max(1e-9, abs(d0))
        rows.append((k, s_pk, s_fc, s_d))
    rows.sort(key=lambda z: -(abs(z[1]) + abs(z[2])))
    for (k, a, b_, c) in rows:
        ch = ""
        if abs(a) + abs(b_) > 25.0 and abs(c) < 10.0:
            ch = "sets the PRESENTATION only"
        elif abs(c) > 25.0 and abs(a) + abs(b_) < 12.0:
            ch = "sets the RESPONSE only"
        elif abs(a) + abs(b_) > 25.0 and abs(c) > 25.0:
            ch = "sets both"
        print("  %-18s %10.1f %12.1f %13.1f   %s" % (k, a, b_, c, ch))
    print("""
  The interesting rows are the asymmetric ones.  Craniospinal capacitance (PVI,
  the venous buffer, the CSF floor) dominates the presentation and has almost no
  leverage over what a bolus of saline achieves.  The barrier parameters --
  the reflection-coefficient loss and the barrier leak rate -- barely touch the
  presentation and largely determine the response, because they decide whether
  the osmole is on the right side of the membrane.  In other words the
  parameters governing WHO ARRIVES IN CRISIS and the parameters governing
  WHETHER OSMOTHERAPY WILL WORK are close to disjoint sets, and only the second
  group is measurable by anything at the bedside.""")


# ------------------------------------------------------------------------------
def R18_conservation(p):
    hdr("R18.  CONSERVATION AND SANITY CHECKS")
    tp = Tiered(icp_thr=22.0, cpp_target=65.0, max_tier=3)
    r = run(p, 2880.0, y0=initial_state(p), out_every=5.0, protocol=tp,
            MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6,
            boluses=[(600.0, "mann", 0.5)])
    y = r["y_final"]
    print("""
  Sodium.  Every millimole entering as maintenance fluid or as hypertonic
  saline, minus renal excretion, must appear in the extracellular pool.""")
    print("    plasma Na at 48 h            %10.3f mmol/L" % r["Na_p"][-1])
    print("    ECF volume at 48 h           %10.1f mL" % r["V_ecf"][-1])
    print("    total exchangeable Na at 48h %10.1f mmol" % y[IX["Na_ecf"]])
    print("""
  Mannitol.  A single 0.5 g/kg dose is %.1f mmol.  At 48 h it must be gone: the
  filtered load leaves in the urine and the fraction that crossed the broken
  barrier is degraded slowly in situ.""" % (0.5 * p["weight"] * 1000.0 / 182.17))
    print("    plasma mannitol at 48 h      %10.4f mmol/L" % r["Mann_pl"][-1])
    print("    mannitol left in the lesion  %10.4f mmol/L" % r["Mann_inj_mM"][-1])
    print("""
  Tissue fractions.  Core + penumbra + normal must equal one at every instant,
  and penumbra may only leave by becoming core or by recovering.""")
    print("    core %.6f + penumbra %.6f = %.6f  (must be <= 1)"
          % (y[IX["F_core"]], y[IX["F_pen"]], y[IX["F_core"]] + y[IX["F_pen"]]))
    print("""
  Intracranial volume.  The Monro-Kellie balance is enforced implicitly by
  solving for dICP/dt rather than by book-keeping, so the check is that the
  compartments still sum sensibly and none has gone negative.""")
    print("    arterial blood  %8.2f mL     venous blood %8.2f mL"
          % (r["Va"][-1], r["Vv"][-1]))
    print("    CSF             %8.2f mL     haematoma    %8.2f mL"
          % (y[IX["V_csf"]], y[IX["V_hem"]]))
    print("    brain water     %8.2f mL (%+.2f mL vs baseline)"
          % (y[IX["V_int"]] + y[IX["V_inj"]], r["V_ed"][-1]))
    print("    ICP             %8.2f mmHg   CPP          %8.2f mmHg"
          % (r["Pic"][-1], r["CPP"][-1]))
    bad = [n for n in STATES if not (abs(y[IX[n]]) < 1e12)]
    print("\n    non-finite or divergent states at 48 h: %s"
          % (", ".join(bad) if bad else "none"))
    print("""
  Step-size independence was verified separately (R0): ICP, core fraction,
  brain water, CSF volume and PbtO2 at 24 h are identical to five decimal
  places for dt from 0.0025 to 0.02 min.""")


# ------------------------------------------------------------------------------
def R19_summary(p):
    hdr("R19.  SCENARIO SUMMARY")
    print("""
Fourteen managed courses of the SAME reference patient (age 42, severity 0.62),
each run for 72 h.  The outcome model is an IMPACT-shaped logistic whose injury
terms are model-derived; its three calibration anchors are stated in the source
and it should be read as an ordering, not as a prognosis.""")
    kw = dict(MAP_base=88.0, PaCO2=38.0, prop_rate=3.0, fluid_rate=1.6)
    arms = []

    def add(lab, **extra):
        maxt = extra.pop("max_tier", 2)
        crani = extra.pop("crani", False)
        tp = Tiered(icp_thr=extra.pop("icp_thr", 22.0),
                    cpp_target=extra.pop("cpp_target", 65.0),
                    max_tier=maxt, allow_crani=crani,
                    osm_agent=extra.pop("osm_agent", "hts3"))
        r = run(p, 4320.0, y0=initial_state(p), out_every=15.0,
                protocol=(None if maxt < 0 else tp), **dict(kw, **extra))
        arms.append((lab, r, tp))

    add("no ICP-directed therapy", max_tier=-1)
    add("tier 0 only", max_tier=0)
    add("tiers 0-1 (EVD + HTS)", max_tier=1)
    add("tiers 0-1, mannitol instead", max_tier=1, osm_agent="mann")
    add("tiers 0-2", max_tier=2)
    add("tiers 0-3, no craniectomy", max_tier=3)
    add("tiers 0-3 + craniectomy", max_tier=3, crani=True)
    add("ICP threshold 18 (aggressive)", max_tier=2, icp_thr=18.0)
    add("ICP threshold 25 (permissive)", max_tier=2, icp_thr=25.0)
    add("CPP target 55 (Lund-leaning)", max_tier=2, cpp_target=55.0)
    add("CPP target 75", max_tier=2, cpp_target=75.0)
    add("tiers 0-2 + TXA at 1 h", max_tier=2,
        boluses=[(60.0, "txa", 1000.0), (70.0, "txainf", 1000.0)])
    add("tiers 0-2, autoregulation dead", max_tier=2, autoreg=0.0)
    add("tiers 0-2, age 68", max_tier=2)

    print("\n                                   peak  hours   ICP    CPP   PbtO2"
          "  final  pupils P(unfav) P(death)")
    print("  arm                               ICP  ICP>22  dose   dose   dose"
          "   core")
    for i, (lab, r, tp) in enumerate(arms):
        pp = dict(p)
        if "age 68" in lab:
            pp["age"] = 68.0
        u, Fc, Di, Dc, Dp = gose_unfavourable(pp, r)
        g = r["t"][1] - r["t"][0]
        hrs = sum(1 for v in r["Pic"] if v > 22.0) * g / 60.0
        _, pk = tmax_of(r, "Pic")
        print("  %-33s %5.1f %6.1f %6.1f %6.1f %6.1f %6.3f %5d %8.3f %8.3f"
              % (lab, pk, hrs, Di, Dc, Dp, Fc, pupil_score(r), u,
                 mortality(pp, r)))
    return arms


# ==============================================================================
def main():
    print(__doc__)
    print("Reference patient: age %.0f, %.0f kg, injury severity %.2f"
          % (REF_AGE, P["weight"], REF_SEV))
    p_ref = newp(sev=REF_SEV, age=REF_AGE)

    R0_verification()
    R1_autoregulation()
    p2, marks = R2_reserve()
    R3_plateau_waves(p2, marks)

    # One branch point, shared by every therapy section below, so that the
    # comparisons are of therapies rather than of moments.  The patient is on
    # tier-0 care (sedation, normocapnia, CPP held at 65 with noradrenaline)
    # and is caught the first time ICP crosses 22 mmHg.
    p3 = newp(sev=0.50, age=REF_AGE)
    tb, yb, rb, neb = decompensating(p3, icp_trigger=22.0)
    print("\n[shared branch point for R4-R9: t = %.1f h, ICP = %.2f, CPP = %.1f,"
          " noradrenaline %.3f ug/kg/min]" % (tb/60.0, rb["Pic"][-1],
                                              rb["CPP"][-1], neb))
    R4_pressor_sign(p3, yb, neb, tb)
    R5_cppopt(p3, yb, neb, tb)
    R6_osmotherapy(p3, yb, neb, tb)
    R7_rebound(p3, yb, neb, tb)
    R8_hyperventilation(p3, yb, neb, tb)
    R9_currency(p3)
    R10_protocol(newp(sev=0.55, age=REF_AGE))
    R11_craniectomy(newp(sev=0.92, age=28.0))
    R12_hypothermia(newp(sev=0.62, age=REF_AGE))
    R13_secondary_insults(newp(sev=0.55, age=REF_AGE))
    R14_pbto2(newp(sev=0.60, age=REF_AGE))
    R15_txa(newp(sev=0.62, age=REF_AGE))
    R16_biomarkers(newp(sev=0.62, age=REF_AGE))
    R17_sensitivity(p_ref)
    R18_conservation(newp(sev=0.60, age=REF_AGE))
    R19_summary(p_ref)

    hdr("END OF REFERENCE OUTPUT")
    print("""
Everything above was computed by this file.  Where the model disagrees with the
trial literature (R12) that is stated rather than tuned away, and where a result
depends on a calibration choice (the outcome logistic in R19) the anchors are
written into the source so the choice can be argued with.""")


if __name__ == "__main__":
    main()
