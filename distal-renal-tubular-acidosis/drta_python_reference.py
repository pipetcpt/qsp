#!/usr/bin/env python3
"""
dRTA QSP model -- dependency-free Python reference implementation.

Purpose: this is the *validation harness* for the mrgsolve model.  Every
equation here is transcribed verbatim into drta_mrgsolve_model.R.  Because no R
runtime is available in this environment, the equations are integrated here
(fixed-step RK4) so that the ODE system is actually run, calibrated against
published trial data, and debugged before it is written out as R.

Time unit: HOURS.  Amounts: mmol or mEq.  Concentrations: mmol/L.

STRUCTURAL THESIS
-----------------
(1) Plasma bicarbonate is a RATIO, not a flux.  The renal acid-excretion
    controller is an INTEGRATOR with a SATURATING ACTUATOR; dRTA is actuator
    saturation.  Once railed, the residual acid gap is carried indefinitely by
    the bone sink, so plasma HCO3 settles wherever the bone dose-response says
    it should -- which is why the same patient reads "HCO3 22.0 (near normal)"
    and "lumbar BMD z -1.1 (clearly abnormal)" at the same visit.
(2) Alkali efficiency is RATE MATCHING, not dose.  Base delivered faster than
    endogenous acid is produced pushes HCO3 over the proximal reabsorptive
    threshold and is wasted as bicarbonaturia.  That waste bites hardest
    exactly AT the therapeutic target, so a schedule change moves the
    RESPONDER RATE far more than it moves the MEAN.
(3) The two endpoints want OPPOSITE kinetics.  Systemic alkalinisation wants
    slow delivery (rate matching).  Citraturia wants FAST delivery, because
    NaDC1 is Tm-limited and a bolus escapes reabsorption.  A two-granule
    formulation (fast citrate + slow bicarbonate) is the only schedule that
    optimises both.
"""
import math

# ----------------------------------------------------------------------------
# state vector layout
# ----------------------------------------------------------------------------
SNAMES = [
    # --- gut / drug ---------------------------------------------------------
    "AG_bicIR",   # 0  IR bicarbonate salt in gut               (mEq)
    "AG_citIR",   # 1  IR citrate salt in gut                   (mmol citrate)
    "AG_bicPR",   # 2  ADV7103 slow bicarbonate granules        (mEq)
    "AG_citPR",   # 3  ADV7103 fast citrate granules            (mmol citrate)
    "AG_KCl",     # 4  KCl in gut                               (mmol K)
    "AG_K",       # 5  K carried by the alkali salts, in gut    (mmol K)
    "CIT_pl",     # 6  plasma citrate                           (mmol/L)
    "AG_hctz",    # 7  thiazide in gut                          (mg)
    "C_hctz",     # 8  thiazide plasma                          (mg/L)
    "AG_vitD",    # 9  cholecalciferol in gut                   (IU)
    "C_25D",      # 10 25-OH vitamin D                          (ng/mL)
    # --- acid-base core ----------------------------------------------------
    "HCO3_e",     # 11 ECF bicarbonate                          (mmol/L)
    "BUF",        # 12 non-bicarbonate buffer base donated      (mEq)
    "PaCO2",      # 13 arterial pCO2                            (mmHg)
    "pHi_PT",     # 14 proximal tubule cell pH
    "pHi_IC",     # 15 alpha-intercalated cell pH
    # --- renal actuator / adaptation ---------------------------------------
    "VH",         # 16 V-ATPase controller output (0..1) SATURATING ACTUATOR
    "PEND",       # 17 pendrin abundance (beta-IC HCO3 secretion), relative
    "NH3P",       # 18 ammoniagenic capacity                    (mEq/h)
    "NDC1",       # 19 NaDC1 citrate-reabsorption abundance, relative
    "NEPH",       # 20 functional nephron fraction (0..1)
    "FIB",        # 21 tubulo-interstitial fibrosis index
    # --- electrolytes / volume --------------------------------------------
    "K_pl",       # 22 plasma potassium                         (mmol/L)
    "KDEF",       # 23 total-body K deficit                     (mmol)
    "Cl_pl",      # 24 plasma chloride                          (mmol/L)
    "VECF",       # 25 ECF volume                               (L)
    "ALDO",       # 26 aldosterone, relative to normal
    "Ca_pl",      # 27 plasma ionised calcium                   (mmol/L)
    "PTH",        # 28 PTH                                      (pg/mL)
    "Pi_pl",      # 29 plasma phosphate                         (mmol/L)
    # --- bone --------------------------------------------------------------
    "BLAB",       # 30 rapidly exchangeable bone base pool      (mEq)
    "BMIN",       # 31 bone mineral mass, fraction of expected
    "OC",         # 32 osteoclast activity, relative
    "OB",         # 33 osteoblast activity, relative
    "bALP",       # 34 bone alkaline phosphatase                (U/L)
    "OSM",        # 35 osteomalacia / unmineralised matrix index
    "BMDz",       # 36 lumbar spine BMD z-score
    "CUMBASE",    # 37 cumulative base withdrawn from bone      (mEq)
    # --- urine / stone -----------------------------------------------------
    "UCa_s",      # 38 smoothed urine calcium                   (mmol/day)
    "UCit_s",     # 39 smoothed urine citrate                   (mmol/day)
    "UpH_s",      # 40 smoothed urine pH (24 h mean)
    "UVol_s",     # 41 smoothed urine volume                    (L/day)
    "SS_s",       # 42 smoothed brushite supersaturation index
    "NC",         # 43 nephrocalcinosis burden
    "STONE",      # 44 stone burden                             (arb.)
    # --- systemic / clinical ----------------------------------------------
    "IGF1",       # 45 IGF-1                                    (ng/mL)
    "Hz",         # 46 height z-score
    "MUS",        # 47 muscle strength index (1 = normal)
    "ADH",        # 48 adherence (0..1)
    "GI",         # 49 GI irritation index
    "HEAR",       # 50 hearing threshold shift                  (dB)
    "TBT",        # 51 cumulative time with HCO3 < 22           (h)
    "AAC",        # 52 integral of max(0, 24 - HCO3)            (mmol/L*h)
    "WASTE",      # 53 cumulative alkali wasted as HCO3 in urine (mEq)
    "GIVEN",      # 54 cumulative alkali equivalents absorbed   (mEq)
    "AG_acid",    # 55 NH4Cl / mineral acid load in gut          (mEq)
]
NS = len(SNAMES)
IX = {n: i for i, n in enumerate(SNAMES)}

# ----------------------------------------------------------------------------
# default parameters
# ----------------------------------------------------------------------------
def default_par():
    return dict(
        # ---- subject ------------------------------------------------------
        BW=30.0,           # kg
        AGE=10.0,          # yr
        BSA=1.05,          # m2
        SEX=0,             # 0 = female, 1 = male
        GFR0=110.0,        # mL/min/1.73 m2 at NEPH = 1
        n_intakes=2.0,
        # ---- diet: net endogenous acid production -------------------------
        NEAP_kg=1.35,      # mEq/kg/day  (adult Western ~1.0; children higher)
        f_basal=0.35,      # fraction of NEAP that is non-meal (constant)
        DIET=1.0,          # multiplier: 0.6 low-acid, 1.6 high-protein
        K_diet_kg=1.15,    # mmol K/kg/day dietary potassium
        F_Kdiet=0.90, f_stool_K=0.10,
        # ---- acid-base setpoints -----------------------------------------
        HCO3_set=24.0,     # controller setpoint (age adjusted in setup)
        THR_gap=1.05,      # proximal reabsorptive threshold above setpoint
        CL_bic_ref=2.6,    # L/h per 1.73 m2, clearance of supra-threshold HCO3 [FITTED]
        FE_leak=0.0018,    # fractional HCO3 leak below threshold
        # ---- respiratory compensation -------------------------------------
        # dPaCO2/dHCO3 ~ 1.2 anchored at the NORMAL point (40 mmHg @ 24 mmol/L).
        # Winter's regression (1.5*HCO3+8) is only valid inside the acidotic
        # range; using it at HCO3 24 puts a healthy subject at pCO2 44-47 and
        # therefore at pH 7.36, which made bone donate base in HEALTH.
        PaCO2_0=40.0, kresp=1.20, HCO3_resp=24.0, tau_resp=3.0,
        # ---- cell pH ------------------------------------------------------
        pHi_PT0=7.20, tau_pHiPT=2.0, gPT=0.055,
        pHi_IC0=7.25, tau_pHiIC=2.0, gIC=0.045,
        # ---- V-ATPase controller (INTEGRATOR + SATURATING ACTUATOR) ------
        # The actuator needs BOTH a fast proportional arm (cell pH / pCO2
        # sensing, minutes-hours) and a slow integral arm (V-ATPase trafficking
        # and transcription, days).  BUG FIX #6: with integral control alone at
        # kI = 0.22/h the integrator traversed its whole 0-1 range twice per day
        # and the controller degenerated into a bang-bang oscillator.
        kI_VH=0.011,       # /h per mmol/L  slow integral (trafficking)
        kP_VH=0.125,       # per mmol/L     fast proportional (cell pH sensing)
        VH_max=1.0,        # actuator ceiling (never exceeded)
        # BUG FIX #9: Jh_max_ref was 5.2 mEq/h, capping maximal NAE at ~118
        # mEq/day.  A healthy adult given an NH4Cl load excretes 300-450
        # mEq/day, so the model could not acidify under load and dropped its
        # HCO3 to 12 -- i.e. every healthy control looked like a dRTA patient on
        # the diagnostic test.  Capacity now allows ~4.5x basal NAE.
        Jh_max_ref=11.0,   # mEq/h per 1.73 m2 maximal distal H+ secretion
        dpH_max=2.95,      # maximal blood->urine pH gradient, healthy
        LES=1.0,           # dRTA lesion: fraction of Jh_max retained (rate defect)
        LES_grad=1.0,      # fraction of dpH_max retained (gradient defect)
        # ---- ammoniagenesis ----------------------------------------------
        NH3P_ref=2.20,     # mEq/h per 1.73 m2 at normal acid-base
                # BUG FIX #10: ammoniagenic capacity was driven ONLY by the plasma HCO3
        # error -- but that error is small precisely BECAUSE the kidney is
        # compensating, so the ammonia arm could never be recruited and the
        # V-ATPase controller railed at 1.0 on a merely high-protein diet.
        # Dietary protein supplies the glutamine, so NEAP itself must drive it.
        tau_NH3=72.0, kNH3_ac=2.40, kNH3_K=0.30, kNH3_load=1.05,
        pK_trap=6.62, s_trap=0.62,
        # ---- pendrin (beta-IC HCO3 secretion) ----------------------------
        PEND0=1.0, tau_pend=96.0, kpend=0.45, Jpend_ref=0.22,
        # ---- phosphate / titratable acid ---------------------------------
        Pi_pl0=1.35, TRP0=0.86, kTRP_PTH=0.16,
        pKa_Pi=6.80, TAoth_ref=0.22,  # other non-NH3 buffers, mEq/h, pKa 5.0
        pKa_oth=5.0,
        # ---- urine CO2 ----------------------------------------------------
        pCO2u=48.0,
        # ---- citrate ------------------------------------------------------
        Vd_cit_kg=0.26, CIT_end=0.50, k_ox=1.40, CIT_pl0=0.12,
        FPE_cit=0.85,      # hepatic first-pass extraction of oral citrate
        Tm_cit_ref=0.895,  # mmol/h per 1.73 m2, NaDC1 Tm at NDC1=1 [FITTED]
        Km_cit_ref=0.28,   # mmol/h per 1.73 m2, half-saturating filtered load
        tau_ndc1=48.0, kndc1=1.9,   # NaDC1 up-regulation by PT cell acidosis
        FEcit_min=0.040,   # floor on fractional excretion of citrate
        pKa3_cit=5.60, K_cacit=0.28,  # L/mmol Ca-citrate complexation
        # ---- drug absorption / release ------------------------------------
        ka_bicIR=1.50, F_bicIR=0.86,
        ka_citIR=1.05, F_citIR=0.94,
        kr_bicPR=0.215,    # ADV7103 bicarbonate granules: ~10-12 h release
        kr_citPR=0.90,     # ADV7103 citrate granules:      ~2-3 h release
        ka_K=1.20, F_K=0.90, f_Kcat=1.0,   # 1.0 = K-based alkali, 0.0 = Na-based
        f_cit_ADV=0.35,    # fraction of ADV7103 alkali equivalents as citrate
        # ---- potassium ----------------------------------------------------
        # plasma K is itself a RATIO: it reads total-body deficit only through a
        # ~300 mmol-per-mmol/L exchange coefficient, so it under-reports stores
        # exactly as plasma HCO3 under-reports the acid gap.
        K_pl0=4.2, kKdef=4.3, tau_K=6.0, kK_sec=0.0, kK_H=0.55,
        kK_Na=0.030,   # extra K+ secretion per mEq/h of Na-based alkali
        kK_aldo=1.20, kK_shift=1.6, kK_hctz=0.45, nK_sec=2.0,
        # ---- volume / aldosterone ----------------------------------------
        tau_aldo=12.0, kaldo_V=3.2, kaldo_K=6.0,
        # ---- calcium / PTH / bone ----------------------------------------
        Ca_pl0=1.22, tau_Ca=8.0, PTH0=32.0, tau_PTH=2.0, kPTH_Ca=48.0,
        kpc_bone=9.5,      # physicochemical dissolution gain  [FITTED]
        kcell_bone=6.0,    # osteoclast-mediated gain          [FITTED]
        pH_ref=7.395,
        BLAB_ref=1400.0, tau_blab=240.0,
        BaseCap_ref=52000.0,  # mEq base per unit BMIN (whole skeleton, 70 kg)
        k_remod=5.7e-6,    # basal remodelling, fraction of BaseCap per hour
        rCaBase=0.035,     # mmol Ca released per mEq bone base [FITTED]
        # ---- non-bicarbonate buffer (apparent HCO3 space) -----------------
        # normal apparent space ~0.5 L/kg = ECF 0.20 + buffer 0.30; the buffer
        # arm EXPANDS as HCO3 falls, reproducing the classical 0.7-1.0 L/kg
        # apparent space of severe metabolic acidosis.
        Cbuf_kg=0.30, kbuf=0.25, Cbuf_exp=0.90, Cbuf_hco3=20.0,
        tau_OC=96.0, kOC_ac=14.0, kOC_PTH=0.55,
        tau_OB=168.0, kOB_ac=6.0,
        bALP0=90.0, tau_bALP=120.0, kbALP=55.0,
        tau_OSM=720.0, kOSM=1.3, kOSM_D=0.9,
        kz_catch=6.0e-6,   # /h  BMD z catch-up                [FITTED]
        kz_loss=1.05e-3,   # BMD z loss per (mEq/h) bone base flux [FITTED]
        # ---- urine calcium ------------------------------------------------
        # NOTE: the Lemann 1999 slope dUCa/dNAE = 0.035 mmol/mEq (PMID 9873210)
        # is deliberately NOT a parameter here -- it is a VALIDATION TARGET.
        # Structurally, calciuria comes from bone Ca release plus acid
        # inhibition of distal TRPV5-mediated Ca reabsorption.
        UCa0=0.055,        # mmol/kg/day baseline
        kUCa_bone=1.0,     # mmol urine Ca per mmol Ca released from bone
        kUCa_NEAP=0.035,   # mmol Ca per mEq acid load (Lemann 1999, PMID 9873210)
        kUCa_pH=50.0,      # mmol/day per pH unit of acidaemia per 70 kg [FITTED]
        kUCa_hctz=0.42, tau_UCa=12.0,
        # ---- urine volume -------------------------------------------------
        UVol_kg=0.035,     # L/kg/day
        kDI=0.55, tau_UVol=12.0,
        # ---- supersaturation / nephrocalcinosis ---------------------------
        Ksp_ref=4.25, tau_SS=12.0,  # normalised: SS = 1 is the brushite threshold
        kNC=1.05e-4, nNC=1.35, kNC_res=1.1e-5,
        kSTONE=7.0e-5, kSTONE_res=2.0e-5,
        kFIB=4.0e-4, tau_FIB=2000.0,
        kNEPH=1.1e-5, NEPH_min=0.12,
        # ---- growth -------------------------------------------------------
        IGF1_ref=220.0, tau_IGF=168.0, K_ac_igf=3.2, h_igf=1.7,
        kh_catch=1.60e-5,  # /h  height z catch-up             [FITTED]
        kh_drag=7.0e-5,    # /h  height z loss when acidotic   [FITTED]
        # ---- muscle / adherence / GI --------------------------------------
        tau_MUS=24.0, kMUS=0.55,
        ADH_ref=0.92, tau_ADH=720.0, kADH_n=0.130, kADH_GI=0.22,
        tau_GI=8.0, kGI=0.042,     # GI irritation per (mEq/h)/kg of IR bolus
        # ---- hearing ------------------------------------------------------
        kHEAR=0.0, HEAR_max=60.0,
        # ---- thiazide / vitamin D ----------------------------------------
        ka_hctz=1.1, ke_hctz=0.088, Vd_hctz=3.6, EC50_hctz=0.09,
        ka_vitD=0.05, ke_25D=0.00035, F_vitD=6.5e-5,
    )


def setup(par):
    """
    Derive age/size-dependent parameters.  MUST BE IDEMPOTENT.

    BUG FIX #7: this function used to overwrite parameters with hardcoded
    constants (e.g. `par["Tm_cit"] = 0.62 * BSA/1.73`), which silently discarded
    every caller override, AND it is invoked twice per simulation, so any
    in-place multiplicative scaling was applied twice.  All BSA/BW scaling now
    reads from separate `*_ref` base constants so re-running setup() is a no-op.
    """
    a = par["AGE"]
    fb = par["BSA"] / 1.73
    fw = par["BW"] / 70.0
    # normal plasma HCO3 setpoint rises with age (infant ~21.5 -> adult ~24.7)
    par["HCO3_set"] = 21.0 + 4.0 * (a / (a + 3.0))
    par["THR_bic"] = par["HCO3_set"] + par["THR_gap"]
    # net endogenous acid production per kg falls with age
    par["NEAP_kg"] = 0.95 + 1.45 * math.exp(-a / 6.5)
    # capacities scale with kidney size (BSA)
    par["Jh_max"] = par["Jh_max_ref"] * fb
    par["NH3P0"] = par["NH3P_ref"] * fb
    par["Tm_cit"] = par["Tm_cit_ref"] * fb
    par["Km_cit"] = par["Km_cit_ref"] * fb
    par["TAoth"] = par["TAoth_ref"] * fb
    par["Jpend"] = par["Jpend_ref"] * fb
    par["GFRh"] = par["GFR0"] * 60.0 / 1000.0 * fb          # L/h
    par["CL_bic"] = par["CL_bic_ref"] * fb
    par["Vd_cit"] = par["Vd_cit_kg"] * par["BW"]
    par["VECF0"] = 0.20 * par["BW"]
    par["BLAB0"] = par["BLAB_ref"] * fw
    par["BaseCap"] = par["BaseCap_ref"] * fw
    par["growth_pot"] = 1.0 / (1.0 + math.exp((a - 15.5) / 1.8))  # 1 child, 0 adult
    par["NEAP_h"] = par["NEAP_kg"] * par["BW"] * par["DIET"] / 24.0
    par["NEAP_ref_h"] = par["NEAP_kg"] * par["BW"] / 24.0   # DIET = 1 reference

    # --- self-consistent baselines -------------------------------------------
    # (a) the bone-dissolution threshold pH_ref must be the pH THIS SUBJECT
    #     reaches in health, or a normal toddler slowly dissolves its skeleton
    pco2_h = min(62.0, max(14.0, par["PaCO2_0"]
                           + par["kresp"] * (par["HCO3_set"] - par["HCO3_resp"])))
    par["pH_ref"] = blood_pH(par["HCO3_set"], pco2_h)
    # (b) endogenous citrate appearance must hold plasma citrate at CIT_pl0
    par["CIT_end"] = (par["k_ox"] * par["Vd_cit"] + par["GFRh"]) * par["CIT_pl0"]
    # (c) basal renal K+ secretion must exactly balance dietary K+ intake, so
    #     any hypokalaemia comes from the LESION and not from tuning
    K_diet_h = par["K_diet_kg"] * par["BW"] / 24.0
    par["kK_sec"] = ((par["F_Kdiet"] - par["f_stool_K"]) * K_diet_h) \
        / (par["K_pl0"] * par["BW"] / 24.0)
    return par


# ----------------------------------------------------------------------------
# dosing schedule
# ----------------------------------------------------------------------------
class Regimen:
    """Alkali regimen.  `kind` in {none, bicIR, citIR, ADV, mixed}."""
    def __init__(self, kind="none", mEq_kg_day=0.0, times=(7.0, 13.0, 19.0),
                 KCl_mmol_day=0.0, hctz_mg=0.0, vitD_IU=0.0, f_cit=None,
                 cation="K", acid_load=(), start_h=0.0, stop_h=1e18):
        self.kind = kind
        self.mEq_kg_day = mEq_kg_day
        self.times = tuple(sorted(times))
        self.KCl = KCl_mmol_day
        self.hctz = hctz_mg
        self.vitD = vitD_IU
        self.f_cit = f_cit
        self.cation = cation
        self.acid_load = tuple(acid_load)   # ((time_h, mEq), ...)
        self.start_h = start_h
        self.stop_h = stop_h

    def n_intakes(self):
        return max(1, len(self.times))


def dose_events(reg, par, t_end):
    """Return sorted list of (time_h, state_index, amount)."""
    ev = []
    if reg.kind == "none" or reg.mEq_kg_day <= 0:
        per_dose = 0.0
    else:
        per_dose = reg.mEq_kg_day * par["BW"] / len(reg.times)
    fc = reg.f_cit
    if fc is None:
        fc = {"bicIR": 0.0, "citIR": 1.0, "ADV": par["f_cit_ADV"],
              "mixed": 0.5, "none": 0.0}[reg.kind]
    day0 = int(reg.start_h // 24)
    nday = int(t_end // 24) + 2
    for d in range(day0, nday):
        for ct in reg.times:
            t = d * 24.0 + ct
            if t < reg.start_h or t > reg.stop_h or t > t_end:
                continue
            if per_dose > 0:
                mEq_cit = per_dose * fc
                mEq_bic = per_dose * (1.0 - fc)
                mmol_cit = mEq_cit / 3.0        # citrate3- -> 3 HCO3
                if reg.kind == "ADV":
                    if mmol_cit > 0:
                        ev.append((t, IX["AG_citPR"], mmol_cit))
                    if mEq_bic > 0:
                        ev.append((t, IX["AG_bicPR"], mEq_bic))
                else:
                    if mmol_cit > 0:
                        ev.append((t, IX["AG_citIR"], mmol_cit))
                    if mEq_bic > 0:
                        ev.append((t, IX["AG_bicIR"], mEq_bic))
                # potassium carried by the salt (K-citrate / KHCO3)
                ev.append((t, IX["AG_K"], per_dose))
            if reg.KCl > 0:
                ev.append((t, IX["AG_KCl"], reg.KCl / len(reg.times)))
            if reg.hctz > 0 and abs(ct - reg.times[0]) < 1e-9:
                ev.append((t, IX["AG_hctz"], reg.hctz))
            if reg.vitD > 0 and abs(ct - reg.times[0]) < 1e-9:
                ev.append((t, IX["AG_vitD"], reg.vitD))
    for (tt, amt) in getattr(reg, "acid_load", ()):
        if tt <= t_end:
            ev.append((tt, IX["AG_acid"], amt))
    ev.sort()
    return ev


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def neap_rate(t, par):
    """Net endogenous acid production, mEq/h, with meal-linked diurnal shape."""
    tot = par["NEAP_kg"] * par["BW"] * par["DIET"] / 24.0
    basal = tot * par["f_basal"]
    meal = tot * (1.0 - par["f_basal"])
    h = t % 24.0
    s = 0.0
    for mt, frac in ((7.5, 0.28), (12.8, 0.34), (19.2, 0.38)):
        d = h - mt
        if d < -12:
            d += 24
        if d > 12:
            d -= 24
        if -1.5 < d < 3.0:                     # raised-cosine over 4.5 h
            s += frac * (1.0 - math.cos(2 * math.pi * (d + 1.5) / 4.5)) / 4.5 * 2.0
    # BUG FIX #8 (found by simulation): the raised-cosine shape integrates to
    # 2.0 h over a day, not 24 h, so the meal-linked arm delivered only 1/12 of
    # its intended acid load and TOTAL NEAP was ~35% of the prescribed value.
    # Every subject therefore looked far less acidotic than the diet implied.
    return basal + meal * s * 12.0


def blood_pH(hco3, paco2):
    hco3 = max(hco3, 0.5)
    paco2 = max(paco2, 8.0)
    return 6.10 + math.log10(hco3 / (0.0301 * paco2))


def _frac_prot(pH, pKa):
    return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def solve_urine_pH(Jh_cap_max, pH_floor, HCO3_del, Pi_load, NH3_avail,
                   Cit_load, TAoth, pKa_Pi, pKa_oth, pKa3_cit, pCO2u,
                   UVol_h, pK_trap, s_trap, pH_bl, guess):
    """
    Solve the luminal proton balance for urine pH.

    F(pH) = Jh_cap(pH) - Demand(pH),  monotonically increasing in pH.
      Jh_cap(pH) = Jh_max * (1 - 10^-(pH - pH_floor))      thermodynamic ceiling
      Demand(pH) = (HCO3 delivered - HCO3 excreted)        base titrated
                 + titratable acid formed
                 + citrate protonated
                 + NH3 trapped as NH4+
    Warm-started bracketed bisection (12 iterations, +/-0.0005 pH).
    """
    lo, hi = max(pH_floor + 1e-4, 4.20), 8.40
    # warm start bracket
    if guess is not None:
        lo = max(lo, guess - 0.45)
        hi = min(hi, guess + 0.45)
        if lo >= hi:
            lo, hi = max(pH_floor + 1e-4, 4.20), 8.40

    def F(pH):
        jh = Jh_cap_max * (1.0 - 10.0 ** (-(pH - pH_floor))) if pH > pH_floor else 0.0
        hco3u = UVol_h * 0.0301 * pCO2u * 10.0 ** (pH - 6.10)
        ta = Pi_load * (_frac_prot(pH, pKa_Pi) - _frac_prot(pH_bl, pKa_Pi))
        ta += TAoth * (_frac_prot(pH, pKa_oth) - _frac_prot(pH_bl, pKa_oth))
        cith = Cit_load * (_frac_prot(pH, pKa3_cit) - _frac_prot(pH_bl, pKa3_cit))
        trap = 1.0 / (1.0 + 10.0 ** ((pH - pK_trap) / s_trap))
        nh4 = NH3_avail * trap
        return jh - ((HCO3_del - hco3u) + ta + cith + nh4)

    flo, fhi = F(lo), F(hi)
    if flo > 0.0:
        lo2 = max(pH_floor + 1e-4, 4.20)
        if lo2 < lo and F(lo2) <= 0.0:
            lo, flo = lo2, F(lo2)
        else:
            return lo
    if fhi < 0.0:
        hi2 = 8.40
        if hi2 > hi and F(hi2) >= 0.0:
            hi, fhi = hi2, F(hi2)
        else:
            return hi
    for _ in range(12):
        mid = 0.5 * (lo + hi)
        if F(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------------------
# right-hand side
# ----------------------------------------------------------------------------
class Model:
    def __init__(self, par):
        self.p = setup(dict(par))
        self._uph = 6.2          # warm start for the urine-pH solve
        self.last = {}

    # -- initial condition ---------------------------------------------------
    def init(self):
        p = self.p
        y = [0.0] * NS
        y[IX["CIT_pl"]] = p["CIT_pl0"]
        y[IX["HCO3_e"]] = p["HCO3_set"]
        y[IX["PaCO2"]] = p["PaCO2_0"]
        y[IX["pHi_PT"]] = p["pHi_PT0"]
        y[IX["pHi_IC"]] = p["pHi_IC0"]
        y[IX["VH"]] = 0.60
        y[IX["PEND"]] = p["PEND0"]
        y[IX["NH3P"]] = p["NH3P0"]
        y[IX["NDC1"]] = 1.0
        y[IX["NEPH"]] = 1.0
        y[IX["K_pl"]] = p["K_pl0"]
        y[IX["Cl_pl"]] = 104.0
        y[IX["VECF"]] = p["VECF0"]
        y[IX["ALDO"]] = 1.0
        y[IX["Ca_pl"]] = p["Ca_pl0"]
        y[IX["PTH"]] = p["PTH0"]
        y[IX["Pi_pl"]] = p["Pi_pl0"]
        y[IX["BLAB"]] = p["BLAB0"]
        y[IX["BMIN"]] = 1.0
        y[IX["OC"]] = 1.0
        y[IX["OB"]] = 1.0
        y[IX["bALP"]] = p["bALP0"]
        y[IX["BMDz"]] = 0.0
        y[IX["UCa_s"]] = p["UCa0"] * p["BW"]
        y[IX["UCit_s"]] = 3.0
        y[IX["UpH_s"]] = 6.0
        y[IX["UVol_s"]] = p["UVol_kg"] * p["BW"]
        y[IX["SS_s"]] = 1.0
        y[IX["IGF1"]] = p["IGF1_ref"]
        y[IX["MUS"]] = 1.0
        y[IX["ADH"]] = p["ADH_ref"]
        y[IX["C_25D"]] = 28.0
        return y

    # -- rhs -----------------------------------------------------------------
    def rhs(self, t, y, store=False):
        p = self.p
        g = y.__getitem__
        d = [0.0] * NS

        HCO3 = max(g(IX["HCO3_e"]), 1.0)
        PaCO2 = max(g(IX["PaCO2"]), 8.0)
        pH_bl = blood_pH(HCO3, PaCO2)
        NEPH = max(g(IX["NEPH"]), p["NEPH_min"])
        GFRh = p["GFRh"] * NEPH
        VECF = max(g(IX["VECF"]), 0.3 * p["VECF0"])

        # ---------------- gut / drug absorption ----------------------------
        adh = min(1.0, max(0.0, g(IX["ADH"])))
        a_bicIR = p["ka_bicIR"] * g(IX["AG_bicIR"])
        a_citIR = p["ka_citIR"] * g(IX["AG_citIR"])
        r_bicPR = p["kr_bicPR"] * g(IX["AG_bicPR"])
        r_citPR = p["kr_citPR"] * g(IX["AG_citPR"])
        d[IX["AG_bicIR"]] = -a_bicIR
        d[IX["AG_citIR"]] = -a_citIR
        d[IX["AG_bicPR"]] = -r_bicPR
        d[IX["AG_citPR"]] = -r_citPR

        bic_abs = p["F_bicIR"] * a_bicIR + p["F_bicIR"] * r_bicPR     # mEq/h
        cit_abs = p["F_citIR"] * (a_citIR + r_citPR)                  # mmol/h
        # NH4Cl / mineral-acid loading test (the diagnostic challenge)
        a_acid = p["ka_bicIR"] * g(IX["AG_acid"])
        d[IX["AG_acid"]] = -a_acid
        acid_abs = 0.95 * a_acid                                      # mEq/h

        # AG_K holds the CATION equivalents carried by the alkali salt.  f_Kcat
        # splits them into K+ and Na+: NaHCO3 delivers Na+ to the distal nephron
        # and therefore AGGRAVATES the K+ wasting it is meant to treat, whereas
        # K-citrate / KHCO3 repletes it.  This is the whole reason dRTA alkali
        # is potassium-based.
        a_salt = p["ka_K"] * g(IX["AG_K"])
        a_kcl = p["ka_K"] * g(IX["AG_KCl"])
        d[IX["AG_K"]] = -a_salt
        d[IX["AG_KCl"]] = -a_kcl
        K_abs_K = p["F_K"] * (p["f_Kcat"] * a_salt + a_kcl)
        Na_alk_h = p["F_K"] * (1.0 - p["f_Kcat"]) * a_salt

        d[IX["AG_hctz"]] = -p["ka_hctz"] * g(IX["AG_hctz"])
        d[IX["C_hctz"]] = (p["ka_hctz"] * g(IX["AG_hctz"]) / p["Vd_hctz"]
                           - p["ke_hctz"] * g(IX["C_hctz"]))
        d[IX["AG_vitD"]] = -p["ka_vitD"] * g(IX["AG_vitD"])
        d[IX["C_25D"]] = (p["F_vitD"] * p["ka_vitD"] * g(IX["AG_vitD"])
                          - p["ke_25D"] * g(IX["C_25D"]))

        # ---------------- citrate disposition ------------------------------
        CIT = max(g(IX["CIT_pl"]), 1e-6)
        cit_filt = GFRh * CIT                                   # mmol/h
        cit_ox = p["k_ox"] * CIT * p["Vd_cit"]                  # mmol/h
        # oral citrate is largely extracted on first pass through the liver and
        # oxidised there; only the escaping fraction reaches the systemic
        # circulation and can be FILTERED.
        cit_fp = p["FPE_cit"] * cit_abs
        cit_sys = (1.0 - p["FPE_cit"]) * cit_abs
        d[IX["CIT_pl"]] = (cit_sys + p["CIT_end"] - cit_ox - cit_filt) / p["Vd_cit"]
        # Every mmol of citrate3- oxidised consumes 3 H+ == yields 3 mEq HCO3,
        # but ONLY drug-derived citrate is an alkali source; endogenous citrate
        # turnover is acid-base neutral.  BUG FIX #5: subtracting the whole
        # endogenous appearance rate CIT_end (rather than the baseline oxidation
        # rate) fed a phantom -57 mEq/day of acid into every healthy subject,
        # which drove the V-ATPase controller onto its ceiling at baseline.
        cit_ox_drug = max(0.0, p["k_ox"] * p["Vd_cit"] * (CIT - p["CIT_pl0"]))
        alk_from_cit = 3.0 * (cit_fp + cit_ox_drug)

        # NaDC1: Tm-limited reabsorption -> a BOLUS escapes reabsorption
        NDC1 = max(g(IX["NDC1"]), 0.05)
        Tm = p["Tm_cit"] * NDC1 * NEPH
        cit_reab = Tm * cit_filt / (p["Km_cit"] + cit_filt)
        # NaDC1 cannot reabsorb the whole filtered load: fractional excretion of
        # citrate has a floor of a few percent even in severe acidosis, so
        # hypocitraturia is profound but never absolute.
        UCit_h = max(p["FEcit_min"] * cit_filt, cit_filt - cit_reab)   # mmol/h

        # ---------------- respiratory compensation -------------------------
        PaCO2_t = min(62.0, max(14.0,
                     p["PaCO2_0"] + p["kresp"] * (HCO3 - p["HCO3_resp"])))
        d[IX["PaCO2"]] = (PaCO2_t - PaCO2) / p["tau_resp"]

        # ---------------- cell pH ------------------------------------------
        pHiPT_t = p["pHi_PT0"] - p["gPT"] * (p["HCO3_set"] - HCO3)
        d[IX["pHi_PT"]] = (pHiPT_t - g(IX["pHi_PT"])) / p["tau_pHiPT"]
        pHiIC_t = p["pHi_IC0"] - p["gIC"] * (p["HCO3_set"] - HCO3)
        d[IX["pHi_IC"]] = (pHiIC_t - g(IX["pHi_IC"])) / p["tau_pHiIC"]

        # ---------------- V-ATPase: INTEGRAL CONTROL, SATURATING ACTUATOR --
        err = p["HCO3_set"] - HCO3
        VH = g(IX["VH"])
        VHraw = VH + p["kP_VH"] * err
        VHc = min(max(VHraw, 0.0), p["VH_max"])       # <-- SATURATING ACTUATOR
        dVH = p["kI_VH"] * err
        # anti-windup: stop integrating once the EFFECTIVE actuator is railed
        if (VHraw >= p["VH_max"] and dVH > 0.0) or (VHraw <= 0.0 and dVH < 0.0):
            dVH = 0.0
        d[IX["VH"]] = dVH

        # ---------------- distal H+ secretion / urine pH -------------------
        UVol = max(g(IX["UVol_s"]), 0.2)
        UVol_h = UVol / 24.0
        # bicarbonate delivered distally
        HCO3_del = (p["CL_bic"] * NEPH * max(0.0, HCO3 - p["THR_bic"])
                    + p["FE_leak"] * GFRh * HCO3
                    + p["Jpend"] * g(IX["PEND"]) * (p["BSA"] / 1.73))
        # phosphate delivered
        TRP = min(0.95, max(0.35, p["TRP0"] - p["kTRP_PTH"] * (g(IX["PTH"]) / p["PTH0"] - 1.0)))
        Pi_load = GFRh * max(g(IX["Pi_pl"]), 0.2) * (1.0 - TRP)
        NH3_avail = max(0.0, g(IX["NH3P"])) * NEPH
        Jh_cap_max = p["Jh_max"] * p["LES"] * VHc * NEPH
        pH_floor = pH_bl - p["dpH_max"] * p["LES_grad"]

        uph = solve_urine_pH(Jh_cap_max, pH_floor, HCO3_del, Pi_load, NH3_avail,
                             UCit_h, p["TAoth"], p["pKa_Pi"], p["pKa_oth"],
                             p["pKa3_cit"], p["pCO2u"], UVol_h,
                             p["pK_trap"], p["s_trap"], pH_bl, self._uph)
        self._uph = uph

        HCO3_u = UVol_h * 0.0301 * p["pCO2u"] * 10.0 ** (uph - 6.10)      # mEq/h
        TA = (Pi_load * (_frac_prot(uph, p["pKa_Pi"]) - _frac_prot(pH_bl, p["pKa_Pi"]))
              + p["TAoth"] * (_frac_prot(uph, p["pKa_oth"]) - _frac_prot(pH_bl, p["pKa_oth"])))
        trap = 1.0 / (1.0 + 10.0 ** ((uph - p["pK_trap"]) / p["s_trap"]))
        NH4_u = NH3_avail * trap
        NAE = TA + NH4_u - HCO3_u                                          # mEq/h

        # ---------------- whole-body acid balance --------------------------
        NEAP = neap_rate(t, p) + acid_abs
        alk_in = bic_abs + alk_from_cit
        # --- SINK 2: non-bicarbonate (intracellular / protein) buffer -------
        # BUG FIX #1 (found by simulation): this flux was subtracted from
        # d(HCO3)/dt instead of added.  With the wrong sign the buffer became a
        # POSITIVE feedback and HCO3_e diverged to 5e9 mmol/L within 14 h.
        # Jbuf > 0 means the buffer DONATES base to the ECF, so it enters
        # d(HCO3)/dt with a plus sign and d(BUF)/dt with a plus sign too.
        Cbuf = p["Cbuf_kg"] * p["BW"] * (
            1.0 + p["Cbuf_exp"] * max(0.0, p["Cbuf_hco3"] - HCO3) / 10.0)
        Jbuf = p["kbuf"] * ((p["HCO3_set"] - HCO3) * Cbuf - g(IX["BUF"]))
        d[IX["BUF"]] = Jbuf
        # --- SINK 3: bone.  RECTIFIED in pH -> convex, so a spiky HCO3
        # profile costs more bone than a flat one of the same MEAN (Jensen).
        acid_pH = max(0.0, p["pH_ref"] - pH_bl)
        blab_av = max(0.0, g(IX["BLAB"])) / p["BLAB0"]
        Jbone = (p["kpc_bone"] * max(0.2, g(IX["BMIN"])) * blab_av
                 + p["kcell_bone"] * max(0.0, g(IX["OC"]))) * acid_pH * (p["BW"] / 70.0)
        d[IX["BLAB"]] = (p["BLAB0"] - g(IX["BLAB"])) / p["tau_blab"] - Jbone
        d[IX["CUMBASE"]] = Jbone
        # basal remodelling cancels at OC = OB = 1
        d[IX["BMIN"]] = (p["k_remod"] * (g(IX["OB"]) - g(IX["OC"]))
                         - Jbone / p["BaseCap"])

        d[IX["HCO3_e"]] = (alk_in + NAE - NEAP + Jbuf + Jbone) / VECF
        d[IX["GIVEN"]] = alk_in
        d[IX["WASTE"]] = HCO3_u

        # ---------------- renal adaptation ---------------------------------
        NH3_t = p["NH3P0"] * (1.0
                              + p["kNH3_load"] * (p["NEAP_h"] / p["NEAP_ref_h"] - 1.0)
                              + p["kNH3_ac"] * min(1.0, max(0.0, err) / 6.0)
                              + p["kNH3_K"] * max(0.0, (p["K_pl0"] - g(IX["K_pl"])) / 1.5))
        d[IX["NH3P"]] = (NH3_t - g(IX["NH3P"])) / p["tau_NH3"]
        PEND_t = p["PEND0"] * (1.0 + p["kpend"] * max(0.0, -err) / 3.0)
        d[IX["PEND"]] = (PEND_t - g(IX["PEND"])) / p["tau_pend"]
        NDC1_t = 1.0 + p["kndc1"] * max(0.0, p["pHi_PT0"] - g(IX["pHi_PT"]))
        d[IX["NDC1"]] = (NDC1_t - g(IX["NDC1"])) / p["tau_ndc1"]

        # ---------------- potassium ----------------------------------------
        # Distal Na+ reabsorption that is NOT electrically matched by H+
        # secretion must be matched by K+ secretion instead: the H+-pump lesion
        # IS a K+-wasting lesion.  Hdef = the unsecreted H+ the cortical
        # collecting duct "should" have moved.
        # Hdef == the H+ the collecting duct FAILED to secrete, expressed as the
        # unmet acid load.  Zero in health by construction, so no parameter has
        # to be re-tuned to keep a normal subject in potassium balance.
        Hdef = max(0.0, p["NEAP_h"] - NAE - alk_in)
        K_diet_h = p["K_diet_kg"] * p["BW"] / 24.0
        K_in = p["F_Kdiet"] * K_diet_h + K_abs_K
        # BUG FIX #12: renal K+ secretion was LINEAR in plasma K+ with a weak
        # aldosterone arm, so plasma K+ was almost undefended -- tripling K+
        # intake with K-citrate drove it to 6.6 mmol/L and Na-based alkali drove
        # it to 1.4.  Real K+ homeostasis has very high loop gain (doubling
        # intake moves plasma K+ by ~0.1-0.3 mmol/L).
        Ksec = (p["kK_sec"] * p["K_pl0"] * (g(IX["K_pl"]) / p["K_pl0"]) ** p["nK_sec"]
                * p["BW"] / 24.0
                * (1.0 + p["kK_aldo"] * (g(IX["ALDO"]) - 1.0)) * NEPH
                + p["kK_H"] * Hdef
                + p["kK_Na"] * Na_alk_h
                + p["kK_hctz"] * g(IX["C_hctz"]) / (p["EC50_hctz"] + g(IX["C_hctz"]))
                * p["BW"] / 24.0 * 0.05)
        d[IX["KDEF"]] = Ksec + p["f_stool_K"] * K_diet_h - K_in
        # plasma K reads the deficit only through a ~300 mmol / (mmol/L) gearing
        Kpl_t = (p["K_pl0"] - g(IX["KDEF"]) / (p["kKdef"] * p["BW"])
                 + p["kK_shift"] * (7.40 - pH_bl))
        d[IX["K_pl"]] = (min(7.5, max(1.4, Kpl_t)) - g(IX["K_pl"])) / p["tau_K"]

        ALDO_t = (1.0 + p["kaldo_V"] * max(0.0, (p["VECF0"] - VECF) / p["VECF0"])
                  + p["kaldo_K"] * max(0.0, (g(IX["K_pl"]) - p["K_pl0"]) / 1.0))
        d[IX["ALDO"]] = (ALDO_t - g(IX["ALDO"])) / p["tau_aldo"]
        d[IX["VECF"]] = (p["VECF0"] * (1.0 - 0.05 * max(0.0, err) / 6.0) - VECF) / 24.0
        d[IX["Cl_pl"]] = ((104.0 + (p["HCO3_set"] - HCO3) * 0.85) - g(IX["Cl_pl"])) / 6.0

        # ---------------- calcium / PTH / phosphate ------------------------
        Ca_bone = p["rCaBase"] * Jbone                              # mmol/h
        UCa_basal = p["UCa0"] * p["BW"] / 24.0
        # Three routes, deliberately separated:
        #  (i)  dietary acid load  -> this is what Lemann's 0.035 mmol/mEq slope
        #       actually measures, and it is present with a NORMAL kidney;
        #  (ii) Ca liberated when bone buffers the unmet acid gap (dRTA-specific);
        #  (iii) direct acid inhibition of distal TRPV5 Ca reabsorption.
        # BUG FIX #11: with only (ii) and (iii) the model reproduced a Lemann
        # slope of 0.0015 mmol/mEq -- 23x too flat -- because systemic pH barely
        # moves in a subject whose kidney can still compensate.
        UCa_h = (UCa_basal
                 + p["kUCa_NEAP"] * (NEAP - alk_in - p["NEAP_ref_h"])
                 + p["kUCa_bone"] * Ca_bone
                 + p["kUCa_pH"] * acid_pH * (p["BW"] / 70.0) / 24.0) \
            * (1.0 - p["kUCa_hctz"] * g(IX["C_hctz"]) / (p["EC50_hctz"] + g(IX["C_hctz"])))
        UCa_h = max(0.0, UCa_h)
        # serum Ca stays normal because the bone Ca is spilled into the urine --
        # dRTA is hypercalciuria WITHOUT hypercalcaemia.
        Ca_t = p["Ca_pl0"] + 0.35 * (Ca_bone - max(0.0, UCa_h - UCa_basal))
        d[IX["Ca_pl"]] = (min(1.55, max(0.85, Ca_t)) - g(IX["Ca_pl"])) / p["tau_Ca"]
        PTH_t = p["PTH0"] * (1.0 + p["kPTH_Ca"] * (p["Ca_pl0"] - g(IX["Ca_pl"]))
                             + 0.35 * max(0.0, (25.0 - g(IX["C_25D"])) / 25.0))
        d[IX["PTH"]] = (max(3.0, PTH_t) - g(IX["PTH"])) / p["tau_PTH"]
        d[IX["Pi_pl"]] = ((p["Pi_pl0"] - 0.10 * (g(IX["PTH"]) / p["PTH0"] - 1.0))
                          - g(IX["Pi_pl"])) / 12.0

        # ---------------- bone cells --------------------------------------
        OC_t = 1.0 + p["kOC_ac"] * acid_pH + p["kOC_PTH"] * (g(IX["PTH"]) / p["PTH0"] - 1.0)
        d[IX["OC"]] = (max(0.2, OC_t) - g(IX["OC"])) / p["tau_OC"]
        OB_t = 1.0 - p["kOB_ac"] * acid_pH
        d[IX["OB"]] = (max(0.15, OB_t) - g(IX["OB"])) / p["tau_OB"]
        bALP_t = p["bALP0"] * (1.0 + p["kbALP"] * acid_pH * 0.5) \
            * (1.0 + 0.5 * p["growth_pot"])
        d[IX["bALP"]] = (bALP_t - g(IX["bALP"])) / p["tau_bALP"]
        OSM_t = (p["kOSM"] * acid_pH * 12.0
                 + p["kOSM_D"] * max(0.0, (20.0 - g(IX["C_25D"])) / 20.0))
        d[IX["OSM"]] = (OSM_t - g(IX["OSM"])) / p["tau_OSM"]
        d[IX["BMDz"]] = (p["kz_catch"] * (0.0 - g(IX["BMDz"]))
                         * (0.35 + 0.65 * p["growth_pot"])
                         - p["kz_loss"] * Jbone / (p["BW"] / 70.0))

        # ---------------- urine ------------------------------------------
        UVol_t = (p["UVol_kg"] * p["BW"]
                  * (1.0 + p["kDI"] * max(0.0, (3.8 - g(IX["K_pl"])) / 1.0)))
        d[IX["UVol_s"]] = (UVol_t - g(IX["UVol_s"])) / p["tau_UVol"]
        d[IX["UCa_s"]] = (UCa_h * 24.0 - g(IX["UCa_s"])) / p["tau_UCa"]
        d[IX["UCit_s"]] = (UCit_h * 24.0 - g(IX["UCit_s"])) / p["tau_UCa"]
        d[IX["UpH_s"]] = (uph - g(IX["UpH_s"])) / 12.0

        # brushite (CaHPO4) supersaturation: Ca(free) x HPO4(2-) / Ksp
        Ca_c = (UCa_h * 24.0) / UVol
        Cit_c = (UCit_h * 24.0) / UVol
        Pi_u_c = (Pi_load * 24.0) / UVol
        Ca_free = Ca_c / (1.0 + p["K_cacit"] * Cit_c)
        HPO4 = Pi_u_c / (1.0 + 10.0 ** (p["pKa_Pi"] - uph))
        SS = Ca_free * HPO4 / p["Ksp_ref"]
        d[IX["SS_s"]] = (SS - g(IX["SS_s"])) / p["tau_SS"]

        drive = max(0.0, SS - 1.0) ** p["nNC"]
        d[IX["NC"]] = p["kNC"] * drive - p["kNC_res"] * g(IX["NC"])
        d[IX["STONE"]] = p["kSTONE"] * drive - p["kSTONE_res"] * g(IX["STONE"])
        FIB_t = 1.0 - math.exp(-p["kFIB"] * max(0.0, g(IX["NC"])) * 1000.0)
        d[IX["FIB"]] = (FIB_t - g(IX["FIB"])) / p["tau_FIB"]
        # nephron loss is driven by the fibrosis LEVEL (bug fix #2: it was
        # driven by the target-minus-state GAP, which vanishes at steady state
        # and made established nephrocalcinosis harmless)
        d[IX["NEPH"]] = -p["kNEPH"] * max(0.0, g(IX["FIB"])) * NEPH

        # ---------------- growth / muscle / adherence ---------------------
        acid_idx = max(0.0, p["HCO3_set"] - HCO3)
        GM = 1.0 / (1.0 + (acid_idx / p["K_ac_igf"]) ** p["h_igf"])
        IGF_t = p["IGF1_ref"] * GM
        d[IX["IGF1"]] = (IGF_t - g(IX["IGF1"])) / p["tau_IGF"]
        gp = p["growth_pot"]
        d[IX["Hz"]] = (p["kh_catch"] * (0.0 - g(IX["Hz"])) * gp * GM
                       - p["kh_drag"] * (1.0 - GM) * gp)
        MUS_t = 1.0 - p["kMUS"] * max(0.0, (3.5 - g(IX["K_pl"])) / 1.5) ** 1.5
        d[IX["MUS"]] = (max(0.05, MUS_t) - g(IX["MUS"])) / p["tau_MUS"]
        # GI irritation from IR alkali bolus (osmotic + CO2 release)
        # BUG FIX #13: the GI-irritation drive was normalised by BW/70 and then
        # divided by 25, so ANY immediate-release bolus saturated it at the
        # min(1, GI) cap and gastrointestinal tolerability became a binary
        # switch rather than a dose-graded penalty.  Drive is now the bolus
        # absorption rate PER KG, which is what a patient actually feels.
        d[IX["GI"]] = p["kGI"] * (a_bicIR + 3.0 * a_citIR) / p["BW"] \
            - g(IX["GI"]) / p["tau_GI"]
        ADH_t = p["ADH_ref"] * math.exp(-p["kADH_n"] * max(0.0, p["n_intakes"] - 2.0)) \
            * (1.0 - p["kADH_GI"] * min(1.0, g(IX["GI"])))
        d[IX["ADH"]] = (ADH_t - g(IX["ADH"])) / p["tau_ADH"]
        d[IX["HEAR"]] = p["kHEAR"] * max(0.0, p["HEAR_max"] - g(IX["HEAR"]))

        # ---------------- exposure integrals ------------------------------
        d[IX["TBT"]] = 1.0 if HCO3 < 22.0 else 0.0
        d[IX["AAC"]] = max(0.0, 24.0 - HCO3)

        if store:
            self.last = dict(pH_bl=pH_bl, uph=uph, NAE=NAE * 24.0, TA=TA * 24.0,
                             NH4=NH4_u * 24.0, HCO3u=HCO3_u * 24.0,
                             NEAP=NEAP * 24.0, alk_in=alk_in * 24.0,
                             Jbone=Jbone, VHc=VHc, SS=SS, GM=GM,
                             UCa=UCa_h * 24.0, UCit=UCit_h * 24.0,
                             Ca_cit=(UCa_h * 24.0 * 40.08) / max(0.02, UCit_h * 24.0 * 192.12),
                             gap=NEAP - alk_in - NAE, adh=adh,
                             FE_HCO3=100.0 * HCO3_u / max(1e-9, GFRh * HCO3),
                             UVol=UVol, Pi_u=Pi_load * 24.0, Hdef=Hdef,
                             cit_filt=cit_filt * 24.0, cit_reab=cit_reab * 24.0)
        return d


# ----------------------------------------------------------------------------
# integrator: RK4 with dose events, adherence applied at dose time
# ----------------------------------------------------------------------------
def simulate(par, reg, t_end, dt=0.25, record_every=6.0, seed=1):
    m = Model(par)
    p = m.p
    p["n_intakes"] = reg.n_intakes() if reg.kind != "none" else 2.0
    p["f_Kcat"] = {"K": 1.0, "Na": 0.0, "NaK": 0.5}[getattr(reg, "cation", "K")]
    y = m.init()
    ev = dose_events(reg, p, t_end)
    ei = 0
    rng = _LCG(seed)
    out = []
    t = 0.0
    nrec = max(1, int(round(record_every / dt)))
    step = 0
    while t < t_end - 1e-9:
        # apply dose events at this instant; ONE adherence draw per dose TIME
        # (bug fix #4: an independent coin per event let a patient take the
        # citrate granules but skip the bicarbonate granules of the same dose)
        if ei < len(ev) and ev[ei][0] <= t + 1e-9:
            take = 1.0 if rng.uniform() < min(1.0, max(0.0, y[IX["ADH"]])) else 0.0
            while ei < len(ev) and ev[ei][0] <= t + 1e-9:
                _, idx, amt = ev[ei]
                y[idx] += amt * take
                ei += 1
        if step % nrec == 0:
            dd = m.rhs(t, y, store=True)
            out.append((t, list(y), dict(m.last)))
        h = dt
        nxt = ev[ei][0] if ei < len(ev) else 1e18
        if t + h > nxt:
            h = max(1e-4, nxt - t)
        k1 = m.rhs(t, y)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(NS)]
        k2 = m.rhs(t + 0.5 * h, y2)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(NS)]
        k3 = m.rhs(t + 0.5 * h, y3)
        y4 = [y[i] + h * k3[i] for i in range(NS)]
        k4 = m.rhs(t + h, y4)
        for i in range(NS):
            y[i] += h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
        # hard non-negativity on amounts / pools
        for i in (0, 1, 2, 3, 4, 5, 7, 9):
            if y[i] < 0.0:
                y[i] = 0.0
        if y[IX["CIT_pl"]] < 1e-8:
            y[IX["CIT_pl"]] = 1e-8
        if y[IX["VH"]] > p["VH_max"]:
            y[IX["VH"]] = p["VH_max"]
        if y[IX["VH"]] < 0.0:
            y[IX["VH"]] = 0.0
        t += h
        step += 1
    dd = m.rhs(t, y, store=True)
    out.append((t, list(y), dict(m.last)))
    return out, m


class _LCG:
    """Tiny deterministic RNG (no numpy)."""
    def __init__(self, seed=1):
        self.s = seed * 2654435761 % 2147483647 or 12345

    def uniform(self):
        self.s = (self.s * 48271) % 2147483647
        return self.s / 2147483647.0

    def normal(self):
        u1 = max(1e-12, self.uniform())
        u2 = self.uniform()
        return math.sqrt(-2 * math.log(u1)) * math.cos(2 * math.pi * u2)
