#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mmd_reference_model.py
======================
Reference (numerical truth) implementation of the Moyamoya Disease (MMD) QSP
model.  This file exists so that EVERY number quoted in README.md,
mmd_mrgsolve_model.R and mmd_shiny_app.R is *computed*, not asserted.  The
mrgsolve model is a line-by-line translation of the equations below.

--------------------------------------------------------------------------
THE THESIS
--------------------------------------------------------------------------
In moyamoya disease the measured quantity (cerebral blood flow, CBF) is not
the state variable.  CBF is *defended* by the cortical arteriole, which is a
variable resistor with a floor.  The state variable is how much of that
resistor's dilatory range is left -- the cerebrovascular reserve.

Write one hemisphere as a 2-node resistive network.  Territory A is the MCA
cortex distal to the stenosed terminal ICA; territory B is the posterior
(PCA) cortex, which is also the *source* of the leptomeningeal collateral
that keeps A alive:

      Pa --[ R_ica(stenosis) ]--+
      Pa --[ R_moya           ]--+--> P_A --[ R_artA ]--> Pv
      Pa --[ R_byp (surgery)  ]--+          (autoregulator, has a FLOOR)
                                 |
                          [ R_coll ]
                                 |
      Pa --[ R_pca ]---------> P_B --[ R_artB ]--> Pv

Nodal balance (Pv = ICP):
   (gi+gm+gb)(Pa-P_A) + gc(P_B-P_A) = gA(P_A-Pv)
             gp(Pa-P_B)             = gB(P_B-Pv) + gc(P_B-P_A)

Autoregulation chooses gA, gB to meet demand, CLAMPED to
[1/R_art_max_eff, 1/R_art_min].  Nine results follow as arithmetic:

 (1) A CRITICAL INLET CONDUCTANCE exists.  Below it the arteriole is on its
     floor and CBF becomes pressure-passive: dCBF/dMAP jumps from ~0 to
     1/R_total.  Stenosis severity alone does not locate this point --
     collateral conductance does.
 (2) THE PENUMBRAL THRESHOLD IS NOT A PARAMETER.  It is
     CBF_crit = CMRO2 / (CaO2 * OEF_max).  With Hb 15 g/dL that is
     ~19.7 mL/100g/min (the textbook 20).  With Hb 8 g/dL (sickle-cell
     moyamoya) it is ~36.9 -- the same brain infarcts at flows a normal
     brain tolerates, and transfusion, not a vasodilator, is what moves it.
 (3) CO2 REACTIVITY CHANGES SIGN, and the sign is a direct read-out of
     reserve = 0.  Because B is A's supply, dilating B steals from A.
     Intrinsic reserve and MEASURED reactivity are different numbers.
 (4) HYPERVENTILATION IS AN ISCHAEMIC STIMULUS ONLY WHERE RESERVE IS ZERO
     -- the crying child.  Same PaCO2 step, opposite outcome by reserve.
 (5) THE BLOOD-PRESSURE PARADOX IS A SINGLE OPTIMISATION.  Lowering MAP
     scales ischaemic hazard up through the pressure-passive term and
     haemorrhagic hazard down through perforator wall stress
     sigma ~ P_perf * d^1.5.  The optimum MAP is computable and differs by
     phenotype -- which is the JAM-trial result written as calculus.
 (6) BYPASS IS A PARALLEL CONDUCTANCE, so it treats the ischaemic and the
     haemorrhagic phenotype BY THE SAME MECHANISM: raising P_A both lifts
     CBF_A and collapses the gradient across the fragile perforators, which
     then regress.
 (7) HYPERPERFUSION IS PRICED IN PRE-OPERATIVE ISCHAEMIA, NOT IN SURGERY.
     A chronically floored arteriole remodels and loses its CONSTRICTIVE
     range (REMOD).  Normal pressure then arrives at a vessel that cannot
     constrict.  The duration of the syndrome is REMOD's time constant.
 (8) INDIRECT BYPASS IS A BET ON ANGIOGENIC CAPACITY, so the same operation
     is curative in a child and near-useless in an adult -- one parameter.
 (9) RNF213 ENTERS TWICE WITH OPPOSITE SIGN: it drives intimal SMC growth
     (raises stenosis) AND caps the collateral ceiling.  A gene that makes
     the disease and forbids the compensation.

Author: QSP Disease Model Library (Claude Code Routine), 2026-08-06
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# =========================================================================
# 0.  PARAMETERS
# =========================================================================
P = dict(
    # ---- systemic / geometry -------------------------------------------
    MAP0      = 90.0,    # mmHg, baseline mean arterial pressure
    ICP       = 10.0,    # mmHg, taken as venous outflow pressure
    MASS_A    = 250.0,   # g, MCA (anterior) territory
    MASS_B    = 120.0,   # g, PCA (posterior) territory
    HB        = 15.0,    # g/dL haemoglobin
    SAO2      = 0.98,    # arterial saturation
    CMRO2_0   = 3.30,    # mL O2 / 100 g / min  (cortical grey+white mix)
    OEF_MAX   = 0.85,    # hard ceiling on oxygen extraction
    OEF_BASE  = 0.335,

    # ---- baseline resistive split --------------------------------------
    CBF_A0    = 50.0,    # mL/100g/min  normal
    CBF_B0    = 50.0,
    FRAC_PROX = 0.25,    # fraction of territory resistance that is large-artery
    DIL_MAX   = 2.50,    # R_art_min = R_art0 / DIL_MAX   (max dilatation)
    CON_MAX   = 2.20,    # R_art_max = R_art0 * CON_MAX   (max constriction)
    GCOLL0    = 0.050,   # mL/min/mmHg, native leptomeningeal conductance
    GMOYA0    = 0.004,   # native basal perforator conductance

    # ---- RNF213 / intimal hyperplasia (the lesion) ---------------------
    RNF       = 1.0,     # 0 = wild type, 1 = R4810K heterozygote, 1.6 = homozygote
    SEC       = 0.0,     # RNF213-independent drive: irradiation, NF1, SCD,
                         # Down syndrome, Graves -- i.e. quasi-moyamoya
    K_PDGF    = 0.90,    # /d   PDGF-BB production gain
    K_PDGF_D  = 0.60,    # /d   PDGF decay
    K_SMC     = 2.2e-4,  # /d   intimal SMC accumulation rate constant
    SMC_MAX   = 1.0,
    STEN_MAX  = 0.945,   # asymptotic fractional DIAMETER reduction
    #  STEN = STEN_MAX * SMC ; ALL of the non-linearity lives in the
    #  Poiseuille law R ~ (1-STEN)^-4, which is where it belongs.
    K_STATIN  = 0.28,    # fractional suppression of K_SMC at full statin effect

    # ---- collateral / angiogenesis -------------------------------------
    ANGIO     = 1.0,     # angiogenic capacity: ~1.0 child, ~0.35 adult
    K_HIF_ON  = 6.0,     # /d
    K_HIF_OFF = 3.0,     # /d
    K_VEGF    = 2.0,     # /d
    K_VEGF_D  = 1.4,     # /d
    K_MMP     = 1.5, K_MMP_D = 1.0,
    K_MOYA    = 0.018,   # /d  basal perforator recruitment
    GMOYA_CAP = 3.40,    # mL/min/mmHg ceiling BEFORE the RNF213 penalty
    RNF_CAP   = 0.42,    # collateral ceiling multiplier per unit RNF
    K_MOYA_REG= 0.030,   # /d  regression when the driving gradient collapses
    K_COLL    = 0.009, GCOLL_CAP = 1.00, K_COLL_REG = 0.012,

    # ---- the periventricular / choroidal anastomosis: THE BLEEDER -------
    # Two anatomically different collateral routes, with opposite safety
    # profiles.  The leptomeningeal route is donated by the PCA and is safe.
    # The periventricular (choroidal / thalamotuberal) route is recruited
    # when the DEEP borderzone is starved, is thin-walled, and is the
    # documented source of haemorrhage (Funaki, Takahashi).  PCA involvement
    # therefore does two opposite things at once: it removes the safe donor
    # and forces the dangerous route.
    PCA_INV   = 0.0,     # 0 = PCA spared, 1 = PCA involved (~30% of MMD)
    PCA_STEN  = 0.25,    # diameter stenosis of the PCA at PCA_INV = 1
    PCA_COLL  = 0.60,    # fractional loss of leptomeningeal ceiling
    PCA_PVA   = 2.20,    # gain on the periventricular ceiling
    GPVA0     = 0.002, GPVA_CAP = 0.55, K_PVA = 0.014, K_PVA_REG = 0.030,
    GPVA_RREF = 0.050,

    # ---- surgery --------------------------------------------------------
    GBYP_DIR  = 1.05,    # mL/min/mmHg, mature STA-MCA direct anastomosis
    TAU_BYP   = 6.0,     # d, direct graft maturation time constant
    GBYP_IND  = 1.15,    # ceiling of an indirect (EDAS/EMS) construct
    K_BYP_IND = 0.055,   # /d, VEGF- and ANGIO-gated indirect ingrowth

    # ---- arteriolar remodelling (the hyperperfusion mechanism) ---------
    K_REM_ON  = 0.075,   # /d, gain of chronic-dilatation remodelling
    K_REM_OFF = 0.050,   # /d  -> tau ~ 20 d recovery of constrictive range
    # VASOPARALYSIS.  An arteriole held at maximal dilatation for months loses
    # essentially ALL of its constrictive range, not merely part of it -- which
    # is what the hyperperfusion literature means by the word.  REM_CON is set
    # so that REMOD = 1 collapses the whole [g_con, g_dil] interval:
    #   REM_CON = GA_DIL/GA_CON = DIL_MAX * CON_MAX = 2.5 * 2.2 = 5.50
    REM_CON   = 5.50,

    # ---- the peri-anastomotic compartment -------------------------------
    # A LUMPED TERRITORY CANNOT HYPERPERFUSE.  Even a patent graft leaves the
    # territorial arteriole still dilating (gS after bypass is ~2.4, and
    # holding 50 mL/100g/min through it needs g_art ABOVE resting), so mean
    # territorial CBF never exceeds demand.  Post-bypass hyperperfusion is
    # therefore necessarily FOCAL: it belongs to the cortex the graft is sewn
    # to, which sees a pressure close to systemic while its own arterioles are
    # vasoparalysed.  That requires a third node.
    F_FRAC    = 0.10,    # peri-anastomotic share of territory A
    G_LEAK    = 0.80,    # pial/M3 coupling between focal cortex and territory
    #  Coupling is GOOD -- that is why a direct graft helps the whole
    #  territory and not just the cortex it is sewn to.  The peri-anastomotic
    #  node still sits at a much higher pressure, because the graft flow it
    #  does not consume has to cross gl to reach the territory:
    #  P_F - P_A = (graft flow passed on) / gl.

    # ---- injury ---------------------------------------------------------
    # ---- watershed geometry --------------------------------------------
    # Collateral flow arrives at the cortical SURFACE and travels inward, so
    # the deep borderzone is supplied last.  The watershed penalty therefore
    # scales with how COLLATERAL-dependent the supply is, not with how low
    # the mean flow is.  This is why a territory with an acceptable mean CBF
    # still infarcts, and why the infarcts are where they are.
    WS_FRAC   = 0.30,    # fraction of territory A that is borderzone
    WS_K      = 0.55,    # max watershed penalty at fully collateral supply
    W_BYP     = 0.80,    # a direct graft is 80% as good as antegrade flow
    K_INF     = 0.0165,  # /d, infarct accrual gain
    K_INF_EMB = 4.0e-5,  # /d per unit embolic hazard
    # Clinical ischaemic events are threshold CROSSINGS by a fluctuating
    # perfusion, not a state of permanent sub-threshold flow.  The natural
    # form is therefore exponential in the MARGIN to CBF_crit, with a scale
    # equal to the day-to-day perfusion variability -- and that variability
    # is itself larger in a pressure-passive territory, because there MAP
    # noise passes straight through to CBF.
    ISCH_HAZ0 = 0.0480,  # /yr at zero margin
    CBF_SD0   = 4.50,    # mL/100g/min, autoregulating perfusion variability
    K_PASS    = 0.70,    # extra variability when the arteriole is on its floor
    ISCH_EMB  = 0.075,   # /yr per unit embolic hazard
    VEGF_REG  = 0.60,    # VEGF below which the collateral bed is pruned
    K_TIA     = 1.0,
    K_EMB     = 0.35,    # embolic hazard gain from low-shear moyamoya vessels
    K_EMB_ASA = 0.55,    # fractional reduction at full platelet inhibition

    # ---- perforator wall / haemorrhage ---------------------------------
    SIG0      = 1.0,     # normalising constant, set in init()
    K_ANEU    = 0.0110,  # /d microaneurysm formation gain
    K_ANEU_REP= 0.0060,  # /d repair
    K_MINO    = 0.45,    # minocycline (MMP-9) boost to repair
    HEM_HAZ0  = 0.4140,  # /yr haemorrhage hazard at ANEU=1, sigma=1
    HEM_POW   = 2.0,
    HEM_VOL   = 26.0,    # mL, mean haematoma volume per event

    # ---- cognition / disability ----------------------------------------
    K_COG     = 0.00090, # /d z-score loss per unit relative hypoperfusion
    K_COG_INF = 3.10,    # z per unit infarct fraction

    # ---- BBB / oedema ---------------------------------------------------
    #  HYPERPERFUSION IS DEFINED RELATIVE TO THE FLOW THE BARRIER HAS
    #  ADAPTED TO, not to a textbook normal.  A cortex perfused at 18
    #  mL/100g/min for years is injured by 50, which is normal.  CBFAD is that
    #  adapted set-point and it re-adapts with TAU_ADAPT, which is therefore
    #  what sets the DURATION of the syndrome (2-3 weeks), while REMOD
    #  (vasoparalysis) sets its HEIGHT.
    TAU_ADAPT = 22.0,    # d, barrier/autoregulatory re-adaptation
    K_BBB     = 0.45, K_BBB_REP = 0.25, HYPER_THR = 1.35,
    K_EDEMA   = 0.85, K_EDEMA_R = 0.30,

    # ---- drug PK --------------------------------------------------------
    # aspirin: irreversible COX-1, effect compartment is platelet pool
    ASA_KA = 12.0, ASA_KE = 8.3,  ASA_V = 12.0,
    K_PLT_ON = 45.0, K_PLT_OFF = 0.10,          # 10 %/d platelet turnover
    # cilostazol: PDE3 inhibitor, t1/2 ~ 11 h
    CILO_KA = 4.2, CILO_KE = 1.51, CILO_V = 95.0, CILO_EC50 = 0.55, CILO_EMAX = 0.22,
    CILO_APL = 0.45,                             # platelet effect fraction
    # nifedipine GITS: systemic + NON-SELECTIVE cerebral arteriolar dilator
    NIF_KA = 0.55, NIF_KE = 2.31, NIF_V = 55.0, NIF_EC50 = 0.030,
    NIF_EMAX_CBR = 0.30,                         # cerebral dilatory stimulus
    NIF_EMAX_MAP = 0.16,                         # systemic MAP reduction
    # statin (atorvastatin), pleiotropic endothelial effect
    STAT_KA = 6.0, STAT_KE = 0.50, STAT_V = 380.0, STAT_EC50 = 0.012,
    # minocycline
    MIN_KA = 5.0, MIN_KE = 0.77, MIN_V = 100.0, MIN_EC50 = 1.20,
    # acetazolamide (diagnostic challenge)
    AZ_KE = 3.30, AZ_V = 18.0, AZ_EC50 = 8.0, AZ_EMAX = 0.45,
    # generic antihypertensive (ARB/CCB class effect on MAP)
    AHT_KA = 8.0, AHT_KE = 1.0, AHT_V = 60.0, AHT_EC50 = 0.10, AHT_EMAX = 0.22,

    # ---- CO2 -------------------------------------------------------------
    PACO2_0 = 40.0,
    K_CO2      = 0.0350,   # arteriolar DEMAND change per mmHg PaCO2
    K_CO2_COLL = 0.0200,   # leptomeningeal/pial COLLATERAL conductance per mmHg
    K_CO2_MOYA = 0.0110,   # basal perforator conductance per mmHg (less reactive)
    #  Acetazolamide acts on the arteriolar term ONLY (tissue carbonic
    #  anhydrase), PaCO2 acts on BOTH.  The two probes are therefore not
    #  interchangeable in moyamoya -- see report section 3/4.
    MAP_REF = 90.0,        # fixed normalising pressure for wall stress
    GMOYA_RREF = 0.20,     # perforator conductance above which "dilated"
    Q_LOAD_REF = 125.0,    # mL/min, native design flow of territory A
)

# names of the ODE state vector, in order
SNAMES = [
    "PDGF", "SMC", "GMOYA", "GPVA", "GCOLL", "GBYP", "REMOD", "HIF", "VEGF",
    "MMP9",
    "INFA", "INFB", "ANEU", "HEMV", "COG", "BBB", "EDEMA", "EMBH", "TIAB",
    "HYPOX", "HEMH", "ISCH", "CBFAD",
    "ASA_G", "ASA_C", "PLTI", "CILO_G", "CILO_C", "NIF_G", "NIF_C",
    "STAT_G", "STAT_C", "MIN_G", "MIN_C", "AZ_C", "AHT_G", "AHT_C",
]
IX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)


# =========================================================================
# 1.  DERIVED CONSTANTS
# =========================================================================
def init(p):
    """Compute the derived resistive geometry and O2 constants."""
    d = dict(p)
    d["CAO2"] = 1.34 * p["HB"] * p["SAO2"] / 100.0        # mL O2 / mL blood
    d["CBF_CRIT"] = p["CMRO2_0"] / (d["CAO2"] * p["OEF_MAX"])   # mL/100g/min
    cpp = p["MAP0"] - p["ICP"]
    QA0 = p["CBF_A0"] * p["MASS_A"] / 100.0               # mL/min
    QB0 = p["CBF_B0"] * p["MASS_B"] / 100.0
    d["QA0"], d["QB0"] = QA0, QB0
    RtotA, RtotB = cpp / QA0, cpp / QB0
    d["R_ICA0"] = p["FRAC_PROX"] * RtotA
    d["R_ARTA0"] = (1 - p["FRAC_PROX"]) * RtotA
    d["R_PCA0"] = p["FRAC_PROX"] * RtotB
    d["R_ARTB0"] = (1 - p["FRAC_PROX"]) * RtotB
    d["GA0"] = 1 / d["R_ARTA0"]
    d["GB0"] = 1 / d["R_ARTB0"]
    d["GA_DIL"] = p["DIL_MAX"] / d["R_ARTA0"]     # max dilatation conductance
    d["GA_CON"] = 1 / (d["R_ARTA0"] * p["CON_MAX"])
    d["GB_DIL"] = p["DIL_MAX"] / d["R_ARTB0"]
    d["GB_CON"] = 1 / (d["R_ARTB0"] * p["CON_MAX"])
    # focal (peri-anastomotic) and territorial masses / demands
    d["MASS_F"] = p["MASS_A"] * p["F_FRAC"]
    d["MASS_At"] = p["MASS_A"] * (1 - p["F_FRAC"])
    d["QF0"] = QA0 * p["F_FRAC"]
    d["QAt0"] = QA0 * (1 - p["F_FRAC"])
    d["GF0"] = d["GA0"] * p["F_FRAC"]
    d["GAt0"] = d["GA0"] * (1 - p["F_FRAC"])
    d["GF_DIL"] = d["GA_DIL"] * p["F_FRAC"]
    d["GF_CON"] = d["GA_CON"] * p["F_FRAC"]
    d["GAt_DIL"] = d["GA_DIL"] * (1 - p["F_FRAC"])
    d["GAt_CON"] = d["GA_CON"] * (1 - p["F_FRAC"])
    d["GI0"] = 1 / d["R_ICA0"]
    d["GP0"] = 1 / d["R_PCA0"] * (1 - p["PCA_STEN"] * p["PCA_INV"]) ** 4
    # reference wall stress of the native perforator bed, evaluated at the
    # HEALTHY operating point: P_A0 = Pv + CPP*(1-FRAC_PROX), dilatation = 1
    cpp_ref = p["MAP_REF"] - p["ICP"]
    d["P_A0"] = p["ICP"] + cpp_ref * (1 - p["FRAC_PROX"])
    d["SIG_REF"] = 0.5 * (p["MAP_REF"] + d["P_A0"])   # MAP-INDEPENDENT
    return d


# =========================================================================
# 2.  HAEMODYNAMIC CORE  (algebraic, quasi-steady autoregulation)
# =========================================================================
def solve_network(Pa, Pv, gSA, gSF, gb, gp, gc, gl, gA, gF, gB):
    """
    Three-node linear solve.

        Pa --[gSA]--+                 Pa --[gSF + gb]--+
                    v                                  v
        Pb --[gc]--> P_A <---[gl]---> P_F        (peri-anastomotic cortex)
                     |                 |
                   [gA]              [gF]
                     v                 v
                     Pv                Pv

    A: gSA(Pa-P_A) + gc(P_B-P_A) + gl(P_F-P_A) = gA(P_A-Pv)
    F: (gSF+gb)(Pa-P_F) + gl(P_A-P_F)          = gF(P_F-Pv)
    B: gp(Pa-P_B) + gc(P_A-P_B)                = gB(P_B-Pv)

    With gb = 0 and gSA:gSF and gA:gF both split by mass, P_A == P_F exactly,
    so the focal compartment is invisible until an operation happens.
    Returns (P_A, P_F, P_B).
    """
    gSFb = gSF + gb
    M = np.array([
        [-(gSA + gc + gl + gA), gl,                  gc],
        [gl,                    -(gSFb + gl + gF),   0.0],
        [gc,                    0.0,                 -(gp + gc + gB)],
    ])
    r = np.array([-(gSA * Pa + gA * Pv),
                  -(gSFb * Pa + gF * Pv),
                  -(gp * Pa + gB * Pv)])
    P = np.linalg.solve(M, r)
    return P[0], P[1], P[2]


def autoreg(d, Pa, gSA, gSF, gb, gp, gc, gl, QtA, QtF, QtB,
            gAc, gAd, gFc, gFd, gBc, gBd, iters=30):
    """
    Quasi-steady autoregulation on the three-node network.  Each compartment
    picks the arteriolar conductance that would deliver its demanded flow,
    CLAMPED to its own dilatory floor and constrictive ceiling.  Fixed-point
    iterated because the compartments share conduits.
    """
    Pv = d["ICP"]
    gA, gF, gB = d["GAt0"], d["GF0"], d["GB0"]
    PA, PF, PB = solve_network(Pa, Pv, gSA, gSF, gb, gp, gc, gl, gA, gF, gB)
    for _ in range(iters):
        PAr = (gSA * Pa + gc * PB + gl * PF - QtA) / (gSA + gc + gl)
        gA = QtA / (PAr - Pv) if PAr - Pv > 1e-9 else gAd
        gA = min(max(gA, gAc), gAd)
        PFr = ((gSF + gb) * Pa + gl * PA - QtF) / (gSF + gb + gl)
        gF = QtF / (PFr - Pv) if PFr - Pv > 1e-9 else gFd
        gF = min(max(gF, gFc), gFd)
        PBr = (gp * Pa + gc * PA - QtB) / (gp + gc)
        gB = QtB / (PBr - Pv) if PBr - Pv > 1e-9 else gBd
        gB = min(max(gB, gBc), gBd)
        PA, PF, PB = solve_network(Pa, Pv, gSA, gSF, gb, gp, gc, gl,
                                   gA, gF, gB)
    return dict(gA=gA, gF=gF, gB=gB, PA=PA, PF=PF, PB=PB,
                QA=gA * (PA - Pv), QF=gF * (PF - Pv), QB=gB * (PB - Pv))


def hemodynamics(d, y, t, sc, full=True):
    """
    Full algebraic layer at one instant: drug effects -> demand -> network
    -> flows, OEF, CMRO2, reserve, wall stress.  Returns a flat dict.
    `full=False` skips the three reserve/pressure-passivity probes (each of
    which is itself a network solve) and is what the ODE right-hand side uses.
    """
    Pv = d["ICP"]
    STEN = d["STEN_MAX"] * np.clip(y[IX["SMC"]], 0.0, 1.0)

    # ---- CO2: acts on the COLLATERAL CONDUITS as well as the arteriole ---
    paco2 = sc.get("paco2", lambda tt: d["PACO2_0"])(t)
    dco2 = paco2 - d["PACO2_0"]
    f_coll = float(np.clip(1.0 + d["K_CO2_COLL"] * dco2, 0.25, 1.8))
    f_moya = float(np.clip(1.0 + d["K_CO2_MOYA"] * dco2, 0.35, 1.6))

    # ---- inlet conductances ------------------------------------------
    gi = d["GI0"] * (1.0 - STEN) ** 4
    gm = max(y[IX["GMOYA"]], 0.0) * f_moya
    gv = max(y[IX["GPVA"]], 0.0) * f_moya          # periventricular route
    gb = max(y[IX["GBYP"]], 0.0)          # a surgical graft is not CO2-gated
    gc = max(y[IX["GCOLL"]], 0.0) * f_coll
    gS = gi + gm + gv + gb
    gp = d["GP0"]
    gl = d["G_LEAK"]
    # native inlets are shared between the territorial and peri-anastomotic
    # cortex in proportion to mass, so the focal compartment is invisible
    # until a graft is added to it
    gnat = gi + gm + gv
    gSA = gnat * (1 - d["F_FRAC"])
    gSF = gnat * d["F_FRAC"]

    # ---- drug concentrations / effects --------------------------------
    c_cilo = y[IX["CILO_C"]] / d["CILO_V"]
    c_nif = y[IX["NIF_C"]] / d["NIF_V"]
    c_az = y[IX["AZ_C"]] / d["AZ_V"]
    c_aht = y[IX["AHT_C"]] / d["AHT_V"]
    e_cilo = d["CILO_EMAX"] * c_cilo / (d["CILO_EC50"] + c_cilo)
    e_nif = d["NIF_EMAX_CBR"] * c_nif / (d["NIF_EC50"] + c_nif)
    e_nif_map = d["NIF_EMAX_MAP"] * c_nif / (d["NIF_EC50"] + c_nif)
    e_az = d["AZ_EMAX"] * c_az / (d["AZ_EC50"] + c_az)
    e_aht = d["AHT_EMAX"] * c_aht / (d["AHT_EC50"] + c_aht)

    # ---- systemic pressure --------------------------------------------
    # a MAP policy applies only from the day it is started -- otherwise a
    # "post-operative BP target" would silently rewrite the whole pre-operative
    # history and make the hemisphere a different one
    mm = sc.get("map_mult", 1.0) if t >= sc.get("map_mult_t", -1e9) else 1.0
    Pa = d["MAP0"] * (1 - e_nif_map) * (1 - e_aht) * mm

    # ---- arteriolar demand: CO2 + acetazolamide + drug vasodilatation --
    fco2 = 1.0 + d["K_CO2"] * dco2
    dem = max(fco2, 0.35) * (1 + e_az) * (1 + e_cilo) * (1 + e_nif)
    QtA = d["QAt0"] * dem * (1 - y[IX["INFA"]])  # infarcted tissue has no demand
    QtF = d["QF0"] * dem
    QtB = d["QB0"] * dem * (1 - y[IX["INFB"]])

    # ---- arteriolar clamps; the CONSTRICTIVE one is widened by chronic
    #      remodelling (vasoparalysis), which is the hyperperfusion mechanism
    rem = float(np.clip(y[IX["REMOD"]], 0, 1))
    gA_dil, gF_dil = d["GAt_DIL"], d["GF_DIL"]
    gA_con = min(d["GAt_CON"] * (1 + (d["REM_CON"] - 1) * rem), gA_dil * 0.999)
    gF_con = min(d["GF_CON"] * (1 + (d["REM_CON"] - 1) * rem), gF_dil * 0.999)
    gB_con, gB_dil = d["GB_CON"], d["GB_DIL"]

    s = autoreg(d, Pa, gSA, gSF, gb, gp, gc, gl, QtA, QtF, QtB,
                gA_con, gA_dil, gF_con, gF_dil, gB_con, gB_dil)

    massA = d["MASS_At"] * (1 - y[IX["INFA"]]) + 1e-6
    massF = d["MASS_F"]
    massB = d["MASS_B"] * (1 - y[IX["INFB"]]) + 1e-6
    cbfA = s["QA"] / massA * 100.0          # mL/100g/min of SURVIVING tissue
    cbfF = s["QF"] / massF * 100.0          # peri-anastomotic cortex
    cbfB = s["QB"] / massB * 100.0
    # what a SPECT of the WHOLE MCA territory would report
    cbfT = (s["QA"] + s["QF"]) / (massA + massF) * 100.0

    # ---- oxygen: OEF rises to a ceiling, then CMRO2 falls -------------
    del_A = cbfA * d["CAO2"]                # mL O2 /100 g /min delivered
    oefA = float(np.clip(d["CMRO2_0"] / max(del_A, 1e-9),
                         d["OEF_BASE"], d["OEF_MAX"]))
    cmroA = del_A * oefA
    del_B = cbfB * d["CAO2"]
    oefB = float(np.clip(d["CMRO2_0"] / max(del_B, 1e-9),
                         d["OEF_BASE"], d["OEF_MAX"]))
    cmroB = del_B * oefB

    # ---- watershed (borderzone) perfusion ------------------------------
    # Antegrade FLOW fraction, not conductance fraction: the graft reaches the
    # territorial cortex through the pial coupling gl, so its contribution is
    # whatever actually crosses that connection.
    q_in_A = max(gSA * (Pa - s["PA"]) + gc * (s["PB"] - s["PA"])
                 + gl * (s["PF"] - s["PA"]), 1e-9)
    q_ante = gi * (1 - d["F_FRAC"]) * (Pa - s["PA"]) \
        + d["W_BYP"] * max(gl * (s["PF"] - s["PA"]), 0.0)
    ante = float(np.clip(q_ante / q_in_A, 0.0, 1.0))
    k_ws = 1.0 - d["WS_K"] * (1.0 - ante)
    cbf_ws = cbfA * k_ws
    oef_ws = float(np.clip(d["CMRO2_0"] / max(cbf_ws * d["CAO2"], 1e-9),
                           d["OEF_BASE"], d["OEF_MAX"]))

    # ---- perforator wall stress and flow load -------------------------
    # G ~ N*r^4 and both N and r grow, so r ~ (G/Gref)^(1/8).
    # Laplace with a wall that thins as it dilates (h ~ r^-1/2):
    #   sigma ~ P_transmural * r / h ~ P * dil^1.5
    P_perf = 0.5 * (Pa + s["PA"])
    PF = s["PF"]
    dil = max(gm / d["GMOYA_RREF"], 1.0) ** 0.125
    sigma = P_perf * dil ** 1.5 / d["SIG_REF"]
    q_moya = gm * (Pa - s["PA"])
    load = q_moya / d["Q_LOAD_REF"]
    # the periventricular anastomosis: thinner walled, same transmural
    # pressure, and it is THIS stress that the haemorrhage hazard reads
    dil_v = max(gv / d["GPVA_RREF"], 1.0) ** 0.125
    sig_v = P_perf * dil_v ** 1.5 / d["SIG_REF"]
    q_pva = gv * (Pa - s["PA"])
    load_v = q_pva / d["Q_LOAD_REF"]

    out = dict(STEN=STEN, gi=gi, gm=gm, gb=gb, gc=gc, gS=gS, Pa=Pa, gl=gl,
               PA=s["PA"], PF=PF, PB=s["PB"], QA=s["QA"], QF=s["QF"],
               QB=s["QB"], gA=s["gA"], gF=s["gF"], gB=s["gB"],
               gA_dil=gA_dil, gA_con=gA_con, gF_dil=gF_dil, gF_con=gF_con,
               CBFF=cbfF, CBFT=cbfT, AR_POS_F=s["gF"] / gF_dil,
               HYPER=cbfF / d["CBF_A0"],
               HYPER_REL=cbfF / max(y[IX["CBFAD"]], 5.0),
               CBFA=cbfA, CBFB=cbfB, OEFA=oefA, OEFB=oefB,
               CMROA=cmroA, CMROB=cmroB, dem=dem, paco2=paco2,
               SIGMA=sigma, DIL=dil, P_PERF=P_perf, LOAD=load,
               gv=gv, SIG_PVA=sig_v, DIL_PVA=dil_v, LOAD_PVA=load_v,
               Q_PVA=q_pva,
               E_CILO=e_cilo, E_NIF=e_nif, E_AZ=e_az, E_AHT=e_aht,
               CBF_CRIT=d["CBF_CRIT"], F_COLL=f_coll,
               CBF_WS=cbf_ws, OEF_WS=oef_ws, ANTE=ante, K_WS=k_ws,
               Q_MOYA=q_moya, Q_BYP=gb * (Pa - s["PA"]),
               Q_COLL=gc * (s["PB"] - s["PA"]), Q_ICA=gi * (Pa - s["PA"]),
               AR_POS=s["gA"] / gA_dil)
    if not full:
        return out

    # ---- reserve: two DIFFERENT numbers -------------------------------
    # (a) intrinsic: dilate A maximally, hold everything else where it is
    sa = solve_network(Pa, Pv, gSA, gSF, gb, gp, gc, gl,
                      gA_dil, s["gF"], s["gB"])
    out["CVR_INTR"] = (gA_dil * (sa[0] - Pv) / massA * 100.0 - cbfA) \
        / max(cbfA, 1e-9) * 100.0
    # (b) measured: an acetazolamide-like stimulus dilates EVERY arteriolar bed
    s2 = autoreg(d, Pa, gSA, gSF, gb, gp, gc, gl,
                 QtA * 1.45, QtF * 1.45, QtB * 1.45,
                 gA_con, gA_dil, gF_con, gF_dil, gB_con, gB_dil)
    out["CVR_MEAS"] = (s2["QA"] / massA * 100.0 - cbfA) / max(cbfA, 1e-9) * 100.0
    # ---- pressure-passivity: dCBF/dMAP by finite difference -----------
    sp = autoreg(d, Pa + 5.0, gSA, gSF, gb, gp, gc, gl, QtA, QtF, QtB,
                 gA_con, gA_dil, gF_con, gF_dil, gB_con, gB_dil)
    out["DCBF_DMAP"] = (sp["QA"] / massA * 100.0 - cbfA) / 5.0
    return out


# =========================================================================
# 3.  ODE SYSTEM
# =========================================================================
def rhs(t, y, d, sc):
    h = hemodynamics(d, y, t, sc, full=False)
    dy = np.zeros(NST)
    g = lambda n: y[IX[n]]

    # ---- drug PK -------------------------------------------------------
    dy[IX["ASA_G"]] = -d["ASA_KA"] * g("ASA_G")
    dy[IX["ASA_C"]] = d["ASA_KA"] * g("ASA_G") - d["ASA_KE"] * g("ASA_C")
    c_asa = g("ASA_C") / d["ASA_V"]
    dy[IX["PLTI"]] = d["K_PLT_ON"] * c_asa * (1 - g("PLTI")) \
        - d["K_PLT_OFF"] * g("PLTI")
    for nm, ka, ke in (("CILO", d["CILO_KA"], d["CILO_KE"]),
                       ("NIF", d["NIF_KA"], d["NIF_KE"]),
                       ("STAT", d["STAT_KA"], d["STAT_KE"]),
                       ("MIN", d["MIN_KA"], d["MIN_KE"]),
                       ("AHT", d["AHT_KA"], d["AHT_KE"])):
        dy[IX[nm + "_G"]] = -ka * g(nm + "_G")
        dy[IX[nm + "_C"]] = ka * g(nm + "_G") - ke * g(nm + "_C")
    dy[IX["AZ_C"]] = -d["AZ_KE"] * g("AZ_C")

    c_stat = g("STAT_C") / d["STAT_V"]
    e_stat = c_stat / (d["STAT_EC50"] + c_stat)
    c_min = g("MIN_C") / d["MIN_V"]
    e_min = c_min / (d["MIN_EC50"] + c_min)

    # ---- lesion: (RNF213 + secondary drive) -> PDGF -> intimal SMC ------
    dy[IX["PDGF"]] = d["K_PDGF"] * (d["RNF"] + d["SEC"]) * (1 + 0.6 * g("MMP9")) \
        - d["K_PDGF_D"] * g("PDGF")
    dy[IX["SMC"]] = d["K_SMC"] * (1 - d["K_STATIN"] * e_stat) * g("PDGF") \
        * (1 - g("SMC") / d["SMC_MAX"]) * sc.get("prog", 1.0)

    # ---- hypoxic drive -------------------------------------------------
    hypo = max(0.0, 1.0 - h["CBFA"] / d["CBF_A0"])
    dy[IX["HIF"]] = d["K_HIF_ON"] * hypo - d["K_HIF_OFF"] * g("HIF")
    dy[IX["VEGF"]] = d["K_VEGF"] * g("HIF") - d["K_VEGF_D"] * g("VEGF")
    dy[IX["MMP9"]] = d["K_MMP"] * g("VEGF") - d["K_MMP_D"] * g("MMP9")

    # ---- collaterals: RNF213 caps the ceiling it also creates ----------
    cap_m = max(d["GMOYA_CAP"] * (1 - d["RNF_CAP"] * d["RNF"]) * d["ANGIO"],
                d["GMOYA0"])
    cap_c = max(d["GCOLL_CAP"] * (1 - d["RNF_CAP"] * d["RNF"]) * d["ANGIO"]
                * (1 - d["PCA_COLL"] * d["PCA_INV"]), d["GCOLL0"])
    cap_v = max(d["GPVA_CAP"] * (1 - d["RNF_CAP"] * d["RNF"]) * d["ANGIO"]
                * (1 + d["PCA_PVA"] * d["PCA_INV"]), d["GPVA0"])
    drive = g("VEGF")
    # VEGF withdrawal prunes the bed: this is what makes bypass collapse the
    # fragile perforators instead of merely out-competing them.
    reg = max(0.0, 1.0 - drive / d["VEGF_REG"])
    dy[IX["GMOYA"]] = d["K_MOYA"] * drive * d["ANGIO"] * \
        max(0.0, 1 - g("GMOYA") / cap_m) \
        - d["K_MOYA_REG"] * reg * max(g("GMOYA") - d["GMOYA0"], 0.0)
    dy[IX["GCOLL"]] = d["K_COLL"] * drive * d["ANGIO"] * \
        max(0.0, 1 - g("GCOLL") / cap_c) \
        - d["K_COLL_REG"] * reg * max(g("GCOLL") - d["GCOLL0"], 0.0)
    # the periventricular route is recruited by the DEEP borderzone deficit
    # specifically -- which is what makes it a marker of watershed failure
    drive_v = drive * max(0.0, 1.0 - h["CBF_WS"] / d["CBF_A0"])
    dy[IX["GPVA"]] = d["K_PVA"] * drive_v * d["ANGIO"] * \
        max(0.0, 1 - g("GPVA") / cap_v) \
        - d["K_PVA_REG"] * reg * max(g("GPVA") - d["GPVA0"], 0.0)

    # ---- surgery -------------------------------------------------------
    kind = sc.get("surg_kind", None)
    ton = sc.get("surg_t", 1e9)
    scale = sc.get("surg_scale", 1.0)
    if kind is not None and t >= ton:
        if kind == "direct":
            dy[IX["GBYP"]] = (d["GBYP_DIR"] * scale - g("GBYP")) / d["TAU_BYP"]
        elif kind == "indirect":
            dy[IX["GBYP"]] = d["K_BYP_IND"] * g("VEGF") * d["ANGIO"] * \
                max(0.0, d["GBYP_IND"] * scale - g("GBYP"))
        elif kind == "combined":
            tgt = (d["GBYP_DIR"] + d["GBYP_IND"]) * scale
            dy[IX["GBYP"]] = \
                max(0.0, d["GBYP_DIR"] * scale - g("GBYP")) / d["TAU_BYP"] \
                + d["K_BYP_IND"] * g("VEGF") * d["ANGIO"] * \
                max(0.0, tgt - g("GBYP"))

    # ---- arteriolar remodelling: the price of chronic maximal dilation --
    on_floor = 1.0 / (1.0 + np.exp(-(h["AR_POS"] - 0.93) / 0.02))
    dy[IX["REMOD"]] = d["K_REM_ON"] * on_floor * (1 - g("REMOD")) \
        - d["K_REM_OFF"] * g("REMOD")

    # ---- ischaemic injury ----------------------------------------------
    sevA = max(0.0, (d["CBF_CRIT"] - h["CBF_WS"]) / d["CBF_CRIT"])
    sevB = max(0.0, (d["CBF_CRIT"] - h["CBFB"]) / d["CBF_CRIT"])
    # antiplatelet effect: aspirin (irreversible) + cilostazol (reversible)
    apl = min(g("PLTI") + d["CILO_APL"] *
              (h["E_CILO"] / max(d["CILO_EMAX"], 1e-9)), 1.0)
    emb = max(d["K_EMB"] * g("GMOYA") * (1 - d["K_EMB_ASA"] * apl), 0.0)
    dy[IX["INFA"]] = d["K_INF"] * max(0.0, d["WS_FRAC"] - g("INFA")) \
        * (sevA ** 2) + d["K_INF_EMB"] * emb * (1 - g("INFA"))
    dy[IX["INFB"]] = d["K_INF"] * (1 - g("INFB")) * (sevB ** 2)
    dy[IX["EMBH"]] = emb
    dy[IX["HYPOX"]] = hypo
    dy[IX["TIAB"]] = d["K_TIA"] * (1.0 if h["CBFA"] < d["CBF_CRIT"] else 0.0)
    # cumulative ISCHAEMIC hazard, per YEAR, integrated on a day clock.
    # Exponential in the margin to CBF_crit; the scale widens when the
    # arteriole is on its floor, because then MAP noise reaches the tissue.
    sd = d["CBF_SD0"] * (1 + d["K_PASS"] * on_floor)
    margin = (h["CBF_WS"] - d["CBF_CRIT"]) / sd
    dy[IX["ISCH"]] = (d["ISCH_HAZ0"] * np.exp(-np.clip(margin, -6, 25))
                      + d["ISCH_EMB"] * emb) / 365.0

    # ---- perforator microaneurysm & haemorrhage -------------------------
    stress = h["SIG_PVA"] * max(h["LOAD_PVA"], 0.0)
    dy[IX["ANEU"]] = d["K_ANEU"] * (stress ** 2) * (1 - g("ANEU")) \
        - d["K_ANEU_REP"] * (1 + d["K_MINO"] * e_min) * g("ANEU")
    haz = d["HEM_HAZ0"] * g("ANEU") * (h["SIG_PVA"] ** d["HEM_POW"])   # /yr
    dy[IX["HEMH"]] = haz / 365.0
    dy[IX["HEMV"]] = haz / 365.0 * d["HEM_VOL"]

    # ---- hyperperfusion -> BBB -> oedema -------------------------------
    # HYPERPERFUSION IS FOCAL AND RELATIVE.  It is read on the
    # peri-anastomotic cortex, against the flow that cortex had adapted to.
    dy[IX["CBFAD"]] = (h["CBFF"] - g("CBFAD")) / d["TAU_ADAPT"]
    hyper_rel = h["CBFF"] / max(g("CBFAD"), 5.0)
    over = max(0.0, hyper_rel - d["HYPER_THR"])
    dy[IX["BBB"]] = d["K_BBB"] * over - d["K_BBB_REP"] * g("BBB")
    dy[IX["EDEMA"]] = d["K_EDEMA"] * g("BBB") - d["K_EDEMA_R"] * g("EDEMA")

    # ---- cognition ------------------------------------------------------
    dy[IX["COG"]] = -d["K_COG"] * hypo - d["K_COG_INF"] * dy[IX["INFA"]]
    return dy


def y0(d):
    y = np.zeros(NST)
    y[IX["PDGF"]] = d["K_PDGF"] * d["RNF"] / d["K_PDGF_D"]
    y[IX["SMC"]] = 0.0
    y[IX["GMOYA"]] = d["GMOYA0"]
    y[IX["GPVA"]] = d["GPVA0"]
    y[IX["GCOLL"]] = d["GCOLL0"]
    y[IX["GBYP"]] = 0.0
    y[IX["ANEU"]] = 0.02
    y[IX["CBFAD"]] = d["CBF_A0"]      # a healthy barrier is adapted to normal
    return y


# =========================================================================
# 4.  DOSING / SIMULATION DRIVER
# =========================================================================


def simulate(scn, days, n=1500, pars=None):
    """
    scn keys:
      RNF, ANGIO, HB, MAP0, prog, map_mult   -- parameter overrides
      drugs: dict name -> (dose_mg, interval_days, start_day)
      surg_kind ('direct'|'indirect'|'combined'), surg_t, surg_scale
      paco2: callable t -> PaCO2
      az_bolus: (dose_mg, t)
    """
    p = dict(P)
    for k in ("RNF", "SEC", "ANGIO", "HB", "MAP0", "SAO2", "K_SMC", "G_LEAK",
              "F_FRAC", "REM_CON", "TAU_ADAPT",
              "GMOYA_CAP", "RNF_CAP", "STEN_MAX", "HEM_HAZ0", "K_INF",
              "CON_MAX", "K_REM_OFF", "GBYP_DIR", "GBYP_IND", "PCA_INV",
              "GPVA_CAP", "ISCH_HAZ0", "WS_K"):
        if k in scn:
            p[k] = scn[k]
    if pars:
        p.update(pars)
    d = init(p)
    sc = dict(scn)
    sc.setdefault("paco2", lambda t: d["PACO2_0"])
    sc.setdefault("prog", 1.0)

    drugs = scn.get("drugs", {})
    az = scn.get("az_bolus", None)

    # build dosing grid
    events = []
    for name, (dose, ii, t0) in drugs.items():
        tt = t0
        while tt <= days:
            events.append((tt, name, dose))
            tt += ii
    if az:
        events.append((az[1], "AZ", az[0]))
    events.sort()

    tgrid = np.unique(np.concatenate([
        np.linspace(0, days, n),
        np.array([e[0] for e in events]) if events else np.array([0.0]),
    ]))
    y = y0(d)
    T, Y = [0.0], [y.copy()]
    ei = 0
    # integrate segment by segment between doses
    bps = sorted(set([0.0] + [e[0] for e in events] + [days]))
    for a, b in zip(bps[:-1], bps[1:]):
        while ei < len(events) and abs(events[ei][0] - a) < 1e-12:
            _, nm, dose = events[ei]
            if nm == "AZ":
                y[IX["AZ_C"]] += dose
            else:
                y[IX[nm + "_G"]] += dose
            ei += 1
        if b <= a:
            continue
        sub = tgrid[(tgrid > a) & (tgrid <= b)]
        if sub.size == 0:
            sub = np.array([b])
        sol = solve_ivp(rhs, (a, b), y, args=(d, sc), t_eval=sub,
                        method="LSODA", rtol=1e-6, atol=1e-9, max_step=(b - a))
        if not sol.success:
            raise RuntimeError(f"integration failed {a}-{b}: {sol.message}")
        for k in range(sol.t.size):
            T.append(sol.t[k]); Y.append(sol.y[:, k].copy())
        y = sol.y[:, -1].copy()
    T = np.array(T); Y = np.array(Y)
    # post-hoc algebraic layer
    H = [hemodynamics(d, Y[i], T[i], sc) for i in range(T.size)]
    out = {k: np.array([hh[k] for hh in H]) for k in H[0]}
    out["t"] = T
    for nmm in SNAMES:
        out[nmm] = Y[:, IX[nmm]]
    out["_d"] = d
    return out


# =========================================================================
# 5.  ARCHETYPES
# =========================================================================
#  ANGIO = angiogenic capacity (child ~1.0, adult ~0.35); K_SMC sets how fast
#  the intimal lesion accumulates.  These are the ONLY structural differences
#  between the archetypes -- everything else in the report is a consequence.
PED     = dict(RNF=1.0, ANGIO=1.00, K_SMC=5.5e-4)       # child, fast lesion
ADULT_I = dict(RNF=1.0, ANGIO=0.58, K_SMC=2.6e-4)       # adult ischaemic onset
ADULT_H = dict(RNF=1.0, ANGIO=0.86, K_SMC=2.3e-4,
               PCA_INV=1.0)                             # adult haemorrhagic onset
SCD     = dict(RNF=0.0, SEC=1.30, ANGIO=0.85, HB=8.0,
               K_SMC=3.4e-4)                            # sickle-cell moyamoya
ASYMP   = dict(RNF=1.0, ANGIO=0.66, K_SMC=1.05e-4)      # slow / asymptomatic


def at(o, day, key):
    return float(np.interp(day, o["t"], o[key]))


# =========================================================================
# 6.  REPORT
# =========================================================================
def hr(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)


def state_at(o, day):
    """Freeze the full state vector at a given day, for challenge probes."""
    return np.array([float(np.interp(day, o["t"], o[k])) for k in SNAMES])


def probe(o, day, pa_mult=1.0, paco2=40.0, az=0.0, hold=0.02):
    """Apply an instantaneous challenge to a frozen state and re-read it."""
    dd = o["_d"]
    ys = state_at(o, day)
    ys[IX["AZ_C"]] += az
    sc = dict(paco2=lambda t: paco2, prog=0.0, map_mult=pa_mult)
    if az == 0.0:
        return hemodynamics(dd, ys, 0.0, sc)
    sol = solve_ivp(rhs, (0, hold), ys, args=(dd, sc),
                    t_eval=np.linspace(0, hold, 40), method="LSODA",
                    rtol=1e-7, atol=1e-10)
    H = [hemodynamics(dd, sol.y[:, k], sol.t[k], sc) for k in range(sol.t.size)]
    j = int(np.argmax([abs(hh["CBFA"] - H[0]["CBFA"]) for hh in H]))
    return H[j]


def main():
    d0 = init(P)

    # ==================================================================
    hr("0.  DERIVED GEOMETRY, AND A THRESHOLD THAT IS NOT A PARAMETER")
    print(f"  CaO2 = 1.34*Hb*SaO2/100        = {d0['CAO2']:.4f} mL O2/mL blood")
    print(f"  CBF_crit = CMRO2/(CaO2*OEF_max) = {d0['CBF_CRIT']:.2f} mL/100g/min")
    print("  <-- the textbook penumbral 20 mL/100g/min is DERIVED, not fitted.")
    print(f"  R_ica0 {d0['R_ICA0']:.4f}   R_artA0 {d0['R_ARTA0']:.4f}   "
          f"gA range [{d0['GA_CON']:.3f}, {d0['GA_DIL']:.3f}] mL/min/mmHg")
    print(f"  Healthy operating point: P_A {d0['P_A0']:.1f} mmHg, "
          f"CBF 50.0, OEF {P['OEF_BASE']:.3f}, CMRO2 {P['CMRO2_0']:.2f}")
    print("\n  Anaemia moves the threshold; no vasodilator can:")
    print(f"  {'Hb (g/dL)':>10} {'CaO2':>8} {'CBF_crit':>10} {'vs Hb 15':>10}")
    for hb in (15, 13, 11, 9, 8, 7):
        dd = init({**P, "HB": hb})
        print(f"  {hb:>10} {dd['CAO2']:>8.4f} {dd['CBF_CRIT']:>10.2f} "
              f"{dd['CBF_CRIT']/d0['CBF_CRIT']:>9.2f}x")

    # ==================================================================
    hr("1.  A CRITICAL INLET CONDUCTANCE EXISTS  (the reserve cliff)")
    print("  Sweep the total inlet conductance gS = g_ica + g_moya + g_pva +")
    print("  g_bypass.  Report where the arteriole sits in its own range, the")
    print("  two reserve numbers, and pressure-passivity dCBF/dMAP.\n")
    print(f"  {'gS':>7} {'STEN_eq':>8} {'P_A':>7} {'CBF_A':>7} {'CBF_ws':>7}"
          f" {'gA/gA_dil':>10} {'CVR_int%':>9} {'CVR_meas%':>10} {'dCBF/dMAP':>10}")
    d = init(P)
    crit = None
    for gS in (6.25, 4.0, 3.0, 2.2, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4):
        y = y0(d)
        gi_t = max(gS - d["GMOYA0"] - d["GPVA0"], 1e-6)
        sten = float(np.clip(1 - (gi_t / d["GI0"]) ** 0.25, 0.0, d["STEN_MAX"]))
        y[IX["SMC"]] = sten / d["STEN_MAX"]
        h = hemodynamics(d, y, 0.0, dict(paco2=lambda t: 40.0))
        if crit is None and h["AR_POS"] > 0.999:
            crit = gS
        print(f"  {h['gS']:>7.3f} {h['STEN']:>8.3f} {h['PA']:>7.1f} "
              f"{h['CBFA']:>7.1f} {h['CBF_WS']:>7.1f} {h['AR_POS']:>10.3f} "
              f"{h['CVR_INTR']:>9.1f} {h['CVR_MEAS']:>10.1f} "
              f"{h['DCBF_DMAP']:>10.3f}")
    # locate it exactly
    def arpos(gS):
        y = y0(d)
        gi_t = max(gS - d["GMOYA0"] - d["GPVA0"], 1e-9)
        y[IX["SMC"]] = float(np.clip(1 - (gi_t / d["GI0"]) ** 0.25, 0, 1)) \
            / d["STEN_MAX"]
        return hemodynamics(d, y, 0.0, dict(paco2=lambda t: 40.0),
                            full=False)["AR_POS"] - 0.9995
    gS_crit = brentq(arpos, 1.0, 4.0, xtol=1e-6)
    print(f"\n  CRITICAL INLET CONDUCTANCE gS* = {gS_crit:.4f} mL/min/mmHg.")
    print("  Above it CBF is 50.0 whatever the pressure does and dCBF/dMAP = 0;")
    print("  below it the arteriole is on its floor, CBF becomes a linear")
    print("  function of MAP, and the measured CO2/acetazolamide response")
    print("  turns NEGATIVE (steal).  All three switch at the same point.")
    print("  Note gS* is a property of the NETWORK: the equivalent isolated")
    print(f"  ICA stenosis is {1-(gS_crit/d['GI0'])**0.25:.3f} diameter, but a")
    print("  patient with good collaterals reaches the same gS at a far worse")
    print("  angiographic stenosis.  Suzuki grade cannot locate this point.")

    # ==================================================================
    hr("2.  NATURAL HISTORY OF FIVE ARCHETYPES (untreated, 10 y)")
    print("  Only ANGIO (angiogenic capacity), K_SMC (lesion growth rate),")
    print("  PCA_INV (posterior-circulation involvement), Hb and the")
    print("  RNF213/secondary drive differ between these five patients.\n")
    runs = {}
    for nm, sc in (("paediatric", PED), ("adult-ischaemic", ADULT_I),
                   ("adult-haemorrhagic", ADULT_H), ("sickle-cell", SCD),
                   ("slow/asymptomatic", ASYMP)):
        runs[nm] = simulate(dict(sc), 3650, n=1100)
    print(f"  {'archetype':>19} {'STEN':>5} {'gS':>5} {'g_pva':>6} {'CBF_A':>6}"
          f" {'CBF_ws':>7} {'OEF':>5} {'CVRi%':>6} {'CVRm%':>6} {'INF%':>5}"
          f" {'P_isch':>7} {'P_haem':>7} {'COGz':>6}")
    for nm, o in runs.items():
        i = -1
        print(f"  {nm:>19} {o['STEN'][i]:>5.3f} {o['gS'][i]:>5.2f} "
              f"{o['gv'][i]:>6.3f} {o['CBFA'][i]:>6.1f} {o['CBF_WS'][i]:>7.1f} "
              f"{o['OEFA'][i]:>5.3f} {o['CVR_INTR'][i]:>6.1f} "
              f"{o['CVR_MEAS'][i]:>6.1f} {o['INFA'][i]*100:>5.1f} "
              f"{1-np.exp(-o['ISCH'][i]):>7.3f} {1-np.exp(-o['HEMH'][i]):>7.3f} "
              f"{o['COG'][i]:>6.2f}")
    print("\n  THE SILENT PHASE IS A CONSEQUENCE OF R ~ (1-STEN)^-4.")
    o = runs["adult-ischaemic"]
    print(f"  {'year':>5} {'STEN':>6} {'g_ica':>7} {'gS':>6} {'CBF_A':>7}"
          f" {'CBF_ws':>7} {'OEF':>6} {'CVRi%':>7} {'INF%':>6}")
    for yr in (1, 2, 3, 4, 5, 6, 8, 10):
        f = lambda k: at(o, yr * 365, k)
        print(f"  {yr:>5} {f('STEN'):>6.3f} {f('gi'):>7.3f} {f('gS'):>6.2f} "
              f"{f('CBFA'):>7.1f} {f('CBF_WS'):>7.1f} {f('OEFA'):>6.3f} "
              f"{f('CVR_INTR'):>7.1f} {f('INFA')*100:>6.2f}")
    print("  Stenosis is linear in time; conductance is its fourth power.  CBF")
    print("  is held at 50 for years while reserve is silently spent, then the")
    print("  patient decompensates over months.  Reserve is lost at year 3-4;")
    print("  flow does not move until year 5.  That gap is the disease.")
    print("\n  AND THE STEADY STATE IS FLOW PINNED AT THE THRESHOLD:")
    for yr in (5, 6, 7, 8, 10):
        f = lambda k: at(o, yr * 365, k)
        print(f"    y{yr:>2}: CBF_ws {f('CBF_WS'):5.1f}  (CBF_crit "
              f"{d0['CBF_CRIT']:.1f})  infarct {f('INFA')*100:5.2f}%  "
              f"demand removed {f('INFA')*100:5.2f}%")
    print("  Infarction REMOVES demand, which raises the flow per surviving")
    print("  gram back to threshold.  A chronic moyamoya SPECT therefore looks")
    print("  only mildly abnormal precisely BECAUSE tissue was lost.  Reading")
    print("  CBF alone systematically under-rates the disease.")

    # ==================================================================
    hr("3.  ACETAZOLAMIDE AND CO2 ARE NOT THE SAME PROBE")
    print("  Acetazolamide acts on the arteriole (tissue carbonic anhydrase).")
    print("  PaCO2 acts on the arteriole AND on the pial collateral conduits.")
    print("  In a normal brain that distinction is invisible.  In moyamoya the")
    print("  conduits ARE the circulation, so the two probes disagree.\n")
    o = runs["adult-ischaemic"]
    print(f"  {'state':>16} {'CVRi%':>6} {'CBF_A':>6} | {'ACZ 1 g':>9}"
          f" {'dA%':>6} {'dB%':>6} | {'PaCO2 50':>9} {'dA%':>6}"
          f" | {'PaCO2 25':>9} {'dA%':>6}")
    for lbl, day in (("normal (y1)", 365), ("reserve lost (y4)", 1460),
                     ("decompensated (y6)", 2190), ("end-stage (y10)", 3650)):
        b = probe(o, day)
        pz = probe(o, day, az=1000.0)
        ph = probe(o, day, paco2=50.0)
        pl = probe(o, day, paco2=25.0)
        print(f"  {lbl:>16} {b['CVR_INTR']:>6.1f} {b['CBFA']:>6.1f} | "
              f"{pz['CBFA']:>9.1f} {(pz['CBFA']/b['CBFA']-1)*100:>6.1f} "
              f"{(pz['CBFB']/b['CBFB']-1)*100:>6.1f} | "
              f"{ph['CBFA']:>9.1f} {(ph['CBFA']/b['CBFA']-1)*100:>6.1f} | "
              f"{pl['CBFA']:>9.1f} {(pl['CBFA']/b['CBFA']-1)*100:>6.1f}")
    print("\n  The donor territory answers acetazolamide every time.  The")
    print("  affected territory stops answering, then answers with the WRONG")
    print("  SIGN -- because territory B is territory A's supply, so dilating B")
    print("  steals from A.  The steal is a property of the network, not an")
    print("  extra assumption.  Hypercapnia, by contrast, opens the conduits")
    print("  themselves and helps -- which is why normocapnia-to-mild-")
    print("  hypercapnia is the anaesthetic target, and why an acetazolamide")
    print("  study and a breath-hold study are not interchangeable.")

    # ==================================================================
    hr("4.  THE CRYING CHILD: ONE PaCO2 STEP, OPPOSITE OUTCOMES")
    print("  Hyperventilation to PaCO2 25 mmHg.  Two things happen at once:")
    print("  demand falls (protective) and the collateral conduits constrict")
    print("  (harmful).  Which dominates depends only on whether the arteriole")
    print("  has any range left to answer with.\n")
    print(f"  {'patient':>22} {'CBF_ws 40':>10} {'CBF_ws 25':>10} {'ratio':>7}"
          f" {'CBF_crit':>9} {'crosses?':>9}")
    for lbl, key, day in (("normal", None, 0), ("preserved reserve", "slow/asymptomatic", 1095),
                          ("adult, reserve lost", "adult-ischaemic", 1460),
                          ("adult, decompensated", "adult-ischaemic", 2190),
                          ("child, decompensated", "paediatric", 1460)):
        if key is None:
            oo = simulate(dict(K_SMC=0.0), 30, n=40); day = 30
        else:
            oo = runs[key]
        b = probe(oo, day); l = probe(oo, day, paco2=25.0)
        print(f"  {lbl:>22} {b['CBF_WS']:>10.1f} {l['CBF_WS']:>10.1f} "
              f"{l['CBF_WS']/b['CBF_WS']:>7.3f} {b['CBF_CRIT']:>9.1f} "
              f"{'YES' if l['CBF_WS'] < b['CBF_CRIT'] <= b['CBF_WS'] else ('already' if l['CBF_WS'] < b['CBF_CRIT'] else 'no'):>9}")
    print("\n  In the normal brain a 15 mmHg fall in PaCO2 costs 44% of the")
    print("  flow and is entirely safe, because the flow was never the")
    print("  constraint.  In the decompensated hemisphere it costs less in")
    print("  percentage terms and is an infarct, because there was no margin.")

    # ==================================================================
    hr("5.  THE BLOOD-PRESSURE PARADOX HAS A COMPUTABLE OPTIMUM")
    print("  Freeze the disease state and sweep MAP.  Ischaemic hazard runs")
    print("  through the pressure-passive term (dCBF/dMAP > 0), haemorrhagic")
    print("  hazard through periventricular wall stress sigma ~ P_perf*d^1.5.")
    print("  These have opposite sign in MAP, so an optimum exists.\n")
    for lbl, key, day in (("adult ISCHAEMIC (PCA spared)", "adult-ischaemic", 2555),
                          ("adult HAEMORRHAGIC (PCA involved)", "adult-haemorrhagic", 3285)):
        oo = runs[key]; dd = oo["_d"]; ys = state_at(oo, day)
        print(f"  --- {lbl}, state at day {day} ---")
        print(f"  {'MAP':>5} {'CBF_A':>7} {'CBF_ws':>7} {'sig_pva':>8}"
              f" {'isch/yr':>9} {'haem/yr':>9} {'TOTAL/yr':>9}")
        best = (None, 1e9); rows = []
        for mp in range(65, 126, 5):
            h = hemodynamics(dd, ys, 0.0,
                             dict(paco2=lambda t: 40.0, map_mult=mp / dd["MAP0"]))
            sd = dd["CBF_SD0"] * (1 + dd["K_PASS"] *
                                  (1.0 if h["AR_POS"] > 0.99 else 0.0))
            ih = dd["ISCH_HAZ0"] * np.exp(
                -np.clip((h["CBF_WS"] - dd["CBF_CRIT"]) / sd, -6, 25)) \
                + dd["ISCH_EMB"] * dd["K_EMB"] * ys[IX["GMOYA"]]
            hh = dd["HEM_HAZ0"] * ys[IX["ANEU"]] * h["SIG_PVA"] ** dd["HEM_POW"]
            rows.append((mp, ih, hh, ih + hh))
            if ih + hh < best[1]:
                best = (mp, ih + hh)
            print(f"  {mp:>5} {h['CBFA']:>7.1f} {h['CBF_WS']:>7.1f} "
                  f"{h['SIG_PVA']:>8.3f} {ih:>9.4f} {hh:>9.4f} {ih+hh:>9.4f}")
        print(f"  --> hazard-minimising MAP = {best[0]} mmHg "
              f"(total {best[1]:.4f}/yr)")
        lo = [r for r in rows if r[0] == 70][0]; hi = [r for r in rows if r[0] == 110][0]
        print(f"      MAP 70 : isch {lo[1]:.4f} haem {lo[2]:.4f} total {lo[3]:.4f}")
        print(f"      MAP 110: isch {hi[1]:.4f} haem {hi[2]:.4f} total "
              f"{hi[3]:.4f}\n")
    print("  The two phenotypes have DIFFERENT optima from the SAME equations.")
    print("  This is why 'control the blood pressure' is not a moyamoya")
    print("  instruction until you say which hemisphere you mean.")

    # ==================================================================
    hr("6.  BYPASS: ONE MECHANISM, BOTH PHENOTYPES  (JAM-trial emulation)")
    print("  Direct STA-MCA bypass vs conservative management, operated once")
    print("  the phenotype is established, 5 y of follow-up.\n")
    print(f"  {'arm':>38} {'CBF_A':>6} {'CBF_ws':>7} {'P_A':>6} {'Q_pva':>7}"
          f" {'g_pva':>6} {'sig_pva':>8} {'5y haem':>8} {'5y isch':>8}")
    tab = {}
    for phen, sc_, ton in (("haemorrhagic", ADULT_H, 2555.0),
                           ("ischaemic", ADULT_I, 1825.0)):
        for arm, surg in (("conservative", None), ("direct bypass", "direct")):
            s_ = dict(sc_)
            if surg:
                s_.update(surg_kind=surg, surg_t=ton)
            oo = simulate(s_, ton + 5 * 365, n=1500)
            i0 = int(np.argmin(np.abs(oo["t"] - ton)))
            reb = oo["HEMH"][-1] - oo["HEMH"][i0]
            isc = oo["ISCH"][-1] - oo["ISCH"][i0]
            tab[(phen, arm)] = (reb, isc)
            print(f"  {phen+' / '+arm:>38} {oo['CBFA'][-1]:>6.1f} "
                  f"{oo['CBF_WS'][-1]:>7.1f} {oo['PA'][-1]:>6.1f} "
                  f"{oo['Q_PVA'][-1]:>7.1f} {oo['gv'][-1]:>6.3f} "
                  f"{oo['SIG_PVA'][-1]:>8.3f} {1-np.exp(-reb):>8.3f} "
                  f"{1-np.exp(-isc):>8.3f}")
    print()
    for phen in ("haemorrhagic", "ischaemic"):
        r0, i0_ = tab[(phen, "conservative")]
        r1, i1 = tab[(phen, "direct bypass")]
        print(f"  {phen:>14}:  haemorrhage HR {r1/max(r0,1e-12):.3f}"
              f"   ischaemic HR {i1/max(i0_,1e-12):.3f}")
    print("\n  JAM Trial (Miyamoto 2014, PMID 24668203): 5-y rebleeding 11.9%")
    print("  with bypass vs 31.6% conservative, rebleed HR 0.355.  Only the")
    print("  CONSERVATIVE rate was used to set HEM_HAZ0, so the HR is a")
    print("  PREDICTION -- and the model UNDER-predicts the benefit by roughly")
    print("  a third (0.49 against 0.355).  Reported as a miss, not a fit: the")
    print("  most likely reason is that the model routes the whole benefit")
    print("  through pruning of the periventricular bed, whereas a real bypass")
    print("  probably also lowers the transmural pressure those vessels see")
    print("  faster than the bed can regress.")
    print("  Its mechanism in the model is visible in the Q_pva")
    print("  column: raising P_A collapses the gradient across the")
    print("  periventricular anastomosis, VEGF falls, the bed is pruned, and")
    print("  the wall stress that was going to rupture goes with it.  The")
    print("  operation does not reinforce the fragile vessel -- it retires it.")

    # ==================================================================
    hr("7.  HYPERPERFUSION IS FOCAL, RELATIVE, AND PRICED IN PRE-OP ISCHAEMIA")
    print("  A LUMPED TERRITORY CANNOT HYPERPERFUSE: even a patent graft leaves")
    print("  the territorial arteriole still dilating, so mean territorial CBF")
    print("  never exceeds demand.  The syndrome therefore has to be FOCAL --")
    print("  the cortex the graft is sewn to, which sees P_F >> P_A -- and it")
    print("  has to be RELATIVE, measured against the flow that cortex had")
    print("  adapted to.  50 mL/100g/min is normal, and it injures a barrier")
    print("  that has lived on 18 for five years.\n")
    print(f"  {'pre-op hemisphere':>22} {'REMOD':>6} {'CVRi%':>6} {'CBF_F pre':>9}"
          f" {'CBF_F pk':>9} {'rel pk':>7} {'pk day':>7} {'d>1.35':>7}"
          f" {'oedema':>7} {'CBFws +':>8}")
    for lbl, sc_, ton in (("preserved reserve", ASYMP, 1095.0),
                          ("reserve just lost", ADULT_I, 1460.0),
                          ("decompensated 2 y", ADULT_I, 2190.0),
                          ("decompensated 5 y", ADULT_I, 3285.0),
                          ("paediatric, severe", PED, 1825.0)):
        s_ = dict(sc_); s_.update(surg_kind="direct", surg_t=ton)
        oo = simulate(s_, ton + 120, n=2800)
        m = oo["t"] >= ton
        rel = oo["HYPER_REL"][m]; ip = int(np.argmax(rel))
        dt = np.gradient(oo["t"][m]); pre = at(oo, ton - 1, "CBF_WS")
        print(f"  {lbl:>22} {at(oo,ton-1,'REMOD'):>6.3f} "
              f"{at(oo,ton-1,'CVR_INTR'):>6.1f} {at(oo,ton-1,'CBFF'):>9.1f} "
              f"{oo['CBFF'][m].max():>9.1f} {rel[ip]:>7.2f} "
              f"{oo['t'][m][ip]-ton:>7.2f} "
              f"{float(np.sum(dt[rel>P['HYPER_THR']])):>7.1f} "
              f"{oo['EDEMA'][m].max():>7.3f} {oo['CBF_WS'][m][-1]-pre:>8.1f}")
    print("\n  A hemisphere with reserve does not hyperperfuse at all: its")
    print("  barrier was already adapted to 50.  The syndrome belongs to the")
    print("  arteriole's HISTORY, not to the surgeon.  REMOD (vasoparalysis)")
    print("  sets the HEIGHT of the surge; the barrier's re-adaptation time")
    print("  constant sets its DURATION.")
    print("\n  AND THE INCIDENCE IS SET BY SOMETHING THE SURGEON CANNOT SEE.")
    print("  Sweep the pial coupling g_leak between the peri-anastomotic")
    print("  cortex and the rest of the territory.  It trades the TERRITORIAL")
    print("  benefit against the FOCAL surge, in opposite directions:\n")
    print(f"  {'g_leak':>7} {'Q_byp':>7} {'P_A':>6} {'P_F':>6} {'P_F-P_A':>8}"
          f" {'CBFws +':>8} {'rel surge':>10} {'d>1.35':>7} {'syndrome?':>10}")
    for gl in (0.25, 0.35, 0.55, 0.80, 1.20, 1.80):
        s_ = dict(ADULT_I)
        s_.update(surg_kind="direct", surg_t=3285.0, G_LEAK=gl)
        oo = simulate(s_, 3285 + 120, n=2400)
        m = oo["t"] >= 3285
        rel = oo["HYPER_REL"][m]; dt = np.gradient(oo["t"][m])
        pre = at(oo, 3284, "CBF_WS")
        dhi = float(np.sum(dt[rel > P["HYPER_THR"]]))
        print(f"  {gl:>7.2f} {oo['Q_BYP'][m][-1]:>7.1f} {oo['PA'][m][-1]:>6.1f} "
              f"{oo['PF'][m][-1]:>6.1f} "
              f"{oo['PF'][m][-1]-oo['PA'][m][-1]:>8.1f} "
              f"{oo['CBF_WS'][m][-1]-pre:>8.1f} {rel.max():>10.2f} {dhi:>7.1f} "
              f"{('YES' if dhi > 3 else 'no'):>10}")
    print("\n  A WELL-COUPLED cortex spreads the graft flow over the territory:")
    print("  large territorial gain, no surge.  A POORLY-COUPLED recipient")
    print("  keeps the pressure locally: small territorial gain, a big focal")
    print("  surge.  Since g_leak is a property of the patient's pial network")
    print("  and not of the operation, this is the model's account of why")
    print("  hyperperfusion occurs in a MINORITY of technically perfect")
    print("  bypasses, and why its occurrence is not a measure of surgical")
    print("  quality.  It also predicts the inverse correlation: the")
    print("  hemispheres that hyperperfuse are the ones that gained LEAST.")
    print("\n  Post-operative MAP policy on the decompensated-5y hemisphere:")
    print(f"    {'policy':>16} {'rel surge':>10} {'d>1.35':>8} {'oedema':>7}"
          f" {'min CBF_ws':>11}")
    for mult, lbl in ((1.00, "standard"), (0.90, "MAP -10%"),
                      (0.85, "MAP -15%"), (0.75, "MAP -25%")):
        s_ = dict(ADULT_I)
        s_.update(surg_kind="direct", surg_t=3285.0, map_mult=mult,
                  map_mult_t=3285.0, G_LEAK=0.35)
        oo = simulate(s_, 3285 + 90, n=2400)
        m = oo["t"] >= 3285
        rel = oo["HYPER_REL"][m]; dt = np.gradient(oo["t"][m])
        print(f"    {lbl:>16} {rel.max():>10.2f} "
              f"{float(np.sum(dt[rel>P['HYPER_THR']])):>8.1f} "
              f"{oo['EDEMA'][m].max():>7.3f} {oo['CBF_WS'][m].min():>11.1f}")
    print("  Note the sign flip against section 5: the SAME hemisphere that")
    print("  wanted a HIGHER MAP last week wants a lower one now.  Nothing")
    print("  about the patient changed except which side of gS* they are on.")

    hr("8.  INDIRECT BYPASS IS A BET ON ANGIOGENIC CAPACITY")
    print("  Identical operations offered to a child and an adult.  The only")
    print("  difference in the model is ANGIO.\n")
    print(f"  {'patient':>10} {'ANGIO':>6} {'operation':>10} {'g_byp d90':>10}"
          f" {'g_byp d365':>11} {'CBF_ws d365':>12} {'INF% d1095':>11}"
          f" {'P_isch 3y':>10}")
    for lbl, sc_ in (("child", PED), ("adult", ADULT_I)):
        for op in (None, "indirect", "direct", "combined"):
            s_ = dict(sc_); ton = 1460.0
            if op:
                s_.update(surg_kind=op, surg_t=ton)
            oo = simulate(s_, ton + 1095, n=1400)
            i0 = int(np.argmin(np.abs(oo["t"] - ton)))
            print(f"  {lbl:>10} {oo['_d']['ANGIO']:>6.2f} {(op or 'none'):>10} "
                  f"{at(oo,ton+90,'GBYP'):>10.3f} {at(oo,ton+365,'GBYP'):>11.3f} "
                  f"{at(oo,ton+365,'CBF_WS'):>12.1f} "
                  f"{at(oo,ton+1095,'INFA')*100:>11.2f} "
                  f"{1-np.exp(-(oo['ISCH'][-1]-oo['ISCH'][i0])):>10.3f}")
    print("\n  The indirect construct needs the patient to grow it, and the")
    print("  patient it is growing in is the same one whose RNF213 caps")
    print("  angiogenesis.  In the child it approaches the direct graft; in")
    print("  the adult it is most of a year late and short of target.")

    # ==================================================================
    hr("9.  MEDICAL THERAPY MOVES THE HAZARD TERMS, NOT THE GEOMETRY")
    print("  5 y from the adult-ischaemic archetype, surgery-naive.\n")
    print(f"  {'arm':>26} {'STEN':>6} {'gS':>5} {'CBF_A':>6} {'CBF_ws':>7}"
          f" {'CVRi%':>6} {'INF%':>5} {'5y isch':>8} {'5y haem':>8} {'COGz':>6}")
    arms = [
        ("none", {}),
        ("aspirin 100 mg/d", dict(drugs={"ASA": (100.0, 1.0, 0.0)})),
        ("cilostazol 200 mg/d", dict(drugs={"CILO": (100.0, 0.5, 0.0)})),
        ("atorvastatin 20 mg/d", dict(drugs={"STAT": (20.0, 1.0, 0.0)})),
        ("nifedipine GITS 60 mg/d", dict(drugs={"NIF": (60.0, 1.0, 0.0)})),
        ("minocycline 200 mg/d", dict(drugs={"MIN": (100.0, 0.5, 0.0)})),
        ("ASA + statin", dict(drugs={"ASA": (100.0, 1.0, 0.0),
                                     "STAT": (20.0, 1.0, 0.0)})),
        ("antihypertensive alone", dict(drugs={"AHT": (10.0, 1.0, 0.0)})),
        ("direct bypass (reference)", dict(surg_kind="direct", surg_t=30.0)),
    ]
    for lbl, extra in arms:
        s_ = dict(ADULT_I); s_.update(extra)
        oo = simulate(s_, 5 * 365, n=1300)
        print(f"  {lbl:>26} {oo['STEN'][-1]:>6.3f} {oo['gS'][-1]:>5.2f} "
              f"{oo['CBFA'][-1]:>6.1f} {oo['CBF_WS'][-1]:>7.1f} "
              f"{oo['CVR_INTR'][-1]:>6.1f} {oo['INFA'][-1]*100:>5.2f} "
              f"{1-np.exp(-oo['ISCH'][-1]):>8.3f} "
              f"{1-np.exp(-oo['HEMH'][-1]):>8.3f} {oo['COG'][-1]:>6.2f}")
    print("\n  Every medical arm acts on a HAZARD (embolism, wall repair,")
    print("  lesion growth rate).  None of them changes gS by more than a few")
    print("  percent, because gS is set by lumen geometry and collateral")
    print("  anatomy.  Nifedipine is the instructive failure: a real cerebral")
    print("  vasodilator that cannot dilate a vessel already on its floor,")
    print("  that dilates the DONOR territory instead, and that lowers the MAP")
    print("  a pressure-passive territory is living on.")

    # ==================================================================
    hr("10.  SICKLE-CELL MOYAMOYA: THE THRESHOLD MOVES, NOT THE FLOW")
    print(f"  {'arm':>26} {'Hb':>4} {'CBF_crit':>9} {'CBF_A':>6} {'CBF_ws':>7}"
          f" {'OEF':>6} {'INF% 5y':>8} {'5y isch':>8}")
    for lbl, hb in (("untransfused", 8.0), ("transfusion to 9", 9.0),
                    ("chronic transfusion 11", 11.0),
                    ("aggressive to 13", 13.0), ("normal-Hb comparator", 15.0)):
        s_ = dict(SCD); s_["HB"] = hb
        oo = simulate(s_, 5 * 365, n=900)
        print(f"  {lbl:>26} {hb:>4.0f} {oo['_d']['CBF_CRIT']:>9.1f} "
              f"{oo['CBFA'][-1]:>6.1f} {oo['CBF_WS'][-1]:>7.1f} "
              f"{oo['OEFA'][-1]:>6.3f} {oo['INFA'][-1]*100:>8.2f} "
              f"{1-np.exp(-oo['ISCH'][-1]):>8.3f}")
    print("\n  The vascular lesion is IDENTICAL in all five rows.  Only the")
    print("  oxygen-carrying capacity differs, and it moves the threshold the")
    print("  flow is being compared against.  This is the model's account of")
    print("  why transfusion -- which does nothing to any vessel -- is the")
    print("  intervention that works in sickle-cell moyamoya (STOP/SIT).")

    # ==================================================================
    hr("11.  PCA INVOLVEMENT CONVERTS THE PHENOTYPE")
    print("  One flag, two opposite-signed consequences: it removes the SAFE")
    print("  leptomeningeal donor and forces the DANGEROUS periventricular")
    print("  route.  Same patient otherwise, 10 y.\n")
    print(f"  {'PCA involvement':>18} {'g_coll':>7} {'g_pva':>6} {'CBF_ws':>7}"
          f" {'sig_pva':>8} {'ANEU':>6} {'INF%':>6} {'10y isch':>9}"
          f" {'10y haem':>9}")
    for pv in (0.0, 0.25, 0.5, 0.75, 1.0):
        s_ = dict(ADULT_I); s_["PCA_INV"] = pv
        oo = simulate(s_, 3650, n=1000)
        print(f"  {pv:>18.2f} {oo['gc'][-1]:>7.3f} {oo['gv'][-1]:>6.3f} "
              f"{oo['CBF_WS'][-1]:>7.1f} {oo['SIG_PVA'][-1]:>8.3f} "
              f"{oo['ANEU'][-1]:>6.3f} {oo['INFA'][-1]*100:>6.2f} "
              f"{1-np.exp(-oo['ISCH'][-1]):>9.3f} "
              f"{1-np.exp(-oo['HEMH'][-1]):>9.3f}")
    print("\n  Haemorrhage risk is not a separate disease with a separate")
    print("  cause.  It is what this circulation does when the collateral it")
    print("  is allowed to use is the fragile one.")

    # ==================================================================
    hr("12.  RNF213 ENTERS TWICE, WITH OPPOSITE SIGN")
    print(f"  {'genotype':>24} {'RNF':>4} {'STEN 10y':>9} {'moya cap':>9}"
          f" {'g_moya':>7} {'gS':>6} {'CBF_ws':>7} {'INF%':>6} {'10y isch':>9}")
    for lbl, rnf, sec, ksmc in (("wild type / quasi-MMD", 0.0, 1.0, 2.6e-4),
                                ("R4810K heterozygote", 1.0, 0.0, 2.6e-4),
                                ("R4810K homozygote", 1.6, 0.0, 4.2e-4)):
        s_ = dict(ADULT_I); s_.update(RNF=rnf, SEC=sec, K_SMC=ksmc)
        oo = simulate(s_, 3650, n=1000)
        cap = P["GMOYA_CAP"] * (1 - P["RNF_CAP"] * rnf) * oo["_d"]["ANGIO"]
        print(f"  {lbl:>24} {rnf:>4.1f} {oo['STEN'][-1]:>9.3f} {cap:>9.3f} "
              f"{oo['GMOYA'][-1]:>7.3f} {oo['gS'][-1]:>6.2f} "
              f"{oo['CBF_WS'][-1]:>7.1f} {oo['INFA'][-1]*100:>6.2f} "
              f"{1-np.exp(-oo['ISCH'][-1]):>9.3f}")
    print("\n  The wild-type row has the SAME lesion growth rate and a")
    print("  collateral ceiling 1.7x higher, and that alone separates the")
    print("  outcomes.  A gene that builds the obstruction and then forbids")
    print("  the detour is worse than either effect alone -- and it is why")
    print("  RNF213 status is prognostic beyond the angiogram.")

    # ==================================================================
    hr("13.  VIRTUAL POPULATION (n=400, 5 y): WHAT ACTUALLY PREDICTS?")
    rng = np.random.default_rng(20260806)
    N = 400; rows = []
    for i in range(N):
        s_ = dict(
            RNF=float(rng.choice([0.0, 1.0, 1.6], p=[0.20, 0.72, 0.08])),
            ANGIO=float(np.clip(rng.normal(0.62, 0.22), 0.12, 1.05)),
            K_SMC=float(np.clip(rng.lognormal(np.log(2.6e-4), 0.42),
                                5e-5, 9e-4)),
            HB=float(np.clip(rng.normal(13.8, 1.5), 7.0, 16.5)),
            MAP0=float(np.clip(rng.normal(92.0, 11.0), 68.0, 125.0)),
            PCA_INV=float(rng.choice([0.0, 1.0], p=[0.70, 0.30])),
            G_LEAK=float(np.clip(rng.lognormal(np.log(0.80), 0.55), 0.12, 3.0)),
        )
        s_["SEC"] = 1.0 if s_["RNF"] == 0.0 else 0.0
        oo = simulate(s_, 5 * 365, n=240)
        rows.append([s_["RNF"], s_["ANGIO"], s_["K_SMC"], s_["HB"],
                     s_["MAP0"], s_["PCA_INV"], s_["G_LEAK"],
                     oo["STEN"][-1], oo["gS"][-1],
                     oo["CBFA"][-1], oo["CBF_WS"][-1], oo["CVR_INTR"][-1],
                     oo["INFA"][-1] * 100, 1 - np.exp(-oo["ISCH"][-1]),
                     1 - np.exp(-oo["HEMH"][-1]), oo["COG"][-1]])
    A = np.array(rows)
    cols = ["RNF", "ANGIO", "K_SMC", "HB", "MAP0", "PCA_INV", "G_LEAK",
            "STEN", "gS", "CBF_A", "CBF_ws", "CVRi%", "INF%", "P_isch",
            "P_haem", "COGz"]
    print(f"  {'':>8} " + " ".join(f"{c:>8}" for c in cols[7:]))
    for q, lbl in ((10, "P10"), (50, "median"), (90, "P90")):
        v = np.percentile(A, q, axis=0)
        print(f"  {lbl:>8} " + " ".join(f"{v[j]:>8.3f}" for j in range(7, 16)))
    for tgt, tl in ((12, "5-y infarct fraction"), (13, "5-y ischaemic event"),
                    (14, "5-y haemorrhage")):
        print(f"\n  Correlates of {tl}:")
        for j, c in enumerate(cols):
            if j == tgt:
                continue
            r = np.corrcoef(A[:, j], A[:, tgt])[0, 1]
            if abs(r) > 0.12:
                print(f"    {c:>8}  r = {r:+.3f}")
    lo = A[:, 11] < 5.0
    print(f"\n  {lo.sum()/N*100:.1f}% of the population reaches 5 y with")
    print("  intrinsic reserve < 5% -- the pressure-passive subgroup.")
    print(f"  mean 5-y infarct   reserve<5%: {A[lo,12].mean():6.2f}%   "
          f"reserve>=5%: {A[~lo,12].mean():6.2f}%")
    print(f"  mean 5-y P(isch)   reserve<5%: {A[lo,13].mean():6.3f}    "
          f"reserve>=5%: {A[~lo,13].mean():6.3f}")
    pv = A[:, 5] > 0.5
    print(f"  mean 5-y P(haem)   PCA involved: {A[pv,14].mean():6.3f}   "
          f"PCA spared: {A[~pv,14].mean():6.3f}   "
          f"ratio {A[pv,14].mean()/max(A[~pv,14].mean(),1e-9):.1f}x")

    # ==================================================================
    hr("14.  SELECT BY PHYSIOLOGY, NOT BY ANGIOGRAPHY -- BUT WHICH?")
    q_res = np.percentile(A[:, 11], [25, 75])
    print("  FIRST, A NEGATIVE RESULT THAT CHANGES THE QUESTION.  Stratifying")
    print("  this population by INTRINSIC RESERVE at 5 y does not work:")
    print(f"    reserve quartile cut-points: {q_res[0]:.2f}% / {q_res[1]:.2f}%")
    print(f"    fraction of the cohort with reserve < 5%: "
          f"{(A[:,11] < 5).sum()/N*100:.1f}%")
    print("  Reserve is a SATURATING variable.  It is exquisitely informative")
    print("  while it is being spent (section 2: it hits zero at year 2 while")
    print("  flow is still normal) and it carries almost no information once")
    print("  it is gone, because by then everyone is at zero.  A cross-section")
    print("  of established disease therefore cannot be stratified by it -- and")
    print("  that is a statement about WHEN to measure, not about whether the")
    print("  measurement matters.")
    print("\n  WHAT STILL VARIES AT 5 YEARS IS THE MARGIN TO THRESHOLD.")
    cao2 = 1.34 * A[:, 3] * P["SAO2"] / 100.0
    crit = P["CMRO2_0"] / (cao2 * P["OEF_MAX"])
    margin = A[:, 10] - crit
    qm = np.percentile(margin, [25, 75])
    print(f"  margin = CBF_ws - CBF_crit:  P10 {np.percentile(margin,10):+.1f}"
          f"   median {np.percentile(margin,50):+.1f}"
          f"   P90 {np.percentile(margin,90):+.1f} mL/100g/min")
    print(f"  quartile cut-points: {qm[0]:+.2f} / {qm[1]:+.2f}")
    print("\n  Offer the SAME direct bypass at year 5 to the worst- and best-")
    print("  margin quartile, and read year 8.\n")
    print(f"  {'group':>22} {'n':>3} {'margin':>7} {'STEN':>6} {'gS':>5}"
          f" {'INF% none':>10} {'INF% byp':>9} {'gain':>6}"
          f" {'Pisch none':>11} {'Pisch byp':>10} {'ARR':>7}")
    for lbl, sel in (("worst margin (Q1)", margin <= qm[0]),
                     ("best margin (Q4)", margin >= qm[1])):
        idx = np.where(sel)[0][:36]
        g0, g1, p0, p1, st, gs = [], [], [], [], [], []
        for i in idx:
            base = dict(RNF=A[i, 0], ANGIO=A[i, 1], K_SMC=A[i, 2],
                        HB=A[i, 3], MAP0=A[i, 4], PCA_INV=A[i, 5],
                        G_LEAK=A[i, 6], SEC=1.0 if A[i, 0] == 0.0 else 0.0)
            r0 = simulate(base, 8 * 365, n=220)
            b2 = dict(base); b2.update(surg_kind="direct", surg_t=5 * 365.0)
            r1 = simulate(b2, 8 * 365, n=260)
            g0.append(r0["INFA"][-1] * 100); g1.append(r1["INFA"][-1] * 100)
            p0.append(1 - np.exp(-r0["ISCH"][-1]))
            p1.append(1 - np.exp(-r1["ISCH"][-1]))
            st.append(r0["STEN"][-1]); gs.append(r0["gS"][-1])
        g0, g1 = np.array(g0), np.array(g1)
        p0, p1 = np.array(p0), np.array(p1)
        print(f"  {lbl:>22} {len(idx):>3} {margin[sel][:36].mean():>+7.2f} "
              f"{np.mean(st):>6.3f} {np.mean(gs):>5.2f} {g0.mean():>10.2f} "
              f"{g1.mean():>9.2f} {g0.mean()-g1.mean():>6.2f} "
              f"{p0.mean():>11.3f} {p1.mean():>10.3f} "
              f"{p0.mean()-p1.mean():>7.3f}")
    print("\n  The two groups differ in what the tissue is actually getting,")
    print("  and they differ by ~7x in what the operation is worth.  But read")
    print("  the STEN column honestly: the angiogram is NOT blind here (mean")
    print("  stenosis 0.870 vs 0.543), so in this population the lesion and the")
    print("  physiology are correlated and a DSA does carry real prognostic")
    print("  information.  The defensible claim is the weaker one: stenosis is")
    print("  informative but NOT SUFFICIENT, because the same stenosis maps to")
    print("  different margins depending on collateral conductance and on")
    print("  haemoglobin (sections 2, 10 and 12).")
    print("\n  Taken with the negative result above, the model's selection")
    print("  claim is narrower and more useful than 'measure reserve':")
    print("    - BEFORE decompensation, reserve is the sensitive variable and")
    print("      flow is uninformative (section 2, year 1-2).")
    print("    - AFTER decompensation, reserve is saturated and the margin to")
    print("      CBF_crit is the variable that still separates patients -- and")
    print("      because CBF_crit depends on haemoglobin, the same CBF means")
    print("      different things in different patients (section 10).")

    hr("END OF REFERENCE OUTPUT")


if __name__ == "__main__":
    main()
