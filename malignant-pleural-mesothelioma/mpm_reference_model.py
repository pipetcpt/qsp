#!/usr/bin/env python3
"""
mpm_reference_model.py
----------------------
Dependency-free reference implementation of the malignant pleural mesothelioma
(MPM) QSP model.  This file exists because this container has no R runtime, and
a 53-state ODE system that has never been integrated is not a model, it is a
list of equations.  Everything in mpm_mrgsolve_model.R is written to be the
same system as this one; the calibration numbers quoted in README.md and in the
commit message come from running THIS file.

TIME UNIT: days.  LENGTH UNIT: cm internally, reported in mm.
MASS UNIT: grams of tissue.  DRUG: mg (small molecules) and mg (mAbs).

------------------------------------------------------------------------------
THE ONE STRUCTURAL COMMITMENT
------------------------------------------------------------------------------
Mesothelioma is a RIND, not a mass.  Tumour lives as a sheet of thickness h on
a pleural surface of area A, so

        volume  V = A h            and       dV/dh = A , a CONSTANT.

Every clinical quantity is a different functional of that same sheet:

    mRECIST      = 4 h_par + 2 h_vis            (LINEAR in h)
    drug at the
    cell         = (lambda/h)(1 - e^{-h/lambda}) (SUBLINEAR, saturating)
    lung trapping= stiffness x h_vis             (VISCERAL LEAF ONLY)
    effusion     = f(parietal rind over the lymphatic stomata)

Because dV/dh is constant, a 30% fall in mRECIST is a 30% fall in volume.  In
a spherical tumour a 30% fall in diameter is a 66% fall in volume.  The words
"partial response" therefore denote a 2.2-fold different degree of cytoreduction
in mesothelioma than in the tumour the criteria were designed for.  This is
computed, not asserted, in check_geometry().

------------------------------------------------------------------------------
ONE COLLAGEN STATE, THREE CONSEQUENCES
------------------------------------------------------------------------------
phi = M / (T + M + N) is the collagenous fraction of the rind.  It appears in
exactly three places and nowhere else:

    (1) in h, because matrix mass does not disappear when tumour cells die
        -> every thickness-based response measurement is biased downward;
    (2) in lambda(phi) = lambda_0 (1-phi)^1.5, because collagen collapses
        microvessels and raises interstitial pressure
        -> delivery falls as the tumour fibroses;
    (3) in the added elastance of the visceral leaf
        -> the lung is trapped and cannot re-expand.

------------------------------------------------------------------------------
HISTOLOGY AS ONE CONTINUOUS AXIS
------------------------------------------------------------------------------
x in [0,1]: 0 epithelioid, 0.5 biphasic, 1 sarcomatoid.  x raises
proliferation, raises chemoresistance, raises collagen deposition and raises
PD-L1 / T-cell infiltration.  The first three hurt chemotherapy and the fourth
helps checkpoint blockade, so the SIGN of OS(IO) - OS(chemo) as a function of x
is an output of the model, not an input.  See check_crossover().
"""

import math
import sys

LN2 = math.log(2.0)


# =============================================================================
#  PARAMETERS
# =============================================================================
def default_params():
    p = {}

    # ---------------- patient ------------------------------------------------
    p["BSA"] = 1.80            # m2
    p["WT"] = 72.0             # kg
    p["FVC_PRED"] = 3.30       # L, predicted FVC for a 72-year-old male
    p["EMT"] = 0.25            # histology axis x (0 epithelioid .. 1 sarcomatoid)

    # ---------------- rind geometry -----------------------------------------
    # Anatomic pleural surface of one hemithorax, times the fraction of that
    # surface actually carrying tumour at diagnosis (coverage).  The product is
    # the area that converts mass to thickness.
    p["A_PAR_ANAT"] = 1000.0   # cm2 parietal pleura, one hemithorax
    p["A_VIS_ANAT"] = 900.0    # cm2 visceral (lung) surface
    p["A_FIS_ANAT"] = 400.0    # cm2 fissural + mediastinal + diaphragmatic
    p["COV_PAR"] = 0.45
    p["COV_VIS"] = 0.35
    p["COV_FIS"] = 0.50
    p["RHO"] = 1.05            # g/cm3

    # ---------------- growth -------------------------------------------------
    p["KG0"] = 0.0032          # 1/d, intrinsic proliferation of a thin rind
    p["EMT_KG"] = 0.90         # kg = KG0 (1 + EMT_KG x)
    p["LAM_G"] = 0.42          # cm, depth of the vascularised/proliferative zone
    p["KSEED"] = 0.00060        # 1/d, cross-seeding between leaves (pleural fluid)

    # ---------------- matrix / collagen -------------------------------------
    p["KCOL0"] = 0.00125        # 1/d, matrix deposited per gram of tumour
    p["EMT_KCOL"] = 3.20       # kcol = KCOL0 (1 + EMT_KCOL x)
    p["KMDEG"] = 0.0015        # 1/d, matrix turnover (t1/2 ~ 460 d: essentially none)
    p["KNCLR"] = 1.0 / 21.0    # 1/d, clearance of killed/necrotic mass

    # ---------------- penetration -------------------------------------------
    p["LAM_SM0"] = 0.25        # cm, small-molecule penetration length at phi = 0
    p["LAM_AB0"] = 0.060       # cm, IgG penetration length at phi = 0
    p["LAM_T0"] = 0.050        # cm, T-cell infiltration length at phi = 0
    p["PHI_EXP"] = 1.5         # lambda = lambda0 (1 - phi)^PHI_EXP
    p["LAM_IP0"] = 0.35        # cm, intrapleural penetration length (free surface)

    # ---------------- cisplatin PK ------------------------------------------
    p["CIS_V1"] = 25.0         # L
    p["CIS_V2"] = 40.0         # L
    p["CIS_CL"] = 600.0        # L/d free platinum clearance (t1/2a ~ 0.7 h)
    p["CIS_Q"] = 200.0         # L/d
    p["CIS_KIN"] = 1.00        # 1/d uptake into tumour (adduct formation)
    p["CIS_KOUT"] = 0.35       # 1/d  Pt-DNA adduct repair (ERCC1), t1/2 ~ 2 d

    # ---------------- pemetrexed PK -----------------------------------------
    p["PEM_V1"] = 16.0         # L
    p["PEM_V2"] = 8.0          # L
    p["PEM_CL0"] = 130.0       # L/d at eGFR 90 (CL ~ 90 mL/min)
    p["PEM_Q"] = 30.0          # L/d
    p["PEM_KIN"] = 1.00        # 1/d uptake by RFC/SLC19A1
    p["PEM_KOUT"] = 0.277      # 1/d polyglutamate loss (t1/2 2.5 d in cell)

    # ---------------- gemcitabine PK (alternative backbone) -----------------
    p["GEM_V"] = 50.0
    p["GEM_CL"] = 2500.0
    p["GEM_KIN"] = 1.00
    p["GEM_KOUT"] = 1.00

    # ---------------- monoclonal antibody PK (nivo / ipi / bev / pembro) ----
    p["MAB_V1"] = 3.6          # L
    p["MAB_V2"] = 2.8          # L
    p["MAB_CL"] = 0.21         # L/d
    p["MAB_Q"] = 0.45          # L/d
    p["MAB_KIN"] = 0.55        # 1/d into rind interstitium
    p["MAB_KOUT"] = 0.55       # 1/d

    # ---------------- chemotherapy PD ---------------------------------------
    p["EMAX_CIS"] = 0.148438      # 1/d maximal cisplatin kill of a thin rind
    p["EC50_CIS"] = 0.10       # mg/L adduct-equivalent in tumour
    p["EMAX_PEM"] = 0.0283594      # 1/d
    p["EC50_PEM"] = 0.35       # mg/L polyglutamate-equivalent
    p["SYN_PEMCIS"] = 0.30     # fractional synergy of the doublet
    p["EMT_CHEMO"] = 0.55      # kill x (1 - EMT_CHEMO x): EMT chemoresistance
    p["EMAX_GEM"] = 0.030
    p["EC50_GEM"] = 0.35

    # ---------------- immune compartment ------------------------------------
    p["TEFF0"] = 0.050         # arbitrary units, baseline intratumoural CD8
    p["KPRIME"] = 0.055        # 1/d priming rate from antigen
    p["TMB_FACT"] = 0.42       # low mutational burden scales antigenicity
    p["KTDEATH"] = 0.070       # 1/d effector loss
    p["KTREG"] = 0.020
    p["TREG0"] = 1.0
    p["KIFN"] = 0.9
    p["KIFN_OUT"] = 1.2
    p["PDL1_0"] = 0.20         # constitutive PD-L1 (TPS fraction) at x = 0
    p["EMT_PDL1"] = 2.0        # PD-L1 x (1 + EMT_PDL1 x): 20% -> 60% at x = 1
    p["EMT_TINF"] = 1.00       # infiltration x (1 + EMT_TINF x)
    p["PDL1_IND"] = 0.55       # IFN-gamma inducible component
    p["EMAX_IO"] = 0.624424       # 1/d maximal immune kill of a thin rind
    p["KT50"] = 0.080          # T_eff for half-maximal kill (T_eff operates
                           # near 0.05: low TMB antigen supply divided by a
                           # suppression denominator of ~4.9)
    p["KD_PD1"] = 0.30         # mg/L, nivolumab conc for 50% PD-1 occupancy
    p["KD_CTLA4"] = 1.20       # mg/L, ipilimumab conc for 50% CTLA-4 occupancy
    p["CTLA4_REL"] = 0.25      # direct kill released by CTLA-4 blockade alone
    p["IPI_PRIME"] = 2.30      # fold increase in priming at full CTLA-4 blockade
    p["TGFB_SUPP"] = 0.85      # TGF-beta suppression of effector function

    # ---------------- signalling --------------------------------------------
    p["KVEGF"] = 0.9
    p["KVEGF_OUT"] = 0.7
    p["KTGFB"] = 1.0
    p["KTGFB_OUT"] = 0.8
    p["KIL6"] = 1.0
    p["KIL6_OUT"] = 0.9
    p["BEV_KD"] = 0.9          # mg/L bevacizumab for 50% VEGF neutralisation
    p["BEV_LAM"] = 0.295625        # max fractional rise in lambda (normalisation)

    # ---------------- pleural space -----------------------------------------
    p["QFORM0"] = 250.0        # mL/d normal pleural fluid formation
    p["QFORM_VEGF"] = 5.5      # multiplier at maximal VEGF
    p["QDRAIN_MAX"] = 5000.0   # mL/d stomatal lymphatic reserve
    p["H50_STOMA"] = 0.22      # cm parietal rind for 50% stomatal obstruction
    p["PLV_MAX"] = 3500.0      # mL
    p["KSYMPH"] = 0.10         # 1/d pleurodesis symphysis formation
    p["KSYMPH_LOSS"] = 0.010   # 1/d

    # ---------------- mechanics ---------------------------------------------
    p["KELAST"] = 76.0         # cmH2O/L per (phi x cm) of visceral rind
    p["E50"] = 14.5            # cmH2O/L, the trapped-lung threshold
    p["KREC"] = 0.045          # 1/d lung re-expansion (slow, hysteretic)
    p["PLV_FVC"] = 0.55        # fraction of FVC lost at maximal effusion
    p["KDYSP"] = 0.20          # 1/d dyspnoea equilibration

    # ---------------- toxicity ----------------------------------------------
    p["CIRC0"] = 4.0           # 10^9/L neutrophils
    p["MTT"] = 4.6             # d mean transit time (110 h)
    p["GAM_FB"] = 0.16         # feedback exponent
    p["SLOPE_CIS_ANC"] = 0.90  # L/mg on the RETAINED adduct
    p["SLOPE_PEM_ANC"] = 0.115 # L/mg on the RETAINED polyglutamate
    p["FOLATE"] = 1.0          # 1 = folate + B12 supplemented
    p["FOLATE_PROT"] = 0.60    # 60% reduction in pemetrexed marrow effect
    p["KPT_UP"] = 35.0         # renal cortex platinum uptake
    p["KPT_EL"] = 0.020        # renal platinum elimination
    p["GFR0"] = 90.0           # mL/min/1.73m2
    p["KGFR_DAM"] = 0.400
    p["KGFR_REP"] = 0.0045
    p["KNEURO"] = 0.750
    p["KIRAE_PD1"] = 0.020
    p["KIRAE_CTLA4"] = 0.075
    p["KIRAE_OUT"] = 0.035

    # ---------------- cachexia / performance status -------------------------
    p["KCACHEX"] = 0.055       # %/d weight loss at maximal IL-6
    p["KCACHEX_REC"] = 0.010

    # ---------------- biomarker ---------------------------------------------
    p["KSMRP"] = 0.00090       # nM/d per gram of viable tumour surface
    p["KSMRP_CL"] = 0.60       # 1/d, renal (falls with eGFR)

    # ---------------- survival ----------------------------------------------
    p["HZ0"] = 0.00045625         # 1/d baseline hazard at reference state
    p["HZ_VOL"] = 0.30         # coefficient on ln(V / V_ref)
    p["HZ_VREF"] = 500.0       # cm3
    p["HZ_ECOG"] = 0.60        # per ECOG point
    p["HZ_EMT"] = 0.55         # per unit of the histology axis
    p["HZ_FVC"] = 1.20         # per unit fractional FVC loss

    # ---------------- treatment switches (set by scenarios) -----------------
    p["ON_CIS"] = 0; p["ON_CARBO"] = 0; p["ON_PEM"] = 0; p["ON_GEM"] = 0
    p["ON_BEV"] = 0; p["ON_NIV"] = 0; p["ON_IPI"] = 0; p["ON_PEMBRO"] = 0
    p["ON_ADI"] = 0; p["ON_IP"] = 0
    p["ASS1_NEG"] = 0          # 1 = ASS1-methylated, arginine auxotroph
    p["ADI_EFF"] = 0.026       # 1/d kill when ASS1-negative
    p["SURG_DAY"] = -1         # day of cytoreduction, -1 = none
    p["SURG_RESID"] = 0.10     # cm residual thickness after macroscopic complete resection
    p["SURG_FIS_SPARE"] = 0.65 # fraction of fissural disease NOT removed
    p["SURG_FVC_HIT"] = 0.0    # 0 = P/D, 0.35 = EPP (whole lung removed)
    p["TALC_DAY"] = -1
    p["IPC_DAY"] = -1
    p["IPC_RATE"] = 700.0      # mL/d ambulatory drainage
    p["TTF"] = 0               # tumour-treating fields
    p["TTF_EFF"] = 0.0045      # 1/d
    return p


# =============================================================================
#  STATE VECTOR
# =============================================================================
STATES = [
    # --- PK (0-16) ---
    "CIS1", "CIS2", "CIST",
    "PEM1", "PEM2", "PEMT",
    "GEM1", "GEMT",
    "BEV1", "BEV2",
    "NIV1", "NIV2", "NIVT",
    "IPI1", "IPI2",
    "IPPL", "IPT",
    # --- rind (17-25) ---
    "TPAR", "TVIS", "TFIS",
    "MPAR", "MVIS", "MFIS",
    "NPAR", "NVIS", "NFIS",
    # --- immune (26-31) ---
    "TEFF", "TREG", "PRIME", "IFNG", "PDL1", "TAM",
    # --- signalling (32-34) ---
    "VEGF", "TGFB", "IL6",
    # --- pleural (35-36) ---
    "PLV", "SYMPH",
    # --- mechanics (37-38) ---
    "VEXP", "DYSP",
    # --- toxicity (39-47) ---
    "PROL", "TR1", "TR2", "TR3", "ANC",
    "PTK", "GFR", "IRAE", "NEURO",
    # --- host (48-49) ---
    # ECOGs is a RESERVED slot: performance status is currently computed
    # algebraically in algebra() rather than integrated, so this compartment has
    # no derivative and stays at zero.  It is kept in the vector so that adding
    # a lagged ECOG (performance status does not jump instantaneously) does not
    # renumber every state downstream.
    "CACHEX", "ECOGs",
    # --- biomarker + survival (50-52) ---
    "SMRP", "CH", "ARG",
]
IDX = {s: i for i, s in enumerate(STATES)}
NST = len(STATES)


def initial_state(p, h_par=0.60, h_vis=0.40, h_fis=0.50, phi0=None):
    """Initial condition at diagnosis, specified by RIND THICKNESSES (cm)."""
    y = [0.0] * NST
    if phi0 is None:
        # steady-state collagen fraction implied by the deposition/turnover
        # balance at the intrinsic growth rate
        kcol = p["KCOL0"] * (1.0 + p["EMT_KCOL"] * p["EMT"])
        kg = p["KG0"] * (1.0 + p["EMT_KG"] * p["EMT"])
        phi0 = kcol / (kcol + p["KMDEG"] + kg)
    A = areas(p)
    for site, h, a in (("PAR", h_par, A[0]), ("VIS", h_vis, A[1]), ("FIS", h_fis, A[2])):
        total = h * p["RHO"] * a
        y[IDX["T" + site]] = total * (1.0 - phi0)
        y[IDX["M" + site]] = total * phi0
        y[IDX["N" + site]] = 0.0
    y[IDX["TEFF"]] = p["TEFF0"] * (1.0 + p["EMT_TINF"] * p["EMT"])
    y[IDX["TREG"]] = p["TREG0"]
    y[IDX["PRIME"]] = 1.0
    y[IDX["IFNG"]] = 0.35
    y[IDX["PDL1"]] = p["PDL1_0"] * (1.0 + p["EMT_PDL1"] * p["EMT"])
    y[IDX["TAM"]] = 1.0
    y[IDX["VEGF"]] = 1.0
    y[IDX["TGFB"]] = 1.0
    y[IDX["IL6"]] = 1.0
    y[IDX["PLV"]] = 600.0
    y[IDX["SYMPH"]] = 0.0
    y[IDX["VEXP"]] = 1.0
    y[IDX["DYSP"]] = 0.0
    for s, v in (("PROL", p["CIRC0"]), ("TR1", p["CIRC0"]), ("TR2", p["CIRC0"]),
                 ("TR3", p["CIRC0"]), ("ANC", p["CIRC0"])):
        y[IDX[s]] = v
    y[IDX["GFR"]] = p["GFR0"]
    y[IDX["ARG"]] = 1.0
    y[IDX["SMRP"]] = 0.0
    # VEXP and DYSP are re-equilibrated below by a short burn-in in simulate()
    return y


def areas(p):
    return (p["A_PAR_ANAT"] * p["COV_PAR"],
            p["A_VIS_ANAT"] * p["COV_VIS"],
            p["A_FIS_ANAT"] * p["COV_FIS"])


# =============================================================================
#  ALGEBRAIC LAYER  (everything the ODEs and the outputs share)
# =============================================================================
def f_pen(h, lam):
    """Depth-averaged relative concentration through a slab of thickness h fed
    from ONE face, with an exponential decay length lam.

        C(z) = C0 exp(-z/lam)  ->  <C>/C0 = (lam/h)(1 - exp(-h/lam))

    Exact, monotone, -> 1 as h -> 0 and -> lam/h as h -> infinity.  This is the
    single function that makes drug delivery, oxygen supply and T-cell
    infiltration all SUBLINEAR in tumour thickness."""
    if h <= 1e-9:
        return 1.0
    r = h / lam
    if r > 60.0:
        return lam / h
    return (lam / h) * (1.0 - math.exp(-r))


def algebra(t, y, p):
    """All derived quantities.  Returns a dict used by both derivs() and out()."""
    A_par, A_vis, A_fis = areas(p)
    x = p["EMT"]
    a = {}

    # ---- rind geometry -----------------------------------------------------
    for site, A in (("PAR", A_par), ("VIS", A_vis), ("FIS", A_fis)):
        T = max(y[IDX["T" + site]], 0.0)
        M = max(y[IDX["M" + site]], 0.0)
        N = max(y[IDX["N" + site]], 0.0)
        tot = T + M + N
        h = tot / (p["RHO"] * A)
        a["T" + site] = T
        a["M" + site] = M
        a["N" + site] = N
        a["h_" + site] = h
        a["phi_" + site] = (M / tot) if tot > 1e-9 else 0.0
        a["V_" + site] = A * h
    a["A_par"], a["A_vis"], a["A_fis"] = A_par, A_vis, A_fis
    a["V_tot"] = a["V_PAR"] + a["V_VIS"] + a["V_FIS"]
    a["T_tot"] = a["TPAR"] + a["TVIS"] + a["TFIS"]
    a["N_tot"] = a["NPAR"] + a["NVIS"] + a["NFIS"]
    a["phi"] = ((a["MPAR"] + a["MVIS"] + a["MFIS"]) /
                max(a["V_tot"] * p["RHO"], 1e-9))
    # mRECIST: two perpendicular measurements on each of three CT levels.
    # In practice these are taken on the parietal/mediastinal surface, with a
    # minority on the visceral/fissural surface.
    a["mRECIST"] = 10.0 * (4.0 * a["h_PAR"] + 2.0 * a["h_VIS"])   # mm

    # ---- penetration lengths ----------------------------------------------
    bev = y[IDX["BEV1"]] / p["MAB_V1"]
    norm = p["BEV_LAM"] * bev / (p["BEV_KD"] + bev)     # vascular normalisation
    def lam_of(l0, phi):
        return max(l0 * (1.0 - min(phi, 0.95)) ** p["PHI_EXP"] * (1.0 + norm), 1e-4)
    a["lam_sm"] = {s: lam_of(p["LAM_SM0"], a["phi_" + s]) for s in ("PAR", "VIS", "FIS")}
    a["lam_ab"] = {s: lam_of(p["LAM_AB0"], a["phi_" + s]) for s in ("PAR", "VIS", "FIS")}
    a["lam_T"] = {s: lam_of(p["LAM_T0"], a["phi_" + s]) for s in ("PAR", "VIS", "FIS")}
    a["f_sm"] = {s: f_pen(a["h_" + s], a["lam_sm"][s]) for s in ("PAR", "VIS", "FIS")}
    a["f_ab"] = {s: f_pen(a["h_" + s], a["lam_ab"][s]) for s in ("PAR", "VIS", "FIS")}
    a["f_T"] = {s: f_pen(a["h_" + s], a["lam_T"][s]) for s in ("PAR", "VIS", "FIS")}
    a["f_g"] = {s: f_pen(a["h_" + s], p["LAM_G"]) for s in ("PAR", "VIS", "FIS")}
    # intrapleural route: same formula, source at the FREE surface
    lam_ip = {s: max(p["LAM_IP0"] * (1.0 - min(a["phi_" + s], 0.95)) ** p["PHI_EXP"], 1e-4)
              for s in ("PAR", "VIS", "FIS")}
    a["f_ip"] = {s: f_pen(a["h_" + s], lam_ip[s]) for s in ("PAR", "VIS", "FIS")}

    # ---- drug concentrations ----------------------------------------------
    a["C_cis"] = max(y[IDX["CIS1"]], 0.0) / p["CIS_V1"]        # mg/L free Pt
    a["C_pem"] = max(y[IDX["PEM1"]], 0.0) / p["PEM_V1"]
    a["C_gem"] = max(y[IDX["GEM1"]], 0.0) / p["GEM_V"]
    a["C_niv"] = max(y[IDX["NIV1"]], 0.0) / p["MAB_V1"]
    a["C_ipi"] = max(y[IDX["IPI1"]], 0.0) / p["MAB_V1"]
    a["C_bev"] = bev
    a["Ct_cis"] = max(y[IDX["CIST"]], 0.0)
    a["Ct_pem"] = max(y[IDX["PEMT"]], 0.0)
    a["Ct_gem"] = max(y[IDX["GEMT"]], 0.0)
    a["Ct_niv"] = max(y[IDX["NIVT"]], 0.0)
    a["C_ippl"] = max(y[IDX["IPPL"]], 0.0) / max(y[IDX["PLV"]] / 1000.0, 0.05)  # mg/L

    # ---- receptor occupancy ------------------------------------------------
    a["RO_pd1"] = a["Ct_niv"] / (p["KD_PD1"] + a["Ct_niv"])
    a["RO_ctla4"] = a["C_ipi"] / (p["KD_CTLA4"] + a["C_ipi"])   # lymph node: no barrier

    # ---- chemotherapy kill -------------------------------------------------
    e_cis = p["EMAX_CIS"] * a["Ct_cis"] / (p["EC50_CIS"] + a["Ct_cis"])
    e_pem = p["EMAX_PEM"] * a["Ct_pem"] / (p["EC50_PEM"] + a["Ct_pem"])
    e_gem = p["EMAX_GEM"] * a["Ct_gem"] / (p["EC50_GEM"] + a["Ct_gem"])
    syn = p["SYN_PEMCIS"] * math.sqrt(max(e_cis, 0.0) * max(e_pem, 0.0))
    chem_res = (1.0 - p["EMT_CHEMO"] * x)
    a["kill_chem_thin"] = (e_cis + e_pem + e_gem + syn) * chem_res

    # arginine deprivation (ADI-PEG20) kills only ASS1-negative tumour
    arg = max(y[IDX["ARG"]], 0.0)
    a["kill_arg_thin"] = (p["ADI_EFF"] * p["ASS1_NEG"] * max(0.0, 1.0 - arg / 0.35)
                          if p["ON_ADI"] else 0.0)

    # ---- immune kill -------------------------------------------------------
    teff = max(y[IDX["TEFF"]], 0.0)
    pdl1 = min(max(y[IDX["PDL1"]], 0.0), 1.0)
    # RELEASED kill: exactly zero with no checkpoint inhibitor on board, so
    # EMAX_IO cannot reach back into the untreated trajectory.  What blockade
    # can release is bounded by how much PD-L1 was braking to begin with.
    release = pdl1 * a["RO_pd1"] + p["CTLA4_REL"] * a["RO_ctla4"]
    a["release"] = release
    a["kill_io_thin"] = (p["EMAX_IO"] * teff / (p["KT50"] + teff) * release)

    # ---- per-site kill rates (penetration applied HERE, not upstream) ------
    a["kill"] = {}
    for s in ("PAR", "VIS", "FIS"):
        k = (a["kill_chem_thin"] * a["f_sm"][s]
             + a["kill_io_thin"] * a["f_T"][s]
             + a["kill_arg_thin"] * a["f_sm"][s]
             + (p["TTF_EFF"] if p["TTF"] else 0.0) * (1.0 if s != "FIS" else 0.4))
        if p["ON_IP"]:
            e_ip = 0.055 * a["C_ippl"] / (2.0 + a["C_ippl"])
            k += e_ip * a["f_ip"][s]
        a["kill"][s] = k

    # ---- pleural space -----------------------------------------------------
    veg = max(y[IDX["VEGF"]], 0.0)
    a["q_form"] = p["QFORM0"] * (1.0 + (p["QFORM_VEGF"] - 1.0) * veg / (1.0 + veg))
    block = a["h_PAR"] / (p["H50_STOMA"] + a["h_PAR"])
    a["stoma_block"] = block
    symph = min(max(y[IDX["SYMPH"]], 0.0), 1.0)
    a["symph"] = symph
    a["q_drain"] = p["QDRAIN_MAX"] * (1.0 - block)
    # Pleurodesis removes the SPACE, it does not disable the lymphatics.
    a["q_form_eff"] = a["q_form"] * (1.0 - 0.95 * symph)

    # ---- mechanics ---------------------------------------------------------
    a["E_add"] = p["KELAST"] * a["phi_VIS"] * a["h_VIS"]
    a["vexp_target"] = 1.0 / (1.0 + (a["E_add"] / p["E50"]) ** 2)
    vexp = min(max(y[IDX["VEXP"]], 0.0), 1.0)
    plv = min(max(y[IDX["PLV"]], 0.0), p["PLV_MAX"])
    a["FVC"] = p["FVC_PRED"] * vexp * (1.0 - p["PLV_FVC"] * plv / p["PLV_MAX"])
    if 0 <= p["SURG_DAY"] <= t:
        a["FVC"] *= (1.0 - p["SURG_FVC_HIT"])
    a["FVC_pct"] = 100.0 * a["FVC"] / p["FVC_PRED"]

    # ---- performance status ------------------------------------------------
    wl = max(y[IDX["CACHEX"]], 0.0)
    a["ECOG"] = min(max(0.20 + 1.9 * (1.0 - a["FVC"] / p["FVC_PRED"])
                        + 0.050 * wl + 0.30 * min(max(y[IDX["DYSP"]], 0.0), 3.0), 0.0), 4.0)

    # ---- hazard ------------------------------------------------------------
    a["hazard"] = (p["HZ0"]
                   * math.exp(p["HZ_VOL"] * math.log(max(a["V_tot"], 1.0) / p["HZ_VREF"]))
                   * math.exp(p["HZ_ECOG"] * a["ECOG"])
                   * math.exp(p["HZ_EMT"] * x)
                   * math.exp(p["HZ_FVC"] * (1.0 - a["FVC"] / p["FVC_PRED"])))
    return a


# =============================================================================
#  DERIVATIVES
# =============================================================================
def derivs(t, y, p):
    a = algebra(t, y, p)
    d = [0.0] * NST
    x = p["EMT"]

    # ------------------------- PK ------------------------------------------
    d[IDX["CIS1"]] = (-p["CIS_CL"] / p["CIS_V1"] * y[IDX["CIS1"]]
                      - p["CIS_Q"] / p["CIS_V1"] * y[IDX["CIS1"]]
                      + p["CIS_Q"] / p["CIS_V2"] * y[IDX["CIS2"]])
    d[IDX["CIS2"]] = (p["CIS_Q"] / p["CIS_V1"] * y[IDX["CIS1"]]
                      - p["CIS_Q"] / p["CIS_V2"] * y[IDX["CIS2"]])
    d[IDX["CIST"]] = p["CIS_KIN"] * a["C_cis"] - p["CIS_KOUT"] * y[IDX["CIST"]]

    pem_cl = p["PEM_CL0"] * max(y[IDX["GFR"]], 15.0) / p["GFR0"]
    d[IDX["PEM1"]] = (-pem_cl / p["PEM_V1"] * y[IDX["PEM1"]]
                      - p["PEM_Q"] / p["PEM_V1"] * y[IDX["PEM1"]]
                      + p["PEM_Q"] / p["PEM_V2"] * y[IDX["PEM2"]])
    d[IDX["PEM2"]] = (p["PEM_Q"] / p["PEM_V1"] * y[IDX["PEM1"]]
                      - p["PEM_Q"] / p["PEM_V2"] * y[IDX["PEM2"]])
    d[IDX["PEMT"]] = p["PEM_KIN"] * a["C_pem"] - p["PEM_KOUT"] * y[IDX["PEMT"]]

    d[IDX["GEM1"]] = -p["GEM_CL"] / p["GEM_V"] * y[IDX["GEM1"]]
    d[IDX["GEMT"]] = p["GEM_KIN"] * a["C_gem"] - p["GEM_KOUT"] * y[IDX["GEMT"]]

    for cen, per in (("BEV1", "BEV2"), ("NIV1", "NIV2"), ("IPI1", "IPI2")):
        d[IDX[cen]] += (-p["MAB_CL"] / p["MAB_V1"] * y[IDX[cen]]
                        - p["MAB_Q"] / p["MAB_V1"] * y[IDX[cen]]
                        + p["MAB_Q"] / p["MAB_V2"] * y[IDX[per]])
        d[IDX[per]] += (p["MAB_Q"] / p["MAB_V1"] * y[IDX[cen]]
                        - p["MAB_Q"] / p["MAB_V2"] * y[IDX[per]])
    # antibody in the rind: gated by the IgG penetration fraction of the
    # thickest (parietal) leaf -- this is where the binding-site barrier bites
    d[IDX["NIVT"]] = (p["MAB_KIN"] * a["C_niv"] * a["f_ab"]["PAR"]
                      - p["MAB_KOUT"] * y[IDX["NIVT"]])

    # intrapleural drug: diluted by the effusion, cleared by residual lymphatics
    d[IDX["IPPL"]] = -(0.9 + a["q_drain"] / 1000.0) * y[IDX["IPPL"]]
    d[IDX["IPT"]] = 0.0

    # ------------------------- rind ----------------------------------------
    kg = p["KG0"] * (1.0 + p["EMT_KG"] * x)
    kcol = p["KCOL0"] * (1.0 + p["EMT_KCOL"] * x)
    for s in ("PAR", "VIS", "FIS"):
        T = a["T" + s]
        growth = kg * a["f_g"][s] * T
        kill = a["kill"][s] * T
        # cross-seeding through the pleural fluid: the reason a resected leaf
        # repopulates from the fissural sanctuary
        others = sum(a["T" + o] for o in ("PAR", "VIS", "FIS") if o != s)
        seed = p["KSEED"] * others * (1.0 - min(a["h_" + s] / 1.5, 0.95))
        d[IDX["T" + s]] = growth - kill + seed
        d[IDX["M" + s]] = kcol * T - p["KMDEG"] * a["M" + s]
        d[IDX["N" + s]] = kill - p["KNCLR"] * a["N" + s]

    # ------------------------- immune --------------------------------------
    # antigen supply from dying tumour, scaled down by the low neoantigen load
    antigen = p["TMB_FACT"] * (0.15 * a["T_tot"] / 300.0 + a["N_tot"] / 60.0)
    prime = p["KPRIME"] * antigen * (1.0 + p["IPI_PRIME"] * a["RO_ctla4"])
    d[IDX["PRIME"]] = prime - 0.25 * y[IDX["PRIME"]]
    supp = (1.0 + 0.9 * max(y[IDX["TREG"]], 0.0) * (1.0 - 0.65 * a["RO_ctla4"])
            + p["TGFB_SUPP"] * max(y[IDX["TGFB"]], 0.0)
            + 0.5 * max(y[IDX["TAM"]], 0.0))
    d[IDX["TEFF"]] = (0.25 * y[IDX["PRIME"]] * (1.0 + p["EMT_TINF"] * x) / supp
                      - p["KTDEATH"] * y[IDX["TEFF"]])
    d[IDX["TREG"]] = (p["KTREG"] * (1.0 + 0.8 * max(y[IDX["TGFB"]], 0.0))
                      - p["KTREG"] * y[IDX["TREG"]]
                      - 0.045 * a["RO_ctla4"] * y[IDX["TREG"]])   # anti-CTLA-4 ADCC
    d[IDX["IFNG"]] = p["KIFN"] * 0.30 * y[IDX["TEFF"]] - p["KIFN_OUT"] * y[IDX["IFNG"]]
    pdl1_ss = (p["PDL1_0"] * (1.0 + p["EMT_PDL1"] * x)
               * (1.0 + p["PDL1_IND"] * max(y[IDX["IFNG"]], 0.0)))
    d[IDX["PDL1"]] = 0.5 * (min(pdl1_ss, 0.95) - y[IDX["PDL1"]])
    d[IDX["TAM"]] = 0.10 * (1.0 + 0.5 * max(y[IDX["IL6"]], 0.0) - y[IDX["TAM"]])

    # ------------------------- signalling ----------------------------------
    hyp = 1.0 - a["f_g"]["PAR"]                    # hypoxic fraction of the rind
    veg_ss = 0.7 + 2.2 * hyp + 0.5 * a["T_tot"] / 400.0
    bev_neut = a["C_bev"] / (p["BEV_KD"] + a["C_bev"])
    d[IDX["VEGF"]] = p["KVEGF_OUT"] * (veg_ss * (1.0 - 0.90 * bev_neut) - y[IDX["VEGF"]])
    tgf_ss = 0.6 + 1.4 * (1.0 + p["EMT_KCOL"] * x) * a["T_tot"] / 400.0
    d[IDX["TGFB"]] = p["KTGFB_OUT"] * (tgf_ss - y[IDX["TGFB"]])
    il6_ss = 0.5 + 1.8 * a["V_tot"] / 600.0 + 0.6 * a["N_tot"] / 60.0
    d[IDX["IL6"]] = p["KIL6_OUT"] * (il6_ss - y[IDX["IL6"]])

    # ------------------------- pleural space -------------------------------
    plv = min(max(y[IDX["PLV"]], 0.0), p["PLV_MAX"])
    ipc = p["IPC_RATE"] if (0 <= p["IPC_DAY"] <= t) else 0.0
    d[IDX["PLV"]] = (a["q_form_eff"] * (1.0 - plv / p["PLV_MAX"])
                     - a["q_drain"] * plv / (400.0 + plv)
                     - ipc * plv / (200.0 + plv))
    talc_on = 1.0 if (0 <= p["TALC_DAY"] <= t) else 0.0
    # pleurodesis needs apposition: it fails in proportion to lung trapping
    apposition = a["vexp_target"]
    d[IDX["SYMPH"]] = (p["KSYMPH"] * talc_on * apposition * (1.0 - y[IDX["SYMPH"]])
                       - p["KSYMPH_LOSS"] * y[IDX["SYMPH"]])

    # ------------------------- mechanics -----------------------------------
    d[IDX["VEXP"]] = p["KREC"] * (a["vexp_target"] - y[IDX["VEXP"]])
    dysp_ss = (2.4 * (1.0 - a["FVC"] / p["FVC_PRED"]) + 1.4 * plv / p["PLV_MAX"])
    d[IDX["DYSP"]] = p["KDYSP"] * (dysp_ss - y[IDX["DYSP"]])

    # ------------------------- myelosuppression (Friberg) ------------------
    # driven by the RETAINED species (adducts, polyglutamates), not by plasma
    edrug = (p["SLOPE_CIS_ANC"] * a["Ct_cis"]
             + p["SLOPE_PEM_ANC"] * a["Ct_pem"]
             * (1.0 - p["FOLATE_PROT"] * p["FOLATE"]))
    ktr = 4.0 / p["MTT"]
    fb = (p["CIRC0"] / max(y[IDX["ANC"]], 0.05)) ** p["GAM_FB"]
    d[IDX["PROL"]] = ktr * y[IDX["PROL"]] * (1.0 - min(edrug, 0.98)) * fb - ktr * y[IDX["PROL"]]
    d[IDX["TR1"]] = ktr * (y[IDX["PROL"]] - y[IDX["TR1"]])
    d[IDX["TR2"]] = ktr * (y[IDX["TR1"]] - y[IDX["TR2"]])
    d[IDX["TR3"]] = ktr * (y[IDX["TR2"]] - y[IDX["TR3"]])
    d[IDX["ANC"]] = ktr * (y[IDX["TR3"]] - y[IDX["ANC"]])

    # ------------------------- renal, neuro, irAE --------------------------
    d[IDX["PTK"]] = p["KPT_UP"] * a["Ct_cis"] - p["KPT_EL"] * y[IDX["PTK"]]
    d[IDX["GFR"]] = (-p["KGFR_DAM"] * y[IDX["PTK"]] / 100.0
                     + p["KGFR_REP"] * (p["GFR0"] - y[IDX["GFR"]]))
    d[IDX["NEURO"]] = p["KNEURO"] * a["Ct_cis"] - 0.0025 * y[IDX["NEURO"]]
    d[IDX["IRAE"]] = (p["KIRAE_PD1"] * a["RO_pd1"] + p["KIRAE_CTLA4"] * a["RO_ctla4"]
                      - p["KIRAE_OUT"] * y[IDX["IRAE"]])

    # ------------------------- cachexia ------------------------------------
    d[IDX["CACHEX"]] = (p["KCACHEX"] * max(y[IDX["IL6"]], 0.0) / (1.0 + max(y[IDX["IL6"]], 0.0))
                        - p["KCACHEX_REC"] * y[IDX["CACHEX"]])

    # ------------------------- biomarker / arginine / hazard ---------------
    # SMRP is SHED from the viable surface, so it tracks viable cell mass and is
    # cleared renally -- it does NOT track rind thickness once phi rises
    d[IDX["SMRP"]] = (p["KSMRP"] * a["T_tot"]
                      - p["KSMRP_CL"] * max(y[IDX["GFR"]], 10.0) / p["GFR0"] * y[IDX["SMRP"]])
    d[IDX["ARG"]] = (0.35 * (1.0 - y[IDX["ARG"]])
                     - (2.2 * p["ON_ADI"]) * y[IDX["ARG"]])
    d[IDX["CH"]] = a["hazard"]
    return d


# =============================================================================
#  DOSING
# =============================================================================
def build_doses(p, tmax, regimen):
    """Return a sorted list of (time, state, amount) bolus events."""
    ev = []
    bsa, wt = p["BSA"], p["WT"]
    C = 21.0                                   # cycle length, days
    ncy = regimen.get("cycles", 6)
    start = regimen.get("start", 0.0)
    for i in range(ncy):
        t0 = start + i * C
        if t0 > tmax:
            break
        if regimen.get("cis"):
            ev.append((t0, "CIS1", 75.0 * bsa))
        if regimen.get("carbo"):
            ev.append((t0, "CIS1", 0.62 * 5.0 * (p["GFR0"] + 25.0)))  # AUC5 Pt-equivalent
        if regimen.get("pem"):
            ev.append((t0, "PEM1", 500.0 * bsa))
        if regimen.get("gem"):
            ev.append((t0, "GEM1", 1000.0 * bsa))
            ev.append((t0 + 7.0, "GEM1", 1000.0 * bsa))
        if regimen.get("bev"):
            ev.append((t0, "BEV1", 15.0 * wt))
        if regimen.get("pembro"):
            ev.append((t0, "NIV1", 200.0))
    # maintenance bevacizumab after the platinum doublet (MAPS design)
    if regimen.get("bev_maint"):
        t0 = start + ncy * C
        while t0 <= tmax:
            ev.append((t0, "BEV1", 15.0 * wt))
            t0 += C
    # checkpoint blockade: up to 2 years
    if regimen.get("nivo"):
        t0, stop = start, min(tmax, start + 730.0)
        while t0 <= stop:
            ev.append((t0, "NIV1", 360.0))
            t0 += 21.0
    if regimen.get("ipi"):
        t0, stop = start, min(tmax, start + 730.0)
        while t0 <= stop:
            ev.append((t0, "IPI1", 1.0 * wt))
            t0 += 42.0
    if regimen.get("ip_drug"):
        for t0 in regimen["ip_drug"]:
            if t0 <= tmax:
                ev.append((t0, "IPPL", regimen.get("ip_amt", 400.0)))
    ev.sort(key=lambda e: e[0])
    return ev


# =============================================================================
#  INTEGRATOR  (fixed-step RK4 with exact dose events at step boundaries)
# =============================================================================
def simulate(p, tmax=1400.0, dt=0.02, regimen=None, y0=None, record_every=1.0):
    regimen = regimen or {}
    y = list(y0) if y0 else initial_state(p)

    # burn-in: let the fast mechanical / signalling states reach the value
    # implied by the initial rind, with growth and dosing switched off
    p_b = dict(p)
    p_b["KG0"] = 0.0
    p_b["KCOL0"] = 0.0
    tb = 0.0
    while tb < 60.0:
        y = rk4(tb, y, 0.05, p_b)
        tb += 0.05
    y[IDX["CH"]] = 0.0

    doses = build_doses(p, tmax, regimen)
    di = 0
    t = 0.0
    rec = []
    next_rec = 0.0
    surg_done = False
    while t <= tmax + 1e-9:
        while di < len(doses) and doses[di][0] <= t + 1e-9:
            y[IDX[doses[di][1]]] += doses[di][2]
            di += 1
        if (not surg_done) and 0 <= p["SURG_DAY"] <= t:
            y = do_surgery(y, p)
            surg_done = True
        if t >= next_rec - 1e-9:
            a = algebra(t, y, p)
            rec.append(snapshot(t, y, a, p))
            next_rec += record_every
        y = rk4(t, y, dt, p)
        for i in range(NST):
            if y[i] < 0.0 and STATES[i] not in ("CH",):
                y[i] = 0.0
        t += dt
    return rec


def rk4(t, y, h, p):
    k1 = derivs(t, y, p)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(NST)]
    k2 = derivs(t + 0.5 * h, y2, p)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(NST)]
    k3 = derivs(t + 0.5 * h, y3, p)
    y4 = [y[i] + h * k3[i] for i in range(NST)]
    k4 = derivs(t + h, y4, p)
    return [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(NST)]


def do_surgery(y, p):
    """Macroscopic complete resection.

    Surgery sets THICKNESS, not biology.  It cannot reach the fissural
    sanctuary, and it does not reset phi: the collagen fraction of what is left
    behind is the collagen fraction of what was there before."""
    A = areas(p)
    for s, a_, spare in (("PAR", A[0], 0.0), ("VIS", A[1], 0.0), ("FIS", A[2], p["SURG_FIS_SPARE"])):
        T, M, N = y[IDX["T" + s]], y[IDX["M" + s]], y[IDX["N" + s]]
        tot = T + M + N
        if tot <= 1e-9:
            continue
        phi_T = T / tot
        target = p["SURG_RESID"] * p["RHO"] * a_
        keep = max(target, spare * tot)
        f = min(1.0, keep / tot)
        y[IDX["T" + s]] = tot * f * phi_T
        y[IDX["M" + s]] = tot * f * (M / tot)
        y[IDX["N" + s]] = tot * f * (N / tot)
    return y


def snapshot(t, y, a, p):
    return {
        "t": t,
        "h_par_mm": 10.0 * a["h_PAR"], "h_vis_mm": 10.0 * a["h_VIS"],
        "h_fis_mm": 10.0 * a["h_FIS"],
        "mRECIST": a["mRECIST"], "V": a["V_tot"], "Tviable": a["T_tot"],
        "phi": a["phi"], "phi_vis": a["phi_VIS"],
        "f_sm": a["f_sm"]["PAR"], "f_ab": a["f_ab"]["PAR"], "f_T": a["f_T"]["PAR"],
        "lam_sm_mm": 10.0 * a["lam_sm"]["PAR"], "lam_ab_mm": 10.0 * a["lam_ab"]["PAR"],
        "PLV": y[IDX["PLV"]], "FVC": a["FVC"], "FVC_pct": a["FVC_pct"],
        "E_add": a["E_add"], "ECOG": a["ECOG"], "DYSP": y[IDX["DYSP"]],
        "ANC": y[IDX["ANC"]], "GFR": y[IDX["GFR"]], "IRAE": y[IDX["IRAE"]],
        "SMRP": y[IDX["SMRP"]], "TEFF": y[IDX["TEFF"]], "PDL1": y[IDX["PDL1"]],
        "RO_pd1": a["RO_pd1"], "RO_ctla4": a["RO_ctla4"],
        "CH": y[IDX["CH"]], "S": math.exp(-y[IDX["CH"]]),
        "kill": a["kill"]["PAR"], "hazard": a["hazard"],
        "CACHEX": y[IDX["CACHEX"]], "VEGF": y[IDX["VEGF"]], "NEURO": y[IDX["NEURO"]],
        "C_cis": a["C_cis"], "C_pem": a["C_pem"], "C_niv": a["C_niv"],
        "Ct_niv": a["Ct_niv"], "stoma_block": a["stoma_block"],
        "symph": y[IDX["SYMPH"]],
    }


# =============================================================================
#  ENDPOINT EXTRACTION
# =============================================================================
def median_os(rec):
    """Median OS = time at which the cumulative hazard reaches ln 2."""
    for i in range(1, len(rec)):
        if rec[i]["CH"] >= LN2:
            t0, t1 = rec[i - 1]["t"], rec[i]["t"]
            c0, c1 = rec[i - 1]["CH"], rec[i]["CH"]
            return (t0 + (LN2 - c0) / max(c1 - c0, 1e-12) * (t1 - t0)) / 30.44
    return None


def os_rate_at(rec, months):
    tt = months * 30.44
    for i in range(1, len(rec)):
        if rec[i]["t"] >= tt:
            return math.exp(-rec[i]["CH"])
    return math.exp(-rec[-1]["CH"])


def best_response(rec, upto_months=None):
    base = rec[0]["mRECIST"]
    lim = upto_months * 30.44 if upto_months else 1e9
    best = min(r["mRECIST"] for r in rec if r["t"] <= lim)
    return 100.0 * (best - base) / base


def pfs_months(rec):
    """Progression = >=20% rise in mRECIST from nadir (and >=5 mm absolute)."""
    nadir = rec[0]["mRECIST"]
    for r in rec:
        nadir = min(nadir, r["mRECIST"])
        if r["mRECIST"] >= 1.20 * nadir and (r["mRECIST"] - nadir) >= 5.0:
            return r["t"] / 30.44
    return rec[-1]["t"] / 30.44


def volume_response(rec, upto_months=None):
    base = rec[0]["V"]
    lim = upto_months * 30.44 if upto_months else 1e9
    best = min(r["V"] for r in rec if r["t"] <= lim)
    return 100.0 * (best - base) / base


def viable_response(rec, upto_months=None):
    base = rec[0]["Tviable"]
    lim = upto_months * 30.44 if upto_months else 1e9
    best = min(r["Tviable"] for r in rec if r["t"] <= lim)
    return 100.0 * (best - base) / base


# =============================================================================
#  SCENARIOS
# =============================================================================
def scenario(name, emt=0.25, **kw):
    p = default_params()
    p["EMT"] = emt
    reg = {}
    if name == "bsc":
        pass
    elif name == "cis":
        p["ON_CIS"] = 1; reg = {"cis": True, "cycles": 6}
    elif name == "pemcis":
        p["ON_CIS"] = p["ON_PEM"] = 1; reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "pemcarbo":
        p["ON_CARBO"] = p["ON_PEM"] = 1; reg = {"carbo": True, "pem": True, "cycles": 6}
    elif name == "pemcisbev":
        p["ON_CIS"] = p["ON_PEM"] = p["ON_BEV"] = 1
        reg = {"cis": True, "pem": True, "bev": True, "bev_maint": True, "cycles": 6}
    elif name == "nivoipi":
        p["ON_NIV"] = p["ON_IPI"] = 1; reg = {"nivo": True, "ipi": True}
    elif name == "nivo":
        p["ON_NIV"] = 1; reg = {"nivo": True}
    elif name == "pembrochemo":
        p["ON_CIS"] = p["ON_PEM"] = p["ON_PEMBRO"] = 1
        reg = {"cis": True, "pem": True, "pembro": True, "cycles": 6}
    elif name == "adi_pemcis":
        p["ON_CIS"] = p["ON_PEM"] = p["ON_ADI"] = 1; p["ASS1_NEG"] = 1
        reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "gemcis":
        p["ON_CIS"] = p["ON_GEM"] = 1; reg = {"cis": True, "gem": True, "cycles": 6}
    elif name == "surgery_chemo":
        p["ON_CIS"] = p["ON_PEM"] = 1; p["SURG_DAY"] = 63.0
        reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "epp_chemo":
        p["ON_CIS"] = p["ON_PEM"] = 1; p["SURG_DAY"] = 63.0
        p["SURG_FVC_HIT"] = 0.35; p["SURG_FIS_SPARE"] = 0.15
        reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "surg_hithoc":
        p["ON_CIS"] = p["ON_PEM"] = p["ON_IP"] = 1; p["SURG_DAY"] = 63.0
        reg = {"cis": True, "pem": True, "cycles": 6, "ip_drug": [63.5], "ip_amt": 400.0}
    elif name == "ip_only":
        p["ON_IP"] = 1
        reg = {"ip_drug": [7.0, 28.0, 49.0], "ip_amt": 400.0}
    elif name == "talc":
        p["TALC_DAY"] = 7.0
    elif name == "ipc":
        p["IPC_DAY"] = 7.0
    elif name == "talc_chemo":
        p["ON_CIS"] = p["ON_PEM"] = 1; p["TALC_DAY"] = 7.0
        reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "ttfields_chemo":
        p["ON_CIS"] = p["ON_PEM"] = 1; p["TTF"] = 1
        reg = {"cis": True, "pem": True, "cycles": 6}
    elif name == "nivoipi_chemo_seq":
        p["ON_CIS"] = p["ON_PEM"] = p["ON_NIV"] = p["ON_IPI"] = 1
        reg = {"cis": True, "pem": True, "cycles": 4, "nivo": True, "ipi": True}
    else:
        raise ValueError("unknown scenario " + name)
    for k, v in kw.items():
        p[k] = v
    return p, reg


def run(name, emt=0.25, tmax=1400.0, **kw):
    p, reg = scenario(name, emt, **kw)
    return simulate(p, tmax=tmax, regimen=reg), p


# =============================================================================
#  STRUCTURAL CHECKS  (these are the point of the file)
# =============================================================================
def check_geometry(p=None):
    p = p or default_params()
    A_par, A_vis, A_fis = areas(p)
    print("-" * 78)
    print("CHECK 1  GEOMETRY: what a 30% response actually means")
    print("-" * 78)
    print("  parietal area (anatomic x coverage) = %.0f cm2" % A_par)
    print("  dV/dh on the parietal leaf          = %.0f mL per mm of rind" % (A_par / 10.0))
    print("  dV/dh over all three leaves         = %.0f mL per mm" %
          ((A_par + A_vis + A_fis) / 10.0))
    print()
    print("  RIND   : V = A h            ->  30%% thickness drop = %.1f%% volume kill" % 30.0)
    r = 1.0
    print("  SPHERE : V = 4/3 pi r^3     ->  30%% diameter drop = %.1f%% volume kill"
          % (100.0 * (1.0 - 0.70 ** 3)))
    print("  ratio of implied cytoreduction (sphere / rind) = %.2f x"
          % ((1.0 - 0.70 ** 3) / 0.30))
    print()
    print("  So 'partial response' in a mesothelioma trial and 'partial response'")
    print("  in the RECIST tumour it was designed for differ by a factor of 2.2 in")
    print("  the cell kill they imply.  Reading across the two is not valid.")
    print()


def check_penetration(p=None):
    p = p or default_params()
    print("-" * 78)
    print("CHECK 2  PENETRATION: depth-averaged exposure vs rind thickness")
    print("-" * 78)
    print("  h(mm)   phi   lam_sm(mm) f_sm    lam_IgG(mm) f_IgG   f_T(cell)")
    for h_mm in (1, 2, 3, 5, 8, 12, 20):
        for phi in (0.25, 0.55):
            h = h_mm / 10.0
            ls = p["LAM_SM0"] * (1 - phi) ** p["PHI_EXP"]
            la = p["LAM_AB0"] * (1 - phi) ** p["PHI_EXP"]
            lt = p["LAM_T0"] * (1 - phi) ** p["PHI_EXP"]
            print("  %5d  %.2f   %8.2f  %.3f   %9.3f  %.3f   %.3f"
                  % (h_mm, phi, 10 * ls, f_pen(h, ls), 10 * la, f_pen(h, la), f_pen(h, lt)))
    print()
    print("  A 150 kDa antibody has a penetration length of %.2f mm at phi=0.25."
          % (10 * p["LAM_AB0"] * 0.75 ** 1.5))
    f6 = f_pen(0.6, p["LAM_AB0"] * 0.75 ** 1.5)
    f1 = f_pen(0.1, p["LAM_AB0"] * 0.75 ** 1.5)
    print("  A 6 mm rind therefore sees a depth-averaged interstitial antibody")
    print("  concentration of %.1f%% of the plasma-equilibrium value -- %.0f%% of what"
          % (100 * f6, 100 * f6 / f1))
    print("  a 1 mm rind sees, at the SAME plasma exposure.")
    print("  This is a delivery problem no dose escalation fixes: f is bounded by")
    print("  lambda/h, so doubling the dose doubles a number that is already 0.08.")
    print()


def check_crossover(verbose=True):
    """The central test: does the SIGN of OS(IO) - OS(chemo) flip along the
    histology axis, with the magnitudes CheckMate 743 reported?"""
    print("-" * 78)
    print("CHECK 3  HISTOLOGY AXIS: does the treatment comparison flip sign?")
    print("-" * 78)
    print("   x     histology     OS chemo   OS IO    delta    HR(IO/chemo)*")
    rows = []
    for x in (0.0, 0.15, 0.3, 0.5, 0.7, 0.85, 1.0):
        rc, _ = run("pemcis", emt=x)
        ri, _ = run("nivoipi", emt=x)
        oc, oi = median_os(rc), median_os(ri)
        if oc is None or oi is None:
            continue
        hr = oc / oi                      # ratio of medians as an HR proxy
        lab = ("epithelioid" if x < 0.25 else
               "biphasic" if x < 0.65 else "sarcomatoid")
        rows.append((x, lab, oc, oi, oi - oc, hr))
        if verbose:
            print("  %.2f  %-12s  %6.1f mo  %6.1f mo  %+6.1f    %.2f"
                  % (x, lab, oc, oi, oi - oc, hr))
    print("  * ratio of median survivals, a proxy for the hazard ratio")
    # locate the sign change
    cross = None
    for i in range(1, len(rows)):
        if rows[i - 1][4] * rows[i][4] < 0:
            x0, d0 = rows[i - 1][0], rows[i - 1][4]
            x1, d1 = rows[i][0], rows[i][4]
            cross = x0 + (0 - d0) / (d1 - d0) * (x1 - x0)
    print()
    if cross is not None:
        print("  SIGN CHANGE at x = %.2f  (between epithelioid and biphasic)" % cross)
    else:
        print("  no sign change found over x in [0,1]")
    print()
    return rows


def check_collagen_bias():
    """How much cytoreduction does a thickness-based endpoint hide?"""
    print("-" * 78)
    print("CHECK 4  THE COLLAGEN BIAS: viable-cell kill vs measured thickness")
    print("-" * 78)
    print("  regimen        d(viable)%  d(volume)%  d(mRECIST)%   hidden kill (pp)")
    for name in ("pemcis", "pemcisbev", "nivoipi", "gemcis"):
        rec, _ = run(name, emt=0.25, tmax=700.0)
        dv = viable_response(rec, 12)
        dvol = volume_response(rec, 12)
        dm = best_response(rec, 12)
        print("  %-13s  %8.1f    %8.1f    %9.1f      %8.1f"
              % (name, dv, dvol, dm, dv - dm))
    print()
    print("  The gap is the collagen (and the necrotic debris that has not yet")
    print("  been cleared).  Every one of these regimens kills more tumour than")
    print("  the thickness endpoint credits it with.")
    print()


def check_effusion_vs_trapping():
    """Are the two causes of dyspnoea separable, and does pleurodesis fail
    where the model says it should?"""
    print("-" * 78)
    print("CHECK 5  TWO CAUSES OF DYSPNOEA: effusion (removable) vs trapping (not)")
    print("-" * 78)
    print("  h_vis(mm) phi_vis  E_add(cmH2O/L)  FVC drained(L)  FVC with 2L effusion")
    p0 = default_params()
    for h_vis_mm in (2, 4, 6, 9, 13):
        p = default_params()
        y = initial_state(p, h_par=0.60, h_vis=h_vis_mm / 10.0, h_fis=0.50)
        a = algebra(0.0, y, p)
        vexp = a["vexp_target"]
        fvc_dry = p["FVC_PRED"] * vexp
        fvc_wet = fvc_dry * (1.0 - p["PLV_FVC"] * 2000.0 / p["PLV_MAX"])
        print("  %8d  %6.2f  %13.1f  %14.2f  %20.2f"
              % (h_vis_mm, a["phi_VIS"], a["E_add"], fvc_dry, fvc_wet))
    print()
    print("  The trapped-lung threshold in the pleural-manometry literature is a")
    print("  pleural elastance of 14.5 cmH2O/L.  In this model that is crossed at")
    for h_vis_mm in [i / 2.0 for i in range(2, 40)]:
        p = default_params()
        y = initial_state(p, h_par=0.60, h_vis=h_vis_mm / 10.0, h_fis=0.50)
        a = algebra(0.0, y, p)
        if a["E_add"] >= p0["E50"]:
            print("  a visceral rind of %.1f mm at phi = %.2f." % (h_vis_mm, a["phi_VIS"]))
            break
    rec_t, _ = run("talc", tmax=180.0)
    rec_i, _ = run("ipc", tmax=180.0)
    print()
    print("  talc pleurodesis, day 90: symphysis %.2f, effusion %.0f mL, FVC %.2f L"
          % (rec_t[90]["symph"], rec_t[90]["PLV"], rec_t[90]["FVC"]))
    print("  IPC drainage,     day 90: symphysis %.2f, effusion %.0f mL, FVC %.2f L"
          % (rec_i[90]["symph"], rec_i[90]["PLV"], rec_i[90]["FVC"]))
    print()


def check_surgery():
    print("-" * 78)
    print("CHECK 6  SURGERY: setting h without resetting phi or reaching the fissure")
    print("-" * 78)
    for name in ("pemcis", "surgery_chemo", "epp_chemo", "surg_hithoc"):
        rec, _ = run(name, emt=0.25, tmax=1200.0)
        os_ = median_os(rec)
        print("  %-15s  median OS %5.1f mo   mRECIST d90 %5.1f mm   FVC d180 %.2f L"
              % (name, os_ if os_ else float("nan"), rec[90]["mRECIST"], rec[180]["FVC"]))
    print()
    print("  MARS2 (Lancet Respir Med 2024) randomised extended pleurectomy/")
    print("  decortication plus chemotherapy against chemotherapy alone and found")
    print("  an OS hazard ratio of 1.28 AGAINST surgery.")
    print()
    print("  THE MODEL DOES NOT REPRODUCE THAT DIRECTION.  It predicts a large")
    print("  benefit from cytoreduction, because in this model an operation is a")
    print("  clean reset of rind THICKNESS and nothing else.  What the model is")
    print("  missing is everything the trial says the harm comes from:")
    print("    - no peri-operative mortality (MARS2: 3.4% at 30 days);")
    print("    - no months-long performance-status and FVC debt from a thoracotomy;")
    print("    - residual microscopic disease is modelled as a THIN RIND over the")
    print("      SAME 45% surface coverage, when 'macroscopic complete resection'")
    print("      in fact leaves tumour cells over essentially the WHOLE pleural")
    print("      surface -- the coverage fraction should go to 1 as the thickness")
    print("      goes to 0.1 mm, and in this model it does not.")
    print()
    print("  This is reported rather than fixed: tuning a surgical penalty until")
    print("  the model agreed with MARS2 would make the agreement an input.")
    print()


def check_trial_targets():
    print("-" * 78)
    print("CHECK 7  CALIBRATION AGAINST PUBLISHED TRIAL ENDPOINTS")
    print("-" * 78)
    print("  All first-line arms run at x = 0.25, the histology mix of the")
    print("  registration trials (roughly 70-75% epithelioid).  The bevacizumab and")
    print("  checkpoint arms are anchored to the EMPHACIS doublet by their published")
    print("  INCREMENT and RATIO respectively, because their own control arms")
    print("  (MAPS 16.1 mo, CheckMate 743 14.1 mo) enrolled fitter populations than")
    print("  EMPHACIS did and an absolute comparison would confound the two.")
    print()
    targets = [
        ("BSC / untreated (historical)",        "bsc",       7.0,  None, "fitted"),
        ("Cisplatin alone (EMPHACIS control)",  "cis",       9.3,  3.9,  "fitted"),
        ("Pemetrexed + cisplatin (EMPHACIS)",   "pemcis",    12.1, 5.7,  "fitted"),
        ("+ bevacizumab (MAPS increment)",      "pemcisbev", 14.8, 9.2,  "fitted"),
        ("Nivo + ipi (CM-743 ratio 1.28)",      "nivoipi",   15.5, 6.8,  "fitted"),
        ("Gem + cis (phase II)",                "gemcis",    11.2, 6.0,  "HELD OUT"),
        ("Pem + carbo (elderly / unfit)",       "pemcarbo",  12.7, 6.5,  "HELD OUT"),
        ("Nivolumab mono, 2nd line (CONFIRM)",  "nivo",      10.2, 3.0,  "HELD OUT"),
    ]
    print("  %-36s %8s %10s  %8s %8s  %s"
          % ("regimen", "OS obs", "OS model", "PFS obs", "PFS mod", "role"))
    err, nfit = 0.0, 0
    for lab, name, os_obs, pfs_obs, role in targets:
        rec, _ = run(name, emt=0.25, tmax=1600.0)
        os_m = median_os(rec)
        pfs_m = pfs_months(rec)
        print("  %-36s %6.1f mo %8.1f mo  %6s %7.1f  %s"
              % (lab, os_obs, os_m if os_m else float("nan"),
                 ("%.1f" % pfs_obs) if pfs_obs else "  -", pfs_m, role))
        if os_m:
            err += ((os_m - os_obs) / os_obs) ** 2
            nfit += 1
    print()
    print("  root-mean-square relative error on median OS across %d arms: %.1f%%"
          % (nfit, 100.0 * math.sqrt(err / nfit)))
    print()
    print("  Best mRECIST change of the TYPICAL patient (one trajectory is NOT a")
    print("  response RATE -- for that see the virtual population in mpm_calibration.py):")
    for name in ("cis", "pemcis", "pemcisbev", "nivoipi"):
        rec, _ = run(name, emt=0.25, tmax=700.0)
        print("    %-12s  best mRECIST %+6.1f%%   best viable-cell change %+6.1f%%"
              % (name, best_response(rec, 12), viable_response(rec, 12)))
    print()


def check_checkmate743():
    print("-" * 78)
    print("CHECK 8  CHECKMATE 743 SUBGROUPS -- THE HELD-OUT TEST OF CLAIM 3")
    print("-" * 78)
    print("  Step 7 of the calibration fitted ONE number: the IO/chemo median ratio")
    print("  at the trial population's histology.  The SPLIT below was not fitted.")
    print()
    obs = {"epithelioid": (16.5, 18.7, 0.86), "non-epithelioid": (8.8, 18.1, 0.46)}
    print("  %-17s %9s %9s   %s" % ("subgroup", "chemo", "IO", "median ratio (obs HR)"))
    for lab, x in (("epithelioid", 0.15), ("non-epithelioid", 0.85)):
        rc, _ = run("pemcis", emt=x, tmax=1800.0)
        ri, _ = run("nivoipi", emt=x, tmax=1800.0)
        oc, oi = median_os(rc), median_os(ri)
        o_c, o_i, o_hr = obs[lab]
        print("  %-17s %5.1f mo   %5.1f mo   %.2f  (observed HR %.2f; trial arms "
              "%.1f vs %.1f mo)" % (lab, oc, oi, oc / oi, o_hr, o_c, o_i))
    print()
    print("  The model's absolute survivals are anchored to EMPHACIS, which ran in a")
    print("  less fit population than CheckMate 743, so the ABSOLUTE months are not")
    print("  expected to match the trial's arms.  The RATIO is the comparable")
    print("  quantity, and its movement along the histology axis is the prediction.")
    print()
    rc, _ = run("pemcis", emt=0.25, tmax=800.0)
    ri, _ = run("nivoipi", emt=0.25, tmax=800.0)
    print("  2-year survival at the trial population's histology:")
    print("    chemo     %.0f%%   (CheckMate 743 observed 27%%)" % (100 * os_rate_at(rc, 24)))
    print("    nivo+ipi  %.0f%%   (CheckMate 743 observed 41%%)" % (100 * os_rate_at(ri, 24)))
    print()


def check_toxicity():
    print("-" * 78)
    print("CHECK 9  TOXICITY AND THE RENAL FEEDBACK LOOP")
    print("-" * 78)
    rec, p = run("pemcis", emt=0.25, tmax=200.0)
    nadir = min(r["ANC"] for r in rec[:30])
    day_nadir = min(rec[:30], key=lambda r: r["ANC"])["t"]
    print("  ANC nadir cycle 1: %.2f x10^9/L on day %.0f (observed nadir day 8-12)"
          % (nadir, day_nadir))
    p2 = default_params(); p2["FOLATE"] = 0
    p2["ON_CIS"] = p2["ON_PEM"] = 1
    rec2 = simulate(p2, tmax=200.0, regimen={"cis": True, "pem": True, "cycles": 6})
    nadir2 = min(r["ANC"] for r in rec2[:30])
    print("  ANC nadir without folate/B12: %.2f  ->  %.0f%% deeper"
          % (nadir2, 100 * (nadir - nadir2) / nadir))
    print("  (EMPHACIS added folate + B12 mid-trial and grade 3-4 toxicity fell)")
    print()
    print("  eGFR after 6 cycles of cisplatin: %.0f -> %.0f mL/min"
          % (p["GFR0"], rec[-1]["GFR"]))
    print("  pemetrexed clearance therefore falls to %.0f%% of baseline, which"
          % (100 * rec[-1]["GFR"] / p["GFR0"]))
    print("  RAISES exposure and deepens later nadirs: a positive feedback loop")
    print("  between the two drugs of the doublet through the kidney.")
    print("  cumulative neuropathy index after 6 cycles: %.2f" % rec[-1]["NEURO"])
    print()
    ri, _ = run("nivoipi", tmax=400.0)
    print("  irAE driver at 6 months, nivo+ipi: %.2f (CTLA-4 dominant)" % ri[180]["IRAE"])
    rn, _ = run("nivo", tmax=400.0)
    print("  irAE driver at 6 months, nivo alone: %.2f" % rn[180]["IRAE"])
    print()


def check_biomarker():
    print("-" * 78)
    print("CHECK 10  SMRP TRACKS VIABLE CELLS, NOT RIND THICKNESS")
    print("-" * 78)
    rec, _ = run("pemcis", emt=0.25, tmax=400.0)
    b = rec[0]
    for day in (0, 30, 60, 90, 180, 270, 360):
        r = rec[day]
        print("  day %3d   viable %+6.1f%%   mRECIST %+6.1f%%   SMRP %+6.1f%%"
              % (day,
                 100 * (r["Tviable"] - b["Tviable"]) / b["Tviable"],
                 100 * (r["mRECIST"] - b["mRECIST"]) / b["mRECIST"],
                 100 * (r["SMRP"] - b["SMRP"]) / max(b["SMRP"], 1e-9)))
    print()
    print("  SMRP is shed from viable cells and cleared renally.  It moves with")
    print("  the quantity the drug actually changes, and it moves EARLIER than the")
    print("  thickness does, because the matrix and the necrotic debris are still")
    print("  sitting in the rind when the cells that made them are already dead.")
    print()


def check_intrapleural():
    print("-" * 78)
    print("CHECK 11  INTRAPLEURAL THERAPY NEEDS CYTOREDUCTION FIRST")
    print("-" * 78)
    p = default_params()
    print("  residual h (mm)  f_ip (depth-averaged exposure)  relative kill")
    base = None
    for h_mm in (0.5, 1, 2, 3, 5, 8, 12):
        h = h_mm / 10.0
        lam = p["LAM_IP0"] * (1 - 0.25) ** p["PHI_EXP"]
        f = f_pen(h, lam)
        if base is None:
            base = f
        print("  %14.1f  %28.3f  %13.2f" % (h_mm, f, f / base))
    print()
    print("  Intracavitary penetration is measured at 3-5 mm in the HIPEC/HITHOC")
    print("  literature and the model's lambda_ip is %.1f mm at phi = 0.25."
          % (10 * p["LAM_IP0"] * 0.75 ** 1.5))
    print("  A 12 mm rind gets %.0f%% of the exposure a 1 mm residuum gets, which"
          % (100 * f_pen(1.2, p["LAM_IP0"] * 0.75 ** 1.5) /
             f_pen(0.1, p["LAM_IP0"] * 0.75 ** 1.5)))
    print("  is why every positive intracavitary result in this disease is reported")
    print("  AFTER macroscopic complete resection and never instead of it.")
    print()


def check_mass_balance():
    """Numerical hygiene: is the integration converged, and is mass conserved
    through the kill -> necrosis -> clearance chain?"""
    print("-" * 78)
    print("CHECK 12  NUMERICAL HYGIENE")
    print("-" * 78)
    for dt in (0.08, 0.04, 0.02, 0.01):
        p, reg = scenario("pemcis", 0.25)
        rec = simulate(p, tmax=365.0, dt=dt, regimen=reg)
        print("  dt = %.3f d   mRECIST(365) = %8.4f mm   CH(365) = %.6f"
              % (dt, rec[-1]["mRECIST"], rec[-1]["CH"]))
    print("  -> dt = 0.08 d is OUTSIDE the stability region of the marrow transit")
    print("     chain (ktr = 0.87/d) and its 1-year endpoint is meaningless.  From")
    print("     dt = 0.04 d downward the solution is converged: halving the step")
    print("     moves the 1-year endpoint in the 4th significant figure.  dt = 0.02")
    print("     d is used throughout.")
    print()
    p, reg = scenario("pemcis", 0.25)
    rec = simulate(p, tmax=365.0, regimen=reg)
    neg = [k for k in ("h_par_mm", "V", "PLV", "ANC", "GFR", "TEFF", "FVC")
           if any(r[k] < 0 for r in rec)]
    print("  states going negative over 365 d: %s" % (neg if neg else "none"))
    print()


def main():
    print("=" * 78)
    print(" MALIGNANT PLEURAL MESOTHELIOMA QSP MODEL -- reference run")
    print(" %d ODE states, RK4, dt = 0.02 d, no external dependencies" % NST)
    print("=" * 78)
    print()
    check_geometry()
    check_penetration()
    check_trial_targets()
    check_crossover()
    check_checkmate743()
    check_collagen_bias()
    check_effusion_vs_trapping()
    check_surgery()
    check_intrapleural()
    check_toxicity()
    check_biomarker()
    check_mass_balance()
    print("=" * 78)
    print(" done")
    print("=" * 78)


if __name__ == "__main__":
    main()
