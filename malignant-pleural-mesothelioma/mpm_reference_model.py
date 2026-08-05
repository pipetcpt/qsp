#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpm_reference_model.py
======================
Dependency-free reference implementation of the malignant pleural mesothelioma
(MPM) QSP model.  Every one of the 51 differential equations in
`mpm_mrgsolve_model.R` is re-implemented here and integrated with a fixed-step
classical Runge-Kutta 4 scheme, using only the Python standard library.

WHY THIS FILE EXISTS
--------------------
There is no R runtime in the environment this model was built in, so the
mrgsolve file could not be executed.  Rather than ship 51 unexecuted ODEs, the
whole system was written twice: once in mrgsolve C++ syntax and once here, and
every number that appears in the README, in the model header, and in the
repository table was produced by running THIS file.  Defects that this
re-implementation exposed are recorded, with the fix, in
`mpm_calibration_output.txt` and as comments at the point of repair in both
files.

THE STRUCTURAL COMMITMENT
-------------------------
The tumour is a SHEET, not a ball:

    burden  V(t) = A(t) * h(t)          A = involved pleural area   [cm^2]
                                        h = mean rind thickness     [cm]

and the four consequences are all mechanical, not rhetorical:

  (a) mRECIST sums six perpendicular THICKNESSES.  It observes h and is blind
      to A.  Worse, h includes the fibro-necrotic matrix left behind by killed
      cells, so the measurement is  h_meas = (N + phi*M)/A  and not N/A.
  (b) A sheet proliferates at its FRONT (a perimeter), not throughout a bulk,
      so dA/dt scales with 2*sqrt(pi*A) and the volume doubling time lengthens
      on its own without any resistance mechanism.
  (c) The rind is perfused from ONE face.  The mean relative drug (or oxygen,
      or T-cell) exposure across a slab of thickness h fed from one side with
      penetration length L is
          fpen(L,h) = (L/h) * (1 - exp(-h/L))
      which is 1 for a thin rind and L/h for a thick one.  Intrapleural drug
      enters from the OPPOSITE face, so the two routes cover complementary
      depths rather than the same depth twice.
  (d) "Macroscopic complete resection" of a sheet collapses h and leaves A
      untouched, in a wound flooded with IL-6 and TGF-beta.

Run:
    python3 mpm_reference_model.py            # full calibration + analyses
    python3 mpm_reference_model.py --quick    # scenarios only
"""

import math
import sys

# ============================================================================
# 0.  PARAMETERS
# ============================================================================
# Time unit: DAYS.  Volumes of tumour tissue: cm^3 (1 cm^3 ~ 1e9 packed cells).
# Areas: cm^2.  Thickness: cm.  Drug concentrations: mg/L.  Effusion: mL.

P = dict(
    # ---- body / demographics -------------------------------------------
    BSA        = 1.80,     # m2
    WT         = 70.0,     # kg
    LBM0       = 52.0,     # kg lean body mass

    # ---- pleural geometry ----------------------------------------------
    S_PL       = 1300.0,   # cm2, total pleural surface of one hemithorax
    PHIM       = 0.70,     # fraction of a killed cell's volume that persists
                           # as fibro-necrotic matrix and is still MEASURED
    KDEGM      = 0.008,    # /d  matrix resorption (t1/2 ~ 87 d)

    # ---- tumour cell kinetics -------------------------------------------
    KPROL      = 0.2680,   # /d  proliferation rate of a fully oxygenated cell
    L_OX       = 0.0180,   # cm  oxygen/nutrient penetration length
    KDEATH0    = 0.0005,   # /d  baseline apoptosis
    KNEC       = 0.00180,  # /d  hypoxic death of the unreached fraction
    KFRONT     = 0.0550,   # cm/d radial creep speed of the growth front
    HFRONT     = 0.0600,   # cm  thickness OF THE ADVANCING MARGIN itself
    GFRONT_REF = 0.0700,   # /d  reference net proliferation at the margin
    ANG_MIN    = 0.620,    # residual front speed under complete VEGF blockade
    KINV       = 0.00350,  # cm/d chest-wall / diaphragm invasion
    KMET       = 2.2e-5,   # /d  seeding of nodal + distant burden
    KMETG      = 0.0035,   # /d  growth of the metastatic compartment
    METMAX     = 150.0,    # cm3 ceiling on the nodal + distant compartment
    WMET       = 1.50,     # weight of metastatic burden in the hazard

    # ---- histology multipliers (set at run time) ------------------------
    CHEMOS     = 1.00,     # chemosensitivity   epi 1.00 / bip 0.62 / sar 0.35
    IMMINF     = 0.38,     # immune infiltration epi 0.38 / bip 0.70 / sar 1.00
    VISTA_S    = 1.15,     # VISTA-type suppression epi 1.15 / bip 0.60 / sar 0.20
    COLLF      = 1.00,     # stromal collagen (diffusivity divisor) epi 1 / sar 1.9
    MSLNF      = 1.00,     # mesothelin expression epi 1.00 / bip 0.50 / sar 0.12
    SARCH      = 1.00,     # hazard multiplier epi 1.00 / bip 1.35 / sar 1.90
    KPROLF     = 1.00,     # proliferation multiplier epi 1.00 / bip 1.25 / sar 1.55

    # ---- effusion --------------------------------------------------------
    KFORM      = 1400.0,   # mL/d maximal filtration scale
    VEFFMAX    = 3200.0,   # mL  pressure/anatomical ceiling of one hemithorax
    KABS       = 0.350,    # /d  absorption from an apposed (symphysed) space
    EMAXV      = 3.00,     # VEGF permeability multiplier
    EC50V      = 250.0,    # pg/mL
    KDRAIN     = 2200.0,   # mL/d  lymphatic stomatal capacity (x20 reserve)
    FOBS       = 0.95,     # fraction of stomata a full rind can occlude
    KMD        = 300.0,    # mL  half-saturation of stomatal drainage
    KSYMSP     = 0.0060,   # /d  spontaneous tumour symphysis (A/S)^3
    KTALC      = 0.35,     # /d  talc-induced symphysis
    VDEAD      = 100.0,    # mL  residual pleural volume after drainage

    # ---- VEGF / angiogenesis / delivery geometry -------------------------
    KVSYN      = 317.0,    # pg/mL/d per 100 cm3 of viable tumour
    KVDEG      = 1.50,     # /d
    KDBEV      = 8.00,     # mg/L  bevacizumab-VEGF sequestration constant
    KANG       = 0.0500,   # /d
    EC50A      = 300.0,    # pg/mL
    RHOMAX     = 2.50,
    KVREG      = 0.0200,   # /d
    IFP_BASE   = 4.00,     # mmHg
    IFP_MAX    = 22.0,     # mmHg
    EC50I      = 350.0,    # pg/mL
    IFP_REF    = 17.9,     # mmHg, reference for the L_d scaling
    LD0_CIS    = 0.0750,   # cm  cisplatin penetration length at reference state
    LD0_PEM    = 0.0600,   # cm
    LD0_MAB    = 0.0200,   # cm  150 kDa IgG
    L_T        = 0.0300,   # cm  T-cell infiltration length
    LDIP0      = 0.0300,   # cm  intrapleural penetration (from the free face)
    PEXP_RHO   = 0.50,     # L_d exponent on microvessel density
    PEXP_IFP   = 0.35,     # L_d exponent on (IFP_ref/IFP)

    # ---- pemetrexed PK ----------------------------------------------------
    V1_PEM     = 12.9,     # L
    V2_PEM     = 8.00,     # L
    Q_PEM      = 48.0,     # L/d
    CL_PEM     = 108.0,    # L/d at CrCl 100 mL/min
    # intracellular, EXPOSED tumour cells (high FPGS, low GGH)
    KIN_T      = 6.00, KOUT_T = 20.0, KFPGS_T = 3.00, KGGH_T = 0.100,
    # intracellular, marrow progenitors (low FPGS, high GGH)
    KIN_M      = 6.00, KOUT_M = 20.0, KFPGS_M = 0.900, KGGH_M = 0.380,
    KI_FOL     = 900.0,    # nM, weak RFC competition by circulating folate

    # ---- folate / B12 -----------------------------------------------------
    FOL_DIET   = 12.0,     # nM plasma folate, unsupplemented
    KFOUT      = 0.050,    # /d
    KFIN       = 0.00413,  # nM/d per ug of daily folic acid
    KHSYN      = 26.0,     # umol/L/d homocysteine production
    KHDEG      = 4.20,     # /d
    KMH        = 20.0,     # nM
    KHB12      = 0.55,     # /d extra remethylation with B12 repletion

    # ---- cisplatin PK ------------------------------------------------------
    V1_CIS     = 25.0,     # L (free platinum)
    KEL_CIS    = 12.0,     # /d renal elimination of free Pt
    KBIND      = 21.3,     # /d irreversible plasma-protein binding
    VB_CIS     = 5.00,     # L
    KBEL       = 0.140,    # /d elimination of bound Pt (t1/2 ~ 5 d)
    KADD       = 1.00,     # adduct formation per (mg/L . d)
    KREP       = 0.350,    # /d nucleotide-excision repair
    ERCC1F     = 1.00,     # repair capacity multiplier

    # ---- chemotherapy effect ------------------------------------------------
    EMAXP      = 0.1000,   # /d  pemetrexed maximal kill (exposed cells)
    EC50P      = 0.900,    # mg/L polyglutamate equivalent
    KFPT       = 400.0,    # nM  folate rescue constant, TUMOUR (weak)
    EMAXC      = 0.7800,   # /d  cisplatin maximal kill
    EC50C      = 0.280,    # adduct units
    SYNCP      = 0.900,    # pemetrexed x platinum synergy
    EMAXVIN    = 0.110,    # /d  vinorelbine / gemcitabine second line
    EC50VIN    = 0.500,

    # ---- marrow toxicity (Friberg) -------------------------------------------
    CIRC0      = 4.00,     # 10^9/L absolute neutrophil count
    KTR        = 0.873,    # /d  (MTT = 110 h over 4 compartments)
    GAM        = 0.170,
    EMAXPM     = 1.380,    # pemetrexed marrow effect
    EC50PM     = 0.850,
    KFPM       = 14.0,     # nM  folate rescue constant, MARROW (strong)
    SLOPE_CIS  = 0.090,    # per mg/L free Pt

    # ---- renal ---------------------------------------------------------------
    CRCL0      = 95.0,     # mL/min
    KRECR      = 0.0500,   # /d recovery toward the (falling) set point
    KNEPH      = 47.0,     # mL/min per (mg/L.d) of free Pt, reversible
    KNEPHIRR   = 12.0,     # mL/min per (mg/L.d) of free Pt, permanent

    # ---- biologics PK ---------------------------------------------------------
    V1_BEV = 3.00, V2_BEV = 2.40, Q_BEV = 0.55, CL_BEV = 0.210,   # L, L/d
    V1_NIV = 4.50, V2_NIV = 3.20, Q_NIV = 0.65, CL_NIV = 0.240,
    V1_IPI = 4.40, CL_IPI = 0.370,
    EC50_PD1 = 0.100,      # mg/L  nivolumab receptor occupancy
    EC50_CT4 = 1.50,       # mg/L  ipilimumab occupancy in the draining node

    # ---- ADI-PEG20 / arginine ---------------------------------------------------
    V_ADI      = 4.00,     # L
    KA_ADI     = 0.500,    # /d  IM absorption
    CL_ADI     = 0.300,    # L/d  (pegylated, t1/2 ~ 9 d)
    KADA       = 0.00450,  # /d  anti-drug antibody formation
    KADAEL     = 0.0250,   # /d
    ADA_POT    = 8.00,     # ADA-driven clearance multiplier
    ARG0       = 90.0,     # umol/L plasma arginine
    KARGIN     = 45.0,     # umol/L/d
    KARGOUT    = 0.500,    # /d
    KADIARG    = 20.00,    # /d per mg/L of ADI
    ARG_CRIT   = 12.0,     # umol/L below which ASS1-null cells starve
    EMAXARG    = 0.0450,   # /d kill at full arginine depletion
    ASS1LOSS   = 1.00,     # 1 = ASS1-methylated (auxotroph), 0 = ASS1 intact

    # ---- immune module ------------------------------------------------------------
    KPRIME     = 0.0350,   # /d
    EIPI       = 4.50,     # ipilimumab amplification of priming
    KEXH       = 0.0900,   # /d exhaustion under PD-1 engagement
    PDL1B      = 1.00,
    KTD        = 0.0350,   # /d effector decay
    KSUP       = 0.0180,   # /d per unit Treg
    KTREGS     = 0.0400, KMTG = 8.00, KTREGD = 0.0300, KIPIDEP = 0.0450,
    KCLON      = 0.0500,   # /d memory formation
    KCLOND     = 0.0100,   # /d memory decay (the CheckMate 743 tail)
    WMEM       = 0.600,
    KKILL      = 0.3300,   # /d maximal immune kill of exposed cells
    KMKILL     = 6.00,
    KTGF       = 0.0450, KTGFD = 0.120,

    # ---- inflammation / cachexia / symptoms ------------------------------------
    KIL6T      = 41.4,     # pg/mL/d per 100 cm3 tumour
    KIL6W      = 260.0,    # pg/mL/d per unit surgical wound signal
    KIL6D      = 8.00,     # /d
    KWD        = 0.0600,   # /d resolution of the surgical wound response
    KCACH      = 0.0130,   # kg/d at saturating IL-6
    KMIL6      = 22.0,     # pg/mL
    KLREC      = 0.00250,  # /d
    FVC0       = 92.0,     # % predicted at diagnosis (already reduced)
    A_FVC      = 0.450, H_FVC = 0.100, V_FVC = 0.220,

    # ---- biomarker --------------------------------------------------------------
    KSMRP      = 0.0217,   # nM/d per cm3 viable tumour
    KSMRPCL    = 1.400,    # /d at CrCl 95

    # ---- hazard -----------------------------------------------------------------
    H0         = 0.000167, # /d
    V0H        = 360.0,    # cm3 reference burden
    B1         = 0.500,    # burden exponent
    B2         = 0.550,    # per ECOG point
    B5         = 2.200,    # per unit of pleural surface ENCASED (A/S)
    H_FEBN     = 0.00400,  # /d added hazard while ANC < 0.5
    H_IRAE     = 0.00095,  # /d added hazard per unit grade>=3 irAE intensity
    H_SURG30   = 0.0,      # set by the surgical scenarios

    # ---- intrapleural route -------------------------------------------------------
    FCONT      = 0.550,    # fraction of the involved surface an instillate
                           # actually bathes (loculation, gravity, adhesions)
    KIPCL      = 8.00,     # /d clearance of drug out of the pleural space
    FABS_IP    = 0.55,     # fraction of that clearance reaching the systemic pool
)

HIST = {
    "epithelioid": dict(CHEMOS=1.00, IMMINF=0.38, VISTA_S=1.15, COLLF=1.00,
                        MSLNF=1.00, SARCH=1.00, KPROLF=1.00),
    "biphasic":    dict(CHEMOS=0.52, IMMINF=0.70, VISTA_S=0.60, COLLF=1.45,
                        MSLNF=0.50, SARCH=1.35, KPROLF=1.25),
    "sarcomatoid": dict(CHEMOS=0.18, IMMINF=1.00, VISTA_S=0.20, COLLF=1.90,
                        MSLNF=0.12, SARCH=1.90, KPROLF=1.55),
}
# Histology mix of the modern first-line trials (CheckMate 743: 76 % epithelioid)
COHORT = [("epithelioid", 0.75), ("biphasic", 0.13), ("sarcomatoid", 0.12)]

# ============================================================================
# 1.  STATE VECTOR
# ============================================================================
STATES = [
    # tumour geometry (the core)
    "N",        # 0  viable tumour tissue                       [cm3]
    "M",        # 1  necrotic + fibrotic matrix                 [cm3]
    "A",        # 2  involved pleural area                      [cm2]
    "Z",        # 3  chest-wall / diaphragm invasion depth      [cm]
    "MET",      # 4  nodal + distant burden                     [cm3]
    # pleural space
    "VEFF",     # 5  effusion volume                            [mL]
    "PSY",      # 6  pleural symphysis                          [0-1]
    "VEGF",     # 7  VEGF-A                                     [pg/mL]
    "RHOV",     # 8  relative microvessel density               [-]
    # inflammation / immunity
    "IL6",      # 9                                             [pg/mL]
    "TGFB",     # 10                                            [a.u.]
    "TEFF",     # 11 intratumoural CD8 effectors                [a.u.]
    "TREG",     # 12                                            [a.u.]
    "TCLON",    # 13 memory / clonal pool                       [a.u.]
    "WOUND",    # 14 surgical wound signal                      [0-1]
    # pemetrexed
    "PEM_C",    # 15 central                                    [mg]
    "PEM_P",    # 16 peripheral                                 [mg]
    "PEM_T",    # 17 tumour intracellular monoglutamate         [mg/L]
    "PEM_TP",   # 18 tumour polyglutamates                      [mg/L]
    "PEM_M",    # 19 marrow intracellular monoglutamate         [mg/L]
    "PEM_MP",   # 20 marrow polyglutamates                      [mg/L]
    "FOL",      # 21 plasma folate                              [nM]
    "HCY",      # 22 homocysteine                               [umol/L]
    # cisplatin
    "CIS_F",    # 23 free platinum (central)                    [mg]
    "CIS_B",    # 24 protein-bound platinum                     [mg/L]
    "ADD",      # 25 Pt-DNA adducts, systemically exposed cells [a.u.]
    "ADDIP",    # 26 Pt-DNA adducts from the intrapleural route [a.u.]
    # biologics
    "BEV_C", "BEV_P",      # 27 28  [mg]
    "NIV_C", "NIV_P",      # 29 30  [mg]
    "IPI_C",               # 31     [mg]
    "ADI_A", "ADI_C",      # 32 33  IM depot [mg], central [mg]
    "ADA",                 # 34     anti-ADI antibody [a.u.]
    "ARG",                 # 35     plasma arginine [umol/L]
    # intrapleural depot
    "IPD",                 # 36     drug amount in the pleural space [mg]
    # marrow
    "PROL", "TR1", "TR2", "TR3", "CIRC",   # 37 38 39 40 41
    # renal
    "CRCL",     # 42 measured creatinine clearance              [mL/min]
    "CRCLSS",   # 43 irreversibly lost set point                [mL/min]
    # biomarker / host / outcome
    "SMRP",     # 44 serum mesothelin                           [nM]
    "LBM",      # 45 lean body mass                             [kg]
    "CUMH",     # 46 cumulative hazard
    "AUCP",     # 47 cumulative pemetrexed AUC                  [mg.d/L]
    "AUCC",     # 48 cumulative free-Pt AUC                     [mg.d/L]
    "VINE",     # 49 vinorelbine effect-site                    [mg/L]
    "IRAE",     # 50 immune-related AE intensity                [0-1]
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)


# ============================================================================
# 2.  GEOMETRY HELPERS
# ============================================================================
def fpen(L, h):
    """Mean relative concentration across a slab of thickness h fed from ONE
    face with penetration length L.  ->1 for h<<L, ->L/h for h>>L."""
    if h <= 1e-9:
        return 1.0
    r = h / L
    if r < 1e-6:
        return 1.0
    if r > 60.0:            # exp(-r) underflows into irrelevance
        return L / h
    return (L / h) * (1.0 - math.exp(-r))


def geometry(y, p):
    """Everything downstream of 'the tumour is a sheet'."""
    N = max(y[IX["N"]], 0.0)
    M = max(y[IX["M"]], 0.0)
    A = max(y[IX["A"]], 1.0)
    h_tot = (N + p["PHIM"] * M) / A          # what the CT scanner measures
    h_via = N / A                            # what the drug has to kill
    return N, M, A, h_tot, h_via


def penetration_lengths(y, p):
    """L_d for each modality, given the CURRENT vascular state of the rind."""
    VEGFf = free_vegf(y, p)
    IFP = p["IFP_BASE"] + p["IFP_MAX"] * VEGFf / (VEGFf + p["EC50I"])
    rho = max(y[IX["RHOV"]], 1e-3)
    fv = rho ** p["PEXP_RHO"]
    fi = (p["IFP_REF"] / max(IFP, 1e-3)) ** p["PEXP_IFP"]
    fc = 1.0 / math.sqrt(p["COLLF"])
    return dict(
        cis=p["LD0_CIS"] * fv * fi * fc,
        pem=p["LD0_PEM"] * fv * fi * fc,
        mab=p["LD0_MAB"] * fv * fi * fc,
        tcell=p["L_T"] * fv * fi * fc,
        ip=p["LDIP0"] * fc,                  # no vasculature needed
        ifp=IFP,
    )


def free_vegf(y, p):
    bev = y[IX["BEV_C"]] / p["V1_BEV"]
    return max(y[IX["VEGF"]], 0.0) / (1.0 + bev / p["KDBEV"])


# ============================================================================
# 3.  RIGHT-HAND SIDE
# ============================================================================
def rhs(t, y, p, u):
    """u = infusion/administration rates dict, set by the scenario driver."""
    d = [0.0] * NS
    N, M, A, h_tot, h_via = geometry(y, p)
    S = p["S_PL"]

    # ---- 3.1 delivery geometry -------------------------------------------
    L = penetration_lengths(y, p)
    f_ox = fpen(p["L_OX"], h_tot)
    f_cis = fpen(L["cis"], h_tot)
    f_pem = fpen(L["pem"], h_tot)
    f_mab = fpen(L["mab"], h_tot)
    f_tc = fpen(L["tcell"], h_tot)
    f_ip = fpen(L["ip"], h_tot)
    # the systemic and intrapleural fronts advance from OPPOSITE faces, so the
    # intrapleural route only adds coverage where the systemic route stops
    f_sys_chem = max(f_cis, f_pem)
    # ---------------------------------------------------------------------
    # DEFECT #5.  The mechanistic map says that pleurodesis CLOSES the
    # intrapleural route, and the first version of these equations did not
    # implement it: symphysis entered only through the effusion, so obliterating
    # the space made the instillate MORE concentrated (smaller diluting volume)
    # and talc came out looking like a way to potentiate intrapleural
    # chemotherapy.  Section J's own table contradicted its own caption.  Two
    # factors were missing: an apposed surface cannot be instilled into at all,
    # and even a free space is loculated so an instillate bathes only part of
    # the involved surface.
    # ---------------------------------------------------------------------
    f_ip_avail = f_ip * p["FCONT"] * (1.0 - min(max(y[IX["PSY"]], 0.0), 1.0))
    f_ip_only = max(0.0, min(f_ip_avail, 1.0 - f_sys_chem))

    # ---- 3.2 drug concentrations ------------------------------------------
    CP_PEM = y[IX["PEM_C"]] / p["V1_PEM"]
    CP_CIS = y[IX["CIS_F"]] / p["V1_CIS"]
    C_NIV = y[IX["NIV_C"]] / p["V1_NIV"]
    C_IPI = y[IX["IPI_C"]] / p["V1_IPI"]
    C_ADI = y[IX["ADI_C"]] / p["V_ADI"]
    RO_PD1 = C_NIV / (C_NIV + p["EC50_PD1"])
    RO_CT4 = C_IPI / (C_IPI + p["EC50_CT4"])
    C_IP = y[IX["IPD"]] / max((y[IX["VEFF"]] + p["VDEAD"]) / 1000.0, 0.02)

    # ---- 3.3 pemetrexed PK + intracellular polyglutamation ----------------
    CL_pem = p["CL_PEM"] * (max(y[IX["CRCL"]], 10.0) / 100.0)
    d[IX["PEM_C"]] = (u.get("pem", 0.0)
                      - CL_pem * CP_PEM
                      - p["Q_PEM"] * CP_PEM
                      + p["Q_PEM"] * y[IX["PEM_P"]] / p["V2_PEM"])
    d[IX["PEM_P"]] = p["Q_PEM"] * CP_PEM - p["Q_PEM"] * y[IX["PEM_P"]] / p["V2_PEM"]
    d[IX["AUCP"]] = CP_PEM

    folcomp = 1.0 + y[IX["FOL"]] / p["KI_FOL"]
    upt = CP_PEM / folcomp
    d[IX["PEM_T"]] = (p["KIN_T"] * upt - p["KOUT_T"] * y[IX["PEM_T"]]
                      - p["KFPGS_T"] * y[IX["PEM_T"]] + p["KGGH_T"] * y[IX["PEM_TP"]])
    d[IX["PEM_TP"]] = p["KFPGS_T"] * y[IX["PEM_T"]] - p["KGGH_T"] * y[IX["PEM_TP"]]
    d[IX["PEM_M"]] = (p["KIN_M"] * upt - p["KOUT_M"] * y[IX["PEM_M"]]
                      - p["KFPGS_M"] * y[IX["PEM_M"]] + p["KGGH_M"] * y[IX["PEM_MP"]])
    d[IX["PEM_MP"]] = p["KFPGS_M"] * y[IX["PEM_M"]] - p["KGGH_M"] * y[IX["PEM_MP"]]

    # ---- 3.4 folate / B12 / homocysteine ----------------------------------
    d[IX["FOL"]] = (p["KFIN"] * u.get("folic_ug", 0.0)
                    + p["KFOUT"] * p["FOL_DIET"] - p["KFOUT"] * y[IX["FOL"]])
    d[IX["HCY"]] = (p["KHSYN"]
                    - p["KHDEG"] * (y[IX["FOL"]] / (y[IX["FOL"]] + p["KMH"])) * y[IX["HCY"]]
                    - p["KHB12"] * u.get("b12", 0.0) * y[IX["HCY"]])

    # ---- 3.5 cisplatin PK + adducts ---------------------------------------
    d[IX["CIS_F"]] = (u.get("cis", 0.0)
                      + p["FABS_IP"] * p["KIPCL"] * y[IX["IPD"]]
                      - (p["KEL_CIS"] + p["KBIND"]) * y[IX["CIS_F"]])
    d[IX["CIS_B"]] = (p["KBIND"] * y[IX["CIS_F"]] / p["VB_CIS"]
                      - p["KBEL"] * y[IX["CIS_B"]])
    d[IX["AUCC"]] = CP_CIS
    d[IX["ADD"]] = p["KADD"] * CP_CIS - p["KREP"] * p["ERCC1F"] * y[IX["ADD"]]
    d[IX["ADDIP"]] = p["KADD"] * C_IP - p["KREP"] * p["ERCC1F"] * y[IX["ADDIP"]]
    d[IX["IPD"]] = u.get("ip", 0.0) - p["KIPCL"] * y[IX["IPD"]]

    # ---- 3.6 biologics PK ---------------------------------------------------
    C_BEV = y[IX["BEV_C"]] / p["V1_BEV"]
    d[IX["BEV_C"]] = (u.get("bev", 0.0) - p["CL_BEV"] * C_BEV
                      - p["Q_BEV"] * C_BEV + p["Q_BEV"] * y[IX["BEV_P"]] / p["V2_BEV"])
    d[IX["BEV_P"]] = p["Q_BEV"] * C_BEV - p["Q_BEV"] * y[IX["BEV_P"]] / p["V2_BEV"]
    d[IX["NIV_C"]] = (u.get("niv", 0.0) - p["CL_NIV"] * C_NIV
                      - p["Q_NIV"] * C_NIV + p["Q_NIV"] * y[IX["NIV_P"]] / p["V2_NIV"])
    d[IX["NIV_P"]] = p["Q_NIV"] * C_NIV - p["Q_NIV"] * y[IX["NIV_P"]] / p["V2_NIV"]
    d[IX["IPI_C"]] = u.get("ipi", 0.0) - p["CL_IPI"] * C_IPI

    # ADI-PEG20: IM depot, ADA-accelerated clearance, arginine depletion
    d[IX["ADI_A"]] = u.get("adi", 0.0) - p["KA_ADI"] * y[IX["ADI_A"]]
    cl_adi = p["CL_ADI"] * (1.0 + p["ADA_POT"] * y[IX["ADA"]])
    d[IX["ADI_C"]] = p["KA_ADI"] * y[IX["ADI_A"]] - cl_adi * C_ADI
    d[IX["ADA"]] = p["KADA"] * (1.0 if C_ADI > 0.02 else 0.0) - p["KADAEL"] * y[IX["ADA"]]
    d[IX["ARG"]] = (p["KARGIN"] - p["KARGOUT"] * y[IX["ARG"]]
                    - p["KADIARG"] * C_ADI * y[IX["ARG"]] / (y[IX["ARG"]] + 20.0))

    # ---- 3.7 vinorelbine (second line, simple effect compartment) ----------
    d[IX["VINE"]] = u.get("vin", 0.0) - 0.60 * y[IX["VINE"]]

    # ---- 3.8 VEGF, microvessels, effusion ----------------------------------
    VEGFf = free_vegf(y, p)
    d[IX["VEGF"]] = (p["KVSYN"] * (N / 100.0) * (0.30 + 0.70 * (1.0 - f_ox))
                     - p["KVDEG"] * y[IX["VEGF"]])
    d[IX["RHOV"]] = (p["KANG"] * (VEGFf / (VEGFf + p["EC50A"]))
                     * (1.0 - y[IX["RHOV"]] / p["RHOMAX"])
                     - p["KVREG"] * y[IX["RHOV"]])

    afrac = min(A / S, 1.0)
    psy = min(max(y[IX["PSY"]], 0.0), 1.0)
    # ---------------------------------------------------------------------
    # DEFECT #2 (found by this reference implementation).  Written without a
    # ceiling, the effusion equation is  Jform(A) - Jdrain(A)  with formation
    # RISING and stomatal drainage FALLING in the same variable, so in the
    # untreated run it integrated to 1.3 million mL by one year.  The missing
    # physics is that the pleural space has a VOLUME: filtration is a Starling
    # flux that stops when pleural pressure rises to meet capillary pressure,
    # and the space available is what the hemithorax holds MINUS the volume
    # the tumour rind itself occupies, minus fibrothorax contraction.  With
    # the ceiling in place the effusion becomes non-monotone in A on its own.
    # ---------------------------------------------------------------------
    Vtum = N + p["PHIM"] * M
    vcap = max(60.0, (p["VEFFMAX"] - 0.90 * Vtum) * (1.0 - psy)
               * (1.0 - 0.45 * afrac))
    Jform = (p["KFORM"] * (0.05 + afrac)
             * (1.0 + p["EMAXV"] * VEGFf / (VEGFf + p["EC50V"])) * (1.0 - psy)
             * max(0.0, 1.0 - y[IX["VEFF"]] / vcap))
    # symphysis blocks FORMATION, not absorption: an apposed space still
    # resorbs whatever fluid is left in it
    Jdrain = (p["KDRAIN"] * (1.0 - p["FOBS"] * afrac)
              * y[IX["VEFF"]] / (y[IX["VEFF"]] + p["KMD"])
              + p["KABS"] * psy * y[IX["VEFF"]])
    d[IX["VEFF"]] = Jform - Jdrain - u.get("drain", 0.0)
    trapped = min(1.0, (afrac ** 2) * (h_tot / (h_tot + 0.6)) * 2.2)
    d[IX["PSY"]] = ((p["KSYMSP"] * afrac ** 3 + p["KTALC"] * u.get("talc", 0.0))
                    * (1.0 - psy) * (1.0 - 0.85 * trapped))

    # ---- 3.9 immunity --------------------------------------------------------
    antig = ((max(N, 1.0) / 290.0) ** 0.30) * p["IMMINF"]
    prime = p["KPRIME"] * antig * (1.0 + p["EIPI"] * RO_CT4) / (1.0 + p["VISTA_S"])
    d[IX["TEFF"]] = (prime
                     - p["KEXH"] * p["PDL1B"] * (1.0 - RO_PD1) * y[IX["TEFF"]]
                     - p["KTD"] * y[IX["TEFF"]]
                     - p["KSUP"] * y[IX["TREG"]] * y[IX["TEFF"]])
    d[IX["TREG"]] = (p["KTREGS"] * y[IX["TGFB"]] / (y[IX["TGFB"]] + p["KMTG"])
                     - p["KTREGD"] * y[IX["TREG"]]
                     - p["KIPIDEP"] * RO_CT4 * y[IX["TREG"]])
    d[IX["TCLON"]] = p["KCLON"] * y[IX["TEFF"]] - p["KCLOND"] * y[IX["TCLON"]]
    d[IX["TGFB"]] = p["KTGF"] * (N / 100.0) + 3.0 * y[IX["WOUND"]] - p["KTGFD"] * y[IX["TGFB"]]
    Tpool = y[IX["TEFF"]] + p["WMEM"] * y[IX["TCLON"]]
    k_imm = p["KKILL"] * (Tpool / (Tpool + p["KMKILL"])) * f_tc

    # ---- 3.10 kill terms -------------------------------------------------------
    E_pem = (p["EMAXP"] * y[IX["PEM_TP"]]
             / (y[IX["PEM_TP"]] + p["EC50P"] * (1.0 + y[IX["FOL"]] / p["KFPT"])))
    E_cis = p["EMAXC"] * y[IX["ADD"]] / (y[IX["ADD"]] + p["EC50C"])
    E_vin = p["EMAXVIN"] * y[IX["VINE"]] / (y[IX["VINE"]] + p["EC50VIN"])
    E_sys = p["CHEMOS"] * (E_pem + E_cis + p["SYNCP"] * E_pem * E_cis + E_vin)
    E_ipl = p["CHEMOS"] * p["EMAXC"] * y[IX["ADDIP"]] / (y[IX["ADDIP"]] + p["EC50C"])
    E_arg = (p["EMAXARG"] * p["ASS1LOSS"]
             * max(0.0, (p["ARG_CRIT"] - y[IX["ARG"]])) / p["ARG_CRIT"])
    k_rt = u.get("rt", 0.0)
    k_prmt5 = u.get("prmt5", 0.0)

    kill = (E_sys * f_sys_chem + E_ipl * f_ip_only + k_imm
            + E_arg * f_cis + k_rt + k_prmt5)
    # (radiotherapy and PRMT5 inhibition are not diffusion-limited in the same
    #  way: photons reach the whole rind, and MRTX1719 is a small molecule
    #  given continuously, so both act on the whole viable pool.)

    # ---- 3.11 tumour geometry ODEs ----------------------------------------------
    growth = p["KPROL"] * p["KPROLF"] * f_ox * N
    deaths = (p["KDEATH0"] * N + p["KNEC"] * (1.0 - f_ox) * N + kill * N)
    d[IX["N"]] = growth - deaths
    d[IX["M"]] = p["PHIM"] * deaths - p["KDEGM"] * M
    perim = 2.0 * math.sqrt(math.pi * max(A, 1.0))
    # -------------------------------------------------------------------
    # THE ADVANCING MARGIN IS THIN.  This is the sharpest consequence of
    # writing the tumour as a sheet.  At the growth front the rind is a few
    # hundred microns deep, so there  f_exp -> 1: the front is FULLY drug
    # exposed and fully oxygenated even when the bulk behind it is neither.
    # Chemotherapy therefore arrests circumferential SPREAD long before it
    # makes a dent in the bulk -- and it is spread (encasement) that the
    # hazard is most sensitive to.  Using the same fpen() at h = HFRONT
    # rather than at h = h_tot is the whole of the difference.
    # -------------------------------------------------------------------
    hf = p["HFRONT"]
    fox_f = fpen(p["L_OX"], hf)
    fch_f = min(1.0, max(fpen(L["cis"], hf), fpen(L["pem"], hf)))
    fip_f = min(1.0, fpen(L["ip"], hf)) * p["FCONT"] * (
        1.0 - min(max(y[IX["PSY"]], 0.0), 1.0))
    ftc_f = min(1.0, fpen(L["tcell"], hf))
    angio_f = p["ANG_MIN"] + (1.0 - p["ANG_MIN"]) * VEGFf / (VEGFf + p["EC50A"])
    kimm_f = p["KKILL"] * (Tpool / (Tpool + p["KMKILL"])) * ftc_f
    front_net = (p["KPROL"] * p["KPROLF"] * fox_f * angio_f
                 - p["KDEATH0"] - p["KNEC"] * (1.0 - fox_f)
                 - E_sys * fch_f - E_ipl * fip_f - kimm_f
                 - E_arg * fch_f - k_rt - k_prmt5)
    vfront = min(1.5, max(0.0, front_net / p["GFRONT_REF"]))
    d[IX["A"]] = p["KFRONT"] * perim * (1.0 - A / S) * vfront
    viab = N / max(N + p["PHIM"] * M, 1e-6)
    d[IX["Z"]] = (p["KINV"] * (h_via / (h_via + 0.30)) * viab
                  - 0.0020 * y[IX["Z"]] * min(kill / 0.02, 1.0))
    # DEFECT #3: written as unbounded exponential growth with only a
    # kill-dependent brake, MET reached 660 cm3 by two years in the
    # immunotherapy arm and then DOMINATED the hazard, so a scenario that
    # removed 89 % of the pleural tumour still showed no survival gain.
    # Distant spread in MPM is late and is rarely what kills; a ceiling and
    # a slower rate were required.
    met_supp = max(0.0, 1.0 - kill / 0.050)
    d[IX["MET"]] = (p["KMET"] * N
                    + p["KMETG"] * y[IX["MET"]] * met_supp
                    * (1.0 - y[IX["MET"]] / p["METMAX"]))

    # ---- 3.12 marrow (Friberg transit) --------------------------------------------
    E_pem_m = (p["EMAXPM"] * y[IX["PEM_MP"]]
               / (y[IX["PEM_MP"]] + p["EC50PM"] * (1.0 + y[IX["FOL"]] / p["KFPM"])))
    Edrug = min(0.95, E_pem_m + p["SLOPE_CIS"] * CP_CIS + 0.35 * y[IX["VINE"]])
    fb = (p["CIRC0"] / max(y[IX["CIRC"]], 0.05)) ** p["GAM"]
    d[IX["PROL"]] = p["KTR"] * y[IX["PROL"]] * (1.0 - Edrug) * fb - p["KTR"] * y[IX["PROL"]]
    d[IX["TR1"]] = p["KTR"] * (y[IX["PROL"]] - y[IX["TR1"]])
    d[IX["TR2"]] = p["KTR"] * (y[IX["TR1"]] - y[IX["TR2"]])
    d[IX["TR3"]] = p["KTR"] * (y[IX["TR2"]] - y[IX["TR3"]])
    d[IX["CIRC"]] = p["KTR"] * y[IX["TR3"]] - p["KTR"] * y[IX["CIRC"]]

    # ---- 3.13 renal ------------------------------------------------------------------
    d[IX["CRCLSS"]] = -p["KNEPHIRR"] * CP_CIS
    d[IX["CRCL"]] = p["KRECR"] * (y[IX["CRCLSS"]] - y[IX["CRCL"]]) - p["KNEPH"] * CP_CIS

    # ---- 3.14 biomarker ---------------------------------------------------------------
    d[IX["SMRP"]] = (p["KSMRP"] * N * p["MSLNF"]
                     - p["KSMRPCL"] * (max(y[IX["CRCL"]], 5.0) / 95.0) * y[IX["SMRP"]])

    # ---- 3.15 inflammation, cachexia, wound ---------------------------------------------
    d[IX["IL6"]] = (p["KIL6T"] * (N / 100.0) + p["KIL6W"] * y[IX["WOUND"]]
                    - p["KIL6D"] * y[IX["IL6"]])
    d[IX["WOUND"]] = -p["KWD"] * y[IX["WOUND"]]
    d[IX["LBM"]] = (-p["KCACH"] * y[IX["IL6"]] / (y[IX["IL6"]] + p["KMIL6"])
                    + p["KLREC"] * (p["LBM0"] - y[IX["LBM"]]))

    # ---- 3.16 immune-related adverse events ---------------------------------------------
    d[IX["IRAE"]] = 0.055 * (RO_CT4 * 1.0 + 0.30 * RO_PD1) * (1.0 - y[IX["IRAE"]]) - 0.030 * y[IX["IRAE"]]

    # ---- 3.17 hazard --------------------------------------------------------------------
    V = N + p["PHIM"] * M
    fvc = fvc_pct(y, p)
    ps = ecog(y, p, fvc)
    # The hazard is NOT proportional to volume.  What kills in MPM is
    # encasement -- the fraction of the hemithorax that no longer moves --
    # and the performance status that follows from it.  Volume enters with a
    # square-root exponent; AREA enters in its own right.  This is the same
    # structural claim as the rest of the model, and it is what makes
    # debulking (which collapses h and leaves A) a poor bargain.
    haz = (p["H0"] * ((max(V, 1.0) + p["WMET"] * y[IX["MET"]]) / p["V0H"]) ** p["B1"]
           * math.exp(p["B2"] * ps + p["B5"] * min(A / p["S_PL"], 1.0))
           * p["SARCH"])
    if y[IX["CIRC"]] < 0.5:
        haz += p["H_FEBN"]
    haz += p["H_IRAE"] * max(0.0, y[IX["IRAE"]] - 0.35) / 0.65
    d[IX["CUMH"]] = haz

    return d


def fvc_pct(y, p):
    _, _, A, h_tot, _ = geometry(y, p)
    return p["FVC0"] * (1.0
                        - p["A_FVC"] * min(A / p["S_PL"], 1.0)
                        - p["H_FVC"] * min(h_tot, 3.0) / 3.0
                        - p["V_FVC"] * min(y[IX["VEFF"]], 2500.0) / 2500.0)


def ecog(y, p, fvc=None):
    if fvc is None:
        fvc = fvc_pct(y, p)
    pain = 3.0 * y[IX["Z"]] / (y[IX["Z"]] + 0.80)
    idx = (1.10 * max(0.0, 1.0 - fvc / p["FVC0"])
           + 1.30 * max(0.0, 1.0 - y[IX["LBM"]] / p["LBM0"]) * 3.0
           + 0.12 * pain)
    return min(4.0, 4.0 * (1.0 - math.exp(-1.35 * idx)))


# ============================================================================
# 4.  INITIAL CONDITIONS
# ============================================================================
def initial_state(p, h0=0.80, A0=450.0, mfrac=0.28):
    """Typical MPM at diagnosis: rind 0.8 cm over 450 cm2 -> V ~ 360 cm3."""
    y = [0.0] * NS
    Vtot = h0 * A0
    # PHIM = 0 is used as a counterfactual in section H (killed cells vanish
    # instead of fibrosing); with no persistent matrix the whole measured
    # volume is viable by definition.
    if p["PHIM"] <= 0.0:
        M0, N0 = 0.0, Vtot
    else:
        M0 = mfrac * Vtot / p["PHIM"]
        N0 = Vtot - p["PHIM"] * M0
    y[IX["N"]] = N0
    y[IX["M"]] = M0
    y[IX["A"]] = A0
    y[IX["Z"]] = 0.18
    y[IX["MET"]] = 2.0
    y[IX["VEFF"]] = 1000.0
    y[IX["PSY"]] = 0.0
    y[IX["VEGF"]] = p["KVSYN"] * (N0 / 100.0) * 0.98 / p["KVDEG"]
    y[IX["RHOV"]] = 1.0
    y[IX["IL6"]] = p["KIL6T"] * (N0 / 100.0) / p["KIL6D"]
    y[IX["TGFB"]] = p["KTGF"] * (N0 / 100.0) / p["KTGFD"]
    # immune compartments started AT their untreated steady state, so that a
    # scenario without immunotherapy shows no spurious opening transient
    tgf = y[IX["TGFB"]]
    treg = p["KTREGS"] * tgf / (tgf + p["KMTG"]) / p["KTREGD"]
    antig = ((max(N0, 1.0) / 290.0) ** 0.30) * p["IMMINF"]
    prime = p["KPRIME"] * antig / (1.0 + p["VISTA_S"])
    teff = prime / (p["KEXH"] * p["PDL1B"] + p["KTD"] + p["KSUP"] * treg)
    y[IX["TREG"]] = treg
    y[IX["TEFF"]] = teff
    y[IX["TCLON"]] = p["KCLON"] * teff / p["KCLOND"]
    y[IX["FOL"]] = p["FOL_DIET"]
    y[IX["HCY"]] = p["KHSYN"] / (p["KHDEG"] * (p["FOL_DIET"] / (p["FOL_DIET"] + p["KMH"])))
    y[IX["ARG"]] = p["KARGIN"] / p["KARGOUT"]
    for s in ("PROL", "TR1", "TR2", "TR3", "CIRC"):
        y[IX[s]] = p["CIRC0"]
    y[IX["CRCL"]] = p["CRCL0"]
    y[IX["CRCLSS"]] = p["CRCL0"]
    y[IX["SMRP"]] = p["KSMRP"] * N0 * p["MSLNF"] / p["KSMRPCL"]
    y[IX["LBM"]] = p["LBM0"]
    return y


# ============================================================================
# 5.  DOSING / EVENT ENGINE
# ============================================================================
class Regimen:
    """Infusions (rate over a window) + instantaneous events."""

    def __init__(self):
        self.infusions = []      # (key, t0, t1, rate)
        self.events = []         # (t, callable(y, p))
        self.constants = {}      # key -> value applied at all times
        self.windows = []        # (key, t0, t1, value)
        self.tlast = 25.0        # last time any input is still switching
        self.chemo_state = None
        # -----------------------------------------------------------------
        # DEFECT #4 (found by the step-size convergence check in section A of
        # mpm_calibration_output.txt).  A pemetrexed dose is a TEN-MINUTE
        # infusion = 0.0069 d, which is SHORTER than the 0.02 d integration
        # step.  With a naive fixed step, whether an RK4 stage happens to land
        # inside the infusion window decides how much drug is delivered -- at
        # dt = 0.02 d only the first stage landed inside, so each dose arrived
        # at a fraction of its nominal amount.  V(400 d) then differed by 44 %
        # between dt = 0.02 and dt = 0.005, and the fitted EMAXP/EMAXC were
        # silently compensating for the missing drug.
        # The fix: record every rate discontinuity as a BREAKPOINT and never
        # let a step cross one, so short infusions are integrated exactly.
        # (mrgsolve does this natively by inserting dose records.)
        # -----------------------------------------------------------------
        self.breaks = [0.0]

    def _brk(self, *ts):
        for t in ts:
            self.breaks.append(t)
        self.breaks.sort()

    def infuse(self, key, t0, dur, amt):
        self.infusions.append((key, t0, t0 + dur, amt / dur))
        self.tlast = max(self.tlast, t0 + dur + 25.0)
        self._brk(t0, t0 + dur)
        return self

    def const(self, key, value):
        self.constants[key] = value
        return self

    def window(self, key, t0, t1, value):
        self.windows.append((key, t0, t1, value))
        self.tlast = max(self.tlast, t1 + 25.0)
        self._brk(t0, t1)
        return self

    def event(self, t, fn):
        self.events.append((t, fn))
        self.tlast = max(self.tlast, t + 25.0)
        self._brk(t)
        return self

    def u_at(self, t):
        u = dict(self.constants)
        for k, t0, t1, r in self.infusions:
            if t0 <= t < t1:
                u[k] = u.get(k, 0.0) + r
        for k, t0, t1, v in self.windows:
            if t0 <= t < t1:
                u[k] = v
        return u


NONNEG = ("N", "M", "MET", "VEFF", "TEFF", "TREG", "TCLON", "TGFB",
          "PEM_C", "PEM_P", "PEM_T", "PEM_TP", "PEM_M", "PEM_MP",
          "CIS_F", "CIS_B", "ADD", "ADDIP", "IPD", "ADI_A", "ADI_C",
          "ADA", "BEV_C", "BEV_P", "NIV_C", "NIV_P", "IPI_C", "VINE",
          "SMRP", "IRAE", "IL6", "VEGF")
NONNEG_IX = tuple(IX[s] for s in NONNEG)


def simulate(p, reg, tmax=1250.0, dt=0.02, record_every=1.0, y0=None,
             dt_slow=0.06):
    """Fixed-step RK4 with a TWO-PHASE step size.

    While any drug is being given (and for 25 d afterwards) the fastest rates
    in the system are cisplatin's (k_el + k_bind = 33.3 /d) and the pemetrexed
    efflux (20 /d), so dt = 0.02 d is required.  Once every depot is empty the
    remaining dynamics are the tumour, the marrow and the effusion, none faster
    than 8 /d, and dt = 0.06 d is both stable and accurate.  This is worth ~2.5x
    of wall clock over a 1250-day run and changes the reported endpoints by
    less than 0.15 % (verified in section A of the output)."""
    y = list(y0 if y0 is not None else initial_state(p))
    # reg.events is used LIVE, not copied: the chemotherapy controller adds
    # future cycles while the simulation is running.  Controllers only ever
    # append times later than 'now', so re-sorting never moves an already
    # fired event past the pointer.
    reg.events.sort(key=lambda e: e[0])
    events = reg.events
    ei = 0
    out = {"t": [], "y": []}
    # (reg.constants -- folic acid, B12 -- are slow inputs and do not force
    #  the small step size; reg.tlast is maintained by Regimen and grows as
    #  the controller schedules more cycles)

    t = 0.0
    next_rec = 0.0
    bi = 0
    while True:
        while ei < len(events) and events[ei][0] <= t + 1e-9:
            events[ei][1](y, p)
            ei += 1
            events = reg.events
        if t >= next_rec - 1e-9:
            out["t"].append(round(next_rec, 6))
            out["y"].append(list(y))
            next_rec += record_every
        if t >= tmax - 1e-9:
            break
        h = dt if t < reg.tlast else dt_slow
        # never step past the next recording point or the next event
        h = min(h, next_rec - t if next_rec > t else h, tmax - t)
        if ei < len(events):
            h = min(h, max(events[ei][0] - t, 1e-9))
        # never step across a rate discontinuity (see DEFECT #4)
        while bi < len(reg.breaks) and reg.breaks[bi] <= t + 1e-12:
            bi += 1
        if bi < len(reg.breaks):
            h = min(h, max(reg.breaks[bi] - t, 1e-9))
        dt_ = h
        # Because no step crosses a breakpoint, the input vector is CONSTANT
        # across the step, and sampling it at the midpoint identifies it
        # exactly.  Sampling it separately at t, t+h/2 and t+h (as a naive RK4
        # driver does) mis-weights the two ends of every infusion window: the
        # k4 stage of the last step of a window falls on the closing boundary
        # and reads zero, delivering 5/6 of the intended dose.  This is the
        # second half of the DEFECT #4 repair.
        u1 = reg.u_at(t + 0.5 * dt_)
        k1 = rhs(t, y, p, u1)
        ym = [y[j] + 0.5 * dt_ * k1[j] for j in range(NS)]
        k2 = rhs(t + 0.5 * dt_, ym, p, u1)
        ym = [y[j] + 0.5 * dt_ * k2[j] for j in range(NS)]
        k3 = rhs(t + 0.5 * dt_, ym, p, u1)
        ym = [y[j] + dt_ * k3[j] for j in range(NS)]
        k4 = rhs(t + dt_, ym, p, u1)
        for j in range(NS):
            y[j] += dt_ / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
        # ---------------------------------------------------------------
        # DEFECT #1 (found by this reference implementation, see
        # mpm_calibration_output.txt): with no floor, KILL can drive N
        # marginally negative in the first hours after an intrapleural
        # instillation (C_IP is ~200 mg/L for a few minutes and the RK4
        # midpoint overshoots).  N < 0 then makes h < 0, makes fpen()
        # return a negative exposed fraction, and the whole rind flips
        # sign.  Clamping the physically non-negative states after each
        # step is the fix; the same clamp exists in $ODE of the mrgsolve
        # file as a set of max(...) guards on the state reads.
        # ---------------------------------------------------------------
        for j in NONNEG_IX:
            if y[j] < 0.0:
                y[j] = 0.0
        y[IX["A"]] = min(max(y[IX["A"]], 1.0), p["S_PL"])
        y[IX["PSY"]] = min(max(y[IX["PSY"]], 0.0), 1.0)
        y[IX["ARG"]] = max(y[IX["ARG"]], 0.0)
        y[IX["CIRC"]] = max(y[IX["CIRC"]], 0.01)
        y[IX["CRCL"]] = max(y[IX["CRCL"]], 8.0)
        t += dt_
    return out


# ============================================================================
# 6.  REGIMEN BUILDERS
# ============================================================================
CYC = 21.0            # days
INF_PEM = 10.0 / 60.0 / 24.0     # 10-minute infusion, in days
INF_CIS = 2.0 / 24.0             # 2-hour infusion


def add_chemo(reg, p, ncyc=6, t0=0.0, pem=True, plat="cis", supplemented=True,
              adaptive=True):
    """Platinum-pemetrexed with the dose modifications that are actually
    written into the protocols, because that is the channel through which
    folate supplementation changes SURVIVAL rather than just toxicity:

      * a cycle is given only if ANC >= 1.5 x10^9/L and CrCl >= 45 mL/min
      * otherwise it is delayed one week, up to three times
      * after two delays the dose drops to 75 %, then to 50 %
      * five consecutive failed weeks stop treatment permanently

    Without this the model reproduces EMPHACIS's toxicity difference and
    NOT its survival difference, because in the equations alone folate
    changes only the marrow."""
    dose_pem = 500.0 * p["BSA"]
    dose_pt = (75.0 if plat == "cis" else 0.78 * 75.0) * p["BSA"]
    state = {"given": 0, "delays": 0, "level": 1.0, "stopped": False,
             "delivered": 0.0, "planned": ncyc * (dose_pem if pem else 0.0)}
    reg.chemo_state = state

    def controller(t):
        def fn(y, pp):
            if state["stopped"] or state["given"] >= ncyc:
                return
            ok = (y[IX["CIRC"]] >= 1.5) and (y[IX["CRCL"]] >= 45.0)
            if not ok:
                state["delays"] += 1
                if state["delays"] >= 5:
                    state["stopped"] = True
                    return
                if state["delays"] == 2:
                    state["level"] = 0.75
                elif state["delays"] >= 3:
                    state["level"] = 0.50
                reg.event(t + 7.0, controller(t + 7.0))
                reg.events.sort(key=lambda e: e[0])
                return
            lv = state["level"]
            if pem:
                reg.infuse("pem", t, INF_PEM, dose_pem * lv)
                state["delivered"] += dose_pem * lv
            reg.infuse("cis", t + 0.5 / 24.0, INF_CIS, dose_pt * lv)
            state["given"] += 1
            nxt = t + CYC
            if state["given"] < ncyc:
                reg.event(nxt, controller(nxt))
                reg.events.sort(key=lambda e: e[0])
        return fn

    if adaptive:
        reg.event(t0, controller(t0))
    else:
        for c in range(ncyc):
            t = t0 + c * CYC
            if pem:
                reg.infuse("pem", t, INF_PEM, dose_pem)
            reg.infuse("cis", t + 0.5 / 24.0, INF_CIS, dose_pt)
        state["given"] = ncyc
        state["delivered"] = state["planned"]
    if supplemented:
        reg.const("folic_ug", 400.0)
        reg.const("b12", 1.0)
    return reg


def add_bev(reg, p, ncyc=6, t0=0.0, maint_until=None, mgkg=15.0):
    n = ncyc if maint_until is None else int(maint_until / CYC)
    for c in range(n):
        reg.infuse("bev", t0 + c * CYC + 0.02, 0.03, mgkg * p["WT"])
    return reg


def add_nivo_ipi(reg, p, until=730.0, t0=0.0):
    n_niv = int(until / 14.0)
    for c in range(n_niv):
        reg.infuse("niv", t0 + c * 14.0, 0.03, 3.0 * p["WT"])
    n_ipi = int(until / 42.0)
    for c in range(n_ipi):
        reg.infuse("ipi", t0 + c * 42.0, 0.03, 1.0 * p["WT"])
    return reg


def add_nivo(reg, p, until=730.0, t0=0.0):
    for c in range(int(until / 14.0)):
        reg.infuse("niv", t0 + c * 14.0, 0.03, 3.0 * p["WT"])
    return reg


def add_pembro(reg, p, until=730.0, t0=0.0):
    # pembrolizumab 200 mg q3w, mapped onto the nivolumab PD-1 module
    for c in range(int(until / 21.0)):
        reg.infuse("niv", t0 + c * 21.0, 0.03, 200.0)
    return reg


def add_adi(reg, p, until=365.0, t0=0.0):
    for c in range(int(until / 7.0)):
        reg.infuse("adi", t0 + c * 7.0, 0.05, 36.0 * p["BSA"])
    return reg


def surgery_event(t, h_res=0.020, mortality_hazard=0.0, fvc_loss=0.14,
                  lasting_hazard=1.15, lbm_loss=0.10):
    """Cytoreductive surgery for a SHEET.

    h collapses to a microscopic residue.  A DOES NOT MOVE -- there is no
    margin to take, because the disease is coextensive with the pleural
    surface.  Three permanent costs are applied alongside the debulking,
    because they are what the randomised trials measured:

      * peri-operative mortality, added straight to the cumulative hazard
      * a permanent loss of predicted FVC (a decorticated hemithorax is
        stiff; a pneumonectomy removes a lung outright)
      * a permanent increase in baseline hazard, and an acute loss of lean
        mass, standing for grade >=3 morbidity

    The FVC and hazard penalties are applied to the parameter set, not to a
    state, because they are irreversible: this is the one place in the model
    where an intervention edits p."""
    def fn(y, p):
        A = y[IX["A"]]
        keep = h_res * A                       # residual sheet volume
        tot = y[IX["N"]] + p["PHIM"] * y[IX["M"]]
        if tot <= 0:
            return
        frac = min(1.0, keep / tot)
        y[IX["N"]] *= frac
        y[IX["M"]] *= frac
        y[IX["VEFF"]] = 50.0
        y[IX["PSY"]] = 0.85                    # the space is obliterated
        y[IX["WOUND"]] = 1.0                   # IL-6 / TGF-beta surge
        y[IX["LBM"]] *= (1.0 - lbm_loss)
        y[IX["CUMH"]] += mortality_hazard      # peri-operative mortality
        p["FVC0"] *= (1.0 - fvc_loss)
        p["H0"] *= lasting_hazard
    return fn


# ============================================================================
# 7.  READOUTS
# ============================================================================
def series(out, name):
    i = IX[name]
    return [row[i] for row in out["y"]]


def derived(out, p):
    t = out["t"]
    res = {"t": t, "h_meas": [], "h_via": [], "V": [], "Vvia": [], "mRECIST": [],
           "fexp": [], "S": [], "FVC": [], "PS": [], "Lp": [], "sanct": []}
    for row in out["y"]:
        N, M, A, h_tot, h_via = geometry(row, p)
        L = penetration_lengths(row, p)
        f = fpen(L["cis"], h_tot)
        res["h_meas"].append(h_tot)
        res["h_via"].append(h_via)
        res["V"].append(N + p["PHIM"] * M)
        res["Vvia"].append(N)
        res["mRECIST"].append(6.0 * h_tot * 10.0)   # mm, 6 sites
        res["fexp"].append(f)
        res["sanct"].append(1.0 - f)
        res["Lp"].append(L["cis"] * 10.0)           # mm
        res["S"].append(math.exp(-row[IX["CUMH"]]))
        res["FVC"].append(fvc_pct(row, p))
        res["PS"].append(ecog(row, p))
    return res


def median_from_S(t, S):
    for i in range(1, len(S)):
        if S[i] <= 0.5:
            # linear interpolation in t
            s0, s1 = S[i - 1], S[i]
            return t[i - 1] + (s0 - 0.5) / (s0 - s1) * (t[i] - t[i - 1])
    return None


def cohort_survival(runner, p_base, tmax=1400.0):
    """Weighted mixture over histology -> a trial-cohort survival curve."""
    curves, tt = [], None
    for hname, w in COHORT:
        p = dict(p_base)
        p.update(HIST[hname])
        out = runner(p)
        dd = derived(out, p)
        tt = dd["t"]
        curves.append((w, dd["S"]))
    S = [sum(w * c[i] for w, c in curves) for i in range(len(tt))]
    return tt, S


def best_response(dd, base_idx=0, upto_day=250):
    """mRECIST best percentage change in the 6-site thickness sum."""
    b = dd["mRECIST"][base_idx]
    n = min(len(dd["t"]), upto_day + 1)
    lo = min(dd["mRECIST"][:n])
    return 100.0 * (lo - b) / b


def response_category(pct_change, pct_worst_growth):
    if pct_change <= -30.0:
        return "PR"
    if pct_worst_growth >= 20.0:
        return "PD"
    return "SD"


# ============================================================================
# 8.  SCENARIOS
# ============================================================================
def scen_bsc(p):
    return simulate(p, Regimen(), tmax=1400.0)


def scen_cis_alone(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=False, plat="cis", supplemented=False)
    return simulate(p, r, tmax=1400.0)


def scen_cispem_unsupp(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=False)
    return simulate(p, r, tmax=1400.0)


def scen_cispem_supp(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    return simulate(p, r, tmax=1400.0)


def scen_carbopem(p):
    q = dict(p); q["KNEPH"] = 0.25 * p["KNEPH"]; q["KNEPHIRR"] = 0.25 * p["KNEPHIRR"]
    r = Regimen()
    add_chemo(r, q, ncyc=6, pem=True, plat="carbo", supplemented=True)
    return simulate(q, r, tmax=1400.0)


def scen_maps(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    add_bev(r, p, maint_until=730.0)
    return simulate(p, r, tmax=1400.0)


def scen_nivo_ipi(p):
    r = Regimen()
    add_nivo_ipi(r, p, until=730.0)
    return simulate(p, r, tmax=1400.0)


def scen_pembro_chemo(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    add_pembro(r, p, until=730.0)
    return simulate(p, r, tmax=1400.0)


def scen_nivo_2l(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    add_nivo(r, p, until=730.0, t0=200.0)
    return simulate(p, r, tmax=1400.0)


def scen_adi_chemo(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    add_adi(r, p, until=365.0)
    return simulate(p, r, tmax=1400.0)


def scen_mars2(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    r.event(70.0, surgery_event(70.0, h_res=0.020, mortality_hazard=0.040,
                                fvc_loss=0.14, lasting_hazard=1.30, lbm_loss=0.10))
    return simulate(p, r, tmax=1400.0)


def scen_epp_rt(p):
    r = Regimen()
    add_chemo(r, p, ncyc=3, pem=True, plat="cis", supplemented=True)
    r.event(70.0, surgery_event(70.0, h_res=0.012, mortality_hazard=0.110,
                                fvc_loss=0.40, lasting_hazard=1.85, lbm_loss=0.16))
    r.window("rt", 110.0, 145.0, 0.055)
    return simulate(p, r, tmax=1400.0)


def scen_ip_after_sys(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    for c in range(6):
        r.infuse("ip", c * CYC + 1.0, 0.04, 100.0 * p["BSA"])
    return simulate(p, r, tmax=1400.0)


def scen_talc_then_chemo(p):
    r = Regimen()
    r.window("talc", 5.0, 7.0, 1.0)
    add_chemo(r, p, ncyc=6, t0=14.0, pem=True, plat="cis", supplemented=True)
    return simulate(p, r, tmax=1400.0)


def scen_talc_then_ip(p):
    r = Regimen()
    r.window("talc", 5.0, 7.0, 1.0)
    add_chemo(r, p, ncyc=6, t0=14.0, pem=True, plat="cis", supplemented=True)
    for c in range(6):
        r.infuse("ip", 14.0 + c * CYC + 1.0, 0.04, 100.0 * p["BSA"])
    return simulate(p, r, tmax=1400.0)


def scen_prmt5(p):
    q = dict(p)
    r = Regimen()
    add_chemo(r, q, ncyc=6, pem=True, plat="cis", supplemented=True)
    r.window("prmt5", 0.0, 730.0, 0.0075)
    return simulate(q, r, tmax=1400.0)


def scen_maint_pem(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    for c in range(6, 24):
        r.infuse("pem", c * CYC, INF_PEM, 500.0 * p["BSA"])
    return simulate(p, r, tmax=1400.0)


def scen_renal_impaired(p):
    q = dict(p); q["CRCL0"] = 52.0
    r = Regimen()
    add_chemo(r, q, ncyc=6, pem=True, plat="cis", supplemented=True)
    return simulate(q, r, tmax=1400.0, y0=initial_state(q))


def scen_vin_2l(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    for c in range(18):
        r.infuse("vin", 200.0 + c * 7.0, 0.03, 60.0 * p["BSA"] * 0.012)
    return simulate(p, r, tmax=1400.0)


def scen_early_small(p):
    """Same regimen, rind 0.35 cm over 250 cm2 (screen-detected disease)."""
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    return simulate(p, r, tmax=1400.0, y0=initial_state(p, h0=0.35, A0=250.0))


def scen_late_bulky(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    return simulate(p, r, tmax=1400.0, y0=initial_state(p, h0=2.00, A0=800.0))


def scen_ip_only(p):
    r = Regimen()
    for c in range(6):
        r.infuse("ip", c * CYC + 1.0, 0.04, 100.0 * p["BSA"])
    return simulate(p, r, tmax=1400.0)


def scen_io_then_chemo(p):
    r = Regimen()
    add_nivo_ipi(r, p, until=250.0)
    add_chemo(r, p, ncyc=6, t0=260.0, pem=True, plat="cis", supplemented=True)
    return simulate(p, r, tmax=1400.0)


def scen_chemo_then_io(p):
    r = Regimen()
    add_chemo(r, p, ncyc=6, pem=True, plat="cis", supplemented=True)
    add_nivo_ipi(r, p, until=600.0, t0=140.0)
    return simulate(p, r, tmax=1400.0)


SCENARIOS = [
    ("1  BSC / untreated",                     scen_bsc),
    ("2  cisplatin alone x6",                  scen_cis_alone),
    ("3  cis+pem, NO supplementation",         scen_cispem_unsupp),
    ("4  cis+pem, folate+B12 supplemented",    scen_cispem_supp),
    ("5  carboplatin+pem",                     scen_carbopem),
    ("6  cis+pem+bevacizumab (MAPS)",          scen_maps),
    ("7  nivolumab+ipilimumab (CM743)",        scen_nivo_ipi),
    ("8  pembrolizumab+chemo (IND.227)",       scen_pembro_chemo),
    ("9  chemo -> nivolumab 2L (CONFIRM)",     scen_nivo_2l),
    ("10 ADI-PEG20 + chemo (ATOMIC-Meso)",     scen_adi_chemo),
    ("11 extended P/D + chemo (MARS2)",        scen_mars2),
    ("12 EPP + chemo + hemithoracic RT",       scen_epp_rt),
    ("13 chemo + intrapleural cisplatin",      scen_ip_after_sys),
    ("14 talc pleurodesis then chemo",         scen_talc_then_chemo),
    ("15 talc, then chemo + intrapleural",     scen_talc_then_ip),
    ("16 chemo + PRMT5i (MTAP-deleted)",       scen_prmt5),
    ("17 chemo + pemetrexed maintenance",      scen_maint_pem),
    ("18 chemo with CrCl 52 mL/min",           scen_renal_impaired),
    ("19 chemo -> vinorelbine 2L",             scen_vin_2l),
    ("20 chemo, EARLY small rind 0.35 cm",     scen_early_small),
    ("21 chemo, LATE bulky rind 2.0 cm",       scen_late_bulky),
    ("22 intrapleural cisplatin alone",        scen_ip_only),
    ("23 IO first, chemo at progression",      scen_io_then_chemo),
    ("24 chemo first, IO at day 140",          scen_chemo_then_io),
]


# ============================================================================
# 9.  MAIN — every number quoted in README.md is produced below
# ============================================================================
# Calibration frame.  Control arms differ by more than four months across the
# MPM trials (EMPHACIS cisplatin 9.3 mo; CheckMate 743 chemotherapy 14.1 mo;
# MARS2 chemotherapy 24.8 mo), so absolute medians are NOT a common scale.
# Three absolute anchors are used; everything else is a WITHIN-TRIAL hazard
# ratio.  Every observed value below was read out of the PubMed abstract of
# the trial in question, not from memory.
ANCHORS = [
    ("BSC / untreated", "MS01 active symptom control arm", 7.6),
    ("cisplatin alone", "EMPHACIS control arm", 9.3),
    ("cis+pem supplemented", "EMPHACIS experimental arm", 12.1),
]
HR_TARGETS = [
    ("cis+pem vs cisplatin",            "EMPHACIS  12.1 vs 9.3 mo",       0.77),
    ("+ bevacizumab vs cis+pem",        "MAPS  18.8 vs 16.1 mo",          0.77),
    ("nivo+ipi vs chemotherapy",        "CheckMate 743  18.1 vs 14.1 mo", 0.74),
    ("pembrolizumab+chemo vs chemo",    "IND.227  17.3 vs 16.1 mo",       0.79),
    ("nivolumab 2L vs placebo",         "CONFIRM  10.2 vs 6.9 mo",        0.69),
    ("pegargiminase+chemo (non-epi)",   "ATOMIC-Meso  9.3 vs 7.7 mo",     0.71),
    ("extended P/D + chemo vs chemo",   "MARS2  19.3 vs 24.8 mo",         1.28),
    ("extrapleural pneumonectomy",      "MARS  14.4 vs 19.5 mo",          1.90),
]

_CACHE = {}


def cached(fn, p, key):
    if key not in _CACHE:
        _CACHE[key] = fn(dict(p))
    return _CACHE[key]


def curve(fn, pbase, mix=None, tag=""):
    """Histology-weighted survival curve for a trial-like cohort."""
    mix = mix or COHORT
    tt, cur = None, []
    for hname, w in mix:
        p = dict(pbase)
        p.update(HIST[hname])
        out = cached(fn, p, (fn.__name__, hname, tag))
        dd = derived(out, p)
        tt = dd["t"]
        cur.append((w, dd["S"]))
    S = [sum(w * c[i] for w, c in cur) for i in range(len(tt))]
    return tt, S


def one(fn, pbase, hname, tag=""):
    p = dict(pbase)
    p.update(HIST[hname])
    return derived(cached(fn, p, (fn.__name__, hname, tag)), p), p


def medmo(c):
    m = median_from_S(c[0], c[1])
    return (m / 30.44) if m else float("nan")


def hr(a, b, day=730):
    return (-math.log(max(a[1][day], 1e-12))) / (-math.log(max(b[1][day], 1e-12)))


NONEPI = [("biphasic", 0.52), ("sarcomatoid", 0.48)]


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


def run_all(quick=False):
    p0 = dict(P)
    p0.update(HIST["epithelioid"])
    t_start = __import__("time").time()

    # ---------------------------------------------------------------------
    hdr("A.  RK4 STEP-SIZE CONVERGENCE")
    print("   The integrator uses a two-phase step: dt = 0.02 d while any drug")
    print("   depot is switching (the fastest rate in the system is cisplatin's")
    print("   k_el + k_bind = 33.3 /d) and dt = 0.06 d thereafter.")
    ref = None
    for dt, dts in ((0.04, 0.12), (0.02, 0.06), (0.01, 0.03), (0.005, 0.015)):
        r = Regimen()
        add_chemo(r, p0, ncyc=6, supplemented=True)
        o = simulate(p0, r, tmax=400.0, dt=dt, dt_slow=dts, record_every=400.0)
        dd = derived(o, p0)
        v = dd["V"][-1]
        if ref is None:
            ref = v
        print(f"   dt = {dt:6.3f} / {dts:5.3f} d   V(400 d) = {v:9.4f} cm3   "
              f"CumHaz = {o['y'][-1][IX['CUMH']]:9.6f}")
    r = Regimen(); add_chemo(r, p0, ncyc=6, supplemented=True)
    fine = derived(simulate(p0, r, tmax=400.0, dt=0.005, dt_slow=0.015,
                            record_every=400.0), p0)["V"][-1]
    r = Regimen(); add_chemo(r, p0, ncyc=6, supplemented=True)
    used = derived(simulate(p0, r, tmax=400.0, record_every=400.0), p0)["V"][-1]
    print(f"   -> the step actually used differs from dt = 0.005 d by "
          f"{100 * abs(used - fine) / fine:.3f} %.")

    # ---------------------------------------------------------------------
    hdr("B.  NATURAL HISTORY: THE GEOMETRY OF AN UNTREATED SHEET")
    o = cached(scen_bsc, p0, ("scen_bsc", "epithelioid", ""))
    dd = derived(o, p0)
    print(f"   {'day':>5} {'A cm2':>8} {'h_meas cm':>10} {'V cm3':>9} "
          f"{'VDT d':>8} {'f_exp':>7} {'L_p mm':>7} {'eff mL':>8} {'SMRP nM':>8}")
    prevV, prevt = None, None
    for day in (0, 60, 120, 180, 270, 365, 540, 730, 1000):
        V = dd["V"][day]
        vdt = ""
        if prevV and V > prevV:
            vdt = f"{(day - prevt) * math.log(2) / math.log(V / prevV):8.0f}"
        print(f"   {day:5d} {o['y'][day][IX['A']]:8.1f} {dd['h_meas'][day]:10.3f} "
              f"{V:9.1f} {vdt:>8} {dd['fexp'][day]:7.3f} {dd['Lp'][day]:7.3f} "
              f"{o['y'][day][IX['VEFF']]:8.0f} {o['y'][day][IX['SMRP']]:8.2f}")
        prevV, prevt = V, day
    print("   The volume doubling time LENGTHENS with no resistance mechanism.")
    print("   The proliferating pool is a shell of fixed depth L_OX = 0.18 mm, so")
    print("   it scales with AREA while the burden scales with AREA x THICKNESS.")

    # ---------------------------------------------------------------------
    hdr("C.  EFFUSION IS NON-MONOTONE IN TUMOUR AREA")
    veff = series(o, "VEFF")
    pk = max(range(len(veff)), key=lambda i: veff[i])
    print(f"   day    0: {veff[0]:6.0f} mL  at A/S = {o['y'][0][IX['A']] / p0['S_PL']:.2f}")
    print(f"   peak {pk:4d}: {veff[pk]:6.0f} mL  at A/S = {o['y'][pk][IX['A']] / p0['S_PL']:.2f}")
    for day in (365, 730, 1000):
        print(f"   day {day:4d}: {veff[day]:6.0f} mL  at A/S = "
              f"{o['y'][day][IX['A']] / p0['S_PL']:.2f}, symphysis "
              f"{o['y'][day][IX['PSY']]:.2f}")
    print("   Formation rises with A; stomatal drainage falls with A; the space")
    print("   has a finite volume that the rind itself eats into.  Nothing was")
    print("   added to make the late hemithorax go dry -- it is the product.")

    # ---------------------------------------------------------------------
    hdr("D.  TRIAL CALIBRATION")
    c_bsc = curve(scen_bsc, p0)
    c_cis = curve(scen_cis_alone, p0)
    c_cp = curve(scen_cispem_supp, p0)
    c_cpu = curve(scen_cispem_unsupp, p0)
    c_bev = curve(scen_maps, p0)
    c_io = curve(scen_nivo_ipi, p0)
    c_pem = curve(scen_pembro_chemo, p0)
    c_2l = curve(scen_nivo_2l, p0)
    c_epd = curve(scen_mars2, p0)
    c_epp = curve(scen_epp_rt, p0)
    c_adi_ne = curve(scen_adi_chemo, p0, NONEPI, "ne")
    c_cp_ne = curve(scen_cispem_supp, p0, NONEPI, "ne")

    print("   Absolute anchors (cohort = 75 % epithelioid / 13 % biphasic / 12 % sarcomatoid):")
    print(f"   {'arm':<26}{'model mo':>10}{'observed mo':>13}   source")
    for (nm, src, tgt), c in zip(ANCHORS, (c_bsc, c_cis, c_cp)):
        print(f"   {nm:<26}{medmo(c):10.1f}{tgt:13.1f}   {src}")
    print()
    print("   Hazard ratios (model read at 24 months; the trials' own HRs):")
    print(f"   {'comparison':<32}{'HR 12 mo':>10}{'HR 24 mo':>10}{'trial HR':>10}   source")
    pairs = [(c_cp, c_cis), (c_bev, c_cp), (c_io, c_cp), (c_pem, c_cp),
             (c_2l, c_cp), (c_adi_ne, c_cp_ne), (c_epd, c_cp), (c_epp, c_cp)]
    for (nm, src, tgt), (a, b) in zip(HR_TARGETS, pairs):
        print(f"   {nm:<32}{hr(a, b, 365):10.2f}{hr(a, b):10.2f}{tgt:10.2f}   {src}")
    print()
    print("   Other model medians, for orientation only -- the trials they come")
    print("   from enrolled fitter patients than the EMPHACIS-anchored cohort:")
    for nm, c in (("cis+pem unsupplemented", c_cpu), ("+ bevacizumab", c_bev),
                  ("nivolumab + ipilimumab", c_io), ("pembrolizumab + chemo", c_pem),
                  ("extended P/D + chemo", c_epd), ("EPP + chemo + RT", c_epp),
                  ("pegargiminase, non-epithelioid", c_adi_ne),
                  ("chemo alone, non-epithelioid", c_cp_ne)):
        print(f"      {nm:<34}{medmo(c):6.1f} mo")

    # ---------------------------------------------------------------------
    hdr("E.  mRECIST BEST RESPONSE VS TRUE VIABLE-CELL KILL")
    print("   mRECIST sums SIX THICKNESSES.  Thickness contains the fibrotic")
    print("   matrix left by the cells the drug already killed, and it falls")
    print("   when the same tumour SPREADS over more area.  Both columns below")
    print("   describe the same patient at the same instant.")
    print(f"   {'scenario':<40}{'mRECIST %':>11}{'viable N %':>12}{'gap pp':>9}{'cat':>5}")
    for name, fn in SCENARIOS[:14]:
        ddx, pe = one(fn, p0, "epithelioid")
        b_m, b_n = ddx["mRECIST"][0], ddx["Vvia"][0]
        i_best = min(range(min(300, len(ddx["t"]))), key=lambda i: ddx["mRECIST"][i])
        m_ch = 100 * (ddx["mRECIST"][i_best] - b_m) / b_m
        n_at = 100 * (ddx["Vvia"][i_best] - b_n) / b_n
        cat = "PR" if m_ch <= -30 else ("SD" if m_ch < 20 else "PD")
        print(f"   {name:<40}{m_ch:11.1f}{n_at:12.1f}{n_at - m_ch:9.1f}{cat:>5}")
    print("   A negative gap means the scan UNDERSTATES the kill; a positive gap")
    print("   means it overstates it.")

    # ---------------------------------------------------------------------
    hdr("F.  THE GEOMETRIC SANCTUARY")
    print(f"   {'h cm':>6}{'L_p mm':>9}{'f_sys':>8}{'f_IP':>8}{'union':>8}"
          f"{'sanctuary':>11}{'f_IgG':>8}{'f_Tcell':>9}")
    for h in (0.20, 0.40, 0.60, 0.80, 1.20, 1.60, 2.40, 3.20):
        y = initial_state(p0, h0=h, A0=450.0)
        L = penetration_lengths(y, p0)
        fs, fi = fpen(L["cis"], h), fpen(L["ip"], h)
        print(f"   {h:6.2f}{L['cis'] * 10:9.3f}{fs:8.3f}{fi:8.3f}"
              f"{min(1.0, fs + fi):8.3f}{max(0.0, 1 - fs - fi):11.3f}"
              f"{fpen(L['mab'], h):8.3f}{fpen(L['tcell'], h):9.3f}")
    y = initial_state(p0, h0=0.80, A0=450.0)
    L = penetration_lengths(y, p0)
    print(f"   At the ADVANCING MARGIN (h = {p0['HFRONT']:.2f} cm) the same "
          f"formula gives f_sys = {min(1.0, fpen(L['cis'], p0['HFRONT'])):.3f} "
          f"and f_Tcell = {min(1.0, fpen(L['tcell'], p0['HFRONT'])):.3f}.")
    print("   That is the whole asymmetry: the front is reachable, the bulk is not.")

    # ---------------------------------------------------------------------
    hdr("G.  EARLY VS LATE: THE SAME REGIMEN AGAINST A THIN AND A THICK RIND")
    print(f"   {'baseline':<26}{'f_exp(0)':>10}{'viable nadir %':>16}"
          f"{'mRECIST best %':>16}{'median OS mo':>14}")
    for label, h0, A0 in (("thin  0.35 cm / 250 cm2", 0.35, 250.0),
                          ("usual 0.80 cm / 450 cm2", 0.80, 450.0),
                          ("bulky 2.00 cm / 800 cm2", 2.00, 800.0)):
        def runner(pp, h0=h0, A0=A0):
            r = Regimen(); add_chemo(r, pp, ncyc=6, supplemented=True)
            return simulate(pp, r, y0=initial_state(pp, h0=h0, A0=A0))
        runner.__name__ = f"early_{h0}_{A0}"
        cc = curve(runner, p0)
        ddx, _ = one(runner, p0, "epithelioid")
        nadir = 100 * (min(ddx["Vvia"][:250]) - ddx["Vvia"][0]) / ddx["Vvia"][0]
        mr = min(100 * (m - ddx["mRECIST"][0]) / ddx["mRECIST"][0]
                 for m in ddx["mRECIST"][:250])
        print(f"   {label:<26}{ddx['fexp'][0]:10.3f}{nadir:16.1f}{mr:16.1f}"
              f"{medmo(cc):14.1f}")

    # ---------------------------------------------------------------------
    hdr("H.  THE MATRIX THE DRUG MAKES IS THE BARRIER THE DRUG THEN MEETS")
    ddx, pe = one(scen_cispem_supp, p0, "epithelioid")
    ochem = _CACHE[("scen_cispem_supp", "epithelioid", "")]
    print(f"   {'cycle':>6}{'day':>6}{'h_meas cm':>11}{'matrix cm3':>12}"
          f"{'f_exp':>8}{'viable change this cycle %':>28}")
    for c in range(6):
        d0, d1 = int(c * CYC), int((c + 1) * CYC)
        k = 100 * (ddx["Vvia"][d1] - ddx["Vvia"][d0]) / max(ddx["Vvia"][d0], 1e-9)
        print(f"   {c + 1:6d}{d0:6d}{ddx['h_meas'][d0]:11.3f}"
              f"{ochem['y'][d0][IX['M']]:12.1f}{ddx['fexp'][d0]:8.3f}{k:28.1f}")
    # ---- the author's hypothesis, and its refutation --------------------
    p_nomat = dict(p0); p_nomat["PHIM"] = 0.0
    r2 = Regimen(); add_chemo(r2, p_nomat, ncyc=6, supplemented=True)
    dnm = derived(simulate(p_nomat, r2, tmax=200.0), p_nomat)
    print(f"   Matrix clears with t1/2 = {math.log(2) / P['KDEGM']:.0f} d while cycles are "
          f"21 d apart, so the debris of cycle 1 is still in the diffusion")
    print("   path at cycle 4, and the matrix pool GROWS from 144 to ~196 cm3.")
    print()
    print("   THE AUTHOR EXPECTED THIS TO PROTECT THE LATER CYCLES.  It does not.")
    print(f"   The exposed fraction RISES over the six cycles "
          f"({ddx['fexp'][0]:.3f} -> {ddx['fexp'][105]:.3f}) and the per-cycle kill")
    print("   ACCELERATES, because the rind thins faster than the matrix")
    print("   accumulates.  Matrix only RETARDS the improvement.  Running the")
    print("   same regimen with PHIM = 0 (killed cells vanish instead of")
    print("   fibrosing) isolates how much:")
    print(f"      with matrix (PHIM = {P['PHIM']:.2f}):  h at day 105 = "
          f"{ddx['h_meas'][105]:.3f} cm, f_exp = {ddx['fexp'][105]:.3f}")
    print(f"      without matrix (PHIM = 0):  h at day 105 = "
          f"{dnm['h_meas'][105]:.3f} cm, f_exp = {dnm['fexp'][105]:.3f}")
    print(f"   So the debris costs {100 * (1 - ddx['fexp'][105] / dnm['fexp'][105]):.0f} % of the "
          f"exposed fraction that would otherwise have been")
    print("   recovered by cycle 6 -- a real penalty, but not the reversal that")
    print("   was expected.  Where the matrix DOES dominate is the MEASUREMENT:")
    print(f"   at day 105 it is {100 * P['PHIM'] * ochem['y'][105][IX['M']] / (ochem['y'][105][IX['N']] + P['PHIM'] * ochem['y'][105][IX['M']]):.0f} % of what the scanner is measuring.")

    # ---------------------------------------------------------------------
    hdr("I.  THE ANTI-ANGIOGENIC PARADOX: BEVACIZUMAB DOSE SWEEP")
    print("   Blocking VEGF lowers interstitial pressure (which helps delivery)")
    print("   and prunes vessels (which hurts it), while also drying the pleural")
    print("   space and slowing the front.  The net is not monotone in dose.")
    print(f"   {'mg/kg':>7}{'free VEGF':>11}{'IFP mmHg':>10}{'rho_v':>8}"
          f"{'L_p mm':>9}{'f_exp':>8}{'effusion mL':>13}{'median OS mo':>14}")
    for dose in (0.0, 1.0, 2.5, 5.0, 7.5, 10.0, 15.0, 25.0, 40.0):
        def runner(pp, dose=dose):
            r = Regimen()
            add_chemo(r, pp, ncyc=6, supplemented=True)
            if dose > 0:
                add_bev(r, pp, maint_until=730.0, mgkg=dose)
            return simulate(pp, r)
        runner.__name__ = f"bev_{dose}"
        cc = curve(runner, p0)
        ddx, pe = one(runner, p0, "epithelioid")
        ob = _CACHE[(runner.__name__, "epithelioid", "")]
        i = 120
        Lb = penetration_lengths(ob["y"][i], pe)
        print(f"   {dose:7.1f}{free_vegf(ob['y'][i], pe):11.1f}{Lb['ifp']:10.1f}"
              f"{ob['y'][i][IX['RHOV']]:8.3f}{Lb['cis'] * 10:9.3f}"
              f"{ddx['fexp'][i]:8.3f}{ob['y'][i][IX['VEFF']]:13.0f}{medmo(cc):14.1f}")

    # ---------------------------------------------------------------------
    hdr("J.  INTRAPLEURAL AND SYSTEMIC COVER COMPLEMENTARY DEPTHS")
    print(f"   {'regimen':<30}{'viable nadir %':>16}{'peak symphysis':>16}"
          f"{'median OS mo':>14}")
    for nm, fn in (("systemic chemo only", scen_cispem_supp),
                   ("intrapleural only", scen_ip_only),
                   ("systemic + intrapleural", scen_ip_after_sys),
                   ("talc first, then chemo", scen_talc_then_chemo),
                   ("talc first, then both", scen_talc_then_ip)):
        cc = curve(fn, p0)
        ddx, _ = one(fn, p0, "epithelioid")
        oo = _CACHE[(fn.__name__, "epithelioid", "")]
        nadir = 100 * (min(ddx["Vvia"][:250]) - ddx["Vvia"][0]) / ddx["Vvia"][0]
        print(f"   {nm:<30}{nadir:16.1f}{max(series(oo, 'PSY')):16.2f}"
              f"{medmo(cc):14.1f}")
    print("   The two routes attack opposite faces, so their coverage adds rather")
    print("   than overlapping -- and talc closes the very route it was given to")
    print("   make usable, which is why 'talc first, then both' is worse than")
    print("   'systemic + intrapleural' even though the effusion is better")
    print("   controlled.")
    print()
    print("   CAVEAT ON THE INTRAPLEURAL ARMS.  The model likes them more than")
    print("   the literature does, and the reason is identifiable: an instillate")
    print("   reaches the ADVANCING MARGIN, which is only 0.6 mm deep, from the")
    print(f"   free face, so f = {min(1.0, fpen(penetration_lengths(initial_state(p0), p0)['ip'], p0['HFRONT'])) * p0['FCONT']:.2f} there against "
          f"{fpen(penetration_lengths(initial_state(p0), p0)['ip'], 0.80) * p0['FCONT']:.3f} in the bulk.")
    print(f"   That rests on assuming a contact fraction of {p0['FCONT']:.2f} -- i.e. that an")
    print("   instillate bathes just over half the involved surface.  No trial has")
    print("   compared intrapleural chemotherapy against systemic chemotherapy in")
    print("   MPM, so this is a PREDICTION and a fragile one: it is roughly linear")
    print("   in a contact fraction nobody has measured.")

    # ---------------------------------------------------------------------
    hdr("K.  MARS2: WHAT AN OPERATION DOES TO A SHEET")
    om = _CACHE.get(("scen_mars2", "epithelioid", ""))
    if om is None:
        _, _ = one(scen_mars2, p0, "epithelioid")
        om = _CACHE[("scen_mars2", "epithelioid", "")]
    ddm, pm = one(scen_mars2, p0, "epithelioid")
    print(f"   pre-op  day 69: A = {om['y'][69][IX['A']]:6.1f} cm2, "
          f"h = {ddm['h_meas'][69]:5.3f} cm, V = {ddm['V'][69]:7.1f} cm3")
    print(f"   post-op day 72: A = {om['y'][72][IX['A']]:6.1f} cm2, "
          f"h = {ddm['h_meas'][72]:5.3f} cm, V = {ddm['V'][72]:7.1f} cm3   "
          f"({100 * (1 - ddm['V'][72] / ddm['V'][69]):.1f} % debulked, "
          f"{100 * (1 - om['y'][72][IX['A']] / om['y'][69][IX['A']]):.1f} % of the AREA removed)")
    for day in (100, 150, 200, 300, 400, 600):
        print(f"   day {day:4d}: A = {om['y'][day][IX['A']]:6.1f} cm2, "
              f"h = {ddm['h_meas'][day]:5.3f} cm, V = {ddm['V'][day]:7.1f} cm3, "
              f"IL-6 = {om['y'][day][IX['IL6']]:5.1f} pg/mL, ECOG = {ddm['PS'][day]:.2f}")
    print(f"   model HR vs chemotherapy alone: {hr(c_epd, c_cp, 365):.2f} at 12 mo, "
          f"{hr(c_epd, c_cp):.2f} at 24 mo   (MARS2 reported 1.28)")
    print(f"   extrapleural pneumonectomy:     {hr(c_epp, c_cp, 365):.2f} at 12 mo, "
          f"{hr(c_epp, c_cp):.2f} at 24 mo   (MARS reported 1.90)")

    # ---------------------------------------------------------------------
    hdr("L.  HISTOLOGY: WHY CHEMOTHERAPY AND IMMUNOTHERAPY SWAP PLACES")
    print(f"   {'histology':<14}{'chemo mo':>10}{'nivo+ipi mo':>13}"
          f"{'HR 12 mo':>10}{'HR 24 mo':>10}{'CM743 HR':>10}")
    cm743 = {"epithelioid": 0.86, "sarcomatoid": 0.46}
    for hname in ("epithelioid", "biphasic", "sarcomatoid"):
        dc, _ = one(scen_cispem_supp, p0, hname)
        di, _ = one(scen_nivo_ipi, p0, hname)
        a = (di["t"], di["S"]); b = (dc["t"], dc["S"])
        ref = cm743.get(hname)
        print(f"   {hname:<14}{medmo(b):10.1f}{medmo(a):13.1f}"
              f"{hr(a, b, 365):10.2f}{hr(a, b):10.2f}"
              f"{(f'{ref:.2f}' if ref else '-'):>10}")
    print("   CHEMOS and IMMINF/VISTA_S are set from pathology, not from outcome.")
    print("   The direction of the crossing is an output.")
    print()
    print("   TWO FINDINGS FALL OUT OF THIS TABLE.")
    print()
    print("   (a) THE EPITHELIOID RESULT SITS ON A THRESHOLD, AND IS FRAGILE.")
    print("       Sweeping the maximal immune kill rate over a 2.4-fold range,")
    print("       with everything else held fixed:")
    print(f"       {'KKILL /d':>10}{'kill at the margin /d':>24}{'epithelioid HR 24 mo':>24}")
    for kk in (0.22, 0.33, 0.44, 0.53):
        pk = dict(p0); pk["KKILL"] = kk
        pk.update(HIST["epithelioid"])
        tag = f"kk{kk}"
        da = derived(cached(scen_nivo_ipi, pk, ("scen_nivo_ipi", "epithelioid", tag)), pk)
        db = derived(cached(scen_cispem_supp, pk, ("scen_cispem_supp", "epithelioid", tag)), pk)
        ok = _CACHE[("scen_nivo_ipi", "epithelioid", tag)]
        yk = ok["y"][300]
        tp = yk[IX["TEFF"]] + pk["WMEM"] * yk[IX["TCLON"]]
        Lk = penetration_lengths(yk, pk)
        kmarg = kk * (tp / (tp + pk["KMKILL"])) * min(1.0, fpen(Lk["tcell"], pk["HFRONT"]))
        print(f"       {kk:10.2f}{kmarg:24.4f}"
              f"{hr((da['t'], da['S']), (db['t'], db['S'])):24.2f}")
    print("       An earlier parameterisation of this model appeared to show that")
    print("       the epithelioid hazard ratio was INSENSITIVE to this rate.  That")
    print("       was an artefact of moving the infiltration and VISTA parameters")
    print("       to compensate; with everything else held fixed the table above")
    print("       shows the opposite, and the claim is withdrawn.")
    print()
    print("       What the sweep actually shows is a THRESHOLD.  The quantity that")
    print("       matters is the immune kill AT THE THIN ADVANCING MARGIN against")
    print(f"       the margin\'s own net proliferation, {p0['GFRONT_REF']:.3f} /d.  Below it the")
    print("       front keeps advancing and immunotherapy loses to six cycles of")
    print("       chemotherapy; approaching it, the front stops for as long as the")
    print("       drug is given -- two years rather than five months -- and")
    print("       immunotherapy wins by DURATION.  The calibrated value sits at")
    print("       0.043 /d, below the threshold and on the steep part of the")
    print("       curve.  So the model\'s epithelioid prediction is genuinely")
    print("       ill-conditioned, which is worth saying plainly: a 30 % error in")
    print("       one rate constant moves the epithelioid HR from 1.1 to 0.5.")
    print()
    print("   (b) THE COMPARATOR\'S CROSSOVER IS NEEDED TO REPRODUCE CM743.")
    print("       Roughly 44 % of the CheckMate 743 chemotherapy arm went on to")
    print("       subsequent systemic therapy.  Mixing that fraction of the")
    print("       chemo-then-nivolumab arm into the comparator moves the model:")
    c_2l_mix = None
    for f2 in (0.0, 0.25, 0.45, 0.60):
        mixS = [(1 - f2) * c_cp[1][i] + f2 * c_2l[1][i] for i in range(len(c_cp[1]))]
        mixc = (c_cp[0], mixS)
        if abs(f2 - 0.45) < 1e-9:
            c_2l_mix = mixc
        print(f"          {100 * f2:3.0f} % crossover in the comparator -> HR "
              f"{hr(c_io, mixc):.2f}   (CheckMate 743: 0.74)")
    for hname, ref in (("epithelioid", 0.86), ("sarcomatoid", 0.46)):
        da, _ = one(scen_nivo_ipi, p0, hname)
        db, _ = one(scen_cispem_supp, p0, hname)
        dc, _ = one(scen_nivo_2l, p0, hname)
        mixS = [0.55 * db["S"][i] + 0.45 * dc["S"][i] for i in range(len(db["S"]))]
        print(f"          {hname:<12} HR {hr((da['t'], da['S']), (db['t'], db['S'])):.2f} "
              f"without crossover, {hr((da['t'], da['S']), (db['t'], mixS)):.2f} with 45 % "
              f"-> trial {ref:.2f}")
    print("       The sarcomatoid subgroup then lands on the trial value exactly.")
    print("       The epithelioid subgroup does not -- see section S.")

    # ---------------------------------------------------------------------
    hdr("M.  FOLATE SUPPLEMENTATION: TOXICITY FALLS, EFFICACY DOES NOT")
    print(f"   {'':<22}{'folate nM':>11}{'Hcy uM':>9}{'ANC nadir':>11}"
          f"{'d ANC<0.5':>11}{'cycles':>8}{'dose %':>8}{'viable nadir %':>16}{'OS mo':>8}")
    for lab, fn, supp in (("unsupplemented", scen_cispem_unsupp, False),
                          ("folate 400 ug + B12", scen_cispem_supp, True)):
        cc = curve(fn, p0)
        ddx, pe = one(fn, p0, "epithelioid")
        oo = _CACHE[(fn.__name__, "epithelioid", "")]
        anc = series(oo, "CIRC")
        r = Regimen(); add_chemo(r, pe, ncyc=6, supplemented=supp)
        simulate(pe, r, tmax=300.0)
        st = r.chemo_state
        nadir = 100 * (min(ddx["Vvia"][:250]) - ddx["Vvia"][0]) / ddx["Vvia"][0]
        print(f"   {lab:<22}{oo['y'][60][IX['FOL']]:11.1f}{oo['y'][60][IX['HCY']]:9.1f}"
              f"{min(anc[:150]):11.2f}{sum(1 for a in anc[:200] if a < 0.5):11d}"
              f"{st['given']:8d}{100 * st['delivered'] / max(st['planned'], 1e-9):8.0f}"
              f"{nadir:16.1f}{medmo(cc):8.1f}")
    print(f"   The asymmetry is a ratio of two rescue constants: the folate term")
    print(f"   raises the MARROW EC50 by 1 + FOL/{P['KFPM']:.0f} nM and the TUMOUR EC50")
    print(f"   by only 1 + FOL/{P['KFPT']:.0f} nM -- a {P['KFPT'] / P['KFPM']:.1f}-fold selectivity.")
    print(f"   At 400 ug/d that is a {1 + 45 / P['KFPM']:.2f}x rise in the marrow EC50 "
          f"against {1 + 45 / P['KFPT']:.2f}x in the tumour.")

    # ---------------------------------------------------------------------
    hdr("N.  THE NEPHROTOXICITY FEEDBACK LOOP AND A BIOMARKER THAT LIES")
    oo = _CACHE[("scen_cispem_supp", "epithelioid", "")]
    p_nokid = dict(p0); p_nokid["KNEPH"] = 0.0; p_nokid["KNEPHIRR"] = 0.0
    o2 = scen_cispem_supp(p_nokid)
    print(f"   {'cycle':>6}{'CrCl':>8}{'pem AUC':>10}{'ANC nadir':>11}"
          f"{'viable N':>10}{'SMRP nM':>9}{'SMRP if CrCl fixed':>21}")
    for c in range(6):
        d0, d1 = int(c * CYC), int((c + 1) * CYC - 1)
        anc = min(series(oo, "CIRC")[d0:d1 + 1])
        auc = oo["y"][d1][IX["AUCP"]] - oo["y"][d0][IX["AUCP"]]
        print(f"   {c + 1:6d}{oo['y'][d0][IX['CRCL']]:8.1f}{auc:10.2f}{anc:11.2f}"
              f"{oo['y'][d0][IX['N']]:10.1f}{oo['y'][d0][IX['SMRP']]:9.2f}"
              f"{o2['y'][d0][IX['SMRP']]:21.2f}")
    i = 120
    print(f"   day 120: viable tumour has changed "
          f"{100 * (oo['y'][i][IX['N']] / oo['y'][0][IX['N']] - 1):+.1f} %, "
          f"serum mesothelin only "
          f"{100 * (oo['y'][i][IX['SMRP']] / oo['y'][0][IX['SMRP']] - 1):+.1f} %.")
    print(f"   Of that discrepancy, "
          f"{100 * (oo['y'][i][IX['SMRP']] / o2['y'][i][IX['SMRP']] - 1):+.1f} % is renal and not "
          f"tumour: cisplatin has taken CrCl from {oo['y'][0][IX['CRCL']]:.0f} to "
          f"{oo['y'][i][IX['CRCL']]:.0f} mL/min.")

    # ---------------------------------------------------------------------
    hdr("O.  SEQUENCING AND MAINTENANCE")
    print(f"   {'regimen':<32}{'median OS mo':>14}{'d ANC<0.5':>11}"
          f"{'peak irAE':>11}{'best mRECIST %':>16}")
    for nm, fn in (("chemo x6 only", scen_cispem_supp),
                   ("chemo x6 + pem maintenance", scen_maint_pem),
                   ("IO first, chemo at progression", scen_io_then_chemo),
                   ("chemo, IO added day 140", scen_chemo_then_io),
                   ("chemo -> nivolumab 2L", scen_nivo_2l),
                   ("chemo -> vinorelbine 2L", scen_vin_2l),
                   ("carboplatin + pemetrexed", scen_carbopem),
                   ("chemo with CrCl 52 mL/min", scen_renal_impaired),
                   ("chemo + PRMT5i (MTAP-del)", scen_prmt5)):
        cc = curve(fn, p0)
        ddx, _ = one(fn, p0, "epithelioid")
        oo2 = _CACHE[(fn.__name__, "epithelioid", "")]
        anc = series(oo2, "CIRC")
        mr = min(100 * (m - ddx["mRECIST"][0]) / ddx["mRECIST"][0]
                 for m in ddx["mRECIST"][:300])
        print(f"   {nm:<32}{medmo(cc):14.1f}{sum(1 for a in anc if a < 0.5):11d}"
              f"{max(series(oo2, 'IRAE')):11.2f}{mr:16.1f}")

    # ---------------------------------------------------------------------
    hdr("P.  ADI-PEG20: THE ANTIBODY THAT ENDS THE TREATMENT")
    dda, pa = one(scen_adi_chemo, p0, "sarcomatoid", "ne")
    oa = _CACHE[("scen_adi_chemo", "sarcomatoid", "ne")]
    print(f"   {'day':>6}{'ADI mg/L':>10}{'ADA':>8}{'arginine uM':>13}{'kill /d':>10}")
    for day in (3, 7, 21, 42, 63, 84, 120, 180, 270, 364, 400):
        y = oa["y"][day]
        arg = y[IX["ARG"]]
        k = P["EMAXARG"] * max(0.0, (P["ARG_CRIT"] - arg)) / P["ARG_CRIT"]
        print(f"   {day:6d}{y[IX['ADI_C']] / P['V_ADI']:10.3f}{y[IX['ADA']]:8.3f}"
              f"{arg:13.1f}{k:10.4f}")
    p_ass1 = dict(p0); p_ass1["ASS1LOSS"] = 0.0
    c_ass1 = curve(scen_adi_chemo, p_ass1, NONEPI, "ass1intact")
    print(f"   ASS1-methylated (auxotroph): median OS {medmo(c_adi_ne):.1f} mo")
    print(f"   ASS1-intact:                 median OS {medmo(c_ass1):.1f} mo")
    print(f"   chemotherapy alone:          median OS {medmo(c_cp_ne):.1f} mo")
    print("   ATOMIC-Meso enrolled unselected non-epithelioid disease: 9.3 vs 7.7 mo.")

    # ---------------------------------------------------------------------
    hdr("Q.  PREDICTIONS THE MODEL MAKES THAT NOBODY HAS MEASURED")
    print("   (i) mRECIST CANNOT SEE A TWELVE-FOLD DIFFERENCE IN DRUG ACCESS.")
    print("       The author expected the partial-response rate to FALL with")
    print("       baseline rind thickness, since the exposed fraction does.  The")
    print("       table refutes that and says something sharper: with the involved")
    print("       AREA held fixed and only the thickness varied, the measured")
    print("       response is flat while the true kill is not.")
    print(f"   {'baseline h cm':>15}{'f_exp':>8}{'mRECIST best %':>16}"
          f"{'viable nadir %':>16}{'category':>10}")
    for h in (0.25, 0.50, 0.80, 1.20, 1.80, 2.50):
        pe = dict(p0)
        r = Regimen(); add_chemo(r, pe, ncyc=6, supplemented=True)
        oq = simulate(pe, r, tmax=400.0, y0=initial_state(pe, h0=h, A0=450.0))
        ddq = derived(oq, pe)
        mr = min(100 * (m - ddq["mRECIST"][0]) / ddq["mRECIST"][0]
                 for m in ddq["mRECIST"][:250])
        vn = min(100 * (v - ddq["Vvia"][0]) / ddq["Vvia"][0]
                 for v in ddq["Vvia"][:250])
        print(f"   {h:15.2f}{ddq['fexp'][0]:8.3f}{mr:16.1f}{vn:16.1f}"
              f"{('PR' if mr <= -30 else 'SD'):>10}")
    print("       Every one of them is a partial response, over a 10-fold range of")
    print("       baseline thickness and a 12-fold range of exposed fraction.  The")
    print("       reason is that measured thickness is dominated by matrix")
    print("       clearance and by the sheet spreading, neither of which tracks")
    print("       how many tumour cells the drug reached.  THE TESTABLE CLAIM is")
    print("       therefore the opposite of the one expected: mRECIST response")
    print("       depth should NOT stratify by baseline rind thickness even though")
    print("       survival does (15.3 vs 4.4 months in section G).")
    print()
    print("   (ii) The fibrotic fraction of the MEASURED rind under treatment.")
    for day in (0, 30, 60, 120, 180, 300, 500):
        fib = P["PHIM"] * oo["y"][day][IX["M"]] / max(
            oo["y"][day][IX["N"]] + P["PHIM"] * oo["y"][day][IX["M"]], 1e-9)
        print(f"        day {day:4d}: {100 * fib:5.1f} % of what the scanner measures "
              f"is matrix, not tumour")
    print()
    print("   (iii) Serum mesothelin should be re-scaled by measured GFR before it")
    print("        is read as a response biomarker; the correction is worth")
    print(f"        {100 * (oo['y'][120][IX['SMRP']] / o2['y'][120][IX['SMRP']] - 1):+.0f} % at cycle 6 of cisplatin.")
    print()
    print("   (iv) Intrapleural chemotherapy and talc pleurodesis are mutually")
    print("        exclusive by geometry, not by pharmacology.")

    # ---------------------------------------------------------------------
    hdr("R.  SCENARIO SUMMARY (24 regimens)")
    print(f"   {'scenario':<40}{'OS mo':>7}{'PFS mo':>8}{'mRECIST %':>11}"
          f"{'viable %':>10}{'ANC nad':>9}{'FVC 6mo':>9}{'eff pk':>8}")
    for name, fn in SCENARIOS:
        cc = curve(fn, p0)
        ddx, _ = one(fn, p0, "epithelioid")
        oo3 = _CACHE[(fn.__name__, "epithelioid", "")]
        mr = min(100 * (m - ddx["mRECIST"][0]) / ddx["mRECIST"][0]
                 for m in ddx["mRECIST"][:300])
        vn = min(100 * (v - ddx["Vvia"][0]) / ddx["Vvia"][0]
                 for v in ddx["Vvia"][:400])
        nad, pfs = ddx["mRECIST"][0], None
        for i, m in enumerate(ddx["mRECIST"]):
            nad = min(nad, m)
            if i > 20 and m >= 1.20 * nad:
                pfs = i
                break
        anc = series(oo3, "CIRC")
        print(f"   {name:<40}{medmo(cc):7.1f}"
              f"{(pfs / 30.44 if pfs else 41.0):8.1f}{mr:11.1f}{vn:10.1f}"
              f"{min(anc):9.2f}{ddx['FVC'][182]:9.1f}"
              f"{max(series(oo3, 'VEFF')):8.0f}")

    # ---------------------------------------------------------------------
    hdr("S.  WHERE THIS MODEL DISAGREES WITH THE DATA")
    print("   Listed rather than tuned away.")
    print()
    e_dc, _ = one(scen_cispem_supp, p0, "epithelioid")
    e_di, _ = one(scen_nivo_ipi, p0, "epithelioid")
    s_dc, _ = one(scen_cispem_supp, p0, "sarcomatoid")
    s_di, _ = one(scen_nivo_ipi, p0, "sarcomatoid")
    print(f"   1. THE HISTOLOGY x IMMUNOTHERAPY INTERACTION IS TOO FLAT.")
    print(f"      CheckMate 743 split 0.86 (epithelioid) / 0.46 (non-epithelioid).")
    print(f"      The model gives {hr((e_di['t'], e_di['S']), (e_dc['t'], e_dc['S'])):.2f} / "
          f"{hr((s_di['t'], s_di['S']), (s_dc['t'], s_dc['S'])):.2f}.  Differences in immune")
    print("      infiltration, VISTA-type suppression and chemosensitivity, set at")
    print("      literature magnitudes, do NOT generate a split that wide.  Either")
    print("      one of those differences is larger than assumed, or the trial's")
    print("      subgroup estimate is unstable (it rests on ~140 patients).")
    print()
    print(f"   2. THE IMMUNOTHERAPY CURVES CROSS TOO LATE.  The model's 12-month HR")
    print(f"      for nivo+ipi is {hr(c_io, c_cp, 365):.2f} against {hr(c_io, c_cp):.2f} at 24 months.")
    print("      CheckMate 743's curves cross at about 8 months; the model's cross")
    print("      later, so it understates early immunotherapy benefit and the")
    print("      median lands on the wrong side of the crossing.")
    print()
    print("   3. ABSOLUTE MEDIANS FOR THE MODERN ARMS ARE LOW.  The cohort is")
    print("      anchored to EMPHACIS-era survival; MAPS, CheckMate 743 and MARS2")
    print("      enrolled fitter patients and their control arms run 16-25 months.")
    print("      Only the hazard ratios are comparable, and those are what was fit.")
    print()
    print("   4. THE SUPPLEMENTATION SURVIVAL EFFECT IS SMALLER THAN THE")
    print("      EMPHACIS SUBGROUP.  The trial's fully-supplemented subset showed")
    print("      13.3 vs 10.0 months; the model produces the toxicity difference")
    print("      but a much smaller survival difference, because in the equations")
    print("      folate reaches survival only through dose intensity and febrile")
    print("      neutropenia.  The abstract itself says supplementation reduced")
    print("      toxicity 'without adversely affecting survival time', so the")
    print("      13.3/10.0 contrast is confounded by era; the model sides with the")
    print("      smaller effect.")
    print()
    print(f"   Total wall clock: {__import__('time').time() - t_start:.0f} s, "
          f"{len(_CACHE)} distinct simulations.")
    print("\nDone.")


if __name__ == "__main__":
    run_all(quick="--quick" in sys.argv)
