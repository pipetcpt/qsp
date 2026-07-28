#!/usr/bin/env python3
# =============================================================================
#  cpp_reference_check.py
#  Central Precocious Puberty (CPP) — independent reference implementation
#  ---------------------------------------------------------------------------
#  A stand-alone, vectorised numpy/RK4 transcription of EXACTLY the same 44
#  ODEs, parameters and dosing regimens that appear in cpp_mrgsolve_model.R.
#  Its purpose is to make every number quoted in README.md reproducible
#  WITHOUT R or mrgsolve, and to let two independent transcriptions of the same
#  equations be cross-checked against each other.
#
#      python3 cpp_reference_check.py        -> cpp_reference_output.txt
#
#  NUMERICAL NOTE.  The gonadotropin and sex-steroid compartments are
#  deliberately COARSE-GRAINED: LH, FSH, E2 and testosterone are represented as
#  pulse-AVERAGED concentrations with relaxation times of 6-12 h rather than
#  their true 20-60 min turnover.  The model resolves days-to-years, not
#  individual GnRH pulses, and pulse FREQUENCY enters explicitly as its own
#  state (PULS) instead.  This keeps a fixed-step RK4 integrator stable and
#  makes the R (lsoda) and python (RK4) solutions agree to <0.1 cm of adult
#  height.  Nothing that matters clinically lives on the sub-6-hour timescale:
#  the agonist flare is set by receptor desensitisation (t1/2 ~1.5 d) and
#  interdose escape by depot dissolution (t1/2 ~15 d), both far slower.
#
#  Nothing here is fitted to individual patient data.
# =============================================================================

import numpy as np

np.seterr(all="ignore")

# -----------------------------------------------------------------------------
#  STATE INDEX  (44 ODEs)
# -----------------------------------------------------------------------------
NAMES = [
    # --- PK / exposure (17) ---
    "LDEPB", "LDEPS", "LA1", "LA2",      # leuprolide depot burst/slow, central, periph
    "TDEPB", "TDEPS", "TA1",             # triptorelin depot burst/slow, central
    "HIMP", "HA1",                       # histrelin implant reservoir, central
    "NDEP", "NA1",                       # nafarelin nasal depot, central
    "XDEP", "XA1",                       # GnRH ANTAGONIST depot, central
    "AIA1",                              # aromatase inhibitor central
    "GHDEP", "GHEFF",                    # rhGH sc depot, effect compartment
    "TAMA1",                             # tamoxifen central
    # --- GnRH receptor (2) ---
    "RS", "RD",
    # --- neuroendocrine (6) ---
    "MKRN3", "KND", "PULS", "LH", "FSH", "INHB",
    # --- gonad / steroid (5) ---
    "FOL", "E2", "TESTO", "DHEAS", "AUTON",
    # --- growth (5) ---
    "GH", "IGF1", "BA", "GPRES", "HT",
    # --- target tissue / body (7) ---
    "UTV", "BST", "ENDO", "BMIZ", "BMDZ", "HF", "E2TRK",
    # --- indices / integrators (2) ---
    "QOL", "CUME2",
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)
assert NST == 44, NST

# -----------------------------------------------------------------------------
#  PARAMETERS  (names/values identical to $PARAM in cpp_mrgsolve_model.R)
# -----------------------------------------------------------------------------
P = dict(
    # ---------------- PK: leuprolide acetate PLGA depot (IM) ----------------
    FBURSTL=0.11, KABL=1.10, KDISL=0.0455,
    CLL=140.0, V1L=60.0, QL=20.0, V2L=40.0,
    EC50L=0.35,
    # ---------------- PK: triptorelin pamoate depot (IM) -------------------
    FBURSTT=0.07, KABT=0.90, KDIST=0.0150,
    CLT=125.0, V1T=55.0, EC50T=0.28,
    # ---------------- PK: histrelin subdermal implant (zero order) ---------
    K0H=65.0, CLH=130.0, V1H=55.0, EC50H=0.12,
    # ---------------- PK: nafarelin nasal spray ----------------------------
    KAN=6.0, FN=0.021, CLN=150.0, V1N=40.0, EC50N=0.55,
    # ---------------- PK: GnRH antagonist (illustrative sc depot) ----------
    KAX=0.030, CLX=95.0, V1X=45.0, EC50X=0.90,
    # ---------------- PK/PD: aromatase inhibitor --------------------------
    CLAI=26.0, VAI=45.0, FAI=0.85, IMAXAI=0.975, IC50AI=3.0,
    # ---------------- PK/PD: rhGH -----------------------------------------
    KAG=1.4, KELG=1.1, VGH=1200.0, EGH_EXO=0.55,
    # ---------------- PK/PD: tamoxifen ------------------------------------
    CLTAM=380.0, VTAM=1100.0, IMAXTAM=0.72, IC50TAM=40.0,
    # ---------------- GnRH receptor dynamics ------------------------------
    KDES=0.45, KDESE=0.02, KREC=0.035, AINT=1.60,
    # ---------------- MKRN3 brake / KNDy pulse generator ------------------
    MK0=1.00, KMK=6.40e-4, KMKSEX=0.75, KMK50=0.35, HMK=4.0,
    KNDMAX=1.15, KNDLES=0.0, KSEC=0.0, TAUKND=60.0, DLK1R=0.10, KLEPB=0.18,
    PULSMIN=6.0, PULSMAX=18.0, TAUP=20.0, PULSREF=12.0,
    # ---------------- gonadotropes ----------------------------------------
    LHB0=0.030, LHMAX=5.00, KSLH=0.62, HLH=3.2, TAULH=0.25,
    FSHB0=0.60, FSHMAX=5.00, KSF=0.45, HFS=1.6, TAUFSH=0.40, KINH=90.0,
    KFB=120.0,
    # ---------------- ovary / steroids ------------------------------------
    KFOL=0.060, KFF=2.60, HFF=3.0, KFOLO=0.030,
    KE2=120.0, KLE=1.20, HLE=1.30, TAUE2=0.50, KPER=1.40, E2FLOOR=0.80,
    KINHB=150.0, KFF2=2.20, TAUI=1.5,
    KAUT=95.0, TAUAUT=40.0, AUTSET=0.0,
    # ---------------- androgens -------------------------------------------
    TB0=4.0, KTG=600.0, KLT=1.60, HLT=1.20, KTA=22.0, TAUT=0.50,
    LEYD=1.0, SEXM=0.0,
    DHM=250.0, DHCA50=11.0, HDH=3.5, TAUDH=60.0,
    # ---------------- GH / IGF-1 ------------------------------------------
    GHB=1.00, EGH=1.80, KGH=25.0, NGH=2.0, TAUGH=1.0,
    KIGF=160.0, TAUIGF=1.0, IGFREF=150.0,
    # ---------------- growth plate / bone age -----------------------------
    GVBASE=5.90, GVSEXM=0.04, PGP=0.45, PIGF=0.50, EE2GV=0.45, KGV=20.0,
    HGV=1.20, ETGV=0.34, KTGV=250.0, HTGV=1.50,
    # bone age is maturation RELATIVE TO THE SAME-AGE NORMAL REFERENCE CHILD,
    # which is what the Greulich-Pyle atlas is; hence a normal child yields
    # dBA/dCA == 1.0 by construction and any deviation is pathology.
    SCBA=1.00, M0=0.30, ME=0.55, KBA=25.0, HBA=1.30,
    MA=0.20, KAND=0.60, MI=0.45, KIG=300.0,
    BAFUS=15.0, BAREF=5.0, XTRA=0.10,
    # prescribed NORMAL reference trajectories (sex-specific) used only as the
    # denominator of the bone-age rate and as the BMD age-reference
    E2N0=1.40, E2NA=54.0, E2NC=11.35, E2NH=11.0,         # girls
    E2N0M=1.60, E2NAM=27.0, E2NCM=12.60, E2NHM=9.0,      # boys
    TN0=5.0, TNA=42.0, TNC=11.40, TNH=8.0,               # girls
    TN0M=6.0, TNAM=480.0, TNCM=12.60, TNHM=8.0,          # boys
    # ---------------- target tissues --------------------------------------
    UTV0=1.20, KUT=10.00, KUTE=25.0, HUT=1.50, TAUUT=30.0,
    KBST=22.0, HBST=2.20, TAUBSTU=90.0, TAUBSTD=240.0,
    KEND=0.030, KENDE=45.0, HEND=2.0, KENDO=0.010, ENDOBLEED=0.503,
    # ---------------- bone mass ------------------------------------------
    KZ=0.55, KBM=25.0, HBM=1.20, KCATCH=0.35, CAVD=0.0, KCAVD=0.30,
    # ---------------- body composition -----------------------------------
    BMIZTGT=0.30, KBMIG=0.25, TAUBMI=365.0,
    # ---------------- symptoms / QoL --------------------------------------
    KHF=6.00, HFSCALE=20.0, TAUHF=7.0, TAUTRK=60.0, TAUQOL=45.0,
)

GIRL_BP = np.array([  # Bayley-Pinneau girls, "average" column, % of adult height
    [6.0, 66.2], [7.0, 69.6], [8.0, 73.0], [9.0, 77.2], [10.0, 80.4],
    [10.5, 82.3], [11.0, 84.4], [11.5, 86.2], [12.0, 88.4], [12.5, 90.6],
    [13.0, 92.2], [13.5, 94.1], [14.0, 95.8], [14.5, 97.4], [15.0, 98.6],
    [15.5, 99.0], [16.0, 99.6], [17.0, 100.0]])
BOY_BP = np.array([
    [7.0, 69.5], [8.0, 72.3], [9.0, 75.2], [10.0, 78.4], [11.0, 81.2],
    [12.0, 84.2], [12.5, 85.8], [13.0, 87.6], [13.5, 90.2], [14.0, 92.7],
    [14.5, 94.8], [15.0, 96.8], [15.5, 97.9], [16.0, 98.6], [16.5, 99.2],
    [17.0, 99.6], [18.0, 100.0]])


def bp_fraction(ba, male=False):
    tab = BOY_BP if male else GIRL_BP
    return np.interp(ba, tab[:, 0], tab[:, 1]) / 100.0


# -----------------------------------------------------------------------------
#  RIGHT-HAND SIDE — y shape (NST, n); every parameter may be scalar or (n,)
# -----------------------------------------------------------------------------
def rhs(t, y, p, ca0, want_aux=False):
    d = np.zeros_like(y)
    CA = ca0 + t / 365.25

    LDEPB, LDEPS, LA1, LA2 = y[0], y[1], y[2], y[3]
    TDEPB, TDEPS, TA1 = y[4], y[5], y[6]
    HIMP, HA1, NDEP, NA1 = y[7], y[8], y[9], y[10]
    XDEP, XA1, AIA1 = y[11], y[12], y[13]
    GHDEP, GHEFF, TAMA1 = y[14], y[15], y[16]
    RS, RD = y[17], y[18]
    MKRN3, KND, PULS, LH, FSH, INHB = y[19], y[20], y[21], y[22], y[23], y[24]
    FOL, E2, TESTO, DHEAS, AUTON = y[25], y[26], y[27], y[28], y[29]
    GH, IGF1, BA, GPRES, HT = y[30], y[31], y[32], y[33], y[34]
    UTV, BST, ENDO, BMIZ, BMDZ, HF, E2TRK = (y[35], y[36], y[37], y[38],
                                             y[39], y[40], y[41])
    QOL = y[42]

    LH = np.maximum(LH, 1e-9)
    FSH = np.maximum(FSH, 1e-9)
    E2 = np.maximum(E2, 1e-9)
    TESTO = np.maximum(TESTO, 1e-9)
    IGF1 = np.maximum(IGF1, 1e-6)

    # ================= 1. PK =================================================
    d[0] = -p["KABL"] * LDEPB
    d[1] = -p["KDISL"] * LDEPS
    d[2] = (p["KABL"] * LDEPB + p["KDISL"] * LDEPS
            - (p["CLL"] / p["V1L"]) * LA1
            - (p["QL"] / p["V1L"]) * LA1 + (p["QL"] / p["V2L"]) * LA2)
    d[3] = (p["QL"] / p["V1L"]) * LA1 - (p["QL"] / p["V2L"]) * LA2
    CLEUP = LA1 / p["V1L"]                                  # ng/mL

    d[4] = -p["KABT"] * TDEPB
    d[5] = -p["KDIST"] * TDEPS
    d[6] = p["KABT"] * TDEPB + p["KDIST"] * TDEPS - (p["CLT"] / p["V1T"]) * TA1
    CTRIP = TA1 / p["V1T"]

    relH = np.where(HIMP > 1.0, p["K0H"], 0.0)              # zero-order implant
    d[7] = -relH
    d[8] = relH - (p["CLH"] / p["V1H"]) * HA1
    CHIST = HA1 / p["V1H"]

    d[9] = -p["KAN"] * NDEP
    d[10] = p["FN"] * p["KAN"] * NDEP - (p["CLN"] / p["V1N"]) * NA1
    CNAF = NA1 / p["V1N"]

    d[11] = -p["KAX"] * XDEP
    d[12] = p["KAX"] * XDEP - (p["CLX"] / p["V1X"]) * XA1
    CANT = XA1 / p["V1X"]

    d[13] = -(p["CLAI"] / p["VAI"]) * AIA1
    CAI = AIA1 / p["VAI"]
    AIEFF = p["IMAXAI"] * CAI / (CAI + p["IC50AI"])

    d[14] = -p["KAG"] * GHDEP
    d[15] = p["KAG"] * GHDEP - p["KELG"] * GHEFF
    CGH = GHEFF / p["VGH"]

    d[16] = -(p["CLTAM"] / p["VTAM"]) * TAMA1
    CTAM = TAMA1 / p["VTAM"]
    TAMEFF = p["IMAXTAM"] * CTAM / (CTAM + p["IC50TAM"])

    # ================= 2. GnRHR occupancy (competitive) & desensitisation ====
    A = (CLEUP / p["EC50L"] + CTRIP / p["EC50T"]
         + CHIST / p["EC50H"] + CNAF / p["EC50N"])
    B = CANT / p["EC50X"]
    den = 1.0 + A + B
    fa = A / den                    # fraction occupied by AGONIST
    fx = B / den                    # fraction occupied by ANTAGONIST
    ffree = np.clip(1.0 - fa - fx, 0.0, 1.0)

    # ================= 3. MKRN3 brake -> KNDy pulse generator ================
    # boys release the brake ~1.8 yr later than girls at the same MK0
    d[19] = -p["KMK"] * (1.0 - (1.0 - p["KMKSEX"]) * p["SEXM"]) * MKRN3
    LEPF = np.clip(1.0 + p["KLEPB"] * BMIZ, 0.55, 1.60)
    RESTR = MKRN3 + p["DLK1R"]
    KNDss = (p["KNDMAX"] * LEPF / (1.0 + (RESTR / p["KMK50"]) ** p["HMK"])
             + p["KNDLES"]
             + p["KSEC"] * BA ** 6 / (BA ** 6 + 10.5 ** 6))
    d[20] = (KNDss - KND) / p["TAUKND"]
    PULSss = p["PULSMIN"] + (p["PULSMAX"] - p["PULSMIN"]) * np.clip(KND, 0.0, 1.2)
    d[21] = (PULSss - PULS) / p["TAUP"]

    # ================= 4. Pituitary: the PULSATILITY DECODER =================
    FBK = 1.0 / (1.0 + (E2 / p["KFB"]) ** 1.5)
    Sendo = KND * (PULS / p["PULSREF"]) ** 0.5 * FBK
    S = RS * (Sendo * ffree + p["AINT"] * fa)

    d[18] = (p["KDES"] * fa + p["KDESE"] * Sendo) * RS - p["KREC"] * RD
    d[17] = -d[18]

    SH = S ** p["HLH"]
    LHss = p["LHB0"] + p["LHMAX"] * SH / (SH + p["KSLH"] ** p["HLH"])
    d[22] = (LHss - LH) / p["TAULH"]
    SF = S ** p["HFS"]
    FSHss = ((p["FSHB0"] + p["FSHMAX"] * SF / (SF + p["KSF"] ** p["HFS"]))
             / (1.0 + INHB / p["KINH"]))
    d[23] = (FSHss - FSH) / p["TAUFSH"]

    # ================= 5. Gonad / steroidogenesis ============================
    fsh_drive = FSH ** p["HFF"] / (FSH ** p["HFF"] + p["KFF"] ** p["HFF"])
    d[25] = p["KFOL"] * fsh_drive * (1.0 - FOL) - p["KFOLO"] * FOL
    d[29] = (p["AUTSET"] - AUTON) / p["TAUAUT"]

    lh_drive = LH ** p["HLE"] / (LH ** p["HLE"] + p["KLE"] ** p["HLE"])
    AROM = FOL * (1.0 - AIEFF)
    E2per = (p["KPER"] * (TESTO / 20.0 + DHEAS / 400.0) * (1.0 - AIEFF)
             * (1.0 + 0.20 * np.maximum(BMIZ, 0.0)))
    # the gonadal-aromatase term is OVARIAN: gate it by sex, so that in boys
    # oestradiol arises only from PERIPHERAL aromatisation of testosterone
    E2ss = (p["E2FLOOR"] + E2per
            + p["KE2"] * (1.0 - p["SEXM"]) * AROM * lh_drive
            + p["KAUT"] * AUTON * (1.0 - AIEFF))
    d[26] = (E2ss - E2) / p["TAUE2"]

    Tss = (p["TB0"]
           + p["KTG"] * (LH ** p["HLT"] / (LH ** p["HLT"] + p["KLT"] ** p["HLT"]))
           * (p["SEXM"] * p["LEYD"] + (1.0 - p["SEXM"]) * 0.08)
           + p["KTA"] * DHEAS / 400.0)
    d[27] = (Tss - TESTO) / p["TAUT"]

    DHss = p["DHM"] * CA ** p["HDH"] / (CA ** p["HDH"] + p["DHCA50"] ** p["HDH"])
    d[28] = (DHss - DHEAS) / p["TAUDH"]

    INHBss = p["KINHB"] * FOL * FSH ** 2 / (FSH ** 2 + p["KFF2"] ** 2)
    d[24] = (INHBss - INHB) / p["TAUI"]

    # ================= 6. GH / IGF-1 — the (+) arm ===========================
    fE2gh = E2 ** p["NGH"] / (E2 ** p["NGH"] + p["KGH"] ** p["NGH"])
    GHss = p["GHB"] * (1.0 + p["EGH"] * fE2gh) + p["EGH_EXO"] * CGH
    d[30] = (GHss - GH) / p["TAUGH"]
    d[31] = (p["KIGF"] * GH - IGF1) / p["TAUIGF"]

    # ================= 7. Growth plate — the (-) arm, IRREVERSIBLE ===========
    GPc = np.clip(GPRES, 0.0, 1.0)
    open_plate = GPc > 1e-6
    fE2gv = E2 ** p["HGV"] / (E2 ** p["HGV"] + p["KGV"] ** p["HGV"])
    fTgv = TESTO ** p["HTGV"] / (TESTO ** p["HTGV"] + p["KTGV"] ** p["HTGV"])
    GV = (p["GVBASE"] * (1.0 + p["GVSEXM"] * p["SEXM"])
          * GPc ** p["PGP"] * (IGF1 / p["IGFREF"]) ** p["PIGF"]
          * (1.0 + p["EE2GV"] * fE2gv + p["ETGV"] * fTgv))
    GV = np.where(open_plate, GV, 0.0)
    d[34] = GV / 365.25

    # ---- the maturation driver, and its NORMAL same-age reference value ----
    hb = p["HBA"]
    fE = lambda x: x ** hb / (x ** hb + p["KBA"] ** hb)
    fA = lambda x: x / (x + p["KAND"])
    fI = lambda x: x / (x + p["KIG"])

    E2plate = E2 * (1.0 - TAMEFF)
    ANDR = DHEAS / 400.0 + TESTO / 400.0
    Rmat = (p["M0"] + p["ME"] * fE(E2plate) + p["MA"] * fA(ANDR)
            + p["MI"] * fI(IGF1))

    E2NORM = np.where(p["SEXM"] > 0.5,
                      p["E2N0M"] + p["E2NAM"] * CA ** p["E2NHM"]
                      / (CA ** p["E2NHM"] + p["E2NCM"] ** p["E2NHM"]),
                      p["E2N0"] + p["E2NA"] * CA ** p["E2NH"]
                      / (CA ** p["E2NH"] + p["E2NC"] ** p["E2NH"]))
    TNORM = np.where(p["SEXM"] > 0.5,
                     p["TN0M"] + p["TNAM"] * CA ** p["TNHM"]
                     / (CA ** p["TNHM"] + p["TNCM"] ** p["TNHM"]),
                     p["TN0"] + p["TNA"] * CA ** p["TNH"]
                     / (CA ** p["TNH"] + p["TNC"] ** p["TNH"]))
    ANDRN = DHss / 400.0 + TNORM / 400.0     # adrenarche is axis-independent
    IGFNORM = p["KIGF"] * p["GHB"] * (1.0 + p["EGH"] * E2NORM ** p["NGH"]
                                      / (E2NORM ** p["NGH"] + p["KGH"] ** p["NGH"]))
    RmatN = (p["M0"] + p["ME"] * fE(E2NORM) + p["MA"] * fA(ANDRN)
             + p["MI"] * fI(IGFNORM))

    fE_BA = fE(E2plate)
    dBA_yr = p["SCBA"] * Rmat / RmatN
    dBA_yr = np.where(open_plate, dBA_yr, 0.0)
    d[32] = dBA_yr / 365.25
    d[33] = np.where(GPRES > 0.0,
                     -(dBA_yr / 365.25) / (p["BAFUS"] - p["BAREF"])
                     * (1.0 + p["XTRA"] * fE_BA), 0.0)

    # ================= 8. Target tissues =====================================
    E2t = E2 * (1.0 - TAMEFF)
    UTVss = p["UTV0"] + p["KUT"] * E2t ** p["HUT"] / (E2t ** p["HUT"] + p["KUTE"] ** p["HUT"])
    d[35] = (UTVss - UTV) / p["TAUUT"]

    BSTss = 1.0 + 4.0 * E2t ** p["HBST"] / (E2t ** p["HBST"] + p["KBST"] ** p["HBST"])
    tauB = np.where(BSTss >= BST, p["TAUBSTU"], p["TAUBSTD"])
    d[36] = (BSTss - BST) / tauB

    fend = E2t ** p["HEND"] / (E2t ** p["HEND"] + p["KENDE"] ** p["HEND"])
    d[37] = p["KEND"] * fend * (1.0 - ENDO) - p["KENDO"] * ENDO

    # ================= 9. Bone mass ==========================================
    fb_E2 = E2 ** p["HBM"] / (E2 ** p["HBM"] + p["KBM"] ** p["HBM"])
    fb_NR = E2NORM ** p["HBM"] / (E2NORM ** p["HBM"] + p["KBM"] ** p["HBM"])
    accr = np.clip((20.0 - CA) / 10.0, 0.0, 1.0)
    # calcium/vitamin D only matters while the child is oestrogen-deprived and
    # the accrual window is still open; without the (supp * accr) gate the term
    # would keep adding indefinitely and lift ADULT bone mass, which it does not
    d[39] = (p["KZ"] * (fb_E2 - fb_NR) + p["KCATCH"] * (-BMDZ) * accr
             + p["KCAVD"] * p["CAVD"] * (fa / (fa + 0.30)) * accr) / 365.25

    # ================= 10. Body composition, symptoms, QoL ===================
    supp = fa / (fa + 0.30)
    d[38] = (p["BMIZTGT"] + p["KBMIG"] * supp - BMIZ) / p["TAUBMI"]

    d[41] = (E2 - E2TRK) / p["TAUTRK"]
    HFss = p["KHF"] * np.clip((E2TRK - E2) / p["HFSCALE"], 0.0, 1.0)
    d[40] = (HFss - HF) / p["TAUHF"]

    QOLss = (1.0 + 2.6 * np.clip((BST - 1.0) / 4.0, 0.0, 1.0)
             + 1.3 * np.clip((BA - CA) / 2.0, 0.0, 1.5)
             + 0.9 * HF / 6.0
             + 1.4 * np.clip((ENDO - 0.60) / 0.30, 0.0, 1.0) * (1.0 if CA < 10.0 else 0.0))
    d[42] = (QOLss - QOL) / p["TAUQOL"]
    d[43] = E2 / 365.25

    if not want_aux:
        return d, None
    return d, dict(GV=GV, S=S, fa=fa, fx=fx, dBA_yr=dBA_yr,
                   CLEUP=CLEUP, CTRIP=CTRIP, CHIST=CHIST, CNAF=CNAF,
                   CANT=CANT, AIEFF=AIEFF, RSo=RS, E2ss=E2ss)


AUXKEYS = ("GV", "S", "fa", "fx", "dBA_yr", "CLEUP", "CTRIP", "CHIST",
           "CNAF", "CANT", "AIEFF")
DEFREC = ("HT", "BA", "E2", "LH", "FSH", "BST", "UTV", "GPRES", "BMDZ",
          "IGF1", "ENDO", "HF", "QOL", "DHEAS", "TESTO", "GH", "FOL",
          "KND", "RS", "CUME2", "BMIZ", "INHB", "PULS", "AUTON")
MINREC = ("HT", "BA", "E2", "LH", "BST", "UTV", "GPRES", "ENDO", "BMDZ")


# -----------------------------------------------------------------------------
def init_state(n, p, ba0=5.0, ht0=108.0, mk0=1.0):
    y = np.zeros((NST, n))
    one = np.ones(n)
    y[17] = 1.0                                     # RS
    y[19] = mk0 * one                               # MKRN3
    y[20] = 0.05
    y[21] = 6.0
    y[22] = 0.10
    y[23] = 1.20
    y[24] = 8.0
    y[25] = 0.12
    y[26] = 4.0
    y[27] = 5.0
    y[28] = 25.0
    y[30] = 1.05
    y[31] = 130.0
    y[32] = ba0 * one
    y[33] = np.clip((p["BAFUS"] - ba0) / (p["BAFUS"] - p["BAREF"]), 0.0, 1.0) * one
    y[34] = ht0 * one
    y[35] = 1.3
    y[36] = 1.0
    y[37] = 0.05
    y[38] = 0.30
    y[41] = 4.0
    y[42] = 1.0
    return y


# -----------------------------------------------------------------------------
#  REGIMEN BUILDERS  (start/stop/dose may be scalars or per-subject arrays)
# -----------------------------------------------------------------------------
def reg_leup(dose_mg, interval_d, start_d, stop_d):
    return dict(kind="depot", amt=np.asarray(dose_mg, float) * 1000.0,
                interval=interval_d, start=start_d, stop=stop_d,
                bi=0, si=1, fburst=P["FBURSTL"])


def reg_trip(dose_mg, interval_d, start_d, stop_d):
    return dict(kind="depot", amt=np.asarray(dose_mg, float) * 1000.0,
                interval=interval_d, start=start_d, stop=stop_d,
                bi=4, si=5, fburst=P["FBURSTT"])


def reg_hist(start_d, stop_d, interval_d=365.0, dose_ug=50000.0):
    return dict(kind="implant", amt=dose_ug, interval=interval_d,
                start=start_d, stop=stop_d, si=7)


def reg_naf(daily_ug, start_d, stop_d, adherence=1.0, seed=1):
    return dict(kind="nasal", amt=np.asarray(daily_ug, float) / 3.0,
                interval=1.0 / 3.0, start=start_d, stop=stop_d, si=9,
                adherence=adherence, seed=seed)


def reg_antag(dose_ug, interval_d, start_d, stop_d):
    return dict(kind="simple", amt=dose_ug, interval=interval_d,
                start=start_d, stop=stop_d, si=11)


def reg_ai(dose_mg, start_d, stop_d):
    return dict(kind="simple", amt=np.asarray(dose_mg, float) * 1000.0 * P["FAI"],
                interval=1.0, start=start_d, stop=stop_d, si=13)


def reg_gh(dose_mg_day, start_d, stop_d):
    return dict(kind="simple", amt=np.asarray(dose_mg_day, float) * 1000.0,
                interval=1.0, start=start_d, stop=stop_d, si=14)


def reg_tam(dose_mg, start_d, stop_d):
    return dict(kind="simple", amt=np.asarray(dose_mg, float) * 1000.0 * 0.4,
                interval=1.0, start=start_d, stop=stop_d, si=16)


# -----------------------------------------------------------------------------
#  SIMULATOR  (RK4, fixed step, vectorised over subjects)
# -----------------------------------------------------------------------------
def simulate(n=1, years=16.0, dt=0.25, p=None, regimens=(), ba0=5.0, ht0=108.0,
             mk0=1.0, ca0=5.0, stride=4, record=DEFREC, extra=AUXKEYS):
    p = dict(P) if p is None else p
    nstep = int(round(years * 365.25 / dt))
    y = init_state(n, p, ba0=ba0, ht0=ht0, mk0=mk0)
    one = np.ones(n)

    regs = []
    for r in regimens:
        r = dict(r)
        r["next"] = np.asarray(r["start"], float) * one
        r["start_a"] = np.asarray(r["start"], float) * one
        r["stop_a"] = np.asarray(r["stop"], float) * one
        r["amt_a"] = np.asarray(r["amt"], float) * one
        r["int_a"] = np.asarray(r["interval"], float) * one
        if r["kind"] == "nasal":
            r["rng"] = np.random.default_rng(r.get("seed", 1))
            r["adh_a"] = np.asarray(r.get("adherence", 1.0), float) * one
        regs.append(r)

    nrec = nstep // stride + 1
    tt = np.zeros(nrec)
    out = {k: np.zeros((nrec, n)) for k in record}
    outx = {k: np.zeros((nrec, n)) for k in extra}

    t = 0.0
    j = 0
    for i in range(nstep + 1):
        for r in regs:
            due = ((t >= r["next"] - 1e-9) & (t < r["stop_a"])
                   & (t >= r["start_a"] - 1e-9))
            if np.any(due):
                if r["kind"] == "depot":
                    y[r["bi"]] += np.where(due, r["amt_a"] * r["fburst"], 0.0)
                    y[r["si"]] += np.where(due, r["amt_a"] * (1 - r["fburst"]), 0.0)
                elif r["kind"] == "implant":
                    y[7] = np.where(due, r["amt_a"], y[7])
                elif r["kind"] == "nasal":
                    take = r["rng"].random(n) < r["adh_a"]
                    y[r["si"]] += np.where(due & take, r["amt_a"], 0.0)
                else:
                    y[r["si"]] += np.where(due, r["amt_a"], 0.0)
                r["next"] = np.where(due, r["next"] + r["int_a"], r["next"])

        need = (i % stride == 0)
        k1, aux = rhs(t, y, p, ca0, want_aux=need)
        if need and j < nrec:
            tt[j] = t
            for k in record:
                out[k][j] = y[IX[k]]
            for k in extra:
                outx[k][j] = aux[k] * one
            j += 1
        if i == nstep:
            break
        k2, _ = rhs(t + dt / 2, y + dt / 2 * k1, p, ca0)
        k3, _ = rhs(t + dt / 2, y + dt / 2 * k2, p, ca0)
        k4, _ = rhs(t + dt, y + dt * k3, p, ca0)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        y[33] = np.clip(y[33], 0.0, 1.0)
        y[17] = np.clip(y[17], 0.0, 1.0)
        y[18] = np.clip(y[18], 0.0, 1.0)
        for k in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16):
            y[k] = np.maximum(y[k], 0.0)
        t = (i + 1) * dt

    return dict(t=tt, ca=ca0 + tt / 365.25, n=n, **out, **outx)


# -----------------------------------------------------------------------------
#  READ-OUTS  (all accept a subject index j)
# -----------------------------------------------------------------------------
def adult_height(r, j=0):
    return r["HT"][-1, j]


def _first_age(r, mask, j):
    idx = np.where(mask)[0]
    return r["ca"][idx[0]] if len(idx) else np.nan


def thelarche_age(r, j=0):
    return _first_age(r, r["BST"][:, j] >= 2.0, j)


def menarche_age(r, j=0):
    ok = ((r["ENDO"][:, j] > 0.60) & (r["UTV"][:, j] > 4.0)
          & (r["BA"][:, j] > 10.0))
    return _first_age(r, ok, j)


def fusion_age(r, j=0):
    return _first_age(r, r["GPRES"][:, j] <= 1e-6, j)


def at_age(r, ca, key, j=0):
    return float(np.interp(ca, r["ca"], r[key][:, j]))


def pah_bp(r, i, j=0, male=False):
    return r["HT"][i, j] / bp_fraction(r["BA"][i, j], male)


def ba_ca_ratio(r, lo, hi, j=0):
    return (at_age(r, hi, "BA", j) - at_age(r, lo, "BA", j)) / (hi - lo)


def gv_mean(r, lo, hi, j=0):
    m = (r["ca"] >= lo) & (r["ca"] <= hi)
    return float(r["GV"][m, j].mean()) if m.any() else np.nan


def run(mk0=1.0, p=None, regs=(), years=16.0, ba0=5.0, ht0=108.0, n=1,
        dt=0.25, stride=4, record=DEFREC, extra=AUXKEYS):
    return simulate(n=n, years=years, dt=dt, p=dict(P) if p is None else p,
                    regimens=regs, ba0=ba0, ht0=ht0, mk0=mk0, stride=stride,
                    record=record, extra=extra)


# =============================================================================
#  REPORT
# =============================================================================
OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def sec(title):
    say()
    say("=" * 84)
    say(title)
    say("=" * 84)


MK_CPP = 0.400        # MKRN3 loss-of-function phenotype used throughout
MK_NORM = 1.00


def main():
    pN = dict(P)

    sec("0.  MODEL STRUCTURE")
    say(f"  ODEs                      : {NST}")
    say(f"  named parameters          : {len(P)}")
    say(f"  integrator                : fixed-step RK4, dt = 0.25 d")

    # =====================================================================
    sec("1.  NATURAL HISTORY — normal puberty vs untreated CPP (girls)")
    resN = run(MK_NORM, pN, years=16.0)
    resC = run(MK_CPP, pN, years=16.0)
    for lab, r in (("NORMAL girl (MKRN3 intact)", resN),
                   ("CPP girl, untreated (MKRN3 LoF)", resC)):
        say()
        say(f"--- {lab} ---")
        say(f"  thelarche (Tanner B2)        {thelarche_age(r):6.2f} y")
        say(f"  menarche                     {menarche_age(r):6.2f} y")
        say(f"  epiphyseal fusion            {fusion_age(r):6.2f} y")
        say(f"  ADULT HEIGHT                 {adult_height(r):6.1f} cm")
        i = int(r["GV"][:, 0].argmax())
        say(f"  peak growth velocity         {r['GV'][i,0]:6.2f} cm/yr at CA {r['ca'][i]:.2f} y")
        say(f"  peak E2                      {r['E2'].max():6.1f} pg/mL")
        say(f"  LH / E2 / UTV at CA 8.0      {at_age(r,8.0,'LH'):.2f} IU/L / "
            f"{at_age(r,8.0,'E2'):.1f} pg/mL / {at_age(r,8.0,'UTV'):.2f} mL")
        say(f"  bone age at CA 8.0           {at_age(r,8.0,'BA'):6.2f} y "
            f"(BA-CA {at_age(r,8.0,'BA')-8.0:+.2f})")
        say(f"  dBA/dCA over CA 7.0-9.0      {ba_ca_ratio(r,7.0,9.0):6.2f}")
        say(f"  IGF-1 at CA 8.0              {at_age(r,8.0,'IGF1'):6.0f} ng/mL")
        say(f"  cumulative E2 exposure       {r['CUME2'][-1,0]:6.0f} pg/mL*yr")
        say(f"  BMD Z at CA 8 / adult        {at_age(r,8.0,'BMDZ'):+.2f} / {r['BMDZ'][-1,0]:+.2f}")
        say(f"  peak psychosocial index      {r['QOL'][:,0].max():6.2f} / 10")
    say()
    say(f"  >>> UNTREATED HEIGHT LOSS = {adult_height(resN)-adult_height(resC):.1f} cm "
        f"({adult_height(resC):.1f} vs {adult_height(resN):.1f} cm)")

    # =====================================================================
    sec("2.  THE SIGN FLIP — height gain vs BONE AGE at start of therapy")
    say("A GnRH agonist suppresses BOTH arms of the oestradiol signal at once.")
    say("Nothing in the model encodes a 'treat before age 8' rule: the rule is")
    say("the OUTPUT of the competition between the (+) growth-velocity arm and")
    say("the (-) plate-reserve arm.")
    starts = np.array([6.6, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0])
    ns = len(starts)
    sd = (starts - 5.0) * 365.25
    rU = run(MK_CPP, pN, years=16.5, n=ns)
    rT = run(MK_CPP, pN, years=16.5, n=ns,
             regs=[reg_leup(11.25, 84.0, sd, sd + 365.25 * 6.0)])
    say()
    say(f"{'start CA':>9} {'start BA':>9} {'GPRES':>7} {'noRx cm':>9} {'Rx cm':>8} "
        f"{'GAIN cm':>9} {'GV on Rx':>9} {'dBA/dCA Rx':>11}")
    bas, gains = [], []
    for j in range(ns):
        ba_s = at_age(rU, starts[j], "BA", j)
        gp_s = at_age(rU, starts[j], "GPRES", j)
        g = adult_height(rT, j) - adult_height(rU, j)
        bas.append(ba_s); gains.append(g)
        hi = min(starts[j] + 2.0, 16.4)
        say(f"{starts[j]:9.1f} {ba_s:9.2f} {gp_s:7.3f} {adult_height(rU,j):9.1f} "
            f"{adult_height(rT,j):8.1f} {g:+9.2f} {gv_mean(rT,starts[j]+0.3,hi,j):9.2f} "
            f"{ba_ca_ratio(rT,starts[j]+0.3,hi,j):11.2f}")
    bas = np.array(bas); gains = np.array(gains)
    say()
    for thr in (4.0, 2.0, 1.0, 0.0):
        if gains.min() < thr < gains.max():
            say(f"  gain falls below {thr:>4.1f} cm at bone age "
                f"{np.interp(-thr, -gains, bas):.2f} y "
                f"(chronological age {np.interp(-thr,-gains,starts):.2f} y)")

    # =====================================================================
    sec("3.  GROWTH VELOCITY IS A MISLEADING MONITOR")
    s0 = (7.5 - 5.0) * 365.25
    rT = run(MK_CPP, pN, years=16.5, stride=2,
             regs=[reg_leup(11.25, 84.0, s0, s0 + 365.25 * 5.0)])
    rU = run(MK_CPP, pN, years=16.5, stride=2)
    say("Treated child (leuprolide 11.25 mg q12wk from CA 7.5):")
    say(f"{'CA':>6} {'GV':>7} {'dBA/dCA':>8} {'BA-CA':>7} {'E2':>6} {'LH':>6} "
        f"{'UTV':>6} {'BP-PAH':>7} {'true-final':>10}")
    for ca in (7.4, 7.6, 8.0, 8.5, 9.0, 10.0, 11.0, 12.0, 12.5):
        i = int(np.argmin(np.abs(rT["ca"] - ca)))
        say(f"{ca:6.1f} {rT['GV'][i,0]:7.2f} {rT['dBA_yr'][i,0]:8.2f} "
            f"{rT['BA'][i,0]-ca:+7.2f} {rT['E2'][i,0]:6.1f} {rT['LH'][i,0]:6.2f} "
            f"{rT['UTV'][i,0]:6.2f} {pah_bp(rT,i):7.1f} {adult_height(rT):10.1f}")
    say()
    say("Untreated, same ages:")
    say(f"{'CA':>6} {'GV':>7} {'dBA/dCA':>8} {'BP-PAH':>7}")
    for ca in (7.6, 8.5, 9.0, 10.0):
        i = int(np.argmin(np.abs(rU["ca"] - ca)))
        say(f"{ca:6.1f} {rU['GV'][i,0]:7.2f} {rU['dBA_yr'][i,0]:8.2f} {pah_bp(rU,i):7.1f}")

    sec("3b. ACROSS SUBJECTS: the model REFUTES the stronger form of this claim")
    say("Hypothesis under test: if growth velocity is a misleading monitor, then")
    say("across subjects a LOWER on-therapy velocity should go with a LARGER")
    say("adult-height gain.  The model says NO, and the reason is instructive:")
    say("between subjects, both quantities are driven by the SAME third variable")
    say("(remaining growth-plate reserve at the start of therapy), so they move")
    say("TOGETHER.  The monitoring trap is therefore purely WITHIN-subject -- the")
    say("step change at the moment therapy starts -- and not a between-subject")
    say("ranking.  This is reported rather than tuned away.")
    rng = np.random.default_rng(7)
    NS = 200
    mk = rng.uniform(0.28, 0.60, NS)
    rU2 = run(mk, pN, years=16.5, n=NS, record=MINREC, extra=("GV", "dBA_yr"))
    t0 = np.array([thelarche_age(rU2, j) for j in range(NS)])
    keep = np.isfinite(t0) & (t0 < 8.0)
    sdp = (np.clip(t0 + 0.4, 6.4, 11.0) - 5.0) * 365.25
    rT2 = run(mk, pN, years=16.5, n=NS, record=MINREC, extra=("GV", "dBA_yr"),
              regs=[reg_leup(11.25, 84.0, sdp, sdp + 365.25 * 5.5)])
    gvs, dbas, gains2 = [], [], []
    for j in np.where(keep)[0]:
        cs = 5.0 + sdp[j] / 365.25
        gvs.append(gv_mean(rT2, cs + 0.5, cs + 2.0, j))
        dbas.append(ba_ca_ratio(rT2, cs + 0.5, cs + 2.0, j))
        gains2.append(adult_height(rT2, j) - adult_height(rU2, j))
    gvs, dbas, gains2 = map(np.array, (gvs, dbas, gains2))
    say(f"n = {len(gvs)} virtual subjects with thelarche < 8.0 y")
    say(f"  corr( mean on-therapy GV  , adult-height GAIN )  r = {np.corrcoef(gvs,gains2)[0,1]:+.3f}")
    say(f"  corr( on-therapy dBA/dCA  , adult-height GAIN )  r = {np.corrcoef(dbas,gains2)[0,1]:+.3f}")
    say(f"  on-therapy GV      mean {gvs.mean():.2f}  range {gvs.min():.2f}-{gvs.max():.2f} cm/yr")
    say(f"  on-therapy dBA/dCA mean {dbas.mean():.2f}  range {dbas.min():.2f}-{dbas.max():.2f}")
    say(f"  adult-height gain  mean {gains2.mean():+.2f}  range {gains2.min():+.2f} to {gains2.max():+.2f} cm")

    # =====================================================================
    sec("4.  THE FLARE — what does it actually cost?")
    s0 = (7.5 - 5.0) * 365.25
    e0 = s0 + 365.25 * 5.0
    rAg = run(MK_CPP, pN, years=16.5, stride=1, regs=[reg_leup(3.75, 28.0, s0, e0)])
    rAn = run(MK_CPP, pN, years=16.5, stride=1, regs=[reg_antag(18000.0, 28.0, s0, e0)])
    rU = run(MK_CPP, pN, years=16.5, stride=1)
    say(f"{'':34} {'agonist':>10} {'antagonist':>12} {'no therapy':>11}")
    m14 = lambda r: (r["t"] >= s0) & (r["t"] <= s0 + 14)
    m30 = lambda r: (r["t"] >= s0) & (r["t"] <= s0 + 30)
    say(f"{'peak LH, first 14 d (IU/L)':34} {rAg['LH'][m14(rAg),0].max():10.2f} "
        f"{rAn['LH'][m14(rAn),0].max():12.2f} {rU['LH'][m14(rU),0].max():11.2f}")
    say(f"{'peak E2, first 14 d (pg/mL)':34} {rAg['E2'][m14(rAg),0].max():10.2f} "
        f"{rAn['E2'][m14(rAn),0].max():12.2f} {rU['E2'][m14(rU),0].max():11.2f}")
    for dd in (7, 14, 28, 56):
        say(f"{'LH at day '+str(dd)+' (IU/L)':34} "
            f"{np.interp(s0+dd, rAg['t'], rAg['LH'][:,0]):10.3f} "
            f"{np.interp(s0+dd, rAn['t'], rAn['LH'][:,0]):12.3f} "
            f"{np.interp(s0+dd, rU['t'], rU['LH'][:,0]):11.3f}")
    say(f"{'peak endometrial state, 30 d':34} {rAg['ENDO'][m30(rAg),0].max():10.3f} "
        f"{rAn['ENDO'][m30(rAn),0].max():12.3f} {rU['ENDO'][m30(rU),0].max():11.3f}")
    say(f"{'bone age at CA 9.0 (y)':34} {at_age(rAg,9.0,'BA'):10.2f} "
        f"{at_age(rAn,9.0,'BA'):12.2f} {at_age(rU,9.0,'BA'):11.2f}")
    say(f"{'ADULT HEIGHT (cm)':34} {adult_height(rAg):10.1f} "
        f"{adult_height(rAn):12.1f} {adult_height(rU):11.1f}")
    say(f"{'gain vs untreated (cm)':34} {adult_height(rAg)-adult_height(rU):+10.2f} "
        f"{adult_height(rAn)-adult_height(rU):+12.2f} {0.0:11.2f}")
    say()
    say(f"  >>> COST OF THE FLARE = {adult_height(rAn)-adult_height(rAg):.2f} cm of adult height")

    NS = 250
    rng = np.random.default_rng(11)
    mk = rng.uniform(0.28, 0.60, NS)
    rU3 = run(mk, pN, years=12.0, n=NS, record=MINREC, extra=())
    t0 = np.array([thelarche_age(rU3, j) for j in range(NS)])
    keep = np.isfinite(t0) & (t0 < 8.2)
    sdp = (np.clip(t0 + 0.6, 6.3, 11.0) - 5.0) * 365.25
    rB = run(mk, pN, years=12.0, n=NS, stride=2, record=MINREC, extra=(),
             regs=[reg_leup(3.75, 28.0, sdp, sdp + 365.25 * 3)])
    bleeds = 0
    used = 0
    pks = []
    for j in np.where(keep)[0]:
        m0 = (rB["t"] >= sdp[j]) & (rB["t"] <= sdp[j] + 18)
        m1 = (rB["t"] > sdp[j] + 18) & (rB["t"] <= sdp[j] + 60)
        used += 1
        pk = rB["ENDO"][m0, j].max()
        pks.append(pk)
        if pk > P["ENDOBLEED"] and (pk - rB["ENDO"][m1, j].min()) / pk > 0.40:
            bleeds += 1
    pks = np.array(pks)
    say(f"  endometrial priming at start of therapy: mean {pks.mean():.3f}, "
        f"range {pks.min():.3f}-{pks.max():.3f} (bleed threshold "
        f"ENDOBLEED = {P['ENDOBLEED']:.2f})")
    say(f"  withdrawal bleeding after 1st depot: {bleeds}/{used} = "
        f"{100.0*bleeds/max(used,1):.1f}%   (reported 5-10%)")
    rBx = run(mk, pN, years=12.0, n=NS, stride=2, record=MINREC, extra=(),
              regs=[reg_antag(18000.0, 28.0, sdp, sdp + 365.25 * 3)])
    bl2 = 0
    for j in np.where(keep)[0]:
        m0 = (rBx["t"] >= sdp[j]) & (rBx["t"] <= sdp[j] + 18)
        m1 = (rBx["t"] > sdp[j] + 18) & (rBx["t"] <= sdp[j] + 60)
        pk = rBx["ENDO"][m0, j].max()
        if pk > P["ENDOBLEED"] and (pk - rBx["ENDO"][m1, j].min()) / pk > 0.40:
            bl2 += 1
    say(f"  SAME cohort on a GnRH ANTAGONIST (no flare): {bl2}/{used} = "
        f"{100.0*bl2/max(used,1):.1f}%  -> in this structure the bleed is")
    say(f"  OESTROGEN WITHDRAWAL from a primed endometrium, not the flare itself.")

    # =====================================================================
    sec("5.  FORMULATION — trough coverage, not potency, is the design variable")
    s0 = (7.4 - 5.0) * 365.25
    e0 = s0 + 365.25 * 5.0
    forms = [
        ("leuprolide 3.75 mg IM q28d", [reg_leup(3.75, 28.0, s0, e0)]),
        ("leuprolide 7.5 mg IM q28d", [reg_leup(7.5, 28.0, s0, e0)]),
        ("leuprolide 11.25 mg IM q12wk", [reg_leup(11.25, 84.0, s0, e0)]),
        ("leuprolide 30 mg IM q12wk", [reg_leup(30.0, 84.0, s0, e0)]),
        ("triptorelin 11.25 mg IM q12wk", [reg_trip(11.25, 84.0, s0, e0)]),
        ("triptorelin 22.5 mg IM q24wk", [reg_trip(22.5, 168.0, s0, e0)]),
        ("histrelin 50 mg implant q12mo", [reg_hist(s0, e0)]),
        ("nafarelin nasal 1800 ug/d adh1.00", [reg_naf(1800.0, s0, e0, 1.00, 3)]),
        ("nafarelin nasal 1800 ug/d adh0.80", [reg_naf(1800.0, s0, e0, 0.80, 3)]),
        ("leuprolide 3.75 mg given LATE q42d", [reg_leup(3.75, 42.0, s0, e0)]),
        ("leuprolide 3.75 mg given LATE q56d", [reg_leup(3.75, 56.0, s0, e0)]),
        ("leuprolide 1.875 mg q28d underdosed", [reg_leup(1.875, 28.0, s0, e0)]),
        ("GnRH antagonist 18 mg q28d", [reg_antag(18000.0, 28.0, s0, e0)]),
    ]
    rU = run(MK_CPP, pN, years=16.5)
    say(f"{'formulation':36} {'meanCp':>7} {'mean fa':>8} {'%t LH>0.5':>10} "
        f"{'%t E2>10':>9} {'AdultHt':>8} {'gain':>7}")
    for lab, regs in forms:
        r = run(MK_CPP, pN, years=16.5, regs=regs, stride=2)
        m = (r["t"] >= s0 + 90) & (r["t"] <= e0)
        cp = (r["CLEUP"] + r["CTRIP"] + r["CHIST"] + r["CNAF"])[m, 0]
        say(f"{lab:36} {cp.mean():7.3f} {r['fa'][m,0].mean():8.3f} "
            f"{100*np.mean(r['LH'][m,0]>0.5):10.1f} {100*np.mean(r['E2'][m,0]>10):9.1f} "
            f"{adult_height(r):8.1f} {adult_height(r)-adult_height(rU):+7.2f}")
    say()
    say("Within-interval profile, leuprolide 3.75 mg q28d, 12th cycle:")
    r = run(MK_CPP, pN, years=16.5, regs=[reg_leup(3.75, 28.0, s0, e0)], stride=1)
    base = s0 + 28 * 11
    say(f"{'day in cycle':>13} {'Cp ng/mL':>9} {'fa':>7} {'RS':>7} {'LH':>7} {'E2':>7}")
    for dd in (0.5, 2, 5, 10, 14, 20, 24, 27.5):
        i = int(np.argmin(np.abs(r["t"] - (base + dd))))
        say(f"{dd:13.1f} {r['CLEUP'][i,0]:9.3f} {r['fa'][i,0]:7.3f} "
            f"{r['RS'][i,0]:7.3f} {r['LH'][i,0]:7.3f} {r['E2'][i,0]:7.2f}")

    # =====================================================================
    sec("6.  ADHERENCE SWEEP (nafarelin nasal) — the hidden failure mode")
    say(f"{'adherence':>10} {'mean fa':>8} {'%t LH>0.5':>10} {'%t E2>10':>9} "
        f"{'AdultHt':>8} {'gain':>7}")
    for adh in (1.00, 0.90, 0.80, 0.70, 0.60, 0.50):
        r = run(MK_CPP, pN, years=16.5, stride=2,
                regs=[reg_naf(1800.0, s0, e0, adh, 5)])
        m = (r["t"] >= s0 + 90) & (r["t"] <= e0)
        say(f"{adh:10.2f} {r['fa'][m,0].mean():8.3f} {100*np.mean(r['LH'][m,0]>0.5):10.1f} "
            f"{100*np.mean(r['E2'][m,0]>10):9.1f} {adult_height(r):8.1f} "
            f"{adult_height(r)-adult_height(rU):+7.2f}")

    # =====================================================================
    sec("7.  GnRH-INDEPENDENT DISEASE (McCune-Albright) — mechanism mismatch")
    pM = dict(pN); pM["AUTSET"] = 0.42; pM["KSEC"] = 0.22
    pAI = dict(pM); pAI["IMAXAI"] = 0.992; pAI["IC50AI"] = 2.2
    sM = (6.0 - 5.0) * 365.25
    eM = sM + 365.25 * 7.0
    rM0 = run(MK_NORM, pM, years=16.5)
    say(f"{'strategy':32} {'meanE2 7-11y':>13} {'BA@CA10':>8} {'AdultHt':>8} {'gain':>7}")
    for lab, regs, pp in [
            ("no treatment", [], pM),
            ("leuprolide 11.25 q12wk", [reg_leup(11.25, 84.0, sM, eM)], pM),
            ("potent AI (letrozole-like)", [reg_ai(2.5, sM, eM)], pAI),
            ("AI + leuprolide", [reg_ai(2.5, sM, eM),
                                 reg_leup(11.25, 84.0, sM, eM)], pAI),
            ("AI + tamoxifen", [reg_ai(2.5, sM, eM), reg_tam(20.0, sM, eM)], pAI)]:
        r = run(MK_NORM, pp, years=16.5, regs=regs)
        m = (r["ca"] >= 7.0) & (r["ca"] <= 11.0)
        say(f"{lab:32} {r['E2'][m,0].mean():13.1f} {at_age(r,10.0,'BA'):8.2f} "
            f"{adult_height(r):8.1f} {adult_height(r)-adult_height(rM0):+7.2f}")
    rC = run(MK_CPP, pN, years=16.5, regs=[reg_leup(11.25, 84.0, sM, eM)])
    rC0 = run(MK_CPP, pN, years=16.5)
    say()
    say(f"  same agonist regimen in TRUE CENTRAL disease: gain "
        f"{adult_height(rC)-adult_height(rC0):+.2f} cm")

    # =====================================================================
    sec("8.  BOYS — aromatase inhibition separates the two oestradiol signs")
    pB = dict(pN); pB["SEXM"] = 1.0; pB["BAFUS"] = 16.8
    rBn = run(MK_NORM, pB, years=18.0, ht0=110.0)
    say(f"  NORMAL boy reference: fusion {fusion_age(rBn):.2f} y, adult height "
        f"{adult_height(rBn):.1f} cm, peak GV {rBn['GV'].max():.2f} cm/yr at CA "
        f"{rBn['ca'][rBn['GV'][:,0].argmax()]:.2f} y")
    sB = (8.0 - 5.0) * 365.25
    eB = sB + 365.25 * 5.0
    rB0 = run(MK_CPP, pB, years=17.5, ht0=108.5)
    say(f"{'strategy':30} {'meanE2':>7} {'meanT':>8} {'BA@CA11':>8} {'AdultHt':>8} "
        f"{'gain':>7} {'BMDZ nadir':>11} {'BMDZ adult':>11}")
    for lab, regs in [
            ("untreated CPP boy", []),
            ("leuprolide 11.25 q12wk", [reg_leup(11.25, 84.0, sB, eB)]),
            ("anastrozole-like AI alone", [reg_ai(1.0, sB, eB)]),
            ("GnRHa + AI", [reg_leup(11.25, 84.0, sB, eB), reg_ai(1.0, sB, eB)])]:
        r = run(MK_CPP, pB, years=17.5, ht0=108.5, regs=regs)
        m = (r["ca"] >= 8.5) & (r["ca"] <= 12.0)
        say(f"{lab:30} {r['E2'][m,0].mean():7.2f} {r['TESTO'][m,0].mean():8.1f} "
            f"{at_age(r,11.0,'BA'):8.2f} {adult_height(r):8.1f} "
            f"{adult_height(r)-adult_height(rB0):+7.2f} {r['BMDZ'][:,0].min():+11.2f} "
            f"{r['BMDZ'][-1,0]:+11.2f}")

    # =====================================================================
    sec("9.  rhGH ADD-ON — and why the model says it is partly self-defeating")
    rU = run(MK_CPP, pN, years=16.5)
    for lab_start, ca_s, dur, wt in (("LATE start (CA 8.6)", 8.6, 4.5, 32.0),
                                     ("EARLY start (CA 7.2)", 7.2, 5.0, 30.0)):
        s = (ca_s - 5.0) * 365.25
        e = s + 365.25 * dur
        say()
        say(f"--- {lab_start} ---")
        say(f"{'strategy':30} {'meanIGF1':>9} {'meanGV':>7} {'dBA/dCA':>8} "
            f"{'AdultHt':>8} {'gain':>7}")
        for lab, regs in (("GnRHa alone", [reg_leup(11.25, 84.0, s, e)]),
                          ("GnRHa + rhGH 0.043 mg/kg/d",
                           [reg_leup(11.25, 84.0, s, e), reg_gh(0.043 * wt, s, e)]),
                          ("rhGH alone", [reg_gh(0.043 * wt, s, e)])):
            r = run(MK_CPP, pN, years=16.5, regs=regs)
            lo, hi = ca_s + 0.3, min(ca_s + 3.0, 16.4)
            m = (r["ca"] >= lo) & (r["ca"] <= hi)
            say(f"{lab:30} {r['IGF1'][m,0].mean():9.0f} {r['GV'][m,0].mean():7.2f} "
                f"{ba_ca_ratio(r,lo,hi):8.2f} {adult_height(r):8.1f} "
                f"{adult_height(r)-adult_height(rU):+7.2f}")

    # =====================================================================
    sec("10. SLOWLY-PROGRESSIVE VARIANT — the over-treatment problem, quantified")
    mks = np.array([0.34, 0.40, 0.46, 0.52, 0.58, 0.66, 0.76])
    rUv = run(mks, pN, years=17.0, n=len(mks))
    t0 = np.array([thelarche_age(rUv, j) for j in range(len(mks))])
    sdv = (np.clip(np.nan_to_num(t0, nan=12.0) + 0.5, 6.3, 11.5) - 5.0) * 365.25
    rTv = run(mks, pN, years=17.0, n=len(mks),
              regs=[reg_leup(11.25, 84.0, sdv, sdv + 365.25 * 5)])
    say(f"{'MK0':>6} {'thelarche':>10} {'peakE2':>7} {'menarche':>9} "
        f"{'noRx cm':>8} {'Rx cm':>7} {'GAIN':>7}")
    for j in range(len(mks)):
        th = t0[j]
        say(f"{mks[j]:6.2f} {th if np.isfinite(th) else float('nan'):10.2f} "
            f"{rUv['E2'][:,j].max():7.1f} {menarche_age(rUv,j):9.2f} "
            f"{adult_height(rUv,j):8.1f} {adult_height(rTv,j):7.1f} "
            f"{adult_height(rTv,j)-adult_height(rUv,j):+7.2f}")

    # =====================================================================
    sec("11. BONE MINERAL DENSITY — the dip is regression of an ADVANCED baseline")
    s = (7.4 - 5.0) * 365.25
    e = s + 365.25 * 4.0
    rT = run(MK_CPP, pN, years=16.5, regs=[reg_leup(11.25, 84.0, s, e)])
    rTC = run(MK_CPP, dict(pN, CAVD=1.0), years=16.5,
              regs=[reg_leup(11.25, 84.0, s, e)])
    rU = run(MK_CPP, pN, years=16.5)
    rN = run(MK_NORM, pN, years=16.5)
    say(f"{'CA (y)':>7} {'normal':>8} {'CPP noRx':>9} {'CPP+GnRHa':>10} {'+Ca/VitD':>10}")
    for ca in (7.0, 7.4, 8.0, 9.0, 10.0, 11.4, 12.0, 14.0, 16.0, 20.4):
        say(f"{ca:7.1f} {at_age(rN,ca,'BMDZ'):+8.2f} {at_age(rU,ca,'BMDZ'):+9.2f} "
            f"{at_age(rT,ca,'BMDZ'):+10.2f} {at_age(rTC,ca,'BMDZ'):+10.2f}")
    i0 = int(np.argmin(np.abs(rT["ca"] - 7.4)))
    k = i0 + int(np.argmin(rT["BMDZ"][i0:, 0]))
    say(f"  nadir on therapy {rT['BMDZ'][k,0]:+.2f} at CA {rT['ca'][k]:.1f} y "
        f"(treatment CA 7.4-11.4)")
    say(f"  final (CA {rT['ca'][-1]:.1f}) : normal {rN['BMDZ'][-1,0]:+.2f} / "
        f"treated {rT['BMDZ'][-1,0]:+.2f} / untreated CPP {rU['BMDZ'][-1,0]:+.2f}")

    # =====================================================================
    sec("12. AXIS RECOVERY AND MENARCHE AFTER STOPPING")
    s = (7.4 - 5.0) * 365.25
    for dur in (3.0, 4.0, 5.0):
        e = s + 365.25 * dur
        r = run(MK_CPP, pN, years=17.0, stride=2, regs=[reg_leup(11.25, 84.0, s, e)])
        ca_stop = 5.0 + e / 365.25
        aft = r["t"] > e
        i1 = np.where(aft & (r["LH"][:, 0] > 1.0))[0]
        i2 = np.where(aft & (r["E2"][:, 0] > 20.0))[0]
        men = menarche_age(r)
        say(f"  duration {dur:.1f} y (stop CA {ca_stop:5.2f}): "
            f"LH>1 IU/L after {((r['t'][i1[0]]-e)/30.44 if len(i1) else float('nan')):5.1f} mo, "
            f"E2>20 after {((r['t'][i2[0]]-e)/30.44 if len(i2) else float('nan')):5.1f} mo, "
            f"menarche {men:5.2f} y = {12*(men-ca_stop):5.1f} mo after stop")
    say()
    say("Receptor-pool recovery after the last depot (11.25 mg q12wk x 4 y):")
    e = s + 365.25 * 4.0
    r = run(MK_CPP, pN, years=17.0, stride=1, regs=[reg_leup(11.25, 84.0, s, e)])
    say(f"{'days after stop':>16} {'Cp':>7} {'fa':>7} {'RS':>7} {'LH':>7} {'E2':>7}")
    for dd in (0, 30, 60, 90, 120, 180, 270):
        i = int(np.argmin(np.abs(r["t"] - (e + dd))))
        say(f"{dd:16d} {r['CLEUP'][i,0]:7.3f} {r['fa'][i,0]:7.3f} {r['RS'][i,0]:7.3f} "
            f"{r['LH'][i,0]:7.2f} {r['E2'][i,0]:7.1f}")

    # =====================================================================
    sec("13. VIRTUAL POPULATION (n=400) — treat-all vs treat-by-bone-age")
    rng = np.random.default_rng(2024)
    NSP = 400
    mk = rng.uniform(0.28, 0.95, NSP)
    ht0 = 108.0 + rng.normal(0, 3.2, NSP)
    pp = dict(pN)
    pp["KMK"] = P["KMK"] * np.exp(rng.normal(0, 0.15, NSP))
    pp["KE2"] = P["KE2"] * np.exp(rng.normal(0, 0.12, NSP))
    pp["GVBASE"] = P["GVBASE"] * np.exp(rng.normal(0, 0.05, NSP))
    rUp = run(mk, pp, years=17.0, n=NSP, ht0=ht0, record=MINREC, extra=())
    t0 = np.array([thelarche_age(rUp, j) for j in range(NSP)])
    sel = np.isfinite(t0) & (t0 < 8.0)
    # REFERRAL DELAY: families and referral pathways vary, so bone age at the
    # moment the decision is actually made spans a wide range.  Without this,
    # every simulated child would be diagnosed with a young bone age and the
    # policy comparison below would be degenerate.
    delay = rng.uniform(0.3, 3.8, NSP)
    ca_dx = np.clip(np.nan_to_num(t0, nan=11.5) + delay, 6.3, 11.8)
    sdp = (ca_dx - 5.0) * 365.25
    rTp = run(mk, pp, years=17.0, n=NSP, ht0=ht0, record=MINREC, extra=(),
              regs=[reg_leup(11.25, 84.0, sdp, sdp + 365.25 * 5)])
    idx = np.where(sel)[0]
    ba_dx = np.array([at_age(rUp, ca_dx[j], "BA", j) for j in idx])
    hu = np.array([adult_height(rUp, j) for j in idx])
    ht = np.array([adult_height(rTp, j) for j in idx])
    gn = ht - hu
    say(f"n screened {NSP}; presenting with thelarche < 8.0 y : {len(idx)}")
    say(f"  bone age at diagnosis   mean {ba_dx.mean():.2f} y  range "
        f"{ba_dx.min():.2f}-{ba_dx.max():.2f}")
    say(f"  untreated adult height  mean {hu.mean():6.1f} cm  sd {hu.std():4.1f}  "
        f"range {hu.min():.1f}-{hu.max():.1f}")
    say(f"  treated   adult height  mean {ht.mean():6.1f} cm  sd {ht.std():4.1f}")
    say(f"  TREAT-ALL mean gain     {gn.mean():+.2f} cm (sd {gn.std():.2f}, "
        f"range {gn.min():+.2f} to {gn.max():+.2f})")
    say(f"  fraction gain <1 cm {100*np.mean(gn<1.0):.1f}% ; gain <0 cm "
        f"{100*np.mean(gn<0.0):.1f}%")
    say()
    say(f"{'policy':28} {'% treated':>10} {'cohort gain':>12} {'gain if treated':>16}")
    for cut in (10.0, 10.5, 11.0, 11.5, 12.0, 99.0):
        s2 = ba_dx < cut
        if s2.sum() == 0:
            continue
        pol = np.where(s2, gn, 0.0)
        lab = "treat all" if cut > 90 else f"treat if BA < {cut:4.1f}"
        say(f"{lab:28} {100*s2.mean():10.1f} {pol.mean():+12.2f} {gn[s2].mean():+16.2f}")
    say()
    say(f"  corr( bone age at diagnosis , gain )  r = {np.corrcoef(ba_dx,gn)[0,1]:+.3f}")
    say(f"  corr( thelarche age , gain )          r = {np.corrcoef(t0[idx],gn)[0,1]:+.3f}")
    say(f"  corr( untreated adult height , gain ) r = {np.corrcoef(hu,gn)[0,1]:+.3f}")

    # =====================================================================
    sec("14. MONITORING — which surrogate separates adequate from inadequate?")
    s = (7.4 - 5.0) * 365.25
    e = s + 365.25 * 5.0
    rG = run(MK_CPP, pN, years=16.5, stride=2, regs=[reg_leup(11.25, 84.0, s, e)])
    rBd = run(MK_CPP, pN, years=16.5, stride=2, regs=[reg_naf(1800.0, s, e, 0.55, 9)])
    say(f"{'surrogate':26} {'adequate':>10} {'inadequate':>12} {'% apart':>9} {'usable':>7}")
    msk = lambda r: (r["t"] >= s + 120) & (r["t"] <= e)
    for lab, key in (("basal LH (IU/L)", "LH"), ("E2 (pg/mL)", "E2"),
                     ("uterine volume (mL)", "UTV"),
                     ("growth velocity (cm/yr)", "GV"), ("dBA/dCA", "dBA_yr"),
                     ("Tanner breast stage", "BST"), ("IGF-1 (ng/mL)", "IGF1")):
        a = rG[key][msk(rG), 0].mean()
        b = rBd[key][msk(rBd), 0].mean()
        rel = abs(b - a) / max(abs(a), 1e-9)
        say(f"{lab:26} {a:10.2f} {b:12.2f} {100*rel:9.0f} "
            f"{'YES' if rel>0.30 else 'no':>7}")
    say()
    say(f"  adult height: adequate {adult_height(rG):.1f} cm vs inadequate "
        f"{adult_height(rBd):.1f} cm (difference "
        f"{adult_height(rG)-adult_height(rBd):.1f} cm)")

    # =====================================================================
    sec("15. LOCAL SENSITIVITY OF THE ADULT-HEIGHT GAIN (+/-20%)")
    s = (7.5 - 5.0) * 365.25
    e = s + 365.25 * 5.0
    keys = ["M0", "ME", "KBA", "MI", "KIG", "EE2GV", "ETGV", "KE2", "KDES",
            "KREC", "AINT", "EC50L", "GVBASE", "PGP", "PIGF", "EGH", "KIGF",
            "MA", "XTRA", "BAFUS", "KFOL", "TAUKND"]
    nk = len(keys)
    pv = {k: np.repeat(P[k], 2 * nk) if not isinstance(P[k], str) else P[k] for k in P}
    pv = dict(P)
    for k in P:
        pv[k] = np.full(2 * nk, float(P[k]))
    for a, kk in enumerate(keys):
        pv[kk][2 * a] = P[kk] * 0.8
        pv[kk][2 * a + 1] = P[kk] * 1.2
    rUs = run(MK_CPP, pv, years=16.5, n=2 * nk, record=MINREC, extra=())
    rTs = run(MK_CPP, pv, years=16.5, n=2 * nk, record=MINREC, extra=(),
              regs=[reg_leup(11.25, 84.0, s, e)])
    g0 = (adult_height(run(MK_CPP, pN, years=16.5,
                           regs=[reg_leup(11.25, 84.0, s, e)]))
          - adult_height(run(MK_CPP, pN, years=16.5)))
    say(f"baseline gain = {g0:+.2f} cm")
    rows = []
    for a, kk in enumerate(keys):
        gl = adult_height(rTs, 2 * a) - adult_height(rUs, 2 * a)
        gh = adult_height(rTs, 2 * a + 1) - adult_height(rUs, 2 * a + 1)
        rows.append((kk, gl, gh, (gh - gl) / 4.0))
    say(f"{'parameter':12} {'-20%':>9} {'+20%':>9} {'d(gain)/10%':>12}")
    for kk, a, b, sv in sorted(rows, key=lambda z: -abs(z[3])):
        say(f"{kk:12} {a:+9.2f} {b:+9.2f} {sv:+12.3f}")

    # =====================================================================
    sec("16. HEAD-LINE SCENARIO TABLE")
    s = (7.4 - 5.0) * 365.25
    e = s + 365.25 * 5.0
    rU = run(MK_CPP, pN, years=16.5)
    rN = run(MK_NORM, pN, years=16.5)
    scen = [
        ("S1  normal puberty reference", rN),
        ("S2  CPP untreated", rU),
        ("S3  CPP + leuprolide 3.75 q28d",
         run(MK_CPP, pN, years=16.5, regs=[reg_leup(3.75, 28.0, s, e)])),
        ("S4  CPP + leuprolide 11.25 q12wk",
         run(MK_CPP, pN, years=16.5, regs=[reg_leup(11.25, 84.0, s, e)])),
        ("S5  CPP + histrelin implant",
         run(MK_CPP, pN, years=16.5, regs=[reg_hist(s, e)])),
        ("S6  CPP + triptorelin 22.5 q24wk",
         run(MK_CPP, pN, years=16.5, regs=[reg_trip(22.5, 168.0, s, e)])),
        ("S7  CPP + nafarelin nasal adh0.80",
         run(MK_CPP, pN, years=16.5, regs=[reg_naf(1800.0, s, e, 0.80, 3)])),
        ("S8  CPP + GnRH antagonist",
         run(MK_CPP, pN, years=16.5, regs=[reg_antag(18000.0, 28.0, s, e)])),
        ("S9  CPP + GnRHa, late start CA 10",
         run(MK_CPP, pN, years=16.5,
             regs=[reg_leup(11.25, 84.0, (10.0-5.0)*365.25, (10.0-5.0)*365.25+365.25*4)])),
        ("S10 CPP + GnRHa + rhGH",
         run(MK_CPP, pN, years=16.5, regs=[reg_leup(11.25, 84.0, s, e),
                                           reg_gh(0.043*30, s, e)])),
    ]
    say(f"{'scenario':36} {'thel':>6} {'menar':>6} {'fusion':>7} {'AdultHt':>8} "
        f"{'gain':>6} {'peakQoL':>8} {'peakHF':>7} {'BMDZnadir':>10}")
    for lab, r in scen:
        say(f"{lab:36} {thelarche_age(r):6.2f} {menarche_age(r):6.2f} "
            f"{fusion_age(r):7.2f} {adult_height(r):8.1f} "
            f"{adult_height(r)-adult_height(rU):+6.2f} {r['QOL'][:,0].max():8.2f} "
            f"{r['HF'][:,0].max():7.2f} {r['BMDZ'][:,0].min():+10.2f}")
    say()
    say(f"  reference normal-girl adult height = {adult_height(rN):.1f} cm")

    # =====================================================================
    sec("17. PK SANITY CHECK")
    for lab, regs, key in (
            ("leuprolide 3.75 mg q28d", [reg_leup(3.75, 28.0, 0.0, 365.0 * 3)], "CLEUP"),
            ("leuprolide 11.25 mg q12wk", [reg_leup(11.25, 84.0, 0.0, 365.0 * 3)], "CLEUP"),
            ("triptorelin 11.25 mg q12wk", [reg_trip(11.25, 84.0, 0.0, 365.0 * 3)], "CTRIP"),
            ("histrelin 50 mg implant", [reg_hist(0.0, 365.0 * 3)], "CHIST"),
            ("nafarelin 1800 ug/d nasal", [reg_naf(1800.0, 0.0, 365.0 * 3)], "CNAF")):
        r = run(MK_CPP, pN, years=4.0, regs=regs, stride=1)
        m = (r["t"] > 200) & (r["t"] < 330)
        say(f"{lab:28} mean Cp {r[key][m,0].mean():7.3f} ng/mL   "
            f"peak {r[key][m,0].max():7.3f}   trough {r[key][m,0].min():7.3f}   "
            f"mean fa {r['fa'][m,0].mean():.3f}")
    r = run(MK_CPP, pN, years=1.5, regs=[reg_ai(1.0, 0.0, 365.0)], stride=1)
    m = (r["t"] > 200) & (r["t"] < 330)
    say(f"{'anastrozole 1 mg/d':28} mean aromatase inhibition "
        f"{100*r['AIEFF'][m,0].mean():.1f}%")

    sec("DONE — every number above is produced by this file alone (numpy + RK4)")


if __name__ == "__main__":
    main()
    with open("cpp_reference_output.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
