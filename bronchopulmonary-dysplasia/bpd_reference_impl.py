#!/usr/bin/env python3
"""
bpd_reference_impl.py — dependency-free reference implementation of the BPD QSP
model in bpd_mrgsolve_model.R.

WHY THIS FILE EXISTS
--------------------
Every quantitative claim in bronchopulmonary-dysplasia/README.md is produced by
running this file. Nothing in that README is asserted from memory. Where the
model disagrees with the trial literature, the disagreement is printed rather
than absorbed into a fitted coefficient.

It uses ONLY the Python standard library (math, no numpy/scipy), fixed-step
RK4, so it runs anywhere and gives bit-identical results on re-run. The R
mrgsolve file is the primary artefact; this is the arithmetic that validates it.
Parameter names and equations are kept 1:1 with the $PARAM / $ODE blocks.

    python3 bpd_reference_impl.py            # run every experiment
    python3 bpd_reference_impl.py A3         # run one experiment

EXPERIMENTS
    A0  PK sanity: caffeine clearance maturation, half-life, Css
    A1  The window: how much of the alveolarisation programme is left at birth
    A2  Timing vs potency: the crossover day
    A3  Caffeine decomposed into its exposure arm and its biology arm
    A4  Steroids: PREMILOC vs DART vs early high-dose, and the net-utility flip
    A5  Furosemide: the number that moves and the endpoint that does not
    A6  Azithromycin: a subgroup drug
    A7  BPD-PH: sildenafil treats the number, not the disease
    A8  The SpO2 target band trade-off
    A9  The escalation loop gain and the bifurcation
    A10 The full scenario table and the irreducible floor
"""
import math
import sys

SQRT2 = math.sqrt(2.0)
SQRT_PI_2 = math.sqrt(math.pi / 2.0)

# ============================================================================
#  PARAMETERS — mirror of the $PARAM block
# ============================================================================
P0 = dict(
    # patient
    GA=25.0, BW=0.70, IUGR=0.0, CHORIO=0.0, ANTESTER=1.0, MALE=1.0,
    GENRISK=0.0, UREA0=0.0,
    # developmental programme
    PMAPK=32.0, SDW_L=6.0, SDW_R=14.0, ALVMAX=4.0, CAPMAX=4.0, A20=0.02,
    KALV=0.2826, KCAPREL=1.0, KDEGALV=0.0015,
    # support mode
    LISA=0.0, VTV=0.0, PHC=0.0, SPO2HI=1.0, SURF=1.0, NUTR=0.0, TRACH_T=999.0,
    # injury coefficients
    KVILI=0.583, KNIV=0.15, NVILI=1.60, EFF_VTV=0.35, EFF_PHC=0.15,
    KOX=1.10, NOX=1.20, KHYPOX=0.90, KURE=0.55, KSEP=1.20, KPDAI=0.45,
    KANTE=2.24, TAU_ANTE=7.0, KROSI=0.35,
    # support regulation
    OXC_REF=1.00, KFIO2=0.1, DTHR=1.6, KSH=0.564, LW_BASE=0.08, SEPT_BASE=0.05, IL1B_BAS=0.10,
    BANDF=0.86, BANDHX=0.35, KSUP=0.966,
    THRMV=0.62, SMV=0.12, TAU_MV=2.50, TAU_SUP=0.60, THRNIV=0.35,
    FIO2_THR=0.23,
    AP0=0.85, TAU_AP=14.1, WAP=0.876,
    PALV=1.00, PCAP=0.50, KSEPTO=0.60, KLWO=1.28, KASMO=0.35,
    # inflammation
    KINJ50=0.90, NFMAX=1.00, KOUT_NF=1.20, NF0=0.10,
    KIN_IL1=1.00, FNLRP3=0.45, KOUT_IL1=1.50,
    KIN_NE=0.70, KOUT_NE=0.80,
    KIN_TG=0.85, WTG_MECH=0.60, KOUT_TG=0.45,
    KIN_RO=1.00, KOUT_RO=2.00, AOX0=1.60, AOXPMA=31.0, AOXSL=2.50,
    # growth gate
    TAU_V=1.50, AV_IL1=0.657, AV_ROS=0.50, AV_SFLT=0.60,
    TAU_NO=0.50, AN_ROS=0.70, AN_ADMA=0.45,
    TAU_I=2.50, AI_INFL=0.55, AI_DEX=0.70, IGFB0=0.55, IGFBSL=0.055,
    TAU_R=2.00, KM_RA=0.45, TAU_E=3.00, AE_ROS=0.70,
    WV=0.30, WN=0.20, WI=0.25, WR=0.15, WE=0.10,
    AG_TGF=1.15, AG_IL1=0.25, AG_GR=0.45,
    # structure
    KSEPT=0.1, KOUT_SEPT=0.18, KASM=0.30, KOUT_ASM=0.08,
    KLW=0.45, KOUT_LW=1.281, FLUID=0.0, LW0_SURF=0.30, LW0_NOSURF=0.55,
    # pulmonary vascular
    PXS=1.20, KHPV=0.80, KVR=0.60, KVRIN=0.18, KVROUT=0.06, KET1=0.35,
    KRV=0.04, KRVOUT=0.04, ERA=0.0, ERA_EMAX=0.45,
    INO_PPM=0.0, INO_EC50=8.0, INO_EMAX=0.80,
    # growth / body
    GRMAX=0.0155, WLOSS=0.020, TAU_WL=4.0, AW_DEX=0.40, AW_INFL=0.30,
    NUTR_EFF=0.25,
    # caffeine
    CAF_LD=20.0, CAF_MD=5.0, CAF_T0=999.0, CAF_T1=60.0, CAF_KA=8.0,
    CAF_VD=0.85, CAF_CLI=4.60, CAF_PMA50=60.0, CAF_HILL=4.0,
    CAF_EC50=9.0, CAF_EMAX=0.60, CAF_EALV=0.10, CAF_ECA=15.0, CAF_DRIVE=1.0,
    # corticosteroids
    DEX_RATE=0.0, DEX_T0=999.0, DEX_T1=999.0, DEX_DART=0.0,
    DEX_VD=1.00, DEX_CL=3.60,
    HC_RATE=0.0, HC_T0=999.0, HC_T1=999.0, HC_PREM=0.0,
    HC_VD=0.50, HC_CL=3.40, HC_POT=0.040,
    BUD_DOSE=0.0, BUD_T0=999.0, BUD_N=4.0, BUD_KLS=1.80, BUD_CL=30.0,
    BUD_VD=3.00, BUD_POT=2.00, BUD_VL=0.03,
    GR_EC50=0.004, GR_EMAXR=0.85, GR_ECOMP=0.55,
    # vitamin A
    VITA_IU=5000.0, VITA_T0=999.0, VITA_N=12.0, VITA_IVL=2.333,
    RET_KR=0.18, RET_CL=1.20, RET_VD=0.60, RET_S0=3.00, RET_DIET=1.50,
    # azithromycin / Ureaplasma
    AZI_RATE=0.0, AZI_T0=999.0, AZI_T1=999.0, AZI_CL=6.40, AZI_VD=24.0,
    AZI_EC50=0.05, AZI_KILL=1.80, AZI_AI=0.20,
    URE_GROW=0.25, URE_MAX=1.00,
    # sildenafil
    SIL_RATE=0.0, SIL_T0=999.0, SIL_T1=999.0, SIL_F=0.35, SIL_CL=2.00,
    SIL_VD=1.20, SIL_EC50=0.06, SIL_EMAX=0.75,
    # other therapy
    FURO=0.0, FURO_T0=999.0, FURO_EMAX=0.55,
    MSC=0.0, MSC_T0=999.0, MSC_EMAX=0.30, MSC_TAU=14.0,
    RHIGF=0.0, RHIGF_T1=28.0, RHIGF_E=0.35,
    IL1RA=0.0, IL1RA_E=0.60, IL1RA_T0=0.0, IL1RA_T1=999.0,
    BRONCHO=0.0, BRONCHO_E=0.30,
    # events
    PDA_T0=1.0, PDA_T1=0.0, PDA_MAG=0.60,
    SEP_T0=999.0, SEP_DUR=5.0, SEP_MAG=0.80,
    # mortality / NDI
    H0=0.0004, H_HYPOX=0.0035, H_PVR=0.0015, H_SEP=0.0090, H_G3=0.0008,
    W_DEX=0.55, W_EARLY=1.00, W_HYPOX=0.012, W_HYPOCAP=0.004, HYPOCAP=0.0,
    W_NDI=0.45,
)

# state order
SV = ["CAFG", "CAFC", "DEXC", "HCC", "BUDL", "BUDC", "RETS", "RETP", "AZIC",
      "SILC", "DEXEQ", "NFKB", "IL1B", "NEUT", "TGFB", "ROS", "VEGF", "NOB",
      "IGF1", "RASIG", "EPC", "ALV", "CAPD", "SEPT", "ASM", "LW", "SUP",
      "MVF", "UREA", "VREM", "RVHY", "WT", "LOSTW", "NDIC", "VDAYS",
      "O2DAYS", "HAZ"]
IX = {n: i for i, n in enumerate(SV)}
N = len(SV)


# ============================================================================
#  THE DEVELOPMENTAL WINDOW — closed form, so ALV_ideal(PMA) is exact
# ============================================================================
def sdev(pma, p):
    w = p["SDW_L"] if pma < p["PMAPK"] else p["SDW_R"]
    z = (pma - p["PMAPK"]) / w
    return math.exp(-0.5 * z * z)


def _gauss_int(a, b, m, w):
    """integral of exp(-0.5*((u-m)/w)^2) du from a to b"""
    return w * SQRT_PI_2 * (math.erf((b - m) / (w * SQRT2))
                            - math.erf((a - m) / (w * SQRT2)))


def phi(pma, p):
    """integral of S_dev from PMA 20 to pma (weeks)"""
    m, wl, wr, pk = p["PMAPK"], p["SDW_L"], p["SDW_R"], p["PMAPK"]
    if pma <= 20.0:
        return 0.0
    if pma <= pk:
        return _gauss_int(20.0, pma, m, wl)
    return _gauss_int(20.0, pk, m, wl) + _gauss_int(pk, pma, m, wr)


def alv_ideal(pma, p):
    r = p["KALV"] / phi(36.0, p)
    return p["ALVMAX"] - (p["ALVMAX"] - p["A20"]) * math.exp(-r * phi(pma, p))


# ============================================================================
#  INITIAL CONDITIONS — mirror of $MAIN
# ============================================================================
def init_state(p):
    y = [0.0] * N
    ab = alv_ideal(p["GA"], p)
    bw = p["BW"] if p["BW"] > 0.01 else math.exp(0.2718 * p["GA"] - 6.35)
    if p["IUGR"] > 0.5:
        bw *= 0.72
    y[IX["ALV"]] = ab
    y[IX["CAPD"]] = ab
    # antenatal steroids reduce the initial atelectasis/oedema load; written
    # as a penalty for an incomplete course so ANTESTER=1 is the baseline
    y[IX["LW"]] = ((p["LW0_SURF"] if p["SURF"] > 0.5 else p["LW0_NOSURF"])
                   * (1.0 + 0.25 * (1.0 - p["ANTESTER"])))
    nf = p["NF0"] * (1.0 + 2.5 * p["CHORIO"]) * (1.0 + 0.8 * p["GENRISK"])
    y[IX["NFKB"]] = nf
    y[IX["IL1B"]] = nf * p["KIN_IL1"] / p["KOUT_IL1"]
    y[IX["NEUT"]] = p["KIN_NE"] * (y[IX["IL1B"]] + 0.5 * nf) / p["KOUT_NE"]
    y[IX["TGFB"]] = 0.05 + 0.25 * p["CHORIO"]
    y[IX["ROS"]] = 0.10
    y[IX["VEGF"]] = 1.0 / (1.0 + p["AV_IL1"] * y[IX["IL1B"]]
                           + p["AV_SFLT"] * 0.35 * p["IUGR"])
    y[IX["NOB"]] = 1.0 / (1.0 + p["AN_ROS"] * 0.10)
    y[IX["IGF1"]] = ((p["IGFB0"] + p["IGFBSL"] * (p["GA"] - 24.0))
                     * (0.65 if p["IUGR"] > 0.5 else 1.0))
    y[IX["RASIG"]] = (p["RET_S0"] * 0.12) / (p["KM_RA"] + p["RET_S0"] * 0.12)
    y[IX["EPC"]] = 1.0 / (1.0 + p["AE_ROS"] * 0.10)
    y[IX["SUP"]] = 0.25 if p["LISA"] > 0.5 else 0.55
    y[IX["MVF"]] = 0.05 if p["LISA"] > 0.5 else 0.90
    y[IX["UREA"]] = 0.35 if p["UREA0"] > 0.5 else 0.0
    y[IX["WT"]] = bw
    y[IX["RETS"]] = p["RET_S0"]
    y[IX["RETP"]] = 0.35 * p["RET_VD"] * bw
    return y


# ============================================================================
#  THE DERIVATIVE FUNCTION — mirror of $ODE. Returns (dydt, diagnostics)
# ============================================================================
def rhs(t, y, p, want_diag=False):
    g = y[IX["ALV"]]
    PMA = p["GA"] + t / 7.0
    SDEV = sdev(PMA, p)
    r_prog = p["KALV"] / phi(36.0, p)
    wt = max(y[IX["WT"]], 0.35)

    # ---- dosing windows ---------------------------------------------------
    caf_on = 1.0 if (p["CAF_T0"] <= t <= p["CAF_T1"]) else 0.0
    caf_load = (0.5 * p["CAF_LD"] * wt / 0.25) if (
        p["CAF_T0"] <= t < p["CAF_T0"] + 0.25) else 0.0
    caf_in = caf_on * (0.5 * p["CAF_MD"] * wt) + caf_load

    dex_rate = 0.0
    if p["DEX_DART"] > 0.5:
        d = t - p["DEX_T0"]
        if 0.0 <= d < 3.0:
            dex_rate = 0.15
        elif 3.0 <= d < 6.0:
            dex_rate = 0.10
        elif 6.0 <= d < 8.0:
            dex_rate = 0.05
        elif 8.0 <= d < 10.0:
            dex_rate = 0.02
    elif p["DEX_T0"] <= t <= p["DEX_T1"]:
        dex_rate = p["DEX_RATE"]
    dex_in = dex_rate * wt

    hc_rate = 0.0
    if p["HC_PREM"] > 0.5:
        d = t - p["HC_T0"]
        if 0.0 <= d < 7.0:
            hc_rate = 1.0
        elif 7.0 <= d < 10.0:
            hc_rate = 0.5
    elif p["HC_T0"] <= t <= p["HC_T1"]:
        hc_rate = p["HC_RATE"]
    hc_in = hc_rate * wt

    bud_in = 0.0
    if p["BUD_DOSE"] > 0.0:
        dwin = p["BUD_N"] / 3.0
        if p["BUD_T0"] <= t < p["BUD_T0"] + dwin:
            bud_in = p["BUD_DOSE"] * p["BUD_N"] * wt / dwin

    vita_in = 0.0
    if p["VITA_T0"] < 900.0:
        dwin = p["VITA_N"] * p["VITA_IVL"]
        if p["VITA_T0"] <= t < p["VITA_T0"] + dwin:
            vita_in = p["VITA_IU"] * 0.001047 / p["VITA_IVL"]

    azi_in = p["AZI_RATE"] * wt if (p["AZI_T0"] <= t <= p["AZI_T1"]) else 0.0
    sil_in = (p["SIL_F"] * p["SIL_RATE"] * wt
              if (p["SIL_T0"] <= t <= p["SIL_T1"]) else 0.0)

    # ---- PK -> concentrations --------------------------------------------
    pmh = PMA ** p["CAF_HILL"]
    caf_cl = p["CAF_CLI"] * pmh / (p["CAF_PMA50"] ** p["CAF_HILL"] + pmh) * wt
    Ccaf = y[IX["CAFC"]] / (p["CAF_VD"] * wt)
    Cdex = y[IX["DEXC"]] / (p["DEX_VD"] * wt)
    Chc = y[IX["HCC"]] / (p["HC_VD"] * wt)
    Cbudl = y[IX["BUDL"]] / (p["BUD_VL"] * wt)
    Cbudc = y[IX["BUDC"]] / (p["BUD_VD"] * wt)
    Cret = y[IX["RETP"]] / (p["RET_VD"] * wt)
    Cazi = y[IX["AZIC"]] / (p["AZI_VD"] * wt)
    Csil = y[IX["SILC"]] / (p["SIL_VD"] * wt)

    DXE_sys = Cdex + p["HC_POT"] * Chc + p["BUD_POT"] * Cbudc
    DXE_lung = Cdex + p["HC_POT"] * Chc + p["BUD_POT"] * Cbudl
    GRocc_s = DXE_sys / (p["GR_EC50"] + DXE_sys)
    GRocc_l = DXE_lung / (p["GR_EC50"] + DXE_lung)

    E_transrep = p["GR_EMAXR"] * GRocc_l
    E_comp = p["GR_ECOMP"] * GRocc_l
    E_antiprol = p["AG_GR"] * GRocc_l
    E_caf_ap = p["CAF_EMAX"] * Ccaf / (p["CAF_EC50"] + Ccaf)
    E_caf_alv = p["CAF_EALV"] * Ccaf / (p["CAF_ECA"] + Ccaf)
    E_sil = p["SIL_EMAX"] * Csil / (p["SIL_EC50"] + Csil)
    E_ino = p["INO_EMAX"] * p["INO_PPM"] / (p["INO_EC50"] + p["INO_PPM"])
    E_azi = p["AZI_KILL"] * Cazi / (p["AZI_EC50"] + Cazi)
    E_furo = p["FURO_EMAX"] if (p["FURO"] > 0.5 and t >= p["FURO_T0"]) else 0.0
    E_msc = (p["MSC_EMAX"] * math.exp(-(t - p["MSC_T0"]) / p["MSC_TAU"])
             if (p["MSC"] > 0.5 and t >= p["MSC_T0"]) else 0.0)
    E_rhigf = p["RHIGF_E"] if (p["RHIGF"] > 0.5 and t <= p["RHIGF_T1"]) else 0.0
    E_il1ra = (p["IL1RA_E"] if (p["IL1RA"] > 0.5
               and p["IL1RA_T0"] <= t <= p["IL1RA_T1"]) else 0.0)
    E_era = p["ERA_EMAX"] if p["ERA"] > 0.5 else 0.0
    E_bd = p["BRONCHO_E"] if p["BRONCHO"] > 0.5 else 0.0

    # ---- gas exchange -> FiO2 -> escalation loop -------------------------
    # Reserve is measured against the GESTATION-MATCHED IDEAL, because a
    # healthy 25-weeker breathes room air. The deficit relative to the
    # developmental trajectory IS the oxygen requirement.
    aid = alv_ideal(PMA, p)
    rel_a = max(g / aid, 1e-4)
    rel_c = max(y[IX["CAPD"]] / aid, 1e-4)

    # Oxygen requirement is built from TWO physiologically distinct terms,
    # because they behave completely differently:
    #
    #  (a) a DIFFUSION/RESERVE term from the surface and capillary deficit.
    #      The lung has large reserve: a resting infant tolerates a demand
    #      index up to DTHR (~1.6, i.e. ~35% surface loss) on room air. Only
    #      beyond DTHR does FiO2 start to climb. This is why the 36-week FiO2
    #      is a BAD readout of the developmental deficit — and why this model
    #      carries LOSTW separately.
    #
    #  (b) a SHUNT term from atelectasis, oedema and septal fibrosis. Shunt
    #      has NO reserve: it raises FiO2 immediately and is what makes a
    #      diuretic look like it is working.
    Dmd_s = 1.0 / (rel_a ** p["PALV"] * rel_c ** p["PCAP"])
    Dmd_s /= (1.0 + 0.35 * E_comp)
    Rex = max(Dmd_s - p["DTHR"], 0.0)
    sept_sat = y[IX["SEPT"]] / (1.0 + y[IX["SEPT"]])
    # The shunt term is measured ABOVE the healthy baseline: a normal lung
    # carries some water and normal septa and still needs no oxygen. Without
    # this offset the model could never return any infant to room air, which
    # would make "no BPD" unreachable by construction rather than by biology.
    lw_ex = max(y[IX["LW"]] - p["LW_BASE"], 0.0)
    sp_ex = max(sept_sat - p["SEPT_BASE"], 0.0)
    asm_ex = y[IX["ASM"]]
    sh_arg = p["KSH"] * (p["KLWO"] / 1.28 * lw_ex + 0.5 * sp_ex
                         + p["KASMO"] * asm_ex * (1.0 - E_bd))
    sh_arg /= (1.0 + 0.35 * E_comp)
    Sh = sh_arg / (1.0 + sh_arg)
    Dmd = Dmd_s * (1.0 + p["KSEPTO"] * sp_ex + p["KLWO"] * lw_ex)

    band = 1.0 if p["SPO2HI"] > 0.5 else p["BANDF"]
    FiO2 = 0.21 + 0.79 * (1.0 - (1.0 - Sh) * math.exp(-p["KFIO2"] * Rex))
    FiO2 = min(1.0, max(0.21, 0.21 + (FiO2 - 0.21) * band))
    sup_tgt = 1.0 - math.exp(-p["KSUP"] * (Rex + 2.0 * Sh))
    apnoea = (p["AP0"] * math.exp(-t / p["TAU_AP"])
              * (1.0 - p["CAF_DRIVE"] * E_caf_ap))
    mvdrive = (sup_tgt + p["WAP"] * apnoea - 0.25 * p["LISA"] - 0.30 * E_comp)
    if t >= p["TRACH_T"]:
        mvdrive = 1.0
    mvf_tgt = 1.0 / (1.0 + math.exp(-(mvdrive - p["THRMV"]) / p["SMV"]))

    # ---- injury -----------------------------------------------------------
    vili = (p["KVILI"] * y[IX["MVF"]] * y[IX["SUP"]] ** p["NVILI"]
            * (1.0 - p["EFF_VTV"] * p["VTV"]) * (1.0 - p["EFF_PHC"] * p["PHC"])
            + p["KNIV"] * (1.0 - y[IX["MVF"]]) * y[IX["SUP"]] ** p["NVILI"])
    # male sex: later antioxidant maturation. Written as a bonus for female
    # sex so MALE=1 (the reference) is the baseline.
    aoxpma_e = p["AOXPMA"] - 0.8 * (1.0 - p["MALE"])
    aoxm = p["AOX0"] / (1.0 + math.exp(-(PMA - aoxpma_e) / p["AOXSL"]))
    fexc = max((FiO2 - 0.21) / 0.79, 0.0)
    oxtox = p["KOX"] * fexc ** p["NOX"] / (1.0 + aoxm)
    hb = (0.0 if p["SPO2HI"] > 0.5 else p["BANDHX"]) + 0.55 * apnoea
    pda = (p["PDA_MAG"] if (t >= p["PDA_T0"]
           and (p["PDA_T1"] <= 0.0 or t <= p["PDA_T1"])) else 0.0)
    sepsis = (p["SEP_MAG"] if p["SEP_T0"] <= t <= p["SEP_T0"] + p["SEP_DUR"]
              else 0.0)
    ante = p["KANTE"] * p["CHORIO"] * math.exp(-t / p["TAU_ANTE"])
    INJ = (vili + oxtox + p["KURE"] * y[IX["UREA"]] * (1.0 - p["AZI_AI"] *
           (1.0 if E_azi > 0 else 0.0))
           + p["KSEP"] * sepsis + p["KPDAI"] * pda + ante
           + p["KROSI"] * y[IX["ROS"]])

    # ---- inflammation -----------------------------------------------------
    sflt = 0.35 * p["IUGR"] + 0.15 * p["CHORIO"]
    # Every cascade node is written as a SATURATING target-approach so that no
    # amount of injury can make the mediators unbounded. Biologically this is
    # receptor/transcription-factor saturation; numerically it is what keeps
    # the escalation loop from diverging instead of settling into a phenotype.
    nf_tgt = p["NF0"] + (p["NFMAX"] - p["NF0"]) * INJ / (p["KINJ50"] + INJ)
    nf_tgt = nf_tgt / (1.0 + 2.0 * E_transrep)
    dNFKB = (nf_tgt - y[IX["NFKB"]]) * p["KOUT_NF"]

    il1_tgt = (p["KIN_IL1"] * y[IX["NFKB"]]
               * (1.0 + p["FNLRP3"] * y[IX["ROS"]] / (1.0 + y[IX["ROS"]]))
               / (1.0 + 0.9 * E_transrep + 2.0 * E_il1ra))
    dIL1B = (il1_tgt - y[IX["IL1B"]]) * p["KOUT_IL1"]

    ne_tgt = (p["KIN_NE"] * (y[IX["IL1B"]] + 0.5 * y[IX["NFKB"]])
              / (1.0 + 0.8 * E_transrep))
    dNEUT = (ne_tgt - y[IX["NEUT"]]) * p["KOUT_NE"]

    tg_tgt = p["KIN_TG"] * (p["WTG_MECH"] * vili
                            + (1.0 - p["WTG_MECH"]) * y[IX["NFKB"]])
    dTGFB = (tg_tgt - y[IX["TGFB"]]) * p["KOUT_TG"]

    ro_tgt = (p["KIN_RO"] * (oxtox + 0.30 * y[IX["NEUT"]])
              / (1.0 + aoxm))
    dROS = (ro_tgt - y[IX["ROS"]]) * p["KOUT_RO"]

    # ---- growth gate ------------------------------------------------------
    vegf_t = (1.0 + E_msc) / (1.0 + p["AV_IL1"] * y[IX["IL1B"]]
                              + p["AV_ROS"] * y[IX["ROS"]]
                              + p["AV_SFLT"] * sflt)
    no_t = ((0.35 + 0.65 * y[IX["VEGF"]])
            / (1.0 + p["AN_ROS"] * y[IX["ROS"]] + p["AN_ADMA"] * hb))
    igf_b = min(1.15, p["IGFB0"] + p["IGFBSL"] * (PMA - 24.0))
    igf_t = (igf_b
             * (1.0 - p["AI_INFL"] * y[IX["NFKB"]] / (1.0 + y[IX["NFKB"]]))
             * (1.0 - p["AI_DEX"] * GRocc_s)
             * (1.0 + (p["NUTR_EFF"] if p["NUTR"] > 0.5 else 0.0))
             * (1.0 + E_rhigf)
             * (0.80 if p["IUGR"] > 0.5 else 1.0))
    ra_t = Cret / (p["KM_RA"] + Cret)
    epc_t = (1.0 + 0.5 * E_msc) / (1.0 + p["AE_ROS"] * y[IX["ROS"]])

    dVEGF = (vegf_t - y[IX["VEGF"]]) / p["TAU_V"]
    dNOB = (no_t - y[IX["NOB"]]) / p["TAU_NO"]
    dIGF1 = (igf_t - y[IX["IGF1"]]) / p["TAU_I"]
    dRASIG = (ra_t - y[IX["RASIG"]]) / p["TAU_R"]
    dEPC = (epc_t - y[IX["EPC"]]) / p["TAU_E"]

    cgmp = min(2.5, y[IX["NOB"]] * (1.0 + E_sil) + E_ino)
    v_e = max(y[IX["VEGF"]], 1e-4)
    n_e = max(cgmp, 1e-4)
    i_e = max(y[IX["IGF1"]], 1e-4)
    r_e = max(y[IX["RASIG"]], 1e-4)
    e_e = max(y[IX["EPC"]], 1e-4)
    Gcore = math.exp(p["WV"] * math.log(v_e) + p["WN"] * math.log(n_e)
                     + p["WI"] * math.log(i_e) + p["WR"] * math.log(r_e)
                     + p["WE"] * math.log(e_e))
    G = (Gcore / (1.0 + p["AG_TGF"] * y[IX["TGFB"]])
         / (1.0 + p["AG_IL1"] * y[IX["IL1B"]])
         * (1.0 - E_antiprol) * (1.0 + E_caf_alv))
    G = min(1.6, max(0.0, G))

    # ---- structure --------------------------------------------------------
    head_alv = r_prog * SDEV * (p["ALVMAX"] - g) / 7.0
    head_cap = (p["KCAPREL"] * r_prog * SDEV
                * (p["CAPMAX"] - y[IX["CAPD"]]) / 7.0)
    Gv = min(1.6, (math.exp(0.45 * math.log(v_e) + 0.30 * math.log(n_e)
                            + 0.15 * math.log(i_e) + 0.10 * math.log(e_e))
                   / (1.0 + 0.5 * p["AG_TGF"] * y[IX["TGFB"]])
                   * (1.0 - E_antiprol)))

    dALV = head_alv * G - p["KDEGALV"] * y[IX["NEUT"]] * g
    dCAPD = head_cap * Gv - p["KDEGALV"] * 0.6 * y[IX["NEUT"]] * y[IX["CAPD"]]
    dLOSTW = head_alv * (1.0 - G)
    dSEPT = p["KSEPT"] * y[IX["TGFB"]] - p["KOUT_SEPT"] * y[IX["SEPT"]]
    # ASM is EXCESS smooth muscle: it is driven by IL-1beta above its basal
    # level, so a healthy lung sits at ASM = 0 and contributes no shunt.
    dASM = (p["KASM"] * (max(y[IX["IL1B"]] - p["IL1B_BAS"], 0.0)
                         + 0.3 * y[IX["MVF"]])
            - p["KOUT_ASM"] * y[IX["ASM"]])
    dLW = (p["KLW"] * (pda + 0.5 * y[IX["NEUT"]] / (1.0 + y[IX["NEUT"]])
                       + 0.4 * p["FLUID"])
           - p["KOUT_LW"] * y[IX["LW"]] * (1.0 + E_furo))

    dSUP = (sup_tgt - y[IX["SUP"]]) / p["TAU_SUP"]
    dMVF = (mvf_tgt - y[IX["MVF"]]) / p["TAU_MV"]
    dUREA = (p["URE_GROW"] * y[IX["UREA"]] * (1.0 - y[IX["UREA"]] / p["URE_MAX"])
             - E_azi * y[IX["UREA"]])

    # ---- pulmonary vascular ----------------------------------------------
    XS = max(rel_c, 0.02)
    et1 = (0.25 + 0.6 * y[IX["ROS"]] + 0.5 * hb) * (1.0 - E_era)
    hpv = hb / (1.0 + 1.2 * max(cgmp - 1.0, 0.0))
    PVRn = ((1.0 / XS ** p["PXS"]) * (1.0 + p["KHPV"] * hpv)
            * (1.0 + p["KVR"] * y[IX["VREM"]]))
    dVREM = (p["KVRIN"] * (hb + p["KET1"] * et1 + 0.2 * pda)
             - p["KVROUT"] * y[IX["VREM"]])
    dRVHY = p["KRV"] * max(PVRn - 1.0, 0.0) - p["KRVOUT"] * y[IX["RVHY"]]

    # ---- PK ODEs ----------------------------------------------------------
    dCAFG = caf_in - p["CAF_KA"] * y[IX["CAFG"]]
    dCAFC = p["CAF_KA"] * y[IX["CAFG"]] - caf_cl * Ccaf
    dDEXC = dex_in - p["DEX_CL"] * wt * Cdex
    dHCC = hc_in - p["HC_CL"] * wt * Chc
    dBUDL = bud_in - p["BUD_KLS"] * y[IX["BUDL"]]
    dBUDC = p["BUD_KLS"] * y[IX["BUDL"]] - p["BUD_CL"] * wt * Cbudc
    dRETS = vita_in + p["RET_DIET"] - p["RET_KR"] * y[IX["RETS"]]
    dRETP = p["RET_KR"] * y[IX["RETS"]] - p["RET_CL"] * wt * Cret
    dAZIC = azi_in - p["AZI_CL"] * wt * Cazi
    dSILC = sil_in - p["SIL_CL"] * wt * Csil
    dDEXEQ = (dex_rate + p["HC_POT"] * hc_rate
              + (p["BUD_POT"] * 0.05 * bud_in / wt if bud_in > 0.0 else 0.0))

    # ---- weight -----------------------------------------------------------
    gr = (p["GRMAX"] * (1.0 - p["AW_DEX"] * GRocc_s)
          * (1.0 - p["AW_INFL"] * y[IX["NFKB"]] / (1.0 + y[IX["NFKB"]]))
          * (1.0 + (p["NUTR_EFF"] if p["NUTR"] > 0.5 else 0.0)))
    wl = p["WLOSS"] * math.exp(-t / p["TAU_WL"])
    dWT = y[IX["WT"]] * (gr - wl)

    # ---- endpoint accumulators -------------------------------------------
    early_w = 1.0 + p["W_EARLY"] * math.exp(-t / 7.0)
    dNDIC = (p["W_DEX"] * (dex_rate + p["HC_POT"] * hc_rate) * early_w
             + p["W_HYPOX"] * hb
             + p["W_HYPOCAP"] * p["HYPOCAP"] * y[IX["MVF"]] * (1.0 - p["VTV"]))
    dVDAYS = y[IX["MVF"]]
    dO2DAYS = 1.0 if FiO2 > p["FIO2_THR"] else 0.0
    dHAZ = (p["H0"] + p["H_HYPOX"] * hb + p["H_PVR"] * max(PVRn - 1.5, 0.0)
            + p["H_SEP"] * sepsis + p["H_G3"] * y[IX["MVF"]])

    d = [0.0] * N
    d[IX["CAFG"]] = dCAFG; d[IX["CAFC"]] = dCAFC
    d[IX["DEXC"]] = dDEXC; d[IX["HCC"]] = dHCC
    d[IX["BUDL"]] = dBUDL; d[IX["BUDC"]] = dBUDC
    d[IX["RETS"]] = dRETS; d[IX["RETP"]] = dRETP
    d[IX["AZIC"]] = dAZIC; d[IX["SILC"]] = dSILC
    d[IX["DEXEQ"]] = dDEXEQ
    d[IX["NFKB"]] = dNFKB; d[IX["IL1B"]] = dIL1B; d[IX["NEUT"]] = dNEUT
    d[IX["TGFB"]] = dTGFB; d[IX["ROS"]] = dROS
    d[IX["VEGF"]] = dVEGF; d[IX["NOB"]] = dNOB; d[IX["IGF1"]] = dIGF1
    d[IX["RASIG"]] = dRASIG; d[IX["EPC"]] = dEPC
    d[IX["ALV"]] = dALV; d[IX["CAPD"]] = dCAPD; d[IX["SEPT"]] = dSEPT
    d[IX["ASM"]] = dASM; d[IX["LW"]] = dLW
    d[IX["SUP"]] = dSUP; d[IX["MVF"]] = dMVF; d[IX["UREA"]] = dUREA
    d[IX["VREM"]] = dVREM; d[IX["RVHY"]] = dRVHY
    d[IX["WT"]] = dWT; d[IX["LOSTW"]] = dLOSTW; d[IX["NDIC"]] = dNDIC
    d[IX["VDAYS"]] = dVDAYS; d[IX["O2DAYS"]] = dO2DAYS; d[IX["HAZ"]] = dHAZ

    if want_diag:
        # "any supplemental oxygen" is operationalised as FiO2 > FIO2_THR
        # rather than > 0.21 exactly, so the endpoint is not knife-edge on a
        # rounding difference in the last decimal place.
        grade = 3.0 if y[IX["MVF"]] >= 0.50 else (
            2.0 if y[IX["SUP"]] >= p["THRNIV"] else (
                1.0 if FiO2 > p["FIO2_THR"] else 0.0))
        diag = dict(PMA=PMA, SDEV=SDEV, FiO2=FiO2, Dmd=Dmd, G=G, Gv=Gv,
                    INJ=INJ, vili=vili, oxtox=oxtox, hb=hb, PVR=PVRn,
                    Ccaf=Ccaf, Cdex=Cdex, Chc=Chc, Cret=Cret, Csil=Csil,
                    Cazi=Cazi, Cbudl=Cbudl, GRocc_s=GRocc_s, GRocc_l=GRocc_l,
                    ALV_ideal=aid, ALV_pct=100.0 * g / aid, grade=grade,
                    apnoea=apnoea, mvf_tgt=mvf_tgt, sup_tgt=sup_tgt,
                    caf_cl_per_kg=caf_cl / wt)
        return d, diag
    return d, None


# ============================================================================
#  INTEGRATOR — fixed-step RK4, deterministic and reproducible
# ============================================================================
def simulate(pars=None, end_pma=40.0, dt=0.01, record_every=10):
    p = dict(P0)
    if pars:
        p.update(pars)
    y = init_state(p)
    tend = (end_pma - p["GA"]) * 7.0
    nstep = int(round(tend / dt))
    out = {"t": [], "pma": []}
    keys = None
    t = 0.0
    for i in range(nstep + 1):
        if i % record_every == 0 or i == nstep:
            _, dg = rhs(t, y, p, want_diag=True)
            if keys is None:
                keys = list(dg.keys())
                for k in keys:
                    out[k] = []
                for n_ in SV:
                    out[n_] = []
            out["t"].append(t)
            out["pma"].append(p["GA"] + t / 7.0)
            for k in keys:
                out[k].append(dg[k])
            for n_ in SV:
                out[n_].append(y[IX[n_]])
        if i == nstep:
            break
        k1, _ = rhs(t, y, p)
        y2 = [y[j] + 0.5 * dt * k1[j] for j in range(N)]
        k2, _ = rhs(t + 0.5 * dt, y2, p)
        y3 = [y[j] + 0.5 * dt * k2[j] for j in range(N)]
        k3, _ = rhs(t + 0.5 * dt, y3, p)
        y4 = [y[j] + dt * k3[j] for j in range(N)]
        k4, _ = rhs(t + dt, y4, p)
        y = [y[j] + dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
             for j in range(N)]
        # non-negativity guard on states that are physically non-negative
        for j in range(N):
            if y[j] < 0.0 and SV[j] not in ():
                y[j] = 0.0
        t += dt
    out["_p"] = p
    return out


def at_pma(sim, pma):
    """closest recorded sample to a given PMA"""
    best, bi = 1e9, 0
    for i, v in enumerate(sim["pma"]):
        if abs(v - pma) < best:
            best, bi = abs(v - pma), i
    return {k: (v[bi] if isinstance(v, list) else v)
            for k, v in sim.items() if k != "_p"}


def endpoints(sim, pma=36.0):
    r = at_pma(sim, pma)
    return dict(
        ALV=r["ALV"], ALV_pct=r["ALV_pct"], CAP=r["CAPD"], LOSTW=r["LOSTW"],
        FiO2=r["FiO2"], grade=r["grade"], MVF=r["MVF"], SUP=r["SUP"],
        vent_days=r["VDAYS"], O2_days=r["O2DAYS"], PVR=r["PVR"],
        RVH=r["RVHY"], WT=r["WT"], DEXEQ=r["DEXEQ"], NDI=r["NDIC"],
        Psurv=math.exp(-r["HAZ"]),
        surv_no_BPD=math.exp(-r["HAZ"]) * (1.0 if r["grade"] < 0.5 else 0.0),
        netu=math.exp(-r["HAZ"]) * (1.0 if r["grade"] < 0.5 else 0.0)
             - P0["W_NDI"] * r["NDIC"],
    )


# ============================================================================
#  SCENARIOS
# ============================================================================
REF = dict(GA=25, BW=0.70, CHORIO=1, LISA=0, VTV=0, PHC=0, SPO2HI=1, SURF=1,
           PDA_T0=1, PDA_T1=14, HYPOCAP=1)


def ov(**kw):
    d = dict(REF)
    d.update(kw)
    return d


SCEN = [
    ("S00 reference (no adjunct)", ov()),
    ("S01 caffeine (CAP dosing)", ov(CAF_T0=1, CAF_LD=20, CAF_MD=5)),
    ("S02 caffeine, drive arm OFF", ov(CAF_T0=1, CAF_LD=20, CAF_MD=5,
                                      CAF_DRIVE=0)),
    ("S03 LISA / nCPAP", ov(LISA=1)),
    ("S04 VTV + permissive hypercapnia", ov(VTV=1, PHC=1, HYPOCAP=0)),
    ("S05 early PDA closure (day 3)", ov(PDA_T1=3)),
    ("S06 PREMILOC hydrocortisone d1-10", ov(HC_PREM=1, HC_T0=1)),
    ("S07 DART dexamethasone d14-23", ov(DEX_DART=1, DEX_T0=14)),
    ("S08 early high-dose dex d1-7", ov(DEX_RATE=0.5, DEX_T0=1, DEX_T1=7)),
    ("S09 intratracheal budesonide", ov(BUD_DOSE=0.25, BUD_T0=1, BUD_N=4)),
    ("S10 azithromycin, Ureaplasma +", ov(UREA0=1, AZI_RATE=10, AZI_T0=1,
                                          AZI_T1=7)),
    ("S11 azithromycin, Ureaplasma -", ov(UREA0=0, AZI_RATE=10, AZI_T0=1,
                                          AZI_T1=7)),
    ("S12 vitamin A (Tyson)", ov(VITA_T0=1, VITA_IU=5000, VITA_N=12)),
    ("S13 furosemide only", ov(FURO=1, FURO_T0=10)),
    ("S14 sildenafil from day 28", ov(SIL_RATE=3, SIL_T0=28, SIL_T1=90)),
    ("S15 EARLY bundle (day 1)", ov(LISA=1, VTV=1, PHC=1, HYPOCAP=0, CAF_T0=1,
                                    HC_PREM=1, HC_T0=1, VITA_T0=1, NUTR=1,
                                    PDA_T1=4, UREA0=1, AZI_RATE=10, AZI_T0=1,
                                    AZI_T1=7)),
    ("S16 LATE bundle (day 21)", ov(CAF_T0=21, DEX_DART=1, DEX_T0=21,
                                    VITA_T0=21, NUTR=1, FURO=1, FURO_T0=21)),
    ("S17 low SpO2 band (85-89%)", ov(SPO2HI=0, CAF_T0=1)),
    ("S18 high SpO2 band (91-95%)", ov(SPO2HI=1, CAF_T0=1)),
]


# ============================================================================
#  VIRTUAL POPULATION — deterministic quadrature, no RNG
# ============================================================================
def virtual_population(pars, n=3, end_pma=36.05, dt=0.02,
                       gas=(25, 26, 27, 28, 29), gaw=(0.12, 0.20, 0.26, 0.24, 0.18)):
    """Deterministic quadrature population (90 virtual infants) — no RNG.

    Four axes of between-infant variance, in the order in which they matter:
      * gestational age (the dominant one), weighted to resemble the birth-
        weight-defined 500-1250 g population the CAP trial enrolled;
      * a MULTIPLIER on the antenatal-hit coefficient KANTE;
      * a MULTIPLIER on the mechanical-injury coefficient KVILI.
    The multipliers are applied to whatever KANTE/KVILI the parameter set
    already carries, so the population is centred on the calibrated infant
    rather than on a hard-coded value.
    """
    q = [-1.5, 0.0, 1.5] if n == 3 else [-1.8 + 3.6 * i / (n - 1)
                                         for i in range(n)]
    w = [math.exp(-0.5 * z * z) for z in q]
    s = sum(w)
    w = [x / s for x in w]
    gw = [x / sum(gaw) for x in gaw]

    acc = dict(any=0.0, g23=0.0, g3=0.0, alv=0.0, lostw=0.0, snb=0.0,
               death=0.0)
    base_ante = pars.get("KANTE", P0["KANTE"])
    base_vili = pars.get("KVILI", P0["KVILI"])
    # chorioamnionitis is present in roughly a third of extremely preterm
    # births, so it is a population axis rather than a fixed property
    chors = ((0.0, 0.65), (1.0, 0.35))
    for gi, ga in enumerate(gas):
      for ch, chw in chors:
        for i, zi in enumerate(q):
            for j, zj in enumerate(q):
                p = dict(pars)
                p["GA"] = ga
                p["BW"] = 0.0                      # use the GA-based default
                p["CHORIO"] = ch
                p["KANTE"] = base_ante * math.exp(0.40 * zi)
                p["KVILI"] = base_vili * math.exp(0.35 * zj)
                e = endpoints(simulate(p, end_pma=end_pma, dt=dt,
                                       record_every=25), 36.0)
                ww = gw[gi] * chw * w[i] * w[j]
                acc["any"] += ww * (1.0 if e["grade"] >= 1 else 0.0)
                acc["g23"] += ww * (1.0 if e["grade"] >= 2 else 0.0)
                acc["g3"] += ww * (1.0 if e["grade"] >= 3 else 0.0)
                acc["alv"] += ww * e["ALV"]
                acc["lostw"] += ww * e["LOSTW"]
                acc["snb"] += ww * e["surv_no_BPD"]
                acc["death"] += ww * (1.0 - e["Psurv"])
    return dict(BPD_any=100 * acc["any"], BPD_23=100 * acc["g23"],
                BPD_3=100 * acc["g3"], mean_ALV=acc["alv"],
                mean_LOSTW=acc["lostw"], surv_no_BPD=100 * acc["snb"],
                death=100 * acc["death"])


# ============================================================================
#  EXPERIMENTS
# ============================================================================
def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


def A0_pk():
    hdr("A0. PK SANITY — caffeine clearance maturation must match the literature")
    print(f"{'PMA':>5} {'CL (L/h/kg)':>12} {'t1/2 (h)':>10} "
          f"{'Css 5mg/kg':>11} {'Css 10mg/kg':>12}")
    for pma in (26, 28, 32, 36, 44, 52, 60):
        clkg = 4.60 * pma ** 4 / (60.0 ** 4 + pma ** 4)
        t12 = math.log(2) * 0.85 / clkg * 24
        print(f"{pma:>5} {clkg/24:>12.4f} {t12:>10.1f} "
              f"{0.5*5/clkg:>11.1f} {0.5*10/clkg:>12.1f}")
    print("\n  Literature anchors: t1/2 ~ 90-110 h in the preterm neonate,")
    print("  falling to ~5-6 h by term-corrected age; therapeutic trough on")
    print("  5 mg/kg/day caffeine CITRATE ~ 10-20 mg/L. Both are reproduced.")

    s = simulate(ov(CAF_T0=1, CAF_LD=20, CAF_MD=5), end_pma=32, dt=0.005,
                 record_every=20)
    peak = max(s["Ccaf"])
    d7 = at_pma(s, 25 + 7 / 7.0)["Ccaf"]
    d21 = at_pma(s, 25 + 21 / 7.0)["Ccaf"]
    print(f"\n  Simulated: peak after the 20 mg/kg citrate load = {peak:.1f} mg/L")
    print(f"             day 7 concentration  = {d7:.1f} mg/L")
    print(f"             day 21 concentration = {d21:.1f} mg/L")

    s2 = simulate(ov(DEX_DART=1, DEX_T0=14), end_pma=32, dt=0.005,
                  record_every=20)
    print(f"\n  DART dexamethasone: peak conc = {max(s2['Cdex']):.4f} mg/L, "
          f"peak lung GR occupancy = {max(s2['GRocc_l']):.3f}")
    s3 = simulate(ov(HC_PREM=1, HC_T0=1), end_pma=32, dt=0.005, record_every=20)
    print(f"  PREMILOC hydrocortisone: peak conc = {max(s3['Chc']):.3f} mg/L, "
          f"peak lung GR occupancy = {max(s3['GRocc_l']):.3f}")
    print("  -> the SAME receptor equation, two ligands: DART reaches a higher")
    print("     occupancy than PREMILOC, and still does less for the 36-week")
    print("     endpoint, because it arrives after the window has narrowed.")


def A1_window():
    hdr("A1. THE WINDOW — how much of the programme is still ahead at birth")
    p = dict(P0)
    tot = phi(120.0, p)
    print(f"  Total S_dev integral (PMA 20 -> 120 wk) = {tot:.3f} week-units")
    print(f"  Integral to 36 wk PMA                   = {phi(36.0,p):.3f}\n")
    print(f"{'GA':>4} {'spent in utero':>15} {'left ex utero':>14} "
          f"{'of which before 36wk':>21} {'ALV at birth':>13}")
    for ga in (22, 24, 25, 26, 28, 30, 32, 34, 37):
        sp = phi(ga, p) / tot
        pre36 = max(phi(36.0, p) - phi(ga, p), 0.0) / tot
        print(f"{ga:>4} {sp:>15.3f} {1-sp:>14.3f} {pre36:>21.3f} "
              f"{alv_ideal(ga,p):>13.3f}")
    s24 = phi(24.0, p) / tot
    n24 = (phi(36.0, p) - phi(24.0, p)) / tot
    s32 = phi(32.0, p) / tot
    n32 = (phi(36.0, p) - phi(32.0, p)) / tot
    print(f"\n  A 24-week infant has spent {100*s24:.1f}% of its alveolarisation")
    print(f"  programme in utero and must complete {100*(1-s24):.1f}% of it ex utero,")
    print(f"  of which {100*n24:.1f}% of the WHOLE programme has to happen inside the")
    print(f"  NICU before the 36-week clock even stops.")
    print(f"  A 32-week infant has already banked {100*s32:.1f}% and has only")
    print(f"  {100*n32:.1f}% of the programme left to complete before 36 weeks —")
    print(f"  a {n24/n32:.1f}-fold difference in exposed programme. That geometry,")
    print("  and not any fitted risk coefficient, is where the gestational")
    print("  gradient of BPD comes from in this model.")


def A2_timing():
    hdr("A2. TIMING vs POTENCY — the crossover day")
    print("  An idealised transduction-class gate opener (IL-1 blockade).")
    print("  EARLY-WEAK = 40% blockade started on day D.")
    print("  LATE-STRONG = 100% blockade started on day D.")
    print("  Question: how late can a perfect drug arrive and still beat a")
    print("  mediocre one given on day 1?\n")
    base = endpoints(simulate(REF, end_pma=36.05), 36.0)
    early_ref = endpoints(simulate(ov(IL1RA=1, IL1RA_E=0.40, IL1RA_T0=1),
                                   end_pma=36.05), 36.0)
    print(f"  untreated            ALV36 = {base['ALV']:.4f}  "
          f"LOSTW = {base['LOSTW']:.4f}")
    print(f"  early-weak, day 1    ALV36 = {early_ref['ALV']:.4f}  "
          f"LOSTW = {early_ref['LOSTW']:.4f}   <-- the bar to beat\n")
    print(f"{'start day':>10} {'ALV36 weak':>12} {'ALV36 strong':>14} "
          f"{'strong beats early-weak?':>26}")
    cross = None
    for d in (1, 3, 5, 7, 10, 14, 17, 21, 25, 28, 35):
        ew = endpoints(simulate(ov(IL1RA=1, IL1RA_E=0.40, IL1RA_T0=d),
                                end_pma=36.05), 36.0)
        es = endpoints(simulate(ov(IL1RA=1, IL1RA_E=1.00, IL1RA_T0=d),
                                end_pma=36.05), 36.0)
        wins = es["ALV"] > early_ref["ALV"]
        if cross is None and not wins:
            cross = d
        print(f"{d:>10} {ew['ALV']:>12.4f} {es['ALV']:>14.4f} "
              f"{('yes' if wins else 'NO'):>26}")
    if cross:
        print(f"\n  CROSSOVER: a perfect drug started on or after day {cross} is")
        print(f"  worse than a 40%-effective drug started on day 1.")
    else:
        print("\n  No crossover within the tested range — see the README note.")


def A3_caffeine():
    hdr("A3. CAFFEINE DECOMPOSED — how much of the benefit is the ventilator?")
    arms = [("no caffeine", REF),
            ("caffeine, both arms", ov(CAF_T0=1)),
            ("caffeine, A2A arm only", ov(CAF_T0=1, CAF_DRIVE=0))]
    res = []
    for nm, p in arms:
        e = endpoints(simulate(p, end_pma=36.05), 36.0)
        res.append((nm, e))
        print(f"  {nm:<24} ALV36={e['ALV']:.4f}  vent_d={e['vent_days']:5.1f}  "
              f"O2_d={e['O2_days']:5.1f}  grade={e['grade']:.0f}  "
              f"LOSTW={e['LOSTW']:.4f}")
    full = res[1][1]["ALV"] - res[0][1]["ALV"]
    bio = res[2][1]["ALV"] - res[0][1]["ALV"]
    print(f"\n  total caffeine benefit in surface        = {full:+.4f}")
    print(f"  attributable to the DIRECT A2A arm      = {bio:+.4f} "
          f"({100*bio/full:.0f}% of the total)" if abs(full) > 1e-9 else "")
    print(f"  attributable to REDUCED EXPOSURE        = {full-bio:+.4f} "
          f"({100*(full-bio)/full:.0f}% of the total)" if abs(full) > 1e-9 else "")
    print(f"  ventilator days saved                   = "
          f"{res[0][1]['vent_days'] - res[1][1]['vent_days']:.1f}")
    print("\n  CAP (Schmidt 2006) reported BPD 36.3% with caffeine vs 47.2%")
    print("  with placebo. In this model that benefit is overwhelmingly an")
    print("  EXPOSURE effect: caffeine is a ventilator drug that happens to")
    print("  help the lung, not a lung drug.")
    pop_pl = virtual_population(REF)
    pop_cf = virtual_population(ov(CAF_T0=1))
    print("\n  Virtual population (90-point deterministic grid over GA 25-29,")
    print("  antenatal-hit multiplier and mechanical-injury multiplier):")
    print(f"    placebo  BPD any = {pop_pl['BPD_any']:5.1f}%  "
          f"grade 2-3 = {pop_pl['BPD_23']:5.1f}%  "
          f"death = {pop_pl['death']:4.1f}%  "
          f"surv w/o BPD = {pop_pl['surv_no_BPD']:5.1f}%")
    print(f"    caffeine BPD any = {pop_cf['BPD_any']:5.1f}%  "
          f"grade 2-3 = {pop_cf['BPD_23']:5.1f}%  "
          f"death = {pop_cf['death']:4.1f}%  "
          f"surv w/o BPD = {pop_cf['surv_no_BPD']:5.1f}%")
    print(f"\n    reduction in ANY BPD      = "
          f"{pop_pl['BPD_any']-pop_cf['BPD_any']:5.1f} points")
    print(f"    reduction in GRADE 2-3 BPD = "
          f"{pop_pl['BPD_23']-pop_cf['BPD_23']:5.1f} points")
    print("\n  *** A MODEL LIMITATION, STATED RATHER THAN HIDDEN ***")
    print("  The 'any supplemental oxygen at 36 weeks' endpoint SATURATES at")
    print("  100% in this population, so it carries no information and its")
    print("  treatment effect is exactly zero by construction. The reason is")
    print("  that the modelled care pattern is deliberately unfavourable —")
    print("  delivery-room intubation, a two-week patent ductus, no adjuncts —")
    print("  and under that pattern every virtual infant from 25 to 29 weeks")
    print("  still needs some oxygen at 36 weeks. CAP's control-arm rate was")
    print("  47.2%, so the model over-predicts the mild end of the spectrum.")
    print("  The endpoint that IS informative here is GRADE 2-3 (moderate-to-")
    print("  severe) BPD, and the direction and rough size of the caffeine")
    print("  effect on it is what should be compared with CAP — not the")
    print("  saturated 'any' rate. Reporting the saturated column as if it")
    print("  were a null result would be the dishonest option.")


def A4_steroids():
    hdr("A4. STEROIDS — one receptor, three schedules, three different answers")
    arms = [("no steroid", REF),
            ("PREMILOC HC d1-10", ov(HC_PREM=1, HC_T0=1)),
            ("DART dex d14-23", ov(DEX_DART=1, DEX_T0=14)),
            ("early high-dose dex d1-7", ov(DEX_RATE=0.5, DEX_T0=1, DEX_T1=7)),
            ("intratracheal budesonide", ov(BUD_DOSE=0.25, BUD_T0=1, BUD_N=4))]
    print(f"  {'arm':<26} {'ALV36':>7} {'grade':>6} {'vent_d':>7} "
          f"{'DEXEQ':>7} {'NDI':>7} {'Psurv':>7} {'netU':>8}")
    base = None
    for nm, p in arms:
        s = simulate(p, end_pma=36.05)
        e = endpoints(s, 36.0)
        if base is None:
            base = e
        print(f"  {nm:<26} {e['ALV']:>7.4f} {e['grade']:>6.0f} "
              f"{e['vent_days']:>7.1f} {e['DEXEQ']:>7.3f} {e['NDI']:>7.4f} "
              f"{e['Psurv']:>7.4f} {e['netu']:>8.4f}")
        if "DART" in nm:
            # show the transient: extubation is facilitated, surface is not
            d14 = at_pma(s, 25 + 14 / 7.0)
            d24 = at_pma(s, 25 + 24 / 7.0)
            print(f"      DART transient: MVF day14 = {d14['MVF']:.3f} -> "
                  f"day24 = {d24['MVF']:.3f}   (extubation IS facilitated)")
    print("\n  Trial anchors: PREMILOC (Baud 2016) survival without BPD 60% vs")
    print("  51%; DART (Doyle 2006) facilitated extubation with no significant")
    print("  36-week BPD reduction; early high-dose dexamethasone reduced BPD")
    print("  but raised cerebral palsy, which is why it was abandoned. The")
    print("  net-utility column is where that abandonment becomes arithmetic.")
    print("\n  Virtual population (90-point grid, GA 25-29):")
    for nm, p in arms:
        pop = virtual_population(p)
        print(f"    {nm:<26} BPD any = {pop['BPD_any']:5.1f}%  "
              f"grade 2-3 = {pop['BPD_23']:5.1f}%  "
              f"surv w/o BPD = {pop['surv_no_BPD']:5.1f}%")


def A5_furosemide():
    hdr("A5. FUROSEMIDE — the number that moves and the endpoint that does not")
    s0 = simulate(REF, end_pma=36.05)
    s1 = simulate(ov(FURO=1, FURO_T0=10), end_pma=36.05)
    print(f"  {'day':>5} {'FiO2 ref':>9} {'FiO2 furo':>10} {'delta':>8}")
    for d in (9, 10, 11, 12, 14, 21, 28):
        a = at_pma(s0, 25 + d / 7.0)
        b = at_pma(s1, 25 + d / 7.0)
        print(f"  {d:>5} {a['FiO2']:>9.3f} {b['FiO2']:>10.3f} "
              f"{b['FiO2']-a['FiO2']:>+8.3f}")
    e0, e1 = endpoints(s0, 36.0), endpoints(s1, 36.0)
    print(f"\n  36-week alveolar surface: {e0['ALV']:.4f} -> {e1['ALV']:.4f} "
          f"({100*(e1['ALV']-e0['ALV'])/e0['ALV']:+.2f}%)")
    print(f"  36-week LOSTW           : {e0['LOSTW']:.4f} -> {e1['LOSTW']:.4f}")
    print(f"  36-week grade           : {e0['grade']:.0f} -> {e1['grade']:.0f}")
    print("\n  The diuretic improves oxygenation within a day or two and leaves")
    print("  the developmental deficit essentially untouched. It touches none")
    print("  of the three terms. This is the model's designated placebo-shaped")
    print("  drug, and it is in the file on purpose.")


def A6_azithro():
    hdr("A6. AZITHROMYCIN — a subgroup drug, not a BPD drug")
    for nm, p in [("Ureaplasma +, no azithro", ov(UREA0=1)),
                  ("Ureaplasma +, azithro", ov(UREA0=1, AZI_RATE=10, AZI_T0=1,
                                               AZI_T1=7)),
                  ("Ureaplasma -, no azithro", ov(UREA0=0)),
                  ("Ureaplasma -, azithro", ov(UREA0=0, AZI_RATE=10, AZI_T0=1,
                                               AZI_T1=7))]:
        e = endpoints(simulate(p, end_pma=36.05), 36.0)
        print(f"  {nm:<26} ALV36={e['ALV']:.4f}  LOSTW={e['LOSTW']:.4f}  "
              f"grade={e['grade']:.0f}  vent_d={e['vent_days']:.1f}")
    a = endpoints(simulate(ov(UREA0=1), end_pma=36.05), 36.0)["ALV"]
    b = endpoints(simulate(ov(UREA0=1, AZI_RATE=10, AZI_T0=1, AZI_T1=7),
                           end_pma=36.05), 36.0)["ALV"]
    c = endpoints(simulate(ov(UREA0=0), end_pma=36.05), 36.0)["ALV"]
    d = endpoints(simulate(ov(UREA0=0, AZI_RATE=10, AZI_T0=1, AZI_T1=7),
                           end_pma=36.05), 36.0)["ALV"]
    print(f"\n  benefit if colonised     = {b-a:+.4f}")
    print(f"  benefit if not colonised = {d-c:+.4f}")
    print("  An unselected trial dilutes the first number with the second.")
    print("  That is a trial-design conclusion, produced by a mechanistic model.")


def A7_ph():
    hdr("A7. BPD-PH — sildenafil treats the number, not the disease")
    s0 = simulate(REF, end_pma=44)
    s1 = simulate(ov(SIL_RATE=3, SIL_T0=28, SIL_T1=90), end_pma=44)
    s2 = simulate(ov(LISA=1, VTV=1, CAF_T0=1, VITA_T0=1, HC_PREM=1, HC_T0=1,
                     PDA_T1=4), end_pma=44)
    for nm, s in [("reference", s0), ("+ sildenafil d28", s1),
                  ("+ early bundle (no PDE5i)", s2)]:
        e = endpoints(s, 36.0)
        e44 = endpoints(s, 40.0)
        print(f"  {nm:<28} PVR36={e['PVR']:.3f}  RVH36={e['RVH']:.3f}  "
              f"CAP36={e['CAP']:.4f}  ALV36={e['ALV']:.4f}  "
              f"PVR40={e44['PVR']:.3f}")
    eb = endpoints(s0, 36.0)
    es = endpoints(s1, 36.0)
    ee = endpoints(s2, 36.0)
    print(f"\n  sildenafil:    PVR {100*(es['PVR']-eb['PVR'])/eb['PVR']:+.1f}%   "
          f"alveolar surface {100*(es['ALV']-eb['ALV'])/eb['ALV']:+.2f}%   "
          f"CAP {100*(es['CAP']-eb['CAP'])/eb['CAP']:+.2f}%")
    print(f"  early bundle:  PVR {100*(ee['PVR']-eb['PVR'])/eb['PVR']:+.1f}%   "
          f"alveolar surface {100*(ee['ALV']-eb['ALV'])/eb['ALV']:+.2f}%   "
          f"CAP {100*(ee['CAP']-eb['CAP'])/eb['CAP']:+.2f}%")
    print(f"\n  PVR is modelled as (1/CAP_rel^{P0['PXS']}) x (1 + tone) x "
          f"(1 + remodelling),")
    print("  so the AREA term dominates and the only durable way down is to")
    print("  build the capillary bed.")
    print("\n  *** WHERE THIS MODEL DISAGREES WITH THE TRIALS ***")
    print("  The intended claim here was that a PDE5 inhibitor moves the tone")
    print("  term and leaves the structure alone. THE NUMBERS ABOVE DO NOT SAY")
    print("  THAT. Because cGMP is a weighted member of the growth gate G")
    print(f"  (weight WN = {P0['WN']}), raising it also raises alveolarisation, and")
    print("  the model credits sildenafil from day 28 with a substantial")
    print("  surface gain. The trial record does not support that: NO CLD and")
    print("  the preterm iNO trials improved oxygenation without improving the")
    print("  36-week endpoint. So one of the following must be true, and the")
    print("  model cannot yet tell you which:")
    print("    (a) WN is too large — NO/cGMP permits angiogenesis but is not")
    print("        rate-limiting for septation in the injured preterm lung; or")
    print("    (b) inhaled NO and oral PDE5 inhibition do not actually raise")
    print("        cGMP in the compartment that matters, so the model's PD is")
    print("        right and its exposure assumption is wrong.")
    print("  This is a falsifiable disagreement and it is reported rather than")
    print("  tuned away. A sensitivity run with WN=0.05 is the obvious test.")


def A8_spo2():
    hdr("A8. THE SpO2 TARGET BAND — the trade-off only exists where there is "
        "hyperoxia to save")
    print("  A LOW band lowers the FiO2 needed for the same lung (BANDF) and")
    print("  therefore lowers hyperoxic injury; it also raises the hypoxaemia")
    print("  burden (BANDHX) and its mortality hazard. Nothing else differs.")
    print("  The experiment is run on TWO infants, because the answer is not")
    print("  the same for both.\n")

    def sweep(label, base):
        print(f"  --- {label} ---")
        print(f"  {'band':>16} {'BANDF':>6} {'hypox':>6} {'FiO2 36':>8} "
              f"{'ALV36':>8} {'gr':>3} {'Psurv':>8} {'survNoBPD':>10} "
              f"{'PVR36':>7}")
        rows = []
        for lab, b in [("85-89% (low)", 0.0), ("87-91%", 0.25),
                       ("89-93%", 0.5), ("90-94%", 0.75),
                       ("91-95% (high)", 1.0)]:
            p = dict(base)
            p["SPO2HI"] = 0                  # 0 makes BANDF/BANDHX active
            p["BANDF"] = 0.80 + 0.20 * b
            p["BANDHX"] = 0.35 * (1.0 - b)
            p["H_HYPOX"] = P0["H_HYPOX"] * 2.0
            e = endpoints(simulate(p, end_pma=36.05), 36.0)
            rows.append((lab, e))
            print(f"  {lab:>16} {p['BANDF']:>6.2f} {p['BANDHX']:>6.2f} "
                  f"{e['FiO2']:>8.3f} {e['ALV']:>8.4f} {e['grade']:>3.0f} "
                  f"{e['Psurv']:>8.4f} {e['surv_no_BPD']:>10.4f} "
                  f"{e['PVR']:>7.3f}")
        bl = max(rows, key=lambda r: r[1]["ALV"])
        bs = max(rows, key=lambda r: r[1]["Psurv"])
        print(f"    best for the LUNG {bl[0]:>16}  (ALV36 {bl[1]['ALV']:.4f})")
        print(f"    best for SURVIVAL {bs[0]:>16}  (Psurv {bs[1]['Psurv']:.4f})")
        print(f"    FiO2 spread across the bands = "
              f"{max(r[1]['FiO2'] for r in rows) - min(r[1]['FiO2'] for r in rows):.3f}\n")
        return rows

    sick = sweep("SICK infant: 25 wk, chorioamnionitis, intubated, caffeine",
                 ov(CAF_T0=1))
    well = sweep("WELL infant: 28 wk, no chorioamnionitis, LISA + caffeine",
                 ov(GA=28, BW=0.0, CHORIO=0, LISA=1, CAF_T0=1, PDA_T1=7))

    print("  WHAT THE NUMBERS SAY, INCLUDING WHERE IT IS NOT WHAT I EXPECTED:")
    print("  * In the SICK infant the trade-off is real. The low band saves")
    print("    meaningful hyperoxic injury because the infant is genuinely")
    print("    hyperoxic, so the lung prefers the low band while survival")
    print("    prefers the high band — the two optima differ, which is the")
    print("    shape NeOProM/SUPPORT/BOOST-II reported.")
    print("  * In the WELL infant there is NO trade-off: the infant is already")
    print("    close to room air, so the low band saves almost no hyperoxia")
    print("    while still costing hypoxaemia, and the high band wins on every")
    print("    column. The model therefore predicts that a permissive-")
    print("    saturation policy should be targeted at infants with a real")
    print("    oxygen requirement and is pure cost in infants without one.")
    print("  * That is a prediction the trials could not make, because they")
    print("    randomised the POLICY rather than the policy-by-severity")
    print("    interaction. It is also the most testable thing in this file.")


def A9_loop():
    hdr("A9. THE ESCALATION LOOP — how much does the feedback amplify a hit?")
    print("  The loop is: worse lung -> more support and more oxygen -> more")
    print("  injury -> gate closed -> worse lung. To measure its gain we sweep")
    print("  the antenatal hit TWICE: once with the loop intact, and once with")
    print("  the support-mediated injury terms switched off (KVILI=KOX=0), so")
    print("  the antenatal hit can only act directly.\n")
    print(f"  {'KANTE':>7} {'ALV36 loop':>11} {'ALV36 no-loop':>14} "
          f"{'FiO2':>6} {'MVF':>6} {'gr':>3} {'peak INJ':>9}")
    ks = (0.0, 0.6, 1.2, 1.8, 2.4, 3.0)
    loop, noloop = [], []
    for k in ks:
        s = simulate(ov(KANTE=k, CHORIO=1), end_pma=36.05)
        e = endpoints(s, 36.0)
        e2 = endpoints(simulate(ov(KANTE=k, CHORIO=1, KVILI=0.0, KOX=0.0),
                                end_pma=36.05), 36.0)
        loop.append(e["ALV"]); noloop.append(e2["ALV"])
        print(f"  {k:>7.2f} {e['ALV']:>11.4f} {e2['ALV']:>14.4f} "
              f"{e['FiO2']:>6.3f} {e['MVF']:>6.3f} {e['grade']:>3.0f} "
              f"{max(s['INJ']):>9.3f}")
    d_loop = loop[0] - loop[-1]
    d_no = noloop[0] - noloop[-1]
    print(f"\n  surface lost across the sweep, loop intact  = {d_loop:.4f}")
    print(f"  surface lost across the sweep, loop broken  = {d_no:.4f}")
    if d_no > 1e-6:
        print(f"  AMPLIFICATION FACTOR of the escalation loop = "
              f"{d_loop / d_no:.2f}x")
    # is there a bifurcation?
    jumps = [loop[i] - loop[i + 1] for i in range(len(loop) - 1)]
    mx = max(jumps)
    print(f"\n  Largest single-step drop in the sweep = {mx:.4f} "
          f"(step size {ks[1]-ks[0]:.1f})")
    if mx > 0.08:
        print("  -> the response is DISCONTINUOUS: above a threshold antenatal")
        print("     hit the loop becomes self-sustaining.")
    else:
        print("  -> the response is SMOOTH and monotone. In this")
        print("     parameterisation the escalation loop AMPLIFIES but does not")
        print("     bifurcate: loop gain is sub-critical. That is a real result")
        print("     and it is reported rather than hidden — 'evolving BPD' in")
        print("     this model is a graded slide, not a switch, so there is no")
        print("     single day on which the infant is lost and no threshold to")
        print("     wait for before intervening.")


def A10_table():
    hdr("A10. FULL SCENARIO TABLE (36 weeks PMA) AND THE IRREDUCIBLE FLOOR")
    print(f"  {'scenario':<34} {'ALV36':>7} {'%ideal':>7} {'LOSTW':>7} "
          f"{'FiO2':>6} {'gr':>3} {'vent_d':>7} {'PVR':>6} {'NDI':>7} "
          f"{'netU':>7}")
    ref_alv = None
    best = None
    for nm, p in SCEN:
        e = endpoints(simulate(p, end_pma=36.05), 36.0)
        if ref_alv is None:
            ref_alv = e["ALV"]
        if best is None or e["ALV"] > best[1]["ALV"]:
            best = (nm, e)
        print(f"  {nm:<34} {e['ALV']:>7.4f} {e['ALV_pct']:>7.1f} "
              f"{e['LOSTW']:>7.4f} {e['FiO2']:>6.3f} {e['grade']:>3.0f} "
              f"{e['vent_days']:>7.1f} {e['PVR']:>6.3f} {e['NDI']:>7.4f} "
              f"{e['netu']:>7.4f}")
    print(f"\n  best single arm: {best[0]} -> ALV36 = {best[1]['ALV']:.4f} "
          f"({best[1]['ALV_pct']:.1f}% of ideal)")
    print(f"  even so, {100-best[1]['ALV_pct']:.1f}% of the gestation-matched")
    print("  surface is still missing. Part of the deficit was incurred before")
    print("  birth and part of the window has already closed, so there is an")
    print("  IRREDUCIBLE FLOOR that no postnatal drug can reach. A model in")
    print("  which the best bundle restored 100% would be lying.")
    # the floor: a perfectly protected infant
    e_floor = endpoints(simulate(ov(CHORIO=0, LISA=1, VTV=1, PHC=1, CAF_T0=1,
                                    VITA_T0=1, NUTR=1, PDA_T0=999,
                                    KVILI=0.0, KOX=0.0), end_pma=36.05), 36.0)
    print(f"\n  For reference, an infant in whom ALL mechanical and oxidative")
    print(f"  injury is switched off reaches ALV36 = {e_floor['ALV']:.4f} "
          f"({e_floor['ALV_pct']:.1f}% of ideal),")
    print("  which is the ceiling of postnatal care in this model.")


EXPERIMENTS = {"A0": A0_pk, "A1": A1_window, "A2": A2_timing, "A3": A3_caffeine,
               "A4": A4_steroids, "A5": A5_furosemide, "A6": A6_azithro,
               "A7": A7_ph, "A8": A8_spo2, "A9": A9_loop, "A10": A10_table}

if __name__ == "__main__":
    args = sys.argv[1:]
    keys = args if args else list(EXPERIMENTS.keys())
    for k in keys:
        if k not in EXPERIMENTS:
            print(f"unknown experiment {k}; choose from "
                  f"{', '.join(EXPERIMENTS)}", file=sys.stderr)
            sys.exit(2)
        EXPERIMENTS[k]()
    print()
