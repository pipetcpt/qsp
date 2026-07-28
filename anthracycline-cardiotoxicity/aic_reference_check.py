#!/usr/bin/env python3
# =============================================================================
# aic_reference_check.py
# Independent numpy transcription of the anthracycline-induced cardiotoxicity
# (AIC) QSP model in aic_mrgsolve_model.R.
#
# WHY THIS FILE EXISTS
# --------------------
# Every number quoted in README.md and in the header of the mrgsolve model is
# an OUTPUT of these equations, not a literature value copied into prose. This
# script is a second, independent implementation (fixed-step RK4, vectorised
# over a virtual population) whose job is to (a) produce those numbers and
# (b) let anyone re-derive them without R/mrgsolve installed.
#
#   python3 aic_reference_check.py            # all analyses (A0..A11)
#   python3 aic_reference_check.py --quick    # index-patient scenarios only
#
# 32 ODE states. Time unit = DAYS. Concentrations: mg/L (plasma),
# arbitrary tissue units normalised to per-cycle reference exposure.
# =============================================================================
import argparse
import numpy as np

# -----------------------------------------------------------------------------
# state index map
# -----------------------------------------------------------------------------
(A1, A2, A3, ALIP, AM, AMP, CHF, CH, CHM, ADEX, T2B, TR1, TR2,
 EACE, EBB, ESTA, EARNI, ESGLT,
 DSB, P53, ROS, MITOD, LIP, MYO, FUNC, FIB, HYP, NH,
 TNI, BNP, AUCH, CUMKILL) = range(32)
NSTATE = 32
SNAMES = ["A1", "A2", "A3", "ALIP", "AM", "AMP", "CHF", "CH", "CHM", "ADEX",
          "T2B", "TR1", "TR2", "EACE", "EBB", "ESTA", "EARNI", "ESGLT",
          "DSB", "P53", "ROS", "MITOD", "LIP", "MYO", "FUNC", "FIB", "HYP",
          "NH", "TNI", "BNP", "AUCH", "CUMKILL"]

# -----------------------------------------------------------------------------
# parameters
# -----------------------------------------------------------------------------
P = dict(
    # --- doxorubicin PK (3-cmt, adult BSA 1.8 m2; AUC 60 mg/m2 = 2.0 ug*h/mL)
    CL=1296.0, V1=60.0, Q2=1200.0, V2=200.0, Q3=600.0, V3=1300.0,
    FM=0.25,                       # fraction of elimination -> doxorubicinol
    CLM=150.0, VM=250.0, QM=200.0, VMP=400.0,
    # --- pegylated liposomal carrier (t1/2 ~3.1 d, 27% released as free drug)
    KEL_LIP=0.22, KREL_LIP=0.11,
    # --- cardiac distribution: TWO pools, deliberately -------------------
    KPINF=600.0, KPHF=25.0,        # fast nuclear/free pool -> tracks PEAKS
    KPIN=3.0,   KPH=40.0,          # slow retained pool     -> tracks AUC
    KPINM=1.5,  KPHM=60.0,         # retained doxorubicinol (accumulates)
    CHFREF=10.0, CHREF=0.25, CHMREF=0.27,
    WM=1.0,                        # doxorubicinol redox potency weight
    # --- dexrazoxane ----------------------------------------------------
    KDEX=6.6, VDEX=25.0, KDEX50=8.0,
    KT2B_SYN=0.35, KT2B_DEG=25.0, KCHEL=2.6,
    # --- trastuzumab ----------------------------------------------------
    KT10=0.075, KT12=0.20, KT21=0.222, VTR=3.0, EC50TR=25.0,
    # --- cardioprotective effect-state onset (titration ~2-3 weeks) -----
    KON=0.20,
    # --- Top2b / DSB arm; P53 is a SLOW genotoxic-stress memory ----------
    KDSB=31.0, TOXN50=2.50, KDSBOUT=0.35, WTR_REP=0.90,
    KP53=0.40, KP53OUT=0.025,
    # --- redox / iron arm (loop gain kept sub-critical) ------------------
    KROS=0.50, WN=0.55, WR=0.45, FEAMP=0.35, WMITO=1.20, WP53R=0.50,
    KROSOUT=1.20,
    KLIP=0.15, KLIPOUT=0.15, WLIPROS=0.35,
    # --- mitochondrial deficit ------------------------------------------
    KMIN=0.012, KMOUT=0.030, WBIOG=1.50,
    # --- irreversible myocyte pool --------------------------------------
    KDROS=4.8e-4, ROS50=1.00, WP53D=1.20, WTRD=0.65, KDNH=5.0e-5,
    KREG=2.0e-5, WTR_REG=0.80,
    # --- reversible functional injury -----------------------------------
    KFIN=1.0e-3, WF_ROS=0.70, WF_TOX=0.30, KFTR=2.2e-3,
    KFOUT=0.012, WFIBF=3.0,
    # --- fibrosis / remodelling -----------------------------------------
    KFIBIN=8.0e-3, WFIB_DEAD=2.5, WFIB_NH=0.30, KFIBOUT=0.010,
    # --- compensatory hypertrophy (the masking term) --------------------
    GH=1.60, HYPMAX=0.32, KH=1.0 / 45.0,
    # --- neurohormonal --------------------------------------------------
    KNH=0.15, WWS_NH=3.0,
    # --- contractility -> LVEF ------------------------------------------
    EF0=62.0, WFIBC=1.20, EFEXP=0.90,
    GLS0=20.5, WG_FUNC=0.90, WG_DEAD=1.60, WG_FIB=0.35,
    # --- biomarkers -----------------------------------------------------
    KTIN=55.0, KTOUT=1.40, WT_ROS=0.35, TNI_BASE=1.5,
    KBIN=40.0, KBOUT=1.20, WB_WS=14.0,
    # --- cardioprotection pharmacology (fractional target engagement) ---
    W_STA_SCAV=0.65, W_BB_SCAV=0.15, W_SGLT_SCAV=0.25,
    W_ACE_NH=0.60, W_BB_NH=0.50, W_ARNI_NH=1.30,
    W_ACE_REC=0.45, W_BB_REC=0.35, W_ARNI_REC=1.00, W_SGLT_REC=0.30,
    W_ACE_FIB=0.25, W_ARNI_FIB=0.50,
)

BSA = 1.8


def hill(x, x50, n):
    xn = np.power(np.maximum(x, 0.0), n)
    return xn / (np.power(x50, n) + xn)


# -----------------------------------------------------------------------------
# right-hand side (vectorised over subjects; y is (NSTATE, nsub))
# -----------------------------------------------------------------------------
def rhs(t, y, p, rate, tgt):
    d = np.zeros_like(y)

    C1 = y[A1] / p["V1"]
    CM = y[AM] / p["VM"]
    CDEX = y[ADEX] / p["VDEX"]
    CTR = y[TR1] / p["VTR"]

    # ---- doxorubicin / carrier / metabolite PK --------------------------
    k10 = p["CL"] / p["V1"]
    k12 = p["Q2"] / p["V1"]
    k21 = p["Q2"] / p["V2"]
    k13 = p["Q3"] / p["V1"]
    k31 = p["Q3"] / p["V3"]
    rel = p["KREL_LIP"] * y[ALIP]
    d[A1] = (-(k10 + k12 + k13) * y[A1] + k21 * y[A2] + k31 * y[A3]
             + rel + rate["dox"])
    d[A2] = k12 * y[A1] - k21 * y[A2]
    d[A3] = k13 * y[A1] - k31 * y[A3]
    d[ALIP] = -p["KEL_LIP"] * y[ALIP] + rate["lipo"]

    klm = p["CLM"] / p["VM"]
    km12 = p["QM"] / p["VM"]
    km21 = p["QM"] / p["VMP"]
    d[AM] = p["FM"] * k10 * y[A1] - (klm + km12) * y[AM] + km21 * y[AMP]
    d[AMP] = km12 * y[AM] - km21 * y[AMP]

    # ---- cardiac pools --------------------------------------------------
    d[CHF] = p["KPINF"] * (C1 - y[CHF] / p["KPHF"])
    d[CH] = p["KPIN"] * (C1 - y[CH] / p["KPH"])
    d[CHM] = p["KPINM"] * (CM - y[CHM] / p["KPHM"])

    TOXN = y[CHF] / p["CHFREF"]                       # nuclear, peak-driven
    TOXR = y[CH] / p["CHREF"] + p["WM"] * y[CHM] / p["CHMREF"]   # redox, AUC

    # ---- dexrazoxane ----------------------------------------------------
    d[ADEX] = -p["KDEX"] * y[ADEX] + rate["dex"]
    FDEX = CDEX / (CDEX + p["KDEX50"])
    # Top2b protein turnover: dexrazoxane drives proteasomal degradation
    # within hours; resynthesis takes days (t1/2 = 2 d).
    d[T2B] = p["KT2B_SYN"] * (1.0 - y[T2B]) - p["KT2B_DEG"] * FDEX * y[T2B]
    TOP2B = y[T2B]
    CHEL = p["KCHEL"] * FDEX

    # ---- trastuzumab ----------------------------------------------------
    d[TR1] = -(p["KT10"] + p["KT12"]) * y[TR1] + p["KT21"] * y[TR2] + rate["tr"]
    d[TR2] = p["KT12"] * y[TR1] - p["KT21"] * y[TR2]
    ETR = CTR / (CTR + p["EC50TR"])

    # ---- cardioprotective effect states ---------------------------------
    for i, key in ((EACE, "ace"), (EBB, "bb"), (ESTA, "sta"),
                   (EARNI, "arni"), (ESGLT, "sglt")):
        d[i] = p["KON"] * (tgt[key] - y[i])

    # ---- Top2b-dependent DNA damage -------------------------------------
    d[DSB] = (p["KDSB"] * TOP2B * hill(TOXN, p["TOXN50"], 2)
              - p["KDSBOUT"] / (1.0 + p["WTR_REP"] * ETR) * y[DSB])
    d[P53] = p["KP53"] * y[DSB] - p["KP53OUT"] * y[P53]

    # ---- labile iron / ROS ----------------------------------------------
    d[LIP] = (p["KLIP"] * (1.0 + p["WLIPROS"] * y[ROS])
              - p["KLIPOUT"] * y[LIP] - CHEL * y[LIP])
    SCAV = 1.0 / (1.0 + p["W_STA_SCAV"] * y[ESTA] + p["W_BB_SCAV"] * y[EBB]
                  + p["W_SGLT_SCAV"] * y[ESGLT])
    d[ROS] = (p["KROS"] * (p["WN"] * hill(TOXN, p["TOXN50"], 2)
                           + p["WR"] * TOXR * (1.0 + p["FEAMP"] * y[LIP])
                           + p["WMITO"] * y[MITOD]
                           + p["WP53R"] * y[P53]) * SCAV
              - p["KROSOUT"] * y[ROS])

    BIOG = 1.0 / (1.0 + p["WBIOG"] * y[P53])
    d[MITOD] = (p["KMIN"] * (0.6 * y[ROS] + 0.4 * y[P53]) * (1.0 - y[MITOD])
                - p["KMOUT"] * BIOG * y[MITOD])

    # ---- neurohormonal tone (blocked fraction) --------------------------
    NHe = y[NH] / (1.0 + p["W_ACE_NH"] * y[EACE] + p["W_BB_NH"] * y[EBB]
                   + p["W_ARNI_NH"] * y[EARNI])

    # ---- irreversible myocyte loss --------------------------------------
    kdeath = (p["KDROS"] * hill(y[ROS], p["ROS50"], 3)
              * (1.0 + p["WP53D"] * y[P53]) * (1.0 + p["WTRD"] * ETR)
              + p["KDNH"] * NHe)
    d[MYO] = -kdeath * y[MYO] + p["KREG"] * (1.0 - y[MYO]) * (1.0 - p["WTR_REG"] * ETR)

    # ---- reversible functional injury ------------------------------------
    RECOV = (p["W_ACE_REC"] * y[EACE] + p["W_BB_REC"] * y[EBB]
             + p["W_ARNI_REC"] * y[EARNI] + p["W_SGLT_REC"] * y[ESGLT])
    d[FUNC] = ((p["KFIN"] * (p["WF_ROS"] * y[ROS] + p["WF_TOX"] * TOXR)
                + p["KFTR"] * ETR) * (1.0 - y[FUNC])
               - p["KFOUT"] * (1.0 + RECOV) / (1.0 + p["WFIBF"] * y[FIB]) * y[FUNC])

    # ---- fibrosis --------------------------------------------------------
    DEAD = 1.0 - y[MYO]
    d[FIB] = (p["KFIBIN"] * (p["WFIB_DEAD"] * DEAD + p["WFIB_NH"] * NHe)
              * (1.0 - y[FIB])
              - p["KFIBOUT"] * (1.0 + p["W_ACE_FIB"] * y[EACE]
                                + p["W_ARNI_FIB"] * y[EARNI]) * y[FIB])

    # ---- compensatory hypertrophy ----------------------------------------
    CONTraw = y[MYO] * (1.0 - y[FUNC])
    DEFICIT = np.maximum(0.0, 1.0 - CONTraw)
    HYPT = np.minimum(p["HYPMAX"], p["GH"] * DEFICIT)
    d[HYP] = p["KH"] * (HYPT - y[HYP])

    CONT = CONTraw * (1.0 + y[HYP]) / (1.0 + p["WFIBC"] * y[FIB])
    LVEF = p["EF0"] * np.power(np.maximum(CONT, 1e-6), p["EFEXP"])
    WS = np.maximum(0.0, (p["EF0"] - LVEF) / p["EF0"]) + 0.5 * y[FUNC]

    d[NH] = p["KNH"] * (p["WWS_NH"] * WS - y[NH])

    # ---- biomarkers -------------------------------------------------------
    d[TNI] = (p["KTIN"] * (1000.0 * kdeath * y[MYO] + p["WT_ROS"] * y[ROS])
              - p["KTOUT"] * y[TNI])
    d[BNP] = p["KBIN"] * (1.0 + p["WB_WS"] * WS) - p["KBOUT"] * y[BNP]

    d[AUCH] = y[CH]
    d[CUMKILL] = kdeath * y[MYO]
    return d


def derived(y, p):
    """Observables from a state array (NSTATE, ...)."""
    CONTraw = y[MYO] * (1.0 - y[FUNC])
    CONT = CONTraw * (1.0 + y[HYP]) / (1.0 + p["WFIBC"] * y[FIB])
    LVEF = p["EF0"] * np.power(np.maximum(CONT, 1e-6), p["EFEXP"])
    GLS = p["GLS0"] * (1.0 - p["WG_FUNC"] * y[FUNC]
                       - p["WG_DEAD"] * (1.0 - y[MYO])
                       - p["WG_FIB"] * y[FIB])
    return dict(LVEF=LVEF, GLS=GLS, TNI=y[TNI] + p["TNI_BASE"], BNP=y[BNP],
                MYO=y[MYO], FUNC=y[FUNC], FIB=y[FIB], ROS=y[ROS],
                HYP=y[HYP], P53=y[P53], MITOD=y[MITOD], CH=y[CH],
                CHF=y[CHF], CHM=y[CHM])


# -----------------------------------------------------------------------------
# regimen construction
# -----------------------------------------------------------------------------
def dox_reg(dose_mg_m2, ncyc, interval, tinf=0.0, start=0.0, lipo=False):
    """IV push (tinf=0) or continuous infusion of `tinf` days per cycle."""
    ev = []
    cmt = "lipo" if lipo else "dox"
    for i in range(ncyc):
        t0 = start + i * interval
        amt = dose_mg_m2 * BSA
        if tinf > 0:
            ev.append(dict(kind="inf", t0=t0, t1=t0 + tinf, rate=amt / tinf,
                           cmt=cmt))
        else:
            ev.append(dict(kind="bolus", t0=t0, amt=amt, cmt=cmt))
    return ev


def dex_reg(dox_mg_m2, ncyc, interval, ratio=10.0, start=0.0):
    """Dexrazoxane 10:1 by rapid IV push immediately before each anthracycline."""
    return [dict(kind="bolus", t0=max(0.0, start + i * interval - 0.02),
                 amt=ratio * dox_mg_m2 * BSA, cmt="dex") for i in range(ncyc)]


def tras_reg(start, weeks=52, wt=70.0):
    """Trastuzumab 8 mg/kg load then 6 mg/kg q3w."""
    ev = [dict(kind="bolus", t0=start, amt=8.0 * wt, cmt="tr")]
    for i in range(1, int(weeks / 3)):
        ev.append(dict(kind="bolus", t0=start + 21.0 * i, amt=6.0 * wt,
                       cmt="tr"))
    return ev


class Scenario:
    def __init__(self, name, events=(), po=None, tend=730.0, label=""):
        self.name = name
        self.events = list(events)
        # po: dict drug -> (start_day, target_engagement)
        self.po = po or {}
        self.tend = tend
        self.label = label or name


# -----------------------------------------------------------------------------
# integrator: fixed-step RK4, vectorised over subjects.
#
# Two segments with different step sizes. While anthracycline / dexrazoxane are
# on board the plasma compartment is stiff (lambda ~52/d) and needs dt = 0.02 d;
# 60 days after the last such dose those compartments are empty, so they are
# pinned at zero and the remaining (slow) system is stepped at dt = 0.10 d.
# Exact per-subject maxima are tracked inside the loop, so reported peaks do not
# depend on the 1-day observation grid.
# -----------------------------------------------------------------------------
PEAKVARS = ("CHF", "CH", "CHM", "ROS", "P53", "TNI", "MITOD", "FUNC", "BNP")
# States pinned to zero in the slow segment: all small-molecule PK plus the
# fast nuclear pool (t1/2 = 0.03 d, hence numerically stiff at dt_slow and
# empty in any case 60 days after the last dose).
PKSTATES = (A1, A2, A3, ALIP, AM, AMP, ADEX, CHF)
CMTS = ("dox", "lipo", "dex", "tr")
DRUGS = ("ace", "bb", "sta", "arni", "sglt")
BOLCMT = dict(dox=A1, lipo=ALIP, dex=ADEX, tr=TR1)


def simulate(sc, p, nsub=1, etas=None, dt=0.02, dt_slow=0.10, obs_dt=1.0):
    n = nsub
    y = np.zeros((NSTATE, n))
    y[MYO] = 1.0
    y[LIP] = 1.0
    y[T2B] = 1.0
    y[BNP] = p["KBIN"] / p["KBOUT"]

    pp = dict(p)
    if etas is not None:
        for k, v in etas.items():
            pp[k] = np.asarray(v)

    nobs = int(round(sc.tend / obs_dt)) + 1
    keys = ("LVEF", "GLS", "TNI", "BNP", "MYO", "FUNC", "FIB", "ROS", "HYP",
            "P53", "MITOD", "CH", "CHF", "CHM")
    out = {k: np.zeros((nobs, n)) for k in keys}
    peak = {k: np.zeros(n) for k in PEAKVARS}
    tobs = np.arange(nobs) * obs_dt
    d0 = derived(y, pp)
    for k in keys:
        out[k][0] = d0[k]

    # last time a stiff (small-molecule) drug is administered
    t_stiff = 0.0
    for e in sc.events:
        if e["cmt"] in ("dox", "lipo", "dex"):
            t_stiff = max(t_stiff, e.get("t1", e["t0"]))
    t_switch = min(sc.tend, np.ceil(t_stiff + 60.0))

    segments = [(0.0, t_switch, dt)]
    if t_switch < sc.tend:
        segments.append((t_switch, sc.tend, dt_slow))

    for (tA, tB, h) in segments:
        nstep = int(round((tB - tA) / h))
        if nstep <= 0:
            continue
        if h != dt:                      # entering the slow segment
            for i in PKSTATES:
                y[i] = 0.0
        tg = tA + np.arange(nstep) * h
        grid = np.concatenate([tg, tg + h / 2, tg + h])
        RT = {c: np.zeros(3 * nstep) for c in CMTS}
        bol = {}
        for e in sc.events:
            if e["kind"] == "inf":
                m = (grid >= e["t0"] - 1e-9) & (grid < e["t1"] - 1e-9)
                RT[e["cmt"]][m] += e["rate"]
            elif tA - 1e-9 <= e["t0"] < tB - 1e-9:
                bol.setdefault(int(round((e["t0"] - tA) / h)), []).append(
                    (BOLCMT[e["cmt"]], e["amt"]))
        GT = {k: np.zeros(3 * nstep) for k in DRUGS}
        for k, (t0, lvl) in sc.po.items():
            GT[k][grid >= t0 - 1e-9] = lvl
        obs_step = int(round(obs_dt / h))

        for st in range(nstep):
            t = tA + st * h
            if st in bol:
                for ist, amt in bol[st]:
                    y[ist] = y[ist] + amt
            r1 = {c: RT[c][st] for c in CMTS}
            r2 = {c: RT[c][nstep + st] for c in CMTS}
            r3 = {c: RT[c][2 * nstep + st] for c in CMTS}
            g1 = {k: GT[k][st] for k in DRUGS}
            g2 = {k: GT[k][nstep + st] for k in DRUGS}
            g3 = {k: GT[k][2 * nstep + st] for k in DRUGS}
            k1 = rhs(t, y, pp, r1, g1)
            k2 = rhs(t + h / 2, y + h / 2 * k1, pp, r2, g2)
            k3 = rhs(t + h / 2, y + h / 2 * k2, pp, r2, g2)
            k4 = rhs(t + h, y + h * k3, pp, r3, g3)
            y = y + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            np.clip(y, 0.0, None, out=y)
            y[MYO] = np.minimum(y[MYO], 1.0)
            y[FUNC] = np.minimum(y[FUNC], 0.98)
            y[FIB] = np.minimum(y[FIB], 0.95)
            y[MITOD] = np.minimum(y[MITOD], 0.98)
            y[T2B] = np.minimum(y[T2B], 1.0)
            if h != dt:
                for i in PKSTATES:
                    y[i] = 0.0
            for k in PEAKVARS:
                idx = globals()[k]
                v = y[idx] + (pp["TNI_BASE"] if k == "TNI" else 0.0)
                np.maximum(peak[k], v, out=peak[k])
            if (st + 1) % obs_step == 0:
                i = int(round((t + h) / obs_dt))
                if i < nobs:
                    dd = derived(y, pp)
                    for k in keys:
                        out[k][i] = dd[k]
    out["t"] = tobs
    out["p"] = pp
    out["peak"] = peak
    return out


# -----------------------------------------------------------------------------
# population sampling
# -----------------------------------------------------------------------------
def sample_pop(n, seed=20260728, risk="standard"):
    rng = np.random.default_rng(seed)
    e = dict(
        KDROS=P["KDROS"] * np.exp(rng.normal(0, 0.55, n)),
        TOXN50=P["TOXN50"] * np.exp(rng.normal(0, 0.30, n)),
        KPIN=P["KPIN"] * np.exp(rng.normal(0, 0.35, n)),
        KPINF=P["KPINF"] * np.exp(rng.normal(0, 0.30, n)),
        KFIN=P["KFIN"] * np.exp(rng.normal(0, 0.40, n)),
        KFIBIN=P["KFIBIN"] * np.exp(rng.normal(0, 0.40, n)),
        EF0=np.clip(P["EF0"] + rng.normal(0, 4.0, n), 52.0, 72.0),
        CL=P["CL"] * np.exp(rng.normal(0, 0.28, n)),
    )
    if risk == "high":          # age >=70, hypertension, prior mediastinal RT
        e["KDROS"] = e["KDROS"] * 1.9
        e["KFIN"] = e["KFIN"] * 1.5
        e["KFIBIN"] = e["KFIBIN"] * 1.6
        e["EF0"] = np.clip(e["EF0"] - 4.0, 45.0, 72.0)
        e["KREG"] = np.full(n, P["KREG"] * 0.4)
    return e


def ctrcd_flags(out, upto=None):
    """CTRCD = LVEF fall >=10 points AND absolute LVEF <50%."""
    ef = out["LVEF"]
    if upto is not None:
        k = int(upto / (out["t"][1] - out["t"][0])) + 1
        ef = ef[:k]
    ef0 = ef[0]
    drop = ef0[None, :] - ef
    ctr = ((drop >= 10.0) & (ef < 50.0)).any(axis=0)
    hf = (ef < 40.0).any(axis=0)
    return ctr, hf


def nadir(out):
    return out["LVEF"].min(axis=0)


# =============================================================================
# scenario library
# =============================================================================
def scenarios():
    S = {}
    S["DOX240"] = Scenario("DOX240", dox_reg(60, 4, 21), tend=730,
                           label="doxorubicin 60 mg/m2 q3w x4 (240 mg/m2)")
    S["DOX360"] = Scenario("DOX360", dox_reg(60, 6, 21), tend=730,
                           label="doxorubicin 60 mg/m2 q3w x6 (360 mg/m2)")
    S["DOX480"] = Scenario("DOX480", dox_reg(60, 8, 21), tend=730,
                           label="doxorubicin 60 mg/m2 q3w x8 (480 mg/m2)")
    S["DOX480_INF72"] = Scenario(
        "DOX480_INF72", dox_reg(60, 8, 21, tinf=3.0), tend=730,
        label="480 mg/m2 given as 72-h continuous infusions")
    S["DOX480_WEEKLY"] = Scenario(
        "DOX480_WEEKLY", dox_reg(20, 24, 7), tend=730,
        label="480 mg/m2 as 20 mg/m2 weekly x24")
    S["DOX480_PLD"] = Scenario(
        "DOX480_PLD", dox_reg(40, 12, 28, lipo=True), tend=730,
        label="pegylated liposomal doxorubicin 40 mg/m2 q4w x12 (480 mg/m2)")
    S["DOX480_DEX"] = Scenario(
        "DOX480_DEX", dox_reg(60, 8, 21) + dex_reg(60, 8, 21), tend=730,
        label="480 mg/m2 + dexrazoxane 10:1")
    S["DOX480_STA"] = Scenario(
        "DOX480_STA", dox_reg(60, 8, 21), po=dict(sta=(0.0, 0.75)), tend=730,
        label="480 mg/m2 + atorvastatin 40 mg from day 0")
    S["DOX480_ACEBB"] = Scenario(
        "DOX480_ACEBB", dox_reg(60, 8, 21),
        po=dict(ace=(0.0, 0.80), bb=(0.0, 0.80)), tend=730,
        label="480 mg/m2 + enalapril 20 mg + carvedilol 25 mg bid from day 0")
    S["DOX480_ALL"] = Scenario(
        "DOX480_ALL", dox_reg(60, 8, 21) + dex_reg(60, 8, 21),
        po=dict(sta=(0.0, 0.75), ace=(0.0, 0.80), bb=(0.0, 0.80)), tend=730,
        label="480 mg/m2 + dexrazoxane + statin + ACEi/BB")
    S["DOX240_TSEQ"] = Scenario(
        "DOX240_TSEQ", dox_reg(60, 4, 21) + tras_reg(84.0), tend=730,
        label="240 mg/m2 then trastuzumab 1 year (sequential, day 84)")
    S["DOX240_TCONC"] = Scenario(
        "DOX240_TCONC", dox_reg(60, 4, 21) + tras_reg(0.0), tend=730,
        label="240 mg/m2 with concurrent trastuzumab 1 year")
    S["TRAS_ALONE"] = Scenario(
        "TRAS_ALONE", tras_reg(0.0), tend=730,
        label="trastuzumab 1 year, no anthracycline")
    return S


# =============================================================================
# analyses
# =============================================================================
def fmt(x, n=2):
    return f"{x:.{n}f}"


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


def a_index_patient(S):
    hdr("A0. INDEX PATIENT (no random effects) — 24-month trajectories")
    print(f"{'scenario':<16} {'EF 0':>6} {'EF 6m':>7} {'EF12m':>7} {'EF24m':>7} "
          f"{'nadir':>6} {'|GLS|12m':>8} {'TnI pk':>7} {'d(pk)':>6} "
          f"{'BNP12m':>7} {'MYO24m':>7} {'FUNC12m':>8} {'FIB24m':>7}")
    res = {}
    for k in ["DOX240", "DOX360", "DOX480", "DOX480_INF72", "DOX480_WEEKLY",
              "DOX480_PLD", "DOX480_DEX", "DOX480_STA", "DOX480_ACEBB",
              "DOX480_ALL", "DOX240_TSEQ", "DOX240_TCONC", "TRAS_ALONE"]:
        o = simulate(S[k], P, nsub=1)
        res[k] = o
        ef = o["LVEF"][:, 0]
        tni = o["TNI"][:, 0]
        tnipk = o["peak"]["TNI"][0]
        print(f"{k:<16} {fmt(ef[0],1):>6} {fmt(ef[182],1):>7} {fmt(ef[365],1):>7} "
              f"{fmt(ef[730],1):>7} {fmt(ef.min(),1):>6} "
              f"{fmt(o['GLS'][365,0],1):>8} {fmt(tnipk,1):>7} "
              f"{int(np.argmax(tni)):>6} {fmt(o['BNP'][365,0],0):>7} "
              f"{fmt(o['MYO'][730,0],3):>7} {fmt(o['FUNC'][365,0],3):>8} "
              f"{fmt(o['FIB'][730,0],3):>7}")
    return res


def a1_dose_response(S):
    hdr("A1. CUMULATIVE-DOSE RESPONSE (n=800 virtual subjects, 60 mg/m2 q3w)")
    n = 800
    etas = sample_pop(n)
    print(f"{'cum dose':>9} {'cycles':>7} {'CTRCD %':>8} {'LVEF<40 %':>10} "
          f"{'mean dEF12m':>12} {'mean nadir':>11}")
    rows = []
    for ncyc in (2, 4, 5, 6, 7, 8, 10):
        sc = Scenario("d", dox_reg(60, ncyc, 21), tend=730)
        o = simulate(sc, P, nsub=n, etas=etas)
        ctr, hf = ctrcd_flags(o)
        dEF = o["LVEF"][0] - o["LVEF"][365]
        rows.append((60 * ncyc, ncyc, 100 * ctr.mean(), 100 * hf.mean(),
                     dEF.mean(), nadir(o).mean()))
        print(f"{60*ncyc:>9} {ncyc:>7} {fmt(100*ctr.mean(),1):>8} "
              f"{fmt(100*hf.mean(),1):>10} {fmt(dEF.mean(),2):>12} "
              f"{fmt(nadir(o).mean(),1):>11}")
    print("\n  Von Hoff 1979 clinical-HF anchors: 300 mg/m2 ~2%, 400 ~3-5%,"
          " 550 ~7-18%.")
    return rows


def a2_leadtime(S):
    hdr("A2. WHAT MOVES FIRST — biomarker/strain/EF lead times (DOX480, n=600)")
    n = 600
    etas = sample_pop(n)
    o = simulate(S["DOX480"], P, nsub=n, etas=etas)
    t = o["t"]
    ef, gls, tni = o["LVEF"], o["GLS"], o["TNI"]

    def first_day(mask):
        idx = np.where(mask.any(axis=0), mask.argmax(axis=0), -1)
        return np.where(idx >= 0, t[np.maximum(idx, 0)], np.nan)

    d_tni = first_day(tni > 14.0)                       # hs-cTnI > 14 ng/L
    d_gls = first_day((gls[0] - gls) / gls[0] >= 0.15)   # relative GLS fall 15%
    d_ef = first_day((ef[0] - ef >= 10.0) & (ef < 50.0))
    ctr, hf = ctrcd_flags(o)
    print(f"  subjects with CTRCD by 24 mo : {100*ctr.mean():.1f}%")
    for nm, dd in (("hs-cTnI > 14 ng/L", d_tni), ("GLS relative fall >15%", d_gls),
                   ("CTRCD (EF)", d_ef)):
        v = dd[~np.isnan(dd)]
        print(f"  {nm:<24} reached by {100*len(v)/n:5.1f}% of subjects, "
              f"median day {np.median(v) if len(v) else float('nan'):.0f}")
    both = (~np.isnan(d_gls)) & (~np.isnan(d_ef))
    print(f"  GLS -> CTRCD lead time (subjects with both, n={both.sum()}): "
          f"median {np.median(d_ef[both]-d_gls[both]):.0f} days")
    both2 = (~np.isnan(d_tni)) & (~np.isnan(d_ef))
    print(f"  cTnI -> CTRCD lead time (n={both2.sum()}): "
          f"median {np.median(d_ef[both2]-d_tni[both2]):.0f} days")
    # positive/negative predictive value of an early troponin
    tni_pos = tni[:170].max(axis=0) > 40.0
    if tni_pos.sum() and (~tni_pos).sum():
        print(f"  cTnI>40 ng/L during chemo (n={tni_pos.sum()}): "
              f"CTRCD in {100*ctr[tni_pos].mean():.1f}% vs "
              f"{100*ctr[~tni_pos].mean():.1f}% if cTnI stayed <=40 "
              f"(NPV {100*(1-ctr[~tni_pos].mean()):.1f}%)")
    return o


def a3_exposure_metric(S):
    hdr("A3. WHICH EXPOSURE METRIC DRIVES INJURY (480 mg/m2, five schedules)")
    print(f"{'schedule':<16} {'peak CHF':>9} {'peak CH':>8} {'AUC_CH':>9} "
          f"{'peak P53':>9} {'peak ROS':>9} {'MYO24m':>8} {'EF12m':>7}")
    for k in ("DOX480", "DOX480_INF72", "DOX480_WEEKLY", "DOX480_PLD",
              "DOX480_DEX"):
        o = simulate(S[k], P, nsub=1)
        print(f"{k:<16} {fmt(o['peak']['CHF'][0],2):>9} {fmt(o['peak']['CH'][0],3):>8} "
              f"{fmt(np.trapezoid(o['CH'][:,0], o['t']),1):>9} "
              f"{fmt(o['peak']['P53'][0],3):>9} {fmt(o['peak']['ROS'][0],3):>9} "
              f"{fmt(o['MYO'][730,0],3):>8} {fmt(o['LVEF'][365,0],1):>7}")
    print("  (identical cumulative dose in every row -> cumulative dose is not"
          " the exposure metric)")


def a4_dex_sparing(S):
    hdr("A4. DEXRAZOXANE DOSE-SPARING EQUIVALENCE (index patient)")
    ref = simulate(Scenario("r", dox_reg(60, 4, 21), tend=730), P, nsub=1)
    target = ref["LVEF"][365, 0]
    print(f"  reference: 240 mg/m2 bolus, no protection -> EF12m = {target:.2f}")
    best = None
    for ncyc in range(4, 15):
        sc = Scenario("d", dox_reg(60, ncyc, 21) + dex_reg(60, ncyc, 21),
                      tend=730)
        o = simulate(sc, P, nsub=1)
        ef = o["LVEF"][365, 0]
        print(f"    + dexrazoxane, {60*ncyc:>3} mg/m2 -> EF12m = {ef:.2f}"
              f"{'   <-- iso-effective' if best is None and ef <= target else ''}")
        if best is None and ef <= target:
            best = 60 * ncyc
    if best:
        print(f"  iso-cardiotoxic dose with dexrazoxane ~ {best} mg/m2 "
              f"({best/240:.1f}x the unprotected 240 mg/m2)")


def a5_reversibility_window(S):
    hdr("A5. THE REVERSIBILITY WINDOW — when HF therapy is started")
    n = 500
    etas = sample_pop(n)
    base = simulate(S["DOX480"], P, nsub=n, etas=etas)
    ctr, _ = ctrcd_flags(base, upto=365)
    print(f"  CTRCD by 12 mo in {100*ctr.mean():.1f}% of {n} subjects; "
          f"analysis restricted to those {ctr.sum()} subjects")
    print(f"{'start day':>10} {'EF at start':>12} {'EF +12mo':>10} "
          f"{'delta':>7} {'recovered %':>12} {'FIB at start':>13}")
    for t0 in (60, 90, 120, 180, 270, 365, 540):
        sc = Scenario("d", dox_reg(60, 8, 21),
                      po=dict(ace=(float(t0), 0.80), bb=(float(t0), 0.80)),
                      tend=float(t0) + 365.0)
        o = simulate(sc, P, nsub=n, etas=etas)
        i0 = t0
        i1 = t0 + 365
        ef0 = o["LVEF"][i0][ctr]
        ef1 = o["LVEF"][i1][ctr]
        base_ef0 = o["LVEF"][0][ctr]
        rec = ((ef1 >= 50.0) | (base_ef0 - ef1 < 5.0)).mean()
        print(f"{t0:>10} {ef0.mean():>12.2f} {ef1.mean():>10.2f} "
              f"{ef1.mean()-ef0.mean():>7.2f} {100*rec:>12.1f} "
              f"{o['FIB'][i0][ctr].mean():>13.3f}")
    print("  Cardinale 2010 anchor: LVEF recovery 64% if enalapril+carvedilol"
          " started <2 mo after detection, ~0% if >6 mo.")


def a6_trastuzumab_interaction(S):
    hdr("A6. ANTHRACYCLINE x TRASTUZUMAB — additive or synergistic?")
    n = 600
    etas = sample_pop(n)
    out = {}
    for k in ("DOX240", "TRAS_ALONE", "DOX240_TSEQ", "DOX240_TCONC"):
        o = simulate(S[k], P, nsub=n, etas=etas)
        ctr, hf = ctrcd_flags(o)
        dEF = (o["LVEF"][0] - o["LVEF"][:456].min(axis=0)).mean()
        # recovery after trastuzumab washout: EF at 24 mo vs nadir
        rec = (o["LVEF"][730] - o["LVEF"][:456].min(axis=0)).mean()
        out[k] = dict(ctr=100 * ctr.mean(), hf=100 * hf.mean(), dEF=dEF,
                      rec=rec, myo=o["MYO"][730].mean())
        print(f"  {k:<14} CTRCD {out[k]['ctr']:5.1f}%  LVEF<40 {out[k]['hf']:4.1f}%"
              f"  max EF fall {dEF:5.2f}  EF regained by 24 mo {rec:5.2f}"
              f"  myocytes left {out[k]['myo']:.3f}")
    a = out["DOX240"]["dEF"]
    b = out["TRAS_ALONE"]["dEF"]
    for k in ("DOX240_TSEQ", "DOX240_TCONC"):
        obs = out[k]["dEF"]
        print(f"  {k}: additive prediction {a+b:.2f}, observed {obs:.2f} -> "
              f"{'SYNERGY' if obs > a+b+0.2 else 'sub-additive/additive'} "
              f"(excess {obs-(a+b):+.2f} EF points)")
    print("  Note EF regained by month 24: trastuzumab-dominated injury is the"
          " reversible kind; anthracycline-dominated injury is not.")


def a7_cbr1(S):
    hdr("A7. CBR1/AKR1C3 METABOLISER STATUS (doxorubicinol formation fraction)")
    n = 500
    print(f"{'FM':>6} {'peak CHM':>9} {'peak ROS':>9} {'EF12m (index)':>14} "
          f"{'CTRCD % (pop)':>14}")
    for fm in (0.12, 0.18, 0.25, 0.34, 0.45):
        p = dict(P); p["FM"] = fm
        o1 = simulate(S["DOX480"], p, nsub=1)
        e = sample_pop(n)
        op = simulate(S["DOX480"], p, nsub=n, etas=e)
        ctr, _ = ctrcd_flags(op)
        print(f"{fm:>6} {o1['CHM'].max():>9.3f} {o1['ROS'].max():>9.3f} "
              f"{o1['LVEF'][365,0]:>14.2f} {100*ctr.mean():>14.1f}")


def a8_population_table(S):
    hdr("A8. VIRTUAL-POPULATION INCIDENCE TABLE (n=800 each, 24 months)")
    n = 800
    etas = sample_pop(n)
    print(f"{'scenario':<16} {'CTRCD %':>8} {'RR vs ref':>10} {'LVEF<40 %':>10} "
          f"{'mean dEF12m':>12} {'myocytes 24m':>13}")
    ref = None
    for k in ("DOX480", "DOX480_INF72", "DOX480_WEEKLY", "DOX480_PLD",
              "DOX480_DEX", "DOX480_STA", "DOX480_ACEBB", "DOX480_ALL",
              "DOX240", "DOX240_TSEQ", "DOX240_TCONC"):
        o = simulate(S[k], P, nsub=n, etas=etas)
        ctr, hf = ctrcd_flags(o)
        if ref is None:
            ref = ctr.mean()
        dEF = (o["LVEF"][0] - o["LVEF"][365]).mean()
        print(f"{k:<16} {100*ctr.mean():>8.1f} {ctr.mean()/ref:>10.2f} "
              f"{100*hf.mean():>10.1f} {dEF:>12.2f} {o['MYO'][730].mean():>13.3f}")
    print("  Cochrane anchors — dexrazoxane RR(HF) 0.29; continuous infusion"
          " RR 0.27; liposomal RR 0.20. STOP-CA: statin 9% vs 22%.")


def a9_prevention_vs_rescue(S):
    hdr("A9. UPSTREAM vs DOWNSTREAM — same drugs given during vs after chemo")
    n = 600
    etas = sample_pop(n)
    print(f"{'strategy':<34} {'CTRCD %':>8} {'EF12m':>7} {'EF24m':>7}")
    combos = [
        ("no protection", {}),
        ("statin during chemo (day 0)", dict(sta=(0.0, 0.75))),
        ("statin after chemo (day 170)", dict(sta=(170.0, 0.75))),
        ("ACEi+BB during chemo (day 0)", dict(ace=(0.0, 0.8), bb=(0.0, 0.8))),
        ("ACEi+BB after chemo (day 170)", dict(ace=(170.0, 0.8), bb=(170.0, 0.8))),
        ("ARNI after chemo (day 170)", dict(arni=(170.0, 0.9))),
        ("SGLT2i during chemo (day 0)", dict(sglt=(0.0, 0.70))),
    ]
    for nm, po in combos:
        sc = Scenario("d", dox_reg(60, 8, 21), po=po, tend=730)
        o = simulate(sc, P, nsub=n, etas=etas)
        ctr, _ = ctrcd_flags(o)
        print(f"{nm:<34} {100*ctr.mean():>8.1f} {o['LVEF'][365].mean():>7.2f} "
              f"{o['LVEF'][730].mean():>7.2f}")
    print("  The statin arm only works while drug is on board (it removes ROS"
          " production); the ACEi/BB/ARNI arms work afterwards (they act on"
          " the reversible functional and fibrotic terms).")


def a10_decomposition(S):
    hdr("A10. WHAT IS THE EF DEFICIT MADE OF (DOX480 index patient)")
    o = simulate(S["DOX480"], P, nsub=1)
    p = P
    print(f"{'day':>5} {'LVEF':>6} {'MYO':>6} {'FUNC':>6} {'FIB':>6} {'HYP':>6} "
          f"{'EF if FUNC=0':>13} {'EF if MYO=1':>12} {'EF if HYP=0':>12}")
    for day in (0, 60, 120, 170, 240, 365, 540, 730):
        MYOv = o["MYO"][day, 0]; FUNCv = o["FUNC"][day, 0]
        FIBv = o["FIB"][day, 0]; HYPv = o["HYP"][day, 0]

        def ef(myo, func, fib, hyp):
            c = myo * (1 - func) * (1 + hyp) / (1 + p["WFIBC"] * fib)
            return p["EF0"] * c ** p["EFEXP"]
        print(f"{day:>5} {ef(MYOv,FUNCv,FIBv,HYPv):>6.2f} {MYOv:>6.3f} "
              f"{FUNCv:>6.3f} {FIBv:>6.3f} {HYPv:>6.3f} "
              f"{ef(MYOv,0,FIBv,HYPv):>13.2f} {ef(1,FUNCv,FIBv,HYPv):>12.2f} "
              f"{ef(MYOv,FUNCv,FIBv,0):>12.2f}")
    print("  'EF if HYP=0' is the EF the patient would have without"
          " compensatory hypertrophy = the size of the mask.")
    # masking: cumulative myocyte loss at the moment EF first drops 5 points
    ef_tr = o["LVEF"][:, 0]
    i5 = np.argmax(ef_tr[0] - ef_tr >= 5.0)
    print(f"  EF has fallen 5 points only at day {i5}, by which time "
          f"{100*(1-o['MYO'][i5,0]):.1f}% of myocytes are already gone "
          f"and |GLS| is down {100*(o['GLS'][0,0]-o['GLS'][i5,0])/o['GLS'][0,0]:.1f}%.")


def a11_highrisk(S):
    hdr("A11. HIGH-RISK PHENOTYPE (age>=70 + hypertension + prior chest RT)")
    n = 700
    for risk in ("standard", "high"):
        e = sample_pop(n, risk=risk)
        for k in ("DOX240", "DOX480", "DOX240_DEXPROT"):
            if k == "DOX240_DEXPROT":
                sc = Scenario("d", dox_reg(60, 4, 21) + dex_reg(60, 4, 21),
                              tend=730)
            else:
                sc = S[k]
            o = simulate(sc, P, nsub=n, etas=e)
            ctr, hf = ctrcd_flags(o)
            print(f"  {risk:<9} {k:<16} CTRCD {100*ctr.mean():5.1f}%   "
                  f"LVEF<40 {100*hf.mean():4.1f}%   "
                  f"EF12m {o['LVEF'][365].mean():5.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    S = scenarios()
    a_index_patient(S)
    if args.quick:
        return
    a1_dose_response(S)
    a2_leadtime(S)
    a3_exposure_metric(S)
    a4_dex_sparing(S)
    a5_reversibility_window(S)
    a6_trastuzumab_interaction(S)
    a7_cbr1(S)
    a8_population_table(S)
    a9_prevention_vs_rescue(S)
    a10_decomposition(S)
    a11_highrisk(S)


if __name__ == "__main__":
    main()
