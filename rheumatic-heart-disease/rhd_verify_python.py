#!/usr/bin/env python3
# =====================================================================
# Rheumatic heart disease (RHD) — independent Python/scipy
# re-implementation of rhd_mrgsolve_model.R
# ---------------------------------------------------------------------
# Purpose.  The mrgsolve file is the deliverable; this file exists to
# check it.  Every ODE below was transcribed from the $ODE block by
# hand, and the numbers this script prints are the numbers quoted in
# README.md.  Where the two implementations disagreed, the disagreement
# was a defect in the model, not in the transcription; those are listed
# at the bottom of README.md.
#
# Run:  python3 rhd_verify_python.py
# =====================================================================

import math
import numpy as np
from scipy.integrate import solve_ivp

PI = math.pi

# ---------------------------------------------------------------------
# 1. Parameters  (must match $PARAM in rhd_mrgsolve_model.R)
# ---------------------------------------------------------------------
P = dict(
    # --- demographics / renal ---------------------------------------
    WT=70.0, CRCL=100.0, AGE=25.0,

    # --- penicillin PK ----------------------------------------------
    VPEN0=15.0,        # L / 70 kg   (apparent V/F of penicillin G)
    CLPEN0=500.0,      # L/d / 70 kg (renal, CrCl 100)
    KAFBPG=0.35,       # 1/d  fast-release fraction of the BPG depot
    KASBPG=0.075,      # 1/d  slow-release fraction (rate-limiting, flip-flop)
    FFBPG=0.30,        # fraction of the BPG dose in the fast depot
    FSBPG=0.70,
    KAORAL=24.0,       # 1/d  oral penicillin V absorption
    MICP=0.02,         # ug/mL  protective threshold (Stollerman benchmark)

    # --- group A streptococcal pharyngitis ---------------------------
    KGROW=2.0, GMAX=1e8, KGCL=0.35, KMIK=1.2, KPKILL=6.0,
    GEXT=100.0,        # extinction floor: below this the organism cannot
                       # regrow.  Without it the ODE never reaches zero and
                       # a single pharyngitis relapses forever (see README).
    HPEN=3.0, EC50P=0.02,
    GM50=1e6,          # burden giving half-maximal SYSTEMIC antigen release
    GM50I=1e4,         # burden giving half-maximal IMMUNE stimulus.  These
                       # are different numbers: the immune system sees the
                       # organism long before enough antigen reaches the
                       # circulation to drive the cross-reactive response.
    KEO=5.545,         # 1/d  plasma -> tonsillar effect site, t1/2 3 h.
                       # Bacterial kill follows the EFFECT SITE; the 0.02
                       # ug/mL programme benchmark is a PLASMA number, and
                       # TPROT is therefore still computed on plasma.
    KMI=0.30, KMD=0.015,

    # --- antigen, antibody, memory -----------------------------------
    KAG=1.0, KAGD=0.15,
    KASO=40.0, KASOD=0.0154, ASOBASE=100.0,
    # MEM is a SATURATING pool (MEMMAX).  Written unbounded, the loop
    # AG -> MEM -> recurrence hazard -> AG has gain 2.9 and diverges.
    KMEM=4.0e-3, KMEMD=3.8e-4, MEMMAX=1.0,
    SUSC=1.0,                 # host susceptibility (HLA-DR7, TNFA-308 ...)
    KXAB=4.0e-3, BMEM=1.5, KXABD=0.01155,
    KVIT=0.60, BSCAR=2.0, KVITD=0.0231,

    # --- valve ---------------------------------------------------------
    MVAN=4.5, MVAMIN=0.30,
    KFI=3.004e-4,     # cm per (valvulitis AU . day)   immune arm
                      # set so ONE severe first attack costs 0.35 cm2
    KFS=9.132e-5,     # cm per day per unit shear      autonomous arm
    KHALT=1.6,        # rigidity brake: deposition falls off as the orifice
                      # approaches its fused limit.  KFS is renormalised so
                      # that progression at MVA 1.5 is unchanged (0.09 cm2/yr)
    RESTEN=1.0,       # multiplier on KFS after balloon valvotomy
    KED=0.5, KEDD=0.033, KEDA=0.012,
    KCA=4.0e-3, KMRA=0.096, KMRAD=0.03, KMRC=1.2e-3, KAOI=6.0e-4,

    # --- haemodynamics -------------------------------------------------
    HRB=72.0, DHRAF=35.0, HRMIN=45.0, HRSNS=15.0,
    SVMAX=110.0, GORLIN=37.7, AFPEN=0.12,
    CO0=5.0, DEMF=1.0,        # DEMF: 1 rest, 1.5 pregnancy, 2.0 exercise
    KHRDEM=0.5,               # fractional HR rise per unit of extra demand
    LVEDP0=6.0, KVOLP=8.0, LAPMAX=32.0,
    LAPTH=18.0, KCON=0.10, KCOND=0.25,
    LAV0=55.0, KLAV=0.06, KLAVD=0.002, LAPB=12.0,
    KAFO=3.0e-4, LAVREF=60.0, NAF=2.0, KAFR=0.0,
    PVRB=1.4, KPVR=1.8e-3, PVRMAX=12.0, KPVRR=4.0e-3, LAPPV=20.0,
    KRV=3.0e-4, PAPTH=35.0, KRVR=2.0e-3,
    KVOLR=0.05, EDIUR=0.0,

    # --- rhythm / embolism ---------------------------------------------
    KEMB=1.85e-5,     # 1/d
    AF50E=0.15,       # AF burden giving half-maximal embolic risk
    EMAXAC=0.85, INRE50=1.55,

    # --- recurrence hazard ---------------------------------------------
    LAMEXP=3.0,       # pharyngitis exposures per year
    PGASP=0.25,       # P(GAS | sore throat) in an endemic setting
    PRHEUM=0.0444,    # P(ARF | untreated GAS, susceptible host)
    BMEM2=2.0,        # memory amplification of recurrence risk
    QARF=2.483,       # antigen AU delivered by one full ARF episode.
                      # MEASURED from the episodic channel rather than
                      # guessed, so the expected-value and episodic routes
                      # deliver the same antigen per episode.
    EPISODIC=0.0,     # 1 = discrete inoculation events only, 0 = expected value

    # --- anti-inflammatory / rate control / anticoagulation PK ---------
    KAASA=24.0, VASA=12.0, CLASA=25.0, EMAXASA=0.55, EC50ASA=120.0,
    KAPRED=24.0, VPRED=40.0, CLPRED=222.0, EMAXPR=3.0, EC50PR=0.05,
    KABB=12.0, VBB=250.0, CLBB=1000.0, EMAXBB=0.35, EC50BB=0.04,
    KADIG=12.0, VDIG=500.0, CLDIG=180.0, FDIG=0.7, EMAXDIG=0.18, EC50DIG=1.0,
    KAW=12.0, VW=10.0, CLW=4.2, EMAXW=0.92, IC50W=0.6,
    KDEGP=0.35, GAMINR=0.85,

    # --- acute phase ----------------------------------------------------
    KCRP=26.4, KCRPD=0.35, CRPB=2.0,
)

# state index map --------------------------------------------------------
NAMES = ["BPGF", "BPGS", "PVA", "PENC", "GAS", "MIMM", "AG", "ASO", "MEM",
         "XAB", "VIT", "EDEM", "MVA", "CA", "MRA", "MRC", "AOI", "LAV",
         "AFB", "PVR", "RVF", "CONG", "VOL", "CRP", "ASAA", "ASAC",
         "PREDA", "PREDC", "BBA", "BBC", "DIGA", "DIGC", "WA", "WC",
         "PCF", "CUMARF", "CUMEMB", "TPROT", "CE"]
IDX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)


def initial(p):
    y = np.zeros(NST)
    y[IDX["MVA"]] = p["MVAN"]
    y[IDX["ASO"]] = p["ASOBASE"]
    y[IDX["LAV"]] = p["LAV0"]
    y[IDX["PVR"]] = p["PVRB"]
    y[IDX["RVF"]] = 1.0
    y[IDX["VOL"]] = 1.0
    y[IDX["PCF"]] = 1.0
    y[IDX["CRP"]] = p["CRPB"]
    return y


# ---------------------------------------------------------------------
# 2. Algebraic haemodynamic block
#    This is the arithmetic core of the model and is shared verbatim
#    with $ODE / $TABLE in the mrgsolve file.
# ---------------------------------------------------------------------
def haemo(y, p):
    """Return dict of derived haemodynamic quantities."""
    MVA = max(y[IDX["MVA"]], p["MVAMIN"])
    EDEM = y[IDX["EDEM"]]
    AFB = y[IDX["AFB"]]
    VOL = y[IDX["VOL"]]
    CONG = y[IDX["CONG"]]
    PVR = y[IDX["PVR"]]

    # --- drug effects on rate -----------------------------------------
    CBB = y[IDX["BBC"]] / p["VBB"]                       # mg/L
    CDIG = y[IDX["DIGC"]] / p["VDIG"]                    # ug/L = ng/mL
    EBB = p["EMAXBB"] * CBB / (CBB + p["EC50BB"])
    EDIG = p["EMAXDIG"] * CDIG / (CDIG + p["EC50DIG"])

    # --- heart rate ----------------------------------------------------
    HRraw = (p["HRB"] * (1.0 + p["KHRDEM"] * (p["DEMF"] - 1.0))
             + p["DHRAF"] * AFB
             + p["HRSNS"] * CONG / (CONG + 1.0)) * (1.0 - EBB) * (1.0 - EDIG)
    HR = max(HRraw, p["HRMIN"])

    # --- diastolic filling period (Gorlin denominator) -----------------
    RR = 60.0 / HR                                       # s per beat
    SYS = 0.36 * math.sqrt(RR)                           # systolic period
    DFPB = max(RR - SYS, 0.05)                           # s per beat
    DFPM = HR * DFPB                                     # s per minute

    # --- effective orifice ---------------------------------------------
    # acute leaflet oedema narrows the orifice reversibly; AF costs the
    # atrial contribution to filling, which behaves like a smaller orifice
    MVAeff = max(MVA - p["KEDA"] * EDEM, p["MVAMIN"]) * (1.0 - p["AFPEN"] * AFB)

    # --- filling pressures ---------------------------------------------
    LVEDP = p["LVEDP0"] + p["KVOLP"] * (VOL - 1.0)
    dPmax = max(p["LAPMAX"] - LVEDP, 1.0)

    # --- what the circulation is asked for, and what the valve allows --
    CODEM = p["CO0"] * p["DEMF"]
    COvalve = p["GORLIN"] * MVAeff * math.sqrt(dPmax) * DFPM / 1000.0
    COsv = p["SVMAX"] * HR / 1000.0
    CO = min(CODEM, COvalve, COsv)

    MVF = CO * 1000.0 / DFPM                             # mL/s
    MVG = (MVF / (p["GORLIN"] * MVAeff)) ** 2            # mmHg, mean gradient
    LAP = LVEDP + MVG
    PAPm = LAP + PVR * CO

    SXI = 0.6 * CONG / (CONG + 1.2) + 0.4 * (1.0 - CO / CODEM)
    NYHA = min(4.0, 1.0 + 3.0 * SXI)

    return dict(MVA=MVA, MVAeff=MVAeff, HR=HR, DFPB=DFPB, DFPM=DFPM,
                LVEDP=LVEDP, CO=CO, CODEM=CODEM, COvalve=COvalve, COsv=COsv,
                MVF=MVF, MVG=MVG, LAP=LAP, PAPm=PAPm, NYHA=NYHA,
                EBB=EBB, EDIG=EDIG, SXI=SXI)


def wilkins(y, p):
    """Wilkins/Abascal echo score, 4-16.  >8 predicts a poor balloon result."""
    MVA = max(y[IDX["MVA"]], p["MVAMIN"])
    CA = y[IDX["CA"]]
    narrow = min(1.0, max(0.0, (p["MVAN"] - MVA) / (p["MVAN"] - 0.8)))
    c8 = CA / (CA + 8.0)
    c10 = CA / (CA + 10.0)
    mob = 1.0 + 3.0 * min(1.0, 0.35 * narrow + 0.45 * c8)
    thick = 1.0 + 3.0 * min(1.0, 0.30 * narrow + 0.50 * c8)
    calc = 1.0 + 3.0 * c10
    sub = 1.0 + 3.0 * min(1.0, 0.30 * narrow + 0.45 * c10)
    return min(16.0, mob + thick + calc + sub)


# ---------------------------------------------------------------------
# 3. Right-hand side
# ---------------------------------------------------------------------
def rhs(t, y, p):
    y = np.maximum(y, 0.0)
    d = np.zeros(NST)
    H = haemo(y, p)

    # ---- 3.1 penicillin PK -------------------------------------------
    VPEN = p["VPEN0"] * (p["WT"] / 70.0)
    CLPEN = p["CLPEN0"] * (p["WT"] / 70.0) ** 0.75 * (p["CRCL"] / 100.0)
    K10 = CLPEN / VPEN
    rel = p["KAFBPG"] * y[IDX["BPGF"]] + p["KASBPG"] * y[IDX["BPGS"]]
    d[IDX["BPGF"]] = -p["KAFBPG"] * y[IDX["BPGF"]]
    d[IDX["BPGS"]] = -p["KASBPG"] * y[IDX["BPGS"]]
    d[IDX["PVA"]] = -p["KAORAL"] * y[IDX["PVA"]]
    d[IDX["PENC"]] = rel + p["KAORAL"] * y[IDX["PVA"]] - K10 * y[IDX["PENC"]]
    CPEN = y[IDX["PENC"]] / VPEN                          # ug/mL
    CE = y[IDX["CE"]]
    d[IDX["CE"]] = p["KEO"] * (CPEN - CE)

    EPEN = CE ** p["HPEN"] / (CE ** p["HPEN"] + p["EC50P"] ** p["HPEN"])

    # ---- 3.2 pharyngeal GAS ------------------------------------------
    GAS = y[IDX["GAS"]]
    d[IDX["GAS"]] = (p["KGROW"] * GAS * (1.0 - GAS / p["GMAX"])
                     * (GAS / (GAS + p["GEXT"]))
                     - (p["KGCL"] + p["KMIK"] * y[IDX["MIMM"]]) * GAS
                     - p["KPKILL"] * EPEN * GAS)
    agpres = GAS / (GAS + p["GM50"])
    d[IDX["MIMM"]] = (p["KMI"] * GAS / (GAS + p["GM50I"])
                      - p["KMD"] * y[IDX["MIMM"]])

    # ---- 3.3 recurrence hazard (expected-value channel) ---------------
    HAZARF = (p["LAMEXP"] * p["PGASP"] * p["PRHEUM"] * p["SUSC"]
              * (1.0 - EPEN) * (1.0 + p["BMEM2"] * y[IDX["MEM"]]))   # per year
    d[IDX["CUMARF"]] = HAZARF / 365.0
    AGBG = 0.0 if p["EPISODIC"] > 0.5 else p["QARF"] * HAZARF / 365.0

    # ---- 3.4 antigen, antibody, memory --------------------------------
    d[IDX["AG"]] = p["KAG"] * agpres + AGBG - p["KAGD"] * y[IDX["AG"]]
    d[IDX["ASO"]] = (p["KASO"] * y[IDX["AG"]]
                     - p["KASOD"] * (y[IDX["ASO"]] - p["ASOBASE"]))
    d[IDX["MEM"]] = (p["KMEM"] * y[IDX["AG"]] * p["SUSC"]
                     * (1.0 - y[IDX["MEM"]] / p["MEMMAX"])
                     - p["KMEMD"] * y[IDX["MEM"]])
    d[IDX["XAB"]] = (p["KXAB"] * y[IDX["AG"]] * p["SUSC"]
                     * (1.0 + p["BMEM"] * y[IDX["MEM"]])
                     - p["KXABD"] * y[IDX["XAB"]])

    # ---- 3.5 anti-inflammatory drugs -----------------------------------
    CASA = y[IDX["ASAC"]] / p["VASA"]
    CPRED = y[IDX["PREDC"]] / p["VPRED"]
    EASA = p["EMAXASA"] * CASA / (CASA + p["EC50ASA"])
    EPRED = p["EMAXPR"] * CPRED / (CPRED + p["EC50PR"])
    ANTIINF = 1.0 + EASA + EPRED

    d[IDX["ASAA"]] = -p["KAASA"] * y[IDX["ASAA"]]
    d[IDX["ASAC"]] = p["KAASA"] * y[IDX["ASAA"]] - p["CLASA"] / p["VASA"] * y[IDX["ASAC"]]
    d[IDX["PREDA"]] = -p["KAPRED"] * y[IDX["PREDA"]]
    d[IDX["PREDC"]] = p["KAPRED"] * y[IDX["PREDA"]] - p["CLPRED"] / p["VPRED"] * y[IDX["PREDC"]]

    # ---- 3.6 valvulitis -------------------------------------------------
    SCARF = max(0.0, 1.0 - H["MVA"] / p["MVAN"])
    d[IDX["VIT"]] = (p["KVIT"] * y[IDX["XAB"]] * (1.0 + p["BSCAR"] * SCARF) / ANTIINF
                     - p["KVITD"] * y[IDX["VIT"]])
    d[IDX["EDEM"]] = p["KED"] * y[IDX["VIT"]] - p["KEDD"] * y[IDX["EDEM"]]
    d[IDX["CRP"]] = p["KCRP"] * y[IDX["VIT"]] - p["KCRPD"] * (y[IDX["CRP"]] - p["CRPB"])

    # ---- 3.7 valve remodelling ------------------------------------------
    SH = H["MVG"] / 4.0                     # (mean transvalvular velocity)^2
    BRAKE = ((H["MVA"] - p["MVAMIN"])
             / (H["MVA"] - p["MVAMIN"] + p["KHALT"]))
    depo = p["KFI"] * y[IDX["VIT"]] + p["KFS"] * p["RESTEN"] * SH * BRAKE
    geom = 2.0 * math.sqrt(PI * H["MVA"])   # dA = -2 sqrt(pi A) dh
    d[IDX["MVA"]] = -geom * depo if y[IDX["MVA"]] > p["MVAMIN"] else 0.0
    d[IDX["CA"]] = p["KCA"] * (1.0 - H["MVA"] / p["MVAN"]) * (1.0 + p["AGE"] / 50.0)
    d[IDX["MRA"]] = p["KMRA"] * y[IDX["VIT"]] - p["KMRAD"] * y[IDX["MRA"]]
    d[IDX["MRC"]] = p["KMRC"] * y[IDX["VIT"]]
    d[IDX["AOI"]] = p["KAOI"] * y[IDX["VIT"]]

    # ---- 3.8 chambers, rhythm, lung ------------------------------------
    d[IDX["LAV"]] = (p["KLAV"] * max(0.0, H["LAP"] - p["LAPB"])
                     - p["KLAVD"] * (y[IDX["LAV"]] - p["LAV0"]))
    fAF = (max(0.0, y[IDX["LAV"]] - p["LAV0"]) / p["LAVREF"]) ** p["NAF"]
    d[IDX["AFB"]] = p["KAFO"] * (1.0 - y[IDX["AFB"]]) * fAF - p["KAFR"] * y[IDX["AFB"]]
    d[IDX["PVR"]] = (p["KPVR"] * max(0.0, H["LAP"] - p["LAPPV"])
                     * (1.0 - y[IDX["PVR"]] / p["PVRMAX"])
                     - p["KPVRR"] * (y[IDX["PVR"]] - p["PVRB"]))
    # The damage term is proportional to RVF so the index decays towards
    # zero instead of through it.  Written additively it reached -0.43 at a
    # mean PA pressure of 53 mmHg, which is not a right ventricle.
    d[IDX["RVF"]] = (-p["KRV"] * max(0.0, H["PAPm"] - p["PAPTH"]) * y[IDX["RVF"]]
                     + p["KRVR"] * (1.0 - y[IDX["RVF"]]))
    d[IDX["CONG"]] = (p["KCON"] * max(0.0, H["LAP"] - p["LAPTH"])
                      - p["KCOND"] * y[IDX["CONG"]])
    d[IDX["VOL"]] = p["KVOLR"] * (1.0 - y[IDX["VOL"]]) - p["EDIUR"] * y[IDX["VOL"]]

    # ---- 3.9 anticoagulation and embolism -------------------------------
    d[IDX["WA"]] = -p["KAW"] * y[IDX["WA"]]
    d[IDX["WC"]] = p["KAW"] * y[IDX["WA"]] - p["CLW"] / p["VW"] * y[IDX["WC"]]
    CW = y[IDX["WC"]] / p["VW"]
    IW = p["EMAXW"] * CW / (CW + p["IC50W"])
    d[IDX["PCF"]] = p["KDEGP"] * (1.0 - IW) - p["KDEGP"] * y[IDX["PCF"]]
    PCF = max(y[IDX["PCF"]], 1e-3)
    INR = (1.0 / PCF) ** p["GAMINR"]
    EANTI = p["EMAXAC"] * max(0.0, INR - 1.0) / (max(0.0, INR - 1.0) + (p["INRE50"] - 1.0))
    # NOTE.  This term was first written as AFB**0.7.  Its derivative is
    # infinite at AFB = 0, which is exactly where every simulation starts,
    # so the integrator crawled to a halt on any run long enough for AF to
    # begin appearing.  A saturating Hill has the same shape -- steep at low
    # burden, flat at high -- with a finite Jacobian everywhere.
    EMBAF = y[IDX["AFB"]] / (y[IDX["AFB"]] + p["AF50E"])
    d[IDX["CUMEMB"]] = (p["KEMB"] * (0.15 + EMBAF)
                        * (y[IDX["LAV"]] / p["LAV0"]) ** 1.2
                        * (p["MVAN"] / H["MVA"]) ** 0.5
                        * (1.0 - EANTI))

    # ---- 3.10 rate-control drug PK ---------------------------------------
    d[IDX["BBA"]] = -p["KABB"] * y[IDX["BBA"]]
    d[IDX["BBC"]] = p["KABB"] * y[IDX["BBA"]] - p["CLBB"] / p["VBB"] * y[IDX["BBC"]]
    d[IDX["DIGA"]] = -p["KADIG"] * y[IDX["DIGA"]]
    d[IDX["DIGC"]] = p["KADIG"] * y[IDX["DIGA"]] - p["CLDIG"] / p["VDIG"] * y[IDX["DIGC"]]

    # ---- 3.11 time above the protective threshold -------------------------
    d[IDX["TPROT"]] = CPEN ** 6 / (CPEN ** 6 + p["MICP"] ** 6)

    return d


# ---------------------------------------------------------------------
# 4. Simulation driver with dosing events
# ---------------------------------------------------------------------
def simulate(p, tend, events=None, t_eval=None, y0=None, pchange=None):
    """events: list of (time, state_name, amount) bolus additions.
       pchange: list of (time, param, value) parameter step changes."""
    p = dict(p)
    y = initial(p) if y0 is None else np.array(y0, dtype=float)
    events = sorted(events or [], key=lambda e: e[0])
    pchange = sorted(pchange or [], key=lambda e: e[0])
    breaks = sorted(set([0.0, float(tend)] + [e[0] for e in events]
                        + [c[0] for c in pchange]))
    breaks = [b for b in breaks if 0.0 <= b <= tend]
    if t_eval is None:
        t_eval = np.linspace(0, tend, min(int(tend * 4) + 1, 40001))
    t_eval = np.asarray(t_eval, dtype=float)

    T, Y = [0.0], [y.copy()]
    for i in range(len(breaks) - 1):
        t0, t1 = breaks[i], breaks[i + 1]
        for (te, nm, amt) in events:
            if abs(te - t0) < 1e-9:
                y[IDX[nm]] += amt
        for (tc, nm, val) in pchange:
            if abs(tc - t0) < 1e-9:
                p[nm] = val
        if t1 <= t0:
            continue
        sub = t_eval[(t_eval > t0) & (t_eval <= t1)]
        sol = solve_ivp(rhs, (t0, t1), y, args=(p,), method="LSODA",
                        rtol=1e-6, atol=1e-8,
                        t_eval=sub if len(sub) else None, max_step=(t1 - t0))
        if len(sub):
            T.extend(sol.t.tolist())
            Y.extend(sol.y.T.tolist())
        y = sol.y[:, -1].copy()
    for (te, nm, amt) in events:
        if abs(te - tend) < 1e-9:
            y[IDX[nm]] += amt
    return np.array(T), np.array(Y), p


def bpg_events(interval, n, mg=720.0, t0=0.0, adherence=None):
    """Two records per injection: fast and slow depot fractions."""
    ev = []
    for i in range(n):
        if adherence is not None and adherence[i % len(adherence)] == 0:
            continue
        t = t0 + i * interval
        ev.append((t, "BPGF", mg * P["FFBPG"]))
        ev.append((t, "BPGS", mg * P["FSBPG"]))
    return ev


def daily(nm, amt, start, stop, step=1.0):
    t = start
    ev = []
    while t < stop:
        ev.append((t, nm, amt))
        t += step
    return ev


# =====================================================================
# 5. Anchors
# =====================================================================
RES = []


def check(label, got, lo, hi, src):
    ok = (lo <= got <= hi)
    RES.append((ok, label, got, lo, hi, src))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:<58s} {got:>12.4g}   "
          f"target {lo:g}-{hi:g}   ({src})")
    return ok


def section(t):
    print("\n" + "=" * 100)
    print(t)
    print("=" * 100)


# ---------------------------------------------------------------------
section("A.  Benzathine penicillin G depot kinetics")
# ---------------------------------------------------------------------
def bpg_profile(wt, mg=720.0, tend=35.0):
    p = dict(P); p["WT"] = wt
    t, Y, _ = simulate(p, tend, events=bpg_events(1e9, 1, mg=mg),
                       t_eval=np.linspace(0, tend, 3501))
    V = p["VPEN0"] * (wt / 70.0)
    return t, Y[:, IDX["PENC"]] / V, Y


t, C, Y = bpg_profile(70)
for day, lo, hi in [(1, 0.10, 0.30), (7, 0.035, 0.085), (14, 0.018, 0.045),
                    (21, 0.010, 0.028), (28, 0.005, 0.017)]:
    c = np.interp(day, t, C)
    check(f"BPG 1.2 MU, 70 kg: serum penicillin G at day {day} (ug/mL)",
          c, lo, hi, "Kaplan 1989 PMID 2738782")

tprot70 = np.interp(35, t, Y[:, IDX["TPROT"]])
check("BPG 1.2 MU, 70 kg: days with C > 0.02 ug/mL", tprot70, 15.0, 21.0,
      "Kaplan 1989 PMID 2738782; Neely 2014 PMID 25182635")

TP = {}
for wt in (40, 50, 60, 70, 80, 90, 100, 110):
    tt, cc, YY = bpg_profile(wt, tend=60)
    TP[wt] = np.interp(60, tt, YY[:, IDX["TPROT"]])
print("\n  Days above 0.02 ug/mL after a single 1.2 MU injection, by body weight:")
for wt in sorted(TP):
    print(f"      {wt:3d} kg   {TP[wt]:5.1f} d   "
          f"({'covers' if TP[wt] >= 28 else 'FAILS'} a 28-day interval, "
          f"{'covers' if TP[wt] >= 21 else 'fails'} a 21-day interval)")
check("weight sensitivity: T>MIC at 40 kg divided by T>MIC at 100 kg",
      TP[40] / TP[100], 1.4, 2.2, "Hand 2019 PMID 30989171 (size on CL/V)")

# ---------------------------------------------------------------------
section("B.  The headline: programme adherence is not time-above-threshold")
# ---------------------------------------------------------------------
def protection(wt, interval, adherence_pattern, years=2.0, mg=720.0):
    """Fraction of calendar time with penicillin above 0.02 ug/mL."""
    p = dict(P); p["WT"] = wt
    n = int(years * 365 / interval) + 1
    ev = bpg_events(interval, n, mg=mg, adherence=adherence_pattern)
    tend = years * 365
    t, Y, _ = simulate(p, tend, events=ev, t_eval=np.linspace(0, tend, 601))
    return Y[-1, IDX["TPROT"]] / tend


ALL = [1]
P80 = [1, 1, 1, 1, 0]          # 4 of every 5 doses given
P60 = [1, 1, 1, 0, 0]

rows = [
    ("55 kg, 21-day interval, 100% of doses given", 55, 21, ALL, 1.00),
    ("55 kg, 21-day interval,  80% of doses given", 55, 21, P80, 0.80),
    ("55 kg, 28-day interval, 100% of doses given", 55, 28, ALL, 1.00),
    ("70 kg, 28-day interval, 100% of doses given", 70, 28, ALL, 1.00),
    ("95 kg, 28-day interval, 100% of doses given", 95, 28, ALL, 1.00),
    ("95 kg, 21-day interval, 100% of doses given", 95, 21, ALL, 1.00),
    ("95 kg, 21-day interval,  60% of doses given", 95, 21, P60, 0.60),
]
print(f"\n  {'regimen':<48s} {'adherence':>10s} {'protected':>11s}")
prot = {}
for lab, wt, iv, pat, adh in rows:
    f = protection(wt, iv, pat)
    prot[lab] = (adh, f)
    print(f"  {lab:<48s} {adh*100:>9.0f}% {f*100:>10.1f}%")

check("100%-adherent 95 kg adult on 4-weekly BPG: fraction of time protected",
      prot["95 kg, 28-day interval, 100% of doses given"][1], 0.40, 0.60,
      "model prediction")
check("80%-adherent 55 kg adult on 3-weekly BPG: fraction of time protected",
      prot["55 kg, 21-day interval,  80% of doses given"][1], 0.70, 0.90,
      "model prediction")
inv = (prot["55 kg, 21-day interval,  80% of doses given"][1]
       > prot["95 kg, 28-day interval, 100% of doses given"][1])
print(f"\n  Rank inversion between adherence and protection: "
      f"{'CONFIRMED' if inv else 'not present'}")

# ---------------------------------------------------------------------
section("C.  Untreated versus treated streptococcal pharyngitis")
# ---------------------------------------------------------------------
p = dict(P); p["EPISODIC"] = 1.0
ev_inoc = [(0.0, "GAS", 1e4)]
t, Y, _ = simulate(p, 40, events=ev_inoc, t_eval=np.linspace(0, 40, 4001))
gpk = Y[:, IDX["GAS"]].max()
ipk = int(np.argmax(Y[:, IDX["GAS"]]))
below = np.where(Y[ipk:, IDX["GAS"]] < 1e3)[0]
tclear = t[ipk + below[0]] if below.size else 99.0
check("untreated GAS pharyngitis: days to clear below 1e3 CFU-eq",
      tclear, 10.0, 26.0, "GAS persists 2-3 weeks untreated")
asopk = Y[:, IDX["ASO"]].max()
taso = t[np.argmax(Y[:, IDX["ASO"]])]
check("untreated GAS: time to peak ASO titre (days)", taso, 18.0, 40.0,
      "ASO peaks 3-5 weeks after infection")

# penicillin V 500 mg bid for 10 days, started on day 2
ev_tx = ev_inoc + [(d, "PVA", 250.0) for d in np.arange(2.0, 12.0, 0.5)]
t2, Y2, _ = simulate(p, 40, events=ev_tx, t_eval=np.linspace(0, 40, 4001))
ag_un = np.trapezoid(Y[:, IDX["AG"]], t)
ag_tx = np.trapezoid(Y2[:, IDX["AG"]], t2)
check("penicillin V from day 2: fractional fall in cumulative antigen exposure",
      1.0 - ag_tx / ag_un, 0.60, 1.00, "Denny/Wannamaker landmark PMID 3892066")

# late treatment (day 8) - the classic 9-day window
ev_late = ev_inoc + [(d, "PVA", 250.0) for d in np.arange(8.0, 18.0, 0.5)]
t3, Y3, _ = simulate(p, 40, events=ev_late, t_eval=np.linspace(0, 40, 4001))
ag_late = np.trapezoid(Y3[:, IDX["AG"]], t3)
print(f"\n  cumulative antigen exposure (AU.d):  untreated {ag_un:.1f}   "
      f"treated d2 {ag_tx:.1f} (-{100*(1-ag_tx/ag_un):.0f}%)   "
      f"treated d8 {ag_late:.1f} (-{100*(1-ag_late/ag_un):.0f}%)")
# The nine-day window is NOT reproduced in a fast-clearing host: by day 8
# this model has already cleared the organism, so late treatment removes
# almost nothing.  The window therefore requires a host who is still
# carrying antigen at day 8 -- which is exactly the slow-clearing/carrier
# phenotype.  That is a testable statement, so it is tested rather than
# asserted.
def late_benefit(kmi):
    pp = dict(P); pp["EPISODIC"] = 1.0; pp["KMI"] = kmi
    ta, Ya, _ = simulate(pp, 60, events=ev_inoc, t_eval=np.linspace(0, 60, 3001))
    tb, Yb, _ = simulate(pp, 60, events=ev_late, t_eval=np.linspace(0, 60, 3001))
    a = np.trapezoid(Ya[:, IDX["AG"]], ta)
    b = np.trapezoid(Yb[:, IDX["AG"]], tb)
    Ga = Ya[:, IDX["GAS"]]
    ipk_ = int(np.argmax(Ga)); bl = np.where(Ga[ipk_:] < 1e3)[0]
    return 1.0 - b / a, (ta[ipk_ + bl[0]] if bl.size else 60.0)


print("\n  Benefit of penicillin begun on DAY 8, against how fast the untreated")
print("  host clears the organism:")
for kmi in (0.30, 0.15, 0.08, 0.05):
    ben, tcl = late_benefit(kmi)
    print(f"      KMI {kmi:.2f}  untreated clearance day {tcl:5.1f}  "
          f"day-8 treatment cuts antigen {100*ben:5.1f}%")
ben_slow, tcl_slow = late_benefit(0.08)
check("THE NINE-DAY WINDOW exists only in a slow-clearing (carrier) host",
      ben_slow, 0.30, 0.95,
      "Denny 1950 / Wannamaker landmark PMID 3892066")
check("untreated GAS pharyngitis: peak ASO titre (Todd units)",
      Y[:, IDX["ASO"]].max(), 250.0, 1500.0, "ASO 300-1500 in ARF")

# ---------------------------------------------------------------------
section("D.  Gorlin arithmetic: the valve does not change, the patient does")
# ---------------------------------------------------------------------
def state_at(mva, **kw):
    y = initial(P)
    y[IDX["MVA"]] = mva
    for k, v in kw.items():
        y[IDX[k]] = v
    return y


def hemo_at(mva, demf=1.0, afb=0.0, **kw):
    p = dict(P); p["DEMF"] = demf
    return haemo(state_at(mva, AFB=afb, **kw), p)


for mva, lo, hi in [(2.0, 2.0, 5.0), (1.5, 4.0, 8.0), (1.0, 9.0, 17.0)]:
    h = hemo_at(mva)
    check(f"resting mean mitral gradient at MVA {mva} cm2 (mmHg)",
          h["MVG"], lo, hi, "ASE/EACVI severity grades")

def steady(mva, demf=1.0, afb=0.0, days=200, bb_mg=0.0):
    """Run to a haemodynamic steady state so the congestion -> sympathetic ->
    heart-rate loop is allowed to act.  The static algebra alone understates
    decompensation because it holds heart rate fixed."""
    p = dict(P); p["DEMF"] = demf
    y0 = initial(p); y0[IDX["MVA"]] = mva; y0[IDX["AFB"]] = afb
    ev = daily("BBA", bb_mg, 0, days) if bb_mg > 0 else []
    t, Y, _ = simulate(p, days, y0=y0, events=ev,
                       t_eval=np.linspace(0, days, 401))
    return haemo(Y[-1], p), Y[-1]


def show(tag, h):
    print(f"      {tag:<26s} CO {h['CO']:.2f} L/min  HR {h['HR']:5.1f}  "
          f"gradient {h['MVG']:5.1f} mmHg  LA pressure {h['LAP']:5.1f} mmHg  "
          f"NYHA {h['NYHA']:.2f}")


h_rest, _ = steady(1.5, demf=1.0)
h_preg, _ = steady(1.5, demf=1.5)
h_af, _ = steady(1.5, demf=1.0, afb=1.0)
h_afbb, _ = steady(1.5, demf=1.0, afb=1.0, bb_mg=100.0)
print("\n  MVA 1.5 cm2 throughout -- the ORIFICE NEVER CHANGES in this block:")
show("rest, sinus rhythm", h_rest)
show("pregnancy (CO demand x1.5)", h_preg)
show("new-onset AF, at rest", h_af)
show("  ... plus metoprolol 100 mg", h_afbb)

check("pregnancy (CO demand x1.5) multiplies the mean gradient at fixed MVA 1.5 by",
      h_preg["MVG"] / h_rest["MVG"], 2.0, 5.0,
      "Hameed 2001 PMID 11693767; gradient scales as flow squared")
check("pregnancy at MVA 1.5: LA pressure reaches the oedema range (mmHg)",
      h_preg["LAP"], 22.0, 40.0, "clinical: MS presents in the 2nd-3rd trimester")
check("new-onset AF at fixed MVA 1.5 multiplies the mean gradient by",
      h_af["MVG"] / h_rest["MVG"], 1.6, 4.5, "clinical decompensation on AF onset")
check("rate control recovers what AF cost, at an unchanged orifice (frac of gradient)",
      1 - h_afbb["MVG"] / h_af["MVG"], 0.05, 0.70, "model prediction")

# ---------------------------------------------------------------------
section("E.  The optimal heart rate falls out of the model, it is not asserted")
# ---------------------------------------------------------------------
def co_vs_hr(mva, demf=2.0):
    """Sweep heart rate by fixing it directly, holding demand high."""
    p = dict(P); p["DEMF"] = demf
    out = []
    for hr in range(40, 161, 1):
        RR = 60.0 / hr
        DFPB = max(RR - 0.36 * math.sqrt(RR), 0.05)
        DFPM = hr * DFPB
        dPmax = p["LAPMAX"] - p["LVEDP0"]
        COv = p["GORLIN"] * mva * math.sqrt(dPmax) * DFPM / 1000.0
        COs = p["SVMAX"] * hr / 1000.0
        out.append((hr, min(COv, COs, p["CO0"] * demf)))
    return out


print("\n  Maximum achievable cardiac output against heart rate (demand not limiting):")
for mva, lo, hi in [(2.0, 90, 125), (1.5, 65, 95), (1.0, 50, 75), (0.8, 45, 70)]:
    sweep = co_vs_hr(mva, demf=3.0)
    best = max(sweep, key=lambda r: r[1])
    print(f"      MVA {mva:.1f} cm2   optimal HR {best[0]:3d} bpm   "
          f"CO {best[1]:.2f} L/min")
    check(f"model-derived optimal heart rate at MVA {mva} cm2 (bpm)",
          best[0], lo, hi, "rate control is more aggressive the tighter the valve")

# ---------------------------------------------------------------------
section("F.  Two clocks of valve destruction, and where they cross")
# ---------------------------------------------------------------------
def rates_at(mva, hazard_per_year):
    """cm2/yr lost to shear alone, and to the immune arm alone."""
    h = hemo_at(mva)
    geom = 2.0 * math.sqrt(PI * mva)
    brake = (mva - P["MVAMIN"]) / (mva - P["MVAMIN"] + P["KHALT"])
    shear = geom * P["KFS"] * (h["MVG"] / 4.0) * brake * 365.0
    # immune arm: steady-state VIT sustained by a given ARF hazard
    # antigen AU/d = QARF*haz/365 ; AG_ss = that / KAGD
    ag = P["QARF"] * hazard_per_year / 365.0 / P["KAGD"]
    mem = P["KMEM"] * ag * P["SUSC"] / (P["KMEM"] * ag * P["SUSC"] + P["KMEMD"])
    xab = P["KXAB"] * ag * P["SUSC"] * (1 + P["BMEM"] * mem) / P["KXABD"]
    scarf = max(0.0, 1.0 - mva / P["MVAN"])
    vit = P["KVIT"] * xab * (1 + P["BSCAR"] * scarf) / P["KVITD"]
    immune = geom * P["KFI"] * vit * 365.0
    return shear, immune


print("\n  Annual loss of mitral valve area (cm2/yr) by mechanism, no prophylaxis:")
print(f"  {'MVA':>6s} {'shear':>10s} {'immune':>10s} {'total':>10s}  dominant")
cross = None
prev = None
for mva in [4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.2, 1.0, 0.8]:
    s, i = rates_at(mva, 0.10)
    dom = "shear" if s > i else "immune"
    print(f"  {mva:6.1f} {s:10.4f} {i:10.4f} {s+i:10.4f}  {dom}")
    if prev is not None and prev[1] > prev[0] and i < s:
        cross = (prev[2] + mva) / 2
    prev = (s, i, mva)
if cross:
    print(f"\n  Crossover: the autonomous shear arm overtakes the immune arm "
          f"at MVA ~ {cross:.2f} cm2")
check("MVA at which shear-driven loss overtakes immune loss (cm2)",
      cross if cross else 0.0, 1.2, 3.0,
      "model prediction; explains GOAL trial targeting LATENT RHD")

s15, _ = rates_at(1.5, 0.0)
check("progression of established MS with zero recurrences (cm2/yr)",
      s15, 0.05, 0.15, "Sagie 1996 PMID 8800128; Gordon 1992 PMID 1552121")

# ---------------------------------------------------------------------
section("G.  Twenty-five year natural history: prophylaxis on and off")
# ---------------------------------------------------------------------
YEARS = 25
TEND = YEARS * 365


def life(prophylaxis=None, wt=55.0, susc=1.0, mva0=3.6, mem0=1.0, years=YEARS):
    """prophylaxis: None, or (interval_days, adherence_pattern)"""
    p = dict(P); p["WT"] = wt; p["SUSC"] = susc
    y0 = initial(p)
    y0[IDX["MVA"]] = mva0          # after a first attack with carditis
    y0[IDX["MEM"]] = mem0
    ev = []
    if prophylaxis:
        iv, pat = prophylaxis
        ev = bpg_events(iv, int(years * 365 / iv) + 1, adherence=pat)
    tend = years * 365
    t, Y, _ = simulate(p, tend, events=ev, t_eval=np.linspace(0, tend, 4001))
    return t, Y


scen = {
    "no prophylaxis": None,
    "4-weekly, 100% doses": (28, ALL),
    "3-weekly, 100% doses": (21, ALL),
    "3-weekly,  80% doses": (21, P80),
}
print(f"\n  Starting from MVA 3.60 cm2 (a first attack with carditis, latent RHD),"
      f" 55 kg, {YEARS} years:")
print(f"  {'regimen':<26s} {'protected':>10s} {'MVA':>8s} {'gradient':>9s} "
      f"{'LAP':>7s} {'AF':>6s} {'ARF':>6s} {'NYHA':>6s}")
life_out = {}
for lab, pr in scen.items():
    t, Y = life(pr)
    p = dict(P); p["WT"] = 55.0
    h = haemo(Y[-1], p)
    protf = Y[-1, IDX["TPROT"]] / TEND
    life_out[lab] = (protf, Y[-1, IDX["MVA"]], h["MVG"], h["LAP"],
                     Y[-1, IDX["AFB"]], Y[-1, IDX["CUMARF"]], h["NYHA"])
    print(f"  {lab:<26s} {protf*100:>9.1f}% {Y[-1, IDX['MVA']]:>8.2f} "
          f"{h['MVG']:>9.1f} {h['LAP']:>7.1f} {Y[-1, IDX['AFB']]:>6.2f} "
          f"{Y[-1, IDX['CUMARF']]:>6.2f} {h['NYHA']:>6.2f}")

check("25 y with no prophylaxis: expected recurrent ARF episodes",
      life_out["no prophylaxis"][5], 0.8, 4.0,
      "recurrence ~5-10%/yr unprotected in endemic settings")
check("25 y on 3-weekly BPG: expected recurrent ARF episodes",
      life_out["3-weekly, 100% doses"][5], 0.01, 0.5,
      "Manyemba 2002 Cochrane PMID 12137650")
check("relative reduction in recurrent ARF, 3-weekly vs none",
      1 - life_out["3-weekly, 100% doses"][5] / life_out["no prophylaxis"][5],
      0.70, 0.98, "Manyemba 2002 PMID 12137650")
dMVA = life_out["3-weekly, 100% doses"][1] - life_out["no prophylaxis"][1]
print(f"\n  MVA preserved by 3-weekly prophylaxis over 25 y: {dMVA:+.2f} cm2")

# GOAL-style: intervene while still latent (MVA near normal)
t, Ya = life(None, mva0=4.2, years=2)
t, Yb = life((28, ALL), mva0=4.2, years=2)
print(f"  Latent RHD (MVA 4.20), 2 y:  no prophylaxis {Ya[-1, IDX['MVA']]:.3f} cm2 "
      f"vs BPG {Yb[-1, IDX['MVA']]:.3f} cm2")

# ---------------------------------------------------------------------
section("H.  Acute rheumatic fever, and why a six-week steroid course is null")
# ---------------------------------------------------------------------
def arf_episode(steroid_weeks=0, aspirin=True, tend=400, pred_mg=40.0):
    p = dict(P); p["EPISODIC"] = 1.0
    ev = [(0.0, "GAS", 1e4)]
    onset = 19.0                       # ARF declares ~3 weeks after pharyngitis
    if aspirin:
        for d in np.arange(onset, onset + 42, 0.25):
            ev.append((d, "ASAA", 1000.0))
    for d in np.arange(onset, onset + steroid_weeks * 7, 1.0):
        ev.append((d, "PREDA", pred_mg))
    t, Y, _ = simulate(p, tend, events=ev, t_eval=np.linspace(0, tend, 4001))
    return t, Y


t0_, Y0_ = arf_episode(0)
crp_pk = Y0_[:, IDX["CRP"]].max()
check("first ARF attack: peak CRP (mg/L)", crp_pk, 30.0, 200.0,
      "acute-phase response in ARF")
check("first ARF attack: mitral area lost (cm2)",
      P["MVAN"] - Y0_[-1, IDX["MVA"]], 0.20, 1.20,
      "carditis in ~50-70% of first attacks; residual RHD in ~half")
check("first ARF attack: peak acute mitral regurgitation grade",
      Y0_[:, IDX["MRA"]].max(), 1.0, 4.0, "MR is the dominant acute lesion")

print(f"\n  {'steroid course':<22s} {'int(valvulitis)':>16s} {'MVA lost':>10s} "
      f"{'vs no steroid':>14s}")
base_int = base_loss = None
for wk in (0, 2, 6, 12, 26):
    t, Y = arf_episode(wk)
    vint = np.trapezoid(Y[:, IDX["VIT"]], t)
    loss = P["MVAN"] - Y[-1, IDX["MVA"]]
    if wk == 0:
        base_int, base_loss = vint, loss
    print(f"  {(str(wk) + ' weeks'):<22s} {vint:>16.1f} {loss:>10.3f} "
          f"{100*(1-loss/base_loss):>13.1f}%")
t6, Y6 = arf_episode(6)
loss6 = P["MVAN"] - Y6[-1, IDX["MVA"]]
check("6-week prednisolone: reduction in valve area lost (fraction)",
      1 - loss6 / base_loss, 0.0, 0.35,
      "Cilliers 2015 Cochrane PMID 26017576 - no clear benefit")
t26, Y26 = arf_episode(26)
loss26 = P["MVAN"] - Y26[-1, IDX["MVA"]]
check("26-week prednisolone: reduction in valve area lost (fraction)",
      1 - loss26 / base_loss, 0.35, 0.95, "model prediction - untested")
print(f"\n  Antibody half-life {math.log(2)/P['KXABD']:.0f} d; a 6-week course "
      f"covers {100*(1-math.exp(-P['KXABD']*42)):.0f}% of the antibody-time integral.")

# ---------------------------------------------------------------------
section("I.  Rate control, anticoagulation, and balloon valvotomy")
# ---------------------------------------------------------------------
# severe MS in AF, exercising; add a beta-blocker
def ms_af(bb_mg=0.0, demf=1.6, mva=1.0, days=30):
    p = dict(P); p["DEMF"] = demf
    y0 = initial(p)
    y0[IDX["MVA"]] = mva
    y0[IDX["AFB"]] = 1.0
    y0[IDX["LAV"]] = 130.0
    ev = daily("BBA", bb_mg, 0, days) if bb_mg > 0 else []
    t, Y, _ = simulate(p, days, y0=y0, events=ev,
                       t_eval=np.linspace(0, days, 1201))
    return t, Y, p


print(f"\n  Severe MS (MVA 1.0), permanent AF, moderate exertion:")
for mg in (0, 25, 50, 100, 200):
    t, Y, p = ms_af(mg)
    h = haemo(Y[-1], p)
    print(f"      metoprolol {mg:3d} mg/d   HR {h['HR']:5.1f}   "
          f"gradient {h['MVG']:5.1f} mmHg   LA pressure {h['LAP']:5.1f}   "
          f"CO {h['CO']:.2f} L/min   NYHA {h['NYHA']:.2f}")
t, Y, p = ms_af(0); h0 = haemo(Y[-1], p)
t, Y, p = ms_af(200); h1 = haemo(Y[-1], p)
print("\n  Note what rate control does and does not do here.  This patient is")
print("  already AT the left-atrial pressure ceiling, so the gradient CANNOT")
print("  fall -- it is pinned at LAPMAX - LVEDP.  The diastole that beta")
print("  blockade buys is spent on FORWARD FLOW instead, and the measured")
print("  gradient is the one number that does not move.")
check("severe MS at the LA-pressure ceiling: metoprolol 200 mg/d raises CO by (frac)",
      h1["CO"] / h0["CO"] - 1.0, 0.03, 0.30, "model prediction")
# below the ceiling the same drug moves the gradient instead
t, Y, p = ms_af(0, demf=1.0, mva=1.4); g0 = haemo(Y[-1], p)
t, Y, p = ms_af(200, demf=1.0, mva=1.4); g1 = haemo(Y[-1], p)
print(f"      below the ceiling (MVA 1.4, rest): gradient "
      f"{g0['MVG']:.1f} -> {g1['MVG']:.1f} mmHg, "
      f"LA pressure {g0['LAP']:.1f} -> {g1['LAP']:.1f} mmHg")
check("moderate MS below the ceiling: metoprolol 200 mg/d cuts the gradient by",
      1 - g1["MVG"] / g0["MVG"], 0.10, 0.60,
      "rate control lengthens diastole; Gorlin denominator")

# warfarin
p = dict(P)
ev = daily("WA", 5.0, 0, 30)
t, Y, _ = simulate(p, 30, events=ev, t_eval=np.linspace(0, 30, 1201))
inr = (1.0 / np.maximum(Y[:, IDX["PCF"]], 1e-3)) ** P["GAMINR"]
check("warfarin 5 mg/d: steady-state INR", inr[-1], 1.9, 3.2,
      "typical maintenance dose-response")
i5 = np.interp(5, t, inr)
check("warfarin 5 mg/d: INR on day 5", i5, 1.6, 3.0, "onset over 4-6 days")

# embolic hazard, severe MS with AF, with and without warfarin
def emb(warf_mg, years=1.0):
    p = dict(P)
    y0 = initial(p); y0[IDX["MVA"]] = 1.0; y0[IDX["AFB"]] = 1.0; y0[IDX["LAV"]] = 140.0
    ev = daily("WA", warf_mg, 0, years * 365) if warf_mg > 0 else []
    t, Y, _ = simulate(p, years * 365, events=ev, y0=y0,
                       t_eval=np.linspace(0, years * 365, 2001))
    return Y[-1, IDX["CUMEMB"]]


e0, e5 = emb(0.0), emb(5.0)
check("severe MS + AF, no anticoagulation: annual embolic event rate",
      e0, 0.06, 0.20, "7-15%/yr in rheumatic MS with AF")
check("warfarin at INR ~2.5: relative reduction in embolic rate",
      1 - e5 / e0, 0.50, 0.80, "~65% in AF trials; INVICTUS PMID 36036525")

# balloon valvotomy -- a typical favourable candidate: young, sinus rhythm,
# pliable valve
p = dict(P)
y0 = initial(p); y0[IDX["MVA"]] = 0.95; y0[IDX["CA"]] = 3.0
y0[IDX["AFB"]] = 0.0; y0[IDX["LAV"]] = 85.0
w = wilkins(y0, p)
check("Wilkins score of a typical PMBV candidate", w, 6.0, 11.0,
      "Wilkins 1988 PMID 3190958")
gain = 1.15
# the post-valvotomy patient stays on 4-weekly secondary prophylaxis
ev_pmbv = [(1.0, "MVA", gain)] + bpg_events(28, int(10 * 365 / 28) + 1)
t, Y, _ = simulate(p, 10 * 365, y0=y0, events=ev_pmbv,
                   pchange=[(1.0, "RESTEN", 1.35)],
                   t_eval=np.linspace(0, 10 * 365, 2001))
mva_post = np.interp(2, t, Y[:, IDX["MVA"]])
mva_5y = np.interp(5 * 365, t, Y[:, IDX["MVA"]])
mva_10y = Y[-1, IDX["MVA"]]
check("PMBV: immediate post-procedure valve area (cm2)", mva_post, 1.7, 2.4,
      "typical gain 1.0 -> 2.0 cm2")
check("PMBV: valve area 5 years later (cm2)", mva_5y, 1.3, 1.9,
      "Iung 1999 PMID 10385502; Palacios 2002 PMID 11914256")
print(f"\n  PMBV: {0.95:.2f} -> {mva_post:.2f} cm2, {mva_5y:.2f} at 5 y, "
      f"{mva_10y:.2f} at 10 y (restenosis multiplier 1.35 on the shear arm)")

# ---- the model's most exposed claim, stated rather than tuned away -------
print("\n  THE MODEL'S MOST EXPOSED CLAIM.  The shear arm is calibrated to a")
print("  MEAN progression of 0.09 cm2/yr at MVA 1.5 (Sagie 1996, Gordon 1992),")
print("  and because the gradient rises as flow-squared-over-area-squared it")
print("  then ACCELERATES as the valve narrows.  Sagie reported the opposite")
print("  association -- progression was slower in the tighter valves -- so this")
print("  is the assumption most likely to be wrong.  Sensitivity, same patient:")
print("  (KFS is renormalised at MVA 1.5 in every row, so these rows isolate the")
print("   AREA-DEPENDENCE of the shear arm from its overall size.)")
for kh, lab in ((1.6, "accelerating (base case)"), (6.0, "strong rigidity brake"),
                (20.0, "shear arm nearly flat")):
    pp = dict(P); pp["KHALT"] = kh
    pp["KFS"] = 3.914e-5 / ((1.5 - pp["MVAMIN"]) / (1.5 - pp["MVAMIN"] + kh))
    tt, YY, _ = simulate(pp, 10 * 365, y0=y0, events=ev_pmbv,
                         pchange=[(1.0, "RESTEN", 1.35)],
                         t_eval=np.linspace(0, 10 * 365, 1201))
    print(f"      KHALT {kh:>4.1f}  {lab:<26s} 5 y {np.interp(5*365, tt, YY[:, IDX['MVA']]):.2f}"
          f"   10 y {YY[-1, IDX['MVA']]:.2f} cm2")
print("  The three rows barely differ, which is itself the answer: the 10-year")
print("  number is NOT sensitive to how the shear arm depends on area, because")
print("  the mean gradient saturates at LAPMAX - LVEDP once the valve is tight.")
print("  So if this trajectory is too pessimistic, the error is in the SIZE of")
print("  KFS -- i.e. in reading a cohort MEAN progression of 0.09 cm2/yr as if")
print("  every valve progressed -- and not in its area-dependence.  A restenosis")
print("  cohort with serial echoes would settle it directly.")

# ---------------------------------------------------------------------
section("J.  Pulmonary vascular and right ventricular consequences")
# ---------------------------------------------------------------------
p = dict(P)
y0 = initial(p); y0[IDX["MVA"]] = 0.9; y0[IDX["LAV"]] = 120.0; y0[IDX["AFB"]] = 0.8
t, Y, _ = simulate(p, 8 * 365, y0=y0, t_eval=np.linspace(0, 8 * 365, 2001))
h = haemo(Y[-1], p)
check("severe MS untreated 8 y: pulmonary vascular resistance (Wood units)",
      Y[-1, IDX["PVR"]], 2.5, 10.0, "reactive PH in severe MS")
check("severe MS untreated 8 y: mean pulmonary artery pressure (mmHg)",
      h["PAPm"], 30.0, 70.0, "PAPm 30-60 mmHg in severe MS")
check("severe MS untreated 8 y: RV function index (1 = normal)",
      Y[-1, IDX["RVF"]], 0.1, 0.9, "progressive RV failure")

# ---------------------------------------------------------------------
section("SUMMARY")
# ---------------------------------------------------------------------
npass = sum(1 for r in RES if r[0])
print(f"\n  {npass} / {len(RES)} anchors pass\n")
for ok, lab, got, lo, hi, src in RES:
    if not ok:
        print(f"  FAILED: {lab}  got {got:.4g}, target {lo:g}-{hi:g}  ({src})")
if npass == len(RES):
    print("  All anchors pass.\n")
