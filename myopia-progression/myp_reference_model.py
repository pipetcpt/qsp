#!/usr/bin/env python3
"""
myp_slow.py -- independent Python re-implementation of the slow (disease) half
of the myopia-progression QSP model.  Dependency-free RK4.

Written BEFORE the mrgsolve file so that every equation is actually integrated
and every number quoted in the README is a measured output.  The fast atropine
PK subsystem lives in pk_atropine.py and is solved separately (it is linear and
autonomous in the disease states); its periodic-steady-state occupancy
summaries enter here through ATRO_X1_*.

TIME UNIT: days.  AGE in years.

DEFECT LOG (things the first integration exposed, all fixed below)
  D1  The effector cascade was written with LINEAR couplings
      (TGFB = 1 - ATG*(RA-1), TIMP2 = 1 + ATI*(TGFB-1)).  With a drive of
      O(1) that drove TGF-beta and then TIMP-2 NEGATIVE, MMP2/TIMP2 hit 49,
      scleral creep hit 125 and the eye reached -246 D in three years.  All
      couplings are now LOG-LINEAR (target = exp(gain*drive), i.e. a power
      law in the upstream variable), which cannot change sign and which is
      also the right functional form for a transcription-mediated cascade.
  D2  CHT_REF (the choroidal thickness that measured axial length is
      referenced to) was set to the analytic baseline, but the 120-day
      settling phase moves the choroid, so t = 0 did not reproduce the
      prescribed baseline refraction (-1.78 D instead of -1.50 D).  CHT_REF
      is now read off AFTER settling.
"""

import math

YR = 365.0

# =========================================================================
# ATROPINE occupancy vs dose (summaries of the linear PK in pk_atropine.py)
# =========================================================================
# The posterior (choroid/sclera) concentration stays far below the
# low-affinity growth-control EC50, so posterior occupancy is nearly LINEAR
# in dose; the two anterior sites are not.  Each site is summarised by
# x1 = C_avg/Kd at 0.01%, with x scaling linearly in dose x adherence.
ATRO_X1_POST = 0.0410     # fitted to the 0.01 / 0.05 / 1% axial-efficacy triplet
ATRO_X1_IRIS = 0.1880     # fitted to LAMP 1-yr photopic pupil diameter
ATRO_X1_CIL = 0.0617      # fitted to LAMP 1-yr accommodative amplitude


def occ(x1, pct, adh=1.0):
    x = x1 * (pct / 0.01) * adh
    return x / (1.0 + x)


P = dict(
    NVIT=1.336, NAQ=1.336, VERTEX=0.012,
    CRAD0=7.80, ACD0=3.60, LT0=3.45, PLENS0=22.60,
    KCRAD=0.002 / YR, KACD=0.012 / YR, KLTAGE=-0.020 / YR, TAUAS=30.0,
    KACDA=0.16, KLTA=-0.14,
    KPL=0.50 / YR, TAUPL=7.0, KMC=0.05, KPLA=0.0 / YR,
    LUXREF=14500.0, TAUDA=2.0, KDA=0.90, PDA=0.35,
    G0=0.254, KPAR=0.35, KGRS=1.20, KETH=0.60, PNEAR=0.50, NEARREF=3.0,
    KH=0.55, DH=1.20, KM=1.40, DM=0.80, WPER=0.70, WCEN=0.30, KSIG=1.00,
    RPR0=-0.30, KRPR=0.25, PRPR=0.80, TAURPR=60.0,
    ALAG0=0.65, KLAGM=0.10, KLAGA=1.60,
    # --- log-linear effector cascade (D1) ---
    TAURA=3.0, ARA=0.55,
    TAUNO=2.0, ANO=0.45,
    TAUTG=14.0, ATG=0.60,
    TAUMM=10.0, AMM=0.70, AMMH=0.40,
    TAUTI=14.0, ATI=0.90,
    KGS=0.060, KGD=0.060, AGD=0.80,
    KCS=0.030, KCD=0.030, ACOD=0.70,
    TAULX=45.0, ALX=0.70, TAUSF=20.0, ASF=0.80, TAUHY=15.0, AHY=0.90,
    TAUCR=20.0, PC1=0.50, PC2=0.35, PC3=0.50, PC4=0.25, PC5=0.30,
    KEXP=0.0771 / YR, TAUAGE=6.0, PEXP=1.50, PIOP=0.50, IOPREF=15.0,
    KTHIN=0.55, KDEGS=0.06, SCT0=1000.0,
    KSTAPH=2.0e-4, ALSTAPH=26.5, KLACQ=8.0e-4,
    CHT0=300.0, TAUCH=5.0, KCH=0.10, KCHA=0.055, KCHR=0.075,
    KCHAL=25.0, TAUCBF=7.0,
    TAURUP=90.0, KRUP=0.55,
    EMAXATR=1.00, ED50ATR=0.070, HATR=1.22, KDIR=0.90,
    PUPD0=4.60, EMAXPD=3.80, AAMP0=13.40,
    TAUOK=5.0, KOKPD=1.80, KOKCR=0.42,
    TAUPBA=1.5, TAUPBC=30.0, EMAXRL=0.62,
    KAMX=6.0, KELMX=2.4, VMX=40.0, EMAXMX=0.38, EC50MX=1.2,
    ADH0=0.92, TAUADH=60.0, KADHP=0.25, KADHN=0.35,
    IOP0=15.5,
    AVI=-2.060, BVI=0.950,
    KMMD=1.1e-5, BMMD=0.95, KRD=2.2e-6, BRD=0.42,
    KCNV=3.0e-6, BCNV=1.05, KGLC=1.4e-5, BGLC=0.30, KCAT=2.0e-5, BCAT=0.22,
)

NAMES = [
    "DA", "RA", "NOSIG", "TGFB", "MMP2", "TIMP2", "SHYP",
    "GAG", "COL1", "LOX", "SFIB", "CREEPS", "SCT", "PSTAPH", "LACQ",
    "CHT", "CBF", "CHOX",
    "ALS", "ACD", "LT", "PLENS", "CRAD", "RPRS",
    "MRUP", "ADH", "AAMPS", "NEARWS",
    "OKEPI", "PBMA", "PBMC", "MXGUT", "MXPLA",
    "IOPS",
    "AXCUM", "SERAUC", "HMMD", "HRD", "HCNV", "HGLC", "HCAT",
]
IX = {n: i for i, n in enumerate(NAMES)}
N = len(NAMES)

POSITIVE = ("DA", "RA", "NOSIG", "TGFB", "MMP2", "TIMP2", "SHYP", "GAG",
            "COL1", "LOX", "SFIB", "CREEPS", "CHT", "CBF", "CHOX", "MRUP",
            "ADH", "AAMPS")


def pw(b, e):
    """Power with a clamped base.

    DEFECT D3: RK4 evaluates the RHS at intermediate stages that are NOT
    guaranteed to respect the positivity of the states, and a negative base
    raised to a fractional exponent silently returns a COMPLEX number in
    Python, which then propagates until an unrelated comparison raises.  In
    the ATOM2 washout runs a transient negative CBF did exactly this.  Every
    power-law coupling in the cascade now goes through this function.
    """
    return (b if b > 1e-6 else 1e-6) ** e


def optics(AL, ACD, LT, PLENS, CRAD, p):
    """Exact paraxial two-thin-lens eye -> spectacle-plane SER (D)."""
    PCORN = (p["NAQ"] - 1.0) / (CRAD / 1000.0)
    d1 = ACD + LT / 2.0
    d2 = max(AL - d1, 5.0)
    L2req = 1000.0 * p["NVIT"] / d2
    L2 = L2req - PLENS
    L1 = L2 / (1.0 + (d1 / 1000.0 / p["NAQ"]) * L2)
    FC = L1 - PCORN
    return FC / (1.0 + p["VERTEX"] * FC), FC, PCORN


def resp(D, p):
    if D > 0.0:
        return p["KH"] * D / (1.0 + D / p["DH"])
    a = -D
    return -p["KM"] * a / (1.0 + a / p["DM"])


def make_rhs(sc, p):
    OUTD = sc.get("OUTD", 1.0); NEARD = sc.get("NEARD", 3.0)
    NPAR = sc.get("NPAR", 2); GRS = sc.get("GRS", 0.70); ETHN = sc.get("ETHN", 1.0)
    TRTPD = sc.get("TRTPD", 0.15); OK_ON = sc.get("OK_ON", 0.0)
    RLRL = sc.get("RLRL", 0.0); ATRO = sc.get("ATRO", 0.0)
    ATRO_STOP = sc.get("ATRO_STOP", 1e9); ATRO_TAPER = sc.get("ATRO_TAPER", 0.0)
    OPT_STOP = sc.get("OPT_STOP", 1e9)
    TRT_START = sc.get("TRT_START", 0.0); FUIRIS_R = sc.get("FUIRIS_R", 1.0)
    AGE0 = sc.get("AGE0", 8.0)
    ADHFIX = sc.get("ADHFIX", None)
    LUXH = OUTD * 10000.0 + (16.0 - OUTD) * 300.0
    GBASE = (p["G0"] * (1.0 + p["KPAR"] * NPAR)
             * (1.0 + p["KGRS"] * (GRS - 0.50))
             * (1.0 + p["KETH"] * ETHN)
             * (NEARD / p["NEARREF"]) ** p["PNEAR"])

    def atro_pct(t):
        if t < TRT_START:
            return 0.0
        if t < ATRO_STOP:
            return ATRO
        if ATRO_TAPER > 0.0:
            return ATRO * max(0.0, 1.0 - (t - ATRO_STOP) / ATRO_TAPER)
        return 0.0

    def opt_on(t):
        return 1.0 if (TRT_START <= t < OPT_STOP) else 0.0

    def rhs(t, y):
        d = [0.0] * N
        g = lambda n: y[IX[n]]
        AGE = AGE0 + t / YR

        pct = atro_pct(t)
        ADH = ADHFIX if ADHFIX is not None else g("ADH")
        OCCP = occ(ATRO_X1_POST, pct, ADH) if pct > 0 else 0.0
        OCCI = occ(ATRO_X1_IRIS * FUIRIS_R, pct, ADH) if pct > 0 else 0.0
        OCCC = occ(ATRO_X1_CIL * FUIRIS_R, pct, ADH) if pct > 0 else 0.0
        # The posterior occupancy is near-linear in dose (the PK run shows
        # C_choroid stays far below the low-affinity EC50), so the drive
        # suppression is written as a hyperbola in DOSE, with ED50ATR
        # directly interpretable as a w/v percentage.
        # DEFECT D4: written first as a plain hyperbola in dose (Hill = 1),
        # which could NOT fit the 0.01 / 0.05 / 1% axial arms simultaneously --
        # the fit drove EMAXATR and KDIR to their bounds and still left the
        # 0.05% arm 9 percentage points short.  The published triplet
        # (12.2 / 51.2 / 97.6% axial reduction at a 1 / 5 / 100-fold dose
        # ratio) implies a Hill slope of ~1.22 in BOTH intervals, so the
        # exponent is now explicit and EMAXATR is FIXED at 1 (structural:
        # complete blockade of the modifiable drive).
        dose_eff = pct * ADH
        if dose_eff > 0:
            dh = pw(dose_eff, p["HATR"])
            EATR = p["EMAXATR"] * dh / (dh + pw(p["ED50ATR"], p["HATR"]))
        else:
            EATR = 0.0
        EMX = p["EMAXMX"] * g("MXPLA") / (p["EC50MX"] + g("MXPLA"))
        ERL = p["EMAXRL"] * (0.35 * g("PBMA") + 0.65 * g("PBMC"))

        CHT = g("CHT"); ALS = g("ALS")
        ALMEAS = ALS - (CHT - sc["_CHT_REF"]) / 1000.0
        CRAD_EFF = g("CRAD") + p["KOKCR"] * g("OKEPI")
        SER, FC, PCORN = optics(ALMEAS, g("ACD"), g("LT"), g("PLENS"), CRAD_EFF, p)
        SER_TRUE, _, _ = optics(ALMEAS, g("ACD"), g("LT"), g("PLENS"), g("CRAD"), p)
        MYO = max(0.0, -SER_TRUE)

        RPRT = p["RPR0"] + p["KRPR"] * pw(MYO, p["PRPR"])
        d[IX["RPRS"]] = (RPRT - g("RPRS")) / p["TAURPR"]
        DPER = g("RPRS") + opt_on(t) * TRTPD - p["KOKPD"] * g("OKEPI")
        AAMP = g("AAMPS")
        ALAG = p["ALAG0"] + p["KLAGM"] * MYO + p["KLAGA"] * OCCC
        SIGRAW = (p["WPER"] * resp(DPER, p)
                  + p["WCEN"] * (NEARD / 12.0) * resp(ALAG, p))

        FDA = 1.0 / (1.0 + p["KDA"] * (g("DA") - 1.0))
        SIG = GBASE + p["KSIG"] * SIGRAW
        SIGEFF = max(SIG * FDA * (1 - EATR) * (1 - EMX) * (1 - ERL) * g("MRUP"),
                     -0.90)

        d[IX["DA"]] = ((LUXH / p["LUXREF"]) ** p["PDA"] - g("DA")) / p["TAUDA"]

        # ---- log-linear cascade: every target is a positive power law ----
        RA = g("RA"); TGFB = g("TGFB"); MMP2 = g("MMP2"); TIMP2 = g("TIMP2")
        SHYP = g("SHYP")
        d[IX["RA"]] = (math.exp(p["ARA"] * SIGEFF) - RA) / p["TAURA"]
        d[IX["NOSIG"]] = (math.exp(-p["ANO"] * SIGEFF) - g("NOSIG")) / p["TAUNO"]
        d[IX["TGFB"]] = (pw(RA, -p["ATG"]) - TGFB) / p["TAUTG"]
        d[IX["MMP2"]] = (pw(RA, p["AMM"]) * pw(SHYP, p["AMMH"]) - MMP2) / p["TAUMM"]
        d[IX["TIMP2"]] = (pw(TGFB, p["ATI"]) - TIMP2) / p["TAUTI"]

        MTR = MMP2 / (TIMP2 if TIMP2 > 1e-4 else 1e-4)
        GAG = g("GAG"); COL1 = g("COL1")
        d[IX["GAG"]] = p["KGS"] * TGFB - p["KGD"] * GAG * pw(MTR, p["AGD"])
        d[IX["COL1"]] = p["KCS"] * TGFB - p["KCD"] * COL1 * pw(MTR, p["ACOD"])
        d[IX["LOX"]] = (pw(TGFB, p["ALX"]) - g("LOX")) / p["TAULX"]
        d[IX["SFIB"]] = (pw(TGFB, p["ASF"]) - g("SFIB")) / p["TAUSF"]
        d[IX["SHYP"]] = (pw(g("CBF"), -p["AHY"]) - SHYP) / p["TAUHY"]

        CRTGT = (pw(MTR, p["PC1"]) * pw(GAG, -p["PC2"]) * pw(COL1, -p["PC3"])
                 * pw(g("LOX"), -p["PC4"]) * pw(SHYP, p["PC5"]))
        d[IX["CREEPS"]] = (CRTGT - g("CREEPS")) / p["TAUCR"]

        PHI = math.exp(-(AGE - 7.0) / p["TAUAGE"])
        # DIRECT scleral limb: ATOM1 measured -0.02 mm of axial change over
        # TWO YEARS on 1% atropine.  The model floor with the myopigenic drive
        # fully removed is the EMMETROPIC growth rate (~0.047 mm/yr at this
        # age), so no amount of drive suppression can reach the ATOM1 value.
        # High-dose atropine must therefore also act DIRECTLY on the sclera
        # (muscarinic receptors are present on scleral fibroblasts), which is
        # the KDIR term.  KDIR is fitted to the 1% arm alone.
        dALS = (p["KEXP"] * PHI * pw(g("CREEPS"), p["PEXP"])
                * pw(g("IOPS") / p["IOPREF"], p["PIOP"])
                * max(1.0 - p["KDIR"] * EATR, 0.02))
        d[IX["ALS"]] = dALS
        d[IX["AXCUM"]] = dALS

        SCT = g("SCT")
        d[IX["SCT"]] = (-p["KTHIN"] * dALS * 1000.0 * (SCT / p["SCT0"])
                        - p["KDEGS"] * max(0.0, MTR - 1.0))
        d[IX["PSTAPH"]] = (p["KSTAPH"] * max(0.0, ALS - p["ALSTAPH"])
                           * (p["SCT0"] / max(SCT, 100.0)))
        d[IX["LACQ"]] = p["KLACQ"] * g("PSTAPH") * max(0.0, ALS - p["ALSTAPH"])

        CHT0E = p["CHT0"] - p["KCHAL"] * max(0.0, ALS - 23.50)
        CHTTGT = CHT0E * (1.0 - p["KCH"] * SIGEFF + p["KCHA"] * EATR
                          + p["KCHR"] * (0.5 * g("PBMA") + 0.5 * g("PBMC")))
        d[IX["CHT"]] = (max(CHTTGT, 40.0) - CHT) / p["TAUCH"]
        d[IX["CBF"]] = (CHT / max(CHT0E, 50.0) - g("CBF")) / p["TAUCBF"]
        d[IX["CHOX"]] = (1.0 / max(g("CBF"), 0.2) - g("CHOX")) / p["TAUCH"]

        FMYO = 1.0 / (1.0 + p["KMC"] * MYO)
        d[IX["PLENS"]] = (-p["KPL"] * math.exp(-(AGE - 7.0) / p["TAUPL"]) * FMYO
                          - p["KPLA"] * EATR)
        d[IX["ACD"]] = ((p["ACD0"] + p["KACD"] * YR * (AGE - 7.0)
                         + p["KACDA"] * OCCC) - g("ACD")) / p["TAUAS"]
        d[IX["LT"]] = ((p["LT0"] + p["KLTAGE"] * YR * (AGE - 7.0)
                        + p["KLTA"] * OCCC) - g("LT")) / p["TAUAS"]
        d[IX["CRAD"]] = p["KCRAD"]

        d[IX["AAMPS"]] = (p["AAMP0"] * (1.0 - OCCC) - AAMP) / 3.0
        PUPD = p["PUPD0"] + p["EMAXPD"] * OCCI
        NEARBLUR = 1.0 - AAMP / p["AAMP0"]
        ADHT = max(p["ADH0"] * (1.0 - p["KADHP"] * OCCI
                                - p["KADHN"] * NEARBLUR), 0.15)
        d[IX["ADH"]] = (ADHT - g("ADH")) / p["TAUADH"]
        d[IX["NEARWS"]] = (NEARD - g("NEARWS")) / 30.0
        d[IX["MRUP"]] = (1.0 + p["KRUP"] * EATR - g("MRUP")) / p["TAURUP"]

        d[IX["OKEPI"]] = (OK_ON * opt_on(t) - g("OKEPI")) / p["TAUOK"]
        d[IX["PBMA"]] = (RLRL * opt_on(t) - g("PBMA")) / p["TAUPBA"]
        d[IX["PBMC"]] = (RLRL * opt_on(t) - g("PBMC")) / p["TAUPBC"]
        d[IX["MXGUT"]] = -p["KAMX"] * g("MXGUT")
        d[IX["MXPLA"]] = p["KAMX"] * g("MXGUT") / p["VMX"] - p["KELMX"] * g("MXPLA")
        d[IX["IOPS"]] = (p["IOP0"] - g("IOPS")) / 30.0

        d[IX["SERAUC"]] = MYO / YR
        d[IX["HMMD"]] = p["KMMD"] * math.exp(p["BMMD"] * (ALMEAS - 26.5))
        d[IX["HRD"]] = p["KRD"] * math.exp(p["BRD"] * MYO)
        d[IX["HCNV"]] = p["KCNV"] * math.exp(p["BCNV"] * max(0.0, ALMEAS - 26.0))
        d[IX["HGLC"]] = p["KGLC"] * math.exp(p["BGLC"] * MYO)
        d[IX["HCAT"]] = p["KCAT"] * math.exp(p["BCAT"] * MYO)

        aux = dict(AGE=AGE, SER=SER, SER_TRUE=SER_TRUE, AL=ALMEAS, ALS=ALS,
                   PUPD=PUPD, AAMP=AAMP, OCCP=OCCP, OCCI=OCCI, OCCC=OCCC,
                   EATR=EATR, ERL=ERL, EMX=EMX, SIGEFF=SIGEFF, CHT=CHT,
                   CREEPS=g("CREEPS"), MTR=MTR, ADH=ADH, ALAG=ALAG, DPER=DPER,
                   PHI=PHI, MRUP=g("MRUP"), SCT=SCT, PLENS=g("PLENS"),
                   ALCR=ALMEAS / g("CRAD"), pct=pct, dALS=dALS * YR)
        return d, aux
    return rhs


def solve_AL0(SER0, p, ACD, LT, PLENS, CRAD):
    lo, hi = 18.0, 34.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        s, _, _ = optics(mid, ACD, LT, PLENS, CRAD, p)
        if s > SER0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def run(sc, years=10.0, dt=0.05, p=P, record_every=91.25):
    ACD = p["ACD0"]; LT = p["LT0"]; PLENS = p["PLENS0"]; CRAD = p["CRAD0"]
    AL0 = sc.get("AL0") or solve_AL0(sc["SER0"], p, ACD, LT, PLENS, CRAD)
    sc["_CHT_REF"] = p["CHT0"] - p["KCHAL"] * max(0.0, AL0 - 23.50)
    y = [0.0] * N
    for n in POSITIVE:
        y[IX[n]] = 1.0
    y[IX["RPRS"]] = p["RPR0"] + p["KRPR"] * max(0.0, -sc["SER0"]) ** p["PRPR"]
    y[IX["SCT"]] = p["SCT0"]; y[IX["CHT"]] = sc["_CHT_REF"]
    y[IX["ALS"]] = AL0; y[IX["ACD"]] = ACD; y[IX["LT"]] = LT
    y[IX["PLENS"]] = PLENS; y[IX["CRAD"]] = CRAD
    y[IX["ADH"]] = p["ADH0"]; y[IX["AAMPS"]] = p["AAMP0"]
    y[IX["NEARWS"]] = sc.get("NEARD", 3.0); y[IX["IOPS"]] = p["IOP0"]

    # --- settle fast states with the eye and its slow optics frozen -------
    # Cached: the settling phase is run with NO treatment, so it depends only
    # on the patient descriptors and the disease parameters, never on the
    # treatment arm.  Without this the calibration loops re-derive the same
    # baseline hundreds of times.
    ckey = (round(p["G0"], 10), round(p["KEXP"], 14), round(p["KCH"], 8),
            sc["SER0"], sc.get("AGE0", 8.0), sc.get("NPAR", 2),
            sc.get("GRS", 0.70), sc.get("ETHN", 1.0), sc.get("OUTD", 1.0),
            sc.get("NEARD", 3.0), round(AL0, 10))
    if ckey in _SETTLE_CACHE:
        y = list(_SETTLE_CACHE[ckey][0])
        sc["_CHT_REF"] = _SETTLE_CACHE[ckey][1]
        return _finish(sc, p, y, AL0, years, dt, record_every)
    sc_settle = dict(sc); sc_settle["TRT_START"] = 1e9; sc_settle["ATRO"] = 0.0
    sc_settle["_CHT_REF"] = sc["_CHT_REF"]
    fs = make_rhs(sc_settle, p)
    frozen = ("ALS", "AXCUM", "SERAUC", "PLENS", "SCT", "ACD", "LT", "CRAD",
              "HMMD", "HRD", "HCNV", "HGLC", "HCAT", "PSTAPH", "LACQ")
    for _ in range(int(200.0 / 0.05)):
        dd, _ = fs(0.0, y)
        for i in range(N):
            if NAMES[i] in frozen:
                continue
            y[i] += 0.05 * dd[i]
            if NAMES[i] in POSITIVE and y[i] < 1e-3:
                y[i] = 1e-3
    # D2: reference the measured axial length to the SETTLED choroid, so that
    # t = 0 reproduces the prescribed baseline refraction exactly.
    sc["_CHT_REF"] = y[IX["CHT"]]
    _SETTLE_CACHE[ckey] = (list(y), sc["_CHT_REF"])
    return _finish(sc, p, y, AL0, years, dt, record_every)


_SETTLE_CACHE = {}


def _finish(sc, p, y, AL0, years, dt, record_every):
    f = make_rhs(sc, p)
    n = int(round(years * YR / dt))
    out, nxt = [], 0.0
    mx = sc.get("MXDOSE", 0.0); last_day = -1
    for k in range(n + 1):
        t = k * dt
        if mx > 0:
            day = int(t)
            if day != last_day and t >= sc.get("TRT_START", 0.0):
                y[IX["MXGUT"]] += mx / 194.19 * 1000.0
                last_day = day
        if t >= nxt - 1e-9:
            _, aux = f(t, y)
            rec = dict(t=t, yr=t / YR); rec.update(aux)
            for nm in ("AXCUM", "SERAUC", "HMMD", "HRD", "HCNV", "HGLC",
                       "HCAT", "PSTAPH", "LACQ", "GAG", "COL1", "LOX",
                       "TIMP2", "MMP2", "TGFB", "RA", "CBF", "DA", "OKEPI",
                       "NOSIG", "SFIB", "CHOX"):
                rec[nm] = y[IX[nm]]
            out.append(rec); nxt += record_every
        if k == n:
            break
        k1, _ = f(t, y)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(N)]
        k2, _ = f(t + 0.5 * dt, y2)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(N)]
        k3, _ = f(t + 0.5 * dt, y3)
        y4 = [y[i] + dt * k3[i] for i in range(N)]
        k4, _ = f(t + dt, y4)
        for i in range(N):
            y[i] += dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
            if NAMES[i] in POSITIVE and y[i] < 1e-3:
                y[i] = 1e-3
    return out, AL0


def lifetime_vi(AL, p=P):
    return 1.0 / (1.0 + math.exp(-(p["AVI"] + p["BVI"] * (AL - 26.0))))


# =========================================================================
#  CALIBRATION: exactly two constants (G0, KEXP) against exactly two targets
# =========================================================================
EMM = dict(AGE0=7.0, SER0=+0.50, NPAR=0, GRS=0.50, ETHN=0.0, OUTD=1.0, NEARD=3.0)
LAMP = dict(AGE0=9.7, SER0=-3.00, NPAR=2, GRS=0.70, ETHN=1.0, OUTD=1.0, NEARD=3.0)


def yr1(sc, p, key="AXCUM"):
    out, _ = run(dict(sc), years=1.0, dt=0.05, p=p, record_every=365.0)
    return out[-1]


def calibrate(verbose=True):
    p = dict(P)
    for it in range(40):
        # KEXP from the emmetrope, G0 from the LAMP placebo arm
        e = yr1(EMM, p)
        p["KEXP"] *= 0.100 / max(e["AXCUM"], 1e-6)
        l = yr1(LAMP, p)
        err = l["AXCUM"] - 0.410
        if abs(err) < 2e-5:
            break
        p["G0"] *= (1.0 - 0.55 * err / 0.410)
    e = yr1(EMM, p); l = yr1(LAMP, p)
    if verbose:
        print("CALIBRATION (2 constants, 2 targets, %d iterations)" % (it + 1))
        print("  G0   = %.5f     KEXP = %.6f mm/yr" % (p["G0"], p["KEXP"] * YR))
        print("  emmetrope age 7 : AL +%.4f mm/yr (target 0.100)   SER %+.3f D/yr"
              % (e["AXCUM"], e["SER"] - 0.50))
        print("  LAMP placebo    : AL +%.4f mm/yr (target 0.410)   SER %+.3f D/yr"
              " (target -0.81)" % (l["AXCUM"], l["SER"] + 3.00))
    return p


if __name__ == "__main__":
    p = calibrate()
    print()
    out, AL0 = run(dict(LAMP), years=3.0, p=p, record_every=182.5)
    print("LAMP placebo trajectory (AL0 = %.3f mm)" % AL0)
    print("%6s %8s %8s %8s %7s %7s %7s %7s %7s" %
          ("yr", "SER", "AL", "CHT", "CREEP", "MTR", "GAG", "COL1", "PLENS"))
    for r in out:
        print("%6.2f %8.3f %8.3f %8.1f %7.3f %7.3f %7.3f %7.3f %7.2f" %
              (r["yr"], r["SER"], r["AL"], r["CHT"], r["CREEPS"], r["MTR"],
               r["GAG"], r["COL1"], r["PLENS"]))
