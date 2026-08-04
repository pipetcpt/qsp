"""
cm_reference_model.py
HIV-associated Cryptococcal Meningitis (CM) QSP model
-- independent, dependency-free Python reference implementation of
   cm_mrgsolve_model.R.

WHY THIS FILE EXISTS
--------------------
Every number quoted in README.md is computed twice: once by mrgsolve
(C++/LSODA) and once here (pure-Python fixed-step RK4, standard library only).
The two implementations share no code, so agreement between them is evidence
that a result follows from the EQUATIONS rather than from one solver's
behaviour.  This file is also the tuning harness: `python3 cm_reference_model.py`
prints the full calibration report (early fungicidal activity per regimen,
intracranial-pressure trajectories, toxicity, 10-week mortality) against the
published trial values it was fitted to.

UNITS
-----
    time            days
    AmB / 5FC / FLU / dexamethasone / sertraline amounts   mg
    concentrations  ug/mL   (= mg/L)
    fungal burden   CFU/mL of CSF
    GXM antigen     ug/mL of CSF
    CSF pressure    mmH2O   (internally cmH2O; 1 cmH2O = 10 mmH2O)
    CSF volumes     mL, flows mL/day
    outflow resistance Rout  cmH2O per (mL/day)
                    (multiply by 1440 for the clinical cmH2O/(mL/min))

THE FOUR STRUCTURAL COMMITMENTS
-------------------------------
1.  Intracranial pressure is NOT a state variable.  It is the Davson residual
    of a volume budget:  ICP = Pss + Pel(Vex + EDEMA), with elastic pressure
    Pel absorbed at a rate Pel/Rout.  Pressure therefore cannot be "treated"
    directly; only formation, resistance or volume can be.

2.  Outflow resistance Rout is driven by CSF GXM (capsular polysaccharide
    obstructing the arachnoid outflow pathway), not by viable yeast.  GXM is
    released BOTH by living yeast (slow trickle) and by dying yeast (whole
    capsule, in a bolus).  Killing the fungus therefore adds antigen.  This
    single sign is the source of most of the model's non-obvious behaviour.

3.  The antigen pool clears ~20x more slowly than the viable pool.  Sterility
    and antigen-negativity are different events weeks apart, so a patient can
    be culture-negative and still be generating pressure and IRIS drive.

4.  Death is a hazard integral with separable terms (burden, pressure,
    perfusion, anaemia, neutropenia, hypokalaemia, IRIS).  Because the terms
    are separable the model can be asked WHICH clock killed a virtual patient,
    and different interventions answer differently.

USAGE
-----
    python3 cm_reference_model.py            # full calibration report
    python3 cm_reference_model.py --brief    # EFA + mortality table only

    import cm_reference_model as M
    r = M.run(M.regimen("ambition"), days=70)
    M.efa(r)                                 # log10 CFU/mL/day over days 1-14
"""

from math import exp, log, log10, sqrt

# ---------------------------------------------------------------------------
# 0.  State vector layout
# ---------------------------------------------------------------------------

SNAMES = [
    # --- amphotericin B ---
    "Ad",     # 0  AmB deoxycholate, central (mg)
    "Ad2",    # 1  AmB deoxycholate, peripheral tissue (mg)
    "Al",     # 2  liposomal AmB, central (mg)
    "Al2",    # 3  liposomal AmB, peripheral (mg)
    "Abr",    # 4  AmB at the CNS effect site (mg)
    "Akid",   # 5  AmB in renal cortex (mg)
    # --- flucytosine ---
    "FCg",    # 6  gut depot (mg)
    "FCc",    # 7  central (mg)
    "FCcsf",  # 8  CSF concentration (ug/mL)
    # --- fluconazole ---
    "FLg",    # 9  gut depot (mg)
    "FLc",    # 10 central (mg)
    "FLcsf",  # 11 CSF concentration (ug/mL)
    # --- adjuncts ---
    "DXg",    # 12 dexamethasone gut (mg)
    "DXc",    # 13 dexamethasone central (mg)
    "SRg",    # 14 sertraline gut (mg)
    "SRc",    # 15 sertraline central (mg)
    "SRbr",   # 16 sertraline brain total (mg)
    "IFNsc",  # 17 exogenous IFN-gamma SC depot (ug)
    # --- fungus ---
    "Fe",     # 18 extracellular, 5FC-susceptible (CFU/mL)
    "Fres",   # 19 extracellular, 5FC-resistant (CFU/mL)
    "Ft",     # 20 phenotypically tolerant / persister (CFU/mL)
    "Fi",     # 21 intracellular (macrophage-resident) (CFU/mL equiv)
    "Fp",     # 22 parenchymal / deep reservoir (CFU/mL equiv)
    "GXM",    # 23 CSF glucuronoxylomannan (ug/mL)
    # --- host immunity ---
    "CD4",    # 24 peripheral CD4 (cells/uL)
    "VL",     # 25 HIV RNA (log10 copies/mL)
    "MAC",    # 26 activated CNS macrophage/microglia index (1 = resting)
    "TH1",    # 27 CNS Th1 effector index
    "IFNG",   # 28 CSF IFN-gamma (pg/mL)
    "PROIN",  # 29 CSF proinflammatory index (TNF-a/IL-6 composite, pg/mL)
    "IL10",   # 30 CSF IL-10 (pg/mL)
    "WBC",    # 31 CSF leucocytes (cells/uL)
    "IMML",   # 32 lagged immune competence (dimensionless)
    "IRISa",  # 33 IRIS activity (dimensionless)
    # --- CNS hydrodynamics ---
    "Vex",    # 34 excess CSF volume (mL)
    "Rout",   # 35 CSF outflow resistance (cmH2O per mL/day)
    "EDEMA",  # 36 cerebral oedema, volume-equivalent (mL)
    "LEAK",   # 37 post-LP dural leak flux (mL/day)
    # --- safety ---
    "GFR",    # 38 (mL/min)
    "Kser",   # 39 serum potassium (mmol/L)
    "Hb",     # 40 haemoglobin (g/dL)
    "ANC",    # 41 absolute neutrophil count (1e9/L)
    "ALT",    # 42 (U/L)
    # --- outcome / bookkeeping ---
    "NEUR",   # 43 neuronal injury index (0-1)
    "HAZ",    # 44 cumulative death hazard
    "AICP",   # 45 integral of (ICP - 250)+ (mmH2O.day)
    "ABR",    # 46 integral of CNS AmB concentration (ug.day/mL)
    "AFC",    # 47 integral of CSF 5FC concentration (ug.day/mL)
    "CLPV",   # 48 cumulative CSF volume drained (mL)
    "CLPG",   # 49 cumulative GXM removed (ug)
    "KILLC",  # 50 cumulative yeast killed (CFU/mL equivalents)
    "ERG",    # 51 fungal membrane ergosterol, relative to untreated (1.0)
    "DIS",    # 52 permanent disability index: the time-integral of injury
]
IX = {n: i for i, n in enumerate(SNAMES)}
NS = len(SNAMES)

# ---------------------------------------------------------------------------
# 1.  Parameters
# ---------------------------------------------------------------------------

P = dict(
    # ---- body ----
    WT=60.0,            # kg
    MAP=85.0,           # mean arterial pressure, mmHg

    # ---- amphotericin B deoxycholate PK (2-cpt, mg / L / day) ----
    VcD=30.0, VpD=210.0, CLD=34.6, QD=40.0,
    # ---- liposomal amphotericin B PK (2-cpt) ----
    VcL=5.0, VpL=12.0, CLL=1.20, QL=3.0,
    # ---- CNS effect site ----
    Vbr=1.30,           # L, CNS distribution volume for AmB
    CLbrD=0.11250,       # L/day of plasma cleared into CNS, deoxycholate
    CLbrL=0.003090,      # L/day, liposomal (per unit PLASMA concentration)
    kbrOut=0.160,       # /day  CNS efflux (t1/2 = 4.3 d)
    # ---- renal cortex ----
    Vkid=0.30,          # L
    CLkdD=0.900,        # L/day into renal cortex, deoxycholate (free drug)
    CLkdL=0.0125,       # L/day, liposomal (liposome-shielded)
    kkdOut=0.055,       # /day

    # ---- flucytosine PK ----
    kaFC=30.0, FFC=0.87, VFC=45.0, CLFC0=156.0,
    fcCSF=0.74,         # CSF:plasma ratio at steady state
    keqFC=12.0,         # /day CSF equilibration

    # ---- fluconazole PK ----
    kaFL=18.0, FFL=0.90, VFL=45.0, CLFL=25.0,
    flCSF=0.80, keqFL=3.0,

    # ---- dexamethasone PK ----
    kaDX=25.0, FDX=0.80, VDX=60.0, CLDX=17.0, KDX=0.006,

    # ---- sertraline PK ----
    kaSR=15.0, FSR=0.44, VSR=1200.0, CLSR=1200.0,
    kSRbr=60.0, kSRbrOut=3.0,   # brain partition (Kp ~ 20)
    fuSRbr=0.015,               # free fraction in brain tissue

    # ---- exogenous IFN-gamma ----
    kaIFN=1.2,          # /day SC absorption
    FIFN=0.6,
    kIFNexo=48.0,       # pg/mL CSF per ug absorbed per day (lumped)

    # ---- fungal growth ----
    g0=0.350,           # /day intrinsic net growth in CSF
    Fmax=1.0e7,         # CFU/mL carrying capacity (nutrient limit)
    gp=0.160,           # /day parenchymal reservoir growth
    Fpmax=3.0e6,
    kSeed=0.0060,       # /day CSF -> parenchyma
    kShed=0.0006,       # /day parenchyma -> CSF

    # ---- amphotericin B pharmacodynamics ----
    # In the achievable CNS range the kill rate is nearly FIRST ORDER in
    # concentration (EC50 >> C).  That is not an aesthetic choice: it is what
    # the observed proportionality between AmB dose and EFA (0.7 vs 1.0 mg/kg,
    # Bicanic 2008) forces.  A saturating fit cannot reproduce a 1.43x dose
    # step producing a 1.37x rate step.  Flagged in README as the model's most
    # exposed PD assumption.
    KmaxA=9.80,         # /day maximum kill rate (~2.7 log10/day, in-vitro range)
    EC50A=2.40,         # ug/mL at the CNS effect site, ergosterol-replete
    hA=1.0,
    # ---- flucytosine PD ----
    KmaxF=0.72, EC50F=22.0, RESFAC=45.0,
    kMut=2.2e-8,        # /day  emergence of 5FC resistance
    # ---- fluconazole / ergosterol PD ----
    # Azole action is routed ENTIRELY through the ergosterol node, so the same
    # parameter that gives fluconazole its activity also costs amphotericin
    # its target.  This is what makes 5FC and FLU non-interchangeable partners.
    kERG=0.85,          # /day membrane ergosterol turnover
    IC50erg=16.0,       # ug/mL CSF fluconazole for 50% ergosterol depletion
    KmaxL=0.344,        # /day kill from ergosterol depletion itself
    gERG=1.0,           # exponent: growth rate scales as ERG^gERG
    aERG=0.5437,          # exponent: AmB EC50 scales as 1/ERG^aERG
    # ---- sertraline PD (from in-vitro MIC) ----
    KmaxSR=0.55, EC50SR=6.0,

    # ---- capsule-mediated drug tolerance ----
    kCap=1.10, KCap=300.0,      # EC50 multiplier from capsule/GXM
    tolPers=0.40,               # residual drug effect on persisters
    fIntra=0.22,                # AmB effect on intracellular yeast
    fParen=0.40,                # drug effect in parenchymal reservoir

    # ---- persister / intracellular exchange ----
    kPer=0.0035, kRev=0.045,
    kPhag=0.075, kEsc=0.085,

    # ---- immune killing ----
    kImm=0.0345, KIFN=180.0,
    kImmI=0.34,                 # intracellular killing by activated macrophage
    kImmP=0.100,

    # ---- GXM antigen ----
    sLive=5.00e-6,      # ug/mL/day per CFU/mL (live shedding)
    sLysis=4.50e-4,     # ug released per CFU killed (whole capsule)
    kGXM=0.055,         # /day clearance (t1/2 = 12.6 d)
    KGsat=800.0,        # ug/mL, saturation of capsule production/shedding

    # ---- host immunity ----
    KAg=400.0,          # ug/mL GXM for half-maximal antigen drive
    KCD4=110.0,         # cells/uL for half-maximal immune competence
    kMACon=1.55, kMACifn=0.30, KIFN2=250.0, kMACoff=0.16,
    kTHon=2.30, kTHoff=0.22,
    kIFNp=210.0, kIFNd=1.60,
    kPRp=95.0, kPRd=0.95, dexPR=0.80,
    kILp=26.0, kILd=0.85,
    kWBCp=0.1330, KWCD4=55.0, kWBCd=0.20, dexWBC=0.85,
    klag=0.055,         # /day lag of the immune-competence tracker
    kIRIS=2.00, KIR=200.0, kIRISoff=0.085,

    # ---- ART ----
    kCD4=0.0062, CD4t=290.0, kVL=0.34,

    # ---- CSF hydrodynamics ----
    Iform0=504.0,       # mL/day (0.35 mL/min)
    fSupp=0.45, KICPf=260.0,
    Pss=8.0,            # cmH2O sagittal sinus pressure
    Pel0=2.80,          # cmH2O elastic pressure at Vex+EDEMA = 0
    Eel=0.0347,         # /mL elastance coefficient
    Rout0=0.005556,     # cmH2O per mL/day  (= 8 cmH2O/(mL/min))
    aG=11.20, KG=60.0, # GXM contribution to Rout (multiplicative)
    aW=0.65, KW2=120.0, # CSF leucocyte contribution
    aI=1.55,            # IRIS contribution
    kRon=0.55, kRoff=0.10,
    kEon=1.35, KE=340.0, KE2=0.55, kEoff=0.17,
    kLeak=1.40,         # /day decay of post-LP dural leak
    Vcsf=150.0,         # mL total CSF (for fractional removal at LP)

    # ---- safety ----
    GFR0=92.0, kGFRrec=0.115, fGFRtox=0.62, KGFRtox=2.30,
    Kser0=4.05, kKrec=0.42, aKamb=1.35, KKamb=2.10,
    Hb0=9.60, kHbrec=0.048, sHbD=0.1150, sHbL=0.00105, sHbF=0.130,
    ANC0=3.10, kANCrec=0.135, sANC=2.50, KFCmar=110.0, hFCmar=1.8,
    ALT0=28.0, kALTrec=0.09, sALTfl=0.75, sALTfc=22.0,

    # ---- injury and hazard ----
    aNicp=0.0130, aNcpp=0.0180, aNpro=0.0100, aNbur=0.0062, aNiris=0.0290,
    aNdex=0.0020,       # steroid myopathy / encephalopathy contribution
    kNrep=0.0550, kDIS=0.0300,
    # ---- hazard thresholds ----
    ANCthr=1.5, Kthr=3.0, GFRthr=60.0, ICPthr=250.0,
    # Hazard coefficients.  TWO-STAGE CALIBRATION (see README):
    # the seven coefficients that induction-therapy trials cannot identify --
    # none of them randomised pressure management, potassium replacement or
    # transfusion -- are FIXED at mechanistically anchored values, and only
    # hB (burden) and hAMS (neurological injury) are fitted, by non-negative
    # least squares, to 17 mortality endpoints from 5 randomised trials.
    h0=0.000500,        # background hazard of advanced HIV
    hB=0.013906,        # per log10 CFU/mL above 3      [FITTED]
    hICP=0.030000,      # per unit of (ICP-250)/250      [fixed]
    hCPP=0.030000,      # per unit of (60-CPP)/20        [fixed, unidentified]
    hAN=0.000400,       # per g/dL fall from Hb0         [fixed]
    hNEUT=0.004000,     # per 1e9/L below ANCthr         [fixed, deliberately small]
    hK=0.001200,        # per mmol/L below Kthr          [fixed]
    hGFR=0.000400,      # per unit of (GFRthr-GFR)/30    [fixed]
    hIRIS=0.008000,     # per unit IRIS activity         [fixed, checked vs COAT]
    hSTER=0.002900,     # per day the CSF culture stays positive  [FITTED]
    hDEX=0.000400,      # steroid-attributable infection/metabolic harm [fixed]
    hAMS=0.004649,      # per unit NEUR x 10             [FITTED]
)

# ---------------------------------------------------------------------------
# 2.  Right-hand side
# ---------------------------------------------------------------------------


def derived(y, p):
    """Algebraic quantities, including ICP (a residual, never a state)."""
    d = {}
    d["Cd"] = y[IX["Ad"]] / p["VcD"]
    d["Cl"] = y[IX["Al"]] / p["VcL"]
    d["Cbr"] = y[IX["Abr"]] / p["Vbr"]
    d["Ckid"] = y[IX["Akid"]] / p["Vkid"]
    d["Cfc"] = y[IX["FCc"]] / p["VFC"]
    d["Cfl"] = y[IX["FLc"]] / p["VFL"]
    d["Cdx"] = y[IX["DXc"]] / p["VDX"]
    d["Csr"] = y[IX["SRc"]] / p["VSR"]
    d["CsrBrFree"] = p["fuSRbr"] * y[IX["SRbr"]] / p["Vbr"]

    # ---- CSF hydrodynamics: the Davson residual --------------------------
    Vload = y[IX["Vex"]] + y[IX["EDEMA"]]
    Pel = p["Pel0"] * exp(p["Eel"] * Vload)
    d["Pel"] = Pel
    d["ICPcm"] = p["Pss"] + Pel
    d["ICP"] = 10.0 * d["ICPcm"]                     # mmH2O
    d["ICPmmHg"] = d["ICP"] / 13.6
    d["CPP"] = p["MAP"] - d["ICPmmHg"]
    d["Iform"] = p["Iform0"] * (1.0 - p["fSupp"] * d["ICP"] / (d["ICP"] + p["KICPf"]))
    d["Iabs"] = Pel / y[IX["Rout"]]

    # ---- fungal burden as the laboratory sees it -------------------------
    d["CFU"] = max(y[IX["Fe"]] + y[IX["Fres"]] + y[IX["Ft"]] + y[IX["Fi"]], 1e-9)
    d["logCFU"] = log10(d["CFU"])
    # CrAg lateral-flow titre: 1:2^n, anchored so 1 ug/mL ~ 1:160
    d["CrAgTitre"] = 160.0 * max(y[IX["GXM"]], 1e-6)

    # ---- drug effects ----------------------------------------------------
    # THE PIVOT OF THE MODEL.  Amphotericin's apparent potency is divided by
    # the ergosterol content of the membrane it has to bind, and fluconazole
    # is the thing that lowers that content.  One node, two opposite signs.
    erg = max(y[IX["ERG"]], 0.12)
    d["erg"] = erg
    cap = 1.0 + p["kCap"] * y[IX["GXM"]] / (y[IX["GXM"]] + p["KCap"])
    d["cap"] = cap
    ec50a = p["EC50A"] * cap / erg ** p["aERG"]
    d["EC50Aeff"] = ec50a
    d["killA"] = (p["KmaxA"] * d["Cbr"] ** p["hA"]
                  / (ec50a ** p["hA"] + d["Cbr"] ** p["hA"]))
    fcsf = y[IX["FCcsf"]]
    d["killF"] = p["KmaxF"] * fcsf / (p["EC50F"] + fcsf)
    d["killFr"] = p["KmaxF"] * fcsf / (p["EC50F"] * p["RESFAC"] + fcsf)
    d["killL"] = p["KmaxL"] * (1.0 - erg)
    d["gsupp"] = erg ** p["gERG"]
    d["killSR"] = (p["KmaxSR"] * d["CsrBrFree"]
                   / (p["EC50SR"] + d["CsrBrFree"]))

    dex = d["Cdx"] / (d["Cdx"] + p["KDX"])
    d["dex"] = dex
    d["killI"] = (p["kImm"] * y[IX["MAC"]]
                  * (1.0 + y[IX["IFNG"]] / (y[IX["IFNG"]] + p["KIFN"]))
                  * (1.0 - 0.75 * dex))

    d["immcomp"] = y[IX["CD4"]] / (y[IX["CD4"]] + p["KCD4"])
    d["antigen"] = y[IX["GXM"]] / (y[IX["GXM"]] + p["KAg"])
    d["fcMarrow"] = (d["Cfc"] ** p["hFCmar"]
                     / (d["Cfc"] ** p["hFCmar"] + p["KFCmar"] ** p["hFCmar"]))
    return d


def rhs(t, y, p, inp, d=None):
    """inp: dict of exogenous infusion / absorption rates for this interval."""
    if d is None:
        d = derived(y, p)
    dy = [0.0] * NS

    # ================= amphotericin B ===================================
    Cd, Cl = d["Cd"], d["Cl"]
    Cp_d2 = y[IX["Ad2"]] / p["VpD"]
    Cp_l2 = y[IX["Al2"]] / p["VpL"]
    dy[IX["Ad"]] = (inp.get("rAmBd", 0.0) - p["CLD"] * Cd
                    - p["QD"] * (Cd - Cp_d2)
                    - p["CLbrD"] * Cd - p["CLkdD"] * Cd)
    dy[IX["Ad2"]] = p["QD"] * (Cd - Cp_d2)
    dy[IX["Al"]] = (inp.get("rAmBl", 0.0) - p["CLL"] * Cl
                    - p["QL"] * (Cl - Cp_l2)
                    - p["CLbrL"] * Cl - p["CLkdL"] * Cl)
    dy[IX["Al2"]] = p["QL"] * (Cl - Cp_l2)
    dy[IX["Abr"]] = (p["CLbrD"] * Cd + p["CLbrL"] * Cl
                     - p["kbrOut"] * y[IX["Abr"]])
    dy[IX["Akid"]] = (p["CLkdD"] * Cd + p["CLkdL"] * Cl
                      - p["kkdOut"] * y[IX["Akid"]])

    # ================= flucytosine =======================================
    # renal clearance scales with GFR -> AmB nephrotoxicity raises 5FC
    clfc = p["CLFC0"] * (0.12 + 0.88 * y[IX["GFR"]] / p["GFR0"])
    dy[IX["FCg"]] = inp.get("rFCg", 0.0) - p["kaFC"] * y[IX["FCg"]]
    dy[IX["FCc"]] = (p["FFC"] * p["kaFC"] * y[IX["FCg"]] - clfc * d["Cfc"])
    dy[IX["FCcsf"]] = p["keqFC"] * (p["fcCSF"] * d["Cfc"] - y[IX["FCcsf"]])

    # ================= fluconazole =======================================
    dy[IX["FLg"]] = inp.get("rFLg", 0.0) - p["kaFL"] * y[IX["FLg"]]
    dy[IX["FLc"]] = (p["FFL"] * p["kaFL"] * y[IX["FLg"]] - p["CLFL"] * d["Cfl"])
    dy[IX["FLcsf"]] = p["keqFL"] * (p["flCSF"] * d["Cfl"] - y[IX["FLcsf"]])

    # ================= adjuncts ==========================================
    dy[IX["DXg"]] = inp.get("rDXg", 0.0) - p["kaDX"] * y[IX["DXg"]]
    dy[IX["DXc"]] = p["FDX"] * p["kaDX"] * y[IX["DXg"]] - p["CLDX"] * d["Cdx"]
    dy[IX["SRg"]] = inp.get("rSRg", 0.0) - p["kaSR"] * y[IX["SRg"]]
    dy[IX["SRc"]] = (p["FSR"] * p["kaSR"] * y[IX["SRg"]] - p["CLSR"] * d["Csr"]
                     - p["kSRbr"] * d["Csr"] + p["kSRbrOut"] * y[IX["SRbr"]])
    dy[IX["SRbr"]] = p["kSRbr"] * d["Csr"] - p["kSRbrOut"] * y[IX["SRbr"]]
    dy[IX["IFNsc"]] = inp.get("rIFN", 0.0) - p["kaIFN"] * y[IX["IFNsc"]]

    # ================= ergosterol ========================================
    ergtar = 1.0 / (1.0 + y[IX["FLcsf"]] / p["IC50erg"])
    dy[IX["ERG"]] = p["kERG"] * (ergtar - y[IX["ERG"]])

    # ================= fungus ============================================
    Ftot = d["CFU"]
    gr = p["g0"] * d["gsupp"] * max(0.0, 1.0 - Ftot / p["Fmax"])
    Fe, Fres, Ft, Fi, Fp = (y[IX["Fe"]], y[IX["Fres"]], y[IX["Ft"]],
                            y[IX["Fi"]], y[IX["Fp"]])

    killE = d["killA"] + d["killF"] + d["killL"] + d["killI"] + d["killSR"]
    killR = d["killA"] + d["killFr"] + d["killL"] + d["killI"] + d["killSR"]
    killT = (p["tolPers"] * (d["killA"] + d["killF"] + d["killL"] + d["killSR"])
             + d["killI"])
    killIn = (p["fIntra"] * d["killA"] + d["killF"] + p["kImmI"] * y[IX["MAC"]]
              * (1.0 - 0.75 * d["dex"]))
    killP = (p["fParen"] * (d["killA"] + d["killF"] + d["killL"])
             + p["kImmP"] * y[IX["MAC"]] * (1.0 - 0.75 * d["dex"]))

    phag = p["kPhag"] * y[IX["MAC"]]
    dy[IX["Fe"]] = (gr * Fe - killE * Fe - phag * Fe + p["kEsc"] * Fi
                    - p["kPer"] * Fe + p["kRev"] * Ft
                    - p["kMut"] * Fe - p["kSeed"] * Fe + p["kShed"] * Fp)
    dy[IX["Fres"]] = (gr * Fres - killR * Fres - phag * Fres + p["kMut"] * Fe
                      - p["kSeed"] * Fres)
    dy[IX["Ft"]] = (0.10 * gr * Ft + p["kPer"] * Fe - p["kRev"] * Ft
                    - killT * Ft)
    dy[IX["Fi"]] = (0.30 * gr * Fi + phag * (Fe + Fres) - p["kEsc"] * Fi
                    - killIn * Fi)
    dy[IX["Fp"]] = (p["gp"] * Fp * max(0.0, 1.0 - Fp / p["Fpmax"])
                    + p["kSeed"] * (Fe + Fres) - p["kShed"] * Fp - killP * Fp)

    killFlux = (killE * Fe + killR * Fres + killT * Ft + killIn * Fi
                + killP * Fp)
    dy[IX["KILLC"]] = killFlux
    gsat = 1.0 / (1.0 + y[IX["GXM"]] / p["KGsat"])
    dy[IX["GXM"]] = ((p["sLive"] * (Fe + Fres + Ft + Fi + Fp)
                      + p["sLysis"] * killFlux) * gsat
                     - p["kGXM"] * y[IX["GXM"]])

    # ================= host immunity =====================================
    art = inp.get("art", 0.0)
    dy[IX["CD4"]] = art * p["kCD4"] * (p["CD4t"] - y[IX["CD4"]])
    dy[IX["VL"]] = -art * p["kVL"] * y[IX["VL"]]

    dy[IX["MAC"]] = (p["kMACon"] * d["antigen"] * (0.30 + 0.70 * d["immcomp"])
                     * (1.0 - 0.80 * d["dex"])
                     + p["kMACifn"] * y[IX["IFNG"]] / (y[IX["IFNG"]] + p["KIFN2"])
                     - p["kMACoff"] * (y[IX["MAC"]] - 1.0))
    dy[IX["TH1"]] = (p["kTHon"] * d["antigen"] * d["immcomp"]
                     * (1.0 - 0.90 * d["dex"]) - p["kTHoff"] * y[IX["TH1"]])
    dy[IX["IFNG"]] = (p["kIFNp"] * y[IX["TH1"]]
                      + p["kIFNexo"] * p["FIFN"] * p["kaIFN"] * y[IX["IFNsc"]]
                      - p["kIFNd"] * y[IX["IFNG"]])
    dy[IX["PROIN"]] = (p["kPRp"] * (y[IX["MAC"]] * d["antigen"]
                                    + 3.0 * y[IX["IRISa"]])
                       - p["kPRd"] * y[IX["PROIN"]]
                       - p["dexPR"] * d["dex"] * y[IX["PROIN"]])
    dy[IX["IL10"]] = p["kILp"] * y[IX["MAC"]] * d["antigen"] - p["kILd"] * y[IX["IL10"]]
    dy[IX["WBC"]] = (p["kWBCp"] * y[IX["PROIN"]]
                     * (y[IX["CD4"]] / (y[IX["CD4"]] + p["KWCD4"]))
                     * (1.0 - p["dexWBC"] * d["dex"])
                     - p["kWBCd"] * y[IX["WBC"]])

    dy[IX["IMML"]] = p["klag"] * (d["immcomp"] - y[IX["IMML"]])
    recov = max(0.0, d["immcomp"] - y[IX["IMML"]])   # proportional to dCD4/dt
    dy[IX["IRISa"]] = (p["kIRIS"] * recov
                       * y[IX["GXM"]] / (y[IX["GXM"]] + p["KIR"])
                       * (1.0 - 0.70 * d["dex"])
                       - p["kIRISoff"] * y[IX["IRISa"]])

    # ================= CNS hydrodynamics =================================
    dy[IX["Vex"]] = d["Iform"] - d["Iabs"] - y[IX["LEAK"]]
    dy[IX["LEAK"]] = -p["kLeak"] * y[IX["LEAK"]]
    Rtar = p["Rout0"] * (1.0
                         + p["aG"] * y[IX["GXM"]] / (y[IX["GXM"]] + p["KG"])
                         + p["aW"] * y[IX["WBC"]] / (y[IX["WBC"]] + p["KW2"])
                         + p["aI"] * y[IX["IRISa"]])
    kR = p["kRon"] if Rtar > y[IX["Rout"]] else p["kRoff"]
    dy[IX["Rout"]] = kR * (Rtar - y[IX["Rout"]])
    dy[IX["EDEMA"]] = (p["kEon"] * (y[IX["PROIN"]] / (y[IX["PROIN"]] + p["KE"])
                                    + 2.0 * y[IX["IRISa"]]
                                    / (y[IX["IRISa"]] + p["KE2"]))
                       - p["kEoff"] * y[IX["EDEMA"]])

    # ================= safety ============================================
    ck = d["Ckid"]
    gtar = p["GFR0"] * (1.0 - p["fGFRtox"] * ck / (ck + p["KGFRtox"]))
    dy[IX["GFR"]] = p["kGFRrec"] * (gtar - y[IX["GFR"]])
    ktar = (p["Kser0"] - p["aKamb"] * ck / (ck + p["KKamb"])
            + inp.get("ksupp", 0.0))
    dy[IX["Kser"]] = p["kKrec"] * (ktar - y[IX["Kser"]])
    dy[IX["Hb"]] = (p["kHbrec"] * (p["Hb0"] - y[IX["Hb"]])
                    - p["sHbD"] * Cd - p["sHbL"] * Cl
                    - p["sHbF"] * d["fcMarrow"])
    # marrow suppression cannot remove precursors that are already gone
    dy[IX["ANC"]] = (p["kANCrec"] * (p["ANC0"] - y[IX["ANC"]])
                     - p["sANC"] * d["fcMarrow"]
                     * y[IX["ANC"]] / (y[IX["ANC"]] + 0.30))
    dy[IX["ALT"]] = (p["kALTrec"] * (p["ALT0"] - y[IX["ALT"]])
                     + p["sALTfl"] * d["Cfl"] + p["sALTfc"] * d["fcMarrow"])

    # ================= injury and hazard =================================
    icpX = max(0.0, d["ICP"] - 250.0) / 100.0
    cppX = max(0.0, 60.0 - d["CPP"]) / 10.0
    burX = max(0.0, d["logCFU"] - 3.0)
    # NEUR is reversible injury; DIS is the fraction of it that does not
    # recover, and is the model's analogue of a 10-week disability endpoint.
    dy[IX["DIS"]] = p["kDIS"] * y[IX["NEUR"]] * (1.0 - y[IX["DIS"]])
    dy[IX["NEUR"]] = ((p["aNicp"] * icpX + p["aNcpp"] * cppX
                       + p["aNpro"] * y[IX["PROIN"]] / (y[IX["PROIN"]] + 400.0)
                       + p["aNbur"] * burX + p["aNiris"] * y[IX["IRISa"]]
                       + p["aNdex"] * d["dex"])
                      * (1.0 - y[IX["NEUR"]]) - p["kNrep"] * y[IX["NEUR"]])

    haz = (p["h0"]
           + p["hB"] * burX
           + p["hICP"] * max(0.0, d["ICP"] - 250.0) / 250.0
           + p["hCPP"] * max(0.0, 60.0 - d["CPP"]) / 20.0
           + p["hAN"] * max(0.0, p["Hb0"] - y[IX["Hb"]])
           + p["hNEUT"] * max(0.0, p["ANCthr"] - y[IX["ANC"]])
           + p["hK"] * max(0.0, p["Kthr"] - y[IX["Kser"]])
           + p["hGFR"] * max(0.0, p["GFRthr"] - y[IX["GFR"]]) / 30.0
           + p["hIRIS"] * y[IX["IRISa"]]
           + (p["hSTER"] if d["CFU"] >= 1.0 else 0.0)
           + p["hDEX"] * d["dex"]
           + p["hAMS"] * y[IX["NEUR"]] * 10.0)
    dy[IX["HAZ"]] = haz
    dy[IX["AICP"]] = max(0.0, d["ICP"] - 250.0)
    dy[IX["ABR"]] = d["Cbr"]
    dy[IX["AFC"]] = y[IX["FCcsf"]]
    return dy


# ---------------------------------------------------------------------------
# 3.  Integrator with events
# ---------------------------------------------------------------------------


def rk4(y, t, h, p, inp):
    k1 = rhs(t, y, p, inp)
    y2 = [y[i] + 0.5 * h * k1[i] for i in range(NS)]
    k2 = rhs(t + 0.5 * h, y2, p, inp)
    y3 = [y[i] + 0.5 * h * k2[i] for i in range(NS)]
    k3 = rhs(t + 0.5 * h, y3, p, inp)
    y4 = [y[i] + h * k3[i] for i in range(NS)]
    k4 = rhs(t + h, y4, p, inp)
    return [y[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
            for i in range(NS)]


def lumbar_puncture(y, p, target_frac=0.50, floor_cm=12.0, max_ml=30.0):
    """Therapeutic LP: drain until closing pressure is the lower of 50% of
    opening pressure and 120 mmH2O, capped at 30 mL.  Removes CSF-resident
    yeast, antigen and leucocytes in proportion to the volume removed."""
    d = derived(y, p)
    Pel = d["Pel"]
    Pel_t = max(target_frac * Pel, floor_cm - p["Pss"], 0.35)
    if Pel_t >= Pel:
        return y, 0.0, 0.0
    dV = min(max_ml, log(Pel / Pel_t) / p["Eel"])
    y = list(y)
    # take the volume out of the free CSF space, not the oedema
    take = min(dV, y[IX["Vex"]] + 40.0)
    y[IX["Vex"]] -= take
    fr = min(0.28, dV / p["Vcsf"])
    gxm_removed = fr * y[IX["GXM"]] * p["Vcsf"]
    for s in ("GXM", "Fe", "Fres", "Ft", "Fi", "WBC"):
        y[IX[s]] *= (1.0 - fr)
    y[IX["LEAK"]] += p["kLeak"] * 190.0        # transient dural leak
    y[IX["CLPV"]] += dV
    y[IX["CLPG"]] += gxm_removed
    return y, dV, gxm_removed


def simulate(y0, schedule, days, p=None, dt=0.0025, record=1.0, lp_times=()):
    """schedule: list of (t_start, t_end, inp_dict) intervals, non-overlapping.
    Returns dict of time series."""
    p = p or P
    y = list(y0)
    # build the break-point grid
    bps = {0.0, float(days)}
    for (a, b, _) in schedule:
        bps.add(float(a)); bps.add(min(float(b), days))
    for tl in lp_times:
        if tl < days:
            bps.add(float(tl))
    grid = sorted(x for x in bps if 0.0 <= x <= days)
    lp_set = sorted(float(x) for x in lp_times if x < days)

    out = {"t": [], "y": []}
    nextrec = 0.0
    t = 0.0

    def cur_inp(tt):
        inp = {}
        for (a, b, d_) in schedule:
            if a - 1e-12 <= tt < b - 1e-12:
                for k, v in d_.items():
                    inp[k] = inp.get(k, 0.0) + v
        return inp

    def snap(tt, yy):
        out["t"].append(tt); out["y"].append(list(yy))

    snap(0.0, y)
    nextrec = record
    for gi in range(len(grid) - 1):
        t0, t1 = grid[gi], grid[gi + 1]
        for tl in lp_set:
            if abs(tl - t0) < 1e-9:
                y, _, _ = lumbar_puncture(y, p)
        inp = cur_inp(0.5 * (t0 + t1))
        n = max(1, int(round((t1 - t0) / dt)))
        h = (t1 - t0) / n
        t = t0
        for _ in range(n):
            y = rk4(y, t, h, p, inp)
            t += h
            if t >= nextrec - 1e-9:
                snap(t, y); nextrec += record
    if out["t"][-1] < days - 1e-9:
        snap(t, y)
    return out


# ---------------------------------------------------------------------------
# 4.  Baseline (the model's own quasi-equilibrium at presentation)
# ---------------------------------------------------------------------------


def baseline(p=None, logCFU=5.0, CD4=25.0, VL=5.4, burn=70.0):
    """Relax the SLOW variables (antigen, resistance, immunity, oedema,
    haematology) to quasi-steady state with the fungal burden CLAMPED at the
    presenting value, then release.  The presenting state is therefore the
    model's own equilibrium given the burden, not a hand-set vector."""
    p = p or P
    y = [0.0] * NS
    F = 10.0 ** logCFU
    y[IX["Fe"]] = F * 0.9930
    y[IX["Ft"]] = F * 0.0020
    y[IX["Fi"]] = F * 0.0050
    y[IX["Fres"]] = F * 1e-6
    y[IX["Fp"]] = F * 0.15
    y[IX["GXM"]] = 40.0
    y[IX["ERG"]] = 1.0
    y[IX["CD4"]] = CD4
    y[IX["VL"]] = VL
    y[IX["MAC"]] = 1.0
    y[IX["TH1"]] = 0.05
    y[IX["IFNG"]] = 5.0
    y[IX["PROIN"]] = 40.0
    y[IX["IL10"]] = 10.0
    y[IX["WBC"]] = 2.0
    y[IX["IMML"]] = CD4 / (CD4 + p["KCD4"])
    y[IX["IRISa"]] = 0.0
    y[IX["Vex"]] = 0.0
    y[IX["Rout"]] = p["Rout0"]
    y[IX["EDEMA"]] = 0.0
    y[IX["GFR"]] = p["GFR0"]
    y[IX["Kser"]] = p["Kser0"]
    y[IX["Hb"]] = p["Hb0"]
    y[IX["ANC"]] = p["ANC0"]
    y[IX["ALT"]] = p["ALT0"]

    clamp = [IX["Fe"], IX["Fres"], IX["Ft"], IX["Fi"], IX["Fp"],
             IX["CD4"], IX["VL"], IX["NEUR"], IX["HAZ"],
             IX["AICP"], IX["ABR"], IX["AFC"], IX["CLPV"], IX["CLPG"],
             IX["KILLC"]]
    keep = {i: y[i] for i in clamp}
    h, t = 0.005, 0.0
    n = int(burn / h)
    for _ in range(n):
        y = rk4(y, t, h, p, {})
        t += h
        for i, v in keep.items():
            y[i] = v
    for i in (IX["NEUR"], IX["DIS"], IX["HAZ"], IX["AICP"], IX["ABR"], IX["AFC"],
              IX["CLPV"], IX["CLPG"], IX["KILLC"]):
        y[i] = 0.0
    return y


# ---------------------------------------------------------------------------
# 5.  Regimens
# ---------------------------------------------------------------------------

def _daily(t0, t1, rate_key, mg_per_day):
    """Continuous-rate approximation of a once-daily oral/IV course."""
    return [(t0, t1, {rate_key: mg_per_day})]


def _amb_d(t0, ndays, mgkg, p=None):
    """AmB deoxycholate, 4-hour infusion once daily."""
    p = p or P
    dose = mgkg * p["WT"]
    inf = 4.0 / 24.0
    out = []
    for k in range(ndays):
        out.append((t0 + k, t0 + k + inf, {"rAmBd": dose / inf}))
    return out


def _amb_l(t0, mgkg, p=None, hours=2.0):
    p = p or P
    dose = mgkg * p["WT"]
    inf = hours / 24.0
    return [(t0, t0 + inf, {"rAmBl": dose / inf})]


def _antifungal(name, p, days):
    """Antifungal / adjunct intervals ONLY.  ART and electrolyte replacement
    are attached exactly once, by regimen(), so that composite regimens
    ("ambition_dex" etc.) cannot double-count them."""
    wt = p["WT"]
    fc = 100.0 * wt                 # flucytosine 100 mg/kg/day
    S = []
    lps = []

    if name == "untreated":
        lab = "No antifungal therapy (natural history)"

    elif name == "flu800":
        S += _daily(0, 14, "rFLg", 800.0) + _daily(14, days, "rFLg", 400.0)
        lab = "Fluconazole 800 mg/d x14d, then 400 mg/d"

    elif name == "flu1200":
        S += _daily(0, 14, "rFLg", 1200.0) + _daily(14, days, "rFLg", 400.0)
        lab = "Fluconazole 1200 mg/d x14d, then 400 mg/d"

    elif name == "oral_acta":
        S += _daily(0, 14, "rFLg", 1200.0) + _daily(0, 14, "rFCg", fc)
        S += _daily(14, days, "rFLg", 800.0)
        lab = "ACTA oral arm: FLU 1200 + 5FC 100 mg/kg x14d"

    elif name == "ambd_alone":
        S += _amb_d(0, 28, 1.0, p) + _daily(28, days, "rFLg", 400.0)
        lab = "AmB-d 1 mg/kg x28d (Day 2013 group 1)"

    elif name == "ambd_flu":
        S += _amb_d(0, 14, 1.0, p) + _daily(0, 14, "rFLg", 800.0)
        S += _daily(14, days, "rFLg", 400.0)
        lab = "AmB-d 1 mg/kg + FLU 800 x14d (Day 2013 group 3)"

    elif name == "ambd_5fc07":
        S += _amb_d(0, 14, 0.7, p) + _daily(0, 14, "rFCg", fc)
        S += _daily(14, days, "rFLg", 400.0)
        lab = "AmB-d 0.7 mg/kg + 5FC x14d (Bicanic 2008 group 1)"

    elif name == "ambd_5fc":
        S += _amb_d(0, 14, 1.0, p) + _daily(0, 14, "rFCg", fc)
        S += _daily(14, days, "rFLg", 400.0)
        lab = "AmB-d 1 mg/kg + 5FC x14d (Day 2013 / Bicanic group 2)"

    elif name == "ambd_5fc_1wk":
        S += _amb_d(0, 7, 1.0, p) + _daily(0, 7, "rFCg", fc)
        S += _daily(7, 14, "rFLg", 1200.0) + _daily(14, days, "rFLg", 800.0)
        lab = "AmB-d 1 mg/kg + 5FC x7d, then FLU 1200 x7d (ACTA/AMBITION control)"

    elif name == "ambition":
        S += _amb_l(0, 10.0, p) + _daily(0, 14, "rFCg", fc)
        S += _daily(0, 14, "rFLg", 1200.0) + _daily(14, days, "rFLg", 800.0)
        lab = "AMBITION: L-AmB 10 mg/kg single dose + 5FC + FLU 1200 x14d"

    elif name == "ambition_lp":
        S, lps, _ = _antifungal("ambition", p, days)
        lps = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 13.0]
        lab = "AMBITION + daily/alternate therapeutic LP x7 (days 1-13)"

    elif name == "ambition_dex":
        S, lps, _ = _antifungal("ambition", p, days)
        # CryptoDex taper: 0.3 mg/kg/d x1wk, 0.2, 0.1, 0.05, then oral taper
        tap = [(0, 7, 0.30), (7, 14, 0.20), (14, 21, 0.10), (21, 28, 0.05),
               (28, 35, 0.03), (35, 42, 0.015)]
        for a, b, mgkg in tap:
            S += _daily(a, b, "rDXg", mgkg * wt)
        lab = "AMBITION + adjunctive dexamethasone (CryptoDex taper)"

    elif name == "ambition_sert":
        S, lps, _ = _antifungal("ambition", p, days)
        S += _daily(0, 14, "rSRg", 400.0) + _daily(14, days, "rSRg", 200.0)
        lab = "AMBITION + sertraline 400 mg/d (ASTRO-CM)"

    elif name == "ambition_ifng":
        S, lps, _ = _antifungal("ambition", p, days)
        for k in (1.0, 3.0, 5.0, 8.0, 10.0, 12.0):
            S.append((k, k + 0.25, {"rIFN": 100.0 / 0.25}))   # 100 ug SC
        lab = "AMBITION + IFN-gamma 100 ug SC x6 doses"

    elif name.startswith("coat"):
        # COAT used AmB-d 0.7-1.0 + fluconazole 800 for 14 days
        S += _amb_d(0, 14, 0.9, p) + _daily(0, 14, "rFLg", 800.0)
        S += _daily(14, days, "rFLg", 400.0)
        lab = "COAT backbone (AmB-d 0.9 + FLU 800 x14d)"

    else:
        raise ValueError(f"unknown regimen {name!r}")

    return S, lps, lab


def regimen(name, p=None, art_day=35.0, lp=(), days=70.0):
    """Return (schedule, lp_times, label)."""
    p = p or P
    S, lps, lab = _antifungal(name, p, days)
    lps = list(lp) + list(lps)
    if art_day is not None and art_day < days:
        S.append((art_day, days, {"art": 1.0}))
        lab += f", ART day {art_day:g}"
    # routine potassium/magnesium replacement whenever amphotericin is given
    if any("rAmBd" in d_ or "rAmBl" in d_ for (_, _, d_) in S):
        S.append((0.0, min(21.0, days), {"ksupp": 0.55}))
    return S, lps, lab


def run(reg, days=70.0, p=None, y0=None, record=0.25, art_day=35.0):
    p = p or P
    if isinstance(reg, str):
        S, lps, lab = regimen(reg, p, art_day=art_day, days=days)
    else:
        S, lps, lab = reg
    y0 = y0 if y0 is not None else baseline(p)
    r = simulate(y0, S, days, p=p, record=record, lp_times=lps)
    r["label"] = lab
    r["p"] = p
    return r


# ---------------------------------------------------------------------------
# 6.  Read-outs
# ---------------------------------------------------------------------------


def get(r, name):
    i = IX[name]
    return [row[i] for row in r["y"]]


def series(r, key):
    """Derived (algebraic) series."""
    p = r["p"]
    return [derived(row, p)[key] for row in r["y"]]


def at(r, t, key):
    """Value of a state or derived quantity at time t (nearest grid point)."""
    ts = r["t"]
    j = min(range(len(ts)), key=lambda i: abs(ts[i] - t))
    if key in IX:
        return r["y"][j][IX[key]]
    return derived(r["y"][j], r["p"])[key]


def efa(r, t0=0.0, t1=14.0):
    """Early fungicidal activity: least-squares slope of log10 CSF CFU/mL."""
    ts, xs, ys = r["t"], [], []
    for j, t in enumerate(ts):
        if t0 - 1e-9 <= t <= t1 + 1e-9:
            xs.append(t)
            ys.append(derived(r["y"][j], r["p"])["logCFU"])
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    den = sum((xs[i] - mx) ** 2 for i in range(n))
    return num / den


def mortality(r, t):
    return 1.0 - exp(-at(r, t, "HAZ"))


def sterile_day(r, thresh=0.0):
    """First day on which the CSF culture is negative (< 1 CFU/mL)."""
    for j, t in enumerate(r["t"]):
        if derived(r["y"][j], r["p"])["logCFU"] < thresh:
            return t
    return None


def peak(r, key, t0=0.0, t1=1e9):
    vals = [(at_v, t) for t, at_v in
            zip(r["t"], (derived(row, r["p"])[key] if key not in IX
                         else row[IX[key]] for row in r["y"]))
            if t0 <= t <= t1]
    return max(vals)


def nadir(r, key, t0=0.0, t1=1e9):
    vals = [(at_v, t) for t, at_v in
            zip(r["t"], (derived(row, r["p"])[key] if key not in IX
                         else row[IX[key]] for row in r["y"]))
            if t0 <= t <= t1]
    return min(vals)


# ---------------------------------------------------------------------------
# 7.  Calibration report
# ---------------------------------------------------------------------------

EFA_TARGETS = [
    ("flu800",       -0.07, "Longley 2008 CID (FLU 800 mg)"),
    ("flu1200",      -0.11, "Nussbaum 2010 CID (FLU 1200 mg)"),
    ("oral_acta",    -0.28, "Nussbaum 2010 CID (FLU 1200 + 5FC)"),
    ("ambd_alone",   -0.31, "Day 2013 NEJM group 1"),
    ("ambd_flu",     -0.32, "Day 2013 NEJM group 3"),
    ("ambd_5fc07",   -0.40, "Bicanic 2008 ratio 0.80 x the 1 mg/kg arm"),
    ("ambd_5fc",     -0.49, "Day 2013 (-0.42) / Bicanic 2008 (-0.56)"),
    ("ambd_5fc_1wk", -0.42, "AMBITION control arm"),
    ("ambition",     -0.40, "AMBITION L-AmB arm"),
]

MORT_TARGETS = {
    "flu800":       (None, 0.600, "fluconazole 800 monotherapy, historical"),
    "flu1200":      (0.300, 0.560, "Nussbaum 2010 FLU-alone arm"),
    "oral_acta":    (0.182, 0.351, "ACTA oral arm"),
    "ambd_alone":   (0.248, 0.436, "Day 2013 group 1 (25/101, 44/101)"),
    "ambd_flu":     (0.197, 0.326, "Day 2013 group 3 (HR 0.78 / 0.71)"),
    "ambd_5fc":     (0.150, 0.300, "Day 2013 group 2 (15/100, 30/100)"),
    "ambd_5fc_1wk": (0.155, 0.270, "ACTA 24.2% / AMBITION control 28.7%"),
    "ambition":     (0.130, 0.248, "AMBITION L-AmB arm"),
}

# Virtual patients.  Each is a DIFFERENT quasi-equilibrium of the same model,
# reached by relaxing the slow variables at a different presenting burden and
# CD4 count -- not a hand-edited state vector.
PHENOTYPES = {
    "median":      dict(logCFU=5.0, CD4=25.0,
                        note="median trial participant"),
    "high_burden": dict(logCFU=5.9, CD4=25.0,
                        note="upper-quartile burden -> the high-pressure patient"),
    "paucicellular": dict(logCFU=5.3, CD4=12.0,
                        note="CD4 12, CSF WBC <5: the COAT high-risk subgroup"),
    "immune_partial": dict(logCFU=4.3, CD4=150.0,
                        note="partially preserved immunity"),
}


def report(brief=False):
    p = P
    y0 = baseline(p)
    d0 = derived(y0, p)
    print("=" * 78)
    print("HIV-ASSOCIATED CRYPTOCOCCAL MENINGITIS -- QSP REFERENCE MODEL")
    print("=" * 78)
    print(f"states: {NS}   parameters: {len(P)}")
    print("\n-- Presenting state (model's own quasi-equilibrium) ------------")
    print(f"  CSF fungal burden      {d0['logCFU']:.2f} log10 CFU/mL")
    print(f"  CSF GXM antigen        {y0[IX['GXM']]:.1f} ug/mL "
          f"(CrAg ~ 1:{int(d0['CrAgTitre']):,})")
    print(f"  CSF opening pressure   {d0['ICP']:.0f} mmH2O "
          f"({d0['ICPcm']:.1f} cmH2O)")
    print(f"  outflow resistance     {y0[IX['Rout']]*1440:.1f} "
          f"cmH2O/(mL/min)  [normal 6-10]")
    print(f"  CSF leucocytes         {y0[IX['WBC']]:.1f} /uL")
    print(f"  CSF IFN-gamma          {y0[IX['IFNG']]:.0f} pg/mL")
    print(f"  CD4                    {y0[IX['CD4']]:.0f} /uL   "
          f"HIV VL {y0[IX['VL']]:.1f} log10")
    print(f"  Hb {y0[IX['Hb']]:.1f} g/dL   ANC {y0[IX['ANC']]:.2f}   "
          f"GFR {y0[IX['GFR']]:.0f}   K {y0[IX['Kser']]:.2f}")

    runs = {}
    print("\n-- Early fungicidal activity, days 0-14 (log10 CFU/mL/day) ----")
    print(f"  {'regimen':16s} {'model':>7s} {'target':>7s}  source")
    for name, tgt, src in EFA_TARGETS:
        r = run(name, days=70.0, p=p, y0=y0)
        runs[name] = r
        e = efa(r)
        flag = "" if abs(e - tgt) < 0.055 else "  <-- off"
        print(f"  {name:16s} {e:7.3f} {tgt:7.2f}  {src}{flag}")

    for extra in ("untreated", "ambition_lp", "ambition_dex",
                  "ambition_sert", "ambition_ifng"):
        runs[extra] = run(extra, days=70.0, p=p, y0=y0)

    print("\n-- Mortality (cumulative hazard -> 1-exp(-H)) -----------------")
    print(f"  {'regimen':16s} {'2wk':>7s} {'tgt':>6s} {'10wk':>7s} {'tgt':>6s}")
    for name, (m2t, m10t, src) in MORT_TARGETS.items():
        r = runs[name]
        m2, m10 = mortality(r, 14.0), mortality(r, 70.0)
        s2 = f"{m2t*100:.1f}" if m2t else "  -"
        print(f"  {name:16s} {m2*100:6.1f}% {s2:>6s} "
              f"{m10*100:6.1f}% {m10t*100:5.1f}   {src}")
    print(f"  {'untreated':16s} {mortality(runs['untreated'],14)*100:6.1f}%"
          f" {'  -':>6s} {mortality(runs['untreated'],70)*100:6.1f}%"
          f"    -   (no antifungal)")

    if brief:
        return runs

    print("\n-- Thesis A: two clocks. Burden vs pressure ------------------")
    print(f"  {'regimen':16s} {'EFA':>6s} {'d14 log':>8s} {'ICP pk':>7s}"
          f" {'day':>5s} {'AUC>250':>8s} {'GXM pk':>7s}")
    for name in ("untreated", "flu1200", "oral_acta", "ambd_5fc",
                 "ambd_5fc_1wk", "ambition", "ambition_lp"):
        r = runs[name]
        ip, ipt = peak(r, "ICP", 0, 42)
        print(f"  {name:16s} {efa(r):6.2f} {at(r,14,'logCFU'):8.2f}"
              f" {ip:7.0f} {ipt:5.1f} {at(r,42,'AICP'):8.0f}"
              f" {peak(r,'GXM',0,42)[0]:7.0f}")

    print("\n  Hazard decomposition at 6 weeks (share of cumulative hazard):")
    for name in ("ambd_5fc", "ambition", "ambition_lp"):
        sh = hazard_shares(runs[name], 42.0)
        tot = sum(sh.values())
        top = sorted(sh.items(), key=lambda kv: -kv[1])[:5]
        s = "  ".join(f"{k} {100*v/tot:.0f}%" for k, v in top)
        print(f"    {name:16s} {s}")

    print("\n-- Thesis B: L-AmB 10 mg/kg vs 7 x AmB-d 1 mg/kg -------------")
    for name in ("ambd_5fc_1wk", "ambition"):
        r = runs[name]
        print(f"  {name:16s} brain AUC0-14 {at(r,14,'ABR'):6.1f} ug.d/mL   "
              f"brain Cmax {peak(r,'Cbr',0,70)[0]:5.2f} ug/mL @d"
              f"{peak(r,'Cbr',0,70)[1]:.1f}   "
              f"T>0.15 ug/mL {t_above(r,'Cbr',0.15):.1f} d")
        print(f"  {'':16s} renal AmB AUC0-14 {trapz(r,'Ckid',0,14):6.1f} "
              f"ug.d/mL   GFR nadir {nadir(r,'GFR',0,42)[0]:.0f}   "
              f"K nadir {nadir(r,'Kser',0,42)[0]:.2f}   "
              f"Hb nadir {nadir(r,'Hb',0,42)[0]:.2f}   "
              f"ANC nadir {nadir(r,'ANC',0,42)[0]:.2f}")

    print("\n-- Thesis C: the flucytosine-clearance trap ------------------")
    for name in ("oral_acta", "ambd_5fc_1wk", "ambd_5fc"):
        r = runs[name]
        print(f"  {name:16s} 5FC CSF AUC0-14 {at(r,14,'AFC'):6.0f} ug.d/mL"
              f"   GFR nadir {nadir(r,'GFR',0,28)[0]:5.0f}"
              f"   peak 5FC plasma {peak(r,'Cfc',0,28)[0]:5.1f} ug/mL"
              f"   ANC nadir {nadir(r,'ANC',0,28)[0]:.2f}")

    print("\n-- Thesis D: dexamethasone --------------------------------------")
    a, dx = runs["ambition"], runs["ambition_dex"]
    print(f"  EFA            {efa(a):6.3f}  ->  {efa(dx):6.3f}"
          f"   (CryptoDex: clearance slower on dexamethasone)")
    print(f"  ICP peak       {peak(a,'ICP',0,42)[0]:6.0f}  ->  "
          f"{peak(dx,'ICP',0,42)[0]:6.0f} mmH2O")
    print(f"  CSF WBC d14    {at(a,14,'WBC'):6.1f}  ->  {at(dx,14,'WBC'):6.1f} /uL")
    print(f"  10-wk death    {mortality(a,70)*100:5.1f}% ->  "
          f"{mortality(dx,70)*100:5.1f}%   (trial 41% -> 47%)")
    print(f"  injury peak    {peak(a,'NEUR',0,70)[0]:6.3f}  ->  "
          f"{peak(dx,'NEUR',0,70)[0]:6.3f}")
    print(f"  disability d70 {at(a,70,'DIS'):6.3f}  ->  {at(dx,70,'DIS'):6.3f}"
          f"   (x{at(dx,70,'DIS')/at(a,70,'DIS'):.2f}; trial good outcome "
          f"25% -> 13%, x0.52)")

    print("\n-- Thesis E: ART timing (COAT backbone) -----------------------")
    print(f"  {'ART day':>8s} {'26wk death':>11s} {'IRIS peak':>10s}"
          f" {'GXM @ART':>9s} {'ICP pk':>7s}")
    for adt in (7, 10, 14, 21, 28, 35, 42, 56):
        r = run("coat", days=182.0, p=p, y0=y0, art_day=float(adt),
                record=0.5)
        print(f"  {adt:8d} {mortality(r,182)*100:10.1f}%"
              f" {peak(r,'IRISa',0,182)[0]:10.3f}"
              f" {at(r,adt,'GXM'):9.0f} {peak(r,'ICP',0,182)[0]:7.0f}")
    print("  (COAT: 45% at 1-2 wks vs 30% at 5 wks, 26-week mortality)")

    print("\n-- Thesis F: sertraline, and why ASTRO-CM had to fail ---------")
    sr = runs["ambition_sert"]
    cs = peak(sr, "CsrBrFree", 0, 14)[0]
    print(f"  free sertraline in brain, peak {cs*1000:.1f} ng/mL "
          f"= {cs/P['EC50SR']*100:.2f}% of EC50 ({P['EC50SR']:.0f} ug/mL)")
    print(f"  EFA {efa(runs['ambition']):.3f} -> {efa(sr):.3f}"
          f"   10-wk death {mortality(runs['ambition'],70)*100:.1f}% -> "
          f"{mortality(sr,70)*100:.1f}%")

    print("\n-- Adjunctive IFN-gamma ---------------------------------------")
    ig = runs["ambition_ifng"]
    print(f"  EFA {efa(runs['ambition']):.3f} -> {efa(ig):.3f}"
          f"   CSF IFN-g peak {peak(ig,'IFNG',0,21)[0]:.0f} pg/mL"
          f"   10-wk death {mortality(ig,70)*100:.1f}%")

    print("\n-- CSF sterility and relapse ---------------------------------")
    for name in ("oral_acta", "ambd_5fc_1wk", "ambition", "ambition_lp"):
        r = runs[name]
        sd = sterile_day(r)
        print(f"  {name:16s} culture-negative day "
              f"{('%.1f' % sd) if sd else '  -  '}"
              f"   GXM at that day "
              f"{(at(r,sd,'GXM') if sd else float('nan')):6.1f} ug/mL"
              f"   GXM d70 {at(r,70,'GXM'):6.1f}"
              f"   parenchymal log {log10(max(at(r,70,'Fp'),1e-9)):5.2f}")

    print("\n-- Virtual patients (each is its own model equilibrium) --------")
    print(f"  {'phenotype':15s} {'logCFU':>7s} {'GXM':>6s} {'OP':>5s} {'WBC':>6s}"
          f" {'CD4':>4s} | {'EFA':>6s} {'10wk':>6s} {'LP 10wk':>8s} {'AUC>250':>8s}")
    for ph, cfg in PHENOTYPES.items():
        yb = baseline(p, logCFU=cfg["logCFU"], CD4=cfg["CD4"])
        db = derived(yb, p)
        ra = run("ambition", days=70.0, p=p, y0=yb, record=0.5)
        rl = run("ambition_lp", days=70.0, p=p, y0=yb, record=0.5)
        print(f"  {ph:15s} {db['logCFU']:7.2f} {yb[IX['GXM']]:6.0f}"
              f" {db['ICP']:5.0f} {yb[IX['WBC']]:6.1f} {yb[IX['CD4']]:4.0f}"
              f" | {efa(ra):6.2f} {mortality(ra,70)*100:5.1f}%"
              f" {mortality(rl,70)*100:7.1f}% {at(ra,70,'AICP'):8.0f}")
    print("  (Rolfes 2014: >=1 therapeutic LP, adjusted RR of 11-day death 0.31)")
    yb = baseline(p, logCFU=PHENOTYPES["high_burden"]["logCFU"], CD4=25.0)
    ra = run("ambition", days=70.0, p=p, y0=yb, record=0.25)
    rl = run("ambition_lp", days=70.0, p=p, y0=yb, record=0.25)
    m11a, m11l = mortality(ra, 11.0), mortality(rl, 11.0)
    print(f"  high-burden phenotype, 11-day mortality {m11a*100:.1f}% -> "
          f"{m11l*100:.1f}%  (RR {m11l/m11a:.2f})")

    print("\n-- Paucicellular subgroup x ART timing (COAT subgroup) --------")
    for ph in ("median", "paucicellular"):
        yb = baseline(p, logCFU=PHENOTYPES[ph]["logCFU"],
                      CD4=PHENOTYPES[ph]["CD4"])
        m = {}
        for adt in (10.0, 35.0):
            rr = run("coat", days=182.0, p=p, y0=yb, art_day=adt, record=0.5)
            m[adt] = (mortality(rr, 182.0), at(rr, 182.0, "HAZ"),
                      peak(rr, "IRISa", 0, 182)[0])
        hr = m[10.0][1] / m[35.0][1]
        print(f"  {ph:15s} early ART {m[10.0][0]*100:5.1f}%  deferred "
              f"{m[35.0][0]*100:5.1f}%   hazard ratio {hr:.2f}"
              f"   IRIS peak {m[10.0][2]:.2f} vs {m[35.0][2]:.2f}")
    print("  (COAT: overall HR 1.73 [1.06-2.82]; CSF WBC <5/uL subgroup "
          "HR 3.87 [1.41-10.58])")

    print("\n-- ACTA partner-drug comparison (flucytosine vs fluconazole) ---")
    m5 = mortality(runs["ambd_5fc"], 70.0)
    mf = mortality(runs["ambd_flu"], 70.0)
    print(f"  AmB + 5FC {m5*100:.1f}%   AmB + FLU {mf*100:.1f}%   "
          f"hazard ratio {at(runs['ambd_5fc'],70,'HAZ')/at(runs['ambd_flu'],70,'HAZ'):.2f}")
    print(f"  (ACTA: 31.1% vs 45.0%, HR 0.62 [0.45-0.84]) -- and in the model "
          f"the ONLY\n   difference is the ergosterol node: EC50(AmB) "
          f"{at(runs['ambd_5fc'],7,'EC50Aeff'):.2f} vs "
          f"{at(runs['ambd_flu'],7,'EC50Aeff'):.2f} ug/mL at day 7")

    print("\n-- Mass / consistency checks ---------------------------------")
    r = runs["ambition"]
    # AmB mass balance for the single liposomal dose
    dose = 10.0 * P["WT"]
    y = r["y"][-1]
    elim = dose - (y[IX["Al"]] + y[IX["Al2"]] + y[IX["Abr"]] + y[IX["Akid"]])
    print(f"  L-AmB dose {dose:.0f} mg; remaining in body at d70 "
          f"{dose-elim:.3f} mg; accounted {100*(elim+(dose-elim))/dose:.6f}%")
    print(f"  LP arm drained {at(runs['ambition_lp'],70,'CLPV'):.0f} mL CSF "
          f"and {at(runs['ambition_lp'],70,'CLPG')/1000:.2f} mg GXM over 7 LPs")
    return runs


def hazard_shares(r, t):
    """Integrate each hazard term separately up to t."""
    p = r["p"]
    keys = ["base", "burden", "ICP", "CPP", "anaemia", "neutropenia",
            "hypoK", "GFR", "IRIS", "injury", "unsterile", "steroid"]
    acc = {k: 0.0 for k in keys}
    ts = r["t"]
    for j in range(len(ts) - 1):
        if ts[j] > t:
            break
        h = min(ts[j + 1], t) - ts[j]
        if h <= 0:
            continue
        y = r["y"][j]
        d = derived(y, p)
        acc["base"] += p["h0"] * h
        acc["burden"] += p["hB"] * max(0.0, d["logCFU"] - 3.0) * h
        acc["ICP"] += p["hICP"] * max(0.0, d["ICP"] - 250.0) / 250.0 * h
        acc["CPP"] += p["hCPP"] * max(0.0, 60.0 - d["CPP"]) / 20.0 * h
        acc["anaemia"] += p["hAN"] * max(0.0, p["Hb0"] - y[IX["Hb"]]) * h
        acc["neutropenia"] += p["hNEUT"] * max(0.0, p["ANCthr"] - y[IX["ANC"]]) * h
        acc["hypoK"] += p["hK"] * max(0.0, p["Kthr"] - y[IX["Kser"]]) * h
        acc["GFR"] += p["hGFR"] * max(0.0, p["GFRthr"] - y[IX["GFR"]]) / 30.0 * h
        acc["IRIS"] += p["hIRIS"] * y[IX["IRISa"]] * h
        acc["unsterile"] += (p["hSTER"] if d["CFU"] >= 1.0 else 0.0) * h
        acc["steroid"] += p["hDEX"] * d["dex"] * h
        acc["injury"] += p["hAMS"] * y[IX["NEUR"]] * 10.0 * h
    return acc


def trapz(r, key, t0, t1):
    ts = r["t"]
    tot = 0.0
    for j in range(len(ts) - 1):
        if ts[j] < t0 - 1e-9 or ts[j] > t1:
            continue
        a = (derived(r["y"][j], r["p"])[key] if key not in IX
             else r["y"][j][IX[key]])
        b = (derived(r["y"][j + 1], r["p"])[key] if key not in IX
             else r["y"][j + 1][IX[key]])
        tot += 0.5 * (a + b) * (ts[j + 1] - ts[j])
    return tot


def t_above(r, key, thresh):
    ts = r["t"]
    tot = 0.0
    for j in range(len(ts) - 1):
        a = (derived(r["y"][j], r["p"])[key] if key not in IX
             else r["y"][j][IX[key]])
        if a >= thresh:
            tot += ts[j + 1] - ts[j]
    return tot


if __name__ == "__main__":
    import sys
    report(brief="--brief" in sys.argv)
