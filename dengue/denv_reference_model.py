#!/usr/bin/env python3
# =============================================================================
#  DENGUE / SEVERE DENGUE (DSS)  --  INDEPENDENT REFERENCE IMPLEMENTATION
# =============================================================================
#  Purpose
#  -------
#  This file re-implements, from the written equations rather than by
#  translating code, the same 45-state ODE system that denv_mrgsolve_model.R
#  encodes.  It exists so that every number quoted in README.md and in the
#  root README table has been produced by a solver rather than asserted, and
#  so that the mrgsolve file has an independent check.
#
#  The organising idea
#  -------------------
#  Dengue is written here as ONE antibody response read through TWO channels
#  of opposite sign.
#
#    channel 1 (protective)  neutralising IgG clears virions -> viraemia falls
#                            -> interferon falls -> the fever breaks.
#    channel 2 (destructive) the SAME IgG binds NS1 and virion -> immune
#                            complexes -> complement -> mast-cell chymase and
#                            cross-reactive T-cell TNF -> the endothelial
#                            glycocalyx is stripped -> the reflection
#                            coefficient sigma falls -> plasma leaks.
#
#  Because the two channels share a driver, defervescence and the onset of
#  plasma leakage are not two events that coincide; they are one event
#  observed twice.  The 24-48 h "critical phase" is then not a clinical rule
#  imposed on the model -- it is the width of the antibody rise.
#
#  The second organising idea is that the leak is self-limiting only through
#  shock.  Filtration is  J_v = Kf * [ (Pc - Pi) - sigma*(PIp - PIi) ]  and
#  the only term that falls fast enough to switch it off is Pc, which falls
#  because the plasma volume feeding it has gone.  Resuscitation therefore
#  restores perfusion by restarting the leak, and the optimum fluid dose is a
#  computed interior minimum rather than a preference.
#
#  Run:  python3 denv_reference_model.py
#  Writes: denv_reference_output.txt, denv_scenario_results.json
# =============================================================================

import json
import math

import numpy as np
from scipy.integrate import solve_ivp

# -----------------------------------------------------------------------------
# 0.  STATE VECTOR  (45 ODEs)
# -----------------------------------------------------------------------------
NAMES = [
    "TGT",    #  0 susceptible FcgR+ target cells (monocyte/macrophage/DC), cells/mL
    "ECL",    #  1 eclipse-phase infected cells, cells/mL
    "INF",    #  2 productively infected cells, cells/mL
    "V",      #  3 plasma viraemia, RNA copies/mL
    "NS1",    #  4 circulating NS1 antigen, ng/mL
    "ABH",    #  5 pre-existing heterotypic IgG, reciprocal PRNT50 titre
    "ABN",    #  6 de-novo (homotypic, neutralising) IgG, reciprocal titre
    "PBL",    #  7 plasmablast / antibody-secreting cell pool, relative
    "CTL",    #  8 activated CD8 T cells, relative
    "IFN",    #  9 type-I interferon, pg/mL
    "TNF",    # 10 TNF-alpha, pg/mL
    "IL10",   # 11 IL-10, pg/mL
    "VEGF",   # 12 free VEGF-A, pg/mL
    "CHYM",   # 13 mast-cell chymase, ng/mL
    "GLX",    # 14 endothelial glycocalyx integrity, 0-1
    "LPX",    # 15 hydraulic conductance multiplier Kf/Kf0
    "VP",     # 16 plasma volume, mL
    "VI",     # 17 interstitial volume, mL
    "PRP",    # 18 plasma protein mass, g
    "PRI",    # 19 interstitial protein mass, g
    "VSER",   # 20 serosal (pleural + peritoneal) effusion volume, mL
    "RBCV",   # 21 red-cell volume, mL
    "PLT",    # 22 platelets, 1e9/L
    "MKC",    # 23 megakaryocyte / marrow platelet output capacity, relative
    "APLT",   # 24 anti-platelet (anti-NS1 cross-reactive) antibody, relative
    "FIB",    # 25 fibrinogen, mg/dL
    "WBC",    # 26 leukocytes, 1e9/L
    "HEP",    # 27 viable hepatocyte fraction, 0-1
    "AST",    # 28 AST, U/L
    "ALT",    # 29 ALT, U/L
    "PGE",    # 30 hypothalamic PGE2, relative
    "TEMP",   # 31 core temperature, degC
    "SVRR",   # 32 systemic vascular resistance multiplier (baroreflex state)
    "HR",     # 33 heart rate, bpm
    "LAC",    # 34 arterial lactate, mmol/L
    "COLL",   # 35 synthetic colloid mass in plasma, g
    "AVD",    # 36 antiviral gut depot, mg
    "AVC",    # 37 antiviral plasma concentration, mg/L
    "APAP",   # 38 paracetamol concentration, mg/L
    "STER",   # 39 methylprednisolone concentration, mg/L
    "CIN",    # 40 cumulative fluid in (IV + oral), mL
    "COUT",   # 41 cumulative urine out, mL
    "NAUC",   # 42 cumulative NS1 exposure, ng/mL*h
    "ICAUC",  # 43 cumulative immune-complex exposure, relative*h
    "ASC",    # 44 antibody-secreting cells (plasmablast -> ASC maturation stage)
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)

# -----------------------------------------------------------------------------
# 1.  PARAMETERS   (per hour unless stated)
# -----------------------------------------------------------------------------
P = dict(
    # ---- host ------------------------------------------------------------
    WT       = 70.0,
    VP0      = 2800.0,    # mL  plasma volume, 40 mL/kg
    VI0      = 10500.0,   # mL  interstitial volume, 150 mL/kg
    RBC0     = 2000.0,    # mL  red-cell volume -> baseline Hct 41.7 %
    TP_P0    = 7.30,      # g/dL plasma total protein
    TP_I0    = 2.00,      # g/dL interstitial total protein (Guyton)

    # ---- target-cell-limited virology ------------------------------------
    TGT0     = 2.0e4,     # cells/mL productively infectable FcgR+ pool
    dT       = 0.010/24,
    KVSAT    = 1.5e7,     # copies/mL: FcgR-mediated uptake saturates
    kECL     = 0.25,      # /h  eclipse exit (4 h)
    dINF0    = 1.0/24,    # /h  infected-cell death, lifespan ~1 d
    kCTLkill = 2.4/24,    # /h  per unit CTL
    BETA0    = 1.19e-9,   # mL/(copy*h)  FcgR-independent entry
    PPROD0   = 1050.0,    # copies/(cell*h)  burst rate, IFN-free
    PINT     = 3.00,      # intrinsic ADE: extra burst size at full opsonisation
    IC50IFN  = 400.0,     # pg/mL  IFN concentration halving viral production
    cV       = 3.5/24,    # /h  virion clearance (t1/2 4.8 h)
    kNeutCl  = 2.0e-4,    # /h per titre unit -- opsonic clearance of virions
    pNS1     = 0.0155,    # (ng/mL)/(cell/mL)/h  NS1 secretion
    cNS1     = 0.098,     # /h  (t1/2 7 h)
    kNS1ab   = 2.2e-4,    # /h per titre unit -- immune-complex clearance

    # ---- antibody-dependent enhancement -----------------------------------
    FCROSS   = 0.12,      # heterotypic IgG is 8.3x less neutralising
    AVIDN    = 3.0,       # homotypic de-novo IgG engages quaternary epitopes and
                          # neutralises at ~3x lower titre (de Alwis 2012, EDE mAbs)
    NT50     = 30.0,      # reciprocal titre for 50 % neutralisation
    HNEUT    = 2.5,       # multi-hit neutralisation
    KOPS     = 4.0,       # reciprocal titre for 50 % FcgR ligation
    PHI      = 14.0,      # FcgR amplification of monocyte entry
    KINT     = 0.85,      # intrinsic ADE strength (IFN suppression)

    # ---- adaptive response -------------------------------------------------
    kPBL     = 0.055, KVPBL = 3.0e6, dPBL = 0.030,
    kPBLM    = 0.185,     # memory B cells answer sooner than naive B cells
    kMAT     = 0.030,     # /h  B cell -> antibody-secreting cell maturation (~33 h);
                          # this lag is what puts the antibody rise, and therefore
                          # BOTH defervescence and the leak, at illness day 4-6
    dASC     = 0.021,
    kAB      = 3.30, dAB = 1.375e-3,        # IgG t1/2 21 d
    dABH     = 6.716e-4,                    # heterotypic waning t1/2 43 d
    kCTL     = 0.058, KICTL = 40.0, dCTL = 0.025,
    MEMT     = 2.4,       # memory (cross-reactive) T-cell amplification
    XAVID    = 0.45,      # cross-reactive T cells kill at 45 % efficiency ...
    XTNF     = 3.1,       # ... and make 3.1x the TNF (original antigenic sin)

    # ---- innate mediators ---------------------------------------------------
    kIFN     = 0.0085, dIFN = 0.10,
    kTNFm    = 34.0, kTNFt = 3.2, dTNF = 0.29,
    kIL10    = 620.0, dIL10 = 0.21,
    kVEGF    = 78.0, dVEGF = 0.14,
    kCHYM    = 13.5, dCHYM = 0.11,

    # ---- endothelium ---------------------------------------------------------
    kGSYN    = 0.0065,    # /h regeneration (t1/2 107 h)
    kGDEG    = 0.082,     # /h maximal stripping.  The glycocalyx is a
                          # STRUCTURE: it integrates the mediator load over
                          # days, which is why the leak lands at illness day
                          # 4-6 and not at the cytokine peak on day 2.
    # Each driver enters as a Hill-2 term: the glycocalyx tolerates a
    # mediator load and then gives way, which is what separates dengue
    # fever from dengue haemorrhagic fever on a continuous driver axis.
    wNS1 = 0.20, KNS1 = 1100.0,  # direct NS1 sialidase/heparanase (Puerta-Guardo 2016)
    wTNF = 0.36, KTNF = 110.0,
    wVEG = 0.12, KVEG = 520.0,
    wCHY = 0.32, KCHY = 11.0,    # mast-cell chymase (Tissera 2017)
    SIG0     = 0.97,      # baseline reflection coefficient
    SIGMIN   = 0.55,      # floor at complete glycocalyx loss
    NSIG     = 1.30,
    ALPHLP   = 1.20,      # Kf can rise 2.2-fold as junctions open
    tauLP    = 3.0,

    # ---- Starling ------------------------------------------------------------
    KF0      = 125.89,    # mL/(h*mmHg): OVERWRITTEN by calibrate() below
    PI0      = -3.0,      # mmHg
    PIMAX    = 2.6,       # mmHg interstitial pressure ceiling
    VSC      = 1600.0,    # mL  interstitial compliance scale
    JL0      = 140.0,     # mL/h baseline lymph flow
    kLYMPH   = 154.0,     # mL/(h*mmHg)
    JLMAX    = 1100.0,    # mL/h lymph pump ceiling
    PSPROT   = 47.0,      # mL/h diffusive protein permeability-surface product
    SIGCOL_F = 1.02,      # HES sieved slightly better than native protein
    PRSYN0   = 1.04,      # g/h total plasma-protein synthesis (25 g/day)
    PRSUP    = 0.55, KPRSUP = 60.0,   # negative acute-phase suppression
    kPRDEG   = 0.00509,   # /h
    # ---- serosal (pleural + peritoneal) compartment ----------------------
    # This is a SECOND Starling bed and it is what makes dengue look like
    # dengue.  Its drainage is not systemic lymph but parietal pleural and
    # diaphragmatic lymphatics, whose maximum is ~0.6 mL/kg/h -- two orders
    # of magnitude below systemic lymph flow.  Fluid therefore accumulates
    # here and is limited only by the pressure it generates.
    KFS      = 19.48,     # mL/(h*mmHg): OVERWRITTEN by calibrate() below
    PSER0    = -5.0,      # mmHg  resting pleural/peritoneal pressure
    CSER     = 150.0,     # mL/mmHg  combined serosal compliance
    TP_S0    = 1.50,      # g/dL  resting serosal protein
    TP_SMAX  = 3.20,      # g/dL  protein of an exudative dengue effusion
    DSER0    = 20.0,      # mL/h  baseline serosal turnover
    VSERB    = 20.0,      # mL    resting pleural + peritoneal fluid volume
    DSERMAX  = 46.0,      # mL/h  ceiling of parietal lymphatic absorption
                          #       (~0.65 mL/kg/h, Miserocchi 1997)

    # ---- haemodynamics --------------------------------------------------------
    SVMAX    = 140.0, VPUN = 1800.0, KSV = 1000.0,
    kINO     = 0.35, kMYOC = 0.28,
    HR0      = 72.0, kHRB = 46.0, kHRT = 7.5, tauHR = 0.10,
    SVR0     = 17.46, BARGAIN = 5.2, SVRMAX = 2.55, tauSVR = 0.16,
    MAPSET   = 88.0,
    kVASOPL  = 0.55, KVASOPL = 140.0,
    CART     = 1.82,      # mL/mmHg -> pulse pressure = SV / CART
    PV0      = 8.0, CVEN = 200.0, FPLAS = 0.583, PVMIN = 1.0,
    RRATIO0  = 7.6, kRRAT = 2.5,

    # ---- renal / intake --------------------------------------------------------
    UO0      = 62.0, KMAPU = 62.0, NMAPU = 3.0,
    kUOVOL   = 4.0,       # volume (pressure) diuresis exponent
    INSENS   = 43.0, ORAL0 = 105.0,
    ANOREX   = 0.12,      # fractional fall in oral intake at peak fever

    # ---- haematology ------------------------------------------------------------
    kPLTP    = 1.05, dPLT = 0.00417,
    kPLTI    = 0.265, kPLTC = 0.030,
    kMKR     = 0.011, kMKS = 0.055, KVMK = 2.0e6,
    kAPLT    = 0.0125, dAPLT = 0.011,
    kFIBC    = 0.020, kFIBS = 0.0090, FIBAPR = 1.45,
    kWBCS    = 0.075, kWBCR = 0.030,
    dPLTTX   = 0.0,       # (transfused platelets share the PLT compartment)

    # ---- liver ---------------------------------------------------------------
    kHEPV    = 0.0042, KVHEP = 1.0e8,
    kHEPS    = 0.020, kHEPAP = 1.1e-4, APAPTH = 22.0,
    kHEPR    = 0.0055,
    kAST     = 7.4e3, kALT = 3.1e3,
    dAST     = 0.041, dALT = 0.0154,

    # ---- fever ---------------------------------------------------------------
    kPGE     = 0.155, KPYR = 62.0, dPGE = 0.115,
    TSETMAX  = 3.10, KPGE = 0.52, tauT = 1.6,
    APAPIC50 = 9.0,

    # ---- drugs ----------------------------------------------------------------
    kaAV = 1.10, kelAV = 0.077, VdAV = 92.0,
    EMAXAV = 0.965, EC50AV = 0.42,
    kelAPAP = 0.315, VdAPAP = 49.0,
    kelSTER = 0.231, VdSTER = 84.0,
    STEREM = 0.72, STEREC = 0.55, STERVIR = 0.70, STERAB = 0.45,
    COLLEL = 0.010, COLLCOP = 5.05,
)


# -----------------------------------------------------------------------------
# 2.  ALGEBRAIC HELPERS
# -----------------------------------------------------------------------------
def calibrate(p):
    """
    Fix the two filtration coefficients by REQUIRING a stationary baseline:
    systemic filtration must equal baseline lymph flow, and serosal
    filtration must equal baseline serosal turnover.  Nothing here is a free
    parameter -- Kf and Kfs are whatever the textbook pressures make them.
    """
    Pa, Pv = p["MAPSET"], p["PV0"]
    r = p["RRATIO0"]
    Pc0 = (Pa + r * Pv) / (1.0 + r)
    PIp0 = cop_plasma(p["TP_P0"])
    PIi0 = cop_interstitium(p["TP_I0"])
    PIs0 = cop_interstitium(p["TP_S0"])
    nfp_sys = (Pc0 - p["PI0"]) - p["SIG0"] * (PIp0 - PIi0)
    nfp_ser = (Pc0 - p["PSER0"]) - p["SIG0"] * (PIp0 - PIs0)
    p["KF0"] = p["JL0"] / nfp_sys
    p["KFS"] = p["DSER0"] / nfp_ser
    p["_Pc0"], p["_NFPsys0"], p["_NFPser0"] = Pc0, nfp_sys, nfp_ser
    return p


def cop_plasma(c):
    """Landis-Pappenheimer colloid osmotic pressure, c = TOTAL protein g/dL."""
    return 2.1 * c + 0.16 * c * c + 0.009 * c * c * c


def cop_interstitium(c):
    """Interstitial COP; calibrated to Guyton's 8 mmHg at 2.0 g/dL."""
    return 3.6 * c + 0.20 * c * c


def ade_factors(abh, abn, p):
    """
    Antibody-dependent enhancement.

    Two occupancies are read off the same antibody pool:

      neutralisation  needs a high, multi-hit stoichiometric occupancy, driven
                      by the NEUTRALISING-EQUIVALENT titre Aeff = ABN+FCROSS*ABH
      opsonisation    needs only one or two IgG per virion to ligate FcgR,
                      driven by the TOTAL titre  Aopz = ABN + ABH.

    With O the opsonised-but-not-neutralised fraction,

        E = (1 - N - O) + PHI*O = (1 - N) * [1 + (PHI-1)*Aopz/(Aopz+KOPS)]

    exactly 1 in a naive host, rising to a maximum at an intermediate titre and
    falling back through 1 once neutralisation takes over.  The bell shape is
    the product of a rising saturating term and a falling Hill term; it is not
    imposed.
    """
    aeff = p["AVIDN"] * abn + p["FCROSS"] * abh
    aopz = abn + abh
    n = (aeff ** p["HNEUT"] / (aeff ** p["HNEUT"] + p["NT50"] ** p["HNEUT"])
         if aeff > 0 else 0.0)
    x = aopz / (aopz + p["KOPS"])
    o = (1.0 - n) * x
    entry = (1.0 - n) * (1.0 + (p["PHI"] - 1.0) * x)
    return entry, o, n


calibrate(P)


def initial_state(p, abh0=0.0, v0=4.0e5):
    """
    t = 0 is FEVER ONSET (illness day 0), not the mosquito bite.  The 4-7 day
    incubation is deliberately outside the model: within-host dengue data begin
    at presentation, and starting every scenario from the same viraemia makes
    the antibody state -- not a different starting point -- the only thing that
    differs between arms.  ECL and INF are seeded at the values consistent with
    exponential growth at rate lambda, so the trajectory is continuous.
    """
    lam = 0.080          # /h, pre-symptomatic growth rate
    y = np.zeros(NST)
    y[IX["TGT"]]  = p["TGT0"]
    y[IX["V"]]    = v0
    y[IX["INF"]]  = (lam + p["cV"]) * v0 / p["PPROD0"]
    y[IX["ECL"]]  = (lam + p["dINF0"]) * y[IX["INF"]] / p["kECL"]
    y[IX["TGT"]]  = p["TGT0"] - y[IX["INF"]] - y[IX["ECL"]]
    febrile = v0 > 0.0
    y[IX["IFN"]]  = 28.0 if febrile else 0.0
    y[IX["PGE"]]  = 0.49 if febrile else 0.0
    y[IX["ABH"]]  = abh0
    y[IX["GLX"]]  = 1.0
    y[IX["LPX"]]  = 1.0
    y[IX["VP"]]   = p["VP0"]
    y[IX["VI"]]   = p["VI0"]
    y[IX["PRP"]]  = p["TP_P0"] / 100.0 * p["VP0"]
    y[IX["PRI"]]  = p["TP_I0"] / 100.0 * p["VI0"]
    y[IX["RBCV"]] = p["RBC0"]
    y[IX["PLT"]]  = 250.0
    y[IX["MKC"]]  = 1.0
    y[IX["FIB"]]  = 300.0
    y[IX["WBC"]]  = 6.5
    y[IX["HEP"]]  = 1.0
    y[IX["AST"]]  = 24.0
    y[IX["ALT"]]  = 22.0
    y[IX["TEMP"]] = 38.4 if febrile else 37.0
    y[IX["SVRR"]] = 1.0
    y[IX["HR"]]   = p["HR0"]
    y[IX["LAC"]]  = 1.0
    y[IX["VSER"]] = p["VSERB"]
    return y


class Regimen:
    """Everything a scenario is allowed to do."""

    def __init__(self, **kw):
        self.crystalloid = kw.get("crystalloid", [])   # (t0, t1, mL/h)
        self.colloid     = kw.get("colloid", [])       # (t0, t1, mL/h of 6 % HES)
        self.albumin     = kw.get("albumin", [])       # (t0, t1, g/h)
        self.oral_extra  = kw.get("oral_extra", [])    # (t0, t1, mL/h)
        self.plt_tx      = kw.get("plt_tx", [])        # (t, 1e9/L)
        self.antiviral   = kw.get("antiviral", None)   # (t0, mg q12h, n)
        self.apap        = kw.get("apap", None)        # (t0, mg q6h, n)
        self.steroid     = kw.get("steroid", None)     # (t0, mg q24h, n)
        self.ns1mab      = kw.get("ns1mab", 0.0)       # extra /h NS1 clearance

    @staticmethod
    def _rate(sched, t):
        return sum(v for t0, t1, v in sched if t0 <= t < t1)

    def crys_rate(self, t): return self._rate(self.crystalloid, t)
    def coll_rate(self, t): return self._rate(self.colloid, t)
    def alb_rate(self, t):  return self._rate(self.albumin, t)
    def oral_rate(self, t): return self._rate(self.oral_extra, t)

    def bolus_times(self, p):
        ev = [(t, IX["PLT"], amt) for t, amt in self.plt_tx]
        if self.antiviral:
            t0, dose, n = self.antiviral
            ev += [(t0 + 12.0 * k, IX["AVD"], dose) for k in range(n)]
        if self.apap:
            t0, dose, n = self.apap
            ev += [(t0 + 6.0 * k, IX["APAP"], dose / p["VdAPAP"]) for k in range(n)]
        if self.steroid:
            t0, dose, n = self.steroid
            ev += [(t0 + 24.0 * k, IX["STER"], dose / p["VdSTER"]) for k in range(n)]
        return sorted(ev, key=lambda e: e[0])


# -----------------------------------------------------------------------------
# 3.  RIGHT-HAND SIDE
# -----------------------------------------------------------------------------
def rhs(t, y, p, reg):
    y = np.maximum(y, 0.0)
    (TGT, ECL, INF, V, NS1, ABH, ABN, PBL, CTL, IFN, TNF, IL10, VEGF, CHYM,
     GLX, LPX, VP, VI, PRP, PRI, VSER, RBCV, PLT, MKC, APLT, FIB, WBC,
     HEP, AST, ALT, PGE, TEMP, SVRR, HR, LAC, COLL, AVD, AVC, APAP, STER,
     CIN, COUT, NAUC, ICAUC, ASC) = y
    VP = max(VP, 300.0)
    VI = max(VI, 3000.0)

    d = np.zeros(NST)

    # ---------- 3.1  ADE ----------------------------------------------------
    entry, opsO, neutN = ade_factors(ABH, ABN, p)
    beta = p["BETA0"] * entry
    memory = 1.0 if ABH > 1.0 else 0.0

    # ---------- 3.2  virology ------------------------------------------------
    av_block = p["EMAXAV"] * AVC / (AVC + p["EC50AV"])
    p_eff = (p["PPROD0"] * (1.0 + p["PINT"] * opsO) / (1.0 + IFN / p["IC50IFN"])
             * (1.0 - av_block))
    infect = beta * V * TGT / (1.0 + V / p["KVSAT"])
    dINF = p["dINF0"] + p["kCTLkill"] * CTL * (p["XAVID"] if memory else 1.0)

    d[IX["TGT"]] = p["dT"] * (p["TGT0"] - TGT) - infect
    d[IX["ECL"]] = infect - p["kECL"] * ECL
    d[IX["INF"]] = p["kECL"] * ECL - dINF * INF
    d[IX["V"]]   = (p_eff * INF - p["cV"] * V
                    - p["kNeutCl"] * (ABN + 0.35 * ABH) * V)
    d[IX["NS1"]] = (p["pNS1"] * (1.0 - av_block) * INF
                    - (p["cNS1"] + reg.ns1mab) * NS1
                    - p["kNS1ab"] * ABN * NS1)

    # ---------- 3.3  adaptive immunity ---------------------------------------
    antigen = V / (V + p["KVPBL"])
    d[IX["PBL"]] = ((p["kPBL"] + p["kPBLM"] * memory) * antigen
                    - (p["dPBL"] + p["kMAT"]) * PBL)
    d[IX["ASC"]] = (p["kMAT"] * PBL * (1.0 - p["STERAB"] * STER / (STER + p["STEREC"]))
                    - p["dASC"] * ASC)
    d[IX["ABN"]] = p["kAB"] * ASC - p["dAB"] * ABN
    d[IX["ABH"]] = -p["dABH"] * ABH
    ster_imm = 1.0 - p["STERVIR"] * STER / (STER + p["STEREC"])
    d[IX["CTL"]] = (p["kCTL"] * INF / (INF + p["KICTL"]) * (1.0 + p["MEMT"] * memory)
                    * ster_imm - p["dCTL"] * CTL)

    # ---------- 3.4  immune complexes and cytokines --------------------------
    # ABN rises while NS1 and V fall, so the product peaks BETWEEN them.
    fAb = ABN / (ABN + 400.0)
    IC = fAb * (NS1 / (NS1 + 240.0) + 0.55 * V / (V + 4.0e6))
    ster_tnf = 1.0 - p["STEREM"] * STER / (STER + p["STEREC"])

    d[IX["IFN"]]  = (p["kIFN"] * INF / (1.0 + p["KINT"] * opsO) - p["dIFN"] * IFN)
    d[IX["TNF"]]  = ((p["kTNFm"] * IC
                      + p["kTNFt"] * CTL * (p["XTNF"] if memory else 1.0)) * ster_tnf
                     - p["dTNF"] * TNF)
    d[IX["IL10"]] = p["kIL10"] * (IC + 0.02 * CTL) - p["dIL10"] * IL10
    d[IX["VEGF"]] = p["kVEGF"] * (IC + 0.006 * TNF) - p["dVEGF"] * VEGF
    d[IX["CHYM"]] = p["kCHYM"] * IC - p["dCHYM"] * CHYM

    # ---------- 3.5  endothelium ---------------------------------------------
    hill2 = lambda x, k: (x * x) / (x * x + k * k)
    damage = (p["wNS1"] * hill2(NS1, p["KNS1"])
              + p["wTNF"] * hill2(TNF, p["KTNF"])
              + p["wVEG"] * hill2(VEGF, p["KVEG"])
              + p["wCHY"] * hill2(CHYM, p["KCHY"]))
    d[IX["GLX"]] = p["kGSYN"] * (1.0 - GLX) - p["kGDEG"] * GLX * damage
    d[IX["LPX"]] = (1.0 + p["ALPHLP"] * (1.0 - GLX) - LPX) / p["tauLP"]

    sigma = p["SIGMIN"] + (p["SIG0"] - p["SIGMIN"]) * GLX ** p["NSIG"]
    sigma_col = min(0.995, sigma * p["SIGCOL_F"])

    # ---------- 3.6  haemodynamics --------------------------------------------
    stressed = max(VP - p["VPUN"], 0.0)
    myoc = 1.0 - p["kMYOC"] * TNF / (TNF + 150.0)
    SV = (p["SVMAX"] * stressed / (p["KSV"] + stressed)
          * (1.0 + p["kINO"] * (SVRR - 1.0)) * myoc)
    CO = SV * HR / 1000.0
    MAP = max(CO * p["SVR0"] * SVRR, 5.0)
    PP = SV / p["CART"]

    ceiling = max(1.0, p["SVRMAX"] * (1.0 - p["kVASOPL"] * TNF / (TNF + p["KVASOPL"])))
    svr_target = min(max(1.0 + p["BARGAIN"] * (p["MAPSET"] - MAP) / p["MAPSET"], 0.85),
                     ceiling)
    d[IX["SVRR"]] = (svr_target - SVRR) / p["tauSVR"]
    hr_target = (p["HR0"] + p["kHRB"] * max(SVRR - 1.0, 0.0)
                 + p["kHRT"] * max(TEMP - 37.0, 0.0))
    d[IX["HR"]] = (hr_target - HR) / p["tauHR"]

    PV_cvp = max(p["PVMIN"], p["PV0"] + (VP - p["VP0"]) / p["FPLAS"] / p["CVEN"])
    rratio = p["RRATIO0"] * (1.0 + p["kRRAT"] * max(SVRR - 1.0, 0.0) / p["SVRMAX"])
    Pc = (MAP + rratio * PV_cvp) / (1.0 + rratio)

    # ---------- 3.7  Starling exchange ------------------------------------------
    Cp   = PRP / VP * 100.0
    Ci   = PRI / VI * 100.0
    Ccol = COLL / VP * 100.0
    PIp  = cop_plasma(Cp) + p["COLLCOP"] * Ccol
    PIi  = cop_interstitium(Ci)
    Pi   = p["PIMAX"] - (p["PIMAX"] - p["PI0"]) * math.exp(-max(VI - p["VI0"], 0.0) / p["VSC"])

    Jv = p["KF0"] * LPX * ((Pc - Pi) - sigma * (PIp - PIi))
    Jl = max(0.0, min(p["JL0"] + p["kLYMPH"] * (Pi - p["PI0"]), p["JLMAX"]))

    # --- the serosal bed: a second Starling exchange with a tiny drain -----
    Pser = p["PSER0"] + max(VSER - p["VSERB"], 0.0) / p["CSER"]
    Cser = p["TP_S0"] + (p["TP_SMAX"] - p["TP_S0"]) * (p["SIG0"] - sigma) \
        / max(p["SIG0"] - p["SIGMIN"], 1e-9)
    PIs = cop_interstitium(Cser)
    Jser = p["KFS"] * LPX * ((Pc - Pser) - sigma * (PIp - PIs))
    Dser = min(p["DSERMAX"], p["DSER0"] * math.sqrt(max(VSER, 0.0) / p["VSERB"]))

    crys = reg.crys_rate(t)
    coll = reg.coll_rate(t)
    albg = reg.alb_rate(t)
    anorexia = 1.0 - p["ANOREX"] * min(max(TEMP - 37.0, 0.0) / 2.5, 1.0)
    oral = (p["ORAL0"] * anorexia + reg.oral_rate(t)) * (0.6 if MAP < 60.0 else 1.0)
    insens = p["INSENS"] * (1.0 + 0.13 * max(TEMP - 37.0, 0.0))
    f_map = ((MAP ** p["NMAPU"] / (MAP ** p["NMAPU"] + p["KMAPU"] ** p["NMAPU"]))
             / (p["MAPSET"] ** p["NMAPU"]
                / (p["MAPSET"] ** p["NMAPU"] + p["KMAPU"] ** p["NMAPU"])))
    f_vol = min(3.0, max(0.15, (VP / p["VP0"]) ** p["kUOVOL"]))
    uo = p["UO0"] * f_map * f_vol

    d[IX["VP"]] = -Jv + Jl - Jser + Dser + crys + coll + oral - uo - insens
    d[IX["VI"]] = Jv - Jl
    d[IX["VSER"]] = Jser - Dser

    Js = (1.0 - sigma) * max(Jv, 0.0) * (PRP / VP) + p["PSPROT"] * (Cp - Ci) / 100.0
    back = Jl * (PRI / VI)
    Jsprot = (Jser - Dser) * Cser / 100.0          # protein following the effusion
    prsyn = p["PRSYN0"] * (1.0 - p["PRSUP"] * TNF / (TNF + p["KPRSUP"]))
    d[IX["PRP"]] = prsyn - Js + back + albg - Jsprot - p["kPRDEG"] * PRP
    d[IX["PRI"]] = Js - back
    d[IX["COLL"]] = (coll * 0.06 - (1.0 - sigma_col) * max(Jv, 0.0) * (COLL / VP)
                     - p["COLLEL"] * COLL)
    d[IX["RBCV"]] = 0.0
    d[IX["CIN"]]  = crys + coll + oral
    d[IX["COUT"]] = uo

    # ---------- 3.8  haematology -------------------------------------------------
    d[IX["MKC"]]  = p["kMKR"] * (1.0 - MKC) - p["kMKS"] * V / (V + p["KVMK"]) * MKC
    d[IX["APLT"]] = p["kAPLT"] * IC * 100.0 - p["dAPLT"] * APLT
    d[IX["PLT"]]  = (p["kPLTP"] * MKC - p["dPLT"] * PLT
                     - p["kPLTI"] * APLT * PLT / 250.0
                     - p["kPLTC"] * (1.0 - GLX) * PLT)
    fib_target = 300.0 * p["FIBAPR"] * (1.0 - 0.55 * (1.0 - HEP))
    d[IX["FIB"]] = p["kFIBS"] * (fib_target - FIB) - p["kFIBC"] * FIB * (1.0 - GLX)
    d[IX["WBC"]] = p["kWBCR"] * (6.5 - WBC) - p["kWBCS"] * WBC * V / (V + p["KVMK"])

    # ---------- 3.9  liver --------------------------------------------------------
    perf_def = max(0.0, 1.0 - CO / 4.6)
    apap_tox = p["kHEPAP"] * max(APAP - p["APAPTH"], 0.0)
    inj = (p["kHEPV"] * V / (V + p["KVHEP"]) + p["kHEPS"] * perf_def ** 2 + apap_tox) * HEP
    d[IX["HEP"]] = p["kHEPR"] * (1.0 - HEP) - inj
    d[IX["AST"]] = p["kAST"] * inj - p["dAST"] * (AST - 24.0)
    d[IX["ALT"]] = p["kALT"] * inj - p["dALT"] * (ALT - 22.0)

    # ---------- 3.10 fever ---------------------------------------------------------
    pyrogen = IFN + 0.06 * TNF
    cox = 1.0 / (1.0 + APAP / p["APAPIC50"])
    d[IX["PGE"]] = p["kPGE"] * pyrogen / (pyrogen + p["KPYR"]) * cox - p["dPGE"] * PGE
    tset = 37.0 + p["TSETMAX"] * PGE / (PGE + p["KPGE"])
    d[IX["TEMP"]] = (tset - TEMP) / p["tauT"]

    # ---------- 3.11 oxygen delivery / lactate --------------------------------------
    Hct = RBCV / (RBCV + VP) * 100.0
    DO2 = CO * (Hct / 3.0) * 1.34 * 0.97 * 10.0 / 1000.0        # L O2/min
    d[IX["LAC"]] = 0.24 + 2.6 * max(0.0, 1.0 - DO2 / 0.72) - 0.24 * LAC * max(HEP, 0.2)

    # ---------- 3.12 drug PK ----------------------------------------------------------
    d[IX["AVD"]]  = -p["kaAV"] * AVD
    d[IX["AVC"]]  = p["kaAV"] * AVD / p["VdAV"] - p["kelAV"] * AVC
    d[IX["APAP"]] = -p["kelAPAP"] * APAP
    d[IX["STER"]] = -p["kelSTER"] * STER

    # ---------- 3.13 exposure integrals -------------------------------------------------
    d[IX["NAUC"]]  = NS1
    d[IX["ICAUC"]] = IC
    return d


# -----------------------------------------------------------------------------
# 4.  DERIVED OUTPUTS
# -----------------------------------------------------------------------------
OBS_KEYS = ["Hct", "sigma", "Jv", "Jl", "netleak", "Pc", "Pi", "PIp", "PIi",
            "MAP", "SBP", "DBP", "PP", "CO", "SV", "SVRr", "TP", "Alb", "Eade",
            "Vser", "shock", "bleed", "fluidbal", "DO2", "Temp", "PLT", "V",
            "NS1", "ABN", "TNF", "IL10", "VEGF", "CHYM", "GLX", "Vp", "Vi",
            "AST", "ALT", "LAC", "HR", "WBC", "FIB", "IC", "IFN", "warnsigns",
            "Jser", "Dser", "Pser", "netser"]


def observe(t, Y, p):
    n = len(t)
    out = {k: np.zeros(n) for k in OBS_KEYS}
    for i in range(n):
        yy = np.maximum(Y[:, i], 0.0)
        VP = max(yy[IX["VP"]], 300.0); VI = max(yy[IX["VI"]], 3000.0)
        PRP, PRI, COLL = yy[IX["PRP"]], yy[IX["PRI"]], yy[IX["COLL"]]
        GLX, LPX, RBCV = yy[IX["GLX"]], yy[IX["LPX"]], yy[IX["RBCV"]]
        SVRR, HRv, TNF = yy[IX["SVRR"]], yy[IX["HR"]], yy[IX["TNF"]]

        sigma = p["SIGMIN"] + (p["SIG0"] - p["SIGMIN"]) * GLX ** p["NSIG"]
        stressed = max(VP - p["VPUN"], 0.0)
        myoc = 1.0 - p["kMYOC"] * TNF / (TNF + 150.0)
        SV = (p["SVMAX"] * stressed / (p["KSV"] + stressed)
              * (1.0 + p["kINO"] * (SVRR - 1.0)) * myoc)
        CO = SV * HRv / 1000.0
        MAP = max(CO * p["SVR0"] * SVRR, 5.0)
        PP = SV / p["CART"]
        PV_cvp = max(p["PVMIN"], p["PV0"] + (VP - p["VP0"]) / p["FPLAS"] / p["CVEN"])
        rratio = p["RRATIO0"] * (1.0 + p["kRRAT"] * max(SVRR - 1.0, 0.0) / p["SVRMAX"])
        Pc = (MAP + rratio * PV_cvp) / (1.0 + rratio)
        Cp, Ci, Ccol = PRP / VP * 100.0, PRI / VI * 100.0, COLL / VP * 100.0
        PIp = cop_plasma(Cp) + p["COLLCOP"] * Ccol
        PIi = cop_interstitium(Ci)
        Pi = p["PIMAX"] - (p["PIMAX"] - p["PI0"]) * math.exp(-max(VI - p["VI0"], 0.0) / p["VSC"])
        Jv = p["KF0"] * LPX * ((Pc - Pi) - sigma * (PIp - PIi))
        Jl = max(0.0, min(p["JL0"] + p["kLYMPH"] * (Pi - p["PI0"]), p["JLMAX"]))
        Pser = p["PSER0"] + max(yy[IX["VSER"]] - p["VSERB"], 0.0) / p["CSER"]
        Cser = p["TP_S0"] + (p["TP_SMAX"] - p["TP_S0"]) * (p["SIG0"] - sigma) \
            / max(p["SIG0"] - p["SIGMIN"], 1e-9)
        Jser = p["KFS"] * LPX * ((Pc - Pser) - sigma * (PIp - cop_interstitium(Cser)))
        Dser = min(p["DSERMAX"], p["DSER0"] * math.sqrt(max(yy[IX["VSER"]], 0.0) / p["VSERB"]))
        Hct = RBCV / (RBCV + VP) * 100.0
        DO2 = CO * (Hct / 3.0) * 1.34 * 0.97 * 10.0 / 1000.0
        entry, _, _ = ade_factors(yy[IX["ABH"]], yy[IX["ABN"]], p)
        fAb = yy[IX["ABN"]] / (yy[IX["ABN"]] + 400.0)
        IC = fAb * (yy[IX["NS1"]] / (yy[IX["NS1"]] + 240.0)
                    + 0.55 * yy[IX["V"]] / (yy[IX["V"]] + 4.0e6))

        # Haemostasis is a PRODUCT of four requirements, so the worst term
        # dominates and raising only one of them does almost nothing.
        f_plt = 1.0 / (1.0 + (max(yy[IX["PLT"]], 1e-6) / 18.0) ** 2.2)
        f_coa = 1.0 / (1.0 + (max(yy[IX["FIB"]], 1e-6) / 90.0) ** 2.5)
        f_ves = (1.0 - GLX) ** 2.0
        f_shk = 1.0 / (1.0 + (MAP / 52.0) ** 5.0)
        bleed = 1.0 - (1.0 - 0.95 * f_plt) * (1.0 - 0.75 * f_coa) \
            * (1.0 - 0.60 * f_ves) * (1.0 - 0.80 * f_shk)

        hct0 = p["RBC0"] / (p["RBC0"] + p["VP0"]) * 100.0
        warn = (int(Hct >= hct0 * 1.20) + int(yy[IX["PLT"]] < 100.0)
                + int(yy[IX["VSER"]] > 250.0) + int(yy[IX["AST"]] > 200.0)
                + int(PP <= 20.0) + int(yy[IX["LAC"]] > 2.5))

        vals = dict(Hct=Hct, sigma=sigma, Jv=Jv, Jl=Jl, netleak=Jv - Jl, Pc=Pc,
                    Pi=Pi, PIp=PIp, PIi=PIi, MAP=MAP, SBP=MAP + 2 * PP / 3,
                    DBP=MAP - PP / 3, PP=PP, CO=CO, SV=SV, SVRr=SVRR, TP=Cp,
                    Alb=0.575 * Cp, Eade=entry, Vser=max(yy[IX["VSER"]] - p["VSERB"], 0.0),
                    shock=1.0 if (PP <= 20.0 or MAP < 60.0) else 0.0, bleed=bleed,
                    fluidbal=yy[IX["CIN"]] - yy[IX["COUT"]], DO2=DO2,
                    Temp=yy[IX["TEMP"]], PLT=yy[IX["PLT"]], V=yy[IX["V"]],
                    NS1=yy[IX["NS1"]], ABN=yy[IX["ABN"]], TNF=TNF,
                    IL10=yy[IX["IL10"]], VEGF=yy[IX["VEGF"]], CHYM=yy[IX["CHYM"]],
                    GLX=GLX, Vp=VP, Vi=VI, AST=yy[IX["AST"]], ALT=yy[IX["ALT"]],
                    LAC=yy[IX["LAC"]], HR=HRv, WBC=yy[IX["WBC"]], FIB=yy[IX["FIB"]],
                    IC=IC, IFN=yy[IX["IFN"]], warnsigns=warn, Jser=Jser,
                    Dser=Dser, Pser=Pser, netser=Jser - Dser)
        for k in OBS_KEYS:
            out[k][i] = vals[k]
    return out


# -----------------------------------------------------------------------------
# 5.  DRIVER
# -----------------------------------------------------------------------------
def simulate(abh0=0.0, reg=None, tend=336.0, p=None, dt=0.5, v0=4.0e5):
    p = calibrate(dict(p)) if p is not None else P
    reg = reg or Regimen()
    y = initial_state(p, abh0, v0)
    events = reg.bolus_times(p)
    tgrid = np.arange(0.0, tend + 1e-9, dt)
    Y = np.zeros((NST, len(tgrid)))
    Y[:, 0] = y
    breaks = sorted(set([0.0] + [e[0] for e in events if 0 < e[0] < tend] + [tend]))
    idx = 0
    for a, b in zip(breaks[:-1], breaks[1:]):
        for te, si, amt in events:
            if abs(te - a) < 1e-9:
                y[si] += amt
        seg = tgrid[(tgrid > a) & (tgrid <= b)]
        if len(seg) == 0:
            continue
        sol = solve_ivp(rhs, (a, b), y, t_eval=seg, args=(p, reg),
                        method="LSODA", rtol=1e-6, atol=1e-8, max_step=1.0)
        if not sol.success:
            raise RuntimeError(sol.message)
        Y[:, idx + 1: idx + 1 + len(seg)] = sol.y
        idx += len(seg)
        y = sol.y[:, -1].copy()
    return tgrid, Y, observe(tgrid, Y, p)


# -----------------------------------------------------------------------------
# 6.  CLINICAL SUMMARY
# -----------------------------------------------------------------------------
def fever_onset(t, obs, thr=38.0):
    w = np.where(obs["Temp"] >= thr)[0]
    return float(t[w[0]]) if len(w) else float("nan")


def defervescence(t, obs, thr=37.8):
    w = np.where(obs["Temp"] >= 38.0)[0]
    if not len(w):
        return float("nan")
    ipk = int(w[0] + np.argmax(obs["Temp"][w[0]:]))
    below = np.where(obs["Temp"][ipk:] < thr)[0]
    return float(t[ipk + below[0]]) if len(below) else float("nan")


def presentation_time(t, obs, n_warn=3):
    """
    When would this patient actually be admitted?  Not at an arbitrary hour:
    at the moment the model itself first satisfies the WHO admission rule of
    three warning signs.  Rounded down to the nearest 6 h ward round.
    """
    w = np.where(obs["warnsigns"] >= n_warn)[0]
    return float(math.floor(t[w[0]] / 6.0) * 6.0) if len(w) else float("nan")


def severity(t, obs):
    """
    Composite severity, 0-1.  Five graded sub-scores, each with a clinical
    threshold below which it contributes nothing, combined with fixed weights.

    The respiratory term is SQUARED: a litre of effusion is tolerated and
    three litres are not, and without that curvature no amount of fluid ever
    looks harmful, which would make the dose-response monotone by
    construction rather than by physiology.
    """
    dt = t[1] - t[0]
    shock = min(1.0, float(np.sum(np.clip(25.0 - obs["PP"], 0, None)) * dt) / 1800.0)
    hypo  = min(1.0, float(np.sum(np.clip(65.0 - obs["MAP"], 0, None)) * dt) / 400.0)
    resp  = min(1.0, max(0.0, (float(np.max(obs["Vser"])) - 400.0) / 2800.0)) ** 2
    bleed = min(1.0, max(0.0, (float(np.max(obs["bleed"])) - 0.15) / 0.55))
    lac   = min(1.0, max(0.0, (float(np.max(obs["LAC"])) - 1.5) / 5.0))
    return float(min(1.0, 0.30 * shock + 0.15 * hypo + 0.25 * resp
                     + 0.15 * bleed + 0.15 * lac))


def severity_parts(t, obs):
    dt = t[1] - t[0]
    return dict(
        shock=min(1.0, float(np.sum(np.clip(25.0 - obs["PP"], 0, None)) * dt) / 1800.0),
        hypotension=min(1.0, float(np.sum(np.clip(65.0 - obs["MAP"], 0, None)) * dt) / 400.0),
        respiratory=min(1.0, max(0.0, (float(np.max(obs["Vser"])) - 400.0) / 2800.0)) ** 2,
        bleeding=min(1.0, max(0.0, (float(np.max(obs["bleed"])) - 0.15) / 0.55)),
        lactate=min(1.0, max(0.0, (float(np.max(obs["LAC"])) - 1.5) / 5.0)))


def summarise(t, Y, obs, p=None):
    p = p or P
    hct0 = p["RBC0"] / (p["RBC0"] + p["VP0"]) * 100.0
    dtg = t[1] - t[0]
    return dict(
        fever_onset_h=fever_onset(t, obs),
        defervescence_h=defervescence(t, obs),
        peak_viraemia_log10=float(np.log10(max(np.max(obs["V"]), 1e-9))),
        t_peak_viraemia_h=float(t[int(np.argmax(obs["V"]))]),
        peak_NS1=float(np.max(obs["NS1"])),
        NS1_AUC=float(Y[IX["NAUC"], -1]),
        IC_AUC=float(Y[IX["ICAUC"], -1]),
        peak_ABN=float(np.max(obs["ABN"])),
        min_sigma=float(np.min(obs["sigma"])),
        min_GLX=float(np.min(obs["GLX"])),
        max_Jv=float(np.max(obs["Jv"])),
        max_netleak=float(np.max(obs["netleak"])),
        max_netser=float(np.max(obs["netser"])),
        max_Hct=float(np.max(obs["Hct"])),
        hct_rise_pct=float((np.max(obs["Hct"]) - hct0) / hct0 * 100.0),
        min_Vp=float(np.min(obs["Vp"])),
        plasma_loss_pct=float((p["VP0"] - np.min(obs["Vp"])) / p["VP0"] * 100.0),
        max_Vser=float(np.max(obs["Vser"])),
        min_PP=float(np.min(obs["PP"])),
        min_MAP=float(np.min(obs["MAP"])),
        min_SBP=float(np.min(obs["SBP"])),
        max_HR=float(np.max(obs["HR"])),
        max_lactate=float(np.max(obs["LAC"])),
        shock_hours=float(np.sum(obs["shock"]) * dtg),
        min_PLT=float(np.min(obs["PLT"])),
        t_PLT_nadir_h=float(t[int(np.argmin(obs["PLT"]))]),
        min_WBC=float(np.min(obs["WBC"])),
        min_FIB=float(np.min(obs["FIB"])),
        min_TP=float(np.min(obs["TP"])),
        min_alb=float(np.min(obs["Alb"])),
        max_AST=float(np.max(obs["AST"])),
        max_ALT=float(np.max(obs["ALT"])),
        max_TNF=float(np.max(obs["TNF"])),
        max_IL10=float(np.max(obs["IL10"])),
        max_VEGF=float(np.max(obs["VEGF"])),
        max_CHYM=float(np.max(obs["CHYM"])),
        max_bleed=float(np.max(obs["bleed"])),
        max_warnsigns=float(np.max(obs["warnsigns"])),
        cum_fluid_in=float(Y[IX["CIN"], -1]),
        net_fluid_balance=float(obs["fluidbal"][-1]),
        respiratory_distress=bool(np.max(obs["Vser"]) >= 1500.0),
        severity_index=float(severity(t, obs)),
    )


# -----------------------------------------------------------------------------
# 7.  PROTOCOLS
# -----------------------------------------------------------------------------
def who_crystalloid(t0, wt=70.0, scale=1.0):
    """WHO 2009 compensated-shock ladder in mL/kg/h, anchored at presentation."""
    ladder = [(0.0, 1.0, 10.0), (1.0, 3.0, 6.0), (3.0, 5.0, 4.0),
              (5.0, 12.0, 2.5), (12.0, 36.0, 1.5), (36.0, 48.0, 1.0)]
    return [(t0 + a, t0 + b, r * wt * scale) for a, b, r in ladder]


def colloid_bolus(t0, wt=70.0, mlkg=15.0, over_h=1.0, n=1, spacing=6.0):
    return [(t0 + spacing * k, t0 + spacing * k + over_h, mlkg * wt / over_h)
            for k in range(n)]


# -----------------------------------------------------------------------------
# 8.  ANALYSES THAT DO NOT NEED THE ODE SOLVER
# -----------------------------------------------------------------------------
def ade_curve(p=None):
    p = p or P
    A = np.logspace(-1, 4.4, 1400)
    E = np.array([ade_factors(a, 0.0, p)[0] for a in A])
    ipk = int(np.argmax(E))
    above = np.where(E > 1.0)[0]
    i = above[-1]
    cross = 10 ** np.interp(1.0, [E[i + 1], E[i]], [np.log10(A[i + 1]), np.log10(A[i])])
    g = lambda a: float(ade_factors(a, 0.0, p)[0])
    return dict(titres=A, E=E, peak_titre=float(A[ipk]), peak_E=float(E[ipk]),
                cross_titre=float(cross),
                E_at_20=g(20), E_at_40=g(40), E_at_80=g(80), E_at_320=g(320),
                E_at_1280=g(1280), E_at_2560=g(2560), E_at_5120=g(5120))


def infant_age_prediction(cord_titre=1280.0, t_half_days=43.0, p=None):
    peak = ade_curve(p)["peak_titre"]
    n_half = math.log2(cord_titre / peak)
    return dict(cord_titre=cord_titre, t_half_days=t_half_days,
                peak_enh_titre=peak, half_lives=n_half,
                age_days=n_half * t_half_days,
                age_months=n_half * t_half_days / 30.44)


def hct_arithmetic(p=None):
    p = p or P
    h0 = p["RBC0"] / (p["RBC0"] + p["VP0"]) * 100.0
    h1 = h0 * 1.20
    vp1 = p["RBC0"] * (100.0 - h1) / h1
    return dict(hct0=h0, hct_threshold=h1, vp_at_threshold=vp1,
                plasma_lost_mL=p["VP0"] - vp1,
                plasma_lost_pct=(p["VP0"] - vp1) / p["VP0"] * 100.0)


def starling_table(p=None):
    p = calibrate(dict(p)) if p else P
    Pc, Pi = p["_Pc0"], p["PI0"]
    PIp, PIi = cop_plasma(p["TP_P0"]), cop_interstitium(p["TP_I0"])
    rows = []
    for s in [0.97, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55]:
        nfp = (Pc - Pi) - s * (PIp - PIi)
        rows.append((s, nfp, p["KF0"] * nfp, p["KF0"] * nfp - p["JL0"]))
    return rows


# -----------------------------------------------------------------------------
# 9.  SCENARIOS
# -----------------------------------------------------------------------------
SEC_TITRE = 55.0
TERT_TITRE = 2600.0


def build_scenarios(tp):
    R = Regimen
    return {
        "S01_primary_naive":            dict(abh0=0.0,        reg=R()),
        "S02_secondary_untreated":      dict(abh0=SEC_TITRE,  reg=R()),
        "S03_secondary_WHO_fluids":     dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp))),
        "S04_under_resuscitation_50":   dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp, scale=0.5))),
        "S05_over_resuscitation_200":   dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp, scale=2.0))),
        "S06_colloid_rescue":           dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp, scale=0.55),
                                                                    colloid=colloid_bolus(tp, mlkg=15.0, n=2))),
        "S07_albumin_rescue":           dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp, scale=0.55),
                                                                    albumin=[(tp, tp + 2.0, 15.0)])),
        "S08_prophylactic_platelets":   dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    plt_tx=[(tp, 38.0), (tp + 24, 38.0), (tp + 48, 38.0)])),
        "S09_steroid_at_presentation":  dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    steroid=(tp, 500.0, 3))),
        "S09b_steroid_at_fever_onset":  dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    steroid=(0.0, 500.0, 3))),
        "S10_antiviral_day2_illness":   dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    antiviral=(132.0, 450.0, 14))),
        "S11_antiviral_at_present":     dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    antiviral=(tp, 450.0, 14))),
        "S12_vaccine_seronegative":     dict(abh0=48.0,       reg=R()),
        "S13_vaccine_seropositive":     dict(abh0=3100.0,     reg=R()),
        "S14_tertiary_high_titre":      dict(abh0=TERT_TITRE, reg=R()),
        "S15_paracetamol_high_dose":    dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp),
                                                                    apap=(120.0, 1000.0, 28))),
        "S16_anti_NS1_mab":             dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp), ns1mab=0.35)),
        "S17_early_oral_rehydration":   dict(abh0=SEC_TITRE,  reg=R(oral_extra=[(tp - 72.0, tp + 96.0, 95.0)])),
        "S18_delayed_presentation_8h":  dict(abh0=SEC_TITRE,  reg=R(crystalloid=who_crystalloid(tp + 8.0))),
    }


def fluid_dose_response(tp, scales, abh0=SEC_TITRE):
    out = []
    for s in scales:
        t, Y, o = simulate(abh0=abh0, reg=Regimen(crystalloid=who_crystalloid(tp, scale=s)))
        r = summarise(t, Y, o); r["scale"] = s
        out.append(r)
    return out


def antiviral_timing(tp, starts, abh0=SEC_TITRE):
    out = []
    for st in starts:
        t, Y, o = simulate(abh0=abh0, reg=Regimen(crystalloid=who_crystalloid(tp),
                                                  antiviral=(st, 450.0, 18)))
        r = summarise(t, Y, o); r["start_h"] = st
        out.append(r)
    return out


def titre_sweep(titres):
    """
    The epidemiological experiment, and it must be run differently from the
    clinical scenarios: here every host receives the SAME MOSQUITO INOCULUM
    and the antibody titre is allowed to decide whether an infection even
    takes hold.  Starting these runs at fever onset would hand a protected
    host a viraemia it could never have reached.
    """
    out = []
    for a in titres:
        t, Y, o = simulate(abh0=a, reg=Regimen(), tend=720.0, v0=50.0)
        r = summarise(t, Y, o); r["titre"] = a
        out.append(r)
    return out


# -----------------------------------------------------------------------------
# 10.  MAIN
# -----------------------------------------------------------------------------
def baseline_selftest():
    """No virus: every state must sit still."""
    t, Y, o = simulate(abh0=0.0, reg=Regimen(), tend=336.0, v0=0.0)
    drift = dict(Vp=o["Vp"][-1] - P["VP0"], Vi=o["Vi"][-1] - P["VI0"],
                 TP=o["TP"][-1] - P["TP_P0"], Hct=o["Hct"][-1] - o["Hct"][0],
                 MAP=o["MAP"][-1] - o["MAP"][0], Jv=o["Jv"][-1] - P["JL0"],
                 PLT=o["PLT"][-1] - 250.0, Temp=o["Temp"][-1] - 37.0)
    return drift


def main():
    out = []
    W = out.append
    W("=" * 79)
    W(" DENGUE / SEVERE DENGUE QSP MODEL -- PYTHON/SciPy REFERENCE RUN")
    W(" 45 ODEs, LSODA -- every number below is solver output, not an assertion")
    W("=" * 79)

    W("")
    W("0. BASELINE SELF-TEST (uninfected host, 14 days)")
    W("-" * 79)
    dr = baseline_selftest()
    W("   " + "  ".join(f"d{k}={v:+.3f}" for k, v in dr.items()))

    ade = ade_curve()
    W("")
    W("A. ANTIBODY-DEPENDENT ENHANCEMENT -- THE ENTRY FACTOR E(titre)")
    W("-" * 79)
    W(f"   peak enhancement       E = {ade['peak_E']:.2f}  at reciprocal titre 1:{ade['peak_titre']:.0f}")
    W(f"   E crosses 1.0 again at   1:{ade['cross_titre']:.0f}   (net protection above this)")
    W(f"   E(1:20) {ade['E_at_20']:.2f}   E(1:40) {ade['E_at_40']:.2f}   "
      f"E(1:80) {ade['E_at_80']:.2f}   E(1:320) {ade['E_at_320']:.2f}")
    W(f"   E(1:1280) {ade['E_at_1280']:.3f}   E(1:2560) {ade['E_at_2560']:.3f}   "
      f"E(1:5120) {ade['E_at_5120']:.4f}")
    W("   Salje 2018 Nature: risk raised at 1:21-1:80, protection above 1:1280.")

    inf = infant_age_prediction()
    W("")
    W("B. INFANT DENGUE -- THE AGE OF PEAK SEVERITY IS A DIVISION PROBLEM")
    W("-" * 79)
    W(f"   cord titre 1:{inf['cord_titre']:.0f}, maternal IgG t1/2 = {inf['t_half_days']:.0f} d")
    W(f"   half-lives down to the peak-enhancement titre 1:{inf['peak_enh_titre']:.0f} = {inf['half_lives']:.2f}")
    W(f"   predicted age of maximum enhancement = {inf['age_days']:.0f} d = {inf['age_months']:.1f} months")
    W("   observed peak of infant DHF: 6-8 months (Kliks 1988; Chau 2009).")
    W("   Two inputs, no infant data: a cord titre and an IgG half-life.")

    ha = hct_arithmetic()
    W("")
    W("C. WHAT A 20 % HAEMATOCRIT RISE ACTUALLY COSTS")
    W("-" * 79)
    W(f"   baseline Hct {ha['hct0']:.1f} %  ->  WHO leak threshold {ha['hct_threshold']:.1f} %")
    W(f"   plasma volume at that haematocrit = {ha['vp_at_threshold']:.0f} mL")
    W(f"   plasma already gone = {ha['plasma_lost_mL']:.0f} mL = {ha['plasma_lost_pct']:.1f} % of the compartment")

    W("")
    W("D. THE STARLING TABLE -- sigma IS THE DISEASE")
    W("-" * 79)
    W("    sigma   net filtration P    J_v (mL/h)   J_v - lymph_baseline")
    for s, nfp, jv, net in starling_table():
        W(f"    {s:4.2f}      {nfp:8.3f} mmHg   {jv:9.0f}      {net:+9.0f}")
    W(f"   Kf = {P['KF0']:.1f} mL/(h*mmHg) is fixed by requiring baseline J_v to")
    W(f"   equal baseline lymph flow ({P['JL0']:.0f} mL/h); lymph ceiling {P['JLMAX']:.0f} mL/h.")

    # anchor runs -----------------------------------------------------------
    t0, Y0, o0 = simulate(abh0=SEC_TITRE, reg=Regimen())
    dv = defervescence(t0, o0)
    tp = presentation_time(t0, o0, n_warn=3)

    W("")
    W("E. SCENARIOS")
    W("-" * 79)
    W(f"   untreated secondary case defervesces at t = {dv:.1f} h (illness day {dv/24:.1f});")
    W(f"   it first shows three WHO warning signs at t = {tp:.0f} h (illness day {tp/24:.1f}),")
    W(f"   and that -- not a chosen hour -- is when every treatment scenario starts.")
    W("")
    scen = build_scenarios(tp)
    results = {}
    hdr = (f"   {'scenario':<30}{'Hct+%':>7}{'PPmin':>7}{'shk_h':>7}{'PLTmin':>7}"
           f"{'eff_mL':>8}{'lact':>6}{'AST':>6}{'sev':>7}")
    W(hdr); W("   " + "-" * (len(hdr) - 3))
    for name, cfg in scen.items():
        t, Y, o = simulate(**cfg)
        s = summarise(t, Y, o)
        results[name] = s
        W(f"   {name:<30}{s['hct_rise_pct']:>7.1f}{s['min_PP']:>7.1f}"
          f"{s['shock_hours']:>7.1f}{s['min_PLT']:>7.0f}{s['max_Vser']:>8.0f}"
          f"{s['max_lactate']:>6.2f}{s['max_AST']:>6.0f}{s['severity_index']:>7.3f}")

    W("")
    W("F. PRIMARY vs SECONDARY -- SAME VIRUS, SAME HOST, ONE TITRE APART")
    W("-" * 79)
    keys = ["fever_onset_h", "defervescence_h", "peak_viraemia_log10", "peak_NS1",
            "NS1_AUC", "IC_AUC", "peak_ABN", "max_TNF", "max_IL10", "max_CHYM",
            "min_GLX", "min_sigma", "max_Jv", "max_netleak", "hct_rise_pct",
            "plasma_loss_pct", "min_TP", "min_PP", "min_MAP", "max_Vser",
            "min_PLT", "min_WBC", "max_AST", "max_lactate", "shock_hours",
            "max_warnsigns", "severity_index"]
    a, b = results["S01_primary_naive"], results["S02_secondary_untreated"]
    W(f"   {'':<22}{'primary':>12}{'secondary':>12}{'ratio':>9}")
    for k in keys:
        va, vb = a[k], b[k]
        r = vb / va if abs(va) > 1e-9 else float('nan')
        W(f"   {k:<22}{va:>12.2f}{vb:>12.2f}{r:>9.2f}")

    W("")
    W("G. FLUID DOSE-RESPONSE -- THE OPTIMUM IS INTERIOR")
    W("-" * 79)
    scales = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
    fdr = fluid_dose_response(tp, scales)
    W(f"   {'x WHO':>6}{'IV mL':>9}{'shock h':>9}{'effusion':>10}{'lactate':>9}"
      f"{'bleed':>8}{'severity':>10}")
    for r in fdr:
        iv = sum((t1 - t0) * v for t0, t1, v in who_crystalloid(tp, scale=r["scale"]))
        W(f"   {r['scale']:>6.2f}{iv:>9.0f}{r['shock_hours']:>9.1f}"
          f"{r['max_Vser']:>10.0f}{r['max_lactate']:>9.2f}{r['max_bleed']:>8.3f}"
          f"{r['severity_index']:>10.3f}")
    best = min(fdr, key=lambda r: r["severity_index"])
    iv_best = sum((t1 - t0) * v for t0, t1, v in who_crystalloid(tp, scale=best["scale"]))
    W(f"   minimum severity at {best['scale']:.2f} x the WHO ladder = {iv_best:.0f} mL IV")

    W("")
    W("H. ANTIVIRAL START TIME -- WHY DAY-5 TRIALS CANNOT WORK")
    W("-" * 79)
    fo_sec = results["S02_secondary_untreated"]["fever_onset_h"]
    starts = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 48.0, 60.0, 72.0, 96.0, 120.0]
    avt = antiviral_timing(tp, starts)
    ref = results["S03_secondary_WHO_fluids"]
    W(f"   {'start h':>8}{'illness d':>11}{'NS1 AUC':>10}{'% untx':>8}{'IC AUC':>9}"
      f"{'Hct+%':>8}{'shock h':>9}{'severity':>10}{'% benefit':>11}")
    for r in avt:
        illd = r["start_h"] / 24.0
        ben = 100.0 * (ref["severity_index"] - r["severity_index"]) / ref["severity_index"]
        W(f"   {r['start_h']:>8.0f}{illd:>11.1f}{r['NS1_AUC']:>10.0f}"
          f"{100*r['NS1_AUC']/ref['NS1_AUC']:>8.1f}{r['IC_AUC']:>9.1f}"
          f"{r['hct_rise_pct']:>8.1f}{r['shock_hours']:>9.1f}"
          f"{r['severity_index']:>10.3f}{ben:>11.1f}")
    W(f"   untreated reference: NS1 AUC {ref['NS1_AUC']:.0f} ng/mL*h, "
      f"severity {ref['severity_index']:.3f}")

    W("")
    W("I. PRE-INFECTION TITRE vs OUTCOME -- A MODEL OUTPUT, NOT AN INPUT")
    W("-" * 79)
    titres = [0.0, 10.0, 20.0, 40.0, 55.0, 80.0, 160.0, 320.0, 640.0, 1280.0,
              2560.0, 5120.0, 10240.0]
    tsw = titre_sweep(titres)
    naive = tsw[0]["severity_index"]
    W(f"   {'titre':>10}{'E':>7}{'peakV':>8}{'NS1pk':>8}{'Hct+%':>7}{'PPmin':>7}"
      f"{'PLTmin':>8}{'severity':>10}{'vs naive':>10}")
    for a_, r in zip(titres, tsw):
        E = ade_factors(a_, 0.0, P)[0]
        W(f"   1:{a_:<8.0f}{E:>7.2f}{r['peak_viraemia_log10']:>8.2f}"
          f"{r['peak_NS1']:>8.0f}{r['hct_rise_pct']:>7.1f}{r['min_PP']:>7.1f}"
          f"{r['min_PLT']:>8.0f}{r['severity_index']:>10.3f}"
          f"{r['severity_index'] - naive:>+10.3f}")

    W("")
    W("J. COLLOID vs CRYSTALLOID -- THE ADVANTAGE IS PROPORTIONAL TO sigma")
    W("-" * 79)
    for mult, lbl in [(0.7, "moderate leak"), (1.0, "severe leak"), (1.3, "profound leak")]:
        pp = dict(P); pp["kGDEG"] = P["kGDEG"] * mult
        tA, YA, oA = simulate(abh0=SEC_TITRE, p=pp, reg=Regimen(crystalloid=who_crystalloid(tp)))
        tB, YB, oB = simulate(abh0=SEC_TITRE, p=pp,
                              reg=Regimen(crystalloid=who_crystalloid(tp, scale=0.55),
                                          colloid=colloid_bolus(tp, mlkg=15.0, n=2)))
        sA, sB = summarise(tA, YA, oA, pp), summarise(tB, YB, oB, pp)
        W(f"   {lbl:<15} sigma_min {sA['min_sigma']:.3f} | crystalloid: shock "
          f"{sA['shock_hours']:5.1f} h sev {sA['severity_index']:.3f} | colloid: shock "
          f"{sB['shock_hours']:5.1f} h sev {sB['severity_index']:.3f} | delta "
          f"{sB['severity_index']-sA['severity_index']:+.3f}")
    W("   Wills 2005 NEJM randomised 512 children and found NO overall difference")
    W("   between Ringer's lactate, dextran-70 and 6 % starch, but a colloid")
    W("   advantage confined to the narrowest-pulse-pressure stratum.  The model")
    W("   reproduces that pattern without being told: the colloid margin is small")
    W("   and grows monotonically with the depth of the sigma defect, because the")
    W("   oncotic term it adds enters the Starling equation multiplied by sigma.")

    W("")
    W("K. HOW MANY mL OF CRYSTALLOID BUY 1 mL OF RETAINED PLASMA")
    W("-" * 79)
    tN, YN, oN = simulate(abh0=SEC_TITRE, reg=Regimen())
    dtg = tN[1] - tN[0]
    i0, i1 = int(tp / dtg), int((tp + 12.0) / dtg)
    vp_ref = oN["Vp"][i1]
    W(f"   {'x WHO':>7}{'IV given 12 h':>15}{'Vp at +12 h':>13}{'gain mL':>10}{'mL per mL':>11}")
    for s in [0.5, 1.0, 1.5, 2.0, 3.0]:
        t, Y, o = simulate(abh0=SEC_TITRE, reg=Regimen(crystalloid=who_crystalloid(tp, scale=s)))
        given = Y[IX["CIN"], i1] - Y[IX["CIN"], i0] - (YN[IX["CIN"], i1] - YN[IX["CIN"], i0])
        gain = o["Vp"][i1] - vp_ref
        W(f"   {s:>7.2f}{given:>15.0f}{o['Vp'][i1]:>13.0f}{gain:>10.0f}"
          f"{given/gain if gain > 1 else float('nan'):>11.1f}")
    pn = dict(P); pn["SIGMIN"] = P["SIG0"]          # endothelium that cannot leak
    tn, Yn, on = simulate(abh0=SEC_TITRE, p=pn, reg=Regimen())
    tn2, Yn2, on2 = simulate(abh0=SEC_TITRE, p=pn,
                             reg=Regimen(crystalloid=who_crystalloid(tp)))
    given_n = Yn2[IX["CIN"], i1] - Yn2[IX["CIN"], i0] - (Yn[IX["CIN"], i1] - Yn[IX["CIN"], i0])
    gain_n = on2["Vp"][i1] - on["Vp"][i1]
    W(f"   control, same patient with sigma held at {P['SIG0']:.2f}: "
      f"{given_n:.0f} mL buys {gain_n:.0f} mL, i.e. {given_n/gain_n:.1f} mL per mL")
    W(f"   (no-IV reference plasma volume at the same instant: {vp_ref:.0f} mL)")

    W("")
    W("L. DEFERVESCENCE AND LEAK ONSET ARE ONE EVENT")
    W("-" * 79)
    t, Y, o = simulate(abh0=SEC_TITRE, reg=Regimen())
    dvt = defervescence(t, o)
    iv_ = int(np.argmax(o["V"])); iab = int(np.argmax(np.gradient(o["ABN"], t)))
    iic = int(np.argmax(o["IC"])); ilk = int(np.argmax(o["netser"]))
    ihc = int(np.argmax(o["Hct"])); ipp = int(np.argmin(o["PP"]))
    ins = int(np.argmax(o["NS1"]))
    for lbl, i in [("peak viraemia", iv_), ("peak NS1", ins),
                   ("steepest rise of neutralising IgG", iab),
                   ("peak immune-complex load", iic)]:
        W(f"   {lbl:<38} t = {t[i]:7.1f} h")
    W(f"   {'defervescence (T < 37.8 C)':<38} t = {dvt:7.1f} h")
    for lbl, i, extra in [("peak serosal leak (J_ser - drainage)", ilk,
                           f"{o['netser'][ilk]:.0f} mL/h"),
                          ("peak haematocrit", ihc, f"{o['Hct'][ihc]:.1f} %"),
                          ("narrowest pulse pressure", ipp, f"{o['PP'][ipp]:.1f} mmHg")]:
        W(f"   {lbl:<38} t = {t[i]:7.1f} h   ({extra})")
    W(f"   -> narrowest pulse pressure minus defervescence = {t[ipp]-dvt:+.1f} h;")
    W(f"      peak haematocrit minus defervescence          = {t[ihc]-dvt:+.1f} h;")
    W(f"      peak serosal filtration rate minus defervescence = {t[ilk]-dvt:+.1f} h.")
    W("   Temperature is computed from hypothalamic PGE2 driven by interferon;")
    W("   the leak is computed from a glycocalyx driven by NS1, TNF, VEGF and")
    W("   chymase.  The two chains share no equation after the antibody rise,")
    W("   yet the clinical nadir and the fever break land within a day of each")
    W("   other, because both are that antibody rise seen from a different side.")

    W("")
    W("L2. WHY PROPHYLACTIC PLATELET TRANSFUSION DOES NOTHING")
    W("-" * 79)
    tA, YA, oA = simulate(abh0=SEC_TITRE, reg=Regimen(crystalloid=who_crystalloid(tp)))
    tB, YB, oB = simulate(abh0=SEC_TITRE,
                          reg=Regimen(crystalloid=who_crystalloid(tp),
                                      plt_tx=[(tp, 38.0), (tp + 24, 38.0), (tp + 48, 38.0)]))
    pa, pb = severity_parts(tA, oA), severity_parts(tB, oB)
    W(f"   platelet nadir      {np.min(oA['PLT']):.0f}  ->  {np.min(oB['PLT']):.0f} x10^9/L "
      f"(+{np.min(oB['PLT'])-np.min(oA['PLT']):.0f} after three pools)")
    W(f"   bleeding index      {np.max(oA['bleed']):.3f}  ->  {np.max(oB['bleed']):.3f}")
    W(f"   bleeding sub-score  {pa['bleeding']:.3f}  ->  {pb['bleeding']:.3f}")
    W(f"   severity            {severity(tA, oA):.3f}  ->  {severity(tB, oB):.3f}")
    W("   Haemostasis is written as a PRODUCT of four requirements -- platelet")
    W("   count, fibrinogen, vessel-wall integrity and perfusion -- so the worst")
    W("   term sets the answer.  At the nadir the vessel-wall term is the worst,")
    W("   and platelets are not the vessel wall.  AAPT (Lye 2017, Lancet) found")
    W("   prophylactic transfusion did not prevent bleeding; the model agrees for")
    W("   a structural reason rather than by fitting the trial.")

    W("")
    W("L3. CORTICOSTEROID START TIME -- THE SAME CLIFF AS THE ANTIVIRAL")
    W("-" * 79)
    W(f"   {'start h':>8}{'illness d':>11}{'peak TNF':>10}{'min GLX':>9}{'Hct+%':>8}"
      f"{'shock h':>9}{'NS1 AUC':>10}{'severity':>10}")
    for st in [0.0, 12.0, 24.0, 36.0, 48.0, 60.0, 72.0, 96.0]:
        t_, Y_, o_ = simulate(abh0=SEC_TITRE,
                              reg=Regimen(crystalloid=who_crystalloid(tp),
                                          steroid=(st, 500.0, 3)))
        s_ = summarise(t_, Y_, o_)
        W(f"   {st:>8.0f}{st/24:>11.1f}{s_['max_TNF']:>10.1f}{s_['min_GLX']:>9.3f}"
          f"{s_['hct_rise_pct']:>8.1f}{s_['shock_hours']:>9.1f}{s_['NS1_AUC']:>10.0f}"
          f"{s_['severity_index']:>10.3f}")
    W(f"   reference without steroid: severity {ref['severity_index']:.3f}, "
      f"NS1 AUC {ref['NS1_AUC']:.0f}")
    W("   Methylprednisolone helps in this model ONLY if it is given before the")
    W("   glycocalyx has integrated its mediator load.  Tam 2012 (Clin Infect Dis)")
    W("   enrolled within 72 h of fever onset and found no benefit; the model")
    W("   agrees at every start time a trial could realistically achieve, and")
    W("   disagrees only in a window that ends before patients seek care.")

    W("")
    W("M. WHICH LIMB OF THE LEAK MATTERS -- LESION-BY-LESION")
    W("-" * 79)
    variants = [("full model", {}),
                ("sigma frozen at 0.97", dict(SIGMIN=P["SIG0"])),
                ("Kf frozen at baseline", dict(ALPHLP=0.0)),
                ("protein sieving off (PS=0)", dict(PSPROT=0.0)),
                ("serosal drain uncapped", dict(DSERMAX=600.0)),
                ("no baroreflex (SVR fixed)", dict(BARGAIN=0.0, SVRMAX=1.0))]
    W(f"   {'variant':<28}{'plasma loss %':>14}{'min TP':>9}{'PPmin':>8}"
      f"{'effusion':>10}{'severity':>10}")
    for lbl, mod in variants:
        pp = dict(P); pp.update(mod)
        tt, YY, oo = simulate(abh0=SEC_TITRE, p=pp, reg=Regimen())
        ss = summarise(tt, YY, oo, pp)
        W(f"   {lbl:<28}{ss['plasma_loss_pct']:>14.1f}{ss['min_TP']:>9.2f}"
          f"{ss['min_PP']:>8.1f}{ss['max_Vser']:>10.0f}{ss['severity_index']:>10.3f}")

    text = "\n".join(out)
    print(text)
    with open("denv_reference_output.txt", "w") as fh:
        fh.write(text + "\n")

    payload = dict(
        baseline_drift=dr,
        ade={k: v for k, v in ade.items() if k not in ("titres", "E")},
        ade_curve=dict(titre=[float(x) for x in ade["titres"][::28]],
                       E=[float(x) for x in ade["E"][::28]]),
        infant=inf, hct=ha, t_present=tp,
        scenarios=results, fluid_dose_response=fdr, antiviral_timing=avt,
        titre_sweep=tsw,
    )
    with open("denv_scenario_results.json", "w") as fh:
        json.dump(payload, fh, indent=1)
    return payload


if __name__ == "__main__":
    main()
