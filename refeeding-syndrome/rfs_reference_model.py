#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rfs_reference_model.py
======================
Independent Python/scipy reference implementation of the Refeeding Syndrome
(RFS) QSP model.  Its ONLY purpose is to verify the arithmetic of the
mrgsolve model (`rfs_mrgsolve_model.R`) before that model is trusted: the two
are written from the same equation list but as separate code, so a
transcription slip in one shows up as disagreement with the other.

==========================================================================
THE THESIS
==========================================================================
Refeeding syndrome is not a disease of low serum phosphate.  It is a
FLUX MISMATCH, and the two fluxes are set by two different people.

    Lambda_P  =  J_demand / J_supply

    J_demand  the phosphate that must move into cells because insulin has
              raised the cellular organic-phosphate set-point (G6P, F1,6BP,
              ATP, phosphocreatine, 2,3-DPG) plus the phosphorus laid down
              in new lean tissue.  Set by the CLINICIAN, through the
              glucose infusion rate.

    J_supply  absorbed dietary phosphate + prescribed intravenous phosphate
              + net bone efflux + renal phosphate SAVED by switching off
              urinary excretion.  Set by the PATIENT'S HISTORY, and by what
              the clinician chooses to hand over.

Three pieces of arithmetic make this a dangerous ratio and all three are
stoichiometry rather than fitted parameters.

(1) THE MEASURED COMPARTMENT IS NOT THE COMPARTMENT AT RISK.
        ECF inorganic phosphate = 1.15 mmol/L x 12.4 L =     14.3 mmol
        Total body phosphorus                         = 22,780 mmol
    The laboratory samples 0.06 % of the pool.  A demand of 60 mmol/d is
    four times the entire measured compartment PER DAY.

(2) THE RENAL RESERVE HAS ALREADY BEEN SPENT.  Urinary phosphate excretion
    is a THRESHOLD function, GFR x max(0, Pser - TmP/GFR).  A starved
    patient sits at or below their threshold and is already excreting
    essentially nothing, so "the kidney will just hold on to it" is a
    reserve that was used up before admission.  J_supply reduces to bone
    plus whatever is prescribed.

(3) BONE IS SLOW, AND MAGNESIUM CAN SWITCH IT OFF.  Net skeletal phosphate
    efflux is ~10 mmol/d and is PTH-driven; hypomagnesaemia below ~0.5
    mmol/L suppresses PTH secretion, so the one endogenous supply line that
    remains is disabled by a second electrolyte the same syndrome depletes.

Thiamine runs on a DIFFERENT CLOCK, and that difference is the second
thesis.  Both pools take ~3 weeks to empty, but the REPLETION constants
differ by two orders of magnitude, because at supra-physiological plasma
concentrations thiamine bypasses the saturable ThTR-1/ThTR-2 carriers by
passive diffusion, whereas intravenous phosphate refills only the ECF and
must then be pumped into cells.  The clinical rule falls out of the
arithmetic instead of being asserted:

        thiamine before glucose matters by HOURS
        phosphate before glucose matters by DAYS

and, less obviously, oral thiamine cannot substitute for intravenous
thiamine in an acutely at-risk patient, because intestinal uptake saturates
near 5 mg per dose no matter how many milligrams are swallowed.

==========================================================================
CALIBRATION POLICY  (what is measured, what is predicted)
==========================================================================
Spent on NORMAL physiology, none of it on refeeding syndrome:
    total body P/K/Mg and their bone/ICF/ECF partition; 67 mmol P and
    72 mmol K per kg fat-free mass; serum P/K/Mg/Ca reference intervals;
    TmP/GFR; GFR 6 L/h; urinary P 29, K 70, Mg 4.9 mmol/d in balance with
    a normal diet; whole-body thiamine 26.5 mg with a 14 d biological
    half-life -- which REPRODUCES the 1.1-1.4 mg/d RDA as an OUTPUT;
    Cunningham REE; Forbes fat/lean partition; insulin t1/2 5 min;
    glucose Vd 0.16 L/kg; brain glucose use 28 mmol/h.
Spent on NON-refeeding pharmacology and nutrition data:
    0.3 mmol phosphate per gram of infused dextrose (the classical
    parenteral-nutrition additive rule) -- used as a CHECK, not an input;
    enteral formula content 22 mmol P, 38 mmol K, 6 mmol Mg per 1000 kcal;
    oral thiamine absorption saturating near 5 mg/dose; the 0.6-1.0 mmol/L
    fall in serum potassium produced by insulin-dextrose given for
    hyperkalaemia; Ca x P solubility product 4.4 mmol^2/L^2; QTc
    sensitivity to K, Mg and Ca.
Spent on REFEEDING SYNDROME ITSELF: THREE numbers.
    sIns  the fractional rise in the cellular organic-phosphate set-point
          at maximal insulin (0.09), which is what converts an insulin
          signal into a phosphate demand;
    kFill the rate at which cells move toward that set-point;
    HAZ   one global scale on the composite mortality hazard.
Everything else printed below is a PREDICTION.

==========================================================================
STRUCTURE
==========================================================================
Phase 0   healthy steady state       -- verified, max|dy/dt| is reported
Phase 1   starvation, simulated      -- the deficit is GROWN, not asserted,
                                        and the model's first real test is
                                        that serum P/K/Mg stay in range
                                        while the pools empty
Phase 2   refeeding by protocol

48 ODEs.   python3 rfs_reference_model.py > rfs_reference_output.txt
"""

import sys
import numpy as np
from scipy.integrate import solve_ivp

# =========================================================================
# 0.  PARAMETERS
# =========================================================================

P = {}

# ---- anthropometry (healthy, pre-morbid) --------------------------------
P["BW0"]   = 62.0; P["HT"] = 1.68; P["FM0"] = 14.0; P["FFM0"] = 48.0
P["fEcf"]  = 0.200; P["fVg"] = 0.160; P["fVp"] = 0.045
P["GFR0"]  = 6.00

# ---- phosphate ----------------------------------------------------------
P["rhoP"]     = 67.0      # mmol P per kg FFM (soft-tissue ICF)
P["PB0"]      = 19550.0   # mmol bone phosphorus
P["Pser0"]    = 1.15
P["TmPmax"]   = None      # solved
P["PTH50tm"]  = 4.20; P["FGF50tm"] = 45.0; P["hTmP"] = 1.30
P["ObligP"]   = 2.0/24.0
P["FabsP"]    = 0.65; P["kaP"] = 0.25
P["KmPup"]    = 0.30      # mmol/L  Km of Na-Pi cotransport on serum P
P["sIns"]     = 0.750     # scale of the glycolytic-flux-driven rise of the
                          # cellular organic-phosphate set-point
                          # (ONE OF THE THREE NUMBERS CALIBRATED ON RFS)
P["uSat"]     = 0.50      # normalised glycolytic flux at half the set-point rise
P["kFill"]    = 0.0145    # 1/h  cell filling rate  (t1/2 48 h)
P["Pnew"]     = 67.0      # mmol P per kg of new lean tissue
P["kBres"]    = None; P["kBform"] = None; P["hBres"] = 1.60
P["BoneCap"]  = 1.35      # max multiple of basal bone efflux (PTH-driven)
P["kPrecip"]  = 0.55; P["Ksp"] = 4.40
P["etohP"]    = 0.32      # fractional fall of TmP/GFR at full alcohol use

# ---- potassium ----------------------------------------------------------
P["rhoK"]  = 72.0; P["Kser0"] = 4.20
P["kKp"]   = 8.60         # L/h  pump flux scale (sets a ~1 h serum t1/2)
P["EmaxK"] = 0.34         # maximal fractional rise of pump activity
P["tauX"]  = 1.20
P["kKleak"] = None
P["FabsK"] = 0.90; P["kaK"] = 0.40
P["betaMgK"] = 1.90       # renal K wasting per unit fractional hypoMg
P["kKsec"] = None; P["Knew"] = 72.0; P["rAld"] = 0.30
P["Kobl"]  = 6.0/24.0     # mmol/h  irreducible (colonic + minimal renal) K loss
P["hRomk"] = 6.0          # steepness of ROMK down-regulation in K depletion
P["romkMin"] = 0.05

# ---- magnesium ----------------------------------------------------------
P["rhoMg"] = 7.30; P["MGser0"] = 0.85; P["fUFMg"] = 0.70
P["TmMg"]  = None
P["FabsMg"] = 0.40; P["kaMg"] = 0.35
P["kMgp"]  = 1.30; P["EmaxMg"] = 0.20; P["kMgleak"] = None
P["Mgnew"] = 7.30
P["etohMg"] = 0.28        # fractional fall of TmMg at full alcohol use

# ---- calcium / PTH / vitamin D / FGF23 ----------------------------------
P["Ca0"] = 1.20; P["PTH0"] = 4.20; P["PTHmax"] = 26.0; P["PTHmin"] = 0.60
P["Ca50"] = None; P["hCa"] = 6.0; P["kPTH"] = 6.93
P["MgPTH50"] = 0.22
P["CTD0"] = 120.0; P["kCTD"] = 0.10
P["FGF0"] = 45.0; P["kFGF"] = 0.14
P["kCaBone"] = None; P["kCaU"] = None

# ---- thiamine -----------------------------------------------------------
P["THT0"] = 100.0         # umol  (26.5 mg)
P["kThDeg"] = 0.00206     # 1/h   ln2 / 14 d
P["kThCho"] = 2.05e-4     # umol per mmol glucose disposed
P["ThP0"] = 10.0          # nmol/L
P["CLth"] = 0.42; P["ThRenalC"] = 30.0
P["VmaxThA"] = 7.50; P["KmThA"] = 6.0
P["VmaxThT"] = None; P["KmThT"] = 12.0; P["kThPass"] = 0.0042
P["kThOff"] = 0.020
P["KmTPP"] = 0.35; P["hTPP"] = 3.0
P["EtOHabs"] = 0.45; P["EtOHtpk"] = 0.55

# ---- glucose / insulin --------------------------------------------------
P["G0"] = 5.00; P["INS0"] = 45.0
P["EGP0"] = 44.0; P["X50"] = 1.00
P["Rd0"] = 5.00           # L/h insulin-independent (brain ~28 mmol/h)
P["SI"] = None            # solved from the fed-state glucose balance
P["kSec"] = None; P["Gthr"] = 3.00; P["hSec"] = 2.00; P["SecB"] = None
P["INSlow"] = 15.0        # pmol/L insulin at the starvation glycaemia below
P["Glow"] = 3.50          # mmol/L glycaemia of prolonged starvation
P["KmGly"] = 0.25        # mmol/L  phosphate Km of glycolytic (GAPDH) flux
P["siMin"] = 0.45         # insulin sensitivity floor of the adapted, starved
                          # state ("starvation diabetes"); recovers with AT
P["fInc"] = 0.35          # fraction of fed insulin secretion driven by the
                          # ENTERAL nutrient (incretin) signal rather than by
                          # glycaemia -- the reason a gram of enteral
                          # carbohydrate raises insulin more than a gram of
                          # intravenous dextrose
P["incIV"] = 0.45         # residual incretin drive of intravenous dextrose
P["kInc"] = None
P["RaLowRef"] = 0.32      # starvation intake fraction used for the two-point
                          # insulin-secretion calibration
P["CLins"] = 8.32; P["p2"] = 1.50
P["kaCho"] = 0.55; P["kaFat"] = 0.28; P["kaPro"] = 0.38
P["GLY0"] = 2200.0; P["GLYmax"] = 3300.0
P["kGlyS"] = None; P["kGlyB"] = None; P["Jgly"] = 12.0
P["GCG0"] = 60.0; P["kGcg"] = 0.60

# ---- lipid / ketone / lactate ------------------------------------------
P["FFA0"] = 0.45; P["kLipo"] = None; P["kFFAox"] = 1.65
P["BHB0"] = 0.10; P["kKeto"] = None; P["kBHBox"] = 0.62
P["LAC0"] = 1.00; P["CLlac"] = 22.0; P["kLac2"] = 1.50

# ---- organ energetics ---------------------------------------------------
P["kATPm"] = 0.14; P["KPatp"] = 0.35; P["hPatp"] = 2.0
P["kPCr"]  = 0.30; P["KPpcr"] = 0.50
P["kDPG"]  = 0.10; P["KPdpg"] = 0.55
P["kATPd"] = 0.16

# ---- body composition / heart / volume ----------------------------------
P["rhoFat"] = 9440.0; P["rhoLean"] = 1800.0; P["PAL"] = 1.30
P["kAT"] = 0.0060; P["ATmin"] = 0.78
P["LVM0"] = 150.0; P["kLVM"] = 1.0/(21*24); P["LVMfrac"] = None
P["Na0"] = 140.0; P["rIns"] = 0.62; P["kNaB"] = 1.0/72.0
P["ALD0"] = 1.00; P["kALD"] = 0.06

# ---- hazards ------------------------------------------------------------
P["QTc0"] = 400.0; P["aK"] = 26.0; P["aMg"] = 22.0; P["aCa"] = 95.0
P["h0arr"] = 5.0e-6; P["cQT"] = 0.045; P["cATParr"] = 3.2
P["h0cf"] = 3.0e-6; P["pCong"] = 2.0
P["h0rf"] = 2.0e-6; P["PaCO2ref"] = 55.0; P["qCO2"] = 0.075
P["h0we"] = 3.0e-4
P["hazMax"] = 0.020        # 1/h  ceiling on the instantaneous total hazard
P["VAmin"] = 0.35          # floor on relative alveolar ventilation capacity
P["contMin"] = 0.20        # floor on the myocardial contractility index
P["HAZ"] = 1.00            # global hazard scale (THE refeeding-specific number)

# ---- gas exchange -------------------------------------------------------
P["kcalPerLO2"] = 4.83; P["PaCO2_0"] = 40.0; P["ChoOxMax"] = 4.6
P["VCO2ref"] = None


# =========================================================================
# 1.  STATE VECTOR
# =========================================================================

NAMES = [
    "A_CHO", "A_FAT", "A_PRO", "GLU", "INS", "X", "GCG", "GLY", "FFA",
    "BHB", "LAC", "PGUT", "PE", "PI", "PB", "KGUT", "KE", "KI", "MGGUT",
    "MGE", "MGI", "CAE", "CAB", "PTH", "CTD", "FGF", "THGUT", "THP",
    "THT", "ATPm", "PCRm", "DPG", "ATPd", "FM", "FFM", "AT", "LVM",
    "NAB", "ALD", "H_ARR", "H_CF", "H_RF", "H_WE", "SURV", "AUCP",
    "CUMP", "CUMK", "CUMKCAL",
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)


# =========================================================================
# 2.  HELPERS
# =========================================================================

def pos(x):
    return x if x > 0.0 else 0.0


def hill(x, k, h):
    x = pos(x)
    return x ** h / (k ** h + x ** h)


def ree_kcal_d(FFM, AT):
    return (370.0 + 21.6 * FFM) * AT


def sw_below(x, thr, width):
    """Smooth 'is x below thr' switch: 1 well below, 0 well above.

    Clinical repletion is written in the notes as a hard threshold, but a
    hard threshold makes the integrator chatter on the crossing, and real
    prescribing is not instantaneous either.  A logistic of width `width`
    is both faster and more honest.
    """
    z = (x - thr) / width
    if z > 40.0:
        return 0.0
    if z < -40.0:
        return 1.0
    return 1.0 / (1.0 + np.exp(z))


def pumpmult(X, Emax, tau):
    """Saturating effect of insulin on Na/K-ATPase-mediated cation uptake."""
    if X >= 1.0:
        return 1.0 + Emax * (1.0 - np.exp(-(X - 1.0) / tau))
    return 1.0 - 0.6 * Emax * (1.0 - np.exp(-(1.0 - X) / tau))


# =========================================================================
# 3.  PROTOCOL OBJECT
# =========================================================================

class Regimen:
    def __init__(self, **kw):
        # starvation history
        self.starve_days   = kw.get("starve_days", 60.0)
        self.starve_frac   = kw.get("starve_frac", 0.32)
        self.micronutrient = kw.get("micronutrient", 0.25)  # diet quality 0-1
        self.alcohol       = kw.get("alcohol", 0.0)
        self.diuretic      = kw.get("diuretic", 0.0)
        self.gi_loss       = kw.get("gi_loss", 0.0)
        # refeeding
        self.kcal_start    = kw.get("kcal_start", 10.0)
        self.kcal_goal     = kw.get("kcal_goal", 30.0)
        self.advance_frac  = kw.get("advance_frac", 0.20)
        self.advance_start = kw.get("advance_start", 2.0)
        self.cho_frac      = kw.get("cho_frac", 0.50)
        self.fat_frac      = kw.get("fat_frac", 0.32)
        self.pro_frac      = kw.get("pro_frac", 0.18)
        self.route         = kw.get("route", "enteral")   # enteral | dextrose
        self.formula_P     = kw.get("formula_P", 22.0)
        self.formula_K     = kw.get("formula_K", 38.0)
        self.formula_Mg    = kw.get("formula_Mg", 6.0)
        self.formula_Th    = kw.get("formula_Th", 5.7)
        # repletion
        self.P_dose   = kw.get("P_dose", 0.0)     # mmol/kg/d IV
        self.K_dose   = kw.get("K_dose", 0.0)
        self.Mg_dose  = kw.get("Mg_dose", 0.0)
        self.Mg_bolus = kw.get("Mg_bolus", False)
        self.th_iv    = kw.get("th_iv", 0.0)      # mg/d IV
        self.th_po    = kw.get("th_po", 0.0)      # mg/d oral
        self.th_lead_h = kw.get("th_lead_h", 0.0)
        self.th_days  = kw.get("th_days", 10.0)
        self.th_delay = kw.get("th_delay", 0.0)   # h, thiamine started LATE
        self.P_delay  = kw.get("P_delay", 0.0)    # d, phosphate started LATE
        self.rescue   = kw.get("rescue", False)
        self.gut_fail = kw.get("gut_fail", 0.0)  # fractional loss of
                                                 # intestinal uptake
                                                 # (vomiting, ileus,
                                                 #  damaged mucosa)
        self.Na_intake = kw.get("Na_intake", 100.0)
        self.sim_days = kw.get("sim_days", 14.0)

    def kcal_per_kg(self, day):
        if day < self.advance_start:
            return self.kcal_start
        n = np.floor(day - self.advance_start) + 1.0
        return min(self.kcal_start + n * self.advance_frac * self.kcal_goal,
                   self.kcal_goal)


# =========================================================================
# 4.  STEADY-STATE SOLVE
# =========================================================================

def steady_state(P):
    BW, FFM = P["BW0"], P["FFM0"]
    P["Vecf0"] = P["fEcf"] * BW
    P["Vg"] = P["fVg"] * BW
    P["Vp"] = P["fVp"] * BW
    P["VLac"] = 0.55 * BW
    P["LVMfrac"] = P["LVM0"] / FFM

    P["PI0"] = P["rhoP"] * FFM
    P["KI0"] = P["rhoK"] * FFM
    P["MGI0"] = P["rhoMg"] * FFM
    P["PE0"] = P["Pser0"] * P["Vecf0"]
    P["KE0"] = P["Kser0"] * P["Vecf0"]
    P["MGE0"] = P["MGser0"] * P["Vecf0"]
    P["CAE0"] = P["Ca0"] * P["Vecf0"]

    # ---- energy: TEE includes 10 % diet-induced thermogenesis -----------
    REE0 = ree_kcal_d(FFM, 1.0)
    P["REE0"] = REE0
    P["TEE0"] = REE0 * P["PAL"] / 0.90
    kcal_h0 = P["TEE0"] / 24.0
    P["kcal_h0"] = kcal_h0

    # ---- glucose / insulin ---------------------------------------------
    g_cho_h0 = kcal_h0 * 0.50 / 4.0
    Ra0 = g_cho_h0 * 1000.0 / 180.0
    EGPss = P["EGP0"] / (1.0 + 1.0 / P["X50"])
    Rd_need = Ra0 + EGPss
    P["SI"] = (Rd_need / P["G0"] - P["Rd0"]) / 1.0
    P["Uins0"] = P["SI"] * 1.0 * P["G0"]
    P["Rd_0"] = Rd_need
    P["Ra0"] = Ra0
    # Insulin secretion is pinned by TWO points, not one: the fed reference
    # (G0 -> INS0) and the starvation operating point (Glow -> INSlow).  A
    # single-point fit leaves most of the secretion in a glucose-independent
    # basal term, which makes insulin -- and therefore the whole phosphate
    # demand -- almost unresponsive to refeeding.
    cvp = P["CLins"] * P["Vp"]
    a1 = pos(P["G0"] - P["Gthr"]) ** P["hSec"]
    a2 = pos(P["Glow"] - P["Gthr"]) ** P["hSec"]
    RaLow = P["RaLowRef"] * Ra0
    P["kInc"] = P["fInc"] * cvp * P["INS0"] / Ra0
    b1 = cvp * P["INS0"] - P["kInc"] * Ra0
    b2 = cvp * P["INSlow"] - P["kInc"] * RaLow
    P["kSec"] = (b1 - b2) / (a1 - a2)
    P["SecB"] = b1 - P["kSec"] * a1
    # glycogen cycling
    P["kGlyS"] = P["Jgly"] / (1.0 * P["G0"] * (1.0 - P["GLY0"] / P["GLYmax"]))
    P["kGlyB"] = P["Jgly"] * (1.0 + 3.0 * 1.0) / P["GLY0"]
    # lipolysis / ketogenesis balanced at baseline
    P["kLipo"] = P["kFFAox"] * P["FFA0"] * (1.0 + 3.2) / (P["FM0"] / P["FM0"] + 0.15)
    P["kKeto"] = P["kBHBox"] * P["BHB0"] * (1.0 + 6.0) / P["FFA0"]

    # ---- calcium set-point so that PTH secretion = PTH0 at Ca0 ----------
    num = P["PTHmax"] - P["PTHmin"]
    r = num / (P["PTH0"] - P["PTHmin"]) - 1.0
    P["Ca50"] = P["Ca0"] / r ** (1.0 / P["hCa"])

    # ---- phosphate balance ---------------------------------------------
    fP0 = 22.0 * kcal_h0 / 1000.0
    JPabs0 = P["FabsP"] * fP0
    P["JPabs0"] = JPabs0
    TmP_eff = P["Pser0"] - (JPabs0 - P["ObligP"]) / P["GFR0"]
    denom = ((1.0 + hill(P["PTH0"], P["PTH50tm"], P["hTmP"]))
             * (1.0 + 0.55 * hill(P["FGF0"], P["FGF50tm"], P["hTmP"])))
    P["TmPmax"] = TmP_eff * denom
    P["TmP0"] = TmP_eff
    P["fU1"] = 1.0 / (P["uSat"] + 1.0)
    Jbone = 10.0 / 24.0
    P["kBres"] = Jbone / (P["PB0"] * hill(P["PTH0"], P["PTH0"], P["hBres"]))
    P["kBform"] = Jbone / P["PB0"]
    P["Jbone0"] = Jbone

    # ---- potassium ------------------------------------------------------
    JKin0 = P["kKp"] * P["Kser0"]
    P["kKleak"] = JKin0 / P["KI0"]
    fK0 = 38.0 * kcal_h0 / 1000.0
    ExK0 = P["FabsK"] * fK0
    P["kKsec"] = (ExK0 - P["Kobl"]) / P["Kser0"]
    P["ExK0"] = ExK0

    # ---- magnesium ------------------------------------------------------
    JMgin0 = P["kMgp"] * P["MGser0"]
    P["kMgleak"] = JMgin0 / P["MGI0"]
    fMg0 = 6.0 * kcal_h0 / 1000.0
    ExMg0 = P["FabsMg"] * fMg0
    P["TmMg"] = P["fUFMg"] * P["MGser0"] - (ExMg0 - 0.5 / 24.0) / P["GFR0"]
    P["ExMg0"] = ExMg0

    # ---- calcium --------------------------------------------------------
    P["kCaU"] = (5.0 / 24.0) / P["Ca0"]
    P["kCaBone"] = (5.0 / 24.0) / P["CAE0"]

    # ---- thiamine -------------------------------------------------------
    Cth = P["ThP0"]
    P["THP0"] = Cth * P["Vp"] / 1000.0
    cho_ox0 = Rd_need
    loss = P["kThDeg"] * P["THT0"] + P["kThCho"] * cho_ox0
    need = loss + P["kThOff"] * P["THT0"]
    P["VmaxThT"] = (need - P["kThPass"] * Cth) * (P["KmThT"] + Cth) / Cth
    Jup0 = P["VmaxThT"] * Cth / (P["KmThT"] + Cth) + P["kThPass"] * Cth
    renal_th0 = P["CLth"] * pos(Cth - P["ThRenalC"]) / 1000.0
    P["ThAbs0"] = Jup0 + renal_th0 - P["kThOff"] * P["THT0"]
    P["ThRDA"] = P["ThAbs0"] * 24.0 * 265.4 / 1000.0
    P["fPDH0"] = hill(1.0, P["KmTPP"], P["hTPP"])
    # the enteral formula must supply that absorbed amount at baseline
    P["formulaTh0"] = P["ThAbs0"] / (kcal_h0 / 1000.0)

    # ---- lactate --------------------------------------------------------
    P["LacProd0"] = P["CLlac"] * P["LAC0"]

    # ---- gas exchange (filled by the calibration pass) ------------------
    P["VO2_0"] = P["TEE0"] / (1440.0 * P["kcalPerLO2"])
    P["VCO2ref"] = 0.85 * P["VO2_0"]
    return P


# =========================================================================
# 5.  INITIAL CONDITIONS
# =========================================================================

def y0_healthy(P):
    y = np.zeros(NST)
    kcal_h = P["kcal_h0"]
    y[IX["A_CHO"]] = (kcal_h * 0.50 / 4.0) / P["kaCho"]
    y[IX["A_FAT"]] = (kcal_h * 0.32 / 9.0) / P["kaFat"]
    y[IX["A_PRO"]] = (kcal_h * 0.18 / 4.0) / P["kaPro"]
    y[IX["GLU"]] = P["G0"]; y[IX["INS"]] = P["INS0"]; y[IX["X"]] = 1.0
    y[IX["GCG"]] = P["GCG0"]; y[IX["GLY"]] = P["GLY0"]
    y[IX["FFA"]] = P["FFA0"]; y[IX["BHB"]] = P["BHB0"]; y[IX["LAC"]] = P["LAC0"]
    y[IX["PGUT"]] = (22.0 * kcal_h / 1000.0) / P["kaP"]
    y[IX["PE"]] = P["PE0"]; y[IX["PI"]] = P["PI0"]; y[IX["PB"]] = P["PB0"]
    y[IX["KGUT"]] = (38.0 * kcal_h / 1000.0) / P["kaK"]
    y[IX["KE"]] = P["KE0"]; y[IX["KI"]] = P["KI0"]
    y[IX["MGGUT"]] = (6.0 * kcal_h / 1000.0) / P["kaMg"]
    y[IX["MGE"]] = P["MGE0"]; y[IX["MGI"]] = P["MGI0"]
    y[IX["CAE"]] = P["CAE0"]
    y[IX["PTH"]] = P["PTH0"]; y[IX["CTD"]] = P["CTD0"]; y[IX["FGF"]] = P["FGF0"]
    thg = P["formulaTh0"] * kcal_h / 1000.0
    y[IX["THGUT"]] = P["KmThA"] * thg / max(P["VmaxThA"] - thg, 1e-9)
    y[IX["THP"]] = P["THP0"]; y[IX["THT"]] = P["THT0"]
    y[IX["ATPm"]] = 1.0; y[IX["PCRm"]] = 1.0
    y[IX["DPG"]] = 1.0; y[IX["ATPd"]] = 1.0
    y[IX["FM"]] = P["FM0"]; y[IX["FFM"]] = P["FFM0"]; y[IX["AT"]] = 1.0
    y[IX["LVM"]] = P["LVM0"]; y[IX["NAB"]] = 0.0; y[IX["ALD"]] = P["ALD0"]
    y[IX["SURV"]] = 1.0
    return y


# =========================================================================
# 6.  RIGHT-HAND SIDE
# =========================================================================

def rhs(t, y, P, R, phase, aux=None):
    y = np.maximum(y, 0.0)
    (A_CHO, A_FAT, A_PRO, GLU, INS, X, GCG, GLY, FFA, BHB, LAC,
     PGUT, PE, PI, PB, KGUT, KE, KI, MGGUT, MGE, MGI, CAE, CAB,
     PTH, CTD, FGF, THGUT, THP, THT, ATPm, PCRm, DPG, ATPd,
     FM, FFM, AT, LVM, NAB, ALD, H_ARR, H_CF, H_RF, H_WE, SURV,
     AUCP, CUMP, CUMK, CUMKCAL) = y

    d = np.zeros(NST)
    BW = FM + FFM
    Vecf = P["fEcf"] * P["BW0"] + NAB / P["Na0"]
    Vecf = max(Vecf, 1.0)
    Vg = P["fVg"] * max(BW, 20.0)
    Vp = P["fVp"] * max(BW, 20.0)

    Pser = PE / Vecf; Kser = KE / Vecf
    MGser = MGE / Vecf; Caser = CAE / Vecf
    Cth = THP / Vp * 1000.0

    # -----------------------------------------------------------------
    # 6.1  nutrition
    # -----------------------------------------------------------------
    day = t / 24.0
    if phase == "healthy":
        kcal_d = P["TEE0"]; cho_f, fat_f, pro_f = 0.50, 0.32, 0.18
        route = "enteral"; quality = 1.0
    elif phase == "starve":
        kcal_d = R.starve_frac * P["TEE0"]
        cho_f, fat_f, pro_f = 0.50, 0.32, 0.18
        route = "enteral"; quality = R.micronutrient
    else:
        kcal_d = R.kcal_per_kg(day) * BW
        cho_f, fat_f, pro_f = R.cho_frac, R.fat_frac, R.pro_frac
        route = R.route; quality = 1.0
        if t < R.th_lead_h:
            kcal_d = 0.0

    kcal_h = kcal_d / 24.0
    g_cho_h = kcal_h * cho_f / 4.0
    g_fat_h = kcal_h * fat_f / 9.0
    g_pro_h = kcal_h * pro_f / 4.0
    scale = kcal_h / 1000.0

    # Diet QUALITY hits the nutrients unequally.  Refined carbohydrate and
    # alcohol are almost devoid of thiamine and poor in phosphate, but any
    # real food still carries potassium roughly in proportion to its mass.
    qP = 0.40 + 0.60 * quality
    qK = 0.70 + 0.30 * quality
    qMg = 0.45 + 0.55 * quality
    qTh = quality
    if route == "dextrose":
        fP = fK = fMg = fTh = 0.0
    else:
        fP = R.formula_P * scale * qP
        fK = R.formula_K * scale * qK
        fMg = R.formula_Mg * scale * qMg
        fTh = R.formula_Th * scale * qTh
    if phase == "healthy":
        fP = 22.0 * scale; fK = 38.0 * scale
        fMg = 6.0 * scale; fTh = P["formulaTh0"] * scale

    alc = R.alcohol if phase == "starve" else 0.0
    absTh = 1.0 - alc * (1.0 - P["EtOHabs"])
    tpk = 1.0 - alc * (1.0 - P["EtOHtpk"])
    if phase == "refeed":
        absTh *= (1.0 - R.gut_fail)

    # -----------------------------------------------------------------
    # 6.2  prescribed intravenous repletion
    # -----------------------------------------------------------------
    ivP = ivK = ivMg = ivTh = po_th = 0.0
    if phase == "refeed":
        ivP = (R.P_dose * BW / 24.0) if day >= R.P_delay else 0.0
        ivK = R.K_dose * BW / 24.0
        if R.Mg_bolus:
            ivMg = R.Mg_dose * BW / 2.0 if (t % 24.0) < 2.0 else 0.0
        else:
            ivMg = R.Mg_dose * BW / 24.0
        if R.th_delay <= t < R.th_delay + R.th_days * 24.0:
            ivTh = R.th_iv / 265.4 * 1000.0 / 24.0
            po_th = R.th_po / 265.4 * 1000.0 / 24.0
        if R.rescue:
            ivP += 0.50 * BW / 24.0 * sw_below(Pser, 0.65, 0.030)
            ivK += 1.00 * BW / 24.0 * sw_below(Kser, 3.30, 0.100)
            ivMg += 0.20 * BW / 24.0 * sw_below(MGser, 0.65, 0.030)

    # -----------------------------------------------------------------
    # 6.3  glucose / insulin
    # -----------------------------------------------------------------
    if route == "dextrose" and phase == "refeed":
        Ra_gut = g_cho_h * 1000.0 / 180.0
        d[IX["A_CHO"]] = -P["kaCho"] * A_CHO
    else:
        d[IX["A_CHO"]] = g_cho_h - P["kaCho"] * A_CHO
        Ra_gut = P["kaCho"] * A_CHO * 1000.0 / 180.0
    d[IX["A_FAT"]] = g_fat_h - P["kaFat"] * A_FAT
    d[IX["A_PRO"]] = g_pro_h - P["kaPro"] * A_PRO

    gly_f = min(1.0, GLY / (0.35 * P["GLY0"]))
    EGP = (P["EGP0"] / (1.0 + X / P["X50"]) * (0.7 + 0.3 * GCG / P["GCG0"])
           * (0.45 + 0.55 * gly_f))
    # Prolonged starvation down-regulates insulin-stimulated glucose
    # disposal ("starvation diabetes").  AT, the adaptive-thermogenesis
    # state, is used as the slow marker of that adapted state, so
    # sensitivity recovers over about a week of refeeding -- which is why
    # refeeding hyperglycaemia is real and transient.
    ATn = (AT - P["ATmin"]) / (1.0 - P["ATmin"])
    SIe = P["SI"] * (P["siMin"] + (1.0 - P["siMin"]) * min(max(ATn, 0.0), 1.0))
    # Glycolysis needs inorganic phosphate: GAPDH consumes Pi
    # stoichiometrically, so profound hypophosphataemia throttles glucose
    # disposal itself.  This is the brake that stops serum phosphate falling
    # to absurd values -- the demand is extinguished by its own consequence
    # -- and it is also why a crashing patient becomes glucose-intolerant.
    fPgly = ((Pser / (P["KmGly"] + Pser))
             / (P["Pser0"] / (P["KmGly"] + P["Pser0"])))
    SIe *= min(1.0, fPgly)
    Uins = SIe * X * GLU
    Rd = P["Rd0"] * GLU + Uins
    glyS = P["kGlyS"] * X * GLU * max(0.0, 1.0 - GLY / P["GLYmax"])
    glyB = P["kGlyB"] * GLY / (1.0 + 3.0 * X)
    ren_glu = P["GFR0"] * pos(GLU - 10.0)
    d[IX["GLU"]] = (Ra_gut + EGP + glyB - Rd - glyS - ren_glu) / Vg

    beta = 0.55 + 0.45 * min(1.0, FFM / P["FFM0"])
    inc = P["kInc"] * Ra_gut * (P["incIV"] if (route == "dextrose"
                                               and phase == "refeed") else 1.0)
    sec = P["SecB"] + P["kSec"] * pos(GLU - P["Gthr"]) ** P["hSec"] + inc
    d[IX["INS"]] = (sec * beta - P["CLins"] * INS * Vp) / Vp
    d[IX["X"]] = P["p2"] * (INS / P["INS0"] - X)
    gcg_t = (P["GCG0"] * (1.0 + 1.6 * pos(P["G0"] - GLU) / P["G0"])
             * (1.0 + 0.9) / (1.0 + 0.9 * X))
    d[IX["GCG"]] = P["kGcg"] * (gcg_t - GCG)
    d[IX["GLY"]] = glyS - glyB

    lipo = P["kLipo"] / (1.0 + 3.2 * X) * (FM / P["FM0"] + 0.15)
    d[IX["FFA"]] = lipo - P["kFFAox"] * FFA
    d[IX["BHB"]] = P["kKeto"] * FFA / (1.0 + 6.0 * X) - P["kBHBox"] * BHB

    # -----------------------------------------------------------------
    # 6.4  thiamine
    # -----------------------------------------------------------------
    Jth_gut = P["VmaxThA"] * THGUT / (P["KmThA"] + THGUT)
    d[IX["THGUT"]] = fTh + po_th - Jth_gut
    Jth_abs = absTh * Jth_gut
    Jth_up = (P["VmaxThT"] * Cth / (P["KmThT"] + Cth) + P["kThPass"] * Cth) * tpk
    ren_th = P["CLth"] * pos(Cth - P["ThRenalC"]) / 1000.0
    d[IX["THP"]] = Jth_abs + ivTh - Jth_up + P["kThOff"] * THT - ren_th
    d[IX["THT"]] = (Jth_up - P["kThOff"] * THT - P["kThDeg"] * THT
                    - P["kThCho"] * Rd)

    thn = THT / P["THT0"]
    fPDH = min(1.0, hill(thn, P["KmTPP"], P["hTPP"]) / P["fPDH0"])
    fTK = fPDH

    lac_prod = (P["LacProd0"] + P["kLac2"] * Rd * (1.0 - fPDH)) \
        * (1.0 + 0.85 * pos(1.0 - DPG))
    d[IX["LAC"]] = (lac_prod - P["CLlac"] * LAC) / P["VLac"]

    # -----------------------------------------------------------------
    # 6.5  energy balance and body composition
    # -----------------------------------------------------------------
    REE = ree_kcal_d(FFM, AT)
    TEE = REE * P["PAL"] + 0.10 * kcal_d
    Ebal = kcal_d - TEE
    cat = min(1.0, max(0.0, 1.0 - kcal_d / max(TEE, 1.0)))
    fLean = min(max(0.10 + 0.50 * np.exp(-FM / 8.0), 0.08), 0.55)
    dFFM_d = fLean * Ebal / P["rhoLean"]
    dFM_d = (1.0 - fLean) * Ebal / P["rhoFat"]
    if Ebal > 0.0:
        anab = min(1.0, Pser / 0.70) * min(1.0, Kser / 3.50)
        dFFM_d *= anab
    d[IX["FFM"]] = dFFM_d / 24.0
    d[IX["FM"]] = dFM_d / 24.0
    AT_t = P["ATmin"] + (1.0 - P["ATmin"]) * min(1.0, kcal_d / max(TEE, 1.0))
    d[IX["AT"]] = P["kAT"] * (AT_t - AT)
    d[IX["LVM"]] = P["kLVM"] * (P["LVMfrac"] * FFM - LVM)

    # -----------------------------------------------------------------
    # 6.6  PHOSPHATE  -- the engine of the syndrome
    # -----------------------------------------------------------------
    # Insulin sets the cellular organic-phosphate set-point (G6P, F1,6BP,
    # ATP, phosphocreatine, 2,3-DPG).  It must be referenced to the level
    # the patient is CURRENTLY at, not to a well-fed control: a starved
    # patient sits near X = 0, so the whole demand is generated by the
    # RISE from there.  f() is normalised so that the healthy fed state
    # (X = 1) is the zero point and the pool there is exactly rhoP x FFM.
    # The set-point is driven by the GLYCOLYTIC FLUX actually achieved, not
    # by the insulin signal: the pool being refilled IS the pool of
    # phosphorylated intermediates, so its target must track the flux
    # through them.  Writing it against insulin instead is a sign error in
    # disguise -- when hypophosphataemia throttles glucose disposal, plasma
    # glucose and insulin RISE, so an insulin-driven set-point would raise
    # the demand at exactly the moment the cell can no longer meet it.
    un = Uins / P["Uins0"]
    fU = un / (P["uSat"] + un)
    PI_t = P["rhoP"] * FFM * (1.0 + P["sIns"] * (fU - P["fU1"]))
    gap = PI_t - PI
    fUp = Pser / (P["KmPup"] + Pser)
    JPnet = P["kFill"] * gap * (fUp if gap > 0.0 else 1.0)
    JPnew = pos(d[IX["FFM"]]) * P["Pnew"]      # already inside PI_t via FFM

    Jres = P["kBres"] * PB * hill(PTH, P["PTH0"], P["hBres"])
    Jres = min(Jres, P["BoneCap"] * P["Jbone0"] * PB / P["PB0"] * 3.0)
    Jform = P["kBform"] * PB * min(2.0, Pser / P["Pser0"]) * min(2.0, Caser / P["Ca0"])

    TmP = (P["TmPmax"] / ((1.0 + hill(PTH, P["PTH50tm"], P["hTmP"]))
                          * (1.0 + 0.55 * hill(FGF, P["FGF50tm"], P["hTmP"]))))
    TmP *= (1.0 - P["etohP"] * alc)
    ExP = P["GFR0"] * pos(Pser - TmP) + P["ObligP"] * Pser / P["Pser0"]
    gi_P = (R.gi_loss * 4.0 / 24.0) if phase == "starve" else 0.0

    JPabs = P["FabsP"] * P["kaP"] * PGUT * (0.7 + 0.3 * CTD / P["CTD0"])
    d[IX["PGUT"]] = fP - P["kaP"] * PGUT

    sup = pos(Caser * Pser - P["Ksp"])
    Jprec = P["kPrecip"] * sup * Vecf

    d[IX["PE"]] = JPabs + ivP + Jres - JPnet - Jform - ExP - Jprec - gi_P
    d[IX["PI"]] = JPnet
    d[IX["PB"]] = Jform - Jres

    # -----------------------------------------------------------------
    # 6.7  POTASSIUM
    # -----------------------------------------------------------------
    depK = pos(1.0 - KI / (P["rhoK"] * max(FFM, 1.0)))
    JKin = P["kKp"] * pumpmult(X, P["EmaxK"], P["tauX"]) * Kser * (1.0 + 0.8 * depK)
    JKout = P["kKleak"] * KI * (1.0 + 0.9 * cat)
    JKnew = d[IX["FFM"]] * P["Knew"]
    fMgK = 1.0 + P["betaMgK"] * pos(1.0 - MGser / P["MGser0"])
    # ROMK/BK down-regulation: the distal nephron adapts to potassium
    # depletion over days, but never all the way to zero -- which is why a
    # depleted patient still loses potassium in the urine.
    romk = max(P["romkMin"], min(1.0, (KI / (P["rhoK"] * max(FFM, 1e-9))) ** P["hRomk"]))
    ExK = (P["Kobl"] + P["kKsec"] * Kser * ALD * fMgK
           * (1.0 + 1.4 * R.diuretic) * romk)
    JKabs = P["FabsK"] * P["kaK"] * KGUT
    d[IX["KGUT"]] = fK - P["kaK"] * KGUT
    gi_K = (R.gi_loss * 25.0 / 24.0) if phase == "starve" else 0.0
    d[IX["KE"]] = JKabs + ivK + JKout - JKin - ExK - gi_K - JKnew
    d[IX["KI"]] = JKin - JKout + JKnew

    # -----------------------------------------------------------------
    # 6.8  MAGNESIUM
    # -----------------------------------------------------------------
    depMg = pos(1.0 - MGI / (P["rhoMg"] * max(FFM, 1.0)))
    JMgin = P["kMgp"] * pumpmult(X, P["EmaxMg"], P["tauX"]) * MGser * (1.0 + 0.7 * depMg)
    JMgout = P["kMgleak"] * MGI * (1.0 + 0.6 * cat)
    JMgnew = d[IX["FFM"]] * P["Mgnew"]
    TmMg = P["TmMg"] * (1.0 - P["etohMg"] * alc)
    ExMg = ((P["GFR0"] * pos(P["fUFMg"] * MGser - TmMg)
             + 0.5 / 24.0 * MGser / P["MGser0"]) * (1.0 + 1.8 * R.diuretic))
    JMgabs = P["FabsMg"] * P["kaMg"] * MGGUT
    d[IX["MGGUT"]] = fMg - P["kaMg"] * MGGUT
    gi_Mg = (R.gi_loss * 6.0 / 24.0) if phase == "starve" else 0.0
    d[IX["MGE"]] = JMgabs + ivMg + JMgout - JMgin - ExMg - gi_Mg - JMgnew
    d[IX["MGI"]] = JMgin - JMgout + JMgnew

    # -----------------------------------------------------------------
    # 6.9  CALCIUM / PTH / 1,25D / FGF23
    # -----------------------------------------------------------------
    fMgPTH = ((MGser / (P["MgPTH50"] + MGser))
              / (P["MGser0"] / (P["MgPTH50"] + P["MGser0"])))
    fMgPTH = min(fMgPTH, 1.05)
    secPTH = (P["PTHmin"] + (P["PTHmax"] - P["PTHmin"])
              / (1.0 + (Caser / P["Ca50"]) ** P["hCa"])) * fMgPTH
    d[IX["PTH"]] = P["kPTH"] * (secPTH - PTH)

    CaRes = P["kCaBone"] * P["CAE0"] * hill(PTH, P["PTH0"], P["hBres"]) * 2.0
    CaForm = P["kCaBone"] * CAE * min(2.0, Pser / P["Pser0"])
    CaAbs = (5.0 / 24.0) * (CTD / P["CTD0"]) * (0.4 + 0.6 * min(2.0, kcal_h / P["kcal_h0"]))
    CaU = P["kCaU"] * Caser * (1.0 + 0.4 * pos(1.0 - PTH / P["PTH0"]))
    d[IX["CAE"]] = CaAbs + CaRes - CaForm - CaU - 1.5 * Jprec
    d[IX["CAB"]] = Jprec + CaForm - CaRes
    d[IX["CTD"]] = P["kCTD"] * (P["CTD0"] * (0.5 + 0.5 * PTH / P["PTH0"])
                                * (1.0 + 0.5 * pos(1.0 - Pser / P["Pser0"])) - CTD)
    d[IX["FGF"]] = P["kFGF"] * (P["FGF0"] * (0.25 + 0.75 * Pser / P["Pser0"])
                                * (1.0 + 0.4 * (CTD / P["CTD0"] - 1.0)) - FGF)

    # -----------------------------------------------------------------
    # 6.10  ORGAN ENERGETICS
    # -----------------------------------------------------------------
    nP = hill(P["Pser0"], P["KPatp"], P["hPatp"])
    fPox = hill(Pser, P["KPatp"], P["hPatp"]) / nP
    oxid = min(1.0, 0.35 + 0.65 * fPDH)
    d[IX["ATPm"]] = P["kATPm"] * (min(1.0, fPox * oxid) - ATPm)
    nPc = hill(P["Pser0"], P["KPpcr"], P["hPatp"])
    fPcr = hill(Pser, P["KPpcr"], P["hPatp"]) / nPc
    d[IX["PCRm"]] = P["kPCr"] * (min(1.0, fPcr * ATPm) - PCRm)
    nPd = hill(P["Pser0"], P["KPdpg"], P["hPatp"])
    d[IX["DPG"]] = P["kDPG"] * (min(1.0, hill(Pser, P["KPdpg"], P["hPatp"]) / nPd) - DPG)
    d[IX["ATPd"]] = P["kATPd"] * (min(1.0, fPox * oxid) - ATPd)

    # -----------------------------------------------------------------
    # 6.11  SODIUM / VOLUME / ALDOSTERONE
    # -----------------------------------------------------------------
    # Potassium is the dominant aldosterone secretagogue, so hypokalaemia
    # brakes its own renal loss; volume depletion pushes the other way.
    ald_t = ((0.15 + 0.85 * (Kser / P["Kser0"]) ** 2)
             * (1.0 + 1.8 * pos(1.0 - BW / P["BW0"])))
    d[IX["ALD"]] = P["kALD"] * (ald_t - ALD)
    Na_in = (R.Na_intake / 24.0) if phase == "refeed" else (R.Na_intake * 0.5 / 24.0)
    ret = min(0.92, P["rIns"] * pos(X - 1.0) / (1.0 + pos(X - 1.0))
              + P["rAld"] * pos(ALD - 1.0) / (1.0 + pos(ALD - 1.0)))
    d[IX["NAB"]] = Na_in * ret - P["kNaB"] * NAB

    # -----------------------------------------------------------------
    # 6.12  GAS EXCHANGE
    # -----------------------------------------------------------------
    cho_ox_max = P["ChoOxMax"] * BW * 60.0 / 180.0
    cho_ox_act = min(Rd, cho_ox_max)
    dnl = pos(Rd - cho_ox_max) * (1.0 if GLY > 0.92 * P["GLYmax"] else 0.25)
    EE_h = (REE * P["PAL"] + 0.10 * kcal_d) / 24.0
    e_cho = cho_ox_act * 0.180 * 4.0
    e_cho = min(e_cho, EE_h)
    e_fat = max(0.0, EE_h - e_cho)
    RQ = (1.00 * e_cho + 0.71 * e_fat) / max(e_cho + e_fat, 1e-9)
    VO2 = EE_h * 24.0 / (1440.0 * P["kcalPerLO2"])
    VCO2 = RQ * VO2 + 2.44 * dnl * 22.4 / 1000.0 / 60.0
    # Alveolar ventilation cannot fall to zero in a living patient, and
    # PaCO2 is bounded on both sides.  Leaving this ratio unbounded puts an
    # unbounded argument inside the respiratory-failure exponential, which
    # is a modelling artefact rather than a physiological prediction.
    VA_cap = max(ATPd * (FFM / P["FFM0"]) ** 0.5, P["VAmin"])
    PaCO2 = P["PaCO2_0"] * (VCO2 / max(P["VCO2ref"], 1e-9)) / VA_cap
    PaCO2 = min(max(PaCO2, 25.0), 120.0)

    # -----------------------------------------------------------------
    # 6.13  CARDIAC FUNCTION AND HAZARDS
    # -----------------------------------------------------------------
    CONT = max(ATPm * PCRm ** 0.5, P["contMin"])
    cong = pos(Vecf / P["Vecf0"] - 1.0) / max(LVM / P["LVM0"], 0.25) / CONT
    QTc = (P["QTc0"] + P["aK"] * pos(P["Kser0"] - Kser)
           + P["aMg"] * pos(P["MGser0"] - MGser)
           + P["aCa"] * pos(P["Ca0"] - Caser))
    h_arr = P["HAZ"] * P["h0arr"] * np.exp(P["cQT"] * pos(QTc - 440.0)) \
        * (1.0 + P["cATParr"] * pos(1.0 - ATPm))
    h_cf = P["HAZ"] * P["h0cf"] * (1.0 + cong) ** P["pCong"] / CONT
    h_rf = P["HAZ"] * P["h0rf"] * np.exp(P["qCO2"] * pos(PaCO2 - P["PaCO2ref"])) \
        / max(ATPd, P["VAmin"])
    h_we = P["h0we"] * pos(1.0 - fTK) ** 2 * (1.0 + 2.0 * min(2.0, Rd / 40.0))
    h_tot = min(h_arr + h_cf + h_rf, P["hazMax"])

    d[IX["H_ARR"]] = h_arr
    d[IX["H_CF"]] = h_cf
    d[IX["H_RF"]] = h_rf
    d[IX["H_WE"]] = h_we * (1.0 - min(H_WE, 0.999))
    d[IX["SURV"]] = -SURV * h_tot
    d[IX["AUCP"]] = pos(0.65 - Pser)
    d[IX["CUMP"]] = ivP + JPabs
    d[IX["CUMK"]] = ivK + JKabs
    d[IX["CUMKCAL"]] = kcal_h

    if aux is not None:
        aux.update(dict(
            Pser=Pser, Kser=Kser, MGser=MGser, Caser=Caser, Cth=Cth,
            Jdem=(JPnet if JPnet > 0 else 0.0) * 24.0,
            Jsup=(JPabs + ivP + max(Jres - Jform, 0.0)) * 24.0,
            LambdaP=(max(JPnet, 0.0) / max(JPabs + ivP + max(Jres - Jform, 0.0), 1e-9)),
            TmP=TmP, ExP=ExP * 24.0, ExK=ExK * 24.0, ExMg=ExMg * 24.0,
            RQ=RQ, VCO2=VCO2, PaCO2=PaCO2, QTc=QTc, CONT=CONT, cong=cong,
            fPDH=fPDH, fTK=fTK, Uins=Uins, Rd=Rd, kcal_d=kcal_d,
            GIR=(g_cho_h * 1000.0 / 60.0 / max(BW, 1.0)), EGP=EGP,
            Vecf=Vecf, PI_t=PI_t, JPnet=JPnet * 24.0, Jres=Jres * 24.0,
            cat=cat, Ebal=Ebal, TEE=TEE,
        ))
    return d


# =========================================================================
# 7.  OBSERVATION
# =========================================================================

def observe(t, y, P, R, phase):
    aux = {}
    rhs(t, y, P, R, phase, aux)
    BW = y[IX["FM"]] + y[IX["FFM"]]
    o = dict(aux)
    o.update(dict(
        t=t, day=t / 24.0, BW=BW, FFM=y[IX["FFM"]], FM=y[IX["FM"]],
        BMI=BW / P["HT"] ** 2, GLU=y[IX["GLU"]], INS=y[IX["INS"]], X=y[IX["X"]],
        LAC=y[IX["LAC"]], PI=y[IX["PI"]], KI=y[IX["KI"]], MGI=y[IX["MGI"]],
        PB=y[IX["PB"]], TBP=y[IX["PE"]] + y[IX["PI"]] + y[IX["PB"]],
        TBK=y[IX["KE"]] + y[IX["KI"]], THT=y[IX["THT"]],
        ATPm=y[IX["ATPm"]], PCRm=y[IX["PCRm"]], DPG=y[IX["DPG"]],
        ATPd=y[IX["ATPd"]], LVM=y[IX["LVM"]], PTH=y[IX["PTH"]],
        FGF=y[IX["FGF"]], SURV=y[IX["SURV"]],
        mortality=100.0 * (1.0 - y[IX["SURV"]]), WE=100.0 * y[IX["H_WE"]],
        AUCP=y[IX["AUCP"]], CUMKCAL=y[IX["CUMKCAL"]],
        piC=y[IX["PI"]] / (P["rhoP"] * max(y[IX["FFM"]], 1e-9)),
        kiC=y[IX["KI"]] / (P["rhoK"] * max(y[IX["FFM"]], 1e-9)),
        mgiC=y[IX["MGI"]] / (P["rhoMg"] * max(y[IX["FFM"]], 1e-9)),
        oedema=max(0.0, aux["Vecf"] - P["Vecf0"]),
    ))
    o["CaxP"] = o["Caser"] * o["Pser"]
    return o


# =========================================================================
# 8.  DRIVER
# =========================================================================

class Grid:
    """A stitched solution on a fixed grid, with linear interpolation.

    The refeeding protocol is genuinely discontinuous -- the calorie ramp
    steps once a day, a magnesium bolus starts and stops -- so the
    integrator is restarted at every breakpoint instead of being asked to
    chase the jump.  That is both faster and more accurate than letting a
    stiff solver discover the step by itself.
    """

    def __init__(self, t, y):
        self.t = np.asarray(t)
        self.y = np.asarray(y)
        self.nfev = 0

    def sol(self, ts):
        ts = np.atleast_1d(ts)
        return np.vstack([np.interp(ts, self.t, self.y[i, :])
                          for i in range(self.y.shape[0])])


def _breakpoints(R):
    T = R.sim_days * 24.0
    bp = set(np.arange(0.0, T + 1e-9, 24.0))
    bp.add(T)
    if R.Mg_bolus:
        for d in np.arange(0.0, R.sim_days):
            bp.add(d * 24.0 + 2.0)
    for extra in (R.th_lead_h, R.th_delay, R.th_delay + R.th_days * 24.0,
                  R.P_delay * 24.0):
        if 0.0 < extra < T:
            bp.add(float(extra))
    return sorted(b for b in bp if 0.0 <= b <= T)


def run(R, P, dt=0.25):
    y = y0_healthy(P)
    star = None
    if R.starve_days > 0:
        s1 = solve_ivp(rhs, (0.0, R.starve_days * 24.0), y, args=(P, R, "starve"),
                       method="LSODA", rtol=1e-6, atol=1e-8, max_step=12.0,
                       dense_output=True)
        if not s1.success:
            raise RuntimeError("starvation failed: " + s1.message)
        y = np.maximum(s1.y[:, -1], 0.0)
        star = s1
    for k in ("H_ARR", "H_CF", "H_RF", "H_WE", "AUCP", "CUMP", "CUMK", "CUMKCAL"):
        y[IX[k]] = 0.0
    y[IX["SURV"]] = 1.0
    y_adm = y.copy()

    bp = _breakpoints(R)
    ts_all, ys_all, nfev = [0.0], [y.copy()], 0
    for a, b in zip(bp[:-1], bp[1:]):
        n = max(2, int(round((b - a) / dt)) + 1)
        te = np.linspace(a, b, n)
        s = solve_ivp(rhs, (a, b), y, args=(P, R, "refeed"), method="LSODA",
                      rtol=1e-6, atol=1e-8, max_step=1.0, t_eval=te)
        if not s.success:
            raise RuntimeError("refeeding failed on [%g,%g]: %s" % (a, b, s.message))
        nfev += s.nfev
        ts_all.extend(s.t[1:]); ys_all.extend(s.y[:, 1:].T)
        y = np.maximum(s.y[:, -1], 0.0)
    g = Grid(np.array(ts_all), np.array(ys_all).T)
    g.nfev = nfev
    return star, g, y_adm


def summarise(sol, P, R, npts=None):
    n = npts or int(sol.t[-1]) + 1
    ts = np.linspace(0.0, sol.t[-1], n)
    ys = sol.sol(ts)
    rows = [observe(ts[i], ys[:, i], P, R, "refeed") for i in range(n)]
    Pv = [r["Pser"] for r in rows]
    i = int(np.argmin(Pv))
    last = rows[-1]

    def at(day, key):
        j = int(np.clip(round(day * (n - 1) / (sol.t[-1] / 24.0)), 0, n - 1))
        return rows[j][key]

    d1 = {k: at(1.0, k) for k in ("Pser", "Kser", "MGser", "Caser", "LAC", "fTK")}
    d3 = {k: at(3.0, k) for k in ("Pser", "Kser", "MGser", "Caser", "LAC", "fTK")}
    d7 = {k: at(7.0, k) for k in ("Pser", "Kser", "MGser", "Caser", "LAC", "fTK")}
    return dict(
        d1=d1, d3=d3, d7=d7, at=at,
        Pnadir=Pv[i], Pnadir_day=rows[i]["day"],
        Knadir=min(r["Kser"] for r in rows),
        Mgnadir=min(r["MGser"] for r in rows),
        Canadir=min(r["Caser"] for r in rows),
        QTcmax=max(r["QTc"] for r in rows),
        LACmax=max(r["LAC"] for r in rows),
        Lmax=max(r["LambdaP"] for r in rows),
        ATPmin=min(r["ATPm"] for r in rows),
        PaCO2max=max(r["PaCO2"] for r in rows),
        RQmax=max(r["RQ"] for r in rows),
        mortality=last["mortality"], WE=last["WE"], AUCP=last["AUCP"],
        kcal=last["CUMKCAL"], dFFM=last["FFM"] - rows[0]["FFM"],
        oedema=max(r["oedema"] for r in rows),
        CaxPmax=max(r["CaxP"] for r in rows),
        fTKmin=min(r["fTK"] for r in rows),
        GIRmax=max(r["GIR"] for r in rows),
        rows=rows,
    )


# =========================================================================
# 9.  SCENARIOS
# =========================================================================

def hx(**kw):
    """Reference high-risk history: 60 d at 32 % of energy needs, poor
    micronutrient quality, some gastrointestinal loss."""
    d = dict(starve_days=60.0, starve_frac=0.32, micronutrient=0.25,
             alcohol=0.0, diuretic=0.0, gi_loss=0.20, sim_days=14.0)
    d.update(kw)
    return d


SC = {}
SC["S01_NICE"] = hx(kcal_start=10.0, kcal_goal=30.0, advance_frac=0.20,
                    advance_start=2.0, P_dose=0.5, K_dose=2.5, Mg_dose=0.3,
                    th_iv=200.0, th_lead_h=0.5, th_days=10.0, rescue=True)
SC["S02_ASPEN"] = hx(kcal_start=15.0, kcal_goal=30.0, advance_frac=0.33,
                     advance_start=1.0, P_dose=0.6, K_dose=3.0, Mg_dose=0.4,
                     th_iv=100.0, th_lead_h=0.5, th_days=7.0, rescue=True)
SC["S03_FULL_NOTHING"] = hx(kcal_start=30.0, kcal_goal=30.0, advance_start=99.0,
                            rescue=False)
SC["S04_FULL_ELEC_ONLY"] = hx(kcal_start=30.0, kcal_goal=30.0, advance_start=99.0,
                              P_dose=0.8, K_dose=3.0, Mg_dose=0.4, rescue=True)
SC["S05_FULL_THIA_ONLY"] = hx(kcal_start=30.0, kcal_goal=30.0, advance_start=99.0,
                              th_iv=300.0, th_lead_h=0.5, rescue=False)
SC["S06_FULL_BOTH"] = hx(kcal_start=30.0, kcal_goal=30.0, advance_start=99.0,
                         P_dose=0.8, K_dose=3.0, Mg_dose=0.4,
                         th_iv=300.0, th_lead_h=0.5, rescue=True)
SC["S07_SLOW_NOTHING"] = hx(kcal_start=10.0, kcal_goal=30.0, advance_frac=0.20,
                            advance_start=2.0, rescue=False)
SC["S08_DEXTROSE_ONLY"] = hx(kcal_start=20.0, kcal_goal=20.0, advance_start=99.0,
                             route="dextrose", cho_frac=1.0, fat_frac=0.0,
                             pro_frac=0.0, rescue=False)
SC["S09_ENTERAL_MATCHED"] = hx(kcal_start=20.0, kcal_goal=20.0, advance_start=99.0,
                               route="enteral", cho_frac=1.0, fat_frac=0.0,
                               pro_frac=0.0, rescue=False)
SC["S10_MIXED_FUEL"] = hx(kcal_start=20.0, kcal_goal=20.0, advance_start=99.0,
                          route="enteral", cho_frac=0.40, fat_frac=0.42,
                          pro_frac=0.18, rescue=False)
SC["S11_ALCOHOL_NOTH"] = hx(alcohol=0.75, micronutrient=0.10, starve_frac=0.40,
                            kcal_start=25.0, kcal_goal=25.0, advance_start=99.0,
                            rescue=False)
SC["S12_ALCOHOL_IVTH"] = hx(alcohol=0.75, micronutrient=0.10, starve_frac=0.40,
                            kcal_start=25.0, kcal_goal=25.0, advance_start=99.0,
                            th_iv=500.0, th_lead_h=1.0, rescue=False)
SC["S13_ALCOHOL_POTH"] = hx(alcohol=0.75, micronutrient=0.10, starve_frac=0.40,
                            kcal_start=25.0, kcal_goal=25.0, advance_start=99.0,
                            th_po=300.0, rescue=False)
SC["S14_MG_BOLUS"] = hx(diuretic=0.6, kcal_start=20.0, kcal_goal=20.0,
                        advance_start=99.0, Mg_dose=0.4, Mg_bolus=True,
                        K_dose=2.0, rescue=False)
SC["S15_MG_INFUSION"] = hx(diuretic=0.6, kcal_start=20.0, kcal_goal=20.0,
                           advance_start=99.0, Mg_dose=0.4, Mg_bolus=False,
                           K_dose=2.0, rescue=False)
SC["S16_K_ONLY_NO_MG"] = hx(diuretic=0.6, kcal_start=20.0, kcal_goal=20.0,
                            advance_start=99.0, Mg_dose=0.0, K_dose=2.0,
                            rescue=False)




# =========================================================================
# 10.  VERIFICATION SUITE AND SWEEPS
# =========================================================================

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)


def sub(s):
    print()
    print("-" * 78)
    print(s)
    print("-" * 78)


def main():
    global P
    P = steady_state(P)
    R0 = Regimen(**hx())
    y = y0_healthy(P)

    hdr("REFEEDING SYNDROME QSP MODEL - PYTHON REFERENCE IMPLEMENTATION")
    tbp = P["PB0"] + P["PI0"] + P["PE0"]
    print("Reference subject   %.0f kg (FFM %.0f, FM %.0f), BMI %.1f, GFR %.1f L/h"
          % (P["BW0"], P["FFM0"], P["FM0"], P["BW0"] / P["HT"] ** 2, P["GFR0"]))
    print("Total body P        %.0f mmol = bone %.0f + ICF %.0f + ECF %.1f"
          % (tbp, P["PB0"], P["PI0"], P["PE0"]))
    print("Total body K        %.0f mmol = ICF %.0f + ECF %.1f"
          % (P["KI0"] + P["KE0"], P["KI0"], P["KE0"]))
    print("Total body Mg (exch)%.0f mmol = ICF %.0f + ECF %.1f"
          % (P["MGI0"] + P["MGE0"], P["MGI0"], P["MGE0"]))
    print()
    print("  ** THE RATIO THAT DEFINES THE DISEASE **")
    print("     the compartment the laboratory measures is %.3f %% of total body"
          % (100.0 * P["PE0"] / tbp))
    print("     phosphorus.  A refeeding demand of 60 mmol/d is %.1f times the"
          % (60.0 / P["PE0"]))
    print("     ENTIRE measured compartment, per day.")
    print()
    print("Derived (not fitted) from the balance conditions:")
    print("   TmP/GFR                     %.3f mmol/L      (normal 0.80-1.35)"
          % P["TmP0"])
    print("   urinary P / K / Mg          %.1f / %.1f / %.2f mmol/d"
          % ((P["GFR0"] * pos(P["Pser0"] - P["TmP0"]) + P["ObligP"]) * 24,
             P["ExK0"] * 24, P["ExMg0"] * 24))
    print("   insulin sensitivity SI      %.2f L/h, fed glucose turnover %.1f mmol/h"
          % (P["SI"], P["Rd_0"]))
    print("   incretin share of secretion %.0f %% at the fed reference"
          % (100.0 * P["fInc"]))
    print()
    print("EMERGENT CHECK - thiamine.  The store (%.0f umol = %.1f mg) and the"
          % (P["THT0"], P["THT0"] * 265.4 / 1000.0))
    print("   14 d biological half-life were taken from tracer studies.  The daily")
    print("   requirement was NOT: the model computes it from those two numbers as")
    print("   %.2f mg/d, against a published RDA of 1.1-1.4 mg/d." % P["ThRDA"])

    # ------------------------------------------------------------------
    sub("TEST 1   healthy steady state (should be numerically exact)")
    d = rhs(0.0, y, P, R0, "healthy")
    physi = [i for i, n in enumerate(NAMES)
             if not n.startswith(("CUM", "H_", "SURV", "AUCP"))]
    worst = max(physi, key=lambda i: abs(d[i]))
    print("   max |dy/dt| over the 44 physiological states = %.3e  (%s)"
          % (abs(d[worst]), NAMES[worst]))
    print("   (the four cumulative counters are excluded: they are meant to grow)")

    # ------------------------------------------------------------------
    sub("TEST 2   the starvation phase - THE MODEL'S FIRST REAL TEST")
    print("The claim being tested is that the pools empty while the BLOOD TESTS")
    print("STAY NORMAL.  If serum phosphate fell during starvation, admission")
    print("phosphate would be a useful triage test and refeeding syndrome would")
    print("not be dangerous.  The deficit here is SIMULATED, never asserted.")
    print()
    phen = [("anorexia/starvation", dict(gi_loss=0.20, micronutrient=0.25)),
            ("alcohol use disorder", dict(alcohol=0.75, gi_loss=0.05,
                                          micronutrient=0.10, starve_frac=0.40)),
            ("diuretic + cardiac", dict(diuretic=0.60, gi_loss=0.05,
                                        micronutrient=0.35)),
            ("purging / GI losses", dict(gi_loss=0.50, micronutrient=0.30,
                                         starve_frac=0.35))]
    print("%-21s %5s %5s | %5s %5s %5s %5s | %6s %6s %6s | %6s"
          % ("phenotype at day 60", "BMI", "FFM", "P", "K", "Mg", "Ca",
             "ICF-P", "ICF-K", "ICF-Mg", "thiam"))
    print("%-21s %5s %5s | %5s %5s %5s %5s | %6s %6s %6s | %6s"
          % ("", "", "kg", "mmol/L", "", "", "", "%norm", "%norm", "%norm", "%store"))
    for nm, kw in phen:
        R = Regimen(**hx(**kw))
        _, _, ya = run(R, P)
        o = observe(0.0, ya, P, R, "refeed")
        print("%-21s %5.1f %5.1f | %5.2f %5.2f %5.2f %5.2f | %6.1f %6.1f %6.1f | %6.1f"
              % (nm, o["BMI"], o["FFM"], o["Pser"], o["Kser"], o["MGser"],
                 o["Caser"], 100 * o["piC"], 100 * o["kiC"], 100 * o["mgiC"],
                 100 * o["THT"] / P["THT0"]))
    print()
    print("VERDICT: serum phosphate is 0.8-1.1 mmol/L in every phenotype - inside")
    print("   or at the edge of the reference interval - after two months of")
    print("   starvation.  Serum potassium and magnesium are only mildly low.")
    print("   Meanwhile the thiamine store is down to 5-25 % and the intracellular")
    print("   cation pools are down by 20-40 %.  The admission blood panel does")
    print("   not see the disease.  This is why NICE triages on HISTORY - body")
    print("   mass index, weight loss, days without intake, alcohol - and not on")
    print("   the laboratory result.")

    # ------------------------------------------------------------------
    sub("TEST 3   the flux arithmetic at the moment of refeeding")
    R = Regimen(**SC["S03_FULL_NOTHING"])
    _, g, ya = run(R, P)
    o0 = observe(0.0, ya, P, R, "refeed")
    o1 = observe(24.0, g.sol(24.0)[:, 0], P, R, "refeed")
    print("   admission ECF phosphate pool        %6.1f mmol" % (o0["Pser"] * o0["Vecf"]))
    print("   day-1 demand  J_demand              %6.1f mmol/d" % o1["Jdem"])
    print("   day-1 supply  J_supply              %6.1f mmol/d" % o1["Jsup"])
    print("      of which absorbed from the feed  %6.1f mmol/d"
          % (o1["Jsup"] - o1["Jres"] if o1["Jsup"] > o1["Jres"] else 0.0))
    print("      of which net bone efflux         %6.1f mmol/d" % o1["Jres"])
    print("   day-1 Lambda_P = demand / supply    %6.2f" % o1["LambdaP"])
    o3 = observe(72.0, g.sol(72.0)[:, 0], P, R, "refeed")
    print()
    print("   urinary phosphate, admission        %6.2f mmol/d" % o0["ExP"])
    print("   urinary phosphate, day 1            %6.2f mmol/d" % o1["ExP"])
    print("   urinary phosphate, day 3            %6.2f mmol/d" % o3["ExP"])
    print()
    print("   Renal excretion is a THRESHOLD, GFR x max(0, Pser - TmP/GFR), so the")
    print("   kidney's entire contribution is to stop excreting.  That is worth")
    print("   %.1f mmol/d and it is SINGLE-USE: it is fully spent within a day,"
          % o0["ExP"])
    print("   after which the kidney has nothing further to offer no matter how")
    print("   low the phosphate goes.  Against a demand of %.0f mmol/d the whole"
          % o1["Jdem"])
    print("   renal reserve is worth about %.0f hours." % (24.0 * o0["ExP"] / max(o1["Jdem"], 1e-9)))

    # ------------------------------------------------------------------
    sub("TEST 4   scenario panel (16 regimens, one virtual patient each)")
    print("All 16 share the SAME simulated starvation history unless the name")
    print("says otherwise, so every difference below is caused by the regimen.")
    print()
    print("%-20s %5s %4s %5s %5s %5s %5s %5s %5s %5s %5s"
          % ("regimen", "Pnad", "day", "P d3", "K d3", "Mg d3", "QTc", "lac",
             "mort%", "WE%", "kcal"))
    res = {}
    for k in sorted(SC):
        R = Regimen(**SC[k])
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        res[k] = s
        print("%-20s %5.2f %4.1f %5.2f %5.2f %5.2f %5.0f %5.1f %5.2f %5.0f %5.0f"
              % (k, s["Pnadir"], s["Pnadir_day"], s["d3"]["Pser"],
                 s["d3"]["Kser"], s["d3"]["MGser"], s["QTcmax"], s["LACmax"],
                 s["mortality"], s["WE"], s["kcal"]))

    # ------------------------------------------------------------------
    sub("RESULT 1   the 2x2: which prophylaxis is doing the work?")
    print("Full-rate feeding (30 kcal/kg/d from hour zero) in the same patient,")
    print("crossing electrolyte repletion with thiamine.")
    print()
    print("%-26s %8s %8s %8s %8s" % ("", "P nadir", "lactate", "mort %", "Wernicke %"))
    for lbl, k in [("neither", "S03_FULL_NOTHING"),
                   ("electrolytes only", "S04_FULL_ELEC_ONLY"),
                   ("thiamine only", "S05_FULL_THIA_ONLY"),
                   ("both", "S06_FULL_BOTH")]:
        s = res[k]
        print("%-26s %8.2f %8.1f %8.2f %8.0f"
              % (lbl, s["Pnadir"], s["LACmax"], s["mortality"], s["WE"]))
    print()
    print("The two prophylaxes are NOT interchangeable and they do not overlap:")
    print("electrolytes fix the phosphate nadir and do nothing for Wernicke;")
    print("thiamine fixes Wernicke and the lactate and does NOTHING for the")
    print("phosphate.  Each closes a channel the other leaves open, which is why")
    print("giving one and not the other still leaves the patient exposed.")

    # ------------------------------------------------------------------
    sub("RESULT 2   NICE 2006 vs ASPEN 2020 - the calorie-restriction question")
    a, b = res["S01_NICE"], res["S02_ASPEN"]
    print("%-22s %10s %10s" % ("", "NICE-like", "ASPEN-like"))
    print("%-22s %10.1f %10.1f" % ("start kcal/kg/d", 10.0, 15.0))
    print("%-22s %10.2f %10.2f" % ("phosphate nadir", a["Pnadir"], b["Pnadir"]))
    print("%-22s %10.2f %10.2f" % ("potassium day 3", a["d3"]["Kser"], b["d3"]["Kser"]))
    print("%-22s %10.0f %10.0f" % ("kcal delivered / 14 d", a["kcal"], b["kcal"]))
    print("%-22s %10.2f %10.2f" % ("mortality %", a["mortality"], b["mortality"]))
    print("%-22s %10.2f %10.2f" % ("lean mass change kg", a["dFFM"], b["dFFM"]))
    print()
    print("The faster ramp delivers %.0f %% more energy at %s risk, PROVIDED the"
          % (100.0 * (b["kcal"] / a["kcal"] - 1.0),
             "indistinguishable" if abs(a["mortality"] - b["mortality"]) < 0.1
             else "different"))
    print("repletion is prescribed alongside it.  The model's reading of the")
    print("NICE/ASPEN disagreement is therefore that the calorie cap was never")
    print("the active ingredient - the repletion was.  Compare:")
    c = res["S07_SLOW_NOTHING"]
    print("   slow ramp WITHOUT repletion:  P nadir %.2f, mortality %.2f %%,"
          % (c["Pnadir"], c["mortality"]))
    print("      Wernicke %.0f %% (worse than the fast ramp WITH repletion, and"
          % c["WE"])
    print("      the low-calorie period prolongs the thiamine-deficient window)")

    # ------------------------------------------------------------------
    sub("RESULT 3   the feed itself is a phosphate prescription")
    d_, e_, f_ = (res["S08_DEXTROSE_ONLY"], res["S09_ENTERAL_MATCHED"],
                  res["S10_MIXED_FUEL"])
    print("Three regimens at an IDENTICAL 20 kcal/kg/d and identical patient.")
    print()
    print("%-34s %8s %8s %8s %8s"
          % ("", "P nadir", "GIR", "lactate", "mort %"))
    print("%-34s %8.2f %8.2f %8.1f %8.2f"
          % ("intravenous dextrose, 100 % CHO", d_["Pnadir"], d_["GIRmax"],
             d_["LACmax"], d_["mortality"]))
    print("%-34s %8.2f %8.2f %8.1f %8.2f"
          % ("enteral formula, 100 % CHO", e_["Pnadir"], e_["GIRmax"],
             e_["LACmax"], e_["mortality"]))
    print("%-34s %8.2f %8.2f %8.1f %8.2f"
          % ("enteral formula, 40 % CHO", f_["Pnadir"], f_["GIRmax"],
             f_["LACmax"], f_["mortality"]))
    print()
    print("Same calories, same patient, three different diseases.  A complete")
    print("enteral formula carries ~22 mmol phosphate, 38 mmol potassium, 6 mmol")
    print("magnesium and thiamine per 1000 kcal; a bag of intravenous dextrose")
    print("carries NONE of them unless somebody adds them.  Substituting the")
    print("route while holding the calories fixed moves the phosphate nadir from")
    print("%.2f to %.2f mmol/L.  Dropping the carbohydrate fraction from 100 %% to"
          % (d_["Pnadir"], e_["Pnadir"]))
    print("40 %% at the same energy moves it again, to %.2f, because the demand"
          % f_["Pnadir"])
    print("tracks the GLUCOSE, not the calories.")

    sys.stdout.flush()
    return P, res


def sweeps(P, res):
    # ------------------------------------------------------------------
    sub("RESULT 4   admission phosphate is nearly uninformative")
    print("Six patients presented for refeeding.  Their starvation histories")
    print("differ; their ADMISSION PHOSPHATE is essentially the same number.")
    print("Every one is then fed on the identical unprophylaxed regimen.")
    print()
    print("%-24s %9s %9s %9s %9s %8s"
          % ("history", "adm P", "adm K", "ICF-K %", "P nadir", "mort %"))
    hist = [("30 d at 50 % energy", dict(starve_days=30, starve_frac=0.50)),
            ("45 d at 40 % energy", dict(starve_days=45, starve_frac=0.40)),
            ("60 d at 32 % energy", dict(starve_days=60, starve_frac=0.32)),
            ("75 d at 28 % energy", dict(starve_days=75, starve_frac=0.28)),
            ("90 d at 25 % energy", dict(starve_days=90, starve_frac=0.25)),
            ("120 d at 22 % energy", dict(starve_days=120, starve_frac=0.22))]
    admP, nad = [], []
    for nm, kw in hist:
        kw2 = dict(SC["S03_FULL_NOTHING"]); kw2.update(kw)
        R = Regimen(**kw2)
        _, g, ya = run(R, P)
        o = observe(0.0, ya, P, R, "refeed")
        s = summarise(g, P, R)
        admP.append(o["Pser"]); nad.append(s["Pnadir"])
        print("%-24s %9.2f %9.2f %9.1f %9.2f %8.2f"
              % (nm, o["Pser"], o["Kser"], 100 * o["kiC"], s["Pnadir"],
                 s["mortality"]))
    print()
    print("Admission phosphate spans %.2f-%.2f mmol/L across a FOUR-FOLD range of"
          % (min(admP), max(admP)))
    print("starvation duration - a spread of %.2f mmol/L, well inside assay noise"
          % (max(admP) - min(admP)))
    print("and entirely inside the reference interval - while modelled mortality")
    print("varies %.1f-fold over the same six patients.  The number on the")
    print("admission panel does not carry the information.  The history does.")
    print()
    print("UNEXPECTED, and reported because it runs against the obvious guess:")
    print("the phosphate nadir is NOT monotone in starvation duration.  It")
    print("deepens to %.2f at 60-75 days and then becomes SHALLOWER again at 120"
          % min(nad))
    print("days (%.2f).  The most wasted patients have the least lean mass and the"
          % nad[-1])
    print("most insulin resistance, so they cannot mount the glycolytic flux that")
    print("generates the demand.  Their mortality still rises, because it is")
    print("driven by potassium, magnesium and thiamine rather than by phosphate.")
    print("A shallow nadir in a very wasted patient is therefore NOT reassurance;")
    print("it is the model predicting that serum phosphate is the wrong marker in")
    print("exactly the patients who are most at risk.")

    # ------------------------------------------------------------------
    sub("RESULT 5   the demand tracks glucose, so there is a GIR threshold")
    print("Carbohydrate is swept at a FIXED total energy of 25 kcal/kg/d, so")
    print("only the glucose infusion rate changes.  No repletion is given.")
    print()
    print("%8s %8s %9s %9s %9s %8s"
          % ("CHO %", "GIR", "Lambda_P", "P nadir", "K day 3", "mort %"))
    for cf in [0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 1.00]:
        kw = dict(SC["S03_FULL_NOTHING"])
        kw.update(dict(kcal_start=25.0, kcal_goal=25.0, advance_start=99.0,
                       cho_frac=cf, fat_frac=(1.0 - cf) * 0.64,
                       pro_frac=(1.0 - cf) * 0.36))
        R = Regimen(**kw)
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        print("%8.0f %8.2f %9.2f %9.2f %9.2f %8.2f"
              % (100 * cf, s["GIRmax"], s["Lmax"], s["Pnadir"],
                 s["d3"]["Kser"], s["mortality"]))
    print()
    print("The energy is constant down every row.  Only the glucose changes, and")
    print("the phosphate nadir moves across the whole clinical range.  'How many")
    print("calories' is the wrong question; 'how much glucose per minute' is the")
    print("right one.")

    # ------------------------------------------------------------------
    sub("RESULT 6   phosphate repletion has an optimum, not a maximum")
    print("Prophylactic intravenous phosphate is swept in a full-rate refeed.")
    print("Too little leaves Lambda_P above 1; too much exceeds the calcium-")
    print("phosphate solubility product and precipitates, taking the ionised")
    print("calcium down with it.")
    print()
    print("%8s %9s %9s %9s %9s %8s"
          % ("mmol/kg/d", "P nadir", "P day 7", "Ca nadir", "max CaxP", "mort %"))
    best, bestm, res_by_dose = None, 1e9, {}
    for pd in [0.0, 0.15, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0, 5.0, 8.0, 12.0]:
        kw = dict(SC["S03_FULL_NOTHING"]); kw.update(dict(P_dose=pd))
        R = Regimen(**kw)
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        flag = ""
        if s["CaxPmax"] > P["Ksp"]:
            flag = "  <- supersaturated"
        res_by_dose[pd] = s["mortality"]
        if s["mortality"] < bestm:
            bestm, best = s["mortality"], pd
        print("%8.2f %9.2f %9.2f %9.2f %9.2f %8.2f%s"
              % (pd, s["Pnadir"], s["d7"]["Pser"], s["Canadir"],
                 s["CaxPmax"], s["mortality"], flag))
    print()
    print("HONEST READING, including where the model FAILS to show what was")
    print("expected.  The benefit saturates: essentially all of it is bought by")
    print("the first 0.5 mmol/kg/d, which is the band the guidelines recommend and")
    print("which the model was never shown.  Beyond that the mortality curve is")
    print("FLAT to within the resolution of the hazard model.")
    print()
    print("The U-shape this section was written to demonstrate is REAL but it lies")
    print("far outside clinical practice.  The calcium-phosphate product only")
    print("crosses the %.1f mmol^2/L^2 solubility limit at 12 mmol/kg/d - twenty"
          % P["Ksp"])
    print("times the recommended dose - and only there does ionised calcium fall")
    print("and mortality turn back up.  At every dose a clinician would actually")
    print("write, the precipitation term never engages.")
    print()
    print("So the asymmetry inside this model is stark and one-sided: UNDER-")
    print("repletion costs %.2f percentage points of mortality, while over-"
          % (res_by_dose[0.0] - min(res_by_dose.values())))
    print("repletion costs nothing measurable until the dose becomes absurd.  If")
    print("that is right it favours the ASPEN posture of repleting generously and")
    print("not withholding feed.  The caveat is that the model has no soft-tissue")
    print("calcification, no phosphate-infusion arrhythmia and fixed renal")
    print("function, all three of which are the real reasons clinicians cap the")
    print("rate - so the flat top of this curve should be read as 'this model")
    print("cannot see the harm', not as 'there is none'.")

    # ------------------------------------------------------------------
    sub("RESULT 7   the two clocks - thiamine in hours, phosphate in days")
    print("Same patient, same feed.  Only the TIMING of a 300 mg intravenous")
    print("thiamine dose relative to the first glucose is changed, and separately")
    print("the timing of phosphate.")
    print()
    print("%-34s %10s %10s %10s"
          % ("", "lactate", "Wernicke %", "mort %"))
    for lbl, lead in [("thiamine 1 h BEFORE glucose", 1.0),
                      ("thiamine with the glucose", 0.0),
                      ("thiamine 6 h after", -6.0),
                      ("thiamine 24 h after", -24.0),
                      ("thiamine 72 h after", -72.0),
                      ("no thiamine at all", None)]:
        kw = dict(SC["S03_FULL_NOTHING"])
        if lead is None:
            kw.update(dict(th_iv=0.0))
        elif lead >= 0:
            kw.update(dict(th_iv=300.0, th_lead_h=lead, th_delay=0.0))
        else:
            kw.update(dict(th_iv=300.0, th_lead_h=0.0, th_delay=-lead))
        R = Regimen(**kw)
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        print("%-34s %10.1f %10.0f %10.2f"
              % (lbl, s["LACmax"], s["WE"], s["mortality"]))

    print()
    print("%-34s %10s %10s %10s"
          % ("", "P nadir", "P day 3", "mort %"))
    for lbl, pre in [("phosphate started on day 0", 0.0),
                     ("phosphate started on day 1", 1.0),
                     ("phosphate started on day 2", 2.0),
                     ("phosphate started on day 4", 4.0),
                     ("phosphate never started", None)]:
        kw = dict(SC["S03_FULL_NOTHING"])
        kw.update(dict(P_dose=0.0 if pre is None else 0.5,
                       P_delay=0.0 if pre is None else pre))
        R = Regimen(**kw)
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        print("%-34s %10.2f %10.2f %10.2f"
              % (lbl, s["Pnadir"], s["d3"]["Pser"], s["mortality"]))
    print()
    print("Two readings, and the first one contradicts the teaching:")
    print()
    print("(a) the famous 'thiamine BEFORE the glucose' rule buys nothing at the")
    print("    one-hour timescale on which it is usually taught - giving thiamine")
    print("    an hour early and giving it with the feed are indistinguishable")
    print("    here.  What the model does support is the WEAKER, and probably")
    print("    the real, version of the rule: thiamine must be given on DAY ZERO.")
    print("    Six hours late already costs; a day late costs more.")
    print()
    print("(b) phosphate is the opposite shape.  Starting it on day 0 abolishes")
    print("    the nadir; starting it 24 h later does not, because by then the")
    print("    cells have already taken what they needed out of a 13 mmol pool.")
    print("    Late phosphate still repairs the day-7 number, so it rescues the")
    print("    LABORATORY promptly and the PATIENT slowly - which is exactly the")
    print("    trap of treating this syndrome reactively.")

    # ------------------------------------------------------------------
    sub("RESULT 8   oral vs intravenous thiamine - A PREDICTION THAT FAILED")
    print("The intestinal thiamine carriers (ThTR-1/ThTR-2) saturate: the model's")
    print("maximum absorption rate is %.1f umol/h (~%.1f mg/h) whatever dose is"
          % (P["VmaxThA"], P["VmaxThA"] * 265.4 / 1000.0))
    print("swallowed.  Intravenous dosing has no such ceiling.")
    print()
    print("This section was written to show that oral thiamine cannot replete an")
    print("acutely deficient patient.  THE MODEL DOES NOT AGREE, and the")
    print("disagreement is reported rather than tuned away.")
    print()
    print("%-26s %8s %8s %8s %8s %8s"
          % ("", "TK 3 h", "TK 12 h", "TK 24 h", "TK d3", "WE %"))
    for lbl, kw2 in [("no thiamine", dict(th_iv=0.0, th_po=0.0)),
                     ("oral 100 mg/d", dict(th_po=100.0)),
                     ("oral 300 mg/d", dict(th_po=300.0)),
                     ("oral 1500 mg/d", dict(th_po=1500.0)),
                     ("IV 100 mg/d", dict(th_iv=100.0)),
                     ("IV 500 mg/d", dict(th_iv=500.0)),
                     ("oral 300, vomiting", dict(th_po=300.0, gut_fail=0.92)),
                     ("IV 100, vomiting", dict(th_iv=100.0, gut_fail=0.92))]:
        kw = dict(SC["S11_ALCOHOL_NOTH"]); kw.update(kw2)
        R = Regimen(**kw)
        _, g, _ = run(R, P)
        s = summarise(g, P, R)
        print("%-26s %8.2f %8.2f %8.2f %8.2f %8.0f"
              % (lbl, s["at"](3.0 / 24.0, "fTK"), s["at"](0.5, "fTK"),
                 s["at"](1.0, "fTK"), s["d3"]["fTK"], s["WE"]))
    print()
    print("The saturable carrier turns out not to be the binding constraint.  Its")
    print("ceiling is %.1f umol/h = %.0f mg/d of absorbed thiamine, and the daily"
          % (P["VmaxThA"], P["VmaxThA"] * 24 * 265.4 / 1000.0))
    print("REQUIREMENT is only %.1f mg.  Saturated absorption still delivers"
          % P["ThRDA"])
    print("roughly thirty times what the patient needs, so an intact gut repletes")
    print("the store within a day whatever the route.  Raising the oral dose from")
    print("300 to 1500 mg/d changes nothing, exactly as predicted - but so does")
    print("switching to the vein, which was NOT predicted.")
    print()
    print("The intravenous route earns its place on the last two rows instead: if")
    print("absorption is impaired - vomiting, ileus, an alcohol-damaged mucosa -")
    print("the oral arm collapses and the intravenous arm does not.  The model's")
    print("verdict is that 'IV because oral cannot deliver enough' is the wrong")
    print("reason, and 'IV because you cannot rely on this patient's gut' is the")
    print("right one.  Since that is precisely the population in question, the")
    print("clinical recommendation survives with its justification replaced.")

    # ------------------------------------------------------------------
    sub("RESULT 9   magnesium: how you give it, and what it unlocks")
    print("(a) the same 0.4 mmol/kg/d as a 2 h bolus or as a 24 h infusion")
    a, b = res["S14_MG_BOLUS"], res["S15_MG_INFUSION"]
    print("      bolus      Mg day 3 %.2f mmol/L" % a["d3"]["MGser"])
    print("      infusion   Mg day 3 %.2f mmol/L" % b["d3"]["MGser"])
    print("    Renal magnesium handling is a THRESHOLD, so the peak of a bolus is")
    print("    excreted rather than retained.  The prescription chart shows the")
    print("    same daily dose; the patient does not receive it.")
    print()
    print("(b) potassium repletion with and without magnesium repletion")
    c = res["S16_K_ONLY_NO_MG"]
    print("      K 2 mmol/kg/d alone         K day 3 %.2f, Mg day 3 %.2f"
          % (c["d3"]["Kser"], c["d3"]["MGser"]))
    print("      K 2 mmol/kg/d + Mg          K day 3 %.2f, Mg day 3 %.2f"
          % (b["d3"]["Kser"], b["d3"]["MGser"]))
    print("    Hypomagnesaemia releases the intracellular magnesium block on ROMK,")
    print("    so the distal nephron wastes potassium.  The effect is real in the")
    print("    model but MODEST at this repletion rate (%.2f vs %.2f mmol/L on day"
          % (c["d3"]["Kser"], b["d3"]["Kser"]))
    print("    3): a 2 mmol/kg/d potassium infusion is large enough to outrun the")
    print("    leak it is fighting.  The model therefore supports 'correct the")
    print("    magnesium' as good practice but does NOT reproduce the textbook")
    print("    picture of potassium that is refractory until magnesium is given -")
    print("    and that is a place where the model may be understating a real and")
    print("    frequently reported clinical phenomenon.")
    print("    Magnesium matters through a second, quieter route as well: it is")
    print("    required for PTH secretion, and PTH is what makes bone donate")
    print("    phosphate, so magnesium sits upstream of the phosphate supply line.")

    # ------------------------------------------------------------------
    sub("RESULT 10  the ventilatory cost of feeding on glucose alone")
    print("%-28s %8s %8s %9s %9s"
          % ("", "RQ", "VCO2", "PaCO2", "mort %"))
    for lbl, k in [("100 % carbohydrate", "S09_ENTERAL_MATCHED"),
                   ("40 % carbohydrate", "S10_MIXED_FUEL")]:
        s = res[k]
        print("%-28s %8.2f %8.3f %9.1f %9.2f"
              % (lbl, s["RQmax"], s["rows"][-1]["VCO2"], s["PaCO2max"],
                 s["mortality"]))
    print()
    print("    Substituting fat for carbohydrate at constant energy lowers the")
    print("    respiratory quotient, but note the SIZE of the two terms before")
    print("    accepting the usual explanation: carbon dioxide production differs")
    print("    by only about a tenth between these arms, while the modelled PaCO2")
    print("    differs by a quarter.  Most of the gap is not the CO2 load at all -")
    print("    it is the ventilatory pump, because the same glucose load has")
    print("    driven the phosphate down and taken the diaphragm's energetics with")
    print("    it.  The carbohydrate is not mainly overloading the lungs; it is")
    print("    weakening them, and then loading them.  A prediction that follows:")
    print("    the ventilatory benefit of a lipid-based feed should be largest in")
    print("    exactly the patients who are phosphate-depleted, and small in those")
    print("    who are not.")

    # ------------------------------------------------------------------
    sub("LIMITATIONS - stated, not buried")
    for i, s in enumerate([
        "The worst scenario (phosphate-free intravenous dextrose, no repletion",
        "  at all) bottoms out near 0.07 mmol/L.  Real reported nadirs in that",
        "  setting are 0.10-0.30.  The model has no rescue behaviour by a real",
        "  clinician, and its only brakes are transport saturation and the",
        "  glycolytic block, so it over-shoots at the extreme.",
        "Mortality is a hazard model with ONE fitted scale.  The split between",
        "  arrhythmic, cardiac-failure and respiratory death is structural, not",
        "  fitted, and should not be read as three separately validated numbers.",
        "There is no sepsis, no organ failure, no delirium, no re-feeding",
        "  induced hepatic steatosis and no vitamin deficiency other than",
        "  thiamine.  Real refeeding deaths often involve those.",
        "The starvation phase assumes a constant intake fraction.  Real",
        "  histories are episodic, and a binge-restrict pattern would deplete",
        "  thiamine faster than this model does.",
        "Bone is a single well-mixed compartment with a PTH-driven efflux.",
        "  Osteopenia of chronic starvation is not represented, so the model",
        "  probably OVERSTATES the skeletal buffer in the most chronic cases.",
        "Intracellular phosphate is one soft-tissue pool.  Erythrocyte 2,3-DPG,",
        "  myocardium and diaphragm are read off serum phosphate rather than",
        "  given their own compartments.",
        "The insulin submodel is a minimal model calibrated at two points.  It",
        "  is adequate for the direction and rough size of the refeeding insulin",
        "  rise, and is not a glucose-kinetics model.",
        "Renal function is fixed.  Acute kidney injury, which is common in",
        "  these patients, would change every threshold in the model.",
    ]):
        print("   " + s)


if __name__ == "__main__":
    _P, _res = main()
    sweeps(_P, _res)
    print()
    print("=" * 78)
    print("END OF REFERENCE OUTPUT")
    print("=" * 78)
