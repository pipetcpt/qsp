#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cgd_kernel.py — the phagosomal oxidant-chemistry kernel of the CGD QSP model.

WHAT THIS FILE IS
-----------------
Chronic granulomatous disease (CGD) is, at bottom, a disease of ONE number:
how many electrons per second the NADPH oxidase (NOX2/CYBB + p22phox/CYBA +
p47phox/NCF1 + p67phox/NCF2 + p40phox/NCF4 + Rac2) pushes across the phagosomal
membrane.  Everything that follows — superoxide, H2O2, HOCl, phagosomal pH,
phagosomal K+, granule protease release, bacterial and fungal killing, and
eventually granuloma, colitis and death — is downstream of that single flux.

This kernel therefore does NOT fit CGD.  It writes down the phagosome of a
NORMAL neutrophil from published rate constants and published concentrations,
and then turns the oxidase flux down by a factor phi and reads off what
happens.  CGD is the phi = 0 limit; hypomorphic CGD, X-linked carrier
mosaicism, post-transplant donor chimerism and lentiviral gene therapy are all
just different ways of writing the same phi.

    ONLY THE NORMAL PHAGOSOME IS CALIBRATED.  THE DISEASE IS PREDICTED.

THE ONE OBJECT THIS KERNEL EXISTS TO PRODUCE
---------------------------------------------
K(phi) — the log10 kill of one ingested organism in 60 minutes as a function of
the fractional oxidase activity of the neutrophil that ate it.  K(phi) is the
ONLY thing that crosses from this file into the whole-body model in
cgd_python_reference.py.  Its SHAPE — how sharp the threshold is — decides
every clinical question in this disease: how much donor chimerism cures a
transplant recipient, how many gene-marked cells a lentiviral product must
deliver, why an X-linked carrier with 20% normal neutrophils is well and one
with 5% is not, and whether the DHR-123 flow assay's single mean-fluorescence
number is a sufficient statistic for the patient in front of you.
(It is not.  See cgd_analysis.py, Section D.)

WHAT IS SPENT AND WHERE
------------------------
FIXED from the literature, never fitted (each cited in cgd_references.md):
    V_phag 1.2 fL, V_bact 0.5 fL                     Winterbourn 2006
    MPO 1.0 mM haem in the phagosome                 Winterbourn 2006
    R_ox 2.0 mM/s superoxide into the phagosome      Winterbourn 2006
    Cl- 0.10 M initial, resupplied                   Chapman 2002
    k1  MPO + H2O2  -> Cpd I          2.3e7 /M/s     Furtmuller 1998
    k2  Cpd I + Cl- -> MPO + HOCl     2.5e4 /M/s     Furtmuller 2000
    k3  Cpd I + AH  -> Cpd II         2.0e7 /M/s     Marquez 1995
    k4  Cpd II + O2.- -> MPO          1.2e6 /M/s     Kettle 2007
    k5  Cpd II + AH -> MPO            1.0e4 /M/s     Marquez 1995
    k6  MPO + O2.-  -> Cpd III        1.1e6 /M/s     Kettle 1988/2007
    k7  Cpd III + O2.- -> MPO + H2O2  1.5e6 /M/s     Kettle 1988/2007
    kdis spontaneous dismutation      5.0e5 /M/s @pH7 Bielski 1985
    kcat bacterial catalase           3.4e7 /M/s     Switala 2002
    kHOCl + Met/Cys                   3.0e7 /M/s     Storkey 2014
    HOCl protein sink 20 mM           phagosomal protein ~100 g/L

CALIBRATED here, and ONLY against measurements made on NORMAL neutrophils —
five numbers, five targets:
    AH0        -> 70% of phagosomal H2O2 ends as HOCl   Winterbourn 2006
    alk_gain   -> peak phagosomal pH 7.80               Segal 1981
    kK_comp    -> peak phagosomal K+ 0.50 M             Reeves 2002
    kox        -> normal PMN kills 95% of S. aureus/60 min
    kprot      -> phi=0 PMN kills 40% of S. aureus/60 min

The last of the five is the only number in this file taken from a CGD
measurement, and it buys exactly one thing: the SIZE of the oxidase-independent
killing arm.  Without it the model would have no way of knowing that a CGD
neutrophil is not simply an inert bag, which it plainly is not — CGD patients
are not neutropenic-equivalent and do not die of Streptococcus.

UNFITTED PREDICTIONS THE KERNEL THEN MAKES (checked in cgd_analysis.py):
    ~1e8 HOCl molecules delivered per ingested bacterium (Winterbourn 2006)
    the ratio of phagosomal MPO to bacterial catalase as competitors for H2O2
    the shape of K(phi), including phi50 and the Hill coefficient
    the whole organism-spectrum ordering of Winkelstein 2000

TWO ARMS THAT ARE CONTESTED IN THE LITERATURE AND ARE THEREFORE SWITCHABLE
---------------------------------------------------------------------------
arms['segal_pH']  Segal 1981 has the normal phagosome alkalinising to ~7.8 and
                  the CGD phagosome staying acid.  Jankowski/Grinstein 2002
                  measured essentially normal acidification in CGD and dispute
                  it.  Base case follows Segal; Section H prices the other one.
arms['protease']  Reeves 2002 has K+ influx, not oxidants, doing the killing by
                  releasing cationic proteases from the sulphated proteoglycan
                  granule matrix.  Segal's critics dispute the magnitude.  Base
                  case includes it; Section H runs without it.
"""

import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

NA = 6.02214076e23

# ----------------------------------------------------------------------------
# 1. FIXED PHYSICAL AND CHEMICAL CONSTANTS  (not fitted)
# ----------------------------------------------------------------------------
FIXED = dict(
    V_phag = 1.2e-15,      # L, phagosome around one coccus     Winterbourn 2006
    V_bact = 0.50e-15,     # L, one S. aureus                   Winterbourn 2006

    R_ox   = 2.0e-3,       # M/s superoxide generation at phi=1 Winterbourn 2006
    tau_ox = 600.0,        # s, decay of the phagosomal burst

    MPO_tot = 1.0e-3,      # M haem                             Winterbourn 2006
    Cl0     = 0.100,       # M                                  Chapman 2002
    kCl_in  = 0.02,        # /s chloride resupply (CLC-3)
    kAH_in  = 0.01,        # /s resupply of 1-electron donors

    k1 = 2.3e7,  k2 = 2.5e4,  k3 = 2.0e7,  k4 = 1.2e6,
    k5 = 1.0e4,  k6 = 1.1e6,  k7 = 1.5e6,  k7b = 0.25,

    kdis_pH7 = 5.0e5,      # /M/s spontaneous dismutation at pH 7.0

    # Superoxide is a charged species and does not cross the bacterial envelope,
    # so bacterial SOD — which protects the organism against its OWN endogenous
    # superoxide — has essentially no access to the phagosomal lumen.  Set to
    # zero in the base case; a non-zero value is a sensitivity, not biology.
    k_sod_b  = 2.0e9,
    SOD_b    = 0.0,

    kcat_bact = 3.4e7,     # /M/s per catalase haem
    kprx_bact = 4.0e7,     # /M/s per peroxiredoxin (AhpC) site
    PRX_b     = 2.0e-5,    # M referenced to phagosome volume

    kHOCl_sink = 3.0e7,    # /M/s HOCl + Met/Cys
    Tgt0       = 2.0e-2,   # M oxidisable protein residues in the phagosome
    f_bact     = 0.30,     # fraction of the HOCl flux intercepted by the organism

    pH_min  = 6.00,        # V-ATPase floor reached when phi = 0   Segal 1981
    pH_rest = 7.20,        # phagosomal pH at t = 0
    kv_atp  = 0.020,       # /s V-ATPase relaxation toward pH_min

    K0      = 0.140,       # M phagosomal K+ at t = 0
    kK_out  = 0.006,       # /s K+ leak back out

    # Non-oxidative killing: TWO arms with OPPOSITE pH dependence.  This is not
    # decoration.  The alkaline serine proteases (elastase, cathepsin G,
    # proteinase 3) need both K+ release and an alkaline phagosome, so CGD
    # closes them twice over; the acid hydrolases (cathepsins B/D/L) and the
    # alpha-defensins HNP1-4 work BETTER in the acid CGD phagosome.  A model
    # with only the alkaline arm predicts that a CGD neutrophil is inert, and a
    # CGD neutrophil is not inert.
    pH_opt_alk  = 8.0,
    pH_mid_acid = 6.5,
    K50_prot    = 0.250,   # M K+ needed to strip cationic proteases off the
    n_prot      = 3.0,     #   anionic proteoglycan matrix (Reeves 2002)
    w_alk       = 0.50,    # split between the two non-oxidative arms (assumed)

    # Doubling time INSIDE a phagosome, not in broth.  A first draft used the
    # 30-minute broth doubling time; at that rate intraphagosomal growth eats a
    # third of the normal neutrophil's whole killing rate, the oxidative arm has
    # to be enormous to overcome it, and switching the non-oxidative arm off
    # then collapsed normal killing from 95% to 7% — an absurdity that exposed
    # the error.  Three hours is the defensible figure for a nutrient-restricted,
    # acidifying phagosome.
    mu_bact = math.log(2) / (3 * 3600.0),  # /s
)

# Calibrated on the NORMAL phagosome only — filled in by calibrate().
CAL = dict(AH0=None, alk_gain=None, kK_comp=None, kox=None, kprot=None,
           prot_norm_alk=1.0, prot_norm_acid=1.0)

TARGETS = dict(
    frac_H2O2_to_HOCl = 0.70,   # Winterbourn 2006
    pH_peak_normal    = 7.80,   # Segal 1981 Nature
    K_peak_normal     = 0.500,  # Reeves 2002 Nature
    kill60_normal     = 0.95,   # fraction of S. aureus killed in 60 min, normal
    kill60_cgd        = 0.40,   # fraction killed in 60 min, phi = 0
)

ARMS_DEFAULT = dict(catalase=True, perox=True, protease=True, segal_pH=True)


# ----------------------------------------------------------------------------
# 2. THE ORGANISM TABLE — and the two columns the textbook collapses into one
# ----------------------------------------------------------------------------
# catalase : intrabacterial catalase haem, M (concentration INSIDE the organism)
# perox    : the organism's OWN H2O2 output, M/s referenced to phagosome volume
# omega    : fraction of this organism's killing that is oxidant-dependent
#
# "Catalase-positive organisms cause CGD infections" is a statement about
# column 1.  The MECHANISM it is supposed to encode — that catalase-negative
# organisms hand the CGD phagosome the H2O2 it cannot make for itself — is a
# statement about column 2.  These are not the same column.  Section B of the
# analysis shows the disease tracks column 2 and is indifferent to column 1.
#
# Peroxigenicity anchor: S. pneumoniae makes ~1 mM H2O2 in 30 min from 1e8
# CFU/mL (Pericone 2000) = 5.6e-18 mol/cell/s, which delivered into a 1.2 fL
# phagosome is 4.6e-3 M/s — larger than the NADPH oxidase's own 1.0e-3 M/s of
# H2O2.  A single pneumococcus arms the phagosome better than NOX2 does.
ORGANISMS = {
    "Staphylococcus aureus":      dict(catalase=5.0e-5, perox=0.0,    omega=0.72, cgd=True),
    "Serratia marcescens":        dict(catalase=4.0e-5, perox=0.0,    omega=0.80, cgd=True),
    "Burkholderia cepacia":       dict(catalase=6.0e-5, perox=0.0,    omega=0.88, cgd=True),
    "Nocardia spp.":              dict(catalase=5.0e-5, perox=0.0,    omega=0.90, cgd=True),
    "Aspergillus fumigatus":      dict(catalase=7.0e-5, perox=0.0,    omega=0.92, cgd=True),
    "Salmonella spp.":            dict(catalase=3.0e-5, perox=0.0,    omega=0.70, cgd=True),
    "Granulibacter bethesdensis": dict(catalase=3.0e-5, perox=0.0,    omega=0.85, cgd=True),
    "Escherichia coli":           dict(catalase=3.0e-5, perox=0.0,    omega=0.45, cgd=False),
    "Candida albicans":           dict(catalase=4.0e-5, perox=0.0,    omega=0.60, cgd=False),
    "Streptococcus pneumoniae":   dict(catalase=0.0,    perox=4.6e-3, omega=0.55, cgd=False),
    "Streptococcus pyogenes":     dict(catalase=0.0,    perox=1.5e-3, omega=0.50, cgd=False),
    "Enterococcus faecalis":      dict(catalase=0.0,    perox=8.0e-4, omega=0.40, cgd=False),
    "Lactobacillus spp.":         dict(catalase=0.0,    perox=3.0e-3, omega=0.35, cgd=False),
    "Haemophilus influenzae":     dict(catalase=2.0e-5, perox=0.0,    omega=0.30, cgd=False),
}

# CGD infection frequency (percentage of patients from whom the organism was
# ever recovered), Winkelstein 2000 (n = 368) with Marciano 2015.  Used ONLY to
# test the model's ordering.  Never fitted.
WINKELSTEIN_FREQ = {
    "Aspergillus fumigatus":     33.0,
    "Staphylococcus aureus":     23.0,
    "Burkholderia cepacia":       7.0,
    "Serratia marcescens":        7.0,
    "Nocardia spp.":              5.0,
    "Salmonella spp.":            4.0,
    "Escherichia coli":           2.0,
    "Candida albicans":           2.0,
    "Streptococcus pneumoniae":   0.8,
    "Streptococcus pyogenes":     0.5,
    "Haemophilus influenzae":     0.5,
}


# ----------------------------------------------------------------------------
# 3. THE PHAGOSOME ODE SYSTEM
# ----------------------------------------------------------------------------
STATES = ["S", "H", "MPO", "CI", "CII", "CIII", "Cl", "AH", "Tgt",
          "pH", "K", "B", "cumHOCl", "cumHOCl_b", "cumH2O2", "cumO2m"]
IDX = {n: i for i, n in enumerate(STATES)}
NSTATE = len(STATES)


def protease_arms(K, pH, F=FIXED):
    """Raw (un-normalised) activities of the two non-oxidative killing arms."""
    pH = min(max(float(pH), 3.0), 11.0)
    K = max(float(K), 0.0)
    gate_K   = K ** F["n_prot"] / (F["K50_prot"] ** F["n_prot"] + K ** F["n_prot"])
    gate_alk = 1.0 / (1.0 + 10.0 ** (F["pH_opt_alk"] - pH))
    P_alk    = gate_K * gate_alk
    P_acid   = 1.0 / (1.0 + 10.0 ** (pH - F["pH_mid_acid"]))
    return P_alk, P_acid


def protease_activity(K, pH, cal, F=FIXED):
    """Total non-oxidative killing capacity, normalised to 1.0 in the NORMAL
    phagosome.  Both arms are normalised separately so that each contributes
    w_alk / (1 - w_alk) of the normal cell's non-oxidative kill."""
    P_alk, P_acid = protease_arms(K, pH, F)
    return (F["w_alk"] * P_alk / cal["prot_norm_alk"] +
            (1.0 - F["w_alk"]) * P_acid / cal["prot_norm_acid"])


def _rhs(t, y, phi, org, F, cal, arms):
    (S, H, MPO, CI, CII, CIII, Cl, AH, Tgt, pH, K, B, _c1, _c2, _c3, _c4) = y

    S = max(S, 0.0);   H = max(H, 0.0);   MPO = max(MPO, 0.0)
    CI = max(CI, 0.0); CII = max(CII, 0.0); CIII = max(CIII, 0.0)
    Cl = max(Cl, 0.0); AH = max(AH, 0.0);  Tgt = max(Tgt, 0.0); B = max(B, 0.0)

    # --- NADPH oxidase: the one number CGD changes --------------------------
    v_ox = phi * F["R_ox"] * math.exp(-t / F["tau_ox"])

    # --- superoxide -> H2O2 --------------------------------------------------
    kdis  = F["kdis_pH7"] * 10.0 ** (7.0 - pH)      # Bielski 1985 pH dependence
    v_dis = kdis * S * S                             # M/s H2O2 from 2 O2.-
    v_sod = F["k_sod_b"] * F["SOD_b"] * S            # off in the base case

    # --- MPO peroxidase cycle ------------------------------------------------
    # k6/k7 together make MPO the phagosome's real superoxide dismutase:
    #   MPO(III) + O2.-      -> Cpd III
    #   Cpd III  + O2.- + 2H+-> MPO(III) + H2O2 + O2
    # i.e. 2 O2.- -> H2O2, catalysed.  Dropping this pair (as an SOD-free model
    # would) starves MPO of the very H2O2 it needs and understates HOCl by ~40x.
    v1 = F["k1"] * MPO * H
    v2 = F["k2"] * CI * Cl
    v3 = F["k3"] * CI * AH
    v4 = F["k4"] * CII * S
    v5 = F["k5"] * CII * AH
    v6 = F["k6"] * MPO * S
    v7 = F["k7"] * CIII * S
    v7b = F["k7b"] * CIII

    # --- the organism's own peroxide chemistry -------------------------------
    cat_eff = org["catalase"] * (F["V_bact"] / F["V_phag"]) if arms["catalase"] else 0.0
    v_cat = F["kcat_bact"] * cat_eff * H
    v_prx = F["kprx_bact"] * F["PRX_b"] * H
    v_bH  = (org["perox"] * min(B, 1.0)) if arms["perox"] else 0.0

    v_HOCl   = v2                          # HOCl is consumed within microseconds
    v_HOCl_b = F["f_bact"] * v_HOCl        # the share that reaches the organism

    dS    = v_ox - 2.0 * v_dis - v_sod - v4 - v6 - v7
    dH    = v_dis + 0.5 * v_sod + v7 + v_bH - v1 - v_cat - v_prx
    dMPO  = -v1 + v2 + v4 + v5 - v6 + v7 + v7b
    dCI   =  v1 - v2 - v3
    dCII  =  v3 - v4 - v5
    dCIII =  v6 - v7 - v7b
    dCl   = -v2 + F["kCl_in"] * (F["Cl0"] - Cl)
    # AH is the pool of one-electron donors that divert Compound I away from
    # chlorination — free tyrosine, urate, ascorbate, and above all the tyrosine
    # residues of the ~100 g/L of protein crammed into the phagosome.  It is
    # held FIXED, not depleted: at 100 g/L, protein tyrosine alone is several mM
    # and the whole 10-minute burst consumes micromolar.  Letting it deplete (as
    # a first draft of this file did) silently removes the entire chlorination /
    # peroxidation competition and pins the HOCl yield at 0.90 regardless of AH0.
    dAH   = 0.0
    dTgt  = -(1.0 - F["f_bact"]) * v_HOCl

    # --- phagosomal pH -------------------------------------------------------
    # Every H2O2 formed from 2 O2.- consumes 2 H+, so the alkalinisation is
    # driven by the TOTAL dismutation flux, spontaneous plus MPO-catalysed.
    if arms["segal_pH"]:
        dpH = cal["alk_gain"] * (v_dis + v7) - F["kv_atp"] * (pH - F["pH_min"])
    else:
        dpH = -F["kv_atp"] * (pH - F["pH_min"])     # Jankowski/Grinstein variant

    # --- phagosomal K+, charge compensation for the electron flux ------------
    dK = cal["kK_comp"] * v_ox - F["kK_out"] * (K - F["K0"])

    # --- killing -------------------------------------------------------------
    P_non = protease_activity(K, pH, cal, F) if arms["protease"] else 0.0
    om = org["omega"]
    h_ox  = cal["kox"]   * (2.0 * om)         * (v_HOCl_b / 1.0e-4)
    h_non = cal["kprot"] * (2.0 * (1.0 - om)) * P_non
    dB = F["mu_bact"] * B - (h_ox + h_non) * B

    return [dS, dH, dMPO, dCI, dCII, dCIII, dCl, dAH, dTgt, dpH, dK, dB,
            v_HOCl, v_HOCl_b, v_dis + 0.5 * v_sod + v7 + v_bH, v_ox]


def run_phagosome(phi, org, t_end=3600.0, cal=None, arms=None, F=FIXED, n_out=181):
    cal = cal if cal is not None else CAL
    arms = dict(ARMS_DEFAULT, **(arms or {}))
    y0 = np.zeros(NSTATE)
    y0[IDX["MPO"]] = F["MPO_tot"]
    y0[IDX["Cl"]]  = F["Cl0"]
    y0[IDX["AH"]]  = cal["AH0"]
    y0[IDX["Tgt"]] = F["Tgt0"]
    y0[IDX["pH"]]  = F["pH_rest"]
    y0[IDX["K"]]   = F["K0"]
    y0[IDX["B"]]   = 1.0                      # one organism, relative units
    sol = solve_ivp(_rhs, (0.0, t_end), y0, args=(phi, org, F, cal, arms),
                    method="LSODA", rtol=1e-7, atol=1e-16,
                    t_eval=np.linspace(0.0, t_end, n_out), first_step=1e-9,
                    max_step=20.0)
    if not sol.success:
        raise RuntimeError("phagosome integration failed: " + sol.message)
    return sol


def summarise(sol, F=FIXED):
    y = sol.y
    o = dict(
        pH_peak      = float(np.max(y[IDX["pH"]])),
        pH_end       = float(y[IDX["pH"]][-1]),
        K_peak       = float(np.max(y[IDX["K"]])),
        S_peak       = float(np.max(y[IDX["S"]])),
        H_peak       = float(np.max(y[IDX["H"]])),
        Cl_min       = float(np.min(y[IDX["Cl"]])),
        CIII_frac    = float(np.max(y[IDX["CIII"]]) / F["MPO_tot"]),
        HOCl_total_M = float(y[IDX["cumHOCl"]][-1]),
        HOCl_bact_M  = float(y[IDX["cumHOCl_b"]][-1]),
        H2O2_total_M = float(y[IDX["cumH2O2"]][-1]),
        O2m_total_M  = float(y[IDX["cumO2m"]][-1]),
        B_end        = float(y[IDX["B"]][-1]),
    )
    o["HOCl_molecules_per_bact"] = o["HOCl_bact_M"] * F["V_phag"] * NA
    o["frac_H2O2_to_HOCl"] = (o["HOCl_total_M"] / o["H2O2_total_M"]
                              if o["H2O2_total_M"] > 0 else float("nan"))
    o["frac_O2m_to_H2O2"] = (2.0 * o["H2O2_total_M"] / o["O2m_total_M"]
                             if o["O2m_total_M"] > 0 else float("nan"))
    return o


def kill_fraction(phi, org, minutes=60.0, cal=None, arms=None, F=FIXED):
    sol = run_phagosome(phi, org, minutes * 60.0, cal, arms, F, n_out=13)
    return 1.0 - min(max(float(sol.y[IDX["B"]][-1]), 0.0), 1.0)


def K_log(phi, org, minutes=60.0, cal=None, arms=None, F=FIXED, floor=1e-14):
    sol = run_phagosome(phi, org, minutes * 60.0, cal, arms, F, n_out=13)
    return -math.log10(max(float(sol.y[IDX["B"]][-1]), floor))


# ----------------------------------------------------------------------------
# 4. CALIBRATION — five numbers, five NORMAL-neutrophil targets
# ----------------------------------------------------------------------------
def calibrate(verbose=True, F=FIXED):
    """Solved in dependency order.

    AH0 first, because the competition between chloride and the low-molecular-
    weight one-electron donors for Compound I sets the whole HOCl budget and
    therefore everything oxidative.  Then pH (which feeds back into the
    dismutation rate constant and into both protease gates), then K+, then the
    two killing constants.  kox and kprot are separable to better than 1%:
    at phi = 0 there is no HOCl at all, so the phi = 0 target sees kprot alone.
    """
    cal = dict(AH0=5e-5, alk_gain=41.0, kK_comp=0.0, kox=0.0, kprot=0.0,
               prot_norm_alk=1.0, prot_norm_acid=1.0)
    ref = ORGANISMS["Staphylococcus aureus"]

    # --- (i)+(ii) AH0 and alk_gain, iterated -------------------------------
    # They are weakly coupled through the pH dependence of the spontaneous
    # dismutation rate constant.  Two passes are enough: the second pass moves
    # AH0 by well under 1%.
    def f_AH(a):
        c = dict(cal); c["AH0"] = a
        s = summarise(run_phagosome(1.0, ref, 600.0, c, None, F, n_out=13), F)
        return s["frac_H2O2_to_HOCl"] - TARGETS["frac_H2O2_to_HOCl"]

    def f_pH(g):
        c = dict(cal); c["alk_gain"] = g
        s = run_phagosome(1.0, ref, 1800.0, c, None, F, n_out=181)
        return float(np.max(s.y[IDX["pH"]])) - TARGETS["pH_peak_normal"]

    for _pass in range(2):
        cal["AH0"] = brentq(f_AH, 1e-8, 1e-2, xtol=1e-14, rtol=1e-12)
        cal["alk_gain"] = brentq(f_pH, 1e-3, 1e7, xtol=1e-8, rtol=1e-10)

    # --- (iii) kK_comp: peak phagosomal K+ 0.50 M ----------------------------
    def f_K(g):
        c = dict(cal); c["kK_comp"] = g
        s = run_phagosome(1.0, ref, 1800.0, c, None, F, n_out=181)
        return float(np.max(s.y[IDX["K"]])) - TARGETS["K_peak_normal"]
    cal["kK_comp"] = brentq(f_K, 1e-4, 1e5, xtol=1e-12, rtol=1e-10)

    # --- normalise each non-oxidative arm to the NORMAL phagosome ------------
    s = run_phagosome(1.0, ref, 3600.0, cal, None, F, n_out=361)
    Kt, pHt = s.y[IDX["K"]], s.y[IDX["pH"]]
    arms_raw = np.array([protease_arms(k, p, F) for k, p in zip(Kt, pHt)])
    cal["prot_norm_alk"]  = float(np.mean(arms_raw[:, 0]))
    cal["prot_norm_acid"] = float(np.mean(arms_raw[:, 1]))

    # --- (iv) kprot from the phi = 0 target (no HOCl exists there) -----------
    def f_prot(kp):
        c = dict(cal); c["kox"] = 0.0; c["kprot"] = kp
        return kill_fraction(0.0, ref, 60.0, c, None, F) - TARGETS["kill60_cgd"]
    cal["kprot"] = brentq(f_prot, 1e-8, 1e3, xtol=1e-14, rtol=1e-12)

    # --- (v) kox from the phi = 1 target -------------------------------------
    def f_ox(ko):
        c = dict(cal); c["kox"] = ko
        return kill_fraction(1.0, ref, 60.0, c, None, F) - TARGETS["kill60_normal"]
    cal["kox"] = brentq(f_ox, 1e-10, 1e4, xtol=1e-16, rtol=1e-12)

    CAL.update(cal)
    if verbose:
        print("  phagosome kernel calibrated on NORMAL neutrophils only:")
        for k in ("AH0", "alk_gain", "kK_comp", "kox", "kprot",
                  "prot_norm_alk", "prot_norm_acid"):
            print(f"    {k:15s} = {cal[k]:.6g}")
    return cal


# ----------------------------------------------------------------------------
# 5. K(phi): the object the whole-body model consumes
# ----------------------------------------------------------------------------
def K_curve(org_name="Staphylococcus aureus", n=41, minutes=60.0,
            cal=None, arms=None, F=FIXED):
    org = ORGANISMS[org_name]
    phis = np.linspace(0.0, 1.0, n)
    return phis, np.array([K_log(p, org, minutes, cal, arms, F) for p in phis])


def fit_hill(phis, Ks):
    """K(phi) = K0 + (Kmax-K0) * phi^h / (phi50^h + phi^h).

    This Hill triple is the ENTIRE interface between the phagosome chemistry and
    the whole-body model, and it is what the mrgsolve file and the Shiny app
    carry instead of the 16-state chemistry.
    """
    from scipy.optimize import least_squares
    K0, Kmax = float(Ks[0]), float(Ks[-1])

    def resid(p):
        phi50, h = p
        pred = K0 + (Kmax - K0) * (phis ** h) / (phi50 ** h + phis ** h + 1e-30)
        return pred - Ks

    # phi50 is allowed ABOVE 1.  A fitted phi50 > 1 is not a failure: it is the
    # statement that killing has NO threshold inside the physiological range and
    # is still accelerating at 100% of normal oxidase activity.  That turns out
    # to be what this kernel says, and it is the reason the clinical threshold
    # in this disease cannot live in the phagosome.  See Section C.
    r = least_squares(resid, [0.5, 1.5], bounds=([1e-4, 0.2], [20.0, 20.0]))
    phi50, h = r.x
    pred = K0 + (Kmax - K0) * (phis ** h) / (phi50 ** h + phis ** h + 1e-30)
    ss = float(np.sum((pred - Ks) ** 2))
    tot = float(np.sum((Ks - np.mean(Ks)) ** 2))
    return dict(K0=K0, Kmax=Kmax, phi50=float(phi50), hill=float(h),
                r2=(1.0 - ss / tot) if tot > 0 else float("nan"),
                rmse=float(np.sqrt(ss / len(Ks))))


def K_hill(phi, hp):
    p = max(float(phi), 0.0)
    return hp["K0"] + (hp["Kmax"] - hp["K0"]) * p ** hp["hill"] / (
        hp["phi50"] ** hp["hill"] + p ** hp["hill"] + 1e-30)


def survival_fraction(phi, hp):
    """Fraction of ingested organisms SURVIVING one 60-min phagosome at phi.

    Population averages must be taken HERE, in survival space, never in the
    log-kill space K(phi).  Averaging log kills across a mixed neutrophil
    population is the arithmetic error that makes an X-linked carrier look like
    a hypomorphic hemizygote.
    """
    return 10.0 ** (-K_hill(phi, hp))


def population_survival(f_normal, phi_residual, hp):
    """A real neutrophil population is a MIXTURE, not a mean.

    f_normal      : fraction of neutrophils with a fully competent oxidase
    phi_residual  : oxidase activity of the remaining (1 - f_normal) cells

    The DHR-123 assay reported as a single mean channel fluorescence or
    stimulation index sees only  f*1 + (1-f)*r.  The phagosome sees the two
    populations separately.  Same number, different disease.
    """
    f = min(max(float(f_normal), 0.0), 1.0)
    return f * survival_fraction(1.0, hp) + (1.0 - f) * survival_fraction(phi_residual, hp)


def dhr_mean(f_normal, phi_residual):
    return f_normal * 1.0 + (1.0 - f_normal) * phi_residual


# ----------------------------------------------------------------------------
# 6. THE CATALASE ARITHMETIC
# ----------------------------------------------------------------------------
def h2o2_competition(org, F=FIXED):
    """Who gets the phagosomal H2O2 — MPO, or the organism's own catalase?
    Pseudo-first-order rate constants, /s.  The catalase dogma, reduced to a
    division."""
    cat_eff = org["catalase"] * (F["V_bact"] / F["V_phag"])
    k_mpo = F["k1"] * F["MPO_tot"]
    k_cat = F["kcat_bact"] * cat_eff
    k_prx = F["kprx_bact"] * F["PRX_b"]
    tot = k_mpo + k_cat + k_prx
    return dict(k_mpo=k_mpo, k_cat=k_cat, k_prx=k_prx,
                frac_to_mpo=k_mpo / tot, frac_to_catalase=k_cat / tot,
                catalase_M_needed_to_halve=k_mpo / (F["kcat_bact"] *
                                                    (F["V_bact"] / F["V_phag"])))


def self_arming_index(org, F=FIXED):
    """How much H2O2 the organism hands the phagosome, relative to what a NORMAL
    oxidase supplies (1.0 = the organism is as good an H2O2 source as NOX2).

    This is what the catalase dogma is actually about, and it is a DIFFERENT
    column of the organism table from catalase content."""
    oxidase_H2O2 = 0.5 * F["R_ox"]       # 2 O2.- -> 1 H2O2
    return org["perox"] / oxidase_H2O2


if __name__ == "__main__":
    calibrate()
    ref = ORGANISMS["Staphylococcus aureus"]
    print()
    for phi in (1.0, 0.5, 0.2, 0.1, 0.05, 0.0):
        s = summarise(run_phagosome(phi, ref), FIXED)
        print(f"phi={phi:4.2f}  pH_peak={s['pH_peak']:5.2f}  K_peak={s['K_peak']:.3f} M  "
              f"O2-={s['S_peak']:.2e}  H2O2={s['H_peak']:.2e}  "
              f"HOCl/bact={s['HOCl_molecules_per_bact']:.2e}  surv={s['B_end']:.4g}")
