#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pws_calibration.py — derived quantities, bifurcation structure and sensitivity
==============================================================================
Everything in this file is either CLOSED-FORM or a root-find on a closed form.
Nothing here is fitted.  It exists so that the claims in README.md can be
checked without running the 64-state simulation, and so that the boundary
between "calibrated" and "derived" is explicit.

    python3 pws_calibration.py > pws_calibration_output.txt
    python3 pws_calibration.py --sens     # also run the local sensitivity sweep
                                          # (this one DOES integrate the ODEs)

CONTENTS
  1  the escape-ratio formalism, swept over PC1/3 activity
  2  the exact conservation result (why a cross-reacting assay is blind)
  3  harmonic-mean sensitivity — which arm is worth repairing, and why
  4  the melanocortin relay CEILING — the bound MC4R agonism cannot pass
  5  ghrelin-arm elasticity — how little hyperghrelinaemia can be worth
  6  the food-seeking bifurcation diagram and its two saddle-node loci
  7  carbetocin: the analytic optimum, and the width of the useful window
  8  DCCR: why the therapeutic index is fixed by two EC50s
  9  the two airway clocks: the peak time in closed form
 10  energy-balance consistency of the reference scaffold
 11  local parameter sensitivity of the primary endpoints (--sens)
"""

import math
import sys

import pws_reference_model as M

P = M.P


def hr(c="-", n=78):
    return c * n


OUT = []
W = OUT.append


# ---------------------------------------------------------------------------
# 1.  THE ESCAPE-RATIO FORMALISM
# ---------------------------------------------------------------------------
def section1():
    W("=" * 78)
    W("PWS QSP MODEL — CALIBRATION AND DERIVED QUANTITIES")
    W("=" * 78)
    W("")
    W("[1] THE ESCAPE-RATIO FORMALISM, SWEPT OVER PC1/3 ACTIVITY")
    W(hr())
    W("  A prohormone P is made at rate S, processed by PC1/3 with pseudo-first-")
    W("  order rate kc*PC13, and lost unprocessed at rate d (secreted as")
    W("  precursor, or degraded).  In the sub-saturating regime:")
    W("")
    W("      dP/dt = S - (kc*PC13 + d) P          =>  P_ss = S / (PC13 + eps)")
    W("      product flux = kc*PC13*P             =>  flux = S / (1 + eps/PC13)")
    W("      eps = d/kc = d*Km/kcat               (the ESCAPE RATIO)")
    W("")
    W("      L(eps, PC13) = 1 - (1+eps)/(1+eps/PC13)     product loss")
    W("      R(eps, PC13) = (1+eps)/(PC13+eps)           precursor accumulation")
    W("")
    W("  Note the limits, which are the whole point:")
    W("      eps -> 0    L -> 0        (fully compensated, precursor piles up)")
    W("      eps -> inf  L -> 1-PC13   (fully proportional, no pile-up)")
    W("")
    branches = [("pro-oxytocin -> OXT", P["EPSPOXT"]),
                ("POMC -> alpha-MSH", P["EPSPOMC"]),
                ("pro-GHRH -> GHRH", P["EPSPGHRH"]),
                ("proinsulin -> insulin", P["EPSPINS"]),
                ("proghrelin -> acyl-ghrelin", P["EPSPGHR"])]
    W("  product loss L (%) as PC1/3 activity falls")
    W("  %-32s %7s %7s %7s %7s %7s %7s"
      % ("branch (eps)", "1.00", "0.80", "0.60", "0.40", "0.30", "0.20"))
    for nm, eps in branches:
        row = "  %-32s" % ("%s (%.2f)" % (nm, eps))
        for pc in (1.0, 0.8, 0.6, 0.4, 0.3, 0.2):
            row += " %7.1f" % (100.0 * M.branch_loss(eps, pc))
        W(row)
    W("")
    W("  precursor accumulation R (x normal) over the same sweep")
    W("  %-32s %7s %7s %7s %7s %7s %7s"
      % ("branch (eps)", "1.00", "0.80", "0.60", "0.40", "0.30", "0.20"))
    for nm, eps in branches:
        row = "  %-32s" % ("%s (%.2f)" % (nm, eps))
        for pc in (1.0, 0.8, 0.6, 0.4, 0.3, 0.2):
            row += " %7.3f" % M.branch_precursor(eps, pc)
        W(row)
    W("")
    W("  The two tables run in OPPOSITE directions down the rows.  At every")
    W("  PC1/3 activity the branch ordering by product loss is the exact")
    W("  reverse of the ordering by precursor accumulation — this is forced by")
    W("  the algebra, not chosen.  Consequence: the precursor:product ratio is")
    W("  a good severity marker only for the branches that are NOT failing.")
    W("")
    pc = M.pc13_activity(P)
    W("  AN EXACT IDENTITY, AND THE MEASUREMENT THEORY THAT FOLLOWS")
    W(hr("."))
    W("  Divide Eq. B by (1 - Eq. A):")
    W("")
    W("      precursor / product  =  (1+eps)/(PC13+eps) * (1+eps/PC13)/(1+eps)")
    W("                           =  (1+eps/PC13)/(PC13+eps)")
    W("                           =  1 / PC13          <-- eps CANCELS")
    W("")
    W("  The precursor:product ratio is 1/PC13 for EVERY branch, exactly,")
    W("  independent of the escape ratio.  Numerically, at PC13 = %.3f:" % pc)
    W("")
    W("  %-32s %14s %14s" % ("branch", "prec/prod", "1/PC13"))
    for nm, eps in branches:
        r = (M.branch_precursor(eps, pc) / (1.0 - M.branch_loss(eps, pc)))
        W("  %-32s %14.10f %14.10f" % (nm, r, 1.0 / pc))
    W("")
    W("  So the model implies a complete and slightly bleak measurement theory:")
    W("")
    W("    precursor + product  =  S            measures SYNTHESIS;")
    W("                                         exactly blind to PC1/3")
    W("    precursor / product  =  1/PC13       measures the CONVERTASE;")
    W("                                         exactly blind to eps, i.e.")
    W("                                         identical for all five branches")
    W("    product alone        =  S/(1+eps/PC13)   the ONLY quantity that")
    W("                                         carries branch information")
    W("")
    W("  The two easy assays are each blind to precisely the thing the other")
    W("  measures, and neither can tell you which branch is failing.  The one")
    W("  informative quantity — absolute bioactive product — is the hardest to")
    W("  measure in a living patient, especially for a neuropeptide released")
    W("  centrally.  This is a testable claim about assay design, not a")
    W("  rhetorical flourish: it says a hyperproinsulinaemia measurement and a")
    W("  pro-oxytocin:oxytocin measurement should give the SAME number in the")
    W("  same patient, and that number grades the convertase, not the")
    W("  phenotype.")


# ---------------------------------------------------------------------------
# 2.  THE CONSERVATION RESULT
# ---------------------------------------------------------------------------
def section2():
    W("")
    W("[2] WHY A CROSS-REACTING ASSAY IS EXACTLY BLIND TO THE LESION")
    W(hr())
    W("  At steady state, mass in = mass out:")
    W("      S = kc*PC13*P  +  d*P  =  (product flux) + (escape flux)")
    W("  So the SUM of the two fluxes equals S for every PC13.  An assay that")
    W("  cannot separate precursor from product measures that sum.")
    W("")
    W("  %-24s %10s %10s %10s %10s"
      % ("branch", "PC13", "product", "escape", "SUM"))
    for nm, eps in (("pro-oxytocin", P["EPSPOXT"]), ("POMC", P["EPSPOMC"]),
                    ("proinsulin", P["EPSPINS"])):
        for pc in (1.0, 0.4, 0.15):
            prod = 1.0 / (1.0 + eps / pc)
            esc = 1.0 - prod
            W("  %-24s %10.2f %10.6f %10.6f %10.6f"
              % (nm if pc == 1.0 else "", pc, prod, esc, prod + esc))
    W("")
    W("  The SUM column is 1.000000 everywhere.  This is not a numerical")
    W("  coincidence and it is not parameter-dependent: it is conservation.")
    W("  It predicts that studies measuring total immunoreactive oxytocin in")
    W("  plasma or CSF in PWS should give inconsistent, near-null results —")
    W("  which is what the literature reports — while bioactive, meal-entrained")
    W("  central release is markedly reduced.  Any trial that uses a")
    W("  cross-reacting oxytocin assay as a target-engagement biomarker is")
    W("  measuring a quantity the model says CANNOT respond.")


# ---------------------------------------------------------------------------
# 3.  HARMONIC-MEAN SENSITIVITY
# ---------------------------------------------------------------------------
def arms_pws():
    """The five satiety arms of a food-secure PWS 12-year-old on GH."""
    return dict(x1=0.2965, x2=0.7073, x3=0.6063, x4=1.0105, x5=1.1593)


def arms_ctrl():
    return dict(x1=0.9985, x2=1.0110, x3=1.0110, x4=1.0110, x5=1.0048)


def sat_of(x):
    w = (P["W1"], P["W2"], P["W3"], P["W4"], P["W5"])
    return 1.0 / sum(wi / x["x%d" % (i + 1)] for i, wi in enumerate(w))


def section3():
    W("")
    W("[3] HARMONIC-MEAN SENSITIVITY — WHICH ARM IS WORTH REPAIRING")
    W(hr())
    W("      SAT = 1 / sum_i (w_i / x_i)")
    W("      dSAT/dx_i = SAT^2 * w_i / x_i^2")
    W("")
    W("  The sensitivity to an arm scales as w_i / x_i^2, so it grows")
    W("  QUADRATICALLY as an arm gets worse.  A harmonic mean does not just")
    W("  refuse to be rescued by its healthy terms — it actively concentrates")
    W("  all available leverage on its worst term.")
    W("")
    x = arms_pws()
    S = sat_of(x)
    W("  PWS @12y, SAT = %.4f" % S)
    W("  %-6s %8s %8s %12s %12s"
      % ("arm", "w_i", "x_i", "dSAT/dx_i", "share of 1/SAT"))
    w = (P["W1"], P["W2"], P["W3"], P["W4"], P["W5"])
    tot = 1.0 / S
    for i, wi in enumerate(w):
        xi = x["x%d" % (i + 1)]
        W("  %-6s %8.2f %8.4f %12.5f %11.1f%%"
          % ("x%d" % (i + 1), wi, xi, S * S * wi / (xi * xi),
             100.0 * (wi / xi) / tot))
    W("")
    W("  the oxytocin arm carries %.1f%% of the resistance to satiety while"
      % (100.0 * (P["W1"] / x["x1"]) / tot))
    W("  holding only %.0f%% of the weight.  That asymmetry is why every other"
      % (100.0 * P["W1"]))
    W("  arm can be repaired for almost nothing:")
    W("")
    W("  %-42s %9s %9s" % ("repair", "SAT", "dSAT"))
    for lab, key, val in (("x2 (vagal) -> 1.0", "x2", 1.0),
                          ("x3 (PYY) -> 1.0", "x3", 1.0),
                          ("x4 (GLP-1) -> 2.0", "x4", 2.0),
                          ("x5 (leptin/insulin) -> 1.5", "x5", 1.5)):
        y = dict(x); y[key] = val
        W("  %-42s %9.4f %+9.4f" % (lab, sat_of(y), sat_of(y) - S))
    y = dict(x); y.update(x2=1.0, x3=1.0, x4=1.0, x5=1.0)
    W("  %-42s %9.4f %+9.4f" % ("ALL FOUR non-oxytocin arms -> 1.0",
                                sat_of(y), sat_of(y) - S))
    for g in (1.5, 2.0, 2.4, 3.0, 3.37):
        y = dict(x); y["x1"] = x["x1"] * g
        W("  %-42s %9.4f %+9.4f" % ("x1 (oxytocin) x %.2f" % g,
                                    sat_of(y), sat_of(y) - S))
    W("")
    W("  control reference SAT = %.4f" % sat_of(arms_ctrl()))


# ---------------------------------------------------------------------------
# 4.  THE MELANOCORTIN RELAY CEILING
# ---------------------------------------------------------------------------
def section4():
    W("")
    W("[4] THE MELANOCORTIN RELAY CEILING")
    W(hr())
    W("      x1 = OXT * relay(MC) * (1 + E_OXTR)")
    W("      relay(MC) = (MC/(K+MC)) * (K+1),   K = KREL = %.2f" % P["KREL"])
    W("")
    W("  relay(1) = 1 by construction, and relay(inf) = K+1 = %.2f."
      % (P["KREL"] + 1.0))
    W("  So MC4R agonism, no matter how potent, cannot multiply the oxytocin")
    W("  arm by more than %.2f — and it starts from PWS's alpha-MSH tone of"
      % (P["KREL"] + 1.0))
    W("  0.609, where relay = %.4f."
      % ((0.609 / (P["KREL"] + 0.609)) * (P["KREL"] + 1.0)))
    W("")
    W("  %-30s %10s %10s %10s" % ("melanocortin input MC", "relay(MC)",
                                  "x1", "SAT"))
    x = arms_pws()
    OXTtone = 0.3281
    for mc in (0.609, 1.0, 2.0, 3.14, 5.0, 10.0, 1e6):
        r = (mc / (P["KREL"] + mc)) * (P["KREL"] + 1.0)
        y = dict(x); y["x1"] = OXTtone * r
        lab = "MC = %.3g" % mc if mc < 1e5 else "MC -> infinity (the CEILING)"
        W("  %-30s %10.4f %10.4f %10.4f" % (lab, r, y["x1"], sat_of(y)))
    W("")
    W("  Setmelanotide 3 mg/d reaches MC = 3.14 in the model, i.e. %.0f%% of the"
      % (100.0 * ((3.14 / (P["KREL"] + 3.14)) * (P["KREL"] + 1.0))
         / (P["KREL"] + 1.0)))
    W("  ceiling already.  There is essentially nothing left to gain from a")
    W("  more potent MC4R agonist, and that is a STRUCTURAL statement:")
    W("  the drug enters ABOVE the block.  It is also why setmelanotide works")
    W("  in POMC and LEPR deficiency (lesion above MC4R) and not in PWS")
    W("  (lesion below it) — the same molecule, two different topologies.")
    W("")
    W("  Contrast the oxytocin-receptor route, which enters BELOW the block")
    W("  and is therefore not bounded by relay():")
    for e in (0.0, 0.46, 0.49, 1.0, 1.4):
        y = dict(x); y["x1"] = OXTtone * ((0.609 / (P["KREL"] + 0.609))
                                         * (P["KREL"] + 1.0)) * (1.0 + e)
        W("    E_OXTR = %.2f  ->  x1 = %.4f, SAT = %.4f" % (e, y["x1"], sat_of(y)))


# ---------------------------------------------------------------------------
# 5.  GHRELIN ELASTICITY
# ---------------------------------------------------------------------------
def section5():
    W("")
    W("[5] GHRELIN-ARM ELASTICITY")
    W(hr())
    K = P["KAG50"]
    W("  occupancy(AG) = AG/(K+AG),  K = %.0f pg/mL" % K)
    W("  elasticity     d ln(occ) / d ln(AG) = K/(K+AG)")
    W("")
    W("  %10s %11s %11s %11s" % ("AG pg/mL", "occupancy", "arm/control",
                                 "elasticity"))
    o0 = 350.0 / (K + 350.0)
    for ag in (150, 250, 350, 500, 700, 870, 1200, 1800):
        occ = ag / (K + ag)
        W("  %10d %11.4f %11.4f %11.4f" % (ag, occ, occ / o0, K / (K + ag)))
    W("")
    W("  At the model's PWS concentration (870 pg/mL) the elasticity is only")
    W("  %.3f: a 50%% fall in acyl-ghrelin moves the receptor arm by %.1f%%."
      % (K / (K + 870.0), 100.0 * ((435.0 / (K + 435.0)) / (870.0 / (K + 870.0)) - 1.0)))
    W("  With the arm's weight in the drive at K_ghr = %.2f, that is a %.1f%%"
      % (P["KGHRD"], 100.0 * P["KGHRD"]
         * abs((435.0 / (K + 435.0)) / (870.0 / (K + 870.0)) - 1.0)))
    W("  change in orexigenic drive.  This is the model's account of three")
    W("  independent negative results (octreotide, livoletide, and")
    W("  ghrelin-lowering in general) with ONE parameter.")
    W("")
    W("  HONESTY NOTE.  K_ghr = %.2f is CALIBRATED to those negative trials,"
      % P["KGHRD"])
    W("  not derived.  What is derived is that a single small value reconciles")
    W("  all three of them simultaneously with a 2-fold hyperghrelinaemia —")
    W("  i.e. the trials are mutually consistent under one number, and any")
    W("  future ghrelin-directed agent is predicted to fail at the same size.")


# ---------------------------------------------------------------------------
# 6.  THE BIFURCATION DIAGRAM
# ---------------------------------------------------------------------------
def seek_rhs(S, DRVe, SAT, selfg):
    sh = S ** P["NSELF"]
    return (P["KSON"] * DRVe * (1.0 - S)
            + selfg * sh / (P["KSHALF"] ** P["NSELF"] + sh) * (1.0 - S)
            - P["KSOFF"] * S * SAT)


def fps(DRVe, SAT, selfg, n=4001):
    out = []
    prev = seek_rhs(0.0, DRVe, SAT, selfg)
    for i in range(1, n):
        S = i / (n - 1.0)
        cur = seek_rhs(S, DRVe, SAT, selfg)
        if (prev < 0.0) != (cur < 0.0):
            lo, hi = (i - 1) / (n - 1.0), S
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                if ((seek_rhs(lo, DRVe, SAT, selfg) < 0.0)
                        != (seek_rhs(mid, DRVe, SAT, selfg) < 0.0)):
                    hi = mid
                else:
                    lo = mid
            Sf = 0.5 * (lo + hi)
            out.append((Sf, "stable" if seek_rhs(Sf + 1e-4, DRVe, SAT, selfg) < 0.0
                        else "SADDLE"))
        prev = cur
    return out


def section6():
    W("")
    W("[6] THE FOOD-SEEKING BIFURCATION DIAGRAM")
    W(hr())
    W("      dS/dt = k_on*DRVe*(1-S) + G*S^4/(K^4+S^4)*(1-S) - k_off*S*SAT")
    W("      k_on = %.3f   G = K_self*(RFL + (1-RFL)*REINF)   K = %.3f"
      % (P["KSON"], P["KSHALF"]))
    W("      k_off = %.2f   K_self = %.2f   RFL = %.2f"
      % (P["KSOFF"], P["KSELF"], P["RFL"]))
    W("")
    DRVi, SAT, G = 2.0774, 0.4565, 3.0677
    W("  PWS @12y on management: DRVi = %.4f, SAT = %.4f, G = %.4f"
      % (DRVi, SAT, G))
    W("")
    W("  (a) sweep the food environment CUE at fixed SAT")
    W("  %8s %9s   %s" % ("CUE", "DRVe", "fixed points"))
    for cue in (0.05, 0.10, 0.20, 0.28, 0.40, 0.50, 0.55, 0.58, 0.60, 0.80, 1.00):
        f = fps(DRVi * cue, SAT, G)
        W("  %8.2f %9.4f   %s"
          % (cue, DRVi * cue,
             "  ".join("%.4f(%s)" % (a, b) for a, b in f)))
    lo, hi = 0.02, 3.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if len(fps(DRVi * mid, SAT, G)) >= 3:
            lo = mid
        else:
            hi = mid
    cuestar = 0.5 * (lo + hi)
    W("")
    W("  LOWER saddle-node (the low state is annihilated) at CUE* = %.4f"
      % cuestar)
    W("")
    W("  (b) sweep satiety at fixed CUE = 0.28 — the DRUG axis.")
    W("      Two reinforcement states matter, and they give two answers:")
    W("        G = %.4f  a patient on the LOW branch (well managed, food"
      % G)
    W("                    mostly matches the ask, so reinforcement is high)")
    W("        G = %.4f  a LATCHED patient under the same management (the ask"
      % 2.61)
    W("                    is large and unmet, so reinforcement is LOWER)")
    W("")
    for Gl, lab in ((G, "low-branch patient"), (2.6100, "latched patient")):
        W("  --- %s, G = %.4f ---" % (lab, Gl))
        W("  %8s   %s" % ("SAT", "fixed points"))
        for sat in (0.40, 0.4565, 0.55, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00):
            f = fps(DRVi * 0.28, sat, Gl)
            W("  %8.4f   %s"
              % (sat, "  ".join("%.4f(%s)" % (a, b) for a, b in f)))
        lo, hi = 0.20, 4.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if len(fps(DRVi * 0.28, mid, Gl)) >= 2:
                lo = mid
            else:
                hi = mid
        satcrit = 0.5 * (lo + hi)
        x = arms_pws()
        others = (P["W2"] / x["x2"] + P["W3"] / x["x3"] + P["W4"] / x["x4"]
                  + P["W5"] / x["x5"])
        need = 1.0 / satcrit - others
        W("  upper saddle-node at SAT_c = %.4f  (%.2f x achieved)"
          % (satcrit, satcrit / SAT))
        if need > 0:
            gain = (P["W1"] / need) / x["x1"]
            W("  required oxytocin-arm gain %.2f x;  carbetocin gives 1.49 x;"
              % gain)
            W("  SHORTFALL %.2f x" % (gain / 1.4924))
        else:
            W("  the oxytocin arm alone cannot reach it (1/SAT >= %.4f)" % others)
        W("")
    W("  Note the direction of that difference: the LATCHED patient is the")
    W("  EASIER target, because being latched under management means being")
    W("  chronically unrewarded, which lowers the self-sensitization gain that")
    W("  holds the state up.  The model therefore predicts that the patients")
    W("  most likely to switch branches on an OXTR agonist are the ones with")
    W("  the highest baseline HQ-CT in a food-secure setting — the same")
    W("  stratum the DCCR subgroup analysis picked out, arrived at from a")
    W("  different direction.")
    W("")
    W("  (c) the cusp — the two loci together")
    W("  %8s %12s %12s" % ("SAT", "CUE* (lower)", "SAT_c reached?"))
    for sat in (0.35, 0.4565, 0.60, 0.75, 0.90):
        lo, hi = 0.01, 3.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if len(fps(DRVi * mid, sat, G)) >= 3:
                lo = mid
            else:
                hi = mid
        n_at_28 = len(fps(DRVi * 0.28, sat, G))
        W("  %8.4f %12.4f %12s"
          % (sat, 0.5 * (lo + hi), "YES (monostable low)" if n_at_28 == 1
             and fps(DRVi * 0.28, sat, G)[0][0] < 0.25 else "no"))
    W("")
    W("  raising satiety WIDENS the range of food environments in which a low")
    W("  state exists.  Drug and environment are therefore not additive")
    W("  interventions on one scale — they move two different boundaries of")
    W("  the same cusp, which is why their combination is not predictable from")
    W("  either alone.")


# ---------------------------------------------------------------------------
# 7.  CARBETOCIN
# ---------------------------------------------------------------------------
def section7():
    W("")
    W("[7] CARBETOCIN — THE ANALYTIC OPTIMUM AND THE WIDTH OF THE WINDOW")
    W(hr())
    E1, K1 = P["ECBMAX"], P["ECB50"]
    E2, K2 = P["EV1AMX"], P["EV1A50"]
    W("      E(C) = E1*C/(K1+C) - E2*C/(K2+C)")
    W("      E1 = %.2f  K1 = %.2f ng/mL   (OXTR agonism)" % (E1, K1))
    W("      E2 = %.2f  K2 = %.2f ng/mL   (V1a cross-activation)" % (E2, K2))
    W("")
    W("      E'(C) = E1*K1/(K1+C)^2 - E2*K2/(K2+C)^2 = 0")
    W("      =>  sqrt(E1K1)/(K1+C) = sqrt(E2K2)/(K2+C)")
    W("      =>  C* = [K1*sqrt(E2K2) - K2*sqrt(E1K1)] / [sqrt(E1K1) - sqrt(E2K2)]")
    a1, a2 = math.sqrt(E1 * K1), math.sqrt(E2 * K2)
    Cs = (K1 * a2 - K2 * a1) / (a1 - a2)

    def E(C):
        return E1 * C / (K1 + C) - E2 * C / (K2 + C)

    W("")
    W("      sqrt(E1K1) = %.6f    sqrt(E2K2) = %.6f" % (a1, a2))
    W("      C*         = %.6f ng/mL" % Cs)
    W("      E(C*)      = %.6f" % E(Cs))
    W("      E(0+) slope = E1/K1 - E2/K2 = %.4f  (positive: it starts to rise)"
      % (E1 / K1 - E2 / K2))
    W("      E(infinity) = E1 - E2 = %.4f  (positive but %.0f%% of the peak)"
      % (E1 - E2, 100.0 * (E1 - E2) / E(Cs)))
    W("")
    W("  A biphasic curve exists at all only because E1/K1 > E2/K2 while")
    W("  E1 - E2 < E(C*): the agonist is the higher-affinity action and the")
    W("  counter-action is the higher-capacity one.  Both conditions hold here,")
    W("  and neither was imposed for that purpose.")
    W("")
    W("  numerical check of the closed form (finite differences):")
    for d in (1e-3, 1e-4, 1e-5):
        W("      dE/dC at C* with h = %.0e : %+.3e" % (d, (E(Cs + d) - E(Cs - d)) / (2 * d)))
    W("")
    W("  dose mapping.  Cavg over a dosing interval = n*D*F/(V*ke)")
    W("      F = %.3f  V = %.1f L  ke = %.1f /d  (n = 3 doses/d)"
      % (P["FCB"], P["VCB"], P["KECB"]))
    W("")
    W("  %10s %11s %10s %10s %10s"
      % ("mg TID", "Cavg ng/mL", "E(C)", "E/E(C*)", "C/C*"))
    for dose in (0.4, 0.8, 1.2, 1.6, 2.4, 3.2, 4.8, 6.4, 9.6, 14.4, 19.2):
        cavg = 3.0 * dose * P["FCB"] / (P["VCB"] * P["KECB"]) * 1000.0
        W("  %10.1f %11.4f %10.4f %10.4f %10.3f"
          % (dose, cavg, E(cavg), E(cavg) / E(Cs), cavg / Cs))
    W("")
    lo, hi = 1e-4, Cs
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if E(mid) < 0.9 * E(Cs):
            lo = mid
        else:
            hi = mid
    c_lo = 0.5 * (lo + hi)
    lo, hi = Cs, 1e4
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if E(mid) > 0.9 * E(Cs):
            lo = mid
        else:
            hi = mid
    c_hi = 0.5 * (lo + hi)
    W("  the 90%%-of-peak window is C in [%.3f, %.3f] ng/mL, i.e. a %.1f-fold"
      % (c_lo, c_hi, c_hi / c_lo))
    W("  concentration range — narrow, and BOUNDED ABOVE.  In dose units that")
    W("  is %.2f to %.2f mg TID."
      % (c_lo * P["VCB"] * P["KECB"] / (3.0 * P["FCB"] * 1000.0),
         c_hi * P["VCB"] * P["KECB"] / (3.0 * P["FCB"] * 1000.0)))
    W("  9.6 mg TID sits at %.0f%% of peak effect — outside the window, on the"
      % (100.0 * E(3.0 * 9.6 * P["FCB"] / (P["VCB"] * P["KECB"]) * 1000.0) / E(Cs)))
    W("  descending limb.  The model therefore reproduces the reported dose")
    W("  ordering (3.2 mg > 9.6 mg) with no dose-specific parameter, and")
    W("  predicts that a dose-ranging study below 3.2 mg should find the")
    W("  effect PRESERVED, not lost.")


# ---------------------------------------------------------------------------
# 8.  DCCR
# ---------------------------------------------------------------------------
def section8():
    W("")
    W("[8] DCCR — A THERAPEUTIC INDEX FIXED BY TWO EC50s")
    W(hr())
    W("  Both actions are the SAME pharmacology (K-ATP opening) in two")
    W("  tissues, so the model cannot give them independent dose-responses:")
    W("      wanted:  AgRP neuron   Emax %.2f  EC50 %.1f ug/mL"
      % (P["EDZMAX"], P["EDZ50"]))
    W("      unwanted: beta cell    Emax %.2f  EC50 %.1f ug/mL"
      % (P["EDZIMX"], P["EDZI50"]))
    W("")
    W("  %10s %10s %10s %10s %10s"
      % ("C ug/mL", "E(AgRP)", "E(bcell)", "ratio", "mg/kg/d"))
    for c in (2.0, 5.0, 10.0, 12.87, 20.0, 30.0, 50.0, 100.0):
        ea = P["EDZMAX"] * c / (P["EDZ50"] + c)
        eb = P["EDZIMX"] * c / (P["EDZI50"] + c)
        dose = c * P["VDZ"] * P["KEDZ"]
        W("  %10.2f %10.4f %10.4f %10.4f %10.2f" % (c, ea, eb, ea / eb, dose))
    W("")
    W("  the ratio E(AgRP)/E(bcell) is monotone and bounded:")
    W("      C -> 0    ratio -> (Emax_a/EC50_a)/(Emax_b/EC50_b) = %.4f"
      % ((P["EDZMAX"] / P["EDZ50"]) / (P["EDZIMX"] / P["EDZI50"])))
    W("      C -> inf  ratio -> Emax_a/Emax_b = %.4f"
      % (P["EDZMAX"] / P["EDZIMX"]))
    W("  so there is no dose at which the wanted action outruns the unwanted")
    W("  one by more than %.2f-fold.  Escalating the dose buys efficacy and"
      % max((P["EDZMAX"] / P["EDZ50"]) / (P["EDZIMX"] / P["EDZI50"]),
            P["EDZMAX"] / P["EDZIMX"]))
    W("  hyperglycaemia in a FIXED proportion.  Only a Kir6.2/SUR1 subtype- or")
    W("  tissue-selective opener changes that ratio, and that is a medicinal-")
    W("  chemistry target the model states quantitatively.")


# ---------------------------------------------------------------------------
# 9.  THE TWO AIRWAY CLOCKS
# ---------------------------------------------------------------------------
def section9():
    W("")
    W("[9] THE TWO AIRWAY CLOCKS — PEAK TIME IN CLOSED FORM")
    W(hr())
    W("  Strip the airway block to its skeleton: a step change in IGF-1 drives")
    W("  an obstructive term up with rate k_f and a protective term up with")
    W("  rate k_s, and AHI is their difference.")
    W("")
    W("      AHI(t) - AHI(0) = A*(1 - exp(-k_f t)) - B*(1 - exp(-k_s t))")
    W("      d/dt = 0  =>  t_peak = ln( (A k_f)/(B k_s) ) / (k_f - k_s)")
    W("")
    kf = P["KLYU"] * 1.0 * (P["LMAX"] - 1.3) + P["KLYD"] * (1.0 + 0.0)
    ks = 1.0 / P["TAURMS"]
    W("  from the parameter block:")
    W("      k_f (lymphoid, up-phase)  ~ %.4f /d   (tau %.0f d)" % (kf, 1.0 / kf))
    W("      k_s (muscle strength)     = 1/TAURMS = %.4f /d   (tau %.0f d)"
      % (ks, P["TAURMS"]))
    W("      adaptation (involution)   = 1/TAUAD  = %.4f /d   (tau %.0f d)"
      % (1.0 / P["TAUAD"], P["TAUAD"]))
    W("")
    W("  %8s %8s %10s %10s" % ("A", "B", "t_peak (d)", "t_peak (wk)"))
    for A, B in ((3.0, 2.0), (3.0, 3.0), (3.0, 5.0), (2.0, 6.0)):
        num = math.log((A * kf) / (B * ks))
        tp = num / (kf - ks)
        W("  %8.1f %8.1f %10.1f %10.1f" % (A, B, tp, tp / 7.0))
    W("")
    W("  Over the whole plausible range of amplitudes the peak sits between")
    W("  weeks 3 and 8, because it is set by the RATIO of two rate constants")
    W("  under a logarithm — the least amplitude-sensitive quantity the model")
    W("  produces.  This is why 'polysomnography at 6-8 weeks' is a robust")
    W("  recommendation and not a convention: the timing barely depends on how")
    W("  big the effect is, only on how fast lymphoid tissue and muscle move.")
    W("")
    W("  The full model (section [7] of pws_reference_output.txt) puts the peak")
    W("  at week 4 with AHI 11.2 -> 13.2 for a typical airway, and 11.2 ->")
    W("  22.1 for a vulnerable one, both returning below baseline by 1 year.")
    W("  Airway vulnerability multiplies the OBSTRUCTIVE terms only, so it")
    W("  scales the window and not the plateau — which is exactly the shape of")
    W("  the reported early-GH mortality signal.")


# ---------------------------------------------------------------------------
# 10.  ENERGY SCAFFOLD CONSISTENCY
# ---------------------------------------------------------------------------
def section10():
    W("")
    W("[10] ENERGY-BALANCE CONSISTENCY OF THE REFERENCE SCAFFOLD")
    W(hr())
    W("  The reference individual must be in energy balance at every age, or")
    W("  the control drifts and every PWS comparison is contaminated.  Check:")
    W("")
    W("  %6s %9s %9s %9s %9s %9s %9s"
      % ("age", "BWnorm", "Schofield", "EIREQ", "PAL", "kcal/cm", "LEPref"))
    for a in (0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 25, 40):
        bw = M.LFMNORM(a) + M.FMNORM(a)
        ree = M.REEREF(a)
        eir = M.EIREQREF(a)
        W("  %6.1f %9.1f %9.0f %9.0f %9.3f %9.1f %9.2f"
          % (a, bw, ree, eir, eir / ree, eir / M.HTNORM(a), M.LEPREF(a)))
    W("")
    W("  PAL is constant at %.3f by construction: EIREQ = REE*(1+PALX*ACTREF)"
      % (M.EIREQREF(10.0) / M.REEREF(10.0)))
    W("  /(1-FDIT) with ACTREF = 1.  The kcal/cm column is the quantity PWS")
    W("  clinics actually prescribe against; the model's reference child sits")
    W("  at %.0f-%.0f kcal/cm, and standard PWS management (about 10-11"
      % (M.EIREQREF(6.0) / M.HTNORM(6.0), M.EIREQREF(14.0) / M.HTNORM(14.0)))
    W("  kcal/cm) is therefore %.0f-%.0f%% of the reference requirement —"
      % (100.0 * 10.0 / (M.EIREQREF(6.0) / M.HTNORM(6.0)),
         100.0 * 11.0 / (M.EIREQREF(14.0) / M.HTNORM(14.0))))
    W("  which is the 60-80%% figure the guidelines quote, recovered rather")
    W("  than assumed.")
    W("")
    W("  PWS expenditure, from composition alone (no PWS-specific REE term):")
    W("  %6s %10s %10s %10s %10s"
      % ("age", "LFM/norm", "ACT", "REE/ref", "TEE/req"))
    for a, lr, act in ((4, 0.84, 0.66), (8, 0.84, 0.66), (12, 0.84, 0.68),
                       (18, 0.85, 0.70)):
        ree = P["WREEL"] * lr + P["WREEF"] * 2.4
        tee = ree * (1.0 + P["PALX"] * act) + 0.10 * 0.85 * 1.80
        W("  %6d %10.2f %10.2f %10.3f %10.3f"
          % (a, lr, act, ree, ree * (1.0 + P["PALX"] * act) / 1.80 + 0.085))
    W("  the 60-70%% of predicted TEE reported in PWS falls out of (low lean")
    W("  mass) x (low activity) with no expenditure parameter of its own.")


# ---------------------------------------------------------------------------
# 11.  LOCAL SENSITIVITY (integrates the ODEs)
# ---------------------------------------------------------------------------
def section11():
    W("")
    W("[11] LOCAL PARAMETER SENSITIVITY OF THE PRIMARY ENDPOINTS")
    W(hr())
    W("  Each parameter is moved +/-20% and the model re-run birth -> 12 y on")
    W("  standard management with growth hormone.  Reported: elasticity")
    W("  d ln(endpoint) / d ln(parameter), central difference.")
    W("")
    keys = ["HQ", "PBF", "HT", "IGF1", "AG", "AHI", "SAT", "SEEK"]
    env = M.ENV_MGMT
    ctl = M.CT(env["avail"], env["cue"], env["titr"])

    def run(**kw):
        return M.run_to(les=1.0, age=12.0, gh_from=1.0, dt=0.25,
                        **dict(env, **kw)).out(ctl)

    base = run()
    W("  baseline: " + "  ".join("%s %.4g" % (k, base[k]) for k in keys))
    W("")
    W("  %-10s" % "parameter" + "".join("%9s" % k for k in keys))
    for par in ("DPC13", "FOXTN", "EPSPOXT", "EPSPOMC", "KREL", "W1",
                "KINHO", "KGHRD", "FGHRC", "KAG50", "KSELF", "KSOFF",
                "KSON", "WLIGF0", "KGIGF", "KTITR", "KLYU", "KAHIF",
                "FPOTGH", "KLEPEI"):
        v0 = P[par]
        try:
            hi = run(**{par: v0 * 1.2})
            lo = run(**{par: v0 * 0.8})
        except Exception as exc:                       # pragma: no cover
            W("  %-10s  (failed: %s)" % (par, exc))
            continue
        row = "  %-10s" % par
        for k in keys:
            if abs(base[k]) < 1e-9:
                row += "%9s" % "-"
            else:
                el = ((hi[k] - lo[k]) / base[k]) / (0.4)
                row += "%9.3f" % el
        W(row)
    W("")
    W("  READ THIS TABLE AS A MAP OF WHAT THE MODEL COMMITS TO.")
    W("  HQ-CT is dominated by the oxytocin branch (FOXTN, EPSPOXT, W1, KINHO)")
    W("  and by the bistability constants, and is nearly INSENSITIVE to the")
    W("  ghrelin parameters (KGHRD, FGHRC, KAG50) — the same orthogonality the")
    W("  drug panel shows, now as a derivative rather than a simulation.")
    W("  %fat and height load onto the somatotropic and titration parameters")
    W("  and not onto the satiety ones.  A model in which one lesion produced")
    W("  one undifferentiated 'severity' could not do this.")


def main():
    section1(); section2(); section3(); section4(); section5()
    section6(); section7(); section8(); section9(); section10()
    if "--sens" in sys.argv:
        section11()
    else:
        W("")
        W("[11] LOCAL PARAMETER SENSITIVITY — skipped (re-run with --sens)")
    W("")
    W("=" * 78)
    W("end of calibration output")
    W("=" * 78)
    print("\n".join(OUT))


if __name__ == "__main__":
    main()
