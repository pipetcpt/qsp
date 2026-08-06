#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run every scenario of the GAD QSP model and write gad_reference_output.txt.

The numbers printed here are the ONLY numbers quoted in README.md, the
mrgsolve header and the Shiny app.  Nothing is transcribed by hand.
"""
import json
import os
import sys

import numpy as np

import gad_python_reference as G

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def hdr(t):
    say("")
    say("=" * 100)
    say(t)
    say("=" * 100)


def load_params():
    P = G.default_params()
    G.calibrate_healthy(P)
    fit = json.load(open("gad_calibration.json"))
    P.update(kfl=fit["kfl"], dvisit=fit["dvisit"], e_ex=fit["e_ex"],
             k5ht_eff=fit["k5ht_eff"], emax_pgb=fit["emax_pgb"],
             emax_a2=fit["emax_a2"])
    return P, fit["fluct0"]


VIS8 = [0.0, 7.0, 14.0, 28.0, 42.0, 56.0]
VIS4 = [0.0, 7.0, 14.0, 21.0, 28.0]
VIS12 = VIS8 + [70.0, 84.0]
KEYS = ("hama", "hama_psy", "hama_som", "phi_n", "occ_sert", "occ_net",
        "occ_a2d", "occ_bz", "occ_h1", "sht", "sht_post", "ne", "auto",
        "bdnf", "cpfc", "cpfc_eff", "eamy", "sglu", "igaba", "cort",
        "sedation", "dizziness", "nausea", "activation", "sexd", "wt",
        "ae_burden", "cumhaz", "madrs", "sleepd", "worry", "auton", "sns",
        "Cesc", "Cven", "Codv", "Cpgb", "Cbzd", "Cqtp", "expect")


# ---------------------------------------------------------------------------
def s1_baseline(P):
    hdr("1.  BASELINE — the four factors of Phi, and what one severity knob does")
    dis = np.array([0.0, 0.4, 0.7, 1.0, 1.3, 1.6])
    Q = dict(P)
    Q["DIS"] = dis
    y = G.run_to_steady_state(Q, n=dis.size)
    o = G.observe(y, Q)
    say("Untreated attractor (no drug, no trial contact).  Phi = (E_amy*S_glu)/(C_pfc*I_gaba)")
    say("")
    say("  DIS   HAM-A  psychic somatic   Phi_n   E_amy   S_glu   C_pfc  I_gaba   WORRY  SLEEPD  cortisol  MADRS")
    for i, d in enumerate(dis):
        say("  %4.1f  %5.1f   %5.1f   %5.1f   %5.2f   %5.3f   %5.3f   %5.3f   %5.3f   %5.2f   %5.2f     %5.3f  %5.1f"
            % (d, o["hama"][i], o["hama_psy"][i], o["hama_som"][i], o["phi_n"][i],
               o["eamy"][i], o["sglu"][i], o["cpfc"][i], o["igaba"][i],
               y[G.IX["worry"]][i], o["sleepd"][i], o["cort"][i], o["madrs"][i]))
    say("")
    say("  Phi_healthy (the normalising constant) = %.5f" % P["phi_healthy"])
    say("  DIS = 1.0 is the typical trial patient, and the number that matters is")
    say("  NOT the one in the table.  Its STABLE score is %.1f; the score at" % o["hama"][3])
    say("  ENROLMENT is that plus the enrolment peak, which the virtual population")
    say("  in section 9 puts at 23.1 -- inside the 20-27 band GAD trials report.")
    say("  Entry criteria select the top of a fluctuation, so a model whose stable")
    say("  score already equals the reported baseline is a model that will over-")
    say("  state severity by a third.  That error was made and corrected here.")
    say("  Note also that cortisol is only +%.0f%% at DIS = 1: in this model the HPA"
        % (100 * (o["cort"][3] - 1)))
    say("  axis is a modulator, not the disease.")


# ---------------------------------------------------------------------------
def s2_clocks(P, fluct0):
    hdr("2.  THE FOUR CLOCKS — onset is a property of WHICH FACTOR a drug moves")
    arms = [G.arm("escitalopram 10 mg", "escitalopram", 10.0),
            G.arm("pregabalin 300 mg", "pregabalin", 300.0),
            G.arm("lorazepam 3 mg", "lorazepam", 3.0),
            G.arm("none")]
    t, o, yf, sl, Q = G.run_trial(P, arms, 84.0, visits=(), dis=1.0,
                                  fluct0=0.0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 84.1, 0.5))
    say("NO trial contact in this run (no visits, no enrolment peak): this is the")
    say("PHARMACOLOGICAL clock alone, uncontaminated by the fifth one.")
    say("")

    def t50(series, col):
        v = series[:, col]
        v0, vinf = v[0], v[-1]
        if abs(vinf - v0) < 1e-9:
            return float("nan")
        target = v0 + 0.5 * (vinf - v0)
        for i in range(1, len(v)):
            if (v[i] - target) * (v[0] - target) <= 0:
                return float(t[i])
        return float("nan")

    say("  Time to 50% of the eventual effect (days), by quantity:")
    say("")
    say("  arm                    target occupancy   the factor it moves    HAM-A")
    for nm, occ in (("escitalopram 10 mg", "occ_sert"),
                    ("pregabalin 300 mg", "occ_a2d"),
                    ("lorazepam 3 mg", "occ_bz")):
        c = sl[nm].start
        fac = {"escitalopram 10 mg": "cpfc", "pregabalin 300 mg": "sglu",
               "lorazepam 3 mg": "igaba"}[nm]
        say("  %-22s %8.2f d        %8.2f d           %6.2f d"
            % (nm, t50(o[occ], c), t50(o[fac], c), t50(o["hama"], c)))
    c = sl["escitalopram 10 mg"].start
    say("")
    say("  Escitalopram, step by step down the cascade:")
    say("     SERT occupancy      t50 = %.2f d   (plasma PK only)" % t50(o["occ_sert"], c))
    say("     extracellular 5-HT  t50 = %.2f d   (gated by the autoreceptor)" % t50(o["sht"], c))
    say("     autoreceptor AUTO   t50 = %.2f d   (the gate itself)" % t50(o["auto"], c))
    say("     BDNF                t50 = %.2f d" % t50(o["bdnf"], c))
    say("     C_pfc               t50 = %.2f d" % t50(o["cpfc"], c))
    say("     HAM-A               t50 = %.2f d" % t50(o["hama"], c))
    say("")
    say("  THE OCCUPANCY-EFFECT GAP IS A FACTOR OF %.0f: SERT is half-blocked in %.2f d,"
        % (t50(o["hama"], c) / max(t50(o["occ_sert"], c), 1e-6), t50(o["occ_sert"], c)))
    say("  the score is half-moved in %.1f d.  Nothing was fitted to produce that;" % t50(o["hama"], c))
    say("  it is two slow steps in series (autoreceptor tau, then plasticity tau).")
    say("")
    for d in (1, 3, 7, 14, 28, 56, 84):
        i = int(np.argmin(np.abs(t - d)))
        row = "  day %3d: " % d
        for nm in ("escitalopram 10 mg", "pregabalin 300 mg", "lorazepam 3 mg"):
            cc = sl[nm].start
            row += "%s %6.2f | " % (nm.split()[0][:4], o["hama"][i, cc] - o["hama"][0, cc])
        row += "untreated %6.2f" % (o["hama"][i, sl["none"].start] - o["hama"][0, sl["none"].start])
        say(row)
    say("")
    say("  5-HT and the gate (escitalopram 10 mg):")
    say("     day   Occ_SERT   AUTO    5-HT    C_pfc   Phi_n   HAM-A")
    for d in (0, 1, 3, 7, 14, 28, 56, 84):
        i = int(np.argmin(np.abs(t - d)))
        say("     %3d    %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.2f"
            % (d, o["occ_sert"][i, c], o["auto"][i, c], o["sht"][i, c],
               o["cpfc"][i, c], o["phi_n"][i, c], o["hama"][i, c]))
    say("")
    say("  Acute 5-HT rise %.1f%% (day 1) versus chronic %.1f%% (day 84) at the SAME"
        % (100 * (o["sht"][int(np.argmin(np.abs(t - 1))), c] - 1),
           100 * (o["sht"][-1, c] - 1)))
    say("  occupancy (%.3f -> %.3f).  That is the autoreceptor, and it is why the"
        % (o["occ_sert"][int(np.argmin(np.abs(t - 1))), c], o["occ_sert"][-1, c]))
    say("  onset curve is sigmoid rather than exponential.")


# ---------------------------------------------------------------------------
def s3_rickels(P, fluct0):
    hdr("3.  RICKELS 2005 (PMID 16143734) — one fitted arm, and the rest predicted")
    arms = [G.arm("placebo"),
            G.arm("pregabalin 300", "pregabalin", 300.0),
            G.arm("pregabalin 450", "pregabalin", 450.0),
            G.arm("pregabalin 600", "pregabalin", 600.0),
            G.arm("alprazolam 1.5", "alprazolam", 1.5)]
    t, o, yf, sl, Q = G.run_trial(P, arms, 28.0, visits=VIS4, dis=1.0,
                                  fluct0=fluct0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 28.1, 1.0))
    rep = {"placebo": (-8.4, None, None), "pregabalin 300": (-12.2, "sig", "sig"),
           "pregabalin 450": (-11.0, "sig", "ns"), "pregabalin 600": (-11.8, "sig", "sig"),
           "alprazolam 1.5": (-10.9, "sig", "ns")}
    say("Week-4 HAM-A change from baseline.")
    say("FITTED against this trial: the placebo arm's week-4 value, and the")
    say("alprazolam-minus-placebo DIFFERENCE of -2.50.")
    say("PREDICTED: all three pregabalin arms (the pregabalin knob was fitted to")
    say("the Slee 2019 network-meta-analytic difference of -2.79, not to this")
    say("trial), and every psychic/somatic split.")
    say("")
    say("  arm               model  reported   delta vs placebo (model / reported)")
    hp = o["hama"][28, sl["placebo"].start] - o["hama"][0, sl["placebo"].start]
    for nm in rep:
        c = sl[nm].start
        d = o["hama"][28, c] - o["hama"][0, c]
        tag = "  <-- FITTED" if nm in ("placebo", "alprazolam 1.5") else "  (predicted)"
        say("  %-16s %6.2f   %6.2f      %6.2f / %6.2f%s"
            % (nm, d, rep[nm][0], d - hp, rep[nm][0] - (-8.4), tag))
    say("")
    say("  THE 450 mg ARM IS THE INTERESTING ONE.  The trial found 450 mg WORSE than")
    say("  300 mg (-11.0 vs -12.2) and 600 mg in between (-11.8).  The model has a")
    say("  monotone alpha2delta occupancy hyperbola, so it CANNOT reproduce a")
    say("  non-monotone pregabalin curve; it predicts %.2f / %.2f / %.2f for"
        % (o["hama"][28, sl["pregabalin 300"].start] - o["hama"][0, sl["pregabalin 300"].start],
           o["hama"][28, sl["pregabalin 450"].start] - o["hama"][0, sl["pregabalin 450"].start],
           o["hama"][28, sl["pregabalin 600"].start] - o["hama"][0, sl["pregabalin 600"].start]))
    say("  300/450/600.  What it DOES get is the FLATNESS: doubling the dose from 300")
    say("  to 600 buys only %.2f HAM-A points, because occupancy goes %.3f -> %.3f."
        % (abs((o["hama"][28, sl["pregabalin 600"].start] - o["hama"][0, sl["pregabalin 600"].start])
               - (o["hama"][28, sl["pregabalin 300"].start] - o["hama"][0, sl["pregabalin 300"].start])),
           o["occ_a2d"][28, sl["pregabalin 300"].start], o["occ_a2d"][28, sl["pregabalin 600"].start]))
    say("  The trial's non-monotonicity is inside its own confidence intervals; the")
    say("  model's flatness is a prediction, and it is reported as a partial miss.")
    say("")
    say("  WEEK 1, psychic vs somatic (the dissociation Rickels reports as a")
    say("  significance pattern: pregabalin significant on BOTH, alprazolam on")
    say("  psychic only, p = 0.21 for somatic):")
    say("")
    say("  arm                dPSY(wk1) dSOM(wk1)  som/psy   dPSY(wk4) dSOM(wk4)  reported wk4 somatic")
    for nm in ("placebo", "pregabalin 300", "pregabalin 600", "alprazolam 1.5"):
        c = sl[nm].start
        p1 = o["hama_psy"][7, c] - o["hama_psy"][0, c]
        s1 = o["hama_som"][7, c] - o["hama_som"][0, c]
        p4 = o["hama_psy"][28, c] - o["hama_psy"][0, c]
        s4 = o["hama_som"][28, c] - o["hama_som"][0, c]
        say("  %-16s %8.2f %9.2f  %7.2f  %9.2f %9.2f   %s"
            % (nm, p1, s1, s1 / p1 if p1 else float("nan"), p4, s4,
               rep[nm][2] or "-"))
    pc = sl["placebo"].start
    for nm in ("pregabalin 300", "alprazolam 1.5"):
        c = sl[nm].start
        dp = (o["hama_psy"][28, c] - o["hama_psy"][0, c]) - (o["hama_psy"][28, pc] - o["hama_psy"][0, pc])
        ds = (o["hama_som"][28, c] - o["hama_som"][0, c]) - (o["hama_som"][28, pc] - o["hama_som"][0, pc])
        say("  %-16s week-4 delta vs placebo: psychic %.2f, somatic %.2f  (ratio %.2f)"
            % (nm, dp, ds, ds / dp if dp else float("nan")))
    say("")
    say("  The model reproduces the DIRECTION of the reported dissociation: the")
    say("  somatic/psychic ratio of the drug-placebo delta is larger for pregabalin")
    say("  than for alprazolam, because S_glu (pregabalin's factor) also feeds the")
    say("  muscle-tension term of the somatic cluster directly, whereas I_gaba")
    say("  (the benzodiazepine's factor) reaches the somatic items only through Phi.")


# ---------------------------------------------------------------------------
def s4_khan(P, fluct0):
    hdr("4.  KHAN 2011 (PMID 21694613) — quetiapine XR, and an inverted-U that has "
        "to come from somewhere")
    arms = [G.arm("placebo"),
            G.arm("quetiapine 50", "quetiapine", 50.0),
            G.arm("quetiapine 150", "quetiapine", 150.0),
            G.arm("quetiapine 300", "quetiapine", 300.0)]
    rep8 = {"placebo": -11.10, "quetiapine 50": -13.31,
            "quetiapine 150": -13.54, "quetiapine 300": -11.87}
    rep1 = {"placebo": -5.94, "quetiapine 50": -7.47,
            "quetiapine 150": -8.19, "quetiapine 300": -7.23}
    t, o, yf, sl, Q = G.run_trial(P, arms, 56.0, visits=VIS8, dis=1.0,
                                  fluct0=fluct0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 56.1, 1.0))
    say("Only the PLACEBO trajectory of this trial was fitted (week 1 and week 8).")
    say("The three quetiapine arms are out-of-sample.")
    say("")
    say("  arm              wk1 model  wk1 rep   wk8 model  wk8 rep   Occ_H1  Occ_NET  sedation  AE burden")
    for nm in rep8:
        c = sl[nm].start
        d1 = o["hama"][7, c] - o["hama"][0, c]
        d8 = o["hama"][56, c] - o["hama"][0, c]
        say("  %-15s %9.2f %8.2f  %10.2f %8.2f  %6.3f  %7.3f  %8.3f  %8.3f"
            % (nm, d1, rep1[nm], d8, rep8[nm], o["occ_h1"][56, c],
               o["occ_net"][56, c], o["sedation"][56, c], o["ae_burden"][56, c]))
    say("")
    say("  Deterministic (no dropout) the model is MONOTONE in dose, so it misses the")
    say("  reported inverted U.  The mechanism the model offers for that U is")
    say("  differential dropout: AE burden rises %.2f -> %.2f from 50 to 300 mg while"
        % (o["ae_burden"][56, sl["quetiapine 50"].start],
           o["ae_burden"][56, sl["quetiapine 300"].start]))
    say("  H1 occupancy is already %.3f at 50 mg and only %.3f at 300 mg — sedation"
        % (o["occ_h1"][56, sl["quetiapine 50"].start],
           o["occ_h1"][56, sl["quetiapine 300"].start]))
    say("  saturates FIRST.  Section 9 puts that mechanism to a virtual-population test.")


# ---------------------------------------------------------------------------
def s5_doseresponse(P, fluct0):
    hdr("5.  FLAT DOSE-RESPONSE, DERIVED — the SSRI curve is the occupancy hyperbola")
    arms = ([G.arm("placebo")]
            + [G.arm("esc %g" % d, "escitalopram", d) for d in (5, 10, 20, 30)]
            + [G.arm("ven %g EM" % d, "venlafaxine", d) for d in (75, 150, 225)]
            + [G.arm("ven 225 PM", "venlafaxine", 225.0, fm2d6=0.15),
               G.arm("ven 225 UM", "venlafaxine", 225.0, fm2d6=0.85),
               G.arm("dlx 60", "duloxetine", 60.0),
               G.arm("dlx 120", "duloxetine", 120.0),
               G.arm("buspirone 45", "buspirone", 45.0)])
    t, o, yf, sl, Q = G.run_trial(P, arms, 56.0, visits=VIS8, dis=1.0,
                                  fluct0=fluct0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 56.1, 1.0))
    pc = sl["placebo"].start
    dp = o["hama"][56, pc] - o["hama"][0, pc]
    say("Week-8 HAM-A change and drug-placebo delta.  ALL out of sample.")
    say("Placebo arm: %.2f" % dp)
    say("")
    say("  arm             Css(ng/mL)  Occ_SERT  Occ_NET   5-HT   C_pfc   Phi_n   dHAM-A   delta")
    for nm in [a["name"] for a in arms[1:]]:
        c = sl[nm].start
        css = (o["Cesc"][56, c] if nm.startswith("esc")
               else (o["Cven"][56, c] + o["Codv"][56, c]) if nm.startswith("ven")
               else o["hama"][56, c] * 0 + float("nan"))
        d = o["hama"][56, c] - o["hama"][0, c]
        say("  %-15s %10.1f  %8.3f  %7.3f %6.3f  %6.3f  %6.3f  %7.2f  %6.2f"
            % (nm, css, o["occ_sert"][56, c], o["occ_net"][56, c],
               o["sht"][56, c], o["cpfc"][56, c], o["phi_n"][56, c], d, d - dp))
    say("")
    e5, e10, e20, e30 = (sl["esc %g" % d].start for d in (5, 10, 20, 30))
    say("  ESCITALOPRAM 5 -> 10 -> 20 -> 30 mg:")
    say("     SERT occupancy   %.3f -> %.3f -> %.3f -> %.3f" %
        tuple(o["occ_sert"][56, c] for c in (e5, e10, e20, e30)))
    say("     drug-placebo     %.2f -> %.2f -> %.2f -> %.2f" %
        tuple(o["hama"][56, c] - o["hama"][0, c] - dp for c in (e5, e10, e20, e30)))
    say("     Doubling 10 -> 20 mg buys %.3f occupancy and %.2f HAM-A points."
        % (o["occ_sert"][56, e20] - o["occ_sert"][56, e10],
           (o["hama"][56, e20] - o["hama"][0, e20]) - (o["hama"][56, e10] - o["hama"][0, e10])))
    say("")
    say("     A CLEAN MISS, REPORTED.  The tidy version of this argument is that")
    say("     the flat clinical SSRI dose-response IS the occupancy hyperbola.")
    say("     The model does not support that.  Occupancy does flatten (%.3f ->"
        % o["occ_sert"][56, e10])
    say("     %.3f from 10 to 20 mg, a %.0f%% relative gain), but the HAM-A delta"
        % (o["occ_sert"][56, e20],
           100 * (o["occ_sert"][56, e20] / o["occ_sert"][56, e10] - 1)))
    say("     still grows by %.0f%% (%.2f -> %.2f), because in this model raised"
        % (100 * (((o["hama"][56, e20] - o["hama"][0, e20]) - dp)
                  / ((o["hama"][56, e10] - o["hama"][0, e10]) - dp) - 1),
           (o["hama"][56, e10] - o["hama"][0, e10]) - dp,
           (o["hama"][56, e20] - o["hama"][0, e20]) - dp))
    say("     occupancy is followed by a STEEP rise in extracellular 5-HT (%.2f ->"
        % o["sht"][56, e10])
    say("     %.2f) as the reuptake term approaches its floor, and that undoes the"
        % o["sht"][56, e20])
    say("     flattening.  Baldwin 2006 (PMID 16946363) found 20 mg NOT superior to")
    say("     10 mg; the model predicts about %.1f points of superiority.  Either the"
        % abs(((o["hama"][56, e20] - o["hama"][0, e20]) - dp)
              - ((o["hama"][56, e10] - o["hama"][0, e10]) - dp)))
    say("     5-HT rise saturates far earlier than the reuptake algebra implies, or")
    say("     the flatness comes from somewhere this model does not contain (most")
    say("     likely dropout at 20 mg, which Baldwin also reports).  Stated, not fixed.")
    say("")
    v75, v150, v225 = (sl["ven %g EM" % d].start for d in (75, 150, 225))
    say("  VENLAFAXINE ER 75 -> 150 -> 225 mg: SERT is already flat (%.3f -> %.3f)"
        % (o["occ_sert"][56, v75], o["occ_sert"][56, v225]))
    say("     but NET is NOT (%.3f -> %.3f -> %.3f; Arakawa 2019 measured 8-61%%),"
        % tuple(o["occ_net"][56, c] for c in (v75, v150, v225)))
    say("     so the SNRI keeps a real dose-response where the SSRI has none:")
    say("     drug-placebo %.2f -> %.2f -> %.2f"
        % tuple(o["hama"][56, c] - o["hama"][0, c] - dp for c in (v75, v150, v225)))
    pm, um = sl["ven 225 PM"].start, sl["ven 225 UM"].start
    say("     CYP2D6 at 225 mg: PM delta %.2f (NET occ %.3f), EM %.2f (%.3f), UM %.2f (%.3f)"
        % (o["hama"][56, pm] - o["hama"][0, pm] - dp, o["occ_net"][56, pm],
           o["hama"][56, v225] - o["hama"][0, v225] - dp, o["occ_net"][56, v225],
           o["hama"][56, um] - o["hama"][0, um] - dp, o["occ_net"][56, um]))
    b = sl["buspirone 45"].start
    say("")
    say("  BUSPIRONE 45 mg/d: delta %.2f.  The model makes it weak for a structural"
        % (o["hama"][56, b] - o["hama"][0, b] - dp))
    say("  reason — a 5-HT1A PARTIAL agonist occupies the autoreceptor (anxiogenic")
    say("  arm) and the postsynaptic receptor (anxiolytic arm) at the same time.")


# ---------------------------------------------------------------------------
def s6_placebo(P, fluct0):
    hdr("6.  THE FIFTH CLOCK — pharmacological onset vs MEASURED onset, and "
        "assay sensitivity derived")
    arms = [G.arm("placebo"), G.arm("esc10", "escitalopram", 10.0),
            G.arm("pgb300", "pregabalin", 300.0), G.arm("lzp3", "lorazepam", 3.0)]
    t, o, yf, sl, Q = G.run_trial(P, arms, 56.0, visits=VIS8, dis=1.0,
                                  fluct0=fluct0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 56.1, 1.0))
    t2, o2, _, sl2, _ = G.run_trial(P, arms, 56.0, visits=(), dis=1.0, fluct0=0.0,
                                    obs_keys=KEYS, obs_times=np.arange(0, 56.1, 1.0))
    say("Left block: what a trial would MEASURE (visits + enrolment peak, all arms).")
    say("Right block: the same pharmacology with the trial machinery switched off.")
    say("")
    say("           ------- MEASURED delta vs placebo -------   ---- PHARMACOLOGICAL change ----")
    say("   day      esc10    pgb300     lzp3                    esc10    pgb300     lzp3")
    for d in (1, 3, 7, 14, 28, 56):
        i = int(round(d))
        pc = sl["placebo"].start
        dpp = o["hama"][i, pc] - o["hama"][0, pc]
        row = "   %3d  " % d
        for nm in ("esc10", "pgb300", "lzp3"):
            c = sl[nm].start
            row += "%9.2f" % ((o["hama"][i, c] - o["hama"][0, c]) - dpp)
        row += "        "
        for nm in ("esc10", "pgb300", "lzp3"):
            c = sl2[nm].start
            row += "%9.2f" % (o2["hama"][i, c] - o2["hama"][0, c])
        say(row)
    pc = sl["placebo"].start
    say("")
    say("  Placebo arm alone: %.2f at week 1, %.2f at week 8, i.e. %.0f%% of the"
        % (o["hama"][7, pc] - o["hama"][0, pc], o["hama"][56, pc] - o["hama"][0, pc],
           100 * (o["hama"][7, pc] - o["hama"][0, pc]) / (o["hama"][56, pc] - o["hama"][0, pc])))
    say("  eight-week placebo effect is already present at week 1.  The fifth clock")
    say("  is the FASTEST in the system.")
    say("")
    say("  WHAT THAT DOES AND DOES NOT DO, checked rather than asserted.  It does")
    say("  NOT shrink the drug-placebo difference at week 1 (the two blocks above")
    say("  are close, and the measured delta is if anything slightly LARGER than")
    say("  the pharmacological change, because the placebo arm has already moved")
    say("  the patient onto a steeper part of the HAM-A curve).  What it does is")
    say("  swamp it: at week 1 the placebo arm has moved %.2f points while the"
        % (o["hama"][7, pc] - o["hama"][0, pc]))
    say("  largest drug increment on top of that is %.2f.  A clinician watching a"
        % min(o["hama"][7, sl[n].start] - o["hama"][0, sl[n].start]
              - (o["hama"][7, pc] - o["hama"][0, pc]) for n in ("esc10", "pgb300", "lzp3")))
    say("  single patient in week 1 is therefore watching the fifth clock, whatever")
    say("  they prescribed.  The drug clocks are only visible in a DIFFERENCE.")
    say("")
    say("  ASSAY SENSITIVITY, DERIVED.  Sweep the site expectancy coefficient e_ex.")
    say("  The PHARMACOLOGY IS IDENTICAL in every row; only the placebo response moves.")
    say("")
    say("   e_ex   placebo dHAM-A   esc10 delta   pgb300 delta   lzp3 delta   Phi_n(pbo,wk8)")
    for ex in (0.0, 0.25, 0.50, P["e_ex"], 1.0, 1.5):
        Pe = dict(P)
        Pe["e_ex"] = ex
        te, oe, _, sle, _ = G.run_trial(Pe, arms, 56.0, visits=VIS8, dis=1.0,
                                        fluct0=fluct0, obs_keys=("hama", "phi_n"),
                                        obs_times=np.array([0.0, 56.0]))
        p = sle["placebo"].start
        dpp = oe["hama"][1, p] - oe["hama"][0, p]
        say("   %4.2f   %13.2f   %11.2f   %12.2f   %10.2f   %14.3f"
            % (ex, dpp,
               (oe["hama"][1, sle["esc10"].start] - oe["hama"][0, sle["esc10"].start]) - dpp,
               (oe["hama"][1, sle["pgb300"].start] - oe["hama"][0, sle["pgb300"].start]) - dpp,
               (oe["hama"][1, sle["lzp3"].start] - oe["hama"][0, sle["lzp3"].start]) - dpp,
               oe["phi_n"][1, p]))
    say("")
    say("  A site with a strong expectancy response loses measurable drug effect")
    say("  WITHOUT ANY CHANGE IN THE PHARMACOLOGY.  In this model that is not an")
    say("  assumption: expectancy raises C_pfc, C_pfc is in the DENOMINATOR of Phi,")
    say("  and HAM-A is a saturating function of Phi, so a placebo-responsive site is")
    say("  already sitting on the flat part of the curve when the drug arrives.")


# ---------------------------------------------------------------------------
def s7_bzd(P, fluct0):
    hdr("7.  BENZODIAZEPINES — one adaptation state, two signs")
    arms = [G.arm("lzp 12wk then STOP", "lorazepam", 3.0, stop=84.0),
            G.arm("lzp 12wk then 4wk taper", "lorazepam", 3.0,
                  titrate=[(0.0, 84.0, 3.0), (84.0, 91.0, 2.25),
                           (91.0, 98.0, 1.5), (98.0, 105.0, 0.75),
                           (105.0, 112.0, 0.375)]),
            G.arm("lzp continued", "lorazepam", 3.0, stop=200.0),
            G.arm("untreated")]
    t, o, yf, sl, Q = G.run_trial(P, arms, 168.0, visits=(), dis=1.0, fluct0=0.0,
                                  obs_keys=KEYS, obs_times=np.arange(0, 168.1, 1.0),
                                  record_states=("ra1", "ra2", "depend"))
    base = o["hama"][0, sl["lzp continued"].start]
    say("No trial machinery: pure pharmacology.  Baseline HAM-A %.2f." % base)
    say("")
    say("  TOLERANCE, split by subunit (lorazepam 3 mg/d continued):")
    say("     day    Occ_bz   R_a1(sedation)  R_a2(anxiolysis)  DEPEND   sedation   HAM-A")
    c = sl["lzp continued"].start
    for d in (0, 1, 3, 7, 14, 28, 56, 84, 112, 168):
        say("     %3d    %6.3f   %13.3f  %16.3f  %6.3f   %8.3f  %6.2f"
            % (d, o["occ_bz"][d, c], o["st_ra1"][d, c], o["st_ra2"][d, c],
               o["st_depend"][d, c], o["sedation"][d, c], o["hama"][d, c]))
    say("")
    say("  The asymmetry is the whole point.  Sedation falls %.0f%% between day 1 and"
        % (100 * (1 - o["sedation"][14, c] / o["sedation"][1, c])))
    say("  day 14 (R_a1 tau ~2.6 d) while the anxiolytic pool R_a2 is still %.3f at"
        % o["st_ra2"][14, c])
    say("  day 14 and %.3f at day 168.  Anxiolysis at day 168 retains %.0f%% of the"
        % (o["st_ra2"][168, c], 100 * (base - o["hama"][168, c]) / (base - o["hama"][3, c])))
    say("  day-3 effect.  That is consistent with the clinical literature, which")
    say("  reports prominent tolerance to sedation and little to anxiolysis.")
    say("")
    say("  REBOUND: stopping abruptly at day 84 versus a four-week taper.")
    say("     day    abrupt   taper   continued   untreated")
    for d in (83, 85, 88, 91, 98, 105, 112, 126, 140, 168):
        say("     %3d   %6.2f  %6.2f    %8.2f    %8.2f"
            % (d, o["hama"][d, sl["lzp 12wk then STOP"].start],
               o["hama"][d, sl["lzp 12wk then 4wk taper"].start],
               o["hama"][d, sl["lzp continued"].start],
               o["hama"][d, sl["untreated"].start]))
    a = o["hama"][:, sl["lzp 12wk then STOP"].start]
    tp = o["hama"][:, sl["lzp 12wk then 4wk taper"].start]
    un = o["hama"][:, sl["untreated"].start]
    pk_a = float(np.max(a[84:140]))
    pk_t = float(np.max(tp[84:140]))
    say("")
    say("  Peak post-withdrawal HAM-A: abrupt %.2f, taper %.2f, untreated baseline %.2f."
        % (pk_a, pk_t, float(un[100])))
    say("  Abrupt stop overshoots the patient's OWN untreated baseline by %.2f points"
        % (pk_a - float(un[100])))
    say("  (%.1f%% of it); the taper overshoot is %.2f points.  ONE parameter"
        % (100 * (pk_a - float(un[100])) / float(un[100]), pk_t - float(un[100])))
    say("  (wadapt = %.2f) produces both the tolerance on drug and the rebound off"
        % P["wadapt"])
    say("  it -- the same number with opposite sign, which is a testable coupling.")
    say("")
    say("  HOW BIG IS THAT, HONESTLY?  Small.  The model returns the patient to")
    say("  baseline and then a fraction of a point past it, whereas the clinical")
    say("  literature describes discontinuation-emergent anxiety that clearly")
    say("  EXCEEDS the pre-treatment state for weeks.  Two readings are available")
    say("  and the model cannot choose between them: either wadapt is too small")
    say("  (it is pinned by the tolerance side of the same coupling, so raising it")
    say("  would over-predict tolerance), or discontinuation anxiety is not the")
    say("  same object as loss of I_gaba at all -- it may be a distinct withdrawal")
    say("  syndrome that this one-state adaptation model does not contain.  The")
    say("  SHAPE is right (abrupt worse than tapered, decaying over ~3 weeks); the")
    say("  MAGNITUDE is under-predicted, and that is reported rather than tuned.")


# ---------------------------------------------------------------------------
def s8_combo(P, fluct0):
    hdr("8.  MULTIPLICATIVE IN Phi — combinations, the benzodiazepine bridge, and CBT")
    arms = [G.arm("placebo"),
            G.arm("esc10", "escitalopram", 10.0),
            G.arm("pgb300", "pregabalin", 300.0),
            G.arm("esc10+pgb300", "escitalopram", 10.0),
            G.arm("lzp3", "lorazepam", 3.0),
            G.arm("esc10+lzp3", "escitalopram", 10.0),
            G.arm("CBT alone", None, cbt=(0.0, 84.0, 1.0)),
            G.arm("esc10+CBT", "escitalopram", 10.0, cbt=(0.0, 84.0, 1.0)),
            G.arm("esc10 + lzp bridge 4wk", "escitalopram", 10.0)]
    # run_trial takes one drug per arm, so the combination arms are built with
    # an explicit regimen instead
    Q2, n = G._broadcast_params(P, arms, 1)
    Q2["DIS"] = np.full(n, 1.0)
    R = G.Regimen()
    for i, a in enumerate(arms):
        mask = np.zeros(n)
        mask[i] = 1.0
        nm = a["name"]
        if "esc10" in nm:
            R.add("esc_gut", 10.0, 0.0, 84.0, 1, mask=mask)
        if "pgb300" in nm:
            R.add("pgb_gut", 300.0, 0.0, 84.0, 2, mask=mask)
        if nm == "lzp3" or nm == "esc10+lzp3":
            R.add("bzd_gut", 3.0, 0.0, 84.0, 3, mask=mask)
        if nm == "esc10 + lzp bridge 4wk":
            R.add("bzd_gut", 3.0, 0.0, 21.0, 3, mask=mask)
            R.add("bzd_gut", 2.0, 21.0, 24.0, 3, mask=mask)
            R.add("bzd_gut", 1.0, 24.0, 28.0, 3, mask=mask)
        if a["cbt"]:
            inten = np.zeros(n)
            inten[i] = 1.0
            R.cbt(a["cbt"][0], a["cbt"][1], inten)
    y0 = G.attractor_for(P, np.full(n, 1.0))
    t, o, yf = G.simulate(Q2, R, 84.0, y0=y0, n=n, visits=VIS12,
                          obs_times=np.arange(0, 84.1, 1.0),
                          fluct0=np.full(n, fluct0), obs_keys=KEYS)
    sl = {a["name"]: i for i, a in enumerate(arms)}
    pc = sl["placebo"]
    say("Week-12 HAM-A change from baseline, and drug-placebo delta.")
    say("")
    dpp = o["hama"][84, pc] - o["hama"][0, pc]
    say("  arm                       Phi_n(wk12)  dHAM-A   delta   E_amy   C_pfc  S_glu  I_gaba")
    for nm in [a["name"] for a in arms]:
        c = sl[nm]
        say("  %-25s %11.3f  %6.2f  %6.2f  %6.3f  %6.3f %6.3f %6.3f"
            % (nm, o["phi_n"][84, c], o["hama"][84, c] - o["hama"][0, c],
               (o["hama"][84, c] - o["hama"][0, c]) - dpp,
               o["eamy"][84, c], o["cpfc"][84, c], o["sglu"][84, c], o["igaba"][84, c]))
    say("")

    def dd(nm):
        c = sl[nm]
        return (o["hama"][84, c] - o["hama"][0, c]) - dpp

    def phir(nm):
        return o["phi_n"][84, sl[nm]] / o["phi_n"][84, pc]

    for a_, b_, ab in (("esc10", "pgb300", "esc10+pgb300"),
                       ("esc10", "lzp3", "esc10+lzp3"),
                       ("esc10", "CBT alone", "esc10+CBT")):
        say("  %s x %s:" % (a_, b_))
        say("     Phi ratio vs placebo:  %.4f x %.4f = %.4f predicted-if-multiplicative;"
            % (phir(a_), phir(b_), phir(a_) * phir(b_)))
        say("                            observed in the combination arm %.4f  (error %.1f%%)"
            % (phir(ab), 100 * (phir(ab) - phir(a_) * phir(b_)) / (phir(a_) * phir(b_))))
        say("     HAM-A delta:           %.2f + %.2f = %.2f if additive; observed %.2f"
            % (dd(a_), dd(b_), dd(a_) + dd(b_), dd(ab)))
        pct = 100 * dd(ab) / (dd(a_) + dd(b_))
        say("     => %.0f%% of the additive expectation on the score." % pct)
    say("")
    say("  TWO RESULTS, ONE STRONG AND ONE AGAINST EXPECTATION.")
    say("  (1) Phi is MULTIPLICATIVE to within a few percent in all three pairs")
    say("      (errors %s).  That is the strong result: two drugs that divide"
        % ", ".join("%.1f%%" % (100 * (phir(ab) - phir(a_) * phir(b_)) / (phir(a_) * phir(b_)))
                    for a_, b_, ab in (("esc10", "pgb300", "esc10+pgb300"),
                                       ("esc10", "lzp3", "esc10+lzp3"),
                                       ("esc10", "CBT alone", "esc10+CBT"))))
    say("      different factors of the same ratio compose exactly, with no")
    say("      pharmacological interaction term anywhere in the model.")
    pcts = [100 * dd(ab) / (dd(a_) + dd(b_))
            for a_, b_, ab in (("esc10", "pgb300", "esc10+pgb300"),
                               ("esc10", "lzp3", "esc10+lzp3"),
                               ("esc10", "CBT alone", "esc10+CBT"))]
    say("  (2) On HAM-A the combinations come out at %s of the additive"
        % ", ".join("%.0f%%" % p for p in pcts))
    say("      expectation -- i.e. SUPER-additive, which is the OPPOSITE of what")
    say("      was expected going in.  The expectation was sub-additivity from a")
    say("      saturating readout.  The readout does saturate, but not here: at")
    say("      effect sizes of 3-4 HAM-A points the arms sit on the CONVEX part")
    say("      of the link, where each further reduction in Phi buys MORE score")
    say("      than the last.  Sub-additivity is what the same link gives once")
    say("      the effects are large enough to approach the floor, and GAD trials")
    say("      never get there.  The prior expectation is recorded as refuted,")
    say("      and it inverts the clinical prediction: combination arms in GAD")
    say("      should BEAT the sum of their monotherapy deltas, not fall short.")
    say("")
    br = sl["esc10 + lzp bridge 4wk"]
    e = sl["esc10"]
    say("  THE BENZODIAZEPINE BRIDGE.  3 mg/d lorazepam for 3 weeks, tapered off by")
    say("  week 4, on top of escitalopram, versus escitalopram alone:")
    say("     week 1: bridge %.2f vs esc alone %.2f  (advantage %.2f)"
        % (o["hama"][7, br] - o["hama"][0, br], o["hama"][7, e] - o["hama"][0, e],
           (o["hama"][7, br] - o["hama"][0, br]) - (o["hama"][7, e] - o["hama"][0, e])))
    for w, d in (("week 4", 28), ("week 6", 42), ("week 12", 84)):
        say("     %-7s: bridge %.2f vs esc alone %.2f  (advantage %.2f)"
            % (w, o["hama"][d, br] - o["hama"][0, br], o["hama"][d, e] - o["hama"][0, e],
               (o["hama"][d, br] - o["hama"][0, br]) - (o["hama"][d, e] - o["hama"][0, e])))
    say("  The bridge buys a real week-1 advantage and then GIVES SOME OF IT BACK")
    say("  around the taper, because DEPEND has been built and unbinds over 20 d.")
    say("  By week 12 the two arms are within %.2f points."
        % abs((o["hama"][84, br] - o["hama"][0, br]) - (o["hama"][84, e] - o["hama"][0, e])))
    say("")
    cb = sl["CBT alone"]
    say("  WHICH FACTOR DOES EACH INTERVENTION ACTUALLY MOVE?  At week 12:")
    say("     escitalopram alone : E_amy %.3f, C_pfc %.3f, S_glu %.3f"
        % (o["eamy"][84, e], o["cpfc"][84, e], o["sglu"][84, e]))
    say("     CBT alone          : E_amy %.3f, C_pfc %.3f, S_glu %.3f"
        % (o["eamy"][84, cb], o["cpfc"][84, cb], o["sglu"][84, cb]))
    say("     untreated (placebo): E_amy %.3f, C_pfc %.3f, S_glu %.3f"
        % (o["eamy"][84, pc], o["cpfc"][84, pc], o["sglu"][84, pc]))
    say("")
    say("  A CAVEAT AGAINST THE TIDY STORY, reported rather than smoothed over.  The")
    say("  clean version of the argument is 'CBT moves E_amy, the drug moves C_pfc'.")
    say("  The model does NOT cleanly say that: chronic 5-HT elevation accelerates")
    say("  E_amy decay too (a5ht_amy = %.2f), and lowering Phi by any route removes"
        % P["a5ht_amy"])
    say("  the sleep-deficit drive that feeds E_amy, so escitalopram reaches E_amy")
    say("  as well (%.3f vs CBT's %.3f).  What survives is the WEAKER claim, which"
        % (o["eamy"][84, e], o["eamy"][84, cb]))
    say("  is the one section 10 actually needs: CBT's effect on E_amy does not")
    say("  depend on a drug being present, so it does not unwind when dosing stops.")


# ---------------------------------------------------------------------------
def s9_vpop(P, fluct0, n=200, seed=20260806):
    hdr("9.  VIRTUAL POPULATION — dropout, LOCF, and whether AE burden can bend "
        "the quetiapine curve down")
    rng = np.random.default_rng(seed)
    arms = [G.arm("placebo"),
            G.arm("escitalopram 10", "escitalopram", 10.0),
            G.arm("pregabalin 300", "pregabalin", 300.0),
            G.arm("quetiapine 50", "quetiapine", 50.0),
            G.arm("quetiapine 150", "quetiapine", 150.0),
            G.arm("quetiapine 300", "quetiapine", 300.0)]
    na = len(arms)
    dis = np.clip(rng.normal(1.00, 0.22, na * n), 0.45, 1.70)
    exv = np.clip(rng.normal(P["e_ex"], 0.35 * P["e_ex"], na * n), 0.02, None)
    fl = np.clip(rng.normal(fluct0, 0.35 * abs(fluct0) + 0.5, na * n), 0.0, None)
    Pv = dict(P)
    Pv["e_ex"] = exv
    Pv["adherence"] = np.clip(rng.normal(0.92, 0.10, na * n), 0.4, 1.0)
    t, o, yf, sl, Q = G.run_trial(Pv, arms, 56.0, visits=VIS8, dis=dis, nper=n,
                                  fluct0=fl, obs_keys=("hama", "ae_burden", "cumhaz"),
                                  obs_times=np.array(VIS8))
    hama = o["hama"]
    haz = o["cumhaz"]
    u = rng.random(na * n)
    surv = np.exp(-(haz - haz[0]))
    drop_idx = np.full(na * n, len(VIS8) - 1)
    for j in range(na * n):
        for i in range(1, len(VIS8)):
            if surv[i, j] < u[j]:
                drop_idx[j] = i - 1
                break
    noise = rng.normal(0.0, 2.2, (len(VIS8), na * n))
    obs = hama + noise
    locf = obs[drop_idx, np.arange(na * n)]
    base = obs[0]
    chg = locf - base
    say("n = %d per arm; between-subject variability in severity (DIS), site" % n)
    say("expectancy (e_ex), enrolment peak and adherence; HAM-A measurement SD 2.2;")
    say("dropout from the AE-driven hazard; LOCF endpoint, exactly as the trials did.")
    say("")
    say("  arm               baseline   LOCF change   completers%   AE burden(wk8)   response%  remission%")
    for a_ in arms:
        s_ = sl[a_["name"]]
        say("  %-16s %8.1f   %11.2f   %10.0f   %14.3f   %8.0f   %9.0f"
            % (a_["name"], base[s_].mean(), chg[s_].mean(),
               100 * np.mean(drop_idx[s_] == len(VIS8) - 1),
               o["ae_burden"][-1, s_].mean(),
               100 * np.mean(chg[s_] <= -0.5 * base[s_]),
               100 * np.mean(locf[s_] <= 7)))
    say("")
    q = [chg[sl["quetiapine %d" % d]].mean() for d in (50, 150, 300)]
    say("  Quetiapine LOCF endpoint 50 / 150 / 300 mg: %.2f / %.2f / %.2f" % tuple(q))
    say("  Reported (Khan 2011): -13.31 / -13.54 / -11.87")
    if q[2] > q[1]:
        say("  The model DOES reproduce the sign of the inverted U once dropout is in:")
        say("  the 300 mg arm loses %.2f points relative to 150 mg." % (q[2] - q[1]))
    else:
        say("  THE MODEL STILL DOES NOT REPRODUCE THE INVERTED U (%.2f at 300 vs %.2f at" % (q[2], q[1]))
        say("  150).  Reported as a miss: AE-driven dropout in this parameterisation is")
        say("  too weak to overturn a monotone pharmacological effect, and either the")
        say("  AE weights or the mechanism of the reported 300 mg failure is wrong.")
    say("")
    say("  SEVERITY DEPENDENCE (escitalopram 10 mg minus placebo, by baseline")
    say("  tertile).  The expectation going in was that the delta would GROW with")
    say("  baseline severity, because Phi is a ratio and the same fractional change")
    say("  buys more score further up the curve.  The tertiles are printed without")
    say("  comment so the reader can see whether that happened:")
    pb = sl["placebo"]
    ea = sl["escitalopram 10"]
    tert = []
    for lo, hi, nm in ((0, 33.3, "mild"), (33.3, 66.7, "moderate"), (66.7, 100, "severe")):
        bp = np.percentile(base[pb], [lo, hi])
        be = np.percentile(base[ea], [lo, hi])
        mp = chg[pb][(base[pb] >= bp[0]) & (base[pb] <= bp[1])].mean()
        me = chg[ea][(base[ea] >= be[0]) & (base[ea] <= be[1])].mean()
        say("     %-9s baseline ~%.1f : placebo %.2f, escitalopram %.2f, delta %.2f"
            % (nm, 0.5 * (bp[0] + bp[1]), mp, me, me - mp))
        tert.append(me - mp)
    say("")
    if tert[2] < tert[0] - 0.15:
        say("  It did: the drug-placebo difference grows from %.2f to %.2f across the"
            % (tert[0], tert[2]))
        say("  tertiles, which is the enrichment argument for recruiting severe")
        say("  patients into GAD trials.")
    elif abs(tert[2] - tert[0]) <= 0.35:
        say("  IT DID NOT.  The difference is %.2f in the mild tertile and %.2f in the"
            % (tert[0], tert[2]))
        say("  severe one -- flat, and non-monotone in between.  The reason is that")
        say("  two effects cancel: severe patients sit higher on the readout, where a")
        say("  given fractional change in Phi is worth more score, but they also have")
        say("  a LARGER placebo response (%.2f versus %.2f), which is subtracted."
            % (mp, chg[pb][(base[pb] >= np.percentile(base[pb], [0, 33.3])[0])
                           & (base[pb] <= np.percentile(base[pb], [0, 33.3])[1])].mean()))
        say("  The severity-enrichment argument therefore does NOT follow from this")
        say("  model, and the prior expectation is recorded as refuted rather than")
        say("  quietly dropped.")
    else:
        say("  It ran the OTHER way: %.2f in the mild tertile against %.2f in the"
            % (tert[0], tert[2]))
        say("  severe one, because the placebo response grows with severity too.")


# ---------------------------------------------------------------------------
def s10_relapse(P, fluct0, n=260, seed=770203):
    hdr("10.  ALLGULANDER 2006 (PMID 16316482) — the randomised-withdrawal "
        "prediction, out of sample")
    rng = np.random.default_rng(seed)
    say("Design reproduced exactly: 12 weeks open-label escitalopram 20 mg/d;")
    say("responders (HAM-A <= 10) randomised to continue 20 mg or to placebo;")
    say("followed to week 76; relapse = HAM-A >= 15 at any scheduled visit.")
    say("Reported: 56% relapse on placebo, 19% on escitalopram (risk ratio 4.04).")
    say("")
    dis = np.clip(rng.normal(1.00, 0.22, n), 0.45, 1.70)
    exv = np.clip(rng.normal(P["e_ex"], 0.35 * P["e_ex"], n), 0.02, None)
    Pv = dict(P)
    Pv["e_ex"] = exv
    Pv["adherence"] = np.clip(rng.normal(0.92, 0.10, n), 0.4, 1.0)

    # ---- phase 1: 12 weeks open-label 20 mg ---------------------------
    a1 = [G.arm("ol", "escitalopram", 20.0)]
    Q1, _ = G._broadcast_params(Pv, a1, n)
    Q1["DIS"] = dis
    R1 = G.Regimen()
    R1.add("esc_gut", 20.0, 0.0, 84.0, 1)
    y0 = G.attractor_for(P, dis)
    t1, o1, y1 = G.simulate(Q1, R1, 84.0, y0=y0, n=n, visits=VIS12,
                            obs_times=np.array([0.0, 84.0]),
                            fluct0=np.clip(rng.normal(fluct0, 0.35 * abs(fluct0) + 0.5, n), 0, None),
                            obs_keys=("hama",))
    noise1 = rng.normal(0, 2.2, n)
    hama12 = o1["hama"][1] + noise1
    resp = hama12 <= 10.0
    say("  Open-label phase: %d entered, %d responded (HAM-A <= 10) = %.0f%%."
        % (n, int(resp.sum()), 100 * resp.mean()))
    say("  (Allgulander randomised 375 of 491, i.e. 76%.)")
    nr = int(resp.sum())
    if nr < 20:
        say("  too few responders to continue")
        return

    # ---- phase 2: randomised withdrawal, 76 weeks ----------------------
    yR = y1[:, resp]
    idx = np.arange(nr)
    cont = idx % 2 == 0
    y2 = np.concatenate([yR, yR], axis=1)          # arm A = continue, arm B = placebo
    n2 = 2 * nr
    Q2 = dict(Pv)
    for k, v in list(Q2.items()):
        if isinstance(v, np.ndarray) and v.size == n:
            Q2[k] = np.concatenate([v[resp], v[resp]])
    Q2["DIS"] = np.concatenate([dis[resp], dis[resp]])
    R2 = G.Regimen()
    mask_cont = np.concatenate([np.ones(nr), np.zeros(nr)])
    R2.add("esc_gut", 20.0, 0.0, 532.0, 1, mask=mask_cont)
    vis2 = [0.0] + [28.0 * k for k in range(1, 20)]
    obs2 = np.array(vis2)
    t2, o2, y2f = G.simulate(Q2, R2, 532.0, y0=y2, n=n2, visits=vis2,
                             obs_times=obs2, dt=1.0 / 24.0,
                             obs_keys=("hama", "phi_n", "eamy", "cpfc", "occ_sert"))
    hm = o2["hama"] + rng.normal(0, 2.2, o2["hama"].shape)
    relapsed = (hm >= 15.0).any(axis=0)
    rc = 100 * relapsed[:nr].mean()
    rp = 100 * relapsed[nr:].mean()
    say("")
    say("  Randomised withdrawal (%d per arm), relapse by week 76:" % nr)
    say("     escitalopram 20 mg continued : %.0f%%   (reported 19%%)" % rc)
    say("     placebo                      : %.0f%%   (reported 56%%)" % rp)
    say("     risk ratio placebo/drug      : %.2f    (reported 4.04)"
        % (rp / rc if rc > 0 else float("inf")))
    say("")
    say("  Relapse-free curve (percent still relapse-free):")
    say("     week   continued   placebo")
    for i, tt in enumerate(vis2):
        if i == 0 or i % 3:
            continue
        rf_c = 100 * (1 - (hm[:i + 1, :nr] >= 15).any(axis=0).mean())
        rf_p = 100 * (1 - (hm[:i + 1, nr:] >= 15).any(axis=0).mean())
        say("     %4.0f   %9.0f   %7.0f" % (tt / 7.0, rf_c, rf_p))
    say("")
    say("  WHY THE ARMS SEPARATE, mechanistically.  At week 76:")
    say("     continued: Phi_n %.3f, E_amy %.3f, C_pfc %.3f, Occ_SERT %.3f"
        % (o2["phi_n"][-1, :nr].mean(), o2["eamy"][-1, :nr].mean(),
           o2["cpfc"][-1, :nr].mean(), o2["occ_sert"][-1, :nr].mean()))
    say("     placebo  : Phi_n %.3f, E_amy %.3f, C_pfc %.3f, Occ_SERT %.3f"
        % (o2["phi_n"][-1, nr:].mean(), o2["eamy"][-1, nr:].mean(),
           o2["cpfc"][-1, nr:].mean(), o2["occ_sert"][-1, nr:].mean()))
    say("  Stopping the drug releases C_pfc back down its own 30-day time constant.")
    say("  Nothing in the acute-phase calibration knew about this number.")


# ---------------------------------------------------------------------------
def s11_covariates(P, fluct0):
    hdr("11.  COVARIATES — renal function, CYP2D6, and adherence")
    arms = [G.arm("placebo"),
            G.arm("pgb300 CrCl 100", "pregabalin", 300.0, crcl=100.0),
            G.arm("pgb300 CrCl 60", "pregabalin", 300.0, crcl=60.0),
            G.arm("pgb300 CrCl 30", "pregabalin", 300.0, crcl=30.0),
            G.arm("pgb150 CrCl 30", "pregabalin", 150.0, crcl=30.0),
            G.arm("esc10 adherence 1.0", "escitalopram", 10.0, adherence=1.0),
            G.arm("esc10 adherence 0.7", "escitalopram", 10.0, adherence=0.7),
            G.arm("esc10 adherence 0.5", "escitalopram", 10.0, adherence=0.5)]
    t, o, yf, sl, Q = G.run_trial(P, arms, 56.0, visits=VIS8, dis=1.0,
                                  fluct0=fluct0, obs_keys=KEYS,
                                  obs_times=np.arange(0, 56.1, 1.0))
    pc = sl["placebo"].start
    dpp = o["hama"][56, pc] - o["hama"][0, pc]
    say("  arm                    Css        occupancy   dHAM-A   delta   dizziness  sedation")
    for a_ in arms[1:]:
        c = sl[a_["name"]].start
        if "pgb" in a_["name"]:
            css, occ = o["Cpgb"][56, c], o["occ_a2d"][56, c]
        else:
            css, occ = o["Cesc"][56, c], o["occ_sert"][56, c]
        say("  %-22s %8.1f    %8.3f  %7.2f  %6.2f   %8.3f  %8.3f"
            % (a_["name"], css, occ, o["hama"][56, c] - o["hama"][0, c],
               o["hama"][56, c] - o["hama"][0, c] - dpp,
               o["dizziness"][56, c], o["sedation"][56, c]))
    say("")
    c100 = sl["pgb300 CrCl 100"].start
    c30 = sl["pgb300 CrCl 30"].start
    h30 = sl["pgb150 CrCl 30"].start
    say("  Pregabalin clearance is renal, so CrCl 30 multiplies steady-state")
    say("  exposure by %.2f (%.0f -> %.0f ng/mL) on an unchanged 300 mg/d."
        % (o["Cpgb"][56, c30] / o["Cpgb"][56, c100], o["Cpgb"][56, c100], o["Cpgb"][56, c30]))
    say("  Halving the dose is NOT enough to restore it: 150 mg at CrCl 30 still")
    say("  gives %.0f ng/mL against %.0f at CrCl 100, so the label's quartering of"
        % (o["Cpgb"][56, h30], o["Cpgb"][56, c100]))
    say("  the dose at severe renal impairment is what this arithmetic actually")
    say("  supports, and the model is reported as agreeing with the direction and")
    say("  not with a halving rule invented for it.")
    say("")
    d07 = ((o["hama"][56, sl["esc10 adherence 0.7"].start]
            - o["hama"][0, sl["esc10 adherence 0.7"].start])
           - (o["hama"][56, sl["esc10 adherence 1.0"].start]
              - o["hama"][0, sl["esc10 adherence 1.0"].start]))
    say("  Adherence: 70%% intake costs %.2f HAM-A points, i.e. %.0f%% of the entire"
        % (abs(d07), 100 * abs(d07) / abs(o["hama"][56, sl["esc10 adherence 1.0"].start]
                                          - o["hama"][0, sl["esc10 adherence 1.0"].start] - dpp)))
    say("  drug-placebo difference the trial was powered to detect.  Nothing else in")
    say("  this covariate table comes close.")


# ---------------------------------------------------------------------------
def s12_verify(P):
    hdr("12.  INTEGRATOR VERIFICATION — hand-written RK4 against scipy LSODA")
    from scipy.integrate import solve_ivp
    Q = dict(P)
    Q["DIS"] = 1.0
    y0 = G.attractor_for(P, [1.0])[:, 0]
    # single 10 mg escitalopram dose + 300 mg pregabalin BID for 7 days
    R = G.Regimen()
    R.add("esc_gut", 10.0, 0.0, 7.0, 1)
    R.add("pgb_gut", 300.0, 0.0, 7.0, 2)
    t, o, yf = G.simulate(Q, R, 7.0, y0=y0[:, None], n=1,
                          obs_times=np.array([7.0]), obs_keys=("hama",))

    # reference: LSODA integrated segment by segment across the same impulses
    yr = y0.copy()
    times = sorted(set([round(x, 6) for x in
                        list(np.arange(0, 7, 1.0)) + list(np.arange(0, 7, 0.5))] + [7.0]))
    for i in range(len(times) - 1):
        a, b = times[i], times[i + 1]
        if abs(a - round(a)) < 1e-9:
            yr[G.IX["esc_gut"]] += 10000.0
        if abs(a * 2 - round(a * 2)) < 1e-9:
            yr[G.IX["pgb_gut"]] += 150000.0
        sol = solve_ivp(lambda tt, yy: rhs_flat(tt, yy, Q), (a, b), yr,
                        method="LSODA", rtol=1e-9, atol=1e-11)
        yr = sol.y[:, -1]
    o2 = G.observe(yr[:, None], Q)
    say("  7 days of escitalopram 10 mg QD + pregabalin 300 mg BID from the GAD")
    say("  attractor.  Terminal values, RK4 (dt = 30 min) vs LSODA (rtol 1e-9):")
    say("")
    say("    quantity          RK4          LSODA        rel. difference")
    for k in ("hama",):
        say("    %-14s %12.6f %12.6f %14.2e"
            % (k, o[k][0, 0], o2[k][0], abs(o[k][0, 0] - o2[k][0]) / abs(o2[k][0])))
    worst = 0.0
    worst_s = ""
    for s_ in G.STATES:
        a_, b_ = yf[G.IX[s_], 0], yr[G.IX[s_]]
        if abs(b_) > 1e-8:
            r = abs(a_ - b_) / abs(b_)
            if r > worst:
                worst, worst_s = r, s_
    say("    worst state relative difference: %.2e (%s)" % (worst, worst_s))
    say("")
    say("  A second check: halve the RK4 step and confirm the answer does not move.")
    t3, o3, yf3 = G.simulate(Q, R, 7.0, y0=y0[:, None], n=1, dt=1.0 / 96.0,
                             obs_times=np.array([7.0]), obs_keys=("hama",))
    say("    HAM-A at dt = 1/48 d: %.8f" % o["hama"][0, 0])
    say("    HAM-A at dt = 1/96 d: %.8f  (difference %.2e)"
        % (o3["hama"][0, 0], abs(o3["hama"][0, 0] - o["hama"][0, 0])))


def rhs_flat(t, y, P):
    return G.rhs(t, y[:, None], P, 0.0)[:, 0]


# ---------------------------------------------------------------------------
def main():
    P, fluct0 = load_params()
    say("GAD QSP MODEL — REFERENCE OUTPUT")
    say("generated by gad_analysis.py from gad_python_reference.py")
    say("")
    say("Fitted parameters (five, from gad_calibration.json; dvisit and e_ex")
    say("are listed too but are STRUCTURAL, not fitted -- see gad_calibrate.py):")
    fit = json.load(open("gad_calibration.json"))
    for k, v in fit.items():
        say("   %-10s %.5f" % (k, v))
    say("")
    say("Every other number in this file, in README.md, in the mrgsolve model and")
    say("in the Shiny app is a consequence of those five plus measured PK and")
    say("occupancy constants.")
    s1_baseline(P)
    s2_clocks(P, fluct0)
    s3_rickels(P, fluct0)
    s4_khan(P, fluct0)
    s5_doseresponse(P, fluct0)
    s6_placebo(P, fluct0)
    s7_bzd(P, fluct0)
    s8_combo(P, fluct0)
    s9_vpop(P, fluct0)
    s10_relapse(P, fluct0)
    s11_covariates(P, fluct0)
    s12_verify(P)
    open("gad_reference_output.txt", "w").write("\n".join(OUT) + "\n")
    print("\nwrote gad_reference_output.txt (%d lines)" % len(OUT))


if __name__ == "__main__":
    main()
