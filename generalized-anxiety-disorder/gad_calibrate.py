#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibrate the GAD QSP model.  FIVE parameters, five numbers.

HOW THIS ENDED UP BEING A SEPARABLE FIT, which is worth recording because the
first attempt was wrong in an instructive way.

The first attempt (kept as calib_stage1_joint.log) fitted fluct0, dvisit,
e_ex, emax_pgb and emax_a2 jointly by least squares.  It converged to a sum of
squared residuals of about 4.3, and essentially the whole of that residual sat
on ONE target: the week-1 placebo response, where the model reached -4.4 to
-4.6 against Khan 2011's reported -5.94.  Adding iterations did not help,
because the miss was structural rather than numerical.  The placebo arm has
two
components: an additive enrolment peak FLUCT that decays with rate kfl, and a
biological expectancy term that has to propagate through WORRY, SLEEPD and
AUTON, each of which has its own time constant of several days.  The
expectancy term therefore CANNOT be fast.  Only FLUCT can be -- and kfl had
been left at a hand-set 0.075/d (tau 13 d), far too slow to put half of the
eight-week placebo effect into the first week, which is what Khan's placebo
arm does.  No amount of moving the other four parameters can repair that.

So kfl became a fitted parameter.  A SECOND thing then went wrong, and it is
worth recording too, because it is a real identifiability result rather than a
bug.  When fluct0, kfl AND dvisit were all fitted to the three placebo
numbers, the optimiser drove dvisit to its lower bound (0.02) and put the
entire placebo response into the additive FLUCT term: fluct0 9.72 points with
tau 5.8 d.  That is a perfectly good fit of three numbers -- and it silently
switches OFF the biological expectancy channel, on which the model's
assay-sensitivity argument depends.

The lesson is that THREE SUMMARY NUMBERS CANNOT SEPARATE a selection artefact
that decays from a real top-down effect that grows.  Both produce a falling
placebo curve.  Fitting them jointly does not discover the split; it just
picks whichever is cheaper.

The split is therefore FIXED STRUCTURALLY rather than fitted:

    dvisit = 0.30   per-visit expectancy increment
    e_ex   = 0.55   how strongly expectancy raises C_pfc

and both are swept explicitly in the assay-sensitivity analysis, where their
influence is the object of study rather than a nuisance parameter.  Any claim
in this repository that depends on the size of the expectancy channel is
therefore a claim about a STRUCTURAL choice, and is labelled as such.

What is left is separable, and is solved as such:

  A.  fluct0, kfl  <- the three PLACEBO targets.  Placebo runs carry no doses,
      so they integrate at dt = 1/24 d instead of 1/48; the placebo curve is
      identical to 4 decimal places at dt = 1/8, 1/12 and 1/24.
  B.  k5ht_eff     <- 1-D solve on the escitalopram-placebo difference.
  C.  emax_pgb     <- 1-D solve on the pregabalin-placebo difference.
  D.  emax_a2      <- 1-D solve on the alprazolam-placebo difference.

FIVE fitted parameters (fluct0, kfl, k5ht_eff, emax_pgb, emax_a2)
against six calibration numbers.  Everything else in
the model is a measured PK or occupancy constant, or structure.
"""
import json

import numpy as np
from scipy.optimize import brentq, least_squares

import gad_python_reference as G

VIS_RICKELS = [0.0, 7.0, 14.0, 21.0, 28.0]
VIS_KHAN = [0.0, 7.0, 14.0, 28.0, 42.0, 56.0]

TGT_PBO_W1 = -5.94      # Khan 2011 (PMID 21694613) placebo, week 1
TGT_PBO_W8 = -11.10     # Khan 2011 placebo, week 8
TGT_PBO_W4 = -8.40      # Rickels 2005 (PMID 16143734) placebo, week 4
# The three EFFICACY knobs are fitted to drug-minus-placebo DIFFERENCES, not
# to absolute change scores.  The difference is the quantity a randomised
# trial actually establishes; the absolute change also contains the placebo
# arm, so fitting it makes every drug parameter absorb the placebo model's
# error.  An earlier version of this file fitted absolute change and the
# consequence was visible immediately: the model matched Rickels' alprazolam
# total of -10.9 exactly while getting its drug-placebo delta wrong by more
# than a factor of two.
TGT_D_ESC = -2.45       # escitalopram vs placebo, Slee 2019 NMA (PMID 30712879)
TGT_D_PGB = -2.79       # pregabalin  vs placebo, Slee 2019 NMA
TGT_D_ALP = -2.50       # alprazolam 1.5 mg vs placebo, Rickels 2005 week 4


def placebo_curve(P, fluct0, visits, tstop, dt=1.0 / 24.0):
    arms = [G.arm("pbo")]
    t, o, yf, sl, Q = G.run_trial(P, arms, tstop, visits=visits, dis=1.0,
                                  fluct0=fluct0, obs_keys=("hama",), dt=dt,
                                  obs_times=np.arange(0, tstop + 1e-9, 1.0))
    return o["hama"][:, 0]


def placebo_resid(x, P0, verbose=False):
    fluct0, kfl = x
    P = dict(P0)
    P.update(kfl=kfl)
    hK = placebo_curve(P, fluct0, VIS_KHAN, 56.0)
    hR = placebo_curve(P, fluct0, VIS_RICKELS, 28.0)
    d1 = hK[7] - hK[0]
    d8 = hK[56] - hK[0]
    d4 = hR[28] - hR[0]
    r = np.array([d1 - TGT_PBO_W1, d4 - TGT_PBO_W4, d8 - TGT_PBO_W8])
    if verbose:
        print("   fluct0 %.4f kfl %.5f | wk1 %.2f (%.2f) wk4 %.2f (%.2f)"
              " wk8 %.2f (%.2f)  SSR %.4f"
              % (fluct0, kfl, d1, TGT_PBO_W1, d4, TGT_PBO_W4,
                 d8, TGT_PBO_W8, float(r @ r)))
    return r


def drug_delta(P, fluct0, drug, mg, tstop=56.0, visits=None, **over):
    """Drug-minus-placebo HAM-A change at tstop, both arms in one run."""
    visits = VIS_KHAN if visits is None else visits
    arms = [G.arm("pbo"), G.arm("drug", drug, mg, **over)]
    t, o, yf, sl, Q = G.run_trial(P, arms, tstop, visits=visits, dis=1.0,
                                  fluct0=fluct0, obs_keys=("hama",),
                                  obs_times=np.array([0.0, tstop]))
    c = sl["drug"].start
    p = sl["pbo"].start
    dd = o["hama"][1, c] - o["hama"][0, c]
    dp = o["hama"][1, p] - o["hama"][0, p]
    return dd - dp, dd, dp


def main():
    P0 = G.default_params()
    G.calibrate_healthy(P0)
    P0["dvisit"] = 0.30
    print("STRUCTURAL, not fitted:  dvisit %.3f   e_ex %.3f"
          % (P0["dvisit"], P0["e_ex"]))
    print("(see the module docstring: three summary numbers cannot separate a")
    print(" decaying selection artefact from a growing top-down effect)")

    # ---- A: the placebo clock -----------------------------------------
    x0 = np.array([6.0, 0.18])
    print("A. placebo clock, start:")
    placebo_resid(x0, P0, verbose=True)
    solA = least_squares(placebo_resid, x0,
                         bounds=(np.array([0.0, 0.01]), np.array([16.0, 2.00])),
                         args=(P0, True), diff_step=0.08,
                         xtol=3e-4, ftol=3e-4, max_nfev=40)
    fluct0, kfl = solA.x
    dvisit = P0["dvisit"]
    print("A. fitted: fluct0 %.4f  kfl %.5f (tau %.1f d)"
          % (fluct0, kfl, np.log(2) / kfl))
    P0["kfl"] = kfl

    def solve(name, key, lo, hi, drug, mg, target, **kw):
        def f(v):
            P = dict(P0)
            P[key] = v
            d, dd, dp = drug_delta(P, fluct0, drug, mg, **kw)
            print("   %s %.4f -> delta %.3f (target %.3f)  [drug %.2f, placebo %.2f]"
                  % (key, v, d, target, dd, dp))
            return d - target
        print("%s:" % name)
        flo, fhi = f(lo), f(hi)
        if flo * fhi < 0:
            return brentq(f, lo, hi, xtol=(hi - lo) * 2e-3)
        print("   WARNING: target not bracketed on [%g, %g]; clamping" % (lo, hi))
        return lo if abs(flo) < abs(fhi) else hi

    # ---- B: the serotonergic efficacy scalar ---------------------------
    k5ht_eff = solve("B. escitalopram 10 mg (serotonergic arm)", "k5ht_eff",
                     0.005, 1.2, "escitalopram", 10.0, TGT_D_ESC)
    P0["k5ht_eff"] = k5ht_eff

    # ---- C: pregabalin --------------------------------------------------
    emax_pgb = solve("C. pregabalin 300 mg", "emax_pgb",
                     0.01, 0.95, "pregabalin", 300.0, TGT_D_PGB)
    P0["emax_pgb"] = emax_pgb

    # ---- D: alprazolam ---------------------------------------------------
    emax_a2 = solve("D. alprazolam 1.5 mg", "emax_a2", 0.05, 11.0,
                    "alprazolam", 1.5, TGT_D_ALP,
                    tstop=28.0, visits=VIS_RICKELS)

    fit = dict(fluct0=float(fluct0), kfl=float(kfl), dvisit=float(dvisit),
               e_ex=float(P0["e_ex"]), k5ht_eff=float(k5ht_eff),
               emax_pgb=float(emax_pgb), emax_a2=float(emax_a2))
    print(json.dumps(fit, indent=2))
    json.dump(fit, open("gad_calibration.json", "w"), indent=2)

    # ---- final residuals -------------------------------------------------
    P = dict(P0)
    P["emax_a2"] = emax_a2
    hK = placebo_curve(P, fluct0, VIS_KHAN, 56.0)
    hR = placebo_curve(P, fluct0, VIS_RICKELS, 28.0)
    de, _, _ = drug_delta(P, fluct0, "escitalopram", 10.0)
    dg, _, _ = drug_delta(P, fluct0, "pregabalin", 300.0)
    da, _, _ = drug_delta(P, fluct0, "alprazolam", 1.5,
                          tstop=28.0, visits=VIS_RICKELS)
    print("\nFINAL FIT vs the six calibration numbers")
    print("  placebo week 1        : %7.2f   target %7.2f" % (hK[7] - hK[0], TGT_PBO_W1))
    print("  placebo week 4        : %7.2f   target %7.2f" % (hR[28] - hR[0], TGT_PBO_W4))
    print("  placebo week 8        : %7.2f   target %7.2f" % (hK[56] - hK[0], TGT_PBO_W8))
    print("  escitalopram 10 delta : %7.2f   target %7.2f" % (de, TGT_D_ESC))
    print("  pregabalin 300 delta  : %7.2f   target %7.2f" % (dg, TGT_D_PGB))
    print("  alprazolam 1.5 delta  : %7.2f   target %7.2f" % (da, TGT_D_ALP))


if __name__ == "__main__":
    main()
