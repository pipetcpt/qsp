#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibration of the HDFN model.  EIGHT numbers are fitted here, in three
independent stages, and every one of them is a published summary statistic.
Everything the analysis script reports afterwards is out of sample.

STAGE 1  (v0, g_pl)        placental conveyor
    target  fetal:maternal total-IgG concentration ratio
            0.075 at 19.5 wk and 1.25 at 39 wk        (Malek 1996, PMID 8955500)

STAGE 2  (kops)            red-cell destruction -- ONE parameter, ONE target
    target  a maternal anti-D of 15 IU/mL (the UK "high risk" band) drives the
            fetus to Hb 0.65 MoM -- i.e. to its first transfusion -- at 26 wk,
            which is the mean gestational age of the first IUT in the cohort
            that measured the decline (Nishie 2012, 26.1 wk)
    NOT a target, therefore a PREDICTION: the 0.40 g/dL/day decline between the
    first and second IUT, and its fall over later intervals.  Kres is fixed
    STRUCTURALLY because it is not identifiable from these data.

STAGE 3  (ksens, fmh_ante) sensitisation and the silent antenatal bleed
    target  16% sensitisation with no prophylaxis at all
    target  1.6% with postpartum 300 ug only -- a residual which in this model
            is almost entirely the UNCOVERED third-trimester silent FMH, so
            fitting it IS measuring that bleed

STAGE 4  (ugt_birth, ugt_t50)  neonatal conjugation
    target  a HEALTHY term newborn peaks at 8 mg/dL total serum bilirubin
            at ~96 h and is below 3 mg/dL by 2 weeks.  No HDFN newborn number
            is used, so every neonatal HDFN result is out of sample.

The fitter is never shown: the 2nd/3rd IUT decline, the MCA-PSV threshold, the
UNITY trial, the residual risk under combined antenatal+postpartum prophylaxis,
the anti-K contrast, or any neonatal HDFN number.
"""
import json
import math
import sys

import numpy as np
from scipy.optimize import brentq, least_squares

import hdfn_python_reference as M

# The log is opened LAZILY.  hdfn_analysis.py imports this module to re-use the
# maternal-only prophylaxis integrator, and an eager open(..., "w") at import
# time silently truncated the calibration log it had just written.
LOG = None


def say(*a):
    global LOG
    s = " ".join(str(x) for x in a)
    print(s)
    if LOG is None:
        LOG = open("calib.log", "a")
    LOG.write(s + "\n")
    LOG.flush()


# ==============================================================================
# STAGE 1 -- the placental conveyor
# ==============================================================================
def fm_ratio(p, gas=(19.5, 39.0)):
    c = M.Ctl()
    r = M.simulate(p, c, ga0=11.0, anti_d_iu=0.0, protocol="none",
                   ga_deliver=41.0, postnatal_days=1.0)
    return [M.at_ga(r, "Cf_ig", g) / M.at_ga(r, "Cig", g) for g in gas]


def stage1():
    say("\n=== STAGE 1: placental conveyor (v0, g_pl) ===")
    say("targets: F/M total IgG = 0.075 at 19.5 wk, 1.25 at 39 wk (Malek 1996)")

    def resid(x):
        p = dict(M.P)
        p["v0"], p["g_pl"] = math.exp(x[0]), x[1]
        got = fm_ratio(p)
        r = [math.log(got[0] / 0.075), math.log(got[1] / 1.25)]
        say("   v0=%.5g g_pl=%.4f -> F/M %.4f %.4f  resid %.4f %.4f"
            % (math.exp(x[0]), x[1], got[0], got[1], r[0], r[1]))
        return r

    sol = least_squares(resid, [math.log(0.002), 0.30], diff_step=0.06,
                        xtol=1e-4, ftol=1e-4)
    v0, g_pl = math.exp(sol.x[0]), sol.x[1]
    say("   FITTED v0 = %.6g /d, g_pl = %.4f /wk" % (v0, g_pl))
    return v0, g_pl


# ==============================================================================
# STAGE 2 -- destruction potency and RES saturation
# ==============================================================================
def severe_case(p, anti_d=15.0, protocol="none", ga_deliver=34.0, ga0=14.0,
                seed=7):
    c = M.Ctl()
    return M.simulate(p, c, ga0=ga0, anti_d_iu=anti_d, protocol=protocol,
                      ga_deliver=ga_deliver, postnatal_days=1.0, seed=seed), c


def ga_at_mom(res, mom=0.65):
    """First gestational age at which Hb falls below `mom` multiples of median."""
    x = res["hbmom"]
    idx = np.where(x < mom)[0]
    if len(idx) == 0:
        return 99.0
    i = idx[0]
    return float(res["ga"][i])


def iut_decline(p, anti_d=15.0, seed=7):
    """Simulate the standard IUT course and return the mean g/dL/day fall in the
    interval after IUT #1, #2 and #3."""
    c = M.Ctl()
    r = M.simulate(p, c, ga0=14.0, anti_d_iu=anti_d, protocol="mca",
                   ga_deliver=35.0, postnatal_days=1.0, seed=seed)
    out = []
    for k, (ga, hbpre, hbpost, n) in enumerate(c.iut_log):
        if k + 1 < len(c.iut_log):
            ga2, hbpre2, _, _ = c.iut_log[k + 1]
            days = (ga2 - ga) * 7.0
            if days > 3:
                out.append((hbpost - hbpre2) / days)
    return out, c, r


def stage2(p0):
    say("\n=== STAGE 2: destruction (kops) -- one parameter, one target ===")
    say("target: anti-D 15 IU/mL -> Hb 0.65 MoM at 26.0 wk")

    def f(lk):
        p = dict(p0)
        p["kops"] = math.exp(lk)
        ga65 = ga_at_mom(severe_case(p)[0], 0.65)
        say("   kops=%.5g -> GA(0.65 MoM) = %.2f wk" % (p["kops"], ga65))
        return ga65 - 26.0

    lk = brentq(f, math.log(0.03), math.log(1.2), xtol=0.01)
    say("   FITTED kops = %.6g" % math.exp(lk))
    return math.exp(lk)


# ==============================================================================
# STAGE 3 -- sensitisation and immune deviation
# ==============================================================================
# Fetomaternal haemorrhage volumes are LOGNORMAL: most are trivial, a small
# tail is enormous, and the whole point of a fixed 300 ug dose is that it
# covers the body of the distribution but not the tail.  Median 0.2 mL of fetal
# whole blood at delivery, sigma 1.82 on the log scale, which puts ~1.6% above
# 10 mL and ~0.3% above 30 mL -- the classic figures.
FMH_MED = 0.20
FMH_SIG = 1.82
FMH_ANTE_FRAC = 0.07   # third-trimester silent FMH as a share of the total load


def fmh_draws(n=20000, seed=3):
    rng = np.random.default_rng(seed)
    return FMH_MED * np.exp(FMH_SIG * rng.standard_normal(n))


def sens_maternal(p, fmh_delivery, fmh_ante, rhig_ante, rhig_pp, delay_d=1.0,
                  days=260.0, dt=0.05):
    """Maternal-side-only integration of the prophylaxis race.  States:
    free fetal D+ RBC, coated fetal D+ RBC, passive anti-D, IM depot, signal.

    fmh volumes are mL of fetal WHOLE blood; antigen-bearing red cells are
    ~0.42 of that.  rhig doses are in ug.
    """
    free = coat = rh = dep = sig = 0.0
    t = 0.0
    # internal clock in days: the silent antenatal bleed at "28 weeks", then
    # delivery 9 weeks (63 d) later, then 6 months of follow-up
    t_ante, t_deliv = 30.0, 93.0
    ante_done = pp_done = False
    while t < days:
        if (not ante_done) and t >= t_ante:
            free += fmh_ante * 0.42
            if rhig_ante > 0:
                dep += rhig_ante * p["ad_iu_per_ug"]
            ante_done = True
        if (not pp_done) and t >= t_deliv:
            free += fmh_delivery * 0.42
            pp_done = True
        if rhig_pp > 0 and abs(t - (t_deliv + delay_d)) < dt / 2:
            dep += rhig_pp * p["ad_iu_per_ug"]
        need = p["iu_per_ml_rbc"] * (free + coat)
        coat_f = rh / (rh + need) if (rh + need) > 0 else 0.0
        to_coat = 4.0 * coat_f * free
        kel = p["k_int"] * (1.0 - p["phi_rescue"] / (1.0 + p["igg0"] / p["K_igg"]))
        d_free = -to_coat - p["k_sen_free"] * free
        d_coat = to_coat - p["k_clear_coat"] * coat
        d_rh = 0.35 * dep - kel * rh - p["k_clear_coat"] * coat * p["iu_per_ml_rbc"]
        d_dep = -0.35 * dep
        # exposure integral, not clearance flux
        ag_prime = free + (1.0 - p["dev"]) * coat
        free += dt * d_free
        coat += dt * d_coat
        rh = max(0.0, rh + dt * d_rh)
        dep += dt * d_dep
        sig += dt * p["ksens"] * ag_prime
        free = max(free, 0.0)
        coat = max(coat, 0.0)
        t += dt
    return 1.0 - math.exp(-sig)


def pop_risk(p, rhig_ante=0.0, rhig_pp=0.0, delay_d=1.0, draws=None):
    if draws is None:
        draws = fmh_draws()
    # thin the draws for speed but keep the tail: stratified sample
    q = np.quantile(draws, np.linspace(0.002, 0.998, 260))
    w = np.ones_like(q) / len(q)
    risk = np.array([sens_maternal(p, v, p["fmh_ante"], rhig_ante, rhig_pp,
                                   delay_d) for v in q])
    return float(np.sum(w * risk))


def stage3(p0):
    say("\n=== STAGE 3: sensitisation (ksens, dev) ===")
    say("targets: 16.0% with no prophylaxis, 1.6% with postpartum 300 ug")
    draws = fmh_draws()

    def resid(x):
        p = dict(p0)
        p["ksens"], p["fmh_ante"] = math.exp(x[0]), math.exp(x[1])
        a = pop_risk(p, 0.0, 0.0, draws=draws)
        b = pop_risk(p, 0.0, 300.0, draws=draws)
        r = [math.log(max(a, 1e-7) / 0.16), math.log(max(b, 1e-7) / 0.016)]
        say("   ksens=%.5g fmh_ante=%.4f mL -> none %.4f  pp %.5f  resid %.3f %.3f"
            % (p["ksens"], p["fmh_ante"], a, b, r[0], r[1]))
        return r

    sol = least_squares(resid, [math.log(0.02), math.log(0.02)], diff_step=0.10,
                        xtol=1e-4, ftol=1e-4)
    ksens, fmh_ante = math.exp(sol.x[0]), math.exp(sol.x[1])
    say("   FITTED ksens = %.6g /(mL.d), fmh_ante = %.5f mL" % (ksens, fmh_ante))
    return ksens, fmh_ante


# ==============================================================================
# STAGE 4 -- neonatal conjugation
# ==============================================================================
def healthy_neonate(p):
    c = M.Ctl()
    r = M.simulate(p, c, ga0=36.0, anti_d_iu=0.0, protocol="none",
                   ga_deliver=40.0, postnatal_days=18.0)
    m = r["born"] > 0
    return r["pna"][m], r["tsb"][m]


def stage4(p0):
    say("\n=== STAGE 4: neonatal conjugation (ugt_birth, ugt_t50) ===")
    say("targets: healthy term TSB peak 8.0 mg/dL at 4.0 d; 2.5 mg/dL at 14 d")

    def resid(x):
        p = dict(p0)
        p["ugt_birth"], p["ugt_t50"] = math.exp(x[0]), math.exp(x[1])
        pna, tsb = healthy_neonate(p)
        i = int(np.argmax(tsb))
        peak, tpeak = tsb[i], pna[i]
        d14 = float(np.interp(14.0, pna, tsb))
        r = [(peak - 8.0) / 1.0, (tpeak - 4.0) / 1.5, (d14 - 2.5) / 1.0]
        say("   ugt0=%.4f t50=%.2f -> peak %.2f at %.2f d, d14 %.2f" %
            (p["ugt_birth"], p["ugt_t50"], peak, tpeak, d14))
        return r

    sol = least_squares(resid, [math.log(0.075), math.log(25.5)],
                        diff_step=0.08, xtol=1e-3, ftol=1e-3)
    return math.exp(sol.x[0]), math.exp(sol.x[1])


# ==============================================================================
def main():
    p = dict(M.P)
    v0, g_pl = stage1()
    p["v0"], p["g_pl"] = v0, g_pl
    kops = stage2(p)
    p["kops"] = kops
    ksens, fmh_ante = stage3(p)
    p["ksens"], p["fmh_ante"] = ksens, fmh_ante
    ugt0, ugt_t50 = stage4(p)
    p["ugt_birth"], p["ugt_t50"] = ugt0, ugt_t50

    fitted = dict(v0=v0, g_pl=g_pl, kops=kops, ksens=ksens, fmh_ante=fmh_ante,
                  ugt_birth=ugt0, ugt_t50=ugt_t50)
    say("\n=== FITTED PARAMETER SET (7 numbers) ===")
    for k, v in fitted.items():
        say("   %-10s = %.6g" % (k, v))
    with open("hdfn_calibration.json", "w") as fh:
        json.dump(fitted, fh, indent=2)
    say("\nwritten to hdfn_calibration.json")


if __name__ == "__main__":
    main()
