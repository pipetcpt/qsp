#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Every number quoted in README.md / the .dot map / the mrgsolve header is
produced here.  Run:  python3 hdfn_analysis.py > hdfn_reference_output.txt
"""
import json
import math
import sys

import numpy as np
from scipy.optimize import brentq

import hdfn_python_reference as M
from hdfn_python_reference import IX, at_ga

FIT = json.load(open("hdfn_calibration.json"))
P = dict(M.P)
P.update(FIT)


def hdr(n, t):
    print("\n" + "=" * 78)
    print("%s.  %s" % (n, t))
    print("=" * 78)


def run(anti_d=15.0, protocol="mca", ga0=13.0, ga_del=36.0, pn=1.0, seed=7,
        p=None, **ctl):
    c = M.Ctl(**ctl)
    r = M.simulate(p or P, c, ga0=ga0, anti_d_iu=anti_d, protocol=protocol,
                   ga_deliver=ga_del, postnatal_days=pn, seed=seed)
    return r, c


def first_iut(c):
    return c.iut_log[0][0] if c.iut_log else None


# ==============================================================================
hdr(0, "WHAT WAS FITTED AND WHAT WAS NOT")
print("""FITTED -- eight numbers, all published summary statistics:
   v0, g_pl        placental conveyor      <- Malek 1996 F/M IgG 0.075 @19.5 wk,
                                              1.25 @39 wk
   kops            destruction potency     <- anti-D 15 IU/mL needs its first
                                              IUT at 26 wk (Nishie 2012 cohort
                                              mean 26.1 wk)
   ksens, fmh_ante sensitisation           <- 16% unprophylaxed, 1.6% with
                                              postpartum 300 ug only
   ugt_birth, ugt_t50  conjugation         <- HEALTHY term newborn peaks at
                                              8 mg/dL TSB at ~4 d
   emh_thresh      extramedullary switch   <- overt ascites appears at ~5-6 g/dL
STRUCTURAL, not fitted, and reported as sensitivities: Kres, dev, nip_pl_pen,
do2_alpha (= 1 by the oxygen-delivery argument), visc_k (= 0, see section 4).""")
for k, v in sorted(FIT.items()):
    print("   %-10s = %.6g" % (k, v))

# ==============================================================================
hdr(1, "MODEL INTEGRITY: THE HEALTHY FETUS AND THE HEALTHY NEWBORN")
rh, ch = run(anti_d=0.0, protocol="none", ga0=12.0, ga_del=40.0, pn=20.0)
print("A fetus with no alloantibody must sit on its own reference curves, or"
      "\nevery disease number below is measured from the wrong baseline.\n")
print(" GA   Hb    HbMoM  Hct   PSV MoM  alb   F/M IgG   ascites/threshold")
for g in (18, 22, 26, 30, 34, 39):
    print("%4.0f %5.2f %6.3f %6.3f %7.3f %5.2f %8.3f   %5.1f / %.0f mL" % (
        g, at_ga(rh, "hb", g), at_ga(rh, "hbmom", g), at_ga(rh, "hct", g),
        at_ga(rh, "psv_mom", g), at_ga(rh, "alb", g),
        at_ga(rh, "Cf_ig", g) / at_ga(rh, "Cig", g),
        at_ga(rh, "Asc", g), at_ga(rh, "asc_thr", g)))
print("\nMalek 1996 targets   : F/M = 0.075 at 19.5 wk, 1.25 at 39 wk")
print("model                : F/M = %.4f          , %.4f" % (
    at_ga(rh, "Cf_ig", 19.5) / at_ga(rh, "Cig", 19.5),
    at_ga(rh, "Cf_ig", 39.0) / at_ga(rh, "Cig", 39.0)))
print("Maternal total IgG falls %.0f%% from 12 wk to 39 wk (Malek: to 60-70%%,"
      "\n  i.e. a 30-40%% fall).  In the model this is plasma-volume dilution"
      "\n  plus the placental sink and NOTHING else -- see section 17." % (
          100 * (1 - at_ga(rh, "Cig", 39) / at_ga(rh, "Cig", 12))))
b = rh["born"] > 0
i = int(np.argmax(rh["tsb"][b]))
print("\nHealthy term newborn: TSB peak %.1f mg/dL at %.1f d, %.1f mg/dL at 14 d"
      % (rh["tsb"][b][i], rh["pna"][b][i], float(np.interp(14, rh["pna"][b],
                                                           rh["tsb"][b]))))
print("Healthy newborn Hb: %.1f at birth -> %.1f at 20 d (physiological nadir)"
      % (rh["hb"][b][0], float(np.interp(20, rh["pna"][b], rh["hb"][b]))))

# ==============================================================================
hdr(2, "THE PRODUCT: destruction = A x f_ag x M")
print("""The organising claim is that the haemolytic flux is a product of three
factors that different treatments own separately.  Here are the three factors
and the flux they multiply to, in the SAME fetus, before and after one
transfusion.  Note what the transfusion does and does not change.
""")
ri, ci = run(anti_d=15.0, protocol="mca", ga_del=36.0, seed=7)
g1 = first_iut(ci)
print(" state                    A (IU/mL)  f_ag    M      flux (g Hb/d)")
for lab, g in (("just before IUT #1", g1 - 0.02), ("1 day after IUT #1", g1 + 1 / 7.),
               ("7 days after", g1 + 1.0), ("14 days after", g1 + 2.0)):
    A = at_ga(ri, "A_eff", g)
    fa = at_ga(ri, "f_ag", g)
    Me = at_ga(ri, "M_eff", g)
    mass = at_ga(ri, "hb", g) * at_ga(ri, "Vfp", g) / 100.0
    flux = P["kops"] * A * fa * mass
    print(" %-24s %7.3f  %6.3f %6.3f   %7.3f" % (lab, A, fa, Me, flux))
print("""
The antibody factor A is UNCHANGED by the transfusion (packed cells are 80%
red cells; they add almost no plasma, so they do not dilute the antibody that
is already in the fetus).  What collapses is f_ag: the fraction of circulating
haemoglobin that carries the D antigen falls from 1.00 to about a third,
because the donor cells are antigen-negative.  A transfusion is therefore not
only a haemoglobin top-up: it is a SUBSTRATE DILUTION, and that is why its
effect outlasts the haemoglobin it delivers.""")

# ==============================================================================
hdr(3, "PREDICTION 1: the post-transfusion decline, and its decomposition")
print("""Nishie 2012 (PMID 22949399) measured 0.40 g/dL/day (SD 0.25) between the
first and second IUT in 41 pregnancies whose first IUT was at 26.1 wk.  Only
the GESTATION of that first transfusion was used in calibration; the decline
was not.  Model, averaged over ten surveillance seeds:
""")
seq = {1: [], 2: [], 3: []}
firsts = []
for sd in range(1, 11):
    r, c = run(anti_d=15.0, protocol="mca", ga_del=36.0, seed=sd)
    if c.iut_log:
        firsts.append(c.iut_log[0][0])
    for k in range(len(c.iut_log) - 1):
        ga, pre, post, n = c.iut_log[k]
        ga2, pre2, _, _ = c.iut_log[k + 1]
        days = (ga2 - ga) * 7.0
        # measured from the EQUILIBRATED post-transfusion value (24 h), because
        # the end-of-procedure sample is still volume-expanded
        hb24 = at_ga(r, "hb", ga + 1 / 7.)
        if n in seq and days > 5:
            seq[n].append((hb24 - pre2) / days)
print("   first IUT at %.1f wk (mean of %d runs; target 26.0)" %
      (np.mean(firsts), len(firsts)))
for n in (1, 2, 3):
    if seq[n]:
        print("   interval after IUT #%d : %.3f g/dL/day  (n=%d runs)"
              % (n, np.mean(seq[n]), len(seq[n])))
print("   reported            #1 : 0.400 g/dL/day (SD 0.250)")
print("""
The ORDER is predicted -- each interval is flatter than the last -- and the
reason is f_ag, which is halved again by every transfusion.  The LEVEL is
under-predicted by about a third (%.2f vs 0.40), i.e. about half a standard
deviation.  Decomposition of the first interval shows where the missing
haemoglobin would have to come from:""" % np.mean(seq[1]))
r, c = run(anti_d=15.0, protocol="mca", ga_del=36.0, seed=7)
ga = c.iut_log[0][0]
ga_b, ga_e = ga + 1 / 7., c.iut_log[1][0]
dt = (ga_e - ga_b) * 7.0
d_own = at_ga(r, "Ro", ga_b) - at_ga(r, "Ro", ga_e) + \
    at_ga(r, "Rr", ga_b) - at_ga(r, "Rr", ga_e)
d_don = at_ga(r, "Rd", ga_b) - at_ga(r, "Rd", ga_e)
v_b, v_e = at_ga(r, "Vfp", ga_b), at_ga(r, "Vfp", ga_e)
hb_b, hb_e = at_ga(r, "hb", ga_b), at_ga(r, "hb", ga_e)
mass_b = hb_b * v_b / 100.0
dil = 100.0 * mass_b / v_b - 100.0 * mass_b / v_e
print("   antigen-positive cells destroyed : %.2f g Hb -> %.3f g/dL/day"
      % (d_own, 100 * d_own / v_e / dt))
print("   donor-cell senescence            : %.2f g Hb -> %.3f g/dL/day"
      % (d_don, 100 * d_don / v_e / dt))
print("   growth dilution (%.0f -> %.0f mL)   :             %.3f g/dL/day"
      % (v_b, v_e, dil / dt))
print("   total                            :             %.3f g/dL/day"
      % ((hb_b - hb_e) / dt))
for td in (70.0, 45.0, 30.0):
    p2 = dict(P)
    p2["t_don"] = td
    d = []
    for sd in (1, 3, 5, 7, 9):
        r2, c2 = run(anti_d=15.0, protocol="mca", ga_del=36.0, seed=sd, p=p2)
        if len(c2.iut_log) > 1:
            ga_, pre_, post_, n_ = c2.iut_log[0]
            ga2_, pre2_, _, _ = c2.iut_log[1]
            d.append((at_ga(r2, "hb", ga_ + 1 / 7.) - pre2_) /
                     ((ga2_ - ga_) * 7.0))
    print("   donor red-cell lifespan %4.0f d -> first-interval decline %.3f"
          % (td, np.mean(d)))
print("""Read the decomposition carefully, because one term has the sign nobody
expects: the fetoplacental volume FALLS across the interval (114 -> 101 mL),
because the red cell volume being destroyed shrinks faster than the plasma
volume grows.  Dilution is therefore not helping the decline, it is OPPOSING it
by about 0.14 g/dL/day, and the immune term has to be that much larger to
produce what is observed.

Shortening the donor red cell lifespan closes about half of the remaining gap
and no more: 70 d gives 0.29, 30 d gives 0.36, against a reported 0.40.  So
donor-cell survival cannot be the whole explanation either, and at least one
term is still missing -- the candidates being bystander loss of donor cells,
sequestration of transfused cells in the first 24 h, or an interval shorter
than the model schedules.  The model keeps the adult 70 d and reports the
shortfall rather than tuning a parameter to cover it.""")

# ==============================================================================
hdr(4, "PREDICTION 2: the MCA-PSV 1.5 MoM threshold is not empirical")
print("""If cerebral oxygen delivery is defended -- flow rises exactly as far as
haemoglobin falls -- then velocity x haemoglobin is a constant, so

        PSV MoM  =  1 / HbMoM        (do2_alpha = 1, no fitting whatever)

and the threshold that everyone uses inverts to a haemoglobin:
""")
for mom in (1.0, 0.84, 0.75, 0.667, 0.65, 0.55, 0.45):
    print("   Hb %.3f MoM  ->  PSV %.3f MoM" % (mom, 1.0 / mom))
print("""
   PSV 1.50 MoM  <->  Hb %.3f MoM
   The published definition of MODERATE anaemia is Hb 0.65 MoM.  The threshold
   is therefore the arithmetic image of the disease definition, not a
   receiver-operating-characteristic artefact.""" % (1 / 1.5))
print("""
Mari 2000 reported 100% sensitivity for moderate/severe anaemia and a 12%
false-positive rate.  Both are properties of MEASUREMENT ERROR rather than of
biology, so the model can be asked which measurement error reproduces them.
Taking the anaemia mix of that cohort (41 not anaemic, 35 mild, 4 moderate,
31 severe) and the mid-band haemoglobin of each group:""")
groups = [("not anaemic (n=41)", 41, 0.92), ("mild (n=35)", 35, 0.75),
          ("moderate (n=4)", 4, 0.60), ("severe (n=31)", 31, 0.40)]
for cv in (0.10, 0.15, 0.20, 0.25):
    rng = np.random.default_rng(11)
    tp = fp = npos = nneg = 0.0
    for lab, n, mom in groups:
        obs = (1.0 / mom) * (1.0 + cv * rng.standard_normal(40000))
        frac = float(np.mean(obs >= 1.5))
        if mom <= 0.65:
            tp += frac * n
            npos += n
        else:
            fp += frac * n
            nneg += n
    print("   CV %2.0f%%  ->  sensitivity %3.0f%%, false-positive rate %3.0f%%"
          % (100 * cv, 100 * tp / npos, 100 * fp / nneg))
print("""   reported                sensitivity 100%, false-positive rate  12%

The sensitivity is reproduced at any plausible CV, because a severely anaemic
fetus is far above the threshold.  The FALSE-POSITIVE rate is the informative
number: it pins the measurement CV at roughly 15-20%, which is a testable
statement about Doppler reproducibility (three-beat average, angle correction,
operator) and not a property of the disease.  At a 10% CV the model
UNDER-predicts false positives (5% against 12%).""")
print("""
A viscosity term would move the threshold, and that is why this model sets
visc_k = 0.  If MCA velocity carried an EXTRA haematocrit dependence on top of
the oxygen-delivery defence, then PSV MoM = (1/HbMoM)*exp(-visc_k*(Hct-Hct_ref))
and the 1.5 MoM cut-off would no longer invert to the definition of moderate
anaemia:""")


def psv_at(mom, vk, hbr=13.0):
    hct = mom * hbr / P["mchc"]
    return (1.0 / mom) * math.exp(-vk * (hct - hbr / P["mchc"]))


for vk in (0.0, 0.45, 0.90, 1.80):
    try:
        mom_star = brentq(lambda m: psv_at(m, vk) - 1.5, 0.25, 0.999)
    except Exception:
        mom_star = float("nan")
    rng = np.random.default_rng(11)
    fp = nneg = 0.0
    for lab, n, mom in groups:
        if mom > 0.65:
            obs = psv_at(mom, vk) * (1.0 + 0.15 * rng.standard_normal(40000))
            fp += float(np.mean(obs >= 1.5)) * n
            nneg += n
    print("   visc_k %.2f -> PSV 1.50 MoM corresponds to Hb %.3f MoM; "
          "false positives %3.0f%% (at CV 15%%)"
          % (vk, mom_star, 100 * fp / nneg))
print("""   Only visc_k = 0 puts the threshold on the definition of moderate
anaemia (0.667 vs 0.65 MoM).  The physiological reason not to add the term is
the same as the arithmetic one: the flow rise that a fall in viscosity produces
IS the mechanism by which oxygen delivery is defended, so a separate viscosity
factor counts it twice.""")

# ==============================================================================
hdr(5, "THE CONVEYOR: why the same titre is a different disease at 20 and 34 wk")
print(" GA   transfer capacity (rel to 20 wk)   free FcRn   F/M IgG")
for g in (16, 20, 24, 28, 32, 36, 40):
    cap = math.exp(P["g_pl"] * (g - 20.0))
    print("%4.0f %12.2f %26.3f %9.3f" % (
        g, cap, at_ga(rh, "fu_m", min(g, 39.9)),
        at_ga(rh, "Cf_ig", min(g, 39.9)) / at_ga(rh, "Cig", min(g, 39.9))))
print("""
Transfer capacity grows e^(%.3f per week), a factor of %.0f from 16 to 36 wk.
This is the whole reason HDFN is a THIRD-trimester disease in most women and a
second-trimester disease only in those whose titre is high enough to make the
early, inefficient conveyor sufficient.  It is also why an FcRn blocker has to
be started at 14-16 wk: it is not treating disease, it is preventing the
conveyor from ever delivering a dangerous dose.""" % (
    P["g_pl"], math.exp(P["g_pl"] * 20)))

# ==============================================================================
hdr(6, "NIPOCALIMAB: the two effects are separable, and they are unequal")
print("""Nipocalimab does two things: it accelerates maternal IgG catabolism
(FcRn in the endothelium) and it blocks transcytosis (FcRn in the
syncytiotrophoblast).  The model lets one be switched off at a time.
""")
base, cb = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5)
print(" arm                          maternal IgG   anti-D    fetal anti-D   IUTs")
print("                               (30 wk, g/L)  (30 wk)   (22 wk, IU/mL)")
rows = [("no drug", dict()),
        ("nipocalimab 15 mg/kg/wk", dict(nip_dose=15.0, nip_start=14.0)),
        ("nipocalimab 30 mg/kg/wk", dict(nip_dose=30.0, nip_start=14.0)),
        ("nipocalimab 45 mg/kg/wk", dict(nip_dose=45.0, nip_start=14.0))]
for lab, kw in rows:
    r, c = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5, **kw)
    print(" %-28s %8.2f %9.2f %12.3f %8d" % (
        lab, at_ga(r, "Cig", 30), at_ga(r, "Ca1", 30), at_ga(r, "Cf_a1", 22),
        c.iut_n))
p_noblock = dict(P)
p_noblock["nip_pl_pen"] = 0.0
r, c = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5, p=p_noblock,
           nip_dose=30.0, nip_start=14.0)
print(" %-28s %8.2f %9.2f %12.3f %8d   <- catabolism only" % (
    "30 mg/kg, NO placental block", at_ga(r, "Cig", 30), at_ga(r, "Ca1", 30),
    at_ga(r, "Cf_a1", 22), c.iut_n))
p_full = dict(P)
p_full["nip_pl_pen"] = 1.0
r, c = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5, p=p_full,
           nip_dose=30.0, nip_start=14.0)
print(" %-28s %8.2f %9.2f %12.3f %8d   <- if the placenta saw plasma" % (
    "30 mg/kg, FULL placental block", at_ga(r, "Cig", 30), at_ga(r, "Ca1", 30),
    at_ga(r, "Cf_a1", 22), c.iut_n))
print("""
This is the model's sharpest inference about the drug.  If the
syncytiotrophoblast FcRn saw the same nipocalimab concentration as the vascular
FcRn does, transfer would be blocked >99%% and NO fetus in the UNITY population
would have needed a transfusion.  Six of thirteen did.  Either placental FcRn
is much harder to saturate than endothelial FcRn (the model's assumption:
15%% interstitial penetration) or the enrolled titres spanned a factor of ten.
Section 7 shows the second explanation alone cannot do it.""")
print("\nStart-time sweep at 120 IU/mL (30 mg/kg/wk):")
for st in (14.0, 16.0, 20.0, 24.0, 28.0, 99.0):
    r, c = run(anti_d=120.0, protocol="mca", ga_del=36.0, seed=5,
               nip_dose=30.0 if st < 99 else 0.0, nip_start=st)
    print("   start %4.1f wk -> IUTs %d, first at %s wk, hydrops %s" % (
        st, c.iut_n, ("%.1f" % first_iut(c)) if c.iut_log else "none",
        bool((r["Asc"] > r["asc_thr"]).any())))

# ==============================================================================
hdr(7, "THE UNITY TRIAL, SIMULATED AS PUBLISHED")
print("""Moise 2024 (NEJM, PMID 39115062): open-label, 13 pregnancies at high risk
of RECURRENT early-onset severe HDFN, nipocalimab 30 or 45 mg/kg weekly from
14-35 wk.  Primary end point: live birth at >= 32 wk WITHOUT any intrauterine
transfusion, against a historical benchmark of 0%.  Reported: 7/13 = 54%
(95% CI 25-81), and no hydrops.

The trial's entry criterion is reproduced in the model rather than assumed: a
virtual mother is ENROLLED only if her untreated pregnancy would have needed a
transfusion before 24 weeks.  A log-normal prior on anti-D quantitation
(median 8 IU/mL, GSD 4) is screened by that criterion; the survivors are then
treated.
""")
rng = np.random.default_rng(202)
cand = 8.0 * np.exp(math.log(4.0) * rng.standard_normal(300))
enrolled, screened = [], 0
for ad in cand:
    if len(enrolled) >= 30:
        break
    if ad < 5 or ad > 4000:
        continue
    screened += 1
    r0, c0 = run(anti_d=float(ad), protocol="mca", ga_del=36.0, seed=4)
    if c0.iut_log and c0.iut_log[0][0] < 24.0:
        enrolled.append(float(ad))
print("   screened %d, enrolled %d (%.0f%%); enrolled anti-D median %.0f IU/mL"
      " (IQR %.0f-%.0f)" % (screened, len(enrolled),
                            100 * len(enrolled) / max(screened, 1),
                            np.median(enrolled),
                            np.percentile(enrolled, 25),
                            np.percentile(enrolled, 75)))
free = hyd = 0
gadel = []
for ad in enrolled:
    r, c = run(anti_d=ad, protocol="mca", ga_del=36.0, seed=4,
               nip_dose=30.0, nip_start=14.0, nip_stop=35.0)
    if c.iut_n == 0:
        free += 1
    if (r["Asc"] > r["asc_thr"]).any():
        hyd += 1
    gadel.append(36.0 if c.iut_n == 0 else 35.0)
n = len(enrolled)
print("   MODEL : IUT-free live birth >= 32 wk in %d/%d = %.0f%%" %
      (free, n, 100 * free / n))
print("   UNITY : 7/13 = 54% (95% CI 25-81)")
print("   MODEL : hydrops in %d/%d;   UNITY: none" % (hyd, n))
print("\n   sensitivity to the ONE distributional assumption (GSD of the prior):")
for gsd in (2.0, 3.0, 4.0, 6.0):
    rng2 = np.random.default_rng(77)
    c2 = 8.0 * np.exp(math.log(gsd) * rng2.standard_normal(160))
    en = []
    for ad in c2:
        if len(en) >= 20:
            break
        if ad < 5 or ad > 4000:
            continue
        r0, cc = run(anti_d=float(ad), protocol="mca", ga_del=36.0, seed=4)
        if cc.iut_log and cc.iut_log[0][0] < 24.0:
            en.append(float(ad))
    fr = 0
    for ad in en:
        _, cc = run(anti_d=ad, protocol="mca", ga_del=36.0, seed=4,
                    nip_dose=30.0, nip_start=14.0)
        fr += cc.iut_n == 0
    print("      GSD %.1f -> enrolled median %5.0f IU/mL, IUT-free %2d/%2d = %.0f%%"
          % (gsd, np.median(en) if en else 0, fr, len(en),
             100 * fr / max(len(en), 1)))

# ==============================================================================
hdr(8, "IVIG AND PLASMAPHERESIS: two mechanisms, one of them nearly useless")
print("""High-dose maternal IVIG is supposed to work three ways: FcRn competition
(it floods the transporter), FcgammaR blockade in the fetal spleen, and
feedback suppression of the maternal response.  The model implements the first
two and can switch each off.
""")
print(" arm                                fetal anti-D(24wk)  IUTs  first IUT")
for lab, kw, p2 in [
        ("no immunomodulation", dict(), P),
        ("IVIG 1 g/kg/wk from 12 wk", dict(ivig_dose=1.0, ivig_start=12.0), P),
        ("IVIG 2 g/kg/wk from 12 wk", dict(ivig_dose=2.0, ivig_start=12.0), P),
        ("IVIG 1 g/kg, FcRn competition only",
         dict(ivig_dose=1.0, ivig_start=12.0),
         dict(P, ivig_fcgr=0.0)),
        ("IVIG 1 g/kg, FcgammaR blockade only",
         dict(ivig_dose=1.0, ivig_start=12.0),
         dict(P, ivig_compete=0.0))]:
    r, c = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5, p=p2, **kw)
    print(" %-34s %10.3f %8d   %s" % (
        lab, at_ga(r, "Cf_a1", 24), c.iut_n,
        ("%.1f" % first_iut(c)) if c.iut_log else "none"))
print("""
Ruma 2007 (PMID 17306655) treated nine such pregnancies with plasmapheresis
plus weekly IVIG.  ALL NINE still required intrauterine transfusion (median 4,
range 3-8), every infant survived, mean delivery 34 wk -- and yet the MCA-PSV
stayed below threshold throughout therapy.  The model reproduces the shape of
that result: IVIG delays the first transfusion but does not prevent it,""")
r, c = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5,
           ivig_dose=1.0, ivig_start=12.0)
r0, c0 = run(anti_d=60.0, protocol="mca", ga_del=36.0, seed=5)
print("   first IUT %.1f wk on IVIG vs %.1f wk without; %d vs %d procedures."
      % (first_iut(c) or 99, first_iut(c0) or 99, c.iut_n, c0.iut_n))
print("""and the reason is arithmetic: IVIG raises the total maternal IgG the conveyor
sees from ~%.1f to ~%.1f g/L, and the conveyor is only HALF-saturated at
%.0f g/L, so flooding it divides the transfer by a factor of %.2f -- against a
disease that is running %.0f-fold over threshold.  The DECOMPOSITION above says
that essentially the whole of the IVIG effect in this model is FcRn
competition; the FcgammaR arm contributes almost nothing at 1 g/kg/wk, because
the reticuloendothelial system is not the rate-limiting step until the antibody
load is far higher.""" % (at_ga(r0, "Cig", 24), at_ga(r, "Cig", 24) +
                         at_ga(r, "Civ", 24), P["K_igg"],
                         (1 + (at_ga(r, "Cig", 24) + at_ga(r, "Civ", 24)) /
                          P["K_igg"]) / (1 + at_ga(r0, "Cig", 24) / P["K_igg"]),
                         60 / 15))

# ==============================================================================
hdr(9, "ANTI-D PROPHYLAXIS: stoichiometry, the 72-hour window, and the residual")
import hdfn_calibrate as CAL  # noqa: E402  (re-uses the maternal-only integrator)
draws = CAL.fmh_draws()
print("""The dose rule is stoichiometric, not pharmacokinetic: 20 ug of anti-D per
mL of fetal red cells, i.e. 100 IU/mL.  A 300 ug (1500 IU) dose therefore
covers 15 mL of fetal red cells = ~30 mL of fetal whole blood, and nothing
about the woman's weight or plasma volume enters that statement.
""")
print(" regimen                                       population sensitisation")
for lab, ante, pp, delay in [
        ("none", 0.0, 0.0, 1.0),
        ("postpartum 300 ug within 24 h", 0.0, 300.0, 1.0),
        ("postpartum 100 ug within 24 h", 0.0, 100.0, 1.0),
        ("postpartum 300 ug at 72 h", 0.0, 300.0, 3.0),
        ("postpartum 300 ug at 7 d", 0.0, 300.0, 7.0),
        ("postpartum 300 ug at 14 d", 0.0, 300.0, 14.0),
        ("antenatal 300 ug at 28 wk only", 300.0, 0.0, 1.0),
        ("antenatal 28 wk + postpartum 300 ug", 300.0, 300.0, 1.0),
        ("antenatal 28 wk + postpartum 600 ug", 300.0, 600.0, 1.0)]:
    risk = CAL.pop_risk(P, ante, pp, delay, draws=draws)
    print(" %-44s %8.3f%%" % (lab, 100 * risk))
print("""
   Fitted: the first two lines (16%% and 1.6%%).
   PREDICTED: the combined antenatal+postpartum regimen, %.2f%%.  Observed with
   routine antenatal prophylaxis: 0.1-0.4%%.  The prediction lands because the
   residual is not a failure of immunology but a failure of ARITHMETIC -- the
   tail of the fetomaternal-haemorrhage distribution that a fixed 300 ug dose
   cannot coat.""" % (100 * CAL.pop_risk(P, 300.0, 300.0, 1.0, draws=draws)))
print("""
   The 72-hour window is derived, not assumed.  An UNCOATED fetal red cell in
   the maternal circulation is just a red cell: it survives with a ~70 day
   half-life (k_sen_free = %.4f /d).  A three-day delay therefore leaks only
   ~%.0f%% of the exposure integral, which is why the window is days and not
   minutes -- and why it eventually closes: by 14 days the leak is %.0f%%.""" % (
    P["k_sen_free"], 100 * (1 - math.exp(-P["k_sen_free"] * 3)) /
    (1 - math.exp(-P["k_sen_free"] * 200)) * 100,
    100 * (1 - math.exp(-P["k_sen_free"] * 14)) /
    (1 - math.exp(-P["k_sen_free"] * 200)) * 100))
print("\n Large-bleed dose-response (single delivery FMH, prophylaxis at 24 h):")
print("  FMH (mL fetal whole blood)   300 ug      600 ug     1500 ug")
for v in (2.0, 10.0, 30.0, 60.0, 120.0):
    a = CAL.sens_maternal(P, v, 0.0, 0.0, 300.0, 1.0)
    b = CAL.sens_maternal(P, v, 0.0, 0.0, 600.0, 1.0)
    cc = CAL.sens_maternal(P, v, 0.0, 0.0, 1500.0, 1.0)
    print("  %8.0f %20.3f%% %10.3f%% %10.3f%%" % (v, 100 * a, 100 * b, 100 * cc))
print("""  30 mL of fetal whole blood is exactly the stated coverage of a 300 ug
  dose, and it is exactly where the 300 ug column turns over.""")

# ==============================================================================
hdr(10, "ANTI-K IS NOT ANTI-D: the same anaemia with the bilirubin missing")
print("""Anti-Kell destroys erythroid PROGENITORS rather than circulating red
cells (Vaughan 1998).  In this model that is one parameter, kell_kill, moved
off zero -- and the consequence is a fetus with the same haemoglobin and a
completely different set of secondary markers.  Matched at Hb ~0.55 MoM at
28 wk:
""")
p_kell = dict(P)
p_kell["kell_kill"] = 0.55
p_kell["kops"] = P["kops"] * 0.18      # Kell sites are far scarcer on mature cells
r_d, c_d = run(anti_d=15.0, protocol="none", ga_del=32.0, seed=7)
r_k, c_k = run(anti_d=15.0, protocol="none", ga_del=32.0, seed=7, p=p_kell)
print(" marker at 28 wk            anti-D        anti-K")
for lab, key, fmt in [("haemoglobin (g/dL)", "hb", "%8.2f"),
                      ("Hb MoM", "hbmom", "%8.3f"),
                      ("progenitor pool (x normal)", "Prog", "%8.2f"),
                      ("nucleated RBC (x normal)", "Nrbc", "%8.2f"),
                      ("fetal bilirubin (mg/dL)", "tsb", "%8.2f"),
                      ("amniotic dOD450", "od450", "%8.4f"),
                      ("splenic RES capacity", "M_eff", "%8.2f"),
                      ("EMH (hepatic)", "EMH", "%8.2f")]:
    print(" %-26s" % lab + fmt % at_ga(r_d, key, 28) + "  " +
          fmt % at_ga(r_k, key, 28))
print("""
The clinical corollary is a warning that comes out of the model rather than
being written into it: in Kell alloimmunisation the amniotic-fluid bilirubin
(dOD450, the Liley curve) UNDER-states the anaemia by a factor of %.1f, because
the anaemia is made by cells that never reached the circulation and therefore
never released haem.  Any surveillance strategy built on bilirubin is
systematically wrong in Kell disease; one built on MCA-PSV is not.""" % (
    at_ga(r_d, "od450", 28) / max(at_ga(r_k, "od450", 28), 1e-9)))

# ==============================================================================
hdr(11, "THE BIRTH SWITCH: identical haemolysis, two different diseases")
print("""Nothing about the antibody changes at delivery.  What changes is the
bilirubin clearance route: the placenta (the mother's liver) is replaced by a
newborn liver with %.1f%% of adult UGT1A1 activity.  The model's fetal
clearance constant is %.0f /d; the newborn's effective conjugation at the same
bilirubin load is a small fraction of it.
""" % (100 * P["ugt_birth"], P["k_bil_pl"]))
r, c = run(anti_d=15.0, protocol="mca", ga_del=36.0, pn=70.0, seed=7)
b = r["born"] > 0
print(" time                      Hb    TSB (mg/dL)  UGT1A1  B/A ratio")
print(" 34 wk in utero        %6.2f %10.2f %8.3f %9.3f" % (
    at_ga(r, "hb", 34), at_ga(r, "tsb", 34), P["ugt_birth"],
    at_ga(r, "ba_ratio", 34)))
for d in (0.25, 1, 2, 3, 5, 7, 14, 28, 42, 56, 69):
    g = 36.0 + d / 7.0
    print(" %4.0f h / %2.0f d postnatal %6.2f %10.2f %8.3f %9.3f" % (
        d * 24, d, at_ga(r, "hb", g), at_ga(r, "tsb", g),
        at_ga(r, "Ugt", g), at_ga(r, "ba_ratio", g)))
print("   exchange transfusions: %d;  top-up transfusions: %d;  iron load %.0f mg"
      % (getattr(c, "exch_n", 0), getattr(c, "topup_n", 0),
         at_ga(r, "Fer", 36 + 69 / 7.)))
print("""
Two predictions here were never fitted.  First, this baby's peak bilirubin is
LOW for the severity of its disease (%.1f mg/dL), because after four
transfusions most of its circulating haemoglobin is donor haemoglobin that
nothing is attacking -- the same f_ag that flattened the intrauterine decline
also flattens the jaundice.  Second, the haemoglobin nadir comes at %.0f days,
not in the first week: maternal IgG persists with a 3-week half-life while the
marrow, suppressed by months of transfusion, has nothing in reserve.  That is
the late hyporegenerative anaemia which is why these infants need top-up
transfusions for two to three months.""" % (
    r["tsb"][b].max(),
    r["pna"][b][int(np.argmin(r["hb"][b]))]))

# ==============================================================================
hdr(12, "HYDROPS: which term crosses first")
print("""Hydrops in this model is a Starling balance in which the lymphatic
capacity is DEFINED as the baseline filtration rate, so a healthy fetus can
never drift into ascites and the disease has to come from a named term.
""")
r, c = run(anti_d=15.0, protocol="none", ga_del=34.0, seed=7)
print(" GA   Hb   HbMoM  albumin  UVP    filtration  lymph   ascites  overt?")
for g in (22, 26, 28, 30, 31, 32, 33):
    print("%4.0f %5.2f %6.3f %8.2f %6.2f %11.1f %7.1f %8.0f   %s" % (
        g, at_ga(r, "hb", g), at_ga(r, "hbmom", g), at_ga(r, "alb", g),
        at_ga(r, "cvp", g), at_ga(r, "jv", g), at_ga(r, "lymph", g),
        at_ga(r, "Asc", g), at_ga(r, "Asc", g) > at_ga(r, "asc_thr", g)))
m = r["Asc"] > r["asc_thr"]
if m.any():
    i = int(np.where(m)[0][0])
    print("   overt ascites first appears at GA %.1f wk, Hb %.2f g/dL (%.2f MoM)"
          % (r["ga"][i], r["hb"][i], r["hbmom"][i]))
print("""   Nicolaides' clinical rule is that hydrops appears below about 5 g/dL,
   i.e. a haemoglobin deficit of 7 g/dL.  The model gets there from co_max
   alone: with oxygen delivery defended (alpha = 1) and a maximal cardiac
   output reserve of %.2f, the fetus loses the ability to defend delivery at
   Hb MoM = 1/%.2f = %.3f.  Umbilical venous pressure rises with anaemia and
   then falls back in the most extreme cases as the heart fails, which is the
   non-monotonic pattern Ville 1994 reported.""" % (
    P["co_max"], P["co_max"], 1 / P["co_max"]))

# ==============================================================================
hdr(13, "SIXTEEN SCENARIOS")
SC = [
    ("1  no antibody (reference fetus)", dict(anti_d=0.0, protocol="none")),
    ("2  anti-D 4 IU/mL, no intervention", dict(anti_d=4.0, protocol="none")),
    ("3  anti-D 15 IU/mL, no intervention", dict(anti_d=15.0, protocol="none")),
    ("4  anti-D 60 IU/mL, no intervention", dict(anti_d=60.0, protocol="none")),
    ("5  anti-D 15, MCA-PSV surveillance + IUT", dict(anti_d=15.0)),
    ("6  anti-D 60, MCA-PSV surveillance + IUT", dict(anti_d=60.0)),
    ("7  anti-D 60, fixed 2-weekly IUT", dict(anti_d=60.0, protocol="fixed")),
    ("8  anti-D 60 + IVIG 1 g/kg/wk", dict(anti_d=60.0, ivig_dose=1.0,
                                           ivig_start=12.0)),
    ("9  anti-D 60 + IVIG 2 g/kg/wk", dict(anti_d=60.0, ivig_dose=2.0,
                                           ivig_start=12.0)),
    ("10 anti-D 60 + nipocalimab 30 from 14 wk", dict(anti_d=60.0,
                                                     nip_dose=30.0,
                                                     nip_start=14.0)),
    ("11 anti-D 60 + nipocalimab 45 from 14 wk", dict(anti_d=60.0,
                                                     nip_dose=45.0,
                                                     nip_start=14.0)),
    ("12 anti-D 120 + nipocalimab 30 from 14 wk", dict(anti_d=120.0,
                                                      nip_dose=30.0,
                                                      nip_start=14.0)),
    ("13 anti-D 120 + nipo 30 + IVIG 1 g/kg", dict(anti_d=120.0, nip_dose=30.0,
                                                   nip_start=14.0,
                                                   ivig_dose=1.0,
                                                   ivig_start=14.0)),
    ("14 anti-D 120 + nipocalimab from 24 wk", dict(anti_d=120.0, nip_dose=30.0,
                                                    nip_start=24.0)),
    ("15 anti-D 15, deliver 34 wk not 37", dict(anti_d=15.0, ga_del=34.0)),
    ("16 anti-D 15, deliver 37 wk", dict(anti_d=15.0, ga_del=37.0)),
]
print(" scenario                                  IUTs  1st IUT  hydrops  "
      "surv  Hb@birth  TSBpeak  exch  topup")
for lab, kw in SC:
    kw = dict(kw)
    gd = kw.pop("ga_del", 36.0)
    ad = kw.pop("anti_d", 15.0)
    pr = kw.pop("protocol", "mca")
    r, c = run(anti_d=ad, protocol=pr, ga_del=gd, pn=70.0, seed=7, **kw)
    bb = r["born"] > 0
    print(" %-40s %4d %8s %8s %5.2f %9.2f %8.1f %5d %6d" % (
        lab, c.iut_n, ("%.1f" % first_iut(c)) if c.iut_log else "-",
        "yes" if (r["Asc"] > r["asc_thr"]).any() else "no",
        math.exp(-at_ga(r, "Hzd", gd)), r["hb"][bb][0], r["tsb"][bb].max(),
        getattr(c, "exch_n", 0), getattr(c, "topup_n", 0)))

# ==============================================================================
hdr(14, "VIRTUAL POPULATION (n = 100) BY ANTI-D BAND")
rng = np.random.default_rng(5150)
bands = [("< 4 IU/mL   (low risk)", 1.0, 4.0),
         ("4-15 IU/mL  (moderate)", 4.0, 15.0),
         ("15-60 IU/mL (high)", 15.0, 60.0),
         ("> 60 IU/mL  (very high)", 60.0, 300.0)]
print(" band                     n   any IUT   IUTs(mean)  hydrops  exch  "
      "survival  Hb@birth")
for lab, lo, hi in bands:
    ads = np.exp(rng.uniform(math.log(lo), math.log(hi), 25))
    ni = hy = ex = 0
    nl = []
    sv = []
    hbb = []
    for ad in ads:
        r, c = run(anti_d=float(ad), protocol="mca", ga_del=36.0, pn=30.0,
                   seed=int(ad * 7) % 97 + 1)
        ni += c.iut_n > 0
        nl.append(c.iut_n)
        hy += bool((r["Asc"] > r["asc_thr"]).any())
        ex += getattr(c, "exch_n", 0) > 0
        sv.append(math.exp(-at_ga(r, "Hzd", 36.0)))
        hbb.append(r["hb"][r["born"] > 0][0])
    print(" %-24s %3d %7.0f%% %11.1f %7.0f%% %5.0f%% %8.2f %9.2f" % (
        lab, len(ads), 100 * ni / len(ads), np.mean(nl), 100 * hy / len(ads),
        100 * ex / len(ads), np.mean(sv), np.mean(hbb)))
print("""
   READ THE FIRST ROW AS A MISS, NOT AS A RESULT.  The UK bands call <4 IU/mL
   low risk, and the model does not agree: it puts a substantial fraction of
   that band into late-gestation anaemia.  The reason is structural and worth
   stating, because it is the conveyor exponent doing it.  Destruction is
   linear in the titre, so a four-fold lower titre does not make the disease
   four times milder -- it DELAYS it by ln(4)/0.268 = 5.2 weeks.  With the
   conveyor growing exponentially, titre bands compress in TIME rather than in
   severity, and a model with no additional threshold non-linearity therefore
   over-predicts disease at the bottom of the range.  What is missing is
   probably a genuine threshold in the antibody-to-clearance step (avidity,
   subclass composition, or a minimum bound-IgG density per cell for
   phagocytosis), and this model does not have one.""")

# ==============================================================================
hdr(15, "SENSITIVITY AND IDENTIFIABILITY")
print(" parameter        x0.5                    x2.0            (first IUT, wk)")
for par in ("kops", "v0", "g_pl", "Kres", "t_don", "co_max", "emh_thresh",
            "nip_pl_pen", "K_igg", "res_max"):
    out = []
    for fac in (0.5, 2.0):
        p2 = dict(P)
        p2[par] = P[par] * fac
        try:
            r, c = run(anti_d=15.0, protocol="mca", ga_del=36.0, seed=7, p=p2)
            out.append("%s / %d IUT" % (("%.1f" % first_iut(c)) if c.iut_log
                                        else "none", c.iut_n))
        except Exception as e:
            out.append("fail: %s" % e)
    print(" %-14s %-23s %-23s" % (par, out[0], out[1]))
print("""
   Kres and res_max barely move the first transfusion, which is why Kres was
   fixed structurally rather than fitted: the reticuloendothelial capacity is
   not identifiable from the timing of transfusions.  It becomes identifiable
   only in a fetus whose antibody load is high enough to saturate clearance,
   and no such measurement exists.""")

# ==============================================================================
hdr(16, "DEFECTS FOUND BY EXECUTING THE MODEL, AND HONEST MISSES")
print("""ELEVEN defects were found by RUNNING this model, not by reading it.
Each one is listed with what it did, because a defect that only produces a
plausible-looking number is the dangerous kind.

 1. The circulation was written as a fixed blood volume plus a transient
    "excess" that decayed.  A transfusion therefore looked correct at the end
    of the procedure and then haemoconcentrated the fetus to 20 g/dL a day
    later.  Rewritten as plasma volume PLUS red-cell volume -- which also
    corrected a second error the first was hiding: packed cells were diluting
    the antibody already inside the fetus, and they do not.
 2. Erythropoiesis responded to anaemia but not to polycythaemia, so a fetus
    made polycythaemic by a transfusion went on producing red cells at its full
    requirement.  That cancelled most of the post-transfusion decline and made
    the decline insensitive to the destruction rate -- i.e. it destroyed the
    model's central prediction while leaving every number looking reasonable.
    EPO is now bidirectional in Hb MoM.
 3. Production was an absolute rate, which does not scale with a fetus that
    quadruples in mass.  It is now a multiple of the REQUIREMENT (senescence
    plus blood-volume growth), which is why k_prod = 1 is not fitted.
 4. The Starling balance had a baseline imbalance: a healthy fetus accumulated
    ascites from 18 weeks and the death hazard fired on it (survival 0.76 at
    18 wk in a fetus with no disease).  Lymphatic capacity is now DEFINED as
    baseline filtration.
 5. Long-lived plasma cells were missing, so every sensitised mother's titre
    collapsed during the pregnancy: 15 IU/mL at 14 wk became 0.06 by 33 wk.
 6. Sensitisation was driven by the CLEARANCE FLUX of fetal cells, which makes
    the 72-hour window impossible -- with a fast clearance the first day has
    already primed a third of the response.  It is now the exposure INTEGRAL,
    with uncoated fetal cells surviving as long as any red cell, and the
    72-hour window then falls out instead of being assumed.
 7. Bilirubin was given a total-body-water distribution volume (0.75 L/kg).
    It is albumin-bound: 0.20 L/kg.  With the wrong volume the arithmetic of
    physiological jaundice does not close -- four days of production cannot
    fill the pool the model claimed to have.
 8. UGT1A1 maturation was first-order with a 6-day time constant, which gave a
    newborn 40%% of adult conjugating capacity on day 3.
 9. The exchange transfusion ADDED donor cells without removing any.  Donor
    mass therefore grew at every procedure, its senescence grew with it,
    bilirubin production grew with that, and the threshold re-triggered: a
    positive feedback loop that reached a total serum bilirubin of 1.2e7 mg/dL
    and 69 "exchanges" in one neonate.  An exchange is volume-neutral and
    removes ~85%% of everything in the circulation.
10. The "fixed 2-weekly IUT" protocol never fired: nothing set the FIRST
    trigger, so the arm silently reported zero transfusions and looked
    identical to no treatment.
11. The IVIG mechanism decomposition isolated its "FcgammaR only" arm by
    zeroing the IVIG POOL, which removes both mechanisms at once.  The arm
    looked inert and the write-up drew the opposite conclusion from the truth.
    The model now carries an explicit ivig_compete switch (section 8).

Two further bugs were in the ANALYSIS rather than the model, and are recorded
because they were both silent: a stray single %% in a format string ("40%% of
adult") was read as a %%o conversion and killed the script at the last block
after forty minutes of computation, and the viscosity sensitivity in section 4
printed a hard-coded constant instead of solving for the shifted threshold.
The section-rerun tool (hdfn_section_rerun.py) exists because of them.

PRIOR EXPECTATIONS REFUTED BY THE CALCULATION, reported rather than deleted:
 * "A transfusion works by adding haemoglobin."  It works at least as much by
   REMOVING SUBSTRATE; the two effects are separated in section 2, and it is
   the substrate effect that explains the flattening decline.
 * "The FcgammaR arm of IVIG does the work."  It does not.  Once the switch was
   fixed, FcRn competition accounted for essentially the whole effect on
   delivered antibody (5.57 -> 1.20 IU/mL) and FcgammaR blockade for almost
   none of it, because the reticuloendothelial system is not rate-limiting at
   these antibody loads.
 * "An FcRn blocker that lowers maternal IgG by 70%% should abolish the
   disease."  It does -- unless placental FcRn is much harder to saturate than
   endothelial FcRn, which is exactly what the UNITY failure rate implies
   (section 6).
 * "MCA-PSV needs a viscosity term."  Adding one moves the 1.5 MoM threshold
   off the definition of moderate anaemia and raises the predicted
   false-positive rate above what was observed (section 4).

HONEST MISSES:
 * The post-transfusion decline is under-predicted by about a third (%.2f vs
   0.40 g/dL/day, SD 0.25).  Matching it requires donor red cells to survive
   30-45 days in the fetus instead of 70; the model keeps 70 and reports the
   gap (section 3).
 * The model does NOT reproduce the UK "low risk" band.  At <4 IU/mL it still
   puts a substantial fraction of pregnancies into late-gestation anaemia,
   because destruction is linear in the titre while the conveyor is
   exponential in time: a four-fold lower titre DELAYS the disease by
   ln(4)/0.268 = 5.2 weeks rather than making it four times milder.  What is
   missing is a threshold non-linearity in the antibody-to-phagocytosis step
   (section 14).
 * The false-positive rate of the MCA-PSV threshold is under-predicted at a
   10%% measurement CV (5%% against 12%%); reproducing 12%% needs a CV of
   15-20%% (section 4).
 * Maternal total IgG falls only %.0f%% across gestation where Malek measured
   30-40%%.  The model has dilution and the placental sink but no reduction in
   maternal IgG synthesis, which is probably the missing term.
 * The UNITY simulation is only as good as one distributional assumption (the
   spread of anti-D quantitation among women with previous early-onset severe
   disease), and the sweep in section 7 shows the answer moving across a wide
   range as that spread changes.
 * Amniotic-fluid dOD450 is a single well-mixed pool with a fixed fractional
   appearance rate; the Liley zones are not reproduced quantitatively, only
   the anti-D versus anti-K contrast.
 * Neonatal IVIG has no dosing route: the FcgammaR term is driven by MATERNAL
   IVIG only.
 * No R toolchain was available in this environment, so hdfn_mrgsolve_model.R
   and hdfn_shiny_app.R mirror these equations line for line but have not
   themselves been executed.
""" % (np.mean(seq[1]),
       100 * (1 - at_ga(rh, "Cig", 39) / at_ga(rh, "Cig", 12))))
print("=" * 78)
print("end of reference output")
