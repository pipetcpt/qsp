#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
icp_calibration.py
==================
Everything in the ICP model that is *fitted* rather than *assumed* is fitted
here, in the open, and every quantitative claim made in README.md is computed
here so that it can be re-derived by re-running this file.

Sections
--------
  A  severity strata: what the four susceptibility vectors produce
  B  the 100 umol/L stillbirth threshold: placental carrier saturation
  C  fitting the three stillbirth-hazard constants to Ovadia 2019
  D  what exponent a TBA-driven hazard would need instead
  E  PITCHES: decomposing the composite primary endpoint
  F  the two-axis dissociation (bile acid axis vs itch axis)
  G  delivery-timing optimisation by bile acid stratum
  H  onset timing, postpartum resolution, twin pregnancy
  I  in-silico ablations that isolate each proposed mechanism

Run:  python3 icp_calibration.py > icp_calibration_output.txt
"""

import math
import copy
import icp_reference_model as M

P = M.P


def trace(gen=None, tdel_ga=39.0, dt=0.02, **kw):
    """Run one scenario and return (per-day rows, endpoint dict, ctl)."""
    ctl = M.make_ctl(gen=gen, tdel_ga=tdel_ga, **kw)
    rows, _ = M.simulate(ctl, dt=dt)
    return rows, M.endpoints(rows, ctl), ctl


def hazard_integral(rows, ctl, n):
    """integral of ARRI^n * (1 + KHY*HYP) dt over GA >= 24 wk, and the window
    length in days.  Trapezoidal on the 1-day record grid."""
    tdel_ga = ctl["ga0"] + ctl["tdel"] / 7.0
    sel = [r for r in rows if 24.0 <= r["ga"] <= tdel_ga + 1e-9]
    I = 0.0
    for k in range(len(sel) - 1):
        a, b = sel[k], sel[k + 1]
        fa = (a["ARRI"] ** n) * (1.0 + P["KHY"] * a["HYP"])
        fb = (b["ARRI"] ** n) * (1.0 + P["KHY"] * b["HYP"])
        I += 0.5 * (fa + fb)
    T = max(0.0, (tdel_ga - 24.0) * 7.0)
    return I, T


def hdr(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------------------
# A. severity strata
# ---------------------------------------------------------------------------
# A fourth, milder vector is needed because the Ovadia "<40 umol/L" band has a
# median around 15-20 umol/L, not 26: using the 26 umol/L vector as the
# representative of that band would push part of the fitted background hazard
# into the bile-acid-dependent term.
MILD1 = dict(gBSEP=0.795, gMDR3=0.695, gSULT=1.72)

STRATA = [("normal pregnancy", None),
          ("ICP, TBA band <40 (repr.)", MILD1),
          ("ICP, TBA band <40 (upper)", M.MILD),
          ("ICP, TBA band 40-99", M.SEVERE),
          ("ICP, TBA band >=100", M.VSEVERE)]

hdr("A. Severity strata produced by the susceptibility vectors (term, 39 wk)")
print(f"{'stratum':<28s} {'TBA':>7s} {'endo':>7s} {'cord':>7s} {'F/M':>5s} "
      f"{'FCL':>7s} {'wbar':>5s} {'ARRI':>7s} {'ALT':>5s} {'VAS':>5s} "
      f"{'onset':>6s}")
A = {}
for name, gen in STRATA:
    rows, e, ctl = trace(gen)
    A[name] = (rows, e, ctl)
    wbar = e["FCL"] / e["CORD"] if e["CORD"] > 0 else 0.0
    ons = "-" if e["onset_ga"] is None else f"{e['onset_ga']:.1f}"
    print(f"{name:<28s} {e['TBA_del']:7.1f} {e['TBA_endo_del']:7.1f} "
          f"{e['CORD']:7.1f} {e['CORD']/e['TBA_del']:5.2f} {e['FCL']:7.2f} "
          f"{wbar:5.3f} {e['ARRI']:7.4f} {e['ALT']:5.0f} {e['VAS']:5.1f} "
          f"{ons:>6s}")
print()
print("wbar = FCL/cord = mean cytotoxicity weight of the fetal pool.  It RISES")
print("with severity because the maternal->fetal diffusive permeability is")
print("higher for the hydrophobic species, so the fetal pool does not merely")
print("get bigger, it gets more hydrophobic.  That is a second, independent")
print("source of supralinearity on top of carrier saturation.")


# ---------------------------------------------------------------------------
# B. placental carrier saturation
# ---------------------------------------------------------------------------
hdr("B. Where the fetal load stops tracking the maternal load, and why")
sweep = [dict(gBSEP=b, gMDR3=m, gSULT=s) for b, m, s in [
    (1.00, 1.00, 1.00), (0.90, 0.85, 1.40), (0.845, 0.76, 1.60),
    (0.795, 0.695, 1.72), (0.735, 0.62, 1.95), (0.695, 0.56, 2.15),
    (0.650, 0.50, 2.40), (0.600, 0.44, 2.65), (0.560, 0.40, 2.85),
    (0.520, 0.35, 3.10), (0.480, 0.31, 3.35), (0.440, 0.27, 3.60)]]
print(f"{'TBA':>8s} {'cord':>8s} {'FCL':>8s} {'JM2F/VP':>8s} {'dlnFCL':>8s} "
      f"{'/dlnTBA':>8s}")
curve = []
for gen in sweep:
    rows, e, ctl = trace(dict(M.GEN0, **gen))
    # placental carrier utilisation at term
    csp = e["TBA_endo_del"]
    util = sum(P["PS"][i] for i in range(4)) / 4.0 * csp * 1.0
    curve.append((e["TBA_del"], e["CORD"], e["FCL"], util))
for k, (t, c, f, u) in enumerate(curve):
    if k == 0:
        print(f"{t:8.1f} {c:8.1f} {f:8.2f} {u/P['VP']:8.2f} {'-':>8s} {'-':>8s}")
        continue
    t0, c0, f0, _ = curve[k - 1]
    sl = (math.log(f) - math.log(f0)) / (math.log(t) - math.log(t0))
    print(f"{t:8.1f} {c:8.1f} {f:8.2f} {u/P['VP']:8.2f} "
          f"{math.log(f)-math.log(f0):8.3f} {sl:8.2f}")
best = max(range(1, len(curve)),
           key=lambda k: (math.log(curve[k][2]) - math.log(curve[k - 1][2]))
           / (math.log(curve[k][0]) - math.log(curve[k - 1][0])))
print()
print(f"steepest local log-log slope d(lnFCL)/d(lnTBA) = "
      f"{(math.log(curve[best][2])-math.log(curve[best-1][2]))/(math.log(curve[best][0])-math.log(curve[best-1][0])):.2f} "
      f"at maternal TBA {curve[best-1][0]:.0f}-{curve[best][0]:.0f} umol/L")
print("The knee is not fitted.  It is where the diffusive maternal->fetal flux")
print("plus fetal synthesis consumes the carrier capacity VP; the position")
print("follows from PS, VP and KMP alone.  What it explains is the LEVEL of")
print("fetal exposure and the fact that the cord:maternal ratio stops falling")
print("and starts rising as the carrier fills (column 2 divided by column 1:")
print("it bottoms out around 0.30 in mild disease and climbs back above 0.40")
print("in the >=100 band).  What it does NOT explain is the bend in the risk")
print("curve -- section D measures that and finds it downstream, in the")
print("myocardium.")


# ---------------------------------------------------------------------------
# C. fitting the stillbirth hazard
# ---------------------------------------------------------------------------
hdr("C. Fitting the three stillbirth-hazard constants (Ovadia 2019 strata)")
# Ovadia 2019 Lancet, IPD meta-analysis, 5269 ICP pregnancies:
#   TBA <40      stillbirth 0.13%
#   TBA 40-99    stillbirth 0.28%
#   TBA >=100    stillbirth 3.44%
TARGET = [0.0013, 0.0028, 0.0344]
reps = [MILD1, M.SEVERE, M.VSEVERE]
traces = []
for gen in reps:
    rows, e, ctl = trace(gen)
    traces.append((rows, ctl, e))

print("representative FCL at term:  " +
      "  ".join(f"{e['FCL']:.2f}" for _, _, e in traces))
print()
print(f"{'HN':>5s} {'HSB0 (/d)':>12s} {'HSBSC (/d)':>12s} "
      f"{'band1':>8s} {'band2':>8s} {'band3':>8s} {'err':>9s}")
bestfit = None
n = 1.0
while n <= 4.001:
    Is, Ts = [], []
    for rows, ctl, _ in traces:
        I, T = hazard_integral(rows, ctl, n)
        Is.append(I)
        Ts.append(T)
    # solve the 2x2 system on bands 1 and 3 exactly; band 2 is the test
    L = [-math.log(1.0 - t) for t in TARGET]
    det = Ts[0] * Is[2] - Ts[2] * Is[0]
    if abs(det) > 1e-30:
        h0 = (L[0] * Is[2] - L[2] * Is[0]) / det
        hs = (L[2] * Ts[0] - L[0] * Ts[2]) / det
        if h0 > 0 and hs > 0:
            pred = [1.0 - math.exp(-(h0 * Ts[k] + hs * Is[k])) for k in range(3)]
            err = abs(pred[1] - TARGET[1]) / TARGET[1]
            if abs(n - round(n * 2) / 2) < 1e-9:
                print(f"{n:5.2f} {h0:12.4e} {hs:12.4e} "
                      + " ".join(f"{100*p:7.3f}%" for p in pred)
                      + f" {100*err:8.2f}%")
            if bestfit is None or err < bestfit[0]:
                bestfit = (err, n, h0, hs, pred)
    n += 0.05
err, n, h0, hs, pred = bestfit
print()
print(f"BEST: HN = {n:.2f}   HSB0 = {h0:.4e} /d   HSBSC = {hs:.4e} /d")
print(f"      predicted  {100*pred[0]:.3f}% / {100*pred[1]:.3f}% / "
      f"{100*pred[2]:.3f}%")
print(f"      target     {100*TARGET[0]:.3f}% / {100*TARGET[1]:.3f}% / "
      f"{100*TARGET[2]:.3f}%   (band-2 residual {100*err:.2f}%)")
print()
print("Three constants, three data points -- so the FIT itself is not evidence.")
print("What is worth reading is the SIZE of the exponent.  HN acts on the")
print("myocardial arrhythmia index ARRI, not on bile acid, and the exponent")
print("that comes out is ~1.5.  Section D asks the same question of the two")
print("upstream variables and shows where the steepness of the clinical")
print("threshold actually lives.")
FIT = dict(HN=n, HSB0=h0, HSBSC=hs)


# ---------------------------------------------------------------------------
# D. the same fit written on maternal TBA instead
# ---------------------------------------------------------------------------
hdr("D. Where the steepness of the 100 umol/L threshold actually lives")
print("The same three outcome numbers are re-fitted three times, changing only")
print("WHICH variable the hazard is written on.  The exponent each version")
print("needs measures how much of the nonlinearity that variable already")
print("carries: a big exponent means the variable is too shallow and the")
print("dose-response is being asked to supply the steepness.")
print()
tba = [e["TBA_del"] for _, _, e in traces]
fcl = [e["FCL"] for _, _, e in traces]
arr = [e["ARRI"] for _, _, e in traces]
L = [-math.log(1.0 - t) for t in TARGET]


def integ_of(rows, ctl, key, n):
    """integral of X^n * (1+KHY*HYP) dt over GA>=24, plus the window."""
    tdel_ga = ctl["ga0"] + ctl["tdel"] / 7.0
    sel = [r for r in rows if 24.0 <= r["ga"] <= tdel_ga + 1e-9]
    I = 0.0
    for k in range(len(sel) - 1):
        a, b = sel[k], sel[k + 1]
        I += 0.5 * ((a[key] ** n) * (1.0 + P["KHY"] * a["HYP"])
                    + (b[key] ** n) * (1.0 + P["KHY"] * b["HYP"]))
    return I, max(0.0, (tdel_ga - 24.0) * 7.0)


def bestexp(key, nmax=12.0):
    best = None
    n = 0.5
    while n <= nmax + 1e-9:
        Is, Ts = [], []
        for rows, ctl, _ in traces:
            I, T = integ_of(rows, ctl, key, n)
            Is.append(I)
            Ts.append(T)
        det = Ts[0] * Is[2] - Ts[2] * Is[0]
        if abs(det) > 1e-60:
            a = (L[0] * Is[2] - L[2] * Is[0]) / det
            b = (L[2] * Ts[0] - L[0] * Ts[2]) / det
            if a > 0 and b > 0:
                p2 = 1.0 - math.exp(-(a * Ts[1] + b * Is[1]))
                e2 = abs(p2 - TARGET[1]) / TARGET[1]
                if best is None or e2 < best[0]:
                    best = (e2, n, a, b, p2)
        n += 0.05
    return best


print(f"{'hazard written on':<40s} {'band1':>8s} {'band2':>8s} {'band3':>8s} "
      f"{'b3/b1':>7s} {'needs n':>8s} {'resid':>7s}")
summary = {}
for label, key, vals in (
        ("maternal total bile acid (the assay)", "TBA", tba),
        ("fetal cytotoxic load FCL", "FCL", fcl),
        ("myocardial arrhythmia index ARRI", "ARRI", arr)):
    b = bestexp(key)
    summary[key] = b
    print(f"{label:<40s} {vals[0]:8.3f} {vals[1]:8.3f} {vals[2]:8.3f} "
          f"{vals[2]/vals[0]:7.1f} {b[1]:8.2f} {100*b[0]:6.2f}%")
print()
print("Reading the table:")
print(f"  * maternal TBA spans {tba[2]/tba[0]:.1f}x across the three bands and")
print(f"    needs an exponent of {summary['TBA'][1]:.2f}.")
print(f"  * the fetal cytotoxic load spans {fcl[2]/fcl[0]:.1f}x -- barely more")
print("    than the maternal assay -- and so needs an exponent of")
print(f"    {summary['FCL'][1]:.2f}, essentially unchanged.  Transport alone")
print("    does NOT explain the threshold.  (What transport does explain is")
print("    the LEVEL and the rising cord:maternal ratio; see section B.)")
print(f"  * the arrhythmia index spans {arr[2]/arr[0]:.1f}x on the same three")
print(f"    pregnancies and needs only {summary['ARRI'][1]:.2f}.")
print()
print("So the amplification is in the MYOCARDIAL RESPONSE, not in the")
print("placenta and not in the dose-response.  Connexin-43 uncoupling is a")
print("threshold process (Hill 1.6 around IC50 ~30 in weighted units) and")
print("calcium overload multiplies it, so a 9-fold change in fetal exposure")
print(f"becomes a {arr[2]/arr[0]:.0f}-fold change in electrical instability.")
print("That matters clinically in a specific way: it predicts the risk is")
print("carried by a step change in myocyte coupling rather than by a graded")
print("increase in 'fetal distress', which is why the events are abrupt and")
print("why cardiotocography does not see them coming.")
print()
print("This is also the model correcting its own author.  The map was drawn")
print("expecting placental carrier saturation to be the source of the")
print("threshold; the fit says it is not, and the ablations in section I")
print("agree.  The saturation is real and it does set how much bile acid")
print("reaches the fetus -- it just is not what makes the risk curve bend.")


# ---------------------------------------------------------------------------
# E. PITCHES
# ---------------------------------------------------------------------------
hdr("E. PITCHES (Chappell 2019): decomposing the composite primary endpoint")
# Trial: 605 women with ICP, UDCA 500 mg bid (titrated) vs placebo.
# Primary composite = perinatal death OR preterm delivery <37 wk
#                     OR neonatal unit admission >= 4 h
#   UDCA     74/298 = 24.8%
#   placebo  85/306 = 27.8%   adjusted RR 0.85 (0.62-1.15), not significant
# Secondary: worst itch VAS improved by 0.7 cm on a 10 cm scale.
# The trial population was mostly mild-moderate ICP, so the model is run on
# the MILD and SEVERE vectors and weighted 70/30, which is roughly the
# severity mix reported in the trial.
def arm(gen, udca):
    rows, e, ctl = trace(gen, udca_mg=udca, udca_start=30.0)
    return e


mix = [(M.MILD, 0.70), (M.SEVERE, 0.30)]
out = {}
for label, dose in (("placebo", 0.0), ("UDCA 1000 mg/d", 1000.0)):
    agg = dict(SB=0.0, PTB=0.0, NICU=0.0, COMP=0.0, TBA=0.0, VAS=0.0, ALT=0.0,
               FCL=0.0)
    for gen, w in mix:
        e = arm(gen, dose)
        agg["SB"] += w * e["SB"]
        agg["PTB"] += w * e["PTB_spont"]
        agg["NICU"] += w * e["NICU"]
        agg["COMP"] += w * e["COMPOSITE"]
        agg["TBA"] += w * e["TBA_del"]
        agg["VAS"] += w * e["VAS"]
        agg["ALT"] += w * e["ALT"]
        agg["FCL"] += w * e["FCL"]
    out[label] = agg
p_, u_ = out["placebo"], out["UDCA 1000 mg/d"]
print(f"{'endpoint':<34s} {'placebo':>10s} {'UDCA':>10s} {'rel.':>9s} "
      f"{'trial':>16s}")
rows_ = [
    ("total bile acid, umol/L", p_["TBA"], u_["TBA"], "falls"),
    ("fetal cytotoxic load (model)", p_["FCL"], u_["FCL"], "not measurable"),
    ("ALT, U/L", p_["ALT"], u_["ALT"], "falls"),
    ("worst itch, 0-10 cm", p_["VAS"], u_["VAS"], "-0.7 cm"),
    ("stillbirth, %", 100 * p_["SB"], 100 * u_["SB"], "1 vs 2 deaths"),
    ("spontaneous PTB <37 wk, %", 100 * p_["PTB"], 100 * u_["PTB"], "-"),
    ("neonatal unit >=4 h, %", 100 * p_["NICU"], 100 * u_["NICU"], "-"),
    ("COMPOSITE, %", 100 * p_["COMP"], 100 * u_["COMP"], "27.8 vs 24.8"),
]
for nm, a, b, tr in rows_:
    rel = f"{100*(b/a-1):+.1f}%" if a > 0 else "-"
    print(f"{nm:<34s} {a:10.2f} {b:10.2f} {rel:>9s} {tr:>16s}")
print()
sb_share = p_["SB"] / p_["COMP"]
print(f"Stillbirth is {100*sb_share:.2f}% of the composite in the placebo arm.")
print("So even abolishing it entirely moves the composite by "
      f"{100*p_['SB']:.2f} percentage points.")
print(f"The model's own stillbirth reduction is "
      f"{100*(1-u_['SB']/p_['SB']):.0f}%, worth "
      f"{100*(p_['SB']-u_['SB']):.3f} points on a composite running at "
      f"{100*p_['COMP']:.0f}%.")
print()
print("That is the whole of PITCHES.  Two of the three components of the")
print("primary endpoint -- preterm delivery and neonatal unit admission --")
print("are set by gestational age at delivery and by the decision to deliver,")
print("neither of which UDCA touches; the model reproduces this because both")
print("are functions of GA at delivery, not of bile acid.  The one component")
print("UDCA can move contributes a fraction of a percentage point.  With 604")
print("women the trial could not have detected its own mechanism, and the")
print("result is therefore consistent with the Ovadia 2021 individual-")
print("participant meta-analysis reporting a stillbirth signal: the two")
print("studies are measuring different things, not disagreeing.")
print()
print("Detectability check -- women per arm needed for 80% power at alpha 0.05")
print("on each component, assuming the model's effect sizes are correct:")


def nper(p1, p2):
    if p1 <= 0 or p2 <= 0 or abs(p1 - p2) < 1e-12:
        return float("inf")
    pbar = 0.5 * (p1 + p2)
    return (1.96 * math.sqrt(2 * pbar * (1 - pbar))
            + 0.84 * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2 \
        / (p1 - p2) ** 2


for nm, a, b in (("composite", p_["COMP"], u_["COMP"]),
                 ("stillbirth alone", p_["SB"], u_["SB"]),
                 ("neonatal unit alone", p_["NICU"], u_["NICU"])):
    v = nper(a, b)
    print(f"   {nm:<22s} n/arm = "
          + ("> 1e7" if v > 1e7 else f"{v:,.0f}"))


# ---------------------------------------------------------------------------
# F. two-axis dissociation
# ---------------------------------------------------------------------------
hdr("F. The bile acid axis and the itch axis respond to different drugs")
DRUGS = [
    ("none", {}),
    ("UDCA 1000 mg/d", dict(udca_mg=1000)),
    ("UDCA 1500 mg/d", dict(udca_mg=1500)),
    ("rifampicin 600 mg/d", dict(rif_mg=600)),
    ("UDCA + rifampicin", dict(udca_mg=1000, rif_mg=600)),
    ("cholestyramine 16 g/d", dict(chol_g=16.0)),
    ("SAMe 1000 mg/d IV", dict(sam_mg=1000)),
    ("IBAT inhibitor (hypoth.)", dict(odev_umol=4.0)),
    ("naltrexone 50 mg/d", dict(ntx_mg=50.0)),
    ("antihistamine 12 mg/d", dict(ah_mg=12.0)),
]
print(f"{'regimen (on the 40-99 stratum)':<30s} {'TBA':>7s} {'dTBA':>7s} "
      f"{'FCL':>7s} {'dFCL':>7s} {'ATX':>6s} {'VAS':>5s} {'dVAS':>6s} "
      f"{'ALT':>5s} {'INR':>5s} {'SB%':>7s}")
base = None
for nm, kw in DRUGS:
    rows, e, ctl = trace(M.SEVERE, **kw)
    if base is None:
        base = e
    print(f"{nm:<30s} {e['TBA_del']:7.1f} "
          f"{100*(e['TBA_del']/base['TBA_del']-1):+6.0f}% {e['FCL']:7.2f} "
          f"{100*(e['FCL']/base['FCL']-1):+6.0f}% {e['ATX']:6.2f} "
          f"{e['VAS']:5.1f} {e['VAS']-base['VAS']:+6.1f} {e['ALT']:5.0f} "
          f"{e['INR']:5.2f} {100*e['SB']:7.3f}")
print()
print("Read the dFCL and dVAS columns together.  Every drug that lowers the")
print("bile acid pool (UDCA, cholestyramine, an IBAT inhibitor) leaves itch")
print("essentially where it was, and the one drug that abolishes itch")
print("(rifampicin) is mid-table on bile acid.  The model does not contain a")
print("rule that says so; it falls out of driving autotaxin from the sex-")
print("steroid load with only a weak, saturating bile-acid term.  The")
print("practical consequence is that a trial powered on itch cannot rank")
print("bile-acid-lowering drugs, and a trial powered on bile acid cannot")
print("rank antipruritics -- which is roughly the state of the ICP")
print("literature.")
print()
print("Note the INR column: cholestyramine is the only agent that lowers")
print("bile acid AND worsens vitamin K status, because it lowers the")
print("intestinal bile acid concentration that fat-soluble vitamin")
print("absorption depends on.  Rifampicin does the same for a different")
print("reason (it shrinks the pool).  Both trade postpartum haemorrhage")
print("risk for bile acid, and the model prices the trade.")


# ---------------------------------------------------------------------------
# G. delivery timing
# ---------------------------------------------------------------------------
hdr("G. Optimal gestational age at delivery, by bile acid stratum")
# Expected-loss framing.  Stillbirth is weighted 1 (loss of the pregnancy);
# the iatrogenic cost of early delivery is the excess respiratory morbidity
# and neonatal unit admission at that GA, weighted by WMORB.  WMORB is the
# only judgement call here and it is reported, not hidden: it is the
# disutility of a neonatal unit admission relative to a stillbirth.
WMORB = 0.055
print(f"(neonatal morbidity weighted at {WMORB:.3f} of a stillbirth)")
print()
print(f"{'stratum':<26s} " + " ".join(f"{g:>7.0f}" for g in
                                      (35, 36, 37, 38, 39, 40)) + "   opt")
for name, gen in (("TBA <40 (repr.)", MILD1), ("TBA ~26", M.MILD),
                  ("TBA 40-99", M.SEVERE), ("TBA >=100", M.VSEVERE),
                  ("TBA >=100 + UDCA", M.VSEVERE)):
    losses, gas = [], []
    for ga in (35, 36, 37, 38, 39, 40):
        kw = dict(udca_mg=1000) if "UDCA" in name else {}
        rows, e, ctl = trace(gen, tdel_ga=float(ga), **kw)
        # re-score stillbirth with the fitted constants
        I, T = hazard_integral(rows, ctl, FIT["HN"])
        sb = 1.0 - math.exp(-(FIT["HSB0"] * T + FIT["HSBSC"] * I))
        loss = sb + WMORB * (e["NICU"] + e["RDS"])
        losses.append(loss)
        gas.append(ga)
    opt = gas[min(range(len(losses)), key=lambda k: losses[k])]
    print(f"{name:<26s} " + " ".join(f"{1000*l:7.2f}" for l in losses)
          + f"   {opt:>3d} wk")
print()
print("Units are milli-expected-losses; only the position of the minimum")
print("matters.")
print()
print("The >=100 band comes out at 37 weeks, which is the direction the")
print("guidelines point (SMFM suggests 36-37 wk at or above 100 umol/L), and")
print("treating that pregnancy pushes the optimum back to 39 weeks -- i.e. in")
print("this model the value of UDCA is not a biochemical endpoint, it is about")
print("two weeks of gestation.  That is a testable and, as far as we can tell,")
print("untested claim.")
print()
print("The 40-99 band does NOT come out where guidelines put it, and the model")
print("should not be read as supporting them there.  It puts the optimum at")
print("39-40 weeks, not 37-38.  The reason is arithmetic and worth stating")
print("plainly: the stillbirth excess Ovadia reports for that band is 0.28%")
print("against 0.13% -- fifteen hundredths of a percentage point -- and no")
print("plausible weighting makes fifteen hundredths of a point outweigh the")
print("neonatal cost of delivering two to three weeks early.  Below is the")
print("neonatal-morbidity weight that WOULD be required:")
print()
print(f"{'WMORB':>8s} {'<40':>6s} {'40-99':>6s} {'>=100':>6s}")
for wm in (0.005, 0.010, 0.020, 0.035, 0.055, 0.080, 0.120):
    line = []
    for gen in (MILD1, M.SEVERE, M.VSEVERE):
        best_ga, best_l = None, None
        for ga in (35, 36, 37, 38, 39, 40):
            rows, e, ctl = trace(gen, tdel_ga=float(ga))
            I, T = hazard_integral(rows, ctl, FIT["HN"])
            sb = 1.0 - math.exp(-(FIT["HSB0"] * T + FIT["HSBSC"] * I))
            l = sb + wm * (e["NICU"] + e["RDS"])
            if best_l is None or l < best_l:
                best_l, best_ga = l, ga
        line.append(best_ga)
    print(f"{wm:8.3f} " + " ".join(f"{g:6d}" for g in line))
print()
print("Only at a weight below about 0.01 -- i.e. only if a neonatal unit")
print("admission is worth less than a hundredth of a stillbirth -- does the")
print("40-99 band move to 37 weeks, and at that weight the <40 band moves")
print("with it, which no guideline recommends.  So either the model is missing")
print("something real in that band, or the 37-38 week recommendation is not")
print("carried by the stillbirth numbers it is usually justified by.  Two")
print("candidates for what the model is missing, both plausible: within-woman")
print("bile acid fluctuates and a woman measured at 70 umol/L may spend days")
print("above 100 between visits (the model runs a smooth trajectory), and the")
print("hazard here is an average over a band whose upper edge is much more")
print("dangerous than its middle.  A version of this model driven by serial")
print("measurements rather than a smooth curve would be the way to test that,")
print("and it is the single most useful extension we can identify.")

# ---------------------------------------------------------------------------
# H. onset, resolution, twins
# ---------------------------------------------------------------------------
hdr("H. Onset timing, postpartum resolution, and twin pregnancy")
print("None of these three is fitted.  Onset and resolution are consequences")
print("of the placental sex-steroid trajectory (E2 and P4-sulfate are ODEs")
print("with placental production switched off at delivery); the twin effect is")
print("a single covariate on that trajectory.")
print()
print(f"{'scenario':<34s} {'onsetGA':>8s} {'TBApk':>7s} {'pkGA':>6s} "
      f"{'<10 by':>8s} {'VASpk':>6s}")
for name, gen, kw in (
        ("normal", None, {}),
        ("ICP <40 stratum", MILD1, {}),
        ("ICP 40-99 stratum", M.SEVERE, {}),
        ("ICP >=100 stratum", M.VSEVERE, {}),
        ("twin, <40 genotype", dict(MILD1, twin=1.55), {}),
        ("twin, 40-99 genotype", dict(M.SEVERE, twin=1.55), {}),
        ("ICP >=100 + UDCA 1000", M.VSEVERE, dict(udca_mg=1000))):
    rows, e, ctl = trace(dict(M.GEN0, **gen) if gen else None, **kw)
    ons = "never" if e["onset_ga"] is None else f"{e['onset_ga']:.1f}"
    res = "n/a" if e["resolve_d"] is None else f"{e['resolve_d']:.0f} d"
    print(f"{name:<34s} {ons:>8s} {e['TBA_peak']:7.1f} "
          f"{e['TBA_peak_ga']:6.1f} {res:>8s} {e['VAS_peak']:6.1f}")
print()
print("The same susceptibility vector that is silent at 20 weeks crosses the")
print("10 umol/L diagnostic line in the third trimester and clears within")
print("about two weeks of delivery, with no parameter changed at either")
print("boundary -- the sex-steroid trajectory does both.  A twin pregnancy")
print("shifts onset earlier and raises the peak on the same genotype, which")
print("is the mechanistic reading of the ~5-fold higher ICP incidence in")
print("twins.")


# ---------------------------------------------------------------------------
# I. ablations
# ---------------------------------------------------------------------------
hdr("I. In-silico ablations: which mechanism carries which result")
saved = copy.deepcopy(P)
W0 = list(M.WTOX)


def ablation_metrics(label, changes=None, weights=None):
    """Re-run the three representative pregnancies under an ablation and
    report the quantities the threshold actually depends on."""
    changes = changes or {}
    for k, v in changes.items():
        P[k] = v
    if weights is not None:
        M.WTOX[:] = weights
    tr, vals = [], []
    for gen in (MILD1, M.SEVERE, M.VSEVERE):
        rows, e, ctl = trace(gen)
        tr.append((rows, ctl, e))
        vals.append(e)
    # required exponent if the hazard is written on ARRI
    globals()["traces"] = tr
    b = bestexp("ARRI")
    bt = bestexp("TBA")
    for k in changes:
        P[k] = saved[k]
    M.WTOX[:] = W0
    return (label, [v["TBA_del"] for v in vals], [v["FCL"] for v in vals],
            [v["ARRI"] for v in vals], b[1], bt[1])


ABL = [
    ablation_metrics("(reference)"),
    ablation_metrics("placental carrier not saturable (V_P x 20)",
                     dict(VP=20 * saved["VP"])),
    ablation_metrics("no maternal-side occupancy (TRANSIN = 0)",
                     dict(TRANSIN=0.0)),
    ablation_metrics("species-uniform diffusive PS (all = 0.77)",
                     dict(PS=[0.77, 0.77, 0.77, 0.77, 0.15])),
    ablation_metrics("uniform cytotoxicity weights (all = 0.5)",
                     weights=[0.5, 0.5, 0.5, 0.5, 0.5]),
    ablation_metrics("no Cx43 uncoupling (IMAXGJ = 0.02)",
                     dict(IMAXGJ=0.02)),
    ablation_metrics("Cx43 uncoupling made linear (HGJ = 1, IC50 x 6)",
                     dict(HGJ=1.0, IC50GJ=6 * saved["IC50GJ"])),
    ablation_metrics("no Ca2+ overload term (ECA = 0)", dict(ECA=0.0)),
]
# restore the reference traces for anything downstream
traces = []
for gen in reps:
    rows, e, ctl = trace(gen)
    traces.append((rows, ctl, e))

print(f"{'ablation':<44s} {'ARRI b1':>8s} {'b3':>8s} {'b3/b1':>7s} "
      f"{'n(ARRI)':>8s} {'n(TBA)':>7s}")
for label, tb, fc, ar, nA, nT in ABL:
    print(f"{label:<44s} {ar[0]:8.4f} {ar[2]:8.4f} {ar[2]/ar[0]:7.1f} "
          f"{nA:8.2f} {nT:7.2f}")
print()
print("Columns to read: b3/b1 is how much electrical instability separates the")
print(">=100 band from the <40 band, and n(ARRI) is the exponent the hazard")
print("then needs.  A mechanism matters here if removing it COLLAPSES b3/b1")
print("and forces n(ARRI) up towards the n(TBA) column.")
print()
print("  * making the placental carrier non-saturable barely moves b3/b1.")
print("    Carrier saturation sets how much bile acid reaches the fetus, not")
print("    how nonlinear the response is.")
print("  * flattening the species-specific diffusive permeabilities, or the")
print("    cytotoxicity weights, likewise leaves the ratio in the same place.")
print("  * what destroys it is LINEARISING connexin-43 uncoupling: b3/b1 falls")
print("    from 37 to 13 and the exponent the hazard needs climbs from 1.55")
print("    towards the bile-acid figure.  Removing the Ca2+ overload term")
print("    does the same thing less strongly.")
print("  * note that merely SCALING the uncoupling down (IMAXGJ = 0.02) leaves")
print("    b3/b1 untouched, because a Hill function scaled by a constant has")
print("    the same shape.  It is the threshold SHAPE that carries the")
print("    result, not the amplitude -- which is a sharper claim, and a")
print("    harder one to satisfy by accident.")
print()
print("The conclusion is specific and falsifiable: the 100 umol/L threshold is")
print("a property of gap-junction uncoupling in fetal myocardium, and it")
print("should therefore track the hydrophobic fraction of the fetal pool")
print("rather than maternal total bile acid.  Cord-blood speciation at")
print("delivery would test it directly; the maternal assay never can.")
print()

# UDCA-specific ablations
hdr("I2. UDCA: which of its mechanisms does what")
print(f"{'ablation':<48s} " + "   ".join(
    f"{d:>6d} mg/d" for d in (500, 1000, 2000)))
print(f"{'  (cells are maternal TBA / fetal cytotoxic load)':<48s}")
for name, changes in (
        ("reference", {}),
        ("no BSEP trafficking effect (ETRAF=0)", dict(ETRAF=0.0)),
        ("no gut 7beta-dehydroxylation to LCA (KDH_UD=0)", dict(KDH_UD=0.0)),
        ("no hepatic LCA sulfation (VLCAS=0)", dict(VLCAS=0.0)),
        ("no placental carrier competition (KMP_UDCA huge)", dict()),
):
    for k, v in changes.items():
        P[k] = v
    if "placental" in name:
        kmp = list(P["KMP"])
        P["KMP"] = kmp[:4] + [1e9]
    rows0, e0, _ = trace(M.SEVERE)
    cells = []
    for dose in (500, 1000, 2000):
        rows, e, ctl = trace(M.SEVERE, udca_mg=dose)
        cells.append(f"{e['TBA_del']:6.1f}/{e['FCL']:5.2f}")
    print(f"{name:<48s} " + "   ".join(cells))
    for k in changes:
        P[k] = saved[k]
    P["KMP"] = list(saved["KMP"])
print()
print("The BSEP trafficking effect is the therapeutic mechanism: remove it and")
print("most of the drug's effect goes with it.  The gut conversion of UDCA to")
print("lithocholate is a real but modest liability at these doses (removing it")
print("improves the drug by a few percentage points, and the gap widens with")
print("dose because the LCA-generating route is first-order in dose while the")
print("hepatic sulfation route that disposes of LCA saturates).  Removing")
print("hepatic LCA sulfation entirely is what makes UDCA look bad -- which is")
print("the model's reading of why UDCA is tolerated in ordinary ICP and is a")
print("genuine concern in the severest cholestasis, where that same enzyme")
print("system is already working against a much larger substrate load.")

print()
print("=" * 78)
print("Fitted constants in the whole model, exhaustively:")
print(f"  HSB0  = {FIT['HSB0']:.4e} /d   background stillbirth hazard")
print(f"  HSBSC = {FIT['HSBSC']:.4e} /d   scale on the arrhythmia index")
print(f"  HN    = {FIT['HN']:.2f}            exponent on the arrhythmia index")
print("  HPT0, HM0, HMSC   background spontaneous-preterm and meconium rates")
print("  VSCALE            cm per unit of the central itch state")
print("  WMORB             neonatal morbidity : stillbirth disutility ratio")
print("Everything else is a transport, binding or turnover constant taken")
print("from the literature ranges cited in icp_references.md, or a scale")
print("chosen so that an uncomplicated pregnancy sits at its normal value.")
print("=" * 78)
