#!/usr/bin/env python3
"""Calibration and derived-quantity analysis for the RILI model.

Two scaling constants in the model are not identifiable from mechanism and
are fitted here against published dose-response data:

  PNI50   the pneumonitis index giving 50% risk of CTCAE grade >=2 RP
  PNISL   the logistic slope on ln(index)

Everything else is either a rate constant taken from the literature or
algebraically pinned (see the [DEFECT n] notes in the model file).

Calibration targets (QUANTEC lung, Marks 2010; Seppenwoolde 2003):
  mean lung dose 13 Gy -> 10% grade >=2      MLD 20 Gy -> 20%
  mean lung dose 24 Gy -> 30%                TD50 (gEUD) ~ 30.8 Gy
  peripheral SBRT                            ~ 9-10% grade >=2

The script then MEASURES three things the model was not told:
  1. the Lyman-Kutcher-Burman volume-effect exponent n, by finding which
     generalised EUD exponent a=1/n collapses the model's own risk onto a
     single curve across 300 random dose-volume histograms
  2. the latency of peak pneumonitis and its dependence on dose
  3. the total dose at which a bin converts to the fibrotic attractor
"""

import math
import random
import sys

sys.path.insert(0, "/home/user/qsp/radiation-induced-lung-injury")
import rili_reference_model as M   # noqa: E402


# ---------------------------------------------------------------------
def peak_index(v, DRX, NFX, tend=400.0, dt=0.10, **kw):
    sim = M.Sim((DRX, NFX, v), M.Regimen(), **kw)
    out, y = M.integrate(sim, tend, dt=dt)
    pk = max(out, key=lambda r: r["PNIe"])
    return pk["PNIe"], pk["t"] - sim.TCOURSE, out, y, sim


# ---------------------------------------------------------------------
# 1.  build a family of plans spanning the clinical MLD range
# ---------------------------------------------------------------------
def scaled_plan(base, k):
    """Scale a DVH toward/away from the low-dose bin to sweep MLD while
    keeping a realistic histogram shape."""
    v = list(base)
    lo = v[0]
    for b in range(1, M.NB):
        v[b] *= k
    v[0] = max(1.0 - sum(v[1:]), 0.0)
    if v[0] < 0.0:
        return None
    s = sum(v)
    return [x / s for x in v], lo


BASE = [0.40, 0.20, 0.12, 0.11, 0.09, 0.08]

print("=" * 96)
print("1.  DOSE-RESPONSE FAMILY  (60 Gy / 30 fx, histogram scaled)")
print("=" * 96)
print("%8s %8s %8s %9s %9s" % ("scale", "MLD", "V20", "PNIpk", "tPk(d)"))
fam = []
for k in (0.15, 0.30, 0.45, 0.60, 0.80, 1.00, 1.20, 1.45, 1.70):
    r = scaled_plan(BASE, k)
    if r is None:
        continue
    v, _ = r
    mld = M.mean_lung_dose(v, 60.0)
    pk, tpk, _, _, _ = peak_index(v, 60.0, 30)
    fam.append((mld, pk, tpk, v))
    print("%8.2f %8.2f %8.1f %9.4f %9.1f"
          % (k, mld, 100 * M.vx(v, 60.0, 20.0), pk, tpk))


# ---------------------------------------------------------------------
# 2.  fit PNI50 / PNISL to QUANTEC
# ---------------------------------------------------------------------
TARGETS = [(13.0, 0.10), (20.0, 0.20), (24.0, 0.30)]


def interp_pk(mld):
    xs = [f[0] for f in fam]
    ys = [f[1] for f in fam]
    if mld <= xs[0]:
        return ys[0]
    if mld >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if mld <= xs[i]:
            w = (mld - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + w * (ys[i] - ys[i - 1])
    return ys[-1]


best = None
for p50 in [0.30 + 0.005 * i for i in range(200)]:
    for sl in [0.20 + 0.01 * i for i in range(60)]:
        sse = 0.0
        for mld, tgt in TARGETS:
            pk = interp_pk(mld)
            z = (math.log(pk) - math.log(p50)) / sl
            pr = 1.0 / (1.0 + math.exp(-z))
            sse += (pr - tgt) ** 2
        if best is None or sse < best[0]:
            best = (sse, p50, sl)
print("\nfitted:  PNI50 = %.4f   PNISL = %.3f   (SSE %.2e)"
      % (best[1], best[2], best[0]))
P50, PSL = best[1], best[2]


def ntcp(pk):
    z = (math.log(pk) - math.log(P50)) / PSL
    return 1.0 / (1.0 + math.exp(-z))


print("\n%8s %10s %10s %10s" % ("MLD", "model", "target", "source"))
for mld, tgt in TARGETS:
    print("%8.1f %9.1f%% %9.1f%%   QUANTEC" % (mld, 100 * ntcp(interp_pk(mld)),
                                               100 * tgt))
for nm, tgt, src, dl in [
        ("sbrt54", 0.095, "SBRT, nominal reserve", 85.0),
        ("sbrt54", 0.095, "SBRT, inoperable COPD reserve", 55.0)]:
    DRX, NFX, v = M.make_plan(nm)
    pk, tpk, _, _, _ = peak_index(v, DRX, NFX, DLCO0=dl)
    print("%8s %9.1f%% %9.1f%%   %s"
          % (nm, 100 * ntcp(pk), 100 * tgt, src))


# ---------------------------------------------------------------------
# 3.  MEASURE the volume-effect exponent the model produces
# ---------------------------------------------------------------------
print("\n" + "=" * 96)
print("2.  DERIVED LYMAN-KUTCHER-BURMAN VOLUME-EFFECT EXPONENT")
print("=" * 96)
print("300 random DVHs at 60 Gy / 30 fx.  For each candidate a = 1/n we")
print("regress ln(peak index) on ln(gEUD_a) and record the residual SD.")
print("The a that minimises it is the exponent the biology implies.")

random.seed(20260804)
pts = []
for _ in range(300):
    w = [random.random() ** 2 for _ in range(M.NB)]
    w[0] += 0.5
    s = sum(w)
    v = [x / s for x in w]
    if M.mean_lung_dose(v, 60.0) < 3.0:
        continue
    pk, _, _, _, _ = peak_index(v, 60.0, 30, tend=300.0, dt=0.15)
    pts.append((v, pk))
print("   usable DVHs: %d" % len(pts))

rows = []
for ai in range(1, 61):
    a = 0.4 + 0.1 * ai
    xs = [math.log(M.gEUD(v, 60.0, a)) for v, _ in pts]
    ys = [math.log(pk) for _, pk in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    b1 = sxy / sxx
    b0 = my - b1 * mx
    resid = [ys[i] - (b0 + b1 * xs[i]) for i in range(n)]
    sd = math.sqrt(sum(r * r for r in resid) / (n - 2))
    rows.append((sd, a, b1))
rows.sort()
sd, abest, slope = rows[0]
print("\n   best a = %.2f  ->  n = 1/a = %.3f   (residual SD %.4f, slope %.3f)"
      % (abest, 1.0 / abest, sd, slope))
print("   published lung n (Seppenwoolde 2003, NSCLC): 0.99")
print("\n   residual SD across candidate exponents:")
for sd, a, _ in sorted(rows, key=lambda r: r[1])[::6]:
    bar = "#" * int(sd * 220)
    print("      a=%5.2f  n=%5.2f  SD=%.4f %s" % (a, 1.0 / a, sd, bar))


# ---------------------------------------------------------------------
# 4.  latency vs dose
# ---------------------------------------------------------------------
print("\n" + "=" * 96)
print("3.  LATENCY IS SET BY AT2 TURNOVER, NOT BY DOSE")
print("=" * 96)
print("%10s %8s %10s %10s %10s" % ("plan", "MLD", "PNIpk", "tPk(d)", "tPk(wk)"))
for k in (0.25, 0.50, 1.00, 1.50):
    v, _ = scaled_plan(BASE, k)
    pk, tpk, _, _, _ = peak_index(v, 60.0, 30)
    print("%10.2f %8.2f %10.4f %10.1f %10.1f"
          % (k, M.mean_lung_dose(v, 60.0), pk, tpk, tpk / 7.0))


print("\n   sensitivity of the peak time to KMIT (1/AT2 turnover):")
for km in (0.0167, 0.0250, 0.0400, 0.0556):
    v, _ = scaled_plan(BASE, 1.0)
    pk, tpk, _, _, _ = peak_index(v, 60.0, 30, params=dict(KMIT=km))
    print("      KMIT=%.4f (tau %4.1f d)  tPk = %5.1f d = %4.1f wk  PNIpk %.4f"
          % (km, 1.0 / km, tpk, tpk / 7.0, pk))


# ---------------------------------------------------------------------
# 5.  the fibrotic conversion threshold in Gy
# ---------------------------------------------------------------------
print("\n" + "=" * 96)
print("4.  DOSE AT WHICH A BIN CONVERTS TO THE FIBROTIC ATTRACTOR")
print("=" * 96)
print("A uniform whole-lung irradiation at 2 Gy/fx, dose swept.  COL(730 d)")
print("above the separatrix at 1.586 means the bin has latched.")
print("%8s %8s %10s %10s %10s %8s" % ("D (Gy)", "NFX", "BED3", "COL730",
                                      "TGFB730", "latched"))
for D in (20, 30, 36, 40, 44, 48, 54, 60, 70):
    NFX = int(round(D / 2.0))
    v = [0.0] * M.NB
    v[5] = 1.0                    # put all volume in the top bin
    DRX = D / M.REP_F[5]
    sim = M.Sim((DRX, NFX, v), M.Regimen())
    out, y = M.integrate(sim, 730.0, dt=0.10)
    col = y[M.sidx(5, "COL")]
    print("%8.1f %8d %10.1f %10.4f %10.4f %8s"
          % (D, NFX, sim.BEDLtot[5], col, y[M.sidx(5, "TGFB")],
             "YES" if col > 1.586 else "-"))
