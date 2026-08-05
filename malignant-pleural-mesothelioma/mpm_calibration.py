#!/usr/bin/env python3
"""
mpm_calibration.py
------------------
Calibrate the MPM QSP model against published trial endpoints, then run every
structural check.  Nothing here is fitted by eye: each free constant is moved
by a one-dimensional bisection against ONE endpoint, in an order chosen so that
each step does not disturb the endpoint fitted before it.

ORDER OF IDENTIFICATION (each line adds exactly one target and one constant)

  1. baseline geometry          -> fixed from CT/manometry literature, not fitted
  2. KG0        <- untreated growth, mRECIST roughly +55% at 6 months
  3. HZ0        <- best supportive care median OS 7.0 months
  4. EMAX_CIS   <- cisplatin monotherapy median OS 9.3 mo  (EMPHACIS control)
  5. EMAX_PEM   <- pemetrexed + cisplatin median OS 12.1 mo (EMPHACIS)
  6. BEV_LAM    <- bevacizumab increment +2.7 mo           (MAPS)
  7. EMAX_IO    <- nivolumab + ipilimumab / chemo median ratio 1.28 (CheckMate 743)

  HELD OUT (never used to fit anything):
    - the CheckMate 743 HISTOLOGY SUBGROUPS (epithelioid ratio 1.13,
      non-epithelioid ratio 2.06).  Step 7 fits ONE number at the trial's
      population histology; the split between the two subgroups is then a
      prediction of the EMT axis.
    - second-line nivolumab monotherapy (CONFIRM, 10.2 vs 6.6 mo)
    - the direction of the MARS2 surgical result
    - ANC nadir day and the folate effect
"""
import math
import sys
import mpm_reference_model as M


TMAX = 1150.0
DT = 0.04


def os_of(name, x, tmax=TMAX, **kw):
    p, reg = M.scenario(name, x, **kw)
    rec = M.simulate(p, tmax=tmax, dt=DT, regimen=reg)
    return M.median_os(rec), rec


def bisect(setter, target, lo, hi, tol=0.05, itmax=18, label=""):
    """Move ONE parameter until ONE endpoint hits its target."""
    flo, fhi = setter(lo), setter(hi)
    if (flo - target) * (fhi - target) > 0:
        print("   ! %s: target %.3f not bracketed by [%.5g -> %.3f, %.5g -> %.3f]"
              % (label, target, lo, flo, hi, fhi))
        return (lo if abs(flo - target) < abs(fhi - target) else hi)
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fm = setter(mid)
        if abs(fm - target) < tol:
            print("   %-12s = %-10.5g  -> %.3f (target %.3f)" % (label, mid, fm, target))
            return mid
        if (flo - target) * (fm - target) <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    print("   %-12s = %-10.5g  -> %.3f (target %.3f)" % (label, mid, fm, target))
    return mid


def calibrate():
    fitted = {}
    print("=" * 78)
    print(" STEPWISE CALIBRATION")
    print("=" * 78)

    # --- 3. baseline hazard from best supportive care ----------------------
    def f_hz(v):
        p, reg = M.scenario("bsc", 0.25)
        p["HZ0"] = v
        return M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)) or 99.0
    fitted["HZ0"] = bisect(f_hz, 7.0, 0.00030, 0.00080, label="HZ0")

    # --- 4. cisplatin monotherapy -----------------------------------------
    def f_cis(v):
        p, reg = M.scenario("cis", 0.25)
        p["HZ0"] = fitted["HZ0"]; p["EMAX_CIS"] = v
        return M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)) or 99.0
    fitted["EMAX_CIS"] = bisect(f_cis, 9.3, 0.05, 0.40, label="EMAX_CIS")

    # --- 5. the doublet ----------------------------------------------------
    def f_pem(v):
        p, reg = M.scenario("pemcis", 0.25)
        p["HZ0"] = fitted["HZ0"]; p["EMAX_CIS"] = fitted["EMAX_CIS"]; p["EMAX_PEM"] = v
        return M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)) or 99.0
    fitted["EMAX_PEM"] = bisect(f_pem, 12.1, 0.005, 0.12, label="EMAX_PEM")

    # --- 6. bevacizumab: a DELIVERY effect, not a kill effect --------------
    base = f_pem(fitted["EMAX_PEM"])
    def f_bev(v):
        p, reg = M.scenario("pemcisbev", 0.25)
        p["HZ0"] = fitted["HZ0"]; p["EMAX_CIS"] = fitted["EMAX_CIS"]
        p["EMAX_PEM"] = fitted["EMAX_PEM"]; p["BEV_LAM"] = v
        return M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)) or 99.0
    fitted["BEV_LAM"] = bisect(f_bev, base + 2.7, 0.02, 1.00, label="BEV_LAM")

    # --- 7. checkpoint blockade -------------------------------------------
    def f_io(v):
        p, reg = M.scenario("nivoipi", 0.25)
        p["HZ0"] = fitted["HZ0"]; p["EMAX_IO"] = v
        return M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)) or 99.0
    fitted["EMAX_IO"] = bisect(f_io, base * 18.1 / 14.1, 0.02, 6.00, label="EMAX_IO")

    print()
    print(" fitted constants:")
    for k, v in fitted.items():
        print("   %-10s = %.6g" % (k, v))
    print()
    print(" %d constants fitted against %d endpoints; every other parameter in"
          % (len(fitted), len(fitted)))
    print(" default_params() is a literature value or a geometric fact.")
    print()
    return fitted


def apply_fit(p, fitted):
    for k, v in fitted.items():
        p[k] = v
    return p


def patch_defaults(fitted):
    """Monkey-patch default_params so every later check uses the fitted set."""
    orig = M.default_params
    def wrapped():
        return apply_fit(orig(), fitted)
    M.default_params = wrapped


# =============================================================================
#  VIRTUAL POPULATION: response RATE is a population quantity, not a trajectory
# =============================================================================
def lcg(seed=20260805):
    s = seed
    while True:
        s = (1103515245 * s + 12345) % (1 << 31)
        yield s / float(1 << 31)


def normal_stream(seed):
    u = lcg(seed)
    while True:
        u1 = max(next(u), 1e-12)
        u2 = next(u)
        r = math.sqrt(-2.0 * math.log(u1))
        yield r * math.cos(2 * math.pi * u2)
        yield r * math.sin(2 * math.pi * u2)


def virtual_population(name, n=40, seed=20260805, tmax=800.0, dt=0.05,
                       emt_mean=0.25, emt_sd=0.22):
    """Inter-individual variability on the four quantities that actually vary
    between patients: histology, growth rate, baseline burden and drug
    sensitivity.  Everything else is held at the typical value."""
    z = normal_stream(seed)
    os_list, resp_list, pfs_list = [], [], []
    for i in range(n):
        p, reg = M.scenario(name, 0.25)
        x = min(max(emt_mean + emt_sd * next(z), 0.0), 1.0)
        p["EMT"] = x
        p["KG0"] *= math.exp(0.35 * next(z))
        p["EMAX_CIS"] *= math.exp(0.45 * next(z))
        p["EMAX_PEM"] *= math.exp(0.45 * next(z))
        p["EMAX_IO"] *= math.exp(0.70 * next(z))
        p["HZ0"] *= math.exp(0.30 * next(z))
        h_par = max(0.25, 0.60 * math.exp(0.35 * next(z)))
        h_vis = max(0.15, 0.40 * math.exp(0.35 * next(z)))
        y0 = M.initial_state(p, h_par=h_par, h_vis=h_vis, h_fis=0.50)
        rec = M.simulate(p, tmax=tmax, dt=dt, regimen=reg, y0=y0)
        o = M.median_os(rec)
        os_list.append(o if o else tmax / 30.44)
        resp_list.append(M.best_response(rec, 12))
        pfs_list.append(M.pfs_months(rec))
    os_list.sort()
    med = os_list[len(os_list) // 2]
    orr = 100.0 * sum(1 for r in resp_list if r <= -30.0) / len(resp_list)
    dcr = 100.0 * sum(1 for r in resp_list if r <= 20.0) / len(resp_list)
    pfs_list.sort()
    return {"medOS": med, "ORR": orr, "DCR": dcr,
            "medPFS": pfs_list[len(pfs_list) // 2],
            "os": os_list, "resp": resp_list}


def population_table(n=40):
    print("=" * 78)
    print(" VIRTUAL POPULATION (n = %d per arm): response RATE is a population" % n)
    print(" quantity and cannot be read off a single trajectory")
    print("=" * 78)
    obs = {
        "cis":        ("Cisplatin alone (EMPHACIS ctrl)", 9.3, 16.7, 3.9),
        "pemcis":     ("Pemetrexed + cisplatin (EMPHACIS)", 12.1, 41.3, 5.7),
        "pemcisbev":  ("Pem + cis + bev (MAPS)", None, 44.4, None),
        "nivoipi":    ("Nivolumab + ipilimumab (CM-743)", None, 40.0, 6.8),
        "gemcis":     ("Gemcitabine + cisplatin", None, 33.0, None),
    }
    print("  %-34s %10s %10s   %8s %8s" % ("arm", "OS model", "OS obs", "ORR mod", "ORR obs"))
    out = {}
    for k, (lab, os_o, orr_o, pfs_o) in obs.items():
        r = virtual_population(k, n=n)
        out[k] = r
        print("  %-34s %8.1f mo %8s   %7.0f%% %7s"
              % (lab, r["medOS"], ("%.1f mo" % os_o) if os_o else "  -",
                 r["ORR"], ("%.0f%%" % orr_o) if orr_o else "  -"))
    print()
    return out


def main():
    fitted = calibrate()
    patch_defaults(fitted)

    M.check_geometry()
    M.check_penetration()
    M.check_trial_targets()
    M.check_crossover()
    M.check_checkmate743()
    M.check_collagen_bias()
    M.check_effusion_vs_trapping()
    M.check_surgery()
    M.check_intrapleural()
    M.check_toxicity()
    M.check_biomarker()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    population_table(n=n)
    M.check_mass_balance()

    print("=" * 78)
    print(" FITTED CONSTANTS TO COPY INTO mpm_reference_model.py / the R model")
    print("=" * 78)
    for k, v in fitted.items():
        print('   p["%s"] = %.6g' % (k, v))
    print()


if __name__ == "__main__":
    main()
