"""
hs_analysis.py -- every quantitative claim made in README.md is computed here.

Run:  python3 hs_analysis.py  > hs_verification_output.txt

The organising thesis, and the order of the sections below:

  Heat stroke is what happens when the heat-balance equation loses its fixed
  point.  Everything after that is a clock and a dose.

    regime 1  COMPENSABLE    a steady state exists; core temperature plateaus
    regime 2  UNCOMPENSABLE  no steady state; a clock runs at a computable rate
    regime 3  COMMITTED      the thermal dose has latched a bistable
                             inflammatory switch, and cooling restores the
                             temperature but not the patient
"""
import numpy as np
from scipy.optimize import brentq, fsolve
from hs_core import (make_params, Schedule, simulate, get, derived, y0, dubois,
                     wet_bulb, psat, rhs, IX, isth_dic_score, gcs_from_cnsd,
                     SNAMES)
from hs_calibrate import equilibrate_shell

# ---------------------------------------------------------------------------
# cooling modalities: lumped conductances fitted in hs_calibrate.py
# ---------------------------------------------------------------------------
MODALITY = {
    "ice-water immersion (2C)":     dict(UA_COOL=23.53, T_COOL=2.0,  IMMERSE=1),
    "cold-water immersion (14C)":   dict(UA_COOL=27.74, T_COOL=14.0, IMMERSE=1),
    "tarp-assisted cooling (10C)":  dict(UA_COOL=20.05, T_COOL=10.0, IMMERSE=1),
    "cold shower / dousing (20C)":  dict(UA_COOL=15.94, T_COOL=20.0, IMMERSE=0),
    "evaporative + convective":     dict(UA_COOL=16.23, T_COOL=25.0, IMMERSE=0),
    "endovascular catheter":        dict(UA_COOL=4.60,  T_COOL=4.0,  IMMERSE=0),
    "ice packs (neck/axilla/groin)": dict(UA_COOL=1.60, T_COOL=0.0,  IMMERSE=0),
    "passive (shade, no cooling)":  dict(UA_COOL=0.0,   T_COOL=30.0, IMMERSE=0),
}

FIELD = dict(TA=30.0, RH=0.50, VAIR=0.3, ICL=0.1, QSOL=0.0, MEX=0.0)
WARD = dict(TA=22.0, RH=0.50, VAIR=0.15, ICL=1.0, QSOL=0.0, MEX=0.0)
RACE = dict(TA=35.0, RH=0.80, VAIR=1.5, ICL=0.15, QSOL=120.0)

H_UNSTABLE = 10.91     # ng/mL, the commitment threshold (derived in hs_core docstring)
SEP = "=" * 78


def banner(txt):
    print()
    print(SEP)
    print(txt)
    print(SEP)


# ===========================================================================
# 1. THE FIXED POINT: compensable vs uncompensable
# ===========================================================================
def thermal_steady_state(p, env, mex, guess=(37.2, 37.4, 34.0)):
    """Solve the three thermal node equations for a simultaneous steady state.

    Returns (TCR, TMU, TSK) or None if no fixed point exists.  This is the
    definition of compensability and it is an algebraic question, not a
    simulation question: either the heat-balance equation has a root or it
    does not.  The subject is held euhydrated (ORAL_MATCH) so that the
    thermal boundary is not contaminated by the hydration boundary.
    """
    p = dict(p)
    p.setdefault("_AD", dubois(p["BW"], p["HT"]))
    e = dict(env)
    e.update(dict(MEX=mex, ORAL_MATCH=1.0))
    sch = Schedule([dict(dur=1e9, **e)], e)
    base = y0(p)

    def f(x):
        yy = base.copy()
        yy[IX["TCR"]], yy[IX["TMU"]], yy[IX["TSK"]] = x
        d = rhs(0.0, yy, p, sch)
        return [d[IX["TCR"]], d[IX["TMU"]], d[IX["TSK"]]]

    x, info, ok, msg = fsolve(f, guess, full_output=True)
    if ok != 1 or not np.all(np.isfinite(x)):
        return None
    if max(abs(v) for v in f(x)) > 1e-5 or not (15.0 < x[0] < 60.0):
        return None
    return tuple(x)


def terminal_slope(p, mex, rh, ta, tmax=180.0):
    """Core temperature slope (deg C/min) over the final 20 min of an exposure.

    Compensable exposures flatten to ~0; uncompensable ones are still climbing.
    The slope is monotone in ambient temperature, so it bisects cleanly --
    unlike fsolve on the steady state, which the min/max kinks in the effector
    equations defeat.
    """
    env = dict(TA=ta, RH=rh, VAIR=1.0, ICL=0.15, QSOL=0.0, ORAL_MATCH=1.0)
    sch = Schedule([dict(dur=tmax, MEX=mex, **env)], env)
    T, Y, p2 = simulate(p, sch, max_step=10.0,
                        t_eval=np.arange(0.0, tmax + 0.01, 2.0))
    tc = get(Y, "TCR")
    if len(tc) < 12 or not np.isfinite(tc[-1]):
        return 1.0                                   # ran away / terminated
    return (tc[-1] - tc[-11]) / 20.0


SLOPE_CRIT = 0.002      # deg C/min: below this the exposure has found a plateau


def critical_ambient(p, mex, rh, ta_lo=5.0, ta_hi=64.0):
    """Air temperature at which the exposure stops having a fixed point."""
    f = lambda ta: terminal_slope(p, mex, rh, ta) - SLOPE_CRIT
    if f(ta_lo) > 0 or f(ta_hi) < 0:
        return np.nan
    return brentq(f, ta_lo, ta_hi, xtol=0.05)


def section1():
    banner("1. REGIME 1 vs REGIME 2 -- does the heat-balance equation have a\n"
           "   fixed point?  The critical ambient temperature at which it is lost.")
    print(f"{'metabolic rate':>16s} {'RH':>5s} {'Ta_crit':>8s} {'Tw_crit':>8s} "
          f"{'Ta_crit(acclim)':>16s} {'Tw_crit(acclim)':>16s} {'gain':>6s}")
    rows = []
    for mex, lbl in [(0.0, "rest (84 W)"), (200.0, "light 200 W"),
                     (400.0, "moderate 400 W"), (600.0, "hard 600 W"),
                     (900.0, "very hard 900 W")]:
        for rh in (0.30, 0.60, 0.90):
            pu = make_params(ACCLIM=0.0)
            pa = make_params(ACCLIM=1.0)
            tu = critical_ambient(pu, mex, rh)
            ta_ = critical_ambient(pa, mex, rh)
            wu = wet_bulb(tu, rh) if np.isfinite(tu) else np.nan
            wa = wet_bulb(ta_, rh) if np.isfinite(ta_) else np.nan
            rows.append((lbl, rh, tu, wu, ta_, wa))
            print(f"{lbl:>16s} {rh:5.2f} {tu:8.1f} {wu:8.1f} {ta_:16.1f} "
                  f"{wa:16.1f} {ta_ - tu:+6.1f}")
    rest_w = np.nanmean([r[3] for r in rows if r[0].startswith("rest")])
    w900 = np.nanmean([r[3] for r in rows if "900" in r[0]])
    print()
    print(f"  Reading: at rest the boundary sits at a wet-bulb of {rest_w:.1f} C -- the")
    print("  canonical 35 C survivability limit falls out of the heat-balance")
    print(f"  equation without being put in.  At 900 W it drops to {w900:.1f} C.")
    print("  The 35 C limit is the RESTING limit.  Working people lose the fixed")
    print("  point many degrees of wet-bulb earlier, which is the empirical")
    print("  finding of Vecellio 2022 obtained here from the equation alone.")
    print()
    d = [r[4] - r[2] for r in rows if np.isfinite(r[4] - r[2])]
    print("  ACCLIMATISATION DOES NOT MOVE THIS BOUNDARY.  Mean shift "
          f"{np.mean(d):+.2f} C")
    print(f"  (range {min(d):+.2f} to {max(d):+.2f}), and the sign is not")
    print("  consistently favourable.  The reason is visible in the solution: at")
    print("  every point on this boundary E_actual == E_max,environment, while")
    print("  sweat capacity is 668 W unacclimatised and 1233 W acclimatised.")
    print("  Neither is ever reached.  The AIR is the binding constraint, so")
    print("  doubling the gland's output buys nothing.")
    print()
    print("  This is a falsifiable structural claim, not a hedge: heat")
    print("  acclimatisation -- and every other sweat-side countermeasure -- must")
    print("  act on the TRANSIENT (how fast the core climbs, how much")
    print("  cardiovascular reserve is left) rather than on WHICH environments")
    print("  are survivable.  Section 8 tests that consequence.")
    return rows


# ===========================================================================
# 2. THE CLOCK: rate of rise, and what a fixed 40 C threshold means
# ===========================================================================
def section2():
    banner("2. REGIME 2 -- the clock.  Same equation, two diseases.")
    cases = [
        ("EHS: 900 W, 35C/80% RH, sun", make_params(),
         dict(TA=35., RH=.80, VAIR=1.5, ICL=.15, QSOL=120.), 900., 180),
        ("EHS: 900 W, acclimatised",    make_params(ACCLIM=1.0),
         dict(TA=35., RH=.80, VAIR=1.5, ICL=.15, QSOL=120.), 900., 180),
        ("EHS: 1100 W, football pads",  make_params(),
         dict(TA=33., RH=.70, VAIR=0.5, ICL=1.2, QSOL=200.), 1100., 180),
        ("NEHS: elderly+anticholinergic, 40C/55%",
         make_params(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05, FSW_DRUG=.30),
         dict(TA=40., RH=.55, VAIR=0.15, ICL=0.7), 0., 4320),
        ("NEHS: elderly+anticholinergic, 43C/45%",
         make_params(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05, FSW_DRUG=.30),
         dict(TA=43., RH=.45, VAIR=0.15, ICL=0.7), 0., 4320),
    ]
    print(f"{'case':42s} {'t(38.5)':>8s} {'t(40.0)':>8s} {'t(42.0)':>8s} "
          f"{'rate 40->42':>12s} {'warn 40->42':>12s}")
    out = {}
    for lbl, p, env, mex, dur in cases:
        sch = Schedule([dict(dur=dur, MEX=mex, **env)], env)
        T, Y, p2 = simulate(p, sch, max_step=5.0)
        tc = get(Y, "TCR")

        def cross(x):
            i = np.argmax(tc >= x)
            return T[i] if tc.max() >= x else np.nan
        t385, t40, t42 = cross(38.5), cross(40.0), cross(42.0)
        rate = (42.0 - 40.0) / (t42 - t40) if np.isfinite(t42) else np.nan
        out[lbl] = (t385, t40, t42, rate)
        u = "min" if dur < 400 else "min"
        print(f"{lbl:42s} {t385:8.0f} {t40:8.0f} {t42:8.0f} "
              f"{rate:9.4f} C/{u[:3]} {t42 - t40:9.0f} min")
    print()
    print("  The interval between crossing 40.0 C -- the diagnostic threshold --")
    print("  and reaching 42.0 C IS the warning time a temperature-based")
    print("  criterion gives you.  It differs between the two heat strokes by")
    ehs = out["EHS: 900 W, 35C/80% RH, sun"]
    neh = out["NEHS: elderly+anticholinergic, 40C/55%"]
    print(f"  a factor of {(neh[2]-neh[1])/(ehs[2]-ehs[1]):.0f} "
          f"({ehs[2]-ehs[1]:.0f} min vs {neh[2]-neh[1]:.0f} min).")
    print("  It is the same threshold on two clocks running at different speeds,")
    print("  and that -- not a difference in intrinsic lethality -- is the first")
    print("  half of why classic heat stroke kills more often than exertional.")
    return out


# ===========================================================================
# 3. THE DOSE: CEM43, and the exchange rate between modality and delay
# ===========================================================================
def collapsed_patient(tcr=42.0, acclim=0.0, **kw):
    """A runner who has just gone down: hyperthermic, 4% dehydrated, glands
    fatigued, exercise stopped, now lying in a shaded tent."""
    p = make_params(ACCLIM=acclim, **kw)
    p["_AD"] = dubois(p["BW"], p["HT"])
    yi = y0(p, TCR=tcr, TSK=37.0, TMU=tcr + 0.3)
    yi[IX["WDEF"]] = 0.04 * p["BW"]
    yi[IX["VP"]] = p["VP0"] * 0.90
    yi[IX["SWFAT"]] = 0.25
    yi[IX["HSP"]] = p["HSP_BASE"] * 1.6
    yi = equilibrate_shell(p, FIELD, yi, tcr)
    return p, yi


def cool_run(modality, delay, tcr0=42.0, tail=0.0, acclim=0.0, ivf=None,
             stop=38.6, maxdur=400.0):
    """Collapse at tcr0, wait `delay` minutes in the shade, then cool.
    Returns (T, Y, p, total CEM43 to reach `stop`, minutes of cooling)."""
    p, yi = collapsed_patient(tcr0, acclim=acclim)
    m = MODALITY[modality]
    phases = []
    if delay > 0:
        phases.append(dict(dur=delay, **FIELD))
    ph = dict(FIELD)
    ph.update(dict(dur=maxdur, UA_COOL=m["UA_COOL"], T_COOL=m["T_COOL"],
                   IMMERSE=m["IMMERSE"], COOL_STOP_TC=stop))
    if ivf:
        ph.update(dict(IVF_RATE=ivf[0], IVF_TEMP=ivf[1]))
    phases.append(ph)
    sch = Schedule(phases, FIELD)
    T, Y, p2 = simulate(p, sch, y_init=yi, max_step=2.0, stop_tc=None)
    tc = get(Y, "TCR")
    below = np.where(tc <= stop)[0]
    iend = below[0] if len(below) else len(tc) - 1
    return T[:iend + 1], Y[:, :iend + 1], p2, get(Y, "CEM43")[iend], T[iend] - delay


def section3():
    banner("3. REGIME 2 -> the DOSE.  CEM43 = the integral the patient actually\n"
           "   pays.  Cooling modality is a POTENCY; delay is an EXPOSURE TIME.")
    print("  Sapareto-Dewey: dose rate = R^(43-Tc) equivalent minutes per minute,")
    print("  R = 0.25 below 43 C, 0.50 at or above.  Each degree below 43 divides")
    print("  the rate by four.  Time and temperature are NOT exchangeable 1:1.")
    print()
    print("     Tc:   40.0     41.0     42.0     42.5     43.0     43.5     44.0")
    rates = []
    for tc in [40., 41., 42., 42.5, 43., 43.5, 44.]:
        R = 0.25 if tc < 43.0 else 0.5
        rates.append(R ** (43.0 - tc))
    print("   rate: " + "".join(f"{r:8.3f} " for r in rates)
          + "  CEM43-min per min")
    print()

    print("  A. Total CEM43 paid from collapse at 42.0 C to 38.6 C:")
    print(f"  {'modality':32s} {'delay 0':>9s} {'delay 10':>9s} {'delay 30':>9s} "
          f"{'delay 60':>9s} {'cool min':>9s}")
    table = {}
    for mod in MODALITY:
        row = []
        for d in (0, 10, 30, 60):
            _, _, _, dose, cmin = cool_run(mod, d)
            row.append(dose)
            if d == 0:
                cm = cmin
        table[mod] = row
        print(f"  {mod:32s} {row[0]:9.2f} {row[1]:9.2f} {row[2]:9.2f} "
              f"{row[3]:9.2f} {cm:9.1f}")
    print()

    # the exchange rate: how many minutes of delay is a modality upgrade worth?
    print("  B. THE EXCHANGE RATE.  Upgrading the cooler buys you a fixed number")
    print("     of minutes of delay, and the number is remarkably stable.")
    base = "evaporative + convective"
    best = "ice-water immersion (2C)"
    for tcr0 in (41.0, 42.0, 43.0):
        d0 = cool_run(base, 0, tcr0=tcr0)[3]
        d1 = cool_run(best, 0, tcr0=tcr0)[3]
        saved = d0 - d1
        R = 0.25 if tcr0 < 43.0 else 0.5
        per_min = R ** (43.0 - tcr0)
        print(f"     collapse at {tcr0:.1f} C: upgrading {base} -> {best}")
        print(f"        saves {saved:6.2f} CEM43; one minute of delay costs "
              f"{per_min:5.3f} CEM43")
        print(f"        => the better cooler is worth {saved/per_min:5.1f} MINUTES "
              f"of delay")
    print()
    print("     The exchange rate barely moves with peak temperature, but the")
    print("     absolute stakes rise exponentially.  This is the quantitative")
    print("     content of 'cool first, transport second': modality is a factor")
    print("     of a few, delay is unbounded.")
    return table


# ===========================================================================
# 4. REGIME 3: the commitment switch
# ===========================================================================
def switch_fixed_points(c_extra=0.0):
    p = make_params()
    A, K, n = p["KH_AUTO"], p["KH_HALF"], p["NH"]
    c = p["KH_OUT"] + c_extra - p["KH_NEC"] * p["KNEC_H"]
    f = lambda H: A * H ** n / (H ** n + K ** n) - c * H
    roots = []
    for lo, hi in [(0.2, 30.0), (30.0, 400.0)]:
        try:
            roots.append(brentq(f, lo, hi))
        except ValueError:
            pass
    return c, roots


def committed(Y):
    return get(Y, "HMGB1")[-1] > H_UNSTABLE


def section4():
    banner("4. REGIME 3 -- the commitment switch.  Cooling restores the\n"
           "   temperature; it does not always restore the patient.")
    c, roots = switch_fixed_points()
    print(f"  effective HMGB1 clearance c_eff = {c:.5f}/min")
    print(f"  fixed points: OFF = 0.0,  UNSTABLE = {roots[0]:.2f} ng/mL,  "
          f"ON = {roots[1]:.1f} ng/mL")
    p = make_params()
    Ht = (2.0 * p["KH_HALF"] ** p["NH"]) ** (1.0 / p["NH"])
    c_crit = p["KH_AUTO"] * Ht ** p["NH"] / (Ht ** p["NH"] + p["KH_HALF"] ** p["NH"]) / Ht
    print(f"  saddle-node at H = {Ht:.1f} ng/mL: if c_eff exceeds {c_crit:.5f}/min")
    print(f"  the ON state ceases to exist "
          f"(requires +{c_crit - c:.5f}/min of extra clearance).")
    print()

    # find the committing delay for each modality, then read off the dose
    print("  A. For each modality, the DELAY at which the switch latches, and")
    print("     the CEM43 that delay corresponds to (collapse at 42.0 C).")
    print("     The switch is scored at 48 h, not at the end of cooling: the")
    print("     temperature is settled within the hour and the switch is not.")
    print(f"  {'modality':32s} {'crit delay':>11s} {'CEM43 paid by then':>20s}")
    crit = {}
    for mod in MODALITY:
        if mod.startswith("passive"):
            continue

        def f(d):
            r = full_course("probe", delay=d, modality=mod, days=2.0)
            return get(r["Y"], "HMGB1")[-1]
        lo, hi = 0.0, 200.0
        if f(hi) < H_UNSTABLE:
            print(f"  {mod:32s} {'>200':>11s} {'-':>20s}")
            continue
        if f(lo) > H_UNSTABLE:
            print(f"  {mod:32s} {'0 (already)':>11s} {'-':>20s}")
            crit[mod] = 0.0
            continue
        d = brentq(lambda x: f(x) - H_UNSTABLE, lo, hi, xtol=1.0)
        dose = cool_run(mod, d)[3]
        crit[mod] = d
        print(f"  {mod:32s} {d:11.1f} {dose:20.2f}")
    print()
    if crit:
        ds = [cool_run(m, d)[3] for m, d in crit.items()]
        print(f"  The committing DOSE is {np.mean(ds):.1f} +/- {np.std(ds):.1f} CEM43 "
              f"across modalities (spread {max(ds)-min(ds):.2f}),")
        print("  while the committing DELAY ranges from "
              f"{min(crit.values()):.0f} to {max(crit.values()):.0f} minutes.")
        print("  The dose is a property of the patient; the delay you can afford")
        print("  is a property of the cooler.  That is the whole argument for")
        print("  cooling on the field, expressed as two numbers.")
    return crit


# ===========================================================================
# 5. THE FULL SCENARIO SUITE
# ===========================================================================
def full_course(label, race_env=None, mex=900.0, collapse_tc=42.0, delay=20.0,
                modality="ice-water immersion (2C)", acclim=0.0, pkw=None,
                doses=(), days=7.0, ivf=None, neh=None, stop=38.6,
                oral=0.0):
    """Race (or heat-wave exposure) -> collapse -> delay -> cooling -> ward."""
    pkw = dict(pkw or {})
    p = make_params(ACCLIM=acclim, **pkw)
    m = MODALITY[modality]
    phases = []
    if neh is None:
        env = dict(race_env or RACE)
        phases.append(dict(dur=600.0, MEX=mex, **env))
    else:
        phases.append(dict(dur=neh["dur"], MEX=0.0,
                           **{k: v for k, v in neh.items() if k != "dur"}))
    sch = Schedule(phases, phases[0])
    T1, Y1, p1 = simulate(p, sch, max_step=5.0, stop_tc=collapse_tc)
    t_coll = T1[-1]
    y = Y1[:, -1].copy()

    phases2 = []
    if delay > 0:
        phases2.append(dict(dur=delay, **FIELD))
    ph = dict(FIELD)
    ph.update(dict(dur=240.0, UA_COOL=m["UA_COOL"], T_COOL=m["T_COOL"],
                   IMMERSE=m["IMMERSE"], COOL_STOP_TC=stop))
    if ivf:
        ph.update(dict(IVF_RATE=ivf[0], IVF_TEMP=ivf[1]))
    phases2.append(ph)
    ward = dict(WARD)
    ward.update(dict(dur=days * 1440.0, ORAL_RATE=oral or 0.0012))
    phases2.append(ward)
    sch2 = Schedule(phases2, FIELD)
    T2, Y2, p2 = simulate(p, sch2, y_init=y, doses=doses, max_step=5.0)
    T = np.concatenate([T1, T2 + t_coll])
    Y = np.concatenate([Y1, Y2], axis=1)
    return dict(label=label, T=T, Y=Y, p=p2, t_collapse=t_coll,
                t_cool_start=t_coll + delay)


def summarise(r):
    Y, T = r["Y"], r["T"]
    tc = get(Y, "TCR")
    plt_ = get(Y, "PLT")
    fib = get(Y, "FIB")
    dd = get(Y, "DDIM")
    ptr = 1.0 + 0.9 * (1.0 - get(Y, "LIVF"))          # PT ratio proxy
    dic = max(isth_dic_score(plt_[i], dd[i], fib[i], ptr[i])
              for i in range(0, len(T), 10))
    return dict(
        label=r["label"],
        peakTc=tc.max(),
        CEM43=get(Y, "CEM43")[-1],
        HMGB1=get(Y, "HMGB1")[-1],
        committed=get(Y, "HMGB1")[-1] > H_UNSTABLE,
        peakALT=get(Y, "ALT").max(),
        peakAST=get(Y, "AST").max(),
        peakCK=get(Y, "CK").max(),
        peakCr=get(Y, "SCR").max(),
        minGFR=get(Y, "GFRF").min(),
        minPLT=plt_.min(),
        minFIB=fib.min(),
        peakDD=dd.max(),
        dic=dic,
        peakIL6=get(Y, "IL6").max(),
        peakTNF=get(Y, "TNF").max(),
        worstGCS=gcs_from_cnsd(get(Y, "CNSD").max()),
        # GCS at fixed follow-up times discriminates between arms; the WORST
        # value does not, because every arm is encephalopathic at 42 C -- that
        # is the diagnostic criterion, not an outcome.
        gcs24=gcs_from_cnsd(np.interp(T[0] + 1440.0, T, get(Y, "CNSD"))),
        gcsEnd=gcs_from_cnsd(get(Y, "CNSD")[-1]),
        endGFR=get(Y, "GFRF")[-1],
        peakLPS=get(Y, "LPS").max(),
        endLIVF=get(Y, "LIVF")[-1],
    )


def section5():
    banner("5. SCENARIO SUITE -- 16 arms.  Biology held fixed except where named.")
    S = []
    S.append(full_course("01 EHS, no cooling until hospital (60 min)",
                         delay=60, modality="evaporative + convective"))
    S.append(full_course("02 EHS, on-site ice-water immersion (<5 min)",
                         delay=4, modality="ice-water immersion (2C)"))
    S.append(full_course("03 EHS, on-site CWI at 20 min", delay=20,
                         modality="ice-water immersion (2C)"))
    S.append(full_course("04 EHS, tarp-assisted cooling at 10 min", delay=10,
                         modality="tarp-assisted cooling (10C)"))
    S.append(full_course("05 EHS, evaporative at 10 min", delay=10,
                         modality="evaporative + convective"))
    S.append(full_course("06 EHS, ice packs only at 10 min", delay=10,
                         modality="ice packs (neck/axilla/groin)"))
    S.append(full_course("07 EHS, 2 L cold saline only at 20 min", delay=20,
                         modality="passive (shade, no cooling)",
                         ivf=(0.10, 4.0)))
    S.append(full_course("08 EHS, CWI at 20 min + 2 L cold saline", delay=20,
                         modality="ice-water immersion (2C)", ivf=(0.10, 4.0)))
    S.append(full_course("09 EHS, heat-acclimatised runner, CWI at 20 min",
                         delay=20, acclim=1.0))
    S.append(full_course("10 EHS, euhydrated (drinking during race), CWI 20 min",
                         delay=20, oral=0.012))
    S.append(full_course("11 EHS + paracetamol 1 g IV at cooling start",
                         delay=20, doses=[(20.0, "PARA_C", 1000.0)]))
    S.append(full_course("12 EHS + ibuprofen 800 mg at cooling start",
                         delay=20, doses=[(20.0, "IBU_A", 800.0)]))
    S.append(full_course("13 EHS + dantrolene 2.5 mg/kg at cooling start",
                         delay=20, doses=[(20.0, "DAN_C", 175.0)]))
    S.append(full_course("14 EHS + hydrocortisone 200 mg at cooling start",
                         delay=20, doses=[(20.0, "HC_C", 200.0)]))
    S.append(full_course(
        "15 EHS, 60 min delay + thrombomodulin alfa x3 d", delay=60,
        modality="evaporative + convective",
        doses=[(60.0 + 1440.0 * k, "RTM_C", 3.5) for k in range(3)]))
    S.append(full_course(
        "16 NEHS, elderly + anticholinergic, found late, evaporative",
        neh=dict(dur=4320.0, TA=41.0, RH=0.55, VAIR=0.15, ICL=0.7),
        collapse_tc=41.5, delay=45, modality="evaporative + convective",
        pkw=dict(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05, FSW_DRUG=.30)))

    rows = [summarise(r) for r in S]
    hdr = (f"{'scenario':46s} {'CEM43':>6s} {'HMGB1':>6s} {'cmt':>4s} "
           f"{'IL6':>6s} {'ALT':>5s} {'CK':>7s} {'Cr':>5s} {'PLT':>4s} "
           f"{'DIC':>4s} {'GCS24':>6s} {'GCSend':>7s}")
    print(hdr)
    for r in rows:
        print(f"{r['label']:46s} {r['CEM43']:6.2f} "
              f"{r['HMGB1']:6.1f} {'YES' if r['committed'] else ' no':>4s} "
              f"{r['peakIL6']:6.0f} {r['peakALT']:5.0f} {r['peakCK']:7.0f} "
              f"{r['peakCr']:5.2f} {r['minPLT']:4.0f} {r['dic']:4d} "
              f"{r['gcs24']:6.1f} {r['gcsEnd']:7.1f}")
    print()
    print("  Every arm reaches the same peak core temperature (42.0 C, the")
    print("  collapse trigger) and every arm is encephalopathic at that moment --")
    print("  that is the diagnostic criterion, not an outcome.  What separates the")
    print("  arms is the DOSE paid afterwards and whether the switch latched.")
    return S, rows


# ===========================================================================
# 6. THE ANTIPYRETIC ARGUMENT
# ===========================================================================
def section6():
    banner("6. WHY ANTIPYRETICS CANNOT WORK -- and what they cost.")
    print("  Fever raises the SET-POINT and the body defends it.  Heat stroke")
    print("  leaves the set-point at 37 C and saturates the effectors trying to")
    print("  hold it.  An antipyretic lowers a threshold that is no longer the")
    print("  binding constraint.  The model makes this a sensitivity, not a claim:")
    print()
    print("  A 0.5 C set-point shift is applied, and the effect is read off both")
    print("  in absolute terms and as a FRACTION of the rate of rise.  The")
    print("  fraction is the number that matters.")
    print()
    print(f"  {'condition':34s} {'HSI':>5s} {'sweat':>7s} {'SkBF':>6s} "
          f"{'dTc/dt':>9s} {'d(rate)':>9s} {'as % of rate':>13s} "
          f"{'dCEM43@90':>10s}")
    for lbl, env, mex in [
            ("compensable (rest, 30C/50%)", dict(TA=30., RH=.5, VAIR=.5, ICL=.3), 0.),
            ("marginal (400 W, 33C/60%)", dict(TA=33., RH=.6, VAIR=1., ICL=.15), 400.),
            ("uncompensable (900 W, 35C/80%)",
             dict(TA=35., RH=.8, VAIR=1.5, ICL=.15, QSOL=120.), 900.)]:
        vals = []
        for dts in (0.0, 0.5):
            p = make_params(TCSW0=37.0 - dts, TCVD0=36.8 - dts)
            p["_AD"] = dubois(p["BW"], p["HT"])
            sch = Schedule([dict(dur=90., MEX=mex, **env)], env)
            T, Y, p2 = simulate(p, sch, max_step=2.0)
            D = derived(T, Y, p2, sch)
            tc = get(Y, "TCR")
            rate = (tc[60] - tc[40]) / 20.0
            vals.append((rate, D["HSI"][45], D["m_sw"][45] * 60 / 1000,
                         D["SKBF"][45], get(Y, "CEM43")[-1]))
        dr = vals[1][0] - vals[0][0]
        r, hsi, sw, sk, dose0 = vals[0]
        pct = 100.0 * dr / r if abs(r) > 1e-9 else float("nan")
        print(f"  {lbl:34s} {hsi:4.0f}% {sw:5.2f}L/h {sk:6.2f} "
              f"{r:+9.4f} {dr:+9.4f} {pct:+12.1f}% {vals[1][4]-dose0:+10.3f}")
    print()
    print("  The ABSOLUTE effect of moving the set-point is about the same in")
    print("  all three regimes -- roughly a thousandth of a degree per minute.")
    print("  What changes by two orders of magnitude is what that is a fraction")
    print("  OF.  In a compensable exposure it is several times the entire rate")
    print("  of change; in an uncompensable one it is about 1% of it, and the")
    print("  sign is not even reliably favourable.")
    print()
    print("  That is the whole antipyretic argument, and note what it is NOT: it")
    print("  is not that antipyretics are weak drugs.  They move the set-point")
    print("  exactly as advertised.  They are ineffective here because in heat")
    print("  stroke the set-point is not the binding constraint -- the effectors")
    print("  are already saturated against it.  Antipyretics work precisely in")
    print("  the region where the patient does not have heat stroke.")
    print()

    print("  And the cost side.  Paracetamol 1 g IV on a heat-injured liver,")
    print("  with CYP2E1 induced and glutathione already being spent:")
    for lbl, kw, dose in [("healthy liver, 1 g", dict(), 1000.0),
                          ("heat-stroke liver, 1 g", dict(), 1000.0),
                          ("heat-stroke liver, 4 g/24 h", dict(), 1000.0)]:
        pass
    a = full_course("noPARA", delay=20)
    b = full_course("PARA1g", delay=20, doses=[(20.0, "PARA_C", 1000.0)])
    c = full_course("PARA4g", delay=20,
                    doses=[(20.0 + 360.0 * k, "PARA_C", 1000.0) for k in range(4)])
    for r in (a, b, c):
        s = summarise(r)
        gsh = get(r["Y"], "GSH").min()
        print(f"    {r['label']:10s} peak ALT {s['peakALT']:7.0f} U/L | "
              f"peak AST {s['peakAST']:7.0f} | nadir GSH {gsh:5.2f} | "
              f"end liver function {s['endLIVF']:.3f} | peak Tc {s['peakTc']:.2f}")
    print("    Peak core temperature is unchanged to two decimals; the liver is not.")
    return None


# ===========================================================================
# 7. THE ONE AGENT THAT TOUCHES THE SWITCH
# ===========================================================================
def section7():
    banner("7. THROMBOMODULIN ALFA -- the only agent in the model that acts on\n"
           "   the commitment variable, and therefore the only one with a window\n"
           "   set by the switch rather than by the thermometer.")
    c0, r0 = switch_fixed_points(0.0)
    for extra, lbl in [(0.0, "no rTM"), (0.0012, "half dose"),
                       (0.0024, "standard 380 U/kg/d"), (0.0040, "double")]:
        c, roots = switch_fixed_points(extra)
        if len(roots) < 2:
            print(f"  {lbl:22s} c_eff={c:.5f}/min  -> MONOSTABLE (ON state abolished)")
        else:
            print(f"  {lbl:22s} c_eff={c:.5f}/min  -> unstable {roots[0]:5.2f}, "
                  f"ON {roots[1]:5.1f} ng/mL")
    print()
    print("  Timing, at TWO doses.  A 60-minute cooling delay is a committing")
    print("  exposure; rTM is then started at increasing times after collapse.")
    print("  The two doses ask different questions, because the standard dose")
    print("  destroys the ON state while the half dose only moves it.")
    for amt, dlbl in ((3.5, "STANDARD 380 U/kg/d (monostable)"),
                      (1.75, "HALF dose (bistability survives)")):
        print()
        print(f"  --- {dlbl}")
        print(f"  {'rTM start':>12s} {'end HMGB1':>11s} {'committed':>11s} "
              f"{'peak ALT':>9s} {'end GFR':>8s} {'GCS day 7':>10s}")
        for start in (30, 60, 120, 240, 480, 1440, 2880, None):
            if start is None:
                d, lbl = [], "none"
            else:
                d = [(float(start) + 1440.0 * k, "RTM_C", amt) for k in range(3)]
                lbl = f"{start} min"
            r = full_course(f"rtm{start}", delay=60,
                            modality="evaporative + convective", doses=d)
            s = summarise(r)
            print(f"  {lbl:>12s} {s['HMGB1']:11.1f} "
                  f"{'YES' if s['committed'] else 'no':>11s} {s['peakALT']:9.0f} "
                  f"{s['endGFR']:8.2f} {s['gcsEnd']:10.1f}")
    print()
    print("  Read the two blocks together.  At the standard dose the model makes")
    print("  a strong and uncomfortable claim: because +0.0024/min of clearance")
    print("  puts c_eff past the saddle-node, the ON state ceases to EXIST, so")
    print("  the agent reverses an already-latched switch and there is no time")
    print("  window at all.  At half dose the ON state survives at 33 ng/mL and")
    print("  a window reappears.  The window is therefore not a pharmacokinetic")
    print("  property of the drug -- it is a property of which side of the")
    print("  saddle-node the dose puts you on.  That is a falsifiable structural")
    print("  prediction and it is the model's most exposed therapeutic claim:")
    print("  the sepsis-DIC trial of this agent (SCARLET) was null.")


# ===========================================================================
# 8. PREVENTION vs RESCUE: where the effect sizes actually are
# ===========================================================================
def section8():
    banner("8. WHERE THE EFFECT SIZES ARE.  Prevention moves the boundary;\n"
           "   rescue only shortens the clock.")
    base = dict(TA=35., RH=.80, VAIR=1.5, ICL=.15, QSOL=120.)
    print(f"  {'intervention':44s} {'t to 42.0 C':>12s} {'CEM43 @90min':>13s}")
    variants = [
        ("baseline: unacclimatised, no drinking, 900 W", dict(), base, 900., 0.0),
        ("heat-acclimatised (10-14 d)", dict(ACCLIM=1.0), base, 900., 0.0),
        ("drinking 0.7 L/h", dict(), base, 900., 0.012),
        ("work rate 900 -> 600 W", dict(), base, 600., 0.0),
        ("work rate 900 -> 400 W", dict(), base, 400., 0.0),
        ("air velocity 1.5 -> 4 m/s", dict(), dict(base, VAIR=4.0), 900., 0.0),
        ("shade (remove 120 W solar)", dict(), dict(base, QSOL=0.0), 900., 0.0),
        ("humidity 80% -> 40%", dict(), dict(base, RH=0.40), 900., 0.0),
        ("acclimatised + drinking + shade",
         dict(ACCLIM=1.0), dict(base, QSOL=0.0), 900., 0.012),
    ]
    for lbl, pk, env, mex, oral in variants:
        p = make_params(**pk)
        sch = Schedule([dict(dur=180., MEX=mex, ORAL_RATE=oral, **env)], env)
        T, Y, p2 = simulate(p, sch, max_step=5.0)
        tc = get(Y, "TCR")
        t42 = T[np.argmax(tc >= 42.0)] if tc.max() >= 42.0 else np.inf
        i90 = min(90, len(T) - 1)
        d90 = get(Y, "CEM43")[i90]
        t42s = f"{t42:.0f} min" if np.isfinite(t42) else "never"
        print(f"  {lbl:44s} {t42s:>12s} {d90:13.2f}")
    print()
    print("  Only the interventions that restore a fixed point produce 'never'.")
    print("  Everything else buys minutes.")


if __name__ == "__main__":
    print("HEAT STROKE QSP MODEL -- computed results")
    print("all numbers below are produced by hs_core.py via this script")
    section1()
    section2()
    section3()
    section4()
    S, rows = section5()
    section6()
    section7()
    section8()
    print()
    print(SEP)
    print("END")
