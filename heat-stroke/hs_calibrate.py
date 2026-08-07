"""
hs_calibrate.py -- fit the lumped cooling conductance UA_COOL for each modality
to published whole-body core cooling rates, and check the thermal core against
independent physiological anchors.

Targets (rectal/oesophageal core cooling rate, deg C/min), taken from the
cooling-modality literature summarised in hs_references.md:

    ice-water immersion (2 C, whole body)   0.22    Casa 2007; Zhang 2015
    cold-water immersion (14 C)             0.17    Zhang 2015; Proulx 2003
    tarp-assisted cooling, ~10 C water      0.14    Hosokawa 2017; Luhring 2016
    cold shower / dousing (20 C)            0.10    Butts 2016
    evaporative + convective (fan + mist)   0.08    Weiner 1980; Gaudio 2016
    endovascular catheter                   0.06    Hoedemaekers 2007
    ice packs to neck/axilla/groin          0.032   McDermott 2009
    passive rest in shade                   0.022   McDermott 2009

The fit is done at Tc = 41.0 C in a standard reference patient (collapsed,
resting, 30 C shaded ambient, 50% RH, still air, unclothed, 3% body-mass water
deficit, partially fatigued sweat glands), which is the state in which the
published rates were actually measured.
"""
import numpy as np
from scipy.optimize import brentq, fsolve
from hs_core import (make_params, Schedule, simulate, get, y0, dubois, IX, rhs)


def equilibrate_shell(p, env, y, tcr):
    """Solve for the muscle and skin temperatures that are in steady state with a
    core held at `tcr`, leaving every other state untouched.

    Without this the first minutes of any cooling run are dominated by the
    core charging or discharging the skin node's heat capacity, which is an
    artefact of the initial condition and not a cooling rate.
    """
    sch = Schedule([dict(dur=1e6, **{k: v for k, v in env.items()})], env)

    def f(x):
        yy = np.array(y, dtype=float)
        yy[IX["TCR"]] = tcr
        yy[IX["TMU"]], yy[IX["TSK"]] = x
        d = rhs(0.0, yy, p, sch)
        return [d[IX["TMU"]], d[IX["TSK"]]]

    sol = fsolve(f, [tcr + 0.3, tcr - 1.0], full_output=True)
    x, _, ok, msg = sol
    if ok != 1:
        raise RuntimeError(f"shell equilibration failed: {msg}")
    y = np.array(y, dtype=float)
    y[IX["TCR"]] = tcr
    y[IX["TMU"]], y[IX["TSK"]] = x
    return y

TARGETS = [
    #  label                         Tmedium  immerse  target rate
    ("ice_water_immersion",              2.0,   1.0,   0.220),
    ("cold_water_immersion_14C",        14.0,   1.0,   0.170),
    ("tarp_assisted_cooling",           10.0,   1.0,   0.140),
    ("cold_shower_dousing",             20.0,   0.0,   0.100),
    ("evaporative_convective",          25.0,   0.0,   0.080),
    ("endovascular_catheter",            4.0,   0.0,   0.060),
    ("ice_packs_neck_axilla_groin",      0.0,   0.0,   0.032),
    ("passive_shade",                   30.0,   0.0,   0.022),
]

# shaded field tent, light shorts, near-still air.  This environment is chosen
# so that the PASSIVE arm reproduces the observed 0.022 C/min, which is the
# comparator against which every published modality rate was measured.
REF_ENV = dict(TA=30.0, RH=0.50, VAIR=0.3, ICL=0.1, QSOL=0.0, MEX=0.0)


def reference_patient():
    p = make_params()
    p["_AD"] = dubois(p["BW"], p["HT"])
    yi = y0(p, TCR=41.0, TSK=36.5, TMU=41.3)
    yi[IX["WDEF"]] = 0.03 * p["BW"]          # 3% body-mass deficit
    yi[IX["VP"]] = p["VP0"] * 0.93
    yi[IX["SWFAT"]] = 0.18
    yi = equilibrate_shell(p, REF_ENV, yi, 41.0)
    return p, yi


def cooling_rate(ua, tmed, immerse, minutes=10.0):
    p, yi = reference_patient()
    sch = Schedule([dict(dur=minutes, MEX=0.0, UA_COOL=ua, T_COOL=tmed,
                         IMMERSE=immerse, COOL_STOP_TC=30.0)], REF_ENV)
    T, Y, p = simulate(p, sch, y_init=yi, t_eval=np.arange(0, minutes + 0.01, 0.5))
    tc = get(Y, "TCR")
    return (tc[0] - tc[-1]) / minutes


def fit():
    out = {}
    print(f"{'modality':32s} {'T_med':>6s} {'imm':>4s} {'target':>7s} "
          f"{'UA fit':>8s} {'achieved':>9s}")
    for name, tmed, imm, tgt in TARGETS:
        try:
            ua = brentq(lambda u: cooling_rate(u, tmed, imm) - tgt, 0.0, 400.0,
                        xtol=1e-3)
        except ValueError:
            lo, hi = cooling_rate(0.0, tmed, imm), cooling_rate(400.0, tmed, imm)
            print(f"{name:32s} UNREACHABLE  passive={lo:.4f} max={hi:.4f} "
                  f"target={tgt}")
            continue
        got = cooling_rate(ua, tmed, imm)
        out[name] = (round(ua, 2), tmed, imm)
        print(f"{name:32s} {tmed:6.1f} {imm:4.0f} {tgt:7.3f} {ua:8.2f} {got:9.4f}")
    return out


if __name__ == "__main__":
    print("=== passive baseline (no device) ===")
    print(f"  passive core cooling in reference patient: "
          f"{cooling_rate(0.0, 30.0, 0.0):.4f} C/min")
    print()
    print("=== fitted lumped conductances ===")
    ua = fit()
    print()
    print("MODALITY_UA = {")
    for k, (u, tm, im) in ua.items():
        print(f'    "{k}": dict(UA_COOL={u}, T_COOL={tm}, IMMERSE={im:.0f}),')
    print("}")
