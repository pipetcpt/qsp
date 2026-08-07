"""Calibration checks and physiological anchors for the thermal core."""
import numpy as np
from hs_core import (make_params, Schedule, simulate, get, derived, dubois,
                     wet_bulb, y0, IX)
from hs_calibrate import equilibrate_shell, reference_patient, REF_ENV, cooling_rate

np.set_printoptions(suppress=True)

print("=" * 74)
print("1. COMPENSABLE vs UNCOMPENSABLE -- the fixed point either exists or not")
print("=" * 74)
for lbl, env, mex, dur in [
        ("rest, 30C / 50% RH, light clothing", dict(TA=30., RH=.5, VAIR=.5, ICL=.5), 0., 240),
        ("walk 300 W, 35C / 40% RH",           dict(TA=35., RH=.4, VAIR=1., ICL=.3), 300., 240),
        ("run 900 W, 35C / 80% RH, sun",       dict(TA=35., RH=.8, VAIR=1.5, ICL=.15, QSOL=120.), 900., 90),
        ("run 900 W, 22C / 40% RH, breeze",    dict(TA=22., RH=.4, VAIR=3., ICL=.15), 900., 240)]:
    p = make_params()
    sch = Schedule([dict(dur=dur, MEX=mex)], env)
    T, Y, p = simulate(p, sch)
    D = derived(T, Y, p, sch)
    tc = get(Y, "TCR")
    plateau = abs(tc[-1] - tc[-31]) < 0.05
    print(f"  {lbl:36s} Tc {tc[0]:.2f}->{tc[-1]:5.2f} "
          f"HSI {D['HSI'][dur//2]:5.0f}%  {'PLATEAU' if plateau else 'RISING'}  "
          f"CEM43 {get(Y,'CEM43')[-1]:6.2f}")

print()
print("=" * 74)
print("2. RATE OF RISE -- the clock, from the same equation")
print("=" * 74)
env = dict(TA=35., RH=.8, VAIR=1.5, ICL=.15, QSOL=120.)
p = make_params()
sch = Schedule([dict(dur=90, MEX=900.)], env)
T, Y, p = simulate(p, sch)
D = derived(T, Y, p, sch)
tc = get(Y, "TCR")
i40, i42 = np.argmax(tc >= 40.), np.argmax(tc >= 42.)
print(f"  EHS runner:  t(40C)={T[i40]:5.1f} min  t(42C)={T[i42]:5.1f} min  "
      f"rate 40->42 = {(tc[i42]-tc[i40])/(T[i42]-T[i40]):.3f} C/min")
print(f"    Hprod {D['H_prod'][30]:.0f} W | Qdry {D['Q_dry'][30]:+.0f} W | "
      f"Emax {D['E_max'][30]:.0f} W | Eact {D['E_act'][30]:.0f} W | "
      f"sweat {D['m_sw'][30]*60/1000:.2f} L/h | HSI {D['HSI'][30]:.0f}%")
print(f"    muscle-core gradient at 30 min: {get(Y,'TMU')[30]-tc[30]:+.2f} C")

print()
print("=" * 74)
print("3. COOLING MODALITY -- fitted conductances, and one prediction")
print("=" * 74)
MODALITY = {
    "ice-water immersion 2C":   dict(UA_COOL=22.08, T_COOL=2.0,  IMMERSE=1),
    "cold-water immersion 14C": dict(UA_COOL=25.93, T_COOL=14.0, IMMERSE=1),
    "tarp-assisted 10C":        dict(UA_COOL=18.76, T_COOL=10.0, IMMERSE=1),
    "cold shower/dousing 20C":  dict(UA_COOL=14.43, T_COOL=20.0, IMMERSE=0),
    "evaporative + convective": dict(UA_COOL=14.58, T_COOL=25.0, IMMERSE=0),
    "endovascular catheter":    dict(UA_COOL=4.10,  T_COOL=4.0,  IMMERSE=0),
    "ice packs neck/axilla":    dict(UA_COOL=1.35,  T_COOL=0.0,  IMMERSE=0),
    "passive shade":            dict(UA_COOL=0.0,   T_COOL=30.0, IMMERSE=0),
}
for lbl, m in MODALITY.items():
    r = cooling_rate(m["UA_COOL"], m["T_COOL"], m["IMMERSE"])
    print(f"  {lbl:28s} UA={m['UA_COOL']:5.2f} W/K  ->  {r:.4f} C/min")
print("  PREDICTION: same tub (UA fitted at 2C = 22.08) filled with 14C water")
print(f"              -> {cooling_rate(22.08, 14.0, 1.0):.4f} C/min "
      f"(observed 0.17; ratio 2C:14C model "
      f"{cooling_rate(22.08,2.,1.)/cooling_rate(22.08,14.,1.):.2f}, observed ~1.16-1.29)")

print()
print("=" * 74)
print("4. COLD INTRAVENOUS SALINE -- a fixed, computable quantity")
print("=" * 74)
for vol, rate in [(1.0, 0.05), (2.0, 0.10), (2.0, 0.067)]:
    p, yi = reference_patient()
    dur = vol / rate
    sch = Schedule([dict(dur=dur, MEX=0., IVF_RATE=rate, IVF_TEMP=4.0),
                    dict(dur=60., MEX=0.)], REF_ENV)
    T, Y, p = simulate(p, sch, y_init=yi)
    tc = get(Y, "TCR")
    passive = 0.0156 * (dur + 60.)
    print(f"  {vol:.0f} L @4C over {dur:.0f} min: nadir {tc.min():.2f} C "
          f"(fall {41.0-tc.min():.2f}), at +60 min {tc[-1]:.2f}; "
          f"passive-only would give {41.0-passive:.2f}")
print("  analytic whole-body equivalent for 2 L: "
      f"{2*4.18*(41-4)/(70*3.470):.2f} C")

print()
print("=" * 74)
print("5. CLASSIC (NON-EXERTIONAL) HEAT STROKE -- the slow clock")
print("=" * 74)
for lbl, ta, rh, kw in [
        ("healthy adult, 38C/45%", 38., .45, {}),
        ("elderly, 40C/55%",       40., .55, dict(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05)),
        ("elderly + anticholinergic, 40C/55%", 40., .55,
         dict(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05, FSW_DRUG=.30)),
        ("elderly + anticholinergic, 43C/45%", 43., .45,
         dict(FSW_AGE=.60, FVD_AGE=.55, MREST=1.05, FSW_DRUG=.30))]:
    p = make_params(**kw)
    env = dict(TA=ta, RH=rh, VAIR=0.15, ICL=0.7)
    sch = Schedule([dict(dur=4320, MEX=0.)], env)
    T, Y, p = simulate(p, sch, max_step=10.0)
    tc = get(Y, "TCR")
    h = lambda n: tc[min(int(n * 60), len(tc) - 1)]
    print(f"  {lbl:38s} 6h {h(6):5.2f}  12h {h(12):5.2f}  24h {h(24):5.2f}  "
          f"48h {h(48):5.2f}  72h {h(72):5.2f}  CEM43 {get(Y,'CEM43')[-1]:7.1f}")

print()
print("=" * 74)
print("6. PSYCHROMETRIC SANITY")
print("=" * 74)
for ta, rh in [(35, 1.0), (40, .5), (45, .3), (35, .5), (31, .8)]:
    print(f"  wet-bulb({ta}C, {rh*100:.0f}%) = {wet_bulb(ta, rh):5.1f} C")
