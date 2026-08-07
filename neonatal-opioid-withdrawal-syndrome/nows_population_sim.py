#!/usr/bin/env python3
"""
Virtual-population runs for the NOWS QSP model.

The deterministic scenarios in nows_verify_python.py describe ONE infant.  The
two trial-level numbers this model is asked to reproduce -- the proportion of
opioid-exposed infants who receive pharmacotherapy, and how that proportion
moves under Eat-Sleep-Console -- are population quantities, so they need a
population.

The point of this file is not to re-fit the trials.  It is the DECOMPOSITION
that the trials cannot do: ESC-NOW changed the assessment instrument AND the
care environment at the same time, so its 52.0% -> 19.5% result cannot be
attributed.  Here the two are separable, because the model applies them through
different parameters:

    ESCMODE  changes WHAT IS MEASURED    (functional failure, not a sum of signs)
    CARE     changes WHAT THE INFANT HAS (arousal drive on the locus coeruleus)
    BF       changes the OPIOID ON BOARD (2-9 ng/mL of methadone from milk)

Four cells therefore run: usual care, ESC criterion alone, care bundle alone,
and both together.

Between-subject variability is placed on the quantities that plausibly carry it
-- maternal dose, gestational age and weight, clearance of each drug, the gain
of the adaptation map, the gain from gap to sign expression, co-exposures, and
the Finnegan threshold itself (inter-rater variance makes the threshold a random
variable, not a constant).

The population median infant is deliberately MILDER than the deterministic base
case: the base case is a moderately severe methadone-exposed term infant, while
about half of an unselected exposed cohort never crosses the treatment
threshold at all.  The severity multiplier POPA carries that difference and is
the one number in this file fitted to a trial-level outcome rather than derived.

Run:  python3 nows_population_sim.py            (n = 200 per cell)
      POPA=0.52 POPS=0.35 python3 nows_population_sim.py

NOT FOR CLINICAL USE.
"""
import copy
import os
import multiprocessing as mp

import numpy as np

import nows_verify_python as M

POPA = float(os.environ.get("POPA", "0.52"))   # cohort severity multiplier
POPS = float(os.environ.get("POPS", "0.35"))   # between-subject SD (log scale)
NSUB = int(os.environ.get("NSUB", "200"))


def draw(rng, arm, care_mode):
    """One virtual infant."""
    p = copy.deepcopy(M.P)
    ga = float(np.clip(rng.normal(38.4, 1.7), 33.0, 41.5))
    p["GA"] = ga
    p["WT0"] = float(np.clip(3.30 * (ga / 39.0) ** 2.6 * np.exp(rng.normal(0, 0.12)),
                             1.4, 4.6))
    if arm == "meth":
        p["MDRUG"] = 1
        p["MDOSE"] = float(np.clip(85 * np.exp(rng.normal(0, 0.42)), 20, 220))
    else:
        p["MDRUG"] = 2
        p["MDOSE"] = float(np.clip(16 * np.exp(rng.normal(0, 0.42)), 2, 40))
    p["CLMREF"] *= float(np.exp(rng.normal(0, 0.30)))
    p["CLDREF"] *= float(np.exp(rng.normal(0, 0.30)))
    p["CLBREF"] *= float(np.exp(rng.normal(0, 0.30)))
    p["AGAIN"] *= float(POPA * np.exp(rng.normal(0, POPS)))
    p["KLC"] *= float(np.exp(rng.normal(0, 0.20)))
    p["FBZD"] = 1.0 if rng.random() < 0.22 else 0.0
    p["FNIC"] = float(rng.random() * 20) if rng.random() < 0.6 else 0.0

    if care_mode == "usual":                 # Finnegan tool, usual environment
        p["CARE"] = float(np.clip(rng.normal(0.50, 0.12), 0.1, 0.9))
        p["ESCMODE"] = 0
        p["BF"] = 1.0 if rng.random() < 0.25 else 0.0
        p["THRSTART"] = float(np.clip(rng.normal(8.0, 1.1), 5, 12))
    elif care_mode == "esc":                 # ESC criterion + full care bundle
        p["CARE"] = float(np.clip(rng.normal(0.86, 0.08), 0.3, 0.99))
        p["ESCMODE"] = 1
        p["BF"] = 1.0 if rng.random() < 0.65 else 0.0
    elif care_mode == "esc_criterion_only":  # ESC criterion, usual environment
        p["CARE"] = float(np.clip(rng.normal(0.50, 0.12), 0.1, 0.9))
        p["ESCMODE"] = 1
        p["BF"] = 1.0 if rng.random() < 0.25 else 0.0
    elif care_mode == "care_only":           # full care bundle, Finnegan tool
        p["CARE"] = float(np.clip(rng.normal(0.86, 0.08), 0.3, 0.99))
        p["ESCMODE"] = 0
        p["BF"] = 1.0 if rng.random() < 0.65 else 0.0
        p["THRSTART"] = float(np.clip(rng.normal(8.0, 1.1), 5, 12))
    return p


def sim(arm, care_mode, n=NSUB, seed=11, tmax=45 * 24):
    rng = np.random.default_rng(seed)
    treated, dur, cum, ready, burden = [], [], [], [], []
    for _ in range(n):
        p = draw(rng, arm, care_mode)
        try:
            o = M.run(p, tmax, derived=False)
            m = M.metrics(p, o)
        except Exception:
            continue
        treated.append(bool(m["treated"]))
        cum.append(m["cum"] * p["WT0"])          # total mg, as the trials report it
        burden.append(m["aucgap"])
        if m["treated"]:
            dur.append(m["days"])
        ready.append(m["ready"] if not np.isnan(m["ready"]) else p["MINOBS"] / 24)
    t = np.array(treated)
    return dict(n=int(len(t)), pct=float(100 * t.mean()),
                dur=float(np.median(dur)) if dur else 0.0,
                cum=float(np.mean(cum)),
                ready=float(np.median(ready)),
                burden=float(np.median(burden)))


CELLS = [("meth", "usual"), ("meth", "esc"), ("meth", "esc_criterion_only"),
         ("meth", "care_only"), ("bup", "usual"), ("bup", "esc")]

LABEL = {"usual": "Finnegan tool, usual care",
         "esc": "ESC criterion + care bundle",
         "esc_criterion_only": "ESC criterion ONLY",
         "care_only": "care bundle ONLY"}


def _job(a):
    return (a, sim(a[0], a[1]))


if __name__ == "__main__":
    print(f"NOWS virtual population — n={NSUB} per cell, POPA={POPA}, POPS={POPS}\n")
    print(f"{'arm':6s} {'cell':30s} {'treated':>8s} {'med days':>9s} "
          f"{'mean mg':>8s} {'med ready':>10s} {'med burden':>11s}")
    print("-" * 80)
    res = {}
    with mp.Pool(min(3, mp.cpu_count())) as po:
        for a, r in po.map(_job, CELLS):
            res[a] = r
            print(f"{a[0]:6s} {LABEL[a[1]]:30s} {r['pct']:7.1f}% {r['dur']:9.1f} "
                  f"{r['cum']:8.2f} {r['ready']:10.1f} {r['burden']:11.1f}")
    print("-" * 80)

    u = res[("meth", "usual")]["pct"]
    e = res[("meth", "esc")]["pct"]
    c = res[("meth", "esc_criterion_only")]["pct"]
    k = res[("meth", "care_only")]["pct"]
    print("\nDECOMPOSITION of the Eat-Sleep-Console effect (percentage points):")
    print(f"  usual care                                  {u:5.1f}%")
    print(f"  changing only the CRITERION                 {c:5.1f}%   "
          f"({u - c:+.1f} pp)")
    print(f"  changing only the CARE ENVIRONMENT          {k:5.1f}%   "
          f"({u - k:+.1f} pp)")
    print(f"  changing both (the trial's intervention)    {e:5.1f}%   "
          f"({u - e:+.1f} pp)")
    print(f"  sum of the single effects {(u - c) + (u - k):.1f} pp vs {u - e:.1f} pp "
          f"together — they overlap and cannot be added.")
    print("\nESC-NOW (Young 2023, NEJM) reported 52.0% -> 19.5% and a time to "
          "medically\nready for discharge of 14.9 -> 8.2 days.  The model is asked "
          "to reproduce the\nSHAPE of that change and the decomposition; the "
          "absolute percentages depend on\nthe cohort severity multiplier POPA, "
          "which is the one number here fitted to a\ntrial-level outcome.")
