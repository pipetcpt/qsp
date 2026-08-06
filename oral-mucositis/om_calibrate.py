#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
om_calibrate.py -- fit the small number of free parameters of the oral
mucositis model, and record exactly which published numbers were spent.

WHAT IS FITTED (stage 1, four parameters, four numbers)
-------------------------------------------------------
  pot_mel   potency of mucosal melphalan on the clonogen pool
  rad_pot   potency of a Gy of mucosal dose on the same pool
  lamS      basal clonogen renewal rate  -- the REGENERATION clock
  greg      gain of the barrier-deficit -> proliferation loop

against the TIMING of severe (WHO >= 3) mucositis in two regimens whose
insults could hardly be more different -- one impulse, one seven-week ramp:

  HDM 200 mg/m2 autologous HSCT   onset  5.5 d, duration  6.0 d
  H&N chemoRT 70 Gy/35 fx + cis   onset 18.0 d, duration 35.0 d

Four numbers, four parameters: the fit is exactly determined, so agreement
here is NOT evidence.  Everything in om_analysis.py beyond section 1 is.

WHAT IS FITTED (stage 2, one parameter, one number)
---------------------------------------------------
  cy_equiv  the melphalan-equivalent alkylator exposure that stands in for
            the etoposide + cyclophosphamide block of the Spielberger 2004
            TBI conditioning regimen, fitted to the PLACEBO arm's median
            grade 3/4 duration of 9 days.  The palifermin arm is then a
            PREDICTION with no further fitting.

WHAT IS FITTED (stage 4, one parameter, one number)
---------------------------------------------------
  pot_5fu   fitted so that bolus 5-FU 425 mg/m2 d1-5 produces 5 days of
            ULCERATIVE (WHO >= 2) stomatitis in the typical patient -- the
            classic Mayo-regimen picture.  5-FU is the only cycle-active
            agent used in the shipped scenarios, so without this number the
            whole S-phase-specific arm of the model (and with it the
            palifermin paradox in a 5-FU setting) would be unconstrained.

WHAT IS FITTED (stage 3, one parameter, one number)
---------------------------------------------------
  pot_cis   fitted to the incremental severe-OM incidence of concurrent
            cisplatin over radiotherapy alone in the virtual population
            (Trotti 2003: 34% RT alone -> 43% chemoradiation).

pot_mtx is NOT constrained by anything here -- methotrexate is carried in
the model for completeness (GVHD prophylaxis in allogeneic transplant) but
appears in none of the shipped scenarios, so its potency is a placeholder and
is labelled as such.

NOTHING ELSE IS FITTED.  Cryotherapy, photobiomodulation, benzydamine,
fractionation, treatment gaps and the palifermin scheduling window are all
run forward from these six numbers.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import least_squares

import om_python_reference as M

SCHEDS = dict(
    crt=lambda **kw: M.sched_chemoRT(**kw),
    hdm=lambda **kw: M.sched_HDM(**kw),
)

LOG = open("calib.log", "w")


def say(s=""):
    # opened lazily: om_analysis.py imports this module for its schedule
    # builders, and an eager open(..., "w") would truncate the calibration
    # log every time the analysis runs.
    global LOG
    print(s)
    if LOG is None:
        LOG = open("calib.log", "a")
    LOG.write(s + "\n")
    LOG.flush()


# ----------------------------------------------------------------------------
# targets
# ----------------------------------------------------------------------------
TARGETS = dict(
    hdm_onset=6.5,     # d after melphalan, WHO>=3
    hdm_dur=6.0,       # d
    crt_onset=19.0,    # d after first fraction, WHO>=3
    crt_dur=38.0,      # d
)
SPIELBERGER_PLACEBO_DUR = 9.0     # d, median grade 3/4 (WHO) duration
TROTTI_RT = 0.34                  # severe OM, RT alone
TROTTI_CRT = 0.43                 # severe OM, concurrent chemoradiation


def run_hdm(p, t_end=45.0, n=901):
    return M.metrics(M.simulate(M.sched_HDM(200.0), p=p, t_end=t_end, n=n))


def run_crt(p, t_end=105.0, n=1051, **kw):
    s = M.sched_chemoRT(**kw)
    return M.metrics(M.simulate(s, p=p, t_end=t_end, n=n))


# ----------------------------------------------------------------------------
# stage 1
# ----------------------------------------------------------------------------
FIT1 = ["pot_mel", "rad_pot", "lamS", "greg"]


def _bisect(f, lo, hi, target, iters=15, increasing=True):
    """monotone 1-D solve"""
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        v = f(mid)
        if (v < target) == increasing:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _nan(v, fb):
    return fb if v != v else v


def solve_inner(lamS, base=None):
    """
    THE NESTING IS DICTATED BY THE SENSITIVITIES, NOT BY TASTE.

    A sweep of the four candidate parameters (recorded in this log) showed:

      * for an IMPULSE insult, severe-OM ONSET is almost potency-independent
        (6.70 -> 6.45 d over a 33% change in pot_mel) but strongly dependent
        on the regeneration GAIN greg (3.90 -> 7.40 d as greg goes 4 -> 14);
      * DURATION carries essentially all of the dose-response
        (5.55 -> 8.60 d over that same 33% potency change).

    That is itself the two-clock claim showing up in the derivatives, and it
    means the well-conditioned assignment is

        greg     <- HDM onset          pot_mel  <- HDM duration
        rad_pot  <- chemoRT onset      lamS     <- chemoRT duration

    An earlier design pinned the potencies with the ONSETS.  Because onset is
    nearly flat in potency there, the bisection ran to potencies 6x higher
    than needed, annihilated the clonogen pool, and produced a 39-day
    non-resolving ulcer.  Getting the nesting wrong was what exposed the
    missing quiescent reserve.
    """
    p = dict(base if base is not None else M.P)
    p["lamS"] = lamS

    def onset_hdm(g):
        q = dict(p); q["greg"] = g
        return _nan(run_hdm(q, t_end=16.0, n=321)["onset_sev"], 0.0)

    def dur_hdm(pot):
        q = dict(p); q["pot_mel"] = pot
        return _nan(run_hdm(q, t_end=50.0, n=1001)["dur_sev"], 0.0)

    # greg and pot_mel are mildly coupled; two sweeps converge them
    for _ in range(2):
        p["greg"] = _bisect(onset_hdm, 1.0, 40.0, TARGETS["hdm_onset"],
                            increasing=True)
        p["pot_mel"] = _bisect(dur_hdm, 5.0, 600.0, TARGETS["hdm_dur"],
                               increasing=True)

    def onset_crt(rp):
        q = dict(p); q["rad_pot"] = rp
        return _nan(run_crt(q, t_end=34.0, n=341)["onset_sev"], 999.0)

    p["rad_pot"] = _bisect(onset_crt, 0.05, 20.0, TARGETS["crt_onset"],
                           increasing=False)
    return p


def stage1():
    say("=" * 78)
    say("STAGE 1 -- four parameters against four timing numbers")
    say("  greg    <- HDM onset       pot_mel <- HDM duration")
    say("  rad_pot <- chemoRT onset   lamS    <- chemoRT duration")
    say("  (assignment chosen from the measured sensitivities, see")
    say("   solve_inner's docstring)")
    say("=" * 78)
    say("")
    say("  %-9s %9s %9s %9s %9s %9s"
        % ("lamS", "greg", "pot_mel", "rad_pot", "crt_dur", "resid"))
    grid = [0.06, 0.09, 0.12, 0.16, 0.21, 0.28, 0.36]
    recs = []
    for lam in grid:
        try:
            p = solve_inner(lam)
            c = run_crt(p)
            cd = _nan(c["dur_sev"], 0.0)
        except Exception as e:
            say("  %-9.3f  FAILED %r" % (lam, e))
            continue
        recs.append((lam, p, cd))
        say("  %-9.3f %9.2f %9.1f %9.3f %9.2f %9.2f"
            % (lam, p["greg"], p["pot_mel"], p["rad_pot"], cd,
               cd - TARGETS["crt_dur"]))

    if not recs:
        raise SystemExit("stage 1 produced no usable points")
    best = min(recs, key=lambda r: abs(r[2] - TARGETS["crt_dur"]))
    # local refinement around the best grid point
    lo = max(best[0] * 0.65, 0.03)
    hi = min(best[0] * 1.55, 0.6)
    say("")
    say("  refining lamS in [%.3f, %.3f]" % (lo, hi))
    for _ in range(6):
        mid = 0.5 * (lo + hi)
        p = solve_inner(mid)
        cd = _nan(run_crt(p)["dur_sev"], 0.0)
        say("    lamS %.4f -> chemoRT duration %.2f" % (mid, cd))
        if cd > TARGETS["crt_dur"]:
            lo = mid          # faster regeneration shortens it
        else:
            hi = mid
    lamS = 0.5 * (lo + hi)
    p = solve_inner(lamS)

    say("")
    for k in FIT1:
        say("    %-10s %.5g" % (k, p[k]))
    h = run_hdm(p, t_end=50.0, n=1001)
    c = run_crt(p)
    say("")
    say("  %-34s %8s %8s" % ("target", "obs", "model"))
    say("  %-34s %8.1f %8.2f" % ("HDM onset WHO>=3 (d)",
                                 TARGETS["hdm_onset"], h["onset_sev"]))
    say("  %-34s %8.1f %8.2f" % ("HDM duration WHO>=3 (d)",
                                 TARGETS["hdm_dur"], h["dur_sev"]))
    say("  %-34s %8.1f %8.2f" % ("chemoRT onset WHO>=3 (d)",
                                 TARGETS["crt_onset"], c["onset_sev"]))
    say("  %-34s %8.1f %8.2f" % ("chemoRT duration WHO>=3 (d)",
                                 TARGETS["crt_dur"], c["dur_sev"]))
    say("")
    say("  UNFITTED consequences of the same four numbers:")
    say("    HDM peak WHO grade            %.0f   (observed: 3-4)"
        % h["peak_who"])
    say("    HDM peak ulcer area fraction  %.3f" % h["peak_area"])
    say("    chemoRT peak WHO grade        %.0f   (observed: 3-4)"
        % c["peak_who"])
    say("    chemoRT severe OM ends        %.1f d after the last fraction"
        % (c["end_sev"] - 46.0))
    say("    HDM ANC nadir                 %.2f x10^9/L" % h["nadir_anc"])
    say("    HDM days on opioid            %.1f d" % h["opidays"])
    return p


# ----------------------------------------------------------------------------
# stage 2 -- the Spielberger conditioning block
# ----------------------------------------------------------------------------
def sched_spielberger(cy_equiv, bsa=1.8, pal_days=None, t0=3.0):
    """
    Spielberger 2004 (NEJM 351:2590) conditioning: fractionated TBI 12 Gy in
    8 fractions over 4 days (days -8..-5), then etoposide 60 mg/kg and
    cyclophosphamide 100 mg/kg (days -4, -2), stem cells day 0.

    The etoposide/cyclophosphamide block is represented by a single
    melphalan-EQUIVALENT alkylator exposure `cy_equiv` (mg/m2), which is the
    ONE number stage 2 fits.  TBI is modelled explicitly with the same
    alpha/beta and the same rad_pot that stage 1 already fixed.
    """
    # t0 shifts the whole conditioning block later so that the three
    # PRE-conditioning palifermin doses (label: 3 consecutive days before
    # conditioning begins) have somewhere to live on a timeline that starts
    # at zero.  It is a pure translation: with no palifermin the trajectory
    # is identical, which is why stage 2's fit is unaffected by it.
    s = {k: (list(v) if isinstance(v, list) else v)
         for k, v in M.ZERO_SCHED.items()}
    rt = []
    # conditioning days 1-4: two 1.5 Gy fractions per day
    for d in [0.0, 1.0, 2.0, 3.0]:
        rt.append((t0 + d, 1.5, 10.0 / 1440.0))
        rt.append((t0 + d + 0.35, 1.5, 10.0 / 1440.0))
    s["rt"] = rt
    s["dose_per_fx"] = 1.5
    amt = cy_equiv * bsa * 1000.0
    dur = 0.5 / 24.0
    s["mel"] = [(t0 + 4.0, t0 + 4.0 + dur, 0.5 * amt / dur),
                (t0 + 6.0, t0 + 6.0 + dur, 0.5 * amt / dur)]
    if pal_days:
        s = M.add_palifermin(s, pal_days)
    return s


def stage2(p):
    say("")
    say("=" * 78)
    say("STAGE 2 -- one parameter (alkylator equivalent of the VP16/Cy block)")
    say("           against ONE number: Spielberger placebo duration = 9 d")
    say("=" * 78)

    def f(ce):
        m = M.metrics(M.simulate(sched_spielberger(ce), p=p,
                                 t_end=45.0, n=1801))
        return (0.0 if m["dur_sev"] != m["dur_sev"] else m["dur_sev"])

    lo, hi = 20.0, 400.0
    for _ in range(38):
        mid = 0.5 * (lo + hi)
        if f(mid) < SPIELBERGER_PLACEBO_DUR:
            lo = mid
        else:
            hi = mid
    ce = 0.5 * (lo + hi)
    say("  cy_equiv = %.1f mg/m2 melphalan-equivalent" % ce)
    say("  placebo arm duration WHO>=3: %.2f d (target 9.0)" % f(ce))
    return ce


# ----------------------------------------------------------------------------
# virtual population
# ----------------------------------------------------------------------------
def vpop(n=400, seed=20260806, cv_sens=0.38, cv_lam=0.30, cv_dcrit=0.07):
    rng = np.random.default_rng(seed)
    return dict(
        sens=np.exp(rng.normal(0, cv_sens, n) - 0.5 * cv_sens ** 2),
        lam=np.exp(rng.normal(0, cv_lam, n) - 0.5 * cv_lam ** 2),
        dcrit=np.exp(rng.normal(0, cv_dcrit, n) - 0.5 * cv_dcrit ** 2),
    )


_POOL = None


def _pool():
    global _POOL
    if _POOL is None:
        _POOL = Pool(processes=max(1, (os.cpu_count() or 2)))
    return _POOL


def _one(arg):
    kind, kw, p, sev = arg
    try:
        kw = dict(kw)
        # pop the horizon BEFORE building the schedule -- it is not a
        # schedule argument.  Leaving it in raised TypeError inside every
        # worker, every worker returned None, and the population silently
        # reported an incidence of 0.000 rather than failing loudly.
        t_end = kw.pop("_t_end", 60.0)
        s = SCHEDS[kind](**kw)
        r = M.simulate(s, p=p, t_end=t_end, n=int(t_end * 10) + 1)
        m = M.metrics(r, sev=sev)
        return (m["incidence"], m["dur_sev"], m["peak_who"], m["peak_area"],
                m["painAUC"], m["opidays"])
    except Exception:
        return None


def run_pop(kind, kw, p, pop, t_end=60.0, sev=3.0):
    """
    Virtual population.  Inter-individual variability is put on exactly three
    handles -- overall cytotoxic sensitivity, regenerative capacity, and the
    barrier threshold -- because those are the three places the structure says
    a patient can differ.  Everything else is shared.
    """
    N = len(pop["sens"])
    args = []
    for i in range(N):
        pi = dict(p)
        pi["sens"] = p["sens"] * float(pop["sens"][i])
        pi["lamS"] = p["lamS"] * float(pop["lam"][i])
        pi["Dcrit"] = p["Dcrit"] * float(pop["dcrit"][i])
        kwi = dict(kw)
        kwi["_t_end"] = t_end
        args.append((kind, kwi, pi, sev))
    out = [o for o in _pool().map(_one, args, chunksize=4) if o is not None]
    n_ok = len(out)
    inc = sum(1 for o in out if o[0]) / max(n_ok, 1)
    g4 = sum(1 for o in out if o[2] >= 4.0) / max(n_ok, 1)
    durs = np.array([o[1] if o[0] else 0.0 for o in out])
    pos = durs[durs > 0]
    return dict(inc=inc, inc4=g4,
                med_dur=float(np.median(pos)) if len(pos) else 0.0,
                mean_dur=float(durs.mean()), durs=durs, n=n_ok,
                mean_painAUC=float(np.mean([o[4] for o in out])),
                mean_opidays=float(np.mean([o[5] for o in out])))


def stage3(p, pop, n_keep=None):
    """
    TWO population parameters against TWO population numbers.

    A subtlety that the first run of this stage got wrong.  Stage 1 calibrated
    the DETERMINISTIC patient to a severe-mucositis DURATION of 38 days -- a
    number that is only defined for patients who actually get severe
    mucositis.  So the stage-1 patient is the median AFFECTED patient, not the
    median enrolled patient, and a population centred on that patient has an
    incidence near 100%, not the 34-43% the trials report.

    The fix is not to re-fit stage 1; it is to admit that the population has a
    LOCATION as well as a spread.  sens_med is that location, and it is fitted
    to the radiotherapy-alone incidence.  pot_cis is then fitted to the
    increment that concurrent cisplatin adds.
    """
    say("")
    say("=" * 78)
    say("STAGE 3 -- two population parameters against two population numbers:")
    say("           Trotti 2003 severe OM, %.0f%% RT alone -> %.0f%% chemoRT"
        % (100 * TROTTI_RT, 100 * TROTTI_CRT))
    say("=" * 78)
    say("")
    say("  (a) sens_med -- the population's LOCATION on the sensitivity axis")
    say("      (the stage-1 patient is the median AFFECTED patient, so the")
    say("       population must sit BELOW it)")

    def inc_rt(sm):
        pp = dict(p)
        pp["sens"] = sm
        return run_pop("crt", dict(cis_mgm2=0.0), pp, pop, t_end=95.0)["inc"]

    lo, hi = 0.20, 1.05
    for _ in range(7):
        mid = 0.5 * (lo + hi)
        v = inc_rt(mid)
        say("      sens_med %.4f -> RT-alone incidence %.3f" % (mid, v))
        if v < TROTTI_RT:
            lo = mid
        else:
            hi = mid
    sens_med = 0.5 * (lo + hi)
    rt_inc = inc_rt(sens_med)
    say("      sens_med = %.4f  (RT-alone incidence %.3f, target %.2f)"
        % (sens_med, rt_inc, TROTTI_RT))

    say("")
    say("  (b) pot_cis -- the increment concurrent cisplatin adds")
    pbase = dict(p)
    pbase["sens"] = sens_med

    def inc_crt(pc):
        pp = dict(pbase)
        pp["pot_cis"] = pc
        return run_pop("crt", dict(cis_mgm2=100.0), pp, pop,
                       t_end=95.0)["inc"]

    grid = [0.0, 400.0, 1200.0, 3000.0]
    incs = []
    for pc in grid:
        v = inc_crt(pc)
        incs.append(v)
        say("      pot_cis %8.1f -> chemoRT incidence %.3f" % (pc, v))
    incs = np.array(incs)
    if incs.max() < TROTTI_CRT:
        pc = grid[-1]
        say("      WARNING: target not reachable on this grid; taking %.1f"
            % pc)
    else:
        pc = float(np.interp(TROTTI_CRT, incs, grid))
    say("      pot_cis = %.4g" % pc)
    return pc, rt_inc, sens_med


# ----------------------------------------------------------------------------
def main():
    p = stage1()
    ce = stage2(p)
    pop = vpop(n=int(sys.argv[1]) if len(sys.argv) > 1 else 120)
    pc, rt_inc = stage3(p, pop)
    p["pot_cis"] = pc
    out = dict(fitted={k: p[k] for k in FIT1},
               cy_equiv=ce, pot_cis=pc, rt_alone_incidence=rt_inc,
               targets=TARGETS,
               spielberger_placebo_dur=SPIELBERGER_PLACEBO_DUR,
               trotti=[TROTTI_RT, TROTTI_CRT])
    json.dump(out, open("om_calibration.json", "w"), indent=2)
    say("")
    say("wrote om_calibration.json")


if __name__ == "__main__":
    main()


# ----------------------------------------------------------------------------
# stage 4 -- the cycle-active agent
# ----------------------------------------------------------------------------
FU_ULCER_DAYS = 5.0     # d of WHO >= 2 stomatitis, bolus 5-FU d1-5


def stage4(p):
    say("")
    say("=" * 78)
    say("STAGE 4 -- one parameter (pot_5fu) against ONE number:")
    say("           bolus 5-FU 425 mg/m2 d1-5 -> %.0f d of WHO>=2 stomatitis"
        % FU_ULCER_DAYS)
    say("=" * 78)

    def f(pot):
        q = dict(p)
        q["pot_5fu"] = pot
        m = M.metrics(M.simulate(M.sched_5FU_bolus(), p=q, t_end=40.0, n=801),
                      sev=2.0)
        return _nan(m["dur_sev"], 0.0)

    pot = _bisect(f, 1.0, 40000.0, FU_ULCER_DAYS, iters=18, increasing=True)
    q = dict(p)
    q["pot_5fu"] = pot
    m2 = M.metrics(M.simulate(M.sched_5FU_bolus(), p=q, t_end=40.0, n=801),
                   sev=2.0)
    m3 = M.metrics(M.simulate(M.sched_5FU_bolus(), p=q, t_end=40.0, n=801),
                   sev=3.0)
    say("  pot_5fu = %.5g" % pot)
    say("  WHO>=2 duration %.2f d (target %.1f)" % (m2["dur_sev"],
                                                    FU_ULCER_DAYS))
    say("  UNFITTED: WHO>=3 duration %.2f d, peak WHO %d, peak area %.3f"
        % (m3["dur_sev"], int(m3["peak_who"]), m3["peak_area"]))
    return pot
