#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
om_analysis.py -- everything the oral mucositis model is actually FOR.

Section 1 replays the calibration.  Sections 2 onwards are consequences: no
parameter is touched after section 1 except where a section explicitly says it
is turning a mechanism OFF in order to decompose it.

Run:  python3 om_analysis.py        (writes om_reference_output.txt)
"""
import json
import os
import sys

import numpy as np

import om_python_reference as M

OUT = open("om_reference_output.txt", "w")


def say(s=""):
    print(s)
    OUT.write(s + "\n")
    OUT.flush()


def hdr(t):
    say("")
    say("=" * 78)
    say(t)
    say("=" * 78)


# ----------------------------------------------------------------------------
# load calibration
# ----------------------------------------------------------------------------
CAL = json.load(open("om_calibration.json"))
P = dict(M.P)
P.update(CAL["fitted"])
P["pot_cis"] = CAL["pot_cis"]
P["pot_5fu"] = CAL.get("pot_5fu", M.P["pot_5fu"])
CY_EQ = CAL["cy_equiv"]


def run(sched, t_end=60.0, n=None, p=None):
    n = n or int(t_end * 10) + 1
    return M.simulate(sched, p=p or P, t_end=t_end, n=n)


def mx(sched, t_end=60.0, p=None, sev=3.0):
    return M.metrics(run(sched, t_end=t_end, p=p), sev=sev)


# ============================================================================
def s1_calibration():
    hdr("1.  WHAT WAS FITTED, AND WHAT IT COST")
    say("")
    say("  EIGHT published numbers were spent, on eight parameters:")
    say("")
    say("    stage 1  pot_mel  rad_pot  lamS  greg   <- 4 timing numbers")
    say("             HDM onset %.1f d / duration %.1f d"
        % (CAL["targets"]["hdm_onset"], CAL["targets"]["hdm_dur"]))
    say("             chemoRT onset %.1f d / duration %.1f d"
        % (CAL["targets"]["crt_onset"], CAL["targets"]["crt_dur"]))
    say("    stage 2  cy_equiv                        <- Spielberger placebo"
        " duration %.0f d" % CAL["spielberger_placebo_dur"])
    say("    stage 3  sens_med                        <- Trotti RT-alone"
        " incidence %.2f" % CAL["trotti"][0])
    say("    stage 3  pot_cis                         <- Trotti chemoRT"
        " incidence %.2f" % CAL["trotti"][1])
    say("    stage 4  pot_5fu                         <- bolus 5-FU"
        " ulcerative stomatitis 5 d")
    say("")
    for k, v in CAL["fitted"].items():
        say("    %-10s = %.5g" % (k, v))
    say("    %-10s = %.5g mg/m2 melphalan-equivalent" % ("cy_equiv", CY_EQ))
    say("    %-10s = %.5g" % ("pot_cis", CAL["pot_cis"]))
    say("    %-10s = %.5g" % ("pot_5fu", P["pot_5fu"]))
    say("    %-10s = %.5g   (population location, not the typical patient)"
        % ("sens_med", CAL.get("sens_med", 1.0)))
    say("")
    say("  pot_mtx is NOT constrained by any of this.  Methotrexate is")
    say("  carried for completeness (GVHD prophylaxis) but appears in none")
    say("  of the shipped scenarios; its potency is a placeholder.")
    say("")
    say("  The fit is exactly determined: eight numbers, eight parameters.")
    say("  Agreement in this section is therefore NOT evidence -- it is")
    say("  arithmetic.  Sections 2-12 are the evidence, because nothing in")
    say("  them was fitted.  Where the fit does NOT close (the HDM onset")
    say("  below) that is reported rather than absorbed.")

    h = mx(M.sched_HDM(200.0), t_end=45.0)
    c = mx(M.sched_chemoRT(), t_end=105.0)
    say("")
    say("  Consequences that were NOT fitted:")
    say("    %-42s %s" % ("HDM peak WHO grade", int(h["peak_who"])))
    say("    %-42s %.3f" % ("HDM peak ulcer area fraction", h["peak_area"]))
    say("    %-42s %.2f" % ("HDM ANC nadir (x10^9/L)", h["nadir_anc"]))
    say("    %-42s %.1f d" % ("HDM days on opioid", h["opidays"]))
    say("    %-42s %s" % ("chemoRT peak WHO grade", int(c["peak_who"])))
    say("    %-42s %.1f d"
        % ("chemoRT severe OM ending after RT end", c["end_sev"] - 46.0))
    return h, c


# ============================================================================
def s2_two_clocks():
    hdr("2.  THE TWO CLOCKS -- the central claim, tested as a 2x2")
    say("")
    say("  CLAIM: onset belongs to the INSULT term, duration to the")
    say("  REGENERATION term.  If true, perturbing one term must move its own")
    say("  endpoint and leave the other's largely alone.")
    say("")
    say("  The insult term is perturbed by scaling `sens` (which multiplies")
    say("  every kill term and nothing else).  The regeneration term is")
    say("  perturbed by scaling `lamS` (which enters only dS/dt renewal).")
    say("")

    for name, sched, t_end, ref_on, ref_du in [
            ("HDM 200 mg/m2", M.sched_HDM(200.0), 45.0, None, None),
            ("chemoRT 70 Gy + cis", M.sched_chemoRT(), 105.0, None, None)]:
        base = mx(sched, t_end=t_end)
        say("  %s" % name)
        say("    %-26s %10s %10s" % ("perturbation", "onset d", "durn d"))
        say("    %-26s %10.2f %10.2f"
            % ("baseline", base["onset_sev"], base["dur_sev"]))
        rows = []
        for lab, key, fac in [("insult  x0.70", "sens", 0.70),
                              ("insult  x1.40", "sens", 1.40),
                              ("regen   x0.70", "lamS", 0.70),
                              ("regen   x1.40", "lamS", 1.40)]:
            p = dict(P)
            p[key] = P.get(key, 1.0) * fac
            m = mx(sched, t_end=t_end, p=p)
            say("    %-26s %10.2f %10.2f"
                % (lab, m["onset_sev"], m["dur_sev"]))
            rows.append((lab, m["onset_sev"], m["dur_sev"]))
        # elasticities
        d_on_ins = abs(rows[1][1] - rows[0][1]) / base["onset_sev"]
        d_du_ins = abs(rows[1][2] - rows[0][2]) / max(base["dur_sev"], 1e-9)
        d_on_reg = abs(rows[3][1] - rows[2][1]) / base["onset_sev"]
        d_du_reg = abs(rows[3][2] - rows[2][2]) / max(base["dur_sev"], 1e-9)
        say("    relative swing over a 2-fold change in the term:")
        say("      insult      -> onset %5.1f%%   duration %5.1f%%"
            % (100 * d_on_ins, 100 * d_du_ins))
        say("      regeneration-> onset %5.1f%%   duration %5.1f%%"
            % (100 * d_on_reg, 100 * d_du_reg))
        say("      SEPARATION RATIO (onset: insult/regen) = %.2f"
            % (d_on_ins / max(d_on_reg, 1e-9)))
        say("      SEPARATION RATIO (durn : regen/insult) = %.2f"
            % (d_du_reg / max(d_du_ins, 1e-9)))
        say("")
    say("  VERDICT -- AND THE CLAIM AS ORIGINALLY WRITTEN DOES NOT SURVIVE")
    say("  INTACT.  The ONSET half holds: in chemoRT a two-fold change in the")
    say("  insult term moves onset ~2.6x more than the same change in the")
    say("  regeneration term.  The DURATION half does NOT: the two terms move")
    say("  duration almost equally (separation ratio ~1.0-1.1).  That is not")
    say("  a numerical accident, it is obvious in hindsight -- duration is the")
    say("  time taken to climb back from a DEPTH, and the insult sets the")
    say("  depth.  The defensible statement is therefore the weaker one:")
    say("")
    say("      ONSET is owned by the insult term.")
    say("      DURATION is shared, and only the regeneration term can")
    say("      shorten it WITHOUT also reducing the insult.")
    say("")
    say("  Section 8 is where that weaker statement earns its keep: it is")
    say("  what separates an agent that can be given to a patient whose")
    say("  cytotoxic dose is fixed from one that cannot.")
    say("")


# ============================================================================
def s3_latency():
    hdr("3.  THE LATENT PERIOD IS 1/k_shed, NOT THE TRANSIT TIME")
    say("")
    say("  For an IMPULSE insult (HDM) the barrier is left with no incoming")
    say("  flux, so it decays at k_shed and ulcerates when it crosses D_crit.")
    say("  The prediction is therefore  t_onset ~ ln(1/D_crit)/k_shed,")
    say("  and it should NOT scale with the transit parameter k_p.")
    say("")
    say("  %-10s %10s %10s %10s" % ("k_shed", "onset d", "pred d", "ratio"))
    for ks in [0.20, 0.25, 0.30, 0.40, 0.50]:
        p = dict(P)
        p["k_shed"] = ks
        m = mx(M.sched_HDM(200.0), t_end=45.0, p=p)
        pred = np.log(1.0 / P["Dcrit"]) / ks
        say("  %-10.2f %10.2f %10.2f %10.2f"
            % (ks, m["onset_sev"], pred, m["onset_sev"] / pred))
    say("")
    say("  %-10s %10s   (transit parameter -- should move onset much less)"
        % ("k_p", "onset d"))
    for kp in [0.45, 0.62, 0.85, 1.10]:
        p = dict(P)
        p["k_p"] = kp
        m = mx(M.sched_HDM(200.0), t_end=45.0, p=p)
        say("  %-10.2f %10.2f" % (kp, m["onset_sev"]))
    say("")
    say("  THE ANALYTIC PREDICTION IS WRONG, AND THE TABLE SHOWS BY HOW")
    say("  MUCH.  ln(1/D_crit)/k_shed under-predicts the onset by 1.6x at")
    say("  k_shed 0.20 and by 3.2x at 0.50, and the model's onset moves only")
    say("  ~22% while the prediction moves 2.5-fold.  The reason is that the")
    say("  closed-form ignores the transit-amplifying pool, which keeps")
    say("  feeding the barrier for several days after the insult and is NOT")
    say("  emptied instantaneously.  What survives is the QUALITATIVE claim,")
    say("  and it survives cleanly: onset tracks the BARRIER LIFETIME")
    say("  parameter and is nearly indifferent to the TRANSIT parameter")
    say("  (k_p 0.45 -> 1.10, a 2.4-fold change, moves onset by 0.3 d).")
    say("")
    say("  For chemoRT the same two sweeps should behave DIFFERENTLY, because")
    say("  there the onset is set by how long the ramp takes to exhaust the")
    say("  clonogen pool, not by barrier lifetime:")
    say("  %-10s %10s" % ("k_shed", "onset d"))
    for ks in [0.20, 0.30, 0.50]:
        p = dict(P)
        p["k_shed"] = ks
        m = mx(M.sched_chemoRT(), t_end=105.0, p=p)
        say("  %-10.2f %10.2f" % (ks, m["onset_sev"]))


# ============================================================================
def s4_cryo():
    hdr("4.  THE CRYOTHERAPY CRITERION -- derived, not looked up")
    say("")
    say("  Ice acts only while it is in the mouth.  So the benefit must be")
    say("")
    say("      benefit  =  (1 - f_cryo) x (fraction of the insult AUC that")
    say("                                  falls inside the ice window)")
    say("")
    say("  and for a drug eliminated with half-life t_half, iced for T from")
    say("  the moment of dosing, that fraction is exactly  1 - 2^(-T/t_half).")
    say("  Everything below is the model computing that, with NO parameter")
    say("  fitted to any cryotherapy trial.")
    say("")

    def auc_mel(sch):
        r = run(sch, t_end=20.0)
        return float(r["AUCmuc"][-1])

    say("  (a) melphalan 200 mg/m2, t_half 1.2 h -- mucosal AUC vs ice window")
    say("      %-12s %12s %12s %12s"
        % ("T_ice (h)", "AUC ratio", "1-2^(-T/th)", "sev durn d"))
    base_auc = auc_mel(M.sched_HDM(200.0))
    base_m = mx(M.sched_HDM(200.0), t_end=45.0)
    say("      %-12s %12.4f %12s %12.2f" % ("none", 1.0, "-",
                                            base_m["dur_sev"]))
    th = 1.2 / 24.0
    for Th in [0.25, 0.5, 1.0, 2.0, 4.0, 6.0]:
        T = Th / 24.0
        s = M.add_cryo(M.sched_HDM(200.0), [(-0.02, T)])
        a = auc_mel(s)
        m = mx(s, t_end=45.0)
        say("      %-12.2f %12.4f %12.4f %12.2f"
            % (Th, a / base_auc, 1 - 2 ** (-T / th), m["dur_sev"]))

    say("")
    say("  (b) DECOMPOSITION: is cryotherapy a FLOW effect or a TEMPERATURE")
    say("      effect?  Run 6 h of ice with each arm disabled in turn.")
    s6 = M.add_cryo(M.sched_HDM(200.0), [(-0.02, 6.0 / 24.0)])
    variants = [
        ("both arms", dict()),
        ("flow arm only  (Q10 = 1)", dict(Q10_cryo=1.0)),
        ("temp arm only  (fQ  = 1)", dict(fQ_cryo=1.0)),
    ]
    say("      %-28s %12s %12s" % ("", "mucosal AUC", "sev durn d"))
    say("      %-28s %12.4f %12.2f" % ("no ice", base_auc, base_m["dur_sev"]))
    for lab, over in variants:
        p = dict(P)
        p.update(over)
        r = M.simulate(s6, p=p, t_end=20.0, n=201)
        m = M.metrics(M.simulate(s6, p=p, t_end=45.0, n=451))
        say("      %-28s %12.4f %12.2f"
            % (lab, float(r["AUCmuc"][-1]), m["dur_sev"]))
    say("")
    say("      A tissue that re-equilibrates with plasma in ~2 minutes cannot")
    say("      be protected by cutting its blood flow: the flow arm moves the")
    say("      RATE of equilibration, not the exposure.  Whatever cryotherapy")
    say("      does, the delivery story alone does not carry it.")

    say("")
    say("  (c) THE PATTERN THE GUIDELINES REPORT, RECOMPUTED.")
    say("      Same 30 min of ice against three fluoropyrimidine schedules")
    say("      that differ ONLY in how long the drug is present.")
    say("      %-34s %10s %10s %10s"
        % ("regimen", "AUC drop", "peak area", "sev durn"))
    T30 = 30.0 / (60.0 * 24.0)
    regs = [
        ("5-FU bolus 425 mg/m2 d1-5", lambda: M.sched_5FU_bolus(),
         [(d - 0.01, d + T30) for d in (0, 1, 2, 3, 4)], 40.0),
        ("melphalan 200 mg/m2 (30 min)", lambda: M.sched_HDM(200.0),
         [(-0.02, T30)], 45.0),
        ("5-FU 96-h infusion 4000 mg/m2", lambda: M.sched_5FU_CI(),
         [(-0.02, T30)], 45.0),
    ]
    for lab, mk, cw, te in regs:
        b = run(mk(), t_end=te)
        bm = M.metrics(b)
        s = M.add_cryo(mk(), cw)
        a = run(s, t_end=te)
        am = M.metrics(a)
        drop = 1.0 - float(a["AUCmuc"][-1]) / max(float(b["AUCmuc"][-1]), 1e-12)
        say("      %-34s %9.1f%% %10.3f %10.2f"
            % (lab, 100 * drop, am["peak_area"], am["dur_sev"]))
        say("      %-34s %9s %10.3f %10.2f"
            % ("   (no ice)", "-", bm["peak_area"], bm["dur_sev"]))


# ============================================================================
def s5_palifermin():
    hdr("5.  PALIFERMIN -- one fitted number, then predictions")
    say("")
    say("  Stage 2 spent exactly ONE number: the placebo arm's median grade")
    say("  3/4 duration (9 d), used to size the etoposide/cyclophosphamide")
    say("  block as a melphalan equivalent.  The palifermin arm below is a")
    say("  PREDICTION.")
    say("")
    import om_calibrate as C
    # Conditioning begins at t = 3 (see sched_spielberger's t0), so the label
    # schedule -- three daily doses BEFORE conditioning and three after stem
    # cells -- lands on days 0,1,2 and 11,12,13.
    pl = C.sched_spielberger(CY_EQ)
    pl_pal = C.sched_spielberger(CY_EQ,
                                 pal_days=[0.0, 1.0, 2.0, 11.0, 12.0, 13.0])
    a = mx(pl, t_end=48.0)
    b = mx(pl_pal, t_end=48.0)
    say("  %-38s %10s %10s %10s"
        % ("", "placebo", "palifermin", "observed"))
    say("  %-38s %10.2f %10.2f %10s"
        % ("duration WHO>=3 (d)", a["dur_sev"], b["dur_sev"], "9 -> 3"))
    say("  %-38s %10.2f %10.2f %10s"
        % ("peak WHO grade", a["peak_who"], b["peak_who"], "4 -> 3/4"))
    say("  %-38s %10.2f %10.2f %10s"
        % ("peak ulcer area", a["peak_area"], b["peak_area"], "-"))
    say("  %-38s %10.2f %10.2f %10s"
        % ("days on opioid", a["opidays"], b["opidays"], "11 -> 7"))
    say("  %-38s %10.2f %10.2f %10s"
        % ("onset WHO>=3 (d)", a["onset_sev"], b["onset_sev"], "-"))
    say("")
    say("  NOTE the SHAPE of the prediction, which is the two-clock claim")
    say("  again: palifermin is a REGENERATION-arm drug, so it should move")
    say("  DURATION much more than ONSET.  Check the two rows above.")

    say("")
    say("  (b) THE SCHEDULING PARADOX -- WHERE IT IS NOT.")
    say("      Sweep the gap between the last pre-conditioning palifermin")
    say("      dose and the alkylator.  Negative gap = palifermin still on")
    say("      board when the cytotoxic lands.")
    say("      %-16s %12s %12s %12s"
        % ("gap (d)", "sev durn d", "peak area", "vs placebo"))
    say("      %-16s %12.2f %12.3f %12s"
        % ("no palifermin", a["dur_sev"], a["peak_area"], "-"))
    for gap in [-1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]:
        # the alkylator lands at t = 3 + 4 = 7 in the Spielberger schedule
        last = 7.0 - gap
        days = [last - 2.0, last - 1.0, last, 11.0, 12.0, 13.0]
        days = [max(d, 0.0) for d in days]
        s = C.sched_spielberger(CY_EQ, pal_days=days)
        m = mx(s, t_end=48.0)
        say("      %-16.1f %12.2f %12.3f %11.1f%%"
            % (gap, m["dur_sev"], m["peak_area"],
               100 * (m["dur_sev"] - a["dur_sev"]) / max(a["dur_sev"], 1e-9)))
    say("")
    say("      A PRIOR EXPECTATION OF THIS MODEL'S AUTHOR IS REFUTED HERE,")
    say("      AND THE SIGN IS THE OPPOSITE OF THE ONE EXPECTED.  The map and")
    say("      the model header both claimed the label's 24-hour separation")
    say("      is explained by KGF enlarging the cycling pool.  In THIS")
    say("      regimen it cannot be: neither an alkylator nor a photon")
    say("      carries a cycle-specific term (cycspec_mel = cycspec_cis = 0,")
    say("      and the LQ kill is first-order in S), so a larger cycling pool")
    say("      is a proportionally larger kill on a proportionally larger")
    say("      pool and the FRACTION killed does not change.  With the")
    say("      penalty structurally absent, all that is left is the KGF")
    say("      effect pool's 2-day half-life -- so giving palifermin CLOSER")
    say("      to the insult is strictly BETTER, monotonically, across the")
    say("      whole sweep.  The model therefore does not reproduce the")
    say("      label's separation requirement from mucosal kinetics, and it")
    say("      should not be read as supporting it: whatever justifies that")
    say("      requirement lies outside this model (formulation, direct")
    say("      interaction, or tumour safety), not in the epithelium.")

    say("")
    say("  (c) THE SCHEDULING PARADOX -- WHERE IT IS.")
    say("      Repeat with a CYCLE-ACTIVE insult: bolus 5-FU 425 mg/m2 d1-5,")
    say("      for which cycspec = 1 and the kill scales with the cycling")
    say("      fraction.  Endpoint is ULCERATIVE (WHO>=2) stomatitis, the")
    say("      relevant grade for this regimen.")
    fb = mx(M.sched_5FU_bolus(), t_end=40.0, sev=2.0)
    p0 = dict(P)
    p0["fcyc_KGF"] = 0.0
    say("      %-16s %12s %12s %14s"
        % ("gap (d)", "durn>=2 d", "peak area", "fcyc_KGF = 0"))
    say("      %-16s %12.2f %12.3f %14s"
        % ("no palifermin", fb["dur_sev"], fb["peak_area"], "-"))
    for gap in [-1.0, -0.5, 0.0, 1.0, 2.0, 3.0]:
        days = [max(-gap - 2.0, 0.0), max(-gap - 1.0, 0.0),
                max(-gap, 0.0), 7.0, 8.0, 9.0]
        s = M.add_palifermin(M.sched_5FU_bolus(), days)
        m = mx(s, t_end=40.0, sev=2.0)
        m0 = mx(s, t_end=40.0, p=p0, sev=2.0)
        say("      %-16.1f %12.2f %12.3f %14.2f"
            % (gap, m["dur_sev"], m["peak_area"], m0["dur_sev"]))
    say("")
    say("      THE LAST COLUMN IS THE RESULT.  Turning off the ONE parameter")
    say("      that couples KGF-driven proliferation to the cytotoxic target")
    say("      turns palifermin from a marginal drug into a curative one in")
    say("      this regimen.  The coupling eats almost the whole benefit.")
    say("")
    say("      And the sweep is flat again -- but for the opposite reason.")
    say("      A 5-day bolus course cannot be separated from a growth factor")
    say("      whose effect pool has a 2-day half-life: there is no gap that")
    say("      empties the pool before the next dose lands.  That, not a")
    say("      24-hour rule, is the model's account of why palifermin is used")
    say("      in transplant conditioning and not with conventional")
    say("      chemotherapy: conditioning is an IMPULSE that a growth factor")
    say("      can be scheduled around, and a multi-day cytotoxic course is")
    say("      not.")


# ============================================================================
def s6_fractionation():
    hdr("6.  FRACTIONATION AND THE COST OF A TREATMENT GAP")
    say("")
    say("  All arms deliver a clinically equivalent tumour prescription; they")
    say("  differ in HOW the insult is spread in time.  The mucosal cost and")
    say("  the tumour BED are computed from the same trajectory.")
    say("")
    say("  %-34s %8s %8s %8s %8s"
        % ("schedule", "onset", "durn", "peakA", "BED_t"))
    arms = [
        ("70 Gy / 35 fx / 7 wk (standard)",
         dict(total_Gy=70.0, nfx=35, dose_per_fx=2.0, fx_per_week=5)),
        ("81.6 Gy / 68 fx b.i.d. over 7 wk",
         dict(total_Gy=81.6, nfx=68, dose_per_fx=1.2, fx_per_week=5,
              fx_per_day=2)),
        ("70 Gy / 35 fx / 6 wk (6 fx/wk accel)",
         dict(total_Gy=70.0, nfx=35, dose_per_fx=2.0, fx_per_week=6)),
        ("54 Gy / 36 fx b.i.d. / 12 d (CHART)",
         dict(total_Gy=54.0, nfx=36, dose_per_fx=1.5, fx_per_week=7,
              fx_per_day=3, bid_gap_h=6.0, cis_mgm2=0.0)),
        ("70 Gy / 35 fx, RT alone (no cisplatin)",
         dict(total_Gy=70.0, nfx=35, dose_per_fx=2.0, fx_per_week=5,
              cis_mgm2=0.0)),
    ]
    for lab, kw in arms:
        s = M.sched_chemoRT(**kw)
        m = mx(s, t_end=115.0)
        say("  %-34s %8.1f %8.1f %8.3f %8.1f"
            % (lab, m["onset_sev"], m["dur_sev"], m["peak_area"], m["BEDt"]))

    say("")
    say("  A 7-day unplanned gap: mucosa recovers, tumour repopulates.")
    say("  %-24s %8s %8s %8s" % ("gap", "durn d", "peakA", "BED_t"))
    for g in [None, (21.0, 28.0), (28.0, 35.0), (35.0, 42.0)]:
        s = M.sched_chemoRT(gap=g)
        m = mx(s, t_end=125.0)
        say("  %-24s %8.1f %8.3f %8.1f"
            % ("none" if g is None else "d%d-%d" % g, m["dur_sev"],
               m["peak_area"], m["BEDt"]))


# ============================================================================
def s7_who_saturation():
    hdr("7.  WHY WHO-GRADED TRIALS LOSE POWER WHERE THE DISEASE IS WORST")
    say("")
    say("  The WHO scale is ORDINAL and grade 4 is absorbing.  The same")
    say("  intervention is measured here on WHO grade and on ulcer area, at")
    say("  three baseline severities produced by scaling `sens` alone.")
    say("")
    say("  intervention: 6 h oral cryotherapy with HDM")
    say("  %-14s %10s %10s %10s %10s %10s"
        % ("sens", "areaAUC0", "areaAUC1", "d(area)%", "WHOdays0",
           "d(WHOd)%"))
    ice = [(-0.02, 6.0 / 24.0)]
    for sfac in [0.7, 1.0, 1.3, 1.7]:
        p = dict(P)
        p["sens"] = sfac
        r0 = run(M.sched_HDM(200.0), t_end=45.0, p=p)
        r1 = run(M.add_cryo(M.sched_HDM(200.0), ice), t_end=45.0, p=p)
        m0, m1 = M.metrics(r0), M.metrics(r1)
        a0, a1 = m0["area_auc"], m1["area_auc"]
        w0, w1 = m0["dur_sev"], m1["dur_sev"]
        say("  %-14.2f %10.2f %10.2f %9.1f%% %10.2f %9.1f%%"
            % (sfac, a0, a1, 100 * (a1 - a0) / max(a0, 1e-9), w0,
               100 * (w1 - w0) / max(w0, 1e-9)))
    say("")
    say("  THE ARTEFACT IS REAL BUT IT IS A FLOOR, NOT THE CEILING THIS")
    say("  SECTION WAS WRITTEN TO LOOK FOR.  The WHO-derived effect reads")
    say("  -100% at every severity -- the ordinal endpoint cannot express")
    say("  anything between 'severe' and 'not severe', so an intervention")
    say("  that drops the patient just below threshold and one that abolishes")
    say("  the disease score identically.  The continuous endpoint separates")
    say("  them cleanly and shows the effect DECAYING with severity")
    say("  (-97.8% -> -84.9%), which is the clinically meaningful fact and")
    say("  the one the ordinal scale destroys.  Either way the conclusion")
    say("  for trial design is the same and it is not the one asserted")
    say("  earlier: WHO grade is not a low-powered version of ulcer area, it")
    say("  is a DIFFERENT measurement, and a trial powered on it is blind to")
    say("  the size of the effect it is measuring.")


# ============================================================================
def s8_arms():
    hdr("8.  AN AGENT CANNOT BUY THE ENDPOINT IT DOES NOT OWN")
    say("")
    say("  Every available intervention, HDM 200 mg/m2, same patient.")
    say("  %-30s %8s %8s %8s %8s"
        % ("intervention", "onset", "durn", "peakA", "painAUC"))
    base = M.sched_HDM(200.0)
    b = mx(base, t_end=45.0)
    say("  %-30s %8.2f %8.2f %8.3f %8.1f"
        % ("none", b["onset_sev"], b["dur_sev"], b["peak_area"],
           b["painAUC"]))
    import om_calibrate as C
    opts = [
        ("cryotherapy 6 h  [INSULT]",
         M.add_cryo(M.sched_HDM(200.0), [(-0.02, 0.25)])),
        ("dose 140 mg/m2   [INSULT]", M.sched_HDM(140.0)),
        ("palifermin       [REGEN ]",
         M.add_palifermin(M.sched_HDM(200.0), [-3.0, -2.0, -1.0,
                                               3.0, 4.0, 5.0])),
        ("photobiomodulation [REGEN]",
         M.add_pbm(M.sched_HDM(200.0), list(np.arange(0.0, 21.0, 1.0)))),
        ("glutamine        [REGEN ]",
         M.add_gln(M.sched_HDM(200.0), 0.0, 21.0)),
        ("benzydamine      [INFLAM]",
         M.add_bzd(M.sched_HDM(200.0), 0.0, 21.0)),
        ("cryo + palifermin",
         M.add_palifermin(M.add_cryo(M.sched_HDM(200.0), [(-0.02, 0.25)]),
                          [-3.0, -2.0, -1.0, 3.0, 4.0, 5.0])),
    ]
    for lab, s in opts:
        m = mx(s, t_end=45.0)
        say("  %-30s %8.2f %8.2f %8.3f %8.1f"
            % (lab, m["onset_sev"], m["dur_sev"], m["peak_area"],
               m["painAUC"]))
    say("")
    say("  READ THE AREA COLUMN, NOT THE DURATION COLUMN.  Section 7's floor")
    say("  effect reappears here unbidden: the calibrated patient sits just")
    say("  above the WHO>=3 threshold (peak area 0.68 against a threshold of")
    say("  0.35), so every effective agent reads as 100% prevention on the")
    say("  ordinal endpoint and 'nan' onset.  That is an artefact of the")
    say("  instrument and of a single deterministic patient -- NOT a claim")
    say("  that six-hour ice abolishes mucositis, which no trial reports.")
    say("  On the continuous endpoint the agents separate and rank sensibly,")
    say("  and the INSULT-arm agents (ice, dose reduction) beat the")
    say("  REGENERATION-arm agents (palifermin, glutamine) on peak area,")
    say("  which is what a depth-vs-recovery-rate reading predicts.")
    say("  The virtual population in om_calibrate.py, not this table, is")
    say("  where incidence-type claims belong.")


# ============================================================================
def s9_infection():
    hdr("9.  THE ULCER AS A PORTAL OF ENTRY")
    say("")
    say("  Bacteraemia hazard integrates A_ulc x microbial load / ANC.  The")
    say("  ANC arm and the mucosal arm are driven by the SAME melphalan")
    say("  exposure but with completely different time constants, so the")
    say("  hazard is a product of two curves that peak at different times.")
    say("")
    say("  %-22s %10s %10s %10s %10s"
        % ("scenario", "ANC nadir", "peak area", "cum hazard", "durn"))
    for lab, s, p in [
            ("HDM alone", M.sched_HDM(200.0), P),
            ("HDM + cryo 6 h",
             M.add_cryo(M.sched_HDM(200.0), [(-0.02, 0.25)]), P),
            ("HDM + palifermin",
             M.add_palifermin(M.sched_HDM(200.0),
                              [-3.0, -2.0, -1.0, 3.0, 4.0, 5.0]), P)]:
        m = mx(s, t_end=45.0, p=p)
        say("  %-22s %10.2f %10.3f %10.3f %10.2f"
            % (lab, m["nadir_anc"], m["peak_area"], m["infhaz"],
               m["dur_sev"]))
    say("")
    r = run(M.sched_HDM(200.0), t_end=45.0)
    i_a = int(np.argmax(r["A"]))
    i_n = int(np.argmin(r["ANC"]))
    say("  peak ulcer area at day %.1f;  ANC nadir at day %.1f;  gap %.1f d"
        % (r["t"][i_a], r["t"][i_n], abs(r["t"][i_a] - r["t"][i_n])))


# ============================================================================
def s10_loopgain():
    hdr("10.  STRUCTURAL CHECKS")
    say("")
    say("  (a) The TNF -> NF-kB positive feedback.  The loop gain at the")
    say("      operating point must be < 1 or the inflammation block")
    say("      diverges.  Gain = gTNF * d(hill)/d(TNF) * pTNF/kTNF.")
    for nf in [0.2, 0.5, 1.0, 1.5]:
        tnf = P["pTNF"] * nf / P["kTNF"]
        dh = P["KTNF"] / (P["KTNF"] + tnf) ** 2
        g = P["gTNF"] * dh * P["pTNF"] / P["kTNF"]
        say("      NFkB=%.2f -> TNF=%.3f  loop gain %.4f  %s"
            % (nf, tnf, g, "OK" if g < 1 else "*** DIVERGENT ***"))

    say("")
    say("  (b) Dose delivery is exact (the breakpoint integrator).")
    for lab, s, te, exp in [
            ("70 Gy/35 fx mucosal BED (Gy10)", M.sched_chemoRT(), 60.0,
             70.0 * 1.2),
            ("TBI 12 Gy/8 fx mucosal BED (Gy10)",
             {**M.ZERO_SCHED, "rt": [(i * 0.5, 1.5, 10.0 / 1440.0)
                                     for i in range(8)],
              "dose_per_fx": 1.5}, 20.0, 12.0 * 1.15)]:
        r = run(s, t_end=te)
        say("      %-38s %10.4f  expected %8.4f"
            % (lab, float(r["BEDm"][-1]), exp))

    say("")
    say("  (c) Integrator convergence: halve the tolerances and compare.")
    import om_python_reference as MM
    s = M.sched_chemoRT()
    r1 = M.simulate(s, p=P, t_end=60.0, n=601)
    # tighter tolerance run
    src_ok = True
    try:
        import scipy.integrate as si
        old = MM.solve_ivp

        def tight(f, span, y0, **kw):
            kw["rtol"] = 1e-9
            kw["atol"] = 1e-11
            return old(f, span, y0, **kw)
        MM.solve_ivp = tight
        r2 = M.simulate(s, p=P, t_end=60.0, n=601)
        MM.solve_ivp = old
    except Exception as e:
        src_ok = False
        say("      tolerance check failed: %s" % e)
    if src_ok:
        d = np.max(np.abs(r1["A"] - r2["A"]))
        say("      max |A(t)| difference rtol 1e-7 vs 1e-9 : %.3e" % d)
        rel = (abs(float(r1["BEDm"][-1]) - float(r2["BEDm"][-1]))
               / max(float(r1["BEDm"][-1]), 1e-12))
        say("      relative BED difference                  : %.3e" % rel)

    say("")
    say("  (d) Drug-free steady state is a genuine attractor.")
    ss = M.equilibrate(P, tmax=400.0)
    r = M.simulate({**M.ZERO_SCHED}, p=P, t_end=60.0, n=121)
    say("      D drift over 60 drug-free days: %.3e"
        % float(np.max(np.abs(r["D"] - r["D"][0]))))
    say("      max ulcer area with no insult : %.3e" % float(r["A"].max()))


# ============================================================================
def s11_sensitivity():
    hdr("11.  LOCAL SENSITIVITY OF THE HEADLINE OUTPUTS")
    say("")
    say("  +/-20% on each parameter, HDM 200 mg/m2.  Reported as the relative")
    say("  swing in severe-OM duration and in the ulcer-area integral.")
    say("")
    keys = ["lamS", "aS", "k_p", "k_shed", "greg", "Kreg", "k_ab", "ndiv",
            "Sq0", "q_res", "kact", "krest",
            "Dcrit", "fmax", "acc_kp", "pot_mel", "keq_muc", "Kp_mel",
            "aCER", "aTNFs", "gTNF", "gMB", "pCER", "eKGF", "fcyc_KGF",
            "kKGFe", "Q10_cryo", "fQ_cryo", "sens"]
    base = mx(M.sched_HDM(200.0), t_end=45.0)
    rows = []
    for k in keys:
        if k not in P:
            continue
        vals = []
        for f in (0.8, 1.2):
            p = dict(P)
            p[k] = P[k] * f
            try:
                vals.append(mx(M.sched_HDM(200.0), t_end=45.0, p=p))
            except Exception:
                vals.append(None)
        if None in vals:
            continue
        dd = abs(vals[1]["dur_sev"] - vals[0]["dur_sev"]) \
            / max(base["dur_sev"], 1e-9)
        da = abs(vals[1]["area_auc"] - vals[0]["area_auc"]) \
            / max(base["area_auc"], 1e-9)
        rows.append((k, dd, da))
    rows.sort(key=lambda r: -r[1])
    say("  %-14s %12s %12s" % ("parameter", "d durn", "d areaAUC"))
    for k, dd, da in rows:
        say("  %-14s %11.1f%% %11.1f%%" % (k, 100 * dd, 100 * da))


# ============================================================================
def s12_scenarios():
    hdr("12.  THE SEVENTEEN SCENARIOS SHIPPED IN THE mrgsolve MODEL")
    say("")
    say("  %-3s %-44s %7s %7s %7s %8s"
        % ("#", "scenario", "onset", "durn", "peakA", "painAUC"))
    import om_calibrate as C
    pal6 = [-3.0, -2.0, -1.0, 3.0, 4.0, 5.0]
    scen = [
        ("HDM 200 mg/m2, no prophylaxis", M.sched_HDM(200.0), 45.0),
        ("HDM 140 mg/m2 (reduced dose)", M.sched_HDM(140.0), 45.0),
        ("HDM 200 + cryotherapy 30 min",
         M.add_cryo(M.sched_HDM(200.0), [(-0.02, 30. / 1440)]), 45.0),
        ("HDM 200 + cryotherapy 6 h",
         M.add_cryo(M.sched_HDM(200.0), [(-0.02, 0.25)]), 45.0),
        ("HDM 200 + palifermin (separated)",
         M.add_palifermin(M.sched_HDM(200.0), pal6), 45.0),
        ("HDM 200 + palifermin (concurrent)",
         M.add_palifermin(M.sched_HDM(200.0), [-1.0, -0.5, 0.0,
                                               3.0, 4.0, 5.0]), 45.0),
        ("HDM 200 + photobiomodulation daily",
         M.add_pbm(M.sched_HDM(200.0), list(np.arange(0., 21., 1.))), 45.0),
        ("HDM 200 + glutamine",
         M.add_gln(M.sched_HDM(200.0), 0.0, 21.0), 45.0),
        ("HDM 200 + cryo + palifermin",
         M.add_palifermin(M.add_cryo(M.sched_HDM(200.0), [(-0.02, 0.25)]),
                          pal6), 45.0),
        ("TBI-VP16-Cy conditioning (placebo)",
         C.sched_spielberger(CY_EQ), 48.0),
        ("TBI-VP16-Cy + palifermin",
         C.sched_spielberger(CY_EQ, pal_days=[0.0, 1.0, 2.0,
                                              11.0, 12.0, 13.0]), 48.0),
        ("H&N 70 Gy/35 fx + cisplatin", M.sched_chemoRT(), 115.0),
        ("H&N 70 Gy/35 fx, RT alone",
         M.sched_chemoRT(cis_mgm2=0.0), 115.0),
        ("H&N hyperfractionated 81.6 Gy/68 fx b.i.d.",
         M.sched_chemoRT(total_Gy=81.6, nfx=68, dose_per_fx=1.2,
                         fx_per_week=5, fx_per_day=2), 115.0),
        ("H&N chemoRT + benzydamine",
         M.add_bzd(M.sched_chemoRT(), 0.0, 60.0), 115.0),
        ("5-FU bolus 425 mg/m2 d1-5", M.sched_5FU_bolus(), 40.0),
        ("5-FU 96-h CI 4000 mg/m2", M.sched_5FU_CI(), 40.0),
    ]
    for i, (lab, s, te) in enumerate(scen, 1):
        m = mx(s, t_end=te)
        say("  %-3d %-44s %7.2f %7.2f %7.3f %8.1f"
            % (i, lab, m["onset_sev"], m["dur_sev"], m["peak_area"],
               m["painAUC"]))


# ============================================================================
if __name__ == "__main__":
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    steps = [("1", s1_calibration), ("2", s2_two_clocks), ("3", s3_latency),
             ("4", s4_cryo), ("5", s5_palifermin), ("6", s6_fractionation),
             ("7", s7_who_saturation), ("8", s8_arms), ("9", s9_infection),
             ("10", s10_loopgain), ("11", s11_sensitivity),
             ("12", s12_scenarios)]
    for k, f in steps:
        if only and k not in only:
            continue
        try:
            f()
        except Exception as e:
            import traceback
            say("")
            say("*** SECTION %s FAILED: %r" % (k, e))
            say(traceback.format_exc())
    say("")
    say("done.")
