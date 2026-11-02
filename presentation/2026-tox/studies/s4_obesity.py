"""Study 4 — Obesity: incretin-therapy titration as a three-way optimization
of weight loss, GI tolerability, and lean-mass (safety), with an elderly /
sarcopenic subgroup.

Drugs: semaglutide (GLP-1 RA) and tirzepatide (GLP-1/GIP). The clinical levers
are the TARGET dose and the ESCALATION schedule. They trade:
  benefit      = total body-weight loss (appetite suppression -> energy deficit)
  tolerability = peak GI adverse-event burden (nausea/vomiting; transient,
                 driven by the RATE of dose escalation, blunted by tolerance)
  safety       = absolute lean (fat-free) mass lost — the axis that matters in
                 older adults, where ~1/3 of excess risk is sarcopenia/frailty

QSP layer:
  * exposure accumulates toward the steady state for the current dose,
  * appetite suppression = Hill(exposure) -> caloric intake,
  * body weight via an energy-balance ODE; weight change partitioned into fat
    and fat-free mass by a Forbes-type rule (age-dependent),
  * GI AE = nausea driven by d(exposure)/dt and level, damped by a tolerance
    state that builds with sustained exposure.

Calibration anchors:
  STEP 1 (Wilding NEJM 2021): semaglutide 2.4 mg -> -14.9% at 68 wk; ~1/3 lean.
  SURMOUNT-1 (Jastreboff NEJM 2022): tirzepatide 5/10/15 mg -> -15/-19.5/-20.9%.
  GI AEs mostly mild-moderate, concentrated in the escalation phase.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import odeint
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as F
F.apply()
import matplotlib.pyplot as plt

OUT = Path("/home/user/qsp/presentation/2026-tox/data"); OUT.mkdir(exist_ok=True, parents=True)

# ---- drug definitions -----------------------------------------------------
# appetite suppression Emax (max fractional intake reduction) and dose EC50
DRUGS = {
    "semaglutide":  dict(Emax=0.285, EC50=0.9,  maxdose=2.4,  color=F.TEAL),
    "tirzepatide":  dict(Emax=0.40,  EC50=6.0,  maxdose=15.0, color=F.NAVY),
}
KE = 1/7.0   # exposure equilibration rate (1/day); long-acting weekly agents

def esc_schedule(drug, target, weeks_step=4, step_frac=None):
    """return a function week->dose (mg), 4-weekly escalation to target."""
    md = DRUGS[drug]["maxdose"]
    if drug == "semaglutide":
        steps = [0.25, 0.5, 1.0, 1.7, 2.4]
    else:
        steps = [2.5, 5.0, 7.5, 10.0, 12.5, 15.0]
    steps = [s for s in steps if s <= target] or [steps[0]]
    if steps[-1] < target: steps.append(target)
    def dose(week):
        idx = min(int(week // weeks_step), len(steps)-1)
        return steps[idx]
    return dose, steps

def simulate(drug, target, age=45, weeks=72, weeks_step=4, BW0=105.0):
    d = DRUGS[drug]
    dose_fn, steps = esc_schedule(drug, target, weeks_step)
    # baseline composition (fraction fat by age/BMI); elderly: less FFM reserve
    fat_frac0 = 0.42 + 0.0015*(age-45)
    FM0 = BW0*fat_frac0; FFM0 = BW0 - FM0
    EE0 = 22.0*FFM0 + 500          # kcal/d resting+activity (Cunningham-ish)
    EI0 = EE0                       # weight-stable baseline
    t = np.arange(0, weeks*7+1)     # days
    # exposure trajectory
    C = np.zeros_like(t, dtype=float)
    for k in range(1, len(t)):
        wk = t[k]/7.0
        Css = dose_fn(wk)
        C[k] = C[k-1] + KE*(Css - C[k-1])
    # appetite suppression fraction = Emax * dose-occupancy
    supp = d["Emax"] * (C/(C + d["EC50"]))
    dose_series = np.array([dose_fn(wk/7.0) for wk in t])
    # ---- energy balance ODE (daily) ----
    FM = np.zeros_like(t, float); FFM = np.zeros_like(t, float); BW = np.zeros_like(t, float)
    FM[0], FFM[0], BW[0] = FM0, FFM0, BW0
    # Forbes partition: fraction of weight change that is FFM
    for k in range(1, len(t)):
        EI = EI0*(1 - supp[k])
        # expenditure falls as mass falls + adaptive thermogenesis
        EE = 22.0*FFM[k-1] + 500*(BW[k-1]/BW0) - 120*supp[k-1]
        bal = EI - EE                      # kcal/d (negative = deficit)
        # Forbes: dFFM/dFM ratio depends on fat mass; elderly lose more lean
        p_lean = (20.0/(20.0+FM[k-1]))      # fraction FFM of change (Forbes-type)
        age_mult = 1 + 0.018*(age-45)       # anabolic resistance -> more lean loss
        p_lean = min(0.6, p_lean*age_mult)
        rho = 9441*(1-p_lean) + 1816*p_lean # kcal/kg composite tissue energy
        dBW = bal/rho
        FFM[k] = FFM[k-1] + dBW*p_lean
        FM[k]  = FM[k-1] + dBW*(1-p_lean)
        BW[k]  = FFM[k] + FM[k]
    # ---- GI adverse events (nausea) with tolerance ----
    # Driver on the 0-1 OCCUPANCY scale so the two drugs are comparable; nausea
    # rises with both the rate of occupancy change (escalation) and its level,
    # and is blunted by a tolerance state that builds with sustained exposure.
    occ = C/(C + d["EC50"])
    docc = np.gradient(occ, t)
    N = np.zeros_like(t, float); Tol = np.zeros_like(t, float)
    for k in range(1, len(t)):
        drive = 22.0*max(docc[k], 0) + 0.9*occ[k]
        Tol[k] = Tol[k-1] + (0.045*occ[k]*(1-Tol[k-1]) - 0.004*Tol[k-1])
        N[k] = N[k-1] + (drive*(1-0.85*Tol[k-1]) - 0.22*N[k-1])
    gi_burden = np.trapezoid(N[t <= 24*7], t[t <= 24*7])   # cumulative over escalation phase
    return dict(t=t/7.0, C=C, occ=occ, supp=supp, BW=BW, FM=FM, FFM=FFM, N=N, Tol=Tol,
                steps=steps, dose_series=dose_series, FM0=FM0, FFM0=FFM0, BW0=BW0,
                pct=(BW0-BW)/BW0*100, lean_lost=FFM0-FFM[-1], fat_lost=FM0-FM[-1],
                final_pct=(BW0-BW[-1])/BW0*100, peakN=N.max(), gi_burden=gi_burden,
                plateauN=N[-1])

# calibrate check
for drug, tgt in [("semaglutide",2.4),("tirzepatide",15.0),("tirzepatide",5.0),("tirzepatide",10.0)]:
    r = simulate(drug, tgt)
    print(f"{drug} {tgt}mg -> {r['final_pct']:.1f}% weight, lean {r['lean_lost']:.1f}kg "
          f"({r['lean_lost']/(r['lean_lost']+r['fat_lost'])*100:.0f}% of loss)")

# ============================================================
# FIG A — titration: exposure, appetite suppression, GI AE (tolerance)
# ============================================================
r_slow = simulate("semaglutide", 2.4, weeks_step=4)
r_fast = simulate("semaglutide", 2.4, weeks_step=2)
fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
ax = axes[0]
ax2 = ax.twinx()
ax2.step(r_slow["t"], r_slow["dose_series"], color=F.GREY, lw=1.6, where="post", zorder=1)
ax2.set_ylabel("dose (mg)", color=F.GREY); ax2.set_ylim(0, 3)
ax2.spines["top"].set_visible(False); ax2.tick_params(colors=F.GREY)
ax.plot(r_slow["t"], r_slow["supp"]*100, color=F.TEAL, lw=2.8, zorder=3,
        label="appetite suppression")
ax.set_zorder(ax2.get_zorder()+1); ax.patch.set_visible(False)
ax.set_xlabel("weeks"); ax.set_ylabel("appetite suppression (%)")
ax.set_title("Exposure builds with each dose step", fontsize=13)
ax.set_ylim(0, 40); ax.set_xlim(0, 72)
ax.legend(loc="lower right", fontsize=10.5); F.despine(ax)
ax = axes[1]
nref_A = max(r_slow["N"].max(), r_fast["N"].max())   # shared scale (both panels)
ax.plot(r_slow["t"], r_slow["N"]/nref_A*100, color=F.GOLD, label="4-weekly escalation")
ax.plot(r_fast["t"], r_fast["N"]/nref_A*100, color=F.RED, label="2-weekly (too fast)")
ax.set_xlabel("weeks"); ax.set_ylabel("nausea intensity (rel.)")
ax.set_title("Nausea tracks escalation, then adapts", fontsize=13)
ax.legend(loc="upper right", fontsize=10); ax.set_xlim(0, 40); F.despine(ax)
fig.subplots_adjust(wspace=0.42)
fig.suptitle("Titration is the handle: faster escalation buys nothing but nausea",
             fontsize=15, fontweight="bold", color=F.NAVY, y=1.03)
F.save(fig, "s4_titration.png")

# ============================================================
# FIG B — efficacy vs tolerability frontier over target dose
# ============================================================
fig, ax = plt.subplots(figsize=(9.2, 5.6))
# global nausea reference across the whole sweep
sweep = {}
for drug in DRUGS:
    md = DRUGS[drug]["maxdose"]
    sweep[drug] = [(tg, simulate(drug, tg)) for tg in np.linspace(md*0.15, md, 10)]
NREF = max(r["gi_burden"] for runs in sweep.values() for _, r in runs)
for drug, mark in [("semaglutide","o"),("tirzepatide","s")]:
    eff = [r["final_pct"] for _, r in sweep[drug]]
    tol = [r["gi_burden"]/NREF*100 for _, r in sweep[drug]]
    ax.plot(tol, eff, "-", color=DRUGS[drug]["color"], marker=mark, ms=7,
            label=drug, lw=2.4)
# mark approved doses
for drug, tg, lab, off in [("semaglutide",2.4,"sema 2.4 mg",(10,-4)),
                           ("tirzepatide",15.0,"tirze 15 mg",(-18,12)),
                           ("tirzepatide",5.0,"tirze 5 mg",(10,-4)),
                           ("tirzepatide",10.0,"tirze 10 mg",(10,-4))]:
    r = simulate(drug, tg)
    ax.scatter([r["gi_burden"]/NREF*100],[r["final_pct"]], s=150, color=DRUGS[drug]["color"],
               edgecolors="white", linewidths=1.5, zorder=5)
    ax.annotate(lab, (r["gi_burden"]/NREF*100, r["final_pct"]), textcoords="offset points",
                xytext=off, fontsize=10.5, fontweight="bold", color=DRUGS[drug]["color"])
ax.set_xlabel("cumulative GI adverse-event burden, escalation phase (relative)")
ax.set_ylabel("total body-weight loss at 72 wk (%)")
ax.set_title("More weight loss costs more GI burden — the efficacy–tolerability frontier", pad=10, fontsize=14)
ax.legend(loc="lower right", title="incretin"); F.despine(ax)
F.save(fig, "s4_frontier.png")

# ============================================================
# FIG C — body composition: fat vs lean, young vs elderly
# ============================================================
fig, ax = plt.subplots(figsize=(9.2, 5.4))
groups = [("Age 45", 45, 0), ("Age 70", 70, 1.2), ("Age 70\n+ muscle-preserving", 70, 2.4)]
x = np.arange(len(groups)); w = 0.55
for i,(lab, age, mode) in enumerate(groups):
    r = simulate("tirzepatide", 15.0, age=age)
    fat, lean = r["fat_lost"], r["lean_lost"]
    if mode == 2.4:  # muscle-preserving analog (resistance ex / slower): shift partition
        shift = lean*0.45; lean -= shift; fat += shift*0.4
    ax.bar(x[i], fat, w, color=F.GOLD, label="fat mass lost" if i==0 else "")
    ax.bar(x[i], lean, w, bottom=fat, color=F.RED, label="lean mass lost" if i==0 else "")
    ax.text(x[i], fat/2, f"{fat:.1f}", ha="center", va="center", color="white", fontweight="bold", fontsize=11)
    ax.text(x[i], fat+lean/2, f"{lean:.1f}", ha="center", va="center", color="white", fontweight="bold", fontsize=11)
    ax.text(x[i], fat+lean+0.4, f"{lean/(fat+lean)*100:.0f}% lean", ha="center", fontsize=10, color=F.INK)
ax.set_xticks(x); ax.set_xticklabels([g[0] for g in groups], fontsize=11)
ax.set_ylabel("mass lost at 72 wk (kg)")
ax.set_title("Same drug, different body-composition cost by age", pad=10)
ax.legend(loc="upper right"); ax.grid(axis="x", visible=False); F.despine(ax)
F.save(fig, "s4_composition.png")

# ============================================================
# FIG D — lean-mass safety across age (the elderly subgroup)
# ============================================================
ages = np.arange(40, 86, 3)
fig, ax = plt.subplots(figsize=(9.2, 5.2))
for drug, tg in [("tirzepatide", 15.0), ("semaglutide", 2.4)]:
    lean = [simulate(drug, tg, age=a)["lean_lost"] for a in ages]
    ax.plot(ages, lean, color=DRUGS[drug]["color"], marker="o", ms=4,
            label=f"{drug} max dose")
# sarcopenia-risk threshold (illustrative): ~8-10% of baseline FFM lost
r0 = simulate("tirzepatide", 15.0, age=70)
thr = 0.09*r0["FFM0"]
ax.axhline(thr, color=F.RED, ls="--", lw=1.6)
ax.text(41, thr+0.15, "sarcopenia-risk zone (>~9% of lean mass)", color=F.RED, fontsize=10.5)
ax.axhspan(thr, ax.get_ylim()[1], color=F.RED, alpha=0.05)
ax.set_xlabel("patient age (years)"); ax.set_ylabel("lean mass lost at 72 wk (kg)")
ax.set_title("Lean-mass loss climbs with age — a dose-reduction signal", pad=10)
ax.legend(loc="lower right"); F.despine(ax)
F.save(fig, "s4_aging.png")

res = {
    "sema24_pct": round(simulate("semaglutide",2.4)["final_pct"],1),
    "tirze15_pct": round(simulate("tirzepatide",15.0)["final_pct"],1),
    "tirze5_pct": round(simulate("tirzepatide",5.0)["final_pct"],1),
    "tirze15_lean_young": round(simulate("tirzepatide",15.0,age=45)["lean_lost"],1),
    "tirze15_lean_old": round(simulate("tirzepatide",15.0,age=70)["lean_lost"],1),
    "lean_frac_young": round(simulate("tirzepatide",15.0,age=45)["lean_lost"]/
                             (simulate("tirzepatide",15.0,age=45)["lean_lost"]+simulate("tirzepatide",15.0,age=45)["fat_lost"])*100,0),
    "lean_frac_old": round(simulate("tirzepatide",15.0,age=70)["lean_lost"]/
                           (simulate("tirzepatide",15.0,age=70)["lean_lost"]+simulate("tirzepatide",15.0,age=70)["fat_lost"])*100,0),
}
(OUT / "s4_results.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
print("S4 done")
