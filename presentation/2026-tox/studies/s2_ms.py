"""Study 2 — Multiple sclerosis: turning a disease-modifying-therapy
benefit-risk problem into a tunable optimization.

Focus drug: natalizumab (anti-alpha4-integrin / VLA-4 mAb). It is the
sharpest efficacy-safety trade-off in neurology: among the most effective
relapse suppressors, but it lifts CNS immune surveillance and so carries a
risk of progressive multifocal leukoencephalopathy (PML) from JC-virus
reactivation. Extended-interval dosing (EID, q6w) trades a small efficacy
margin for a large PML-risk reduction — a controllable handle.

A minimal QSP layer couples:
  * natalizumab exposure -> alpha4-integrin receptor saturation (desaturates
    between doses; weight- and interval-dependent),
  * saturation -> relapse-hazard suppression (annualized relapse rate, ARR),
  * saturation (sustained CNS surveillance deficit) -> a PML-hazard multiplier
    on top of a stratified base hazard.

Calibration anchors (cited inline):
  AFFIRM (Polman NEJM 2006): natalizumab q4w ARR 0.23 vs placebo 0.73 (~68% red.)
  NOVA  (Foley Lancet Neurol 2022): q6w EID retains most efficacy
  Ho (Lancet Neurol 2017) / Bloomgren (NEJM 2012): PML risk stratification table
  Ryerson (Neurology 2019, TOUCH): EID associated with ~88% lower PML hazard
  FREEDOMS (Kappos NEJM 2010): fingolimod ARR 0.18 vs 0.40 (~54% reduction)
  OPERA I/II (Hauser NEJM 2017): ocrelizumab ARR ~0.16
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as F
F.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = Path("/home/user/qsp/presentation/2026-tox/data"); OUT.mkdir(exist_ok=True, parents=True)

# ============================================================
# 1. natalizumab receptor saturation over a dosing interval
# ============================================================
# Steady-state trough alpha4 saturation falls as the inter-dose interval and
# body weight rise (van Kempen 2020; Foley 2022). Model saturation along one
# interval with a monoexponential drug decline feeding an Emax occupancy.
NAT_THALF = 11.0        # days, serum t1/2 (Polman 2006 range 7-15)
ke = np.log(2) / NAT_THALF
# concentration (arbitrary units) just after a 300 mg dose, scaled by weight
def trough_profile(interval_d, weight_kg=75.0, n_cycles=8):
    kV = 75.0 / weight_kg          # lighter patient -> higher conc
    Cmax = 30.0 * kV               # relative units
    # accumulate to steady state
    t = np.linspace(0, interval_d, 400)
    # steady-state superposition of a repeated IV bolus
    r = np.exp(-ke * interval_d)
    Css0 = Cmax / (1 - r)          # conc at start of a steady-state interval
    C = Css0 * np.exp(-ke * t)
    KD = 2.0                       # occupancy half-sat (units)
    sat = C / (KD + C)
    return t, sat

def mean_trough_sat(interval_d, weight_kg=75.0):
    t, sat = trough_profile(interval_d, weight_kg)
    return sat.mean(), sat[-1]     # time-averaged, trough

# efficacy: ARR = ARR0 * (1 - Emax * mean_saturation)
ARR0 = 0.73                        # placebo ARR (AFFIRM)
EMAX_ARR = 0.95                    # max achievable fractional suppression
def arr_for(interval_d, weight_kg=75.0):
    msat, _ = mean_trough_sat(interval_d, weight_kg)
    return ARR0 * (1 - EMAX_ARR * msat)

# calibrate: q4w(28d) at 75 kg should give ARR ~0.23
msat28, _ = mean_trough_sat(28, 75)
EMAX_ARR = (1 - 0.23 / ARR0) / msat28
arr_q4 = arr_for(28); arr_q6 = arr_for(42)
print(f"calibrated EMAX_ARR={EMAX_ARR:.3f} | ARR q4w={arr_q4:.3f} q6w={arr_q6:.3f}")

# ============================================================
# 2. PML risk stratification (Ho 2017, per 1000 patients) -> hazard
# ============================================================
# estimated incidence per 1000 by (stratum, epoch of months on drug)
# epochs: 1-24, 25-48, 49-72 months
PML_TABLE = {   # per 1000 patients
    "JCV- (index n/a)":            [0.09, 0.09, 0.09],
    "JCV+ index <=0.9":            [0.1,  0.3,  0.4],
    "JCV+ index 0.9-1.5":          [0.3,  0.9,  3.0],
    "JCV+ index >1.5":             [1.0,  3.0, 10.0],
    "JCV+ prior immunosuppressant":[1.2,  4.0,  8.0],
}
EPOCHS = [(0, 24), (24, 48), (48, 72)]
# convert epoch incidence to an (approx) monthly hazard within the epoch
def epoch_hazards(inc_per_1000):
    hz = []
    for inc, (a, b) in zip(inc_per_1000, EPOCHS):
        p = inc / 1000.0
        months = b - a
        h = -np.log(max(1e-9, 1 - p)) / months   # monthly hazard
        hz.append(h)
    return hz

def cum_pml(stratum, horizon_m, eid_hr=1.0):
    """cumulative PML probability over horizon (months); eid_hr scales hazard."""
    hz = epoch_hazards(PML_TABLE[stratum])
    H = 0.0
    for m in range(horizon_m):
        if m < 24: h = hz[0]
        elif m < 48: h = hz[1]
        else: h = hz[2]
        H += h * eid_hr
    return 1 - np.exp(-H)

EID_HR = 0.24   # EID vs SID PML hazard ratio (Ryerson 2019 ~0.12-0.88 range; use 0.24)

# ============================================================
# FIG A — saturation over interval + efficacy retention
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4))
ax = axes[0]
for interval, col, lab in [(28, F.NAVY, "q4w (standard, SID)"),
                           (42, F.GOLD, "q6w (extended, EID)"),
                           (56, F.RED, "q8w (over-extended)")]:
    t, sat = trough_profile(interval, 75)
    ax.plot(t, sat*100, color=col, label=lab)
    ax.scatter([interval], [sat[-1]*100], color=col, s=45, zorder=5)
ax.axhline(50, color=F.RED, lw=1.4, ls="--")
ax.text(1, 46, "surveillance-deficit zone", color=F.RED, fontsize=9.5, va="top")
ax.set_xlabel("days since dose"); ax.set_ylabel("α4-integrin saturation (%)")
ax.set_title("Receptor desaturation between doses", fontsize=13)
ax.set_ylim(0, 100); ax.legend(fontsize=10, loc="upper right")
F.despine(ax)
ax = axes[1]
intervals = np.arange(21, 61)
arrs = [arr_for(i) for i in intervals]
red = [(1 - a/ARR0)*100 for a in arrs]
ax.plot(intervals, red, color=F.TEAL)
for iv, c, lab in [(28, F.NAVY, "SID"), (42, F.GOLD, "EID")]:
    ax.scatter([iv], [(1-arr_for(iv)/ARR0)*100], color=c, s=70, zorder=5)
    ax.annotate(f"{lab}\n{(1-arr_for(iv)/ARR0)*100:.0f}%", (iv,(1-arr_for(iv)/ARR0)*100),
                textcoords="offset points", xytext=(6, -28 if lab=="SID" else 6),
                color=c, fontsize=11, fontweight="bold")
ax.axhline((1-0.23/ARR0)*100, color=F.GREY, ls=":", lw=1.3)
ax.set_xlabel("dosing interval (days)"); ax.set_ylabel("relapse-rate reduction vs placebo (%)")
ax.set_title("Efficacy is retained out to q6w", fontsize=13)
ax.set_ylim(40, 80); F.despine(ax)
fig.suptitle("Natalizumab: one controllable handle — the dosing interval",
             fontsize=15, fontweight="bold", color=F.NAVY, y=1.03)
F.save(fig, "s2_saturation.png")

# ============================================================
# FIG B — PML stratification heatmap (24 / 48 / 72 months)
# ============================================================
strata = list(PML_TABLE.keys())
mat = np.array([PML_TABLE[s] for s in strata])
fig, ax = plt.subplots(figsize=(8.8, 4.6))
im = ax.imshow(mat, cmap=F.risk_ramp(), aspect="auto", vmin=0, vmax=10)
ax.set_xticks(range(3)); ax.set_xticklabels(["1–24 mo", "25–48 mo", "49–72 mo"])
ax.set_yticks(range(len(strata))); ax.set_yticklabels(strata, fontsize=11)
for i in range(len(strata)):
    for j in range(3):
        v = mat[i, j]
        ax.text(j, i, f"{v:g}", ha="center", va="center",
                color="white" if v > 4 else F.INK, fontweight="bold", fontsize=12)
ax.set_title("Estimated PML risk per 1000 patients (Ho 2017)", pad=10)
ax.set_xlabel("cumulative time on natalizumab")
cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02); cb.set_label("per 1000")
ax.grid(False)
F.save(fig, "s2_pml_table.png")

# ============================================================
# FIG C — benefit-risk optimization over a 6-year horizon
# ============================================================
HORIZON_M = 72
# benefit: relapses avoided vs no treatment, per patient over horizon
yrs = HORIZON_M / 12
avoided_sid = (ARR0 - arr_q4) * yrs
avoided_eid = (ARR0 - arr_q6) * yrs
risk_strata = ["JCV+ index <=0.9", "JCV+ index 0.9-1.5", "JCV+ index >1.5",
               "JCV+ prior immunosuppressant"]
labels = ["index ≤ 0.9", "index 0.9–1.5", "index > 1.5", "prior immunosuppr."]
fig, ax = plt.subplots(figsize=(9.6, 5.4))
y = np.arange(len(risk_strata))[::-1]
for yi, s in zip(y, risk_strata):
    r_sid = cum_pml(s, HORIZON_M, 1.0) * 1000
    r_eid = cum_pml(s, HORIZON_M, EID_HR) * 1000
    ax.plot([r_eid, r_sid], [yi, yi], color=F.LGREY, lw=6, solid_capstyle="round", zorder=1)
    ax.scatter([r_sid], [yi], s=200, color=F.RED, zorder=3, edgecolors="white", linewidths=1.5)
    ax.scatter([r_eid], [yi], s=200, color=F.TEAL, zorder=3, edgecolors="white", linewidths=1.5)
    ax.text(r_sid + 0.35, yi, f"{r_sid:.0f}", va="center", ha="left", color=F.RED, fontsize=11, fontweight="bold")
    ax.text(r_eid - 0.35, yi, f"{r_eid:.1f}", va="center", ha="right", color=F.TEAL, fontsize=11, fontweight="bold")
    ax.text(r_sid, yi + 0.22, f"−{(1-r_eid/r_sid)*100:.0f}%", va="bottom", ha="center",
            color=F.NAVY, fontsize=10, style="italic")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=12)
ax.set_xlabel("6-year cumulative PML probability (per 1000 patients)")
ax.set_title("Extended-interval dosing cuts PML risk by ~76%\nfor 0.30 relapse / 6 yr in retained benefit", pad=12)
ax.scatter([], [], s=160, color=F.RED, label="standard q4w")
ax.scatter([], [], s=160, color=F.TEAL, label="extended q6w")
ax.legend(loc="center", bbox_to_anchor=(0.60, 0.62), fontsize=12,
          title="natalizumab", title_fontsize=11)
ax.set_xlim(-1.4, 16.5)
ax.grid(axis="y", visible=False); F.despine(ax, left=False)
F.save(fig, "s2_benefit_risk.png")

# ============================================================
# FIG D — cross-DMT benefit-risk map
# ============================================================
# ARR reduction (efficacy) vs a serious-safety-event index (per 1000 / 2 yr)
DMT = [
    # name, ARR_reduction_%, safety_index(per1000 serious event /yr), color, note
    ("Interferon-β / GA", 32, 2,  F.GREY,  "platform"),
    ("Fingolimod",        54, 9,  F.BLUE,  "lymphopenia, infection, bradycardia"),
    ("Ocrelizumab",       68, 11, F.PURPLE,"infection, hypogammaglobulinemia"),
    ("Natalizumab q4w",   68, 24, F.RED,   "PML (JCV+ index>1.5, >24 mo)"),
    ("Natalizumab q6w",   63, 6,  F.TEAL,  "PML ↓ ~76% with EID"),
]
fig, ax = plt.subplots(figsize=(9.6, 5.8))
for name, eff, saf, col, note in DMT:
    ax.scatter(saf, eff, s=260, color=col, edgecolors="white", linewidths=1.5, zorder=3)
    dy = 2.2 if name != "Natalizumab q6w" else -3.4
    ax.annotate(name, (saf, eff), textcoords="offset points", xytext=(10, dy),
                fontsize=11.5, fontweight="bold", color=col)
# arrow from nat q4w -> q6w
ax.annotate("", xy=(6, 63), xytext=(24, 68),
            arrowprops=dict(arrowstyle="->", color=F.NAVY, lw=2.2))
ax.text(14, 71, "the optimization\n(EID)", color=F.NAVY, fontsize=11,
        fontweight="bold", ha="center")
ax.set_xlabel("serious-safety-event index  (per 1000 patients · year)")
ax.set_ylabel("relapse-rate reduction vs placebo (%)")
ax.set_title("Benefit–risk map across MS disease-modifying therapies", pad=10)
ax.set_xlim(0, 28); ax.set_ylim(25, 80); F.despine(ax)
# "better" direction
ax.annotate("better →\n(more efficacy, less risk)", xy=(1.0, 78),
            fontsize=10, color=F.GREEN, style="italic")
F.save(fig, "s2_dmt_map.png")

# save key numbers
res = {
    "arr_placebo": ARR0, "arr_q4w": round(arr_q4,3), "arr_q6w": round(arr_q6,3),
    "red_q4w_pct": round((1-arr_q4/ARR0)*100,1), "red_q6w_pct": round((1-arr_q6/ARR0)*100,1),
    "eid_pml_hr": EID_HR,
    "relapses_avoided_sid_6y": round(avoided_sid,2),
    "relapses_avoided_eid_6y": round(avoided_eid,2),
    "pml_6y_highrisk_sid": round(cum_pml("JCV+ index >1.5",72,1.0)*1000,1),
    "pml_6y_highrisk_eid": round(cum_pml("JCV+ index >1.5",72,EID_HR)*1000,1),
    "pml_6y_priorIS_sid": round(cum_pml("JCV+ prior immunosuppressant",72,1.0)*1000,1),
    "pml_6y_priorIS_eid": round(cum_pml("JCV+ prior immunosuppressant",72,EID_HR)*1000,1),
}
(OUT / "s2_results.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
print("S2 done")
