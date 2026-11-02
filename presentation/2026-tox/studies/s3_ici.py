"""Study 3 — Oncology: checkpoint-inhibitor efficacy vs immune-related
toxicity as a tunable optimization (metastatic melanoma).

Nivolumab (anti-PD-1) + ipilimumab (anti-CTLA-4) is the most effective and
the most toxic regimen in melanoma. The two antibodies are not interchangeable
handles: CTLA-4 blockade (ipilimumab) drives much of the autoimmune toxicity,
while anti-tumor efficacy is PD-1-dominated once priming has occurred. So the
dose RATIO is a controllable handle on the benefit-risk frontier — exactly what
CheckMate 511 demonstrated (NIVO3/IPI1 halved grade 3-4 toxicity at preserved
efficacy vs NIVO1/IPI3).

A QSP layer links receptor occupancy -> two effector pools (anti-tumor,
self-reactive) with Bliss-type synergy -> ORR and grade 3-4 irAE rate, and a
dynamic tumor + irAE time course. Occupancy-to-outcome maps are calibrated to:

  CheckMate 067 (Larkin NEJM 2015/2019): ORR nivo+ipi 58%, nivo 45%, ipi 19%;
      grade 3-4 treatment-related AE: combo 59%, nivo 23%, ipi 28%.
  CheckMate 511 (Lebbe JCO 2019): NIVO3/IPI1 grade 3-4 34% vs NIVO1/IPI3 48%;
      ORR 46% vs 51% (not significantly different).
  Weber (J Clin Oncol 2017): irAE onset timing by organ system.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.integrate import odeint
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as F
F.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = Path("/home/user/qsp/presentation/2026-tox/data"); OUT.mkdir(exist_ok=True, parents=True)

# ---- receptor occupancy ---------------------------------------------------
# Key biology: ipilimumab's contribution to anti-tumor PRIMING saturates at a
# low dose (efficacy flat from 1 to 3 mg/kg), while its contribution to
# autoimmune toxicity keeps climbing super-linearly with dose. So efficacy and
# toxicity read CTLA-4 occupancy through two different curves.
def nivo_occ(d): return d / (d + 0.3)            # PD-1 saturates at low mg/kg
def ipi_eff(d):  return d / (d + 0.30)           # priming nearly saturated by 1 mg/kg
def ipi_tox(d):  return d**1.8 / (d**1.8 + 2.6)  # toxicity keeps rising super-linearly
def ipi_occ(d):  return ipi_eff(d)               # back-compat alias

# regimens (nivo mg/kg, ipi mg/kg)
REG = {
    "nivolumab 3 (mono)":      (3, 0),
    "ipilimumab 3 (mono)":     (0, 3),
    "NIVO1 / IPI3  (CM511)":   (1, 3),
    "NIVO3 / IPI1  (CM511)":   (3, 1),
}
# calibration targets: ORR %, grade3-4 % — monotherapy arms from CheckMate 067,
# combination arms from the head-to-head CheckMate 511 (one consistent comparison)
TARG = {
    "nivolumab 3 (mono)":    (45, 23),   # CM067
    "ipilimumab 3 (mono)":   (19, 28),   # CM067
    "NIVO1 / IPI3  (CM511)": (51, 48),   # CM511 standard combo
    "NIVO3 / IPI1  (CM511)": (46, 34),   # CM511 reduced-ipilimumab combo
}

names = list(REG)
no = np.array([nivo_occ(REG[n][0]) for n in names])
ie = np.array([ipi_eff(REG[n][1]) for n in names])    # efficacy occupancy
it = np.array([ipi_tox(REG[n][1]) for n in names])    # toxicity occupancy
orr_t = np.array([TARG[n][0] for n in names])
g34_t = np.array([TARG[n][1] for n in names])

# anti-tumor effector score -> ORR via Emax/Hill
def orr_model(p):
    wn, wi, ws, EC = p
    T = wn*no + wi*ie + ws*no*ie
    return 75.0 * T / (EC + T)
fit_o = least_squares(lambda p: orr_model(p) - orr_t, [1.0, 0.4, 1.2, 0.8],
                      bounds=([0,0,0,0.05],[5,5,5,5])).x

# self-reactive effector -> grade3-4 irAE (ipi-tox-weighted) via Emax
G34_EMAX = 95.0
def g34_model(p):
    rn, ri, rs, k = p
    R = rn*no + ri*it + rs*no*it
    return G34_EMAX * R / (k + R)
w_g = np.array([1, 1, 1.6, 1.6])   # weight the two combos (the headline result)
fit_g = least_squares(lambda p: (g34_model(p) - g34_t)*w_g, [0.3, 1.0, 0.8, 0.6],
                      bounds=([0,0,0,0.05],[5,5,5,5])).x

print("ORR fit", np.round(fit_o,3), "pred", np.round(orr_model(fit_o),1), "targ", orr_t)
print("g34 fit", np.round(fit_g,3), "pred", np.round(g34_model(fit_g),1), "targ", g34_t)

def eval_regimen(nd, idose):
    n = nivo_occ(nd)
    wn, wi, ws, EC = fit_o
    rn, ri, rs, k = fit_g
    T = wn*n + wi*ipi_eff(idose) + ws*n*ipi_eff(idose)
    R = rn*n + ri*ipi_tox(idose) + rs*n*ipi_tox(idose)
    return 75.0*T/(EC+T), G34_EMAX*R/(k+R)

# ============================================================
# FIG A — the efficacy / toxicity Pareto frontier over dose ratio
# ============================================================
fig, ax = plt.subplots(figsize=(9.6, 6.0))
# background cloud: grid of regimens
grid = []
for nd in np.linspace(0.3, 3, 16):
    for idose in np.linspace(0, 3, 16):
        o, g = eval_regimen(nd, idose)
        grid.append((g, o, nd, idose))
grid = np.array(grid)
sc = ax.scatter(grid[:,0], grid[:,1], c=grid[:,3], cmap=F.risk_ramp(),
                s=26, alpha=0.6, edgecolors="none", zorder=1)
cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02); cb.set_label("ipilimumab dose (mg/kg)")
# Pareto frontier (max ORR for each toxicity level)
order = grid[grid[:,0].argsort()]
front_x, front_y, best = [], [], -1
for g, o, nd, idose in order:
    if o > best:
        best = o; front_x.append(g); front_y.append(o)
ax.plot(front_x, front_y, color=F.NAVY, lw=2.6, zorder=2, label="efficacy–safety frontier")
# real regimens
cols = {"nivolumab 3 (mono)": F.BLUE, "ipilimumab 3 (mono)": F.GOLD,
        "NIVO1 / IPI3  (CM511)": F.RED, "NIVO3 / IPI1  (CM511)": F.TEAL}
# per-regimen label offset (dx, dy in points) to avoid overlap
LBL = {"nivolumab 3 (mono)": "nivolumab mono",
       "ipilimumab 3 (mono)": "ipilimumab mono",
       "NIVO1 / IPI3  (CM511)": "NIVO1 / IPI3",
       "NIVO3 / IPI1  (CM511)": "NIVO3 / IPI1"}
OFF = {
    "nivolumab 3 (mono)":    (12, 10),
    "ipilimumab 3 (mono)":   (12, 10),
    "NIVO1 / IPI3  (CM511)": (-70, -26),
    "NIVO3 / IPI1  (CM511)": (-150, 12),
}
for n in names:
    o, g = eval_regimen(*REG[n])
    ax.scatter([g], [o], s=240, color=cols[n], edgecolors="white", linewidths=1.6, zorder=4)
    dx, dy = OFF[n]
    ax.annotate(LBL[n], (g, o), textcoords="offset points", xytext=(dx, dy),
                fontsize=11, fontweight="bold", color=cols[n])
# arrow CM067 -> CM511
o1, g1 = eval_regimen(*REG["NIVO1 / IPI3  (CM511)"])
o2, g2 = eval_regimen(*REG["NIVO3 / IPI1  (CM511)"])
ax.annotate("", xy=(g2, o2), xytext=(g1, o1),
            arrowprops=dict(arrowstyle="->", color=F.NAVY, lw=2.4))
ax.text((g1+g2)/2-1, (o1+o2)/2-11, "CheckMate 511: shift ipilimumab down\n−%d pts grade 3–4, efficacy held" % (g1-g2),
        fontsize=10.5, color=F.NAVY, fontweight="bold")
ax.set_xlabel("grade 3–4 immune-related AE rate (%)")
ax.set_ylabel("objective response rate (%)")
ax.set_title("Checkpoint-inhibitor dosing: one frontier, two antibodies", pad=12)
ax.legend(loc="lower right"); F.despine(ax)
F.save(fig, "s3_pareto.png")

# ============================================================
# FIG B — dynamic tumor + irAE time course (combo)
# ============================================================
# states: Ttum (tumor burden, rel), Etum (anti-tumor effector), Eself (self-reactive)
def rhs(y, t, nd, idose):
    Ttum, Etum, Eself = y
    n, i = nivo_occ(nd), ipi_occ(idose)
    drive_t = 0.9*n + 0.4*i + 1.3*n*i
    drive_s = 0.3*n + 1.0*i + 0.8*n*i
    dEtum = 0.20*drive_t - 0.08*Etum
    dEself = 0.12*drive_s - 0.05*Eself      # slower build -> delayed irAE
    r, K, kill = 0.045, 100.0, 0.9
    dT = r*Ttum*(1 - Ttum/K) - kill*Etum*Ttum/(Ttum+20)
    return [dT, dEtum, dEself]
t = np.linspace(0, 52, 400)   # weeks
sol = odeint(rhs, [60, 0, 0], t, args=REG["NIVO3 / IPI1  (CM511)"])
Ttum, Etum, Eself = sol.T
fig, ax = plt.subplots(figsize=(9.4, 5.2))
ax.plot(t, Ttum/60*100, color=F.NAVY, label="tumor burden (% baseline)")
ax.axhline(70, color=F.GREY, ls=":", lw=1.3)
ax.text(0.5, 72, "−30% = partial response", color=F.GREY, fontsize=10)
ax2 = ax.twinx()
ax2.plot(t, Eself, color=F.RED, label="self-reactive effector (irAE driver)")
ax2.plot(t, Etum, color=F.TEAL, ls="--", label="anti-tumor effector")
ax2.set_ylabel("effector pool (relative)", color=F.INK)
ax.set_xlabel("weeks on therapy"); ax.set_ylabel("tumor burden (% baseline)")
ax.set_title("One regimen, two clocks: response builds, toxicity follows", pad=10)
l1, la1 = ax.get_legend_handles_labels(); l2, la2 = ax2.get_legend_handles_labels()
ax.legend(l1+l2, la1+la2, loc="center right", fontsize=10.5)
ax.set_ylim(0, 110); ax2.set_ylim(0, 3.2)
F.despine(ax, left=True, bottom=True); ax2.spines["top"].set_visible(False)
F.save(fig, "s3_timecourse.png")

# ============================================================
# FIG C — irAE onset timing by organ (Weber 2017 pattern)
# ============================================================
organs = [("Skin (rash, pruritus)", 3.5, 2.0, F.GOLD),
          ("GI (colitis, diarrhea)", 6.5, 2.5, F.RED),
          ("Hepatic (transaminitis)", 7.5, 2.5, F.PURPLE),
          ("Pulmonary (pneumonitis)", 8.5, 3.5, F.BLUE),
          ("Endocrine (thyroid, hypophysitis)", 10.5, 3.5, F.TEAL),
          ("Renal (nephritis)", 12.0, 3.5, F.NAVY)]
fig, ax = plt.subplots(figsize=(9.6, 5.0))
y = np.arange(len(organs))[::-1]
for yi, (name, med, sd, col) in zip(y, organs):
    ax.barh(yi, 2*sd, left=med-sd, height=0.5, color=col, alpha=0.75)
    ax.scatter([med], [yi], color=col, s=60, zorder=4, edgecolors="white")
    ax.text(med+sd+0.3, yi, f"~{med:g} wk", va="center", fontsize=10.5, color=F.INK)
ax.set_yticks(y); ax.set_yticklabels([o[0] for o in organs], fontsize=11)
ax.set_xlabel("weeks from start of checkpoint therapy to irAE onset")
ax.set_title("The kinetic signature: irAEs arrive organ by organ", pad=10)
ax.set_xlim(0, 18); ax.grid(axis="y", visible=False); F.despine(ax, left=False)
F.save(fig, "s3_organ_timing.png")

res = {name: {"orr": round(eval_regimen(*REG[name])[0],1),
              "g34": round(eval_regimen(*REG[name])[1],1)} for name in names}
res["cm511_tox_drop_pts"] = round(eval_regimen(1,3)[1]-eval_regimen(3,1)[1],1)
res["cm511_orr_change_pts"] = round(eval_regimen(1,3)[0]-eval_regimen(3,1)[0],1)
(OUT / "s3_results.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
print("S3 done")
