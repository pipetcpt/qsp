"""Study 1 figures: the QSP library as a safety atlas."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as F
F.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

D = Path("/home/user/qsp/presentation/2026-tox/data")
models = json.load(open(D / "library_models.json"))
summ = json.load(open(D / "library_summary.json"))
N = 416  # repository headline count

# English labels for therapeutic areas
KO2EN = {
    "종양": "Oncology", "신경": "Neurology", "내분비·대사": "Endocrine/metabolic",
    "소화기·간담도": "GI/hepatobiliary", "혈액": "Hematology", "심혈관": "Cardiovascular",
    "신장·비뇨": "Renal/urologic", "희귀·유전": "Rare/genetic", "감염": "Infectious",
    "호흡기": "Respiratory", "자가면역·류마티스": "Autoimmune/rheumatic",
    "정신·신경": "Psychiatric", "피부": "Dermatologic", "안과": "Ophthalmic",
    "혈관염": "Vasculitis", "근골격·신경근육": "Musculoskeletal/NM",
    "부인·생식": "Gynecologic/repro", "신생아·소아": "Neonatal/pediatric",
    "중독·환경": "Toxicology/environ.", "기타": "Other", "면역결핍": "Immunodeficiency",
    "알레르기": "Allergy",
}

# ========================================================================
# FIG 1 — milestone growth trajectory (documented in KSMB/KBCS decks)
# ========================================================================
fig, ax = plt.subplots(figsize=(9.2, 5.0))
milestones = [("KSMB 2026\n(special session)", 191),
              ("KBCS 2026\n(invited lecture)", 258),
              ("KSOT 2026\n(this talk)", N)]
x = np.arange(len(milestones))
vals = [m[1] for m in milestones]
ax.plot(x, vals, "-o", color=F.NAVY, lw=3.0, ms=12, zorder=3,
        markerfacecolor="white", markeredgewidth=3, markeredgecolor=F.NAVY)
for xi, (lab, v) in zip(x, milestones):
    ax.annotate(f"{v}", (xi, v), textcoords="offset points", xytext=(0, 14),
                ha="center", fontsize=17, fontweight="bold", color=F.NAVY)
ax.set_xticks(x)
ax.set_xticklabels([m[0] for m in milestones], fontsize=12)
ax.set_ylabel("Cumulative disease QSP models")
ax.set_ylim(0, N * 1.18)
ax.set_title("One LLM-built model per day — a library that compounds", pad=14)
# cadence annotation
ax.annotate("+1 model / day,\nautonomous & quality-gated",
            xy=(1.0, 258), xytext=(0.15, 360),
            fontsize=12, color=F.TEAL, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=F.TEAL, lw=1.6))
F.despine(ax)
F.save(fig, "s1_growth.png")

# ========================================================================
# FIG 2 — therapeutic-area coverage with safety-axis fraction
# ========================================================================
tot = summ["total_by_cat"]; saf = summ["safety_by_cat"]
cats = sorted(tot, key=lambda k: tot[k], reverse=True)
cats = [c for c in cats if tot[c] >= 2]
y = np.arange(len(cats))[::-1]
totv = [tot[c] for c in cats]
safv = [saf.get(c, 0) for c in cats]
fig, ax = plt.subplots(figsize=(9.6, 6.6))
ax.barh(y, totv, color=F.LGREY, height=0.72, label="models in area")
ax.barh(y, safv, color=F.TEAL, height=0.72, label="carry an explicit safety / toxicity axis")
for yi, t, sv in zip(y, totv, safv):
    ax.text(t + 0.6, yi, str(t), va="center", ha="left", fontsize=10.5, color=F.INK)
ax.set_yticks(y)
ax.set_yticklabels([KO2EN.get(c, c) for c in cats], fontsize=11)
ax.set_xlabel("Number of disease models")
ax.set_xlim(0, max(totv) * 1.12)
ax.set_title("22 therapeutic areas — three-quarters already carry a safety axis", pad=12)
ax.legend(loc="lower right", fontsize=11)
ax.grid(axis="y", visible=False)
F.despine(ax, left=False)
F.save(fig, "s1_areas.png")

# ========================================================================
# FIG 3 — quality-gate distributions (ODE states, references)
# ========================================================================
odes = np.array([m["odes"] for m in models])
refs = np.array([m["refs"] for m in models])
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.3))
ax = axes[0]
ax.hist(odes, bins=np.arange(10, 86, 4), color=F.BLUE, edgecolor="white")
ax.axvline(15, color=F.RED, lw=2.2, ls="--")
ax.text(15.6, ax.get_ylim()[1]*0.92, "gate ≥ 15", color=F.RED, fontsize=11, fontweight="bold")
ax.set_xlabel("ODE compartments per model"); ax.set_ylabel("models")
ax.set_title(f"Model size (median {int(np.median(odes))})", fontsize=13)
F.despine(ax)
ax = axes[1]
ax.hist(refs, bins=np.arange(0, 230, 12), color=F.TEAL, edgecolor="white")
ax.axvline(30, color=F.RED, lw=2.2, ls="--")
ax.text(33, ax.get_ylim()[1]*0.92, "gate ≥ 30", color=F.RED, fontsize=11, fontweight="bold")
ax.set_xlabel("PubMed references per model"); ax.set_ylabel("models")
ax.set_title(f"Literature grounding (median {int(np.median(refs))})", fontsize=13)
F.despine(ax)
fig.suptitle("Every model clears hard quality gates — size and grounding",
             fontsize=15, fontweight="bold", color=F.NAVY, y=1.03)
F.save(fig, "s1_quality.png")

# ========================================================================
# FIG 4 — the atlas: each model a point (map size vs ODE size), sized by refs
# ========================================================================
top_areas = sorted(tot, key=lambda k: tot[k], reverse=True)[:7]
colmap = {c: F.CYCLE[i] for i, c in enumerate(top_areas)}
fig, ax = plt.subplots(figsize=(9.4, 6.0))
for m in models:
    c = colmap.get(m["category"], F.LGREY)
    ax.scatter(m["odes"], m["nodes"], s=8 + m["refs"]*0.9, c=c,
               alpha=0.62, edgecolors="white", linewidths=0.4, zorder=2)
ax.set_xlabel("ODE compartments  (PK + disease dynamics)")
ax.set_ylabel("Mechanistic-map nodes")
ax.set_title("416 models, one consistent schema — the atlas", pad=12)
handles = [Patch(color=colmap[c], label=KO2EN.get(c, c)) for c in top_areas]
handles.append(Patch(color=F.LGREY, label="other areas"))
ax.legend(handles=handles, loc="upper left", fontsize=10, ncol=2)
ax.set_xlim(10, 86)
F.despine(ax)
# size legend note
ax.text(0.99, 0.02, "marker area ~ reference count", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=10, color=F.GREY, style="italic")
F.save(fig, "s1_atlas.png")

print("S1 figures done")
