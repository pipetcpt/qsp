"""Merigolix flagship figure — the estrogen-threshold trade-off as a tunable
optimization (the named example in the abstract).

Oral GnRH-receptor antagonist for endometriosis: suppresses LH/FSH -> estradiol
(E2) -> lesion/pain regression (benefit), but too much suppression -> bone-mineral
-density loss and vasomotor symptoms (risk). The 'estrogen threshold' hypothesis
says aim for a partial E2 band (~20-40 pg/mL). A 2-compartment PK + HPG-axis PD
QSP model turns dose selection into an optimization over that band.
Anchors (library Shiny run): E2 trough 2.6 pg/mL; pain 3.1->1.1; lesion 18.5->5.3 mm;
BMD -1.22%; hot flashes 5.4->7.9/day.
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as F
F.apply()
import matplotlib.pyplot as plt

# dose -> E2 trough (Hill suppression of baseline ~ 50 pg/mL follicular mean)
dose = np.linspace(0, 240, 200)   # mg/day (relative scale)
E2_0 = 55.0
E2 = E2_0 * (1 - 0.97*dose/(dose + 45))     # suppression Hill
# benefit (pain, lesion) and risk (BMD loss, hot flashes) as functions of E2
def pain(E2):   return 1.0 + 2.2*(E2/E2_0)**0.7          # lower E2 -> less pain
def lesion(E2): return 4.0 + 15.0*(E2/E2_0)**0.9
def bmd_loss(E2): return 2.6*(1 - E2/E2_0)**1.6          # lower E2 -> more BMD loss (%)
def flashes(E2):  return 4.5 + 7.0*(1 - E2/E2_0)**1.4

fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.7))
# left: dose -> E2 with target band
ax = axes[0]
ax.plot(dose, E2, color=F.NAVY, lw=3)
ax.axhspan(20, 40, color="#E4F1EB", zorder=0)
ax.text(150, 30, "target E2 band\n(~20–40 pg/mL)", color=F.GREEN, fontsize=10.5, fontweight="bold", ha="center")
ax.set_xlabel("merigolix dose (relative)"); ax.set_ylabel("estradiol trough (pg/mL)")
ax.set_title("Dose sets the degree of E2 suppression", fontsize=13)
ax.set_ylim(0, 56); F.despine(ax)
# right: benefit & risk vs E2 (x reversed: more suppression to the right)
ax = axes[1]
E2x = np.linspace(55, 2, 200)
axB = ax
axB.plot(E2x, pain(E2x), color=F.TEAL, lw=2.6, label="pelvic pain (NRS)")
axB.plot(E2x, lesion(E2x)/5, color=F.GREEN, lw=2.6, ls="--", label="lesion size (÷5 mm)")
axR = ax.twinx()
axR.plot(E2x, bmd_loss(E2x), color=F.RED, lw=2.6, label="BMD loss (%)")
axR.plot(E2x, flashes(E2x), color=F.GOLD, lw=2.6, ls="--", label="hot flashes (/day)")
ax.axvspan(20, 40, color="#E4F1EB", zorder=0)
ax.invert_xaxis()
ax.set_xlabel("estradiol trough (pg/mL)   →  more suppression")
ax.set_ylabel("benefit  (pain · lesion)", color=F.INK)
axR.set_ylabel("risk  (BMD loss · hot flashes)", color=F.INK)
ax.set_title("Benefit and risk cross inside the window", fontsize=13)
l1,la1 = axB.get_legend_handles_labels(); l2,la2 = axR.get_legend_handles_labels()
ax.legend(l1+l2, la1+la2, loc="upper center", fontsize=9.5, ncol=2)
ax.set_ylim(0, 5.5); axR.set_ylim(0, 12); axR.spines["top"].set_visible(False)
F.despine(ax, left=True, bottom=True)
fig.suptitle("Merigolix: QSP turns the estrogen-threshold trade-off into an optimization",
             fontsize=14.5, fontweight="bold", color=F.NAVY, y=1.03)
fig.subplots_adjust(wspace=0.42)
F.save(fig, "s6_merigolix.png")
print("S6 merigolix done")
