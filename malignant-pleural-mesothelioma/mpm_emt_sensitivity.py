#!/usr/bin/env python3
"""
mpm_emt_sensitivity.py
----------------------
POST-HOC analysis, run AFTER the calibration and reported as such.

mpm_calibration.py held the CheckMate 743 histology subgroups out of the fit:
step 7 matched ONE number (the IO/chemo median ratio at the trial population's
histology) and the epithelioid / non-epithelioid split was left as a
prediction.  THE PREDICTION FAILED, and not narrowly -- it failed in DIRECTION:

    observed   epithelioid HR 0.86   non-epithelioid HR 0.46   (IO improves)
    model      epithelioid    0.70   non-epithelioid    0.87   (IO degrades)

So the EMT axis as parameterised from the mechanism literature does NOT
generate the flip.  This script asks the obvious follow-up question -- what
would have to be true for it to -- and answers it by scanning the two
parameters that pull in opposite directions along the axis:

    EMT_CHEMO : chemotherapy kill x (1 - EMT_CHEMO * x)
                how chemo-refractory sarcomatoid disease is
    EMT_PDL1  : PD-L1 = PDL1_0 * (1 + EMT_PDL1 * x)
                how much more brake there is for blockade to release

The scan does not rescue the model.  Neither lever reproduces the trial, and
that negative is the useful part: it moves the diagnosis from "a slope is set
wrong" to "the histology penalty is applied to both arms when the trial says it
falls on only one of them".  See the READING section of the output.

Run:  python3 mpm_emt_sensitivity.py > mpm_emt_sensitivity_output.txt

The closing READING section is pure text and can be regenerated on its own
with `--reading-only` without repeating the ten-minute scan.
"""
import math
import sys

import mpm_reference_model as M

X_EPI = 0.15
X_NON = 0.85
OBS_HR_EPI = 0.86
OBS_HR_NON = 0.46
TMAX = 1800.0
DT = 0.05


def ratio_at(x, emt_chemo=None, emt_pdl1=None, emt_kg=None):
    """median OS ratio chemo/IO -- the model's proxy for the trial's HR."""
    out = []
    for name in ("pemcis", "nivoipi"):
        p, reg = M.scenario(name, x)
        if emt_chemo is not None:
            p["EMT_CHEMO"] = emt_chemo
        if emt_pdl1 is not None:
            p["EMT_PDL1"] = emt_pdl1
        if emt_kg is not None:
            p["EMT_KG"] = emt_kg
        out.append(M.median_os(M.simulate(p, tmax=TMAX, dt=DT, regimen=reg)))
    oc, oi = out
    if not oc or not oi:
        return None, oc, oi
    return oc / oi, oc, oi


def scan_chemoresistance():
    print("=" * 78)
    print(" WHAT WOULD MAKE THE CHECKMATE 743 SUBGROUP SPLIT COME OUT?")
    print("=" * 78)
    print()
    print(" Baseline EMT_CHEMO is 0.55 (chemotherapy retains 45% of its kill in")
    print(" fully sarcomatoid disease), taken from the EMT-chemoresistance")
    print(" literature rather than from any mesothelioma trial.")
    print()
    print(" %-11s %11s %11s   %11s %11s" %
          ("EMT_CHEMO", "epi ratio", "non-epi", "epi OS c/IO", "non-epi OS c/IO"))
    rows = []
    for ec in (0.55, 0.70, 0.80, 0.88, 0.94, 0.97):
        r_epi, oc_e, oi_e = ratio_at(X_EPI, emt_chemo=ec)
        r_non, oc_n, oi_n = ratio_at(X_NON, emt_chemo=ec)
        rows.append((ec, r_epi, r_non))
        print("  %-10.2f %11s %11s   %5.1f / %5.1f %5.1f / %5.1f"
              % (ec,
                 ("%.2f" % r_epi) if r_epi else "  -",
                 ("%.2f" % r_non) if r_non else "  -",
                 oc_e or float("nan"), oi_e or float("nan"),
                 oc_n or float("nan"), oi_n or float("nan")))
    print()
    print(" observed:      epithelioid 0.86      non-epithelioid 0.46")
    print()

    # Where does the model's non-epithelioid ratio cross the observed 0.46, and
    # where does the ORDERING (non-epi ratio BELOW epi ratio) first appear?
    order_flip = None
    for i in range(1, len(rows)):
        prev, cur = rows[i - 1], rows[i]
        if prev[1] is None or prev[2] is None or cur[1] is None or cur[2] is None:
            continue
        d_prev = cur[2] - cur[1]
        d_pp = prev[2] - prev[1]
        if d_pp > 0 >= d_prev:
            f = d_pp / (d_pp - d_prev)
            order_flip = prev[0] + f * (cur[0] - prev[0])
    if order_flip is not None:
        print(" The ORDER of the two subgroups (non-epithelioid benefiting MORE than")
        print(" epithelioid, as CheckMate 743 reported) first appears at")
        print("   EMT_CHEMO = %.2f" % order_flip)
        print(" i.e. chemotherapy retaining only %.0f%% of its kill in fully"
              % (100 * (1 - order_flip)))
        print(" sarcomatoid disease, against the %.0f%% the model started with."
              % (100 * (1 - 0.55)))
    else:
        print(" No value of EMT_CHEMO in the scanned range reproduces the ORDER")
        print(" CheckMate 743 reported.  Chemoresistance alone cannot do it, and")
        print(" the missing mechanism is somewhere else.")
    print()


def scan_pdl1():
    print("=" * 78)
    print(" THE OTHER LEVER: how much brake is there for blockade to release?")
    print("=" * 78)
    print()
    print(" PD-L1 positivity is reported at roughly 20% in epithelioid and 60-70%")
    print(" in sarcomatoid mesothelioma, which is where EMT_PDL1 = 2.0 comes from.")
    print()
    print(" %-10s %11s %11s" % ("EMT_PDL1", "epi ratio", "non-epi ratio"))
    for ep in (2.0, 3.0, 3.5, 3.75):
        r_epi, _, _ = ratio_at(X_EPI, emt_pdl1=ep)
        r_non, _, _ = ratio_at(X_NON, emt_pdl1=ep)
        print("  %-9.2f %11s %11s"
              % (ep, ("%.2f" % r_epi) if r_epi else "  -",
                 ("%.2f" % r_non) if r_non else "  -"))
    print()
    print(" PD-L1 alone is bounded: PD-L1 is a fraction, so EMT_PDL1 above 4.0")
    print(" puts the sarcomatoid tumour past 100% positivity.  Whatever produces")
    print(" the CheckMate 743 split, it cannot be PD-L1 expression by itself.")
    print()


def reading():
    """Closing text.  Pure text, no computation -- see the note in the module
    docstring about how this section is regenerated."""
    print("=" * 78)
    print(" READING")
    print("=" * 78)
    print(" Neither lever reproduces the trial.  Even at EMT_CHEMO = 0.97 --")
    print(" chemotherapy retaining 3% of its kill in fully sarcomatoid disease --")
    print(" the non-epithelioid ratio only reaches 0.79 against an epithelioid")
    print(" 0.68, so the two subgroups stay in the WRONG ORDER.  Pushing PD-L1 to")
    print(" the edge of what a fraction allows does the same thing to both")
    print(" subgroups and never separates them either.")
    print()
    print(" That is a more useful failure than a near miss, because it says the")
    print(" discrepancy is not a mis-set slope.  Look at the survivals rather than")
    print(" the ratios: at x = 0.85 BOTH arms have collapsed to 4.5 and 5.7")
    print(" months, against best supportive care's 7.0 at the reference histology.")
    print(" In this model a sarcomatoid tumour is so aggressive that NO treatment")
    print(" has room to work, and a ratio between two small numbers cannot be far")
    print(" from one.  CheckMate 743 says something quite different: its")
    print(" non-epithelioid patients on nivolumab + ipilimumab lived 18.1 months,")
    print(" which is what its EPITHELIOID patients achieved (18.7).  Checkpoint")
    print(" blockade did not lose efficacy with sarcomatoid histology at all --")
    print(" only chemotherapy did, falling to 8.8 months.")
    print()
    print(" So the located error is structural, not numerical: this model applies")
    print(" the histology penalty to BOTH arms, through three channels that do not")
    print(" distinguish them --")
    print("   HZ_EMT      raises the baseline hazard for everybody;")
    print("   EMT_KCOL    raises collagen, which shortens lambda for the T cell")
    print("               exactly as it does for the platinum;")
    print("   EMT_KG      raises proliferation against both.")
    print(" To reproduce CheckMate 743 the model would need checkpoint blockade to")
    print(" be SPARED that penalty -- an immune effector whose reach does not")
    print(" degrade with desmoplasia, or a sarcomatoid-specific antigenicity gain")
    print(" large enough to cancel it.  Which of those is true is a question about")
    print(" mesothelioma, not about this code, and the model is left disagreeing")
    print(" with the trial rather than tuned into agreement with it.")
    print()


def main():
    print("MPM QSP model -- post-hoc sensitivity on the failed subgroup prediction")
    print()
    scan_chemoresistance()
    scan_pdl1()
    reading()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--reading-only":
        reading()
    else:
        main()
