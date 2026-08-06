#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mkmrgsolve.py -- write om_mrgsolve_model.R from the template plus the values
that om_calibrate.py actually fitted.

The point of generating rather than hand-editing is that the R model and the
Python reference cannot drift apart: the fitted numbers live in exactly one
place (om_calibration.json) and both consumers read them from there.

It also handles the ONE unit correction that the two implementations treat
differently.  The Python driver renormalises the barrier to D = 1 at runtime
by scaling aS; mrgsolve has no such step, so the corrected value has to be
baked in.  At the drug-free attractor

    S = S0,  P1 = P2 = P3 = aS.S0/k_p,  D = k_p.P3/k_shed = aS.S0/k_shed

so D = 1 requires exactly  aS = k_shed / S0.  Any other value would move
D_crit off the steady state it is defined as a fraction of, and the two
implementations would silently disagree about what "ulcerated" means.
"""
import json

CAL = json.load(open("om_calibration.json"))
F = CAL["fitted"]

import om_python_reference as M

k_shed = M.P["k_shed"]
S0 = M.P["S0"]
aS_corr = k_shed / S0

# The 5-FU potency IS fitted (stage 4: bolus 5-FU d1-5 -> 5 d of ulcerative
# stomatitis).  It is read from the calibration rather than from the module
# default, which is only a starting value.
pot_5fu = CAL["pot_5fu"]

subs = {
    "__LAMS__": "%.6g" % F["lamS"],
    "__GREG__": "%.6g" % F["greg"],
    "__POTMEL__": "%.6g" % F["pot_mel"],
    "__RADPOT__": "%.6g" % F["rad_pot"],
    "__POTCIS__": "%.6g" % CAL["pot_cis"],
    "__CYEQ__": "%.6g" % CAL["cy_equiv"],
    "__AS__": "%.6g" % aS_corr,
}

src = open("om_mrgsolve_template.R").read()
for k, v in subs.items():
    src = src.replace(k, v)
if "pot_5fu: 12000" not in src:
    raise SystemExit("pot_5fu placeholder not found in the template")
src = src.replace("pot_5fu: 12000",
                  "pot_5fu: %-8.6g" % pot_5fu)

missing = [k for k in subs if k in src]
if missing:
    raise SystemExit("unsubstituted placeholders: %s" % missing)

open("om_mrgsolve_model.R", "w").write(src)
print("wrote om_mrgsolve_model.R")
for k, v in subs.items():
    print("  %-12s %s" % (k, v))
