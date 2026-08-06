#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Emit gad_mrgsolve_model.R from gad_mrgsolve_template.R.

Two things are substituted so that the R model can never drift out of sync
with the Python reference:

  1. the five FITTED parameters, read from gad_calibration.json;
  2. the $MAIN initial-condition block, which encodes the untreated attractor
     as a cubic in DIS fitted to the reference's own relaxation.
"""
import json

import numpy as np

import gad_python_reference as G

RNAME = {
    "sht": "SHT", "auto": "AUTO", "ne": "NE", "a2auto": "A2AUTO",
    "bdnf": "BDNF", "cpfc": "CPFC", "eamy": "EAMY", "traf": "TRAF",
    "ra2": "RA2", "ra1": "RA1", "depend": "DEPEND",
    "crh": "CRH", "acth": "ACTH", "cort": "CORT", "gr": "GR",
    "sns": "SNS", "auton": "AUTON", "sleepd": "SLEEPD", "worry": "WORRY",
    "expect": "EXPECT", "rnau": "RNAU", "rdizz": "RDIZZ", "rh1": "RH1",
    "sexd": "SEXD", "wt": "WT", "ract": "RACT", "madrs": "MADRSS",
    "cumhaz": "CUMHAZ",
}


def main():
    P = G.default_params()
    G.calibrate_healthy(P)
    grid, Y = G.attractor_grid(P)
    m = (grid >= 0.40) & (grid <= 1.75)
    x = grid[m]

    lines = []
    worst = 0.0
    worst_s = ""
    for s, rn in RNAME.items():
        yv = Y[G.IX[s]][m]
        c = np.polyfit(x, yv, 3)
        pred = np.polyval(c, x)
        rel = float(np.max(np.abs(pred - yv) / np.maximum(np.abs(yv), 1e-6)))
        if rel > worst:
            worst, worst_s = rel, s
        lines.append("%-9s = %+ .10e*d3 %+ .10e*d2 %+ .10e*d1 %+ .10e;"
                     % (rn + "_0", c[0], c[1], c[2], c[3]))
    lines.append("")
    lines.append("// worst relative error of these cubics: %.2e (on %s)" % (worst, worst_s))
    init = "\n".join(lines)

    fit = json.load(open("gad_calibration.json"))
    tpl = open("gad_mrgsolve_template.R").read()
    tpl = tpl.replace("@@INITBLOCK@@", init)
    for k in ("emax_pgb", "emax_a2", "e_ex", "dvisit", "fluct0"):
        tpl = tpl.replace("@@%s@@" % k, "%.6f" % fit[k])
    tpl = tpl.replace("phi_healthy  : 0.98484849",
                      "phi_healthy  : %.8f" % P["phi_healthy"])
    tpl = tpl.replace("sglu_healthy : 1.00000000",
                      "sglu_healthy : %.8f" % P["sglu_healthy"])
    open("gad_mrgsolve_model.R", "w").write(tpl)
    print("wrote gad_mrgsolve_model.R (worst attractor cubic error %.2e on %s)"
          % (worst, worst_s))


if __name__ == "__main__":
    main()
