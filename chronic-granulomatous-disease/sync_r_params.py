#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sync_r_params.py — push calibrated numbers from cgd_calibration.json into
the R files, so that they cannot drift away from the Python implementation
that was actually executed.

Two kinds of number are synced:

  1. The AUTO-SYNC block of cgd_mrgsolve_model.R's $PARAM section — the four
     Hill parameters of K(phi), which are the ENTIRE interface between the
     phagosome kernel and the whole-body model.  If these drift, the R model
     is simulating a different disease from the Python one.

  2. The calibrated constants block of cgd_shiny_app.R.

Nothing here is clever.  It exists because hand-copying seven numbers between
four files is exactly the sort of thing that silently goes wrong, and because
a reviewer should be able to see that the R model's parameters have a
machine-checkable provenance rather than a remembered one.

    python3 sync_r_params.py            # write
    python3 sync_r_params.py --check    # verify only, non-zero exit on drift
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CAL = os.path.join(HERE, "cgd_calibration.json")
RMOD = os.path.join(HERE, "cgd_mrgsolve_model.R")
RSHINY = os.path.join(HERE, "cgd_shiny_app.R")

BEGIN = "// ---- BEGIN AUTO-SYNC (written by sync_r_params.py) ------------------------"
END = "// ---- END AUTO-SYNC ---------------------------------------------------------"

SH_BEGIN = "# ---- BEGIN AUTO-SYNC (written by sync_r_params.py) -----------------------"
SH_END = "# ---- END AUTO-SYNC -------------------------------------------------------"


def build_blocks(cal):
    hp = cal["Kphi_hill"]
    wb = cal["whole_body"]
    dv = cal["derived"]

    r = [BEGIN]
    r.append(f"K0      : {hp['K0']:.5f} : log10 kill at phi=0, from the phagosome kernel")
    r.append(f"Kmax    : {hp['Kmax']:.5f} : log10 kill at phi=1, from the phagosome kernel")
    r.append(f"phi50   : {hp['phi50']:.5f} : Hill midpoint of K(phi)")
    r.append(f"hillK   : {hp['hill']:.5f} : Hill coefficient of K(phi)")
    r.append(END)

    s = [SH_BEGIN]
    s.append("CAL <- list(")
    s.append(f"  K0        = {hp['K0']:.6f},   # log10 kill at phi = 0 (the whole non-oxidative arm)")
    s.append(f"  Kmax      = {hp['Kmax']:.6f},   # log10 kill at phi = 1")
    s.append(f"  phi50     = {hp['phi50']:.6f},   # Hill midpoint: NOT a threshold, it is above 0.5")
    s.append(f"  hillK     = {hp['hill']:.6f},")
    s.append(f"  lambda_b  = {wb['lambda_b']:.6f},   # bacterial exposures / patient-year")
    s.append(f"  muN_b     = {wb['muN_b']:.6f},   # log10 mean bacterial inoculum")
    s.append(f"  lambda_f  = {wb['lambda_f']:.6f},   # fungal exposures / patient-year")
    s.append(f"  muN_f     = {wb['muN_f']:.6f},   # log10 mean conidial inoculum")
    s.append(f"  sdN       = {wb['sdN']:.6f},")
    s.append(f"  CFR_B     = {wb['CFR_B']:.6f},   # case fatality, bacterial (FIXED, not fitted)")
    s.append(f"  CFR_F     = {wb['CFR_F']:.6f},   # case fatality, invasive aspergillosis (FIXED)")
    s.append(f"  Ncrit_cgd = {dv['logNcrit_cgd_null']:.6f},   # log10, untreated X-CGD null")
    s.append(f"  Ncrit_hlt = {dv['logNcrit_healthy']:.6f},   # log10, healthy control")
    s.append(f"  E_TMP     = {dv['E_TMP']:.6f},   # /d bacterial kill at steady-state co-trimoxazole")
    s.append(f"  E_ITZ     = {dv['E_ITZ']:.6f},   # /d hyphal kill at steady-state itraconazole")
    s.append(f"  E_IFN     = {dv['E_IFN']:.6f},   # fold boost of the oxidase-independent arm")
    s.append(f"  dhr_xover = {dv['dhr_crossover']:.6f}    # DHR mean where mosaic/uniform ordering flips")
    s.append(")")
    s.append(SH_END)
    return "\n".join(r), "\n".join(s)


def splice(path, begin, end, block):
    txt = open(path).read()
    if begin not in txt or end not in txt:
        raise SystemExit(f"{os.path.basename(path)}: AUTO-SYNC markers not found")
    i = txt.index(begin)
    j = txt.index(end) + len(end)
    return txt[:i] + block + txt[j:], txt


def main():
    check = "--check" in sys.argv
    cal = json.load(open(CAL))
    rblock, sblock = build_blocks(cal)

    ok = True
    for path, (b, e, blk) in ((RMOD, (BEGIN, END, rblock)),
                              (RSHINY, (SH_BEGIN, SH_END, sblock))):
        if not os.path.exists(path):
            print(f"  SKIP {os.path.basename(path)} (missing)")
            continue
        new, old = splice(path, b, e, blk)
        if new == old:
            print(f"  OK   {os.path.basename(path)} already in sync")
        elif check:
            print(f"  DRIFT {os.path.basename(path)} differs from cgd_calibration.json")
            ok = False
        else:
            open(path, "w").write(new)
            print(f"  WROTE {os.path.basename(path)}")
    if check and not ok:
        sys.exit(1)
    print("PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
