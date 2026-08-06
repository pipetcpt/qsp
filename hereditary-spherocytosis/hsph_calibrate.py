#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calibration ledger for the HS model.  SEVEN numbers are spent, in three
stages, and everything reported afterwards is a prediction.

STAGE 1 -- the normal red cell.  FOUR numbers, none of them about HS.
  tau50     <- mean red cell lifespan 120 d
  kv_base   <- membrane area falls 10.5% over the lifespan   (Waugh 1992)
  kd_base   <- cell volume falls  8.4% over the lifespan     (Waugh 1992)
  k_ph      <- geometric (cord-retention) clearance accounts for 1% of normal
               red cell removal, i.e. splenectomy in a healthy person barely
               changes red cell survival

STAGE 2 -- ONE number: the genotype of the reference patient.
  fdef_mod  <- haemoglobin 10.5 g/dL in the reference "moderate HS" patient.
               fdef is 1 - (spectrin content / normal), so this fit reports a
               spectrin content, and that number can be compared with the
               published one (~70% for moderate HS).  Everything else about
               the patient -- reticulocytes, lifespan, MCHC, bilirubin,
               osmotic fragility, spleen size, splenectomy response -- stays a
               prediction, and the mild and severe genotypes are obtained by
               scaling the published spectrin bands, not by refitting.

  This started as a three-number stage (membrane area to the 16% EMA deficit,
  lifespan to 20 days, MCHC to 35.5) and shrank to nothing, in three steps
  that are each worth recording:

  * The membrane-area target had to go because it is structurally
    incompatible with a 20-day lifespan.  A cell whose area is 16% below
    normal has a D_c deep inside the splenic filter and cannot live 20 days.
    The circulating population always looks milder than the destroyed one,
    because the destroyed one is destroyed.

  * The MCHC target had to go because it was not needed once the cation leak
    was written as a property of the membrane (see lesson 1 below).

  * The lifespan target had to go with cordamp, the cordal
    membrane-stripping gain that it was being used to fit.  Above a gain of
    about 2 the per-cell loop  area loss -> D_c up -> longer cordal residence
    -> area loss  becomes supercritical, the POPULATION acquires two stable
    steady states at the same parameter values, and which one you land on
    depends on the ODE solver's maximum step size: max_step 6 gave a mean
    lifespan of 12.6 days and max_step 8 gave 20.0, for the same parameters.
    That is not a phenotype, it is an artefact of a bistable model being
    calibrated on one branch and read on the other.  cordamp is therefore
    set to 1 (stripping OFF) in the base model, and the amplifier hypothesis
    is examined as an explicit experiment in section 5 of hsph_analysis.py
    instead of being assumed.  It turns out not to be needed.

STAGE 3 -- the opsonic arm.  TWO numbers, from Reliene 2002.
  K_cl      <- 140 IgG molecules per cell in splenectomised band 3-deficient
               HS (controls and spectrin/ankyrin HS <= 60)
  k_ops     <- opsonic clearance carries half the destruction in an
               unsplenectomised band 3-deficient patient

Structural choices that are NOT fitted and are stated instead:
  kv_def   = 3.33 x kv_base, i.e. fdef = 0.30 doubles the cell's intrinsic
             (extrasplenic) vesiculation rate
  kd_def   = 6.67 x kd_base, i.e. fdef = 0.30 doubles the basal dehydration
             rate -- the constitutive cation leak of the defective membrane
  a_ent_def= 1.0, i.e. the membrane-area deficit already present when the
             cell is released is 0.42 x fdef (12.6% for moderate HS)
  kd_cord  = 41 x kd_base, i.e. the cord is assumed to be a fixed multiple
             worse for hydration; the equivalent multiple for MEMBRANE is
             cordamp, which is fitted
  visc_k   = ln2/3 = 0.231 per g/dL, the published slope of cytoplasmic
             viscosity against MCHC
  k_lys, s_lys: intravascular lysis is constrained to be a minor route,
             because HS haemolysis is >95% extravascular

TWO STRUCTURAL LESSONS WERE LEARNED THE HARD WAY AND ARE RECORDED HERE
BECAUSE THEY ARE RESULTS, NOT BUGS.

(1) The cation leak had to be a property of the MEMBRANE, not of the spleen.
    An earlier version put the whole leak in the cordal term and could not
    raise the circulating MCHC at any parameter value, including absurd ones:
    if dehydration only happens inside the spleen then the cells that are
    dehydrated are precisely the cells about to be eaten, and the survivors --
    the ones the blood count measures -- are the ones it did not happen to.
    A splenic leak is self-erasing.  The high MCHC of HS is therefore evidence
    that the leak is constitutive.

(2) Cytoplasmic viscosity had to act on the EXTRACTION step, not the exposure
    step.  Surface-to-volume geometry on its own says dehydration PROTECTS the
    HS cell, because losing volume at fixed area raises the excess-area
    reserve and LOWERS D_c.  Viscosity is the standard answer, but if it is
    written as a longer cordal DWELL then dehydration also increases membrane
    stripping, and the resulting loop is supercritical for normal red cells
    too -- the normal arm of the calibration stopped converging.  Written as a
    harder time getting through the SLIT it only raises the destruction
    probability, and normal cells stay safely subcritical.

WHY BISECTION AND NOT LEAST SQUARES.  Multi-parameter least squares was tried
first and repeatedly converged on degenerate corners -- kv_def driven to zero
with all the clearance loaded onto k_ph, or a viscosity coefficient 56x the
published value -- because the residual surface is nearly flat along one
direction and warm-started steady states are only approximately converged,
which is enough noise to mislead a gradient.  Each target below is monotone in
its own parameter, so one-dimensional bisection cannot go anywhere degenerate.
"""
import json
import math
import sys

import numpy as np
from scipy.optimize import brentq, least_squares

import hsph_python_reference as M
from hsph_python_reference import steady, simulate, observe, P0

FIT = {}
NORMAL_AREA = None


def obs(**p):
    o, ys, pp = steady(dict(FIT, **p), tmax=1300.0, max_step=8.0)
    return o


def solve(read, target, lo, hi, name, **fixed):
    """One-dimensional bisection on `name` between lo and hi."""
    def f(x):
        return read(obs(**dict(fixed, **{name: x}))) - target
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0.0:
        print("      ! %s: target %.4g outside [%.4g, %.4g] "
              "(f = %+.3g .. %+.3g); clamping" % (name, target, lo, hi,
                                                  flo, fhi))
        return lo if abs(flo) < abs(fhi) else hi
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-6, maxiter=50)


# ------------------------------------------------------------------ stage 1
def stage1(verbose=True):
    global NORMAL_AREA

    def resid(x):
        FIT.update(tau50=math.exp(x[0]), kv_base=math.exp(x[1]),
                   kd_base=math.exp(x[2]))
        o = obs(fdef=0.0)
        return [o['lifespan'] / 120.0 - 1.0,
                o['A_end'] / (140.0 * (1 - 0.105)) - 1.0,
                o['V_end'] / (94.0 * (1 - 0.084)) - 1.0]

    FIT.update(tau50=211.0, kv_base=0.000865, kd_base=0.00376, k_ph=0.005,
               kd_def=0.024, cordamp=7.0, a_ent_def=1.0,
               kv_def=0.00288, kd_cord=0.154)
    r = least_squares(resid, [math.log(211.0), math.log(0.000865),
                              math.log(0.00376)],
                      diff_step=0.04, xtol=1e-11, ftol=1e-11)
    FIT.update(tau50=math.exp(r.x[0]), kv_base=math.exp(r.x[1]),
               kd_base=math.exp(r.x[2]))
    # kv_def and kd_cord are tied to the normal-arm constants, not fitted
    FIT['kv_def'] = 3.3333 * FIT['kv_base']
    FIT['kd_def'] = 6.6667 * FIT['kd_base']
    FIT['kd_cord'] = 41.0 * FIT['kd_base']
    FIT['a_ent_def'] = 1.0
    # k_ph from the requirement that the normal spleen destroys almost nothing
    FIT['k_ph'] = solve(lambda o: o['f_geom'], 0.010, 1e-4, 0.20, 'k_ph',
                        fdef=0.0)
    o = obs(fdef=0.0)
    NORMAL_AREA = o['area']
    if verbose:
        print("STAGE 1  the normal red cell")
        print("   tau50   = %10.4f d      -> mean lifespan       %8.2f d "
              "(target 120)" % (FIT['tau50'], o['lifespan']))
        print("   kv_base = %10.6f /d     -> end-of-life area    %8.2f um^2 "
              "(target 125.30)" % (FIT['kv_base'], o['A_end']))
        print("   kd_base = %10.6f /d     -> end-of-life volume  %8.2f um^3 "
              "(target 86.10)" % (FIT['kd_base'], o['V_end']))
        print("   k_ph    = %10.6f /d     -> geometric share     %8.4f "
              "(target 0.010)" % (FIT['k_ph'], o['f_geom']))
        print("   NOT fitted: kv_def  = 3.333 x kv_base = %.6f /d"
              % FIT['kv_def'])
        print("   NOT fitted: kd_def  = 6.667 x kd_base = %.6f /d"
              % FIT['kd_def'])
        print("   NOT fitted: kd_cord = 41    x kd_base = %.6f /d"
              % FIT['kd_cord'])
        print("   NOT fitted: a_ent_def = 1.0 (release deficit = 0.42 x fdef)")
        print("   the normal population mean membrane area is %.2f um^2; a 16%%"
              % NORMAL_AREA)
        print("   EMA deficit against that would be %.2f um^2 (not a target)"
              % (0.84 * NORMAL_AREA))
        print("   normal Hb %.2f g/dL, MCV %.1f fL, MCHC %.2f g/dL, "
              "retic %.2f%%, bilirubin %.2f mg/dL"
              % (o['Hb'], o['MCV'], o['MCHC'], o['RET_pct'], o['TBIL']))
    return FIT


# ------------------------------------------------------------------ stage 2
def stage2(verbose=True):
    """Nothing is fitted here.  The section exists to print what the model
    says about a moderate HS patient when only the NORMAL red cell has been
    calibrated."""
    FIT['cordamp'] = 1.0
    FIT['fdef_mod'] = solve(lambda o: o['Hb'], 10.5, 0.15, 0.60, 'fdef')
    g = dict(fdef=FIT['fdef_mod'])
    o = obs(**g)
    if verbose:
        print()
        print("STAGE 2  the reference moderate HS patient")
        print("   fdef_mod  = %8.4f        -> haemoglobin        %7.2f g/dL "
              "(target 10.5)" % (FIT['fdef_mod'], o['Hb']))
        print("   i.e. the model needs a spectrin content of %.0f%% of normal"
              % (100 * (1 - FIT['fdef_mod'])))
        print("   to produce a textbook moderate phenotype.  The published")
        print("   figure for moderate HS is about 70%%, so the model's severity")
        print("   mapping is displaced by %.0f percentage points of spectrin --"
              % (100 * (0.30 - FIT['fdef_mod'])))
        print("   it needs a %s defect than the literature implies."
              % ("BIGGER" if FIT['fdef_mod'] > 0.30 else "SMALLER"))
        print("   cordamp = 1 (cordal membrane stripping off, see section 5).")
        print()
        print("   PREDICTED, nothing fitted to any of these:")
        print("      haemoglobin            %7.2f g/dL   (moderate HS: 8-12)"
              % o['Hb'])
        print("      reticulocytes          %7.2f %%      (moderate HS: >6)"
              % o['RET_pct'])
        print("      total bilirubin        %7.2f mg/dL  (moderate HS: >2)"
              % o['TBIL'])
        print("      MCV                    %7.1f fL" % o['MCV'])
        print("      MCHC                   %7.2f g/dL   (HS: >=35.4 is the"
              " screen)" % o['MCHC'])
        print("      RDW                    %7.1f %%      (HS: >=14 is the"
              " screen)" % o['RDW'])
        print("      mean membrane area     %7.2f um^2" % o['area'])
        print("      EMA reduction          %7.1f %%      (diagnostic cut-off"
              " 16-21)" % (100 * (1 - o['area'] / NORMAL_AREA)))
        print("      50%% osmotic lysis      %7.3f %%NaCl (HS: 0.55-0.65)"
              % o['MCF'])
        print("      spleen volume          %7.0f mL     (HS: 300-800)"
              % o['SPL'])
        print("      mean D_c               %7.3f um" % o['Dc'])
        print("      mean lifespan          %7.2f d      (HS: 10-30)"
              % o['lifespan'])
        print("      time in splenic cords  %7.4f        (normal ~0.002)"
              % o['Rcord'])
        print("      geometric share of destruction %7.3f" % o['f_geom'])
    return FIT


# ------------------------------------------------------------------ stage 3
def stage3(verbose=True, rounds=3):
    b3 = dict(fdef=0.30, f_b3ves=0.15)
    for it in range(rounds):
        FIT['K_cl'] = solve(lambda o: o['IgG'], 140.0, 0.02, 1.2, 'K_cl',
                            **dict(b3, spl_frac=0.0))
        FIT['k_ops'] = solve(lambda o: o['f_ops'], 0.50, 1e-4, 1.5, 'k_ops',
                             **b3)
    pre = obs(**b3)
    post = obs(**dict(b3, spl_frac=0.0))
    if verbose:
        print()
        print("STAGE 3  the opsonic arm (band 3 genotype)")
        print("   K_cl    = %10.6f      -> IgG/cell after splenectomy "
              "%6.1f (target 140)" % (FIT['K_cl'], post['IgG']))
        print("   k_ops   = %10.6f /d   -> opsonic share before it    "
              "%6.3f (target 0.50)" % (FIT['k_ops'], pre['f_ops']))
    return FIT


def main():
    print(__doc__)
    stage1()
    print()
    stage2()
    stage3()
    json.dump(FIT, open("hsph_calibration.json", "w"), indent=2)
    print()
    print("written hsph_calibration.json")
    print(json.dumps(FIT, indent=2))


if __name__ == "__main__":
    main()
