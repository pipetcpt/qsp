#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hereditary spherocytosis -- executable reference implementation.

72 ODEs.  This file is the model of record: the mrgsolve (.R) file mirrors it
equation for equation.  Everything here is run, not asserted.

------------------------------------------------------------------ THE IDEA
A red cell is two numbers, membrane area S and volume V, and every clinical
observation in HS is a function of those two numbers (see hsph_geometry.py).
The disease is a LOOP in that pair:

    area loss -> D_c rises -> longer residence in the splenic cords
              -> more area loss

and the loop is closed inside one organ.  Two things follow that are not
obvious from the standard "membrane protein defect -> spherocyte -> spleen
eats it" cartoon:

  1. The spleen is not (only) the SINK, it is the AMPLIFIER.  Splenectomy
     removes both, and the model can separate them by switching off cordal
     membrane-stripping while leaving splenic phagocytosis intact, and vice
     versa.  Only one of the two switches should carry the clinical benefit.

  2. Area loss is only spherocytising if the volume underneath is NOT lost
     with it.  Normal senescence loses area 10.5% and volume 8.4% (Waugh
     1992, rabbit, 50 d) -- "little change in sphericity", in the paper's own
     words.  So what makes an HS vesicle pathogenic is not that there are
     more of them, but WHAT IS IN THEM.  The model carries the vesicle's
     haemoglobin content as an explicit parameter and can be run with HS
     vesiculation rates and normal vesicle contents.

A THIRD mechanism is needed because HS clearance is not all geometric.
Reliene 2002 (Blood 100:2208) found that splenectomy helps spectrin/ankyrin-
deficient cells far more than band 3-deficient cells, and that band
3-deficient cells carry up to 140 IgG molecules per cell while
spectrin/ankyrin-deficient cells carry <=60 (control level).  Their proposed
mechanism is an accounting rule: spectrin/ankyrin-deficient cells shed band 3
IN the vesicles, so band-3 surface density -- and with it the clusters that
low-affinity natural anti-band-3 antibody needs for bivalent binding -- falls
as fast as area does.  Band 3-deficient cells have a relative skeletal excess,
do not shed band 3, and so their band-3 density RISES as area falls.  The
model therefore tracks band 3 copies per cell as a fifth cohort state, and
the genotype enters as ONE number: the band-3 content of the shed vesicle.
Opsonic clearance is shared between spleen and liver; geometric clearance is
splenic only.  That asymmetry is what makes splenectomy genotype-dependent.

------------------------------------------------------------- STATE VECTOR
 9 red cell age cohorts x 5 states               45
   N   cells                       (1e12/L)
   NA  N * membrane area           (1e12 um^2/L)
   NV  N * volume                  (1e12 um^3/L = mL/L)
   NH  N * haemoglobin mass        (1e12 pg/L)
   NB  N * band 3 copies           (1e12 * 1e6 /L)
 erythropoiesis                                    5   EPO PROG ERB RETM RETB
 transfused donor cells                            1   NDON
 spleen                                            3   SPLV CORD MAC
 haem catabolism / hepatobiliary                   7   HPT FHB BILU BILC BILE STONE LDH
 iron                                              4   HEPC FERR FELIV FESPL
 drug PK/PD                                        5   MGUT MCEN MPER ATP DPG
 insults / cofactors                               2   PARVO FOL
                                                  --
                                                  72
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

from hsph_geometry import (d_crit, v_sphere, sphericity, sphericity_index,
                           osmotic_lysis_point, v_solid, C_ISO, NACL_PER_MOSM)

# --------------------------------------------------------------- cohorts
AGE_EDGES = np.array([0.0, 4.0, 9.0, 16.0, 26.0, 42.0, 64.0, 92.0, 130.0, 200.0])
NC = len(AGE_EDGES) - 1                       # 9
DTAU = np.diff(AGE_EDGES)                     # bin widths, days
ADV = 1.0 / DTAU                              # advection rate out of each bin
AGE_MID = 0.5 * (AGE_EDGES[:-1] + AGE_EDGES[1:])

# state layout
iN, iNA, iNV, iNH, iNB = (np.arange(NC) + k * NC for k in range(5))
o = 5 * NC
(iEPO, iPROG, iERB, iRETM, iRETB, iNDON, iSPLV, iCORD, iMAC, iHPT, iFHB,
 iBILU, iBILC, iBILE, iSTONE, iLDH, iHEPC, iFERR, iFELIV, iFESPL, iMGUT,
 iMCEN, iMPER, iATP, iDPG, iPARVO, iFOL) = range(o, o + 27)
NSTATE = o + 27                               # 72

# ------------------------------------------------------------- parameters
P0 = dict(
    # ---- whole body
    BV=5.0,                    # blood volume, L
    # ---- inherited defect (genotype)
    fdef=0.0,                  # vertical-linkage deficit, 0 = normal
    a_ent_def=0.30,            # fraction of the deficit already present at
                               # reticulocyte release (rest is acquired)
    f_b3ves=1.00,              # band 3 content of the shed vesicle, relative
                               # to the parent membrane. 1.0 = spectrin/ankyrin
                               # deficiency (band 3 leaves with the vesicle);
                               # ~0.15 = band 3 deficiency (it does not).
    f_hbves=0.05,              # haemoglobin content of the shed vesicle,
                               # relative to parent MCHC.  HS vesicles are
                               # described as haemoglobin-free.
    # ---- reticulocyte entry geometry
    A0=140.0,                  # um^2 membrane area at release
    V0=94.0,                   # um^3 (mature cell at release; the
                               #       population mean lands near 90)
    H0=30.0,                   # pg haemoglobin
    B0=1.20,                   # 1e6 band 3 copies
    V_RET_MULT=1.22,           # reticulocyte volume excess
    # ---- membrane loss
    kv_base=0.00084,           # /day fractional area loss, normal, no cords
    kv_def=0.0180,             # /day extra per unit fdef
    cordamp=1.0,               # membrane stripping inside the cord, as a
                               #   multiple of the circulating rate.  1.0 = OFF
                               #   in the base model.  See section 5 of
                               #   hsph_analysis.py: switching this on makes
                               #   the per-cell loop supercritical, and above
                               #   about 2 the population acquires TWO stable
                               #   steady states at the same parameters, so
                               #   which one you get depends on the
                               #   integrator's step size.  It is also not
                               #   needed: the base model reproduces the HS
                               #   phenotype without it.
    p_atp=1.5,                 # ATP exponent on vesiculation
    r_ves3=0.025,              # um, vesicle volume/area ratio = radius/3
    # ---- dehydration (K+ loss: KCC + Gardos)
    kd_base=0.00270,           # /day  [CALIBRATED: Waugh volume loss]
    kd_def=0.0,                # /day extra per unit fdef -- a CONSTITUTIVE
                               #   cation leak of the defective membrane,
                               #   present everywhere in the circulation.
                               #   [CALIBRATED to MCHC].  It cannot be moved
                               #   into the cordal term: a leak acquired only
                               #   inside the spleen is self-erasing, because
                               #   the cells that acquire it are exactly the
                               #   cells that are then destroyed, so it never
                               #   appears in the circulating MCHC.
    kd_cord=0.108,             # /day extra at full cordal residence.  NOT
                               #   fitted: set to (cordamp-1)*kd_base, i.e.
                               #   the cord is assumed to be as much worse for
                               #   hydration as it is for membrane.
    MCHC_MAX=41.0,             # g/dL, densest circulating red cell; sets the
                               #   volume floor for the cation leak
    # ---- splenic geometry filter
    f_pass0=72.0,              # splenic passes per day per cell
    spl_flow_exp=0.35,         # flow ~ mass^0.35, not mass^1
    spl_flow_cap=1.60,         # and never more than 1.6x normal
    D50=3.270,                 # um, D_c at which half the splenic passes are
                               #     diverted into the cordal (slow) route
    wD=0.090,                  # um, steepness of that filter  [CALIBRATED]
    tau0=0.00104,              # d (1.5 min) cordal dwell of a normal cell
    Dc_ref=2.856,              # um, D_c of the reticulocyte-fresh normal cell
    w_esc=0.300,               # um, e-folding of dwell with D_c
    visc_k=0.231,              # 1/(g/dL): cytoplasmic viscosity doubles per
                               #   3 g/dL of MCHC above 33 (Chien).  NOT
                               #   fitted.  Without this term the S/V geometry
                               #   alone predicts that dehydration PROTECTS
                               #   the cell, because losing volume at fixed
                               #   area lowers D_c.  It does not: a dense cell
                               #   is slow, and slow is what kills it.
    MCHC_ref=33.0,             # g/dL, reference for the viscosity law
    R_MAX=0.60,                # cap on cordal time fraction
    k_ph=0.0108,               # /day phagocytic hazard while held in the
                               #     cord  [CALIBRATED]
    # ---- opsonic arm
    Bden_ref=1.20 / 140.0,     # band 3 copies per um^2, normal
    IgG_bg=45.0,               # molecules/cell, background (Reliene control)
    IgG_span=180.0,
    K_cl=0.28,                 # band 3 density excess at half-maximal IgG
                               #    [CALIBRATED to Reliene's 140 IgG/cell]
    IgG_thr=70.0,              # molecules/cell below which nothing is eaten
    IgG50=38.0,                # EXCESS over threshold for half-max phagocytosis
    m_ops=3.0,
    k_ops=0.0333,              # /day  [CALIBRATED]
    spl_ops_cap=1.5,           # splenic opsonic capacity scales with mass
                               #     only up to 1.5x
    w_spl_ops=0.72,            # share of opsonic clearance that is splenic
    w_liv_ops=0.28,            # ... and hepatic (survives splenectomy)
    # ---- non-splenic geometric clearance (hepatic sinusoid)
    k_liv_geo=0.035, D_liv=3.45, w_liv=0.16,
    # ---- intravascular lysis of near-spherical cells
    k_lys=0.06, s_lys=0.970,   # intravascular lysis is a MINOR route in HS
                               #   (<5% of destruction); constrained, not fitted
    # ---- senescence (molecular tagging, not geometry)
    k_sen=0.80, tau50=118.0, m_sen=8.0,
    # ---- erythropoiesis
    k_epo=24.0, EPO_a=4.17, EPO_b=0.211, EPO_norm=10.0,
    Emax_mar=8.0, K_epo=48.0,
    k_prog=1 / 3.0, k_erb=1 / 4.0, k_retm=1 / 2.5, tau_ret0=1.0,
    prod0=0.0417,              # 1e12/L/d baseline red cell production
    shift_max=1.6,             # marrow retic release acceleration
    # ---- spleen
    SPL_base=150.0, SPL_max=1200.0, k_spl=1 / 220.0, g_spl=0.146,
    spl_frac=1.0,              # 1 = intact; set by splenectomy scenarios
    k_regrow=0.0,              # /day remnant regrowth (partial splenectomy)
    k_mac=1 / 14.0, a_mac=0.8, b_mac=2.0,
    f_sen_spl=0.30, f_lys_spl=0.0, W0=0.0125,
    # ---- haem catabolism
    BR_PER_G=34.0,             # mg bilirubin per g haemoglobin
    Vd_bil=5.0, Vmax_ugt=12750.0, Km_ugt=30.0, ugt_f=1.0,
    k_bilc=20.0, k_bile=1.0, f_ucb_bile=0.012, k_stone=2.3e-6,
    k_hp_syn=14.0, k_hp_deg=0.14, k_hp_bind=9.0,
    k_fhb=8.0, LDH_per=1.0, k_ldh=2.2, LDH0=170.0,
    # ---- iron
    k_hepc=0.5, HEPC0=1.0, erfe_k=0.45, k_ferr=1 / 30.0,
    fe_per_cell=1.05,          # mg iron per 1e12 cells (0.34 mg/mL RBC)
    k_fe_abs=0.0012, k_fe_liv=1 / 400.0, tx_fe=200.0, k_chel=0.0,
    # ---- mitapivat
    ka_m=1.6, CL_m=4.2, Vc_m=38.0, Q_m=3.0, Vp_m=25.0,
    EC50_m=1150.0, Emax_atp=0.62, Emax_dpg=0.45, k_atp=2.0, k_dpg=1.2,
    dose_m=0.0, tau_dose=0.5,
    # ---- insults / cofactors
    parvo_t=-1.0, parvo_dur=8.0, parvo_supp=0.97,
    fol_ok=1.0, k_fol=0.05,
    # ---- transfusion
    tx_start=-1.0, tx_interval=28.0, tx_units=2.0, tx_cells=2.0,
    k_don=1 / 55.0,
)


# ------------------------------------------------------------ helper maps
def _dc_vec(A, V):
    Vs = A ** 1.5 / (6.0 * math.sqrt(math.pi))
    s = np.clip(V / Vs, 1e-9, 1.0)
    th = np.arccos(-s)
    return 2.0 * np.sqrt(A / math.pi) * np.cos(th / 3.0 - 2.0 * math.pi / 3.0), Vs, s


def _safe(num, den):
    return num / np.maximum(den, 1e-12)


def initial_state(p):
    """Healthy-ish start; the calibration always integrates to steady state
    so the initial condition only has to be in the basin."""
    y = np.zeros(NSTATE)
    n = p['prod0'] * 1.0
    for i in range(NC):
        y[iN[i]] = n * DTAU[i] * math.exp(-AGE_MID[i] / 110.0) / 1.0
    y[iN] *= 5.0 / max(y[iN].sum(), 1e-9)
    A_ent = p['A0'] * (1.0 - p['a_ent_def'] * p['fdef'] * 0.42)
    y[iNA] = y[iN] * A_ent
    y[iNV] = y[iN] * p['V0']
    y[iNH] = y[iN] * p['H0']
    y[iNB] = y[iN] * p['B0'] * (1.0 - p['f_b3ves']
                                * (1.0 - A_ent / p['A0']))
    y[iEPO] = 10.0
    y[iPROG] = p['prod0'] / p['k_prog']
    y[iERB] = p['prod0'] / p['k_erb']
    y[iRETM] = p['prod0'] / p['k_retm']
    y[iRETB] = p['prod0'] * p['tau_ret0']
    y[iSPLV] = p['SPL_base']
    y[iMAC] = 1.0
    y[iHPT] = 100.0
    y[iBILU] = 0.6
    y[iBILC] = 0.15
    y[iLDH] = p['LDH0']
    y[iHEPC] = p['HEPC0']
    y[iFERR] = 90.0
    y[iFELIV] = 0.8
    y[iFESPL] = 300.0
    y[iATP] = 1.0
    y[iDPG] = 1.0
    y[iFOL] = 1.0
    return y


# ------------------------------------------------------------------ RHS
def rhs(t, y, p):
    y = np.maximum(y, 0.0)
    dy = np.zeros_like(y)

    N = y[iN]
    Nsafe = np.maximum(N, 1e-14)
    A = _safe(y[iNA], Nsafe)
    V = _safe(y[iNV], Nsafe)
    H = _safe(y[iNH], Nsafe)
    B = _safe(y[iNB], Nsafe)
    A = np.maximum(A, 20.0)
    V = np.maximum(V, 15.0)

    # ---- geometry
    Dc, Vsph, s = _dc_vec(A, V)

    # ---- spleen: functional mass fraction
    spl_rel = p['spl_frac'] * y[iSPLV] / p['SPL_base']
    # Splenic blood flow does not scale linearly with splenic mass -- flow
    # per gram falls in a big spleen -- and it matters that it does not:
    # with linear scaling the loop  hazard -> workload -> mass -> flow ->
    # hazard  has a gain above one, the spleen runs to its cap, and the
    # steady state stops being unique (the same parameters gave a mean
    # lifespan of 15.5, 18.6 or 30.8 days depending only on the
    # integrator's maximum step).
    f_pass = p['f_pass0'] * min(spl_rel ** p['spl_flow_exp'],
                               p['spl_flow_cap'])

    # ---- cordal filter
    p_slow = 1.0 / (1.0 + np.exp(-(Dc - p['D50']) / p['wD']))
    # Cytoplasmic viscosity enters the EXTRACTION step, not the exposure
    # step.  A dense cell is not held in the cord for longer -- it is
    # held in the SLIT for longer, and that is where it is eaten.
    # Putting viscosity on the dwell time instead couples dehydration
    # back into membrane stripping and makes the loop supercritical
    # for normal red cells as well as HS ones.
    visc = np.exp(np.clip(p['visc_k'] * (100.0 * H / V - p['MCHC_ref']),
                          -3.0, 4.0))
    tau_c = p['tau0'] * np.exp(np.clip((Dc - p['Dc_ref']) / p['w_esc'],
                                       -8, 8))
    Rcord = np.minimum(p['R_MAX'], f_pass * p_slow * tau_c)

    # ---- ATP (mitapivat and cordal metabolic stress)
    atp = max(y[iATP], 0.25)
    atp_f = atp ** (-p['p_atp'])

    # ---- membrane loss and dehydration
    kv = (p['kv_base'] + p['kv_def'] * p['fdef']) \
        * (1.0 + (p['cordamp'] - 1.0) * Rcord) * atp_f
    vA = A * kv                                   # um^2/day
    Vfloor = 100.0 * H / p['MCHC_MAX']
    kd = (p['kd_base'] + p['kd_def'] * p['fdef']
          + p['kd_cord'] * Rcord) / max(atp, 0.3)
    vV = vA * p['r_ves3'] + np.maximum(V - Vfloor, 0.0) * kd
    vH = vA * p['r_ves3'] * (H / V) * p['f_hbves']
    vB = (vA / A) * B * p['f_b3ves']

    # ---- clearance
    #  (1) geometric: retention in the cord, then phagocytosis there
    pi_dest = 1.0 - np.exp(-p['k_ph'] * y[iMAC] * tau_c * visc)
    h_geom = f_pass * p_slow * pi_dest
    #  (2) opsonic: band 3 clustering -> natural anti-band-3 IgG -> FcgammaR
    rho = (B / A) / p['Bden_ref']
    ce = np.maximum(rho - 1.0, 0.0)
    IgG = p['IgG_bg'] + p['IgG_span'] * ce ** 2 / (p['K_cl'] ** 2 + ce ** 2)
    fops = (np.maximum(IgG - p['IgG_thr'], 0.0) / p['IgG50']) ** p['m_ops']
    h_ops = p['k_ops'] * (p['w_spl_ops'] * min(spl_rel, p['spl_ops_cap'])
                          * y[iMAC] + p['w_liv_ops']) \
        * fops / (1.0 + fops)
    #  (3) hepatic geometric (survives splenectomy)
    h_liv = p['k_liv_geo'] / (1.0 + np.exp(-(Dc - p['D_liv']) / p['w_liv']))
    #  (4) intravascular lysis of the nearly-spherical
    h_lys = p['k_lys'] * np.maximum(s - p['s_lys'], 0.0) / (1.0 - p['s_lys'])
    #  (5) senescence: molecular tagging, geometry-independent
    h_sen = p['k_sen'] * (AGE_MID / p['tau50']) ** p['m_sen']

    h = h_geom + h_ops + h_liv + h_lys + h_sen

    # ---- reticulocyte input
    tau_ret = p['tau_ret0'] * (1.0 + 0.9 * min(y[iEPO] / p['EPO_norm'] / 12.0, 1.7))
    prod_in = y[iRETB] / tau_ret                  # 1e12/L/d entering cohort 1
    A_ent = p['A0'] * (1.0 - p['a_ent_def'] * p['fdef'] * 0.42)

    # ---- cohort ODEs (upwind advection of the extensive quantities)
    inflowN = np.empty(NC)
    inflowN[0] = prod_in
    inflowN[1:] = ADV[:-1] * N[:-1]
    dy[iN] = inflowN - ADV * N - h * N
    # Band 3 at release must lose the same fraction the MEMBRANE has already
    # lost by release, weighted by how much band 3 the shed vesicle carries.
    # Getting this wrong (entering every genotype with the full B0 on a
    # already-deficient membrane) gave every genotype an elevated band 3
    # SURFACE DENSITY at birth, so even the spectrin/ankyrin arm -- which by
    # construction should carry control-level IgG -- came out at 107
    # molecules per cell, and the genotype x splenectomy prediction inverted.
    B_ent = p['B0'] * (1.0 - p['f_b3ves'] * (1.0 - A_ent / p['A0']))
    for idx, mean, sink, ent in ((iNA, A, vA, A_ent),
                                 (iNV, V, vV, p['V0']),
                                 (iNH, H, vH, p['H0']),
                                 (iNB, B, vB, B_ent)):
        X = y[idx]
        inflow = np.empty(NC)
        inflow[0] = prod_in * ent
        inflow[1:] = ADV[:-1] * X[:-1]
        dy[idx] = inflow - ADV * X - h * X - N * sink

    # ---- blood counts (reticulocytes count as red cells)
    RBC = N.sum() + y[iRETB] + y[iNDON]
    HBmass = (y[iNH].sum() + y[iRETB] * p['H0'] + y[iNDON] * 30.0)   # 1e12 pg/L
    Hb = HBmass / 10.0                                               # g/dL
    Vol = (y[iNV].sum() + y[iRETB] * p['V0'] * p['V_RET_MULT']
           + y[iNDON] * 88.0)
    Hct = Vol / 1000.0

    # ---- erythropoiesis
    EPO_t = 10.0 ** (p['EPO_a'] - p['EPO_b'] * max(Hb, 2.0))
    dy[iEPO] = p['k_epo'] * (EPO_t - y[iEPO])
    e = max(y[iEPO] - p['EPO_norm'], 0.0)
    stim = 1.0 + (p['Emax_mar'] - 1.0) * e / (p['K_epo'] + e)
    parvo = 1.0 - p['parvo_supp'] * y[iPARVO]
    folf = y[iFOL] / (0.35 + y[iFOL]) / (1.0 / 1.35)
    fef = 1.0
    dy[iPROG] = p['prod0'] * stim * max(parvo, 0.0) * min(folf, 1.0) * fef \
        - p['k_prog'] * y[iPROG]
    dy[iERB] = p['k_prog'] * y[iPROG] - p['k_erb'] * y[iERB]
    shift = 1.0 + (p['shift_max'] - 1.0) * e / (p['K_epo'] + e)
    dy[iRETM] = p['k_erb'] * y[iERB] - p['k_retm'] * shift * y[iRETM]
    dy[iRETB] = p['k_retm'] * shift * y[iRETM] - prod_in
    dy[iFOL] = p['k_fol'] * (p['fol_ok'] - y[iFOL]) \
        - 0.010 * max(stim - 1.0, 0.0) * y[iFOL]

    # ---- transfusion (donor cells: normal geometry, no defect)
    dy[iNDON] = -p['k_don'] * y[iNDON]
    if p['tx_start'] >= 0.0 and t >= p['tx_start']:
        ph = (t - p['tx_start']) % p['tx_interval']
        if ph < p['tau_dose']:
            dy[iNDON] += p['tx_units'] * p['tx_cells'] / p['tau_dose'] / p['BV']

    # ---- destruction fluxes
    destN = (h * N).sum() + ADV[-1] * N[-1]                # 1e12 cells/L/d
    destHb = (h * y[iNH]).sum() + ADV[-1] * y[iNH][-1]     # 1e12 pg/L/d = g/L/d
    destHb += p['k_don'] * y[iNDON] * 30.0
    destHb_g = destHb * p['BV']                            # g/day
    lysHb_g = ((h_lys * y[iNH]).sum()) * p['BV']

    # ---- spleen mass and macrophage activation
    W = ((h_geom + p['f_sen_spl'] * h_sen + p['f_lys_spl'] * h_lys) * N).sum() \
        + (p['w_spl_ops'] * min(spl_rel, 3.0)
           / (p['w_spl_ops'] + p['w_liv_ops'])) * (h_ops * N).sum()
    W0 = p['W0']
    load = W / W0
    SPL_t = min(p['SPL_base'] * (1.0 + p['g_spl'] * max(load - 1.0, 0.0)),
                p['SPL_max'])
    dy[iSPLV] = p['k_spl'] * (SPL_t - y[iSPLV]) * p['spl_frac'] \
        + p['k_regrow'] * y[iSPLV] * max(1.0 - y[iSPLV] / max(SPL_t, 1.0), 0.0)
    dy[iCORD] = (Rcord * N).sum() - y[iCORD]        # fast: cord pool tracks
    dy[iMAC] = p['k_mac'] * (1.0 + p['a_mac'] * max(load - 1.0, 0.0)
                             / (1.0 + max(load - 1.0, 0.0) / p['b_mac'])
                             - y[iMAC])

    # ---- haem catabolism / bilirubin
    BR_prod = p['BR_PER_G'] * destHb_g                     # mg/day
    dy[iFHB] = lysHb_g * 1000.0 / p['BV'] / 10.0 - p['k_fhb'] * y[iFHB]
    dy[iHPT] = p['k_hp_syn'] - p['k_hp_deg'] * y[iHPT] \
        - p['k_hp_bind'] * y[iFHB] * y[iHPT] / (10.0 + y[iHPT])
    conj = p['Vmax_ugt'] * p['ugt_f'] * y[iBILU] / (p['Km_ugt'] + y[iBILU])
    dy[iBILU] = (BR_prod - conj) / (10.0 * p['Vd_bil'])
    dy[iBILC] = conj / (10.0 * p['Vd_bil']) - p['k_bilc'] * y[iBILC]
    dy[iBILE] = p['k_bilc'] * y[iBILC] * 10.0 * p['Vd_bil'] - p['k_bile'] * y[iBILE]
    ucb_bile = p['f_ucb_bile'] * y[iBILE] * (1.0 + 0.5 * (y[iBILU] / 0.6 - 1.0))
    dy[iSTONE] = p['k_stone'] * max(ucb_bile, 0.0) * (1.0 - min(y[iSTONE], 1.0))
    dy[iLDH] = p['k_ldh'] * (p['LDH0'] + 210.0 * (destN / (5.0 / 120.0) - 1.0)
                             + 900.0 * lysHb_g - y[iLDH])

    # ---- iron
    erfe = max(stim - 1.0, 0.0)
    dy[iHEPC] = p['k_hepc'] * (p['HEPC0'] * (1.0 + 0.9 * (y[iFELIV] / 0.8 - 1.0))
                               / (1.0 + p['erfe_k'] * erfe) - y[iHEPC])
    abs_fe = p['k_fe_abs'] * 1000.0 / max(y[iHEPC], 0.05)          # mg/day
    tx_fe = 0.0
    if p['tx_start'] >= 0.0 and t >= p['tx_start']:
        tx_fe = p['tx_units'] * p['tx_fe'] / p['tx_interval']
    dy[iFELIV] = (abs_fe + 0.55 * tx_fe - p['k_chel'] * y[iFELIV] * 12.0
                  - 1.0) * p['k_fe_liv']
    dy[iFESPL] = 0.45 * tx_fe + 0.02 * destN * p['fe_per_cell'] * p['BV'] \
        - 0.02 * y[iFESPL] - p['k_chel'] * y[iFESPL] * 0.6
    dy[iFERR] = p['k_ferr'] * (30.0 + 120.0 * y[iFELIV] / 0.8 - y[iFERR])

    # ---- mitapivat PK/PD
    dose = 0.0
    if p['dose_m'] > 0.0:
        ph = t % 0.5
        if ph < 0.02:
            dose = p['dose_m'] / 0.02
    dy[iMGUT] = dose - p['ka_m'] * 24.0 * y[iMGUT]
    conc = y[iMCEN] / p['Vc_m'] * 1000.0                      # ng/mL
    dy[iMCEN] = p['ka_m'] * 24.0 * y[iMGUT] \
        - p['CL_m'] * 24.0 * y[iMCEN] / p['Vc_m'] \
        - p['Q_m'] * 24.0 * (y[iMCEN] / p['Vc_m'] - y[iMPER] / p['Vp_m'])
    dy[iMPER] = p['Q_m'] * 24.0 * (y[iMCEN] / p['Vc_m'] - y[iMPER] / p['Vp_m'])
    drv = conc / (p['EC50_m'] + conc)
    dy[iATP] = p['k_atp'] * (1.0 + p['Emax_atp'] * drv
                             - 0.10 * min((Rcord * N).sum() / max(RBC, 1e-6), 1.0)
                             - y[iATP])
    dy[iDPG] = p['k_dpg'] * (1.0 - p['Emax_dpg'] * drv - y[iDPG])

    # ---- parvovirus B19
    if p['parvo_t'] >= 0.0:
        on = 1.0 if (p['parvo_t'] <= t < p['parvo_t'] + p['parvo_dur']) else 0.0
        dy[iPARVO] = 3.0 * (on - y[iPARVO])
    else:
        dy[iPARVO] = -3.0 * y[iPARVO]

    return dy


# ------------------------------------------------------------- observables
def observe(t, y, p):
    y = np.maximum(y, 0.0)
    N = y[iN]
    Ns = np.maximum(N, 1e-14)
    A = np.maximum(_safe(y[iNA], Ns), 20.0)
    V = np.maximum(_safe(y[iNV], Ns), 15.0)
    H = _safe(y[iNH], Ns)
    B = _safe(y[iNB], Ns)
    Dc, Vsph, s = _dc_vec(A, V)

    RET = y[iRETB]
    RBC = N.sum() + RET + y[iNDON]
    Hb = (y[iNH].sum() + RET * p['H0'] + y[iNDON] * 30.0) / 10.0
    Vol = (y[iNV].sum() + RET * p['V0'] * p['V_RET_MULT'] + y[iNDON] * 88.0)
    Hct = Vol / 1000.0
    MCV = Vol / max(RBC, 1e-9)
    MCH = (y[iNH].sum() + RET * p['H0'] + y[iNDON] * 30.0) / max(RBC, 1e-9)
    MCHC = 100.0 * MCH / max(MCV, 1e-9)

    # RDW: CV of the volume distribution over cohorts + retics + donor
    vv = np.concatenate([V, [p['V0'] * p['V_RET_MULT']], [88.0]])
    ww = np.concatenate([N, [RET], [y[iNDON]]])
    ww = ww / max(ww.sum(), 1e-12)
    mu = (ww * vv).sum()
    sd = math.sqrt(max((ww * (vv - mu) ** 2).sum(), 0.0))
    RDW = 100.0 * sd / max(mu, 1e-9)

    spl_rel = p['spl_frac'] * y[iSPLV] / p['SPL_base']
    # Splenic blood flow does not scale linearly with splenic mass -- flow
    # per gram falls in a big spleen -- and it matters that it does not:
    # with linear scaling the loop  hazard -> workload -> mass -> flow ->
    # hazard  has a gain above one, the spleen runs to its cap, and the
    # steady state stops being unique (the same parameters gave a mean
    # lifespan of 15.5, 18.6 or 30.8 days depending only on the
    # integrator's maximum step).
    f_pass = p['f_pass0'] * min(spl_rel ** p['spl_flow_exp'],
                               p['spl_flow_cap'])
    p_slow = 1.0 / (1.0 + np.exp(-(Dc - p['D50']) / p['wD']))
    # Cytoplasmic viscosity enters the EXTRACTION step, not the exposure
    # step.  A dense cell is not held in the cord for longer -- it is
    # held in the SLIT for longer, and that is where it is eaten.
    # Putting viscosity on the dwell time instead couples dehydration
    # back into membrane stripping and makes the loop supercritical
    # for normal red cells as well as HS ones.
    visc = np.exp(np.clip(p['visc_k'] * (100.0 * H / V - p['MCHC_ref']),
                          -3.0, 4.0))
    tau_c = p['tau0'] * np.exp(np.clip((Dc - p['Dc_ref']) / p['w_esc'],
                                       -8, 8))
    Rcord = np.minimum(p['R_MAX'], f_pass * p_slow * tau_c)
    pi_dest = 1.0 - np.exp(-p['k_ph'] * y[iMAC] * tau_c * visc)
    h_geom = f_pass * p_slow * pi_dest
    rho = (B / A) / p['Bden_ref']
    ce = np.maximum(rho - 1.0, 0.0)
    IgG = p['IgG_bg'] + p['IgG_span'] * ce ** 2 / (p['K_cl'] ** 2 + ce ** 2)
    fops = (np.maximum(IgG - p['IgG_thr'], 0.0) / p['IgG50']) ** p['m_ops']
    h_ops = p['k_ops'] * (p['w_spl_ops'] * min(spl_rel, p['spl_ops_cap'])
                          * y[iMAC] + p['w_liv_ops']) \
        * fops / (1.0 + fops)
    h_liv = p['k_liv_geo'] / (1.0 + np.exp(-(Dc - p['D_liv']) / p['w_liv']))
    h_lys = p['k_lys'] * np.maximum(s - p['s_lys'], 0.0) / (1.0 - p['s_lys'])
    h_sen = p['k_sen'] * (AGE_MID / p['tau50']) ** p['m_sen']
    h = h_geom + h_ops + h_liv + h_lys + h_sen

    destN = (h * N).sum() + ADV[-1] * N[-1]
    lifespan = N.sum() / max(destN, 1e-12)
    w = N / max(N.sum(), 1e-12)
    # geometry of the cells at the moment they are destroyed -- the model's
    # counterpart of a density-fractionated "senescent" population
    wd = h * N
    wd[-1] += ADV[-1] * N[-1]
    wd = wd / max(wd.sum(), 1e-12)
    A_end = float((wd * A).sum())
    V_end = float((wd * V).sum())
    age_death = float((wd * AGE_MID).sum())

    # osmotic fragility of the circulating population
    mcf = np.array([osmotic_lysis_point(A[i], V[i], 100.0 * H[i] / V[i])[1]
                    for i in range(NC)])
    mcf_pop = float((w * mcf).sum())

    return dict(
        t=t, Hb=Hb, Hct=100 * Hct, RBC=RBC, MCV=MCV, MCH=MCH, MCHC=MCHC,
        RDW=RDW, RET_pct=100.0 * RET / max(RBC, 1e-9), RET_abs=1000.0 * RET,
        EPO=y[iEPO], lifespan=lifespan, destN=destN,
        Dc=float((w * Dc).sum()), sph=float((w * s).sum()),
        SI=float((w * (36 * math.pi * V ** 2) ** (1 / 3.) / A).sum()),
        area=float((w * A).sum()), vol=float((w * V).sum()),
        EMA=float((w * A).sum()) / p['A0'],
        MCF=mcf_pop,
        IgG=float((w * IgG).sum()), A_end=A_end, V_end=V_end,
        age_death=age_death,
        Rcord=float((w * Rcord).sum()),
        f_geom=float((h_geom * N).sum() / max(destN, 1e-12)),
        f_ops=float((h_ops * N).sum() / max(destN, 1e-12)),
        f_sen=float((h_sen * N).sum() / max(destN, 1e-12)),
        f_liv=float((h_liv * N).sum() / max(destN, 1e-12)),
        f_lys=float((h_lys * N).sum() / max(destN, 1e-12)),
        TBIL=y[iBILU] + y[iBILC], IBIL=y[iBILU], DBIL=y[iBILC],
        LDH=y[iLDH], HPT=y[iHPT], SPL=p['spl_frac'] * y[iSPLV],
        MAC=y[iMAC],
        STONE=100.0 * y[iSTONE], FERR=y[iFERR], FELIV=y[iFELIV],
        ATP=y[iATP], DPG=y[iDPG], MITA=y[iMCEN] / p['Vc_m'] * 1000.0,
        BRprod=p['BR_PER_G'] * (((h * y[iNH]).sum() + ADV[-1] * y[iNH][-1])
                                * p['BV']),
        NDON=y[iNDON], FOL=y[iFOL],
    )


# ---------------------------------------------------------------- drivers
def simulate(p, tmax=400.0, y0=None, n=401, rtol=1e-7, atol=1e-9,
             max_step=1.0):
    pp = dict(P0)
    pp.update(p)
    if y0 is None:
        y0 = initial_state(pp)
    te = np.linspace(0.0, tmax, n)
    sol = solve_ivp(rhs, (0.0, tmax), y0, args=(pp,), method='LSODA',
                    t_eval=te, rtol=rtol, atol=atol, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, pp


def steady(p, tmax=900.0, max_step=6.0):
    sol, pp = simulate(p, tmax=tmax, n=61, max_step=max_step)
    return observe(sol.t[-1], sol.y[:, -1], pp), sol.y[:, -1], pp
