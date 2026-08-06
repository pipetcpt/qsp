#!/usr/bin/env python3
# =====================================================================
#  csdh_verify.py -- independent Python/scipy re-implementation of the
#  chronic subdural haematoma (cSDH) QSP model.
#
#  PURPOSE.  The mrgsolve file csdh_mrgsolve_model.R is the deliverable;
#  this file exists to CHECK it.  Every ODE here was written from the
#  same equation sheet but typed independently, integrated with LSODA,
#  and compared against the calibration targets below.  Defects found
#  this way are listed in README.md.
#
#  THE MODEL IN ONE PARAGRAPH.  A chronic subdural haematoma is not a
#  clot that is failing to resorb.  It is a SECRETING ORGAN: a
#  neomembrane whose immature capillaries filter plasma into the
#  subdural space faster than meningeal lymphatics can take it away.
#  So the volume is a FIXED POINT of
#         dV/dt = J_exudation + J_rebleed - J_absorption
#  and every therapy is a statement about which of those three terms it
#  touches.  Surgery empties the STOCK and leaves the membrane.
#  Embolisation starves the membrane and leaves the stock.  That single
#  distinction is what the 2024-25 randomised trials are arguing about.
#
#  Units: time = days, volume = mL, drug amounts = mg, pressures = mmHg.
# =====================================================================

import numpy as np
from scipy.integrate import solve_ivp


def sp(x, s=0.7):
    """Softplus: a smooth max(0, x).  The hard kinks in the pressure
    terms made the system stiff enough that LSODA stalled indefinitely;
    smoothing them is the fix, not loosening the tolerances."""
    return s * np.log1p(np.exp(np.clip(x / s, -60.0, 60.0)))

# ---------------------------------------------------------------------
# 1.  PARAMETERS
# ---------------------------------------------------------------------
P = dict(

    # ---- patient descriptors -------------------------------------
    AGE      = 76.0,     # years
    WT       = 70.0,     # kg
    # V_RES = intracranial reserve (atrophy) volume in mL.  THE key
    # patient covariate: it appears exactly twice in the model, with
    # opposite clinical signs (see section 8).
    V_RES    = 30.0,

    # ---- geometry -------------------------------------------------
    # The neomembrane has an EXTENT of its own.  Draining the fluid does
    # not un-build the membrane, so area must be a state, not a function
    # of the current volume -- otherwise a thin post-operative film has
    # almost no exchange area and no drained cavity can ever refill.
    A_MAX    = 180.0,    # cm2, area of one cerebral convexity
    K_A      = 22.0,     # mL, volume at which the collection has spread halfway
    K_SPREAD = 0.35,     # /d spreading rate
    # Apposition is a PERMANENT healing event, not a pause.  Sustained
    # collapse lets the two membrane leaves FUSE, and that is the only way a
    # 48 h drain can buy a durable benefit rather than a ~10 day delay.
    K_FUSEX  = 0.022,    # /d loss of membrane extent once layers appose
    D_CAP    = 35.0,     # mm, max plausible thickness
    D_MEANF  = 0.45,     # mean/max thickness ratio for a crescent
    P_MIN    = 1.5,      # mmHg, floor on intracranial pressure
    V_FLOOR  = 0.5,      # mL, floor on cavity volume in washout terms

    # ---- Starling exchange across the outer neomembrane -----------
    P_CAP    = 25.0,     # membrane capillary hydrostatic pressure, mmHg
    PI_PL    = 25.0,     # plasma colloid osmotic pressure, mmHg
    PI_COEF  = 3.6,      # mmHg per g/dL of cavity protein
    LP_IMM   = 0.240641,     # mL/(d * 100cm2 * mmHg), immature neovessel
    LP_MAT   = 0.06,     # ... mature vessel: ~12x tighter
    SIG_IMM  = 0.25,     # reflection coefficient, immature (fenestrated)
    SIG_MAT  = 0.85,     # ... mature
    PROT_PL  = 0.070,    # plasma protein, g/mL
    PS_PROT  = 3.0,      # mL/d per 100cm2, diffusive protein equilibration

    # ---- clearance of the collection ------------------------------
    K_LYMPH  = 0.981762,     # mL/d per 100 cm2 (meningeal lymphatic route)
    K_ARACH  = 0.085,    # /d (bulk, arachnoid granulation route)
    SEPT_ABS = 0.55,     # septation penalty on absorption
    P_SINK   = 1.0,      # mmHg, dural lymphatic / deep cervical node sink
    P_REF    = 10.0,     # mmHg, normalisation of the absorption gradient
    K_DRAIN  = 3.0,      # /d active subdural drain
    P_DRAIN  = 2.0,      # mmHg, height the drain is set at

    # ---- rebleed from the fragile membrane ------------------------
    K_RB     = 3.70,     # mL/(d * 100cm2) per unit N_CAP at full fragility
    MMA_FLOOR= 0.371842,     # fraction of membrane supply NOT from the MMA
                         # (ophthalmic / occipital / falcine collaterals, and
                         #  incomplete embolisation) -- at 0.15 embolisation
                         #  was absolutely curative in every patient
    W_IL8    = 0.15,     # autocrine IL-8 amplification of macrophage drive
    A_MMP    = 0.60,     # MMP-9 contribution to fragility
    KM_MMP   = 220.0,    # ng/mL
    A_ANG    = 0.35,     # Ang-2 contribution to fragility
    KM_ANG   = 12.0,     # ng/mL
    K_RECAN  = 0.004,    # /d MMA recanalisation / collateral reconstitution

    # ---- haemoglobin / density ------------------------------------
    HB_BLOOD = 150.0,    # mg/mL whole blood haemoglobin
    K_HBDEG  = 0.080,    # /d haem breakdown & clearance

    # ---- fibrinolysis loop ----------------------------------------
    K_TPA    = 6500.0,   # ng/d per 100cm2 per unit N_CAP
    KDEG_TPA = 3.0,      # /d
    K_PAI    = 0.055,    # /d per (ng/mL) PAI-1
    PAI0     = 22.0,     # ng/mL basal PAI-1 in cavity
    K_PAISYN = 2.2,      # ng/mL/d
    K_PAIDEG = 0.10,     # /d
    K_PLS    = 0.052,    # plasmin activity per (ng/mL tPA) per day
    KEL_PLS  = 1.4,      # /d
    IC50_TXA = 5.0,      # mg/L, TXA inhibition of plasminogen activation
    K_FIBF   = 3.2,      # mg fibrin per mL of extravasated blood
    K_LYS    = 0.55,     # /d per unit plasmin
    K_FDPY   = 1.0,      # mg FDP per mg fibrin lysed
    KDEG_FDP = 0.12,     # /d
    K_FDP    = 150.0,    # ug/mL, FDP conc halving haemostatic competence
    IC50_DOAC= 110.0,    # ng/mL apixaban in cavity

    # ---- angiogenesis / membrane biology --------------------------
    K_SPROUT = 0.410,    # /d
    N_MAX    = 2.50,     # carrying capacity of vessel density
    # Maturation is not a fixed rate.  An ACTIVE haematoma keeps its own
    # vessels immature (VEGF/Ang-2 antagonise pericyte stabilisation), so
    # the membrane stays leaky as long as the haem switch is lit.  Quench
    # the switch and the same vessels mature, sigma rises, Lp falls and
    # the exudation collapses.  This is also exactly where atorvastatin
    # acts, so its effect and spontaneous healing share one term.
    K_MAT    = 0.0600,   # /d intrinsic immature -> mature conversion
    MAT_SUPP = 0.920,    # fraction of maturation suppressed by active disease
    K_MMPDES = 0.100,    # /d MMP-9-driven loss of MATURE vessels (remodelling)
    KM_MMPD  = 220.0,    # ng/mL
    KD_CAP   = 0.2000,   # /d immature vessel regression on VEGF withdrawal
    KD_MAT   = 0.0200,   # /d mature vessel loss
    EMAX_ATV = 6.0,      # atorvastatin fold-increase of maturation
    EC50_ATV = 1.6,      # ug/L effect-site atorvastatin
    K_VEGF   = 480.0,    # pg/mL/d per unit N_MAC
    KEL_VEGF = 4.0,      # /d
    VEGF_REF = 120.0,    # pg/mL normalisation
    EMAX_DEXV= 2.30,     # dexamethasone suppression of VEGF transcription
    EC50_DEX = 6.0,      # ug/L
    DEX_LP   = 0.60,     # dexamethasone direct anti-permeability weight
    HB50     = 0.899831,      # g/dL cavity haemoglobin at half-maximal drive
    HB_HILL  = 2.0,      # steepness of the haem -> inflammation switch
    D50_PAT  = 5.0,      # mm, cavity thickness at half patency
    N_PAT    = 2.0,      # steepness of the apposition term
    K_FUSE   = 0.079721,    # /d membrane involution when layers are apposed
    K_MAC    = 1.820,    # /d macrophage / eosinophil recruitment
    KD_MAC   = 0.350,    # /d
    MAC_MAX  = 2.20,
    K_ANG2   = 3.4, KEL_ANG2 = 0.35,
    K_MMP9   = 62.0, KEL_MMP9 = 0.30, ATV_MMP = 0.45,
    K_IL6    = 210.0, KEL_IL6 = 1.2,
    K_IL8    = 175.0, KEL_IL8 = 1.0,
    IL6_REF  = 900.0, IL8_REF = 800.0,
    K_MEMO   = 0.045, KD_MEMO = 0.012,   # outer membrane growth
    K_SEPT   = 0.020, KD_SEPT = 0.004,   # septation

    # ---- brain mechanics ------------------------------------------
    ICP0     = 8.0,      # mmHg
    K_PV     = 0.022,    # /mL pressure-volume exponent
    COMP_F   = 0.60,     # fraction of excess volume absorbed by deformation
    K_COMP   = 0.20,     # /d loading rate (fast)
    K_REEXP0 = 0.094087,    # /d unloading rate (slow) in a non-atrophic brain
    LAM_ATR  = 26.0,     # mL, atrophy constant slowing re-expansion
    MLS_MAX  = 25.0,     # mm
    K_MLS    = 45.0,     # mL
    K_CBF    = 0.055, KREC_CBF = 0.10,
    K_NEUR   = 0.0075, KREP_NEUR= 0.010,

    # ---- clinical -------------------------------------------------
    K_SYMP   = 0.42, KREC_SYMP= 0.30,
    K_COGL   = 0.055, KREC_COG = 0.030,
    LAM_REOP = 0.001731,    # /d peak reoperation hazard
    D_REOP   = 10.0,     # mm thickness at which reoperation is considered
    S_REOP   = 2.2,      # mm steepness
    LAM_THR  = 0.00045,  # /d thromboembolic hazard off anticoagulation
    LAM_DTH  = 0.00016,  # /d baseline mortality

    # ---- steroid toxicity -----------------------------------------
    GLU0     = 5.4, K_GLU = 3.6, KEL_GLU = 1.1,
    K_INF    = 0.055, KEL_INF = 0.035,
    K_MYO    = 0.085, KEL_MYO = 0.045,

    # ---- PK: dexamethasone (2-cmt + effect) -----------------------
    DEX_KA = 6.0, DEX_F = 0.80, DEX_V1 = 35.0, DEX_V2 = 30.0,
    DEX_Q = 360.0, DEX_CL = 250.0, DEX_KE0 = 1.2,
    # ---- PK: atorvastatin -----------------------------------------
    ATV_KA = 6.0, ATV_F = 0.14, ATV_V1 = 400.0, ATV_V2 = 5200.0,
    ATV_Q = 900.0, ATV_CL = 8000.0, ATV_KE0 = 0.08,
    # ---- PK: tranexamic acid --------------------------------------
    TXA_V1 = 12.0, TXA_V2 = 18.0, TXA_Q = 288.0, TXA_CL = 180.0,
    TXA_KDEGH = 0.05,
    # ---- PK: apixaban ---------------------------------------------
    DOAC_V1 = 21.0, DOAC_V2 = 12.0, DOAC_Q = 60.0, DOAC_CL = 79.0,
    DOAC_KDEGH = 0.05,
)

NST = 50
(IA_DEXD, IA_DEXC, IA_DEXP, ICE_DEX,
 IA_ATVD, IA_ATVC, IA_ATVP, ICE_ATV,
 IA_TXAC, IA_TXAP, IA_TXAH,
 IA_DOAC, IA_DOAP, IA_DOAH,
 IV_HEM, IM_HB, IM_PROT, IM_FIB, IM_FDP, IM_TPA, IC_PLS, IC_PAI,
 IN_CAP, IN_MAT, IC_VEGF, IC_ANG2, IC_MMP9, IC_IL6, IC_IL8, IN_MAC,
 IM_MEMO, IM_SEPT, IP_MMA, IF_DRN, IN_PC, IF_PLT, IC_FBG,
 IV_COMP, IX_MLSI, IX_CBF, IN_NEUR,
 IS_SYMP, IS_COG, IGLU, IX_INF, IX_MYO,
 IH_REOP, IH_THR, IH_DTH, IA_EXT) = range(NST)


# ---------------------------------------------------------------------
# 2.  ALGEBRAIC LAYER  (everything the ODEs need that is not a state)
# ---------------------------------------------------------------------
def algebra(y, p, infusions=None):
    """Return a dict of derived quantities. Kept separate so that the
    verification harness can interrogate the intermediate fluxes."""
    a = {}
    V = max(y[IV_HEM], 1e-4)

    # -- geometry ---------------------------------------------------
    # Area comes from the MEMBRANE's extent; thickness is then whatever
    # the fluid volume implies over that area.
    A_mem = max(p['A_MAX'] * min(1.0, max(y[IA_EXT], 1e-3)), 4.0)  # cm2
    A100 = A_mem / 100.0
    dmean_cm = V / A_mem                                          # cm
    dmax = min(p['D_CAP'], 10.0 * dmean_cm / p['D_MEANF'])        # mm
    a['A_mem'], a['A100'], a['dmax'] = A_mem, A100, dmax
    # Cavity patency.  A subdural space is not a rigid container: once the
    # brain re-apposes the dura the two membrane layers touch, filtration
    # has nowhere to go, the neovessels are compressed and the leaves
    # fuse.  Without this a 0.5 mL post-drainage film was still refilled
    # at 16 mL/day, so no operation could ever succeed.
    a['f_pat'] = dmax**p['N_PAT'] / (dmax**p['N_PAT'] + p['D50_PAT']**p['N_PAT'])

    # -- drug effect sites -----------------------------------------
    ce_dex = y[ICE_DEX]
    ce_atv = y[ICE_ATV]
    a['E_DEXV'] = p['EMAX_DEXV'] * ce_dex / (p['EC50_DEX'] + ce_dex)
    a['E_ATV']  = p['EMAX_ATV']  * ce_atv / (p['EC50_ATV'] + ce_atv)
    f_dex_vegf  = 1.0 / (1.0 + a['E_DEXV'])
    f_dex_lp    = 1.0 / (1.0 + p['DEX_LP'] * a['E_DEXV'])
    a['f_dex_vegf'] = f_dex_vegf
    a['f_dex_lp']   = f_dex_lp

    # -- vessel population ------------------------------------------
    ncap, nmat = max(y[IN_CAP], 0.0), max(y[IN_MAT], 0.0)
    ntot = ncap + nmat + 1e-9
    fmat = nmat / ntot
    # pericyte coverage tracks maturity but is accelerated by statin
    pc = min(1.0, max(0.0, y[IN_PC]))
    sig_eff = p['SIG_IMM'] + (p['SIG_MAT'] - p['SIG_IMM']) * (0.5 * fmat + 0.5 * pc)
    Lp_eff  = (p['LP_IMM'] * ncap + p['LP_MAT'] * nmat) * f_dex_lp
    a['sig_eff'], a['Lp_eff'], a['fmat'] = sig_eff, Lp_eff, fmat

    # -- cavity composition ------------------------------------------
    C_PROT = (y[IM_PROT] / V) * 100.0            # g/dL   (M_PROT in g)
    C_HB   = (y[IM_HB]   / V) / 10.0             # g/dL   (M_HB in mg)
    C_FDP  = (y[IM_FDP]  / V) * 1000.0           # ug/mL  (M_FDP in mg)
    C_TPA  = (y[IM_TPA]  / V)                    # ng/mL  (M_TPA in ng)
    C_TXAH = (y[IA_TXAH] / V) * 1000.0           # mg/L   (A_TXAH in mg)
    C_DOAH = (y[IA_DOAH] / V) * 1e6              # ng/mL  (A_DOAH in mg)
    a.update(C_PROT=C_PROT, C_HB=C_HB, C_FDP=C_FDP, C_TPA=C_TPA,
             C_TXAH=C_TXAH, C_DOAH=C_DOAH)
    a['HU'] = 10.0 + 5.5 * C_HB + 1.5 * C_PROT
    hb_n = C_HB / p['HB50']
    a['haem_drive'] = hb_n**p['HB_HILL'] / (1.0 + hb_n**p['HB_HILL'])

    # -- brain mechanics ---------------------------------------------
    R = p['V_RES'] + max(y[IV_COMP], 0.0)        # effective reserve
    excess = sp(V - R, 1.5)
    a['R'] = R
    a['excess'] = excess
    # ICP is TWO-SIDED.  When the brain has not re-expanded into a
    # drained cavity the pressure is SUB-baseline; that is why the space
    # persists, why absorption stalls, and why a suction drain does not
    # simply empty the head.  Flooring ICP at baseline let a 48 h drain
    # take a 78 mL haematoma to 0.06 mL.
    a['ICP'] = p['P_MIN'] + p['ICP0'] * np.exp(
        np.clip(p['K_PV'] * (V - R), -30.0, 30.0))
    a['MLS'] = p['MLS_MAX'] * excess / (excess + p['K_MLS'])
    a['space'] = sp(R - V, 1.5)                 # room available to refill

    # -- Starling fluxes ----------------------------------------------
    pi_hem = p['PI_COEF'] * C_PROT
    dP = (p['P_CAP'] - a['ICP']) - sig_eff * (p['PI_PL'] - pi_hem)
    a['pi_hem'] = pi_hem
    a['dP_net'] = dP
    J_ex = Lp_eff * A100 * sp(dP, 0.8) * a['f_pat']
    a['J_ex'] = J_ex

    # -- haemostatic competence ---------------------------------------
    f_fdp  = 1.0 / (1.0 + C_FDP / p['K_FDP'])
    f_doac = 1.0 / (1.0 + C_DOAH / p['IC50_DOAC'])
    f_plt  = min(1.0, max(0.0, y[IF_PLT]))
    f_fbg  = y[IC_FBG] / (y[IC_FBG] + 90.0)
    H = f_fdp * f_doac * f_plt * (0.45 + 0.55 * f_fbg)
    a['H'] = H
    a['f_fdp'], a['f_doac'] = f_fdp, f_doac

    # -- rebleed -------------------------------------------------------
    frag = 1.0 + p['A_MMP'] * y[IC_MMP9] / (p['KM_MMP'] + y[IC_MMP9]) \
               + p['A_ANG'] * y[IC_ANG2] / (p['KM_ANG'] + y[IC_ANG2])
    supply = p['MMA_FLOOR'] + (1.0 - p['MMA_FLOOR']) * min(1.0, max(0.0, y[IP_MMA]))
    J_rb = p['K_RB'] * ncap * A100 * frag * (1.0 - H) * supply * a['f_pat']
    a['frag'], a['supply'], a['J_rb'] = frag, supply, J_rb

    # -- absorption -----------------------------------------------------
    # Efflux to the dural lymphatics / arachnoid granulations is
    # PRESSURE-DRIVEN, not merely kinetic.  Without this term a 3 mL
    # collection clears in a day (so the model resolved every incident
    # haematoma before its membrane could form) and a 140 mL one was
    # never pushed back hard enough.
    p_drive = sp(a['ICP'] - p['P_SINK'], 0.8) / p['P_REF']
    # The area-borne term is multiplied by patency: a collapsed cavity has
    # no fluid to clear.  Omitting this let absorption keep running on a
    # cavity that no longer existed and drove the volume to -6 mL.
    J_abs = (p['K_LYMPH'] * A100 * a['f_pat'] + p['K_ARACH'] * V) * p_drive / \
            (1.0 + p['SEPT_ABS'] * y[IM_SEPT])
    a['p_drive'] = p_drive
    J_drn = (p['K_DRAIN'] * min(1.0, max(0.0, y[IF_DRN])) * V
             * sp(a['ICP'] - p['P_DRAIN'], 0.8) / p['P_REF'])
    a['J_abs'], a['J_drn'] = J_abs, J_drn

    # -- plasma drug concentrations --------------------------------------
    a['C_DEX']  = y[IA_DEXC] / p['DEX_V1'] * 1000.0     # ug/L
    a['C_ATV']  = y[IA_ATVC] / p['ATV_V1'] * 1000.0     # ug/L
    a['C_TXA']  = y[IA_TXAC] / p['TXA_V1']              # mg/L
    a['C_DOAC'] = y[IA_DOAC] / p['DOAC_V1'] * 1000.0    # ng/mL
    return a


# ---------------------------------------------------------------------
# 3.  RIGHT-HAND SIDE
# ---------------------------------------------------------------------
def rhs(t, y, p, rate):
    """rate = dict of zero-order input rates (mg/d) set by the scenario."""
    y = np.maximum(y, 0.0)
    a = algebra(y, p)
    d = np.zeros(NST)
    V = max(y[IV_HEM], 1e-4)

    # ---------------- PK -------------------------------------------
    d[IA_DEXD] = rate.get('dex', 0.0) * p['DEX_F'] - p['DEX_KA'] * y[IA_DEXD]
    d[IA_DEXC] = (p['DEX_KA'] * y[IA_DEXD]
                  - p['DEX_CL'] * y[IA_DEXC] / p['DEX_V1']
                  - p['DEX_Q'] * (y[IA_DEXC] / p['DEX_V1'] - y[IA_DEXP] / p['DEX_V2']))
    d[IA_DEXP] = p['DEX_Q'] * (y[IA_DEXC] / p['DEX_V1'] - y[IA_DEXP] / p['DEX_V2'])
    d[ICE_DEX] = p['DEX_KE0'] * (a['C_DEX'] - y[ICE_DEX])

    d[IA_ATVD] = rate.get('atv', 0.0) * p['ATV_F'] - p['ATV_KA'] * y[IA_ATVD]
    d[IA_ATVC] = (p['ATV_KA'] * y[IA_ATVD]
                  - p['ATV_CL'] * y[IA_ATVC] / p['ATV_V1']
                  - p['ATV_Q'] * (y[IA_ATVC] / p['ATV_V1'] - y[IA_ATVP] / p['ATV_V2']))
    d[IA_ATVP] = p['ATV_Q'] * (y[IA_ATVC] / p['ATV_V1'] - y[IA_ATVP] / p['ATV_V2'])
    d[ICE_ATV] = p['ATV_KE0'] * (a['C_ATV'] - y[ICE_ATV])

    d[IA_TXAC] = (rate.get('txa', 0.0)
                  - p['TXA_CL'] * y[IA_TXAC] / p['TXA_V1']
                  - p['TXA_Q'] * (y[IA_TXAC] / p['TXA_V1'] - y[IA_TXAP] / p['TXA_V2'])
                  - (a['J_ex'] + a['J_rb']) * a['C_TXA'] / 1000.0)
    d[IA_TXAP] = p['TXA_Q'] * (y[IA_TXAC] / p['TXA_V1'] - y[IA_TXAP] / p['TXA_V2'])
    # THE POINT: tranexamic acid reaches its site of action by being
    # carried in with the exudate, not by diffusion.  Its cavity
    # half-life is therefore V/J_abs (~6 days), not its 3 h plasma one.
    d[IA_TXAH] = ((a['J_ex'] + a['J_rb']) * a['C_TXA'] / 1000.0
                  - (a['J_abs'] + a['J_drn']) * (y[IA_TXAH] / V)
                  - p['TXA_KDEGH'] * y[IA_TXAH])

    d[IA_DOAC] = (rate.get('doac', 0.0)
                  - p['DOAC_CL'] * y[IA_DOAC] / p['DOAC_V1']
                  - p['DOAC_Q'] * (y[IA_DOAC] / p['DOAC_V1'] - y[IA_DOAP] / p['DOAC_V2'])
                  - (a['J_ex'] + a['J_rb']) * a['C_DOAC'] / 1e6)
    d[IA_DOAP] = p['DOAC_Q'] * (y[IA_DOAC] / p['DOAC_V1'] - y[IA_DOAP] / p['DOAC_V2'])
    d[IA_DOAH] = ((a['J_ex'] + a['J_rb']) * a['C_DOAC'] / 1e6
                  - (a['J_abs'] + a['J_drn']) * (y[IA_DOAH] / V)
                  - p['DOAC_KDEGH'] * y[IA_DOAH])

    # ---------------- the volume balance ---------------------------
    d[IV_HEM] = a['J_ex'] + a['J_rb'] - a['J_abs'] - a['J_drn']

    # ---------------- cavity contents ------------------------------
    # V is FLOORED here.  Washout of the cavity solutes is a fractional rate
    # with the volume in the denominator, so as a successfully treated cavity
    # empties this term diverges and the solute equations become violently
    # stiff -- the embolisation arms simply hung the integrator.  Below about
    # half a millilitre there is no cavity left to have a concentration in.
    out_frac = (a['J_abs'] + a['J_drn']) / max(V, p['V_FLOOR'])
    d[IM_HB]   = (p['HB_BLOOD'] * a['J_rb']
                  - p['K_HBDEG'] * y[IM_HB] - out_frac * y[IM_HB])
    # Convective protein entry is (1-sigma)*C_plasma; a diffusive term is
    # added because the neomembrane is not a dialysis bag.  Together with
    # J_ex + J_rb = J_abs at steady volume these BOUND the cavity protein
    # below plasma protein, which is what stops the oncotic term from
    # becoming a driver of its own filtration.
    d[IM_PROT] = (a['J_ex'] * (1.0 - a['sig_eff']) * p['PROT_PL']
                  + a['J_rb'] * p['PROT_PL']
                  + p['PS_PROT'] * a['A100'] * (p['PROT_PL'] - y[IM_PROT] / V)
                  - out_frac * y[IM_PROT])

    # fibrin is laid down by each rebleed in proportion to how well the
    # patient can actually clot, then eaten by plasmin
    fib_lys = p['K_LYS'] * y[IC_PLS] * y[IM_FIB]
    d[IM_FIB] = (p['K_FIBF'] * a['J_rb'] * a['H'] * (y[IC_FBG] / 300.0)
                 - fib_lys - out_frac * y[IM_FIB])
    d[IM_FDP] = (p['K_FDPY'] * fib_lys - p['KDEG_FDP'] * y[IM_FDP]
                 - out_frac * y[IM_FDP])
    d[IM_TPA] = (p['K_TPA'] * y[IN_CAP] * a['A100']
                 - (p['KDEG_TPA'] + p['K_PAI'] * y[IC_PAI]) * y[IM_TPA]
                 - out_frac * y[IM_TPA])
    f_txa = 1.0 / (1.0 + a['C_TXAH'] / p['IC50_TXA'])
    d[IC_PLS] = p['K_PLS'] * a['C_TPA'] * f_txa - p['KEL_PLS'] * y[IC_PLS]
    d[IC_PAI] = p['K_PAISYN'] * (1.0 + 0.4 * y[IN_MAC]) - p['K_PAIDEG'] * y[IC_PAI]

    # ---------------- membrane biology ------------------------------
    vegf_n = y[IC_VEGF] / p['VEGF_REF']
    k_mat_eff = (p['K_MAT'] * (1.0 + a['E_ATV'])
                 * (1.0 - p['MAT_SUPP'] * a['haem_drive']))
    d[IN_CAP] = (p['K_SPROUT'] * vegf_n * a['supply']
                 * max(0.0, 1.0 - (y[IN_CAP] + y[IN_MAT]) / p['N_MAX'])
                 - k_mat_eff * y[IN_CAP] - p['KD_CAP'] * y[IN_CAP]
                 - p['K_FUSE'] * (1.0 - a['f_pat']) * y[IN_CAP])
    mmp_des = p['K_MMPDES'] * y[IC_MMP9] / (p['KM_MMPD'] + y[IC_MMP9])
    d[IN_MAT] = k_mat_eff * y[IN_CAP] - (p['KD_MAT'] + mmp_des) * y[IN_MAT]
    d[IN_PC]  = 0.010 * (1.0 + a['E_ATV']) * (a['fmat'] - y[IN_PC])

    d[IC_VEGF] = (p['K_VEGF'] * y[IN_MAC] * a['f_dex_vegf']
                  - p['KEL_VEGF'] * y[IC_VEGF])
    # THE SWITCH.  Haem is the stimulus that recruits the macrophages that
    # make the VEGF that builds the leaky vessels that bleed and release
    # more haem.  That loop is closed, so its steepness decides whether
    # the disease is bistable.  With a LINEAR haem term the model has one
    # fixed point at 113 mL and surgery cannot cure anybody -- every
    # scenario refilled to the untreated volume and recurrence was 100%
    # by construction.  A Hill exponent of 2 gives the loop a threshold:
    # wash the haem out and it extinguishes, leave some behind and it
    # relights.  This is the arithmetic content of "irrigate until clear".
    # The IL-8 amplification is GATED ON HAEM.  Written as a free autocrine
    # term it is a second positive feedback with no external input, so the
    # inflammatory state never switches off: macrophages settled at 0.63 with
    # the haem drive at 0.12, the membrane always relit, surgery never cured
    # anybody and MMA embolisation was the only thing that worked (P(reop)
    # 0.0007 vs a reported 0.041).
    il8_n = y[IC_IL8] / p['IL8_REF']
    d[IN_MAC] = (p['K_MAC'] * a['haem_drive'] * (1.0 + p['W_IL8'] * il8_n)
                 * a['f_dex_vegf']
                 * max(0.0, 1.0 - y[IN_MAC] / p['MAC_MAX'])
                 - p['KD_MAC'] * y[IN_MAC])
    d[IC_ANG2] = p['K_ANG2'] * y[IN_CAP] * a['f_dex_vegf'] - p['KEL_ANG2'] * y[IC_ANG2]
    d[IC_MMP9] = (p['K_MMP9'] * y[IN_MAC] * a['f_dex_vegf']
                  / (1.0 + p['ATV_MMP'] * a['E_ATV']) - p['KEL_MMP9'] * y[IC_MMP9])
    d[IC_IL6]  = p['K_IL6'] * y[IN_MAC] * a['f_dex_vegf'] - p['KEL_IL6'] * y[IC_IL6]
    d[IC_IL8]  = (p['K_IL8'] * y[IN_MAC] * a['haem_drive'] * a['f_dex_vegf']
                  - p['KEL_IL8'] * y[IC_IL8])

    d[IM_MEMO] = p['K_MEMO'] * y[IN_MAC] * (1.0 - y[IM_MEMO] / 2.0) - p['KD_MEMO'] * y[IM_MEMO]
    d[IM_SEPT] = p['K_SEPT'] * y[IM_MEMO] * (1.0 - y[IM_SEPT] / 2.0) - p['KD_SEPT'] * y[IM_SEPT]

    # Spreading and fusing are blended with a sigmoid rather than switched
    # with an if/else.  A hard branch on a state variable makes the solver
    # chatter on the switching surface: with the discontinuous form, a late
    # embolisation (scenario 07) stalled LSODA indefinitely.
    ext_target = V / (V + p['K_A'])
    w_ext = 1.0 / (1.0 + np.exp(-np.clip((ext_target - y[IA_EXT]) / 0.01,
                                         -60.0, 60.0)))
    d[IA_EXT] = (w_ext * p['K_SPREAD'] * (ext_target - y[IA_EXT])
                 - (1.0 - w_ext) * p['K_FUSEX'] * (1.0 - a['f_pat']) * y[IA_EXT])

    d[IP_MMA]  = p['K_RECAN'] * (1.0 - y[IP_MMA])
    d[IF_DRN]  = 0.0
    d[IF_PLT]  = 0.05 * (1.0 - y[IF_PLT])
    d[IC_FBG]  = 0.25 * (300.0 - y[IC_FBG])

    # ---------------- brain mechanics -------------------------------
    # Loads fast, unloads slowly, and unloads more slowly the more
    # atrophic the brain.  This asymmetry is the whole recurrence story.
    target = p['COMP_F'] * max(0.0, V - p['V_RES'])
    k_reexp = p['K_REEXP0'] * np.exp(-p['V_RES'] / p['LAM_ATR'])
    w_cmp = 1.0 / (1.0 + np.exp(-np.clip((target - y[IV_COMP]) / 0.05,
                                         -60.0, 60.0)))
    d[IV_COMP] = (w_cmp * p['K_COMP'] * (target - y[IV_COMP])
                  - (1.0 - w_cmp) * k_reexp * (y[IV_COMP] - target))

    d[IX_MLSI] = a['MLS'] / 10.0
    d[IX_CBF]  = p['K_CBF'] * (a['MLS'] / 10.0) - p['KREC_CBF'] * y[IX_CBF]
    d[IN_NEUR] = -p['K_NEUR'] * y[IX_CBF] + p['KREP_NEUR'] * (1.0 - y[IN_NEUR])

    # ---------------- clinical ---------------------------------------
    drive = 0.55 * (a['MLS'] / 10.0) + 0.45 * (1.0 - y[IN_NEUR]) / 0.3
    d[IS_SYMP] = p['K_SYMP'] * drive - p['KREC_SYMP'] * y[IS_SYMP]
    d[IS_COG]  = p['K_COGL'] * drive - p['KREC_COG'] * y[IS_COG]

    d[IGLU]  = p['K_GLU'] * a['E_DEXV'] / p['EMAX_DEXV'] - p['KEL_GLU'] * (y[IGLU] - p['GLU0'])
    d[IX_INF] = p['K_INF'] * a['E_DEXV'] - p['KEL_INF'] * y[IX_INF]
    d[IX_MYO] = p['K_MYO'] * a['E_DEXV'] - p['KEL_MYO'] * y[IX_MYO]

    # ---------------- hazards -----------------------------------------
    reop_geo = 1.0 / (1.0 + np.exp(-(a['dmax'] - p['D_REOP']) / p['S_REOP']))
    reop_sym = y[IS_SYMP] / (y[IS_SYMP] + 0.55)
    d[IH_REOP] = p['LAM_REOP'] * reop_geo * reop_sym
    # thromboembolic hazard is carried by the patient who is OFF their
    # anticoagulant; it is the price of the haemostatic safety above
    d[IH_THR]  = p['LAM_THR'] / (1.0 + a['C_DOAC'] / 30.0)
    d[IH_DTH]  = p['LAM_DTH'] * (1.0 + 1.6 * y[IX_INF] + 0.9 * (a['MLS'] / 10.0))
    return d


# ---------------------------------------------------------------------
# 4.  INITIAL CONDITIONS
# ---------------------------------------------------------------------
def y0_chronic(p, V0=None):
    """A mature, symptomatic cSDH sitting near its own fixed point."""
    y = np.zeros(NST)
    y[IV_HEM]  = 78.0 if V0 is None else V0
    y[IM_HB]   = 2600.0
    y[IM_PROT] = 4.2
    y[IM_FIB]  = 42.0
    y[IM_FDP]  = 15.0
    y[IM_TPA]  = 2100.0
    y[IC_PLS]  = 0.9
    y[IC_PAI]  = P['PAI0']
    y[IN_CAP]  = 1.10
    y[IN_MAT]  = 0.35
    y[IC_VEGF] = 150.0
    y[IC_ANG2] = 11.0
    y[IC_MMP9] = 230.0
    y[IC_IL6]  = 950.0
    y[IC_IL8]  = 820.0
    y[IN_MAC]  = 1.05
    y[IM_MEMO] = 1.10
    y[IM_SEPT] = 0.45
    y[IP_MMA]  = 1.0
    y[IF_DRN]  = 0.0
    y[IN_PC]   = 0.24
    y[IF_PLT]  = 1.0
    y[IC_FBG]  = 300.0
    y[IV_COMP] = P['COMP_F'] * max(0.0, y[IV_HEM] - p['V_RES'])
    y[IN_NEUR] = 0.86
    y[IS_SYMP] = 0.95
    y[IS_COG]  = 0.9
    y[IGLU]    = p['GLU0']
    y[IA_EXT]  = y[IV_HEM] / (y[IV_HEM] + p['K_A'])
    return y


def y0_fresh(p, V0=12.0):
    """A thin subacute subdural collection a week or two after a minor
    head injury: blood, no neomembrane, no neovessels yet."""
    y = np.zeros(NST)
    y[IV_HEM]  = V0
    y[IM_HB]   = p['HB_BLOOD'] * V0        # frank blood
    y[IM_PROT] = p['PROT_PL'] * V0
    y[IM_FIB]  = 60.0 * V0 / 12.0
    y[IM_FDP]  = 0.5
    y[IM_TPA]  = 20.0
    y[IC_PLS]  = 0.05
    y[IC_PAI]  = p['K_PAISYN'] / p['K_PAIDEG']
    y[IN_CAP]  = 0.04
    y[IN_MAT]  = 0.01
    y[IC_VEGF] = 20.0
    y[IC_ANG2] = 1.0
    y[IC_MMP9] = 25.0
    y[IC_IL6]  = 60.0
    y[IC_IL8]  = 50.0
    y[IN_MAC]  = 0.06
    y[IM_MEMO] = 0.02
    y[IM_SEPT] = 0.00
    y[IP_MMA]  = 1.0
    y[IF_DRN]  = 0.0
    y[IN_PC]   = 0.20
    y[IF_PLT]  = 1.0
    y[IC_FBG]  = 300.0
    y[IV_COMP] = p['COMP_F'] * max(0.0, V0 - p['V_RES'])
    y[IN_NEUR] = 1.0
    y[IGLU]    = p['GLU0']
    y[IA_EXT]  = V0 / (V0 + p['K_A'])
    return y


def run_in(p, V_pres=78.0, tmax=400.0, V0=12.0, doac_rate=0.0):
    """Grow the haematoma from a fresh thin collection and stop when it
    reaches the presentation volume.  This replaces an assumed initial
    condition with a MODEL OUTPUT, and it makes the latency from injury
    to presentation a falsifiable prediction (clinically 3-8 weeks).

    Returns (state_at_presentation, latency_days, reached?)."""
    y = y0_fresh(p, V0)
    rate = {'doac': doac_rate} if doac_rate else {}

    def hit(t, yy, p, rate):
        return yy[IV_HEM] - V_pres
    hit.terminal, hit.direction = True, 1

    sol = solve_ivp(rhs, (0.0, tmax), y, args=(p, rate), method='LSODA',
                    rtol=1e-6, atol=1e-8, max_step=2.0, events=hit)
    yend = np.maximum(sol.y[:, -1], 0.0)
    # The hazard clocks and the shift-exposure integral must start at
    # PRESENTATION.  Carrying the run-in's accumulation forward gave every
    # arm, treated or not, a 58% reoperation probability at t = 0.
    for i in (IH_REOP, IH_THR, IH_DTH, IX_MLSI):
        yend[i] = 0.0
    reached = len(sol.t_events[0]) > 0
    return yend, float(sol.t[-1]), reached


def run_in_trace(p, tmax=180.0, V0=12.0):
    """Same run-in but returning the whole trajectory, for the
    natural-history plots and for the latency validation."""
    y = y0_fresh(p, V0)
    te = np.arange(0.0, tmax + 1e-9, 0.25)
    sol = solve_ivp(rhs, (0.0, tmax), y, args=(p, {}), method='LSODA',
                    rtol=1e-6, atol=1e-8, max_step=2.0, t_eval=te)
    return sol.t, sol.y.T


# ---------------------------------------------------------------------
# 5.  SCENARIO ENGINE
# ---------------------------------------------------------------------
class Scn:
    """A scenario is a list of (time, action) plus infusion windows."""
    def __init__(self, name, tend=180.0, p=None, V_pres=78.0, doac_runin=0.0):
        self.name, self.tend = name, tend
        self.p = dict(P)
        if p: self.p.update(p)
        self.V_pres = V_pres
        self.doac_runin = doac_runin
        self.latency = None
        self.events = []      # (t, fn(y, p))
        self.windows = []     # (t0, t1, drug, rate mg/d)

    def add(self, t, fn):  self.events.append((t, fn)); return self
    def infuse(self, t0, t1, drug, rate): self.windows.append((t0, t1, drug, rate)); return self

    # ---- standard actions --------------------------------------------
    def surgery(self, t, residual=0.22, drain_days=0.0, wash=0.15):
        """Burr-hole craniostomy.  Empties the STOCK; the membrane, the
        vessel population and the un-re-expanded brain are all left
        exactly as they were.  That is why it recurs.

        `wash` is the irrigation: it scales the cavity CONCENTRATIONS, not
        just the amounts.  Getting this wrong (scaling amounts and volume
        by the same factor, so that every concentration was preserved)
        left the haem stimulus untouched and made irrigation a no-op."""
        def do(y, p):
            y[IV_HEM] *= residual
            for i in (IM_HB, IM_PROT, IM_FIB, IM_FDP, IM_TPA, IA_TXAH, IA_DOAH):
                y[i] *= residual * wash
            y[IM_SEPT] *= 0.6            # irrigation breaks some loculi
            y[IF_DRN] = 1.0 if drain_days > 0 else 0.0
        self.add(t, do)
        if drain_days > 0:
            self.add(t + drain_days, lambda y, p: y.__setitem__(IF_DRN, 0.0))
        return self

    def embolise(self, t, residual_perf=0.05):
        """MMA embolisation.  Leaves the STOCK untouched and starves the
        membrane.  Effect is therefore an integral over weeks."""
        self.add(t, lambda y, p: y.__setitem__(IP_MMA, residual_perf))
        return self

    def dex(self, t0=0.0, days=14.0, dose=16.0, taper=True):
        if taper:
            self.infuse(t0, t0 + 8, 'dex', dose)
            self.infuse(t0 + 8, t0 + 11, 'dex', dose * 0.5)
            self.infuse(t0 + 11, t0 + days, 'dex', dose * 0.25)
        else:
            self.infuse(t0, t0 + days, 'dex', dose)
        return self

    def atorva(self, t0=0.0, days=56.0, dose=20.0):
        return self.infuse(t0, t0 + days, 'atv', dose)

    def txa(self, t0=0.0, days=90.0, dose=750.0, F=0.34):
        return self.infuse(t0, t0 + days, 'txa', dose * F)

    def doac(self, t0, t1, dose=10.0, F=0.5):
        return self.infuse(t0, t1, 'doac', dose * F)

    # ---- run ------------------------------------------------------------
    def run(self, dt=0.25):
        p = self.p
        # t = 0 is PRESENTATION, reached by growing the haematoma from a
        # fresh thin collection -- so the membrane, the vessel population,
        # the fibrinolytic milieu and the brain's compression deficit all
        # arrive with histories that are consistent with each other.
        y, self.latency, ok = run_in(p, self.V_pres, doac_rate=self.doac_runin)
        if not ok:
            raise RuntimeError(f"{self.name}: run-in never reached V_pres")
        tmarks = sorted(set([0.0, self.tend]
                            + [t for t, _ in self.events]
                            + [w[0] for w in self.windows]
                            + [w[1] for w in self.windows]))
        tmarks = [t for t in tmarks if 0.0 <= t <= self.tend]
        T, Y = [0.0], [y.copy()]
        for k in range(len(tmarks) - 1):
            t0, t1 = tmarks[k], tmarks[k + 1]
            for te, fn in self.events:
                if abs(te - t0) < 1e-9:
                    fn(y, p)
            rate = {}
            for (w0, w1, drug, r) in self.windows:
                if w0 <= t0 + 1e-9 and t1 <= w1 + 1e-9:
                    rate[drug] = rate.get(drug, 0.0) + r
            if t1 - t0 < 1e-9:
                continue
            teval = np.arange(t0, t1 + 1e-9, dt)
            if teval[-1] < t1 - 1e-9:
                teval = np.append(teval, t1)
            sol = solve_ivp(rhs, (t0, t1), y, args=(p, rate), t_eval=teval,
                            method='LSODA', rtol=1e-6, atol=1e-8, max_step=2.0)
            if not sol.success:
                raise RuntimeError(f"{self.name}: integration failed at t={t0}")
            for j in range(1, sol.y.shape[1]):
                T.append(sol.t[j]); Y.append(sol.y[:, j].copy())
            y = np.maximum(sol.y[:, -1].copy(), 0.0)
        for te, fn in self.events:
            if abs(te - self.tend) < 1e-9:
                fn(y, p)
        return np.array(T), np.array(Y), p


# ---------------------------------------------------------------------
# 6.  OUTPUT HELPERS
# ---------------------------------------------------------------------
def series(T, Y, p, key):
    return np.array([algebra(Y[i], p)[key] for i in range(len(T))])

def at(T, X, t):
    return float(np.interp(t, T, X))

def reop_prob(T, Y, t):
    return 1.0 - np.exp(-at(T, Y[:, IH_REOP], t))
