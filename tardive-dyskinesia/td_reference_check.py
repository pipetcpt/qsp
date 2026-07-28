#!/usr/bin/env python3
# =============================================================================
# Tardive Dyskinesia (TD) QSP model — independent numerical reference
# 지연성 운동이상증 QSP 모델 — mrgsolve 모델의 독립 검증 구현
#
# This file is a numpy/scipy transcription of the SAME equations that appear in
# td_mrgsolve_model.R ($ODE block).  It exists so that every number quoted in
# README.md can be re-derived without R / mrgsolve installed.
#
#   python3 td_reference_check.py            # run all analyses, print tables
#   python3 td_reference_check.py --quick    # scenarios only
#
# 40 states (19 PK/exposure + 21 disease / basal-ganglia / endpoint).
# Time unit = DAYS.  Concentrations = ng/mL.  Dosing is zero-order daily input
# (DOSE mg/day delivered continuously) — identical in the R model, so the two
# implementations are numerically comparable.
# =============================================================================
import sys
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------------
# state index map
# ----------------------------------------------------------------------------
S = dict(
    AP_gut=0, AP_cen=1, AP_per=2, AP_lai=3,
    VAL_gut=4, VAL_cen=5, NBI_cen=6, NBI_per=7,
    DTB_gut=8, HTB_cen=9, HTB_per=10,
    CLZ_gut=11, CLZ_cen=12,
    AMA_gut=13, AMA_cen=14,
    ACH_gut=15, ACH_cen=16,
    GKB=17, BONT=18,
    CE_AP=19,
    DA_CYT=20, DA_VES=21, DA_SYN=22, DEPOL=23,
    RUP=24, ROS=25, SDAM=26,
    IND=27, GPE=28, STN=29, GPI=30, THAL=31,
    AIMS=32, PSYCH=33, PARK=34, DEPR=35, QTC=36,
    ADHER=37, DOSEMULT=38, CUMOCC=39,
)
NS = len(S)

# ----------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------
P0 = dict(
    # --- antipsychotic (risperidone active-moiety equivalents) -------------
    KA_AP=20.0, CL_AP=120.0, V2_AP=100.0, V3_AP=80.0, Q_AP=30.0,
    KA_LAI=0.16, F_LAI=0.0,          # fraction of AP dose given as LAI depot
    KE0_AP=12.0,                     # effect-site equilibration (striatum)
    EC50_D2=12.0, HILL_D2=1.25,
    # --- clozapine ---------------------------------------------------------
    KA_CLZ=15.0, CL_CLZ=900.0, V_CLZ=500.0,
    EC50_D2_CLZ=900.0, EC50_CLZ_EFF=250.0,
    # --- valbenazine -> NBI-98782 -----------------------------------------
    KA_VAL=10.0, CL_VAL=500.0, V_VAL=500.0, FM_VAL=0.30,
    CL_NBI=1000.0, V_NBI=900.0, Q_NBI=200.0, V3_NBI=600.0,
    EC50_VMAT_NBI=26.0,
    # --- deutetrabenazine -> (alpha+beta)-HTBZ ------------------------------
    KA_DTB=12.0, FM_DTB=0.50,
    CL_HTB=1200.0, V_HTB=800.0, Q_HTB=150.0, V3_HTB=500.0,
    EC50_VMAT_HTB=19.0,
    CYP2D6=1.0,                      # CL multiplier for NBI + HTBZ
    # --- amantadine / anticholinergic / ginkgo / botulinum ------------------
    KA_AMA=8.0, CL_AMA=400.0, V_AMA=350.0, EC50_AMA=400.0,
    KA_ACH=10.0, CL_ACH=1000.0, V_ACH=800.0, EC50_ACH=1.5,
    KIN_GKB=0.143, KOUT_GKB=0.143, ANTIOX_MAX=0.85,
    KOUT_BONT=0.0111, BONT_MAX=0.35,
    # --- dopamine handling -------------------------------------------------
    SYN0=1.0, A_AUTO=0.60, KVES=5.0, KMAO=3.0, KREL=2.0,
    KDAT=8.0, KDIFF=1.0, KBACK=7.2,
    A_FIRE=0.50, KIN_DEP=0.0167, KOUT_DEP=0.0167, DEPOL_MAX=0.45,
    DA_SYN0=0.13889,                 # analytic baseline synaptic DA
    DA_CYT0=0.25,                    # analytic baseline cytosolic DA
    # --- D2 supersensitivity (RUP) ----------------------------------------
    KIN_R=0.0035, KOUT_R=0.003333, OCC50_R=0.70, HILL_R=4.0,
    W_ACH_R=0.50, W_SDAM_R=0.80,
    # --- oxidative stress --------------------------------------------------
    KIN_ROS=0.03333, KOUT_ROS=0.03333,
    W_CYT_ROS=0.80, W_OCC_ROS=0.50, W_RUP_ROS=0.25, ROS_CONST=0.0,
    W_S_ROS=2.20, S_ROS50=0.50, HILL_S_ROS=4.0,
    RISK_FGA=1.0, W_AMA_ROS=0.30,
    # --- structural latch (SDAM) ------------------------------------------
    KIN_S=0.0025, KOUT_S=0.00030,
    W_ROS_S=0.55, W_RUP_S=0.45, S50=0.85, HILL_S=8.0,
    # --- postsynaptic D2 stimulation -> basal ganglia ----------------------
    G_RUP=2.20, MASK=0.55, G_ACH=0.50, G_AMA=0.35, G_S=1.20,
    K_IND=0.80, K_DIR=0.50,
    KBG=0.50, K_THAL_GPI=0.80,
    AIMS_MAX=26.0, TH50=0.55, HILL_TH=2.0, K_AIMS=0.1429,
    # --- psychiatric / motor / mood endpoints ------------------------------
    OCC50_P=0.58, HILL_P=3.0, E_MAX_P=1.25, W_SUPER_P=0.20, W_CLZ_P=0.75,
    K_PSY=0.0333,
    OCC50_PK=0.78, HILL_PK=4.0, VOCC50_PK=0.70, HILL_VPK=3.0,
    W_RUP_PK=0.30, K_PARK=0.0714,
    VOCC50_DEP=0.75, HILL_VDEP=2.0, W_VMAT_DEP=0.60, W_OCC_DEP=0.20, K_DEPR=0.0476,
    QT_NBI=0.080, QT_HTB=0.250, K_QTC=1.0,
    W_PARK_AD=0.50, W_DEPR_AD=0.40, W_PSY_AD=0.30, K_ADHER=0.0333,
    # --- clinician dose-escalation policy ---------------------------------
    ESC_ON=0.0, K_ESC=0.004, PSY_TARGET=0.20, DOSEMULT_MAX=2.0,
    # --- host risk modifiers ----------------------------------------------
    AGE=58.0, ESTROGEN=0.0, DM_RISK=0.0, GEN_RISK=1.0,
    # --- regimen -----------------------------------------------------------
    DOSE_AP=4.0, DOSE_AP2=4.0, TSW_AP=1e6,
    DOSE_CLZ=0.0, TSTART_CLZ=1e6,
    DOSE_VAL=0.0, TSTART_VAL=1e6, TSTOP_VAL=1e6,
    DOSE_DTB=0.0, TSTART_DTB=1e6, TSTOP_DTB=1e6,
    DOSE_AMA=0.0, TSTART_AMA=1e6,
    DOSE_ACH=0.0, TSTART_ACH=1e6,
    GKB_ON=0.0, TSTART_GKB=1e6,
    BONT_ON=0.0, TSTART_BONT=1e6, BONT_INT=90.0,
    TAPER=1.0,                       # days over which regimen changes ramp
)


def hill(x, x50, n):
    x = max(x, 0.0)
    xn = x ** n
    return xn / (xn + x50 ** n) if (xn + x50 ** n) > 0 else 0.0


def step(t, t0, tau=1.0):
    """smooth unit step at t0 (same form in the R model)"""
    return 0.5 * (1.0 + np.tanh((t - t0) / tau))


# ----------------------------------------------------------------------------
# right-hand side
# ----------------------------------------------------------------------------
def rhs(t, y, p):
    d = np.zeros(NS)
    g = lambda k: y[S[k]]

    # ---------------- regimen drivers -------------------------------------
    tau = p['TAPER']
    ap_dose = p['DOSE_AP'] + (p['DOSE_AP2'] - p['DOSE_AP']) * step(t, p['TSW_AP'], tau)
    ap_dose *= g('DOSEMULT')
    adh = min(max(g('ADHER'), 0.0), 1.0)
    ap_oral = ap_dose * (1.0 - p['F_LAI']) * adh
    ap_depot = ap_dose * p['F_LAI']              # LAI: adherence-independent
    clz_dose = p['DOSE_CLZ'] * step(t, p['TSTART_CLZ'], tau) * adh
    val_dose = p['DOSE_VAL'] * (step(t, p['TSTART_VAL'], tau) - step(t, p['TSTOP_VAL'], tau))
    dtb_dose = p['DOSE_DTB'] * (step(t, p['TSTART_DTB'], tau) - step(t, p['TSTOP_DTB'], tau))
    ama_dose = p['DOSE_AMA'] * step(t, p['TSTART_AMA'], tau)
    ach_dose = p['DOSE_ACH'] * step(t, p['TSTART_ACH'], tau)
    gkb_in = p['GKB_ON'] * step(t, p['TSTART_GKB'], tau)
    bont_in = p['BONT_ON'] * step(t, p['TSTART_BONT'], tau)

    # ---------------- PK ---------------------------------------------------
    C_AP = g('AP_cen') / p['V2_AP'] * 1000.0        # mg/L -> ng/mL
    C_APp = g('AP_per') / p['V3_AP'] * 1000.0
    d[S['AP_gut']] = ap_oral - p['KA_AP'] * g('AP_gut')
    d[S['AP_lai']] = ap_depot - p['KA_LAI'] * g('AP_lai')
    d[S['AP_cen']] = (p['KA_AP'] * g('AP_gut') + p['KA_LAI'] * g('AP_lai')
                      - p['CL_AP'] * C_AP / 1000.0
                      - p['Q_AP'] * (C_AP - C_APp) / 1000.0)
    d[S['AP_per']] = p['Q_AP'] * (C_AP - C_APp) / 1000.0

    C_CLZ = g('CLZ_cen') / p['V_CLZ'] * 1000.0
    d[S['CLZ_gut']] = clz_dose - p['KA_CLZ'] * g('CLZ_gut')
    d[S['CLZ_cen']] = p['KA_CLZ'] * g('CLZ_gut') - p['CL_CLZ'] * C_CLZ / 1000.0

    C_VAL = g('VAL_cen') / p['V_VAL'] * 1000.0
    C_NBI = g('NBI_cen') / p['V_NBI'] * 1000.0
    C_NBIp = g('NBI_per') / p['V3_NBI'] * 1000.0
    cl_nbi = p['CL_NBI'] * p['CYP2D6']
    d[S['VAL_gut']] = val_dose - p['KA_VAL'] * g('VAL_gut')
    d[S['VAL_cen']] = p['KA_VAL'] * g('VAL_gut') - p['CL_VAL'] * C_VAL / 1000.0
    d[S['NBI_cen']] = (p['FM_VAL'] * p['CL_VAL'] * C_VAL / 1000.0
                       - cl_nbi * C_NBI / 1000.0
                       - p['Q_NBI'] * (C_NBI - C_NBIp) / 1000.0)
    d[S['NBI_per']] = p['Q_NBI'] * (C_NBI - C_NBIp) / 1000.0

    C_HTB = g('HTB_cen') / p['V_HTB'] * 1000.0
    C_HTBp = g('HTB_per') / p['V3_HTB'] * 1000.0
    cl_htb = p['CL_HTB'] * p['CYP2D6']
    d[S['DTB_gut']] = dtb_dose - p['KA_DTB'] * g('DTB_gut')
    d[S['HTB_cen']] = (p['FM_DTB'] * p['KA_DTB'] * g('DTB_gut')
                       - cl_htb * C_HTB / 1000.0
                       - p['Q_HTB'] * (C_HTB - C_HTBp) / 1000.0)
    d[S['HTB_per']] = p['Q_HTB'] * (C_HTB - C_HTBp) / 1000.0

    C_AMA = g('AMA_cen') / p['V_AMA'] * 1000.0
    d[S['AMA_gut']] = ama_dose - p['KA_AMA'] * g('AMA_gut')
    d[S['AMA_cen']] = p['KA_AMA'] * g('AMA_gut') - p['CL_AMA'] * C_AMA / 1000.0

    C_ACH = g('ACH_cen') / p['V_ACH'] * 1000.0
    d[S['ACH_gut']] = ach_dose - p['KA_ACH'] * g('ACH_gut')
    d[S['ACH_cen']] = p['KA_ACH'] * g('ACH_gut') - p['CL_ACH'] * C_ACH / 1000.0

    d[S['GKB']] = p['KIN_GKB'] * gkb_in - p['KOUT_GKB'] * g('GKB')
    d[S['BONT']] = bont_in / p['BONT_INT'] - p['KOUT_BONT'] * g('BONT')

    # ---------------- receptor occupancies --------------------------------
    d[S['CE_AP']] = p['KE0_AP'] * (C_AP - g('CE_AP'))
    occ_ap = hill(g('CE_AP'), p['EC50_D2'], p['HILL_D2'])
    occ_clz = hill(C_CLZ, p['EC50_D2_CLZ'], 1.0)
    OCC = min(occ_ap + occ_clz * (1.0 - occ_ap), 0.995)      # combined D2 block
    EFF_CLZ = hill(C_CLZ, p['EC50_CLZ_EFF'], 2.0)
    OCCV = min(hill(C_NBI, p['EC50_VMAT_NBI'], 1.0)
               + hill(C_HTB, p['EC50_VMAT_HTB'], 1.0) * (1 - hill(C_NBI, p['EC50_VMAT_NBI'], 1.0)),
               0.98)
    E_AMA = hill(C_AMA, p['EC50_AMA'], 1.0)
    E_ACH = hill(C_ACH, p['EC50_ACH'], 1.0)
    ANTIOX = p['ANTIOX_MAX'] * min(g('GKB'), 1.0)
    E_BONT = p['BONT_MAX'] * min(g('BONT'), 1.0)

    # ---------------- dopamine pools --------------------------------------
    syn = p['SYN0'] * (1.0 + p['A_AUTO'] * OCC)
    fire = max(1.0 + p['A_FIRE'] * OCC - g('DEPOL'), 0.05)
    kves = p['KVES'] * (1.0 - OCCV)
    d[S['DA_CYT']] = (syn + p['KBACK'] * g('DA_SYN')
                      - (kves + p['KMAO']) * g('DA_CYT'))
    d[S['DA_VES']] = kves * g('DA_CYT') - p['KREL'] * fire * g('DA_VES')
    d[S['DA_SYN']] = (p['KREL'] * fire * g('DA_VES')
                      - (p['KDAT'] + p['KDIFF']) * g('DA_SYN'))
    d[S['DEPOL']] = (p['KIN_DEP'] * p['DEPOL_MAX'] * hill(OCC, 0.75, 4.0)
                     - p['KOUT_DEP'] * g('DEPOL'))

    DA_N = g('DA_SYN') / p['DA_SYN0']
    CYT_N = g('DA_CYT') / p['DA_CYT0']

    # ---------------- plasticity: supersensitivity, ROS, structural latch --
    risk_age = 1.0 + 0.015 * max(p['AGE'] - 40.0, 0.0)
    RISKMOD = risk_age * p['GEN_RISK'] * (1.0 + 0.25 * p['DM_RISK']) \
        * (1.0 - 0.25 * p['ESTROGEN'])
    d[S['RUP']] = (p['KIN_R'] * hill(OCC, p['OCC50_R'], p['HILL_R'])
                   * RISKMOD * (1.0 + p['W_ACH_R'] * E_ACH)
                   - p['KOUT_R'] * g('RUP') / (1.0 + p['W_SDAM_R'] * g('SDAM')))

    ros_drive = (p['W_CYT_ROS'] * max(CYT_N ** 1.5 - 1.0, 0.0)
                 + p['W_OCC_ROS'] * hill(OCC, 0.70, 4.0) * p['RISK_FGA']
                 + p['W_RUP_ROS'] * g('RUP')
                 + p['W_S_ROS'] * hill(y[S['SDAM']], p['S_ROS50'],
                                      p['HILL_S_ROS'])
                 + 0.25 * max(p['AGE'] - 40.0, 0.0) / 30.0
                 + 0.20 * p['DM_RISK'] + p['ROS_CONST'])
    d[S['ROS']] = (p['KIN_ROS'] * ros_drive
                   - p['KOUT_ROS'] * (1.0 + ANTIOX + p['W_AMA_ROS'] * E_AMA)
                   * g('ROS'))

    drive_s = p['W_ROS_S'] * g('ROS') + p['W_RUP_S'] * g('RUP')
    d[S['SDAM']] = (p['KIN_S'] * hill(drive_s, p['S50'], p['HILL_S'])
                    * (1.0 - g('SDAM')) - p['KOUT_S'] * g('SDAM'))

    # ---------------- postsynaptic D2 stimulation -------------------------
    D2STIM = (DA_N * (1.0 + p['G_RUP'] * g('RUP')) * (1.0 - p['MASK'] * OCC)
              * (1.0 + p['G_ACH'] * E_ACH) / (1.0 + p['G_AMA'] * E_AMA))
    EXC = max(D2STIM - 1.0, 0.0) + p['G_S'] * g('SDAM')

    # ---------------- basal-ganglia loop ---------------------------------
    ind_t = 1.0 / (1.0 + p['K_IND'] * EXC)
    d[S['IND']] = p['KBG'] * (ind_t - g('IND'))
    gpe_t = 2.0 - g('IND')
    d[S['GPE']] = p['KBG'] * (gpe_t - g('GPE'))
    stn_t = 1.0 / (0.5 + 0.5 * g('GPE'))
    d[S['STN']] = p['KBG'] * (stn_t - g('STN'))
    gpi_t = 0.55 * g('STN') + 0.45 / (1.0 + p['K_DIR'] * EXC)
    d[S['GPI']] = p['KBG'] * (gpi_t - g('GPI'))
    thal_t = 1.0 / (1.0 - p['K_THAL_GPI'] + p['K_THAL_GPI'] * g('GPI'))
    d[S['THAL']] = p['KBG'] * (thal_t - g('THAL'))

    aims_t = (p['AIMS_MAX'] * hill(g('THAL') - 1.0, p['TH50'], p['HILL_TH'])
              * (1.0 - E_BONT))
    d[S['AIMS']] = p['K_AIMS'] * (aims_t - g('AIMS'))

    # ---------------- other endpoints -------------------------------------
    control = (p['E_MAX_P'] * hill(OCC, p['OCC50_P'], p['HILL_P'])
               / (1.0 + p['W_SUPER_P'] * g('RUP') ** 1.5)
               + p['W_CLZ_P'] * EFF_CLZ)
    psy_t = max(1.0 - min(control, 1.0), 0.0)
    d[S['PSYCH']] = p['K_PSY'] * (psy_t - g('PSYCH'))

    park_t = max(1.2 * hill(OCC, p['OCC50_PK'], p['HILL_PK'])
                 + 0.9 * hill(OCCV, p['VOCC50_PK'], p['HILL_VPK'])
                 - p['W_RUP_PK'] * g('RUP'), 0.0)
    d[S['PARK']] = p['K_PARK'] * (park_t - g('PARK'))

    depr_t = (p['W_VMAT_DEP'] * hill(OCCV, p['VOCC50_DEP'], p['HILL_VDEP'])
              + p['W_OCC_DEP'] * hill(OCC, 0.80, 4.0))
    d[S['DEPR']] = p['K_DEPR'] * (depr_t - g('DEPR'))

    qtc_t = p['QT_NBI'] * C_NBI + p['QT_HTB'] * C_HTB
    d[S['QTC']] = p['K_QTC'] * (qtc_t - g('QTC'))

    adh_t = 1.0 / (1.0 + p['W_PARK_AD'] * g('PARK') + p['W_DEPR_AD'] * g('DEPR')
                   + p['W_PSY_AD'] * g('PSYCH'))
    d[S['ADHER']] = p['K_ADHER'] * (adh_t - g('ADHER'))

    esc = p['ESC_ON'] * p['K_ESC'] * max(g('PSYCH') - p['PSY_TARGET'], 0.0)
    d[S['DOSEMULT']] = esc * (p['DOSEMULT_MAX'] - g('DOSEMULT'))

    d[S['CUMOCC']] = OCC
    return d


def y0(p):
    y = np.zeros(NS)
    y[S['DA_CYT']] = p['DA_CYT0']
    y[S['DA_VES']] = 0.625
    y[S['DA_SYN']] = p['DA_SYN0']
    for k in ('IND', 'GPE', 'STN', 'GPI', 'THAL', 'ADHER', 'DOSEMULT'):
        y[S[k]] = 1.0
    return y


def run(days=1825, out=None, **kw):
    p = dict(P0)
    p.update(kw)
    t_eval = np.arange(0, days + 1, 1.0) if out is None else np.asarray(out, float)
    sol = solve_ivp(rhs, (0.0, float(days)), y0(p), args=(p,), method='LSODA',
                    rtol=1e-6, atol=1e-9, t_eval=t_eval, max_step=8.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y, p


def derived(t, Y, p):
    """recompute the algebraic read-outs along a solution"""
    C_AP = Y[S['AP_cen']] / p['V2_AP'] * 1000.0
    C_NBI = Y[S['NBI_cen']] / p['V_NBI'] * 1000.0
    C_HTB = Y[S['HTB_cen']] / p['V_HTB'] * 1000.0
    C_CLZ = Y[S['CLZ_cen']] / p['V_CLZ'] * 1000.0
    oa = np.array([hill(c, p['EC50_D2'], p['HILL_D2']) for c in Y[S['CE_AP']]])
    oc = np.array([hill(c, p['EC50_D2_CLZ'], 1.0) for c in C_CLZ])
    OCC = np.minimum(oa + oc * (1 - oa), 0.995)
    on = np.array([hill(c, p['EC50_VMAT_NBI'], 1.0) for c in C_NBI])
    oh = np.array([hill(c, p['EC50_VMAT_HTB'], 1.0) for c in C_HTB])
    OCCV = np.minimum(on + oh * (1 - on), 0.98)
    return dict(C_AP=C_AP, C_NBI=C_NBI, C_HTB=C_HTB, C_CLZ=C_CLZ,
                OCC=OCC, OCCV=OCCV, DA_N=Y[S['DA_SYN']] / p['DA_SYN0'],
                AIMS=Y[S['AIMS']], RUP=Y[S['RUP']], SDAM=Y[S['SDAM']],
                ROS=Y[S['ROS']], PSYCH=Y[S['PSYCH']], PARK=Y[S['PARK']],
                DEPR=Y[S['DEPR']], QTC=Y[S['QTC']], ADHER=Y[S['ADHER']],
                CUMOCC=Y[S['CUMOCC']], DOSEMULT=Y[S['DOSEMULT']])


def at(t, v, day):
    return float(np.interp(day, t, v))


# =============================================================================
# 1. scenarios
# =============================================================================
SCEN = {
    'S1 natural history (risp 4 mg, 5 y)':
        dict(),
    'S2 FGA-equivalent high occupancy (risp-eq 8 mg)':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6),
    'S3 50% dose reduction at y2':
        dict(DOSE_AP=8.0, DOSE_AP2=4.0, TSW_AP=730.0, RISK_FGA=1.6),
    'S4 full withdrawal at y2':
        dict(DOSE_AP=8.0, DOSE_AP2=0.0, TSW_AP=730.0, RISK_FGA=1.6),
    'S5 switch to clozapine at y2':
        dict(DOSE_AP=8.0, DOSE_AP2=0.0, TSW_AP=730.0, RISK_FGA=1.6,
             DOSE_CLZ=350.0, TSTART_CLZ=730.0),
    'S6 valbenazine 80 mg from y2':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6,
             DOSE_VAL=80.0, TSTART_VAL=730.0),
    'S7 deutetrabenazine 36 mg from y2':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6,
             DOSE_DTB=36.0, TSTART_DTB=730.0),
    'S8 benztropine add-on from y2':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6,
             DOSE_ACH=2.0, TSTART_ACH=730.0),
    'S9 ginkgo+amantadine from y2':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6,
             GKB_ON=1.0, TSTART_GKB=730.0, DOSE_AMA=200.0, TSTART_AMA=730.0),
    'S10 clinician escalation policy':
        dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6, ESC_ON=1.0),
}


def scenarios():
    print('\n' + '=' * 108)
    print('1. TREATMENT SCENARIOS  (5-year simulation; AIMS = dyskinesia '
          'subscore 0-26)')
    print('=' * 108)
    print(f"{'scenario':46s}{'AIMS':>7s}{'AIMS':>7s}{'RUP':>7s}"
          f"{'SDAM':>7s}{'OCC':>7s}{'PSY':>7s}{'PARK':>7s}{'DEPR':>7s}"
          f"{'QTc':>7s}")
    print(f"{'':46s}{'y2':>7s}{'y5':>7s}{'y5':>7s}{'y5':>7s}{'y5':>7s}"
          f"{'y5':>7s}{'y5':>7s}{'y5':>7s}{'y5':>7s}")
    print('-' * 108)
    res = {}
    for name, kw in SCEN.items():
        t, Y, p = run(1825, **kw)
        D = derived(t, Y, p)
        res[name] = (t, D)
        print(f"{name:46s}{at(t,D['AIMS'],730):7.2f}{at(t,D['AIMS'],1825):7.2f}"
              f"{at(t,D['RUP'],1825):7.3f}{at(t,D['SDAM'],1825):7.3f}"
              f"{at(t,D['OCC'],1825):7.3f}{at(t,D['PSYCH'],1825):7.3f}"
              f"{at(t,D['PARK'],1825):7.3f}{at(t,D['DEPR'],1825):7.3f}"
              f"{at(t,D['QTC'],1825):7.2f}")
    return res


# =============================================================================
# 2. cumulative-exposure threshold: dose x duration at matched exposure
# =============================================================================
def exposure_threshold():
    print('\n' + '=' * 108)
    print('2. IS IT DOSE OR IS IT TIME?  matched cumulative occupancy-days, '
          'different dose/duration splits')
    print('=' * 108)
    print(f"{'dose (risp-eq mg/d)':>22s}{'duration (d)':>14s}"
          f"{'OCC':>8s}{'occ-days':>10s}{'drive-days':>12s}{'SDAM end':>10s}"
          f"{'AIMS end':>10s}{'AIMS stop+6y':>13s}{'outcome':>14s}")
    print('-' * 108)
    rows = []
    for dose, dur in [(16.0, 300), (8.0, 460), (6.0, 560), (4.0, 700),
                      (3.0, 900), (2.0, 1500), (1.5, 1825)]:
        t, Y, p = run(dur + 2190, DOSE_AP=dose, DOSE_AP2=0.0,
                      TSW_AP=float(dur), RISK_FGA=1.6)
        D = derived(t, Y, p)
        occ = at(t, D['OCC'], dur - 5)
        occd = at(t, D['CUMOCC'], dur)
        m = t <= dur
        drv = float(np.trapezoid([hill(o, p['OCC50_R'], p['HILL_R'])
                                  for o in D['OCC'][m]], t[m]))
        sd_end = at(t, D['SDAM'], dur)
        aims_end = at(t, D['AIMS'], dur)
        aims5 = at(t, D['AIMS'], dur + 2190)
        out = 'PERSISTENT' if aims5 > 2.0 else 'resolved'
        rows.append((dose, dur, occ, occd, drv, sd_end, aims_end, aims5, out))
        print(f"{dose:22.1f}{dur:14d}{occ:8.3f}{occd:10.1f}{drv:12.1f}"
              f"{sd_end:10.3f}{aims_end:10.2f}{aims5:13.2f}{out:>14s}")
    return rows


# =============================================================================
# 3. reversibility window: when does withdrawal stop working?
# =============================================================================
def reversibility_window():
    print('\n' + '=' * 108)
    print('3. REVERSIBILITY WINDOW — withdrawal at increasing exposure '
          'duration, all followed 6 y past stop')
    print('=' * 108)
    print(f"{'stop day':>10s}{'SDAM@stop':>11s}{'RUP@stop':>10s}"
          f"{'AIMS@stop':>11s}{'AIMS peak':>11s}{'peak day':>10s}"
          f"{'AIMS +2y':>10s}{'AIMS +6y':>10s}{'outcome':>13s}")
    print('-' * 108)
    rows = []
    for stop in [60, 90, 120, 150, 180, 210, 240, 260, 270, 275, 280, 285,
                 290, 300, 330, 365, 460, 550, 730, 1095, 1460]:
        t, Y, p = run(stop + 2190, DOSE_AP=8.0, DOSE_AP2=0.0,
                      TSW_AP=float(stop), RISK_FGA=1.6)
        D = derived(t, Y, p)
        m = t >= stop
        ip = int(np.argmax(D['AIMS'][m]))
        rows.append((stop, at(t, D['SDAM'], stop), at(t, D['RUP'], stop),
                     at(t, D['AIMS'], stop), float(D['AIMS'][m][ip]),
                     float(t[m][ip] - stop), at(t, D['AIMS'], stop + 730),
                     at(t, D['AIMS'], stop + 2190)))
        r = rows[-1]
        out = 'PERSISTENT' if r[7] > 2.0 else 'resolved'
        print(f"{r[0]:10d}{r[1]:11.3f}{r[2]:10.3f}{r[3]:11.2f}{r[4]:11.2f}"
              f"{r[5]:10.0f}{r[6]:10.2f}{r[7]:10.2f}{out:>13s}")
    return rows


# =============================================================================
# 4. the withdrawal paradox and its crossover time
# =============================================================================
def withdrawal_crossover():
    print('\n' + '=' * 108)
    print('4. WITHDRAWAL PARADOX — AIMS after dose change at day 730 '
          '(reference = continue unchanged)')
    print('=' * 108)
    print(f"{'strategy':30s}{'AIMS d730':>10s}{'d744':>8s}{'d760':>8s}"
          f"{'d820':>8s}{'d1095':>8s}{'d1825':>8s}{'d2555':>8s}"
          f"{'peak':>8s}{'crossover d':>13s}")
    print('-' * 108)
    tR, YR, pR = run(2555, DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
    DR = derived(tR, YR, pR)
    rows = []
    for name, kw in [
        ('continue 8 mg', dict(DOSE_AP2=8.0)),
        ('reduce to 6 mg', dict(DOSE_AP2=6.0)),
        ('reduce to 4 mg', dict(DOSE_AP2=4.0)),
        ('reduce to 2 mg', dict(DOSE_AP2=2.0)),
        ('stop', dict(DOSE_AP2=0.0)),
        ('stop + clozapine 350', dict(DOSE_AP2=0.0, DOSE_CLZ=350.0,
                                      TSTART_CLZ=730.0)),
    ]:
        t, Y, p = run(2555, DOSE_AP=8.0, TSW_AP=730.0, RISK_FGA=1.6, **kw)
        D = derived(t, Y, p)
        m = t >= 730
        dif = D['AIMS'][m] - DR['AIMS'][m]
        cross = next((float(t[m][i] - 730) for i in range(len(dif))
                      if dif[i] < 0), float('nan'))
        rows.append((name, [at(t, D['AIMS'], d) for d in
                            (730, 744, 760, 820, 1095, 1825, 2555)],
                     float(D['AIMS'][m].max()), cross))
        r = rows[-1]
        print(f"{r[0]:30s}" + ''.join(f"{v:8.2f}" for v in
              [r[1][0]] + r[1][1:]) + f"{r[2]:8.2f}{r[3]:13.0f}")
    return rows


# =============================================================================
# 5. VMAT2 dose-response and the CYP2D6 axis
# =============================================================================
def vmat2_dose_response():
    print('\n' + '=' * 108)
    print('5. VMAT2 INHIBITOR DOSE-RESPONSE (added at y2 to ongoing 8 mg; '
          'read out at week 6 = day 772)')
    print('=' * 108)
    print(f"{'drug / dose':28s}{'Cave':>8s}{'VMAT2occ':>10s}{'DA_syn':>8s}"
          f"{'AIMS':>8s}{'dAIMS':>8s}{'%':>7s}{'PARK':>7s}{'DEPR':>7s}"
          f"{'QTc':>7s}{'TI':>7s}{'ROS y5':>8s}{'SDAM y5':>9s}{'AIMS y5':>9s}")
    print('-' * 108)
    t0, Y0_, p0 = run(1825, DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
    D0 = derived(t0, Y0_, p0)
    base = at(t0, D0['AIMS'], 730)
    ref772 = at(t0, D0['AIMS'], 772)
    ref_y5 = (at(t0, D0['ROS'], 1825), at(t0, D0['SDAM'], 1825),
              at(t0, D0['AIMS'], 1825))
    rows = []
    grid = ([('valbenazine', d, dict(DOSE_VAL=float(d), TSTART_VAL=730.0))
             for d in (20, 40, 60, 80, 120)]
            + [('deutetrabenazine', d, dict(DOSE_DTB=float(d),
                                            TSTART_DTB=730.0))
               for d in (12, 24, 36, 48, 72)])
    for drug, dose, kw in grid:
        t, Y, p = run(1825, DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6, **kw)
        D = derived(t, Y, p)
        c = at(t, D['C_NBI'] if drug == 'valbenazine' else D['C_HTB'], 772)
        a = at(t, D['AIMS'], 772)
        da = a - base
        pk = at(t, D['PARK'], 772)
        dp = at(t, D['DEPR'], 772)
        ti = (-da) / max(pk + dp, 1e-6)
        rows.append((drug, dose, c, at(t, D['OCCV'], 772),
                     at(t, D['DA_N'], 772), a, da, 100 * da / base, pk, dp,
                     at(t, D['QTC'], 772), ti, at(t, D['ROS'], 1825),
                     at(t, D['SDAM'], 1825), at(t, D['AIMS'], 1825)))
        r = rows[-1]
        print(f"{drug + ' ' + str(dose) + ' mg':28s}{r[2]:8.1f}{r[3]:10.3f}"
              f"{r[4]:8.3f}{r[5]:8.2f}{r[6]:8.2f}{r[7]:7.1f}{r[8]:7.3f}"
              f"{r[9]:7.3f}{r[10]:7.2f}{r[11]:7.2f}{r[12]:8.3f}{r[13]:9.3f}"
              f"{r[14]:9.2f}")
    print(f"  reference (no VMAT2 inhibitor): AIMS day730 = {base:.2f}, "
          f"day772 = {ref772:.2f}, ROS y5 = {ref_y5[0]:.3f}, "
          f"SDAM y5 = {ref_y5[1]:.3f}, AIMS y5 = {ref_y5[2]:.2f}")
    return rows


def cyp2d6_panel():
    print('\n' + '=' * 108)
    print('6. CYP2D6 PHENOTYPE x VMAT2 INHIBITOR (day 772 read-out, '
          'label doses)')
    print('=' * 108)
    print(f"{'phenotype (CL mult)':26s}{'drug':20s}{'Cave':>8s}"
          f"{'VMAT2occ':>10s}{'AIMS':>8s}{'dAIMS%':>8s}{'PARK':>7s}"
          f"{'DEPR':>7s}{'QTc':>7s}")
    print('-' * 108)
    t0, Y0_, p0 = run(800, DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
    base = at(t0, derived(t0, Y0_, p0)['AIMS'], 730)
    rows = []
    for ph, mult in [('UM (1.6x)', 1.6), ('normal (1.0x)', 1.0),
                     ('IM (0.7x)', 0.7), ('PM (0.5x)', 0.5),
                     ('PM + CYP3A4 inh (0.35x)', 0.35)]:
        for drug, kw in [('valbenazine 80', dict(DOSE_VAL=80.0,
                                                 TSTART_VAL=730.0)),
                         ('deutetrab. 36', dict(DOSE_DTB=36.0,
                                                TSTART_DTB=730.0))]:
            t, Y, p = run(800, DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6,
                          CYP2D6=mult, **kw)
            D = derived(t, Y, p)
            c = at(t, D['C_NBI'] if 'val' in drug else D['C_HTB'], 772)
            a = at(t, D['AIMS'], 772)
            rows.append((ph, drug, c, at(t, D['OCCV'], 772), a,
                         100 * (a - base) / base, at(t, D['PARK'], 772),
                         at(t, D['DEPR'], 772), at(t, D['QTC'], 772)))
            r = rows[-1]
            print(f"{ph:26s}{drug:20s}{r[2]:8.1f}{r[3]:10.3f}{r[4]:8.2f}"
                  f"{r[5]:8.1f}{r[6]:7.3f}{r[7]:7.3f}{r[8]:7.2f}")
    return rows


# =============================================================================
# 7. suppression vs. disease modification, and the washout test
# =============================================================================
def suppression_vs_modification():
    print('\n' + '=' * 108)
    print('7. SUPPRESSION vs DISEASE MODIFICATION — 2 y of therapy from y2, '
          'then washout at y4')
    print('=' * 108)
    print(f"{'arm':34s}{'AIMS y2':>9s}{'AIMS y4':>9s}{'RUP y4':>8s}"
          f"{'SDAM y4':>9s}{'AIMS y4+8w':>11s}{'AIMS y6':>9s}"
          f"{'rebound vs y2':>14s}")
    print('-' * 108)
    arms = {
        'no treatment': dict(),
        'valbenazine 80 (stop y4)': dict(DOSE_VAL=80.0, TSTART_VAL=730.0,
                                         TSTOP_VAL=1460.0),
        'deutetrab. 36 (stop y4)': dict(DOSE_DTB=36.0, TSTART_DTB=730.0,
                                        TSTOP_DTB=1460.0),
        'clozapine switch y2': dict(DOSE_AP2=0.0, TSW_AP=730.0,
                                    DOSE_CLZ=350.0, TSTART_CLZ=730.0),
        'clozapine + valbenaz. (stop y4)': dict(DOSE_AP2=0.0, TSW_AP=730.0,
                                                DOSE_CLZ=350.0,
                                                TSTART_CLZ=730.0,
                                                DOSE_VAL=80.0,
                                                TSTART_VAL=730.0,
                                                TSTOP_VAL=1460.0),
    }
    rows = []
    for name, kw in arms.items():
        base_kw = dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
        base_kw.update(kw)
        t, Y, p = run(2190, **base_kw)
        D = derived(t, Y, p)
        rows.append((name, at(t, D['AIMS'], 730), at(t, D['AIMS'], 1460),
                     at(t, D['RUP'], 1460), at(t, D['SDAM'], 1460),
                     at(t, D['AIMS'], 1516), at(t, D['AIMS'], 2190),
                     at(t, D['AIMS'], 1516) - at(t, D['AIMS'], 730)))
        r = rows[-1]
        print(f"{r[0]:34s}{r[1]:9.2f}{r[2]:9.2f}{r[3]:8.3f}{r[4]:9.3f}"
              f"{r[5]:11.2f}{r[6]:9.2f}{r[7]:14.2f}")
    return rows


# =============================================================================
# 8. combination interaction: does suppression help modification?
# =============================================================================
def combination_interaction():
    print('\n' + '=' * 108)
    print('8. INTERACTION — clozapine switch x valbenazine (2x2, effect = '
          'AIMS change from day 730)')
    print('=' * 108)
    arms = {}
    for clz in (0, 1):
        for val in (0, 1):
            kw = dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
            if clz:
                kw.update(DOSE_AP2=0.0, TSW_AP=730.0, DOSE_CLZ=350.0,
                          TSTART_CLZ=730.0)
            if val:
                kw.update(DOSE_VAL=80.0, TSTART_VAL=730.0)
            t, Y, p = run(2555, **kw)
            D = derived(t, Y, p)
            arms[(clz, val)] = {d: at(t, D['AIMS'], d)
                                for d in (730, 772, 1095, 1825, 2555)}
    print(f"{'endpoint':12s}{'neither':>10s}{'valben':>10s}{'clozap':>10s}"
          f"{'both':>10s}{'E_val':>9s}{'E_clz':>9s}{'interaction':>13s}")
    print('-' * 108)
    rows = []
    for d in (772, 1095, 1825, 2555):
        b = arms[(0, 0)][730]
        e00 = arms[(0, 0)][d] - b
        e01 = arms[(0, 1)][d] - b
        e10 = arms[(1, 0)][d] - b
        e11 = arms[(1, 1)][d] - b
        inter = e11 - (e01 + e10 - e00)
        rows.append((d, e00, e01, e10, e11, e01 - e00, e10 - e00, inter))
        print(f"{'day ' + str(d):12s}{e00:10.2f}{e01:10.2f}{e10:10.2f}"
              f"{e11:10.2f}{e01 - e00:9.2f}{e10 - e00:9.2f}{inter:13.2f}")
    return rows


# =============================================================================
# 9. the opposed-levers frontier: AIMS vs psychosis control
# =============================================================================
def opposed_levers():
    print('\n' + '=' * 108)
    print('9. OPPOSED LEVERS — the two ways to lower AIMS trade against '
          'different things (day 1095)')
    print('=' * 108)
    print(f"{'strategy from y2':34s}{'OCC':>7s}{'VMAT2':>7s}{'AIMS':>8s}"
          f"{'PSYCH':>8s}{'PARK':>7s}{'DEPR':>7s}{'ADHER':>7s}{'RUP':>7s}"
          f"{'SDAM':>7s}")
    print('-' * 108)
    strat = [
        ('continue 8 mg', dict()),
        ('lower D2 block: 4 mg', dict(DOSE_AP2=4.0, TSW_AP=730.0)),
        ('lower D2 block: 2 mg', dict(DOSE_AP2=2.0, TSW_AP=730.0)),
        ('lower D2 block: stop', dict(DOSE_AP2=0.0, TSW_AP=730.0)),
        ('clozapine 350 switch', dict(DOSE_AP2=0.0, TSW_AP=730.0,
                                      DOSE_CLZ=350.0, TSTART_CLZ=730.0)),
        ('lower DA supply: VBZ 40', dict(DOSE_VAL=40.0, TSTART_VAL=730.0)),
        ('lower DA supply: VBZ 80', dict(DOSE_VAL=80.0, TSTART_VAL=730.0)),
        ('lower DA supply: VBZ 120', dict(DOSE_VAL=120.0, TSTART_VAL=730.0)),
        ('both: 4 mg + VBZ 80', dict(DOSE_AP2=4.0, TSW_AP=730.0,
                                     DOSE_VAL=80.0, TSTART_VAL=730.0)),
    ]
    rows = []
    for name, kw in strat:
        base = dict(DOSE_AP=8.0, DOSE_AP2=8.0, RISK_FGA=1.6)
        base.update(kw)
        t, Y, p = run(1200, **base)
        D = derived(t, Y, p)
        rows.append((name, at(t, D['OCC'], 1095), at(t, D['OCCV'], 1095),
                     at(t, D['AIMS'], 1095), at(t, D['PSYCH'], 1095),
                     at(t, D['PARK'], 1095), at(t, D['DEPR'], 1095),
                     at(t, D['ADHER'], 1095), at(t, D['RUP'], 1095),
                     at(t, D['SDAM'], 1095)))
        r = rows[-1]
        print(f"{r[0]:34s}{r[1]:7.3f}{r[2]:7.3f}{r[3]:8.2f}{r[4]:8.3f}"
              f"{r[5]:7.3f}{r[6]:7.3f}{r[7]:7.3f}{r[8]:7.3f}{r[9]:7.3f}")
    return rows


# =============================================================================
# 10. host risk-factor scan
# =============================================================================
def risk_scan():
    print('\n' + '=' * 108)
    print('10. HOST RISK MODIFIERS — identical regimen (8 mg risp-eq, 5 y), '
          'different host')
    print('=' * 108)
    print(f"{'host':44s}{'AIMS y1':>9s}{'AIMS y3':>9s}{'AIMS y5':>9s}"
          f"{'RUP y5':>8s}{'SDAM y5':>9s}{'latch day':>11s}")
    print('-' * 108)
    hosts = [
        ('25 y, SGA', dict(AGE=25.0, RISK_FGA=1.0)),
        ('45 y, SGA', dict(AGE=45.0, RISK_FGA=1.0)),
        ('58 y, SGA (index patient, SGA)', dict(AGE=58.0, RISK_FGA=1.0)),
        ('58 y, FGA (index patient)', dict(AGE=58.0, RISK_FGA=1.6)),
        ('25 y, FGA', dict(AGE=25.0, RISK_FGA=1.6)),
        ('75 y, FGA', dict(AGE=75.0, RISK_FGA=1.6)),
        ('58 y, FGA, diabetes', dict(AGE=58.0, RISK_FGA=1.6, DM_RISK=1.0)),
        ('58 y, FGA, estrogen-replete', dict(AGE=58.0, RISK_FGA=1.6,
                                             ESTROGEN=1.0)),
        ('58 y, FGA, high-risk genotype', dict(AGE=58.0, RISK_FGA=1.6,
                                               GEN_RISK=1.3)),
        ('58 y, FGA + benztropine 2 mg', dict(AGE=58.0, RISK_FGA=1.6,
                                              DOSE_ACH=2.0, TSTART_ACH=0.0)),
        ('58 y, FGA + ginkgo EGb761', dict(AGE=58.0, RISK_FGA=1.6,
                                          GKB_ON=1.0, TSTART_GKB=0.0)),
    ]
    rows = []
    for name, kw in hosts:
        base = dict(DOSE_AP=8.0, DOSE_AP2=8.0)
        base.update(kw)
        t, Y, p = run(1825, **base)
        D = derived(t, Y, p)
        idx = np.where(D['SDAM'] > 0.49)[0]
        latch = float(t[idx[0]]) if len(idx) else float('nan')
        rows.append((name, at(t, D['AIMS'], 365), at(t, D['AIMS'], 1095),
                     at(t, D['AIMS'], 1825), at(t, D['RUP'], 1825),
                     at(t, D['SDAM'], 1825), latch))
        r = rows[-1]
        print(f"{r[0]:44s}{r[1]:9.2f}{r[2]:9.2f}{r[3]:9.2f}{r[4]:8.3f}"
              f"{r[5]:9.3f}{r[6]:11.0f}")
    return rows


# =============================================================================
# 11. bistability of the structural latch (SDAM nullcline)
# =============================================================================
def latch_bistability():
    print('\n' + '=' * 108)
    print('11. STRUCTURAL LATCH — dSDAM/dt with NO external drive '
          '(drug-free), i.e. self-sustaining or not')
    print('=' * 108)
    p = dict(P0)
    print(f"{'SDAM':>8s}{'drive':>9s}{'Hill':>9s}{'dSDAM/dt (1/d)':>17s}"
          f"{'direction':>12s}")
    print('-' * 108)
    rows = []
    for s in np.arange(0.0, 1.0001, 0.05):
        # drug-free: RUP -> 0, so ROS is sustained only by the damage feedback
        ros = p['W_S_ROS'] * hill(s, p['S_ROS50'], p['HILL_S_ROS'])
        drive = p['W_ROS_S'] * ros
        h = hill(drive, p['S50'], p['HILL_S'])
        ds = p['KIN_S'] * h * (1 - s) - p['KOUT_S'] * s
        rows.append((s, drive, h, ds))
        print(f"{s:8.2f}{drive:9.3f}{h:9.4f}{ds:17.3e}"
              f"{('up' if ds > 0 else 'down'):>12s}")
    sign = [np.sign(r[3]) for r in rows]
    thr = [i for i in range(1, len(rows)) if sign[i] > 0 and sign[i - 1] <= 0]
    if thr:
        i0 = thr[0]
        top = [i for i in range(i0 + 1, len(rows)) if sign[i] < 0]
        print(f"\n  lower stable state: SDAM = 0 (drug-free repair wins)")
        print(f"  unstable threshold  : SDAM between "
              f"{rows[i0 - 1][0]:.2f} and {rows[i0][0]:.2f}")
        if top:
            print(f"  upper stable state  : SDAM between "
                  f"{rows[top[0] - 1][0]:.2f} and {rows[top[0]][0]:.2f} "
                  f"-> self-sustaining structural damage (irreversible TD)")
    else:
        print("\n  (no interior root found: system is monostable)")
    return rows


def main():
    quick = '--quick' in sys.argv
    scenarios()
    if quick:
        return
    exposure_threshold()
    reversibility_window()
    withdrawal_crossover()
    vmat2_dose_response()
    cyp2d6_panel()
    suppression_vs_modification()
    combination_interaction()
    opposed_levers()
    risk_scan()
    latch_bistability()
    print('\n' + '=' * 108)
    print('done — all numbers above are model outputs, not literature values.')
    print('=' * 108)


if __name__ == '__main__':
    main()
