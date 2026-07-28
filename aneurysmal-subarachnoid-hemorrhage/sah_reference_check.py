#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sah_reference_check.py
======================
Independent vectorised numpy/RK4 transcription of the aSAH-DCI QSP model
(sah_mrgsolve_model.R).  Its only job is to COMPUTE every number quoted in
README.md / sah_references.md rather than let them be asserted.

The model in one sentence
------------------------
Delayed cerebral ischaemia after aneurysmal subarachnoid haemorrhage is not
posed here as "vasospasm -> ischaemia", but as ONE SHARED ARTERIOLAR
VASODILATORY RESERVE with FOUR COMPETING CONSUMERS:

    reserve   = RA0 - RAMIN                        (arteriolar resistance range)
    consumers = dRL   large-artery spasm           (the only one angiography sees)
                dRMt  microvascular / pericyte tone
                dRMc  capillary microthrombosis (capillary drop-out)
                dCPP  perfusion-pressure loss (ICP up, MAP down)
                dSD   spreading-depolarisation inverse coupling

    CBF = CPP / (RL + RA + RM),  RA = AREG*clamp(CPP/CBFtarget - RL - RM) +
                                      (1-AREG)*RA0

Ischaemia begins only when the SUM of the demands exhausts the reserve.  So the
map from any single consumer to the endpoint is nonlinear and context-dependent,
and blocking one consumer helps only the subpopulation in which that consumer
was the MARGINAL one.  Nothing in the model contains a "day 4-10" switch: the
DCI time window is the convolution of clot lysis, CSF oxyhaemoglobin scavenging
(haptoglobin genotype dependent) and inducible HO-1.

38 ODEs.  Time unit = DAYS.  RK4, dt = 0.01 d (14.4 min).

Run:  python3 sah_reference_check.py            (writes stdout; ~1-2 min)
"""

import sys
import numpy as np

rng_seed = 20260728
N_SUB = 900              # virtual subjects per arm
T_END = 21.0             # days simulated
DT = 0.010               # integration step (days)
REC_EVERY = 10           # record every 0.1 d

# =============================================================================
# 0. STATE INDEX
# =============================================================================
S = dict(
    # --- drug PK / exposure (8) ---
    NIMG=0, NIMC=1, CLAZ=2, CILO=3, STAT=4, MILR=5, NICA=6, KETA=7,
    # --- clot / haemoglobin arm (7) ---
    CLOT=8, IVH=9, OXYHB=10, HP=11, HEME=12, HO1=13, BOX=14,
    # --- effectors (6) ---
    ET1=15, NOB=16, ROS=17, RHOK=18, INFL=19, PAI=20,
    # --- vascular mechanics (5) ---
    SPASM=21, STRUCT=22, RMIC=23, MTHR=24, AREG=25,
    # --- spreading depolarisation / neuronal (3) ---
    SDSUS=26, SDBUR=27, SDCUM=28,
    # --- perfusion / injury (7) ---
    EDEMA=29, HYDRO=30, ICP=31, MAP=32, OGD=33, INFVOL=34, EBI=35,
    # --- systemic (2) ---
    NAS=36, FLUID=37,
    # --- erythrocyte pool undergoing haemolysis (1) ---
    RBCL=38,
)
NS = len(S)

# =============================================================================
# 1. FIXED SYSTEM PARAMETERS
# =============================================================================
P = dict(
    # --- haemorrhage clearance -------------------------------------------------
    KMAT=0.250,      # clot -> erythrocytes entering haemolysis (1/d)
    KMATV=0.185,     # same, intraventricular compartment (1/d)
    KHEM=0.230,      # haemolysis, RBC pool -> free oxyHb (1/d)
    KDRAIN=0.036,    # extra physical clearance per unit of drainage intensity (1/d)
    KREL=9.20,       # oxyHb released per unit haemolysed (RU)
    KHPCL=1.40,      # haptoglobin-mediated oxyHb scavenging (1/d)
    KMHP=0.45,       # Km, haptoglobin saturation
    KMAC=0.55,       # CD163/macrophage uptake baseline (1/d)
    KCSF=0.22,       # bulk CSF clearance of free oxyHb (1/d)
    KHPIN=0.60, KHPOUT=0.60, KHPUSE=0.85,
    KHEME=0.85,      # oxyHb -> free heme (1/d)
    KHOX=0.95,       # HO-1 dependent heme degradation (1/d)
    KHO1IN=1.30, KHO1OUT=0.36,    # HO-1 induction, t1/2 1.9 d (the lag)
    KBOXIN=0.55, KBOXOUT=0.42,    # bilirubin oxidation products
    # --- effectors -------------------------------------------------------------
    KETIN=0.80, KETOUT=0.80,      # endothelin-1
    TAUNO=0.30,                   # NO bioavailability time constant (d)
    KNOHB=0.85, KNOROS=0.35,      # oxyHb and ROS scavenging weights
    KROSIN=0.60, KROSOUT=1.05,
    TAURHO=1.00,                  # Rho-kinase Ca-sensitisation (slow)
    KINFIN=0.40, KINFOUT=0.55,
    KPAIIN=0.45, KPAIOUT=0.60,
    # --- vascular --------------------------------------------------------------
    KSPIN=0.550, KSPOUT=0.62, SPMAX=0.82,
    KSTIN=0.045, KSTOUT=0.045,    # structural remodelling (slow, sticky)
    TAURM=0.50, KRMIC=1.15,       # microvascular tone time constant / gain
    KMTIN=0.0480, KMTOUT=0.130,   # microthrombosis on/off
    KAREGR=0.16, KAREGL=0.35,     # autoregulation recovery / loss
    # --- spreading depolarisation ---------------------------------------------
    TAUSDS=0.40, SDMAX=7.0, KSDB=0.60,
    KINV=0.55,                    # SD inverse-coupling gain on RA
    # --- perfusion / injury ---------------------------------------------------
    KEDIN=0.30, KEDOUT=0.35,
    KHYIN=0.22, KHYOUT=0.22, KEVD=2.4,
    TAUICP=0.30, TAUMAP=0.20,
    KREP=1.55,                    # oxygen-debt repair (1/d)
    KINFV=0.42, OGD50=0.300, OGDN=4.0, VRISK=190.0,
    KEBIREC=0.030,
    KNAD=6.50, KNAR=0.42,
    KFLIN=0.56, KFLOUT=0.50,
    # --- haemodynamic constants ----------------------------------------------
    CBF0=50.0, CPP0=80.0,
    RL0=0.250, RA0=0.850, RM0=0.500,     # sum = 1.60 = CVR0 (CPP0/CBF0)
    FSEG=0.330,      # fraction of large-artery resistance in the spastic segment
    RAMIN=0.200, RAMAX=2.20,
    OEFMAX=0.85, OEF0=0.35,
    KCTH=1.75,       # penalty on maximum O2 extraction per unit heterogeneity
    # --- drug PK ---------------------------------------------------------------
    KANIM=36.0, KELNIM=9.79, VNIM=100.0, FNIM=0.13,   # nimodipine (oral)
    KELCLZ=8.32, VCLZ=20.0,                            # clazosentan
    KELCIL=1.512, VCIL=90.0, FCIL=0.90,                # cilostazol
    KELMIL=6.65, VMIL=30.0,                            # milrinone
    KELNIC=1.00, VNIC=0.50,                            # intrathecal nicardipine
    KELKET=1.00, VKET=1.00,                            # ketamine (normalised)
    TAUSTAT=1.00,                                      # statin effect site
    # --- drug PD ---------------------------------------------------------------
    EC50NIM=30.0,  EMXNIM_SP=0.16, EMXNIM_RM=0.45, EMXNIM_SD=0.55, EMXNIM_NP=0.50,
    DMAP_NIM_PO=5.0, DMAP_NIM_IV=15.0,
    EC50CLZ=400.0, EMXCLZ_ET=0.70, EMXCLZ_RM=0.08, DMAP_CLZ=6.0, DHGB_CLZ=0.80,
    EC50CIL=600.0, EMXCIL_MT=0.70, EMXCIL_RM=0.28, EMXCIL_SP=0.10,
    EC50MIL=150.0, EMXMIL_RM=0.34, EMXMIL_SP=0.14, DMAP_MIL=-4.0,
    EC50NIC=2000.0, EMXNIC_SP=0.55, EMXNIC_RM=0.22,
    EC50KET=1000.0, EMXKET_SD=0.65,
    EMXSTAT_NO=0.03, EMXSTAT_INF=0.04,
    # --- outcome model ---------------------------------------------------------
    B0=-4.10, B_EBI=3.15, B_INF=1.08, B_HYD=1.05, B_AGE=0.036,
    B_HARM=1.55, B_NA=0.45, B_REB=2.50, B_SD=0.85,
)

# =============================================================================
# 2. VIRTUAL POPULATION
# =============================================================================
def make_population(n=N_SUB, seed=rng_seed):
    r = np.random.default_rng(seed)
    p = {}
    p['AGE'] = np.clip(r.normal(56, 12, n), 25, 88)
    p['WFNS'] = r.choice([1, 2, 3, 4, 5], size=n, p=[.34, .16, .10, .24, .16])
    p['MFI'] = r.choice([1, 2, 3, 4], size=n, p=[.10, .20, .34, .36])
    p['HP22'] = r.binomial(1, 0.36, n).astype(float)
    p['MAP0'] = np.clip(r.normal(96, 11, n), 70, 132)
    p['HGB0'] = np.clip(r.normal(12.6, 1.7, n), 7.5, 17.0)
    p['THRP'] = r.lognormal(0.0, 0.90, n)          # microthrombosis propensity
    p['SPSN'] = r.lognormal(0.0, 0.45, n)          # large-artery vasoreactivity
    p['RMSN'] = r.lognormal(0.0, 0.30, n)          # microvascular reactivity
    p['AREG0'] = np.clip(r.normal(0.86, 0.11, n) - 0.055 * (p['WFNS'] - 1), 0.10, 0.99)
    p['EBI0'] = np.clip(0.08 + 0.155 * (p['WFNS'] - 1) + r.normal(0, 0.07, n), 0.0, 1.0)
    p['RA0i'] = P['RA0'] * r.lognormal(0.0, 0.12, n)
    p['RAMINi'] = P['RAMIN'] * r.lognormal(0.0, 0.15, n)
    p['CLF'] = r.lognormal(0.0, 0.30, n)           # shared clearance factor
    p['VRISKi'] = 42.0 * r.lognormal(0.0, 0.65, n)  # tissue volume at risk (mL)
    p['NAP'] = r.lognormal(0.0, 0.42, n)            # natriuresis / SIADH propensity
    p['FLP'] = r.lognormal(-0.5 * 0.65 ** 2, 0.65, n)  # fluid-retention propensity (mean 1)
    # clot burden by modified Fisher grade: 1 thin, 2 thin+IVH, 3 thick, 4 thick+IVH
    clot_map = np.array([0.35, 0.45, 1.00, 1.06])
    ivh_map = np.array([0.00, 0.45, 0.00, 0.60])
    p['CLOT0'] = clot_map[p['MFI'] - 1] * r.lognormal(0.0, 0.16, n)
    p['IVH0'] = ivh_map[p['MFI'] - 1] * r.lognormal(0.0, 0.20, n)
    p['EVD'] = (p['IVH0'] > 0.30).astype(float)
    p['REB'] = r.binomial(1, 0.045, n).astype(float)   # rebleed before securing
    p['n'] = n
    return p


def cao2_rel(hgb):
    """Arterial O2 content relative to the population reference (Hb 12.6 g/dL)."""
    ref = 12.6 * 1.34 * 0.97 + 0.003 * 95.0
    return (hgb * 1.34 * 0.97 + 0.003 * 95.0) / ref


# =============================================================================
# 3. SCENARIOS
# =============================================================================
def scen(name, **kw):
    d = dict(name=name,
             nim=0.0, nim_iv=0.0,      # oral mg q4h  /  IV mg/h
             clz=0.0,                   # mg/h
             cil=0.0,                   # mg bid
             stat=0.0,                  # 0/1 (simvastatin 40 mg)
             mil=0.0,                   # ug/kg/min
             nic=0.0,                   # mg/d intrathecal
             ket=0.0,                   # 0/1
             drain=1.0,                 # clot-clearance multiplier
             drain_hb=1.0,              # CSF free-Hb clearance multiplier
             maph=0.0,                  # induced hypertension, mmHg added d4-14
             hgb_target=0.0,            # transfusion threshold (g/dL)
             ko_large=0.0, ko_micro=0.0, ko_thromb=0.0, ko_sd=0.0, ko_areg=0.0,
             t_start=1.0, t_stop=14.0)
    d.update(kw)
    return d


SOC = dict(nim=60.0, t_start=0.5, t_stop=21.0)      # oral nimodipine = standard care

SCENARIOS = [
    scen('S1  supportive only (no nimodipine)'),
    scen('S2  SoC: oral nimodipine 60 mg q4h x21d', **SOC),
    scen('S3  IV nimodipine 2 mg/h x14d', nim_iv=2.0, t_start=0.5, t_stop=14.0),
    scen('S4  SoC + clazosentan 15 mg/h d1-14', clz=15.0, **SOC),
    scen('S5  SoC + clazosentan 5 mg/h d1-14', clz=5.0, **SOC),
    scen('S6  SoC + cilostazol 100 mg bid', cil=100.0, **SOC),
    scen('S7  SoC + simvastatin 40 mg qd', stat=1.0, **SOC),
    scen('S8  SoC + early lumbar drainage', drain=2.20, drain_hb=1.85, **SOC),
    scen('S9  SoC + intrathecal nicardipine implant', nic=4.0, **SOC),
    scen('S10 SoC + induced hypertension (MAP +20)', maph=20.0, **SOC),
    scen('S11 SoC + milrinone 0.5 ug/kg/min d4-14', mil=0.5, **SOC),
    scen('S12 SoC + ketamine (SD suppression)', ket=1.0, **SOC),
    scen('S13 SoC + transfusion to Hb 10 g/dL', hgb_target=10.0, **SOC),
    scen('S14 SoC + clazosentan + cilostazol + drainage',
         clz=15.0, cil=100.0, drain=2.20, drain_hb=1.85, **SOC),
]

KNOCKOUTS = [
    scen('K0  SoC (reference)', **SOC),
    scen('K1  perfect large-artery block', ko_large=1.0, **SOC),
    scen('K2  perfect microvascular-tone block', ko_micro=1.0, **SOC),
    scen('K3  perfect microthrombosis block', ko_thromb=1.0, **SOC),
    scen('K4  perfect SD block', ko_sd=1.0, **SOC),
    scen('K5  autoregulation preserved', ko_areg=1.0, **SOC),
    scen('K12 large + micro', ko_large=1.0, ko_micro=1.0, **SOC),
    scen('K13 large + thrombosis', ko_large=1.0, ko_thromb=1.0, **SOC),
    scen('K23 micro + thrombosis', ko_micro=1.0, ko_thromb=1.0, **SOC),
    scen('K1234 all four consumers blocked',
         ko_large=1.0, ko_micro=1.0, ko_thromb=1.0, ko_sd=1.0, **SOC),
]


# =============================================================================
# 4. RIGHT-HAND SIDE
# =============================================================================
def hill(x, x50, n):
    xn = np.power(np.maximum(x, 0.0), n)
    return xn / (xn + x50 ** n)


def perfusion(y, pop, sc, extra):
    """Algebraic haemodynamics.  Returns dict of derived quantities."""
    SPT = np.clip(y[S['SPASM']] + y[S['STRUCT']], 0.0, 0.85)
    if sc['ko_large']:
        SPT = np.zeros_like(SPT)
    RMICv = np.maximum(y[S['RMIC']], 1.0)
    if sc['ko_micro']:
        RMICv = np.ones_like(RMICv)
    MT = np.clip(y[S['MTHR']], 0.0, 0.75)
    if sc['ko_thromb']:
        MT = np.zeros_like(MT)

    # Only the spastic SEGMENT of the large artery obeys r^-4; the rest of the
    # conducting-artery resistance is unaffected.  This is why a 33% angiographic
    # narrowing is a modest resistance load and a 50% one is not.
    RL = P['RL0'] * ((1.0 - P['FSEG'])
                     + P['FSEG'] / np.power(np.clip(1.0 - SPT, 0.15, 1.0), 4.0))
    RMtone = P['RM0'] * (RMICv - 1.0)
    RMthr = P['RM0'] * RMICv * (1.0 / np.clip(1.0 - MT, 0.25, 1.0) - 1.0)
    RM = P['RM0'] + RMtone + RMthr

    ICPv = y[S['ICP']]
    MAPv = y[S['MAP']]
    CPP = np.clip(MAPv - ICPv, 20.0, 200.0)

    SDfrac = np.clip(y[S['SDBUR']] / P['SDMAX'], 0.0, 1.0)
    if sc['ko_sd']:
        SDfrac = np.zeros_like(SDfrac)
    demand = 1.0 + 0.30 * SDfrac

    AREGv = np.ones_like(MAPv) if sc['ko_areg'] else np.clip(y[S['AREG']], 0.0, 1.0)

    CBFt = P['CBF0'] * demand
    RAdes = np.clip(CPP / CBFt - RL - RM, pop['RAMINi'], P['RAMAX'])
    RA = AREGv * RAdes + (1.0 - AREGv) * pop['RA0i']
    # SD inverse neurovascular coupling: in compromised tissue SD CONSTRICTS
    RA = RA * (1.0 + P['KINV'] * SDfrac * np.clip(y[S['SDSUS']], 0, 1))
    RA = np.clip(RA, pop['RAMINi'], P['RAMAX'])

    CBF = CPP / (RL + RA + RM)
    CaO2 = extra['cao2']
    DO2 = CBF * CaO2
    CMRO2 = P['CBF0'] * P['OEF0'] * demand
    # Capillary transit-time heterogeneity (Jespersen & Ostergaard): microvascular
    # constriction and capillary plugging DEGRADE the maximum extractable oxygen
    # fraction, so a CBF that would be adequate in a homogeneous bed is not.
    # This is the term that makes the microcirculation a first-class consumer.
    CTH = (RMICv - 1.0) + 2.5 * MT
    oefmax = P['OEFMAX'] / (1.0 + P['KCTH'] * CTH)
    OEF = np.minimum(oefmax, CMRO2 / np.maximum(DO2, 1e-6))
    ISCH = np.clip(1.0 - DO2 * oefmax / np.maximum(CMRO2, 1e-6), 0.0, 1.0)
    PBTO2 = 25.0 * (DO2 * (1.0 - OEF)) / (P['CBF0'] * 1.0 * (1.0 - P['OEF0']))
    TCD = 60.0 * (CBF / P['CBF0']) / np.power(np.clip(1.0 - SPT, 0.20, 1.0), 2.0)

    # reserve bookkeeping (units of arteriolar resistance)
    reserve = np.maximum(pop['RA0i'] - pop['RAMINi'], 1e-6)
    d_large = RL - P['RL0']
    d_micro = RMtone
    d_thr = RMthr
    d_cpp = np.maximum(0.0, (P['CPP0'] - CPP)) / P['CBF0']
    d_sd = RA * P['KINV'] * SDfrac * np.clip(y[S['SDSUS']], 0, 1)
    # the ischaemic margin: maximum achievable flow / flow actually required
    cbf_max = CPP / (RL + pop['RAMINi'] + RM)
    cbf_thr = CMRO2 / np.maximum(CaO2 * oefmax, 1e-6)
    margin = cbf_max / np.maximum(cbf_thr, 1e-6)
    # extra flow demanded purely by the loss of extraction efficiency
    d_cth = CMRO2 / np.maximum(CaO2, 1e-6) * (1.0 / np.maximum(oefmax, 1e-6)
                                              - 1.0 / P['OEFMAX'])
    return dict(SPT=SPT, RL=RL, RA=RA, RM=RM, CPP=CPP, CBF=CBF, ISCH=ISCH,
                PBTO2=PBTO2, TCD=TCD, OEF=OEF, demand=demand, SDfrac=SDfrac,
                reserve=reserve, d_large=d_large, d_micro=d_micro, d_thr=d_thr,
                d_cpp=d_cpp, d_sd=d_sd, d_cth=d_cth, margin=margin,
                cbf_max=cbf_max, cbf_thr=cbf_thr, CTH=CTH, oefmax=oefmax)


def drug_effects(t, y, pop, sc):
    on = (t >= sc['t_start']) & (t <= sc['t_stop'])
    e = {}
    e['NIM'] = y[S['NIMC']] / (y[S['NIMC']] + P['EC50NIM'])
    e['CLZ'] = y[S['CLAZ']] / (y[S['CLAZ']] + P['EC50CLZ'])
    e['CIL'] = y[S['CILO']] / (y[S['CILO']] + P['EC50CIL'])
    e['MIL'] = y[S['MILR']] / (y[S['MILR']] + P['EC50MIL'])
    e['NIC'] = y[S['NICA']] / (y[S['NICA']] + P['EC50NIC'])
    e['KET'] = y[S['KETA']] / (y[S['KETA']] + P['EC50KET'])
    e['STAT'] = np.clip(y[S['STAT']], 0.0, 1.0)
    e['on'] = on
    return e


def rhs(t, y, pop, sc, extra):
    d = np.zeros_like(y)
    e = drug_effects(t, y, pop, sc)
    h = perfusion(y, pop, sc, extra)
    ISCH = h['ISCH']

    # ---------------- drug PK -------------------------------------------------
    cl = pop['CLF']
    d[S['NIMG']] = -P['KANIM'] * y[S['NIMG']]
    d[S['NIMC']] = (P['KANIM'] * P['FNIM'] * y[S['NIMG']] / P['VNIM'] * 1000.0
                    - P['KELNIM'] * cl * y[S['NIMC']])
    if sc['nim_iv'] > 0:
        rate_iv = sc['nim_iv'] * 24.0 * ((t >= sc['t_start']) & (t <= sc['t_stop']))
        d[S['NIMC']] = d[S['NIMC']] + rate_iv * 1000.0 / P['VNIM']
    d[S['CLAZ']] = -P['KELCLZ'] * cl * y[S['CLAZ']]
    if sc['clz'] > 0:
        r_clz = sc['clz'] * 24.0 * 1000.0 / P['VCLZ'] * ((t >= sc['t_start']) & (t <= sc['t_stop']))
        d[S['CLAZ']] = d[S['CLAZ']] + r_clz
    d[S['CILO']] = -P['KELCIL'] * cl * y[S['CILO']]
    if sc['cil'] > 0:
        r_cil = sc['cil'] * 2.0 * P['FCIL'] * 1000.0 / P['VCIL'] * ((t >= sc['t_start']) & (t <= sc['t_stop']))
        d[S['CILO']] = d[S['CILO']] + r_cil
    d[S['MILR']] = -P['KELMIL'] * cl * y[S['MILR']]
    if sc['mil'] > 0:
        onm = (t >= 4.0) & (t <= sc['t_stop'])
        r_mil = sc['mil'] * 70.0 * 60.0 * 24.0 / 1000.0 * 1000.0 / P['VMIL'] * onm
        d[S['MILR']] = d[S['MILR']] + r_mil
    d[S['NICA']] = -P['KELNIC'] * y[S['NICA']]
    if sc['nic'] > 0:
        r_nic = sc['nic'] * 1000.0 / P['VNIC'] * ((t >= sc['t_start']) & (t <= sc['t_stop']))
        d[S['NICA']] = d[S['NICA']] + r_nic
    d[S['KETA']] = -P['KELKET'] * y[S['KETA']]
    if sc['ket'] > 0:
        d[S['KETA']] = d[S['KETA']] + 1480.0 * P['KELKET'] * ((t >= sc['t_start']) & (t <= sc['t_stop']))
    d[S['STAT']] = (sc['stat'] * ((t >= sc['t_start']) & (t <= sc['t_stop'])) - y[S['STAT']]) / P['TAUSTAT']

    # ---------------- clot / haemoglobin -------------------------------------
    # Clot mass does not release oxyHb directly: subarachnoid erythrocytes must
    # first enter haemolysis.  This CLOT -> RBCL -> oxyHb transit is the ONLY
    # reason the model has a day-5-to-8 window; there is no explicit switch.
    # CSF drainage / cisternal irrigation does NOT speed haemolysis: it physically
    # REMOVES erythrocytes before they lyse, i.e. it truncates the input function.
    krem = P['KDRAIN'] * np.maximum(sc['drain'] - 1.0, 0.0)
    mat = P['KMAT'] * y[S['CLOT']]
    matv = P['KMATV'] * y[S['IVH']]
    hemo = P['KHEM'] * y[S['RBCL']]
    d[S['CLOT']] = -mat - krem * y[S['CLOT']]
    d[S['IVH']] = -matv - krem * y[S['IVH']]
    d[S['RBCL']] = mat + 0.70 * matv - hemo - krem * y[S['RBCL']]
    effhp = 1.0 - 0.40 * pop['HP22']                     # Hp2-2 scavenges worse
    hpav = y[S['HP']] / (y[S['HP']] + P['KMHP'])
    clr_hb = (P['KHPCL'] * effhp * hpav
              + P['KMAC'] * (1.0 + 3.0 * y[S['HO1']])
              + P['KCSF'] * sc['drain_hb'])
    d[S['OXYHB']] = P['KREL'] * hemo - clr_hb * y[S['OXYHB']]
    d[S['HP']] = (P['KHPIN'] * (1.0 + 0.8 * y[S['INFL']]) - P['KHPOUT'] * y[S['HP']]
                  - P['KHPUSE'] * hpav * y[S['OXYHB']])
    d[S['HEME']] = P['KHEME'] * y[S['OXYHB']] - P['KHOX'] * (0.20 + y[S['HO1']]) * y[S['HEME']]
    drv_ho1 = y[S['HEME']] + 0.5 * y[S['OXYHB']]
    d[S['HO1']] = P['KHO1IN'] * drv_ho1 / (1.0 + drv_ho1) - P['KHO1OUT'] * y[S['HO1']]
    d[S['BOX']] = (P['KBOXIN'] * P['KHOX'] * (0.20 + y[S['HO1']]) * y[S['HEME']]
                   - P['KBOXOUT'] * y[S['BOX']])

    # ---------------- effectors ----------------------------------------------
    d[S['ET1']] = (P['KETIN'] * (0.10 + y[S['OXYHB']] + 0.30 * y[S['ROS']]
                                 + 0.25 * y[S['INFL']]) - P['KETOUT'] * y[S['ET1']])
    ETeff = y[S['ET1']] * (1.0 - P['EMXCLZ_ET'] * e['CLZ'])
    nob_t = (1.0 + P['EMXSTAT_NO'] * e['STAT']) / (
        1.0 + P['KNOHB'] * y[S['OXYHB']] + P['KNOROS'] * y[S['ROS']])
    d[S['NOB']] = (nob_t - y[S['NOB']]) / P['TAUNO']
    d[S['ROS']] = (P['KROSIN'] * (0.10 + 0.90 * (y[S['HEME']] + 0.40 * y[S['BOX']]
                                                  + 0.50 * y[S['INFL']]))
                   - P['KROSOUT'] * y[S['ROS']])
    rhok_t = 0.55 * ETeff + 0.30 * y[S['ROS']] + 0.35 * (1.0 - np.clip(y[S['NOB']], 0, 1))
    d[S['RHOK']] = (rhok_t - y[S['RHOK']]) / P['TAURHO']
    d[S['INFL']] = (P['KINFIN'] * (y[S['OXYHB']] + 0.25 * y[S['IVH']] + 0.30 * y[S['EBI']])
                    - P['KINFOUT'] * y[S['INFL']] * (1.0 + P['EMXSTAT_INF'] * e['STAT']))
    d[S['PAI']] = (P['KPAIIN'] * (0.10 + y[S['INFL']] + 0.50 * y[S['OXYHB']])
                   - P['KPAIOUT'] * y[S['PAI']])

    # ---------------- vascular -----------------------------------------------
    vasod_sp = np.clip(P['EMXNIM_SP'] * e['NIM'] + P['EMXNIC_SP'] * e['NIC']
                       + P['EMXMIL_SP'] * e['MIL'] + P['EMXCIL_SP'] * e['CIL'], 0.0, 0.85)
    drive_sp = pop['SPSN'] * (ETeff + 0.55 * y[S['RHOK']]
                              + 0.35 * (1.0 - np.clip(y[S['NOB']], 0, 1))
                              + 0.25 * y[S['BOX']]) * (1.0 - vasod_sp)
    d[S['SPASM']] = (P['KSPIN'] * drive_sp * (P['SPMAX'] - y[S['SPASM']])
                     - P['KSPOUT'] * y[S['SPASM']])
    d[S['STRUCT']] = P['KSTIN'] * y[S['SPASM']] * y[S['RHOK']] - P['KSTOUT'] * y[S['STRUCT']]
    vasod_rm = np.clip(P['EMXNIM_RM'] * e['NIM'] + P['EMXCLZ_RM'] * e['CLZ']
                       + P['EMXCIL_RM'] * e['CIL'] + P['EMXMIL_RM'] * e['MIL']
                       + P['EMXNIC_RM'] * e['NIC'], 0.0, 0.80)
    rmic_t = 1.0 + pop['RMSN'] * P['KRMIC'] * (
        0.24 * ETeff + 0.14 * (1.0 - np.clip(y[S['NOB']], 0, 1))
        + 0.10 * y[S['ROS']] + 0.08 * y[S['INFL']]) * (1.0 - vasod_rm)
    d[S['RMIC']] = (rmic_t - y[S['RMIC']]) / P['TAURM']
    anti = np.clip(P['EMXCIL_MT'] * e['CIL'], 0.0, 0.85)
    d[S['MTHR']] = (P['KMTIN'] * y[S['PAI']] * (1.0 + 0.50 * y[S['INFL']]) * pop['THRP']
                    * (1.0 - y[S['MTHR']]) * (1.0 - anti) - P['KMTOUT'] * y[S['MTHR']])
    d[S['AREG']] = (P['KAREGR'] * (pop['AREG0'] - y[S['AREG']])
                    - P['KAREGL'] * y[S['AREG']] * (0.90 * ISCH + 0.10 * y[S['ROS']]
                                                    + 0.07 * y[S['INFL']] + 0.25 * y[S['EBI']]))

    # ---------------- spreading depolarisation -------------------------------
    sds_t = np.clip(0.10 + 1.25 * ISCH + 0.25 * (1.0 - np.clip(y[S['NOB']], 0, 1))
                    + 0.35 * y[S['EBI']]
                    + 0.40 * y[S['INFVOL']] / np.maximum(pop['VRISKi'], 1.0), 0.0, 1.0)
    d[S['SDSUS']] = (sds_t - y[S['SDSUS']]) / P['TAUSDS']
    sdrate = (P['SDMAX'] * np.power(np.clip(y[S['SDSUS']], 0, 1), 2.0)
              * (1.0 - P['EMXNIM_SD'] * e['NIM'] - P['EMXKET_SD'] * e['KET']))
    sdrate = np.maximum(sdrate, 0.0)
    d[S['SDBUR']] = sdrate - P['KSDB'] * y[S['SDBUR']] / 1.0
    d[S['SDCUM']] = sdrate

    # ---------------- perfusion / injury -------------------------------------
    d[S['EDEMA']] = (P['KEDIN'] * (ISCH + 0.40 * y[S['INFL']]
                                   + 0.50 * y[S['INFVOL']] / np.maximum(pop['VRISKi'], 1.0))
                     - P['KEDOUT'] * y[S['EDEMA']])
    d[S['HYDRO']] = (P['KHYIN'] * (y[S['IVH']] + 0.30 * y[S['INFL']])
                     - P['KHYOUT'] * y[S['HYDRO']] * (1.0 + P['KEVD'] * pop['EVD']))
    icp_t = 8.0 + 13.0 * y[S['HYDRO']] + 11.0 * y[S['EDEMA']] + 5.0 * y[S['IVH']]
    d[S['ICP']] = (icp_t - y[S['ICP']]) / P['TAUICP']
    onh = (t >= 4.0) & (t <= 14.0)
    map_t = (pop['MAP0'] + sc['maph'] * onh
             - (P['DMAP_NIM_PO'] * e['NIM'] if sc['nim'] > 0 else 0.0)
             - (P['DMAP_NIM_IV'] * e['NIM'] if sc['nim_iv'] > 0 else 0.0)
             - P['DMAP_CLZ'] * e['CLZ']
             - P['DMAP_MIL'] * e['MIL'])
    d[S['MAP']] = (map_t - y[S['MAP']]) / P['TAUMAP']
    d[S['OGD']] = ISCH - P['KREP'] * y[S['OGD']]
    npro = 1.0 - P['EMXNIM_NP'] * e['NIM']
    d[S['INFVOL']] = (P['KINFV'] * np.maximum(pop['VRISKi'] - y[S['INFVOL']], 0.0)
                      * hill(y[S['OGD']], P['OGD50'], P['OGDN']) * npro)
    d[S['EBI']] = -P['KEBIREC'] * y[S['EBI']]
    d[S['NAS']] = (-P['KNAD'] * pop['NAP'] * (0.50 * y[S['IVH']] + 0.40 * y[S['INFL']]
                                              + 0.30 * y[S['EBI']])
                   + P['KNAR'] * (140.0 - y[S['NAS']]))
    d[S['FLUID']] = (P['KFLIN'] * pop['FLP'] * e['CLZ']
                     + 0.25 * pop['FLP'] * (sc['maph'] > 0) * onh
                     - P['KFLOUT'] * y[S['FLUID']])
    return d


# =============================================================================
# 5. INTEGRATOR
# =============================================================================
def simulate(pop, sc):
    n = pop['n']
    y = np.zeros((NS, n))
    y[S['CLOT']] = pop['CLOT0']
    y[S['IVH']] = pop['IVH0']
    y[S['HP']] = 1.0
    y[S['NOB']] = 1.0
    y[S['RMIC']] = 1.0
    y[S['AREG']] = pop['AREG0']
    y[S['ICP']] = 8.0 + 13.0 * 0.1 + 5.0 * pop['IVH0']
    y[S['MAP']] = pop['MAP0']
    y[S['NAS']] = 139.0
    y[S['EBI']] = pop['EBI0']

    hgb = pop['HGB0'].copy()
    if sc['hgb_target'] > 0:
        hgb = np.maximum(hgb, sc['hgb_target'])
    hgb = hgb - P['DHGB_CLZ'] * (sc['clz'] > 0)
    extra = dict(cao2=cao2_rel(hgb))

    nstep = int(round(T_END / DT))
    # oral dose schedules (impulses)
    dose_nim = np.arange(sc['t_start'], sc['t_stop'] + 1e-9, 1.0 / 6.0) if sc['nim'] > 0 else np.array([])
    inim = 0

    nrec = nstep // REC_EVERY + 1
    rec = {k: np.zeros((nrec, n)) for k in
           ('t', 'CBF', 'PBTO2', 'ISCH', 'SPT', 'TCD', 'INFVOL', 'OXYHB', 'ET1',
            'MTHR', 'RMIC', 'AREG', 'SDCUM', 'CPP', 'd_large', 'd_micro', 'd_thr',
            'd_cpp', 'd_sd', 'd_cth', 'reserve', 'RA', 'NAS', 'FLUID', 'ICP', 'HP',
            'HO1', 'margin', 'CTH', 'oefmax', 'cbf_max')}
    tISCH = np.zeros(n)      # cumulative days with ISCH > 0.15 after day 3
    inf_d3 = None
    ri = 0
    t = 0.0
    for i in range(nstep + 1):
        if i % REC_EVERY == 0 and ri < nrec:
            h = perfusion(y, pop, sc, extra)
            rec['t'][ri] = t
            for k in ('CBF', 'PBTO2', 'ISCH', 'SPT', 'TCD', 'CPP', 'RA',
                      'd_large', 'd_micro', 'd_thr', 'd_cpp', 'd_sd', 'd_cth',
                      'reserve', 'margin', 'CTH', 'oefmax', 'cbf_max'):
                rec[k][ri] = h[k]
            for k, idx in (('INFVOL', 'INFVOL'), ('OXYHB', 'OXYHB'), ('ET1', 'ET1'),
                           ('MTHR', 'MTHR'), ('RMIC', 'RMIC'), ('AREG', 'AREG'),
                           ('SDCUM', 'SDCUM'), ('NAS', 'NAS'), ('FLUID', 'FLUID'),
                           ('ICP', 'ICP'), ('HP', 'HP'), ('HO1', 'HO1')):
                rec[k][ri] = y[S[idx]]
            ri += 1
        if i == nstep:
            break
        # oral impulses
        while inim < len(dose_nim) and dose_nim[inim] < t + DT:
            y[S['NIMG']] += sc['nim']             # mg into the gut depot
            inim += 1
        if abs(t - 3.0) < DT / 2:
            inf_d3 = y[S['INFVOL']].copy()
        k1 = rhs(t, y, pop, sc, extra)
        k2 = rhs(t + DT / 2, y + DT / 2 * k1, pop, sc, extra)
        k3 = rhs(t + DT / 2, y + DT / 2 * k2, pop, sc, extra)
        k4 = rhs(t + DT, y + DT * k3, pop, sc, extra)
        y = y + DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        y = np.maximum(y, 0.0)
        t += DT
        if t > 3.0:
            h = perfusion(y, pop, sc, extra)
            tISCH += (h['ISCH'] > 0.15) * DT
    if inf_d3 is None:
        inf_d3 = np.zeros(n)

    out = dict(rec=rec, y=y, tISCH=tISCH, inf_d3=inf_d3, pop=pop, sc=sc)
    out.update(endpoints(out))
    return out


# =============================================================================
# 6. ENDPOINTS
# =============================================================================
def endpoints(o):
    rec, y, pop = o['rec'], o['y'], o['pop']
    tt = rec['t'][:, 0]
    w510 = (tt >= 4.0) & (tt <= 11.0)
    spt_max = rec['SPT'][w510].max(axis=0)
    tcd_max = rec['TCD'][w510].max(axis=0)
    pbt_min = rec['PBTO2'][tt >= 3.0].min(axis=0)
    dci_inf = y[S['INFVOL']] - o['inf_d3']
    dci = ((dci_inf >= 3.0) | (o['tISCH'] >= 0.25)).astype(float)
    ang_any = (spt_max >= 0.25).astype(float)
    ang_ms = (spt_max >= 0.33).astype(float)      # CONSCIOUS-1 "moderate-severe"
    ang_sev = (spt_max >= 0.50).astype(float)
    hyd_max = rec['ICP'][tt >= 2.0].max(axis=0)
    na_min = rec['NAS'].min(axis=0)
    fluid_max = rec['FLUID'].max(axis=0)
    harm = 0.85 * np.maximum(0.0, fluid_max - 0.35)
    ebim = pop['EBI0']
    logit = (P['B0'] + P['B_EBI'] * ebim + P['B_INF'] * np.log1p(y[S['INFVOL']])
             + P['B_HYD'] * np.clip((hyd_max - 15.0) / 15.0, 0, 2)
             + P['B_AGE'] * (pop['AGE'] - 56.0) + P['B_HARM'] * harm
             + P['B_SD'] * np.log1p(y[S['SDCUM']]) / 3.0
             + P['B_NA'] * (na_min < 130.0) + P['B_REB'] * pop['REB'])
    ppoor = 1.0 / (1.0 + np.exp(-logit))
    return dict(spt_max=spt_max, tcd_max=tcd_max, pbt_min=pbt_min, dci=dci,
                dci_inf=dci_inf, ang_any=ang_any, ang_ms=ang_ms, ang_sev=ang_sev,
                infvol=y[S['INFVOL']], ppoor=ppoor, na_min=na_min,
                fluid_max=fluid_max, sdcum=y[S['SDCUM']], areg_min=rec['AREG'].min(axis=0),
                icp_max=hyd_max)


# =============================================================================
# 7. REPORTING HELPERS
# =============================================================================
def pct(x):
    return 100.0 * float(np.mean(x))


def head(title):
    print('\n' + '=' * 78)
    print(title)
    print('=' * 78)


def main():
    np.set_printoptions(suppress=True)
    pop = make_population()
    head('POPULATION (N = %d)' % pop['n'])
    print('  age        mean %.1f  (SD %.1f)' % (pop['AGE'].mean(), pop['AGE'].std()))
    print('  WFNS 1-2   %.1f%%   WFNS 4-5 %.1f%%' %
          (pct(pop['WFNS'] <= 2), pct(pop['WFNS'] >= 4)))
    print('  mFisher 3-4 %.1f%%  IVH %.1f%%  EVD %.1f%%' %
          (pct(pop['MFI'] >= 3), pct(pop['IVH0'] > 0), pct(pop['EVD'])))
    print('  Hp2-2      %.1f%%   Hb mean %.1f g/dL' % (pct(pop['HP22']), pop['HGB0'].mean()))

    res = {}
    for sc in SCENARIOS:
        res[sc['name']] = simulate(pop, sc)
        sys.stderr.write('.')
        sys.stderr.flush()
    sys.stderr.write('\n')

    ref = res['S2  SoC: oral nimodipine 60 mg q4h x21d']

    head('R1  SCENARIO TABLE  (N = %d each, same virtual subjects)' % pop['n'])
    print('%-46s %7s %7s %7s %7s %8s %8s %8s' %
          ('scenario', 'angANY%', 'angMS%', 'angSEV%', 'DCI%', 'infDCI', 'PbtO2', 'poor%'))
    for sc in SCENARIOS:
        o = res[sc['name']]
        dm = o['dci'] > 0.5
        print('%-46s %7.1f %7.1f %7.1f %7.1f %8.1f %8.1f %8.1f' %
              (sc['name'], pct(o['ang_any']), pct(o['ang_ms']), pct(o['ang_sev']),
               pct(o['dci']),
               np.median(o['infvol'][dm]) if dm.sum() else 0.0,
               np.median(o['pbt_min']), pct(o['ppoor'])))
    print('  (angMS = calibre loss >=33%, angSEV >=50%;  infDCI = median infarct'
          ' volume in mL AMONG DCI patients;  PbtO2 = median nadir, mmHg;'
          '  poor = mean P(mRS 4-6))')

    head('R2  RISK RATIOS vs SoC (S2)')
    print('%-46s %8s %8s %8s %8s' % ('scenario', 'RR angMS', 'RR DCI', 'RR poor', 'd infarct'))
    for sc in SCENARIOS:
        o = res[sc['name']]
        rr_ms = np.mean(o['ang_ms']) / np.mean(ref['ang_ms'])
        rr_dci = np.mean(o['dci']) / np.mean(ref['dci'])
        rr_p = np.mean(o['ppoor']) / np.mean(ref['ppoor'])
        print('%-46s %8.2f %8.2f %8.2f %8.1f' %
              (sc['name'], rr_ms, rr_dci, rr_p, np.median(o['infvol']) - np.median(ref['infvol'])))

    head('R3  THE DISSOCIATION: angiographic spasm vs tissue outcome (SoC arm)')
    o = ref
    x, v = o['spt_max'], o['infvol']
    r = np.corrcoef(x, v)[0, 1]
    print('  corr(min-caliber loss, final infarct volume)  r = %.3f   r2 = %.3f' % (r, r * r))
    r2 = np.corrcoef(o['tcd_max'], v)[0, 1]
    print('  corr(peak TCD velocity  , final infarct volume) r = %.3f   r2 = %.3f' % (r2, r2 * r2))
    dci_no_ms = np.mean(o['dci'] * (1 - o['ang_ms'])) / np.mean(o['dci'])
    ms_no_dci = np.mean(o['ang_ms'] * (1 - o['dci'])) / max(np.mean(o['ang_ms']), 1e-9)
    print('  DCI WITHOUT moderate-severe angiographic spasm : %.1f%% of all DCI' % (100 * dci_no_ms))
    print('  moderate-severe spasm WITHOUT DCI              : %.1f%% of all such patients' % (100 * ms_no_dci))
    for thr in (120, 160, 200):
        sens = np.mean((o['tcd_max'] > thr) * o['dci']) / np.mean(o['dci'])
        spec = np.mean((o['tcd_max'] <= thr) * (1 - o['dci'])) / np.mean(1 - o['dci'])
        print('  TCD > %3d cm/s as a DCI test: sensitivity %.1f%%  specificity %.1f%%'
              % (thr, 100 * sens, 100 * spec))
    for thr in (20, 15):
        sens = np.mean((o['pbt_min'] < thr) * o['dci']) / np.mean(o['dci'])
        spec = np.mean((o['pbt_min'] >= thr) * (1 - o['dci'])) / np.mean(1 - o['dci'])
        print('  PbtO2 < %2d mmHg as a DCI test: sensitivity %.1f%%  specificity %.1f%%'
              % (thr, 100 * sens, 100 * spec))

    head('R4  RESERVE DECOMPOSITION  (fraction of arteriolar reserve consumed)')
    rec = ref['rec']
    tt = rec['t'][:, 0]
    print('  --- share of total demand (population mean, per day) ---')
    for day in (2, 4, 6, 8, 10, 14):
        j = int(np.argmin(np.abs(tt - day)))
        comp = np.array([np.mean(rec['d_large'][j]), np.mean(rec['d_micro'][j]),
                         np.mean(rec['d_thr'][j]), np.mean(rec['d_cpp'][j]),
                         np.mean(rec['d_sd'][j])])
        sh = 100 * comp / comp.sum()
        print('  day %2d  large %5.1f%%  microtone %5.1f%%  thrombosis %5.1f%%  '
              'CPP %5.1f%%  SD %5.1f%%   demand/reserve %.2f' %
              (day, sh[0], sh[1], sh[2], sh[3], sh[4],
               np.mean((rec['d_large'][j] + rec['d_micro'][j] + rec['d_thr'][j]
                        + rec['d_cpp'][j] + rec['d_sd'][j]) / rec['reserve'][j])))
    print('  large-artery share of TOTAL cerebrovascular resistance:')
    for day in (0, 6):
        j = int(np.argmin(np.abs(tt - day)))
        rl = P['RL0'] * ((1 - P['FSEG']) + P['FSEG'] / np.power(np.clip(1 - rec['SPT'][j], 0.15, 1), 4))
        tot_r = rec['CPP'][j] / rec['CBF'][j]
        print('     day %d: %.1f%%' % (day, 100 * np.mean(rl / tot_r)))

    head('R5  THE EMERGENT TIME WINDOW (no day-4 switch exists in the model)')
    for lab, mask in (('all', np.ones(pop['n'], bool)),
                      ('mFisher 1-2', pop['MFI'] <= 2),
                      ('mFisher 3-4', pop['MFI'] >= 3),
                      ('Hp1-1/1-2', pop['HP22'] == 0),
                      ('Hp2-2', pop['HP22'] == 1)):
        hb = rec['OXYHB'][:, mask].mean(axis=1)
        et = rec['ET1'][:, mask].mean(axis=1)
        sp = rec['SPT'][:, mask].mean(axis=1)
        isc = rec['ISCH'][:, mask].mean(axis=1)
        print('  %-12s oxyHb peak d%4.1f (%.2f RU) | ET-1 peak d%4.1f | caliber nadir d%4.1f '
              '(%.0f%% loss) | ischaemia peak d%4.1f' %
              (lab, tt[hb.argmax()], hb.max(), tt[et.argmax()], tt[sp.argmax()],
               100 * sp.max(), tt[isc.argmax()]))
    print('  DCI incidence by mFisher grade:  ' + '  '.join(
        'mF%d %.1f%%' % (g, pct(ref['dci'][pop['MFI'] == g])) for g in (1, 2, 3, 4)))
    print('  DCI incidence by Hp genotype  :  Hp1-1/1-2 %.1f%%   Hp2-2 %.1f%%' %
          (pct(ref['dci'][pop['HP22'] == 0]), pct(ref['dci'][pop['HP22'] == 1])))
    print('  DCI incidence by autoregulation tertile (AREG0):  ' + '  '.join(
        '%s %.1f%%' % (lab, pct(ref['dci'][m])) for lab, m in
        zip(('low', 'mid', 'high'),
            (pop['AREG0'] <= np.quantile(pop['AREG0'], 1 / 3),
             (pop['AREG0'] > np.quantile(pop['AREG0'], 1 / 3)) &
             (pop['AREG0'] < np.quantile(pop['AREG0'], 2 / 3)),
             pop['AREG0'] >= np.quantile(pop['AREG0'], 2 / 3)))))

    head('R6  PATH-KNOCKOUT EXPERIMENT (perfect, unattainable single-target drugs)')
    kres = {}
    for sc in KNOCKOUTS:
        kres[sc['name']] = simulate(pop, sc)
        sys.stderr.write('.')
        sys.stderr.flush()
    sys.stderr.write('\n')
    k0 = kres['K0  SoC (reference)']
    print('%-40s %8s %8s %8s %9s' % ('knockout', 'DCI%', 'RR DCI', 'infarct', 'RR poor'))
    for sc in KNOCKOUTS:
        o = kres[sc['name']]
        print('%-40s %8.1f %8.2f %8.1f %9.2f' %
              (sc['name'], pct(o['dci']), np.mean(o['dci']) / np.mean(k0['dci']),
               np.median(o['infvol']), np.mean(o['ppoor']) / np.mean(k0['ppoor'])))
    a = 1 - np.mean(kres['K1  perfect large-artery block']['dci']) / np.mean(k0['dci'])
    b = 1 - np.mean(kres['K3  perfect microthrombosis block']['dci']) / np.mean(k0['dci'])
    ab = 1 - np.mean(kres['K13 large + thrombosis']['dci']) / np.mean(k0['dci'])
    print('  interaction (large + thrombosis): singly %.3f + %.3f = %.3f, jointly %.3f'
          '  -> %+.3f (SUB-additive: each perfect block already removes most of the'
          ' DCI, so the pair hits the ceiling)' % (a, b, a + b, ab, ab - (a + b)))
    # the outcome ceiling: abolishing DCI entirely does NOT abolish poor outcome
    kall = kres['K1234 all four consumers blocked']
    print('  OUTCOME CEILING: with all four consumers perfectly blocked, DCI falls to '
          '%.1f%% but poor outcome only falls to %.1f%% of baseline (RR %.2f).' %
          (pct(kall['dci']), 100 * np.mean(kall['ppoor']) / np.mean(k0['ppoor']),
           np.mean(kall['ppoor']) / np.mean(k0['ppoor'])))
    print('  => %.0f%% of the poor-outcome risk is early brain injury and non-DCI'
          ' channels, which is the hard ceiling on ANY DCI-directed drug.' %
          (100 * np.mean(kall['ppoor']) / np.mean(k0['ppoor'])))
    # what a REAL drug achieves on the consumer it targets, vs a perfect block
    clz = res['S4  SoC + clazosentan 15 mg/h d1-14']
    print('  real vs perfect on the SAME consumer (large artery):')
    print('     perfect block   : DCI RR %.2f   poor RR %.2f' %
          (np.mean(kres['K1  perfect large-artery block']['dci']) / np.mean(k0['dci']),
           np.mean(kres['K1  perfect large-artery block']['ppoor']) / np.mean(k0['ppoor'])))
    print('     clazosentan     : DCI RR %.2f   poor RR %.2f' %
          (np.mean(clz['dci']) / np.mean(ref['dci']),
           np.mean(clz['ppoor']) / np.mean(ref['ppoor'])))

    head('R7  INDUCED HYPERTENSION AND MILRINONE DEPEND ON AUTOREGULATION STATE')
    lo = pop['AREG0'] <= np.quantile(pop['AREG0'], 0.33)
    hi = pop['AREG0'] >= np.quantile(pop['AREG0'], 0.67)
    for nm in ('S10 SoC + induced hypertension (MAP +20)', 'S11 SoC + milrinone 0.5 ug/kg/min d4-14'):
        o = res[nm]
        print('  %-44s DCI RR: impaired-autoreg %.2f | intact-autoreg %.2f' %
              (nm, np.mean(o['dci'][lo]) / np.mean(ref['dci'][lo]),
               np.mean(o['dci'][hi]) / np.mean(ref['dci'][hi])))

    head('R8  CLAZOSENTAN: WHERE THE ANGIOGRAPHIC WIN IS LOST')
    o = res['S4  SoC + clazosentan 15 mg/h d1-14']
    print('  moderate-severe angiographic spasm  %.1f%% -> %.1f%%   (RR %.2f)' %
          (pct(ref['ang_ms']), pct(o['ang_ms']),
           np.mean(o['ang_ms']) / np.mean(ref['ang_ms'])))
    print('  peak caliber loss (median)          %.0f%% -> %.0f%%' %
          (100 * np.median(ref['spt_max']), 100 * np.median(o['spt_max'])))
    print('  DCI                                 %.1f%% -> %.1f%%   (RR %.2f)' %
          (pct(ref['dci']), pct(o['dci']), np.mean(o['dci']) / np.mean(ref['dci'])))
    print('  poor outcome                        %.1f%% -> %.1f%%   (RR %.2f)' %
          (pct(ref['ppoor']), pct(o['ppoor']), np.mean(o['ppoor']) / np.mean(ref['ppoor'])))
    print('  fluid-retention burden (median max) %.2f -> %.2f  (IQR %.2f-%.2f on '
          'clazosentan)' %
          (np.median(ref['fluid_max']), np.median(o['fluid_max']),
           np.quantile(o['fluid_max'], .25), np.quantile(o['fluid_max'], .75)))
    for thr in (0.45, 1.0, 1.5):
        print('     fraction above FLUID %.2f : %.1f%% -> %.1f%%' %
              (thr, pct(ref['fluid_max'] > thr), pct(o['fluid_max'] > thr)))
    # decomposition of the null: benefit-only and harm-only counterfactuals
    sc_b_ = scen('clz benefit only', clz=15.0, **SOC)
    P_save = (P['DMAP_CLZ'], P['DHGB_CLZ'], P['KFLIN'])
    P['DMAP_CLZ'], P['DHGB_CLZ'], P['KFLIN'] = 0.0, 0.0, 0.0
    o_b = simulate(pop, sc_b_)
    P['DMAP_CLZ'], P['DHGB_CLZ'], P['KFLIN'] = P_save
    print('  counterfactual, harm channel switched OFF: poor outcome %.1f%% (RR %.2f)' %
          (pct(o_b['ppoor']), np.mean(o_b['ppoor']) / np.mean(ref['ppoor'])))
    print('  => the angiographic effect is real; the endpoint effect is eaten by '
          'redundancy first and harm second.')

    head('R9  MONITORING LEAD TIMES (SoC arm, DCI patients only)')
    m = ref['dci'] > 0.5
    tt = rec['t'][:, 0]

    def first_cross(arr, thr, above=True, after=1.0):
        idx = np.full(arr.shape[1], np.nan)
        ok = (arr > thr) if above else (arr < thr)
        ok = ok & (tt[:, None] > after)
        any_ = ok.any(axis=0)
        idx[any_] = tt[ok.argmax(axis=0)][any_]
        return idx
    t_tcd = first_cross(rec['TCD'][:, m], 120.0, True)
    t_pbt = first_cross(rec['PBTO2'][:, m], 20.0, False)
    t_inf = first_cross(np.diff(rec['INFVOL'][:, m], axis=0, prepend=rec['INFVOL'][:1, m]), 0.05, True)
    for lab, arr in (('TCD > 120 cm/s', t_tcd), ('PbtO2 < 20 mmHg', t_pbt),
                     ('infarct starts', t_inf)):
        v = arr[~np.isnan(arr)]
        print('  %-18s median day %.1f  (IQR %.1f-%.1f, n=%d)' %
              (lab, np.median(v), np.quantile(v, .25), np.quantile(v, .75), len(v)))

    head('R10 SODIUM AND THE OTHER SYSTEMIC CHANNELS (SoC arm)')
    print('  hyponatraemia < 135 mmol/L  %.1f%%   < 130 mmol/L %.1f%%' %
          (pct(ref['na_min'] < 135), pct(ref['na_min'] < 130)))
    print('  ICP > 20 mmHg               %.1f%%   (EVD in %.1f%%)' %
          (pct(ref['icp_max'] > 20), pct(pop['EVD'])))
    print('  cumulative SD count, median %.0f  (DCI %.0f vs no-DCI %.0f)' %
          (np.median(ref['sdcum']), np.median(ref['sdcum'][m]),
           np.median(ref['sdcum'][~m])))
    print('  autoregulation nadir, median %.2f  (DCI %.2f vs no-DCI %.2f)' %
          (np.median(ref['areg_min']), np.median(ref['areg_min'][m]),
           np.median(ref['areg_min'][~m])))

    head('R11 ANAEMIA / OXYGEN CONTENT AS A CONSUMER OF THE SAME NODE')
    for lab, m2 in (('Hb < 10 g/dL', pop['HGB0'] < 10), ('Hb 10-13', (pop['HGB0'] >= 10) & (pop['HGB0'] < 13)),
                    ('Hb >= 13', pop['HGB0'] >= 13)):
        print('  %-14s n=%3d  DCI %.1f%%  median infarct %.1f mL' %
              (lab, int(m2.sum()), pct(ref['dci'][m2]), np.median(ref['infvol'][m2])))
    ot = res['S13 SoC + transfusion to Hb 10 g/dL']
    m3 = pop['HGB0'] < 10
    print('  transfusing only the anaemic subgroup (n=%d): DCI %.1f%% -> %.1f%% (RR %.2f)' %
          (int(m3.sum()), pct(ref['dci'][m3]), pct(ot['dci'][m3]),
           np.mean(ot['dci'][m3]) / max(np.mean(ref['dci'][m3]), 1e-9)))

    head('R12 NIMODIPINE: ROUTE MATTERS BECAUSE THE HARM CHANNEL IS CPP')
    s1, s2, s3 = res['S1  supportive only (no nimodipine)'], ref, res['S3  IV nimodipine 2 mg/h x14d']
    print('  no nimodipine   : angMS %.1f%%  DCI %.1f%%  poor %.1f%%' %
          (pct(s1['ang_ms']), pct(s1['dci']), pct(s1['ppoor'])))
    print('  oral nimodipine : angMS %.1f%%  DCI %.1f%%  poor %.1f%%  (RR poor %.2f)' %
          (pct(s2['ang_ms']), pct(s2['dci']), pct(s2['ppoor']),
           np.mean(s2['ppoor']) / np.mean(s1['ppoor'])))
    print('  IV nimodipine   : angMS %.1f%%  DCI %.1f%%  poor %.1f%%  (RR poor %.2f)' %
          (pct(s3['ang_ms']), pct(s3['dci']), pct(s3['ppoor']),
           np.mean(s3['ppoor']) / np.mean(s1['ppoor'])))
    print('  mean CPP day 4-11: no-nim %.1f | oral %.1f | IV %.1f mmHg' %
          tuple(np.mean(r['rec']['CPP'][(rec['t'][:, 0] >= 4) & (rec['t'][:, 0] <= 11)])
                for r in (s1, s2, s3)))

    head('R13 CALIBRATION TARGETS vs MODEL (SoC arm unless noted)')
    rows = [
        ('any angiographic vasospasm', '60-70%', '%.1f%%' % pct(ref['ang_any'])),
        ('moderate-severe vasospasm (CONSCIOUS-1 placebo 66%)', '~66%', '%.1f%%' % pct(ref['ang_ms'])),
        ('DCI incidence', '~30%', '%.1f%%' % pct(ref['dci'])),
        ('poor outcome mRS 4-6 at 90 d', '30-35%', '%.1f%%' % pct(ref['ppoor'])),
        ('hyponatraemia < 135', '30-50%', '%.1f%%' % pct(ref['na_min'] < 135)),
        ('clazosentan RR mod-sev spasm (CONSCIOUS-1)', '0.35',
         '%.2f' % (np.mean(res['S4  SoC + clazosentan 15 mg/h d1-14']['ang_ms']) / np.mean(ref['ang_ms']))),
        ('clazosentan RR poor outcome (CONSCIOUS-2)', '~1.05',
         '%.2f' % (np.mean(res['S4  SoC + clazosentan 15 mg/h d1-14']['ppoor']) / np.mean(ref['ppoor']))),
        ('nimodipine RR poor outcome (Cochrane 0.67)', '0.67',
         '%.2f' % (np.mean(ref['ppoor']) / np.mean(res['S1  supportive only (no nimodipine)']['ppoor']))),
        ('cilostazol RR DCI (meta ~0.47)', '0.47',
         '%.2f' % (np.mean(res['S6  SoC + cilostazol 100 mg bid']['dci']) / np.mean(ref['dci']))),
        ('simvastatin RR poor outcome (STASH null)', '~1.0',
         '%.2f' % (np.mean(res['S7  SoC + simvastatin 40 mg qd']['ppoor']) / np.mean(ref['ppoor']))),
        ('lumbar drainage RR poor outcome (EARLYDRAIN 0.76)', '0.76',
         '%.2f' % (np.mean(res['S8  SoC + early lumbar drainage']['ppoor']) / np.mean(ref['ppoor']))),
    ]
    for a_, b_, c_ in rows:
        print('  %-52s target %-8s model %s' % (a_, b_, c_))

    print('\nDone.')


if __name__ == '__main__':
    main()
