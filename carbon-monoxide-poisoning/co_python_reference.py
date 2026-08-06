"""=============================================================================
 Carbon Monoxide Poisoning QSP model -- INDEPENDENT Python/scipy implementation
=============================================================================

This file is the cross-check for co_mrgsolve_model.R.  It is not a port: it was
written first, executed, and only then transcribed into mrgsolve, so that the R
model is a transcription of something that has actually been run.  Running this
file regenerates co_reference_output.txt, which is the source of every number
quoted in README.md.

    python3 co_python_reference.py > co_reference_output.txt

The cross-check found and fixed nine defects during development.  They are listed
in README.md under "What the independent implementation caught"; the four most
consequential are commented in place below:

  * the tissue CO driving force (the Haldane free tension is a lung-exchange
    construct and using it here made hyperbaric oxygen pump CO INTO muscle);
  * adduct clearance on a myelin-turnover rather than a days timescale (without
    which delayed sequelae were structurally unreachable at every dose);
  * the strict unknown-parameter check in P(), which exists because a silent
    CN_rate/CN_dose mismatch had been running the whole fire-smoke scenario with
    no cyanide in it, and which then immediately caught a second such divergence;
  * the non-uniform output grid in run(), because a uniform grid over a 90-day
    horizon samples every ~108 min, misses every acute peak, and had inverted a
    conclusion about the timing of an antidote.

Structure of this file:
  1. parameters, helpers and the 45-ODE right-hand side
  2. closed-form derivations (results 1, 2, 3, 9, measurement model)
  3. the scenario / population / sensitivity driver

Requires numpy and scipy only.
=============================================================================
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sys as _sys

# Independent Python/scipy reference implementation of the CO poisoning QSP model.
# 
# Mirrors co_mrgsolve_model.R state-for-state and equation-for-equation.  Written
# first, so that the mrgsolve file is a transcription of something that has been
# executed, and so the numbers in the README are computed.
# 
# State vector (45):
#   0  Cenv     room CO (ppm)
#   1  ACO      CO bound to haemoglobin (mL CO, whole body)
#   2  Pbr      brain tissue CO tension (mmHg x 1e3, i.e. milli-mmHg)
#   3  Pht      heart tissue CO tension
#   4  Pmu      skeletal muscle CO tension  (the slow myoglobin reservoir)
#   5  Prest    other tissue CO tension
#   6  MbHt     cardiac myoglobin CO-occupied fraction
#   7  MbMu     skeletal myoglobin CO-occupied fraction
#   8  CcOb     brain cytochrome c oxidase CO-inhibited fraction
#   9  CcOh     cardiac cytochrome c oxidase CO-inhibited fraction
#  10  ATPb     brain energy charge (fraction of baseline)
#  11  ATPh     cardiac energy charge
#  12  Lac      arterial lactate (mM)
#  13  CBF      cerebral blood flow (mL/100g/min)
#  14  PaCO2    arterial CO2 tension (mmHg)
#  15  XO       xanthine oxidase converted fraction
#  16  ROS      brain ROS index
#  17  NOx      NO / peroxynitrite index
#  18  Neut     adherent neutrophils (index)
#  19  MPO      myeloperoxidase activity (index)
#  20  LPO      lipid peroxidation product (index)
#  21  MBPad    MBP adduct burden (fraction of myelin modified)
#  22  Tcell    autoreactive lymphocyte clone (index)
#  23  Micro    activated microglia
#  24  Demy     demyelination burden
#  25  HO1      heme oxygenase-1 induction
#  26  HIF      HIF-1 alpha
#  27  Necb     brain necrosis burden (globus pallidus / white matter)
#  28  Edema    cerebral oedema
#  29  ICP      intracranial pressure (mmHg)
#  30  TnI      cardiac troponin I (ng/mL)
#  31  EF       LV ejection fraction (fraction)
#  32  CK       creatine kinase (U/L)
#  33  Cr       serum creatinine (mg/dL)
#  34  Cog      neurocognitive composite (fraction of baseline)
#  35  HazD     cumulative DNS hazard
#  36  CNbl     blood cyanide (uM)
#  37  CNtis    tissue cyanide (uM)
#  38  OHCbl    hydroxocobalamin (uM)
#  39  NACc     N-acetylcysteine central (mg/L)
#  40  NACp     NAC peripheral (mg/L)
#  41  GSH      brain glutathione (fraction of baseline)
#  42  Oxy      oxypurinol (mg/L)
#  43  FetCO    fetal COHb fraction
#  44  ATPgp    watershed (globus pallidus) energy charge


NS = 45
PH2O = 47.0

# ------------------------------------------------------------------ parameters
def P(**kw):
    p = dict(
        # --- subject
        Hb=15.0, Vb=5500.0, WT=70.0,
        # --- CFK
        DL=25.0, VA=4200.0, MHald=245.0, Vco_endo=0.007, shunt=0.72,
        # --- exposure / environment
        Rsrc=0.0, Vroom=30000.0, ACH=0.5, texp=0.0, ppm_fix=-1.0,
        # --- therapy
        FiO2=0.21, ATA=1.0, thbo_start=1e9, thbo_dur=90.0, nhbo=1, hbo_gap=480.0,
        to2_start=1e9, to2_stop=1e9, FiO2_trt=0.85, VAtrt=4200.0, FiCO2=0.0,
        ppm_amb=1.5,
        # --- tissue CO kinetics (min)  brain has no Mb, muscle is the reservoir
        tau_br=9.0, tau_ht=45.0, tau_mu=320.0, tau_rest=120.0,
        # --- tissue CO partition and myoglobin
        kappa=1750.0, Mmb=36.0,
        # --- cytochrome c oxidase
        Kc=0.55, Ko=0.55, kon_cco=0.10, koff_cco=0.0022, CcO_base=0.0095,
        Kc_cn=6.0,
        # --- oxygen transport
        P50=26.8, nHill=2.7, CMRO2=3.30, alphaCO2=0.75,
        CBF0=50.0, tau_cbf=1.5, PaCO2_0=40.0,
        # --- metabolism
        tau_atp=6.0, Km_o2=12.0, PtO2ref=33.6, ATPthr=0.92, f_wshed=0.55,
        kLac=0.55, kLacEl=0.030, Lac0=1.0,
        # --- oxidative cascade
        kXOon=0.028, kXOoff=0.0016, kROS=0.185, kROSel=0.10, Kros_o2=45.0,
        kNO=0.225, kNOel=0.09,
        kadh=0.0075, kdeadh=0.0075, hbo_adh=0.82, tau_hbo_adh=420.0,
        kMPO=0.0045, kMPOel=0.0045,
        kLPO=0.00114, kLPOel=0.0021, wMPO=0.85,
        kAd=0.00050, kSpread=0.0004, kAdClr=4.0e-5,
        # --- adaptive immunity (the DNS engine)
        Tprol=0.00120, Tdeath=0.00057, theta=0.450, nT=6.0, Tmax=1.0, T0=0.004,
        kMicro=0.0010, kMicroEl=0.0012,
        kDemy=6.6e-5, kDemyDir=2.4e-5, kRepair=4.5e-5,
        # --- adaptive / protective
        kHO1=0.0016, kHO1el=0.00035, kHIF=0.006, kHIFel=0.0016,
        # --- necrosis / oedema / ICP
        ATPcrit=0.42, kNec=0.0075, kEdema=0.020, kEdemaEl=0.0018,
        ICP0=10.0, kICP=42.0,
        # --- cardiac
        ATPcrit_h=0.72, kTrop=1.5, kTropEl=0.0060, kEF=0.85, kEFrec=0.0016, EF0=0.62,
        # --- muscle / kidney
        CK0=100.0, kCK=300.0, MbCK=0.15, kCKel=0.0009,
        Cr0=0.9, kCr=0.00025, kCrEl=0.0013,
        # --- outcome weights
        wDemy=0.62, wNec=0.55, tau_cog=4320.0, kHaz=0.02,
        # --- cyanide (fire smoke)
        CN_rate=0.0, CN_tau=25.0, CL_cn=0.55, V_cn=25.0, k_rhod=0.0075,
        k_cnt=0.045, lam_cn=2.2, k_ohcbl=0.085,
        # --- drug PK
        ohcbl_dose=0.0, ohcbl_t=1e9, CL_ohc=0.30, V_ohc=18.0,
        nac_dose=0.0, nac_t=1e9, CL_nac=13.0, V1_nac=12.0, Q_nac=8.0, V2_nac=20.0,
        kGSH=0.0032, kGSHox=0.0075, EC50_nac=42.0,
        allo_dose=0.0, allo_t=1e9, CL_oxy=1.2, V_oxy=25.0, IC50_oxy=6.0,
        # --- fetal
        fet_ratio=1.8, fet_tau=250.0, fet=0.0,
        # --- misc switches
        hypothermia=0.0,
    )
    bad = [k for k in kw if k not in p]
    if bad:
        raise KeyError(f"unknown parameter(s): {bad}")   # silent typos are how the
        # fire-smoke scenario once ran with no cyanide in it at all
    p.update(kw)
    return p

def cap_of(p):   return 1.34*p['Hb']/100.0      # mL O2 per mL blood at full saturation

# ------------------------------------------------------------- helper functions
def hill(x, k, n):
    xn = max(x, 0.0)**n
    return xn/(xn + k**n)

def sat_o2(po2, p50, n):
    po2 = max(po2, 1e-9)
    return po2**n/(po2**n + p50**n)

def po2_of_sat(s, p50, n):
    s = min(max(s, 1e-9), 1-1e-9)
    return p50*(s/(1.0-s))**(1.0/n)

def therapy(t, p):
    """Return (FiO2, ATA, VA, PaCO2_target, in_chamber, since_hbo)."""
    inO2 = (t >= p['to2_start']) and (t < p['to2_stop'])
    FiO2 = p['FiO2_trt'] if inO2 else p['FiO2']
    VA   = p['VAtrt'] if inO2 else p['VA']
    ata, chamber, since = 1.0, 0.0, 1e9
    for i in range(int(p['nhbo'])):
        t0 = p['thbo_start'] + i*p['hbo_gap']
        if t >= t0:
            since = min(since, t - t0)
        if t0 <= t < t0 + p['thbo_dur']:
            ata, chamber, FiO2 = p['ATA'], 1.0, 1.0
    return FiO2, ata, VA, chamber, since

def PcO2_of(FiO2, ata, PaCO2, shunt):
    """Mean pulmonary capillary O2 tension.  `shunt` folds venous admixture and
    mask leak into one efficiency factor, calibrated ONCE on Weaver's 74 min."""
    PA = FiO2*(760.0*ata - PH2O) - PaCO2/0.8
    if FiO2 <= 0.25:
        return 100.0                       # room air, by construction
    return shunt*PA

# --------------------------------------------------------------------- the RHS
def rhs(t, y, p):
    (Cenv, ACO, Pbr, Pht, Pmu, Prest, MbHt, MbMu, CcOb, CcOh, ATPb, ATPh, Lac,
     CBF, PaCO2, XO, ROS, NOx, Neut, MPO, LPO, MBPad, Tcell, Micro, Demy, HO1,
     HIF, Necb, Edema, ICP, TnI, EF, CK, Cr, Cog, HazD, CNbl, CNtis, OHCbl,
     NACc, NACp, GSH, Oxy, FetCO, ATPgp) = y
    d = np.zeros(NS)

    FiO2, ata, VA, chamber, since_hbo = therapy(t, p)
    PL   = 760.0*ata - PH2O
    cap  = cap_of(p)

    # ---- 0. environment -----------------------------------------------------
    src = p['Rsrc'] if t < p['texp'] else 0.0
    d[0] = src/p['Vroom']*1e6 - p['ACH']/60.0*Cenv
    Cenv_eff = p['ppm_fix'] if p['ppm_fix'] >= 0.0 and t < p['texp'] else Cenv
    PIco = (0.0 if chamber > 0.5 else (Cenv_eff + p['ppm_amb'])*1e-6*PL)

    # ---- CO2 / ventilation --------------------------------------------------
    # inspired CO2 (carbogen) offsets the hypocapnia of hyperventilation
    PaCO2_ss = p['PaCO2_0']*p['VA']/max(VA, 500.0) + p['FiCO2']*(760.0*ata-PH2O)*0.55
    d[14] = (PaCO2_ss - PaCO2)/2.0

    # ---- 1. CFK: haemoglobin CO --------------------------------------------
    COHb  = ACO/p['Vb']                        # mL CO / mL blood
    FCOHb = min(max(COHb/cap, 0.0), 0.98)
    O2Hb  = max(cap - COHb, 1e-9)
    PcO2  = PcO2_of(FiO2, ata, PaCO2, p['shunt'])
    Bres  = 1.0/p['DL'] + PL/max(VA, 500.0)
    Vco   = p['Vco_endo']*(1.0 + 1.9*HO1)      # HO-1 makes CO: a positive feedback
    d[1]  = Vco - (COHb*PcO2/(p['MHald']*O2Hb) - PIco)/Bres

    PCOb = COHb*PcO2/(p['MHald']*O2Hb)         # mmHg, Haldane free tension (lung only)
    Ptgt = p['kappa']*FCOHb                    # milli-mmHg, set by blood CO CONTENT

    # ---- 2-5. tissue CO (perfusion limited, milli-mmHg) --------------------
    d[2] = (Ptgt - Pbr  )/p['tau_br']
    d[3] = (Ptgt - Pht  )/p['tau_ht']
    d[4] = (Ptgt - Pmu  )/p['tau_mu']
    d[5] = (Ptgt - Prest)/p['tau_rest']

    # ---- oxygen transport ---------------------------------------------------
    PaO2  = FiO2*(760.0*ata - PH2O) - PaCO2/0.8 - 10.0
    SaO2  = sat_o2(PaO2, p['P50'], p['nHill'])
    CaO2  = 1.34*p['Hb']*SaO2*(1.0 - FCOHb) + 0.003*max(PaO2, 0.0)
    P50e  = p['P50']*(1.0 - p['alphaCO2']*FCOHb)      # Haldane left shift

    # cerebral autoregulation: flow defends delivery, and responds to CO2
    fco2  = (0.28 + 0.72/(1.0+np.exp(-(PaCO2-32.0)/6.0))) / \
            (0.28 + 0.72/(1.0+np.exp(-(p['PaCO2_0']-32.0)/6.0)))
    DO2n  = CBF*CaO2/(p['CBF0']*1.34*p['Hb']*0.97)
    CBFss = p['CBF0']*fco2*(1.0 + 1.25*max(0.0, 1.0 - DO2n))
    CBFss = min(CBFss, 2.4*p['CBF0'])
    d[13] = (CBFss - CBF)/p['tau_cbf']

    # brain tissue (= end-capillary) PO2 by Fick on the left-shifted curve
    CMRO2_dem = p['CMRO2']*(1.0 - 0.35*p['hypothermia'])
    CcO_tot_b = 1.0 - (1.0-CcOb)*(1.0 - hill(CNtis, p['Kc_cn'], 1.0))
    CvO2 = CaO2 - CMRO2_dem*100.0/max(CBF, 1.0)
    Sv   = min(max((CvO2 - 0.003*40.0)/max(1.34*p['Hb']*(1.0-FCOHb), 1e-6), 1e-6), 1-1e-6)
    PtO2 = max(po2_of_sat(Sv, P50e, p['nHill']), 0.05)

    # ---- 6-7. myoglobin -----------------------------------------------------
    MbHt_ss = 1.0/(1.0 + PtO2/max(p['Mmb']*Pht*1e-3, 1e-9))
    MbMu_ss = 1.0/(1.0 + 40.0 /max(p['Mmb']*Pmu*1e-3, 1e-9))   # muscle PO2 ~ 40
    d[6] = (MbHt_ss - MbHt)/6.0
    d[7] = (MbMu_ss - MbMu)/25.0

    # ---- 8-9. cytochrome c oxidase -----------------------------------------
    # CO binds only the REDUCED a3 haem, so O2 protects: drive ~ PCO/(Kc(1+PtO2/Ko))
    drive_b = (Pbr*1e-3)/(p['Kc']*(1.0 + PtO2/p['Ko']))
    drive_h = (Pht*1e-3)/(p['Kc']*(1.0 + 40.0 /p['Ko']))
    d[8] = p['kon_cco']*drive_b*(1.0 - CcOb) - p['koff_cco']*(1.0 + 2.6*chamber)*CcOb
    d[9] = p['kon_cco']*drive_h*(1.0 - CcOh) - p['koff_cco']*(1.0 + 2.6*chamber)*CcOh

    # ---- 10-12. energetics --------------------------------------------------
    fo2_ref  = p['PtO2ref']/(p['PtO2ref'] + p['Km_o2'])
    fo2      = PtO2/(PtO2 + p['Km_o2'])
    supply_b = (1.0 - CcO_tot_b)/(1.0 - p['CcO_base'])*min(1.0, fo2/fo2_ref)
    supply_b = min(supply_b, 1.0)
    d[10] = (supply_b - ATPb)/p['tau_atp']
    # the watershed is not handicapped at baseline; it is STEEPER, because it sits
    # lower on its own O2 supply curve, so the same fractional fall in PtO2 costs more
    fo2_gp     = (p['f_wshed']*PtO2)/(p['f_wshed']*PtO2 + p['Km_o2'])
    fo2_gp_ref = (p['f_wshed']*p['PtO2ref'])/(p['f_wshed']*p['PtO2ref'] + p['Km_o2'])
    supply_gp  = min(1.0, (1.0 - CcO_tot_b)/(1.0 - p['CcO_base'])*min(1.0, fo2_gp/fo2_gp_ref))
    d[44] = (supply_gp - ATPgp)/p['tau_atp']
    CcO_tot_h = 1.0 - (1.0-CcOh)*(1.0 - hill(CNtis, p['Kc_cn'], 1.0))
    supply_h  = min(1.0, (1.0 - CcO_tot_h)/(1.0 - p['CcO_base'])*(1.0 - 0.55*MbHt))
    d[11] = (supply_h - ATPh)/p['tau_atp']
    d[12] = p['kLac']*(max(0.0, p['ATPthr']-supply_b) + max(0.0, p['ATPthr']-supply_h)) \
            - p['kLacEl']*(Lac - p['Lac0'])

    # ---- 15-17. oxidative burst --------------------------------------------
    # XDH -> XO conversion needs energy failure; ROS then needs O2 to be restored
    Edef  = max(0.0, p['ATPthr'] - ATPb)          # exactly 0 in health
    d[15] = p['kXOon']*Edef*(1.0-XO) - p['kXOoff']*XO
    o2avail = PtO2/(PtO2 + 8.0) + 0.55*chamber
    I_oxy   = 1.0 - Oxy/(Oxy + p['IC50_oxy'])
    d[16] = p['kROS']*XO*o2avail*I_oxy*(2.0-GSH) - p['kROSel']*ROS
    d[17] = p['kNO']*(max(0.0, FCOHb-0.01) + 0.5*max(0.0, CcOb-p['CcO_base'])) \
            - p['kNOel']*NOx

    # ---- 18-20. neutrophil -> myeloperoxidase -> lipid peroxidation ---------
    # HBO blocks beta-2 integrin adhesion directly (Thom): NOT an oxygen-delivery effect
    hbo_adh = p['hbo_adh']*np.exp(-since_hbo/p['tau_hbo_adh']) if since_hbo < 1e8 else 0.0
    d[18] = p['kadh']*(ROS + 0.6*NOx)*(1.0 - min(hbo_adh, 0.95)) - p['kdeadh']*Neut
    d[19] = p['kMPO']*Neut - p['kMPOel']*MPO
    d[20] = p['kLPO']*(ROS + p['wMPO']*MPO)*(2.0-GSH) - p['kLPOel']*LPO

    # ---- 21. the antigen ----------------------------------------------------
    d[21] = (p['kAd']*LPO + p['kSpread']*Demy)*(1.0 - MBPad) - p['kAdClr']*MBPad

    # ---- 22. autoreactive clone: growth only above an adduct threshold ------
    d[22] = p['Tprol']*Tcell*hill(MBPad, p['theta'], p['nT'])*(1.0 - Tcell/p['Tmax']) \
            - p['Tdeath']*(Tcell - p['T0'])
    d[23] = p['kMicro']*(max(0.0, Tcell-p['T0']) + 0.4*LPO) - p['kMicroEl']*Micro

    # ---- 24. demyelination --------------------------------------------------
    d[24] = (p['kDemy']*max(0.0, Tcell-p['T0']) + p['kDemyDir']*Micro)*(1.0-Demy) \
            - p['kRepair']*Demy

    # ---- 25-26. adaptive ----------------------------------------------------
    d[25] = p['kHO1']*(ROS + 4.0*max(0.0, FCOHb-0.01)) - p['kHO1el']*HO1
    d[26] = p['kHIF']*Edef - p['kHIFel']*HIF

    # ---- 27-29. necrosis, oedema, ICP --------------------------------------
    d[27] = p['kNec']*max(0.0, p['ATPcrit']-ATPgp)*(1.0-Necb)
    d[28] = p['kEdema']*(Necb + 0.35*max(0.0, p['ATPcrit']-ATPgp)) - p['kEdemaEl']*Edema
    d[29] = (p['ICP0'] + p['kICP']*Edema - ICP)/8.0

    # ---- 30-31. heart -------------------------------------------------------
    d[30] = p['kTrop']*max(0.0, p['ATPcrit_h']-ATPh) - p['kTropEl']*TnI
    EFss  = p['EF0']*(1.0 - p['kEF']*max(0.0, p['ATPcrit_h']-ATPh))
    d[31] = (EFss - EF)/25.0 + p['kEFrec']*(p['EF0']-EF)

    # ---- 32-33. muscle and kidney ------------------------------------------
    d[32] = p['kCK']*max(0.0, MbMu-p['MbCK']) - p['kCKel']*(CK - p['CK0'])
    d[33] = p['kCr']*max(0.0, CK-5000.0)/1000.0 - p['kCrEl']*(Cr - p['Cr0'])

    # ---- 34-35. outcome -----------------------------------------------------
    Cogss = 1.0 - p['wDemy']*Demy - p['wNec']*Necb
    d[34] = (Cogss - Cog)/p['tau_cog']
    d[35] = p['kHaz']*Demy + 0.5*p['kHaz']*Necb

    # ---- 36-38. cyanide -----------------------------------------------------
    cn_in = p['CN_rate'] if t < p['CN_tau'] else 0.0
    cnb  = max(CNbl, 0.0)
    ohc  = max(OHCbl, 0.0)
    scav = p['k_ohcbl']*ohc*cnb
    d[36] = cn_in - (p['CL_cn']/p['V_cn'])*cnb - p['k_rhod']*cnb \
            - scav - p['k_cnt']*(cnb - CNtis/p['lam_cn'])
    d[37] = p['k_cnt']*(cnb - CNtis/p['lam_cn'])*(p['V_cn']/40.0) - 0.010*CNtis
    d[38] = -(p['CL_ohc']/p['V_ohc'])*ohc - scav

    # ---- 39-42. drugs -------------------------------------------------------
    d[39] = -(p['CL_nac']/p['V1_nac'])*NACc - (p['Q_nac']/p['V1_nac'])*NACc \
            + (p['Q_nac']/p['V2_nac'])*NACp
    d[40] =  (p['Q_nac']/p['V2_nac'])*NACc - (p['Q_nac']/p['V2_nac'])*NACp
    d[41] = p['kGSH']*(1.0 + 1.6*NACc/(NACc + p['EC50_nac']))*(1.0-GSH) \
            - p['kGSHox']*ROS*GSH
    d[42] = -(p['CL_oxy']/p['V_oxy'])*Oxy

    # ---- 43. fetus ----------------------------------------------------------
    d[43] = (p['fet_ratio']*FCOHb - FetCO)/p['fet_tau']
    return d

# ------------------------------------------------------------------- init/solve
def cohb_equilibrium(p):
    """COHb fraction at which endogenous production balances pulmonary loss."""
    cap = cap_of(p)
    B   = 1.0/p['DL'] + (760.0-PH2O)/p['VA']
    PIa = p['ppm_amb']*1e-6*(760.0-PH2O)
    f   = lambda F: p['Vco_endo'] - (100.0*F/(p['MHald']*(1.0-F)) - PIa)/B
    return brentq(f, 1e-9, 0.5)

def y0_of(p):
    y = np.zeros(NS)
    y[1]  = cohb_equilibrium(p)*cap_of(p)*p['Vb']
    y[10] = y[11] = y[44] = 1.0          # ATP
    y[12] = p['Lac0']
    y[13] = p['CBF0']
    y[14] = p['PaCO2_0']
    y[22] = p['T0']
    y[29] = p['ICP0']
    y[31] = p['EF0']
    y[32] = p['CK0']
    y[33] = p['Cr0']
    y[34] = 1.0
    y[41] = 1.0
    return y

def run(p, tend=60*24*45, npt=4000, y0=None, doses=None):
    """doses: list of (time_min, state_index, amount) bolus additions."""
    y = y0_of(p) if y0 is None else np.array(y0, float)
    doses = sorted(doses or [], key=lambda z: z[0])
    # Peak statistics are read off this grid, so it must be DENSE where the fast
    # variables move.  A uniform grid over a 90-day horizon samples every ~108 min
    # and silently misses the acute peaks entirely -- which once made a well-timed
    # antidote look as though it RAISED tissue cyanide (an artefact of the missed
    # peak in the comparator arm, not of the drug).
    tgrid = np.unique(np.concatenate([
        np.linspace(0.0, 720.0, 481),            # first 12 h, every 1.5 min
        np.linspace(720.0, 4320.0, 241),         # to day 3, every 15 min
        np.linspace(4320.0, tend, max(npt, 2)),  # the long tail
        np.array([d[0] for d in doses] or [0.0])]))
    tgrid = tgrid[tgrid <= tend]
    out_t, out_y, t0 = [], [], 0.0
    breaks = sorted(set([0.0] + [d[0] for d in doses] + [tend]))
    for a, b in zip(breaks[:-1], breaks[1:]):
        for dt, idx, amt in doses:
            if abs(dt-a) < 1e-9: y[idx] += amt
        sub = tgrid[(tgrid >= a) & (tgrid <= b)]
        if len(sub) < 2: sub = np.array([a, b])
        s = solve_ivp(rhs, (a, b), y, args=(p,), t_eval=sub, method='LSODA',
                      rtol=1e-7, atol=1e-9, max_step=30.0)
        out_t.append(s.t); out_y.append(s.y); y = s.y[:, -1].copy()
    return np.concatenate(out_t), np.hstack(out_y)

def cohb(y, p):  return y[1]/(cap_of(p)*p['Vb'])


# =============================================================================
# SELF-ALIAS: the two blocks below were developed as separate scripts
# =============================================================================
C = _sys.modules[__name__]


# =============================================================================
# CLOSED-FORM DERIVATIONS (results 1, 2, 3, 9 and the measurement model)
# =============================================================================

def cfk_only(t, A, PcO2, p):
    """The CFK sub-system alone, extracted from the full model's dxdt_ACO."""
    COHb = A[0]/p['Vb']
    O2Hb = max(C.cap_of(p) - COHb, 1e-9)
    B    = 1.0/p['DL'] + p['PL_'] /max(p['VA_'], 500.0)
    return [p['Vco_endo'] - (COHb*PcO2/(p['MHald']*O2Hb))/B]

def halflife(PcO2, ata=1.0, VA=None, F0=0.30, **kw):
    p = C.P(**kw)
    p['PL_'] = 760.0*ata - PH2O
    p['VA_'] = VA if VA else p['VA']
    A0 = F0*C.cap_of(p)*p['Vb']; tgt = 0.5*A0
    ev = lambda t, y, *a: y[0]-tgt; ev.terminal, ev.direction = True, -1
    s = solve_ivp(cfk_only, (0, 60000), [A0], args=(PcO2, p), events=ev,
                  rtol=1e-11, atol=1e-13, max_step=10.0)
    return s.t_events[0][0] if len(s.t_events[0]) else np.inf

p0 = C.P()
B  = 1.0/p0['DL'] + (760.0-PH2O)/p0['VA']
print("="*80)
print("RESULT 1.  What sets the CO elimination half-life, and the hyperbaric floor")
print("="*80)
print(f"  B = 1/DL + PL/VA = {1/p0['DL']:.5f} + {(760.0-PH2O)/p0['VA']:.5f} = {B:.5f} mmHg.min/mL")
print(f"  ventilation is {100*((760.0-PH2O)/p0['VA'])/B:.1f}% of the CO transfer resistance,"
      f" membrane diffusion {100*(1/p0['DL'])/B:.1f}%")
print()
print("   delivery                     PcO2 (mmHg)   model t1/2   observed")
for nm, PcO2, obs in [("room air, FiO2 0.21", 100.0, "320 min"),
                      ("nasal cannula 6 L/min", 250.0, "-"),
                      ("non-rebreather mask", 432.0, "74 min"),
                      ("intubated, FiO2 1.0", 560.0, "-")]:
    print(f"   {nm:28s} {PcO2:8.0f}      {halflife(PcO2):8.1f} min   {obs}")
print("   (Weaver 2000, Chest: room air 320 min, tight 100% O2 mask 74 min)")
print()
print("   -- hyperbaric: pressure raises the driving force AND the resistance --")
for ata in [1.5, 2.0, 2.4, 2.8, 3.0, 6.0, 20.0]:
    PcO2 = 760.0*ata - PH2O - 50.0
    Bh   = 1.0/p0['DL'] + (760.0*ata-PH2O)/p0['VA']
    print(f"   HBO {ata:5.1f} ATA  PcO2 {PcO2:8.0f}  B {Bh:.4f} ({Bh/B:5.2f}x air)"
          f"   t1/2 = {halflife(PcO2, ata=ata):6.1f} min")
O2Hb  = C.cap_of(p0)*(1-0.30)
kfl   = p0['VA']/(p0['MHald']*O2Hb*p0['Vb'])
tfl   = np.log(2)/kfl
print()
print("   ANALYTIC FLOOR.  As ATA -> inf, PcO2 -> 760*ATA and PL/VA -> 760*ATA/VA, so")
print("      PcO2/B -> VA   and   k -> VA/(M [O2Hb] Vb)")
print(f"      = {p0['VA']:.0f}/({p0['MHald']:.0f} x {O2Hb:.4f} x {p0['Vb']:.0f})"
      f" = {kfl:.5f} /min   ->   t1/2_min = {tfl:.1f} min")
print("      Depends ONLY on alveolar ventilation, haemoglobin and blood volume.")
print(f"      Observed ~20 min at 2.5-3 ATA is {100*(1-20/tfl):.0f}% BELOW this floor:")
print("      CFK is FALSIFIED at hyperbaric pressure.  Reported, not fitted away.")
print("      Corollary: inside a chamber the only remaining lever is ventilation.")
for VAL in [4.2, 6, 8, 10, 12]:
    k = VAL*1000.0/(p0['MHald']*O2Hb*p0['Vb'])
    print(f"         VA {VAL:5.1f} L/min -> floor t1/2 {np.log(2)/k:5.1f} min")

print()
print("="*80)
print("RESULT 2.  HBO works before it has removed any CO")
print("="*80)
def diss(ata, FiO2=1.0, PaCO2=40.0, Aa=10.0):
    PAO2 = FiO2*(760.0*ata - PH2O) - PaCO2/0.8
    return 0.003*(PAO2-Aa), PAO2-Aa
avdo2, cmro2_av = 4.6, 6.3
for ata in [1.0, 2.0, 2.4, 2.8, 3.0]:
    d, pao2 = diss(ata)
    print(f"  {ata:.1f} ATA 100% O2: PaO2 {pao2:7.0f} mmHg  dissolved {d:5.2f} mL/dL"
          f"  = {d/avdo2:5.2f}x whole-body A-V diff, {d/cmro2_av:5.2f}x cerebral")
astar = brentq(lambda a: diss(a)[0] - cmro2_av, 1.0, 5.0)
print(f"  -> dissolved O2 alone covers the entire cerebral demand at {astar:.2f} ATA.")
print("     The therapeutic effect PRECEDES the pharmacokinetics.")

print()
print("="*80)
print("RESULT 3.  CO is worse than the anaemia with the identical arterial O2 content")
print("="*80)
P50, nH, alpha = p0['P50'], p0['nHill'], p0['alphaCO2']
inv  = lambda s, p50: p50*(s/(1.0-s))**(1.0/nH)
CaO2 = lambda Hb, F: 1.34*Hb*0.97*(1.0-F) + 0.003*95.0
def pvo2(Hb, F, ext=5.0):
    cap = 1.34*Hb*(1.0-F); Cv = CaO2(Hb, F) - ext
    Sv  = min(max((Cv - 0.12)/cap, 1e-6), 1-1e-6)
    return inv(Sv, P50*(1.0 - alpha*F))
print("  COHb   CaO2   equivalent   P50eff   PvO2 with CO   PvO2 of that   deficit")
print("          mL/dL  anaemia Hb   mmHg       mmHg          anaemia       mmHg")
for F in [0.0, 0.10, 0.20, 0.30, 0.40, 0.50]:
    ca = CaO2(15.0, F)
    Hbeq = brentq(lambda h: CaO2(h, 0.0) - ca, 0.5, 25.0)
    print(f"  {100*F:4.0f}% {ca:6.2f}   {Hbeq:6.2f}    {P50*(1-alpha*F):6.2f}   "
          f"{pvo2(15.0,F):9.1f}     {pvo2(Hbeq,0.0):9.1f}     {pvo2(Hbeq,0.0)-pvo2(15.0,F):6.1f}")
print("  -> at COHb 40% the content equals Hb 9.0 g/dL, a level nobody transfuses")
print("     urgently, yet the tissue PO2 is ~9 mmHg lower.  Two hits from one ligand.")

print()
print("="*80)
print("RESULT 9.  The fetal compartment sets the treatment duration")
print("="*80)
def matfet(t, y, PcO2, p):
    A, Ff = y
    dA = cfk_only(t, [A], PcO2, p)[0]
    Fm = A/(C.cap_of(p)*p['Vb'])
    return [dA, (p['fet_ratio']*Fm - Ff)/p['fet_tau']]
def cross(tt, v, thr=0.05):
    ix = np.where(v < thr)[0]
    return tt[ix[0]] if len(ix) else np.nan
for nm, PcO2 in [("room air", 100.0), ("non-rebreather mask", 432.0),
                 ("intubated FiO2 1.0", 560.0)]:
    p = C.P(); p['PL_'] = 760.0-PH2O; p['VA_'] = p['VA']
    F0 = 0.30
    s = solve_ivp(matfet, (0, 6000), [F0*C.cap_of(p)*p['Vb'], p['fet_ratio']*F0],
                  args=(PcO2, p), rtol=1e-10, atol=1e-12, dense_output=True, max_step=5.0)
    tt = np.linspace(0, 6000, 60000)
    Fm = s.sol(tt)[0]/(C.cap_of(p)*p['Vb']); Ff = s.sol(tt)[1]
    tm, tf = cross(tt, Fm), cross(tt, Ff)
    print(f"  {nm:22s} maternal COHb<5% at {tm:6.1f} min, FETAL at {tf:6.1f} min -> {tf/tm:4.2f}x")
print("  -> the bedside 'treat for ~5x the maternal clearance time' rule is DERIVED.")
print("     Note the direction: the better the maternal treatment, the LARGER the gap.")

print()
print("="*80)
print("MEASUREMENT MODEL.  Pulse oximetry reports the poison as if it were oxygen")
print("="*80)
E6 = {'O2': 0.081, 'HH': 0.845, 'CO': 0.083}
E9 = {'O2': 0.290, 'HH': 0.170, 'CO': 0.010}
def spo2(FCO):
    fO2 = (1-FCO)*0.97; fHH = (1-FCO)*0.03
    R = (fO2*E6['O2']+fHH*E6['HH']+FCO*E6['CO'])/(fO2*E9['O2']+fHH*E9['HH']+FCO*E9['CO'])
    f = lambda s: (s*E6['O2']+(1-s)*E6['HH'])/(s*E9['O2']+(1-s)*E9['HH']) - R
    return brentq(f, 1e-9, 1-1e-9), fO2
print("   true SaO2   COHb    SpO2 displayed   saturation gap")
for FCO in [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60]:
    sp, tru = spo2(FCO)
    print(f"   {100*tru:7.1f}%  {100*FCO:5.1f}%     {100*sp:7.1f}%       {100*(sp-tru):6.1f} points")
print("   -> SpO2 tracks (O2Hb + COHb): the more CO the patient carries, the more")
print("      reassuring the monitor becomes.  The gap IS the COHb.")


# =============================================================================
# SCENARIO, VIRTUAL POPULATION AND SENSITIVITY DRIVER
# =============================================================================
D = 1440.0
PPM = {'mild': (538., 60), 'moderate': (1409., 60), 'severe': (2298., 60),
       'critical': (4198., 45)}
PPM_S, PPM_C = PPM['severe'][0], PPM['critical'][0]
MASK = dict(to2_start=90., to2_stop=450.)

def sim(ppm, mins, tend=90*D, npt=700, doses=None, **kw):
    p = C.P(ppm_fix=ppm, texp=mins, **kw)
    t, y = C.run(p, tend=tend, npt=npt, doses=doses)
    F = np.array([C.cohb(y[:, i], p) for i in range(y.shape[1])])
    return p, t, y, F

def Q(t, y, F):
    return dict(COHb=100*F.max(), CcO=100*y[8].max(), ATP=y[10].min(), GP=y[44].min(),
                Lac=y[12].max(), LPO=y[20].max(), MBP=y[21].max(), Tc=y[22].max(),
                Demy=y[24, -1], Nec=y[27, -1], TnI=y[30].max(), CK=y[32].max(),
                Cog=y[34, -1])

print("#"*80)
print("#  CARBON MONOXIDE POISONING QSP MODEL -- reference output")
print("#  45 ODEs, independent Python/scipy implementation")
print("#  Re-run after three analysis-pipeline fixes: a uniform output grid that")
print("#  undersampled every acute peak, `shunt` promoted from a function default to a")
print("#  real parameter, and a CN_rate/CN_dose name mismatch that had been silently")
print("#  running the fire-smoke scenario with no cyanide in it.")
print("#"*80)

print("\n" + "="*80); print("EXPOSURE SCALE"); print("="*80)
for k, (ppm, mins) in PPM.items():
    p, t, y, F = sim(ppm, mins, tend=float(mins), npt=50)
    print(f"  {k:9s} {ppm:7.0f} ppm x {mins:3d} min  ->  peak COHb {100*F.max():5.1f}%")
print("  OSHA 8-h PEL 50 ppm | NIOSH IDLH 1200 ppm | engine in a closed garage >30000 ppm")

print("\n" + "="*80); print("TREATMENT SCENARIOS"); print("="*80)
S = [
 ("01 severe, no treatment",         PPM_S, 60, {}, None),
 ("02 severe, nasal cannula",        PPM_S, 60, dict(to2_start=90., to2_stop=450., FiO2_trt=0.44), None),
 ("03 severe, O2 mask 6 h",          PPM_S, 60, MASK, None),
 ("04 severe, O2 mask 24 h",         PPM_S, 60, dict(to2_start=90., to2_stop=90+1440.), None),
 ("05 severe, intubated+hypervent",  PPM_S, 60, dict(to2_start=90., to2_stop=450., VAtrt=8000.), None),
 ("06 severe, carbogen 5% CO2",      PPM_S, 60, dict(to2_start=90., to2_stop=450., VAtrt=8000., FiCO2=0.05), None),
 ("07 severe, HBO x1 @2h 3ATA",      PPM_S, 60, dict(**MASK, thbo_start=120., ATA=3.0, nhbo=1), None),
 ("08 severe, HBO x3 @2h (Weaver)",  PPM_S, 60, dict(to2_start=90., to2_stop=90+2880., thbo_start=120., ATA=3.0, nhbo=3), None),
 ("09 severe, HBO x3 @20h",          PPM_S, 60, dict(to2_start=90., to2_stop=90+2880., thbo_start=1200., ATA=3.0, nhbo=3), None),
 ("10 severe, HBO x3 @2h 2.0ATA",    PPM_S, 60, dict(to2_start=90., to2_stop=90+2880., thbo_start=120., ATA=2.0, nhbo=3), None),
 ("11 critical, O2 mask",            PPM_C, 45, dict(to2_start=60., to2_stop=60+1440.), None),
 ("12 critical, HBO x3 early",       PPM_C, 45, dict(to2_start=60., to2_stop=60+2880., thbo_start=90., ATA=3.0, nhbo=3), None),
 ("13 severe + anaemia Hb 9",        PPM_S, 60, dict(**MASK, Hb=9.0), None),
 ("14 severe + NAC 150 mg/kg",       PPM_S, 60, MASK, (90.0, 39, 150.0*70/12.0)),
]
hdr = (f"{'scenario':32s}{'COHb%':>7s}{'CcO%':>6s}{'ATPb':>6s}{'ATPgp':>6s}{'Lac':>6s}"
       f"{'LPO':>6s}{'MBPad':>7s}{'Demy':>7s}{'Nec':>6s}{'TnI':>6s}{'CK':>7s}{'Cog':>6s}{'DNS':>5s}")
print(hdr); print("-"*len(hdr))
for lab, ppm, mins, kw, dose in S:
    p, t, y, F = sim(ppm, mins, doses=[dose] if dose else None, **kw)
    q = Q(t, y, F)
    print(f"{lab:32s}{q['COHb']:7.1f}{q['CcO']:6.1f}{q['ATP']:6.3f}{q['GP']:6.3f}"
          f"{q['Lac']:6.1f}{q['LPO']:6.2f}{q['MBP']:7.3f}{q['Demy']:7.3f}{q['Nec']:6.3f}"
          f"{q['TnI']:6.1f}{q['CK']:7.0f}{q['Cog']:6.3f}"
          f"{'YES' if q['Demy']>0.05 else '.':>5s}")

print("\n" + "="*80)
print("RESULT A.  Two occupancies, two clocks -- COHb is the one that is measured")
print("="*80)
p, t, y, F = sim(PPM_S, 60, npt=1500, **MASK)
def tfrac(sig, fr):
    i0 = int(np.argmax(sig)); ix = np.where(sig[i0:] < fr*sig.max())[0]
    return t[i0+ix[0]] if len(ix) else np.nan
for nm, sig in [("COHb (blood, measured)", F), ("brain tissue CO", y[2]),
                ("cardiac myoglobin-CO", y[6]), ("skeletal muscle CO", y[4]),
                ("brain cytochrome c oxidase", y[8])]:
    print(f"  {nm:28s} peak {sig.max():9.4g}  t(50%) {tfrac(sig,.5):7.1f} min"
          f"  t(10%) {tfrac(sig,.1):8.1f} min")
i0 = int(np.argmax(F)); ix = np.where(F[i0:] < 0.05)[0]; tN = t[i0+ix[0]]
print(f"  -> COHb reaches the 'normal' 5% at {tN:.0f} min ({tN/60:.1f} h); brain cytochrome c")
print(f"     oxidase is still {100*np.interp(tN,t,y[8]):.1f}% inhibited at that moment, and reaches 10% of")
print(f"     its own peak only at {tfrac(y[8],.1):.0f} min -- {tfrac(y[8],.1)/tfrac(F,.1):.1f}x later than COHb does.")

print("\n" + "="*80)
print("RESULT B.  DNS is a bistable switch with a computable threshold")
print("="*80)
p0 = C.P(); Hc = p0['Tdeath']/p0['Tprol']
Mc = p0['theta']*(Hc/(1.0-Hc))**(1.0/p0['nT'])
print(f"  expands iff Tprol*H(MBPad) > Tdeath => H > {Hc:.4f};  MBPad_crit = {Mc:.5f}")
print("    ppm   peak COHb   peak LPO   peak MBPad   /crit   Tcell   Demy   Cog    DNS")
for ppm in [700, 1100, 1500, 1700, 1900, 2300, 3200, 4000]:
    p, t, y, F = sim(ppm, 60, **MASK); q = Q(t, y, F)
    print(f"  {ppm:5.0f}   {q['COHb']:7.1f}%   {q['LPO']:8.2f}   {q['MBP']:9.3f}  "
          f"{q['MBP']/Mc:6.2f}x  {q['Tc']:6.3f} {q['Demy']:6.3f} {q['Cog']:6.3f}   "
          f"{'YES' if q['Demy']>0.05 else 'no':>4s}")



# =============================================================================
# VIRTUAL POPULATION (RESULT C).  This block dominates the runtime:
# 150 subjects x 90-day horizon x 7 arms.
# =============================================================================
# ---------------------------------------------------------------- virtual population
print("\n" + "="*80)
print("RESULT C.  Virtual population: the clinical incidence curve is a threshold DISTRIBUTION")
print("="*80)
rng = np.random.default_rng(20260806)
N = 150
def vpop(ppm, mins, n=N, **kw):
    hits, cogs = 0, []
    for i in range(n):
        kv = dict(
            Hb      = float(np.clip(rng.normal(14.2, 1.7), 8.0, 18.0)),
            VA      = float(np.clip(rng.normal(4200, 500), 3000, 5600)),
            kAd     = float(p0['kAd']*np.exp(rng.normal(0, 0.30))),
            Tprol   = float(p0['Tprol']*np.exp(rng.normal(0, 0.22))),
            kon_cco = float(p0['kon_cco']*np.exp(rng.normal(0, 0.25))),
            CBF0    = float(np.clip(rng.normal(50, 7), 32, 68)),
            kXOon   = float(p0['kXOon']*np.exp(rng.normal(0, 0.30))),
        )
        p, t, y, F = sim(ppm, mins, tend=90*D, npt=700, **kv, **kw)
        hits += int(y[24, -1] > 0.05); cogs.append(y[34, -1])
    return 100.0*hits/n, float(np.mean(cogs))
print("   exposure           n    DNS incidence   mean cognitive score at 90 d")
for lab, ppm, mins in [("mild  (COHb 10%)", PPM['mild'][0], 60),
                       ("moderate (COHb 25%)", PPM['moderate'][0], 60),
                       ("severe  (COHb 40%)", PPM_S, 60),
                       ("critical (COHb 55%)", PPM_C, 45)]:
    inc, cg = vpop(ppm, mins, **MASK)
    print(f"   {lab:20s} {N:4d}   {inc:8.1f}%        {cg:.3f}")
print("  observed: DNS in 10-40% of symptomatic patients (Weaver 2007; Pepe 2011; Choi 1983)")
print("\n   same populations, with early HBO x3:")
for lab, ppm, mins in [("moderate (COHb 25%)", PPM['moderate'][0], 60), ("severe  (COHb 40%)", PPM_S, 60),
                       ("critical (COHb 55%)", PPM_C, 45)]:
    inc, cg = vpop(ppm, mins, to2_start=90., to2_stop=90+2880.,
                   thbo_start=120., ATA=3.0, nhbo=3)
    print(f"   {lab:20s} {N:4d}   {inc:8.1f}%        {cg:.3f}")

print("\n" + "="*80)
print("RESULT D.  The HBO window is set by the adduct, not the carboxyhaemoglobin")
print("="*80)
p, t, y, F = sim(PPM_C, 45, to2_start=60., to2_stop=60+2880.); qb = Q(t, y, F)
print(f"   {'HBO start':>11s} {'COHb then':>10s} {'peak MBPad':>11s} {'Demy':>7s} {'Cog':>7s}")
print(f"   {'none':>11s} {'-':>10s} {qb['MBP']:11.3f} {qb['Demy']:7.3f} {qb['Cog']:7.3f}")
for th in [65, 90, 120, 150, 180, 210, 240, 300, 360, 480, 720, 1440]:
    p, t, y, F = sim(PPM_C, 45, to2_start=60., to2_stop=60+2880.,
                     thbo_start=float(th), ATA=3.0, nhbo=3)
    q = Q(t, y, F)
    print(f"   {th:8.0f} min {100*np.interp(th,t,F):9.1f}% {q['MBP']:11.3f} "
          f"{q['Demy']:7.3f} {q['Cog']:7.3f}")

print("\n" + "="*80)
print("RESULT E.  HBO's actions are REDUNDANT, not additive")
print("="*80)
base = dict(to2_start=60., to2_stop=60+2880., thbo_start=90., ATA=3.0, nhbo=3)
p, t, y, F = sim(PPM_C, 45, to2_start=60., to2_stop=60+2880.); qn = Q(t, y, F)
print(f"   normobaric oxygen only                        Cog {qn['Cog']:.3f}  MBPad {qn['MBP']:.3f}")
for lab, kw in [("full HBO x3 at 3.0 ATA", {}),
                ("3.0 ATA, anti-adhesion action REMOVED", dict(hbo_adh=0.0)),
                ("1.35 ATA, anti-adhesion KEPT", dict(ATA=1.35)),
                ("1.35 ATA AND no anti-adhesion", dict(ATA=1.35, hbo_adh=0.0)),
                ("2.0 ATA", dict(ATA=2.0)),
                ("x1 session instead of x3", dict(nhbo=1)),
                ("x6 sessions instead of x3", dict(nhbo=6))]:
    kk = dict(base); kk.update(kw)
    p, t, y, F = sim(PPM_C, 45, **kk); q = Q(t, y, F)
    print(f"   {lab:45s} Cog {q['Cog']:.3f}  MBPad {q['MBP']:.3f}")
print("  -> either action alone is SUFFICIENT; only removing both loses the benefit.")
print("     Session number is irrelevant in the model: x1 = x3 = x6.")

print("\n" + "="*80)
print("RESULT F.  Fire smoke: CO and cyanide converge on the SAME enzyme")
print("="*80)
print("   arm                                      CcO_CO% CcO_tot% ATPb   Lac   MBPad   Cog")
for lab, ppm, cn, extra, dose in [
   ("CO alone (COHb ~30%)",           1700., 0.0, {}, None),
   ("CN alone (2 uM/min x 25 min)",      0., 2.0, {}, None),
   ("CO + CN (fire smoke)",           1700., 2.0, {}, None),
   ("CO + CN + OHCbl 5 g at 5 min",   1700., 2.0, {}, (5.0, 38, 205.0)),
   ("CO + CN + OHCbl 5 g at 25 min",  1700., 2.0, {}, (25.0, 38, 205.0)),
   ("CO + CN + OHCbl 5 g at 50 min",  1700., 2.0, {}, (50.0, 38, 205.0)),
   ("CO + CN + early HBO x3 (no OHCbl)", 1700., 2.0,
        dict(thbo_start=120., ATA=3.0, nhbo=3), None)]:
    kk = dict(to2_start=90., to2_stop=90+1440.); kk.update(extra)
    p, t, y, F = sim(ppm, 60, CN_rate=cn, doses=[dose] if dose else None, **kk)
    ct = y[37].max(); fcn = ct/(ct + p['Kc_cn'])
    tot = 1-(1-y[8].max())*(1-fcn)
    print(f"   {lab:40s} {100*y[8].max():7.1f} {100*tot:8.1f} {y[10].min():5.3f} "
          f"{y[12].max():5.1f} {y[21].max():6.3f} {y[34,-1]:6.3f}")
print("  -> the occupancies compose on ONE axis, 1-(1-fCO)(1-fCN).  Oxygen cannot displace")
print("     cyanide and cobalamin cannot displace CO, and the antidote's value collapses")
print("     with delay: full at 5 min, none by 50.")

print("\n" + "="*80)
print("RESULT G.  Sensitivity of the 90-day cognitive score (severe, O2 mask)")
print("="*80)
p, t, y, F = sim(PPM_S, 60, **MASK); ref = Q(t, y, F)['Cog']
sens = []
for k in ['kappa','kon_cco','koff_cco','Kc','Ko','ATPthr','kXOon','kROS','kadh','kMPO',
          'kLPO','kAd','kSpread','kAdClr','theta','Tprol','Tdeath','kDemy','Hb','VA',
          'DL','MHald','f_wshed','Km_o2','CBF0','hbo_adh','CMRO2','kRepair','shunt']:
    v = p0[k]; o = []
    for f in (0.8, 1.25):
        pp, tt, yy, FF = sim(PPM_S, 60, npt=600, **{k: v*f}, **MASK)
        o.append(Q(tt, yy, FF)['Cog'])
    sens.append((k, abs(o[1]-o[0])/max(1e-9, abs(ref)), o[0], o[1]))
sens.sort(key=lambda z: -z[1])
print(f"   {'parameter':12s} {'|dCog|/Cog':>11s}   {'-20%':>7s} {'+25%':>7s}")
for k, s, lo, hi in sens[:18]:
    print(f"   {k:12s} {s:11.4f}   {lo:7.3f} {hi:7.3f}")
print(f"   (reference Cog = {ref:.3f})")
zero = [k for k, s, _, _ in sens if s < 1e-6]
print(f"   no influence on the 90-day score: {', '.join(zero) if zero else 'none'}")
