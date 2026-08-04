#!/usr/bin/env python3
# =============================================================================
# essential-tremor / et_verify.py
#
# Independent, dependency-free re-implementation of the Essential Tremor QSP
# model (et_mrgsolve_model.R).  Pure-Python RK4.  No numpy, no scipy.
#
# Purpose: this file is NOT a convenience wrapper around the R model.  It is a
# from-scratch second implementation whose job is to DISAGREE with the R model
# if either one is wrong.  Every number quoted in README.md and in the commit
# message comes from running this file.
#
#   python3 et_verify.py            # run all scenarios, print tables
#   python3 et_verify.py --quick    # coarser dt, for a fast smoke test
#
# Core formulation (see et_qsp_model.dot cluster ⑧):
#   a tremor is a LIMIT CYCLE, not a level.
#     G_total = G0*(1+PROG) * [ w_C*Phi_C + w_P*Phi_P ]
#     Phi_C   = (a_O*phi_olive + a_R*phi_cblthal) * phi_thal * phi_ctx
#     Phi_P   = phi_spindle * phi_nmj
#     mu      = G_total - 1                      (bifurcation parameter)
#     dr/dt   = (1/tauA) * ( r*(mu - beta*r^2) + eps )
#     r*      = sqrt(mu/beta)  for mu>0 ;  ~0 for mu<0
#   Amplitude comes from mu (gain).  Frequency comes from tau_loop (delay).
#   The therapeutic CEILING comes from topology: the olive is one of two
#   PARALLEL branches; the Vim relay is a SERIES element.
# =============================================================================

import math, sys

PI = math.pi
EXP = math.exp
SQ = math.sqrt
LOG10 = math.log10

# =============================================================================
# PARAMETERS
# =============================================================================
def P0():
    p = {}

    # ---- oscillator core -----------------------------------------------------
    p['G0']    = 1.60     # untreated loop gain, moderate ET (mild 1.15, severe 6.0)
    p['WC']    = 0.60     # central loop weight
    p['WP']    = 0.40     # peripheral stretch-reflex loop weight  (WC+WP=1)
    p['AO']    = 0.35     # olivary share of the central drive  <-- the T-type ceiling
    p['AR']    = 0.65     # cerebello-thalamic resonance share   (AO+AR=1)
    p['BETA']  = 1.0      # Hopf saturation coefficient
    p['TAUA']  = 0.35     # h, tremor-envelope time constant (~20 min)
    p['EPS']   = 1.0e-4   # noise seeding of the limit cycle
    p['PROG0'] = 0.0
    p['KPROG'] = 5.71e-6  # /h  -> G0 doubles in 20 y
    p['KEXC']  = 0.0      # hypothesis switch: oscillation -> cerebellar injury

    # ---- effector gain shares (which effectors oscillate at all) ------------
    p['HDG']   = 0.55     # head/neck effector gain share (0.85 = head-tremor phenotype)
    p['VXG']   = 0.45     # laryngeal effector gain share
    p['WC_HD'] = 0.85     # neck loop is mostly central...
    p['WP_HD'] = 0.15     # ...so beta2 spindle block helps it much less

    # ---- limb / effector biomechanics --------------------------------------
    p['J0']    = 0.00256  # kg m^2, hand about wrist
    p['KST']   = 7.30     # N m /rad, wrist stiffness  -> f0 = 8.50 Hz unloaded
    p['ZETA0'] = 0.35     # damping RATIO, held constant under load.
    # defect found by this file: holding the damping COEFFICIENT fixed made
    # zeta fall as 1/sqrt(J), so +500 g sharpened the resonance and DOUBLED ET
    # amplitude (1.91 -> 3.70 cm).  Loading provokes co-contraction, which
    # raises B with K; a constant ratio is the defensible choice.
    p['BD']    = 0.0955   # N m s /rad (documentation only; see ZETA0)
    p['MLOAD'] = 0.0      # kg added mass (weight-loading test)
    p['LLOAD'] = 0.10     # m, moment arm of the added mass
    p['KAMP']  = 1.902    # cm per unit r per unit |H|  (calibration)
    p['KNZ']   = 0.021    # cm, physiological-tremor noise gain
    p['KCAT']  = 12.0     # adrenergic amplification of physiological tremor
    p['FNZ']   = 0.55     # beta2-dependent fraction of resting physiological tremor
    p['J0_HD'] = 0.0250   # kg m^2, head about C1-C2
    p['KST_HD']= 8.80     # -> f0_head = 2.99 Hz
    p['BD_HD'] = 0.329    # -> zeta_head = 0.35
    p['KAMP_HD']= 3.10    # degrees per unit r per unit |H|
    p['NECKBTX']= 0.0

    # ---- loop delay / frequency --------------------------------------------
    p['TAU0']  = 0.18182  # s, total loop delay -> 5.50 Hz
    p['AGE']   = 60.0
    p['KAGE']  = 0.110    # per decade above 50 -> -0.055 Hz/y (Elble longitudinal)
    p['KPC']   = 0.42     # Purkinje loss lengthens the loop
    p['KLOADN']= 0.012    # mass loading barely touches the CENTRAL delay
    p['KVISD'] = 0.055    # visuomotor feedback adds delay during intention

    # ---- peripheral loop / beta2 -------------------------------------------
    p['FB2']   = 0.60     # beta2-dependent fraction of peripheral loop gain
    p['ADR0']  = 0.20     # nM resting adrenaline equivalent
    # defect found by this file: with KD_AG=30 the resting spindle beta2
    # occupancy was 0.0066 while carrying 60% of peripheral loop gain, so RAG
    # had ~150-fold headroom and a thyrotoxic "EPT" patient came out with mu=1.41
    # -- i.e. a limit cycle, which is exactly what EPT is NOT.  Putting the
    # resting system at 36% occupancy caps RAG at 1/OCC_AG0 = 2.75.
    p['KD_AG'] = 0.35     # nM, adrenaline at spindle beta2
    p['STRESS']= 0.0      # multiplier on adrenergic drive
    p['CAFF']  = 0.0
    p['THYRO'] = 0.0
    p['SALB']  = 0.0      # nM beta2-agonist equivalent
    p['SF']    = 3.0      # NMJ safety factor
    p['HN']    = 2.5      # supralinearity of release on SNAP-25

    # ---- disease / thalamic terms ------------------------------------------
    p['DNDIS0']= 0.0
    p['PCINT0']= 1.0
    p['KDEGPC']= 1.427e-6 # /h  -> PC integrity 1.00 -> 0.75 over 20 y
    p['KG_OL'] = 0.70     # GABA-A potentiation reaching the olive
    p['KREBF'] = 2.20     # acute counter-adaptation -> rebound gain
    p['KREBS'] = 1.20     # chronic counter-adaptation -> baseline shift
    # defect found by this file: with a single symmetric tau the adaptation
    # never outlived the ethanol, so the rebound was +0.1% -- i.e. absent.
    # Acute tolerance rises fast (Mellanby) and decays slowly; that asymmetry
    # IS the rebound.
    p['KAF']   = 0.90; p['TAUF_ON'] = 1.0; p['TAUF_OFF'] = 5.0
    p['KAS']   = 0.90; p['TAUS'] = 720.0      # slow adaptation (30 d)
    p['HARM']  = 0.0      # harmaline drive on the olive
    p['KRR']   = 3.70     # thalamic re-routing gain (habituation)
    p['REROUTE_MAX'] = 0.45
    p['TAU_RR']= 13140.0  # h  (1.5 y)

    # ---- Vim lesion / DBS ---------------------------------------------------
    p['LESION']= 0.0; p['VLES'] = 0.0
    p['V50L']  = 45.0; p['HL'] = 2.5          # mm^3, lesion efficacy
    p['V50A']  = 260.0; p['AXL'] = 1.00       # mm^3, lesion ataxia
    p['DBSON'] = 0.0; p['FSTIM'] = 130.0; p['VTA'] = 250.0
    p['F50D']  = 80.0; p['HD'] = 4.0; p['V50D'] = 90.0
    p['EMAX_DBS'] = 0.90
    p['KENT']  = 0.30; p['F50E'] = 60.0       # low-frequency entrainment
    p['BILAT'] = 0.0

    # ---- propranolol --------------------------------------------------------
    p['MW_PRP']=259.3; p['FU_PRP']=0.10
    p['KA_PRP']=0.35; p['V1_PRP']=250.0; p['V2_PRP']=300.0
    p['CL_PRP']=60.0; p['Q_PRP']=60.0
    p['F0_PRP']=0.20; p['FMX_PRP']=0.25; p['FD50_PRP']=150.0
    p['KPUU_PRP']=1.50; p['KE0_PRP']=1.20
    p['KI_PRP_B1']=1.80; p['KI_PRP_B2']=0.60
    p['EMAX_PRPC']=0.25; p['EC50_PRPC']=150.0   # non-beta central component (nM)
    p['KB2REG']=0.55; p['TAU_B2REG']=336.0

    # ---- atenolol (beta1-selective, low Kp: the negative control) ----------
    p['MW_ATN']=266.3; p['FU_ATN']=0.95
    p['KA_ATN']=0.80; p['V_ATN']=70.0; p['CL_ATN']=10.0; p['F_ATN']=0.50
    p['KI_ATN_B1']=30.0; p['KI_ATN_B2']=1000.0

    # ---- nadolol (non-selective, peripherally restricted) ------------------
    p['MW_NAD']=309.4; p['FU_NAD']=0.75
    p['KA_NAD']=0.40; p['V_NAD']=140.0; p['CL_NAD']=18.0; p['F_NAD']=0.32
    p['KI_NAD_B1']=3.00; p['KI_NAD_B2']=1.20

    # ---- primidone / phenobarbital / PEMA ----------------------------------
    p['MW_PRM']=218.25; p['MW_PB']=232.24; p['MW_PEM']=190.20
    p['KA_PRM']=1.50; p['F_PRM']=0.92; p['V_PRM']=40.0; p['CL_PRM']=2.80
    p['FM_PB']=0.25; p['FM_PEM']=0.45
    p['V_PB']=45.0;  p['CL_PB']=0.32
    p['V_PEM']=45.0; p['CL_PEM']=2.00
    p['KP_PRM']=0.90; p['KP_PB']=0.70; p['KE0_PRM']=1.00; p['KE0_PB']=0.35
    p['EMAX_PRM']=0.62; p['EC50_PRM']=28.0     # uM, GABA-A (parent = active moiety)
    p['EMAX_PB'] =0.55; p['EC50_PB'] =220.0    # uM
    p['EMAX_PEM']=0.10; p['EC50_PEM']=400.0
    p['EMAX_NACH']=0.12; p['EC50_NACH']=30.0   # cortical Na-channel term
    p['KIND']=1.30; p['TAU_IND']=336.0

    # ---- topiramate / gabapentin / benzodiazepine ---------------------------
    p['MW_TOP']=339.4; p['KA_TOP']=1.00; p['F_TOP']=0.80
    p['V_TOP']=60.0; p['CL_TOP']=1.40; p['KP_TOP']=0.90
    p['EMAX_TOPC']=0.28; p['EC50_TOPC']=14.0
    p['EMAX_TOPG']=0.15; p['EC50_TOPG']=30.0
    p['MW_GBP']=171.2; p['KA_GBP']=1.20; p['V_GBP']=60.0; p['CL_GBP']=10.0
    p['KP_GBP']=0.15; p['EMAX_GBP']=0.18; p['EC50_GBP']=3.0
    p['EMAX_BZD']=0.30; p['BZDLEV']=0.0

    # ---- ethanol / 1-octanol ------------------------------------------------
    p['V_ETH']=47.6; p['KA_ETH']=4.00; p['VMAX_ETH']=7.14; p['KM_ETH']=0.08
    p['KE0_ETH']=6.00; p['KP_ETH']=1.00
    p['EMAX_ETHG']=0.65; p['EC50_ETHG']=0.45   # g/L, GABA-A
    p['EMAX_ETHO']=0.75; p['EC50_ETHO']=0.30   # g/L, olivary decoupling
    p['EMAX_ETHC']=0.10; p['EC50_ETHC']=0.90   # g/L, cortical
    p['KPCTX_ETH']=1.00; p['EC50_INT']=0.80
    p['MW_OCT']=130.2; p['KA_OCT']=2.00; p['V_OCT']=50.0; p['CL_OCT']=25.0
    p['KP_OCT']=1.00
    p['EMAX_OCTG']=0.35; p['EC50_OCTG']=120.0  # uM
    p['EMAX_OCTO']=0.60; p['EC50_OCTO']=60.0   # uM
    p['KPCTX_OCT']=0.15; p['EC50_INTO']=250.0

    # ---- T-type calcium blocker --------------------------------------------
    p['MW_TTB']=400.0; p['KA_TTB']=1.20; p['F_TTB']=0.70
    p['V_TTB']=200.0; p['CL_TTB']=15.0; p['FU_TTB']=0.15
    p['KPUU_TTB']=1.50; p['KE0_TTB']=0.80
    p['EMAX_TT']=0.85; p['IC50_TT']=80.0       # nM free brain (state-dependent)

    # ---- botulinum toxin ----------------------------------------------------
    p['KDEG_BTX']=0.35; p['KCL_BTX']=0.00328   # per U.h
    p['KR_SNAP']=3.21e-4                        # SNAP-25 recovery, t1/2 90 d
    p['FSPILL']=0.15                            # guided 0.15 / unguided 0.45

    # ---- organ systems ------------------------------------------------------
    p['HR0']=72.0; p['EHR']=0.31; p['SBP0']=132.0; p['ESBP']=14.0
    p['FEV10']=3.20; p['ASTHMA']=0.0; p['EFEV_A']=0.30; p['EFEV_N']=0.04
    p['TAU_FEV']=2.0
    p['EC50_SED_PB']=90.0; p['EC50_SED_PRM']=45.0
    p['TAU_SEDTOL']=240.0; p['KSEDTOL']=0.80
    p['HCO30']=24.0; p['EHCO3']=5.0; p['EC50_HCO3']=20.0; p['TAU_HCO3']=48.0
    p['BW0']=78.0; p['EBW']=0.09; p['EC50_BW']=25.0; p['TAU_BW']=1440.0
    p['ECOG']=45.0; p['EC50_COG']=30.0; p['TAU_COG']=336.0
    p['BMD0']=1.00; p['KBMD']=1.9e-6
    p['ALT0']=25.0; p['EALT']=45.0; p['TAU_ALT']=720.0
    p['ETHCHR']=0.0

    # ---- task / state -------------------------------------------------------
    p['TASK_INT']=0.0; p['KVIS']=0.22; p['FATIGUE']=0.0; p['KFAT']=0.10
    p['ASYM']=0.80        # non-dominant / dominant amplitude ratio
    p['HDBASE']=0.0       # head-tremor phenotype flag handled via HDG
    return p

# state index map -------------------------------------------------------------
SN = ['A_PRPG','A_PRPC','A_PRPP','C_PRPB','A_ATNG','A_ATNC','A_NADG','A_NADC',
      'A_PRMG','A_PRMC','A_PBC','A_PEMC','C_PRMB','C_PBB',
      'A_TOPG','A_TOPC','A_GBPG','A_GBPC',
      'A_ETHG','A_ETHC','C_ETHB','A_OCTG','A_OCTC',
      'A_TTBG','A_TTBC','C_TTBB','A_BTXT','A_BTXG',
      'SNAPT','SNAPG','ADAPTF','ADAPTS','B2REG','IND','SEDTOL',
      'PCINT','DNDIS','REROUTE','PROG',
      'R_UL','R_HD','R_VX',
      'FEV1','HCO3','BW','BMD','ALT','COG']
IX = {n:i for i,n in enumerate(SN)}
NS = len(SN)

def init_state(p):
    y = [0.0]*NS
    y[IX['SNAPT']]=1.0; y[IX['SNAPG']]=1.0
    y[IX['PCINT']]=p['PCINT0']; y[IX['DNDIS']]=p['DNDIS0']
    y[IX['PROG']]=p['PROG0']
    y[IX['FEV1']]=p['FEV10']; y[IX['HCO3']]=p['HCO30']
    y[IX['BW']]=p['BW0']; y[IX['BMD']]=p['BMD0']; y[IX['ALT']]=p['ALT0']
    # seat the tremor envelope at its untreated fixed point
    a = alg(y, p)
    y[IX['R_UL']] = SQ(max(a['MU'],0.0)/p['BETA']) if a['MU']>0 else 1e-3
    y[IX['R_HD']] = SQ(max(a['MU_HD'],0.0)/p['BETA']) if a['MU_HD']>0 else 1e-3
    y[IX['R_VX']] = SQ(max(a['MU_VX'],0.0)/p['BETA']) if a['MU_VX']>0 else 1e-3
    return y

# =============================================================================
# ALGEBRA  (mirrors $MAIN + the derived block of the mrgsolve model)
# =============================================================================
def hmag(f, f0, z):
    r = f/f0
    return 1.0/SQ((1.0-r*r)**2 + (2.0*z*r)**2)

def alg(y, p):
    a = {}
    g = y.__getitem__
    # ---- concentrations ----------------------------------------------------
    C_PRP = g(IX['A_PRPC'])/p['V1_PRP']
    CF_PRP = C_PRP*p['FU_PRP']*1e6/p['MW_PRP']              # nM free
    CB_PRP = g(IX['C_PRPB'])                                 # nM free brain
    C_ATN = g(IX['A_ATNC'])/p['V_ATN']
    CF_ATN = C_ATN*p['FU_ATN']*1e6/p['MW_ATN']
    C_NAD = g(IX['A_NADC'])/p['V_NAD']
    CF_NAD = C_NAD*p['FU_NAD']*1e6/p['MW_NAD']
    C_PRM = g(IX['A_PRMC'])/p['V_PRM']; CU_PRM = C_PRM*1e3/p['MW_PRM']   # uM
    C_PB  = g(IX['A_PBC'])/p['V_PB'];   CU_PB  = C_PB*1e3/p['MW_PB']
    C_PEM = g(IX['A_PEMC'])/p['V_PEM']; CU_PEM = C_PEM*1e3/p['MW_PEM']
    CB_PRM = g(IX['C_PRMB']); CB_PB = g(IX['C_PBB'])
    C_TOP = g(IX['A_TOPC'])/p['V_TOP']; CU_TOP = C_TOP*1e3/p['MW_TOP']
    CB_TOP = p['KP_TOP']*CU_TOP
    C_GBP = g(IX['A_GBPC'])/p['V_GBP']; CU_GBP = C_GBP*1e3/p['MW_GBP']
    CB_GBP = p['KP_GBP']*CU_GBP
    C_ETH = g(IX['A_ETHC'])/p['V_ETH']; CB_ETH = g(IX['C_ETHB'])
    C_OCT = g(IX['A_OCTC'])/p['V_OCT']; CU_OCT = C_OCT*1e3/p['MW_OCT']
    CB_OCT = p['KP_OCT']*CU_OCT
    C_TTB = g(IX['A_TTBC'])/p['V_TTB']; CB_TTB = g(IX['C_TTBB'])
    a.update(C_PRP=C_PRP, CF_PRP=CF_PRP, CB_PRP=CB_PRP, C_ATN=C_ATN,
             CF_ATN=CF_ATN, C_NAD=C_NAD, CF_NAD=CF_NAD, C_PRM=C_PRM,
             CU_PRM=CU_PRM, C_PB=C_PB, CU_PB=CU_PB, CU_PEM=CU_PEM,
             CB_PRM=CB_PRM, CB_PB=CB_PB, C_TOP=C_TOP, CB_TOP=CB_TOP,
             CB_GBP=CB_GBP, C_ETH=C_ETH, CB_ETH=CB_ETH, CB_OCT=CB_OCT,
             C_TTB=C_TTB, CB_TTB=CB_TTB)

    # ---- beta-adrenoceptor occupancies ------------------------------------
    IQ_B2 = CF_PRP/p['KI_PRP_B2'] + CF_ATN/p['KI_ATN_B2'] + CF_NAD/p['KI_NAD_B2']
    IQ_B1 = CF_PRP/p['KI_PRP_B1'] + CF_ATN/p['KI_ATN_B1'] + CF_NAD/p['KI_NAD_B1']
    OCCB2 = IQ_B2/(1.0+IQ_B2); OCCB1 = IQ_B1/(1.0+IQ_B1)
    # defect found by this file: putting B2REG on the agonist made atenolol make
    # tremor WORSE on treatment (+2.5%).  Up-regulation is a receptor-number
    # effect, so it multiplies the beta2-mediated GAIN term and only shows up as
    # withdrawal rebound once the antagonist has washed out.
    AG = p['ADR0']*(1.0+p['STRESS']+p['CAFF']+p['THYRO']) + p['SALB']
    OCC_AG  = AG/(p['KD_AG']*(1.0+IQ_B2) + AG)
    OCC_AG0 = p['ADR0']/(p['KD_AG'] + p['ADR0'])
    RAG = OCC_AG/OCC_AG0
    a.update(IQ_B2=IQ_B2, OCCB2=OCCB2, OCCB1=OCCB1, OCC_AG=OCC_AG, RAG=RAG)

    # ---- peripheral loop --------------------------------------------------
    PHI_SPIN = 1.0 - p['FB2'] + p['FB2']*(1.0+g(IX['B2REG']))*RAG
    den = 1.0 - EXP(-p['SF'])
    PHI_NMJ = (1.0-EXP(-p['SF']*g(IX['SNAPT'])**p['HN']))/den
    GRIP    = 100.0*(1.0-EXP(-p['SF']*g(IX['SNAPG'])**p['HN']))/den
    PHI_P   = PHI_SPIN*PHI_NMJ
    a.update(PHI_SPIN=PHI_SPIN, PHI_NMJ=PHI_NMJ, PHI_P=PHI_P, GRIP=GRIP)

    # ---- GABA-A potentiation ----------------------------------------------
    P_PRM = p['EMAX_PRM']*CB_PRM/(p['EC50_PRM']+CB_PRM)
    P_PB  = p['EMAX_PB'] *CB_PB /(p['EC50_PB'] +CB_PB)
    P_PEM = p['EMAX_PEM']*CU_PEM/(p['EC50_PEM']+CU_PEM)
    P_ETH = p['EMAX_ETHG']*CB_ETH/(p['EC50_ETHG']+CB_ETH)
    P_OCT = p['EMAX_OCTG']*CB_OCT/(p['EC50_OCTG']+CB_OCT)
    P_TOP = p['EMAX_TOPG']*CB_TOP/(p['EC50_TOPG']+CB_TOP)
    P_BZD = p['EMAX_BZD']*p['BZDLEV']
    P_RAW = min(0.92, P_PRM+P_PB+P_PEM+P_ETH+P_OCT+P_TOP+P_BZD)
    P_EFF = P_RAW/(1.0+g(IX['ADAPTF'])+g(IX['ADAPTS']))
    # The rebound is the UNOPPOSED part of the adaptation.  While the drug is
    # present the adaptation is chasing it from below and cancels; once the drug
    # clears, what is left is excess excitability.  Writing REB as
    # 1+KREBF*ADAPTF instead destroyed the acute suppression (PHI_CBL went
    # above 1 at peak ethanol) -- defect found by this file.
    P_ETHR = p['EMAX_ETHG']*CB_ETH/(p['EC50_ETHG']+CB_ETH)
    REB   = 1.0 + p['KREBF']*max(0.0, g(IX['ADAPTF'])-p['KAF']*P_RAW) \
                + p['KREBS']*max(0.0, g(IX['ADAPTS'])-p['KAS']*P_ETHR)
    a.update(P_PRM=P_PRM, P_PB=P_PB, P_RAW=P_RAW, P_EFF=P_EFF, REB=REB)

    # ---- olivary branch (PARALLEL) ---------------------------------------
    BLK_TT = p['EMAX_TT']*CB_TTB/(p['IC50_TT']+CB_TTB)
    ETH_OL = p['EMAX_ETHO']*CB_ETH/(p['EC50_ETHO']+CB_ETH) \
           + p['EMAX_OCTO']*CB_OCT/(p['EC50_OCTO']+CB_OCT)
    ETH_OL = min(ETH_OL, 0.95)
    PHI_OL = (1.0-BLK_TT)*(1.0-ETH_OL)*(1.0-p['KG_OL']*P_EFF) \
             *(1.0+p['KG_OL']*(REB-1.0))*(1.0+p['HARM'])
    # ---- cerebello-thalamic branch (PARALLEL) ----------------------------
    PHI_CBL = (1.0+g(IX['DNDIS']))*REB*(1.0-P_EFF) \
              *(1.0-p['EMAX_PRPC']*CB_PRP/(p['EC50_PRPC']+CB_PRP))
    # ---- thalamic relay (SERIES) -----------------------------------------
    VL = p['VLES']
    LES_EFF = p['LESION']*(VL**p['HL']/(p['V50L']**p['HL']+VL**p['HL'])) if VL>0 else 0.0
    DBSF = p['FSTIM']**p['HD']/(p['F50D']**p['HD']+p['FSTIM']**p['HD'])
    DBSV = p['VTA']/(p['V50D']+p['VTA'])
    DBSB = p['EMAX_DBS']*DBSF*DBSV*p['DBSON']
    ENTR = 1.0 + p['KENT']*max(0.0,(p['F50E']-p['FSTIM'])/p['F50E'])*p['DBSON']
    PHI_TH = max(0.02, (1.0-0.92*LES_EFF)*(1.0-DBSB)*ENTR
                       *(1.0+p['KRR']*g(IX['REROUTE'])))
    # ---- cortical ---------------------------------------------------------
    PHI_CTX = (1.0-p['EMAX_TOPC']*CB_TOP/(p['EC50_TOPC']+CB_TOP)) \
            * (1.0-p['EMAX_NACH']*CB_PRM/(p['EC50_NACH']+CB_PRM)) \
            * (1.0-p['EMAX_GBP']*CB_GBP/(p['EC50_GBP']+CB_GBP)) \
            * (1.0-p['EMAX_ETHC']*CB_ETH/(p['EC50_ETHC']+CB_ETH)) \
            * (1.0+p['KFAT']*p['FATIGUE'])
    PHI_C = (p['AO']*PHI_OL + p['AR']*PHI_CBL)*PHI_TH*PHI_CTX
    a.update(BLK_TT=BLK_TT, PHI_OL=PHI_OL, PHI_CBL=PHI_CBL, PHI_TH=PHI_TH,
             PHI_CTX=PHI_CTX, PHI_C=PHI_C, LES_EFF=LES_EFF, DBSB=DBSB)

    # ---- loop gain and bifurcation ---------------------------------------
    G0T = p['G0']*(1.0+g(IX['PROG']))
    GTOT = G0T*(p['WC']*PHI_C + p['WP']*PHI_P)*(1.0+p['KVIS']*p['TASK_INT'])
    MU = GTOT - 1.0
    GHD = G0T*p['HDG']*(p['WC_HD']*PHI_C + p['WP_HD']*PHI_P) \
          *(1.0-0.55*p['NECKBTX'])
    MU_HD = GHD - 1.0
    GVX = G0T*p['VXG']*(0.90*PHI_C + 0.10*PHI_P)
    MU_VX = GVX - 1.0
    a.update(G0T=G0T, GTOT=GTOT, MU=MU, MU_HD=MU_HD, MU_VX=MU_VX)

    # ---- frequency and mechanics -----------------------------------------
    JREL = p['MLOAD']*p['LLOAD']**2/p['J0']
    TAU_L = p['TAU0']*(1.0 + p['KAGE']*(p['AGE']-50.0)/10.0
                        + p['KPC']*(1.0-g(IX['PCINT']))
                        + p['KVISD']*p['TASK_INT'] + p['KLOADN']*JREL)
    FNEUR = 1.0/TAU_L
    JTOT = p['J0'] + p['MLOAD']*p['LLOAD']**2
    F0 = (1.0/(2.0*PI))*SQ(p['KST']/JTOT)
    ZETA = p['ZETA0']
    HN_ = hmag(FNEUR, F0, ZETA); HP_ = hmag(F0, F0, ZETA)
    # defect found by this file: NZ = 1+KCAT*(RAG-1) went to -0.73 on
    # propranolol, i.e. a negative physiological tremor.  Physiological tremor
    # has a beta2-independent component; only FNZ of it is adrenergic.
    NZ = (1.0 + p['KCAT']*(RAG-1.0)) if RAG > 1.0 else max(0.30, 1.0-p['FNZ']*(1.0-RAG))
    A_PHYS = p['KNZ']*NZ*HP_
    A_LC = p['KAMP']*g(IX['R_UL'])*HN_
    A_UL = SQ(A_LC*A_LC + A_PHYS*A_PHYS)
    F_OBS = (A_LC*FNEUR + A_PHYS*F0)/(A_LC+A_PHYS)
    F0_HD = (1.0/(2.0*PI))*SQ(p['KST_HD']/p['J0_HD'])
    Z_HD = p['ZETA0']
    A_HD = p['KAMP_HD']*g(IX['R_HD'])*hmag(FNEUR, F0_HD, Z_HD)
    A_VX = 6.0*g(IX['R_VX'])
    a.update(FNEUR=FNEUR, F0=F0, ZETA=ZETA, A_LC=A_LC, A_PHYS=A_PHYS,
             A_UL=A_UL, F_OBS=F_OBS, A_HD=A_HD, A_VX=A_VX, NZ=NZ)

    # ---- Elble log transform -> rating scales ----------------------------
    def T(x): return min(4.0, max(0.0, 2.0+2.0*LOG10(max(x,1e-6))))
    TR = T(A_UL); TL = T(A_UL*p['ASYM'])
    THD = min(4.0, max(0.0, 2.0+2.0*LOG10(max(A_HD/3.0,1e-6))))
    TVX = min(4.0, max(0.0, 1.2*LOG10(max(A_VX,1e-6))+1.5)) if A_VX>0 else 0.0
    UL = 3.0*TR + 3.0*TL
    SPI = min(4.0,max(0.0,0.85*TR+0.15)) + min(4.0,max(0.0,0.85*TL+0.15))
    HW = min(4.0,max(0.0,0.90*TR))
    LL = 2.0*min(4.0,max(0.0,0.55*TR))
    FACE = min(4.0,0.5*THD); TONG = min(4.0,0.5*THD)
    STAND = min(4.0,max(0.0,0.60*TR+0.40*THD))
    TETRAS_PS = UL+SPI+HW+LL+THD+FACE+TONG+TVX+STAND
    GRIPLOSS = (100.0-GRIP)/100.0
    ADLf = min(1.0, max(0.0, 0.24*TR + 0.35*GRIPLOSS))
    TETRAS_ADL = 48.0*ADLf
    FTM = 2.25*TETRAS_PS
    SPIRAL = min(4.0,max(0.0,0.85*TR+0.15))
    BAINF = min(10.0, 2.5*TR)
    a.update(T_UL=TR, TETRAS_PS=TETRAS_PS, TETRAS_ADL=TETRAS_ADL, FTM=FTM,
             SPIRAL=SPIRAL, BAINF=BAINF)

    # ---- organ systems ----------------------------------------------------
    HR = p['HR0']*(1.0-p['EHR']*OCCB1)
    SBP = p['SBP0']-p['ESBP']*OCCB1
    FEV_SS = p['FEV10']*(1.0-p['EFEV_A']*p['ASTHMA']*OCCB2-p['EFEV_N']*OCCB2)
    SED_RAW = 100.0*min(1.0, 0.55*CB_PB/(p['EC50_SED_PB']+CB_PB)
                          + 0.45*CB_PRM/(p['EC50_SED_PRM']+CB_PRM)
                          + 0.60*CB_ETH/(0.90+CB_ETH)
                          + 0.30*CB_TOP/(30.0+CB_TOP)
                          + 0.35*CB_GBP/(4.0+CB_GBP)
                          + 0.25*CB_PRP/(300.0+CB_PRP))
    SED = SED_RAW/(1.0+g(IX['SEDTOL']))
    ATAX = 100.0*min(1.0, 0.50*CB_PB/(120.0+CB_PB)
                        + 0.60*CB_ETH/(1.20+CB_ETH)
                        + p['AXL']*(VL/(p['V50A']+VL) if VL>0 else 0.0)
                        + 0.35*DBSV*p['DBSON']*(1.0+p['BILAT'])/2.0)
    INTOX = 100.0*min(1.0, CB_ETH*p['KPCTX_ETH']/(p['EC50_INT']+CB_ETH*p['KPCTX_ETH'])
                        + CB_OCT*p['KPCTX_OCT']/(p['EC50_INTO']+CB_OCT*p['KPCTX_OCT']))
    QUEST = min(100.0, 100.0*(0.52*ADLf + 0.13*SED/100.0 + 0.12*ATAX/100.0
                    + 0.09*g(IX['COG'])/100.0 + 0.08*GRIPLOSS + 0.06*INTOX/100.0))
    a.update(HR=HR, SBP=SBP, FEV_SS=FEV_SS, SED_RAW=SED_RAW, SED=SED,
             ATAX=ATAX, INTOX=INTOX, QUEST=QUEST, GRIPLOSS=GRIPLOSS)
    return a

# =============================================================================
# DERIVATIVES
# =============================================================================
def derivs(t, y, p, a=None):
    if a is None: a = alg(y, p)
    d = [0.0]*NS
    g = y.__getitem__

    # propranolol 2-cpt + effect site
    d[IX['A_PRPG']] = -p['KA_PRP']*g(IX['A_PRPG'])
    d[IX['A_PRPC']] = p['KA_PRP']*g(IX['A_PRPG'])*p['_FPRP'] \
                    - p['CL_PRP']/p['V1_PRP']*g(IX['A_PRPC']) \
                    - p['Q_PRP']*(g(IX['A_PRPC'])/p['V1_PRP']-g(IX['A_PRPP'])/p['V2_PRP'])
    d[IX['A_PRPP']] = p['Q_PRP']*(g(IX['A_PRPC'])/p['V1_PRP']-g(IX['A_PRPP'])/p['V2_PRP'])
    d[IX['C_PRPB']] = p['KE0_PRP']*(p['KPUU_PRP']*a['CF_PRP']-g(IX['C_PRPB']))
    # atenolol / nadolol
    d[IX['A_ATNG']] = -p['KA_ATN']*g(IX['A_ATNG'])
    d[IX['A_ATNC']] = p['KA_ATN']*g(IX['A_ATNG'])*p['F_ATN'] - p['CL_ATN']/p['V_ATN']*g(IX['A_ATNC'])
    d[IX['A_NADG']] = -p['KA_NAD']*g(IX['A_NADG'])
    d[IX['A_NADC']] = p['KA_NAD']*g(IX['A_NADG'])*p['F_NAD'] - p['CL_NAD']/p['V_NAD']*g(IX['A_NADC'])
    # primidone -> PB, PEMA
    kel = p['CL_PRM']/p['V_PRM']
    d[IX['A_PRMG']] = -p['KA_PRM']*g(IX['A_PRMG'])
    d[IX['A_PRMC']] = p['KA_PRM']*g(IX['A_PRMG'])*p['F_PRM'] - kel*g(IX['A_PRMC'])
    d[IX['A_PBC']]  = p['FM_PB']*kel*g(IX['A_PRMC'])*(p['MW_PB']/p['MW_PRM']) \
                    - p['CL_PB']/p['V_PB']*g(IX['A_PBC'])
    d[IX['A_PEMC']] = p['FM_PEM']*kel*g(IX['A_PRMC'])*(p['MW_PEM']/p['MW_PRM']) \
                    - p['CL_PEM']/p['V_PEM']*g(IX['A_PEMC'])
    d[IX['C_PRMB']] = p['KE0_PRM']*(p['KP_PRM']*a['CU_PRM']-g(IX['C_PRMB']))
    d[IX['C_PBB']]  = p['KE0_PB'] *(p['KP_PB'] *a['CU_PB'] -g(IX['C_PBB']))
    # topiramate (cleared faster when induced) / gabapentin
    d[IX['A_TOPG']] = -p['KA_TOP']*g(IX['A_TOPG'])
    d[IX['A_TOPC']] = p['KA_TOP']*g(IX['A_TOPG'])*p['F_TOP'] \
                    - p['CL_TOP']*(1.0+p['KIND']*g(IX['IND']))/p['V_TOP']*g(IX['A_TOPC'])
    d[IX['A_GBPG']] = -p['KA_GBP']*g(IX['A_GBPG'])
    d[IX['A_GBPC']] = p['KA_GBP']*g(IX['A_GBPG']) - p['CL_GBP']/p['V_GBP']*g(IX['A_GBPC'])
    # ethanol (Michaelis-Menten) / octanol
    d[IX['A_ETHG']] = -p['KA_ETH']*g(IX['A_ETHG'])
    d[IX['A_ETHC']] = p['KA_ETH']*g(IX['A_ETHG']) \
                    - p['VMAX_ETH']*a['C_ETH']/(p['KM_ETH']+a['C_ETH'])
    d[IX['C_ETHB']] = p['KE0_ETH']*(p['KP_ETH']*a['C_ETH']-g(IX['C_ETHB']))
    d[IX['A_OCTG']] = -p['KA_OCT']*g(IX['A_OCTG'])
    d[IX['A_OCTC']] = p['KA_OCT']*g(IX['A_OCTG']) - p['CL_OCT']/p['V_OCT']*g(IX['A_OCTC'])
    # T-type blocker
    d[IX['A_TTBG']] = -p['KA_TTB']*g(IX['A_TTBG'])
    d[IX['A_TTBC']] = p['KA_TTB']*g(IX['A_TTBG'])*p['F_TTB'] - p['CL_TTB']/p['V_TTB']*g(IX['A_TTBC'])
    CF_TTB = a['C_TTB']*p['FU_TTB']*1e6/p['MW_TTB']
    d[IX['C_TTBB']] = p['KE0_TTB']*(p['KPUU_TTB']*CF_TTB-g(IX['C_TTBB']))
    # botulinum depots + SNAP-25 pools
    d[IX['A_BTXT']] = -p['KDEG_BTX']*g(IX['A_BTXT'])
    d[IX['A_BTXG']] = -p['KDEG_BTX']*g(IX['A_BTXG'])
    d[IX['SNAPT']] = p['KR_SNAP']*(1.0-g(IX['SNAPT'])) - p['KCL_BTX']*g(IX['A_BTXT'])*g(IX['SNAPT'])
    d[IX['SNAPG']] = p['KR_SNAP']*(1.0-g(IX['SNAPG'])) - p['KCL_BTX']*g(IX['A_BTXG'])*g(IX['SNAPG'])
    # adaptation states
    _tgt = p['KAF']*a['P_RAW']
    d[IX['ADAPTF']] = (_tgt - g(IX['ADAPTF']))/(p['TAUF_ON'] if _tgt > g(IX['ADAPTF'])
                                                else p['TAUF_OFF'])
    P_ETHONLY = p['EMAX_ETHG']*a['CB_ETH']/(p['EC50_ETHG']+a['CB_ETH'])
    d[IX['ADAPTS']] = (p['KAS']*P_ETHONLY - g(IX['ADAPTS']))/p['TAUS']
    d[IX['B2REG']]  = (p['KB2REG']*a['OCCB2'] - g(IX['B2REG']))/p['TAU_B2REG']
    d[IX['IND']]    = (a['CU_PB']/(a['CU_PB']+60.0) - g(IX['IND']))/p['TAU_IND']
    d[IX['SEDTOL']] = (p['KSEDTOL']*a['SED_RAW']/100.0 - g(IX['SEDTOL']))/p['TAU_SEDTOL']
    # disease
    rsq = min(g(IX['R_UL'])**2, 1e6)
    d[IX['PCINT']] = -p['KDEGPC']*(1.0+p['KEXC']*rsq)
    d[IX['DNDIS']] = (2.4*(1.0-g(IX['PCINT'])) - g(IX['DNDIS']))/8760.0
    supp = 1.0 if (a['LES_EFF']>0.2 or a['DBSB']>0.2) else 0.0
    d[IX['REROUTE']] = (p['REROUTE_MAX']*supp - g(IX['REROUTE']))/p['TAU_RR']
    d[IX['PROG']] = p['KPROG']*(1.0+p['KEXC']*rsq)
    # tremor envelopes (Hopf normal form, slow-envelope reduction)
    # NOTE (defect found by this file, 1st run): with dr/dt = (1/tauA) r (mu-b r^2)
    # the linearised relaxation time near r* is tauA/(2 mu), so it SHRINKS as the
    # patient gets more severe.  At G0=6 (mu=5) that is 0.035 h and RK4 at any
    # step usable over 24 weeks diverges (G0=12 returned NaN).  Dividing by
    # (1+|mu|) leaves the fixed point r*=sqrt(mu/beta) untouched and caps the
    # relaxation time at tauA/2 for all severities.
    iT = 1.0/p['TAUA']; B = p['BETA']; e = p['EPS']
    rU=g(IX['R_UL']); rH=g(IX['R_HD']); rV=g(IX['R_VX'])
    d[IX['R_UL']] = iT*(rU*(a['MU']   - B*rU*rU) + e)/(1.0+abs(a['MU']))
    d[IX['R_HD']] = iT*(rH*(a['MU_HD']- B*rH*rH) + e)/(1.0+abs(a['MU_HD']))
    d[IX['R_VX']] = iT*(rV*(a['MU_VX']- B*rV*rV) + e)/(1.0+abs(a['MU_VX']))
    # organs
    d[IX['FEV1']] = (a['FEV_SS']-g(IX['FEV1']))/p['TAU_FEV']
    HCSS = p['HCO30']-p['EHCO3']*a['CB_TOP']/(p['EC50_HCO3']+a['CB_TOP'])
    d[IX['HCO3']] = (HCSS-g(IX['HCO3']))/p['TAU_HCO3']
    BWSS = p['BW0']*(1.0-p['EBW']*a['CB_TOP']/(p['EC50_BW']+a['CB_TOP']))
    d[IX['BW']] = (BWSS-g(IX['BW']))/p['TAU_BW']
    d[IX['BMD']] = -p['KBMD']*g(IX['IND'])*g(IX['BMD'])
    ALTSS = p['ALT0']+p['EALT']*p['ETHCHR']
    d[IX['ALT']] = (ALTSS-g(IX['ALT']))/p['TAU_ALT']
    CGSS = p['ECOG']*a['CB_TOP']/(p['EC50_COG']+a['CB_TOP'])
    d[IX['COG']] = (CGSS-g(IX['COG']))/p['TAU_COG']
    return d

# =============================================================================
# INTEGRATOR
# =============================================================================
def run(p, tend, dt=0.1, events=(), sample=None, y0=None, qss=None):
    """events: list of (time, state_name, amount). sample: list of output times.

    qss: quasi-steady-state the three tremor envelopes.  Defect found by this
    file: multi-year runs need a coarse dt, but the envelope relaxation time is
    TAUA/2 = 0.175 h, so at dt = 1 h the RK4 step is 5.7x the stability limit
    and the envelope oscillated numerically -- the 20-year progression run
    returned a NON-MONOTONIC amplitude (2.16 cm at yr 2, 2.09 at yr 5) while mu
    rose monotonically from 0.61 to 2.26.  For any dt above ~0.2 h the envelope
    is quasi-static anyway, so pin it to its fixed point instead.  mrgsolve is
    not affected (LSODA adapts its own step); this is a verifier-side fix.
    """
    if qss is None: qss = (dt >= 0.20)
    p = dict(p)
    p['_FPRP'] = p['F0_PRP']  # set per-dose below
    y = init_state(p) if y0 is None else list(y0)
    ev = sorted(events, key=lambda e: e[0]); ei = 0
    if sample is None: sample = [tend]
    sample = sorted(sample); si = 0
    out = []
    t = 0.0
    n = int(round(tend/dt))
    for k in range(n+1):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            _, nm, amt = ev[ei]
            if nm == 'A_PRPG':
                p['_FPRP'] = p['F0_PRP'] + p['FMX_PRP']*amt/(p['FD50_PRP']+amt)
            y[IX[nm]] += amt; ei += 1
        while si < len(sample) and sample[si] <= t + 1e-9:
            a = alg(y, p); a['t'] = t
            for nm in ('R_UL','SNAPT','SNAPG','PCINT','PROG','FEV1','HCO3',
                       'BW','BMD','ALT','COG','ADAPTF','ADAPTS','IND','REROUTE'):
                a[nm] = y[IX[nm]]
            out.append(a); si += 1
        if k == n: break
        k1 = derivs(t, y, p)
        y2 = [y[i]+0.5*dt*k1[i] for i in range(NS)]
        k2 = derivs(t+0.5*dt, y2, p)
        y3 = [y[i]+0.5*dt*k2[i] for i in range(NS)]
        k3 = derivs(t+0.5*dt, y3, p)
        y4 = [y[i]+dt*k3[i] for i in range(NS)]
        k4 = derivs(t+dt, y4, p)
        for i in range(NS):
            y[i] += dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
        y[IX['SNAPT']] = min(1.0, max(1e-4, y[IX['SNAPT']]))
        y[IX['SNAPG']] = min(1.0, max(1e-4, y[IX['SNAPG']]))
        y[IX['R_UL']] = max(0.0, y[IX['R_UL']])
        y[IX['R_HD']] = max(0.0, y[IX['R_HD']])
        y[IX['R_VX']] = max(0.0, y[IX['R_VX']])
        t = (k+1)*dt
        if qss:
            aq = alg(y, p)
            y[IX['R_UL']] = SQ(aq['MU']   /p['BETA']) if aq['MU']    > 0 else 1e-3
            y[IX['R_HD']] = SQ(aq['MU_HD']/p['BETA']) if aq['MU_HD'] > 0 else 1e-3
            y[IX['R_VX']] = SQ(aq['MU_VX']/p['BETA']) if aq['MU_VX'] > 0 else 1e-3
    return out, y

def run_avg(p, tend, dt=0.1, events=(), extra_last=24.0, qss=None):
    """Run to tend and return the MEAN of every output over the final dosing
    interval.  Sampling a single trough (the first version of this file did)
    understates a once-daily drug's effect: propranolol 160 mg read -31% at
    trough vs -41% on the interval mean.  Trial assessments are not troughs."""
    smp = [tend-extra_last + i*1.0 for i in range(int(extra_last)+1)]
    smp = [x for x in smp if x >= 0.0]
    o, y = run(p, tend, dt, events=events, sample=smp, qss=qss)
    keys = [k for k in o[-1] if isinstance(o[-1][k], float)]
    m = {k: sum(z[k] for z in o)/len(o) for k in keys}
    m['_trough'] = o[-1]
    return m, y


def qd(dose, tend, state, interval=24.0, start=0.0):
    """repeated dosing event list"""
    ev = []; t = start
    while t < tend - 1e-9:
        ev.append((t, state, dose)); t += interval
    return ev

# =============================================================================
# FAST-OSCILLATOR CROSS-CHECK
#   The slow-envelope reduction dr/dt = (1/tauA) r (mu - beta r^2) is only
#   legitimate if r* really equals the RMS amplitude of the underlying
#   oscillator.  Integrate the full 2-D Hopf oscillator at 5.5 Hz and compare.
# =============================================================================
def fast_oscillator_check(mu, f, tau_s=0.5, beta=1.0, secs=60.0, dt=2e-5):
    # defect found by this file: the first version used tauA in HOURS as the
    # envelope timescale, so over a 40 s window the oscillator never grew and
    # the check returned ~0.068 for every mu.  The physical envelope relaxes in
    # a few cycles; tau_s = 0.5 s.  TAUA in the reduced model is a numerical
    # choice, not this timescale.
    w = 2.0*PI*f
    x, yv = 0.02, 0.0
    n = int(secs/dt); s2 = 0.0; m = 0
    tstart = int(0.6*n)
    for k in range(n):
        r2 = x*x+yv*yv
        lam = (mu - beta*r2)/tau_s
        dx = lam*x - w*yv
        dy = lam*yv + w*x
        x += dt*dx; yv += dt*dy
        if k > tstart:
            s2 += x*x; m += 1
    rms = SQ(s2/m)
    return rms*SQ(2.0), SQ(max(mu,0.0)/beta)   # peak (=rms*sqrt2) vs r*

# =============================================================================
# SCENARIOS
# =============================================================================
W = 24*7.0
def base(**kw):
    p = P0(); p.update(kw); return p

def summarise(o, label):
    return dict(label=label, A=o['A_UL'], TP=o['TETRAS_PS'], MU=o['MU'],
                G=o['GTOT'], F=o['F_OBS'], HR=o['HR'], SED=o['SED'],
                GRIP=o['GRIP'], Q=o['QUEST'], ATAX=o['ATAX'], INTOX=o['INTOX'],
                FEV=o['FEV1'], ADL=o['TETRAS_ADL'])

def pct(new, old): return 100.0*(new-old)/old

def main():
    quick = '--quick' in sys.argv
    DT = 0.25 if quick else 0.10
    DTF = 0.02 if quick else 0.005
    R = {}
    print('='*78)
    print(' ESSENTIAL TREMOR QSP MODEL — independent Python verification')
    print('='*78)

    # -- envelope reduction validity -----------------------------------------
    print('\n[0] Slow-envelope reduction vs full 2-D oscillator')
    print('     mu      peak(full osc)   r*(envelope)   rel.err')
    for mu in (0.15, 0.60, 1.40, 5.00):
        pk, rs = fast_oscillator_check(mu, 5.5)
        print('   %5.2f      %10.4f     %10.4f     %6.2f%%'
              % (mu, pk, rs, 100*(pk-rs)/rs))

    # -- 1. untreated baselines by severity ----------------------------------
    print('\n[1] Untreated baselines (24 weeks, no drug)')
    print('   %-22s %8s %8s %8s %8s %8s' % ('phenotype','G','mu','A(cm)','TETRAS','f(Hz)'))
    for nm, G0 in (('mild G0=1.15',1.15),('moderate G0=1.60',1.60),
                   ('severe G0=6.0',6.00),('very severe G0=12',12.0)):
        m,_ = run_avg(base(G0=G0), 8*W, DT)
        R['base_'+nm.split()[0]] = m
        print('   %-22s %8.3f %8.3f %8.3f %8.1f %8.2f'
              % (nm, m['GTOT'], m['MU'], m['A_UL'], m['TETRAS_PS'], m['F_OBS']))
    B = R['base_moderate']

    # -- 2. beta blockers -----------------------------------------------------
    print('\n[2] beta-blockade, 24 weeks (moderate ET, baseline A=%.3f cm, TETRAS=%.1f)'
          % (B['A_UL'], B['TETRAS_PS']))
    print('   %-26s %7s %7s %7s %8s %7s %7s' %
          ('regimen','A(cm)','dA%','TETRAS','dTETRAS','HR','OCCb2'))
    betas = [('propranolol LA 60 mg',   'A_PRPG', 60.0),
             ('propranolol LA 120 mg',  'A_PRPG',120.0),
             ('propranolol LA 160 mg',  'A_PRPG',160.0),
             ('propranolol LA 240 mg',  'A_PRPG',240.0),
             ('propranolol LA 320 mg',  'A_PRPG',320.0),
             ('atenolol 100 mg',        'A_ATNG',100.0),
             ('nadolol 120 mg',         'A_NADG',120.0)]
    for nm, st, d in betas:
        z,_ = run_avg(base(), 8*W, DT, events=qd(d, 8*W, st)); R[nm] = z
        print('   %-26s %7.3f %7.1f %7.1f %8.2f %7.1f %7.3f'
              % (nm, z['A_UL'], pct(z['A_UL'],B['A_UL']), z['TETRAS_PS'],
                 z['TETRAS_PS']-B['TETRAS_PS'], z['HR'], z['OCCB2']))

    # -- 3. the log-scale artefact -------------------------------------------
    z = R['propranolol LA 160 mg']
    print('\n[3] The log-scale artefact (propranolol 160 mg)')
    print('   accelerometric amplitude : %+.1f %%' % pct(z['A_UL'],B['A_UL']))
    print('   TETRAS performance score : %+.2f points (%+.1f %% of baseline)'
          % (z['TETRAS_PS']-B['TETRAS_PS'], pct(z['TETRAS_PS'],B['TETRAS_PS'])))
    print('   single upper-limb item   : %+.2f points' % (z['T_UL']-B['T_UL']))
    print('   amplitude ratio needed for a 1.00-point item change: %.2f-fold'
          % (10**0.5))

    # -- 4. severity-dependence of response (the sqrt law) -------------------
    print('\n[4] Same drug, different severity — propranolol LA 160 mg')
    print('   %-14s %8s %8s %8s %9s %9s' % ('phenotype','A0','A_on','dA%','TETRAS0','dTETRAS'))
    for nm, G0 in (('mild',1.15),('moderate',1.60),('severe',6.00),('v.severe',12.0)):
        b = R['base_'+('mild' if nm=='mild' else 'moderate' if nm=='moderate'
                       else 'severe' if nm=='severe' else 'very')]
        z,_ = run_avg(base(G0=G0), 8*W, DT, events=qd(160.0, 8*W,'A_PRPG'))
        print('   %-14s %8.3f %8.3f %8.1f %9.1f %9.2f'
              % (nm, b['A_UL'], z['A_UL'], pct(z['A_UL'],b['A_UL']),
                 b['TETRAS_PS'], z['TETRAS_PS']-b['TETRAS_PS']))
        R['sev_'+nm] = z

    # -- 5. primidone ---------------------------------------------------------
    print('\n[5] Primidone — the parent molecule is the active moiety')
    print('   %-22s %8s %8s %8s %8s %8s %7s' %
          ('regimen','A(cm)','dA%','P_prim','P_PB','%parent','SED'))
    for nm, d in (('primidone 250 mg',250.0),('primidone 500 mg',500.0),
                  ('primidone 750 mg',750.0)):
        z,_ = run_avg(base(), 8*W, DT,
                      events=qd(d/3.0, 8*W,'A_PRMG', interval=8.0)); R[nm] = z
        fp = 100*z['P_PRM']/(z['P_PRM']+z['P_PB']+1e-12)
        print('   %-22s %8.3f %8.1f %8.4f %8.4f %8.1f %7.1f'
              % (nm, z['A_UL'], pct(z['A_UL'],B['A_UL']), z['P_PRM'],
                 z['P_PB'], fp, z['SED']))
    # day-1 vs steady state
    o,_ = run(base(), 8*W, DT, events=qd(250.0/3.0, 8*W,'A_PRMG', interval=8.0),
              sample=[2.0, 6.0, 12.0, 24.0, 72.0, 336.0, 8*W])
    print('   time-course on 250 mg/d (P_prim / P_PB / A):')
    for z in o:
        print('     t=%7.1f h   P_prim %.4f   P_PB %.4f   A %.3f cm  (dA %+.1f%%)'
              % (z['t'], z['P_PRM'], z['P_PB'], z['A_UL'], pct(z['A_UL'],B['A_UL'])))
    # phenobarbital level needed to match
    tgt = R['primidone 250 mg']['P_RAW']
    need_uM = P0()['EC50_PB']*tgt/max(P0()['EMAX_PB']-tgt,1e-9)
    print('   phenobarbital brain conc needed to match primidone 250 mg: %.0f uM'
          % need_uM)
    print('   -> plasma total %.0f uM = %.1f mg/L (therapeutic range 10-40 mg/L)'
          % (need_uM/P0()['KP_PB'], need_uM/P0()['KP_PB']*P0()['MW_PB']/1e3))

    # -- 6. combination -------------------------------------------------------
    print('\n[6] Combination: additive in gain, sublinear in amplitude')
    ev = qd(160.0, 8*W,'A_PRPG') + qd(250.0/3.0, 8*W,'A_PRMG', interval=8.0)
    C,_ = run_avg(base(), 8*W, DT, events=ev); R['combo']=C
    dA_p = pct(R['propranolol LA 160 mg']['A_UL'], B['A_UL'])
    dA_r = pct(R['primidone 250 mg']['A_UL'], B['A_UL'])
    dA_c = pct(C['A_UL'], B['A_UL'])
    print('   propranolol alone %+.1f%% | primidone alone %+.1f%% | sum %+.1f%%'
          % (dA_p, dA_r, dA_p+dA_r))
    print('   observed combination %+.1f%%  -> EXCESS over the sum %+.1f points'
          % (dA_c, dA_c-(dA_p+dA_r)))
    print('   (sqrt(mu) is CONCAVE, so pushing mu toward 0 makes the amplitude')
    print('    map steepest exactly where the second drug acts -> supra-additive)')
    print('   gain: G %.3f -> %.3f (prop) / %.3f (prim) / %.3f (both)'
          % (B['GTOT'], R['propranolol LA 160 mg']['GTOT'],
             R['primidone 250 mg']['GTOT'], C['GTOT']))
    R['_combo_nums'] = (dA_p, dA_r, dA_c)
    return R, DT, DTF, B

_MAINOUT = None


# =============================================================================
# PART 2 — ethanol, T-type ceiling, surgery, botulinum, diagnostics, progression
# =============================================================================
def main2(R, DT, DTF, B):
    print('\n[7] Ethanol 2 standard drinks (28 g) — suppression then rebound')
    p = base(); ev = [(0.0,'A_ETHG',28.0)]
    smp = [0,0.25,0.5,0.75,1.0,1.5,2,3,4,5,6,8,10,12,16,24,36]
    o,_ = run(p, 40.0, DTF, events=ev, sample=smp)
    print('   %6s %8s %8s %8s %8s %8s %8s' %
          ('t(h)','BAC(g/L)','P_eff','ADAPTF','G','A(cm)','dA%'))
    A0 = o[0]['A_UL']
    nadir=(0,1e9); reb=(0,-1e9)
    for z in o:
        print('   %6.2f %8.3f %8.3f %8.3f %8.3f %8.3f %+8.1f'
              % (z['t'], z['C_ETH'], z['P_EFF'], z['ADAPTF'], z['GTOT'],
                 z['A_UL'], pct(z['A_UL'],A0)))
        if z['A_UL']<nadir[1]: nadir=(z['t'],z['A_UL'])
        if z['t']>4 and z['A_UL']>reb[1]: reb=(z['t'],z['A_UL'])
    print('   nadir  %.3f cm (%+.1f%%) at t=%.2f h' % (nadir[1],pct(nadir[1],A0),nadir[0]))
    print('   rebound peak %.3f cm (%+.1f%%) at t=%.1f h' % (reb[1],pct(reb[1],A0),reb[0]))

    print('\n[8] Ethanol dose-dependence of the rebound')
    print('   %8s %10s %10s %10s %10s' % ('drinks','nadir dA%','t_nadir','rebound dA%','INTOX'))
    for nd in (1,2,3,4):
        o,_ = run(base(), 48.0, DTF, events=[(0.0,'A_ETHG',14.0*nd)],
                  sample=[i*0.25 for i in range(193)])
        A0=o[0]['A_UL']
        na=min(o,key=lambda z:z['A_UL']); rb=max([z for z in o if z['t']>6],key=lambda z:z['A_UL'])
        ix=max(o,key=lambda z:z['INTOX'])
        print('   %8d %10.1f %10.2f %10.1f %10.1f'
              % (nd, pct(na['A_UL'],A0), na['t'], pct(rb['A_UL'],A0), ix['INTOX']))

    print('\n[9] Daily ethanol 3 drinks x 90 d — the self-medication trap')
    ev = qd(42.0, 90*24.0, 'A_ETHG', interval=24.0, start=19.0)
    smp = [24*d+18.0 for d in (1,7,14,30,60,89)] + [24*d+3.0 for d in (1,7,30,89)]
    o,_ = run(base(), 91*24.0, 0.02, events=ev, sample=sorted(smp))
    print('   %8s %8s %8s %8s %8s %8s' % ('day','clock','ADAPTS','G','A(cm)','dA%'))
    for z in o:
        print('   %8.0f %8.1f %8.3f %8.3f %8.3f %+8.1f'
              % (z['t']//24, z['t']%24, z['ADAPTS'], z['GTOT'], z['A_UL'],
                 pct(z['A_UL'],B['A_UL'])))

    print('\n[10] 1-octanol vs ethanol — tremor benefit per unit intoxication')
    print('   %-22s %9s %9s %9s' % ('agent','peak dA%','peak INTOX','ratio'))
    for nm, st, d in (('ethanol 2 drinks','A_ETHG',28.0),
                      ('1-octanol 8 mg/kg','A_OCTG',560.0),
                      ('1-octanol 16 mg/kg','A_OCTG',1120.0)):
        o,_ = run(base(), 24.0, DTF, events=[(0.0,st,d)],
                  sample=[i*0.25 for i in range(97)])
        na=min(o,key=lambda z:z['A_UL']); ix=max(o,key=lambda z:z['INTOX'])
        dA=pct(na['A_UL'],B['A_UL'])
        print('   %-22s %9.1f %9.1f %9.2f' % (nm, dA, ix['INTOX'], abs(dA)/max(ix['INTOX'],0.01)))

    print('\n[11] T-type calcium blocker — the ceiling is TOPOLOGICAL (a_O)')
    print('   %-34s %8s %8s %8s %8s' % ('setting','phi_ol','Phi_C','G','dA%'))
    hum,_ = run_avg(base(), 4*W, DT, events=qd(100.0, 4*W,'A_TTBG'))
    print('   %-34s %8.3f %8.3f %8.3f %8.1f'
          % ('human ET, 100 mg/d', hum['PHI_OL'], hum['PHI_C'], hum['GTOT'],
             pct(hum['A_UL'],B['A_UL'])))
    hum2,_ = run_avg(base(EMAX_TT=1.0, IC50_TT=1.0), 4*W, DT, events=qd(100.0,4*W,'A_TTBG'))
    print('   %-34s %8.3f %8.3f %8.3f %8.1f'
          % ('human ET, PERFECT Cav3 block', hum2['PHI_OL'], hum2['PHI_C'],
             hum2['GTOT'], pct(hum2['A_UL'],B['A_UL'])))
    # harmaline model: the oscillation is purely olivary and purely central
    pr = base(AO=1.0, AR=0.0, WC=0.95, WP=0.05, G0=2.2, HARM=0.0)
    rb,_ = run_avg(pr, 4*W, DT); R['rat_base']=rb
    rt,_ = run_avg(pr, 4*W, DT, events=qd(100.0, 4*W,'A_TTBG'))
    print('   %-34s %8.3f %8.3f %8.3f %8.1f'
          % ('harmaline rat (a_O=1), 100 mg/d', rt['PHI_OL'], rt['PHI_C'],
             rt['GTOT'], pct(rt['A_UL'],rb['A_UL'])))
    print('   -> same drug, same occupancy, ONE parameter (a_O) different')
    print('   ceiling sweep over a_O at perfect Cav3 block:')
    for ao in (0.20,0.35,0.50,0.65,0.80,1.00):
        pb2 = base(AO=ao, AR=1.0-ao); b2,_ = run_avg(pb2, 2*W, DT)
        z2,_ = run_avg(base(AO=ao, AR=1.0-ao, EMAX_TT=1.0, IC50_TT=1.0),
                       2*W, DT, events=qd(100.0,2*W,'A_TTBG'))
        print('     a_O=%.2f  ceiling dA = %+6.1f %%' % (ao, pct(z2['A_UL'],b2['A_UL'])))

    print('\n[12] MRgFUS thalamotomy — one lesion volume, two outcomes')
    print('   %8s %8s %8s %8s %8s %8s' % ('V(mm3)','LES_EFF','G','A(cm)','dA%','ATAX'))
    for v in (20,40,60,90,120,180,250,400):
        z,_ = run_avg(base(LESION=1.0, VLES=float(v)), 4*W, DT)
        print('   %8d %8.3f %8.3f %8.3f %8.1f %8.1f'
              % (v, z['LES_EFF'], z['GTOT'], z['A_UL'], pct(z['A_UL'],B['A_UL']), z['ATAX']))

    print('\n[13] Vim DBS — frequency is not a dial, it is a switch')
    print('   %8s %8s %8s %8s %8s %8s' % ('f(Hz)','DBSB','PHI_TH','G','A(cm)','dA%'))
    for f in (10,30,50,60,80,100,130,185,250):
        z,_ = run_avg(base(DBSON=1.0, FSTIM=float(f), VTA=250.0), 4*W, DT)
        print('   %8d %8.3f %8.3f %8.3f %8.3f %8.1f'
              % (f, z['DBSB'], z['PHI_TH'], z['GTOT'], z['A_UL'], pct(z['A_UL'],B['A_UL'])))

    print('\n[14] Habituation — re-routing capacity decides who relapses')
    print('   %10s %8s %8s %8s %8s' % ('REROUTE_max','yr','REROUTE','A(cm)','dA%'))
    for rmax in (0.45, 0.75, 0.95, 1.00):
        for yr in (0.0, 1.0, 3.0, 5.0):
            z,_ = run_avg(base(LESION=1.0, VLES=120.0, REROUTE_MAX=rmax),
                          max(4*W, yr*8760.0), 2.0 if yr>0 else DT)
            print('   %10.2f %8.1f %8.3f %8.3f %8.1f'
                  % (rmax, yr, z['REROUTE'], z['A_UL'], pct(z['A_UL'],B['A_UL'])))

    print('\n[15] Botulinum toxin — the window is spatial precision, not dose')
    print('   %-30s %8s %8s %8s %8s %8s' %
          ('injection','SNAP_T','SNAP_G','dA%','grip%','QUEST'))
    for nm, dose, fs in (('100 U, US/EMG-guided (f=0.15)',100.0,0.15),
                         ('100 U, unguided (f=0.45)',100.0,0.45),
                         ('50 U, unguided (f=0.45)',50.0,0.45),
                         ('150 U, guided (f=0.15)',150.0,0.15)):
        pp = base(FSPILL=fs)
        ev = [(0.0,'A_BTXT',dose*(1.0-fs)), (0.0,'A_BTXG',dose*fs)]
        o,_ = run(pp, 6*W, DT, events=ev, sample=[i*24.0 for i in range(1,43)])
        z = min(o, key=lambda q: q['A_UL'])
        print('   %-30s %8.3f %8.3f %8.1f %8.1f %8.1f'
              % (nm, z['SNAPT'], z['SNAPG'], pct(z['A_UL'],B['A_UL']),
                 z['GRIP'], z['QUEST']))

    print('\n[16] Topiramate and gabapentin')
    for nm, st, d, iv in (('topiramate 200 mg/d','A_TOPG',100.0,12.0),
                          ('topiramate 400 mg/d','A_TOPG',200.0,12.0),
                          ('gabapentin 1800 mg/d','A_GBPG',600.0,8.0)):
        z,_ = run_avg(base(), 8*W, DT, events=qd(d, 8*W, st, interval=iv))
        print('   %-22s dA %+6.1f %%  TETRAS %+5.2f  HCO3 %.1f  BW %.1f kg  COG %.0f/100'
              % (nm, pct(z['A_UL'],B['A_UL']), z['TETRAS_PS']-B['TETRAS_PS'],
                 z['HCO3'], z['BW'], z['COG']))

    print('\n[17] Weight-loading test — ET vs enhanced physiological tremor')
    print('   %-28s %8s %8s %8s %8s' % ('condition','MU','f0(Hz)','f_obs','A(cm)'))
    _LT = {}
    for nm, pp in (('ET, unloaded', base()),
                   ('ET, +500 g', base(MLOAD=0.5)),
                   ('EPT (thyrotoxic), unloaded', base(G0=0.42, THYRO=9.0)),
                   ('EPT (thyrotoxic), +500 g', base(G0=0.42, THYRO=9.0, MLOAD=0.5)),
                   ('EPT on propranolol 160', base(G0=0.42, THYRO=9.0))):
        ev = qd(160.0, 4*W,'A_PRPG') if 'propranolol' in nm else ()
        z,_ = run_avg(pp, 4*W, DT, events=ev)
        _LT[nm] = z
        print('   %-28s %8.3f %8.2f %8.2f %8.3f'
              % (nm, z['MU'], z['F0'], z['F_OBS'], z['A_UL']))
    print('   -> loading moves the ET peak by %+.2f Hz and the EPT peak by %+.2f Hz'
          % (_LT['ET, +500 g']['F_OBS'] - _LT['ET, unloaded']['F_OBS'],
             _LT['EPT (thyrotoxic), +500 g']['F_OBS']
             - _LT['EPT (thyrotoxic), unloaded']['F_OBS']))
    print('      (the standard clinical discriminator, reproduced from the equations)')

    print('\n[18] Head tremor — emerges when the neck effector crosses G=1')
    print('   %-26s %8s %8s %8s %8s' % ('phenotype / drug','MU_HD','A_HD(deg)','A_UL','dA_HD%'))
    for nm, hdg, g0 in (('HDG=0.55, G0=1.6',0.55,1.6),('HDG=0.85, G0=1.6',0.85,1.6),
                        ('HDG=0.85, G0=2.4',0.85,2.4),('HDG=0.85, G0=6.0',0.85,6.0)):
        b,_ = run_avg(base(HDG=hdg, G0=g0), 4*W, DT)
        z,_ = run_avg(base(HDG=hdg, G0=g0), 4*W, DT, events=qd(160.0,4*W,'A_PRPG'))
        dh = pct(z['A_HD'], b['A_HD']) if b['A_HD']>1e-3 else 0.0
        print('   %-26s %8.3f %8.3f %8.3f %8.1f'
              % (nm, b['MU_HD'], b['A_HD'], b['A_UL'], dh))
        if hdg==0.85 and g0==2.4:
            print('     arm on propranolol %+.1f%% vs head %+.1f%%  <- same drug'
                  % (pct(z['A_UL'],b['A_UL']), dh))

    print('\n[19] Adrenergic challenge (performance stress + caffeine)')
    print('   %-30s %8s %8s %8s' % ('condition','G','A(cm)','dA% vs rest'))
    for nm, pp, ev in (('rest', base(), ()),
                       ('stress x4 + caffeine', base(STRESS=3.0, CAFF=1.0), ()),
                       ('stress + propranolol 160', base(STRESS=3.0, CAFF=1.0),
                        qd(160.0,4*W,'A_PRPG')),
                       ('rest + propranolol 160', base(), qd(160.0,4*W,'A_PRPG'))):
        z,_ = run_avg(pp, 4*W, DT, events=ev)
        print('   %-30s %8.3f %8.3f %8.1f'
              % (nm, z['GTOT'], z['A_UL'], pct(z['A_UL'],B['A_UL'])))

    print('\n[20] The contraindication IS the efficacy (same beta2 occupancy)')
    print('   %-26s %8s %8s %8s %8s' % ('drug','OCCb2','dA%','FEV1(L)','dFEV1%'))
    for nm, st, d in (('propranolol 160','A_PRPG',160.0),('nadolol 120','A_NADG',120.0),
                      ('atenolol 100','A_ATNG',100.0)):
        z,_ = run_avg(base(ASTHMA=1.0), 4*W, DT, events=qd(d,4*W,st))
        print('   %-26s %8.3f %8.1f %8.3f %8.1f'
              % (nm, z['OCCB2'], pct(z['A_UL'],B['A_UL']), z['FEV1'],
                 pct(z['FEV1'], P0()['FEV10'])))

    print('\n[21] Progression over 20 years (sqrt law -> decelerating growth)')
    print('   %6s %8s %8s %8s %8s %8s' % ('yr','G','mu','A(cm)','TETRAS','dA/yr%'))
    prev=None
    for yr in (0,1,2,5,10,15,20):
        z,_ = run_avg(base(), max(4*W, yr*8760.0), 2.0 if yr>0 else DT)
        rate = 0.0 if prev is None or yr==0 else 100*(z['A_UL']/prev-1)/max(yr-pyr,1)
        print('   %6d %8.3f %8.3f %8.3f %8.1f %8.1f'
              % (yr, z['GTOT'], z['MU'], z['A_UL'], z['TETRAS_PS'], rate))
        prev=z['A_UL']; pyr=yr

    print('\n[21b] Is early suppression disease-modifying?  (KEXC hypothesis,')
    print('      default OFF -- reported as a hypothesis, NOT as a result)')
    for nm, kx, ev2 in (('KEXC 0, untreated', 0.0, ()),
                        ('KEXC 0.3, untreated', 0.30, ()),
                        ('KEXC 0.3, propranolol 10 y', 0.30,
                         qd(160.0, 10*8760.0,'A_PRPG'))):
        z,_ = run_avg(base(KEXC=kx), 10*8760.0, 2.0, events=ev2)
        print('   %-30s PROG %.3f  G %.3f  A %.3f cm  TETRAS %.1f'
              % (nm, z['PROG'], z['GTOT'], z['A_UL'], z['TETRAS_PS']))

    print('\n[22] Isolating propranolol\'s non-beta central component')
    for d in (160.0, 320.0):
        a1,_ = run_avg(base(), 4*W, DT, events=qd(d,4*W,'A_PRPG'))
        a2,_ = run_avg(base(EMAX_PRPC=0.0), 4*W, DT, events=qd(d,4*W,'A_PRPG'))
        print('   %.0f mg: with central %+.1f%% | peripheral only %+.1f%% -> central adds %.1f pts'
              % (d, pct(a1['A_UL'],B['A_UL']), pct(a2['A_UL'],B['A_UL']),
                 pct(a1['A_UL'],B['A_UL'])-pct(a2['A_UL'],B['A_UL'])))
    print('\n' + '='*78)

if __name__ == '__main__':
    _R, _DT, _DTF, _B = main()
    main2(_R, _DT, _DTF, _B)
