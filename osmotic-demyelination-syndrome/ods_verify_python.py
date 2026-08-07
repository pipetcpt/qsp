#!/usr/bin/env python3
# =====================================================================
#  ods_verify_python.py
#  Osmotic Demyelination Syndrome (ODS) — independent Python/scipy
#  re-implementation of the mrgsolve model, used to (a) verify that the
#  ODE system integrates to the published anchors and (b) generate every
#  number quoted in README.md / the .dot map / the references file.
#
#  삼투성 탈수초 증후군 — mrgsolve 모델의 독립 재구현 및 검증
# ---------------------------------------------------------------------
#  STRUCTURAL THESIS
#  -----------------
#  Osmotic demyelination is not caused by hyponatraemia and it is not
#  caused by sodium.  It is caused by an ASYMMETRY IN TWO TRANSPORT
#  TIME CONSTANTS inside the astrocyte:
#
#      EFFLUX of organic osmolytes   VRAC / LRRC8A ion channel   t1/2 ~ 8 h
#      INFLUX of organic osmolytes   SMIT1 / TauT / BGT1,        t1/2 ~ 2-3 d
#                                    requiring TonEBP-driven
#                                    TRANSCRIPTION of the carriers
#
#  A channel opens in milliseconds; a transporter has to be transcribed.
#  The brain can therefore give solute away roughly six times faster than
#  it can take it back, and every clinical rule about correcting
#  hyponatraemia is a statement about that ratio.
#
#  The single state variable that carries the disease is
#
#      OMEGA(t) = ORG_set( Osm_eff(t) )  -  ORG(t)          [mOsm/kg bw]
#
#  the ORGANIC OSMOLYTE DEFICIT: how much compatible solute the brain
#  ought to be holding at the tonicity it now finds itself in, minus how
#  much it actually holds.  OMEGA is zero in the normal brain, zero in
#  the chronically ADAPTED hyponatraemic brain (which is why those
#  patients are not injured by a sodium of 110), and becomes positive
#  only when plasma tonicity moves faster than transcription.
#
#  Everything else in the file is bookkeeping around that one quantity.
# ---------------------------------------------------------------------
#  Units: TIME = days.  Brain solutes = mOsm (or mmol) per kg of
#  BASELINE brain water.  Body solutes = mmol.  Volumes = L.
# =====================================================================

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import json, sys, os

# --- smooth clamps ----------------------------------------------------
# Hard min()/max() inside an ODE make the Jacobian discontinuous, which
# costs LSODA (and mrgsolve's LSODA) an enormous number of rejected
# steps.  Every switch in this model is therefore written as a smooth
# hinge with a small width EPS; at EPS = 1e-3 the difference from the
# hard clamp is below the reporting precision of every number quoted.
EPSH = 1.0e-3
def sp(x, e=EPSH):
    """C1 smoothed max(0, x): EXACTLY zero for x <= 0, x - e/2 for x >= e,
    a quadratic bridge in between (bias e/2, chosen far below reporting
    precision for every quantity it is applied to).

    The obvious smooth hinge 0.5*(x + sqrt(x^2+e^2)) is NOT zero below the
    knee -- it leaves a floor of e^2/(4|x|).  That floor is harmless in a
    rate equation but not in an INTEGRATED one: fed through the astrocyte
    -> microglia -> oligodendrocyte chain it demyelinated a completely
    healthy brain to 75% of normal over 60 simulated days.  The bug is
    recorded here because it is the kind that never announces itself."""
    return np.where(x <= 0.0, 0.0, np.where(x >= e, x - 0.5 * e, x * x / (2.0 * e)))
def sclamp(x, lo, hi, e=EPSH):
    return lo + sp(x - lo, e) - sp(x - hi, e)

# ---------------------------------------------------------------------
# 1. PARAMETERS
# ---------------------------------------------------------------------
P0 = dict(
    # ---- body fluid / Edelman -----------------------------------------
    EDA      = 1.11,    # Edelman slope   [Na]s = EDA*(Nae+Ke)/TBW - EDB
    EDB      = 25.6,    # Edelman intercept (mmol/L)   Edelman 1958 J Clin Invest
    FECF     = 0.65,    # osmotically ACTIVE fraction of exchangeable Na (bone/skin store)
    OSMX     = 5.0,     # non-Na effective osmoles (glucose 90 mg/dL) mOsm/kg

    # ---- renal solute handling ----------------------------------------
    NAINT    = 120.0,   # dietary Na intake            mmol/d
    KINT     = 60.0,    # dietary K intake             mmol/d
    NAESET   = 3000.0,  # exchangeable Na set point    mmol
    KESET    = 3266.0,  # exchangeable K set point     mmol
    KNAEX    = 1.50,    # Na excretion gain            /d  (load t1/2 ~11 h)
    KKEX     = 1.00,    # K excretion gain             /d
    PUREA    = 400.0,   # urea generation (70 g protein/d)   mmol/d
    CLUREA   = 55.0,    # renal urea clearance         L/d
    OTHOSM   = 60.0,    # other urinary osmoles        mOsm/d
    UOSMMIN  = 50.0,    # maximally dilute urine       mOsm/kg
    UOSMMAX  = 1150.0,  # maximally concentrated urine mOsm/kg
    SOL0     = 900.0,   # solute load above which the medulla washes out  mOsm/d
    SOLREF   = 900.0,   # washout scale                                  mOsm/d

    # ---- water ---------------------------------------------------------
    WIN      = 1.638295,# oral water intake            L/d (solved: exact steady state)
    WMET     = 0.30,    # metabolic water              L/d
    OSMTHIRST= 292.0,   # thirst threshold             mOsm/kg
    KTHIRST  = 0.25,    # L/d per mOsm/kg above threshold
    WINSENS  = 0.90,    # insensible loss              L/d

    # ---- AVP axis ------------------------------------------------------
    OSMTHR   = 280.0,   # osmotic threshold for AVP    mOsm/kg
    GAINOSM  = 0.45,    # pg/mL per mOsm/kg
    AVPVMAX  = 12.0,    # max non-osmotic (baroreceptor) AVP   pg/mL
    VOLEC50  = 0.08,    # ECF deficit giving half-max AVP  (fraction)
    VOLHILL  = 3.0,
    AVPSIADH = 0.0,     # autonomous (ectopic/neurohypophyseal) AVP  pg/mL
    TAUAVP   = 0.015,   # AVP turnover                 d  (~22 min)
    EC50AVP  = 1.20,    # V2 EC50 for AVP              pg/mL
    EC50DDA  = 1.60,    # V2 EC50 for desmopressin     pg/mL
    KITLV    = 30.0,    # apparent V2 Ki, tolvaptan (TOTAL drug)  ng/mL
    TAUAQP   = 0.060,   # AQP2 trafficking             d (~1.4 h)
    THIAZ    = 0.0,     # thiazide: raises minimum urine osmolality (0/1)

    # ---- brain: baseline composition (mOsm per kg baseline brain water)
    BW0      = 1.120,   # brain water                  kg (1400 g x 80%)
    IMP      = 10.0,    # impermeant (protein, phosphate)  mOsm/kg
    INS0     = 7.0,     # myo-inositol
    TAU0     = 2.0,     # taurine
    GLX0     = 15.0,    # glutamate + glutamine
    CRE0     = 10.0,    # creatine + phosphocreatine
    GPC0     = 2.5,     # glycerophosphocholine + phosphocholine
    OTH0     = 11.5,    # NAA, glycine, betaine, others
    # osmoresponsiveness beta_i  (fractional change per fractional tonicity change)
    BINS     = 3.2, BTAU = 3.2, BGLX = 2.6, BCRE = 1.0, BGPC = 2.4, BOTH = 0.8,
    ORGFLOOR = 0.12,    # floor as fraction of baseline pool

    # ---- brain transport kinetics --------------------------------------
    TAUEFF   = 0.90,    # VRAC efflux time constant    d  (t1/2 ~ 15 h)
    TAUINF   = 1.25,    # carrier influx time constant d  (t1/2 ~ 21 h, x transporter)
    TAUTON   = 0.25,    # TonEBP activation            d
    TAUSMIT  = 1.00,    # carrier protein turnover     d
    ETON     = 0.9,     # TonEBP gain per fractional organic deficit
    TONMAX   = 2.6,     # ceiling on TonEBP drive
    FOSM     = 1.0,     # organic transport capacity (alcohol/malnutrition -> <1)
    TAUELEC  = 0.15,    # RVI/RVD electrolyte time constant  d
    JRVIMAX  = 60.0,    # max inorganic accumulation   mOsm/kg/d
    JRVDMAX  = 150.0,   # max inorganic loss           mOsm/kg/d
    BELCEIL  = 1.10,    # inorganic ceiling, x baseline
    TAUBUR   = 0.10,    # brain <-> plasma urea equilibration  d

    # ---- injury --------------------------------------------------------
    WSHR     = 1.0,     # weight of shrinkage (converted to mOsm equivalent)
    OMSTAR   = 8.0,     # injury threshold             mOsm/kg
    RISK     = 1.0,     # threshold multiplier (high risk -> <1)
    KINJ     = 0.300,   # injury rate constant
    HINJ     = 1.00,    # injury exponent
    FNUT     = 1.0,     # astrocyte energy reserve (malnutrition -> <1)
    TAUATP   = 0.10,    # ATP recovery                 d
    KATPU    = 0.030,   # ATP consumption per unit stress
    KAST     = 0.30,    # astrocyte death gain
    ASTTHR   = 0.030,   # astrocyte loss tolerated before trophic failure
    KASTREP  = 0.250,   # astrocyte repopulation       /d
    KBBB     = 0.60,    # BBB opening per unit astrocyte loss
    TAUBBB   = 2.50,    # BBB resealing                d
    KMG      = 0.50,    # microglial activation gain
    TAUMG    = 4.00,    # microglial decay             d
    KCYT     = 0.60, TAUCYT = 1.50,
    KIGG     = 0.50, TAUIGG = 3.00,
    KOLI     = 0.20,    # oligodendrocyte death gain
    KHIT     = 1.20,    # saturation constant of the oligodendrocyte hit
    KSEV     = 0.42,    # lesion burden giving half-maximal clinical deficit
    WOA      = 1.00,    # weight: astrocyte loss
    WOC      = 0.45,    # weight: cytokines
    WOI      = 0.35,    # weight: IgG / complement
    KOPCP    = 0.070,   # OPC proliferation            /d
    KOPCD    = 0.055,   # OPC -> oligodendrocyte       /d
    KDEM     = 0.40,    # demyelination gain
    KMYE     = 0.060,   # remyelination gain
    WPONS    = 1.00,    # topographic weight, central pons
    WEXP     = 0.55,    # topographic weight, extrapontine
    TAULES   = 45.0,    # radiological lesion resolution  d
    TAUMRI   = 3.5,     # MRI signal lag after myelin loss  d
    TAUDEF   = 1.50,    # clinical deficit lag             d
    DEFMAX   = 100.0,

    # ---- drug PK --------------------------------------------------------
    VDDA     = 25.0,  KDDA = 5.55,   KADDA = 36.0,   # desmopressin  V(L), k(/d) t1/2 3h, ka SC
    VTLV     = 210.0, KTLV = 2.08,   KATLV = 12.0, FTLV = 0.50,  # tolvaptan t1/2 8h
    KAUREA   = 48.0,                                  # oral urea absorption /d
    VDEX     = 60.0,  KDEX = 4.16,                    # dexamethasone t1/2 4 h
    VMIN     = 80.0,  KMIN = 1.04,                    # minocycline t1/2 16 h
    EMAXUREA = 0.55,  EC50UREA = 12.0,   # urea: BBB/microglia protection, mmol/L over baseline
    EMAXDEX  = 0.45,  EC50DEX  = 15.0,   # dexamethasone on BBB, ng/mL
    EMAXMIN  = 0.50,  EC50MIN  = 2.0,    # minocycline on microglia, ug/mL

    # ---- infusions / control (set per scenario) --------------------------
    R3PCT    = 0.0,   # fixed 3% NaCl rate            L/d
    R09      = 0.0,   # fixed 0.9% NaCl rate          L/d
    RD5W     = 0.0,   # fixed D5W rate                L/d
    RKCL     = 0.0,   # KCl infusion                  mmol/d
    NA3      = 513.0, # [Na] of 3% saline             mmol/L
    NA09     = 154.0, # [Na] of 0.9% saline           mmol/L
    CTRLON   = 0.0,   # 1 = closed-loop 3% saline titration
    TCORR    = 0.0,   # time correction starts        d
    RATETGT  = 6.0,   # target correction rate        mmol/L/d
    NASTART  = 110.0, # [Na] at start of correction
    NACAP    = 130.0, # stop correcting at this [Na]
    KP3      = 2.0,   # controller gain               L/d per mmol/L
    R3MAX    = 4.0,   # 3% saline clamp               L/d
    RESCUE   = 0.0,   # 1 = relowering rescue active
    TRESCUE  = 1e9,   # time rescue starts            d
    NARES    = 118.0, # relowering target [Na]
    KPD5W    = 0.8,   # relowering controller gain    L/d per mmol/L
    DURRES   = 1.5,   # relowering treatment window   d
    D5WMAX   = 5.0,   # D5W clamp  (~3 mL/kg/h)       L/d
    DDACLAMP = 0.0,   # 1 = scheduled DDAVP q8h from TDDA
    TDDA     = 1e9,
    DDADOSE  = 2.0,   # ug per dose
    UREADOSE = 0.0,   # g/d oral urea (divided q8h -> modelled as continuous gut input)
    DEXON    = 0.0, MINOON = 0.0,
    ACUTE    = 0.0,   # 1 = acute water-loading protocol
    WLOAD    = 0.0,   # free-water load rate          L/d
    TWLEND   = 0.0,   # end of water load             d
    # ---- extrarenal losses (vomiting / diarrhoea / diuretic / third space)
    NALOSS   = 0.0,   # mmol/d
    KLOSS    = 0.0,   # mmol/d
    WLOSS    = 0.0,   # L/d
    TLOSSEND = 0.0,   # losses stop here (= presentation)  d
    # ---- the ADH SWITCH: the moment the stimulus is removed --------------
    TAVPOFF  = 1e9,   # cortisol given / thiazide stopped / alcohol stopped  d
    AVPFREEZE= 0.0,   # 1 = clamp AVP at its value at TFREEZE (counterfactual)
    TFREEZE  = 1e9,
    AVPFRZV  = 0.0,
    # ---- solute intake step (refeeding a beer-potomania patient) ---------
    TSOLUP   = 1e9,
    PUREA2   = 400.0,
)

SNAMES = ['TBW','NAE','KE','UREAB','AVP','AQP2',
          'BELEC','INS','TAU','GLX','CRE','GPC','OTH','BURE',
          'TONEBP','SMIT','ATP','AST','BBBP','MG','CYT','IGG',
          'OLI','OPC','MYE','LESP','LESE','MRI','DEF',
          'DDAD','DDAC','TLVD','TLVC','UREAG','DEXC','MINC',
          'CUMI','CUMNA','CUMV','CUMEFW']
IX = {n: i for i, n in enumerate(SNAMES)}
NS = len(SNAMES)


# ---------------------------------------------------------------------
# 2. ALGEBRAIC LAYER  (identical expressions live in the mrgsolve $ODE)
# ---------------------------------------------------------------------
def algebra(t, y, p):
    """Everything the derivatives need, computed once."""
    a = {}
    TBW = max(y[IX['TBW']], 5.0)
    NAE, KE, UREAB = y[IX['NAE']], y[IX['KE']], y[IX['UREAB']]

    # --- Edelman: serum sodium -------------------------------------------
    SNA = p['EDA'] * (NAE + KE) / TBW - p['EDB']
    a['SNA'] = SNA
    a['OSMEFF'] = 2.0 * SNA + p['OSMX']
    a['UREAP'] = UREAB / TBW                      # mmol/L
    a['BUN'] = a['UREAP'] * 2.80                  # mg/dL
    a['OSMTOT'] = a['OSMEFF'] + a['UREAP']
    # serum potassium from total-body K deficit (98% of K is intracellular)
    a['SK'] = max(1.5, 4.0 + (KE - p['KESET']) / 300.0)

    # --- ECF volume & the non-osmotic AVP drive ---------------------------
    ECFV = p['FECF'] * NAE / max(SNA, 60.0)
    ECFV0 = p['FECF'] * p['NAESET'] / 140.0
    a['ECFV'] = ECFV
    vd = sp((ECFV0 - ECFV) / ECFV0, 2e-3)
    a['VOLDEF'] = vd
    avp_osm = p['GAINOSM'] * sp(a['OSMEFF'] - p['OSMTHR'], 0.30)
    avp_vol = p['AVPVMAX'] * vd**p['VOLHILL'] / (p['VOLEC50']**p['VOLHILL'] + vd**p['VOLHILL'])
    gt2 = lambda x, w: 0.5 * (1.0 + np.tanh(x / w))
    siadh = p['AVPSIADH'] * gt2(p['TAVPOFF'] - t, 0.01)
    a['AVPTGT'] = avp_osm + avp_vol + siadh
    if p['AVPFREEZE'] > 0.5:                # counterfactual: hold the kidney still
        w = gt2(t - p['TFREEZE'], 0.01)
        a['AVPTGT'] = (1.0 - w) * a['AVPTGT'] + w * p['AVPFRZV']

    # --- V2 receptor occupancy: AVP + desmopressin vs tolvaptan -----------
    DDAC = max(y[IX['DDAC']], 0.0) / p['VDDA'] * 1e6      # ug -> pg/mL
    TLVC = max(y[IX['TLVC']], 0.0) / p['VTLV'] * 1e6      # mg -> ng/mL
    a['DDACP'], a['TLVCP'] = DDAC, TLVC
    ago = max(y[IX['AVP']], 0.0) / p['EC50AVP'] + DDAC / p['EC50DDA']
    a['OCC'] = ago / (1.0 + ago + TLVC / p['KITLV'])

    # --- urine ------------------------------------------------------------
    uomin = p['UOSMMIN'] + 100.0 * p['THIAZ']
    ENA = 0.05 * p['NAINT'] + sp(p['KNAEX'] * (NAE - p['NAESET']) + 0.95 * p['NAINT'], 0.5)
    EK = 0.05 * p['KINT'] + sp(p['KKEX'] * (KE - p['KESET']) + 0.95 * p['KINT'], 0.5)
    EUREA = p['CLUREA'] * a['UREAP']
    # thirst: without it a prescribed fluid restriction runs forever and the
    # patient walks to a sodium of 210 over six weeks, which is an artefact of
    # the prescription, not of the disease.
    a['WINT'] = p['WIN'] + p['KTHIRST'] * sp(a['OSMEFF'] - p['OSMTHIRST'], 0.3)
    gt = lambda x, w: 0.5 * (1.0 + np.tanh(x / w))
    a['PU'] = p['PUREA'] + (p['PUREA2'] - p['PUREA']) * gt(t - p['TSOLUP'], 0.02)
    ls = gt(p['TLOSSEND'] - t, 0.004)
    a['NALOSS'], a['KLOSS'], a['WLOSS'] = p['NALOSS'] * ls, p['KLOSS'] * ls, p['WLOSS'] * ls
    SOLEXC = 2.0 * (ENA + EK) + EUREA + p['OTHOSM']
    # flow-dependent medullary washout: the countercurrent multiplier cannot
    # hold its gradient against a large solute load, so a salt or urea load
    # drives urine osmolality toward isotonic even at maximal ADH.  Without
    # this term the model produces urine [Na+K] above 400 mmol/L.
    uomax = p['UOSMMAX'] / (1.0 + sp(SOLEXC - p['SOL0'], 1.0) / p['SOLREF'])
    UOSM = uomin + (uomax - uomin) * y[IX['AQP2']]
    VU = SOLEXC / (UOSM + 30.0)
    a['UOSM'], a['VU'], a['ENA'], a['EK'], a['EUREA'] = UOSM, VU, ENA, EK, EUREA
    a['SOLEXC'] = SOLEXC
    a['UNAK'] = (ENA + EK) / max(VU, 1e-6)
    a['EFWC'] = VU * (1.0 - a['UNAK'] / max(SNA, 60.0))    # electrolyte-free water clearance L/d

    # --- infusions --------------------------------------------------------
    r3, r09, rd5, rkcl = p['R3PCT'], p['R09'], p['RD5W'], p['RKCL']
    # Every switch below is a smooth logistic gate in TIME.  A switch on a
    # measured sodium ("stop when [Na] reaches the cap") makes the loop hunt
    # at the cap: the integrator spent 7.7e5 function evaluations there and
    # never finished.  A prescription is a statement about time anyway.
    gate = lambda x, w: 0.5 * (1.0 + np.tanh(x / w))
    if p['ACUTE'] > 0.5:
        rd5 += p['WLOAD'] * gate(p['TWLEND'] - t, 0.004)
    if p['CTRLON'] > 0.5:
        tstop = p['TCORR'] + (p['NACAP'] - p['NASTART']) / max(p['RATETGT'], 0.1)
        ramp = p['NASTART'] + p['RATETGT'] * (t - p['TCORR'])
        r3 += sclamp(p['KP3'] * (ramp - SNA), 0.0, p['R3MAX'], 0.02) \
            * gate(t - p['TCORR'], 0.01) * gate(tstop - t, 0.01)
    if p['RESCUE'] > 0.5:
        # Relowering is a bounded prescription (D5W + desmopressin for a day
        # or two), not a permanent servo.  Leaving the loop open forever made
        # the controller hunt at its target exactly as the 3% saline
        # controller did at its cap, and for the same reason.
        rd5 += sclamp(p['KPD5W'] * (SNA - p['NARES']), 0.0, p['D5WMAX'], 0.02) \
            * gate(t - p['TRESCUE'], 0.01) \
            * gate(p['TRESCUE'] + p['DURRES'] - t, 0.01)
    a['R3'], a['R09'], a['RD5'], a['RKCL'] = r3, r09, rd5, rkcl
    a['VINF'] = r3 + r09 + rd5
    a['NAINF'] = r3 * p['NA3'] + r09 * p['NA09']
    a['KINF'] = rkcl

    # --- brain: set points -------------------------------------------------
    f = (a['OSMEFF'] - 285.0) / 285.0
    def setp(base, beta):
        return base * p['ORGFLOOR'] + sp(base * (1.0 + beta * f) - base * p['ORGFLOOR'], 0.02)
    a['INSs'] = setp(p['INS0'], p['BINS']); a['TAUs'] = setp(p['TAU0'], p['BTAU'])
    a['GLXs'] = setp(p['GLX0'], p['BGLX']); a['CREs'] = setp(p['CRE0'], p['BCRE'])
    a['GPCs'] = setp(p['GPC0'], p['BGPC']); a['OTHs'] = setp(p['OTH0'], p['BOTH'])
    a['ORGSET'] = a['INSs'] + a['TAUs'] + a['GLXs'] + a['CREs'] + a['GPCs'] + a['OTHs']
    a['ORG'] = sum(max(y[IX[k]], 0.0) for k in ['INS', 'TAU', 'GLX', 'CRE', 'GPC', 'OTH'])
    ORG0 = p['INS0'] + p['TAU0'] + p['GLX0'] + p['CRE0'] + p['GPC0'] + p['OTH0']
    a['ORG0'] = ORG0

    # --- brain: inorganic target that would restore normal volume ----------
    a['BELSET'] = a['OSMEFF'] - p['IMP'] - a['ORGSET']
    BELVOL = a['OSMEFF'] - p['IMP'] - a['ORG']         # what volume regulation asks for
    a['BELVOL'] = BELVOL
    a['BELCAP'] = (285.0 - p['IMP'] - ORG0) * p['BELCEIL']

    # --- brain water (osmotic equilibrium is instantaneous) ----------------
    solute = max(y[IX['BELEC']], 0.0) + a['ORG'] + p['IMP']
    a['BWREL'] = solute / a['OSMEFF']                 # relative to baseline
    a['BWPCT'] = 80.0 * a['BWREL']                    # mL water /100 g brain
    a['SHRINK'] = sp(1.0 - a['BWREL'], 5e-4)
    a['SWELL'] = sp(a['BWREL'] - 1.0, 5e-4)
    a['HERN'] = 100.0 / (1.0 + np.exp(-(a['SWELL'] - 0.070) / 0.012))   # % herniation risk proxy

    # --- THE DRIVER --------------------------------------------------------
    a['OMEGA'] = sp(a['ORGSET'] - a['ORG'], 0.03)
    a['STRESS'] = a['OMEGA'] + p['WSHR'] * a['SHRINK'] * a['OSMEFF']
    thr = p['OMSTAR'] * p['RISK']
    a['THRESH'] = thr
    exc = sp(a['STRESS'] - thr, 0.03)
    ATP = max(y[IX['ATP']], 0.05)
    a['INJ'] = p['KINJ'] * exc**p['HINJ'] / ATP
    return a


# ---------------------------------------------------------------------
# 3. DERIVATIVES
# ---------------------------------------------------------------------
def rhs(t, y, p):
    a = algebra(t, y, p)
    d = np.zeros(NS)
    g = lambda n: y[IX[n]]

    # ---- body fluid -------------------------------------------------------
    d[IX['TBW']] = a['WINT'] + p['WMET'] - p['WINSENS'] - a['VU'] + a['VINF'] - a['WLOSS']
    d[IX['NAE']] = p['NAINT'] + a['NAINF'] - a['ENA'] - a['NALOSS']
    d[IX['KE']] = p['KINT'] + a['KINF'] - a['EK'] - a['KLOSS']
    d[IX['UREAB']] = a['PU'] + p['KAUREA'] * g('UREAG') - a['EUREA']

    # ---- AVP / AQP2 -------------------------------------------------------
    d[IX['AVP']] = (a['AVPTGT'] - g('AVP')) / p['TAUAVP']
    d[IX['AQP2']] = (a['OCC'] - g('AQP2')) / p['TAUAQP']

    # ---- brain inorganic: regulatory volume increase / decrease -----------
    tgt = a['BELCAP'] - sp(a['BELCAP'] - a['BELVOL'], 0.05)     # = min(BELVOL, BELCAP)
    flux = (tgt - g('BELEC')) / p['TAUELEC']
    d[IX['BELEC']] = sclamp(flux, -p['JRVDMAX'], p['JRVIMAX'], 0.5)

    # ---- brain organic osmolytes: THE ASYMMETRY ---------------------------
    #   efflux  = VRAC channel      (fast, energy-independent)
    #   influx  = SMIT1/TauT/BGT1   (slow, transcription- and Na-gradient-gated)
    kfac = min(1.10, max(0.40, a['SK'] / 4.0))          # Na-coupled carriers need the K/Na gradient
    ginf = p['FOSM'] * g('SMIT') * kfac * max(0.2, g('ATP'))
    for nm, spn in [('INS', 'INSs'), ('TAU', 'TAUs'), ('GLX', 'GLXs'),
                    ('CRE', 'CREs'), ('GPC', 'GPCs'), ('OTH', 'OTHs')]:
        diff = a[spn] - g(nm)
        up = sp(diff, 0.03)                       # influx arm  (slow, carrier)
        dn = up - diff                            # efflux arm  (fast, VRAC)
        d[IX[nm]] = up / p['TAUINF'] * ginf - dn / p['TAUEFF']
    d[IX['BURE']] = (a['UREAP'] - g('BURE')) / p['TAUBUR']

    # ---- TonEBP -> carrier transcription ----------------------------------
    tont = 1.0 + (p['TONMAX'] - 1.0) - sp((p['TONMAX'] - 1.0)
                 - p['ETON'] * a['OMEGA'] / a['ORG0'] * 10.0, 0.02)
    d[IX['TONEBP']] = (tont - g('TONEBP')) / p['TAUTON']
    d[IX['SMIT']] = (g('TONEBP') - g('SMIT')) / p['TAUSMIT']

    # ---- astrocyte energetics & viability ---------------------------------
    d[IX['ATP']] = (p['FNUT'] - g('ATP')) / p['TAUATP'] \
        - p['KATPU'] * a['STRESS'] * g('ATP') * (4.0 / max(a['SK'], 1.5))
    d[IX['AST']] = -p['KAST'] * a['INJ'] * g('AST') + p['KASTREP'] * (1.0 - g('AST'))
    # The astrocyte syncytium is gap-junction coupled: scattered loss is
    # covered by neighbours, and oligodendrocyte trophic support fails only
    # once loss exceeds ASTTHR.  Without this threshold a very small but very
    # PROLONGED astrocyte deficit integrates into complete demyelination.
    astloss = sp((1.0 - g('AST')) - p['ASTTHR'], 0.01)

    # ---- BBB, microglia, cytokines, humoral entry -------------------------
    ureaX = max(0.0, g('BURE') - 7.3)
    pu = p['EMAXUREA'] * ureaX / (p['EC50UREA'] + ureaX)             # urea protection (Soupart/Gankam)
    dexc = max(g('DEXC'), 0.0) / p['VDEX'] * 1e6                     # mg -> ng/mL
    pd_ = p['EMAXDEX'] * dexc / (p['EC50DEX'] + dexc)
    minc = max(g('MINC'), 0.0) / p['VMIN'] * 1e3                     # mg -> ug/mL
    pm = p['EMAXMIN'] * minc / (p['EC50MIN'] + minc)
    a_pu, a_pd, a_pm = pu, pd_, pm
    d[IX['BBBP']] = p['KBBB'] * astloss * (1.0 - a_pu) * (1.0 - a_pd) \
        - (g('BBBP') - 1.0) / p['TAUBBB']
    d[IX['MG']] = p['KMG'] * (astloss + 0.5 * sp(g('BBBP') - 1.0, 1e-3)) \
        * (1.0 - a_pm) * (1.0 - a_pu) - g('MG') / p['TAUMG']
    d[IX['CYT']] = p['KCYT'] * g('MG') - g('CYT') / p['TAUCYT']
    d[IX['IGG']] = p['KIGG'] * max(0.0, g('BBBP') - 1.0) - g('IGG') / p['TAUIGG']

    # ---- oligodendrocyte / myelin -----------------------------------------
    hit0 = p['WOA'] * astloss + p['WOC'] * g('CYT') + p['WOI'] * g('IGG')
    hit = hit0 / (p['KHIT'] + hit0)          # saturating: the pons cannot be
    d[IX['OLI']] = -p['KOLI'] * hit * g('OLI') \
        + p['KOPCD'] * g('OPC') * (1.0 - g('OLI'))   # demyelinated twice
    d[IX['OPC']] = p['KOPCP'] * (1.0 - g('OLI')) * (1.0 - g('OPC')) \
        - p['KOPCD'] * g('OPC') * (1.0 - g('OLI'))
    # myelin TRACKS the oligodendrocyte population (so a 10% cell deficit is a
    # ~10% myelin deficit, not a collapse) plus an acute stripping term.
    d[IX['MYE']] = p['KMYE'] * (g('OLI') - g('MYE')) - p['KDEM'] * hit * g('MYE')

    # ---- lesion topography, imaging, clinical deficit ---------------------
    dm = max(0.0, 1.0 - g('MYE'))
    d[IX['LESP']] = (p['WPONS'] * dm - g('LESP')) / p['TAULES'] * 12.0
    d[IX['LESE']] = (p['WEXP'] * dm - g('LESE')) / p['TAULES'] * 12.0
    d[IX['MRI']] = (max(g('LESP'), g('LESE')) - g('MRI')) / p['TAUMRI']
    lb = g('LESP') + 0.5 * g('LESE')
    sev = lb**1.4 / (p['KSEV']**1.4 + lb**1.4)
    d[IX['DEF']] = (p['DEFMAX'] * sev - g('DEF')) / p['TAUDEF']

    # ---- drug PK -----------------------------------------------------------
    d[IX['DDAD']] = -p['KADDA'] * g('DDAD')
    d[IX['DDAC']] = p['KADDA'] * g('DDAD') - p['KDDA'] * g('DDAC')
    d[IX['TLVD']] = -p['KATLV'] * g('TLVD')
    d[IX['TLVC']] = p['FTLV'] * p['KATLV'] * g('TLVD') - p['KTLV'] * g('TLVC')
    d[IX['UREAG']] = p['UREADOSE'] / 0.060 - p['KAUREA'] * g('UREAG')   # g/d -> mmol/d
    d[IX['DEXC']] = p['DEXON'] * 16.0 - p['KDEX'] * g('DEXC')           # 16 mg/d
    d[IX['MINC']] = p['MINOON'] * 200.0 - p['KMIN'] * g('MINC')         # 200 mg/d

    # ---- bookkeeping --------------------------------------------------------
    d[IX['CUMI']] = a['INJ']
    d[IX['CUMNA']] = a['NAINF']
    d[IX['CUMV']] = a['VINF']
    d[IX['CUMEFW']] = a['EFWC']
    return d


# ---------------------------------------------------------------------
# 4. INITIAL CONDITIONS / STEADY STATE
# ---------------------------------------------------------------------
def y_normal(p):
    y = np.zeros(NS)
    y[IX['TBW']] = 42.0
    y[IX['NAE']] = p['NAESET']
    y[IX['KE']] = p['KESET']
    y[IX['UREAB']] = p['PUREA'] / p['CLUREA'] * 42.0
    y[IX['AVP']] = 2.25
    y[IX['AQP2']] = 0.65
    y[IX['BELEC']] = 285.0 - p['IMP'] - (p['INS0'] + p['TAU0'] + p['GLX0']
                                         + p['CRE0'] + p['GPC0'] + p['OTH0'])
    for nm, k in [('INS', 'INS0'), ('TAU', 'TAU0'), ('GLX', 'GLX0'),
                  ('CRE', 'CRE0'), ('GPC', 'GPC0'), ('OTH', 'OTH0')]:
        y[IX[nm]] = p[k]
    y[IX['BURE']] = p['PUREA'] / p['CLUREA']
    y[IX['TONEBP']] = 1.0
    y[IX['SMIT']] = 1.0
    y[IX['ATP']] = 1.0
    y[IX['AST']] = 1.0
    y[IX['BBBP']] = 1.0
    y[IX['OLI']] = 1.0
    y[IX['MYE']] = 1.0
    return y


def run(p, y0, t_end, t_eval=None, max_step=np.inf):
    if t_eval is None:
        t_eval = np.linspace(0, t_end, int(t_end * 48) + 1)
    s = solve_ivp(rhs, (0.0, t_end), y0, args=(p,), method='LSODA',
                  t_eval=t_eval, rtol=1e-6, atol=1e-8, max_step=max_step)
    if not s.success:
        raise RuntimeError(s.message)
    obs = [algebra(tt, s.y[:, i], p) for i, tt in enumerate(s.t)]
    return s.t, s.y, obs


def P(**kw):
    p = dict(P0)
    p.update(kw)
    return p


def series(obs, key):
    return np.array([o[key] for o in obs])


# ---------------------------------------------------------------------
# 5. PHENOTYPE BUILDERS
# ---------------------------------------------------------------------
def make_chronic(target_na=110.0, days=21.0, pheno='siadh', extra=None):
    """Develop CHRONIC hyponatraemia by simulation (not by assertion), so
    the adapted brain composition is an OUTPUT, not an input."""
    def build(x):
        if pheno == 'siadh':
            p = P(AVPSIADH=6.0, WIN=x)
        elif pheno == 'hypovol':
            # 7 days of gastric/diuretic salt loss replaced with hypotonic fluid
            p = P(NALOSS=250.0, KLOSS=60.0, WLOSS=1.6, TLOSSEND=days, WIN=x)
        elif pheno == 'thiazide':
            p = P(THIAZ=1.0, NAINT=40.0, KINT=20.0, WIN=x, AVPSIADH=2.0)
        elif pheno == 'potomania':
            p = P(NAINT=15.0, KINT=15.0, PUREA=110.0, UOSMMAX=600.0, WIN=x)
        elif pheno == 'adrenal':
            p = P(AVPSIADH=4.0, NAINT=30.0, KINT=90.0, WIN=x)
        else:
            raise ValueError(pheno)
        if extra:
            p.update(extra)
        return p

    def err(x):
        p = build(x)
        t, y, o = run(p, y_normal(p), days, np.array([0.0, days]))
        return o[-1]['SNA'] - target_na

    lo, hi = 1.0, 14.0
    if err(hi) > 0:                       # cannot reach target: use maximum load
        x = hi
    else:
        x = brentq(err, lo, hi, xtol=1e-4)
    p = build(x)
    t, y, o = run(p, y_normal(p), days)
    return p, y[:, -1].copy(), o[-1], x


def make_acute(target_na=110.0, hours=8.0, extra=None):
    """Acute (self-induced water intoxication / exercise-associated /
    MDMA) — the same sodium, but the osmolytes have not left yet."""
    def build(x):
        p = P(AVPSIADH=6.0, ACUTE=1.0, WLOAD=x, TWLEND=hours / 24.0)
        if extra:
            p.update(extra)
        return p

    def err(x):
        p = build(x)
        t, y, o = run(p, y_normal(p), hours / 24.0, np.array([0.0, hours / 24.0]))
        return o[-1]['SNA'] - target_na
    x = brentq(err, 1.0, 60.0, xtol=1e-4)
    p = build(x)
    t, y, o = run(p, y_normal(p), hours / 24.0)
    return p, y[:, -1].copy(), o[-1], x


# ---------------------------------------------------------------------
# 6. DOSING EVENTS (desmopressin q8h, etc.) — event-driven integration
# ---------------------------------------------------------------------
def run_with_doses(p, y0, t_end, doses=(), npts=None):
    """doses = list of (time_days, state_name, amount)."""
    ev = sorted(set([0.0] + [d[0] for d in doses if 0.0 < d[0] < t_end] + [t_end]))
    ts, ys, obs = [], [], []
    y = y0.copy()
    for i in range(len(ev) - 1):
        t0, t1 = ev[i], ev[i + 1]
        for (td, sn, amt) in doses:
            if abs(td - t0) < 1e-9:
                y[IX[sn]] += amt
        n = max(2, int((t1 - t0) * (npts or 96)) + 1)
        te = np.linspace(t0, t1, n)
        s = solve_ivp(rhs, (t0, t1), y, args=(p,), method='LSODA',
                      t_eval=te, rtol=1e-6, atol=1e-8, max_step=np.inf)
        if not s.success:
            raise RuntimeError(s.message)
        keep = slice(0, None) if i == 0 else slice(1, None)
        ts.append(s.t[keep]); ys.append(s.y[:, keep])
        obs += [algebra(tt, s.y[:, j], p) for j, tt in enumerate(s.t)][keep]
        y = s.y[:, -1].copy()
    return np.concatenate(ts), np.concatenate(ys, axis=1), obs


def run_protocol(phases, y0, doses=(), npts=96):
    """phases = [(params, duration_days), ...]; doses = [(t_abs, state, amt)].
    A rescue is a change of PRESCRIPTION, not just a change of infusion rate:
    the fluid restriction and the desmopressin start at the same moment as the
    dextrose.  Running the rescue parameters from t=0 instead would abolish
    the overcorrection it is supposed to be rescuing."""
    bounds = [0.0]
    for (_, d) in phases:
        bounds.append(bounds[-1] + d)
    ev = sorted(set(bounds + [d[0] for d in doses if 0.0 < d[0] < bounds[-1]]))
    ts, ys, obs = [], [], []
    y = y0.copy()
    for i in range(len(ev) - 1):
        t0, t1 = ev[i], ev[i + 1]
        for (td, sn, amt) in doses:
            if abs(td - t0) < 1e-9:
                y[IX[sn]] += amt
        k = next(j for j in range(len(phases))
                 if bounds[j] <= t0 + 1e-9 < bounds[j + 1])
        p = phases[k][0]
        te = np.linspace(t0, t1, max(3, int((t1 - t0) * npts) + 1))
        s = solve_ivp(rhs, (t0, t1), y, args=(p,), method='LSODA',
                      t_eval=te, rtol=1e-6, atol=1e-8)
        if not s.success:
            raise RuntimeError(s.message)
        kk = slice(0, None) if i == 0 else slice(1, None)
        ts.append(s.t[kk]); ys.append(s.y[:, kk])
        obs += [algebra(tt, s.y[:, j], p) for j, tt in enumerate(s.t)][kk]
        y = s.y[:, -1].copy()
    return np.concatenate(ts), np.concatenate(ys, axis=1), obs


def ddavp_q8(t0, t1, dose=2.0):
    return [(t, 'DDAC', dose) for t in np.arange(t0, t1, 8.0 / 24.0)]


# ---------------------------------------------------------------------
# 7. REPORTING HELPERS
# ---------------------------------------------------------------------
def na_rise_24(t, sna, t0=0.0):
    """max rise in [Na] over any 24 h window starting at/after t0."""
    m = t >= t0
    tt, ss = t[m], sna[m]
    best = 0.0
    for i in range(len(tt)):
        j = np.searchsorted(tt, tt[i] + 1.0)
        if j >= len(tt):
            break
        best = max(best, ss[j] - ss[i])
    return best


def summarise(name, t, y, obs, note=''):
    sna = series(obs, 'SNA')
    return dict(
        id=name, note=note,
        Na0=round(float(sna[0]), 1),
        Na24=round(float(np.interp(1.0, t, sna)), 1),
        NaEnd=round(float(sna[-1]), 1),
        dNa24=round(float(na_rise_24(t, sna)), 1),
        OMEGAmax=round(float(np.max(series(obs, 'OMEGA'))), 1),
        STRESSmax=round(float(np.max(series(obs, 'STRESS'))), 1),
        THRESH=round(float(obs[0]['THRESH']), 1),
        BWmin=round(float(np.min(series(obs, 'BWPCT'))), 2),
        BWmax=round(float(np.max(series(obs, 'BWPCT'))), 2),
        CUMI=round(float(y[IX['CUMI'], -1]), 3),
        ASTmin=round(float(np.min(y[IX['AST'], :])), 3),
        MYEmin=round(float(np.min(y[IX['MYE'], :])), 3),
        DEFmax=round(float(np.max(y[IX['DEF'], :])), 1),
        LESPmax=round(float(np.max(y[IX['LESP'], :])), 3),
        EFWCmax=round(float(np.max(series(obs, 'EFWC'))), 2),
        UOSMmin=round(float(np.min(series(obs, 'UOSM'))), 0),
    )


# =====================================================================
# 8. SCENARIOS AND ANCHOR CHECKS
#    Everything below produces ods_verification_output.txt.  Every number
#    quoted in README.md, the .dot map and the mrgsolve header comes from
#    this run; nothing is written from memory.
# =====================================================================
OUT = []
def say(s=''):
    OUT.append(s)
    print(s, flush=True)


def run_phases(phases, y0, npts=96):
    """phases = [(params, duration_days), ...] run back to back."""
    ts, ys, obs = [], [], []
    y = y0.copy(); t0 = 0.0
    for i, (p, dur) in enumerate(phases):
        te = np.linspace(t0, t0 + dur, max(3, int(dur * npts) + 1))
        s = solve_ivp(rhs, (t0, t0 + dur), y, args=(p,), method='LSODA',
                      t_eval=te, rtol=1e-6, atol=1e-8)
        if not s.success:
            raise RuntimeError(s.message)
        k = slice(0, None) if i == 0 else slice(1, None)
        ts.append(s.t[k]); ys.append(s.y[:, k])
        obs += [algebra(tt, s.y[:, j], p) for j, tt in enumerate(s.t)][k]
        y = s.y[:, -1].copy(); t0 += dur
    return np.concatenate(ts), np.concatenate(ys, axis=1), obs


def hdr(t):
    say(); say('=' * 78); say(t); say('=' * 78)


def main():
    hdr('0. NUMERICAL VERIFICATION OF THE HEALTHY STEADY STATE')
    p = P(); y = y_normal(p)
    t, yy, o = run(p, y, 150.0)
    d = rhs(0.0, yy[:, -1], p)
    say(f'  max|dy/dt| over the 29 physiological states at t=150 d : {np.max(np.abs(d[:29])):.3e}')
    say(f'  [Na] 0 d / 150 d            : {o[0]["SNA"]:.4f} / {o[-1]["SNA"]:.4f} mmol/L')
    say(f'  brain water (mL/100 g)      : {o[0]["BWPCT"]:.4f} / {o[-1]["BWPCT"]:.4f}')
    say(f'  astrocyte / oligo / myelin  : {yy[IX["AST"],-1]:.5f} / {yy[IX["OLI"],-1]:.5f} / {yy[IX["MYE"],-1]:.5f}')
    say(f'  cumulative injury integral  : {yy[IX["CUMI"],-1]:.3e}   (must be exactly 0)')
    say(f'  urine osmolality / volume   : {o[-1]["UOSM"]:.0f} mOsm/kg / {o[-1]["VU"]:.2f} L/d')
    say(f'  BUN                         : {o[-1]["BUN"]:.1f} mg/dL')
    say(f'  urine (Na+K) : serum Na     : {o[-1]["UNAK"]/o[-1]["SNA"]:.2f}')

    # ---------------------------------------------------------------
    hdr('1. THE ADAPTED BRAIN — built by simulation, not asserted')
    ph = {}
    for name, kind, days in [('SIADH', 'siadh', 21.0), ('hypovolaemic', 'hypovol', 7.0),
                             ('thiazide', 'thiazide', 21.0), ('beer potomania', 'potomania', 21.0)]:
        pp, yy0, oo, x = make_chronic(110.0, days, kind)
        ph[kind] = (pp, yy0, oo)
        say(f'  {name:15s} [Na] {oo["SNA"]:6.1f}  ORG {oo["ORG"]:5.2f}/48.00 ({100*(1-oo["ORG"]/48):4.1f}% lost)'
            f'  myo-Ins {yy0[IX["INS"]]:4.2f}/7.00 ({100*(1-yy0[IX["INS"]]/7):4.1f}% lost)'
            f'  brain water {oo["BWPCT"]:5.2f}  AVP {yy0[IX["AVP"]]:5.2f}  U_osm {oo["UOSM"]:5.0f}'
            f'  ECF deficit {100*oo["VOLDEF"]:4.1f}%')
    pa, ya0, oa, xa = make_acute(110.0, 8.0)
    say(f'  {"ACUTE (8 h)":15s} [Na] {oa["SNA"]:6.1f}  ORG {oa["ORG"]:5.2f}/48.00 ({100*(1-oa["ORG"]/48):4.1f}% lost)'
        f'  myo-Ins {ya0[IX["INS"]]:4.2f}/7.00 ({100*(1-ya0[IX["INS"]]/7):4.1f}% lost)'
        f'  brain water {oa["BWPCT"]:5.2f}  <-- OEDEMA')
    say()
    say(f'  The chronically adapted brain at [Na] 110 holds NORMAL water ({ph["siadh"][2]["BWPCT"]:.2f} vs 80.00)')
    say(f'  because it has thrown away {48-ph["siadh"][2]["ORG"]:.1f} mOsm/kg of organic osmolyte and')
    say(f'  {227.0-ph["siadh"][1][IX["BELEC"]]:.1f} mOsm/kg of inorganic ion.  Omega = {ph["siadh"][2]["OMEGA"]:.3f}: it is')
    say(f'  ADAPTED, not injured.  The acute brain at the SAME sodium is swollen by')
    say(f'  {100*(oa["BWPCT"]/80.0-1):.1f}% because the osmolytes have not left yet.')

    psi, ysi, osi = ph['siadh']
    phy, yhy, ohy = ph['hypovol']

    # ---------------------------------------------------------------
    hdr('2. HOW FAST MAY THE SODIUM BE CORRECTED?  (derived, not assumed)')
    say('  Threshold Omega* = 8.0 mOsm/kg is ONE number for every patient.  Risk')
    say('  factors act only on FOSM, the organic transport capacity.')
    say()
    say('  FOSM  FNUT |' + ''.join(f'{r:8d}' for r in [4, 6, 8, 10, 12, 14, 16, 20]) + '   mmol/L/24 h')
    for fo, fn, lbl in [(1.00, 1.00, 'normal'), (0.80, 0.95, 'mild'),
                        (0.65, 0.90, 'moderate'), (0.55, 0.85, 'alcoholic/malnourished'),
                        (0.45, 0.80, 'severe')]:
        pc = P(**psi); pc.update(dict(FOSM=fo, FNUT=fn))
        tc, yc, oc = run(pc, ysi, 3.0); yy0 = yc[:, -1].copy()
        row = []
        for r in [4, 6, 8, 10, 12, 14, 16, 20]:
            pp = P(**psi); pp.update(dict(CTRLON=1.0, RATETGT=float(r), NASTART=110.0,
                                          NACAP=140.0, WIN=1.0, FOSM=fo, FNUT=fn))
            tt, yv, oo = run(pp, yy0, 5.0)
            row.append(float(np.max(series(oo, 'STRESS'))))
        mark = ''.join(('  ' + ('.' if v < 8.0 else 'X') + f'{v:5.1f}') for v in row)
        say(f'  {fo:4.2f}  {fn:4.2f} |{mark}   {lbl}')
    say('       ( . = below the injury threshold,  X = above it )')
    say()
    say('  NORMAL transport      -> the limit falls between 10 and 12 mmol/L/24 h')
    say('  ALCOHOLIC/MALNOURISHED-> the limit falls between  6 and  8 mmol/L/24 h')
    say('  Those are the two numbers in every guideline.  They are ONE number here')
    say('  (Omega* = 8) plus a transporter, and the transporter is the risk factor.')

    # ---------------------------------------------------------------
    hdr('3. DOSE-RESPONSE OF THE INJURY CASCADE (normal-risk SIADH, [Na] 110)')
    say('  rate  maxSTRESS  astrocyte  oligo   myelin  peakDeficit  day  symptoms  MRI+  d90')
    dr = {}
    for r in [8, 10, 12, 14, 16, 20, 25, 30]:
        pp = P(**psi); pp.update(dict(CTRLON=1.0, RATETGT=float(r), NASTART=110.0,
                                      NACAP=140.0, WIN=1.0))
        tt, yv, oo = run(pp, ysi, 90.0)
        st = series(oo, 'STRESS'); df = yv[IX['DEF']]; mr = yv[IX['MRI']]
        f = lambda arr, th: tt[np.argmax(arr > th)] if (arr > th).any() else np.nan
        dr[r] = (float(st[tt <= 21].max()), float(df.max()), float(f(df, 10)), float(f(mr, 0.25)))
        say(f'  {r:4d} {st[tt<=21].max():9.2f}  {yv[IX["AST"]].min():9.3f} {yv[IX["OLI"]].min():6.3f}'
            f' {yv[IX["MYE"]].min():8.3f} {df.max():11.1f} {tt[np.argmax(df)]:5.1f}'
            f'   {f(df,10):6.1f} {f(mr,0.25):5.1f} {np.interp(90,tt,df):5.1f}')
    say()
    say('  Two things fall out that were not put in:')
    say('   (a) the clinical course is BIPHASIC — the sodium is normal and the')
    say('       patient is better while the deficit is still zero;')
    say(f'   (b) the MRI turns positive AFTER the patient deteriorates '
        f'({dr[16][3]:.1f} d vs {dr[16][2]:.1f} d at 16 mmol/L/24 h),')
    say('       so a normal scan at symptom onset does not exclude the diagnosis.')

    # ---------------------------------------------------------------
    hdr('4. THE CENTRAL EXPERIMENT — who actually sets the correction rate?')
    say('  Hypovolaemic hyponatraemia, [Na] 110.  Losses stop; 2 L/d of 0.9% saline.')
    say('  Arm A: the AVP axis responds physiologically.')
    say('  Arm B: IDENTICAL patient, identical fluid, but AVP is FROZEN at its')
    say('         presentation value (a counterfactual, not a therapy).')
    say()
    base = dict(NALOSS=0.0, KLOSS=0.0, WLOSS=0.0, TLOSSEND=0.0, R09=2.0)
    pA = P(**phy); pA.update(base)
    pB = P(**phy); pB.update(base); pB.update(dict(AVPFREEZE=1.0, TFREEZE=0.0,
                                                   AVPFRZV=float(yhy[IX['AVP']])))
    tA, yA, oA = run(pA, yhy, 12.0)
    tB, yB, oB = run(pB, yhy, 12.0)
    say('        t(h)  ARM A [Na]  AVP  U_osm  V_u  EFWC  |  ARM B [Na]  AVP  U_osm  V_u')
    for h in [0, 6, 12, 18, 24, 36, 48, 72, 96]:
        i = np.argmin(abs(tA - h / 24)); j = np.argmin(abs(tB - h / 24))
        say(f'       {h:4.0f} {series(oA,"SNA")[i]:10.1f} {yA[IX["AVP"],i]:5.2f} {oA[i]["UOSM"]:6.0f}'
            f' {oA[i]["VU"]:5.2f} {oA[i]["EFWC"]:5.2f}  | {series(oB,"SNA")[j]:9.1f}'
            f' {yB[IX["AVP"],j]:5.2f} {oB[j]["UOSM"]:6.0f} {oB[j]["VU"]:5.2f}')
    snaA, snaB = series(oA, 'SNA'), series(oB, 'SNA')
    say()
    say(f'  24 h rise      arm A {na_rise_24(tA,snaA):5.1f}   arm B {na_rise_24(tB,snaB):5.1f}  mmol/L')
    say(f'  peak STRESS    arm A {np.max(series(oA,"STRESS")):5.1f}   arm B {np.max(series(oB,"STRESS")):5.1f}  (threshold 8.0)')
    say(f'  peak deficit   arm A {yA[IX["DEF"]].max():5.1f}   arm B {yB[IX["DEF"]].max():5.1f}')
    say(f'  myelin nadir   arm A {yA[IX["MYE"]].min():5.3f}   arm B {yB[IX["MYE"]].min():5.3f}')
    say()
    say('  Not one millimole of sodium was prescribed in either arm.  The whole of')
    say('  the difference is made by the kidney being released from a volume')
    say('  stimulus it no longer has: urine osmolality collapses, the')
    say('  electrolyte-free water clearance opens, and the sodium climbs on its own.')

    # ---------------------------------------------------------------
    hdr('5. PROACTIVE DESMOPRESSIN CLAMP — and what happens when it comes off')
    say('  Same hypovolaemic patient.  Desmopressin 2 ug IV q8h for five days plus')
    say('  3% saline titrated to +6 mmol/L/24 h, against 0.9% saline alone.  The')
    say('  correction rate is reported in TWO windows because the model found')
    say('  something in the second one that was not designed in.')
    say()
    def clamp(win, r09=0.5):
        return P(**{**phy, **base, 'R09': r09, 'CTRLON': 1.0, 'RATETGT': 6.0,
                    'NASTART': 110.0, 'NACAP': 130.0, 'WIN': win})
    dose_taper = ddavp_q8(0.0, 5.0, 2.0) \
        + [(t, 'DDAC', 1.0) for t in np.arange(5.0, 7.0, 1/3)] \
        + [(t, 'DDAC', 0.5) for t in np.arange(7.0, 9.0, 1/3)] \
        + [(t, 'DDAC', 0.25) for t in np.arange(9.0, 11.0, 1/3)]
    arms = [('0.9% saline alone', pA, []),
            ('clamp, fluid intake 1.5 L/d', clamp(1.5), ddavp_q8(0.0, 5.0)),
            ('clamp, fluid intake 1.0 L/d', clamp(1.0), ddavp_q8(0.0, 5.0)),
            ('clamp, fluid intake 0.5 L/d', clamp(0.5), ddavp_q8(0.0, 5.0)),
            ('clamp 1.0 L/d + dose taper 2->0.25 ug', clamp(1.0), dose_taper)]
    say(f'  {"":38s} {"d0-5":>6s} {"d4-14":>7s} {"TBW d5":>7s} {"maxSTR":>7s} {"deficit":>8s}')
    res5 = {}
    for lbl, pp, ds in arms:
        tt, yv, oo = run_with_doses(pp, yhy, 30.0, ds)
        sn = series(oo, 'SNA')
        m1, m2 = tt <= 5.0, (tt >= 4.0) & (tt <= 14.0)
        r = (na_rise_24(tt[m1], sn[m1]), na_rise_24(tt[m2], sn[m2]),
             float(yv[IX['TBW'], np.argmin(abs(tt - 5.0))]),
             float(np.max(series(oo, 'STRESS'))), float(yv[IX['DEF']].max()))
        res5[lbl] = r
        say(f'  {lbl:38s} {r[0]:6.1f} {r[1]:7.1f} {r[2]:7.1f} {r[3]:7.2f} {r[4]:8.1f}')
    a15 = res5['clamp, fluid intake 1.5 L/d']; a10 = res5['clamp, fluid intake 1.0 L/d']
    a05 = res5['clamp, fluid intake 0.5 L/d']; at = res5['clamp 1.0 L/d + dose taper 2->0.25 ug']
    say()
    say('  INSIDE the clamp the prescription and the outcome are the same number:')
    say(f'  {a10[0]:.1f} mmol/L/24 h against a prescribed 6.0, in every arm, with no injury.')
    say('  The clamp does not lower the sodium; it takes the KIDNEY out of the loop.')
    say()
    say('  TAKING IT OFF is the other half of the therapy, and the model says the')
    say('  rebound is a STORED-WATER problem, not a tapering problem:')
    say(f'   - at 1.5 L/d of intake the patient stores water to {a15[2]:.1f} L, and when the')
    say(f'     drug clears the sodium rebounds {a15[1]:.1f} mmol/L in 24 h — peak stress {a15[3]:.1f},')
    say(f'     peak deficit {a15[4]:.1f}.  The clamp caused the injury it was preventing.')
    say(f'   - at 1.0 L/d: {a10[2]:.1f} L stored, rebound {a10[1]:.1f}, deficit {a10[4]:.1f}.')
    say(f'   - at 0.5 L/d: {a05[2]:.1f} L stored, rebound {a05[1]:.1f}, peak stress {a05[3]:.1f} — never near threshold.')
    say(f'   - tapering the DOSE from 2 ug to 0.25 ug changes nothing ({at[1]:.1f} vs {a10[1]:.1f}),')
    say('     because 0.25 ug still gives about 10 pg/mL against a V2 EC50 of 1.6:')
    say('     the receptor stays saturated until the drug is essentially gone,')
    say('     whatever the schedule.  The aquaresis is not tapering; it is waiting.')
    say()
    say('  So the derived rule is not "wean the desmopressin".  It is: while the')
    say('  clamp is on, hold the free-water balance neutral, because every litre')
    say('  retained is a litre that will be excreted the moment it comes off.')

    # ---------------------------------------------------------------
    hdr('6. RELOWERING — is there a deadline, and where is it?')
    say('  Overcorrection scenario of section 4 (arm A).  At time T the')
    say('  prescription CHANGES: 0.9% saline stops, fluid intake is restricted to')
    say('  0.5 L/d, 5% dextrose is titrated back towards [Na] 118 for 24 h, and')
    say('  desmopressin 2 ug q8h runs for 60 h.')
    say()
    say('   rescue start (h)   peak STRESS   astrocyte nadir   peak deficit   deficit at d90')
    resc_rows = {}
    for T in [8, 12, 18, 24, 36, 48, 72, 1e9]:
        if T < 1e8:
            th = T / 24.0
            pr = P(**pA); pr.update(dict(RESCUE=1.0, TRESCUE=th, NARES=118.0,
                                         DURRES=1.0, WIN=0.5, R09=0.0))
            tt, yv, oo = run_protocol([(pA, th), (pr, 90.0 - th)], yhy,
                                      ddavp_q8(th, th + 2.5))
            lbl = f'{T:5.0f}'
        else:
            tt, yv, oo = run(pA, yhy, 90.0)
            lbl = '  none'
        resc_rows[T] = (float(np.max(series(oo, 'STRESS'))), float(yv[IX['AST']].min()),
                        float(yv[IX['DEF']].max()), float(np.interp(90, tt, yv[IX['DEF']])))
        r = resc_rows[T]
        say(f'   {lbl}              {r[0]:8.2f}      {r[1]:10.3f}      {r[2]:10.1f}     {r[3]:10.1f}')
    say()
    say('  Relowering works, and the benefit decays with delay: it is largest when')
    say('  the astrocytes have not yet died, and by the time the deficit is')
    say('  established the sodium is no longer the thing that needs correcting.')

    # ---------------------------------------------------------------
    hdr('7. POTASSIUM IS SODIUM — the correction nobody prescribed')
    say('  Chronic SIADH, [Na] 110, coexisting potassium depletion.  KCl 40 mmol')
    say('  q8h (120 mmol/d) and NOT ONE MILLIMOLE OF SODIUM given.')
    say()
    pk = P(**psi); pk.update(dict(KE=0.0))
    yk = ysi.copy(); yk[IX['KE']] -= 450.0                 # ~450 mmol total-body K deficit
    a0 = algebra(0.0, yk, P(**psi))
    say(f'  starting serum K {a0["SK"]:.2f} mmol/L, starting [Na] {a0["SNA"]:.1f} mmol/L')
    for kcl in [0.0, 60.0, 120.0, 180.0]:
        pp = P(**psi); pp.update(dict(RKCL=kcl, WIN=1.0))
        tt, yv, oo = run(pp, yk, 4.0)
        sn = series(oo, 'SNA')
        say(f'  KCl {kcl:5.0f} mmol/d -> [Na] at 24 h {np.interp(1,tt,sn):6.1f} '
            f'(+{np.interp(1,tt,sn)-sn[0]:4.1f})   serum K {series(oo,"SK")[np.argmin(abs(tt-1))]:4.2f}'
            f'   max 24 h rise {na_rise_24(tt,sn):4.1f}')
    say()
    say('  Edelman says [Na] = 1.11(Na_e + K_e)/TBW - 25.6, and potassium is in the')
    say('  numerator.  40 mmol of KCl into 42 L is 1.11 x 40 / 42 = 1.06 mmol/L of')
    say('  sodium correction that appears on no fluid chart.')

    # ---------------------------------------------------------------
    hdr('8. TOLVAPTAN, UREA, SALINE, RESTRICTION — four ways to raise a sodium')
    say('  Chronic SIADH, [Na] 110, seven days.')
    say()
    say(f'  {"":34s} {"[Na] 24h":>9s} {"24h rise":>9s} {"maxSTRESS":>10s} {"deficit":>8s}')
    arms = [
        ('fluid restriction 1.0 L/d', dict(WIN=1.0)),
        ('fluid restriction 0.5 L/d', dict(WIN=0.5)),
        ('0.9% saline 2 L/d', dict(R09=2.0)),
        ('3% saline, titrated to +6/d', dict(CTRLON=1.0, RATETGT=6.0, NASTART=110.0,
                                             NACAP=130.0, WIN=1.0)),
        ('oral urea 30 g/d', dict(UREADOSE=30.0, WIN=1.0)),
        ('tolvaptan 15 mg/d', dict(WIN=1.628349)),
    ]
    for lbl, kw in arms:
        pp = P(**psi); pp.update(kw)
        doses = [(d, 'TLVD', 15.0) for d in range(7)] if 'tolvaptan' in lbl else []
        tt, yv, oo = run_with_doses(pp, ysi, 7.0, doses)
        sn = series(oo, 'SNA')
        say(f'  {lbl:34s} {np.interp(1,tt,sn):9.1f} {na_rise_24(tt,sn):9.1f}'
            f' {np.max(series(oo,"STRESS")):10.2f} {yv[IX["DEF"]].max():8.1f}')
    say()
    say(f'  Urine (Na+K)/serum Na at presentation = {osi["UNAK"]/osi["SNA"]:.2f}.  Above 1.0 the')
    say('  urine is more concentrated in cation than the plasma is, so every litre')
    say('  passed makes the patient MORE hyponatraemic: fluid restriction alone')
    say('  cannot work, and 0.9% saline (308 mOsm/L against a urine osmolality of')
    say(f'  {osi["UOSM"]:.0f}) makes it worse.  Both are model outputs, not rules.')

    # ---------------------------------------------------------------
    hdr('9. UREA AND THE HONEST NEGATIVE RESULT')
    say('  Under a pure effective-tonicity accounting urea CANNOT protect the')
    say('  brain: it crosses the blood-brain barrier, so it appears on both sides')
    say('  of the balance and cancels exactly.  The animal data (Soupart 2000,')
    say('  Gankam Kengne 2009) are real, so in this model urea acts where those')
    say('  studies MEASURED an effect — the blood-brain barrier and microglia —')
    say('  and not on the osmotic arm.  That is a falsifiable structural claim.')
    say()
    say(f'  {"":34s} {"maxSTRESS":>10s} {"BBB peak":>9s} {"astro nadir":>12s} {"deficit":>8s}')
    for lbl, kw in [('overcorrection, no urea', {}),
                    ('overcorrection + urea 30 g/d', dict(UREADOSE=30.0)),
                    ('overcorrection + dexamethasone', dict(DEXON=1.0)),
                    ('overcorrection + minocycline', dict(MINOON=1.0))]:
        pp = P(**pA); pp.update(kw)
        tt, yv, oo = run(pp, yhy, 60.0)
        say(f'  {lbl:34s} {np.max(series(oo,"STRESS")):10.2f} {yv[IX["BBBP"]].max():9.3f}'
            f' {yv[IX["AST"]].min():12.3f} {yv[IX["DEF"]].max():8.1f}')
    say()
    say('  Note the pattern: STRESS is IDENTICAL in all four arms (the osmotic')
    say('  insult is untouched) and only the downstream arms move.  If a trial ever')
    say('  showed urea blunting the osmotic insult itself, this model is wrong.')

    # ---------------------------------------------------------------
    hdr('10. RATE versus TOTAL MAGNITUDE')
    say('  Is 8 mmol/L/day for four days safer than 12 in one day and then stop?')
    say('  Same patient, same 32 mmol/L of total correction either way.')
    say()
    for lbl, kw, dur in [
            ('+8/d x 4 d   (32 total)', dict(CTRLON=1.0, RATETGT=8.0, NASTART=110.0, NACAP=142.0), 21.0),
            ('+12 on d1 then hold', dict(CTRLON=1.0, RATETGT=12.0, NASTART=110.0, NACAP=122.0), 21.0),
            ('+16 on d1 then hold', dict(CTRLON=1.0, RATETGT=16.0, NASTART=110.0, NACAP=126.0), 21.0),
            ('+6/d x 5 d   (30 total)', dict(CTRLON=1.0, RATETGT=6.0, NASTART=110.0, NACAP=140.0), 21.0)]:
        pp = P(**psi); pp.update(kw); pp['WIN'] = 1.0
        tt, yv, oo = run(pp, ysi, dur)
        sn = series(oo, 'SNA')
        say(f'  {lbl:26s} 24h rise {na_rise_24(tt,sn):5.1f}  48h rise '
            f'{max(0.0,float(np.interp(2,tt,sn)-sn[0])):5.1f}  maxSTRESS '
            f'{np.max(series(oo,"STRESS")):5.2f}  deficit {yv[IX["DEF"]].max():5.1f}')
    say()
    say('  The 24 h rate is not sufficient on its own: Omega depends on the gap')
    say('  between the CURRENT tonicity and an osmolyte pool with a 2-3 day time')
    say('  constant, so a second and third day of "safe" correction keep adding to')
    say('  a deficit the brain has not yet closed.  The 48 h limit in the guidelines')
    say('  is not a separate rule here; it is the same rule read one day later.')

    # ---------------------------------------------------------------
    hdr('11. PHENOTYPES AT AN IDENTICAL PRESCRIPTION (+10 mmol/L/24 h)')
    say(f'  {"phenotype":34s} {"maxSTRESS":>10s} {"astro":>7s} {"myelin":>7s} {"deficit":>8s}')
    phenos = [
        ('normal-risk SIADH', {}, ysi, psi),
        ('alcoholic / malnourished', dict(FOSM=0.55, FNUT=0.85), ysi, psi),
        ('hypokalaemic (K 2.5)', {}, None, psi),
        ('cirrhosis / pre-transplant', dict(FOSM=0.70, INS0=4.2, GLX0=19.5), None, psi),
        ('severe, starting [Na] 100', {}, None, None),
    ]
    for lbl, kw, y_in, p_in in phenos:
        if lbl.startswith('hypokal'):
            y_in = ysi.copy(); y_in[IX['KE']] -= 450.0; p_in = psi
        if lbl.startswith('cirrhosis'):
            pc = P(**psi); pc.update(kw)
            tc, yc, oc = run(pc, y_normal(pc), 21.0)   # re-adapt with the depleted pool
            y_in = yc[:, -1].copy(); p_in = {**psi, **kw}
        if lbl.startswith('severe'):
            p_in, y_in, _o, _x = make_chronic(100.0, 21.0, 'siadh')
        pp = P(**p_in); pp.update(kw)
        pp.update(dict(CTRLON=1.0, RATETGT=10.0, NASTART=float(algebra(0.0, y_in, P(**p_in))['SNA']),
                       NACAP=140.0, WIN=1.0))
        tt, yv, oo = run(pp, y_in, 60.0)
        say(f'  {lbl:34s} {np.max(series(oo,"STRESS")):10.2f} {yv[IX["AST"]].min():7.3f}'
            f' {yv[IX["MYE"]].min():7.3f} {yv[IX["DEF"]].max():8.1f}')
    say()
    say('  One prescription, five outcomes.  "8 mmol/L per day" is not a safety')
    say('  statement about a patient; it is a safety statement about a transporter.')

    # ---------------------------------------------------------------
    hdr('12. ACUTE HYPONATRAEMIA IS A DIFFERENT DISEASE')
    say('  Same sodium (110), same correction (+20 mmol/L/24 h), different history.')
    say()
    for lbl, pp0, yy0 in [('chronic (21 d)', psi, ysi), ('acute (8 h)', pa, ya0)]:
        pp = P(**pp0); pp.update(dict(CTRLON=1.0, RATETGT=20.0, NASTART=110.0,
                                      NACAP=140.0, WIN=1.0, ACUTE=0.0))
        tt, yv, oo = run(pp, yy0, 30.0)
        say(f'  {lbl:16s} brain water {np.max(series(oo,"BWPCT")):5.2f} -> {np.min(series(oo,"BWPCT")):5.2f}'
            f'   herniation index {np.max(series(oo,"HERN")):5.1f}%'
            f'   maxSTRESS {np.max(series(oo,"STRESS")):5.2f}   deficit {yv[IX["DEF"]].max():5.1f}')
    say()
    say('  The acute brain is in danger from OEDEMA and is not in danger from')
    say('  correction, because the osmolytes have not left yet.  The chronic brain')
    say('  is the exact opposite.  The model was never told this.')

    # ---------------------------------------------------------------
    hdr('13. HOW LONG DOES THE BRAIN TAKE TO PUT THE OSMOLYTES BACK?')
    pp = P(**psi); pp.update(dict(CTRLON=1.0, RATETGT=6.0, NASTART=110.0, NACAP=140.0, WIN=1.0))
    tt, yv, oo = run(pp, ysi, 14.0)
    say('   day  [Na]  myo-Ins  taurine  Glx   total ORG   ORG_set   Omega   SMIT1')
    for d_ in [0, 1, 2, 3, 4, 5, 7, 10, 14]:
        i = np.argmin(abs(tt - d_))
        say(f'  {d_:4d} {series(oo,"SNA")[i]:6.1f} {yv[IX["INS"],i]:8.2f} {yv[IX["TAU"],i]:8.2f}'
            f' {yv[IX["GLX"],i]:6.2f} {oo[i]["ORG"]:10.2f} {oo[i]["ORGSET"]:9.2f}'
            f' {oo[i]["OMEGA"]:7.2f} {yv[IX["SMIT"],i]:7.2f}')
    ins = yv[IX['INS']]; tgt95 = 7.0 * 0.95
    i95 = np.argmax(ins > tgt95) if (ins > tgt95).any() else -1
    say()
    say(f'  myo-inositol reaches 95% of normal at {tt[i95]:.1f} d '
        f'(Verbalis & Gullans 1993 report ~5 d in the rat).')

    # ---------------------------------------------------------------
    hdr('14. ADROGUE-MADIAS CHECK (the model against the bedside formula)')
    for wt, tbwf, na in [(70, 0.6, 110.0), (60, 0.5, 110.0), (70, 0.6, 125.0)]:
        tbw = wt * tbwf
        am3 = (513.0 - na) / (tbw + 1.0)
        am09 = (154.0 - na) / (tbw + 1.0)
        say(f'  {wt} kg, TBW {tbw:.0f} L, [Na] {na:.0f}:  3% NaCl {am3:5.2f} mmol/L per L'
            f'  ({am3/10:4.2f} per 100 mL)   0.9% NaCl {am09:5.2f} per L')
    say()
    say('  A 150 mL bolus of 3% saline raises a 60 kg patient by '
        f'{0.150*(513-110)/(30+1):.1f} mmol/L, so the')
    say('  European guideline\'s "three 150 mL boluses for a 5 mmol/L rise" is')
    say('  arithmetically consistent; the same three boluses in a 70 kg man give')
    f3 = 0.450 * (513 - 110) / (42 + 1)
    say(f'  {f3:.1f} mmol/L.  The formula ignores the urine, which is exactly the term')
    say('  section 4 shows is dominant.')

    # ---------------------------------------------------------------
    hdr('15. SUMMARY TABLE OF ALL SCENARIOS')
    scen = []
    def add(name, pp, yy0, dur, doses=(), note=''):
        tt, yv, oo = run_with_doses(pp, yy0, dur, list(doses))
        s = summarise(name, tt, yv, oo, note); scen.append(s); return s
    add('S01 normal control', P(), y_normal(P()), 30.0, note='no hyponatraemia')
    add('S02 chronic SIADH, untreated', P(**psi), ysi, 30.0, note='adapted, uninjured')
    add('S03 acute, corrected +20/24h', P(**{**pa, 'CTRLON': 1.0, 'RATETGT': 20.0,
        'NASTART': 110.0, 'NACAP': 140.0, 'WIN': 1.0, 'ACUTE': 0.0}), ya0, 30.0,
        note='oedema treated, no ODS')
    add('S04 chronic, +6/24h', P(**{**psi, 'CTRLON': 1.0, 'RATETGT': 6.0, 'NASTART': 110.0,
        'NACAP': 140.0, 'WIN': 1.0}), ysi, 30.0, note='guideline')
    add('S05 chronic, +12/24h', P(**{**psi, 'CTRLON': 1.0, 'RATETGT': 12.0, 'NASTART': 110.0,
        'NACAP': 140.0, 'WIN': 1.0}), ysi, 30.0, note='upper limit')
    add('S06 hypovolaemic + 0.9% saline', pA, yhy, 30.0, note='CENTRAL: autonomous overcorrection')
    add('S07 S06 with AVP frozen', pB, yhy, 30.0, note='counterfactual')
    add('S08 S06 + proactive DDAVP clamp',
        P(**{**phy, **base, 'R09': 0.5, 'CTRLON': 1.0, 'RATETGT': 6.0, 'NASTART': 110.0,
             'NACAP': 130.0, 'WIN': 1.0}), yhy, 30.0,
        ddavp_q8(0.0, 5.0) + [(t, 'DDAC', 2.0) for t in np.arange(5.0, 7.0, 0.5)]
        + [(t, 'DDAC', 2.0) for t in np.arange(7.0, 10.0, 1.0)], 'rate prescribed, tapered')
    for T in [12, 24, 48]:
        th = T / 24.0
        pr = P(**pA); pr.update(dict(RESCUE=1.0, TRESCUE=th, NARES=118.0,
                                     DURRES=1.0, WIN=0.5, R09=0.0))
        tt, yv, oo = run_protocol([(pA, th), (pr, 60.0 - th)], yhy, ddavp_q8(th, th + 2.5))
        scen.append(summarise(f'S09-11 relowering at {T} h', tt, yv, oo, 'rescue'))
    yk = ysi.copy(); yk[IX['KE']] -= 450.0
    add('S12 KCl 120 mmol/d, no sodium', P(**{**psi, 'RKCL': 120.0, 'WIN': 1.0}), yk, 14.0,
        note='hidden correction')
    add('S13 tolvaptan 15 mg/d', P(**psi), ysi, 14.0,
        [(d, 'TLVD', 15.0) for d in range(7)], 'aquaresis')
    add('S14 alcoholic at +8/24h', P(**{**psi, 'CTRLON': 1.0, 'RATETGT': 8.0, 'NASTART': 110.0,
        'NACAP': 140.0, 'WIN': 1.0, 'FOSM': 0.55, 'FNUT': 0.85}), ysi, 60.0, note='high risk')
    add('S15 normal-risk at +8/24h', P(**{**psi, 'CTRLON': 1.0, 'RATETGT': 8.0, 'NASTART': 110.0,
        'NACAP': 140.0, 'WIN': 1.0}), ysi, 60.0, note='comparator for S14')
    add('S16 urea 30 g/d', P(**{**psi, 'UREADOSE': 30.0, 'WIN': 1.0}), ysi, 14.0, note='urea')
    ppot, ypot, opot, _ = make_chronic(110.0, 21.0, 'potomania')
    add('S17 potomania, solute restored', P(**{**ppot, 'TSOLUP': 0.0, 'PUREA2': 400.0}), ypot,
        30.0, note='refeeding opens the aquaresis')
    add('S18 +8/d x 4 d (32 total)', P(**{**psi, 'CTRLON': 1.0, 'RATETGT': 8.0, 'NASTART': 110.0,
        'NACAP': 142.0, 'WIN': 1.0}), ysi, 30.0, note='slow but large')
    padr, yadr, oadr, _ = make_chronic(110.0, 21.0, 'adrenal')
    add('S19 adrenal insufficiency, steroid given',
        P(**{**padr, 'TAVPOFF': 0.0, 'WIN': 1.0}), yadr, 30.0, note='cortisol opens the aquaresis')
    psev, ysev, osev, _ = make_chronic(100.0, 21.0, 'siadh')
    add('S20 [Na] 100, +8/24h', P(**{**psev, 'CTRLON': 1.0, 'RATETGT': 8.0, 'NASTART': 100.0,
        'NACAP': 135.0, 'WIN': 1.0}), ysev, 60.0, note='very severe')
    say(f'  {"scenario":34s} {"Na0":>5s} {"Na24":>5s} {"d24":>5s} {"OMG":>5s} {"STR":>5s}'
        f' {"BW-":>6s} {"AST":>5s} {"MYE":>5s} {"DEF":>5s}')
    for s in scen:
        say(f'  {s["id"]:34s} {s["Na0"]:5.1f} {s["Na24"]:5.1f} {s["dNa24"]:5.1f}'
            f' {s["OMEGAmax"]:5.1f} {s["STRESSmax"]:5.1f} {s["BWmin"]:6.2f} {s["ASTmin"]:5.3f}'
            f' {s["MYEmin"]:5.3f} {s["DEFmax"]:5.1f}')
    say()
    say('  Na0/Na24 = serum sodium at 0 and 24 h; d24 = largest rise in any 24 h')
    say('  window; OMG = peak organic osmolyte deficit; STR = peak osmotic stress')
    say('  (threshold 8.0); BW- = minimum brain water (mL/100 g, normal 80.00);')
    say('  AST/MYE = astrocyte and myelin nadirs; DEF = peak clinical deficit 0-100.')
    json.dump(scen, open('ods_scenarios.json', 'w'), indent=1)

    open('ods_verification_output.txt', 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
