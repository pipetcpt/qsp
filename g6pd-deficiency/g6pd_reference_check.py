#!/usr/bin/env python3
"""
g6pd_reference_check.py
=======================

An INDEPENDENT re-implementation of `g6pd_mrgsolve_model.R` in pure Python
(no numpy, no scipy, no R) with a fixed-step RK4 integrator.

Why this file exists
--------------------
Every number quoted in `README.md` and in the scenario comments of the R
model is produced here. Two implementations written from the same equations
but with different integrators, different languages and different state
orderings will not agree by accident, so running both is a real check on
the model rather than a restatement of it.

It also means the model's claims are reproducible by anyone with a Python 3
interpreter and nothing else installed:

    python3 g6pd_reference_check.py            # scenario table
    python3 g6pd_reference_check.py --curves   # + the age-activity table

Runtime is about 3-6 minutes for the full scenario set (the 250-day drug-free
lead-in that every run needs in order to start from its OWN steady state is
most of it).

What is deliberately NOT here
-----------------------------
This script covers the mechanism. It does not reproduce the R model's
convenience layer (mrgsolve event objects, the mosaic mixer, the Shiny
bindings). Where the two differ in a way that matters, the R model is the
reference and this file is the check.
"""

import math
import sys

LN2 = math.log(2.0)

# =====================================================================
#  PARAMETERS  (identical to $PARAM in g6pd_mrgsolve_model.R)
# =====================================================================
P = dict(
    # --- the variant: these two numbers are the genotype ---
    E0=1.00, TAU=62.0, FRET=1.50,
    # --- red cell age structure ---
    LSPAN=120.0, TRET=1.20, MCH=29.70, BV=5.00, VPLAS=3.00, VRBCL=2.00,
    # --- redox capacity and damage ---
    VMAXOX=60.0, OXBASE=0.060, KHZ=1.00, KPIT=0.50,
    KEVMAX=0.80, HZ50=0.10, NEV=3.0,
    KIVMAX=8.00, OMEG50=20.0, NIV=3.0, OXREF=0.02, SPLEEN=1.0,
    # --- erythropoiesis ---
    HB0=15.0, EPO0=10.0, GEPO=5.40, EPOMAX=6000.0,
    PROD0=0.041667, EMAXERY=5.00, EC50ERY=60.0,
    TTM0=5.00, FSHIFT=0.60, FERY=1.00,
    # --- primaquine ---
    FPQ=0.96, KAPQ=36.0, VPQ=240.0, CLPQ=576.0, KMPQ=2.10, AS2D6=1.0, SPQ=6.00,
    # --- tafenoquine ---
    FTQ=0.90, KATQ=12.0, VTQ=1600.0, CLTQ=73.9, STQ=6.00,
    # --- dapsone ---
    FDP=0.90, KADP=12.0, VDP=70.0, CLDP=46.6,
    FMNOH=0.05, FNAT2=1.0, SDAP=1.50, KMETDAP=1000.0,
    # --- rasburicase / urate ---
    VRBX=6.0, CLRBX=4.75, KCATRBX=3.30, KMUR=0.10,
    VUR=35.0, KGENUR=3.00, CLUR=9.00, SRBX=0.005,
    # --- methylene blue ---
    VMB=1400.0, CLMB=4400.0, KMBRED=540.0, SMB=8.00,
    # --- fava ---
    FFV=0.30, KAFV=12.0, VFV=40.0, CLFV=333.0, SFV=0.80,
    # --- infection ---
    INFON=0.0, OXINF=0.15,
    # --- methaemoglobin ---
    KCYB5=12.0, FMETBAS=12.0, KMETOX=73.0,
    # --- bilirubin ---
    BILHB=34.0, VBIL=12.0, CLBIL=2.92, FUGT=1.0,
    NEO=0.0, UGTMAT0=0.03, KUGTMAT=0.025, PNA0=0.0, KENT=0.0, PHOTO=0.0,
    ALBM=0.60, KABIL=7.0e7,
    # --- intravascular handling ---
    HP0=1.30, KSYNHP=0.20, KBINDHP=200.0, STOIHP=1.50, CLRENHB=60.0,
    KTUBI=0.005, KTUBR=0.35, KTUB50=1.00,
    CRGEN=1000.0, CLCR=100.0, VCRD=42.0,
)

NBIN = 12
EREF = 0.5593   # normal-subject age-weighted mean activity (see the R model)

# state layout
iR, iZ, iM, iRET = 0, 12, 24, 27
iPQG, iPQC, iOXF, iTQG, iTQC, iDPG, iDPC = 28, 29, 30, 31, 32, 33, 34
iRBX, iURA, iMB, iFVG, iFVC = 35, 36, 37, 38, 39
iMET, iUCB, iFHB, iHPT, iTUB, iCRE = 40, 41, 42, 43, 44, 45
NSTATE = 46


def pos(x):
    return x if x > 0.0 else 0.0


def ebins(p):
    wb = p['LSPAN'] / NBIN
    return [p['E0'] * math.exp(-LN2 * wb * (i + 0.5) / p['TAU'])
            for i in range(NBIN)]


def deriv(t, y, p, Eb, Vb):
    d = [0.0] * NSTATE
    wb = p['LSPAN'] / NBIN
    kage = NBIN / p['LSPAN']

    # ---- 0. population summary -------------------------------------
    Rmat = sum(pos(y[iR + i]) for i in range(NBIN))
    Rret = pos(y[iRET])
    Rtot = Rmat + Rret
    Esum = sum(pos(y[iR + i]) * Eb[i] for i in range(NBIN)) \
        + Rret * p['E0'] * p['FRET']
    meanE = Esum / Rtot if Rtot > 1e-9 else 0.0
    nadphav = min(1.0, meanE / EREF)

    # ---- 1. drug concentrations ------------------------------------
    CPQ = y[iPQC] / p['VPQ']
    CTQ = y[iTQC] / p['VTQ']
    CDP = y[iDPC] / p['VDP']
    CRBX = y[iRBX] / p['VRBX']
    CMB = y[iMB] / p['VMB']
    CFV = y[iFVC] / p['VFV']
    CNHOH = p['FNAT2'] * p['FMNOH'] * CDP

    urc = pos(y[iURA]) / p['VUR']
    v_uox = p['KCATRBX'] * CRBX * urc / (p['KMUR'] + urc)
    h2o2 = v_uox * p['VUR']                       # mmol/d, whole body

    mb_red = CMB * nadphav
    mb_unred = CMB * (1.0 - nadphav)

    # ---- 2. total oxidant flux --------------------------------------
    OX = (p['OXBASE']
          + p['SPQ'] * pos(y[iOXF])
          + p['STQ'] * CTQ
          + p['SDAP'] * CNHOH
          + p['SRBX'] * h2o2 / p['VRBCL']
          + p['SMB'] * mb_unred
          + p['SFV'] * CFV
          + p['OXINF'] * p['INFON'])

    # ---- 3. per-bin redox, damage, removal ---------------------------
    hz50n = p['HZ50'] ** p['NEV']
    omegn = p['OMEG50'] ** p['NIV']
    lysEV = lysIV = 0.0
    for i in range(NBIN):
        ri = pos(y[iR + i])
        zi = pos(y[iZ + i])
        unbuf = pos(OX - Vb[i])
        d[iZ + i] = p['KHZ'] * unbuf - p['KPIT'] * zi

        zev = zi ** p['NEV'] if zi > 0 else 0.0
        kev = p['KEVMAX'] * p['SPLEEN'] * zev / (hz50n + zev) if zi > 0 else 0.0
        omeg = OX / (Vb[i] + p['OXREF'])
        omn = omeg ** p['NIV']
        kiv = p['KIVMAX'] * omn / (omegn + omn) if zi > p['HZ50'] else 0.0

        lysEV += kev * ri
        lysIV += kiv * ri
        inflow = (Rret / p['TRET']) if i == 0 else kage * pos(y[iR + i - 1])
        d[iR + i] = inflow - kage * ri - (kev + kiv) * ri
    senesce = kage * pos(y[iR + NBIN - 1])

    # ---- 4. erythropoiesis ------------------------------------------
    Hb = Rtot * p['MCH'] / 10.0
    epo = min(p['EPOMAX'], p['EPO0'] * (p['HB0'] / max(Hb, 1.0)) ** p['GEPO'])
    dE = pos(epo - p['EPO0'])
    drive = dE / (dE + p['EC50ERY'])
    prod = p['PROD0'] * p['FERY'] * (1.0 + p['EMAXERY'] * drive)
    if epo < p['EPO0']:
        prod = p['PROD0'] * p['FERY'] * (epo / p['EPO0'])
    kM = 3.0 / (p['TTM0'] * (1.0 - p['FSHIFT'] * drive))
    tretb = p['TRET'] * (1.0 + drive)
    d[iM] = prod - kM * y[iM]
    d[iM + 1] = kM * y[iM] - kM * y[iM + 1]
    d[iM + 2] = kM * y[iM + 1] - kM * y[iM + 2]
    d[iRET] = kM * y[iM + 2] - Rret / tretb

    # ---- 5. methaemoglobin ------------------------------------------
    kred = p['KCYB5'] + p['KMBRED'] * mb_red
    mform = (p['FMETBAS'] + p['KMETDAP'] * CNHOH
             + p['KMETOX'] * OX * (1.0 - nadphav))
    d[iMET] = mform * (1.0 - pos(y[iMET]) / 100.0) - kred * pos(y[iMET])

    # ---- 6. bilirubin -------------------------------------------------
    hb_g_d = (lysEV + lysIV + senesce) * p['BV'] * p['MCH']
    bilprod = p['BILHB'] * hb_g_d
    ugtmat = 1.0
    if p['NEO'] > 0.5:
        ugtmat = p['UGTMAT0'] + (1.0 - p['UGTMAT0']) * \
            (1.0 - math.exp(-p['KUGTMAT'] * (p['PNA0'] + t)))
    clb = p['CLBIL'] * p['FUGT'] * ugtmat + p['PHOTO']
    d[iUCB] = bilprod / (p['VBIL'] * 10.0) + p['KENT'] - clb * pos(y[iUCB])

    # ---- 7. free Hb / haptoglobin / tubule / creatinine ---------------
    iv_g_d = lysIV * p['BV'] * p['MCH']
    fhb, hp = pos(y[iFHB]), pos(y[iHPT])
    bind = p['KBINDHP'] * hp * fhb
    filt = p['CLRENHB'] * fhb
    d[iFHB] = iv_g_d / p['VPLAS'] - bind - filt
    d[iHPT] = p['KSYNHP'] * (p['HP0'] - hp) - p['STOIHP'] * bind
    d[iTUB] = p['KTUBI'] * filt - p['KTUBR'] * pos(y[iTUB])
    gfrf = 1.0 / (1.0 + pos(y[iTUB]) / p['KTUB50'])
    d[iCRE] = (p['CRGEN'] / 10.0 - p['CLCR'] * gfrf * pos(y[iCRE])) / p['VCRD']

    # ---- 8. drug PK ---------------------------------------------------
    d[iPQG] = -p['KAPQ'] * y[iPQG]
    d[iPQC] = p['KAPQ'] * y[iPQG] - (p['CLPQ'] / p['VPQ']) * y[iPQC]
    d[iOXF] = p['KMPQ'] * (p['AS2D6'] * CPQ - pos(y[iOXF]))
    d[iTQG] = -p['KATQ'] * y[iTQG]
    d[iTQC] = p['KATQ'] * y[iTQG] - (p['CLTQ'] / p['VTQ']) * y[iTQC]
    d[iDPG] = -p['KADP'] * y[iDPG]
    d[iDPC] = p['KADP'] * y[iDPG] - (p['CLDP'] / p['VDP']) * y[iDPC]
    d[iRBX] = -(p['CLRBX'] / p['VRBX']) * y[iRBX]
    d[iURA] = p['KGENUR'] - (p['CLUR'] / p['VUR']) * pos(y[iURA]) - h2o2
    d[iMB] = -(p['CLMB'] / p['VMB']) * y[iMB]
    d[iFVG] = -p['KAFV'] * y[iFVG]
    d[iFVC] = p['KAFV'] * y[iFVG] - (p['CLFV'] / p['VFV']) * y[iFVC]
    return d


# =====================================================================
#  VARIANTS
# =====================================================================
VARIANTS = {
    'normal':        (1.00, 62, 'B (wild type)'),
    'Aplus':         (0.90, 45, 'A+ (A376G)'),
    'Aminus':        (0.55, 13, 'A- (A376G+G202A)'),
    'Mahidol':       (0.35, 20, 'Mahidol (G487A)'),
    'Mediterranean': (0.05, 12, 'Mediterranean (C563T)'),
    'Canton':        (0.08, 13, 'Canton (G1376T)'),
    'ClassI':        (0.02,  7, 'Class I (CNSHA)'),
}


def patient(variant='normal', wt=70, ugt='6/6', cyp2d6=1.0, nat2='fast',
            spleen=1.0, marrow=1.0, neonate=False, pna=0.0):
    e0, tau, _ = VARIANTS[variant]
    p = dict(P)
    p.update(E0=e0, TAU=tau, AS2D6=cyp2d6,
             FNAT2=1.5 if nat2 == 'slow' else 1.0,
             FUGT={'6/6': 1.0, '6/7': 0.70, '7/7': 0.35}[ugt],
             SPLEEN=spleen, FERY=marrow,
             BV=0.0714 * wt, VPLAS=0.0429 * wt, VRBCL=0.0286 * wt)
    if neonate:
        p.update(NEO=1.0, PNA0=pna, ALBM=0.45,
                 BV=0.085 * wt, VPLAS=0.048 * wt, VRBCL=0.037 * wt,
                 VBIL=0.55 * wt, HB0=17.0, LSPAN=80.0, MCH=34.0,
                 PROD0=0.041667 * (120.0 / 80.0), KENT=0.8, OXBASE=0.15,
                 CRGEN=1000.0 * wt / 70, CLCR=100.0 * wt / 70, VCRD=42.0 * wt / 70)
    return p


def astar(p, ox):
    """The closed form the R model reports as ASTAR."""
    oxeff = ox - p['KPIT'] * p['HZ50']
    if oxeff <= 0:
        return p['LSPAN']
    if p['E0'] * p['VMAXOX'] <= oxeff:
        return 0.0
    return min(p['LSPAN'],
               (p['TAU'] / LN2) * math.log(p['E0'] * p['VMAXOX'] / oxeff))


def bfree_nM(tsb, albm, ka):
    """One-site bilirubin-albumin binding, solved exactly."""
    btot = tsb * 10.0 / 584.7 / 1000.0
    b = ka * (albm / 1000.0 - btot) + 1.0
    return 1e9 * (-b + math.sqrt(b * b + 4.0 * ka * btot)) / (2.0 * ka)


# =====================================================================
#  SIMULATION
# =====================================================================
def simulate(p, doses=(), days=60, lead=250.0, h=0.005, record=0.25):
    """doses: iterable of (time_after_t0, state_index, amount)."""
    wb = p['LSPAN'] / NBIN
    Eb = ebins(p)
    Vb = [p['VMAXOX'] * e for e in Eb]

    y = [0.0] * NSTATE
    for i in range(NBIN):
        y[iR + i] = p['PROD0'] * wb
    for k in range(3):
        y[iM + k] = p['PROD0'] * p['TTM0'] / 3.0
    y[iRET] = p['PROD0'] * p['TRET']
    y[iMET] = p['FMETBAS'] / p['KCYB5']
    y[iUCB] = 0.60
    y[iHPT] = p['HP0']
    y[iCRE] = 1.00
    y[iURA] = p.get('URATE0', p['KGENUR'] * p['VUR'] / p['CLUR'])

    sched = sorted((lead + t, idx, amt) for t, idx, amt in doses)
    si = 0
    out = []
    n = int(round((lead + days) / h))
    nrec = max(1, int(round(record / h)))
    for k in range(n):
        t = k * h
        while si < len(sched) and sched[si][0] <= t + h / 2:
            _, idx, amt = sched[si]
            y[idx] += amt
            si += 1
        if k % nrec == 0 and t >= lead - 1.0:
            out.append(snapshot(t - lead, y, p, Eb))
        k1 = deriv(t, y, p, Eb, Vb)
        y1 = [y[i] + h / 2 * k1[i] for i in range(NSTATE)]
        k2 = deriv(t + h / 2, y1, p, Eb, Vb)
        y2 = [y[i] + h / 2 * k2[i] for i in range(NSTATE)]
        k3 = deriv(t + h / 2, y2, p, Eb, Vb)
        y3 = [y[i] + h * k3[i] for i in range(NSTATE)]
        k4 = deriv(t + h, y3, p, Eb, Vb)
        y = [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(NSTATE)]
    out.append(snapshot(n * h - lead, y, p, Eb))
    return out


def snapshot(t, y, p, Eb):
    Rmat = sum(pos(y[iR + i]) for i in range(NBIN))
    Rret = pos(y[iRET])
    Rtot = Rmat + Rret
    Esum = sum(pos(y[iR + i]) * Eb[i] for i in range(NBIN)) \
        + Rret * p['E0'] * p['FRET']
    meanE = Esum / Rtot if Rtot > 1e-9 else 0.0
    nadphav = min(1.0, meanE / EREF)

    CTQ = y[iTQC] / p['VTQ']
    CDP = y[iDPC] / p['VDP']
    CMB = y[iMB] / p['VMB']
    CFV = y[iFVC] / p['VFV']
    urc = pos(y[iURA]) / p['VUR']
    v_uox = p['KCATRBX'] * (y[iRBX] / p['VRBX']) * urc / (p['KMUR'] + urc)
    OX = (p['OXBASE'] + p['SPQ'] * pos(y[iOXF]) + p['STQ'] * CTQ
          + p['SDAP'] * p['FNAT2'] * p['FMNOH'] * CDP
          + p['SRBX'] * (v_uox * p['VUR']) / p['VRBCL']
          + p['SMB'] * CMB * (1.0 - nadphav)
          + p['SFV'] * CFV + p['OXINF'] * p['INFON'])
    a = astar(p, OX)
    tsb = pos(y[iUCB])
    return dict(t=t, HB=Rtot * p['MCH'] / 10.0,
                RETPCT=100.0 * Rret / Rtot if Rtot > 1e-9 else 0.0,
                G6PDPCT=100.0 * meanE / EREF, OX=OX, ASTAR=a,
                ATRISK=100.0 * (1.0 - a / p['LSPAN']),
                METPCT=pos(y[iMET]), TSB=tsb,
                BFREE=bfree_nM(tsb, p['ALBM'], p['KABIL']),
                HPTG=pos(y[iHPT]), FHB=pos(y[iFHB]),
                CREA=pos(y[iCRE]), URATE=pos(y[iURA]) / p['VUR'] * 168.1 / 10.0)


def summarise(name, out):
    base = min(out, key=lambda r: abs(r['t']))
    post = [r for r in out if r['t'] >= -1e-9]
    nad = min(post, key=lambda r: r['HB'])
    drop = 100.0 * (base['HB'] - nad['HB']) / base['HB']
    end = post[-1]
    recov = (100.0 * (end['HB'] - nad['HB']) / (base['HB'] - nad['HB'])
             if base['HB'] - nad['HB'] > 0.05 else float('nan'))
    dt = post[1]['t'] - post[0]['t']
    return dict(scenario=name, Hb_base=base['HB'], Hb_nadir=nad['HB'],
                drop=drop, nadir_d=nad['t'], Hb_end=end['HB'], recov=recov,
                astar_min=min(r['ASTAR'] for r in post),
                atrisk_max=max(r['ATRISK'] for r in post),
                retic_max=max(r['RETPCT'] for r in post),
                assay_base=base['G6PDPCT'],
                assay_max=max(r['G6PDPCT'] for r in post),
                methb_max=max(r['METPCT'] for r in post),
                tsb_max=max(r['TSB'] for r in post),
                bfree_max=max(r['BFREE'] for r in post),
                hapto_min=min(r['HPTG'] for r in post),
                fhb_max=max(r['FHB'] for r in post),
                hburia_d=sum(dt for r in post if r['FHB'] > 0.10),
                crea_max=max(r['CREA'] for r in post))


# ---- dose builders (state index + amount, times relative to exposure) ----
def pq_daily(dose, days, f=0.96):
    return [(float(d), iPQG, f * dose) for d in range(days)]


def pq_weekly(dose, weeks, f=0.96):
    return [(7.0 * w, iPQG, f * dose) for w in range(weeks)]


def tafenoquine(dose, f=0.90):
    return [(0.0, iTQG, f * dose)]


def dapsone(dose, days, f=0.90):
    return [(float(d), iDPG, f * dose) for d in range(days)]


def rasburicase(mg, days):
    return [(float(d), iRBX, mg) for d in range(days)]


def methylene_blue(mg, at):
    return [(float(at), iMB, mg)]


def fava(mg, f=0.30):
    return [(0.0, iFVG, f * mg)]


# =====================================================================
#  SCENARIOS
# =====================================================================
def run_all(quick=False):
    rows = []
    D = 60 if not quick else 30

    def go(name, p, doses=(), days=D, lead=250.0, rec=0.25, h=0.005):
        out = simulate(p, doses, days=days, lead=lead, h=h, record=rec)
        rows.append(summarise(name, out))
        r = rows[-1]
        print(f"  {name:52s} nadir {r['Hb_nadir']:5.2f} "
              f"({r['drop']:5.1f}% at d{r['nadir_d']:5.1f})  end {r['Hb_end']:5.2f}")
        sys.stdout.flush()
        return out

    print("\n[1-3] the same prescription in three genotypes")
    go("1  A- - primaquine 30 mg daily x60 d",
       patient('Aminus'), pq_daily(30, 60))
    go("2  Mediterranean - primaquine 30 mg daily x60 d",
       patient('Mediterranean'), pq_daily(30, 60))
    go("3  normal - primaquine 30 mg daily x60 d",
       patient('normal'), pq_daily(30, 60))

    print("\n[4] the half-life is the toxicology")
    go("4a A- - tafenoquine 300 mg single dose",
       patient('Aminus'), tafenoquine(300))
    go("4b A- - primaquine 15 mg daily x14 d",
       patient('Aminus'), pq_daily(15, 14))

    print("\n[5] weekly dosing is a different mechanism")
    go("5  A- - primaquine 45 mg WEEKLY x8",
       patient('Aminus'), pq_weekly(45, 8), days=70)

    print("\n[6] favism")
    go("6  Mediterranean child 20 kg - fava bean meal",
       patient('Mediterranean', wt=20), fava(1000), days=40, rec=0.05)

    print("\n[7] dapsone, then the antidote")
    mbruns = {}
    for tag, var in (("7a normal", 'normal'), ("7b A-", 'Aminus')):
        mbruns[tag] = go(f"{tag} - dapsone 100 mg + methylene blue on d30",
                         patient(var), dapsone(100, 45) + methylene_blue(140, 30),
                         days=45, rec=0.02)
    mb_report(mbruns)

    print("\n[8] rasburicase: the tumour sets the oxidant dose")
    for tag, var, ur, gen in (
            ("8a normal - rasburicase, HIGH urate", 'normal', 41.7, 90.0),
            ("8b A- - rasburicase, HIGH urate", 'Aminus', 41.7, 90.0),
            ("8c Mediterranean - rasburicase, HIGH urate", 'Mediterranean', 41.7, 90.0),
            ("8d Mediterranean - rasburicase, LOW urate", 'Mediterranean', 16.7, 20.0)):
        p = patient(var); p['URATE0'] = ur; p['KGENUR'] = gen
        go(tag, p, rasburicase(14, 3), days=30, rec=0.05)
    p = patient('Mediterranean'); p['URATE0'] = 41.7; p['KGENUR'] = 90.0
    go("8e Mediterranean - allopurinol instead (no uricase)", p, (), days=30, rec=0.05)

    print("\n[9] the assay trap")
    go("9  A- - 10-day course, followed to 120 d",
       patient('Aminus'), pq_daily(30, 10), days=120)

    print("\n[10] neonatal jaundice: a PRODUCT of two hits")
    for tag, var, ugt in (("10a neither hit", 'normal', '6/6'),
                          ("10b G6PD only", 'Mediterranean', '6/6'),
                          ("10c UGT1A1 7/7 only", 'normal', '7/7'),
                          ("10d BOTH hits", 'Mediterranean', '7/7')):
        go(f"{tag} - neonate 3 kg, 12 d",
           patient(var, wt=3, ugt=ugt, neonate=True, pna=0.5),
           (), days=12, lead=0.0, rec=0.05)

    print("\n[11] remove the rescue")
    go("11 A- + aplastic crisis - primaquine 30 mg daily",
       patient('Aminus', marrow=0.15), pq_daily(30, 60))

    print("\n[12] CYP2D6 poor metaboliser: no bioactivation, no haemolysis")
    go("12 A- CYP2D6 PM - primaquine 30 mg daily",
       patient('Aminus', cyp2d6=0.0), pq_daily(30, 60))
    return rows


def mb_report(runs):
    """The summary table reports a MAXIMUM, which is exactly the wrong
    statistic for an antidote. Print what happens in the hours AFTER the
    methylene blue goes in."""
    print("\n  --- methylene blue 2 mg/kg given on day 30 ---")
    print(f"    {'':22s}" + "".join(f"{('d'+str(t)):>9s}"
          for t in ('29.9', '30.1', '30.5', '31.0', '32.0', '35.0')))
    for tag, out in runs.items():
        def at(tt):
            return min(out, key=lambda r: abs(r['t'] - tt))
        row = "".join(f"{at(t)['METPCT']:9.1f}"
                      for t in (29.9, 30.1, 30.5, 31.0, 32.0, 35.0))
        print(f"    {tag+' MetHb %':22s}{row}")
    for tag, out in runs.items():
        def at(tt):
            return min(out, key=lambda r: abs(r['t'] - tt))
        row = "".join(f"{at(t)['OX']:9.2f}"
                      for t in (29.9, 30.1, 30.5, 31.0, 32.0, 35.0))
        print(f"    {tag+' oxidant flux':22s}{row}")
    n, a = runs.get("7a normal"), runs.get("7b A-")
    if n and a:
        def mn(o, tt):
            return min(o, key=lambda r: abs(r['t'] - tt))
        print(f"    -> normal: MetHb {mn(n,29.9)['METPCT']:.1f}% -> "
              f"{mn(n,30.5)['METPCT']:.1f}% in 12 h (the antidote works)")
        print(f"    -> A-    : MetHb {mn(a,29.9)['METPCT']:.1f}% -> "
              f"{mn(a,30.5)['METPCT']:.1f}% in 12 h, and the oxidant flux rises "
              f"{mn(a,29.9)['OX']:.2f} -> {mn(a,30.1)['OX']:.2f}")
        print("       (no NADPH -> no leuco-MB -> the dye stays oxidised "
              "and becomes the insult)")


def age_table():
    print("\n" + "=" * 100)
    print("TABLE 1.  The genotype IS the age-activity curve.  E(a) = E0*exp(-ln2*a/TAU)")
    print("          No simulation here at all -- this is arithmetic on two numbers.")
    print("=" * 100)
    ages = [0, 20, 40, 60, 80, 100, 120]
    oxs = [("baseline", 0.06), ("PQ 30 mg", 0.372), ("fava", 6.0)]
    hdr = f"{'variant':24s} {'E0':>5s} {'TAU':>4s} " + \
        " ".join(f"{'E@'+str(a):>8s}" for a in ages) + " " + \
        " ".join(f"{'a*|'+n:>12s}" for n, _ in oxs)
    print(hdr); print("-" * len(hdr))
    for key, (e0, tau, lab) in VARIANTS.items():
        p = dict(P); p['E0'] = e0; p['TAU'] = tau
        es = " ".join(f"{e0*math.exp(-LN2*a/tau):8.2e}" for a in ages)
        aa = " ".join(f"{astar(p, o):12.1f}" for _, o in oxs)
        print(f"{lab:24s} {e0:5.2f} {tau:4.0f} {es} {aa}")
    print("\nRead the a* columns: at the primaquine-30 mg load a normal cell of any")
    print("age is safe (120 = nothing at risk), an A- cell is at risk past ~87 days")
    print("(the oldest quarter), and a Mediterranean cell past ~39 days (two thirds).")
    print("At a fava load not even a Mediterranean reticulocyte is safe (a* = 0).")


def main():
    if '--curves' in sys.argv:
        age_table()
    quick = '--quick' in sys.argv
    print("\n" + "=" * 100)
    print("TABLE 2.  Scenario results (independent RK4 re-implementation)")
    print("=" * 100)
    rows = run_all(quick=quick)
    print("\n" + "=" * 158)
    hdr = (f"{'scenario':52s}{'Hb0':>6s}{'nadir':>7s}{'drop%':>7s}{'day':>6s}"
           f"{'end':>7s}{'recov%':>8s}{'a*min':>7s}{'risk%':>7s}{'ret%':>6s}"
           f"{'assay0':>8s}{'assayX':>8s}{'MetHb':>7s}{'TSB':>6s}{'Bf':>6s}"
           f"{'Hp':>6s}{'fHb':>7s}{'HbUd':>6s}{'Cr':>6s}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        rc = "n/a" if r['recov'] != r['recov'] else f"{r['recov']:.0f}"
        print(f"{r['scenario']:52s}{r['Hb_base']:6.2f}{r['Hb_nadir']:7.2f}"
              f"{r['drop']:7.1f}{r['nadir_d']:6.1f}{r['Hb_end']:7.2f}{rc:>8s}"
              f"{r['astar_min']:7.1f}{r['atrisk_max']:7.1f}{r['retic_max']:6.1f}"
              f"{r['assay_base']:8.1f}{r['assay_max']:8.1f}{r['methb_max']:7.1f}"
              f"{r['tsb_max']:6.2f}{r['bfree_max']:6.1f}{r['hapto_min']:6.2f}"
              f"{r['fhb_max']:7.3f}{r['hburia_d']:6.1f}{r['crea_max']:6.2f}")
    print("=" * 158)
    print("Hb0 baseline g/dL | drop% from own baseline | recov% of the deficit "
          "regained by the end of follow-up |")
    print("a*min critical age (d) | risk% of red cell mass older than a* | "
          "assay0/assayX the LABORATORY reading |")
    print("Bf free bilirubin nM | Hp haptoglobin g/L | fHb free plasma Hb g/L "
          "| HbUd days of haemoglobinuria")


if __name__ == '__main__':
    main()
