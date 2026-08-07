#!/usr/bin/env python3
"""
lep_verify_python.py — independent re-implementation of the leprosy QSP model
(lep_mrgsolve_model.R) in Python/scipy, used to (a) tune the parameters against
published anchors and (b) check that the mrgsolve equations do what the README
claims they do.

Run:  python3 lep_verify_python.py            (prints the anchor table)
      python3 lep_verify_python.py --tune     (prints extra diagnostics)

Nothing here is imported by the R model; it exists so that every number quoted
in README.md has been produced by an implementation written from the equations
rather than copied from the model that generated it.
"""
import argparse
import math
import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# parameters (days, mg, L; bacilli in units of 1e6 organisms per gram)
# ---------------------------------------------------------------------------
P = dict(
    # ---- rifampicin -------------------------------------------------------
    FRIF=0.90, KARIF=36.0, VRIF=50.0, CLRIF=240.0,
    KENZ=0.14, EMXIND=1.20, EC50IND=1.00,
    KMAXR=12.0, EC50R=0.30,
    # ---- dapsone ----------------------------------------------------------
    FDAP=0.93, KADAP=24.0, VDAP=75.0, CLDAP=43.2, FINDD=0.35,
    KMAXD=0.068, EC50D=0.50,
    FNOH=0.12, KNOH=12.0,
    KMET=6774.0, KRED=12.0, GRED=1.0,
    KINHB=0.10833, KOUTHB=0.008333, HEMO=18.0, GHEM=1.0,
    # ---- clofazimine ------------------------------------------------------
    FCLO=0.50, KACLO=6.0, VCLO1=100.0, VCLO2=3000.0, QCLO=100.0, CLCLO=43.0,
    A2REF=2100.0,
    KMAXC=0.090, EC50C=2.00,
    IMAXC=0.45, IC50C=0.80,
    KPIGIN=0.42, KPIGOUT=0.00385,
    # ---- prednisolone -----------------------------------------------------
    FPDN=0.85, KAPDN=48.0, VPDN=40.0, CLPDN=216.0,
    IC50P=0.020, IMAXP=0.90,
    KHPA=0.0476, HPAMAX=0.85,
    # ---- thalidomide ------------------------------------------------------
    FTHA=0.90, KATHA=14.4, VTHA=90.0, CLTHA=240.0,
    IMAXT=0.70, IC50T=0.50,
    # ---- ofloxacin + minocycline (lumped "ROM") ---------------------------
    FROM=0.95, KAROM=36.0, VROM=100.0, CLROM=192.0,
    KMAXO=0.150, EC50O=1.00,
    # ---- bacterial dynamics ----------------------------------------------
    MUMAX=0.05545,          # ln2 / 12.5 d generation time
    BMAX=190.0,             # carrying capacity (1e6 /g)
    KNAT=0.01043,           # natural death of growing bacilli (1/d)
    FNATP=0.096,            # persister death as a fraction of KNAT
    KCLR0=0.00168,          # dead-bacillus degradation, live-suppressed (1/d)
    CLRBOOST=2.00,          # de-repression of degradation once bacilli are dead
    KSUB=50.0,              # live burden at which suppression is half-relieved
    KGP=2.50e-8, KPG=0.0015, # growing <-> persister switching
    PTOL=0.002,             # persister drug tolerance multiplier
    FPERS=1e-5,             # persister fraction at diagnosis
    FRES=0.0,               # resistant fraction at diagnosis
    RRIF=0.0, RDAP=1.0, RCLO=1.0, RROM=1.0,   # resistant-clone susceptibility
    KHOST=2.00,             # maximal cell-mediated killing (1/d)
    KSAT=1.00,              # live burden that half-saturates granuloma containment
    PIMM=0.02,              # persister tolerance of cell-mediated killing
    # ---- antigen ----------------------------------------------------------
    YB=0.10, YD=0.90, KAGCL=0.35,
    FNERVE=0.05, KAGNCL=0.06,
    # ---- humoral / reaction ----------------------------------------------
    KINAB=1.00, KAB50=3.00, KOUTAB=0.0077,
    KFIC=0.0010, KCIC=2.00,
    KOUTT=12.0, STIC=1.20, SAG=1.50, KAGT=9.00,
    KINN=0.50, KOUTN=0.50,
    WTNF=0.60, WNEU=1.00, E50=2.00, HENL=4.0,
    # ---- cellular immunity -----------------------------------------------
    SPEC=0.05, CMIMAX=1.0, KCMI=0.0500, BSUP=1.00,
    KTREG=0.02, AG50T=4.50, TREGSUP=0.50,
    # ---- type 1 reaction --------------------------------------------------
    KT1R=0.1429, KAGN=0.30, PT1R=0.60,
    T1R50=12.0, HT1R=5.0,   # nerve damage needs a REACTION, not mere activity
    # ---- nerve ------------------------------------------------------------
    KDAM=0.70, WT1=1.00, WE2=0.50, ENLNEU=0.40,
    KREC=0.01667, SREC=2.00, KFIX=0.005556,
    # ---- misc -------------------------------------------------------------
    KLES=0.0714, KLES1=30.0,
    KALT=0.1429, ALT0=25.0, HRIF=0.50,
    TISSG=100.0,            # grams of bacilliferous tissue
    NCRIT=1.60e6,            # viable organisms per "escape" (Poisson relapse)
)

S = ['RIFG', 'RIFC', 'ENZ', 'DAPG', 'DAPC', 'NOH', 'METHB', 'HB',
     'CLOG', 'CLO1', 'CLO2', 'PDNG', 'PDNC', 'THAG', 'THAC', 'ROMG', 'ROMC',
     'BG', 'BP', 'BR', 'BD', 'AG', 'AGN', 'AB', 'IC', 'TNF', 'NEU',
     'CMI', 'TREG', 'T1R', 'NFIR', 'NFIP', 'LES', 'PIG', 'HPA', 'ALT',
     'ENLC', 'AGC']
IX = {n: i for i, n in enumerate(S)}
NST = len(S)


def derived(y, p):
    """Everything the ODEs and the output table share."""
    d = {}
    d['CR'] = y[IX['RIFC']] / p['VRIF']
    d['CD'] = y[IX['DAPC']] / p['VDAP']
    d['CN'] = y[IX['NOH']] / p['VDAP']
    d['CC'] = y[IX['CLO1']] / p['VCLO1']
    d['CLOE'] = y[IX['CLO2']] / p['A2REF']
    d['CP'] = y[IX['PDNC']] / p['VPDN']
    d['CT'] = y[IX['THAC']] / p['VTHA']
    d['CO'] = y[IX['ROMC']] / p['VROM']
    d['GROCC'] = d['CP'] / (p['IC50P'] + d['CP'])

    d['kR'] = p['KMAXR'] * d['CR'] / (p['EC50R'] + d['CR'])
    d['kD'] = p['KMAXD'] * d['CD'] / (p['EC50D'] + d['CD'])
    d['kC'] = p['KMAXC'] * d['CLOE'] / (p['EC50C'] + d['CLOE'])
    d['kO'] = p['KMAXO'] * d['CO'] / (p['EC50O'] + d['CO'])
    _blive0 = max(y[IX['BG']] + y[IX['BP']] + y[IX['BR']], 0.0)
    d['kHOST'] = p['KHOST'] * y[IX['CMI']] * p['KSAT'] / (p['KSAT'] + _blive0)

    BG, BP, BR, BD = y[IX['BG']], y[IX['BP']], y[IX['BR']], y[IX['BD']]
    d['BLIVE'] = max(BG + BP + BR, 0.0)
    d['BTOT'] = d['BLIVE'] + max(BD, 0.0)

    d['killS'] = d['kR'] + d['kD'] + d['kC'] + d['kO'] + d['kHOST']
    d['killR'] = (p['RRIF'] * d['kR'] + p['RDAP'] * d['kD'] +
                  p['RCLO'] * d['kC'] + p['RROM'] * d['kO'] + d['kHOST'])
    d['killP'] = p['PTOL'] * (d['kR'] + d['kD'] + d['kC'] + d['kO']) + \
        p['PIMM'] * d['kHOST']

    # degradation of dead bacilli is de-repressed once the live burden falls:
    # live organisms actively block phagosome maturation.
    d['KCLR'] = p['KCLR0'] * (1.0 + p['CLRBOOST'] *
                              (1.0 - d['BLIVE'] / (d['BLIVE'] + p['KSUB'])))

    d['MU'] = max(p['MUMAX'] * (1.0 - d['BLIVE'] / p['BMAX']), 0.0)

    d['DEATHFLUX'] = (p['KNAT'] * BG + p['FNATP'] * p['KNAT'] * BP + p['KNAT'] * BR
                      + d['killS'] * BG + d['killP'] * BP + d['killR'] * BR)
    d['RELEASE'] = p['YB'] * d['DEATHFLUX'] + p['YD'] * d['KCLR'] * max(BD, 0.0)

    # inhibition of TNF transcription by the three anti-inflammatory arms
    d['IPDN'] = 1.0 - p['IMAXP'] * d['GROCC']
    d['ITHA'] = 1.0 - p['IMAXT'] * d['CT'] / (p['IC50T'] + d['CT'])
    d['ICLO'] = 1.0 - p['IMAXC'] * d['CLOE'] / (p['IC50C'] + d['CLOE'])

    Edrive = p['WTNF'] * max(y[IX['TNF']] - 1.0, 0.0) + p['WNEU'] * max(y[IX['NEU']], 0.0)
    eh = Edrive ** p['HENL']
    d['ENL'] = 100.0 * eh / (p['E50'] ** p['HENL'] + eh)
    d['BI'] = math.log10(max(d['BTOT'], 1e-30)) + 3.0
    d['MI'] = 100.0 * d['BLIVE'] / max(d['BTOT'], 1e-30)
    return d


def rhs(t, y, p):
    d = derived(y, p)
    dy = np.zeros(NST)
    g = lambda n: y[IX[n]]

    # ---- rifampicin -------------------------------------------------------
    dy[IX['RIFG']] = -p['KARIF'] * g('RIFG')
    dy[IX['RIFC']] = p['KARIF'] * g('RIFG') * p['FRIF'] - \
        (p['CLRIF'] * g('ENZ') / p['VRIF']) * g('RIFC')
    dy[IX['ENZ']] = p['KENZ'] * (1.0 + p['EMXIND'] * d['CR'] /
                                 (p['EC50IND'] + d['CR']) - g('ENZ'))

    # ---- dapsone ----------------------------------------------------------
    CLD = p['CLDAP'] * (1.0 + p['FINDD'] * (g('ENZ') - 1.0))
    dy[IX['DAPG']] = -p['KADAP'] * g('DAPG')
    dy[IX['DAPC']] = p['KADAP'] * g('DAPG') * p['FDAP'] - (CLD / p['VDAP']) * g('DAPC')
    dy[IX['NOH']] = p['FNOH'] * CLD * d['CD'] - p['KNOH'] * g('NOH')
    dy[IX['METHB']] = p['KMET'] * d['CN'] - p['KRED'] * p['GRED'] * g('METHB')
    dy[IX['HB']] = p['KINHB'] - p['KOUTHB'] * g('HB') * \
        (1.0 + p['HEMO'] * d['CN'] * p['GHEM'])

    # ---- clofazimine ------------------------------------------------------
    dy[IX['CLOG']] = -p['KACLO'] * g('CLOG')
    dy[IX['CLO1']] = (p['KACLO'] * g('CLOG') * p['FCLO']
                      - (p['CLCLO'] / p['VCLO1']) * g('CLO1')
                      - (p['QCLO'] / p['VCLO1']) * g('CLO1')
                      + (p['QCLO'] / p['VCLO2']) * g('CLO2'))
    dy[IX['CLO2']] = (p['QCLO'] / p['VCLO1']) * g('CLO1') - \
        (p['QCLO'] / p['VCLO2']) * g('CLO2')

    # ---- prednisolone / thalidomide / ROM ---------------------------------
    dy[IX['PDNG']] = -p['KAPDN'] * g('PDNG')
    dy[IX['PDNC']] = p['KAPDN'] * g('PDNG') * p['FPDN'] - (p['CLPDN'] / p['VPDN']) * g('PDNC')
    dy[IX['THAG']] = -p['KATHA'] * g('THAG')
    dy[IX['THAC']] = p['KATHA'] * g('THAG') * p['FTHA'] - (p['CLTHA'] / p['VTHA']) * g('THAC')
    dy[IX['ROMG']] = -p['KAROM'] * g('ROMG')
    dy[IX['ROMC']] = p['KAROM'] * g('ROMG') * p['FROM'] - (p['CLROM'] / p['VROM']) * g('ROMC')

    # ---- bacterial populations -------------------------------------------
    BG, BP, BR, BD = g('BG'), g('BP'), g('BR'), g('BD')
    dy[IX['BG']] = (d['MU'] * BG - p['KNAT'] * BG - d['killS'] * BG
                    - p['KGP'] * BG + p['KPG'] * BP)
    dy[IX['BP']] = (p['KGP'] * BG - p['KPG'] * BP
                    - p['FNATP'] * p['KNAT'] * BP - d['killP'] * BP)
    dy[IX['BR']] = d['MU'] * BR - p['KNAT'] * BR - d['killR'] * BR
    dy[IX['BD']] = d['DEATHFLUX'] - d['KCLR'] * BD

    # ---- antigen ----------------------------------------------------------
    dy[IX['AG']] = d['RELEASE'] - p['KAGCL'] * g('AG')
    dy[IX['AGN']] = p['FNERVE'] * d['RELEASE'] - p['KAGNCL'] * g('AGN')
    dy[IX['AGC']] = d['RELEASE']

    # ---- humoral arm ------------------------------------------------------
    dy[IX['AB']] = p['KINAB'] * g('AG') / (p['KAB50'] + g('AG')) - p['KOUTAB'] * g('AB')
    dy[IX['IC']] = p['KFIC'] * g('AG') * g('AB') - p['KCIC'] * g('IC')

    # ---- TNF and neutrophils ---------------------------------------------
    drive = p['STIC'] * g('IC') + p['SAG'] * g('AG') / (p['KAGT'] + g('AG'))
    dy[IX['TNF']] = p['KOUTT'] * ((1.0 + drive) * d['IPDN'] * d['ITHA'] * d['ICLO']
                                  - g('TNF'))
    dy[IX['NEU']] = p['KINN'] * (g('IC') + 0.4 * max(g('TNF') - 1.0, 0.0)) * \
        d['IPDN'] * d['ITHA'] - p['KOUTN'] * g('NEU')

    # ---- cellular immunity (anergy lifts as antigen falls) ----------------
    # anergy is imposed by LIVE organisms (PGL-1, M2 polarisation, IL-10),
    # so it lifts within days of an effective bactericide - which is why the
    # reversal reaction is a treatment-onset phenomenon.
    cmi_t = p['CMIMAX'] * p['SPEC'] / (1.0 + d['BLIVE'] / p['BSUP']) * \
        (1.0 - p['TREGSUP'] * g('TREG'))
    dy[IX['CMI']] = p['KCMI'] * (cmi_t - g('CMI'))
    dy[IX['TREG']] = p['KTREG'] * (g('AG') / (p['AG50T'] + g('AG')) - g('TREG'))

    # ---- type 1 reaction --------------------------------------------------
    t1_t = 100.0 * g('CMI') * (g('AGN') / (p['KAGN'] + g('AGN'))) * \
        (1.0 - p['PT1R'] * d['GROCC'])
    dy[IX['T1R']] = p['KT1R'] * (t1_t - g('T1R'))

    # ---- nerve: a reversible pool draining into an irreversible sink ------
    th = g('T1R') ** p['HT1R']
    dmg = (p['WT1'] * th / (p['T1R50'] ** p['HT1R'] + th)
           + p['WE2'] * p['ENLNEU'] * d['ENL'] / 100.0)
    room = max(100.0 - g('NFIR') - g('NFIP'), 0.0) / 100.0
    dy[IX['NFIR']] = (p['KDAM'] * dmg * room
                      - p['KREC'] * (1.0 + p['SREC'] * d['GROCC']) * g('NFIR')
                      - p['KFIX'] * g('NFIR'))
    dy[IX['NFIP']] = p['KFIX'] * g('NFIR')

    # ---- descriptive states ----------------------------------------------
    les_t = 100.0 * (0.5 * d['BLIVE'] / (d['BLIVE'] + p['KLES1'])
                     + 0.3 * g('T1R') / 100.0 + 0.2 * d['ENL'] / 100.0)
    dy[IX['LES']] = p['KLES'] * (les_t - g('LES'))
    dy[IX['PIG']] = p['KPIGIN'] * d['CLOE'] - p['KPIGOUT'] * g('PIG')
    dy[IX['HPA']] = p['KHPA'] * ((1.0 - p['HPAMAX'] * d['GROCC']) - g('HPA'))
    dy[IX['ALT']] = p['KALT'] * (p['ALT0'] * (1.0 + p['HRIF'] * d['CR']) - g('ALT'))
    dy[IX['ENLC']] = d['ENL'] / 100.0
    return dy


# ---------------------------------------------------------------------------
# initial condition: a patient who has been incubating for years, at the
# quasi-steady state implied by the parameters (found by relaxation).
# ---------------------------------------------------------------------------
def initial_state(p, B0):
    y = np.zeros(NST)
    y[IX['ENZ']] = 1.0
    y[IX['HB']] = p['KINHB'] / p['KOUTHB']
    y[IX['HPA']] = 1.0
    y[IX['ALT']] = p['ALT0']
    y[IX['BG']] = B0 * (1.0 - p['FPERS'] - p['FRES'])
    y[IX['BP']] = B0 * p['FPERS']
    y[IX['BR']] = B0 * p['FRES']
    kclr = p['KCLR0'] * (1.0 + p['CLRBOOST'] * (1.0 - B0 / (B0 + p['KSUB'])))
    y[IX['BD']] = p['KNAT'] * B0 / kclr
    rel = p['YB'] * p['KNAT'] * B0 + p['YD'] * kclr * y[IX['BD']]
    y[IX['AG']] = rel / p['KAGCL']
    y[IX['AGN']] = p['FNERVE'] * rel / p['KAGNCL']
    y[IX['AB']] = (p['KINAB'] * y[IX['AG']] / (p['KAB50'] + y[IX['AG']])) / p['KOUTAB']
    y[IX['IC']] = p['KFIC'] * y[IX['AG']] * y[IX['AB']] / p['KCIC']
    y[IX['TREG']] = y[IX['AG']] / (p['AG50T'] + y[IX['AG']])
    y[IX['CMI']] = p['CMIMAX'] * p['SPEC'] / (1.0 + B0 / p['BSUP']) * \
        (1.0 - p['TREGSUP'] * y[IX['TREG']])
    drive = p['STIC'] * y[IX['IC']] + p['SAG'] * y[IX['AG']] / (p['KAGT'] + y[IX['AG']])
    y[IX['TNF']] = 1.0 + drive
    y[IX['NEU']] = p['KINN'] * (y[IX['IC']] + 0.4 * drive) / p['KOUTN']
    y[IX['T1R']] = 100.0 * y[IX['CMI']] * y[IX['AGN']] / (p['KAGN'] + y[IX['AGN']])
    return y


# ---------------------------------------------------------------------------
# dosing
# ---------------------------------------------------------------------------
def sched(cmt, amt, start, ii, n):
    return [(start + ii * k, cmt, amt) for k in range(int(n))]


def mdt_mb(months=12, start=0.0, clo_load=False, rom=False):
    """WHO multibacillary MDT: rifampicin 600 monthly, clofazimine 300 monthly
    + 50 daily, dapsone 100 daily."""
    ev = []
    ev += sched('RIFG', 600, start, 30, months)
    ev += sched('CLOG', 300, start, 30, months)
    ev += sched('CLOG', 50, start, 1, 30 * months)
    ev += sched('DAPG', 100, start, 1, 30 * months)
    if clo_load:                      # 300 mg daily for the first 30 days
        ev += sched('CLOG', 250, start, 1, 30)
    if rom:                           # accelerated-kill arm
        ev += sched('ROMG', 400, start, 1, 30 * months)
    return ev


def mdt_pb(months=6, start=0.0):
    ev = sched('RIFG', 600, start, 30, months)
    ev += sched('DAPG', 100, start, 1, 30 * months)
    return ev


def pred_taper(start, weeks=20, top=40.0, floor=5.0, step_wk=2):
    """Prednisolone taper. weeks=20, floor=5 reproduces the standard WHO
    20-week course for a type 1 reaction; longer courses hold a maintenance
    dose once the floor is reached."""
    ev = []
    dose = top
    for w in range(int(weeks)):
        ev += sched('PDNG', dose, start + 7 * w, 1, 7)
        if (w + 1) % step_wk == 0:
            dose = max(dose - 5.0, floor)
    return ev


def thal_course(start, days=84, mg=300.0):
    return sched('THAG', mg, start, 1, days)


def simulate(events, tend, p, B0, t_eval=None, y0=None):
    ev = sorted(events, key=lambda e: e[0])
    y = initial_state(p, B0) if y0 is None else y0.copy()
    if t_eval is None:
        t_eval = np.arange(0.0, tend + 0.5, 1.0)
    times = sorted(set([0.0] + [e[0] for e in ev if 0 <= e[0] <= tend] + [tend]))
    T, Y = [], []
    for a, b in zip(times[:-1], times[1:]):
        for (te, cmt, amt) in ev:
            if abs(te - a) < 1e-9:
                y[IX[cmt]] += amt
        pts = t_eval[(t_eval >= a) & (t_eval < b)]
        sol = solve_ivp(rhs, (a, b), y, args=(p,), method='LSODA',
                        rtol=1e-7, atol=1e-10, max_step=1.0,
                        t_eval=np.unique(np.concatenate([[a], pts, [b]])))
        T.append(sol.t[:-1]); Y.append(sol.y[:, :-1])
        y = sol.y[:, -1].copy()
    for (te, cmt, amt) in ev:
        if abs(te - times[-1]) < 1e-9:
            y[IX[cmt]] += amt
    T.append(np.array([times[-1]])); Y.append(y.reshape(-1, 1))
    return np.concatenate(T), np.concatenate(Y, axis=1)


def table(T, Y, p):
    out = {k: [] for k in ('BI', 'MI', 'ENL', 'T1R', 'NFIR', 'NFIP', 'NFI',
                           'BLIVE', 'BTOT', 'AGRATE', 'CLOE', 'METHB', 'HB',
                           'PIG', 'TNF', 'CR', 'CD', 'CC', 'CP', 'HPA')}
    for j in range(Y.shape[1]):
        d = derived(Y[:, j], p)
        out['BI'].append(d['BI']); out['MI'].append(d['MI'])
        out['ENL'].append(d['ENL']); out['T1R'].append(Y[IX['T1R'], j])
        out['NFIR'].append(Y[IX['NFIR'], j]); out['NFIP'].append(Y[IX['NFIP'], j])
        out['NFI'].append(Y[IX['NFIR'], j] + Y[IX['NFIP'], j])
        out['BLIVE'].append(d['BLIVE']); out['BTOT'].append(d['BTOT'])
        out['AGRATE'].append(d['RELEASE']); out['CLOE'].append(d['CLOE'])
        out['METHB'].append(Y[IX['METHB'], j]); out['HB'].append(Y[IX['HB'], j])
        out['PIG'].append(Y[IX['PIG'], j]); out['TNF'].append(Y[IX['TNF'], j])
        out['CR'].append(d['CR']); out['CD'].append(d['CD'])
        out['CC'].append(d['CC']); out['CP'].append(d['CP'])
        out['HPA'].append(Y[IX['HPA'], j])
    return {k: np.array(v) for k, v in out.items()}


def relapse_risk(blive, p):
    """Poisson single-escape interpretation of the residual viable burden."""
    n = blive * 1e6 * p['TISSG']
    return 1.0 - math.exp(-n / p['NCRIT'])


def par(**kw):
    q = dict(P); q.update(kw); return q


# ---------------------------------------------------------------------------
CHECKS = []


def check(name, value, lo, hi, unit='', note=''):
    ok = (value >= lo) and (value <= hi)
    CHECKS.append((name, value, lo, hi, unit, ok, note))
    return ok


def main(tune=False):
    LL = par(SPEC=0.05)     # lepromatous pole, BI ~6 at diagnosis
    BL = par(SPEC=0.35)     # borderline, BI ~5
    TT = par(SPEC=0.95)     # tuberculoid pole, smear negative
    B_LL, B_BL, B_TT = 150.0, 20.0, 1e-3

    # ---------------- A. drug PK anchors ----------------------------------
    p = par()
    T, Y = simulate(sched('RIFG', 600, 0, 30, 1), 2.0, p, 150.0,
                    t_eval=np.arange(0, 2.001, 0.005))
    cr = Y[IX['RIFC']] / p['VRIF']
    check('rifampicin Cmax after 600 mg', cr.max(), 6.0, 14.0, 'mg/L',
          'observed 8-12 mg/L')
    half = np.interp(0.5, [0, 1], [0, 1])
    i = np.argmax(cr)
    t_half = T[i:][np.argmin(np.abs(cr[i:] - cr.max() / 2))] - T[i]
    check('rifampicin terminal t1/2', t_half * 24, 2.5, 5.0, 'h', 'observed 3-4 h')

    T, Y = simulate(sched('DAPG', 100, 0, 1, 60), 60.0, p, 150.0)
    cd = Y[IX['DAPC']] / p['VDAP']
    check('dapsone Css (100 mg/d)', cd[-1], 1.2, 3.2, 'mg/L', 'observed 1.5-3 ug/mL')
    check('methaemoglobin on dapsone 100 mg/d', Y[IX['METHB']][-1], 4.0, 12.0, '%',
          'observed 5-12%')
    check('haemoglobin fall on dapsone 100 mg/d', 13.0 - Y[IX['HB']][-1], 0.7, 2.0,
          'g/dL', 'observed ~1-1.5 g/dL')

    g6 = par(GRED=0.35, GHEM=6.0)
    T, Y = simulate(sched('DAPG', 100, 0, 1, 60), 60.0, g6, 150.0)
    check('methaemoglobin, G6PD-deficient', Y[IX['METHB']][-1], 12.0, 35.0, '%')
    check('haemoglobin, G6PD-deficient', Y[IX['HB']][-1], 6.0, 10.0, 'g/dL')

    # clofazimine accumulation: fraction of steady state at 30 and 180 days
    T, Y = simulate(mdt_mb(24), 720.0, p, 150.0)
    cloe = Y[IX['CLO2']] / p['A2REF']
    ss = cloe[-1]
    check('clofazimine plasma Css', (Y[IX['CLO1']] / p['VCLO1'])[540], 0.4, 1.2,
          'mg/L', 'observed 0.5-0.9 ug/mL')
    check('clofazimine tissue depot at 30 d', 100 * cloe[30] / ss, 15.0, 40.0, '% of SS',
          'one 70-day half-life is ~1 month of loading')
    check('clofazimine tissue depot at 180 d', 100 * cloe[180] / ss, 70.0, 95.0, '% of SS')
    check('clofazimine pigmentation at 12 mo', Y[IX['PIG']][360], 40.0, 90.0, 'index')

    # ---------------- B. natural history ----------------------------------
    y0 = initial_state(LL, B_LL)
    T, Y = simulate([], 1825.0, LL, B_LL)
    tb = table(T, Y, LL)
    check('untreated LL bacterial index (5 y)', tb['BI'][-1], 5.4, 6.6, 'BI')
    check('untreated LL morphological index', tb['MI'][-1], 5.0, 30.0, '%',
          'observed 5-25% at diagnosis')

    T, Y = simulate([], 1825.0, TT, B_TT)
    tb = table(T, Y, TT)
    check('untreated TT: self-clearing burden', tb['BI'][-1], -30.0, 1.0, 'BI',
          'tuberculoid pole is smear negative')

    # ---------------- C. the two clocks -----------------------------------
    # C1. BI decline on MDT-MB
    T, Y = simulate(mdt_mb(12), 1095.0, LL, B_LL)
    mdt = table(T, Y, LL)
    dBI = mdt['BI'][0] - mdt['BI'][365]
    check('BI fall in year 1 of MDT-MB', dBI, 0.5, 1.2, 'log10/yr',
          'observed 0.6-1.0 log/yr')
    check('MI at 21 days of MDT-MB', mdt['MI'][21], 0.0, 1.0, '%',
          'observed: MI zero within 2-3 weeks of rifampicin')
    check('BI still positive 3 years after starting MDT', mdt['BI'][1095],
          2.0, 5.0, 'BI', 'smears stay positive for years after cure')

    # C2. the discrimination experiment: rifampicin-containing vs dapsone alone
    T, Yd = simulate(sched('DAPG', 100, 0, 1, 365), 365.0, LL, B_LL)
    dds = table(T, Yd, LL)
    dBI_dds = dds['BI'][0] - dds['BI'][365]
    check('BI fall in year 1, dapsone monotherapy', dBI_dds, 0.5, 1.4, 'log10/yr')
    check('BI gap MDT minus dapsone at 1 y', abs(dBI - dBI_dds), 0.0, 0.35, 'log10',
          'the BI cannot tell the two regimens apart')
    mi_gap = dds['MI'][21] - mdt['MI'][21]
    check('MI gap at 21 d (dapsone minus MDT)', mi_gap, 5.0, 30.0, '%',
          'the MI can')
    check('viable-burden gap at 21 d (log10)',
          math.log10(dds['BLIVE'][21] / max(mdt['BLIVE'][21], 1e-30)), 2.0, 12.0,
          'log10', 'the quantity the BI is blind to')
    check('MI reaches zero on dapsone alone by', float(np.argmax(dds['MI'] < 1.0)),
          100.0, 250.0, 'd', 'observed ~5 months')

    # ---------------- D. ENL: the derivative, not the level ---------------
    SLOWKILL = par(SPEC=0.05, KMAXR=0.45)   # same regimen, 20x weaker bactericide
    T, Ys = simulate(mdt_mb(12), 730.0, LL, B_LL)
    slow = table(T, Ys, LL)
    T, Yf = simulate(mdt_mb(12), 730.0, SLOWKILL, B_LL)
    fast = table(T, Yf, SLOWKILL)
    check('day of peak antigen liberation on MDT',
          float(np.argmax(slow['AGRATE'])), 0.0, 120.0, 'd',
          'the antigen surge is early, and it decays with the dead-bacillus pool')
    check('antigen liberation rate at 12 mo, as % of its peak',
          100 * slow['AGRATE'][360] / slow['AGRATE'].max(), 5.0, 40.0, '%')
    check('MDT raises the antigen liberation rate by',
          slow['AGRATE'].max() / slow['AGRATE'][0], 2.0, 4.0, 'x',
          'killing does not create antigen, it un-blocks its degradation')
    _, Ys2 = simulate(mdt_mb(12), 2500.0, LL, B_LL)
    _, Yf2 = simulate(mdt_mb(12), 2500.0, SLOWKILL, B_LL)
    ag_slow = Ys2[IX['AGC']][-1]
    ag_fast = Yf2[IX['AGC']][-1]
    check('cumulative antigen, fast vs 20x slower kill',
          100 * abs(ag_fast - ag_slow) / ag_slow, 0.0, 20.0, '% difference',
          'CONSERVATION: area set by B0 (residual gap = extra replication in the slow arm)')
    check('ENL peak, fast vs 20x slower bactericide',
          slow['ENL'].max() / max(fast['ENL'].max(), 1e-9), 1.15, 8.0, 'x',
          'the peak is NOT conserved')
    check('BI at 12 mo, fast vs 20x slower bactericide',
          abs(fast['BI'][360] - slow['BI'][360]), 0.0, 0.30, 'log10',
          'a 20-fold kill-rate difference is invisible on the BI')

    # ENL is a lepromatous-pole event and T1R a borderline-pole event
    _, Yll = simulate(mdt_mb(12), 730.0, LL, B_LL); tll = table(_, Yll, LL)
    _, Ybl = simulate(mdt_mb(12), 730.0, BL, B_BL); tbl = table(_, Ybl, BL)
    _, Ytt = simulate(mdt_pb(6), 730.0, TT, B_TT); ttt = table(_, Ytt, TT)
    check('peak ENL score, lepromatous patient', tll['ENL'].max(), 25.0, 100.0, '/100',
          'ENL belongs to the lepromatous pole')
    check('peak ENL score, tuberculoid patient', ttt['ENL'].max(), 0.0, 3.0, '/100',
          'and essentially never occurs at the tuberculoid pole')
    check('peak T1R activity, borderline vs lepromatous',
          tbl['T1R'].max() / max(tll['T1R'].max(), 1e-9), 1.3, 50.0, 'x',
          'reversal reaction belongs to the borderline zone')

    # ---------------- E. clofazimine timing -------------------------------
    T, Yn = simulate(mdt_mb(12), 730.0, LL, B_LL)
    T, Yl = simulate(mdt_mb(12, clo_load=True), 730.0, LL, B_LL)
    tn, tl = table(T, Yn, LL), table(T, Yl, LL)
    check('clofazimine depot at 30 d, loaded vs standard',
          tl['CLOE'][30] / tn['CLOE'][30], 2.5, 7.0, 'x')
    check('ENL burden reduction from a 1-month clofazimine load',
          100 * (Yn[IX['ENLC']][-1] - Yl[IX['ENLC']][-1]) / Yn[IX['ENLC']][-1],
          10.0, 75.0, '%', 'PREDICTION - never tested in a trial')

    # ---------------- F. the nerve window ---------------------------------
    # A borderline patient whose reversal reaction begins at day ~20 of MDT.
    # Two separate questions: how LONG must the steroid course be, and how
    # much does DELAY cost?
    base = mdt_mb(12)
    T0 = 20.0

    def nfi_perm(ev, end=1200.0):
        _, Y = simulate(ev, end, BL, B_BL)
        return float(Y[IX['NFIP']][-1])

    _, Yr = simulate(base, 1200.0, BL, B_BL)
    tr = table(_, Yr, BL)
    check('day of peak reversal-reaction activity', float(np.argmax(tr['T1R'])),
          20.0, 180.0, 'd', 'observed: reactions cluster in the first 6 months')
    never = nfi_perm(base)
    check('permanent nerve deficit, no steroid', never, 5.0, 40.0, 'points')

    dur = {w: nfi_perm(base + pred_taper(T0, weeks=w, floor=10.0))
           for w in (20, 32, 52)}
    check('permanent deficit, 20-week WHO taper at onset', dur[20], 0.0, never, 'points')
    check('deficit ratio, 20-week course vs 52-week course',
          dur[20] / max(dur[52], 1e-9), 1.3, 6.0, 'x',
          'Cochrane: a 20-week course leaves residual deficit; longer is better')
    check('duration monotonicity (52 < 32 < 20 weeks)',
          1.0 if dur[52] < dur[32] < dur[20] else 0.0, 1.0, 1.0, 'bool')

    res = {d: nfi_perm(base + pred_taper(T0 + d, weeks=52, floor=10.0))
           for d in (0, 30, 90, 180, 365)}
    for delay, lo, hi in ((30, 1.05, 2.5), (90, 1.8, 5.5), (180, 3.0, 9.0)):
        check(f'permanent deficit ratio, {delay} d delay vs immediate',
              res[delay] / max(res[0], 1e-9), lo, hi, 'x')
    check('deficit after a 1-year delay, as % of never treating',
          100 * res[365] / never, 85.0, 105.0, '%',
          'past the reaction, steroids have nothing left to save')
    # half-life of the salvageable deficit
    sal = {d: (never - v) / (never - res[0]) for d, v in res.items()}
    xs = sorted(sal); ys = [sal[x] for x in xs]
    thalf = float(np.interp(0.5, [-y for y in ys][::-1], xs[::-1])) if ys[0] > 0.5 else 0.0
    thalf = float(np.interp(-0.5, [-y for y in ys], xs))
    check('half-life of the salvageable deficit', thalf, 45.0, 220.0, 'd',
          'THE WINDOW: half of what steroids could have saved is gone by then')

    # ---------------- G. relapse ------------------------------------------
    out = {}
    for months, key in ((6, 'U-MDT 6 mo'), (12, 'MDT 12 mo'), (24, 'MDT 24 mo')):
        T, Y = simulate(mdt_mb(months), 30 * months + 1, LL, B_LL)
        blive = derived(Y[:, -1], LL)['BLIVE']
        out[key] = 100 * relapse_risk(blive, LL)
    check('relapse risk, MDT 12 months (BI 6)', out['MDT 12 mo'], 0.3, 4.0, '%',
          'observed 0.7-3% at 10 y')
    check('relapse risk, uniform MDT 6 months (BI 6)', out['U-MDT 6 mo'], 1.0, 12.0, '%',
          'observed ~2-3% at 5 y, higher than 12 mo')
    check('relapse risk ratio, 6 mo vs 12 mo', out['U-MDT 6 mo'] / max(out['MDT 12 mo'], 1e-9),
          1.5, 20.0, 'x')
    check('relapse risk, MDT 24 months', out['MDT 24 mo'], 0.0, out['MDT 12 mo'], '%')

    # low-BI patient on 6 months: should be safe
    T, Y = simulate(mdt_pb(6), 181, par(SPEC=0.60), B_TT)
    check('relapse risk, PB patient after 6 months', 100 * relapse_risk(
        derived(Y[:, -1], par(SPEC=0.60))['BLIVE'], P), 0.0, 1.0, '%')

    # ---------------- H. containment threshold and resistance -------------
    # Cell-mediated killing at a collapsed burden is KHOST*SPEC; net growth of
    # a replicating clone is MUMAX - KNAT.  The model therefore has a sharp
    # immunological set-point below which no residual clone can be contained.
    spec_star = (P['MUMAX'] - P['KNAT']) / P['KHOST']
    check('containment threshold (immunological set-point)', spec_star,
          0.010, 0.050, 'SPEC',
          'below it a residual clone regrows; above it the patient stays cured')
    lo_spec, hi_spec = 0.001, 0.20
    for _ in range(24):                       # bisection on the actual ODEs
        mid = 0.5 * (lo_spec + hi_spec)
        pm = par(SPEC=mid)
        _, Ym = simulate(mdt_mb(12), 3000.0, pm, B_LL)
        grew = derived(Ym[:, -1], pm)['BLIVE'] > derived(Ym[:, 400], pm)['BLIVE']
        if grew:
            lo_spec = mid
        else:
            hi_spec = mid
    check('bisection on the ODEs agrees with the algebraic threshold',
          100 * abs(0.5 * (lo_spec + hi_spec) - spec_star) / spec_star, 0.0, 12.0, '%')

    # Historical experiment: dapsone monotherapy in a polar (non-upgrading)
    # patient, with a pre-existing folP1 mutant at 1e-6.
    anerg = par(SPEC=0.010, FRES=1e-6, RRIF=1.0, RDAP=0.0, RCLO=1.0, RROM=1.0)
    T, Y = simulate(sched('DAPG', 100, 0, 1, 3650), 3650.0, anerg, B_LL)
    check('dapsone monotherapy: resistant clone at 10 y (log10/g)',
          math.log10(max(Y[IX['BR']][-1], 1e-30) * 1e6), 5.0, 12.0, 'log10/g',
          'secondary dapsone resistance, the reason MDT exists')
    check('dapsone monotherapy: susceptible clone at 10 y (log10/g)',
          math.log10(max(Y[IX['BG']][-1], 1e-30) * 1e6), -30.0, 4.0, 'log10/g')
    T, Y = simulate(mdt_mb(12), 3650.0, anerg, B_LL)
    check('MDT in the same patient: resistant clone at 10 y (log10/g)',
          math.log10(max(Y[IX['BR']][-1], 1e-30) * 1e6), -30.0, 5.0, 'log10/g',
          'the companion drugs cover the dapsone-resistant mutant')

    # ---------------- I. single-dose rifampicin PEP -----------------------
    T, Y = simulate(sched('RIFG', 600, 0, 30, 1), 7.0, LL, 1e-6)
    check('SDR-PEP: log10 kill of a subclinical inoculum at 7 d',
          math.log10(1e-6 / max(derived(Y[:, -1], LL)['BLIVE'], 1e-30)), 2.5, 8.0,
          'log10', 'one 600 mg dose kills ~99.9% of viable M. leprae')

    # ---------------- report ----------------------------------------------
    w = max(len(c[0]) for c in CHECKS)
    npass = sum(1 for c in CHECKS if c[5])
    print('=' * (w + 62))
    print('LEPROSY QSP MODEL — independent Python re-implementation'.center(w + 62))
    print('=' * (w + 62))
    for name, val, lo, hi, unit, ok, note in CHECKS:
        flag = 'PASS' if ok else 'FAIL'
        print(f'{name:<{w}}  {val:>12.4g} {unit:<10} [{lo:g}, {hi:g}]  {flag}'
              + (f'   {note}' if note else ''))
    print('-' * (w + 62))
    print(f'{npass}/{len(CHECKS)} anchors pass')

    if tune:
        print('\n--- extra diagnostics -------------------------------------')
        T, Y = simulate(mdt_mb(12), 1095.0, LL, B_LL)
        tb = table(T, Y, LL)
        for d in (0, 30, 90, 180, 365, 730, 1095):
            print(f'  day {d:>5}: BI {tb["BI"][d]:5.2f}  MI {tb["MI"][d]:7.3f}%  '
                  f'ENL {tb["ENL"][d]:5.1f}  AGrate {tb["AGRATE"][d]:6.2f}  '
                  f'CLOE {tb["CLOE"][d]:5.2f}  NFI {tb["NFI"][d]:5.1f}')
        print('  nerve window:', {k: round(v, 2) for k, v in res.items()})
        print('  relapse risk:', {k: round(v, 2) for k, v in out.items()})
    return 0 if npass == len(CHECKS) else 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tune', action='store_true')
    a = ap.parse_args()
    raise SystemExit(main(a.tune))
