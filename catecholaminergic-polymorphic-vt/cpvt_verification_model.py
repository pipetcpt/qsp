#!/usr/bin/env python3
"""
CPVT QSP prototype v2 -- independent re-implementation used to CALIBRATE and
VERIFY the mrgsolve model before writing it.

Structure (fixed after v1 exposed two defects):
  * every process that happens ONCE PER BEAT (ICaL entry, systolic SR release,
    SERCA uptake, NCX extrusion) is written per beat and multiplied by HR;
    only the leak and the spontaneous waves are per-minute processes.  v1
    scaled release with HR but not uptake, which made cytosolic Ca run away to
    770 uM and made the SR *deplete* during exercise.
  * the SOICR wave rate is a HARD threshold in X = CaJSR_free / Cth.  v1 used
    a sigmoid centred on X = 1, which gave a wild-type heart 11 waves/min at
    X = 0.75 and therefore exercise VT in a normal subject.
"""
import numpy as np
from scipy.integrate import solve_ivp

P = dict(
    # ---------------- sympathetic drive (time in MINUTES) ----------------
    KNE=6.0, KUPT=1.0, KEPI=3.0, KCLE=0.35,
    FLCSD=0.0, ELCSD=0.60,
    EC50NE=3.0, EC50EPI=1.5, FEPI=0.60,
    # ---------------- heart rate ----------------
    HRV0=55.0, HRVEX=45.0, HRMX=230.0, KMHR=4.0, TAUHR=0.35,
    # ---------------- cAMP / PKA ----------------
    KAC=4.0, KPDE=4.0, KPKAON=1.0, KPKAOFF=1.0, KMPKA=0.50,
    KPI1=3.0, KPI1OFF=1.5, TAUPH=0.5,
    KPLB=0.28, KPRYR=0.32, KPICA=0.30,
    # ---------------- Ca cycling: PER BEAT (uM cytosol / beat) ----------
    GCAL=12.0, AICA=0.90,
    VUP=46.2, KMUP=0.35, ASERCA=1.00,
    VNCX=62.5, KMNCX=0.60, ENCX=0.55,
    GREL=0.0398, KRELH=1000.0, HREL=2.0,
    GLEAK=0.045,                     # per MINUTE
    VTR=15.0, VNSR=0.055, VJSR=0.0057, BETAI=0.045,
    KMITI=1.5, KMITO=1.0, KMMIT=0.60,
    BCSQ=14000.0, KCSQ=630.0, KCSQTO=0.0020, CSQSET=1.0,
    # ---------------- RyR2 / SOICR ----------------
    CTH0=2330.0, GRYR=0.0, GCSQI=0.0,
    BPRYR=0.60, BCAMK=0.55,
    WMAX=110.0, KW=0.15, HW=2.0,     # hard-threshold wave rate
    FRELW=0.25, AFRW=2.0,            # fraction of the store a wave releases
    KREF=0.40, TAUREC=0.030, FREC=1.0,
    # ---------------- DAD -> triggered beat ----------------
    GDADV=1.64, GK1=1.0, TAUV=0.05,
    VREQ0=22.0, MTRIG=6.0, PNAV=1.5, KO=4.2,
    KVT=22.0, NVT_H=4.0,
    # ---------------- flecainide ----------------
    KONNA=0.50, KOFFNA=2.00, NAMIN=0.35,
    IC50FLR=4.0, EFLR=0.25, AQRS=0.55, QRS0=95.0,
    # ---------------- chronic ----------------
    KCMON=0.020, KCMOFF=0.010, CAI0=0.377,
    KFIB=2e-5, KFIBD=1e-5,
    # ---------------- hazard / ICD ----------------
    KHAZ=0.010, KSHK=0.90, TAUSHK=1.0, ESHK=6.0, ICD=0.0,
    # ---------------- PK ----------------
    KANAD=1.2, FNAD=0.30, VCNAD=150.0, VPNAD=150.0, QNAD=20.0, CLNAD=12.0,
    EC50NAD=20.0,
    KAMET=1.5, FMET=0.40, VCMET=250.0, VPMET=500.0, QMET=90.0, CLMET=150.0,
    EC50MET=45.0, CYP2D6=1.0,
    KABIS=1.0, FBIS=0.90, VCBIS=230.0, CLBIS=14.5, EC50BIS=12.0,
    KAFLE=0.9, FFLE=0.90, VCFLE=350.0, VPFLE=250.0, QFLE=40.0, CLFLE=35.0,
    FUFLE=0.60, MWFLE=414.0,
    KAVER=1.2, FVER=0.22, VCVER=300.0, CLVER=60.0, IC50VER=80.0, EMXVER=0.55,
    # reference (baseline) phospho-levels, filled by settle()
    PICA0=0.0, PPLB0=0.0, PRYR0=0.0, CFLEfree=0.0,
)

NAMES = """NEJ EPI HR CAMP PKA PI1 PPLB PRYR PICA
CAI CANSR CAJSRT CAMT RREF CSQ VDAD NAAV PVCR VTB NPVC NVT HZRD NSHOCK SYNCH
CAMK CAINT FIBR QRSD FATG SHKD
NADD NADC NADP METD METC METP BISD BISC FLED FLEC FLEP VERD VERC EXPO""".split()
IDX = {n: i for i, n in enumerate(NAMES)}
NS = len(NAMES)


H60 = 60.0          # hours -> minutes


def csq_free(total, bmax, kd):
    b = total - bmax - kd
    return 0.5 * (b + np.sqrt(b * b + 4.0 * kd * total))


def rhs(t, y, p, drive):
    g = lambda n: y[IDX[n]]
    d = np.zeros(NS)
    EXER, EMOT, EPIINF = drive(t)

    # -------- sympathetic drive -------------------------------------------
    SNS = min(1.0, 0.10 + 0.90 * EXER + 0.55 * EMOT)
    NEJ = max(g('NEJ'), 0.0)
    d[IDX['NEJ']] = p['KNE'] * SNS * (1 - p['FLCSD'] * p['ELCSD']) - p['KUPT'] * NEJ
    EPI = max(g('EPI'), 0.0)
    d[IDX['EPI']] = (p['KEPI'] * SNS * SNS + EPIINF + p['ESHK'] * g('SHKD')
                     - p['KCLE'] * EPI)

    # -------- competitive beta1 blockade ----------------------------------
    CNAD = max(g('NADC'), 0.0) / p['VCNAD'] * 1e3
    CMET = max(g('METC'), 0.0) / p['VCMET'] * 1e3
    CBIS = max(g('BISC'), 0.0) / p['VCBIS'] * 1e3
    BSHIFT = 1.0 + CNAD / p['EC50NAD'] + CMET / p['EC50MET'] + CBIS / p['EC50BIS']
    AGON = NEJ / p['EC50NE'] + p['FEPI'] * EPI / p['EC50EPI']
    SYMP = AGON / BSHIFT
    B1ACT = SYMP / (1.0 + SYMP)

    # -------- heart rate ---------------------------------------------------
    HRV = p['HRV0'] + p['HRVEX'] * EXER
    E = SYMP / (SYMP + p['KMHR'])
    d[IDX['HR']] = (HRV + (p['HRMX'] - HRV) * E - g('HR')) / p['TAUHR']
    HR = min(max(g('HR'), 30.0), 260.0)

    # -------- cAMP / PKA / phospho-targets --------------------------------
    d[IDX['CAMP']] = p['KAC'] * B1ACT - p['KPDE'] * g('CAMP')
    CAMP = max(g('CAMP'), 0.0)
    fpka = CAMP**2 / (CAMP**2 + p['KMPKA']**2)
    d[IDX['PKA']] = p['KPKAON'] * fpka * (1 - g('PKA')) - p['KPKAOFF'] * g('PKA')
    PKA = max(g('PKA'), 0.0)
    d[IDX['PI1']] = p['KPI1'] * PKA * (1 - g('PI1')) - p['KPI1OFF'] * g('PI1')
    PP1 = 1.0 / (1.0 + 1.5 * g('PI1'))
    for nm, km in (('PPLB', p['KPLB']), ('PRYR', p['KPRYR']), ('PICA', p['KPICA'])):
        tgt = PKA**2 / (PKA**2 + (km * PP1)**2)
        d[IDX[nm]] = (tgt - g(nm)) / p['TAUPH']

    # -------- Ca cycling ---------------------------------------------------
    CAI = max(g('CAI'), 1e-4)
    CANSR = max(g('CANSR'), 1.0)
    CAJSRT = max(g('CAJSRT'), 1.0)
    CSQ = min(max(g('CSQ'), 0.01), 2.0)
    d[IDX['CSQ']] = p['KCSQTO'] * (p['CSQSET'] - CSQ)
    CAJSRF = csq_free(CAJSRT, p['BCSQ'] * CSQ, p['KCSQ'])

    CVER = max(g('VERC'), 0.0) / p['VCVER'] * 1e3
    fICa = max((1 + p['AICA'] * (g('PICA') - p['PICA0']))
               * (1 - p['EMXVER'] * CVER / (CVER + p['IC50VER'])), 0.05)
    fSERCA = 1 + p['ASERCA'] * (g('PPLB') - p['PPLB0'])

    RREF = min(max(g('RREF'), 0.0), 0.98)
    fload = (CAJSRF / p['KRELH'])**p['HREL'] / (1 + (CAJSRF / p['KRELH'])**p['HREL'])

    # per-beat fluxes x HR
    JCAL = p['GCAL'] * fICa * HR
    JUP = p['VUP'] * fSERCA * CAI**2 / (CAI**2 + p['KMUP']**2) * HR
    JNCX = (p['VNCX'] * CAI**3 / (CAI**3 + p['KMNCX']**3)
            * 60.0 * (HR / 60.0)**p['ENCX'])
    JREL = p['GREL'] * (1 - RREF) * CAJSRF * fload * fICa * HR
    JMITI = p['KMITI'] * CAI**2 / (CAI**2 + p['KMMIT']**2) * HR
    JMITO = p['KMITO'] * max(g('CAMT'), 0.0) / 100.0 * HR
    # per-minute
    JLEAK = p['GLEAK'] * CAJSRF * (1 + p['BCAMK'] * g('CAMK'))

    # -------- SOICR: hard threshold ---------------------------------------
    CTHR = (p['CTH0'] * (1 - p['GRYR']) * (1 - p['GCSQI'])
            * (1 + p['EFLR'] * p['CFLEfree'] / (p['CFLEfree'] + p['IC50FLR']))
            / (1 + p['BPRYR'] * (g('PRYR') - p['PRYR0'])
               + p['BCAMK'] * g('CAMK')))
    X = CAJSRF / CTHR
    xs = max(X - 1.0, 0.0)
    # RyR2 refractoriness sets a CEILING on how often a wave can recur
    # (1/tau_recovery); it does NOT scale with heart rate.  Writing it as
    # a (1-RREF) multiplier -- as v2 did -- handed heart rate a spurious
    # ANTI-arrhythmic channel that made bradycardia and LCSD proarrhythmic.
    WCAP = 1.0 / (p['TAUREC'] * p['FREC'])
    WAVE = min(p['WMAX'] * xs**p['HW'] / (p['KW']**p['HW'] + xs**p['HW']), WCAP)
    FRW = p['FRELW'] * (1 + p['AFRW'] * max(1.0 - CSQ, 0.0))
    QRELc = FRW * CAJSRT * p['VJSR']          # uM of cytosol per wave
    JSPONT = WAVE * QRELc

    d[IDX['CAI']] = p['BETAI'] * (JCAL + JREL + JLEAK + JSPONT
                                  - JUP - JNCX - JMITI + JMITO)
    d[IDX['CAMT']] = (JMITI - JMITO) * 0.20
    JTR = p['VTR'] * (CANSR - CAJSRF)
    d[IDX['CANSR']] = (JUP - JTR) / p['VNSR']
    d[IDX['CAJSRT']] = (JTR - JREL - JLEAK - JSPONT) / p['VJSR']
    d[IDX['RREF']] = (p['KREF'] * (JREL + JSPONT) / 1000.0 * (1 - RREF)
                      - RREF / (p['TAUREC'] * p['FREC']))

    # -------- DAD -> triggered beat ---------------------------------------
    NAAV = min(max(g('NAAV'), p['NAMIN']), 1.0)
    d[IDX['NAAV']] = (p['KOFFNA'] * (1 - NAAV)
                      - p['KONNA'] * p['CFLEfree'] * (HR / 60.0) * NAAV)
    SINK = p['GK1'] * np.sqrt(p['KO'] / 5.4)
    VDADT = p['GDADV'] * QRELc / SINK * (1.0 if WAVE > 1e-9 else 0.0)
    d[IDX['VDAD']] = (VDADT - g('VDAD')) / p['TAUV']
    VDAD = max(g('VDAD'), 0.0)
    VREQ = p['VREQ0'] / NAAV**p['PNAV']
    PTRIG = VDAD**p['MTRIG'] / (VDAD**p['MTRIG'] + VREQ**p['MTRIG'])
    d[IDX['PVCR']] = (min(WAVE * PTRIG, HR) - g('PVCR')) / 0.10
    PVCR = max(g('PVCR'), 0.0)
    r = (PVCR / p['KVT'])**p['NVT_H']
    d[IDX['VTB']] = (r / (1 + r) - g('VTB')) / 0.10
    d[IDX['NPVC']] = PVCR
    d[IDX['NVT']] = g('VTB') * 2.0
    d[IDX['HZRD']] = p['KHAZ'] * g('VTB')
    d[IDX['SYNCH']] = p['KHAZ'] * 0.45 * g('VTB')
    dsh = p['ICD'] * p['KSHK'] * g('VTB')
    d[IDX['NSHOCK']] = dsh
    d[IDX['SHKD']] = dsh - g('SHKD') / p['TAUSHK']

    # -------- chronic / readouts ------------------------------------------
    d[IDX['CAMK']] = (p['KCMON'] * max(CAI / p['CAI0'] - 1.0, 0.0) * (1 - g('CAMK'))
                      - p['KCMOFF'] * g('CAMK'))
    d[IDX['CAINT']] = max(CAI / p['CAI0'] - 1.0, 0.0)
    d[IDX['FIBR']] = p['KFIB'] * g('CAMK') - p['KFIBD'] * g('FIBR')
    d[IDX['QRSD']] = (p['QRS0'] * (1 + p['AQRS'] * (1 - NAAV)) - g('QRSD')) / 0.5
    d[IDX['FATG']] = ((1 - E / 0.60) - g('FATG')) / 5.0
    d[IDX['EXPO']] = EXER

    # -------- PK -----------------------------------------------------------
    for dep, cen, per, ka, F, VC, VP, Q, CL in (
        ('NADD', 'NADC', 'NADP', 'KANAD', 'FNAD', 'VCNAD', 'VPNAD', 'QNAD', 'CLNAD'),
        ('METD', 'METC', 'METP', 'KAMET', 'FMET', 'VCMET', 'VPMET', 'QMET', 'CLMET'),
        ('FLED', 'FLEC', 'FLEP', 'KAFLE', 'FFLE', 'VCFLE', 'VPFLE', 'QFLE', 'CLFLE'),
    ):
        # NB: ka, CL and Q are given in 1/h and L/h; the model runs in
        # MINUTES, hence the /H60 conversion (this was a real defect in v2).
        cl = p[CL] * (p['CYP2D6'] if cen == 'METC' else 1.0) / H60
        ain = p[ka] / H60 * g(dep)
        d[IDX[dep]] = -ain
        d[IDX[cen]] = (p[F] * ain - cl / p[VC] * g(cen)
                       - p[Q] / H60 * (g(cen) / p[VC] - g(per) / p[VP]))
        d[IDX[per]] = p[Q] / H60 * (g(cen) / p[VC] - g(per) / p[VP])
    d[IDX['BISD']] = -p['KABIS'] / H60 * g('BISD')
    d[IDX['BISC']] = (p['FBIS'] * p['KABIS'] / H60 * g('BISD')
                      - p['CLBIS'] / H60 / p['VCBIS'] * g('BISC'))
    d[IDX['VERD']] = -p['KAVER'] / H60 * g('VERD')
    d[IDX['VERC']] = (p['FVER'] * p['KAVER'] / H60 * g('VERD')
                      - p['CLVER'] / H60 / p['VCVER'] * g('VERC'))
    return d


REST = lambda t: (0.0, 0.0, 0.0)


def y_init(p):
    y = np.zeros(NS)
    y[IDX['HR']] = 65; y[IDX['CAI']] = 0.377; y[IDX['CANSR']] = 1180
    y[IDX['CAJSRT']] = 10000; y[IDX['CSQ']] = p['CSQSET']
    y[IDX['NAAV']] = 1.0; y[IDX['CAMT']] = 100; y[IDX['QRSD']] = 95
    return y


_REFS = None


def reference_phospho(minutes=800.0):
    """Basal phospho-levels of a NORMAL, drug-free, un-denervated heart.

    These are FIXED CONSTANTS of the model, not per-scenario quantities.
    Recomputing them for each parameter set (as v2 did) is a real defect: an
    intervention that lowers resting sympathetic tone -- LCSD, a beta-blocker
    present at settle time -- also lowers its own reference, so the model then
    measures delta-PKA from a lowered floor and the intervention comes out
    PROARRHYTHMIC.  In v2 this made LCSD raise peak PVC rate from 32 to 47/min.
    """
    global _REFS
    if _REFS is None:
        p = dict(P)
        p.update(GRYR=0.0, GCSQI=0.0, FLCSD=0.0, CSQSET=1.0, CFLEfree=0.0,
                 PICA0=0.0, PPLB0=0.0, PRYR0=0.0)
        y = y_init(p)
        y = solve_ivp(rhs, (0, minutes), y, args=(p, REST), method='LSODA',
                      rtol=1e-8, atol=1e-10).y[:, -1]
        _REFS = (y[IDX['PICA']], y[IDX['PPLB']], y[IDX['PRYR']])
    return _REFS


def settle(p, minutes=800.0):
    """Settle to the resting state of THIS parameter set, with the phospho
    reference points held at normal-physiology values."""
    p = dict(p)
    p['PICA0'], p['PPLB0'], p['PRYR0'] = reference_phospho()
    y = y_init(p)
    y = solve_ivp(rhs, (0, minutes), y, args=(p, REST), method='LSODA',
                  rtol=1e-8, atol=1e-10).y[:, -1]
    for k in ('NPVC', 'NVT', 'HZRD', 'NSHOCK', 'SYNCH', 'CAINT', 'EXPO'):
        y[IDX[k]] = 0.0
    return y, p


def run(p, y0, protocol, tmax, npts=1500, doses=()):
    p = dict(p)
    times = sorted(set([0.0] + [dd[0] for dd in doses] + [tmax]))
    y = y0.copy(); T, Y = [], []

    def f(t, yy, pp, dr):
        pp = dict(pp)
        pp['CFLEfree'] = (max(yy[IDX['FLEC']], 0) / pp['VCFLE'] * 1e3
                          * pp['FUFLE'] / pp['MWFLE'])
        return rhs(t, yy, pp, dr)

    for a, b in zip(times[:-1], times[1:]):
        for dt_, nm, amt in doses:
            if abs(dt_ - a) < 1e-9:
                y[IDX[nm]] += amt
        n = max(25, int(npts * (b - a) / tmax))
        s = solve_ivp(f, (a, b), y, args=(p, protocol),
                      t_eval=np.linspace(a, b, n), method='LSODA',
                      rtol=1e-7, atol=1e-9)
        T.append(s.t); Y.append(s.y); y = s.y[:, -1]
    return np.concatenate(T), np.concatenate(Y, axis=1), p


def exercise(t0=10.0, ramp=12.0, rec=12.0):
    def dr(t):
        if t < t0:
            return (0.0, 0.0, 0.0)
        if t < t0 + ramp:
            return ((t - t0) / ramp, 0.0, 0.0)
        return (max(0.0, 1.0 - (t - t0 - ramp) / rec), 0.0, 0.0)
    return dr


def genotype(p, g):
    p = dict(p)
    if g == 'CPVT1':
        p['GRYR'] = 0.30; p['FREC'] = 0.60
    elif g == 'CPVT1s':          # severe / homozygous-like
        p['GRYR'] = 0.40; p['FREC'] = 0.50
    elif g == 'CPVT2':
        p['CSQSET'] = 0.45; p['GCSQI'] = 0.28; p['FREC'] = 0.45
    return p


def metrics(T, Y, p):
    hr, pv, vt = Y[IDX['HR']], Y[IDX['PVCR']], Y[IDX['VTB']]
    on = np.where(pv > 1.0)[0]
    onv = np.where(vt > 0.5)[0]
    ve = 0
    if pv.max() > 0.5: ve = 1
    if pv.max() > 5: ve = 2
    if pv.max() > 12: ve = 3
    if vt.max() > 0.5: ve = 4
    return dict(HRpk=hr.max(), PVCpk=pv.max(), VTpk=vt.max(),
                HRon=hr[on[0]] if len(on) else np.nan,
                HRvt=hr[onv[0]] if len(onv) else np.nan,
                NPVC=Y[IDX['NPVC']][-1], VE=ve,
                minVT=float(np.trapezoid(vt, T)))


def show(tag, T, Y, p):
    m = metrics(T, Y, p)
    print(f"{tag:44s} HRpk={m['HRpk']:5.0f} PVC={m['PVCpk']:7.2f}/min "
          f"VT={m['VTpk']:5.3f} HR@PVC={m['HRon']:5.0f} HR@VT={m['HRvt']:5.0f} "
          f"nPVC={m['NPVC']:7.0f} VE={m['VE']} VTmin={m['minVT']:5.2f}")
    return m


if __name__ == '__main__':
    for gt in ('WT', 'CPVT1', 'CPVT1s', 'CPVT2'):
        p = genotype(P, gt)
        y0, p = settle(p)
        F0 = csq_free(y0[IDX['CAJSRT']], p['BCSQ'] * y0[IDX['CSQ']], p['KCSQ'])
        CT0 = p['CTH0'] * (1 - p['GRYR']) * (1 - p['GCSQI'])
        print(f"--- {gt:7s} rest: HR={y0[IDX['HR']]:5.1f} CAI={y0[IDX['CAI']]:.3f} "
              f"NSR={y0[IDX['CANSR']]:6.0f} JSRt={y0[IDX['CAJSRT']]:6.0f} "
              f"JSRf={F0:6.0f} Cth={CT0:6.0f} X={F0/CT0:5.3f} "
              f"PKA={y0[IDX['PKA']]:.3f} PVC={y0[IDX['PVCR']]:.3f}")
        T, Y, p = run(p, y0, exercise(), 40.0)
        show(gt + ' exercise test', T, Y, p)
        for tt in (10, 13, 15, 17, 19, 22):
            i = int(np.argmin(abs(T - tt)))
            F = csq_free(Y[IDX['CAJSRT']][i], p['BCSQ'] * Y[IDX['CSQ']][i], p['KCSQ'])
            CT = (p['CTH0'] * (1 - p['GRYR']) * (1 - p['GCSQI'])
                  / (1 + p['BPRYR'] * (Y[IDX['PRYR']][i] - p['PRYR0'])
                     + p['BCAMK'] * Y[IDX['CAMK']][i]))
            print(f"      t={tt:3.0f} HR={Y[IDX['HR']][i]:6.1f} CAI={Y[IDX['CAI']][i]:.3f}"
                  f" JSRf={F:6.0f} Cth={CT:6.0f} X={F/CT:5.3f} "
                  f"VDAD={Y[IDX['VDAD']][i]:5.1f} PVC={Y[IDX['PVCR']][i]:7.2f} "
                  f"VT={Y[IDX['VTB']][i]:5.3f}")
