"""Reference (non-mrgsolve) transcription of the CRPS QSP model.

Purpose: every quantitative claim in README.md and in the header of
crps_mrgsolve_model.R was produced by THIS file, so the numbers can be
re-derived with nothing but numpy + scipy (no R, no mrgsolve, no compiler).
The equations, parameter values and dosing regimens are the same as in
crps_mrgsolve_model.R; if the two ever disagree, the mrgsolve file is the
model of record and this file is the thing to fix.

    python3 crps_reference_check.py            # run every analysis
    python3 crps_reference_check.py nat win    # run selected analyses

Analyses: nat (natural history), kf (fear-avoidance trait bifurcation),
inj (insult-magnitude scan), win (therapeutic window), arm (arm
decomposition), ket (ketamine washout), dt (dose vs timing), ph (arm
ordering by phenotype), bone (neridronate), scs (SCS habituation).
"""
import numpy as np
from scipy.integrate import solve_ivp

P = dict(
    # ketamine PK
    KET_V1=22., KET_V2=200., KET_Q=60., KET_CL=85., KET_FM=0.80,
    NORKET_V=70., NORKET_CL=40., NORKET_POT=0.30,
    PRED_KA=2.0, PRED_F=0.80, PRED_V=40., PRED_CL=9.0,
    NER_V=20., NER_CLR=6.0, NER_KBONE=4.0, NER_KOFF=2e-5,
    GBP_VMAX=45., GBP_KM=420., GBP_F=0.90, GBP_V=58., GBP_CL=7.5,
    AMT_KA=1.0, AMT_F=0.50, AMT_V=900., AMT_CL=45.,
    NAC_KA=1.2, NAC_F=0.10, NAC_V=35., NAC_CL=15.,
    IVIG_V=5.0, IVIG_KEL=0.0014,
    # forcing
    INJ_AMP=1.0, INJ_TAU=400.0, INJ_T0=0., IMMOB_DUR=720., IMMOB_AMP=0.25,
    # fast node
    KIN_NP=0.02, KOUT_NP=0.05, NP_SENS_FB=0.3,
    KIN_CYT=0.02, KOUT_CYT=0.045, W_NP_CYT=1.2,
    KIN_NGF=0.01, KOUT_NGF=0.015,
    KIN_EDEMA=0.02, KOUT_EDEMA=0.03,
    KIN_ROS=0.05, KOUT_ROS=0.09, W_HYP_ROS=0.70,
    KIN_AAB=0.0015, KOUT_AAB=0.0018, W_AAB_PS=0.30,
    KIN_A1=0.0025, KOUT_A1=0.006, W_NGF_A1=0.80, SYMP_TONE=1.0,
    # vaso
    PERF_BASE=1.0, KPERF=0.03, W_NP_PERF=3.0, W_A1_PERF=1.4,
    W_ROS_PERF=0.40, W_SS_PERF=0.5, KIN_HYP=0.05, KOUT_HYP=0.04, W_DIS_HYP=0.2,
    TEMP_GAIN=3.0,
    # periph sens
    KIN_PS=0.04, KOUT_PS=0.07, W_NP_PS=0.35, W_CYT_PS=0.40, W_ROS_PS=0.30,
    W_A1_PS=0.45, W_HYP_PS=0.35, W_BONE_PS=0.30,
    # central ring
    KIN_SS=0.012, KOUT_SS=0.01, W_PS_SS=1.6,
    KIN_GLIA=0.0016, KOUT_GLIA=0.0008, W_GLIA_SS=0.4,
    W_GLIA_SELF=1.3, GLIA50=0.36, HILL_GLIA=16.0,
    DINH_BASE=1.0, KDINH=0.004, W_SS_DINH=0.55,
    KIN_CTX=0.0025, KOUT_CTX=0.0015, W_SS_CTX=0.45, W_DIS_CTX=0.9,
    KIN_DIS=0.006, KOUT_DIS=0.004, KFEAR=1.0, W_ROM_DIS=0.15,
    PAIN50_DRIVE=2.5, PAIN50_FEAR=5.0, HILL_FEAR=3.0, DIS50_CTX=0.45, HILL_CTX=3.0,
    # pain / rom
    PAIN_MAX=10., KPAIN=0.25, W_PS_PAIN=2.2, W_SS_PAIN=3.6,
    W_CTX_PAIN=3.2, W_DIS_PAIN=1.6,
    KROM=0.010, ROM_MAX=1.0, W_EDEMA_ROM=0.35, W_PAIN_ROM=0.055,
    W_DIS_ROM=0.45,
    # bone
    KIN_OC=0.008, OC_BASE=1.0, W_NP_OC=0.70, W_CYT_OC=0.50, W_DIS_OC=0.35,
    KBMD=0.00016, BMD_REC=2.5, BMD_BASE=1.0, KCTX=0.05, CTX_BASE=0.35, CTX_GAIN=0.30,
    # PD
    EC50_KET=130., EMAX_KET=0.85, KET_GLIA_E=0.30,
    EC50_PRED=25., EMAX_PRED=0.75, PRED_EDEMA=0.55,
    EMAX_NER=0.80, BONE50_NER=25.,
    EC50_GBP=4000., EMAX_GBP=0.35,
    EC50_AMT=60., EMAX_AMT=0.60,
    EC50_NAC=900., EMAX_NAC=0.55,
    EC50_IVIG=6.0, EMAX_IVIG=0.65,
    VASODIL=0., VASODIL_E=0.45,
    A_KET=0.75, A_GBP=0.30, A_AMT=0.25, A_SCS=0.80, A_MAX=0.85,
    # interventions
    REHAB=0., REHAB_T0=0., REHAB_DUR=4380., REHAB_CTX=0.85, REHAB_DIS=0.90,
    REHAB_FEAR=0.60,
    SCS_ON=0., SCS_T0=1e9, SCS_EFF=0.55, SCS_HAB_K=8e-5,
    SYMPBLOCK=0., SB_T0=1e9, SB_EFF=0.60, SB_TAU=336.,
)

NAMES = ['KET_C1','KET_C2','NORKET','PRED_GUT','PRED_C','NER_C','NER_BONE',
         'GBP_GUT','GBP_C','AMT_GUT','AMT_C','NAC_GUT','NAC_C','IVIG_C',
         'NP','CYT','NGF','EDEMA','ROS','AAB','ALPHA1','PSENS','PERF','HYPOX',
         'SSENS','GLIA','DINH','CORTEX','DISUSE','PAIN','ROM','OC','BMD','CTXI']
IDX = {n: i for i, n in enumerate(NAMES)}


def y0(p):
    y = np.zeros(len(NAMES))
    y[IDX['PERF']] = p['PERF_BASE']
    y[IDX['DINH']] = p['DINH_BASE']
    y[IDX['ROM']] = p['ROM_MAX']
    y[IDX['OC']] = p['OC_BASE']
    y[IDX['BMD']] = p['BMD_BASE']
    y[IDX['CTXI']] = p['CTX_BASE']
    return y


def hill(x, x50, n):
    x = max(x, 0.0)
    return x**n / (x50**n + x**n)


def rhs(t, y, p, infusions):
    (KET_C1, KET_C2, NORKET, PRED_GUT, PRED_C, NER_C, NER_BONE, GBP_GUT,
     GBP_C, AMT_GUT, AMT_C, NAC_GUT, NAC_C, IVIG_C, NP, CYT, NGF, EDEMA,
     ROS, AAB, ALPHA1, PSENS, PERF, HYPOX, SSENS, GLIA, DINH, CORTEX,
     DISUSE, PAIN, ROM, OC, BMD, CTXI) = y
    d = np.zeros_like(y)

    ket_c = 1000.*KET_C1/p['KET_V1']
    nork_c = 1000.*NORKET/p['NORKET_V']
    ket_eff = ket_c + p['NORKET_POT']*nork_c
    pred_c = 1000.*PRED_C/p['PRED_V']
    gbp_c = 1000.*GBP_C/p['GBP_V']
    amt_c = 1000.*AMT_C/p['AMT_V']
    nac_c = 1000.*NAC_C/p['NAC_V']
    ivig_c = IVIG_C/p['IVIG_V']

    E_ket = p['EMAX_KET']*ket_eff/(p['EC50_KET']+ket_eff)
    E_ketgl = p['KET_GLIA_E']*ket_eff/(p['EC50_KET']+ket_eff)
    E_pred = p['EMAX_PRED']*pred_c/(p['EC50_PRED']+pred_c)
    E_predE = p['PRED_EDEMA']*pred_c/(p['EC50_PRED']+pred_c)
    E_gbp = p['EMAX_GBP']*gbp_c/(p['EC50_GBP']+gbp_c)
    E_amt = p['EMAX_AMT']*amt_c/(p['EC50_AMT']+amt_c)
    E_nac = p['EMAX_NAC']*nac_c/(p['EC50_NAC']+nac_c)
    E_ivig = p['EMAX_IVIG']*ivig_c/(p['EC50_IVIG']+ivig_c)
    E_ner = p['EMAX_NER']*NER_BONE/(p['BONE50_NER']+NER_BONE)

    t_inj = t - p['INJ_T0']
    INJ = p['INJ_AMP']*np.exp(-t_inj/p['INJ_TAU']) if t_inj >= 0 else 0.
    IMMOB = p['IMMOB_AMP'] if (p['INJ_T0'] <= t < p['INJ_T0']+p['IMMOB_DUR']) else 0.
    rehab = p['REHAB'] if (p['REHAB'] > 0 and p['REHAB_T0'] <= t
                           < p['REHAB_T0']+p['REHAB_DUR']) else 0.
    scs_t = t - p['SCS_T0']
    scs = p['SCS_EFF']*np.exp(-p['SCS_HAB_K']*scs_t) if (p['SCS_ON'] > 0 and scs_t >= 0) else 0.
    sb_t = t - p['SB_T0']
    sb = (p['SB_EFF']*np.exp(-sb_t/(p['SB_TAU']/3.))
          if (p['SYMPBLOCK'] > 0 and 0 <= sb_t < p['SB_TAU']) else 0.)

    # ---- PK ----
    rate_ket = infusions.get('KET_C1', lambda tt: 0.)(t)
    d[IDX['KET_C1']] = (rate_ket + p['KET_Q']*KET_C2/p['KET_V2']
                        - (p['KET_Q']+p['KET_CL'])*KET_C1/p['KET_V1'])
    d[IDX['KET_C2']] = p['KET_Q']*(KET_C1/p['KET_V1'] - KET_C2/p['KET_V2'])
    d[IDX['NORKET']] = (p['KET_FM']*p['KET_CL']*KET_C1/p['KET_V1']
                        - p['NORKET_CL']*NORKET/p['NORKET_V'])
    d[IDX['PRED_GUT']] = -p['PRED_KA']*PRED_GUT
    d[IDX['PRED_C']] = (p['PRED_F']*p['PRED_KA']*PRED_GUT
                        - p['PRED_CL']*PRED_C/p['PRED_V'])
    d[IDX['NER_C']] = -(p['NER_CLR']+p['NER_KBONE'])*NER_C/p['NER_V'] + p['NER_KOFF']*NER_BONE
    d[IDX['NER_BONE']] = p['NER_KBONE']*NER_C/p['NER_V'] - p['NER_KOFF']*NER_BONE
    abs_gbp = p['GBP_VMAX']*GBP_GUT/(p['GBP_KM']+GBP_GUT)
    d[IDX['GBP_GUT']] = -abs_gbp
    d[IDX['GBP_C']] = p['GBP_F']*abs_gbp - p['GBP_CL']*GBP_C/p['GBP_V']
    d[IDX['AMT_GUT']] = -p['AMT_KA']*AMT_GUT
    d[IDX['AMT_C']] = p['AMT_F']*p['AMT_KA']*AMT_GUT - p['AMT_CL']*AMT_C/p['AMT_V']
    d[IDX['NAC_GUT']] = -p['NAC_KA']*NAC_GUT
    d[IDX['NAC_C']] = p['NAC_F']*p['NAC_KA']*NAC_GUT - p['NAC_CL']*NAC_C/p['NAC_V']
    d[IDX['IVIG_C']] = -p['IVIG_KEL']*IVIG_C

    # ---- fast peripheral node ----
    d[IDX['NP']] = p['KIN_NP']*(INJ + p['NP_SENS_FB']*PSENS)*(1-NP) - p['KOUT_NP']*NP
    d[IDX['CYT']] = p['KIN_CYT']*(0.6*INJ + p['W_NP_CYT']*NP)*(1-E_pred)*(1-CYT) - p['KOUT_CYT']*CYT
    d[IDX['NGF']] = p['KIN_NGF']*(0.5*INJ + 0.8*CYT)*(1-NGF) - p['KOUT_NGF']*NGF
    d[IDX['EDEMA']] = (p['KIN_EDEMA']*(NP+0.5*CYT)*(1-E_predE)*(1-EDEMA)
                       - p['KOUT_EDEMA']*EDEMA*(1+0.5*rehab))
    d[IDX['ROS']] = (p['KIN_ROS']*(CYT + p['W_HYP_ROS']*HYPOX)*(1-ROS)
                     - p['KOUT_ROS']*ROS*(1+3*E_nac))
    d[IDX['AAB']] = p['KIN_AAB']*(INJ+0.4*CYT)*(1-AAB) - p['KOUT_AAB']*AAB
    d[IDX['ALPHA1']] = (p['KIN_A1']*(p['W_NGF_A1']*NGF+0.3*INJ)*p['SYMP_TONE']*(1-sb)*(1-ALPHA1)
                        - p['KOUT_A1']*ALPHA1)
    aab_drive = p['W_AAB_PS']*AAB*(1-E_ivig)
    bone_drive = p['W_BONE_PS']*max(OC-p['OC_BASE'], 0.)
    ps_in = (p['W_NP_PS']*(NP+0.6*NGF) + p['W_CYT_PS']*CYT + p['W_ROS_PS']*ROS
             + p['W_A1_PS']*ALPHA1*p['SYMP_TONE']*(1-sb) + p['W_HYP_PS']*HYPOX
             + aab_drive + bone_drive)
    d[IDX['PSENS']] = p['KIN_PS']*ps_in*(1-PSENS) - p['KOUT_PS']*PSENS

    # ---- vasomotor ----
    perf_target = (p['PERF_BASE'] + p['W_NP_PERF']*NP
                   - (p['W_A1_PERF']*ALPHA1*p['SYMP_TONE'] + p['W_ROS_PERF']*ROS
                      + p['W_SS_PERF']*SSENS*p['SYMP_TONE'])
                   * (1 - (p['VASODIL_E'] if p['VASODIL'] > 0 else 0.)))
    d[IDX['PERF']] = p['KPERF']*(perf_target-PERF)
    perf_def = max(p['PERF_BASE']-PERF, 0.)
    d[IDX['HYPOX']] = p['KIN_HYP']*(perf_def + p['W_DIS_HYP']*DISUSE)*(1-HYPOX) - p['KOUT_HYP']*HYPOX

    # ---- central ring ----
    spinal_gain = ((p['W_PS_SS']*PSENS*(1+p['W_GLIA_SS']*GLIA)
                    + p['W_GLIA_SELF']*hill(GLIA, p['GLIA50'], p['HILL_GLIA']))
                   * (1-E_ket)*(1-E_gbp)*(1-scs))
    d[IDX['SSENS']] = p['KIN_SS']*spinal_gain*(1-SSENS) - p['KOUT_SS']*SSENS*DINH
    d[IDX['GLIA']] = p['KIN_GLIA']*SSENS*(1-E_ketgl)*(1-GLIA) - p['KOUT_GLIA']*GLIA
    dinh_target = p['DINH_BASE']*(1 - p['W_SS_DINH']*SSENS/(1+SSENS))*(1+E_amt)
    d[IDX['DINH']] = p['KDINH']*(dinh_target-DINH)
    d[IDX['CORTEX']] = (p['KIN_CTX']*(p['W_SS_CTX']*SSENS
                        + p['W_DIS_CTX']*hill(DISUSE, p['DIS50_CTX'], p['HILL_CTX']))*(1-CORTEX)
                        - (p['KOUT_CTX'] + p['REHAB_CTX']*p['KIN_CTX']*rehab*8.)*CORTEX)
    fear_gain = p['KFEAR']*(1 - p['REHAB_FEAR']*rehab)
    dis_in = (fear_gain*hill(PAIN, p['PAIN50_FEAR'], p['HILL_FEAR'])
              + p['W_ROM_DIS']*(p['ROM_MAX']-ROM) + IMMOB)
    d[IDX['DISUSE']] = (p['KIN_DIS']*dis_in*(1-DISUSE)
                        - (p['KOUT_DIS'] + p['REHAB_DIS']*p['KIN_DIS']*rehab*6.)*DISUSE)

    pain_drive = (p['W_PS_PAIN']*PSENS + p['W_SS_PAIN']*SSENS
                  + p['W_CTX_PAIN']*CORTEX + p['W_DIS_PAIN']*DISUSE)
    acute = min(p['A_KET']*E_ket + p['A_GBP']*E_gbp + p['A_AMT']*E_amt
                + p['A_SCS']*scs, p['A_MAX'])
    pain_target = (p['PAIN_MAX']*pain_drive/(p['PAIN50_DRIVE']+pain_drive))*(1-acute)
    d[IDX['PAIN']] = p['KPAIN']*(pain_target-PAIN)
    rom_target = p['ROM_MAX']*(1 - p['W_EDEMA_ROM']*EDEMA
                               - p['W_PAIN_ROM']*PAIN
                               - p['W_DIS_ROM']*DISUSE)
    rom_target = max(rom_target, 0.05)
    d[IDX['ROM']] = p['KROM']*(rom_target-ROM)*(1+2*rehab)

    oc_target = (p['OC_BASE']*(1 + p['W_NP_OC']*NP + p['W_CYT_OC']*CYT
                               + p['W_DIS_OC']*DISUSE)*(1-E_ner))
    d[IDX['OC']] = p['KIN_OC']*(oc_target-OC)
    d[IDX['CTXI']] = p['KCTX']*(p['CTX_BASE'] + p['CTX_GAIN']*(OC-p['OC_BASE']) - CTXI)
    d[IDX['BMD']] = (-p['KBMD']*max(OC-p['OC_BASE'], 0.)*BMD
                     + p['KBMD']*p['BMD_REC']*(1-DISUSE)*(p['BMD_BASE']-BMD))
    return d


# ---------------------------------------------------------------- dosing
def bolus_times(start, ii, addl):
    return [start + i*ii for i in range(addl+1)]


def run(pars=None, doses=None, ket=None, end=2*8760, npts=None):
    """doses: list of (cmt, amt, times). ket: (start, hours, rate_max)."""
    p = dict(P)
    if pars:
        p.update(pars)
    infus = {}
    if ket:
        start, hours, rmax = ket
        steps = [0.25, 0.5, 0.75, 1.0]
        dur = hours/len(steps)
        def kr(t, start=start, dur=dur, steps=steps, rmax=rmax):
            if t < start or t >= start+hours:
                return 0.
            i = int((t-start)//dur)
            i = min(i, len(steps)-1)
            return rmax*steps[i]
        infus['KET_C1'] = kr
    events = []
    for (cmt, amt, times) in (doses or []):
        for tt in times:
            events.append((tt, cmt, amt))
    if ket:
        events += [(ket[0], None, None), (ket[0]+ket[1], None, None)]
    events.sort(key=lambda e: e[0])
    breaks = sorted(set([e[0] for e in events] + [0., p['IMMOB_DUR'],
                    p['REHAB_T0'], p['REHAB_T0']+p['REHAB_DUR'],
                    p['SB_T0'], p['SB_T0']+p['SB_TAU'], p['SCS_T0'], end]))
    breaks = [b for b in breaks if 0 <= b <= end]
    y = y0(p)
    ts, ys = [0.], [y.copy()]
    cur = 0.
    for b in breaks + [end]:
        # apply events at time cur
        for (tt, cmt, amt) in events:
            if cmt is not None and abs(tt-cur) < 1e-9:
                y[IDX[cmt]] += amt
        if b <= cur:
            continue
        n = max(int((b-cur)/6)+2, 3)
        teval = np.linspace(cur, b, n)
        sol = solve_ivp(rhs, (cur, b), y, args=(p, infus), method='LSODA',
                        t_eval=teval, rtol=1e-6, atol=1e-9, max_step=6.)
        ts.extend(sol.t[1:]); ys.extend(sol.y[:, 1:].T)
        y = sol.y[:, -1].copy()
        cur = b
    T = np.array(ts)
    Y = np.array(ys)
    return T, Y, p


def get(T, Y, name, day=None):
    col = Y[:, IDX[name]]
    if day is None:
        return col[-1]
    i = int(np.argmin(np.abs(T-day*24)))
    return col[i]


def css(T, Y, p, day=None):
    pain = get(T, Y, 'PAIN', day)
    perf = get(T, Y, 'PERF', day)
    ed = get(T, Y, 'EDEMA', day)
    rom = get(T, Y, 'ROM', day)
    ta = p['TEMP_GAIN']*(perf-p['PERF_BASE'])
    c = (4*pain/p['PAIN_MAX'] + 4*min(abs(ta)/3., 1.)
         + 4*ed/(0.8+ed) + 4*(p['ROM_MAX']-rom)/p['ROM_MAX'])
    return min(c, 16.)


# ---------------------------------------------------------------- regimens
def pred(start):
    return [('PRED_GUT', 40., bolus_times(start, 24, 13)),
            ('PRED_GUT', 20., bolus_times(start+14*24, 24, 6)),
            ('PRED_GUT', 10., bolus_times(start+21*24, 24, 6))]


def nac(start, days=180):
    return [('NAC_GUT', 600., bolus_times(start, 8, int(days*3)-1))]


def ner(start):
    return [('NER_C', 100., bolus_times(start, 72, 3))]


def gbp(start, days=360):
    return [('GBP_GUT', 600., bolus_times(start, 8, int(days*3)-1))]


def amt(start, days=360):
    return [('AMT_GUT', 50., bolus_times(start, 24, days-1))]


if __name__ == '__main__':
    import sys
    T, Y, p = run(end=3*8760)
    for d in [7, 30, 90, 180, 365, 730, 1095]:
        print(f"day {d:5d}  NRS {get(T,Y,'PAIN',d):5.2f}  CSS {css(T,Y,p,d):5.2f} "
              f"PSENS {get(T,Y,'PSENS',d):5.2f} SSENS {get(T,Y,'SSENS',d):5.2f} "
              f"GLIA {get(T,Y,'GLIA',d):5.2f} CTX {get(T,Y,'CORTEX',d):5.2f} "
              f"DIS {get(T,Y,'DISUSE',d):5.2f} ROM {get(T,Y,'ROM',d):5.2f} "
              f"TA {3*(get(T,Y,'PERF',d)-1):5.2f} BMD {get(T,Y,'BMD',d):5.3f}")


## ==========================================================================
## Verification battery
## ==========================================================================

B = dict(P)   # calibrated parameters are the module defaults
D = 24


def R(extra=None, doses=None, ket=None, end=3*8760):
    p = dict(B)
    if extra:
        p.update(extra)
    return run(pars=p, doses=doses, ket=ket, end=end)


def pkg(start, arms=('pred', 'nac', 'rehab')):
    """literature multimodal package starting at `start` hours"""
    ex, ds = {}, []
    if 'rehab' in arms:
        ex.update(REHAB=1., REHAB_T0=float(start), REHAB_DUR=180*D)
    if 'pred' in arms:
        ds += pred(start)
    if 'nac' in arms:
        ds += nac(start, 180)
    return ex, ds


def head(t):
    print("\n" + "="*78 + f"\n{t}\n" + "="*78)


def natural():
    head("1. NATURAL HISTORY (distal radius fracture, 30-day cast, vulnerable trait)")
    T, Y, p = R()
    print(f"{'day':>5} {'NRS':>5} {'CSS':>5} {'TA_C':>6} {'PSENS':>6} {'SSENS':>6} {'GLIA':>6} "
          f"{'CORTEX':>7} {'DISUSE':>7} {'ROM':>5} {'BMD':>5} {'CTXI':>5} {'EDEMA':>6}")
    for d in [3, 7, 14, 21, 30, 45, 60, 90, 180, 365, 730, 1095]:
        print(f"{d:5d} {get(T,Y,'PAIN',d):5.2f} {css(T,Y,p,d):5.2f} {3*(get(T,Y,'PERF',d)-1):+6.2f} "
              f"{get(T,Y,'PSENS',d):6.3f} {get(T,Y,'SSENS',d):6.3f} {get(T,Y,'GLIA',d):6.3f} "
              f"{get(T,Y,'CORTEX',d):7.3f} {get(T,Y,'DISUSE',d):7.3f} {get(T,Y,'ROM',d):5.2f} "
              f"{get(T,Y,'BMD',d):5.3f} {get(T,Y,'CTXI',d):5.2f} {get(T,Y,'EDEMA',d):6.3f}")


def _kf(k):
    T, Y, p = R(dict(KFEAR=k))
    return (k, round(float(get(T, Y, 'PAIN')), 2), round(float(get(T, Y, 'GLIA')), 3),
            round(float(get(T, Y, 'CORTEX')), 3), round(float(get(T, Y, 'DISUSE')), 3),
            round(css(T, Y, p), 2), round(float(get(T, Y, 'BMD'))), )


def kfear_scan():
    head("2. BIFURCATION IN THE FEAR-AVOIDANCE TRAIT (KFEAR), 3-year endpoint")
    ks = [0.2, 0.35, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0, 1.2, 1.4]
    res = list(map(_kf, ks))
    print(f"{'KFEAR':>6} {'NRS':>6} {'GLIA':>6} {'CORTEX':>7} {'DISUSE':>7} {'CSS':>6}")
    for k, n, g, c, di, cs, _ in res:
        print(f"{k:6.2f} {n:6.2f} {g:6.3f} {c:7.3f} {di:7.3f} {cs:6.2f}")
    return res


def _inj(a):
    T, Y, p = R(dict(INJ_AMP=a))
    return (a, round(float(get(T, Y, 'PAIN')), 2), round(float(get(T, Y, 'GLIA')), 3),
            round(float(max(Y[:, IDX['PAIN']])), 2), round(css(T, Y, p), 2))


def inj_scan():
    head("3. INSULT-MAGNITUDE SCAN (three attractors), 3-year endpoint")
    aa = [0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 0.9, 0.95, 1.0, 1.2, 1.5]
    res = list(map(_inj, aa))
    print(f"{'INJ_AMP':>8} {'NRS_3y':>7} {'GLIA':>6} {'peakNRS':>8} {'CSS':>6}")
    for a, n, g, pk, cs in res:
        print(f"{a:8.2f} {n:7.2f} {g:6.3f} {pk:8.2f} {cs:6.2f}")
    return res


def _win(d):
    ex, ds = pkg(d*D)
    T, Y, p = R(ex, doses=ds)
    return (d, round(float(get(T, Y, 'PAIN')), 2), round(css(T, Y, p), 2),
            round(float(get(T, Y, 'GLIA')), 3), round(float(get(T, Y, 'CORTEX')), 3),
            round(float(get(T, Y, 'BMD')), 3), round(float(get(T, Y, 'ROM')), 3))


def window():
    head("4. THERAPEUTIC WINDOW: identical package (prednisolone taper + NAC + rehab)")
    ds = [3, 7, 14, 21, 30, 45, 60, 75, 90, 120, 180, 240, 365]
    res = list(map(_win, ds))
    print(f"{'start_d':>8} {'NRS_3y':>7} {'CSS':>6} {'GLIA':>6} {'CORTEX':>7} {'BMD':>6} {'ROM':>6}")
    for d, n, c, g, cx, b, r in res:
        print(f"{d:8d} {n:7.2f} {c:6.2f} {g:6.3f} {cx:7.3f} {b:6.3f} {r:6.3f}")
    return res


def _armwin(args):
    d, arms = args
    ex, ds = pkg(d*D, arms)
    T, Y, p = R(ex, doses=ds)
    return (d, arms, round(float(get(T, Y, 'PAIN')), 2))


def window_by_arm():
    head("5. WHICH ARM CARRIES THE WINDOW? (single-arm packages)")
    combos = [(d, a) for d in [7, 30, 60, 120, 240]
              for a in [('pred',), ('nac',), ('rehab',), ('pred', 'nac'), ('pred', 'nac', 'rehab')]]
    res = list(map(_armwin, combos))
    arms = [('pred',), ('nac',), ('rehab',), ('pred', 'nac'), ('pred', 'nac', 'rehab')]
    print(f"{'start_d':>8} " + " ".join(f"{'+'.join(a):>18}" for a in arms))
    for d in [7, 30, 60, 120, 240]:
        row = {a: v for (dd, a, v) in res if dd == d}
        print(f"{d:8d} " + " ".join(f"{row[a]:18.2f}" for a in arms))
    return res


def ketamine():
    head("6. KETAMINE 100-h INFUSION: acute analgesia vs disease modification")
    st = 30*D
    runs = {}
    runs['untreated'] = R(end=2*8760)
    runs['ketamine'] = R(ket=(st, 100, 22), end=2*8760)
    ex, _ = pkg(st, ('rehab',))
    runs['rehab'] = R(ex, end=2*8760)
    runs['ket+rehab'] = R(ex, ket=(st, 100, 22), end=2*8760)
    ks = list(runs)
    print(f"{'day':>5} " + " ".join(f"{k:>11s}" for k in ks))
    for d in [29, 31, 32, 34, 37, 44, 58, 79, 107, 210, 395, 700]:
        print(f"{d:5d} " + " ".join(f"{get(*runs[k][:2],'PAIN',d):11.2f}" for k in ks))
    T, Y, p = runs['ketamine']
    print(f"ketamine Cmax {1000*max(Y[:,IDX['KET_C1']])/p['KET_V1']:.0f} ng/mL, "
          f"norketamine Cmax {1000*max(Y[:,IDX['NORKET']])/p['NORKET_V']:.0f} ng/mL")
    v = {k: float(get(*runs[k][:2], 'PAIN', 700)) for k in ks}
    ek = v['untreated']-v['ketamine']; er = v['untreated']-v['rehab']
    ekr = v['untreated']-v['ket+rehab']
    print(f"700-day effects: ketamine {ek:+.2f}, rehab {er:+.2f}, combination {ekr:+.2f}, "
          f"interaction {ekr-(ek+er):+.2f}")
    return runs


def _dt(args):
    rate, delay = args
    T, Y, p = R(ket=(delay*D, 100, rate), end=2*8760)
    return (rate, delay, round(float(get(T, Y, 'PAIN')), 2),
            round(1000*float(max(Y[:, IDX['KET_C1']]))/p['KET_V1'], 0))


def dose_timing():
    head("7. DOSE vs TIMING (ketamine, 2-year NRS) — plus rehab-coupled variant")
    rates = [5.5, 11, 22, 44, 88]
    delays = [14, 30, 60, 120, 240]
    res = list(map(_dt, [(r, d) for r in rates for d in delays]))
    print("ketamine alone:")
    print(f"{'mg/h':>6} " + " ".join(f"{'d'+str(d):>8}" for d in delays) + "   Cmax(ng/mL)")
    for r in rates:
        row = {d: v for (rr, d, v, c) in res if rr == r}
        cmax = [c for (rr, d, v, c) in res if rr == r][0]
        print(f"{r:6.1f} " + " ".join(f"{row[d]:8.2f}" for d in delays) + f"   {cmax:.0f}")
    return res


def _ph(args):
    label, extra, start, arm = args
    ex = dict(extra); ds = None
    if arm == 'steroid': ds = pred(start)
    elif arm == 'antiox': ds = nac(start, 180)
    elif arm == 'vasodil': ex['VASODIL'] = 1.
    elif arm == 'bisphos': ds = ner(start)
    elif arm == 'gabapentin': ds = gbp(start, 360)
    elif arm == 'amitript': ds = amt(start, 360)
    elif arm == 'rehab': ex.update(REHAB=1., REHAB_T0=float(start), REHAB_DUR=180*D)
    elif arm == 'sympblock': ex.update(SYMPBLOCK=1., SB_T0=float(start))
    elif arm == 'ivig': ds = [('IVIG_C', 35., [start + i*28*D for i in range(6)])]
    T, Y, p = R(ex, doses=ds, end=2*8760)
    return (label, arm, round(float(get(T, Y, 'PAIN')), 2), round(css(T, Y, p), 2))


def phenotype():
    head("8. INTERVENTION ORDERING BY PHENOTYPE / TIMING (delta NRS at 2 years)")
    arms = ['none', 'steroid', 'antiox', 'vasodil', 'bisphos', 'gabapentin',
            'amitript', 'rehab', 'sympblock', 'ivig']
    conds = [('warm-early d14 (SYMP 0.8)', dict(SYMP_TONE=0.8), 14*D),
             ('cold-late d240 (SYMP 1.7)', dict(SYMP_TONE=1.7), 240*D)]
    jobs = [(c[0], c[1], c[2], a) for c in conds for a in arms]
    res = list(map(_ph, jobs))
    for label, extra, start in conds:
        base = [v for (l, a, v, cs) in res if l == label and a == 'none'][0]
        print(f"\n{label}: untreated 2-y NRS {base:.2f}")
        rows = [(a, base-v, cs) for (l, a, v, cs) in res if l == label and a != 'none']
        rows.sort(key=lambda r: -r[1])
        for i, (a, dv, cs) in enumerate(rows, 1):
            print(f"   {i:2d}. {a:12s} dNRS {dv:+6.2f}   CSS {cs:5.2f}")
    return res


def bone():
    head("9. BONE AXIS: neridronate 100 mg IV x4 (Varenna-style), early vs late")
    runs = {'untreated': R(end=8760),
            'nerid_d30': R(doses=ner(30*D), end=8760),
            'nerid_d120': R(doses=ner(120*D), end=8760)}
    print(f"{'day':>5} " + " ".join(f"{k+'_NRS':>16s}" for k in runs))
    for d in [30, 40, 70, 120, 160, 250, 365]:
        print(f"{d:5d} " + " ".join(f"{get(*runs[k][:2],'PAIN',d):16.2f}" for k in runs))
    print(f"{'day':>5} {'CTXI_untr':>10} {'CTXI_ner':>9} {'BMD_untr':>9} {'BMD_ner':>8} {'NER_bone_mg':>12}")
    for d in [30, 40, 70, 120, 250, 365]:
        print(f"{d:5d} {get(*runs['untreated'][:2],'CTXI',d):10.3f} "
              f"{get(*runs['nerid_d30'][:2],'CTXI',d):9.3f} "
              f"{get(*runs['untreated'][:2],'BMD',d):9.3f} "
              f"{get(*runs['nerid_d30'][:2],'BMD',d):8.3f} "
              f"{get(*runs['nerid_d30'][:2],'NER_BONE',d):12.1f}")
    return runs


def scs():
    head("10. SCS AT MONTH 6 WITH HABITUATION (Kemler 2000 vs 2008)")
    runs = {'untreated': R(end=5*8760), 'scs': R(dict(SCS_ON=1., SCS_T0=180.*D), end=5*8760)}
    print(f"{'month':>6} {'NRS_untreated':>14} {'NRS_SCS':>9} {'delta':>7}")
    for m in [6, 7, 9, 12, 18, 24, 36, 48, 60]:
        d = m*30
        a = float(get(*runs['untreated'][:2], 'PAIN', d))
        b = float(get(*runs['scs'][:2], 'PAIN', d))
        print(f"{m:6d} {a:14.2f} {b:9.2f} {a-b:+7.2f}")
    return runs



if __name__ == '__main__':
    import sys
    fns = dict(nat=natural, kf=kfear_scan, inj=inj_scan, win=window,
               arm=window_by_arm, ket=ketamine, dt=dose_timing, ph=phenotype,
               bone=bone, scs=scs)
    which = sys.argv[1:] or ['all']
    for name, fn in fns.items():
        if 'all' in which or name in which:
            fn()
