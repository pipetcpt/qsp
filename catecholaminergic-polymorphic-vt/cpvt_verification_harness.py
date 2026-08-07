#!/usr/bin/env python3
"""Scenario suite / verification harness for the CPVT QSP model."""
import numpy as np
from cpvt_verification_model import P, IDX, NS, rhs, settle, run, exercise, csq_free, metrics
from cpvt_verification_pk import pk_ss

ROWS = []


def scenario(name, gt=(('GRYR', 0.24), ('FREC', 0.60)), drugs=(), offset_h=24.0,
             pmod=None, protocol=None, tmax=40.0):
    """drugs = ((drug, dose_mg, tau_h), ...) at steady state; the test starts
    `offset_h` after the last dose."""
    p = dict(P)
    for k, v in gt:
        p[k] = v
    if pmod:
        p.update(pmod)
    y0, p = settle(p)
    conc = {}
    for spec in drugs:
        dg, dose, tau_h = spec[0], spec[1], spec[2]
        # Each drug is read at ITS OWN dosing-interval position: asking for
        # "24 h after the nadolol dose" means a bid drug is 12 h after its own
        # last dose, not 24 h (which would be two missed doses).  A 4th element
        # overrides this with an explicit offset -- that is how a genuinely
        # MISSED dose is expressed.
        off = spec[3] if len(spec) > 3 else (
            offset_h if offset_h <= tau_h else (offset_h % tau_h or tau_h))
        amts, c = pk_ss(dg, dose, tau_h, off, p)
        for nm, a in amts.items():
            y0[IDX[nm]] = a
        conc[dg] = c
    prot = protocol if protocol else exercise(t0=10, ramp=12, rec=12)
    T, Y, p = run(p, y0, prot, tmax, npts=1600)
    m = metrics(T, Y, p)
    m['name'] = name
    m['conc'] = conc
    m['QRS'] = Y[IDX['QRSD']].max()
    m['NAAV'] = Y[IDX['NAAV']].min()
    m['HZ'] = Y[IDX['HZRD']][-1]
    m['SHK'] = Y[IDX['NSHOCK']][-1]
    m['EPIpk'] = Y[IDX['EPI']].max()
    ROWS.append(m)
    cs = ' '.join(f"{k}={v:6.1f}" for k, v in conc.items())
    print(f"{name:47s} HRpk={m['HRpk']:5.0f} PVC={m['PVCpk']:7.2f} "
          f"VT={m['VTpk']:5.3f} VTmin={m['minVT']:5.2f} nPVC={m['NPVC']:6.0f} "
          f"VE={m['VE']} QRS={m['QRS']:5.0f} {cs}")
    return m, T, Y, p


NORM = (('GRYR', 0.0), ('FREC', 1.0))
CP1 = (('GRYR', 0.24), ('FREC', 0.60))
CP1S = (('GRYR', 0.30), ('FREC', 0.55))
CP1M = (('GRYR', 0.18), ('FREC', 0.70))
CP2 = (('GRYR', 0.0), ('FREC', 0.45))
CP2P = dict(CSQSET=0.45, GCSQI=0.28)

