"""Steady-state PK by integrating ONLY the linear PK subsystem (cheap),
so the stiff 44-state physiology model only has to run the 40-min test."""
import numpy as np
from scipy.integrate import solve_ivp

SPEC = {   # drug: (ka, F, VC, VP, Q, CL, (dep, cen, per))
    'nad': ('KANAD', 'FNAD', 'VCNAD', 'VPNAD', 'QNAD', 'CLNAD', ('NADD', 'NADC', 'NADP')),
    'met': ('KAMET', 'FMET', 'VCMET', 'VPMET', 'QMET', 'CLMET', ('METD', 'METC', 'METP')),
    'fle': ('KAFLE', 'FFLE', 'VCFLE', 'VPFLE', 'QFLE', 'CLFLE', ('FLED', 'FLEC', 'FLEP')),
    'bis': ('KABIS', 'FBIS', 'VCBIS', None, None, 'CLBIS', ('BISD', 'BISC', None)),
    'ver': ('KAVER', 'FVER', 'VCVER', None, None, 'CLVER', ('VERD', 'VERC', None)),
}


def pk_ss(drug, dose_mg, tau_h, offset_h, p, ndays=18):
    kaK, FK, VCK, VPK, QK, CLK, cmts = SPEC[drug]
    ka, F, VC = p[kaK] / 60.0, p[FK], p[VCK]
    CL = p[CLK] * (p['CYP2D6'] if drug == 'met' else 1.0) / 60.0
    two = VPK is not None
    VP_ = p[VPK] if two else 1.0
    Q_ = (p[QK] / 60.0) if two else 0.0

    def f(t, y):
        return np.array([
            -ka * y[0],
            F * ka * y[0] - CL / VC * y[1] - Q_ * (y[1] / VC - y[2] / VP_),
            Q_ * (y[1] / VC - y[2] / VP_)])

    tau = tau_h * 60.0
    n = max(2, int(ndays * 1440.0 / tau))
    y = np.zeros(3)
    for _ in range(n):
        y[0] += dose_mg
        y = solve_ivp(f, (0, tau), y, method='LSODA', rtol=1e-9, atol=1e-12).y[:, -1]
    y[0] += dose_mg
    y = solve_ivp(f, (0, offset_h * 60.0), y, method='LSODA',
                  rtol=1e-9, atol=1e-12).y[:, -1]
    out = {cmts[0]: y[0], cmts[1]: y[1]}
    if two:
        out[cmts[2]] = y[2]
    return out, y[1] / VC * 1e3
