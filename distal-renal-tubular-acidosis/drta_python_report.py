#!/usr/bin/env python3
"""Shared 24-h summary reporter for the dRTA model."""
from drta import *


def summarize(out, m, last_h=24.0):
    """Mean / min / max over the final `last_h` hours."""
    t_end = out[-1][0]
    sel = [(t, y, L) for (t, y, L) in out if t >= t_end - last_h - 1e-9]
    n = len(sel)
    def mS(k):  return sum(L[k] for _, _, L in sel) / n
    def mY(k):  return sum(y[IX[k]] for _, y, _ in sel) / n
    def mnY(k): return min(y[IX[k]] for _, y, _ in sel)
    def mxY(k): return max(y[IX[k]] for _, y, _ in sel)
    def mnS(k): return min(L[k] for _, _, L in sel)
    def mxS(k): return max(L[k] for _, _, L in sel)
    y_end = out[-1][1]
    r = dict(
        HCO3=mY("HCO3_e"), HCO3min=mnY("HCO3_e"), HCO3max=mxY("HCO3_e"),
        pH=mS("pH_bl"), pCO2=mY("PaCO2"),
        uPH=mS("uph"), uPHmin=mnS("uph"), uPHmax=mxS("uph"),
        NAE=mS("NAE"), TA=mS("TA"), NH4=mS("NH4"), HCO3u=mS("HCO3u"),
        NEAP=mS("NEAP"), alk=mS("alk_in"), gap=mS("gap") * 24.0,
        VH=mS("VHc"), VHmax=mxS("VHc"),
        Jbone=mS("Jbone") * 24.0, CUMBASE=y_end[IX["CUMBASE"]],
        K=mY("K_pl"), Kmin=mnY("K_pl"), KDEF=y_end[IX["KDEF"]],
        Cl=mY("Cl_pl"),
        UCa=mS("UCa"), UCit=mS("UCit"), UVol=mS("UVol"),
        CaCit=mS("Ca_cit"), SS=mS("SS"), FE=mS("FE_HCO3"),
        BMDz=y_end[IX["BMDz"]], Hz=y_end[IX["Hz"]], BMIN=y_end[IX["BMIN"]],
        NC=y_end[IX["NC"]], STONE=y_end[IX["STONE"]], FIB=y_end[IX["FIB"]],
        NEPH=y_end[IX["NEPH"]], eGFR=m.p["GFR0"] * y_end[IX["NEPH"]],
        IGF1=y_end[IX["IGF1"]], MUS=y_end[IX["MUS"]], ADH=y_end[IX["ADH"]],
        bALP=y_end[IX["bALP"]], PTH=y_end[IX["PTH"]], OSM=y_end[IX["OSM"]],
        Pi=y_end[IX["Pi_pl"]], Ca=y_end[IX["Ca_pl"]],
        TBT=y_end[IX["TBT"]], AAC=y_end[IX["AAC"]],
        WASTE=y_end[IX["WASTE"]], GIVEN=y_end[IX["GIVEN"]],
        GM=mS("GM"), NDC1=mY("NDC1"),
    )
    r["waste_frac"] = r["WASTE"] / r["GIVEN"] if r["GIVEN"] > 0 else 0.0
    # fraction of the last 24 h with HCO3 below the responder threshold
    thr = 22.0
    r["frac_below"] = sum(1 for _, y, _ in sel if y[IX["HCO3_e"]] < thr) / n
    return r


def line(label, r):
    print(f"{label:30s} HCO3 {r['HCO3']:5.1f} ({r['HCO3min']:4.1f}-{r['HCO3max']:4.1f}) "
          f"pH {r['pH']:.3f} uPH {r['uPH']:.2f} NAE {r['NAE']:5.1f} "
          f"gap {r['gap']:6.1f} VH {r['VH']:.2f} K {r['K']:4.2f} Cl {r['Cl']:5.1f} "
          f"UCa {r['UCa']:5.2f} UCit {r['UCit']:5.2f} Ca/Cit {r['CaCit']:5.2f} "
          f"SS {r['SS']:5.2f} bone {r['Jbone']:6.1f}")


def run(par_upd, reg, t_end, dt=0.25, rec=1.0, seed=1):
    p = default_par()
    p.update(par_upd)
    out, m = simulate(p, reg, t_end=t_end, dt=dt, record_every=rec, seed=seed)
    return summarize(out, m), out, m


# canonical phenotypes -------------------------------------------------------
ADULT = dict(BW=70., AGE=35., BSA=1.80, GFR0=105.)
CHILD = dict(BW=30., AGE=10., BSA=1.05, GFR0=115.)
INFANT = dict(BW=12., AGE=2., BSA=0.54, GFR0=105.)

# dRTA lesion severities: LES = retained H+-pump Vmax, LES_grad = retained
# maximal blood->urine pH gradient
LES_COMPLETE   = dict(LES=0.18, LES_grad=0.58)   # complete dRTA
LES_SEVERE     = dict(LES=0.12, LES_grad=0.52)   # severe infantile (ATP6V0A4)
LES_INCOMPLETE = dict(LES=0.34, LES_grad=0.70)   # incomplete dRTA (normal HCO3)
