#!/usr/bin/env python3
"""Produce every number quoted in README.md, from the verified ODE system.

    python3 csdh_report.py > csdh_reference_output.txt
"""
import numpy as np
from csdh_verify import *

W = 78


def hdr(t):
    print("\n" + "=" * W); print(t); print("=" * W)


def sub(t):
    print("\n--- " + t + " " + "-" * max(0, W - 5 - len(t)))


# ---------------------------------------------------------------------
hdr("1.  NATURAL HISTORY OF AN INCIDENT cSDH  (a model OUTPUT, not an input)")
T, Y = run_in_trace(P, tmax=340.0)
print(f"{'day':>5} {'V(mL)':>8} {'d(mm)':>7} {'A(cm2)':>7} {'MLS':>6} {'ICP':>6} "
      f"{'J_ex':>6} {'J_rb':>6} {'J_abs':>6} {'haem':>6} {'NCAP':>6} {'HU':>6}")
for t in [0, 7, 14, 21, 28, 42, 60, 90, 150, 250, 340]:
    i = np.argmin(abs(T - t)); a = algebra(Y[i], P)
    print(f"{t:5.0f} {Y[i,IV_HEM]:8.2f} {a['dmax']:7.1f} {a['A_mem']:7.1f} "
          f"{a['MLS']:6.2f} {a['ICP']:6.1f} {a['J_ex']:6.2f} {a['J_rb']:6.2f} "
          f"{a['J_abs']:6.2f} {a['haem_drive']:6.3f} {Y[i,IN_CAP]:6.3f} {a['HU']:6.1f}")
_, lat, _ = run_in(P, 78.0, tmax=500.0)
print(f"\nlatency injury -> 78 mL presentation : {lat:.1f} days "
      f"(clinically 3-8 weeks)")
i = np.argmax(Y[:, IV_HEM])
print(f"plateau volume                       : {Y[-1,IV_HEM]:.1f} mL")
a = algebra(Y[-1], P)
print(f"at plateau: thickness {a['dmax']:.1f} mm, MLS {a['MLS']:.2f} mm, "
      f"ICP {a['ICP']:.1f} mmHg")
print(f"           cavity protein {a['C_PROT']:.2f} g/dL (reported 4-8), "
      f"HU {a['HU']:.0f} (isodense)")
print(f"           cavity tPA {a['C_TPA']:.1f} ng/mL, FDP {a['C_FDP']:.0f} ug/mL, "
      f"H {a['H']:.3f}")
print(f"           J_ex {a['J_ex']:.2f} + J_rb {a['J_rb']:.2f} "
      f"- J_abs {a['J_abs']:.2f} = {a['J_ex']+a['J_rb']-a['J_abs']:+.3f} mL/day")

# ---------------------------------------------------------------------
hdr("2.  THE 16 SCENARIOS")


def mkscn(tend=180):
    S = []
    S.append(Scn("01 conservative / no treatment", tend))
    S.append(Scn("02 burr-hole, no drain", tend).surgery(0, 0.25))
    S.append(Scn("03 burr-hole + drain 48h", tend).surgery(0, 0.25, drain_days=2.))
    S.append(Scn("04 burr-hole + drain, poor irrigation", tend)
             .surgery(0, 0.25, drain_days=2., wash=0.55))
    S.append(Scn("05 MMA embolisation alone", tend).embolise(0))
    S.append(Scn("06 burr-hole + drain + MMAE", tend)
             .surgery(0, 0.25, drain_days=2.).embolise(0))
    S.append(Scn("07 burr-hole + drain, MMAE at day 30", tend)
             .surgery(0, 0.25, drain_days=2.).embolise(30))
    S.append(Scn("08 dexamethasone alone (Dex-CSDH)", tend).dex())
    S.append(Scn("09 surgery + dexamethasone", tend)
             .surgery(0, 0.25, drain_days=2.).dex())
    S.append(Scn("10 atorvastatin 20 mg x 8 wk (ATOCH)", tend).atorva())
    S.append(Scn("11 atorvastatin + low-dose dex (ATOCH II)", tend)
             .atorva().dex(dose=2.25, days=35))
    S.append(Scn("12 tranexamic acid 750 mg/d", tend).txa())
    S.append(Scn("13 surgery + tranexamic acid", tend)
             .surgery(0, 0.25, drain_days=2.).txa())
    S.append(Scn("14 surgery + DOAC resumed day 7", tend)
             .surgery(0, 0.25, drain_days=2.).doac(7, tend))
    S.append(Scn("15 surgery + DOAC resumed day 30", tend)
             .surgery(0, 0.25, drain_days=2.).doac(30, tend))
    S.append(Scn("16 surgery + drain + MMAE + atorvastatin", tend)
             .surgery(0, 0.25, drain_days=2.).embolise(0).atorva())
    return S


print(f"{'scenario':<42} {'V90':>7} {'V180':>7} {'d180':>6} {'P(reop)':>8} "
      f"{'P(fav)':>7} {'P(thr)':>7}")
RES = {}
for s in mkscn():
    T, Y, p = s.run(dt=1.0)
    RES[s.name[:2]] = (T, Y, p)
    mrs = (0.55 * Y[:, IS_COG] + 0.45 * Y[:, IS_SYMP] + 0.9 * Y[:, IX_MYO]
           + 0.7 * Y[:, IX_INF] + 0.35 * (1 - Y[:, IN_NEUR]) / 0.3)
    pfav = 1.0 / (1.0 + np.exp((mrs - 1.55) / 0.42))
    print(f"{s.name:<42} {at(T,Y[:,IV_HEM],90):7.1f} {at(T,Y[:,IV_HEM],180):7.1f} "
          f"{at(T,series(T,Y,p,'dmax'),180):6.1f} "
          f"{1-np.exp(-at(T,Y[:,IH_REOP],180)):8.3f} "
          f"{at(T,pfav,180):7.3f} {1-np.exp(-at(T,Y[:,IH_THR],180)):7.4f}")

# ---------------------------------------------------------------------
hdr("3.  FLUX LEDGER: which term does each therapy touch?")
print("At day 30, relative to untreated, in mL/day:")
print(f"{'arm':<34} {'J_ex':>8} {'J_rb':>8} {'J_abs':>8} {'net':>8} {'Lp':>7} {'sigma':>7}")
for key, nm in [("01", "untreated"), ("03", "burr-hole + drain"),
                ("06", "+ MMA embolisation"), ("09", "+ dexamethasone"),
                ("10", "atorvastatin (conservative)"),
                ("13", "+ tranexamic acid")]:
    T, Y, p = RES[key]
    i = np.argmin(abs(T - 30))
    a = algebra(Y[i], p)
    print(f"{nm:<34} {a['J_ex']:8.2f} {a['J_rb']:8.2f} {a['J_abs']:8.2f} "
          f"{a['J_ex']+a['J_rb']-a['J_abs']:8.2f} {a['Lp_eff']:7.3f} "
          f"{a['sig_eff']:7.3f}")

# ---------------------------------------------------------------------
hdr("4.  BISTABILITY: the haem switch, probed by irrigation quality")
print("Same operation, same drain, only the irrigation differs.")
print(f"{'wash':>7} {'V30':>8} {'V90':>8} {'V180':>8} {'haem30':>8} {'NCAP30':>8} "
      f"{'P(reop)':>8}")
for w in [0.05, 0.10, 0.15, 0.25, 0.40, 0.55, 0.75, 1.00]:
    s = Scn(f"w{w}", 180).surgery(0, 0.25, drain_days=2., wash=w)
    T, Y, p = s.run(dt=1.0)
    hd = series(T, Y, p, 'haem_drive')
    print(f"{w:7.2f} {at(T,Y[:,IV_HEM],30):8.2f} {at(T,Y[:,IV_HEM],90):8.2f} "
          f"{at(T,Y[:,IV_HEM],180):8.2f} {at(T,hd,30):8.3f} "
          f"{at(T,Y[:,IN_CAP],30):8.3f} {1-np.exp(-at(T,Y[:,IH_REOP],180)):8.3f}")

sub("and with the Hill exponent set to 1 (switch removed)")
print(f"{'wash':>7} {'V180 (Hill=2)':>15} {'V180 (Hill=1)':>15}")
for w in [0.05, 0.25, 1.00]:
    a2 = Scn("h2", 180).surgery(0, 0.25, drain_days=2., wash=w).run(dt=1.5)
    a1 = Scn("h1", 180, p=dict(HB_HILL=1.0)).surgery(
        0, 0.25, drain_days=2., wash=w).run(dt=1.5)
    print(f"{w:7.2f} {at(a2[0],a2[1][:,IV_HEM],180):15.2f} "
          f"{at(a1[0],a1[1][:,IV_HEM],180):15.2f}")

# ---------------------------------------------------------------------
hdr("5.  V_RES APPEARS TWICE, WITH OPPOSITE SIGNS")
print("Atrophy buffers the pressure (protective) AND slows re-expansion (harmful).")
print(f"{'V_RES':>7} {'latency':>8} {'V@present':>10} {'MLS@pres':>9} {'k_reexp':>9} "
      f"{'Vcomp30':>9} {'V90':>8} {'P(reop)':>8}")
for vr in [5, 12, 20, 30, 40, 50, 60]:
    p = dict(P); p['V_RES'] = vr
    # presentation is defined by SYMPTOMS, so let the threshold be a shift, not
    # a volume: find the volume at which MLS first reaches 8 mm
    Tn, Yn = run_in_trace(p, tmax=400.)
    mls = series(Tn, Yn, p, 'MLS')
    j = np.argmax(mls >= 8.0) if (mls >= 8.0).any() else len(mls) - 1
    vpres = max(Yn[j, IV_HEM], 20.0)
    s = Scn(f"vr{vr}", 180, p=p, V_pres=vpres).surgery(0, 0.25, drain_days=2.)
    T, Y, pp = s.run(dt=1.0)
    kr = p['K_REEXP0'] * np.exp(-vr / p['LAM_ATR'])
    print(f"{vr:7.0f} {s.latency:8.1f} {vpres:10.1f} "
          f"{at(T,series(T,Y,pp,'MLS'),0):9.2f} {kr:9.4f} "
          f"{at(T,Y[:,IV_COMP],30):9.1f} {at(T,Y[:,IV_HEM],90):8.1f} "
          f"{1-np.exp(-at(T,Y[:,IH_REOP],180)):8.3f}")

# ---------------------------------------------------------------------
hdr("6.  MMA EMBOLISATION IS AN INTEGRAL, NOT AN EVENT")
Ta, Ya, pa = RES["03"]
Tb, Yb, pb = RES["06"]
print("Surgery+drain vs surgery+drain+MMAE. The absolute risk difference has to")
print("GROW with follow-up because the therapy removes the SOURCE, not the STOCK.")
print(f"{'day':>6} {'V ctrl':>9} {'V MMAE':>9} {'NCAP ctrl':>10} {'NCAP MMAE':>10} "
      f"{'P ctrl':>8} {'P MMAE':>8} {'ARD':>8}")
for t in [7, 14, 30, 60, 90, 120, 180]:
    pc = 1 - np.exp(-at(Ta, Ya[:, IH_REOP], t))
    pm = 1 - np.exp(-at(Tb, Yb[:, IH_REOP], t))
    print(f"{t:6.0f} {at(Ta,Ya[:,IV_HEM],t):9.2f} {at(Tb,Yb[:,IV_HEM],t):9.2f} "
          f"{at(Ta,Ya[:,IN_CAP],t):10.3f} {at(Tb,Yb[:,IN_CAP],t):10.3f} "
          f"{pc:8.3f} {pm:8.3f} {pc-pm:8.3f}")

sub("embolisation timing sweep (day of MMAE after the operation)")
print(f"{'MMAE day':>9} {'V90':>8} {'V180':>8} {'P(reop)180':>11}")
for td in [0, 7, 14, 30, 60, 90, -1]:
    if td < 0:
        s = Scn("none", 180).surgery(0, 0.25, drain_days=2.)
    else:
        s = Scn(f"e{td}", 180).surgery(0, 0.25, drain_days=2.).embolise(td)
    T, Y, p = s.run(dt=1.5)
    lbl = "never" if td < 0 else str(td)
    print(f"{lbl:>9} {at(T,Y[:,IV_HEM],90):8.2f} {at(T,Y[:,IV_HEM],180):8.2f} "
          f"{1-np.exp(-at(T,Y[:,IH_REOP],180)):11.3f}")

# ---------------------------------------------------------------------
hdr("7.  TRANEXAMIC ACID: 3 h IN PLASMA, ~6 DAYS AT THE SITE OF ACTION")
s = Scn("txa", 120).surgery(0, 0.25, drain_days=2.).txa(days=120)
T, Y, p = s.run(dt=0.25)
print(f"{'day':>7} {'C plasma (mg/L)':>16} {'C cavity (mg/L)':>16} {'ratio':>8} "
      f"{'V/J_abs (d)':>12}")
for t in [0.5, 1, 3, 7, 14, 21, 30, 60, 90, 120]:
    i = np.argmin(abs(T - t)); a = algebra(Y[i], p)
    cp = a['C_TXA']
    print(f"{t:7.1f} {cp:16.3f} {a['C_TXAH']:16.3f} "
          f"{(a['C_TXAH']/cp if cp>1e-9 else np.nan):8.2f} "
          f"{Y[i,IV_HEM]/max(a['J_abs'],1e-6):12.2f}")

# ---------------------------------------------------------------------
hdr("8.  DEXAMETHASONE: ONE DRUG, TWO SIGNS")


def pfav_of(T, Y):
    mrs = (0.55 * Y[:, IS_COG] + 0.45 * Y[:, IS_SYMP] + 0.9 * Y[:, IX_MYO]
           + 0.7 * Y[:, IX_INF] + 0.35 * (1 - Y[:, IN_NEUR]) / 0.3)
    return 1.0 / (1.0 + np.exp((mrs - 1.55) / 0.42))


Ta, Ya, pa = RES["03"]; Tb, Yb, pb = RES["09"]
print(f"{'endpoint':<34} {'placebo':>10} {'dex':>10}")
print(f"{'P(reoperation) at 180 d':<34} "
      f"{1-np.exp(-at(Ta,Ya[:,IH_REOP],180)):10.3f} "
      f"{1-np.exp(-at(Tb,Yb[:,IH_REOP],180)):10.3f}   <- helps locally")
print(f"{'P(favourable, mRS 0-3) at 180 d':<34} "
      f"{at(Ta,pfav_of(Ta,Ya),180):10.3f} {at(Tb,pfav_of(Tb,Yb),180):10.3f}"
      f"   <- harms globally")
print(f"{'peak glucose (mmol/L)':<34} {Ya[:,IGLU].max():10.2f} "
      f"{Yb[:,IGLU].max():10.2f}")
print(f"{'peak infection burden':<34} {Ya[:,IX_INF].max():10.3f} "
      f"{Yb[:,IX_INF].max():10.3f}")
print(f"{'peak myopathy burden':<34} {Ya[:,IX_MYO].max():10.3f} "
      f"{Yb[:,IX_MYO].max():10.3f}")
print(f"{'volume at 90 d (mL)':<34} {at(Ta,Ya[:,IV_HEM],90):10.2f} "
      f"{at(Tb,Yb[:,IV_HEM],90):10.2f}")

# ---------------------------------------------------------------------
hdr("9.  THE FIBRINOLYTIC LOOP L2")
Ta, Ya, pa = RES["01"]; Tb, Yb, pb = RES["12"]
print(f"{'quantity at day 60':<30} {'no TXA':>12} {'TXA':>12}")
for lbl, key in [("cavity tPA (ng/mL)", 'C_TPA'), ("plasmin (a.u.)", None),
                 ("cavity FDP (ug/mL)", 'C_FDP'),
                 ("haemostatic competence H", 'H'), ("J_rb (mL/day)", 'J_rb')]:
    if key is None:
        va, vb = at(Ta, Ya[:, IC_PLS], 60), at(Tb, Yb[:, IC_PLS], 60)
    else:
        va = at(Ta, series(Ta, Ya, pa, key), 60)
        vb = at(Tb, series(Tb, Yb, pb, key), 60)
    print(f"{lbl:<30} {va:12.3f} {vb:12.3f}")
Ha = at(Ta, series(Ta, Ya, pa, 'H'), 60); Hb = at(Tb, series(Tb, Yb, pb, 'H'), 60)
print(f"\nrebleed amplification 1/H : {1/Ha:.2f}x  ->  {1/Hb:.2f}x")

# ---------------------------------------------------------------------
hdr("10.  ANTICOAGULATION: THE TRADE IS EXPLICIT")
print(f"{'resumption day':>15} {'V180':>8} {'P(reop)':>9} {'P(thromb)':>11} "
      f"{'sum':>8}")
for td in [3, 7, 14, 30, 60, 180]:
    s = Scn(f"r{td}", 180).surgery(0, 0.25, drain_days=2.).doac(td, 180)
    T, Y, p = s.run(dt=1.5)
    pr = 1 - np.exp(-at(T, Y[:, IH_REOP], 180))
    pt = 1 - np.exp(-at(T, Y[:, IH_THR], 180))
    lbl = "never" if td >= 180 else str(td)
    print(f"{lbl:>15} {at(T,Y[:,IV_HEM],180):8.2f} {pr:9.3f} {pt:11.4f} "
          f"{pr+pt:8.3f}")

# ---------------------------------------------------------------------
hdr("11.  TRIAL LEDGER  (model vs published)")
Td, Yd, pd_ = RES["03"]; Tn, Yn2, pn = RES["02"]
Tm, Ym, pm_ = RES["06"]; Tx, Yx, px = RES["09"]
Tc, Yc, pc_ = RES["01"]; Tv, Yv, pv = RES["10"]
rows = [
 ("recurrence 90 d, burr-hole + drain",
  1 - np.exp(-at(Td, Yd[:, IH_REOP], 90)), "0.093", "Santarius 2009 Lancet"),
 ("recurrence 90 d, burr-hole no drain",
  1 - np.exp(-at(Tn, Yn2[:, IH_REOP], 90)), "0.240", "Santarius 2009 Lancet"),
 ("treatment failure 180 d, surgery alone",
  1 - np.exp(-at(Td, Yd[:, IH_REOP], 180)), "0.113", "EMBOLISE 2024 NEJM"),
 ("treatment failure 180 d, + MMAE",
  1 - np.exp(-at(Tm, Ym[:, IH_REOP], 180)), "0.041", "EMBOLISE 2024 NEJM"),
 ("repeat surgery, surgery + dexamethasone",
  1 - np.exp(-at(Tx, Yx[:, IH_REOP], 180)), "0.017", "Dex-CSDH 2020 NEJM"),
 ("favourable outcome, dexamethasone",
  at(Tx, pfav_of(Tx, Yx), 180), "0.839", "Dex-CSDH 2020 NEJM"),
 ("favourable outcome, placebo",
  at(Td, pfav_of(Td, Yd), 180), "0.903", "Dex-CSDH 2020 NEJM"),
 ("conversion to surgery 56 d, atorvastatin",
  1 - np.exp(-at(Tv, Yv[:, IH_REOP], 56)), "0.112", "ATOCH 2018 JAMA Neurol"),
 ("conversion to surgery 56 d, placebo",
  1 - np.exp(-at(Tc, Yc[:, IH_REOP], 56)), "0.235", "ATOCH 2018 JAMA Neurol"),
 ("latency injury -> presentation (days)", lat, "21-56", "clinical range"),
 ("cavity protein at plateau (g/dL)",
  algebra(Y_pl := run_in_trace(P, tmax=340.)[1][-1], P)['C_PROT'], "4-8",
  "cSDH fluid analyses"),
]
print(f"{'endpoint':<44} {'model':>8} {'published':>10}  source")
for nm, mv, pv_, src in rows:
    print(f"{nm:<44} {mv:8.3f} {pv_:>10}  {src}")

# ---------------------------------------------------------------------
hdr("12.  LOCAL SENSITIVITY  (+/- 20%)")
base = Scn("b", 180).surgery(0, 0.25, drain_days=2.)
Tb0, Yb0, pb0 = base.run(dt=1.5)
b90 = at(Tb0, Yb0[:, IV_HEM], 90)
br = 1 - np.exp(-at(Tb0, Yb0[:, IH_REOP], 180))
keys = ['LP_IMM', 'SIG_IMM', 'K_LYMPH', 'K_RB', 'K_FDP', 'K_TPA', 'HB50',
        'K_MAC', 'KD_CAP', 'K_MAT', 'MAT_SUPP', 'V_RES', 'K_REEXP0',
        'P_CAP', 'K_PV', 'D50_PAT', 'K_SPROUT', 'K_A', 'MMA_FLOOR', 'KD_MAC']
print(f"{'parameter':<12} {'dV90 (%)':>10} {'dP(reop) (%)':>13}")
sens = []
for k in keys:
    try:
        hi = Scn("h", 180, p={k: P[k] * 1.2}).surgery(0, 0.25, drain_days=2.).run(dt=1.5)
        lo = Scn("l", 180, p={k: P[k] * 0.8}).surgery(0, 0.25, drain_days=2.).run(dt=1.5)
        dv = (at(hi[0], hi[1][:, IV_HEM], 90) - at(lo[0], lo[1][:, IV_HEM], 90)) / max(b90, 1e-6) * 100
        dr = ((1 - np.exp(-at(hi[0], hi[1][:, IH_REOP], 180)))
              - (1 - np.exp(-at(lo[0], lo[1][:, IH_REOP], 180)))) / max(br, 1e-6) * 100
    except Exception:
        dv = dr = float('nan')
    sens.append((k, dv, dr))
    print(f"{k:<12} {dv:10.1f} {dr:13.1f}")
print("\nranked by |dP(reop)|:")
for k, dv, dr in sorted(sens, key=lambda r: -abs(r[2]))[:8]:
    print(f"   {k:<12} {dr:8.1f}%")

# ---------------------------------------------------------------------
hdr("13.  INTEGRITY CHECKS")
s = Scn("chk", 180).surgery(0, 0.25, drain_days=2.).atorva().txa()
T, Y, p = s.run(dt=0.5)
names = ['A_DEXD','A_DEXC','A_DEXP','CE_DEX','A_ATVD','A_ATVC','A_ATVP','CE_ATV',
         'A_TXAC','A_TXAP','A_TXAH','A_DOAC','A_DOAP','A_DOAH','V_HEM','M_HB',
         'M_PROT','M_FIB','M_FDP','M_TPA','C_PLS','C_PAI','N_CAP','N_MAT',
         'C_VEGF','C_ANG2','C_MMP9','C_IL6','C_IL8','N_MAC','M_MEMO','M_SEPT',
         'P_MMA','F_DRN','N_PC','F_PLT','C_FBG','V_COMP','X_MLSI','X_CBF',
         'N_NEUR','S_SYMP','S_COG','GLU','X_INF','X_MYO','H_REOP','H_THR',
         'H_DTH','A_EXT']
bad = [names[i] for i in range(NST) if Y[:, i].min() < -1e-6]
print("negative states:", ", ".join(bad) if bad else "none")
prot = series(T, Y, p, 'C_PROT')
print(f"cavity protein max {prot.max():.2f} g/dL  (must stay below plasma 7.0)")
icp = series(T, Y, p, 'ICP')
print(f"ICP range {icp.min():.2f} - {icp.max():.2f} mmHg")
fp = series(T, Y, p, 'f_pat')
print(f"patency range {fp.min():.3f} - {fp.max():.3f}")
H = series(T, Y, p, 'H')
print(f"haemostatic competence range {H.min():.3f} - {H.max():.3f}")
print(f"\ncalibrated parameter values actually used:")
for k in ['LP_IMM','K_LYMPH','K_MAC','KD_MAC','K_SPROUT','KD_CAP','LAM_REOP',
          'K_REEXP0','MMA_FLOOR','HB50','K_FUSE','HB_HILL','MAT_SUPP']:
    print(f"   {k:<12} = {P[k]:.6f}")
print("\n" + "=" * W)
print("All numbers above are produced by csdh_verify.py, the independent")
print("scipy re-implementation of the 50-ODE system in csdh_mrgsolve_model.R.")
print("=" * W)
