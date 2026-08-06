#!/usr/bin/env python3
# Pure-Python RK4 reference re-implementation of the CHI QSP model.
# Dependency-free (no numpy/scipy) so it can be run anywhere as an
# independent check on the mrgsolve model in chi_mrgsolve_model.R.
#
#   python3 chi_verify_run.py > chi_verification_output.txt
#
# Cross-validating these 36 ODEs against the R model exposed five real
# defects (basal amino-acid drive, insulin-only glycogen synthesis,
# nifedipine/sirolimus unit mismatch, 1000x glucagon dosing error, and
# diazoxide written as force-opening channels rather than shifting the
# channel ATP/ADP set-point). See README section 6.
"""
Reference re-implementation of the CHI QSP model (pure Python, RK4).
Used to (a) calibrate the handful of normal-physiology parameters and
(b) verify every number quoted in the README / R model.
"""
import math

# ----------------------------------------------------------------- parameters
P = dict(
    # --- beta-cell metabolic signal (ATP/ADP proxy) ---
    Rbas=1.0, Rmax=8.0, KG=8.0, hG=1.7,
    kAA=2.0, KAA=0.5, aGDH=3.0,        # GDH / leucine anaplerosis
    kMCT1=0.0,                          # pyruvate/lactate entry (SLC16A1)
    # --- K_ATP channel ---
    R50=2.0, nR=2.5, gKmax=12.0,
    EK=-75.0, Eleak=-20.0, gleak=1.0,
    # --- Ca / electrical ---
    V50=-50.0, kV=4.0, tauCa=0.02,
    # --- amplifying pathway ---
    ampmin=0.35, KA=7.0, hA=1.5, acamp=0.6,
    # --- secretion ---
    Smax_ins=1290.0, Sleak=0.0015,
    kelI=8.318,                         # insulin t1/2 5 min
    aCP=0.0333, kelCP=1.386,            # C-peptide t1/2 30 min
    kx=6.0,                             # insulin action tau 10 min
    # --- SSTR2 / octreotide ---
    gGIRK=1.2, EC50oct=1.0, kdesen=0.010, krecov=0.020,
    # --- glucose fluxes (mg/kg/min) ---
    Vg=0.22, UbrMax=6.44, Kbr=40.0, acbf=0.35, Gcbf=70.0,
    CMRreq=4.00, VketMax=3.0, Kket=6.0, VlacMax=0.8, Klac=4.0,
    Uii0=2.0, KUii=90.0, kren=0.02,
    Uidmax=20.0, KI=40.0, nI=2.0, G0id=90.0,
    kglyMax=3.60, KGLY=800.0, kgngMax=3.00,
    IC50gly=30.0, IC50gng=45.0, ng=2.0,
    ksyn=2.5, KsynI=10.0, Ksyng=100.0, GLYmax=4000.0,
    kabs=0.012,
    # --- lipolysis / ketogenesis ---
    kLipo=7.25, IC50lip=12.0, nl=1.8, kFFAox=10.0, aepiL=0.5,
    kket=0.4167, IC50ket=15.0, nk=2.0, kketox=1.5, agcgK=3.0, akg=2.5,
    # --- counterregulation ---
    kgcgS=563.0, kgcgE=5.199, agcg=3.0, Ghalf=70.0, kgs=8.0,
    IC50a=60.0, GCG0=100.0, KGCG=400.0, ashift=6.0,
    kepiS=695.0, kepiE=13.90, aepi=8.0, Gepi=60.0, kes=6.0, EPI0=50.0,
    kcortS=3.70, kcortE=0.462, acort=2.0, CORT0=8.0,
    # --- ammonia (GDH) ---
    kNH3=40.0, kNH3e=1.0, aNH3=0.917,
    kAAel=1.2, AA0=0.12,
    klacp=1.30, klacc=1.0, LAC0=1.3,
    # --- drug PK ---
    ka_dzx=0.60, Vc_dzx=0.25, CL_dzx=0.00866, Q_dzx=0.010, Vp_dzx=0.20,
    Emax_dzx=0.95, EC50_dzx=25.0, kdzx=0.60,
    ka_oct=1.80, V_oct=0.27, CL_oct=0.110,
    ka_gcg=6.0, V_gcg=0.25, fgcgconv=1000.0,
    Vc_ers=0.055, CL_ers=0.000160, Q_ers=0.00020, Vp_ers=0.045,
    Emax_ers=0.70, EC50_ers=45.0,
    ka_sir=1.0, V_sir=6.0, CL_sir=0.070, Emax_sir=0.35, EC50_sir=0.012,
    V_nif=1.2, CL_nif=0.35, Emax_nif=0.60, EC50_nif=1.20,
    CL_ex9=1.20, V_ex9=0.20, Emax_ex9=0.55, EC50_ex9=30.0, EX9inf=0.0,
    # --- disease axis ---
    g_ab=1.0, w_ab=0.0, dens_a=1.0, sGTP=1.0, KGshift=1.0, mct1=0.0,
    # --- beta-cell mass ---
    BMASS0=1.0, kprol=0.0010, kapo=0.00050, Bmax=1.3,
    # --- outcome ---
    kdev=1.0,
    # --- infusions / dosing (set per scenario) ---
    GIR_fix=0.0, feed=0.0, GCG_inf=0.0, Gtarget=70.0,
    Kp=0.10, Ki=0.60, loop=0.0, GIRmax=40.0,
)

S = ["GLU", "GLUi", "INS", "X", "CPEP", "GLY", "FFA", "BOHB", "LAC", "GCG",
     "EPI", "CORT", "NH3", "CABn", "CABa", "CAMP", "BMASS", "AA", "GGUT",
     "DZXg", "DZX", "DZXp", "OCTs", "OCT", "ROCT", "GCGs", "ERS", "ERSp",
     "SIRg", "SIR", "NIF", "EX9", "AUCHYPO", "TFUEL", "DEV", "GIRi"]
IX = {n: i for i, n in enumerate(S)}


def sig_low(x, x50, k):
    """1 at low x, 0 at high x."""
    return 1.0 / (1.0 + math.exp((x - x50) / k))


def derived(y, p):
    """All algebraic quantities. Returns dict."""
    GLU, GLUi, INS, X, CPEP, GLY, FFA, BOHB, LAC, GCG, EPI, CORT, NH3, \
        CABn, CABa, CAMP, BMASS, AA, GGUT, DZXg, DZX, DZXp, OCTs, OCT, ROCT, \
        GCGs, ERS, ERSp, SIRg, SIR, NIF, EX9, AUC, TF, DEV, GIRi = y
    d = {}
    Gb = max(GLU, 1e-6) / 18.016
    fGDH = p["aGDH"] * (1.0 - p["sGTP"])
    gdh = p["kAA"] * (AA / (p["KAA"] + AA)) * (1.0 + fGDH)
    KGe = p["KG"] * p["KGshift"]
    Rt = (p["Rbas"] + (p["Rmax"] - p["Rbas"]) * Gb ** p["hG"]
          / (KGe ** p["hG"] + Gb ** p["hG"]) + gdh
          + p["kMCT1"] * p["mct1"] * LAC)
    Fdzx = p["Emax_dzx"] * DZX / (p["EC50_dzx"] + DZX)
    R50e = p["R50"] * (1.0 + p["kdzx"] * Fdzx)   # diazoxide shifts the set-point
    PoATP = 1.0 / (1.0 + (Rt / p["R50"]) ** p["nR"])
    Po = 1.0 / (1.0 + (Rt / R50e) ** p["nR"])
    Foct = (OCT / (p["EC50oct"] + OCT)) * ROCT
    GGIRK = p["gGIRK"] * Foct
    GKn = p["gKmax"] * 1.0 * Po
    GKa = p["gKmax"] * p["g_ab"] * Po
    Vn = ((GKn + GGIRK) * p["EK"] + p["gleak"] * p["Eleak"]) / (GKn + GGIRK + p["gleak"])
    Va = ((GKa + GGIRK) * p["EK"] + p["gleak"] * p["Eleak"]) / (GKa + GGIRK + p["gleak"])
    Fnif = p["Emax_nif"] * NIF / (p["EC50_nif"] + NIF)
    fCan = (1.0 - Fnif) / (1.0 + math.exp(-(Vn - p["V50"]) / p["kV"]))
    fCaa = (1.0 - Fnif) / (1.0 + math.exp(-(Va - p["V50"]) / p["kV"]))
    amp = p["ampmin"] + (1 - p["ampmin"]) * Gb ** p["hA"] / (p["KA"] ** p["hA"] + Gb ** p["hA"])
    Fsir = p["Emax_sir"] * SIR / (p["EC50_sir"] + SIR)
    fcamp = 1.0 + p["acamp"] * (CAMP - 1.0)
    mass = (1.0 - p["w_ab"]) * CABn + p["w_ab"] * p["dens_a"] * CABa
    Ssec = p["Smax_ins"] * BMASS * (mass * amp * fcamp + p["Sleak"]) * (1.0 - Fsir)

    Fers = p["Emax_ers"] * ERS / (p["EC50_ers"] + ERS)
    Xeff = X

    # --- glucose fluxes ---
    fgcg = GCG / p["GCG0"]
    Egcg = GCG / (p["KGCG"] + GCG)              # glucagon receptor occupancy
    Egcg0 = p["GCG0"] / (p["KGCG"] + p["GCG0"])  # occupancy at basal glucagon
    # PKA (glucagon) and PP1 (insulin) compete on phosphorylase: raising cAMP
    # shifts the insulin IC50 for glycogenolysis to the right.
    IC50g_eff = p["IC50gly"] * (1.0 + p["ashift"] * max(0.0, Egcg - Egcg0)
                                / max(1e-9, 1.0 - Egcg0))
    finsG = 1.0 / (1.0 + (Xeff / IC50g_eff) ** p["ng"])
    finsN = 1.0 / (1.0 + (Xeff / p["IC50gng"]) ** p["ng"])
    Vgly = (p["kglyMax"] * (GLY / (p["KGLY"] + GLY)) * (Egcg / Egcg0)
            * (1.0 + 0.5 * (EPI / p["EPI0"] - 1.0)) * finsG)
    Vgng = (p["kgngMax"] * (0.4 + 0.6 * min(fgcg, 3.0)) * finsN
            * (1.0 + 0.25 * (CORT / p["CORT0"] - 1.0)) * (max(LAC, 0.1) / p["LAC0"]) ** 0.3)
    EGP = max(Vgly + Vgng, 0.0)
    Ra = p["kabs"] * GGUT
    Ubr = (p["UbrMax"] * GLU / (p["Kbr"] + GLU)
           * (1.0 + p["acbf"] * max(0.0, (p["Gcbf"] - GLU)) / p["Gcbf"]))
    Uket = p["VketMax"] * BOHB / (p["Kket"] + BOHB)
    Ulac = p["VlacMax"] * LAC / (p["Klac"] + LAC)
    FUEL = Ubr + Uket + Ulac
    fuelr = FUEL / p["CMRreq"]
    Uii = p["Uii0"] * GLU / (p["KUii"] + GLU)
    Uren = p["kren"] * max(0.0, GLU - 180.0)
    Uid = (p["Uidmax"] * Xeff ** p["nI"] / (p["KI"] ** p["nI"] + Xeff ** p["nI"])
           * (GLU / p["G0id"]))
    GIR = p["GIR_fix"]
    if p["loop"] > 0.5:
        GIR = max(0.0, min(p["GIRmax"],
                           p["Kp"] * (p["Gtarget"] - GLU) + GIRi))
    d.update(Gb=Gb, Rt=Rt, Po=Po, PoATP=PoATP, Vn=Vn, Va=Va, fCan=fCan, fCaa=fCaa,
             amp=amp, Ssec=Ssec, EGP=EGP, Vgly=Vgly, Vgng=Vgng, Ra=Ra, Ubr=Ubr,
             Uket=Uket, Ulac=Ulac, FUEL=FUEL, fuelr=fuelr, Uii=Uii, Uid=Uid,
             Uren=Uren, GIR=GIR, fgcg=fgcg, Fdzx=Fdzx, Foct=Foct, Fers=Fers,
             fGDH=fGDH, GGIRK=GGIRK, finsG=finsG, finsN=finsN, mass=mass,
             Egcg=Egcg, IC50g_eff=IC50g_eff, R50e=R50e)
    return d


def rhs(t, y, p):
    dd = derived(y, p)
    (GLU, GLUi, INS, X, CPEP, GLY, FFA, BOHB, LAC, GCG, EPI, CORT, NH3,
     CABn, CABa, CAMP, BMASS, AA, GGUT, DZXg, DZX, DZXp, OCTs, OCT, ROCT,
     GCGs, ERS, ERSp, SIRg, SIR, NIF, EX9, AUC, TF, DEV, GIRi) = y
    p_ = p
    kconv = 6.0 / p["Vg"]
    dy = [0.0] * len(y)
    Fers = dd["Fers"]

    dy[IX["GLU"]] = kconv * (dd["EGP"] + dd["Ra"] + dd["GIR"]
                             - dd["Uii"] - dd["Uid"] - dd["Uren"] - dd["Ubr"])
    dy[IX["GLUi"]] = (GLU - GLUi) / 0.17
    dy[IX["INS"]] = dd["Ssec"] - p["kelI"] * INS
    dy[IX["X"]] = p["kx"] * (INS * (1.0 - Fers) - X)
    dy[IX["CPEP"]] = p["aCP"] * dd["Ssec"] - p["kelCP"] * CPEP
    # hepatic glycogen synthesis is driven by PORTAL glucose delivery (Ra),
    # gated by insulin -- so it is ~0 in a true fast even though basal insulin
    # is non-zero.  This is what makes fasting tolerance a real prediction.
    Vsyn = (p["ksyn"] * dd["Ra"] * (X / (p["KsynI"] + X))
            * (GLU / (p["Ksyng"] + GLU)) * max(0.0, 1.0 - GLY / p["GLYmax"]))
    dy[IX["GLY"]] = 60.0 * (Vsyn - dd["Vgly"])
    fl = 1.0 / (1.0 + (X / p["IC50lip"]) ** p["nl"])
    dy[IX["FFA"]] = (p["kLipo"] * fl * (1.0 + p["aepiL"] * (EPI / p["EPI0"] - 1.0))
                     - p["kFFAox"] * FFA)
    fk = 1.0 / (1.0 + (X / p["IC50ket"]) ** p["nk"])
    switch = 1.0 + p["akg"] * max(0.0, 1.0 - GLY / 2000.0)
    dy[IX["BOHB"]] = (p["kket"] * FFA * fk * (1.0 + p["agcgK"] * (dd["fgcg"] - 1.0))
                      * switch - p["kketox"] * BOHB)
    ex = 0.0
    dy[IX["LAC"]] = p["klacp"] * (1.0 + ex) * (GLU / 90.0) - p["klacc"] * LAC
    falpha = 1.0 / (1.0 + X / p["IC50a"])
    dy[IX["GCG"]] = (p["kgcgS"] * (1.0 + p["agcg"] * sig_low(GLU, p["Ghalf"], p["kgs"]))
                     * falpha - p["kgcgE"] * GCG + p["fgcgconv"] * (p["ka_gcg"] * GCGs + p["GCG_inf"]) / p["V_gcg"])
    dy[IX["EPI"]] = (p["kepiS"] * (1.0 + p["aepi"] * sig_low(GLU, p["Gepi"], p["kes"]))
                     - p["kepiE"] * EPI)
    dy[IX["CORT"]] = (p["kcortS"] * (1.0 + p["acort"] * sig_low(GLU, 55.0, 6.0))
                      - p["kcortE"] * CORT)
    dy[IX["NH3"]] = p["kNH3"] * (1.0 + p["aNH3"] * dd["fGDH"]) - p["kNH3e"] * NH3
    dy[IX["CABn"]] = (dd["fCan"] - CABn) / p["tauCa"]
    dy[IX["CABa"]] = (dd["fCaa"] - CABa) / p["tauCa"]
    Fex9 = p["Emax_ex9"] * EX9 / (p["EC50_ex9"] + EX9)
    dy[IX["CAMP"]] = (1.0 - Fex9 - 0.55 * dd["Foct"] - CAMP) / 0.05
    dy[IX["BMASS"]] = (p["kprol"] * BMASS * max(0.0, 1.0 - BMASS / p["Bmax"])
                       - p["kapo"] * BMASS)
    dy[IX["AA"]] = -p["kAAel"] * (AA - p["AA0"])
    dy[IX["GGUT"]] = p["feed"] - 60.0 * p["kabs"] * GGUT
    # ---- drugs ----
    dy[IX["DZXg"]] = -p["ka_dzx"] * DZXg
    dy[IX["DZX"]] = (0.9 * p["ka_dzx"] * DZXg / p["Vc_dzx"]
                     - p["CL_dzx"] * DZX / p["Vc_dzx"]
                     - p["Q_dzx"] * (DZX - DZXp) / p["Vc_dzx"])
    dy[IX["DZXp"]] = p["Q_dzx"] * (DZX - DZXp) / p["Vp_dzx"]
    dy[IX["OCTs"]] = -p["ka_oct"] * OCTs
    dy[IX["OCT"]] = p["ka_oct"] * OCTs / p["V_oct"] - p["CL_oct"] * OCT / p["V_oct"]
    dy[IX["ROCT"]] = p["krecov"] * (1.0 - ROCT) - p["kdesen"] * ROCT * dd["Foct"]
    dy[IX["GCGs"]] = -p["ka_gcg"] * GCGs
    dy[IX["ERS"]] = (-p["CL_ers"] * ERS / p["Vc_ers"]
                     - p["Q_ers"] * (ERS - ERSp) / p["Vc_ers"])
    dy[IX["ERSp"]] = p["Q_ers"] * (ERS - ERSp) / p["Vp_ers"]
    dy[IX["SIRg"]] = -p["ka_sir"] * SIRg
    dy[IX["SIR"]] = p["ka_sir"] * SIRg / p["V_sir"] - p["CL_sir"] * SIR / p["V_sir"]
    dy[IX["NIF"]] = -p["CL_nif"] * NIF / p["V_nif"]
    dy[IX["EX9"]] = (p["EX9inf"] - p["CL_ex9"] * EX9) / p["V_ex9"]
    # ---- outcome integrals ----
    dy[IX["AUCHYPO"]] = max(0.0, 70.0 - GLU)
    dy[IX["TFUEL"]] = sig_low(dd["fuelr"], 1.0, 0.02)
    dy[IX["DEV"]] = p["kdev"] * max(0.0, 1.0 - dd["fuelr"])
    if p["loop"] > 0.5:
        dy[IX["GIRi"]] = p["Ki"] * (p["Gtarget"] - GLU)
    return dy


def y0(p):
    y = [0.0] * len(S)
    y[IX["GLU"]] = 75.0
    y[IX["GLUi"]] = 75.0
    y[IX["INS"]] = 5.0
    y[IX["X"]] = 5.0
    y[IX["CPEP"]] = 1.0
    y[IX["GLY"]] = 2000.0
    y[IX["FFA"]] = 0.60
    y[IX["BOHB"]] = 0.15
    y[IX["LAC"]] = 1.30
    y[IX["GCG"]] = 100.0
    y[IX["EPI"]] = 50.0
    y[IX["CORT"]] = 8.0
    y[IX["NH3"]] = 40.0
    y[IX["CABn"]] = 0.04
    y[IX["CABa"]] = 0.04
    y[IX["CAMP"]] = 1.0
    y[IX["BMASS"]] = p["BMASS0"]
    y[IX["AA"]] = p["AA0"]
    y[IX["ROCT"]] = 1.0
    return y


def run(p, tend, dt=0.001, y_init=None, events=None, record=None):
    y = y_init[:] if y_init else y0(p)
    t = 0.0
    n = int(round(tend / dt))
    out = []
    ev = sorted(events or [], key=lambda e: e[0])
    ei = 0
    for i in range(n):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            _, name, amt, mode = ev[ei]
            if mode == "add":
                y[IX[name]] += amt
            else:
                y[IX[name]] = amt
            ei += 1
        k1 = rhs(t, y, p)
        y2 = [y[j] + 0.5 * dt * k1[j] for j in range(len(y))]
        k2 = rhs(t + dt / 2, y2, p)
        y3 = [y[j] + 0.5 * dt * k2[j] for j in range(len(y))]
        k3 = rhs(t + dt / 2, y3, p)
        y4 = [y[j] + dt * k3[j] for j in range(len(y))]
        k4 = rhs(t + dt, y4, p)
        y = [y[j] + dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j])
             for j in range(len(y))]
        for j in (IX["GLU"], IX["INS"], IX["X"], IX["GLY"], IX["FFA"],
                  IX["BOHB"], IX["CPEP"], IX["NH3"]):
            if y[j] < 0:
                y[j] = 0.0
        t += dt
        if record and i % record == 0:
            out.append((t, y[:], derived(y, p)))
    return t, y, out


def P2(**kw):
    q = dict(P)
    q.update(kw)
    return q
