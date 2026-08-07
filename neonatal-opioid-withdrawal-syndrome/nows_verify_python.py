#!/usr/bin/env python3
"""
Independent Python/scipy re-implementation of the NOWS QSP model, used to
verify nows_mrgsolve_model.R.

WHY THIS FILE EXISTS
--------------------
The mrgsolve model cannot be executed in the environment this library is built
in (no R toolchain), so every ODE, every algebraic block and every parameter in
nows_mrgsolve_model.R is re-implemented here from the same equations and run
against a list of numeric anchors.  The file does two things:

  1. cross-checks all parameters against the $PARAM block of the R file by
     parsing it directly, so the two files cannot silently drift apart; and
  2. runs 48 anchors -- pharmacokinetic, structural, and clinical -- and
     prints PASS/FAIL for each.

Run:  python3 nows_verify_python.py            (scenarios + anchors)
      python3 nows_verify_python.py --anchors  (anchors only)

NOT FOR CLINICAL USE.
"""
import numpy as np
from scipy.integrate import solve_ivp

P = dict(
    WT0=3.30, GA=39.0,
    MDRUG=1, MDOSE=90.0, FBZD=0.0, FNIC=0.0, CARE=0.55, BF=0.0,
    ADFORCE=-1.0,
    PMA50U=54.2, HILLU=3.92, PMA50R=47.7, HILLR=3.40,
    PMA50C=52.0, HILLC=3.20, PNAIND=0.35, TIND=240.0,
    FORAL=0.30, KAM=1.10, CLMREF=0.95, V1MREF=3.40, V2MREF=5.44, QMREF=0.60,
    FM6G=0.30, CLGREF=0.10, VGREF=1.70, KE0M=0.45, KE0G=0.25,
    RPM=1.0, RPG=0.035,
    FDOR=0.75, KAD=0.35, CLDREF=0.737, VDREF=17.0, KE0D=0.30, RPD=3.0,
    CORD_D=0.55, KMPD=3.60, RIDD=0.025,
    FBSL=0.55, KAB=0.60, CLBREF=1.60, VBREF=40.0, KE0B=0.25,
    RPB=1.0, EMAXB=0.80, BADRV=0.40,
    FNORB=0.35, CLNREF=0.28, VNREF=15.0, KE0N=0.20, RPN=0.30,
    CORD_B=0.65, KMPB=0.28, CNB0=1.50, RIDB=0.005,
    FCLO=0.85, KACL=1.50, CLCREF=0.80, VCREF=10.2, KE0C=0.60,
    EC50A2=1.60, EMAXA2=0.34, KTACH=0.0100, KTOFF=0.0040,
    FPHB=0.90, KAP=0.80, CLPREF=0.0212, VPREF=3.06, EC50P=22.0, EMAXP=0.55,
    EC50MU=12.0, EC50BU=0.22, EMAXMU=1.0, KSHIFT=3.00,
    KRDOWN=0.0040, KRREC=0.0060,
    AGAIN=1.123, EA50A=0.160, GA50A=33.5, HGAA=8.0,
    ATMAX=0.180, EAT50=0.350, KATON=0.0167, KADOFF=0.002406, AMAX=1.0,
    LCB0=0.70, KENV=0.55, KLC=11.9, KCG=0.30, KINNE=1.10, KOUTNE=1.00,
    KINC=0.35, KOUTC=0.35, KINA=0.25, KOUTA=0.25, KING=0.12, KOUTG=0.12,
    WBZ=1.60, WNIC=0.45, KSEDO=1.30,
    SCNS=1.20, SANS=0.80, SGI=0.60, KFENV=1.55, KFS=0.35,
    THE=3.30, SE=0.55, THS=3.00, SS=0.55, THC=3.60, SC=0.55, KCARE=1.30,
    TAUE=6.0, TAUS=6.0, TAUC=4.0,
    KCALMAX=115.0, KCALBASE=58.0, KACT=6.0, KCALGAIN=5200.0,
    KELBZ=0.0069, KBZON=0.020, KBZOFF=0.0050, KELNIC=0.0058,
    FTARGET=6.0, THRSTART=8.0, KUP=0.0030, RUP=0.020, DMAX=1.00, DMIN=0.04,
    TAPMODE=1, TAPFRAC=0.10, TAPR=0.10/24.0, KTRK=0.05, GTOL=0.10,
    MINOBS=120.0, STABREQ=48.0, KLATCH=1.0, DINIT=0.32, KLM=0.50,
    ESCFAIL=0.95, TRTDRUG=1, CDOSE=0.0, PDOSE=0.0, PLOAD=0.0, ESCMODE=0,
)

S = ["AGUT_M","AC_M","AP_M","CE_M","A_G","CE_G","AGUT_D","AC_D","CE_D",
     "AGUT_B","AC_B","CE_B","AC_N","CE_N","AGUT_C","AC_C","CE_C","TACH",
     "AGUT_P","AC_P","AD","AT","RMU","NE","CNS","ANS","GI","SLP","EAT","CONS",
     "WT","BZD","BZW","NICW","CUMM","CUMD","SEIZH","AUCGAP","AUCSED",
     "FNAS_S","DOSE","TRTON","STAB","DSTAB"]
IDX = {n: i for i, n in enumerate(S)}
N = len(S)


def mat(pma, p50, hill):
    return pma**hill / (p50**hill + pma**hill)


def fetal_conc(p):
    if p["MDRUG"] == 1:
        return dict(cd=p["CORD_D"]*p["KMPD"]*p["MDOSE"], cb=0.0, cn=0.0)
    if p["MDRUG"] == 2:
        cb = p["CORD_B"]*p["KMPB"]*p["MDOSE"]
        return dict(cd=0.0, cb=cb, cn=p["CNB0"]*cb)
    if p["MDRUG"] == 3:
        return dict(cd=40.0, cb=0.0, cn=0.0)
    return dict(cd=0.0, cb=0.0, cn=0.0)


def effects(cem, ceg, ced, ceb, cen, A, rmu, p):
    UF = p["RPM"]*cem + p["RPG"]*ceg + p["RPD"]*ced + p["RPN"]*cen
    UB = p["RPB"]*ceb
    ec50 = p["EC50MU"]*(1.0 + p["KSHIFT"]*A)
    ec50b = p["EC50BU"]*(1.0 + p["KSHIFT"]*A)
    xf = UF/ec50
    xb = UB/ec50b
    den = 1.0 + xf + xb
    thf, thb = xf/den, xb/den
    EMU = rmu*(p["EMAXMU"]*thf + p["EMAXB"]*thb)
    EAD = rmu*(p["EMAXMU"]*thf + p["EMAXB"]*p["BADRV"]*thb)
    return EMU, EAD, thf, thb, ec50


def init_state(p):
    y = np.zeros(N)
    WT = p["WT0"]; y[IDX["WT"]] = WT
    c = fetal_conc(p)
    FGA = mat(p["GA"], p["GA50A"], p["HGAA"])
    A, r = 0.5, 0.8
    for _ in range(500):
        EMU, EAD, _, _, _ = effects(0, 0, c["cd"], c["cb"], c["cn"], A, r, p)
        r = p["KRREC"]/(p["KRREC"] + p["KRDOWN"]*EMU)
        A = p["AGAIN"]*FGA*EAD/(p["EA50A"] + EAD)
    EMU, EAD, _, _, _ = effects(0, 0, c["cd"], c["cb"], c["cn"], A, r, p)
    AT0 = p["ATMAX"]*EAD/(p["EAT50"] + EAD)
    AD0 = max(A - AT0, 0.0)
    if p["ADFORCE"] >= 0:
        AD0 = p["ADFORCE"]
    y[IDX["AD"]] = AD0; y[IDX["AT"]] = AT0; y[IDX["RMU"]] = r
    WS = WT/3.4
    if c["cd"] > 0:
        y[IDX["AC_D"]] = c["cd"]/1000.0*p["VDREF"]*WS*(0.25 if p["MDRUG"] == 3 else 1.0)
        y[IDX["CE_D"]] = c["cd"]
    if c["cb"] > 0:
        y[IDX["AC_B"]] = c["cb"]/1000.0*p["VBREF"]*WS; y[IDX["CE_B"]] = c["cb"]
        y[IDX["AC_N"]] = c["cn"]/1000.0*p["VNREF"]*WS; y[IDX["CE_N"]] = c["cn"]
    lc = p["LCB0"]*(1.0 + p["KENV"]*(1.0-p["CARE"]))
    ne = p["KINNE"]/p["KOUTNE"]*lc
    y[IDX["NE"]] = ne
    y[IDX["CNS"]] = ne; y[IDX["ANS"]] = ne; y[IDX["GI"]] = ne
    y[IDX["EAT"]] = 1.0/(1.0+np.exp((ne-p["THE"])/p["SE"]))
    y[IDX["SLP"]] = 1.0/(1.0+np.exp((ne-p["THS"])/p["SS"]))
    y[IDX["CONS"]] = 1.0/(1.0+np.exp((ne-p["THC"]-p["KCARE"]*p["CARE"])/p["SC"]))
    y[IDX["FNAS_S"]] = (p["SCNS"]+p["SANS"]+p["SGI"])*ne + p["KFENV"]*(1-p["CARE"])
    if p["FBZD"] > 0:
        y[IDX["BZD"]] = 1.0
    if p["FNIC"] > 0:
        y[IDX["NICW"]] = min(1.0, p["FNIC"]/20.0)
    y[IDX["AC_P"]] = p["PLOAD"]*WT
    return y


def rhs(t, y, p, rec=None):
    g = {n: y[i] for i, n in enumerate(S)}
    WT = max(g["WT"], 0.5)
    PMA = p["GA"] + t/168.0
    SZ = (WT/3.4)**0.75
    WS = WT/3.4
    RUGT = mat(PMA, p["PMA50U"], p["HILLU"])/mat(40.0, p["PMA50U"], p["HILLU"])
    RGFR = mat(PMA, p["PMA50R"], p["HILLR"])/mat(40.0, p["PMA50R"], p["HILLR"])
    RCYP = (mat(PMA, p["PMA50C"], p["HILLC"])/mat(40.0, p["PMA50C"], p["HILLC"])
            * (1.0 + p["PNAIND"]*(1.0-np.exp(-t/p["TIND"]))))

    CLM = p["CLMREF"]*SZ*RUGT; V1M = p["V1MREF"]*WS; V2M = p["V2MREF"]*WS
    QM = p["QMREF"]*SZ
    CLG = p["CLGREF"]*SZ*RGFR; VG = p["VGREF"]*WS
    CLD = p["CLDREF"]*SZ*RCYP; VD = p["VDREF"]*WS
    CLB = p["CLBREF"]*SZ*RCYP; VB = p["VBREF"]*WS
    CLN = p["CLNREF"]*SZ*RGFR; VN = p["VNREF"]*WS
    CLC = p["CLCREF"]*SZ*(0.4+0.6*RGFR); VC = p["VCREF"]*WS
    CLP = p["CLPREF"]*SZ*(0.5+0.5*RCYP); VP = p["VPREF"]*WS

    CM = 1000.0*g["AC_M"]/V1M; CPM = 1000.0*g["AP_M"]/V2M
    CG = 1000.0*g["A_G"]/VG; CD = 1000.0*g["AC_D"]/VD
    CB = 1000.0*g["AC_B"]/VB; CN = 1000.0*g["AC_N"]/VN
    CC = 1000.0*g["AC_C"]/VC; CP = g["AC_P"]/VP

    DOSE = max(g["DOSE"], 0.0)
    rate_m = rate_d = rate_b = 0.0
    if p["TRTDRUG"] == 1:
        rate_m = DOSE*WT/24.0
    elif p["TRTDRUG"] == 2:
        rate_d = DOSE*WT/24.0/6.0
    else:
        rate_b = DOSE*WT/24.0/20.0
    rate_c = p["CDOSE"]*WT/1000.0/24.0
    rate_p = p["PDOSE"]*WT/24.0
    bm_d = p["RIDD"]*(p["MDOSE"]/70.0)*WT/24.0*p["BF"] if p["MDRUG"] == 1 else 0.0
    bm_b = p["RIDB"]*(p["MDOSE"]/70.0)*WT/24.0*p["BF"] if p["MDRUG"] == 2 else 0.0

    dAGUT_M = -p["KAM"]*g["AGUT_M"] + rate_m
    dAC_M = (p["FORAL"]*p["KAM"]*g["AGUT_M"] - CLM*CM/1000.0
             - QM*CM/1000.0 + QM*CPM/1000.0)
    dAP_M = QM*CM/1000.0 - QM*CPM/1000.0
    dCE_M = p["KE0M"]*(CM - g["CE_M"])
    dA_G = p["FM6G"]*1.62*CLM*CM/1000.0 - CLG*CG/1000.0
    dCE_G = p["KE0G"]*(CG - g["CE_G"])
    dAGUT_D = -p["KAD"]*g["AGUT_D"] + rate_d + bm_d
    dAC_D = p["FDOR"]*p["KAD"]*g["AGUT_D"] - CLD*CD/1000.0
    dCE_D = p["KE0D"]*(CD - g["CE_D"])
    dAGUT_B = -p["KAB"]*g["AGUT_B"] + rate_b + bm_b
    dAC_B = p["FBSL"]*p["KAB"]*g["AGUT_B"] - CLB*CB/1000.0
    dCE_B = p["KE0B"]*(CB - g["CE_B"])
    dAC_N = p["FNORB"]*CLB*CB/1000.0 - CLN*CN/1000.0
    dCE_N = p["KE0N"]*(CN - g["CE_N"])
    dAGUT_C = -p["KACL"]*g["AGUT_C"] + rate_c
    dAC_C = p["FCLO"]*p["KACL"]*g["AGUT_C"] - CLC*CC/1000.0
    dCE_C = p["KE0C"]*(CC - g["CE_C"])
    dTACH = p["KTACH"]*(CC/(p["EC50A2"]+CC))*(1.0-g["TACH"]) - p["KTOFF"]*g["TACH"]
    dAGUT_P = -p["KAP"]*g["AGUT_P"] + rate_p
    dAC_P = p["FPHB"]*p["KAP"]*g["AGUT_P"] - CLP*CP

    A = min(g["AD"] + g["AT"], p["AMAX"])
    EMU, EAD, thf, thb, ec50 = effects(g["CE_M"], g["CE_G"], g["CE_D"],
                                       g["CE_B"], g["CE_N"], A, g["RMU"], p)
    dRMU = p["KRREC"]*(1.0-g["RMU"]) - p["KRDOWN"]*EMU*g["RMU"]
    EA2 = p["EMAXA2"]*g["CE_C"]/(p["EC50A2"]*(1.0+1.2*g["TACH"]) + g["CE_C"])
    ITONE = EMU + EA2 - EMU*EA2

    dAD = -p["KADOFF"]*g["AD"]
    dAT = p["KATON"]*(p["ATMAX"]*EAD/(p["EAT50"]+EAD) - g["AT"])

    GAP = A - ITONE
    GP = max(GAP, 0.0); GN = max(-GAP, 0.0)
    SED = min(p["EMAXP"]*CP/(p["EC50P"]+CP) + p["KSEDO"]*GN, 0.95)
    LC = (p["LCB0"]*(1.0 + p["KENV"]*(1.0-p["CARE"]))
          + p["KLC"]*GP*(1.0 - p["KCG"]*p["CARE"]))
    dNE = p["KINNE"]*LC - p["KOUTNE"]*g["NE"]

    dBZD = -p["KELBZ"]*g["BZD"]
    dBZW = (p["KBZON"]*p["FBZD"]*(1.0-g["BZD"])*(1.0-g["BZW"])
            - p["KBZOFF"]*g["BZW"])
    dNICW = -p["KELNIC"]*g["NICW"]
    dCNS = p["KINC"]*(g["NE"]+p["WBZ"]*g["BZW"])*(1.0-0.80*SED) - p["KOUTC"]*g["CNS"]
    dANS = p["KINA"]*(g["NE"]+p["WNIC"]*g["NICW"])*(1.0-0.50*SED) - p["KOUTA"]*g["ANS"]
    dGI = p["KING"]*g["NE"]*(1.0-0.30*SED) - p["KOUTG"]*g["GI"]
    FRAW = (p["SCNS"]*g["CNS"] + p["SANS"]*g["ANS"] + p["SGI"]*g["GI"]
            + p["KFENV"]*(1.0-p["CARE"]))
    dFNAS = p["KFS"]*(FRAW - g["FNAS_S"])

    Etar = (1.0/(1.0+np.exp((g["NE"]-p["THE"])/p["SE"])))*(1.0-0.60*SED)
    Star = (1.0/(1.0+np.exp((g["NE"]-p["THS"])/p["SS"])))*(1.0-0.20*SED)+0.20*SED
    Ctar = 1.0/(1.0+np.exp((g["NE"]-p["THC"]-p["KCARE"]*p["CARE"])/p["SC"]))
    dEAT = (Etar-g["EAT"])/p["TAUE"]
    dSLP = (Star-g["SLP"])/p["TAUS"]
    dCONS = (Ctar-g["CONS"])/p["TAUC"]
    dWT = WT*(p["KCALMAX"]*g["EAT"] - p["KCALBASE"] - p["KACT"]*g["NE"])/p["KCALGAIN"]/24.0

    if p["ESCMODE"] == 0:
        sig = g["FNAS_S"]; thr = p["THRSTART"]; tgt = p["FTARGET"]
    else:
        fail = 3.0 - (g["EAT"]+g["SLP"]+g["CONS"])
        sig = 4.0 + 12.0*max(fail-p["ESCFAIL"], 0.0); thr = 8.0; tgt = 6.0
    dTRTON = p["KLATCH"]*(1.0-g["TRTON"]) if sig > thr else 0.0
    djump = p["DINIT"]*dTRTON
    ERR = sig - tgt
    dSTAB = 0.0
    dDOSE = 0.0
    if g["TRTON"] > 0.5:
        if ERR > 0:
            dDOSE = min(p["KUP"]*ERR, p["RUP"]); dSTAB = -0.5*g["STAB"]
        else:
            dSTAB = 1.0
            if g["STAB"] > p["STABREQ"]:
                if p["TAPMODE"] == 0:
                    dDOSE = -p["TAPFRAC"]*DOSE/24.0
                elif p["TAPMODE"] == 1:
                    dDOSE = -p["TAPFRAC"]*g["DSTAB"]/24.0
                elif p["TAPMODE"] == 2:
                    dDOSE = -p["TAPFRAC"]*p["DINIT"]/24.0
                else:
                    atar = max(A - p["GTOL"], 0.0)
                    ctar = ec50*atar/max(p["EMAXMU"]*g["RMU"]-atar, 1e-3)
                    dDOSE = p["KTRK"]*(ctar/44.7 - DOSE)
    dDOSE += djump
    if DOSE <= p["DMIN"] and dDOSE < 0:
        dDOSE = 0.0
    if DOSE >= p["DMAX"] and dDOSE > 0:
        dDOSE = 0.0

    dDSTAB = p["KLM"]*(DOSE-g["DSTAB"]) if DOSE > g["DSTAB"] else 0.0
    dCUMM = (rate_m + 6.0*rate_d + 20.0*rate_b)/WT
    dCUMD = 1.0/24.0 if DOSE > p["DMIN"] else 0.0
    dSEIZH = 0.004*max(GP-0.55, 0.0)**2
    d = np.zeros(N)
    vals = dict(AGUT_M=dAGUT_M, AC_M=dAC_M, AP_M=dAP_M, CE_M=dCE_M, A_G=dA_G,
                CE_G=dCE_G, AGUT_D=dAGUT_D, AC_D=dAC_D, CE_D=dCE_D,
                AGUT_B=dAGUT_B, AC_B=dAC_B, CE_B=dCE_B, AC_N=dAC_N, CE_N=dCE_N,
                AGUT_C=dAGUT_C, AC_C=dAC_C, CE_C=dCE_C, TACH=dTACH,
                AGUT_P=dAGUT_P, AC_P=dAC_P, AD=dAD, AT=dAT, RMU=dRMU, NE=dNE,
                CNS=dCNS, ANS=dANS, GI=dGI, SLP=dSLP, EAT=dEAT, CONS=dCONS,
                WT=dWT, BZD=dBZD, BZW=dBZW, NICW=dNICW, CUMM=dCUMM, CUMD=dCUMD,
                SEIZH=dSEIZH, AUCGAP=GP, AUCSED=GN, FNAS_S=dFNAS, DOSE=dDOSE,
                TRTON=dTRTON, STAB=dSTAB, DSTAB=dDSTAB)
    for k, v in vals.items():
        d[IDX[k]] = v
    if rec is not None:
        rec.update(dict(EMU=EMU, EAD=EAD, EA2=EA2, ITONE=ITONE, GAP=GAP, A=A,
                        FRAW=FRAW, CM=CM, CD=CD, CB=CB, CN=CN, CG=CG, CC=CC,
                        CP=CP, SED=SED, LC=LC, CLM=CLM, RUGT=RUGT, EC50=ec50,
                        thf=thf, thb=thb, sig=sig))
    return d


DER = ["EMU","EAD","EA2","ITONE","GAP","A","FRAW","CM","CD","CB","CN","CG",
       "CC","CP","SED","LC","CLM","RUGT","EC50","thf","thb","sig"]


def run(p, tmax=45*24, derived=True):
    y0 = init_state(p)
    tt = np.arange(0, tmax+0.5, 1.0)
    sol = solve_ivp(rhs, (0, tmax), y0, args=(p,), t_eval=tt,
                    method="LSODA", rtol=1e-6, atol=1e-9, max_step=1.0)
    out = {n_: sol.y[i] for i, n_ in enumerate(S)}
    out["t"] = sol.t
    if not derived:
        for k in DER:
            out[k] = np.zeros(len(sol.t))
        out["GAP"] = out["A"] = np.zeros(len(sol.t))
        return out
    der = {k: [] for k in DER}
    for j in range(len(sol.t)):
        rec = {}
        rhs(sol.t[j], sol.y[:, j], p, rec)
        for k in DER:
            der[k].append(rec[k])
    for k in DER:
        out[k] = np.array(der[k])
    return out


def metrics(p, o):
    t = o["t"]
    trt = o["DOSE"] > p["DMIN"]
    treated = o["TRTON"].max() > 0.5
    dur = trt.sum()*1.0/24.0
    stop = t[np.max(np.where(trt)[0])]/24.0 if trt.any() else 0.0
    ok = (o["EAT"] > 0.55) & (o["SLP"] > 0.55) & (o["CONS"] > 0.55) & (~trt)
    ready = np.nan
    run_ = 0
    i0 = int(np.max(np.where(trt)[0])) if trt.any() else 0
    i0 = max(i0, int(p["MINOBS"]))
    for i in range(i0, len(t)):
        run_ = run_+1 if ok[i] else 0
        if run_ >= 48:
            ready = t[i]/24.0
            break
    return dict(A0=o["A"][0], peakF=o["FNAS_S"].max(),
                tpeak=t[o["FNAS_S"].argmax()]/24.0, treated=treated,
                days=dur, stop=stop, cum=o["CUMM"][-1], ready=ready,
                wt=o["WT"][-1], maxdose=o["DOSE"].max(),
                aucgap=o["AUCGAP"][-1], aucsed=o["AUCSED"][-1],
                gapstop=o["GAP"][int(min(len(t)-1, stop*24))],
                peakgap=o["GAP"].max())


def show(p, label="", tmax=45*24):
    o = run(p, tmax); m = metrics(p, o)
    print(f"{label:32s} A0={m['A0']:.3f} pkF={m['peakF']:5.1f}@d{m['tpeak']:4.1f} "
          f"trt={'Y' if m['treated'] else 'N'} Dmx={m['maxdose']:.2f} "
          f"days={m['days']:5.1f} stop=d{m['stop']:4.1f} rdy=d{m['ready']:5.1f} "
          f"cum={m['cum']:5.2f} wt={m['wt']:.2f}")
    return o, m


# ==========================================================================
#  PARAMETER CROSS-CHECK AGAINST nows_mrgsolve_model.R
# ==========================================================================
import os
import re


def parse_r_params(path=None):
    """Read the $PARAM @annotated block out of the mrgsolve file."""
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "nows_mrgsolve_model.R")
    if not os.path.exists(path):
        return None
    src = open(path).read()
    m = re.search(r"\$PARAM @annotated(.*?)\$CMT", src, re.S)
    if not m:
        return None
    out = {}
    for line in m.group(1).splitlines():
        mm = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*:\s*([-\d.eE+]+)\s*:", line)
        if mm:
            out[mm.group(1)] = float(mm.group(2))
    return out


def check_params():
    rp = parse_r_params()
    if rp is None:
        print("  [skip] nows_mrgsolve_model.R not found next to this file")
        return 0, 0
    ok = bad = 0
    only_py = sorted(set(P) - set(rp))
    only_r = sorted(set(rp) - set(P))
    for k in sorted(set(P) & set(rp)):
        a, b = float(P[k]), float(rp[k])
        if abs(a - b) <= 1e-6 + 1e-4 * max(abs(a), abs(b)):
            ok += 1
        else:
            bad += 1
            print(f"  MISMATCH {k}: python {a} vs R {b}")
    for k in only_py:
        if k in ("TAPR",):          # python-only convenience aliases
            continue
        bad += 1
        print(f"  ONLY IN PYTHON: {k}")
    for k in only_r:
        bad += 1
        print(f"  ONLY IN R: {k}")
    print(f"  parameters matched: {ok}   discrepancies: {bad}")
    return ok, bad


# ==========================================================================
#  ANCHORS
# ==========================================================================
RESULTS = []


def anchor(name, value, lo, hi, unit=""):
    ok = (value >= lo) and (value <= hi)
    RESULTS.append((name, value, lo, hi, unit, ok))
    return ok


def terminal_thalf(t, c):
    import numpy as _np
    c = _np.asarray(c)
    i0 = int(len(c) * 0.55)
    i1 = int(len(c) * 0.85)
    sel = slice(i0, i1)
    y = _np.log(_np.maximum(c[sel], 1e-12))
    x = _np.asarray(t)[sel]
    k = -_np.polyfit(x, y, 1)[0]
    return _np.log(2.0) / k if k > 0 else float("inf")


def run_anchors():
    import copy
    import numpy as _np
    b = copy.deepcopy(P)

    def v(**kw):
        q = copy.deepcopy(b); q.update(kw); return q

    # ---------- maturation arithmetic ----------
    fu40 = mat(40.0, b["PMA50U"], b["HILLU"])
    anchor("UGT2B7 capacity at PMA 40 wk (fraction of adult)", fu40, 0.225, 0.242)
    anchor("morphine CL ratio, PMA 44 vs 40 wk",
           mat(44., b["PMA50U"], b["HILLU"]) / fu40, 1.28, 1.34)
    anchor("morphine CL ratio, PMA 38 vs 34 wk",
           mat(38., b["PMA50U"], b["HILLU"]) / mat(34., b["PMA50U"], b["HILLU"]),
           1.40, 1.48)
    anchor("adaptation maturity FGA(34)/FGA(39)",
           mat(34., b["GA50A"], b["HGAA"]) / mat(39., b["GA50A"], b["HGAA"]),
           0.66, 0.72)

    # ---------- single-drug PK ----------
    anchor("morphine CL at term (L/h/kg)", b["CLMREF"] / 3.4, 0.25, 0.31, "L/h/kg")
    anchor("morphine Vss at term (L/kg)",
           (b["V1MREF"] + b["V2MREF"]) / 3.4, 2.3, 2.9, "L/kg")
    anchor("methadone CL at term (L/h/kg)", b["CLDREF"] / 3.4, 0.19, 0.24, "L/h/kg")
    anchor("methadone t1/2 at term (h)",
           0.693 * (b["VDREF"] / 3.4) / (b["CLDREF"] / 3.4), 14., 20., "h")
    anchor("M6G t1/2 at term (h)",
           0.693 * (b["VGREF"] / 3.4) / (b["CLGREF"] / 3.4), 10., 15., "h")
    anchor("clonidine t1/2 at term (h)",
           0.693 * (b["VCREF"] / 3.4) / (b["CLCREF"] / 3.4), 7., 13., "h")
    anchor("phenobarbital t1/2 at term (h)",
           0.693 * (b["VPREF"] / 3.4) / (b["CLPREF"] / 3.4), 85., 125., "h")
    anchor("cord methadone on 90 mg/day (ng/mL)",
           b["CORD_D"] * b["KMPD"] * 90., 140., 210., "ng/mL")
    anchor("cord buprenorphine on 16 mg/day (ng/mL)",
           b["CORD_B"] * b["KMPB"] * 16., 2.0, 4.0, "ng/mL")

    # ---------- washout half-life measured from the simulation ----------
    o = run(v(THRSTART=999), 8 * 24)
    anchor("simulated neonatal methadone t1/2 (h)",
           terminal_thalf(o["t"], o["CD"]), 13., 21., "h")

    # ---------- treatment concentrations ----------
    ob = run(b)
    on = ob["DOSE"] > b["DMIN"]
    cm = ob["CM"][on]
    anchor("morphine plasma during treatment, median (ng/mL)",
           float(_np.median(cm)), 12., 45., "ng/mL")
    anchor("morphine plasma during treatment, max (ng/mL)",
           float(cm.max()), 20., 60., "ng/mL")
    oc = run(v(CDOSE=6.0))
    anchor("clonidine plasma on 6 ug/kg/day (ng/mL)",
           float(_np.median(oc["CC"][int(10 * 24):])), 0.5, 1.3, "ng/mL")
    op = run(v(PDOSE=5.0, PLOAD=20.0))
    anchor("phenobarbital plasma on 5 mg/kg/day (mg/L)",
           float(_np.median(op["CP"][int(10 * 24):])), 12., 40., "mg/L")

    # ---------- structural invariants ----------
    anchor("ITONE never exceeds 1 (GIRK union)", float(ob["ITONE"].max()), 0.0, 1.0)
    anchor("mu + alpha-2 are sub-additive",
           float((oc["EMU"] + oc["EA2"] - oc["ITONE"]).min()), -1e-9, 1.0)
    anchor("receptor availability stays in (0,1]", float(ob["RMU"].min()), 0.05, 1.0)
    anchor("no negative drug amounts", float(min(ob["AC_M"].min(), ob["AC_D"].min(),
                                                 ob["A_G"].min())), -1e-6, 1e9)
    anchor("dose never exceeds DMAX", float(ob["DOSE"].max()), 0.0, b["DMAX"] + 1e-6)
    dose_int = float(_np.trapezoid(ob["DOSE"], ob["t"]) / 24.0)
    anchor("cumulative dose matches the dose integral (ratio)",
           float(ob["CUMM"][-1] / dose_int), 0.97, 1.03)
    anchor("weight stays positive", float(ob["WT"].min()), 1.0, 10.0, "kg")

    # ---------- adaptation block ----------
    a0 = {d: run(v(MDOSE=d, THRSTART=999), 24)["A"][0] for d in (40, 90, 160)}
    anchor("A0 on maternal methadone 90 mg/day", a0[90], 0.66, 0.73)
    anchor("A0 spread across 40-160 mg/day (relative)",
           (a0[160] - a0[40]) / a0[90], 0.0, 0.05)
    a0b = run(v(MDRUG=2, MDOSE=16, THRSTART=999), 24)["A"][0]
    anchor("A0 buprenorphine / A0 methadone", a0b / a0[90], 0.55, 0.80)

    # ---------- untreated natural history ----------
    ou = run(v(THRSTART=999))
    anchor("untreated peak Finnegan, methadone", float(ou["FNAS_S"].max()), 13.5, 18.)
    anchor("day of untreated peak", float(ou["t"][ou["FNAS_S"].argmax()] / 24), 3., 8., "d")
    anchor("baseline score in a well infant", float(ou["FNAS_S"][0]), 1.5, 4.5)
    t8 = ou["t"][ou["FNAS_S"] > 8][0] / 24
    anchor("onset day, methadone (score > 8)", float(t8), 2.0, 4.0, "d")
    os_ = run(v(MDRUG=3, THRSTART=999))
    anchor("onset day, short-acting opioid", float(os_["t"][os_["FNAS_S"] > 8][0] / 24),
           0.3, 1.6, "d")
    obu = run(v(MDRUG=2, MDOSE=16, THRSTART=999))
    anchor("onset day, buprenorphine", float(obu["t"][obu["FNAS_S"] > 8][0] / 24),
           2.0, 5.5, "d")
    anchor("untreated peak Finnegan, buprenorphine", float(obu["FNAS_S"].max()), 10., 14.)
    anchor("maternal-dose shift in onset day (40 -> 160 mg/day)",
           float(run(v(MDOSE=160, THRSTART=999))["t"][run(v(MDOSE=160, THRSTART=999))["FNAS_S"] > 8][0] / 24
                 - run(v(MDOSE=40, THRSTART=999))["t"][run(v(MDOSE=40, THRSTART=999))["FNAS_S"] > 8][0] / 24),
           0.7, 2.0, "d")

    # ---------- treated course ----------
    mb = metrics(b, ob)
    anchor("treatment days, methadone-exposed reference", mb["days"], 17., 30., "d")
    anchor("peak dose, methadone-exposed reference", mb["maxdose"], 0.55, 1.00, "mg/kg/d")
    anchor("weight gain per day on treatment (g/day)",
           (mb["wt"] - b["WT0"]) * 1000 / 45, 22., 42., "g/d")
    anchor("weight change per day, untreated (g/day)",
           (metrics(v(THRSTART=999), ou)["wt"] - b["WT0"]) * 1000 / 45, -15., 12., "g/d")

    mbu = metrics(v(MDRUG=2, MDOSE=16), run(v(MDRUG=2, MDOSE=16)))
    anchor("buprenorphine/methadone treatment-day ratio",
           mbu["days"] / mb["days"], 0.50, 0.85)
    anchor("buprenorphine/methadone cumulative-morphine ratio",
           mbu["cum"] / mb["cum"], 0.30, 0.65)

    # ---------- preterm ----------
    mpt = metrics(v(GA=34, WT0=2.10), run(v(GA=34, WT0=2.10)))
    mptT = metrics(v(GA=34, WT0=2.10, ADFORCE=0.578),
                   run(v(GA=34, WT0=2.10, ADFORCE=0.578)))
    anchor("preterm treatment days below term", mpt["days"] / mb["days"], 0.2, 0.85)
    anchor("preterm over-sedation exceeds term (ratio)",
           mpt["aucsed"] / max(mb["aucsed"], 1e-6), 5.0, 1e6)
    anchor("preterm with TERM adaptation is close to term (day ratio)",
           mptT["days"] / mb["days"], 0.80, 1.05)

    # ---------- weaning rules ----------
    W = {k: metrics(v(TAPMODE=k), run(v(TAPMODE=k))) for k in (0, 1, 2)}
    Wa = metrics(v(TAPMODE=3, GTOL=0.15), run(v(TAPMODE=3, GTOL=0.15)))
    anchor("weaning rules span > 8 treatment days",
           max(w["days"] for w in W.values()) - min(w["days"] for w in W.values()),
           8., 40., "d")
    anchor("slower wean buys lower withdrawal burden",
           W[0]["aucgap"] / W[1]["aucgap"], 0.20, 0.80)
    exch = (W[0]["days"] - W[1]["days"]) / (W[1]["aucgap"] - W[0]["aucgap"])
    anchor("exchange rate, treatment days per unit of burden avoided",
           exch, 0.20, 0.70, "d/unit")
    anchor("A-tracking oracle does not dominate the stabilisation rule",
           1.0 if (Wa["days"] >= W[1]["days"] or Wa["aucgap"] >= W[1]["aucgap"]) else 0.0,
           0.5, 1.5)

    # ---------- self-perpetuation ----------
    mnp = metrics(v(ATMAX=1e-6), run(v(ATMAX=1e-6)))
    anchor("days added by the treatment's own re-induced adaptation",
           mb["days"] - mnp["days"], 1.5, 8.0, "d")

    # ---------- adjuncts ----------
    mcl = metrics(v(CDOSE=6.0), run(v(CDOSE=6.0)))
    anchor("clonidine 6 ug/kg/day shortens treatment (fraction)",
           1 - mcl["days"] / mb["days"], 0.15, 0.45)
    mph = metrics(v(PDOSE=5, PLOAD=20), run(v(PDOSE=5, PLOAD=20)))
    anchor("phenobarbital shortens the course (pure opioid)",
           1 - mph["days"] / mb["days"], 0.10, 0.40)
    anchor("phenobarbital RAISES cumulative withdrawal burden",
           mph["aucgap"] / mb["aucgap"], 1.05, 2.0)
    anchor("phenobarbital costs weight gain (fraction of control)",
           (mph["wt"] - b["WT0"]) / (mb["wt"] - b["WT0"]), 0.10, 0.60)
    mbz = metrics(v(FBZD=1), run(v(FBZD=1)))
    anchor("benzodiazepine co-exposure lengthens treatment (fraction)",
           mbz["days"] / mb["days"] - 1, 0.15, 0.80)
    mbzp = metrics(v(FBZD=1, PDOSE=5, PLOAD=20), run(v(FBZD=1, PDOSE=5, PLOAD=20)))
    anchor("phenobarbital helps more in polysubstance than in pure opioid",
           (mbz["days"] - mbzp["days"]) - (mb["days"] - mph["days"]), 0.5, 12., "d")

    # ---------- care and ESC ----------
    mroom = metrics(v(CARE=0.90), run(v(CARE=0.90)))
    mesc = metrics(v(CARE=0.90, ESCMODE=1), run(v(CARE=0.90, ESCMODE=1)))
    anchor("rooming-in shortens treatment (fraction)", 1 - mroom["days"] / mb["days"],
           0.10, 0.40)
    anchor("ESC shortens treatment further (fraction vs Finnegan)",
           1 - mesc["days"] / mb["days"], 0.35, 0.70)
    anchor("ESC raises the withdrawal burden the infant actually bears",
           mesc["aucgap"] / mb["aucgap"], 1.2, 2.5)
    anchor("ESC burden stays below the untreated burden",
           mesc["aucgap"] / metrics(v(THRSTART=999), ou)["aucgap"], 0.2, 0.85)
    mbf = metrics(v(BF=1), run(v(BF=1)))
    anchor("breastfeeding shortens treatment (fraction)", 1 - mbf["days"] / mb["days"],
           0.15, 0.55)
    return RESULTS


def print_anchors(res):
    npass = sum(1 for r in res if r[5])
    print(f"\n{'ANCHOR':62s} {'VALUE':>10s}  {'RANGE':>18s}  RESULT")
    print("-" * 106)
    for name, val, lo, hi, unit, ok in res:
        print(f"{name:62s} {val:10.4g}  [{lo:8.3g},{hi:8.3g}]  {'PASS' if ok else 'FAIL'}")
    print("-" * 106)
    print(f"{npass}/{len(res)} anchors pass")
    return npass, len(res)


SCENARIOS = [
    ("methadone-exposed, untreated",        dict(THRSTART=999)),
    ("methadone-exposed + oral morphine",   dict()),
    ("maternal methadone 40 mg/day",        dict(MDOSE=40)),
    ("maternal methadone 160 mg/day",       dict(MDOSE=160)),
    ("buprenorphine-exposed",               dict(MDRUG=2, MDOSE=16)),
    ("buprenorphine-exposed, untreated",    dict(MDRUG=2, MDOSE=16, THRSTART=999)),
    ("short-acting opioid",                 dict(MDRUG=3)),
    ("preterm 34 wk, own adaptation",       dict(GA=34, WT0=2.10)),
    ("preterm 34 wk, TERM adaptation",      dict(GA=34, WT0=2.10, ADFORCE=0.578)),
    ("wean -10%/day of current dose",       dict(TAPMODE=0)),
    ("wean -10%/day of standard dose",      dict(TAPMODE=2)),
    ("wean by tracking A (oracle)",         dict(TAPMODE=3, GTOL=0.15)),
    ("titrate to a score of 3",             dict(FTARGET=3.0)),
    ("titrate to a score of 7.8",           dict(FTARGET=7.8)),
    ("+ clonidine 6 ug/kg/day",             dict(CDOSE=6)),
    ("+ clonidine 12 ug/kg/day",            dict(CDOSE=12)),
    ("+ phenobarbital, pure opioid",        dict(PDOSE=5, PLOAD=20)),
    ("benzodiazepine co-exposure",          dict(FBZD=1)),
    ("benzodiazepine + phenobarbital",      dict(FBZD=1, PDOSE=5, PLOAD=20)),
    ("rooming-in (CARE 0.90)",              dict(CARE=0.90)),
    ("Eat-Sleep-Console (CARE 0.90)",       dict(CARE=0.90, ESCMODE=1)),
    ("breastfeeding",                       dict(BF=1)),
    ("sublingual buprenorphine treatment",  dict(TRTDRUG=3)),
    ("oral methadone treatment",            dict(TRTDRUG=2)),
    ("no re-inducible adaptation pool",     dict(ATMAX=1e-6)),
]


def run_scenarios():
    import copy
    print(f"{'SCENARIO':38s} {'A0':>5s} {'peakF':>6s} {'day':>5s} {'trt':>4s} "
          f"{'Dmax':>5s} {'days':>6s} {'ready':>6s} {'mg/kg':>6s} {'AUCgap':>7s} "
          f"{'AUCsed':>7s} {'dWT g':>6s}")
    print("-" * 118)
    for nm, kw in SCENARIOS:
        p = copy.deepcopy(P); p.update(kw)
        o = run(p); m = metrics(p, o)
        print(f"{nm:38s} {m['A0']:5.3f} {m['peakF']:6.1f} {m['tpeak']:5.1f} "
              f"{'Y' if m['treated'] else 'N':>4s} {m['maxdose']:5.2f} "
              f"{m['days']:6.1f} {m['ready']:6.1f} {m['cum']:6.2f} "
              f"{m['aucgap']:7.1f} {m['aucsed']:7.1f} "
              f"{(m['wt']-p['WT0'])*1000:6.0f}")


if __name__ == "__main__":
    import sys
    print("=" * 118)
    print("NOWS QSP model — independent verification of nows_mrgsolve_model.R")
    print("=" * 118)
    print("\n[1] PARAMETER CROSS-CHECK  (python dict vs the $PARAM block of the R file)")
    ok, bad = check_params()
    if "--anchors" not in sys.argv:
        print("\n[2] SCENARIOS")
        run_scenarios()
    print("\n[3] ANCHORS")
    res = run_anchors()
    npass, ntot = print_anchors(res)
    sys.exit(0 if (bad == 0 and npass == ntot) else 1)
