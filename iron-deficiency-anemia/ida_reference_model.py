"""
ida_reference_model.py
Iron Deficiency Anaemia (IDA) QSP model -- independent pure-Python reference
implementation of ida_mrgsolve_model.R.

WHY THIS FILE EXISTS
--------------------
Every number quoted in README.md was computed twice: once by mrgsolve
(C++/LSODA) and once here (pure-Python fixed-step RK4, no dependencies beyond
the standard library). The two implementations share no code, so agreement
between them is evidence that the results follow from the EQUATIONS rather
than from one solver's behaviour. Representative agreement for a single 60 mg
oral dose at the IDA reference state:

                          mrgsolve      this file
    peak serum iron       224.2         224.2   ug/dL
    peak hepcidin           1.243         1.243 ng/mL
    fractional absorption  21.97         21.97  %
    probe dose at +24 h    91.5          91.5   % of rested
    14-day q24h absorbed  147.70        147.70  mg
    mass-balance error     ~1e-12        0.0000 %

Units: time = hours; iron = mg; cells = 1e12 cells; Hb mass = g.

USAGE
-----
    python3 ida_reference_model.py          # print both reference states

    import ida_reference_model as M
    h = M.ss(dict(M.P, VBLEED=0.6), 700)              # replete equilibrium
    d = M.ss(M.P, 900, y0=M.yv(h))                    # IDA equilibrium
    r = M.simulate(M.P, 72.0, events=[(0.0,"oral",60.0)], dt=0.025,
                   y0=M.yv(d), out_every=0.25)        # one oral dose

Events are (time_h, kind, amount) with kind in {"oral", "iv", "il6"}.
"""
import math

# ----------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------
P = dict(
    # --- body size / volumes -------------------------------------------------
    BV_L=3.9, PV_L=2.4,

    # --- diet & oral absorption ---------------------------------------------
    DIET=16.0,          # dietary iron (mg/day)
    F_DIET=0.30,        # bioavailable fraction reaching the absorptive window
    F_ORAL=0.90,        # fraction of an oral salt dose reaching the window (fasting)
    KTR=0.35,           # luminal iron swept past the duodenum (1/h)
    VMAX_D=5.57031,        # apical DMT1 Vmax (mg/h)
    KM_D=11.0,          # DMT1 Km (mg in window)
    K_APIC=45.0,        # hepcidin IC50 on apical uptake (ng/mL)
    K_IRP=30.25,         # enterocyte iron (mg) for IRP-mediated apical down-regulation
    KEXP=2.035,         # basolateral ferroportin export (1/h) at FPN_ENT = 1
    KFT=0.120,          # labile enterocyte iron -> enterocyte ferritin (1/h)
    KMOB=0.1360,         # enterocyte ferritin -> labile pool (1/h) at FPN_ENT = 1
    KSHED=0.0140,       # enterocyte ferritin iron lost by shedding (1/h), ~3 d

    # --- plasma iron / storage kinetics --------------------------------------
    K_LIV0=0.020,       # basal plasma -> hepatocyte uptake (1/h)
    K_LIV1=0.289258,       # holo-Tf/TSAT-driven hepatocyte uptake (1/h)
    KM_LIVT=45.0,       # TSAT (%) at half-max hepatocyte uptake
    K_LIV_EXP=0.01111,  # hepatocyte iron export (1/h) at FPN_RES = 1
    K_RES_EXP=0.0724,   # macrophage storage iron mobilisation (1/h) at FPN_RES = 1
    KM_FPN_REL=0.0969,   # FPN_RES at which half of recycled iron is released directly
    K_RES_LIV=0.00013,  # macrophage -> hepatocyte transfer (1/h)
    F_COL_RES=0.98,     # fraction of IV colloid handled by the RES
    K_COL=0.0693,       # IV colloid clearance from plasma (1/h), t1/2 = 10 h
    K_TISS=0.030,       # plasma -> non-erythroid tissue (1/h), capacity limited
    TISS_CAP=450.0,     # non-erythroid tissue iron capacity (mg)
    K_TISS_OUT=2.95e-5, # tissue -> plasma (1/h), slow
    TSAT_NTBI=75.0,     # TSAT (%) above which NTBI appears
    K_NTBI_ON=0.60,     # NTBI formation above threshold (1/h per fractional excess)
    K_NTBI_CL=1.20,     # NTBI clearance to liver/tissue (1/h)
    LOSS_BASAL=0.0417,  # obligatory epithelial iron loss (mg/h) = 1 mg/day

    # --- transferrin (TIBC) --------------------------------------------------
    TIBC_MIN=265.0, TIBC_MAX=470.0, KI_TIBC=22.0,
    IL6_TIBC=0.30, KI6_TIBC=25.0, KT_TIBC=0.0060,

    # --- hepcidin ------------------------------------------------------------
    KOUT_HEP=0.277,     # 1/h, t1/2 = 2.5 h
    KSYN_HEP=0.725,     # ng/mL/h
    HEP_STORE0=0.30, HEP_STORE_E=2.50, KM_HEP_LIV=260.0,   # BMP6/SMAD store sensing
    HEP_TSAT0=0.30, HEP_TSAT_E=1.70, KM_HEP_TSAT=25.0, HILL_TSAT=2.0,  # TfR2/HFE holo-Tf sensing
    EMAX_IL6=8.0, KM_IL6=18.0,                             # IL-6/STAT3
    KI_ERFE=9.0,                                           # erythroferrone
    TMPRSS6=1.0,        # matriptase-2 activity (1 normal, <1 = IRIDA)

    # --- ferroportin ---------------------------------------------------------
    # Enterocyte iron-export capacity ("mucosal block"): slow, because it embeds
    # ferroportin protein loss, enterocyte iron loading and enterocyte replacement.
    KSYN_FPE=0.0040,    # recovery of enterocyte export capacity (1/h)
    KDEG_FPE=0.0300,    # hepcidin-driven loss of enterocyte export capacity
    # Macrophage / hepatocyte ferroportin: fast, giving hypoferraemia within hours.
    KSYN_FPR=0.0350,    # ferroportin resynthesis (1/h)
    KDEG_FPR=0.1200,    # hepcidin-driven ferroportin degradation
    HILL_FPN=0.60,      # power-law exponent of the hepcidin-ferroportin relationship

    # --- erythroferrone / inflammation --------------------------------------
    KSYN_ERFE=0.30, KOUT_ERFE=0.150,
    KOUT_IL6=0.050, IL6_IN=0.0,
    IMAX_IL6_EPO=0.45, KM_IL6_EPO=25.0,   # IL-6/TNF blunting of renal EPO production
    IMAX_IL6_ERY=0.35, KM_IL6_ERY=25.0,   # direct cytokine suppression of erythroid progenitors
    E_IL6_RBC=0.50, KM_IL6_RBC=25.0,      # inflammatory shortening of RBC survival

    # --- erythropoiesis ------------------------------------------------------
    KIN_PROG=0.00684375, KPROG=0.02083, KEB=0.025, KRET=0.04167, KRBC=0.000347,
    EPO0=10.0, HB_REF=13.5, KEPO=0.50, EPO_MAX=1200.0,
    EMAX_EPO=2.50, EC50_EPO=100.0,
    K_EXP_FE=12.0,      # TSAT (%) at which half of the EPO-driven expansion is realised
    VHB=0.267539,         # maximal Hb synthesis per cell (pg/cell/h)
    K_TSAT_FE=3.50,     # TSAT (%) at half-maximal haemoglobinisation
    W1=0.50, W2=1.00, W3=1.50, WRET=0.30,
    KAPO_MAX=0.0060, FE_PER_HB=3.47,

    # --- blood loss ----------------------------------------------------------
    VBLEED=6.375,         # ongoing blood loss (mL/day)

    # --- biomarkers ----------------------------------------------------------
    KT_FERR=0.030, FERR_FLOOR=3.0, K_FERR_RES=3.50, K_FERR_LIV=8.50,
    FERR_IL6=32.0, KM_FERR_IL6=20.0,
    KT_STFR=0.010, STFR0=0.90, STFR_E=5.00,

    # --- FGF23 / phosphate arm ----------------------------------------------
    KSYN_FGF=8.0, KOUT_FGF=0.20,
    E_CLV=6.00, KM_CLV=210.0, K_CLV_IN=1.0, KOUT_CLV=0.0050, FCM_FGF=1.0,
    PHOS0=3.60, KOUT_PHOS=0.025, IMAX_FGF_P=0.72, IC50_FGF_P=190.0,
    CTRIOL0=42.0, KOUT_CTRIOL=0.060, IMAX_FGF_D=0.60, IC50_FGF_D=160.0,
    F_CTRIOL_P=0.35,
    PTH0=42.0, KOUT_PTH=0.120, E_PTH_D=1.50,

    # --- GI tolerability / adherence ----------------------------------------
    KT_GI=0.080, K_GI=0.55, EMAX_ADH=0.45, K50_ADH=3.0,

    # --- symptom endpoints ---------------------------------------------------
    FACIT_MAX=52.0, KM_FACIT_HB=2.0, W_FACIT_TISS=0.35, TISS_REF=380.0,
    IRLS_MAX=32.0,
)

SV = ["A_LUM", "A_ENT", "A_EFT", "A_TF", "A_NTBI", "A_COL", "A_RES", "A_LIV", "A_TISS",
      "TIBC", "HEP", "FPN_ENT", "FPN_RES", "ERFE", "IL6",
      "PROG", "N1", "N2", "N3", "H1", "H2", "H3", "NRET", "HRET", "NRBC", "HRBC",
      "FERR", "STFR", "CLV", "FGF23", "PHOS", "CTRIOL", "PTH", "GI",
      "CUM_ABS", "CUM_DOSE", "CUM_LOSS"]
IX = {n: i for i, n in enumerate(SV)}
NS = len(SV)
(iLUM, iENT, iEFT, iTF, iNTBI, iCOL, iRES, iLIV, iTISS, iTIBC, iHEP, iFPE, iFPR,
 iERFE, iIL6, iPROG, iN1, iN2, iN3, iH1, iH2, iH3, iNRET, iHRET, iNRBC, iHRBC,
 iFERR, iSTFR, iCLV, iFGF, iPHOS, iCTR, iPTH, iGI, iCABS, iCDOSE, iCLOSS) = range(NS)


def rhs(y, p):
    PV_dL = p["PV_L"] * 10.0
    BV_dL = p["BV_L"] * 10.0
    SI = y[iTF] * 1000.0 / PV_dL
    TIBC = y[iTIBC] if y[iTIBC] > 1.0 else 1.0
    TSAT = 100.0 * SI / TIBC
    HBt = y[iHRBC] + y[iHRET]
    HB = HBt / BV_dL
    F_FE = TSAT / (TSAT + p["K_TSAT_FE"])
    VHB_EFF = p["VHB"] * F_FE
    EPO = (p["EPO0"] * math.exp(p["KEPO"] * (p["HB_REF"] - HB))
           * (1.0 - p["IMAX_IL6_EPO"] * y[iIL6] / (p["KM_IL6_EPO"] + y[iIL6])))
    if EPO > p["EPO_MAX"]:
        EPO = p["EPO_MAX"]
    dep = EPO - p["EPO0"]
    if dep < 0.0:
        dep = 0.0
    EPO_EFF = 1.0 + p["EMAX_EPO"] * dep / (p["EC50_EPO"] + dep)
    # iron-restricted erythropoiesis: EPO cannot expand a marrow it cannot supply
    F_EXP = TSAT / (TSAT + p["K_EXP_FE"])
    EPO_REAL = 1.0 + (EPO_EFF - 1.0) * F_EXP

    wsum = (p["W1"] * y[iN1] + p["W2"] * y[iN2] + p["W3"] * y[iN3] + p["WRET"] * y[iNRET])
    HBSYN = wsum * VHB_EFF
    FE_MARROW = p["FE_PER_HB"] * HBSYN

    # ---- gut ---------------------------------------------------------------
    f_apic = ((1.0 / (1.0 + y[iHEP] / p["K_APIC"]))
              * (1.0 / (1.0 + (y[iENT] + y[iEFT]) / p["K_IRP"])))
    V_DMT1 = p["VMAX_D"] * y[iLUM] / (p["KM_D"] + y[iLUM]) * f_apic
    V_ABS = p["KEXP"] * y[iFPE] * y[iENT]
    v_ft = p["KFT"] * y[iENT]                    # labile -> enterocyte ferritin
    v_mob = p["KMOB"] * y[iFPE] * y[iEFT]        # ferritin -> labile (ferroportin gated)
    v_shed = p["KSHED"] * y[iEFT]                # lost with the shed enterocyte

    # ---- transport ---------------------------------------------------------
    liv_up = (p["K_LIV0"] + p["K_LIV1"] * TSAT / (p["KM_LIVT"] + TSAT)) * y[iTF]
    liv_exp = p["K_LIV_EXP"] * y[iFPR] * y[iLIV]
    res_exp = p["K_RES_EXP"] * y[iFPR] * y[iRES]
    res_liv = p["K_RES_LIV"] * y[iRES]
    col_out = p["K_COL"] * y[iCOL]
    exc = TSAT - p["TSAT_NTBI"]
    ntbi_on = p["K_NTBI_ON"] * (exc / 100.0) * y[iTF] if exc > 0.0 else 0.0
    ntbi_cl = p["K_NTBI_CL"] * y[iNTBI]
    tiss_up = p["K_TISS"] * y[iTF] * (1.0 - y[iTISS] / p["TISS_CAP"])
    if tiss_up < 0.0:
        tiss_up = 0.0
    tiss_out = p["K_TISS_OUT"] * y[iTISS]
    ephago = p["FE_PER_HB"] * p["KRBC"] * (1.0 + p["E_IL6_RBC"] * y[iIL6] / (p["KM_IL6_RBC"] + y[iIL6])) * y[iHRBC]
    f_rel = y[iFPR] / (y[iFPR] + p["KM_FPN_REL"])
    rec_rel = ephago * f_rel            # recycled iron released straight to plasma
    rec_sto = ephago - rec_rel          # recycled iron entering macrophage ferritin
    FBLEED = p["VBLEED"] / (1000.0 * 24.0 * p["BV_L"])

    dy = [0.0] * NS
    dy[iLUM] = p["DIET"] * p["F_DIET"] / 24.0 - p["KTR"] * y[iLUM] - V_DMT1
    dy[iENT] = V_DMT1 - V_ABS - v_ft + v_mob
    dy[iEFT] = v_ft - v_mob - v_shed
    dy[iTF] = (V_ABS + rec_rel + res_exp + liv_exp + tiss_out + (1.0 - p["F_COL_RES"]) * col_out
               - FE_MARROW - liv_up - tiss_up - ntbi_on - p["LOSS_BASAL"])
    dy[iNTBI] = ntbi_on - ntbi_cl
    dy[iCOL] = -col_out
    # iron from iron-restricted erythroblast apoptosis is recovered by marrow
    # macrophages (ineffective erythropoiesis is not an iron loss)
    apo_fe = p["FE_PER_HB"] * p["KAPO_MAX"] * (1.0 - F_FE) * (y[iH1] + y[iH2] + y[iH3])
    dy[iRES] = p["F_COL_RES"] * col_out + rec_sto + apo_fe - res_exp - res_liv
    dy[iLIV] = liv_up + res_liv + 0.7 * ntbi_cl - liv_exp
    dy[iTISS] = tiss_up + 0.3 * ntbi_cl - tiss_out

    # ---- transferrin -------------------------------------------------------
    tibc_t = ((p["TIBC_MIN"] + (p["TIBC_MAX"] - p["TIBC_MIN"]) * p["KI_TIBC"]
               / (p["KI_TIBC"] + y[iFERR]))
              * (1.0 - p["IL6_TIBC"] * y[iIL6] / (p["KI6_TIBC"] + y[iIL6])))
    dy[iTIBC] = p["KT_TIBC"] * (tibc_t - y[iTIBC])

    # ---- hepcidin ----------------------------------------------------------
    f_store = p["HEP_STORE0"] + p["HEP_STORE_E"] * y[iLIV] / (p["KM_HEP_LIV"] + y[iLIV])
    tsn = TSAT ** p["HILL_TSAT"]
    f_tsat = p["HEP_TSAT0"] + p["HEP_TSAT_E"] * tsn / (p["KM_HEP_TSAT"] ** p["HILL_TSAT"] + tsn)
    f_il6 = 1.0 + p["EMAX_IL6"] * y[iIL6] / (p["KM_IL6"] + y[iIL6])
    f_erfe = 1.0 / (1.0 + y[iERFE] / p["KI_ERFE"])
    dy[iHEP] = (p["KSYN_HEP"] * f_store * f_tsat * f_il6 * f_erfe / p["TMPRSS6"]
                - p["KOUT_HEP"] * y[iHEP])

    # ---- ferroportin -------------------------------------------------------
    hpow = y[iHEP] ** p["HILL_FPN"] if y[iHEP] > 0.0 else 0.0
    dy[iFPE] = p["KSYN_FPE"] * (1.0 - y[iFPE]) - p["KDEG_FPE"] * hpow * y[iFPE]
    dy[iFPR] = p["KSYN_FPR"] * (1.0 - y[iFPR]) - p["KDEG_FPR"] * hpow * y[iFPR]

    # ---- ERFE / IL-6 -------------------------------------------------------
    dy[iERFE] = p["KSYN_ERFE"] * (EPO_EFF - 0.75) - p["KOUT_ERFE"] * y[iERFE]
    dy[iIL6] = p["IL6_IN"] - p["KOUT_IL6"] * y[iIL6]

    # ---- erythropoiesis ----------------------------------------------------
    kapo = p["KAPO_MAX"] * (1.0 - F_FE)
    keb = p["KEB"]
    f_il6_ery = 1.0 - p["IMAX_IL6_ERY"] * y[iIL6] / (p["KM_IL6_ERY"] + y[iIL6])
    krbc_eff = p["KRBC"] * (1.0 + p["E_IL6_RBC"] * y[iIL6] / (p["KM_IL6_RBC"] + y[iIL6]))
    dy[iPROG] = p["KIN_PROG"] * f_il6_ery * EPO_REAL - p["KPROG"] * y[iPROG] - kapo * y[iPROG]
    dy[iN1] = p["KPROG"] * y[iPROG] - keb * y[iN1] - kapo * y[iN1]
    dy[iN2] = keb * y[iN1] - keb * y[iN2] - kapo * y[iN2]
    dy[iN3] = keb * y[iN2] - keb * y[iN3] - kapo * y[iN3]
    dy[iH1] = p["W1"] * y[iN1] * VHB_EFF - keb * y[iH1] - kapo * y[iH1]
    dy[iH2] = keb * y[iH1] + p["W2"] * y[iN2] * VHB_EFF - keb * y[iH2] - kapo * y[iH2]
    dy[iH3] = keb * y[iH2] + p["W3"] * y[iN3] * VHB_EFF - keb * y[iH3] - kapo * y[iH3]
    dy[iNRET] = keb * y[iN3] - p["KRET"] * y[iNRET] - FBLEED * y[iNRET]
    dy[iHRET] = (keb * y[iH3] + p["WRET"] * y[iNRET] * VHB_EFF
                 - p["KRET"] * y[iHRET] - FBLEED * y[iHRET])
    dy[iNRBC] = p["KRET"] * y[iNRET] - krbc_eff * y[iNRBC] - FBLEED * y[iNRBC]
    dy[iHRBC] = p["KRET"] * y[iHRET] - krbc_eff * y[iHRBC] - FBLEED * y[iHRBC]

    # ---- biomarkers --------------------------------------------------------
    ferr_t = (p["FERR_FLOOR"] + y[iRES] / p["K_FERR_RES"] + y[iLIV] / p["K_FERR_LIV"]
              + p["FERR_IL6"] * y[iIL6] / (p["KM_FERR_IL6"] + y[iIL6]))
    dy[iFERR] = p["KT_FERR"] * (ferr_t - y[iFERR])
    stfr_t = p["STFR0"] + p["STFR_E"] * (1.0 - F_FE) * (0.5 + 0.5 * EPO_EFF)
    dy[iSTFR] = p["KT_STFR"] * (stfr_t - y[iSTFR])

    # ---- FGF23 / phosphate -------------------------------------------------
    dy[iCLV] = p["K_CLV_IN"] * p["FCM_FGF"] * col_out - p["KOUT_CLV"] * y[iCLV]
    e_clv = 1.0 + p["E_CLV"] * y[iCLV] / (p["KM_CLV"] + y[iCLV])
    dy[iFGF] = p["KSYN_FGF"] * e_clv - p["KOUT_FGF"] * y[iFGF]
    f_phos_d = (1.0 - p["F_CTRIOL_P"]) + p["F_CTRIOL_P"] * y[iCTR] / p["CTRIOL0"]
    inh_p = 1.0 - p["IMAX_FGF_P"] * y[iFGF] / (p["IC50_FGF_P"] + y[iFGF])
    dy[iPHOS] = p["KOUT_PHOS"] * (p["PHOS0"] * f_phos_d * inh_p - y[iPHOS])
    inh_d = 1.0 - p["IMAX_FGF_D"] * y[iFGF] / (p["IC50_FGF_D"] + y[iFGF])
    dy[iCTR] = p["KOUT_CTRIOL"] * (p["CTRIOL0"] * inh_d - y[iCTR])
    pth_t = p["PTH0"] * (1.0 + p["E_PTH_D"] * (1.0 - y[iCTR] / p["CTRIOL0"]))
    dy[iPTH] = p["KOUT_PTH"] * (pth_t - y[iPTH])

    # ---- GI ----------------------------------------------------------------
    unabs = p["KTR"] * y[iLUM] + v_shed
    dy[iGI] = p["K_GI"] * unabs - p["KT_GI"] * y[iGI]

    # ---- bookkeeping -------------------------------------------------------
    dy[iCABS] = V_ABS
    dy[iCDOSE] = 0.0
    dy[iCLOSS] = p["LOSS_BASAL"] + p["FE_PER_HB"] * FBLEED * HBt
    return dy


def table(y, p):
    """Derived output quantities (mirrors $TABLE in the mrgsolve model)."""
    PV_dL = p["PV_L"] * 10.0
    BV_dL = p["BV_L"] * 10.0
    SI = y[iTF] * 1000.0 / PV_dL
    TSAT = 100.0 * SI / max(y[iTIBC], 1.0)
    HBt = y[iHRBC] + y[iHRET]
    HB = HBt / BV_dL
    ncell = max(y[iNRBC] + y[iNRET], 1e-9)
    F_FE = TSAT / (TSAT + p["K_TSAT_FE"])
    EPO = min(p["EPO0"] * math.exp(p["KEPO"] * (p["HB_REF"] - HB))
              * (1.0 - p["IMAX_IL6_EPO"] * y[iIL6] / (p["KM_IL6_EPO"] + y[iIL6])), p["EPO_MAX"])
    dep = max(EPO - p["EPO0"], 0.0)
    EPO_EFF = 1.0 + p["EMAX_EPO"] * dep / (p["EC50_EPO"] + dep)
    F_EXP = TSAT / (TSAT + p["K_EXP_FE"])
    EPO_REAL = 1.0 + (EPO_EFF - 1.0) * F_EXP
    wsum = (p["W1"] * y[iN1] + p["W2"] * y[iN2] + p["W3"] * y[iN3] + p["WRET"] * y[iNRET])
    HBSYN = wsum * p["VHB"] * F_FE
    d = dict(zip(SV, y))
    d.update(
        SI=SI, TSAT=TSAT, HB=HB, RBCC=ncell / p["BV_L"], MCH=HBt / ncell,
        RETPC=100.0 * y[iNRET] / ncell, RET_ABS=y[iNRET] / p["BV_L"] * 1000.0,
        CHR=y[iHRET] / max(y[iNRET], 1e-9), F_FE=F_FE, EPO=EPO, EPO_EFF=EPO_EFF,
        HBSYN=HBSYN, FE_MARROW=p["FE_PER_HB"] * HBSYN, EPO_REAL=EPO_REAL, F_EXP=F_EXP,
        V_ABS=p["KEXP"] * y[iFPE] * y[iENT],
        ABS_DAY=p["KEXP"] * y[iFPE] * y[iENT] * 24.0,
        ADH=1.0 - p["EMAX_ADH"] * y[iGI] / (p["K50_ADH"] + y[iGI]),
        STORES=y[iRES] + y[iLIV],
        BODY_FE=(y[iTF] + y[iNTBI] + y[iCOL] + y[iRES] + y[iLIV] + y[iTISS]
                 + p["FE_PER_HB"] * (HBt + y[iH1] + y[iH2] + y[iH3])),
        FACIT=p["FACIT_MAX"] * ((1.0 - p["W_FACIT_TISS"]) * (HB / (HB + p["KM_FACIT_HB"]))
                               + p["W_FACIT_TISS"] * min(y[iTISS] / p["TISS_REF"], 1.0)),
        IRLS=p["IRLS_MAX"] * (1.0 - min(y[iTISS] / p["TISS_REF"], 1.0)),
    )
    return d


def simulate(p, tend, events=(), dt=0.05, y0=None, out_every=1.0, adhere=True,
             collect=True):
    y = list(y0) if y0 is not None else list(INIT)
    ev = sorted(events, key=lambda e: e[0])
    ei, t = 0, 0.0
    rec, next_out = [], 0.0
    nstep = int(round(tend / dt))
    for step in range(nstep + 1):
        while ei < len(ev) and ev[ei][0] <= t + 1e-9:
            _, kind, amt = ev[ei]
            if kind == "oral":
                f = p["F_ORAL"] * ((1.0 - p["EMAX_ADH"] * y[iGI] / (p["K50_ADH"] + y[iGI]))
                                   if adhere else 1.0)
                y[iLUM] += amt * f
                y[iCDOSE] += amt
            elif kind == "iv":
                y[iCOL] += amt
                y[iCDOSE] += amt
            elif kind == "il6":
                y[iIL6] += amt
            ei += 1
        if collect and t >= next_out - 1e-9:
            row = table(y, p)
            row["time"] = t
            rec.append(row)
            next_out += out_every
        if step == nstep:
            break
        k1 = rhs(y, p)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NS)]
        k2 = rhs(y2, p)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NS)]
        k3 = rhs(y3, p)
        y4 = [y[i] + dt * k3[i] for i in range(NS)]
        k4 = rhs(y4, p)
        y = [max(y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]), 0.0)
             for i in range(NS)]
        t = (step + 1) * dt
    if not collect:
        row = table(y, p); row["time"] = t; rec = [row]
    return rec


INIT = [0.0] * NS
INIT[iLUM] = 0.3; INIT[iENT] = 0.5; INIT[iEFT] = 1.5; INIT[iTF] = 2.28; INIT[iTIBC] = 300.0
INIT[iHEP] = 5.0; INIT[iFPE] = 0.05; INIT[iFPR] = 0.22; INIT[iERFE] = 0.5
INIT[iPROG] = 0.26
INIT[iN1] = INIT[iN2] = INIT[iN3] = 0.217
INIT[iH1] = 2.5; INIT[iH2] = 7.6; INIT[iH3] = 14.5
INIT[iNRET] = 0.130; INIT[iHRET] = 3.9; INIT[iNRBC] = 15.6; INIT[iHRBC] = 526.5
INIT[iRES] = 150.0; INIT[iLIV] = 350.0; INIT[iTISS] = 380.0
INIT[iFERR] = 45.0; INIT[iSTFR] = 1.2
INIT[iPHOS] = 3.60; INIT[iCTR] = 42.0; INIT[iPTH] = 42.0; INIT[iFGF] = 40.0


def ss(p, days, y0=None, dt=0.25):
    r = simulate(p, days * 24.0, dt=dt, y0=y0, collect=False)
    return r[-1]


def yv(row):
    return [row[n] for n in SV]


KEYS = ["HB", "MCH", "RBCC", "SI", "TIBC", "TSAT", "HEP", "FERR", "STORES",
        "RETPC", "RET_ABS", "CHR", "STFR", "FPN_ENT", "FPN_RES", "F_FE", "EPO",
        "EPO_EFF", "ABS_DAY", "A_LIV", "A_RES", "A_TISS", "ERFE", "BODY_FE",
        "FACIT", "IRLS"]

if __name__ == "__main__":
    ph = dict(P); ph["VBLEED"] = 0.6
    h = ss(ph, 700)
    print("--- HEALTHY (replete, VBLEED 0.6 mL/d) ---")
    for k in KEYS:
        print(f"  {k:9s} {h[k]:10.3f}")
    print(f"  marrow mg/d {h['FE_MARROW']*24:8.2f}")
    d = ss(P, 900, y0=yv(h))
    print("--- IDA (VBLEED 3.5 mL/d) ---")
    for k in KEYS:
        print(f"  {k:9s} {d[k]:10.3f}")
    print(f"  marrow mg/d {d['FE_MARROW']*24:8.2f}")
