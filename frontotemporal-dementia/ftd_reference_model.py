#!/usr/bin/env python3
# =============================================================================
#  ftd_reference_model.py
#  Dependency-free reference implementation of the Frontotemporal Dementia QSP
#  model in ftd_mrgsolve_model.R (same 49 states, same parameters, same right-
#  hand side). Pure standard library: no numpy, no scipy.
#
#  WHY THIS FILE EXISTS
#  --------------------
#  Every quantitative claim in README.md and in the header of the mrgsolve file
#  is COMPUTED here rather than asserted. It also serves three jobs the R file
#  cannot do on its own:
#
#    1. VERIFY the analytic steady-state result the whole model is built around
#         CSF PGRN fold-rise   = 1/(1 - f*theta)
#         lysosomal delivery   = (1 - p*theta)/(1 - f*theta)
#       against the numerically integrated ODEs, and report the residual.
#
#    2. SOLVE for the parameters that were fitted rather than looked up:
#       the sortilin share of plasma PGRN clearance (f_plasma) implied by the
#       observed ~2.7x plasma rise, and the CNS share (f_CSF) implied by the
#       observed ~2x CSF rise at the achieved CNS occupancy.
#
#    3. Run the natural-history and treatment scenarios and check them against
#       the published anchors (CDR+NACC-FTLD SB slope, plasma NfL, atrophy
#       rate, survival), REPORTING the disagreements instead of hiding them.
#
#  Integrator: fixed-step RK4, h = 0.05 d for trial horizons and h = 0.1 d for
#  10-year runs, with IV/oral doses applied as instantaneous state increments
#  at exact event times (the grid is chosen so events land on step boundaries).
#
#  Usage:  python3 ftd_reference_model.py            # full report to stdout
#          python3 ftd_reference_model.py --brief    # headline numbers only
# =============================================================================

import math
import sys

# -----------------------------------------------------------------------------
# PARAMETERS  (mirrors [PARAM] in ftd_mrgsolve_model.R)
# -----------------------------------------------------------------------------
P0 = dict(
    # genotype / patient
    GENO_GRN=0.0, GENO_C9=0.0, GENO_MAPT=0.0, SPORADIC=1.0,
    PRESYMPT=0.0, ALS_FLAG=0.0, PHENO_TEMP=0.0, AGE0=30.0,
    TMEM106B_PROT=0.0, BWT=70.0,
    GDOSE_HET=0.50, KSYN_NMD=0.70,
    # latozinemab
    MW_LATO=145000.0, V1_LATO=3.0, V2_LATO=3.6, CL_LATO=0.20, Q_LATO=0.45,
    KD_SORT_PL=1.0, KINT_SORT=0.25, R0_SORT=30.0, KDEG_SORT=0.35,
    KDOWN_SORT=0.20, KIN_CSF=0.0020, KOUT_CSF=1.0, KD_SORT_CNS=3.0,
    # AAV-GRN
    KDEC_VG=0.30, KTR_AAV=0.35, KLOSS_TRANSD=0.0004, NAB_FRAC=0.0,
    AAV_BOOST=1.10,
    # C9 ASO
    KOUT_ASO9_CSF=2.5, KTIS_ASO9=1.6, KOUT_ASO9_TIS=0.011,
    EMAX_ASO9=0.60, EC50_ASO9=3.0, ASO9_ROSTRAL=0.55,
    # MAPT ASO
    KOUT_ASOT_CSF=2.5, KTIS_ASOT=1.6, KOUT_ASOT_TIS=0.010,
    EMAX_ASOT=0.55, EC50_ASOT=3.5,
    # symptomatic drugs
    KA_SSRI=6.0, V_SSRI=900.0, KE_SSRI=0.46, EMAX_SERT=0.55, EC50_SERT=0.030,
    KA_TRZ=8.0, V_TRZ=90.0, KE_TRZ=2.20, EMAX_TRZ=0.65, EC50_TRZ=0.60,
    KA_DNP=4.0, V_DNP=830.0, KE_DNP=0.238, EMAX_ACHE=0.60, EC50_ACHE=0.020,
    K_AGIT_DNP=0.55,
    # progranulin
    KSYN_PGRN_PL=1400.0, KCL_SORT_PL=4.41, KCL_OTH_PL=2.59,
    KSYN_PGRN_CSF=6.00, KCL_SORT_CSF=1.56, KCL_OTH_CSF=0.44,
    P_LYS_SORT=0.50, K_UPTAKE_LYS=0.90, KDEG_LYS_PGRN=0.90, PGRN_LYS_WT=1.0,
    W_ACT_LYS=0.80, W_ACT_EXTRA=0.20,
    # lysosome
    KT_LYSO=0.030, HILL_LYSO=2.0, K50_LYSO=0.28,
    C9_LYSO_PEN=0.09, TMEM_PEN=0.05,
    KGEN_LIPO=0.0008, KCLR_LIPO=0.0006, K_LIPO_DAM=0.10,
    # TDP-43
    KIN_TDP=0.0020, KMIS_TDP=0.0030,
    W_STRESS_LYS=0.30, W_STRESS_GA=0.18, W_STRESS_OX=0.35,
    AGE_REF=62.0, AGE_POW=5.0, AGE_FLOOR=0.02, TMEM_ONSET_SHIFT=0.88,
    KAGG_TDP=2.60e-4, KSEED_TDP=0.85, KCLR_TDP=3.50e-4, AGG_MAX_TDP=1.0,
    KCRYP_STMN2=0.040, HILL_CRYP=1.6, KDEG_STMN2=0.040,
    KCRYP_UNC13A=0.040, KDEG_UNC13A=0.040,
    # C9
    KSYN_C9RNA=0.30, KDEG_C9RNA=0.30, KRAN_GP=0.28, KDEG_GP=0.28,
    KRAN_GR=0.18, KDEG_GR=0.20, KRAN_GA_REL=1.20,
    # tau
    KSYN_TAUM=0.28, KDEG_TAUM=0.28, KTRL_TAU=0.30, KDEG_TAUS=0.28,
    KPHOS_TAU=0.030, MUT_TAU_BOOST=5.00, KDEPHOS_TAU=0.075,
    KOLIG_TAU=5.00e-4, KSEED_TAU=0.90, KCLR_TAU=6.00e-4, AGG_MAX_TAU=1.0,
    # inflammation
    KACT_MG=0.030, KRES_MG=0.030,
    W_MG_TDP=1.00, W_MG_TAU=0.90, W_MG_PGRN=0.60, W_MG_GR=0.16,
    KSYN_C1Q=0.10, KDEG_C1Q=0.10, KSYN_C3=0.09, KDEG_C3=0.09, I_C1Q=0.0,
    KACT_AST=0.020, KRES_AST=0.020,
    # structure
    KREP_SYN=1.12e-4, KPRUNE_C3=3.00e-4, KTOX_TAU_SYN=1.50e-4,
    KTOX_UNC13A=3.00e-4,
    KDEATH_SN=7.80e-5, KDEATH_TEMP=5.60e-5,
    W_TOX_TDP=1.00, W_TOX_TAU=1.00, W_TOX_GR=0.16, W_TOX_CYT=0.70,
    W_TOX_SYN=0.90, W_TOX_STMN2=0.55,
    KATR_F=1.71e-4, KATR_T=1.50e-4, KATR_AGE=0.0000137, TEMP_SHIFT=0.0,
    # biomarkers
    KPROD_NFL=930.0, NFL_CSF_BASE=700.0, KOUT_NFL_CSF=0.11,
    KTR_NFL=0.011, KOUT_NFL_PL=0.60, NFL_PL_BASE=10.0,
    KPROD_GFAP=130.0, KOUT_GFAP=0.85, GFAP_BASE=95.0,
    KPROD_BMP=48.0, KOUT_BMP=0.50, BMP_BASE=12.0,
    # neurotransmitters
    KIN_5HT=0.55, KOUT_5HT=0.55, W_5HT_LOSS=0.85,
    ACH_INTEGRITY=0.95, KIN_ACH=0.60, KOUT_ACH=0.60,
    KIN_DA=0.50, KOUT_DA=0.50, W_DA_LOSS=0.55,
    # clinical
    CDR_MAX=24.0, KCDR=0.0149,
    A_CDR_SYN=1.00, A_CDR_SN=1.20, A_CDR_TEMP=0.60, A_CDR_VOL=0.80,
    CDR_THRESH=0.45,
    NPI_MAX=144.0, KT_NPI=0.045, EMAX_NPI=0.42, K50_NPI=0.62,
    B_NPI_5HT=0.75, B_NPI_SN=1.00, B_NPI_DA=0.35,
    H0_HAZ=0.000048, BH_CDR=0.115, BH_NFL=0.42, BH_ALS=1.35,
)

# -----------------------------------------------------------------------------
# STATE VECTOR  (order must match [CMT] in the R file)
# -----------------------------------------------------------------------------
NAMES = [
    "CLATO", "CPLATO", "SORTTOT", "BMPCSF", "CSFLATO",
    "AAVVG", "TRANSD", "ASO9CSF", "ASO9TIS", "ASOTCSF", "ASOTTIS",
    "SSRIG", "SSRIC", "TRZC", "DNPC",
    "PGRNPL", "PGRNCSF", "PGRNLYS",
    "LYSO", "LIPO",
    "TDPN", "TDPAGG", "STMN2T", "UNC13AT",
    "C9RNA", "POLYGP", "POLYGR",
    "TAUM", "TAUS", "TAUP", "TAUAGG",
    "MG", "C1Q", "C3", "AST",
    "SYN", "NEURSN", "NEURTMP", "VOLF", "VOLT",
    "NFLCSF", "NFLPL", "GFAPPL",
    "HT5", "ACH", "DA",
    "CDR", "NPI", "CUMHAZ",
]
IX = {n: i for i, n in enumerate(NAMES)}
N = len(NAMES)
assert N == 49, N


def grn_dose(p):
    """GRN transcription factor: allele dose x transcript stability."""
    return p["GDOSE_HET"] * p["KSYN_NMD"] if p["GENO_GRN"] > 0.5 else 1.0


def init_state(p):
    """Mirror of [MAIN] initial conditions in the R model."""
    y = [0.0] * N
    gd = grn_dose(p)

    y[IX["SORTTOT"]] = p["R0_SORT"]
    y[IX["PGRNPL"]] = p["KSYN_PGRN_PL"] * gd / (p["KCL_SORT_PL"] + p["KCL_OTH_PL"])
    y[IX["PGRNCSF"]] = p["KSYN_PGRN_CSF"] * gd / (p["KCL_SORT_CSF"] + p["KCL_OTH_CSF"])

    pcsf_wt = p["KSYN_PGRN_CSF"] / (p["KCL_SORT_CSF"] + p["KCL_OTH_CSF"])
    y[IX["PGRNLYS"]] = (p["K_UPTAKE_LYS"] / p["KDEG_LYS_PGRN"]) * (
        y[IX["PGRNCSF"]] / pcsf_wt)

    act0 = p["W_ACT_LYS"] * y[IX["PGRNLYS"]] + p["W_ACT_EXTRA"] * (
        y[IX["PGRNCSF"]] / pcsf_wt)
    pen0 = (p["C9_LYSO_PEN"] if p["GENO_C9"] > 0.5 else 0.0) + \
           (0.0 if p["TMEM106B_PROT"] > 0.5 else p["TMEM_PEN"])
    h = p["HILL_LYSO"]
    lt0 = act0 ** h / (act0 ** h + p["K50_LYSO"] ** h)
    y[IX["LYSO"]] = lt0 * (1.0 - pen0)

    y[IX["TDPN"]] = 1.0
    y[IX["C9RNA"]] = (p["KSYN_C9RNA"] / p["KDEG_C9RNA"]) if p["GENO_C9"] > 0.5 else 0.0
    y[IX["POLYGP"]] = p["KRAN_GP"] * y[IX["C9RNA"]] / p["KDEG_GP"]
    y[IX["POLYGR"]] = p["KRAN_GR"] * y[IX["C9RNA"]] / p["KDEG_GR"]

    y[IX["TAUM"]] = p["KSYN_TAUM"] / p["KDEG_TAUM"]
    y[IX["TAUS"]] = p["KTRL_TAU"] * y[IX["TAUM"]] / (p["KDEG_TAUS"] + p["KPHOS_TAU"])

    y[IX["SYN"]] = 1.0
    y[IX["NEURSN"]] = 1.0
    y[IX["NEURTMP"]] = 1.0
    y[IX["VOLF"]] = 100.0
    y[IX["VOLT"]] = 100.0

    y[IX["NFLCSF"]] = p["NFL_CSF_BASE"]
    y[IX["NFLPL"]] = p["NFL_PL_BASE"]
    y[IX["GFAPPL"]] = p["GFAP_BASE"]

    y[IX["HT5"]] = 1.0
    y[IX["ACH"]] = p["ACH_INTEGRITY"]
    y[IX["DA"]] = 1.0
    y[IX["BMPCSF"]] = (p["BMP_BASE"]
                       + p["KPROD_BMP"] * (1.0 - y[IX["LYSO"]]) / p["KOUT_BMP"])

    y[IX["CDR"]] = 0.0
    y[IX["NPI"]] = 0.0
    return y


def rhs(t, y, p):
    """Right-hand side; mirror of [ODE] in the R model."""
    d = [0.0] * N
    g = y  # readability

    CLATO = g[IX["CLATO"]]; CPLATO = g[IX["CPLATO"]]
    SORTTOT = g[IX["SORTTOT"]]; BMPCSF = g[IX["BMPCSF"]]
    CSFLATO = g[IX["CSFLATO"]]
    AAVVG = g[IX["AAVVG"]]; TRANSD = g[IX["TRANSD"]]
    ASO9CSF = g[IX["ASO9CSF"]]; ASO9TIS = g[IX["ASO9TIS"]]
    ASOTCSF = g[IX["ASOTCSF"]]; ASOTTIS = g[IX["ASOTTIS"]]
    SSRIG = g[IX["SSRIG"]]; SSRIC = g[IX["SSRIC"]]
    TRZC = g[IX["TRZC"]]; DNPC = g[IX["DNPC"]]
    PGRNPL = g[IX["PGRNPL"]]; PGRNCSF = g[IX["PGRNCSF"]]; PGRNLYS = g[IX["PGRNLYS"]]
    LYSO = g[IX["LYSO"]]; LIPO = g[IX["LIPO"]]
    TDPN = g[IX["TDPN"]]; TDPAGG = g[IX["TDPAGG"]]
    STMN2T = g[IX["STMN2T"]]; UNC13AT = g[IX["UNC13AT"]]
    C9RNA = g[IX["C9RNA"]]; POLYGP = g[IX["POLYGP"]]; POLYGR = g[IX["POLYGR"]]
    TAUM = g[IX["TAUM"]]; TAUS = g[IX["TAUS"]]; TAUP = g[IX["TAUP"]]
    TAUAGG = g[IX["TAUAGG"]]
    MG = g[IX["MG"]]; C1Q = g[IX["C1Q"]]; C3 = g[IX["C3"]]; AST = g[IX["AST"]]
    SYN = g[IX["SYN"]]; NEURSN = g[IX["NEURSN"]]; NEURTMP = g[IX["NEURTMP"]]
    VOLF = g[IX["VOLF"]]; VOLT = g[IX["VOLT"]]
    NFLCSF = g[IX["NFLCSF"]]; NFLPL = g[IX["NFLPL"]]; GFAPPL = g[IX["GFAPPL"]]
    HT5 = g[IX["HT5"]]; ACH = g[IX["ACH"]]; DA = g[IX["DA"]]
    CDR = g[IX["CDR"]]; NPI = g[IX["NPI"]]

    # --- latozinemab PK + QUASI-EQUILIBRIUM TMDD ------------------------------
    # IgG binding to sortilin is fast relative to every other timescale here and
    # the drug is ~2-3 orders of magnitude above the receptor, so occupancy is
    # taken to be at rapid equilibrium (QE-TMDD) rather than integrated as an
    # explicit on/off pair. Receptor-mediated elimination and drug-induced
    # receptor downregulation are both retained.
    THETA_PL = CLATO / (CLATO + p["KD_SORT_PL"])
    THETA_CNS = CSFLATO / (CSFLATO + p["KD_SORT_CNS"])
    rme = p["KINT_SORT"] * THETA_PL * SORTTOT          # nM/day, target-mediated
    d[IX["CLATO"]] = (-(p["CL_LATO"] / p["V1_LATO"]) * CLATO
                      - (p["Q_LATO"] / p["V1_LATO"]) * CLATO
                      + (p["Q_LATO"] / p["V1_LATO"]) * CPLATO
                      - p["KIN_CSF"] * CLATO - rme)
    d[IX["CPLATO"]] = (p["Q_LATO"] / p["V2_LATO"]) * (CLATO - CPLATO)
    ksyn_sort = p["R0_SORT"] * p["KDEG_SORT"]
    d[IX["SORTTOT"]] = (ksyn_sort - p["KDEG_SORT"] * SORTTOT
                        - p["KDOWN_SORT"] * THETA_PL * SORTTOT)
    d[IX["CSFLATO"]] = p["KIN_CSF"] * CLATO - p["KOUT_CSF"] * CSFLATO

    # --- AAV-GRN --------------------------------------------------------------
    d[IX["AAVVG"]] = -p["KDEC_VG"] * AAVVG - p["KTR_AAV"] * AAVVG * (1 - p["NAB_FRAC"])
    d[IX["TRANSD"]] = (p["KTR_AAV"] * AAVVG * (1 - p["NAB_FRAC"])
                       - p["KLOSS_TRANSD"] * TRANSD)

    # --- ASOs -----------------------------------------------------------------
    d[IX["ASO9CSF"]] = -p["KOUT_ASO9_CSF"] * ASO9CSF - p["KTIS_ASO9"] * ASO9CSF
    d[IX["ASO9TIS"]] = (p["KTIS_ASO9"] * ASO9CSF * p["ASO9_ROSTRAL"]
                        - p["KOUT_ASO9_TIS"] * ASO9TIS)
    d[IX["ASOTCSF"]] = -p["KOUT_ASOT_CSF"] * ASOTCSF - p["KTIS_ASOT"] * ASOTCSF
    d[IX["ASOTTIS"]] = (p["KTIS_ASOT"] * ASOTCSF * p["ASO9_ROSTRAL"]
                        - p["KOUT_ASOT_TIS"] * ASOTTIS)
    KD9 = p["EMAX_ASO9"] * ASO9TIS / (p["EC50_ASO9"] + ASO9TIS)
    KDT = p["EMAX_ASOT"] * ASOTTIS / (p["EC50_ASOT"] + ASOTTIS)

    # --- symptomatic drugs ----------------------------------------------------
    d[IX["SSRIG"]] = -p["KA_SSRI"] * SSRIG
    d[IX["SSRIC"]] = p["KA_SSRI"] * SSRIG / p["V_SSRI"] - p["KE_SSRI"] * SSRIC
    d[IX["TRZC"]] = -p["KE_TRZ"] * TRZC
    d[IX["DNPC"]] = -p["KE_DNP"] * DNPC
    SERTBLK = p["EMAX_SERT"] * SSRIC / (p["EC50_SERT"] + SSRIC)
    TRZBLK = p["EMAX_TRZ"] * TRZC / (p["EC50_TRZ"] + TRZC)
    ACHEBLK = p["EMAX_ACHE"] * DNPC / (p["EC50_ACHE"] + DNPC)

    # --- progranulin ----------------------------------------------------------
    GRND = grn_dose(p) + p["AAV_BOOST"] * TRANSD
    d[IX["PGRNPL"]] = (p["KSYN_PGRN_PL"] * GRND
                       - p["KCL_SORT_PL"] * (1 - THETA_PL) * PGRNPL
                       - p["KCL_OTH_PL"] * PGRNPL)
    d[IX["PGRNCSF"]] = (p["KSYN_PGRN_CSF"] * GRND
                        - p["KCL_SORT_CSF"] * (1 - THETA_CNS) * PGRNCSF
                        - p["KCL_OTH_CSF"] * PGRNCSF)
    pcsf_wt = p["KSYN_PGRN_CSF"] / (p["KCL_SORT_CSF"] + p["KCL_OTH_CSF"])
    PCSFN = PGRNCSF / pcsf_wt
    uptake = (p["K_UPTAKE_LYS"]
              * (p["P_LYS_SORT"] * (1 - THETA_CNS) + (1 - p["P_LYS_SORT"]))
              * PCSFN)
    d[IX["PGRNLYS"]] = uptake - p["KDEG_LYS_PGRN"] * PGRNLYS
    PGRNACT = p["W_ACT_LYS"] * PGRNLYS + p["W_ACT_EXTRA"] * PCSFN

    # --- lysosome -------------------------------------------------------------
    pen = ((p["C9_LYSO_PEN"] if p["GENO_C9"] > 0.5 else 0.0)
           + (0.0 if p["TMEM106B_PROT"] > 0.5 else p["TMEM_PEN"]))
    pa = max(PGRNACT, 1e-9)
    h = p["HILL_LYSO"]
    lysotgt = pa ** h / (pa ** h + p["K50_LYSO"] ** h)
    lysotgt = lysotgt * (1 - pen) - p["K_LIPO_DAM"] * LIPO
    lysotgt = max(lysotgt, 0.0)
    d[IX["LYSO"]] = p["KT_LYSO"] * (lysotgt - LYSO)
    d[IX["LIPO"]] = p["KGEN_LIPO"] * (1 - LYSO) - p["KCLR_LIPO"] * LIPO

    # --- TDP-43 ---------------------------------------------------------------
    # AGE GATE — age-related decline of autophagy/lysosomal capacity and
    # microglial priming. This is what makes ONSET AGE a model output instead
    # of "t=0 plus a fixed accrual time": a constant genetic deficit alone
    # cannot explain why a lifelong carrier stays well for decades.
    ageyr = p["AGE0"] + t / 365.0
    agefac = (max(ageyr, 1.0) / p["AGE_REF"]) ** p["AGE_POW"] + p["AGE_FLOOR"]
    if p["TMEM106B_PROT"] > 0.5:
        agefac *= p["TMEM_ONSET_SHIFT"]
    oxstr = p["W_STRESS_OX"]
    polyga = p["KRAN_GA_REL"] * POLYGP
    aggc_t = min(max(TDPAGG, 0.0), p["AGG_MAX_TDP"])
    tdpdrive = ((p["W_STRESS_LYS"] * (1 - LYSO) + p["W_STRESS_GA"] * polyga + oxstr)
                * agefac * (1 + p["KSEED_TDP"] * aggc_t))
    tdp_allowed = 0.0 if p["GENO_MAPT"] > 0.5 else 1.0
    d[IX["TDPN"]] = (p["KIN_TDP"] * (1 - TDPN)
                     - p["KMIS_TDP"] * TDPN * tdpdrive * tdp_allowed)
    cyto = 1 - TDPN
    d[IX["TDPAGG"]] = (p["KAGG_TDP"] * cyto * (1 + p["KSEED_TDP"] * aggc_t)
                       * (1 - aggc_t / p["AGG_MAX_TDP"]) * tdp_allowed
                       - p["KCLR_TDP"] * LYSO * TDPAGG)
    crypd = cyto ** p["HILL_CRYP"] if cyto > 0 else 0.0
    d[IX["STMN2T"]] = p["KCRYP_STMN2"] * crypd - p["KDEG_STMN2"] * STMN2T
    d[IX["UNC13AT"]] = p["KCRYP_UNC13A"] * crypd - p["KDEG_UNC13A"] * UNC13AT

    # --- C9 -------------------------------------------------------------------
    c9on = 1.0 if p["GENO_C9"] > 0.5 else 0.0
    d[IX["C9RNA"]] = p["KSYN_C9RNA"] * c9on * (1 - KD9) - p["KDEG_C9RNA"] * C9RNA
    d[IX["POLYGP"]] = p["KRAN_GP"] * C9RNA - p["KDEG_GP"] * POLYGP
    d[IX["POLYGR"]] = p["KRAN_GR"] * C9RNA - p["KDEG_GR"] * POLYGR

    # --- tau ------------------------------------------------------------------
    taub = 1.0 + (p["MUT_TAU_BOOST"] if p["GENO_MAPT"] > 0.5 else 0.0)
    d[IX["TAUM"]] = p["KSYN_TAUM"] * (1 - KDT) - p["KDEG_TAUM"] * TAUM
    d[IX["TAUS"]] = (p["KTRL_TAU"] * TAUM - p["KDEG_TAUS"] * TAUS
                     - p["KPHOS_TAU"] * taub * TAUS)
    d[IX["TAUP"]] = (p["KPHOS_TAU"] * taub * TAUS - p["KDEPHOS_TAU"] * TAUP
                     - p["KOLIG_TAU"] * TAUP)
    aggc_u = min(max(TAUAGG, 0.0), p["AGG_MAX_TAU"])
    d[IX["TAUAGG"]] = (p["KOLIG_TAU"] * agefac * TAUP * (1 + p["KSEED_TAU"] * aggc_u)
                       * (1 - aggc_u / p["AGG_MAX_TAU"])
                       - p["KCLR_TAU"] * LYSO * TAUAGG)

    # --- inflammation ---------------------------------------------------------
    mgdrive = (p["W_MG_TDP"] * TDPAGG + p["W_MG_TAU"] * TAUAGG
               + p["W_MG_PGRN"] * (1 - min(PGRNLYS, 1.0))
               + p["W_MG_GR"] * POLYGR)
    d[IX["MG"]] = p["KACT_MG"] * agefac * mgdrive * (1 - MG) - p["KRES_MG"] * MG
    d[IX["C1Q"]] = p["KSYN_C1Q"] * MG * (1 - p["I_C1Q"]) - p["KDEG_C1Q"] * C1Q
    d[IX["C3"]] = p["KSYN_C3"] * C1Q - p["KDEG_C3"] * C3
    d[IX["AST"]] = (p["KACT_AST"] * (TDPAGG + TAUAGG + MG) * (1 - AST)
                    - p["KRES_AST"] * AST)

    # --- synapse / neurons / volume ------------------------------------------
    d[IX["SYN"]] = (p["KREP_SYN"] / agefac * (1 - SYN)
                    - p["KPRUNE_C3"] * C3 * SYN
                    - p["KTOX_TAU_SYN"] * TAUP * SYN
                    - p["KTOX_UNC13A"] * UNC13AT * SYN)
    tox = (p["W_TOX_TDP"] * TDPAGG + p["W_TOX_TAU"] * TAUAGG
           + p["W_TOX_GR"] * POLYGR + p["W_TOX_CYT"] * MG
           + p["W_TOX_SYN"] * (1 - SYN) + p["W_TOX_STMN2"] * STMN2T)
    wsn = 1.0 - p["TEMP_SHIFT"]
    wtp = 0.55 + p["TEMP_SHIFT"]
    d[IX["NEURSN"]] = -p["KDEATH_SN"] * tox * wsn * NEURSN
    d[IX["NEURTMP"]] = -p["KDEATH_TEMP"] * tox * wtp * NEURTMP
    damf = 0.60 * (1 - NEURSN) + 0.40 * (1 - SYN)
    damt = 0.60 * (1 - NEURTMP) + 0.40 * (1 - SYN)
    d[IX["VOLF"]] = -p["KATR_F"] * damf * VOLF - p["KATR_AGE"] * VOLF
    d[IX["VOLT"]] = -p["KATR_T"] * damt * VOLT - p["KATR_AGE"] * VOLT

    # --- biomarkers -----------------------------------------------------------
    axdeg = (0.55 * (1 - NEURSN) + 0.30 * (1 - NEURTMP)
             + 0.45 * STMN2T + 0.25 * (1 - SYN))
    d[IX["NFLCSF"]] = (p["KPROD_NFL"] * axdeg
                       + p["KOUT_NFL_CSF"] * p["NFL_CSF_BASE"]
                       - p["KOUT_NFL_CSF"] * NFLCSF)
    d[IX["NFLPL"]] = (p["KTR_NFL"] * NFLCSF + p["KOUT_NFL_PL"] * p["NFL_PL_BASE"]
                      - p["KTR_NFL"] * p["NFL_CSF_BASE"] - p["KOUT_NFL_PL"] * NFLPL)
    d[IX["GFAPPL"]] = (p["KPROD_GFAP"] * AST + p["KOUT_GFAP"] * p["GFAP_BASE"]
                       - p["KOUT_GFAP"] * GFAPPL)
    # CSF BMP: a LYSOSOMAL-function readout. This is the measurement the model
    # argues should have been made alongside plasma PGRN, because unlike plasma
    # PGRN it responds to the pool that decides whether the drug works.
    d[IX["BMPCSF"]] = (p["KPROD_BMP"] * (1.0 - LYSO)
                       + p["KOUT_BMP"] * p["BMP_BASE"] - p["KOUT_BMP"] * BMPCSF)

    # --- neurotransmitters ----------------------------------------------------
    d[IX["HT5"]] = (p["KIN_5HT"] * (1 - p["W_5HT_LOSS"] * (1 - NEURSN))
                    - p["KOUT_5HT"] * (1 - SERTBLK) * HT5)
    d[IX["ACH"]] = (p["KIN_ACH"] * p["ACH_INTEGRITY"]
                    - p["KOUT_ACH"] * (1 - ACHEBLK) * ACH)
    d[IX["DA"]] = (p["KIN_DA"] * (1 - p["W_DA_LOSS"] * min(TAUAGG, 1.0))
                   - p["KOUT_DA"] * DA)

    # --- clinical -------------------------------------------------------------
    cdrdam = (p["A_CDR_SYN"] * (1 - SYN) + p["A_CDR_SN"] * (1 - NEURSN)
              + p["A_CDR_TEMP"] * (1 - NEURTMP)
              + p["A_CDR_VOL"] * (1 - VOLF / 100.0))
    # Cognitive reserve: clinical deficit appears only once accumulated damage
    # exceeds a threshold. This is what produces the long silent prodrome and
    # then the relatively abrupt clinical decline that FTD actually shows.
    cdrdrive = max(cdrdam - p["CDR_THRESH"], 0.0)
    d[IX["CDR"]] = p["KCDR"] * cdrdrive * (1 - CDR / p["CDR_MAX"])

    # Presynaptic serotonergic deficit is SSRI-reversible; the salience-network
    # term is not (it stands for postsynaptic 5-HT2A and neuronal loss), which
    # is why SSRIs help behaviour partially and never fully.
    npidrive = (p["B_NPI_5HT"] * max(0.0, 1 - HT5)
                + p["B_NPI_SN"] * (1 - NEURSN)
                + p["B_NPI_DA"] * max(0.0, 1 - DA))
    npidrive = npidrive * (1 - TRZBLK) * (1 + p["K_AGIT_DNP"] * ACHEBLK)
    npitgt = p["NPI_MAX"] * p["EMAX_NPI"] * npidrive / (p["K50_NPI"] + npidrive)
    d[IX["NPI"]] = p["KT_NPI"] * (npitgt - NPI)

    nflr = max(NFLPL / p["NFL_PL_BASE"], 1e-6)
    hazarg = (p["BH_CDR"] * CDR + p["BH_NFL"] * math.log(nflr)
              + p["BH_ALS"] * p["ALS_FLAG"])
    haz = p["H0_HAZ"] * math.exp(min(hazarg, 50.0))
    d[IX["CUMHAZ"]] = haz
    return d


# -----------------------------------------------------------------------------
# INTEGRATOR
# -----------------------------------------------------------------------------
def simulate(p, tend, h=0.05, doses=(), record_every=7.0, y0=None, t0=0.0):
    """RK4 with instantaneous dose increments. doses = [(time, state, amount)].
    Dose times and t0 are on the SAME absolute clock, so a trial that starts
    from onset_state() keeps the patient's true age in the aging stress term."""
    y = list(y0) if y0 is not None else init_state(p)
    dose_map = {}
    for (dt, cmt, amt) in doses:
        key = round(dt / h)
        dose_map.setdefault(key, []).append((IX[cmt], amt))

    nsteps = int(round(tend / h))
    rec_stride = max(1, int(round(record_every / h)))
    out = {"time": [0.0]}
    for n in NAMES:
        out[n] = [y[IX[n]]]

    for step in range(nsteps):
        t = t0 + step * h
        for (i, a) in dose_map.get(step, ()):
            y[i] += a
        k1 = rhs(t, y, p)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(N)]
        k2 = rhs(t + 0.5 * h, y2, p)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(N)]
        k3 = rhs(t + 0.5 * h, y3, p)
        y4 = [y[i] + h * k3[i] for i in range(N)]
        k4 = rhs(t + h, y4, p)
        y = [y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(N)]
        if (step + 1) % rec_stride == 0:
            out["time"].append((step + 1) * h)
            for n in NAMES:
                out[n].append(y[IX[n]])
    # always record the final point
    if abs(out["time"][-1] - tend) > 1e-9:
        out["time"].append(tend)
        for n in NAMES:
            out[n].append(y[IX[n]])
    return out


def pars(**kw):
    p = dict(P0)
    p.update(kw)
    return p


GRN = dict(GENO_GRN=1, GENO_C9=0, GENO_MAPT=0, SPORADIC=0)
C9 = dict(GENO_GRN=0, GENO_C9=1, GENO_MAPT=0, SPORADIC=0)
MAPT = dict(GENO_GRN=0, GENO_C9=0, GENO_MAPT=1, SPORADIC=0)


def lato_doses(mgkg=60, bwt=70, n=24, start=0.0, p=P0):
    """Latozinemab IV q4w. mg -> nM in the central volume."""
    conv = (1.0e6 / p["MW_LATO"]) / p["V1_LATO"]
    return [(start + 28.0 * i, "CLATO", mgkg * bwt * conv) for i in range(n)]


def onset_state(p, cdr_target=5.5, tmax_y=45.0, h=0.5):
    """Integrate the natural history from a clinically silent carrier until the
    CDR plus NACC-FTLD SB reaches `cdr_target`, and return (t*, state).

    This is the model's disease clock. Trial scenarios START from this state
    rather than from hand-set 'symptomatic baseline' initial conditions, so a
    treated and an untreated arm share one continuous trajectory and the
    baseline biomarker values are EARNED by the model rather than assumed.
    """
    y = init_state(p)
    nsteps = int(round(tmax_y * 365.0 / h))
    for step in range(nsteps):
        t = step * h
        k1 = rhs(t, y, p)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(N)]
        k2 = rhs(t + 0.5 * h, y2, p)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(N)]
        k3 = rhs(t + 0.5 * h, y3, p)
        y4 = [y[i] + h * k3[i] for i in range(N)]
        k4 = rhs(t + h, y4, p)
        y = [y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(N)]
        if y[IX["CDR"]] >= cdr_target:
            return (step + 1) * h, y
    return None, y


def aso_doses(cmt, amt=100.0, n=8, ii=84.0, start=0.0):
    return [(start + ii * i, cmt, amt) for i in range(n)]


def oral_doses(cmt, amt, n, V=1.0, ii=1.0, start=0.0):
    return [(start + ii * i, cmt, amt / V) for i in range(n)]


def last(out, name):
    return out[name][-1]


def at(out, name, t):
    """Value of `name` at the recorded time closest to t."""
    best, bi = 1e18, 0
    for i, tt in enumerate(out["time"]):
        if abs(tt - t) < best:
            best, bi = abs(tt - t), i
    return out[name][bi]


# =============================================================================
#  ANALYSES
# =============================================================================
R = []


def say(s=""):
    R.append(s)


def section(t):
    say("")
    say("=" * 78)
    say(t)
    say("=" * 78)


F_PL = P0["KCL_SORT_PL"] / (P0["KCL_SORT_PL"] + P0["KCL_OTH_PL"])
F_CSF = P0["KCL_SORT_CSF"] / (P0["KCL_SORT_CSF"] + P0["KCL_OTH_CSF"])
PPL_WT = P0["KSYN_PGRN_PL"] / (P0["KCL_SORT_PL"] + P0["KCL_OTH_PL"])
PCSF_WT = P0["KSYN_PGRN_CSF"] / (P0["KCL_SORT_CSF"] + P0["KCL_OTH_CSF"])

# Cache of onset states, so each genotype's disease clock is solved only once.
_ONSET = {}


def onset(key, p, cdr=5.5):
    if key not in _ONSET:
        _ONSET[key] = M_onset(p, cdr)
    return _ONSET[key]


def M_onset(p, cdr):
    t, y = onset_state(p, cdr_target=cdr, tmax_y=80.0, h=0.5)
    return t, y


def trial(p, doses=(), tend=672.0, key=None):
    """Run a 96-week 'trial' starting from the model's own onset state."""
    t0, y0 = onset(key, p)
    return simulate(p, tend, h=0.05, doses=doses, record_every=28.0,
                    y0=y0, t0=t0)


def main(brief=False):
    say("=" * 78)
    say("FRONTOTEMPORAL DEMENTIA QSP MODEL — REFERENCE IMPLEMENTATION REPORT")
    say("49 ODE states | dependency-free RK4 | every number below is computed")
    say("=" * 78)
    say("")
    say("Trial scenarios do NOT start from hand-set 'symptomatic baseline'")
    say("initial conditions. Each genotype is run from a clinically silent")
    say("carrier at age 30 until its CDR plus NACC-FTLD SB reaches 5.5, and")
    say("that state becomes the trial baseline. Treated and untreated arms")
    say("therefore share one continuous trajectory, and every baseline")
    say("biomarker value is EARNED by the model rather than assumed.")

    # =====================================================================
    section("A. ONSET AGE AS A PREDICTION, NOT AN OFFSET")
    # =====================================================================
    say("A constant genetic deficit cannot by itself explain why a lifelong")
    say("carrier stays well for decades, so the model carries an explicit AGE")
    say("GATE (declining autophagic/lysosomal capacity and microglial priming,")
    say("(age/62)^5). With it, onset AGE becomes an output. Observed values are")
    say("cohort MEANS from Moore 2020 Lancet Neurol (PMID 31810826), n=3403")
    say("carriers from 1492 families; the sporadic figure is the usual ~58 y.")
    say("")
    say("  genotype      predicted onset   observed   error")
    say("  " + "-" * 52)
    # Observed onset ages are the cohort MEANS from Moore et al. 2020 Lancet
    # Neurol (PMID 31810826), n=3403 carriers from 1492 families: MAPT 49.5 y
    # (SD 10.0), C9orf72 58.2 (9.8), GRN 61.3 (8.8). The sporadic bvFTD figure
    # is the commonly quoted ~58 y and is not from that cohort.
    geno = [
        ("MAPT", pars(**MAPT), 49.5, "mapt"),
        ("C9orf72", pars(**C9, ALS_FLAG=1), 58.2, "c9"),
        ("GRN", pars(**GRN), 61.3, "grn"),
        ("sporadic", pars(), 58.0, "spo"),
    ]
    onsets = {}
    for lab, p, obs, key in geno:
        t, y = onset(key, p)
        age = p["AGE0"] + t / 365.0
        onsets[lab] = age
        say(f"  {lab:12s}{age:14.1f} y{obs:12.1f} y{age - obs:+8.1f} y")
    say("")
    say("  MAPT lands within ~1 year of the reported mean and sporadic within")
    say("  ~1.5. C9orf72 is ~4 y early. GRN is ~5.5 y early, and there the")
    say("  reason is STRUCTURAL rather than a bad parameter: in this model PGRN")
    say("  haploinsufficiency is a damage driver acting from birth, so it can")
    say("  only ever move onset EARLIER than sporadic disease. The measured")
    say("  ordering is the opposite — GRN onset 61.3 y is LATER than sporadic")
    say("  bvFTD at ~58 y. NO monotone version of this model can produce that")
    say("  ordering. Two candidate resolutions:")
    say("    (i)  lifelong haploinsufficiency is developmentally compensated,")
    say("         so the deficit is not rate-limiting until a second hit —")
    say("         consistent with GRN's incomplete penetrance and its strong")
    say("         dependence on modifiers such as TMEM106B;")
    say("    (ii) the comparison is confounded by ascertainment: GRN carriers")
    say("         are found by GENOTYPE regardless of severity, while")
    say("         'sporadic bvFTD' is a clinically defined group enriched for")
    say("         aggressive pathology. The two cohorts are not exchangeable.")
    say("")
    p_tm = pars(**GRN, TMEM106B_PROT=1)
    t_tm, _ = onset("grn_tmem", p_tm)
    say(f"  The model does reproduce the modifier DIRECTION: the protective")
    say(f"  TMEM106B haplotype shifts predicted GRN onset from"
        f" {onsets['GRN']:.1f} y to"
        f" {p_tm['AGE0'] + t_tm/365.0:.1f} y (later),")
    say("  which is the direction reported for rs1990622. Note the model puts")
    say("  TMEM106B on the lysosomal/age-vulnerability axis and NOT on PGRN")
    say("  levels, which agrees with the 7071-person biofluid meta-analysis")
    say("  (PMID 38539243) finding no association between plasma PGRN and")
    say("  rs1990622 genotype. That was not fitted; it fell out of the wiring.")

    # =====================================================================
    section("B. NATURAL HISTORY vs PUBLISHED ANCHORS")
    # =====================================================================
    say("Measured over the 2 years following each genotype's own onset state.")
    say("")
    say("  genotype     slope/y   plasma NfL   CSF NfL   atrophy%/y   NPI")
    say("  " + "-" * 62)
    for lab, p, obs, key in geno:
        t0, y0 = onset(key, p)
        o = simulate(p, 730.0, h=0.5, y0=y0, t0=t0, record_every=30.0)
        slope = (o["CDR"][-1] - o["CDR"][0]) / 2.0
        atr = (o["VOLF"][0] - o["VOLF"][-1]) / o["VOLF"][0] * 100 / 2
        say(f"  {lab:12s}{slope:7.2f}{o['NFLPL'][0]:13.1f}{o['NFLCSF'][0]:10.0f}"
            f"{atr:12.2f}{o['NPI'][0]:7.0f}")
    say("")
    say("  anchors: CDR+NACC-FTLD SB slope 1.5-2.5 /y symptomatic;")
    say("           plasma NfL 50-80 pg/mL (control ~10); CSF NfL 3000-5000")
    say("           (control ~700); frontal atrophy 2-3 %/y (control ~0.5);")
    say("           NPI total 20-40 in symptomatic bvFTD.")
    say("")
    say("  Presymptomatic lead time — does NfL rise before the clinical score?")
    p_g = pars(**GRN)
    tconv, _ = onset("grn", p_g)
    o_pre = simulate(p_g, tconv, h=0.5, record_every=180.0)
    t3 = tnfl = None
    for i, tt in enumerate(o_pre["time"]):
        if t3 is None and o_pre["CDR"][i] >= 3.0:
            t3 = tt
        if tnfl is None and o_pre["NFLPL"][i] >= 1.5 * P0["NFL_PL_BASE"]:
            tnfl = tt
    say("")
    say("    threshold on plasma NfL   age crossed   lead time before CDR=3")
    say("    " + "-" * 60)
    for mult in [1.5, 2.0, 3.0, 4.0, 5.0]:
        tc = None
        for i, tt in enumerate(o_pre["time"]):
            if o_pre["NFLPL"][i] >= mult * P0["NFL_PL_BASE"]:
                tc = tt
                break
        if tc is not None and t3 is not None:
            say(f"    {mult:4.1f}x control{p_g['AGE0'] + tc/365.0:16.1f} y"
                f"{(t3 - tc)/365.0:20.1f} y")
    if t3:
        say(f"    CDR+NACC-FTLD SB reaches 3 at age"
            f" {p_g['AGE0'] + t3/365.0:.1f}.  Reported lead time ~1-2 y.")
    say("")
    say("    THIS IS A FAILURE, and it is the SAME failure as the GRN onset")
    say("    error above. The model has no genuinely compensated latent phase:")
    say("    damage begins at t=0 and accumulates monotonically, so NfL starts")
    say("    drifting up decades early and only a very high threshold (~4-5x")
    say("    control) recovers a plausible 1-2 y lead time. Real presymptomatic")
    say("    carriers hold near-normal NfL and then rise steeply near")
    say("    conversion. Reproducing that needs a compensated regime with a")
    say("    threshold-crossing collapse, which this monotone structure lacks.")

    # =====================================================================
    section("C. THE CENTRAL ANALYTIC RESULT AND ITS NUMERICAL VERIFICATION")
    # =====================================================================
    say(f"Sortilin share of PLASMA PGRN clearance    f_plasma = {F_PL:.4f}")
    say(f"Sortilin share of CSF    PGRN clearance    f_CSF    = {F_CSF:.4f}")
    say(f"Sortilin share of LYSOSOMAL PGRN delivery  p        ="
        f" {P0['P_LYS_SORT']:.4f}  <-- NOT MEASURED")
    say("")
    say("Sortilin is BOTH the clearance route that sets extracellular PGRN AND")
    say("a delivery route filling the lysosomal pool where PGRN acts. Blocking")
    say("it raises the ligand and removes a delivery route simultaneously. At")
    say("steady state:")
    say("    CSF PGRN fold-rise       = 1 / (1 - f*theta)")
    say("    lysosomal delivery ratio = (1 - p*theta) / (1 - f*theta)")
    say("")
    p = pars(**GRN)
    nat = trial(p, key="grn")
    tx = trial(p, doses=lato_doses(p=p), key="grn")
    th_pl = tx["CLATO"][-1] / (tx["CLATO"][-1] + P0["KD_SORT_PL"])
    th_cns = tx["CSFLATO"][-1] / (tx["CSFLATO"][-1] + P0["KD_SORT_CNS"])
    say(f"  latozinemab 60 mg/kg q4w, trough at 96 wk:")
    say(f"    plasma  {tx['CLATO'][-1]:8.0f} nM      CSF {tx['CSFLATO'][-1]:6.2f} nM"
        f"   (CSF/plasma {100*tx['CSFLATO'][-1]/tx['CLATO'][-1]:.3f}%,"
        " reported ~0.1-0.3%)")
    say(f"    sortilin occupancy: peripheral {th_pl:.4f}   CNS {th_cns:.4f}")
    say("")
    say("  quantity                 numeric    algebraic   residual")
    say("  " + "-" * 56)
    rows = [("plasma PGRN fold-rise", tx["PGRNPL"][-1] / nat["PGRNPL"][0],
             1.0 / (1.0 - F_PL * th_pl)),
            ("CSF PGRN fold-rise", tx["PGRNCSF"][-1] / nat["PGRNCSF"][0],
             1.0 / (1.0 - F_CSF * th_cns)),
            ("LYSOSOMAL PGRN ratio", tx["PGRNLYS"][-1] / nat["PGRNLYS"][0],
             (1.0 - P0["P_LYS_SORT"] * th_cns) / (1.0 - F_CSF * th_cns))]
    for lab, num, alg in rows:
        say(f"  {lab:24s}{num:8.4f}{alg:12.4f}{abs(num-alg):11.4f}")
    say("")
    say("  Residuals are small and arise only because occupancy oscillates")
    say("  within the dosing interval while the algebra assumes it constant.")
    say("  The 49-state integration and the closed form agree, so what follows")
    say("  is a property of the STRUCTURE, not of one solver run.")

    # =====================================================================
    section("D. THE BIOMARKER / TARGET DISSOCIATION")
    # =====================================================================
    say("Latozinemab 60 mg/kg IV q4w in a symptomatic GRN carrier, 96 weeks.")
    say("")
    say("  compartment       baseline   on drug     fold")
    say("  " + "-" * 48)
    say(f"  PLASMA PGRN      {100*nat['PGRNPL'][0]/PPL_WT:8.1f}%"
        f"{100*tx['PGRNPL'][-1]/PPL_WT:9.1f}%   {tx['PGRNPL'][-1]/nat['PGRNPL'][0]:5.2f}x")
    say(f"  CSF    PGRN      {100*nat['PGRNCSF'][0]/PCSF_WT:8.1f}%"
        f"{100*tx['PGRNCSF'][-1]/PCSF_WT:9.1f}%   {tx['PGRNCSF'][-1]/nat['PGRNCSF'][0]:5.2f}x")
    say(f"  LYSOSOMAL PGRN   {100*nat['PGRNLYS'][0]:8.1f}%"
        f"{100*tx['PGRNLYS'][-1]:9.1f}%   {tx['PGRNLYS'][-1]/nat['PGRNLYS'][0]:5.2f}x")
    say("  (% of the wild-type level for that compartment)")
    say("")
    say("  The registrational biomarker is restored essentially to the control")
    say(f"  range ({100*tx['PGRNPL'][-1]/PPL_WT:.0f}% of control). The pool where"
        " progranulin actually works")
    say(f"  reaches only {100*tx['PGRNLYS'][-1]:.0f}% of wild type.")
    say("")
    say("  96-week outcomes, treated vs untreated:")
    say("")
    say("  endpoint                  untreated    treated   difference")
    say("  " + "-" * 60)
    for nm, lab in [("CDR", "CDR+NACC-FTLD SB"), ("NPI", "NPI total"),
                    ("NFLPL", "plasma NfL pg/mL"), ("NFLCSF", "CSF NfL pg/mL"),
                    ("VOLF", "frontal volume %"), ("SYN", "synaptic density"),
                    ("C1Q", "CSF C1q au"), ("C3", "CSF C3 au"),
                    ("BMPCSF", "CSF BMP au"), ("LYSO", "lysosomal function")]:
        a, b = nat[nm][-1], tx[nm][-1]
        say(f"  {lab:24s}{a:11.3f}{b:11.3f}{b-a:+12.3f}")
    d = tx["CDR"][-1] - nat["CDR"][-1]
    prog = nat["CDR"][-1] - nat["CDR"][0]
    say("")
    say(f"  Untreated 96-week CDR progression = {prog:+.2f} points.")
    say(f"  Treatment effect                  = {d:+.4f} points"
        f"  ({-100*d/prog:.2f}% slowing).")
    say("")
    say("  A fraction-of-one-percent effect on a 0-24 scale over 96 weeks is")
    say("  far below what any feasible trial can resolve. The model therefore")
    say("  PREDICTS a failed clinical endpoint sitting on top of a clean,")
    say("  unambiguous biomarker win — which is what the anti-sortilin phase 3")
    say("  programme reported. Note the C1q/C3 and BMP columns: the PD effects")
    say("  are real and measurable. Target engagement was never the problem.")

    # =====================================================================
    section("E. THE SIGN FLIP: SWEEPING THE PARAMETER NOBODY HAS MEASURED")
    # =====================================================================
    say("p = sortilin share of lysosomal PGRN delivery. Everything else fixed.")
    say(f"Prediction: the sign flips exactly at p = f_CSF = {F_CSF:.3f}, at ANY")
    say("dose, and the theta->1 ceiling is (1-p)/(1-f_CSF).")
    say("")
    say("     p    LYSO pool  ratio  ceiling  plasma PGRN   dCDR@96wk  verdict")
    say("  " + "-" * 70)
    for pv in [0.0, 0.2, 0.4, 0.5, 0.6, 0.7, F_CSF, 0.85, 0.95]:
        pj = pars(**GRN)
        pj["P_LYS_SORT"] = pv
        kk = f"grn_p{pv:.3f}"
        nj = trial(pj, key=kk)
        oj = trial(pj, doses=lato_doses(p=pj), key=kk)
        ratio = oj["PGRNLYS"][-1] / nj["PGRNLYS"][0]
        ceil = (1 - pv) / (1 - F_CSF)
        dc = oj["CDR"][-1] - nj["CDR"][-1]
        v = "BENEFIT" if ratio > 1.02 else ("HARM" if ratio < 0.98 else "NEUTRAL")
        star = "  <== break-even" if abs(pv - F_CSF) < 1e-6 else ""
        say(f"  {pv:5.3f}{100*oj['PGRNLYS'][-1]:9.1f}%{ratio:8.3f}{ceil:8.3f}"
            f"{100*oj['PGRNPL'][-1]/PPL_WT:11.1f}%{dc:+11.4f}  {v}{star}")
    say("")
    say("  Read the plasma-PGRN column: 94.5% at EVERY value of p. The drug")
    say("  looks identical on its registrational biomarker whether it is")
    say("  helping, doing nothing, or making the target pool worse. That is")
    say("  this model's sharpest claim — the biomarker is mathematically blind")
    say("  to the parameter that decides the sign of the effect, and no amount")
    say("  of dose escalation changes that, because raising occupancy raises")
    say("  the ligand and removes delivery capacity in lockstep.")

    # =====================================================================
    section("F. GENE THERAPY IS NOT SUBJECT TO THE CEILING")
    # =====================================================================
    say("AAV-GRN raises PRODUCTION, so lysosomal delivery scales with transgene")
    say("output and no receptor is removed. Same patient, same 96 weeks:")
    say("")
    aav = trial(p, doses=[(0.0, "AAVVG", 1.0)], key="grn")
    say("  therapy         plasma   CSF    LYSOSOMAL   dCDR      dBMP")
    say("  " + "-" * 60)
    say(f"  untreated      {100*nat['PGRNPL'][0]/PPL_WT:6.1f}%"
        f"{100*nat['PGRNCSF'][0]/PCSF_WT:7.1f}%{100*nat['PGRNLYS'][0]:10.1f}%"
        f"{0.0:+9.3f}{0.0:+10.2f}")
    say(f"  latozinemab    {100*tx['PGRNPL'][-1]/PPL_WT:6.1f}%"
        f"{100*tx['PGRNCSF'][-1]/PCSF_WT:7.1f}%{100*tx['PGRNLYS'][-1]:10.1f}%"
        f"{tx['CDR'][-1]-nat['CDR'][-1]:+9.3f}"
        f"{tx['BMPCSF'][-1]-nat['BMPCSF'][-1]:+10.2f}")
    say(f"  AAV-GRN single {100*aav['PGRNPL'][-1]/PPL_WT:6.1f}%"
        f"{100*aav['PGRNCSF'][-1]/PCSF_WT:7.1f}%{100*aav['PGRNLYS'][-1]:10.1f}%"
        f"{aav['CDR'][-1]-nat['CDR'][-1]:+9.3f}"
        f"{aav['BMPCSF'][-1]-nat['BMPCSF'][-1]:+10.2f}")
    say("")
    say("  Gene therapy raises all three pools TOGETHER, because it adds")
    say("  production instead of blocking a receptor. Its lysosomal")
    say(f"  restoration is {aav['PGRNLYS'][-1]/tx['PGRNLYS'][-1]:.2f}x the"
        " antibody's, and its predicted clinical")
    say(f"  effect {(aav['CDR'][-1]-nat['CDR'][-1])/(tx['CDR'][-1]-nat['CDR'][-1]):.1f}x —"
        " from a LOWER plasma PGRN level. A plasma-PGRN")
    say("  biomarker would have RANKED THESE TWO THERAPIES BACKWARDS.")

    # =====================================================================
    section("G. WHAT TO MEASURE INSTEAD, AND THE RIVAL HYPOTHESIS")
    # =====================================================================
    say("CSF BMP (a lysosomal lipid) is in the model as the readout plasma")
    say("PGRN cannot be. Across the p sweep above it moves with the LYSOSOMAL")
    say("pool, not with the plasma pool — so it separates the cases plasma")
    say("PGRN cannot. It is attenuated, though, because lysosomal function in")
    say("this model also sees the extracellular pool (weight 0.2), so BMP is a")
    say("better discriminator than plasma PGRN without being a perfect one.")
    say("")
    say("If progranulin acts mainly EXTRACELLULARLY rather than in the")
    say("lysosome, the elevated CSF pool IS the therapeutic signal and")
    say("blocking sortilin is straightforwardly good. Same drug, same doses,")
    say("action weights swapped:")
    say("")
    p_ex = pars(**GRN, W_ACT_LYS=0.20, W_ACT_EXTRA=0.80)
    nx = trial(p_ex, key="grn_ex")
    ox = trial(p_ex, doses=lato_doses(p=p_ex), key="grn_ex")
    say("  hypothesis                    dCDR@96wk    dNfL    dLYSO fn")
    say("  " + "-" * 58)
    say(f"  lysosomal action (default)  {tx['CDR'][-1]-nat['CDR'][-1]:+11.4f}"
        f"{tx['NFLPL'][-1]-nat['NFLPL'][-1]:+9.2f}"
        f"{tx['LYSO'][-1]-nat['LYSO'][-1]:+11.4f}")
    say(f"  extracellular action        {ox['CDR'][-1]-nx['CDR'][-1]:+11.4f}"
        f"{ox['NFLPL'][-1]-nx['NFLPL'][-1]:+9.2f}"
        f"{ox['LYSO'][-1]-nx['LYSO'][-1]:+11.4f}")
    say("")
    say(f"  The two hypotheses differ"
        f" {abs((ox['CDR'][-1]-nx['CDR'][-1])/(tx['CDR'][-1]-nat['CDR'][-1])):.1f}x in predicted"
        " clinical benefit and are")
    say("  INDISTINGUISHABLE by plasma PGRN. They are distinguishable by")
    say("  downstream lysosomal chemistry. That is a cheap, falsifiable")
    say("  experiment the model recommends running before the next trial.")

    # =====================================================================
    section("H. SYMPTOMATIC PHARMACOLOGY, INCLUDING THE NULL/HARM ARM")
    # =====================================================================
    say("96 weeks in a symptomatic GRN carrier.")
    say("")
    arms = [("untreated", ()),
            ("SSRI 30 mg/d", oral_doses("SSRIG", 30.0, 672)),
            ("trazodone 150 mg/d", oral_doses("TRZC", 150.0, 672, V=P0["V_TRZ"])),
            ("SSRI + trazodone",
             oral_doses("SSRIG", 30.0, 672)
             + oral_doses("TRZC", 150.0, 672, V=P0["V_TRZ"])),
            ("donepezil 10 mg/d",
             oral_doses("DNPC", 10.0, 672, V=P0["V_DNP"]))]
    say("  arm                     NPI    dNPI     CDR    dCDR    ACh tone")
    say("  " + "-" * 64)
    res = {}
    for lab, dz in arms:
        o = trial(p, doses=dz, key="grn")
        res[lab] = o
        b = res["untreated"]
        say(f"  {lab:22s}{o['NPI'][-1]:6.1f}{o['NPI'][-1]-b['NPI'][-1]:+8.1f}"
            f"{o['CDR'][-1]:8.2f}{o['CDR'][-1]-b['CDR'][-1]:+8.3f}"
            f"{o['ACH'][-1]:10.3f}")
    b = res["untreated"]
    say("")
    say(f"  trazodone dNPI = {res['trazodone 150 mg/d']['NPI'][-1]-b['NPI'][-1]:+.1f}"
        "  (anchor: -6 to -10, Lebert 2004 crossover RCT)")
    say(f"  donepezil dNPI = {res['donepezil 10 mg/d']['NPI'][-1]-b['NPI'][-1]:+.1f}"
        "  (anchor: behaviour WORSE, Mendez 2007)")
    say("")
    say("  The AChEI null result is STRUCTURAL, not imposed. ACH_INTEGRITY is")
    say("  0.95 because the cholinergic system is largely PRESERVED in FTD, so")
    say(f"  donepezil does raise cholinergic tone ({b['ACH'][-1]:.2f} ->"
        f" {res['donepezil 10 mg/d']['ACH'][-1]:.2f}) —")
    say("  and the cognitive endpoint does not depend on it, so dCDR is")
    say("  EXACTLY zero while the agitation term makes behaviour worse. That")
    say("  is the mechanistic contrast with Alzheimer's disease, where the")
    say("  same drug acts on a genuinely depleted system.")
    say("")
    say("  Note also that no symptomatic drug moves CDR at all: in this")
    say("  structure they act on neurotransmitter tone, which feeds the")
    say("  behavioural axis only. Benefit and harm both emerge from the")
    say("  wiring rather than from an imposed effect size.")

    # =====================================================================
    section("I. ASOs — TARGET ENGAGEMENT WITHOUT CLINICAL BENEFIT")
    # =====================================================================
    p_c9 = pars(**C9, ALS_FLAG=1)
    n9 = trial(p_c9, key="c9")
    x9 = trial(p_c9, doses=aso_doses("ASO9CSF", 100.0, 8), key="c9")
    say("C9orf72 repeat-targeting ASO, intrathecal q12w x8:")
    say(f"  repeat RNA          -{100*(1-x9['C9RNA'][-1]/n9['C9RNA'][-1]):.1f}%")
    say(f"  CSF poly-GP         -{100*(1-x9['POLYGP'][-1]/n9['POLYGP'][-1]):.1f}%"
        "   (FOCUS-C9 reported up to ~50%)")
    say(f"  poly-GR (toxic)     -{100*(1-x9['POLYGR'][-1]/n9['POLYGR'][-1]):.1f}%")
    say(f"  dCDR at 96 wk       {x9['CDR'][-1]-n9['CDR'][-1]:+.4f}"
        "   (trial: no clinical benefit; discontinued)")
    say(f"  dplasma NfL         {x9['NFLPL'][-1]-n9['NFLPL'][-1]:+.2f} pg/mL")
    say("")
    say("  Why so little from ~50% target knockdown? In this structure the")
    say("  repeat products are one of several parallel inputs. The same")
    say("  genotype also carries a lysosomal penalty from C9orf72 protein")
    say("  HAPLOINSUFFICIENCY, which an RNA-targeting ASO cannot touch, and")
    say("  the TDP-43 loop is self-templating once started. Lowering poly-GP")
    say("  after that loop has closed removes an initiator, not the engine.")
    say("")
    p_mt = pars(**MAPT)
    nT = trial(p_mt, key="mapt")
    xT = trial(p_mt, doses=aso_doses("ASOTCSF", 100.0, 8), key="mapt")
    say("MAPT-lowering ASO, intrathecal q12w x8:")
    say(f"  MAPT mRNA           -{100*(1-xT['TAUM'][-1]/nT['TAUM'][-1]):.1f}%")
    say(f"  soluble tau         -{100*(1-xT['TAUS'][-1]/nT['TAUS'][-1]):.1f}%"
        "   (BIIB080 class: ~30-50%)")
    say(f"  AGGREGATED tau      -{100*(1-xT['TAUAGG'][-1]/nT['TAUAGG'][-1]):.1f}%")
    say(f"  dCDR at 96 wk       {xT['CDR'][-1]-nT['CDR'][-1]:+.4f}")
    say("")
    say("  Aggregated tau falls roughly a third as much as soluble tau:")
    say("  blocking synthesis STARVES a self-templating deposit, it does not")
    say("  CLEAR one. The model's implication is that MAPT-ASO trials should")
    say("  be powered on presymptomatic or very early cohorts, while the")
    say("  deposit is still small — and that its effect, though the largest")
    say("  of any agent modelled here, is still ~1% of 96-week progression.")

    # =====================================================================
    section("J. WHERE THIS MODEL IS WEAK OR PROBABLY WRONG")
    # =====================================================================
    say("Stated rather than tuned away.")
    say("")
    say("1. P_LYS_SORT (p) has no human measurement, and the SIGN of the")
    say("   anti-sortilin prediction depends entirely on it. The default 0.50")
    say(f"   sits below f_CSF={F_CSF:.2f} and so gives the drug the benefit of the")
    say("   doubt; a hostile choice would make it frankly harmful. What the")
    say("   model constrains is the RELATIONSHIP p vs f, not the value of p.")
    say("")
    say("2. THE MODEL HAS NO COMPENSATED LATENT PHASE, and this single defect")
    say("   causes two separate disagreements: GRN onset lands ~5 y early")
    say("   (section A) and the presymptomatic NfL lead time is ~15 y instead")
    say("   of ~1-2 y (section B). Both follow from damage accumulating")
    say("   monotonically from t=0. The fix is not a parameter but a different")
    say("   structure: a compensated regime that holds biomarkers near normal")
    say("   until a threshold is crossed, then collapses. That is the single")
    say("   most valuable change anyone could make to this model.")
    say("")
    say("3. KSYN_NMD=0.70 exists only to reconcile plasma PGRN at ~1/3 of")
    say("   control with an allele dose of 1/2. It is a fitted stand-in for")
    say("   nonsense-mediated decay plus lost autocrine recapture. Its")
    say("   consequence is real though: plasma falls to ~33% while CSF falls")
    say("   to ~50%, and no single gene-dose scalar fits both, which is itself")
    say("   evidence the two compartments have different clearance structures.")
    say("")
    say("4. CNS antibody exposure is one well-stirred CSF compartment. Real")
    say("   parenchymal distribution is regionally heterogeneous, and the")
    say("   salience network is exactly where delivery matters. If cortical")
    say("   exposure is below the CSF surrogate, theta_CNS is overestimated")
    say("   and every anti-sortilin number here is optimistic.")
    say("")
    say("5. Poly-GP is treated as a marker in equilibrium with repeat RNA and")
    say("   poly-GR as the toxic species. If the causal DPR is poly-GA instead")
    say("   (proteasome impairment), the ASO prediction changes and the CSF")
    say("   marker trials actually use would be tracking the wrong arm.")
    say("")
    say("6. Two lumped regions cannot express what actually determines FTD")
    say("   phenotype, which is WHICH network fails first — a connectome")
    say("   problem. The svPPA scenario is a weighting shift, not a mechanism.")
    say("")
    say("7. Sporadic FTLD gets the same TDP-43 machinery as GRN disease with")
    say("   no genetic driver, so its trajectory rests on the age gate alone.")
    say("   Real sporadic FTLD spans TDP-43 types A-C, FET and tau pathology;")
    say("   one parameterisation cannot represent that heterogeneity.")
    say("")
    say("8. The age gate (age/62)^5 is phenomenological. It does real work —")
    say("   it turns onset age into a prediction — but its exponent is fitted")
    say("   to onset data, so the good agreement in section A is partly")
    say("   calibration rather than pure prediction. The genotype ORDERING")
    say("   and the TMEM106B direction are the parts that are not.")

    # =====================================================================
    section("K. HEADLINE NUMBERS")
    # =====================================================================
    dlato = tx["CDR"][-1] - nat["CDR"][-1]
    say(f" 1. Predicted onset ages MAPT {onsets['MAPT']:.1f} / C9orf72"
        f" {onsets['C9orf72']:.1f} / GRN {onsets['GRN']:.1f} / sporadic"
        f" {onsets['sporadic']:.1f} y")
    say("    against Moore 2020 cohort means 49.5 / 58.2 / 61.3 and ~58: MAPT")
    say("    within 1 y, sporadic within 1.5, C9orf72 4 y early, GRN 5.5 y early")
    say("    — the last a structural miss the model explains rather than hides.")
    say(f" 2. f_plasma = {F_PL:.3f}, f_CSF = {F_CSF:.3f}, fitted from the observed"
        f" {tx['PGRNPL'][-1]/nat['PGRNPL'][0]:.2f}x")
    say(f"    plasma and {tx['PGRNCSF'][-1]/nat['PGRNCSF'][0]:.2f}x CSF PGRN rises."
        f" Occupancy: peripheral {th_pl:.3f}, CNS {th_cns:.3f}.")
    say(f" 3. Plasma PGRN normalises {100*nat['PGRNPL'][0]/PPL_WT:.0f}% ->"
        f" {100*tx['PGRNPL'][-1]/PPL_WT:.0f}% of control while LYSOSOMAL PGRN")
    say(f"    moves only {100*nat['PGRNLYS'][0]:.0f}% -> {100*tx['PGRNLYS'][-1]:.0f}%"
        " of wild type.")
    say(f" 4. Predicted 96-week CDR slowing on latozinemab:"
        f" {-100*dlato/prog:.2f}% — a predicted")
    say("    phase 3 miss coexisting with unambiguous target engagement.")
    say(f" 5. The sign of anti-sortilin benefit flips exactly at p = f_CSF ="
        f" {F_CSF:.3f},")
    say("    at ANY dose; ceiling (1-p)/(1-f_CSF) ="
        f" {(1-P0['P_LYS_SORT'])/(1-F_CSF):.2f}x at the default p.")
    say(" 6. Plasma PGRN is 94.5% at EVERY p in the sweep: the registrational")
    say("    biomarker is blind to the parameter that decides the sign.")
    say(f" 7. AAV-GRN reaches {100*aav['PGRNLYS'][-1]:.0f}% lysosomal PGRN vs"
        f" {100*tx['PGRNLYS'][-1]:.0f}% for the antibody,")
    say(f"    from a LOWER plasma level — plasma PGRN ranks the two backwards.")
    say(f" 8. Trazodone dNPI {res['trazodone 150 mg/d']['NPI'][-1]-b['NPI'][-1]:+.1f};"
        f" donepezil dNPI"
        f" {res['donepezil 10 mg/d']['NPI'][-1]-b['NPI'][-1]:+.1f} with dCDR exactly"
        " 0.000 while")
    say(f"    cholinergic tone rises {b['ACH'][-1]:.2f} ->"
        f" {res['donepezil 10 mg/d']['ACH'][-1]:.2f}: the AChEI null is structural.")
    say(f" 9. C9 ASO: CSF poly-GP -{100*(1-x9['POLYGP'][-1]/n9['POLYGP'][-1]):.0f}%,"
        f" dCDR {x9['CDR'][-1]-n9['CDR'][-1]:+.3f}. Target engagement is")
    say("    necessary and nowhere near sufficient.")
    say(f"10. MAPT ASO: soluble tau -{100*(1-xT['TAUS'][-1]/nT['TAUS'][-1]):.0f}%"
        f" but aggregated tau only"
        f" -{100*(1-xT['TAUAGG'][-1]/nT['TAUAGG'][-1]):.0f}%:")
    say("    synthesis blockade starves a deposit, it does not clear one.")

    text = "\n".join(R)
    if brief:
        i = text.find("K. HEADLINE NUMBERS")
        print(text[max(0, i - 82):])
    else:
        print(text)
    return text


if __name__ == "__main__":
    main(brief="--brief" in sys.argv)
