#!/usr/bin/env python3
# ============================================================================
#  MOGAD QSP model — dependency-free reference implementation
#  ---------------------------------------------------------------------------
#  This file is the executable twin of `mogad_mrgsolve_model.R`. It encodes the
#  SAME 34 ODEs with the SAME parameter values, integrated with fixed-step RK4
#  using nothing but the Python standard library.
#
#  Why it exists: every number quoted in README.md and in the calibration notes
#  of the R model is produced by running this file, so the claims are checkable
#  without an R/mrgsolve installation:
#
#      python3 mogad_reference_impl.py
#
#  It is a verification artefact, not the primary deliverable — the mrgsolve
#  model is. If the two ever disagree, the R file is authoritative for
#  structure and this file is authoritative for the quoted numbers.
# ============================================================================
from __future__ import annotations
import math

# ---------------------------------------------------------------------------
# Parameters (identical names/values to $PARAM in mogad_mrgsolve_model.R)
# ---------------------------------------------------------------------------
P = dict(
    # --- patient / phenotype ------------------------------------------------
    SEV=1.0,          # attack-drive scaling (0.5 mild .. 1.5 severe)
    SITE_ON=1.00,     # topography weight, optic nerve
    SITE_SC=0.45,     # topography weight, spinal cord
    AQP4_MODE=0.0,    # 1 = AQP4-IgG NMOSD comparator arm (astrocytopathy on)
    AQP4_OLD=25.0,    # extra oligodendrocyte-death coefficient when AQP4_MODE=1
    AQP4_GFAP=30.0,   # extra astrocytic GFAP release when AQP4_MODE=1
    BYST_GFAP=0.6,    # bystander GFAP release in MOGAD (small)

    # --- glucocorticoids ----------------------------------------------------
    KEL_MP=6.6,       # 1/d, methylprednisolone (t1/2 ~2.5 h)
    MPPOT=1.20,       # MP -> prednisone-equivalent potency factor
    KA_PRED=8.0,      # 1/d oral absorption
    F_PRED=0.85,      # oral bioavailability
    KEL_PRED=5.5,     # 1/d (t1/2 ~3 h)
    PRED0=0.0,        # mg/d, starting oral prednisone dose
    TAPER_DAYS=1.0,   # d, linear taper length of the oral prednisone wean
    PRED_START=0.0,   # d, day oral prednisone begins
    KIN_GC=0.70, KOUT_GC=0.70, GC50=4.70,
    IMAX_GC_PB=0.85,  # max plasmablast suppression by GC (reversible)
    IMAX_GC_IL6=0.60, IMAX_GC_T=0.55, IMAX_GC_BBB=0.60, IMAX_GC_M=0.35,

    # --- IVIG ---------------------------------------------------------------
    KEL_IVIG=0.040, EC50_IVIG=4.0,
    FRN_IVIG=0.45,    # max fractional acceleration of IgG catabolism
    KNEUT=1.35,       # effector-level neutralisation potency (FcgR block etc.)

    # --- rituximab ----------------------------------------------------------
    KEL_RTX=0.075, EC50_RTX=0.25, KDEP_RTX=0.60,
    # --- IL-6R blocker (tocilizumab / satralizumab) -------------------------
    KEL_TCZ=0.035, EC50_TCZ=0.15, IMAX_TCZ=0.85,
    # --- antimetabolite (MMF / azathioprine) --------------------------------
    KEL_MMF=0.35, EC50_MMF=0.60, IMAX_MMF=0.35, MMF_RATE=0.0,
    # --- FcRn inhibitor -----------------------------------------------------
    KEL_FCRN=0.30, EC50_FCRN=0.40, FRN_INH=2.30,
    # --- C5 inhibitor -------------------------------------------------------
    KEL_C5I=0.06, EC50_C5I=0.35, IMAX_C5I=0.95, IMAX_C5H=0.35,
    # --- plasma exchange ----------------------------------------------------
    KOUT_PLEX=4.0, KPLEX=0.90,
    # --- attack trigger -----------------------------------------------------
    KOUT_TRIG=0.30, ATRIG=1.5, BTRIG=6.0, CTRIG=25.0, DTRIG=5.0, ETRIG=70.0,

    # --- B lineage / antibody ----------------------------------------------
    KIN_MBC=0.008, KOUT_MBC=0.008,
    KGEN_PB=0.030435, KOUT_PB=0.140,
    FRAC_ESC=0.60,    # CD20-negative (rituximab-proof) share of PB generation
    KSYN_AB=0.099, KEL_AB=0.033,
    KSYN_IgG=0.330, FLOOR_IgG=0.40,
    EMAX_IL6=6.0, EC50_IL6=2.0,

    # --- cytokines / T cells / barrier -------------------------------------
    KIN_IL6=8.0, KOUT_IL6=4.0, KT_IL6=0.50,
    KIN_T=0.080, KOUT_T=0.080,
    KIN_BBB=0.00918, KOUT_BBB=0.50, KIL6_BBB=2.0, EC50_IL6B=6.0, KT_BBB=1.20,

    # --- CNS compartment & effectors ---------------------------------------
    KIN_CNS=0.050, KOUT_CNS=0.250,
    KIN_COMP=1.0, KOUT_COMP=1.0,
    KIN_MPHI=0.30, KOUT_MPHI=0.30, KC5A=2.0,
    KCDC=1.50, KADCP=8.0, KDIR=0.50,
    INJMAX=0.90, OPS50=1.20, HI=2.50,

    # --- tissue -------------------------------------------------------------
    KREMY=0.035,      # remyelination rate at full oligodendrocyte survival
    KOLD=0.08,        # oligodendrocyte death coefficient (SMALL in MOGAD)
    KOLREP=0.010,
    KAX=0.22, AXTHR=0.06,      # acute axonal loss above an injury threshold
    KAX2=0.010, DEMYTHR=0.35,  # chronic denuded-axon degeneration
    KIN_ED=3.0, KOUT_ED=0.15,

    # --- biomarkers ---------------------------------------------------------
    KIN_NFL=900.0, KBASE_NFL=0.64, KOUT_NFL=0.080,
    KIN_GFAP=110.0, KBASE_GFAP=6.40, KOUT_GFAP=0.080,
    KIN_CSF=25.0, KOUT_CSF=0.12, KM_CSF=1.0,

    # --- relapse hazard -----------------------------------------------------
    LAM0=0.00213,     # 1/d, maximal hazard (ARR_max = 0.78/y)
    AB50=1.80, HAZ_H=3.0,
    KAMP=6.0, EC50_AMP=1.2, POW_BBB=0.70,
    IL6_BASE=3.0, BBB_BASE=0.05,

    # --- toxicity -----------------------------------------------------------
    KBMD=2.70e-4, KBMDREC=1.0e-4,

    # --- endpoint mapping ---------------------------------------------------
    WM=2.00, WA=3.50, AXRES=0.15, RNFL0=98.0, TITER_REF=100.0,
)

# state order
S = ["MP_C", "PRED_G", "PRED_C", "GC_EFF", "IVIG_C", "RTX_C", "TCZ_C", "MMF_C",
     "FCRN_C", "C5I_C", "PLEX_D", "TRIG", "MBC", "PB", "IL6", "TCELL", "AB",
     "TOTIgG", "BBB", "AB_CNS", "COMP", "MPHI", "MYEL_ON", "AXON_ON",
     "MYEL_SC", "AXON_SC", "OL", "EDEMA", "NFL", "GFAP", "CSFC",
     "CUMHAZ", "CUMSTER", "BMD"]
IX = {n: i for i, n in enumerate(S)}


def init_state(p):
    y = [0.0] * len(S)
    y[IX["MBC"]] = 1.0
    y[IX["PB"]] = 1.0
    y[IX["IL6"]] = p["IL6_BASE"]
    y[IX["TCELL"]] = 1.0
    y[IX["AB"]] = 3.0
    y[IX["TOTIgG"]] = 10.0
    y[IX["BBB"]] = p["BBB_BASE"]
    y[IX["AB_CNS"]] = 0.030
    y[IX["COMP"]] = 0.030
    y[IX["MPHI"]] = 0.090
    y[IX["MYEL_ON"]] = 0.97
    y[IX["AXON_ON"]] = 1.0
    y[IX["MYEL_SC"]] = 0.97
    y[IX["AXON_SC"]] = 1.0
    y[IX["OL"]] = 1.0
    y[IX["NFL"]] = 8.0
    y[IX["GFAP"]] = 80.0
    y[IX["CSFC"]] = 3.0
    y[IX["BMD"]] = 1.0
    return y


def hillp(x, k, h):
    x = max(x, 0.0)
    xn = x ** h
    return xn / (k ** h + xn)


def derivs(t, y, p):
    (MP_C, PRED_G, PRED_C, GC_EFF, IVIG_C, RTX_C, TCZ_C, MMF_C, FCRN_C, C5I_C,
     PLEX_D, TRIG, MBC, PB, IL6, TCELL, AB, TOTIgG, BBB, AB_CNS, COMP, MPHI,
     MYEL_ON, AXON_ON, MYEL_SC, AXON_SC, OL, EDEMA, NFL, GFAP, CSFC,
     CUMHAZ, CUMSTER, BMD) = y
    d = [0.0] * len(S)

    # ---- glucocorticoid PK/PD ---------------------------------------------
    trel = t - p["PRED_START"]
    if p["PRED0"] > 0 and 0.0 <= trel < p["TAPER_DAYS"]:
        PRED_IN = p["PRED0"] * (1.0 - trel / p["TAPER_DAYS"])
    else:
        PRED_IN = 0.0
    d[IX["MP_C"]] = -p["KEL_MP"] * MP_C
    d[IX["PRED_G"]] = PRED_IN - p["KA_PRED"] * PRED_G
    d[IX["PRED_C"]] = p["KA_PRED"] * PRED_G * p["F_PRED"] - p["KEL_PRED"] * PRED_C
    PEQ = PRED_C + p["MPPOT"] * MP_C
    d[IX["GC_EFF"]] = p["KIN_GC"] * PEQ - p["KOUT_GC"] * GC_EFF
    GC_occ = GC_EFF / (p["GC50"] + GC_EFF)

    # ---- other drug PK ----------------------------------------------------
    d[IX["IVIG_C"]] = -p["KEL_IVIG"] * IVIG_C
    d[IX["RTX_C"]] = -p["KEL_RTX"] * RTX_C
    d[IX["TCZ_C"]] = -p["KEL_TCZ"] * TCZ_C
    d[IX["MMF_C"]] = p["MMF_RATE"] - p["KEL_MMF"] * MMF_C
    d[IX["FCRN_C"]] = -p["KEL_FCRN"] * FCRN_C
    d[IX["C5I_C"]] = -p["KEL_C5I"] * C5I_C
    d[IX["PLEX_D"]] = -p["KOUT_PLEX"] * PLEX_D
    d[IX["TRIG"]] = -p["KOUT_TRIG"] * TRIG

    RTX_occ = RTX_C / (p["EC50_RTX"] + RTX_C)
    IL6R_block = p["IMAX_TCZ"] * TCZ_C / (p["EC50_TCZ"] + TCZ_C)
    MMF_eff = p["IMAX_MMF"] * MMF_C / (p["EC50_MMF"] + MMF_C)
    FCRN_occ = FCRN_C / (p["EC50_FCRN"] + FCRN_C)
    C5I_block = p["IMAX_C5I"] * C5I_C / (p["EC50_C5I"] + C5I_C)
    IVIG_sat = IVIG_C / (p["EC50_IVIG"] + IVIG_C)

    # ---- B lineage / antibody --------------------------------------------
    IL6_sig = IL6 * (1.0 - IL6R_block)
    FIL6 = 1.0 + p["EMAX_IL6"] * IL6_sig / (p["EC50_IL6"] + IL6_sig)
    GCsup = p["IMAX_GC_PB"] * GC_occ

    d[IX["MBC"]] = (p["KIN_MBC"] * (1.0 - MMF_eff) * (1.0 + p["ATRIG"] * TRIG)
                    - p["KOUT_MBC"] * MBC - p["KDEP_RTX"] * RTX_occ * MBC)
    gen = ((1.0 - p["FRAC_ESC"]) * MBC + p["FRAC_ESC"])
    d[IX["PB"]] = (p["KGEN_PB"] * gen * FIL6 * (1.0 - GCsup)
                   * (1.0 - 0.5 * MMF_eff) * (1.0 + p["BTRIG"] * TRIG)
                   - p["KOUT_PB"] * PB)

    FcRn_acc = p["FRN_IVIG"] * IVIG_sat + p["FRN_INH"] * FCRN_occ
    d[IX["AB"]] = (p["KSYN_AB"] * PB - p["KEL_AB"] * AB * (1.0 + FcRn_acc)
                   - p["KPLEX"] * PLEX_D * AB)
    d[IX["TOTIgG"]] = (p["KSYN_IgG"] * (p["FLOOR_IgG"] + (1 - p["FLOOR_IgG"]) * MBC)
                       - p["KEL_AB"] * TOTIgG * (1.0 + FcRn_acc)
                       - p["KPLEX"] * PLEX_D * TOTIgG)

    d[IX["IL6"]] = (p["KIN_IL6"] * (1.0 + p["CTRIG"] * TRIG + p["KT_IL6"] * TCELL)
                    * (1.0 - p["IMAX_GC_IL6"] * GC_occ) - p["KOUT_IL6"] * IL6)
    d[IX["TCELL"]] = (p["KIN_T"] * (1.0 + p["DTRIG"] * TRIG)
                      * (1.0 - p["IMAX_GC_T"] * GC_occ) * (1.0 - MMF_eff)
                      - p["KOUT_T"] * TCELL)

    bbb_drive = (1.0 + p["ETRIG"] * TRIG
                 + p["KIL6_BBB"] * IL6_sig / (p["EC50_IL6B"] + IL6_sig)
                 + p["KT_BBB"] * TCELL)
    d[IX["BBB"]] = (p["KIN_BBB"] * bbb_drive * (1.0 - BBB)
                    - p["KOUT_BBB"] * BBB * (1.0 + p["IMAX_GC_BBB"] * GC_occ))

    # ---- CNS transit & effectors -----------------------------------------
    AB_eff = AB / (1.0 + p["KNEUT"] * IVIG_sat)
    d[IX["AB_CNS"]] = p["KIN_CNS"] * BBB * AB_eff - p["KOUT_CNS"] * AB_CNS
    OPS = AB_CNS
    d[IX["COMP"]] = p["KIN_COMP"] * OPS * (1.0 - C5I_block) - p["KOUT_COMP"] * COMP
    d[IX["MPHI"]] = (p["KIN_MPHI"] * (OPS + p["KC5A"] * COMP)
                     * (1.0 - p["IMAX_GC_M"] * GC_occ) - p["KOUT_MPHI"] * MPHI)

    OPS_eff = p["KCDC"] * COMP + p["KADCP"] * MPHI * OPS + p["KDIR"] * OPS
    INJ = p["SEV"] * p["INJMAX"] * hillp(OPS_eff, p["OPS50"], p["HI"])
    INJ_ON = p["SITE_ON"] * INJ
    INJ_SC = p["SITE_SC"] * INJ

    # ---- tissue -----------------------------------------------------------
    d[IX["MYEL_ON"]] = -INJ_ON * MYEL_ON + p["KREMY"] * OL * (1.0 - MYEL_ON)
    d[IX["MYEL_SC"]] = -INJ_SC * MYEL_SC + p["KREMY"] * OL * (1.0 - MYEL_SC)

    ax_on = (p["KAX"] * max(INJ_ON - p["AXTHR"], 0.0) * AXON_ON
             + p["KAX2"] * max(1.0 - MYEL_ON - p["DEMYTHR"], 0.0) * AXON_ON)
    ax_sc = (p["KAX"] * max(INJ_SC - p["AXTHR"], 0.0) * AXON_SC
             + p["KAX2"] * max(1.0 - MYEL_SC - p["DEMYTHR"], 0.0) * AXON_SC)
    d[IX["AXON_ON"]] = -ax_on
    d[IX["AXON_SC"]] = -ax_sc

    d[IX["OL"]] = (-p["KOLD"] * (1.0 + p["AQP4_MODE"] * p["AQP4_OLD"]) * INJ * OL
                   + p["KOLREP"] * (1.0 - OL))
    d[IX["EDEMA"]] = p["KIN_ED"] * INJ - p["KOUT_ED"] * EDEMA * (1.0 + 2.0 * GC_occ)

    # ---- biomarkers -------------------------------------------------------
    d[IX["NFL"]] = p["KIN_NFL"] * (ax_on + ax_sc) + p["KBASE_NFL"] - p["KOUT_NFL"] * NFL
    ASTRO_INJ = (p["BYST_GFAP"] + p["AQP4_MODE"] * p["AQP4_GFAP"]) * INJ
    d[IX["GFAP"]] = p["KIN_GFAP"] * ASTRO_INJ + p["KBASE_GFAP"] - p["KOUT_GFAP"] * GFAP
    d[IX["CSFC"]] = p["KIN_CSF"] * (BBB + p["KM_CSF"] * MPHI) - p["KOUT_CSF"] * CSFC

    # ---- relapse hazard ---------------------------------------------------
    Hill_AB = hillp(AB_eff, p["AB50"], p["HAZ_H"])
    amp_il6 = ((1.0 + p["KAMP"] * IL6_sig / (p["EC50_AMP"] + IL6_sig))
               / (1.0 + p["KAMP"] * p["IL6_BASE"] / (p["EC50_AMP"] + p["IL6_BASE"])))
    amp_bbb = (max(BBB, 1e-9) / p["BBB_BASE"]) ** p["POW_BBB"]
    AMP = amp_il6 * amp_bbb * (1.0 - p["IMAX_C5H"] * C5I_block)
    HAZ = p["LAM0"] * Hill_AB * AMP
    d[IX["CUMHAZ"]] = HAZ

    # ---- toxicity ---------------------------------------------------------
    d[IX["CUMSTER"]] = (p["KEL_PRED"] * PRED_C + p["KEL_MP"] * MP_C) / 1000.0
    d[IX["BMD"]] = -p["KBMD"] * GC_occ + p["KBMDREC"] * (1.0 - BMD)
    return d


def outputs(t, y, p):
    """Derived (algebraic) outputs — mirrors $TABLE in the R model."""
    g = dict(zip(S, y))
    GC_occ = g["GC_EFF"] / (p["GC50"] + g["GC_EFF"])
    IL6R_block = p["IMAX_TCZ"] * g["TCZ_C"] / (p["EC50_TCZ"] + g["TCZ_C"])
    IVIG_sat = g["IVIG_C"] / (p["EC50_IVIG"] + g["IVIG_C"])
    C5I_block = p["IMAX_C5I"] * g["C5I_C"] / (p["EC50_C5I"] + g["C5I_C"])
    AB_eff = g["AB"] / (1.0 + p["KNEUT"] * IVIG_sat)
    IL6_sig = g["IL6"] * (1.0 - IL6R_block)
    Hill_AB = hillp(AB_eff, p["AB50"], p["HAZ_H"])
    amp_il6 = ((1.0 + p["KAMP"] * IL6_sig / (p["EC50_AMP"] + IL6_sig))
               / (1.0 + p["KAMP"] * p["IL6_BASE"] / (p["EC50_AMP"] + p["IL6_BASE"])))
    amp_bbb = (max(g["BBB"], 1e-9) / p["BBB_BASE"]) ** p["POW_BBB"]
    AMP = amp_il6 * amp_bbb * (1.0 - p["IMAX_C5H"] * C5I_block)
    HAZ = p["LAM0"] * Hill_AB * AMP
    # Deficit = reversible conduction block (demyelination) + irreversible axonal
    # loss ABOVE a functional reserve. The reserve term is what lets MOGAD show
    # marked RNFL thinning together with good visual recovery (commitment 2).
    logMAR = min(3.0, p["WM"] * (1 - g["MYEL_ON"])
                 + p["WA"] * max(0.0, (1 - g["AXON_ON"]) - p["AXRES"]))
    edss = min(9.5, 4.5 * (1 - g["MYEL_SC"])
               + 9.0 * max(0.0, (1 - g["AXON_SC"]) - 0.05)
               + 0.48 * min(logMAR, 2.0))
    return dict(
        TIME=t, TITER=p["TITER_REF"] * g["AB"], AB=g["AB"], AB_EFF=AB_eff,
        GC_OCC=GC_occ, ARR=HAZ * 365.0, HAZ=HAZ,
        PREL=1.0 - math.exp(-g["CUMHAZ"]), CUMHAZ=g["CUMHAZ"],
        LOGMAR=logMAR, SNELLEN_DEN=20.0 * (10 ** logMAR), EDSS=edss,
        MYEL_ON=g["MYEL_ON"], AXON_ON=g["AXON_ON"], MYEL_SC=g["MYEL_SC"],
        AXON_SC=g["AXON_SC"], OL=g["OL"], BBB=g["BBB"], NFL=g["NFL"],
        GFAP=g["GFAP"], CSFC=g["CSFC"], RNFL=p["RNFL0"] * g["AXON_ON"],
        IGG_TOT=g["TOTIgG"] + g["IVIG_C"], IL6_OBS=g["IL6"] * (1 + 3 * IL6R_block),
        PB=g["PB"], MBC=g["MBC"], EDEMA=g["EDEMA"], CUMSTER=g["CUMSTER"],
        BMD=g["BMD"], IVIG_C=g["IVIG_C"], INJ_PROXY=g["AB_CNS"],
    )


def simulate(par_over=None, doses=None, tend=365.0, dt=0.01, record_every=1.0):
    """RK4 with instantaneous bolus events. doses = [(time, state, amount), ...]"""
    p = dict(P)
    if par_over:
        p.update(par_over)
    y = init_state(p)
    doses = sorted(doses or [], key=lambda z: z[0])
    di = 0
    rec = []
    n = int(round(tend / dt))
    next_rec = 0.0
    for k in range(n + 1):
        t = k * dt
        while di < len(doses) and doses[di][0] <= t + 1e-9:
            y[IX[doses[di][1]]] += doses[di][2]
            di += 1
        if t >= next_rec - 1e-9:
            rec.append(outputs(t, y, p))
            next_rec += record_every
        if k == n:
            break
        k1 = derivs(t, y, p)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(len(S))]
        k2 = derivs(t + 0.5 * dt, y2, p)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(len(S))]
        k3 = derivs(t + 0.5 * dt, y3, p)
        y4 = [y[i] + dt * k3[i] for i in range(len(S))]
        k4 = derivs(t + dt, y4, p)
        y = [y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
             for i in range(len(S))]
        # numeric guards (fractions stay in [0,1])
        for nm in ("MYEL_ON", "MYEL_SC", "AXON_ON", "AXON_SC", "OL", "BBB"):
            y[IX[nm]] = min(max(y[IX[nm]], 0.0), 1.0)
    return rec


# ---------------------------------------------------------------------------
# Dosing helpers
# ---------------------------------------------------------------------------
def ivmp(day0=4.0, n=5):
    """IV methylprednisolone 1 g daily x n -> 1250 mg prednisone-equivalent."""
    return [(day0 + i, "MP_C", 1250.0) for i in range(n)]


def plex(day0=2.0, n=5, spacing=2.0):
    return [(day0 + i * spacing, "PLEX_D", 1.0) for i in range(n)]


def ivig_acute(day0=1.0):
    """2 g/kg divided over 2 days ~ +24 g/L total."""
    return [(day0, "IVIG_C", 12.0), (day0 + 1, "IVIG_C", 12.0)]


def ivig_maint(start=0.0, tend=365.0, interval=28.0, gkg=1.0):
    return [(t, "IVIG_C", 12.0 * gkg)
            for t in frange(start, tend, interval)]


def rtx(start=0.0, tend=365.0):
    """1 g x2 two weeks apart, then 1 g every 6 months."""
    ev = [(start, "RTX_C", 1.0), (start + 14, "RTX_C", 1.0)]
    t = start + 182.0
    while t <= tend:
        ev.append((t, "RTX_C", 1.0))
        t += 182.0
    return ev


def il6rb(start=0.0, tend=365.0, interval=28.0):
    return [(t, "TCZ_C", 1.0) for t in frange(start, tend, interval)]


def fcrn(start=0.0, tend=365.0, interval=7.0):
    return [(t, "FCRN_C", 1.0) for t in frange(start, tend, interval)]


def c5i(start=0.0, tend=365.0, interval=14.0):
    return [(t, "C5I_C", 1.0) for t in frange(start, tend, interval)]


def trigger(day=0.0, amt=1.0):
    return [(day, "TRIG", amt)]


def frange(a, b, step):
    out, t = [], a
    while t <= b + 1e-9:
        out.append(t)
        t += step
    return out


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------
def base_attack():
    """Index attack at day 0 (trigger drives BBB opening + PB burst)."""
    return trigger(0.0, 1.0)


SCENARIOS = {
    "S1 Untreated index attack": dict(par={}, dose=base_attack()),
    "S2 IVMP x5, SHORT taper (60 mg, 28 d; 840 mg total)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=28.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5)),
    "S3 IVMP x5, LONG taper (18.7 mg, 90 d; 840 mg total)":
        dict(par=dict(PRED0=18.67, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5)),
    "S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5)),
    "S5 IVMP + PLEX 5x (severe attack)":
        dict(par=dict(SEV=1.35, PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + plex(6, 5, 2)),
    "S6 IVMP + acute IVIG 2 g/kg":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + ivig_acute(5)),
    "S7 Maintenance IVIG 1 g/kg q4wk":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + ivig_maint(28, 365, 28, 1.0)),
    "S8 Maintenance IVIG 0.4 g/kg q8wk (under-dosed)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + ivig_maint(28, 365, 56, 0.4)),
    "S9 Rituximab (1 g x2, then q6mo)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + rtx(28, 365)),
    "S10 IL-6R blockade (tocilizumab/satralizumab q4wk)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + il6rb(28, 365, 28)),
    "S11 MMF/azathioprine (continuous)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0, MMF_RATE=0.42),
             dose=base_attack() + ivmp(4, 5)),
    "S12 FcRn inhibitor (weekly SC)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + fcrn(28, 365, 7)),
    "S13 C5 inhibitor (q2wk)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + c5i(28, 365, 14)),
    "S14 IVIG + IL-6R blockade (combination)":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + ivig_maint(28, 365, 28, 1.0)
             + il6rb(28, 365, 28)),
    "S15 AQP4-NMOSD comparator (same attack, same IVMP)":
        dict(par=dict(AQP4_MODE=1.0, PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5)),
    "S16 Intercurrent infection at d200 on maintenance IVIG":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + ivig_maint(28, 365, 28, 1.0)
             + trigger(200.0, 1.0)),
    "S17 Intercurrent infection at d200, no maintenance":
        dict(par=dict(PRED0=60.0, TAPER_DAYS=90.0, PRED_START=9.0),
             dose=base_attack() + ivmp(4, 5) + trigger(200.0, 1.0)),
}


def summarise(rec):
    """Key readouts. 'Steady' ARR is taken at day 330 (post-taper, on
    maintenance) so it reflects the maintained state, not the attack."""
    by_t = {round(r["TIME"], 3): r for r in rec}

    def at(t):
        return by_t[min(by_t, key=lambda k: abs(k - t))]

    nadir = min(rec, key=lambda r: r["MYEL_ON"])
    peaklog = max(r["LOGMAR"] for r in rec)
    return dict(
        logMAR_nadir=peaklog,
        logMAR_d180=at(180)["LOGMAR"],
        logMAR_d365=at(365)["LOGMAR"],
        VA_den_d365=at(365)["SNELLEN_DEN"],
        MYEL_nadir=nadir["MYEL_ON"],
        AXON_ON_d365=at(365)["AXON_ON"],
        OL_min=min(r["OL"] for r in rec),
        EDSS_nadir=max(r["EDSS"] for r in rec),
        EDSS_d365=at(365)["EDSS"],
        NFL_peak=max(r["NFL"] for r in rec),
        GFAP_peak=max(r["GFAP"] for r in rec),
        CSF_peak=max(r["CSFC"] for r in rec),
        titer_d180=at(180)["TITER"],
        titer_d330=at(330)["TITER"],
        ARR_d330=at(330)["ARR"],
        ARR_mean_90_365=sum(r["HAZ"] for r in rec if 90 <= r["TIME"] <= 365)
        / max(1, len([r for r in rec if 90 <= r["TIME"] <= 365])) * 365,
        Prelapse_365=at(365)["PREL"],
        Prelapse_180=at(180)["PREL"],
        steroid_g=at(365)["CUMSTER"],
        BMD_d365=at(365)["BMD"],
        IgG_d330=at(330)["IGG_TOT"],
    )


def main():
    print("=" * 108)
    print("MOGAD QSP reference implementation — scenario summary (365-day horizon)")
    print("=" * 108)
    hdr = (f"{'scenario':<52}{'ARR@330':>8}{'P(rel)1y':>9}{'logMAR':>8}"
           f"{'logMAR':>8}{'EDSS':>6}{'titer':>7}{'NfL':>7}{'GFAP':>7}{'ster':>6}")
    print(hdr)
    print(f"{'':<52}{'':>8}{'':>9}{'nadir':>8}{'d365':>8}{'d365':>6}"
          f"{'d330':>7}{'peak':>7}{'peak':>7}{'g':>6}")
    print("-" * 108)
    results = {}
    for name, cfg in SCENARIOS.items():
        rec = simulate(cfg["par"], cfg["dose"], tend=365.0)
        s = summarise(rec)
        results[name] = s
        print(f"{name:<52}{s['ARR_d330']:>8.2f}{s['Prelapse_365']:>9.2f}"
              f"{s['logMAR_nadir']:>8.2f}{s['logMAR_d365']:>8.2f}"
              f"{s['EDSS_d365']:>6.1f}{s['titer_d330']:>7.0f}"
              f"{s['NFL_peak']:>7.0f}{s['GFAP_peak']:>7.0f}{s['steroid_g']:>6.2f}")
    print("-" * 108)

    # --- commitment 2: oligodendrocyte survival sets the recovery ceiling ---
    print("\nCommitment 2 — recovery ceiling is set by oligodendrocyte survival")
    for nm in ("S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)",
               "S15 AQP4-NMOSD comparator (same attack, same IVMP)"):
        s = results[nm]
        print(f"  {nm[:56]:<58} OL_min={s['OL_min']:.2f}  "
              f"logMAR nadir={s['logMAR_nadir']:.2f} -> d365={s['logMAR_d365']:.2f}  "
              f"final VA 20/{s['VA_den_d365']:.0f}  GFAP peak={s['GFAP_peak']:.0f}")

    # --- commitment 3: taper duration ---------------------------------------
    print("\nNEGATIVE RESULT — matched-dose taper duration barely matters "
          "(840 mg prednisone in both arms)")
    for nm in ("S2 IVMP x5, SHORT taper (60 mg, 28 d; 840 mg total)",
               "S3 IVMP x5, LONG taper (18.7 mg, 90 d; 840 mg total)"):
        s = results[nm]
        print(f"  {nm[:56]:<58} P(relapse) 6 mo={s['Prelapse_180']:.3f}  "
              f"1 y={s['Prelapse_365']:.3f}  cum steroid={s['steroid_g']:.2f} g")

    # --- taper-duration sweep ----------------------------------------------
    print("\nTaper-duration sweep at matched total dose (840 mg prednisone-equivalent)")
    print(f"  {'taper (d)':>10}{'start mg/d':>12}{'P(rel) 6mo':>12}{'P(rel) 1y':>11}"
          f"{'ARR@180':>9}{'BMD d365':>10}")
    for td in (14, 28, 45, 60, 90, 120, 180):
        p0 = 2 * 840.0 / td
        rec = simulate(dict(PRED0=p0, TAPER_DAYS=float(td), PRED_START=9.0),
                       base_attack() + ivmp(4, 5), tend=365.0)
        s = summarise(rec)
        arr180 = [r["ARR"] for r in rec if abs(r["TIME"] - 180) < 0.6][0]
        print(f"  {td:>10}{p0:>12.1f}{s['Prelapse_180']:>12.3f}"
              f"{s['Prelapse_365']:>11.3f}{arr180:>9.2f}{s['BMD_d365']:>10.3f}")

    # --- commitment 1/3: where does each drug act? -------------------------
    print("\nCommitments 1 & 3 — TITRE REDUCTION DOES NOT PREDICT RELAPSE "
          "PROTECTION; the node of action does")
    print(f"  {'therapy':<34}{'titre @330 (1:x)':>18}{'titre vs none':>15}"
          f"{'ARR@330':>9}{'ARR vs none':>13}")
    ref = results["S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)"]
    rows = [("none (steroid taper only)", "S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)"),
            ("maintenance IVIG 1 g/kg", "S7 Maintenance IVIG 1 g/kg q4wk"),
            ("rituximab", "S9 Rituximab (1 g x2, then q6mo)"),
            ("IL-6R blockade", "S10 IL-6R blockade (tocilizumab/satralizumab q4wk)"),
            ("MMF / azathioprine", "S11 MMF/azathioprine (continuous)"),
            ("FcRn inhibitor", "S12 FcRn inhibitor (weekly SC)"),
            ("C5 inhibitor", "S13 C5 inhibitor (q2wk)"),
            ("IVIG + IL-6R blockade", "S14 IVIG + IL-6R blockade (combination)")]
    for lbl, key in rows:
        s = results[key]
        print(f"  {lbl:<34}{s['titer_d330']:>18.0f}"
              f"{s['titer_d330']/ref['titer_d330']:>15.2f}"
              f"{s['ARR_d330']:>9.2f}{s['ARR_d330']/ref['ARR_d330']:>13.2f}")

    # --- observed-vs-model table ------------------------------------------
    print("\nCalibration check against published cohorts")
    obs = [("Relapsing MOGAD, no preventive therapy", 0.64, "0.58-0.70",
            results["S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)"]["ARR_d330"]),
           ("Maintenance IVIG >=1 g/kg q4wk", 0.22, "0.15-0.32",
            results["S7 Maintenance IVIG 1 g/kg q4wk"]["ARR_d330"]),
           ("IL-6R blockade", 0.09, "0.06-0.14",
            results["S10 IL-6R blockade (tocilizumab/satralizumab q4wk)"]["ARR_d330"]),
           ("Rituximab", 0.59, "0.59-0.63",
            results["S9 Rituximab (1 g x2, then q6mo)"]["ARR_d330"]),
           ("Mycophenolate", 0.67, "0.67-0.84",
            results["S11 MMF/azathioprine (continuous)"]["ARR_d330"])]
    print(f"  {'cohort':<40}{'observed ARR':>14}{'(range)':>12}{'model ARR':>11}")
    for lbl, o, rng, m in obs:
        print(f"  {lbl:<40}{o:>14.2f}{rng:>12}{m:>11.2f}")

    print("\nOther readouts quoted in README.md")
    for lbl, key in [("steroid taper only", "S4 IVMP x5, guideline taper (60 mg, 90 d; 2700 mg total)"),
                     ("rituximab", "S9 Rituximab (1 g x2, then q6mo)"),
                     ("FcRn inhibitor", "S12 FcRn inhibitor (weekly SC)"),
                     ("maintenance IVIG", "S7 Maintenance IVIG 1 g/kg q4wk"),
                     ("AQP4 comparator", "S15 AQP4-NMOSD comparator (same attack, same IVMP)"),
                     ("acute IVIG add-on", "S6 IVMP + acute IVIG 2 g/kg"),
                     ("PLEX, severe attack", "S5 IVMP + PLEX 5x (severe attack)")]:
        s2 = results[key]
        print(f"  {lbl:<24} CSF peak={s2['CSF_peak']:>5.0f}/uL  "
              f"total IgG d330={s2['IgG_d330']:>5.1f} g/L  "
              f"BMD d365={s2['BMD_d365']:.3f}  RNFL d365="
              f"{98*s2['AXON_ON_d365']:.0f} um  OL_min={s2['OL_min']:.2f}")
    print()


if __name__ == "__main__":
    main()
