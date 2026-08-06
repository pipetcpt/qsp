#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cgd_python_reference.py — the executable reference implementation of the
whole-body CGD QSP model.  The mrgsolve file (cgd_mrgsolve_model.R) and the
Shiny app mirror THIS file equation for equation; where they disagree, this
file is right, because this one has been run.

STRUCTURE
---------
Three layers, and the interface between them is deliberately narrow:

  LAYER 1  cgd_kernel.py — phagosomal chemistry.  Exports exactly one object,
           the survival of an ingested organism as a function of the oxidase
           activity phi of the neutrophil that ate it.  Nothing else crosses.

  LAYER 2  THE FOCUS (this file, focus_rhs) — one infection focus: extracellular
           organisms, organisms sheltering inside phagocytes that failed to kill
           them, recruited neutrophils, and antimicrobial drug effect.  Its
           output is the CRITICAL INOCULUM N_crit: the smallest number of
           organisms that grows instead of being cleared.

  LAYER 3  THE PATIENT (this file, patient_rhs) — 53 ODEs: granulopoiesis,
           the focus, Aspergillus, the cytokine network, granuloma and colitis,
           tryptophan/kynurenine, seven drug PK models, donor chimerism, gene
           therapy, and the survival hazard.

WHY N_crit IS THE CENTRE OF THE MODEL
--------------------------------------
A deterministic ODE cannot produce "0.83 infections per patient-year"; it
produces either clearance or sepsis.  The standard dodge is to bolt a fitted
hazard function onto the side of the model, which then does all the work and
learns nothing.  This model instead computes, from the focus equations, the
critical inoculum N_crit — the sharp deterministic boundary between clearance
and establishment — and then lets a FIXED, patient-independent exposure
distribution (how often you inhale/inoculate how many organisms) decide how
often that boundary is crossed:

    infection rate = lambda * P( log10 N0 > log10 N_crit )

lambda and the mean inoculum are calibrated ONCE, on two targets: the bacterial
infection rate of an untreated X-linked CGD null patient (1.90/patient-year,
Margolis 1990) and of a healthy control (~0.02/patient-year).  After that,
EVERY other rate in this model is a prediction: the effect of co-trimoxazole,
of itraconazole, of interferon gamma, of 12% versus 25% donor chimerism, of a
hypomorphic mutation, of being a female carrier.  None of them is fitted.

WHAT IS CALIBRATED IN THIS FILE (eight numbers, all listed here)
-----------------------------------------------------------------
    lambda_b, mu_N_b   bacterial exposure: rate and mean inoculum
    lambda_f, mu_N_f   fungal exposure: rate and mean inoculum
    k_h_bact, k_h_fung mortality per unit established infection burden
    amp_IL1            how much the missing oxidase amplifies IL-1beta
    k_gran             granuloma formation gain

Everything else is fixed from the literature or from Layer 1.
"""

import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import norm

import cgd_kernel as KER

HERE = os.path.dirname(os.path.abspath(__file__))

# ============================================================================
# PARAMETERS
# ============================================================================
P = dict(

    # ---- granulopoiesis and neutrophil kinetics ---------------------------
    # Dancey 1976 / Cartwright 1964: marrow mitotic pool ~2.1e9/kg, post-mitotic
    # reserve ~5.6e9/kg, circulating ANC 4000/uL with a blood half-life of ~7 h
    # and a marrow transit time of ~6 d.
    tau_prog   = 6.0,        # d, mitotic pool transit
    tau_res    = 5.0,        # d, post-mitotic reserve transit
    kout_ANC   = math.log(2) / (7.0 / 24.0),   # /d, blood disappearance
    ANC0       = 4000.0,     # /uL
    GCSF0      = 0.025,      # ng/mL basal
    kin_GCSF   = 1.0,        # /d
    kout_GCSF  = 24.0,       # /d
    GCSF_EC50  = 0.20,       # ng/mL, effect on marrow release
    GCSF_Emax  = 6.0,
    # G-CSF is cleared by the neutrophils themselves (receptor-mediated), which
    # is what closes the feedback loop.
    kcl_GCSF_N = 4.0e-3,     # /d per (cell/uL)

    MONO0      = 400.0,      # /uL
    kout_MONO  = 1.0,        # /d

    # ---- recruitment to the focus -----------------------------------------
    k_rec      = 40.0,       # /d, maximal neutrophil recruitment
    Krec       = 50.0,       # CFU, half-maximal recruitment
    kout_NEUT  = 1.2,        # /d, neutrophil death at the focus
    k_rec_mono = 4.0,
    kout_MAC   = 0.10,       # /d

    # ---- the focus: bacteria ----------------------------------------------
    mu_b       = 4.0,        # /d net growth of S. aureus in tissue
    Bmax       = 1.0e10,     # CFU carrying capacity of one focus
    k_up       = 3.0e-2,     # /d per (cell/uL): phagocytic uptake rate
    k_release  = 0.50,       # /d, organisms released when the phagocyte dies
    mu_bi      = 0.10,       # /d replication inside a failed phagosome
    Km_up      = 1.0e7,      # CFU, saturation of uptake
    N_resident = 200.0,      # cells/uL-equivalent of resident phagocytes

    # RETURN FACTOR.  Each phagocytic uptake returns  s_enc * k_release /
    # (k_release - mu_bi)  organisms to the extracellular pool.  At phi = 1 that
    # is 0.06 and phagocytosis is almost perfectly absorbing; at phi = 0 it is
    # 0.75, so a CGD neutrophil removes only a quarter of what it swallows and
    # hands the rest back alive.  If this factor is allowed above 1 the
    # phagocyte becomes a net AMPLIFIER, every CGD patient is septic from birth,
    # and the model is wrong; an early draft of this file had exactly that, with
    # mu_bi = 0.35/d, and produced a critical inoculum of one organism.

    # Necrosis.  Without it the focus outcome is independent of inoculum size,
    # every exposure has the same fate, and there is no critical inoculum to
    # find.  With it, a focus that grows past ~1e5 organisms before the
    # neutrophils arrive walls itself off from both the neutrophils and the
    # antibiotic — which is the clinical fact that a CGD abscess has to be
    # drained and cannot be talked out of existence with co-trimoxazole.
    k_nec      = 3.0,        # /d
    Knec       = 1.0e5,      # CFU
    kout_nec   = 0.15,       # /d

    # ---- the focus: Aspergillus -------------------------------------------
    # Conidia are handled by macrophages (oxidase-dependent); only neutrophils
    # damage hyphae, and hyphal damage is the most oxidant-dependent killing in
    # the whole model (omega = 0.92).  This asymmetry is why CGD is above all a
    # mould disease and why the fungal arm responds to chimerism differently
    # from the bacterial arm.
    k_germ     = 0.35,       # /d conidium -> hypha (slow in vivo)
    mu_h       = 2.2,        # /d hyphal extension
    Hmax       = 1.0e9,
    k_up_c     = 0.25,       # /d per (cell/uL) macrophage conidial uptake
    k_dam_h    = 3.0e-2,     # /d per (cell/uL) neutrophil hyphal damage
    k_GM       = 3.0e-8,     # galactomannan index per hypha per day
    kout_GM    = 1.4,        # /d

    # ---- cytokines (pg/mL unless stated) -----------------------------------
    IL1B0 = 2.0,   kout_IL1B = 24.0,
    IL180 = 60.0,  kout_IL18 = 12.0,
    TNF0  = 5.0,   kout_TNF  = 40.0,
    IL60  = 2.0,   kout_IL6  = 24.0,
    IL170 = 5.0,   kout_IL17 = 6.0,
    IL100 = 4.0,   kout_IL10 = 24.0,
    IFNg0 = 3.0,   kout_IFNg = 30.0,
    CRP0  = 1.0,   kout_CRP  = 0.37,     # /d, CRP t1/2 19 h
    k_CRP_IL6 = 0.185,

    k_IL1_stim = 0.030,      # per CFU-equivalent of stimulus
    k_IL1_apop = 0.60,
    k_TNF_stim = 0.020,
    k_IL6_TNF  = 0.35,
    k_IL6_IL1  = 0.55,
    k_IL17_IL1 = 0.020,
    k_IL10_eff = 0.40,       # IL-10 released by efferocytosis of apoptotic PMN
    k_IFNg_IL18 = 0.010,

    # ---- the efferocytosis defect ------------------------------------------
    # Phosphatidylserine externalisation on the apoptotic neutrophil is itself
    # an oxidation-dependent event, so the CGD neutrophil dies without being
    # recognised.  The uncleared corpse is both an IL-1beta stimulus and the
    # missing source of the IL-10 that would have switched the lesion off.
    # This is the single mechanism that makes CGD inflammation self-sustaining
    # in the absence of any organism.
    k_apop     = 0.55,       # /d apoptotic neutrophils generated at the focus
    N_tissue0  = 50.0,       # resident tissue phagocytes turning over everywhere
    gut_stim   = 1.0,        # constant translocated microbial stimulus (everyone)
    k_effero   = 3.0,        # /d clearance, at full oxidase competence
    effero_min = 0.25,       # residual clearance at phi = 0

    # ---- granuloma and colitis ---------------------------------------------
    kout_GRAN  = 0.020,      # /d
    IL1_gran_EC50 = 30.0,    # pg/mL
    k_fib      = 0.010,
    kout_FIB   = 0.0015,
    k_COL      = 0.320,
    kout_COL   = 0.10,

    # ---- tryptophan / kynurenine -------------------------------------------
    # Romani 2008 proposed that NADPH oxidase is required for IDO activity and
    # that the tryptophan pathway, not the oxidant, is the real effector.  The
    # claim did not replicate (De Ravin 2010).  It is carried here as an
    # explicitly switchable arm, default OFF, so that the cost of believing it
    # can be priced rather than argued about.
    TRP0 = 55.0, kin_TRP = 55.0, kout_TRP = 1.0,
    k_IDO = 0.35, IDO_IFNg_EC50 = 20.0, kout_KYN = 3.0,

    # ---- NETs ---------------------------------------------------------------
    k_NET = 0.30, kout_NET = 1.0,

    # ---- PHARMACOKINETICS ---------------------------------------------------
    # Reference patient: 30 kg, BSA 1.1 m2 (CGD is diagnosed in childhood).
    # CONVENTION, and it is worth stating because getting it wrong is the most
    # common way a QSP model quietly reports nonsense: every drug amount is in
    # mg and every volume in L, so every concentration is mg/L.  The only
    # exception is interferon gamma-1b, whose amounts are in ug and whose
    # volume is in L, so its concentration is ug/L = ng/mL.

    # co-trimoxazole, trimethoprim 5 mg/kg/d divided BID
    TMP_F = 0.90, TMP_ka = 12.0, TMP_Vc = 30.0, TMP_Vp = 18.0,
    TMP_Q = 30.0, TMP_CL = 82.0,        # L/d, gives t1/2 ~10 h
    SMX_F = 0.85, SMX_ka = 10.0, SMX_V = 6.3, SMX_CL = 10.5,
    TMP_fu = 0.55, SMX_fu = 0.30,

    # itraconazole 200 mg BID capsule + the active hydroxy metabolite
    ITZ_F = 0.55, ITZ_ka = 4.0, ITZ_Vc = 120.0, ITZ_Vp = 250.0,
    ITZ_Q = 150.0, ITZ_CL = 220.0,      # L/d
    OHITZ_fm = 0.75, OHITZ_V = 300.0, OHITZ_CL = 150.0,

    # voriconazole, saturable (CYP2C19) elimination
    VOR_F = 0.95, VOR_ka = 12.0, VOR_V = 138.0,
    VOR_Vmax = 1400.0, VOR_Km = 3.0,    # mg/d, mg/L

    # interferon gamma-1b 50 ug/m2 SC three times weekly  (ug and ng/mL)
    IFN_F = 0.90, IFN_ka = 6.0, IFN_V = 30.0, IFN_CL = 100.0,

    # prednisolone
    PRED_F = 0.85, PRED_ka = 20.0, PRED_V = 18.0, PRED_CL = 100.0, PRED_fu = 0.25,

    # anakinra 2 mg/kg/d SC
    ANA_F = 0.95, ANA_ka = 6.0, ANA_V = 18.0, ANA_CL = 60.0,

    # ---- PHARMACODYNAMICS ---------------------------------------------------
    TMP_Emax = 6.0,  TMP_EC50 = 0.50,    # /d, mg/L free TMP, on S. aureus
    ITZ_Emax = 2.4,  ITZ_EC50 = 0.25,    # /d, on hyphal extension
    VOR_Emax = 3.4,  VOR_EC50 = 0.50,

    # Interferon gamma-1b.  The 1991 trial found a 67% reduction in serious
    # infection WITHOUT a reproducible rise in superoxide production, so this
    # model is not allowed to let it touch phi.  It acts only on the
    # oxidase-INDEPENDENT arm and on recruitment.  The consequence is a
    # falsifiable prediction: its benefit should be roughly FLAT across residual
    # ROS, unlike itraconazole's.  Section G tests that.
    IFN_Emax_nonox = 1.9, IFN_EC50 = 0.30,   # ng/mL
    IFN_Emax_rec   = 0.7,

    PRED_Emax_cyt = 0.85, PRED_EC50 = 12.0,  # ng/mL FREE prednisolone
    PRED_Emax_mig = 0.45,
    ANA_Imax = 0.92, ANA_IC50 = 250.0,       # ng/mL

    # ---- transplant and gene therapy ---------------------------------------
    k_engraft   = 0.12,      # /d approach to the final chimerism plateau
    k_vcn_loss  = 0.0011,    # /d slow loss of gene-marked myelopoiesis
    gvhd_hazard = 0.0,       # set per scenario

    # ---- mortality ----------------------------------------------------------
    h_base = 0.0008,         # /yr baseline, converted to /d below
    kout_ALT = 0.35, ALT0 = 25.0,
)

# The eight numbers calibrated in this file.  Placeholders until calibrate().
CALW = dict(lambda_b=None, muN_b=None, lambda_f=None, muN_f=None,
            k_h_bact=None, k_h_fung=None, amp_IL1=None, k_gran=None)

SDN = 1.0        # log10 SD of the inoculum size distribution, FIXED not fitted

WHOLE_TARGETS = dict(
    rate_bact_cgd_null   = 1.90,   # /patient-year, no prophylaxis   Margolis 1990
    rate_bact_healthy    = 0.02,   # /patient-year
    rate_fung_cgd_null   = 0.19,   # /patient-year                   Winkelstein 2000
    rate_fung_healthy    = 0.0005,
    mortality_cgd_null   = 0.050,  # /patient-year, pre-prophylaxis era
    colitis_index_cgd    = 5.0,    # 0-10 activity index in established disease
    gran_index_cgd       = 1.0,    # normalised granuloma burden
)


# ============================================================================
# GENOTYPES — every one of them is just a way of writing (f_normal, phi_res)
# ============================================================================
GENOTYPES = {
    # name                       f_normal  phi_res   note
    "healthy":                   (1.00, 1.00),
    "X-CGD null (gp91phox-)":    (0.00, 0.00),
    "X-CGD hypomorph 5%":        (0.00, 0.05),
    "X-CGD hypomorph 20%":       (0.00, 0.20),
    "p47phox-deficient (AR)":    (0.00, 0.12),   # AR-CGD retains measurable ROS
    "p67phox-deficient (AR)":    (0.00, 0.03),
    "X-CGD carrier 50%":         (0.50, 0.00),
    "X-CGD carrier 20%":         (0.20, 0.00),
    "X-CGD carrier 5%":          (0.05, 0.00),
}


def dhr_index(f_normal, phi_res):
    """What the DHR-123 flow assay reports: one number for the whole
    population.  Two very different marrows can produce the same one."""
    return KER.dhr_mean(f_normal, phi_res)


# ============================================================================
# LAYER 2 — THE FOCUS, AND THE CRITICAL INOCULUM
# ============================================================================
def _sat(x, k):
    return x / (x + k)

# ----------------------------------------------------------------------------
# EXTINCTION.  An ODE has no notion of "one organism", so a state left at 1e-30
# by the integrator's round-off is still multiplied by mu_b = 4/d, and over a
# 180-day simulation exp(720) turns it into 1e10 CFU.  A first draft of this
# file did exactly that and reported that an untreated CGD patient who was never
# inoculated with anything develops a 1e10 CFU abscess.  Every replicating
# compartment therefore carries the factor X/(X+1): below one organism there is
# nothing left to divide.
# ----------------------------------------------------------------------------
def _alive(x):
    return x / (x + 1.0)




def focus_rhs(t, y, ctx):
    """One infection focus.  Four states:
        B    extracellular organisms
        Bi   organisms sheltering inside phagocytes that ingested them and
             failed to kill them
        N    recruited neutrophils
        NEC  necrotic / walled-off tissue

    Bi is the whole reason mosaicism and uniform residual activity are not the
    same disease.  An organism ingested by an oxidase-NULL neutrophil is not
    merely un-killed: it is hidden from every competent neutrophil in the
    lesion, it goes on dividing, and it is released alive when its host cell
    dies.  An organism ingested by a uniformly WEAK neutrophil is at least
    being killed slowly the whole time.  A DHR mean cannot tell these apart.
    These equations can.

    NEC is what makes the outcome depend on the SIZE of the inoculum rather
    than only on the patient.  A focus that outruns neutrophil recruitment
    walls itself off from the neutrophils and from the antibiotic alike.
    """
    B, Bi, N, NEC = np.maximum(y, 0.0)
    p = ctx["P"]

    pen = 1.0 / (1.0 + NEC)                     # delivery into the lesion
    NT = (N + p["N_resident"]) * pen
    uptake = p["k_up"] * NT * B / (1.0 + B / p["Km_up"])

    dB = (p["mu_b"] * B * _alive(B) * (1.0 - B / p["Bmax"]) - uptake
          + p["k_release"] * Bi - ctx["drug_kill_b"] * pen * B)
    dBi = uptake * ctx["s_enc"] + p["mu_bi"] * Bi * _alive(Bi) - p["k_release"] * Bi
    dN = (p["k_rec"] * ctx["rec_gain"] * _sat(B + Bi, p["Krec"])
          * ctx["ANC"] / p["ANC0"] * 20.0 - p["kout_NEUT"] * N)
    dNEC = p["k_nec"] * _sat(B + Bi, p["Knec"]) - p["kout_nec"] * NEC
    return [dB, dBi, dN, dNEC]


SERIOUS_B = 1.0e7      # CFU at one focus = a clinically serious bacterial infection
SERIOUS_H = 1.0e6      # hyphal units = invasive aspergillosis


def focus_outcome(N0, ctx, days=45.0):
    """Integrate one inoculum.

    The endpoint is deliberately NOT sterilisation.  A CGD patient carries
    smouldering, subclinical foci more or less permanently, and demanding that
    the model drive every one of them below a single organism would classify a
    harmless 40-CFU nodule and a lung abscess as the same event.  The endpoint
    that the trials actually counted, and that this model therefore counts, is
    a SERIOUS infection: a focus that reaches 1e7 organisms."""
    sol = solve_ivp(focus_rhs, (0.0, days), [float(N0), 0.0, 0.0, 0.0], args=(ctx,),
                    method="LSODA", rtol=1e-7, atol=1e-6,
                    t_eval=np.linspace(0, days, 181))
    tot = sol.y[0] + sol.y[1]
    peak = float(np.max(tot))
    return dict(cleared=bool(peak < SERIOUS_B), peak=peak,
                end=float(tot[-1]), nec=float(np.max(sol.y[3])))


def boosted_survivals(f_normal, phi_res, hp, nonox_boost=1.0):
    """Survival of an ingested organism in a competent and in a mutant cell,
    with interferon gamma-1b's effect on the oxidase-INDEPENDENT arm applied.

    K0 = K(phi = 0) is, by construction, the whole non-oxidative log-kill,
    because at phi = 0 no HOCl is produced at all.  Boosting that arm by a
    factor b therefore ADDS (b - 1) * K0 log to every cell, competent or not.
    It is additive in log-kill and identical for both subpopulations, which is
    exactly the mechanistic claim the 1991 trial forces on us: interferon did
    not restore superoxide, so whatever it does must be phi-independent.

    A first draft wrote this as s ** (1 / b), which — since s < 1 and 1/b < 1 —
    RAISED survival, and the model duly reported that interferon gamma more
    than doubles the infection rate in CGD.
    """
    extra = (max(float(nonox_boost), 0.0) - 1.0) * hp["K0"]
    s1 = 10.0 ** (-(KER.K_hill(1.0, hp) + extra))
    sr = 10.0 ** (-(KER.K_hill(phi_res, hp) + extra))
    return min(max(s1, 1e-12), 1.0), min(max(sr, 1e-12), 1.0)


def make_ctx(f_normal, phi_res, hp, ANC=None, drug_kill_b=0.0, rec_gain=1.0,
             nonox_boost=1.0, p=None):
    """Assemble everything the focus needs.  s_enc is where LAYER 1 enters."""
    p = p or P
    s1, sr = boosted_survivals(f_normal, phi_res, hp, nonox_boost)
    s_enc = f_normal * s1 + (1.0 - f_normal) * sr
    return dict(P=p, s_enc=float(min(max(s_enc, 1e-12), 1.0)),
                ANC=ANC if ANC is not None else p["ANC0"],
                drug_kill_b=drug_kill_b, rec_gain=rec_gain)


def critical_inoculum(ctx, lo=1.0, hi=1e9, days=45.0):
    """The sharp deterministic boundary: the smallest inoculum that establishes,
    returned as log10 CFU.

    +inf  no inoculum up to 1e9 produces a serious infection (protected)
    -inf  a single organism produces a serious infection (no defence at all)
    """
    if focus_outcome(hi, ctx, days)["cleared"]:
        return float("inf")
    if not focus_outcome(lo, ctx, days)["cleared"]:
        return -float("inf")
    a, b = math.log10(lo), math.log10(hi)
    for _ in range(40):
        m = 0.5 * (a + b)
        if focus_outcome(10.0 ** m, ctx, days)["cleared"]:
            a = m
        else:
            b = m
        if b - a < 1e-3:
            break
    return 0.5 * (a + b)


def infection_rate(logNcrit, lam, muN, sdn=SDN):
    """Exposure crosses the boundary this often."""
    if not np.isfinite(logNcrit):
        return 0.0 if logNcrit > 0 else lam
    return float(lam * norm.sf((logNcrit - muN) / sdn))


# ---- the fungal focus, same idea, different cells --------------------------
def fungal_focus_rhs(t, y, ctx):
    """The fungal focus is NOT the bacterial focus with different constants.
    A conidium is ingestible and is handled inside a macrophage phagosome, so
    it sees exactly the kernel's s(phi).  A germinated hypha is far too large
    to ingest and can only be attacked from the outside, by oxidants and NETs
    released onto its surface — the least salvageable form of killing there is.
    That asymmetry, not any fitted constant, is why Aspergillus is the leading
    cause of death in CGD and why the fungal arm answers to chimerism on a
    different scale from the bacterial arm."""
    C, H, N, M, NEC = np.maximum(y, 0.0)
    p = ctx["P"]
    pen = 1.0 / (1.0 + NEC)
    NT = (N + p["N_resident"]) * pen
    # A conidium has two exits and only one of them is protective: it is
    # ingested by a macrophage (and then survives with probability s_conid,
    # straight out of the kernel), or it escapes and germinates outright.
    ingest = p["k_up_c"] * M * pen
    germ = (ingest * ctx["s_conid"] + p["k_germ"]) * C
    dC = -(ingest + p["k_germ"]) * C
    dH = (germ + p["mu_h"] * H * _alive(H) * (1.0 - H / p["Hmax"])
          - p["k_dam_h"] * NT * H * ctx["hyphal_competence"]
          - ctx["drug_kill_f"] * pen * H)
    dN = (p["k_rec"] * ctx["rec_gain"] * _sat(H + C, p["Krec"])
          * ctx["ANC"] / p["ANC0"] * 20.0 - p["kout_NEUT"] * N)
    dM = (p["k_rec_mono"] * _sat(C + H, p["Krec"]) * ctx["MONO"] / p["MONO0"] * 20.0
          - p["kout_MAC"] * M)
    dNEC = p["k_nec"] * _sat(H, p["Knec"]) - p["kout_nec"] * NEC
    return [dC, dH, dN, dM, dNEC]


def fungal_outcome(N0, ctx, days=45.0):
    sol = solve_ivp(fungal_focus_rhs, (0.0, days), [float(N0), 0.0, 0.0, 20.0, 0.0],
                    args=(ctx,), method="LSODA", rtol=1e-7, atol=1e-6,
                    t_eval=np.linspace(0, days, 181))
    peak = float(np.max(sol.y[1]))
    return dict(cleared=bool(peak < SERIOUS_H), peak=peak,
                conid_end=float(sol.y[0][-1]))


def make_fctx(f_normal, phi_res, hp, ANC=None, MONO=None, drug_kill_f=0.0,
              rec_gain=1.0, nonox_boost=1.0, p=None):
    p = p or P
    org = KER.ORGANISMS["Aspergillus fumigatus"]
    # Conidial survival inside the macrophage uses the same kernel, with
    # Aspergillus' own omega (0.92) — the most oxidant-dependent organism in the
    # table.  Hyphal competence is the neutrophil's ability to damage a target
    # far too large to ingest, which is essentially pure oxidant + NET work.
    s1, sr = boosted_survivals(f_normal, phi_res, hp, nonox_boost)
    s_conid = f_normal * s1 + (1.0 - f_normal) * sr
    hyph = f_normal * 1.0 + (1.0 - f_normal) * (0.10 + 0.90 * phi_res)
    return dict(P=p, s_conid=float(min(max(s_conid, 1e-12), 1.0)),
                hyphal_competence=float(hyph) * nonox_boost,
                ANC=ANC if ANC is not None else p["ANC0"],
                MONO=MONO if MONO is not None else p["MONO0"],
                drug_kill_f=drug_kill_f, rec_gain=rec_gain)


def critical_inoculum_f(fctx, lo=1.0, hi=1e9, days=45.0):
    if fungal_outcome(hi, fctx, days)["cleared"]:
        return float("inf")
    if not fungal_outcome(lo, fctx, days)["cleared"]:
        return -float("inf")
    a, b = math.log10(lo), math.log10(hi)
    for _ in range(40):
        m = 0.5 * (a + b)
        if fungal_outcome(10.0 ** m, fctx, days)["cleared"]:
            a = m
        else:
            b = m
        if b - a < 1e-3:
            break
    return 0.5 * (a + b)


# ============================================================================
# LAYER 3 — THE PATIENT: 53 ODEs
# ============================================================================
SV = ["PROG", "RES", "ANC", "GCSF", "MONO", "NEUT_T", "MAC", "EPI",
      "BACT", "BACT_i", "NEC", "CONID", "HYPH", "GM",
      "IL1B", "IL18", "TNF", "IL6", "IL17", "IL10", "IFNg", "CRP",
      "APOP", "NET", "GRAN", "FIB", "COL", "TRP", "KYN",
      "TMP_g", "TMP_c", "TMP_p", "SMX_g", "SMX_c",
      "ITZ_g", "ITZ_c", "ITZ_p", "OHITZ",
      "VOR_g", "VOR_c", "IFN_sc", "IFN_c",
      "PRED_g", "PRED_c", "ANA_sc", "ANA_c",
      "CHIM", "VCN", "CUMBACT", "CUMFUNG", "SURV", "STEROID", "ALT"]
S = {n: i for i, n in enumerate(SV)}
NS = len(SV)
assert NS == 53, NS


def patient_rhs(t, y, R):
    p = R["P"]
    c = R["CAL"]
    hp = R["hp"]
    v = np.maximum(y, 0.0)
    g = lambda n: v[S[n]]

    # ---------------- drug concentrations ----------------------------------
    # All amounts mg, all volumes L, so all concentrations mg/L.  Interferon
    # gamma-1b alone is in ug and L, so its concentration is ug/L = ng/mL.
    TMP  = g("TMP_c") / p["TMP_Vc"]                     # mg/L
    TMPf = TMP * p["TMP_fu"]
    TMPp = g("TMP_p") / p["TMP_Vp"]
    SMX  = g("SMX_c") / p["SMX_V"]
    ITZ  = g("ITZ_c") / p["ITZ_Vc"]
    ITZp = g("ITZ_p") / p["ITZ_Vp"]
    OHI  = g("OHITZ") / p["OHITZ_V"]
    VOR  = g("VOR_c") / p["VOR_V"]
    IFNd = g("IFN_c") / p["IFN_V"]                      # ng/mL
    PREDf = g("PRED_c") / p["PRED_V"] * 1000.0 * p["PRED_fu"]   # ng/mL free
    ANA  = g("ANA_c") / p["ANA_V"] * 1000.0             # ng/mL

    E_TMP  = p["TMP_Emax"] * TMPf / (p["TMP_EC50"] + TMPf)
    E_AZOL = (p["ITZ_Emax"] * (ITZ + 0.6 * OHI) / (p["ITZ_EC50"] + ITZ + 0.6 * OHI)
              + p["VOR_Emax"] * VOR / (p["VOR_EC50"] + VOR))
    E_IFNn = 1.0 + p["IFN_Emax_nonox"] * IFNd / (p["IFN_EC50"] + IFNd)
    E_IFNr = 1.0 + p["IFN_Emax_rec"] * IFNd / (p["IFN_EC50"] + IFNd)
    I_PRED = p["PRED_Emax_cyt"] * PREDf / (p["PRED_EC50"] + PREDf)
    I_PMIG = p["PRED_Emax_mig"] * PREDf / (p["PRED_EC50"] + PREDf)
    I_ANA  = p["ANA_Imax"] * ANA / (p["ANA_IC50"] + ANA)

    # ---------------- the oxidase state of the marrow -----------------------
    # Three sources of oxidase-competent neutrophils, all of them the same
    # variable as far as the phagosome is concerned: what you were born with,
    # what a donor gave you, and what a vector put back.
    f_norm = min(1.0, R["f_normal"] + g("CHIM") + g("VCN"))
    phi_res = R["phi_res"]
    s1, sr = boosted_survivals(f_norm, phi_res, hp, E_IFNn)
    s_enc = min(max(f_norm * s1 + (1.0 - f_norm) * sr, 1e-12), 1.0)
    # phi_eff is the population-average oxidase competence used by the arms that
    # genuinely respond to a mean (efferocytosis, inflammasome restraint), as
    # opposed to killing, which responds to the distribution.
    phi_eff = f_norm + (1.0 - f_norm) * phi_res

    # ---------------- granulopoiesis ----------------------------------------
    GCSF = g("GCSF")
    rel = 1.0 + p["GCSF_Emax"] * GCSF / (p["GCSF_EC50"] + GCSF)
    prod0 = p["ANC0"] * p["kout_ANC"]
    dPROG = prod0 - g("PROG") / p["tau_prog"]
    dRES  = g("PROG") / p["tau_prog"] - g("RES") / p["tau_res"] * rel / R["rel0"]
    dANC  = (g("RES") / p["tau_res"] * rel / R["rel0"] - p["kout_ANC"] * g("ANC")
             - p["k_rec"] * E_IFNr * (1.0 - I_PMIG)
             * _sat(g("BACT") + g("BACT_i") + g("HYPH"), p["Krec"]) * g("ANC") / p["ANC0"] * 20.0)
    dGCSF = (p["kin_GCSF"] * (1.0 + 3.0 * _sat(g("IL1B") + g("TNF"), 40.0))
             - p["kout_GCSF"] * GCSF - p["kcl_GCSF_N"] * GCSF * g("ANC"))
    dMONO = p["MONO0"] * p["kout_MONO"] * (1.0 + 1.5 * _sat(g("IL6"), 50.0)) \
        - p["kout_MONO"] * g("MONO")

    # ---------------- the focus ---------------------------------------------
    pen = 1.0 / (1.0 + g("NEC"))
    NT = (g("NEUT_T") + p["N_resident"]) * pen
    uptake = p["k_up"] * NT * g("BACT") / (1.0 + g("BACT") / p["Km_up"])
    dBACT = (p["mu_b"] * g("BACT") * _alive(g("BACT")) * (1.0 - g("BACT") / p["Bmax"]) - uptake
             + p["k_release"] * g("BACT_i") - E_TMP * pen * g("BACT"))
    dBACTi = (uptake * s_enc + p["mu_bi"] * g("BACT_i") * _alive(g("BACT_i"))
              - p["k_release"] * g("BACT_i"))
    dNEC = (p["k_nec"] * _sat(g("BACT") + g("BACT_i") + g("HYPH"), p["Knec"])
            - p["kout_nec"] * g("NEC"))

    hyph_comp = (f_norm + (1.0 - f_norm) * (0.10 + 0.90 * phi_res)) * E_IFNn
    s_conid = s_enc
    dCONID = -p["k_up_c"] * g("MAC") * pen * g("CONID") - p["k_germ"] * g("CONID")
    dHYPH = (p["k_germ"] * g("CONID") * s_conid
             + p["mu_h"] * g("HYPH") * _alive(g("HYPH")) * (1.0 - g("HYPH") / p["Hmax"])
             - p["k_dam_h"] * NT * g("HYPH") * hyph_comp - E_AZOL * pen * g("HYPH"))
    dGM = p["k_GM"] * g("HYPH") - p["kout_GM"] * g("GM")

    dNEUT_T = (p["k_rec"] * E_IFNr * (1.0 - I_PMIG)
               * _sat(g("BACT") + g("BACT_i") + g("HYPH"), p["Krec"]) * g("ANC") / p["ANC0"]
               * 20.0 - p["kout_NEUT"] * g("NEUT_T"))
    dMAC = (p["k_rec_mono"] * _sat(g("CONID") + g("HYPH") + g("BACT"), p["Krec"])
            * g("MONO") / p["MONO0"] * 20.0 - p["kout_MAC"] * g("MAC"))

    # ---------------- efferocytosis: the resolution defect -------------------
    eff = p["effero_min"] + (1.0 - p["effero_min"]) * phi_eff
    dAPOP = (p["k_apop"] * (g("NEUT_T") + g("MAC") + p["N_tissue0"])
             - p["k_effero"] * eff * g("APOP"))
    dNET = p["k_NET"] * g("NEUT_T") * (0.25 + 0.75 * phi_eff) - p["kout_NET"] * g("NET")

    # ---------------- cytokines ----------------------------------------------
    # amp_IL1 encodes the fact that NADPH-oxidase-derived ROS RESTRAINS the
    # NLRP3 inflammasome and drives LC3-associated phagocytosis; without it the
    # same stimulus makes several-fold more IL-1beta.  This one gain is what
    # makes the disease inflammatory rather than merely infectious, and it is
    # the reason anakinra works on colitis while antibiotics do not.
    amp = 1.0 + c["amp_IL1"] * (1.0 - phi_eff)
    stim = (math.log10(1.0 + g("BACT") + g("BACT_i") + g("HYPH") + g("CONID"))
            + p["gut_stim"])
    dIL1B = (p["IL1B0"] * p["kout_IL1B"]
             + amp * (p["k_IL1_stim"] * stim * 100.0 + p["k_IL1_apop"] * g("APOP"))
             * (1.0 - I_PRED) - p["kout_IL1B"] * g("IL1B"))
    dIL18 = (p["IL180"] * p["kout_IL18"] + amp * p["k_IL1_stim"] * stim * 120.0
             * (1.0 - I_PRED) - p["kout_IL18"] * g("IL18"))
    dTNF = (p["TNF0"] * p["kout_TNF"] + p["k_TNF_stim"] * stim * 100.0 * (1.0 - I_PRED)
            - p["kout_TNF"] * g("TNF"))
    IL1_free = g("IL1B") * (1.0 - I_ANA)
    dIL6 = (p["IL60"] * p["kout_IL6"] + (p["k_IL6_TNF"] * g("TNF")
            + p["k_IL6_IL1"] * IL1_free) * (1.0 - I_PRED) - p["kout_IL6"] * g("IL6"))
    dIL17 = (p["IL170"] * p["kout_IL17"] + p["k_IL17_IL1"] * IL1_free * 100.0
             * (1.0 - I_PRED) - p["kout_IL17"] * g("IL17"))
    dIL10 = (p["IL100"] * p["kout_IL10"]
             + p["k_IL10_eff"] * p["k_effero"] * eff * g("APOP")
             - p["kout_IL10"] * g("IL10"))
    dIFNg = (p["IFNg0"] * p["kout_IFNg"] + p["k_IFNg_IL18"] * g("IL18") * 10.0
             - p["kout_IFNg"] * g("IFNg"))
    dCRP = p["k_CRP_IL6"] * g("IL6") - p["kout_CRP"] * g("CRP")

    # ---------------- granuloma, fibrosis, colitis ---------------------------
    dIL1 = max(IL1_free - p["IL1B0"], 0.0)          # INCREMENT, see note
    dIL17x = max(g("IL17") - p["IL170"], 0.0)
    dEPI = 0.35 * _sat(dIL1 + dIL17x, 60.0) * g("MAC") - 0.02 * g("EPI")
    dGRAN = (c["k_gran"] * _sat(dIL1, p["IL1_gran_EC50"]) * (1.0 + 0.01 * g("EPI"))
             - p["kout_GRAN"] * g("GRAN") * (1.0 + 2.0 * I_PRED + 2.0 * I_ANA))
    dFIB = p["k_fib"] * g("GRAN") - p["kout_FIB"] * g("FIB")
    dCOL = (p["k_COL"] * (_sat(dIL1, 25.0) * 10.0 + 0.02 * g("GRAN"))
            - p["kout_COL"] * g("COL") * (1.0 + 3.0 * I_PRED + 3.5 * I_ANA))

    # ---------------- tryptophan / kynurenine (switchable arm) ---------------
    ido = p["k_IDO"] * _sat(g("IFNg"), p["IDO_IFNg_EC50"]) * (
        phi_eff if R["arm_ido_needs_nox"] else 1.0)
    dTRP = p["kin_TRP"] - p["kout_TRP"] * g("TRP") - ido * g("TRP")
    dKYN = ido * g("TRP") - p["kout_KYN"] * g("KYN")

    # ---------------- PK -----------------------------------------------------
    dTMPg = -p["TMP_ka"] * g("TMP_g")
    dTMPc = (p["TMP_F"] * p["TMP_ka"] * g("TMP_g") - p["TMP_CL"] * TMP
             - p["TMP_Q"] * (TMP - TMPp))
    dTMPp = p["TMP_Q"] * (TMP - TMPp)
    dSMXg = -p["SMX_ka"] * g("SMX_g")
    dSMXc = p["SMX_F"] * p["SMX_ka"] * g("SMX_g") - p["SMX_CL"] * SMX

    dITZg = -p["ITZ_ka"] * g("ITZ_g")
    dITZc = (p["ITZ_F"] * p["ITZ_ka"] * g("ITZ_g") - p["ITZ_CL"] * ITZ
             - p["ITZ_Q"] * (ITZ - ITZp))
    dITZp = p["ITZ_Q"] * (ITZ - ITZp)
    dOHI = p["OHITZ_fm"] * p["ITZ_CL"] * ITZ - p["OHITZ_CL"] * OHI

    dVORg = -p["VOR_ka"] * g("VOR_g")
    dVORc = (p["VOR_F"] * p["VOR_ka"] * g("VOR_g")
             - p["VOR_Vmax"] * VOR / (p["VOR_Km"] + VOR))
    dIFNsc = -p["IFN_ka"] * g("IFN_sc")
    dIFNc = p["IFN_F"] * p["IFN_ka"] * g("IFN_sc") - p["IFN_CL"] * IFNd
    dPREDg = -p["PRED_ka"] * g("PRED_g")
    dPREDc = (p["PRED_F"] * p["PRED_ka"] * g("PRED_g")
              - p["PRED_CL"] * g("PRED_c") / p["PRED_V"])
    dANAsc = -p["ANA_ka"] * g("ANA_sc")
    dANAc = (p["ANA_F"] * p["ANA_ka"] * g("ANA_sc")
             - p["ANA_CL"] * g("ANA_c") / p["ANA_V"])

    # ---------------- corrected myelopoiesis ---------------------------------
    dCHIM = p["k_engraft"] * (R["chim_target"] - g("CHIM"))
    dVCN = p["k_engraft"] * (R["vcn_target"] - g("VCN")) - p["k_vcn_loss"] * g("VCN")

    # ---------------- outcomes ------------------------------------------------
    burden_b = math.log10(1.0 + g("BACT") + g("BACT_i"))
    burden_f = math.log10(1.0 + g("HYPH"))
    h = (p["h_base"] / 365.0
         + c["k_h_bact"] * max(burden_b - 4.0, 0.0)
         + c["k_h_fung"] * max(burden_f - 3.0, 0.0)
         + 0.00002 * g("FIB")
         + R["extra_hazard"] / 365.0)
    dSURV = -h * g("SURV")
    dSTEROID = g("PRED_c") / p["PRED_V"] * p["PRED_V"]   # mg-days of exposure
    tox = 0.30 * (ITZ + OHI) + 0.25 * VOR + 0.002 * (TMP + SMX)
    dALT = p["kout_ALT"] * p["ALT0"] * (1.0 + tox) - p["kout_ALT"] * g("ALT")

    d = np.zeros(NS)
    for nm, val in [
        ("PROG", dPROG), ("RES", dRES), ("ANC", dANC), ("GCSF", dGCSF),
        ("MONO", dMONO), ("NEUT_T", dNEUT_T), ("MAC", dMAC), ("EPI", dEPI),
        ("BACT", dBACT), ("BACT_i", dBACTi), ("NEC", dNEC),
        ("CONID", dCONID), ("HYPH", dHYPH),
        ("GM", dGM), ("IL1B", dIL1B), ("IL18", dIL18), ("TNF", dTNF),
        ("IL6", dIL6), ("IL17", dIL17), ("IL10", dIL10), ("IFNg", dIFNg),
        ("CRP", dCRP), ("APOP", dAPOP), ("NET", dNET), ("GRAN", dGRAN),
        ("FIB", dFIB), ("COL", dCOL), ("TRP", dTRP), ("KYN", dKYN),
        ("TMP_g", dTMPg), ("TMP_c", dTMPc), ("TMP_p", dTMPp),
        ("SMX_g", dSMXg), ("SMX_c", dSMXc),
        ("ITZ_g", dITZg), ("ITZ_c", dITZc), ("ITZ_p", dITZp), ("OHITZ", dOHI),
        ("VOR_g", dVORg), ("VOR_c", dVORc), ("IFN_sc", dIFNsc), ("IFN_c", dIFNc),
        ("PRED_g", dPREDg), ("PRED_c", dPREDc), ("ANA_sc", dANAsc), ("ANA_c", dANAc),
        ("CHIM", dCHIM), ("VCN", dVCN),
        ("CUMBACT", max(burden_b - 4.0, 0.0)), ("CUMFUNG", max(burden_f - 3.0, 0.0)),
        ("SURV", dSURV), ("STEROID", dSTEROID), ("ALT", dALT),
    ]:
        d[S[nm]] = val
    return d


def initial_state(P_=None):
    p = P_ or P
    y = np.zeros(NS)
    y[S["PROG"]] = p["ANC0"] * p["kout_ANC"] * p["tau_prog"]
    y[S["RES"]]  = p["ANC0"] * p["kout_ANC"] * p["tau_res"]
    y[S["ANC"]]  = p["ANC0"]
    y[S["GCSF"]] = p["GCSF0"]
    y[S["MONO"]] = p["MONO0"]
    y[S["IL1B"]] = p["IL1B0"]; y[S["IL18"]] = p["IL180"]
    y[S["TNF"]]  = p["TNF0"];  y[S["IL6"]]  = p["IL60"]
    y[S["IL17"]] = p["IL170"]; y[S["IL10"]] = p["IL100"]
    y[S["IFNg"]] = p["IFNg0"]; y[S["CRP"]]  = p["CRP0"]
    y[S["TRP"]]  = p["TRP0"]
    y[S["SURV"]] = 1.0
    y[S["ALT"]]  = 25.0
    return y


# ---- dosing ---------------------------------------------------------------
def build_doses(regimen, days, wt=30.0, bsa=1.1):
    """Return a list of (time_d, state_name, amount) dosing events."""
    ev = []
    def rep(start, interval, n, name, amt):
        for i in range(int(n)):
            ev.append((start + i * interval, name, amt))
    nd = int(days)
    if regimen.get("tmpsmx"):
        rep(0.0, 0.5, nd * 2, "TMP_g", 5.0 * wt / 2.0)
        rep(0.0, 0.5, nd * 2, "SMX_g", 25.0 * wt / 2.0)
    if regimen.get("itra"):
        rep(0.0, 0.5, nd * 2, "ITZ_g", 100.0 if wt < 40 else 200.0)
    if regimen.get("vori"):
        rep(0.0, 0.5, nd * 2, "VOR_g", 9.0 * wt / 2.0 * 2)
    if regimen.get("ifng"):
        # 50 ug/m2 SC three times weekly, on days 0, 2, 4 of each week
        for w in range(int(days / 7) + 1):
            for off in (0.0, 2.0, 4.0):
                tt = w * 7.0 + off
                if tt <= days:
                    ev.append((tt, "IFN_sc", 50.0 * bsa))
    if regimen.get("pred"):
        rep(regimen.get("pred_start", 0.0), 1.0,
            regimen.get("pred_days", nd), "PRED_g", regimen.get("pred_dose", 1.0) * wt)
    if regimen.get("anakinra"):
        rep(regimen.get("ana_start", 0.0), 1.0,
            regimen.get("ana_days", nd), "ANA_sc", 2.0 * wt)
    ev.sort(key=lambda e: e[0])
    return ev


def simulate(genotype="X-CGD null (gp91phox-)", days=365.0, regimen=None,
             hp=None, cal=None, seed_bact=0.0, seed_conid=0.0,
             chim_target=0.0, vcn_target=0.0, extra_hazard=0.0,
             arm_ido_needs_nox=False, p=None, n_out=None, f_normal=None,
             phi_res=None, wt=30.0):
    p = p or P
    cal = cal or CALW
    regimen = regimen or {}
    if f_normal is None or phi_res is None:
        f_normal, phi_res = GENOTYPES[genotype]
    GCSF = p["GCSF0"]
    rel0 = 1.0 + p["GCSF_Emax"] * GCSF / (p["GCSF_EC50"] + GCSF)
    R = dict(P=p, CAL=cal, hp=hp, f_normal=f_normal, phi_res=phi_res,
             chim_target=chim_target, vcn_target=vcn_target,
             extra_hazard=extra_hazard, arm_ido_needs_nox=arm_ido_needs_nox,
             rel0=rel0)

    y = initial_state(p)
    y[S["BACT"]] = seed_bact
    y[S["CONID"]] = seed_conid

    ev = build_doses(regimen, days, wt=wt)
    n_out = n_out or int(days) + 1
    t_grid = np.linspace(0.0, days, n_out)
    times = sorted(set([0.0, days] + [e[0] for e in ev] + list(t_grid)))
    out_t, out_y = [0.0], [y.copy()]
    ptr = 0
    for i in range(len(times) - 1):
        t0, t1 = times[i], times[i + 1]
        while ptr < len(ev) and abs(ev[ptr][0] - t0) < 1e-9:
            y[S[ev[ptr][1]]] += ev[ptr][2]
            ptr += 1
        if t1 - t0 < 1e-12:
            continue
        sol = solve_ivp(patient_rhs, (t0, t1), y, args=(R,), method="LSODA",
                        rtol=1e-6, atol=1e-8, max_step=0.25)
        if not sol.success:
            raise RuntimeError(f"patient integration failed at t={t0}: {sol.message}")
        y = np.maximum(sol.y[:, -1], 0.0)
        for nm in ("BACT", "BACT_i", "CONID", "HYPH"):
            if y[S[nm]] < 1e-9:
                y[S[nm]] = 0.0
        out_t.append(t1); out_y.append(y.copy())
    T = np.array(out_t); Y = np.array(out_y).T
    idx = np.searchsorted(T, t_grid)
    idx = np.clip(idx, 0, len(T) - 1)
    return dict(t=T[idx], y=Y[:, idx], S=S, R=R, full_t=T, full_y=Y)


def get(res, name):
    return res["y"][S[name]]
