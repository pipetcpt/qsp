#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Niemann-Pick disease type C (NPC) — QSP reference implementation
니만-피크병 C형 · 파이썬 참조 구현 (mrgsolve 모델의 수치적 쌍둥이)
================================================================================

WHY THIS FILE EXISTS
--------------------
`npc_mrgsolve_model.R` is the deliverable; this file is its independently
written twin.  Every parameter, every right-hand side and every derived output
is duplicated here so that (a) the numbers quoted in README.md are produced by
code that actually ran, and (b) an error in one implementation shows up as a
disagreement rather than as a plausible-looking curve.  Anything reported in
the README as a model result is printed by this script.

THE THREE STRUCTURAL CLAIMS BEING TESTED
----------------------------------------
(1) COMPARTMENT CLAIM.  Plasma cholestane-3beta,5alpha,6beta-triol is generated
    by SATURABLE oxidation of stored sterol in the VISCERAL pool.  Both the
    patient mean (88.31 ng/mL) and the control mean (5.97 ng/mL) of
    PMID 33228797 are hit exactly by a two-point calibration.  The model then
    PREDICTS, rather than assumes, that triol is a poor severity marker,
    because in patients it sits at ~82% of its own ceiling.

(2) RESERVE CLAIM.  Neurological function is split as
        FUNC = 1 - D_rev - D_irr
    with D_rev reversible on a ~2-week time constant and D_irr monotone
    non-decreasing.  A symptomatic drug moves D_rev only; a disease-modifying
    drug lowers dD_irr/dt.

(3) DESIGN CLAIM.  Consequence of (2): a 12-week crossover on SARA can only see
    D_rev, and a 12-month change-from-baseline on NPCCSS cannot separate a
    constant offset from a slope change.  Both approved drugs were tested with
    the design matched to their own assumed mechanism.  This script runs each
    trial design against BOTH mechanisms and reports what each design would
    have concluded.

CALIBRATION TARGETS (all from cited, individually verified publications)
-----------------------------------------------------------------------
  T1  plasma triol, NPC patients          88.31 ng/mL   PMID 33228797
  T2  plasma triol, controls               5.97 ng/mL   PMID 33228797
  T3  5-domain NPCCSS natural progression  1.5 pt/yr    PMID 33228797
  T4  17-domain NPCCSS progression         2.7-2.9 /yr  PMID 33228797
  T5  arimoclomol 12-mo 5D difference     -1.40 pt      PMID 34418116
  T6  arimoclomol 12-mo R4D difference    -1.70 pt      PMID 40520915
  T7  arimoclomol placebo arm 5D           +2.11 pt/yr  PMID 34418116
  T8  levacetylleucine 12-wk SARA diff    -1.28 pt      PMID 38294974
  T9  triol vs 5D NPCCSS Spearman rho      0.265        PMID 33228797
  T10 adrabetadex open-label progressors   7/14 vs 21/21  PMID 28803710

Run:  python3 npc_reference_model.py
      python3 npc_reference_model.py --quick     (coarser dt, ~4x faster)
================================================================================
"""

import math
import sys
import argparse
from collections import OrderedDict

import numpy as np

# =============================================================================
# 0.  PARAMETERS  —  this dict is the single source of truth.  The mrgsolve
#     $PARAM block is generated from it (see --emit-param) and must match.
# =============================================================================

P = OrderedDict([

    # ---- 1. genotype -> NPC1 protein --------------------------------------
    ("ksyn_npc1",   1.00),   # NPC1 synthesis, arbitrary units/d
    ("kerad",       2.00),   # ER-associated degradation rate, /d
    ("kexit0",      1.00),   # max ER exit rate of correctly folded NPC1, /d
    ("kdegL",       0.15),   # lysosomal NPC1 turnover, /d (t1/2 ~4.6 d)
    ("EChsp",       0.80),   # HSP70 fold-change above 1 giving half-max folding gain
    ("Emax_fold",   0.300),  # max fractional closure of the folding gap by HSP70
                             # (gives a 2.8x rise in functional NPC1 -- see README:
                             # even the arithmetic CEILING of this mechanism cannot
                             # reach the published arimoclomol effect size)
    ("f_leak",      0.010),  # NPC1-independent lysosomal cholesterol leak (fraction of WT)

    # ---- 2. cholesterol compartments --------------------------------------
    ("Jin_v",       1.00),   # visceral lysosomal cholesterol influx, units/d
    ("Vmax_v",      3.00),   # visceral NPC1-mediated egress Vmax
    ("Km_v",        1.00),   # visceral egress Km
    # The CNS pool runs on a NEURONAL timescale.  Brain cholesterol turns over
    # with a half-life of months-to-years (24S-hydroxycholesterol efflux is the
    # only quantitatively important route), whereas the hepatic pool turns over
    # in days.  Jin_c, Vmax_c and kbas_chol_c are the visceral values divided by
    # tau_cns = 25, which leaves every CNS STEADY STATE unchanged and makes the
    # APPROACH to it 25x slower.  This one factor is what creates the years-long
    # latent period between birth and neurological onset.
    ("tau_cns",    25.00),   # documentation only; the divisions are baked in below
    ("Jin_c",     0.01400),  # = 0.35 / 25
    ("Vmax_c",    0.04200),  # = 1.05 / 25
    ("Km_c",        1.00),
    ("kbas_chol",   0.020),  # visceral non-specific cholesterol clearance, /d
    ("kbas_chol_c",0.00080), # = 0.020 / 25  (CNS)

    # ---- 3. glycosphingolipid ---------------------------------------------
    ("ksyn_gsl_v",  0.60),   # visceral GSL synthesis (UGCG flux)
    ("ksyn_gsl_c", 0.01200), # CNS GSL synthesis  (= 0.30 / 25)
    ("kdeg_gsl",    0.60),   # visceral GSL catabolism, /d (x hydrolase activity)
    ("kdeg_gsl_c", 0.02400), # CNS GSL catabolism (= 0.60 / 25)

    # ---- 4. sphingosine / lysosomal calcium -------------------------------
    ("ksph",        1.00),   # sphingosine generation, units/d
    ("kexp_sph",    4.00),   # sphingosine export, /d
    ("fmin_sph",    0.12),   # NPC1-independent floor of sphingosine export
    ("Ksph",        1.15),   # sphingosine IC50 on acidic Ca2+ store refilling
    ("Ktr_chol",   60.00),   # cholesterol IC50 on TRPML1
    ("kca",         5.00),   # Ca2+ store equilibration rate, /d

    # ---- 5. hydrolase capacity --------------------------------------------
    ("Kph",        14.00),   # cholesterol load raising lysosomal pH (CNS units)
    ("hill_ph",     1.50),
    ("e_hyd_hsp",   0.20),   # HSP70 gain on hydrolase activity
    ("khyd",        3.00),   # hydrolase equilibration, /d

    # ---- 6. autophagy backlog ---------------------------------------------
    ("kaut_in",     1.00),
    ("kaut_out",    2.50),

    # ---- 7. mitochondria / ROS --------------------------------------------
    ("Kmito",      22.00),   # cholesterol load halving mitochondrial function
    ("kmito",       0.50),
    ("aros",        2.50),   # ROS gain per unit lost mitochondrial function
    ("kros",        4.00),

    # ---- 8. Purkinje cell pools -------------------------------------------
    ("S_thresh",    0.200),  # stress absorbed without any damage accrual
    ("kon_pc",      0.020),  # healthy -> stressed, per unit stress per day
    ("koff_pc",     0.010),  # stressed -> healthy (recovery)
    ("Krec",        0.50),   # stress that halves recovery
    # ---- THE RESERVE, WRITTEN AS AN INTEGRAL --------------------------------
    # DAM = integral of stress over time, with NO repair term.  Purkinje death
    # is gated on DAM crossing D_reserve.  Two consequences fall out of this one
    # choice, and both are published facts the model was not fitted to:
    #   * latency to onset is inversely proportional to stress, so residual NPC1
    #     activity spreads onset over decades (juvenile ~8 yr, adult ~25+ yr);
    #   * the post-gate death rate SATURATES, so the slope after onset is nearly
    #     independent of the age of onset -- which is exactly the headline of
    #     Yanjanin 2010 (PMID 19415691), "linear clinical progression,
    #     independent of age of onset".
    # A repair term would make the crossing time depend only LOGARITHMICALLY on
    # stress, collapsing the onset range; that version was tried and rejected.
    ("D_reserve", 6736.55),  # stress x days of damage the cerebellum absorbs
    ("hill_gate",   4.00),   # sharpness of the reserve threshold
    ("kdie",    1.13537e-4),  # stressed -> dead, /d once the gate is open
    # developmental vulnerability: the same biochemical insult costs more during
    # active cerebellar development and myelination (PMID 15217094).  This is a
    # PATIENT-LEVEL covariate, not a genotype consequence -- residual NPC1
    # activity alone does not determine the onset form, and siblings with
    # identical genotypes can differ by a decade.
    ("v_dev",       0.600),
    ("tau_dev",     2.500),  # years
    ("w_chol",      0.55),   # stress weights (normalised, healthy S = 0)
    ("w_gsl",       0.22),
    ("w_ros",       0.30),
    ("w_aut",       0.18),
    ("w_infl",      0.60),
    ("w_ca",        0.35),

    # ---- 9. neuroinflammation ---------------------------------------------
    ("a_chol",      0.35),
    ("a_pc",        0.90),
    ("a_dead",      0.55),
    ("kinfl",       0.050),

    # ---- 10. synapses / cerebellar volume ---------------------------------
    ("ksyn_rep",    0.020),
    ("kprune",      0.020),
    ("g_vol_pc",    0.85),
    ("g_vol_syn",   0.25),
    ("kcbl",        0.0040),

    # ---- 11. function decomposition (THE structural core) -----------------
    # Reversible dysfunction is a GRADED, saturating function of the CURRENT
    # biochemical stress (which already contains inflammation, ROS, calcium and
    # autophagic backlog).  It is deliberately NOT a function of how many cells
    # are in a "stressed" pool: a saturating cell-count would make D_rev
    # saturate within weeks of birth, which no NPC patient does.
    ("cap_rev",     0.335),  # ceiling on reversible dysfunction
    ("K_rev",       6.000),  # stress giving half of cap_rev
    ("cr_pcs",      0.000),  # deliberately zero: a saturating cell-count would
                             # make D_rev plateau within weeks of birth
    ("krev",        0.060),  # D_rev time constant, /d  (t1/2 ~11.6 d)
    ("g_irr_pc",    1.60),   # irreversible from dead PC
    ("g_irr_syn",   0.00),   # synapse density and cerebellar volume are READOUTS
                             # of Purkinje loss here, not independent drivers
    ("g_irr_cbl",   0.00),   # cerebellar volume is a READOUT of PC loss, not an
                             # independent driver -- including both double counts
    ("kirr",        0.010),  # D_irr tracking rate, /d (one-sided)

    # ---- 12. clinical scale mapping ---------------------------------------
    # SARA scores what the patient can do RIGHT NOW under cerebellar control, so
    # it is weighted heavily towards D_rev.  The NPCCSS domains (ambulation,
    # speech, swallow, fine motor, cognition) are milestone-like and are
    # weighted towards D_irr.  This asymmetry is not cosmetic: it is why a
    # symptomatic drug can win a 12-week SARA endpoint while a rate endpoint
    # needs NPCCSS and a year (see section 9).
    ("sara_max",   40.00),
    ("sara_a",      1.80),   # weight of D_rev in SARA
    ("sara_b",      1.00),   # weight of D_irr in SARA
    ("sara_p",      1.00),   # linear: NPC progression is reported as linear
    ("sara_k",     30.3141),# scale (calibrated: SARA = 15.91 at trial entry)
    ("n5_max",     25.00),
    ("n5_a",        0.25),   # weight of D_rev in NPCCSS
    ("n5_b",        1.00),
    ("n5_p",        1.00),
    ("n5_k",       30.00),
    ("n4_frac",     0.80),   # 4-domain is 4/5 of the domains ...
    ("n4_gain",     1.16),   # ... but rescored for linearity (PMID 40520915)
    ("n17_k",       1.95),   # 17-domain = n17_k x 5-domain + hearing/seizure/etc.
    ("sacc_rev",    0.40),   # saccade velocity sensitivity to D_rev
    ("sacc_irr",    1.60),   # ... and to D_irr
    ("swal_a",      0.55),
    ("swal_b",      1.55),

    # ---- 13. survival hazard ----------------------------------------------
    ("h0",        2.0e-5),   # baseline hazard, /d
    ("h_swal",  1.21875e-4), # hazard from dysphagia (aspiration; calibrated to
                             # a median survival of 26 yr in the juvenile form)
    ("h_liver",   4.0e-4),   # hazard from infantile liver disease

    # ---- 14. biomarkers ---------------------------------------------------
    ("ktri",       54.15),   # triol generation scale (calibrated, T1+T2)
    ("Ktri",        8.35),   # triol generation SATURATION constant (calibrated)
    ("ktri_el",     0.50),   # triol elimination, /d
    ("kppcs",       6.20),   # lysoSM-509 / PPCS generation
    ("Kppcs",      10.00),
    ("b_ppcs_gsl",  0.45),
    ("kppcs_el",    0.35),
    ("ktcg",       12.00),   # bile acid B (TCG)
    ("Ktcg",        9.00),
    ("ktcg_el",     0.60),
    ("knfl",      420.00),   # CSF NfL from Purkinje death flux
    ("knfl_base",   1.20),
    ("knfl_el",     0.14),
    ("kcalb",     900.00),   # CSF calbindin from Purkinje death flux
    ("kcalb_el",    0.30),

    # ---- 15. miglustat PK/PD ----------------------------------------------
    ("mig_ka",      8.00),   # /d
    ("mig_F",       0.97),
    ("mig_V",      90.00),   # L
    ("mig_CL",    230.00),   # L/d  -> t1/2 6.5 h
    ("mig_Q",      20.00),
    ("mig_Vp",     60.00),
    ("mig_Kp_br",   0.45),   # brain:plasma partition
    ("mig_ke0",     1.00),   # brain effect-site equilibration, /d
    ("mig_IC50",   30.00),   # uM, UGCG inhibition
    ("mig_Emax",    0.85),
    ("mig_MW",    219.28),

    # ---- 16. arimoclomol PK/PD --------------------------------------------
    ("ari_ka",     12.00),
    ("ari_F",       0.60),
    ("ari_V",     150.00),   # L
    ("ari_CL",    624.00),   # L/d -> t1/2 ~4 h
    ("ari_Kp_csf",  0.15),
    ("ari_ke0",     2.00),
    ("ari_EC50",  120.00),   # ng/mL in CSF
    ("ari_Emax",    1.60),   # max HSP70 fold-increase above baseline
    ("ari_stress0", 0.25),   # co-inducer: fraction of Emax available without stress

    # ---- 17. levacetylleucine PK/PD ---------------------------------------
    ("nal_ka",     30.00),
    ("nal_F",       1.00),
    ("nal_V",     253.00),   # L   (label: Vss/F 253 L)
    ("nal_CL",   3336.00),   # L/d (label: 139 L/h)
    ("nal_ke0",     0.15),   # effect-site equilibration, /d (t1/2 ~4.6 d)
    ("nal_EC50",    0.80),   # mg/L at effect site
    ("nal_Emax_sym",0.30779),  # max fractional reduction of D_rev target
    ("nal_e_gluc",  0.120),  # max gain on brain glucose metabolism readout
    ("nal_Emax_dm", 0.080),  # max fractional reduction of CNS lipid influx

    # ---- 18. adrabetadex (2-HPbCD) PK/PD ----------------------------------
    ("cd_Vcsf",     0.15),   # L
    ("cd_kout",     4.16),   # CSF -> systemic, /d (t1/2 4 h)
    ("cd_kin_br",   0.60),   # CSF -> brain ECF/intracellular, /d
    ("cd_kout_br",  1.50),   # brain ECF elimination, /d
    ("cd_Vsys",    18.00),   # L
    ("cd_klsys",   24.00),   # /d systemic elimination (renal, rapid)
    ("cd_kext",    3.85e-4), # cholesterol extraction, per (mg/L) per day.  Anchored so
                             # that intrathecal dosing NORMALISES the neuronal pool
                             # (~1.2x wild type) rather than emptying it.  The first
                             # value tried, 0.0130, drove CNS cholesterol to 0.06x
                             # wild type -- arithmetically fine, since stress is
                             # floored at zero, and biologically absurd.  Efficacy is
                             # almost unchanged because the stress term saturates
                             # once the pool is back in the normal range.
    ("cd_koto",   1.94e-5),  # outer-hair-cell loss, per (mg/L in CSF) per day.
                             # Anchored to a ~25 dB high-frequency threshold shift at
                             # 12 months, which is the clinical order of magnitude;
                             # the first value tried destroyed every outer hair cell
                             # within a few months.
    ("cd_hear_max",55.00),   # dB threshold shift at total OHC loss
    ("cd_hear_p",   0.70),
])

# genotype table: (f_null fraction of alleles, theta0 foldability, f_npc2)
GENOTYPES = {
    "WT":            dict(f_null=0.00, theta0=1.000, f_npc2=1.00),
    "I1061T/I1061T": dict(f_null=0.00, theta0=0.055, f_npc2=1.00),  # juvenile
    "I1061T/null":   dict(f_null=0.50, theta0=0.055, f_npc2=1.00),  # late infantile
    "null/null":     dict(f_null=1.00, theta0=0.055, f_npc2=1.00),  # perinatal
    "mild/mild":     dict(f_null=0.00, theta0=0.322, f_npc2=1.00),  # adult onset
    "NPC2":          dict(f_null=0.00, theta0=1.000, f_npc2=0.03),  # NPC2 disease
}

# =============================================================================
# 1.  STATE VECTOR
# =============================================================================

STATES = [
    # drug PK (13)
    "MIG_GUT", "MIG_CEN", "MIG_PER", "MIG_BR",
    "ARI_GUT", "ARI_CEN", "ARI_CSF",
    "NAL_GUT", "NAL_CEN", "NAL_BR",
    "CD_CSF", "CD_BR", "CD_SYS",
    # NPC1 protein (2)
    "NPC1_ER", "NPC1_L",
    # visceral (2)
    "CHOL_V", "GSL_V",
    # CNS lipids (4)
    "CHOL_C", "GSL_C", "SPH", "CA_LY",
    # lysosomal / metabolic function (4)
    "HYD", "AUTOPH", "MITO", "ROS",
    # Purkinje pools + network (5)
    "PC", "PC_S", "PC_LOST", "INFL", "SYN",
    "CBL", "DAM",
    # function decomposition (2)
    "D_REV", "D_IRR",
    # biomarkers (5)
    "TRIOL", "PPCS", "TCG", "NFL", "CALB",
    # cochlea (1)
    "OHC",
    # cumulative hazard (1)
    "CUMHAZ",
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)

# integer index constants (I_CHOL_C etc.) -- the RHS is called ~10^7 times, so
# every dict lookup removed from it matters.
for _s, _i in IX.items():
    globals()["I_" + _s] = _i


class Par:
    """Parameter bundle with attribute access (fast) instead of dict lookup."""
    __slots__ = list(P.keys())

    def __init__(self, d):
        for k in P:
            setattr(self, k, float(d[k]))

    def copy_with(self, **kw):
        d = {k: getattr(self, k) for k in P}
        d.update(kw)
        return Par(d)

    def as_dict(self):
        return {k: getattr(self, k) for k in P}


# =============================================================================
# 2.  REFERENCE (WILD-TYPE) STEADY STATE
#     Every normalised quantity in the RHS is divided by its WT value so that
#     a healthy system sits at stress S = 0 exactly.
# =============================================================================

def wt_reference(p):
    """Analytic wild-type steady state of the untreated lipid/lysosome core."""
    g = GENOTYPES["WT"]
    # NPC1
    theta = g["theta0"]
    er = p.ksyn_npc1 * (1 - g["f_null"]) / (p.kexit0 * theta + p.kerad)
    npc1L = er * p.kexit0 * theta / p.kdegL
    f_npc1 = 1.0                      # by definition of the reference
    f_eg = 1.0 * g["f_npc2"] + p.f_leak

    def chol_ss(Jin, Vmax, Km, kbas):
        # Jin = Vmax*f_eg*x/(Km+x) + kbas*x   ->  quadratic in x
        a = kbas
        b = Vmax * f_eg + kbas * Km - Jin
        c = -Jin * Km
        return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)

    chol_v = chol_ss(p.Jin_v, p.Vmax_v, p.Km_v, p.kbas_chol)
    chol_c = chol_ss(p.Jin_c, p.Vmax_c, p.Km_c, p.kbas_chol_c)

    sph = p.ksph / (p.kexp_sph * (p.fmin_sph + (1 - p.fmin_sph) * f_npc1))
    ca = (1.0 / (1.0 + (sph / p.Ksph) ** 2)) * (1.0 / (1.0 + chol_c / p.Ktr_chol))
    hyd = (1.0 / (1.0 + (chol_c / p.Kph) ** p.hill_ph)) * (0.45 + 0.55 * ca)
    gsl_c = p.ksyn_gsl_c / (p.kdeg_gsl_c * hyd)
    gsl_v = p.ksyn_gsl_v / (p.kdeg_gsl * hyd)
    aut = p.kaut_in / (p.kaut_out * hyd * (0.5 + 0.5 * ca))
    mito = 1.0 / (1.0 + (chol_c / p.Kmito) ** 2)
    ros = 1.0 + p.aros * (1.0 / mito - 1.0)

    return dict(NPC1_ER=er, NPC1_L=npc1L, CHOL_V=chol_v, CHOL_C=chol_c,
                SPH=sph, CA=ca, HYD=hyd, GSL_C=gsl_c, GSL_V=gsl_v,
                AUT=aut, MITO=mito, ROS=ros, NPC1_L_WT=npc1L)


PAR = Par(P)
WT = wt_reference(PAR)


# =============================================================================
# 3.  DOSING
# =============================================================================

class Regimen:
    """A drug regimen: which compartment, dose amount, interval, start, stop."""

    def __init__(self, target, amt, ii, start=0.0, stop=1e9):
        self.target = target      # state name receiving the dose
        self.amt = amt            # mg
        self.ii = ii              # days between doses
        self.start = start
        self.stop = stop

    def doses_in(self, t0, t1):
        """Total amount delivered in the window (t0, t1]."""
        if self.amt == 0 or t1 <= self.start or t0 >= self.stop:
            return 0.0
        lo = max(t0, self.start)
        hi = min(t1, self.stop)
        if hi <= lo:
            return 0.0
        n0 = math.floor((lo - self.start) / self.ii + 1e-9)
        n1 = math.floor((hi - self.start) / self.ii + 1e-9)
        n = n1 - n0
        # include the dose exactly at self.start when the window opens there
        if lo <= self.start + 1e-9 < hi or (t0 <= self.start < t1):
            n += 1
        return n * self.amt


def mig_regimen(mg_per_dose=200.0, start=0.0, stop=1e9):
    return Regimen("MIG_GUT", mg_per_dose, 1.0 / 3.0, start, stop)


def ari_regimen(mg_per_dose=124.0, start=0.0, stop=1e9):
    return Regimen("ARI_GUT", mg_per_dose, 1.0 / 3.0, start, stop)


def nal_regimen(g_per_day=4.0, start=0.0, stop=1e9):
    return Regimen("NAL_GUT", g_per_day * 1000.0 / 3.0, 1.0 / 3.0, start, stop)


def cd_regimen(mg=900.0, weeks=2.0, start=0.0, stop=1e9):
    return Regimen("CD_CSF", mg, 7.0 * weeks, start, stop)


# =============================================================================
# 4.  RIGHT-HAND SIDE
# =============================================================================

def rhs(t, y, ctx):
    p = ctx["pr"]                 # Par object: attribute access
    g = ctx["geno"]
    # plain Python floats, not a numpy array: this RHS is evaluated ~10^7 times
    # per full report and numpy scalar arithmetic is ~4x slower than float.
    d = [0.0] * NS

    # ---------------- drug PK ---------------------------------------------
    mig_gut, mig_cen, mig_per, mig_br = y[I_MIG_GUT], y[I_MIG_CEN], y[I_MIG_PER], y[I_MIG_BR]
    a_mig = p.mig_ka * mig_gut
    cl_mig = p.mig_CL / p.mig_V
    q_in = p.mig_Q * (mig_cen / p.mig_V - mig_per / p.mig_Vp)
    d[I_MIG_GUT] = -a_mig
    d[I_MIG_CEN] = p.mig_F * a_mig - cl_mig * mig_cen - q_in
    d[I_MIG_PER] = q_in
    C_mig_pl = mig_cen / p.mig_V / p.mig_MW * 1000.0            # uM
    d[I_MIG_BR] = p.mig_ke0 * (p.mig_Kp_br * C_mig_pl - mig_br)
    C_mig_br = mig_br                                                 # uM

    ari_gut, ari_cen, ari_csf = y[I_ARI_GUT], y[I_ARI_CEN], y[I_ARI_CSF]
    a_ari = p.ari_ka * ari_gut
    d[I_ARI_GUT] = -a_ari
    d[I_ARI_CEN] = p.ari_F * a_ari - p.ari_CL / p.ari_V * ari_cen
    C_ari_pl = ari_cen / p.ari_V * 1e6 / 1000.0                    # ng/mL
    d[I_ARI_CSF] = p.ari_ke0 * (p.ari_Kp_csf * C_ari_pl - ari_csf)
    C_ari_csf = ari_csf

    nal_gut, nal_cen, nal_br = y[I_NAL_GUT], y[I_NAL_CEN], y[I_NAL_BR]
    a_nal = p.nal_ka * nal_gut
    d[I_NAL_GUT] = -a_nal
    d[I_NAL_CEN] = p.nal_F * a_nal - p.nal_CL / p.nal_V * nal_cen
    C_nal_pl = nal_cen / p.nal_V                                   # mg/L
    d[I_NAL_BR] = p.nal_ke0 * (C_nal_pl - nal_br)
    C_nal_br = nal_br

    cd_csf, cd_br, cd_sys = y[I_CD_CSF], y[I_CD_BR], y[I_CD_SYS]
    d[I_CD_CSF] = -(p.cd_kout + p.cd_kin_br) * cd_csf
    d[I_CD_BR] = p.cd_kin_br * cd_csf - p.cd_kout_br * cd_br
    d[I_CD_SYS] = p.cd_kout * cd_csf - p.cd_klsys * cd_sys
    C_cd_csf = cd_csf / p.cd_Vcsf                                  # mg/L
    C_cd_br = cd_br / p.cd_Vcsf                                    # mg/L (same nominal volume)
    C_cd_sys = cd_sys / p.cd_Vsys

    # ---------------- drug effects ----------------------------------------
    # miglustat: UGCG inhibition, separately in viscera (plasma) and CNS (brain)
    I_mig_v = p.mig_Emax * C_mig_pl / (p.mig_IC50 + C_mig_pl)
    I_mig_c = p.mig_Emax * C_mig_br / (p.mig_IC50 + C_mig_br)

    # arimoclomol: HSF1 CO-inducer -> amplification requires cellular stress.
    # Stress proxy = normalised CNS/visceral lipid load, capped at 1.
    chol_c, chol_v = y[I_CHOL_C], y[I_CHOL_V]
    stress_cell = min(1.0, p.ari_stress0 +
                      (1 - p.ari_stress0) * min(1.0, (chol_c / WT["CHOL_C"] - 1.0) / 12.0))
    HSP70 = 1.0 + p.ari_Emax * stress_cell * C_ari_csf / (p.ari_EC50 + C_ari_csf)

    # levacetylleucine: three separable actions
    E_nal = C_nal_br / (p.nal_EC50 + C_nal_br)
    E_nal_sym = p.nal_Emax_sym * E_nal * ctx["nal_sym_on"]
    E_nal_gluc = p.nal_e_gluc * E_nal * ctx["nal_sym_on"]
    E_nal_dm = p.nal_Emax_dm * E_nal * ctx["nal_dm_on"]

    # ---------------- NPC1 protein ----------------------------------------
    hsp_drive = max(0.0, HSP70 - 1.0)
    theta = g["theta0"] + (1 - g["theta0"]) * p.Emax_fold * \
        hsp_drive / (p.EChsp + hsp_drive)
    if theta > 1.0:            # a folding yield is a fraction, hard bound
        theta = 1.0
    er, npc1L = y[I_NPC1_ER], y[I_NPC1_L]
    kexit = p.kexit0 * theta
    d[I_NPC1_ER] = p.ksyn_npc1 * (1 - g["f_null"]) - er * (kexit + p.kerad)
    d[I_NPC1_L] = er * kexit - p.kdegL * npc1L
    f_npc1 = npc1L / WT["NPC1_L_WT"]
    f_eg = f_npc1 * g["f_npc2"] + p.f_leak

    # ---------------- cholesterol pools -----------------------------------
    eg_v = p.Vmax_v * f_eg * chol_v / (p.Km_v + chol_v)
    eg_c = p.Vmax_c * f_eg * chol_c / (p.Km_c + chol_c)
    d[I_CHOL_V] = (p.Jin_v - eg_v - p.kbas_chol * chol_v
                       - p.cd_kext * C_cd_sys * chol_v)
    d[I_CHOL_C] = (p.Jin_c * (1 - E_nal_dm) - eg_c - p.kbas_chol_c * chol_c
                       - p.cd_kext * C_cd_br * chol_c)

    # ---------------- sphingosine, calcium, hydrolases ---------------------
    sph = y[I_SPH]
    d[I_SPH] = p.ksph - p.kexp_sph * (p.fmin_sph +
                                                (1 - p.fmin_sph) * f_npc1) * sph
    ca = y[I_CA_LY]
    ca_ss = (1.0 / (1.0 + (sph / p.Ksph) ** 2)) * (1.0 / (1.0 + chol_c / p.Ktr_chol))
    d[I_CA_LY] = p.kca * (ca_ss - ca)

    hyd = y[I_HYD]
    hyd_ss = ((1.0 / (1.0 + (chol_c / p.Kph) ** p.hill_ph))
              * (1.0 + p.e_hyd_hsp * max(0.0, HSP70 - 1.0))
              * (0.45 + 0.55 * ca))
    d[I_HYD] = p.khyd * (hyd_ss - hyd)

    # ---------------- glycosphingolipids ----------------------------------
    gsl_v, gsl_c = y[I_GSL_V], y[I_GSL_C]
    d[I_GSL_V] = p.ksyn_gsl_v * (1 - I_mig_v) - p.kdeg_gsl * gsl_v * hyd
    d[I_GSL_C] = p.ksyn_gsl_c * (1 - I_mig_c) - p.kdeg_gsl_c * gsl_c * hyd

    # ---------------- autophagy, mitochondria, ROS ------------------------
    aut = y[I_AUTOPH]
    prot = hyd * (0.5 + 0.5 * ca)
    d[I_AUTOPH] = p.kaut_in - p.kaut_out * prot * aut

    mito = y[I_MITO]
    mito_ss = (1.0 / (1.0 + (chol_c / p.Kmito) ** 2)) * (1.0 + E_nal_gluc)
    d[I_MITO] = p.kmito * (mito_ss - mito)
    ros = y[I_ROS]
    ros_ss = 1.0 + p.aros * max(0.0, 1.0 / max(mito, 1e-6) - 1.0)
    d[I_ROS] = p.kros * (ros_ss - ros)

    # ---------------- Purkinje cell stress and pools ----------------------
    infl = y[I_INFL]
    stress = (p.w_chol * (chol_c / WT["CHOL_C"] - 1.0) / 10.0
              + p.w_gsl * (gsl_c / WT["GSL_C"] - 1.0)
              + p.w_ros * (ros - WT["ROS"])
              + p.w_aut * (aut / WT["AUT"] - 1.0)
              + p.w_infl * infl
              + p.w_ca * (1.0 - ca / WT["CA"]))
    # FUNCTIONAL RESERVE: only the stress above the threshold does damage.
    stress = max(0.0, stress - p.S_thresh)

    # developmental vulnerability, largest in the first years of life: the same
    # biochemical insult costs MORE while the cerebellum is still being built and
    # myelinated (PMID 15217094).  It multiplies BOTH the rate at which the
    # reserve is spent and the death rate once the gate is open -- which is what
    # lets it move the age of ONSET and not only the slope afterwards.
    vuln = 1.0 + p.v_dev * math.exp(-(t / 365.25) / p.tau_dev)

    # cumulative damage: a pure integral of weighted stress, NO repair term
    d[I_DAM] = stress * vuln
    ru = y[I_DAM] / p.D_reserve
    run = ru ** p.hill_gate
    gate = run / (1.0 + run)

    pc, pcs, pcl = y[I_PC], y[I_PC_S], y[I_PC_LOST]
    to_stress = p.kon_pc * stress * pc
    recover = p.koff_pc * pcs / (1.0 + stress / p.Krec)
    die = p.kdie * pcs * gate * vuln
    d[I_PC] = -to_stress + recover
    d[I_PC_S] = to_stress - recover - die
    d[I_PC_LOST] = die

    infl_ss = (p.a_chol * (chol_c / WT["CHOL_C"] - 1.0) / 25.0
               + p.a_pc * pcs + p.a_dead * pcl)
    d[I_INFL] = p.kinfl * (infl_ss - infl)

    syn = y[I_SYN]
    d[I_SYN] = p.ksyn_rep * ((pc + pcs) - syn) - p.kprune * infl * syn

    cbl = y[I_CBL]
    cbl_t = 1.0 - p.g_vol_pc * pcl - p.g_vol_syn * (1.0 - syn)
    d[I_CBL] = p.kcbl * (cbl_t - cbl)

    # ---------------- function decomposition ------------------------------
    glucm = mito / WT["MITO"]
    tgt_rev = (p.cap_rev * stress / (p.K_rev + stress) + p.cr_pcs * pcs)
    tgt_rev = min(0.90, max(0.0, tgt_rev)) * (1.0 - E_nal_sym)
    d[I_D_REV] = p.krev * (tgt_rev - y[I_D_REV])

    # LINEAR in accumulated cell loss: NPC progression on the NPCCSS is reported
    # as linear (PMID 19415691), and a saturating form cannot produce that.
    tgt_irr = min(0.95, p.g_irr_pc * pcl + p.g_irr_syn * (1.0 - syn)
                  + p.g_irr_cbl * (1.0 - cbl))
    d[I_D_IRR] = p.kirr * max(0.0, tgt_irr - y[I_D_IRR])

    # ---------------- biomarkers ------------------------------------------
    d[I_TRIOL] = p.ktri * chol_v / (p.Ktri + chol_v) - p.ktri_el * y[I_TRIOL]
    d[I_PPCS] = (p.kppcs * (chol_v + p.b_ppcs_gsl * gsl_v) /
                     (p.Kppcs + chol_v + p.b_ppcs_gsl * gsl_v)
                     - p.kppcs_el * y[I_PPCS])
    d[I_TCG] = p.ktcg * chol_v / (p.Ktcg + chol_v) - p.ktcg_el * y[I_TCG]
    d[I_NFL] = p.knfl * die + p.knfl_base - p.knfl_el * y[I_NFL]
    d[I_CALB] = p.kcalb * die - p.kcalb_el * y[I_CALB]

    # ---------------- cochlea (cyclodextrin ototoxicity) -------------------
    ohc = y[I_OHC]
    d[I_OHC] = -p.cd_koto * C_cd_csf * ohc

    # ---------------- survival hazard -------------------------------------
    swal = clinical_swallow(y, p)
    haz = p.h0 + p.h_swal * (swal / 5.0) ** 2 + p.h_liver * ctx["liver_risk"]
    d[I_CUMHAZ] = haz

    return d


# =============================================================================
# 5.  DERIVED (ALGEBRAIC) OUTPUTS
# =============================================================================

def clinical_sara(y, p=PAR):
    x = p.sara_a * y[I_D_REV] + p.sara_b * y[I_D_IRR]
    return min(p.sara_max, p.sara_k * max(0.0, x) ** p.sara_p)


def clinical_n5(y, p=PAR):
    x = p.n5_a * y[I_D_REV] + p.n5_b * y[I_D_IRR]
    return min(p.n5_max, p.n5_k * max(0.0, x) ** p.n5_p)


def clinical_n4(y, p=PAR):
    """R4DNPCCSS: 4 of the 5 domains, rescored for linearity (PMID 40520915)."""
    return min(20.0, clinical_n5(y, p) * p.n4_frac * p.n4_gain)


def clinical_n17(y, p=PAR):
    """17-domain adds hearing, seizures, cataplexy, behaviour, ... (max 61)."""
    base = clinical_n5(y, p) * p.n17_k
    hearing_pts = 5.0 * min(1.0, hearing_shift(y, p) / p.cd_hear_max)
    return min(61.0, base + hearing_pts)


def clinical_swallow(y, p=PAR):
    x = p.swal_a * y[I_D_REV] + p.swal_b * y[I_D_IRR]
    return min(5.0, 5.0 * min(1.0, x))


def saccade_velocity(y, p=PAR):
    """Peak horizontal saccade velocity, % of normal."""
    return 100.0 * math.exp(-(p.sacc_rev * y[I_D_REV] +
                              p.sacc_irr * y[I_D_IRR]))


def hearing_shift(y, p=PAR):
    lost = 1.0 - y[I_OHC]
    if lost <= 0.0:
        return 0.0
    return p.cd_hear_max * lost ** p.cd_hear_p


def survival_prob(y, p=PAR):
    return math.exp(-y[I_CUMHAZ])


DERIVED = OrderedDict([
    ("SARA", clinical_sara), ("NPCCSS5", clinical_n5), ("NPCCSS4", clinical_n4),
    ("NPCCSS17", clinical_n17), ("SWALLOW", clinical_swallow),
    ("SACCADE", saccade_velocity), ("HEARING_dB", hearing_shift),
    ("SURV", survival_prob),
])


# =============================================================================
# 6.  INITIAL CONDITIONS AND INTEGRATOR
# =============================================================================

def initial_state(genotype, p=PAR):
    """
    Birth state.  Cholesterol/GSL pools start at the WT steady state (fetal
    circulation and maternal handling keep the newborn near normal) and are
    then allowed to fill in at the genotype-determined rate.  The NPC1 protein
    pools start at their genotype steady state.
    """
    g = GENOTYPES[genotype]
    y = [0.0] * NS
    theta = g["theta0"]
    er = p.ksyn_npc1 * (1 - g["f_null"]) / (p.kexit0 * theta + p.kerad)
    y[I_NPC1_ER] = er
    y[I_NPC1_L] = er * p.kexit0 * theta / p.kdegL
    y[I_CHOL_V] = WT["CHOL_V"]
    y[I_CHOL_C] = WT["CHOL_C"]
    y[I_GSL_V] = WT["GSL_V"]
    y[I_GSL_C] = WT["GSL_C"]
    y[I_SPH] = WT["SPH"]
    y[I_CA_LY] = WT["CA"]
    y[I_HYD] = WT["HYD"]
    y[I_AUTOPH] = WT["AUT"]
    y[I_MITO] = WT["MITO"]
    y[I_ROS] = WT["ROS"]
    y[I_PC] = 1.0
    y[I_SYN] = 1.0
    y[I_CBL] = 1.0
    y[I_OHC] = 1.0
    y[I_TRIOL] = p.ktri * WT["CHOL_V"] / (p.Ktri + WT["CHOL_V"]) / p.ktri_el
    y[I_PPCS] = (p.kppcs * (WT["CHOL_V"] + p.b_ppcs_gsl * WT["GSL_V"]) /
                     (p.Kppcs + WT["CHOL_V"] + p.b_ppcs_gsl * WT["GSL_V"]) /
                     p.kppcs_el)
    y[I_TCG] = p.ktcg * WT["CHOL_V"] / (p.Ktcg + WT["CHOL_V"]) / p.ktcg_el
    y[I_NFL] = p.knfl_base / p.knfl_el
    return y


def simulate(genotype="I1061T/I1061T", years=30.0, regimens=(), dt=0.01,
             record_every=1.0, nal_sym_on=1.0, nal_dm_on=1.0,
             liver_risk=0.0, y0=None, p=PAR, t0=0.0):
    """Fixed-step RK4 with dosing applied as instantaneous bolus at grid times."""
    ctx = dict(pr=p, geno=GENOTYPES[genotype], nal_sym_on=nal_sym_on,
               nal_dm_on=nal_dm_on, liver_risk=liver_risk)
    y = initial_state(genotype, p) if y0 is None else [float(v) for v in y0]
    tend = t0 + years * 365.25
    nsteps = int(round((tend - t0) / dt))
    rec_stride = max(1, int(round(record_every / dt)))

    ts, ys = [t0], [list(y)]
    t = t0
    h6 = dt / 6.0
    hh = dt / 2.0
    for k in range(nsteps):
        # doses at the START of this step
        for r in regimens:
            amt = r.doses_in(t - dt / 2, t + dt / 2)
            if amt:
                y[IX[r.target]] += amt
        k1 = rhs(t, y, ctx)
        k2 = rhs(t + hh, [a + hh * b for a, b in zip(y, k1)], ctx)
        k3 = rhs(t + hh, [a + hh * b for a, b in zip(y, k2)], ctx)
        k4 = rhs(t + dt, [a + dt * b for a, b in zip(y, k3)], ctx)
        # RK4 combination, with hard non-negativity on non-negative pools
        y = [v if v > 0.0 else 0.0 for v in
             (u + h6 * (a + 2.0 * b + 2.0 * c + e)
              for u, a, b, c, e in zip(y, k1, k2, k3, k4))]
        t += dt
        if (k + 1) % rec_stride == 0:
            ts.append(t)
            ys.append(list(y))
    # ALWAYS record the final state.  Without this, ys[-1] is the last state on
    # the recording grid, which can be a whole `record_every` earlier than the
    # requested end -- and every burn-in that takes ys[-1] as a trial-entry state
    # is then silently shifted back in time.  (This was a real defect: with
    # record_every = 400 d it shifted trial entry from age 13.0 to age 12.05 and
    # moved every treatment-scenario number in the report.)
    if not ts or ts[-1] < t - 1e-9:
        ts.append(t)
        ys.append(list(y))
    return np.array(ts), np.array(ys)


def out_table(ts, ys, p=PAR):
    """Attach derived outputs to a trajectory."""
    o = {"t": ts, "years": ts / 365.25}
    for s in STATES:
        o[s] = ys[:, IX[s]]
    for name, fn in DERIVED.items():
        o[name] = np.array([fn(row, p) for row in ys])
    o["f_NPC1"] = ys[:, I_NPC1_L] / WT["NPC1_L_WT"]
    o["CHOL_C_fold"] = ys[:, I_CHOL_C] / WT["CHOL_C"]
    o["CHOL_V_fold"] = ys[:, I_CHOL_V] / WT["CHOL_V"]
    o["FUNC"] = np.maximum(0.0, 1.0 - ys[:, I_D_REV] - ys[:, I_D_IRR])
    return o


def at(o, key, year):
    return float(np.interp(year, o["years"], o[key]))


# =============================================================================
# 7.  CALIBRATION HELPERS
# =============================================================================

def secant_solve(f, x0, x1, target, tol=1e-4, maxit=25):
    """Simple secant root find for f(x) = target."""
    f0, f1 = f(x0) - target, f(x1) - target
    for _ in range(maxit):
        if abs(f1) < tol:
            return x1
        if abs(f1 - f0) < 1e-14:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = max(1e-6, x2)
        x0, f0, x1 = x1, f1, x2
        f1 = f(x1) - target
    return x1


# =============================================================================
# 8.  SCENARIOS
# =============================================================================

def natural_history(genotype, years=35.0, dt=0.02, liver_risk=0.0):
    ts, ys = simulate(genotype, years=years, dt=dt, record_every=7.0,
                      liver_risk=liver_risk)
    return out_table(ts, ys)


def onset_age(o, thresh=3.0):
    """Age (years) at which NPCCSS5 first exceeds `thresh` -- clinical onset."""
    idx = np.where(o["NPCCSS5"] >= thresh)[0]
    return float(o["years"][idx[0]]) if len(idx) else float("nan")


def median_survival(o):
    idx = np.where(o["CUMHAZ"] >= math.log(2.0))[0]
    return float(o["years"][idx[0]]) if len(idx) else float("nan")


def burn_in(genotype, years, dt=0.02):
    """Run untreated to `years` of age and return the state, for trial entry."""
    ts, ys = simulate(genotype, years=years, dt=dt, record_every=30.0)
    return ys[-1].copy()


# =============================================================================
# 9.  MAIN REPORT
# =============================================================================

def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="coarser dt")
    ap.add_argument("--emit-param", action="store_true",
                    help="print the mrgsolve $PARAM block and exit")
    args = ap.parse_args()
    if args.emit_param:
        for k, v in P.items():
            print(f"{k} = {v},")
        return

    DTS = 0.30 if args.quick else 0.15      # drug-free (slow) integration step
    DTD = 0.06 if args.quick else 0.03      # with drugs on board

    def run(geno, years, regs=(), y0=None, t0=0.0, vd=None, pp=None,
            rec=14.0, dt=None, **kw):
        pp = (pp or PAR)
        if vd is not None:
            pp = pp.copy_with(v_dev=vd)
        dt = dt or (DTD if regs else DTS)
        ts, ys = simulate(geno, years=years, regimens=list(regs), dt=dt,
                          record_every=rec, y0=y0, t0=t0, p=pp, **kw)
        return out_table(ts, ys, pp)

    def entry_state(geno, age, vd):
        return simulate(geno, years=age, dt=DTS, record_every=400.0,
                        p=PAR.copy_with(v_dev=vd))[1][-1]

    # ---------------------------------------------------------------- 1
    hr("1. WILD-TYPE REFERENCE STATE (every normalisation anchors here)")
    for k in ["CHOL_V", "CHOL_C", "SPH", "CA", "HYD", "GSL_C", "AUT", "MITO", "ROS"]:
        print(f"   WT {k:8s} = {WT[k]:.4f}")
    print("   A healthy system therefore sits at cerebellar stress = 0 exactly.")

    # ---------------------------------------------------------------- 2
    hr("2. ONSET ARCHETYPES = GENOTYPE x DEVELOPMENTAL VULNERABILITY")
    print("   Residual NPC1 activity alone does NOT fix the onset form: the storage")
    print("   phenotype SATURATES, so I1061T/null and I1061T/I1061T differ 2-fold in")
    print("   residual activity and hardly at all in steady-state load.  v_dev is")
    print("   carried as a PATIENT covariate (PMID 15217094), not derived.\n")
    print(f"   {'archetype':18s} {'genotype':16s} {'v_dev':>6s} {'f_NPC1':>8s} "
          f"{'CHOL_C x':>9s} {'onset':>6s} {'surv':>6s}")
    ARCH = [("perinatal", "null/null", 2.5, 1.0, 6.0),
            ("early-infantile", "I1061T/null", 2.5, 0.0, 25.0),
            ("late-infantile", "I1061T/null", 1.2, 0.0, 25.0),
            ("juvenile", "I1061T/I1061T", 0.6, 0.0, 40.0),
            ("adolescent/adult", "mild/mild", 0.2, 0.0, 60.0),
            ("NPC2 disease", "NPC2", 1.0, 0.0, 30.0)]
    arch_out = {}
    for name, g, vd, lr, yrs in ARCH:
        o = run(g, yrs, vd=vd, liver_risk=lr)
        arch_out[name] = o
        print(f"   {name:18s} {g:16s} {vd:6.1f} {at(o,'f_NPC1',min(5,yrs-1)):8.4f} "
              f"{at(o,'CHOL_C_fold',min(10,yrs-1)):9.1f} "
              f"{onset_age(o):6.1f} {median_survival(o):6.1f}")
    o_j = arch_out["juvenile"]
    o_wt = run("WT", 45.0, vd=0.6)

    # ---------------------------------------------------------------- 3
    hr("3. BIOMARKER CALIBRATION (T1, T2) AND WHAT IT PREDICTS")
    tri_pat, tri_ctl = at(o_j, "TRIOL", 10.0), at(o_wt, "TRIOL", 10.0)
    print(f"   plasma triol, NPC       {tri_pat:8.2f} ng/mL   target 88.31")
    print(f"   plasma triol, control   {tri_ctl:8.2f} ng/mL   target  5.97")
    print(f"   fold-change             {tri_pat/tri_ctl:8.2f} x       target 14.79")
    cv = at(o_j, "CHOL_V", 10.0)
    print(f"\n   The two-point fit forces a SATURATION constant Ktri = {PAR.Ktri:.2f} while the")
    print(f"   patient's visceral pool sits at CHOL_V = {cv:.1f}, i.e. triol generation is at")
    print(f"   {100*cv/(PAR.Ktri+cv):.1f}% of its own ceiling.  The weak severity correlation is")
    print( "   therefore PREDICTED, not fitted -- see section 8.")

    # ---------------------------------------------------------------- 4
    hr("4. NATURAL HISTORY OF THE JUVENILE PHENOTYPE (T3; T4 is validation)")
    print(f"   {'age':>5s} {'NPCCSS5':>8s} {'NPCCSS17':>9s} {'SARA':>6s} {'D_rev':>7s} "
          f"{'D_irr':>7s} {'PC_lost':>8s} {'DAM/res':>8s} {'sacc%':>6s} {'CBL':>6s}")
    for yr in [2, 4, 6, 8, 10, 12, 13, 16, 20, 25, 30]:
        print(f"   {yr:5d} {at(o_j,'NPCCSS5',yr):8.2f} {at(o_j,'NPCCSS17',yr):9.2f} "
              f"{at(o_j,'SARA',yr):6.2f} {at(o_j,'D_REV',yr):7.3f} {at(o_j,'D_IRR',yr):7.3f} "
              f"{at(o_j,'PC_LOST',yr):8.3f} {at(o_j,'DAM',yr)/PAR.D_reserve:8.2f} "
              f"{at(o_j,'SACCADE',yr):6.1f} {at(o_j,'CBL',yr):6.3f}")
    sl5 = at(o_j, "NPCCSS5", 13) - at(o_j, "NPCCSS5", 12)
    sl17 = at(o_j, "NPCCSS17", 13) - at(o_j, "NPCCSS17", 12)
    print(f"\n   5-domain slope, age 12-13    {sl5:6.3f} pt/yr   target 1.50   (T3, fitted)")
    print(f"   17-domain slope, age 12-13   {sl17:6.3f} pt/yr   target 2.7-2.9 (T4, VALIDATION)")
    print(f"   SARA slope, age 13-14        "
          f"{at(o_j,'SARA',14)-at(o_j,'SARA',13):6.3f} pt/yr   (derived, no published value)")
    print(f"   median survival              {median_survival(o_j):6.1f} yr")
    print( "\n   ONSET-INDEPENDENT SLOPE (T14, structural validation of the damage integral):")
    for name in ["late-infantile", "juvenile", "adolescent/adult"]:
        o = arch_out[name]
        a0 = onset_age(o)
        if a0 == a0:
            print(f"      {name:18s} onset {a0:5.1f} yr -> slope over the next 5 yr "
                  f"{(at(o,'NPCCSS5',a0+5)-at(o,'NPCCSS5',a0))/5:5.2f} pt/yr")

    # ---------------------------------------------------------------- 5
    hr("5. CONSERVATION AND MONOTONICITY CHECKS")
    tot = o_j["PC"] + o_j["PC_S"] + o_j["PC_LOST"]
    print(f"   Purkinje pool conservation  max|PC+PC_S+PC_LOST-1| = "
          f"{float(np.max(np.abs(tot-1.0))):.3e}")
    print(f"   D_IRR monotone              min step = {float(np.min(np.diff(o_j['D_IRR']))):.3e}")
    print(f"   PC_LOST monotone            min step = {float(np.min(np.diff(o_j['PC_LOST']))):.3e}")
    print(f"   DAM monotone                min step = {float(np.min(np.diff(o_j['DAM']))):.3e}")
    print(f"   WT healthy at 45 yr         NPCCSS5 = {at(o_wt,'NPCCSS5',45):.4f}, "
          f"PC_LOST = {at(o_wt,'PC_LOST',45):.2e}, DAM/res = "
          f"{at(o_wt,'DAM',45)/PAR.D_reserve:.3f}")

    # ---------------------------------------------------------------- 6
    hr("6. DRUG EXPOSURE (steady state, adult dosing)")
    y13 = entry_state("I1061T/I1061T", 13.0, 0.6)
    for label, reg, kind, lit in [
        ("miglustat 200 mg tid", [mig_regimen()], "mig", "Css 11-14 uM"),
        ("arimoclomol 124 mg tid", [ari_regimen()], "ari", "t1/2 ~4 h, linear"),
        ("levacetylleucine 4 g/d", [nal_regimen()], "nal", "label Cmax 8.3 ug/mL"),
        ("IT adrabetadex 900 mg", [cd_regimen()], "cd", "CSF t1/2 ~4 h")]:
        o = run("I1061T/I1061T", 0.25, reg, y0=y13, t0=13*365.25, vd=0.6,
                dt=0.002, rec=0.01)
        n = len(o["t"]) // 4
        if kind == "mig":
            c = o["MIG_CEN"][-n:] / PAR.mig_V / PAR.mig_MW * 1000.0
            extra = (f"brain {float(np.mean(o['MIG_BR'][-n:])):.2f} uM -> "
                     f"CNS UGCG inhib "
                     f"{100*PAR.mig_Emax*float(np.mean(o['MIG_BR'][-n:]))/(PAR.mig_IC50+float(np.mean(o['MIG_BR'][-n:]))):.1f}% "
                     f"vs visceral "
                     f"{100*PAR.mig_Emax*float(np.mean(c))/(PAR.mig_IC50+float(np.mean(c))):.1f}%")
        elif kind == "ari":
            c = o["ARI_CEN"][-n:] / PAR.ari_V * 1000.0
            extra = (f"CSF {float(np.mean(o['ARI_CSF'][-n:])):.1f} ng/mL -> HSP70 "
                     f"{float(np.mean(o['HSP70'][-n:]) if 'HSP70' in o else 0):.2f}x"
                     if 'HSP70' in o else
                     f"CSF {float(np.mean(o['ARI_CSF'][-n:])):.1f} ng/mL")
        elif kind == "nal":
            c = o["NAL_CEN"][-n:] / PAR.nal_V
            extra = f"effect site {float(np.mean(o['NAL_BR'][-n:])):.3f} mg/L"
        else:
            c = o["CD_CSF"][-n:] / PAR.cd_Vcsf
            extra = f"peak CSF {float(np.max(c)):.0f} mg/L"
        print(f"   {label:24s} Cavg {float(np.mean(c)):8.2f}  Cmax {float(np.max(c)):8.2f}"
              f"   [{lit}]")
        print(f"   {'':24s} {extra}")

    # ---------------------------------------------------------------- 7
    hr("7. TREATMENT SCENARIOS  (juvenile, entry at age 13 where SARA = 15.91)")
    SCEN = OrderedDict([
        ("S01 untreated", []),
        ("S06 miglustat", [mig_regimen()]),
        ("S07 levacetylleucine", [nal_regimen()]),
        ("S08 arimoclomol", [ari_regimen()]),
        ("S09 arimoclomol + miglustat", [ari_regimen(), mig_regimen()]),
        ("S10 levacetylleucine + miglustat", [nal_regimen(), mig_regimen()]),
        ("S11 triple oral", [nal_regimen(), mig_regimen(), ari_regimen()]),
        ("S12 IT cyclodextrin 900 q2wk", [cd_regimen()]),
        ("S13 IT cyclodextrin + miglustat", [cd_regimen(), mig_regimen()]),
    ])
    print(f"   {'scenario':34s} {'dN5/yr':>7s} {'dN4/yr':>7s} {'dSARA 12wk':>11s} "
          f"{'D_rev':>7s} {'D_irr':>7s} {'triol':>7s} {'dB':>5s}")
    res = {}
    for name, reg in SCEN.items():
        o = run("I1061T/I1061T", 3.0, reg, y0=y13, t0=13*365.25, vd=0.6, rec=3.5)
        res[name] = o
        wk12 = 13.0 + 12*7/365.25
        print(f"   {name:34s} {at(o,'NPCCSS5',14)-at(o,'NPCCSS5',13):7.3f} "
              f"{at(o,'NPCCSS4',14)-at(o,'NPCCSS4',13):7.3f} "
              f"{at(o,'SARA',wk12)-at(o,'SARA',13):11.3f} "
              f"{at(o,'D_REV',14):7.3f} {at(o,'D_IRR',14):7.3f} "
              f"{at(o,'TRIOL',14):7.1f} {at(o,'HEARING_dB',16):5.1f}")

    # ---------------------------------------------------------------- 8
    hr("8. THE COMPARTMENT CLAIM: triol vs severity across a virtual cohort (T9)")
    # The published rho comes from 36 patients aged 2-18 (PMID 33228797), i.e. a
    # PAEDIATRIC cohort with essentially no adult-onset genotypes.  The virtual
    # cohort is matched to that enrolment; including adult-onset patients would
    # widen the triol range and inflate the correlation, which is a different
    # comparison and is reported separately below.
    def cohort(names, ages):
        rows = []
        for nm in names:
            o = arch_out[nm]
            for age in ages:
                if age < o["years"][-1]:
                    n5 = at(o, "NPCCSS5", age)
                    if n5 > 0.2:
                        rows.append((at(o, "TRIOL", age), n5))
        return rows
    AGES = [2, 4, 6, 8, 10, 12, 14, 16, 18]
    rows = cohort(["early-infantile", "late-infantile", "juvenile"], AGES)
    rows_wide = cohort(["early-infantile", "late-infantile", "juvenile",
                        "adolescent/adult"], AGES + [22, 26, 30, 35])
    tri = np.array([r[0] for r in rows]); n5v = np.array([r[1] for r in rows])
    triW = np.array([r[0] for r in rows_wide]); n5W = np.array([r[1] for r in rows_wide])

    def spearman(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        ra -= ra.mean(); rb -= rb.mean()
        return float(ra @ rb / math.sqrt((ra @ ra) * (rb @ rb)))
    rho = spearman(tri, n5v)
    rhoW = spearman(triW, n5W)
    print(f"   study-matched cohort (ages 2-18, paediatric genotypes) n = {len(rows)}")
    print(f"   Spearman rho(triol, NPCCSS5) = {rho:.3f}      published rho = 0.265 (T9)")
    print(f"   triol spans {tri.min():.1f}-{tri.max():.1f} ng/mL while NPCCSS5 spans "
          f"{n5v.min():.1f}-{n5v.max():.1f}")
    print(f"\n   cohort WIDENED to include adult-onset disease   n = {len(rows_wide)}")
    print(f"   Spearman rho = {rhoW:.3f}, triol spans {triW.min():.1f}-{triW.max():.1f} ng/mL")
    print( "   No noise term is used anywhere.  In the study-matched cohort the triol")
    print( "   range is narrow BECAUSE generation is saturated, so ranking patients by")
    print( "   triol barely ranks them by severity.  Widening the cohort to adult-onset")
    print( "   disease -- where the visceral pool is far below saturation -- restores the")
    print( "   correlation.  The biomarker is informative exactly where it is not")
    print( "   saturated, which is not where the trials enrol.")

    # ---------------------------------------------------------------- 9
    hr("9. THE DESIGN CLAIM: each published design against BOTH mechanisms")
    def bench(sym, dm, dm_scale, years, key, regs=None):
        pp = PAR.copy_with(nal_Emax_dm=PAR.nal_Emax_dm * dm_scale, v_dev=0.6)
        regs = [nal_regimen()] if regs is None else regs
        o = run("I1061T/I1061T", years, regs, y0=y13, t0=13*365.25, pp=pp,
                rec=3.5, nal_sym_on=sym, nal_dm_on=dm)
        o0 = run("I1061T/I1061T", years, [], y0=y13, t0=13*365.25, pp=pp, rec=3.5)
        te = 13.0 + years
        return (at(o, key, te) - at(o, key, 13.0)), (at(o0, key, te) - at(o0, key, 13.0))

    sym12, pbo12 = bench(1.0, 0.0, 1.0, 1.0, "NPCCSS5")
    sym_benefit = sym12 - pbo12
    print(f"   symptomatic-only drug, 12-month NPCCSS5 difference {sym_benefit:+.3f} pt")
    scale = secant_solve(lambda s: bench(0.0, 1.0, s, 1.0, "NPCCSS5")[0]
                         - bench(0.0, 1.0, s, 1.0, "NPCCSS5")[1],
                         1.0, 5.0, sym_benefit, tol=3e-3, maxit=10)
    dm12, _ = bench(0.0, 1.0, scale, 1.0, "NPCCSS5")
    print(f"   disease-modifying drug tuned to match (dm potency x{scale:.2f}) "
          f"{dm12-pbo12:+.3f} pt")
    wk = 12*7/365.25
    s_sym, s_pbo = bench(1.0, 0.0, 1.0, wk, "SARA")
    s_dm, _ = bench(0.0, 1.0, scale, wk, "SARA")
    print(f"\n   --- DESIGN A: 12-week crossover on SARA (IB1001-301; T8 = -1.28) ---")
    print(f"   placebo arm            {s_pbo:+.3f} pt")
    print(f"   symptomatic drug       {s_sym:+.3f} pt  -> difference {s_sym-s_pbo:+.3f}")
    print(f"   disease-modifying drug {s_dm:+.3f} pt  -> difference {s_dm-s_pbo:+.3f}")
    r = abs(s_sym-s_pbo)/max(1e-9, abs(s_dm-s_pbo))
    print(f"   SENSITIVITY RATIO      {r:.1f} x in favour of the symptomatic mechanism")
    print(f"\n   --- DESIGN B: 12-month parallel on NPCCSS5 (arimoclomol; T5 = -1.40) ---")
    print(f"   symptomatic drug       {sym_benefit:+.3f} pt")
    print(f"   disease-modifying drug {dm12-pbo12:+.3f} pt  (matched by construction)")
    print( "   Design B cannot tell them apart; design A can only see the offset.")
    print( "   NEITHER published design identifies the split.  Only withdrawal does.")

    # ---------------------------------------------------------------- 10
    hr("10. WITHDRAWAL: the only design that separates the two components")
    o = run("I1061T/I1061T", 2.0, [nal_regimen(stop=14*365.25)], y0=y13,
            t0=13*365.25, vd=0.6, rec=3.5)
    o0 = run("I1061T/I1061T", 2.0, [], y0=y13, t0=13*365.25, vd=0.6, rec=3.5)
    print(f"   {'age':>6s} {'SARA':>7s} {'untr':>7s} {'D_rev':>7s} {'D_irr':>7s}  phase")
    for yr, ph in [(13.0, "baseline"), (13.25, "on drug 3 mo"), (14.0, "on drug 12 mo"),
                   (14.10, "36 d off"), (14.25, "90 d off"), (15.0, "12 mo off")]:
        print(f"   {yr:6.2f} {at(o,'SARA',yr):7.3f} {at(o0,'SARA',yr):7.3f} "
              f"{at(o,'D_REV',yr):7.3f} {at(o,'D_IRR',yr):7.3f}  {ph}")
    print(f"\n   SARA rebound in the 90 d after stopping {at(o,'SARA',14.25)-at(o,'SARA',14.0):+.3f} pt")
    print(f"   residual gap vs untreated at 12 mo off    "
          f"{at(o,'SARA',15.0)-at(o0,'SARA',15.0):+.3f} pt  (the disease-modifying part)")

    # ---------------------------------------------------------------- 11
    hr("11. LEVACETYLLEUCINE vs THE PUBLISHED LONG-TERM EXTENSION (T12, T13)")
    print("   The model reproduces the 12-week randomised difference exactly, and then")
    print("   FAILS to reproduce the single-arm extension.  Reported, not tuned away.\n")
    print(f"   {'window':8s} {'drug (from baseline)':>21s} {'untreated':>10s} "
          f"{'difference':>11s}   published")
    for yrs, lab, pub in [(wk, "12 wk", "-1.28 vs placebo"),
                          (1.0, "12 mo", "-1.88 from baseline"),
                          (1.5, "18 mo", "-1.64 from baseline")]:
        o = run("I1061T/I1061T", yrs, [nal_regimen()], y0=y13, t0=13*365.25,
                vd=0.6, rec=3.5)
        o0 = run("I1061T/I1061T", yrs, [], y0=y13, t0=13*365.25, vd=0.6, rec=3.5)
        te = 13.0 + yrs
        a = at(o, "SARA", te) - at(o, "SARA", 13.0)
        b = at(o0, "SARA", te) - at(o0, "SARA", 13.0)
        print(f"   {lab:8s} {a:21.3f} {b:10.3f} {a-b:11.3f}   {pub}")
    print("\n   Implied untreated SARA slope if the drug is purely symptomatic with the")
    print("   trial's own within-arm offset of -1.97 pt:")
    for t, obs in [(1.0, 1.88), (1.5, 1.64)]:
        print(f"      at {t:.1f} yr, observed -{obs:.2f} from baseline -> implied slope "
              f"{(1.97-obs)/t:+.3f} pt/yr")
    print(f"   The model's own derived untreated SARA slope is "
          f"{at(o_j,'SARA',14)-at(o_j,'SARA',13):+.3f} pt/yr -- a 7-16x disagreement.")
    print("   There is no published untreated SARA slope for NPC, so this is left open.")

    # ---------------------------------------------------------------- 12
    hr("12. ARIMOCLOMOL: the effect size is gated by the RESERVE, not by the dose")
    print(f"   {'Emax_fold':>10s} {'f_NPC1 off':>11s} {'f_NPC1 on':>10s} {'fold':>6s} "
          f"{'dN5 untr':>9s} {'dN5 ari':>8s} {'diff':>7s} {'reduction':>10s}")
    for ef in [0.15, 0.30, 0.60, 1.20, 3.40, 60.0]:
        pp = PAR.copy_with(Emax_fold=ef, v_dev=0.6)
        o = run("I1061T/I1061T", 1.0, [ari_regimen()], y0=y13, t0=13*365.25, pp=pp, rec=7.0)
        o0 = run("I1061T/I1061T", 1.0, [], y0=y13, t0=13*365.25, pp=pp, rec=7.0)
        a = at(o, "NPCCSS5", 14) - at(o, "NPCCSS5", 13)
        b = at(o0, "NPCCSS5", 14) - at(o0, "NPCCSS5", 13)
        f0, f1 = at(o0, "f_NPC1", 14), at(o, "f_NPC1", 14)
        print(f"   {ef:10.2f} {f0:11.4f} {f1:10.4f} {f1/f0:6.2f} {b:9.3f} {a:8.3f} "
              f"{a-b:7.3f} {100*(1-a/b):9.1f}%")
    print(f"\n   Published (PMID 34418116): -1.40 pt, 65% of progression, "
          f"95% CI -2.76 to -0.03.")
    print( "   Even at the arithmetic CEILING of this mechanism (theta = 1, i.e. folding")
    print( "   fully corrected) the model removes only ~48% of progression in a patient")
    print( "   who entered at age 13.  Now vary ONLY the age at entry:\n")
    print(f"   {'entry age':>10s} {'DAM/reserve':>12s} {'dN5 untr':>9s} {'dN5 ari':>8s} "
          f"{'diff':>7s} {'reduction':>10s}")
    for age in [5.0, 8.0, 10.0, 13.0, 16.0]:
        ya = entry_state("I1061T/I1061T", age, 0.6)
        o = run("I1061T/I1061T", 1.0, [ari_regimen()], y0=ya, t0=age*365.25, vd=0.6, rec=7.0)
        o0 = run("I1061T/I1061T", 1.0, [], y0=ya, t0=age*365.25, vd=0.6, rec=7.0)
        a = at(o, "NPCCSS5", age+1) - at(o, "NPCCSS5", age)
        b = at(o0, "NPCCSS5", age+1) - at(o0, "NPCCSS5", age)
        print(f"   {age:10.1f} {at(o0,'DAM',age)/PAR.D_reserve:12.2f} {b:9.3f} {a:8.3f} "
              f"{a-b:7.3f} {100*(1-a/b) if abs(b)>1e-6 else float('nan'):9.1f}%")
    print( "\n   FALSIFIABLE PREDICTION: arimoclomol's effect should be strongly modified")
    print( "   by reserve status at entry -- large before the gate opens, near-absent")
    print( "   after.  The trial enrolled ages 2-18, so its point estimate is dominated")
    print( "   by its youngest patients, and its wide CI is what that looks like.")

    # ---------------------------------------------------------------- 13
    hr("13. GENOTYPE x MECHANISM: HSP70 needs rescuable NPC1; cyclodextrin does not")
    print(f"   {'genotype':16s} {'f_NPC1 off':>11s} {'f_NPC1 on':>10s} "
          f"{'ari benefit':>12s} {'CD benefit':>11s}")
    for gname, age in [("I1061T/I1061T", 8.0), ("I1061T/null", 6.0),
                       ("mild/mild", 27.0), ("null/null", 0.8), ("NPC2", 8.0)]:
        vd = 2.5 if gname == "null/null" else 0.6
        ya = entry_state(gname, age, vd)
        o0 = run(gname, 1.0, [], y0=ya, t0=age*365.25, vd=vd, rec=7.0)
        b = at(o0, "NPCCSS5", age+1) - at(o0, "NPCCSS5", age)
        oa = run(gname, 1.0, [ari_regimen()], y0=ya, t0=age*365.25, vd=vd, rec=7.0)
        oc = run(gname, 1.0, [cd_regimen()], y0=ya, t0=age*365.25, vd=vd, rec=7.0)
        da = (at(oa, "NPCCSS5", age+1) - at(oa, "NPCCSS5", age)) - b
        dc = (at(oc, "NPCCSS5", age+1) - at(oc, "NPCCSS5", age)) - b
        print(f"   {gname:16s} {at(o0,'f_NPC1',age+1):11.4f} {at(oa,'f_NPC1',age+1):10.4f} "
              f"{da:12.3f} {dc:11.3f}")
    print("\n   Arimoclomol is inert in null/null and in NPC2 disease by construction of")
    print("   the mechanism, not by a switch.  Cyclodextrin keeps its effect in both.")

    # ---------------------------------------------------------------- 14
    hr("14. MIGLUSTAT x ARIMOCLOMOL: multiplicative on one flux (T15)")
    eff = {}
    for name, reg in [("neither", []), ("miglustat", [mig_regimen()]),
                      ("arimoclomol", [ari_regimen()]),
                      ("both", [ari_regimen(), mig_regimen()])]:
        o = run("I1061T/I1061T", 1.0, reg, y0=y13, t0=13*365.25, vd=0.6, rec=7.0)
        eff[name] = at(o, "NPCCSS5", 14) - at(o, "NPCCSS5", 13)
        print(f"   {name:12s} 12-month NPCCSS5 progression {eff[name]:6.3f} pt")
    a = eff["neither"] - eff["arimoclomol"]
    mg = eff["neither"] - eff["miglustat"]
    bo = eff["neither"] - eff["both"]
    print(f"\n   arimoclomol alone saves {a:6.3f};  miglustat alone {mg:6.3f};  "
          f"additive {a+mg:6.3f}")
    print(f"   model, both together    {bo:6.3f}  -> "
          f"{'SUPRA' if bo > a+mg else 'sub'}-additive by {bo-(a+mg):+.3f} pt")
    print(f"   arimoclomol saves {bo-mg:6.3f} pt on a miglustat background vs "
          f"{a:6.3f} pt alone")
    print( "\n   THE PREDICTED SYNERGY DOES NOT APPEAR.  It was expected, because HSP70")
    print( "   raises the folding YIELD while miglustat lowers the LOAD on the same flux,")
    print( "   so the two should multiply.  They do not, for a reason the model makes")
    print( "   plain: miglustat's own CNS effect is only 12.5% UGCG inhibition, worth")
    print(f"   {mg:.3f} pt/yr, so there is almost nothing for arimoclomol to interact with.")
    print( "   Published: -1.40 overall but -2.06 in the miglustat stratum (PMID 34418116).")
    print( "   The model does NOT reproduce that stratum difference and says so.")

    # ---------------------------------------------------------------- 15
    hr("15. ADRABETADEX: mechanism-inseparable ototoxicity eats its own benefit (T10)")
    ocd = run("I1061T/I1061T", 1.5, [cd_regimen()], y0=y13, t0=13*365.25, vd=0.6, rec=7.0)
    o0 = run("I1061T/I1061T", 1.5, [], y0=y13, t0=13*365.25, vd=0.6, rec=7.0)
    print(f"   {'metric':38s} {'untreated':>10s} {'IT-CD':>9s}")
    for lbl, key in [("CNS cholesterol (x WT) at 12 mo", "CHOL_C_fold"),
                     ("CSF calbindin at 12 mo", "CALB"),
                     ("CSF NfL at 12 mo", "NFL"),
                     ("plasma triol at 12 mo", "TRIOL"),
                     ("hearing threshold shift, dB, 12 mo", "HEARING_dB")]:
        print(f"   {lbl:38s} {at(o0,key,14):10.2f} {at(ocd,key,14):9.2f}")
    d5 = (at(ocd,"NPCCSS5",14)-at(ocd,"NPCCSS5",13)) - (at(o0,"NPCCSS5",14)-at(o0,"NPCCSS5",13))
    d17 = (at(ocd,"NPCCSS17",14)-at(ocd,"NPCCSS17",13)) - (at(o0,"NPCCSS17",14)-at(o0,"NPCCSS17",13))
    nh = lambda oo, y: at(oo, "NPCCSS5", y) * PAR.n17_k
    d17n = (nh(ocd,14)-nh(ocd,13)) - (nh(o0,14)-nh(o0,13))
    print(f"\n   5-domain benefit (no hearing domain)      {d5:+.3f} pt")
    print(f"   17-domain benefit WITHOUT hearing scored  {d17n:+.3f} pt")
    print(f"   17-domain benefit WITH hearing scored     {d17:+.3f} pt")
    if abs(d17n) > 1e-9:
        if d17 * d17n < 0:
            print(f"   the hearing domain does not merely erase the benefit, it REVERSES")
            print(f"   the sign: {d17n:+.2f} pt of benefit becomes {d17:+.2f} pt of harm")
        else:
            print(f"   the drug's own ototoxicity consumes {100*(1-d17/d17n):.0f}% of its "
                  f"measured benefit")
    print( "   This is the model-level account of Ory 2017 scoring 'NSS minus hearing'")
    print( "   (PMID 28803710), and of the randomised sham-controlled phase 2b/3 -- with")
    print( "   hearing in the composite and no open-label expectancy -- finding nothing.")

    # ---------------------------------------------------------------- 16
    hr("16. EARLY vs LATE START: the reserve is an integral, so waiting costs")
    ref = run("I1061T/I1061T", 20.0, [], vd=0.6, rec=14.0)
    ref20 = at(ref, "NPCCSS5", 20.0)
    print(f"   {'start age':>10s} {'DAM/res at start':>17s} {'yrs on drug':>12s} "
          f"{'NPCCSS5@20':>11s} {'D_rev@20':>9s} {'D_irr@20':>9s} {'pts saved':>10s}")
    for start in [2.0, 5.0, 8.0, 12.0]:
        ya = entry_state("I1061T/I1061T", start, 0.6)
        regs = [nal_regimen(start=start*365.25), mig_regimen(start=start*365.25),
                ari_regimen(start=start*365.25)]
        o = run("I1061T/I1061T", 20.0-start, regs, y0=ya, t0=start*365.25, vd=0.6, rec=14.0)
        print(f"   {start:10.1f} {at(ref,'DAM',start)/PAR.D_reserve:17.2f} "
              f"{20.0-start:12.1f} {at(o,'NPCCSS5',20):11.2f} {at(o,'D_REV',20):9.3f} "
              f"{at(o,'D_IRR',20):9.3f} {ref20-at(o,'NPCCSS5',20):10.2f}")
    print(f"   {'untreated':>10s} {'-':>17s} {0.0:12.1f} {ref20:11.2f} "
          f"{at(ref,'D_REV',20):9.3f} {at(ref,'D_IRR',20):9.3f} {0.0:10.2f}")

    # ---------------------------------------------------------------- 17
    hr("17. MODEL vs PUBLISHED TARGETS")
    rows = [("T1  plasma triol, patients (ng/mL)", 88.31, tri_pat, "fit"),
            ("T2  plasma triol, controls (ng/mL)", 5.97, tri_ctl, "fit"),
            ("T3  5-domain NPCCSS slope (pt/yr)", 1.50, sl5, "fit"),
            ("T4  17-domain NPCCSS slope (pt/yr)", 2.80, sl17, "VALID"),
            ("T8  levacetylleucine 12-wk SARA diff", -1.28, s_sym - s_pbo, "fit"),
            ("T9  Spearman rho triol vs NPCCSS5", 0.265, rho, "VALID"),
            ("T11 SARA at trial entry", 15.91, at(o_j, "SARA", 13.0), "fit")]
    print(f"   {'target':38s} {'role':>6s} {'published':>10s} {'model':>9s} {'err %':>8s}")
    for lbl, pub, mod, role in rows:
        print(f"   {lbl:38s} {role:>6s} {pub:10.3f} {mod:9.3f} "
              f"{100*(mod-pub)/pub:7.1f}%")
    print("\n   The MISSES are T5/T6 (section 12), T12/T13 (section 11) and T15")
    print("   (section 14).  None of them was tuned away; see README.md section 7.")
    print("\nDone.\n")


if __name__ == "__main__":
    main()
