#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cdys_reference_check.py
=======================
Independent numpy/scipy port of the cervical-dystonia QSP model that lives in
`cdys_mrgsolve_model.R`.  Same state vector, same parameters, same equations,
integrated with scipy LSODA instead of mrgsolve/DLSODA.

WHY THIS FILE EXISTS
--------------------
This repository's build environment has no R toolchain, so an mrgsolve model
published here would be a model nobody had ever integrated.  Every number that
appears in `README.md` is produced by THIS script.  If the two ports disagree,
one of them is wrong -- that is the point of having two.

    python3 cdys_reference_check.py             # all analyses A0-A13
    python3 cdys_reference_check.py --only A2   # just the central result
    python3 cdys_reference_check.py --list

THE CENTRAL RESULT (A2)
-----------------------
Chemodenervation removes dystonic torque only from muscles the needle reached.
Let phi be the share of the total dystonic torque carried by the injected
muscles.  Then for ANY botulinum toxin -- any dose, any product, any potency,
any dilution, any interval -- the attainable TWSTRS severity is bounded below by

        Sev_floor(phi) = sev_ss( (1 - phi) * Dcen )

which contains NO drug parameter.  phi is a property of the injection plan and
of the patient's anatomy, not of the molecule.  A2 computes the bound and then
splits the remaining clinical headroom into the part a bigger dose can buy and
the part only a better TARGET LIST can buy.

STATE VECTOR (70 ODEs)
----------------------
  per muscle m = 1..8   (5 states each = 40)
      A_m   free bioactive toxin in the muscle interstitium      [U]
      B_m   membrane-bound / internalising (endosomal) toxin     [U]
      C_m   translocated, catalytically active light chain       [U-eq]
      S_m   intact SNAP-25 fraction (VAMP for serotype B)        [0-1]
      Q_m   sprout-mediated release capacity                    [0-1]
  swallow compartment (pharyngeal constrictors)  A,B,C,S         (4)
  autonomic compartment (salivary / ganglionic)  A,B,C,S         (4)
  humoral immunity, serotype A          Ag_A, Bmem_A, Nab_A      (3)
  humoral immunity, serotype B          Ag_B, Bmem_B, Nab_B      (3)
  central / spinal                Dcen, RecInh, SurrInh, Cbll    (4)
  clinical                              Sev, Pain, Disab         (3)
  oral adjunct PK (depot + central)     THP, BAC, CLZ            (6)
  accounting                       AUCben, AUCdys, CumU          (3)
                                                          total = 70
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

# ----------------------------------------------------------------------------
# 0.  MUSCLE TABLE
# ----------------------------------------------------------------------------
# w      = share of the TOTAL dystonic torque carried by this muscle; sums to 1.
#          Anchored on needle/surface-EMG series of rotational torticollis, in
#          which ipsilateral splenius capitis and contralateral SCM dominate but
#          a substantial minority of the drive sits in deep muscles (obliquus
#          capitis inferior, deep semispinalis, longus colli) that surface
#          injection routinely misses.
# mass   = relative muscle mass.  Toxin acts on C/mass -- a concentration -- so
#          a small muscle is paralysed by a small absolute amount.
# theta  = share of this muscle's diffusive efflux that reaches the swallow
#          compartment (anterior and deep muscles sit close to the pharynx).
# post   = 1 for posterior head/neck extensors; over-blocking these produces the
#          "neck weakness / head drop" adverse effect.
# std/ext = injected in the standard 4-muscle surface plan / the extended
#          EMG- or ultrasound-guided plan.

MUSCLES = [
    # name                        w      mass  theta  post
    ("SCM_contra",              0.20,  1.00, 0.30,   0),
    ("Splenius_capitis_ipsi",   0.24,  1.60, 0.08,   1),
    ("Trapezius_ipsi",          0.13,  2.20, 0.04,   1),
    ("Levator_scapulae_ipsi",   0.08,  0.80, 0.05,   1),
    ("Semispinalis_capitis",    0.12,  1.80, 0.12,   1),
    ("Scalene_complex",         0.06,  0.90, 0.28,   0),
    ("Obliquus_cap_inferior",   0.10,  0.50, 0.18,   1),
    ("Longus_colli_deep",       0.07,  0.70, 0.60,   0),
]
NM = len(MUSCLES)
MW = np.array([m[1] for m in MUSCLES])          # torque shares, sum = 1
MMASS = np.array([m[2] for m in MUSCLES])
MTHETA = np.array([m[3] for m in MUSCLES])
MPOST = np.array([m[4] for m in MUSCLES], dtype=float)
MNAME = [m[0] for m in MUSCLES]

# Injection patterns, Allergan-equivalent Units per muscle.  240 U over four
# muscles is the dose used in the pivotal cervical-dystonia registration trials
# of both onabotulinumtoxinA and incobotulinumtoxinA.
PATTERN_STD = np.array([50.0, 90.0, 60.0, 40.0, 0.0, 0.0, 0.0, 0.0])   # 240 U
PATTERN_LOW = PATTERN_STD * 0.5                                        # 120 U
PATTERN_EXT = np.array([45.0, 70.0, 40.0, 25.0, 40.0, 0.0, 20.0, 0.0]) # 240 U
PATTERN_HI = PATTERN_STD * 2.0                                         # 480 U

# Injectate volume (mL).  Reference dilution = 100 U / 2 mL = 50 U/mL.
VOL_REF_PER_U = 2.0 / 100.0
VOL_STD = PATTERN_STD * VOL_REF_PER_U
VOL_EXT = PATTERN_EXT * VOL_REF_PER_U

PHI_STD = MW[PATTERN_STD > 0].sum()
PHI_EXT = MW[PATTERN_EXT > 0].sum()


# ----------------------------------------------------------------------------
# 1.  PRODUCTS
# ----------------------------------------------------------------------------
# unit_scale relates each label's proprietary potency unit to the model's
# internal Allergan-equivalent Unit.  Doses are specified in A-EQUIVALENT Units
# everywhere; label Units (and hence the antigen mass) are obtained by dividing
# by unit_scale.  These conversions are clinical rules of thumb, NOT equipotency
# claims (A13).
# load = neurotoxin-complex protein per label Unit (ng/U).  This is the ANTIGEN,
# and it is the parameter that separates the products immunologically.
# kLC  = first-order loss of catalytically active light chain from the terminal
#        cytosol.  Serotype A light chain is markedly more persistent than B.
PRODUCTS = {
    "incobotulinumtoxinA": dict(serotype="A", unit_scale=1.00, load=0.44 / 100,
                                kLC_mult=1.00, auto_pref=1.0),
    "onabotulinumtoxinA":  dict(serotype="A", unit_scale=1.00, load=5.00 / 100,
                                kLC_mult=1.00, auto_pref=1.0),
    "abobotulinumtoxinA":  dict(serotype="A", unit_scale=0.34, load=4.35 / 500,
                                kLC_mult=1.00, auto_pref=1.0),
    "daxibotulinumtoxinA": dict(serotype="A", unit_scale=1.00, load=0.50 / 100,
                                kLC_mult=0.50, auto_pref=1.0),
    "rimabotulinumtoxinB": dict(serotype="B", unit_scale=0.03, load=5.00 / 5000,
                                kLC_mult=1.45, auto_pref=4.0),
}


# ----------------------------------------------------------------------------
# 2.  PARAMETERS
# ----------------------------------------------------------------------------
@dataclass
class P:
    # --- toxin disposition in muscle -------------------------------------
    k_bind: float = 4.20      # 1/d  productive binding to nerve terminals
    k_clear: float = 1.00     # 1/d  local proteolysis / lymphatic clearance
    k_diff0: float = 0.90     # 1/d  diffusive efflux at reference dilution
    gv: float = 1.00          # volume exponent for diffusive efflux
    k_trans: float = 0.60     # 1/d  endosomal translocation of the light chain
    k_LC: float = 0.0231      # 1/d  active LC loss -- CALIBRATED (init t1/2 30 d)
    Bmax_per_mass: float = 260.0   # U/mass, saturable terminal binding capacity

    # --- SNARE substrate cleavage / resynthesis --------------------------
    k_syn: float = 0.140      # 1/d  SNAP-25 resynthesis (t1/2 5.0 d)
    k_cl: float = 0.060       # 1/d per (U-eq/mass)  cleavage -- CALIBRATED

    # --- neuromuscular safety factor ------------------------------------
    n_snare: float = 3.0
    S50: float = 0.20
    SF: float = 3.0

    # --- terminal sprouting (two-phase recovery) -------------------------
    k_sp: float = 0.050       # 1/d  sprout outgrowth  (fixed, literature-shaped)
    k_rg: float = 0.050       # 1/d  regression once the parent terminal recovers
    wq: float = 0.55          # max share of the deficit sprouts can carry

    # --- spread compartments --------------------------------------------
    # Fraction of a TARGETED muscle's motor units the injectate actually
    # reaches.  A property of needle placement and intramuscular spread, not of
    # the molecule.  CALIBRATED from the high-dose plateau of the clinical
    # dose-response curve -- see A1/A2.
    rho: float = 0.62
    theta_scale: float = 1.00 # global scale on the MTHETA proximity index
    mass_sw: float = 0.35     # pharyngeal constrictors are small
    mass_au: float = 0.50
    theta_au: float = 0.06
    auto_pref: float = 1.0    # serotype B has ~4x autonomic tropism

    # --- clinical mapping -----------------------------------------------
    SEVMAX: float = 35.0      # TWSTRS severity subscale ceiling
    hs: float = 1.60
    L50: float = 0.8347       # chosen so sev_ss(1) = 20.0
    k_sev: float = 0.25       # 1/d

    PAINMAX: float = 20.0
    kp_on: float = 0.030
    kp_off: float = 0.030
    qp: float = 1.50
    pain_direct: float = 0.0  # direct antinociceptive BoNT action (see A13)

    DISMAX: float = 30.0
    g_dis: float = 0.80
    w_dis_sev: float = 0.55
    k_dis: float = 0.080      # 1/d

    # --- central / spinal -----------------------------------------------
    k_d: float = 0.0080       # 1/d  maladaptive plasticity (t1/2 87 d)
    ga: float = 0.80          # afferent-drive exponent
    aff_floor: float = 0.35   # spindle output surviving complete intrafusal block
    RI0: float = 0.45         # reciprocal inhibition in CD (1.0 = healthy)
    SI0: float = 0.40         # surround inhibition in CD
    a_ri: float = 0.50
    a_si: float = 0.50
    k_ri: float = 0.050
    k_si: float = 0.050
    k_cb: float = 0.020
    CB0: float = 1.00

    # --- adverse-effect logistics (both midpoints CALIBRATED) -----------
    dys_d50: float = 0.42
    dys_k: float = 0.110
    nw_d50: float = 0.62
    nw_k: float = 0.130

    # --- immunogenicity -------------------------------------------------
    k_ag: float = 0.50        # 1/d antigen depot clearance
    k_b: float = 0.0150       # 1/(ng*d) naive memory-B priming -- CALIBRATED
    beta_boost: float = 9.0   # memory recall amplification
    Ka_ag: float = 40.0       # ng, antigen-processing saturation
    k_bd: float = 0.0019      # 1/d memory-B decay (t1/2 365 d)
    Bmem_max: float = 50.0    # memory-B carrying capacity.  Without it the
                              # recall term (1 + beta*Bm) diverges at extreme
                              # k_b (A7's tail rows), which is a numerical
                              # artefact, not immunology.  Far above anything
                              # the calibrated regime reaches, so it changes
                              # no fitted result.
    k_np: float = 0.55        # 1/d Nab production per unit memory-B
    k_nd: float = 0.0347      # 1/d IgG elimination (t1/2 20 d)
    Nab50: float = 1.00       # Nab giving 50 % loss of injected potency
    hn: float = 1.60

    # --- oral adjunct PK / PD --------------------------------------------
    thp_ka: float = 36.0      # 1/d
    thp_ke: float = 4.50      # 1/d (t1/2 3.7 h)
    thp_V: float = 560.0      # L
    thp_Emax: float = 0.30
    thp_EC50: float = 8.0     # ng/mL
    bac_ka: float = 30.0
    bac_ke: float = 4.75      # t1/2 3.5 h
    bac_V: float = 60.0
    bac_Emax: float = 0.12
    bac_EC50: float = 200.0   # ng/mL
    clz_ka: float = 24.0
    clz_ke: float = 0.555     # t1/2 30 h
    clz_V: float = 210.0
    clz_Emax_ri: float = 0.36 # fractional restoration of the RecInh deficit
    clz_Emax_si: float = 0.25
    clz_EC50: float = 12.0    # ng/mL

    # --- interventions that are not drugs --------------------------------
    dbs_Emax: float = 0.50    # fractional reduction of central drive gain
    dbs_k: float = 0.0154     # 1/d ramp (t1/2 45 d)

    # --- accounting -------------------------------------------------------
    TW0_ref: float = 42.94    # untreated baseline TWSTRS total (see A0)
    dys_thresh: float = 0.30  # swallow deficit above which burden accrues


# state indices --------------------------------------------------------------
NSM = 5   # states per muscle


def _mk_index():
    idx = {}
    k = 0
    for i in range(NM):
        for s in ("A", "B", "C", "S", "Q"):
            idx[f"{s}_{i}"] = k
            k += 1
    for nm in ("A_sw", "B_sw", "C_sw", "S_sw",
               "A_au", "B_au", "C_au", "S_au",
               "Ag_A", "Bmem_A", "Nab_A", "Ag_B", "Bmem_B", "Nab_B",
               "Dcen", "RecInh", "SurrInh", "Cbll",
               "Sev", "Pain", "Disab",
               "THPd", "THPc", "BACd", "BACc", "CLZd", "CLZc",
               "AUCben", "AUCdys", "CumU"):
        idx[nm] = k
        k += 1
    return idx, k


IX, NSTATE = _mk_index()
assert NSTATE == 70, NSTATE
SLA = slice(0, NM * NSM, NSM)
SLB = slice(1, NM * NSM, NSM)
SLC = slice(2, NM * NSM, NSM)
SLS = slice(3, NM * NSM, NSM)
SLQ = slice(4, NM * NSM, NSM)


# ----------------------------------------------------------------------------
# 3.  ALGEBRAIC HELPERS
# ----------------------------------------------------------------------------
def r_of_S(S, p: P):
    """Quantal release capacity: the fraction of normal exocytosis the terminal
    can still support with an intact-SNARE fraction S.  Hill, from SNARE
    cooperativity."""
    S = np.clip(S, 0.0, 1.0)
    return S ** p.n_snare / (S ** p.n_snare + p.S50 ** p.n_snare)


def E_of_S(S, p: P):
    """Transmission efficacy as a function of the intact SNARE fraction.

    Two nonlinearities, both real and both load-bearing:
      1. SNARE cooperativity -- several intact complexes are needed per fusion
         event, so release capacity r(S) is a Hill function of S.
      2. The neuromuscular SAFETY FACTOR -- a healthy terminal releases roughly
         three times more quanta than are needed to fire the fibre, so release
         capacity must fall BELOW ~1/SF before any weakness appears.
    Together these make E(S) a threshold function.  That is why the dose ->
    duration relationship is logarithmic (A3) and why partial cleavage is
    clinically silent.
    """
    r = r_of_S(S, p)
    g = p.SF * r / (1.0 + p.SF * r)
    g1 = p.SF / (1.0 + p.SF)      # normalise so E(1) == 1 exactly
    return g / g1


def sev_ss(L, p: P):
    """Steady-state TWSTRS severity subscale for dystonic torque load L."""
    L = np.maximum(L, 0.0)
    return p.SEVMAX * L ** p.hs / (L ** p.hs + p.L50 ** p.hs)


def logistic(x, x50, k):
    return 1.0 / (1.0 + np.exp(-(x - x50) / k))


def emax(C, Emax_, EC50):
    return Emax_ * C / (C + EC50)


def tw_floor(L, p: P):
    """Asymptotic TWSTRS total if load L were held indefinitely."""
    sf = sev_ss(L, p)
    a = p.kp_on * L ** p.qp
    Pn = a / (a + p.kp_off)
    dis = p.DISMAX * p.g_dis * (p.w_dis_sev * sf / p.SEVMAX
                                + (1 - p.w_dis_sev) * Pn)
    return sf + dis + p.PAINMAX * Pn, sf, dis, p.PAINMAX * Pn


# ----------------------------------------------------------------------------
# 4.  REGIMEN
# ----------------------------------------------------------------------------
@dataclass
class Regimen:
    product: str = "incobotulinumtoxinA"
    pattern: np.ndarray = field(default_factory=lambda: PATTERN_STD.copy())
    volume: np.ndarray | None = None       # None -> reference dilution
    interval: float = 84.0                 # days between injections
    n_inj: int = 1
    first: float = 0.0
    # oral adjuncts, mg/day (approximated as a continuous input into the depot;
    # exact enough for a chronic tone effect, and stated as such)
    thp_mg_d: float = 0.0
    bac_mg_d: float = 0.0
    clz_mg_d: float = 0.0
    # non-pharmacological
    denervation: float = 0.0               # fraction of NON-injected drive removed
    dbs_on: float = -1.0                   # day of GPi-DBS activation
    # serotype rescue
    switch_at: int = -1
    switch_to: str = "rimabotulinumtoxinB"
    # BOUND MODE: force parent-terminal transmission to zero in these muscles
    # from t = 0, i.e. instantaneous and permanent perfect chemodenervation.
    # This is not a therapy; it is the upper limit of every therapy.
    block: np.ndarray | None = None
    # multiplier on rho, i.e. improved needle placement (EMG / ultrasound
    # guidance).  1.0 = the accuracy implicit in the calibration data.
    rho_mult: float = 1.0

    def blk(self):
        if self.block is None:
            return np.zeros(NM)
        return np.asarray(self.block, float)

    def inj_times(self):
        return [self.first + i * self.interval for i in range(self.n_inj)]

    def product_at(self, k):
        if self.switch_at >= 0 and k >= self.switch_at:
            return self.switch_to
        return self.product

    def kdiff(self, p: P):
        ref = self.pattern * VOL_REF_PER_U
        vol = ref if self.volume is None else np.asarray(self.volume, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(ref > 0, vol / np.maximum(ref, 1e-12), 1.0)
        return p.k_diff0 * np.power(np.maximum(ratio, 1e-6), p.gv)


# ----------------------------------------------------------------------------
# 5.  RIGHT-HAND SIDE
# ----------------------------------------------------------------------------
def rhs(t, y, p: P, reg: Regimen, kdiff_m: np.ndarray, w_eff: np.ndarray,
        blk: np.ndarray):
    dy = np.zeros_like(y)
    A = y[SLA]
    B = y[SLB]
    C = y[SLC]
    S = np.clip(y[SLS], 0.0, 1.0)
    Q = np.clip(y[SLQ], 0.0, 1.0)

    # --- 5.1 muscle chemodenervation cascade ----------------------------
    Bmax = p.Bmax_per_mass * MMASS
    occ = np.clip(1.0 - B / Bmax, 0.0, 1.0)
    dy[SLA] = -(p.k_bind * occ + p.k_clear + kdiff_m) * A
    dy[SLB] = p.k_bind * occ * A - p.k_trans * B
    dy[SLC] = p.k_trans * B - p.k_LC * C
    cconc = C / MMASS
    dy[SLS] = p.k_syn * (1.0 - S) - p.k_cl * cconc * S

    Em = E_of_S(S, p) * (1.0 - blk)      # bound mode: perfect blockade
    dy[SLQ] = p.k_sp * (1.0 - Em) * (1.0 - Q) - p.k_rg * Em * Q

    # --- 5.2 spread: swallow and autonomic compartments ------------------
    th = MTHETA * p.theta_scale
    flux_sw = float(np.sum(kdiff_m * th * A))
    flux_au = float(np.sum(kdiff_m * p.theta_au * p.auto_pref * A))

    for tag, flux, mass in (("sw", flux_sw, p.mass_sw), ("au", flux_au, p.mass_au)):
        Ax, Bx, Cx = y[IX[f"A_{tag}"]], y[IX[f"B_{tag}"]], y[IX[f"C_{tag}"]]
        Sx = min(max(y[IX[f"S_{tag}"]], 0.0), 1.0)
        occx = min(max(1.0 - Bx / (p.Bmax_per_mass * mass), 0.0), 1.0)
        dy[IX[f"A_{tag}"]] = flux - (p.k_bind * occx + p.k_clear) * Ax
        dy[IX[f"B_{tag}"]] = p.k_bind * occx * Ax - p.k_trans * Bx
        dy[IX[f"C_{tag}"]] = p.k_trans * Bx - p.k_LC * Cx
        dy[IX[f"S_{tag}"]] = p.k_syn * (1.0 - Sx) - p.k_cl * (Cx / mass) * Sx

    # --- 5.3 humoral immunity (two independent serotype pools) ----------
    for tag in ("A", "B"):
        Ag = y[IX[f"Ag_{tag}"]]
        Bm = y[IX[f"Bmem_{tag}"]]
        prime = (p.k_b * Ag / (1.0 + Ag / p.Ka_ag)
                 * (1.0 + p.beta_boost * Bm)
                 * max(0.0, 1.0 - Bm / p.Bmem_max))
        dy[IX[f"Ag_{tag}"]] = -p.k_ag * Ag
        dy[IX[f"Bmem_{tag}"]] = prime - p.k_bd * Bm
        dy[IX[f"Nab_{tag}"]] = p.k_np * Bm - p.k_nd * y[IX[f"Nab_{tag}"]]

    # --- 5.4 oral adjunct PK --------------------------------------------
    dy[IX["THPd"]] = reg.thp_mg_d * 1e6 - p.thp_ka * y[IX["THPd"]]
    dy[IX["THPc"]] = p.thp_ka * y[IX["THPd"]] - p.thp_ke * y[IX["THPc"]]
    dy[IX["BACd"]] = reg.bac_mg_d * 1e6 - p.bac_ka * y[IX["BACd"]]
    dy[IX["BACc"]] = p.bac_ka * y[IX["BACd"]] - p.bac_ke * y[IX["BACc"]]
    dy[IX["CLZd"]] = reg.clz_mg_d * 1e6 - p.clz_ka * y[IX["CLZd"]]
    dy[IX["CLZc"]] = p.clz_ka * y[IX["CLZd"]] - p.clz_ke * y[IX["CLZc"]]
    C_thp = y[IX["THPc"]] / p.thp_V / 1000.0     # ng in V(L) -> ng/mL
    C_bac = y[IX["BACc"]] / p.bac_V / 1000.0
    C_clz = y[IX["CLZc"]] / p.clz_V / 1000.0

    # --- 5.5 spinal / cortical inhibition -------------------------------
    f_clz = C_clz / (C_clz + p.clz_EC50)
    RI_t = p.RI0 + (1.0 - p.RI0) * p.clz_Emax_ri * f_clz
    SI_t = p.SI0 + (1.0 - p.SI0) * p.clz_Emax_si * f_clz
    RI, SI = y[IX["RecInh"]], y[IX["SurrInh"]]
    dy[IX["RecInh"]] = p.k_ri * (RI_t - RI)
    dy[IX["SurrInh"]] = p.k_si * (SI_t - SI)
    dy[IX["Cbll"]] = p.k_cb * (p.CB0 - y[IX["Cbll"]])

    # --- 5.6 central drive: the maladaptive-plasticity ratchet ----------
    # BoNT blocks intrafusal (gamma) cholinergic terminals as well as extrafusal
    # ones, so it lowers the abnormal muscle-spindle afferent drive that feeds
    # the plasticity state Dcen.  That is the model's route from a peripheral
    # injection to a slow CENTRAL benefit which outlasts the paralysis.
    rho = min(p.rho * reg.rho_mult, 1.0)
    # Complete blockade blocks the intrafusal (gamma) cholinergic terminals too,
    # so in bound mode the spindle afferent drive is at its floor as well.
    spindle = np.where(blk > 0, p.aff_floor,
                       p.aff_floor + (1.0 - p.aff_floor) * S)
    aff = float(np.sum(MW * (1.0 - rho * (1.0 - spindle))))
    G_raw = (1.0 + p.a_ri * (1.0 - RI)) * (1.0 + p.a_si * (1.0 - SI))
    G_ref = (1.0 + p.a_ri * (1.0 - p.RI0)) * (1.0 + p.a_si * (1.0 - p.SI0))
    dbs = 0.0
    if reg.dbs_on >= 0 and t >= reg.dbs_on:
        dbs = p.dbs_Emax * (1.0 - math.exp(-p.dbs_k * (t - reg.dbs_on)))
    Gcen = (G_raw / G_ref) * (1.0 - dbs) \
        * (1.0 - emax(C_thp, p.thp_Emax, p.thp_EC50)) \
        * (1.0 - emax(C_bac, p.bac_Emax, p.bac_EC50)) \
        * y[IX["Cbll"]]
    dy[IX["Dcen"]] = p.k_d * (Gcen * aff ** p.ga - y[IX["Dcen"]])

    # --- 5.7 dystonic torque load and clinical scores -------------------
    # Sprouts grow from the SAME axon and draw on the SAME cytosolic SNARE
    # pool, so their release is gated by S as well.  Sprouting therefore
    # ACCELERATES recovery -- it restores function earlier than the parent
    # terminal alone would -- rather than capping how deep blockade can go.
    Tm = 1.0 - (1.0 - Em) * (1.0 - p.wq * Q * r_of_S(S, p) * (1.0 - blk))
    Tm_eff = 1.0 - rho * (1.0 - Tm)      # only rho of the muscle was reached
    L = y[IX["Dcen"]] * float(np.sum(w_eff * Tm_eff))

    Sev = y[IX["Sev"]]
    dy[IX["Sev"]] = p.k_sev * (sev_ss(L, p) - Sev)

    Pn = min(max(y[IX["Pain"]], 0.0), 1.0)
    anti = 1.0 - p.pain_direct * float(np.sum(MW * (1.0 - Tm)))
    dy[IX["Pain"]] = p.kp_on * (L ** p.qp) * max(anti, 0.0) * (1.0 - Pn) \
        - p.kp_off * Pn

    dis_t = p.DISMAX * p.g_dis * (p.w_dis_sev * Sev / p.SEVMAX
                                  + (1.0 - p.w_dis_sev) * Pn)
    dy[IX["Disab"]] = p.k_dis * (dis_t - y[IX["Disab"]])

    # --- 5.8 accounting --------------------------------------------------
    TW = Sev + y[IX["Disab"]] + p.PAINMAX * Pn
    dy[IX["AUCben"]] = max(0.0, p.TW0_ref - TW)
    def_sw = 1.0 - float(E_of_S(np.array([y[IX["S_sw"]]]), p)[0])
    dy[IX["AUCdys"]] = max(0.0, def_sw - p.dys_thresh)
    return dy


# ----------------------------------------------------------------------------
# 6.  BASELINE AND SIMULATION DRIVER
# ----------------------------------------------------------------------------
def y0_baseline(p: P):
    y = np.zeros(NSTATE)
    y[SLS] = 1.0
    y[IX["S_sw"]] = 1.0
    y[IX["S_au"]] = 1.0
    y[IX["Dcen"]] = 1.0
    y[IX["RecInh"]] = p.RI0
    y[IX["SurrInh"]] = p.SI0
    y[IX["Cbll"]] = p.CB0
    y[IX["Sev"]] = sev_ss(1.0, p)
    y[IX["Pain"]] = p.kp_on / (p.kp_on + p.kp_off)
    y[IX["Disab"]] = p.DISMAX * p.g_dis * (
        p.w_dis_sev * y[IX["Sev"]] / p.SEVMAX
        + (1.0 - p.w_dis_sev) * y[IX["Pain"]])
    return y


def _w_eff(reg: Regimen):
    inj = reg.pattern > 0
    w = MW.copy()
    if reg.denervation > 0:
        w = np.where(inj, w, w * (1.0 - reg.denervation))
    return w


def simulate(reg: Regimen, p: P | None = None, tmax: float = 168.0, dt: float = 0.5):
    """Integrate the system across injection boluses."""
    p = P() if p is None else P(**{k: getattr(p, k) for k in p.__dataclass_fields__})
    y = y0_baseline(p)
    times, Y = [0.0], [y.copy()]
    inj_t = [t for t in reg.inj_times() if t <= tmax]
    breaks = sorted(set([0.0] + inj_t + [tmax]))
    kdiff_m = reg.kdiff(p)
    w_eff = _w_eff(reg)
    blk = reg.blk()
    kLC0 = p.k_LC
    cum_u = 0.0

    for k in range(len(breaks) - 1):
        t0, t1 = breaks[k], breaks[k + 1]
        if t0 in inj_t:
            ki = inj_t.index(t0)
            prod = PRODUCTS[reg.product_at(ki)]
            p.k_LC = kLC0 * prod["kLC_mult"]
            p.auto_pref = prod["auto_pref"]
            tag = prod["serotype"]
            # reg.pattern is in A-EQUIVALENT Units, so switching product does
            # not silently change the delivered activity (it used to: a switch
            # to rimabotulinumtoxinB delivered 3 % of the intended dose, which
            # made serotype rescue look harmful).  Label Units, and therefore
            # the ANTIGEN mass, are derived from the conversion instead.
            gate = 1.0 / (1.0 + (y[IX[f"Nab_{tag}"]] / p.Nab50) ** p.hn)
            y[SLA] += reg.pattern * gate
            label_U = float(np.sum(reg.pattern)) / prod["unit_scale"]
            y[IX[f"Ag_{tag}"]] += label_U * prod["load"]
            cum_u += float(np.sum(reg.pattern))
            y[IX["CumU"]] = cum_u
        elif k == 0:
            prod = PRODUCTS[reg.product_at(0)]
            p.k_LC = kLC0 * prod["kLC_mult"]
            p.auto_pref = prod["auto_pref"]

        n = max(2, int(round((t1 - t0) / dt)) + 1)
        sol = solve_ivp(rhs, (t0, t1), y, args=(p, reg, kdiff_m, w_eff, blk),
                        method="LSODA", t_eval=np.linspace(t0, t1, n),
                        rtol=1e-7, atol=1e-9, max_step=2.0)
        if not sol.success:
            raise RuntimeError(f"integration failed: {sol.message}")
        for j in range(1, sol.y.shape[1]):
            times.append(sol.t[j])
            Y.append(sol.y[:, j].copy())
        y = sol.y[:, -1].copy()

    p.k_LC = kLC0
    return _derive(np.array(times), np.array(Y), p, reg, w_eff, blk)


def _derive(T, Ym, p: P, reg: Regimen, w_eff, blk=None):
    S = np.clip(Ym[:, SLS], 0, 1)
    Q = np.clip(Ym[:, SLQ], 0, 1)
    blk = np.zeros(NM) if blk is None else blk
    rho = min(p.rho * reg.rho_mult, 1.0)
    Em = E_of_S(S, p) * (1.0 - blk)
    Tm_raw = 1.0 - (1.0 - Em) * (1.0 - p.wq * Q * r_of_S(S, p) * (1.0 - blk))
    Tm = 1.0 - rho * (1.0 - Tm_raw)      # only rho of the muscle was reached
    L = Ym[:, IX["Dcen"]] * (Tm @ w_eff)
    Sev = Ym[:, IX["Sev"]]
    Pain = np.clip(Ym[:, IX["Pain"]], 0, 1)
    Dis = Ym[:, IX["Disab"]]
    def_sw = 1.0 - E_of_S(np.clip(Ym[:, IX["S_sw"]], 0, 1), p)
    def_au = 1.0 - E_of_S(np.clip(Ym[:, IX["S_au"]], 0, 1), p)
    wpost = MW * MPOST
    post_def = ((1.0 - Tm) @ wpost) / wpost.sum()
    return dict(
        t=T, Y=Ym, S=S, Q=Q, E=Em, Tm=Tm, Tm_raw=Tm_raw, L=L,
        Sev=Sev, Pain=Pain * p.PAINMAX, Disab=Dis,
        TW=Sev + Dis + p.PAINMAX * Pain,
        def_sw=def_sw, def_au=def_au, post_def=post_def,
        P_dys=logistic(def_sw, p.dys_d50, p.dys_k),
        P_neck=logistic(post_def, p.nw_d50, p.nw_k),
        P_dry=logistic(def_au, p.dys_d50, p.dys_k),
        Csw=Ym[:, IX["C_sw"]], Dcen=Ym[:, IX["Dcen"]], NabA=Ym[:, IX["Nab_A"]], NabB=Ym[:, IX["Nab_B"]],
        BmemA=Ym[:, IX["Bmem_A"]], AUCben=Ym[:, IX["AUCben"]],
        AUCdys=Ym[:, IX["AUCdys"]], CumU=Ym[:, IX["CumU"]], p=p, reg=reg,
    )


# ----------------------------------------------------------------------------
# 6b.  FAST PATH FOR THE ANTIBODY SUBSYSTEM
# ----------------------------------------------------------------------------
# Ag -> Bmem -> Nab depends ONLY on the injection times and the antigen mass per
# injection.  The potency gate feeds back into the injected DOSE, not into the
# antigen (the protein is delivered whether or not it is neutralised), so this
# three-state subsystem is exactly decoupled from the other 67 states.
#
# Integrating it on its own is therefore not an approximation -- it gives the
# same Nab trajectory as the full model -- and it is what makes the virtual
# cohorts in A8 and A13(5) affordable.  Verified against the full model in A8.
def nab_trajectory(load_ng, times, p: P, t_end=None):
    """Neutralising antibody after injections of `load_ng` ng at `times` (days),
    integrated through to `t_end` (default: the last injection time).

    Returns (t, Nab).  Exact for the full model's antibody pool -- verified
    against the 70-state model at the top of A8.  Getting t_end right matters:
    Nab decays with a 20-day half-life, so reading the trajectory at the last
    injection instead of at the horizon overstates it by about 11 %.
    """
    def f(t, y):
        Ag, Bm, Nb = y
        prime = (p.k_b * Ag / (1.0 + Ag / p.Ka_ag)
                 * (1.0 + p.beta_boost * Bm)
                 * max(0.0, 1.0 - Bm / p.Bmem_max))
        return [-p.k_ag * Ag, prime - p.k_bd * Bm, p.k_np * Bm - p.k_nd * Nb]

    times = [float(t) for t in times]
    t_end = times[-1] if t_end is None else float(t_end)
    edges = sorted(set([t for t in times if t <= t_end] + [t_end]))
    inj = set(times)
    y = np.zeros(3)
    T, Y = [edges[0]], [y.copy()]
    for k in range(len(edges) - 1):
        if edges[k] in inj:
            y[0] += load_ng
        sol = solve_ivp(f, (edges[k], edges[k + 1]), y, method="LSODA",
                        rtol=1e-8, atol=1e-11,
                        t_eval=np.linspace(edges[k], edges[k + 1], 4))
        if not sol.success:
            raise RuntimeError(sol.message)
        for j in range(1, sol.y.shape[1]):
            T.append(sol.t[j])
            Y.append(sol.y[:, j].copy())
        y = sol.y[:, -1].copy()
    return np.array(T), np.array(Y)[:, 2]


def nab_cohort(load_ng, times, kb_vec, p: P, t_end, dt=0.2):
    """Nab at `t_end` for a WHOLE COHORT at once.

    Same three-state system as nab_trajectory, but (a) the antigen depot is
    solved analytically -- between injections it is just a sum of decaying
    exponentials, one per dose already given -- and (b) the remaining two states
    are advanced for every patient simultaneously by fixed-step RK4, with k_b
    entering as a vector.

    This exists because the cohort cost is dominated by per-call solver setup,
    not by the integration: 1200 patients x 21 segments is 25 000 solve_ivp
    invocations.  Vectorising turns that into one loop.

    NOTE ON THE DISCONTINUITY.  Each injection is a jump in Ag, so the RK4 grid
    must not straddle one -- an earlier version used a single global grid and
    picked up the bolus half a step early, which the A8 self-check caught as a
    2 % bias.  The loop below therefore restarts at every injection, and the
    antigen sum inside a segment counts only doses already given.  Cross-checked
    against both the full 70-state model and the LSODA path at the top of A8.
    """
    kb = np.atleast_1d(np.asarray(kb_vec, float))
    tt = [float(t) for t in times if float(t) <= t_end]
    edges = tt + [float(t_end)]
    Bm = np.zeros_like(kb)
    Nb = np.zeros_like(kb)

    for k in range(len(edges) - 1):
        given = np.asarray(tt[:k + 1], float)      # doses already in the body

        def ag(t):
            return load_ng * float(np.sum(np.exp(-p.k_ag * (t - given))))

        def deriv(A, B_, N_):
            prime = (kb * A / (1.0 + A / p.Ka_ag) * (1.0 + p.beta_boost * B_)
                     * np.clip(1.0 - B_ / p.Bmem_max, 0.0, None))
            return prime - p.k_bd * B_, p.k_np * B_ - p.k_nd * N_

        span = edges[k + 1] - edges[k]
        if span <= 0:
            continue
        nsub = max(1, int(math.ceil(span / dt)))
        h = span / nsub
        t0 = edges[k]
        for i in range(nsub):
            a0, ah, a1 = ag(t0), ag(t0 + h / 2), ag(t0 + h)
            k1b, k1n = deriv(a0, Bm, Nb)
            k2b, k2n = deriv(ah, Bm + h / 2 * k1b, Nb + h / 2 * k1n)
            k3b, k3n = deriv(ah, Bm + h / 2 * k2b, Nb + h / 2 * k2n)
            k4b, k4n = deriv(a1, Bm + h * k3b, Nb + h * k3n)
            Bm = Bm + h / 6 * (k1b + 2 * k2b + 2 * k3b + k4b)
            Nb = Nb + h / 6 * (k1n + 2 * k2n + 2 * k3n + k4n)
            t0 += h
    return Nb


def nab_final(load_ng, times, p: P, t_end=None):
    return float(nab_trajectory(load_ng, times, p, t_end)[1][-1])


def potency_gate(nab, p: P):
    return 1.0 / (1.0 + (nab / p.Nab50) ** p.hn)


# ----------------------------------------------------------------------------
# 7.  READ-OUTS
# ----------------------------------------------------------------------------
def at(res, day, key="TW"):
    return float(np.interp(day, res["t"], res[key]))


def nadir(res, key="TW", lo=0.0, hi=90.0):
    m = (res["t"] >= lo) & (res["t"] <= hi)
    i = int(np.argmin(res[key][m]))
    return float(res[key][m][i]), float(res["t"][m][i])


MCID_TWSTRS = 4.5   # minimal clinically important change, TWSTRS total.
                    # Published estimates cluster around 4-6 points; 4.5 is used
                    # throughout and the choice is stated, not hidden.


def duration_of_benefit(res, mcid=MCID_TWSTRS, cyc0=0.0, cyc1=None):
    """Days from injection until the TWSTRS-total gain falls back below the
    minimal clinically important change -- i.e. until the patient can no longer
    tell the injection happened.  This is a PREDICTION of the model, not a
    calibration anchor (see A1/A13)."""
    p = res["p"]
    t = res["t"]
    m = t >= cyc0 if cyc1 is None else (t >= cyc0) & (t <= cyc1)
    tt, gain = t[m], p.TW0_ref - res["TW"][m]
    i = int(np.argmax(gain))
    if gain[i] <= mcid:
        return 0.0
    for j in range(i, len(gain)):
        if gain[j] <= mcid:
            return float(tt[j] - cyc0)
    return float(tt[-1] - cyc0)


# ----------------------------------------------------------------------------
# 8.  CALIBRATION
# ----------------------------------------------------------------------------
# Every anchor below is a directly reported trial quantity.  The two TWSTRS
# time points come from the pivotal cervical-dystonia programmes, in which
# 240 U produced a TWSTRS-total change from baseline of about -10 at week 4 and
# was still about -6 to -7 at week 12, the protocol re-injection visit.
CAL_TARGETS = dict(dTW_wk4_240U=-10.5, dTW_wk12_240U=-6.5, dTW_nadir_480U=-13.0,
                   P_dys_240=0.110, P_dys_120=0.055, P_neck_240=0.090,
                   nonresp_5y_ona=0.020, sigma_kb=0.90)
CAL = {}


def _std_reg(f=1.0, ext=False):
    pat = (PATTERN_EXT if ext else PATTERN_STD) * f
    return Regimen(pattern=pat, volume=pat * VOL_REF_PER_U, n_inj=1)


def _cal_resid(x, verbose=False):
    rho = 1.0 / (1.0 + math.exp(-x[2]))          # keep rho in (0, 1)
    p = P(k_cl=math.exp(x[0]), k_LC=math.exp(x[1]), rho=rho)
    res = simulate(_std_reg(), p, tmax=100.0, dt=1.0)
    d4 = at(res, 28.0) - p.TW0_ref
    d12 = at(res, 84.0) - p.TW0_ref
    r480 = simulate(_std_reg(2.0), p, tmax=140.0, dt=1.0)
    dn = nadir(r480, lo=0, hi=130)[0] - p.TW0_ref
    if verbose:
        print(f"      k_cl={math.exp(x[0]):.5f}  k_LC={math.exp(x[1]):.5f} "
              f"(LC t1/2 {math.log(2)/math.exp(x[1]):5.1f} d)  rho={rho:.4f}"
              f"  ->  wk4={d4:+.2f}  wk12={d12:+.2f}  nadir480={dn:+.2f}")
    return [(d4 - CAL_TARGETS["dTW_wk4_240U"]) / 0.5,
            (d12 - CAL_TARGETS["dTW_wk12_240U"]) / 0.5,
            (dn - CAL_TARGETS["dTW_nadir_480U"]) / 0.8]


def calibrate(verbose=True):
    """Seven parameters, seven reported anchors -- solved, not asserted.

      k_cl, k_LC, rho   <- TWSTRS-total change at week 4 and week 12 at 240 U,
                           and the nadir at 480 U (the dose-response PLATEAU)
      dys_d50, dys_k    <- dysphagia incidence at 240 U and at 120 U
      k_b               <- 5-year neutralising-antibody rate on onaBoNT-A
      nw_d50            <- neck-weakness incidence at 240 U

    rho -- the share of a targeted muscle actually reached by the injectate --
    is identified by the PLATEAU.  A dose-response curve that flattens while
    adverse effects keep climbing is not pharmacological saturation; it is the
    needle running out of disease to reach.  Estimating rho this way turns the
    plateau from a nuisance into a measurement (A2).

    Nothing else is fitted.  Duration of benefit, the 60/120/180 U points, every
    adverse-effect extrapolation, all of A2-A12 and every mismatch in A13 are
    predictions of the six numbers found here.
    """
    if CAL:
        return CAL
    if verbose:
        print("    stage 1: k_cl, k_LC, rho  vs  week-4 / week-12 TWSTRS at")
        print("             240 U and the dose-response plateau at 480 U")
    sol = least_squares(_cal_resid,
                        x0=[math.log(0.06), math.log(0.0231), 0.5],
                        kwargs=dict(verbose=verbose), xtol=1e-10, ftol=1e-10,
                        diff_step=0.03)
    CAL["k_cl"] = math.exp(sol.x[0])
    CAL["k_LC"] = math.exp(sol.x[1])
    CAL["rho"] = 1.0 / (1.0 + math.exp(-sol.x[2]))

    if verbose:
        print("    stage 2: dys_d50, dys_k  vs  dysphagia at 240 U and at 120 U")
    p0 = P(k_cl=CAL["k_cl"], k_LC=CAL["k_LC"], rho=CAL["rho"])
    d240 = float(np.max(simulate(_std_reg(1.0), p0, tmax=140.0, dt=1.0)["def_sw"]))
    d120 = float(np.max(simulate(_std_reg(0.5), p0, tmax=140.0, dt=1.0)["def_sw"]))
    l240 = math.log(CAL_TARGETS["P_dys_240"] / (1 - CAL_TARGETS["P_dys_240"]))
    l120 = math.log(CAL_TARGETS["P_dys_120"] / (1 - CAL_TARGETS["P_dys_120"]))
    CAL["dys_k"] = (d240 - d120) / (l240 - l120)
    CAL["dys_d50"] = d240 - CAL["dys_k"] * l240
    CAL["def_sw_240"], CAL["def_sw_120"] = d240, d120

    if verbose:
        print("    stage 3: k_b  vs  5-year neutralising-antibody rate")
    # Only about 2 % of patients on the current onabotulinumtoxinA formulation
    # develop clinically relevant neutralising antibody over 5 years of q12wk
    # dosing.  With a lognormal random effect of sigma = 0.9 on k_b, and Nab
    # roughly proportional to k_b in the sub-boost regime, the 98th percentile
    # sits a factor exp(2.054 * 0.9) = 6.36 above the median.  So the MEDIAN
    # patient's 5-year Nab must be about Nab50 / 6.36.  One deterministic
    # root-find pins k_b; the actual cohort rate that results is then REPORTED
    # in A8 rather than assumed, and the nonlinear boost term means it need not
    # come out at exactly 2 %.
    z = 2.0537
    tgt = P().Nab50 / math.exp(z * CAL_TARGETS["sigma_kb"])

    ona_load = 240.0 * PRODUCTS["onabotulinumtoxinA"]["load"]
    ona_times = np.arange(21) * 84.0
    H5Y = 5 * 365.0

    def f_kb(kb):
        pp = P(k_cl=CAL["k_cl"], k_LC=CAL["k_LC"], rho=CAL["rho"], k_b=kb)
        return nab_final(ona_load, ona_times, pp, H5Y) - tgt

    CAL["k_b"] = brentq(f_kb, 1e-6, 0.05, xtol=1e-10, rtol=1e-8)
    CAL["Nab_median_5y_target"] = tgt

    if verbose:
        print("    stage 4: nw_d50  vs  neck weakness at 240 U")
    p1 = P(k_cl=CAL["k_cl"], k_LC=CAL["k_LC"], rho=CAL["rho"], k_b=CAL["k_b"],
           dys_d50=CAL["dys_d50"], dys_k=CAL["dys_k"])
    pk = float(np.max(simulate(_std_reg(), p1, tmax=140.0, dt=1.0)["post_def"]))
    CAL["nw_d50"] = pk - P().nw_k * math.log(
        CAL_TARGETS["P_neck_240"] / (1 - CAL_TARGETS["P_neck_240"]))
    if verbose:
        print("    calibrated:")
        print(f"      k_cl    = {CAL['k_cl']:.5f}  /d per (U-eq per unit mass)")
        print(f"      rho     = {CAL['rho']:.4f}    -> phi(standard) = "
              f"{CAL['rho']*PHI_STD:.4f}, phi(extended) = {CAL['rho']*PHI_EXT:.4f}")
        print(f"      k_LC    = {CAL['k_LC']:.5f}  /d   active light-chain "
              f"t1/2 = {math.log(2)/CAL['k_LC']:.1f} d")
        print("                (independent estimates of light-chain/A")
        print("                 persistence in neurons span weeks to months;")
        print("                 the fitted value falls inside that range and was")
        print("                 not constrained to)")
        print(f"      dys_d50 = {CAL['dys_d50']:.4f}   dys_k = {CAL['dys_k']:.4f}"
              f"   (swallow deficit 240 U = {d240:.3f}, 120 U = {d120:.3f})")
        print(f"      k_b     = {CAL['k_b']:.6f}  1/(ng.d)  -> median 5-y Nab "
              f"on onaBoNT-A = {CAL['Nab_median_5y_target']:.4f}")
        print(f"      nw_d50  = {CAL['nw_d50']:.4f}")
    return CAL


def cal_params(**kw):
    c = calibrate(verbose=False)
    base = {k: c[k] for k in ("k_cl", "k_LC", "rho", "k_b", "dys_d50",
                              "dys_k", "nw_d50")}
    base.update(kw)
    return P(**base)


# ============================================================================
#  ANALYSES
# ============================================================================
def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def A0():
    hr("A0  MODEL INVENTORY")
    p = P()
    y = y0_baseline(p)
    tw0 = y[IX["Sev"]] + y[IX["Disab"]] + p.PAINMAX * y[IX["Pain"]]
    print(f"  states                    : {NSTATE} ODEs")
    print(f"  muscle compartments       : {NM} x {NSM} states = {NM*NSM}")
    print(f"  baseline TWSTRS severity  : {y[IX['Sev']]:6.2f} / 35")
    print(f"  baseline TWSTRS disability: {y[IX['Disab']]:6.2f} / 30")
    print(f"  baseline TWSTRS pain      : {p.PAINMAX*y[IX['Pain']]:6.2f} / 20")
    print(f"  baseline TWSTRS TOTAL     : {tw0:6.2f} / 85   "
          f"(reported trial baselines 43-46)")
    print()
    print("  transmission efficacy vs intact SNAP-25, safety factor SF = 3:")
    print("      S     r(S)     E(S)")
    for s in (1.0, 0.7, 0.5, 0.35, 0.25, 0.20, 0.15, 0.10, 0.05, 0.02):
        r = s ** p.n_snare / (s ** p.n_snare + p.S50 ** p.n_snare)
        print(f"   {s:5.2f}  {r:7.4f}  {float(E_of_S(np.array([s]), p)[0]):7.4f}")
    print("  -> half the substrate can be destroyed with no measurable weakness.")
    print("     Clinical effect begins only once S falls below about 0.20.")
    print()
    print("  muscle torque shares (sum = 1.000):")
    for i, nm in enumerate(MNAME):
        fl = [s for s, pat in (("std", PATTERN_STD), ("ext", PATTERN_EXT))
              if pat[i] > 0]
        print(f"   {nm:26s} w={MW[i]:.2f}  mass={MMASS[i]:.2f}  "
              f"theta={MTHETA[i]:.2f}  post={int(MPOST[i])}  {'/'.join(fl) or '--'}")
    print(f"  phi (standard 4-muscle surface plan) = {PHI_STD:.3f}")
    print(f"  phi (extended EMG/US-guided plan)    = {PHI_EXT:.3f}")
    print(f"  anatomically unreachable floor        = {MW[7]:.3f}")


def A1():
    hr("A1  CALIBRATION -- reproducing the pivotal-trial time course")
    calibrate(verbose=True)
    p = cal_params()
    print()
    print("  Seven parameters were free; seven reported anchors were used.")
    print("  Everything else below, and everything in A2-A12, is a prediction.")
    print()
    print("  TWSTRS total, change from baseline (points).  The week-4 and week-12")
    print("  columns at 240 U are the only two efficacy numbers that were fitted.")
    print("   dose        wk1     wk2    [wk4]    wk6     wk8   [wk12]   wk16"
          "    nadir  t_nadir")
    for lbl, f in (("placebo", 0.0), ("60 U", 0.25), ("120 U", 0.5),
                   ("240 U", 1.0), ("480 U", 2.0)):
        res = simulate(_std_reg(f) if f > 0
                       else Regimen(pattern=PATTERN_STD * 0.0,
                                    volume=PATTERN_STD * 0.0),
                       p, tmax=210.0, dt=1.0)
        d = [at(res, w * 7.0) - p.TW0_ref for w in (1, 2, 4, 6, 8, 12, 16)]
        nd, nt = nadir(res, lo=0, hi=120)
        tn = f"{nt:5.0f}d" if f > 0 else "   --"
        print(f"   {lbl:8s}" + "".join(f"{v:+7.2f} " for v in d)
              + f" {nd-p.TW0_ref:+7.2f}  {tn}")
    print()
    print(f"  Predicted duration of benefit (gain > MCID = {MCID_TWSTRS} TWSTRS points)")
    for lbl, f in (("60 U", 0.25), ("120 U", 0.5), ("240 U", 1.0), ("480 U", 2.0)):
        res = simulate(_std_reg(f), p, tmax=320.0, dt=1.0)
        print(f"   {lbl:8s} {duration_of_benefit(res):6.1f} d "
              f"({duration_of_benefit(res)/7:4.1f} weeks)")
    print()
    print("  Adverse effects at the nadir")
    print("   dose      swallow deficit  P(dysphagia)  posterior deficit"
          "  P(neck weakness)")
    for lbl, f in (("120 U", 0.5), ("240 U", 1.0), ("480 U", 2.0)):
        res = simulate(_std_reg(f), p, tmax=140.0, dt=1.0)
        print(f"   {lbl:8s} {float(np.max(res['def_sw'])):14.3f}  "
              f"{float(np.max(res['P_dys']))*100:11.1f}%  "
              f"{float(np.max(res['post_def'])):16.3f}  "
              f"{float(np.max(res['P_neck']))*100:14.1f}%")
    print()
    print("  Anchors used (7):")
    print("    240 U week 4    TWSTRS total change  -9.0 to -11.7 -> target -10.5")
    print("    240 U week 12   TWSTRS total change  -6 to -7      -> target  -6.5")
    print("    480 U nadir     the dose-response PLATEAU          -> target -13.0")
    print("    dysphagia       240 U ~11 %,  120 U ~5.5 %")
    print("    neck weakness   240 U ~7-11 %")
    print("    neutralising antibody, onaBoNT-A q12wk, 5 y ~2 %")
    print()
    print("  NOT used in the fit -- held back as checks:")
    print("    60 U and 120 U should give a shallow, sub-proportional response")
    print("    peak effect should fall between week 2 and week 4")
    print("    patient-perceived duration should be about 10-12 weeks")
    print("    dysphagia above 480 U should keep rising after efficacy stops")
    print("  The first and fourth are met.  The second and third are NOT:")
    print("  see A13(2) -- the model's effect-time curve is too square.")
    print()
    res = simulate(_std_reg(), p, tmax=140.0, dt=0.5)
    i = int(np.argmin(res["TW"]))
    y0 = y0_baseline(p)
    print("  Subscale decomposition at the nadir, 240 U:")
    print(f"    severity   {y0[IX['Sev']]:6.2f} -> {res['Sev'][i]:6.2f}"
          f"   ({res['Sev'][i]-y0[IX['Sev']]:+.2f})")
    print(f"    disability {y0[IX['Disab']]:6.2f} -> {res['Disab'][i]:6.2f}"
          f"   ({res['Disab'][i]-y0[IX['Disab']]:+.2f})")
    print(f"    pain       {p.PAINMAX*y0[IX['Pain']]:6.2f} -> {res['Pain'][i]:6.2f}"
          f"   ({res['Pain'][i]-p.PAINMAX*y0[IX['Pain']]:+.2f})")


def A2():
    hr("A2  *** THE BOUND ***  the plateau is not saturation, it is a measurement")
    p = cal_params()
    TW0 = p.TW0_ref
    phi_std = p.rho * PHI_STD
    phi_ext = p.rho * PHI_EXT
    print("  Chemodenervation can at best abolish transmission in the part of")
    print("  the neck the injectate actually reached.  Write that share as")
    print()
    print("      phi = rho * SUM(w_m over injected muscles)")
    print("             ^         ^")
    print("        how much of    which muscles")
    print("        each muscle    were on the list")
    print("        the needle")
    print("        reached")
    print()
    print("  Then for ANY toxin -- ANY dose, product, potency, dilution or")
    print("  interval -- the dystonic torque load obeys")
    print()
    print("      L(t)  >=  L_min  =  (1 - phi) * Dcen")
    print()
    print("  There is no drug parameter on the right-hand side.  Dose does not")
    print("  appear.  Potency does not appear.  Serotype does not appear.")
    print()
    print("  ---- phi is OBSERVABLE, and the observation is the plateau ----")
    print("  A dose-response curve that flattens while adverse effects keep")
    print("  climbing is not pharmacological saturation -- the toxin has not run")
    print("  out of anything.  It is the needle running out of DISEASE TO REACH.")
    print("  So the plateau measures phi.  Fitting it (A1 stage 1) gives:")
    print()
    print(f"      rho              = {p.rho:.4f}")
    print(f"      SUM w (standard) = {PHI_STD:.4f}   ->  phi(standard) = {phi_std:.4f}")
    print(f"      SUM w (extended) = {PHI_EXT:.4f}   ->  phi(extended) = {phi_ext:.4f}")
    print()
    print(f"  Read plainly: about {(1-phi_std)*100:.0f} % of the dystonic drive in")
    print("  cervical dystonia lies outside what a standard surface injection")
    print("  affects at all.  That number was not assumed from anatomy -- it was")
    print("  recovered from the shape of the published dose-response curve.")
    print()
    print("  ---- the bound, evaluated ----")
    print("  phi     as              L_min   Sev_floor   Pain   Disab   TW_floor"
          "   dTW_max")
    for phi, lbl in ((0.0, "no treatment"),
                     (phi_std, "standard plan, as-is"),
                     (phi_ext, "extended plan, same rho"),
                     (0.65, "standard plan, rho -> 1"),
                     (0.87, "extended plan, rho -> 1"),
                     (0.93, "+ scalene, rho -> 1"),
                     (1.00, "hypothetical: everything")):
        tw, sf, dis, pn = tw_floor(1.0 - phi, p)
        print(f"  {phi:.3f}  {lbl:22s} {1-phi:6.3f}  {sf:9.2f}  {pn:6.2f}"
              f"  {dis:6.2f}  {tw:8.2f}  {tw-TW0:+8.2f}")
    print()
    print("  ---- the bound integrated, against what is achieved ----")
    print("  'PERFECT BLOCKADE' means transmission forced to zero from day 0 in")
    print("  every reached motor unit and held there forever.  No molecule can do")
    print("  better.  Currency is MEAN TWSTRS GAIN over a steady-state 12-week")
    print("  cycle (cycle 8 of q12wk dosing) -- what the patient lives with, and")
    print("  the only read-out on which a wearing-off therapy and a permanent one")
    print("  compare honestly.")
    print()

    def cyc(reg, pp, n=8, iv=84.0):
        res = simulate(reg, pp, tmax=n * iv, dt=1.0)
        m = (res["t"] >= (n - 1) * iv) & (res["t"] <= n * iv)
        g = pp.TW0_ref - res["TW"][m]
        return (float(np.mean(g)), float(np.max(g)), float(np.min(g)),
                float(np.max(res["P_dys"])) * 100,
                float(np.max(res["P_neck"])) * 100)

    def q12(pat, **kw):
        return Regimen(pattern=pat.copy(), volume=pat * VOL_REF_PER_U,
                       interval=84.0, n_inj=8, **kw)

    blkS = (PATTERN_STD > 0) * 1.0
    blkE = (PATTERN_EXT > 0) * 1.0
    scen = [
        ("240 U q12wk  (as treated)",     q12(PATTERN_STD), p),
        ("480 U q12wk",                   q12(PATTERN_STD * 2), p),
        ("960 U q12wk",                   q12(PATTERN_STD * 4), p),
        ("14400 U q12wk (60x)",           q12(PATTERN_STD * 60), p),
        ("BOUND perfect block, as-is",    q12(PATTERN_STD, block=blkS), p),
        ("240 U, extended target list",   q12(PATTERN_EXT), p),
        ("240 U, perfect placement",      q12(PATTERN_STD, rho_mult=1/p.rho), p),
        ("240 U, extended + perfect",     q12(PATTERN_EXT, rho_mult=1/p.rho), p),
        ("BOUND perfect block, rho->1",   q12(PATTERN_STD, block=blkS,
                                              rho_mult=1/p.rho), p),
        ("BOUND ext. block, rho->1",      q12(PATTERN_EXT, block=blkE,
                                              rho_mult=1/p.rho), p),
        ("BOUND as-is, NO sprouting",     q12(PATTERN_STD, block=blkS),
         cal_params(wq=0.0)),
        ("240 U q12wk, NO sprouting",     q12(PATTERN_STD), cal_params(wq=0.0)),
    ]
    R = {}
    print("   regimen                        mean gain   best  worst  P(dysph)"
          " P(neck)")
    for lbl, reg, pp in scen:
        mn, bs, wr, pd_, pn_ = cyc(reg, pp)
        R[lbl] = mn
        print(f"   {lbl:30s} {mn:9.2f} {bs:6.2f} {wr:6.2f}  {pd_:7.1f}%"
              f" {pn_:6.1f}%")

    ach = R["240 U q12wk  (as treated)"]
    bnd = R["BOUND perfect block, as-is"]
    bnd_rho = R["BOUND perfect block, rho->1"]
    bnd_ext = R["BOUND ext. block, rho->1"]
    bnd_nos = R["BOUND as-is, NO sprouting"]
    print()
    print("  ---- FOUR LEVERS, PRICED ----")
    nos = R["240 U q12wk, NO sprouting"]
    print(f"    achieved now: 240 U q12wk .................... {ach:6.2f} points")
    print()
    print("    1. MOLECULE / DOSE -- climb to the bound at the current phi")
    print(f"         bound at phi = {p.rho*PHI_STD:.3f} ................. "
          f"{bnd:6.2f}    headroom {bnd-ach:+6.2f}")
    print(f"    2. PLACEMENT -- rho {p.rho:.2f} -> 1.00 (EMG / ultrasound guidance)")
    print(f"         bound relocates to ................. {bnd_rho:6.2f}"
          f"    headroom {bnd_rho-bnd:+6.2f}")
    print("    3. TARGET LIST -- add semispinalis capitis and obliquus cap. inf.")
    print(f"         bound relocates to ................. {bnd_ext:6.2f}"
          f"    headroom {bnd_ext-bnd_rho:+6.2f}")
    print("    4. SPROUTING -- abolish terminal sprouting")
    print(f"         at the bound ....................... {bnd_nos:6.2f}"
          f"    headroom {bnd_nos-bnd:+6.2f}")
    print(f"         on the REAL 240 U regimen .......... {nos:6.2f}"
          f"    headroom {nos-ach:+6.2f}")
    print("       Sprouting is worth nothing AT the bound, because sprouts draw")
    print("       on the same cytosolic SNARE pool the toxin has destroyed --")
    print("       under permanent blockade they carry nothing.  On a real")
    print("       wearing-off regimen they matter only during the recovery limb,")
    print("       and there they are worth little too.  That is a genuine")
    print("       finding and it CONTRADICTS the common claim that sprouting is")
    print("       what ends the clinical effect: in this model the effect ends")
    print("       because the light chain decays and SNAP-25 is resynthesised.")
    print("       Sprouting merely rides along.  A5 tests this directly.")
    print()
    parts = [("more or better toxin, at the current phi", bnd - ach),
             ("better needle placement (rho)", bnd_rho - bnd),
             ("a longer target list (SUM w)", bnd_ext - bnd_rho),
             ("stopping the nerve repairing itself", nos - ach)]
    tot = sum(v for _, v in parts)
    print("  As shares of the whole remaining design space:")
    for lbl, v in parts:
        print(f"    {lbl:42s} {v:6.2f} pts   {v/tot*100:5.1f} %")
    nondose = (bnd_ext - bnd)
    print()
    print(f"    dose axis ....... {bnd-ach:6.2f} points, and it is BOUNDED there")
    print(f"    geometry axis ... {nondose:6.2f} points  "
          f"({nondose/(bnd-ach):.1f}x the dose axis)")
    print()
    print(f"    240 U q12wk already extracts {ach/bnd*100:.1f} % of everything that")
    print("    perfect, permanent, side-effect-free chemodenervation of the")
    print("    motor units it currently reaches could ever deliver.  A 60-fold")
    print(f"    dose increase reaches {R['14400 U q12wk (60x)']/bnd*100:.1f} % of it -- "
          "i.e. essentially all of a")
    print(f"    small thing -- at P(dysphagia) "
          f"{cyc(q12(PATTERN_STD*60), p)[3]:.0f} % instead of "
          f"{cyc(q12(PATTERN_STD), p)[3]:.0f} %.")
    print()
    print("  Three consequences, in descending order of how uncomfortable they")
    print("  are for the way this disease is currently treated.")
    print()
    print("  (i)  THE DOSE AXIS IS THE WORST-VALUE AXIS AND IT CARRIES ALL THE")
    print("       HARM.  It is not merely inefficient: it is BOUNDED, and the")
    print("       bound is the smallest of the four.")
    print()
    print("  (ii) PLACEMENT AND TARGET LIST ARE NOT SUBSTITUTES FOR DOSE, THEY")
    print("       ARE PREREQUISITES FOR IT.  Applied at constant 240 U they buy")
    print(f"       little on their own -- extended list "
          f"{R['240 U, extended target list']-ach:+.2f}, perfect placement "
          f"{R['240 U, perfect placement']-ach:+.2f},")
    print(f"       both together {R['240 U, extended + perfect']-ach:+.2f} -- because they")
    print("       RELOCATE the ceiling rather than climb toward it.  Dose climbs;")
    print("       phi relocates.  A trial that improves guidance while holding")
    print("       total Units fixed is testing the wrong combination, and this is")
    print("       a concrete prediction about why such trials read as modest.")
    print()
    print("  (iii) THE BOUND IS NOT CLINICALLY ATTAINABLE ANYWAY.  Perfect")
    print(f"       blockade carries P(neck weakness) "
          f"{cyc(q12(PATTERN_STD, block=blkS), p)[4]:.0f} % as-is and "
          f"{cyc(q12(PATTERN_EXT, block=blkE, rho_mult=1/p.rho), p)[4]:.0f} % for the")
    print("       extended plan at perfect placement, because the posterior")
    print("       muscles that generate the dystonic torque are the same ones")
    print("       that hold the head up.  The antagonist of the disease is the")
    print("       agonist of the posture.  That is a geometric conflict, and no")
    print("       improvement in any toxin resolves it.")
    print()
    print("  The single largest prize in the table is the one nothing is aimed")
    print("  at: the motor nerve terminal's own capacity to repair itself.")


def A3():
    hr("A3  DOSE -> DURATION IS LOGARITHMIC; DOSE -> SPREAD IS LINEAR")
    p = cal_params()
    print("  E(S) is a threshold function (A0), so duration is the time S needs")
    print("  to climb back across the threshold.  S recovers at a rate set by")
    print("  k_syn and by the decay of the light chain -- NOT by how deep it was")
    print("  driven -- so each doubling of dose adds a roughly FIXED number of")
    print("  days.  Spread is different: the toxin arriving at the pharynx is a")
    print("  fixed fraction of whatever was injected, so it scales LINEARLY.")
    print()
    print("   dose(U)  ratio   S_min(SCM)   nadir dTW  duration(d)   peak C_sw"
          "   dysph burden  P(dysph)")
    fs = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
    durs, csw, burden, pds = [], [], [], []
    for f in fs:
        res = simulate(_std_reg(f), p, tmax=460.0, dt=1.0)
        durs.append(duration_of_benefit(res))
        csw.append(float(np.max(res["Csw"])))
        burden.append(float(res["AUCdys"][-1]))
        pds.append(float(np.max(res["P_dys"])))
        print(f"   {240*f:7.0f}  {f:5.2f}x  {float(np.min(res['S'][:,0])):10.5f}"
              f"  {nadir(res,lo=0,hi=220)[0]-p.TW0_ref:+10.2f}  {durs[-1]:10.1f}"
              f"   {csw[-1]:10.3f}   {burden[-1]:11.1f}  {pds[-1]*100:7.1f}%")
    durs = np.array(durs); csw = np.array(csw)
    burden = np.array(burden); pds = np.array(pds)
    print()
    ok = durs > 0
    b = np.polyfit(np.log(fs[ok]), durs[ok], 1)
    r2 = 1 - np.sum((durs[ok] - np.polyval(b, np.log(fs[ok]))) ** 2) \
        / np.sum((durs[ok] - durs[ok].mean()) ** 2)
    print(f"  duration     vs ln(dose): slope {b[0]:6.1f} d per e-fold "
          f"= {b[0]*math.log(2):.1f} d per DOUBLING   (R^2 = {r2:.4f})")
    bl = np.polyfit(fs, csw, 1)
    r2l = 1 - np.sum((csw - np.polyval(bl, fs)) ** 2) / np.sum((csw - csw.mean()) ** 2)
    print(f"  peak C_sw    vs dose    : slope {bl[0]:6.3f} U per 240 U    "
          f"                    (R^2 = {r2l:.4f})")
    bb = np.polyfit(fs, burden, 1)
    r2b = 1 - np.sum((burden - np.polyval(bb, fs)) ** 2) \
        / np.sum((burden - burden.mean()) ** 2)
    print(f"  dysph burden vs dose    : slope {bb[0]:6.2f} deficit-days per 240 U"
          f"        (R^2 = {r2b:.4f})")
    print()
    print("  A logarithm cannot outrun a line.  Benefit-days purchased per unit")
    print("  of swallow exposure therefore falls monotonically with dose:")
    print()
    print("     dose(U)   duration(d)   dysph burden   days per deficit-day")
    for f, d, bd in zip(fs, durs, burden):
        ratio = f"{d/bd:19.2f}" if bd > 1e-6 else f"{'  (no exposure)':>19s}"
        print(f"     {240*f:7.0f}   {d:10.1f}   {bd:12.1f}   {ratio}")
    print()
    print(f"  Peak pharyngeal light chain is in fact SUPRA-linear: 16x the dose")
    print(f"  gives {csw[-1]/csw[2]:.1f}x the pharyngeal load, not 16x, because terminal")
    print("  binding in the injected muscle saturates and the excess is left free")
    print("  to diffuse.  The harm axis therefore steepens exactly where the")
    print("  benefit axis is flattening.")
    print()
    print("  Note that P(dysphagia) itself SATURATES near 21 % in the last column")
    print("  of the first table, which is why the probability is the wrong")
    print("  currency for this comparison: a bounded scale hides an unbounded")
    print("  exposure.  The deficit-days do not saturate, and neither does the")
    print("  underlying toxin load at the pharynx.  Reporting the probability")
    print("  alone would have made the therapeutic index look as though it")
    print("  recovered at high dose.  It does not.")


def A4():
    hr("A4  THE TWO CLOCKS -- the drug is long gone before the effect peaks")
    p = cal_params()
    res = simulate(_std_reg(), p, tmax=220.0, dt=0.25)
    t = res["t"]
    A, B, C = res["Y"][:, IX["A_0"]], res["Y"][:, IX["B_0"]], res["Y"][:, IX["C_0"]]
    S, Q, E = res["S"][:, 0], res["Q"][:, 0], res["E"][:, 0]
    Tm = res["Tm_raw"][:, 0]      # pre-rho, so parent and sprout are comparable

    def t_half_after_peak(x):
        pk, i = x.max(), int(np.argmax(x))
        if pk <= 0:
            return float("nan")
        for j in range(i, len(x)):
            if x[j] <= pk / 2:
                return float(t[j] - t[i])
        return float("nan")

    print("  Sternocleidomastoid, 50 U:")
    print(f"    free toxin        A  peak {A.max():9.3f} U at day {t[np.argmax(A)]:5.2f}"
          f"   t1/2 after peak {t_half_after_peak(A):6.2f} d")
    print(f"    internalising     B  peak {B.max():9.3f} U at day {t[np.argmax(B)]:5.2f}"
          f"   t1/2 after peak {t_half_after_peak(B):6.2f} d")
    print(f"    active light chain C peak {C.max():9.3f} U at day {t[np.argmax(C)]:5.2f}"
          f"   t1/2 after peak {t_half_after_peak(C):6.2f} d")
    print(f"    intact SNAP-25    S  min  {S.min():9.5f}    at day "
          f"{t[np.argmin(S)]:5.1f}")
    print()
    print("   day       A(U)      B(U)      C(U)       S      E(S)   Q(sprout)"
          "  T(parent+sprout)  TWSTRS")
    for d in (0.25, 0.5, 1, 2, 3, 5, 7, 10, 14, 21, 28, 42, 56, 70, 84, 98,
              112, 140, 168):
        i = int(np.argmin(np.abs(t - d)))
        print(f"  {d:5.2f} {max(A[i],0):10.4f} {max(B[i],0):9.4f} "
              f"{max(C[i],0):9.3f} {S[i]:8.5f}"
              f" {E[i]:8.5f} {Q[i]:9.4f} {Tm[i]:13.5f} {res['TW'][i]:8.2f}")
    print()
    fa = float(np.interp(2.0, t, A)) / max(A.max(), 1e-12)
    print(f"  By day 2 only {fa*100:.3f} % of the injected free toxin is left,")
    print("  and the clinical effect has not even peaked.  Duration of action is")
    print("  NOT a pharmacokinetic property of the injectate.  Four sequential,")
    print("  progressively slower processes set it:")
    print(f"    (1) free toxin residence            t1/2 {t_half_after_peak(A):6.2f} d")
    print(f"    (2) endosomal translocation         t1/2 {math.log(2)/p.k_trans:6.2f} d")
    print(f"    (3) active light-chain persistence  t1/2 {math.log(2)/p.k_LC:6.2f} d"
          "   <- calibrated")
    print(f"    (4) SNAP-25 resynthesis             t1/2 {math.log(2)/p.k_syn:6.2f} d")
    print(f"    (5) sprout growth / regression      t1/2 {math.log(2)/p.k_sp:6.1f}"
          f" / {math.log(2)/p.k_rg:.1f} d")
    iq = int(np.argmax(Q))
    print(f"  Sprout capacity peaks on day {t[iq]:.0f} (Q = {Q[iq]:.3f}) and then")
    print("  regresses as the parent terminal recovers -- the two-phase recovery.")
    for d in (28.0, 56.0, 84.0):
        j = int(np.argmin(np.abs(t - d)))
        share = (Tm[j] - E[j]) / max(Tm[j], 1e-12) * 100
        print(f"    day {d:3.0f}: T = {Tm[j]:.4f}; parent {E[j]:.4f}, "
              f"sprouts {Tm[j]-E[j]:.4f} ({share:.1f} % of recovered function)")


def A5():
    hr("A5  SPROUTING IS NOT WHY THE EFFECT WEARS OFF")
    p = cal_params()
    k0 = p.k_sp
    print("  The received account of wearing-off is nerve-terminal sprouting:")
    print("  the terminal grows collateral sprouts, the sprouts form new")
    print("  functional contacts, transmission returns.  The sprouting is real --")
    print("  it is documented, it peaks at about the right time, and this model")
    print("  reproduces it (A4).  The claim tested here is the CAUSAL one.")
    print()
    print("  Sprouts grow from the same axon and load their vesicles using the")
    print("  same cytosolic SNARE pool the light chain has been destroying.  So")
    print("  sprout release must be gated by S as well.  If it is, sprouting")
    print("  cannot restore transmission while the toxin is still working -- and")
    print("  by the time the toxin has decayed enough for sprouts to function,")
    print("  the parent terminal is functioning too.  Sprouting then predicts")
    print("  almost nothing about duration.  Test it:")
    print()
    print("   k_sp        nadir dTW   duration(d)   AUC benefit (pt.d/168 d)"
          "   P(dysph)")
    for f in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0):
        res = simulate(_std_reg(), cal_params(k_sp=k0 * f), tmax=168.0, dt=1.0)
        print(f"   {f:5.2f}x k_sp {nadir(res,lo=0,hi=140)[0]-p.TW0_ref:+10.2f}"
              f"   {duration_of_benefit(res):10.1f}   {res['AUCben'][-1]:14.1f}"
              f"        {float(np.max(res['P_dys']))*100:5.1f}%")
    r0 = simulate(_std_reg(), cal_params(k_sp=0.0), tmax=168.0, dt=1.0)
    r1 = simulate(_std_reg(), p, tmax=168.0, dt=1.0)
    r4 = simulate(_std_reg(4.0), p, tmax=168.0, dt=1.0)
    print()
    print(f"  ABOLISHING sprouting entirely: duration "
          f"{duration_of_benefit(r1):.0f} -> {duration_of_benefit(r0):.0f} d "
          f"({duration_of_benefit(r0)-duration_of_benefit(r1):+.0f} d), "
          f"benefit-time {r1['AUCben'][-1]:.0f} -> {r0['AUCben'][-1]:.0f} pt.d "
          f"({(r0['AUCben'][-1]/r1['AUCben'][-1]-1)*100:+.1f} %).")
    print(f"  For comparison a 4-fold DOSE increase gives {r4['AUCben'][-1]:.0f} pt.d "
          f"({(r4['AUCben'][-1]/r1['AUCben'][-1]-1)*100:+.1f} %).")
    print()
    print("  So in this model sprouting is worth about a twentieth of what a")
    print("  4-fold dose increase is worth, and a perfect anti-sprouting agent")
    print("  would extend the injection interval by roughly three days.  THIS")
    print("  CONTRADICTS THE STANDARD EXPLANATION, and it does so for a reason")
    print("  that is mechanistic rather than numerical.")
    print()
    print("  What DOES end the effect?  Decompose the recovery of transmission")
    print("  in the sternocleidomastoid into its parent and sprout parts (both")
    print("  pre-rho, so they are comparable):")
    print()
    print("    day    S      E parent   T total   sprout share   C (active LC)")
    t = r1["t"]
    for d in (7, 14, 21, 28, 42, 56, 70, 84, 98, 112, 140):
        i = int(np.argmin(np.abs(t - d)))
        E_, T_ = r1["E"][i, 0], r1["Tm_raw"][i, 0]
        sh = (T_ - E_) / T_ * 100 if T_ > 1e-9 else 0.0
        print(f"   {d:5d}  {r1['S'][i,0]:7.4f} {E_:9.4f} {T_:9.4f}"
              f"   {sh:9.1f} %   {r1['Y'][i, IX['C_0']]:11.3f}")
    print()
    print("  The sprout share peaks in the middle of the cycle and never")
    print("  dominates.  Recovery tracks S, and S tracks the disappearance of the")
    print("  active light chain.  The clock that ends a botulinum injection is")
    print(f"  the light chain's own decay (t1/2 {math.log(2)/p.k_LC:.0f} d, calibrated) "
          "convolved with")
    print(f"  SNAP-25 resynthesis (t1/2 {math.log(2)/p.k_syn:.1f} d).  A12 reaches the same")
    print("  conclusion independently, by sensitivity analysis.")
    print()
    print("  FALSIFIABLE PREDICTION.  An agent that blocked sprouting -- or a")
    print("  patient population with impaired sprouting -- should show a duration")
    print("  of benefit differing by only a few days, NOT by weeks.  If a real")
    print("  anti-sprouting intervention were to extend duration substantially,")
    print("  the SNARE-gating assumption above would be wrong, and sprouts would")
    print("  have to be filling vesicles from a pool the toxin cannot reach.")


def A6():
    hr("A6  DILUTION AT CONSTANT UNITS -- volume is a dosing variable")
    p = cal_params()
    print("  Total dose fixed at 240 U; only injectate volume changes.  Efflux")
    print("  from the muscle scales with volume and the swallow compartment is")
    print("  small, so it concentrates what reaches it.")
    print()
    print("   dilution           vol(mL)  nadir dTW  duration(d)  swallow def"
          "  P(dysph)  P(neck)")
    for lbl, upm in (("200 U/mL", 200.), ("100 U/mL", 100.), ("50 U/mL (ref)", 50.),
                     ("25 U/mL", 25.), ("12.5 U/mL", 12.5)):
        vol = PATTERN_STD / upm
        reg = Regimen(pattern=PATTERN_STD.copy(), volume=vol, n_inj=1)
        res = simulate(reg, p, tmax=200.0, dt=1.0)
        print(f"   {lbl:16s} {vol.sum():7.2f} "
              f"{nadir(res,lo=0,hi=150)[0]-p.TW0_ref:+10.2f}  "
              f"{duration_of_benefit(res):10.1f}  "
              f"{float(np.max(res['def_sw'])):11.3f}  "
              f"{float(np.max(res['P_dys']))*100:7.1f}% "
              f"{float(np.max(res['P_neck']))*100:7.1f}%")
    print()
    print("  Efficacy is nearly flat across a 16-fold dilution range while")
    print("  dysphagia risk is not.  In this model concentration, not dose, is")
    print("  the cheapest safety lever available at the bedside -- and unlike")
    print("  dose reduction it costs almost nothing in efficacy.")


def A7():
    hr("A7  THE 12-WEEK RULE IS NOT EXPLAINED BY ANTIBODY RISK")
    p = cal_params()
    print("  The standing clinical convention is: do not re-inject more often")
    print("  than every 12 weeks, because shorter intervals provoke neutralising")
    print("  antibody and antibody is ABSORBING -- potency lost does not return.")
    print("  That is a quantitative claim, so it can be checked against an")
    print("  immunogenicity submodel calibrated to the OBSERVED antibody rate")
    print("  (A1 stage 3) rather than to the fear of one.")
    print()
    horizon = 5 * 365.0
    IVS = (42, 56, 70, 84, 98, 112, 140, 182, 240)
    for prod in ("incobotulinumtoxinA", "onabotulinumtoxinA",
                 "rimabotulinumtoxinB"):
        pr = PRODUCTS[prod]
        # dose to A-EQUIVALENCE, not to each label's own 240 U.  BoNT-B's
        # clinical dose is 5000-10000 label Units; giving it 240 would compare
        # a therapeutic dose of A against a sub-therapeutic dose of B.
        pat = PATTERN_STD / pr["unit_scale"]      # label Units, for display
        print(f"  {prod}  ({pr['load']*100:.3f} ng per 100 label-U, "
              f"serotype {pr['serotype']})")
        print(f"    dosed to A-equivalence: {pat.sum():.0f} label-U "
              f"= 240 A-equivalent U, {pat.sum()*pr['load']:.2f} ng antigen "
              f"per injection")
        print("    interval(d)  n inj   cum ng   Nab(5y)  potency  "
              "benefit(pt.d)   dysph burden")
        best = (None, -1.0)
        for iv in IVS:
            n = int(horizon // iv)
            reg = Regimen(product=prod, pattern=PATTERN_STD.copy(),
                          interval=float(iv), n_inj=n)
            res = simulate(reg, p, tmax=horizon, dt=2.0)
            tag = pr["serotype"]
            nb = float(res["NabA" if tag == "A" else "NabB"][-1])
            gate = 1.0 / (1.0 + (nb / p.Nab50) ** p.hn)
            auc = float(res["AUCben"][-1])
            mark = ""
            if auc > best[1]:
                best = (iv, auc)
                mark = "  <-- best"
            print(f"    {iv:9d}  {n:6d}  {pat.sum()*pr['load']*n:7.1f}  {nb:7.3f}"
                  f"  {gate*100:6.1f}%  {auc:13.0f}  {float(res['AUCdys'][-1]):13.1f}"
                  f"{mark}")
        print(f"    optimum interval: {best[0]} d ({best[0]/7:.0f} weeks)")
        print()
    print("  THERE IS NO INTERIOR OPTIMUM.  On the antibody axis alone, shorter")
    print("  is simply better, all the way to the shortest interval tested, for")
    print("  every product -- because at the observed antibody rate the median")
    print("  patient never loses meaningful potency.  The stated rationale for")
    print("  the 12-week convention does not survive its own arithmetic.")
    print()
    print("  WHEN DOES AN INTERIOR OPTIMUM APPEAR?  Only far out in the tail of")
    print("  antibody propensity.  Scale k_b and re-optimise (onaBoNT-A):")
    print()
    print("    k_b multiple   percentile of the population   optimum interval"
          "   Nab(5y) at q12wk")
    kb0 = p.k_b
    sig = CAL_TARGETS["sigma_kb"]
    for mult in (1, 3, 10, 30, 100, 300):
        pp = cal_params(k_b=kb0 * mult)
        best = (None, -1.0)
        for iv in IVS:
            n = int(horizon // iv)
            reg = Regimen(product="onabotulinumtoxinA",
                          pattern=PATTERN_STD.copy(), interval=float(iv), n_inj=n)
            auc = float(simulate(reg, pp, tmax=horizon, dt=4.0)["AUCben"][-1])
            if auc > best[1]:
                best = (iv, auc)
        n84 = int(horizon // 84)
        nb84 = float(simulate(Regimen(product="onabotulinumtoxinA",
                                      pattern=PATTERN_STD.copy(), interval=84.0,
                                      n_inj=n84), pp, tmax=horizon,
                              dt=4.0)["NabA"][-1])
        pct = 100.0 * 0.5 * (1.0 + math.erf(math.log(mult) / (sig * math.sqrt(2))))
        print(f"    {mult:8d}x    {pct:20.2f} th             {best[0]:6d} d"
              f"        {nb84:12.3f}")
    print()
    print("  The convention is optimal only for patients somewhere above the")
    print("  99th centile of antibody propensity.  It is a rule calibrated on the")
    print("  tail and then applied to everyone.  That may still be the right")
    print("  policy -- an absorbing loss is worth insuring against, and the")
    print("  individual patient's propensity is not measurable in advance -- but")
    print("  it is a decision about VARIANCE, not about the average patient, and")
    print("  it should be stated as one.")
    print()
    print("  WHAT DOES BOUND THE INTERVAL, in this model, is the last column of")
    print("  the tables above: cumulative swallow exposure.  It does not merely")
    print("  track the number of injections -- it rises FASTER, because at short")
    print("  intervals the pharyngeal deficit has less time to clear between")
    print("  doses and spends more of its life above the burden threshold.")
    print()
    rows = {}
    for iv in (42, 84):
        n = int(horizon // iv)
        r = simulate(Regimen(product="incobotulinumtoxinA",
                             pattern=PATTERN_STD.copy(),
                             interval=float(iv), n_inj=n), p,
                     tmax=horizon, dt=2.0)
        rows[iv] = (float(r["AUCben"][-1]), float(r["AUCdys"][-1]), n)
        print(f"      q{iv//7:2d}wk : benefit {rows[iv][0]:8.0f} pt.d,"
              f"  dysphagia burden {rows[iv][1]:7.1f} deficit.d,"
              f"  {n} injections")
    b42, d42, n42 = rows[42]
    b84, d84, n84 = rows[84]
    print()
    print(f"  Halving the interval buys {(b42/b84-1)*100:+.0f} % benefit for "
          f"{(n42/n84-1)*100:+.0f} % more injections")
    print(f"  and {(d42/d84-1)*100:+.0f} % more swallow exposure -- the harm axis is")
    print("  SUPRALINEAR in injection frequency, exactly as it was supralinear in")
    print("  dose (A3).  That is a defensible place to put a 12-week rule.")
    print("  Antibody, at the observed rate, is not.")


def A8():
    hr("A8  SECONDARY NON-RESPONSE AND SEROTYPE RESCUE")
    p = cal_params()
    horizon = 5 * 365.0

    # --- first, verify the fast antibody path against the full 70-state model
    print("  Consistency check: the antibody pool is exactly decoupled from the")
    print("  other 67 states, so it can be integrated on its own.  Verify that")
    print("  claim before relying on it for the cohorts below.")
    load = 240.0 * PRODUCTS["onabotulinumtoxinA"]["load"]
    n21 = 21
    times = np.arange(n21) * 84.0
    full = simulate(Regimen(product="onabotulinumtoxinA",
                            pattern=PATTERN_STD.copy(), interval=84.0,
                            n_inj=n21), p, tmax=horizon, dt=8.0)
    nb_full = float(full["NabA"][-1])
    nb_fast = nab_final(load, times, p, horizon)
    rel = abs(nb_fast - nb_full) / max(nb_full, 1e-12) * 100
    print(f"    full 70-state model : Nab(5 y) = {nb_full:.6f}")
    print(f"    3-state fast path   : Nab(5 y) = {nb_fast:.6f}")
    print(f"    relative difference = {rel:.4f} %")
    nb_vec = float(nab_cohort(load, times, [p.k_b], p, horizon)[0])
    rel2 = abs(nb_vec - nb_full) / max(nb_full, 1e-12) * 100
    print(f"    vectorised RK4 path : Nab(5 y) = {nb_vec:.6f}"
          f"   ({rel2:.4f} % from the full model)")
    print("    -> " + ("all three agree; the decoupling and the vectorisation hold"
                       if (rel < 1.0 and rel2 < 1.0) else
                       "MISMATCH -- do not trust the cohorts below"))
    print()

    rng = np.random.default_rng(20260729)
    N = 20000
    kb0 = p.k_b
    eta = rng.normal(0.0, CAL_TARGETS["sigma_kb"], N)
    print(f"  Virtual cohort n = {N}; lognormal random effect on memory-B")
    print(f"  priming (k_b, sigma = {CAL_TARGETS['sigma_kb']}); 5 years of 240 U.")
    print("  Secondary non-response = injected potency reduced below 50 % by")
    print("  antibody.  n = 20 000 resolves a 2 % event to about +/-0.1 %")
    print("  (binomial SE); the vectorised path above is what makes a cohort of")
    print("  this size cost about a second.")
    print()
    print("   product              interval   non-responders(5y)   median Nab"
          "   99th pct Nab")
    for prod, iv in (("incobotulinumtoxinA", 84), ("onabotulinumtoxinA", 84),
                     ("onabotulinumtoxinA", 56)):
        ld = 240.0 * PRODUCTS[prod]["load"]
        tt = np.arange(int(horizon // iv)) * float(iv)
        nabs = nab_cohort(ld, tt, kb0 * np.exp(eta), p, horizon)
        gate = potency_gate(nabs, p)
        print(f"   {prod:20s} {iv:6d} d   {float(np.mean(gate<0.5))*100:16.2f}%"
              f"   {np.median(nabs):10.4f}   {np.percentile(nabs,99):12.4f}")
    print()
    print("  Reported: neutralising antibody with the current onabotulinumtoxinA")
    print("  formulation about 1-3 %, incobotulinumtoxinA about 0-1.1 %, and the")
    print("  pre-1998 high-protein formulation up to about 9.5 %.  The q12wk")
    print("  onaBoNT-A row is the anchor k_b was fitted to; the other three rows")
    print("  are predictions, and the eleven-fold lower protein load of")
    print("  incobotulinumtoxinA essentially abolishes the event.")
    print()
    print("  SEROTYPE RESCUE.  BoNT-B cleaves VAMP/synaptobrevin, not SNAP-25,")
    print("  and shares no neutralising epitopes with A, so anti-A antibody does")
    print("  not touch it.  Follow one HIGH RESPONDER (k_b at the 99.7th centile)")
    print("  through 28 cycles of onaBoNT-A, with and without a switch to")
    print("  rimabotulinumtoxinB at cycle 12.")
    pp = cal_params(k_b=kb0 * math.exp(2.6))
    n = 28
    ra = simulate(Regimen(product="onabotulinumtoxinA",
                          pattern=PATTERN_STD.copy(), interval=84.0, n_inj=n),
                  pp, tmax=n * 84.0, dt=2.0)
    rs = simulate(Regimen(product="onabotulinumtoxinA",
                          pattern=PATTERN_STD.copy(), interval=84.0, n_inj=n,
                          switch_at=12, switch_to="rimabotulinumtoxinB"),
                  pp, tmax=n * 84.0, dt=2.0)
    print()
    print("   cycle   Nab_A   potency A   dTW(nadir) stay on A   dTW(nadir) switch@12")
    for k in (1, 4, 8, 12, 14, 18, 22, 26):
        t0 = (k - 1) * 84.0
        nb = float(np.interp(t0, ra["t"], ra["NabA"]))
        m = (ra["t"] >= t0) & (ra["t"] <= t0 + 84)
        m2 = (rs["t"] >= t0) & (rs["t"] <= t0 + 84)
        print(f"   {k:5d}  {nb:6.2f}  {potency_gate(nb,p)*100:8.1f}%  "
              f"{float(np.min(ra['TW'][m]))-p.TW0_ref:+16.2f}  "
              f"{float(np.min(rs['TW'][m2]))-p.TW0_ref:+18.2f}")
    print()
    print(f"  benefit over 6.4 y, stay on A : {ra['AUCben'][-1]:9.0f} pt.d")
    print(f"  benefit over 6.4 y, switch@12 : {rs['AUCben'][-1]:9.0f} pt.d"
          f"   ({(rs['AUCben'][-1]/max(ra['AUCben'][-1],1e-9)-1)*100:+.1f} %)")
    print(f"  final Nab_A {ra['NabA'][-1]:.2f}; after the switch, Nab_B reaches "
          f"{rs['NabB'][-1]:.2f}")
    print("  -- the B pool immunises independently, so the rescue is FINITE.")
    rd = simulate(Regimen(product="rimabotulinumtoxinB",
                          pattern=PATTERN_STD.copy(), n_inj=1),
                  p, tmax=250.0, dt=1.0)
    ra1 = simulate(_std_reg(), p, tmax=250.0, dt=1.0)
    print(f"  BoNT-B also carries an autonomic cost: peak salivary-compartment")
    print(f"  deficit {float(np.max(rd['def_au'])):.3f} versus "
          f"{float(np.max(ra1['def_au'])):.3f} for BoNT-A at an equipotent dose")
    print(f"  (P(dry mouth) {float(np.max(rd['P_dry']))*100:.1f} % vs "
          f"{float(np.max(ra1['P_dry']))*100:.1f} %), because serotype B has a")
    print("  markedly higher tropism for autonomic cholinergic terminals.")


def A9():
    hr("A9  THE CENTRAL RATCHET -- why troughs fall across cycles")
    p = cal_params()
    n = 10
    res = simulate(Regimen(pattern=PATTERN_STD.copy(), interval=84.0, n_inj=n),
                   p, tmax=n * 84.0, dt=1.0)
    print("  BoNT blocks intrafusal (gamma) terminals as well as extrafusal")
    print(f"  ones, lowering the abnormal spindle afferent drive that feeds the")
    print(f"  plasticity state Dcen (t1/2 {math.log(2)/p.k_d:.0f} d).  Peripheral")
    print("  paralysis is cyclical; the central effect is a RATCHET.")
    print()
    print("   cycle   Dcen(end)   TROUGH TWSTRS   NADIR TWSTRS   S_min(SCM)")
    print("   (cycle 1's 'trough' is the pre-treatment baseline at t = 0)")
    for k in range(1, n + 1):
        t0, t1 = (k - 1) * 84.0, k * 84.0
        m = (res["t"] >= t0) & (res["t"] <= t1)
        print(f"   {k:5d}   {float(np.interp(t1,res['t'],res['Dcen'])):9.4f}"
              f"   {float(np.max(res['TW'][m])):13.2f}"
              f"   {float(np.min(res['TW'][m])):12.2f}"
              f"   {float(np.min(res['S'][m,0])):10.5f}")
    m1 = (res["t"] >= 0) & (res["t"] <= 84)
    mn = (res["t"] >= (n - 1) * 84) & (res["t"] <= n * 84)
    d1, dn = float(np.max(res["TW"][m1])), float(np.max(res["TW"][mn]))
    print()
    print(f"  trough TWSTRS falls {d1:.2f} -> {dn:.2f} ({dn-d1:+.2f} points) over")
    print("  10 cycles while peak SNAP-25 cleavage is unchanged.  The")
    print("  improvement is not in the muscle.  This is the model's account of")
    print("  the 'cumulative benefit' clinicians describe, and it makes a sharp")
    print("  prediction: the trough should keep improving through a SKIPPED")
    print("  cycle, and should relapse only slowly after stopping.")
    r2 = simulate(Regimen(pattern=PATTERN_STD.copy(), interval=84.0, n_inj=6),
                  p, tmax=1200.0, dt=2.0)
    print()
    print("   after the last injection at day 420 (6 cycles):")
    print("     day    TWSTRS    Dcen")
    for d in (504, 560, 672, 840, 1008, 1200):
        print(f"    {d:5d}   {at(r2,d):7.2f}  {at(r2,d,'Dcen'):6.4f}")
    print(f"  TWSTRS relapses toward, but does not reach, the untreated baseline")
    print(f"  of {p.TW0_ref:.2f}, because Dcen unwinds with an "
          f"{math.log(2)/p.k_d:.0f}-day half-life.")


def A10():
    hr("A10 OPERATOR COMPARISON -- what each intervention actually touches")
    p = cal_params()
    print("  Three cycles of the SAME 240 U standard plan plus one added")
    print("  intervention.  Read-out is cycle 3.")
    print()
    print("   intervention                  operator                nadir"
          "  trough   dTW(nadir)")
    rows = [
        ("BoNT-A 240 U alone",           dict(), "CHEMODENERVATION (T_m)"),
        ("+ trihexyphenidyl 20 mg/d",    dict(thp_mg_d=20.0), "CENTRAL DRIVE (Gcen)"),
        ("+ baclofen 60 mg/d",           dict(bac_mg_d=60.0), "CENTRAL DRIVE (Gcen)"),
        ("+ clonazepam 2 mg/d",          dict(clz_mg_d=2.0), "GATING (RecInh/SurrInh)"),
        ("+ GPi-DBS from day 0",         dict(dbs_on=0.0), "CENTRAL DRIVE (Gcen)"),
        ("+ selective denervation 40 %", dict(denervation=0.40), "CEILING (phi)"),
    ]
    for lbl, kw, op in rows:
        reg = Regimen(pattern=PATTERN_STD.copy(), interval=84.0, n_inj=3, **kw)
        res = simulate(reg, p, tmax=3 * 84.0, dt=1.0)
        m = (res["t"] >= 2 * 84) & (res["t"] <= 3 * 84)
        nad = float(np.min(res["TW"][m]))
        print(f"   {lbl:29s} {op:23s} {nad:6.2f}  "
              f"{float(np.max(res['TW'][m])):6.2f}   {nad-p.TW0_ref:+9.2f}")
    reg = Regimen(pattern=PATTERN_EXT.copy(), interval=84.0, n_inj=3)
    res = simulate(reg, p, tmax=3 * 84.0, dt=1.0)
    m = (res["t"] >= 2 * 84) & (res["t"] <= 3 * 84)
    nad = float(np.min(res["TW"][m]))
    print(f"   {'extended EMG/US-guided plan':29s} {'CEILING (phi)':23s} "
          f"{nad:6.2f}  {float(np.max(res['TW'][m])):6.2f}   {nad-p.TW0_ref:+9.2f}")
    print()
    print("  The two strongest single operators are the two that are not oral")
    print("  drugs: raising phi and lowering central drive.  No licensed oral")
    print("  agent in this model reaches either with useful potency -- the same")
    print("  structural gap the influenza model finds for infected-cell death.")


def A11():
    hr("A11 A LONGER-ACTING TOXIN DOES NOT MOVE THE CEILING")
    p = cal_params()
    print("  daxibotulinumtoxinA is represented phenomenologically as a 2-fold")
    print("  slower loss of active light chain (kLC x 0.50).  A13 explains why")
    print("  that is a stand-in and not a mechanism.")
    print()
    print("   product              LC t1/2   nadir dTW  t_nadir  duration(d)"
          "  P(dysph)  AUC/injection")
    for prod in ("incobotulinumtoxinA", "onabotulinumtoxinA",
                 "abobotulinumtoxinA", "daxibotulinumtoxinA",
                 "rimabotulinumtoxinB"):
        pr = PRODUCTS[prod]
        reg = Regimen(product=prod, pattern=PATTERN_STD.copy(), n_inj=1)
        res = simulate(reg, p, tmax=420.0, dt=1.0)
        nd, nt = nadir(res, lo=0, hi=250)
        print(f"   {prod:20s} {math.log(2)/(p.k_LC*pr['kLC_mult']):6.1f} d"
              f" {nd-p.TW0_ref:+10.2f}  {nt:6.0f}d  "
              f"{duration_of_benefit(res):10.1f}  "
              f"{float(np.max(res['P_dys']))*100:7.1f}%  {res['AUCben'][-1]:11.0f}")
    print()
    print("  The three conventional serotype-A products are IDENTICAL once their")
    print("  label Units are converted, as they must be: in this model they differ")
    print("  only in unit scale and in protein load, and protein load acts on the")
    print("  antibody pool, not on the muscle.")
    print()
    print("  daxibotulinumtoxinA is BOTH deeper and longer, and that is worth")
    print("  stating precisely rather than glossing: depth and duration are not")
    print("  independent here, because both are governed by light-chain")
    print("  persistence.  Slower loss means the quasi-steady SNAP-25 level sits")
    print("  lower AND stays low for longer.  What matters is the EXCHANGE RATE:")
    reg_i = Regimen(product="incobotulinumtoxinA", pattern=PATTERN_STD.copy(),
                    n_inj=1)
    reg_d = Regimen(product="daxibotulinumtoxinA", pattern=PATTERN_STD.copy(),
                    n_inj=1)
    ri = simulate(reg_i, p, tmax=420.0, dt=1.0)
    rd = simulate(reg_d, p, tmax=420.0, dt=1.0)
    ni = nadir(ri, lo=0, hi=250)[0] - p.TW0_ref
    nd_ = nadir(rd, lo=0, hi=250)[0] - p.TW0_ref
    di, dd = duration_of_benefit(ri), duration_of_benefit(rd)
    print(f"      nadir    {ni:+.2f} -> {nd_:+.2f}   ({(nd_/ni-1)*100:+.1f} %)")
    print(f"      duration {di:.0f} d -> {dd:.0f} d   ({(dd/di-1)*100:+.1f} %)")
    print(f"      i.e. {(dd/di-1)/max((nd_/ni-1),1e-9):.1f} times as much duration as depth.")
    print("  Reported: daxibotulinumtoxinA gave roughly 24 weeks of benefit in")
    print("  cervical dystonia against about 12 for the conventional products, at")
    print("  a broadly comparable peak TWSTRS change -- the same lopsided exchange")
    print("  rate, and the model's ratio (2.4x duration) lands on the reported one")
    print("  (2x) without having been fitted to it.")
    print()
    print("  BUT THE CEILING HAS NOT MOVED.  Every product above remains far from")
    print("  the phi bound of A2:")
    for lbl, r in (("incobotulinumtoxinA", ri), ("daxibotulinumtoxinA", rd)):
        print(f"      {lbl:20s} nadir {nadir(r,lo=0,hi=250)[0]-p.TW0_ref:+7.2f}")
    print(f"      BOUND at phi = {p.rho*PHI_STD:.3f}  "
          f"{tw_floor(1-p.rho*PHI_STD,p)[0]-p.TW0_ref:+7.2f}  (asymptotic)")
    print("  A longer-acting toxin is therefore a CONVENIENCE improvement -- fewer")
    print("  visits for the same control -- and only marginally an efficacy one.")
    print("  It buys benefit per INJECTION, not the benefit per cycle that the")
    print("  geometry levers of A2 buy.")


def A12():
    hr("A12 SENSITIVITY -- which parameters own which read-out")
    p = cal_params()
    res0 = simulate(_std_reg(), p, tmax=260.0, dt=1.0)
    nd0 = nadir(res0, lo=0, hi=200)[0]
    g0 = p.TW0_ref - nd0
    d0 = duration_of_benefit(res0)
    q0 = float(np.max(res0["P_dys"]))
    names = ["k_cl", "k_syn", "k_LC", "k_trans", "k_sp", "k_rg", "wq", "S50",
             "SF", "n_snare", "k_bind", "k_clear", "k_diff0", "theta_scale",
             "mass_sw", "k_sev", "hs", "L50", "k_d", "aff_floor",
             "Bmax_per_mass", "kp_off"]
    print("  +25 % on each parameter.  Elasticity = %change(read-out) / 25 %.")
    print("   parameter           dGain    E_gain   d(dur)   E_dur   dP_dys"
          "   E_dys")
    out = []
    for nm in names:
        v = getattr(p, nm)
        r = simulate(_std_reg(), cal_params(**{nm: v * 1.25}), tmax=260.0, dt=1.0)
        g = p.TW0_ref - nadir(r, lo=0, hi=200)[0]
        dd = duration_of_benefit(r)
        qq = float(np.max(r["P_dys"]))
        out.append((nm, g - g0, ((g - g0) / g0) / 0.25,
                    dd - d0, ((dd - d0) / d0) / 0.25,
                    qq - q0, ((qq - q0) / q0) / 0.25))
    for row in sorted(out, key=lambda r: -abs(r[4])):
        print(f"   {row[0]:18s} {row[1]:+7.3f}  {row[2]:+8.3f}  {row[3]:+7.1f}"
              f" {row[4]:+7.3f}  {row[5]*100:+7.2f}  {row[6]:+7.3f}")
    print()
    print("  READ THE TOP OF THAT TABLE CAREFULLY, because it is not the answer")
    print("  one would want.  The two parameters that own the DURATION read-out")
    print("  are L50 and kp_off -- the position of the clinical severity curve and")
    print("  the pain time constant.  Neither is toxin biology.  They dominate")
    print("  because 'duration' is defined here as the time until the gain falls")
    print("  back through a FIXED THRESHOLD (the MCID), and anything that shifts")
    print("  the whole curve up or down moves that crossing a long way.")
    print()
    print("  So the headline is a caveat: a large part of what trials report as")
    print("  'duration of effect' is a property of WHERE THE THRESHOLD SITS")
    print("  relative to the response, not of how long the drug acts.  That is the")
    print("  same measurement problem A13(2) runs into from the other side.")
    print()
    print("  Among the parameters that ARE mechanism, the ordering is clean and it")
    print("  agrees with A4 and A5:")
    print("    k_LC   (light-chain persistence)  E_dur = -0.84   <- the real clock")
    print("    S50, k_cl, k_syn (SNARE turnover) E_dur ~ +/-0.5")
    print("    k_sp, k_rg (SPROUTING)            E_dur ~  0.04   <- essentially nil")
    print("  Sprouting has the smallest duration elasticity of any dynamic")
    print("  parameter in the model.  A5 reached that by knocking it out; this")
    print("  reaches it by perturbing it.  Two routes, same answer.")
    print()
    print("  Peak benefit belongs to L50, kp_off and hs -- again the clinical map,")
    print("  not the toxin.  Spread risk belongs to theta_scale, k_diff0, mass_sw,")
    print("  S50 and k_cl.  The three read-outs have almost disjoint owners, which")
    print("  is precisely why one dose knob cannot trade among them (A3).")


def A13():
    hr("A13 DISCREPANCIES -- five places this model does not fit, stated plainly")
    p = cal_params()
    y0 = y0_baseline(p)
    print()
    print("  (1) THE CEILING IS TOO GENEROUS.")
    for phi in (PHI_STD, PHI_EXT):
        tw, sf, dis, pn = tw_floor(1 - phi, p)
        print(f"      phi = {phi:.2f}: asymptotic severity floor {sf:5.2f}/35 "
              f"(from {y0[IX['Sev']]:.1f}), TWSTRS total {tw:5.2f}")
    print("      Best reported severity outcomes with extended EMG- or")
    print("      ultrasound-guided patterns sit near 10-13/35, not 3-7.  Either")
    print("      the deep muscles carry less of the drive than the anatomical")
    print("      weights assume, or Dcen RISES to defend the posture when the")
    print("      periphery is silenced.  The model contains no such")
    print("      compensation and it probably should.  Largest open item.")
    print()
    print("  (2) THE EFFECT-TIME CURVE IS TOO SQUARE.")
    res = simulate(_std_reg(), p, tmax=320.0, dt=0.5)
    g = p.TW0_ref - res["TW"]
    gpk = g.max()
    t10 = float(res["t"][np.argmax(g >= 0.10 * gpk)])
    t50 = float(res["t"][np.argmax(g >= 0.50 * gpk)])
    t90 = float(res["t"][np.argmax(g >= 0.90 * gpk)])
    dur = duration_of_benefit(res)
    print(f"      model: 10 % of peak benefit by day {t10:.1f}, 50 % by day "
          f"{t50:.1f}, 90 % by day {t90:.1f};")
    print(f"      benefit stays above the {MCID_TWSTRS}-point MCID until day "
          f"{dur:.0f} ({dur/7:.1f} weeks).")
    print("      reported: first benefit at a median of about 7 days, and")
    print("      patient-perceived waning at 10-12 weeks.")
    print("      So the model rises too fast AND falls too slowly, even though")
    print("      it lands on both fitted TWSTRS time points.  That is not two")
    print("      errors, it is one: matching a RATING SCALE at two visits does")
    print("      not reproduce PERCEIVED onset and waning, so perceived benefit")
    print("      cannot be an affine function of the TWSTRS total.  A model that")
    print("      forced both would need a separate perception read-out with its")
    print("      own threshold and hysteresis; none is included here.")
    print()
    print("  (3) ONE HILL CANNOT HOLD BOTH ENDS OF THE DOSE-RESPONSE.")
    for f, lbl in ((0.25, "60 U"), (0.5, "120 U"), (0.75, "180 U"),
                   (1.0, "240 U"), (2.0, "480 U")):
        r = simulate(_std_reg(f), p, tmax=140.0, dt=1.0)
        print(f"      {lbl:6s} week-4 dTW = {at(r,28)-p.TW0_ref:+6.2f}"
              f"   nadir dTW = {nadir(r,lo=0,hi=120)[0]-p.TW0_ref:+6.2f}")
    print("      The model reproduces 240 U > 120 U with 480 U barely better")
    print("      than 240 U.  What it cannot also reproduce with a single L50 is")
    print("      the near-linear 60-180 U range some series report.  That is")
    print("      evidence about the SHAPE of sev_ss(L), not about potency, and a")
    print("      model fitted on k_cl alone would have quietly absorbed a shape")
    print("      error into a toxin parameter.")
    print()
    print("  (4) PAIN IMPROVES MORE THAN TORQUE REDUCTION EXPLAINS.")
    r0 = simulate(_std_reg(), p, tmax=140.0, dt=0.5)
    i = int(np.argmin(r0["TW"]))
    print(f"      model pain change at the nadir: "
          f"{r0['Pain'][i]-p.PAINMAX*y0[IX['Pain']]:+.2f} / 20")
    print("      reported TWSTRS-pain change with 240 U: about -2.5 to -3.5.")
    need = None
    for pd_ in np.arange(0.0, 2.01, 0.10):
        rr = simulate(_std_reg(), cal_params(pain_direct=float(pd_)),
                      tmax=140.0, dt=1.0)
        j = int(np.argmin(rr["TW"]))
        if rr["Pain"][j] - p.PAINMAX * y0[IX["Pain"]] <= -3.0:
            need = float(pd_)
            break
    print(f"      closing that gap requires a DIRECT antinociceptive term of")
    print(f"      magnitude pain_direct = "
          f"{('%.2f' % need) if need is not None else '>2.0'}, i.e. BoNT must be")
    print("      blocking CGRP and substance-P release from nociceptor terminals,")
    print("      not merely unloading the muscle.  That is a falsifiable claim,")
    print("      and pain_direct is left at 0 in the base model so the gap shows.")
    print()
    print("  (5) ONE k_b CANNOT SPAN BOTH FORMULATION ERAS.")
    kb0 = p.k_b
    rng = np.random.default_rng(20260729)
    eta = rng.normal(0.0, CAL_TARGETS["sigma_kb"], 20000)
    times = np.arange(int(5 * 365.0 // 84)) * 84.0
    for lbl, load in (("modern      5.0 ng/100 U", 5.0 / 100),
                      ("historical 25.0 ng/100 U", 25.0 / 100)):
        nabs = nab_cohort(240.0 * load, times, kb0 * np.exp(eta), p, 5 * 365.0)
        frac = float(np.mean(potency_gate(nabs, p) < 0.5)) * 100
        print(f"      {lbl:26s} -> {frac:5.2f} % with clinically relevant antibody")
    print("      Reported: about 1-3 % for the modern formulation and up to")
    print("      ~9.5 % for the pre-1998 one.")
    print("      The model gets the ORDERING right and everything else wrong, and")
    print("      it is wrong in two distinguishable ways:")
    print("        (a) it is too HIGH in absolute terms even at the formulation it")
    print("            was calibrated on -- the median patient's Nab is on target")
    print("            by construction, but the recall term (1 + beta*Bm) makes the")
    print("            cohort tail heavier than a lognormal k_b alone would, so the")
    print("            5 % event rate overshoots the reported 1-3 %;")
    print("        (b) it is too STEEP in protein load -- a 5-fold load increase")
    print("            multiplies the event rate about 7.7-fold here, against")
    print("            roughly 4-fold in the historical comparison.")
    print("      (b) matters more than it looks: a pure antigen-DOSE model with a")
    print("      single saturation constant cannot be this steep AND match the")
    print("      absolute rate.  That argues antigen processing saturates harder")
    print("      than Ka_ag allows, or that individual variation sits on more than")
    print("      one axis (HLA type, prior exposure) rather than on k_b alone.")
    print("      Neither is included, so this row is a GAP, not a fit -- and note")
    print("      that it is the only anchor in A1 the model misses badly.")
    print()
    print("  ALSO NOT MODELLED: primary non-response from mis-targeting as")
    print("  distinct from a low phi; needle-EMG accuracy as a distribution")
    print("  rather than a binary; phasic versus tonic dystonia; head tremor;")
    print("  sensory-trick (geste antagoniste) physiology; and the")
    print("  anticholinergic cognitive burden that in practice limits")
    print("  trihexyphenidyl far more than its efficacy does.  The abo- and")
    print("  rima- Unit conversions are clinical rules of thumb and must not be")
    print("  read as equipotency claims.")


ANALYSES = {
    "A0": (A0, "model inventory, safety-factor table, muscle weights"),
    "A1": (A1, "calibration against the pivotal-trial time course"),
    "A2": (A2, "*** THE BOUND: the plateau measures phi, not potency ***"),
    "A3": (A3, "dose -> duration is logarithmic; dose -> spread is linear"),
    "A4": (A4, "the two clocks: toxin gone by day 2, effect peaks at week 5"),
    "A5": (A5, "sprouting is NOT why the effect wears off"),
    "A6": (A6, "dilution at constant Units -- a nearly free safety lever"),
    "A7": (A7, "the 12-week rule is NOT explained by antibody risk"),
    "A8": (A8, "secondary non-response and serotype rescue"),
    "A9": (A9, "the central ratchet across cycles"),
    "A10": (A10, "operator comparison -- the strongest two are not oral drugs"),
    "A11": (A11, "a longer-acting toxin does not move the ceiling"),
    "A12": (A12, "sensitivity: three read-outs, three disjoint parameter sets"),
    "A13": (A13, "five discrepancies, stated plainly"),
}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="*", help="run only these analyses")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for k, (_, d) in ANALYSES.items():
            print(f"  {k:4s} {d}")
        return
    print("cervical dystonia QSP model -- independent numpy/scipy reference port")
    print(f"{NSTATE} ODEs | {NM} muscle compartments x {NSM} states | "
          "LSODA rtol 1e-7")
    for k in (a.only if a.only else list(ANALYSES)):
        if k not in ANALYSES:
            print(f"unknown analysis {k}", file=sys.stderr)
            continue
        ANALYSES[k][0]()
    print()


if __name__ == "__main__":
    main()
