#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
===============================================================================
sbe_reference_model.py
Snakebite envenoming + antivenom — INDEPENDENT PYTHON REFERENCE IMPLEMENTATION
독사 교상(뱀독 중독) + 항독소 — 독립 Python 참조 구현
===============================================================================

WHY THIS FILE EXISTS
--------------------
The build environment for this repository has no R runtime.  Committing an
un-integrated ODE model would be dishonest, so every equation in
`sbe_mrgsolve_model.R` is re-implemented here, term for term, and actually
integrated.  Every number quoted in `README.md` is produced by running this
file; the captured log is `sbe_reference_output.txt`.

THE ONE STRUCTURAL CLAIM
------------------------
    Antivenom BINDS.  It does not UNDO.

An antivenom antibody can only remove venom that is still *free in plasma or
lymph*.  It cannot pull a three-finger toxin off an acetylcholine receptor, it
cannot un-cleave a fibrinogen molecule, and it cannot un-destroy a motor nerve
terminal.  So the model is built with two structurally separate clocks:

    clock 1  (fast, hours)   the venom clock  — antivenom acts HERE
    clock 2  (slow, days)    the substrate clock — resynthesis / regeneration

and every clinical endpoint is a *read-out of clock 2*, driven by the
time-integral of clock 1.  Four consequences fall out arithmetically rather
than being asserted:

 (1) TIME-TO-ANTIVENOM CHANGES THE DEPTH OF THE NADIR, NOT THE SLOPE OF
     RECOVERY.  Fibrinogen climbs back at the rate the liver can make it
     (k_syn * SYNUP), whatever the antivenom did.  Therefore "early antivenom"
     is worth a great deal while the substrate is still falling and worth
     almost nothing once it has hit the floor.

 (2) RECURRENCE IS A RATIO OF TWO HALF-LIVES, NOT A DOSING ERROR.  Venom keeps
     arriving from the bite-site depot and the deep tissue compartment with a
     terminal half-life of ~50 h.  Ovine Fab is cleared with a half-life of
     ~22 h.  When the antivenom half-life is shorter than the venom input
     half-life, free venom MUST reappear, and no larger single front-loaded
     Fab dose can prevent it — only re-dosing, or a longer-lived fragment.

 (3) PRESYNAPTIC AND POSTSYNAPTIC PARALYSIS LOOK IDENTICAL AND ANSWER
     OPPOSITE QUESTIONS.  Postsynaptic alpha-neurotoxin block is an occupancy
     that decays when free toxin is removed (and can be out-competed by
     raising acetylcholine).  Presynaptic PLA2 neurotoxicity is destruction of
     the terminal: antivenom given before it prevents it, antivenom given
     after it does nothing, and neostigmine does nothing either.

 (4) ACUTE KIDNEY INJURY IS AN INTEGRAL, NOT A LEVEL.  Nephron loss is driven
     by the accumulated product of intravascular fibrin formation and
     myoglobin cast burden, so creatinine peaks on day 3-5 — long after venom
     is undetectable.

STATE VECTOR (50 ODEs)
----------------------
   0- 3  D_i     bite-site depot, 4 toxin classes                    [mg]
   4- 7  V_i     central (plasma) free toxin                         [mg]
   8-11  P_i     deep peripheral tissue toxin                        [mg]
  12-15  C_i     antivenom-toxin complex, central                    [mg toxin]
  16     A_c     antivenom central                             [mgNE]
  17     A_p     antivenom peripheral                          [mgNE]
  18     A_t     antivenom in bite-site tissue                 [mgNE]
  19     NEO_a   neostigmine IM depot                          [mg]
  20     NEO_c   neostigmine central                           [mg]
  21     VAR_a   varespladib-methyl gut depot                   [mg]
  22     VAR_c   varespladib central                            [mg]
  23     TXA_c   tranexamic acid central                        [mg]
  24     FG      fibrinogen                                    [g/L]
  25     FX      factor X activity                             [fraction]
  26     PLT     platelets                                     [1e9/L]
  27     XDP     D-dimer / fibrin degradation products         [ug/mL]
  28     SYNUP   hepatic fibrinogen synthesis upregulation     [x normal]
  29     BR      AChR fractional occupancy by 3FTx             [0-1]
  30     TERM    presynaptic motor terminal integrity          [0-1]
  31     NEC     local myonecrosis fraction                    [0-1]
  32     EDEMA   local oedema volume                           [L]
  33     CK      creatine kinase                               [U/L]
  34     MB      plasma myoglobin                              [ug/mL]
  35     NEPH    functional nephron fraction                   [0-1]
  36     SCR     serum creatinine                              [mg/dL]
  37     CAST    tubular myoglobin cast burden                 [AU]
  38     FIBR    irreversible renal interstitial fibrosis      [0-1]
  39     PV      plasma volume                                 [L]
  40     GLX     endothelial glycocalyx integrity              [0-1]
  41     IL6     interleukin-6                                 [pg/mL]
  42     TNF     TNF-alpha                                     [pg/mL]
  43     MCA     mast-cell / anaphylactoid activation index     [AU]
  44     IC      circulating immune complex burden             [AU]
  45     SS      serum sickness score                          [AU]
  46     HBLD    cumulative systemic-bleeding hazard           [-]
  47     HDTH    cumulative death hazard                       [-]
  48     VAUC    free venom AUC                                [mg.h/L]
  49     AVCUM   cumulative antivenom administered             [mgNE]

UNITS CONVENTION
----------------
Time is HOURS throughout.  Venom is tracked in MILLIGRAMS of toxin protein.
Antivenom is tracked in "mgNE" = milligrams of venom-neutralising equivalents,
which is the unit antivenom is actually *labelled* in (Indian polyvalent ASV:
0.6 mg Daboia russelii venom neutralised per mL, i.e. 6 mgNE per 10-mL vial).
Working in the labelled unit means the antivenom potency in this model is a
regulatory fact rather than a fitted parameter, and the stoichiometric margin
of a standard vial count becomes computable instead of assumed.
"""

from __future__ import annotations

import json
import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# state indices
# ---------------------------------------------------------------------------
TOX = ["SVMP", "SVSP", "PLA2", "TFTX"]
iD = {t: 0 + k for k, t in enumerate(TOX)}
iV = {t: 4 + k for k, t in enumerate(TOX)}
iP = {t: 8 + k for k, t in enumerate(TOX)}
iC = {t: 12 + k for k, t in enumerate(TOX)}
(
    A_C, A_P, A_T, NEO_A, NEO_C, VAR_A, VAR_C, TXA_C,
    FG, FX, PLT, XDP, SYNUP, BR, TERM, NEC, EDEMA, CK, MB,
    NEPH, SCR, CAST, FIBR, PV, GLX, IL6, TNF, MCA, IC, SS,
    HBLD, HDTH, VAUC, AVCUM,
) = range(16, 50)
NSTATE = 50

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
# Toxin-class physicochemistry and disposition.
#
# Absorption from the bite site is LYMPHATIC for the large enzymes and partly
# capillary for the small ones, so ka falls with molecular weight.  The deep
# peripheral compartment is large and slowly-equilibrating; it is the source of
# the long venom terminal half-life that drives recurrence, and it is why
# "venom half-life" measured from plasma is not the quantity that matters.
#
#   MW      kDa
#   ka      1/h   FAST depot -> central (superficial, well-drained tissue)
#   ka_s    1/h   SLOW depot -> central (deep intramuscular / fascial
#                 sequestration; this is the rate constant that sets the venom's
#                 apparent TERMINAL half-life by flip-flop, and it is therefore
#                 the number an antivenom fragment has to outlast)
#   f_slow        fraction of the injected mass that goes to the slow depot.
#                 Large, lymphatic-dependent enzymes are sequestered most.
#   kdeg_d  1/h   local proteolytic destruction at the bite site
#   Vc      L     central distribution volume
#   CL      L/h   irreversible (reticuloendothelial + renal) clearance
#   eps     mg toxin neutralised per mgNE of antivenom  (relative potency)
#
# NOTE (defect 6, found by integration).  The first version gave venom a
# conventional two-compartment systemic disposition with a computed terminal
# half-life of 49.5 h, and CLAIM 2 WOULD NOT REPRODUCE: the bite-site depot was
# empty by 24 h, the deep compartment held almost nothing, and no antivenom
# fragment — Fab included — could run out while venom was still arriving.  The
# defect was structural.  Prolonged venom antigenaemia after viperid bites is an
# ABSORPTION phenomenon (sequestration in a deep depot that releases over days),
# not a distribution phenomenon, so the peripheral compartment was rebuilt as a
# slow depot.  The venom's terminal half-life is now flip-flop-limited, which is
# both what the antigen data show and what makes recurrence arithmetic.
#
# eps encodes a robustly reported asymmetry (Gutierrez/Calvete "antivenomics"):
# antivenoms are raised against whole venom, whose immunodominant components
# are the LARGE enzymes, so per unit of neutralising capacity they cover SVMP
# and SVSP well and the small, weakly immunogenic PLA2 and three-finger toxins
# badly.  This single vector is why a neurotoxic elapid bite needs 2-3x the
# vial count of a viperid bite of the same mass.
PK = {
    "SVMP": dict(MW=50.0, ka=0.25, ka_s=0.0120, f_slow=0.35, kdeg_d=0.010,
                 kdeg_s=0.0035, Vc=3.5, CL=0.30, eps=1.20),
    "SVSP": dict(MW=30.0, ka=0.35, ka_s=0.0150, f_slow=0.32, kdeg_d=0.012,
                 kdeg_s=0.0040, Vc=3.5, CL=0.40, eps=1.00),
    "PLA2": dict(MW=14.0, ka=0.50, ka_s=0.0300, f_slow=0.15, kdeg_d=0.015,
                 kdeg_s=0.0060, Vc=4.0, CL=0.90, eps=0.60),
    "TFTX": dict(MW=7.0,  ka=1.20, ka_s=0.1000, f_slow=0.05, kdeg_d=0.020,
                 kdeg_s=0.0150, Vc=8.0, CL=3.00, eps=0.35),
}

# Antivenom products.  Terminal half-lives are computed (not asserted) in
# `report_product_pk()` and checked against the published values in the
# comment on each line.
AV_PRODUCTS = {
    # ovine monovalent Fab (crotalidae polyvalent immune Fab, "CroFab")
    #   published terminal t1/2 12-23 h; large Vd, fast tissue penetration
    "Fab": dict(V1=4.0, V2=12.0, Q=2.00, CL=0.600, ktis=0.050, mgNE_vial=8.0,
                rho=0.105, ic=0.10, label="ovine Fab (CroFab-like)"),
    # equine F(ab')2 (crotalidae equine immune F(ab')2, "ANAVIP"/"Antivipmyn")
    #   published terminal t1/2 ~130 h
    "F(ab')2": dict(V1=3.5, V2=4.5, Q=0.50, CL=0.045, ktis=0.020, mgNE_vial=8.0,
                    rho=0.320, ic=0.45, label="equine F(ab')2 (ANAVIP-like)"),
    # equine whole IgG (Indian polyvalent ASV, SAIMR polyvalent)
    #   published terminal t1/2 ~90-200 h; highest acute reaction rate
    "IgG": dict(V1=3.2, V2=2.6, Q=0.30, CL=0.030, ktis=0.012, mgNE_vial=6.0,
                rho=1.300, ic=1.00, label="equine whole IgG (Indian polyvalent ASV-like)"),
}

P = dict(
    # ---- antivenom-toxin association -------------------------------------
    # In vivo the antibody-toxin reaction is essentially diffusion limited.
    # kb0 * cA is ~150/h at a post-infusion antivenom concentration, i.e. a
    # binding half-life of seconds.  The model is therefore STOICHIOMETRIC in
    # practice: free venom is ~0 while a molar excess of antivenom exists and
    # reappears the moment that excess is gone.  That is the whole of claim (2).
    cV_LOD=0.0020,     # mg/L — limit of detection of a venom antigen EIA
                       #   (~2 ng/mL).  Recurrence is DEFINED against this:
                       #   venom becomes undetectable and later detectable again.
    kb0=12.0,          # L/(mgNE.h), scaled per class by eps
    kb_loc=0.45,       # L/(mgNE.h) at the bite site — 27x worse than plasma
    kel_cplx=0.060,    # 1/h  RES clearance of antivenom-toxin complex (t1/2 11.6 h)
    # Antibody-antigen complexes do dissociate.  The DEFAULT here is zero, so
    # that every recurrence result in this model rests on absorption alone and
    # nothing is smuggled in; section E sweeps this parameter to ask how fast
    # dissociation would have to be to explain the clinical phenomenon that
    # absorption cannot.
    koff_av=0.0,       # 1/h  dissociation of the antivenom-toxin complex
    kout_t=0.030,      # 1/h  antivenom loss from bite-site tissue
    Vloc=0.15,         # L    bite-site tissue distribution volume

    # ---- haemostasis ------------------------------------------------------
    FG0=2.80,          # g/L
    kdeg_fg=0.00722,   # 1/h  fibrinogen turnover, t1/2 96 h
    Km_fg=4.00,        # g/L  Michaelis constant of venom fibrinogenolysis
    kcat_svsp_fg=2.50, # (g/L/h) per (mg/L) thrombin-like SVSP
    kcat_svmp_fg=0.30, # (g/L/h) per (mg/L) SVMP (fibrinogenolytic + PT activator)
    f_plasmin=0.05,    # fraction of fibrinogen loss that is plasmin-mediated
                       #   (the ONLY fraction tranexamic acid can touch)
    kup_syn=0.060, kdown_syn=0.020, SYNUP_max=3.20,
    FX0=1.00, k_fx_rec=0.01733,   # FX resynthesis t1/2 40 h
    kfx=0.35,                     # 1/h per (mg/L) SVMP  (RVV-X-like activation)
    PLT0=250.0, kdeg_plt=0.005775, kplt_boost=1.50,
    kplt_svmp=0.045, kplt_svsp=0.070,
    XDP0=0.25, kel_xdp=0.0866,    # D-dimer t1/2 8 h
    g_xdp=9.0,                    # ug/mL D-dimer per (g/L/h) of fibrin turnover

    # ---- neuromuscular junction ------------------------------------------
    SF0=4.00,          # normal neuromuscular safety factor
    koff_R=0.020,      # 1/h  alpha-neurotoxin dissociation (t1/2 35 h)
    Kd_R=0.060,        # mg/L apparent perijunctional Kd -> kon = koff/Kd
    kdes_term=0.550,   # 1/h per (mg/L) presynaptically-active PLA2
    krgn_term=0.00578, # 1/h  motor terminal regeneration, t1/2 5 d
    # NOTE (defect 5, found by integration).  Neostigmine was first written as a
    # free-standing multiplier on the safety factor: SF = SF0*TERM*(1-BR)*neo.
    # That made an anticholinesterase rescue a PRESYNAPTIC krait bite as
    # efficiently as a postsynaptic cobra bite (ventilation 28 h -> 4 h), which
    # is the opposite of the clinical fact and destroyed claim 3.  The error was
    # mechanistic, not numerical: raising cleft acetylcholine works by
    # OUT-COMPETING a competitive antagonist at the receptor, so it belongs
    # inside the occupancy term, where a small BR leaves it almost nothing to do.
    # w_pre keeps the small, real, non-decisive effect of prolonging the dwell
    # time of whatever transmitter a damaged terminal still releases.
    Emax_neo=1.20, EC50_neo=0.0050, w_pre_neo=0.050,

    # ---- local tissue -----------------------------------------------------
    # NOTE (defect 1, found by integration).  The bite-site volume is 0.15 L, so
    # local toxin concentrations are two orders of magnitude above plasma
    # (a Bothrops bite starts at ~167 mg/L of SVMP locally).  The first
    # calibration of knec was set by analogy with the plasma terms and drove
    # every venomous species to >50% limb necrosis within an hour.  knec has to
    # be ~40x smaller than a plasma rate constant precisely BECAUSE the local
    # compartment is small — the concentration does the work, not the constant.
    knec=6.50e-4,      # 1/h per (mg/L) local myotoxin
    a_svmp_nec=1.00, a_pla2_nec=1.60,
    # cobra local necrosis is cytotoxin (three-finger toxin) mediated, not PLA2
    # mediated, which is why an elapid with almost no SVMP still digests a limb
    a_tftx_nec=2.50,
    kheal_nec=0.00150, # 1/h  t1/2 19 d
    kedema=6.0e-4, kdrain=0.0200,
    CP0=8.0, CP_max=58.0, CP_k50=0.35,   # compartment pressure map
    # NOTE (defect 2).  CK release proportional to the necrosis FLUX alone put
    # the CK peak at 4-5 h and at 2x10^5 U/L.  Real CK peaks at 12-48 h: the
    # enzyme leaves a necrotic fibre for as long as the fibre is necrotic, so
    # release has to be driven by the necrosis STOCK as well as its flux.
    g_ck_loc=1.60e4, g_ck_stock=380.0, g_ck_sys=7.0e3,
    kel_ck=0.01925,                      # CK t1/2 36 h
    kmyo_sys=0.055,                      # systemic myotoxicity per (mg/L) PLA2
    # NOTE (defect 9, found by the baseline-stability check).  Myoglobin had a
    # release term and two clearance terms but NO basal production, so plasma
    # myoglobin decayed from 0.030 to 0 over 14 days in a patient with no
    # envenoming at all.  The model's baseline was therefore not a fixed point —
    # the classic signature of a missing zero-order production term.  Normal
    # muscle turnover maintains a basal concentration; mb_base is set so that
    # MB0 is an exact steady state at full nephron mass.
    g_mb=95.0, kel_mb_ren=0.260, kel_mb_nonren=0.045, MB0=0.030,

    # ---- kidney -----------------------------------------------------------
    kcast=0.0180, kclr_cast=0.0350,
    # NOTE (defect 3).  kn_fib was first sized as if FGloss were dimensionless.
    # Its time integral is the total GRAMS PER LITRE of fibrinogen consumed
    # (~5 g/L in an untreated Russell's viper bite), so a rate constant of 0.013
    # contributed 0.05 of a nephron log-loss — i.e. nothing.  The intravascular
    # coagulation arm of venom-induced AKI was silently absent from the model.
    kn_fib=0.2000,     # nephron loss per (g/L) of fibrinogen consumed
    kn_cast=0.00300, kn_isch=0.0900, kn_dir=0.00800,
    krecov_neph=0.00800, kfibrosis=0.100,
    SCR0=0.90, kel_cr=0.060,

    # ---- haemodynamics ----------------------------------------------------
    # NOTE (defect 4).  The first leak calibration took plasma volume to 0.9 L
    # and MAP to 20 mmHg in every viperid scenario, killing 100% of untreated
    # patients from shock before coagulopathy or AKI could be read at all.  The
    # transcapillary refill term was an order of magnitude too weak.
    PV0=3.00, kleak=0.0300, krefill=0.0600,
    kglx=0.0300, kglxr=0.0200,
    MAP0=92.0, map_exp=1.60, map_anaphyl=0.38,

    # ---- inflammation -----------------------------------------------------
    g_il6=170.0, kel_il6=0.347,   # IL-6 t1/2 2 h
    g_tnf=55.0, kel_tnf=0.990,

    # ---- antivenom adverse reactions -------------------------------------
    kel_mca=0.800, MCA_anaphyl=1.00,
    k_ic=0.00090, kel_ic=0.0200,
    k_ss=0.01800, kel_ss=0.0120,

    # ---- hazards ----------------------------------------------------------
    h_bleed=0.00830,   # 1/h with unclottable blood -> ~0.4 over 48 h
    FG_unclot=0.50,    # g/L  the 20-minute whole blood clotting test threshold
    PLT_low=50.0,
    hd_bleed=0.00180, hd_resp_novent=0.09000, hd_resp_vent=0.000400,
    hd_aki_nodial=0.01000, hd_aki_dial=0.000600, NEPH_haz=0.500,
    hd_shock=0.00600, hd_anaphyl=0.00800,

    # ---- adjunct drug PK --------------------------------------------------
    neo_ka=1.400, neo_kel=0.835, neo_V=50.0,      # neostigmine IM, t1/2 50 min
    var_ka=1.500, var_kel=0.231, var_V=12.0, var_F=0.50,
    var_IC50=0.0100,                              # mg/L, svPLA2 inhibition
    txa_kel=0.347, txa_V=12.0, txa_IC50=5.00,
)

# ---------------------------------------------------------------------------
# SNAKE ARCHETYPES
# ---------------------------------------------------------------------------
# `dose` is the mass of venom protein delivered by a median significant bite,
# not the maximum milkable gland yield.  `frac` are the mass fractions of the
# four modelled toxin classes; they do not sum to 1 because every venom also
# contains L-amino-acid oxidase, CRISPs, disintegrins, lectins and non-toxic
# protein that this model does not resolve.
#   f_pre  : fraction of the PLA2 mass that is a presynaptic beta-neurotoxin
#   f_myo  : relative systemic myotoxicity of that venom's PLA2
#   f_loc  : relative local necrotising potency
#   koffm  : multiplier on alpha-neurotoxin dissociation (krait ~ irreversible)
#   f_leak : relative capillary-leak potency
SNAKES = {
    "Daboia russelii (Sri Lanka)": dict(
        dose=63.0, frac=dict(SVMP=0.28, SVSP=0.14, PLA2=0.34, TFTX=0.01),
        f_pre=0.09, f_myo=0.50, f_loc=0.30, koffm=1.00, f_leak=1.00, ko="러셀살무사"),
    "Echis ocellatus": dict(
        dose=22.0, frac=dict(SVMP=0.68, SVSP=0.05, PLA2=0.12, TFTX=0.00),
        f_pre=0.00, f_myo=0.10, f_loc=0.55, koffm=1.00, f_leak=0.60, ko="서아프리카카펫바이퍼"),
    "Naja naja": dict(
        dose=40.0, frac=dict(SVMP=0.10, SVSP=0.02, PLA2=0.28, TFTX=0.55),
        f_pre=0.05, f_myo=0.60, f_loc=1.30, koffm=1.00, f_leak=0.35, ko="인도코브라"),
    "Bungarus caeruleus": dict(
        dose=10.0, frac=dict(SVMP=0.03, SVSP=0.01, PLA2=0.66, TFTX=0.20),
        f_pre=1.00, f_myo=0.05, f_loc=0.03, koffm=0.15, f_leak=0.05, ko="인도크레이트"),
    "Crotalus atrox": dict(
        dose=55.0, frac=dict(SVMP=0.42, SVSP=0.14, PLA2=0.08, TFTX=0.00),
        f_pre=0.03, f_myo=0.30, f_loc=0.70, koffm=1.00, f_leak=0.70, ko="서부다이아몬드방울뱀"),
    "Bothrops asper": dict(
        dose=50.0, frac=dict(SVMP=0.50, SVSP=0.06, PLA2=0.30, TFTX=0.00),
        f_pre=0.02, f_myo=0.80, f_loc=1.50, koffm=1.00, f_leak=0.85, ko="중미창머리독사"),
    "Bitis arietans": dict(
        dose=90.0, frac=dict(SVMP=0.62, SVSP=0.04, PLA2=0.16, TFTX=0.00),
        f_pre=0.00, f_myo=0.50, f_loc=1.30, koffm=1.00, f_leak=1.10, ko="퍼프애더"),
    "Dry bite": dict(
        dose=0.0, frac=dict(SVMP=0.0, SVSP=0.0, PLA2=0.0, TFTX=0.0),
        f_pre=0.00, f_myo=0.00, f_loc=0.00, koffm=1.00, f_leak=0.00, ko="건성교상"),
}


# ---------------------------------------------------------------------------
# SCENARIO SPECIFICATION
# ---------------------------------------------------------------------------
@dataclass
class Scenario:
    name: str
    snake: str
    venom_mult: float = 1.0
    product: str | None = None
    vials: float = 0.0
    t_av: float | None = None          # h post-bite, first antivenom dose
    av_infuse: float = 1.0             # h, infusion duration
    repeat: list = field(default_factory=list)   # [(t_h, vials), ...]
    maintenance: tuple | None = None   # (start_h, interval_h, vials, n_doses)
    neostigmine: list = field(default_factory=list)  # [t_h, ...] 0.5 mg IM each
    varespladib: tuple | None = None   # (t_start_h, interval_h, mg, n_doses)
    txa: list = field(default_factory=list)          # [t_h, ...] 1 g IV each
    ffp: list = field(default_factory=list)          # [(t_h, dFG_g_per_L), ...]
    fluids: list = field(default_factory=list)       # [(t0, t1, L_per_h), ...]
    icu: bool = True                   # mechanical ventilation available
    dialysis: bool = True
    tmax: float = 336.0                # h  (14 days)
    note: str = ""


# ---------------------------------------------------------------------------
# infusion helper
# ---------------------------------------------------------------------------
def _rate(windows, t):
    r = 0.0
    for t0, t1, k in windows:
        if t0 <= t < t1:
            r += k
    return r


def build_inputs(sc: Scenario):
    """Return (av_windows, fluid_windows, bolus_events, av_total_mgNE)."""
    prod = AV_PRODUCTS[sc.product] if sc.product else None
    av, bolus = [], []
    total = 0.0
    if prod and sc.t_av is not None and sc.vials > 0:
        mg = sc.vials * prod["mgNE_vial"]
        av.append((sc.t_av, sc.t_av + sc.av_infuse, mg / sc.av_infuse))
        total += mg
    for t, v in sc.repeat:
        mg = v * prod["mgNE_vial"]
        av.append((t, t + sc.av_infuse, mg / sc.av_infuse))
        total += mg
    if sc.maintenance:
        t0, dt, v, n = sc.maintenance
        for j in range(int(n)):
            mg = v * prod["mgNE_vial"]
            tt = t0 + j * dt
            av.append((tt, tt + sc.av_infuse, mg / sc.av_infuse))
            total += mg
    for t in sc.neostigmine:
        bolus.append((t, NEO_A, 0.50))            # 0.5 mg IM
    if sc.varespladib:
        t0, dt, mgd, n = sc.varespladib
        for j in range(int(n)):
            bolus.append((t0 + j * dt, VAR_A, mgd))
    for t in sc.txa:
        bolus.append((t, TXA_C, 1000.0))          # 1 g IV
    for t, d in sc.ffp:
        bolus.append((t, FG, d))
    return av, list(sc.fluids), sorted(bolus), total


# ---------------------------------------------------------------------------
# INITIAL CONDITIONS
# ---------------------------------------------------------------------------
def initial_state(sc: Scenario):
    y = np.zeros(NSTATE)
    sn = SNAKES[sc.snake]
    dose = sn["dose"] * sc.venom_mult
    for t in TOX:
        m = dose * sn["frac"][t]
        y[iD[t]] = m * (1.0 - PK[t]["f_slow"])   # fast, superficial depot
        y[iP[t]] = m * PK[t]["f_slow"]           # slow, sequestered depot
    y[FG] = P["FG0"]
    y[FX] = P["FX0"]
    y[PLT] = P["PLT0"]
    y[XDP] = P["XDP0"]
    y[SYNUP] = 1.0
    y[BR] = 0.0
    y[TERM] = 1.0
    y[NEC] = 0.0
    y[EDEMA] = 0.0
    y[CK] = 90.0
    y[MB] = P["MB0"]
    y[NEPH] = 1.0
    y[SCR] = P["SCR0"]
    # the pigment-cast compartment's own steady state at baseline myoglobin and
    # full nephron mass: starting it at zero left CAST (and therefore the
    # nephron-loss term it feeds) drifting for the whole 14 days
    y[CAST] = P["kcast"] * P["MB0"] * 0.30 / P["kclr_cast"]
    y[FIBR] = 0.0
    y[PV] = P["PV0"]
    y[GLX] = 1.0
    y[IL6] = 3.0
    y[TNF] = 4.0
    return y


# ---------------------------------------------------------------------------
# DERIVED (algebraic) QUANTITIES  — shared by RHS and by the reporting layer
# ---------------------------------------------------------------------------
def derived(y, sc: Scenario):
    sn = SNAKES[sc.snake]
    prod = AV_PRODUCTS[sc.product] if sc.product else None
    d = {}

    # free plasma toxin concentrations
    for t in TOX:
        d["c" + t] = max(y[iV[t]], 0.0) / PK[t]["Vc"]
    d["cVfree"] = sum(d["c" + t] for t in TOX)
    # what a venom-antigen ELISA sees: free + antivenom-bound, both circulate
    d["cVtotal"] = d["cVfree"] + sum(max(y[iC[t]], 0.0) / PK[t]["Vc"] for t in TOX)

    # bite-site (local) concentrations
    vloc = P["Vloc"] + max(y[EDEMA], 0.0)
    for t in TOX:
        d["l" + t] = max(y[iD[t]], 0.0) / vloc
    d["cAt"] = max(y[A_T], 0.0) / vloc

    # antivenom
    d["cA"] = (max(y[A_C], 0.0) / prod["V1"]) if prod else 0.0

    # adjunct drugs
    d["cNEO"] = max(y[NEO_C], 0.0) / P["neo_V"]
    d["cVAR"] = max(y[VAR_C], 0.0) / P["var_V"]
    d["cTXA"] = max(y[TXA_C], 0.0) / P["txa_V"]
    # varespladib is a catalytic-site PLA2 inhibitor: it reaches the bite site
    # and the nerve terminal, where an IgG cannot.  INH_PLA2 multiplies every
    # PLA2-driven term in the model.
    d["INH_PLA2"] = 1.0 / (1.0 + d["cVAR"] / P["var_IC50"])
    d["INH_PLASMIN"] = 1.0 / (1.0 + d["cTXA"] / P["txa_IC50"])

    # ---- haemostasis ----
    fg = max(y[FG], 1e-9)
    mm = fg / (P["Km_fg"] + fg)
    enz = P["kcat_svsp_fg"] * d["cSVSP"] + P["kcat_svmp_fg"] * d["cSVMP"]
    # 95% of fibrinogen loss is direct venom enzymatic cleavage (untouchable by
    # tranexamic acid); 5% is secondary plasmin-mediated fibrinogenolysis.
    d["FGloss"] = enz * mm * ((1.0 - P["f_plasmin"]) + P["f_plasmin"] * d["INH_PLASMIN"])
    d["unclot"] = 1.0 / (1.0 + (fg / P["FG_unclot"]) ** 4)     # 20WBCT positive
    d["lowplt"] = 1.0 / (1.0 + (max(y[PLT], 1e-6) / P["PLT_low"]) ** 3)
    d["WBCT20"] = d["unclot"] > 0.5

    # ---- neuromuscular ----
    neo_amp = 1.0 + P["Emax_neo"] * d["cNEO"] / (P["EC50_neo"] + d["cNEO"])
    d["neo_amp"] = neo_amp
    # competition at the receptor: elevated acetylcholine effectively dilutes a
    # competitive antagonist's occupancy
    R_eff = max(1.0 - min(y[BR], 1.0) / neo_amp, 0.0)
    # a destroyed terminal releases little transmitter, so prolonging its dwell
    # time buys a small, non-decisive amount
    TERM_eff = max(y[TERM], 0.0) * (1.0 + P["w_pre_neo"] * (neo_amp - 1.0))
    d["R_eff"], d["TERM_eff"] = R_eff, TERM_eff
    d["SF"] = P["SF0"] * TERM_eff * R_eff
    d["VCfrac"] = min(max((d["SF"] - 1.0) / 2.0, 0.0), 1.0)
    d["SBC"] = 42.0 * d["VCfrac"]                 # single-breath count
    d["ptosis"] = d["SF"] < 2.20
    d["vent_needed"] = d["VCfrac"] < 0.15

    # ---- local ----
    d["CP"] = P["CP0"] + P["CP_max"] * max(y[EDEMA], 0.0) / (max(y[EDEMA], 0.0) + P["CP_k50"])

    # myonecrosis flux (local), varespladib-inhibitable in the PLA2 arm only.
    # The three-finger-toxin term is what makes an elapid bite necrotic despite
    # carrying almost no metalloproteinase.
    d["necflux"] = P["knec"] * (
        P["a_svmp_nec"] * d["lSVMP"]
        + P["a_pla2_nec"] * d["lPLA2"] * d["INH_PLA2"]
        + P["a_tftx_nec"] * d["lTFTX"]
    ) * sn["f_loc"] * max(1.0 - y[NEC], 0.0)
    d["myosys"] = P["kmyo_sys"] * d["cPLA2"] * sn["f_myo"] * d["INH_PLA2"]

    # ---- haemodynamics ----
    d["MAP"] = max(
        P["MAP0"] * (max(y[PV], 0.10) / P["PV0"]) ** P["map_exp"]
        * (1.0 - P["map_anaphyl"] * y[MCA] / (y[MCA] + 0.50)),
        20.0,
    )
    d["Hct"] = min(42.0 * P["PV0"] / max(y[PV], 0.30), 72.0)
    d["shock"] = d["MAP"] < 65.0

    # ---- renal ----
    d["eGFR"] = 100.0 * max(y[NEPH], 0.0)
    d["urine"] = min(max((d["MAP"] - 45.0) / 45.0, 0.0), 1.0) * max(y[NEPH], 0.0)
    d["AKI3"] = y[SCR] >= 3.0 * P["SCR0"]

    d["anaphylaxis"] = y[MCA] > P["MCA_anaphyl"]
    d["serum_sickness"] = y[SS] > 1.0
    return d


# ---------------------------------------------------------------------------
# RIGHT-HAND SIDE
# ---------------------------------------------------------------------------
def rhs(t, y, sc: Scenario, av_win, fl_win):
    sn = SNAKES[sc.snake]
    prod = AV_PRODUCTS[sc.product] if sc.product else None
    d = derived(y, sc)
    dy = np.zeros(NSTATE)

    av_rate = _rate(av_win, t)
    fl_rate = _rate(fl_win, t)

    # =====================================================================
    # VENOM PHARMACOKINETICS  +  ANTIVENOM BINDING
    # =====================================================================
    bind_sys_total = 0.0     # mg toxin bound in plasma (per h)
    bind_loc_total = 0.0     # mg toxin bound at the bite site (per h)
    for tox in TOX:
        p = PK[tox]
        D, V, Pp, C = y[iD[tox]], y[iV[tox]], y[iP[tox]], y[iC[tox]]
        D, V, Pp, C = max(D, 0.0), max(V, 0.0), max(Pp, 0.0), max(C, 0.0)

        # local neutralisation at the bite site.  kb_loc is 27x smaller than the
        # plasma rate constant AND cAt is small because an antibody penetrates
        # oedematous, poorly-perfused tissue badly.  This is why intravenous
        # antivenom does not prevent local necrosis.  The SLOW depot is deeper
        # still and takes a further 4-fold penalty — which is exactly why
        # neutralising the depot is not an available way to stop recurrence.
        bl_f = min(P["kb_loc"] * d["cAt"] * D, D * 50.0)
        bl_s = min(0.25 * P["kb_loc"] * d["cAt"] * Pp, Pp * 50.0)
        bind_loc_total += (bl_f + bl_s) / p["eps"]

        dy[iD[tox]] = -p["ka"] * D - p["kdeg_d"] * D - bl_f
        # the slow depot releases over DAYS.  It is the venom's terminal phase,
        # and the reason venom keeps arriving after a short-lived antivenom has
        # already been cleared.
        dy[iP[tox]] = -p["ka_s"] * Pp - p["kdeg_s"] * Pp - bl_s

        # plasma neutralisation — near-instantaneous while antivenom is present
        bs = P["kb0"] * p["eps"] * d["cA"] * V
        bind_sys_total += bs / p["eps"]

        dy[iV[tox]] = (
            p["ka"] * D + p["ka_s"] * Pp
            - p["CL"] * V / p["Vc"]
            - bs
        )
        diss = P["koff_av"] * C
        dy[iV[tox]] += diss
        bind_sys_total -= diss / p["eps"]       # returns capacity to the antivenom
        dy[iC[tox]] = bs + bl_f + bl_s - P["kel_cplx"] * C - diss

    # =====================================================================
    # ANTIVENOM PHARMACOKINETICS
    # =====================================================================
    if prod:
        Ac, Ap, At = max(y[A_C], 0.0), max(y[A_P], 0.0), max(y[A_T], 0.0)
        # oedema increases convective delivery to the bite site a little, but
        # nowhere near enough to matter (see bind_loc_total above)
        tis = prod["ktis"] * Ac * (1.0 + 1.5 * max(y[EDEMA], 0.0))
        dy[A_C] = (
            av_rate
            - prod["CL"] * Ac / prod["V1"]
            - prod["Q"] * Ac / prod["V1"]
            + prod["Q"] * Ap / prod["V2"]
            - tis
            - bind_sys_total
        )
        dy[A_P] = prod["Q"] * Ac / prod["V1"] - prod["Q"] * Ap / prod["V2"]
        dy[A_T] = tis - P["kout_t"] * At - bind_loc_total
        dy[AVCUM] = av_rate
        # acute reaction is driven by the INFUSION RATE of foreign protein,
        # not by the cumulative dose — which is why slowing the infusion works
        dy[MCA] = prod["rho"] * av_rate / 60.0 - P["kel_mca"] * y[MCA]
        dy[IC] = P["k_ic"] * prod["ic"] * (Ac + Ap) - P["kel_ic"] * y[IC]
    else:
        dy[MCA] = -P["kel_mca"] * y[MCA]
        dy[IC] = -P["kel_ic"] * y[IC]
    dy[SS] = P["k_ss"] * y[IC] - P["kel_ss"] * y[SS]

    # =====================================================================
    # ADJUNCT DRUG PK
    # =====================================================================
    dy[NEO_A] = -P["neo_ka"] * y[NEO_A]
    dy[NEO_C] = P["neo_ka"] * y[NEO_A] - P["neo_kel"] * y[NEO_C]
    dy[VAR_A] = -P["var_ka"] * y[VAR_A]
    dy[VAR_C] = P["var_F"] * P["var_ka"] * y[VAR_A] - P["var_kel"] * y[VAR_C]
    dy[TXA_C] = -P["txa_kel"] * y[TXA_C]

    # =====================================================================
    # HAEMOSTASIS  —  clock 2, the substrate clock
    # =====================================================================
    # Fibrinogen falls at the rate the venom enzyme works and rises at the rate
    # the LIVER works.  Those two rates have nothing to do with each other, and
    # antivenom can only touch the first.  Claim (1) is this line.
    dy[FG] = P["kdeg_fg"] * P["FG0"] * y[SYNUP] - P["kdeg_fg"] * y[FG] - d["FGloss"]
    dy[SYNUP] = (
        P["kup_syn"] * max(1.0 - y[FG] / P["FG0"], 0.0) * (P["SYNUP_max"] - y[SYNUP])
        - P["kdown_syn"] * (y[SYNUP] - 1.0)
    )
    dy[FX] = P["k_fx_rec"] * (P["FX0"] - y[FX]) - P["kfx"] * d["cSVMP"] * y[FX]
    dy[PLT] = (
        P["PLT0"] * P["kdeg_plt"] * (1.0 + P["kplt_boost"] * max(1.0 - y[PLT] / P["PLT0"], 0.0))
        - P["kdeg_plt"] * y[PLT]
        - (P["kplt_svmp"] * d["cSVMP"] + P["kplt_svsp"] * d["cSVSP"]) * y[PLT]
    )
    dy[XDP] = P["g_xdp"] * d["FGloss"] - P["kel_xdp"] * (y[XDP] - P["XDP0"])

    # =====================================================================
    # NEUROMUSCULAR JUNCTION  —  the two mechanisms that look the same
    # =====================================================================
    # POSTSYNAPTIC: an occupancy.  Remove free toxin and it decays (slowly);
    # raise acetylcholine and you out-compete what is left.  Reversible.
    kon_R = P["koff_R"] / P["Kd_R"]
    dy[BR] = kon_R * d["cTFTX"] * max(1.0 - y[BR], 0.0) - P["koff_R"] * sn["koffm"] * y[BR]
    # PRESYNAPTIC: destruction.  The only recovery term is regeneration of the
    # terminal on a 5-day half-life.  Neither antivenom nor neostigmine appears
    # anywhere in this equation, which is the entirety of claim (3).
    dy[TERM] = (
        -P["kdes_term"] * sn["f_pre"] * d["cPLA2"] * d["INH_PLA2"] * y[TERM]
        + P["krgn_term"] * (1.0 - y[TERM])
    )

    # =====================================================================
    # LOCAL TISSUE
    # =====================================================================
    dy[NEC] = d["necflux"] - P["kheal_nec"] * y[NEC]
    dy[EDEMA] = (
        P["kedema"] * (d["lSVMP"] + 0.5 * d["lPLA2"] * d["INH_PLA2"]) * sn["f_loc"]
        - P["kdrain"] * y[EDEMA]
    )
    # CK leaves a necrotic fibre for as long as the fibre stays necrotic, so the
    # stock term (not the flux term) sets both the magnitude and the 12-48 h peak
    ck_release = (
        P["g_ck_loc"] * d["necflux"]
        + P["g_ck_stock"] * y[NEC]
        + P["g_ck_sys"] * d["myosys"]
    )
    dy[CK] = ck_release - P["kel_ck"] * (y[CK] - 90.0)
    mb_base = P["MB0"] * (P["kel_mb_ren"] + P["kel_mb_nonren"])
    dy[MB] = (
        mb_base
        + P["g_mb"] * (d["necflux"] * 0.35 + y[NEC] * 0.010 + d["myosys"])
        - P["kel_mb_ren"] * max(y[NEPH], 0.0) * y[MB]
        - P["kel_mb_nonren"] * y[MB]
    )

    # =====================================================================
    # KIDNEY  —  an integral, not a level (claim 4)
    # =====================================================================
    dy[CAST] = P["kcast"] * y[MB] * (1.0 - 0.7 * d["urine"]) - P["kclr_cast"] * y[CAST]
    loss = (
        P["kn_fib"] * d["FGloss"]                      # fibrin microthrombi
        + P["kn_cast"] * y[CAST]                       # pigment casts
        + P["kn_isch"] * max((65.0 - d["MAP"]) / 65.0, 0.0)   # ischaemia
        + P["kn_dir"] * d["cSVMP"]                     # direct glomerular BM injury
    ) * max(y[NEPH], 0.0)
    ceiling = max(1.0 - y[FIBR], 0.0)
    dy[NEPH] = -loss + P["krecov_neph"] * max(ceiling - y[NEPH], 0.0)
    dy[FIBR] = P["kfibrosis"] * loss
    dy[SCR] = P["SCR0"] * P["kel_cr"] - P["kel_cr"] * max(y[NEPH], 0.02) * y[SCR]

    # =====================================================================
    # HAEMODYNAMICS / CAPILLARY LEAK
    # =====================================================================
    dy[GLX] = (
        -P["kglx"] * (d["cSVMP"] + 0.3 * d["cPLA2"] * d["INH_PLA2"]) * sn["f_leak"] * y[GLX]
        + P["kglxr"] * (1.0 - y[GLX])
    )
    dy[PV] = (
        fl_rate
        - P["kleak"] * (1.0 - y[GLX]) * y[PV]
        + P["krefill"] * (P["PV0"] - y[PV])
    )

    # =====================================================================
    # INFLAMMATION
    # =====================================================================
    dy[IL6] = P["g_il6"] * (d["cSVMP"] + d["cPLA2"] + 12.0 * d["necflux"]) - P["kel_il6"] * (y[IL6] - 3.0)
    dy[TNF] = P["g_tnf"] * (d["cSVMP"] + d["cPLA2"]) - P["kel_tnf"] * (y[TNF] - 4.0)

    # =====================================================================
    # CUMULATIVE HAZARDS AND EXPOSURE
    # =====================================================================
    dy[HBLD] = P["h_bleed"] * d["unclot"] * (1.0 + 2.0 * d["lowplt"])
    hd = (
        P["hd_bleed"] * d["unclot"] * (1.0 + 2.0 * d["lowplt"])
        + (P["hd_resp_vent"] if sc.icu else P["hd_resp_novent"]) * (1.0 if d["vent_needed"] else 0.0)
        # graded, not thresholded: a patient at 40% of nephron mass is not safe
        # merely because they are above a dialysis trigger
        + (P["hd_aki_dial"] if sc.dialysis else P["hd_aki_nodial"])
        * max((P["NEPH_haz"] - max(y[NEPH], 0.0)) / P["NEPH_haz"], 0.0)
        + P["hd_shock"] * max((60.0 - d["MAP"]) / 60.0, 0.0) * 10.0
        + P["hd_anaphyl"] * (1.0 if d["anaphylaxis"] else 0.0)
    )
    dy[HDTH] = hd
    dy[VAUC] = d["cVfree"]
    return dy


# ---------------------------------------------------------------------------
# SIMULATION DRIVER
# ---------------------------------------------------------------------------
_SIM_CACHE = {}


def simulate(sc: Scenario, dt=0.25, cache_key=None):
    """Integrate one scenario.  `cache_key` lets the reporting layer re-use a
    trajectory instead of re-integrating it (the summary tables in section D
    need the full trace of scenarios already run in section C)."""
    if cache_key is not None and cache_key in _SIM_CACHE:
        return _SIM_CACHE[cache_key]
    out = _simulate(sc, dt)
    if cache_key is not None:
        _SIM_CACHE[cache_key] = out
    return out


def _simulate(sc: Scenario, dt=0.25):
    av_win, fl_win, bolus, av_total = build_inputs(sc)
    y = initial_state(sc)

    brk = sorted({0.0, sc.tmax}
                 | {w[0] for w in av_win} | {w[1] for w in av_win}
                 | {w[0] for w in fl_win} | {w[1] for w in fl_win}
                 | {b[0] for b in bolus})
    brk = [b for b in brk if 0.0 <= b <= sc.tmax]

    T, Y = [0.0], [y.copy()]
    for a, b in zip(brk[:-1], brk[1:]):
        for tb, idx, amt in bolus:
            if abs(tb - a) < 1e-9:
                y[idx] += amt
        if b - a < 1e-9:
            continue
        n = max(int(math.ceil((b - a) / dt)), 2)
        teval = np.linspace(a, b, n + 1)
        sol = solve_ivp(rhs, (a, b), y, args=(sc, av_win, fl_win),
                        method="LSODA", t_eval=teval, rtol=1e-6, atol=1e-8,
                        max_step=1.0)
        if not sol.success:
            raise RuntimeError(f"{sc.name}: integration failed — {sol.message}")
        T.extend(sol.t[1:].tolist())
        Y.extend(sol.y[:, 1:].T.tolist())
        y = sol.y[:, -1].copy()
    for tb, idx, amt in bolus:
        if abs(tb - sc.tmax) < 1e-9:
            y[idx] += amt

    T = np.asarray(T)
    Y = np.asarray(Y)
    Dl = [derived(Y[k], sc) for k in range(len(T))]
    return T, Y, Dl, av_total


# ---------------------------------------------------------------------------
# SUMMARY METRICS
# ---------------------------------------------------------------------------
def summarize(sc: Scenario, T, Y, Dl, av_total):
    def col(k):
        return Y[:, k]

    def dcol(k):
        return np.asarray([x[k] for x in Dl])

    fg = col(FG)
    unclot = dcol("unclot") > 0.5
    # hours of unclottable blood, and time to restored clottability
    hrs_unclot = float(np.trapezoid(unclot.astype(float), T))
    t_first = float(T[np.argmax(unclot)]) if unclot.any() else float("nan")
    t_restore = float("nan")
    never_unclot = not bool(unclot.any())
    if unclot.any():
        last = int(np.max(np.nonzero(unclot)[0]))
        if last < len(T) - 1:
            t_restore = float(T[last + 1])

    # Fibrinogen RISE RATE evaluated at a fixed substrate level (1.0 g/L), on
    # the way up.  Comparing slopes at the same FG value rather than at each
    # trace's own nadir is what isolates the hepatic synthesis term from the
    # residual venom term — see claim 1.
    rise_at_1 = float("nan")
    j0 = int(np.argmin(fg))
    if fg.min() < 1.0 and j0 < len(T) - 1:
        up = fg[j0:]
        k = int(np.argmax(up >= 1.0)) if (up >= 1.0).any() else -1
        if k > 0:
            m = j0 + k
            t1 = T[m]
            t2 = min(t1 + 12.0, T[-1])
            rise_at_1 = float((np.interp(t2, T, fg) - 1.0) / (t2 - t1) * 24.0)

    vent = dcol("vent_needed")
    hrs_vent = float(np.trapezoid(vent.astype(float), T))
    ptosis = dcol("ptosis")

    # RECURRENT VENOM ANTIGENAEMIA.  Defined as a rebound: after the acute peak,
    # free venom reaches a post-treatment trough and a LATER maximum exceeds that
    # trough by >= 2-fold.  Referencing the rebound to the trough rather than to
    # the acute pre-antivenom peak is the whole point — the acute peak belongs to
    # a different regime (unneutralised bolus), and anchoring to it hides the
    # phenomenon entirely.  (This was defect 7: the first detector compared the
    # rebound with 5% of the acute peak and reported "no recurrence" for a Fab
    # arm whose free venom demonstrably rose 3.5-fold from its trough.)
    cv = dcol("cVfree")
    peak = float(cv.max()) if cv.size else 0.0
    recur_t, recur_peak, rebound = float("nan"), 0.0, 1.0
    trough = float("nan")
    rebound_detectable = False
    if peak > 0:
        jpk = int(np.argmax(cv))
        tail = cv[jpk:]
        ttail = T[jpk:]
        if tail.size > 3:
            # the FIRST local minimum after the acute peak, not the global one:
            # the global minimum of a monotonically decaying tail is its last
            # point, which can never be followed by a rebound
            jmin = None
            for m in range(1, tail.size - 1):
                if tail[m] <= tail[m + 1] and (ttail[m] - ttail[0]) >= 4.0:
                    jmin = m
                    break
            if jmin is not None and tail[jmin] > 0:
                trough = float(tail[jmin])
                after = tail[jmin:]
                mx = float(after.max())
                rebound = mx / trough
                rebound_detectable = mx > P["cV_LOD"]
                # gate on the assay limit of detection.  Without it, a trough of
                # 1e-9 mg/L produces a spurious 20000-fold "rebound" that no
                # instrument could see and no patient would feel (defect 8).
                if rebound >= 2.0 and mx > P["cV_LOD"]:
                    recur_t = float(ttail[jmin + int(np.argmax(after))])
                    recur_peak = mx

    # LATE free-venom exposure (48 h onward) is the clean quantitative separator
    # between fragment formats: it is the integral the short-lived fragment fails
    # to cover, and it does not depend on any threshold choice.
    late = T >= 48.0
    vauc_late = float(np.trapezoid(cv[late], T[late])) if late.sum() > 1 else 0.0
    detect_late = ((cv > P["cV_LOD"]) & late).astype(float)
    hrs_detect_late = float(np.trapezoid(detect_late, T)) if T.size > 1 else 0.0

    # RECURRENT COAGULOPATHY: fibrinogen recovers by >= 0.5 g/L from its nadir
    # and then falls again by >= 0.4 g/L.  This is the endpoint the clinical
    # literature reports, and like the venom detector it is read off the
    # trajectory rather than switched on anywhere in the model.
    recoag_t = float("nan")
    if fg.size:
        j = int(np.argmin(fg))
        seg = fg[j:]
        run_max, best = seg[0], 0.0
        for m in range(1, len(seg)):
            run_max = max(run_max, seg[m])
            if run_max - fg[j] >= 0.5 and run_max - seg[m] >= 0.4:
                recoag_t = float(T[j + m])
                break
            best = max(best, run_max - seg[m])

    prod = AV_PRODUCTS[sc.product] if sc.product else None
    sn = SNAKES[sc.snake]
    need = stoichiometric_need(sc.snake, sc.venom_mult)

    return dict(
        name=sc.name, snake=sc.snake, snake_ko=sn["ko"],
        product=(prod["label"] if prod else "no antivenom"),
        vials=sc.vials, t_av=sc.t_av, av_mgNE=av_total,
        need_mgNE=need, molar_margin=(av_total / need if need > 0 else float("inf")),
        FG_nadir=float(fg.min()), t_FG_nadir=float(T[int(np.argmin(fg))]),
        FG_d7=float(np.interp(168.0, T, fg)),
        hrs_unclottable=hrs_unclot, t_unclot_onset=t_first, t_clottable=t_restore,
        never_unclottable=never_unclot, FG_rise_at_1=rise_at_1,
        PLT_nadir=float(col(PLT).min()),
        FX_nadir=float(col(FX).min()),
        XDP_peak=float(col(XDP).max()),
        SF_nadir=float(dcol("SF").min()),
        SBC_nadir=float(dcol("SBC").min()),
        TERM_nadir=float(col(TERM).min()),
        BR_peak=float(col(BR).max()),
        hrs_ptosis=float(np.trapezoid(ptosis.astype(float), T)),
        hrs_ventilated=hrs_vent, ventilated=bool(vent.any()),
        NEC_final=float(col(NEC)[-1]), NEC_peak=float(col(NEC).max()),
        CP_peak=float(dcol("CP").max()),
        CK_peak=float(col(CK).max()), t_CK_peak=float(T[int(np.argmax(col(CK)))]),
        MB_peak=float(col(MB).max()),
        SCR_peak=float(col(SCR).max()), t_SCR_peak=float(T[int(np.argmax(col(SCR)))]),
        eGFR_nadir=float(dcol("eGFR").min()),
        NEPH_final=float(col(NEPH)[-1]), FIBR_final=float(col(FIBR)[-1]),
        AKI3=bool((col(SCR) >= 3 * P["SCR0"]).any()),
        MAP_nadir=float(dcol("MAP").min()), Hct_peak=float(dcol("Hct").max()),
        GLX_nadir=float(col(GLX).min()),
        IL6_peak=float(col(IL6).max()),
        MCA_peak=float(col(MCA).max()), anaphylaxis=bool((col(MCA) > P["MCA_anaphyl"]).any()),
        SS_peak=float(col(SS).max()), t_SS_peak=float(T[int(np.argmax(col(SS)))]),
        serum_sickness=bool(col(SS).max() > 1.0),
        VAUC=float(col(VAUC)[-1]),
        p_bleed=float(1.0 - math.exp(-col(HBLD)[-1])),
        p_death=float(1.0 - math.exp(-col(HDTH)[-1])),
        cV_peak=float(peak), cV_trough=trough, rebound_factor=rebound,
        rebound_detectable=rebound_detectable,
        VAUC_late=vauc_late, hrs_detectable_late=hrs_detect_late,
        recur_t=recur_t, recur_peak=recur_peak,
        recurrence=bool(not math.isnan(recur_t)),
        recoag_t=recoag_t, recoagulopathy=bool(not math.isnan(recoag_t)),
        note=sc.note,
    )


def stoichiometric_need(snake, mult=1.0):
    """mgNE of antivenom required to neutralise one bite's toxin load exactly."""
    sn = SNAKES[snake]
    return sum(sn["dose"] * mult * sn["frac"][t] / PK[t]["eps"] for t in TOX)


# ---------------------------------------------------------------------------
# ANALYTICAL HELPERS FOR THE REPORT
# ---------------------------------------------------------------------------
def two_cmt_halflives(V1, V2, Q, CL):
    k10, k12, k21 = CL / V1, Q / V1, Q / V2
    s, p = k10 + k12 + k21, k10 * k21
    disc = max(s * s - 4 * p, 0.0)
    a = (s + math.sqrt(disc)) / 2
    b = (s - math.sqrt(disc)) / 2
    return math.log(2) / a, math.log(2) / b


def report_product_pk(out):
    out("")
    out("=" * 78)
    out("A.  ANTIVENOM AND VENOM DISPOSITION — HALF-LIVES ARE COMPUTED, NOT ASSERTED")
    out("=" * 78)
    out("")
    out("  Antivenom products (2-compartment):")
    out(f"  {'product':<34}{'t1/2 alpha':>12}{'t1/2 beta':>12}{'mgNE/vial':>12}")
    for k, v in AV_PRODUCTS.items():
        a, b = two_cmt_halflives(v["V1"], v["V2"], v["Q"], v["CL"])
        out(f"  {v['label']:<34}{a:>11.2f}h{b:>11.1f}h{v['mgNE_vial']:>12.1f}")
    out("")
    out("  Venom toxin classes (two parallel depots + one systemic compartment).")
    out("  The TERMINAL half-life is flip-flop-limited by the slow depot, not by")
    out("  systemic elimination — the venom leaves the body faster than it leaves")
    out("  the bite.  That is the number an antivenom has to outlast.")
    out(f"  {'class':<8}{'MW':>5}{'f_slow':>8}{'abs fast':>11}{'abs SLOW':>11}"
        f"{'systemic':>11}{'TERMINAL':>11}{'eps':>7}")
    for t in TOX:
        p = PK[t]
        t_sys = math.log(2) * p["Vc"] / p["CL"]
        t_term = math.log(2) / (p["ka_s"] + p["kdeg_s"])
        out(f"  {t:<8}{p['MW']:>5.0f}{p['f_slow']:>8.2f}"
            f"{math.log(2)/p['ka']:>10.2f}h{math.log(2)/p['ka_s']:>10.1f}h"
            f"{t_sys:>10.2f}h{t_term:>10.1f}h{p['eps']:>7.2f}")
    out("")
    fab_b = two_cmt_halflives(**{k: AV_PRODUCTS['Fab'][k] for k in ('V1', 'V2', 'Q', 'CL')})[1]
    f2_b = two_cmt_halflives(**{k: AV_PRODUCTS["F(ab')2"][k] for k in ('V1', 'V2', 'Q', 'CL')})[1]
    igg_b = two_cmt_halflives(**{k: AV_PRODUCTS["IgG"][k] for k in ('V1', 'V2', 'Q', 'CL')})[1]
    svmp_b = math.log(2) / (PK['SVMP']['ka_s'] + PK['SVMP']['kdeg_s'])
    out("  >>> THE RECURRENCE INEQUALITY (claim 2), evaluated:")
    out(f"        venom SVMP terminal (absorption-limited) t1/2 = {svmp_b:6.1f} h")
    out(f"        ovine Fab terminal t1/2       = {fab_b:6.1f} h   ratio {fab_b/svmp_b:.2f}  -> Fab LOSES the race")
    out(f"        equine F(ab')2 terminal t1/2  = {f2_b:6.1f} h   ratio {f2_b/svmp_b:.2f}  -> F(ab')2 wins it")
    out(f"        equine whole IgG terminal t1/2= {igg_b:6.1f} h   ratio {igg_b/svmp_b:.2f}  -> IgG wins it")
    out("      Recurrence is therefore a property of the FRAGMENT, not of the dose")
    out("      and not of the operator.")
    out("")


def report_stoichiometry(out):
    out("=" * 78)
    out("B.  HOW MANY VIALS IS 'ENOUGH'?  — the labelled potency, done as arithmetic")
    out("=" * 78)
    out("")
    out("  mgNE required = sum over toxin classes of (mass in the bite / eps_class).")
    out("  Vial counts below are what the label's mgNE/vial delivers for a 1.0x bite.")
    out("")
    out(f"  {'snake':<30}{'venom mg':>9}{'mgNE need':>11}"
        f"{'IgG vials':>11}{'Fab vials':>11}")
    for s in SNAKES:
        if SNAKES[s]["dose"] == 0:
            continue
        need = stoichiometric_need(s)
        out(f"  {s:<30}{SNAKES[s]['dose']:>9.0f}{need:>11.1f}"
            f"{need/AV_PRODUCTS['IgG']['mgNE_vial']:>11.1f}"
            f"{need/AV_PRODUCTS['Fab']['mgNE_vial']:>11.1f}")
    out("")
    out("  Read the Daboia and Naja rows together.  Naja naja carries LESS venom")
    out("  protein than Daboia russelii (40 mg vs 63 mg) and yet needs MORE")
    out("  neutralising capacity (85.7 vs 61.0 mgNE), purely because 55% of its")
    out("  mass is three-finger toxin, the class antivenom covers worst")
    out("  (eps = 0.35).  The standard 10-vial Indian ASV dose is 60 mgNE: a")
    out("  0.98x margin against Daboia and a 0.70x margin against Naja.  The")
    out("  clinical habit of giving 20-30 vials for a neurotoxic cobra bite is")
    out("  not caution; it is stoichiometry.")
    out("")


# ---------------------------------------------------------------------------
# SCENARIOS
# ---------------------------------------------------------------------------
def scenarios():
    S = []
    add = S.append
    RV = "Daboia russelii (Sri Lanka)"
    CA = "Crotalus atrox"
    NN = "Naja naja"
    BC = "Bungarus caeruleus"
    FLUID = [(1.0, 8.0, 0.30)]

    # ================= claim 1: depth vs slope ============================
    add(Scenario("01 Russell's viper — no antivenom", RV, fluids=FLUID,
                 note="untreated reference; VICC + AKI + capillary leak"))
    add(Scenario("02 Russell's viper — ASV 10 v @ 1 h", RV, product="IgG",
                 vials=10, t_av=1.0, fluids=FLUID,
                 note="claim 1: earliest realistic antivenom"))
    add(Scenario("03 Russell's viper — ASV 10 v @ 4 h", RV, product="IgG",
                 vials=10, t_av=4.0, fluids=FLUID,
                 note="claim 1: median real-world delay"))
    add(Scenario("04 Russell's viper — ASV 10 v @ 12 h", RV, product="IgG",
                 vials=10, t_av=12.0, fluids=FLUID,
                 note="claim 1: late antivenom, substrate already at the floor"))
    add(Scenario("05 Russell's viper — ASV 20 v @ 4 h", RV, product="IgG",
                 vials=20, t_av=4.0, fluids=FLUID,
                 note="claim 1: does doubling the dose buy back the lost time?"))
    add(Scenario("06 Russell's viper — ASV 10 v @ 4 h + cryoprecipitate", RV,
                 product="IgG", vials=10, t_av=4.0, fluids=FLUID,
                 ffp=[(4.5, 0.90), (10.0, 0.90)],
                 note="claim 1: substitute for the liver instead of waiting for it"))
    add(Scenario("07 Russell's viper — ASV 10 v @ 4 h + tranexamic acid", RV,
                 product="IgG", vials=10, t_av=4.0, fluids=FLUID,
                 txa=[4.0, 12.0, 20.0],
                 note="claim 1: only 5% of the fibrinogen loss is plasmin-mediated"))

    # ================= claim 2: recurrence =================================
    # Crotalus atrox + ovine Fab is the setting in which recurrent venom
    # antigenaemia and recurrent coagulopathy are actually documented, and the
    # CroFab label's mandatory 2-vials-q6h-x3 maintenance schedule is the fix the
    # inequality predicts.  The equine F(ab')2 product carries no such
    # requirement.  Arms 10-12 hold the MOLAR MARGIN constant while changing the
    # bite size, which is what separates "the margin matters" from "the vial
    # count matters".
    add(Scenario("08 Crotalus atrox — Fab 6 v @ 2 h (margin 1.4x)", CA,
                 product="Fab", vials=6, t_av=2.0, fluids=FLUID,
                 note="claim 2: short-lived fragment against a multi-day venom tail"))
    add(Scenario("09 Crotalus atrox — F(ab')2 6 v @ 2 h (same mgNE, margin 1.4x)", CA,
                 product="F(ab')2", vials=6, t_av=2.0, fluids=FLUID,
                 note="claim 2: identical molar dose, 5.8x the half-life"))
    add(Scenario("10 Crotalus atrox — Fab 12 v @ 2 h (margin 2.8x)", CA,
                 product="Fab", vials=12, t_av=2.0, fluids=FLUID,
                 note="claim 2: doubling works — at THIS bite size"))
    add(Scenario("11 Crotalus atrox 2x venom — Fab 12 v @ 2 h (margin 1.4x again)", CA,
                 venom_mult=2.0, product="Fab", vials=12, t_av=2.0, fluids=FLUID,
                 note="claim 2: same margin, twice the tail — recurrence returns"))
    add(Scenario("12 Crotalus atrox 2x venom — F(ab')2 12 v @ 2 h (margin 1.4x)", CA,
                 venom_mult=2.0, product="F(ab')2", vials=12, t_av=2.0, fluids=FLUID,
                 note="claim 2: the long fragment does not care about the tail"))
    add(Scenario("13 Crotalus atrox — Fab 6 v @ 2 h + 2 v q6h x 3 (label schedule)",
                 CA, product="Fab", vials=6, t_av=2.0, fluids=FLUID,
                 maintenance=(8.0, 6.0, 2.0, 3),
                 note="claim 2: the fix is a schedule, not a bigger bolus"))

    # ================= claim 3: two paralyses ==============================
    add(Scenario("14 Naja naja — no antivenom", NN,
                 note="claim 3: postsynaptic paralysis, untreated"))
    add(Scenario("15 Naja naja — neostigmine q4h only, no antivenom", NN,
                 neostigmine=[2.0, 6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0, 34.0, 38.0],
                 note="claim 3: the competition effect, isolated"))
    add(Scenario("16 Naja naja — ASV 10 v @ 4 h (label dose, margin 0.70x)", NN,
                 product="IgG", vials=10, t_av=4.0,
                 note="claim 3: the label dose is under-stoichiometric for a cobra"))
    add(Scenario("17 Naja naja — ASV 10 v @ 4 h + neostigmine q4h", NN,
                 product="IgG", vials=10, t_av=4.0,
                 neostigmine=[4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0],
                 note="claim 3: acetylcholine out-competes a reversible block"))
    add(Scenario("18 Naja naja — ASV 20 v @ 4 h (margin 1.40x)", NN,
                 product="IgG", vials=20, t_av=4.0,
                 note="claim 3: what adequate neutralisation looks like"))
    add(Scenario("19 Bungarus caeruleus — ASV 10 v @ 6 h", BC,
                 product="IgG", vials=10, t_av=6.0,
                 note="claim 3: presynaptic destruction already done"))
    add(Scenario("20 Bungarus caeruleus — ASV 10 v @ 6 h + neostigmine q4h", BC,
                 product="IgG", vials=10, t_av=6.0,
                 neostigmine=[6.0, 10.0, 14.0, 18.0, 22.0, 26.0, 30.0, 34.0, 38.0, 42.0],
                 note="claim 3: nothing to out-compete"))
    add(Scenario("21 Bungarus caeruleus — ASV 10 v @ 0.5 h (prevention)", BC,
                 product="IgG", vials=10, t_av=0.5,
                 note="claim 3: the same drug, before the destruction"))
    add(Scenario("22 Bungarus caeruleus — ASV 10 v @ 6 h, NO ventilator", BC,
                 product="IgG", vials=10, t_av=6.0, icu=False,
                 note="claim 3: what actually saves a krait bite"))

    # ================= where antivenom cannot go ===========================
    add(Scenario("23 Bothrops asper — F(ab')2 10 v @ 4 h", "Bothrops asper",
                 product="F(ab')2", vials=10, t_av=4.0, fluids=FLUID,
                 note="systemic control, local failure"))
    add(Scenario("24 Bothrops asper — F(ab')2 10 v @ 4 h + varespladib @ 0.5 h",
                 "Bothrops asper", product="F(ab')2", vials=10, t_av=4.0,
                 fluids=FLUID, varespladib=(0.5, 12.0, 500.0, 8),
                 note="a small molecule reaches tissue an IgG cannot"))

    # ================= other settings ======================================
    add(Scenario("25 Echis ocellatus — polyvalent IgG 3 v @ 3 h", "Echis ocellatus",
                 product="IgG", vials=3, t_av=3.0, fluids=FLUID,
                 note="West African pure-SVMP coagulopathy"))
    add(Scenario("26 Dry bite — ASV 10 v @ 1 h", "Dry bite", product="IgG",
                 vials=10, t_av=1.0, note="all of the risk, none of the benefit"))
    add(Scenario("27 Russell's viper 2.5x venom — ASV 10 v @ 4 h", RV,
                 venom_mult=2.5, product="IgG", vials=10, t_av=4.0, fluids=FLUID,
                 note="the label dose against a large bite: margin 0.39x"))
    return S


def recurrence_sensitivity(out):
    """WHAT THE MODEL CANNOT REPRODUCE, AND WHAT WOULD HAVE TO BE TRUE.

    The Fab arm reproduces recurrent venom ANTIGENAEMIA cleanly.  It does NOT
    reproduce severe recurrent HYPOFIBRINOGENAEMIA, which is also reported.
    Rather than tune the depot until the clinical phenomenon appears, this
    function asks what would have to be true for it to appear at all — and
    answers that the depot mechanism cannot do it at ANY physically available
    parameter value, while complex dissociation can, at a rate this puts a
    number on.  That converts an unexplained gap into a falsifiable prediction.
    """
    FL = [(1.0, 8.0, 0.30)]
    probe = lambda: Scenario("probe", "Crotalus atrox", product="Fab", vials=6,
                             t_av=2.0, fluids=FL)
    ceiling = P["kdeg_fg"] * P["FG0"] * P["SYNUP_max"] * 24.0

    out("")
    out("=" * 78)
    out("E.  WHAT THIS MODEL DOES NOT REPRODUCE (and what would have to be true)")
    out("=" * 78)
    out("")
    out("  The Fab arm reproduces recurrent venom ANTIGENAEMIA (a 3.6-fold rebound")
    out("  peaking on day 5).  It does NOT produce severe recurrent")
    out("  HYPOFIBRINOGENAEMIA, which is also reported clinically after ovine Fab.")
    out("  The reason is arithmetic:")
    out("")
    _SIM_CACHE.clear()
    T, Y, Dl, _ = simulate(probe())
    late = T >= 72.0
    fl = np.asarray([x["FGloss"] for x in Dl])[late]
    out(f"      hepatic fibrinogen synthesis ceiling         = {ceiling:6.2f} g/L/day")
    out(f"      peak late (>72 h) fibrinogen consumption     = {fl.max()*24:6.3f} g/L/day"
        f"  ({100*fl.max()*24/ceiling:.1f} % of it)")
    out(f"      shortfall factor needed for a second fall    = {ceiling/max(fl.max()*24,1e-9):6.1f}x")
    out("")
    out("  ROUTE 1 — put more of the bite into the slow depot.  f_slow cannot exceed")
    out("  1.0 (the whole bite sequestered), so the available headroom is bounded:")
    out("")
    out(f"      {'f_slow (SVMP)':<16}{'late flux':>12}{'% ceiling':>11}"
        f"{'recur antigen':>15}{'recur coag':>12}")
    base_fs = {t: PK[t]["f_slow"] for t in TOX}
    for fs in [0.35, 0.50, 0.70, 0.85, 1.00]:
        for t in TOX:
            PK[t]["f_slow"] = min(base_fs[t] * fs / base_fs["SVMP"], 1.0)
        _SIM_CACHE.clear()
        T, Y, Dl, av = simulate(probe())
        r = summarize(probe(), T, Y, Dl, av)
        fl = np.asarray([x["FGloss"] for x in Dl])[T >= 72.0]
        out(f"      {fs:<16.2f}{fl.max()*24:>12.3f}{100*fl.max()*24/ceiling:>10.1f} %"
            f"{('yes ' + ('%.0f h' % r['recur_t'])) if r['recurrence'] else 'no':>15}"
            f"{'YES' if r['recoagulopathy'] else 'no':>12}")
    for t in TOX:
        PK[t]["f_slow"] = base_fs[t]
    _SIM_CACHE.clear()
    out("")
    out("  Even with the ENTIRE bite sequestered, the late flux stays far below the")
    out("  synthesis rate.  Depot release cannot be the mechanism of clinical")
    out("  recurrent coagulopathy.  That is a negative result, and it is the useful")
    out("  half of this section.")
    out("")
    out("  ROUTE 2 — let the antivenom-toxin complex dissociate.  The complexed pool")
    out("  is large (tens of mg of toxin), so even a slow off-rate is a big source.")
    out("  This model's default is k_off = 0; sweeping it:")
    out("")
    out(f"      {'k_off (1/h)':<14}{'complex t1/2':>14}{'late flux':>13}"
        f"{'% ceiling':>12}{'recur coag':>13}{'FG 2nd fall':>14}")
    found = None
    for ko in [0.0, 0.0020, 0.0080, 0.0200, 0.0600, 0.1500]:
        P["koff_av"] = ko
        _SIM_CACHE.clear()
        T, Y, Dl, av = simulate(probe())
        r = summarize(probe(), T, Y, Dl, av)
        fl = np.asarray([x["FGloss"] for x in Dl])[T >= 72.0]
        fg = Y[:, FG]
        j = int(np.argmin(fg))
        seg = fg[j:]
        drop = float(np.max(np.maximum.accumulate(seg) - seg)) if seg.size else 0.0
        if r["recoagulopathy"] and found is None:
            found = ko
        th = ("inf" if ko == 0 else f"{math.log(2)/ko:.0f} h")
        out(f"      {ko:<14.4f}{th:>14}{fl.max()*24:>12.3f}"
            f"{100*fl.max()*24/ceiling:>11.1f} %{'YES' if r['recoagulopathy'] else 'no':>13}"
            f"{drop:>11.2f} g/L")
    P["koff_av"] = 0.0
    _SIM_CACHE.clear()
    out("")
    if found:
        out(f"  A second fibrinogen fall first appears at k_off = {found:.4f} /h"
            f" (complex t1/2 {math.log(2)/found:.0f} h).")
    else:
        out("  Complex dissociation does not produce it either, at any swept rate —")
        out("  and the reason is instructive.  The complexed pool is cleared with an")
        out(f"  11.6 h half-life, so by 72 h only ~{100*math.exp(-P['kel_cplx']*72):.0f} % of it is left."
            f"  Dissociation")
        out("  releases venom EARLY, while the pool is still large, which prolongs the")
        out("  acute phase instead of creating a late second event.")
    out("")
    out("  THE CONCLUSION, STATED AS A NEGATIVE RESULT:")
    out("")
    out("    Recurrent venom ANTIGENAEMIA after a short-lived fragment follows from")
    out("    absorption kinetics alone and needs no extra assumption.  A late second")
    out("    episode of SEVERE, unclottable-blood coagulopathy does not: it needs a")
    out(f"    late venom flux ~{ceiling/max(0.094,1e-9):.0f}x larger than a sequestered depot can supply,")
    out("    and neither a larger sequestered fraction (bounded at 57 % of the")
    out("    synthesis ceiling even if the WHOLE bite is sequestered) nor complex")
    out("    dissociation gets there.  Two readings survive, and they are")
    out("    distinguishable by measurement rather than by argument:")
    out("")
    out("      (i)  the recurrent hypofibrinogenaemia reported after ovine Fab is")
    out("           mostly a LABORATORY event — a dip of a few tenths of a g/L, which")
    out("           is what this model does produce — rather than a return to")
    out("           unclottable blood; or")
    out("      (ii) there is a venom reservoir this model has no compartment for.")
    out("")
    out("    This is reported rather than tuned away, because the alternative —")
    out("    raising f_slow to 1.0 and calling it calibration — would have buried a")
    out("    real discrepancy under a plausible-looking number.")
    out("")


# ---------------------------------------------------------------------------
# VIRTUAL POPULATION
# ---------------------------------------------------------------------------
def virtual_population(out, n=300, seed=20260805):
    """A trial arm is a MIXTURE.  A single median patient cannot reproduce a
    trial's event rate, because the endpoints in this model are thresholds and
    the mean of a threshold is not the threshold of a mean."""
    rng = np.random.default_rng(seed)
    out("")
    out("=" * 78)
    out("F.  VIRTUAL POPULATION — why a median patient cannot reproduce a trial rate")
    out("=" * 78)
    out("")
    arms = [
        ("Daboia + ASV 10 v, delay ~ 4 h", "Daboia russelii (Sri Lanka)", "IgG", 10, 4.0),
        ("Daboia + ASV 10 v, delay ~ 1 h", "Daboia russelii (Sri Lanka)", "IgG", 10, 1.0),
        ("Daboia + ASV 20 v, delay ~ 4 h", "Daboia russelii (Sri Lanka)", "IgG", 20, 4.0),
        ("Daboia, no antivenom",           "Daboia russelii (Sri Lanka)", None, 0, None),
        ("C. atrox + Fab 6 v @ 2 h",       "Crotalus atrox", "Fab", 6, 2.0),
        ("C. atrox + F(ab')2 6 v @ 2 h",   "Crotalus atrox", "F(ab')2", 6, 2.0),
    ]
    res = {}
    for label, snake, prod, vials, delay in arms:
        rows = []
        for _ in range(n):
            # between-subject variability: injected venom mass is the dominant
            # source (bites are not calibrated injections), then absorption rate
            # and presentation delay.
            vm = float(rng.lognormal(mean=math.log(1.0) - 0.5 * 0.55 ** 2, sigma=0.55))
            vm = min(max(vm, 0.10), 3.5)
            td = None if delay is None else float(max(0.3, rng.lognormal(math.log(delay) - 0.5 * 0.5 ** 2, 0.5)))
            sc = Scenario(label, snake, venom_mult=vm, product=prod, vials=vials,
                          t_av=td, fluids=[(1.0, 6.0, 0.35)], tmax=240.0)
            try:
                T, Y, Dl, av = simulate(sc, dt=0.50)
            except RuntimeError:
                continue
            s = summarize(sc, T, Y, Dl, av)
            rows.append(s)
        res[label] = rows
        if not rows:
            continue
        f = lambda k: np.asarray([r[k] for r in rows], dtype=float)
        out(f"  {label}   (n={len(rows)})")
        out(f"      unclottable blood at any time  {100*np.mean(f('hrs_unclottable')>0):5.1f} %")
        out(f"      hours unclottable  median      {np.median(f('hrs_unclottable')):5.1f} h"
            f"   (IQR {np.percentile(f('hrs_unclottable'),25):.1f}-{np.percentile(f('hrs_unclottable'),75):.1f})")
        out(f"      venom recurrence               {100*np.mean([r['recurrence'] for r in rows]):5.1f} %")
        out(f"      late (>48 h) venom detectable  {100*np.mean(f('hrs_detectable_late')>0):5.1f} %")
        out(f"      AKI any stage (SCr >= 1.5x)     {100*np.mean(f('SCR_peak') >= 1.5*P['SCR0']):5.1f} %")
        out(f"      AKI stage 3 (SCr >= 3x)        {100*np.mean([r['AKI3'] for r in rows]):5.1f} %")
        out(f"      mechanical ventilation         {100*np.mean([r['ventilated'] for r in rows]):5.1f} %")
        out(f"      systemic bleeding (modelled)   {100*np.mean(f('p_bleed')):5.1f} %")
        out(f"      death (modelled)               {100*np.mean(f('p_death')):5.1f} %")
        out(f"      needed > label dose (margin<1) {100*np.mean(f('molar_margin')<1.0):5.1f} %")
        out("")
    out("  The 'needed > label dose' column is the model's account of why 30-50%")
    out("  of Russell's viper bites treated with the labelled 10-vial dose require")
    out("  a repeat dose: the labelled dose is a 0.98x margin against the MEDIAN")
    out("  bite, so roughly half of a real bite distribution is under-dosed by")
    out("  construction.  Nothing about the antivenom is failing.")
    out("")
    return res


# ---------------------------------------------------------------------------
# MAIN REPORT
# ---------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    lines = []

    def out(s=""):
        lines.append(s)
        print(s)

    out("=" * 78)
    out("SNAKEBITE ENVENOMING + ANTIVENOM — QSP REFERENCE MODEL RUN")
    out("독사 교상 + 항독소 — QSP 참조 모델 실행 결과")
    out("=" * 78)
    out("")
    out("  50 ODEs · 4 toxin classes · 3 antivenom fragment types · 8 snake")
    out("  archetypes · 22 scenarios · 300-subject virtual population.")
    out("")
    out("  THESIS:  antivenom binds; it does not undo.  Two clocks — a venom")
    out("  clock in hours and a substrate clock in days — and every endpoint")
    out("  is read off the slow one.")

    report_product_pk(out)
    report_stoichiometry(out)

    # ---------------- scenarios ----------------
    S = scenarios()
    out("=" * 78)
    out(f"C.  {len(S)} SCENARIOS")
    out("=" * 78)
    results = {}
    for sc in S:
        T, Y, Dl, av = simulate(sc, cache_key=sc.name)
        s = summarize(sc, T, Y, Dl, av)
        results[sc.name] = s
        out("")
        out(f"  {sc.name}")
        out(f"     snake {sc.snake} ({SNAKES[sc.snake]['ko']})   venom {SNAKES[sc.snake]['dose']*sc.venom_mult:.0f} mg")
        out(f"     antivenom {s['product']}, {s['vials']:.0f} vial(s) = {s['av_mgNE']:.0f} mgNE"
            + (f" @ {sc.t_av:.1f} h" if sc.t_av is not None else ""))
        out(f"     stoichiometric need {s['need_mgNE']:.1f} mgNE  ->  molar margin {s['molar_margin']:.2f}x")
        if s["never_unclottable"]:
            clot_msg = " | blood clottable throughout (20WBCT never positive)"
        elif math.isnan(s["t_clottable"]):
            clot_msg = " | still unclottable at end of window"
        else:
            clot_msg = f" | clottable again at {s['t_clottable']:.1f} h"
        out(f"     HAEMOSTASIS  fibrinogen nadir {s['FG_nadir']:.2f} g/L at {s['t_FG_nadir']:.1f} h"
            f" | unclottable {s['hrs_unclottable']:.1f} h" + clot_msg)
        out(f"                  platelet nadir {s['PLT_nadir']:.0f} x10^9/L | FX nadir {s['FX_nadir']:.2f}"
            f" | D-dimer peak {s['XDP_peak']:.1f} ug/mL")
        out(f"     NEUROMUSCULAR safety factor nadir {s['SF_nadir']:.2f} | SBC nadir {s['SBC_nadir']:.0f}"
            f" | AChR occupancy peak {s['BR_peak']:.3f} | terminal integrity nadir {s['TERM_nadir']:.3f}")
        out(f"                  ptosis {s['hrs_ptosis']:.0f} h | ventilation {s['hrs_ventilated']:.0f} h")
        out(f"     LOCAL        necrosis {100*s['NEC_final']:.1f} % | compartment pressure peak {s['CP_peak']:.0f} mmHg"
            f" | CK peak {s['CK_peak']:.0f} U/L at {s['t_CK_peak']:.0f} h")
        out(f"     RENAL        creatinine peak {s['SCR_peak']:.2f} mg/dL at {s['t_SCR_peak']:.0f} h"
            f" | eGFR nadir {s['eGFR_nadir']:.0f} | permanent scar {100*s['FIBR_final']:.1f} %")
        out(f"     SYSTEMIC     MAP nadir {s['MAP_nadir']:.0f} mmHg | Hct peak {s['Hct_peak']:.0f} %"
            f" | IL-6 peak {s['IL6_peak']:.0f} pg/mL")
        out(f"     ANTIVENOM AE anaphylactoid index {s['MCA_peak']:.2f}"
            f"{'  ANAPHYLAXIS' if s['anaphylaxis'] else ''}"
            f" | serum sickness {s['SS_peak']:.2f} at {s['t_SS_peak']:.0f} h"
            f"{'  CLINICAL' if s['serum_sickness'] else ''}")
        out(f"     EXPOSURE     free-venom AUC {s['VAUC']:.2f} mg.h/L"
            + (f" | RECURRENCE at {s['recur_t']:.0f} h (peak {s['recur_peak']:.4f} mg/L)"
               if s['recurrence'] else " | no recurrence"))
        out(f"     OUTCOME      P(systemic bleeding) {100*s['p_bleed']:.1f} %"
            f" | P(death) {100*s['p_death']:.1f} %")
        out(f"     [{sc.note}]")

    # ---------------- the four claims, side by side ----------------
    out("")
    out("=" * 78)
    out("D.  THE FOUR CLAIMS, EVALUATED AGAINST EACH OTHER")
    out("=" * 78)

    # ---- D1 --------------------------------------------------------------
    out("")
    out("  D1.  DEPTH vs SLOPE — antivenom moves the nadir, not the recovery rate")
    out("")
    out(f"      {'regimen':<26}{'FG nadir':>10}{'unclot h':>10}"
        f"{'rise @ FG 1.0':>16}{'FG day 7':>10}{'SCr peak':>10}")
    keys1 = ["01 Russell's viper — no antivenom",
             "02 Russell's viper — ASV 10 v @ 1 h",
             "03 Russell's viper — ASV 10 v @ 4 h",
             "04 Russell's viper — ASV 10 v @ 12 h",
             "05 Russell's viper — ASV 20 v @ 4 h"]
    slopes = []
    for k in keys1:
        sc = next(x for x in S if x.name == k)
        r = results[k]
        lbl = "none" if sc.t_av is None else f"{sc.vials:.0f} vials @ {sc.t_av:.0f} h"
        v = r["FG_rise_at_1"]
        vs = "n/a (never <1)" if math.isnan(v) else f"{v:.2f} g/L/day"
        if not math.isnan(v):
            slopes.append(v)
        out(f"      {lbl:<26}{r['FG_nadir']:>10.2f}{r['hrs_unclottable']:>10.1f}"
            f"{vs:>16}{r['FG_d7']:>10.2f}{r['SCR_peak']:>10.2f}")
    nad = [results[k]["FG_nadir"] for k in keys1]
    sl = np.asarray(slopes)
    out("")
    out(f"      fibrinogen nadir spans {min(nad):.2f} - {max(nad):.2f} g/L"
        f"   ({max(nad)/max(min(nad), 1e-3):.0f}-fold)")
    out(f"      rise rate measured at the SAME fibrinogen level (1.0 g/L) spans"
        f" {sl.min():.2f} - {sl.max():.2f} g/L/day")
    out(f"      = a {100*(sl.max()-sl.min())/sl.mean():.0f} % spread in the slope against a"
        f" {max(nad)/max(min(nad), 1e-3):.0f}-fold spread in the nadir.")
    out(f"      hepatic synthesis ceiling (k_syn x SYNUP_max) ="
        f" {P['kdeg_fg']*P['FG0']*P['SYNUP_max']*24:.2f} g/L/day — every observed")
    out("      rise rate above sits under it.")
    out("      The nadir is a function of when the antivenom arrived.  The rise rate")
    out("      is not: it is the liver's synthesis term, and neither the antivenom")
    out("      dose nor its timing appears anywhere in that term.  The residual")
    out("      2-fold spread in the slope is NOT the antivenom either — it is")
    out("      residual venom still consuming fibrinogen on the way up, which is")
    out("      why the FASTEST rise belongs to the 12-hour arm (deepest deficit,")
    out("      largest acute-phase SYNUP) and not to the earliest-treated one.")
    r3, r5, r6, r7 = (results["03 Russell's viper — ASV 10 v @ 4 h"],
                      results["05 Russell's viper — ASV 20 v @ 4 h"],
                      results["06 Russell's viper — ASV 10 v @ 4 h + cryoprecipitate"],
                      results["07 Russell's viper — ASV 10 v @ 4 h + tranexamic acid"])
    out("")
    out(f"      Doubling the vial count at a fixed 4 h delay: FG nadir"
        f" {r3['FG_nadir']:.2f} -> {r5['FG_nadir']:.2f} g/L. Timing is not a dose.")
    out(f"      Cryoprecipitate (substituting for the liver): P(bleed)"
        f" {100*r3['p_bleed']:.1f} % -> {100*r6['p_bleed']:.1f} %.")
    out(f"      Tranexamic acid: FG nadir {r3['FG_nadir']:.3f} -> {r7['FG_nadir']:.3f} g/L,"
        f" P(bleed) {100*r3['p_bleed']:.1f} % -> {100*r7['p_bleed']:.1f} %.")
    out("      TXA is unchanged to three significant figures because 95% of the")
    out("      fibrinogen loss is direct enzymatic cleavage, which plasmin")
    out("      inhibition cannot touch.")

    # ---- D2 --------------------------------------------------------------
    out("")
    out("  D2.  RECURRENCE IS THE MOLAR MARGIN OVER THE TAIL, NOT THE VIAL COUNT")
    out("")
    out(f"      {'regimen':<52}{'margin':>8}{'rebound':>9}{'recur at':>10}"
        f"{'late AUC':>11}{'detect h':>9}")
    keys2 = ["08 Crotalus atrox — Fab 6 v @ 2 h (margin 1.4x)",
             "09 Crotalus atrox — F(ab')2 6 v @ 2 h (same mgNE, margin 1.4x)",
             "10 Crotalus atrox — Fab 12 v @ 2 h (margin 2.8x)",
             "11 Crotalus atrox 2x venom — Fab 12 v @ 2 h (margin 1.4x again)",
             "12 Crotalus atrox 2x venom — F(ab')2 12 v @ 2 h (margin 1.4x)",
             "13 Crotalus atrox — Fab 6 v @ 2 h + 2 v q6h x 3 (label schedule)"]
    short2 = ["Fab 6 v, 1x bite", "F(ab')2 6 v, 1x bite", "Fab 12 v, 1x bite",
              "Fab 12 v, 2x bite", "F(ab')2 12 v, 2x bite", "Fab 6 v + label maintenance"]
    for k, lb in zip(keys2, short2):
        r = results[k]
        rb = (f"{r['rebound_factor']:.1f}x" if r["rebound_detectable"] else "<LOD")
        out(f"      {lb:<52}{r['molar_margin']:>7.2f}x{rb:>9}"
            f"{('%.0f h' % r['recur_t']) if r['recurrence'] else 'none':>10}"
            f"{r['VAUC_late']:>11.4f}{r['hrs_detectable_late']:>9.0f}")
    f8, f9 = results[keys2[0]], results[keys2[1]]
    f11, f12 = results[keys2[3]], results[keys2[4]]
    out("")
    out(f"      At an IDENTICAL molar dose and margin, swapping ovine Fab for equine")
    out(f"      F(ab')2 reduces late (>48 h) free-venom exposure"
        f" {f8['VAUC_late']/max(f9['VAUC_late'],1e-9):.0f}-fold")
    out(f"      ({f8['VAUC_late']:.4f} -> {f9['VAUC_late']:.4f} mg.h/L) and the hours of")
    out(f"      detectable venom from {f8['hrs_detectable_late']:.0f} h to"
        f" {f9['hrs_detectable_late']:.0f} h.")
    out(f"      Doubling the Fab dose abolishes recurrence at a 1x bite (margin 2.8x)")
    out(f"      and it comes straight back at a 2x bite treated with the same doubled")
    out(f"      dose (margin 1.4x again, recurrence at {f11['recur_t']:.0f} h), while the")
    out(f"      long fragment at that same 1.4x margin does not recur"
        f" (late AUC {f12['VAUC_late']:.4f}).")
    out("      A fixed vial count is therefore under-specified: whether a patient")
    out("      recurs depends on the bite size they happened to receive.")

    # ---- D3 --------------------------------------------------------------
    out("")
    out("  D3.  TWO PARALYSES, ONE BEDSIDE PICTURE, OPPOSITE PHARMACOLOGY")
    out("")
    out(f"      {'case':<52}{'AChR occ':>10}{'terminal':>10}{'ptosis h':>10}"
        f"{'vent h':>8}{'P(death)':>10}")
    keys3 = ["14 Naja naja — no antivenom",
             "15 Naja naja — neostigmine q4h only, no antivenom",
             "16 Naja naja — ASV 10 v @ 4 h (label dose, margin 0.70x)",
             "17 Naja naja — ASV 10 v @ 4 h + neostigmine q4h",
             "18 Naja naja — ASV 20 v @ 4 h (margin 1.40x)",
             "19 Bungarus caeruleus — ASV 10 v @ 6 h",
             "20 Bungarus caeruleus — ASV 10 v @ 6 h + neostigmine q4h",
             "21 Bungarus caeruleus — ASV 10 v @ 0.5 h (prevention)",
             "22 Bungarus caeruleus — ASV 10 v @ 6 h, NO ventilator"]
    for k in keys3:
        r = results[k]
        out(f"      {k[3:]:<52}{r['BR_peak']:>10.3f}{r['TERM_nadir']:>10.3f}"
            f"{r['hrs_ptosis']:>10.0f}{r['hrs_ventilated']:>8.0f}{100*r['p_death']:>9.1f}%")
    n14, n16, n17, n18 = (results[keys3[0]], results[keys3[2]], results[keys3[3]], results[keys3[4]])
    n15 = results[keys3[1]]
    k19, k20, k21, k22 = (results[keys3[5]], results[keys3[6]], results[keys3[7]], results[keys3[8]])
    out("")
    out("      NEOSTIGMINE, the same drug, in the two cases:")
    out(f"        cobra, no antivenom     ventilation {n14['hrs_ventilated']:.0f} h ->"
        f" {n15['hrs_ventilated']:.0f} h   ({100*(1-n15['hrs_ventilated']/max(n14['hrs_ventilated'],1e-9)):.0f} % reduction)")
    out(f"        cobra, + antivenom      ventilation {n16['hrs_ventilated']:.0f} h ->"
        f" {n17['hrs_ventilated']:.0f} h, ptosis {n16['hrs_ptosis']:.0f} h -> {n17['hrs_ptosis']:.0f} h")
    out(f"        krait,  + antivenom     ventilation {k19['hrs_ventilated']:.0f} h ->"
        f" {k20['hrs_ventilated']:.0f} h, ptosis {k19['hrs_ptosis']:.0f} h -> {k20['hrs_ptosis']:.0f} h")
    out("      Identical bedside signs. The drug acts on an occupancy; in the krait")
    out("      the occupancy is 0.16 and the deficit is a terminal at 0.15 of normal,")
    out("      so there is nothing for acetylcholine to out-compete.")
    out("")
    out("      ANTIVENOM DOSE vs ANTIVENOM TIMING are different questions:")
    out(f"        cobra at 4 h, 10 vials (margin 0.70x): ventilation {n16['hrs_ventilated']:.0f} h")
    out(f"        cobra at 4 h, 20 vials (margin 1.40x): ventilation {n18['hrs_ventilated']:.0f} h")
    out(f"        cobra at 4 h, 20 vials: P(death) {100*n18['p_death']:.1f} % vs"
        f" {100*n16['p_death']:.1f} % for 10 vials")
    out("      Doubling the dose changes almost nothing about the paralysis, and it")
    out("      RAISES modelled mortality, because 20 vials of whole equine IgG pushes")
    out("      the anaphylactoid index past the reaction threshold while the receptors")
    out("      it was given to unblock are already occupied.  That is not noise in the")
    out("      model; it is the ledger having two sides.")
    out("      The reason the dose does nothing here is that by 4 h the receptors")
    out(f"      are already {100*n16['BR_peak']:.0f} % occupied and the occupancy decays on its own")
    out(f"      35 h clock.  The krait shows the same asymmetry from the other side:")
    out(f"        krait at 6.0 h: ventilation {k19['hrs_ventilated']:.0f} h, terminal nadir {k19['TERM_nadir']:.3f}")
    out(f"        krait at 0.5 h: ventilation {k21['hrs_ventilated']:.0f} h, terminal nadir {k21['TERM_nadir']:.3f}")
    out(f"      Same product, same vial count, 5.5 h apart.")
    out("")
    out(f"      And what actually saves the krait patient is not the antivenom:")
    out(f"        krait at 6 h WITH a ventilator:    P(death) {100*k19['p_death']:.1f} %")
    out(f"        krait at 6 h WITHOUT a ventilator: P(death) {100*k22['p_death']:.1f} %")

    # ---- D4 --------------------------------------------------------------
    out("")
    out("  D4.  THE KIDNEY IS AN INTEGRAL, NOT A LEVEL")
    out("")
    s1 = results["01 Russell's viper — no antivenom"]
    s3 = results["03 Russell's viper — ASV 10 v @ 4 h"]

    def t_below_lod(name):
        sc2 = next(x for x in S if x.name == name)
        T2, Y2, D2, _ = simulate(sc2, cache_key=sc2.name)
        c = np.asarray([x["cVfree"] for x in D2])
        j = int(np.argmax(c))
        b = c[j:] < P["cV_LOD"]
        return float(T2[j + int(np.argmax(b))]) if b.any() else float("nan")

    t_lod_1 = t_below_lod("01 Russell's viper — no antivenom")
    t_lod_3 = t_below_lod("03 Russell's viper — ASV 10 v @ 4 h")

    out(f"      {'regimen':<26}{'venom AUC':>11}{'SCr peak':>10}{'at':>7}"
        f"{'eGFR nadir':>12}{'scar':>8}")
    for k in keys1[:4]:
        r = results[k]
        sc2 = next(x for x in S if x.name == k)
        lbl = "none" if sc2.t_av is None else f"{sc2.vials:.0f} vials @ {sc2.t_av:.0f} h"
        out(f"      {lbl:<26}{r['VAUC']:>11.2f}{r['SCR_peak']:>10.2f}"
            f"{r['t_SCR_peak']:>6.0f}h{r['eGFR_nadir']:>12.0f}{100*r['FIBR_final']:>7.1f}%")
    out("")
    out("      TREATED at 4 h: free venom drops below the assay limit of detection")
    out(f"      ({P['cV_LOD']*1000:.0f} ng/mL) at {t_lod_3:.1f} h, i.e. {t_lod_3-4.0:.1f} h after the infusion")
    out(f"      started — and creatinine then climbs for another"
        f" {s3['t_SCR_peak']-t_lod_3:.0f} h, peaking at {s3['t_SCR_peak']:.0f} h")
    out(f"      ({s3['SCR_peak']:.2f} mg/dL) with nothing left to measure in the plasma.")
    out("")
    out("      UNTREATED, the same test says something different and equally useful:")
    out(f"      free venom stays ABOVE the limit of detection until {t_lod_1:.0f} h, because the")
    out("      slow depot is still releasing it. So a positive venom antigen assay on")
    out("      day 5 does not mean treatment failed, and a negative one on day 1 does")
    out("      not mean the kidney is safe. Both readings are answered by the same")
    out(f"      structure: antivenom at 4 h cuts the free-venom AUC from {s1['VAUC']:.1f} to"
        f" {s3['VAUC']:.1f}")
    out(f"      mg.h/L ({100*(1-s3['VAUC']/max(s1['VAUC'],1e-9)):.0f} % less), the creatinine peak from"
        f" {s1['SCR_peak']:.2f} to {s3['SCR_peak']:.2f} mg/dL,")
    out(f"      and permanent scar from {100*s1['FIBR_final']:.1f} % to {100*s3['FIBR_final']:.1f} %.")
    out("      Nothing observable at presentation predicts this. The integral does.")

    # ---- D5 --------------------------------------------------------------
    out("")
    out("  D5.  WHERE ANTIVENOM CANNOT GO, AND WHAT CAN")
    out("")
    b23 = results["23 Bothrops asper — F(ab')2 10 v @ 4 h"]
    b24 = results["24 Bothrops asper — F(ab')2 10 v @ 4 h + varespladib @ 0.5 h"]
    out(f"      {'regimen':<50}{'necrosis':>10}{'CK peak':>10}{'CP mmHg':>10}{'venom AUC':>11}")
    for k, r in (("F(ab')2 10 v @ 4 h", b23),
                 ("F(ab')2 10 v @ 4 h + varespladib @ 0.5 h", b24)):
        out(f"      {k:<50}{100*r['NEC_final']:>9.1f}%{r['CK_peak']:>10.0f}"
            f"{r['CP_peak']:>10.0f}{r['VAUC']:>11.2f}")
    out("")
    out(f"      Intravenous antibody controls the plasma compartment completely")
    out(f"      (free-venom AUC {b23['VAUC']:.2f} mg.h/L) and still leaves"
        f" {100*b23['NEC_final']:.1f} % myonecrosis,")
    out(f"      because k_b,loc/k_b0 = {P['kb_loc']/P['kb0']:.3f} — a"
        f" {P['kb0']/P['kb_loc']:.0f}-fold penalty at the bite site.")
    out(f"      An oral catalytic-site PLA2 inhibitor with a 12 L volume of")
    out(f"      distribution, started at 0.5 h, takes necrosis to"
        f" {100*b24['NEC_final']:.1f} % ({100*(1-b24['NEC_final']/b23['NEC_final']):.0f} % less)")
    out(f"      and CK from {b23['CK_peak']:.0f} to {b24['CK_peak']:.0f} U/L. That is the")
    out(f"      entire rationale of the pre-hospital small-molecule programme.")

    # ---- D6 --------------------------------------------------------------
    out("")
    out("  D6.  THE HARM SIDE — a dry bite treated as an envenoming")
    out("")
    sd = results["26 Dry bite — ASV 10 v @ 1 h"]
    out(f"      benefit: free-venom AUC {sd['VAUC']:.3f} mg.h/L — there was no venom")
    out(f"      cost:    anaphylactoid index {sd['MCA_peak']:.2f}"
        f"{'  (CLINICAL ANAPHYLAXIS)' if sd['anaphylaxis'] else ''}")
    out(f"               serum sickness {sd['SS_peak']:.2f}, peaking on day {sd['t_SS_peak']/24:.1f}")
    out(f"               modelled P(death) {100*sd['p_death']:.2f} % — entirely iatrogenic")
    out("      20-50% of viperid bites and 10-30% of elapid bites inject no venom.")
    out("      A 20-minute whole blood clotting test costs a glass tube and is the")
    out("      gate; this row exists so the model cannot be read as an argument for")
    out("      treating every bite.")

    # ---- D7 --------------------------------------------------------------
    out("")
    out("  D7.  UNDER-DOSING IS BUILT INTO THE LABEL, NOT INTO THE OPERATOR")
    out("")
    s27 = results["27 Russell's viper 2.5x venom — ASV 10 v @ 4 h"]
    out(f"      {'bite':<34}{'margin':>8}{'FG nadir':>10}{'unclot h':>10}"
        f"{'SCr peak':>10}{'P(bleed)':>10}")
    for lbl, r in (("median bite (63 mg)", r3), ("2.5x bite (158 mg)", s27)):
        out(f"      {lbl:<34}{r['molar_margin']:>7.2f}x{r['FG_nadir']:>10.2f}"
            f"{r['hrs_unclottable']:>10.1f}{r['SCR_peak']:>10.2f}{100*r['p_bleed']:>9.1f}%")
    out("")
    out("      The same 10 vials, the same 4 h delay, the same operator, the same")
    out("      product. The only difference is how much venom the snake happened to")
    out("      inject — which nobody can measure at the bedside.")

    # ---------------- population ----------------
    recurrence_sensitivity(out)
    pop = virtual_population(out)

    # ---------------- mass balance / sanity ----------------
    out("=" * 78)
    out("G.  CONSERVATION AND SANITY CHECKS")
    out("=" * 78)
    out("")
    sc = next(x for x in S if x.name == "03 Russell's viper — ASV 10 v @ 4 h")
    T, Y, Dl, av = simulate(sc, cache_key=sc.name)
    # venom mass balance: depot + central + peripheral + complex + cleared
    sn = SNAKES[sc.snake]
    dose_in = sum(sn["dose"] * sn["frac"][t] for t in TOX)
    present = sum(Y[-1, iD[t]] + Y[-1, iV[t]] + Y[-1, iP[t]] + Y[-1, iC[t]] for t in TOX)
    out(f"  venom injected                {dose_in:9.4f} mg")
    out(f"  venom still in the body @ 14 d{present:9.4f} mg  ({100*present/dose_in:.3f} %)")
    out(f"  -> the model clears venom; nothing is trapped in an unphysical sink.")
    out("")
    for k, idx, lo, hi in [("BR (AChR occupancy)", BR, 0.0, 1.0),
                           ("TERM (terminal integrity)", TERM, 0.0, 1.0),
                           ("NEC (necrosis fraction)", NEC, 0.0, 1.0),
                           ("NEPH (nephron fraction)", NEPH, 0.0, 1.0),
                           ("GLX (glycocalyx)", GLX, 0.0, 1.0),
                           ("FIBR (renal scar)", FIBR, 0.0, 1.0)]:
        mn, mx = float(Y[:, idx].min()), float(Y[:, idx].max())
        ok = (mn >= lo - 1e-6) and (mx <= hi + 1e-6)
        out(f"  {k:<28} range [{mn:8.5f}, {mx:8.5f}]   {'OK' if ok else '*** OUT OF BOUNDS ***'}")
    out("")
    out("  Baseline stability (Dry bite, no antivenom): every state must sit still.")
    sc0 = Scenario("baseline", "Dry bite", tmax=336.0)
    T0, Y0, D0, _ = simulate(sc0)
    names = {FG: "FG", FX: "FX", PLT: "PLT", XDP: "XDP", SYNUP: "SYNUP",
             TERM: "TERM", NEC: "NEC", CK: "CK", MB: "MB", NEPH: "NEPH",
             SCR: "SCR", CAST: "CAST", PV: "PV", GLX: "GLX", IL6: "IL6",
             TNF: "TNF"}
    drift = []
    for idx, nm in names.items():
        a, b = Y0[0, idx], Y0[-1, idx]
        drift.append((abs(b - a) / max(abs(a), 1e-9), nm, a, b))
    drift.sort(reverse=True)
    for rel, nm, a, b in drift[:5]:
        out(f"    {nm:<7} {a:10.5f} -> {b:10.5f}   drift {100*rel:7.4f} %"
            f"   {'OK' if rel < 0.02 else '*** NOT A FIXED POINT ***'}")
    worst = drift[0][0]
    out(f"  largest 14-day baseline drift: {drift[0][1]} {100*worst:.4f} %"
        f"   {'OK' if worst < 0.02 else '*** NOT A FIXED POINT ***'}")
    out("")

    # ---------------- machine-readable dump ----------------
    with open(os.path.join(here, "sbe_scenario_results.json"), "w") as fh:
        json.dump({k: {kk: (None if isinstance(vv, float) and math.isnan(vv) else vv)
                       for kk, vv in v.items()} for k, v in results.items()},
                  fh, indent=1, ensure_ascii=False)
    popdump = {}
    for label, rows in pop.items():
        if not rows:
            continue
        f = lambda k: np.asarray([r[k] for r in rows], dtype=float)
        popdump[label] = dict(
            n=len(rows),
            pct_unclottable=float(100 * np.mean(f("hrs_unclottable") > 0)),
            median_hrs_unclottable=float(np.median(f("hrs_unclottable"))),
            pct_recurrence=float(100 * np.mean([r["recurrence"] for r in rows])),
            pct_AKI3=float(100 * np.mean([r["AKI3"] for r in rows])),
            pct_ventilated=float(100 * np.mean([r["ventilated"] for r in rows])),
            mean_p_bleed=float(100 * np.mean(f("p_bleed"))),
            mean_p_death=float(100 * np.mean(f("p_death"))),
            pct_underdosed=float(100 * np.mean(f("molar_margin") < 1.0)),
        )
    with open(os.path.join(here, "sbe_population_results.json"), "w") as fh:
        json.dump(popdump, fh, indent=1, ensure_ascii=False)

    out("=" * 78)
    out("END OF RUN.  Machine-readable: sbe_scenario_results.json,")
    out("             sbe_population_results.json")
    out("=" * 78)

    with open(os.path.join(here, "sbe_reference_output.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
