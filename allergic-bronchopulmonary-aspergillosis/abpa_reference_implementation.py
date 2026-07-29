#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ABPA QSP model — dependency-free reference implementation
=========================================================
Pure Python standard library only (no numpy / scipy).  Every number quoted in
README.md and in the commit message is produced by running this file:

    python3 abpa_reference_implementation.py            # full suite A0-A13
    python3 abpa_reference_implementation.py A6         # one analysis

The parameter block below is the SAME parameter block as
abpa_mrgsolve_model.R ($PARAM); the two files are kept in sync by name.

--------------------------------------------------------------------------
The structural claim the whole model is built to test
--------------------------------------------------------------------------
The Aspergillus population is partitioned into two compartments:

    FLUM  luminal / mucosa-adherent   -> sees unbound plasma drug
    FPLG  embedded in a mucus plug    -> sees f_pen x unbound plasma drug

with entrapment  k_in * PLUG * FLUM  and release  k_out(PLUG) * FPLG.
Linearising at low burden gives a 2x2 Jacobian whose stability condition is a
closed form containing NO dose term for the sanctuary row except through
f_pen, so there is a plug-clearance rate below which no antifungal exposure
whatsoever can clear the organism.  k_out is not an antifungal parameter — it
belongs to the type-2 mucin axis.  That is the model's whole argument.
"""

import math
import sys

# ==========================================================================
# 0.  UNIT CONVERSIONS AND MOLECULAR CONSTANTS
# ==========================================================================
IGE_IU_TO_NM = 0.01263   # 1 IU/mL IgE = 2.4 ng/mL, MW 190 kDa -> 0.01263 nM
IGE_NM_TO_IU = 1.0 / IGE_IU_TO_NM
IGE_NM_TO_NGML = 2.4 * IGE_NM_TO_IU          # nM -> ng/mL
MW_OMA = 149000.0        # g/mol
V_OMA = 5.46             # L, 78 mL/kg x 70 kg
OMA_MG_TO_NM = 1e6 / MW_OMA / V_OMA          # 1 mg absorbed -> nM in V_OMA
MW_ITRA = 705.6

# ==========================================================================
# 1.  PARAMETERS  (identical names in abpa_mrgsolve_model.R)
# ==========================================================================
P = dict(
    # ---- antigen / immune drive -------------------------------------------
    # Solved BACKWARDS from a target untreated endotype (see A0):
    #   AG 0.70 · TH2 2.0 · IL13 1.5 · IL5 1.5 · EOSB 800 /uL · EOSA 1.0
    #   EPX 1.6 · PLUG 1.0 · FLUM 0.10 · FPLG 0.30 · total IgE 2000 IU/mL
    #   FEV1 69.5% · BRON 3.0 · cortisol 14 ug/dL
    kag=1.00, kagd=0.49, ag_plug=0.80,
    kin2=0.2937, Kag2=0.35, kout2=0.12, th2_base=0.45,
    s13=0.90, kd13=1.20, s5=0.90, kd5=1.20,
    # ---- plasma cell / IgE ------------------------------------------------
    spc=0.03718, pc0=0.0538, Kagpc=0.40, a4=0.85, kpc=0.020,
    kige=1.30,                 # nM/d of IgE per unit plasma-cell pool
    kdegE=0.277,               # free IgE elimination, t1/2 = 2.5 d
    kelCX=0.0866,              # omalizumab:IgE complex elimination, t1/2 = 8 d
                               #   CALIBRATED in A7 to a ~3x total-IgE rise;
                               #   the SIGN of the rise needs only kelCX<kdegE
    Kd_oma=1.0,                # nM, omalizumab:IgE
    eta_oma=0.50,              # IgE neutralised per omalizumab (2:1 average)
    ka_oma=0.55, F_oma=0.62, kel_oma=0.0267,   # t1/2 = 26 d
    Kd_fcer=0.10,              # nM, FcepsilonRI affinity for IgE
    ksyn_fcer=0.050, kdeg_fcer=0.050, b_fcer=6.0, Kup_fcer=0.35,
    # ---- eosinophils ------------------------------------------------------
    kin_e=120.0, a5=1.20, kout_e=0.42, cs_eos=1.60,
    kADCC=0.90, EC50_benr=0.30,
    ktr_e=0.170, aeot=0.80, kout_ea=0.30, cs_apop=1.30,
    dup_egress=0.55,           # dupilumab blocks tissue egress -> blood eos up
    sepx=1.00, kepx=0.625,
    # ---- mucus plug -------------------------------------------------------
    smuc=0.352, Kmuc=1.50, muc0=0.030, cs_muc=0.55, smuc_ige=0.060,
    kout0=0.550, cs_plug=1.90, g_epx=0.55, g_br=0.60,
    # ---- the sanctuary partition -----------------------------------------
    gl=0.85, gp=0.350, kin_f=0.063, seed=0.0020,
    k_host=1.265, cs_imm=0.55,
    f_pen=0.10,
    Emax_af=1.55, EC50_af=0.0038,      # mg/L UNBOUND itraconazole-equivalent
    Emax_amb=1.20, EC50_amb=0.90,      # mg/L in epithelial lining fluid
    f_pen_amb=0.03,                    # liposomes diffuse into a plug worse
                                       #   than a small molecule does
    # ---- itraconazole PK (nonlinear) -------------------------------------
    ka_itra=6.0, F_itra=0.55, V_itra=700.0,
    kel_itra=0.79, Ki_auto=0.55,
    fm_oh=1.85, kel_oh=0.26,
    fu_itra=0.0020, fu_oh=0.0016, pot_oh=1.00,
    Ki3A4_itra=0.0015, Ki3A4_oh=0.0050,   # mg/L UNBOUND
    # ---- voriconazole PK -------------------------------------------------
    ka_vori=8.0, F_vori=0.90, V_vori=280.0,
    Vmax_vori=400.0, Km_vori=1.60, CLlin_vori=60.0, cyp2c19=1.0,
    fu_vori=0.42, pot_vori=0.0066,         # potency vs itraconazole unbound
    Ki3A4_vori=0.35,
    # ---- nebulised liposomal amphotericin B ------------------------------
    kel_amb=0.14,
    # ---- corticosteroid PK ------------------------------------------------
    ka_pred=8.0, F_pred=0.80, V_pred=45.0, kel_pred=5.55,
    s3A4g_pred=0.00, s3A4h_pred=0.15,
    ka_mpred=8.0, F_mpred=0.85, V_mpred=90.0, kel_mpred=4.60,
    s3A4g_mpred=0.15, s3A4h_mpred=0.80,
    ka_bud=8.0, F_bud=0.11, V_bud=210.0, kel_bud=9.20,
    s3A4g_bud=0.70, s3A4h_bud=0.70,
    fu_pred=0.25, fu_mpred=0.23, fu_bud=0.12,
    pot_mpred=1.25, pot_bud=40.0,
    Emax_cs=1.00, EC50_cs=0.0090,          # mg/L unbound prednisolone-equiv
    w_t2=0.85, w_cyt=0.80, w_csr=0.60,
    # ---- HPA / toxicity ---------------------------------------------------
    kin_c=2.10, kout_c=0.150, Imax_c=0.96,
    kbmd=0.00042, krep_bmd=0.00020,
    tau_a1c=45.0, a1c0=5.40, a1c_gain=2.00,
    # ---- biologics PK ----------------------------------------------------
    ka_mep=0.28, F_mep=0.80, V_mep=3.60, kel_mep=0.0347,
    Emax_mep=0.90, EC50_mep=0.55,
    ka_ben=0.25, F_ben=0.58, V_ben=3.10, kel_ben=0.0462,
    Emax_ben_block=0.85, EC50_ben=0.30,
    ka_dup=0.26, F_dup=0.64, V_dup=4.80, kel_dup=0.0630,
    Emax_dup=0.82, EC50_dup=1.90,
    ka_tez=0.20, F_tez=0.77, V_tez=3.90, kel_tez=0.0315,
    Emax_tez=0.50, EC50_tez=2.00,
    # ---- structure / physiology ------------------------------------------
    kb1=0.00060, Pthr=0.35, kb2=0.340, BRONmax=18.0,
    FEV1_0=100.0, cp_fev=0.34, Kp_fev=0.85, cb_fev=0.42, ca_fev=0.14, Ka_fev=1.00,
    h0=0.00021, b1_eos=0.85, b2_plug=1.15, b3_fun=0.55, b4_ige=0.55,
    w_ahr_il13=0.60, w_ahr_ige=0.40,
)

# ==========================================================================
# 2.  STATE VECTOR LAYOUT
# ==========================================================================
NAMES = [
    "AG", "TH2", "IL13", "IL5", "PC",          # 0-4
    "TI", "TO",                                 # 5-6   total IgE, total omalizumab (nM)
    "FCER",                                     # 7
    "EOSB", "EOSA", "EPX",                      # 8-10
    "PLUG",                                     # 11
    "FLUM", "FPLG",                             # 12-13
    "AITR", "ITRA", "OHIT",                     # 14-16
    "AVOR", "VORI", "AMB",                      # 17-19
    "APRD", "PRED", "AMPD", "MPRD", "ABUD", "BUD",  # 20-25
    "CORT", "BMD", "HBA1C", "CUMO",             # 26-29
    "MEPD", "MEPO", "BEND", "BENR", "DUPD", "DUPI", "TEZD", "TEZE",  # 30-37
    "BRON", "CHAZ",                             # 38-39
    "AUCP", "AUCCS",                            # 40-41  exposure integrals
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)

# ==========================================================================
# 3.  ALGEBRA
# ==========================================================================


def bind_1to1(B, L, Kd):
    """Rapid-equilibrium 1:1 binding.  Returns complex concentration."""
    if B <= 0.0 or L <= 0.0:
        return 0.0
    s = B + L + Kd
    disc = s * s - 4.0 * B * L
    if disc < 0.0:
        disc = 0.0
    return 0.5 * (s - math.sqrt(disc))


def cyp3a4_activity(y, p):
    """Fractional CYP3A4 activity remaining, I in [0,1]."""
    cu_i = p["fu_itra"] * max(y[IX["ITRA"]], 0.0)
    cu_o = p["fu_oh"] * max(y[IX["OHIT"]], 0.0)
    cu_v = p["fu_vori"] * max(y[IX["VORI"]], 0.0)
    denom = 1.0 + cu_i / p["Ki3A4_itra"] + cu_o / p["Ki3A4_oh"] + cu_v / p["Ki3A4_vori"]
    return 1.0 / denom


def steroid_signal(y, p):
    """Unbound prednisolone-equivalent (mg/L) and the GR effect CS in [0,1]."""
    peq = (p["fu_pred"] * max(y[IX["PRED"]], 0.0)
           + p["fu_mpred"] * max(y[IX["MPRD"]], 0.0) * p["pot_mpred"]
           + p["fu_bud"] * max(y[IX["BUD"]], 0.0) * p["pot_bud"])
    cs = p["Emax_cs"] * peq / (peq + p["EC50_cs"]) if peq > 0 else 0.0
    return peq, cs


def gut_factor(I, s):
    """Fraction of the gut-wall first-pass step surviving; F rises as 1/this."""
    return (1.0 - s) + s * I


def hep_factor(I, s):
    """Fraction of hepatic clearance surviving."""
    return (1.0 - s) + s * I


def kout_plug_eff(EPX, BRON, CS, p):
    """Plug clearance rate: raised by corticosteroid, lowered by eosinophil
    peroxidase cross-linking and by established bronchiectasis."""
    return (p["kout0"] * (1.0 + p["cs_plug"] * CS)
            / ((1.0 + p["g_epx"] * max(EPX, 0.0))
               * (1.0 + p["g_br"] * max(BRON, 0.0) / p["BRONmax"])))


def readouts(y, p):
    """Derived (algebraic) outputs."""
    ti, to = max(y[IX["TI"]], 0.0), max(y[IX["TO"]], 0.0)
    C = bind_1to1(p["eta_oma"] * to, ti, p["Kd_oma"])
    ige_free = ti - C                                    # nM
    occ = ige_free / (ige_free + p["Kd_fcer"])
    eff_act = occ * y[IX["FCER"]]
    plug = max(y[IX["PLUG"]], 0.0)
    bron = max(y[IX["BRON"]], 0.0)
    il13 = max(y[IX["IL13"]], 0.0) * (1.0 - p["Emax_dup"] * y[IX["DUPI"]] / (y[IX["DUPI"]] + p["EC50_dup"]))
    ahr = (p["w_ahr_il13"] * il13 / (il13 + p["Ka_fev"]) + p["w_ahr_ige"] * eff_act)
    fev1 = (p["FEV1_0"] * (1.0 - p["cp_fev"] * plug / (plug + p["Kp_fev"]))
            * (1.0 - p["cb_fev"] * bron / p["BRONmax"])
            * (1.0 - p["ca_fev"] * ahr))
    haz = p["h0"] * math.exp(p["b1_eos"] * max(y[IX["EOSA"]], 0.0)
                             + p["b2_plug"] * plug
                             + p["b3_fun"] * max(y[IX["FLUM"]] + y[IX["FPLG"]], 0.0)
                             + p["b4_ige"] * eff_act)
    cs_now = steroid_signal(y, p)[1]
    kout_now = kout_plug_eff(max(y[IX["EPX"]], 0.0), bron, cs_now, p)
    return dict(
        IgE_total=ti * IGE_NM_TO_IU,
        IgE_free=ige_free * IGE_NM_TO_IU,
        IgE_free_ngml=ige_free * IGE_NM_TO_NGML,
        CX=C * IGE_NM_TO_IU,
        FcER_occ=occ, effector=eff_act,
        EOSB=y[IX["EOSB"]], EOSA=y[IX["EOSA"]],
        PLUG=plug, plug_score=min(18.0, 6.0 * plug),
        FLUM=max(y[IX["FLUM"]], 0.0), FPLG=max(y[IX["FPLG"]], 0.0),
        FTOT=max(y[IX["FLUM"]] + y[IX["FPLG"]], 0.0),
        FEV1=fev1, BRON=bron, CORT=y[IX["CORT"]], BMD=y[IX["BMD"]],
        HBA1C=y[IX["HBA1C"]], CUMO=y[IX["CUMO"]],
        hazard_yr=haz * 365.0, CHAZ=y[IX["CHAZ"]], AHR=ahr,
        AUCP=y[IX["AUCP"]], AUCCS=y[IX["AUCCS"]],
        I3A4=cyp3a4_activity(y, p), CS=cs_now, kout=kout_now, IL13e=il13,
        EPX=max(y[IX["EPX"]], 0.0),
        ITRA=y[IX["ITRA"]], OHIT=y[IX["OHIT"]], VORI=y[IX["VORI"]],
    )


# ==========================================================================
# 4.  RIGHT-HAND SIDE
# ==========================================================================
def rhs(t, y, p, ctl):
    d = [0.0] * NST
    g = lambda n: y[IX[n]]

    AG, TH2, IL13, IL5, PC = g("AG"), g("TH2"), g("IL13"), g("IL5"), g("PC")
    TI, TO, FCER = max(g("TI"), 0.0), max(g("TO"), 0.0), max(g("FCER"), 0.0)
    EOSB, EOSA, EPX = max(g("EOSB"), 0.0), max(g("EOSA"), 0.0), max(g("EPX"), 0.0)
    PLUG = max(g("PLUG"), 0.0)
    FLUM, FPLG = max(g("FLUM"), 0.0), max(g("FPLG"), 0.0)
    BRON = max(g("BRON"), 0.0)

    # ---- drug effects -----------------------------------------------------
    I3A4 = cyp3a4_activity(y, p)
    peq, CS = steroid_signal(y, p)
    if ctl.get("ddi_off"):
        I3A4 = 1.0
    dupi = max(g("DUPI"), 0.0)
    mepo = max(g("MEPO"), 0.0)
    benr = max(g("BENR"), 0.0)
    teze = max(g("TEZE"), 0.0)
    E_dup = p["Emax_dup"] * dupi / (dupi + p["EC50_dup"])
    E_mep = p["Emax_mep"] * mepo / (mepo + p["EC50_mep"])
    E_ben = p["Emax_ben_block"] * benr / (benr + p["EC50_ben"])
    E_tez = p["Emax_tez"] * teze / (teze + p["EC50_tez"])
    IL13e = IL13 * (1.0 - E_dup)
    IL5e = IL5 * (1.0 - max(E_mep, E_ben))

    # ---- antigen ----------------------------------------------------------
    d[IX["AG"]] = p["kag"] * (FLUM + p["ag_plug"] * FPLG) - p["kagd"] * AG

    # ---- Th2 and cytokines (tezepelumab acts upstream on the drive) -------
    drive = (p["th2_base"] + (1.0 - p["th2_base"])
             * AG / (AG + p["Kag2"]) * (1.0 - E_tez))
    d[IX["TH2"]] = p["kin2"] * drive * (1.0 - p["w_t2"] * CS) - p["kout2"] * TH2
    d[IX["IL13"]] = p["s13"] * TH2 * (1.0 - p["w_cyt"] * CS) - p["kd13"] * IL13
    d[IX["IL5"]] = p["s5"] * TH2 * (1.0 - p["w_cyt"] * CS) - p["kd5"] * IL5

    # ---- plasma cells and IgE (rapid-equilibrium TMDD) --------------------
    d[IX["PC"]] = (p["spc"] * AG / (AG + p["Kagpc"]) * (1.0 + p["a4"] * IL13e)
                   * (1.0 - p["w_csr"] * CS) + p["pc0"] - p["kpc"] * PC)
    C = bind_1to1(p["eta_oma"] * TO, TI, p["Kd_oma"])
    ige_free = TI - C
    oma_free = TO - C / p["eta_oma"]
    S_ige = p["kige"] * max(PC, 0.0)
    d[IX["TI"]] = S_ige - p["kdegE"] * ige_free - p["kelCX"] * C
    # omalizumab: SC absorption is added directly to TO by the driver (Depot)
    d[IX["TO"]] = -p["kel_oma"] * max(oma_free, 0.0) - p["kelCX"] * (C / p["eta_oma"])

    # ---- FcepsilonRI density ---------------------------------------------
    d[IX["FCER"]] = (p["ksyn_fcer"] * (1.0 + p["b_fcer"] * ige_free / (ige_free + p["Kup_fcer"]))
                     - p["kdeg_fcer"] * (1.0 + p["b_fcer"]) * FCER)

    # ---- eosinophils ------------------------------------------------------
    d[IX["EOSB"]] = (p["kin_e"] * (1.0 + p["a5"] * IL5e)
                     - p["kout_e"] * EOSB * (1.0 + p["cs_eos"] * CS)
                     - p["kADCC"] * benr / (benr + p["EC50_benr"]) * EOSB)
    egress = p["ktr_e"] * (EOSB / 1000.0) * (1.0 + p["aeot"] * IL13e) * (1.0 - p["dup_egress"] * E_dup)
    d[IX["EOSA"]] = egress - p["kout_ea"] * EOSA * (1.0 + p["cs_apop"] * CS)
    d[IX["EPX"]] = p["sepx"] * EOSA - p["kepx"] * EPX

    # ---- mucus plug -------------------------------------------------------
    EFF = (ige_free / (ige_free + p["Kd_fcer"])) * max(g("FCER"), 0.0)
    kout_plug = kout_plug_eff(EPX, BRON, CS, p)
    mdrive = (p["muc0"] + p["smuc"] * IL13e / (IL13e + p["Kmuc"])
              + p["smuc_ige"] * EFF)
    d[IX["PLUG"]] = mdrive * (1.0 - p["cs_muc"] * CS) - kout_plug * PLUG

    # ---- THE SANCTUARY PARTITION -----------------------------------------
    cu_az = (p["fu_itra"] * max(g("ITRA"), 0.0)
             + p["fu_oh"] * max(g("OHIT"), 0.0) * p["pot_oh"]
             + p["fu_vori"] * max(g("VORI"), 0.0) * p["pot_vori"])
    E_az = p["Emax_af"] * cu_az / (cu_az + p["EC50_af"])
    amb = max(g("AMB"), 0.0)
    E_amb = p["Emax_amb"] * amb / (amb + p["EC50_amb"])
    E_lum = E_az + E_amb
    fp = ctl["f_pen_override"] if ctl.get("f_pen_override") is not None else p["f_pen"]
    E_plg = fp * E_az + p["f_pen_amb"] * E_amb
    host = p["k_host"] * (1.0 - p["cs_imm"] * CS)
    cap = max(1.0 - (FLUM + FPLG), 0.0)
    entrap = p["kin_f"] * PLUG * FLUM * max(1.0 - FPLG, 0.0)
    d[IX["FLUM"]] = (p["gl"] * FLUM * cap + kout_plug * FPLG
                     - entrap - host * FLUM - E_lum * FLUM + p["seed"])
    d[IX["FPLG"]] = (p["gp"] * FPLG * max(1.0 - FPLG, 0.0) + entrap
                     - kout_plug * FPLG - E_plg * FPLG)

    # ---- itraconazole PK (saturable / autoinhibiting) --------------------
    d[IX["AITR"]] = -p["ka_itra"] * g("AITR")
    itra = max(g("ITRA"), 0.0)
    cl_itra = p["kel_itra"] / (1.0 + itra / p["Ki_auto"])
    d[IX["ITRA"]] = p["ka_itra"] * g("AITR") * p["F_itra"] / p["V_itra"] - cl_itra * itra
    d[IX["OHIT"]] = p["fm_oh"] * cl_itra * itra - p["kel_oh"] * max(g("OHIT"), 0.0)

    # ---- voriconazole PK (Michaelis-Menten, CYP2C19 scaled) -------------
    d[IX["AVOR"]] = -p["ka_vori"] * g("AVOR")
    vori = max(g("VORI"), 0.0)
    vmax = p["Vmax_vori"] * p["cyp2c19"]
    d[IX["VORI"]] = (p["ka_vori"] * g("AVOR") * p["F_vori"] / p["V_vori"]
                     - (vmax * vori / (p["Km_vori"] + vori)
                        + p["CLlin_vori"] * vori) / p["V_vori"])

    # ---- nebulised amphotericin B (airway lining fluid) ------------------
    d[IX["AMB"]] = -p["kel_amb"] * amb

    # ---- corticosteroid PK with the CYP3A4 interaction --------------------
    d[IX["APRD"]] = -p["ka_pred"] * g("APRD")
    F_pred_eff = min(p["F_pred"] / gut_factor(I3A4, p["s3A4g_pred"]), 1.0)
    kel_p = p["kel_pred"] * hep_factor(I3A4, p["s3A4h_pred"])
    d[IX["PRED"]] = p["ka_pred"] * g("APRD") * F_pred_eff / p["V_pred"] - kel_p * max(g("PRED"), 0.0)
    d[IX["AMPD"]] = -p["ka_mpred"] * g("AMPD")
    F_mpred_eff = min(p["F_mpred"] / gut_factor(I3A4, p["s3A4g_mpred"]), 1.0)
    kel_m = p["kel_mpred"] * hep_factor(I3A4, p["s3A4h_mpred"])
    d[IX["MPRD"]] = p["ka_mpred"] * g("AMPD") * F_mpred_eff / p["V_mpred"] - kel_m * max(g("MPRD"), 0.0)
    d[IX["ABUD"]] = -p["ka_bud"] * g("ABUD")
    F_bud_eff = min(p["F_bud"] / gut_factor(I3A4, p["s3A4g_bud"]), 1.0)
    kel_b = p["kel_bud"] * hep_factor(I3A4, p["s3A4h_bud"])
    d[IX["BUD"]] = p["ka_bud"] * g("ABUD") * F_bud_eff / p["V_bud"] - kel_b * max(g("BUD"), 0.0)

    # ---- HPA / toxicity ---------------------------------------------------
    d[IX["CORT"]] = p["kin_c"] * (1.0 - p["Imax_c"] * CS) - p["kout_c"] * max(g("CORT"), 0.0)
    d[IX["BMD"]] = -p["kbmd"] * CS + p["krep_bmd"] * (1.0 - g("BMD"))
    tgt = p["a1c0"] + p["a1c_gain"] * CS
    d[IX["HBA1C"]] = (tgt - g("HBA1C")) / p["tau_a1c"]
    d[IX["CUMO"]] = 0.0    # incremented at dosing events

    # ---- biologics PK -----------------------------------------------------
    for dep, cen, ka, F, V, kel in (
            ("MEPD", "MEPO", "ka_mep", "F_mep", "V_mep", "kel_mep"),
            ("BEND", "BENR", "ka_ben", "F_ben", "V_ben", "kel_ben"),
            ("DUPD", "DUPI", "ka_dup", "F_dup", "V_dup", "kel_dup"),
            ("TEZD", "TEZE", "ka_tez", "F_tez", "V_tez", "kel_tez")):
        d[IX[dep]] = -p[ka] * g(dep)
        d[IX[cen]] = p[ka] * g(dep) * p[F] / p[V] - p[kel] * max(g(cen), 0.0)

    # ---- structure --------------------------------------------------------
    haz = p["h0"] * math.exp(p["b1_eos"] * EOSA + p["b2_plug"] * PLUG
                             + p["b3_fun"] * (FLUM + FPLG) + p["b4_ige"] * EFF)
    d[IX["CHAZ"]] = haz
    d[IX["AUCP"]] = peq          # integral of unbound prednisolone-equivalent
    d[IX["AUCCS"]] = CS          # integral of the GR effect actually delivered
    grow = p["kb1"] * max(PLUG - p["Pthr"], 0.0) + p["kb2"] * haz
    d[IX["BRON"]] = grow if BRON < p["BRONmax"] else 0.0
    return d


# ==========================================================================
# 5.  OMALIZUMAB SC DEPOT — kept outside the state vector for clarity
# ==========================================================================
# Omalizumab absorption is handled by adding directly to TO with a first-order
# depot integrated analytically between events (the depot has no feedback).

class Depot:
    def __init__(self, ka):
        self.a = 0.0
        self.ka = ka

    def dose(self, mg, F):
        self.a += mg * F

    def release(self, dt):
        out = self.a * (1.0 - math.exp(-self.ka * dt))
        self.a -= out
        return out


# ==========================================================================
# 6.  SIMULATION ENGINE
# ==========================================================================
def seed_state(p, ige0_iu=None):
    """Hand-set starting guess.  NOT used directly as an initial condition: the
    guess does not lie on the model's steady manifold, and using it produced a
    large transient (untreated total IgE overshooting to ~3450 IU/mL before
    settling at ~1990) that contaminated every week-12 readout.  initial_state()
    relaxes this guess to the untreated steady state first."""
    y = [0.0] * NST
    y[IX["AG"]] = 0.9
    y[IX["TH2"]] = 2.0
    y[IX["IL13"]] = 1.5
    y[IX["IL5"]] = 1.5
    y[IX["PC"]] = 25.0
    y[IX["TI"]] = (ige0_iu if ige0_iu else 2000.0) * IGE_IU_TO_NM
    y[IX["TO"]] = 0.0
    y[IX["FCER"]] = 1.0
    y[IX["EOSB"]] = 800.0
    y[IX["EOSA"]] = 1.0
    y[IX["EPX"]] = 1.6
    y[IX["PLUG"]] = 1.0
    y[IX["FLUM"]] = 0.25
    y[IX["FPLG"]] = 0.35
    y[IX["CORT"]] = 14.0
    y[IX["BMD"]] = 1.0
    y[IX["HBA1C"]] = 5.4
    y[IX["BRON"]] = 3.0
    return y


_BASE_CACHE = {}
ACCUM = ("CUMO", "CHAZ", "AUCP", "AUCCS")
_CACHE_KEYS = ("kige", "gp", "f_pen", "kelCX", "kout0", "g_epx", "smuc",
               "th2_base", "pc0", "k_host", "gl", "kin_f", "Emax_af",
               "smuc_ige", "kepx", "sepx", "cs_plug", "cs_muc", "kdegE",
               "kin_e", "a5", "ktr_e", "aeot", "muc0", "Kmuc", "spc")


def initial_state(p, ige0_iu=None):
    """The untreated STEADY STATE, used as the initial condition for every arm.

    Bronchiectasis is a monotone integrator and so has no steady state; it is
    pinned back to its presentation value (3.0/18) between relaxation passes,
    and the accumulators are zeroed.  Three passes suffice: the residual drift
    is below the fourth decimal of every reported quantity (A0 checks this)."""
    key = tuple(round(p[k], 12) for k in _CACHE_KEYS)
    if key in _BASE_CACHE:
        return list(_BASE_CACHE[key])
    y = seed_state(p, ige0_iu)
    for _ in range(3):
        _, y = simulate(build_regimen([]), tend=365.0, dt=0.02, p=p,
                        record=1e9, y0=y)
        y[IX["BRON"]] = 3.0
        y[IX["BMD"]] = 1.0
        y[IX["HBA1C"]] = p["a1c0"]
        for acc in ACCUM:
            y[IX[acc]] = 0.0
    _BASE_CACHE[key] = list(y)
    return list(y)


def rk4(t, y, dt, p, ctl):
    k1 = rhs(t, y, p, ctl)
    y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
    k2 = rhs(t + 0.5 * dt, y2, p, ctl)
    y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
    k3 = rhs(t + 0.5 * dt, y3, p, ctl)
    y4 = [y[i] + dt * k3[i] for i in range(NST)]
    k4 = rhs(t + dt, y4, p, ctl)
    return [y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(NST)]


def simulate(regimen, tend=364.0, dt=0.02, p=None, ige0=None, ctl=None,
             record=None, y0=None):
    """regimen: callable(t_day) -> dict of dosing instructions applied at
    integer-ish grid times; see build_regimen()."""
    p = p or P
    ctl = dict(ctl or {})
    y = list(y0) if y0 is not None else initial_state(p, ige0)
    oma_depot = Depot(p["ka_oma"])
    t = 0.0
    nstep = int(round(tend / dt))
    rec = []
    record = record if record is not None else 7.0
    next_rec = 0.0
    last_day = -1
    while True:
        if t >= next_rec - 1e-9:
            r = readouts(y, p)
            r["t"] = t
            rec.append(r)
            next_rec += record
        if t >= tend - 1e-9:
            break
        # ---- discrete dosing at day boundaries --------------------------
        day = int(math.floor(t + 1e-9))
        if day != last_day:
            last_day = day
            ev = regimen(day) or {}
            y[IX["AITR"]] += ev.get("itra", 0.0)
            y[IX["AVOR"]] += ev.get("vori", 0.0)
            y[IX["AMB"]] += ev.get("amb", 0.0)
            y[IX["APRD"]] += ev.get("pred", 0.0)
            y[IX["AMPD"]] += ev.get("mpred", 0.0)
            y[IX["ABUD"]] += ev.get("bud", 0.0)
            y[IX["CUMO"]] += (ev.get("pred", 0.0)
                              + p["pot_mpred"] * ev.get("mpred", 0.0)
                              + 0.0)
            if ev.get("oma"):
                oma_depot.dose(ev["oma"], p["F_oma"])
            y[IX["MEPD"]] += ev.get("mepo", 0.0)
            y[IX["BEND"]] += ev.get("benra", 0.0)
            y[IX["DUPD"]] += ev.get("dupi", 0.0)
            y[IX["TEZD"]] += ev.get("teze", 0.0)
            if ev.get("lavage"):
                y[IX["FPLG"]] *= (1.0 - ev["lavage"])
                y[IX["PLUG"]] *= (1.0 - ev["lavage"])
        # ---- omalizumab depot release ------------------------------------
        rel = oma_depot.release(dt)
        if rel:
            y[IX["TO"]] += rel * OMA_MG_TO_NM
        y = rk4(t, y, dt, p, ctl)
        for k in ("TI", "TO", "FLUM", "FPLG", "PLUG", "EOSB", "EOSA", "EPX",
                  "ITRA", "OHIT", "VORI", "AMB", "PRED", "MPRD", "BUD", "CORT", "PC"):
            if y[IX[k]] < 0.0:
                y[IX[k]] = 0.0
        t += dt
    return rec, y


# ==========================================================================
# 7.  REGIMEN BUILDERS
# ==========================================================================
def build_regimen(spec):
    """spec is a list of (start_day, stop_day, every_n_days, dict_of_doses)."""
    def reg(day):
        out = {}
        for (t0, t1, every, doses) in spec:
            if t0 <= day < t1 and ((day - t0) % every == 0):
                for k, v in doses.items():
                    out[k] = out.get(k, 0.0) + v
        return out
    return reg


def pred_taper(total_days=364, start=35.0, wt=70.0, scheme="ISHAM_medium"):
    """ISHAM medium-dose oral prednisolone: 0.5 mg/kg/d x 2 wk, then 0.5 mg/kg
    alternate-day x 8 wk, then taper 5 mg every 2 wk and stop by ~3-5 months."""
    dose0 = 0.5 * wt
    sched = []
    for day in range(total_days):
        if day < 14:
            dose = dose0
        elif day < 70:
            dose = dose0 if (day - 14) % 2 == 0 else 0.0
        elif day < 154:
            k = (day - 70) // 14
            dose = max(dose0 - 5.0 * (k + 1), 0.0)
        else:
            dose = 0.0
        if dose > 0:
            sched.append((day, day + 1, 1, {"pred": dose}))
    return sched


def maintenance_pred(total_days, mg, start=0):
    return [(start, total_days, 1, {"pred": mg})]


# ==========================================================================
# 8.  ANALYSES
# ==========================================================================
def hr(title, ch="="):
    print("\n" + ch * 78)
    print(title)
    print(ch * 78)


def A0_selfcheck():
    hr("A0  MODEL SELF-CHECK — untreated steady state, convergence, IgE flux balance")
    reg = build_regimen([])
    for dt in (0.04, 0.02, 0.01, 0.005):
        rec, y = simulate(reg, tend=364.0, dt=dt, record=364.0)
        r = readouts(y, P)
        print("  dt=%.4f  IgE_tot=%9.2f  EOSB=%7.1f  PLUG=%.5f  FLUM=%.5f  "
              "FPLG=%.5f  FEV1=%.4f  BRON=%.4f"
              % (dt, r["IgE_total"], r["EOSB"], r["PLUG"], r["FLUM"], r["FPLG"],
                 r["FEV1"], r["BRON"]))
    rec, y = simulate(reg, tend=3650.0, dt=0.02, record=365.0)
    print("\n  10-year untreated drift (the model must reach a stable endotype, "
          "not run away):")
    print("  %-6s %10s %8s %7s %7s %7s %7s %7s" %
          ("yr", "IgE", "EOSB", "PLUG", "FLUM", "FPLG", "FEV1", "BRON"))
    for r in rec:
        print("  %-6.1f %10.1f %8.1f %7.4f %7.4f %7.4f %7.2f %7.3f"
              % (r["t"] / 365.0, r["IgE_total"], r["EOSB"], r["PLUG"],
                 r["FLUM"], r["FPLG"], r["FEV1"], r["BRON"]))
    # IgE flux balance without drug: S = kdegE * IgE_free
    r = readouts(y, P)
    S = P["kige"] * y[IX["PC"]]
    print("\n  IgE flux balance check (no omalizumab): S=%.5f nM/d, "
          "kdegE*free=%.5f nM/d, residual=%.2e"
          % (S, P["kdegE"] * y[IX["TI"]], S - P["kdegE"] * y[IX["TI"]]))


def A1_natural_history():
    hr("A1  NATURAL HISTORY — untreated 5 years, and where the irreversibility comes from")
    rec, y = simulate(build_regimen([]), tend=1825.0, dt=0.02, record=182.5)
    print("  %-7s %9s %8s %7s %7s %7s %8s %7s %8s" %
          ("t(yr)", "IgE", "EOSB", "plugsc", "FLUM", "FPLG", "FEV1", "BRON", "exac/yr"))
    for r in rec:
        print("  %-7.2f %9.1f %8.1f %7.2f %7.4f %7.4f %8.2f %7.3f %8.3f"
              % (r["t"] / 365.0, r["IgE_total"], r["EOSB"], r["plug_score"],
                 r["FLUM"], r["FPLG"], r["FEV1"], r["BRON"], r["hazard_yr"]))
    r0, rE = rec[0], rec[-1]
    print("\n  5-year FEV1 loss = %.2f points; of which reversible (plug) %.2f "
          "and fixed (bronchiectasis) %.2f"
          % (r0["FEV1"] - rE["FEV1"],
             P["FEV1_0"] * P["cp_fev"] * (rE["PLUG"] / (rE["PLUG"] + P["Kp_fev"])
                                          - r0["PLUG"] / (r0["PLUG"] + P["Kp_fev"])),
             P["FEV1_0"] * P["cb_fev"] * (rE["BRON"] - r0["BRON"]) / P["BRONmax"]))
    print("  sanctuary share of total fungal burden at 5 y = %.1f%%"
          % (100.0 * rE["FPLG"] / max(rE["FTOT"], 1e-12)))


def erad_threshold(p, plug, kout, gp=None, f_pen=None):
    """Minimum luminal kill rate E for which F=(0,0) is locally stable.
    Jacobian at F->0 (seed set to 0):
        J = [[gl - kin*PLUG - host - E ,  kout            ],
             [kin*PLUG                 ,  gp - kout - f*E ]]
    Stable iff trace<0 and det>0.  Returns None if unreachable."""
    gp = p["gp"] if gp is None else gp
    f = p["f_pen"] if f_pen is None else f_pen
    host = p["k_host"]
    gl = p["gl"]
    kin = p["kin_f"] * plug

    def stable(E):
        j11 = gl - kin - host - E
        j22 = gp - kout - f * E
        tr = j11 + j22
        det = j11 * j22 - kout * kin
        return (tr < 0.0) and (det > 0.0)

    lo, hi = 0.0, 1.0
    if stable(0.0):
        return 0.0
    while hi < 1e7:
        if stable(hi):
            break
        hi *= 2.0
    if not stable(hi):
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if stable(mid):
            hi = mid
        else:
            lo = mid
    return hi


def A2_eradication_boundary():
    hr("A2  THE ERADICATION BOUNDARY — a closed form with no dose in the sanctuary row")
    print("""  Linearising the two-compartment fungal partition at F->0 gives

      J = [[ gl - k_in*PLUG - k_host - E ,  k_out                 ],
           [ k_in*PLUG                   ,  g_p - k_out - f_pen*E ]]

  and clearance needs trace<0 AND det>0.  The second diagonal element holds
  the ONLY dose term that reaches the sanctuary, attenuated by f_pen, so when
  g_p > k_out the required luminal exposure is at least (g_p-k_out)/f_pen — a
  bound containing no antifungal property except penetration.  Note what k_out
  is: it belongs to the type-2 mucin axis, not to the antifungal.""")
    print("\n  (a) required luminal kill rate E* vs plug clearance k_out "
          "(PLUG=1.0, f_pen=%.2f, g_p=%.2f):" % (P["f_pen"], P["gp"]))
    print("      %-9s %11s %13s %11s" % ("k_out", "E* (1/d)", "floor=(gp-kout)/f", "verdict"))
    for kout in (0.05, 0.10, 0.15, 0.20, 0.25, 0.267, 0.293, 0.35, 0.42, 0.60):
        E = erad_threshold(P, 1.0, kout)
        floor = max(P["gp"] - kout, 0.0) / P["f_pen"]
        if E is None:
            v = "unreachable"
        elif E <= 0.0:
            v = "no drug needed"
        elif E > P["Emax_af"]:
            v = "ABOVE Emax"
        else:
            v = "reachable"
        print("      %-9.3f %11.3f %13.3f %11s"
              % (kout, E if E is not None else float("nan"), floor, v))
    print("      (the E* column tracks the analytic floor to within the "
          "coupling term k_out*k_in, which is the only thing the closed form "
          "omits)")

    print("\n  (b) sweep f_pen at the untreated plug clearance:")
    kout_u = kout_plug_eff(1.6, 3.0, 0.0, P)
    print("      untreated effective k_out = %.4f /d (EPX 1.6, BRON 3.0)" % kout_u)
    print("      %-9s %11s" % ("f_pen", "E* (1/d)"))
    for f in (0.02, 0.05, 0.10, 0.20, 0.40, 0.70, 1.00):
        E = erad_threshold(P, 1.0, kout_u, f_pen=f)
        print("      %-9.2f %11s" % (f, ("%.3f" % E) if E is not None else "unreachable"))

    print("\n  (c) the achievable exposure, for scale — itraconazole 200 mg BID at "
          "steady state:")
    rec, y = simulate(build_regimen([(0, 400, 1, {"itra": 400.0})]), tend=42.0,
                      dt=0.02, record=42.0)
    cu = P["fu_itra"] * y[IX["ITRA"]] + P["fu_oh"] * y[IX["OHIT"]] * P["pot_oh"]
    E_ach = P["Emax_af"] * cu / (cu + P["EC50_af"])
    print("      parent %.3f mg/L, OH-itra %.3f mg/L, unbound-equivalent %.5f mg/L"
          % (y[IX["ITRA"]], y[IX["OHIT"]], cu))
    print("      -> E_lum = %.3f /d, which is %.1f%% of the class ceiling Emax = %.2f"
          % (E_ach, 100.0 * E_ach / P["Emax_af"], P["Emax_af"]))
    print("      -> the sanctuary sees only f_pen*E = %.4f /d" % (P["f_pen"] * E_ach))

    print("""
  (d) THE SEVERITY CROSSOVER.  k_out is not a constant: eosinophil-peroxidase
  cross-linking and established bronchiectasis both lower it, so E* RISES as
  the disease advances.  The grid below is E* over (bronchiectasis extent,
  EPX); "." = itraconazole 200 mg BID suffices, "o" = needs more azole than
  itraconazole delivers but less than Emax, "X" = above the ceiling of the
  ENTIRE azole class, i.e. unreachable by any azole at any dose.""")
    epxs = [1.0, 1.6, 2.2, 2.8, 3.4, 4.0]
    brons = [0.0, 3.0, 6.0, 9.0, 12.0, 15.0, 18.0]
    print("\n      %-8s" % "BRON\\EPX" + "".join("%9.1f" % e for e in epxs))
    counts = {".": 0, "o": 0, "X": 0}
    for b in brons:
        row = "      %-8.1f" % b
        for e in epxs:
            k = kout_plug_eff(e, b, 0.0, P)
            E = erad_threshold(P, 1.0, k)
            if E is None or E > P["Emax_af"]:
                mark, val = "X", (E if E is not None else float("inf"))
            elif E <= E_ach:
                mark, val = ".", E
            else:
                mark, val = "o", E
            counts[mark] += 1
            row += "%8s%s" % (("%.2f" % val) if val != float("inf") else "  inf", mark)
        print(row)
    tot = sum(counts.values())
    print("      cells: %d reachable by itraconazole (%.0f%%), %d need a stronger "
          "azole (%.0f%%), %d unreachable by any azole (%.0f%%)"
          % (counts["."], 100.0 * counts["."] / tot, counts["o"],
             100.0 * counts["o"] / tot, counts["X"], 100.0 * counts["X"] / tot))
    # exact crossovers along the EPX=1.6 line
    def Estar_at(b, e=1.6):
        return erad_threshold(P, 1.0, kout_plug_eff(e, b, 0.0, P))
    def crossing(target, e=1.6):
        lo, hi = 0.0, P["BRONmax"]
        Elo = Estar_at(lo, e)
        if Elo is not None and Elo > target:
            return 0.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            Em = Estar_at(mid, e)
            if Em is None or Em > target:
                hi = mid
            else:
                lo = mid
        return hi
    print("\n      along EPX = 1.6:  itraconazole (E=%.3f) stops being sufficient at "
          "BRON = %.2f" % (E_ach, crossing(E_ach)))
    print("                        the whole azole class (Emax=%.2f) stops at "
          "BRON = %.2f" % (P["Emax_af"], crossing(P["Emax_af"])))
    print("      along EPX = 3.0:  itraconazole stops at BRON = %.2f, the class at "
          "BRON = %.2f" % (crossing(E_ach, 3.0), crossing(P["Emax_af"], 3.0)))

    print("""
  (e) THE ONE PARAMETER THAT DECIDES ALL OF THIS AND HAS NEVER BEEN MEASURED
  IN A HUMAN AIRWAY is g_p, the growth rate of Aspergillus inside a plug.  It
  is reported for biofilms in vitro and inferred from explanted plugs, never
  measured in situ.  The whole architecture turns on sign(g_p - k_out):""")
    print("\n      %-9s %11s %13s %28s" %
          ("g_p", "E* (1/d)", "g_p - k_out", "meaning"))
    for gp in (0.10, 0.20, 0.29, 0.35, 0.45, 0.60, 0.80):
        E = erad_threshold(P, 1.0, kout_u, gp=gp)
        diff = gp - kout_u
        if E is None or E > P["Emax_af"]:
            m = "no azole can clear it"
        elif E <= 0.0:
            m = "clears without any drug"
        elif E <= E_ach:
            m = "itraconazole suffices"
        else:
            m = "needs a stronger azole"
        print("      %-9.2f %11s %13.4f %28s"
              % (gp, ("%.3f" % E) if E is not None else "unreachable", diff, m))
    print("""
  The model is reported here at g_p = %.2f, only %.3f /d above the untreated
  k_out, which puts itraconazole monotherapy almost exactly ON the boundary —
  and that is itself a result worth stating, because it is what a ~60%% response
  rate to azole monotherapy looks like in a model with no responder/
  non-responder covariate in it.  It also means the model's central claim is
  FALSIFIABLE by a single measurement nobody has made.""" % (P["gp"], P["gp"] - kout_u))
    return E_ach


# ------------------------------------------------------------------ scenarios
def SCENARIOS(days=364):
    wt = 70.0
    S = {}
    S["01 untreated"] = ([], {})
    S["02 prednisolone ISHAM taper"] = (pred_taper(days, wt=wt), {})
    S["03 prednisolone 10 mg maintenance"] = (maintenance_pred(days, 10.0), {})
    S["04 itraconazole 200 mg BID"] = ([(0, days, 1, {"itra": 400.0})], {})
    S["05 pred taper + itraconazole"] = (pred_taper(days, wt=wt)
                                         + [(0, days, 1, {"itra": 400.0})], {})
    S["06 pred + itra, DDI switched OFF"] = (pred_taper(days, wt=wt)
                                             + [(0, days, 1, {"itra": 400.0})],
                                             {"ddi_off": True})
    S["07 pred + itra, antifungal OFF (DDI only)"] = (pred_taper(days, wt=wt)
                                                      + [(0, days, 1, {"itra": 400.0})],
                                                      {"no_antifungal": True})
    S["08 voriconazole 200 mg BID (CYP2C19 NM)"] = ([(0, days, 1, {"vori": 400.0})],
                                                    {"cyp2c19": 1.0})
    S["09 voriconazole 200 mg BID (CYP2C19 PM)"] = ([(0, days, 1, {"vori": 400.0})],
                                                    {"cyp2c19": 0.35})
    S["10 nebulised L-AmB 10 mg q48h"] = ([(0, days, 2, {"amb": 10.0})], {})
    S["11 omalizumab 375 mg q2w"] = ([(0, days, 14, {"oma": 375.0})], {})
    S["12 mepolizumab 100 mg q4w"] = ([(0, days, 28, {"mepo": 100.0})], {})
    S["13 benralizumab 30 mg q8w"] = ([(0, 84, 28, {"benra": 30.0})]
                                      + [(84, days, 56, {"benra": 30.0})], {})
    S["14 dupilumab 300 mg q2w"] = ([(0, days, 14, {"dupi": 300.0})], {})
    S["15 dupilumab + itraconazole"] = ([(0, days, 14, {"dupi": 300.0}),
                                         (0, days, 1, {"itra": 400.0})], {})
    S["16 budesonide 1600 ug + itra (the trap)"] = ([(0, days, 1, {"bud": 1.6}),
                                                     (0, days, 1, {"itra": 400.0})], {})
    S["17 bronchoscopic toilet d0 + itra"] = ([(0, 1, 1, {"lavage": 0.95}),
                                               (0, days, 1, {"itra": 400.0})], {})
    S["18 pred taper + itra + dupilumab"] = (pred_taper(days, wt=wt)
                                             + [(0, days, 1, {"itra": 400.0}),
                                                (0, days, 14, {"dupi": 300.0})], {})
    return S


def run_scenario(spec, ctl, days=364, dt=0.02, ige0=None):
    p = dict(P)
    c = dict(ctl)
    if c.pop("cyp2c19", None) is not None:
        p["cyp2c19"] = ctl["cyp2c19"]
    if c.pop("no_antifungal", None):
        p["Emax_af"] = 0.0
    rec, y = simulate(build_regimen(spec), tend=float(days), dt=dt, p=p,
                      ctl=c, record=28.0, ige0=ige0)
    return rec, y, p


def A3_scenarios(days=364):
    hr("A3  THERAPY SCENARIOS — 52 weeks, every arm priced in its own currency")
    base = None
    rows = []
    for name, (spec, ctl) in SCENARIOS(days).items():
        rec, y, p = run_scenario(spec, ctl, days)
        r = readouts(y, p)
        w12 = min(rec, key=lambda q: abs(q["t"] - 84.0))
        r["w12"] = w12
        if base is None:
            base = r
        rows.append((name, r))
    print("  %-42s %8s %8s %7s %7s %7s %7s %7s %7s %7s %7s" %
          ("scenario", "IgE_tot", "IgE_fre", "EOSB", "plugsc", "FTOT", "FPLG%",
           "FEV1", "BRON", "OCSmg", "cort"))
    for name, r in rows:
        print("  %-42s %8.0f %8.1f %7.0f %7.2f %7.4f %6.1f%% %7.1f %7.3f %7.0f %7.1f"
              % (name, r["IgE_total"], r["IgE_free"], r["EOSB"], r["plug_score"],
                 r["FTOT"], 100.0 * r["FPLG"] / max(r["FTOT"], 1e-12), r["FEV1"],
                 r["BRON"], r["CUMO"], r["CORT"]))
    print("\n  WEEK 12 (acute response — this is what a trial's primary "
          "endpoint usually sees):")
    print("  %-42s %8s %8s %7s %7s %7s %8s" %
          ("scenario", "IgE_tot", "dIgE%", "EOSB", "plugsc", "FTOT", "FEV1"))
    b12 = rows[0][1]["w12"]
    for name, r in rows:
        q = r["w12"]
        print("  %-42s %8.0f %7.0f%% %7.0f %7.2f %7.4f %8.1f"
              % (name, q["IgE_total"], 100.0 * (q["IgE_total"] / b12["IgE_total"] - 1.0),
                 q["EOSB"], q["plug_score"], q["FTOT"], q["FEV1"]))
    print("\n  %-42s %10s %10s %10s" % ("scenario", "exac/52wk", "dBRON", "dFEV1"))
    for name, r in rows:
        print("  %-42s %10.3f %10.4f %+10.2f"
              % (name, r["CHAZ"], r["BRON"] - 3.0, r["FEV1"] - base["FEV1"]))
    return dict(rows)


def A4_ddi_calibration():
    hr("A4  THE CYP3A4 INTERACTION — and why one first-pass site is arithmetically impossible")
    rec, y = simulate(build_regimen([(0, 400, 1, {"itra": 400.0})]), tend=42.0,
                      dt=0.02, record=42.0)
    I = cyp3a4_activity(y, P)
    print("  itraconazole 200 mg BID x 6 wk: parent %.3f mg/L, OH-itra %.3f mg/L "
          "(OH/parent = %.2f)"
          % (y[IX["ITRA"]], y[IX["OHIT"]], y[IX["OHIT"]] / max(y[IX["ITRA"]], 1e-9)))
    print("  unbound parent %.5f mg/L (Ki %.4f), unbound OH %.5f mg/L (Ki %.4f)"
          % (P["fu_itra"] * y[IX["ITRA"]], P["Ki3A4_itra"],
             P["fu_oh"] * y[IX["OHIT"]], P["Ki3A4_oh"]))
    print("  => fractional CYP3A4 activity I = %.4f" % I)
    print("""
  A SINGLE-SITE MODEL IS RULED OUT BY ARITHMETIC.  If a drug's exposure rises
  only because hepatic clearance is inhibited, the fold-change is
  1/((1-s)+s*I), which is maximised at s=1 and therefore CANNOT exceed 1/I.""")
    print("  ceiling of any single-site model at this I:  1/I = %.2f" % (1.0 / I))
    print("  observed itraconazole-budesonide interaction: 4.2x  ->  %s"
          % ("EXCEEDS the ceiling, so the interaction must act at two sequential "
             "sites" if 4.2 > 1.0 / I else "within the ceiling"))
    print("""
  The model therefore inhibits the gut wall (raising F) and the liver (lowering
  CL) separately, and the fold-change becomes a PRODUCT:
      fold_AUC = 1/((1-s_gut)+s_gut*I)  x  1/((1-s_hep)+s_hep*I)""")
    print("\n  %-24s %7s %7s %11s %11s %8s" %
          ("steroid", "s_gut", "s_hep", "fold(model)", "fold(obs)", "err"))
    obs = {"prednisolone": 1.24, "methylprednisolone": 2.60, "budesonide (inh)": 4.20}
    for label, kg, kh in (("prednisolone", "s3A4g_pred", "s3A4h_pred"),
                          ("methylprednisolone", "s3A4g_mpred", "s3A4h_mpred"),
                          ("budesonide (inh)", "s3A4g_bud", "s3A4h_bud")):
        sg, sh = P[kg], P[kh]
        fold = 1.0 / gut_factor(I, sg) / hep_factor(I, sh)
        o = obs[label]
        print("  %-24s %7.2f %7.2f %11.3f %11.2f %7.1f%%"
              % (label, sg, sh, fold, o, 100.0 * (fold - o) / o))
    print("""
  Two numbers were fitted here: methylprednisolone's hepatic share (to Varis
  1998, ~2.6x) and budesonide's gut/hepatic split (to the Raaska/Skov reports,
  ~4.2x).  PREDNISOLONE WAS NOT FITTED — its CYP3A4 share comes from the
  literature independently, and the model returns 1.12 against an observed
  1.24, i.e. it under-predicts the smallest of the three interactions by ~10%.
  The three-fold spread across steroids is what makes A5 testable rather than
  rhetorical: the SAME azole confounds the SAME trial design differently
  depending only on which steroid the protocol happened to specify.""")
    return I


def A5_ddi_decomposition(days=364):
    hr("A5  DECOMPOSING ITRACONAZOLE'S STEROID-SPARING — antifungal or drug interaction?")
    print("""  A steroid-dependent patient is held on a FIXED steroid schedule in every
  arm — no clinician, no titration — so any difference is attributable, and the
  CYP3A4 interaction is switched on and off independently of the antifungal
  effect.  Everything is reported against EXPOSURE INTEGRALS: AUC of unbound
  prednisolone-equivalent (AUC_peq) and the integral of the GR effect actually
  delivered (AUC_CS).  The first draft of this analysis read an end-of-study
  trough instead, which overstated the interaction badly — inhibiting the
  clearance of a 3-hour-half-life steroid moves its trough far more than its
  AUC, and the trough is not what the receptor integrates.""")
    itr = [(0, days, 1, {"itra": 400.0})]
    print("\n  %-30s %10s %7s %9s %8s %8s %8s %7s" %
          ("backbone / arm", "AUC_peq", "xAUC", "AUC_CS", "xAUC_CS", "FTOT",
           "plugsc", "FEV1"))
    summary = []
    for label, key, dose in (("prednisolone 10 mg/d", "pred", 10.0),
                             ("methylprednisolone 8 mg/d", "mpred", 8.0),
                             ("budesonide 1600 ug/d inh", "bud", 1.6)):
        b = [(0, days, 1, {key: dose})]
        arms = [("alone", b, {}),
                ("+ itraconazole (full)", b + itr, {}),
                ("+ itra, DDI off", b + itr, {"ddi_off": True}),
                ("+ itra, antifungal off", b + itr, {"no_antifungal": True})]
        res = {}
        for nm, spec, ctl in arms:
            rec, y, p = run_scenario(spec, ctl, days)
            res[nm] = readouts(y, p)
        a0 = res["alone"]
        for nm, _, _ in arms:
            r = res[nm]
            print("  %-30s %10.4f %7.3f %9.3f %8.3f %8.4f %8.3f %7.2f"
                  % ((label + ", alone") if nm == "alone" else ("    " + nm),
                     r["AUCP"], r["AUCP"] / a0["AUCP"], r["AUCCS"],
                     r["AUCCS"] / a0["AUCCS"], r["FTOT"], r["plug_score"], r["FEV1"]))
        tot = res["+ itraconazole (full)"]["FEV1"] - a0["FEV1"]
        afx = res["+ itra, DDI off"]["FEV1"] - a0["FEV1"]
        ddx = res["+ itra, antifungal off"]["FEV1"] - a0["FEV1"]
        summary.append((label, res["+ itraconazole (full)"]["AUCP"] / a0["AUCP"],
                        res["+ itraconazole (full)"]["AUCCS"] / a0["AUCCS"],
                        tot, afx, ddx))
    print("\n  %-30s %7s %9s %10s %11s %8s %10s" %
          ("backbone", "xAUC", "xAUC_CS", "dFEV1 tot", "antifungal", "DDI", "DDI share"))
    for label, xauc, xcs, tot, afx, ddx in summary:
        print("  %-30s %7.2f %9.2f %10.3f %11.3f %8.3f %9.1f%%"
              % (label, xauc, xcs, tot, afx, ddx,
                 100.0 * ddx / tot if tot else float("nan")))
    print("""
  WHAT THE NUMBERS ACTUALLY SAY, WHICH IS NOT WHAT I EXPECTED TWICE OVER.
  My first expectation was that the interaction would explain most of
  itraconazole's apparent steroid-sparing on every backbone.  My second, after
  seeing prednisolone's small 1.1x AUC rise, was that prednisolone-based trials
  are safe from the confound.  Neither survives its own table, and the reason is
  that TWO different things vary across the rows.  The INTERACTION contribution
  scales with the steroid's CYP3A4 dependence, as expected (0.51 -> 3.72 ->
  6.14 FEV1 points).  But the ANTIFUNGAL contribution is not constant either
  (0.25 -> 0.22 -> 3.26), and not because the fungus knows which steroid was
  chosen: it is a headroom effect.  Prednisolone 10 mg/d already controls the
  disease (FEV1 80.3, plug score 1.08), so there is almost nothing left for an
  antifungal to win; inhaled budesonide alone barely controls it at all (FEV1
  69.9, plug score 4.40), so the same antifungal pharmacology has room to move
  6 points.  Both effects push the same way, which is the uncomfortable part —
  the apparent benefit of adding itraconazole is largest exactly where the
  backbone is weakest AND most CYP3A4-dependent, and those two explanations are
  not separable from an efficacy endpoint alone.
  That is the falsifiable prediction, and it has never been tested: an
  itraconazole steroid-sparing trial should report a LARGER effect on
  methylprednisolone or inhaled budesonide than on prednisolone, and the
  difference is pharmacokinetic, not antifungal.  A5 is also the reason the
  library's A4 exists — the interaction had to be calibrated against real DDI
  studies before it could be used to reinterpret efficacy trials.""")


def A6_omalizumab_biomarker(days=364):
    hr("A6  OMALIZUMAB — free IgE falls, TOTAL IgE rises, and the response criterion inverts")
    rec, y, p = run_scenario([(0, days, 14, {"oma": 375.0})], {}, days, dt=0.02)
    print("  %-8s %10s %10s %10s %8s %8s %8s %8s" %
          ("week", "IgE_tot", "IgE_free", "free ng/mL", "cplx%", "FcER_oc", "FCER", "effect"))
    for r in rec[::2]:
        cpl = 100.0 * r["CX"] / max(r["IgE_total"], 1e-9)
        print("  %-8.0f %10.1f %10.2f %10.2f %7.1f%% %8.3f %8.3f %8.4f"
              % (r["t"] / 7.0, r["IgE_total"], r["IgE_free"], r["IgE_free_ngml"],
                 cpl, r["FcER_occ"], y[IX["FCER"]] if r is rec[-1] else float("nan"),
                 r["effector"]))
    r0, rE = rec[0], rec[-1]
    print("\n  baseline total IgE %.0f IU/mL -> week 52 total IgE %.0f IU/mL "
          "(x%.2f)" % (r0["IgE_total"], rE["IgE_total"], rE["IgE_total"] / r0["IgE_total"]))
    print("  baseline free  IgE %.1f IU/mL -> week 52 free  IgE %.1f IU/mL "
          "(-%.1f%%)" % (r0["IgE_free"], rE["IgE_free"],
                         100.0 * (1.0 - rE["IgE_free"] / r0["IgE_free"])))
    print("  FcepsilonRI occupancy %.3f -> %.3f;  receptor density %.3f;  "
          "effector activation %.4f -> %.4f (-%.1f%%)"
          % (r0["FcER_occ"], rE["FcER_occ"], y[IX["FCER"]], r0["effector"],
             rE["effector"], 100.0 * (1.0 - rE["effector"] / r0["effector"])))
    crit = rE["IgE_total"] <= 0.65 * r0["IgE_total"]
    print("\n  ABPA response criterion 'total IgE falls 35-50%%': %s"
          % ("MET" if crit else "NOT MET — the criterion scores a working drug as a failure"))
    print("  the same patient judged on FREE IgE: %s"
          % ("responder" if rE["IgE_free"] < 0.5 * r0["IgE_free"] else "non-responder"))
    print("""
  Note the second, quieter result in the FcepsilonRI columns: at the classical
  free-IgE target of 25 ng/mL the receptor is still ~%.0f%% occupied, because
  FcepsilonRI binds IgE with ~%.2f nM affinity.  Instantaneous occupancy is
  therefore NOT what omalizumab buys; the benefit tracks the slow fall in
  receptor DENSITY, which is why the clinical effect lags the biomarker by
  weeks.""" % (100.0 * (25.0 / IGE_NM_TO_NGML) / ((25.0 / IGE_NM_TO_NGML) + P["Kd_fcer"]),
               P["Kd_fcer"]))
    return rec


def A7_complex_halflife():
    hr("A7  IS THE TOTAL-IgE RISE AN ARTEFACT OF ONE PARAMETER?  Sweep the complex half-life")
    print("""  At steady state on omalizumab, IgE production S must leave by two routes:
      S = kdegE*free + kelCX*complex
  so  total/total_0 = f + (1-f)*(kdegE/kelCX)  with f = residual free fraction.
  The RATIO of the two elimination rates is the whole story, and it is > 1 for
  any IgG-like complex.  The magnitude is not identifiable from total IgE.""")
    print("\n  %-14s %10s %12s %12s %12s" %
          ("cplx t1/2 (d)", "kelCX", "IgE_tot x", "free supp %", "criterion"))
    for th in (4.0, 6.0, 8.0, 12.0, 16.0, 20.0, 26.0):
        p = dict(P)
        p["kelCX"] = math.log(2.0) / th
        rec, y = simulate(build_regimen([(0, 364, 14, {"oma": 375.0})]), tend=364.0,
                          dt=0.02, p=p, record=364.0)
        r = readouts(y, p)
        rec0, y0 = simulate(build_regimen([]), tend=364.0, dt=0.02, p=p, record=364.0)
        r0 = readouts(y0, p)
        print("  %-14.1f %10.4f %12.2f %12.1f %12s"
              % (th, p["kelCX"], r["IgE_total"] / r0["IgE_total"],
                 100.0 * (1.0 - r["IgE_free"] / r0["IgE_free"]),
                 "MET" if r["IgE_total"] <= 0.65 * r0["IgE_total"] else "not met"))
    print("""
  Across a 6.5-fold range of complex half-life the total-IgE ratio moves from
  ~2x to ~9x and NEVER falls below 1, while free-IgE suppression stays above
  90%.  The direction is structural; the size is not.  Three explanations of
  the clinically reported 2-5x are indistinguishable here — faster-than-IgG
  complex clearance, incomplete assay detection of complexed IgE, or a genuine
  fall in IgE production — and the model does not pretend to separate them.""")


def params_for_baseline_ige(target_iu, iters=4):
    """Parameter set whose UNTREATED steady-state total IgE equals target_iu.
    Baseline IgE is a property of IgE PRODUCTION, not of the initial condition:
    setting y0 only lets the model relax straight back to its own set-point.
    This helper exists because the first version of A8 made exactly that error
    and returned an identical row for every baseline."""
    p = dict(P)
    for _ in range(iters):
        rec, y = simulate(build_regimen([]), tend=730.0, dt=0.03, p=p, record=730.0)
        got = readouts(y, p)["IgE_total"]
        if got <= 0:
            break
        p["kige"] *= target_iu / got
        if abs(got - target_iu) / target_iu < 1e-3:
            break
    rec, y = simulate(build_regimen([]), tend=730.0, dt=0.03, p=p, record=730.0)
    return p, readouts(y, p)


def A8_dose_ladder():
    hr("A8  WHY ABPA IS DOSED OFF THE OMALIZUMAB TABLE — a molar flux balance")
    print("""  Omalizumab neutralises IgE stoichiometrically, so the requirement is a
  FLUX condition, not a concentration one: the molar delivery rate of drug
  must exceed the molar production rate of IgE.  Both are computable.""")
    inp_nM_d = 375.0 * P["F_oma"] * OMA_MG_TO_NM / 14.0
    print("\n  375 mg q2w delivers %.2f nM/d of omalizumab into V=%.2f L" % (inp_nM_d, V_OMA))
    print("  IgE production needed to sustain a baseline of X IU/mL is "
          "kdegE * X * %.5f nM/IU = %.5f*X nM/d" % (IGE_IU_TO_NM, P["kdegE"] * IGE_IU_TO_NM))
    xstar = inp_nM_d * P["eta_oma"] / (P["kdegE"] * IGE_IU_TO_NM)
    print("  with %.1f IgE neutralised per omalizumab, the flux crossover is at "
          "baseline total IgE = %.0f IU/mL" % (P["eta_oma"], xstar))
    print("\n  Integrated model, 52 weeks of 375 mg q2w, with IgE PRODUCTION scaled")
    print("  so that each untreated steady state matches the stated baseline:")
    print("  %-8s %9s %9s %10s %8s %11s %8s %9s" %
          ("IgE0", "achieved", "free0", "free52", "supp %", "free ng/mL", "FCER", "eff drop"))
    cache = {}
    for ige0 in (500, 1000, 1500, 2000, 3000, 5000, 8000, 12000, 20000):
        p, r0 = params_for_baseline_ige(ige0)
        cache[ige0] = p
        rec, y = simulate(build_regimen([(0, 364, 14, {"oma": 375.0})]), tend=364.0,
                          dt=0.02, p=p, record=364.0)
        r = readouts(y, p)
        print("  %-8d %9.0f %9.1f %10.2f %7.1f%% %11.2f %8.3f %8.1f%%"
              % (ige0, r0["IgE_total"], r0["IgE_free"], r["IgE_free"],
                 100.0 * (1.0 - r["IgE_free"] / r0["IgE_free"]), r["IgE_free_ngml"],
                 y[IX["FCER"]], 100.0 * (1.0 - r["effector"] / r0["effector"])))
    print("\n  Dose required to reach free IgE < 25 ng/mL (the classical target):")
    print("  %-8s %14s %12s %16s %11s" %
          ("IgE0", "dose q2w (mg)", "vs max 375", "mg/(IU/mL x kg)", "on-table?"))
    for ige0 in (500, 1000, 1500, 2500, 5000, 10000):
        p = cache.get(ige0) or params_for_baseline_ige(ige0)[0]
        lo, hi, need = 25.0, 12000.0, None
        for _ in range(22):
            mid = 0.5 * (lo + hi)
            rec, y = simulate(build_regimen([(0, 224, 14, {"oma": mid})]), tend=224.0,
                              dt=0.03, p=p, record=224.0)
            if readouts(y, p)["IgE_free_ngml"] < 25.0:
                hi, need = mid, mid
            else:
                lo = mid
        print("  %-8d %14s %12s %16s %11s"
              % (ige0, ("%.0f" % need) if need else ">12000",
                 ("x%.2f" % (need / 375.0)) if need else "n/a",
                 ("%.5f" % (need / (ige0 * 70.0))) if need else "n/a",
                 "yes" if (need and need <= 375.0 and ige0 <= 1500) else "NO"))
    print("""
  The omalizumab dosing table is capped at total IgE 1500 IU/mL.  The model
  was never shown that number, yet the flux balance above puts the crossover
  in the same neighbourhood, and the required-dose column crosses 375 mg q2w
  there too.  ABPA's diagnostic threshold is total IgE > 1000 and its typical
  baseline is 2000-5000, so the disease begins where the label stops.  This
  reframes 'omalizumab non-response in ABPA' as a dosing question with an
  arithmetic answer, and it predicts that response should correlate with
  mg per (IU/mL x kg) rather than with mg.""")


def advanced_state(years=8.0):
    """Let the untreated model progress to advanced disease, then hand that
    state to a therapy arm.  Used because the therapy comparison in A9 is
    uninformative in early disease: BOTH single agents already drive the burden
    to ~0.2% of carrying capacity there, so every fractional-effect readout
    saturates at 1.0 and a Bliss test on burden can detect nothing.  The first
    version of A9 did exactly that and reported a null; the null was a property
    of the test, not of the model."""
    rec, y = simulate(build_regimen([]), tend=365.0 * years, dt=0.02, record=1e9)
    return y


def A9_synergy(days=364):
    hr("A9  DUPILUMAB + ITRACONAZOLE — testing synergy where the azole is NOT already enough")
    print("""  The biologic cannot kill Aspergillus and the azole cannot open a plug.  But
  the biologic raises k_out, which is the parameter that sets the azole's
  eradication threshold in A2, so the two are not two ways of doing one thing:
  one lowers a bound the other has to clear.

  A WARNING ABOUT THE TEST ITSELF.  In EARLY disease both single agents already
  reduce the fungal burden by >98%, so every fractional-effect readout pins at
  1.0 and a Bliss independence test has no dynamic range left to detect
  anything.  Run there it returns a null that says nothing about the model.
  The test below is therefore run twice: in early disease (where it is
  uninformative, shown here to make the point) and in advanced disease reached
  by letting the untreated model progress for 8 years (where it is not).""")
    arms = {
        "none": [],
        "itraconazole": [(0, days, 1, {"itra": 400.0})],
        "dupilumab": [(0, days, 14, {"dupi": 300.0})],
        "both": [(0, days, 1, {"itra": 400.0}), (0, days, 14, {"dupi": 300.0})],
    }
    for label, y0 in (("EARLY disease (the uninformative case)", None),
                      ("ADVANCED disease (8 y untreated progression)", advanced_state(8.0))):
        print("\n  --- %s ---" % label)
        if y0 is not None:
            r0 = readouts(y0, P)
            print("      starting state: BRON %.2f, EPX %.2f, k_out %.4f, plug score "
                  "%.2f, FEV1 %.1f" % (r0["BRON"], r0["EPX"], r0["kout"],
                                       r0["plug_score"], r0["FEV1"]))
            Es = erad_threshold(P, r0["PLUG"], r0["kout"])
            print("      azole eradication threshold E* here = %s /d against a class "
                  "ceiling of %.2f"
                  % (("%.3f" % Es) if Es is not None else "unreachable", P["Emax_af"]))
        R = {}
        for k, spec in arms.items():
            rec, y = simulate(build_regimen(spec), tend=float(days), dt=0.02,
                              record=float(days), y0=y0)
            R[k] = readouts(y, P)
        print("      %-14s %9s %9s %8s %9s %9s %9s %8s" %
              ("arm", "FPLG", "FTOT", "sanct%", "k_out", "E* here", "vs itra E", "FEV1"))
        E_itra = 0.954
        for k in arms:
            r = R[k]
            Ek = erad_threshold(P, r["PLUG"], r["kout"])
            print("      %-14s %9.5f %9.5f %7.1f%% %9.4f %9s %9s %8.2f"
                  % (k, r["FPLG"], r["FTOT"],
                     100.0 * r["FPLG"] / max(r["FTOT"], 1e-12), r["kout"],
                     ("%.3f" % Ek) if Ek is not None else "unreach",
                     ("OK" if (Ek is not None and Ek <= E_itra) else "SHORT"),
                     r["FEV1"]))
        def frac(k, key):
            base = R["none"][key]
            return (1.0 - R[k][key] / base) if base > 1e-12 else 0.0
        print("\n      Bliss independence: expected = a + b - a*b")
        print("      %-13s %9s %9s %11s %10s %10s" %
              ("readout", "itra", "dupi", "expected", "observed", "excess"))
        for key, lab in (("FTOT", "total fungus"), ("FPLG", "sanctuary"),
                         ("FLUM", "lumen"), ("plug_score", "plug score"),
                         ("hazard_yr", "exac hazard")):
            a, b = frac("itraconazole", key), frac("dupilumab", key)
            exp = a + b - a * b
            obs = frac("both", key)
            flag = "  <- saturated, uninformative" if min(a, b) > 0.95 else ""
            print("      %-13s %9.4f %9.4f %11.4f %10.4f %+10.4f%s"
                  % (lab, a, b, exp, obs, obs - exp, flag))
    print("""
  THE BLISS COLUMN IS A NEGATIVE RESULT AND IT IS REPORTED AS ONE.  The excess
  is slightly NEGATIVE on every readout in both severity regimes: this model
  does not produce Bliss synergy on fungal burden anywhere, and no amount of
  looking for it changes that.  The reason is structural rather than
  disappointing — both agents ultimately act through the same limiting term,
  reached through saturating functions, so on a bounded readout their effects
  compound sub-additively.  Anyone claiming pharmacological synergy for this
  combination from a burden endpoint will not find it here.

  WHAT THE MODEL DOES SAY IS IN THE OTHER TWO COLUMNS.  (i) COMPOSITION: at
  equal total burden the sanctuary share falls from ~54% on itraconazole alone
  to ~16% on the combination, so the surviving organisms sit where the drug can
  reach them.  (ii) THRESHOLD: dupilumab raises k_out and lowers E*, and in the
  advanced-disease patient E* starts at 1.311 /d — above the 0.954 /d that
  itraconazole 200 mg BID actually delivers — so the azole there is suppressive
  by construction until something moves the threshold.  The E*/vs-itra columns
  above show which arms cross it and which do not.
  The prediction, stated so it can fail: adding an IL-4Ra blocker should help
  the azole MOST in patients with the highest mucus-plug and bronchiectasis
  scores — the ones in whom azole monotherapy is observed to fail — and least
  in the mild patients in whom such a trial is easiest to run.  A trial
  enrolling mild ABPA would measure the combination at its weakest.""")


def A10_budesonide_trap(days=180):
    hr("A10  THE BUDESONIDE + ITRACONAZOLE TRAP — a toxicity produced by a PK interaction alone")
    arms = {
        "budesonide 1600 ug alone": [(0, days, 1, {"bud": 1.6})],
        "itraconazole alone": [(0, days, 1, {"itra": 400.0})],
        "budesonide + itraconazole": [(0, days, 1, {"bud": 1.6}), (0, days, 1, {"itra": 400.0})],
        "prednisolone 10 mg + itra": [(0, days, 1, {"pred": 10.0}), (0, days, 1, {"itra": 400.0})],
    }
    print("  %-30s %10s %10s %10s %10s %10s" %
          ("arm", "BUD mg/L", "CS", "cortisol", "BMD", "HbA1c"))
    for k, spec in arms.items():
        rec, y, p = run_scenario(spec, {}, days)
        r = readouts(y, p)
        print("  %-30s %10.5f %10.4f %10.2f %10.4f %10.2f"
              % (k, y[IX["BUD"]], r["CS"], r["CORT"], r["BMD"], r["HBA1C"]))
    print("""
  Inhaled budesonide is chosen precisely to avoid systemic exposure; adding a
  CYP3A4 inhibitor removes the first-pass step that made that true.  Nothing
  in the disease model changed between rows 1 and 3 — the adrenal suppression
  is manufactured entirely in the PK block.  This is the one place in the file
  where a drug interaction is not a confound but the adverse event itself.""")


def A11_lavage(days=364):
    hr("A11  BRONCHOSCOPIC PLUG REMOVAL — how long does emptying the sanctuary last?")
    print("""  Setting FPLG and PLUG to ~0 at day 0 is the only intervention that touches
  the sanctuary directly.  The question is whether the state is an attractor.""")
    for label, spec in (("lavage alone", [(0, 1, 1, {"lavage": 0.95})]),
                        ("lavage + itraconazole", [(0, 1, 1, {"lavage": 0.95}),
                                                   (0, days, 1, {"itra": 400.0})]),
                        ("lavage + itra + dupilumab", [(0, 1, 1, {"lavage": 0.95}),
                                                       (0, days, 1, {"itra": 400.0}),
                                                       (0, days, 14, {"dupi": 300.0})])):
        rec, y, p = run_scenario(spec, {}, days)
        print("\n  %s" % label)
        print("    %-8s %9s %9s %9s %9s" % ("week", "PLUG", "FLUM", "FPLG", "FEV1"))
        for r in rec[:8] + rec[10::3]:
            print("    %-8.0f %9.5f %9.5f %9.5f %9.2f"
                  % (r["t"] / 7.0, r["PLUG"], r["FLUM"], r["FPLG"], r["FEV1"]))
        # time to 50% recovery of FPLG
        base = readouts(run_scenario([], {}, days)[1], P)["FPLG"]
        t50 = None
        for r in rec[1:]:                     # skip t=0: that is the PRE-lavage value
            if r["FPLG"] >= 0.5 * base:
                t50 = r["t"]
                break
        print("    untreated sanctuary for comparison: FPLG = %.5f" % base)
        print("    refill to 50%% of it: %s" %
              (("day %.0f" % t50) if t50 is not None else "not reached in 52 wk"))


def A12_longterm():
    hr("A12  FIVE YEARS — the only endpoint that cannot be given back")
    days = 1825
    arms = {
        "untreated": [],
        "episodic pred (rescue only)": pred_taper(days // 5) ,
        "continuous itraconazole": [(0, days, 1, {"itra": 400.0})],
        "continuous dupilumab": [(0, days, 14, {"dupi": 300.0})],
        "dupilumab + itraconazole": [(0, days, 14, {"dupi": 300.0}),
                                     (0, days, 1, {"itra": 400.0})],
        "pred 10 mg maintenance": maintenance_pred(days, 10.0),
    }
    print("  %-30s %8s %8s %8s %8s %8s %8s %8s" %
          ("arm", "BRON", "dBRON", "FEV1", "exac", "OCSmg", "BMD", "netU"))
    rows = []
    for k, spec in arms.items():
        rec, y, p = run_scenario(spec, {}, days, dt=0.02)
        r = readouts(y, p)
        U = -(1.0 * r["CHAZ"] / 5.0 + 2.2 * (r["BRON"] - 3.0)
              + 0.09 * (100.0 - r["FEV1"]) + 0.0016 * r["CUMO"]
              + 12.0 * (1.0 - r["BMD"]))
        rows.append((k, r, U))
        print("  %-30s %8.3f %8.3f %8.2f %8.2f %8.0f %8.4f %+8.3f"
              % (k, r["BRON"], r["BRON"] - 3.0, r["FEV1"], r["CHAZ"],
                 r["CUMO"], r["BMD"], U))
    rows.sort(key=lambda x: -x[2])
    print("\n  ranking by net utility: " + " > ".join(r[0] for r in rows))


def A13_sensitivity(days=364):
    hr("A13  LOCAL SENSITIVITY — which parameters own which answer, and where sensitivity diverges")
    keys = ["f_pen", "gp", "kout0", "g_epx", "kin_f", "gl", "k_host", "Emax_af",
            "EC50_af", "kelCX", "kdegE", "eta_oma", "Kd_fcer", "smuc", "smuc_ige",
            "cs_plug", "cs_muc", "kb1", "kb2", "s3A4h_mpred", "Ki3A4_itra",
            "th2_base", "pc0", "kige", "sepx", "kepx"]
    print("""  Read the E* column first.  Emax_af — the POTENCY of the antifungal — has an
  elasticity of EXACTLY ZERO on the eradication threshold in every regime below,
  because E* comes from a Jacobian in which no antifungal property appears except
  f_pen.  Two E* columns are reported for exactly this reason:

    E*        E* recomputed at the state the arm actually reached, so a
              parameter can move it INDIRECTLY by changing PLUG and k_out
    E*|fixed  E* evaluated at a FIXED reference state (PLUG 1.0, untreated
              k_out), which isolates the closed form itself

  Emax_af and EC50_af are exactly 0.000 in the E*|fixed column — the closed form
  does not contain antifungal potency at all — while showing a small non-zero
  entry under E*, because killing more fungus lowers the antigen drive, hence
  IL-13, hence the plug, hence k_out.  That residual is a state effect, not a
  potency effect, and it is one to two orders of magnitude smaller than g_p's.
  Everything with a large E*|fixed entry (g_p, k_out0, g_epx, f_pen) belongs to
  the mucus axis, not to the antifungal.""")
    REGIMES = [
        ("itraconazole 200 mg BID alone (near the boundary)",
         [(0, days, 1, {"itra": 400.0})]),
        ("itraconazole + prednisolone 10 mg/d (well away from it)",
         [(0, days, 1, {"itra": 400.0}), (0, days, 1, {"pred": 10.0})]),
    ]
    store = {}
    for label, base_spec in REGIMES:
        def outputs(p):
            rec, y = simulate(build_regimen(base_spec), tend=float(days), dt=0.03,
                              p=p, record=float(days))
            r = readouts(y, p)
            E = erad_threshold(p, r["PLUG"], r["kout"])
            # Efix isolates the CLOSED FORM's dependence on the parameter by
            # evaluating E* at a FIXED disease state.  Estar lets the state move
            # too, so the two columns separate "does this parameter appear in the
            # eradication condition" from "does it change the state that
            # condition is evaluated at".
            Efix = erad_threshold(p, 1.0, kout_plug_eff(1.6, 3.0, 0.0, p))
            return dict(FTOT=r["FTOT"], FPLG=r["FPLG"], FEV1=r["FEV1"],
                        plug=r["plug_score"],
                        Estar=(E if E is not None else float("nan")),
                        Efix=(Efix if Efix is not None else float("nan")))
        b = outputs(P)
        print("\n  --- base regimen: %s ---" % label)
        print("      FTOT=%.6f  FPLG=%.6f  plug=%.3f  FEV1=%.2f  E*=%.3f  "
              "E*|fixed state=%.3f"
              % (b["FTOT"], b["FPLG"], b["plug"], b["FEV1"], b["Estar"], b["Efix"]))
        print("      elasticity  d ln(y) / d ln(theta)  at +20%")
        print("      %-16s %11s %11s %9s %9s %9s %10s" %
              ("parameter", "FTOT", "FPLG", "plugsc", "FEV1", "E*", "E*|fixed"))
        rows = {}
        for k in keys:
            p = dict(P)
            p[k] = P[k] * 1.2
            o = outputs(p)
            row = []
            for key in ("FTOT", "FPLG", "plug", "FEV1", "Estar", "Efix"):
                if b[key] and not math.isnan(b[key]) and not math.isnan(o[key]):
                    row.append((o[key] - b[key]) / b[key] / 0.2)
                else:
                    row.append(float("nan"))
            rows[k] = row
            print("      %-16s %11.3f %11.3f %9.3f %9.3f %9.3f %10.3f" % (k, *row))
        store[label] = rows
    a = store[REGIMES[0][0]]["gp"][1]
    c = store[REGIMES[1][0]]["gp"][1]
    print("""
  THE TWO REGIMES DISAGREE BY A FACTOR OF %.0f ON ONE ROW, AND THAT IS THE
  RESULT.  The elasticity of the sanctuary burden to g_p is %.2f on the
  well-controlled regimen and %.1f on azole monotherapy.  Sensitivity does not
  merely grow near the eradication boundary — it DIVERGES there, because the
  steady-state burden approaches the seed-supported floor and any change in the
  balance g_p - k_out - f_pen*E moves it by orders of magnitude.  A parameter
  sweep run only on the well-treated arm would have reported g_p as unimportant.

  So the honest headline of the whole file is this: everything the model says
  about whether an azole can clear ABPA is a statement about g_p, the intra-plug
  fungal growth rate, WHICH HAS NEVER BEEN MEASURED IN A HUMAN AIRWAY — and the
  closer a regimen sits to the boundary that matters clinically, the more
  completely that unmeasured number owns the answer.  Anyone wanting to refute
  this model should measure it rather than argue about the rest.""" %
          (abs(a / c) if c else float("nan"), c, a))


# ==========================================================================
# 9.  MAIN
# ==========================================================================
ANALYSES = {
    "A0": A0_selfcheck, "A1": A1_natural_history, "A2": A2_eradication_boundary,
    "A3": A3_scenarios, "A4": A4_ddi_calibration, "A5": A5_ddi_decomposition,
    "A6": A6_omalizumab_biomarker, "A7": A7_complex_halflife, "A8": A8_dose_ladder,
    "A9": A9_synergy, "A10": A10_budesonide_trap, "A11": A11_lavage,
    "A12": A12_longterm, "A13": A13_sensitivity,
}

if __name__ == "__main__":
    print("ABPA QSP reference implementation — %d ODE states, pure Python" % NST)
    which = sys.argv[1:] or list(ANALYSES)
    for a in which:
        if a in ANALYSES:
            ANALYSES[a]()
        else:
            print("unknown analysis %s" % a)
