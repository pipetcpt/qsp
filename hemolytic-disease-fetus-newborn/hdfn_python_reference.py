#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hemolytic Disease of the Fetus and Newborn (HDFN) -- executable reference model.

This file is the arithmetic of record.  `hdfn_mrgsolve_model.R` mirrors it
equation for equation; if the two ever disagree, THIS file is right, because
this is the one that was run.

--------------------------------------------------------------------------------
THE ONE SENTENCE
--------------------------------------------------------------------------------
Fetal red cells are destroyed at a rate that is a PRODUCT of three factors --

      destruction  =  A(t)          x   f_ag(t)        x   M(t)
                      antibody          antigen-positive    reticuloendothelial
                      delivered TO      FRACTION of the      clearance capacity
                      the fetus         circulating mass

-- and each therapeutic class in this disease owns exactly ONE of the three:

  * FcRn blockade (nipocalimab), plasmapheresis, high-dose IVIG competition
        -> factor A      (they change what arrives, not what it does)
  * intrauterine transfusion (IUT)
        -> factor f_ag   (donor cells are antigen-NEGATIVE, so transfusing does
                          not merely add haemoglobin, it DILUTES the substrate)
  * high-dose IVIG (FcgammaR blockade), splenic capacity
        -> factor M

The clinically famous number that falls out of factor 2 without being fitted is
the slowing of the post-transfusion haemoglobin decline: ~0.4 g/dL/day between
the first and second IUT, less between later ones.  Nobody has to fit that
sequence.  Once the fetal (antigen-positive) red-cell mass has been diluted
~2-fold by the first transfusion and ~4-fold by the second, the immune term has
less substrate to work on, and what remains of the decline is donor-cell
senescence plus growth dilution.  Section 3 of the analysis decomposes it.

--------------------------------------------------------------------------------
WHAT IS FITTED (six numbers) AND WHAT IS PREDICTED (everything else)
--------------------------------------------------------------------------------
FITTED
  1-2. v0, g_pl (placental conveyor)  <- Malek 1996 fetal:maternal IgG ratio
                                         0.075 at 19.5 wk and 1.25 at 39 wk
  3-4. kops, Kres (destruction)       <- 15 IU/mL anti-D reaches Hb 0.65 MoM by
                                         ~26 wk (UK quantitation risk bands);
                                         post-IUT decline 0.40 g/dL/day
                                         (Nishie 2012)
  5-6. ksens, dev (immune deviation)  <- 16% sensitisation without prophylaxis,
                                         1.6% with postpartum 300 ug only

PREDICTED (never shown to the fitter)
  * decline after the 2nd and 3rd IUT (target 0.30 / 0.20 g/dL/day)
  * that MCA-PSV 1.5 MoM is the arithmetic image of Hb 0.65 MoM
  * the 12% false-positive rate of that threshold
  * 0.1-0.2% residual sensitisation on antenatal + postpartum prophylaxis
  * the UNITY nipocalimab result (7/13 = 54% IUT-free live birth >= 32 wk)
  * the 20 ug anti-D per mL fetal RBC stoichiometry and the 72 h window
  * anti-K anaemia with LOW bilirubin, anti-D anaemia with high bilirubin
  * why a fetus with 8 g/dL of haemolysis is not jaundiced and a newborn is
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

# ==============================================================================
# 0.  GESTATIONAL REFERENCE FUNCTIONS
# ==============================================================================
# Hadlock 50th-centile estimated fetal weight, fitted here (log-cubic) so the
# growth curve in the model is the published curve and not an invention.
_EFW_GA = np.array([16, 20, 24, 28, 32, 36, 40], float)
_EFW_G = np.array([146.0, 331.0, 670.0, 1150.0, 1810.0, 2620.0, 3400.0])
_EFW_COEF = np.polyfit(_EFW_GA, np.log(_EFW_G), 3)


def efw_ref(ga):
    """Median estimated fetal weight (g) at gestational age ga (weeks)."""
    return float(np.exp(np.polyval(_EFW_COEF, np.clip(ga, 14.0, 42.0))))


def hb_ref(ga):
    """Median fetal haemoglobin (g/dL).  Nicolaides: 10.9 at 17 wk -> 14.5 at 40."""
    return 10.9 + 0.1565 * (np.clip(ga, 17.0, 42.0) - 17.0)


def hb_ref_pn(ga_birth, pna_d):
    """Postnatal haemoglobin reference (g/dL).  Cord Hb is ~14% above the
    intrauterine venous value; it then falls to the physiological nadir
    (11.0 g/dL at term, 9.0 at 32 wk) with a ~22 d time constant."""
    hb0 = hb_ref(ga_birth) * 1.14
    nadir = 9.0 + 0.25 * (min(ga_birth, 40.0) - 32.0)
    return nadir + (hb0 - nadir) * math.exp(-max(pna_d, 0.0) / 22.0)


def alb_ref(ga, p):
    """Median fetal/neonatal albumin (g/dL): ~1.9 at 20 wk, ~3.3 at term."""
    return p["alb0"] + p["alb_slope"] * (min(ga, 40.0) - 20.0)


def ugt_ref(pna, p):
    """UGT1A1 activity as a fraction of adult.  ~7.5% at birth, half-mature at
    ~26 d, adult by ~3 months -- the clock that makes neonatal, but not fetal,
    haemolysis a bilirubin problem."""
    x = (max(pna, 0.0) / p["ugt_t50"]) ** p["ugt_n"]
    return p["ugt_birth"] + (1.0 - p["ugt_birth"]) * x / (1.0 + x)


def psv_med(ga):
    """Median MCA peak systolic velocity (cm/s), Mari 2000: exp(2.31+0.046*GA)."""
    return math.exp(2.31 + 0.046 * np.clip(ga, 16.0, 42.0))


def vm_plasma(ga):
    """Maternal plasma volume (L): 2.5 L pre-pregnant -> +45% by the third
    trimester.  This term alone reproduces the observed fall of maternal total
    IgG to 60-70% of early-pregnancy values (Malek 1996) by dilution."""
    return 2.5 * (1.0 + 0.45 / (1.0 + math.exp(-(ga - 24.0) / 4.5)))


# ==============================================================================
# 1.  PARAMETERS
# ==============================================================================
P = dict(
    # ---- maternal IgG homeostasis (amounts in g, concentrations in g/L) ------
    igg0=10.0,          # g/L maternal total IgG at 12 wk
    k_int=0.120,        # /d pinocytotic internalisation of IgG
    phi_rescue=0.967,   # fraction of internalised IgG rescued when FcRn free
    K_igg=30.0,         # g/L IgG concentration at half-saturation of FcRn
    K_nip=8.0,          # mg/L nipocalimab at half-saturation of FcRn
    # ---- anti-D humoral response (amounts in IU) -----------------------------
    ad_iu_per_ug=5.0,   # 1 ug anti-D = 5 IU  (300 ug = 1500 IU)
    k_pc=0.05,          # /d plasma-cell turnover
    sec_pc=1.0,         # IU/d per unit plasma cell (scaled by titre init)
    k_bm=0.002,         # /d memory B-cell decay
    k_prime=0.60,       # per mL fetal RBC processed -> memory B-cell priming
    k_boost=0.35,       # /d differentiation of memory B -> plasma cells on Ag
    # ---- fetomaternal haemorrhage / prophylaxis ------------------------------
    k_sen_free=0.0099,  # /d loss of UNCOATED fetal D+ RBC from the maternal
    #                     circulation.  This is ~ln2/70 d, i.e. ordinary red
    #                     cell senescence, because an uncoated fetal cell is
    #                     just a red cell.  IT IS THE WHOLE REASON THE 72-HOUR
    #                     WINDOW IS 72 HOURS AND NOT 72 MINUTES: the antigen
    #                     sits there for months, so a three-day delay leaks
    #                     only ~3% of the exposure integral.
    k_clear_coat=6.0,   # /d clearance of RhIG-coated fetal D+ RBC (spleen)
    iu_per_ml_rbc=100.0,  # IU anti-D needed to coat 1 mL fetal RBC (=20 ug/mL)
    ksens=0.030,        # FITTED (set by calibrate)
    dev=0.98,           # STRUCTURAL: efficiency with which the coated-cell
    #                     clearance route avoids priming.  Not identifiable
    #                     from the two population risks, so it is fixed and
    #                     its influence is reported as a sensitivity.
    fmh_ante=0.016,     # FITTED mL of fetal whole blood in third-trimester
    #                     silent fetomaternal haemorrhage
    # ---- nipocalimab PK (2-cpt + saturable FcRn route) ----------------------
    nip_V1=3.0,         # L central (maternal plasma)
    nip_V2=4.0,         # L peripheral
    nip_Q=0.80,         # L/d intercompartmental
    nip_CLlin=0.85,     # L/d non-saturable clearance.  An FcRn blocker cannot
    #                     be rescued by the receptor it blocks, so its own
    #                     clearance is high and its half-life short (~5 d).
    nip_CLfcrn=0.25,    # L/d additional FcRn-mediated clearance at low conc
    nip_pl_pen=0.15,    # nipocalimab concentration in the placental
    #                     interstitium as a fraction of maternal plasma.  This
    #                     is why the drug can suppress maternal IgG by 70% and
    #                     still leave a fetus that needs transfusion: the
    #                     endothelial FcRn that sets catabolism sees the plasma
    #                     concentration, the syncytiotrophoblast FcRn that
    #                     moves antibody to the fetus sees a seventh of it.
    # ---- IVIG ---------------------------------------------------------------
    ivig_f=1.0,         # bioavailability of IV IgG into the maternal pool
    ivig_compete=1.0,   # switch: does IVIG-derived IgG compete for FcRn?  Set
    #                     to 0 to isolate the FcgammaR mechanism.  Zeroing the
    #                     IVIG POOL instead removes both mechanisms at once,
    #                     which is how the first version of this decomposition
    #                     reached the wrong conclusion.
    ivig_fcgr=0.55,     # max fractional blockade of fetal FcgammaR at 1 g/kg/wk
    ivig_K=6.0,         # g/L IVIG-derived IgG for half-maximal FcgammaR blockade
    # ---- placental conveyor -------------------------------------------------
    v0=0.020,           # FITTED transfer capacity at 20 wk
    g_pl=0.220,         # FITTED /wk exponential growth of transfer capacity
    sub3=0.60,          # IgG3 transfer relative to IgG1 (Malek: IgG1>IgG4>IgG3>IgG2)
    # ---- fetal IgG kinetics -------------------------------------------------
    vd_igg_f=2.0,       # fetal IgG distribution volume as multiple of plasma vol
    k_int_f=0.10,       # /d fetal IgG internalisation
    # ---- red cells ----------------------------------------------------------
    kbv=0.105,          # mL fetoplacental blood per g EFW at normal Hct
    #                   plasma volume is DERIVED as kbv*(1-Hct_ref(GA)) so that
    #                   a fetus with reference haemoglobin has exactly the
    #                   reference blood volume -- there is no free parameter here.
    #                     The circulation is written as plasma volume PLUS red
    #                     cell volume, not as a fixed blood volume, because a
    #                   The circulation is written as plasma volume PLUS red
    #                   cell volume, not as a fixed blood volume, because a
    #                   transfusion of packed cells adds cell volume and only
    #                   20% plasma -- which is why it does NOT dilute the
    #                   maternal antibody already inside the fetus.
    k_vol=0.60,         # /d regulation of plasma volume toward its target
    vpl_expand=0.50,    # fractional plasma-volume expansion at zero Hct
    mchc=32.8,          # g/dL mean corpuscular haemoglobin concentration
    bv_neo=92.0,        # mL blood per kg after birth (kbv*EFW/(bv_neo*wt)
    #                     = 1.14, which IS the cord-vs-intrauterine Hb step)
    hct_hb=3.05,        # Hct(%) per g/dL Hb (i.e. Hb = Hct/3.05)
    kops=0.020,         # FITTED /d per (IU/mL x site-density) opsonisation rate
    Kres=0.30,          # STRUCTURAL g opsonised Hb at half-maximal RES
    #                     clearance.  Not identifiable from the calibration
    #                     targets; fixed so that the splenic clearance is about
    #                     half-saturated in severe disease, and its influence
    #                     is reported as a sensitivity instead of fitted.
    vmax_res=1.0,       # g Hb/d per unit RES mass, cleared when unsaturated
    k_res_grow=0.10,    # /d expansion of splenic/hepatic macrophage mass
    res_max=4.0,        # x normal RES capacity attainable
    t_own=75.0,         # d fetal red-cell lifespan (senescence)
    t_don=70.0,         # d donor red-cell lifespan in the fetus/newborn
    k_mat=0.5,          # /d reticulocyte -> mature red cell
    site_D=1.0,         # relative D antigen density (R1r ~ 1.0; R2R2 ~ 3.0)
    pot3=1.7,           # ADCC potency of IgG3 relative to IgG1 per IU
    # ---- erythropoiesis -----------------------------------------------------
    epo0=15.0,          # mU/mL fetal EPO at normal Hb
    epo_n=5.1,          # EPO varies as (1/HbMoM)^n -- BIDIRECTIONAL, so that a
    #                     fetus made polycythaemic by a transfusion SUPPRESSES
    #                     its own erythropoiesis.  n = 5.1 gives the observed
    #                     ~35-fold EPO rise at Hb 0.5 MoM.
    epo_cap=40.0,       # maximal fold rise of EPO
    epo_pow=1.2,        # exponent linking EPO to progenitor recruitment
    k_epo=6.0,          # /d EPO elimination
    k_prog=0.55,        # /d progenitor recruitment at EPO0
    prog_max=5.0,       # maximal MARROW progenitor expansion (x normal)
    emh_prod=0.60,      # red cells made per unit hepatic extramedullary pool
    k_prog_d=0.481,     # /d progenitor loss (set so Prog=1 is the
    #                     steady state at normal EPO -- not a free parameter)
    k_prod=1.00,        # production as a multiple of the REQUIREMENT
    #                     (senescence + growth) at Prog = 1
    k_emh=0.09,         # /d recruitment of hepatic extramedullary erythropoiesis
    k_emh_off=0.06,     # /d resolution
    emh_thresh=4.0,     # progenitor expansion above which EMH is recruited
    kell_kill=0.0,      # /d per (IU/mL) killing of PROGENITORS (anti-K scenario)
    # ---- liver, albumin, hydrops -------------------------------------------
    alb0=1.90,          # g/dL fetal albumin at 20 wk ...
    alb_slope=0.070,    # ... rising 0.07 g/dL per week to ~3.3 at term
    k_alb_cat=0.115,    # /d albumin catabolism
    emh_alb=0.40,       # maximal fractional suppression of albumin synthesis
    k_portal=0.80,      # mmHg umbilical venous pressure per unit EMH
    cvp0=4.5,           # mmHg baseline umbilical venous pressure
    lp_s=6.0,           # mL/d/mmHg/kg fetal capillary filtration coefficient
    pi_frac=0.55,       # interstitial oncotic pressure as a fraction of plasma
    asc_compl=150.0,    # mL of ascites per mmHg of interstitial pressure
    p_i=1.0,            # mmHg baseline interstitial hydrostatic pressure
    sigma=0.85,         # oncotic reflection coefficient
    alb_pi=4.5,         # mmHg plasma oncotic pressure per g/dL albumin
    #                     (fetal colloid osmotic pressure ~11-14 mmHg)
    k_lymph_p=0.50,     # /mmHg rise in lymphatic return per mmHg of
    #                     interstitial pressure -- the "safety factor"
    lymph_cvp=0.16,     # /mmHg suppression of lymph return by venous pressure
    k_perm=0.35,        # increase in Lp*S per unit hypoxia index
    asc_hydrops=25.0,   # mL ascites that is sonographically overt at 26 wk
    # ---- oxygen, cardiac ----------------------------------------------------
    co_ref=450.0,       # mL/min/kg combined fetal cardiac output
    sao2=0.55,          # fetal arterial saturation
    do2_alpha=1.0,      # exponent of cerebral flow compensation (DERIVED, =1)
    visc_k=0.0,         # viscosity contribution to MCA velocity.  SET TO ZERO
    #                     ON PURPOSE: see analysis 4.  A non-zero value
    #                     double-counts, because the flow rise that viscosity
    #                     causes is the very mechanism by which oxygen
    #                     delivery is defended, and that is already the
    #                     do2_alpha term.
    co_max=2.20,        # maximal fetal cardiac-output reserve (x baseline).
    #                     With do2_alpha = 1 this single number FIXES the
    #                     haemoglobin at which the fetus can no longer defend
    #                     oxygen delivery: HbMoM = 1/2.2 = 0.455, i.e. ~5.7
    #                     g/dL at 26 wk -- which is where hydrops is observed.
    k_fail=0.22,        # /d accumulation of cardiac decompensation
    k_recov=0.12,       # /d recovery of cardiac reserve
    # ---- bilirubin ----------------------------------------------------------
    bil_per_g=34.0,     # mg bilirubin per g haemoglobin catabolised
    bil_shunt=2.0,      # mg/kg/d 'shunt' bilirubin (ineffective erythropoiesis,
    #                     haem turnover other than circulating red cells)
    k_bil_pl=14.0,      # /d fetal bilirubin clearance ACROSS THE PLACENTA
    #                     (set so normal fetal plasma bilirubin is ~0.5 mg/dL)
    vd_bp=0.09,         # L/kg rapidly mixing (plasma+equilibrating) bilirubin
    vd_bt=0.11,         # L/kg slowly mixing tissue bilirubin.  The SUM, 0.20
    #                     L/kg, is what a double-volume exchange transfusion
    #                     implies: it removes ~85% of the plasma pool but only
    #                     ~45% of the body pool, which is the rebound.
    k_bx=10.0,          # /d plasma <-> tissue bilirubin exchange
    k_amn=0.55,         # /d appearance of fetal bilirubin in amniotic fluid
    k_amn_out=0.42,     # /d amniotic-fluid bilirubin turnover
    od450_per_mg=0.052, # dOD450 units per mg/dL amniotic bilirubin equivalent
    ugt_birth=0.075,    # UGT1A1 activity at birth (fraction of adult) FITTED
    ugt_t50=25.5,       # d for half-maturation of UGT1A1 FITTED
    ugt_n=1.5,          # Hill coefficient of UGT1A1 maturation
    km_bil=8.0,         # mg/dL Michaelis constant of hepatic conjugation
    vmax_ugt=95.0,      # mg/d/kg conjugation at full (adult) activity
    k_ehc=0.22,         # /d enterohepatic recirculation
    k_gut=0.55,         # /d faecal loss of gut bilirubin
    k_photo=0.022,      # /d per (uW/cm2/nm) photoisomerisation of plasma bilirubin
    #                     (30 uW/cm2/nm -> ~0.4 mg/dL/h fall at TSB 15)
    k_iso_out=1.6,      # /d biliary loss of photoisomers
    k_iso_back=0.35,    # /d reversion of photoisomer to native bilirubin
    alb_bind=0.72,      # mg bilirubin bound per g/dL albumin (molar 1:1 => 8.5)
    bind_thresh=0.65,   # free-fraction-weighted B/A ratio at neurotoxic risk
    # ---- hazards ------------------------------------------------------------
    h0=0.00012,         # /d background fetal loss
    h_hydrops=0.010,    # /d additional hazard with overt hydrops
    h_severe=0.020,     # /d additional hazard at Hb < 0.35 MoM
    h_acid=0.030,       # /d additional hazard at pH < 7.15
    h_iut=0.012,        # per-procedure loss (Zwiers 2017: 1.2%/procedure)
    # ---- fetal acid-base ----------------------------------------------------
    do2_crit=0.55,      # fraction of reference DO2 below which lactate rises
    k_lac=2.2,          # mmol/L/d lactate accumulation at zero reserve
    k_lac_out=1.1,      # /d lactate clearance
    ph0=7.35,
    ph_lac=0.035,       # pH units per mmol/L lactate
)

# state vector -----------------------------------------------------------------
SV = [
    "Mig",    # 0  maternal total IgG (g)
    "Ma1",    # 1  maternal anti-D IgG1 (IU)
    "Ma3",    # 2  maternal anti-D IgG3 (IU)
    "Mbm",    # 3  maternal memory B-cell pool (rel)
    "Mpc",    # 4  maternal anti-D plasma cells (rel)
    "Mfree",  # 5  uncoated fetal D+ RBC in maternal blood (mL)
    "Mcoat",  # 6  RhIG-coated fetal D+ RBC in maternal blood (mL)
    "Mrhig",  # 7  passive (prophylactic) anti-D in mother (IU)
    "Mdep",   # 8  intramuscular RhIG depot (IU)
    "Sens",   # 9  cumulative primary-sensitisation signal (-log survival)
    "Nc",     # 10 nipocalimab central (mg)
    "Np",     # 11 nipocalimab peripheral (mg)
    "Miv",    # 12 IVIG-derived maternal IgG (g)
    "Fig",    # 13 fetal total IgG (g)
    "Fa1",    # 14 fetal anti-D IgG1 (IU)
    "Fa3",    # 15 fetal anti-D IgG3 (IU)
    "Fnip",   # 16 fetal nipocalimab (mg)
    "EFW",    # 17 estimated fetal weight (g)
    "Vpl",    # 18 fetal/neonatal PLASMA volume (mL)
    "Ro",     # 19 own antigen-POSITIVE red-cell Hb mass (g)
    "Rd",     # 20 donor antigen-NEGATIVE red-cell Hb mass (g)
    "Rr",     # 21 own reticulocyte Hb mass (g)
    "Rop",    # 22 opsonised red-cell Hb mass awaiting clearance (g)
    "Prog",   # 23 erythroid progenitor pool (rel)
    "EMH",    # 24 hepatic extramedullary erythropoiesis (rel)
    "Epo",    # 25 fetal erythropoietin (mU/mL)
    "Res",    # 26 reticuloendothelial (macrophage) capacity (rel)
    "Alb",    # 27 fetal/neonatal albumin CONCENTRATION (g/dL)
    "Asc",    # 28 ascites + interstitial oedema volume (mL)
    "Card",   # 29 cardiac decompensation index (0-1)
    "Lac",    # 30 lactate (mmol/L)
    "Bil",    # 31 plasma bilirubin (mg)
    "Bex",    # 32 extravascular tissue bilirubin (mg)
    "Bgut",   # 33 intestinal bilirubin available for recirculation (mg)
    "Biso",   # 34 photoisomer pool (mg)
    "Amn",    # 35 amniotic-fluid bilirubin equivalent (mg)
    "Ugt",    # 36 UGT1A1 activity (fraction of adult)
    "Bind",   # 37 cumulative free-bilirubin neurotoxic exposure (index-d)
    "Hzd",    # 38 cumulative perinatal-loss hazard
    "Fer",    # 39 transfusional iron load (mg Fe)
    "Nrbc",   # 40 circulating nucleated red cells (rel)
    "Wgt",    # 41 postnatal weight (kg) after birth
    "Adose",  # 42 cumulative anti-D delivered to fetus (IU, bookkeeping)
    "Dest",   # 43 cumulative antigen-positive Hb destroyed (g, bookkeeping)
    "SenC",   # 44 cumulative senescent Hb loss (g, bookkeeping)
]
IX = {k: i for i, k in enumerate(SV)}
NS = len(SV)


class Ctl:
    """Everything the clinician (or the trial protocol) does."""

    def __init__(self, **kw):
        self.born = False
        self.ga_birth = None
        self.nip_dose = 0.0        # mg/kg/week
        self.nip_start = 99.0      # weeks
        self.nip_stop = 35.0
        self.ivig_dose = 0.0       # g/kg/week
        self.ivig_start = 99.0
        self.photo = 0.0           # uW/cm2/nm
        self.mat_wt = 70.0         # kg
        self.iut_n = 0
        self.iut_log = []          # (ga, hb_pre, hb_post, volume)
        self.ehc_on = 0.0          # enterohepatic circulation (neonate feeding)
        self.__dict__.update(kw)


# ==============================================================================
# 2.  ALGEBRAIC READOUTS
# ==============================================================================
def readouts(t, y, p, c):
    ga = t / 7.0
    born = c.born and (c.ga_birth is not None) and (ga >= c.ga_birth)
    # ---- maternal ------------------------------------------------------------
    Vm = vm_plasma(min(ga, c.ga_birth if c.ga_birth else ga))
    Cig = y[IX["Mig"]] / Vm                      # g/L
    Civ = y[IX["Miv"]] / Vm                      # g/L
    Cnip = y[IX["Nc"]] / p["nip_V1"]             # mg/L
    Ca1 = y[IX["Ma1"]] / (Vm * 1000.0)           # IU/mL
    Ca3 = y[IX["Ma3"]] / (Vm * 1000.0)           # IU/mL
    Crh = y[IX["Mrhig"]] / (Vm * 1000.0)         # IU/mL
    # free FcRn fraction (mother) -- IgG, IVIG and nipocalimab all compete
    fu_m = 1.0 / (1.0 + (Cig + p["ivig_compete"] * Civ) / p["K_igg"]
                  + Cnip / p["K_nip"])
    # ---- fetus / newborn -----------------------------------------------------
    efw = max(y[IX["EFW"]], 80.0)
    wt = y[IX["Wgt"]] if born else efw / 1000.0
    pna = (ga - c.ga_birth) * 7.0 if born else 0.0
    Vrbc = (y[IX["Ro"]] + y[IX["Rd"]] + y[IX["Rr"]]) / (p["mchc"] / 100.0)
    Vpl = max(y[IX["Vpl"]], 1.0)
    Vfp = Vpl + Vrbc                             # mL whole blood
    hbmass = y[IX["Ro"]] + y[IX["Rd"]] + y[IX["Rr"]]
    hb = 100.0 * hbmass / Vfp                    # g/dL
    hct = Vrbc / Vfp                             # fraction
    Vd_igg = p["vd_igg_f"] * Vpl / 1000.0        # L
    Cf_ig = y[IX["Fig"]] / max(Vd_igg, 1e-6)     # g/L
    Cf_a1 = y[IX["Fa1"]] / max(Vd_igg * 1000.0, 1e-6)   # IU/mL
    Cf_a3 = y[IX["Fa3"]] / max(Vd_igg * 1000.0, 1e-6)
    Cf_nip = y[IX["Fnip"]] / max(Vd_igg, 1e-6)
    alb = y[IX["Alb"]]                           # g/dL, tracked directly
    # ---- factor 1: antibody actually present in the fetus -------------------
    A_eff = (Cf_a1 + p["pot3"] * Cf_a3) * p["site_D"]
    # ---- factor 2: antigen-positive fraction --------------------------------
    f_ag = (y[IX["Ro"]] + y[IX["Rr"]]) / max(hbmass, 1e-9)
    # ---- factor 3: clearance capacity ---------------------------------------
    fcgr_block = p["ivig_fcgr"] * Civ / (p["ivig_K"] + Civ)
    M_eff = y[IX["Res"]] * (1.0 - fcgr_block)
    # ---- oxygen delivery and the Doppler readout ----------------------------
    hbr = hb_ref_pn(c.ga_birth, pna) if born else hb_ref(ga)
    hbmom = hb / hbr
    # cerebral flow is what the Doppler sees, and it is defended to hold
    # oxygen delivery constant -- it is NOT capped by the cardiac reserve,
    # because the brain is preferentially perfused.
    demand_brain = (1.0 / max(hbmom, 0.12)) ** p["do2_alpha"]
    # total cardiac output IS capped; the gap between what the fetus needs and
    # what the heart can deliver is what decompensates it.
    demand = min(demand_brain, p["co_max"])
    co = p["co_ref"] * demand
    do2 = demand * hbmom                      # relative oxygen delivery
    do2_abs = p["co_ref"] * demand * hb * 1.34 * p["sao2"] / 100.0  # mL O2/min/kg
    hct_ref = hbr / p["mchc"]
    visc = math.exp(-p["visc_k"] * (hct - hct_ref))
    psv = psv_med(ga) * demand_brain * visc
    psv_mom = psv / psv_med(ga)
    # ---- hydrops -------------------------------------------------------------
    cvp = p["cvp0"] + p["k_portal"] * y[IX["EMH"]] + 6.0 * y[IX["Card"]]
    pi_p = p["alb_pi"] * alb
    hypox = max(0.0, 1.0 - do2 / p["do2_crit"]) if do2 < p["do2_crit"] else 0.0
    kg = max(wt, 0.05)
    lp = p["lp_s"] * kg * (1.0 + p["k_perm"] * hypox)
    pc = 0.55 * cvp + 8.0
    p_i_eff = p["p_i"] + y[IX["Asc"]] / p["asc_compl"]
    drive = (pc - p_i_eff) - p["sigma"] * pi_p * (1.0 - p["pi_frac"])
    jv = lp * drive
    # THE LYMPHATIC CAPACITY IS DEFINED AS THE BASELINE FILTRATION.  Writing it
    # this way means the model makes no claim about the absolute filtration rate
    # -- only about DEPARTURES from balance -- so a healthy fetus never drifts
    # into ascites, and hydrops is forced to come from hypoalbuminaemia, venous
    # pressure or permeability rather than from a parameter mismatch.
    pc0 = 0.55 * p["cvp0"] + 8.0
    pi_ref = p["alb_pi"] * alb_ref(ga if not born else c.ga_birth, p)
    jv_ref = p["lp_s"] * kg * ((pc0 - p["p_i"]) -
                               p["sigma"] * pi_ref * (1.0 - p["pi_frac"]))
    lymph = jv_ref * (1.0 + p["k_lymph_p"] * (p_i_eff - p["p_i"])) * \
        math.exp(-p["lymph_cvp"] * max(0.0, cvp - p["cvp0"]))
    asc_thr = p["asc_hydrops"] * efw / 900.0
    # ---- bilirubin -----------------------------------------------------------
    vd_b = p["vd_bp"] * wt                       # L rapidly mixing pool
    tsb = y[IX["Bil"]] / max(vd_b * 10.0, 1e-9)  # mg/dL  (mg / (L*10 dL))
    ctis = y[IX["Bex"]] / max(p["vd_bt"] * wt * 10.0, 1e-9)
    bcap = p["alb_bind"] * alb * 10.0 * vd_b     # mg bindable
    bfree_ratio = max(0.0, (y[IX["Bil"]] - bcap)) / max(y[IX["Bil"]], 1e-9)
    ba_ratio = tsb / max(alb, 0.3) / 8.5         # molar bilirubin/albumin ratio
    od450 = p["od450_per_mg"] * y[IX["Amn"]] / max(0.4 * efw / 1000.0, 0.05)
    return dict(ga=ga, born=born, Vm=Vm, Cig=Cig, Civ=Civ, Cnip=Cnip, Ca1=Ca1,
                Ca3=Ca3, Crh=Crh, fu_m=fu_m, efw=efw, wt=wt, Vfp=Vfp, hb=hb,
                hct=hct, Vpl=Vpl, Vrbc=Vrbc, Cf_ig=Cf_ig, Cf_a1=Cf_a1, Cf_a3=Cf_a3,
                Cf_nip=Cf_nip, alb=alb, A_eff=A_eff, f_ag=f_ag, M_eff=M_eff,
                fcgr_block=fcgr_block, hbr=hbr, hbmom=hbmom, co=co, do2=do2,
                psv=psv, psv_mom=psv_mom, do2_abs=do2_abs, cvp=cvp, jv=jv, lymph=lymph,
                jv_ref=jv_ref, drive=drive, asc_thr=asc_thr, hydrops=y[IX["Asc"]] > asc_thr, tsb=tsb,
                pna=pna, demand_brain=demand_brain,
                bfree_ratio=bfree_ratio, ba_ratio=ba_ratio, od450=od450, ctis=ctis,
                hypox=hypox, demand=demand, hbmass=hbmass, vd_b=vd_b)


# ==============================================================================
# 3.  RIGHT-HAND SIDE
# ==============================================================================
def rhs(t, y, p, c):
    y = np.maximum(y, 0.0)
    r = readouts(t, y, p, c)
    d = np.zeros(NS)
    ga, born = r["ga"], r["born"]

    # ---------------- maternal IgG and anti-D --------------------------------
    kel_igg = p["k_int"] * (1.0 - p["phi_rescue"] * r["fu_m"])
    syn_igg = p["igg0"] * vm_plasma(12.0) * p["k_int"] * \
        (1.0 - p["phi_rescue"] / (1.0 + p["igg0"] / p["K_igg"]))
    d[IX["Mig"]] = syn_igg - kel_igg * y[IX["Mig"]]
    d[IX["Miv"]] = -kel_igg * y[IX["Miv"]]
    # anti-D: secretion by plasma cells, catabolism shared with bulk IgG
    d[IX["Ma1"]] = p["sec_pc"] * y[IX["Mpc"]] - kel_igg * y[IX["Ma1"]]
    d[IX["Ma3"]] = 0.18 * p["sec_pc"] * y[IX["Mpc"]] - kel_igg * y[IX["Ma3"]]
    d[IX["Mrhig"]] = 0.35 * y[IX["Mdep"]] - kel_igg * y[IX["Mrhig"]] \
        - p["k_clear_coat"] * y[IX["Mcoat"]] * p["iu_per_ml_rbc"]
    d[IX["Mdep"]] = -0.35 * y[IX["Mdep"]]

    # ---------------- fetomaternal haemorrhage / sensitisation ---------------
    # coating fraction: available passive + active anti-D vs antigen on offer
    need = p["iu_per_ml_rbc"] * (y[IX["Mfree"]] + y[IX["Mcoat"]])
    avail = y[IX["Mrhig"]] + y[IX["Ma1"]] + y[IX["Ma3"]]
    coat_f = avail / (avail + need) if (avail + need) > 0 else 0.0
    to_coat = 4.0 * coat_f * y[IX["Mfree"]]
    d[IX["Mfree"]] = -to_coat - p["k_sen_free"] * y[IX["Mfree"]]
    d[IX["Mcoat"]] = to_coat - p["k_clear_coat"] * y[IX["Mcoat"]]
    # PRIMING IS PROPORTIONAL TO THE EXPOSURE INTEGRAL, not to a clearance
    # flux: what sensitises a woman is antigen sitting in her circulation
    # available to antigen-presenting cells.  Coated cells are cleared by the
    # splenic route within hours and (with efficiency dev) do not prime.
    ag_prime = y[IX["Mfree"]] + (1.0 - p["dev"]) * y[IX["Mcoat"]]
    d[IX["Sens"]] = p["ksens"] * ag_prime
    d[IX["Mbm"]] = p["k_prime"] * ag_prime * (1.0 - y[IX["Mbm"]] / 5.0) \
        - p["k_bm"] * y[IX["Mbm"]]
    # Long-lived marrow plasma cells hold an alloantibody titre for years; a
    # fetomaternal haemorrhage (including the one every IUT causes) boosts the
    # pool, which then decays BACK TO ITS SET POINT rather than to zero.  A
    # first-order decay to zero would make every sensitised mother's titre
    # collapse during the pregnancy, which is not what is observed.
    mset = getattr(c, "mpc_set", 0.0)
    d[IX["Mpc"]] = p["k_boost"] * y[IX["Mbm"]] * ag_prime \
        - p["k_pc"] * (y[IX["Mpc"]] - mset)

    # ---------------- nipocalimab PK ----------------------------------------
    cl_nip = p["nip_CLlin"] + p["nip_CLfcrn"] * p["K_nip"] / (p["K_nip"] + r["Cnip"])
    d[IX["Nc"]] = -cl_nip * r["Cnip"] - p["nip_Q"] * (r["Cnip"] - y[IX["Np"]] / p["nip_V2"])
    d[IX["Np"]] = p["nip_Q"] * (r["Cnip"] - y[IX["Np"]] / p["nip_V2"])

    # ---------------- the placental conveyor ---------------------------------
    if not born:
        cap = p["v0"] * math.exp(p["g_pl"] * (ga - 20.0))
        fu_pl = 1.0 / (1.0 + (r["Cig"] + p["ivig_compete"] * r["Civ"]) / p["K_igg"]
                       + p["nip_pl_pen"] * r["Cnip"] / p["K_nip"])
        J_ig = cap * fu_pl * r["Cig"]                          # g/d
        J_iv = cap * fu_pl * r["Civ"]
        J_a1 = cap * fu_pl * r["Ca1"] * 1000.0                 # IU/d
        J_a3 = p["sub3"] * cap * fu_pl * r["Ca3"] * 1000.0
        J_rh = cap * fu_pl * r["Crh"] * 1000.0
        J_nip = 0.55 * cap * fu_pl * r["Cnip"]                 # mg/d
    else:
        fu_pl = 0.0
        J_ig = J_iv = J_a1 = J_a3 = J_rh = J_nip = 0.0
    # mass balance: the conveyor is a SINK on the maternal side.  Leaving this
    # out overstates maternal titres late in gestation, when the fetal IgG pool
    # has grown to ~20% of the maternal one.
    d[IX["Mig"]] -= J_ig
    d[IX["Miv"]] -= J_iv
    d[IX["Ma1"]] -= J_a1
    d[IX["Ma3"]] -= J_a3
    d[IX["Mrhig"]] -= J_rh
    d[IX["Nc"]] -= J_nip
    kel_f = p["k_int_f"] * (1.0 - p["phi_rescue"] /
                            (1.0 + r["Cf_ig"] / p["K_igg"] + r["Cf_nip"] / p["K_nip"]))
    d[IX["Fig"]] = J_ig + J_iv - kel_f * y[IX["Fig"]]
    d[IX["Fa1"]] = J_a1 + J_rh - kel_f * y[IX["Fa1"]]
    d[IX["Fa3"]] = J_a3 - kel_f * y[IX["Fa3"]]
    d[IX["Fnip"]] = J_nip - (p["nip_CLlin"] / p["nip_V1"] * 1.5) * y[IX["Fnip"]]
    d[IX["Adose"]] = J_a1 + J_a3

    # ---------------- growth -------------------------------------------------
    if not born:
        target = efw_ref(ga)
        # growth is slowed by severe anaemia/hydrops
        d[IX["EFW"]] = (target - y[IX["EFW"]]) * 0.35 + \
            (efw_ref(ga + 0.01) - target) / 0.07 * (1.0 - 0.45 * y[IX["Card"]])
        d[IX["Wgt"]] = 0.0
    else:
        d[IX["EFW"]] = 0.0
        d[IX["Wgt"]] = 0.018 * y[IX["Wgt"]] * (1.0 - 0.3 * y[IX["Card"]])
    # plasma volume is REGULATED, and it expands as the haematocrit falls --
    # which is why a severely anaemic fetus is not simply a smaller circulation
    hct_ref = r["hbr"] / p["mchc"]
    vb_ref = p["bv_neo"] * r["wt"] if born else p["kbv"] * r["efw"]
    vpl_t = vb_ref * (1.0 - hct_ref) * \
        (1.0 + p["vpl_expand"] * max(0.0, 1.0 - r["hct"] / hct_ref))
    d[IX["Vpl"]] = p["k_vol"] * (vpl_t - y[IX["Vpl"]])

    # ---------------- THE PRODUCT: red-cell destruction ----------------------
    ops = p["kops"] * r["A_eff"] * (y[IX["Ro"]] + y[IX["Rr"]])       # g Hb/d
    clr = p["vmax_res"] * r["M_eff"] * y[IX["Rop"]] / (p["Kres"] + y[IX["Rop"]] + 1e-12)
    sen_own = math.log(2.0) / p["t_own"] * y[IX["Ro"]]
    sen_don = math.log(2.0) / p["t_don"] * y[IX["Rd"]]
    # PRODUCTION IS WRITTEN AS A MULTIPLE OF THE REQUIREMENT, not as an
    # absolute rate: at Prog = 1 the fetus exactly replaces senescence and
    # exactly keeps up with the growth of its own blood volume.  This is why
    # k_prod = 1 is not a fitted number.
    if born:
        mref = hb_ref_pn(c.ga_birth, r["pna"]) * p["bv_neo"] * r["wt"] / 100.0
        dmref = (hb_ref_pn(c.ga_birth, r["pna"] + 0.05) * p["bv_neo"] *
                 r["wt"] / 100.0 - mref) / 0.05
    else:
        mref = hb_ref(ga) * p["kbv"] * efw_ref(ga) / 100.0
        dmref = (hb_ref(ga + 0.01) * p["kbv"] * efw_ref(ga + 0.01) / 100.0
                 - mref) / 0.07
    req = math.log(2.0) / p["t_own"] * mref + max(dmref, 0.0)
    prod = p["k_prod"] * (y[IX["Prog"]] + p["emh_prod"] * y[IX["EMH"]]) * req
    d[IX["Rr"]] = prod - p["k_mat"] * y[IX["Rr"]] - ops * y[IX["Rr"]] / \
        max(y[IX["Ro"]] + y[IX["Rr"]], 1e-9)
    d[IX["Ro"]] = p["k_mat"] * y[IX["Rr"]] - sen_own - ops * y[IX["Ro"]] / \
        max(y[IX["Ro"]] + y[IX["Rr"]], 1e-9)
    d[IX["Rd"]] = -sen_don
    d[IX["Rop"]] = ops - clr
    d[IX["Dest"]] = clr
    d[IX["SenC"]] = sen_own + sen_don
    d[IX["Res"]] = p["k_res_grow"] * (y[IX["Rop"]] / (p["Kres"] + y[IX["Rop"]] + 1e-12)) * \
        (p["res_max"] - y[IX["Res"]]) - 0.05 * (y[IX["Res"]] - 1.0)

    # ---------------- erythropoiesis ----------------------------------------
    epo_t = p["epo0"] * min((1.0 / max(r["hbmom"], 0.20)) ** p["epo_n"], p["epo_cap"])
    d[IX["Epo"]] = p["k_epo"] * (epo_t - y[IX["Epo"]])
    epo_drive = (y[IX["Epo"]] / p["epo0"]) ** p["epo_pow"]
    kill_prog = p["kell_kill"] * r["A_eff"] * y[IX["Prog"]]
    d[IX["Prog"]] = p["k_prog"] * epo_drive * (1.0 - y[IX["Prog"]] / p["prog_max"]) \
        - p["k_prog_d"] * y[IX["Prog"]] - kill_prog
    d[IX["EMH"]] = p["k_emh"] * max(0.0, y[IX["Prog"]] - p["emh_thresh"]) \
        - p["k_emh_off"] * y[IX["EMH"]]
    d[IX["Nrbc"]] = 0.6 * (y[IX["Prog"]] - y[IX["Nrbc"]])

    # ---------------- liver, albumin, hydrops -------------------------------
    alb_target = alb_ref(ga if not born else c.ga_birth, p) * \
        (1.0 - p["emh_alb"] * y[IX["EMH"]] / (1.0 + y[IX["EMH"]]))
    d[IX["Alb"]] = p["k_alb_cat"] * (alb_target - y[IX["Alb"]])
    if not born:
        d[IX["Asc"]] = max(r["jv"] - r["lymph"], -0.25 * y[IX["Asc"]])
    else:
        d[IX["Asc"]] = -0.25 * y[IX["Asc"]]
    d[IX["Card"]] = p["k_fail"] * max(0.0, r["demand_brain"] / p["co_max"] - 1.0) * \
        (1.0 - y[IX["Card"]]) - p["k_recov"] * y[IX["Card"]]
    d[IX["Lac"]] = p["k_lac"] * r["hypox"] - p["k_lac_out"] * y[IX["Lac"]]

    # ---------------- bilirubin ---------------------------------------------
    bil_prod = p["bil_per_g"] * (clr + sen_own + sen_don) + p["bil_shunt"] * r["wt"]
    if not born:
        bil_out = p["k_bil_pl"] * y[IX["Bil"]]
        amn_in = p["k_amn"] * y[IX["Bil"]] * 0.06
        conj = 0.0
        photo = 0.0
        ehc = 0.0
    else:
        conj = p["vmax_ugt"] * r["wt"] * y[IX["Ugt"]] * r["tsb"] / (r["tsb"] + p["km_bil"])
        bil_out = conj
        amn_in = 0.0
        photo = p["k_photo"] * c.photo * y[IX["Bil"]]
        ehc = p["k_ehc"] * y[IX["Bgut"]] * c.ehc_on
    # plasma <-> extravascular exchange: the reason phototherapy rebounds and
    # the reason an exchange transfusion removes less bilirubin than the
    # plasma pool suggests.
    xchg = p["k_bx"] * (r["tsb"] - r["ctis"]) * p["vd_bp"] * r["wt"] * 10.0
    d[IX["Bil"]] = bil_prod - bil_out - photo + p["k_iso_back"] * y[IX["Biso"]] \
        + ehc - xchg
    d[IX["Bex"]] = xchg
    d[IX["Biso"]] = photo - p["k_iso_out"] * y[IX["Biso"]] - p["k_iso_back"] * y[IX["Biso"]]
    d[IX["Bgut"]] = conj - p["k_gut"] * y[IX["Bgut"]] - ehc
    d[IX["Amn"]] = amn_in - p["k_amn_out"] * y[IX["Amn"]]
    d[IX["Ugt"]] = 3.0 * (ugt_ref(r["pna"], p) - y[IX["Ugt"]]) if born else 0.0
    d[IX["Bind"]] = max(0.0, r["ba_ratio"] - p["bind_thresh"])

    # ---------------- hazard -------------------------------------------------
    ph = p["ph0"] - p["ph_lac"] * y[IX["Lac"]]
    d[IX["Hzd"]] = p["h0"] + (p["h_hydrops"] if r["hydrops"] else 0.0) \
        + (p["h_severe"] if r["hbmom"] < 0.35 else 0.0) \
        + (p["h_acid"] if ph < 7.15 else 0.0)
    d[IX["Fer"]] = 0.0
    return d


# ==============================================================================
# 4.  INITIAL CONDITIONS
# ==============================================================================
def init_state(p, ga0=12.0, anti_d_iu=0.0, hb_frac=1.0):
    y = np.zeros(NS)
    Vm = vm_plasma(ga0)
    y[IX["Mig"]] = p["igg0"] * Vm
    y[IX["Ma1"]] = anti_d_iu * Vm * 1000.0 * 0.85
    y[IX["Ma3"]] = anti_d_iu * Vm * 1000.0 * 0.15
    # plasma-cell pool that sustains this titre at steady state
    kel = p["k_int"] * (1.0 - p["phi_rescue"] / (1.0 + p["igg0"] / p["K_igg"]))
    y[IX["Mpc"]] = kel * y[IX["Ma1"]] / p["sec_pc"] if anti_d_iu > 0 else 0.0
    y[IX["Mbm"]] = 2.0 if anti_d_iu > 0 else 0.0
    y[IX["EFW"]] = efw_ref(ga0)
    efw = y[IX["EFW"]]
    Vfp = p["kbv"] * efw
    hb = hb_ref(ga0) * hb_frac
    mass = hb * Vfp / 100.0
    y[IX["Ro"]] = mass * 0.92
    y[IX["Rr"]] = mass * 0.08
    y[IX["Vpl"]] = Vfp - mass / (p["mchc"] / 100.0)
    y[IX["Prog"]] = 1.0
    y[IX["Epo"]] = p["epo0"]
    y[IX["Res"]] = 1.0
    hct = hb * p["hct_hb"] / 100.0
    Vpl = Vfp * (1.0 - hct)
    y[IX["Alb"]] = alb_ref(ga0, p)
    y[IX["Ugt"]] = p["ugt_birth"]
    y[IX["Bil"]] = 0.4 * p["vd_bp"] * (y[IX["EFW"]] / 1000.0) * 10.0
    y[IX["Bex"]] = 0.4 * p["vd_bt"] * (y[IX["EFW"]] / 1000.0) * 10.0
    y[IX["Nrbc"]] = 1.0
    y[IX["Wgt"]] = efw / 1000.0
    y[IX["Fig"]] = 0.0
    return y


# ==============================================================================
# 5.  IUT ARITHMETIC
# ==============================================================================
def iut_volume(p, y, c, target_hct=0.45, donor_hct=0.80, gamax=None):
    """Volume (mL) of donor red cells needed to reach target haematocrit.

    Mass balance, not a nomogram:  (Vfp*Hct0 + V*Hct_don) / (Vfp + V) = Hct_t
    """
    r = readouts(0.0, y, p, c)
    if r["hct"] >= target_hct:
        return 0.0
    # (Vrbc + Hct_don*V) / (Vb + V) = Hct_target
    V = (target_hct * r["Vfp"] - r["Vrbc"]) / (donor_hct - target_hct)
    Vfp = r["Vfp"]
    # a single procedure does not expand the fetoplacental circulation by more
    # than ~60% (Nishie 2012 reported a mean expansion of 51%)
    return float(min(V, 0.60 * Vfp))


def do_iut(y, p, c, target_hct=0.45, donor_hct=0.80):
    V = iut_volume(p, y, c, target_hct, donor_hct)
    if V <= 0:
        return 0.0
    y[IX["Rd"]] += donor_hct * V * p["mchc"] / 100.0
    y[IX["Vpl"]] += (1.0 - donor_hct) * V
    y[IX["Fer"]] += V * 1.08          # ~1.08 mg Fe per mL of packed red cells
    # a transfusion is also a fetomaternal haemorrhage: it boosts the mother
    y[IX["Mfree"]] += 0.35
    c.iut_n += 1
    y[IX["Hzd"]] += p["h_iut"]
    return V


# ==============================================================================
# 6.  SIMULATION DRIVER
# ==============================================================================
def simulate(p, c, ga0=12.0, ga_end=None, anti_d_iu=6.0, hb_frac=1.0,
             protocol="mca", ga_deliver=37.0, postnatal_days=90.0,
             iut_target=0.45, mca_cut=1.50, record_every=0.5, y0=None,
             fmh_events=(), rhig_events=(), seed=None):
    """Integrate from ga0 (weeks) through delivery and the neonatal period.

    protocol: "none"  -- observe only
              "mca"   -- weekly MCA-PSV, cordocentesis+IUT if >= mca_cut MoM
              "fixed" -- IUT every 14 d once the first is triggered
    """
    rng = np.random.default_rng(seed if seed is not None else 12345)
    y = init_state(p, ga0, anti_d_iu, hb_frac) if y0 is None else y0.copy()
    c.mpc_set = float(y[IX["Mpc"]])
    t = ga0 * 7.0
    t_birth = ga_deliver * 7.0
    t_end = t_birth + postnatal_days
    rows = []
    next_iut = None
    fmh = {round(g, 3): v for g, v in fmh_events}
    rhig = {round(g, 3): v for g, v in rhig_events}
    step = 1.0
    while t < t_end - 1e-9:
        ga = t / 7.0
        # ---- discrete weekly / daily actions --------------------------------
        for g, v in list(fmh.items()):
            if ga >= g:
                y[IX["Mfree"]] += v
                fmh.pop(g)
        for g, v in list(rhig.items()):
            if ga >= g:
                y[IX["Mdep"]] += v * p["ad_iu_per_ug"]
                rhig.pop(g)
        if (not c.born) and c.nip_dose > 0 and c.nip_start <= ga <= c.nip_stop:
            if abs((ga - c.nip_start) * 7.0 % 7.0) < 0.5:
                y[IX["Nc"]] += c.nip_dose * c.mat_wt
        if (not c.born) and c.ivig_dose > 0 and ga >= c.ivig_start:
            if abs((ga - c.ivig_start) * 7.0 % 7.0) < 0.5:
                y[IX["Miv"]] += c.ivig_dose * c.mat_wt * p["ivig_f"]
        # ---- obstetric surveillance and IUT ---------------------------------
        if (not c.born) and protocol != "none" and 17.5 <= ga <= 35.0:
            r = readouts(t, y, p, c)
            trigger = False
            weekly = (ga - ga0) * 7.0 % 7.0 < 0.5
            if protocol == "mca" and weekly:
                psv_obs = r["psv_mom"] * (1.0 + 0.10 * rng.standard_normal())
                trigger = psv_obs >= mca_cut
            elif protocol == "fixed" and weekly and next_iut is None:
                # "fixed" means: once anaemia is established, transfuse on the
                # calendar rather than on the Doppler.  Something still has to
                # start the course, and that something is the first cordocentesis
                # showing Hb below 0.85 MoM.
                trigger = r["hb"] < 0.85 * r["hbr"]
            if next_iut is not None and ga >= next_iut:
                trigger = True
            scheduled = next_iut is not None and ga >= next_iut
            if trigger:
                hb_obs = r["hb"] * (1.0 + 0.03 * rng.standard_normal())
                # a FIRST transfusion needs the cordocentesis to confirm
                # anaemia; a SCHEDULED repeat is given on the calendar, which
                # is what fetal-therapy centres actually do, because MCA-PSV
                # loses accuracy once donor cells are in the circulation
                if scheduled or protocol == "fixed" or hb_obs < 0.75 * r["hbr"]:
                    do_iut(y, p, c, target_hct=iut_target)
                    r2 = readouts(t, y, p, c)
                    c.iut_log.append((ga, r["hb"], r2["hb"], c.iut_n))
                    next_iut = ga + (2.0 if (protocol == "fixed" or c.iut_n == 1)
                                     else 3.0)
                    if next_iut > 34.5:
                        next_iut = None
        # ---- birth -----------------------------------------------------------
        if (not c.born) and t >= t_birth - 1e-9:
            c.born = True
            c.ga_birth = ga_deliver
            c.ehc_on = 1.0
            r = readouts(t, y, p, c)
            y[IX["Wgt"]] = r["efw"] / 1000.0
        # ---- neonatal management --------------------------------------------
        if c.born:
            r = readouts(t, y, p, c)
            hours = (t - t_birth) * 24.0
            thr_photo = 5.0 + min(hours, 96.0) * 0.10
            thr_exch = 12.0 + min(hours, 96.0) * 0.09
            c.photo = 30.0 if r["tsb"] > thr_photo else 0.0
            last_ex = getattr(c, "t_exch", -99.0)
            if r["tsb"] > thr_exch and (t - last_ex) > 0.5:
                # A double-volume exchange is VOLUME-NEUTRAL: it removes about
                # 85% of everything in the circulation -- bilirubin, antibody,
                # the infant's own cells AND any donor cells already there --
                # and replaces the red cell volume with donor cells.  Adding
                # donor cells without removing the old ones (the first version
                # of this block) makes the donor pool, its senescence, and
                # therefore bilirubin production grow with every procedure: a
                # positive feedback loop that ran to TSB 1.2e7 mg/dL.
                fx = 0.85
                y[IX["Bil"]] *= (1.0 - fx)
                y[IX["Bex"]] *= (1.0 - 0.35)      # tissue pool re-equilibrates
                y[IX["Fa1"]] *= (1.0 - fx)
                y[IX["Fa3"]] *= (1.0 - fx)
                y[IX["Fig"]] *= (1.0 - fx)
                keep = 1.0 - fx
                y[IX["Ro"]] *= keep
                y[IX["Rr"]] *= keep
                y[IX["Rd"]] *= keep
                # restore the red cell volume to a post-exchange Hct of 0.50
                mass_now = y[IX["Ro"]] + y[IX["Rd"]] + y[IX["Rr"]]
                vpl = max(y[IX["Vpl"]], 1.0)
                mass_t = 0.50 / (1.0 - 0.50) * vpl * p["mchc"] / 100.0
                if mass_t > mass_now:
                    y[IX["Rd"]] += mass_t - mass_now
                    y[IX["Fer"]] += (mass_t - mass_now) / (p["mchc"] / 100.0) * 1.08
                c.__dict__.setdefault("exch_n", 0)
                c.exch_n += 1
                c.t_exch = t
            if r["hb"] < 7.5 and hours > 24.0:
                y[IX["Rd"]] += 0.60 * 0.15 * r["Vfp"] * p["mchc"] / 100.0
                y[IX["Vpl"]] += 0.40 * 0.15 * r["Vfp"]
                c.__dict__.setdefault("topup_n", 0)
                c.topup_n += 1
        # ---- integrate -------------------------------------------------------
        h = min(step, t_end - t)
        sol = solve_ivp(rhs, (t, t + h), y, args=(p, c), method="LSODA",
                        rtol=1e-6, atol=1e-9, max_step=h)
        if not sol.success:
            raise RuntimeError("integration failed at GA %.2f: %s" % (ga, sol.message))
        y = np.maximum(sol.y[:, -1], 0.0)
        t += h
        r = readouts(t, y, p, c)
        rows.append(dict(t=t, ga=t / 7.0, **{k: r[k] for k in
                    ("hb", "hbmom", "psv_mom", "tsb", "od450", "Cf_a1", "Ca1",
                     "Cig", "Civ", "Cnip", "Cf_nip", "f_ag", "A_eff", "alb",
                     "hct", "efw", "hydrops", "cvp", "do2", "Cf_ig", "M_eff",
                     "ba_ratio", "hbr", "pna", "demand_brain", "wt", "Vfp",
                     "Vpl", "Vrbc", "jv", "lymph", "Cf_a3", "asc_thr",
                     "cvp", "M_eff", "Crh", "fu_m", "bfree_ratio",
                     "do2_abs")},
                    Asc=y[IX["Asc"]], Ro=y[IX["Ro"]], Rd=y[IX["Rd"]],
                    Rr=y[IX["Rr"]], Rop=y[IX["Rop"]], Prog=y[IX["Prog"]],
                    EMH=y[IX["EMH"]], Epo=y[IX["Epo"]], Res=y[IX["Res"]],
                    Sens=y[IX["Sens"]], Hzd=y[IX["Hzd"]], Lac=y[IX["Lac"]],
                    Amn=y[IX["Amn"]], Dest=y[IX["Dest"]], Bind=y[IX["Bind"]],
                    Fer=y[IX["Fer"]], iut_n=c.iut_n, Nrbc=y[IX["Nrbc"]],
                    Card=y[IX["Card"]], Ugt=y[IX["Ugt"]], Bil=y[IX["Bil"]],
                    Mfree=y[IX["Mfree"]], Mcoat=y[IX["Mcoat"]], Mpc=y[IX["Mpc"]],
                    Fig=y[IX["Fig"]], born=1.0 if c.born else 0.0))
    out = {k: np.array([row[k] for row in rows], float) for k in rows[0]}
    out["_ctl"] = c
    out["_y"] = y
    return out


def at_ga(res, key, ga):
    return float(np.interp(ga, res["ga"], res[key]))
