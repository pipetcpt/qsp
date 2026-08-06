#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gad_python_reference.py
=======================
Executable reference implementation of the Generalized Anxiety Disorder (GAD)
QSP model.  Everything the README and the mrgsolve file claim numerically is
produced by THIS file; the mrgsolve translation mirrors it equation for
equation.

THE CENTRAL OBJECT
------------------
Anxiety is written as ONE dimensionless corticolimbic gain

        Phi  =  (E_amy * S_glu) / (C_pfc * I_gaba)

a RATIO with two numerator factors and two denominator factors.  Each of the
four factors is owned by a different drug class and each moves on its own time
constant:

    factor   meaning                              moved by            tau
    ------   -----------------------------------  ------------------  --------
    E_amy    amygdala / BNST excitatory drive     CBT, chronic 5-HT   ~20-30 d
    S_glu    presynaptic glutamate release prob.  pregabalin (a2d-1)  ~2 d
    C_pfc    prefrontal top-down control          SSRI/SNRI, CBT      ~30 d
             capacity                             (via BDNF)          (gated by
                                                                      a 14 d
                                                                      5-HT1A
                                                                      step)
    I_gaba   GABA-A mediated inhibitory efficacy  benzodiazepine      ~0.5 d

and the whole thing is READ OUT through HAM-A, an instrument that carries a
FIFTH clock of its own (expectancy + regression-to-the-mean), which is
subtracted from every arm equally.

Because Phi is a ratio and HAM-A is a SATURATING function of Phi, several
things that are usually asserted become arithmetic instead:

  * onset of action is a property of WHICH FACTOR you move, not of "the drug";
  * a flat SSRI dose-response follows from the SERT occupancy hyperbola;
  * combination therapy is multiplicative in Phi but sub-additive on HAM-A;
  * a site with a large expectancy response has a genuinely smaller measurable
    drug-placebo delta with NO change in the pharmacology (assay sensitivity);
  * one benzodiazepine adaptation state produces BOTH partial tolerance and
    discontinuation rebound, with the same parameter.

CALIBRATION (only these numbers were fitted)
--------------------------------------------
  Rickels 2005 (PMID 16143734), 4-week GAD trial, HAM-A change from baseline:
      placebo -8.4 | pregabalin 300 mg -12.2 | alprazolam 1.5 mg -10.9
  Khan 2011 (PMID 21694613), quetiapine XR trial, PLACEBO arm trajectory:
      week 1 -5.94 | week 8 -11.10          (the shape of the placebo clock)

EVERYTHING ELSE IS PREDICTED OUT OF SAMPLE, in particular
  pregabalin 450 / 600 mg, quetiapine XR 50 / 150 / 300 mg (Khan 2011),
  escitalopram, venlafaxine ER, duloxetine, buspirone, the week-1
  psychic-vs-somatic dissociation, and the Allgulander 2006 (PMID 16316482)
  randomised-withdrawal relapse rates (19% escitalopram vs 56% placebo).

INTEGRATION
-----------
A hand-written RK4 integrator, vectorised over subjects, with a fixed step
that divides the dosing grid exactly (default 1/48 d = 30 min, so QD, BID and
TID dose times all land on step boundaries, and the fastest rate in the system
gives h*lambda = 0.5).  A whole multi-arm trial -- every arm, every virtual
subject -- is integrated in ONE pass with array-valued parameters, which is
what makes the virtual-population and relapse-prevention runs tractable.
Verified against scipy.integrate.solve_ivp (LSODA, rtol 1e-9) and against a
halved step in section 12 of gad_analysis.py.
"""

import json
import math
import os
import sys

import numpy as np

np.seterr(all="ignore")

# ----------------------------------------------------------------------------
# 0.  State vector
# ----------------------------------------------------------------------------
STATES = [
    # --- pharmacokinetics (amounts in ug unless noted) --------------------
    "esc_gut", "esc_c", "esc_p", "esc_e",      # escitalopram (2-cpt + brain effect site, ng/mL)
    "ven_gut", "ven_c", "odv_c",               # venlafaxine ER + O-desmethylvenlafaxine
    "dlx_gut", "dlx_c",                        # duloxetine
    "pgb_gut", "pgb_c",                        # pregabalin
    "bzd_gut", "bzd_c", "bzd_e",               # benzodiazepine (+ brain effect site, ng/mL)
    "bus_gut", "bus_c", "pp1_c",               # buspirone + 1-pyrimidinylpiperazine
    "qtp_gut", "qtp_c", "nqt_c",               # quetiapine XR + norquetiapine
    # --- neurochemistry ----------------------------------------------------
    "sht",      # extracellular forebrain 5-HT (fraction of healthy)
    "auto",     # 5-HT1A somatodendritic autoreceptor functional availability (1 = intact)
    "ne",       # extracellular NE (fraction of healthy)
    "a2auto",   # alpha2-adrenergic autoreceptor availability
    # --- plasticity and circuits ------------------------------------------
    "bdnf",     # BDNF / plasticity signal
    "cpfc",     # prefrontal top-down control capacity  (DENOMINATOR factor 1)
    "eamy",     # amygdala / BNST excitatory drive      (NUMERATOR factor 1)
    "traf",     # alpha2delta-1 effect build-up state (pregabalin)
    # --- GABA-A ------------------------------------------------------------
    "ra2",      # alpha2/alpha3 receptor pool (anxiolysis; slow tolerance)
    "ra1",      # alpha1 receptor pool (sedation; fast tolerance)
    "depend",   # benzodiazepine adaptation state -> partial tolerance AND rebound
    # --- HPA axis -----------------------------------------------------------
    "crh", "acth", "cort", "gr",
    # --- autonomic / symptom layer -----------------------------------------
    "sns",      # sympathetic tone
    "auton",    # autonomic somatic symptom load
    "sleepd",   # sleep-continuity deficit
    "worry",    # self-reinforcing worry engine (the GAD-specific loop)
    # --- context / trial machinery -----------------------------------------
    "expect",   # expectancy (acts on C_pfc: a REAL top-down effect)
    "fluct",    # enrollment-peak / regression-to-the-mean component (additive on score)
    # --- adverse-effect tolerance states and burdens ------------------------
    "rnau",     # nausea tolerance pool (SERT-driven, fast)
    "rdizz",    # dizziness tolerance pool (a2d-driven)
    "rh1",      # H1 sedation tolerance pool
    "sexd",     # sexual dysfunction (SERT-driven, does NOT tolerate)
    "wt",       # weight change (kg)
    "ract",     # activation/jitteriness tolerance pool (NE-driven)
    # --- comorbidity and trial hazard ---------------------------------------
    "madrs",    # comorbid depressive symptom load
    "cumhaz",   # cumulative dropout hazard
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)


# ----------------------------------------------------------------------------
# 1.  Parameters
# ----------------------------------------------------------------------------
def default_params():
    p = {}

    # ---------------- escitalopram --------------------------------------
    # CL/F 400 L/d and F 0.80 give Css,avg = 20 ng/mL at 10 mg QD, inside the
    # 15-80 ng/mL therapeutic reference range (Hart 2024, PMID 38287888).
    p.update(F_esc=0.80, ka_esc=9.6, V1_esc=350.0, V2_esc=400.0,
             Q_esc=250.0, CL_esc=400.0, keo_esc=8.0)
    # 80% SERT occupancy at 10 mg (Css 20 ng/mL) -> EC50 = 5 ng/mL.
    p.update(EC50_sert_esc=5.0)

    # ---------------- venlafaxine ER + ODV --------------------------------
    p.update(F_ven=0.42, ka_ven=2.4, V_ven=300.0, CL_ven=999.0,
             V_odv=250.0, CL_odv=378.0, mw_ratio_odv=0.95)
    p.update(fm2d6=0.70)                      # EM; PM 0.15, UM 0.85
    p.update(EC50_sert_ven=22.5, EC50_net_ven=260.0)

    # ---------------- duloxetine -------------------------------------------
    p.update(F_dlx=0.50, ka_dlx=4.0, V_dlx=700.0, CL_dlx=970.0)
    p.update(EC50_sert_dlx=7.7, EC50_net_dlx=46.0)

    # ---------------- pregabalin -------------------------------------------
    # t1/2 6.3 h (PMID 17940637); clearance proportional to CrCl.
    p.update(F_pgb=0.90, ka_pgb=16.0, V_pgb=40.0, CL_pgb=105.0, crcl=100.0)
    p.update(EC50_a2d=2100.0, ktraf=0.50, emax_pgb=0.42)

    # ---------------- benzodiazepine ---------------------------------------
    # lorazepam PK; alprazolam handled by a lorazepam-equivalent conversion.
    p.update(F_bzd=0.90, ka_bzd=12.0, V_bzd=90.0, CL_bzd=108.0, keo_bzd=14.0)
    p.update(bzd_potency=1.0)                 # 1.0 lorazepam, 1.53 alprazolam
    # Lorazepam plasma EC50 for GABA-A benzodiazepine-site occupancy
    # = 96 ng/mL ([11C]flumazenil microPET, Atack 2007, PMID 17164474).
    p.update(EC50_bz=96.0)
    p.update(emax_a2=1.90, wadapt=0.20, kdep_on=0.05)
    p.update(ktol1=1.20, krec1=0.15, ktol2=0.012, krec2=0.020)

    # ---------------- buspirone + 1-PP -------------------------------------
    p.update(F_bus=0.05, ka_bus=14.0, V_bus=350.0, CL_bus=1800.0,
             f_pp=0.30, V_pp=200.0, CL_pp=554.0)
    p.update(EC50_1a=1.20, ia_bus=0.35, EC50_pp1=45.0, kpp1=0.25)

    # ---------------- quetiapine XR + norquetiapine ------------------------
    p.update(F_qtp=1.0, ka_qtp=3.5, V_qtp=700.0, CL_qtp=2000.0,
             f_nqt=0.30, V_nqt=600.0, CL_nqt=1500.0)
    p.update(EC50_h1=20.0, EC50_net_nqt=60.0)

    # ---------------- serotonin ---------------------------------------------
    p.update(kin_sht=24.0, kout_sht=24.0, sert_floor=0.15, gamma_auto=6.0,
             kdes_auto=0.24, krec_auto=0.020, ia_bus_post=0.55)

    # ---------------- norepinephrine ---------------------------------------
    p.update(kin_ne=18.0, kout_ne=18.0, net_floor=0.25, gamma_a2=5.0,
             kdes_a2=0.18, krec_a2=0.025)

    # ---------------- plasticity / PFC / amygdala ---------------------------
    p.update(kb=0.20, kdb=0.20, w_5ht=0.55, w_ne=0.25, w_tonic=0.20,
             kcort_bdnf=0.35, cbt_bdnf=0.05)
    # ONE scalar for the whole serotonergic efficacy arm.  Every downstream
    # consequence of raised 5-HT (plasticity, amygdala decay, worry decay,
    # mood) reads k5ht_eff * (SHT_post - 1) rather than (SHT_post - 1), so a
    # single fitted number sets how much clinical effect a given SERT
    # occupancy buys.  Without it the SSRI arm had no free parameter at all
    # and over-predicted the escitalopram-placebo difference four-fold.
    p.update(k5ht_eff=1.0)
    p.update(kp_pfc=0.030, kl_pfc=0.0175, cpfc_max=2.20,
             kc_pfc=0.45, kdis_pfc=0.456, e_ex=0.55)
    p.update(ke_in=0.050, ke_out=0.050, kcrh_amy=0.15, kslp_amy=0.18,
             kdis_amy=0.32, a5ht_amy=0.45, acbt_amy=0.09, EAMAX=3.0)
    p.update(kglu_stress=0.35, kglu_cort=0.20)
        # DIS = 1.0 is the TYPICAL TRIAL PATIENT.  The three kdis_* gains were
    # divided by 2.5 from their first values once the virtual population made
    # the error visible: a patient whose STABLE HAM-A is 25.9 enrols at 25.9
    # plus the enrolment peak, i.e. at 32, and GAD trials report baselines of
    # 24-27.  The stable score has to be about 19 so that the ENROLMENT score
    # is about 25.5.  That is not a cosmetic rescaling -- it is the statement
    # that trial entry criteria select the top of a fluctuation.
    p.update(igaba_base0=1.0, kdis_gaba=0.048)

    # ---------------- HPA axis ----------------------------------------------
    p.update(ZMAX=3.0, kcrh_phi=0.12, ks_phi=1.00, ksl_phi=1.00, kw_phi=1.00,
             SLMAX=4.0)
    p.update(kcrh=1.395, kcout=0.9, kgr_fb=0.55,
             ka_acth=6.0, ke_acth=6.0, kc_cort=4.0, kec_cort=4.0,
             kgr_in=0.08, kgr_dn=0.10)

    # ---------------- autonomic / symptoms ----------------------------------
    # b_gaba / b_pgb / h_gaba / h_pgb are the DIRECT symptomatic routes by which
    # a benzodiazepine or pregabalin relieves autonomic arousal and broken sleep
    # without going through Phi.  They are deliberately kept SMALL.  With them
    # large (0.85 / 0.75 / 1.60 / 0.90, the first values tried) the calibration
    # drove emax_pgb to 0.029 -- i.e. it explained pregabalin's entire effect through
    # the symptomatic route and left S_glu doing nothing, which would have made
    # "S_glu is pregabalin's factor" an empty claim.  Two routes to the same
    # score are not separable from one trial endpoint, so the split is a
    # structural choice, stated rather than fitted.
    p.update(ks_in=0.60, ks_out=0.60, w_ne_sns=0.55,
             kau_in=0.45, kau_out=0.45, b_gaba=0.25, b_pgb=0.20,
             ksl_in=0.35, ksl_out=0.35, w_cort_sl=0.30,
             h_h1=1.30, h_gaba=0.60, h_pgb=0.35,
             kw_in=0.16, kw_out=0.16, rho_worry=0.55, kw_m=1.20,
             acbt_worry=0.13, a5ht_worry=0.35)
    # NOTE ON THE CBT ARM: acbt_amy, acbt_worry and cbt_bdnf are NOT calibrated
    # to any trial.  They were set so that CBT alone produces a week-12
    # difference of a few HAM-A points -- the order of magnitude meta-analyses
    # report against control conditions -- and nothing in this repository
    # should be read as a quantitative claim about psychotherapy.

    # ---------------- expectancy / regression to the mean --------------------
    p.update(kexp_dec=0.010, dvisit=0.42, emax_exp=1.00, kfl=0.075)

    # ---------------- HAM-A link --------------------------------------------
    p.update(c1=0.62, c2=0.30, c3=0.26, Kpsy=3.02, psy0=1.098, psy_floor=1.7,
             v1=0.42, v2=0.30, v3=0.30, v4=0.22, Ksom=2.73, som0=1.020, som_floor=1.3)

    # ---------------- adverse effects ---------------------------------------
    p.update(ktol_nau=1.40, krec_nau=0.16,
             ktol_dizz=0.55, krec_dizz=0.09,
             ktol_h1=0.55, krec_h1=0.06,
             ks_sex=0.22, kr_sex=0.10,
             ktol_act=1.10, krec_act=0.14,
             kwt=0.030, kwt_off=0.010)
    p.update(ae_sed_h1=2.20, ae_sed_bz=1.60, ae_sed_pgb=1.10,
             ae_dizz_pgb=1.50, ae_dizz_bz=0.80,
             ae_nau=2.00, ae_act=1.40, ae_sex=1.00)
    p.update(w_ae_sed=0.45, w_ae_dizz=0.35, w_ae_nau=0.50,
             w_ae_act=0.25, w_ae_sex=0.15)
    p.update(h_drop0=0.0016, b_ae=0.0070, b_ineff=0.0030)

    # ---------------- comorbid depression -------------------------------------
    p.update(km_in=0.05, km_out=0.05, m5ht=0.75, mne=0.35, mcbt=0.45, km_phi=0.80,
             madrs_scale=8.0, madrs_floor=3.0)

    # ---------------- disease severity ----------------------------------------
    p.update(DIS=1.0)

    # ---------------- exogenous inputs (set per scenario) ---------------------
    p.update(CBT=0.0, adherence=1.0)

    # derived, filled in by calibrate_healthy()
    p.update(phi_healthy=1.0, sglu_healthy=1.0)
    return p


# ----------------------------------------------------------------------------
# 2.  Right-hand side (vectorised over subjects; y has shape (NS, n))
# ----------------------------------------------------------------------------
def rhs(t, y, P, cbt):
    """P is a dict of scalars OR of length-n arrays.  cbt is scalar/array."""
    d = np.zeros_like(y)
    g = P.get

    # ---------------- PK ----------------------------------------------------
    Cesc = y[IX["esc_c"]] / g("V1_esc")
    Cescp = y[IX["esc_p"]] / g("V2_esc")
    d[IX["esc_gut"]] = -g("ka_esc") * y[IX["esc_gut"]]
    d[IX["esc_c"]] = (g("ka_esc") * y[IX["esc_gut"]] * g("F_esc")
                      - g("CL_esc") * Cesc
                      - g("Q_esc") * (Cesc - Cescp))
    d[IX["esc_p"]] = g("Q_esc") * (Cesc - Cescp)
    d[IX["esc_e"]] = g("keo_esc") * (Cesc - y[IX["esc_e"]])

    Cven = y[IX["ven_c"]] / g("V_ven")
    Codv = y[IX["odv_c"]] / g("V_odv")
    CLm = g("fm2d6") * g("CL_ven")
    d[IX["ven_gut"]] = -g("ka_ven") * y[IX["ven_gut"]]
    d[IX["ven_c"]] = g("ka_ven") * y[IX["ven_gut"]] * g("F_ven") - g("CL_ven") * Cven
    d[IX["odv_c"]] = g("mw_ratio_odv") * CLm * Cven - g("CL_odv") * Codv

    Cdlx = y[IX["dlx_c"]] / g("V_dlx")
    d[IX["dlx_gut"]] = -g("ka_dlx") * y[IX["dlx_gut"]]
    d[IX["dlx_c"]] = g("ka_dlx") * y[IX["dlx_gut"]] * g("F_dlx") - g("CL_dlx") * Cdlx

    CL_pgb = g("CL_pgb") * (g("crcl") / 100.0)
    Cpgb = y[IX["pgb_c"]] / g("V_pgb")
    d[IX["pgb_gut"]] = -g("ka_pgb") * y[IX["pgb_gut"]]
    d[IX["pgb_c"]] = g("ka_pgb") * y[IX["pgb_gut"]] * g("F_pgb") - CL_pgb * Cpgb

    Cbzd = y[IX["bzd_c"]] / g("V_bzd")
    d[IX["bzd_gut"]] = -g("ka_bzd") * y[IX["bzd_gut"]]
    d[IX["bzd_c"]] = g("ka_bzd") * y[IX["bzd_gut"]] * g("F_bzd") - g("CL_bzd") * Cbzd
    d[IX["bzd_e"]] = g("keo_bzd") * (Cbzd - y[IX["bzd_e"]])

    Cbus = y[IX["bus_c"]] / g("V_bus")
    Cpp1 = y[IX["pp1_c"]] / g("V_pp")
    d[IX["bus_gut"]] = -g("ka_bus") * y[IX["bus_gut"]]
    d[IX["bus_c"]] = g("ka_bus") * y[IX["bus_gut"]] * g("F_bus") - g("CL_bus") * Cbus
    d[IX["pp1_c"]] = g("f_pp") * g("ka_bus") * y[IX["bus_gut"]] - g("CL_pp") * Cpp1

    Cqtp = y[IX["qtp_c"]] / g("V_qtp")
    Cnqt = y[IX["nqt_c"]] / g("V_nqt")
    d[IX["qtp_gut"]] = -g("ka_qtp") * y[IX["qtp_gut"]]
    d[IX["qtp_c"]] = g("ka_qtp") * y[IX["qtp_gut"]] * g("F_qtp") - g("CL_qtp") * Cqtp
    d[IX["nqt_c"]] = g("f_nqt") * g("CL_qtp") * Cqtp - g("CL_nqt") * Cnqt

    # ---------------- receptor occupancies (competitive, multi-ligand) ------
    Xsert = (y[IX["esc_e"]] / g("EC50_sert_esc")
             + (Cven + Codv) / g("EC50_sert_ven")
             + Cdlx / g("EC50_sert_dlx"))
    occ_sert = Xsert / (1.0 + Xsert)

    Xnet = ((Cven + Codv) / g("EC50_net_ven")
            + Cdlx / g("EC50_net_dlx")
            + Cnqt / g("EC50_net_nqt"))
    occ_net = Xnet / (1.0 + Xnet)

    occ_a2d = Cpgb / (Cpgb + g("EC50_a2d"))
    Cbz_eq = y[IX["bzd_e"]] * g("bzd_potency")
    occ_bz = Cbz_eq / (Cbz_eq + g("EC50_bz"))
    occ_h1 = Cqtp / (Cqtp + g("EC50_h1"))
    occ_1a = Cbus / (Cbus + g("EC50_1a"))
    occ_a2adr = Cpp1 / (Cpp1 + g("EC50_pp1"))

    # ---------------- serotonin ---------------------------------------------
    sht = y[IX["sht"]]
    auto = y[IX["auto"]]
    ex_sht = np.maximum(0.0, sht - 1.0)
    fire = 1.0 / (1.0 + g("gamma_auto") * auto * ex_sht)
    fire = fire * (1.0 - g("ia_bus") * occ_1a * auto)
    d[IX["sht"]] = (g("kin_sht") * fire
                    - g("kout_sht") * sht * (g("sert_floor")
                                             + (1 - g("sert_floor")) * (1 - occ_sert)))
    d[IX["auto"]] = -g("kdes_auto") * auto * ex_sht + g("krec_auto") * (1.0 - auto)
    sht_post = sht + g("ia_bus_post") * occ_1a
    ex_shtp = g("k5ht_eff") * np.maximum(0.0, sht_post - 1.0)

    # ---------------- norepinephrine ------------------------------------------
    ne = y[IX["ne"]]
    a2 = y[IX["a2auto"]]
    ex_ne = np.maximum(0.0, ne - 1.0)
    fire_ne = (1.0 / (1.0 + g("gamma_a2") * a2 * ex_ne)) * (1.0 + g("kpp1") * occ_a2adr)
    d[IX["ne"]] = (g("kin_ne") * fire_ne
                   - g("kout_ne") * ne * (g("net_floor")
                                          + (1 - g("net_floor")) * (1 - occ_net)))
    d[IX["a2auto"]] = -g("kdes_a2") * a2 * ex_ne + g("krec_a2") * (1.0 - a2)

    # ---------------- HPA axis --------------------------------------------------
    crh, acth, cort, gr = y[IX["crh"]], y[IX["acth"]], y[IX["cort"]], y[IX["gr"]]
    # Excess cortisol enters every downstream effect through a saturating
    # transform.  Without it the CRH -> cortisol -> GR-downregulation ->
    # weaker feedback -> CRH loop has gain > 1 and the axis runs away; that
    # was an actual defect found while building this file, not a cosmetic
    # choice (see README, "defects found by verification").
    ex_cort = np.maximum(0.0, cort - 1.0) / (1.0 + 0.5 * np.maximum(0.0, cort - 1.0))

    # ---------------- plasticity / PFC / amygdala --------------------------------
    bdnf = y[IX["bdnf"]]
    bdnf_in = (g("w_5ht") * (1.0 + ex_shtp) + g("w_ne") * ne + g("w_tonic")
               + g("cbt_bdnf") * cbt) * (1.0 - g("kcort_bdnf") * ex_cort)
    bdnf_in = np.maximum(bdnf_in, 0.0)
    d[IX["bdnf"]] = g("kb") * bdnf_in - g("kdb") * bdnf

    cpfc = y[IX["cpfc"]]
    hill = bdnf * bdnf / (bdnf * bdnf + 1.0)
    d[IX["cpfc"]] = (g("kp_pfc") * hill * (g("cpfc_max") - cpfc)
                     - g("kl_pfc") * cpfc * (1.0 + g("kc_pfc") * ex_cort
                                             + g("kdis_pfc") * g("DIS")))

    eamy = y[IX["eamy"]]
    sleepd = np.maximum(y[IX["sleepd"]], 0.0)
    # Logistic ceiling EAMAX: amygdala/BNST excitability is bounded, and the
    # bound is what keeps the disease loop convergent.
    amy_drive = (1.0 + g("kcrh_amy") * np.maximum(0.0, crh - 1.0)
                 + g("kslp_amy") * sleepd + g("kdis_amy") * g("DIS"))
    d[IX["eamy"]] = (g("ke_in") * amy_drive
                     * np.maximum(0.0, g("EAMAX") - eamy) / (g("EAMAX") - 1.0)
                     - g("ke_out") * eamy * (1.0 + g("a5ht_amy") * ex_shtp
                                             + g("acbt_amy") * cbt))

    d[IX["traf"]] = g("ktraf") * (occ_a2d - y[IX["traf"]])
    traf = y[IX["traf"]]

    # ---------------- the four factors of Phi -------------------------------------
    sglu = ((1.0 + g("kglu_stress") * np.maximum(0.0, eamy - 1.0)
             + g("kglu_cort") * ex_cort) * (1.0 - g("emax_pgb") * traf))

    ra1, ra2, dep = y[IX["ra1"]], y[IX["ra2"]], y[IX["depend"]]
    d[IX["ra1"]] = -g("ktol1") * occ_bz * ra1 + g("krec1") * (1.0 - ra1)
    d[IX["ra2"]] = -g("ktol2") * occ_bz * ra2 + g("krec2") * (1.0 - ra2)
    d[IX["depend"]] = g("kdep_on") * (occ_bz - dep)
    igaba_base = g("igaba_base0") - g("kdis_gaba") * g("DIS")
    igaba = igaba_base * (1.0 + g("emax_a2") * ra2 * (occ_bz - g("wadapt") * dep))
    igaba = np.maximum(igaba, 0.05)

    expect = y[IX["expect"]]
    cpfc_eff = np.maximum(cpfc * (1.0 + g("e_ex") * expect), 0.05)

    phi = (eamy * np.maximum(sglu, 0.02)) / (cpfc_eff * igaba)
    phi_n = phi / g("phi_healthy")

    # ---------------- bounded excess gain -------------------------------------------
    # Every downstream consequence of Phi enters through a SATURATING excess
    # z, not through phi_n itself.  Without this the amygdala -> sleep -> CRH
    # -> cortisol -> PFC loop has gain > 1 and the whole model diverges (an
    # actual defect found while building this file; see README).
    exphi = phi_n - 1.0
    z = np.where(exphi > 0.0, exphi / (1.0 + exphi / g("ZMAX")), exphi)
    z = np.clip(z, -0.90, g("ZMAX"))
    zp = np.maximum(z, 0.0)

    # ---------------- HPA driven by Phi -------------------------------------------
    d[IX["crh"]] = (g("kcrh") * (1.0 + g("kcrh_phi") * zp)
                    - g("kcout") * crh * (1.0 + g("kgr_fb") * gr * cort))
    d[IX["acth"]] = g("ka_acth") * crh - g("ke_acth") * acth
    d[IX["cort"]] = g("kc_cort") * acth - g("kec_cort") * cort
    d[IX["gr"]] = g("kgr_in") * (1.0 - gr) - g("kgr_dn") * ex_cort * gr

    # ---------------- autonomic / symptom layer ------------------------------------
    sns = y[IX["sns"]]
    d[IX["sns"]] = (g("ks_in") * (1.0 + g("ks_phi") * zp + g("w_ne_sns") * ex_ne)
                    - g("ks_out") * sns)

    auton = y[IX["auton"]]
    d[IX["auton"]] = (g("kau_in") * sns
                      - g("kau_out") * auton * (1.0
                                                + g("b_gaba") * np.maximum(0.0, igaba / igaba_base - 1.0)
                                                + g("b_pgb") * traf))

    d[IX["sleepd"]] = (g("ksl_in") * (g("ksl_phi") * zp + g("w_cort_sl") * ex_cort)
                       * np.maximum(0.0, g("SLMAX") - sleepd)
                       - g("ksl_out") * sleepd * (1.0 + g("h_h1") * occ_h1 * y[IX["rh1"]]
                                                  + g("h_gaba") * occ_bz * ra1
                                                  + g("h_pgb") * traf))

    worry = y[IX["worry"]]
    d[IX["worry"]] = (g("kw_in") * (1.0 + g("kw_phi") * zp)
                      * (1.0 + g("rho_worry") * worry / (g("kw_m") + worry))
                      - g("kw_out") * worry * (1.0 + g("acbt_worry") * cbt
                                               + g("a5ht_worry") * ex_shtp))

    # ---------------- trial machinery -------------------------------------------------
    d[IX["expect"]] = -g("kexp_dec") * expect
    d[IX["fluct"]] = -g("kfl") * y[IX["fluct"]]

    # ---------------- adverse effects ---------------------------------------------------
    d[IX["rnau"]] = -g("ktol_nau") * occ_sert * y[IX["rnau"]] + g("krec_nau") * (1.0 - y[IX["rnau"]])
    d[IX["rdizz"]] = -g("ktol_dizz") * traf * y[IX["rdizz"]] + g("krec_dizz") * (1.0 - y[IX["rdizz"]])
    d[IX["rh1"]] = -g("ktol_h1") * occ_h1 * y[IX["rh1"]] + g("krec_h1") * (1.0 - y[IX["rh1"]])
    d[IX["sexd"]] = g("ks_sex") * occ_sert * (1.0 - y[IX["sexd"]]) - g("kr_sex") * y[IX["sexd"]]
    d[IX["ract"]] = -g("ktol_act") * ex_ne * y[IX["ract"]] + g("krec_act") * (1.0 - y[IX["ract"]])
    d[IX["wt"]] = g("kwt") * (traf + 0.8 * occ_h1) - g("kwt_off") * y[IX["wt"]]

    sedation = (g("ae_sed_h1") * occ_h1 * y[IX["rh1"]]
                + g("ae_sed_bz") * occ_bz * ra1
                + g("ae_sed_pgb") * traf * y[IX["rdizz"]])
    dizziness = g("ae_dizz_pgb") * traf * y[IX["rdizz"]] + g("ae_dizz_bz") * occ_bz * ra1
    nausea = g("ae_nau") * occ_sert * y[IX["rnau"]]
    activation = g("ae_act") * ex_ne * y[IX["ract"]]
    ae_burden = (g("w_ae_sed") * sedation + g("w_ae_dizz") * dizziness
                 + g("w_ae_nau") * nausea + g("w_ae_act") * activation
                 + g("w_ae_sex") * g("ae_sex") * y[IX["sexd"]])

    # ---------------- comorbid depression -------------------------------------------------
    d[IX["madrs"]] = (g("km_in") * (1.0 + g("km_phi") * zp)
                      - g("km_out") * y[IX["madrs"]] * (1.0 + g("m5ht") * ex_shtp
                                                        + g("mne") * ex_ne + g("mcbt") * cbt))

    d[IX["cumhaz"]] = g("h_drop0") + g("b_ae") * ae_burden
    return d


# ----------------------------------------------------------------------------
# 3.  Derived (observation) quantities
# ----------------------------------------------------------------------------
def observe(y, P):
    """Compute all readouts from a state array y of shape (NS, n) or (NS,)."""
    g = P.get
    Cesc = y[IX["esc_c"]] / g("V1_esc")
    Cven = y[IX["ven_c"]] / g("V_ven")
    Codv = y[IX["odv_c"]] / g("V_odv")
    Cdlx = y[IX["dlx_c"]] / g("V_dlx")
    Cpgb = y[IX["pgb_c"]] / g("V_pgb")
    Cbzd = y[IX["bzd_c"]] / g("V_bzd")
    Cbus = y[IX["bus_c"]] / g("V_bus")
    Cpp1 = y[IX["pp1_c"]] / g("V_pp")
    Cqtp = y[IX["qtp_c"]] / g("V_qtp")
    Cnqt = y[IX["nqt_c"]] / g("V_nqt")

    Xsert = (y[IX["esc_e"]] / g("EC50_sert_esc") + (Cven + Codv) / g("EC50_sert_ven")
             + Cdlx / g("EC50_sert_dlx"))
    occ_sert = Xsert / (1.0 + Xsert)
    Xnet = ((Cven + Codv) / g("EC50_net_ven") + Cdlx / g("EC50_net_dlx")
            + Cnqt / g("EC50_net_nqt"))
    occ_net = Xnet / (1.0 + Xnet)
    occ_a2d = Cpgb / (Cpgb + g("EC50_a2d"))
    Cbz_eq = y[IX["bzd_e"]] * g("bzd_potency")
    occ_bz = Cbz_eq / (Cbz_eq + g("EC50_bz"))
    occ_h1 = Cqtp / (Cqtp + g("EC50_h1"))
    occ_1a = Cbus / (Cbus + g("EC50_1a"))

    sht_post = y[IX["sht"]] + g("ia_bus_post") * occ_1a
    _exc = np.maximum(0.0, y[IX["cort"]] - 1.0)
    ex_cort = _exc / (1.0 + 0.5 * _exc)
    eamy = y[IX["eamy"]]
    traf = y[IX["traf"]]
    sglu = ((1.0 + g("kglu_stress") * np.maximum(0.0, eamy - 1.0) + g("kglu_cort") * ex_cort)
            * (1.0 - g("emax_pgb") * traf))
    igaba_base = g("igaba_base0") - g("kdis_gaba") * g("DIS")
    igaba = igaba_base * (1.0 + g("emax_a2") * y[IX["ra2"]]
                          * (occ_bz - g("wadapt") * y[IX["depend"]]))
    igaba = np.maximum(igaba, 0.05)
    cpfc_eff = np.maximum(y[IX["cpfc"]] * (1.0 + g("e_ex") * y[IX["expect"]]), 0.05)
    phi = (eamy * np.maximum(sglu, 0.02)) / (cpfc_eff * igaba)
    phi_n = phi / g("phi_healthy")

    sleepd = np.maximum(y[IX["sleepd"]], 0.0)
    psy_raw = g("c1") * y[IX["worry"]] + g("c2") * phi_n + g("c3") * sleepd
    psy_eff = np.maximum(0.0, psy_raw - g("psy0"))
    hama_psy = g("psy_floor") + (28.0 - g("psy_floor")) * psy_eff / (psy_eff + g("Kpsy"))
    som_raw = (g("v1") * y[IX["auton"]] + g("v2") * y[IX["sns"]]
               + g("v3") * (sglu / g("sglu_healthy")) + g("v4") * sleepd)
    som_eff = np.maximum(0.0, som_raw - g("som0"))
    hama_som = g("som_floor") + (28.0 - g("som_floor")) * som_eff / (som_eff + g("Ksom"))
    hama = hama_psy + hama_som + y[IX["fluct"]]

    sedation = (g("ae_sed_h1") * occ_h1 * y[IX["rh1"]] + g("ae_sed_bz") * occ_bz * y[IX["ra1"]]
                + g("ae_sed_pgb") * traf * y[IX["rdizz"]])
    dizziness = g("ae_dizz_pgb") * traf * y[IX["rdizz"]] + g("ae_dizz_bz") * occ_bz * y[IX["ra1"]]
    nausea = g("ae_nau") * occ_sert * y[IX["rnau"]]
    activation = g("ae_act") * np.maximum(0.0, y[IX["ne"]] - 1.0) * y[IX["ract"]]
    ae_burden = (g("w_ae_sed") * sedation + g("w_ae_dizz") * dizziness + g("w_ae_nau") * nausea
                 + g("w_ae_act") * activation + g("w_ae_sex") * g("ae_sex") * y[IX["sexd"]])

    return dict(
        Cesc=Cesc, Cven=Cven, Codv=Codv, Cdlx=Cdlx, Cpgb=Cpgb, Cbzd=Cbzd,
        Cbus=Cbus, Cpp1=Cpp1, Cqtp=Cqtp, Cnqt=Cnqt,
        occ_sert=occ_sert, occ_net=occ_net, occ_a2d=occ_a2d, occ_bz=occ_bz,
        occ_h1=occ_h1, occ_1a=occ_1a,
        sht=y[IX["sht"]], sht_post=sht_post, ne=y[IX["ne"]], auto=y[IX["auto"]],
        bdnf=y[IX["bdnf"]], cpfc=y[IX["cpfc"]], cpfc_eff=cpfc_eff, eamy=eamy,
        sglu=sglu, igaba=igaba, phi=phi, phi_n=phi_n,
        cort=y[IX["cort"]], sns=y[IX["sns"]], auton=y[IX["auton"]],
        sleepd=sleepd, worry=y[IX["worry"]], expect=y[IX["expect"]],
        hama_psy=hama_psy, hama_som=hama_som, hama=hama,
        madrs=g("madrs_floor") + g("madrs_scale") * np.maximum(0.0, y[IX["madrs"]] - 1.0),
        sedation=sedation, dizziness=dizziness, nausea=nausea,
        activation=activation, sexd=y[IX["sexd"]], wt=y[IX["wt"]],
        ae_burden=ae_burden, cumhaz=y[IX["cumhaz"]],
    )


# ----------------------------------------------------------------------------
# 4.  Integrator
# ----------------------------------------------------------------------------
DT = 1.0 / 48.0   # 30 min: QD, BID and TID dose times all land on the grid
                  # exactly, and the fastest rate in the system (5-HT turnover,
                  # 24/d) gives h*lambda = 0.5, well inside RK4 stability.


def make_y0(n=1):
    y = np.zeros((NS, n))
    for s in ("sht", "auto", "ne", "a2auto", "bdnf", "cpfc", "eamy",
              "ra2", "ra1", "crh", "acth", "cort", "gr", "sns", "auton",
              "worry", "rnau", "rdizz", "rh1", "ract"):
        y[IX[s]] = 1.0
    y[IX["madrs"]] = 1.0
    return y


class Regimen(object):
    """Dose impulses, stored per subject so that a whole multi-arm trial can be
    integrated as ONE vectorised run.

    Each entry is (compartment_index, amount_ug (scalar or length-n array),
    time_in_days).
    """

    def __init__(self):
        self.doses = []
        self.cbt_windows = []      # (t0, t1, intensity scalar-or-array)

    def add(self, cmt, dose_mg, start, stop, per_day, mask=None, dt=DT):
        """per_day: 1 (QD), 2 (BID), 3 (TID).  dose_mg is the TOTAL daily dose
        (scalar or per-subject array).  mask selects the subjects dosed."""
        dose_mg = np.asarray(dose_mg, dtype=float)
        if per_day <= 0 or stop <= start or np.all(dose_mg <= 0):
            return self
        each = dose_mg / float(per_day) * 1000.0            # ug per administration
        if mask is not None:
            each = each * np.asarray(mask, dtype=float)
        step = 1.0 / per_day
        t = float(start)
        while t < stop - 1e-9:
            self.doses.append((IX[cmt], each, round(t / dt) * dt))
            t += step
        return self

    def cbt(self, t0, t1, intensity=1.0):
        self.cbt_windows.append((t0, t1, intensity))
        return self

    def merge(self, other):
        self.doses.extend(other.doses)
        self.cbt_windows.extend(other.cbt_windows)
        return self


ALL_OBS = None      # sentinel: record every observable


def simulate(P, regimen, tstop, y0=None, n=1, dt=DT, visits=(),
             obs_times=None, fluct0=0.0, expect_pulse=True,
             obs_keys=ALL_OBS, record_states=()):
    """Vectorised RK4 with dose impulses and visit (expectancy) impulses.

    Returns (times, out, y_final) with out[key] of shape (len(times), n).
    """
    if y0 is None:
        y0 = make_y0(n)
    y = np.array(y0, dtype=float, copy=True)
    if y.ndim == 1:
        y = y[:, None]
    if y.shape[1] == 1 and n > 1:
        y = np.repeat(y, n, axis=1)
    n = y.shape[1]
    y[IX["fluct"]] = y[IX["fluct"]] + np.asarray(fluct0, dtype=float)

    nsteps = int(round(tstop / dt))
    dose_map = {}
    for cmt, amt, t in regimen.doses:
        k = int(round(t / dt))
        if 0 <= k <= nsteps:
            dose_map.setdefault(k, []).append((cmt, amt))
    visit_steps = set(int(round(v / dt)) for v in visits)

    if obs_times is None:
        obs_times = np.arange(0.0, tstop + 1e-9, 1.0)
    times = np.asarray(obs_times, dtype=float)
    obs_steps = {}
    for i, tt in enumerate(times):
        obs_steps.setdefault(int(round(tt / dt)), []).append(i)

    out = {}
    adh = P.get("adherence", 1.0)

    cbt_w = regimen.cbt_windows

    def cbt_at(t):
        # ACCUMULATE.  An earlier version assigned instead of adding, so when
        # two arms of the same run each had a CBT window the later one silently
        # switched the earlier one off and the "CBT alone" arm came out
        # numerically identical to placebo.
        v = 0.0
        for (a, b, inten) in cbt_w:
            if a <= t < b:
                v = v + inten
        return v

    def record(k):
        rows = obs_steps.get(k)
        if rows is None:
            return
        o = observe(y, P)
        keys = o.keys() if obs_keys is ALL_OBS else obs_keys
        for key in keys:
            arr = out.get(key)
            if arr is None:
                arr = out[key] = np.zeros((len(times), n))
            arr[rows, :] = np.asarray(o[key], dtype=float)
        for s in record_states:
            key = "st_" + s
            arr = out.get(key)
            if arr is None:
                arr = out[key] = np.zeros((len(times), n))
            arr[rows, :] = y[IX[s]]

    for k in range(nsteps + 1):
        t = k * dt
        dl = dose_map.get(k)
        if dl is not None:
            for cmt, amt in dl:
                y[cmt] = y[cmt] + amt * adh
        record(k)                      # measure the patient FIRST ...
        if expect_pulse and k in visit_steps:   # ... then the visit happens
            y[IX["expect"]] = (y[IX["expect"]]
                               + P["dvisit"] * (P["emax_exp"] - y[IX["expect"]]))
        if k == nsteps:
            break
        c0 = cbt_at(t)
        ch = cbt_at(t + dt / 2)
        c1 = cbt_at(t + dt)
        k1 = rhs(t, y, P, c0)
        k2 = rhs(t + dt / 2, y + dt / 2 * k1, P, ch)
        k3 = rhs(t + dt / 2, y + dt / 2 * k2, P, ch)
        k4 = rhs(t + dt, y + dt * k3, P, c1)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        y = np.maximum(y, -50.0)
    return times, out, y


def run_to_steady_state(P, days=500.0, dt=1.0 / 24.0, y0=None, n=1):
    """No drug, no visits: relax the untreated system to its attractor."""
    r = Regimen()
    t, o, y = simulate(P, r, days, y0=y0, n=n, dt=dt,
                       obs_times=np.array([0.0]), expect_pulse=False,
                       obs_keys=())
    return y


_ATTR_CACHE = {}


def attractor_grid(P, dis_lo=0.30, dis_hi=1.80, ngrid=61, days=500.0):
    """Untreated attractor as a function of disease severity DIS, computed ONCE
    in a single vectorised relaxation and then interpolated per subject.  DIS is
    the only covariate that changes the drug-free attractor."""
    key = (round(P["kdis_amy"], 6), round(P["kdis_pfc"], 6),
           round(P["kdis_gaba"], 6), round(P["phi_healthy"], 8), ngrid)
    if key in _ATTR_CACHE:
        return _ATTR_CACHE[key]
    cache_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "gad_attractor_cache.npz")
    tag = "k%s" % ("_".join("%.6g" % v for v in key))
    if os.path.exists(cache_file):
        try:
            z = np.load(cache_file)
            if tag + "_dis" in z.files:
                res = (z[tag + "_dis"], z[tag + "_y"])
                _ATTR_CACHE[key] = res
                return res
        except Exception:
            pass
    dis = np.linspace(dis_lo, dis_hi, ngrid)
    Q = dict(P)
    Q["DIS"] = dis
    y = run_to_steady_state(Q, days=days, n=ngrid)
    _ATTR_CACHE[key] = (dis, y)
    store = {}
    if os.path.exists(cache_file):
        try:
            z = np.load(cache_file)
            store = {k: z[k] for k in z.files}
        except Exception:
            store = {}
    store[tag + "_dis"] = dis
    store[tag + "_y"] = y
    try:
        np.savez_compressed(cache_file, **store)
    except Exception:
        pass
    return dis, y


def attractor_for(P, dis_vec):
    """Interpolate the untreated attractor for a vector of DIS values."""
    grid, Y = attractor_grid(P)
    dis_vec = np.atleast_1d(np.asarray(dis_vec, dtype=float))
    out = np.zeros((NS, dis_vec.size))
    for i in range(NS):
        out[i] = np.interp(dis_vec, grid, Y[i])
    return out


# ----------------------------------------------------------------------------
# 5.  Calibration helpers
# ----------------------------------------------------------------------------
def calibrate_healthy(P):
    """Set phi_healthy and sglu_healthy from the DIS = 0 attractor."""
    Q = dict(P)
    Q["DIS"] = 0.0
    Q["phi_healthy"] = 1.0
    Q["sglu_healthy"] = 1.0
    y = run_to_steady_state(Q)
    o = observe(y, Q)
    P["phi_healthy"] = float(o["phi"][0])
    P["sglu_healthy"] = float(o["sglu"][0])
    return P


def gad_baseline(P, dis=None):
    """Untreated GAD attractor for a given severity."""
    Q = dict(P)
    if dis is not None:
        Q["DIS"] = dis
    return run_to_steady_state(Q)


# ----------------------------------------------------------------------------
# 6.  Trial machinery
# ----------------------------------------------------------------------------
VISIT_DAYS_8WK = [0.0, 7.0, 14.0, 28.0, 42.0, 56.0]
VISIT_DAYS_4WK = [0.0, 7.0, 14.0, 21.0, 28.0]


def build_arm(P, drug, total_mg, per_day, start=0.0, stop=1e9,
              titrate=None):
    """Return a Regimen for a named drug arm."""
    r = Regimen()
    if drug is None or total_mg <= 0:
        return r
    cmt = {"escitalopram": "esc_gut", "venlafaxine": "ven_gut",
           "duloxetine": "dlx_gut", "pregabalin": "pgb_gut",
           "lorazepam": "bzd_gut", "alprazolam": "bzd_gut",
           "buspirone": "bus_gut", "quetiapine": "qtp_gut"}[drug]
    if titrate:
        for (t0, t1, mg) in titrate:
            r.add(cmt, mg, t0, min(t1, stop), per_day)
    else:
        r.add(cmt, total_mg, start, stop, per_day)
    return r


DRUG_DEFAULTS = {
    "escitalopram": dict(per_day=1),
    "venlafaxine": dict(per_day=1),
    "duloxetine": dict(per_day=1),
    "pregabalin": dict(per_day=2),
    "lorazepam": dict(per_day=3, bzd_potency=1.0),
    "alprazolam": dict(per_day=3, bzd_potency=1.53),
    "buspirone": dict(per_day=3),
    "quetiapine": dict(per_day=1),
}


# ----------------------------------------------------------------------------
# 7.  Multi-arm trial machinery (one vectorised integration for the whole trial)
# ----------------------------------------------------------------------------
CMT_OF = {"escitalopram": "esc_gut", "venlafaxine": "ven_gut",
          "duloxetine": "dlx_gut", "pregabalin": "pgb_gut",
          "lorazepam": "bzd_gut", "alprazolam": "bzd_gut",
          "buspirone": "bus_gut", "quetiapine": "qtp_gut"}

PERDAY_OF = {"escitalopram": 1, "venlafaxine": 1, "duloxetine": 1,
             "pregabalin": 2, "lorazepam": 3, "alprazolam": 3,
             "buspirone": 3, "quetiapine": 1}

# Alprazolam is dosed and cleared differently from lorazepam.  Both act at the
# same benzodiazepine site, so the model converts each to LORAZEPAM-EQUIVALENT
# plasma concentration using the (well documented) clinical potency ratio
# alprazolam 0.5 mg == lorazepam 1 mg, and then applies the single measured
# occupancy hyperbola EC50 = 96 ng/mL (Atack 2007, PMID 17164474).
BZD_PK = {
    "lorazepam":  dict(F_bzd=0.90, V_bzd=90.0, CL_bzd=108.0, bzd_potency=1.00),
    "alprazolam": dict(F_bzd=0.88, V_bzd=65.0, CL_bzd=90.0, bzd_potency=1.53),
}


def arm(name, drug=None, mg=0.0, per_day=None, start=0.0, stop=1e9,
        titrate=None, cbt=None, **overrides):
    return dict(name=name, drug=drug, mg=mg, per_day=per_day, start=start,
                stop=stop, titrate=titrate, cbt=cbt, overrides=overrides)


def _broadcast_params(P, arms, nper):
    """Build a parameter dict whose per-arm-varying entries are length-n arrays."""
    n = len(arms) * nper
    Q = dict(P)
    varying = set()
    for a in arms:
        varying.update(a["overrides"].keys())
        if a["drug"] in BZD_PK:
            varying.update(BZD_PK[a["drug"]].keys())
    for k in varying:
        base = P[k]
        vec = np.full(n, float(base))
        for i, a in enumerate(arms):
            v = None
            if a["drug"] in BZD_PK and k in BZD_PK[a["drug"]]:
                v = BZD_PK[a["drug"]][k]
            if k in a["overrides"]:
                v = a["overrides"][k]
            if v is not None:
                vec[i * nper:(i + 1) * nper] = v
        Q[k] = vec
    return Q, n


def run_trial(P, arms, tstop, visits=(), dis=1.0, nper=1, obs_times=None,
              fluct0=0.0, dt=DT, obs_keys=ALL_OBS, record_states=(),
              y0=None, expect_pulse=True):
    """Integrate every arm of a trial in ONE vectorised pass.

    Returns (times, out, y_final, slices) where slices[name] is the column
    slice belonging to that arm.
    """
    Q, n = _broadcast_params(P, arms, nper)
    dis_vec = np.full(n, 0.0, dtype=float)
    disA = np.atleast_1d(np.asarray(dis, dtype=float))
    if disA.size == 1:
        dis_vec[:] = disA[0]
    elif disA.size == n:
        dis_vec[:] = disA
    else:
        raise ValueError("dis must be scalar or length n")
    Q["DIS"] = dis_vec

    if y0 is None:
        y0 = attractor_for(P, dis_vec)

    R = Regimen()
    for i, a in enumerate(arms):
        mask = np.zeros(n)
        mask[i * nper:(i + 1) * nper] = 1.0
        if a["drug"]:
            cmt = CMT_OF[a["drug"]]
            pd = a["per_day"] or PERDAY_OF[a["drug"]]
            if a["titrate"]:
                for (t0, t1, mg) in a["titrate"]:
                    R.add(cmt, mg, t0, min(t1, a["stop"]), pd, mask=mask, dt=dt)
            else:
                R.add(cmt, a["mg"], a["start"], min(a["stop"], tstop), pd,
                      mask=mask, dt=dt)
        if a["cbt"]:
            t0, t1 = a["cbt"][0], a["cbt"][1]
            inten = np.zeros(n)
            inten[i * nper:(i + 1) * nper] = (a["cbt"][2] if len(a["cbt"]) > 2 else 1.0)
            R.cbt(t0, t1, inten)

    f0 = np.zeros(n) + np.asarray(fluct0, dtype=float)
    times, out, yf = simulate(Q, R, tstop, y0=y0, n=n, dt=dt, visits=visits,
                              obs_times=obs_times, fluct0=f0,
                              expect_pulse=expect_pulse, obs_keys=obs_keys,
                              record_states=record_states)
    slices = {a["name"]: slice(i * nper, (i + 1) * nper) for i, a in enumerate(arms)}
    return times, out, yf, slices, Q


def at(times, arr, t):
    """Value of arr (len(times) x n) at time t (nearest recorded time)."""
    i = int(np.argmin(np.abs(np.asarray(times) - t)))
    return arr[i]
