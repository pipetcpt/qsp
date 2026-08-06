#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tap_python_reference.py — independent reference implementation of the
toxic-alcohol-poisoning (methanol / ethylene glycol) QSP model.

WHY THIS FILE EXISTS
--------------------
The mrgsolve model (tap_mrgsolve_model.R) is the deliverable; this file is the
INDEPENDENT re-implementation used to verify it.  Every ODE is written here
from the same equations but in a different language, with a different
integrator (scipy LSODA), and the two are required to agree.  The numbers
quoted in README.md are produced by THIS file, so they can be reproduced
without an R installation:

    python3 tap_python_reference.py            # all sections
    python3 tap_python_reference.py --section 3

Units throughout: TIME = hours, AMOUNTS = mmol, CONCENTRATIONS = mmol/L (mM),
volumes = L, clearances = L/h, PaCO2 = mmHg, GFR = fraction of normal.
Clinical mg/dL values are produced only at the reporting boundary.
"""
from __future__ import annotations

import argparse
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# molecular weights (g/mol) and mg/dL <-> mM conversions
# ---------------------------------------------------------------------------
MW = dict(meoh=32.042, form=46.025, eg=62.068, glyc=76.051, oxal=90.034,
          etoh=46.068, fom=82.104)


def mgdl_to_mM(x, sp):
    return x * 10.0 / MW[sp]


def mM_to_mgdl(x, sp):
    return x * MW[sp] / 10.0


# ---------------------------------------------------------------------------
# state vector
# ---------------------------------------------------------------------------
SNAMES = [
    # --- methanol arm -----------------------------------------------------
    "A_GUTM",    # 0  gut methanol (mmol)
    "A_MEOH1",   # 1  central methanol (mmol)
    "A_MEOH2",   # 2  peripheral methanol (mmol)
    "A_FORM1",   # 3  central formate (mmol)
    "A_FORM2",   # 4  peripheral formate (mmol)
    "A_FCNS",    # 5  CNS formate (mmol)
    "A_FVIT",    # 6  vitreous/retinal formate (mmol)
    "THF",       # 7  relative hepatic tetrahydrofolate pool (1 = normal)
    # --- ethylene glycol arm ---------------------------------------------
    "A_GUTE",    # 8  gut ethylene glycol (mmol)
    "A_EG1",     # 9  central EG (mmol)
    "A_EG2",     # 10 peripheral EG (mmol)
    "A_GLYC1",   # 11 central glycolate (mmol)
    "A_GLYC2",   # 12 peripheral glycolate (mmol)
    "A_GLX",     # 13 glyoxylate (mmol)
    "A_OXAL",    # 14 plasma/ECF oxalate (mmol)
    "A_CAOXK",   # 15 renal CaOx deposited (mmol, cumulative)
    "A_CAOXT",   # 16 extrarenal CaOx deposited (mmol, cumulative)
    # --- antidotes and cofactors -----------------------------------------
    "A_GUTF",    # 17 gut fomepizole
    "A_FOM1",    # 18 central fomepizole
    "A_FOM2",    # 19 peripheral fomepizole
    "A_GUTETH",  # 20 gut ethanol
    "A_ETOH1",   # 21 central ethanol
    "A_ETOH2",   # 22 peripheral ethanol
    "FOLEF",     # 23 folinic-acid effect on Vmax_THF (0 = none)
    "THIEF",     # 24 thiamine effect on AST branch
    "PYREF",     # 25 pyridoxine effect on ALT branch
    # --- acid-base / electrolytes ----------------------------------------
    "HCO3",      # 26 plasma bicarbonate (mM)
    "LACT",      # 27 plasma L-lactate (mM)
    "PACO2",     # 28 arterial PaCO2 (mmHg)
    "CAION",     # 29 ionised calcium, pH-uncorrected pool (mM)
    # --- organ state ------------------------------------------------------
    "GFRF",      # 30 GFR as fraction of the patient's own baseline
    "PTINJ",     # 31 proximal tubular injury (0-1)
    "ATPC",      # 32 CNS ATP / energy charge (1 = normal)
    "PUT",       # 33 putaminal injury index (0-1)
    "OPTIC",     # 34 optic/retinal injury index (0-1)
    # --- cumulative bookkeeping and hazards -------------------------------
    "AUC_FCNS",  # 35 mM*h of CNS formate
    "AUC_GLYC",  # 36 mM*h of plasma glycolate
    "HAZ_MORT",  # 37 cumulative mortality hazard
    "HAZ_BLIND", # 38 cumulative blindness hazard
    "CUM_OX_M",  # 39 cumulative methanol oxidised (mmol)
    "CUM_OX_E",  # 40 cumulative EG oxidised (mmol)
    "CUM_HD_M",  # 41 methanol removed by dialysis (mmol)
    "CUM_HD_F",  # 42 formate removed by dialysis (mmol)
    "CUM_HD_E",  # 43 EG removed by dialysis (mmol)
    "CUM_HD_G",  # 44 glycolate removed by dialysis (mmol)
    "CUM_HD_FOM",# 45 fomepizole removed by dialysis (mmol)
    "CUM_UR_F",  # 46 formate excreted in urine (mmol)
    "CUM_UR_G",  # 47 glycolate excreted in urine (mmol)
    "CUM_UR_OX", # 48 oxalate excreted in urine (mmol)
    "CUM_BIC",   # 49 exogenous NaHCO3 delivered (infusion + dialysate + push)
]
IX = {n: i for i, n in enumerate(SNAMES)}
NS = len(SNAMES)


# ---------------------------------------------------------------------------
# parameters
# ---------------------------------------------------------------------------
def default_params(wt=70.0, gfr_base=1.0):
    p = {}
    p["WT"] = wt
    p["GFR0"] = 7.2 * gfr_base           # L/h  (120 mL/min = 7.2 L/h at 1.0)

    # -- distribution volumes (L) ------------------------------------------
    p["V1M"] = 0.35 * wt                 # methanol central
    p["V2M"] = 0.25 * wt                 # methanol peripheral (total 0.60 L/kg)
    p["V1E"] = 0.38 * wt                 # EG central
    p["V2E"] = 0.27 * wt                 # EG peripheral (total 0.65 L/kg)
    p["V1F"] = 0.30 * wt                 # formate central
    p["V2F"] = 0.20 * wt                 # formate peripheral (total 0.50 L/kg)
    p["V1G"] = 0.28 * wt                 # glycolate central
    p["V2G"] = 0.20 * wt                 # glycolate peripheral
    p["VCNS"] = 0.020 * wt               # brain water
    p["VVIT"] = 0.05                     # vitreous + retina (L, fixed)
    p["VECF"] = 0.20 * wt                # ECF, for Ca and oxalate
    p["VOX"] = 0.25 * wt                 # oxalate distribution
    p["V1FOM"] = 0.42 * wt               # fomepizole central
    p["V2FOM"] = 0.28 * wt               # fomepizole peripheral (0.70 L/kg)
    p["V1ET"] = 0.35 * wt
    p["V2ET"] = 0.25 * wt
    p["VLACT"] = 0.35 * wt

    # -- intercompartmental clearances (L/h) -------------------------------
    p["QM"] = 30.0
    p["QE"] = 30.0
    p["QF"] = 8.0                        # formate exchanges slowly -> rebound
    p["QG"] = 8.0
    p["QFOM"] = 25.0
    p["QET"] = 40.0

    # -- absorption --------------------------------------------------------
    p["KA"] = 1.8                        # /h, alcohols
    p["KAFOM"] = 2.5
    p["SLOWE"] = 0.65                    # ethanol slows gastric emptying
    p["KSLOWE"] = 15.0                   # mM ethanol for half of that effect

    # -- ADH: shared enzyme, competing substrates + competitive inhibitor --
    #    v_i = Vmax_i*(S_i/Km_i) / (1 + sum_j S_j/Km_j + F/Ki)
    p["VMAX_M"] = 95.0                   # mmol/h methanol oxidation capacity
    p["KM_M"] = 8.0                      # mM
    p["VMAX_E"] = 52.0                   # mmol/h EG oxidation capacity
    p["KM_E"] = 6.0                      # mM
    p["VMAX_ET"] = 150.0                 # mmol/h ethanol oxidation capacity
    p["KM_ET"] = 1.0                     # mM
    p["KI_FOM"] = 1.5e-4                 # mM (0.15 uM)
    p["FCAT_M"] = 0.06                   # catalase/CYP2E1 route, methanol
    p["FCAT_E"] = 0.01                   # EG oxidation is essentially all ADH

    # -- formate disposition ----------------------------------------------
    p["VMAX_THF"] = 18.0                 # mmol/h, folate-dependent oxidation
    p["KM_THF"] = 0.5                    # mM
    p["KTHF_USE"] = 0.0018               # THF consumed per mmol formate
    p["KTHF_REG"] = 0.030                # /h THF regeneration toward 1
    p["CLF_R0"] = 0.45                   # renal formate CL scale (x GFR frac)
    p["FR_MAX"] = 0.99                   # max fractional tubular reabsorption
    p["KFR"] = 2.5e-4                    # f_HA half-constant for reabsorption
    p["PKA_F"] = 3.75                    # formic acid
    p["PKA_G"] = 3.83                    # glycolic acid

    # -- CNS / vitreous transfer (only un-ionised acid crosses) ------------
    p["PSF"] = 900.0                     # L/h * f_HA  -> brain
    p["PSV"] = 8.0                       # L/h * f_HA  -> vitreous
    p["PHB_OFF"] = -0.40                 # brain pH = plasma pH + offset
    p["PHB_SLOPE"] = 0.50                # brain pH tracks 50% of plasma shift
    p["PHV_OFF"] = -0.30

    # -- glyoxylate branch -------------------------------------------------
    p["K_GLX"] = 4.0                     # /h  glyoxylate turnover
    p["KM_GLXO"] = 3.0                   # mM  glycolate oxidase Km
    p["VMAX_GLXO"] = 5.0                 # mmol/h - SLOW, which is exactly why
    #                                      glycolate is what accumulates
    p["BR_OX"] = 0.25                    # oxalate branch (the MINOR route)
    p["BR_TH"] = 0.45                    # AST/thiamine branch
    p["BR_PY"] = 0.30                    # ALT/pyridoxine branch
    p["ETH_TH"] = 1.10                   # extra AST branch on thiamine
    p["ETH_PY"] = 1.00                   # extra ALT branch on pyridoxine
    p["FRAC_EG_FORM"] = 0.04             # minor EG -> formate route

    # -- oxalate / calcium oxalate ----------------------------------------
    p["KSP"] = 2.32e-3                   # mM^2 (2.32e-9 M^2, CaOx monohydrate)
    p["CFTUB"] = 50.0                    # tubular concentration factor
    p["SSCRIT_K"] = 8.0                  # metastable limit, tubular fluid
    p["SSCRIT_T"] = 4.0                  # metastable limit, soft tissue
    p["FP_MAX"] = 0.60                   # max fraction of the filtered oxalate
    #                                      load that precipitates in the tubule
    p["KSS_K"] = 250.0                   # tubular SS excess for half of FP_MAX
    p["KDEP_T"] = 0.10                   # /h soft-tissue precipitation rate
    p["KSS_T"] = 20.0
    p["BUFCA"] = 0.40                    # fraction of Ca loss seen as ionised
    p["KCA_REP"] = 0.055                 # /h bone/PTH/gut restoration
    p["CA0"] = 1.18                      # mM ionised calcium baseline
    p["DCA_DPH"] = 0.36                  # mM ionised Ca fall per pH unit rise

    # -- acid-base ---------------------------------------------------------
    p["NA"] = 140.0
    p["CL"] = 104.0
    p["HCO30"] = 24.0
    p["AG0"] = 12.0
    p["TAU_CO2"] = 0.40                  # h, respiratory response lag
    p["PACO20"] = 40.0
    p["KHCO3_REN"] = 0.010               # /h renal HCO3 regeneration (NH4+ capacity)
    # Bicarbonate transfer across the dialyser is NOT at the small-solute
    # clearance: plasma equilibrates with the dialysate in minutes and the rest
    # of the (very large) buffer space refills slowly behind it.  Written at
    # 14.4 L/h the first version delivered >300 mmol/h and drove every dialysed
    # patient to HCO3 47 and PaCO2 79.
    p["CL_HD_BIC"] = 3.0                 # L/h effective, whole-body
    p["HCO3_DIAL"] = 35.0                # dialysate bicarbonate (mM)
    p["PH_TARGET"] = 7.35                # bicarbonate is titrated, not poured
    p["PH_TAU"] = 0.03
    p["VNA"] = 0.60                      # L/kg over which the Na load spreads
    p["KNA_OUT"] = 0.014                 # /h renal washout of the alkali load

    # -- lactate -----------------------------------------------------------
    p["LACT0"] = 1.0
    p["KLACT_OUT"] = 0.60                # /h
    p["KLACT_ATP"] = 5.0                 # anaerobic production gain
    p["KLACT_ETOH"] = 0.020              # NADH-driven rise per mM ethanol

    # -- mitochondrial / CNS injury ---------------------------------------
    p["KI_CCO"] = 12.0                   # mM formate at cytochrome c oxidase
    p["TAU_ATP"] = 0.50                  # h
    p["THR_PUT"] = 0.25                  # ATP deficit threshold for putamen
    p["K_PUT"] = 0.55                    # /h
    p["K_PUT_REC"] = 0.004
    p["THR_OPT_ATP"] = 0.18              # retinal ATP deficit that starts injury
    p["K_OPT"] = 0.40                    # /h at full deficit
    p["K_OPT_REC"] = 0.006
    p["LOGMAR_MAX"] = 2.0

    # -- renal injury ------------------------------------------------------
    p["K_PT_GALD"] = 0.0025              # per (mmol/h) of EG oxidation flux
    p["K_PT_REC"] = 0.040                # /h
    p["VE_THR"] = 6.0                    # mmol/h EG flux the tubule tolerates
    p["K_PT_XTAL"] = 0.35                # per (mmol/h) of renal CaOx deposition
    p["DEP_THR"] = 0.70                  # mmol/h CaOx the tubule can pass
    p["K_GFR_INJ"] = 0.22                # /h at PTINJ = 1
    p["K_GFR_REP"] = 0.035               # /h
    p["GFR_MIN"] = 0.05                  # anuric floor

    # -- sedation / respiratory failure -----------------------------------
    p["KSED_M"] = 60.0                   # mM methanol for unit sedation score
    p["KSED_E"] = 55.0
    p["KSED_ET"] = 45.0
    p["KSED_PH"] = 0.55                  # pH units below 7.2 -> unit score
    p["KSED_ATP"] = 0.45
    p["PACO2_FAIL"] = 65.0               # mmHg target when fully obtunded
    p["COMA_RESP"] = 0.72                # coma index at which failure sets in

    # -- hazards -----------------------------------------------------------
    p["H_PH"] = 0.16                     # /h per (7.20 - pH)^2 unit
    p["H_PUT"] = 0.030                   # /h per unit PUT
    p["H_RESP"] = 0.055                  # /h per unit respiratory failure
    p["H_BLIND"] = 0.85                  # /h per unit OPTIC above threshold
    p["BLIND_THR"] = 0.18
    p["PH_FATAL"] = 6.60                 # below this the run is non-survivable
    p["NOACID"] = 0.0                    # 1 = the metabolite is a KETONE
    #                                      (isopropanol control: same ADH flux,
    #                                       same osmolal gap, no acid at all)

    # -- fomepizole PK -----------------------------------------------------
    p["VMAX_FOM"] = 0.90                 # mmol/h
    p["KM_FOM"] = 0.05                   # mM
    p["AUTOIND"] = 0.60                  # fractional Vmax rise
    p["T_AUTO"] = 30.0                   # h at which autoinduction is half-on
    p["TAU_AUTO"] = 8.0

    # -- cofactor effect kinetics -----------------------------------------
    p["KOUT_FOL"] = 0.14                 # /h decay of folinate effect
    p["EMAX_FOL"] = 0.45                 # fractional rise of VMAX_THF
    p["KOUT_TH"] = 0.10
    p["KOUT_PY"] = 0.10

    # -- extracorporeal ----------------------------------------------------
    p["CL_HD"] = 14.4                    # L/h  (240 mL/min) small solutes
    p["CL_HD_FOM"] = 12.0                # L/h  fomepizole
    p["CL_HD_ET"] = 13.0                 # L/h  ethanol
    p["CL_CRRT"] = 2.4                   # L/h  (40 mL/min)
    return p


def y0_state(p):
    y = np.zeros(NS)
    y[IX["THF"]] = 1.0
    y[IX["HCO3"]] = p["HCO30"]
    y[IX["LACT"]] = p["LACT0"]
    y[IX["PACO2"]] = p["PACO20"]
    y[IX["CAION"]] = p["CA0"]
    y[IX["GFRF"]] = 1.0
    y[IX["ATPC"]] = 1.0
    y[IX["A_OXAL"]] = 0.002 * p["VOX"]   # ~2 uM endogenous plasma oxalate
    return y


# ---------------------------------------------------------------------------
# helper physiology
# ---------------------------------------------------------------------------
def f_HA(pH, pKa):
    """Fraction of a weak acid present as the un-ionised, membrane-permeant
    species.  This is the single most load-bearing function in the model."""
    return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def bicarb_space(hco3, wt):
    """Apparent bicarbonate space, L.  Garella/Fernandez-type expansion:
    space (L/kg) = 0.40 + 2.6/[HCO3].  The space GROWS as acidosis deepens,
    which is why plasma HCO3 decelerates on its way down."""
    h = max(hco3, 3.0)
    return (0.40 + 2.6 / h) * wt


def blood_pH(hco3, paco2):
    hco3 = max(hco3, 0.5)
    paco2 = max(paco2, 8.0)
    return 6.10 + math.log10(hco3 / (0.0301 * paco2))


def urine_pH(hco3):
    base = 5.0 + 3.0 * (hco3 / (hco3 + 20.0))
    # renal bicarbonate threshold: bicarbonaturia switches the urine alkaline
    extra = 1.30 / (1.0 + math.exp(-(hco3 - 27.0) / 1.2))
    return min(8.0, base + extra)


# ---------------------------------------------------------------------------
# right-hand side
# ---------------------------------------------------------------------------
def rhs(t, y, p, rx):
    y = np.maximum(y, 0.0)
    d = np.zeros(NS)

    # ---- concentrations -------------------------------------------------
    C_M = y[IX["A_MEOH1"]] / p["V1M"]
    C_M2 = y[IX["A_MEOH2"]] / p["V2M"]
    C_E = y[IX["A_EG1"]] / p["V1E"]
    C_E2 = y[IX["A_EG2"]] / p["V2E"]
    C_F = y[IX["A_FORM1"]] / p["V1F"]
    C_F2 = y[IX["A_FORM2"]] / p["V2F"]
    C_FC = y[IX["A_FCNS"]] / p["VCNS"]
    C_FV = y[IX["A_FVIT"]] / p["VVIT"]
    C_G = y[IX["A_GLYC1"]] / p["V1G"]
    C_G2 = y[IX["A_GLYC2"]] / p["V2G"]
    C_GX = y[IX["A_GLX"]] / p["VECF"]
    C_OX = y[IX["A_OXAL"]] / p["VOX"]
    C_FOM = y[IX["A_FOM1"]] / p["V1FOM"]
    C_FOM2 = y[IX["A_FOM2"]] / p["V2FOM"]
    C_ET = y[IX["A_ETOH1"]] / p["V1ET"]
    C_ET2 = y[IX["A_ETOH2"]] / p["V2ET"]

    hco3 = max(y[IX["HCO3"]], 0.5)
    paco2 = max(y[IX["PACO2"]], 8.0)
    pH = blood_pH(hco3, paco2)
    pHu = urine_pH(hco3)
    gfrf = max(y[IX["GFRF"]], 0.02)
    gfr = p["GFR0"] * gfrf

    # ---- ionised calcium as measured (pH-corrected) ---------------------
    ca_meas = max(0.30, y[IX["CAION"]] - p["DCA_DPH"] * (pH - 7.40))

    # ---- therapy switches -----------------------------------------------
    hd = rx["hd"](t)                 # 0 = off, 1 = IHD, 0.5 flag for CRRT
    cl_hd = (p["CL_HD"] if hd == 1 else (p["CL_CRRT"] if hd == 2 else 0.0))
    cl_hd_fom = (p["CL_HD_FOM"] if hd == 1 else (p["CL_CRRT"] if hd == 2 else 0.0))
    cl_hd_et = (p["CL_HD_ET"] if hd == 1 else (p["CL_CRRT"] if hd == 2 else 0.0))
    bic_inf = rx["bicarb"](t)        # mmol/h of NaHCO3
    fom_inf = rx["fom_inf"](t)       # mmol/h
    eth_inf = rx["eth_inf"](t)       # mmol/h
    ca_inf = rx["ca_inf"](t)         # mmol/h elemental Ca
    fol_inf = rx["fol"](t)           # arbitrary units/h
    thi_inf = rx["thi"](t)
    pyr_inf = rx["pyr"](t)

    # =====================================================================
    # 1 · ADH — one enzyme, three competing substrates, one inhibitor
    # =====================================================================
    sM = C_M / p["KM_M"]
    sE = C_E / p["KM_E"]
    sET = C_ET / p["KM_ET"]
    inh = C_FOM / p["KI_FOM"]
    den = 1.0 + sM + sE + sET + inh

    v_M = p["VMAX_M"] * sM / den
    v_E = p["VMAX_E"] * sE / den
    v_ET = p["VMAX_ET"] * sET / den
    # non-ADH routes are NOT blocked by fomepizole - the reason blockade is
    # never absolute even at enormous inhibitor concentrations
    v_M += p["FCAT_M"] * p["VMAX_M"] * C_M / (p["KM_M"] * 6.0 + C_M)
    v_E += p["FCAT_E"] * p["VMAX_E"] * C_E / (p["KM_E"] * 6.0 + C_E)

    # =====================================================================
    # 2 · absorption and parent-alcohol disposition
    # =====================================================================
    ka = p["KA"] / (1.0 + p["SLOWE"] * C_ET / (p["KSLOWE"] + C_ET))
    abs_M = ka * y[IX["A_GUTM"]]
    abs_E = ka * y[IX["A_GUTE"]]
    d[IX["A_GUTM"]] = -abs_M
    d[IX["A_GUTE"]] = -abs_E

    q_M = p["QM"] * (C_M - C_M2)
    q_E = p["QE"] * (C_E - C_E2)
    cl_ren_M = 0.050 * gfr                     # near-nil: reabsorbed with water
    cl_pulm_M = 0.15                           # L/h, blood:air 1350:1
    cl_ren_E = 0.26 * gfr                      # 20-30% of dose unchanged

    d[IX["A_MEOH1"]] = (abs_M - q_M - v_M
                        - (cl_ren_M + cl_pulm_M) * C_M - cl_hd * C_M)
    d[IX["A_MEOH2"]] = q_M
    d[IX["A_EG1"]] = (abs_E - q_E - v_E - cl_ren_E * C_E - cl_hd * C_E)
    d[IX["A_EG2"]] = q_E

    # =====================================================================
    # 3 · formate: production, folate-dependent oxidation, renal, CNS
    # =====================================================================
    vmax_thf = p["VMAX_THF"] * (1.0 + p["EMAX_FOL"] * y[IX["FOLEF"]])
    ox_form = vmax_thf * y[IX["THF"]] * C_F / (p["KM_THF"] + C_F)

    # renal handling: non-ionic back-diffusion of the un-ionised acid, so the
    # reabsorbed fraction is a function of URINE pH
    fha_u = f_HA(pHu, p["PKA_F"])
    fr_F = p["FR_MAX"] * fha_u / (fha_u + p["KFR"])
    cl_ren_F = gfr * (1.0 - fr_F) * p["CLF_R0"]
    exc_F = cl_ren_F * C_F

    fha_p = f_HA(pH, p["PKA_F"])
    phc = max(pH, 6.60)          # the offsets are not measured below this
    pH_b = 7.00 + p["PHB_SLOPE"] * (phc - 7.40)
    pH_v = 7.10 + p["PHB_SLOPE"] * (phc - 7.40)
    fha_b = f_HA(pH_b, p["PKA_F"])
    fha_v = f_HA(pH_v, p["PKA_F"])
    j_cns = p["PSF"] * (fha_p * C_F - fha_b * C_FC)
    j_vit = p["PSV"] * (fha_p * C_F - fha_v * C_FV)

    q_F = p["QF"] * (C_F - C_F2)
    form_from_eg = p["FRAC_EG_FORM"] * v_E
    acidM = 0.0 if p["NOACID"] > 0.5 else 1.0

    d[IX["A_FORM1"]] = (acidM * (v_M + form_from_eg) - q_F - ox_form - exc_F
                        - j_cns - j_vit - cl_hd * C_F)
    d[IX["A_FORM2"]] = q_F
    d[IX["A_FCNS"]] = j_cns
    d[IX["A_FVIT"]] = j_vit
    d[IX["THF"]] = (-p["KTHF_USE"] * ox_form
                    + p["KTHF_REG"] * (1.0 - y[IX["THF"]]))

    # =====================================================================
    # 4 · glycolate, glyoxylate branch ratio, oxalate
    # =====================================================================
    q_G = p["QG"] * (C_G - C_G2)
    v_glxo = p["VMAX_GLXO"] * C_G / (p["KM_GLXO"] + C_G)   # glycolate oxidase
    fha_ug = f_HA(pHu, p["PKA_G"])
    fr_G = p["FR_MAX"] * fha_ug / (fha_ug + p["KFR"])
    cl_ren_G = gfr * (1.0 - fr_G)
    exc_G = cl_ren_G * C_G

    d[IX["A_GLYC1"]] = (v_E * (1.0 - p["FRAC_EG_FORM"]) - q_G - v_glxo
                        - exc_G - cl_hd * C_G)
    d[IX["A_GLYC2"]] = q_G

    # branch ratio out of glyoxylate: oxalate vs the two cofactor-dependent
    # detoxifying routes.  Thiamine/pyridoxine change the RATIO, not the flux.
    b_ox = p["BR_OX"]
    b_th = p["BR_TH"] * (1.0 + p["ETH_TH"] * y[IX["THIEF"]])
    b_py = p["BR_PY"] * (1.0 + p["ETH_PY"] * y[IX["PYREF"]])
    btot = b_ox + b_th + b_py
    out_glx = p["K_GLX"] * y[IX["A_GLX"]]
    d[IX["A_GLX"]] = v_glxo - out_glx
    to_oxal = out_glx * b_ox / btot

    # Calcium oxalate supersaturation - the model's only threshold.  Crystal
    # formation cannot exceed the rate at which oxalate is DELIVERED, so the
    # tubular deposition is written as a saturating FRACTION of the filtered
    # load.  (Writing it as a rate proportional to the supersaturation excess
    # let the first version deposit 50 mmol/h of CaOx and drive every simulated
    # patient's ionised calcium to zero within hours.)
    ss_p = ca_meas * C_OX / p["KSP"]
    ss_k = ss_p * p["CFTUB"]
    ex_k = max(0.0, ss_k - p["SSCRIT_K"])
    ex_t = max(0.0, ss_p - p["SSCRIT_T"])
    filt_ox = gfr * C_OX
    fprec = p["FP_MAX"] * ex_k / (ex_k + p["KSS_K"])
    dep_k = fprec * filt_ox
    exc_OX = (1.0 - fprec) * filt_ox
    dep_t = p["KDEP_T"] * ex_t / (ex_t + p["KSS_T"]) * y[IX["A_OXAL"]]

    d[IX["A_OXAL"]] = to_oxal - dep_k - dep_t - exc_OX
    d[IX["A_CAOXK"]] = dep_k
    d[IX["A_CAOXT"]] = dep_t

    # =====================================================================
    # 5 · antidote PK
    # =====================================================================
    autoi = 1.0 + p["AUTOIND"] / (1.0 + math.exp(-(t - p["T_AUTO"]) / p["TAU_AUTO"]))
    el_fom = p["VMAX_FOM"] * autoi * C_FOM / (p["KM_FOM"] + C_FOM)
    abs_F = p["KAFOM"] * y[IX["A_GUTF"]]
    q_FOM = p["QFOM"] * (C_FOM - C_FOM2)
    d[IX["A_GUTF"]] = -abs_F
    d[IX["A_FOM1"]] = abs_F + fom_inf - q_FOM - el_fom - cl_hd_fom * C_FOM
    d[IX["A_FOM2"]] = q_FOM

    abs_ET = ka * y[IX["A_GUTETH"]]
    q_ET = p["QET"] * (C_ET - C_ET2)
    d[IX["A_GUTETH"]] = -abs_ET
    d[IX["A_ETOH1"]] = abs_ET + eth_inf - q_ET - v_ET - cl_hd_et * C_ET
    d[IX["A_ETOH2"]] = q_ET

    d[IX["FOLEF"]] = fol_inf - p["KOUT_FOL"] * y[IX["FOLEF"]]
    d[IX["THIEF"]] = thi_inf - p["KOUT_TH"] * y[IX["THIEF"]]
    d[IX["PYREF"]] = pyr_inf - p["KOUT_PY"] * y[IX["PYREF"]]

    # =====================================================================
    # 6 · acid-base
    # =====================================================================
    # Titratable acid, mmol/h, counted as ORGANIC ANION EQUIVALENTS RETAINED.
    #   +1 per formate produced           (formic acid)
    #   +1 per glycolate produced         (glycolic acid)
    #   +1 more per oxalate produced      (glyoxylate 1- -> oxalate 2-)
    #   -1 per formate oxidised to CO2    (bicarbonate regenerated)
    #   -1 per anion excreted in urine or removed in dialysate (with NH4+)
    #   -1 per glyoxylate routed to glycine / alpha-hydroxy-beta-ketoadipate
    #   Oxalate that PRECIPITATES as CaOx is deliberately NOT credited back:
    #   its proton has already left as CO2, so the bicarbonate deficit stays.
    neutralised = out_glx * (b_th + b_py) / btot
    acid_in = (acidM * (v_M + form_from_eg)
               + v_E * (1.0 - p["FRAC_EG_FORM"]) + to_oxal)
    acid_out = (ox_form + exc_F + exc_G + cl_hd * (C_F + C_G)
                + neutralised + 2.0 * exc_OX)
    net_acid = acid_in - acid_out
    lact_prod = (p["KLACT_OUT"] * p["LACT0"] * p["VLACT"]
                 + p["KLACT_ATP"] * max(0.0, 1.0 - y[IX["ATPC"]]) * p["VLACT"]
                 + p["KLACT_ETOH"] * C_ET * p["VLACT"])
    lact_out = p["KLACT_OUT"] * y[IX["LACT"]] * p["VLACT"]
    d[IX["LACT"]] = (lact_prod - lact_out) / p["VLACT"]
    net_acid += (lact_prod - lact_out)

    Vb = bicarb_space(hco3, p["WT"])
    # bicarbonate is titrated to a pH target, the way it is actually given
    _e = min(50.0, max(-50.0, (pH - p["PH_TARGET"]) / p["PH_TAU"]))
    bic_eff = bic_inf / (1.0 + math.exp(_e))
    hd_bic = p["CL_HD_BIC"] * (p["HCO3_DIAL"] - hco3) if hd else 0.0
    ren_bic = p["KHCO3_REN"] * gfrf * (p["HCO30"] - hco3) * Vb   # mmol/h, signed
    d[IX["HCO3"]] = (-net_acid + bic_eff + hd_bic + ren_bic) / Vb
    # Bicarbonaturia takes the sodium with it, so the exogenous alkali load and
    # the sodium load stay coupled; without this the osmolal gap drifted NEGATIVE
    # after a large bicarbonate infusion.  Written first as min(0, ren_bic) it
    # introduced a derivative kink exactly at HCO3 = 24, which LSODA's numerical
    # Jacobian turned into a small negative excursion of the cumulative state;
    # a smooth first-order washout is both cleaner and closer to the physiology.
    d[IX["CUM_BIC"]] = (bic_eff + max(0.0, hd_bic)
                        - p["KNA_OUT"] * gfrf * y[IX["CUM_BIC"]])

    # respiratory compensation, with failure as obtundation deepens
    sed = (C_M / p["KSED_M"] + C_E / p["KSED_E"] + C_ET / p["KSED_ET"]
           + max(0.0, 7.20 - pH) / p["KSED_PH"]
           + max(0.0, 1.0 - y[IX["ATPC"]]) / p["KSED_ATP"])
    coma = sed / (1.0 + sed)
    respf = 1.0 / (1.0 + math.exp(-(coma - p["COMA_RESP"]) / 0.06))
    # Winter's formula is a METABOLIC ACIDOSIS rule.  Extrapolated above a
    # normal bicarbonate it predicted PaCO2 79 mmHg; metabolic alkalosis
    # compensates at only ~0.7 mmHg per mM.
    if hco3 <= p["HCO30"]:
        tgt_comp = max(12.0, 1.5 * hco3 + 8.0)
    else:
        tgt_comp = min(55.0, 40.0 + 0.7 * (hco3 - p["HCO30"]))
    tgt = tgt_comp * (1.0 - respf) + p["PACO2_FAIL"] * respf
    d[IX["PACO2"]] = (tgt - paco2) / p["TAU_CO2"]

    # calcium
    d[IX["CAION"]] = (-(dep_k + dep_t) * p["BUFCA"] / p["VECF"]
                      + ca_inf * p["BUFCA"] / p["VECF"]
                      + p["KCA_REP"] * (p["CA0"] - y[IX["CAION"]]))

    # =====================================================================
    # 7 · organ injury
    # =====================================================================
    # glycolaldehyde/glyoxylate toxicity is a THRESHOLD phenomenon: the tubule
    # detoxifies a low aldehyde flux (GSH, aldehyde dehydrogenases) and fails
    # above it.  Without this threshold the residual non-ADH flux that survives
    # fomepizole slowly destroys the kidney of every simulated patient.
    hl_ald = v_E ** 3 / (v_E ** 3 + p["VE_THR"] ** 3)
    hl_xt = dep_k ** 3 / (dep_k ** 3 + p["DEP_THR"] ** 3)
    inj_drive = (p["K_PT_GALD"] * hl_ald * v_E
                 + p["K_PT_XTAL"] * hl_xt * dep_k)
    d[IX["PTINJ"]] = inj_drive * (1.0 - y[IX["PTINJ"]]) - p["K_PT_REC"] * y[IX["PTINJ"]]
    d[IX["GFRF"]] = (-p["K_GFR_INJ"] * y[IX["PTINJ"]] ** 2 * (gfrf - p["GFR_MIN"])
                     + p["K_GFR_REP"] * (1.0 - gfrf) * (1.0 - 0.70 * y[IX["PTINJ"]]))

    atp_tgt = 1.0 / (1.0 + C_FC / p["KI_CCO"] + C_GX / 8.0)
    d[IX["ATPC"]] = (atp_tgt - y[IX["ATPC"]]) / p["TAU_ATP"]

    deficit = max(0.0, (1.0 - y[IX["ATPC"]]) - p["THR_PUT"])
    d[IX["PUT"]] = p["K_PUT"] * deficit * (1.0 - y[IX["PUT"]]) - p["K_PUT_REC"] * y[IX["PUT"]]
    # the retina fails for the same reason the putamen does - loss of oxidative
    # phosphorylation - so the driver is a RETINAL ATP deficit, which makes
    # visual loss a threshold function of vitreous formate rather than a linear
    # one.  A linear driver made every prolonged low-level exposure blind.
    atp_ret = 1.0 / (1.0 + C_FV / p["KI_CCO"])
    def_opt = max(0.0, (1.0 - atp_ret) - p["THR_OPT_ATP"])
    d[IX["OPTIC"]] = (p["K_OPT"] * def_opt * (1.0 - y[IX["OPTIC"]])
                      - p["K_OPT_REC"] * y[IX["OPTIC"]])

    # =====================================================================
    # 8 · hazards and bookkeeping
    # =====================================================================
    d[IX["AUC_FCNS"]] = C_FC
    d[IX["AUC_GLYC"]] = C_G
    d[IX["HAZ_MORT"]] = (p["H_PH"] * max(0.0, 7.20 - pH) ** 2 * 25.0
                         + p["H_PUT"] * y[IX["PUT"]]
                         + p["H_RESP"] * respf)
    d[IX["HAZ_BLIND"]] = p["H_BLIND"] * max(0.0, y[IX["OPTIC"]] - p["BLIND_THR"])
    d[IX["CUM_OX_M"]] = v_M
    d[IX["CUM_OX_E"]] = v_E
    d[IX["CUM_HD_M"]] = cl_hd * C_M
    d[IX["CUM_HD_F"]] = cl_hd * C_F
    d[IX["CUM_HD_E"]] = cl_hd * C_E
    d[IX["CUM_HD_G"]] = cl_hd * C_G
    d[IX["CUM_HD_FOM"]] = cl_hd_fom * C_FOM
    d[IX["CUM_UR_F"]] = exc_F
    d[IX["CUM_UR_G"]] = exc_G
    d[IX["CUM_UR_OX"]] = exc_OX
    return d


# ---------------------------------------------------------------------------
# derived outputs
# ---------------------------------------------------------------------------
def outputs(t, y, p, rx):
    """Everything the clinician actually looks at."""
    o = {}
    C_M = y[IX["A_MEOH1"]] / p["V1M"]
    C_E = y[IX["A_EG1"]] / p["V1E"]
    C_F = y[IX["A_FORM1"]] / p["V1F"]
    C_G = y[IX["A_GLYC1"]] / p["V1G"]
    C_ET = y[IX["A_ETOH1"]] / p["V1ET"]
    C_FOM = y[IX["A_FOM1"]] / p["V1FOM"]
    C_OX = y[IX["A_OXAL"]] / p["VOX"]
    C_FC = y[IX["A_FCNS"]] / p["VCNS"]
    C_FV = y[IX["A_FVIT"]] / p["VVIT"]
    hco3 = max(y[IX["HCO3"]], 0.5)
    paco2 = max(y[IX["PACO2"]], 8.0)
    pH = blood_pH(hco3, paco2)

    o["t"] = t
    o["MEOH_mgdL"] = mM_to_mgdl(C_M, "meoh")
    o["EG_mgdL"] = mM_to_mgdl(C_E, "eg")
    o["ETOH_mgdL"] = mM_to_mgdl(C_ET, "etoh")
    o["FORM_mM"] = C_F
    o["FORM_mgL"] = C_F * MW["form"]
    o["GLYC_mM"] = C_G
    o["FOM_ugmL"] = C_FOM * MW["fom"]
    o["FCNS_mM"] = C_FC
    o["FVIT_mM"] = C_FV
    o["OXAL_uM"] = C_OX * 1000.0
    o["pH"] = pH
    o["HCO3"] = hco3
    o["PACO2"] = paco2
    o["pH_urine"] = urine_pH(hco3)
    o["LACT"] = y[IX["LACT"]]

    # --- the two gaps ----------------------------------------------------
    # Exogenous NaHCO3 raises sodium as well as bicarbonate.  Ignoring that made
    # the anion gap go to -11 on a bicarbonate infusion.  With the sodium
    # counted, the anion gap is almost UNCHANGED by bicarbonate therapy - which
    # is the correct and clinically useful behaviour: alkali fixes the pH
    # without erasing the diagnostic gap.
    dNa = y[IX["CUM_BIC"]] / (p["VNA"] * p["WT"])
    o["NA"] = p["NA"] + dNa
    o["AG"] = o["NA"] - p["CL"] - hco3
    # osmolal gap = measured - calculated.  Unmeasured solutes add to measured
    # only; bicarbonate is measured but absent from the calculated formula, so
    # it enters with +1; sodium enters measured with +1 and calculated with +2.
    unmeasured = (C_M + C_E + C_ET + C_F + C_G + max(0.0, y[IX["LACT"]] - p["LACT0"])
                  + C_OX)
    o["OG"] = unmeasured + (hco3 - p["HCO30"]) - dNa
    o["OG_etoh_corr"] = o["OG"] - C_ET
    o["GAPSUM"] = o["OG"] + (o["AG"] - p["AG0"])
    o["OG_over_AG"] = o["OG"] / max(0.5, o["AG"] - p["AG0"])

    # --- ADH blockade arithmetic ----------------------------------------
    sM = C_M / p["KM_M"]
    sE = C_E / p["KM_E"]
    sET = C_ET / p["KM_ET"]
    inh = C_FOM / p["KI_FOM"]
    o["ADH_den"] = 1.0 + sM + sE + sET + inh
    base = 1.0 + sM + sE
    o["INH_FACTOR"] = o["ADH_den"] / max(base, 1e-9)
    o["FLUX_M"] = p["VMAX_M"] * sM / o["ADH_den"]
    o["FLUX_E"] = p["VMAX_E"] * sE / o["ADH_den"]
    o["FLUX_FRAC_LEFT"] = 1.0 / o["INH_FACTOR"]

    # --- ion trapping ----------------------------------------------------
    o["fHA_plasma"] = f_HA(pH, p["PKA_F"])
    o["fHA_ratio_vs_740"] = o["fHA_plasma"] / f_HA(7.40, p["PKA_F"])
    pH_b = 7.00 + p["PHB_SLOPE"] * (pH - 7.40)
    o["CNS_eq_ratio"] = o["fHA_plasma"] / f_HA(pH_b, p["PKA_F"])

    # --- calcium, crystals, kidney ---------------------------------------
    ca_meas = max(0.30, y[IX["CAION"]] - p["DCA_DPH"] * (pH - 7.40))
    o["CA_ION"] = ca_meas
    o["CA_total_mgdL"] = ca_meas / 0.5 * 4.008          # ~50% ionised
    o["SS_plasma"] = ca_meas * C_OX / p["KSP"]
    o["SS_tubular"] = o["SS_plasma"] * p["CFTUB"]
    o["CAOX_KID"] = y[IX["A_CAOXK"]]
    o["GFR_pct"] = 100.0 * y[IX["GFRF"]]
    o["PTINJ"] = y[IX["PTINJ"]]
    o["QTc_ms"] = 400.0 + 95.0 * max(0.0, p["CA0"] - ca_meas)

    # --- CNS / eye / outcome ---------------------------------------------
    o["ATP"] = y[IX["ATPC"]]
    o["PUT"] = y[IX["PUT"]]
    o["OPTIC"] = y[IX["OPTIC"]]
    o["logMAR"] = p["LOGMAR_MAX"] * y[IX["OPTIC"]]
    sed = (C_M / p["KSED_M"] + C_E / p["KSED_E"] + C_ET / p["KSED_ET"]
           + max(0.0, 7.20 - pH) / p["KSED_PH"]
           + max(0.0, 1.0 - y[IX["ATPC"]]) / p["KSED_ATP"])
    o["COMA"] = sed / (1.0 + sed)
    o["P_DEATH"] = 1.0 - math.exp(-y[IX["HAZ_MORT"]])
    o["P_BLIND"] = 1.0 - math.exp(-y[IX["HAZ_BLIND"]])
    o["CUM_OX_M"] = y[IX["CUM_OX_M"]]
    o["CUM_OX_E"] = y[IX["CUM_OX_E"]]
    o["CUM_HD_M"] = y[IX["CUM_HD_M"]]
    o["CUM_HD_F"] = y[IX["CUM_HD_F"]]
    o["CUM_HD_FOM"] = y[IX["CUM_HD_FOM"]]
    o["THF"] = y[IX["THF"]]

    # --- diagnostic curiosities -----------------------------------------
    # glycolate cross-reacts on some lactate-oxidase point-of-care electrodes
    o["LACT_POC"] = y[IX["LACT"]] + 0.62 * C_G
    o["LACT_GAP"] = o["LACT_POC"] - y[IX["LACT"]]
    return o


# ---------------------------------------------------------------------------
# scenario machinery
# ---------------------------------------------------------------------------
def zero(t):
    return 0.0


def step_fn(windows, level=1.0):
    """Return f(t) = level inside any (start, stop) window, else 0."""
    def f(t):
        for a, b in windows:
            if a <= t < b:
                return level
        return 0.0
    return f


def hd_fn(windows, mode=1):
    def f(t):
        for a, b in windows:
            if a <= t < b:
                return mode
        return 0
    return f


class Scenario:
    def __init__(self, name, **kw):
        self.name = name
        self.wt = kw.get("wt", 70.0)
        self.gfr_base = kw.get("gfr_base", 1.0)
        self.meoh_g_per_kg = kw.get("meoh_g_per_kg", 0.0)
        self.eg_mL = kw.get("eg_mL", 0.0)          # mL of pure ethylene glycol
        self.etoh_g_per_kg = kw.get("etoh_g_per_kg", 0.0)   # co-ingested
        self.tend = kw.get("tend", 96.0)
        # therapies
        self.fom = kw.get("fom", None)   # dict(start=, bolus_mgkg=, q=, maint=, hd_q=)
        self.eth = kw.get("eth", None)   # dict(start=, load_gkg=, rate_mgkgh=)
        self.bicarb = kw.get("bicarb", None)  # dict(start=, stop=, rate_mmol_h=)
        self.hd = kw.get("hd", [])       # list of (start, stop)
        self.crrt = kw.get("crrt", [])
        self.folinate = kw.get("folinate", None)
        self.thiamine = kw.get("thiamine", None)
        self.pyridoxine = kw.get("pyridoxine", None)
        self.calcium = kw.get("calcium", None)
        self.fom_hd_boost = kw.get("fom_hd_boost", True)
        self.noacid = kw.get("noacid", 0.0)
        self.notes = kw.get("notes", "")

    # ------------------------------------------------------------------
    def build(self):
        p = default_params(self.wt, self.gfr_base)
        p["NOACID"] = self.noacid
        y = y0_state(p)
        boluses = []            # (time, index, amount mmol)

        if self.meoh_g_per_kg > 0:
            y[IX["A_GUTM"]] = self.meoh_g_per_kg * self.wt * 1000.0 / MW["meoh"]
        if self.eg_mL > 0:
            y[IX["A_GUTE"]] = self.eg_mL * 1.1132 * 1000.0 / MW["eg"]
        if self.etoh_g_per_kg > 0:
            y[IX["A_GUTETH"]] = self.etoh_g_per_kg * self.wt * 1000.0 / MW["etoh"]

        rx = dict(hd=lambda t: 0, bicarb=zero, fom_inf=zero, eth_inf=zero,
                  ca_inf=zero, fol=zero, thi=zero, pyr=zero)

        # ---- haemodialysis / CRRT windows ----
        hdw = list(self.hd)
        crw = list(self.crrt)
        if hdw and crw:
            def hdf(t):
                for a, b in hdw:
                    if a <= t < b:
                        return 1
                for a, b in crw:
                    if a <= t < b:
                        return 2
                return 0
            rx["hd"] = hdf
        elif hdw:
            rx["hd"] = hd_fn(hdw, 1)
        elif crw:
            rx["hd"] = hd_fn(crw, 2)

        # ---- fomepizole ----
        if self.fom:
            f = self.fom
            t0 = f["start"]
            load = f.get("bolus_mgkg", 15.0) * self.wt / MW["fom"]
            boluses.append((t0, IX["A_FOM1"], load))
            q = f.get("q", 12.0)
            maint = f.get("maint_mgkg", 10.0) * self.wt / MW["fom"]
            hq = f.get("hd_q", 4.0)
            tt = t0 + q
            n = 0
            while tt < self.tend:
                on_hd = any(a <= tt < b for a, b in hdw)
                dose = maint * (1.5 if n >= 4 else 1.0)   # 15 mg/kg after 4 doses
                boluses.append((tt, IX["A_FOM1"], dose))
                n += 1
                tt += (hq if (on_hd and self.fom_hd_boost) else q)
            # explicit q4h top-ups through dialysis when the schedule is honoured
            if self.fom_hd_boost:
                for a, b in hdw:
                    tt = a + hq
                    while tt < b:
                        boluses.append((tt, IX["A_FOM1"], maint))
                        tt += hq

        # ---- ethanol as antidote ----
        if self.eth:
            e = self.eth
            boluses.append((e["start"], IX["A_ETOH1"],
                            e.get("load_gkg", 0.7) * self.wt * 1000.0 / MW["etoh"]))
            rate = e.get("rate_mgkgh", 100.0) * self.wt / MW["etoh"]   # mmol/h
            hd_mult = e.get("hd_mult", 2.5)
            stop = e.get("stop", self.tend)
            st = e["start"]

            def ethf(t, st=st, stop=stop, rate=rate, hd_mult=hd_mult):
                if not (st <= t < stop):
                    return 0.0
                on = any(a <= t < b for a, b in hdw)
                return rate * (hd_mult if on else 1.0)
            rx["eth_inf"] = ethf

        # ---- bicarbonate ----
        if self.bicarb:
            b = self.bicarb
            rx["bicarb"] = step_fn([(b["start"], b.get("stop", self.tend))],
                                   b.get("rate_mmol_h", 25.0))
            if b.get("push_mmol", 0):
                boluses.append((b["start"], IX["HCO3"],
                                b["push_mmol"] / bicarb_space(24.0, self.wt)))
                boluses.append((b["start"], IX["CUM_BIC"], b["push_mmol"]))

        # ---- cofactors ----
        if self.folinate:
            f = self.folinate
            rx["fol"] = step_fn([(f["start"], f.get("stop", self.tend))],
                                f.get("rate", 0.20))
        if self.thiamine:
            f = self.thiamine
            rx["thi"] = step_fn([(f["start"], f.get("stop", self.tend))],
                                f.get("rate", 0.16))
        if self.pyridoxine:
            f = self.pyridoxine
            rx["pyr"] = step_fn([(f["start"], f.get("stop", self.tend))],
                                f.get("rate", 0.16))
        if self.calcium:
            c = self.calcium
            rx["ca_inf"] = step_fn([(c["start"], c.get("stop", self.tend))],
                                   c.get("rate_mmol_h", 2.0))
        return p, y, rx, sorted(boluses)


def simulate(sc, dt=0.05):
    p, y, rx, boluses = sc.build()
    # break points: bolus times, therapy edges
    brk = {0.0, sc.tend}
    for tb, _, _ in boluses:
        brk.add(tb)
    for a, b in list(sc.hd) + list(sc.crrt):
        brk.add(a); brk.add(b)
    for d in (sc.bicarb, sc.eth, sc.folinate, sc.thiamine, sc.pyridoxine, sc.calcium):
        if d:
            brk.add(d["start"])
            if d.get("stop"):
                brk.add(d["stop"])
    brk = sorted(x for x in brk if 0.0 <= x <= sc.tend)

    ts, ys = [0.0], [y.copy()]
    cur = y.copy()
    for i in range(len(brk) - 1):
        t0, t1 = brk[i], brk[i + 1]
        for tb, ix, amt in boluses:
            if abs(tb - t0) < 1e-9:
                cur[ix] += amt
        if t1 - t0 < 1e-9:
            continue
        n = max(2, int(math.ceil((t1 - t0) / dt)) + 1)
        teval = np.linspace(t0, t1, n)
        sol = solve_ivp(rhs, (t0, t1), cur, args=(p, rx), method="LSODA",
                        t_eval=teval, rtol=1e-7, atol=1e-9, max_step=0.25)
        if not sol.success:
            raise RuntimeError(f"{sc.name}: integration failed at {t0}-{t1}: {sol.message}")
        for k in range(1, sol.y.shape[1]):
            ts.append(sol.t[k]); ys.append(sol.y[:, k].copy())
        cur = sol.y[:, -1].copy()
    # final bolus at tend is irrelevant
    T = np.array(ts)
    Y = np.array(ys)
    O = [outputs(T[i], Y[i], p, rx) for i in range(len(T))]
    # A model that keeps integrating past pH 6.6 reports numbers (pH 5.5, CNS
    # formate 120 mM) that describe a corpse, not a patient.  Truncate at the
    # first non-survivable point and mark it.
    fatal = None
    for i, o in enumerate(O):
        if o["pH"] < p["PH_FATAL"]:
            fatal = T[i]
            break
    if fatal is not None:
        k = int(np.searchsorted(T, fatal, side="right"))
        T, Y, O = T[:k], Y[:k], O[:k]
        for o in O:
            o["FATAL_AT"] = fatal
    return p, T, Y, O


def series(O, key):
    return np.array([o[key] for o in O])


def at(T, O, key, tt):
    i = int(np.argmin(np.abs(T - tt)))
    return O[i][key]


def summarise(sc, T, O):
    ph = series(O, "pH")
    return dict(
        name=sc.name,
        t_end=float(T[-1]),
        fatal=float(O[-1].get("FATAL_AT", float("nan"))),
        pH_min=float(ph.min()),
        t_pH_min=float(T[int(np.argmin(ph))]),
        HCO3_min=float(series(O, "HCO3").min()),
        AG_max=float(series(O, "AG").max()),
        OG_max=float(series(O, "OG").max()),
        FORM_max=float(series(O, "FORM_mM").max()),
        FCNS_max=float(series(O, "FCNS_mM").max()),
        FVIT_max=float(series(O, "FVIT_mM").max()),
        GLYC_max=float(series(O, "GLYC_mM").max()),
        OXAL_max_uM=float(series(O, "OXAL_uM").max()),
        CA_min=float(series(O, "CA_ION").min()),
        GFR_min=float(series(O, "GFR_pct").min()),
        CAOXK=float(series(O, "CAOX_KID").max()),
        logMAR=float(series(O, "logMAR")[-1]),
        PUT=float(series(O, "PUT")[-1]),
        P_death=float(series(O, "P_DEATH")[-1]),
        P_blind=float(series(O, "P_BLIND")[-1]),
        OX_M=float(series(O, "CUM_OX_M")[-1]),
        OX_E=float(series(O, "CUM_OX_E")[-1]),
        HD_M=float(series(O, "CUM_HD_M")[-1]),
        HD_F=float(series(O, "CUM_HD_F")[-1]),
        HD_FOM=float(series(O, "CUM_HD_FOM")[-1]),
    )


# ---------------------------------------------------------------------------
# scenario library (mirrors the mrgsolve file)
# ---------------------------------------------------------------------------
def scenarios():
    S = {}
    # --------------------------- methanol ------------------------------
    S["M1_untreated"] = Scenario(
        "M1 · methanol 0.7 g/kg, untreated", meoh_g_per_kg=0.7, tend=72,
        notes="the natural history the antidote era abolished")
    S["M2_fom_early"] = Scenario(
        "M2 · methanol 0.7 g/kg, fomepizole at 2 h", meoh_g_per_kg=0.7, tend=96,
        fom=dict(start=2.0))
    S["M3_fom_late"] = Scenario(
        "M3 · methanol 0.7 g/kg, fomepizole at 14 h", meoh_g_per_kg=0.7, tend=96,
        fom=dict(start=14.0))
    S["M4_fom_hd"] = Scenario(
        "M4 · methanol 0.7 g/kg, fomepizole 8 h + HD 9-15 h (antidote re-dosed q4h)",
        meoh_g_per_kg=0.7, tend=96, fom=dict(start=8.0), hd=[(9.0, 15.0)],
        bicarb=dict(start=8.0, stop=26.0, rate_mmol_h=30.0, push_mmol=100.0))
    S["M5_hd_only"] = Scenario(
        "M5 · methanol 0.7 g/kg, HD alone at 15 h (no antidote)",
        meoh_g_per_kg=0.7, tend=72, hd=[(15.0, 21.0)])
    S["M6_coingest"] = Scenario(
        "M6 · methanol 0.7 g/kg + ethanol 0.8 g/kg co-ingested, untreated",
        meoh_g_per_kg=0.7, etoh_g_per_kg=0.8, tend=96)
    S["M7_ethanol_rx"] = Scenario(
        "M7 · methanol 0.7 g/kg, ethanol antidote from 4 h",
        meoh_g_per_kg=0.7, tend=96,
        eth=dict(start=4.0, load_gkg=0.7, rate_mgkgh=100.0, stop=60.0))
    S["M8_bicarb_only"] = Scenario(
        "M8 · methanol 0.7 g/kg, bicarbonate only from 8 h",
        meoh_g_per_kg=0.7, tend=72,
        bicarb=dict(start=8.0, stop=48.0, rate_mmol_h=35.0, push_mmol=100.0))
    S["M9_full"] = Scenario(
        "M9 · methanol 0.7 g/kg, full care at 6 h (fomepizole+bicarb+HD)",
        meoh_g_per_kg=0.7, tend=96, fom=dict(start=6.0), hd=[(7.0, 13.0)],
        bicarb=dict(start=6.0, stop=24.0, rate_mmol_h=30.0, push_mmol=100.0),
        folinate=dict(start=6.0, stop=48.0))
    S["M10_folinate"] = Scenario(
        "M10 · methanol 0.7 g/kg, fomepizole 14 h + folinic acid",
        meoh_g_per_kg=0.7, tend=96, fom=dict(start=14.0),
        folinate=dict(start=14.0, stop=60.0))
    S["M11_fom_hd_nodose"] = Scenario(
        "M11 · as M4 but fomepizole NOT re-dosed during dialysis",
        meoh_g_per_kg=0.7, tend=96, fom=dict(start=8.0), hd=[(9.0, 15.0)],
        bicarb=dict(start=8.0, stop=26.0, rate_mmol_h=30.0, push_mmol=100.0),
        fom_hd_boost=False)
    S["M12_massive"] = Scenario(
        "M12 · methanol 1.5 g/kg massive, fomepizole+HD at 4 h",
        meoh_g_per_kg=1.5, tend=120, fom=dict(start=4.0),
        hd=[(5.0, 13.0), (20.0, 26.0)],
        bicarb=dict(start=4.0, stop=30.0, rate_mmol_h=40.0, push_mmol=150.0))
    S["M13_crrt"] = Scenario(
        "M13 · as M4 but CRRT (40 mL/min) for 30 h instead of 6 h of IHD",
        meoh_g_per_kg=0.7, tend=120, fom=dict(start=8.0), crrt=[(9.0, 39.0)],
        bicarb=dict(start=8.0, stop=34.0, rate_mmol_h=30.0, push_mmol=100.0))
    # ------------------------ ethylene glycol ---------------------------
    S["M14_eth_hd_up"] = Scenario(
        "M14 · methanol 0.7 g/kg, ETHANOL antidote 8 h + HD 9-15 h, rate x2.5 on HD",
        meoh_g_per_kg=0.7, tend=96,
        eth=dict(start=8.0, load_gkg=0.7, rate_mgkgh=100.0, stop=60.0, hd_mult=2.5),
        hd=[(9.0, 15.0)],
        bicarb=dict(start=8.0, stop=26.0, rate_mmol_h=30.0, push_mmol=100.0))
    S["M15_eth_hd_flat"] = Scenario(
        "M15 · as M14 but the ethanol infusion rate is NOT increased on dialysis",
        meoh_g_per_kg=0.7, tend=96,
        eth=dict(start=8.0, load_gkg=0.7, rate_mgkgh=100.0, stop=60.0, hd_mult=1.0),
        hd=[(9.0, 15.0)],
        bicarb=dict(start=8.0, stop=26.0, rate_mmol_h=30.0, push_mmol=100.0))
    S["M16_crrt_nodose"] = Scenario(
        "M16 · as M13 (30 h CRRT) but fomepizole NOT re-dosed for the circuit",
        meoh_g_per_kg=0.7, tend=120, fom=dict(start=8.0), crrt=[(9.0, 39.0)],
        bicarb=dict(start=8.0, stop=34.0, rate_mmol_h=30.0, push_mmol=100.0),
        fom_hd_boost=False)
    S["M17_isopropanol"] = Scenario(
        "M17 · CONTROL: same ADH flux, but the metabolite is a KETONE",
        meoh_g_per_kg=0.7, tend=72, noacid=1.0)
    S["E1_untreated"] = Scenario(
        "E1 · ethylene glycol 90 mL, untreated", eg_mL=90.0, tend=96)
    S["E2_fom_early"] = Scenario(
        "E2 · EG 90 mL, fomepizole at 2 h, no dialysis", eg_mL=90.0, tend=96,
        fom=dict(start=2.0))
    S["E3_fom_late_hd"] = Scenario(
        "E3 · EG 90 mL, fomepizole 12 h + HD 13-19 h", eg_mL=90.0, tend=96,
        fom=dict(start=12.0), hd=[(13.0, 19.0)],
        bicarb=dict(start=12.0, stop=30.0, rate_mmol_h=30.0, push_mmol=100.0))
    S["E4_cofactors"] = Scenario(
        "E4 · EG 90 mL, thiamine+pyridoxine only from 2 h", eg_mL=90.0, tend=96,
        thiamine=dict(start=2.0, stop=60.0), pyridoxine=dict(start=2.0, stop=60.0))
    S["E5_fom_cofactors"] = Scenario(
        "E5 · EG 90 mL, fomepizole 2 h + thiamine/pyridoxine", eg_mL=90.0,
        tend=96, fom=dict(start=2.0), thiamine=dict(start=2.0, stop=60.0),
        pyridoxine=dict(start=2.0, stop=60.0))
    S["E6_ckd"] = Scenario(
        "E6 · EG 90 mL in CKD (GFR 30%), fomepizole at 2 h, no dialysis",
        eg_mL=90.0, tend=120, gfr_base=0.30, fom=dict(start=2.0))
    S["E7_ethanol_rx"] = Scenario(
        "E7 · EG 90 mL, ethanol antidote from 3 h + HD", eg_mL=90.0, tend=96,
        eth=dict(start=3.0, load_gkg=0.7, rate_mgkgh=100.0, stop=40.0),
        hd=[(6.0, 12.0)],
        bicarb=dict(start=3.0, stop=24.0, rate_mmol_h=30.0, push_mmol=100.0))
    S["E8_calcium"] = Scenario(
        "E8 · EG 90 mL, fomepizole 12 h + HD + empiric calcium loading",
        eg_mL=90.0, tend=96, fom=dict(start=12.0), hd=[(13.0, 19.0)],
        bicarb=dict(start=12.0, stop=30.0, rate_mmol_h=30.0, push_mmol=100.0),
        calcium=dict(start=12.0, stop=36.0, rate_mmol_h=3.0))
    return S


# ===========================================================================
#                          VERIFICATION SECTIONS
# ===========================================================================
def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


def sec1_kinetic_constants():
    hdr("1 · ADH BLOCKADE ARITHMETIC — how complete is each antidote, really?")
    p = default_params()
    print("Competing-substrate form:  v_i = Vmax_i*(S_i/Km_i) / (1 + SUM_j S_j/Km_j + F/Ki)")
    print("Km (mM): ethanol %.1f  <  ethylene glycol %.1f  <  methanol %.1f"
          % (p["KM_ET"], p["KM_E"], p["KM_M"]))
    print("Ki fomepizole = %.2g mM (%.2f uM)\n" % (p["KI_FOM"], p["KI_FOM"] * 1000))
    print("%-28s %10s %10s %10s %12s" %
          ("condition", "denom", "inh.factor", "flux left", "flux mmol/h"))
    for label, meoh_mgdl, etoh_mgdl, fom_ugml in [
            ("methanol 100, nothing", 100, 0, 0),
            ("methanol 100, ethanol 100", 100, 100, 0),
            ("methanol 100, ethanol 150", 100, 150, 0),
            ("methanol 100, fome 10 ug/mL", 100, 0, 10),
            ("methanol 100, fome 6 ug/mL", 100, 0, 6),
            ("methanol 20, nothing", 20, 0, 0),
            ("methanol 20, ethanol 100", 20, 100, 0),
            ("methanol 20, fome 10 ug/mL", 20, 0, 10),
            ("methanol 400, ethanol 100", 400, 100, 0),
            ("methanol 400, fome 10 ug/mL", 400, 0, 10)]:
        M = mgdl_to_mM(meoh_mgdl, "meoh")
        E = mgdl_to_mM(etoh_mgdl, "etoh")
        F = fom_ugml / MW["fom"] * 1.0          # ug/mL -> umol/mL = mM
        sM = M / p["KM_M"]
        den = 1 + sM + E / p["KM_ET"] + F / p["KI_FOM"]
        base = 1 + sM
        v = p["VMAX_M"] * sM / den
        print("%-28s %10.1f %10.1f %9.3f%% %12.2f"
              % (label, den, den / base, 100 * base / den, v))
    print("""
READ THIS TABLE, NOT THE SLOGAN.  A competitive block is diluted by its own
substrate, so 'fomepizole is 1000x ethanol' is false at the bedside.  At a
methanol of 100 mg/dL, an ethanol infusion held at target still leaves ~18% of
the acid-forming flux running; fomepizole at 10 ug/mL leaves 0.6%.  The honest
ratio is ~30x, not ~1000x — and at a methanol of 400 mg/dL even fomepizole
leaves 2.4%, which is the arithmetic behind 'blockade is never absolute'.""")


def sec2_ion_trapping():
    hdr("2 · ION TRAPPING — why pH, not concentration, is the prognosticator")
    p = default_params()
    print("f_HA(pH) = 1/(1 + 10^(pH - 3.75))   [fraction of formate that can cross a membrane]\n")
    print("%8s %14s %10s %14s %10s %12s"
          % ("pH", "f_HA plasma", "vs 7.40", "f_HA brain", "vs 7.40", "eq CNS/plasma"))
    ref_p = f_HA(7.40, p["PKA_F"])
    ref_b = f_HA(7.00, p["PKA_F"])
    for ph in [7.45, 7.40, 7.30, 7.20, 7.10, 7.00, 6.90, 6.80, 6.70]:
        phb = 7.00 + p["PHB_SLOPE"] * (ph - 7.40)
        fp, fb = f_HA(ph, p["PKA_F"]), f_HA(phb, p["PKA_F"])
        print("%8.2f %14.3e %10.2f %14.3e %10.2f %12.3f"
              % (ph, fp, fp / ref_p, fb, fb / ref_b, fp / fb))
    print("""
TWO SEPARATE EFFECTS, SAME EQUATION, SAME SIGN.  Falling from pH 7.40 to 6.80
multiplies the ENTRY RATE CONSTANT into the brain by 3.98 and multiplies the
EQUILIBRIUM CNS/plasma formate ratio by 1.78.  At identical plasma formate the
CNS burden therefore roughly doubles-to-triples across a range of pH that
clinicians routinely treat as a number to be corrected rather than a mechanism.
Bicarbonate runs both terms backwards, and simultaneously raises the urinary
f_ion so renal formate clearance rises: one intervention, three terms.""")


def sec3_bicarb_and_urine():
    hdr("3 · BICARBONATE — the same equation seen from the kidney")
    p = default_params()
    print("Tubular reabsorption is non-ionic diffusion of the un-ionised acid, so")
    print("the reabsorbed fraction follows f_HA(URINE pH).\n")
    print("%10s %10s %12s %14s %14s"
          % ("HCO3", "urine pH", "f_HA urine", "reabsorbed", "renal CL (L/h)"))
    for h in [4, 8, 12, 16, 20, 24, 28, 32]:
        pu = urine_pH(h)
        fu = f_HA(pu, p["PKA_F"])
        fr = p["FR_MAX"] * fu / (fu + p["KFR"])
        print("%10.0f %10.2f %12.3e %13.1f%% %14.2f"
              % (h, pu, fu, 100 * fr, p["GFR0"] * (1 - fr)))
    print("""
The clinically reachable span is a 4-to-5-fold change in renal formate
clearance.  Note the direction of the trap: the SICKER the patient (low HCO3,
acid urine) the MORE formate the kidney gives back — a second self-amplifying
loop, and the reason 'the acidosis is only a number' is wrong.""")


def sec4_bicarb_space():
    hdr("4 · WHY PLASMA BICARBONATE DECELERATES ON ITS WAY DOWN")
    print("Apparent bicarbonate space = (0.40 + 2.6/[HCO3]) L/kg  (70 kg)\n")
    print("%10s %16s %22s" % ("HCO3 (mM)", "space (L)", "space / Vd(methanol) 42 L"))
    for h in [24, 20, 16, 12, 10, 8, 6, 4]:
        v = bicarb_space(h, 70.0)
        print("%10.0f %16.1f %22.3f" % (h, v, v / 42.0))
    print("""
This is the quantitative content of the two-gap invariant.  The osmolal gap
falls by (converted mmol)/Vd_methanol = /42 L.  The anion gap rises by
(titrated mmol)/bicarbonate space.  At a normal HCO3 the space is 35.6 L, so
the AG rises ~18% FASTER than the OG falls and the SUM drifts upward.  As the
acidosis deepens the space expands past 42 L, the AG rise decelerates, and the
sum comes back down.  The invariant is therefore good to roughly +/-20% over
the clinically relevant range — which is stated, not assumed, and is measured
in section 5.""")


def sec5_two_gap_clock():
    hdr("5 · THE TWO-GAP CLOCK, MEASURED IN THE SIMULATION")
    S = scenarios()
    sc = S["M1_untreated"]
    p, T, Y, O = simulate(sc)
    dose = sc.meoh_g_per_kg * sc.wt * 1000.0 / MW["meoh"]
    print("methanol %.1f g/kg = %.0f mmol ingested; untreated\n" % (sc.meoh_g_per_kg, dose))
    print("%6s %10s %8s %8s %8s %9s %9s %8s"
          % ("t (h)", "MeOH", "formate", "OG", "AG-12", "OG/(AG-12)", "gap sum", "pH"))
    print("%6s %10s %8s %8s %8s %9s %9s %8s"
          % ("", "mg/dL", "mM", "mOsm", "mEq", "ratio", "mOsm+mEq", ""))
    for tt in [0.5, 1, 2, 4, 6, 8, 12, 16, 20, 24, 30, 36, 48]:
        i = int(np.argmin(np.abs(T - tt)))
        o = O[i]
        print("%6.1f %10.1f %8.2f %8.1f %8.1f %9.2f %9.1f %8.3f"
              % (tt, o["MEOH_mgdL"], o["FORM_mM"], o["OG"], o["AG"] - 12,
                 o["OG_over_AG"], o["GAPSUM"], o["pH"]))
    gs = series(O, "GAPSUM")
    early = gs[np.argmin(np.abs(T - 1.0))]
    print("\ngap sum at 1 h  = %.1f      (dose/Vd_methanol = %.1f)"
          % (early, dose / p["V1M"] / (1 + p["V2M"] / p["V1M"])))
    j = int(np.argmax(gs))
    print("gap sum peak    = %.1f at t = %.1f h   (+%.0f%% vs 1 h)"
          % (gs[j], T[j], 100 * (gs[j] / early - 1)))
    print("gap sum at 48 h = %.1f            (%.0f%% of the 1 h value)"
          % (gs[-1], 100 * gs[-1] / early))
    print("""
THREE THINGS THE TABLE SAYS, ONE OF WHICH CONTRADICTS THE TEACHING.

(1) THE CLOCK WORKS, over the interval that matters.  OG/(AG-12) falls from 37
    at half an hour to 0.65 by 24 h - four orders of magnitude of ordering
    information about TIME SINCE INGESTION, with no reference to severity.  It
    is NOT monotone forever: once the reaction is over, formate is cleared, the
    anion gap closes, and the ratio turns back up (0.65 -> 1.26 by 48 h).  So
    the ratio dates the ingestion only on the way down, and a ratio near 1 is
    ambiguous between "half-way through" and "finished and recovering".  The
    disambiguator is the SUM, which is still large in the first case.

(2) THE SUM IS A FLOOR ON THE DOSE, NOT AN ESTIMATE OF IT.  At 1 h it reads
    35.0 against a true 36.4 mmol/L of ingested methanol - 96% recovery.  It
    peaks 4% above that (the bicarbonate-space effect of section 4) and then
    decays to 31% of its initial value by 48 h as formate is oxidised and
    excreted.  Read early it is quantitative; read late it under-reads.

(3) THE OSMOLAL GAP DOES NOT RETURN TO ZERO, AND THIS IS DERIVED.  It bottoms
    at 12.2 mOsm with the methanol essentially all gone.  Textbooks say the
    osmolal gap normalises as the anion gap widens because formate replaces
    bicarbonate one-for-one - but that is only true if the two share a
    distribution volume, and they do not.  Formate spreads through 35 L while
    the bicarbonate it titrates is buffered over 45-73 L, so plasma formate
    rises MORE than plasma bicarbonate falls (this is the same fact as the
    delta-ratio >1 of organic acidosis) and the difference shows up as a
    residual osmolal gap.  A methanol patient 24 h out is therefore expected to
    have a MODERATELY raised osmolal gap and a wide anion gap simultaneously,
    and a clinician who takes the persisting gap as evidence of unabsorbed
    alcohol will reach for the wrong therapy.""")


def sec6_natural_history():
    hdr("6 · SCENARIO LIBRARY — every run, one line each")
    S = scenarios()
    rows = []
    for k, sc in S.items():
        p, T, Y, O = simulate(sc)
        rows.append((k, summarise(sc, T, O)))
    print("%-20s %6s %5s %5s %6s %6s %6s %5s %5s %5s %5s %5s %5s %5s"
          % ("scenario", "pHmin", "HCO3", "AG", "form", "FCNS", "glyc", "CaOx",
             "iCa", "GFR%", "PUT", "logMR", "P(d)", "P(b)"))
    for k, s in rows:
        print("%-20s %6.3f %5.1f %5.0f %6.2f %6.2f %6.2f %5.1f %5.2f %5.0f %5.2f %5.2f %5.2f %5.2f%s"
              % (k, s["pH_min"], s["HCO3_min"], s["AG_max"], s["FORM_max"],
                 s["FCNS_max"], s["GLYC_max"], s["CAOXK"], s["CA_min"],
                 s["GFR_min"], s["PUT"], s["logMAR"], s["P_death"], s["P_blind"],
                 ("   DIED %.0fh" % s["fatal"]) if s["fatal"] == s["fatal"] else ""))
    return dict(rows)


def sec7_fomepizole_timing():
    hdr("7 · THE COST OF DELAY — fomepizole start time, swept")
    print("methanol 0.7 g/kg, fomepizole only (no dialysis, no bicarbonate)\n")
    print("%8s %8s %8s %8s %8s %8s %8s"
          % ("start h", "pH min", "form max", "FCNS max", "logMAR", "P(die)", "%dose ox"))
    dose = 0.7 * 70 * 1000 / MW["meoh"]
    for t0 in [0, 1, 2, 4, 6, 8, 10, 12, 14, 18, 24, 30]:
        sc = Scenario("x", meoh_g_per_kg=0.7, tend=96, fom=dict(start=float(t0)))
        p, T, Y, O = simulate(sc)
        s = summarise(sc, T, O)
        print("%8.0f %8.3f %8.2f %8.2f %8.2f %8.2f %8.0f"
              % (t0, s["pH_min"], s["FORM_max"], s["FCNS_max"], s["logMAR"],
                 s["P_death"], 100 * s["OX_M"] / dose))
    print("""
The 'antidote window' is not a rule, it is the shape of this column.  Because
the parent alcohol is consumed by a saturable enzyme, the fraction already
oxidised grows almost linearly for the first several hours and then plateaus.
Two features are worth naming.  First, there is no cliff: the outcome columns
degrade smoothly, so 'the antidote window' is a marketing phrase for a slope.
Second, the LAST column saturates before the outcome columns do -- between 18
and 30 h the fraction of the dose already oxidised barely moves (79 -> 84%), and
yet the difference between those two presentations is the difference between a
survivor and a death.  What is still changing after the reaction has finished is
not how much acid was made but how long it stays, and nothing about fomepizole
addresses that.  This is the argument for dialysing the late presenter, and it
is why 'the level was only 40 mg/dL, so we held the antidote' is the wrong
sentence: by then the level is a measure of what is NOT going to hurt them.""")


def sec8_eg_vs_meoh_blocked():
    hdr("8 · SAME ANTIDOTE, OPPOSITE MANAGEMENT — the blocked-state half-life")
    print("Both alcohols, ADH fully blocked from t = 2 h, no dialysis.\n")
    for label, kw in [("methanol 0.7 g/kg", dict(meoh_g_per_kg=0.7)),
                      ("ethylene glycol 90 mL", dict(eg_mL=90.0)),
                      ("EG 90 mL, GFR 30%", dict(eg_mL=90.0, gfr_base=0.30))]:
        sc = Scenario("x", tend=160, fom=dict(start=2.0), **kw)
        p, T, Y, O = simulate(sc)
        key = "MEOH_mgdL" if kw.get("meoh_g_per_kg") else "EG_mgdL"
        c = series(O, key)
        i0 = int(np.argmax(c))
        # log-linear slope over the blocked decay phase
        m = (T > T[i0] + 4) & (c > 5)
        if m.sum() > 10:
            sl = np.polyfit(T[m], np.log(c[m]), 1)[0]
            thalf = -math.log(2) / sl
        else:
            thalf = float("nan")
        i20 = np.where((c < 20) & (T > T[i0]))[0]
        t20 = T[i20[0]] if len(i20) else float("nan")
        print("%-24s peak %6.1f mg/dL   blocked t1/2 = %5.1f h   "
              "<20 mg/dL at t = %5.1f h" % (label, c[i0], thalf, t20))
    print("""
This is statement (3) of the thesis as a number.  ADH blockade converts the
poison into a reservoir; whether that is a cure or a holding action is decided
by ONE clearance term.  Methanol has no renal exit, so a blocked patient needs
days of antidote or a dialysis machine.  Ethylene glycol leaks out through the
glomerulus, so blockade alone finishes the job — UNLESS the poison has already
taken the kidney, at which point EG starts behaving like methanol.  The third
row is the whole argument for dialysing the EG patient who is already in AKI.""")


def sec9_hd_antidote_loop():
    hdr("9 · THE LOOP THAT ISN'T — dialysis removes the antidote, but which one?")
    print("""The map draws haemodialysis removing fomepizole as a self-defeating
loop.  The model was built to display that loop.  It does not: the loop is
real in the equations and negligible in the patient, and the reason is one
ratio.  What matters is not how fast dialysis clears the antidote but how much
MARGIN there is between the therapeutic concentration and the concentration at
which the block fails -- i.e. C_therapeutic / K_inhibition.
""")
    p = default_params()
    fom_ther = 10.0 / MW["fom"]                       # 10 ug/mL in mM
    eth_ther = mgdl_to_mM(100.0, "etoh")              # 100 mg/dL in mM
    print("%-14s %14s %14s %14s" % ("antidote", "therapeutic", "K_inh (mM)", "margin C/K"))
    print("%-14s %14.4g %14.4g %14.0f"
          % ("fomepizole", fom_ther, p["KI_FOM"], fom_ther / p["KI_FOM"]))
    print("%-14s %14.4g %14.4g %14.0f"
          % ("ethanol", eth_ther, p["KM_ET"], eth_ther / p["KM_ET"]))
    kel_fom = p["CL_HD_FOM"] / p["V1FOM"]
    kel_eth = p["CL_HD_ET"] / p["V1ET"]
    print("\ndialytic rate constants: fomepizole %.3f /h (t1/2 %.1f h), "
          "ethanol %.3f /h (t1/2 %.1f h)"
          % (kel_fom, math.log(2) / kel_fom, kel_eth, math.log(2) / kel_eth))
    print("hours of uninterrupted dialysis needed to spend the whole margin:")
    print("   fomepizole  %5.1f h" % (math.log(fom_ther / p["KI_FOM"]) / kel_fom))
    print("   ethanol     %5.1f h" % (math.log(eth_ther / p["KM_ET"]) / kel_eth))
    print("""
A conventional 4-6 h session cannot spend fomepizole's margin; it can spend
most of ethanol's.  So the instruction 'increase the antidote during dialysis'
is a comfortable safety margin for fomepizole and a quantitative necessity for
ethanol -- the opposite of how the two are usually taught.  Below, the model is
asked to confirm it.
""")
    S = scenarios()
    pairs = [("fomepizole, 6 h IHD", "M4_fom_hd", "M11_fom_hd_nodose"),
             ("ethanol, 6 h IHD", "M14_eth_hd_up", "M15_eth_hd_flat"),
             ("fomepizole, 30 h CRRT", "M13_crrt", "M16_crrt_nodose")]
    print("Metrics that can SEE inside the session: total methanol oxidised over")
    print("the whole admission, and the CNS formate exposure integral.  (pH nadir")
    print("and peak formate are blind here - both occur before dialysis starts.)\n")
    print("%-24s %-14s %10s %10s %9s %8s"
          % ("comparison", "arm", "MeOH ox", "AUC FCNS", "logMAR", "P(die)"))
    print("%-24s %-14s %10s %10s %9s %8s"
          % ("", "", "mmol", "mM*h", "", ""))
    store = {}
    for lab, ka, kb in pairs:
        base = None
        for tag, k in (("re-dosed", ka), ("NOT re-dosed", kb)):
            sc = S[k]
            _, T, Y, O = simulate(sc)
            sm = summarise(sc, T, O)
            store[k] = (T, O, sm)
            auc = float(Y[-1][IX["AUC_FCNS"]])
            if base is None:
                base = (sm["OX_M"], auc)
                extra = ""
            else:
                extra = ("   (%+.1f%% oxidised, %+.1f%% exposure)"
                         % (100 * (sm["OX_M"] / base[0] - 1),
                            100 * (auc / base[1] - 1)))
            print("%-24s %-14s %10.1f %10.1f %9.3f %8.3f%s"
                  % (lab if tag == "re-dosed" else "", tag, sm["OX_M"], auc,
                     sm["logMAR"], sm["P_death"], extra))
    print("\nInside the session — inhibitor concentration and the flux it leaves:")
    print("%6s | %10s %8s %8s | %8s %8s %8s %8s"
          % ("t (h)", "fome", "F/Ki", "flux", "EtOH up", "flux", "EtOH flat", "flux"))
    Ta, Oa, _ = store["M11_fom_hd_nodose"]
    Tb, Ob, _ = store["M14_eth_hd_up"]
    Tc, Oc, _ = store["M15_eth_hd_flat"]
    for tt in [8.0, 9.0, 10, 11, 12, 13, 14, 14.9, 16, 18, 24]:
        ia = int(np.argmin(np.abs(Ta - tt)))
        ib = int(np.argmin(np.abs(Tb - tt)))
        ic = int(np.argmin(np.abs(Tc - tt)))
        print("%6.1f | %10.3f %8.0f %8.3f | %8.1f %8.3f %8.1f %8.3f"
              % (tt, Oa[ia]["FOM_ugmL"],
                 Oa[ia]["FOM_ugmL"] / MW["fom"] / 1.5e-4, Oa[ia]["FLUX_M"],
                 Ob[ib]["ETOH_mgdL"], Ob[ib]["FLUX_M"],
                 Oc[ic]["ETOH_mgdL"], Oc[ic]["FLUX_M"]))
    print("""
CONCLUSION, AGAINST THE AUTHOR'S PRIOR — THREE PARTS.

(a) The fomepizole arm is flat.  Six hours of dialysis with no top-up leaves the
    plasma concentration at 173x its Ki; the flux moves from 0.132 to 0.174
    mmol/h against an unblocked 67, and every outcome column is identical.  The
    q4h instruction is a comfortable margin, not a knife-edge.

(b) The ethanol arm is NOT flat, and the instantaneous flux is where you see it:
    at the end of the session the un-boosted arm is oxidising 3.74 mmol/h
    against 1.36 mmol/h in the boosted arm, a 2.8-fold difference, because the
    ethanol has fallen 121 -> 27 mg/dL.  Ethanol's therapeutic margin is 22, not
    812, and a single session spends most of it.

(c) The aggregate effect is nevertheless small (+2.4% of the dose oxidised),
    and the reason is worth stating because it was not designed in: DIALYSIS
    PARTLY PROTECTS AGAINST ITS OWN ANTIDOTE REMOVAL, because it is removing the
    substrate at the same time.  The two terms in the same session pull in
    opposite directions on the same numerator.

    Second-order and honestly reported: P(death) is marginally LOWER in the
    un-boosted ethanol arm (0.143 vs 0.155), because a 2.5x infusion rate buys
    a small amount of blockade at the cost of a larger amount of sedation, and
    sedation sits in the respiratory-failure hazard.  The model is not asserting
    that under-dosing ethanol is good; it is showing that ethanol's benefit and
    harm live in the same dose, which fomepizole's do not.

Cluster 17 of the map keeps the loop drawn, with this caveat attached rather
than quietly redrawn.""")


def sec10_crystal_threshold():
    hdr("10 · THE ONLY THRESHOLD IN THE MODEL — CaOx supersaturation")
    p = default_params()
    print("SS = [Ca2+][oxalate]/Ksp,  Ksp(CaOx monohydrate) = 2.32e-9 M^2")
    print("tubular fluid concentrates the filtrate ~%.0fx\n" % p["CFTUB"])
    print("%12s %10s %12s %14s %14s"
          % ("oxalate uM", "Ca2+ mM", "SS plasma", "SS tubular", "deposition?"))
    for ox in [2, 5, 10, 20, 40, 80, 160]:
        ssp = 1.10 * (ox / 1000.0) / p["KSP"]
        sst = ssp * p["CFTUB"]
        tag = ("tubule + tissue" if ssp > p["SSCRIT_T"] else
               ("tubule only" if sst > p["SSCRIT_K"] else "none"))
        print("%12.0f %10.2f %12.2f %14.1f %14s" % (ox, 1.10, ssp, sst, tag))
    print("\nRATE-DEPENDENCE, MEASURED: same total dose, different flux.")
    print("%-40s %10s %10s %10s %10s"
          % ("", "peak ox", "CaOx kid", "GFR min", "iCa min"))
    for lab, kw in [("EG 90 mL untreated", {}),
                    ("EG 90 mL, fomepizole at 2 h", dict(fom=dict(start=2.0))),
                    ("EG 90 mL, fomepizole at 6 h", dict(fom=dict(start=6.0))),
                    ("EG 90 mL, fomepizole at 12 h", dict(fom=dict(start=12.0)))]:
        sc = Scenario("x", eg_mL=90.0, tend=120, **kw)
        p2, T, Y, O = simulate(sc)
        s = summarise(sc, T, O)
        print("%-40s %10.1f %10.2f %10.0f %10.3f"
              % (lab, s["OXAL_max_uM"], s["CAOXK"], s["GFR_min"], s["CA_min"]))
    print("""
There is no equivalent table for methanol, and that asymmetry is the point.
Formate injures in proportion to its integral, so an equal dose delivered
slowly does equal harm.  Oxalate injures only above a solubility product, so an
equal dose delivered slowly does LESS harm — the same total oxalate simply
leaves in the urine.  This is why 'fomepizole bought time' is a weaker
statement for methanol than for ethylene glycol.""")


def sec11_cofactor_nonadditivity():
    hdr("11 · A PREDICTED NON-ADDITIVITY — thiamine/pyridoxine vs fomepizole")
    rows = []
    for lab, kw in [
            ("nothing", {}),
            ("thiamine + pyridoxine", dict(thiamine=dict(start=2.0, stop=72.0),
                                           pyridoxine=dict(start=2.0, stop=72.0))),
            ("fomepizole 2 h", dict(fom=dict(start=2.0))),
            ("fomepizole 2 h + both cofactors",
             dict(fom=dict(start=2.0), thiamine=dict(start=2.0, stop=72.0),
                  pyridoxine=dict(start=2.0, stop=72.0))),
            ("fomepizole 12 h", dict(fom=dict(start=12.0))),
            ("fomepizole 12 h + both cofactors",
             dict(fom=dict(start=12.0), thiamine=dict(start=2.0, stop=72.0),
                  pyridoxine=dict(start=2.0, stop=72.0)))]:
        sc = Scenario("x", eg_mL=90.0, tend=120, **kw)
        p, T, Y, O = simulate(sc)
        rows.append((lab, summarise(sc, T, O)))
    print("%-36s %10s %10s %10s %8s"
          % ("EG 90 mL +", "peak ox uM", "CaOx mmol", "GFR min %", "pH min"))
    for lab, s in rows:
        print("%-36s %10.1f %10.2f %10.0f %8.3f"
              % (lab, s["OXAL_max_uM"], s["CAOXK"], s["GFR_min"], s["pH_min"]))
    b = {l: s for l, s in rows}
    print("\n%-28s %14s %14s" % ("cofactor benefit on CaOx", "absolute mmol", "relative %"))
    for a, c, lab in [("nothing", "thiamine + pyridoxine", "alone"),
                      ("fomepizole 2 h", "fomepizole 2 h + both cofactors",
                       "on top of fome 2 h"),
                      ("fomepizole 12 h", "fomepizole 12 h + both cofactors",
                       "on top of fome 12 h")]:
        d = b[a]["CAOXK"] - b[c]["CAOXK"]
        print("%-28s %14.2f %13.1f%%"
              % (lab, d, 100 * d / max(b[a]["CAOXK"], 1e-9)))
    print("""
THE PREDICTION SURVIVES IN ABSOLUTE TERMS AND FAILS IN RELATIVE TERMS, AND THE
DIFFERENCE MATTERS.  Thiamine and pyridoxine widen two escape routes out of the
glyoxylate node, so what they own is a RATIO applied to a flux.  The RELATIVE
benefit is therefore roughly scale-free -- about 36-52% off the crystal burden
whatever else is running -- which is why the percentage column looks flat and
why a percentage is the wrong way to read this.  The ABSOLUTE benefit tracks the
residual flux exactly: 16 mmol of crystal prevented with no antidote, 15 mmol
alongside a LATE antidote, and 2.4 mmol alongside a prompt one.

So the corrected statement is narrower and more useful than the one the map
implies.  The cofactors are not redundant to fomepizole in general; they are
redundant to PROMPT fomepizole, because prompt fomepizole is the only thing that
leaves no glyoxylate to redirect.  In the late presenter -- the patient for whom
the antidote is already mostly beside the point -- two cheap and harmless drugs
retain essentially their full standalone effect.  A trial that randomises
thiamine and pyridoxine in patients who all receive fomepizole within two hours
is designed to find nothing; the same trial in late presenters is not.

Note also what the cofactors do NOT fix: GFR moves only 26 -> 29% in the
untreated arm, because most of that injury is glycolaldehyde acting directly on
the tubule, upstream of the branch point they operate on.""")


def sec12_isopropanol_control():
    hdr("12 · THE CONTROL — an alcohol whose metabolite is not an acid")
    print("""Isopropanol is not simulated as a separate compartment.  It is
simulated by DELETING ONE TERM: scenario M17 runs the identical ADH flux but the
metabolite is a ketone (acetone), so no proton and no anion appear.  Everything
downstream of the acid must vanish while the OSMOLAL GAP survives untouched --
which is the falsifiable content of "the metabolite is the poison".
""")
    S = scenarios()
    a, b = S["M1_untreated"], S["M17_isopropanol"]
    _, T1, Y1, O1 = simulate(a)
    _, T2, Y2, O2 = simulate(b)
    s1, s2 = summarise(a, T1, O1), summarise(b, T2, O2)
    print("%-34s %14s %16s" % ("", "acid-forming", "ketone-forming"))
    rows = [("peak osmolal gap (mOsm)", "OG_max", 1),
            ("peak anion gap (mEq)", "AG_max", 1),
            ("pH nadir", "pH_min", 3),
            ("bicarbonate nadir (mM)", "HCO3_min", 1),
            ("peak plasma formate (mM)", "FORM_max", 2),
            ("peak CNS formate (mM)", "FCNS_max", 2),
            ("putaminal injury index", "PUT", 3),
            ("final logMAR", "logMAR", 2),
            ("P(death)", "P_death", 3),
            ("methanol oxidised (mmol)", "OX_M", 0)]
    for lab, k, dg in rows:
        print(("%-34s %14." + str(dg) + "f %16." + str(dg) + "f")
              % (lab, s1[k], s2[k]))
    print("""
Same alcohol, same ADH flux (%.0f vs %.0f mmol oxidised), same osmolal gap, one
term deleted -- and the disease disappears.  That is isopropanol poisoning: a
large osmolal gap, a drunk patient, ketones, and a NORMAL anion gap, managed
with a bed and time.  Any model of methanol that does not do this when the acid
term is removed is modelling the wrong molecule.

The second control lives in cluster 20 of the map and needs no code: propylene
glycol raises BOTH gaps, but the acid is ordinary L-lactate and the source is
the patient's own sedative infusion.  The gaps are not specific.""" %
          (s1["OX_M"], s2["OX_M"]))


def sec13_prognostic_ranking():
    hdr("13 · pH vs CONCENTRATION AS A PREDICTOR — a derived, not fitted, result")
    print("""Sixteen untreated/partly-treated methanol runs spanning dose and
presentation.  For each, record admission-equivalent methanol, admission
formate, admission pH, and the final outcome.  Then ask which admission
variable orders the outcome better.  Nothing in the model was fitted to make
pH win; pH wins because it sits inside the delivery term.
""")
    rows = []
    for dose in [0.4, 0.6, 0.8, 1.0]:
        for t_adm in [6.0, 12.0, 18.0, 24.0]:
            sc = Scenario("x", meoh_g_per_kg=dose, tend=96,
                          fom=dict(start=t_adm), bicarb=None)
            p, T, Y, O = simulate(sc)
            i = int(np.argmin(np.abs(T - t_adm)))
            s = summarise(sc, T, O)
            rows.append(dict(dose=dose, t=t_adm, meoh=O[i]["MEOH_mgdL"],
                             form=O[i]["FORM_mM"], pH=O[i]["pH"],
                             out=s["P_death"], logmar=s["logMAR"],
                             died=(s["fatal"] == s["fatal"])))
    print("%6s %6s %9s %9s %8s %9s %8s %6s"
          % ("g/kg", "t_adm", "MeOH", "formate", "pH", "P(die)", "logMAR", "died"))
    for r in rows:
        print("%6.2f %6.0f %9.1f %9.2f %8.3f %9.3f %8.2f %6s"
              % (r["dose"], r["t"], r["meoh"], r["form"], r["pH"], r["out"],
                 r["logmar"], "yes" if r["died"] else ""))
    print("\n(Rows flagged 'died' were truncated at the first non-survivable pH,")
    print(" so their admission values repeat and their injury integrals stop early;")
    print(" the correlations are therefore reported twice, with and without them.)")

    def spearman(a, b):
        ra = np.argsort(np.argsort(a)).astype(float)
        rb = np.argsort(np.argsort(b)).astype(float)
        return float(np.corrcoef(ra, rb)[0, 1])
    for subset, tag in [(rows, "all 16 runs"),
                        ([r for r in rows if not r["died"]], "survivable runs only")]:
        out = np.array([r["out"] for r in subset])
        print("\nSpearman rho with P(death) — %s (n=%d):" % (tag, len(subset)))
        for k, lab in [("meoh", "admission methanol (mg/dL)"),
                       ("form", "admission formate (mM)"),
                       ("pH", "admission pH")]:
            v = np.array([r[k] for r in subset])
            print("   %-32s rho = %+.3f" % (lab, spearman(v, out)))
    print("""
THE CLINICAL CLAIM IS CONFIRMED; THE AUTHOR'S STRONGER VERSION IS NOT.

Confirmed, and strongly: the ADMISSION ALCOHOL LEVEL is useless and runs the
WRONG WAY (rho = -0.13).  That is not noise, it is structure -- a high methanol
means the reaction has not happened yet, which is exactly the patient the
antidote saves.  Any triage rule that escalates on the alcohol level alone is
sorting patients by how salvageable they are and treating that as severity.

Not confirmed: 'pH beats everything'.  Plasma FORMATE is the better single
predictor here (rho = +0.94 vs pH's -0.84), and it should be, because formate is
the proximate cause and pH is one step downstream of it.  The reason pH rules
the bedside is not that it is the better variable; it is that formate cannot be
measured in the four hours during which the decision has to be made.  pH is a
real-time surrogate for formate that costs nothing and arrives in ninety
seconds, and its residual advantage over formate is confined to the delivery
term -- two patients with the same formate are not equally poisoned, and the
more acidotic one is worse.  Stated that carefully, the claim is defensible;
stated as 'pH is the best predictor', it is an artefact of assay availability.""")


def sec14_conservation_checks():
    hdr("14 · CONSERVATION AND SANITY CHECKS (the model must not leak)")
    S = scenarios()
    bad = 0
    print("%-24s %14s %14s %10s %8s"
          % ("scenario", "MeOH in (mmol)", "accounted", "err %", "neg?"))
    for k in ["M1_untreated", "M4_fom_hd", "M9_full", "M12_massive"]:
        sc = S[k]
        p, T, Y, O = simulate(sc)
        din = sc.meoh_g_per_kg * sc.wt * 1000.0 / MW["meoh"]
        y = Y[-1]
        remaining = (y[IX["A_GUTM"]] + y[IX["A_MEOH1"]] + y[IX["A_MEOH2"]])
        oxidised = y[IX["CUM_OX_M"]]
        dial = y[IX["CUM_HD_M"]]
        # renal + pulmonary loss is the residual; recompute by integration
        acc = remaining + oxidised + dial
        # renal/pulm methanol loss
        cM = np.array([o["MEOH_mgdL"] for o in O]) * 10 / MW["meoh"]
        gf = np.array([o["GFR_pct"] for o in O]) / 100.0
        loss = np.trapezoid((0.050 * p["GFR0"] * gf + 0.15) * cM, T)
        acc += loss
        err = 100 * (acc - din) / din
        neg = bool((Y < -1e-6).any())
        bad += abs(err) > 1.0 or neg
        print("%-24s %14.1f %14.1f %10.3f %8s"
              % (k, din, acc, err, "YES" if neg else "no"))
    print("\nEG mass balance:")
    print("%-24s %14s %14s %10s" % ("scenario", "EG in (mmol)", "accounted", "err %"))
    for k in ["E1_untreated", "E3_fom_late_hd", "E6_ckd"]:
        sc = S[k]
        p, T, Y, O = simulate(sc)
        din = sc.eg_mL * 1.1132 * 1000.0 / MW["eg"]
        y = Y[-1]
        cE = np.array([o["EG_mgdL"] for o in O]) * 10 / MW["eg"]
        gf = np.array([o["GFR_pct"] for o in O]) / 100.0
        loss = np.trapezoid(0.26 * p["GFR0"] * gf * cE, T)
        acc = (y[IX["A_GUTE"]] + y[IX["A_EG1"]] + y[IX["A_EG2"]]
               + y[IX["CUM_OX_E"]] + y[IX["CUM_HD_E"]] + loss)
        err = 100 * (acc - din) / din
        bad += abs(err) > 1.0
        print("%-24s %14.1f %14.1f %10.3f" % (k, din, acc, err))
    print("\n%s" % ("ALL CHECKS PASS" if bad == 0 else "*** %d CHECK(S) FAILED ***" % bad))
    return bad


def sec15_reference_trajectory():
    hdr("15 · REFERENCE TRAJECTORY FOR mrgsolve CROSS-CHECK (M4)")
    sc = scenarios()["M4_fom_hd"]
    p, T, Y, O = simulate(sc)
    keys = ["MEOH_mgdL", "FORM_mM", "FCNS_mM", "FVIT_mM", "FOM_ugmL", "HCO3",
            "pH", "PACO2", "NA", "AG", "OG", "LACT", "ATP", "PUT", "OPTIC",
            "logMAR", "P_DEATH", "GFR_pct", "CA_ION"]
    print("t," + ",".join(keys))
    for tt in [0, 2, 4, 8, 12, 14, 15, 18, 21, 24, 30, 36, 48, 60, 72, 96]:
        i = int(np.argmin(np.abs(T - tt)))
        print("%.1f," % T[i] + ",".join("%.5g" % O[i][k] for k in keys))
    print("""
These 16 rows are the acceptance test for tap_mrgsolve_model.R.  The R model is
required to reproduce every column to 3 significant figures; the tolerance is
tight deliberately, because the two implementations share only the equations.""")


SECTIONS = [
    (1, sec1_kinetic_constants),
    (2, sec2_ion_trapping),
    (3, sec3_bicarb_and_urine),
    (4, sec4_bicarb_space),
    (5, sec5_two_gap_clock),
    (6, sec6_natural_history),
    (7, sec7_fomepizole_timing),
    (8, sec8_eg_vs_meoh_blocked),
    (9, sec9_hd_antidote_loop),
    (10, sec10_crystal_threshold),
    (11, sec11_cofactor_nonadditivity),
    (12, sec12_isopropanol_control),
    (13, sec13_prognostic_ranking),
    (14, sec14_conservation_checks),
    (15, sec15_reference_trajectory),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", type=int, action="append", default=None)
    a = ap.parse_args()
    want = set(a.section) if a.section else None
    for n, fn in SECTIONS:
        if want is None or n in want:
            fn()
    print()


if __name__ == "__main__":
    main()
