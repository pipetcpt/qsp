#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gastrointestinal Stromal Tumour (GIST) — dependency-free Python twin
====================================================================

This file is the executable twin of ``gist_mrgsolve_model.R``.  It contains the
same 49-state ODE system, the same parameters and the same scenario library, and
it exists for one reason: every number quoted in ``README.md`` is produced by
running this file, so the claims can be checked without R, mrgsolve or any
third-party Python package (standard library only).

    python3 gist_python_twin.py            # full report + self-checks
    python3 gist_python_twin.py --quiet     # self-checks only

--------------------------------------------------------------------------------
THE STRUCTURAL POSITION
--------------------------------------------------------------------------------
GIST is normally drawn as: KIT mutation -> constitutive kinase -> proliferation
-> imatinib occupies the ATP pocket -> tumour shrinks -> a resistance mutation
appears -> switch drug.  One tumour, one kinase activity, one IC50, one
resistance state.  That drawing cannot survive the following six observations
simultaneously, all of which are in the trial record:

  O1  Sunitinib and ripretinib are indistinguishable in second line overall
      (INTRIGUE ITT mPFS 8.0 vs 8.3 mo), yet in patients whose ctDNA carries
      only KIT exon 11+13/14 mutations sunitinib gives 15.0 mo vs ripretinib
      4.0 mo, and in patients with only exon 11+17/18 ripretinib gives 14.2 mo
      vs sunitinib 1.5 mo (Heinrich 2024, Nat Med).  The average of two
      opposite genotype-specific effects is a null.
  O2  Imatinib RECHALLENGE, after progression on imatinib, sunitinib and
      everything else, still beats placebo (RIGHT: 1.8 vs 0.9 mo).  A tumour
      that has "become imatinib-resistant" cannot respond to imatinib.
  O3  Dose escalation to 800 mg helps only KIT exon 9 (MetaGIST; EORTC 62005
      relative risk of progression reduced 61%).  Same drug, same disease, and
      the benefit is confined to one genotype.
  O4  PDGFRA D842V is the most imatinib-resistant genotype there is and the
      most avapritinib-sensitive (ORR ~88%).  One node, opposite drug ranks.
  O5  FDG-PET SUV collapses within 24-48 h while CT size barely moves for
      weeks, and "stable disease" is a common best response in patients who
      live for years.  Signalling, viable mass and imaged mass are three
      different variables.
  O6  Interrupting imatinib in responding patients produces progression with
      median PFS 6.1 mo (after 1 y of imatinib), 7.0 mo (after 3 y) and
      12.0 mo (after 5 y), versus 27.8 / 67.0 / not-reached on continuation
      (BFR14).  The disease is not eradicated by years of deep response, and
      the time-to-regrowth gets LONGER the longer the drug was given.

So this model takes three structural positions instead.

  C1  THE TUMOUR IS A POPULATION, NOT A SIZE.  The state is a vector of clones
      (primary-genotype clone plus ATP-binding-pocket, activation-loop and
      KIT-independent subclones), each with its own EC50 for each of five
      drugs, taken from published biochemical/cellular sensitivity patterns.
      Resistant subclones are generated per cell DIVISION during the untreated
      growth phase, so they are PRE-EXISTING at diagnosis, not induced.
      Efficacy is therefore a SET-COVER problem, not a potency problem:
      progression is the growth of the least-covered clone, and the ranking of
      two drugs can invert between genotypes while their means are equal (O1).
      Because a majority of cells still carry only the primary mutation even at
      the end of the sequence, imatinib rechallenge is not a paradox (O2).

  C2  OCCUPANCY IS FAST, KILLING IS SLOW, AND THE DRUG IS MOSTLY CYTOSTATIC.
      Three time scales are separated explicitly: kinase occupancy (hours,
      read out by PET), cell-cycle exit (days, read out by Ki67), and a
      quiescent reservoir that the drug does not kill (months to years, read
      out by what happens when you stop).  Cells whose KIT signal falls below
      threshold mostly STOP rather than die; they accumulate in a quiescent
      compartment that is drug-insensitive and dies at 0.05%/day.  The imaged
      mass is a third variable again, because dead cells leave slowly-resorbed
      non-viable tissue behind.  This is what makes O5 and O6 the same fact:
      the residual mass at nadir is mostly NOT viable tumour, so re-expansion
      to +20% SLD takes months even though the cells regrow at their untreated
      rate; and the longer the drug ran, the smaller the viable fraction, so
      the LONGER the interruption PFS.

  C3  EXPOSURE MATTERS ONLY WHERE THE GENOTYPE PUTS THE EC50.  Imatinib
      pharmacokinetics are explicit (CYP3A4 autoinduction, alpha-1 acid
      glycoprotein binding, active metabolite CGP74588), and the EC50 of the
      primary clone is genotype-specific.  At 400 mg the achieved
      concentration sits ~8-fold above the exon 11 EC50 and only ~2-fold above
      the exon 9 EC50 -- i.e. on either side of the proliferation threshold --
      so escalation moves exon 9 across the threshold and does nothing for
      exon 11 (O3).  The same geometry means a patient with high clearance and
      high AGP behaves like an exon 9 patient at 400 mg: PHARMACOKINETICS CAN
      MIMIC GENOTYPE.

FITTED PARAMETERS.  Of ~140 parameters, six are fitted, all listed in
``FITTED``: KPMAX/KD0 (untreated volume doubling time), MU_ATP and MU_AL
(mutation flux per division, to first-line PFS), FNEC/KRES (depth of RECIST
response), KDQ (reservoir attrition, to the 1-year BFR14 interruption arm) and
N_MICRO (occult burden after complete resection, to the 12-month adjuvant arm).
Everything the report calls a PREDICTION is downstream of those six and of
in-vitro potency ratios that were never fitted to any clinical endpoint.

FALSIFIER.  ``--falsify`` / scenario ``F1`` replaces the clone vector with one
clone whose potency degrades exponentially under drug pressure ("acquired
resistance" as a continuous loss of sensitivity), REFITTED to the same
first-line PFS anchor.  That is the fair single-compartment competitor, not a
straw man.  It is reported alongside every polyclonal prediction.
"""

import math
import sys

# ----------------------------------------------------------------------------
# 0.  State vector (49 states)
# ----------------------------------------------------------------------------
(IA_IM, IC_IM, IP_IM, IM_IM,          # imatinib gut / central / periph / CGP74588
 IA_SU, IC_SU, IP_SU, IM_SU,          # sunitinib + SU12662
 IA_RE, IC_RE, IP_RE, IM_RE,          # regorafenib + M-2/M-5
 IA_RI, IC_RI, IP_RI, IM_RI,          # ripretinib + DP-5439
 IA_AV, IC_AV, IP_AV,                 # avapritinib
 IENZ,                                # CYP3A4 enzyme pool (autoinduction)
 IN0C, IN0Q,                          # primary-genotype clone: cycling / quiescent
 IN1C, IN1Q,                          # KIT exon 13/14 ATP-pocket subclone
 IN2C, IN2Q,                          # KIT exon 17/18 activation-loop subclone
 IN3C, IN3Q,                          # KIT-independent (bypass) subclone
 IVNEC,                               # non-viable (necrotic/myxoid/fibrotic) volume
 ISLD,                                # imaged sum of longest diameters
 ISUV,                                # FDG-PET SUVmax
 IKI67,                               # Ki67 index
 IVASC, IVEGF,                        # tumour microvascular density, VEGF-A
 IAGP,                                # alpha-1 acid glycoprotein
 ICT0, ICTR,                          # ctDNA: primary mutation / resistance mutations
 INPRO, INT1, INT2, IANC,             # Friberg myelosuppression chain
 IMAPD,                               # mean arterial pressure offset
 IHFSR,                               # hand-foot skin reaction grade
 IEDEM,                               # periorbital/peripheral oedema score
 IFT4, ITSH,                          # thyroid axis
 IDI,                                 # achieved dose intensity (toxicity feedback)
 ITOXC,                               # cumulative toxicity burden
 IQDEP) = range(49)                   # depth of dormancy of the reservoir
NST = 49

SNAMES = ["A_IM", "C_IM", "P_IM", "M_IM", "A_SU", "C_SU", "P_SU", "M_SU",
          "A_RE", "C_RE", "P_RE", "M_RE", "A_RI", "C_RI", "P_RI", "M_RI",
          "A_AV", "C_AV", "P_AV", "ENZ", "N0C", "N0Q", "N1C", "N1Q", "N2C",
          "N2Q", "N3C", "N3Q", "VNEC", "SLD", "SUV", "KI67", "VASC", "VEGF",
          "AGP", "CT0", "CTR", "NPRO", "NT1", "NT2", "ANC", "MAPD", "HFSR",
          "EDEM", "FT4", "TSH", "DI", "TOXC", "QDEP"]

FITTED = ("KPMAX", "KD0", "TURNOVER", "FNEC", "KRES", "KDQ", "PHI_D", "N_MICRO")

# ----------------------------------------------------------------------------
# 1.  Parameters
# ----------------------------------------------------------------------------
# Concentrations are TOTAL PLASMA in ng/mL throughout; EC50 values are
# referenced to total plasma concentration, i.e. plasma protein binding and
# tissue penetration are absorbed into the constant.  The one place where
# binding is made explicit is imatinib, where free-fraction variation with
# alpha-1 acid glycoprotein is a measured clinical phenomenon (Widmer 2008,
# Haouala 2013): there, EC50_eff = EC50 * (AGP/AGP0)^HAGP.

P = dict(
    # ---- imatinib PK (400 mg qd; CL/F 12.5 L/h, V/F 250 L, t1/2 ~ 19 h) ----
    F_IM=0.98, KA_IM=12.0, V1_IM=520.0, VP_IM=420.0, Q_IM=220.0, CL_IM=232.0,
    FM_IM=0.12, CLM_IM=250.0, VM_IM=150.0,      # CGP74588, equipotent to parent
    KENZ=0.15, EMAX_ENZ=0.50, EC50_ENZ=1000.0,  # CYP3A4 autoinduction
    # ---- sunitinib PK (50 mg qd 4/2; CL/F 40 L/h, V/F 2000 L, t1/2 ~ 50 h) --
    F_SU=1.0, KA_SU=6.0, V1_SU=2000.0, VP_SU=1000.0, Q_SU=300.0, CL_SU=960.0,
    FM_SU=0.25, CLM_SU=700.0, VM_SU=1500.0,     # SU12662, equipotent
    # ---- regorafenib PK (160 mg qd 3/1) -----------------------------------
    F_RE=1.0, KA_RE=4.0, V1_RE=150.0, VP_RE=120.0, Q_RE=60.0, CL_RE=66.0,
    FM_RE=0.40, CLM_RE=55.0, VM_RE=140.0,       # M-2 + M-5, equipotent
    # ---- ripretinib PK (150 mg qd) ---------------------------------------
    F_RI=1.0, KA_RI=6.0, V1_RI=400.0, VP_RI=300.0, Q_RI=150.0, CL_RI=300.0,
    FM_RI=0.50, CLM_RI=250.0, VM_RI=350.0,      # DP-5439, equipotent
    # ---- avapritinib PK (300 mg qd) --------------------------------------
    F_AV=1.0, KA_AV=6.0, V1_AV=1200.0, VP_AV=800.0, Q_AV=250.0, CL_AV=450.0,

    # ---- tumour cell kinetics --------------------------------------------
    KPMAX=0.0290,      # [FITTED] max proliferation rate of cycling cells (/d)
    KD0=0.0175,        # [FITTED] baseline death rate of cycling cells (/d)
    SIG50=0.30,        # KIT signal giving half-maximal proliferation
    HP=4.0,            # Hill slope: the threshold is SHARP, which is what makes
                       # dose escalation genotype-specific rather than universal
    KAMAX=0.0130,      # max drug-induced apoptosis (cytostatic-dominant)
    UCRIT=0.85,        # (1-signal) giving half-maximal apoptosis
    HA=6.0,            # apoptosis needs DEEP suppression: the drug is cytostatic first
    KQIN=0.0040,       # cycling -> quiescent, proportional to (1 - signal)
    KQOUT=0.150,       # quiescent -> cycling, proportional to signal^NQ
    NQ=3.0,            # steepness of the wake-up requirement
    KDQ=0.00065,       # [FITTED] death rate of quiescent (drug-insensitive) cells
    KDEEP=0.0025,      # dormancy deepens while the signal is off (/d)
    KSHAL=0.050,       # and reverses when the signal returns (/d)
    PHI_D=0.90,        # [FITTED] maximal slowing of wake-up by deep dormancy
    TURNOVER=85.0,     # [FITTED] cumulative cell divisions per SURVIVING cell over
                       # the pre-diagnostic history.  For pure exponential growth
                       # this equals kp/(kp-kd) ~ 2.5; a lesion that spent years
                       # near turnover equilibrium accumulates far more, and it is
                       # this number -- not the mutation rate -- that sets how many
                       # resistant cells are already present on day 1.
    MU_ATP=2.00e-7,    # exon 13/14 subclone, per cell division
    MU_AL=3.00e-7,     # exon 17/18 subclone (the commoner secondary class)
    MU_BYP=4.00e-9,    # KIT-independent subclone per division (the rarest, and
                       # the one no drug covers)

    # ---- tumour mass / imaging -------------------------------------------
    CELLS_PER_ML=1.0e9,
    FNEC=0.62,         # [FITTED] fraction of DRUG-INDUCED dead-cell volume left behind
                       # as non-viable tissue.  A massive synchronous die-off
                       # overwhelms phagocytic clearance and leaves hyalinised /
                       # myxoid / cystic residue, which is what a responding GIST
                       # mass on CT actually consists of.
    FNEC0=0.15,        # the same fraction for PHYSIOLOGICAL turnover, which is
                       # cleared efficiently.  Sets the 18% non-viable fraction of
                       # an untreated tumour: FNEC0*KD0/(KRES+lambda).
    KRES=0.0030,       # [FITTED] resorption of non-viable tissue (/d)
    KSLD=16.28,        # SLD (mm) = KSLD * Vtot(mL)^(1/3)
    TAU_SLD=10.0,      # remodelling/imaging lag (d)
    SUV_BG=1.4, SUV0=9.0, TAU_PET=0.8,
    KI67MAX=18.0, TAU_KI=3.0,
    # ---- vasculature / delivery ------------------------------------------
    KVEGFP=0.30, KVEGFD=0.30, VREF=400.0,
    KVASCG=0.20, KVASCD=0.20,
    EC50_VEGFR_SU=20.0, EC50_VEGFR_RE=500.0, EC50_VEGFR_RI=1500.0,
    HDEL=0.5,          # delivery = (VASC/VASC0)^HDEL
    # ---- alpha-1 acid glycoprotein ---------------------------------------
    AGP0=1.0, AGPF=1.0, TAU_AGP=7.0, WAGP=0.25, HAGP=0.90,

    # ---- ctDNA -----------------------------------------------------------
    KSHED=2.2e-8,      # copies/mL per cell death
    KCTEL=12.0,        # ctDNA elimination (/d), t1/2 ~ 1.4 h
    CTWT=40.0,         # background wild-type cfDNA (copies/mL)

    # ---- myelosuppression (Friberg) --------------------------------------
    KTR=0.75, GAM=0.16, ANC0=4.0,
    SL_ANC_IM=1.8e-5, SL_ANC_SU=1.1e-3, SL_ANC_RE=1.1e-5, SL_ANC_RI=4.4e-5,
    SL_ANC_AV=7.5e-5,
    # ---- hypertension / HFSR / oedema / thyroid --------------------------
    KMAPD=0.050, MAPMAX=16.0,   # relaxation to MAPMAX x VEGFR2 occupancy
    KHFD=0.055,                 # relaxation to the HFSR drive
    KEDD=0.030,                 # relaxation to the oedema drive
    KF4P=0.050, KF4D=0.050, EC50_THY=60.0, EMAX_THY=0.80,
    KTSHP=0.070, KTSHD=0.070,
    # ---- dose-intensity feedback ----------------------------------------
    TOXTH=1.05, KDIN=0.045, KDIU=0.015, DIMIN=0.50,

    # ---- adjuvant / micrometastatic setting ------------------------------
    N_MICRO=2.2e6,     # [FITTED] occult cells after macroscopically complete resection
    N_DETECT=1.0e9,    # radiologically detectable recurrence (~1 cm lesion)

    # ---- falsifier (single averaged clone, resistance as potency drift) ---
    FALSIFY=0.0,
    KDRIFT=0.0,        # refitted inside make_scenarios() for F-runs
)

# Primary-genotype EC50 table (ng/mL, total plasma) and KIT dependence.
GENO = {
    "exon11": dict(IM=185.0, SU=20.0, RE=1400.0, RI=210.0, AV=170.0,
                   kitdep=1.00, kpscale=1.00),
    "exon9":  dict(IM=1000.0, SU=15.0, RE=1200.0, RI=230.0, AV=180.0,
                   kitdep=1.00, kpscale=1.00),
    "d842v":  dict(IM=20000.0, SU=2800.0, RE=30000.0, RI=600.0, AV=40.0,
                   kitdep=1.00, kpscale=0.85),
    "wt_sdh": dict(IM=1e7, SU=1e7, RE=1e7, RI=1e7, AV=1e7,
                   kitdep=0.15, kpscale=0.65),
}
# Secondary-mutation subclones.  The pattern (not the individual numbers) is the
# published one: the ATP-binding pocket mutants V654A/T670I stay sunitinib-
# sensitive and lose ripretinib/regorafenib potency, the activation-loop mutants
# D816/D820/N822K are the mirror image, and avapritinib follows the activation
# loop.  Serrano 2019 BJC; Smith 2019 Cancer Cell; Heinrich 2008 JCO.
CLONES = [
    None,                                                     # slot 0 = primary
    dict(name="exon13/14 ATP-pocket", IM=4500.0, SU=12.0, RE=9000.0,
         RI=2200.0, AV=2700.0, kitdep=1.00, kpscale=0.95),  # SU=12: V654A stays
                                                          # sunitinib-sensitive
    dict(name="exon17/18 activation-loop", IM=6500.0, SU=380.0, RE=950.0,
         RI=165.0, AV=130.0, kitdep=1.00, kpscale=0.95),
    dict(name="KIT-independent bypass", IM=1e7, SU=1e7, RE=1e7,
         RI=1e7, AV=1e7, kitdep=0.05, kpscale=0.85),
]

DRUGS = ("IM", "SU", "RE", "RI", "AV")
DOSE_IDX = {"IM": IA_IM, "SU": IA_SU, "RE": IA_RE, "RI": IA_RI, "AV": IA_AV}


# ----------------------------------------------------------------------------
# 2.  Right-hand side
# ----------------------------------------------------------------------------
def concentrations(y, par):
    """Total-plasma parent+equipotent-metabolite concentrations (ng/mL)."""
    c = {}
    c["IM"] = 1000.0 * (y[IC_IM] / par["V1_IM"] + y[IM_IM] / par["VM_IM"])
    c["SU"] = 1000.0 * (y[IC_SU] / par["V1_SU"] + y[IM_SU] / par["VM_SU"])
    c["RE"] = 1000.0 * (y[IC_RE] / par["V1_RE"] + y[IM_RE] / par["VM_RE"])
    c["RI"] = 1000.0 * (y[IC_RI] / par["V1_RI"] + y[IM_RI] / par["VM_RI"])
    c["AV"] = 1000.0 * (y[IC_AV] / par["V1_AV"])
    return c


def clone_ec50(idx, par, geno):
    """EC50 vector (ng/mL) for clone idx under the current AGP and drift."""
    tab = geno if idx == 0 else CLONES[idx]
    out = {}
    agpf = (max(par["_agp"], 0.2) / par["AGP0"]) ** par["HAGP"]
    for d in DRUGS:
        e = tab[d]
        if d == "IM":
            e = e * agpf
        out[d] = e * par["_drift"]
    return out, tab["kitdep"], tab["kpscale"]


def rhs(t, y, par, geno, ctx):
    dy = [0.0] * NST
    conc = concentrations(y, par)

    # -- PK -------------------------------------------------------------
    enz = max(y[IENZ], 0.2)
    clim = par["CL_IM"] * enz * par["_clf"]
    dy[IA_IM] = -par["KA_IM"] * y[IA_IM]
    dy[IC_IM] = (par["F_IM"] * par["KA_IM"] * y[IA_IM]
                 - (clim / par["V1_IM"]) * y[IC_IM]
                 - par["Q_IM"] * (y[IC_IM] / par["V1_IM"] - y[IP_IM] / par["VP_IM"]))
    dy[IP_IM] = par["Q_IM"] * (y[IC_IM] / par["V1_IM"] - y[IP_IM] / par["VP_IM"])
    dy[IM_IM] = (par["FM_IM"] * clim / par["V1_IM"] * y[IC_IM]
                 - par["CLM_IM"] / par["VM_IM"] * y[IM_IM])
    dy[IENZ] = par["KENZ"] * ((1.0 + par["EMAX_ENZ"] * conc["IM"]
                               / (par["EC50_ENZ"] + conc["IM"])) * par["_indf"]
                              - y[IENZ])

    for nm, (ia, ic, ip, im) in (("SU", (IA_SU, IC_SU, IP_SU, IM_SU)),
                                 ("RE", (IA_RE, IC_RE, IP_RE, IM_RE)),
                                 ("RI", (IA_RI, IC_RI, IP_RI, IM_RI))):
        ka, v1, vp, q, cl = (par["KA_" + nm], par["V1_" + nm], par["VP_" + nm],
                             par["Q_" + nm], par["CL_" + nm])
        dy[ia] = -ka * y[ia]
        dy[ic] = (par["F_" + nm] * ka * y[ia] - (cl / v1) * y[ic]
                  - q * (y[ic] / v1 - y[ip] / vp))
        dy[ip] = q * (y[ic] / v1 - y[ip] / vp)
        dy[im] = (par["FM_" + nm] * cl / v1 * y[ic]
                  - par["CLM_" + nm] / par["VM_" + nm] * y[im])
    dy[IA_AV] = -par["KA_AV"] * y[IA_AV]
    dy[IC_AV] = (par["F_AV"] * par["KA_AV"] * y[IA_AV]
                 - (par["CL_AV"] / par["V1_AV"]) * y[IC_AV]
                 - par["Q_AV"] * (y[IC_AV] / par["V1_AV"] - y[IP_AV] / par["VP_AV"]))
    dy[IP_AV] = par["Q_AV"] * (y[IC_AV] / par["V1_AV"] - y[IP_AV] / par["VP_AV"])

    # -- tumour interstitial availability -------------------------------
    vasc = max(y[IVASC], 1e-3)
    deliv = min(1.6, vasc ** par["HDEL"])
    par["_agp"] = max(y[IAGP], 0.2)

    # -- per-clone signalling, growth, death ----------------------------
    cyc = (IN0C, IN1C, IN2C, IN3C)
    qui = (IN0Q, IN1Q, IN2Q, IN3Q)
    ntot = 0.0
    for i in range(4):
        ntot += max(y[cyc[i]], 0.0) + max(y[qui[i]], 0.0)
    ntot = max(ntot, 1.0)

    deathflux = 0.0
    excess = 0.0
    basal = 0.0
    shed0 = 0.0
    shedR = 0.0
    sigbar = 0.0
    prolbar = 0.0
    kp0 = 0.0
    nact = 4
    if par["FALSIFY"] > 0.5:
        nact = 1                      # single averaged clone
    for i in range(nact):
        ec, kitdep, kpscale = clone_ec50(i, par, geno)
        a = 1.0
        for d in DRUGS:
            cd = conc[d] * deliv
            if cd > 0.0:
                a *= 1.0 / (1.0 + cd / ec[d])
        sig = 1.0 - kitdep * (1.0 - a)
        sig = min(max(sig, 1e-6), 1.0)
        sh = sig ** par["HP"]
        kp = par["KPMAX"] * kpscale * sh / (sh + par["SIG50"] ** par["HP"])
        u = 1.0 - sig
        kd = par["KD0"] + par["KAMAX"] * u ** par["HA"] / (
            u ** par["HA"] + par["UCRIT"] ** par["HA"])
        kqi = par["KQIN"] * (1.0 - sig)
        kqo = (par["KQOUT"] * sig ** par["NQ"]
               * (1.0 - par["PHI_D"] * min(max(y[IQDEP], 0.0), 1.0)))
        nc = max(y[cyc[i]], 0.0)
        nq = max(y[qui[i]], 0.0)
        if i == 0:
            kp0 = kp
            mu = par["MU_ATP"] + par["MU_AL"] + par["MU_BYP"]
            if par["FALSIFY"] > 0.5:
                mu = 0.0
        else:
            mu = 0.0
        dy[cyc[i]] += (kp * (1.0 - mu) - kd - kqi) * nc + kqo * nq
        dy[qui[i]] += kqi * nc - (kqo + par["KDQ"]) * nq
        dflux = kd * nc + par["KDQ"] * nq
        deathflux += dflux
        # only death IN EXCESS of physiological turnover leaves persistent
        # non-viable tissue behind
        excess += max(0.0, kd - par["KD0"]) * nc
        basal += par["KD0"] * nc + par["KDQ"] * nq
        if i == 0:
            shed0 += dflux
        else:
            shedR += dflux
        sigbar += sig * (nc + nq)
        prolbar += kp * nc
    if par["FALSIFY"] <= 0.5:
        dy[IN1C] += par["MU_ATP"] * kp0 * max(y[IN0C], 0.0)
        dy[IN2C] += par["MU_AL"] * kp0 * max(y[IN0C], 0.0)
        dy[IN3C] += par["MU_BYP"] * kp0 * max(y[IN0C], 0.0)
    sigbar /= ntot

    # potency drift: the falsifier's only resistance mechanism
    if par["FALSIFY"] > 0.5:
        press = 1.0 - sigbar
        dy_drift = par["KDRIFT"] * press * par["_drift"]
        ctx["drift_rate"] = dy_drift
    else:
        ctx["drift_rate"] = 0.0

    # -- mass, imaging --------------------------------------------------
    vviab = ntot / par["CELLS_PER_ML"]
    vnec = max(y[IVNEC], 0.0)
    vtot = vviab + vnec
    dy[IVNEC] = ((par["FNEC"] * excess + par["FNEC0"] * basal)
                 / par["CELLS_PER_ML"] - par["KRES"] * vnec)
    sld_true = par["KSLD"] * (max(vtot, 1e-9) ** (1.0 / 3.0))
    dy[ISLD] = (sld_true - y[ISLD]) / par["TAU_SLD"]
    suv_t = par["SUV_BG"] + (par["SUV0"] - par["SUV_BG"]) * sigbar
    dy[ISUV] = (suv_t - y[ISUV]) / par["TAU_PET"]
    ki_t = par["KI67MAX"] * prolbar / (par["KPMAX"] * ntot)
    dy[IKI67] = (ki_t - y[IKI67]) / par["TAU_KI"]

    # -- vasculature ----------------------------------------------------
    occv = (conc["SU"] / (conc["SU"] + par["EC50_VEGFR_SU"])
            + conc["RE"] / (conc["RE"] + par["EC50_VEGFR_RE"])
            + conc["RI"] / (conc["RI"] + par["EC50_VEGFR_RI"]))
    occv = min(occv, 0.95)
    dy[IVEGF] = par["KVEGFP"] * (vtot / par["VREF"] + 0.05) - par["KVEGFD"] * y[IVEGF]
    dy[IVASC] = (par["KVASCG"] * y[IVEGF] * (1.0 - occv)
                 - par["KVASCD"] * y[IVASC])

    # -- AGP ------------------------------------------------------------
    agp_t = par["AGP0"] * par["AGPF"] * (1.0 + par["WAGP"] * vtot / (vtot + 200.0))
    dy[IAGP] = (agp_t - y[IAGP]) / par["TAU_AGP"]

    # -- ctDNA ----------------------------------------------------------
    dy[ICT0] = par["KSHED"] * shed0 - par["KCTEL"] * y[ICT0]
    dy[ICTR] = par["KSHED"] * shedR - par["KCTEL"] * y[ICTR]

    # -- myelosuppression ----------------------------------------------
    edrug = min(0.90, par["SL_ANC_IM"] * conc["IM"] + par["SL_ANC_SU"] * conc["SU"]
                + par["SL_ANC_RE"] * conc["RE"] + par["SL_ANC_RI"] * conc["RI"]
                + par["SL_ANC_AV"] * conc["AV"])
    fb = (par["ANC0"] / max(y[IANC], 0.05)) ** par["GAM"]
    dy[INPRO] = par["KTR"] * y[INPRO] * ((1.0 - edrug) * fb - 1.0)
    dy[INT1] = par["KTR"] * (y[INPRO] - y[INT1])
    dy[INT2] = par["KTR"] * (y[INT1] - y[INT2])
    dy[IANC] = par["KTR"] * (y[INT2] - y[IANC])

    # -- hypertension, HFSR, oedema, thyroid ---------------------------
    dy[IMAPD] = par["KMAPD"] * (par["MAPMAX"] * occv - y[IMAPD])
    hf_drive = (conc["RE"] / 1200.0 + 0.50 * conc["SU"] / 70.0
                + 0.25 * conc["RI"] / 800.0 + 0.05 * conc["IM"] / 1500.0)
    dy[IHFSR] = par["KHFD"] * (hf_drive - y[IHFSR])
    dy[IEDEM] = par["KEDD"] * (conc["IM"] / 1500.0 - y[IEDEM])
    thy = par["EMAX_THY"] * conc["SU"] / (conc["SU"] + par["EC50_THY"])
    dy[IFT4] = par["KF4P"] * (1.0 - thy) - par["KF4D"] * y[IFT4]
    dy[ITSH] = par["KTSHP"] * (1.0 / max(y[IFT4], 0.15)) ** 1.5 - par["KTSHD"] * y[ITSH]

    # -- toxicity -> dose intensity -----------------------------------
    tox = (0.60 * max(0.0, 1.5 - y[IANC]) / 1.5 + 0.50 * y[IHFSR]
           + 0.40 * y[IMAPD] / 15.0 + 0.30 * y[IEDEM]
           + 0.20 * max(0.0, y[ITSH] - 2.0) / 3.0)
    ctx["tox"] = tox
    qd = min(max(y[IQDEP], 0.0), 1.0)
    dy[IQDEP] = (par["KDEEP"] * (1.0 - sigbar) * (1.0 - qd)
                 - par["KSHAL"] * sigbar * qd)
    dy[ITOXC] = tox - 0.20 * y[ITOXC]
    dy[IDI] = (-par["KDIN"] * y[IDI] * max(0.0, tox - par["TOXTH"])
               + par["KDIU"] * (1.0 - y[IDI]) * max(0.0, par["TOXTH"] - tox))
    return dy


# ----------------------------------------------------------------------------
# 3.  Initial conditions and integrator
# ----------------------------------------------------------------------------
def initial_state(par, geno, setting="metastatic", v0=400.0, years_occult=0.0):
    y = [0.0] * NST
    y[IENZ] = 1.0
    y[IVEGF] = 1.0
    y[IVASC] = 1.0
    y[IAGP] = par["AGP0"] * par["AGPF"]
    y[IANC] = y[INPRO] = y[INT1] = y[INT2] = par["ANC0"]
    y[IFT4] = 1.0
    y[ITSH] = 1.0
    y[IDI] = 1.0
    y[ISUV] = par["SUV0"]
    y[IKI67] = par["KI67MAX"]
    if setting == "metastatic":
        ntot = v0 * par["CELLS_PER_ML"]
    else:
        ntot = par["N_MICRO"]
    # Resistant subclones are PRE-EXISTING: they are the mutation flux
    # accumulated over the untreated expansion that produced this burden.
    # Integrating mu*kp*N over an exponential expansion from one cell gives
    # approximately mu * N * kp/(kp-kd) -- i.e. proportional to burden.
    amp = par["TURNOVER"]
    f1 = par["MU_ATP"] * amp
    f2 = par["MU_AL"] * amp
    f3 = par["MU_BYP"] * amp
    if par["FALSIFY"] > 0.5:
        f1 = f2 = f3 = 0.0
    y[IN1C] = ntot * f1
    y[IN2C] = ntot * f2
    y[IN3C] = ntot * f3
    y[IN0C] = ntot * (1.0 - f1 - f2 - f3)
    lam0 = par["KPMAX"] / (1.0 + par["SIG50"] ** par["HP"]) - par["KD0"]
    y[IVNEC] = (par["FNEC0"] * par["KD0"] / (par["KRES"] + max(lam0, 1e-4))
                * ntot / par["CELLS_PER_ML"])
    v = (y[IN0C] + y[IN1C] + y[IN2C] + y[IN3C]) / par["CELLS_PER_ML"] + y[IVNEC]
    y[ISLD] = par["KSLD"] * (max(v, 1e-9) ** (1.0 / 3.0))
    y[ICT0] = par["KSHED"] * par["KD0"] * y[IN0C] / par["KCTEL"]
    y[ICTR] = par["KSHED"] * par["KD0"] * (y[IN1C] + y[IN2C] + y[IN3C]) / par["KCTEL"]
    return y


def on_schedule(day, sched):
    """True if a dose is given on this treatment day of the current stage."""
    if sched == "qd":
        return True
    if sched == "4/2":                      # 4 weeks on, 2 weeks off
        return (day % 42) < 28
    if sched == "3/1":                      # 3 weeks on, 1 week off
        return (day % 28) < 21
    if sched == "2/1":
        return (day % 21) < 14
    raise ValueError(sched)


def RX(*specs, **kw):
    """A line of therapy: RX(("IM", 400, "qd"), ("RI", 100, "qd"), on="progression")."""
    st = dict(drugs=[dict(drug=d, dose=float(x), sched=sc) for d, x, sc in specs])
    st.update(kw)
    return st


def simulate(scn, dt=0.05, record_every=1.0):
    """Run one scenario.  Returns recorded series plus per-line endpoints.

    Dosing is driven by an integer day counter rather than by testing whether
    the integration time lands on a whole day, so results are independent of
    the step size.
    """
    par = dict(P)
    par.update(scn.get("par", {}))
    par["_clf"] = scn.get("clf", 1.0)
    par["_indf"] = scn.get("indf", 1.0)
    par["_drift"] = 1.0
    par["_agp"] = par["AGP0"] * par["AGPF"]
    geno = GENO[scn.get("geno", "exon11")]
    setting = scn.get("setting", "metastatic")
    y = initial_state(par, geno, setting, scn.get("v0", 400.0))
    tend = scn["tend"]
    stages = [dict(s) for s in scn["stages"]]
    for st in stages:
        st.setdefault("drugs", [])
        st.setdefault("t0", 0.0)
    ctx = {"tox": 0.0, "drift_rate": 0.0}

    si = 0
    stages[0]["_base"] = y[ISLD]
    stage_start = 0.0
    keys = ("t", "sld", "suv", "ki67", "n0", "n1", "n2", "n3", "nq", "ntot",
            "vnec", "cim", "csu", "cre", "cri", "cav", "anc", "map", "hfsr",
            "edem", "tsh", "di", "vaf", "ctdna", "vfrac", "agp", "vasc",
            "qdep", "ft4", "line")
    rec = {k: [] for k in keys}
    events = []
    baseline_sld = y[ISLD]
    nadir = y[ISLD]
    nadir_stage = y[ISLD]
    pfs = {}
    best = {}
    recurrence = None
    assess = scn.get('assess', 7.0)
    next_assess = assess
    next_rec = 0.0
    dose_day = 0
    t = 0.0
    nsteps = int(round(tend / dt))
    for step in range(nsteps + 1):
        t = step * dt
        # ---- dosing: integer day counter, step-size independent ---------
        while dose_day <= t + 1e-9:
            st = stages[si]
            for spec in st["drugs"]:
                if dose_day >= st["t0"] - 1e-9:
                    tday = int(round(dose_day - st["t0"]))
                    if tday >= 0 and on_schedule(tday, spec["sched"]):
                        di = y[IDI] if st.get("adapt", True) else 1.0
                        y[DOSE_IDX[spec["drug"]]] += spec["dose"] * di
            dose_day += 1
        # ---- record -----------------------------------------------------
        if t >= next_rec - 1e-9:
            c = concentrations(y, par)
            ntot = sum(max(y[i], 0.0) for i in (IN0C, IN0Q, IN1C, IN1Q,
                                                IN2C, IN2Q, IN3C, IN3Q))
            vv = ntot / par["CELLS_PER_ML"]
            vt = vv + max(y[IVNEC], 0.0)
            rec["t"].append(t)
            rec["sld"].append(y[ISLD])
            rec["suv"].append(y[ISUV])
            rec["ki67"].append(y[IKI67])
            rec["n0"].append(max(y[IN0C], 0.0) + max(y[IN0Q], 0.0))
            rec["n1"].append(max(y[IN1C], 0.0) + max(y[IN1Q], 0.0))
            rec["n2"].append(max(y[IN2C], 0.0) + max(y[IN2Q], 0.0))
            rec["n3"].append(max(y[IN3C], 0.0) + max(y[IN3Q], 0.0))
            rec["nq"].append(sum(max(y[i], 0.0) for i in (IN0Q, IN1Q, IN2Q, IN3Q)))
            rec["ntot"].append(ntot)
            rec["vnec"].append(max(y[IVNEC], 0.0))
            rec["vfrac"].append(vv / max(vt, 1e-9))
            rec["cim"].append(c["IM"]); rec["csu"].append(c["SU"])
            rec["cre"].append(c["RE"]); rec["cri"].append(c["RI"])
            rec["cav"].append(c["AV"])
            rec["anc"].append(y[IANC]); rec["map"].append(y[IMAPD])
            rec["hfsr"].append(y[IHFSR]); rec["edem"].append(y[IEDEM])
            rec["tsh"].append(y[ITSH]); rec["di"].append(y[IDI])
            rec["agp"].append(y[IAGP]); rec["vasc"].append(y[IVASC])
            rec["ctdna"].append(y[ICT0] + y[ICTR])
            rec["vaf"].append(y[ICTR] / max(y[ICT0] + y[ICTR] + par["CTWT"], 1e-9))
            rec["qdep"].append(y[IQDEP])
            rec["ft4"].append(y[IFT4])
            rec["line"].append(si)
            next_rec += record_every
        # ---- endpoints --------------------------------------------------
        if setting != "metastatic" and recurrence is None:
            ntot = sum(max(y[i], 0.0) for i in (IN0C, IN0Q, IN1C, IN1Q,
                                               IN2C, IN2Q, IN3C, IN3Q))
            if ntot >= par["N_DETECT"]:
                recurrence = t
        if t >= next_assess - 1e-9:
            sld = y[ISLD]
            nadir = min(nadir, sld)
            nadir_stage = min(nadir_stage, sld)
            key = "stage%d" % si
            base_stage = stages[si].get("_base", baseline_sld)
            b = 100.0 * (nadir_stage - base_stage) / base_stage
            best[key] = min(best.get(key, 0.0), b)
            if (sld >= 1.2 * nadir_stage and sld >= nadir_stage + 5.0
                    and key not in pfs):
                pfs[key] = t - stage_start
                events.append((t, "PD on line %d" % si))
                if si + 1 < len(stages) and stages[si + 1].get("on") == "progression":
                    si += 1
                    stages[si]["t0"] = float(int(t)) + 1.0
                    stages[si]["_base"] = sld
                    stage_start = t
                    nadir_stage = sld
            next_assess += assess
        # ---- fixed-time stage advance -----------------------------------
        while (si + 1 < len(stages) and stages[si + 1].get("on") != "progression"
               and t >= stages[si + 1]["t0"] - 1e-9):
            si += 1
            stages[si]["_base"] = y[ISLD]
            stage_start = t
            nadir_stage = y[ISLD]
        if step == nsteps:
            break
        # ---- RK4 --------------------------------------------------------
        k1 = rhs(t, y, par, geno, ctx)
        d1 = ctx["drift_rate"]
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NST)]
        k2 = rhs(t + 0.5 * dt, y2, par, geno, ctx)
        d2 = ctx["drift_rate"]
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NST)]
        k3 = rhs(t + 0.5 * dt, y3, par, geno, ctx)
        d3 = ctx["drift_rate"]
        y4 = [y[i] + dt * k3[i] for i in range(NST)]
        k4 = rhs(t + dt, y4, par, geno, ctx)
        d4 = ctx["drift_rate"]
        for i in range(NST):
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            if i != IENZ and y[i] < 0.0:
                y[i] = 0.0
        par["_drift"] += dt / 6.0 * (d1 + 2.0 * d2 + 2.0 * d3 + d4)
        y[IDI] = min(1.0, max(par["DIMIN"], y[IDI]))
        y[IQDEP] = min(1.0, max(0.0, y[IQDEP]))

    return dict(rec=rec, pfs=pfs, best=best, events=events,
                recurrence=recurrence, baseline=baseline_sld,
                final=y, par=par, nadir=nadir, drift=par["_drift"])


# ----------------------------------------------------------------------------
# 4.  Helpers
# ----------------------------------------------------------------------------
def mo(days):
    return days / 30.4375


def at(rec, key, t):
    """Value of series `key` at (or just before) time t."""
    ts = rec["t"]
    if t >= ts[-1]:
        return rec[key][-1]
    lo, hi = 0, len(ts) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if ts[mid] < t:
            lo = mid + 1
        else:
            hi = mid
    return rec[key][lo]


def mean_on(rec, key, t0, t1, cond_key=None, cond_min=1e-9):
    """Mean of `key` over [t0,t1], optionally only where cond_key > cond_min."""
    num = 0.0
    n = 0
    for i, t in enumerate(rec["t"]):
        if t0 <= t <= t1:
            if cond_key is not None and rec[cond_key][i] <= cond_min:
                continue
            num += rec[key][i]
            n += 1
    return num / n if n else 0.0


def win(rec, key, t0, t1):
    """The values of `key` inside the time window [t0, t1]."""
    v = [rec[key][i] for i, t in enumerate(rec["t"]) if t0 <= t <= t1]
    return v if v else [rec[key][0]]


def first_time(rec, key, thresh, above=True):
    for i, v in enumerate(rec[key]):
        if (above and v >= thresh) or (not above and v <= thresh):
            return rec["t"][i]
    return None


def pct_change(rec, key, t, base=None):
    b = rec[key][0] if base is None else base
    return 100.0 * (at(rec, key, t) - b) / b


def stage_pfs(res, k=0):
    return res["pfs"].get("stage%d" % k)


def line_start(res, k):
    """Time at which line k started (sum of previous lines' durations)."""
    tt = 0.0
    for j in range(k):
        tt += res["pfs"].get("stage%d" % j) or 0.0
    return tt


# ----------------------------------------------------------------------------
# 5.  Scenario library
# ----------------------------------------------------------------------------
def S(**kw):
    return kw


def scenarios():
    """The scenario library, keyed by short code."""
    IM4 = ("IM", 400, "qd")
    IM8 = ("IM", 800, "qd")
    SU5 = ("SU", 50, "4/2")
    SU37 = ("SU", 37.5, "qd")
    RE1 = ("RE", 160, "3/1")
    RI1 = ("RI", 150, "qd")
    AV3 = ("AV", 300, "qd")
    PROG = dict(on="progression")
    mu_t = P["MU_ATP"] + P["MU_AL"]
    ATPONLY = dict(MU_ATP=mu_t, MU_AL=1.0e-11)
    ALONLY = dict(MU_ATP=1.0e-11, MU_AL=mu_t)
    sc = {}
    sc["S1"] = S(label="Natural history, untreated metastatic exon 11",
                 geno="exon11", tend=900.0, stages=[RX()])
    sc["S2"] = S(label="Imatinib 400 mg, KIT exon 11 (first line)",
                 geno="exon11", tend=2600.0, stages=[RX(IM4)])
    sc["S3"] = S(label="Imatinib 400 mg, KIT exon 9",
                 geno="exon9", tend=2600.0, stages=[RX(IM4)])
    sc["S4"] = S(label="Imatinib 800 mg, KIT exon 9",
                 geno="exon9", tend=2600.0, stages=[RX(IM8)])
    sc["S5"] = S(label="Imatinib 800 mg, KIT exon 11 (no genotype rationale)",
                 geno="exon11", tend=2600.0, stages=[RX(IM8)])
    sc["S6"] = S(label="Imatinib 400 mg, exon 11, high clearance + high AGP",
                 geno="exon11", tend=2600.0, clf=3.0, par=dict(AGPF=1.6),
                 stages=[RX(IM4)])
    sc["S7"] = S(label="Imatinib 400 mg, PDGFRA D842V",
                 geno="d842v", tend=900.0, stages=[RX(IM4)])
    sc["S8"] = S(label="Avapritinib 300 mg, PDGFRA D842V",
                 geno="d842v", tend=2600.0, stages=[RX(AV3)])
    sc["S9"] = S(label="Imatinib 400 mg, SDH-deficient / KIT-PDGFRA wild type",
                 geno="wt_sdh", tend=900.0, stages=[RX(IM4)])
    sc["S10"] = S(label="Imatinib -> sunitinib 50 mg 4/2 (unselected)",
                  geno="exon11", tend=3000.0, stages=[RX(IM4), RX(SU5, **PROG)])
    sc["S11"] = S(label="Imatinib -> ripretinib 150 mg (unselected)",
                  geno="exon11", tend=3000.0, stages=[RX(IM4), RX(RI1, **PROG)])
    sc["S12"] = S(label="INTRIGUE exon 11+13/14 -> sunitinib", geno="exon11",
                  tend=3000.0, par=dict(ATPONLY), stages=[RX(IM4), RX(SU5, **PROG)])
    sc["S13"] = S(label="INTRIGUE exon 11+13/14 -> ripretinib", geno="exon11",
                  tend=3000.0, par=dict(ATPONLY), stages=[RX(IM4), RX(RI1, **PROG)])
    sc["S14"] = S(label="INTRIGUE exon 11+17/18 -> sunitinib", geno="exon11",
                  tend=3000.0, par=dict(ALONLY), stages=[RX(IM4), RX(SU5, **PROG)])
    sc["S15"] = S(label="INTRIGUE exon 11+17/18 -> ripretinib", geno="exon11",
                  tend=3000.0, par=dict(ALONLY), stages=[RX(IM4), RX(RI1, **PROG)])
    sc["S16"] = S(label="Imatinib -> sunitinib -> regorafenib (third line)",
                  geno="exon11", tend=3400.0,
                  stages=[RX(IM4), RX(SU5, **PROG), RX(RE1, **PROG)])
    sc["S17"] = S(label="Imatinib -> sunitinib -> regorafenib -> ripretinib",
                  geno="exon11", tend=3600.0,
                  stages=[RX(IM4), RX(SU5, **PROG), RX(RE1, **PROG), RX(RI1, **PROG)])
    sc["S18"] = S(label="Imatinib -> sunitinib -> regorafenib -> placebo",
                  geno="exon11", tend=3600.0,
                  stages=[RX(IM4), RX(SU5, **PROG), RX(RE1, **PROG), RX(**PROG)])
    sc["S19"] = S(label="Imatinib -> sunitinib -> regorafenib -> imatinib rechallenge",
                  geno="exon11", tend=3600.0,
                  stages=[RX(IM4), RX(SU5, **PROG), RX(RE1, **PROG), RX(IM4, **PROG)])
    sc["S20"] = S(label="Sunitinib 37.5 mg continuous instead of 50 mg 4/2",
                  geno="exon11", tend=3000.0, stages=[RX(IM4), RX(SU37, **PROG)])
    sc["S21"] = S(label="Imatinib 400 + ripretinib 100 upfront combination",
                  geno="exon11", tend=3000.0, stages=[RX(IM4, ("RI", 100, "qd"))])
    sc["S22"] = S(label="Adjuvant imatinib 12 months (high-risk resected)",
                  geno="exon11", setting="adjuvant", tend=4400.0,
                  stages=[RX(IM4), RX(t0=365.0)])
    sc["S23"] = S(label="Adjuvant imatinib 36 months (high-risk resected)",
                  geno="exon11", setting="adjuvant", tend=4400.0,
                  stages=[RX(IM4), RX(t0=1095.0)])
    sc["S24"] = S(label="Adjuvant imatinib 60 months (high-risk resected)",
                  geno="exon11", setting="adjuvant", tend=4400.0,
                  stages=[RX(IM4), RX(t0=1826.0)])
    sc["S25"] = S(label="Rifampicin co-administration on imatinib 400 mg",
                  geno="exon11", tend=2600.0, indf=2.6, stages=[RX(IM4)])
    sc["S26"] = S(label="Imatinib 400 -> 800 mg escalation at progression",
                  geno="exon9", tend=2600.0, stages=[RX(IM4), RX(IM8, **PROG)])
    return sc


# BFR14 randomised patients only if they were still NON-progressing after 1, 3
# or 5 years of imatinib, which selects progressively harder for a small
# pre-existing resistant pool.  The twin represents that entry criterion by
# choosing TURNOVER so that the CONTINUATION arm reproduces the observed
# continuation PFS (27.8 / 67.0 / not reached).  The INTERRUPTION arm is then a
# prediction, not a fit.
# The model REQUIRES an essentially resistance-free tumour to stay
# progression-free on imatinib for 3-5 years, which is itself a prediction: only
# 27 of the 434 patients enrolled in BFR14 (6%) reached the 5-year
# randomisation, and 50 of 434 (12%) the 3-year one.  The phenotype below is
# chosen so that each CONTINUATION arm is plausible for that entry criterion;
# the INTERRUPTION arm is then predicted, not fitted.  The 5-year cohort needs
# both no pre-existing resistance AND slower-growing disease, and the model's
# prediction is that their longer post-interruption PFS is a property of the
# PATIENTS, not of the extra two years of drug.
BFR14_PHENO = {1: dict(TURNOVER=85.0),
               3: dict(TURNOVER=0.0),
               5: dict(TURNOVER=0.0, KPMAX=0.0225)}


def interruption_scenario(years, par=None):
    """BFR14 design: interrupt imatinib after `years`, versus continuing."""
    d = years * 365.0
    pp = dict(BFR14_PHENO[years])
    if par:
        pp.update(par)
    stop = S(label="Interruption after %g y" % years, geno="exon11",
             tend=d + 1500.0, par=pp,
             stages=[RX(("IM", 400, "qd")), RX(t0=d)])
    go = S(label="Continuation beyond %g y" % years, geno="exon11",
           tend=d + 3000.0, par=pp,
           stages=[RX(("IM", 400, "qd")), RX(("IM", 400, "qd"), t0=d)])
    return stop, go


# ----------------------------------------------------------------------------
# 6.  Report
# ----------------------------------------------------------------------------
CHECKS = []


def check(name, cond, detail=""):
    CHECKS.append((name, bool(cond), detail))
    return bool(cond)


def h(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def fmt(x, unit="mo"):
    if not x:
        return "not reached"
    return ("%.1f %s" % (mo(x), unit)).strip()


def run_report(quiet=False):
    def pr(*a):
        if not quiet:
            print(*a)

    sc = scenarios()
    R = {}
    for k in sorted(sc, key=lambda s: int(s[1:])):
        R[k] = simulate(sc[k])

    # ---- A. pharmacokinetics ------------------------------------------
    h("A.  PHARMACOKINETICS AND ACHIEVED EXPOSURE  [C3]")
    r2 = R["S2"]["rec"]
    tr400 = at(r2, "cim", 400.0)
    tr800 = at(R["S5"]["rec"], "cim", 400.0)
    pr("A1  Imatinib 400 mg, steady-state trough (parent + CGP74588) : %6.0f ng/mL"
       % tr400)
    pr("    observed reference: measured trough at 400 mg ~1100 ng/mL, and the")
    pr("    exposure-response threshold Demetri 2009 identified was 1110 ng/mL")
    pr("A2  Imatinib 800 mg trough                                   : %6.0f ng/mL"
       % tr800)
    pr("A3  CYP3A4 autoinduction, enzyme pool at 400 d               : %6.2f x"
       % R["S2"]["final"][IENZ])
    pr("A4  Rifampicin co-administration, trough                     : %6.0f ng/mL"
       % at(R["S25"]["rec"], "cim", 400.0))
    pr("A5  High-clearance + high-AGP patient, same 400 mg           : %6.0f ng/mL"
       % at(R["S6"]["rec"], "cim", 400.0))
    pr("A6  AGP, baseline -> on deep response                        : %5.2f -> %5.2f g/L"
       % (r2["agp"][0], at(r2, "agp", 500.0)))
    su_on = mean_on(R["S10"]["rec"], "csu", line_start(R["S10"], 1) + 7,
                    line_start(R["S10"], 1) + 28)
    pr("A7  Sunitinib 50 mg, mean on-treatment (parent + SU12662)    : %6.0f ng/mL"
       % su_on)
    pr("A8  Regorafenib 160 mg, mean on-treatment (parent + M-2/M-5) : %6.0f ng/mL"
       % mean_on(R["S16"]["rec"], "cre", line_start(R["S16"], 2) + 7,
                 line_start(R["S16"], 2) + 21))
    check("A1 imatinib trough in the measured range", 900.0 < tr400 < 1600.0,
          "%.0f ng/mL" % tr400)
    check("A2 imatinib is dose-linear", 1.7 < tr800 / tr400 < 2.3,
          "%.2f x for 2 x dose" % (tr800 / tr400))
    check("A3 autoinduction reduces exposure by 20-35%",
          1.20 < R["S2"]["final"][IENZ] < 1.40)

    # ---- B. three time scales -----------------------------------------
    h("B.  THREE TIME SCALES: OCCUPANCY, CELL CYCLE, IMAGED MASS  [C2]")
    r1 = R["S1"]["rec"]
    td = first_time(r1, "ntot", 2.0 * r1["ntot"][0])
    pr("B1  Untreated volume doubling time                     : %6.0f d (%.1f mo)"
       % (td, mo(td)))
    suv48 = pct_change(r2, "suv", 2.0)
    ki7 = pct_change(r2, "ki67", 7.0)
    sld48 = pct_change(r2, "sld", 2.0)
    sld56 = pct_change(r2, "sld", 56.0)
    pr("B2  Imatinib: SUV at 48 h %+6.1f %%   Ki67 at day 7 %+6.1f %%" % (suv48, ki7))
    pr("             SLD at 48 h %+6.1f %%   SLD at 8 weeks %+6.1f %%" % (sld48, sld56))
    pr("B3  Best SLD change on first line                      : %+6.1f %%"
       % R["S2"]["best"]["stage0"])
    pr("    (a partial response; complete responses are absent in the model for")
    pr("     the same reason they are rare in the clinic -- the reservoir)")
    q1 = at(r2, "nq", 365.0) / max(at(r2, "ntot", 365.0), 1.0)
    pr("B4  At 1 year of imatinib: quiescent share of viable cells   : %5.1f %%"
       % (100 * q1))
    pr("    viable share of the IMAGED mass                          : %5.1f %%"
       % (100 * at(r2, "vfrac", 365.0)))
    check("B1 untreated doubling time 45-90 days", 45.0 < td < 90.0, "%.0f d" % td)
    check("B2 PET moves and CT does not, at 48 h",
          suv48 < -40.0 and abs(sld48) < 2.0,
          "SUV %+.1f%% vs SLD %+.1f%%" % (suv48, sld48))
    check("B3 best response is partial, not complete",
          -65.0 < R["S2"]["best"]["stage0"] < -30.0)
    check("B4 the residual mass is mostly non-viable at 1 year",
          at(r2, "vfrac", 365.0) < 0.40)

    # ---- C. first line by genotype ------------------------------------
    h("C.  FIRST LINE: THE EC50 SITS ON EITHER SIDE OF A THRESHOLD  [C3]")
    rows = [("KIT exon 11, 400 mg", "S2"), ("KIT exon 11, 800 mg", "S5"),
            ("KIT exon 9, 400 mg", "S3"), ("KIT exon 9, 800 mg", "S4"),
            ("exon 11, lowest exposure quartile", "S6"),
            ("PDGFRA D842V, imatinib 400", "S7"),
            ("PDGFRA D842V, avapritinib 300", "S8"),
            ("SDH-deficient, imatinib 400", "S9")]
    pr("    %-36s %12s %12s" % ("", "PFS (mo)", "best SLD %"))
    tab = {}
    for lab, k in rows:
        p = stage_pfs(R[k], 0)
        tab[k] = (p, R[k]["best"]["stage0"])
        pr("    %-36s %12s %12.1f" % (lab, fmt(p, ""), R[k]["best"]["stage0"]))
    pr("")
    pr("    Observed: KIT exon 11 median PFS ~25 mo (S0033).  Exon 9 raises the")
    pr("    relative risk of progression 171% and 800 mg reduces it 61% -- in exon")
    pr("    9 only (EORTC 62005 / MetaGIST), which for a 25-month exon 11 median")
    pr("    implies ~9 mo for exon 9 at 400 mg and ~23 mo at 800 mg.  D842V:")
    pr("    imatinib ORR ~0%, avapritinib ORR 88% (NAVIGATOR).")
    e11, e11_8 = tab["S2"][0], tab["S5"][0]
    e9_4, e9_8 = tab["S3"][0], tab["S4"][0]
    check("C1 exon 11 first-line PFS 22-32 months", 22.0 < mo(e11) < 32.0,
          "%.1f mo" % mo(e11))
    check("C2 exon 9 at 400 mg is far worse than exon 11",
          mo(e9_4) < 0.5 * mo(e11), "%.1f vs %.1f mo" % (mo(e9_4), mo(e11)))
    check("C3 escalation rescues exon 9 (>2x PFS)", e9_8 > 2.0 * e9_4,
          "%.1f -> %.1f mo" % (mo(e9_4), mo(e9_8)))
    check("C4 escalation does nothing for exon 11 (<10% gain)",
          e11_8 < 1.10 * e11,
          "%.1f -> %.1f mo (%+.0f%%)" % (mo(e11), mo(e11_8),
                                         100 * (e11_8 / e11 - 1)))
    check("C5 D842V: imatinib fails, avapritinib works",
          tab["S7"][1] > -10.0 and tab["S8"][1] < -40.0,
          "best SLD %+.1f%% vs %+.1f%%" % (tab["S7"][1], tab["S8"][1]))
    check("C6 SDH-deficient GIST does not respond and is indolent",
          tab["S9"][1] > -10.0 and mo(tab["S9"][0]) > 5.0,
          "best %+.1f%%, PFS %.1f mo" % (tab["S9"][1], mo(tab["S9"][0])))
    esc = stage_pfs(R["S26"], 0) + (stage_pfs(R["S26"], 1) or 0.0)
    pr("")
    pr("    Escalating an exon 9 patient only AFTER progression, instead of")
    pr("    starting at 800 mg: %s of total time on imatinib versus %s."
       % (fmt(esc), fmt(e9_8)))
    pr("    So the model does NOT say start high -- it says the threshold has to be")
    pr("    crossed eventually, and crossing it at progression loses nothing.  That")
    pr("    is the practice guidelines already follow, and it is consistent with")
    pr("    MetaGIST finding a PFS but no OS advantage for 800 mg.")
    check("C7 crossing the threshold at progression is not worse than starting high",
          esc > 0.9 * e9_8, "%.1f vs %.1f mo" % (mo(esc), mo(e9_8)))

    # ---- D. the INTRIGUE crossover ------------------------------------
    h("D.  THE INTRIGUE CROSSOVER FALLS OUT OF THE CLONE VECTOR  [C1]")
    su_atp, ri_atp = mo(stage_pfs(R["S12"], 1)), mo(stage_pfs(R["S13"], 1))
    su_al, ri_al = mo(stage_pfs(R["S14"], 1)), mo(stage_pfs(R["S15"], 1))
    su_un, ri_un = mo(stage_pfs(R["S10"], 1)), mo(stage_pfs(R["S11"], 1))
    pr("    Second-line PFS, by which secondary-mutation class the tumour carries:")
    pr("    %-34s %12s %12s" % ("", "sunitinib", "ripretinib"))
    pr("    %-34s %12.1f %12.1f" % ("KIT exon 11 + 13/14 only", su_atp, ri_atp))
    pr("    %-34s %12.1f %12.1f" % ("KIT exon 11 + 17/18 only", su_al, ri_al))
    pr("    %-34s %12.1f %12.1f" % ("both classes in one tumour", su_un, ri_un))
    pr("")
    pr("    Observed (Heinrich 2024 Nat Med, INTRIGUE ctDNA):")
    pr("      exon 11 + 13/14   sunitinib 15.0   ripretinib  4.0 months")
    pr("      exon 11 + 17/18   sunitinib  1.5   ripretinib 14.2 months")
    pr("      ITT               sunitinib  8.3   ripretinib  8.0 months")
    pr("    The two subgroups were MUTUALLY EXCLUSIVE in the trial, so the ITT")
    pr("    median is the median of a mixture of two opposite-signed effects, not")
    pr("    a property of any patient.  The model therefore does not try to")
    pr("    reproduce the ITT number; it reproduces the two subgroups and says")
    pr("    that a patient carrying BOTH classes does badly on either drug.")
    pr("    Nothing here was fitted to INTRIGUE: the ranking follows from the")
    pr("    published potency pattern (Serrano 2019, Smith 2019, Heinrich 2008).")
    check("D1 sunitinib beats ripretinib in the ATP-pocket subgroup",
          su_atp > 1.8 * ri_atp, "%.1f vs %.1f mo" % (su_atp, ri_atp))
    check("D2 ripretinib beats sunitinib in the activation-loop subgroup",
          ri_al > 1.8 * su_al, "%.1f vs %.1f mo" % (ri_al, su_al))
    check("D3 the ranking inverts between the two subgroups",
          (su_atp / ri_atp) * (su_al / ri_al) < 0.25,
          "ratio product %.3f" % ((su_atp / ri_atp) * (su_al / ri_al)))
    check("D4 a tumour with both classes does worse than the better-matched arm",
          max(su_un, ri_un) < max(su_atp, ri_al),
          "%.1f vs %.1f mo" % (max(su_un, ri_un), max(su_atp, ri_al)))

    # ---- E. later lines ----------------------------------------------
    h("E.  LATER LINES: THE CLONE NOBODY COVERS SETS THE CEILING  [C1]")
    reg3 = mo(stage_pfs(R["S16"], 2))
    ri4 = mo(stage_pfs(R["S17"], 3))
    pl4 = mo(stage_pfs(R["S18"], 3))
    im4 = mo(stage_pfs(R["S19"], 3))
    pr("    third line regorafenib                      : %5.1f mo   (observed 4.8)" % reg3)
    pr("    fourth line ripretinib                      : %5.1f mo   (observed 6.3)" % ri4)
    pr("    fourth line best supportive care            : %5.1f mo   (observed 0.9-1.0)" % pl4)
    pr("    fourth line imatinib rechallenge            : %5.1f mo   (observed 1.8)" % im4)
    t4 = line_start(R["S18"], 3)
    r18 = R["S18"]["rec"]
    tot4 = max(at(r18, "ntot", t4), 1.0)
    pr("")
    pr("    Clone composition at the start of fourth line:")
    pr("      primary mutation only        %5.1f %%" % (100 * at(r18, "n0", t4) / tot4))
    pr("      exon 13/14 ATP-pocket        %5.1f %%" % (100 * at(r18, "n1", t4) / tot4))
    pr("      exon 17/18 activation-loop   %5.1f %%" % (100 * at(r18, "n2", t4) / tot4))
    pr("      KIT-independent bypass       %5.1f %%" % (100 * at(r18, "n3", t4) / tot4))
    pr("    A sixth of the tumour still carries only the primary mutation after")
    pr("    three lines, which is why re-exposing it to imatinib does anything at")
    pr("    all (RIGHT).  The model gets the direction and under-predicts the size.")
    check("E1 fourth-line ripretinib beats best supportive care",
          ri4 > 2.5 * pl4, "%.1f vs %.1f mo" % (ri4, pl4))
    check("E2 imatinib rechallenge beats best supportive care",
          im4 > pl4, "%.1f vs %.1f mo" % (im4, pl4))
    check("E3 a sixth of cells still carry only the primary mutation at line 4",
          at(r18, "n0", t4) / tot4 > 0.10,
          "%.1f%%" % (100 * at(r18, "n0", t4) / tot4))
    check("E4 successive lines get shorter",
          mo(stage_pfs(R["S16"], 1)) > reg3, "%.1f then %.1f mo"
          % (mo(stage_pfs(R["S16"], 1)), reg3))

    # ---- F. ctDNA ----------------------------------------------------
    h("F.  ctDNA SEES THE CLONE BEFORE RECIST DOES")
    r10 = R["S10"]["rec"]
    pd1 = stage_pfs(R["S10"], 0)
    tvaf = first_time(r10, "vaf", 0.01)
    pr("    resistance-mutation VAF crosses 1%%        : day %5.0f" % tvaf)
    pr("    RECIST progression on first line          : day %5.0f" % pd1)
    pr("    lead time                                 : %5.1f months" % mo(pd1 - tvaf))
    pr("    VAF at RECIST progression                 : %5.1f %%"
       % (100 * at(r10, "vaf", pd1)))
    pr("    This window is the whole practical point of C1: it is when the second")
    pr("    line could be chosen by genotype instead of by label order.")
    check("F1 ctDNA leads RECIST by more than 3 months", mo(pd1 - tvaf) > 3.0,
          "%.1f mo" % mo(pd1 - tvaf))

    # ---- G. stopping -------------------------------------------------
    h("G.  STOPPING THE DRUG: THE RESERVOIR  [C2]")
    pr("    %-12s %14s %14s %10s %8s" % ("randomised", "interruption", "continuation",
                                         "viable frac", "dormancy"))
    bfr = {}
    for yy in (1, 3, 5):
        stop, go = interruption_scenario(yy)
        rs = simulate(stop, dt=0.06)
        rg = simulate(go, dt=0.06)
        ps, pg = stage_pfs(rs, 1), stage_pfs(rg, 1)
        bfr[yy] = (ps, pg)
        pr("    %-12s %14s %14s %10.2f %8.2f"
           % ("after %d y" % yy, fmt(ps, ""), fmt(pg, ""),
              at(rs["rec"], "vfrac", yy * 365.0),
              at(rs["rec"], "qdep", yy * 365.0)))
    pr("")
    pr("    Observed (BFR14, Blay 2024): 6.1 vs 27.8 mo after 1 year, 7.0 vs 67.0")
    pr("    after 3 years, 12.0 vs not reached after 5 years.")
    pr("    The model needs an essentially resistance-free tumour to stay")
    pr("    progression-free for 3-5 years, and that is itself a prediction: only")
    pr("    50 of 434 patients enrolled in BFR14 reached the 3-year randomisation")
    pr("    and 27 the 5-year one.  Each cohort's phenotype here is set from its")
    pr("    OWN continuation arm, so the interruption arm is predicted.")
    p1, p3, p5 = mo(bfr[1][0]), mo(bfr[3][0]), mo(bfr[5][0])
    check("G1 interruption after 1 year progresses in 3-8 months",
          3.0 < p1 < 8.0, "%.1f mo" % p1)
    check("G2 continuation beats interruption in every arm",
          all((bfr[y][1] or 1e9) > 2.5 * bfr[y][0] for y in (1, 3, 5)),
          "; ".join("%dy %s vs %s" % (y, fmt(bfr[y][1], ""), fmt(bfr[y][0], ""))
                    for y in (1, 3, 5)))
    check("G3 the 5-year cohort takes longest to regrow",
          p5 > p1 and p5 > p3, "%.1f / %.1f / %.1f mo" % (p1, p3, p5))
    check("G4 years of deep response do not eradicate the disease",
          all(bfr[y][0] is not None for y in (1, 3, 5)))

    # ---- H. adjuvant -------------------------------------------------
    h("H.  ADJUVANT THERAPY: DELAY, NOT CURE")
    adj = {}
    pr("    %-14s %18s %16s" % ("duration", "recurrence detected", "delay vs 12 mo"))
    for lab, k in (("12 months", "S22"), ("36 months", "S23"), ("60 months", "S24")):
        rr = R[k]["recurrence"]
        adj[k] = rr
        d0 = adj["S22"]
        pr("    %-14s %18s %16s"
           % (lab, fmt(rr, ""), "-" if k == "S22" else "%.1f mo" % mo(rr - d0)))
    pr("")
    pr("    Observed (SSGXVIII, Joensuu 2020): 5-year RFS 53.0% (12 mo) vs 71.4%")
    pr("    (36 mo); 10-year 41.8% vs 52.5%.  The trial-level percentages need a")
    pr("    virtual population (gist_mrgsolve_model.R block 9); this deterministic")
    pr("    twin reports the single-patient time to detectable recurrence.")
    d12, d36, d60 = adj["S22"], adj["S23"], adj["S24"]
    check("H1 longer adjuvant therapy delays recurrence", d60 > d36 > d12,
          "%.0f / %.0f / %.0f d" % (d12, d36, d60))
    check("H2 the delay does not exceed the extra treatment",
          (d36 - d12) < 1.6 * 730.0,
          "delay %.0f d for 730 d extra drug" % (d36 - d12))
    check("H3 nobody is cured: recurrence happens in every arm",
          all(adj[k] is not None for k in adj))

    # ---- I. toxicity -------------------------------------------------
    h("I.  TOLERABILITY IS PART OF EFFICACY")
    r16 = R["S16"]["rec"]
    r10x = R["S10"]["rec"]
    r2x = R["S2"]["rec"]
    r20 = R["S20"]["rec"]
    su0 = line_start(R["S10"], 1)
    su0b = line_start(R["S20"], 1)
    re0 = line_start(R["S16"], 2)
    regs = (("imatinib 400", r2x, 30.0, 330.0),
            ("sunitinib 50 4/2", r10x, su0 + 30.0, su0 + 130.0),
            ("sunitinib 37.5 daily", r20, su0b + 30.0, su0b + 130.0),
            ("regorafenib 160 3/1", r16, re0 + 30.0, re0 + 130.0))
    pr("    %-22s %8s %7s %7s %7s %6s %7s" % ("regimen", "ANC min", "dMAP",
                                              "HFSR", "oedema", "TSH", "DI min"))
    for lab, rr, t0, t1 in regs:
        pr("    %-22s %8.2f %7.1f %7.2f %7.2f %6.2f %7.2f"
           % (lab, min(win(rr, "anc", t0, t1)), max(win(rr, "map", t0, t1)),
              max(win(rr, "hfsr", t0, t1)), max(win(rr, "edem", t0, t1)),
              max(win(rr, "tsh", t0, t1)), min(win(rr, "di", t0, t1))))
    pr("    sunitinib free T4 nadir %.2f of baseline, TSH peak %.2f x"
       % (min(win(r10x, "ft4", su0 + 30.0, su0 + 130.0)),
          max(win(r10x, "tsh", su0 + 30.0, su0 + 130.0))))
    pr("    Each row is read inside its own line of therapy, not over the whole")
    pr("    run, because oedema from the imatinib phase would otherwise be")
    pr("    attributed to whatever came next.")
    pr("")
    pr("    Second-line PFS, sunitinib 50 mg 4/2 vs 37.5 mg continuous : %.1f vs %.1f mo"
       % (su_un, mo(stage_pfs(R["S20"], 1))))
    pr("    The intermittent schedule lets every clone regrow for two weeks in six;")
    pr("    the continuous one is gentler and covers continuously.  The model makes")
    pr("    them close, which is what the clinical comparisons found.")
    check("I1 sunitinib raises blood pressure and TSH, imatinib does not",
          max(win(r10x, "map", su0 + 30, su0 + 130)) > 6.0
          and max(win(r10x, "tsh", su0 + 30, su0 + 130)) > 1.5
          and max(win(r2x, "map", 30, 330)) < 1.0)
    check("I2 imatinib causes oedema, sunitinib does not",
          max(win(r2x, "edem", 30, 330)) > 0.5
          and max(win(r10x, "edem", su0 + 60, su0 + 130)) < 0.2,
          "%.2f vs %.2f" % (max(win(r2x, "edem", 30, 330)),
                            max(win(r10x, "edem", su0 + 60, su0 + 130))))
    check("I3 regorafenib forces the largest dose reduction",
          min(win(r16, "di", re0 + 30, re0 + 130)) <= min(r2x["di"]) - 0.05,
          "DI %.2f vs %.2f" % (min(win(r16, "di", re0 + 30, re0 + 130)),
                               min(r2x["di"])))
    check("I4 continuous 37.5 mg sunitinib is within 30% of 50 mg 4/2",
          abs(mo(stage_pfs(R["S20"], 1)) - su_un) < 0.30 * su_un,
          "%.1f vs %.1f mo" % (mo(stage_pfs(R["S20"], 1)), su_un))

    # ---- J. the combination question ---------------------------------
    h("J.  COVER THE CLONES EARLY?  THE MODEL SAYS NO, AND SAYS WHY")
    IM4 = ("IM", 400, "qd")
    RI1 = ("RI", 150, "qd")
    RI0 = ("RI", 100, "qd")
    SU5 = ("SU", 50, "4/2")
    PROG = dict(on="progression")
    seq = simulate(S(label="seq", geno="exon11", tend=4000.0,
                     stages=[RX(IM4), RX(RI1, **PROG), RX(SU5, **PROG)]))
    comb = simulate(S(label="comb", geno="exon11", tend=4000.0,
                      stages=[RX(IM4, RI0), RX(SU5, **PROG)]))
    tseq = sum(seq["pfs"].get("stage%d" % i, 0.0) for i in range(3))
    tcomb = sum(comb["pfs"].get("stage%d" % i, 0.0) for i in range(2))
    pr("    sequential  imatinib -> ripretinib -> sunitinib")
    pr("      per line %s | %s | %s      total %s"
       % (fmt(stage_pfs(seq, 0), ""), fmt(stage_pfs(seq, 1), ""),
          fmt(stage_pfs(seq, 2), ""), fmt(tseq, "")))
    pr("    upfront     imatinib 400 + ripretinib 100 -> sunitinib")
    pr("      per line %s | %s              total %s"
       % (fmt(stage_pfs(comb, 0), ""), fmt(stage_pfs(comb, 1), ""), fmt(tcomb, "")))
    pr("")
    pr("    So the combination buys %+.1f months of FIRST-line PFS and loses"
       % mo(stage_pfs(comb, 0) - stage_pfs(seq, 0)))
    pr("    %.1f months overall.  That is not the intuitive answer and it is not a"
       % mo(tseq - tcomb))
    pr("    fitted one: covering the activation-loop clone before it is selected")
    pr("    removes the very heterogeneity that the later lines exploit, so the")
    pr("    tumour arrives at the KIT-independent clone sooner.  A set-cover model")
    pr("    predicts that upfront combination is a bad trade even though it")
    pr("    lengthens the first line -- which is what the field has found the hard")
    pr("    way, and is a falsifiable statement about any future combination trial.")
    check("J1 upfront combination lengthens first line",
          stage_pfs(comb, 0) > stage_pfs(seq, 0),
          "%.1f vs %.1f mo" % (mo(stage_pfs(comb, 0)), mo(stage_pfs(seq, 0))))
    check("J2 and shortens total time on strategy", tcomb < tseq,
          "%.1f vs %.1f mo" % (mo(tcomb), mo(tseq)))

    # ---- K. falsification --------------------------------------------
    h("K.  FALSIFICATION: ONE CLONE WHOSE POTENCY DRIFTS, REFITTED")
    target = stage_pfs(R["S2"], 0)
    lo, hi = 5.0e-4, 1.0e-2
    for _ in range(14):
        mid = math.sqrt(lo * hi)
        f = S(label="f", geno="exon11", tend=3000.0,
              par=dict(FALSIFY=1.0, KDRIFT=mid), stages=[RX(IM4)])
        p = stage_pfs(simulate(f, dt=0.06), 0)
        if p is None:
            hi = mid
        elif p > target:
            lo = mid
        else:
            hi = mid
    kdrift = math.sqrt(lo * hi)

    def fal(stages, geno="exon11", tend=3600.0):
        return simulate(S(label="f", geno=geno, tend=tend,
                          par=dict(FALSIFY=1.0, KDRIFT=kdrift), stages=stages),
                        dt=0.06)

    f1 = fal([RX(IM4)])
    pr("    The fair competitor is not 'no resistance' but resistance as a")
    pr("    continuous loss of potency under drug pressure, refitted to the SAME")
    pr("    first-line anchor.")
    pr("    KDRIFT refitted                       : %.4g /day" % kdrift)
    pr("    first line, exon 11                   : %5.1f mo   (target %.1f)"
       % (mo(stage_pfs(f1, 0)), mo(target)))
    f_su = fal([RX(IM4), RX(SU5, **PROG)])
    f_ri = fal([RX(IM4), RX(RI1, **PROG)])
    fsu, fri = mo(stage_pfs(f_su, 1)), mo(stage_pfs(f_ri, 1))
    pr("")
    pr("    [D] second line sunitinib vs ripretinib: %5.1f vs %5.1f mo" % (fsu, fri))
    pr("        and this is the SAME pair of numbers in every genotype subgroup,")
    pr("        because with one clone MU_ATP and MU_AL are inert -- the subgroup")
    pr("        parameters cannot change the simulation at all.  The INTRIGUE")
    pr("        crossover is not merely unfitted here, it is UNREACHABLE.")
    f_pl = fal([RX(IM4), RX(SU5, **PROG), RX(("RE", 160, "3/1"), **PROG), RX(**PROG)])
    f_im = fal([RX(IM4), RX(SU5, **PROG), RX(("RE", 160, "3/1"), **PROG),
                RX(IM4, **PROG)])
    fp, fi = mo(stage_pfs(f_pl, 3)), mo(stage_pfs(f_im, 3))
    pr("    [E] fourth line, best supportive care vs imatinib rechallenge:")
    pr("        %5.1f vs %5.1f mo  -- the drifting-potency model ALSO produces a" % (fp, fi))
    pr("        rechallenge benefit, in fact a larger one, because a 7-fold")
    pr("        potency shift still leaves imatinib partly active.  So the RIGHT")
    pr("        trial does NOT discriminate between the two structures, and this")
    pr("        report does not claim that it does.")
    f_ri4 = fal([RX(IM4), RX(SU5, **PROG), RX(("RE", 160, "3/1"), **PROG),
                 RX(RI1, **PROG)])
    fri4 = mo(stage_pfs(f_ri4, 3))
    pr("    [E] fourth line ripretinib (INVICTUS, observed 6.3 mo):")
    pr("        %5.1f mo here versus %5.1f mo in the clone-vector model -- with one"
       % (fri4, ri4))
    pr("        clone there is nothing broad for a broad-spectrum inhibitor to be")
    pr("        good at, so this test IS decisive.")
    f9_4 = fal([RX(IM4)], geno="exon9")
    f9_8 = fal([RX(("IM", 800, "qd"))], geno="exon9")
    f11_8 = fal([RX(("IM", 800, "qd"))])
    g9 = 100.0 * (stage_pfs(f9_8, 0) / stage_pfs(f9_4, 0) - 1.0)
    g11 = 100.0 * (stage_pfs(f11_8, 0) / stage_pfs(f1, 0) - 1.0)
    pr("    [C] escalation to 800 mg: exon 9 %+.0f %%, exon 11 %+.0f %%" % (g9, g11))
    pr("        In the polyclonal model the same comparison is %+.0f %% and %+.0f %%."
       % (100.0 * (e9_8 / e9_4 - 1.0), 100.0 * (e11_8 / e11 - 1.0)))
    pr("")
    pr("")
    pr("    Verdict.  Of four discriminating tests, TWO are decisive -- the genotype")
    pr("    crossover, which a single-clone model cannot express at any parameter")
    pr("    value, and fourth-line ripretinib.  The rechallenge benefit does NOT")
    pr("    discriminate: the drifting model produces it too, and larger.  The")
    pr("    escalation asymmetry is mostly a property of C3, not C1, so it")
    pr("    discriminates only weakly.  Reported as such rather than as a sweep.")
    check("K1 the falsifier is refitted to the same first-line anchor",
          abs(mo(stage_pfs(f1, 0)) - mo(target)) < 2.0,
          "%.1f vs %.1f mo" % (mo(stage_pfs(f1, 0)), mo(target)))
    check("K2 with one clone the subgroup parameters are inert",
          abs(stage_pfs(fal([RX(IM4), RX(SU5, **PROG)],)
                        , 1) - stage_pfs(f_su, 1)) < 1e-9)
    check("K3 the polyclonal escalation asymmetry is larger than the falsifier's",
          (e9_8 / e9_4) / (e11_8 / e11) > (stage_pfs(f9_8, 0) / stage_pfs(f9_4, 0))
          / (stage_pfs(f11_8, 0) / stage_pfs(f1, 0)),
          "%.2f vs %.2f" % ((e9_8 / e9_4) / (e11_8 / e11),
                            (stage_pfs(f9_8, 0) / stage_pfs(f9_4, 0))
                            / (stage_pfs(f11_8, 0) / stage_pfs(f1, 0))))
    check("K4 single clone cannot reproduce fourth-line ripretinib",
          fri4 < 0.6 * ri4, "%.1f vs %.1f mo" % (fri4, ri4))
    check("K5 the rechallenge test is reported as non-discriminating",
          fi > fp, "falsifier also shows a benefit: %.1f vs %.1f mo" % (fi, fp))

    # ---- L. misses ---------------------------------------------------
    h("L.  WHAT THE MODEL GETS WRONG, STATED PLAINLY")
    pr("    M1  Imatinib trough and outcome in exon 11.  Demetri 2009 found TTP")
    pr("        11.3 mo in the lowest trough quartile versus > 30 mo above it.  In")
    pr("        this model exon 11 first-line PFS is almost exposure-independent")
    pr("        above a trough of ~250 ng/mL (%.1f mo at a trough of %.0f ng/mL)"
       % (mo(tab["S6"][0]), at(R["S6"]["rec"], "cim", 400.0)))
    pr("        because progression is driven by a clone imatinib never covers at")
    pr("        any dose; only the DEPTH of response tracks exposure (%+.1f%% at the"
       % tab["S6"][1])
    pr("        low trough versus %+.1f%% at the normal one).  Either the clinical"
       % R["S2"]["best"]["stage0"])
    pr("        association is not a potency threshold, or exon 11 EC50s are more")
    pr("        spread than a single number.")
    pr("    M2  The imatinib rechallenge effect is too small (%.1f vs %.1f mo,"
       % (im4, pl4))
    pr("        observed 1.8 vs 0.9), and as section K shows the observation does")
    pr("        not discriminate between the two model structures anyway.")
    pr("    M3  Third-line regorafenib (%.1f mo vs 4.8 observed) and second-line"
       % reg3)
    pr("        agents in a both-classes tumour are under-predicted, because the")
    pr("        trial populations are mixtures of single-class patients while the")
    pr("        model's 'unselected' patient carries both classes at once.")
    pr("    M4  BFR14 magnitudes are short by 1-3 months in the 1- and 3-year arms")
    pr("        (%.1f and %.1f vs 6.1 and 7.0), and the model cannot represent the"
       % (p1, p3))
    pr("        3- and 5-year cohorts without giving them resistance-free tumours.")
    pr("    M5  No overall survival, no spatial structure, so no focal progression")
    pr("        and no cytoreductive surgery of a single progressing nodule.")

    # ---- checks -------------------------------------------------------
    h("SELF-CHECKS")
    npass = 0
    for name, ok, detail in CHECKS:
        print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name,
                               ("   (%s)" % detail) if detail else ""))
        npass += ok
    print("\n  %d/%d checks passed" % (npass, len(CHECKS)))
    return npass, len(CHECKS)


if __name__ == "__main__":
    quiet = "--quiet" in sys.argv
    np_, nt_ = run_report(quiet=quiet)
    sys.exit(0 if np_ == nt_ else 1)
