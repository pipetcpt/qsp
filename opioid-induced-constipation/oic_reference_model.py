#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oic_reference_model.py
======================
Opioid-Induced Constipation (OIC) — independent Python/scipy implementation of
the 51-ODE QSP model that is shipped in this directory as
`oic_mrgsolve_model.R`.

WHY THIS FILE EXISTS
--------------------
The container this model was built in has **no R runtime**.  Committing an
unexecuted ODE model would be dishonest.  So every equation was implemented
here, integrated for real, and the R file was transcribed from *this* — not the
other way round.  Every number quoted in `README.md` is produced by running
this file (`python3 oic_reference_model.py`), and its stdout is committed
verbatim as `oic_reference_output.txt`.

THE ONE STRUCTURAL CLAIM
------------------------
A PAMORA's therapeutic window is **not** a potency.  It is the ratio of two
receptor occupancies computed at the *same* plasma concentration, and that
ratio is set by P-glycoprotein at the blood-brain barrier:

      SI  =  OCC_ant(enteric plexus)  /  OCC_ant(CNS)

Because the enteric plexus is outside the BBB and the CNS is not, and because
both sites hold the *same* receptor competing with the *same* agonist, SI
collapses (in the low-occupancy limit) to 1/Kp_uu of the antagonist, corrected
by the agonist's own brain partitioning.  Nothing about "selectivity for
peripheral mu receptors" needs to be asserted; it is division.

What follows arithmetically, and is computed in analyses A-J below:

  A. SI ranks naloxone << naloxegol < methylnaltrexone < naldemedine, spanning
     149-fold from Kp_uu alone -- while the MOR binding Ki spans only 25-fold,
     in the OPPOSITE order.
  B. A P-gp inhibitor does not merely raise AUC, it *rotates the ratio*.  The
     label's "reduce the dose" advice repairs a CYP3A4 interaction and does not
     repair a P-gp one: it surrenders efficacy without recovering the margin.

TWO OF THE AUTHOR'S PRIOR HYPOTHESES WERE REFUTED BY RUNNING THE MODEL, AND
BOTH ARE REPORTED AS SUCH RATHER THAN QUIETLY REPLACED:

  C. Predicted: escalating the opioid dose to hold analgesia would progressively
     WORSEN constipation.  It does not.  Enteric transduction is already
     saturated at an ordinary analgesic dose, so SBM/week is flat from 60 to
     200 mg/day (1.49 -> 1.34).  What the tolerance asymmetry actually costs the
     patient is on the other side of the ledger: the escalation buys no
     analgesia (pain pinned at 5.04) while the gut keeps every milligram.  This
     also explains an epidemiological fact a dose-response account cannot --
     that OIC prevalence correlates poorly with opioid dose.
  D. Predicted: the transit-hydration feedback loop is the amplifier that makes
     a PAMORA's clinical effect exceed the receptor occupancy it buys.  It is
     not.  Disabling that loop entirely moves the drug fold-response from 2.89
     to 2.91.  The amplification is carried by the ANORECTAL brakes -- the
     rectal urge threshold (fold 1.60) and sphincter tone (1.88) -- i.e. it is
     the product of several independent brakes released together by one
     receptor, dominated by the END of the colon.  That predicts responders
     split by anorectal rather than transit phenotype, which is testable.

  H. Cmax-vs-Cavg is invisible to a weekly-frequency endpoint (splitting the
     dose mildly HELPS both drugs) and visible only in the acute 4-h laxation
     endpoint -- which is precisely the endpoint methylnaltrexone was licensed
     on and naloxegol and naldemedine were not.

Units: time h, drug amounts mg (peptides ug), colonic solids g, water mL,
osmoles mmol, receptor concentrations nM.

Run:  python3 oic_reference_model.py
"""

import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

np.seterr(all="ignore")

# ---------------------------------------------------------------------------
# 0.  STATE VECTOR  (51 ODEs)
# ---------------------------------------------------------------------------
NAMES = [
    # --- opioid PK (5) ---
    "AOP_DEP", "AOP_CEN", "AOP_PER", "COP_BR", "AOP_LUM",
    # --- antagonist (PAMORA) PK (5) ---
    "APAM_DEP", "APAM_CEN", "APAM_PER", "CPAM_BR", "APAM_LUM",
    # --- adjunct laxatives / prokinetic (6) ---
    "APEG_LUM", "ALAC_LUM", "ALUB_LUM", "ALIN_LUM", "APRO_CEN", "APRO_PER",
    # --- receptor trafficking (3) ---
    "RG_AV", "RC_AV", "ARR",
    # --- enteric signalling (5) ---
    "CAMP", "ACH", "NOVIP", "CGMP", "CLC2",
    # --- motility (2) ---
    "HAPC", "TONE",
    # --- colonic content: 4 segments x (solids, water, osmoles) (12) ---
    "S1", "S2", "S3", "S4",
    "W1", "W2", "W3", "W4",
    "O1", "O2", "O3", "O4",
    # --- symptoms (5) ---
    "BSFS", "STRAIN", "DIST", "PACSYM", "QOL",
    # --- CNS / analgesia / withdrawal (4) ---
    "PAIN", "TOLA", "WD", "NAUSEA",
    # --- counters & safety (4) ---
    "CUM_SBM", "CUM_CSBM", "CUM_RESC", "IMPACT",
]
IX = {n: i for i, n in enumerate(NAMES)}
NST = len(NAMES)
assert NST == 51, NST


# ---------------------------------------------------------------------------
# 1.  PARAMETERS
# ---------------------------------------------------------------------------
def base_params():
    p = {}

    # ===================== OPIOID (index drug = oxycodone) =================
    # Oxycodone: MW 315.4, F~0.75, CL ~45 L/h, V1 90 L, V2 120 L, Q 40 L/h,
    # fu ~0.55.  Kp_uu ~3 -- oxycodone is one of the few opioids with ACTIVE
    # influx at the BBB (Bostrom/Hammarlund-Udenaes).  That is not decoration:
    # it is the term that makes an antagonist's brain competition *harder*,
    # i.e. it improves PAMORA safety (see analysis A, column "SI").
    p["MWOP"] = 315.4
    p["FOP"] = 0.75
    p["KAOP"] = 1.2
    p["CLOP"] = 45.0
    p["V1OP"] = 90.0
    p["V2OP"] = 120.0
    p["QOP"] = 40.0
    p["FUOP"] = 0.55
    p["KPUUOP"] = 3.0
    p["KEOOP"] = 0.60
    p["KIOP"] = 30.0          # nM, MOR Ki (unbound)
    p["FLUMOP"] = 0.25        # fraction of an ORAL opioid dose reaching colonic lumen
    p["KLUMOUT"] = 0.15       # /h luminal opioid loss (absorption + degradation)

    # ===================== ANTAGONIST (set per drug below) =================
    p.update(ANTAGONISTS["naloxegol"])

    # ===================== PROKINETIC (prucalopride, 5-HT4) ================
    p["MWPRO"] = 367.9
    p["FPRO"] = 0.90
    p["KAPRO"] = 1.0
    p["CLPRO"] = 18.0
    p["V1PRO"] = 320.0
    p["V2PRO"] = 300.0
    p["QPRO"] = 25.0
    p["EC50PRO"] = 8.0        # nM at enteric 5-HT4
    p["EMAXPRO"] = 0.55       # max fractional rise of enteric cAMP

    # ===================== GUT-WALL / LUMEN GEOMETRY =======================
    p["QPORT"] = 60.0         # L/h portal blood flow (pre-systemic gut-wall term)
    p["VLUM"] = 0.35          # L effective colonic luminal volume for plexus access
    # DEFECT #5: at FLUMACC 0.04 the LUMINAL opioid term delivered ~1500 nM to
    # the plexus against 19 nM unbound in plasma -- an 80-fold swamp that
    # saturated enteric occupancy and destroyed all opioid DOSE-dependence
    # (oxycodone 30 and 60 mg/day gave 2.60 and 2.58 SBM/wk).  Oral opioids do
    # expose the plexus above plasma, but not by two orders of magnitude.
    p["FLUMACC"] = 0.002      # fraction of luminal conc that reaches myenteric plexus
    # DEFECT #7: the luminal->plexus access factor was shared by every drug,
    # which made ORAL methylnaltrexone 450 mg a miracle drug (SBM 11.0/wk, a
    # complete normalisation the trials do not show).  A permanently charged
    # quaternary ammonium cannot cross the mucosa the way a neutral lipophilic
    # molecule can -- that inability IS why the drug is peripherally
    # restricted.  It must therefore be penalised on the LUMINAL route too,
    # not only at the BBB.  Carried per drug as FLUMACCX.
    p["FLUMACCX"] = 1.0       # per-drug multiplier on luminal mucosal access

    # ===================== RECEPTOR TRAFFICKING ============================
    # THE tolerance asymmetry.  KDESC/KDESG = 7.5.  This single ratio is why
    # analgesic tolerance develops and constipation tolerance does not.
    p["KINRC"] = 0.0060       # /h CNS MOR resensitisation
    p["KDESC"] = 0.0045       # /h CNS MOR desensitisation (per unit occupancy)
    p["KINRG"] = 0.0060       # /h enteric MOR resensitisation
    p["KDESG"] = 0.00060      # /h enteric MOR desensitisation <-- 7.5x smaller
    p["KONARR"] = 0.012
    p["KOFFARR"] = 0.010
    p["BARR"] = 1.20          # beta-arrestin amplification of CNS desensitisation

    # ===================== TRANSDUCTION (operational agonism) ==============
    # Occupancy -> effect.  OCC50G < OCC50C: the gut arm turns on at LOWER
    # occupancy than the analgesic arm, which is why constipation appears from
    # the first dose while analgesia needs titration.
    p["OCC50G"] = 0.30
    p["HG"] = 1.40
    p["OCC50C"] = 0.30
    p["HC"] = 2.00

    # ===================== ENTERIC SIGNALLING ==============================
    p["KSCAMP"] = 2.50        # /h
    p["KDCAMP"] = 2.50        # /h
    p["EMAXMOR"] = 0.82       # max fractional cAMP suppression by enteric MOR
    p["KACH"] = 1.50          # /h
    p["HACH"] = 1.60          # cAMP -> ACh release exponent
    p["KNO"] = 0.30
    p["ANO"] = 0.50           # MOR -> inhibitory NO/VIP tone
    p["KTONE"] = 1.20
    p["ATONE"] = 0.80         # MOR -> non-propulsive segmental tone
    p["NTONE"] = 0.85         # tone exponent damping propulsion
    p["ANOVIP"] = 0.50

    # secretagogues
    p["KCGMP"] = 0.50
    p["EC50LIN"] = 55.0       # ug in lumen, half-max GC-C activation
    p["EMAXCGMP"] = 2.60
    p["KCLC2"] = 0.60
    p["EC50LUB"] = 14.0       # ug in lumen, half-max ClC-2 activation
    p["EMAXCLC2"] = 1.00
    # methadone blocks ClC-2 -- the reason lubiprostone failed in methadone
    # users while naldemedine/naloxegol did not.  IC50 in opioid-occupancy units.
    p["ICLC2_METH"] = 0.35
    p["ECLC2"] = 2.30         # ClC-2 -> secretion gain
    # DEFECT #3 (found by failing to reproduce an endpoint, not by reading):
    # the first build had enteric MOR act on MOTILITY only.  It could then
    # not produce methylnaltrexone's defining result -- laxation within 4 h --
    # topping out at P=0.19 even with COMPLETE enteric blockade, because
    # stool hydration cannot change in 4 h if the only route to it is transit
    # time.  Enteric MOR also sits on submucosal secretomotor neurons and
    # suppresses active anion/fluid secretion.  Blocking it raises luminal
    # water within minutes, which is the fast arm of acute laxation.
    p["ESECMOR"] = 0.55       # max fractional suppression of secretion by MOR
    p["ECGMP"] = 1.05         # cGMP -> CFTR secretion gain

    # ===================== MOTILITY ========================================
    p["HAPC0"] = 0.28         # /h baseline high-amplitude propagating contractions
    p["KH"] = 1.50
    p["KMACH"] = 0.85
    p["HDIST"] = 0.45         # colonic filling gates HAPC

    # ===================== COLONIC CONTENT =================================
    p["SIN"] = 1.46           # g dry solids/h entering caecum (35 g/day)
    p["WIN"] = 62.5           # mL/h water entering caecum (1500 mL/day)
    # DEFECT #1 (found by integrating, not by reading): water absorption was
    # written first-order in luminal water, so the colon could dry a segment
    # to w->0, WFUN->0, propulsion->0, and solids accumulated WITHOUT BOUND --
    # in the HEALTHY arm (Stot 1064 g, SBM 0.02/wk).  Physiologically the
    # mucosa cannot extract water that is bound to the solid phase: absorption
    # must stop at w_min = WBIND/(1+WBIND).  With WBIND 1.60 that floor is
    # 0.615, which is also the measured water content of hard (Bristol 1-2)
    # stool.  This term is what makes the transit-hydration loop an AMPLIFIER
    # rather than a runaway.
    p["WBIND"] = 1.60         # mL water bound per g dry solids (unabsorbable)
    # PEG 3350 retains water by hydrogen bonding as well as colligatively.
    # Hydration number ~2.5 H2O per ethylene-oxide unit, 3350/44 = 76 units
    # => 76 x 2.5 x 18 / 3350 = 1.02 g water per g polymer.  Entered as an
    # osmole-equivalent so it uses the same unabsorbable-water machinery.
    p["WPEG"] = 1.02          # mL H-bonded water per g PEG 3350
    p["KABS"] = [2.050, 0.137, 0.211, 0.0279]  # /h water absorption per segment
    p["VSEC"] = [3.2, 2.2, 1.5, 0.9]           # mL/h baseline secretion per segment
    p["KPROP"] = [0.232, 0.278, 0.510, 0.000]  # /h base propulsion (4th unused)
    p["W50"] = 0.680          # stool water fraction at half-max propulsive efficiency
    p["GW"] = 3.5             # steepness of the hydration->propulsion term
    p["OSMPL"] = 0.290        # mmol osmoles per mL isotonic water
    p["KFERM"] = 0.075        # /h lactulose fermentation rate (colonic bacteria)
    p["FERMAMP"] = 3.60       # osmole amplification: 1 lactulose -> ~3.6 SCFA/gas osmoles
    p["FGAS"] = 0.30          # fraction of fermentation osmoles that are gas (bloating)

    # ===================== DEFECATION ======================================
    p["KDEF"] = 0.415         # /h maximal urge/defaecation rate
    p["VDEF"] = 170.0         # g rectosigmoid load at half-max urge
    # straining can expel hard stool at the anus even though the colon
    # cannot propel it, so the anorectal hydration gate is milder than the
    # colonic one (W50D < W50, GWD < GW).
    p["W50D"] = 0.550
    p["GWD"] = 2.50
    p["HDEF"] = 2.60
    p["HHAPC"] = 0.50
    p["AANO"] = 3.00          # MOR -> internal anal sphincter tone / blunted RAIR
    # DEFECT #4: with only a sphincter-tone brake the chronic arm could not
    # be pushed below ~2.3 SBM/wk without an implausible AANO, because a
    # heavily loaded rectum kept driving the urge term up.  The missing
    # mechanism is opioid-induced RECTAL HYPOSENSITIVITY: enteric/spinal MOR
    # raises the distension volume at which urge is perceived.  That is the
    # "reduced urge / incomplete evacuation" complaint patients actually
    # report, and it is what a PAMORA restores.  It enters as a shift of the
    # urge threshold itself, not as a gain on the rate.
    p["AVSENS"] = 2.50        # MOR-driven rise in the rectal urge threshold
    p["FEVAC"] = 0.85         # fraction of rectosigmoid content voided per event

    # --------- rescue laxative (KODIAC / COMPOSE protocol: bisacodyl or
    # enema permitted after 72 h with no bowel movement).  DEFECT #2: without
    # this the model has no way to shed the retained load, so the OIC arm
    # integrates to a colon holding >2 kg of solids.  It is also required for
    # ENDPOINT FIDELITY: an "SBM" is by definition a bowel movement with NO
    # rescue laxative in the preceding 24 h, so rescue-driven evacuations must
    # empty the colon WITHOUT being counted.
    p["KRESC"] = 0.012        # /h maximal rescue-laxative-driven BM rate (~2/wk cap)
    p["RTHR"] = 0.0139        # /h spontaneous rate below which rescue is taken (=1/72 h)
    p["FEVACR"] = 0.75        # rectosigmoid fraction voided per rescue BM
    p["FEVACR3"] = 0.35       # descending-colon fraction recruited by rescue
    p["RESCUE_ON"] = 1.0      # set 0 to disable rescue (mechanistic experiments)

    # ===================== SYMPTOMS ========================================
    p["KBSFS"] = 0.30
    p["WB50"] = 0.766
    p["HB"] = 9.0
    p["KSTR"] = 0.25
    p["W50S"] = 0.660         # water fraction at half-max straining relief
    p["HS"] = 8.0
    p["KDIST"] = 0.10
    p["KSYM"] = 0.12
    p["KQOL"] = 0.045
    p["DIST0"] = 900.0        # g total colonic content giving distension score 2

    # ===================== CNS / PAIN / WITHDRAWAL =========================
    p["PAIN0"] = 7.2          # untreated chronic pain NRS
    p["PAINMAX"] = 6.0        # max opioid analgesia (NRS units)
    p["KPAIN"] = 0.25
    p["KTOLA"] = 0.0016       # /h non-receptor counter-adaptation
    p["KTOLAR"] = 0.0022
    p["ATOLA"] = 1.4          # counter-adaptation weight on pain
    p["APAINGI"] = 0.45       # abdominal pain contribution of constipation to NRS
    p["WDMAX"] = 34.0         # COWS scale ceiling used here
    p["WD50"] = 0.145         # CNS antagonist occupancy at half-max withdrawal
    p["HWD"] = 2.6
    p["KWD"] = 0.9
    p["DEPEND"] = 1.0         # physical dependence multiplier
    p["KNAUS"] = 0.30
    p["ENAUSLUB"] = 0.9       # lubiprostone systemic nausea
    p["ENAUSDIST"] = 0.45

    # ===================== SAFETY ==========================================
    p["KIMP"] = 0.010         # impaction hazard accrual

    # ===================== DOSING (filled by scenarios) ====================
    p["op_dose"] = 0.0        # mg per administration
    p["op_int"] = 12.0        # h
    p["pam_dose"] = 0.0
    p["pam_int"] = 24.0
    p["pro_dose"] = 0.0
    p["pro_int"] = 24.0
    p["peg_dose"] = 0.0       # g/day PEG3350
    p["lac_dose"] = 0.0       # g/day lactulose
    p["lub_dose"] = 0.0       # ug per dose
    p["lub_int"] = 12.0
    p["lin_dose"] = 0.0       # ug per dose
    p["lin_int"] = 24.0

    # ===================== MODIFIERS =======================================
    p["PGPINH"] = 1.0         # multiplier on antagonist Kp_uu (P-gp inhibition)
    p["CYP3A4INH"] = 1.0      # divisor on antagonist clearance
    p["METHADONE"] = 0.0      # 1 = index opioid is methadone (ClC-2 block)
    p["FLUIDX"] = 1.0         # patient fluid intake multiplier (on WIN)
    p["ABSX"] = 1.0           # patient colonic water-absorption multiplier
    return p


# Per-antagonist parameter blocks.
#   KIANT   = MOR binding Ki (nM, unbound) -- literature
#   KPUU    = unbound brain / unbound plasma -- the selectivity-determining term
#   KIGUT   = OPERATIONAL enteric potency (nM).  This is the ONLY fitted
#             quantity per drug (see analysis G).  When KIGUT == KIANT the
#             plasma-unbound binding constant fully explains the trial.
ANTAGONISTS = {
    "none": dict(
        ANT_NAME="none", MWANT=400.0, FANT=0.0, KAANT=1.0, CLANT=20.0,
        V1ANT=80.0, V2ANT=80.0, QANT=10.0, FUANT=1.0, KPUU=0.0, KEOANT=0.5,
        KIANT=1e9, KIGUT=1e9, FLUMANT=0.0, ORALANT=1, FLUMACCX=1.0,
    ),
    # ---- naloxegol (PEGylated naloxol; P-gp substrate by design) ----------
    "naloxegol": dict(
        ANT_NAME="naloxegol", MWANT=651.8, FANT=0.50, KAANT=1.0, CLANT=22.0,
        V1ANT=68.0, V2ANT=100.0, QANT=15.0, FUANT=0.96, KPUU=0.020,
        KEOANT=0.50, KIANT=7.4, KIGUT=2.6851, FLUMANT=0.35, ORALANT=1,
        FLUMACCX=1.0,
    ),
    # ---- naldemedine -----------------------------------------------------
    "naldemedine": dict(
        ANT_NAME="naldemedine", MWANT=570.6, FANT=0.55, KAANT=2.5, CLANT=8.5,
        V1ANT=155.0, V2ANT=90.0, QANT=12.0, FUANT=0.065, KPUU=0.012,
        KEOANT=0.50, KIANT=0.34, KIGUT=0.00821, FLUMANT=0.30, ORALANT=1,
        FLUMACCX=1.0,
    ),
    # ---- methylnaltrexone, SC (quaternary amine: permanently charged) -----
    "methylnaltrexone_sc": dict(
        ANT_NAME="methylnaltrexone SC", MWANT=356.5, FANT=0.82, KAANT=3.0,
        CLANT=17.6, V1ANT=84.0, V2ANT=60.0, QANT=20.0, FUANT=0.885,
        KPUU=0.005, KEOANT=0.50, KIANT=28.0, KIGUT=0.7945, FLUMANT=0.02,
        ORALANT=0, FLUMACCX=0.0015,
    ),
    # ---- methylnaltrexone, oral 450 mg (works from the LUMEN) -------------
    "methylnaltrexone_po": dict(
        ANT_NAME="methylnaltrexone PO", MWANT=356.5, FANT=0.012, KAANT=0.9,
        CLANT=17.6, V1ANT=84.0, V2ANT=60.0, QANT=20.0, FUANT=0.885,
        KPUU=0.005, KEOANT=0.50, KIANT=28.0, KIGUT=0.7945, FLUMANT=0.85,
        ORALANT=1, FLUMACCX=0.0015,
    ),
    # ---- naloxone PO (the negative control: Kp_uu ~ 1) -------------------
    "naloxone_po": dict(
        ANT_NAME="naloxone PO", MWANT=327.4, FANT=0.02, KAANT=1.5, CLANT=90.0,
        V1ANT=180.0, V2ANT=120.0, QANT=30.0, FUANT=0.55, KPUU=1.00,
        KEOANT=0.80, KIANT=1.1, KIGUT=1.1, FLUMANT=0.75, ORALANT=1,
        FLUMACCX=1.0,
    ),
}


def set_antagonist(p, name):
    p = dict(p)
    p.update(ANTAGONISTS[name])
    return p


# ---------------------------------------------------------------------------
# 2.  HELPERS
# ---------------------------------------------------------------------------
def hill(x, x50, n):
    x = max(x, 0.0)
    if x <= 0.0:
        return 0.0
    xn = x ** n
    return xn / (xn + x50 ** n)


def wfrac(W, S):
    tot = W + S
    if tot <= 1e-9:
        return 0.75
    return min(max(W / tot, 0.0), 1.0)


# ---------------------------------------------------------------------------
# 3.  RIGHT-HAND SIDE
# ---------------------------------------------------------------------------
def rhs(t, y, p):
    y = np.maximum(y, 0.0)
    d = np.zeros(NST)
    g = IX

    # ---------------- opioid PK ----------------
    Cop = y[g["AOP_CEN"]] / p["V1OP"]                    # mg/L
    Cop_nM = Cop / p["MWOP"] * 1e6
    Aop_u = p["FUOP"] * Cop_nM                            # unbound plasma nM
    abs_op = p["KAOP"] * y[g["AOP_DEP"]]
    d[g["AOP_DEP"]] = -abs_op
    d[g["AOP_CEN"]] = (p["FOP"] * abs_op
                       - (p["CLOP"] / p["V1OP"]) * y[g["AOP_CEN"]]
                       - (p["QOP"] / p["V1OP"]) * y[g["AOP_CEN"]]
                       + (p["QOP"] / p["V2OP"]) * y[g["AOP_PER"]])
    d[g["AOP_PER"]] = ((p["QOP"] / p["V1OP"]) * y[g["AOP_CEN"]]
                       - (p["QOP"] / p["V2OP"]) * y[g["AOP_PER"]])
    d[g["COP_BR"]] = p["KEOOP"] * (p["KPUUOP"] * Aop_u - y[g["COP_BR"]])
    d[g["AOP_LUM"]] = p["FLUMOP"] * abs_op - p["KLUMOUT"] * y[g["AOP_LUM"]]

    # pre-systemic gut-wall increment (absorption flux / portal flow)
    gw_op = p["FUOP"] * (p["FOP"] * abs_op / p["QPORT"]) / p["MWOP"] * 1e6
    lum_op = (y[g["AOP_LUM"]] / p["VLUM"]) / p["MWOP"] * 1e6 * p["FLUMACC"]
    A_gut = Aop_u + gw_op + lum_op
    A_cns = y[g["COP_BR"]]

    # ---------------- antagonist PK ----------------
    CL_ant = p["CLANT"] / p["CYP3A4INH"]
    Cpam = y[g["APAM_CEN"]] / p["V1ANT"]
    Cpam_nM = Cpam / p["MWANT"] * 1e6
    Bp_u = p["FUANT"] * Cpam_nM
    abs_pam = p["KAANT"] * y[g["APAM_DEP"]]
    d[g["APAM_DEP"]] = -abs_pam
    d[g["APAM_CEN"]] = (p["FANT"] * abs_pam
                        - (CL_ant / p["V1ANT"]) * y[g["APAM_CEN"]]
                        - (p["QANT"] / p["V1ANT"]) * y[g["APAM_CEN"]]
                        + (p["QANT"] / p["V2ANT"]) * y[g["APAM_PER"]])
    d[g["APAM_PER"]] = ((p["QANT"] / p["V1ANT"]) * y[g["APAM_CEN"]]
                        - (p["QANT"] / p["V2ANT"]) * y[g["APAM_PER"]])
    # P-gp inhibition acts on Kp_uu, NOT on clearance.  That separation is the
    # whole of analysis B.
    kpuu_eff = min(p["KPUU"] * p["PGPINH"], 1.0)
    d[g["CPAM_BR"]] = p["KEOANT"] * (kpuu_eff * Bp_u - y[g["CPAM_BR"]])
    d[g["APAM_LUM"]] = (p["FLUMANT"] * abs_pam * p["ORALANT"]
                        - p["KLUMOUT"] * y[g["APAM_LUM"]])

    gw_pam = p["FUANT"] * (p["FANT"] * abs_pam * p["ORALANT"] / p["QPORT"]) \
        / p["MWANT"] * 1e6
    lum_pam = (y[g["APAM_LUM"]] / p["VLUM"]) / p["MWANT"] * 1e6 * p["FLUMACC"] \
        * p["FLUMACCX"]
    B_gut_free = Bp_u + gw_pam + lum_pam
    B_cns = y[g["CPAM_BR"]]

    # ---------------- competitive occupancy ----------------
    # enteric plexus: operational potency KIGUT (fitted, one per drug)
    rg_a = A_gut / p["KIOP"]
    rg_b = B_gut_free / p["KIGUT"]
    den_g = 1.0 + rg_a + rg_b
    OCCg_ag = rg_a / den_g
    OCCg_ant = rg_b / den_g
    # CNS: binding Ki (unfitted) -- the safety side must not borrow a fit
    rc_a = A_cns / p["KIOP"]
    rc_b = B_cns / p["KIANT"]
    den_c = 1.0 + rc_a + rc_b
    OCCc_ag = rc_a / den_c
    OCCc_ant = rc_b / den_c

    SIGg = hill(OCCg_ag * y[g["RG_AV"]], p["OCC50G"], p["HG"])
    SIGc = hill(OCCc_ag * y[g["RC_AV"]], p["OCC50C"], p["HC"])

    # ---------------- receptor trafficking ----------------
    d[g["RG_AV"]] = (p["KINRG"] * (1.0 - y[g["RG_AV"]])
                     - p["KDESG"] * OCCg_ag * y[g["RG_AV"]])
    d[g["RC_AV"]] = (p["KINRC"] * (1.0 - y[g["RC_AV"]])
                     - p["KDESC"] * OCCc_ag * (1.0 + p["BARR"] * y[g["ARR"]])
                     * y[g["RC_AV"]])
    d[g["ARR"]] = p["KONARR"] * OCCc_ag - p["KOFFARR"] * y[g["ARR"]]

    # ---------------- prucalopride (5-HT4, Gs on the SAME adenylyl cyclase) --
    abs_pro = p["KAPRO"] * 0.0  # depot handled as direct central input below
    Cpro = y[g["APRO_CEN"]] / p["V1PRO"] / p["MWPRO"] * 1e6   # nM
    E5HT4 = p["EMAXPRO"] * Cpro / (Cpro + p["EC50PRO"])
    d[g["APRO_CEN"]] = (-(p["CLPRO"] / p["V1PRO"]) * y[g["APRO_CEN"]]
                        - (p["QPRO"] / p["V1PRO"]) * y[g["APRO_CEN"]]
                        + (p["QPRO"] / p["V2PRO"]) * y[g["APRO_PER"]] + abs_pro)
    d[g["APRO_PER"]] = ((p["QPRO"] / p["V1PRO"]) * y[g["APRO_CEN"]]
                        - (p["QPRO"] / p["V2PRO"]) * y[g["APRO_PER"]])

    # ---------------- enteric second messengers ----------------
    d[g["CAMP"]] = (p["KSCAMP"] * (1.0 - p["EMAXMOR"] * SIGg + E5HT4)
                    - p["KDCAMP"] * y[g["CAMP"]])
    d[g["ACH"]] = p["KACH"] * (y[g["CAMP"]] ** p["HACH"] - y[g["ACH"]])
    d[g["NOVIP"]] = p["KNO"] * (1.0 + p["ANO"] * SIGg - y[g["NOVIP"]])

    # secretagogues (luminal, non-absorbed)
    d[g["ALUB_LUM"]] = -p["KLUMOUT"] * y[g["ALUB_LUM"]]
    d[g["ALIN_LUM"]] = -p["KLUMOUT"] * y[g["ALIN_LUM"]]
    # methadone blocks ClC-2 -- a DIRECT channel block, not a MOR effect
    clc2_block = 1.0
    if p["METHADONE"] > 0.5:
        clc2_block = 1.0 / (1.0 + OCCg_ag / p["ICLC2_METH"])
    clc2_drive = (p["EMAXCLC2"] * hill(y[g["ALUB_LUM"]], p["EC50LUB"], 1.0)
                  * clc2_block)
    d[g["CLC2"]] = p["KCLC2"] * (clc2_drive - y[g["CLC2"]])
    cgmp_drive = 1.0 + p["EMAXCGMP"] * hill(y[g["ALIN_LUM"]], p["EC50LIN"], 1.0)
    d[g["CGMP"]] = p["KCGMP"] * (cgmp_drive - y[g["CGMP"]])

    # ---------------- motility ----------------
    Stot = sum(y[g[k]] for k in ("S1", "S2", "S3", "S4"))
    Wtot = sum(y[g[k]] for k in ("W1", "W2", "W3", "W4"))
    fill_gate = (max(Stot + Wtot, 1.0) / 300.0) ** p["HDIST"]
    ach_term = (y[g["ACH"]] / (y[g["ACH"]] + p["KMACH"])) / (1.0 / (1.0 + p["KMACH"]))
    hapc_drive = (p["HAPC0"] * ach_term * fill_gate
                  / (1.0 + p["ANOVIP"] * max(y[g["NOVIP"]] - 1.0, 0.0)))
    d[g["HAPC"]] = p["KH"] * (hapc_drive - y[g["HAPC"]])
    d[g["TONE"]] = p["KTONE"] * (1.0 + p["ATONE"] * SIGg - y[g["TONE"]])

    hapc_rel = max(y[g["HAPC"]], 1e-6) / p["HAPC0"]
    tone_div = max(y[g["TONE"]], 1e-3) ** p["NTONE"]

    # ---------------- colonic content ----------------
    S = [y[g["S1"]], y[g["S2"]], y[g["S3"]], y[g["S4"]]]
    W = [y[g["W1"]], y[g["W2"]], y[g["W3"]], y[g["W4"]]]
    O = [y[g["O1"]], y[g["O2"]], y[g["O3"]], y[g["O4"]]]

    kprop = []
    wf = []
    for i in range(4):
        w = wfrac(W[i], S[i])
        wf.append(w)
        kprop.append(p["KPROP"][i] * hapc_rel * hill(w, p["W50"], p["GW"]) / tone_div)

    sec_gain = ((1.0 - p["ESECMOR"] * SIGg)
                + p["ECLC2"] * y[g["CLC2"]]
                + p["ECGMP"] * max(y[g["CGMP"]] - 1.0, 0.0))

    # osmotic laxatives enter segment 1
    d[g["APEG_LUM"]] = -kprop[0] * y[g["APEG_LUM"]]
    ferm = p["KFERM"] * y[g["ALAC_LUM"]]
    d[g["ALAC_LUM"]] = -ferm - kprop[0] * y[g["ALAC_LUM"]]
    # osmoles delivered to segment 1: unfermented lactulose + PEG + fermentation
    osm_in = ferm * p["FERMAMP"]

    # PEG and lactulose themselves are osmotically active while present
    O_extra = [y[g["APEG_LUM"]] * (1.0 + 3.350 * p["WPEG"] * p["OSMPL"])
               + y[g["ALAC_LUM"]], 0.0, 0.0, 0.0]

    # ---------------- defecation ----------------
    load4 = S[3] + W[3]
    vdef_eff = p["VDEF"] * (1.0 + p["AVSENS"] * SIGg)
    Rdef = (p["KDEF"] * hill(load4, vdef_eff, p["HDEF"]) * hill(wf[3], p["W50D"], p["GWD"])
            * (hapc_rel ** p["HHAPC"]) / (1.0 + p["AANO"] * SIGg))
    # rescue laxative: taken when the SPONTANEOUS rate falls below 1 per 72 h
    # DEFECT #6: a piecewise-linear rescue trigger put a genuine DISCONTINUITY
    # in the steady state (SBM jumped 2.44 -> 1.57 between AANO 3.0 and 3.2)
    # because rescue drains the rectum, which lowers the urge term, which
    # lowers the spontaneous rate further.  Smooth the trigger so the residual
    # steepness is the BIOLOGY (the transit-hydration loop) and not the rule.
    resc_need = 1.0 - hill(Rdef, p["RTHR"], 2.5)
    Rresc = p["KRESC"] * p["RESCUE_ON"] * resc_need
    evac4 = p["FEVAC"] * Rdef + p["FEVACR"] * Rresc
    evac3 = p["FEVACR3"] * Rresc

    # outflow rate constants first (segment 3 is also recruited by rescue),
    # THEN inflows -- otherwise the mass that rescue pulls out of the
    # descending colon never arrives in the rectosigmoid and is destroyed.
    kout = [kprop[0], kprop[1], kprop[2] + evac3, evac4]
    inS = [p["SIN"], kout[0] * S[0], kout[1] * S[1], kout[2] * S[2]]
    inW = [p["WIN"] * p["FLUIDX"], kout[0] * W[0], kout[1] * W[1], kout[2] * W[2]]
    inO = [osm_in, kout[0] * O[0], kout[1] * O[1], kout[2] * O[2]]

    for i in range(4):
        Wosm = (O[i] + O_extra[i]) / p["OSMPL"]
        # absorbable water = free water above BOTH the osmotically obliged
        # volume and the water bound to the solid phase (DEFECT #1)
        wfree = W[i] - max(Wosm, p["WBIND"] * S[i])
        absw = p["KABS"][i] * p["ABSX"] * max(wfree, 0.0)
        secw = p["VSEC"][i] * sec_gain
        out_i = kout[i]
        d[g[f"S{i+1}"]] = inS[i] - out_i * S[i]
        d[g[f"W{i+1}"]] = inW[i] + secw - absw - out_i * W[i]
        d[g[f"O{i+1}"]] = inO[i] - out_i * O[i]

    d[g["CUM_SBM"]] = Rdef
    d[g["CUM_RESC"]] = Rresc
    pcomp = hill(wf[3], p["W50"], p["GW"]) * min(hapc_rel, 1.5)
    d[g["CUM_CSBM"]] = Rdef * min(pcomp, 1.0)

    # ---------------- symptoms ----------------
    bsfs_t = 1.0 + 6.0 * hill(wf[3], p["WB50"], p["HB"])
    d[g["BSFS"]] = p["KBSFS"] * (bsfs_t - y[g["BSFS"]])
    strain_t = 4.0 * (1.0 - hill(wf[3], p["W50S"], p["HS"]))
    d[g["STRAIN"]] = p["KSTR"] * (min(strain_t, 4.0) - y[g["STRAIN"]])
    gas = p["FGAS"] * sum(O)
    dist_t = 4.0 * hill(Stot + Wtot + 3.0 * gas, p["DIST0"], 2.0)
    d[g["DIST"]] = p["KDIST"] * (dist_t - y[g["DIST"]])
    sbm_wk = 168.0 * Rdef
    freq_sym = 4.0 * (1.0 - hill(sbm_wk, 3.5, 2.0))
    sym_t = 0.34 * freq_sym + 0.28 * y[g["STRAIN"]] + 0.22 * y[g["DIST"]] \
        + 0.16 * (4.0 * (1.0 - hill(y[g["BSFS"]], 3.0, 4.0)))
    d[g["PACSYM"]] = p["KSYM"] * (sym_t - y[g["PACSYM"]])
    d[g["QOL"]] = p["KQOL"] * (0.85 * y[g["PACSYM"]] - y[g["QOL"]])

    # ---------------- CNS ----------------
    tol_pen = p["ATOLA"] * y[g["TOLA"]]
    pain_t = (p["PAIN0"] - p["PAINMAX"] * SIGc + tol_pen
              + p["APAINGI"] * y[g["DIST"]] * 0.5)
    d[g["PAIN"]] = p["KPAIN"] * (min(max(pain_t, 0.0), 10.0) - y[g["PAIN"]])
    d[g["TOLA"]] = p["KTOLA"] * OCCc_ag - p["KTOLAR"] * y[g["TOLA"]]
    wd_t = p["WDMAX"] * p["DEPEND"] * hill(OCCc_ant, p["WD50"], p["HWD"])
    d[g["WD"]] = p["KWD"] * (wd_t - y[g["WD"]])
    naus_t = (p["ENAUSLUB"] * hill(y[g["ALUB_LUM"]], p["EC50LUB"] * 2.0, 1.0)
              + p["ENAUSDIST"] * y[g["DIST"]] * 0.5)
    d[g["NAUSEA"]] = p["KNAUS"] * (min(naus_t, 3.0) - y[g["NAUSEA"]])

    # ---------------- safety ----------------
    d[g["IMPACT"]] = p["KIMP"] * max(1.0 - hill(sbm_wk, 1.6, 3.0), 0.0) \
        * max(1.0 - hill(wf[3], 0.62, 8.0), 0.0)

    return d


# ---------------------------------------------------------------------------
# 4.  SIMULATION DRIVER
# ---------------------------------------------------------------------------
def init_state(p):
    y = np.zeros(NST)
    g = IX
    y[g["RG_AV"]] = 1.0
    y[g["RC_AV"]] = 1.0
    y[g["CAMP"]] = 1.0
    y[g["ACH"]] = 1.0
    y[g["NOVIP"]] = 1.0
    y[g["CGMP"]] = 1.0
    y[g["HAPC"]] = p["HAPC0"]
    y[g["TONE"]] = 1.0
    y[g["S1"]], y[g["S2"]], y[g["S3"]], y[g["S4"]] = 11.7, 10.2, 5.8, 27.5
    y[g["W1"]], y[g["W2"]], y[g["W3"]], y[g["W4"]] = 46.8, 36.2, 18.5, 82.5
    y[g["BSFS"]] = 4.0
    y[g["PAIN"]] = p["PAIN0"]
    return y


def dose_events(p, tmax):
    """Return sorted list of (time, state_index, amount)."""
    ev = []
    g = IX
    if p["op_dose"] > 0 and p["op_int"] > 0:
        t = 0.0
        while t < tmax:
            ev.append((t, g["AOP_DEP"], p["op_dose"]))
            t += p["op_int"]
    if p["pam_dose"] > 0 and p["pam_int"] > 0:
        t = p.get("pam_start", 0.0)
        while t < tmax:
            ev.append((t, g["APAM_DEP"], p["pam_dose"]))
            t += p["pam_int"]
    if p["pro_dose"] > 0:
        t = p.get("pro_start", 0.0)
        while t < tmax:
            ev.append((t, g["APRO_CEN"], p["pro_dose"] * p["FPRO"]))
            t += p["pro_int"]
    if p["peg_dose"] > 0:
        # PEG3350 MW 3350 -> mmol osmoles
        t = 0.0
        while t < tmax:
            ev.append((t, g["APEG_LUM"], p["peg_dose"] / 3.350))
            t += 24.0
    if p["lac_dose"] > 0:
        # lactulose MW 342.3 -> mmol
        t = 0.0
        while t < tmax:
            ev.append((t, g["ALAC_LUM"], p["lac_dose"] / 0.3423))
            t += 24.0
    if p["lub_dose"] > 0:
        t = 0.0
        while t < tmax:
            ev.append((t, g["ALUB_LUM"], p["lub_dose"]))
            t += p["lub_int"]
    if p["lin_dose"] > 0:
        t = 0.0
        while t < tmax:
            ev.append((t, g["ALIN_LUM"], p["lin_dose"]))
            t += p["lin_int"]
    ev.sort(key=lambda e: e[0])
    return ev


def simulate(p, tmax, y0=None, nout_per_h=1.0):
    """Integrate with bolus dose events.  Returns (t, Y) with Y (nt, NST)."""
    y = init_state(p) if y0 is None else np.array(y0, dtype=float)
    ev = dose_events(p, tmax)
    # break points
    times = sorted(set([0.0] + [e[0] for e in ev] + [tmax]))
    T, Y = [], []
    for k in range(len(times) - 1):
        t0, t1 = times[k], times[k + 1]
        for (te, idx, amt) in ev:
            if abs(te - t0) < 1e-9:
                y[idx] += amt
        if t1 <= t0:
            continue
        n = max(2, int(math.ceil((t1 - t0) * nout_per_h)) + 1)
        teval = np.linspace(t0, t1, n)
        sol = solve_ivp(rhs, (t0, t1), y, args=(p,), method="LSODA",
                        t_eval=teval, rtol=1e-6, atol=1e-9, max_step=1.0)
        if not sol.success:
            raise RuntimeError(f"integration failed at t={t0}: {sol.message}")
        T.append(sol.t if k == 0 else sol.t[1:])
        Y.append(sol.y.T if k == 0 else sol.y.T[1:])
        y = sol.y[:, -1].copy()
    for (te, idx, amt) in ev:
        if abs(te - tmax) < 1e-9:
            y[idx] += amt
    return np.concatenate(T), np.vstack(Y)


# ---------------------------------------------------------------------------
# 5.  READOUTS
# ---------------------------------------------------------------------------
def readouts(t, Y, p, window_h=168.0):
    """Summarise the LAST `window_h` hours."""
    g = IX
    m = t >= (t[-1] - window_h)
    tt, YY = t[m], Y[m]
    dt = tt[-1] - tt[0]
    sbm = (YY[-1, g["CUM_SBM"]] - YY[0, g["CUM_SBM"]]) / dt * 168.0
    csbm = (YY[-1, g["CUM_CSBM"]] - YY[0, g["CUM_CSBM"]]) / dt * 168.0
    resc = (YY[-1, g["CUM_RESC"]] - YY[0, g["CUM_RESC"]]) / dt * 168.0

    # recompute mechanistic algebra at the final point
    yf = Y[-1]
    Cop_nM = (yf[g["AOP_CEN"]] / p["V1OP"]) / p["MWOP"] * 1e6
    Aop_u = p["FUOP"] * Cop_nM
    Cpam_nM = (yf[g["APAM_CEN"]] / p["V1ANT"]) / p["MWANT"] * 1e6
    Bp_u = p["FUANT"] * Cpam_nM
    lum_op = (yf[g["AOP_LUM"]] / p["VLUM"]) / p["MWOP"] * 1e6 * p["FLUMACC"]
    lum_pam = (yf[g["APAM_LUM"]] / p["VLUM"]) / p["MWANT"] * 1e6 * p["FLUMACC"] \
        * p["FLUMACCX"]
    A_gut, B_gut = Aop_u + lum_op, Bp_u + lum_pam
    rg_a, rg_b = A_gut / p["KIOP"], B_gut / p["KIGUT"]
    dg = 1 + rg_a + rg_b
    rc_a, rc_b = yf[g["COP_BR"]] / p["KIOP"], yf[g["CPAM_BR"]] / p["KIANT"]
    dc = 1 + rc_a + rc_b

    # transit time: total colonic content / output flux
    Stot = sum(yf[g[k]] for k in ("S1", "S2", "S3", "S4"))
    Wtot = sum(yf[g[k]] for k in ("W1", "W2", "W3", "W4"))
    ctt = (Stot) / max(p["SIN"], 1e-9)

    return dict(
        SBM_wk=sbm, CSBM_wk=csbm, RESC_wk=resc,
        BSFS=yf[g["BSFS"]], STRAIN=yf[g["STRAIN"]], DIST=yf[g["DIST"]],
        PACSYM=yf[g["PACSYM"]], QOL=yf[g["QOL"]],
        PAIN=yf[g["PAIN"]], WD=yf[g["WD"]], NAUSEA=yf[g["NAUSEA"]],
        OCCg_ag=rg_a / dg, OCCg_ant=rg_b / dg,
        OCCc_ag=rc_a / dc, OCCc_ant=rc_b / dc,
        RG_AV=yf[g["RG_AV"]], RC_AV=yf[g["RC_AV"]],
        HAPC=yf[g["HAPC"]], TONE=yf[g["TONE"]], CAMP=yf[g["CAMP"]],
        w4=wfrac(yf[g["W4"]], yf[g["S4"]]),
        w1=wfrac(yf[g["W1"]], yf[g["S1"]]),
        Stot=Stot, Wtot=Wtot, CTT_h=ctt,
        IMPACT=yf[g["IMPACT"]],
        Cop_u_nM=Aop_u, Cpam_u_nM=Bp_u, Cpam_br_nM=yf[g["CPAM_BR"]],
        Cop_br_nM=yf[g["COP_BR"]],
    )


# ---------------------------------------------------------------------------
# 6.  SCENARIOS
# ---------------------------------------------------------------------------
OPIOID_STD = dict(op_dose=30.0, op_int=12.0)     # oxycodone 60 mg/day

SCENARIOS = [
    ("S01 healthy (no opioid)",            "none",                dict(op_dose=0.0)),
    ("S02 OIC untreated (oxycodone 60/d)", "none",                dict()),
    ("S03 + PEG3350 17 g/d",               "none",                dict(peg_dose=17.0)),
    ("S04 + lactulose 20 g/d",             "none",                dict(lac_dose=20.0)),
    ("S05 + lubiprostone 24 ug bid",       "none",                dict(lub_dose=24.0, lub_int=12.0)),
    ("S06 + linaclotide 145 ug qd",        "none",                dict(lin_dose=145.0)),
    ("S07 + prucalopride 2 mg qd",         "none",                dict(pro_dose=2.0)),
    ("S08 + naloxegol 12.5 mg qd",         "naloxegol",           dict(pam_dose=12.5)),
    ("S09 + naloxegol 25 mg qd",           "naloxegol",           dict(pam_dose=25.0)),
    ("S10 + naldemedine 0.2 mg qd",        "naldemedine",         dict(pam_dose=0.2)),
    ("S11 + methylnaltrexone 12 mg SC qod","methylnaltrexone_sc", dict(pam_dose=12.0, pam_int=48.0)),
    ("S12 + methylnaltrexone 450 mg PO qd","methylnaltrexone_po", dict(pam_dose=450.0)),
    ("S13 + naloxone 20 mg PO tid",        "naloxone_po",         dict(pam_dose=20.0, pam_int=8.0)),
    ("S14 + naloxegol 25 + PEG 17 g",      "naloxegol",           dict(pam_dose=25.0, peg_dose=17.0)),
    ("S15 naloxegol 25 + strong P-gp inh", "naloxegol",           dict(pam_dose=25.0, PGPINH=10.0, CYP3A4INH=3.4)),
    ("S16 methadone 60/d + lubiprostone",  "none",                dict(lub_dose=24.0, lub_int=12.0, METHADONE=1.0)),
]


def run_scenario(name, drug, mods, days=84):
    p = base_params()
    p = set_antagonist(p, drug)
    p.update(OPIOID_STD)
    p.update(mods)
    t, Y = simulate(p, 24 * days)
    r = readouts(t, Y, p)
    r["name"] = name
    r["drug"] = p["ANT_NAME"]
    return r, p, t, Y


# ---------------------------------------------------------------------------
# 7.  ANALYSES
# ---------------------------------------------------------------------------
def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def analysis_A():
    hr("A.  SELECTIVITY INDEX -- computed, not asserted\n"
       "    SI = enteric antagonist occupancy / CNS antagonist occupancy,\n"
       "    both at the SAME plasma concentration in the SAME patient.")
    rows = []
    specs = [("naloxone_po", 20.0, 8.0), ("naloxegol", 25.0, 24.0),
             ("naldemedine", 0.2, 24.0), ("methylnaltrexone_sc", 12.0, 24.0)]
    print(f"{'drug':22} {'Kp_uu':>7} {'1/Kp_uu':>8} {'OCC_gut':>8} {'OCC_cns':>9} "
          f"{'SI':>8} {'COWS':>6} {'dPain':>7} {'SBM/wk':>7}")
    base, _, _, _ = run_scenario("ref", "none", {})
    for drug, dose, ivl in specs:
        r, p, _, _ = run_scenario(drug, drug, dict(pam_dose=dose, pam_int=ivl))
        si = r["OCCg_ant"] / max(r["OCCc_ant"], 1e-12)
        rows.append((p["ANT_NAME"], p["KPUU"], si, r))
        print(f"{p['ANT_NAME']:22} {p['KPUU']:7.3f} {1/p['KPUU']:8.1f} "
              f"{r['OCCg_ant']:8.3f} {r['OCCc_ant']:9.5f} {si:8.1f} "
              f"{r['WD']:6.2f} {r['PAIN']-base['PAIN']:+7.2f} {r['SBM_wk']:7.2f}")
    print("\n  The SI column is not a fitted quantity and not an assumption: it is")
    print("  one division.  It spans %.0f-fold across these four drugs while the"
          % (max(x[2] for x in rows) / min(x[2] for x in rows)))
    print("  MOR binding Ki spans only 25-fold, and in the OPPOSITE order.")
    return rows


def analysis_B():
    hr("B.  A P-gp INHIBITOR DOES NOT RAISE EXPOSURE -- IT ROTATES THE RATIO.\n"
       "    Naloxegol 25 mg qd, perturbed three ways, with and without the\n"
       "    label-directed dose reduction to 12.5 mg.")
    cases = [
        ("no interaction",                dict(), 25.0),
        ("CYP3A4 inhibitor only (AUCx3.4)", dict(CYP3A4INH=3.4), 25.0),
        ("CYP3A4 inhibitor, dose 12.5",   dict(CYP3A4INH=3.4), 12.5),
        ("P-gp inhibitor only (Kp_uu x10)", dict(PGPINH=10.0), 25.0),
        ("P-gp inhibitor, dose 12.5",     dict(PGPINH=10.0), 12.5),
        ("dual CYP3A4+P-gp (verapamil-like)", dict(CYP3A4INH=3.4, PGPINH=10.0), 25.0),
        ("dual, dose reduced to 12.5",    dict(CYP3A4INH=3.4, PGPINH=10.0), 12.5),
    ]
    print(f"{'case':36} {'Cu,pl nM':>9} {'OCC_gut':>8} {'OCC_cns':>9} {'SI':>7} "
          f"{'COWS':>6} {'SBM/wk':>7}")
    for lab, mods, dose in cases:
        m = dict(mods); m["pam_dose"] = dose
        r, p, _, _ = run_scenario(lab, "naloxegol", m)
        si = r["OCCg_ant"] / max(r["OCCc_ant"], 1e-12)
        print(f"{lab:36} {r['Cpam_u_nM']:9.1f} {r['OCCg_ant']:8.3f} "
              f"{r['OCCc_ant']:9.5f} {si:7.1f} {r['WD']:6.2f} {r['SBM_wk']:7.2f}")
    print("""
  SI is dose-invariant only in the LOW-occupancy limit.  At label doses the
  enteric receptor is already 57% occupied, so halving the dose does move SI --
  but read the two perturbations against each other and they behave differently
  in kind, not just in degree:

    CYP3A4 inhibition is a pure EXPOSURE shift.  Kp_uu never moved.  Halving the
    dose walks both occupancies back down the same curve and recovers most of
    the lost margin (SI 17.5 -> 31.3, COWS 1.79 -> 0.33) at a cost in efficacy
    the patient can afford (SBM 7.10 -> 5.26, still above the 4.30 baseline).

    P-gp inhibition moves Kp_uu itself.  The plasma concentration is UNCHANGED
    (6.9 nM in both rows), enteric occupancy is UNCHANGED (0.569), and the only
    thing that moved is the brain.  Halving the dose now buys almost no safety
    (SI 9.8 -> 13.3, still 7-fold below baseline) and pays for it with efficacy
    the patient needed (SBM 4.30 -> 3.11, i.e. back below the response
    threshold).  The dose reduction is treating the wrong variable.

  The dual case is the clinically common one -- verapamil, diltiazem and
  quinidine all do both -- and it is the one where the label's advice is least
  adequate: COWS 30.4 falling only to 23.5 on the reduced dose.""")


def analysis_C():
    hr("C.  TOLERANCE IS ASYMMETRIC -- but the consequence is NOT what a\n"
       "    naive reading predicts.  Equilibrate 12 weeks on oxycodone 60 mg/day,\n"
       "    then titrate weekly for 24 weeks trying to hold pain NRS at 4.0.")
    p = base_params()
    p = set_antagonist(p, "none")
    p.update(OPIOID_STD)
    t, Y = simulate(p, 24 * 84)          # equilibrate BEFORE the titration
    y = Y[-1].copy()
    dose = 30.0
    target = 4.0
    print(f"{'week':>5} {'oxy mg/d':>9} {'OCC_cns':>8} {'RC_avail':>9} "
          f"{'pain':>6} {'OCC_gut':>8} {'RG_avail':>9} {'SBM/wk':>7} {'BSFS':>6}")
    trace = []
    for wk in range(0, 25):
        p["op_dose"] = dose
        t, Y = simulate(p, 168.0, y0=y)
        r = readouts(t, Y, p, window_h=168.0)
        y = Y[-1].copy()
        trace.append((wk, dose * 2, r))
        if wk % 4 == 0 or wk == 24:
            print(f"{wk:5d} {dose*2:9.1f} {r['OCCc_ag']:8.3f} {r['RC_AV']:9.3f} "
                  f"{r['PAIN']:6.2f} {r['OCCg_ag']:8.3f} {r['RG_AV']:9.3f} "
                  f"{r['SBM_wk']:7.2f} {r['BSFS']:6.2f}")
        if r["PAIN"] > target + 0.25:
            dose *= 1.22
        elif r["PAIN"] < target - 0.25:
            dose *= 0.94
        dose = min(dose, 100.0)          # cap at 200 mg/day, a real-world ceiling
    w0, wN = trace[0][2], trace[-1][2]
    print(f"""
  Read the columns, not the intuition.  Over 24 weeks the dose rises
  {trace[-1][1]/trace[0][1]:.1f}-fold.  CNS receptor availability falls
  {w0['RC_AV']:.3f} -> {wN['RC_AV']:.3f} while enteric availability holds at
  {wN['RG_AV']:.3f}: the asymmetry is real and it is {(1-wN['RC_AV'])/(1-wN['RG_AV']):.0f}-fold.

  But the gut endpoint does NOT deteriorate.  Enteric occupancy climbs
  {w0['OCCg_ag']:.3f} -> {wN['OCCg_ag']:.3f} and SBM/week goes
  {w0['SBM_wk']:.2f} -> {wN['SBM_wk']:.2f}.  This model therefore does NOT
  support the story that dose escalation progressively worsens OIC.  It says
  something more specific and more testable: enteric transduction is already
  saturated at an ordinary analgesic dose, so the constipation a patient has at
  60 mg/day is essentially the constipation they will have at 200 mg/day.

  That is the model's explanation for an epidemiological fact a simple
  dose-response account cannot handle: OIC prevalence correlates poorly with
  opioid dose.

  The clinically important asymmetry is on the OTHER side of the ledger.  The
  controller never reaches its target: pain sits at {wN['PAIN']:.2f} against a
  goal of {target:.1f} through the entire escalation, and stops improving once
  the dose passes ~{trace[8][1]:.0f} mg/day.  So the escalation buys no
  analgesia and the gut keeps every milligram of it.  That -- not a worsening
  bowel score -- is what tolerance asymmetry costs the patient, and it is the
  clinical situation a PAMORA exists to break.""")


def analysis_D():
    hr("D.  WHERE DOES THE AMPLIFICATION COME FROM?  A PAMORA wins ~0.57 of\n"
       "    enteric occupancy and buys a ~2.9-fold change in bowel frequency.\n"
       "    This analysis asks which mechanism carries that gain -- and refutes\n"
       "    the author's prior answer.")
    base, _, _, _ = run_scenario("ref", "none", {})
    print(f"{'naloxegol mg':>13} {'OCC_ant':>8} {'SIGg':>7} {'SBM/wk':>7} "
          f"{'x baseline':>11} {'w4':>7} {'CTT h':>7}")
    pts = []
    for dose in [0.0, 3.125, 6.25, 12.5, 25.0, 50.0, 100.0]:
        r = base if dose == 0 else run_scenario("d", "naloxegol",
                                                dict(pam_dose=dose))[0]
        pts.append((dose, r["OCCg_ant"], r["SBM_wk"], r["w4"]))
        print(f"{dose:13.3f} {r['OCCg_ant']:8.3f} {r['OCCg_ag']:7.3f} "
              f"{r['SBM_wk']:7.2f} {r['SBM_wk']/base['SBM_wk']:11.2f} "
              f"{r['w4']:7.3f} {r['CTT_h']:7.0f}")
    (_, o1, s1, _), (_, o2, s2, _) = pts[3], pts[4]
    print(f"\n  d ln(SBM) / d ln(enteric blockade) between 12.5 and 25 mg = "
          f"{(math.log(s2)-math.log(s1))/(math.log(o2)-math.log(o1)):.2f}")

    print("""
  Enteric MOR does not act on one thing.  It acts on FIVE brakes at once:
  propulsive drive (via cAMP/ACh), segmental tone, inhibitory NO/VIP tone,
  anal sphincter tone, and the rectal urge threshold -- plus a secretion arm
  that feeds the hydration term.  Disable each in turn and re-measure the SAME
  0 -> 25 mg naloxegol comparison.  Whichever brake's removal FLATTENS the
  fold-change is the one carrying the amplification.
""")
    print(f"{'brake disabled':44} {'untreated':>10} {'nalox 25':>9} {'fold':>6}")
    variants = [
        ("(none -- full model)", {}),
        ("hydration -> propulsion loop (GW=GWD=0)", dict(GW=0.0, GWD=0.0)),
        ("rectal urge threshold shift (AVSENS=0)", dict(AVSENS=0.0)),
        ("anal sphincter tone (AANO=0)", dict(AANO=0.0)),
        ("segmental muscle tone (ATONE=0)", dict(ATONE=0.0)),
        ("inhibitory NO/VIP tone (ANO=0)", dict(ANO=0.0)),
        ("secretion suppression (ESECMOR=0)", dict(ESECMOR=0.0)),
        ("cAMP -> ACh -> HAPC drive (EMAXMOR=0)", dict(EMAXMOR=0.0)),
    ]
    rows = []
    for lab, mods in variants:
        p0 = base_params(); p0 = set_antagonist(p0, "none")
        p0.update(OPIOID_STD); p0.update(mods)
        t, Y = simulate(p0, 24 * 84); u = readouts(t, Y, p0)["SBM_wk"]
        p1 = base_params(); p1 = set_antagonist(p1, "naloxegol")
        p1.update(OPIOID_STD); p1.update(mods); p1["pam_dose"] = 25.0
        t, Y = simulate(p1, 24 * 84); a = readouts(t, Y, p1)["SBM_wk"]
        rows.append((lab, u, a, a / max(u, 1e-9)))
        print(f"{lab:44} {u:10.2f} {a:9.2f} {a/max(u,1e-9):6.2f}")
    full = rows[0][3]
    ranked = sorted(rows[1:], key=lambda r: r[3])
    print(f"""
  RESULT -- and it contradicts the hypothesis this analysis was written to
  confirm.  Removing the transit-hydration feedback loop entirely moves the
  fold-response from {full:.2f} to {rows[1][3]:.2f}.  That loop is NOT the
  amplifier.  It sets the absolute LEVEL of stool hydration -- untreated
  frequency drops from {rows[0][1]:.2f} to {rows[1][1]:.2f} when it is removed
  -- but it carries none of the drug RESPONSE.

  Ranked by how much removing them flattens the response:
""")
    for lab, u, a_, f in ranked:
        print(f"    fold {f:5.2f}   {lab}")
    top1, top2 = ranked[0], ranked[1]
    print(f"""
  The two that matter are ANORECTAL, not colonic: the rectal urge threshold
  ({top1[3]:.2f}) and internal anal sphincter tone ({top2[3]:.2f}).  Colonic
  propulsion (the cAMP -> ACh -> HAPC arm) is third at
  {[r for r in ranked if 'EMAXMOR' in r[0]][0][3]:.2f}, and segmental tone,
  inhibitory NO/VIP tone and the secretion arm contribute essentially nothing
  to the fold-change at all.

  So the amplification is not a feedback gain.  It is the product of several
  independent brakes released together by one receptor, and it is dominated by
  the brakes at the END of the colon.  Two consequences worth testing:

    1. A PAMORA should help an OIC patient whose problem is a lost urge and a
       tight sphincter far more than one with true colonic inertia -- the model
       predicts a responder split along anorectal rather than transit
       phenotype, which anorectal manometry could check.
    2. A drug that removes the cause acts on all the brakes at once and so
       MULTIPLIES; a secretagogue or an osmotic agent compensates on one
       downstream arm and so ADDS.  That is the model's structural reason for
       PAMORAs outperforming laxatives by more than their potencies suggest.""")


def analysis_E():
    hr("E.  OSMOLE ARITHMETIC.  PEG 3350 and lactulose are both called\n"
       "    'osmotic laxatives'.  Per gram they are not remotely comparable,\n"
       "    and the reason is a molecular weight and a bacterium.")
    for lab, g, mw, ferm in [("PEG 3350 17 g", 17.0, 3350.0, False),
                             ("lactulose 20 g", 20.0, 342.3, True),
                             ("lactulose 40 g", 40.0, 342.3, True)]:
        mmol = g / mw * 1000.0
        obliged = mmol / 0.290
        amp = 3.6 if ferm else 1.0
        print(f"  {lab:16} {mmol:7.2f} mmol -> {obliged:7.1f} mL obliged; "
              f"fermented x{amp:.1f} -> {obliged*amp:7.1f} mL")
    print("\n  Simulated, in the same OIC patient:")
    print(f"{'arm':28} {'SBM/wk':>7} {'BSFS':>6} {'w4':>7} {'distension':>11} {'rescue/wk':>10}")
    for lab, mods in [("untreated", {}),
                      ("PEG 3350 17 g/d", dict(peg_dose=17.0)),
                      ("PEG 3350 34 g/d", dict(peg_dose=34.0)),
                      ("lactulose 20 g/d", dict(lac_dose=20.0)),
                      ("lactulose 40 g/d", dict(lac_dose=40.0))]:
        r, _, _, _ = run_scenario(lab, "none", mods)
        print(f"{lab:28} {r['SBM_wk']:7.2f} {r['BSFS']:6.2f} {r['w4']:7.3f} "
              f"{r['DIST']:11.2f} {r['RESC_wk']:10.2f}")
    print("\n  PEG is dosed in grams of a 3350 Da polymer, so 17 g is 5.1 mmol and")
    print("  obliges ~17 mL.  Lactulose is 342 Da AND is cleaved by colonic flora")
    print("  into osmotically active fragments.  The distension column is the")
    print("  price: the same fermentation that makes the osmoles makes the gas.")


def analysis_F():
    hr("F.  WHY LUBIPROSTONE FAILED IN METHADONE PATIENTS.\n"
       "    Methadone blocks ClC-2 directly -- a channel block, not a MOR effect.\n"
       "    So the drug and the antidote collide on the same protein.")
    print(f"{'arm':40} {'SBM/wk':>7} {'BSFS':>6} {'w4':>7} {'ClC-2':>7}")
    for lab, mods in [
        ("oxycodone, untreated", {}),
        ("oxycodone + lubiprostone 24 ug bid", dict(lub_dose=24.0, lub_int=12.0)),
        ("oxycodone + linaclotide 145 ug qd", dict(lin_dose=145.0)),
        ("methadone, untreated", dict(METHADONE=1.0)),
        ("methadone + lubiprostone 24 ug bid", dict(METHADONE=1.0, lub_dose=24.0, lub_int=12.0)),
        ("methadone + linaclotide 145 ug qd", dict(METHADONE=1.0, lin_dose=145.0)),
        ("methadone + naldemedine 0.2 mg", dict(METHADONE=1.0)),
    ]:
        drug = "naldemedine" if "naldemedine" in lab else "none"
        m = dict(mods)
        if drug == "naldemedine":
            m["pam_dose"] = 0.2
        r, _, t, Y = run_scenario(lab, drug, m)
        clc2 = Y[-1, IX["CLC2"]]
        print(f"{lab:40} {r['SBM_wk']:7.2f} {r['BSFS']:6.2f} {r['w4']:7.3f} {clc2:7.3f}")
    print("\n  Linaclotide reaches the same CFTR-mediated secretion through GC-C /")
    print("  cGMP and is untouched by the block.  A PAMORA is untouched twice over,")
    print("  because it removes the cause rather than compensating downstream.")


def analysis_G():
    hr("G.  CALIBRATION PROVENANCE -- exactly what was fitted, and to what.\n"
       "    THREE parameters were fitted in this entire model.  Everything else\n"
       "    is a literature value or a structural assumption.")
    print("""
  parameter  drug                fitted to                                 value
  ---------  ------------------  ----------------------------------------  --------
  KIGUT      naloxegol           KODIAC-04/05 25 mg arm, mean SBM/wk 4.3    2.685 nM
  KIGUT      naldemedine         COMPOSE-1/2 0.2 mg arm, mean SBM/wk 5.0    0.00821 nM
  KIGUT      methylnaltrexone    Thomas 2008, placebo-subtracted 4-h        0.795 nM
                                 laxation difference of 33 pp

  Each is an OPERATIONAL enteric potency.  Compare it with the drug's own
  published MOR binding Ki -- a number the model never fitted:
""")
    print(f"  {'drug':22} {'Ki (lit, nM)':>13} {'KIGUT (fitted)':>15} {'ratio':>8}")
    for d, kg in [("naloxegol", 2.6851), ("naldemedine", 0.00821),
                  ("methylnaltrexone_sc", 0.7945)]:
        ki = ANTAGONISTS[d]["KIANT"]
        print(f"  {ANTAGONISTS[d]['ANT_NAME']:22} {ki:13.3f} {kg:15.5f} {kg/ki:8.4f}")
    print("""
  All three need MORE enteric potency than plasma-unbound binding predicts
  (2.8x, 41x, 35x).  The model does NOT explain the spread, and it is worth
  being explicit that it does not:

    - a single shared error in the AGONIST anchor (KIOP) would move all three
      ratios by the SAME factor.  They differ by 15-fold, so that is not it.
    - plasma protein binding does not order them either: naloxegol (fu 0.96)
      needs the smallest correction, but methylnaltrexone (fu 0.885) needs a
      correction as large as naldemedine's (fu 0.065).

  Candidates left open: accumulation of antagonist in the gut wall beyond the
  pre-systemic term modelled here, active uptake into enteric neurons, error in
  the published CL/F values used to set plasma exposure, and the fact that the
  three trials enrolled different populations against endpoints defined over
  different windows.  Any per-drug SBM prediction inherits this uncertainty.
  The SELECTIVITY results (analyses A and B) do not: they are ratios in which
  KIGUT cancels.""")


def analysis_H():
    hr("H.  PEAK-DRIVEN vs EXPOSURE-DRIVEN -- and which ENDPOINT can see it.\n"
       "    Same total daily dose, split differently, read two ways: the chronic\n"
       "    weekly frequency, and the acute probability of laxation within 4 h.")
    # a steady, untreated OIC patient to give the first dose to
    p0 = base_params(); p0 = set_antagonist(p0, "none"); p0.update(OPIOID_STD)
    t, Y = simulate(p0, 24 * 84); y0 = Y[-1].copy()

    def acute(drug, dose):
        p = base_params(); p = set_antagonist(p, drug)
        p.update(OPIOID_STD); p["pam_dose"] = dose; p["pam_int"] = 1e9
        t2, Y2 = simulate(p, 4.0, y0=y0, nout_per_h=8)
        return 1 - math.exp(-(Y2[-1, IX["CUM_SBM"]] - Y2[0, IX["CUM_SBM"]]))

    print(f"{'regimen':44} {'SBM/wk':>7} {'P(lax<4h)':>10}")
    print(f"{'(no antagonist)':44} {'':>7} {acute('none', 0.0):10.3f}")
    for lab, drug, dose, ivl in [
        ("naloxegol 25 mg qd", "naloxegol", 25.0, 24.0),
        ("naloxegol 12.5 mg bid (same 25 mg/d)", "naloxegol", 12.5, 12.0),
        ("naloxegol 6.25 mg qid (same 25 mg/d)", "naloxegol", 6.25, 6.0),
        ("methylnaltrexone 12 mg SC qd", "methylnaltrexone_sc", 12.0, 24.0),
        ("methylnaltrexone 6 mg SC bid (same 12)", "methylnaltrexone_sc", 6.0, 12.0),
        ("methylnaltrexone 3 mg SC qid (same 12)", "methylnaltrexone_sc", 3.0, 6.0),
    ]:
        r, _, _, _ = run_scenario(lab, drug, dict(pam_dose=dose, pam_int=ivl))
        print(f"{lab:44} {r['SBM_wk']:7.2f} {acute(drug, dose):10.3f}")
    print("""
  The chronic column refutes the framing this analysis started with: splitting
  the dose does not punish methylnaltrexone, it mildly HELPS both drugs,
  because a weekly frequency is an average and smoother exposure raises the
  average occupancy.  Cmax-vs-Cavg is invisible to that endpoint.

  The acute column is where the distinction lives.  It is also the endpoint
  methylnaltrexone was actually licensed on (rescue-free laxation within 4 h),
  and the one naloxegol and naldemedine were NOT licensed on -- which, read
  this way, looks less like a marketing choice than like each sponsor choosing
  the endpoint their molecule's exposure profile could win.""")


def analysis_I():
    hr("I.  ALL SCENARIOS")
    print(f"{'scenario':38} {'SBM/wk':>7} {'CSBM':>6} {'BSFS':>6} {'strain':>7} "
          f"{'PAC-SYM':>8} {'pain':>6} {'COWS':>6} {'resc/wk':>8}")
    out = []
    for name, drug, mods in SCENARIOS:
        r, _, _, _ = run_scenario(name, drug, mods)
        out.append(r)
        print(f"{name:38} {r['SBM_wk']:7.2f} {r['CSBM_wk']:6.2f} {r['BSFS']:6.2f} "
              f"{r['STRAIN']:7.2f} {r['PACSYM']:8.2f} {r['PAIN']:6.2f} "
              f"{r['WD']:6.2f} {r['RESC_wk']:8.2f}")
    return out


def _make_pop(n, rng, spread=1.0):
    pop = []
    for _ in range(n):
        pop.append(dict(
            op_dose=float(np.exp(rng.normal(math.log(30.0), 0.45 * spread))),
            KIOP=float(30.0 * np.exp(rng.normal(0.0, 0.30 * spread))),
            ABSX=float(np.exp(rng.normal(0.0, 0.18 * spread))),
            FLUIDX=float(np.exp(rng.normal(0.0, 0.16 * spread))),
            SIN=float(1.46 * np.exp(rng.normal(0.0, 0.22 * spread))),
            AVSENS=float(2.5 * np.exp(rng.normal(0.0, 0.30 * spread))),
            KDEF=float(0.415 * np.exp(rng.normal(0.0, 0.22 * spread))),
        ))
    return pop


def _run_arm(pop, drug, mods, weeks):
    sbms = []
    for ind in pop:
        p = base_params(); p = set_antagonist(p, drug)
        p.update(OPIOID_STD); p.update(ind); p.update(mods)
        try:
            t, Y = simulate(p, 24 * 7 * weeks)
            sbms.append(readouts(t, Y, p)["SBM_wk"])
        except Exception:
            sbms.append(float("nan"))
    a = np.array(sbms, dtype=float)
    return a[~np.isnan(a)]


def analysis_J2(n=60, weeks=8, seed=20260806, spread=3.2):
    hr(f"J2. VIRTUAL POPULATION WITH DISPERSION CALIBRATED TO THE TRIAL\n"
       f"    PLACEBO ARM (n={n}, every parameter CV inflated {spread}-fold).\n"
       "    In J the untreated arm produced 0% responders against a trial\n"
       "    placebo rate of 29%.  Rather than leave that as a bare failure,\n"
       "    ask what dispersion WOULD reproduce it, then re-read the drug arms.")
    rng = np.random.default_rng(seed)
    pop = _make_pop(n, rng, spread)
    arms = [("untreated", "none", {}),
            ("naloxegol 25 mg", "naloxegol", dict(pam_dose=25.0)),
            ("naldemedine 0.2 mg", "naldemedine", dict(pam_dose=0.2))]
    res = {}
    for lab, drug, mods in arms:
        a = _run_arm(pop, drug, mods, weeks)
        res[lab] = a
        print(f"  {lab:26} median SBM/wk {np.median(a):5.2f}   "
              f"responder (>=3/wk) {100*np.mean(a >= 3.0):5.1f}%")
    u = 100 * np.mean(res["untreated"] >= 3.0)
    print()
    for lab, trial_abs, trial_delta in [("naloxegol 25 mg", 44.4, 12.7),
                                        ("naldemedine 0.2 mg", 47.6, 16.0)]:
        r = 100 * np.mean(res[lab] >= 3.0)
        print(f"  {lab:26} model {r:5.1f}% (delta {r-u:+5.1f} pp)   "
              f"vs trial {trial_abs:.1f}% (delta {trial_delta:+.1f} pp)")
    print(f"""
  Reading this honestly: the {spread}-fold inflation is NOT a mechanistic
  finding, it is a measurement of what the model is missing.  Every individual
  CV it implies is implausible on its own -- an opioid-dose CV of
  {0.45*spread:.2f} on the log scale is not a real prescribing distribution.
  What the number says is that the mechanistic parameters varied here carry
  only about 1/{spread:.0f} of the between-patient variance that the trials
  actually contain, and that the rest sits in things this model has no state
  for: adherence, diet, mobility, comorbid anorectal disease, concurrent
  constipating drugs, and the placebo response itself.

  With the dispersion forced to match, the ACTIVE arms are predictions again --
  and they still overshoot.  The model reproduces the typical patient's mean
  SBM/week (that was fitted) but not the shape of the response distribution.
  Do not use this model to predict responder rates.""")
    return res


def analysis_J(n=60, weeks=8, seed=20260806):
    hr(f"J.  VIRTUAL POPULATION (n={n}) -- responder rate, >=3 SBM/week.\n"
       "    NOTE: this model has NO placebo mechanism.  The untreated arm is\n"
       "    mechanistically untreated, so it is NOT comparable to a trial\n"
       "    placebo arm (KODIAC 29%, COMPOSE 34%).  Only the DIFFERENCE between\n"
       "    an active arm and the untreated arm here can be read against the\n"
       "    placebo-subtracted trial difference.")
    rng = np.random.default_rng(seed)
    pop = []
    for _ in range(n):
        pop.append(dict(
            op_dose=float(np.exp(rng.normal(math.log(30.0), 0.45))),
            KIOP=float(30.0 * np.exp(rng.normal(0.0, 0.30))),
            ABSX=float(np.exp(rng.normal(0.0, 0.18))),
            FLUIDX=float(np.exp(rng.normal(0.0, 0.16))),
            SIN=float(1.46 * np.exp(rng.normal(0.0, 0.22))),
            AVSENS=float(2.5 * np.exp(rng.normal(0.0, 0.30))),
            KDEF=float(0.415 * np.exp(rng.normal(0.0, 0.22))),
        ))
    arms = [("untreated", "none", {}),
            ("naloxegol 25 mg", "naloxegol", dict(pam_dose=25.0)),
            ("naldemedine 0.2 mg", "naldemedine", dict(pam_dose=0.2)),
            ("PEG 3350 17 g", "none", dict(peg_dose=17.0)),
            ("lactulose 20 g", "none", dict(lac_dose=20.0)),
            ("naloxegol 25 + lactulose", "naloxegol", dict(pam_dose=25.0, lac_dose=20.0))]
    res = {}
    for lab, drug, mods in arms:
        sbms = []
        for ind in pop:
            p = base_params(); p = set_antagonist(p, drug)
            p.update(OPIOID_STD); p.update(ind); p.update(mods)
            try:
                t, Y = simulate(p, 24 * 7 * weeks)
                sbms.append(readouts(t, Y, p)["SBM_wk"])
            except Exception:
                sbms.append(float("nan"))
        a = np.array(sbms, dtype=float)
        a = a[~np.isnan(a)]
        res[lab] = a
        print(f"  {lab:26} median SBM/wk {np.median(a):5.2f}   "
              f"responder (>=3/wk) {100*np.mean(a >= 3.0):5.1f}%")
    u = 100 * np.mean(res["untreated"] >= 3.0)
    print()
    for lab, trial in [("naloxegol 25 mg", 12.7), ("naldemedine 0.2 mg", 16.0)]:
        d = 100 * np.mean(res[lab] >= 3.0) - u
        print(f"  {lab:26} model delta {d:+5.1f} pp   vs trial placebo-subtracted "
              f"{trial:+.1f} pp")
    a, b = res["naloxegol 25 mg"], res["lactulose 20 g"]
    c = res["naloxegol 25 + lactulose"]
    da, db, dc = (100*np.mean(a>=3)-u), (100*np.mean(b>=3)-u), (100*np.mean(c>=3)-u)
    print(f"\n  combination arithmetic (responder pp over untreated):")
    print(f"    naloxegol alone {da:+.1f} pp;  lactulose alone {db:+.1f} pp;  "
          f"sum {da+db:+.1f} pp;  ACTUAL combination {dc:+.1f} pp")
    print(f"    -> the two act on DIFFERENT arms of the same loop (receptor vs")
    print(f"       water), so the loop amplifies them together rather than twice.")
    print(f"\n  PEG 3350 17 g arm: {100*np.mean(res['PEG 3350 17 g']>=3)-u:+.1f} pp over "
          f"untreated -- see the PEG discussion in analysis E.")
    return res


def main():
    print(__doc__)
    print(f"\nstate variables: {NST}")
    analysis_G()
    analysis_A()
    analysis_B()
    analysis_C()
    analysis_D()
    analysis_E()
    analysis_F()
    analysis_H()
    analysis_I()
    res = analysis_J()
    analysis_J2()
    out = {k: dict(median=float(np.median(v)),
                   responder_pct=float(100 * np.mean(v >= 3.0)),
                   n=int(v.size)) for k, v in res.items()}
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "oic_population_results.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\n[wrote oic_population_results.json]")


if __name__ == "__main__":
    main()
