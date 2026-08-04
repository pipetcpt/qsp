#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 PATENT DUCTUS ARTERIOSUS OF PREMATURITY  --  QSP REFERENCE IMPLEMENTATION
================================================================================

 Independent Python re-implementation of ``pda_mrgsolve_model.R``.  Its purpose
 is NOT to be the deliverable model -- the mrgsolve file is -- but to be a
 second pair of hands: every structural equation is written twice, once here in
 SciPy/LSODA and once in mrgsolve's C++, and the two are compared numerically.
 Divergence means one of them is wrong, and in this project that comparison has
 repeatedly been the thing that found the bug.

 --------------------------------------------------------------------------
 WHY THIS DISEASE NEEDS A QSP MODEL AT ALL
 --------------------------------------------------------------------------
 The ductus arteriosus is not a lesion.  It is a normal fetal vessel held open
 by prostaglandin E2 (placental production + immature pulmonary clearance)
 acting on EP4 -> Gs -> cAMP in ductal smooth muscle.  After birth the
 placental PGE2 source disappears, pulmonary 15-hydroxyprostaglandin
 dehydrogenase (15-PGDH) starts destroying what is left, and rising PaO2 drives
 the smooth muscle to contract.  In a term infant that sequence closes the duct
 in 24-72 h and then *remodels* it permanently.  In a 25-week infant it often
 does neither.

 Every therapeutic controversy in this field is quantitative, not qualitative:

  1. COX inhibitors close the duct.  They do not improve the outcomes the duct
     is blamed for.  Baby-OSCAR (2024) closed ducts with early ibuprofen and
     got 69.4% vs 63.5% death-or-moderate/severe-BPD -- numerically WORSE.
     BeNeDuctus (2023) found expectant management non-inferior.  A model that
     cannot reproduce "the drug works and the patient does not benefit" is not
     modelling this disease.

  2. Almost every preterm duct closes eventually without any drug.  Semberova
     (2017) followed 297 infants <27 wk with a non-intervention policy: 94-98%
     closed spontaneously before discharge, median 71 days in the <26 wk group.
     So the treatment effect on any fixed-time closure endpoint is bounded
     above by the untreated non-closure fraction at that time -- a competing
     -risk ceiling that shrinks as the endpoint moves later.  This model makes
     that ceiling an output rather than an assumption.

  3. The three drugs are not three doses of one drug.  Indomethacin and
     ibuprofen compete with arachidonate in the cyclooxygenase channel.
     Acetaminophen reduces the ferryl-protoporphyrin radical at the physically
     separate PEROXIDASE site, and is therefore competitive with peroxide, not
     with substrate.  One structural consequence: acetaminophen's potency is a
     function of ductal peroxide tone, so it should degrade under
     chorioamnionitis/sepsis while ibuprofen does not.  That is a falsifiable
     prediction produced by writing the enzyme correctly.

  4. A constricted duct is not a closed duct.  Permanent closure requires
     ductal WALL hypoxia (lumen occlusion cuts off luminal O2 diffusion) ->
     VEGF/TGF-beta -> neointimal cushion.  The preterm ductal wall is thin
     enough to stay oxygenated by diffusion even when constricted, so
     remodelling never triggers and the duct reopens when drug clears.  This
     one mechanism explains the ~25% reopening rate in <28 wk infants and its
     near-absence at term, and it makes closure BISTABLE rather than graded.

  5. Efficacy and toxicity are the same molecular event in another organ.
     Renal PGE2 maintains GFR in the immature kidney; gut mucosal PGE2
     maintains the mucosal barrier; platelet COX-1 makes thromboxane.  You
     cannot separate them by dose -- only by organ-specific prostanoid
     dependence and by the non-COX vasoconstriction that indomethacin has and
     ibuprofen does not.

 --------------------------------------------------------------------------
 WHAT IS FITTED AND WHAT IS PREDICTED  (read this before trusting a number)
 --------------------------------------------------------------------------
 FITTED -- 9 parameters, 11 targets, in three separable stages:
   * NET50, NETW, TAUSYN0, KTAUSYN, KINVGA, KEP4, KREMOD
       jointly to 8 gestational spontaneous-closure times AND 3 treated
       closure times.  Jointly is essential: fitting the threshold to
       untreated closure alone puts it beyond the reach of any drug.
   * KI_IBU, KI_IND, IC50_APAP -> comparable day-7 closure across the drugs
   * H_BPD0 (level) and B_BPD_COX (arm separation) -> Baby-OSCAR, both arms.
       B_BUR is deliberately NOT fitted, which is what makes the 16-point
       burden-versus-harm gap an interpretable quantity rather than an
       artefact of a free parameter.
 PREDICTED (no parameter touched afterwards):
   * reopening, its GA dependence, and the early-vs-late corollary
   * loss of acetaminophen effect under high peroxide tone
   * the failure of LATE and echo-confirmed treatment to reduce burden at all
   * BeNeDuctus, PDA-TOLERATE, TIPP
   * renal / cerebral / mesenteric flow separation between the two NSAIDs
   * the complete structural failure of any dose below ~25.5 weeks

 THE MODEL CONTRADICTS THE HYPOTHESIS IT WAS BUILT TO TEST.  It was written
 expecting to show that targeted late treatment captures most of the benefit
 for a fraction of the exposure.  It says the opposite: PDA burden accrues
 almost entirely in the FIRST WEEK, so treating at day 7 or day 10 removes
 essentially no burden (-9 and -0 units against -227 for day 2) while
 delivering the COX exposure in full.  Targeting needs early PREDICTION, not
 late confirmation.

 Units: time h (postnatal), concentrations mg/L (total) and uM (unbound),
 flows mL/min/kg, pressures mmHg, resistances mmHg/(mL/min/kg).

 Run:  python3 pda_reference_model.py            (all scenarios + calibration)
       python3 pda_reference_model.py --quick     (skip population runs)
================================================================================
"""

import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

# ==============================================================================
#  SECTION 1.  STATE VECTOR
# ==============================================================================
# 50 ODEs.  Names are shared verbatim with the mrgsolve model so that the two
# implementations can be diffed compartment by compartment.

SNAMES = [
    # ---- drug PK (0-8) -------------------------------------------------------
    "IBU1",    # 0  ibuprofen central          (mg/kg)
    "IBU2",    # 1  ibuprofen peripheral       (mg/kg)
    "IBUE",    # 2  ibuprofen ductal effect site (mg/L equivalent)
    "IND1",    # 3  indomethacin central       (mg/kg)
    "IND2",    # 4  indomethacin peripheral    (mg/kg)
    "INDE",    # 5  indomethacin effect site (slowly-reversible COX binding)
    "APAPG",   # 6  acetaminophen gut depot    (mg/kg)
    "APAP1",   # 7  acetaminophen central      (mg/kg)
    "APAP2",   # 8  acetaminophen peripheral   (mg/kg)
    # ---- hepatic safety (9-10) ----------------------------------------------
    "GSH",     # 9  hepatic glutathione (fraction of normal)
    "ALT",     # 10 alanine aminotransferase (U/L)
    # ---- prostanoid system (11-18) ------------------------------------------
    "PGE2D",   # 11 ductal tissue PGE2 (fraction of fetal baseline)
    "PGE2P",   # 12 circulating PGE2 (pg/mL)
    "PGDH",    # 13 pulmonary 15-PGDH capacity (fraction of term)
    "CAMP",    # 14 ductal SMC cAMP (fraction of fetal baseline)
    "NOD",     # 15 ductal NO tone (fraction of fetal baseline)
    "EP4",     # 16 ductal EP4 receptor density (fraction of fetal baseline)
    "PEROX",   # 17 ductal wall peroxide tone (fraction of normal)
    "CASMC",   # 18 ductal SMC cytosolic Ca (fraction of maximal)
    # ---- ductus structure (19-23) -------------------------------------------
    "TONE",    # 19 ductal smooth-muscle tone (0-1)
    "REMOD",   # 20 anatomic remodelling / neointimal cushion (0-1)
    "WALLTH",  # 21 ductal wall thickness index (fraction of term)
    "VEGFD",   # 22 ductal wall VEGF (fraction of normal)
    "TGFBD",   # 23 ductal wall TGF-beta1 (fraction of normal)
    # ---- haemodynamics (24-32) ----------------------------------------------
    "PVR",     # 24 pulmonary vascular resistance
    "SVR",     # 25 systemic vascular resistance
    "PLA",     # 26 left atrial pressure (mmHg)
    "LVDIL",   # 27 LV dilation / volume-load index (fraction above baseline)
    "EVLW",    # 28 extravascular lung water (mL/kg)
    "CRS",     # 29 respiratory system compliance (mL/cmH2O/kg)
    "PDABUR",  # 30 PDA burden = integral of significant shunt (mL/kg/min * day)
    "MESDEF",  # 31 cumulative mesenteric flow deficit (fraction * day)
    "CERDEF",  # 32 cumulative cerebral flow deficit (fraction * day)
    # ---- kidney (33-37) -----------------------------------------------------
    "GFR",     # 33 glomerular filtration rate (mL/min/kg)
    "SCR",     # 34 serum creatinine (mg/dL)
    "UOCUM",   # 35 cumulative urine output (mL/kg)
    "FLUID",   # 36 net fluid balance (mL/kg)
    "PGE2K",   # 37 renal medullary/afferent PGE2 (fraction of normal)
    # ---- gut, platelet, bilirubin (38-42) -----------------------------------
    "PGE2G",   # 38 gut mucosal PGE2 (fraction of normal)
    "TXA2",    # 39 platelet thromboxane A2 capacity (fraction of normal)
    "BTIME",   # 40 bleeding time index (fraction of normal)
    "TBILI",   # 41 total serum bilirubin (mg/dL)
    "BFREE",   # 42 free (unbound) bilirubin index (fraction of normal)
    # ---- lung development and outcome hazards (43-49) -----------------------
    "ALV",     # 43 alveolarisation index (fraction of gestation-matched)
    "VENTIX",  # 44 cumulative respiratory severity (support index * day)
    "COXIX",   # 45 cumulative systemic COX-inhibition exposure (fraction*day)
    "HBPD",    # 46 cumulative hazard, death or moderate/severe BPD
    "HNEC",    # 46 cumulative hazard, NEC stage >= 2
    "HIVH",    # 47 cumulative hazard, severe (grade 3-4) IVH
    "HSIP",    # 48 cumulative hazard, spontaneous intestinal perforation
    "HDEATH",  # 49 cumulative hazard, death before 36 wk PMA
]
IX = {n: i for i, n in enumerate(SNAMES)}
NST = len(SNAMES)


# ==============================================================================
#  SECTION 2.  PARAMETERS
# ==============================================================================

def default_params():
    p = {}

    # ------------------------------------------------------------------ patient
    p["GA"] = 26.0          # gestational age at birth (weeks)
    p["BW"] = 0.80          # birth weight (kg)
    p["PAO2"] = 55.0        # arterial PO2 (mmHg) -- typical targeted range
    p["SEPSIS"] = 0.0       # 0/1 chorioamnionitis / early-onset sepsis flag
    p["HCORT"] = 0.0        # 0/1 concomitant early hydrocortisone
    p["ANTESTER"] = 1.0     # 0/1 complete antenatal corticosteroid course

    # --------------------------------------------------- ibuprofen PK (per kg)
    # Aranda 1997; Van Overmeire 2001; Hirt 2008 popPK.  Reported preterm
    # values are wide AND mutually inconsistent: Aranda's own CL 2.06 mL/h/kg
    # with Vd 0.062 L/kg implies t1/2 21 h, not the 30.5 h reported alongside
    # them, and a Vd that small implies a Cmax of 160 mg/L after 10 mg/kg,
    # which no assay report supports.  We therefore adopt a SELF-CONSISTENT
    # set inside the reported envelope (t1/2 20-43 h, CL 1.4-4.5 mL/h/kg,
    # Vss 0.06-0.25 L/kg) and verify the simulated t1/2 and Cmax land in it.
    # Only the UNBOUND concentration reaches COX, so the absolute total-
    # concentration scale is absorbed into the fitted Ki anyway.
    p["V1_IBU"] = 0.140     # L/kg
    p["V2_IBU"] = 0.080     # L/kg
    p["Q_IBU"] = 0.012      # L/h/kg
    p["CL_IBU0"] = 0.00420  # L/h/kg at PNA 0 -> t1/2 ~ 36 h
    p["CLMAT_IBU"] = 0.20   # per day fractional increase in CL (maturation)
    p["FU_IBU"] = 0.0120    # unbound fraction (preterm; 6-12x adult, and
                            # further raised by hypoalbuminaemia + bilirubin)
    p["MW_IBU"] = 206.3     # g/mol
    p["KE0_IBU"] = 0.35     # 1/h ductal effect-site equilibration

    # ------------------------------------------------ indomethacin PK (per kg)
    p["V1_IND"] = 0.250     # L/kg
    p["V2_IND"] = 0.120     # L/kg
    p["Q_IND"] = 0.025      # L/h/kg
    p["CL_IND0"] = 0.01160  # L/h/kg  -> t1/2 ~ 22 h
    p["CLMAT_IND"] = 0.12   # per day
    p["FU_IND"] = 0.010     # unbound fraction
    p["MW_IND"] = 357.8
    p["KE0_IND"] = 0.12     # slow: models slowly-reversible tight COX binding

    # ----------------------------------------------- acetaminophen PK (per kg)
    p["V1_APAP"] = 0.60     # L/kg
    p["V2_APAP"] = 0.40     # L/kg
    p["Q_APAP"] = 0.30      # L/h/kg
    p["CL_APAP0"] = 0.120   # L/h/kg -> t1/2 ~ 6 h (preterm; adult ~0.35)
    p["CLMAT_APAP"] = 0.06  # per day
    p["KA_APAP"] = 0.60     # 1/h enteral absorption
    p["F_APAP"] = 0.85      # enteral bioavailability
    p["FU_APAP"] = 0.85
    p["MW_APAP"] = 151.2

    # ---- acetaminophen hepatic safety (NAPQI / glutathione) -----------------
    p["FNAPQI"] = 0.030     # fraction of clearance via CYP oxidation (neonate)
    p["KGSH"] = 0.020       # 1/h glutathione resynthesis
    p["GSHCAP"] = 1.0
    p["KALT"] = 0.60        # ALT production gain per unit unquenched NAPQI
    p["KALTOUT"] = 0.0173   # 1/h ALT elimination (t1/2 40 h)
    p["ALT0"] = 20.0        # U/L baseline

    # -------------------------------------------------- COX enzymology (ductus)
    # Unbound-drug potencies at the ductal cyclooxygenase.  Indomethacin and
    # ibuprofen are competitive in the ARACHIDONATE CHANNEL; acetaminophen acts
    # at the PEROXIDASE site and is written multiplicatively with them.
    p["KI_IBU"] = 0.220     # uM, unbound, arachidonate-channel (apparent,
                            # in-vivo, at the ductal COX -- FITTED)
    p["KI_IND"] = 0.00083   # uM, unbound (FITTED)
    p["IC50_APAP"] = 6.0   # uM, unbound, at basal peroxide tone (FITTED)
    p["KPEROX"] = 1.00      # peroxide tone that doubles acetaminophen IC50

    # -------------------------------------------------- prostanoid dynamics
    # Local ductal PGE2 synthesis is not static: ductal COX-2 expression and
    # prostanoid production fall postnatally, and they fall FASTER at higher
    # gestational age.  Without this term the untreated duct never closes,
    # which is the single most important thing the model must get right --
    # spontaneous closure is the competing process that bounds every trial.
    p["KSYN_PGE2D"] = 0.90  # 1/h ductal PGE2 synthesis (normalised)
    p["KDEG_PGE2D"] = 0.90  # 1/h ductal PGE2 degradation
    p["FSYN_INF"] = 0.18    # asymptotic fraction of fetal synthesis capacity
    p["TAUSYN0"] = 200.0    # h, postnatal decay time constant at 26 wk (FITTED)
    p["KTAUSYN"] = 0.140    # per wk: faster involution at higher GA (FITTED)
    p["KIN_PGE2D"] = 0.055  # ductal uptake gain from plasma (normalised)
    p["PGE2P0"] = 700.0     # pg/mL fetal circulating PGE2
    p["RPROD_PGE2P"] = 22.0 # pg/mL/h postnatal extra-placental production
    p["KPUL_PGE2"] = 0.055  # 1/h pulmonary clearance gain (x PGDH x flow)
    p["KOTH_PGE2"] = 0.010  # 1/h non-pulmonary clearance
    p["PGDH0_GA"] = 22.0    # wk at which pulmonary 15-PGDH ~ 0
    p["PGDH_SL"] = 0.070    # per wk gestational gain
    p["KPGDH"] = 0.0060     # 1/h postnatal 15-PGDH induction
    p["KCAMP"] = 2.0        # 1/h cAMP turnover
    p["KNO"] = 0.10         # 1/h NO tone turnover
    p["NOGA"] = 0.055       # per wk decline of ductal NO dependence with GA
    p["KNODEC"] = 0.0300    # 1/h postnatal loss of ductal NO tone (t1/2 ~23 h)
    p["KEP4"] = 0.00600     # 1/h postnatal EP4 down-regulation at 26 wk
    p["EP4GA"] = 0.075      # per wk lower EP4 density with advancing GA
    p["EP4INF"] = 0.30      # asymptotic EP4 fraction
    # The whole postnatal involution PROGRAMME -- EP4 loss, NO loss, contractile
    # maturation, prostanoid-synthesis involution -- runs slower the more
    # immature the infant.  One shared gestational scaling reproduces the very
    # wide spread of spontaneous closure times (days at term, >10 weeks at
    # 24 wk) without a separate parameter per gestation.
    p["KINVGA"] = 0.140     # per wk acceleration of involution kinetics (FITTED)
    p["KPEROXIN"] = 0.030   # 1/h peroxide turnover
    p["PEROXSEP"] = 3.2     # peroxide tone multiplier with sepsis/chorio

    # -------------------------------------------------- ductal contraction
    # Tone is written as a SIGMOID of a net constrictor-minus-dilator drive
    # rather than as a product of gains.  A product form makes the achievable
    # tone the arithmetic product of every immaturity factor, so a 26-week
    # duct could never close at any drug exposure -- which is wrong: it can,
    # and does, in ~55% of treated infants.  The sigmoid puts the gestational
    # limit where it physically belongs instead: in TMAXGA, the maximal
    # occlusion a thin immature muscular wall can generate.
    p["TMAX_GA50"] = 21.5   # wk, half-maximal achievable occlusion
    p["TMAX_SL"] = 0.55     # per wk
    # The weights matter structurally, not just numerically.  An earlier
    # version had GMAT 0.42 and slow NO decay, so most of the postnatal rise
    # in net drive came from contractile maturation and loss of NO -- two
    # terms NO DRUG TOUCHES.  Fitting NET50 to untreated closure times then
    # placed the threshold above anything PGE2 suppression could reach, and
    # every drug in the model became inert (peak tone 0.02 at 92% ductal COX
    # inhibition).  The fix is not a bigger dose: it is recognising that the
    # dominant driver of postnatal ductal constriction IS the decline of the
    # PGE2/EP4 drive, which is exactly the axis the drugs act on.
    p["GO2"] = 1.00         # weight of O2-driven constriction in net drive
    p["GMAT"] = 0.10        # weight of postnatal contractile maturation
    p["GPGE2"] = 1.00       # weight of PGE2/EP4 relaxation
    p["GNO"] = 0.35         # weight of NO relaxation
    p["NET50"] = 0.350      # net drive at half-maximal tone (FITTED)
    p["NETW"] = 0.0600       # steepness of the net-drive->tone map (FITTED)
    p["TAU_MAT"] = 300.0    # h, postnatal contractile maturation at 26 wk
    p["KTONE"] = 0.28       # 1/h tone equilibration
    p["P50_O2"] = 30.0      # mmHg PaO2 for half-maximal O2-driven tone at term
    p["P50_GA"] = 0.130     # per wk: lower GA -> higher P50 (less O2-sensitive)
    p["EMAX_PGE2"] = 0.97   # maximal relaxation attributable to PGE2/EP4
    p["EC50_PGE2"] = 0.300  # normalised PGE2*EP4 for half-maximal relaxation
    p["EMAX_NO"] = 0.55     # maximal relaxation attributable to NO
    p["EC50_NO"] = 0.45
    p["ALPHA_D"] = 0.965    # fraction of diameter abolished at full tone
    p["DMAX_GA"] = 0.135    # mm per wk: larger duct at higher GA
    p["DMAX_INT"] = -1.35   # mm intercept -> 26 wk duct 2.16 mm
    p["DCLOSE"] = 0.30      # mm, diameter below which shunt is trivial

    # -------------------------------------------------- remodelling (bistable)
    # Kajino 2002 (Pediatr Res) is the key experiment: after constriction it is
    # VASA VASORUM HYPOPERFUSION that produces medial hypoxia and hence
    # anatomic remodelling.  So the media has two oxygen supplies with
    # opposite thickness dependence:
    #   * DIFFUSION from the lumen and adventitial surface -- adequate when the
    #     wall is THIN, falling off with the square of diffusion distance;
    #   * VASA VASORUM -- present only when the wall is THICK, and ABOLISHED
    #     by the constriction itself.
    # A thin preterm wall therefore stays oxygenated while constricted (no
    # hypoxia -> no cushion -> reopening on drug washout); a thick term wall
    # loses its vasa and goes profoundly hypoxic (permanent closure).  Because
    # the wall also thickens postnatally, the same duct constricted at day 20
    # remodels where it would not have at day 2 -- so the model predicts that
    # EARLY treatment reopens more often than late treatment, which is what is
    # observed.  None of this was fitted to reopening data.
    p["WALLTH_GA"] = 0.075  # per wk ductal wall thickness
    p["WALLTH_INT"] = -1.20
    p["KWALL"] = 0.0025     # 1/h postnatal wall thickening
    p["FLUM_O2"] = 0.80     # diffusional O2 supply weight (thin-wall dominant)
    p["FVASA_O2"] = 0.20    # vasa vasorum O2 supply weight (thick-wall only)
    p["KDIFF"] = 2.60       # inverse diffusion length (per unit wall index)
    p["O2CRIT"] = 16.0      # mmHg wall PO2 below which hypoxic signalling fires
    p["KVEGF"] = 0.045      # 1/h VEGF induction per mmHg hypoxia
    p["KVEGFOUT"] = 0.030
    p["KTGFB"] = 0.020
    p["KTGFBOUT"] = 0.015
    p["KREMOD"] = 0.00060    # 1/h remodelling gain
    p["KUNREMOD"] = 0.00004 # 1/h (near-irreversible)
    # REMODELLING COMPETENCE is itself gestationally immature, and this term is
    # what reproduces the sharpest feature in the data.  Semberova's medians
    # jump 5.5-fold -- 71 d below 26+0 wk versus 13 d at 26+0-27+6 -- across
    # about one week of gestation.  Above the occlusion threshold (~25.5 wk)
    # tone alone closes the duct, which is fast; below it, closure REQUIRES a
    # neointimal cushion.  Without this factor the model built that cushion
    # almost as fast in a 24-week duct as in a 30-week one and closed it in
    # ~33 d instead of ~85, and neither KREMOD nor KWALL fixed it (scanning
    # KREMOD over 6.5-fold moved closure only 30->38 d, because remodelling is
    # not rate-limiting once hypoxia is on -- the missing physiology was the
    # duct's CAPACITY to remodel at all).  The very preterm ductus has sparse
    # vasa vasorum and immature endothelial/SMC migration and hyaluronan
    # machinery, so it cannot build the cushion.
    p["KREMGA50"] = 25.4    # wk, half-maximal remodelling competence
    p["KREMGA_SL"] = 1.20   # per wk

    # -------------------------------------------------- haemodynamics
    p["PVR0"] = 0.115       # mmHg/(mL/min/kg) at birth
    p["PVRINF"] = 0.015     # asymptotic
    p["KPVR"] = 0.022       # 1/h postnatal PVR fall
    p["PVR_EVLW"] = 0.00120 # PVR rise per mL/kg lung water
    p["SVR0"] = 0.165       # -> mean BP 30 mmHg at Qs 180 mL/min/kg
    p["KSVR"] = 0.004
    p["SVRINF"] = 0.205
    p["QS_TARGET"] = 180.0  # mL/min/kg desired systemic flow
    p["QMAX_LV"] = 400.0    # mL/min/kg maximal LV output (preterm reserve)
    p["QMAX_GA"] = 20.0     # mL/min/kg per wk of GA added to LV reserve
    p["PLA0"] = 4.0         # mmHg
    p["KPLA"] = 0.35        # 1/h
    p["PLA_Q"] = 0.022      # mmHg per mL/min/kg of pulmonary flow
    # Ductal resistance ~ 1/d^4 (Poiseuille for a short tube).  KRDUCT is set
    # so that a 2 mm duct gives R ~ 0.07, i.e. ~200 mL/min/kg of shunt at a
    # 14 mmHg mean gradient, and a 1 mm duct gives a trivial 13 mL/min/kg.
    # The 4th power is what makes "moderate-to-large >=1.5 mm" a real
    # threshold rather than an arbitrary cut-off.
    p["KRDUCT"] = 1.10      # mmHg/(mL/min/kg) * mm^4  (ductal resistance gain)
    p["RDUCT_MIN"] = 0.0060 # non-restrictive floor
    p["KLVDIL"] = 0.0035    # 1/h LV dilation per unit shunt fraction
    p["KLVDILOUT"] = 0.010
    p["PDIA0"] = 26.0       # mmHg baseline diastolic pressure
    p["KDIA"] = 0.040       # mmHg fall per mL/min/kg of diastolic runoff
    p["QSIG"] = 100.0       # mL/min/kg shunt defining "haemodynamically significant"

    # -------------------------------------------------- lung water / mechanics
    p["KEVLW"] = 0.0035     # mL/kg/h per mL/min/kg of excess pulmonary flow
    p["KEVLWOUT"] = 0.025   # 1/h lymphatic clearance
    p["EVLW0"] = 6.0        # mL/kg
    p["CRS0"] = 0.75        # mL/cmH2O/kg at 26 wk
    p["CRS_GA"] = 0.10      # per wk
    p["KCRS"] = 0.030       # 1/h compliance recovery
    p["CRS_EVLW"] = 0.030   # compliance loss per mL/kg lung water

    # -------------------------------------------------- kidney
    p["GFR0"] = 0.62        # mL/min/kg at 26 wk
    p["GFR_GA"] = 0.115     # per wk
    p["KGFRMAT"] = 0.0045   # 1/h postnatal GFR maturation
    p["GFR_PGE2"] = 0.52    # fraction of GFR that is renal-PGE2 dependent
    p["KPGE2K"] = 0.25      # 1/h renal PGE2 turnover
    p["KI_IBU_K"] = 2.20    # uM unbound (renal COX-2 slightly less sensitive)
    p["KI_IND_K"] = 0.0150
    p["IC50_APAP_K"] = 260.0 # renal peroxide tone is high -> weak APAP effect
    # Creatinine: neonatal SCr starts at the maternal value and falls over
    # 1-2 weeks; a drug-induced GFR fall shows up as a blunted fall or a rise.
    p["CRPROD"] = 0.00360   # mg/dL per h creatinine production
    p["KCREL"] = 0.01290    # per h per (mL/min/kg) elimination gain
    p["SCR0"] = 0.95        # mg/dL at birth (maternal)
    p["KUO"] = 3.80         # mL/kg/h urine per mL/min/kg GFR
    p["FLUIDIN"] = 3.5      # mL/kg/h total fluid intake

    # -------------------------------------------------- gut / platelet / bili
    p["KPGE2G"] = 0.20      # 1/h gut PGE2 turnover
    p["KI_IBU_G"] = 1.80
    p["KI_IND_G"] = 0.0130
    p["IC50_APAP_G"] = 200.0
    p["KTXA2"] = 0.010      # 1/h platelet pool turnover (platelet lifespan)
    p["KI_IBU_PLT"] = 1.10  # platelet COX-1
    p["KI_IND_PLT"] = 0.0060
    p["KBTIME"] = 0.15
    p["TBILI0"] = 4.0       # mg/dL
    p["KBILIIN"] = 0.16     # mg/dL/h production
    p["KBILIOUT"] = 0.030   # 1/h conjugation+phototherapy
    p["BILIDISP_IBU"] = 0.024 # free-bilirubin rise per mg/L total ibuprofen

    # -------------------------------------------------- non-COX vasoconstriction
    # The single parameter that separates the two NSAIDs clinically.
    # Indomethacin reduces cerebral, mesenteric and renal flow velocity acutely
    # (Patel 2000; Mosca 1997; Pezzati 1999); ibuprofen essentially does not.
    p["VC_IND"] = 0.30      # maximal extra fractional flow reduction
    p["VC_IBU"] = 0.030
    p["VC50_IND"] = 0.35    # mg/L total indomethacin for half-maximal effect
    p["VC50_IBU"] = 30.0    # mg/L total ibuprofen

    # -------------------------------------------------- lung development
    p["ALV0"] = 1.0
    p["KALV"] = 0.00060     # 1/h alveolarisation
    p["ALV_VENT"] = 0.55    # arrest per unit ventilator severity
    p["ALV_EVLW"] = 0.018   # arrest per mL/kg lung water

    # -------------------------------------------------- outcome hazards
    # Hazards accumulate to 36 wk PMA (the BPD endpoint).  Only B_GA, B_BUR and
    # B_VENT are free and fitted to the Baby-OSCAR arms; every other
    # coefficient is fixed a priori from the epidemiology, because fitting ten
    # coefficients to two proportions would be unidentifiable theatre.
    # THE DRUG-HARM TERM, and why it has to exist.
    # The burden mechanism alone makes early ibuprofen BENEFICIAL: it lowers the
    # shunt, lowers PDA burden, lowers the BPD hazard.  Baby-OSCAR found the
    # opposite -- 69.2% vs 63.5%, and death 13.6% vs 10.3%.  So EITHER burden
    # barely matters, OR the drug harms something the burden term does not see.
    # This model takes the second option and names it: COX-2-derived
    # prostaglandins are required for normal alveolar septation, and transient
    # renal impairment adds fluid to an already wet lung.  B_BPD_COX is the
    # magnitude of that harm per unit of cumulative systemic COX inhibition, and
    # it is FITTED to the arm SEPARATION while B_BUR is fixed a priori.  The
    # model therefore answers a question the trial cannot: GIVEN the burden
    # reduction the drug actually achieves, how large must a competing harm be
    # to cancel it?  That number is the useful output here, not the fit itself.
    p["B_BPD_COX"] = 0.5996  # per unit cumulative COX exposure [FITTED]
    p["H_BPD0"] = 0.000130608  # 1/h baseline
    p["B_GA"] = 0.335       # per wk below 28 (FITTED)
    p["B_BUR"] = 0.00055    # per unit PDA burden (FITTED)
    p["B_VENT"] = 0.052     # per unit cumulative ventilator severity (FITTED)
    p["H_NEC0"] = 0.0000090
    p["B_NEC_MES"] = 3.6    # per unit cumulative mesenteric deficit
    p["B_NEC_PGE2G"] = 1.5  # per unit gut PGE2 suppression
    p["H_IVH0"] = 0.000320
    p["B_IVH_GA"] = 0.40
    p["B_IVH_BT"] = 0.55    # per unit bleeding-time prolongation
    p["B_IVH_CER"] = 1.10   # per unit cerebral flow instability
    p["IVH_IND_PROT"] = 0.42 # germinal-matrix maturation effect of indomethacin
    p["H_SIP0"] = 0.0000075
    p["B_SIP_IND"] = 2.6
    p["B_SIP_HC"] = 2.2
    p["H_DTH0"] = 0.0000260
    p["B_DTH_GA"] = 0.45
    p["B_DTH_QS"] = 1.30    # per unit systemic flow deficit
    p["IVH_WINDOW"] = 96.0  # h -- severe IVH essentially occurs in the first 4 d

    # -------------------------------------------------- dosing infusion rates
    # Continuous-infusion rates.  These two exist only in the Python
    # implementation: mrgsolve expresses infusions through the `rate` column of
    # the event object, so there is no matching $PARAM entry.  They are the only
    # two entries the Python/R parameter diff is allowed to report as
    # Python-only; anything else appearing there is a transcription bug.
    p["RIN_IBU"] = 0.0      # mg/kg/h ibuprofen infusion (continuous scenarios)
    p["RIN_IND"] = 0.0
    return p


# ==============================================================================
#  SECTION 3.  ALGEBRAIC LAYER
# ==============================================================================

def derived(t, y, p):
    """All algebraic quantities.  Returns dict; shared with the ODE routine so
    that outputs and derivatives can never disagree."""
    d = {}
    PNAd = t / 24.0
    d["PNAd"] = PNAd
    GA = p["GA"]

    # ---------------------------------------------------------------- drug conc
    d["CIBU"] = max(0.0, y[IX["IBU1"]]) / p["V1_IBU"]        # mg/L total
    d["CIND"] = max(0.0, y[IX["IND1"]]) / p["V1_IND"]
    d["CAPAP"] = max(0.0, y[IX["APAP1"]]) / p["V1_APAP"]
    # unbound, uM
    d["UIBU"] = d["CIBU"] * p["FU_IBU"] * 1000.0 / p["MW_IBU"]
    d["UIND"] = d["CIND"] * p["FU_IND"] * 1000.0 / p["MW_IND"]
    d["UAPAP"] = d["CAPAP"] * p["FU_APAP"] * 1000.0 / p["MW_APAP"]
    # ductal effect-site unbound concentrations (uM)
    d["UIBUE"] = max(0.0, y[IX["IBUE"]]) * p["FU_IBU"] * 1000.0 / p["MW_IBU"]
    d["UINDE"] = max(0.0, y[IX["INDE"]]) * p["FU_IND"] * 1000.0 / p["MW_IND"]

    # ------------------------------------------------------- COX inhibition
    # Arachidonate-channel competitive block (ibuprofen + indomethacin,
    # additive in occupancy space because they compete for the same site).
    perox = max(1e-6, y[IX["PEROX"]])

    def cox_block(u_ibu, u_ind, u_apap, ki_ibu, ki_ind, ic50_apap):
        occ = u_ibu / ki_ibu + u_ind / ki_ind
        i_chan = occ / (1.0 + occ)
        # peroxidase-site reduction: acetaminophen competes with peroxide, so
        # its apparent IC50 scales with peroxide tone.
        ic50_eff = ic50_apap * (1.0 + perox / p["KPEROX"])
        i_perox = u_apap / (u_apap + ic50_eff)
        # the two sites are in series on one catalytic cycle -> multiplicative
        return 1.0 - (1.0 - i_chan) * (1.0 - i_perox), i_chan, i_perox

    d["ICOX_D"], d["ICHAN_D"], d["IPEROX_D"] = cox_block(
        d["UIBUE"], d["UINDE"], d["UAPAP"], p["KI_IBU"], p["KI_IND"], p["IC50_APAP"])
    d["ICOX_K"], _, _ = cox_block(
        d["UIBU"], d["UIND"], d["UAPAP"], p["KI_IBU_K"], p["KI_IND_K"], p["IC50_APAP_K"])
    d["ICOX_G"], _, _ = cox_block(
        d["UIBU"], d["UIND"], d["UAPAP"], p["KI_IBU_G"], p["KI_IND_G"], p["IC50_APAP_G"])
    occ_plt = d["UIBU"] / p["KI_IBU_PLT"] + d["UIND"] / p["KI_IND_PLT"]
    d["ICOX_PLT"] = occ_plt / (1.0 + occ_plt)

    # --------------------------------------------------- ductal relaxant drive
    ep4 = max(0.0, y[IX["EP4"]])
    pge2d = max(0.0, y[IX["PGE2D"]])
    drive = pge2d * ep4
    d["RELAX_PGE2"] = p["EMAX_PGE2"] * drive / (p["EC50_PGE2"] + drive)
    no = max(0.0, y[IX["NOD"]])
    d["RELAX_NO"] = p["EMAX_NO"] * no / (p["EC50_NO"] + no)

    # ------------------------------------------------------- O2-driven tone
    p50 = p["P50_O2"] * math.exp(p["P50_GA"] * (28.0 - GA))
    d["P50"] = p50
    pao2 = p["PAO2"]
    d["O2GAIN"] = pao2 ** 2 / (pao2 ** 2 + p50 ** 2)

    # ---------------- maximal achievable occlusion (muscular wall maturity) ---
    # This, not the drive, is where extreme prematurity puts its hard limit:
    # a thin immature medial layer cannot generate complete luminal occlusion.
    # At 24 wk the ceiling leaves a 0.44 mm residual lumen -- above the 0.30 mm
    # closure threshold -- so NO exposure of ANY drug closes that duct.  At
    # 26 wk the ceiling leaves 0.24 mm and closure becomes achievable.  The
    # steep clinical GA gradient in treatment success falls out of this rather
    # than being fitted to it.
    d["TMAX"] = 1.0 / (1.0 + math.exp(-(GA - p["TMAX_GA50"]) * p["TMAX_SL"]))

    # ------------------------- net constrictor drive -> tone (sigmoid) -------
    d["FINV"] = math.exp(p["KINVGA"] * (GA - 26.0))    # involution speed factor
    d["MAT"] = 1.0 - math.exp(-t * d["FINV"] / p["TAU_MAT"])
    d["NET"] = (p["GO2"] * d["O2GAIN"] + p["GMAT"] * d["MAT"]
                - p["GPGE2"] * d["RELAX_PGE2"] - p["GNO"] * d["RELAX_NO"])
    d["TONE_TGT"] = d["TMAX"] / (
        1.0 + math.exp(-(d["NET"] - p["NET50"]) / p["NETW"]))
    d["TONE_TGT"] = min(1.0, max(0.0, d["TONE_TGT"]))

    # ------------------------------------------------------- ductal geometry
    dmax = max(0.6, p["DMAX_INT"] + p["DMAX_GA"] * GA)
    d["DMAX"] = dmax
    tone = min(1.0, max(0.0, y[IX["TONE"]]))
    remod = min(1.0, max(0.0, y[IX["REMOD"]]))
    d["DDUCT"] = dmax * (1.0 - remod) * (1.0 - p["ALPHA_D"] * tone)
    d["CLOSED"] = 1.0 if d["DDUCT"] < p["DCLOSE"] else 0.0

    # ------------------------------------------------------------ haemodynamics
    dd = max(0.05, d["DDUCT"])
    d["RDUCT"] = p["KRDUCT"] / (dd ** 4) + p["RDUCT_MIN"]
    PVR = max(0.004, y[IX["PVR"]])
    SVR = max(0.05, y[IX["SVR"]])
    PLA = max(0.0, y[IX["PLA"]])
    qmax = p["QMAX_LV"] + p["QMAX_GA"] * (GA - 26.0)

    # Solve the 3-equation loop:
    #   Pao = Qs*SVR ; Ppa = Qp*PVR + Pla ; Qsh = (Pao-Ppa)/Rduct ; Qp = Qs+Qsh
    denom = d["RDUCT"] + PVR
    qs = p["QS_TARGET"]
    qsh = (qs * (SVR - PVR) - PLA) / denom
    qsh = max(0.0, qsh)
    qp = qs + qsh
    if qp > qmax:                       # LV reserve exhausted -> systemic steal
        qs = (qmax + PLA / denom) / (1.0 + (SVR - PVR) / denom)
        qs = max(40.0, qs)
        qsh = max(0.0, (qs * (SVR - PVR) - PLA) / denom)
        qp = qs + qsh
    d["QSH"] = qsh
    d["QP"] = qp
    d["QS"] = qs
    d["PAO"] = qs * SVR
    d["PPA"] = qp * PVR + PLA
    d["QPQS"] = qp / max(1.0, qs)
    d["SHFRAC"] = qsh / max(1.0, qp)
    d["QSDEF"] = max(0.0, (p["QS_TARGET"] - qs) / p["QS_TARGET"])

    # diastolic runoff and organ perfusion
    d["PDIA"] = max(8.0, p["PDIA0"] + 0.05 * (GA - 26.0) * 4.0
                    - p["KDIA"] * qsh)
    vc_ind = p["VC_IND"] * d["CIND"] / (d["CIND"] + p["VC50_IND"])
    vc_ibu = p["VC_IBU"] * d["CIBU"] / (d["CIBU"] + p["VC50_IBU"])
    d["VCDRUG"] = vc_ind + vc_ibu
    perf = (d["PDIA"] / p["PDIA0"]) * (1.0 - d["VCDRUG"])
    d["QCERREL"] = max(0.05, perf)
    d["QMESREL"] = max(0.05, perf * (1.0 - 0.25 * d["SHFRAC"]))
    d["QRENREL"] = max(0.05, perf)

    # ------------------------------------------------------------ lung mechanics
    evlw = max(0.0, y[IX["EVLW"]])
    crs = max(0.05, y[IX["CRS"]])
    crs_ref = p["CRS0"] + p["CRS_GA"] * (GA - 26.0)
    d["CRSREL"] = crs / crs_ref
    # ventilator severity index: worse compliance and more lung water -> more support
    d["VENTSEV"] = max(0.0, (1.0 / max(0.15, d["CRSREL"]) - 1.0)
                       + 0.030 * max(0.0, evlw - p["EVLW0"]))
    d["FIO2"] = min(1.0, 0.21 + 0.42 * d["VENTSEV"])

    # ------------------------------------------------------------------ kidney
    gfr = max(0.02, y[IX["GFR"]])
    d["UO"] = p["KUO"] * gfr * 60.0 / 60.0     # mL/kg/h
    d["SCR"] = y[IX["SCR"]]

    # --------------------------------------------------------- unbound APAP tox
    d["NAPQI"] = p["FNAPQI"] * p["CL_APAP0"] * d["CAPAP"]   # mg/kg/h flux
    return d


# ==============================================================================
#  SECTION 4.  ODE SYSTEM
# ==============================================================================

def rhs(t, y, p):
    a = derived(t, y, p)
    dy = np.zeros(NST)
    GA = p["GA"]
    PNAd = a["PNAd"]

    # ------------------------------------------------------------------ PK
    cl_ibu = p["CL_IBU0"] * (1.0 + p["CLMAT_IBU"] * PNAd)
    cl_ind = p["CL_IND0"] * (1.0 + p["CLMAT_IND"] * PNAd)
    cl_apap = p["CL_APAP0"] * (1.0 + p["CLMAT_APAP"] * PNAd)

    c1i = y[IX["IBU1"]] / p["V1_IBU"]
    c2i = y[IX["IBU2"]] / p["V2_IBU"]
    dy[IX["IBU1"]] = p["RIN_IBU"] - cl_ibu * c1i - p["Q_IBU"] * (c1i - c2i)
    dy[IX["IBU2"]] = p["Q_IBU"] * (c1i - c2i)
    dy[IX["IBUE"]] = p["KE0_IBU"] * (a["CIBU"] - y[IX["IBUE"]])

    c1n = y[IX["IND1"]] / p["V1_IND"]
    c2n = y[IX["IND2"]] / p["V2_IND"]
    dy[IX["IND1"]] = p["RIN_IND"] - cl_ind * c1n - p["Q_IND"] * (c1n - c2n)
    dy[IX["IND2"]] = p["Q_IND"] * (c1n - c2n)
    dy[IX["INDE"]] = p["KE0_IND"] * (a["CIND"] - y[IX["INDE"]])

    c1a = y[IX["APAP1"]] / p["V1_APAP"]
    c2a = y[IX["APAP2"]] / p["V2_APAP"]
    dy[IX["APAPG"]] = -p["KA_APAP"] * y[IX["APAPG"]]
    dy[IX["APAP1"]] = (p["F_APAP"] * p["KA_APAP"] * y[IX["APAPG"]]
                       - cl_apap * c1a - p["Q_APAP"] * (c1a - c2a))
    dy[IX["APAP2"]] = p["Q_APAP"] * (c1a - c2a)

    # ---------------------------------------------------- hepatic safety
    gsh = max(0.0, y[IX["GSH"]])
    napqi = a["NAPQI"]
    quench = napqi * gsh
    dy[IX["GSH"]] = p["KGSH"] * (p["GSHCAP"] - gsh) - 0.55 * quench
    unquenched = napqi * (1.0 - gsh)
    dy[IX["ALT"]] = (p["KALT"] * unquenched * 100.0
                     - p["KALTOUT"] * (y[IX["ALT"]] - p["ALT0"]))

    # ---------------------------------------------------- prostanoid system
    # circulating PGE2: placental source is gone at t=0 (initial condition),
    # clearance grows with pulmonary blood flow x 15-PGDH capacity.
    pgdh = max(0.0, y[IX["PGDH"]])
    qpfrac = a["QP"] / 200.0
    dy[IX["PGE2P"]] = (p["RPROD_PGE2P"] * (1.0 - a["ICOX_G"])
                       - (p["KPUL_PGE2"] * pgdh * qpfrac + p["KOTH_PGE2"])
                       * y[IX["PGE2P"]])
    dy[IX["PGDH"]] = p["KPGDH"] * (1.0 - pgdh)

    # ductal tissue PGE2.  FSYN is the postnatal involution of local ductal
    # prostanoid synthesis capacity -- faster at higher GA.  This is the term
    # that makes spontaneous closure happen, and it is the reason the treatment
    # effect on any fixed-time closure endpoint has a ceiling.
    tau_syn = p["TAUSYN0"] / math.exp(p["KTAUSYN"] * (GA - 26.0))
    fsyn = p["FSYN_INF"] + (1.0 - p["FSYN_INF"]) * math.exp(-t / tau_syn)
    dy[IX["PGE2D"]] = (p["KSYN_PGE2D"] * fsyn * (1.0 - a["ICOX_D"])
                       - p["KDEG_PGE2D"] * y[IX["PGE2D"]]
                       + p["KIN_PGE2D"] * (y[IX["PGE2P"]] / p["PGE2P0"])
                       * (1.0 - a["ICOX_D"]))

    # cAMP and Ca are fast readouts of the PGE2/EP4 axis (kept as states so the
    # map's signalling cluster corresponds to real variables)
    dy[IX["CAMP"]] = p["KCAMP"] * (a["RELAX_PGE2"] / max(1e-6, p["EMAX_PGE2"])
                                   - y[IX["CAMP"]])
    dy[IX["CASMC"]] = 2.0 * (y[IX["TONE"]] - y[IX["CASMC"]])

    # NO tone: lower at higher GA, and lost postnatally.  This is the "second
    # dilator" that makes COX inhibition alone insufficient in the immature
    # duct and motivates the (so far unsuccessful) NOS-inhibition idea.
    finv = a["FINV"]
    no_tgt = (math.exp(-p["NOGA"] * (GA - 24.0))
              * math.exp(-p["KNODEC"] * finv * t))
    dy[IX["NOD"]] = p["KNO"] * (no_tgt - y[IX["NOD"]])

    # EP4 down-regulates postnatally toward EP4INF of its gestation-specific
    # starting density
    ep4_tgt = math.exp(-p["EP4GA"] * (GA - 24.0))
    dy[IX["EP4"]] = -p["KEP4"] * finv * (y[IX["EP4"]] - p["EP4INF"] * ep4_tgt)

    # peroxide tone (raises acetaminophen's apparent IC50)
    perox_tgt = 1.0 + (p["PEROXSEP"] - 1.0) * p["SEPSIS"]
    dy[IX["PEROX"]] = p["KPEROXIN"] * (perox_tgt - y[IX["PEROX"]])

    # ---------------------------------------------------- ductal contraction
    dy[IX["TONE"]] = p["KTONE"] * (a["TONE_TGT"] - y[IX["TONE"]])

    # ---------------------------------------------------- wall O2 and remodelling
    wall = max(0.05, y[IX["WALLTH"]])
    wall_ref = max(0.05, p["WALLTH_INT"] + p["WALLTH_GA"] * GA)
    dy[IX["WALLTH"]] = p["KWALL"] * (1.0 - wall)
    lumen_frac = a["DDUCT"] / a["DMAX"]
    # Two supplies with OPPOSITE wall-thickness dependence (see parameter
    # comments): diffusion, adequate only across a thin wall, and vasa
    # vasorum, present only in a thick wall and shut off by constriction.
    wr = wall_ref * wall
    diff_supply = p["FLUM_O2"] * (0.35 + 0.65 * lumen_frac) \
        / (1.0 + (p["KDIFF"] * wr) ** 2)
    vasa_supply = p["FVASA_O2"] * wr * (1.0 - min(1.0, y[IX["TONE"]]))
    wallo2 = p["PAO2"] * (diff_supply + vasa_supply)
    hyp = max(0.0, p["O2CRIT"] - wallo2)
    dy[IX["VEGFD"]] = p["KVEGF"] * hyp - p["KVEGFOUT"] * y[IX["VEGFD"]]
    dy[IX["TGFBD"]] = (p["KTGFB"] * y[IX["VEGFD"]]
                       - p["KTGFBOUT"] * y[IX["TGFBD"]])
    remod = min(1.0, max(0.0, y[IX["REMOD"]]))
    frem = 1.0 / (1.0 + math.exp(-(GA - p["KREMGA50"]) * p["KREMGA_SL"]))
    dy[IX["REMOD"]] = (p["KREMOD"] * frem * y[IX["TGFBD"]] * (1.0 - remod)
                       - p["KUNREMOD"] * remod)

    # ------------------------------------------------------------ haemodynamics
    pvr_tgt = (p["PVRINF"] + p["PVR_EVLW"] * max(0.0, y[IX["EVLW"]] - p["EVLW0"]))
    dy[IX["PVR"]] = -p["KPVR"] * (y[IX["PVR"]] - pvr_tgt)
    dy[IX["SVR"]] = p["KSVR"] * (p["SVRINF"] - y[IX["SVR"]])
    dy[IX["PLA"]] = p["KPLA"] * (p["PLA0"] + p["PLA_Q"] * a["QSH"] - y[IX["PLA"]])
    dy[IX["LVDIL"]] = (p["KLVDIL"] * a["SHFRAC"] * 10.0
                       - p["KLVDILOUT"] * y[IX["LVDIL"]])

    # ------------------------------------------------------- lung water/mechanics
    excess_qp = max(0.0, a["QP"] - 200.0)
    dy[IX["EVLW"]] = (p["KEVLW"] * excess_qp * (1.0 + 0.5 * y[IX["PLA"]] / 4.0)
                      - p["KEVLWOUT"] * (y[IX["EVLW"]] - p["EVLW0"]))
    crs_ref = p["CRS0"] + p["CRS_GA"] * (GA - 26.0)
    crs_tgt = crs_ref * (1.0 - min(0.80, p["CRS_EVLW"]
                                   * max(0.0, y[IX["EVLW"]] - p["EVLW0"])))
    crs_tgt *= (1.0 + 0.10 * p["ANTESTER"])
    dy[IX["CRS"]] = p["KCRS"] * (crs_tgt - y[IX["CRS"]])

    # ---------------------------------------------------------------- burden
    # PDA burden: time integral of shunt ABOVE the haemodynamic-significance
    # threshold, in (mL/min/kg)*day.  This -- not closure at day 7 -- is the
    # variable the outcome hazards see.
    dy[IX["PDABUR"]] = max(0.0, a["QSH"] - p["QSIG"]) / 24.0
    dy[IX["MESDEF"]] = max(0.0, 1.0 - a["QMESREL"]) / 24.0
    dy[IX["CERDEF"]] = max(0.0, 1.0 - a["QCERREL"]) / 24.0

    # ------------------------------------------------------------------ kidney
    dy[IX["PGE2K"]] = p["KPGE2K"] * ((1.0 - a["ICOX_K"]) - y[IX["PGE2K"]])
    gfr_ref = (p["GFR0"] + p["GFR_GA"] * (GA - 26.0)) * (1.0 + p["KGFRMAT"] * t)
    pge2k = max(0.0, y[IX["PGE2K"]])
    gfr_tgt = gfr_ref * ((1.0 - p["GFR_PGE2"]) + p["GFR_PGE2"] * pge2k) \
        * (0.55 + 0.45 * a["PDIA"] / p["PDIA0"]) * (1.0 - 0.5 * a["VCDRUG"])
    dy[IX["GFR"]] = 0.25 * (gfr_tgt - y[IX["GFR"]])
    # creatinine: production constant, elimination proportional to GFR
    dy[IX["SCR"]] = p["CRPROD"] - p["KCREL"] * y[IX["GFR"]] * y[IX["SCR"]]
    dy[IX["UOCUM"]] = a["UO"]
    dy[IX["FLUID"]] = p["FLUIDIN"] - a["UO"] - 2.2

    # -------------------------------------------------------- gut, platelet, bili
    dy[IX["PGE2G"]] = p["KPGE2G"] * ((1.0 - a["ICOX_G"]) - y[IX["PGE2G"]])
    dy[IX["TXA2"]] = p["KTXA2"] * ((1.0 - a["ICOX_PLT"]) - y[IX["TXA2"]])
    bt_tgt = 1.0 + 1.9 * (1.0 - max(0.0, y[IX["TXA2"]]))
    dy[IX["BTIME"]] = p["KBTIME"] * (bt_tgt - y[IX["BTIME"]])
    dy[IX["TBILI"]] = p["KBILIIN"] - p["KBILIOUT"] * y[IX["TBILI"]]
    bfree_tgt = (y[IX["TBILI"]] / p["TBILI0"]) * (
        1.0 + p["BILIDISP_IBU"] * a["CIBU"])
    dy[IX["BFREE"]] = 0.5 * (bfree_tgt - y[IX["BFREE"]])

    # ------------------------------------------------------ lung development
    arrest = (p["ALV_VENT"] * a["VENTSEV"]
              + p["ALV_EVLW"] * max(0.0, y[IX["EVLW"]] - p["EVLW0"]))
    dy[IX["ALV"]] = p["KALV"] * max(0.0, 1.0 - arrest) - 0.0
    dy[IX["VENTIX"]] = a["VENTSEV"] / 24.0

    # ------------------------------------------------------------ hazards
    # systemic COX-inhibition exposure, in fraction*days
    dy[IX["COXIX"]] = a["ICOX_G"] / 24.0

    ga_pen = math.exp(p["B_GA"] * max(0.0, 28.0 - GA))
    h_bpd = (p["H_BPD0"] * ga_pen
             * (1.0 + p["B_BUR"] * y[IX["PDABUR"]])
             * (1.0 + p["B_VENT"] * y[IX["VENTIX"]])
             * (1.0 + p["B_BPD_COX"] * y[IX["COXIX"]])
             * (1.0 - 0.12 * p["ANTESTER"]))
    dy[IX["HBPD"]] = h_bpd

    dy[IX["HNEC"]] = (p["H_NEC0"] * math.exp(0.30 * max(0.0, 28.0 - GA))
                      * (1.0 + p["B_NEC_MES"] * y[IX["MESDEF"]])
                      * (1.0 + p["B_NEC_PGE2G"] * (1.0 - max(0.0, y[IX["PGE2G"]]))))

    # severe IVH: a first-96-h phenomenon.  Two opposing drug terms --
    # platelet COX-1 inhibition prolongs bleeding time (harm) while
    # indomethacin matures germinal-matrix vessels (benefit, TIPP).
    ivh_gate = 1.0 if t < p["IVH_WINDOW"] else math.exp(-(t - p["IVH_WINDOW"]) / 48.0)
    ind_prot = 1.0 - p["IVH_IND_PROT"] * (y[IX["INDE"]] / (y[IX["INDE"]] + 0.25))
    dy[IX["HIVH"]] = (p["H_IVH0"] * ivh_gate
                      * math.exp(p["B_IVH_GA"] * max(0.0, 28.0 - GA))
                      * (1.0 + p["B_IVH_BT"] * max(0.0, y[IX["BTIME"]] - 1.0))
                      * (1.0 + p["B_IVH_CER"] * a["QSDEF"])
                      * max(0.15, ind_prot))

    dy[IX["HSIP"]] = (p["H_SIP0"] * math.exp(0.30 * max(0.0, 28.0 - GA))
                      * (1.0 + p["B_SIP_IND"] * (y[IX["INDE"]] / (y[IX["INDE"]] + 0.3)))
                      * (1.0 + p["B_SIP_HC"] * p["HCORT"]))

    dy[IX["HDEATH"]] = (p["H_DTH0"] * math.exp(p["B_DTH_GA"] * max(0.0, 28.0 - GA))
                        * (1.0 + p["B_DTH_QS"] * a["QSDEF"])
                        * (1.0 + 0.6 * a["VENTSEV"]))
    return dy


# ==============================================================================
#  SECTION 5.  INITIAL CONDITIONS
# ==============================================================================

def init_state(p):
    y = np.zeros(NST)
    GA = p["GA"]
    y[IX["GSH"]] = 1.0
    y[IX["ALT"]] = p["ALT0"]
    y[IX["PGE2D"]] = 1.0
    y[IX["PGE2P"]] = p["PGE2P0"]
    y[IX["PGDH"]] = max(0.02, min(1.0, p["PGDH_SL"] * (GA - p["PGDH0_GA"])))
    y[IX["CAMP"]] = 1.0
    y[IX["NOD"]] = math.exp(-p["NOGA"] * (GA - 24.0))
    y[IX["EP4"]] = math.exp(-p["EP4GA"] * (GA - 24.0))
    y[IX["PEROX"]] = 1.0 + (p["PEROXSEP"] - 1.0) * p["SEPSIS"]
    y[IX["CASMC"]] = 0.05
    y[IX["TONE"]] = 0.02
    y[IX["REMOD"]] = 0.0
    y[IX["WALLTH"]] = 0.35
    y[IX["VEGFD"]] = 0.0
    y[IX["TGFBD"]] = 0.0
    y[IX["PVR"]] = p["PVR0"]
    y[IX["SVR"]] = p["SVR0"]
    y[IX["PLA"]] = p["PLA0"]
    y[IX["LVDIL"]] = 0.0
    y[IX["EVLW"]] = p["EVLW0"]
    y[IX["CRS"]] = (p["CRS0"] + p["CRS_GA"] * (GA - 26.0)) * (1.0 + 0.10 * p["ANTESTER"])
    y[IX["GFR"]] = p["GFR0"] + p["GFR_GA"] * (GA - 26.0)
    y[IX["SCR"]] = p["SCR0"]
    y[IX["PGE2K"]] = 1.0
    y[IX["PGE2G"]] = 1.0
    y[IX["TXA2"]] = 1.0
    y[IX["BTIME"]] = 1.0
    y[IX["TBILI"]] = p["TBILI0"]
    y[IX["BFREE"]] = 1.0
    y[IX["ALV"]] = p["ALV0"]
    return y


# ==============================================================================
#  SECTION 6.  SIMULATION DRIVER (bolus events by interval splitting)
# ==============================================================================

def simulate(p, doses, tend_h, dt=1.0):
    """doses: list of (time_h, compartment_name, amount_mg_per_kg) bolus events,
    or (time_h, 'RIN_IBU', rate) to change an infusion rate."""
    y = init_state(p)
    p = dict(p)
    events = sorted(doses, key=lambda e: e[0])
    tgrid = np.arange(0.0, tend_h + dt / 2, dt)
    out = {n: np.zeros(len(tgrid)) for n in SNAMES}
    extra_names = ["DDUCT", "QSH", "QP", "QS", "QPQS", "PAO", "PPA", "PDIA",
                   "CIBU", "CIND", "CAPAP", "UIBU", "UIND", "UAPAP",
                   "ICOX_D", "ICOX_K", "ICOX_G", "ICOX_PLT", "ICHAN_D",
                   "IPEROX_D", "TONE_TGT", "RELAX_PGE2", "RELAX_NO", "O2GAIN",
                   "TMAX", "CLOSED", "VENTSEV", "FIO2", "UO", "QCERREL",
                   "QMESREL", "QRENREL", "SHFRAC", "QSDEF", "CRSREL", "VCDRUG"]
    for n in extra_names:
        out[n] = np.zeros(len(tgrid))

    # breakpoints = dose times inside the horizon
    bps = sorted(set([0.0] + [e[0] for e in events if 0 < e[0] < tend_h] + [tend_h]))
    ti = 0
    ycur = y.copy()
    # apply any t=0 events
    for e in events:
        if abs(e[0]) < 1e-9:
            if e[1] in IX:
                ycur[IX[e[1]]] += e[2]
            else:
                p[e[1]] = e[2]

    for k in range(len(bps) - 1):
        t0, t1 = bps[k], bps[k + 1]
        seg = tgrid[(tgrid >= t0 - 1e-9) & (tgrid <= t1 + 1e-9)]
        if len(seg) == 0 or seg[0] > t0 + 1e-9:
            seg = np.concatenate(([t0], seg))
        sol = solve_ivp(rhs, (t0, t1), ycur, t_eval=seg, args=(p,),
                        method="LSODA", rtol=1e-6, atol=1e-9, max_step=2.0)
        if not sol.success:
            raise RuntimeError("integration failed: " + sol.message)
        for j, tt in enumerate(sol.t):
            idx = int(round(tt / dt))
            if 0 <= idx < len(tgrid) and abs(tgrid[idx] - tt) < dt / 2:
                yy = sol.y[:, j]
                for n in SNAMES:
                    out[n][idx] = yy[IX[n]]
                a = derived(tt, yy, p)
                for n in extra_names:
                    out[n][idx] = a[n]
        ycur = sol.y[:, -1].copy()
        for e in events:
            if abs(e[0] - t1) < 1e-9:
                if e[1] in IX:
                    ycur[IX[e[1]]] += e[2]
                else:
                    p[e[1]] = e[2]
    out["time"] = tgrid
    out["day"] = tgrid / 24.0
    # event probabilities from cumulative hazards
    for h, nm in [("HBPD", "P_BPD"), ("HNEC", "P_NEC"), ("HIVH", "P_IVH"),
                  ("HSIP", "P_SIP"), ("HDEATH", "P_DEATH")]:
        out[nm] = 1.0 - np.exp(-out[h])
    return out


# ==============================================================================
#  SECTION 7.  DOSING REGIMENS
# ==============================================================================

def rx_ibuprofen(start_h, hi=False):
    """IV ibuprofen lysine 10-5-5 mg/kg q24h (standard) or 20-10-10 (high)."""
    d0, d1 = (20.0, 10.0) if hi else (10.0, 5.0)
    return [(start_h, "IBU1", d0),
            (start_h + 24.0, "IBU1", d1),
            (start_h + 48.0, "IBU1", d1)]


def rx_ibuprofen_infusion(start_h):
    """Continuous infusion: 20 mg/kg over 24 h then 10 mg/kg/24 h x 2."""
    return [(start_h, "RIN_IBU", 20.0 / 24.0),
            (start_h + 24.0, "RIN_IBU", 10.0 / 24.0),
            (start_h + 72.0, "RIN_IBU", 0.0)]


def rx_indomethacin(start_h, dose=0.2, maint=0.1, n=3, interval=24.0):
    ev = [(start_h, "IND1", dose)]
    for i in range(1, n):
        ev.append((start_h + i * interval, "IND1", maint))
    return ev


def rx_acetaminophen(start_h, dose=15.0, interval=6.0, days=3, iv=True):
    ev = []
    n = int(days * 24 / interval)
    cmt = "APAP1" if iv else "APAPG"
    for i in range(n):
        ev.append((start_h + i * interval, cmt, dose))
    return ev


# ==============================================================================
#  SECTION 8.  READOUT HELPERS
# ==============================================================================

def closure_day(out, thresh_mm=None, p=None):
    """First day at which ductal diameter drops below the closure threshold."""
    th = thresh_mm if thresh_mm is not None else (p or default_params())["DCLOSE"]
    dd = out["DDUCT"]
    idx = np.where(dd < th)[0]
    return out["day"][idx[0]] if len(idx) else None


def sustained_closure(out, p, at_day, need_days=2.0):
    """Closed at `at_day` AND still closed need_days later (guards against
    transient constriction being scored as closure)."""
    dd = out["DDUCT"]
    day = out["day"]
    i0 = int(np.argmin(np.abs(day - at_day)))
    i1 = int(np.argmin(np.abs(day - (at_day + need_days))))
    return bool(dd[i0] < p["DCLOSE"] and dd[min(i1, len(dd) - 1)] < p["DCLOSE"])


def reopened(out, p):
    """Reached closure at some point, then diameter recovered above threshold."""
    dd = out["DDUCT"]
    idx = np.where(dd < p["DCLOSE"])[0]
    if not len(idx):
        return False
    return bool(np.any(dd[idx[0]:] > p["DCLOSE"] * 1.6))


def at_day(out, name, d):
    i = int(np.argmin(np.abs(out["day"] - d)))
    return out[name][i]


def auc(out, name, d0, d1):
    day = out["day"]
    m = (day >= d0) & (day <= d1)
    return float(np.trapezoid(out[name][m], day[m]))


# ==============================================================================
#  SECTION 9.  SCENARIOS
# ==============================================================================

def build_scenarios():
    """16 scenarios.  Each returns (label, params, dose events, horizon_days)."""
    S = []

    def base(**kw):
        p = default_params()
        p.update(kw)
        return p

    S.append(("S1  Expectant management, 26 wk",
              base(GA=26.0), [], 90))
    S.append(("S2  Early ibuprofen d2 (Baby-OSCAR), 26 wk",
              base(GA=26.0), rx_ibuprofen(48.0), 90))
    S.append(("S3  Late ibuprofen d7, standard dose, 26 wk",
              base(GA=26.0), rx_ibuprofen(168.0), 90))
    S.append(("S4  Late ibuprofen d7, HIGH dose 20-10-10, 26 wk",
              base(GA=26.0), rx_ibuprofen(168.0, hi=True), 90))
    S.append(("S5  Continuous ibuprofen infusion d2, 26 wk",
              base(GA=26.0), rx_ibuprofen_infusion(48.0), 90))
    S.append(("S6  Indomethacin 0.2-0.1-0.1 d2, 26 wk",
              base(GA=26.0), rx_indomethacin(48.0), 90))
    S.append(("S7  Prophylactic indomethacin 0.1 q24 x3 from 8 h (TIPP)",
              base(GA=26.0), rx_indomethacin(8.0, dose=0.1, maint=0.1, n=3), 90))
    S.append(("S8  IV acetaminophen 15 q6h x3 d, d2, 26 wk",
              base(GA=26.0), rx_acetaminophen(48.0), 90))
    S.append(("S9  Acetaminophen d2 WITH chorioamnionitis (high peroxide)",
              base(GA=26.0, SEPSIS=1.0), rx_acetaminophen(48.0), 90))
    S.append(("S10 Ibuprofen d2 WITH chorioamnionitis (control for S9)",
              base(GA=26.0, SEPSIS=1.0), rx_ibuprofen(48.0), 90))
    S.append(("S11 Ibuprofen + acetaminophen combination d2",
              base(GA=26.0), rx_ibuprofen(48.0) + rx_acetaminophen(48.0), 90))
    S.append(("S12 Ibuprofen d2, second course d7 after reopening",
              base(GA=26.0), rx_ibuprofen(48.0) + rx_ibuprofen(168.0), 90))
    S.append(("S13 Indomethacin d2 + early hydrocortisone (SIP risk)",
              base(GA=26.0, HCORT=1.0), rx_indomethacin(48.0), 90))
    S.append(("S14 Ibuprofen d2, extreme prematurity 24 wk",
              base(GA=24.0, BW=0.62), rx_ibuprofen(48.0), 90))
    S.append(("S15 Ibuprofen d2, 29 wk",
              base(GA=29.0, BW=1.20), rx_ibuprofen(48.0), 90))
    S.append(("S16 Targeted: treat at d10 only (echo-guided), 26 wk",
              base(GA=26.0), rx_ibuprofen(240.0, hi=True), 90))
    return S


# ==============================================================================
#  SECTION 10.  MAIN
# ==============================================================================

def hline(c="="):
    print(c * 78)


def run_all(quick=False):
    p0 = default_params()

    hline()
    print(" PDA-OF-PREMATURITY QSP MODEL -- PYTHON REFERENCE IMPLEMENTATION")
    print(" %d ODE states" % NST)
    hline()

    # ---------------------------------------------------------------- PK check
    print("\n[1] DRUG EXPOSURE CHECKS (26 wk, 0.8 kg)\n")
    o = simulate(p0, rx_ibuprofen(48.0), 24 * 12, dt=0.25)
    cmax = o["CIBU"].max()
    i48 = int(48 / 0.25)
    # terminal half-life after last dose
    seg = o["CIBU"][int(96 / 0.25):int(240 / 0.25)]
    tt = o["time"][int(96 / 0.25):int(240 / 0.25)]
    m = seg > 0.01
    kel = -np.polyfit(tt[m], np.log(seg[m]), 1)[0]
    vss = p0["V1_IBU"] + p0["V2_IBU"]
    t12_instant = math.log(2) * vss / p0["CL_IBU0"]
    print("  Ibuprofen 10-5-5 mg/kg IV q24h from 48 h")
    print("    Cmax(total)            %8.1f mg/L" % cmax)
    print("    C trough @ 24 h        %8.1f mg/L" % o["CIBU"][int(72 / 0.25)])
    print("    unbound Cmax           %8.2f uM   (fu = %.3f)"
          % (o["UIBU"].max(), p0["FU_IBU"]))
    print("    t1/2 at PNA 0          %8.1f h    (reported range 20-43 h)"
          % t12_instant)
    print("    APPARENT terminal t1/2 %8.1f h" % (math.log(2) / kel))
    print("      The apparent value is SHORTER than the instantaneous one, and")
    print("      that is not a bug: clearance matures ~20%%/day, so a log-linear")
    print("      slope fitted across days 4-10 is steepened by the maturation")
    print("      occurring inside the sampling window.  Published preterm")
    print("      ibuprofen half-lives estimated without a maturation term are")
    print("      biased the same way, which is one reason the reported ranges")
    print("      are so wide and mutually inconsistent.")
    print("    AUC0-inf               %8.0f mg*h/L" % auc(o, "CIBU", 2, 12))
    print("    peak ductal COX inhib  %8.1f %%" % (100 * o["ICOX_D"].max()))

    o = simulate(p0, rx_indomethacin(48.0), 24 * 12, dt=0.25)
    print("\n  Indomethacin 0.2-0.1-0.1 mg/kg IV q24h from 48 h")
    print("    Cmax(total)            %8.3f mg/L" % o["CIND"].max())
    print("    unbound Cmax           %8.4f uM" % o["UIND"].max())
    print("    peak ductal COX inhib  %8.1f %%" % (100 * o["ICOX_D"].max()))

    o = simulate(p0, rx_acetaminophen(48.0), 24 * 12, dt=0.25)
    print("\n  Acetaminophen 15 mg/kg IV q6h x 72 h from 48 h")
    print("    Cmax(total)            %8.1f mg/L" % o["CAPAP"].max())
    print("    C trough (steady)      %8.1f mg/L" % o["CAPAP"][int(116 / 0.25)])
    print("    unbound Cmax           %8.0f uM" % o["UAPAP"].max())
    print("    peak ductal COX inhib  %8.1f %%" % (100 * o["ICOX_D"].max()))
    print("    peak ALT               %8.0f U/L" % o["ALT"].max())

    # ------------------------------------------- clearance maturation (thesis 4)
    print("\n[2] CLEARANCE MATURATION: THE SAME REGIMEN AT d2 VS d7\n")
    a = simulate(p0, rx_ibuprofen(48.0), 24 * 12, dt=0.25)
    b = simulate(p0, rx_ibuprofen(168.0), 24 * 16, dt=0.25)
    auc_a = auc(a, "CIBU", 2.0, 8.0)
    auc_b = auc(b, "CIBU", 7.0, 13.0)
    print("  AUC over the 6 d following the first dose")
    print("    start d2   %8.0f mg*h/L   peak COX inhibition %5.1f %%"
          % (auc_a, 100 * a["ICOX_D"].max()))
    print("    start d7   %8.0f mg*h/L   peak COX inhibition %5.1f %%"
          % (auc_b, 100 * b["ICOX_D"].max()))
    print("    ratio d7/d2                %6.2f" % (auc_b / auc_a))
    c = simulate(p0, rx_ibuprofen(168.0, hi=True), 24 * 16, dt=0.25)
    print("    start d7, HIGH dose %6.0f mg*h/L  (restores exposure: ratio %.2f)"
          % (auc(c, "CIBU", 7.0, 13.0), auc(c, "CIBU", 7.0, 13.0) / auc_a))

    # --------------------------------------------------- spontaneous closure
    print("\n[3] SPONTANEOUS CLOSURE WITHOUT ANY DRUG (calibration target:")
    print("    Semberova 2017 -- median 71 d at <26 wk, faster with rising GA)\n")
    print("    GA(wk)  closure day   d(day3)mm  Qsh(day3)  Qp:Qs(d3)")
    for ga in [24, 25, 26, 27, 28, 30, 32]:
        p = default_params()
        p["GA"] = float(ga)
        p["BW"] = 0.5 + 0.11 * (ga - 23)
        o = simulate(p, [], 24 * 130, dt=2.0)
        cd = closure_day(o, p=p)
        print("     %4d   %10s   %8.2f   %8.0f   %7.2f"
              % (ga, ("%.0f" % cd) if cd else ">130",
                 at_day(o, "DDUCT", 3), at_day(o, "QSH", 3),
                 at_day(o, "QPQS", 3)))

    # ------------------------------------------------------------- scenarios
    # Reporting rule that matters: a scenario's "closure" is only counted as
    # TREATMENT-ATTRIBUTABLE if it happens within 8 days of the first dose.
    # Otherwise the table would credit every drug with the spontaneous closure
    # that was going to happen anyway -- which is precisely the reasoning error
    # this whole model exists to expose.  Drug effects on continuous readouts
    # are likewise reported as differences from a GA-matched no-drug control.
    print("\n[4] SIXTEEN TREATMENT SCENARIOS")
    print("    'close' = treatment-attributable closure (within 8 d of dose 1);")
    print("    'spont' = the same patient's closure day with NO drug;")
    print("    dSCr / dBUR / dCOMP are differences from that no-drug control.\n")
    hdr = ("  %-42s %6s %6s %5s %7s %6s %7s" %
           ("scenario", "close", "spont", "reop", "dBUR", "dSCr", "dCOMP"))
    print(hdr)
    print("  " + "-" * len(hdr))
    results = {}
    ctrl_cache = {}
    for label, p, doses, days in build_scenarios():
        o = simulate(p, doses, 24 * days, dt=1.0)
        results[label] = (p, o)
        key = (p["GA"], p["SEPSIS"], p["HCORT"], p["ANTESTER"], p["PAO2"])
        if key not in ctrl_cache:
            ctrl_cache[key] = simulate(p, [], 24 * days, dt=1.0)
        c = ctrl_cache[key]
        cd = closure_day(o, p=p)
        cd0 = closure_day(c, p=p)
        t0 = min([e[0] for e in doses], default=None)
        attrib = (cd is not None and t0 is not None
                  and cd <= t0 / 24.0 + 8.0)
        # dSCr must be the ELEMENTWISE peak difference from the control, not the
        # difference of two window maxima: neonatal creatinine falls
        # monotonically from the maternal value, so both maxima are simply the
        # birth value and the difference of maxima is identically zero however
        # nephrotoxic the drug is.  The drug's signature is a BLUNTED fall.
        w = slice(0, int(24 * 12))
        dscr = float(np.max(o["SCR"][w] - c["SCR"][w]))
        print("  %-42s %6s %6s %5s %+7.0f %+6.2f %+7.1f"
              % (label[:42],
                 ("d%.1f" % cd) if attrib else "no",
                 ("d%.0f" % cd0) if cd0 else ">end",
                 "YES" if reopened(o, p) else "-",
                 at_day(o, "PDABUR", 36) - at_day(c, "PDABUR", 36),
                 dscr,
                 100 * (at_day(o, "P_BPD", 70) - at_day(c, "P_BPD", 70))))

    # ------------------------------------------------- drug-class separation
    print("\n[5] THE TWO NSAIDs DIVERGE IN THE ORGANS, NOT AT THE DUCT")
    print("    (indomethacin has non-COX vasoconstriction; ibuprofen does not)\n")
    for lbl, key in [("ibuprofen d2", "S2  Early ibuprofen d2 (Baby-OSCAR), 26 wk"),
                     ("indomethacin d2", "S6  Indomethacin 0.2-0.1-0.1 d2, 26 wk"),
                     ("acetaminophen d2", "S8  IV acetaminophen 15 q6h x3 d, d2, 26 wk")]:
        p, o = results[key]
        w = slice(48, 48 + 96)
        print("    %-16s cerebral flow nadir %5.0f %%   mesenteric nadir %5.0f %%"
              % (lbl, 100 * o["QCERREL"][w].min(), 100 * o["QMESREL"][w].min()))
        print("    %-16s GFR nadir %5.0f %% of pre-dose   peak SCr %4.2f mg/dL"
              "   urine nadir %4.1f mL/kg/h"
              % ("", 100 * o["GFR"][w].min() / o["GFR"][47],
                 o["SCR"][w].max(), o["UO"][w].min()))
        print("    %-16s bleeding-time index peak %4.2f   platelet COX %5.1f %%"
              % ("", o["BTIME"][w].max(), 100 * o["ICOX_PLT"][w].max()))

    # ---------------------------------------------- peroxide thesis (thesis 3)
    print("\n[6] PREDICTION: ACETAMINOPHEN IS PEROXIDE-SENSITIVE, IBUPROFEN IS NOT\n")
    for lbl, k in [("APAP, no sepsis", "S8  IV acetaminophen 15 q6h x3 d, d2, 26 wk"),
                   ("APAP, chorioamnionitis", "S9  Acetaminophen d2 WITH chorioamnionitis (high peroxide)"),
                   ("IBU, no sepsis", "S2  Early ibuprofen d2 (Baby-OSCAR), 26 wk"),
                   ("IBU, chorioamnionitis", "S10 Ibuprofen d2 WITH chorioamnionitis (control for S9)")]:
        p, o = results[k]
        w = slice(48, 48 + 120)
        # Attributable closure only: the drug gets no credit for the
        # spontaneous closure that was going to happen at day 17 regardless.
        cd = closure_day(o, p=p)
        att = cd is not None and cd <= 10.0
        print("    %-24s peak ductal COX inhibition %5.1f %%   tone peak %4.2f"
              "   closed %s"
              % (lbl, 100 * o["ICOX_D"][w].max(), o["TONE"][w].max(),
                 ("yes, d%.1f" % cd) if att else "NO"))

    # --------------------------------------------------- bistability (thesis 4)
    print("\n[7] WHY PRETERM DUCTS REOPEN: WALL O2 GATES REMODELLING\n")
    print("    GA(wk)  wall idx  wallPO2@peak-tone  neointima d20  reopened")
    for ga in [24, 26, 28, 30, 34, 38]:
        p = default_params()
        p["GA"] = float(ga)
        p["BW"] = 0.5 + 0.11 * (ga - 23)
        o = simulate(p, rx_ibuprofen(48.0), 24 * 40, dt=1.0)
        wall_ref = p["WALLTH_INT"] + p["WALLTH_GA"] * ga
        i = int(np.argmax(o["TONE"]))
        lf = o["DDUCT"][i] / (p["DMAX_INT"] + p["DMAX_GA"] * ga)
        w = o["WALLTH"][i]
        po2 = p["PAO2"] * (p["FLUM_O2"] * lf + p["FVASA_O2"]) / (1 + 2 * wall_ref * w)
        print("     %4d   %8.2f   %17.1f   %13.2f   %8s"
              % (ga, wall_ref, po2, at_day(o, "REMOD", 20),
                 "yes" if reopened(o, p) else "no"))

    if quick:
        print("\n(--quick: population runs skipped)")
        return results

    # ------------------------------------------------- population / trial recon
    print("\n[8] VIRTUAL TRIAL RECONSTRUCTION")
    print("    N=%d per arm, log-normal IIV on ductal tone capacity, PGE2 drive,")
    print("    ibuprofen CL and lung immaturity.\n")
    trial_reconstruction(n=140)
    return results


# ==============================================================================
#  SECTION 11.  VIRTUAL POPULATION
# ==============================================================================

# Inter-individual variability, as log-normal CVs on the parameters that
# actually decide an individual's fate.  DMAX_INT and NET50 matter most: the
# duct's unconstricted size and the constrictor threshold together determine
# whether a given infant sits above or below the structural occlusion limit,
# which is why a population contains 24-week responders and 28-week failures
# even though the deterministic model says 24 wk cannot close and 28 wk can.
IIV = {"DMAX_INT": 0.09,     # unconstricted ductal diameter
       "NET50": 0.10,        # constrictor threshold
       "KSYN_PGE2D": 0.16,   # ductal PGE2 drive
       "CL_IBU0": 0.32,      # ibuprofen clearance
       "CRS0": 0.14,         # lung immaturity
       "EC50_PGE2": 0.22,    # EP4 sensitivity
       "H_BPD0": 0.20}       # unmodelled BPD risk


def sample_subject(rng, ga_dist=None, **over):
    p = default_params()
    if ga_dist is not None:
        p["GA"] = float(rng.choice(ga_dist[0], p=ga_dist[1]))
    p["BW"] = max(0.42, 0.5 + 0.11 * (p["GA"] - 23) + rng.normal(0, 0.09))
    for k, cv in IIV.items():
        p[k] = p[k] * math.exp(rng.normal(0, cv))
    p.update(over)
    return p


def arm(rng, n, regimen_fn, horizon_d=80, ga_dist=None, require_large=True,
        **over):
    """Simulate an arm.  require_large mimics the trials' echocardiographic
    entry criterion (large PDA on the screening echo) by rejecting subjects
    whose duct is already small at the screening time."""
    rows = []
    tries = 0
    while len(rows) < n and tries < n * 6:
        tries += 1
        p = sample_subject(rng, ga_dist=ga_dist, **over)
        doses = regimen_fn(p)
        o = simulate(p, doses, 24 * horizon_d, dt=2.0)
        if require_large and at_day(o, "DDUCT", 2.0) < 1.5:
            continue
        rows.append((p, o))
    return rows


def summarise(rows, tag):
    n = len(rows)
    cl3w = np.mean([1.0 if sustained_closure(p, p, 0) else 0 for p, _ in []]) if False else None
    closed_3w = np.mean([1.0 if (closure_day(o, p=p) is not None
                                 and closure_day(o, p=p) <= 21) else 0.0
                         for p, o in rows])
    reop = [reopened(o, p) for p, o in rows]
    bur = np.mean([at_day(o, "PDABUR", 36) for _, o in rows])
    p_bpd = np.mean([at_day(o, "P_BPD", 70) for _, o in rows])
    p_dth = np.mean([at_day(o, "P_DEATH", 70) for _, o in rows])
    p_nec = np.mean([at_day(o, "P_NEC", 70) for _, o in rows])
    p_ivh = np.mean([at_day(o, "P_IVH", 70) for _, o in rows])
    # composite: death OR moderate/severe BPD  (independent-hazard combination)
    comp = np.mean([1.0 - (1.0 - at_day(o, "P_BPD", 70)) * (1.0 - at_day(o, "P_DEATH", 70))
                    for _, o in rows])
    comp3 = np.mean([1.0 - (1 - at_day(o, "P_BPD", 70)) * (1 - at_day(o, "P_DEATH", 70))
                     * (1 - at_day(o, "P_NEC", 70)) for _, o in rows])
    print("    %-30s n=%3d  closed<=3wk %5.1f%%  reopen %5.1f%%  burden %5.0f"
          % (tag, n, 100 * closed_3w, 100 * np.mean(reop), bur))
    print("    %-30s        death/BPD %5.1f%%   +NEC %5.1f%%   IVH %5.1f%%"
          % ("", 100 * comp, 100 * comp3, 100 * p_ivh))
    return dict(n=n, closed=closed_3w, reopen=np.mean(reop), burden=bur,
                comp=comp, comp3=comp3, ivh=p_ivh, nec=p_nec, death=p_dth,
                bpd=p_bpd)


def trial_reconstruction(n=140):
    rng = np.random.default_rng(20260804)
    # Baby-OSCAR enrolled 23-28 wk with a large PDA within 72 h.
    ga_bo = (np.array([23., 24., 25., 26., 27., 28.]),
             np.array([0.06, 0.15, 0.20, 0.22, 0.20, 0.17]))

    print("  Baby-OSCAR (Gupta 2024 NEJM; 23-28 wk, large PDA <72 h)")
    print("    observed: death or moderate/severe BPD 69.4% ibuprofen"
          " vs 63.5% expectant")
    a = arm(rng, n, lambda p: [], ga_dist=ga_bo)
    b = arm(rng, n, lambda p: rx_ibuprofen(48.0), ga_dist=ga_bo)
    ra = summarise(a, "expectant")
    rb = summarise(b, "early ibuprofen d2")
    print("    modelled risk difference (ibu - expectant): %+.1f points"
          % (100 * (rb["comp"] - ra["comp"])))

    print("\n  BeNeDuctus (Hundscheid 2023 NEJM; <28 wk, composite NEC/BPD/death)")
    print("    observed: 46.5% expectant vs 63.5% early ibuprofen (expectant"
          " non-inferior)")
    print("    modelled: %.1f%% expectant vs %.1f%% ibuprofen  (diff %+.1f)"
          % (100 * ra["comp3"], 100 * rb["comp3"],
             100 * (rb["comp3"] - ra["comp3"])))

    print("\n  TIPP (Schmidt 2001 NEJM; prophylactic indomethacin in ELBW)")
    print("    observed: severe IVH 9% vs 13%; no change in death/disability")
    c = arm(rng, max(60, n // 2), lambda p: rx_indomethacin(8.0, 0.1, 0.1, 3),
            ga_dist=ga_bo)
    rc = summarise(c, "prophylactic indomethacin")
    print("    modelled severe IVH: %.1f%% prophylaxis vs %.1f%% expectant"
          % (100 * rc["ivh"], 100 * ra["ivh"]))
    print("    modelled death/BPD:  %.1f%% prophylaxis vs %.1f%% expectant"
          % (100 * rc["comp"], 100 * ra["comp"]))

    print("\n  Closure rate by gestational age, standard ibuprofen from d2")
    print("    (pooled literature: ~70-80%% at 28-30 wk, ~30-50%% below 26 wk)")
    for ga in [24, 26, 28, 30]:
        rows = arm(rng, 45, lambda p: rx_ibuprofen(48.0), horizon_d=30,
                   GA=float(ga), require_large=False)
        cl = np.mean([1.0 if (closure_day(o, p=p) is not None
                              and closure_day(o, p=p) <= 8) else 0.0
                      for p, o in rows])
        rp = np.mean([reopened(o, p) for p, o in rows])
        print("      %d wk   closure by d8 %5.1f%%   reopening %5.1f%%"
              % (ga, 100 * cl, 100 * rp))


if __name__ == "__main__":
    run_all(quick="--quick" in sys.argv)
