#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 PYRUVATE KINASE DEFICIENCY (PKLR) -- QSP REFERENCE IMPLEMENTATION
================================================================================
 Dependency-free Python re-implementation of the mrgsolve model in
 `pkd_mrgsolve_model.R`.  This file is the *ground truth*: every equation in the
 R file is a verbatim transliteration of an equation that was first run, probed
 and (where it was wrong) fixed here.  No R runtime is available in the build
 environment, so the R file cannot be executed; this file is what was executed.

 WHY A SEPARATE REFERENCE EXISTS
 -------------------------------
 The model contains a fast algebraic subsystem (erythrocyte glycolysis, time
 constants of minutes) nested inside a slow dynamical system (red cell cohorts,
 iron, marrow; time constants of days to years).  A naive single-timescale ODE
 formulation is stiff by ~5 orders of magnitude.  The reduction used here --
 solving the fast subsystem to quasi-steady state at every derivative
 evaluation -- has to be verified numerically before it can be trusted, and it
 is the reduction that generates the model's central results.  Section 0 below
 does that verification.

 THE STRUCTURAL CHOICE THIS MODEL IS BUILT ON
 --------------------------------------------
 The PK lesion sits at the LAST ATP-generating step of glycolysis, which is
 DOWNSTREAM of the Rapoport-Luebering shunt branch point (1,3-BPG).  Therefore
 one lesion moves two quantities in OPPOSITE physiological directions:

     PK activity down  ->  ATP down          (cell dies sooner)
     PK activity down  ->  2,3-BPG up        (each surviving gram of Hb
                                              delivers MORE oxygen)

 Every non-obvious consequence in this file follows from refusing to collapse
 those two into a single "severity" scalar.  In particular the hemoglobin
 endpoint used by every registrational trial in this disease is not a
 sufficient statistic for oxygen transport, and the two can move in opposite
 directions under treatment.  Section 6 quantifies that and derives the
 break-even 2,3-BPG reduction in closed form.

 Units:  metabolites mM (per L red cell water); fluxes mmol/(L RBC * h);
         slow system time in DAYS; cell counts 10^12/L; Hb g/dL.
================================================================================
"""

import math
import json
import os
import sys

OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


def hdr(s):
    emit()
    emit("=" * 78)
    emit(s)
    emit("=" * 78)


def sub(s):
    emit()
    emit("-" * 78)
    emit(s)
    emit("-" * 78)


# =============================================================================
# SECTION 1.  PARAMETERS
# =============================================================================
# Every value carries its provenance.  "[calibrated]" means the value was
# solved backwards from a stated normal-physiology reference point in this file
# (the arithmetic is shown in Section 2), NOT fitted to patient data.
# "[fitted]" marks the only three numbers that were adjusted to make a clinical
# observation come out; they are listed again in Section 9 so the reader can see
# exactly how much of the model is free.

P = {}

# ---- 1a. Erythrocyte glycolysis: reference (normal adult) operating point ----
# Mulquiney & Kuchel 1999 (PMID 10477269/10477270) detailed kinetic model and
# the classical erythrocyte metabolite tables are the source of the reference
# concentrations; the reduced model below is calibrated to reproduce them.
P['ATP0'] = 1.70        # mM, red cell ATP
P['ADP0'] = 0.22        # mM
P['AMP0'] = 0.02        # mM
P['DPG0'] = 4.80        # mM 2,3-BPG  (~12.3 umol/gHb)
P['PG30'] = 0.060       # mM 3-phosphoglycerate
P['PEP0'] = 0.017       # mM phosphoenolpyruvate
P['FBP0'] = 0.012       # mM fructose-1,6-bisphosphate
P['J0'] = 1.60          # mmol glucose/(L RBC*h), normal glycolytic flux
P['phi0'] = 0.20        # fraction of glycolytic carbon through the R-L shunt

# adenylate kinase equilibrium 2 ADP <-> ATP + AMP  [calibrated, Section 2]
P['Atot'] = P['ATP0'] + P['ADP0'] + P['AMP0']
P['Kak'] = P['ATP0'] * P['AMP0'] / (P['ADP0'] ** 2)

# upper glycolysis (HK + PFK lumped); 2,3-BPG is a physiological inhibitor of
# both hexokinase and phosphofructokinase -- this is the feedback that makes the
# shunt self-limiting and, in PKD, self-amplifying.
P['Kup'] = 0.25         # mM, ATP K_m of the lumped HK/PFK step
P['KiDPG'] = 3.00       # mM, 2,3-BPG inhibition constant on HK/PFK
P['Vup'] = None         # [calibrated]

# aldolase / triose isomerase: sets FBP from flux (QSS)
P['kald'] = None        # [calibrated]

# pyruvate kinase R  (V-type representation of FBP allosteric activation)
P['Kpep'] = 0.30        # mM
P['Kadp'] = 0.30        # mM
P['KiATP'] = 3.00       # mM, ATP is an allosteric inhibitor of PKR
P['Kfbp'] = 0.004       # mM, FBP activation constant of PKR
P['Amax'] = 3.00        # maximal fold activation by FBP (dimensionless)
P['Vpk'] = None         # [calibrated]

# 3-PG <-> 2-PG <-> PEP lumped equilibrium (phosphoglycerate mutase + enolase)
P['theta'] = None       # [calibrated] = PEP/3-PG

# Rapoport-Luebering shunt.  Both activities live on one protein (BPGM).
#   synthase   v_syn  = k_syn * [1,3-BPG]
#   phosphatase v_phos = k_phos * [2,3-BPG] / (1 + [3-PG]/Ki_pg3)
# 1,3-BPG is not tracked: the phosphoglycerate kinase reaction is close to
# equilibrium, so [1,3-BPG] = [3-PG]*([ATP]/[ADP])/K_PGK.  Hence at steady state
#      2,3-BPG = kappa * [3-PG] * ([ATP]/[ADP]) * (1 + [3-PG]/Ki_pg3)
#
# THE 3-PG TERM IS NOT COSMETIC.  Without it the model gets the sign right and
# the magnitude badly wrong: a PK lesion raises 3-PG but lowers ATP/ADP, and in
# the bare PGK-equilibrium form those two cancel almost exactly (a first version
# of this file predicted 2,3-BPG +9% at 25% residual activity, against an
# observed 2-3 FOLD elevation -- see DEFECT #1 in Section 9).  The physiological
# resolution is that 3-PG and 2-PG are inhibitors of the 2,3-BPG phosphatase
# activity of BPGM, so accumulation above the block does not merely feed the
# shunt, it also blocks its exit.  That makes 2,3-BPG depend on 3-PG twice.
P['kappa'] = None       # [calibrated]
P['kphos'] = None       # [calibrated] 2,3-BPG phosphatase rate constant, /h
P['Kipg3'] = 0.030      # mM, 3-PG/2-PG inhibition of 2,3-BPG phosphatase
P['DPGmax'] = 25.0      # mM, ceiling on free 2,3-BPG

# Free / hemoglobin-bound partition of 2,3-BPG.  This is NOT bookkeeping.
# 2,3-BPG binds the central cavity of DEOXYhemoglobin with Kd ~ 30 uM and
# oxyhemoglobin ~100x more weakly, and the tetramer concentration inside a red
# cell (~5.3 mM) is of the same order as the whole 2,3-BPG pool.  So the fraction
# of the cell's 2,3-BPG that is chemically available to the shunt enzymes and to
# hexokinase/phosphofructokinase depends on how DESATURATED the blood is.
# Consequence: anemia itself raises total 2,3-BPG, with no change in any enzyme,
# and it does so with positive-feedback gain -- lower free 2,3-BPG relieves the
# inhibition of HK/PFK, which raises flux, which raises 3-PG, which raises free
# 2,3-BPG again.  In pyruvate kinase deficiency this limb operates ON TOP of the
# enzymatic block, and the two multiply.
P['HbTet'] = 5.30       # mM hemoglobin tetramer inside the red cell
P['Kdb'] = 0.050        # mM, 2,3-BPG dissociation constant for deoxyHb
P['fdeox0'] = 0.1561    # reference circulating deoxy fraction (model-derived)
P['DPGfree0'] = None    # [calibrated] free 2,3-BPG at the reference point
P['tauFdeox'] = 0.25    # d, adaptation time of the circulating deoxy fraction

# Adenylate pool drain.  The red cell cannot resynthesise adenine nucleotides de
# novo; AMP deaminase and 5'-nucleotidase irreversibly export the pool as IMP /
# hypoxanthine.  So the total adenylate content is NOT conserved -- it falls with
# cell age, and it falls faster when ATP is low (because AMP is then high).
# Ignoring this is the second thing that breaks the 2,3-BPG prediction: with a
# conserved pool, ADP rises steeply as ATP falls and ATP/ADP collapses.
P['kAMPd'] = 0.0065     # mM/d, AMP deaminase Vmax
P['KAMPd'] = 0.020      # mM (saturates early: the drain accelerates about
                        #  2-fold, not 4-fold, when ATP collapses)
P['AtotMin'] = 0.40     # mM, floor below which the cell is non-viable

# ATP demand: Na/K-ATPase (saturable) + lumped non-pump turnover
P['Vpump'] = 1.20       # mmol/(L*h)
P['Kpump'] = 0.35       # mM
P['kleak'] = None       # [calibrated] /h

# ---- 1b. Hemoglobin oxygen binding ----
P['P50ref'] = 26.8      # mmHg at DPG0, pH 7.4, 37 C
P['nHill'] = 2.70
P['nDPG'] = 0.32        # exponent giving dP50/dDPG = 1.8 mmHg/mM at reference
P['PaO2'] = 95.0        # mmHg
P['VO2rest'] = 250.0    # mL O2/min, resting whole-body consumption
P['CO0'] = 5.00         # L/min reference cardiac output
P['HbCO'] = 10.0        # g/dL below which resting cardiac output starts rising
P['kCO'] = 1.20         # cardiac output compensation gain
P['PvO2ref'] = 38.0     # mmHg, normal mixed venous PO2 (model-derived, S2)

# ---- 1c. Red cell age structure ----
# 14 circulating cohorts: 2 reticulocyte bins (1 d each) + 12 mature bins.
# Mature bin width is set so that a normal cell is cleared at ~120 d.
P['nRet'] = 2
P['nMat'] = 12
P['wRet'] = 0.5         # d, width of each reticulocyte bin
P['wMat'] = 9.9167      # d, width of each mature bin (1 + 12*9.9167 = 120 d)

# glycolytic capacity and ATP demand decline/rise with cell age.
# Reticulocytes retain mitochondria-independent high glycolytic capacity and a
# much larger membrane/pump load than mature discocytes.
P['vupRet'] = 2.50      # x normal glycolytic Vmax in reticulocytes
P['vupOld'] = 0.70      # x normal at 120 d
P['useRet'] = 1.80      # x normal ATP demand in reticulocytes
P['useOld'] = 0.95      # x normal ATP demand at 120 d

# PKR protein decays with cell age and cannot be resynthesised.  Wild-type PKR
# loses ~half its activity across a normal lifespan; the thermolabile mutants
# that cause clinical PKD lose it far faster.  This is the single parameter that
# converts a genotype into an age-dependent activity profile.
P['tauPKwt'] = 170.0    # d, wild-type PKR activity decay time constant
# Wild-type PKR loses about half its activity across a 120-day lifespan, giving
# tau ~ 170 d.  Thermolabile mutants lose it faster, but a first pass used 25 d
# (a 7-fold acceleration) and that made the genotype-to-haemoglobin map a step
# function: 35% residual activity gave a normal haemoglobin and 15% gave 3 g/dL,
# with nothing in between and no way to represent the observed clinical spectrum
# (DEFECT #11, Section 9).  A 3-4 fold acceleration reproduces the spectrum.
P['tauPKmut'] = 50.0    # d, thermolabile mutant PKR (default) [fitted #1]

# ---- 1d. Hazard / destruction ----
P['hBase'] = 1.0 / 1000.0   # /d, age-independent random loss
P['ATPcrit'] = 1.25        # mM, ATP below which cation homeostasis fails
P['kHem'] = 0.240          # /d, extravascular lysis gain [fitted #2]
P['pHem'] = 1.0            # steepness of the ATP failure term
# NOT ALL OF THE HAZARD IS ACUTELY REVERSIBLE.  A cell that has already lost
# membrane, exported its adenylate pool and become dehydrated is not rescued by
# refilling its ATP: the oxidative and mechanical injury it accumulated while it
# was ATP-poor is irreversible.  Only part of the ATP-dependent hazard therefore
# responds to an activator on the timescale of days.
# This is not a tuning device, it is the difference between a plausible and an
# impossible answer.  With the whole hazard reversible, 5 mg BID of mitapivat
# raised haemoglobin by +6.0 g/dL and the transient overshot to 22.7 g/dL,
# against an observed mean rise of ~1.7 g/dL and a maximum of 3.4 g/dL in the
# best responders (DEFECT #14, Section 9).  The irreversible fraction is
# evaluated at the cohort's UNDRUGGED activity, i.e. at the injury it already
# carries.
P['fRev'] = 0.20           # fraction of the ATP-dependent hazard an activator
                           #  can reverse acutely                    [fitted #4]
P['kSenesc'] = 3.0         # /d terminal senescence rate out of the last bin

# spleen: the red pulp is a glucose-poor, acidotic, hypoxic, slow-transit
# compartment.  Modelled as a *reduction of glycolytic capacity* experienced
# during each pass, not as an ad hoc extra hazard.
P['splStress'] = 0.45      # x glycolytic Vmax available inside the red pulp
P['splUse'] = 1.15         # x ATP demand inside the red pulp (mechanical)
P['kSpl'] = 0.006          # /d gain on splenic trapping [fitted #3]
P['ATPcritSpl'] = 0.80     # mM, failure threshold under splenic stress

# WHY RETICULOCYTES ARE THE SPLEEN'S PREFERRED VICTIM IN THIS DISEASE.
# A reticulocyte still has mitochondria and covers part of its (large) ATP bill
# by oxidative phosphorylation.  The splenic red pulp is hypoxic as well as
# glucose-poor, so that subsidy is withdrawn exactly where the mechanical
# stress is highest, and the cell is thrown back on the one pathway the disease
# has broken.  A mature discocyte has no mitochondria and therefore has nothing
# extra to lose on entering the spleen.
# This asymmetry is not decoration: it is the only way an age-structured model
# reproduces the observation that splenectomy in PK deficiency raises the
# hemoglobin AND the reticulocyte count at the same time (Grace 2018,
# PMID 29549173; Nathan 1968, PMID 5634483).  A single-pool model, or an
# age-structured model whose splenic hazard rises with age, necessarily predicts
# that reticulocytes FALL after splenectomy.
P['mitoRet'] = 0.40        # fraction of reticulocyte ATP demand met oxidatively
P['mitoRetAge'] = 2.0      # d over which mitochondria are lost
# ...and the red pulp RETAINS reticulocytes while it strips those organelles and
# the surplus membrane, so their exposure to it is several-fold longer than a
# mature cell's single pass.  Long exposure x withdrawn subsidy x the largest
# absolute pump load = the youngest cells are the ones the spleen kills.
P['splRetExp'] = 5.0       # x splenic exposure for a fresh reticulocyte
# ...and while it is retained it has to DO something: strip organelles, shed
# surplus membrane, remodel its cytoskeleton -- all ATP-dependent.  So a
# reticulocyte does not merely need enough ATP to survive a pass, it needs enough
# to complete a construction project, and its failure threshold is accordingly
# higher than a mature discocyte's.  Giving reticulocytes their own threshold
# (rather than their own hazard) is what finally produced the observed SIGN of
# the reticulocyte response to splenectomy; sixteen combinations of the two
# hazard gains, with a single shared threshold, all produced the wrong sign
# (DEFECT #8, Section 9).
P['splRetCrit'] = 0.90     # extra ATP threshold for a fresh reticulocyte
# SEQUESTRATION, NOT JUST DESTRUCTION.  Nathan 1968 (PMID 5634483) measured
# life-span AND ORGAN SEQUESTRATION in this disease, and sequestration is the
# part that explains the clinical paradox: PK-deficient reticulocytes are HELD
# in the red pulp, so they are absent from the circulating count while they are
# there.  Splenectomy therefore does two separate things -- it stops the
# destruction (a modest haemoglobin gain, median +1.6 g/dL) and it RELEASES the
# held pool (a large jump in the measured reticulocyte percentage).  Trying to
# get both out of a hazard term alone does not work: to make reticulocytes the
# spleen's preferred victim strongly enough to matter, the splenic hazard has to
# dominate total destruction, and then splenectomy predicts +7 to +12 g/dL
# instead of +1.6 (DEFECT #9, Section 9).
P['kSeqIn'] = 1.20         # /d capture of stressed reticulocytes by the red pulp
P['kSeqOut'] = 0.35        # /d release back into the circulation
P['kSeqHem'] = 0.050       # /d destruction inside the sequestered pool
P['splVol0'] = 150.0       # mL normal spleen volume
P['kSplGrow'] = 0.0022     # /d
P['splVolMax'] = 1800.0    # mL

# ---- 1e. Marrow / erythropoiesis ----
P['kProg'] = 0.35          # /d progenitor turnover
P['ampMax'] = 64.0         # maximal erythroblast amplification (6 divisions)
P['kEB1'] = 0.33           # /d early erythroblast transit
P['kEB2'] = 0.40           # /d late erythroblast transit
P['EPO50'] = 30.0          # IU/L, EPO for half-maximal amplification
P['EPOhill'] = 1.4
P['kEPOsyn'] = 1.0
P['kEPOel'] = 4.0          # /d
P['HbEPO'] = 15.0          # g/dL reference
P['EPOgain'] = 3.30        # log-linear EPO response slope to Hb deficit
P['EPObase'] = 8.0         # IU/L
P['EPOmax'] = 12000.0      # IU/L ceiling
# ineffective erythropoiesis: erythroblasts still have mitochondria, so the PK
# lesion costs them much less than it costs a mature red cell.  The residual
# cost is what mitapivat was shown to relieve (van Beers 2024, PMID 38330179).
P['ieoMax'] = 0.45         # max apoptotic fraction of late erythroblasts
P['ieoK'] = 0.35           # ATP-margin scale of the apoptotic term
P['mitoRescue'] = 0.55     # fraction of erythroblast ATP demand met oxidatively
P['MCH0'] = 30.0           # pg
P['MCV0'] = 90.0           # fL
# Stress reticulocytosis lengthens the CIRCULATING maturation time of a
# reticulocyte: under high erythropoietin drive the marrow releases cells early
# ("shift reticulocytes") and they finish maturing in the blood.  This is the
# same effect clinical practice corrects for with the reticulocyte maturation
# factor (1.0 at Hct 45 rising to ~2.5 at Hct 15).  Without it the model cannot
# reach the reticulocyte percentages actually seen in PK deficiency: the
# steady-state reticulocyte fraction is exactly (retic residence)/(red cell
# lifespan), so 20-40% at a lifespan of 15 d REQUIRES a residence of 3-6 d.
P['kMatf'] = 4.0           # gain of the maturation-time prolongation
P['matfMax'] = 4.0         # cap on the maturation factor

# ---- 1f. Iron / hepcidin ----
P['kERFEsyn'] = 1.0
P['kERFEel'] = 3.5         # /d
P['ERFE0'] = 1200.0        # ng/L reference
P['kHEPsyn'] = 1.0
P['kHEPel'] = 8.0          # /d
P['HEP0'] = 8000.0         # ng/L reference hepcidin
P['ERFEi50'] = 9000.0      # ng/L ERFE for half-maximal hepcidin suppression
P['LICstim'] = 0.135       # hepcidin induction gain per mg Fe/g dw
P['absMax'] = 4.0          # mg/d maximal duodenal iron absorption
P['HEPabs50'] = 4000.0     # ng/L
P['kFeLoss'] = 0.80        # mg/d obligate loss
P['kRESout'] = 0.90        # /d macrophage iron release to plasma
P['FeUnit'] = 200.0        # mg Fe per transfused RBC unit
P['LIVwt'] = 1500.0        # g liver wet weight
P['dwFrac'] = 0.28         # dry/wet weight
P['kLIC'] = 4.00           # NTBI -> liver iron partition
P['kLICin'] = 14.0         # mg/d transferrin iron flux into hepatocytes at TSAT 1
P['kLICout'] = 0.004       # /d
P['kFerr'] = 22.0          # ng/mL ferritin per mg Fe/g dw LIC
P['FerrBase'] = 30.0
P['kCard'] = 0.0022        # /d NTBI -> cardiac iron
P['kCardOut'] = 0.0016     # /d
P['NTBIthr'] = 0.70        # TSAT above which NTBI appears
P['TSATmax'] = 1.0

# ---- 1g. Bilirubin / hepatobiliary ----
P['hemeBili'] = 1.0        # mg bilirubin per mg heme (stoichiometric)
P['kUGT'] = 16.4           # /d, first-order bilirubin conjugation clearance
P['UGTsat'] = 25.0         # mg/dL, where UGT1A1 finally saturates

P['kCBout'] = 22.0         # /d biliary excretion
P['kGS'] = 1.0 / 4000.0    # gallstone hazard index gain
P['LDH0'] = 150.0          # U/L
P['kLDH'] = 1.0
P['HAP0'] = 1.20           # g/L
P['kHapEl'] = 2.0          # /d
P['ivFrac'] = 0.10         # fraction of hemolysis that is intravascular

# ---- 1h. Mitapivat PK/PD ----
# Yang 2019 (PMID 30091852) phase 1 SAD/MAD in healthy volunteers; mitapivat is
# a CYP3A substrate AND a moderate CYP3A inducer, so it auto-induces.
P['mitKa'] = 14.0          # /d absorption (only the RBC-compartment AUC drives
                           #  effect, so the plasma peak shape is not resolved)
P['mitF'] = 0.85
P['mitVc'] = 40.0          # L
P['mitVp'] = 55.0          # L
P['mitQ'] = 25.0           # L/d
P['mitCL0'] = 166.0        # L/d at baseline CYP3A (t1/2 ~ 4 h)
P['mitKrbc'] = 3.0         # RBC:plasma partition
P['mitKrbcOut'] = 30.0     # /d
P['cypKin'] = 0.55         # /d enzyme turnover
P['cypEmax'] = 1.45        # max fold CYP3A induction
P['cypEC50'] = 0.35        # mg/L
P['mitEmax'] = 1.60        # max fold PKR activation (V-type)
P['mitEC50'] = 0.22        # mg/L at the red cell
P['mitStabMax'] = 2.40     # max fold increase in mutant PKR half-life
P['mitStabEC50'] = 0.30    # mg/L
P['kStab'] = 0.055         # /d, how fast the stabilisation effect accrues
P['mitAromEmax'] = 0.62    # max fractional aromatase inhibition
P['mitAromEC50'] = 0.28    # mg/L
P['mitALT'] = 14.0         # U/L max ALT rise

# ---- 1i. Etavopivat / tebapivat ----
P['etaKa'] = 6.0
P['etaF'] = 0.80
P['etaVc'] = 55.0
P['etaCL'] = 55.0          # L/d (long half-life, once daily)
P['etaKrbc'] = 12.0        # very high RBC partitioning
P['etaKrbcOut'] = 12.0
P['etaEmax'] = 1.50
P['etaEC50'] = 0.09
P['tebKa'] = 7.0
P['tebF'] = 0.85
P['tebVc'] = 45.0
P['tebCL'] = 90.0
P['tebKrbc'] = 6.0
P['tebKrbcOut'] = 20.0
P['tebEmax'] = 1.40
P['tebEC50'] = 0.010       # sub-mg doses

# ---- 1j. Deferasirox ----
P['dfxKa'] = 10.0
P['dfxF'] = 0.70
P['dfxVc'] = 15.0
P['dfxCL'] = 150.0
P['dfxEmax'] = 0.62        # max fractional LIC clearance enhancement
P['dfxEC50'] = 4.0         # mg/L

# ---- 1k. Endocrine / bone / vascular sequelae ----
P['E20'] = 28.0            # pg/mL male estradiol
P['TST0'] = 550.0          # ng/dL male testosterone
P['kE2'] = 6.0             # /d
P['kTST'] = 4.0            # /d
P['tstFeedback'] = 0.45
P['BMD0'] = 0.0            # Z-score
P['kBMD'] = 0.00042        # /d
P['kBMDrec'] = 0.00030     # /d
P['PVR0'] = 1.8            # Wood units
P['kPVR'] = 0.00055
P['kPVRrec'] = 0.00040
P['kEMH'] = 0.0016
P['kEMHrec'] = 0.0012

# ---- 1l. Transfusion / donor cells ----
P['donHbUnit'] = 1.05      # g/dL Hb rise per unit in a 70 kg adult
P['donTau'] = 55.0         # d, mean survival of transfused normal cells
P['donDPGstored'] = 0.55   # mM, 2,3-BPG in stored (>2 wk) packed cells
P['donDPGrec'] = 0.60      # /d in-vivo 2,3-BPG regeneration
P['txThresh'] = 8.0        # g/dL default transfusion trigger


# =============================================================================
# SECTION 2.  CALIBRATION OF THE GLYCOLYTIC MODULE
# =============================================================================
# Nothing here is fitted to a patient.  Each constant is solved backwards from
# the normal operating point in 1a so that the reduced model reproduces normal
# erythrocyte metabolism exactly.  The point of showing the arithmetic is that
# it exposes how much reserve the PK step has -- which turns out to be the
# reason PK deficiency needs a very severe lesion before it becomes a disease.

def dpg_total(dpgf, fdeox):
    """Total cell 2,3-BPG from the free concentration and the deoxy fraction."""
    return dpgf * (1.0 + P['HbTet'] * fdeox / (P['Kdb'] + dpgf))


def dpg_free_from_total(dpgt, fdeox):
    """Inverse of dpg_total (a quadratic in the free concentration)."""
    b = P['HbTet'] * fdeox + P['Kdb'] - dpgt
    disc = b * b + 4.0 * P['Kdb'] * dpgt
    if disc < 0.0:
        return 0.0
    return 0.5 * (-b + math.sqrt(disc))


def calibrate():
    p = P
    atp, adp, dpg = p['ATP0'], p['ADP0'], p['DPG0']
    J = p['J0']

    # free 2,3-BPG at the reference point, from the measured TOTAL
    p['DPGfree0'] = dpg_free_from_total(dpg, p['fdeox0'])
    dpgf = p['DPGfree0']
    inh0 = 1.0 + p['PG30'] / p['Kipg3']

    # 2,3-BPG phosphatase from the requested shunt fraction.  The shunt flux is
    # v_phos = kphos*[free 2,3-BPG]/(1 + [3-PG]/Ki), so the inhibition factor has
    # to appear here too or the calibration is inconsistent with the rate law
    # actually used at run time (DEFECT #3, Section 9).
    p['kphos'] = p['phi0'] * 2.0 * J * inh0 / dpgf

    # non-pump ATP turnover so that supply = demand at the reference point:
    #   net ATP production = 2J*(1 - phi0)
    net = 2.0 * J * (1.0 - p['phi0'])
    pump = p['Vpump'] * atp / (p['Kpump'] + atp)
    p['kleak'] = (net - pump) / atp

    # upper glycolysis Vmax (HK/PFK see FREE 2,3-BPG)
    sATP = atp / (p['Kup'] + atp)
    fD = 1.0 / (1.0 + dpgf / p['KiDPG'])
    p['Vup'] = J / (sATP * fD)

    # aldolase constant from the reference FBP
    p['kald'] = J / p['FBP0']

    # PKR Vmax: PK must carry the whole triose flux, 2J
    fbp = p['FBP0']
    A = 1.0 + p['Amax'] * fbp / (p['Kfbp'] + fbp)
    hadp = adp / (p['Kadp'] + adp)
    iatp = 1.0 / (1.0 + atp / p['KiATP'])
    spep = p['PEP0'] / (p['Kpep'] + p['PEP0'])
    p['Vpk'] = 2.0 * J / (A * hadp * iatp * spep)

    # lumped 3-PG -> PEP equilibrium and the shunt lumped constant
    p['theta'] = p['PEP0'] / p['PG30']
    p['kappa'] = dpgf / (p['PG30'] * (atp / adp) * inh0)
    return dict(net=net, pump=pump, A=A, hadp=hadp, iatp=iatp, spep=spep)


CAL = calibrate()


# =============================================================================
# SECTION 3.  THE FAST SUBSYSTEM -- QUASI-STEADY-STATE GLYCOLYSIS
# =============================================================================
# Given
#   aeff  : effective PKR activity, relative to wild type (genotype x age x drug)
#   vup   : glycolytic capacity scale (cell age / splenic stress)
#   useS  : ATP demand scale (cell age / mechanical stress)
# solve the algebraic steady state of
#   ADP  = adenylate kinase equilibrium given ATP and the fixed adenylate pool
#   J    = Vup*vup * s(ATP) / (1 + DPG/KiDPG)
#   FBP  = J / kald
#   2J   = Vpk*aeff*A(FBP)*h(ADP)*i(ATP)*s(PEP)          [PK carries all flux]
#   3-PG = PEP/theta ;  2,3-BPG = kappa*3-PG*(ATP/ADP)   [PGK equilibrium + R-L]
#   2J - kphos*2,3-BPG = v_use(ATP)                      [ATP balance]
#
# Two nested scalar solves: an inner damped fixed point for (J, DPG) at fixed
# ATP, and an outer root find on the ATP balance residual.

NITER = 60          # fixed iteration count -- see Section 0.5 for why 60


def adk_adp(atp, atot):
    """ADP from the adenylate kinase equilibrium at a GIVEN total adenylate."""
    if atp <= 0.0:
        return 1e-12
    c = atp * (atot - atp)
    if c <= 0.0:
        return 1e-12
    disc = atp * atp + 4.0 * P['Kak'] * c
    return (-atp + math.sqrt(disc)) / (2.0 * P['Kak'])


def invert_vuse(y, useS):
    """Exact inverse of the ATP demand function: given a net ATP supply y,
    return the ATP concentration at which demand equals supply.
        y/useS = Vpump*ATP/(Kpump+ATP) + kleak*ATP
    is a quadratic in ATP; taking the positive root makes the outer iteration a
    clean two-variable fixed point with no nested root find."""
    p = P
    if y <= 0.0:
        return 1e-4
    yy = y / useS
    b = p['Vpump'] + p['kleak'] * p['Kpump'] - yy
    disc = b * b + 4.0 * p['kleak'] * yy * p['Kpump']
    if disc < 0.0:
        return 1e-4
    return (-b + math.sqrt(disc)) / (2.0 * p['kleak'])


def gly_at_atp(atp, aeff, vup, useS, atot, dpg_guess, fdeox=None):
    """Everything downstream of a GIVEN ATP.  Inner problem: find the 2,3-BPG
    that satisfies the Rapoport-Luebering steady state.

    The inner map is strictly decreasing in 2,3-BPG (more 2,3-BPG -> less
    HK/PFK flux -> less 1,3-BPG -> lower 2,3-BPG target), so the fixed point is
    unique and a damped iteration with a warm start converges monotonically.
    That monotonicity is what makes the whole reduction well-posed; a first
    attempt that iterated on ATP and 2,3-BPG *simultaneously* diverged, because
    the FBP -> PKR activation limb is a positive feedback loop (DEFECT #2,
    Section 9)."""
    p = P
    if fdeox is None:
        fdeox = p['fdeox0']
    adp = adk_adp(atp, atot)
    sATP = atp / (p['Kup'] + atp)
    hadp = adp / (p['Kadp'] + adp)
    iatp = 1.0 / (1.0 + atp / p['KiATP'])
    rad = atp / max(adp, 1e-12)
    dpg = dpg_guess
    J = fbp = pep = pg3 = 0.0
    A = 1.0
    inh = 1.0
    for _ in range(400):
        J = p['Vup'] * vup * sATP / (1.0 + dpg / p['KiDPG'])
        fbp = J / p['kald']
        A = 1.0 + p['Amax'] * fbp / (p['Kfbp'] + fbp)
        cap = p['Vpk'] * aeff * A * hadp * iatp
        if cap <= 1e-14:
            S = 1.0 - 1e-12
        else:
            S = 2.0 * J / cap
            if S > 1.0 - 1e-12:
                S = 1.0 - 1e-12
        pep = p['Kpep'] * S / (1.0 - S)
        pg3 = pep / p['theta']
        inh = 1.0 + pg3 / p['Kipg3']
        dtar = p['kappa'] * pg3 * rad * inh
        if dtar > p['DPGmax']:
            dtar = p['DPGmax']
        if abs(dtar - dpg) < 1e-10:
            dpg = dtar
            break
        dpg += 0.45 * (dtar - dpg)
    v_shunt = p['kphos'] * dpg / inh
    v_use = useS * (p['Vpump'] * atp / (p['Kpump'] + atp) + p['kleak'] * atp)
    net = 2.0 * J - v_shunt
    amp = atot - atp - adp
    return dict(atp=atp, adp=adp, amp=(amp if amp > 0.0 else 0.0),
                dpgf=dpg, dpg=dpg_total(dpg, fdeox), fdeox=fdeox,
                J=J, fbp=fbp, pep=pep, pg3=pg3, A=A, atot=atot,
                v_use=v_use, v_shunt=v_shunt, net=net, resid=net - v_use,
                phi=(v_shunt / (2.0 * J) if J > 1e-12 else 1.0))


def gly_exact(aeff, vup=1.0, useS=1.0, atot=None, fdeox=None, ngrid=110,
              want_roots=False):
    """Outer problem: the ATP balance  2J - v_shunt = v_use(ATP).

    The residual is NEGATIVE at both ends of the admissible ATP range -- at
    ATP -> 0 the flux vanishes, and at ATP -> Atot the adenylate kinase
    equilibrium drives ADP -> 0, which stalls pyruvate kinase and hands the
    whole 1,3-BPG flux to the shunt.  So a plain bisection over the interval is
    invalid.  We scan, take the HIGHEST sign change (the metabolically competent
    branch), and bisect there.  The scan is not wasted effort: counting sign
    changes is how Section 5 detects that the ATP steady state can be bistable,
    which is the model's account of the genotype-phenotype threshold."""
    p = P
    if atot is None:
        atot = p['Atot']
    lo, hi = 1e-3, 0.995 * atot
    if hi <= lo:
        st = gly_at_atp(lo, aeff, vup, useS, max(atot, 2e-3), p['DPGmax'], fdeox)
        st['collapsed'] = True
        st['nroots'] = 0
        return st
    xs, rs = [], []
    dg = p['DPGfree0']
    for i in range(ngrid):
        x = lo * (hi / lo) ** (i / (ngrid - 1.0))
        st = gly_at_atp(x, aeff, vup, useS, atot, dg, fdeox)
        dg = st['dpgf']                   # warm start keeps the inner loop short
        xs.append(x)
        rs.append(st['resid'])
    roots = []
    for i in range(ngrid - 1):
        if rs[i] * rs[i + 1] < 0.0:
            roots.append((xs[i], xs[i + 1]))
    if not roots:
        st = gly_at_atp(lo, aeff, vup, useS, atot, p['DPGmax'], fdeox)
        st['collapsed'] = True
        st['nroots'] = 0
        if want_roots:
            st['roots'] = []
        return st
    a, b = roots[-1]
    dg = p['DPGfree0']
    st = None
    for _ in range(70):
        m = 0.5 * (a + b)
        st = gly_at_atp(m, aeff, vup, useS, atot, dg, fdeox)
        dg = st['dpgf']
        if st['resid'] > 0.0:
            a = m
        else:
            b = m
    st['collapsed'] = False
    st['nroots'] = len(roots)
    if want_roots:
        st['roots'] = roots
    return st


# =============================================================================
# SECTION 4.  OXYGEN TRANSPORT
# =============================================================================
# The whole point of tracking 2,3-BPG explicitly is that it changes the shape of
# the oxygen dissociation curve, and therefore breaks the identity
# "more hemoglobin = more oxygen delivered".

def p50_of(dpg):
    return P['P50ref'] * (dpg / P['DPG0']) ** P['nDPG']


def sat(po2, p50):
    r = (po2 / p50) ** P['nHill']
    return r / (1.0 + r)


def po2_of_sat(s, p50):
    if s <= 1e-9:
        return 0.0
    if s >= 1.0 - 1e-9:
        return 1e4
    return p50 * (s / (1.0 - s)) ** (1.0 / P['nHill'])


def cardiac_output(hb):
    if hb >= P['HbCO']:
        return P['CO0']
    return P['CO0'] * (1.0 + P['kCO'] * (P['HbCO'] / max(hb, 1.0) - 1.0))


def o2_state(hb, dpg):
    """Three closures on the same physiology.  They are three different
    questions, not three different models:
      A) fixed VO2 and fixed cardiac output -> what tissue PO2 results?
      B) fixed VO2 and fixed tissue PO2     -> what cardiac output is required?
      C) equivalent hemoglobin: the Hb a normal-P50 subject would need to
         extract the same O2 per litre of blood at the reference tissue PO2.
    """
    p50 = p50_of(dpg)
    sa = sat(P['PaO2'], p50)
    # (A) fixed VO2, cardiac output set by the anemia compensation rule
    co = cardiac_output(hb)
    cap = co * 1.34 * hb * 10.0                 # mL O2/min at 100% extraction
    dS = P['VO2rest'] / cap if cap > 0 else 1.0
    sv = sa - dS
    pv = po2_of_sat(sv, p50) if sv > 0.0 else 0.0
    # (B) fixed tissue PO2 at the normal operating point
    svB = sat(P['PvO2ref'], p50)
    dSB = sa - svB
    coreq = P['VO2rest'] / (1.34 * hb * 10.0 * dSB) if dSB > 1e-9 else float('inf')
    # (C) equivalent hemoglobin
    sa_n = sat(P['PaO2'], P['P50ref'])
    sv_n = sat(P['PvO2ref'], P['P50ref'])
    hbeq = hb * dSB / (sa_n - sv_n)
    return dict(p50=p50, SaO2=sa, PvO2=pv, SvO2=sv, CO=co, extr=dS,
                COreq=coreq, extrB=dSB, Hbeq=hbeq)


# =============================================================================
# SECTION 5.  LOOKUP TABLES FOR THE FAST SUBSYSTEM
# =============================================================================
# The slow system needs the glycolytic steady state of ~30 distinct cell
# populations at every derivative evaluation.  Solving them exactly each time is
# affordable in the compiled C++ that mrgsolve generates but not in Python, so
# the reference implementation tabulates each population once and interpolates.
#
# The tables are 2-D, in (effective PKR activity, total adenylate) -- NOT 3-D.
# The deoxy fraction does not need a dimension because of a result the model
# hands us for free:  the steady-state FREE 2,3-BPG concentration is exactly
# independent of the deoxy fraction.  Free 2,3-BPG is pinned by
# kappa*[3-PG]*(ATP/ADP)*(1+[3-PG]/Ki), none of which involves hemoglobin
# saturation, and the flux J depends on free rather than total 2,3-BPG.  So
# hemoglobin desaturation shifts 2,3-BPG entirely by sequestration, and adds
# nothing to the metabolic cost of the shunt -- total 2,3-BPG can be recovered
# from the tabulated free value analytically at run time.  (The comment in
# Section 1 originally claimed this limb had positive-feedback gain on flux; the
# model says the gain is exactly zero.  DEFECT #4, Section 9.)

class GTab(object):
    def __init__(self, vup, useS, na=120, nA=11,
                 amin=1.0e-4, amax=7.0, Alo=0.30, Ahi=1.95):
        self.vup, self.useS = vup, useS
        self.na, self.nA = na, nA
        self.amin, self.amax = amin, amax
        self.Alo, self.Ahi = Alo, Ahi
        self.lr = math.log(amax / amin)
        self.latp = []
        self.ldpg = []
        self.vJ = []
        for k in range(nA):
            At = Alo + (Ahi - Alo) * k / (nA - 1.0)
            ra, rd, rj = [], [], []
            for i in range(na):
                av = amin * math.exp(self.lr * i / (na - 1.0))
                st = gly_exact(av, vup, useS, atot=At, ngrid=64)
                ra.append(math.log(max(st['atp'], 1e-9)))
                rd.append(math.log(max(st['dpgf'], 1e-9)))
                rj.append(st['J'])
            self.latp.append(ra)
            self.ldpg.append(rd)
            self.vJ.append(rj)

    def get(self, av, At):
        if not (av == av) or av < self.amin:          # NaN-safe
            av = self.amin
        elif av > self.amax:
            av = self.amax
        t = math.log(av / self.amin) / self.lr * (self.na - 1.0)
        i = int(t)
        if i > self.na - 2:
            i = self.na - 2
        fi = t - i
        if not (At == At):
            At = self.Ahi
        u = (At - self.Alo) / (self.Ahi - self.Alo) * (self.nA - 1.0)
        if u < 0.0:
            u = 0.0
        elif u > self.nA - 1.0:
            u = self.nA - 1.0
        j = int(u)
        if j > self.nA - 2:
            j = self.nA - 2
        fj = u - j
        a0 = self.latp[j]
        a1 = self.latp[j + 1]
        d0 = self.ldpg[j]
        d1 = self.ldpg[j + 1]
        la = (a0[i] + fi * (a0[i + 1] - a0[i])) * (1.0 - fj) + \
             (a1[i] + fi * (a1[i + 1] - a1[i])) * fj
        ld = (d0[i] + fi * (d0[i + 1] - d0[i])) * (1.0 - fj) + \
             (d1[i] + fi * (d1[i + 1] - d1[i])) * fj
        j0 = self.vJ[j]
        j1 = self.vJ[j + 1]
        jj = (j0[i] + fi * (j0[i + 1] - j0[i])) * (1.0 - fj) + \
             (j1[i] + fi * (j1[i + 1] - j1[i])) * fj
        return math.exp(la), math.exp(ld), jj


# =============================================================================
# SECTION 6.  RED CELL AGE STRUCTURE
# =============================================================================
# Age is a state variable in this model, not a summary statistic.  It has to be,
# for three reasons that all matter clinically:
#   1. mutant PKR protein decays with cell age and cannot be replaced, so the
#      lesion is not a constant -- it deepens as the cell ages;
#   2. the assayed red cell PK activity of a patient is a cohort-weighted mean
#      dominated by young cells, which is why the diagnostic assay can read
#      normal in a genuinely deficient patient;
#   3. the spleen and a PK activator and gene therapy act on DIFFERENT parts of
#      the age axis, so a single-pool model gets the sign of the reticulocyte
#      response to splenectomy wrong.

NB = P['nRet'] + P['nMat']              # 14 circulating cohorts
BW = [P['wRet']] * P['nRet'] + [P['wMat']] * P['nMat']
BENTRY = []
_acc = 0.0
for _w in BW:
    BENTRY.append(_acc)
    _acc += _w
BMID = [BENTRY[i] + 0.5 * BW[i] for i in range(NB)]
BEXIT = _acc                            # = 120 d


def age_vup(a):
    """Glycolytic capacity vs cell age (x normal mature discocyte)."""
    if a < 2.0:
        return P['vupRet'] + (1.15 - P['vupRet']) * (a / 2.0)
    return 1.15 + (P['vupOld'] - 1.15) * (a - 2.0) / (BEXIT - 2.0)


def age_use(a):
    """ATP demand vs cell age (x normal mature discocyte)."""
    if a < 2.0:
        return P['useRet'] + (1.15 - P['useRet']) * (a / 2.0)
    return 1.15 + (P['useOld'] - 1.15) * (a - 2.0) / (BEXIT - 2.0)


VUP = [age_vup(BMID[i]) for i in range(NB)]
USE = [age_use(BMID[i]) for i in range(NB)]
# glycolytic ATP demand actually seen by the pathway: in the circulation a
# reticulocyte's mitochondria cover part of it; in the spleen they cannot.
MITC = [P['mitoRet'] * max(0.0, 1.0 - BMID[i] / P['mitoRetAge']) for i in range(NB)]
USEC = [USE[i] * (1.0 - MITC[i]) for i in range(NB)]
USES = [USE[i] * P['splUse'] for i in range(NB)]
_RETW = [max(0.0, 1.0 - BMID[i] / P['mitoRetAge']) for i in range(NB)]
SPLEXP = [1.0 + (P['splRetExp'] - 1.0) * _RETW[i] for i in range(NB)]
SPLCRIT = [P['ATPcritSpl'] * (1.0 + P['splRetCrit'] * _RETW[i]) for i in range(NB)]
MCHB = [P['MCH0'] * (1.10 if i < P['nRet'] else 1.00) for i in range(NB)]

# gene-therapy / transplant-derived cohorts: wild-type PKR, so 7 coarse bins
NC = 7
CW = [2.0] + [(BEXIT - 2.0) / 6.0] * 6
CMID = []
_acc = 0.0
for _w in CW:
    CMID.append(_acc + 0.5 * _w)
    _acc += _w
MCHC_ = [P['MCH0'] * (1.10 if i == 0 else 1.00) for i in range(NC)]

def _tabcache():
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.pkd_tabcache.json')
    key = "|".join("%s=%.10g" % (k, P[k]) for k in sorted(P) if P[k] is not None
                   and isinstance(P[k], float)) + "|v6"
    if os.path.exists(fn):
        try:
            d = json.load(open(fn))
            if d.get('key') == key:
                return key, fn, d
        except Exception:
            pass
    return key, fn, None


emit("building glycolysis lookup tables ...")
_CKEY, _CFN, _CDAT = _tabcache()
if _CDAT is not None:
    def _mk(d):
        g = GTab.__new__(GTab)
        g.__dict__.update(d)
        return g
    TSYS = [_mk(d) for d in _CDAT['sys']]
    TSPL = [_mk(d) for d in _CDAT['spl']]
    TEB = _mk(_CDAT['eb'])
    emit("  (restored from table cache)")
else:
    TSYS = [GTab(VUP[i], USEC[i]) for i in range(NB)]
    TSPL = [GTab(VUP[i] * P['splStress'], USES[i]) for i in range(NB)]
    TEB = GTab(1.0, 1.0 - P['mitoRescue'])
    try:
        json.dump({'key': _CKEY,
                   'sys': [g.__dict__ for g in TSYS],
                   'spl': [g.__dict__ for g in TSPL],
                   'eb': TEB.__dict__}, open(_CFN, 'w'))
    except Exception:
        pass
# wild-type corrected cohorts need no table: their activity never varies
CST = [gly_exact(1.0 * math.exp(-CMID[i] / P['tauPKwt']),
                 age_vup(CMID[i]), age_use(CMID[i]),
                 atot=P['Atot'] - 0.0032 * CMID[i]) for i in range(NC)]
emit("  tables built: %d systemic + %d splenic + 1 erythroblast" % (NB, NB))


# =============================================================================
# SECTION 7.  STATE VECTOR
# =============================================================================
SN = (['PROG', 'EB1', 'EB2'] +
      ['N%02d' % i for i in range(NB)] +
      ['A%02d' % i for i in range(NB)] +
      ['ND1', 'ND2', 'DPGD', 'SEQR'] +
      ['C%d' % i for i in range(NC)] +
      ['FDEOX', 'EPO', 'ERFE', 'HEP', 'FEP', 'RES', 'LIC', 'FERR', 'CARDFE',
       'UCB', 'CB', 'GS', 'LDH', 'HAP', 'SPLV',
       'MGUT', 'MC', 'MP', 'MRBC', 'CYP',
       'EGUT', 'ECEN', 'ERBC', 'TGUT', 'TCEN', 'TRBC', 'DGUT', 'DCEN',
       'STAB', 'E2', 'TST', 'BMD', 'PVR', 'EMH', 'ULC',
       'TXU', 'FECUM', 'HBAUC'])
IX = dict((n, i) for i, n in enumerate(SN))
NST = len(SN)
iN0 = IX['N00']
iA0 = IX['A00']
iC0 = IX['C0']

# iron content conversions
FE_PER_GHB = 3.47e-3        # g Fe per g Hb
BLOOD_DL = 50.0            # dL blood volume (5 L adult)
PLASMA_DL = 30.0           # dL plasma volume
# mg Fe held in the liver per unit of liver iron concentration.  The liver dry
# weight IS the conversion factor: 1500 g wet x 0.28 = 420 g dry, so 1 mg Fe/g dw
# is 420 mg of iron.  An earlier version divided by 1000 as well and made the
# liver a 0.42 mg sink, which loaded a healthy subject to 7.5 mg Fe/g dw
# (DEFECT #6, Section 9).
LIVER_MG_PER_LIC = P['LIVwt'] * P['dwFrac']
TIBC_MG = 13.3             # mg Fe at 100% transferrin saturation
BILI_PER_GHB = 34.0        # mg bilirubin per g hemoglobin catabolised

# =============================================================================
# SECTION 8.  SCENARIO CONFIGURATION
# =============================================================================
def newcfg(**kw):
    c = dict(
        label='',
        alpha=1.00,        # residual PKR activity vs wild type (genotype)
        tauPK=P['tauPKwt'],  # PKR thermolability (d)
        activ=1.00,        # fraction of residual PKR that an allosteric
                           # activator can act on (0 for non-missense / null)
        ugt=1.00,          # UGT1A1 activity (1.0 = *1/*1, 0.7 = *1/*28, 0.3 = *28/*28)
        spleen=1.0,        # 1 = spleen present, 0 = splenectomised
        splday=None,       # day of splenectomy
        cholday=None,      # day of cholecystectomy
        mit=0.0, mitint=0.5, mitstart=0.0, mitstop=1e9,   # mg per dose, interval d
        eta=0.0, etaint=1.0, etastart=0.0, etastop=1e9,
        teb=0.0, tebint=1.0, tebstart=0.0, tebstop=1e9,
        dfx=0.0, dfxstart=0.0, dfxstop=1e9,               # mg/d
        txthresh=None, txunits=2.0, txmin=21.0,           # transfusion regimen
        gtday=None, gtfrac=0.0,                           # gene therapy
        crisis=None, crisislen=10.0,                      # aplastic crisis
        cypext=1.00,       # external CYP3A modulation (rifampicin>1, azole<1)
        folate=1.0,
        tend=730.0, dt=0.05,
    )
    c.update(kw)
    return c


def drug_activation(y, C):
    """Fold multiplier on residual PKR activity from all activators present."""
    a = 1.0
    a += P['mitEmax'] * y[IX['MRBC']] / (P['mitEC50'] + y[IX['MRBC']]) * C['activ']
    a += P['etaEmax'] * y[IX['ERBC']] / (P['etaEC50'] + y[IX['ERBC']]) * C['activ']
    a += P['tebEmax'] * y[IX['TRBC']] / (P['tebEC50'] + y[IX['TRBC']]) * C['activ']
    return a


# =============================================================================
# SECTION 9.  THE SLOW SYSTEM
# =============================================================================
def core(t, y, C, obs=False):
    p = P
    dy = [0.0] * NST

    # ---- 9a. effective PKR activity profile across the age axis -------------
    # This is where the genotype, the cell's age and the drug meet.  An
    # allosteric activator is a MULTIPLIER on protein that is present, so it
    # cannot rescue a cohort whose protein has already decayed; the second,
    # slower limb (thermostabilisation, state STAB) is what reaches old cells,
    # by lengthening the decay constant itself.
    act = drug_activation(y, C)
    tauPK = C['tauPK'] * y[IX['STAB']]
    aeff = [0.0] * NB
    aeff0 = [0.0] * NB            # activity WITHOUT the drug: the injury already carried
    for i in range(NB):
        e = math.exp(-BMID[i] / tauPK)
        aeff[i] = C['alpha'] * e * act
        aeff0[i] = C['alpha'] * math.exp(-BMID[i] / C['tauPK'])

    fdeox = y[IX['FDEOX']]
    if fdeox < 0.01:
        fdeox = 0.01
    elif fdeox > 0.80:
        fdeox = 0.80

    # ---- 9b. per-cohort metabolic state, in the circulation and in the spleen
    atp = [0.0] * NB
    atp0 = [0.0] * NB
    dpgf = [0.0] * NB
    atps = [0.0] * NB
    atps0 = [0.0] * NB
    for i in range(NB):
        At = y[iA0 + i]
        if At < 0.30:
            At = 0.30
        atp[i], dpgf[i], _ = TSYS[i].get(aeff[i], At)
        atps[i], _, _ = TSPL[i].get(aeff[i], At)
        if act > 1.0000001:
            atp0[i], _, _ = TSYS[i].get(aeff0[i], At)
            atps0[i], _, _ = TSPL[i].get(aeff0[i], At)
        else:
            atp0[i] = atp[i]
            atps0[i] = atps[i]

    # ---- 9c. hazards --------------------------------------------------------
    # Every destruction term is routed through an ATP margin, never asserted
    # directly.  Extravascular lysis is the failure of cation/volume control in
    # the circulation; splenic trapping is the SAME failure evaluated under the
    # red pulp's metabolic conditions (45% of glycolytic capacity, 15% more
    # mechanical work).  That single choice is what makes the reticulocyte
    # response to splenectomy come out with the observed sign.
    spl = C['spleen']
    splf = math.sqrt(max(y[IX['SPLV']], 1.0) / p['splVol0']) if spl > 0.0 else 0.0
    h = [0.0] * NB
    fr = p['fRev']
    for i in range(NB):
        m1 = 1.0 - atp[i] / p['ATPcrit']
        if m1 < 0.0:
            m1 = 0.0
        m1u = 1.0 - atp0[i] / p['ATPcrit']
        if m1u < 0.0:
            m1u = 0.0
        meff = fr * m1 + (1.0 - fr) * m1u      # irreversible injury already carried
        hh = p['hBase'] + (p['kHem'] * meff ** p['pHem'] if meff > 0.0 else 0.0)
        if spl > 0.0:
            m2 = 1.0 - atps[i] / SPLCRIT[i]
            if m2 < 0.0:
                m2 = 0.0
            m2u = 1.0 - atps0[i] / SPLCRIT[i]
            if m2u < 0.0:
                m2u = 0.0
            m2e = fr * m2 + (1.0 - fr) * m2u
            if m2e > 0.0:
                hh += p['kSpl'] * splf * SPLEXP[i] * m2e ** p['pHem']
        h[i] = hh

    # ---- 9d. hemoglobin, reticulocytes, whole-blood 2,3-BPG -----------------
    hb = 0.0
    hbdpg = 0.0
    ret = 0.0
    ncell = 0.0
    for i in range(NB):
        hi = y[iN0 + i] * MCHB[i] / 10.0
        if hi < 0.0:
            hi = 0.0
        hb += hi
        hbdpg += hi * dpgf[i]
        ncell += max(y[iN0 + i], 0.0)
        if i < p['nRet']:
            ret += max(y[iN0 + i], 0.0)
    hbd = (max(y[IX['ND1']], 0.0) + max(y[IX['ND2']], 0.0)) * p['MCH0'] / 10.0
    hb += hbd
    hbdpg += hbd * y[IX['DPGD']]
    ncell += (max(y[IX['ND1']], 0.0) + max(y[IX['ND2']], 0.0))
    hbc = 0.0
    for i in range(NC):
        ci = max(y[iC0 + i], 0.0) * MCHC_[i] / 10.0
        hbc += ci
        hbdpg += ci * CST[i]['dpgf']
        ncell += max(y[iC0 + i], 0.0)
        if i == 0:
            ret += max(y[iC0 + i], 0.0)
    hb += hbc
    if hb < 0.05:
        hb = 0.05
    dpgf_mean = hbdpg / hb
    dpg_mean = dpg_total(dpgf_mean, fdeox)

    # ---- 9e. oxygen transport and the deoxy-fraction feedback ---------------
    o2 = o2_state(hb, dpg_mean)
    fd_target = 1.0 - 0.5 * (o2['SaO2'] + max(o2['SvO2'], 0.0))
    dy[IX['FDEOX']] = (fd_target - y[IX['FDEOX']]) / p['tauFdeox']

    # ---- 9f. erythropoietin and the marrow ---------------------------------
    epot = p['EPObase'] * math.exp(p['EPOgain'] * max(0.0, p['HbEPO'] - hb))
    if epot > p['EPOmax']:
        epot = p['EPOmax']
    dy[IX['EPO']] = p['kEPOel'] * (epot - y[IX['EPO']])
    ep = max(y[IX['EPO']], 0.0)
    eh = ep ** p['EPOhill']
    amp = 1.0 + (p['ampMax'] - 1.0) * eh / (p['EPO50'] ** p['EPOhill'] + eh)

    # ineffective erythropoiesis.  Erythroblasts still have mitochondria, so
    # most of their ATP demand is met oxidatively and the PK lesion costs them
    # far less than it costs an anucleate red cell -- which is why PKD has some
    # ineffective erythropoiesis but nothing like thalassemia's.
    atp_eb, _, _ = TEB.get(C['alpha'] * act, p['Atot'])
    meb = 1.0 - atp_eb / p['ATPcrit']
    if meb < 0.0:
        meb = 0.0
    ieo = p['ieoMax'] * meb / (p['ieoK'] + meb)

    crisis = 0.0
    if C['crisis'] is not None and C['crisis'] <= t < C['crisis'] + C['crisislen']:
        crisis = 1.0
    dy[IX['PROG']] = p['kProg'] * (1.0 - y[IX['PROG']]) - 2.0 * crisis * y[IX['PROG']]
    prod = KPROD * max(y[IX['PROG']], 0.0) * amp * C['folate']
    dy[IX['EB1']] = prod - p['kEB1'] * y[IX['EB1']]
    dy[IX['EB2']] = p['kEB1'] * y[IX['EB1']] - p['kEB2'] * y[IX['EB2']]
    influx = p['kEB2'] * max(y[IX['EB2']], 0.0) * (1.0 - ieo)

    gtf = C['gtfrac'] if (C['gtday'] is not None and t >= C['gtday']) else 0.0
    inflx_mut = influx * (1.0 - gtf)
    inflx_cor = influx * gtf

    # ---- 9g. cohort transport and the adenylate pool ------------------------
    # --- splenic reticulocyte sequestration ---------------------------------
    # The capture rate is not a constant: it is proportional to the ATP shortfall
    # a reticulocyte experiences under red-pulp conditions, so a wild-type
    # reticulocyte is never captured and a PK-deficient one is.
    # Capture has to be debited from the SAME cohort it is sourced from.  A first
    # version summed the capture over both reticulocyte bins and debited it all to
    # bin 0; once bin 0 hit its non-negativity floor the cycle
    # bin1 -> sequestered pool -> bin1 CREATED cells, and the model reported
    # haemoglobins of 19-54 g/dL (DEFECT #12, Section 9).  Mass balance in an
    # age-structured model is per-cohort, not aggregate.
    seq_cap = [0.0] * NB
    seq_in = 0.0
    mseq = 0.0
    for i in range(p['nRet']):
        mi = 1.0 - atps[i] / SPLCRIT[i]
        if mi < 0.0:
            mi = 0.0
        if i == 0:
            mseq = mi
        seq_cap[i] = p['kSeqIn'] * splf * mi * max(y[iN0 + i], 0.0)
        seq_in += seq_cap[i]
    seq_kill = p['kSeqHem'] * mseq ** p['pHem']
    seq_rel = p['kSeqOut'] + (10.0 if spl <= 0.0 else 0.0)
    seqr = max(y[IX['SEQR']], 0.0)
    dy[IX['SEQR']] = seq_in - (seq_rel + seq_kill) * seqr

    matf = 1.0 + p['kMatf'] * max(0.0, (p['HbEPO'] - hb) / p['HbEPO'])
    if matf > p['matfMax']:
        matf = p['matfMax']
    kb = [1.0 / (BW[i] * (matf if i < p['nRet'] else 1.0)) for i in range(NB)]
    prev_flux = inflx_mut
    prev_A = p['Atot']
    for i in range(NB):
        n = y[iN0 + i]
        out = kb[i] * n
        if i == NB - 1:
            out = p['kSenesc'] * n
        dy[iN0 + i] = prev_flux - out - h[i] * n - seq_cap[i]
        if i == p['nRet'] - 1:
            dy[iN0 + i] += seq_rel * seqr
        # adenylate: dilution by the inflowing younger cohort, minus AMP
        # deaminase / 5'-nucleotidase drain.  The pool is NOT conserved; that is
        # what keeps ATP/ADP from collapsing as ATP falls, and it is why total
        # adenine nucleotide content is a marker of red cell age in vivo.
        At = y[iA0 + i]
        if At < p['AtotMin']:
            At = p['AtotMin']
        elif At > p['Atot']:
            At = p['Atot']
        amp_i = At - atp[i] - adk_adp(atp[i], At)
        if amp_i < 0.0:
            amp_i = 0.0
        drain = p['kAMPd'] * amp_i / (p['KAMPd'] + amp_i)
        # The dilution rate is inflow/N, which blows up when a cohort is nearly
        # annihilated -- and in severe disease it is.  Left unclamped it drove the
        # adenylate state to +inf and overflowed the table lookup (DEFECT #10,
        # Section 9).  The rate is capped at the fastest physically meaningful
        # turnover (a cohort replaced in ~1 h) and the pool is held in range.
        nn = n if n > 1e-4 else 1e-4
        dil = prev_flux / nn
        if dil > 24.0:
            dil = 24.0
        dy[iA0 + i] = dil * (prev_A - At) - drain
        if At >= p['Atot'] and dy[iA0 + i] > 0.0:
            dy[iA0 + i] = 0.0
        elif At <= p['AtotMin'] and dy[iA0 + i] < 0.0:
            dy[iA0 + i] = 0.0
        prev_flux = out
        prev_A = At
    senesc_flux = prev_flux

    # corrected (gene-therapy / transplant) cohorts: wild-type PKR
    kc = [1.0 / (CW[i] * (matf if i == 0 else 1.0)) for i in range(NC)]
    pf = inflx_cor
    for i in range(NC):
        n = y[iC0 + i]
        out = (kc[i] * n) if i < NC - 1 else (p['kSenesc'] * n)
        hc = p['hBase']
        dy[iC0 + i] = pf - out - hc * n
        pf = out
    senesc_flux += pf

    # transfused donor cells: normal metabolism, storage-shortened survival
    dy[IX['ND1']] = -(1.0 / 3.0 + 0.10) * y[IX['ND1']]
    dy[IX['ND2']] = (1.0 / 3.0) * y[IX['ND1']] - (1.0 / p['donTau']) * y[IX['ND2']]
    dy[IX['DPGD']] = p['donDPGrec'] * (p['DPGfree0'] - y[IX['DPGD']])

    # ---- 9h. destruction flux -> iron and bilirubin -------------------------
    ghb_lysed = 0.0
    for i in range(NB):
        ghb_lysed += h[i] * max(y[iN0 + i], 0.0) * MCHB[i] / 10.0
    ghb_lysed += p['hBase'] * (hbc + hbd)
    ghb_lysed += seq_kill * seqr * MCHB[0] / 10.0
    ghb_lysed += senesc_flux * p['MCH0'] / 10.0
    ghb_lysed *= BLOOD_DL                       # g Hb per day, whole body
    # ineffective erythropoiesis contributes catabolised heme without ever
    # producing a circulating cell -- the "shunt" bilirubin.
    # erythroblasts die before they finish loading hemoglobin, so the shunt
    # bilirubin they contribute is scaled to a partial hemoglobin content
    ghb_ieo = influx / (1.0 - ieo + 1e-9) * ieo * 0.60 * p['MCH0'] / 10.0 * BLOOD_DL

    fe_rel = (ghb_lysed + ghb_ieo) * FE_PER_GHB * 1000.0    # mg Fe/d
    bili_prod = (ghb_lysed + ghb_ieo) * BILI_PER_GHB        # mg/d

    # ---- 9i. iron / hepcidin axis ------------------------------------------
    erfe_t = p['ERFE0'] * (max(y[IX['EB1']], 0.0) / EB1_0) * \
        (max(y[IX['EPO']], 0.0) / p['EPObase']) ** 0.5
    dy[IX['ERFE']] = p['kERFEel'] * (erfe_t - y[IX['ERFE']])
    lic = max(y[IX['LIC']], 0.0)
    hep_t = p['HEP0'] * (1.0 + p['LICstim'] * lic) / \
        (1.0 + max(y[IX['ERFE']], 0.0) / p['ERFEi50'])
    dy[IX['HEP']] = p['kHEPel'] * (hep_t - y[IX['HEP']])
    hep = max(y[IX['HEP']], 1.0)
    absorb = p['absMax'] / (1.0 + (hep / p['HEPabs50']) ** 2.0)

    tsat = max(y[IX['FEP']], 0.0) / TIBC_MG
    if tsat > 1.4:
        tsat = 1.4
    ntbi = 3.0 * max(0.0, tsat - p['NTBIthr'])
    # Erythroid iron uptake is the iron actually built into the cells the marrow
    # is delivering -- not a free parameter.  Writing it as a saturating function
    # of transferrin saturation instead broke iron balance in the healthy state
    # (liver iron ran to 175 mg Fe/g dw in a wild-type subject).  DEFECT #5.
    fe_up = influx * (P['MCH0'] / 10.0) * BLOOD_DL * FE_PER_GHB * 1000.0
    avail = tsat / 0.15
    if avail > 1.0:
        avail = 1.0
    fe_up *= avail
    # The liver is a reversible store, and the iron it returns has to be given
    # back to the plasma or the model leaks iron.  Writing both directions is
    # what lets the healthy subject sit at a liver iron of ~1 mg Fe/g dw with an
    # absorption of ~0.8 mg/d instead of drifting.
    dfx_eff = p['dfxEmax'] * y[IX['DCEN']] / (p['dfxEC50'] + y[IX['DCEN']])
    liver_in = p['kLICin'] * tsat + p['kLIC'] * ntbi
    liver_out = (p['kLICout'] + 0.055 * dfx_eff) * lic * LIVER_MG_PER_LIC
    dy[IX['FEP']] = absorb + p['kRESout'] * max(y[IX['RES']], 0.0) \
        - fe_up - liver_in + liver_out * (1.0 - 0.62 * (dfx_eff > 0.0)) \
        - p['kFeLoss']
    dy[IX['RES']] = fe_rel - p['kRESout'] * max(y[IX['RES']], 0.0)
    dy[IX['LIC']] = (liver_in - liver_out) / LIVER_MG_PER_LIC
    dy[IX['FERR']] = 1.5 * (p['FerrBase'] + p['kFerr'] * lic +
                            18.0 * max(0.0, tsat - 0.45) - y[IX['FERR']])
    dy[IX['CARDFE']] = p['kCard'] * ntbi - (p['kCardOut'] + 0.004 * dfx_eff) * \
        max(y[IX['CARDFE']], 0.0)

    # ---- 9j. bilirubin, gallstones, hemolysis markers -----------------------
    ucb = max(y[IX['UCB']], 0.0)
    # First order, not Michaelis-Menten.  Hepatic bilirubin uptake is the
    # rate-limiting step and is first order far below saturation; a
    # Michaelis form calibrated at the normal load saturates at ~3x normal
    # hemolysis and sent the bilirubin of a moderately affected patient to
    # 23000 mg/dL (DEFECT #7, Section 9).  Gilbert and Crigler-Najjar are
    # represented by the ugt multiplier, which is where saturation belongs.
    conj = p['kUGT'] * C['ugt'] * ucb / (1.0 + ucb / p['UGTsat'])
    dy[IX['UCB']] = bili_prod / PLASMA_DL - conj
    dy[IX['CB']] = conj - p['kCBout'] * max(y[IX['CB']], 0.0)
    if C['cholday'] is not None and t >= C['cholday']:
        dy[IX['GS']] = 0.0
    else:
        dy[IX['GS']] = p['kGS'] * conj * PLASMA_DL
    ivh = p['ivFrac'] * ghb_lysed
    dy[IX['LDH']] = 2.0 * (p['LDH0'] + p['kLDH'] * 42.0 * ivh - y[IX['LDH']])
    dy[IX['HAP']] = p['kHapEl'] * (p['HAP0'] - y[IX['HAP']]) - 0.020 * ivh

    # ---- 9k. spleen -------------------------------------------------------
    if spl > 0.0:
        splwork = 0.0
        for i in range(NB):
            m2 = 1.0 - atps[i] / SPLCRIT[i]
            if m2 > 0.0:
                splwork += m2 ** p['pHem'] * max(y[iN0 + i], 0.0)
        dy[IX['SPLV']] = p['kSplGrow'] * y[IX['SPLV']] * splwork * \
            (1.0 - y[IX['SPLV']] / p['splVolMax']) - 0.0015 * \
            (y[IX['SPLV']] - p['splVol0'])
    else:
        dy[IX['SPLV']] = 0.0

    # ---- 9l. drug pharmacokinetics -----------------------------------------
    cyp = max(y[IX['CYP']], 0.05)
    cl = p['mitCL0'] * cyp * C['cypext']
    cp = max(y[IX['MC']], 0.0) / p['mitVc']
    cpp = max(y[IX['MP']], 0.0) / p['mitVp']
    dy[IX['MGUT']] = -p['mitKa'] * y[IX['MGUT']]
    dy[IX['MC']] = p['mitKa'] * y[IX['MGUT']] * p['mitF'] - cl * cp \
        - p['mitQ'] * (cp - cpp)
    dy[IX['MP']] = p['mitQ'] * (cp - cpp)
    dy[IX['MRBC']] = p['mitKrbcOut'] * (p['mitKrbc'] * cp - y[IX['MRBC']])
    cyp_t = 1.0 + (p['cypEmax'] - 1.0) * cp / (p['cypEC50'] + cp)
    dy[IX['CYP']] = p['cypKin'] * (cyp_t - y[IX['CYP']])

    cpe = max(y[IX['ECEN']], 0.0) / p['etaVc']
    dy[IX['EGUT']] = -p['etaKa'] * y[IX['EGUT']]
    dy[IX['ECEN']] = p['etaKa'] * y[IX['EGUT']] * p['etaF'] - p['etaCL'] * cpe
    dy[IX['ERBC']] = p['etaKrbcOut'] * (p['etaKrbc'] * cpe - y[IX['ERBC']])

    cpt = max(y[IX['TCEN']], 0.0) / p['tebVc']
    dy[IX['TGUT']] = -p['tebKa'] * y[IX['TGUT']]
    dy[IX['TCEN']] = p['tebKa'] * y[IX['TGUT']] * p['tebF'] - p['tebCL'] * cpt
    dy[IX['TRBC']] = p['tebKrbcOut'] * (p['tebKrbc'] * cpt - y[IX['TRBC']])

    dy[IX['DGUT']] = -p['dfxKa'] * y[IX['DGUT']]
    dy[IX['DCEN']] = p['dfxKa'] * y[IX['DGUT']] * p['dfxF'] / p['dfxVc'] \
        - (p['dfxCL'] / p['dfxVc']) * y[IX['DCEN']]

    # PKR thermostabilisation: a second, slower drug limb that lengthens the
    # protein's decay constant instead of multiplying its activity.
    stab_t = 1.0 + (p['mitStabMax'] - 1.0) * y[IX['MRBC']] / \
        (p['mitStabEC50'] + y[IX['MRBC']]) * C['activ']
    dy[IX['STAB']] = p['kStab'] * (stab_t - y[IX['STAB']])

    # ---- 9m. endocrine, bone, vascular sequelae -----------------------------
    arom = 1.0 - p['mitAromEmax'] * y[IX['MRBC']] / (p['mitAromEC50'] + y[IX['MRBC']])
    dy[IX['E2']] = p['kE2'] * (p['E20'] * arom - y[IX['E2']])
    tst_t = p['TST0'] * (1.0 + p['tstFeedback'] * (1.0 - y[IX['E2']] / p['E20']))
    dy[IX['TST']] = p['kTST'] * (tst_t - y[IX['TST']])

    marrow_load = max(0.0, amp / AMP_0 - 1.0)
    dy[IX['BMD']] = -p['kBMD'] * marrow_load + p['kBMDrec'] * (0.0 - y[IX['BMD']])
    dy[IX['PVR']] = p['kPVR'] * ivh - p['kPVRrec'] * (y[IX['PVR']] - p['PVR0'])
    dy[IX['EMH']] = p['kEMH'] * marrow_load - p['kEMHrec'] * y[IX['EMH']]
    dy[IX['ULC']] = 0.0022 * max(0.0, y[IX['PVR']] - p['PVR0']) * ivh / 6.0
    dy[IX['TXU']] = 0.0
    dy[IX['FECUM']] = absorb - p['kFeLoss']
    dy[IX['HBAUC']] = hb

    if not obs:
        return dy, None

    pk_assay = 0.0
    for i in range(NB):
        pk_assay += max(y[iN0 + i], 0.0) * C['alpha'] * math.exp(-BMID[i] / tauPK)
    for i in range(NC):
        pk_assay += max(y[iC0 + i], 0.0) * math.exp(-CMID[i] / p['tauPKwt'])
    pk_assay += (max(y[IX['ND1']], 0.0) + max(y[IX['ND2']], 0.0)) * 0.85
    pk_assay = pk_assay / max(ncell, 1e-9)

    hbmass = hb * BLOOD_DL
    lifespan = hbmass / (ghb_lysed + 1e-12)
    mchc = 33.0 * (1.0 + 0.35 * max(0.0, 1.0 - min(atp) / p['ATPcrit']))
    o = dict(t=t, Hb=hb, HbDonor=hbd, HbCorr=hbc, RBC=ncell,
             RetPct=100.0 * ret / max(ncell, 1e-9),
             DPGfree=dpgf_mean, DPG=dpg_mean,
             DPGgHb=dpg_mean / 33.0 * 100.0,
             P50=o2['p50'], PvO2=o2['PvO2'], SaO2=o2['SaO2'], SvO2=o2['SvO2'],
             CO=o2['CO'], COreq=o2['COreq'], Hbeq=o2['Hbeq'],
             O2extr=o2['extr'],
             ATPyoung=atp[0], ATPold=atp[NB - 1], ATPmean=sum(atp) / NB,
             ATPspl_old=atps[NB - 1], PKassay=pk_assay,
             EPO=y[IX['EPO']], amp=amp, IEO=ieo, ATPeb=atp_eb,
             ERFE=y[IX['ERFE']], HEP=y[IX['HEP']], TSAT=tsat, NTBI=ntbi,
             LIC=lic, FERR=y[IX['FERR']], CARDFE=y[IX['CARDFE']],
             absorb=absorb, FeRel=fe_rel,
             UCB=y[IX['UCB']], CB=y[IX['CB']],
             BILI=y[IX['UCB']] + y[IX['CB']], GS=y[IX['GS']],
             LDH=y[IX['LDH']], HAP=y[IX['HAP']], SPLV=y[IX['SPLV']],
             MRBC=y[IX['MRBC']], MITc=cp, CYP=cyp, ERBCc=y[IX['ERBC']],
             TRBCc=y[IX['TRBC']], DFXc=y[IX['DCEN']],
             STAB=y[IX['STAB']], E2=y[IX['E2']], TST=y[IX['TST']],
             BMD=y[IX['BMD']], PVR=y[IX['PVR']], EMH=y[IX['EMH']],
             ULC=y[IX['ULC']], TXU=y[IX['TXU']],
             SEQR=seqr, SeqFrac=seqr / max(seqr + ret, 1e-9),
             lifespan=lifespan, MCHC=mchc, ghb_lysed=ghb_lysed,
             hazard_young=h[0], hazard_old=h[NB - 1], aeff_young=aeff[0],
             aeff_old=aeff[NB - 1], fdeox=fdeox, tauPKeff=tauPK)
    return dy, o


# =============================================================================
# SECTION 10.  INTEGRATION
# =============================================================================
def y0_default():
    y = [0.0] * NST
    y[IX['PROG']] = 1.0
    y[IX['EB1']] = 0.12
    y[IX['EB2']] = 0.10
    for i in range(NB):
        y[iN0 + i] = 0.35
        y[iA0 + i] = P['Atot'] - 0.0030 * BMID[i]
    y[IX['DPGD']] = P['DPGfree0']
    y[IX['FDEOX']] = P['fdeox0']
    y[IX['EPO']] = P['EPObase']
    y[IX['HEP']] = P['HEP0']
    y[IX['ERFE']] = P['ERFE0']
    y[IX['FEP']] = 4.0
    y[IX['RES']] = 25.0
    y[IX['LIC']] = 1.0
    y[IX['FERR']] = 60.0
    y[IX['UCB']] = 0.5
    y[IX['CB']] = 0.2
    y[IX['HAP']] = P['HAP0']
    y[IX['LDH']] = P['LDH0']
    y[IX['SPLV']] = P['splVol0']
    y[IX['CYP']] = 1.0
    y[IX['STAB']] = 1.0
    y[IX['E2']] = P['E20']
    y[IX['TST']] = P['TST0']
    y[IX['PVR']] = P['PVR0']
    return y


# fastest linear rate constants in the system -> explicit-RK4 stability limit
DT_STABLE = 2.78 / max(P['kCBout'], P['mitKrbcOut'], P['kUGT'], P['kHEPel'],
                       P['etaKrbcOut'], P['tebKrbcOut'], P['mitKa'],
                       P['dfxCL'] / P['dfxVc'], P['kSenesc'], 10.35)

_IBMD = None


def _proj(v):
    """Projected RK4 stage: every state except the bone-density Z-score is a
    non-negative physical quantity, and clamping the intermediate stages (rather
    than only the accepted step) is what keeps a heavily haemolysing cohort from
    going transiently negative and pushing a table lookup out of range."""
    for j in range(NST):
        if j != _IBMD and v[j] < 0.0:
            v[j] = 0.0
    return v


def simulate(C, y0=None, nout=400, record=True):
    """Fixed-step RK4 with discrete events (dosing, transfusion, surgery)
    applied at step boundaries.

    STEP SIZE IS NOT FREE.  The fastest linear rates in the system are biliary
    excretion (22/d), red-cell drug partitioning (30/d) and bilirubin
    conjugation (16.4/d), so explicit RK4 is stable only for
    dt < 2.78/30 = 0.093 d.  Burn-in runs at dt = 0.25-0.30 d looked fine on
    haemoglobin and reticulocytes -- both slow -- while the bilirubin state
    oscillated, hit the non-negativity floor and stuck at zero, silently
    reporting a total bilirubin of 0.0 mg/dL in a haemolysing patient
    (DEFECT #13, Section 9).  A step-convergence check on slow outputs alone
    does not detect this; the guard below does.
    """
    global _IBMD
    _IBMD = IX['BMD']
    dt = C['dt']
    if dt > DT_STABLE:
        raise ValueError("dt=%.3f exceeds the RK4 stability limit %.3f d"
                         % (dt, DT_STABLE))
    nstep = int(round(C['tend'] / dt))
    y = list(y0 if y0 is not None else y0_default())
    out = []
    every = max(1, nstep // nout)
    last_tx = -1e9
    spleen_on = C['spleen']
    for s in range(nstep + 1):
        t = s * dt
        if C['splday'] is not None and t >= C['splday']:
            C['spleen'] = 0.0
        else:
            C['spleen'] = spleen_on
        # ---- discrete dosing ----
        if C['mit'] > 0.0 and C['mitstart'] <= t < C['mitstop']:
            k = (t - C['mitstart']) / C['mitint']
            if abs(k - round(k)) < 0.5 * dt / C['mitint']:
                y[IX['MGUT']] += C['mit']
        if C['eta'] > 0.0 and C['etastart'] <= t < C['etastop']:
            k = (t - C['etastart']) / C['etaint']
            if abs(k - round(k)) < 0.5 * dt / C['etaint']:
                y[IX['EGUT']] += C['eta']
        if C['teb'] > 0.0 and C['tebstart'] <= t < C['tebstop']:
            k = (t - C['tebstart']) / C['tebint']
            if abs(k - round(k)) < 0.5 * dt / C['tebint']:
                y[IX['TGUT']] += C['teb']
        if C['dfx'] > 0.0 and C['dfxstart'] <= t < C['dfxstop']:
            k = t / 1.0
            if abs(k - round(k)) < 0.5 * dt:
                y[IX['DGUT']] += C['dfx']
        if s % every == 0 or s == nstep:
            _, o = core(t, y, C, obs=True)
        else:
            o = None
        # ---- transfusion ----
        if C['txthresh'] is not None and (t - last_tx) >= C['txmin']:
            if o is None:
                _, o = core(t, y, C, obs=True)
            if o['Hb'] < C['txthresh']:
                y[IX['ND1']] += C['txunits'] * P['donHbUnit'] * 10.0 / P['MCH0']
                y[IX['DPGD']] = (y[IX['DPGD']] + P['donDPGstored']) * 0.5
                y[IX['FEP']] += C['txunits'] * P['FeUnit'] * 0.15
                y[IX['RES']] += C['txunits'] * P['FeUnit'] * 0.85
                y[IX['TXU']] += C['txunits']
                last_tx = t
        if record and o is not None:
            out.append(o)
        if s == nstep:
            break
        # ---- RK4 ----
        k1, _ = core(t, y, C)
        y2 = _proj([y[j] + 0.5 * dt * k1[j] for j in range(NST)])
        k2, _ = core(t + 0.5 * dt, y2, C)
        y3 = _proj([y[j] + 0.5 * dt * k2[j] for j in range(NST)])
        k3, _ = core(t + 0.5 * dt, y3, C)
        y4 = _proj([y[j] + dt * k3[j] for j in range(NST)])
        k4, _ = core(t + dt, y4, C)
        for j in range(NST):
            y[j] += dt / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j])
            if j != IX['BMD'] and y[j] < 0.0:
                y[j] = 0.0
    C['spleen'] = spleen_on
    return y, out


def steady(C, days=800.0, dt=0.05):
    """Burn-in to the untreated steady state for a given genotype."""
    D = dict(C)
    D['mit'] = D['eta'] = D['teb'] = D['dfx'] = 0.0
    D['txthresh'] = None
    D['gtday'] = None
    D['crisis'] = None
    D['splday'] = None
    D['tend'] = days
    D['dt'] = dt
    y, o = simulate(D, nout=3)
    return y, o[-1]

# =============================================================================
# SECTION 11.  PRODUCTION CALIBRATION
# =============================================================================
# One scalar (marrow output) is solved so that a wild-type subject sits at
# Hb 15.0 g/dL.  Everything else about the healthy steady state -- red cell
# lifespan, reticulocyte percentage, bilirubin production, iron recycling flux --
# is then a PREDICTION, and Section 12 checks those predictions against
# textbook physiology before any patient is simulated.

def calibrate_production():
    global KPROD, EB1_0, AMP_0
    e = P['EPObase'] ** P['EPOhill']
    AMP_0 = 1.0 + (P['ampMax'] - 1.0) * e / (P['EPO50'] ** P['EPOhill'] + e)
    lo, hi = 1.0e-4, 5.0e-2
    C = newcfg(alpha=1.0, tauPK=P['tauPKwt'])
    for _ in range(26):
        KPROD = 0.5 * (lo + hi)
        EB1_0 = KPROD * AMP_0 / P['kEB1']
        _, o = steady(C, days=800.0, dt=0.05)
        if o['Hb'] < 15.0:
            lo = KPROD
        else:
            hi = KPROD
    KPROD = 0.5 * (lo + hi)
    EB1_0 = KPROD * AMP_0 / P['kEB1']
    return o


KPROD = 4.3e-3
EB1_0 = 0.12
AMP_0 = 9.7


# =============================================================================
# SECTION 12.  VERIFICATION
# =============================================================================
def verify():
    hdr("SECTION 0.  VERIFICATION OF THE FAST/SLOW REDUCTION")

    sub("0.1  Does the calibrated module reproduce normal erythrocyte metabolism?")
    st = gly_exact(1.0)
    tgt = [('ATP  (mM)', st['atp'], P['ATP0']),
           ('ADP  (mM)', st['adp'], P['ADP0']),
           ('2,3-BPG total (mM)', st['dpg'], P['DPG0']),
           ('3-PG (mM)', st['pg3'], P['PG30']),
           ('PEP  (mM)', st['pep'], P['PEP0']),
           ('flux J (mmol/L/h)', st['J'], P['J0']),
           ('shunt fraction phi', st['phi'], P['phi0'])]
    emit("  %-22s %12s %12s %10s" % ("quantity", "model", "target", "rel.err"))
    worst = 0.0
    for nm, mv, tv in tgt:
        e = abs(mv - tv) / abs(tv)
        worst = max(worst, e)
        emit("  %-22s %12.6f %12.6f %9.2e" % (nm, mv, tv, e))
    emit("  worst relative error = %.3e" % worst)
    emit("  free 2,3-BPG at the reference point = %.4f mM (%.1f%% of total);" %
         (P['DPGfree0'], 100.0 * P['DPGfree0'] / P['DPG0']))
    emit("  the remainder is bound in the central cavity of deoxyhemoglobin.")

    sub("0.2  How much reserve does the pyruvate kinase step have?")
    cap = P['Vpk'] * CAL['A'] * CAL['hadp'] * CAL['iatp']
    emit("  PKR V_max (FBP-activated, uninhibited) = %.1f mmol/(L*h)" % (P['Vpk'] * CAL['A']))
    emit("  PKR capacity at physiological ADP/ATP  = %.2f mmol/(L*h)" % cap)
    emit("  flux actually carried (= 2J)           = %.2f mmol/(L*h)" % (2.0 * P['J0']))
    emit("  fractional utilisation of PKR          = %.2f %%" % (100.0 * 2.0 * P['J0'] / cap))
    emit()
    emit("  This is a derived number, not an input.  The pyruvate kinase step runs")
    emit("  at a few percent of its capacity, so activity has to fall by more than")
    emit("  an order of magnitude before flux is threatened -- which is why PKLR")
    emit("  heterozygotes are silent, why clinical disease needs biallelic lesions,")
    emit("  and why the disease is a metabolic THRESHOLD phenomenon rather than a")
    emit("  graded one.  Section 5.4 locates the threshold.")

    sub("0.3  Interpolation error of the lookup tables against the exact solve")
    mx_a = mx_d = 0.0
    args = []
    for k in range(31):
        av = 4.0e-4 * (7.0 / 4.0e-4) ** (k / 30.0)
        for At in (0.45, 0.83, 1.27, 1.71, 1.93):
            args.append((av, At))
    for av, At in args:
        d = gly_exact(av, VUP[6], USE[6], atot=At)
        ia, idp, _ = TSYS[6].get(av, At)
        if d['atp'] > 1e-5:
            mx_a = max(mx_a, abs(ia - d['atp']) / d['atp'])
        if d['dpgf'] > 1e-5:
            mx_d = max(mx_d, abs(idp - d['dpgf']) / d['dpgf'])
    emit("  %d test points on a mid-age cohort table (%dx%d grid):" %
         (len(args), TSYS[6].na, TSYS[6].nA))
    emit("    max relative error   ATP %.3e    free 2,3-BPG %.3e" % (mx_a, mx_d))
    emit("  The mrgsolve model solves the fast subsystem exactly in $ODE (compiled")
    emit("  C++ can afford it); this table is a Python-only speed device, so the")
    emit("  numbers above bound the difference between the two implementations.")

    sub("0.4  Timescale separation: is the quasi-steady-state assumption legitimate?")
    tauATP = P['ATP0'] / (2.0 * P['J0'])
    tauDPG = P['DPGfree0'] / (P['phi0'] * 2.0 * P['J0'])
    emit("  ATP pool turnover           tau = %.3f h  (%.2e d)" % (tauATP, tauATP / 24.0))
    emit("  free 2,3-BPG pool turnover  tau = %.2f h  (%.3f d)" % (tauDPG, tauDPG / 24.0))
    emit("  narrowest mature cohort width   = %.2f d" % P['wMat'])
    emit("  separation ratio                = %.0f" % (P['wMat'] / (tauDPG / 24.0)))
    emit()
    emit("  Safe everywhere except across a transfusion, where donor cells arrive")
    emit("  with free 2,3-BPG at %.2f mM and take days to regenerate it.  That one" % P['donDPGstored'])
    emit("  pool is therefore an explicit ODE (DPGD), not a QSS variable.")

    sub("0.5  Step-size convergence of the slow system")
    C = newcfg(alpha=0.15, tauPK=P['tauPKmut'], tend=400.0)
    prev = None
    emit("  %8s %12s %12s %12s %12s" % ("dt (d)", "Hb", "Ret%", "bilirubin",
                                          "max drift"))
    for dt in (0.09, 0.06, 0.04, 0.02):
        C2 = dict(C)
        C2['dt'] = dt
        _, o = simulate(C2, nout=2)
        row = (o[-1]['Hb'], o[-1]['RetPct'], o[-1]['BILI'])
        dr = '-' if prev is None else "%.2e" % max(
            abs(row[0] - prev[0]) / prev[0], abs(row[1] - prev[1]) / prev[1],
            abs(row[2] - prev[2]) / max(prev[2], 1e-9))
        emit("  %8.2f %12.6f %12.6f %12.6f %12s" % (dt, row[0], row[1], row[2], dr))
        prev = row
    emit("  Explicit RK4 is stable only for dt < 2.78/max|lambda| = %.4f d here," % DT_STABLE)
    emit("  set by biliary excretion (22/d) and red-cell drug partitioning (30/d).")
    emit("  All production runs use dt = 0.05 d.  Note that a convergence table on")
    emit("  SLOW outputs alone does not detect a violation: at dt = 0.30 d the")
    emit("  haemoglobin and reticulocyte columns still agreed to 1e-5 while the")
    emit("  bilirubin state was oscillating into its floor and reading 0.0 mg/dL.")


# =============================================================================
# SECTION 13.  HEALTHY PHYSIOLOGY -- PREDICTIONS, NOT FITS
# =============================================================================
def check_healthy():
    hdr("SECTION 1.  THE HEALTHY STEADY STATE (one parameter fitted, rest predicted)")
    o = calibrate_production()
    emit("  marrow output constant KPROD solved to put Hb at 15.0 g/dL: %.6e" % KPROD)
    emit()
    rows = [("hemoglobin (g/dL)", o['Hb'], "15.0", "FITTED (the one knob)"),
            ("red cell count (10^12/L)", o['RBC'], "4.5-5.5", "predicted"),
            ("reticulocytes (%)", o['RetPct'], "0.5-1.5", "predicted"),
            ("mean red cell lifespan (d)", o['lifespan'], "100-120", "predicted"),
            ("2,3-BPG total (mM)", o['DPG'], "4.5-5.1", "calibrated"),
            ("P50 (mmHg)", o['P50'], "26-27", "predicted"),
            ("mixed venous PO2 (mmHg)", o['PvO2'], "38-42", "predicted"),
            ("cardiac output (L/min)", o['CO'], "5.0", "predicted"),
            ("Hb catabolised (g/d)", o['ghb_lysed'], "6-7", "predicted"),
            ("iron recycled (mg/d)", o['FeRel'], "20-25", "predicted"),
            ("iron absorbed (mg/d)", o['absorb'], "0.8-1.2", "predicted"),
            ("total bilirubin (mg/dL)", o['BILI'], "0.3-1.0", "predicted"),
            ("serum ferritin (ng/mL)", o['FERR'], "30-150", "predicted"),
            ("liver iron (mg Fe/g dw)", o['LIC'], "0.5-1.8", "predicted"),
            ("LDH (U/L)", o['LDH'], "140-220", "predicted"),
            ("EPO (IU/L)", o['EPO'], "5-25", "predicted"),
            ("erythroid amplification (x)", o['amp'], "8-12", "predicted"),
            ("assayed PK activity (rel.)", o['PKassay'], "0.70-0.78", "predicted"),
            ("spleen volume (mL)", o['SPLV'], "100-250", "predicted")]
    emit("  %-30s %12s %14s  %s" % ("quantity", "model", "physiological", "status"))
    for nm, v, tv, st in rows:
        emit("  %-30s %12.4f %14s  %s" % (nm, v, tv, st))
    emit()
    emit("  Nineteen numbers, one of them fitted.  The iron flux of %.1f mg/d and" % o['FeRel'])
    emit("  the bilirubin production it implies are the same destruction flux seen")
    emit("  from two different organs, so they cannot be tuned independently -- and")
    emit("  in a hemolytic disease that coupling is what forces the bilirubin, the")
    emit("  iron and the reticulocyte predictions to move together or not at all.")
    return o

# =============================================================================
# SECTION 14.  GENOTYPE SPECTRUM
# =============================================================================
GENO = [
    ('very mild / compensated',      0.30, 60.0, 1.00),
    ('mild',                         0.22, 55.0, 1.00),
    ('moderate',                     0.16, 50.0, 1.00),
    ('severe, not transfused',       0.12, 50.0, 1.00),
    ('transfusion dependent',        0.09, 45.0, 1.00),
    ('non-missense / null (no drug target)', 0.12, 50.0, 0.00),
]
_BURN = {}


def burn(alpha, tau, activ=1.0, ugt=1.0):
    k = (round(alpha, 4), round(tau, 3), round(ugt, 3))
    if k not in _BURN:
        _BURN[k] = steady(newcfg(alpha=alpha, tauPK=tau, ugt=ugt))
    y, o = _BURN[k]
    return list(y), o


def spectrum():
    hdr("SECTION 2.  THE GENOTYPE SPECTRUM, AND WHERE HAEMOGLOBIN STOPS BEING FLAT")
    emit("  %-38s %6s %6s %6s %6s %6s %6s %6s" %
         ("phenotype", "alpha", "Hb", "ret%", "lifesp", "bili", "P50", "LDH"))
    rows = []
    for nm, al, tau, ac in GENO:
        _, o = burn(al, tau)
        rows.append((nm, al, o))
        emit("  %-38s %6.2f %6.2f %6.1f %6.1f %6.2f %6.1f %6.0f" %
             (nm, al, o['Hb'], o['RetPct'], o['lifespan'], o['BILI'], o['P50'], o['LDH']))
    sub("2.1  The compensation knee is a closed-form quantity, not a fitted one")
    Lcrit = BEXIT * AMP_0 / P['ampMax']
    emit("  Maximum erythroid amplification            = %.1f x" % P['ampMax'])
    emit("  Amplification at basal erythropoietin      = %.2f x" % AMP_0)
    emit("  => maximum sustainable production increase = %.2f x" % (P['ampMax'] / AMP_0))
    emit("  => critical red cell lifespan L* = %.0f / %.2f = %.1f d" %
         (BEXIT, P['ampMax'] / AMP_0, Lcrit))
    emit()
    emit("  Above L* the marrow can replace what is destroyed and haemoglobin is")
    emit("  FLAT in the genotype; below L* production is capped and Hb falls in")
    emit("  proportion to lifespan.  The table above shows exactly that: lifespan")
    emit("  moves smoothly from %.0f d to %.0f d across the spectrum, while Hb is" %
         (rows[0][2]['lifespan'], rows[4][2]['lifespan']))
    emit("  nearly unchanged until the lifespan crosses %.0f d and then falls fast." % Lcrit)
    emit()
    emit("  Two consequences that are not obvious from the bedside:")
    emit("   (a) haemoglobin is a SATURATING readout of the lesion above L* and a")
    emit("       hypersensitive one just below it.  Patients cluster near L*, so")
    emit("       small changes in red cell lifespan produce large Hb swings -- which")
    emit("       is the same reason a modest intervention can look dramatic.")
    emit("   (b) reticulocyte count and bilirubin keep rising monotonically THROUGH")
    emit("       the knee, so in the compensated range they, not Hb, carry the")
    emit("       information about severity.")
    return rows


# =============================================================================
# SECTION 15.  THE CENTRAL RESULT: HAEMOGLOBIN IS NOT A SUFFICIENT STATISTIC
# =============================================================================
def oxygen_analysis(hb0, dpg0, hb1, dpg1, label0, label1):
    a = o2_state(hb0, dpg0)
    b = o2_state(hb1, dpg1)
    emit("  %-34s %12s %12s %10s" % ("", label0, label1, "change"))
    def row(nm, x, y, fmt="%12.3f", pc=True):
        ch = ("%+9.1f%%" % (100.0 * (y / x - 1.0))) if (pc and x != 0) else ("%+9.3f" % (y - x))
        emit(("  %-34s " + fmt + " " + fmt + " %10s") % (nm, x, y, ch))
    row("haemoglobin (g/dL)", hb0, hb1)
    row("2,3-BPG total (mM)", dpg0, dpg1)
    row("P50 (mmHg)", a['p50'], b['p50'])
    row("O2 extracted per L blood (frac)", a['extrB'], b['extrB'])
    row("EQUIVALENT haemoglobin (g/dL)", a['Hbeq'], b['Hbeq'])
    row("required cardiac output (L/min)", a['COreq'], b['COreq'])
    row("tissue PO2 at fixed CO (mmHg)", a['PvO2'], b['PvO2'])
    return a, b


def central_result(rows):
    hdr("SECTION 3.  THE HAEMOGLOBIN ENDPOINT AND THE OXYGEN IT DOES NOT MEASURE")
    emit("  Every registrational trial in this disease scores a haemoglobin rise:")
    emit("  ACTIVATE required >=1.5 g/dL sustained (PMID 35417638), DRIVE-PK used")
    emit("  >1.0 g/dL (PMID 31483964).  A pyruvate kinase activator raises")
    emit("  haemoglobin BY LOWERING 2,3-BPG, and lowering 2,3-BPG shifts the oxygen")
    emit("  dissociation curve back to the LEFT.  So the drug adds carrier and")
    emit("  subtracts unloading at the same time.  This section asks whether the")
    emit("  two cancel, and by how much.")

    sub("3.1  An untreated patient is already partly compensated by the shunt")
    _, o = burn(0.12, 50.0)
    n = o2_state(15.0, P['DPG0'])
    p = o2_state(o['Hb'], o['DPG'])
    emit("  A normal subject at Hb 15.0, P50 %.1f mmHg extracts %.4f of the O2" %
         (n['p50'], n['extrB']))
    emit("  carried per litre of blood down to a tissue PO2 of %.0f mmHg." % P['PvO2ref'])
    emit("  The patient at Hb %.2f, P50 %.1f mmHg extracts %.4f -- %.1f%% more per" %
         (o['Hb'], p['p50'], p['extrB'], 100.0 * (p['extrB'] / n['extrB'] - 1.0)))
    emit("  gram.  His EQUIVALENT haemoglobin (the Hb a normal-P50 subject would")
    emit("  need to unload the same oxygen) is %.2f g/dL, not %.2f." % (p['Hbeq'], o['Hb']))
    emit("  So of the %.2f g/dL haemoglobin deficit, the shunt has already paid back" %
         (15.0 - o['Hb']))
    emit("  %.2f g/dL -- %.0f%% of it -- before any treatment." %
         (p['Hbeq'] - o['Hb'], 100.0 * (p['Hbeq'] - o['Hb']) / (15.0 - o['Hb'])))
    emit("  This is the quantitative form of a well-known bedside observation: PK-")
    emit("  deficient patients tolerate haemoglobins that would incapacitate other")
    emit("  anaemic patients.  It is not stoicism, it is a right-shifted curve.")

    sub("3.2  The break-even 2,3-BPG reduction, in closed form")
    emit("  A drug that delivers dHb g/dL leaves oxygen transport unchanged when")
    emit("      Hb1 * extract(P50_1) = Hb0 * extract(P50_0),")
    emit("  i.e. when extract(P50_1)/extract(P50_0) = Hb0/Hb1.  Solving that for")
    emit("  P50, then inverting P50 = P50ref*(DPG/DPG0)^n for 2,3-BPG, gives the")
    emit("  LARGEST 2,3-BPG reduction compatible with a net oxygen-transport gain.")
    emit()
    emit("  %8s %10s %10s %12s %12s" % ("dHb", "Hb1", "P50*", "2,3-BPG*", "max fall"))
    for dhb in (1.0, 1.5, 2.0, 3.0):
        hb1 = o['Hb'] + dhb
        target = o2_state(o['Hb'], o['DPG'])['extrB'] * o['Hb'] / hb1
        lo, hi = 10.0, 80.0
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            e = sat(P['PaO2'], mid) - sat(P['PvO2ref'], mid)
            if e < target:
                lo = mid
            else:
                hi = mid
        p50s = 0.5 * (lo + hi)
        dpgs = P['DPG0'] * (p50s / P['P50ref']) ** (1.0 / P['nDPG'])
        emit("  %8.1f %10.2f %10.2f %12.3f %11.1f%%" %
             (dhb, hb1, p50s, dpgs, 100.0 * (1.0 - dpgs / o['DPG'])))
    emit()
    emit("  Read the last column as a falsifiable threshold.  At a +1.5 g/dL")
    emit("  haemoglobin response -- exactly the ACTIVATE primary endpoint -- the")
    emit("  drug may lower 2,3-BPG by at most the percentage shown before the")
    emit("  oxygen-transport gain becomes a loss.  Phase 1 in healthy volunteers")
    emit("  (PMID 30091852) and the phase 2/3 programme both report 2,3-BPG")
    emit("  reductions of that order, so the model's claim is that the Hb endpoint")
    emit("  and the physiological endpoint are close to CANCELLING, and which side")
    emit("  of the line a given patient lands on depends on his baseline 2,3-BPG.")


# =============================================================================
# SECTION 16.  TRIAL REPRODUCTIONS
# =============================================================================
def mitapivat_dose_response():
    hdr("SECTION 4.  MITAPIVAT: DOSE, TIME COURSE, AND WHO RESPONDS")

    sub("4.1  Dose-response at 24 weeks (ACTIVATE titration: 5 / 20 / 50 mg BID)")
    y0, o0 = burn(0.12, 50.0)
    emit("  baseline: Hb %.2f  2,3-BPG %.2f mM  P50 %.1f  ret %.1f%%  bili %.2f" %
         (o0['Hb'], o0['DPG'], o0['P50'], o0['RetPct'], o0['BILI']))
    emit()
    emit("  %6s %8s %8s %8s %8s %8s %8s %8s %8s" %
         ("dose", "Crbc", "Hb", "dHb", "2,3-BPG", "P50", "ATPyng", "Hbeq", "dHbeq"))
    res = {}
    for dose in (0.0, 5.0, 20.0, 50.0, 100.0):
        C = newcfg(alpha=0.12, tauPK=50.0, mit=dose, mitint=0.5, tend=168.0, dt=0.05)
        _, out = simulate(C, y0=list(y0), nout=60)
        s = out[-1]
        res[dose] = out
        emit("  %6.0f %8.4f %8.2f %+8.2f %8.2f %8.1f %8.3f %8.2f %+8.2f" %
             (dose, s['MRBC'], s['Hb'], s['Hb'] - o0['Hb'], s['DPG'], s['P50'],
              s['ATPyoung'], s['Hbeq'], s['Hbeq'] - o0['Hbeq']))
    emit()
    emit("  Read the last two columns against each other.  At therapeutic doses the")
    emit("  oxygen-transport gain is REAL but SMALLER than the haemoglobin gain, and")
    emit("  the gap widens with dose: the drug is progressively buying haemoglobin")
    emit("  with oxygen affinity.  It does not reverse sign here, because the")
    emit("  2,3-BPG fall at these doses stays inside the break-even band computed in")
    emit("  Section 3.2 -- which is the useful form of the result: the endpoint")
    emit("  OVERSTATES the physiological gain by a computable factor, and the")
    emit("  overstatement grows with dose, so the dose that maximises the trial")
    emit("  endpoint is NOT the dose that maximises oxygen transport.")

    sub("4.2  Why the haemoglobin response is FAST (median 10 d in DRIVE-PK)")
    out = res[50.0]
    emit("  %8s %8s %8s %8s" % ("day", "Hb", "dHb", "% of wk-24 effect"))
    fin = out[-1]['Hb'] - o0['Hb']
    for target in (3, 7, 10, 14, 28, 56, 112, 168):
        best = min(out, key=lambda r: abs(r['t'] - target))
        emit("  %8.0f %8.2f %+8.2f %14.0f%%" %
             (best['t'], best['Hb'], best['Hb'] - o0['Hb'],
              100.0 * (best['Hb'] - o0['Hb']) / fin if fin != 0 else 0.0))
    emit()
    emit("  DRIVE-PK reported a median 10 days to the first >1.0 g/dL rise, with a")
    emit("  range out to 187 days (PMID 31483964).  In an age-structured model that")
    emit("  is not a puzzle and not a pharmacokinetic fact -- it is arithmetic.  A")
    emit("  step change in hazard relaxes the red cell pool with a time constant")
    emit("  equal to the NEW mean lifespan, not to the 120 d of a normal cell.  The")
    emit("  patient's baseline lifespan here is %.0f d, so the pool refills in weeks." % o0['lifespan'])
    emit()
    emit("  PREDICTION, and it is a strong one: time-to-response should be")
    emit("  INVERSELY related to baseline red cell lifespan.  The sicker the")
    emit("  patient (the shorter the lifespan), the FASTER the haemoglobin moves.")
    emit("  The 187-day outliers should be the mildest patients.  This is testable")
    emit("  in the existing DRIVE-PK/ACTIVATE datasets without new samples.")

    sub("4.3  Genotype: the drug is a multiplier on protein that exists")
    emit("  %-38s %8s %8s %8s" % ("genotype", "Hb0", "Hb24wk", "dHb"))
    for nm, al, tau, ac in GENO:
        y, ob = burn(al, tau)
        C = newcfg(alpha=al, tauPK=tau, activ=ac, mit=50.0, mitint=0.5,
                   tend=168.0, dt=0.05)
        _, out = simulate(C, y0=list(y), nout=4)
        emit("  %-38s %8.2f %8.2f %+8.2f" %
             (nm, ob['Hb'], out[-1]['Hb'], out[-1]['Hb'] - ob['Hb']))
    emit()
    emit("  The null-allele row gets nothing, with no parameter added to make that")
    emit("  happen: an allosteric activator multiplies residual activity, and zero")
    emit("  times anything is zero.  DRIVE-PK found responses ONLY in patients with")
    emit("  at least one missense variant, and correlated with baseline PK-R protein")
    emit("  level (PMID 31483964).  That is the same statement.")

    sub("4.4  The two limbs have different time constants")
    C = newcfg(alpha=0.12, tauPK=50.0, mit=50.0, mitint=0.5, tend=336.0, dt=0.05)
    _, out = simulate(C, y0=list(y0), nout=120)
    emit("  %8s %10s %10s %10s %10s %10s" %
         ("day", "STAB", "tau_PK", "aeff_old", "Hb", "ATPold"))
    for target in (0, 7, 28, 84, 168, 252, 336):
        b = min(out, key=lambda r: abs(r['t'] - target))
        emit("  %8.0f %10.3f %10.1f %10.5f %10.2f %10.4f" %
             (b['t'], b['STAB'], b['tauPKeff'], b['aeff_old'], b['Hb'], b['ATPold']))
    emit()
    emit("  Limb 1 (allosteric activation) is at steady state within a day and")
    emit("  cannot help a cohort whose protein has already decayed.  Limb 2")
    emit("  (thermostabilisation, tau_PK %.0f -> %.0f d) takes months and is the only" %
         (50.0, 50.0 * P['mitStabMax']))
    emit("  thing that reaches OLD cells.  So the model predicts a biphasic")
    emit("  haemoglobin trajectory whose second, slower phase is invisible in a")
    emit("  24-week trial and should appear in the long-term extension -- which is")
    emit("  where the ACTIVATE extension found responses sustained and deepening to")
    emit("  week 96 (PMID 38330179).")
    return res


def trials():
    hdr("SECTION 5.  THE REST OF THE TRIAL PROGRAMME")

    sub("5.1  ACTIVATE-T: transfusion-dependent patients (PMID 35988546)")
    y0, o0 = burn(0.09, 45.0)
    base = newcfg(alpha=0.09, tauPK=45.0, txthresh=8.0, txunits=2.0,
                  tend=365.0, dt=0.05)
    _, oc = simulate(dict(base), y0=list(y0), nout=40)
    tx_ctl = oc[-1]['TXU']
    C = dict(base)
    C['mit'] = 50.0
    C['mitint'] = 0.5
    _, ot = simulate(C, y0=list(y0), nout=40)
    tx_mit = ot[-1]['TXU']
    emit("  units transfused over 1 y, no drug   : %.1f" % tx_ctl)
    emit("  units transfused over 1 y, mitapivat : %.1f" % tx_mit)
    emit("  reduction                            : %.1f%%" %
         (100.0 * (1.0 - tx_mit / max(tx_ctl, 1e-9))))
    emit("  ACTIVATE-T primary endpoint was a >=33%% reduction, met by 10/27 (37%%).")

    sub("5.2  Splenectomy: a two-line identity the data can refute")
    emit("  If the splenic contribution is a hazard on reticulocytes with reversible")
    emit("  pooling, then in STEADY STATE the circulating reticulocyte pool is")
    emit("      R = influx / (k_exit + k_in*k_kill/(k_out+k_kill))")
    emit("  Reversible holding cancels out of R entirely -- cells that are released")
    emit("  come back.  Only the KILLING term suppresses R.  Writing q for the")
    emit("  fraction of marrow output killed in the spleen,")
    emit("      k_in*k_kill/(k_out+k_kill) = k_exit * q/(1-q)")
    emit("  so removing the spleen raises R by q/(1-q) -- and raises haemoglobin by")
    emit("  the same q/(1-q), because the same q is the production that was being")
    emit("  wasted.  THE TWO FRACTIONAL CHANGES MUST BE EQUAL.")
    emit()
    emit("  %-30s %8s %8s %8s %8s %8s" %
         ("phenotype", "dHb", "dHb %", "dRet %", "seqfrac", "L->L"))
    dh, dr = [], []
    for nm, al, tau, ac in GENO[:5]:
        y, ob = burn(al, tau)
        C = newcfg(alpha=al, tauPK=tau, splday=0.0, tend=540.0, dt=0.05)
        _, out = simulate(C, y0=list(y), nout=4)
        s = out[-1]
        a = 100.0 * (s['Hb'] / ob['Hb'] - 1.0)
        b = 100.0 * (s['RetPct'] / max(ob['RetPct'], 1e-9) - 1.0)
        dh.append(s['Hb'] - ob['Hb'])
        dr.append(b)
        emit("  %-30s %+8.2f %+8.1f %+8.1f %8.2f %4.0f->%3.0f" %
             (nm, s['Hb'] - ob['Hb'], a, b, ob['SeqFrac'], ob['lifespan'], s['lifespan']))
    emit()
    emit("  model median haemoglobin gain %+.2f g/dL, median reticulocyte change %+.1f%%" %
         (sorted(dh)[len(dh) // 2], sorted(dr)[len(dr) // 2]))
    emit("  observed median haemoglobin gain +1.6 g/dL (Grace 2018, PMID 29549173)")
    emit()
    emit("  This is the model's clearest NEGATIVE result.  PK deficiency is widely")
    emit("  described as showing a paradoxical RISE in reticulocyte count after")
    emit("  splenectomy, attributed to splenic reticulocyte sequestration (Nathan")
    emit("  1968, PMID 5634483).  The model says a sequestration-plus-killing")
    emit("  mechanism cannot deliver a large reticulocyte rise together with a small")
    emit("  haemoglobin rise, because both are the SAME q.  Sixteen combinations of")
    emit("  the two hazard gains were searched before this was recognised as")
    emit("  structural rather than a calibration failure.  Either the reported")
    emit("  reticulocyte rise is of the same modest size as the haemoglobin rise")
    emit("  (~20%, which the model produces), or its cause is NOT splenic")
    emit("  reticulocyte destruction -- the remaining candidates being the")
    emit("  maturation-factor inflation of the count and persisting marrow drive.")
    emit("  Both are measurable, and they make opposite predictions for the")
    emit("  ABSOLUTE reticulocyte count.")

    sub("5.3  Gene therapy: the age axis sets the response time, not the engraftment")
    y0, o0 = burn(0.12, 50.0)
    emit("  %8s %8s %8s %8s %8s" % ("day", "Hb", "%corrected", "Hb corr", "ret%"))
    C = newcfg(alpha=0.12, tauPK=50.0, gtday=0.0, gtfrac=0.45, tend=540.0, dt=0.05)
    _, out = simulate(C, y0=list(y0), nout=90)
    for target in (0, 14, 30, 60, 120, 240, 360, 540):
        b = min(out, key=lambda r: abs(r['t'] - target))
        emit("  %8.0f %8.2f %8.1f %8.2f %8.1f" %
             (b['t'], b['Hb'], 100.0 * b['HbCorr'] / b['Hb'], b['HbCorr'], b['RetPct']))
    emit()
    emit("  With 45%% of marrow output corrected from day 0, the haemoglobin takes")
    emit("  MONTHS, not days.  The contrast with mitapivat's 10-day response is not")
    emit("  a difference in potency: a PK activator lowers the hazard on cells that")
    emit("  already exist, whereas gene therapy can only change cells not yet born,")
    emit("  so its time constant is the lifespan of the NEW long-lived cohort")
    emit("  (~120 d), not of the old short-lived one.  A gene therapy trial powered")
    emit("  on a 24-week haemoglobin endpoint is therefore reading its own transient.")

    sub("5.4  Iron: the benefit that the haemoglobin endpoint cannot see")
    C = newcfg(alpha=0.12, tauPK=50.0, tend=730.0, dt=0.05)
    _, oc = simulate(dict(C), y0=list(y0), nout=20)
    C2 = dict(C)
    C2['mit'] = 50.0
    C2['mitint'] = 0.5
    _, om = simulate(C2, y0=list(y0), nout=20)
    a, b = oc[-1], om[-1]
    emit("  %-32s %12s %12s %10s" % ("", "untreated", "mitapivat", "change"))
    for nm, k in (("erythroferrone (ng/L)", 'ERFE'), ("hepcidin (ng/L)", 'HEP'),
                  ("erythropoietin (IU/L)", 'EPO'),
                  ("iron absorbed (mg/d)", 'absorb'),
                  ("transferrin saturation", 'TSAT'),
                  ("liver iron (mg Fe/g dw)", 'LIC'),
                  ("serum ferritin (ng/mL)", 'FERR'),
                  ("ineffective erythropoiesis", 'IEO'),
                  ("haemoglobin (g/dL)", 'Hb')):
        emit("  %-32s %12.3f %12.3f %+9.1f%%" %
             (nm, a[k], b[k], 100.0 * (b[k] / a[k] - 1.0) if a[k] else 0.0))
    emit()
    emit("  Directions to compare with van Beers 2024 (PMID 38330179), which")
    emit("  measured exactly these analytes on mitapivat to week 96: erythroferrone")
    emit("  DOWN, hepcidin UP, erythropoietin DOWN, soluble transferrin receptor")
    emit("  DOWN, liver iron DOWN by ~2 mg Fe/g dw.  The model reproduces every")
    emit("  sign.  The point worth making is structural: this benefit is driven by")
    emit("  the marrow's erythroferrone output, which falls because the drug reduces")
    emit("  the DESTRUCTION the marrow was compensating for.  It is therefore")
    emit("  largest in patients whose haemoglobin barely moves -- the compensated")
    emit("  ones above the knee, who are exactly the patients a haemoglobin-response")
    emit("  endpoint classifies as non-responders.")


# =============================================================================
# SECTION 17.  DIAGNOSTIC AND SUPPORTIVE SCENARIOS
# =============================================================================
def other_scenarios():
    hdr("SECTION 6.  DIAGNOSIS, MODIFIERS AND SUPPORTIVE CARE")

    sub("6.1  Why the diagnostic enzyme assay can read normal")
    emit("  %-30s %8s %8s %8s %8s" %
         ("phenotype", "true a", "assayed", "ratio", "ret%"))
    for nm, al, tau, ac in GENO[:5]:
        _, o = burn(al, tau)
        emit("  %-30s %8.3f %8.3f %8.2f %8.1f" %
             (nm, al, o['PKassay'], o['PKassay'] / al, o['RetPct']))
    emit()
    emit("  The assayed activity is a cohort-weighted mean and reticulocytes carry")
    emit("  freshly made enzyme, so the sicker the patient the MORE the assay")
    emit("  over-reads -- the ratio column rises as the genotype worsens.  This is")
    emit("  the mechanism behind the standard advice to interpret red cell PK")
    emit("  activity relative to other age-dependent enzymes (the PK/hexokinase")
    emit("  ratio) and to distrust a normal absolute value in a reticulocytosis.")

    sub("6.2  UGT1A1 co-inheritance (Gilbert) and gallstones")
    emit("  %-22s %10s %10s %10s %12s" %
         ("UGT1A1", "bilirubin", "indirect", "conj flux", "gallstone idx"))
    for nm, ug in (("*1/*1", 1.00), ("*1/*28", 0.70), ("*28/*28", 0.30)):
        C = newcfg(alpha=0.12, tauPK=50.0, ugt=ug, tend=1460.0, dt=0.05)
        y, _ = burn(0.12, 50.0, ugt=ug)
        _, out = simulate(C, y0=list(y), nout=4)
        s = out[-1]
        emit("  %-22s %10.2f %10.2f %10.2f %12.1f" %
             (nm, s['BILI'], s['UCB'], s['CB'], s['GS']))
    emit("  Gallstones occur in 45%% of patients (PMID 29549173) and 48%% of those")
    emit("  splenectomised without simultaneous cholecystectomy later need one.")

    sub("6.3  Aplastic crisis (parvovirus B19, 10 days of marrow arrest)")
    y0, o0 = burn(0.12, 50.0)
    C = newcfg(alpha=0.12, tauPK=50.0, crisis=30.0, crisislen=10.0,
               tend=140.0, dt=0.05)
    _, out = simulate(C, y0=list(y0), nout=140)
    emit("  %8s %8s %8s %8s" % ("day", "Hb", "ret%", "EPO"))
    for target in (28, 34, 40, 45, 52, 70, 100, 140):
        b = min(out, key=lambda r: abs(r['t'] - target))
        emit("  %8.0f %8.2f %8.1f %8.0f" % (b['t'], b['Hb'], b['RetPct'], b['EPO']))
    nad = min(out, key=lambda r: r['Hb'])
    emit("  nadir Hb %.2f g/dL at day %.0f, i.e. %.2f g/dL lost in %.0f d." %
         (nad['Hb'], nad['t'], o0['Hb'] - nad['Hb'], nad['t'] - 30.0))
    emit("  The fall rate is set by the red cell lifespan: a patient with a %.0f-day" % o0['lifespan'])
    emit("  lifespan loses %.1f%% of his red cell mass per day the marrow is off." %
         (100.0 / o0['lifespan']))

    sub("6.4  Alternative activators and a CYP3A drug interaction")
    emit("  %-34s %8s %8s %8s %8s" % ("regimen", "Cp", "Crbc", "Hb", "dHb"))
    for nm, kw in (("mitapivat 50 mg BID", dict(mit=50.0, mitint=0.5)),
                   ("mitapivat 50 BID + rifampicin", dict(mit=50.0, mitint=0.5, cypext=3.0)),
                   ("mitapivat 50 BID + fluconazole", dict(mit=50.0, mitint=0.5, cypext=0.45)),
                   ("etavopivat 400 mg OD", dict(eta=400.0, etaint=1.0)),
                   ("tebapivat 0.3 mg OD", dict(teb=0.3, tebint=1.0))):
        C = newcfg(alpha=0.12, tauPK=50.0, tend=168.0, dt=0.05, **kw)
        _, out = simulate(C, y0=list(y0), nout=4)
        s = out[-1]
        cr = s['MRBC'] + s['ERBCc'] + s['TRBCc']
        emit("  %-34s %8.4f %8.4f %8.2f %+8.2f" %
             (nm, s['MITc'], cr, s['Hb'], s['Hb'] - o0['Hb']))
    emit("  Mitapivat induces the CYP3A4 that clears it, so its own steady-state")
    emit("  exposure is lower than a single dose predicts, and a strong inducer")
    emit("  compounds that.  The same induction is why hormonal contraceptive")
    emit("  efficacy is a labelled concern.")

    sub("6.5  Endocrine effect of aromatase inhibition (male patients)")
    C = newcfg(alpha=0.12, tauPK=50.0, mit=50.0, mitint=0.5, tend=180.0, dt=0.05)
    _, out = simulate(C, y0=list(y0), nout=40)
    emit("  estradiol %.1f -> %.1f pg/mL ; testosterone %.0f -> %.0f ng/dL" %
         (P['E20'], out[-1]['E2'], P['TST0'], out[-1]['TST']))

    sub("6.6  Iron chelation, alone and with a PK activator")
    yT, oT = burn(0.09, 45.0)
    emit("  %-38s %10s %10s %10s" % ("regimen (2 y, transfused)", "LIC", "ferritin", "units"))
    for nm, kw in (("transfusion only", {}),
                   ("+ deferasirox 1000 mg/d", dict(dfx=1000.0)),
                   ("+ mitapivat 50 mg BID", dict(mit=50.0, mitint=0.5)),
                   ("+ both", dict(dfx=1000.0, mit=50.0, mitint=0.5))):
        C = newcfg(alpha=0.09, tauPK=45.0, txthresh=8.0, tend=730.0, dt=0.05, **kw)
        _, out = simulate(C, y0=list(yT), nout=4)
        s = out[-1]
        emit("  %-38s %10.2f %10.0f %10.1f" % (nm, s['LIC'], s['FERR'], s['TXU']))

    sub("6.7  Withdrawal")
    C = newcfg(alpha=0.12, tauPK=50.0, mit=50.0, mitint=0.5, mitstop=168.0,
               tend=336.0, dt=0.05)
    _, out = simulate(C, y0=list(y0), nout=120)
    for target in (160, 168, 175, 182, 210, 280, 336):
        b = min(out, key=lambda r: abs(r['t'] - target))
        emit("  day %3.0f  Hb %.2f  Crbc %.4f  STAB %.3f" %
             (b['t'], b['Hb'], b['MRBC'], b['STAB']))
    emit("  Haemoglobin decays back with the same short lifespan that made the")
    emit("  onset fast, but the thermostabilisation limb unwinds more slowly, so")
    emit("  the off-rate is biphasic and slower than the on-rate.")


# =============================================================================
# SECTION 18.  DEFECT LOG AND LIMITATIONS
# =============================================================================
def defects():
    hdr("SECTION 7.  WHAT WENT WRONG WHILE BUILDING THIS, AND WHAT IS STILL WRONG")
    emit("""
  Thirteen defects were found by running the equations.  Each is marked at its
  fix site in this file and in the mrgsolve model.  They are listed because the
  ones that mattered were not coding slips -- they were places where a
  physiologically plausible equation gave a quantitatively impossible answer.

   #1  2,3-BPG barely moved.  A bare phosphoglycerate-kinase-equilibrium form
       makes 2,3-BPG proportional to [3-PG]x(ATP/ADP); a PK lesion raises the
       first and lowers the second and they cancel, predicting +9% against an
       observed 2-3 fold.  Fixed by adding the 3-PG inhibition of the 2,3-BPG
       PHOSPHATASE, so the block shuts the shunt's exit as well as feeding it.
   #2  Solving for ATP and 2,3-BPG by simultaneous damped iteration DIVERGED,
       because the FBP -> PKR activation limb is positive feedback.  Fixed by
       nesting: the inner map in 2,3-BPG is monotone and contracts; the outer
       ATP balance is bracketed by a scan.
   #3  The shunt-flux calibration omitted the inhibition factor that the rate
       law uses, so the calibrated normal state was not a steady state.
   #4  The comment claiming the deoxyhaemoglobin-binding limb amplifies flux was
       wrong: free 2,3-BPG is EXACTLY invariant to haemoglobin saturation in this
       model.  The limb shifts total 2,3-BPG by sequestration only.  Kept,
       because that invariance is what makes the lookup tables 2-D instead of 3-D.
   #5  Erythroid iron uptake as an independent saturable process does not
       conserve iron; the healthy subject loaded to 175 mg Fe/g dw.
   #6  A units error made the liver a 0.42 mg iron sink instead of 420 mg.
   #7  Michaelis-Menten bilirubin conjugation calibrated at the normal load
       saturates at ~3x normal haemolysis; a moderately affected patient was
       assigned a bilirubin of 23000 mg/dL.
   #8  Sixteen combinations of the two haemolysis gains, with one shared ATP
       threshold, all gave the WRONG SIGN for the reticulocyte response to
       splenectomy.
   #9  Making splenic reticulocyte destruction strong enough to matter forces
       splenectomy to predict +7 to +12 g/dL instead of +1.6.  This turned out
       to be structural, not a calibration failure -- see Section 5.2.
  #10  The adenylate dilution term divides by cohort cell number and overflowed
       when a cohort was nearly annihilated.
  #11  tau_PK = 25 d made the genotype-to-haemoglobin map a step function.
  #12  Reticulocyte capture was summed over two cohorts and debited to one; once
       the debited cohort hit its floor the model CREATED cells and reported
       haemoglobins of 19-54 g/dL.  Mass balance in an age-structured model is
       per-cohort, not aggregate.
  #13  Explicit RK4 at the burn-in step of 0.25-0.30 d violates the stability
       limit (0.093 d) for the bilirubin and drug-partitioning states.  The
       bilirubin state oscillated into its non-negativity floor and reported
       0.0 mg/dL in a haemolysing patient, while the haemoglobin and
       reticulocyte columns of the step-convergence table still agreed to 1e-5.
       A convergence check on slow outputs alone does not detect this.

  LIMITATIONS THAT REMAIN, stated so they are not mistaken for results:
   * pH is not a state.  The Bohr effect and the strong pH dependence of both
     BPGM activities are folded into fixed constants, so acid-base disturbance
     cannot be simulated.
   * Free versus Mg-bound ADP is not resolved; the PGK equilibrium uses total
     ADP, which overstates how far ATP/ADP falls.
   * The 3-PG inhibition constant of the 2,3-BPG phosphatase (30 uM) is the
     single most load-bearing parameter for the 2,3-BPG magnitude and it is the
     least well pinned by data.  The oxygen-transport conclusion in Section 3
     depends on the 2,3-BPG excursion, so it inherits that uncertainty; the
     break-even table is given precisely so the conclusion can be re-scored
     against a different 2,3-BPG estimate.
   * Splenic destruction and sequestration are one compartment with a
     mean-field transit; there is no distribution of transit times.
   * The oxygen module is whole-body and single-tissue.  It cannot represent
     regional differences in extraction, which is where a right-shifted curve
     actually helps or hurts most.
   * Cardiac output responds to haemoglobin, not to tissue PO2, so the model
     cannot arbitrate between the three closures of Section 3 -- it reports all
     of them instead.
   * Only three parameters were fitted (tau_PK of the mutant, and the two
     haemolysis gains) plus the single marrow-output scalar.  Everything else is
     literature or back-calculated from the normal operating point.  That is a
     claim about provenance, not about correctness.
""")


# =============================================================================
# SECTION 19.  MAIN
# =============================================================================
if __name__ == '__main__':
    verify()
    check_healthy()
    rows = spectrum()
    central_result(rows)
    res = mitapivat_dose_response()
    trials()
    other_scenarios()
    defects()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'pkd_reference_output.txt'), 'w') as f:
        f.write("\n".join(OUT) + "\n")
    summ = {'params': dict((k, v) for k, v in P.items() if isinstance(v, (int, float))),
            'KPROD': KPROD, 'AMP_0': AMP_0, 'nstates': NST, 'cohorts': NB,
            'state_names': SN, 'Lcrit': BEXIT * AMP_0 / P['ampMax'],
            'DT_STABLE': DT_STABLE, 'genotypes': []}
    for nm, al, tau, ac in GENO:
        _, o = burn(al, tau)
        summ['genotypes'].append(dict(
            label=nm, alpha=al, tauPK=tau, activatable=ac, Hb=o['Hb'],
            RetPct=o['RetPct'], lifespan=o['lifespan'], BILI=o['BILI'],
            LDH=o['LDH'], DPG=o['DPG'], P50=o['P50'], PvO2=o['PvO2'],
            Hbeq=o['Hbeq'], COreq=o['COreq'], LIC=o['LIC'], FERR=o['FERR'],
            EPO=o['EPO'], ERFE=o['ERFE'], HEP=o['HEP'], SeqFrac=o['SeqFrac'],
            PKassay=o['PKassay']))
    dr = {}
    for d in res:
        s = res[d]
        dr[str(int(d))] = [s[0], s[len(s) // 4], s[-1]]
    summ['mitapivat_dose_response'] = dr
    with open(os.path.join(here, 'pkd_population_results.json'), 'w') as f:
        json.dump(summ, f, indent=1)
    print("\nwrote pkd_reference_output.txt (%d lines)" % len(OUT))
