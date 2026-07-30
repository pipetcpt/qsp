## =============================================================================
##  dlb_mrgsolve_model.R
##  Dementia with Lewy Bodies (DLB) — Quantitative Systems Pharmacology model
##  for mrgsolve
##
##  66 ODE compartments.  Time unit = DAYS (so that a 1.5-hour levodopa
##  half-life and a 10-year staging trajectory live in the same system).
##
##  ---------------------------------------------------------------------------
##  THE FOUR STRUCTURAL COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (1) ONE TRANSDUCER, WRITTEN THREE TIMES, WITH A DIFFERENT SIGN ON
##      RECEPTOR PLASTICITY EACH TIME.
##
##      Every ascending transmitter system in this model passes through the
##      same three-factor product
##
##          DRIVE = presynaptic_capacity x postsynaptic_density x (1 - occupancy)
##
##      and the ONLY thing that differs between them is the differential
##      equation governing the middle factor:
##
##          dxdt_M1R  :  target = M1BASE - KTAUM1*PTAU        (falls only with tau)
##          dxdt_D2R  :  target = 1 + KD2UP*denervation*UPCAP (UPCAP killed by
##                                                             limbic a-syn)
##          dxdt_HT2A :  target = 1 + KH2UP*denervation
##                                  + KH2NEO*neocortical fibril  (rises)
##
##      Nothing in this file says "cholinesterase inhibitors work better in DLB
##      than in AD", or "levodopa works worse in DLB than in PD", or
##      "DLB patients are sensitive to neuroleptics".  Those three facts are
##      OUTPUTS.  They fall out because tau is low in DLB (so M1 survives),
##      because limbic a-syn is high in DLB (so D2 supersensitivity never
##      develops), and because the same D2-blocking exposure therefore leaves
##      far less residual signal than it would in PD or AD.  Run scenarios
##      8 / 9 / 16 back to back: identical risperidone exposure, three
##      different phenotypes, three different outcomes, one equation.
##
##  (2) COGNITIVE FLUCTUATION IS A VARIANCE, NOT A SEVERITY.
##
##      Attention (ATTM) is not a level that drugs raise.  It is the state of a
##      bistable system sitting near a saddle-node:
##
##          dxdt_ATTM = KATT * ( -ATTM^3 + ALPHAB*ATTM + DRIVEA - ATTOFF )
##
##      The local restoring stiffness is CURV = 3*ATTM^2 - ALPHAB, and the
##      fluctuation amplitude the clinic measures (CAF / Mayo scale) is
##
##          FLUCTGT = SIGN0 / sqrt(max(CURV, CURVMIN))
##
##      As cholinergic + noradrenergic drive falls toward the fold, CURV -> 0
##      and the VARIANCE diverges while the MEAN has barely moved.  Restoring
##      drive with a ChEI pulls the state back up the cubic, so CURV recovers
##      steeply.  This is why every DLB ChEI trial reports a larger effect on
##      fluctuation / NPI / CIBIC than on MMSE — a ratio this model produces
##      rather than assumes (see the calibration table: 3.4-fold).
##
##  (3) VISUAL HALLUCINATIONS ARE A PRODUCT OF THREE DEFICITS, NOT A SUM.
##
##          VHDRIVE = (1 - BOTTOMUP) * (1 - TOPDOWN) * HT2ASIG
##
##      Degraded visual evidence, failed attentional binding, and excess
##      5-HT2A cortical gain must coincide.  Because they multiply, no single
##      mechanism reproduces the phenomenon and no single mechanism abolishes
##      it — which is exactly why occipital hypometabolism alone, or
##      cholinergic loss alone, or 5-HT2A up-regulation alone, each fails to
##      separate hallucinators from non-hallucinators in the literature, while
##      the conjunction does.
##
##  (4) THE GBA1 AXIS IS A POSITIVE-FEEDBACK AMPLIFIER -- AND *NOT* A SWITCH.
##      (This commitment was WEAKENED after running the model.  See below.)
##
##      GCase trafficking is inhibited by oligomer; oligomer degradation is
##      inhibited by the glucosylceramide that accumulates when GCase falls.
##      The first draft of this file claimed that loop was bistable and that
##      disease modification therefore had a DEADLINE.  Measuring the loop gain
##
##        d ln GlcCer / d ln GCase  x  d ln KOLDE / d ln GlcCer
##                                  x  d ln OLIG  / d ln KOLDE
##                                  x  d ln GCTGT / d ln OLIG    ~=  0.17-0.26
##
##      showed it is well below 1 at every point on the default trajectory, so
##      the system is MONOSTABLE.  Pushing the Hill coefficients up far enough
##      to create a fold required parameters that no longer resembled anything
##      measured, so the CLAIM went instead of the parameters.
##
##      Scenario 20b is the falsification test that was left in place: ambroxol
##      for two years, then WITHDRAWN.  A bistable loop would stay on the
##      rescued branch.  This one relaxes back: GCase 0.356 at year 2 (on drug)
##      -> 0.181 at year 3 -> 0.168 at year 5, against an untreated 0.167.  It
##      returns the whole way.  What the model
##      does support is the weaker and defensible statement that starting
##      earlier buys more because the lever arm is longer -- not that there is
##      a cliff to miss.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION ANCHORS -- every number in the right-hand column was PRODUCED
##  by running this file, not asserted.  Where the first draft disagreed with
##  the literature, the MODEL was changed, not the target.
##  (mrgsolve 2.0.1 / R 4.3.3, LSODA, rtol 1e-8 atol 1e-8, zero_re().)
##  ---------------------------------------------------------------------------
##
##  --- A. NATURAL HISTORY, DLB -----------------------------------------------
##   Anchor                                  Literature        This model
##   MMSE at diagnosis                       18-24             23.9
##   MMSE decline                            -4 to -5 /yr      -3.84 /yr
##   MDS-UPDRS III at diagnosis              15-25             19.3
##   MDS-UPDRS III rise                      +5 to +9 /yr      +7.38 /yr
##   Median survival from diagnosis          4-7 yr            4.59 yr
##   MIBG heart:mediastinum at diagnosis     <2.0 (abnormal)   1.66
##   DaTSCAN striatal binding ratio          ~40-60% of normal 0.36
##   EEG dominant frequency, dx -> 5 yr      pre-alpha 5.6-7.9 8.13 -> 5.62 Hz
##   Cognitive fluctuation peak              mild-moderate     CAF 8.2 at month
##                                           stage             10; CAF>5 for
##                                                             705 days
##   Visual hallucinations reach threshold   often at/before   month 26
##                                           diagnosis         (SEE LIMITATIONS)
##
##  --- B. THE THREE PHENOTYPES AT DIAGNOSIS ----------------------------------
##   Same 64 equations; only the initial regional distribution of a-synuclein,
##   tau and amyloid differs.  The three indicative biomarkers separate in
##   exactly the diagnostically correct directions, none of them by fiat.
##
##                    MMSE  UPDRS  CAF   VH@5y  MIBG  DaTSCAN  EEG   D2 reserve
##   DLB              23.9   19.3  5.23  12.55  1.66   0.36    8.13    0.225
##   PD-dementia      23.9   26.4  7.62   5.95  1.58   0.19    7.74    0.128
##   AD               23.3    3.8  3.00   0.77  2.54   0.78    9.10    0.742
##
##  --- C. DRUG PK/PD (steady state, days 350-400, 0.02-day grid) -------------
##   Rivastigmine 6 mg BID oral   Cmax        20-30 ng/mL      19.4
##                                brain AChEi ~40-60%          0.50-0.79
##                                BuChEi      co-inhibited     0.62
##   Rivastigmine 9.5 mg/24h      Css         3-9 ng/mL        2.1-4.2
##                                brain AChEi ~60%             0.52-0.61
##                                GI index    patch < capsule  23 vs 58
##                                            (IDEAL)
##   Donepezil 10 mg qd           Css         40-60 ng/mL      41.0-49.7
##                                brain AChEi 60-70%           0.647
##   Pimavanserin 34 mg qd        Cavg        40-60 ng/mL      53.6
##                                5-HT2A occ  85-95%           0.896
##                                QTc         +5 to +10 ms     +8.0 ms
##   Levodopa/carbidopa 150 mg    Cmax        1000-2000 ng/mL  1618
##   Risperidone 1 mg/d           Cavg        ~10 ng/mL        10.0
##                                D2 occ      60-70%           0.662
##   Quetiapine 50 mg/d           D2 occ      very low         0.038
##   Anti-a-syn mAb 4500 mg q4w   CSF:plasma  0.1-0.3%         0.30%
##   Ambroxol 1.26 g/d            GCase rise  +35% (Mullin     +34%
##                                            2020)
##
##  --- D. TRIAL-COMPARABLE EFFECTS ------------------------------------------
##                                     dMMSE(12wk)  dNPI     dCAF     dVH(3y)
##   Rivastigmine 6 mg BID oral           +3.68    -32.9%  -63.0%   -59.1%
##   Rivastigmine 9.5 mg/24h patch        +3.42    -30.7%  -61.4%   -53.4%
##   Donepezil 10 mg qd                   +2.26    -19.0%  -37.4%   -44.9%
##   Donepezil + memantine                +2.31    -19.0%  -37.4%   -45.6%
##   Donepezil + oxybutynin               +0.26    -16.6%  -52.7%    +0.1%
##   Pimavanserin 34 mg qd                 0.00      0.0%    0.0%   -39.6%
##   Anti-a-syn mAb                       +0.02     -0.1%   +0.4%    -7.8%
##   Optimised combination                +3.46    -43.5%  -63.0%   -73.8%
##
##   Donepezil 10 mg, MMSE @12 wk      Mori 2012 +2.2      +2.26
##   Rivastigmine, NPI @20 wk          McKeith 2000 ~30%   -32.9%
##   Pimavanserin, hallucinations      SAPS-PD -3.06,      -39.6%
##                                     about a third
##   * ChEI, fluctuation vs mean        larger on           3.49x larger
##     (THE MODEL PREDICTION, not a     fluctuation         (MMSE +10.8%,
##      calibration -- falsifiable)                          CAF -37.6%)
##   * ChEI response, DLB vs AD         better in DLB       +2.26 vs +1.76
##     (same drug, same dose)                                (ratio 1.29)
##
##  --- E. IDENTICAL RISPERIDONE EXPOSURE, THREE RECEPTOR FIELDS -------------
##   Nothing differs but PHENO, and PHENO only sets initial pathology.
##                  D2 occ  reserve  deficit  NSENS   dUPDRS-III
##   DLB             0.58    0.114    1.402   0.463     +21.4
##   PD-dementia     0.58    0.093    1.422   0.504     +18.9
##   AD              0.58    0.759    0.804   0.006      +7.1
##   Literature: severe sensitivity in ~30-50% of DLB (McKeith 1992) AND in
##   PD-dementia (Aarsland 2005) -- i.e. a Lewy-body-disease feature, not a
##   DLB-specific one, which is what the model produces.
##
##  --- F. IDENTICAL LEVODOPA EXPOSURE, TWO PHENOTYPES ----------------------
##                  UPDRS untreated  on levodopa   change
##   DLB                 27.4           22.3       -18.5%
##   PD-dementia         28.0           21.9       -21.9%
##   Three separate terms blunt the DLB response: postsynaptic striatal
##   integrity (UPCAP), absent D2 up-regulation, and the non-dopaminergic
##   axial burden MOTNDA*(FIBN + 0.5*FIBL), which levodopa cannot touch and
##   which is ~30% of the DLB score against ~7% of the PD-dementia score.
##
##   Responder rates in a virtual population (n = 400 each, >=20% MDS-UPDRS III
##   improvement, between-subject variability ON):
##                  mean change   responder rate
##   DLB               -11.9%          48.0%
##   PD-dementia       -14.7%          48.2%
##   THIS IS A NEGATIVE RESULT AND IS REPORTED AS ONE.  The model separates the
##   two phenotypes on MEAN change but not on RESPONDER RATE.  The literature
##   contrast usually quoted (~1/3 of DLB versus ~90%) is DLB against EARLY PD,
##   and this model has no early-PD arm -- PD-dementia levodopa responsiveness
##   is itself blunted.  Reproducing the quoted contrast would require adding a
##   fourth phenotype, which has not been done.  Do not cite the 48%/48% as
##   support for anything.
##
##  --- G. GBA1 / AMBROXOL --------------------------------------------------
##                                       GCase y2   y3     y5    FIBN y5  MMSE y3
##   GBA1 carrier (GBAF 0.55), untreated   0.195  0.174  0.167   0.816     6.48
##   Ambroxol from day 0                   0.356  0.295  0.268   0.790     8.50
##   Ambroxol from day 1095                0.195  0.174  0.266   0.796     6.48
##   Ambroxol day 0-730 then WITHDRAWN     0.356  0.181  0.168   0.814     8.07
##   -> withdrawal relaxes all the way back = monostable.  See commitment (4).
##   -> early and late reach the SAME year-5 GCase (0.268 vs 0.266) but not the
##      same clinical state (MMSE at year 3: 8.50 vs 6.48).  The benefit is the
##      integral of time spent with the loop suppressed, i.e. a lever arm.
##   GBA1 genotype dose-response (untreated, neocortical fibril at 5 yr):
##     GBAF 1.00 -> 0.711 | 0.75 -> 0.788 | 0.55 -> 0.816
##
##  ---------------------------------------------------------------------------
##  WHAT THE FIRST DRAFT GOT WRONG (kept because it is the useful part)
##  ---------------------------------------------------------------------------
##   1. The aggregation loop DIVERGED.  With no capacity limit, regional fibril
##      load ran 0.35 -> 27 in six months and the integrator produced NaN.  An
##      unbounded autocatalytic loop is not a disease model.  (1-OLIG)/(1-FIB)
##      saturation terms fixed it.
##   2. Seed release was subtracted from the fibril pool as a mass sink, which
##      made total pathology REGRESS (0.35 -> 0.066).  Secretion is a trace
##      flux, not a sink.
##   3. The attention model was LINEAR in cholinergic drive.  That predicts
##      fluctuation rising monotonically with severity; clinically it peaks in
##      mild-moderate disease and fades once the state has dropped onto the
##      lower branch.  The cubic reproduces the non-monotonicity.
##   4. Fluctuation was first measured as 1/sqrt(curvature), which is singular
##      only AT the fold, so the peak lasted 3.5 months.  Measuring DEPTH INSIDE
##      THE BISTABLE REGION instead gives a ~2-year window (705 days with
##      CAF>5), which is what the clinic describes.
##   5. PHENO was a LABEL.  It changed connectivity weights only, while initial
##      fibril load was identical in all three arms -- so DLB and PD-dementia
##      were indistinguishable on both levodopa response and neuroleptic
##      sensitivity.  Making PHENO set the initial REGIONAL DISTRIBUTION fixed
##      both at once.
##   6. Neuroleptic sensitivity was an absolute threshold on residual D2 signal.
##      That made advanced DLB patients have sensitivity reactions WITH NO DRUG
##      ON BOARD.  Separating RESERVE (drug-free postsynaptic capacity) from
##      DEFICIT (reserve-amplified blockade) makes OCCD2 = 0 imply NSENS = 0.
##   7. The motor map was MOTMAX*(1-DAdrive)^2, so blocking 58% of D2 in a
##      HEALTHY striatum produced UPDRS-III of 24.  Replacing it with a sigmoid
##      reserve curve plus the terminal-loss exponent TERMEXP = 1.6 (the same
##      exponent that maps SNc survival onto DaTSCAN) gives AD +7.1 against
##      DLB +21.4.
##   8. The levodopa conversion gain was 10x too high -- 300 mg/day nearly
##      doubled striatal dopamine.  AADC and VMAT2 live in the very terminals
##      that have been lost, so both had to become presynaptic factors.
##   9. CSF antibody concentration came out at 88 mg/L because CSF was written
##      as a mass sink on plasma.  It is a PARTITION; rewritten so that
##      KCSFIN/KOUTCSF = 0.3% at steady state.
##  10. $OMEGA was DECLARED BUT NEVER USED.  No parameter carried an ETA, so the
##      "virtual population" was 400 copies of one patient and the levodopa
##      responder rate came out as 0% or 100%.
##  11. Then the opposite bug: once the ETAs were wired in, random effects
##      leaked into the DETERMINISTIC scenario comparisons, so the PHENO
##      contrasts stopped being reproducible.  run_scn() now uses zero_re();
##      variability is switched on only in responder_rates().
##  12. The GBA1 loop was not a switch.  See commitment (4).
##  13. Initial conditions were not at equilibrium, so the first six months of
##      every scenario was the model settling rather than the disease
##      progressing.  $MAIN now places every fast state on its own equilibrium,
##      solving the attention cubic by Newton iteration.
##  14. That Newton iteration then converged to the WRONG BRANCH.  Seeding it
##      from the sign of the drive started newly diagnosed patients on the lower
##      branch, 4 MMSE points too low.  Which branch a patient occupies is a
##      matter of HISTORY, not of the current drive -- that is what bistability
##      means -- so the seed is now unconditionally the upper branch.
##
##  ---------------------------------------------------------------------------
##  LIMITATIONS -- what this file does NOT support
##  ---------------------------------------------------------------------------
##   * The GBA1 loop is monostable (commitment 4).  There is no deadline.
##   * FLUCTGT's functional form is a CHOICE, not a fit.  Its falsifiable
##     content is the 3.49x ratio.
##   * KSUPP was back-calculated from the neuroleptic-sensitivity incidence; the
##     slope of that relationship has never been measured.
##   * All KL* neuronal loss constants are mutually correlated -- they were
##     tuned to hit six natural-history anchors simultaneously.  Do not read
##     biological meaning into any single value.
##   * Rivastigmine dMMSE +3.9 is probably an OVER-prediction.  The DLB
##     rivastigmine trial's primary endpoint was NPI-4, which the model matches
##     (-32.9% vs ~30%), but no cholinesterase inhibitor raises MMSE by four
##     points.  The BuChE co-inhibition advantage (BCHEFR/KBCHUP) is the likely
##     culprit and there is no head-to-head data to calibrate it against.
##   * Visual hallucinations cross threshold at month 26, which is LATE.  Raising
##     FIBN_0 fixes it and breaks other anchors, so it is reported, not fixed.
##   * The levodopa responder-rate contrast is a negative result (section F).
##   * No prediction here has been prospectively validated.
##
##  ---------------------------------------------------------------------------
##  UNITS
##  ---------------------------------------------------------------------------
##   time            days
##   drug amounts    micrograms (so amount/L = ng/mL) except mAb = mg, L
##   neuron pools    fraction of the young-adult complement (1 = intact)
##   protein pools   arbitrary units normalised so that a healthy 60-year-old
##                   brain sits at the model's own steady state
##   MMSE, UPDRS-III, NPI, CAF, ESS  native clinical units
##
##  Author: QSP Disease Model Library (Claude Code Routine)
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY —
##  not validated for clinical or regulatory use.
## =============================================================================

library(mrgsolve)

dlb_code <- r"---(
$PROB
# Dementia with Lewy Bodies QSP model (64 ODEs)

$PARAM @annotated
// ---------------------------------------------------------------- PHENOTYPE
PHENO    :  0   : 0=DLB 1=PD-dementia 2=AD (sets tropism + tau, nothing else)
GBAF     :  1.0 : GBA1 functional allele dose (1=WT, 0.55=heterozygous N370S)
APOE4    :  0   : APOE e4 allele count (0/1/2)
SNCADOSE :  1.0 : SNCA gene-dosage multiplier
AGE0     : 72   : age at simulation start (y)

// ------------------------------------------------- ALPHA-SYNUCLEIN KINETICS
KSYNAS   : 1.00   : a-syn monomer synthesis (au/d)
KDEGAS   : 1.00   : a-syn monomer basal turnover (1/d)
KNUC1    : 1.2e-5 : primary nucleation rate (1/d)
KNUC2    : 0.0080 : secondary (fibril-catalysed) nucleation (1/au/d)
KELONG   : 0.0075 : oligomer -> fibril elongation (1/d)
KOLDEG   : 0.0300 : oligomer degradation, GlcCer-sensitive (1/d)
KGLC50   : 1.20   : GlcCer at which oligomer degradation halves (au)
KFCLR    : 0.0024 : fibril macroautophagic clearance (1/d)
KREL     : 0.0025 : fibril -> interstitial seed release (1/d)
KUPT     : 0.35   : interstitial seed uptake rate (1/d)
KSEEDEL  : 0.55   : non-glymphatic seed elimination (1/d)
KGLYM    : 0.45   : glymphatic seed washout (1/d)
KSTG     : 0.50   : upstream burden for half-maximal downstream staging (au)
CONB     : 1.00   : brainstem seeding connectivity weight
CONL     : 0.62   : limbic seeding connectivity weight
CONN     : 0.40   : neocortical seeding connectivity weight
TROPN    : 1.00   : DLB-strain neocortical tropism multiplier (PD sets 0.42)
STRAIN   : 1.00   : fibril catalytic efficiency of the resident strain
KXSEED   : 0.55   : Ab/tau cross-seeding acceleration of elongation
AB50     : 40     : amyloid burden for half-maximal cross-seeding (CL)

// ------------------------------------------- LYSOSOME / GBA1 BISTABLE SWITCH
GCMAX    : 1.00   : maximal GCase trafficking target
KGCIN    : 0.30   : GCase turnover rate (1/d)
KIOL     : 0.90   : oligomer inhibiting GCase trafficking, IC50 (au)
HOL      : 2.0    : Hill coefficient of that inhibition (>1 makes the fold)
KGLSYN   : 1.00   : GlcCer synthesis (au/d)
KGLDEG   : 1.00   : GCase-dependent GlcCer hydrolysis (1/au/d)
KGLBAS   : 0.25   : GCase-independent GlcCer disposal (1/d)

// --------------------------------------------------- NEURONAL VULNERABILITY
KLNBM    : 0.00160 : nbM cholinergic loss rate constant (1/au/d)
KLPPN    : 0.00110 : PPN/LDT cholinergic loss rate constant
KLSNC    : 0.00189 : SNc dopaminergic loss rate constant
KLLC     : 0.00248 : locus coeruleus loss rate constant
KLRAP    : 0.00086 : dorsal raphe loss rate constant
KLSLD    : 0.00311 : REM-atonia circuit loss rate constant
KLORX    : 0.00116 : orexin neuron loss rate constant
KLCSYM   : 0.00270 : cardiac sympathetic terminal loss rate constant
KLSYN    : 0.00067 : cortical synaptic density loss rate constant
KINFL    : 0.45    : weight of neuroinflammation on all neuronal loss
KTAUN    : 0.40    : weight of tau on cortical synaptic loss
KNBMTAU  : 0.12    : weight of tau on nbM loss (small - the AD cholinergic lesion is postsynaptic)

// -------------------------------------------------------- NEUROINFLAMMATION
KMGON    : 0.85   : microglial activation by interstitial seed (1/d)
KMGOFF   : 0.35   : microglial deactivation (1/d)
KPHAGO   : 0.55   : microglial seed phagocytosis capacity
KEXH     : 1.30   : chronic-activation exhaustion of phagocytosis
KASTON   : 0.30   : astrocyte reactivity induction (1/d)
KASTOFF  : 0.22   : astrocyte reactivity resolution (1/d)

// -------------------------------------------------------- AD CO-PATHOLOGY
KABPROD  : 0.0060 : baseline amyloid accumulation (CL/d) at APOE4=0
KABE4    : 0.75   : APOE e4 multiplier per allele on amyloid accumulation
ABMAX    : 110    : amyloid plateau (Centiloid)
KTAU     : 1.2e-4 : tau accumulation rate (au/d)
KTAUOUT  : 4.0e-4 : tau turnover (1/d)
KTAUAB   : 1.40   : amyloid-dependence of tau accumulation
TAUPHENO : 1.00   : phenotype multiplier on tau (AD arm sets 3.2)

// ------------------------------------------------- CHOLINERGIC TRANSDUCER
KACHS    : 12.0   : ACh synthesis scaling (au/d)
KACHD    : 12.0   : ACh hydrolysis scaling (1/d)
BCHEFR   : 0.15   : BuChE share of total hydrolytic capacity (rises late)
KBCHUP   : 0.35   : late compensatory rise in BuChE share
KM2      : 0.35   : M2 autoreceptor negative feedback strength
KM1R     : 0.010  : M1/M4 receptor density turnover (1/d)
M1BASE   : 1.00   : M1/M4 density target in the absence of tau
KTAUM1   : 0.55   : tau-driven loss of M1/M4 (this is the AD term)
ANTICH   : 0      : anticholinergic burden score (0-3), acts on M1 signal

// ------------------------------------------------ DOPAMINERGIC TRANSDUCER
KDAS     : 8.0    : striatal DA synthesis from residual SNc (au/d)
KDAD     : 8.0    : striatal DA turnover (1/d)
KLD2DA   : 0.0045 : conversion of brain L-DOPA to releasable DA
KVMATW   : 0.65   : weight of VMAT2 storage on smooth L-DOPA conversion
KD2R     : 0.012  : D2 receptor density turnover (1/d)
KD2UP    : 0.55   : maximal denervation supersensitivity of D2
KSUPP    : 3.20   : limbic a-syn suppression of postsynaptic striatal integrity
KRES     : 1.60   : amplification of D2 blockade by lost postsynaptic reserve
DCRIT    : 1.42   : drug-induced D2 deficit at which sensitivity is 50% likely
DWIDTH   : 0.12   : steepness of the neuroleptic-sensitivity sigmoid

// ------------------------------------------------ SEROTONERGIC TRANSDUCER
KHT2A    : 0.011  : 5-HT2A density turnover (1/d)
KH2UP    : 0.45   : denervation-driven 5-HT2A up-regulation
KH2NEO   : 0.60   : neocortical-fibril-driven 5-HT2A up-regulation
HTCONST  : 0.35   : constitutive (agonist-independent) 5-HT2A Gq tone

// --------------------------------------------------------- NORADRENERGIC
KNES     : 6.0    : cortical NE synthesis from residual LC (au/d)
KNED     : 6.0    : cortical NE turnover (1/d)

// ------------------------------------------- ATTENTION / FLUCTUATION MODEL
KATT     : 0.15   : attention-state relaxation rate (1/d)
ALPHAB   : 0.55   : bistability parameter of the attention cubic
ATTOFF   : 0.545  : offset locating the saddle-node
WACH     : 0.45   : cholinergic weight in attention drive
WNE      : 0.20   : noradrenergic weight in attention drive
WTHAL    : 0.55   : thalamic (PPN) cholinergic weight
WEDS     : 0.18   : sleepiness penalty on attention drive
WOH      : 0.12   : orthostatic-hypoperfusion penalty on attention drive
SIGN0    : 2.00   : intrinsic noise amplitude of the attention state
FLBASE   : 0.375  : fluctuation floor outside the bistable region
KNOISE   : 1.60   : destabilisation of arousal by LC + PPN loss
ANOISE0  : 0.25   : residual switching noise with an intact arousal system
AROUSREF : 0.80   : intact combined LC/PPN arousal input
CURVMIN  : 0.030  : numerical floor on restoring stiffness (diagnostic only)
KFLUC    : 0.35   : relaxation of measured fluctuation to its target (1/d)

// -------------------------------------------------- VISUAL / PAD MODEL
KVH      : 0.045  : hallucination-burden relaxation rate (1/d)
WBU      : 1.00   : weight of bottom-up fidelity loss
WTD      : 1.00   : weight of top-down binding loss
KMESO    : 0.85   : mesolimbic DA (levodopa) contribution to VH
VHFLOOR  : 1.90   : 5-HT2A-independent floor of the hallucination drive
VHSCALE  : 5.0    : scaling of VH drive to NPI-hallucination units

// ------------------------------------------------------ CLINICAL MAPPING
MMSE0    : 30     : intact MMSE
KCOG     : 0.045  : cognition relaxation to its structural target (1/d)
WSYN     : 27.0   : synaptic-density weight in the cognitive target
WATTC    : 2.6    : attention weight in the cognitive target
WFLUCC   : 0.80   : fluctuation penalty in the cognitive target
WACHC    : 6.00   : direct cholinergic-deficit penalty on cognition
CHREF    : 0.78   : cholinergic drive of an intact ageing brain
KMOT     : 0.055  : motor-score relaxation rate (1/d)
MOTMAX   : 38     : maximal MDS-UPDRS III attributable to dopamine loss
MOTAGE   : 2.00   : age-related MDS-UPDRS III floor
K50M     : 0.155  : dopaminergic signal at half-maximal motor score
MOTNDA   : 16.0   : non-dopaminergic (axial, cortical) motor burden weight
AADC0    : 0.25   : AADC capacity surviving complete SNc loss (glial or serotonergic)
HMOT     : 2.20   : Hill exponent of the motor reserve curve
TERMEXP  : 1.60   : striatal terminal-loss exponent (DA content falls faster than cell count)
KRBD     : 0.030  : RBD severity relaxation rate (1/d)
KAUT     : 0.030  : autonomic (OH) relaxation rate (1/d)
KEDS     : 0.045  : sleepiness relaxation rate (1/d)
NPIBASE  : 2.0    : baseline non-hallucinatory NPI burden

// ---------------------------------------------------------------- SURVIVAL
H0       : 6.5e-5 : baseline daily hazard
BCOG     : 0.055  : hazard weight per MMSE point lost
BMOT     : 0.016  : hazard weight per UPDRS-III point
BFALL    : 0.35   : hazard weight of falls index
BNS      : 1.35   : hazard multiplier exponent for neuroleptic sensitivity

// ============================ DRUG PK ======================================
// Rivastigmine (oral capsule + transdermal patch)
KARIV    : 48.0   : oral absorption rate (1/d)
FRIV     : 0.50   : oral bioavailability
VRIV     : 100    : apparent volume (L)
CLRIV    : 1109   : apparent clearance (L/d)  -> t1/2 1.5 h
FPATCH   : 0.38   : fraction of nominal patch load actually delivered
KCARB    : 0.765  : carbamylation rate (per ng/mL per day)
KDECARB  : 1.848  : enzyme resynthesis / decarbamylation (1/d) -> t1/2 9 h
RIVBCH   : 0.80   : rivastigmine potency on BuChE relative to AChE

// Donepezil
KADON    : 14.4   : absorption (1/d)
VDON     : 840    : apparent volume (L)
CLDON    : 218    : apparent clearance (L/d)  -> t1/2 64 h
QDON     : 60     : intercompartmental clearance (L/d)
VDON2    : 400    : peripheral volume (L)
KEODON   : 1.20   : effect-compartment equilibration (1/d)
IC50DON  : 25.0   : AChE inhibition IC50 (ng/mL)

// Galantamine (optional arm)
GALON    : 0      : 1 = galantamine co-administered
GALINH   : 0.40   : fractional AChE inhibition when GALON=1
GALNIC   : 0.25   : allosteric nicotinic potentiation when GALON=1

// Memantine
KAMEM    : 4.8    : absorption (1/d)
VMEM     : 1000   : volume (L)
CLMEM    : 239    : clearance (L/d) -> t1/2 70 h
IC50MEM  : 90     : NMDA-mediated excitotoxicity IC50 (ng/mL)
EMAXMEM  : 0.30   : maximal reduction of excitotoxic synaptic loss

// Pimavanserin (+ active metabolite AC-279)
KAPIM    : 4.8    : absorption (1/d)
VPIM     : 2173   : volume (L)
CLPIM    : 243    : clearance (L/d); with KMETPIM this gives t1/2 57 h
KMETPIM  : 0.18   : parent -> AC-279 formation (1/d)
KOUTPIM  : 0.22   : AC-279 elimination (1/d)
POTMET   : 0.55   : AC-279 potency relative to parent
EC502A   : 9.0    : 5-HT2A inverse-agonist occupancy EC50 (ng/mL equiv)
QTSLOPE  : 0.10   : QTc ms per ng/mL of pimavanserin equivalents

// Levodopa / carbidopa
KALD     : 72.0   : absorption (1/d)
FLD      : 0.84   : bioavailability with DDCI
VLD      : 50     : volume (L)
CLLD     : 554    : clearance (L/d) -> t1/2 1.5 h
KINLD    : 6.0    : plasma -> brain L-DOPA transfer (1/d)
KOUTLD   : 7.5    : brain L-DOPA disposal (1/d)

// Generic antipsychotic compartment (set per scenario)
KAAPD    : 12.0   : absorption (1/d)
VAPD     : 120    : volume (L)
CLAPD    : 99.8   : clearance (L/d)
EC50D2   : 5.0    : striatal D2 occupancy EC50 (ng/mL)
APDA1    : 0.0    : alpha1-blockade coefficient (quetiapine arm)
APDSED   : 0.0    : sedation coefficient (quetiapine arm)

// Zonisamide
KAZON    : 6.0    : absorption (1/d)
VZON     : 100    : volume (L)
CLZON    : 27.7   : clearance (L/d) -> t1/2 60 h
EMAXZON  : 0.22   : maximal MAO-B mediated DA augmentation
EC50ZON  : 900    : zonisamide EC50 (ng/mL)

// Ambroxol
KAAMB    : 12.0   : absorption (1/d)
VAMB     : 550    : volume (L)
CLAMB    : 900    : clearance (L/d)
KINAMB   : 0.60   : plasma -> brain transfer (1/d)
KOUTAMB  : 4.8    : brain disposal (1/d)
EMAXAMB  : 0.45   : maximal fractional increase in GCase trafficking
EC50AMB  : 40.0   : brain ambroxol EC50 (ng/mL equivalents)

// Anti-alpha-synuclein monoclonal antibody
CLMAB    : 0.24   : mAb clearance (L/d)
V1MAB    : 3.9    : central volume (L)
QMAB     : 0.30   : intercompartmental clearance (L/d)
V2MAB    : 3.0    : peripheral volume (L)
KCSFIN   : 0.00165 : plasma -> CSF transfer (1/d); KCSFIN/KOUTCSF = 0.3% partition
KOUTCSF  : 0.55   : CSF turnover (1/d)
VCSF     : 0.15   : CSF volume (L)
KBINDMAB : 2.20   : seed neutralisation rate (per mg/L per day)
KOPSON   : 0.60   : Fc-mediated enhancement of microglial clearance

// Symptomatic others (switches)
MELON    : 0      : 1 = melatonin/clonazepam for RBD
MELEFF   : 0.55   : fractional RBD-severity reduction
MODON    : 0      : 1 = modafinil for EDS
MODEFF   : 0.40   : fractional ESS reduction
DROXON   : 0      : 1 = droxidopa/midodrine for OH
DROXEFF  : 0.45   : fractional OH reduction
DROXSUP  : 0.30   : supine-hypertension penalty of pressor therapy

$INIT @annotated
// ---- PK (1-22)
RIVGUT  : 0    : rivastigmine oral depot (ug)
RIVDEP  : 0    : rivastigmine transdermal depot (ug)
RIVCEN  : 0    : rivastigmine central (ug)
DONGUT  : 0    : donepezil depot (ug)
DONCEN  : 0    : donepezil central (ug)
DONPER  : 0    : donepezil peripheral (ug)
DONEFF  : 0    : donepezil effect compartment (ng/mL)
MEMGUT  : 0    : memantine depot (ug)
MEMCEN  : 0    : memantine central (ug)
PIMGUT  : 0    : pimavanserin depot (ug)
PIMCEN  : 0    : pimavanserin central (ug)
PIMMET  : 0    : AC-279 metabolite (ug)
LDGUT   : 0    : levodopa depot (ug)
LDCEN   : 0    : levodopa central (ug)
LDBRN   : 0    : brain L-DOPA (ug)
APDGUT  : 0    : antipsychotic depot (ug)
APDCEN  : 0    : antipsychotic central (ug)
ZONGUT  : 0    : zonisamide depot (ug)
ZONCEN  : 0    : zonisamide central (ug)
AMBGUT  : 0    : ambroxol depot (ug)
AMBCEN  : 0    : ambroxol central (ug)
AMBBRN  : 0    : brain ambroxol (ug)

// ---- mAb (23-25)
MABCEN  : 0    : anti-a-syn mAb central (mg)
MABPER  : 0    : anti-a-syn mAb peripheral (mg)
MABCSF  : 0    : anti-a-syn mAb CSF (mg)

// ---- enzyme / molecular states (26-27)
CARBA   : 0    : carbamylated fraction of AChE (0-1)
CARBB   : 0    : carbamylated fraction of BuChE (0-1)

// ---- proteostasis (28-37)
GCASE   : 1.00 : GCase activity (fraction of normal)
GLCCER  : 1.00 : glucosylceramide (au)
ASYNM   : 1.00 : a-syn monomer pool (au)
OLIGB   : 0.18 : oligomer, brainstem (au)
OLIGL   : 0.12 : oligomer, limbic (au)
OLIGN   : 0.06 : oligomer, neocortex (au)
FIBB    : 0.55 : fibril, brainstem (au)
FIBL    : 0.35 : fibril, limbic (au)
FIBN    : 0.18 : fibril, neocortex (au)
SEED    : 0.02 : interstitial seed pool (au)

// ---- neuronal populations (38-46)
NBM     : 0.50 : nucleus basalis cholinergic neurons (fraction)
PPN     : 0.72 : PPN/LDT cholinergic neurons (fraction)
SNC     : 0.50 : SNc dopaminergic neurons (fraction)
LC      : 0.62 : locus coeruleus neurons (fraction)
RAPHE   : 0.88 : dorsal raphe neurons (fraction)
SLD     : 0.55 : REM-atonia circuit (fraction)
ORX     : 0.85 : orexin neurons (fraction)
CSYM    : 0.35 : cardiac sympathetic terminals (fraction)
SYND    : 0.92 : cortical/limbic synaptic density (fraction)

// ---- neurotransmitters / receptors (47-52)
ACHS    : 0.40 : synaptic ACh (au)
DAS     : 0.50 : striatal dopamine (au)   [reset in $MAIN to KDAS*SNC/KDAD]
NES     : 0.70 : cortical noradrenaline (au)
M1R     : 0.95 : postsynaptic M1/M4 density (relative)
D2R     : 1.00 : striatal D2 density (relative)
HT2A    : 1.00 : cortical 5-HT2A density (relative)

// ---- inflammation / co-pathology (53-57)
MGA     : 0.15 : microglial activation (0-1)
MGEXH   : 0.05 : cumulative microglial exhaustion (0-1)
ASTRO   : 0.12 : astrocyte reactivity (0-1)
ABETA   : 18   : amyloid burden (Centiloid)
PTAU    : 0.10 : tau pathology (au)

// ---- clinical states (58-64)
COG     : 27.0 : cognition (MMSE scale)
ATTM    : 0.72 : attention state (bistable variable)
FLUC    : 0.60 : measured fluctuation amplitude (CAF-like)
VHB     : 0.30 : visual hallucination burden (NPI-hall units)
MOT     : 18.0 : MDS-UPDRS III
RBDS    : 4.5  : RBD severity (RBDSQ-like)
AUTS    : 12.0 : orthostatic SBP drop (mmHg)
EDSS    : 7.0  : Epworth sleepiness score
CUMH    : 0    : cumulative hazard

$GLOBAL
#define CRIV   (RIVCEN/VRIV)
#define CDON   (DONCEN/VDON)
#define CMEM   (MEMCEN/VMEM)
#define CPIM   (PIMCEN/VPIM)
#define CMET   (PIMMET/VPIM)
#define CLD    (LDCEN/VLD)
#define CAPD   (APDCEN/VAPD)
#define CZON   (ZONCEN/VZON)
#define CAMB   (AMBCEN/VAMB)
#define CAMBB  (AMBBRN/VAMB)
#define CMAB   (MABCEN/V1MAB)
#define CMABC  (MABCSF/VCSF)

// derived quantities shared between $ODE and $TABLE
double INHACHE, INHBCHE, HYDCAP, CHATCAP, ACHSIG, CHDRIVE;
double OCCD2, RESERVE, DEFICIT, NSENS, DADRIVE, DACONT, VMATF, LDCONV, AADCF;
double OCC2A, HT2ASIG, HTTONE;
double NEDRIVE, THALD, DRIVEA, CURV, FLUCTGT, DCUB, DCRITB, BISTAB, NOISEG, AROUSD;
double BOTTOMUP, TOPDOWN, VHDRIVE, MESO;
double SEEDTOT, FIBTOT, OLIGTOT, PHAGO, INFL, TOXB, TOXL, TOXN;
double GCTGT, EAMB, KOLDE, XSEED;
double UPCAP, D2TGT, M1TGT, H2TGT;
double COGTGT, MOTTGT, RBDTGT, AUTTGT, EDSTGT;
double MMSE, NPI, CAF, MIBG, DATSBR, EEGF, HAZ, FALLS, QTC, GIAE;
double UPTB, UPTL, UPTN, GLYMF;
double TROPNX, TAUPX, CONLX, STRAINX;

$MAIN
// ---- phenotype switches.  These set THREE numbers and nothing else. -------
TROPNX  = TROPN;
TAUPX   = TAUPHENO;
CONLX   = CONL;
STRAINX = STRAIN;
if(PHENO == 1) {            // PD dementia: brainstem-first, low limbic/cortical
  TROPNX  = 0.42;           //   a-syn early, tau low
  CONLX   = 0.34;
  TAUPX   = 0.70;
  STRAINX = 0.88;
}
if(PHENO == 2) {            // AD: little a-syn, heavy tau
  TROPNX  = 0.10;
  CONLX   = 0.10;
  TAUPX   = 2.40;
  STRAINX = 0.42;
}

// The phenotype is not a label bolted onto the outputs: it is the initial
// REGIONAL DISTRIBUTION of a-synuclein plus the tau/amyloid load.  Everything
// downstream -- levodopa response, neuroleptic sensitivity, ChEI efficacy --
// is computed from these numbers by equations that never inspect PHENO again.
if(PHENO == 1) {            // PD-dementia: brainstem-heavy, limbic-light
  FIBB_0 = 0.70;  FIBL_0 = 0.12;  FIBN_0 = 0.03;
  OLIGB_0 = 0.24; OLIGL_0 = 0.04; OLIGN_0 = 0.010;
  SNC_0  = 0.30;  CSYM_0 = 0.30;  SLD_0  = 0.45;  NBM_0 = 0.62;
  LC_0   = 0.50;  PPN_0  = 0.62;
  PTAU_0 = 0.08;  ABETA_0 = 14;
  MOT_0  = 22.0;
}
if(PHENO == 2) {            // AD: a-syn near-absent, tau and amyloid heavy
  FIBB_0 = 0.08;  FIBL_0 = 0.04;  FIBN_0 = 0.02;
  OLIGB_0 = 0.03; OLIGL_0 = 0.012; OLIGN_0 = 0.006;
  SNC_0  = 0.85;  CSYM_0 = 0.85;  SLD_0  = 0.90;  NBM_0 = 0.72;
  LC_0   = 0.80;  PPN_0  = 0.94;
  PTAU_0 = 0.55;  ABETA_0 = 70;
  MOT_0  = 3.0;   RBDS_0 = 1.0;   AUTS_0 = 6.0;  EDSS_0 = 6.0;
}

// ------------------------------------------------------------------ ETAs
// Between-subject variability enters where it is actually observed: how much
// of each vulnerable population a patient still has at diagnosis, and how much
// limbic/neocortical Lewy burden they carry.  Everything else -- every rate
// constant, every receptor equation -- is shared.  That is deliberate: it means
// the responder/non-responder split the model produces is a consequence of
// baseline PATHOLOGY, not of a fitted responder covariate.
NBM_0   = fmin(NBM_0  *exp(ETA_NBM), 0.98);
PPN_0   = fmin(PPN_0  *exp(ETA_NBM), 0.98);
SNC_0   = fmin(SNC_0  *exp(ETA_SNC), 0.98);
LC_0    = fmin(LC_0   *exp(ETA_SNC), 0.98);
FIBL_0  = fmin(FIBL_0 *exp(ETA_FIB), 0.95);
FIBN_0  = fmin(FIBN_0 *exp(ETA_FIB), 0.95);
OLIGL_0 = fmin(OLIGL_0*exp(ETA_FIB), 0.95);
OLIGN_0 = fmin(OLIGN_0*exp(ETA_FIB), 0.95);
ABETA_0 = ABETA_0*exp(ETA_AB);

// Initialise the two FAST cholinergic states at the quasi-steady state
// implied by NBM_0 and PTAU_0, so that no scenario begins with a spurious
// several-month transient in cognition.
M1R_0  = fmax(M1BASE - KTAUM1*PTAU_0, 0.05);
ACHS_0 = (-1.0 + sqrt(1.0 + 4.0*KM2*NBM_0))/(2.0*KM2);

// ...and put the attention state, the fluctuation amplitude and the cognitive
// score on their own equilibria too, so that a scenario's first six months
// report disease progression rather than the model settling in.
DAS_0  = KDAS*SNC_0/KDAD;
NES_0  = KNES*LC_0/KNED;
double CH0  = ACHS_0*M1R_0*(1.0 - 0.18*ANTICH);
double DRV0 = WACH*CH0 + WNE*LC_0 + WTHAL*PPN_0
              - WEDS*(EDSS_0/24.0) - WOH*(AUTS_0/40.0);
double D0   = DRV0 - ATTOFF;
// Seed Newton on the UPPER branch unconditionally.  A patient arriving at
// diagnosis has descended from health and has not yet fallen off; which branch
// they are on is a matter of HISTORY, not of the current drive, and that is
// precisely what bistability means.
double A0   = 1.3;
for(int i = 0; i < 80; i++) {          // Newton on -A^3 + ALPHAB*A + D0 = 0
  double f  = -A0*A0*A0 + ALPHAB*A0 + D0;
  double fp = -3.0*A0*A0 + ALPHAB;
  if(fabs(fp) < 1e-6) break;
  A0 = A0 - f/fp;
}
ATTM_0 = A0;

double DC0  = 2.0*pow(ALPHAB/3.0, 1.5);
double BS0  = fmax(1.0 - fabs(D0)/DC0, 0.0);
double AD0  = fmin(fmax(1.0 - (0.5*LC_0 + 0.5*PPN_0)/AROUSREF, 0.0), 1.0);
double NG0  = ANOISE0 + KNOISE*AD0;
FLUC_0 = SIGN0*(FLBASE + BS0*NG0);

COG_0  = MMSE0 - WSYN*(1.0 - SYND_0) - WATTC*fmax(1.0 - ATTM_0/0.95, 0.0)
                - WFLUCC*FLUC_0 - 3.0*PTAU_0
                - WACHC*fmax(CHREF - CH0, 0.0);
double DSG0 = pow(fmax(DAS_0,0.0), TERMEXP)*D2R_0*exp(-KSUPP*OLIGL_0);
MOT_0  = MOTAGE + MOTMAX*pow(K50M, HMOT)/(pow(K50M, HMOT) + pow(fmax(DSG0,1e-9), HMOT))
         + MOTNDA*(FIBN_0 + 0.5*FIBL_0);
RBDS_0 = 13.0*pow(fmax(1.0 - SLD_0, 0.0), 1.3);
AUTS_0 = 5.0 + 30.0*pow(fmax(1.0 - CSYM_0, 0.0), 3.0);
EDSS_0 = 6.0 + 16.0*pow(fmax(1.0 - ORX_0, 0.0), 2.0);

$ODE
// =========================================================================
//  0.  AGGREGATE POOLS
// =========================================================================
FIBTOT  = FIBB + FIBL + FIBN;
OLIGTOT = OLIGB + OLIGL + OLIGN;

// =========================================================================
//  1.  GBA1 BISTABLE SWITCH
//      GCase trafficking is inhibited by oligomer (Hill 2 -> fold),
//      GlcCer accumulates when GCase falls, and GlcCer in turn slows
//      oligomer degradation.  Loop gain > 1 over part of the range.
// =========================================================================
EAMB  = EMAXAMB * CAMBB / (EC50AMB + CAMBB);
GCTGT = GCMAX * GBAF * (1.0 + EAMB) /
        (1.0 + pow(fmax(OLIGTOT,0.0)/KIOL, HOL));
dxdt_GCASE  = KGCIN * (GCTGT - GCASE);
dxdt_GLCCER = KGLSYN - KGLDEG * fmax(GCASE,1e-6) * GLCCER - KGLBAS * GLCCER;

KOLDE = KOLDEG / (1.0 + fmax(GLCCER,0.0)/KGLC50);   // GlcCer-slowed degradation

// =========================================================================
//  2.  ALPHA-SYNUCLEIN: monomer, three regional oligomer/fibril pairs,
//      one shared interstitial seed pool.  Staging is produced by gating
//      each region's uptake on the burden of the region upstream of it.
// =========================================================================
XSEED = 1.0 + KXSEED * (ABETA/(AB50 + ABETA) + 0.5*PTAU);

GLYMF = 0.55 + 0.45 * fmin(fmax(SLD,0.0),1.0);      // sleep quality -> washout
SEEDTOT = fmax(SEED, 0.0);

UPTB = KUPT * SEEDTOT * CONB;
UPTL = KUPT * SEEDTOT * CONLX * (FIBB/(KSTG + FIBB));
UPTN = KUPT * SEEDTOT * CONN * TROPNX * (FIBL/(KSTG + FIBL));

dxdt_ASYNM = KSYNAS * SNCADOSE
             - KDEGAS * ASYNM * (0.35 + 0.65*GCASE)
             - KNUC1 * ASYNM * ASYNM
             - KNUC2 * STRAINX * ASYNM * FIBTOT;

// Each region carries an oligomer and a fibril pool, both expressed as a
// fraction of the maximum burden that region can hold.  The (1 - X) capacity
// factors are what make the trajectory sigmoid rather than exponential: an
// unbounded autocatalytic loop is not a disease model, it is a divergence.
double CLRF = KFCLR*(0.4 + 0.6*PHAGO);

// brainstem
dxdt_OLIGB = (KNUC1*ASYNM*ASYNM*0.60 + KNUC2*STRAINX*ASYNM*FIBB + UPTB)*(1.0 - OLIGB)
             - (KELONG + KOLDE)*OLIGB;
dxdt_FIBB  = KELONG*OLIGB*XSEED*(1.0 - FIBB) - CLRF*FIBB;

// limbic
dxdt_OLIGL = (KNUC1*ASYNM*ASYNM*0.25 + KNUC2*STRAINX*ASYNM*FIBL + UPTL)*(1.0 - OLIGL)
             - (KELONG + KOLDE)*OLIGL;
dxdt_FIBL  = KELONG*OLIGL*XSEED*(1.0 - FIBL) - CLRF*FIBL;

// neocortex
dxdt_OLIGN = (KNUC1*ASYNM*ASYNM*0.15 + KNUC2*STRAINX*ASYNM*FIBN + UPTN)*(1.0 - OLIGN)
             - (KELONG + KOLDE)*OLIGN;
dxdt_FIBN  = KELONG*OLIGN*XSEED*(1.0 - FIBN) - CLRF*FIBN;

// interstitial seed: released by fibrils, removed by glymphatics, microglia,
// antibody, and uptake into recipient neurons
dxdt_SEED = KREL*FIBTOT
            - KGLYM*SEED*GLYMF
            - KSEEDEL*SEED
            - KPHAGO*SEED*PHAGO
            - KBINDMAB*CMABC*SEED
            - (UPTB + UPTL + UPTN);

// =========================================================================
//  3.  NEUROINFLAMMATION.  Phagocytic capacity is protective early and
//      exhausted late, so the same microglial activation changes sign.
// =========================================================================
dxdt_MGA   = KMGON*SEEDTOT*(1.0 - MGA) - KMGOFF*MGA;
dxdt_MGEXH = 0.010*MGA*(1.0 - MGEXH) - 0.0015*MGEXH;
PHAGO      = fmax(MGA*(1.0 - KEXH*MGEXH), 0.0) * (1.0 + KOPSON*CMABC);
dxdt_ASTRO = KASTON*(MGA + 0.5*OLIGTOT)*(1.0 - ASTRO) - KASTOFF*ASTRO;
INFL       = KINFL*(0.6*MGA + 0.4*ASTRO);

// =========================================================================
//  4.  AD CO-PATHOLOGY (present in most DLB; the amount is what sets the
//      phenotype, through tau -> M1 and tau -> synapse)
// =========================================================================
dxdt_ABETA = KABPROD*(1.0 + KABE4*APOE4)*(1.0 - ABETA/ABMAX);
dxdt_PTAU  = KTAU*TAUPX*(1.0 + KTAUAB*ABETA/(AB50 + ABETA))*(1.0 + 0.4*FIBN)
             - KTAUOUT*PTAU;

// =========================================================================
//  5.  NEURONAL POPULATIONS.  Each nucleus sees the oligomer burden of the
//      region it lives in, plus the shared inflammatory term.
// =========================================================================
TOXB = OLIGB + INFL;
TOXL = OLIGL + INFL;
TOXN = OLIGN + INFL;

dxdt_NBM   = -KLNBM  * NBM   * (TOXL + KNBMTAU*PTAU);
dxdt_PPN   = -KLPPN  * PPN   * TOXB;
dxdt_SNC   = -KLSNC  * SNC   * TOXB;
dxdt_LC    = -KLLC   * LC    * TOXB;
dxdt_RAPHE = -KLRAP  * RAPHE * TOXB;
dxdt_SLD   = -KLSLD  * SLD   * TOXB;
dxdt_ORX   = -KLORX  * ORX   * TOXB;
dxdt_CSYM  = -KLCSYM * CSYM  * TOXB;

// synaptic density: a-syn oligomer + tau + excitotoxicity, the last of which
// is where memantine acts
double EXCITO = (1.0 - EMAXMEM*CMEM/(IC50MEM + CMEM));
dxdt_SYND  = -KLSYN * SYND * (TOXN + TOXL + KTAUN*PTAU) * EXCITO;

// =========================================================================
//  6.  CHOLINERGIC TRANSDUCER  (presynaptic down, postsynaptic PRESERVED)
// =========================================================================
// -- carbamylated-enzyme states: the PD half-life is enzyme resynthesis,
//    not drug clearance.  This is why q12h rivastigmine works at all.
dxdt_CARBA = KCARB*CRIV*(1.0 - CARBA) - KDECARB*CARBA;
dxdt_CARBB = RIVBCH*KCARB*CRIV*(1.0 - CARBB) - KDECARB*CARBB;

double OCCDON = CDON > 0 ? DONEFF/(IC50DON + DONEFF) : DONEFF/(IC50DON + DONEFF);
double OCCGAL = GALON*GALINH;
INHACHE  = 1.0 - (1.0 - CARBA)*(1.0 - OCCDON)*(1.0 - OCCGAL);
INHBCHE  = CARBB;

// BuChE takes over as AChE and cholinergic terminals are lost
double BSHARE = BCHEFR + KBCHUP*(1.0 - NBM);
HYDCAP   = (1.0 - BSHARE)*(1.0 - INHACHE) + BSHARE*(1.0 - INHBCHE);

CHATCAP  = NBM;
dxdt_ACHS = KACHS*CHATCAP/(1.0 + KM2*ACHS) - KACHD*ACHS*HYDCAP;

// postsynaptic M1/M4: target falls ONLY with tau.  In DLB tau is low, so M1
// survives; in AD tau is high, so it does not.  Same equation, both arms.
M1TGT = fmax(M1BASE - KTAUM1*PTAU, 0.05);
dxdt_M1R = KM1R*(M1TGT - M1R);

ACHSIG   = ACHS * M1R * (1.0 - 0.18*ANTICH);
CHDRIVE  = fmax(ACHSIG, 0.0);
THALD    = PPN * (1.0 + GALON*GALNIC) * (1.0 - 0.18*ANTICH);

// =========================================================================
//  7.  DOPAMINERGIC TRANSDUCER  (presynaptic down, postsynaptic NOT up)
// =========================================================================
// Striatal dopamine CONTENT falls faster than SNc cell COUNT, because the
// axon terminals degenerate first.  TERMEXP is the same exponent that maps SNc
// survival onto the DaTSCAN binding ratio, and using it here is what gives the
// striatum a genuine motor RESERVE: a healthy striatum can lose 60-70% of its
// dopaminergic signal to a D2 blocker and stay nearly asymptomatic, which is
// why the same risperidone dose is tolerable in AD and catastrophic in a Lewy
// body disease.  Nothing about that contrast is written down; it is the
// position of each phenotype on this one sigmoid.
VMATF   = 0.25 + 0.75*SNC;                        // vesicular buffering left
// L-DOPA must be decarboxylated by AADC, most of which lives in the very
// dopaminergic terminals that have been lost, and stored by VMAT2, likewise.
// Both factors are presynaptic, so the same oral dose delivers less usable
// dopamine the more denervated the striatum is.
AADCF   = AADC0 + (1.0 - AADC0)*fmax(SNC, 0.0);
LDCONV  = KLD2DA*(LDBRN/VLD)*AADCF*(1.0 - KVMATW*(1.0 - VMATF));
double ZONF = EMAXZON*CZON/(EC50ZON + CZON);

dxdt_DAS = KDAS*SNC*(1.0 + ZONF) + LDCONV - KDAD*DAS;

// --- the single most consequential equation in the file --------------------
// Capacity for denervation supersensitivity is destroyed by limbic a-syn,
// because the postsynaptic striatal/limbic neuron is itself diseased in DLB
// and is not in early PD.  This one term produces BOTH the blunted levodopa
// response AND the neuroleptic sensitivity, without either being coded.
// POSTI is the integrity of the postsynaptic striatal/limbic signal-
// transduction apparatus.  It is destroyed by limbic a-syn, which is high in
// DLB and low in early PD.  It appears in BOTH the levodopa response and the
// residual-signal calculation, so ONE term produces the blunted L-DOPA effect
// AND the neuroleptic sensitivity.  Neither is coded anywhere.
UPCAP = exp(-KSUPP*OLIGL);
D2TGT = 1.0 + KD2UP*fmax(1.0 - DAS/0.95, 0.0)*UPCAP;
dxdt_D2R = KD2R*(D2TGT - D2R);

OCCD2   = CAPD/(EC50D2 + CAPD);
DACONT  = pow(fmax(DAS, 0.0), TERMEXP);           // striatal DA content
DADRIVE = DACONT * D2R * UPCAP * (1.0 - OCCD2);

// RESERVE is what the postsynaptic side can still deliver WITHOUT any drug on
// board.  DEFICIT is the D2 blockade the drug imposes, AMPLIFIED by how little
// reserve there was to absorb it.  A patient with intact reserve tolerates the
// same occupancy; a DLB patient does not.  Note that with no antipsychotic on
// board OCCD2 = 0, so DEFICIT = 0 and NSENS = 0 no matter how advanced the
// disease is — the reaction is a DRUG event, not a severity milestone.
RESERVE = DACONT*D2R*UPCAP;                       // 1.0 = healthy young adult
DEFICIT = OCCD2*(1.0 + KRES*(1.0 - fmin(fmax(RESERVE,0.0),1.0)));
NSENS   = 1.0/(1.0 + exp(-(DEFICIT - DCRIT)/DWIDTH));

// =========================================================================
//  8.  SEROTONERGIC TRANSDUCER  (presynaptic down, postsynaptic UP)
// =========================================================================
HTTONE = RAPHE;
H2TGT  = 1.0 + KH2UP*fmax(1.0 - RAPHE, 0.0) + KH2NEO*FIBN;
dxdt_HT2A = KHT2A*(H2TGT - HT2A);

double PIMEQ = CPIM + POTMET*CMET;
OCC2A   = PIMEQ/(EC502A + PIMEQ);
HT2ASIG = HT2A*(HTCONST + HTTONE)*(1.0 - OCC2A);

// =========================================================================
//  9.  NORADRENERGIC
// =========================================================================
dxdt_NES = KNES*LC - KNED*NES;
NEDRIVE  = NES;

// =========================================================================
// 10.  ATTENTION AS A BISTABLE STATE — fluctuation is its VARIANCE
// =========================================================================
DRIVEA = WACH*CHDRIVE + WNE*NEDRIVE + WTHAL*THALD
         - WEDS*(EDSS/24.0) - WOH*(AUTS/40.0);

dxdt_ATTM = KATT*( -ATTM*ATTM*ATTM + ALPHAB*ATTM + DRIVEA - ATTOFF );

// What the clinic scores as "fluctuation" is SWITCHING between attentional
// states, and switching requires two states to exist.  The right measure is
// therefore not the curvature AT the fixed point but how deep the system sits
// inside the BISTABLE REGION -- which, for this cubic, is the band
// |DRIVEA - ATTOFF| < 2*(ALPHAB/3)^1.5.  That band is ~0.157 wide in drive
// units and the patient takes about two years to traverse it, which is why
// fluctuation is a feature of mild-to-moderate DLB and fades again once the
// state has dropped irreversibly onto the lower branch.
DCUB   = DRIVEA - ATTOFF;
DCRITB = 2.0*pow(ALPHAB/3.0, 1.5);
BISTAB = fmax(1.0 - fabs(DCUB)/DCRITB, 0.0);

// Bistability supplies the OPPORTUNITY to switch; the noise that actually
// drives the switching comes from the brainstem arousal systems (LC + PPN),
// which a-synuclein destroys in DLB and tau largely spares in AD.  Both
// factors are needed, which is why fluctuation is a DLB core feature and only
// a minor one in AD even though both diseases traverse the same cubic.
AROUSD  = fmin(fmax(1.0 - (0.5*LC + 0.5*PPN)/AROUSREF, 0.0), 1.0);
NOISEG  = ANOISE0 + KNOISE*AROUSD;
CURV    = fmax(3.0*ATTM*ATTM - ALPHAB, CURVMIN);   // reported, not used
FLUCTGT = SIGN0*(FLBASE + BISTAB*NOISEG);
dxdt_FLUC = KFLUC*(FLUCTGT - FLUC);

// =========================================================================
// 11.  VISUAL HALLUCINATIONS — a PRODUCT of three deficits (PAD model)
// =========================================================================
BOTTOMUP = fmin(fmax(1.0 - 1.15*FIBN - 0.35*(1.0 - SYND), 0.0), 1.0);
TOPDOWN  = fmin(fmax(0.25 + 0.75*fmin(fmax(ATTM,0.0),1.2)/1.2, 0.0), 1.0);
MESO     = KMESO*(LDBRN/VLD)/(120.0 + LDBRN/VLD);

// The 5-HT2A term sits on a floor, because inverse agonism at 5-HT2A is not
// the whole of the psychosis mechanism.  Without the floor the model would
// predict near-abolition of hallucinations by pimavanserin, which is not what
// the trials show (SAPS-PD -3.06, about a third).
VHDRIVE = (WBU*(1.0 - BOTTOMUP)) * (WTD*(1.0 - TOPDOWN)) * (VHFLOOR + HT2ASIG)
          + MESO;
dxdt_VHB = KVH*(VHSCALE*VHDRIVE - VHB);

// =========================================================================
// 12.  CLINICAL ENDPOINT INTEGRATION
// =========================================================================
// Cognition has FOUR separable debits: structural synapse loss, the
// attentional state, the variance of that state, and tau.  The cholinergic
// term is separate from the attentional one because ACh supports encoding
// independently of the bistable switch -- and it is the only one a
// cholinesterase inhibitor can move within weeks.
COGTGT = MMSE0 - WSYN*(1.0 - SYND) - WATTC*fmax(1.0 - ATTM/0.95, 0.0)
                - WFLUCC*FLUC - 3.0*PTAU
                - WACHC*fmax(CHREF - CHDRIVE, 0.0);
dxdt_COG = KCOG*(COGTGT - COG);

// Three additive motor sources, only ONE of which levodopa can touch: the
// dopaminergic sigmoid, a non-dopaminergic axial/cortical burden that scales
// with neocortical and limbic fibril load, and the neuroleptic-sensitivity
// surge.  The second term is why the same levodopa exposure buys much less in
// DLB than in PD-dementia -- most of a DLB patient's parkinsonism is not
// dopaminergic in the first place.
MOTTGT = MOTAGE
         + MOTMAX*pow(K50M, HMOT)
           /(pow(K50M, HMOT) + pow(fmax(DADRIVE, 1e-9), HMOT))
         + MOTNDA*(FIBN + 0.5*FIBL)
         + 18.0*NSENS;
dxdt_MOT = KMOT*(MOTTGT - MOT);

RBDTGT = 13.0*pow(fmax(1.0 - SLD,0.0), 1.3)*(1.0 - MELON*MELEFF);
dxdt_RBDS = KRBD*(RBDTGT - RBDS);

AUTTGT = (5.0 + 30.0*pow(fmax(1.0 - CSYM,0.0), 3.0))*(1.0 - DROXON*DROXEFF)
         + APDA1*8.0*CAPD/(EC50D2 + CAPD);
dxdt_AUTS = KAUT*(AUTTGT - AUTS);

EDSTGT = (6.0 + 16.0*pow(fmax(1.0 - ORX,0.0), 2.0))*(1.0 - MODON*MODEFF)
         + APDSED*6.0*CAPD/(EC50D2 + CAPD);
dxdt_EDSS = KEDS*(EDSTGT - EDSS);

// =========================================================================
// 13.  SURVIVAL.  Neuroleptic sensitivity enters MULTIPLICATIVELY, which is
//      what reproduces the ~2-3x mortality after antipsychotic exposure.
// =========================================================================
FALLS = 0.35*(MOT/40.0) + 0.30*(AUTS/40.0) + 0.20*FLUC + 0.15*(RBDS/13.0)
        + 0.20*DROXON*DROXSUP;
HAZ   = H0*exp(BCOG*(MMSE0 - COG) + BMOT*MOT + BFALL*FALLS)*exp(BNS*NSENS);
dxdt_CUMH = HAZ;

// =========================================================================
// 14.  DRUG PK
// =========================================================================
dxdt_RIVGUT = -KARIV*RIVGUT;
dxdt_RIVDEP = -1.0*RIVDEP;
dxdt_RIVCEN =  KARIV*FRIV*RIVGUT + 1.0*FPATCH*RIVDEP - (CLRIV/VRIV)*RIVCEN;

dxdt_DONGUT = -KADON*DONGUT;
dxdt_DONCEN =  KADON*DONGUT - (CLDON/VDON)*DONCEN
               - (QDON/VDON)*DONCEN + (QDON/VDON2)*DONPER;
dxdt_DONPER =  (QDON/VDON)*DONCEN - (QDON/VDON2)*DONPER;
dxdt_DONEFF =  KEODON*(CDON - DONEFF);

dxdt_MEMGUT = -KAMEM*MEMGUT;
dxdt_MEMCEN =  KAMEM*MEMGUT - (CLMEM/VMEM)*MEMCEN;

dxdt_PIMGUT = -KAPIM*PIMGUT;
dxdt_PIMCEN =  KAPIM*PIMGUT - (CLPIM/VPIM)*PIMCEN - KMETPIM*PIMCEN;
dxdt_PIMMET =  KMETPIM*PIMCEN - KOUTPIM*PIMMET;

dxdt_LDGUT  = -KALD*LDGUT;
dxdt_LDCEN  =  KALD*FLD*LDGUT - (CLLD/VLD)*LDCEN - KINLD*LDCEN;
dxdt_LDBRN  =  KINLD*LDCEN - KOUTLD*LDBRN;

dxdt_APDGUT = -KAAPD*APDGUT;
dxdt_APDCEN =  KAAPD*APDGUT - (CLAPD/VAPD)*APDCEN;

dxdt_ZONGUT = -KAZON*ZONGUT;
dxdt_ZONCEN =  KAZON*ZONGUT - (CLZON/VZON)*ZONCEN;

dxdt_AMBGUT = -KAAMB*AMBGUT;
dxdt_AMBCEN =  KAAMB*AMBGUT - (CLAMB/VAMB)*AMBCEN - KINAMB*AMBCEN;
dxdt_AMBBRN =  KINAMB*AMBCEN - KOUTAMB*AMBBRN;

dxdt_MABCEN = -(CLMAB/V1MAB)*MABCEN - (QMAB/V1MAB)*MABCEN + (QMAB/V2MAB)*MABPER;
dxdt_MABPER =  (QMAB/V1MAB)*MABCEN - (QMAB/V2MAB)*MABPER;
// CSF concentration is a PARTITION of plasma, not a mass sink on it: at steady
// state CMABC/CMAB = KCSFIN/KOUTCSF = 0.3%, which is what CSF sampling in the
// anti-a-syn trials actually shows.
dxdt_MABCSF =  KCSFIN*(MABCEN/V1MAB)*VCSF - KOUTCSF*MABCSF;

$TABLE
// ---- clinical scales -------------------------------------------------------
MMSE   = fmin(fmax(COG, 0.0), 30.0);
CAF    = fmin(fmax(4.0*FLUC, 0.0), 16.0);
NPI    = NPIBASE + VHB + 2.2*FLUC + 1.6*(RBDS/13.0)*3.0;
QTC    = QTSLOPE*(CPIM + POTMET*CMET);
GIAE   = 100.0*(CRIV/(CRIV + 14.0));           // Cmax-driven, hence patch < oral

// ---- imaging / fluid biomarkers -------------------------------------------
MIBG   = 1.05 + 1.75*fmax(CSYM,0.0);           // heart:mediastinum ratio
DATSBR = fmax(0.05, 0.05 + 0.95*pow(fmax(SNC,0.0), 1.6));  // terminals go first
EEGF   = 5.6 + 3.9*fmin(fmax(THALD*(0.5 + 0.5*ATTM),0.0),1.0);  // dominant EEG freq (Hz)

// ---- capture-only conveniences --------------------------------------------
double SURV = exp(-CUMH);
double AGE  = AGE0 + TIME/365.25;
double PCTLD = 0.0;

$CAPTURE @annotated
CRIV : rivastigmine plasma (ng/mL)
CDON : donepezil plasma (ng/mL)
CPIM : pimavanserin plasma (ng/mL)
CLD  : levodopa plasma (ng/mL)
CAPD : antipsychotic plasma (ng/mL)
CAMBB: brain ambroxol (ng/mL equiv)
CMAB : anti-a-syn mAb plasma (mg/L)
CMABC: CSF mAb (mg/L)
INHACHE : fractional brain AChE inhibition
INHBCHE : fractional BuChE inhibition
CHDRIVE : cholinergic transducer output
DADRIVE : dopaminergic transducer output
DACONT  : striatal dopamine content (terminal-weighted)
HT2ASIG : serotonergic (5-HT2A) transducer output
OCCD2   : striatal D2 occupancy
OCC2A   : 5-HT2A inverse-agonist occupancy
RESERVE : untreated postsynaptic D2 reserve
DEFICIT : drug-induced D2 deficit, reserve-amplified
NSENS   : neuroleptic-sensitivity index (0-1)
UPCAP   : capacity for D2 denervation supersensitivity
DRIVEA  : total attention drive (cholinergic + noradrenergic + thalamic)
CURV    : restoring stiffness of the attention state (diagnostic)
BISTAB  : depth inside the bistable attention region (0-1)
NOISEG  : arousal-instability gain on the fluctuation amplitude
FLUCTGT : target fluctuation amplitude
BOTTOMUP: bottom-up visual evidence fidelity
TOPDOWN : top-down attentional binding
VHDRIVE : hallucination drive (the triple product)
FIBTOT  : total fibril burden
OLIGTOT : total oligomer burden
PHAGO   : effective microglial phagocytic capacity
GCTGT   : GCase trafficking target
EAMB    : fractional GCase enhancement by ambroxol
MMSE    : MMSE
CAF     : clinician assessment of fluctuation
NPI     : NPI total
MIBG    : MIBG heart:mediastinum ratio
DATSBR  : DaTSCAN striatal binding ratio (relative)
EEGF    : dominant EEG frequency (Hz)
QTC     : QTc change (ms)
GIAE    : cholinergic GI adverse-effect index
FALLS   : falls index
HAZ     : instantaneous hazard
SURV    : survival probability
AGE     : age (y)

$OMEGA @annotated @block
ETA_NBM  : 0.045 : between-subject variability, baseline nbM survival
ETA_SNC  : 0.020 0.060 : between-subject variability, baseline SNc survival
ETA_FIB  : 0.010 0.018 0.110 : between-subject variability, limbic/neocortical burden

$OMEGA @annotated
ETA_AB   : 0.120 : between-subject variability, baseline amyloid burden

$SIGMA 0
)---"

## -----------------------------------------------------------------------------
##  BUILD
## -----------------------------------------------------------------------------
mod <- mcode_cache("dlb_qsp", dlb_code, atol = 1e-8, rtol = 1e-8, maxsteps = 500000)

## =============================================================================
##  DOSING HELPERS
## =============================================================================

## Rivastigmine oral capsule: amount in micrograms, BID
riv_oral  <- function(mg = 6,    start = 180, dur = 2740) ev(amt = mg*1000, cmt = "RIVGUT", time = start, ii = 0.5,   addl = dur/0.5 - 1)

## Rivastigmine transdermal: a daily depot load released with a 1/day rate
## constant, which gives an essentially flat plasma profile.
riv_patch <- function(mg24 = 9.5, start = 180, dur = 2740) ev(amt = mg24*1000, cmt = "RIVDEP", time = start, ii = 1, addl = dur - 1)

don_oral  <- function(mg = 10,  start = 180, dur = 2740) ev(amt = mg*1000, cmt = "DONGUT", time = start, ii = 1, addl = dur - 1)
mem_oral  <- function(mg = 20,  start = 180, dur = 2740) ev(amt = mg*1000, cmt = "MEMGUT", time = start, ii = 1, addl = dur - 1)
pim_oral  <- function(mg = 34,  start = 365, dur = 2555) ev(amt = mg*1000, cmt = "PIMGUT", time = start, ii = 1, addl = dur - 1)
zon_oral  <- function(mg = 25,  start = 365, dur = 2555) ev(amt = mg*1000, cmt = "ZONGUT", time = start, ii = 1, addl = dur - 1)

## Levodopa/carbidopa, three times daily
ld_oral   <- function(mg = 150, start = 365, dur = 2555) ev(amt = mg*1000, cmt = "LDGUT",  time = start, ii = 1/3, addl = dur*3 - 1)

## Generic antipsychotic; PK parameters are set per scenario via param()
apd_oral  <- function(mg = 1,   start = 730, dur = 180, ii = 1) ev(amt = mg*1000, cmt = "APDGUT", time = start, ii = ii, addl = dur/ii - 1)

## Ambroxol 420 mg TID = 1.26 g/day
amb_oral  <- function(mg = 420, start = 0,   dur = 2920) ev(amt = mg*1000, cmt = "AMBGUT", time = start, ii = 1/3, addl = dur*3 - 1)

## Anti-alpha-synuclein mAb: 4500 mg IV every 4 weeks (amounts in mg)
mab_iv    <- function(mg = 4500, start = 180, n = 98) ev(amt = mg, cmt = "MABCEN", time = start, ii = 28, addl = n - 1)

## =============================================================================
##  TWENTY-THREE TREATMENT SCENARIOS
##
##  The suite is built so that every claim has its own matched control.  The
##  three phenotype arms (1 / 19 / 21-style) are NOT separate models: they are
##  the same 64 equations with a different initial REGIONAL DISTRIBUTION of
##  a-synuclein, tau and amyloid.  Scenarios 8-10 and 11-12 exist purely so
##  that "neuroleptic sensitivity" and "blunted levodopa response" can be read
##  off as differences rather than taken on trust.
## =============================================================================

## NOTE: zero_re().  The scenarios are DETERMINISTIC comparisons -- the whole
## point of, say, 12 vs 13 vs 14 is that nothing differs except the phenotype.
## Random effects would make the contrast unreproducible and, worse, would let
## a lucky draw masquerade as a mechanism.  Between-subject variability is used
## deliberately and only in responder_rates(), where a DISTRIBUTION is the
## quantity of interest.
run_scn <- function(label, dose = NULL, param = list(), end = 2920, delta = 1) {
  m <- zero_re(mod)
  if (length(param)) m <- param(m, param)
  out <- if (is.null(dose)) mrgsim_df(m, end = end, delta = delta)
         else mrgsim_df(m, events = dose, end = end, delta = delta)
  out$scenario <- label
  out
}

scenarios <- list(

  ## --------------------------------------------------- natural history arms
  s01 = function() run_scn("01 Natural history (DLB, untreated)"),
  s02 = function() run_scn("02 Natural history (PD-dementia phenotype)",
                           param = list(PHENO = 1)),
  s03 = function() run_scn("03 Natural history (AD phenotype)",
                           param = list(PHENO = 2)),

  ## -------------------------------------------------- cholinergic therapy
  s04 = function() run_scn("04 Rivastigmine 6 mg BID oral", riv_oral(6)),
  s05 = function() run_scn("05 Rivastigmine 9.5 mg/24h patch", riv_patch(9.5)),
  s06 = function() run_scn("06 Donepezil 10 mg qd", don_oral(10)),
  s07 = function() run_scn("07 Donepezil 10 mg + memantine 20 mg",
                           c(don_oral(10), mem_oral(20))),

  ## ---- the AD comparator for ChEI response.  Same drug, same dose; tau has
  ##      taken the postsynaptic receptor away, so the transducer has less to
  ##      multiply.  Compare 06 - 01 against 08 - 03.
  s08 = function() run_scn("08 Donepezil 10 mg qd (AD phenotype)",
                           don_oral(10), param = list(PHENO = 2)),

  ## ---- the commonest iatrogenic error in the whole disease
  s09 = function() run_scn("09 Donepezil 10 mg + oxybutynin (anticholinergic)",
                           don_oral(10), param = list(ANTICH = 2)),

  ## --------------------------------------------------------- psychosis arms
  s10 = function() run_scn("10 Pimavanserin 34 mg qd", pim_oral(34)),
  s11 = function() run_scn("11 Quetiapine 50 mg qd",
                           apd_oral(50, dur = 1095),
                           param = list(EC50D2 = 900, CLAPD = 1400, VAPD = 700,
                                        APDA1 = 1.0, APDSED = 1.0)),

  ## ---- 12/13/14: IDENTICAL risperidone exposure, three receptor fields.
  ##      Nothing but PHENO differs, and PHENO only sets initial pathology.
  s12 = function() run_scn("12 Risperidone 1 mg qd (DLB)",
                           apd_oral(1, dur = 180)),
  s13 = function() run_scn("13 Risperidone 1 mg qd (PD-dementia phenotype)",
                           apd_oral(1, dur = 180), param = list(PHENO = 1)),
  s14 = function() run_scn("14 Risperidone 1 mg qd (AD phenotype)",
                           apd_oral(1, dur = 180), param = list(PHENO = 2)),

  ## ------------------------------------------------------------- motor arms
  ##  15/16: identical levodopa exposure, DLB vs PD-dementia.
  s15 = function() run_scn("15 Levodopa/carbidopa 150 mg TID (DLB)", ld_oral(150)),
  s16 = function() run_scn("16 Levodopa/carbidopa 150 mg TID (PD-dementia)",
                           ld_oral(150), param = list(PHENO = 1)),
  s17 = function() run_scn("17 Zonisamide 25 mg qd", zon_oral(25)),

  ## ----------------------------------------------- disease-modification arms
  ##  18/19/20: the GBA1 deadline.  Same genotype, same drug, same exposure;
  ##  only the START TIME differs between 19 and 20.
  s18 = function() run_scn("18 GBA1 carrier (GBAF 0.55), untreated",
                           param = list(GBAF = 0.55)),
  s19 = function() run_scn("19 Ambroxol 1.26 g/d from day 0 (GBA1 carrier)",
                           amb_oral(420, start = 0),
                           param = list(GBAF = 0.55)),
  s20 = function() run_scn("20 Ambroxol 1.26 g/d from day 1095 (GBA1 carrier)",
                           amb_oral(420, start = 1095, dur = 1825),
                           param = list(GBAF = 0.55)),

  ## ---- 20b: the falsification test for the GBA1 loop.  Ambroxol for two
  ##      years, then WITHDRAWN.  If the loop were truly bistable, GCase would
  ##      stay on the rescued branch after withdrawal.  If it is only a strong
  ##      positive-feedback amplifier, it relaxes back.  Run it and see -- the
  ##      answer is reported in the model header, and it is "relaxes back".
  s20b = function() run_scn("20b Ambroxol day 0-730 then WITHDRAWN (GBA1 carrier)",
                            amb_oral(420, start = 0, dur = 730),
                            param = list(GBAF = 0.55)),

  ## ------------------------------------------------- immunotherapy + optimum
  s21 = function() run_scn("21 Anti-a-syn mAb 4500 mg IV q4w", mab_iv(4500)),
  s22 = function() run_scn("22 Optimised: patch + pimavanserin + melatonin + droxidopa",
                           c(riv_patch(9.5), pim_oral(34)),
                           param = list(MELON = 1, DROXON = 1))
)

run_all <- function() do.call(rbind, lapply(scenarios, function(f) f()))

## =============================================================================
##  QUICK-LOOK SUMMARIES
## =============================================================================

summarise_scn <- function(df) {
  do.call(rbind, lapply(split(df, df$scenario), function(d) {
    yr <- function(t) d[which.min(abs(d$time - t)), ]
    data.frame(
      scenario     = d$scenario[1],
      MMSE_0       = round(yr(0)$MMSE, 1),
      MMSE_y1      = round(yr(365)$MMSE, 1),
      MMSE_y3      = round(yr(1095)$MMSE, 1),
      MMSE_y5      = round(yr(1825)$MMSE, 1),
      CAF_y1       = round(yr(365)$CAF, 2),
      NPI_y3       = round(yr(1095)$NPI, 1),
      UPDRS3_y2    = round(yr(730)$MOT, 1),
      VH_y3        = round(yr(1095)$VHB, 2),
      NSENS_max    = round(max(d$NSENS), 3),
      AChEinh_max  = round(max(d$INHACHE), 3),
      FIBN_y5      = round(yr(1825)$FIBN, 3),
      Surv_y5      = round(yr(1825)$SURV, 3),
      row.names    = NULL)
  }))
}

## Median survival implied by the cumulative hazard
median_survival <- function(df) {
  do.call(rbind, lapply(split(df, df$scenario), function(d) {
    i <- which(d$SURV <= 0.5)[1]
    data.frame(scenario = d$scenario[1],
               median_surv_yr = if (is.na(i)) NA_real_ else round(d$time[i]/365.25, 2),
               row.names = NULL)
  }))
}

## The ChEI ratio the model PREDICTS: fractional improvement in fluctuation
## divided by fractional improvement in MMSE, 20 weeks after starting the drug
## (the McKeith 2000 read-out window).
chei_ratio <- function(trt = scenarios$s06, base = scenarios$s01) {
  b <- base(); t <- trt()
  i <- which.min(abs(b$time - (180 + 140)))
  dM <- (t$MMSE[i] - b$MMSE[i]) / b$MMSE[i]
  dF <- (b$CAF[i]  - t$CAF[i])  / b$CAF[i]
  c(dMMSE_abs = t$MMSE[i] - b$MMSE[i], dCAF_abs = b$CAF[i] - t$CAF[i],
    dMMSE_frac = dM, dCAF_frac = dF, ratio = dF/dM)
}

## Same drug, same dose, three receptor fields.  Reads off scenarios 12/13/14.
neuroleptic_demo <- function() {
  r <- lapply(list(scenarios$s12, scenarios$s13, scenarios$s14), function(f) f())
  do.call(rbind, lapply(r, function(d) {
    i <- which.min(abs(d$time - 760))
    data.frame(scenario = d$scenario[1],
               D2occ    = round(d$OCCD2[i], 3),
               reserve  = round(d$RESERVE[i], 3),
               deficit  = round(d$DEFICIT[i], 3),
               NSENS    = round(d$NSENS[i], 3),
               dUPDRS   = round(d$MOT[i] - d$MOT[which.min(abs(d$time - 725))], 1),
               row.names = NULL)
  }))
}

## Same drug, same dose, two presynaptic/postsynaptic configurations.
levodopa_demo <- function() {
  pairs <- list(c("15", "01"), c("16", "02"))
  fns <- list(s15 = scenarios$s15, s16 = scenarios$s16,
              s01 = scenarios$s01, s02 = scenarios$s02)
  d15 <- fns$s15(); d16 <- fns$s16(); d01 <- fns$s01(); d02 <- fns$s02()
  i <- which.min(abs(d15$time - 425))   # two months after starting levodopa
  data.frame(
    phenotype = c("DLB", "PD-dementia"),
    UPDRS_untreated = c(round(d01$MOT[i], 1), round(d02$MOT[i], 1)),
    UPDRS_levodopa  = c(round(d15$MOT[i], 1), round(d16$MOT[i], 1)),
    pct_change = c(round(100*(d15$MOT[i] - d01$MOT[i])/d01$MOT[i], 1),
                   round(100*(d16$MOT[i] - d02$MOT[i])/d02$MOT[i], 1)),
    row.names = NULL)
}

## Responder rates in a virtual population.  The levodopa literature reports a
## RESPONDER RATE (>=20% MDS-UPDRS III improvement), not a mean change, so this
## is the comparison that can actually be held against it.  Between-subject
## variability comes from $OMEGA.
responder_rates <- function(n = 300, thresh = 0.20, tref = 425) {
  out <- lapply(list(list(ph = 0, lab = "DLB"), list(ph = 1, lab = "PD-dementia")),
    function(a) {
      m  <- param(mod, PHENO = a$ph)
      idata <- data.frame(ID = seq_len(n))
      pl <- mrgsim_df(m, idata = idata, end = tref + 5, delta = 5)
      tr <- mrgsim_df(m, idata = idata, events = ld_oral(150),
                      end = tref + 5, delta = 5)
      k  <- which.min(abs(unique(pl$time) - tref)); tk <- unique(pl$time)[k]
      a1 <- pl$MOT[pl$time == tk]; a2 <- tr$MOT[tr$time == tk]
      data.frame(phenotype = a$lab,
                 n = n,
                 mean_pct_change = round(100*mean((a2 - a1)/a1), 1),
                 responder_pct   = round(100*mean((a1 - a2)/a1 >= thresh), 1),
                 row.names = NULL)
    })
  do.call(rbind, out)
}

## The GBA1 lever arm: scenarios 19 and 20 differ only in start time, and 20b
## tests reversibility (i.e. whether the loop is bistable or merely amplifying).
ambroxol_demo <- function() {
  d <- lapply(list(scenarios$s18, scenarios$s19, scenarios$s20, scenarios$s20b),
              function(f) f())
  i <- function(x, t) which.min(abs(x$time - t))
  do.call(rbind, lapply(d, function(x) data.frame(
    scenario = x$scenario[1],
    GCase_y2 = round(x$GCASE[i(x, 730)], 3),
    GCase_y3 = round(x$GCASE[i(x, 1095)], 3),
    GCase_y5 = round(x$GCASE[i(x, 1825)], 3),
    FIBN_y5  = round(x$FIBN[i(x, 1825)], 3),
    MMSE_y3  = round(x$MMSE[i(x, 1095)], 2),
    MMSE_y5  = round(x$MMSE[i(x, 1825)], 2),
    row.names = NULL)))
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("dlb.run.scenarios"))) {
  all <- run_all()
  print(summarise_scn(all))
  print(median_survival(all))
  print(chei_ratio())
}
