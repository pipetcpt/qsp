## ============================================================================
##  Progressive Supranuclear Palsy (PSP) -- QSP model for mrgsolve
##  ---------------------------------------------------------------------------
##  psp_mrgsolve_model.R
##
##  71 ODEs | 8 anatomical regions | 16 agents | 29 scenarios
##  Time unit: DAYS.  t = 0 is the biological onset of the tauopathy, NOT the
##  clinic visit.  Trials enrol at t = TENROL (default 1278 d ~ 3.5 y), which is
##  where the model reproduces the observed enrolment PSPRS of ~38.
##
##  ---------------------------------------------------------------------------
##  WHY THIS MODEL IS BUILT THE WAY IT IS
##  ---------------------------------------------------------------------------
##  Two PSP programmes reported target-engagement numbers that point the
##  OPPOSITE way from their clinical results:
##
##    gosuranemab (PASSPORT, n=486)   CSF unbound N-terminal tau  -98%
##                                    PSPRS wk52  10.4 vs 10.6  (p = 0.85)
##                                    CSF NfL     unchanged
##    NIO752 (MAPT ASO, phase 1, n=59) CSF total tau / p-tau181   -20%
##                                    CSF NfL     stabilised (placebo +40%)
##
##  A model with ONE "tau lowering" axis cannot contain both facts.  This model
##  therefore separates tau into pools that differ by 10^3-10^4 in size:
##
##      intraneuronal soluble monomer   TMON   ~2000 nM   (the engine)
##      regional aggregate load         A1..A8 0-1        (the pathology)
##      ISF N-terminal fragments        ETN    ~1 nM      (the assay analyte)
##      ISF seed-competent assemblies   ETS    ~0.02 nM   (the causal species)
##      CSF                             CTN, CTS ~10 pM
##
##  SEVEN STRUCTURAL CLAIMS, each of which is a line of code, not a fitted knob
##  ---------------------------------------------------------------------------
##  1. POOL SIZE.  F_EXT = ETN/TMON ~ 5e-4.  Any extracellular-only mechanism
##     is bounded by F_EXT x PHI_ACC at ANY dose or affinity.
##
##  2. ENGAGEMENT PARADOX.  The -98% is not fitted.  2000 mg q4w with V1 = 3 L,
##     t1/2 = 25 d gives Cp,avg ~ 300 mg/L = 2000 nM; Kp,brain = 0.15% gives
##     AISF ~ 3 nM; a 0.06 nM Kd gives free fraction 1/(1+50) = 0.02, i.e. -98%.
##     Real PK x real affinity => the observed number.  98% ENGAGEMENT IS ITSELF
##     THE EVIDENCE THAT THE POOL WAS SMALL ENOUGH TO SATURATE.
##
##  3. EPITOPE.  ETN (N-terminal fragments; abundant; NOT seeding-competent) and
##     ETS (mid-domain/MTBR assemblies; scarce; N-terminus already cleaved off by
##     caspase-3 / calpain / AEP) are two states with two Kd's.  An N-terminal
##     mAb has KD_NTAB = 0.06 nM and KD_STAB = 1e6 nM.  One epitope parameter
##     produces -98% of the measured species and 0% of the causal one.
##
##  4. GEOMETRY.  Progression is sequential recruitment of 8 nuclei through a
##     connectome (matrix W).  Each region's aggregation is LOGISTIC and
##     saturates in ~1.5-2 y.  A sum of saturating sigmoids is linear only if
##     regions ignite at a constant rate -- so the famous ~10-11 PSPRS pts/yr
##     linearity is EVIDENCE OF A FRONT, and it means every already-involved
##     nucleus sits at plateau where d/dk of the plateau is ~0.  This is why
##     tideglusib, davunetide and the OGA class are structurally inert: they
##     move k, and k multiplies a saturated quantity.
##
##  5. ELASTICITY, not a square-root law.  The continuum Fisher-KPP result
##     v ~ sqrt(D|k|M) is the LOWER bound on the sensitivity of the clinical
##     slope to monomer knockdown; a discrete pulled front on a sparse graph
##     gives v ~ k|dx / ln(k/D), i.e. nearly LINEAR, as the upper bound.  The
##     model does not assert either -- it MEASURES the elasticity
##     e = dln(slope)/dln(M) and reports it, because that one number decides
##     whether the Preserve ASO dose is sufficient.  (fn psp_elasticity())
##
##  6. ASO GRADIENT.  Intrathecal ASO exposure is cortex >> deep nuclei, but the
##     wave front sits in deep nuclei while the CSF biomarker is dominated by
##     cortex.  Parameter ASO_DEEP (0.4) makes the CSF readout OVERSTATE the
##     effect on the front.  The model therefore predicts a Preserve effect
##     SMALLER than the -20% CSF signal suggests, and says by how much.
##
##  7. TIMING.  Aggregate load does not cause the score; it causes a loss rate,
##     which integrates into neuron count, which crosses a regional reserve
##     threshold, which moves the score through two transit compartments.  Two
##     integrations and one threshold sit between drug and endpoint.  Blocking
##     100% of aggregation for 52 weeks from enrolment is worth a couple of
##     points; the same block 3 y earlier is worth most of the disease.
##     (scenarios S24 / S25 make this an in-model counterfactual)
##
##  ---------------------------------------------------------------------------
##  CALIBRATION ANCHORS (see psp_validate() at the bottom for the full table)
##  ---------------------------------------------------------------------------
##    PSPRS at enrolment (t = 3.5 y)                        ~38      pts
##    PSPRS 52-week change, placebo                         10.6     pts
##    Natural-history slope (Golbe 2007)                    11.3     pts/y
##    Gosuranemab CSF unbound N-terminal tau                -98      %
##    Gosuranemab PSPRS 52-wk change                        10.4     pts
##    NIO752 CSF total tau / p-tau181                       -20      %
##    NIO752 CSF NfL                                        ~0 vs +40 %
##    Median survival from onset, PSP-RS                     7.9     y
##    Median survival from onset, PSP-P                      ~10     y
##    PSP-P PSPRS slope                                      ~5      pts/y
##    Vertical saccade peak velocity at enrolment          150-200   deg/s
##    Midbrain area at enrolment                           70-85     mm^2
##    MRPI at enrolment                                    >13.55
##
##  ---------------------------------------------------------------------------
##  REQUIREMENTS:  R >= 4.1, mrgsolve >= 1.0 (developed on mrgsolve 2.0.1)
##  USAGE:         source("psp_mrgsolve_model.R"); psp_validate()
## ============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

# =============================================================================
#  MODEL CODE
# =============================================================================

psp_code <- r"---(

$PROB
# Progressive Supranuclear Palsy -- 4R tauopathy travelling-wave QSP model
- 71 ODEs, 8 regions, 16 agents
- time in DAYS from biological onset

$PARAM @annotated
// ---------------------------------------------------------------------------
// GENETICS / CONSTITUTION
// ---------------------------------------------------------------------------
GMAPT   :  2    : MAPT H1 haplotype dose 0/1/2 (raises MAPT expression)
FH1     :  0.10 : fractional MAPT mRNA increase per H1 allele
TRIM11F :  1.0  : TRIM11 disaggregase capacity multiplier (1 = normal)

// ---------------------------------------------------------------------------
// MAPT TRANSCRIPTION AND MONOMER (the ONLY node upstream of the wave)
// ---------------------------------------------------------------------------
KDEG_M  :  0.50 : MAPT mRNA turnover (1/day)
KSYN_T  :  100  : tau monomer synthesis per unit mRNA (nM/day)
// KSYN_T0 defines an UNMOVING monomer reference: the growth term must be
// compared against a fixed baseline, or a change in expression is invisible.
KSYN_T0 :  100  : nominal tau synthesis rate defining the monomer reference
MAPT_KD :  0    : imposed fractional MAPT knockdown (elasticity probe, 0-1)
TKD     :  0    : time at which the imposed knockdown switches on (days)
KDEG_T  :  0.05 : tau monomer degradation (1/day, t1/2 ~14 d)
KREL_T  :  5e-4 : fractional monomer release to ISF (1/day)
TMON0   :  2000 : reference intraneuronal soluble monomer (nM)

// ---------------------------------------------------------------------------
// POST-TRANSLATIONAL MODIFICATION
// ---------------------------------------------------------------------------
KDEG_G  :  0.20 : GSK-3beta turnover (1/day)
WCYTO_G :  0.60 : cytokine drive on GSK-3beta tone
KOGT    :  0.50 : O-GlcNAc addition (1/day)
KOGA    :  3.00 : O-GlcNAcase removal (1/day)
KP300   :  0.20 : p300/CBP tau acetylation (1/day)
KDEAC   :  1.80 : SIRT1/HDAC6 deacetylation (1/day)
KMTS    :  0.30 : microtubule-binding equilibration rate (1/day)
KD_MT   :  1.00 : phospho-tone at which MT-bound fraction halves
WP_AGG  :  1.20 : weight of hyperphosphorylation on aggregation propensity
WA_AGG  :  1.50 : weight of acetylation on aggregation propensity
WO_AGG  :  2.00 : weight of O-GlcNAc (anti-aggregant)
WM_AGG  :  0.80 : weight of MT sequestration (protective)

// ---------------------------------------------------------------------------
// PROTEOSTASIS
// ---------------------------------------------------------------------------
KPG     :  0.10 : progranulin equilibration (1/day)
KAU     :  0.08 : autophagy-lysosome capacity equilibration (1/day)
WAGG_LYS:  0.35 : lysosomal impairment per unit mean aggregate load
AUTO_MIN:  0.25 : floor on clearance capacity

// ---------------------------------------------------------------------------
// AGGREGATION AND THE WAVE
// ---------------------------------------------------------------------------
KAGG    :  0.0055635 : intrinsic templated-growth rate constant (1/day)
KDIS    :  6e-4   : clearance-driven disaggregation (1/day)
FDIS    :  1.0    : fold-increase in disaggregation after TDIS (counterfactual)
TDIS    :  1e9    : time at which FDIS switches on (days)
KTR     :  0.0031841 : inter-regional transfer coefficient (effective seed units)
// (the remainder travels by exosome / tunnelling nanotube / synaptic microdomain)
PHI_ACC :  0.45   : fraction of transfer flux exposed to bulk ISF
ASEED   :  0.0109512 : aggregate load in the origin nucleus at t = 0
ORIG_RS :  1      : 1 = PSP-RS origin (STN + midbrain), 0 = PSP-P origin (SNc)
FBLOCK  :  0      : in-model counterfactual: fraction of growth blocked
TBLOCK  :  1e9    : time at which FBLOCK switches on (days)
ASO_DEEP:  0.40   : fraction of ASO knockdown realised at the deep-nucleus front

// regional templating propensity (4R vulnerability)
KP1 : 1.20 : STN
KP2 : 1.10 : GPi/SNr
KP3 : 0.95 : SNc
KP4 : 1.30 : midbrain tegmentum (riMLF / rip)
KP5 : 1.20 : pedunculopontine nucleus
KP6 : 0.70 : dentate nucleus
KP7 : 0.5997 : frontal cortex
KP8 : 0.80 : locus coeruleus / bulbar

// connectome weights W[receiver][donor]
W21 : 0.90 : STN  -> GPi/SNr
W31 : 0.50 : STN  -> SNc
W12 : 0.30 : GPi  -> STN
W42 : 0.80 : GPi  -> midbrain tegmentum
W23 : 0.30 : SNc  -> GPi
W54 : 0.80 : midbrain -> PPN
W45 : 0.30 : PPN -> midbrain
W65 : 0.40 : PPN -> dentate
W74 : 0.3505 : midbrain -> frontal
W76 : 0.30 : dentate -> frontal (via thalamus)
W85 : 0.50 : PPN -> locus coeruleus / bulbar

// ---------------------------------------------------------------------------
// NEURONAL LOSS
// ---------------------------------------------------------------------------
KDEG_N  :  3.3371e-4 : aggregate-driven neuronal loss rate (1/day)
KAGE_N  :  1.0e-5 : age-related loss (1/day)
WCY_N   :  0.50   : cytokine contribution to neuronal toxicity
WIF_N   :  0.35   : type-I interferon contribution
WLP_N   :  0.30   : lipid-peroxidation contribution
KT1 : 1.00 : regional toxicity multiplier STN
KT2 : 1.00 : GPi/SNr
KT3 : 0.95 : SNc
KT4 : 1.5465 : midbrain tegmentum
KT5 : 1.05 : PPN
KT6 : 0.80 : dentate
KT7 : 0.75 : frontal cortex
KT8 : 0.90 : LC / bulbar

// ---------------------------------------------------------------------------
// EXTRACELLULAR TAU
// ---------------------------------------------------------------------------
KRELN   :  0.45  : ISF N-terminal fragment appearance (nM/day per unit source)
WAGG_N  :  1.40  : aggregate contribution to N-fragment release
KOUT_N  :  0.60  : ISF free N-fragment clearance (1/day)
KRELS   :  0.011 : ISF seed-competent appearance (nM/day per unit source)
WAGG_S  :  3.00  : aggregate contribution to seed release
KOUT_S  :  0.60  : ISF free seed clearance (1/day)
// An IgG:tau complex leaves the ISF SLOWER than free fragmented tau, not
// faster: it is larger and its exit is bulk-flow / FcRn limited.  This is the
// standard target-stabilisation sink, and it is why an anti-tau mAb drives
// UNBOUND tau down while TOTAL tau goes UP -- two readouts, opposite signs,
// one mechanism.
KOUT_CX :  0.15  : clearance of antibody:tau complex (1/day, 4x SLOWER than free)
KDEATH_T:  6.0   : release from dying neurons (nM per unit death rate)
KISF_CSF:  0.35  : ISF -> CSF transfer (1/day)
RVOL    :  0.030 : ISF:CSF volume ratio factor
KCSF_OUT:  4.0   : CSF bulk turnover (1/day, ~500 mL/d over 150 mL)

// ---------------------------------------------------------------------------
// GLIA / INFLAMMATION
// ---------------------------------------------------------------------------
MG0     :  0.05 : baseline activated-microglia fraction
KMG     :  0.05 : microglial equilibration (1/day)
WMG_A   :  0.55 : aggregate drive on microglial activation
WMG_D   :  120  : neuronal-death drive on microglial activation
WMG_IF  :  0.30 : interferon drive on microglial activation
CYTO0   :  0.10 : baseline cytokine tone
KCY     :  0.15 : cytokine equilibration (1/day)
WCY_MG  :  1.20 : microglia -> cytokine gain
AS0     :  0.08 : baseline reactive-astrocyte fraction
KAS     :  0.05 : astrocyte equilibration (1/day)
WAS_C   :  0.60 : cytokine -> astrocyte gain
WAS_A   :  0.45 : aggregate -> astrocyte gain (tufted astrocytes)
IFN0    :  0.05 : baseline type-I interferon tone
KIF     :  0.08 : interferon equilibration (1/day)
WIF_A   :  0.40 : aggregate -> LINE-1 -> interferon gain
LPO0    :  0.10 : baseline lipid peroxidation
KLP     :  0.10 : lipid-peroxidation equilibration (1/day)
WLP_A   :  0.35 : aggregate -> lipid peroxidation gain

// ---------------------------------------------------------------------------
// FLUID AND IMAGING BIOMARKERS
// ---------------------------------------------------------------------------
NFL_IN  :  199.8 : baseline CSF NfL input (pg/mL/day)
// NfL comes overwhelmingly from large myelinated axons, so the source is
// volume-weighted, NOT the plain mean death rate.  The late-igniting frontal
// white matter dwarfs the STN, which is why CSF NfL rises steeply exactly
// through the symptomatic phase while the mean death rate is already flat.
VN1 : 0.02 : NfL axonal-volume weight, STN
VN2 : 0.07 : GPi/SNr
VN3 : 0.04 : SNc
VN4 : 0.12 : midbrain tegmentum
VN5 : 0.08 : PPN
VN6 : 0.12 : dentate / superior cerebellar peduncle
VN7 : 0.45 : frontal cortex and its white matter
VN8 : 0.10 : LC / bulbar
KNFL_S  :  6.000e5 : CSF NfL released per unit weighted neuronal death rate
KNFL_A  :  718.6 : CSF NfL from ongoing (non-lethal) axonal injury, volume-weighted
KNFL_O  :  0.30 : CSF NfL elimination (1/day)
KP_NFL  :  4.58e-4 : CSF -> plasma NfL transfer (1/day)
KNFLP_O :  0.020 : plasma NfL elimination (1/day, t1/2 ~35 d)
GFAP_IN :  3.6  : baseline plasma GFAP input (pg/mL/day)
KGF     :  260  : astrocyte -> plasma GFAP gain
KGF_O   :  0.10 : plasma GFAP elimination (1/day)
MIDA_MX :  125  : healthy midbrain area (mm2)
KMIDA   :  0.010 : midbrain area equilibration (1/day)
MIDA_FL :  0.50 : residual midbrain area fraction at complete R4 loss
SCP_MX  :  3.40 : healthy superior cerebellar peduncle width (mm)
KSCP    :  0.008 : SCP equilibration (1/day)
SCP_FL  :  0.68 : residual SCP fraction at complete R6 loss
PONS_A  :  480  : pons area (mm2, assumed spared)
MCP_W   :  8.00 : middle cerebellar peduncle width (mm, assumed spared)

// ---------------------------------------------------------------------------
// CLINICAL SCORE -- reserve thresholds, Hill slopes, PSPRS item weights
// ---------------------------------------------------------------------------
TH1 : 0.45 : reserve threshold STN
TH2 : 0.40 : GPi/SNr
TH3 : 0.50 : SNc
TH4 : 0.15 : midbrain tegmentum (LOWEST reserve -> earliest sign)
TH5 : 0.25 : PPN
TH6 : 0.50 : dentate
TH7 : 0.55 : frontal cortex
TH8 : 0.45 : LC / bulbar
HL1 : 4.0 : Hill slope STN
HL2 : 4.0 : GPi/SNr
HL3 : 4.0 : SNc
HL4 : 2.0 : midbrain tegmentum (shallow -> gradual, unthresholded-like)
HL5 : 3.0 : PPN
HL6 : 4.0 : dentate
HL7 : 3.0 : frontal cortex
HL8 : 3.0 : LC / bulbar
WS1 :  8 : PSPRS weight STN
WS2 : 12 : GPi/SNr
WS3 : 10 : SNc (limb motor)
WS4 : 14 : midbrain tegmentum (oculomotor)
WS5 : 22 : PPN (gait / midline -- dominant contributor)
WS6 :  6 : dentate
WS7 : 14 : frontal cortex (mentation + history)
WS8 : 14 : LC / bulbar
TSC : 1.0013 : global scale on all regional reserve thresholds
HSC : 1.6312 : global scale on all reserve Hill slopes
KLAG    : 0.0055 : committed-damage transit rate (1/day, 2 compartments)
PSPRS_MX: 100 : PSPRS ceiling
VSV0    : 450 : healthy vertical saccade peak velocity (deg/s)
VSV_FL  : 40  : floor on vertical saccade peak velocity (deg/s)
GAM_SAC : 1.0 : exponent linking riMLF survival to saccade velocity
KDYS    : 0.010 : dysphagia equilibration (1/day)

// ---------------------------------------------------------------------------
// SURVIVAL (competing risks)
// ---------------------------------------------------------------------------
HP0     : 3.1619e-5 : baseline aspiration-pneumonia hazard (1/day)
BDYS    : 0.95   : log-hazard per unit dysphagia score
HO0     : 2.9369e-5 : baseline other-cause hazard (1/day)
BPSP    : 1.10   : log-hazard per PSPRS/50
PEG     : 0      : 1 = PEG gastrostomy in place
PEG_RR  : 0.45   : relative aspiration hazard with PEG

// ---------------------------------------------------------------------------
// ANTI-TAU mAb PK AND EPITOPE (gosuranemab / tilavonemab / mid-domain probe)
// ---------------------------------------------------------------------------
V1_AB   : 3.0   : mAb central volume (L)
V2_AB   : 3.6   : mAb peripheral volume (L)
CL_AB   : 0.085 : mAb clearance (L/day, t1/2 ~ 25 d)
Q_AB    : 0.45  : mAb intercompartmental clearance (L/day)
KP_AB   : 0.0015: brain:plasma partition for IgG (0.15%)
KBBB    : 0.50  : ISF antibody equilibration (1/day)
KCSFE   : 0.60  : ISF <-> CSF antibody equilibration (1/day)
MW_AB   : 150   : mAb molecular weight (kDa)
KD_NTAB : 0.0344 : Kd for ISF N-terminal fragments (nM)  <- N-terminal mAb
KD_STAB : 1e6   : Kd for ISF seed-competent species (nM) <- ~no binding
KVAX_D  : 0.0116: AADvac1 titre decay (1/day, t1/2 ~60 d)
KVAX_AB : 0.020 : titre unit -> ISF anti-tau IgG (nM per unit)

// ---------------------------------------------------------------------------
// MAPT ASO PK (NIO752)
// ---------------------------------------------------------------------------
KSC_ASO : 0.35  : CSF -> brain tissue distribution (1/day)
KCL_ASO : 1.20  : CSF ASO elimination (1/day)
KTIS_ASO: 0.0092: brain tissue ASO elimination (1/day, t1/2 ~75 d)
EMAX_ASO: 0.85  : maximal fractional MAPT mRNA knockdown
IC50_ASO: 1.15  : tissue ASO amount at half-maximal knockdown (mg-equivalents)

// ---------------------------------------------------------------------------
// SMALL MOLECULES AND SYMPTOMATIC AGENTS
// ---------------------------------------------------------------------------
KA_TID : 12.0 : tideglusib absorption (1/day)
F_TID  : 0.50 : tideglusib bioavailability
V_TID  : 200  : tideglusib effect-site volume (L)
KE_TID : 1.00 : tideglusib effect-site elimination (1/day, irreversible-target surrogate)
EMX_TID: 0.75 : maximal fractional GSK-3beta inactivation
IC5_TID: 1.20 : tideglusib EC50 (mg/L)

KA_OGA : 12.0 : OGA inhibitor absorption (1/day)
F_OGA  : 0.80 : OGA inhibitor bioavailability
V_OGA  : 100  : OGA inhibitor volume (L)
KE_OGA : 0.70 : OGA inhibitor elimination (1/day)
EMX_OGA: 0.90 : maximal fractional OGA inhibition
IC5_OGA: 0.10 : OGA inhibitor EC50 (mg/L)

KA_DAV : 6.0  : davunetide nasal-to-CNS transfer (1/day)
F_DAV  : 0.030: davunetide CNS bioavailability
V_DAV  : 50   : davunetide volume (L)
KE_DAV : 4.0  : davunetide elimination (1/day)
EMX_DAV: 0.35 : maximal increase in MT-bound tau fraction
IC5_DAV: 0.004: davunetide EC50 (mg/L)

KA_SAL : 6.0  : salsalate absorption (1/day)
F_SAL  : 0.90 : salsalate bioavailability
V_SAL  : 12   : salicylate volume (L)
KE_SAL : 0.90 : salicylate elimination (1/day)
EMX_SAL: 0.55 : maximal fractional p300 inhibition
IC5_SAL: 140  : salicylate EC50 (mg/L)

V_FAS  : 60   : fasudil volume (L)
KE_FAS : 6.0  : fasudil elimination (1/day)
EMX_FAS: 0.30 : maximal fractional increase in autophagic capacity
IC5_FAS: 0.15 : fasudil EC50 (mg/L)

KA_TPN : 12.0 : TPN-101 absorption (1/day)
F_TPN  : 0.60 : TPN-101 bioavailability
V_TPN  : 150  : TPN-101 volume (L)
KE_TPN : 1.40 : TPN-101 elimination (1/day)
EMX_TPN: 0.70 : maximal fractional LINE-1 RT / interferon suppression
IC5_TPN: 0.50 : TPN-101 EC50 (mg/L)

KA_LM  : 12.0 : LM11A-31 absorption (1/day)
F_LM   : 0.40 : LM11A-31 bioavailability
V_LM   : 250  : LM11A-31 volume (L)
KE_LM  : 2.80 : LM11A-31 elimination (1/day)
EMX_LM : 0.75 : maximal fractional suppression of tau acetylation
IC5_LM : 0.60 : LM11A-31 EC50 (mg/L)

KA_EZE : 12.0 : ezeprogind absorption (1/day)
F_EZE  : 0.55 : ezeprogind bioavailability
V_EZE  : 100  : ezeprogind volume (L)
KE_EZE : 2.00 : ezeprogind elimination (1/day)
EMX_EZE: 0.60 : maximal fractional progranulin increase
IC5_EZE: 0.25 : ezeprogind EC50 (mg/L)

KA_RIL : 12.0 : riluzole absorption (1/day)
F_RIL  : 0.60 : riluzole bioavailability
V_RIL  : 250  : riluzole volume (L)
KE_RIL : 1.40 : riluzole elimination (1/day)
EMX_RIL: 0.12 : maximal fractional reduction in excitotoxic neuronal loss
IC5_RIL: 0.20 : riluzole EC50 (mg/L)

KA_LDP : 24.0 : levodopa absorption (1/day)
F_LDP  : 0.80 : levodopa bioavailability
V_LDP  : 50   : levodopa volume (L)
KE_LDP : 13.0 : levodopa elimination (1/day)
EMX_LDP: 0.22 : maximal symptomatic relief of the nigral score component
IC5_LDP: 1.50 : levodopa EC50 (mg/L)

V_ZOL  : 40   : zolpidem volume (L)
KE_ZOL : 6.90 : zolpidem elimination (1/day)
EMX_ZOL: 0.10 : maximal symptomatic relief of the pallidal score component
IC5_ZOL: 0.08 : zolpidem EC50 (mg/L)

$CMT @annotated
// regional aggregate load (0-1 of regional maximum)
A1 : STN aggregate load
A2 : GPi/SNr aggregate load
A3 : SNc aggregate load
A4 : midbrain tegmentum aggregate load
A5 : PPN aggregate load
A6 : dentate aggregate load
A7 : frontal cortex aggregate load
A8 : LC/bulbar aggregate load
// regional surviving-neuron fraction
N1 : STN surviving fraction
N2 : GPi/SNr surviving fraction
N3 : SNc surviving fraction
N4 : midbrain tegmentum surviving fraction
N5 : PPN surviving fraction
N6 : dentate surviving fraction
N7 : frontal cortex surviving fraction
N8 : LC/bulbar surviving fraction
// molecular
MRNA  : MAPT mRNA (relative)
TMON  : intraneuronal soluble tau monomer (nM)
GSK   : GSK-3beta activity (relative)
OGLC  : O-GlcNAc-tau fraction
ACT   : acetyl-tau fraction
MTS   : microtubule-bound tau fraction
PGRNS : progranulin (relative)
AUTO  : autophagy-lysosome clearance capacity (relative)
// extracellular
ETN : ISF N-terminal tau fragments (nM)
ETS : ISF seed-competent tau (nM)
CTN : CSF N-terminal tau (nM)
CTS : CSF seed-competent tau (nM)
// glia
MGA  : activated microglia fraction
CYTOK: cytokine tone (relative)
ASTR : reactive astrocyte fraction
IFNA : type-I interferon tone (relative)
LPOX : lipid peroxidation (relative)
// biomarkers
NFLC  : CSF NfL (pg/mL)
NFLP  : plasma NfL (pg/mL)
GFAPP : plasma GFAP (pg/mL)
MIDA  : midbrain area (mm2)
SCPW  : superior cerebellar peduncle width (mm)
// clinical
LAG1 : committed-damage transit 1 (PSPRS points)
LAG2 : committed-damage transit 2 (PSPRS points)
DYS  : dysphagia severity (0-4)
CIFP : cumulative incidence of pneumonia death
CIFO : cumulative incidence of other-cause death
SURV : overall survival probability
// mAb PK
MABC : mAb central (mg)
MABP : mAb peripheral (mg)
AISF : ISF antibody (nM)
ACSF : CSF antibody (nM)
// ASO PK
ASOC : CSF ASO (mg)
ASOT : brain tissue ASO (mg)
// small molecules
TIDG : tideglusib gut (mg)
TIDC : tideglusib effect site (mg)
OGAG : OGA inhibitor gut (mg)
OGAC : OGA inhibitor central (mg)
DAVD : davunetide nasal depot (mg)
DAVC : davunetide CNS (mg)
SALG : salsalate gut (mg)
SALC : salicylate central (mg)
FASC : fasudil central (mg)
TPNG : TPN-101 gut (mg)
TPNC : TPN-101 central (mg)
LMG  : LM11A-31 gut (mg)
LMC  : LM11A-31 central (mg)
EZEG : ezeprogind gut (mg)
EZEC : ezeprogind central (mg)
RILG : riluzole gut (mg)
RILC : riluzole central (mg)
LDPG : levodopa gut (mg)
LDPC : levodopa central (mg)
ZOLC : zolpidem central (mg)
VTIT : AADvac1 anti-tau titre (arbitrary units)

$GLOBAL
#define MRNAS  (1.0 + FH1*GMAPT)
// Baseline acetyl-tau fraction and the clearance factor it implies.  These must
// appear in BOTH the initial condition and the monomer reference, or the model
// starts off its own steady state and the growth term sees a phantom 2-5%
// "knockdown" of the wrong sign.
#define ACT0   (KP300/(KP300 + KDEAC))
#define CLR0   (1.0 - 0.50*ACT0)
#define TMONS  (KSYN_T*MRNAS/(KDEG_T*CLR0 + KREL_T))
#define TMREF  (KSYN_T0*MRNAS/(KDEG_T*CLR0 + KREL_T))

// reserve-gated regional score contribution
double resv(double loss, double thr, double hl) {
  double L = loss; if (L < 1e-12) L = 1e-12;
  double a = pow(L, hl);
  double b = pow(thr, hl);
  return a/(a+b);
}
// bounded Emax
double emx(double e, double c, double ic50) {
  if (c <= 0.0) return 0.0;
  return e*c/(ic50 + c);
}

$MAIN
// ---- initial conditions -----------------------------------------------------
// PSP-RS: the front is seeded in STN and the midbrain tegmentum.
// PSP-P : the front is seeded in SNc, whose connectome outflow is weaker,
//         which is the ONLY difference between the two phenotypes here.
A1_0 = (ORIG_RS > 0.5) ? ASEED       : 0.0;
A4_0 = (ORIG_RS > 0.5) ? 0.30*ASEED  : 0.0;
A3_0 = (ORIG_RS > 0.5) ? 0.0         : ASEED;
A2_0 = 0.0; A5_0 = 0.0; A6_0 = 0.0; A7_0 = 0.0; A8_0 = 0.0;

N1_0 = 1.0; N2_0 = 1.0; N3_0 = 1.0; N4_0 = 1.0;
N5_0 = 1.0; N6_0 = 1.0; N7_0 = 1.0; N8_0 = 1.0;

MRNA_0  = MRNAS;
TMON_0  = TMONS;
GSK_0   = 1.0;
OGLC_0  = KOGT/(KOGT + KOGA);
ACT_0   = ACT0;
MTS_0   = 1.0/(1.0 + 1.0/KD_MT);
PGRNS_0 = 1.0;
AUTO_0  = 1.0;

// analytic steady states of the two extracellular pools at zero pathology
ETN_0 = KRELN/KOUT_N;
ETS_0 = KRELS/KOUT_S;
CTN_0 = KISF_CSF*RVOL*ETN_0/KCSF_OUT;
CTS_0 = KISF_CSF*RVOL*ETS_0/KCSF_OUT;

MGA_0   = MG0;
CYTOK_0 = CYTO0;
ASTR_0  = AS0;
IFNA_0  = IFN0;
LPOX_0  = LPO0;

NFLC_0  = NFL_IN/KNFL_O;
NFLP_0  = KP_NFL*NFLC_0/KNFLP_O;
GFAPP_0 = (GFAP_IN + KGF*AS0)/KGF_O;
MIDA_0  = MIDA_MX;
SCPW_0  = SCP_MX;

LAG1_0 = 0.0; LAG2_0 = 0.0; DYS_0 = 0.0;
CIFP_0 = 0.0; CIFO_0 = 0.0; SURV_0 = 1.0;

$ODE
// ===========================================================================
//  DERIVED QUANTITIES
// ===========================================================================
double AUT = fmax(AUTO, AUTO_MIN);

// -- drug effects ------------------------------------------------------------
double CP_AB   = MABC/V1_AB;                        // mg/L
double ITIDE   = emx(EMX_TID, TIDC/V_TID,  IC5_TID);
double IOGA    = emx(EMX_OGA, OGAC/V_OGA,  IC5_OGA);
double EDAVU   = emx(EMX_DAV, DAVC/V_DAV,  IC5_DAV);
double ISALS   = emx(EMX_SAL, SALC/V_SAL,  IC5_SAL);
double EFASU   = emx(EMX_FAS, FASC/V_FAS,  IC5_FAS);
double ITPN    = emx(EMX_TPN, TPNC/V_TPN,  IC5_TPN);
double ILM11   = emx(EMX_LM,  LMC/V_LM,    IC5_LM );
double EEZE    = emx(EMX_EZE, EZEC/V_EZE,  IC5_EZE);
double IRILU   = emx(EMX_RIL, RILC/V_RIL,  IC5_RIL);
double IASO    = emx(EMAX_ASO, ASOT,       IC50_ASO);
double MKD     = (SOLVERTIME >= TKD) ? MAPT_KD : 0.0;

// -- rapid-equilibrium antibody binding (the whole of AXIS 2 and AXIS 3) -----
// The antibody is in vast molar excess over both ISF tau pools, which is why
// rapid equilibrium is the right approximation AND why 98% engagement of the
// N-terminal pool is achievable at 0.15% brain partition.
double FRN = 1.0/(1.0 + AISF/KD_NTAB);   // free fraction of ETN
double FRS = 1.0/(1.0 + AISF/KD_STAB);   // free fraction of ETS
double FRNC= 1.0/(1.0 + ACSF/KD_NTAB);   // free fraction in CSF

// -- PTM state -> aggregation propensity ------------------------------------
double PHOS  = GSK*(1.0 - 0.60*(OGLC/(KOGT/(KOGT+KOGA))-1.0)*0.5);
if (PHOS < 0.05) PHOS = 0.05;
double MTS_T = fmin(0.95, (1.0/(1.0 + PHOS/KD_MT))*(1.0 + EDAVU));

double PHI0  = (1.0 + WP_AGG*1.0)*(1.0 + WA_AGG*ACT0)
             / ((1.0 + WO_AGG*(KOGT/(KOGT+KOGA)))*(1.0 + WM_AGG*(1.0/(1.0+1.0/KD_MT))));
double PHI   = (1.0 + WP_AGG*PHOS)*(1.0 + WA_AGG*ACT)
             / ((1.0 + WO_AGG*OGLC)*(1.0 + WM_AGG*MTS));
double PHI_N = PHI/PHI0;                 // 1 at healthy baseline

// -- monomer seen BY THE FRONT (AXIS 6: ASO reaches cortex, not deep nuclei) -
//  TMON/TMREF mixes two different things: (i) the DRUG-attributable reduction
//  in MAPT expression, which is subject to the intrathecal depth gradient, and
//  (ii) a PHYSIOLOGICAL rise driven by aggregate-impaired autophagy, which is
//  not.  Scaling the product by ASO_DEEP would silently scale the autophagy
//  feedback too, so the two are separated explicitly.
double RASO  = (1.0 - IASO)*(1.0 - MKD);               // drug-attributable
if (RASO < 1e-6) RASO = 1e-6;
double RPHYS = (TMON/TMREF)/RASO;                      // everything else
double MFRONT= RPHYS*(1.0 - ASO_DEEP*(1.0 - RASO));    // what the front sees
if (MFRONT < 0.05) MFRONT = 0.05;

// -- counterfactual growth block (in-model, not a drug) ---------------------
double BLK = (SOLVERTIME >= TBLOCK) ? (1.0 - FBLOCK) : 1.0;

double KG = KAGG*PHI_N*MFRONT*BLK;
double KDX = KDIS*((SOLVERTIME >= TDIS) ? FDIS : 1.0);

// -- accessible fraction of the inter-regional transfer flux ----------------
// An extracellular antibody can only intercept the part of the transfer flux
// that is exposed to bulk ISF.  The rest travels by exosome, tunnelling
// nanotube or a synaptic microdomain a 150 kDa IgG at 1-10 nM cannot police.
double ACCESS = (1.0 - PHI_ACC) + PHI_ACC*FRS;

// -- sources of aggregate for the connectome (living neurons only) ----------
double S1 = A1*N1;
double S2 = A2*N2;
double S3 = A3*N3;
double S4 = A4*N4;
double S5 = A5*N5;
double S6 = A6*N6;
double S7 = A7*N7;
double S8 = A8*N8;

double SD1 = KTR*ACCESS*( W12*S2 );
double SD2 = KTR*ACCESS*( W21*S1 + W23*S3 );
double SD3 = KTR*ACCESS*( W31*S1 );
double SD4 = KTR*ACCESS*( W42*S2 + W45*S5 );
double SD5 = KTR*ACCESS*( W54*S4 );
double SD6 = KTR*ACCESS*( W65*S5 );
double SD7 = KTR*ACCESS*( W74*S4 + W76*S6 );
double SD8 = KTR*ACCESS*( W85*S5 );

// ===========================================================================
//  1. REGIONAL AGGREGATION -- logistic (saturating) growth + connectome seeding
//     The saturation is the reason local-kinetic drugs are inert (AXIS 4).
// ===========================================================================
dxdt_A1 = KG*KP1*(A1 + SD1)*(1.0 - A1) - KDX*AUT*TRIM11F*A1;
dxdt_A2 = KG*KP2*(A2 + SD2)*(1.0 - A2) - KDX*AUT*TRIM11F*A2;
dxdt_A3 = KG*KP3*(A3 + SD3)*(1.0 - A3) - KDX*AUT*TRIM11F*A3;
dxdt_A4 = KG*KP4*(A4 + SD4)*(1.0 - A4) - KDX*AUT*TRIM11F*A4;
dxdt_A5 = KG*KP5*(A5 + SD5)*(1.0 - A5) - KDX*AUT*TRIM11F*A5;
dxdt_A6 = KG*KP6*(A6 + SD6)*(1.0 - A6) - KDX*AUT*TRIM11F*A6;
dxdt_A7 = KG*KP7*(A7 + SD7)*(1.0 - A7) - KDX*AUT*TRIM11F*A7;
dxdt_A8 = KG*KP8*(A8 + SD8)*(1.0 - A8) - KDX*AUT*TRIM11F*A8;

double ABAR = (A1+A2+A3+A4+A5+A6+A7+A8)/8.0;

// ===========================================================================
//  2. NEURONAL LOSS -- the first integration between drug and endpoint
// ===========================================================================
double GTOX = WCY_N*(CYTOK-CYTO0) + WIF_N*(IFNA-IFN0) + WLP_N*(LPOX-LPO0);
double RN   = KDEG_N*(1.0 - IRILU);

double T1 = KT1*(A1 + GTOX);
double T2 = KT2*(A2 + GTOX);
double T3 = KT3*(A3 + GTOX);
double T4 = KT4*(A4 + GTOX);
double T5 = KT5*(A5 + GTOX);
double T6 = KT6*(A6 + GTOX);
double T7 = KT7*(A7 + GTOX);
double T8 = KT8*(A8 + GTOX);

double D1 = (KAGE_N + RN*fmax(T1,0.0))*N1;
double D2 = (KAGE_N + RN*fmax(T2,0.0))*N2;
double D3 = (KAGE_N + RN*fmax(T3,0.0))*N3;
double D4 = (KAGE_N + RN*fmax(T4,0.0))*N4;
double D5 = (KAGE_N + RN*fmax(T5,0.0))*N5;
double D6 = (KAGE_N + RN*fmax(T6,0.0))*N6;
double D7 = (KAGE_N + RN*fmax(T7,0.0))*N7;
double D8 = (KAGE_N + RN*fmax(T8,0.0))*N8;

dxdt_N1 = -D1; dxdt_N2 = -D2; dxdt_N3 = -D3; dxdt_N4 = -D4;
dxdt_N5 = -D5; dxdt_N6 = -D6; dxdt_N7 = -D7; dxdt_N8 = -D8;

double DEATH = (D1+D2+D3+D4+D5+D6+D7+D8)/8.0;
// volume-weighted axonal degeneration rate -> the NfL source
double DAXON = VN1*D1 + VN2*D2 + VN3*D3 + VN4*D4
             + VN5*D5 + VN6*D6 + VN7*D7 + VN8*D8;
// ongoing NON-lethal axonal injury in regions the front is still entering:
// this is the term that makes CSF NfL RISE through the symptomatic phase
double AAXON = VN1*A1*N1 + VN2*A2*N2 + VN3*A3*N3 + VN4*A4*N4
             + VN5*A5*N5 + VN6*A6*N6 + VN7*A7*N7 + VN8*A8*N8;

// ===========================================================================
//  3. TAU SYNTHESIS, PTM AND PROTEOSTASIS
// ===========================================================================
dxdt_MRNA = KDEG_M*( MRNAS*(1.0 - IASO)*(1.0 - MKD) - MRNA );

double CLR_T = pow(AUT, 0.5)*(1.0 - 0.50*ACT);
dxdt_TMON = KSYN_T*MRNA - KDEG_T*CLR_T*TMON - KREL_T*TMON;

dxdt_GSK  = KDEG_G*( (1.0 + WCYTO_G*(CYTOK - CYTO0)) - GSK*(1.0 + ITIDE/(1.0-fmin(ITIDE,0.95))) );
dxdt_OGLC = KOGT*(1.0 - OGLC) - KOGA*(1.0 - IOGA)*OGLC;
dxdt_ACT  = KP300*(1.0 - ISALS)*(1.0 - ILM11)*(1.0 - ACT) - KDEAC*ACT;
dxdt_MTS  = KMTS*(MTS_T - MTS);
dxdt_PGRNS= KPG*( (1.0 + EEZE) - PGRNS );

double AUTO_T = fmax(AUTO_MIN,
                 pow(PGRNS,0.5)*(1.0 + EFASU)*(1.0 - WAGG_LYS*ABAR));
dxdt_AUTO = KAU*(AUTO_T - AUTO);

// ===========================================================================
//  4. EXTRACELLULAR TAU -- the exhaust.  Two species, two Kd's (AXIS 3)
//     Complexes are cleared 4x faster than free tau, so TOTAL tau falls too;
//     the assay that reported -98% measured the FREE fraction of ETN.
// ===========================================================================
double SRC_N = KRELN*(TMON/TMREF + WAGG_N*ABAR) + KDEATH_T*DEATH;
double SRC_S = KRELS*(TMON/TMREF + WAGG_S*ABAR) + 0.02*KDEATH_T*DEATH;

dxdt_ETN = SRC_N - KOUT_N*FRN*ETN - KOUT_CX*(1.0-FRN)*ETN;
dxdt_ETS = SRC_S - KOUT_S*FRS*ETS - KOUT_CX*(1.0-FRS)*ETS;
dxdt_CTN = KISF_CSF*RVOL*ETN - KCSF_OUT*CTN;
dxdt_CTS = KISF_CSF*RVOL*ETS - KCSF_OUT*CTS;

// ===========================================================================
//  5. GLIA AND INNATE IMMUNITY (amplifier, not initiator)
// ===========================================================================
dxdt_MGA  = KMG*( MG0 + WMG_A*ABAR + WMG_D*DEATH + WMG_IF*(IFNA-IFN0) - MGA );
dxdt_CYTOK= KCY*( CYTO0 + WCY_MG*(MGA - MG0) - CYTOK );
dxdt_ASTR = KAS*( AS0 + WAS_C*(CYTOK - CYTO0) + WAS_A*ABAR - ASTR );
dxdt_IFNA = KIF*( (IFN0 + WIF_A*ABAR)*(1.0 - ITPN) - IFNA );
dxdt_LPOX = KLP*( LPO0 + WLP_A*ABAR - LPOX );

// ===========================================================================
//  6. BIOMARKERS
// ===========================================================================
dxdt_NFLC  = NFL_IN + KNFL_S*DAXON + KNFL_A*AAXON - KNFL_O*NFLC;
dxdt_NFLP  = KP_NFL*NFLC - KNFLP_O*NFLP;
dxdt_GFAPP = GFAP_IN + KGF*(ASTR - AS0) - KGF_O*GFAPP;
dxdt_MIDA  = KMIDA*( MIDA_MX*(MIDA_FL + (1.0-MIDA_FL)*(0.5*N4 + 0.3*N1 + 0.2*N5)) - MIDA );
dxdt_SCPW  = KSCP *( SCP_MX *(SCP_FL  + (1.0-SCP_FL )*N6) - SCPW );

// ===========================================================================
//  7. CLINICAL SCORE -- reserve gating + two transit delays (AXIS 7)
// ===========================================================================
double R1 = resv(1.0-N1, TH1*TSC, HL1*HSC);
double R2 = resv(1.0-N2, TH2*TSC, HL2*HSC);
double R3 = resv(1.0-N3, TH3*TSC, HL3*HSC);
double R4 = resv(1.0-N4, TH4*TSC, HL4*HSC);
double R5 = resv(1.0-N5, TH5*TSC, HL5*HSC);
double R6 = resv(1.0-N6, TH6*TSC, HL6*HSC);
double R7 = resv(1.0-N7, TH7*TSC, HL7*HSC);
double R8 = resv(1.0-N8, TH8*TSC, HL8*HSC);

// symptomatic agents act on the SCORE, not on the pathology
double SYM_LD = emx(EMX_LDP, LDPC/V_LDP, IC5_LDP);   // levodopa -> nigral part
double SYM_ZP = emx(EMX_ZOL, ZOLC/V_ZOL, IC5_ZOL);   // zolpidem -> pallidal part

double RAW = WS1*R1 + WS2*R2*(1.0-SYM_ZP) + WS3*R3*(1.0-SYM_LD) + WS4*R4
           + WS5*R5 + WS6*R6 + WS7*R7 + WS8*R8;

dxdt_LAG1 = KLAG*(RAW  - LAG1);
dxdt_LAG2 = KLAG*(LAG1 - LAG2);

dxdt_DYS  = KDYS*( 4.0*R8 - DYS );

// ===========================================================================
//  8. SURVIVAL -- competing risks written properly:  dCIF/dt = h * S
// ===========================================================================
double PSPRSn = fmin(PSPRS_MX, LAG2);
double HPN = HP0*exp(BDYS*DYS)*( (PEG > 0.5) ? PEG_RR : 1.0 );
double HOT = HO0*exp(BPSP*PSPRSn/50.0);
dxdt_CIFP = HPN*SURV;
dxdt_CIFO = HOT*SURV;
dxdt_SURV = -(HPN + HOT)*SURV;

// ===========================================================================
//  9. PHARMACOKINETICS
// ===========================================================================
// anti-tau mAb: 2-compartment IgG + ISF/CSF equilibration at Kp = 0.15%
dxdt_MABC = -(CL_AB/V1_AB)*MABC - (Q_AB/V1_AB)*MABC + (Q_AB/V2_AB)*MABP;
dxdt_MABP =  (Q_AB/V1_AB)*MABC - (Q_AB/V2_AB)*MABP;
// mg/L -> nM for a 150 kDa IgG: x1e6/150000 = x6.6667
double AISF_T = KP_AB*(CP_AB*1e6/(MW_AB*1000.0)) + KVAX_AB*VTIT;
dxdt_AISF = KBBB*(AISF_T - AISF);
dxdt_ACSF = KCSFE*(AISF - ACSF);
dxdt_VTIT = -KVAX_D*VTIT;

// MAPT ASO: intrathecal, months-long tissue residence
dxdt_ASOC = -(KSC_ASO + KCL_ASO)*ASOC;
dxdt_ASOT =  KSC_ASO*ASOC - KTIS_ASO*ASOT;

// small molecules
dxdt_TIDG = -KA_TID*TIDG;
dxdt_TIDC =  F_TID*KA_TID*TIDG - KE_TID*TIDC;
dxdt_OGAG = -KA_OGA*OGAG;
dxdt_OGAC =  F_OGA*KA_OGA*OGAG - KE_OGA*OGAC;
dxdt_DAVD = -KA_DAV*DAVD;
dxdt_DAVC =  F_DAV*KA_DAV*DAVD - KE_DAV*DAVC;
dxdt_SALG = -KA_SAL*SALG;
dxdt_SALC =  F_SAL*KA_SAL*SALG - KE_SAL*SALC;
dxdt_FASC = -KE_FAS*FASC;
dxdt_TPNG = -KA_TPN*TPNG;
dxdt_TPNC =  F_TPN*KA_TPN*TPNG - KE_TPN*TPNC;
dxdt_LMG  = -KA_LM*LMG;
dxdt_LMC  =  F_LM*KA_LM*LMG - KE_LM*LMC;
dxdt_EZEG = -KA_EZE*EZEG;
dxdt_EZEC =  F_EZE*KA_EZE*EZEG - KE_EZE*EZEC;
dxdt_RILG = -KA_RIL*RILG;
dxdt_RILC =  F_RIL*KA_RIL*RILG - KE_RIL*RILC;
dxdt_LDPG = -KA_LDP*LDPG;
dxdt_LDPC =  F_LDP*KA_LDP*LDPG - KE_LDP*LDPC;
dxdt_ZOLC = -KE_ZOL*ZOLC;

$TABLE
// ---- clinical -------------------------------------------------------------
double PSPRS  = fmin(PSPRS_MX, LAG2);
// mPSPRS-10 / 15-item modification: a lower-noise subset dominated by the
// gait/midline, bulbar and oculomotor items (rescaled to the same range)
double MP10   = fmin(100.0, LAG2*1.02);
double SEADL  = fmax(0.0, 100.0 - 1.05*PSPRS);
// vertical saccade peak velocity: an essentially UNTHRESHOLDED readout of
// riMLF burst-neuron survival, which is why it is the earliest sign and why
// it floors while the disease continues (AXIS 5)
double VSV    = fmax(VSV_FL, VSV0*pow(fmax(N4,0.0), GAM_SAC));
double MRPI   = (PONS_A/MIDA)*(MCP_W/SCPW);

// ---- fluid biomarkers -----------------------------------------------------
double FRN_C  = 1.0/(1.0 + ACSF/KD_NTAB);
double CSF_NUB = CTN*FRN_C;          // "CSF unbound N-terminal tau" (PASSPORT)
double CSF_TT  = CTN + CTS;          // CSF total tau (mid-region assay)
double CSF_P181= 0.62*CTN + 1.1*CTS; // CSF p-tau181 (phospho-epitope weighted)
double CSF_MTBR= CTS;                // MTBR-tau243 tracks the causal species
double CSF_SEED= CTS*(1.0/(1.0 + ACSF/KD_STAB));  // seeding activity
double AB_ISF  = AISF;

// ---- pathology summaries --------------------------------------------------
double ATOT   = (A1+A2+A3+A4+A5+A6+A7+A8)/8.0;
double NTOT   = (N1+N2+N3+N4+N5+N6+N7+N8)/8.0;
double F_EXT  = ETN/TMON;            // the pool-size ratio (AXIS 1)
double KNOCK  = 100.0*(1.0 - TMON/TMREF);

$CAPTURE @annotated
PSPRS  : PSPRS total score (0-100)
MP10   : mPSPRS-10 modified score
SEADL  : Schwab & England ADL (%)
VSV    : vertical saccade peak velocity (deg/s)
MRPI   : MR parkinsonism index
CSF_NUB: CSF unbound N-terminal tau (nM)
CSF_TT : CSF total tau (nM)
CSF_P181: CSF p-tau181 (nM-equivalents)
CSF_MTBR: CSF MTBR-tau243 (nM)
CSF_SEED: CSF seeding activity (nM-equivalents)
AB_ISF : ISF antibody concentration (nM)
ATOT   : mean regional aggregate load
NTOT   : mean regional surviving-neuron fraction
F_EXT  : ISF tau / intraneuronal tau (pool-size ratio)
KNOCK  : monomer knockdown (%)
)---"

psp_mod <- mrgsolve::mcode_cache("psp_qsp", psp_code, soloc = tempdir())

# =============================================================================
#  SCENARIO LIBRARY
# =============================================================================
#  All trial scenarios enrol at TENROL (default 3.5 y from biological onset),
#  which is where the model reproduces the observed enrolment PSPRS of ~38.
# =============================================================================

# TENROL is NOT assumed: it is the time at which the natural-history run first
# reaches the observed PASSPORT/ARISE enrolment severity of PSPRS = 38.  Fitting
# it rather than fixing it turns the pre-clinical duration into a model OUTPUT:
# ~7.0 y of tauopathy precede trial-entry severity, of which ~3.5 y precede any
# symptom at all.  Every trial scenario therefore starts here.
TENROL <- 2583      # days from biological onset to trial enrolment (7.07 y)

# Published survival is counted from SYMPTOM onset, not from the first tau
# aggregate.  PASSPORT/ARISE enrolled a median ~3.5 y after symptom onset, so
# the model's implied symptom onset is TENROL - 3.5 y.  Everything reported
# "from onset" below is therefore reported from TSYMPT, and the ~3.6 y of
# pre-symptomatic tauopathy that precedes it is a model OUTPUT, not an input.
TSYMPT <- TENROL - 1278   # 1305 d = 3.57 y of pre-symptomatic disease
WK52   <- 364       # one year of treatment
WK72   <- 504       # Preserve trial duration

.ev <- function(...) mrgsolve::ev(...)

## ---- helper: repeated dosing block ----------------------------------------
dose_block <- function(cmt, amt, ii, n, start = TENROL) {
  mrgsolve::ev(time = start, amt = amt, cmt = cmt, ii = ii, addl = max(0, n - 1))
}

psp_scenarios <- function() {
  list(

    # --- natural history ---------------------------------------------------
    S00 = list(label = "S00 | Natural history, PSP-RS (placebo)",
               par = list(), ev = NULL),

    S01 = list(label = "S01 | Natural history, PSP-P (nigral wave origin)",
               par = list(ORIG_RS = 0), ev = NULL),

    # --- anti-tau monoclonal antibodies -----------------------------------
    S02 = list(label = "S02 | Gosuranemab 2000 mg IV q4w x 52 wk (PASSPORT)",
               par = list(KD_NTAB = 0.0344, KD_STAB = 1e6),
               ev  = dose_block("MABC", 2000, 28, 13)),

    S03 = list(label = "S03 | Tilavonemab 4000 mg IV q4w x 52 wk (ARISE)",
               par = list(KD_NTAB = 0.0450, KD_STAB = 1e6),
               ev  = dose_block("MABC", 4000, 28, 13)),

    S04 = list(label = "S04 | Mid-domain/MTBR mAb 2000 mg IV q4w x 52 wk (probe)",
               par = list(KD_NTAB = 1e6, KD_STAB = 0.0344),
               ev  = dose_block("MABC", 2000, 28, 13)),

    S05 = list(label = "S05 | Mid-domain mAb with PHI_ACC = 1 (all transfer accessible)",
               par = list(KD_NTAB = 1e6, KD_STAB = 0.0344, PHI_ACC = 1.0),
               ev  = dose_block("MABC", 2000, 28, 13)),

    S06 = list(label = "S06 | NIO752 ASO, phase-1-like: 4 doses q4w (calibration)",
               par = list(),
               ev  = dose_block("ASOC", 1.9, 28, 4)),

    S07 = list(label = "S07 | NIO752 ASO, Preserve-like: q12w x 72 wk",
               par = list(),
               ev  = dose_block("ASOC", 1.9, 84, 6)),

    S08 = list(label = "S08 | MAPT ASO 4x dose q12w x 72 wk (elasticity probe)",
               par = list(),
               ev  = dose_block("ASOC", 7.6, 84, 6)),

    S09 = list(label = "S09 | MAPT ASO started 3 y before SYMPTOM onset",
               par = list(),
               ev  = dose_block("ASOC", 1.9, 84, 28, start = TSYMPT - 1095)),

    # --- local-kinetic agents: the three structural nulls ------------------
    S10 = list(label = "S10 | Tideglusib 800 mg PO qd x 52 wk (TAUROS)",
               par = list(),
               ev  = dose_block("TIDG", 800, 1, 364)),

    S11 = list(label = "S11 | OGA inhibitor 30 mg PO qd x 52 wk",
               par = list(),
               ev  = dose_block("OGAG", 30, 1, 364)),

    S12 = list(label = "S12 | Davunetide 30 mg intranasal bid x 52 wk (AL-108-231)",
               par = list(),
               ev  = dose_block("DAVD", 30, 0.5, 728)),

    # --- other mechanisms --------------------------------------------------
    S13 = list(label = "S13 | Salsalate 1500 mg PO bid x 52 wk",
               par = list(),
               ev  = dose_block("SALG", 1500, 0.5, 728)),

    S14 = list(label = "S14 | Fasudil 40 mg IV bid x 52 wk (ROCK inhibition)",
               par = list(),
               ev  = dose_block("FASC", 40, 0.5, 728)),

    S15 = list(label = "S15 | TPN-101 400 mg PO qd x 52 wk (LINE-1 RT)",
               par = list(),
               ev  = dose_block("TPNG", 400, 1, 364)),

    S16 = list(label = "S16 | LM11A-31 400 mg PO bid x 52 wk (p75NTR)",
               par = list(),
               ev  = dose_block("LMG", 400, 0.5, 728)),

    S17 = list(label = "S17 | Ezeprogind 60 mg PO bid x 52 wk (progranulin)",
               par = list(),
               ev  = dose_block("EZEG", 60, 0.5, 728)),

    S18 = list(label = "S18 | AADvac1 active vaccine, 6 doses then boosters",
               par = list(KD_NTAB = 0.0344, KD_STAB = 1e6),
               ev  = dose_block("VTIT", 60, 28, 9)),

    S19 = list(label = "S19 | Riluzole 50 mg PO bid x 52 wk (NNIPPS)",
               par = list(),
               ev  = dose_block("RILG", 50, 0.5, 728)),

    # --- symptomatic -------------------------------------------------------
    S20 = list(label = "S20 | Levodopa/carbidopa 250 mg tid x 52 wk",
               par = list(),
               ev  = dose_block("LDPG", 250, 1/3, 1092)),

    S21 = list(label = "S21 | Zolpidem 10 mg qd x 52 wk",
               par = list(),
               ev  = dose_block("ZOLC", 10, 1, 364)),

    # --- combinations ------------------------------------------------------
    S22 = list(label = "S22 | MAPT ASO + mid-domain mAb x 72 wk",
               par = list(KD_NTAB = 1e6, KD_STAB = 0.0344),
               ev  = dose_block("ASOC", 1.9, 84, 6) +
                     dose_block("MABC", 2000, 28, 18)),

    S23 = list(label = "S23 | MAPT ASO + ezeprogind x 72 wk (supply down, clearance up)",
               par = list(),
               ev  = dose_block("ASOC", 1.9, 84, 6) +
                     dose_block("EZEG", 60, 0.5, 1008)),

    # --- in-model counterfactuals: the timing trap ------------------------
    S24 = list(label = "S24 | COUNTERFACTUAL: 100% aggregation block from enrolment",
               par = list(FBLOCK = 1.0, TBLOCK = TENROL), ev = NULL),

    S25 = list(label = "S25 | COUNTERFACTUAL: 100% block 3 y before symptom onset",
               par = list(FBLOCK = 1.0, TBLOCK = TSYMPT - 1095), ev = NULL),

    S26 = list(label = "S26 | PEG gastrostomy at enrolment (hazard modifier only)",
               par = list(PEG = 1), ev = NULL),

    S27 = list(label = "S27 | Triple, started 3 y before symptom onset",
               par = list(KD_NTAB = 1e6, KD_STAB = 0.0344),
               ev  = dose_block("ASOC", 7.6, 84, 32, start = TSYMPT - 1095) +
                     dose_block("MABC", 2000, 28, 96, start = TSYMPT - 1095) +
                     dose_block("EZEG", 60,   0.5, 5840, start = TSYMPT - 1095)),

    S28 = list(label = "S28 | COUNTERFACTUAL: 50x disaggregation from enrolment",
               par = list(FDIS = 50, TDIS = TENROL), ev = NULL)
  )
}

# =============================================================================
#  SIMULATION DRIVERS
# =============================================================================

#' Run one scenario
#' @param id scenario id, e.g. "S02"
#' @param end simulation end time (days from biological onset)
#' @param delta output grid (days)
psp_run <- function(id = "S00", end = 5475, delta = 7, extra_par = list()) {
  sc  <- psp_scenarios()[[id]]
  if (is.null(sc)) stop("unknown scenario: ", id)
  m   <- psp_mod
  par <- modifyList(sc$par, extra_par)
  if (length(par)) m <- mrgsolve::param(m, par)
  out <- if (is.null(sc$ev)) {
    mrgsolve::mrgsim(m, end = end, delta = delta, atol = 1e-10, rtol = 1e-8,
                     maxsteps = 200000)
  } else {
    mrgsolve::mrgsim(m, events = sc$ev, end = end, delta = delta,
                     atol = 1e-10, rtol = 1e-8, maxsteps = 200000)
  }
  d <- as.data.frame(out)
  d$scenario <- id
  d$label    <- sc$label
  d
}

#' Run every scenario and return one long data frame
psp_run_all <- function(end = 5475, delta = 7) {
  ids <- names(psp_scenarios())
  do.call(rbind, lapply(ids, function(i) psp_run(i, end = end, delta = delta)))
}

## ---- small readers ---------------------------------------------------------
.at <- function(d, t, col) {
  i <- which.min(abs(d$time - t)); d[[col]][i]
}
.median_surv <- function(d) {
  i <- which(d$SURV <= 0.5)
  if (!length(i)) return(NA_real_)
  d$time[i[1]]/365.25
}
## time at which a run first reaches a given PSPRS -- used so that PSP-P is
## compared with PSP-RS at MATCHED SEVERITY rather than at matched calendar time
.t_at_psprs <- function(d, target = 38) {
  i <- which(d$PSPRS >= target)
  if (!length(i)) return(NA_real_)
  d$time[i[1]]
}

# =============================================================================
#  VALIDATION HARNESS
# =============================================================================
#  Every row is a number published by someone else.  The model was fitted to a
#  SUBSET (marked FIT); the rest are held out.  Failures are reported, not
#  repaired -- see the FAILURES section of README.md.
# =============================================================================

psp_validate <- function(verbose = TRUE) {

  end <- 5475
  rs  <- psp_run("S00", end = end)          # natural history PSP-RS
  pp  <- psp_run("S01", end = end)          # natural history PSP-P
  go  <- psp_run("S02", end = end)          # gosuranemab
  ti  <- psp_run("S03", end = end)          # tilavonemab
  md  <- psp_run("S04", end = end)          # mid-domain mAb
  a1  <- psp_run("S06", end = end)          # ASO phase-1-like
  pv  <- psp_run("S07", end = end)          # ASO Preserve-like
  tg  <- psp_run("S10", end = end)          # tideglusib
  og  <- psp_run("S11", end = end)          # OGA inhibitor
  dv  <- psp_run("S12", end = end)          # davunetide
  cf1 <- psp_run("S24", end = end)          # counterfactual from enrolment
  cf2 <- psp_run("S25", end = end)          # counterfactual pre-symptomatic

  T0 <- TENROL; T1 <- TENROL + WK52

  ch <- function(d) .at(d, T1, "PSPRS") - .at(d, T0, "PSPRS")

  rows <- list()
  add <- function(what, obs, pred, unit = "", tag = "held-out") {
    rows[[length(rows)+1]] <<- data.frame(
      quantity = what, observed = obs, predicted = pred, unit = unit,
      logratio = if (is.na(obs) || is.na(pred) || obs <= 0 || pred <= 0) NA_real_
                 else abs(log10(pred/obs)),
      status = tag, stringsAsFactors = FALSE)
  }

  ## --- disease progression -------------------------------------------------
  add("PSPRS at enrolment (3.5 y after symptom onset)", 38.0, .at(rs, T0, "PSPRS"), "pts", "FIT")
  add("PSPRS 52-wk change, placebo (PASSPORT)", 10.6, ch(rs), "pts", "FIT")
  add("Natural-history slope (Golbe 2007)", 11.3,
      (.at(rs, T0+730, "PSPRS") - .at(rs, T0, "PSPRS"))/2, "pts/y", "held-out")
  ## PSP-P is compared at MATCHED SEVERITY (its own PSPRS=38 crossing), which is
  ## how the observational cohorts report it -- not at matched calendar time.
  tp <- .t_at_psprs(pp, 38)
  add("PSP-P PSPRS slope at matched severity", 5.0,
      .at(pp, tp + WK52, "PSPRS") - .at(pp, tp, "PSPRS"), "pts/y", "held-out")
  add("Median survival from onset, PSP-RS", 7.9,
      .median_surv(rs) - TSYMPT/365.25, "y", "FIT")
  add("Median survival from onset, PSP-P", 10.0,
      .median_surv(pp) - (tp - 1278)/365.25, "y", "held-out")
  add("PSP-P delay to trial-entry severity vs PSP-RS", NA,
      (tp - TENROL)/365.25, "y", "prediction")

  ## --- the two target-engagement numbers -----------------------------------
  gn0 <- .at(rs, T1, "CSF_NUB"); gn1 <- .at(go, T1, "CSF_NUB")
  add("Gosuranemab CSF unbound N-terminal tau change", -98.0,
      100*(gn1/gn0 - 1), "%", "FIT")
  add("Gosuranemab PSPRS 52-wk change", 10.4, ch(go), "pts", "held-out")
  add("Tilavonemab PSPRS 52-wk change (futility)", 10.6, ch(ti), "pts", "held-out")

  at0 <- .at(rs, T0+112, "CSF_TT"); at1 <- .at(a1, T0+112, "CSF_TT")
  add("NIO752 CSF total tau change (phase 1)", -20.0,
      100*(at1/at0 - 1), "%", "FIT")
  ap0 <- .at(rs, T0+112, "CSF_P181"); ap1 <- .at(a1, T0+112, "CSF_P181")
  add("NIO752 CSF p-tau181 change", -20.0, 100*(ap1/ap0 - 1), "%", "held-out")

  n0 <- .at(rs, T0, "NFLC")
  add("CSF NfL rise over 1 y, placebo", 40.0,
      100*(.at(rs, T1, "NFLC")/n0 - 1), "%", "FIT")
  add("CSF NfL rise over 1 y, ASO-treated", 0.0,
      100*(.at(a1, T1, "NFLC")/n0 - 1), "%", "held-out")

  ## --- the three structural nulls ------------------------------------------
  add("Tideglusib PSPRS 52-wk change (TAUROS null)", 10.6, ch(tg), "pts", "held-out")
  add("Davunetide PSPRS 52-wk change (AL-108-231 null)", 11.0, ch(dv), "pts", "held-out")
  add("OGA inhibitor PSPRS 52-wk change (predicted null)", NA, ch(og), "pts", "prediction")

  ## --- imaging and oculomotor ----------------------------------------------
  add("Midbrain area at enrolment", 78.0, .at(rs, T0, "MIDA"), "mm2", "held-out")
  add("MRPI at enrolment (cut-off 13.55)", 15.0, .at(rs, T0, "MRPI"), "-", "held-out")
  add("Vertical saccade peak velocity at enrolment", 175.0, .at(rs, T0, "VSV"), "deg/s", "held-out")

  ## --- the pool-size ratio that carries AXIS 1 -----------------------------
  add("ISF tau / intraneuronal tau (F_EXT)", 5.0e-4, .at(rs, T0, "F_EXT"), "-", "structural")
  add("ISF antibody concentration at steady state (expect 1-10 nM)", NA,
      .at(go, T1, "AB_ISF"), "nM", "structural")
  add("Gosuranemab CSF TOTAL tau change (sink effect, expect UP)", NA,
      100*(.at(go, T1, "CSF_TT")/.at(rs, T1, "CSF_TT") - 1), "%", "structural")

  ## --- the timing trap (in-model counterfactuals) --------------------------
  add("Perfect block from enrolment: 52-wk PSPRS benefit", NA,
      ch(rs) - ch(cf1), "pts", "prediction")
  add("Perfect block 3 y before symptom onset: PSPRS benefit at enrol+52wk", NA,
      .at(rs, T1, "PSPRS") - .at(cf2, T1, "PSPRS"), "pts", "prediction")
  add("Mid-domain mAb: 52-wk PSPRS benefit", NA, ch(rs) - ch(md), "pts", "prediction")
  cl <- psp_run("S28", end = end)
  add("50x DISAGGREGATION from enrolment: 52-wk PSPRS benefit", NA,
      ch(rs) - ch(cl), "pts", "prediction")
  add("ASO Preserve-like: 72-wk PSPRS benefit", NA,
      (.at(rs, TENROL+WK72, "PSPRS") - .at(rs, T0, "PSPRS")) -
      (.at(pv, TENROL+WK72, "PSPRS") - .at(pv, T0, "PSPRS")), "pts", "prediction")

  res <- do.call(rbind, rows)
  if (verbose) {
    cat("\n=== PSP QSP model -- validation against published anchors ===\n\n")
    pr <- res
    pr$observed  <- ifelse(is.na(pr$observed),  "--", formatC(pr$observed,  format = "g", digits = 4))
    pr$predicted <- formatC(pr$predicted, format = "g", digits = 4)
    pr$logratio  <- ifelse(is.na(pr$logratio), "--", formatC(pr$logratio, format = "f", digits = 3))
    print(pr, row.names = FALSE, right = FALSE)
    lr <- res$logratio[!is.na(res$logratio)]
    cat("\nmedian |log10(pred/obs)| over", length(lr), "quantitative anchors:",
        formatC(median(lr), format = "f", digits = 3), "\n\n")
  }
  invisible(res)
}

# =============================================================================
#  ELASTICITY OF THE CLINICAL SLOPE TO MONOMER KNOCKDOWN
# =============================================================================
#  The single number that decides whether an ASO dose is sufficient.
#  Continuum Fisher-KPP predicts e = 0.5; a discrete pulled front on a sparse
#  graph predicts e -> 1.  The model MEASURES it rather than assuming either.
# =============================================================================

psp_elasticity <- function(kd = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),
                           horizon = WK72, verbose = TRUE) {
  Y <- 365.25
  runkd <- function(k, tkd, end) {
    m <- mrgsolve::param(psp_mod, MAPT_KD = k/100, TKD = tkd, ASO_DEEP = 1.0)
    as.data.frame(mrgsolve::mrgsim(m, end = end, delta = 7, atol = 1e-10,
                                   rtol = 1e-8, maxsteps = 200000))
  }
  ## ---- (a) TREATMENT elasticity ------------------------------------------
  ##  knockdown imposed AT ENROLMENT, slope over the following `horizon`.
  ##  This is the quantity the Preserve trial will actually measure.
  ## ---- (b) PREVENTION elasticity -----------------------------------------
  ##  knockdown imposed from biological onset; read out as the DELAY in reaching
  ##  trial-entry severity (PSPRS = 38).  This is the quantity that matters to
  ##  the patient, and the two are not remotely the same number.
  out <- data.frame()
  b_tx <- runkd(0, TENROL, TENROL + horizon + 7)
  s0   <- .at(b_tx, TENROL + horizon, "PSPRS") - .at(b_tx, TENROL, "PSPRS")
  for (k in kd) {
    dtx <- runkd(k, TENROL, TENROL + horizon + 7)
    sl  <- .at(dtx, TENROL + horizon, "PSPRS") - .at(dtx, TENROL, "PSPRS")
    dpv <- runkd(k, 0, 12000)
    i   <- which(dpv$PSPRS >= 38)
    t38 <- if (length(i)) dpv$time[i[1]] else NA_real_
    out <- rbind(out, data.frame(
      knockdown_pct = k,
      tx_slope      = round(sl, 3),
      tx_rel        = round(sl/s0, 4),
      prev_t38_y    = round(t38/Y, 2),
      prev_delay_y  = round((t38 - TENROL)/Y, 2)))
  }
  el <- function(v) {
    ok <- out$knockdown_pct %in% c(10, 30)
    x  <- log(1 - out$knockdown_pct[ok]/100); y <- log(v[ok])
    diff(y)/diff(x)
  }
  e_tx <- el(out$tx_rel)
  if (verbose) {
    cat("\n=== Elasticity of PSP progression to MAPT knockdown ===\n")
    cat("tx_*   : knockdown STARTED AT ENROLMENT (what Preserve measures)\n")
    cat("prev_* : knockdown from biological onset (what matters to the patient)\n\n")
    print(out, row.names = FALSE)
    cat("\ntreatment-phase elasticity dln(slope)/dln(M) around 20% knockdown:",
        formatC(e_tx, format = "f", digits = 3), "\n")
    cat("continuum Fisher-KPP predicts 0.50 | discrete pulled front predicts ~1.0\n")
    hv <- try(approx(out$tx_rel, out$knockdown_pct, xout = 0.5)$y, silent = TRUE)
    cat("knockdown needed to halve the TREATMENT-phase slope:",
        if (inherits(hv, "try-error") || is.na(hv)) "not achievable at any knockdown"
        else paste0(formatC(hv, format = "f", digits = 1), "%"), "\n")
    cat("delay in reaching trial-entry severity at 50% lifelong knockdown:",
        formatC(out$prev_delay_y[out$knockdown_pct == 50], format = "f", digits = 2),
        "years\n\n")
  }
  invisible(list(table = out, elasticity_treatment = e_tx))
}

# =============================================================================
#  SCENARIO COMPARISON TABLE
# =============================================================================

psp_compare <- function(horizon = WK52, verbose = TRUE) {
  ids  <- names(psp_scenarios())
  base <- psp_run("S00", end = TENROL + horizon + 7)
  b0   <- .at(base, TENROL, "PSPRS")
  bch  <- .at(base, TENROL + horizon, "PSPRS") - b0
  rows <- do.call(rbind, lapply(ids, function(i) {
    d  <- psp_run(i, end = 5475)
    ch <- .at(d, TENROL + horizon, "PSPRS") - .at(d, TENROL, "PSPRS")
    ## Change-from-baseline FLATTERS the placebo arm whenever an intervention
    ## delays the disease: a delayed patient sits lower on the curve at the
    ## fixed calendar enrolment time and therefore has MORE room to change.
    ## The absolute score is reported alongside it so the artefact is visible.
    data.frame(
      id       = i,
      label    = substr(psp_scenarios()[[i]]$label, 1, 58),
      dPSPRS   = round(ch, 2),
      vs_pbo   = round(ch - bch, 2),
      absPSPRS = round(.at(d, TENROL + horizon, "PSPRS"), 2),
      abs_vs_pbo = round(.at(d, TENROL + horizon, "PSPRS") -
                         .at(base, TENROL + horizon, "PSPRS"), 2),
      CSF_NUB  = signif(.at(d, TENROL + horizon, "CSF_NUB"), 3),
      CSF_TT   = signif(.at(d, TENROL + horizon, "CSF_TT"), 3),
      NfL_CSF  = round(.at(d, TENROL + horizon, "NFLC"), 0),
      VSV      = round(.at(d, TENROL + horizon, "VSV"), 0),
      medSurv  = round(.median_surv(d), 2),
      stringsAsFactors = FALSE)
  }))
  if (verbose) {
    cat("\n=== PSP scenario comparison at", horizon, "days of treatment ===\n")
    cat("dPSPRS/vs_pbo  = change from enrolment (the trial endpoint)\n")
    cat("absPSPRS/abs_* = absolute score (what the patient actually has)\n\n")
    print(rows, row.names = FALSE, right = FALSE)
    cat("\n")
  }
  invisible(rows)
}

# =============================================================================
#  IF RUN NON-INTERACTIVELY: print the whole numerical suite
# =============================================================================
if (!interactive() && identical(environment(), globalenv())) {
  if (isTRUE(as.logical(Sys.getenv("PSP_RUN_SUITE", "FALSE")))) {
    psp_validate()
    psp_elasticity()
    psp_compare()
  }
}
