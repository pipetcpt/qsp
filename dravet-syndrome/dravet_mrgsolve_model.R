# =============================================================================
#  Dravet Syndrome (SCN1A Developmental & Epileptic Encephalopathy)
#  Quantitative Systems Pharmacology model — mrgsolve implementation
# =============================================================================
#
#  44 ODEs. Drug PK for clobazam (+ the active metabolite norclobazam),
#  stiripentol, cannabidiol (+ 7-OH-CBD), fenfluramine (+ norfenfluramine),
#  valproate, a lamotrigine-like sodium-channel blocker, rescue diazepam and
#  an intrathecal SCN1A-directed antisense oligonucleotide; a CYP turnover
#  layer that carries the drug-drug interactions; target biology from SCN1A
#  transcript through Nav1.1 to interneuron firing capacity; an E:I balance
#  that generates the seizure hazard; and clinical read-outs for seizure
#  frequency, status epilepticus, development, somnolence, weight,
#  transaminases, valve exposure and SUDEP.
#
#  ---------------------------------------------------------------------------
#  THE FOUR STRUCTURAL COMMITMENTS
#
#  (1) The lesion is in INHIBITION. Nav1.1 is expressed preferentially in
#      GABAergic interneurons. The excitatory compartment runs on Nav1.2/1.6
#      and is left intact. ALLELE = 0.5 in Dravet, 1.0 in a healthy control;
#      that single number is the disease.
#
#  (2) Interneuron firing capacity is a STEEP threshold function of sodium
#      reserve (HILL_INT = 6), pyramidal excitability a SHALLOW one
#      (HILL_PYR = 1.5). The healthy operating point sits far above threshold,
#      the Dravet operating point just above it. Consequence: sodium-channel
#      blockers are anticonvulsant in a normal brain and proconvulsant in
#      Dravet. That sign flip is an OUTPUT. Nothing in this file tests for
#      "Dravet" when computing a drug effect.
#
#  (3) CLOBAZAM IS TWO DRUGS. Norclobazam is cleared by polymorphic CYP2C19
#      and reaches 5-20x the parent, so every CYP2C19 inhibitor is an indirect
#      GABAergic drug. Stiripentol and cannabidiol both inhibit CYP2C19 AND
#      have direct targets. PK_ROUTE and PD_ROUTE switch the two paths
#      independently, which is how the decomposition is done.
#      Stiripentol's own GABA-A action is given a SEPARATE saturable term
#      (GAIN_STP) because it binds an alpha3-preferring site, not the
#      benzodiazepine site. Folding them together would make any GABAergic
#      add-on structurally unable to help a patient already on clobazam.
#
#  (4) EFFICACY AND SEDATION SHARE THE NORCLOBAZAM NODE, so the PK route and
#      the direct route have different therapeutic indices.
#
#  ---------------------------------------------------------------------------
#  CALIBRATION (see dravet_references.md for the PMIDs)
#
#  Fitted, 8 parameters, each to a named published quantity:
#    KI_CBD_2C19, KI_7OH_2C19  <- norclobazam x6.0 on cannabidiol (Geffrey 2015)
#    KI_CBD_3A4                <- clobazam parent x1.6            (Geffrey 2015)
#    KI_STP_2C19               <- norclobazam x2.5 on stiripentol (Jullien 2015)
#    EMAX_PAM                  <- clobazam withdrawal doubles seizure frequency
#    EMAX_CBD                  <- GWPCARE1/2, cannabidiol 20 mg/kg/day
#    EMAX_5HT                  <- Study 1, fenfluramine 0.7 mg/kg/day
#    EMAX_STP_SITE             <- STICLO, stiripentol 50 mg/kg/day
#
#  Pinned by an INDEPENDENT observation and not fitted to any of the drugs
#  being decomposed:
#    EC50_PAM_NCLB             <- Hashi 2015: norclobazam near 1100 ng/mL
#                                 already associates with >=90% seizure
#                                 control, i.e. the benzodiazepine-site
#                                 response is near-saturated BELOW the
#                                 concentration a standard clobazam dose
#                                 reaches. The central conclusion of this
#                                 model rests on this one number.
#
#  Reproduced without being fitted: fenfluramine 0.2 mg/kg/day, cannabidiol
#  10 mg/kg/day, Study 1504 under trial-faithful enrichment, the
#  sodium-channel-blocker sign flip, the febrile susceptibility ratio, the
#  antisense ceiling, and the CYP2C19 genotype-independence of stiripentol
#  response that Kouga 2015 observed.
#
#  A dependency-free Python twin of this model, dravet_reference_impl.py,
#  reproduces every number quoted in README.md.
#
#  Units: time in DAYS. Amounts in mg/kg, volumes in L/kg, clearances in
#  L/day/kg, so concentrations come out in mg/L. Doses are mg/kg per
#  administration.
# =============================================================================

library(mrgsolve)

dravet_code <- '
$PARAM @annotated
// ---------------- clobazam -------------------------------------------------
F_CLB    :  0.90 : Clobazam oral bioavailability
KA_CLB   : 12.0  : Clobazam absorption rate (1/day)
V_CLB    :  1.50 : Clobazam central volume (L/kg)
CL_CLB   :  2.25 : Clobazam total clearance (L/day/kg)
Q_CLB    :  1.10 : Clobazam intercompartmental clearance (L/day/kg)
VP_CLB   :  1.30 : Clobazam peripheral volume (L/kg)
FM_NCLB  :  0.85 : Fraction of clobazam demethylated to norclobazam
V_NCLB   :  0.863 : Norclobazam volume (L/kg)
CL_NCLB  :  0.239 : Norclobazam CYP2C19 clearance (L/day/kg)
Q_NCLB   :  0.30 : Norclobazam intercompartmental clearance (L/day/kg)
VP_NCLB  :  1.00 : Norclobazam peripheral volume (L/kg)
// ---------------- stiripentol (Michaelis-Menten) ---------------------------
F_STP    :  0.90 : Stiripentol bioavailability
KA_STP   : 10.0  : Stiripentol absorption rate (1/day)
V_STP    :  1.50 : Stiripentol volume (L/kg)
VMAX_STP : 90.0  : Stiripentol Vmax (mg/kg/day)
KM_STP   : 10.0  : Stiripentol Km (mg/L)
// ---------------- cannabidiol ---------------------------------------------
F_CBD    :  0.10 : Cannabidiol bioavailability
KA_CBD   :  8.0  : Cannabidiol absorption rate (1/day)
V_CBD    : 20.0  : Cannabidiol central volume (L/kg)
CL_CBD   :  5.70 : Cannabidiol clearance (L/day/kg)
Q_CBD    :  6.0  : Cannabidiol intercompartmental clearance (L/day/kg)
VP_CBD   : 25.0  : Cannabidiol peripheral volume (L/kg)
FM_7OH   :  0.40 : Fraction of cannabidiol to 7-OH-CBD
V_7OH    :  6.0  : 7-OH-CBD volume (L/kg)
CL_7OH   :  2.67 : 7-OH-CBD clearance (L/day/kg)
V_7COOH  :  4.0  : 7-COOH-CBD volume (L/kg)
CL_7COOH :  3.0  : 7-COOH-CBD clearance (L/day/kg)
// ---------------- fenfluramine --------------------------------------------
F_FFA    :  0.75 : Fenfluramine bioavailability
KA_FFA   : 10.0  : Fenfluramine absorption rate (1/day)
V_FFA    : 11.5  : Fenfluramine central volume (L/kg)
CL_FFA   :  9.55 : Fenfluramine clearance (L/day/kg)
Q_FFA    :  3.0  : Fenfluramine intercompartmental clearance (L/day/kg)
VP_FFA   : 14.0  : Fenfluramine peripheral volume (L/kg)
FM_NOR   :  0.35 : Fraction to norfenfluramine
V_NOR    : 15.5  : Norfenfluramine volume (L/kg)
CL_NOR   :  7.35 : Norfenfluramine clearance (L/day/kg)
FRAC_CYP_FFA : 0.70 : Share of fenfluramine clearance via CYP1A2/2B6/2D6
// ---------------- valproate ----------------------------------------------
F_VPA    :  1.00 : Valproate bioavailability
KA_VPA   :  8.0  : Valproate absorption rate (1/day)
V_VPA    :  0.31 : Valproate volume (L/kg)
CL_VPA   :  0.43 : Valproate clearance (L/day/kg)
// ---------------- sodium-channel blocker (lamotrigine-like) --------------
F_NVB    :  0.98 : Blocker bioavailability
KA_NVB   :  8.0  : Blocker absorption rate (1/day)
V_NVB    :  1.77 : Blocker volume (L/kg)
CL_NVB   :  1.22 : Blocker clearance (L/day/kg)
FRAC_UGT_NVB : 0.75 : Share of blocker clearance via UGT1A4
// ---------------- rescue benzodiazepine ----------------------------------
V_DZP    :  2.0  : Diazepam volume (L/kg)
CL_DZP   :  1.4  : Diazepam clearance (L/day/kg)
// ---------------- enzyme turnover & inhibition ---------------------------
KE_ENZ      : 0.50   : Enzyme turnover rate (1/day)
KI_STP_2C19 : 4.956  : Stiripentol Ki, CYP2C19 (mg/L)
KI_CBD_2C19 : 0.0742 : Cannabidiol Ki, CYP2C19 (mg/L)
KI_7OH_2C19 : 0.1113 : 7-OH-CBD Ki, CYP2C19 (mg/L)
KI_STP_3A4  : 8.0    : Stiripentol Ki, CYP3A4 (mg/L)
KI_CBD_3A4  : 0.671  : Cannabidiol Ki, CYP3A4 (mg/L)
KI_STP_1A2  : 12.0   : Stiripentol Ki, CYP1A2/2B6/2D6 (mg/L)
KI_VPA_UGT  : 45.0   : Valproate Ki, UGT1A4 (mg/L)
KI_CBD_UGT  : 0.60   : Cannabidiol Ki, UGT (mg/L)
E2C19_GENO  : 1.00   : CYP2C19 activity: UM 1.6 NM 1.0 IM 0.55 PM 0.20
// ---------------- SCN1A / Nav1.1 -----------------------------------------
ALLELE   : 0.50 : Functional SCN1A allele dose (0.5 Dravet 1.0 control)
PE_FRAC  : 0.33 : Poison-exon 20N share of transcripts
ASO_EMAX : 0.85 : Maximum achievable poison-exon skipping
ASO_EC50 : 25.0 : Brain ASO amount for half-maximal skipping
K_NAV    : 0.50 : Nav1.1 turnover (1/day)
SENS_T   : 0.110 : Mutant-selective thermal loss of Nav1.1 per degC
SENS_T_WT : 0.010 : Wild-type thermal sensitivity per degC
// ---------------- the interneuron / pyramidal asymmetry -----------------
EC50_INT : 0.45 : Interneuron sodium reserve at half-maximal firing capacity
HILL_INT : 6.0  : Interneuron threshold steepness
EC50_PYR : 0.35 : Pyramidal sodium reserve at half-maximal excitability
HILL_PYR : 1.5  : Pyramidal steepness
IC50_NVB_INT : 9.33 : Blocker IC50 at interneuron Nav (fast-spiking)
IC50_NVB_PYR : 37.3 : Blocker IC50 at pyramidal Nav
IC50_NVB_PROP : 3.0 : Blocker IC50 for propagation block
EMAX_NVB_PROP : 0.55 : Maximum propagation block
// ---------------- GABA-A pharmacology -----------------------------------
EC50_PAM_CLB  : 0.25 : Clobazam potency at the benzodiazepine site (mg/L)
EC50_PAM_NCLB : 0.45 : Norclobazam potency (mg/L) -- Hashi 2015
EC50_PAM_DZP  : 0.30 : Diazepam potency (mg/L)
EMAX_PAM      : 0.65396 : Maximum benzodiazepine-site inhibitory gain
EC50_PAM_STP  : 10.0 : Stiripentol potency at its alpha3 site (mg/L)
EMAX_STP_SITE : 1.82397 : Maximum stiripentol-site inhibitory gain
FTOL     : 0.35 : Maximum fractional receptor tolerance
KTOL     : 0.01155 : Tolerance rate (1/day), t1/2 ~60 d
PAM_TOL  : 2.0  : PAM load for half-maximal tolerance
EMAX_VPA : 0.45 : Maximum valproate inhibitory gain
EC50_VPA : 60.0 : Valproate EC50 (mg/L)
// ---------------- non-GABAergic protection ------------------------------
EMAX_5HT : 6.37925 : Maximum serotonergic protection
EC50_5HT : 1.00 : Serotonergic EC50 (mg/L); >> Css, so response is linear
W_NORFFA : 0.70 : Norfenfluramine potency weight
K_HT5    : 1.20 : Serotonergic tone turnover (1/day)
EMAX_S1  : 0.10 : Maximum sigma-1 contribution
EC50_S1  : 0.055 : Sigma-1 EC50 (mg/L)
EMAX_CBD : 0.11082 : Maximum cannabidiol direct protection
EC50_CBD_PD : 0.05 : Cannabidiol EC50 (mg/L); saturated at 10 mg/kg/day
W_7OH    : 0.50 : 7-OH-CBD potency weight
K_ADEN   : 1.00 : Adenosine tone turnover (1/day)
EMAX_ADEN : 0.0 : Adenosine contribution (deliberately zero; see README)
EC50_ADEN : 0.30 : Adenosine EC50 (mg/L)
EMAX_KETO : 0.22 : Ketogenic diet contribution
KETO     : 0.0  : Ketogenic diet switch
// ---------------- seizure hazard ----------------------------------------
HAZ0     : 1.03669 : Hazard scale, set so baseline MCSF is 15/month
EI0      : 1.00 : E:I normalisation
NH       : 2.20 : Hazard exponent in E:I
LAM_MAX  : 8.0  : Ceiling on countable convulsive seizures per day
SE_ESC   : 6.0  : Escalation into status as the countable rate saturates
K_KIND   : 0.0011 : Kindling accumulation
K_KIND_OFF : 0.010 : Kindling decay
K_REMOD  : 1.2e-4 : Chronic interneuron attrition
F_SE     : 0.020 : Status epilepticus share of convulsive events
// ---------------- triggers ----------------------------------------------
K_INF    : 0.35 : Infection resolution (1/day)
K_TEMP   : 2.0  : Temperature equilibration (1/day)
TEMP_GAIN : 2.2 : degC per unit infection burden
TRIG_GAIN : 1.6 : Hazard multiplier per degC above 37
// ---------------- development -------------------------------------------
DQ0      : 100.0 : Starting developmental quotient
K_DQ_SZ  : 0.052 : DQ loss per unit seizure rate
K_DQ_SED : 0.028 : DQ loss per unit somnolence
DQ_FLOOR : 22.0 : Floor on DQ
AGE0     : 1.0  : Age at simulation start (years)
// ---------------- safety ------------------------------------------------
EMAX_SOMN : 0.62 : Maximum somnolence
PAM50_SOMN : 2.60 : PAM load for half-maximal somnolence
K_SOMN   : 0.35 : Somnolence equilibration (1/day)
K_WGT_FFA : 0.40 : Weight loss per mg/L fenfluramine
K_WGT_STP : 0.012 : Weight loss per 10 mg/L stiripentol
K_WGT_REC : 0.020 : Weight recovery (1/day)
ALT_BASE : 25.0 : Baseline ALT (U/L)
K_ALT    : 1.9  : ALT response to cannabidiol
K_ALT_VPA : 2.4 : Valproate amplification of the ALT signal
K_ALT_OUT : 0.045 : ALT turnover (1/day)
K_VLV    : 0.010 : Valve index accumulation per mg/L norfenfluramine
K_VLV_OFF : 0.020 : Valve index resolution
K_SUDEP  : 1.15e-4 : SUDEP hazard per unit seizure rate
W_SUDEP_SE : 6.0 : SUDEP weight on status epilepticus
K_FILT   : 0.25 : Observation filter (1/day)
// ---------------- decomposition switches --------------------------------
PK_ROUTE : 1.0 : 1 = CYP inhibition allowed
PD_ROUTE : 1.0 : 1 = direct target engagement allowed
NAIVE    : 0.0 : 1 = start benzodiazepine-naive
BG_CLB   : 0.5 : Background clobazam dose for the tolerance initialisation
BG_VPA   : 30.0 : Background valproate dose

$CMT @annotated
CLB_G   : Clobazam gut depot (mg/kg)
CLB_C   : Clobazam central (mg/kg)
CLB_P   : Clobazam peripheral (mg/kg)
NCLB_C  : Norclobazam central (mg/kg)
NCLB_P  : Norclobazam peripheral (mg/kg)
STP_G   : Stiripentol gut depot (mg/kg)
STP_C   : Stiripentol central (mg/kg)
CBD_G   : Cannabidiol gut depot (mg/kg)
CBD_C   : Cannabidiol central (mg/kg)
CBD_P   : Cannabidiol peripheral (mg/kg)
CBD_7OH : 7-OH-CBD (mg/kg)
CBD_7COOH : 7-COOH-CBD (mg/kg)
FFA_G   : Fenfluramine gut depot (mg/kg)
FFA_C   : Fenfluramine central (mg/kg)
FFA_P   : Fenfluramine peripheral (mg/kg)
FFA_NOR : Norfenfluramine (mg/kg)
VPA_G   : Valproate gut depot (mg/kg)
VPA_C   : Valproate central (mg/kg)
NVB_G   : Sodium-channel blocker gut depot (mg/kg)
NVB_C   : Sodium-channel blocker central (mg/kg)
DZP_C   : Rescue diazepam central (mg/kg)
E2C19   : CYP2C19 relative activity
E3A4    : CYP3A4 relative activity
E1A2    : CYP1A2/2B6/2D6 relative activity
EUGT    : UGT1A4 relative activity
ASO_CSF : Antisense oligonucleotide in CSF (mg)
ASO_BR  : Antisense oligonucleotide in brain (mg)
NAV_INT : Nav1.1 function in interneurons (normalised)
RGABA   : GABA-A receptor availability (tolerance state)
HT5     : Serotonergic tone (mg/L equivalent)
ADEN    : Adenosine tone (mg/L equivalent)
KIND    : Kindling / use-dependent excitability
REMOD   : Chronic interneuron attrition
TCORE   : Core temperature (degC)
INFECT  : Infection burden
BURD    : Cumulative convulsive seizures
SECNT   : Cumulative status epilepticus episodes
DQ      : Developmental quotient
SOMN    : Somnolence (0-1)
WGT     : Weight z-score change
SUDEPH  : Cumulative SUDEP hazard
ALT     : Alanine aminotransferase (U/L)
VLV     : Valve regurgitation index
MCSF_F  : Filtered monthly convulsive seizure frequency

$GLOBAL
#define DPM 30.4375

$MAIN
F_CLB_G = F_CLB;
F_STP_G = F_STP;
F_CBD_G = F_CBD;
F_FFA_G = F_FFA;
F_VPA_G = F_VPA;
F_NVB_G = F_NVB;

E2C19_0 = E2C19_GENO;
E3A4_0  = 1.0;
E1A2_0  = 1.0;
EUGT_0  = 1.0;
NAV_INT_0 = ALLELE;
TCORE_0 = 37.0;
DQ_0    = DQ0;
ALT_0   = ALT_BASE;

// Benzodiazepine tolerance is initialised at the steady state implied by the
// BACKGROUND clobazam dose, because patients entering an add-on trial have
// been taking clobazam for months. Starting at 1.0 instead would give every
// arm a spurious upward drift in seizure frequency across the trial.
double cclb0  = BG_CLB * F_CLB / CL_CLB;
double cnclb0 = FM_NCLB * BG_CLB * F_CLB / (CL_NCLB * (E2C19_GENO > 1e-6 ? E2C19_GENO : 1e-6));
double pam0   = cclb0 / EC50_PAM_CLB + cnclb0 / EC50_PAM_NCLB;
RGABA_0 = (NAIVE >= 0.5) ? 1.0 : (1.0 - FTOL * pam0 / (pam0 + PAM_TOL));

$ODE
// ---- concentrations -----------------------------------------------------
double c_clb  = CLB_C   / V_CLB;
double c_nclb = NCLB_C  / V_NCLB;
double c_stp  = STP_C   / V_STP;
double c_cbd  = CBD_C   / V_CBD;
double c_7oh  = CBD_7OH / V_7OH;
double c_ffa  = FFA_C   / V_FFA;
double c_nor  = FFA_NOR / V_NOR;
double c_vpa  = VPA_C   / V_VPA;
double c_nvb  = NVB_C   / V_NVB;
double c_dzp  = DZP_C   / V_DZP;

// ---- enzyme turnover: this is where the PK route of the DDI lives -------
double inh2c19 = PK_ROUTE * (c_stp / KI_STP_2C19
                           + c_cbd / KI_CBD_2C19
                           + c_7oh / KI_7OH_2C19);
double e2c19ss = E2C19_GENO / (1.0 + inh2c19);
double e3a4ss  = 1.0 / (1.0 + PK_ROUTE * (c_stp / KI_STP_3A4 + c_cbd / KI_CBD_3A4));
double e1a2ss  = 1.0 / (1.0 + PK_ROUTE * c_stp / KI_STP_1A2);
double eugtss  = 1.0 / (1.0 + c_vpa / KI_VPA_UGT + c_cbd / KI_CBD_UGT);
dxdt_E2C19 = KE_ENZ * (e2c19ss - E2C19);
dxdt_E3A4  = KE_ENZ * (e3a4ss  - E3A4);
dxdt_E1A2  = KE_ENZ * (e1a2ss  - E1A2);
dxdt_EUGT  = KE_ENZ * (eugtss  - EUGT);

// ---- PK ------------------------------------------------------------------
double ka_clb  = KA_CLB * CLB_G;
double cl_form = FM_NCLB * CL_CLB * E3A4;
double cl_oth  = (1.0 - FM_NCLB) * CL_CLB;
double q_clb   = Q_CLB * (c_clb - CLB_P / VP_CLB);
dxdt_CLB_G = -ka_clb;
dxdt_CLB_C = ka_clb - (cl_form + cl_oth) * c_clb - q_clb;
dxdt_CLB_P = q_clb;

double q_nclb = Q_NCLB * (c_nclb - NCLB_P / VP_NCLB);
dxdt_NCLB_C = cl_form * c_clb - CL_NCLB * E2C19 * c_nclb - q_nclb;
dxdt_NCLB_P = q_nclb;

double ka_stp = KA_STP * STP_G;
dxdt_STP_G = -ka_stp;
dxdt_STP_C = ka_stp - VMAX_STP * c_stp / (KM_STP + c_stp);

double ka_cbd = KA_CBD * CBD_G;
double q_cbd  = Q_CBD * (c_cbd - CBD_P / VP_CBD);
double cl_7ohf = FM_7OH * CL_CBD * E2C19;
double cl_cbdo = (1.0 - FM_7OH) * CL_CBD;
dxdt_CBD_G = -ka_cbd;
dxdt_CBD_C = ka_cbd - (cl_7ohf + cl_cbdo) * c_cbd - q_cbd;
dxdt_CBD_P = q_cbd;
dxdt_CBD_7OH = cl_7ohf * c_cbd - CL_7OH * c_7oh;
dxdt_CBD_7COOH = CL_7OH * c_7oh - CL_7COOH * CBD_7COOH / V_7COOH;

double ka_ffa = KA_FFA * FFA_G;
double cl_ffa = CL_FFA * (1.0 - FRAC_CYP_FFA + FRAC_CYP_FFA * E1A2);
double q_ffa  = Q_FFA * (c_ffa - FFA_P / VP_FFA);
dxdt_FFA_G = -ka_ffa;
dxdt_FFA_C = ka_ffa - cl_ffa * c_ffa - q_ffa;
dxdt_FFA_P = q_ffa;
dxdt_FFA_NOR = FM_NOR * cl_ffa * c_ffa - CL_NOR * c_nor;

double ka_vpa = KA_VPA * VPA_G;
dxdt_VPA_G = -ka_vpa;
dxdt_VPA_C = ka_vpa - CL_VPA * c_vpa;

double ka_nvb = KA_NVB * NVB_G;
double cl_nvb = CL_NVB * (1.0 - FRAC_UGT_NVB + FRAC_UGT_NVB * EUGT);
dxdt_NVB_G = -ka_nvb;
dxdt_NVB_C = ka_nvb - cl_nvb * c_nvb;

dxdt_DZP_C = -CL_DZP * c_dzp;

dxdt_ASO_CSF = -0.9 * ASO_CSF;
dxdt_ASO_BR  = 0.9 * 0.45 * ASO_CSF - 0.0116 * ASO_BR;

// ---- target biology ------------------------------------------------------
// Poison-exon skipping restores productive splicing, but only from the intact
// wild-type allele, so there is a hard arithmetic ceiling at ALLELE/(1-PE).
double aso_eff = ASO_EMAX * ASO_BR / (ASO_EC50 + ASO_BR);
double pe      = PE_FRAC * (1.0 - aso_eff);
double nav_tgt = (ALLELE * (1.0 - pe)) / (1.0 - PE_FRAC);
dxdt_NAV_INT = K_NAV * (nav_tgt - NAV_INT);

double sens = (ALLELE < 0.99) ? SENS_T : SENS_T_WT;
double fT   = 1.0 - sens * ((TCORE > 37.0) ? (TCORE - 37.0) : 0.0);
if(fT < 0.05) fT = 0.05;

// Use-dependent sodium-channel block: MORE occupancy where firing is faster.
double occ_int = c_nvb / (c_nvb + IC50_NVB_INT);
double occ_pyr = c_nvb / (c_nvb + IC50_NVB_PYR);
double na_int  = NAV_INT * (1.0 - occ_int) * fT;
double na_pyr  = 1.0 * (1.0 - occ_pyr);
if(na_int < 1e-9) na_int = 1e-9;
if(na_pyr < 1e-9) na_pyr = 1e-9;

double ni  = pow(na_int, HILL_INT);
double cap_int = ni / (ni + pow(EC50_INT, HILL_INT));
double np_ = pow(na_pyr, HILL_PYR);
double cap_pyr = np_ / (np_ + pow(EC50_PYR, HILL_PYR));

// Benzodiazepine site: saturates, and Hashi 2015 puts the saturation BELOW
// the concentration a standard clobazam dose already reaches.
double pam = c_clb / EC50_PAM_CLB + c_nclb / EC50_PAM_NCLB + c_dzp / EC50_PAM_DZP;
double gain_pam = 1.0 + EMAX_PAM * RGABA * pam / (1.0 + pam);
double rg_ss = 1.0 - FTOL * pam / (pam + PAM_TOL);
dxdt_RGABA = KTOL * (rg_ss - RGABA);

// Stiripentol binds a DISTINCT alpha3-preferring site, so it gets its own
// saturable term rather than competing for an already-saturated one.
double pam_stp  = PD_ROUTE * c_stp / EC50_PAM_STP;
double gain_stp = 1.0 + EMAX_STP_SITE * pam_stp / (1.0 + pam_stp);

double gain_vpa = 1.0 + EMAX_VPA * c_vpa / (EC50_VPA + c_vpa);

// Non-GABAergic protection.
double ffa_drive = c_ffa + W_NORFFA * c_nor;
dxdt_HT5 = K_HT5 * (ffa_drive - HT5);
double e_5ht = EMAX_5HT * HT5 / (EC50_5HT + HT5);
double e_s1  = EMAX_S1 * c_ffa / (EC50_S1 + c_ffa);
double cbd_drive = PD_ROUTE * (c_cbd + W_7OH * c_7oh);
double e_cbd = EMAX_CBD * cbd_drive / (EC50_CBD_PD + cbd_drive);
dxdt_ADEN = K_ADEN * (cbd_drive - ADEN);
double e_aden = EMAX_ADEN * ADEN / (EC50_ADEN + ADEN);
double prot_ind = 1.0 + e_5ht + e_s1 + e_cbd + e_aden + EMAX_KETO * KETO;

// ---- E:I balance and seizure hazard -------------------------------------
double inh = cap_int * gain_pam * gain_stp * gain_vpa * (1.0 - REMOD);
if(inh < 1e-9) inh = 1e-9;
double exc = cap_pyr / prot_ind;
double ei  = exc / inh;

double prop = 1.0 - EMAX_NVB_PROP * c_nvb / (c_nvb + IC50_NVB_PROP);
double trig = 1.0 + TRIG_GAIN * ((TCORE > 37.0) ? (TCORE - 37.0) : 0.0);

double raw = HAZ0 * pow(ei / EI0, NH) * (1.0 + KIND) * trig * prop;
// The countable convulsive seizure rate has a physiological ceiling; past it
// the discrete events merge and the correct description is status epilepticus.
double lam = LAM_MAX * raw / (raw + LAM_MAX);
double sat = lam / LAM_MAX;
double lam_se = F_SE * lam * trig * (1.0 + SE_ESC * sat);

dxdt_KIND  = K_KIND * lam - K_KIND_OFF * KIND;
dxdt_REMOD = K_REMOD * lam * (1.0 - REMOD);

dxdt_INFECT = -K_INF * INFECT;
dxdt_TCORE  = K_TEMP * (37.0 + TEMP_GAIN * INFECT - TCORE);

dxdt_BURD  = lam;
dxdt_SECNT = lam_se;

double som_ss = EMAX_SOMN * pam / (PAM50_SOMN + pam);
dxdt_SOMN = K_SOMN * (som_ss - SOMN);

double age  = AGE0 + SOLVERTIME / 365.25;
double agew = exp(-0.22 * ((age > 1.0) ? (age - 1.0) : 0.0));
dxdt_DQ = (DQ <= DQ_FLOOR) ? 0.0
          : -(K_DQ_SZ * lam + K_DQ_SED * SOMN) * agew;

dxdt_WGT = -(K_WGT_FFA * c_ffa + K_WGT_STP * c_stp / 10.0) - K_WGT_REC * WGT;
dxdt_SUDEPH = K_SUDEP * (lam + W_SUDEP_SE * lam_se);

double alt_in = ALT_BASE * K_ALT_OUT
              * (1.0 + K_ALT * (c_cbd + 0.4 * c_7oh)
                     * (1.0 + K_ALT_VPA * c_vpa / 70.0));
dxdt_ALT = alt_in - K_ALT_OUT * ALT;
dxdt_VLV = K_VLV * c_nor - K_VLV_OFF * VLV;

dxdt_MCSF_F = K_FILT * (lam * DPM - MCSF_F);

$TABLE
double C_CLB  = CLB_C / V_CLB;
double C_NCLB = NCLB_C / V_NCLB;
double C_STP  = STP_C / V_STP;
double C_CBD  = CBD_C / V_CBD;
double C_7OH  = CBD_7OH / V_7OH;
double C_FFA  = FFA_C / V_FFA;
double C_NOR  = FFA_NOR / V_NOR;
double C_VPA  = VPA_C / V_VPA;
double C_NVB  = NVB_C / V_NVB;
double RATIO_NCLB = (C_CLB > 1e-12) ? C_NCLB / C_CLB : 0.0;

double occI = C_NVB / (C_NVB + IC50_NVB_INT);
double occP = C_NVB / (C_NVB + IC50_NVB_PYR);
double fTt  = 1.0 - ((ALLELE < 0.99) ? SENS_T : SENS_T_WT)
                    * ((TCORE > 37.0) ? (TCORE - 37.0) : 0.0);
if(fTt < 0.05) fTt = 0.05;
double naI = NAV_INT * (1.0 - occI) * fTt;
double naP = 1.0 * (1.0 - occP);
if(naI < 1e-9) naI = 1e-9;
double CAP_INT = pow(naI, HILL_INT) / (pow(naI, HILL_INT) + pow(EC50_INT, HILL_INT));
double CAP_PYR = pow(naP, HILL_PYR) / (pow(naP, HILL_PYR) + pow(EC50_PYR, HILL_PYR));

double PAM = C_CLB / EC50_PAM_CLB + C_NCLB / EC50_PAM_NCLB + (DZP_C / V_DZP) / EC50_PAM_DZP;
double GAIN_PAM = 1.0 + EMAX_PAM * RGABA * PAM / (1.0 + PAM);
double PSTP = PD_ROUTE * C_STP / EC50_PAM_STP;
double GAIN_STP = 1.0 + EMAX_STP_SITE * PSTP / (1.0 + PSTP);
double GAIN_VPA = 1.0 + EMAX_VPA * C_VPA / (EC50_VPA + C_VPA);
double PROT = 1.0 + EMAX_5HT * HT5 / (EC50_5HT + HT5)
                  + EMAX_S1 * C_FFA / (EC50_S1 + C_FFA)
                  + EMAX_CBD * (PD_ROUTE * (C_CBD + W_7OH * C_7OH))
                    / (EC50_CBD_PD + PD_ROUTE * (C_CBD + W_7OH * C_7OH))
                  + EMAX_KETO * KETO;
double INH = CAP_INT * GAIN_PAM * GAIN_STP * GAIN_VPA * (1.0 - REMOD);
if(INH < 1e-9) INH = 1e-9;
double EI = (CAP_PYR / PROT) / INH;
double PROPB = 1.0 - EMAX_NVB_PROP * C_NVB / (C_NVB + IC50_NVB_PROP);
double TRIGF = 1.0 + TRIG_GAIN * ((TCORE > 37.0) ? (TCORE - 37.0) : 0.0);
double RAW = HAZ0 * pow(EI / EI0, NH) * (1.0 + KIND) * TRIGF * PROPB;
double LAM = LAM_MAX * RAW / (RAW + LAM_MAX);
double MCSF = LAM * DPM;
double SE_YR = F_SE * LAM * TRIGF * (1.0 + SE_ESC * LAM / LAM_MAX) * 365.25;
double SUDEP_CUM = 100.0 * (1.0 - exp(-SUDEPH));

$CAPTURE C_CLB C_NCLB RATIO_NCLB C_STP C_CBD C_7OH C_FFA C_NOR C_VPA C_NVB
$CAPTURE CAP_INT CAP_PYR PAM GAIN_PAM GAIN_STP INH EI MCSF SE_YR SUDEP_CUM
'

mod <- mcode("dravet", dravet_code)

# =============================================================================
#  DOSING HELPERS
# =============================================================================
#  All doses are mg/kg per day; `nd` is the number of administrations per day.
#  ADDON_START is day 42, matching the 4-week baseline observation used for the
#  trial-style read-out below.
# =============================================================================
BASE_START <- 14; BASE_END <- 42
TRT_START  <- 56; TRT_END  <- 140
ADDON_START <- 42; TITR <- 7

dose_seq <- function(cmt, mgkg_day, start = 0, until = 200, nd = 2) {
  if (mgkg_day <= 0 || until <= start) return(NULL)
  n <- max(1, floor((until - start) * nd))
  data.frame(ID = 1, time = start + seq_len(n) / nd - 1 / nd,
             amt = mgkg_day / nd, cmt = cmt, evid = 1)
}

bind_ev <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(data.frame(ID = 1, time = 0, amt = 0,
                                        cmt = "CLB_G", evid = 1)[0, ])
  d <- do.call(rbind, parts)
  d[order(d$time), ]
}

background_reg <- function(clb = 0.5, vpa = 30, until = 200) {
  bind_ev(dose_seq("CLB_G", clb, 0, until, 2),
          dose_seq("VPA_G", vpa, 0, until, 2))
}

titrated <- function(cmt, dose, nd = 2, until = 200) {
  bind_ev(dose_seq(cmt, dose * 0.5, ADDON_START, ADDON_START + TITR, nd),
          dose_seq(cmt, dose, ADDON_START + TITR, until, nd))
}

addon <- function(what, dose, until = 200) {
  switch(what,
    none = NULL,
    ffa  = titrated("FFA_G", dose, 2, until),
    cbd  = titrated("CBD_G", dose, 2, until),
    stp  = titrated("STP_G", dose, 3, until),
    ltg  = dose_seq("NVB_G", dose, ADDON_START, until, 2),
    aso  = data.frame(ID = 1, time = ADDON_START + c(0, 90, 180),
                      amt = dose, cmt = "ASO_CSF", evid = 1),
    stop("unknown add-on"))
}

# -----------------------------------------------------------------------------
#  Trial-style read-out: median monthly convulsive seizure frequency over a
#  4-week baseline observation window, then over the maintenance window after
#  titration -- the structure the registration trials actually used.
# -----------------------------------------------------------------------------
run_arm <- function(mod, arm = "none", dose = 0, clb = 0.5, vpa = 30,
                    tend = 140, extra = NULL, param = list()) {
  reg <- bind_ev(background_reg(clb, vpa, tend), addon(arm, dose, tend), extra)
  m <- param(mod, c(list(BG_CLB = clb, BG_VPA = vpa), param))
  # delta must be fine relative to the dosing interval: BID dosing with a
  # 2 h absorption half-life and an 11 h parent half-life means a 6 h output
  # grid ALIASES the peaks and understates every average concentration.
  out <- mrgsim_d(m, reg, end = tend, delta = 0.02, atol = 1e-9, rtol = 1e-9,
                  hmax = 0.1)
  as.data.frame(out)
}

pct_reduction <- function(mod, arm = "none", dose = 0, clb = 0.5, vpa = 30,
                          tend = 140, ...) {
  d <- run_arm(mod, arm, dose, clb, vpa, tend, ...)
  cum <- function(tt) approx(d$time, d$BURD, xout = tt, ties = "ordered")$y
  base <- (cum(BASE_END) - cum(BASE_START)) / (BASE_END - BASE_START) * 30.4375
  trt  <- (cum(TRT_END) - cum(TRT_START)) / (TRT_END - TRT_START) * 30.4375
  c(baseline = base, treated = trt, reduction = 100 * (base - trt) / base)
}

# -----------------------------------------------------------------------------
#  SCENARIOS (21)
# -----------------------------------------------------------------------------
scenarios <- list(
  # --- background and monotherapy
  list(id = 1,  label = "Valproate 30 alone",                arm = "none", dose = 0,   clb = 0.0),
  list(id = 2,  label = "Valproate + clobazam 0.25",         arm = "none", dose = 0,   clb = 0.25),
  list(id = 3,  label = "Valproate + clobazam 0.5 (reference)", arm = "none", dose = 0, clb = 0.5),
  list(id = 4,  label = "Valproate + clobazam 1.0",          arm = "none", dose = 0,   clb = 1.0),
  # --- the three approved add-ons
  list(id = 5,  label = "+ stiripentol 50 mg/kg/day",        arm = "stp",  dose = 50),
  list(id = 6,  label = "+ cannabidiol 10 mg/kg/day",        arm = "cbd",  dose = 10),
  list(id = 7,  label = "+ cannabidiol 20 mg/kg/day",        arm = "cbd",  dose = 20),
  list(id = 8,  label = "+ fenfluramine 0.2 mg/kg/day",      arm = "ffa",  dose = 0.2),
  list(id = 9,  label = "+ fenfluramine 0.7 mg/kg/day",      arm = "ffa",  dose = 0.7),
  # --- combination and the dose cap
  list(id = 10, label = "+ stiripentol + fenfluramine 0.4",  arm = "stp",  dose = 50,
       extra = quote(titrated("FFA_G", 0.4, 2, 140))),
  # --- decomposition: structural nulls
  list(id = 11, label = "stiripentol, PD route deleted",     arm = "stp",  dose = 50,
       param = list(PD_ROUTE = 0)),
  list(id = 12, label = "stiripentol, PK route deleted",     arm = "stp",  dose = 50,
       param = list(PK_ROUTE = 0)),
  list(id = 13, label = "cannabidiol 20, PD route deleted",  arm = "cbd",  dose = 20,
       param = list(PD_ROUTE = 0)),
  list(id = 14, label = "cannabidiol 20, PK route deleted",  arm = "cbd",  dose = 20,
       param = list(PK_ROUTE = 0)),
  # --- pharmacogenetics
  list(id = 15, label = "stiripentol, CYP2C19 poor metaboliser", arm = "stp", dose = 50,
       param = list(E2C19_GENO = 0.20)),
  list(id = 16, label = "cannabidiol 20, CYP2C19 poor metaboliser", arm = "cbd", dose = 20,
       param = list(E2C19_GENO = 0.20)),
  # --- the contraindicated arm, in both hosts
  list(id = 17, label = "lamotrigine 5 mg/kg/day, DRAVET",   arm = "ltg",  dose = 5),
  list(id = 18, label = "lamotrigine 5 mg/kg/day, HEALTHY CONTROL", arm = "ltg", dose = 5,
       param = list(ALLELE = 1.0)),
  # --- disease modification
  list(id = 19, label = "zorevunersen 60 mg intrathecal q90d", arm = "aso", dose = 60,
       tend = 400),
  # --- non-pharmacological and dose sparing
  list(id = 20, label = "ketogenic diet",                    arm = "none", dose = 0,
       param = list(KETO = 1.0)),
  list(id = 21, label = "cannabidiol 20 + clobazam cut to 0.25", arm = "cbd", dose = 20,
       clb = 0.25)
)

run_scenarios <- function(mod) {
  rows <- lapply(scenarios, function(s) {
    clb  <- if (is.null(s$clb)) 0.5 else s$clb
    tend <- if (is.null(s$tend)) 140 else s$tend
    par  <- if (is.null(s$param)) list() else s$param
    ex   <- if (is.null(s$extra)) NULL else eval(s$extra)
    r <- pct_reduction(mod, s$arm, s$dose, clb = clb, tend = tend,
                       extra = ex, param = par)
    data.frame(id = s$id, scenario = s$label,
               baseline_MCSF = round(r[["baseline"]], 2),
               treated_MCSF  = round(r[["treated"]], 2),
               pct_reduction = round(r[["reduction"]], 1))
  })
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
#  BETWEEN-PATIENT VARIABILITY (mirrors ETA_SPEC in the Python twin)
#  Applied by resampling parameters per subject rather than through $OMEGA,
#  so that the same log-normal spread is used in both implementations.
# -----------------------------------------------------------------------------
virtual_population <- function(mod, arm = "none", dose = 0, n = 150,
                               clb = 0.5, seed = 20260729) {
  set.seed(seed)
  cv <- c(CL_CLB = 0.30, CL_NCLB = 0.40, VMAX_STP = 0.35, CL_CBD = 0.45,
          CL_FFA = 0.30, HAZ0 = 0.75, EMAX_5HT = 0.45, EMAX_CBD = 0.35,
          EMAX_PAM = 0.20, EMAX_STP_SITE = 0.55, EC50_INT = 0.07)
  base <- as.list(param(mod))
  geno <- sample(c(1.0, 0.55, 0.20), n, replace = TRUE, prob = c(0.62, 0.30, 0.08))
  reds <- numeric(n)
  for (i in seq_len(n)) {
    pl <- lapply(names(cv), function(k) base[[k]] * exp(rnorm(1, 0, cv[[k]])))
    names(pl) <- names(cv)
    pl$E2C19_GENO <- geno[i]
    reds[i] <- pct_reduction(mod, arm, dose, clb = clb, param = pl)[["reduction"]]
  }
  list(median = median(reds),
       r50 = 100 * mean(reds >= 50), r75 = 100 * mean(reds >= 75),
       worse = 100 * mean(reds < 0), reds = reds)
}

# -----------------------------------------------------------------------------
#  DIAGNOSTICS — each asserts a published anchor or a structural claim
# -----------------------------------------------------------------------------
diagnostics <- function(mod) {
  chk <- function(name, value, target, tol, note = "") {
    ok <- abs(value - target) <= tol
    cat(sprintf("%-46s %10.3f  target %8.3f +/-%6.3f  %s %s\n",
                name, value, target, tol, if (ok) "PASS" else "FAIL", note))
    ok
  }
  d0 <- run_arm(mod, "none", 0)
  tail0 <- subset(d0, time >= 126)
  res <- c(
    chk("baseline MCSF (/month)",
        pct_reduction(mod, "none", 0)[["baseline"]], 15.0, 0.6),
    chk("clobazam Css (mg/L)", mean(tail0$C_CLB), 0.200, 0.02),
    chk("norclobazam Css (mg/L)", mean(tail0$C_NCLB), 1.600, 0.15),
    chk("N-CLB : CLB ratio", mean(tail0$RATIO_NCLB), 8.0, 1.0),
    chk("valproate Css (mg/L)", mean(tail0$C_VPA), 69.7, 4.0),
    chk("stiripentol Css (mg/L)",
        mean(subset(run_arm(mod, "stp", 50), time >= 126)$C_STP), 10.0, 1.5,
        "target 8-12"),
    chk("fenfluramine Css (ng/mL)",
        1000 * mean(subset(run_arm(mod, "ffa", 0.7), time >= 126)$C_FFA),
        55.0, 8.0),
    chk("norclobazam x on cannabidiol 20",
        mean(subset(run_arm(mod, "cbd", 20), time >= 126)$C_NCLB) /
          mean(tail0$C_NCLB), 6.0, 0.8, "Geffrey 2015"),
    chk("clobazam x on cannabidiol 20",
        mean(subset(run_arm(mod, "cbd", 20), time >= 126)$C_CLB) /
          mean(tail0$C_CLB), 1.60, 0.25, "Geffrey 2015"),
    chk("norclobazam x on stiripentol",
        mean(subset(run_arm(mod, "stp", 50), time >= 126)$C_NCLB) /
          mean(tail0$C_NCLB), 2.50, 0.4, "Jullien 2015"),
    chk("stiripentol reduction (%)",
        pct_reduction(mod, "stp", 50)[["reduction"]], 76.0, 4.0, "STICLO"),
    chk("fenfluramine 0.7 reduction (%)",
        pct_reduction(mod, "ffa", 0.7)[["reduction"]], 55.7, 4.0, "Study 1"),
    chk("cannabidiol 20 reduction (%)",
        pct_reduction(mod, "cbd", 20)[["reduction"]], 22.2, 4.0, "GWPCARE1/2"),
    chk("fenfluramine 0.2 reduction (%) [PREDICTION]",
        pct_reduction(mod, "ffa", 0.2)[["reduction"]], 23.1, 5.0, "Study 1"),
    chk("cannabidiol 10 reduction (%) [PREDICTION]",
        pct_reduction(mod, "cbd", 10)[["reduction"]], 21.8, 5.0, "GWPCARE2"),
    chk("Dravet interneuron capacity",
        mean(tail0$CAP_INT), 0.653, 0.02, "structural"),
    chk("healthy interneuron capacity",
        mean(subset(run_arm(mod, "none", 0, param = list(ALLELE = 1.0)),
                    time >= 126)$CAP_INT), 0.992, 0.02, "structural")
  )
  cat(sprintf("\n%d of %d diagnostics passed\n", sum(res), length(res)))
  invisible(res)
}

# -----------------------------------------------------------------------------
#  THE CENTRAL DECOMPOSITION
# -----------------------------------------------------------------------------
decompose <- function(mod) {
  f <- function(arm, dose) {
    both <- pct_reduction(mod, arm, dose)[["reduction"]]
    pk   <- pct_reduction(mod, arm, dose, param = list(PD_ROUTE = 0))[["reduction"]]
    pd   <- pct_reduction(mod, arm, dose, param = list(PK_ROUTE = 0))[["reduction"]]
    c(both = both, pk_only = pk, pd_only = pd, pk_share = 100 * pk / both)
  }
  rbind(`stiripentol 50` = f("stp", 50), `cannabidiol 20` = f("cbd", 20))
}

# -----------------------------------------------------------------------------
#  THE SODIUM-CHANNEL-BLOCKER SIGN FLIP
#  Identical drug, identical dose, identical equations. Only ALLELE differs.
# -----------------------------------------------------------------------------
navblock_signflip <- function(mod, dose = 5) {
  out <- lapply(c(Dravet = 0.5, healthy = 1.0), function(a) {
    r <- pct_reduction(mod, "ltg", dose, param = list(ALLELE = a))
    c(baseline = r[["baseline"]], on_drug = r[["treated"]],
      pct_change = -r[["reduction"]])
  })
  do.call(rbind, out)
}

if (identical(environment(), globalenv()) && !interactive()) {
  cat("\n=== DIAGNOSTICS ===\n");        diagnostics(mod)
  cat("\n=== SCENARIOS ===\n");          print(run_scenarios(mod), row.names = FALSE)
  cat("\n=== DECOMPOSITION ===\n");      print(round(decompose(mod), 1))
  cat("\n=== SODIUM-CHANNEL BLOCKER ===\n"); print(round(navblock_signflip(mod), 2))
}
