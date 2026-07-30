$PROB
 =============================================================================
 Gastrointestinal Stromal Tumour (GIST) — QSP / PK-PD model for mrgsolve
 =============================================================================
 49 compartments · 5 drugs with active metabolites · clone-resolved KIT/PDGFRA
 signalling · quiescent reservoir · imaged-vs-viable mass · organ toxicity ·
 26 therapeutic scenarios.

 Companion files
   gist_qsp_model.dot / .svg / .png   mechanistic map
   gist_python_twin.py                dependency-free executable twin; it is
                                      the numerical reference for every value
                                      quoted in README.md and it self-checks
   gist_scenarios.R                   R driver: the 26 scenarios and endpoints
   gist_references.md                 the literature the parameters come from

 ---------------------------------------------------------------------------
 WHY THIS MODEL IS BUILT THE WAY IT IS
 ---------------------------------------------------------------------------
 The usual drawing of GIST is a chain: KIT mutation -> constitutive kinase ->
 proliferation -> imatinib occupies the ATP pocket -> tumour shrinks -> a
 resistance mutation appears -> switch drug.  One tumour, one kinase activity,
 one IC50, one resistance state.  Six things in the trial record cannot all be
 true at once in that drawing:

  O1 Sunitinib and ripretinib tie in second line (INTRIGUE ITT 8.3 vs 8.0 mo),
     but split by ctDNA genotype: KIT exon 11+13/14 -> sunitinib 15.0 vs
     ripretinib 4.0 mo; exon 11+17/18 -> ripretinib 14.2 vs sunitinib 1.5 mo.
  O2 Imatinib RECHALLENGE after failure of imatinib, sunitinib and everything
     else still beats placebo (RIGHT, 1.8 vs 0.9 mo).
  O3 800 mg helps only KIT exon 9 (MetaGIST; EORTC 62005 risk of progression
     reduced 61%).
  O4 PDGFRA D842V is the most imatinib-resistant and the most
     avapritinib-sensitive genotype (NAVIGATOR ORR ~88%).
  O5 FDG-PET SUV collapses in 24-48 h while CT size barely moves for weeks.
  O6 Interrupting imatinib in responders gives PFS 6.1 / 7.0 / 12.0 months
     (after 1 / 3 / 5 years) versus 27.8 / 67.0 / not reached on continuation
     (BFR14) — years of deep response do not eradicate the disease, and the
     time to regrowth gets LONGER the longer the drug was given.

 Three structural positions are taken instead:

  C1 THE TUMOUR IS A POPULATION, NOT A SIZE.  Four clones (primary genotype,
     exon 13/14 ATP-binding-pocket, exon 17/18 activation-loop,
     KIT-independent bypass), each cycling and quiescent, each with its own
     EC50 for each of five drugs.  Subclones are generated per cell DIVISION,
     overwhelmingly during the pre-diagnostic expansion, so they are
     PRE-EXISTING rather than induced.  Efficacy is a SET-COVER problem:
     progression is the growth of the least-covered clone.  Two drugs can
     therefore swap rank between genotypes while their means are equal (O1),
     and a majority clone carrying only the primary mutation survives to the
     end of the sequence, which is why rechallenge works (O2).

  C2 OCCUPANCY IS FAST, KILLING IS SLOW, AND THE DRUG IS MOSTLY CYTOSTATIC.
     Occupancy (hours, read by PET) -> cell-cycle exit (days, read by Ki67)
     -> a quiescent reservoir the drug cannot kill (months-years, read by what
     happens when you stop).  Wake-up from quiescence needs mitogenic signal
     steeply (SIG^3) and is slowed further by dormancy depth.  The imaged mass
     is a third variable: a synchronous drug-induced die-off leaves
     slowly-resorbed non-viable tissue, so the residual mass on CT is mostly
     not viable tumour.  O5 and O6 are then the same fact.

  C3 EXPOSURE MATTERS WHERE THE GENOTYPE PUTS THE EC50.  Imatinib PK is
     explicit (CYP3A4 autoinduction, AGP binding, CGP74588).  At 400 mg the
     achieved concentration sits ~8-fold above the exon 11 EC50 and ~1.5-fold
     above the exon 9 EC50, i.e. on opposite sides of the proliferation
     threshold, so escalation moves exon 9 across it and does nothing for
     exon 11 (O3).

 FITTED PARAMETERS: KPMAX, KD0 (untreated doubling time), TURNOVER
 (pre-diagnostic divisions per surviving cell -> pre-existing resistant pool,
 fitted to first-line PFS), FNEC/KRES (depth of RECIST response), KDQ and
 PHI_D (reservoir attrition and dormancy depth, fitted to BFR14).  Eight of
 ~150.  Drug potency RATIOS across clones are taken from published
 biochemical/cellular data and were never fitted to a clinical endpoint —
 which is what makes the INTRIGUE crossover a prediction.

 Concentrations are TOTAL PLASMA in ng/mL and EC50 values are referenced to
 total plasma concentration, so plasma protein binding and tissue penetration
 are absorbed into the constant.  The exception is imatinib, where AGP
 binding is explicit: EC50_eff = EC50 * (AGP/AGP0)^HAGP.

 Usage
   library(mrgsolve)
   mod <- mread("gist_mrgsolve_model.R")
   out <- mrgsim(mod, ev(amt = 400, ii = 1, addl = 900, cmt = "A_IM"),
                 end = 900, delta = 1)
   plot(out, SLD + SUV + N_TOT + VAF_R ~ time)
   source("gist_scenarios.R")      # all 26 scenarios and their endpoints
 =============================================================================

# Gastrointestinal stromal tumour: clone-resolved QSP model
# 49 compartments, 5 drugs, 26 scenarios.  See header for the argument.

// -----------------------------------------------------------------------------
$PARAM @annotated
// ---- genotype switch --------------------------------------------------------
GENO   : 1    : Primary genotype (1 KIT exon 11, 2 KIT exon 9, 3 PDGFRA D842V, 4 SDH-deficient/WT)
SETTING: 1    : 1 = metastatic, 2 = adjuvant (occult residual disease)
V0     : 400  : Baseline viable tumour volume, metastatic setting (mL)

// ---- imatinib PK (400 mg qd; CL/F ~ 10 L/h pre-induction, V/F 520 L) --------
F_IM   : 0.98 : Imatinib bioavailability (-)
KA_IM  : 12.0 : Imatinib absorption rate (1/day)
V1_IM  : 520  : Imatinib central volume (L)
VP_IM  : 420  : Imatinib peripheral volume (L)
Q_IM   : 220  : Imatinib intercompartmental clearance (L/day)
CL_IM  : 232  : Imatinib clearance before autoinduction (L/day)
FM_IM  : 0.12 : Fraction of imatinib clearance forming CGP74588 (-)
CLM_IM : 250  : CGP74588 clearance (L/day)
VM_IM  : 150  : CGP74588 volume (L)
KENZ   : 0.15 : CYP3A4 pool turnover (1/day)
EMAXENZ: 0.50 : Maximal CYP3A4 autoinduction by imatinib (-)
EC50ENZ: 1000 : Imatinib concentration for half-maximal autoinduction (ng/mL)
INDF   : 1.0  : External CYP3A4 induction factor (rifampicin ~2.6, azole < 1)
CLF    : 1.0  : Clearance multiplier for the individual (PK variability)

// ---- sunitinib PK (50 mg qd 4/2; t1/2 ~ 50 h) -------------------------------
F_SU   : 1.0  : Sunitinib bioavailability (-)
KA_SU  : 6.0  : Sunitinib absorption rate (1/day)
V1_SU  : 2000 : Sunitinib central volume (L)
VP_SU  : 1000 : Sunitinib peripheral volume (L)
Q_SU   : 300  : Sunitinib intercompartmental clearance (L/day)
CL_SU  : 960  : Sunitinib clearance (L/day)
FM_SU  : 0.25 : Fraction forming SU12662 (-)
CLM_SU : 700  : SU12662 clearance (L/day)
VM_SU  : 1500 : SU12662 volume (L)

// ---- regorafenib PK (160 mg qd 3/1) ----------------------------------------
F_RE   : 1.0  : Regorafenib bioavailability (-)
KA_RE  : 4.0  : Regorafenib absorption rate (1/day)
V1_RE  : 150  : Regorafenib central volume (L)
VP_RE  : 120  : Regorafenib peripheral volume (L)
Q_RE   : 60   : Regorafenib intercompartmental clearance (L/day)
CL_RE  : 66   : Regorafenib clearance (L/day)
FM_RE  : 0.40 : Fraction forming M-2 + M-5 (-)
CLM_RE : 55   : M-2/M-5 clearance (L/day)
VM_RE  : 140  : M-2/M-5 volume (L)

// ---- ripretinib PK (150 mg qd) ---------------------------------------------
F_RI   : 1.0  : Ripretinib bioavailability (-)
KA_RI  : 6.0  : Ripretinib absorption rate (1/day)
V1_RI  : 400  : Ripretinib central volume (L)
VP_RI  : 300  : Ripretinib peripheral volume (L)
Q_RI   : 150  : Ripretinib intercompartmental clearance (L/day)
CL_RI  : 300  : Ripretinib clearance (L/day)
FM_RI  : 0.50 : Fraction forming DP-5439 (-)
CLM_RI : 250  : DP-5439 clearance (L/day)
VM_RI  : 350  : DP-5439 volume (L)

// ---- avapritinib PK (300 mg qd) --------------------------------------------
F_AV   : 1.0  : Avapritinib bioavailability (-)
KA_AV  : 6.0  : Avapritinib absorption rate (1/day)
V1_AV  : 1200 : Avapritinib central volume (L)
VP_AV  : 800  : Avapritinib peripheral volume (L)
Q_AV   : 250  : Avapritinib intercompartmental clearance (L/day)
CL_AV  : 450  : Avapritinib clearance (L/day)

// ---- clone x drug potency matrix (EC50, ng/mL total plasma) -----------------
// Primary clone, selected by GENO in $MAIN.
E11_IM : 185  : Imatinib EC50, KIT exon 11 (ng/mL)
E11_SU : 20   : Sunitinib EC50, KIT exon 11 (ng/mL)
E11_RE : 1400 : Regorafenib EC50, KIT exon 11 (ng/mL)
E11_RI : 210  : Ripretinib EC50, KIT exon 11 (ng/mL)
E11_AV : 170  : Avapritinib EC50, KIT exon 11 (ng/mL)
E9_IM  : 1000 : Imatinib EC50, KIT exon 9 (ng/mL)
E9_SU  : 15   : Sunitinib EC50, KIT exon 9 (ng/mL)
E9_RE  : 1200 : Regorafenib EC50, KIT exon 9 (ng/mL)
E9_RI  : 230  : Ripretinib EC50, KIT exon 9 (ng/mL)
E9_AV  : 180  : Avapritinib EC50, KIT exon 9 (ng/mL)
D842_IM: 20000: Imatinib EC50, PDGFRA D842V (ng/mL)
D842_SU: 2800 : Sunitinib EC50, PDGFRA D842V (ng/mL)
D842_RE: 30000: Regorafenib EC50, PDGFRA D842V (ng/mL)
D842_RI: 600  : Ripretinib EC50, PDGFRA D842V (ng/mL)
D842_AV: 40   : Avapritinib EC50, PDGFRA D842V (ng/mL)
// Secondary-mutation subclones.  The PATTERN is the published one: the
// ATP-binding-pocket mutants V654A/T670I keep sunitinib sensitivity and lose
// ripretinib/regorafenib potency; the activation-loop mutants D816/D820/N822K
// are the mirror image; avapritinib follows the activation loop.
C1_IM  : 4500 : Imatinib EC50, exon 13/14 ATP-pocket subclone (ng/mL)
C1_SU  : 12   : Sunitinib EC50, exon 13/14 subclone (ng/mL)
C1_RE  : 9000 : Regorafenib EC50, exon 13/14 subclone (ng/mL)
C1_RI  : 2200 : Ripretinib EC50, exon 13/14 subclone (ng/mL)
C1_AV  : 2700 : Avapritinib EC50, exon 13/14 subclone (ng/mL)
C2_IM  : 6500 : Imatinib EC50, exon 17/18 activation-loop subclone (ng/mL)
C2_SU  : 380  : Sunitinib EC50, exon 17/18 subclone (ng/mL)
C2_RE  : 950  : Regorafenib EC50, exon 17/18 subclone (ng/mL)
C2_RI  : 165  : Ripretinib EC50, exon 17/18 subclone (ng/mL)
C2_AV  : 130  : Avapritinib EC50, exon 17/18 subclone (ng/mL)
C3_EC  : 1e7  : EC50 of the KIT-independent bypass subclone, all drugs (ng/mL)
KITDEP1: 1.0  : KIT dependence of the ATP-pocket subclone (-)
KITDEP2: 1.0  : KIT dependence of the activation-loop subclone (-)
KITDEP3: 0.05 : KIT dependence of the bypass subclone (-)
KPS1   : 0.95 : Proliferation scale, ATP-pocket subclone (-)
KPS2   : 0.95 : Proliferation scale, activation-loop subclone (-)
KPS3   : 0.85 : Proliferation scale, bypass subclone (-)

// ---- tumour cell kinetics ---------------------------------------------------
KPMAX  : 0.0290 : FITTED maximal proliferation rate of cycling cells (1/day)
KD0    : 0.0175 : FITTED baseline death rate of cycling cells (1/day)
SIG50  : 0.30   : KIT signal giving half-maximal proliferation (-)
HP     : 4.0    : Hill slope, signal to proliferation (-)
KAMAX  : 0.0130 : Maximal drug-induced apoptosis (1/day) -- capped: cytostatic first
UCRIT  : 0.85   : (1 - signal) giving half-maximal apoptosis (-)
HA     : 6.0    : Hill slope, suppression to apoptosis (-)
KQIN   : 0.0040 : Cycling to quiescent, proportional to (1 - signal) (1/day)
KQOUT  : 0.150  : Quiescent to cycling, proportional to signal^NQ (1/day)
NQ     : 3.0    : Steepness of the wake-up requirement (-)
KDQ    : 0.00065: FITTED death rate of quiescent, drug-insensitive cells (1/day)
KDEEP  : 0.0025 : Dormancy deepens while the signal is off (1/day)
KSHAL  : 0.050  : Dormancy reverses when the signal returns (1/day)
PHI_D  : 0.90   : FITTED maximal slowing of wake-up by deep dormancy (-)
MU_ATP : 2.0e-7 : Exon 13/14 subclone per cell division (-)
MU_AL  : 3.0e-7 : Exon 17/18 subclone per cell division (-)
MU_BYP : 4.0e-9 : KIT-independent subclone per cell division (-)
TURNOVER: 85    : FITTED cumulative divisions per surviving cell pre-diagnosis (-)

// ---- mass and imaging -------------------------------------------------------
CELLSML: 1e9    : Cells per mL of viable tumour (-)
FNEC   : 0.62   : FITTED fraction of drug-induced dead-cell volume left non-viable (-)
FNEC0  : 0.15   : Same fraction for physiological turnover, efficiently cleared (-)
KRES   : 0.0030 : FITTED resorption of non-viable tissue (1/day)
KSLD   : 16.28  : SLD (mm) per cube root of total volume (mm/mL^(1/3))
TAUSLD : 10.0   : Imaging and remodelling lag (day)
SUV_BG : 1.4    : Background FDG SUV (-)
SUV0   : 9.0    : Untreated tumour SUVmax (-)
TAUPET : 0.8    : PET response lag (day)
KI67MX : 18.0   : Untreated Ki67 index (%)
TAUKI  : 3.0    : Ki67 lag (day)

// ---- vasculature and delivery ----------------------------------------------
KVEGFP : 0.30 : VEGF-A production per unit relative burden (1/day)
KVEGFD : 0.30 : VEGF-A degradation (1/day)
VREF   : 400  : Reference tumour volume for VEGF drive (mL)
KVASCG : 0.20 : Microvessel formation rate (1/day)
KVASCD : 0.20 : Microvessel loss rate (1/day)
ECVR_SU: 20   : Sunitinib EC50 for VEGFR2 (ng/mL)
ECVR_RE: 500  : Regorafenib EC50 for VEGFR2 (ng/mL)
ECVR_RI: 1500 : Ripretinib EC50 for VEGFR2 (ng/mL)
HDEL   : 0.5  : Delivery exponent on relative microvascular density (-)

// ---- alpha-1 acid glycoprotein ---------------------------------------------
AGP0   : 1.0  : Reference AGP (g/L)
AGPF   : 1.0  : Individual AGP multiplier (-)
TAUAGP : 7.0  : AGP turnover (day)
WAGP   : 0.25 : Maximal fractional AGP rise with tumour burden (-)
HAGP   : 0.90 : Exponent linking AGP to apparent imatinib EC50 (-)

// ---- ctDNA -----------------------------------------------------------------
KSHED  : 2.2e-8 : ctDNA copies/mL released per cell death (-)
KCTEL  : 12.0   : ctDNA elimination (1/day)
CTWT   : 40     : Background wild-type cfDNA (copies/mL)

// ---- haematology and toxicity ----------------------------------------------
KTR    : 0.75   : Myeloid transit rate (1/day)
GAM    : 0.16   : Neutrophil feedback exponent (-)
ANC0   : 4.0    : Baseline ANC (10^9/L)
SLANCIM: 1.8e-5 : Imatinib slope on myeloid proliferation (mL/ng)
SLANCSU: 1.1e-3 : Sunitinib slope on myeloid proliferation (mL/ng)
SLANCRE: 1.1e-5 : Regorafenib slope on myeloid proliferation (mL/ng)
SLANCRI: 4.4e-5 : Ripretinib slope on myeloid proliferation (mL/ng)
SLANCAV: 7.5e-5 : Avapritinib slope on myeloid proliferation (mL/ng)
KMAPD  : 0.050  : Blood-pressure relaxation rate (1/day)
MAPMAX : 16.0   : Maximal MAP rise at full VEGFR2 occupancy (mmHg)
KHFD   : 0.055  : Hand-foot skin reaction relaxation rate (1/day)
KEDD   : 0.030  : Oedema relaxation rate (1/day)
KF4P   : 0.050  : Free T4 production (1/day)
KF4D   : 0.050  : Free T4 elimination (1/day)
EC50THY: 60     : Sunitinib EC50 for thyroid capillary regression (ng/mL)
EMAXTHY: 0.80   : Maximal suppression of T4 production (-)
KTSHP  : 0.070  : TSH production (1/day)
KTSHD  : 0.070  : TSH elimination (1/day)
TOXTH  : 1.05   : Composite toxicity index tolerated without dose reduction (-)
KDIN   : 0.045  : Rate of dose reduction above the threshold (1/day)
KDIU   : 0.015  : Rate of dose re-escalation below it (1/day)
DIMIN  : 0.50   : Minimum achieved dose intensity (-)
DIFIX  : 0      : 1 = hold dose intensity at 1 (turn the feedback off)

// ---- adjuvant setting ------------------------------------------------------
N_MICRO: 2.2e6  : FITTED occult cells after macroscopically complete resection (-)
N_DET  : 1.0e9  : Radiologically detectable recurrence (cells)

// -----------------------------------------------------------------------------
$CMT @annotated
A_IM : Imatinib gut (mg)
C_IM : Imatinib central (mg)
P_IM : Imatinib peripheral (mg)
M_IM : CGP74588 central (mg)
A_SU : Sunitinib gut (mg)
C_SU : Sunitinib central (mg)
P_SU : Sunitinib peripheral (mg)
M_SU : SU12662 central (mg)
A_RE : Regorafenib gut (mg)
C_RE : Regorafenib central (mg)
P_RE : Regorafenib peripheral (mg)
M_RE : M-2/M-5 central (mg)
A_RI : Ripretinib gut (mg)
C_RI : Ripretinib central (mg)
P_RI : Ripretinib peripheral (mg)
M_RI : DP-5439 central (mg)
A_AV : Avapritinib gut (mg)
C_AV : Avapritinib central (mg)
P_AV : Avapritinib peripheral (mg)
ENZ  : CYP3A4 enzyme pool (relative)
N0C  : Primary-genotype clone, cycling (cells)
N0Q  : Primary-genotype clone, quiescent (cells)
N1C  : Exon 13/14 ATP-pocket subclone, cycling (cells)
N1Q  : Exon 13/14 subclone, quiescent (cells)
N2C  : Exon 17/18 activation-loop subclone, cycling (cells)
N2Q  : Exon 17/18 subclone, quiescent (cells)
N3C  : KIT-independent bypass subclone, cycling (cells)
N3Q  : KIT-independent subclone, quiescent (cells)
VNEC : Non-viable tumour tissue (mL)
SLD  : Imaged sum of longest diameters (mm)
SUV  : FDG-PET SUVmax (-)
KI67 : Ki67 index (%)
VASC : Relative microvascular density (-)
VEGF : VEGF-A (relative)
AGP  : Alpha-1 acid glycoprotein (g/L)
CT0  : ctDNA, primary mutation (copies/mL)
CTR  : ctDNA, resistance mutations (copies/mL)
NPRO : Myeloid proliferating pool (10^9/L)
NT1  : Myeloid transit 1 (10^9/L)
NT2  : Myeloid transit 2 (10^9/L)
ANC  : Absolute neutrophil count (10^9/L)
MAPD : Mean arterial pressure offset (mmHg)
HFSR : Hand-foot skin reaction grade (-)
EDEM : Periorbital/peripheral oedema score (-)
FT4  : Free T4 (relative)
TSH  : TSH (relative)
DI   : Achieved dose intensity (-)
TOXC : Cumulative toxicity burden (-)
QDEP : Depth of dormancy of the reservoir (-)

// -----------------------------------------------------------------------------
$GLOBAL
#define CIM   (1000.0*(C_IM/V1_IM + M_IM/VM_IM))
#define CSU   (1000.0*(C_SU/V1_SU + M_SU/VM_SU))
#define CRE   (1000.0*(C_RE/V1_RE + M_RE/VM_RE))
#define CRI   (1000.0*(C_RI/V1_RI + M_RI/VM_RI))
#define CAV   (1000.0*(C_AV/V1_AV))
#define NTOT  (N0C+N0Q+N1C+N1Q+N2C+N2Q+N3C+N3Q)

// primary-clone EC50 vector, filled in $MAIN from GENO
double e0im, e0su, e0re, e0ri, e0av, kitdep0, kps0;

// coverage of one clone by the current five concentrations, and the resulting
// relative KIT pathway output.  Competitive and independent, so the coverages
// multiply; kitdep < 1 means part of the growth signal does not come from KIT.
double covr(double cim, double csu, double cre, double cri, double cav,
            double eim, double esu, double ere, double eri, double eav,
            double kitdep) {
  double a = 1.0;
  if (cim > 0) a /= (1.0 + cim/eim);
  if (csu > 0) a /= (1.0 + csu/esu);
  if (cre > 0) a /= (1.0 + cre/ere);
  if (cri > 0) a /= (1.0 + cri/eri);
  if (cav > 0) a /= (1.0 + cav/eav);
  double s = 1.0 - kitdep*(1.0 - a);
  if (s < 1e-6) s = 1e-6;
  if (s > 1.0)  s = 1.0;
  return s;
}
double kprol(double s, double kpmax, double kps, double sig50, double hp) {
  double sh = pow(s, hp);
  return kpmax*kps*sh/(sh + pow(sig50, hp));
}
double kdeath(double s, double kd0, double kamax, double ucrit, double ha) {
  double u = 1.0 - s;
  return kd0 + kamax*pow(u, ha)/(pow(u, ha) + pow(ucrit, ha));
}

// -----------------------------------------------------------------------------
$MAIN
// ---- primary-genotype EC50 vector -----------------------------------------
if (GENO == 2) {                       // KIT exon 9
  e0im = E9_IM; e0su = E9_SU; e0re = E9_RE; e0ri = E9_RI; e0av = E9_AV;
  kitdep0 = 1.0; kps0 = 1.00;
} else if (GENO == 3) {                // PDGFRA D842V
  e0im = D842_IM; e0su = D842_SU; e0re = D842_RE; e0ri = D842_RI; e0av = D842_AV;
  kitdep0 = 1.0; kps0 = 0.85;
} else if (GENO == 4) {                // SDH-deficient / KIT-PDGFRA wild type
  e0im = 1e7; e0su = 1e7; e0re = 1e7; e0ri = 1e7; e0av = 1e7;
  kitdep0 = 0.15; kps0 = 0.65;         // indolent and largely KIT-independent
} else {                               // KIT exon 11 (default)
  e0im = E11_IM; e0su = E11_SU; e0re = E11_RE; e0ri = E11_RI; e0av = E11_AV;
  kitdep0 = 1.0; kps0 = 1.00;
}

// ---- initial conditions ---------------------------------------------------
double n_start = (SETTING == 2) ? N_MICRO : V0*CELLSML;
// Resistant subclones are PRE-EXISTING: the mutation flux integrated over the
// pre-diagnostic history, which is TURNOVER divisions per surviving cell.
double f1 = MU_ATP*TURNOVER;
double f2 = MU_AL*TURNOVER;
double f3 = MU_BYP*TURNOVER;
N1C_0 = n_start*f1;
N2C_0 = n_start*f2;
N3C_0 = n_start*f3;
N0C_0 = n_start*(1.0 - f1 - f2 - f3);
double lam0 = KPMAX/(1.0 + pow(SIG50, HP)) - KD0;
VNEC_0 = FNEC0*KD0/(KRES + (lam0 > 1e-4 ? lam0 : 1e-4))*n_start/CELLSML;
double v_0 = n_start/CELLSML + VNEC_0;
SLD_0  = KSLD*pow(v_0, 1.0/3.0);
SUV_0  = SUV0;
KI67_0 = KI67MX;
ENZ_0  = 1.0;
VASC_0 = 1.0;
VEGF_0 = 1.0;
AGP_0  = AGP0*AGPF;
NPRO_0 = ANC0; NT1_0 = ANC0; NT2_0 = ANC0; ANC_0 = ANC0;
FT4_0  = 1.0; TSH_0 = 1.0; DI_0 = 1.0;
CT0_0  = KSHED*KD0*N0C_0/KCTEL;
CTR_0  = KSHED*KD0*(N1C_0 + N2C_0 + N3C_0)/KCTEL;

// -----------------------------------------------------------------------------
$ODE
// =========================== pharmacokinetics ==============================
double clim = CL_IM*ENZ*CLF;
dxdt_A_IM = -KA_IM*A_IM;
dxdt_C_IM = F_IM*KA_IM*A_IM - (clim/V1_IM)*C_IM - Q_IM*(C_IM/V1_IM - P_IM/VP_IM);
dxdt_P_IM = Q_IM*(C_IM/V1_IM - P_IM/VP_IM);
dxdt_M_IM = FM_IM*clim/V1_IM*C_IM - CLM_IM/VM_IM*M_IM;
dxdt_ENZ  = KENZ*((1.0 + EMAXENZ*CIM/(EC50ENZ + CIM))*INDF - ENZ);

dxdt_A_SU = -KA_SU*A_SU;
dxdt_C_SU = F_SU*KA_SU*A_SU - (CL_SU/V1_SU)*C_SU - Q_SU*(C_SU/V1_SU - P_SU/VP_SU);
dxdt_P_SU = Q_SU*(C_SU/V1_SU - P_SU/VP_SU);
dxdt_M_SU = FM_SU*CL_SU/V1_SU*C_SU - CLM_SU/VM_SU*M_SU;

dxdt_A_RE = -KA_RE*A_RE;
dxdt_C_RE = F_RE*KA_RE*A_RE - (CL_RE/V1_RE)*C_RE - Q_RE*(C_RE/V1_RE - P_RE/VP_RE);
dxdt_P_RE = Q_RE*(C_RE/V1_RE - P_RE/VP_RE);
dxdt_M_RE = FM_RE*CL_RE/V1_RE*C_RE - CLM_RE/VM_RE*M_RE;

dxdt_A_RI = -KA_RI*A_RI;
dxdt_C_RI = F_RI*KA_RI*A_RI - (CL_RI/V1_RI)*C_RI - Q_RI*(C_RI/V1_RI - P_RI/VP_RI);
dxdt_P_RI = Q_RI*(C_RI/V1_RI - P_RI/VP_RI);
dxdt_M_RI = FM_RI*CL_RI/V1_RI*C_RI - CLM_RI/VM_RI*M_RI;

dxdt_A_AV = -KA_AV*A_AV;
dxdt_C_AV = F_AV*KA_AV*A_AV - (CL_AV/V1_AV)*C_AV - Q_AV*(C_AV/V1_AV - P_AV/VP_AV);
dxdt_P_AV = Q_AV*(C_AV/V1_AV - P_AV/VP_AV);

// ===================== tumour interstitial availability ====================
double vasc  = (VASC > 1e-3) ? VASC : 1e-3;
double deliv = pow(vasc, HDEL);
if (deliv > 1.6) deliv = 1.6;
double cim = CIM*deliv, csu = CSU*deliv, cre = CRE*deliv;
double cri = CRI*deliv, cav = CAV*deliv;

// AGP binding raises the APPARENT imatinib EC50 of every clone [C3]
double agpf = pow((AGP > 0.2 ? AGP : 0.2)/AGP0, HAGP);

// =============== per-clone signalling, growth and death [C1][C2] ===========
double s0 = covr(cim, csu, cre, cri, cav, e0im*agpf, e0su, e0re, e0ri, e0av, kitdep0);
double s1 = covr(cim, csu, cre, cri, cav, C1_IM*agpf, C1_SU, C1_RE, C1_RI, C1_AV, KITDEP1);
double s2 = covr(cim, csu, cre, cri, cav, C2_IM*agpf, C2_SU, C2_RE, C2_RI, C2_AV, KITDEP2);
double s3 = covr(cim, csu, cre, cri, cav, C3_EC, C3_EC, C3_EC, C3_EC, C3_EC, KITDEP3);

double kp0 = kprol(s0, KPMAX, kps0, SIG50, HP);
double kp1 = kprol(s1, KPMAX, KPS1, SIG50, HP);
double kp2 = kprol(s2, KPMAX, KPS2, SIG50, HP);
double kp3 = kprol(s3, KPMAX, KPS3, SIG50, HP);
double kd0c = kdeath(s0, KD0, KAMAX, UCRIT, HA);
double kd1c = kdeath(s1, KD0, KAMAX, UCRIT, HA);
double kd2c = kdeath(s2, KD0, KAMAX, UCRIT, HA);
double kd3c = kdeath(s3, KD0, KAMAX, UCRIT, HA);

// exit from quiescence needs mitogenic signal steeply, and is slowed further by
// dormancy depth: this is what keeps the reservoir intact for years [C2]
double qd  = (QDEP < 0) ? 0 : ((QDEP > 1) ? 1 : QDEP);
double wake = KQOUT*(1.0 - PHI_D*qd);
double qi0 = KQIN*(1.0 - s0), qo0 = wake*pow(s0, NQ);
double qi1 = KQIN*(1.0 - s1), qo1 = wake*pow(s1, NQ);
double qi2 = KQIN*(1.0 - s2), qo2 = wake*pow(s2, NQ);
double qi3 = KQIN*(1.0 - s3), qo3 = wake*pow(s3, NQ);

double mu = MU_ATP + MU_AL + MU_BYP;
dxdt_N0C = (kp0*(1.0 - mu) - kd0c - qi0)*N0C + qo0*N0Q;
dxdt_N0Q = qi0*N0C - (qo0 + KDQ)*N0Q;
dxdt_N1C = (kp1 - kd1c - qi1)*N1C + qo1*N1Q + MU_ATP*kp0*N0C;
dxdt_N1Q = qi1*N1C - (qo1 + KDQ)*N1Q;
dxdt_N2C = (kp2 - kd2c - qi2)*N2C + qo2*N2Q + MU_AL*kp0*N0C;
dxdt_N2Q = qi2*N2C - (qo2 + KDQ)*N2Q;
dxdt_N3C = (kp3 - kd3c - qi3)*N3C + qo3*N3Q + MU_BYP*kp0*N0C;
dxdt_N3Q = qi3*N3C - (qo3 + KDQ)*N3Q;

// ======================= mass, imaging, ctDNA [C2] =========================
double ntot   = (NTOT > 1.0) ? NTOT : 1.0;
double dflux0 = kd0c*N0C + KDQ*N0Q;
double dfluxR = kd1c*N1C + KDQ*N1Q + kd2c*N2C + KDQ*N2Q + kd3c*N3C + KDQ*N3Q;
// only death IN EXCESS of physiological turnover leaves persistent non-viable
// tissue: a synchronous die-off overwhelms clearance, a trickle does not
double excess = fmax(0.0, kd0c - KD0)*N0C + fmax(0.0, kd1c - KD0)*N1C
              + fmax(0.0, kd2c - KD0)*N2C + fmax(0.0, kd3c - KD0)*N3C;
double basal  = KD0*(N0C + N1C + N2C + N3C) + KDQ*(N0Q + N1Q + N2Q + N3Q);
double vviab  = ntot/CELLSML;
double vnec   = (VNEC > 0) ? VNEC : 0;
double vtot   = vviab + vnec;
dxdt_VNEC = (FNEC*excess + FNEC0*basal)/CELLSML - KRES*vnec;
dxdt_SLD  = (KSLD*pow((vtot > 1e-9 ? vtot : 1e-9), 1.0/3.0) - SLD)/TAUSLD;

double sigbar = (s0*(N0C+N0Q) + s1*(N1C+N1Q) + s2*(N2C+N2Q) + s3*(N3C+N3Q))/ntot;
dxdt_SUV  = (SUV_BG + (SUV0 - SUV_BG)*sigbar - SUV)/TAUPET;
dxdt_KI67 = (KI67MX*(kp0*N0C + kp1*N1C + kp2*N2C + kp3*N3C)/(KPMAX*ntot) - KI67)/TAUKI;
dxdt_CT0  = KSHED*dflux0 - KCTEL*CT0;
dxdt_CTR  = KSHED*dfluxR - KCTEL*CTR;

// ======================= vasculature and delivery ==========================
double occv = CSU/(CSU + ECVR_SU) + CRE/(CRE + ECVR_RE) + CRI/(CRI + ECVR_RI);
if (occv > 0.95) occv = 0.95;
dxdt_VEGF = KVEGFP*(vtot/VREF + 0.05) - KVEGFD*VEGF;
dxdt_VASC = KVASCG*VEGF*(1.0 - occv) - KVASCD*VASC;
dxdt_AGP  = (AGP0*AGPF*(1.0 + WAGP*vtot/(vtot + 200.0)) - AGP)/TAUAGP;

// ============================== toxicity ===================================
double edrug = SLANCIM*CIM + SLANCSU*CSU + SLANCRE*CRE + SLANCRI*CRI + SLANCAV*CAV;
if (edrug > 0.90) edrug = 0.90;
double ancx = (ANC > 0.05) ? ANC : 0.05;
dxdt_NPRO = KTR*NPRO*((1.0 - edrug)*pow(ANC0/ancx, GAM) - 1.0);
dxdt_NT1  = KTR*(NPRO - NT1);
dxdt_NT2  = KTR*(NT1 - NT2);
dxdt_ANC  = KTR*(NT2 - ANC);
dxdt_MAPD = KMAPD*(MAPMAX*occv - MAPD);
double hfdrive = CRE/1200.0 + 0.50*CSU/70.0 + 0.25*CRI/800.0 + 0.05*CIM/1500.0;
dxdt_HFSR = KHFD*(hfdrive - HFSR);
dxdt_EDEM = KEDD*(CIM/1500.0 - EDEM);
dxdt_FT4  = KF4P*(1.0 - EMAXTHY*CSU/(CSU + EC50THY)) - KF4D*FT4;
dxdt_TSH  = KTSHP*pow(1.0/((FT4 > 0.15) ? FT4 : 0.15), 1.5) - KTSHD*TSH;

double tox = 0.60*fmax(0.0, 1.5 - ANC)/1.5 + 0.50*HFSR + 0.40*MAPD/15.0
           + 0.30*EDEM + 0.20*fmax(0.0, TSH - 2.0)/3.0;
dxdt_TOXC = tox - 0.20*TOXC;
// toxicity -> dose reduction -> LESS COVERAGE -> faster selection.  This loop
// is the reason a drug's tolerability is part of its efficacy.
dxdt_DI = (DIFIX > 0.5) ? 0.0
        : (-KDIN*DI*fmax(0.0, tox - TOXTH) + KDIU*(1.0 - DI)*fmax(0.0, TOXTH - tox));

// dormancy depth follows the time-averaged signal, symmetrically but slowly
dxdt_QDEP = KDEEP*(1.0 - sigbar)*(1.0 - qd) - KSHAL*sigbar*qd;

// -----------------------------------------------------------------------------
$TABLE
double C_IMAT  = CIM;                        // ng/mL, parent + CGP74588
double C_SUNI  = CSU;
double C_REGO  = CRE;
double C_RIPR  = CRI;
double C_AVAP  = CAV;
double N_TOT   = NTOT;
double N_PRIM  = N0C + N0Q;
double N_ATP   = N1C + N1Q;
double N_LOOP  = N2C + N2Q;
double N_BYP   = N3C + N3Q;
double N_QUIES = N0Q + N1Q + N2Q + N3Q;
double V_VIAB  = N_TOT/CELLSML;
double V_TOT   = V_VIAB + VNEC;
double V_FRAC  = V_VIAB/((V_TOT > 1e-9) ? V_TOT : 1e-9);   // viable fraction of the imaged mass
double VAF_R   = CTR/(CT0 + CTR + CTWT);                   // resistance-mutation VAF
double CTDNA   = CT0 + CTR;
double MAP_TOT = 92.0 + MAPD;
double DETECT  = (N_TOT >= N_DET) ? 1 : 0;                 // adjuvant: radiological recurrence

$CAPTURE @annotated
C_IMAT  : Imatinib + CGP74588 (ng/mL)
C_SUNI  : Sunitinib + SU12662 (ng/mL)
C_REGO  : Regorafenib + M-2/M-5 (ng/mL)
C_RIPR  : Ripretinib + DP-5439 (ng/mL)
C_AVAP  : Avapritinib (ng/mL)
N_TOT   : Total viable tumour cells (-)
N_PRIM  : Primary-genotype clone (cells)
N_ATP   : Exon 13/14 ATP-pocket subclone (cells)
N_LOOP  : Exon 17/18 activation-loop subclone (cells)
N_BYP   : KIT-independent bypass subclone (cells)
N_QUIES : Quiescent reservoir (cells)
V_VIAB  : Viable tumour volume (mL)
V_TOT   : Imaged tumour volume (mL)
V_FRAC  : Viable fraction of the imaged mass (-)
VAF_R   : Resistance-mutation ctDNA VAF (-)
CTDNA   : Total ctDNA (copies/mL)
MAP_TOT : Mean arterial pressure (mmHg)
DETECT  : Recurrence detectable (adjuvant setting) (0/1)
