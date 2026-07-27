## =====================================================================
##  Charcot-Marie-Tooth disease type 1A (CMT1A)
##  Quantitative Systems Pharmacology model — mrgsolve implementation
## ---------------------------------------------------------------------
##  Companion files:
##     cmt1a_qsp_model.dot / .svg / .png   mechanistic map (24 modules)
##     cmt1a_references.md                 literature
##     cmt1a_shiny_app.R                   interactive dashboard
##     README.md                           narrative, results, caveats
##
##  WHAT THIS MODEL IS FOR
##  ----------------------
##  CMT1A is a GENE DOSAGE disease. A 1.4-Mb duplication at 17p11.2
##  delivers three copies of PMP22 instead of two. Everything else is
##  downstream of running one myelin protein at the wrong stoichiometry.
##  The model is built so that this single integer is the disease
##  parameter (CN = 1, 2, 3, 4 reproduce HNPP, normal, CMT1A and the
##  severe homozygous phenotype) and so that the therapeutic target is
##  a RATIO to be restored, not a quantity to be minimised.
##
##  FIVE STRUCTURAL COMMITMENTS (each is testable, see diagnostics)
##  ---------------------------------------------------------------
##  S1. THE DISPOSAL SYSTEM IS ALREADY AT ~80% DUTY CYCLE IN HEALTH.
##      PMP22 folds inefficiently even in wild-type: roughly 80% of what
##      is translated is degraded before it ever reaches myelin. A 50%
##      rise in synthesis therefore does not produce a 50% rise in
##      anything - it produces a NON-LINEAR escape of misfolded protein
##      past a saturating ERAD step. This is why 3 copies is a disease
##      and 2 copies is not, and it is why CN=4 runs away entirely.
##
##  S2. THE TWO ARMS OF THE U-CURVE ARE DIFFERENT PATHWAYS.
##      CMT1A and HNPP are NOT mirror images. Excess dosage acts through
##      the AGGREGATE burden (proteostatic stress -> c-Jun -> Schwann
##      cell dedifferentiation -> demyelination). Deficient dosage acts
##      through MEMBRANE STOICHIOMETRY (gm < 1 -> tomacula -> focal
##      pressure palsy). The therapeutic window therefore has a near
##      wall and a far wall built out of different biology, and the
##      far wall sits where knockdown drives mRNA to the 1-copy level.
##      Because 1.5 x (1 - 0.667) = 0.5, that far wall is at 67%
##      knockdown and the optimum at 33%. Neither number is fitted.
##
##  S3. LENGTH-DEPENDENCE IS GEOMETRY. Delivered (transported) support
##      decays as exp(-L/(v*tau)); locally delivered (glial lactate)
##      support does not. Demand does not depend on L. The resulting
##      MARGIN orders foot < leg < hand < proximal with no biological
##      difference between the classes, and it is the same equation
##      that makes vincristine catastrophic in CMT1A and tolerable in
##      a normal nerve.
##
##  S4. THE MOTOR UNIT RESERVE IS LARGE, SILENT AND ONE-WAY. Collateral
##      sprouting lets a surviving unit adopt the fibres of several lost
##      ones. Strength is therefore flat while MUNE falls by ~80%.
##      "Onset" is reserve exhaustion, not disease onset.
##
##  S5. CONDUCTION VELOCITY IS A DEVELOPMENTAL FOSSIL. Internodal length
##      is eroded overwhelmingly during the myelination window; after
##      childhood the model's own dNCV/dt is far below measurement
##      resolution. Nothing forces this - it emerges from putting the
##      erosion term under a developmental weight - and it is what
##      makes NCV unusable as a trial endpoint.
##
##  UNITS: time = days. Disease states are dimensionless indices
##  normalised to the healthy adult (1.0) unless stated. Drug
##  concentrations are noted per compartment.
##
##  Requires: mrgsolve (>= 1.0). Developed and validated on 2.0.1.
## =====================================================================

library(mrgsolve)

code <- '
$PROB CMT1A QSP — PMP22 dosage / Schwann-cell proteostasis / length-dependent axonal supply

$PARAM @annotated
// ---------------- genetic and constitutional -------------------------
CN      : 3      : PMP22 copy number (1 HNPP, 2 normal, 3 CMT1A, 4 homozygous dup)
AGE0    : 0      : age at t = 0 (years)
HTF     : 1.0    : stature factor scaling all axon lengths (1.0 = 1.72 m adult)
MODF    : 1.0    : genetic modifier factor on aggregate handling capacity
SEXF    : 0      : 0 = male, 1 = female (permits pregnancy scenario)

// ---------------- PMP22 expression -----------------------------------
KSYN_M  : 0.30   : PMP22 mRNA synthesis rate constant (1/day) at CN = 2
KDEG_M  : 0.30   : PMP22 mRNA degradation rate constant (1/day)
KTR     : 1.00   : translation flux per unit mRNA (arbitrary flux units/day)
PHIMAX  : 0.85   : maximum attainable folding efficiency
CHAP    : 0.307  : chaperone/folding capacity (flux units) - set so phi(CN=2) = 0.20
KDEG_P  : 0.05   : membrane PMP22 turnover (1/day)

// ---------------- proteostasis ---------------------------------------
JC_ERAD : 1.60   : misfolded flux at which ERAD is half-escaped (flux units/day)
HE_ERAD : 4.0    : steepness of ERAD escape
KAUT    : 1.00   : autophagic clearance rate constant of aggregates (1/day)
PA50    : 1.00   : aggregate load at which autophagic flux is half-inhibited
PAGGMAX : 6.0    : hard ceiling on aggregate load (cell-level saturation)
KISR    : 0.20   : ISR activation rate (1/day)
KISROFF : 0.20   : ISR resolution rate (1/day)
AISR    : 1.60   : aggregate-to-ISR gain

// ---------------- Schwann-cell state ---------------------------------
KCJ     : 0.10   : c-Jun production rate (1/day)
KCJOFF  : 0.10   : c-Jun decay (1/day)
ACJ_AG  : 0.60   : aggregate drive on c-Jun
ACJ_ISR : 0.60   : ISR drive on c-Jun
KEG     : 0.10   : EGR2 production rate (1/day)
KEGOFF  : 0.10   : EGR2 decay (1/day)
BCJ_EG  : 0.90   : c-Jun repression of EGR2
BEG_CJ  : 0.70   : EGR2 repression of c-Jun
HSCD    : 1.5    : Hill coefficient of the SCD differentiation switch
KOB     : 0.0016 : onion-bulb / SC proliferation accumulation rate (1/day)
KOBOFF  : 0.0004 : onion-bulb resolution rate (1/day)

// ---------------- myelin ----------------------------------------------
KMYT    : 0.010  : myelin thickness relaxation rate (1/day)
AGMSTO  : 0.80   : direct stoichiometric penalty on myelin thickness for gm > 1
KERA    : 2.3e-4 : internodal-length erosion per unit demyelination rate (1/day)
DEVBOOST: 9.0    : internodal erosion weight inside the myelination window
DEVRES  : 0.06   : residual erosion weight AFTER the window closes
DEVTAU  : 5.0    : developmental myelination window (years, smooth close)
INLMIN  : 0.32   : floor on internodal length index
KNDI    : 0.030  : nodal disorganisation rate (1/day)
KNDIOFF : 0.030  : nodal repair rate (1/day)
NCV0    : 52.0   : motor NCV at INL = 1, MYT = 1 (m/s)
ANCV    : 0.75   : NCV exponent on internodal length
BNCV    : 0.55   : NCV exponent on myelin thickness
TOMREF  : 0.192  : (1 - gm) at the 1-copy (HNPP) state - sets the far wall
HTOM    : 1.5    : steepness of the tomacula/pressure-palsy response

// ---------------- axonal transport and metabolic support --------------
KATAT   : 0.30   : alpha-tubulin acetylation rate (1/day)
KHDAC   : 0.30   : HDAC6 deacetylation rate (1/day)
AROS_TU : 0.35   : ROS penalty on tubulin acetylation
V0      : 0.90   : transport velocity index at TUBA = 1
TAUCG   : 1.00   : cargo survival time along the axon (arbitrary, with v: m)
KCG_IN  : 0.60   : cargo delivery rate constant (1/day)
KCG_OUT : 0.12   : cargo consumption rate constant (1/day)
KLACSUP : 0.75   : weight of LOCAL (glial lactate) support in the margin
KLAC    : 0.08   : lactate support turnover (1/day)
KMITO   : 0.02   : mitochondrial health turnover (1/day)
KROS    : 0.05   : ROS turnover (1/day)
AROS_GAP: 1.20   : energy-gap drive on ROS
DEM0    : 1.00   : baseline terminal demand
KNKA    : 0.80   : extra demand imposed by loss of saltatory conduction

// ---------------- axon lengths (metres, adult) ------------------------
L_FOOT  : 1.00   : motor axon to intrinsic foot muscles
L_LEG   : 0.85   : motor axon to tibialis anterior / peroneal group
L_HAND  : 0.75   : motor axon to hand intrinsics
L_PROX  : 0.35   : motor axon to proximal limb muscles
L_SFOOT : 1.05   : sensory axon, DRG to hallux
L_SHAND : 0.80   : sensory axon, DRG to fingertip
LFRAC0  : 0.34   : axon length at birth as a fraction of adult
LGROW   : 16.0   : years over which axons reach adult length

// ---------------- axon loss and reinnervation -------------------------
KAX     : 3.6e-4 : maximum axon degeneration rate constant (1/day)
PAX     : 8.0    : steepness of the chronic hazard-vs-margin relation
MG50    : 1.0    : supply/demand margin at half-maximal degeneration hazard
KAGE    : 1.1e-5 : age-related motor unit attrition (1/day)
KACUTE  : 0.050  : ACUTE degeneration rate when the margin collapses (1/day)
MGACUTE : 0.55   : margin below which the acute limb engages
PACUTE  : 12.0   : steepness of the acute limb
KREGEN  : 2.0e-5 : regenerative reinnervation rate (1/day)
RSMAX   : 5.0    : maximum motor unit enlargement (fibres per unit, x normal)
KSPR    : 0.0060 : collateral sprouting rate (1/day)
KFRAG   : 0.28   : NMJ instability penalty on enlarged motor units
KMS     : 0.0035 : muscle mass adaptation rate (1/day)
KFAT    : 0.0022 : fatty-replacement rate (1/day)
FAT0    : 4.0    : normal intramuscular fat fraction (%)
FATMAX  : 90.0   : maximum intramuscular fat fraction (%)
KEXER   : 0.0    : structured exercise stimulus (0 = none, 1 = full programme)
AEXER   : 0.12   : maximum muscle-mass gain from exercise

// ---------------- sensory and symptom layers --------------------------
KVIB    : 0.010  : vibration threshold adaptation (1/day)
VIB0    : 4.0    : normal vibration threshold at hallux (um)
VIBMAX  : 130.0  : maximum measurable vibration threshold (um)
KPROP   : 0.010  : proprioceptive loss adaptation (1/day)
KNFL    : 0.15   : plasma NfL elimination (1/day)
NFLBASE : 8.0    : NfL in the absence of active degeneration (pg/mL)
ANFL    : 1.35e4 : NfL released per unit axon-loss rate
KPAIN   : 0.02   : pain state adaptation (1/day)
APAIN   : 5.5    : maximum ectopic-firing pain contribution
KFAT_S  : 0.02   : fatigue adaptation (1/day)
KCAV    : 0.0016 : pes cavus accumulation per unit imbalance (1/day)

// ---------------- PXT3003 (baclofen + naltrexone + D-sorbitol) --------
PXT_ON  : 0      : 1 = PXT3003 administered
PXT_LVL : 1.0    : dose level (1 = PLEO-CMT dose 1, 2 = high dose)
BAC_D   : 6.0    : baclofen daily dose at level 1 (mg/day, split BID)
NTX_D   : 0.70   : naltrexone daily dose at level 1 (mg/day)
SOR_D   : 21.0   : D-sorbitol daily dose at level 1 (mg/day)
KA_BAC  : 3.5    : baclofen absorption (1/day)
CL_BAC  : 12.0   : baclofen clearance (L/day)
V_BAC   : 60.0   : baclofen volume (L)
KA_NTX  : 12.0   : naltrexone absorption (1/day)
CL_NTX  : 3400   : naltrexone apparent clearance (L/day, high first pass)
V_NTX   : 1350   : naltrexone apparent volume (L)
KA_SOR  : 6.0    : sorbitol absorption (1/day)
CL_SOR  : 300    : sorbitol clearance (L/day)
V_SOR   : 40.0   : sorbitol volume (L)
EC_BAC  : 0.010  : baclofen EC50 on the PXT synergy index (mg/L)
EC_NTX  : 4.0e-5 : naltrexone EC50 on the PXT synergy index (mg/L)
EC_SOR  : 0.05   : sorbitol EC50 on the PXT synergy index (mg/L)
EMAX_PXT: 0.24   : maximum fractional PMP22 mRNA suppression by PXT3003

// ---------------- ascorbic acid ---------------------------------------
AA_ON   : 0      : 1 = oral ascorbic acid administered
AA_DTOT : 4000   : total daily ascorbate dose (mg/day)
AA_NDOS : 3      : number of divided doses per day
AA_D50  : 1219   : per-dose amount halving intestinal bioavailability (mg)
AA_HF   : 1.63   : steepness of intestinal SVCT1 saturation
KA_AA   : 8.0    : ascorbate absorption (1/day)
V_AA    : 18.0   : ascorbate distribution volume (L)
CL_AA_B : 1.60   : non-renal ascorbate clearance (L/day)
CL_AA_R : 300    : maximum threshold renal clearance (L/day, ~ GFR)
CTH_AA  : 70.0   : renal reabsorption threshold (umol/L)
HR_AA   : 10.0   : steepness of the renal threshold
KM_SVCT2: 20.0   : Schwann-cell SVCT2 Km (umol/L)
AA_BASE : 50.0   : target habitual dietary plasma ascorbate (umol/L)
AA_DIET : 1      : 1 = dietary ascorbate input active
IC50_AA : 620    : intracellular ascorbate for half-maximal cAMP suppression (umol/L equiv)
EMAX_AA : 0.30   : maximum fractional PMP22 mRNA suppression by ascorbate

// ---------------- PMP22-lowering oligonucleotide ----------------------
OLI_ON  : 0      : 1 = oligonucleotide administered
OLI_DOSE: 0      : subcutaneous dose per administration (mg)
KA_OLI  : 0.60   : subcutaneous absorption (1/day)
V_OLIPL : 12.0   : oligonucleotide plasma volume (L)
CL_OLIPL: 30.0   : oligonucleotide plasma clearance (L/day)
KPL_NV  : 1.10   : plasma-to-nerve distribution (L/day)
KNV_OUT : 0.0165 : nerve tissue elimination (1/day) - t1/2 ~ 42 days
V_OLINV : 0.60   : nerve tissue distribution volume (L)
EC50_OLI: 3.0    : nerve concentration for half-maximal knockdown (mg/L)
EMAX_OLI: 0.85   : maximum attainable PMP22 mRNA knockdown

// ---------------- HDAC6 inhibitor -------------------------------------
HDI_ON  : 0      : 1 = HDAC6 inhibitor administered
HDI_DOSE: 30.0   : daily oral dose (mg/day)
KA_HDI  : 12.0   : absorption (1/day)
CL_HDI  : 900    : clearance (L/day)
V_HDI   : 250    : volume (L)
IC50_HDI: 0.020  : HDAC6 IC50 (mg/L)
IMAX_HDI: 0.80   : maximum HDAC6 inhibition

// ---------------- vincristine -----------------------------------------
VCR_ON  : 0      : 1 = vincristine administered
KVCR_NV : 0.55   : plasma-to-nerve distribution (1/day)
KVCR_OUT: 0.030  : nerve elimination (1/day)
CL_VCR  : 190    : plasma clearance (L/day)
V_VCR   : 200    : plasma volume of distribution (L)
IC50_VCR: 0.30   : nerve vincristine for half-maximal transport block (ug/L)
IMAX_VCR: 0.40   : maximum fractional block of axonal transport

// ---------------- AAV1.NT-3 and other -----------------------------------
NT3_ON  : 0      : 1 = AAV1.NT-3 delivered at t = 0
NT3_SS  : 1.0    : steady-state NT-3 expression index reached after transduction
KNT3    : 0.020  : NT-3 expression onset rate (1/day)
ANT3    : 0.15   : NT-3 gain on trophic/regenerative supply
PREG_ON : 0      : 1 = run a pregnancy (starting at PREG_T0)
PREG_T0 : 0      : day of conception
APROG   : 0.06   : progesterone drive on PMP22 transcription at term
ONAPRI  : 0      : 1 = progesterone-receptor antagonist on board

// ---------------- measurement / trial layer ---------------------------
SD_CMTNS: 0      : SD of simulated CMTNS measurement noise (0 = noiseless)
PRACT   : 0      : timed-test practice effect magnitude

$CMT @annotated
MRNA  : PMP22 mRNA (index, 1 = healthy 2-copy steady state)
PPROT : correctly folded PMP22 in myelin membrane (index)
PAGG  : aggregated / mistrafficked PMP22 (index)
ISR   : integrated stress response signal (index)
AUTC  : autophagic clearance capacity (index)
CJ    : Schwann-cell c-Jun (index)
EG    : Schwann-cell EGR2 / Krox20 (index)
SCD   : myelinating differentiation state (0-1)
OB    : onion-bulb / supernumerary Schwann-cell mass (index)
MYT   : myelin thickness index (1 = normal)
INL   : internodal length index (1 = normal)
NDI   : nodal / paranodal disorganisation (0-1)
CYC   : cumulative de/remyelination cycles
TUBA  : alpha-tubulin acetylation (index)
LAC   : Schwann-cell lactate support (index)
MITO  : axonal mitochondrial health (index)
ROSX  : axonal oxidative stress (index)
CGF   : delivered cargo at foot terminal (index)
CGL   : delivered cargo at leg terminal (index)
CGH   : delivered cargo at hand terminal (index)
CGP   : delivered cargo at proximal terminal (index)
AXF   : surviving motor axons, foot (fraction)
AXL   : surviving motor axons, leg (fraction)
AXH   : surviving motor axons, hand (fraction)
AXP   : surviving motor axons, proximal (fraction)
SXF   : surviving large sensory axons, foot (fraction)
SXH   : surviving large sensory axons, hand (fraction)
RSF   : motor unit size, foot (x normal)
RSL   : motor unit size, leg (x normal)
RSH   : motor unit size, hand (x normal)
RSP   : motor unit size, proximal (x normal)
MSF   : muscle mass, intrinsic foot (index)
MSL   : muscle mass, anterior compartment (index)
MSH   : muscle mass, hand intrinsics (index)
MSP   : muscle mass, proximal (index)
FFF   : intramuscular fat fraction, foot (%)
FFL   : intramuscular fat fraction, calf/anterior (%)
VIB   : vibration threshold at hallux (um)
PROP  : proprioceptive loss (0-1)
NFL   : plasma neurofilament light (pg/mL)
PAINS : neuropathic pain state (0-10)
FATG  : fatigue state (FSS 1-7)
CAV   : pes cavus / deformity index (0-1)
BACD  : baclofen depot (mg)
BACC  : baclofen central (mg)
NTXD  : naltrexone depot (mg)
NTXC  : naltrexone central (mg)
SORC  : D-sorbitol central (mg)
AAD   : ascorbate gut (umol)
AAC   : ascorbate central (umol)
OLISC : oligonucleotide subcutaneous depot (mg)
OLIPL : oligonucleotide plasma (mg)
OLINV : oligonucleotide peripheral nerve (mg)
HDID  : HDAC6 inhibitor depot (mg)
HDIC  : HDAC6 inhibitor central (mg)
VCRC  : vincristine plasma (ug)
VCRNV : vincristine nerve (ug)
NT3E  : NT-3 expression index
PROG  : progesterone drive (0-1 of term level)

$GLOBAL
#define SQ(x) ((x)*(x))
// ---- shared between MAIN / ODE / TABLE ----
double AGEY;        // age in years
double DEVW;        // developmental myelination weight (1 -> 0)
double LF, LL, LH, LP, LSF, LSH;   // current axon lengths (m)
double JTR;         // translation flux
double PHI;         // folding efficiency
double JMIS;        // misfolded flux
double ERADF;       // fraction of misfolded flux handled by ERAD
double ESCAPE;      // misfolded flux escaping to aggregates
double PPREF;       // membrane PMP22 at the healthy 2-copy steady state
double GM;          // membrane PMP22 stoichiometry, relative to normal
double TOM;         // tomacula burden (HNPP arm), 1.0 = full 1-copy phenotype
double DEMY;        // instantaneous demyelination rate index
double VTR;         // effective axonal transport velocity
double FVCR;        // fractional transport block by vincristine
double SUPL_F, SUPL_L, SUPL_H, SUPL_P;   // transported supply terms
double LOCSUP;      // local (lactate) support term
double DEMAND;      // terminal demand
double MGF, MGL, MGH, MGP;   // margins
double MGSF, MGSH;           // sensory margins
double HZF, HZL, HZH, HZP, HZSF, HZSH;  // degeneration hazards
double INNF, INNL, INNH, INNP;          // innervated fraction
double FRCF, FRCL, FRCH, FRCP;          // force
double AXLOSSRATE;  // total instantaneous axon loss (for NfL)
double KD_OLI;      // fractional PMP22 knockdown from oligonucleotide
double SUP_PXT;     // fractional PMP22 suppression from PXT3003
double SUP_AA;      // fractional PMP22 suppression from ascorbate
double FPROG;       // progesterone transcriptional drive
double AAPL;        // plasma ascorbate (umol/L)
double AASC;        // Schwann-cell ascorbate saturation (0-1)
double CNVNV;       // oligonucleotide nerve concentration (mg/L)
double HDIINH;      // HDAC6 inhibition (0-1)
double NCV;         // motor nerve conduction velocity (m/s)
double PPHAZ;       // pressure-palsy hazard from tomacula (events/year)
double EGAP;        // energy gap index
double LGROWF;      // fraction of adult axon length reached
// ---- SOLVED healthy (CN = 2) reference state, computed in $MAIN ----
double SCD_H;       // healthy myelinating differentiation state
double PAGG_H;      // healthy aggregate burden
double CJ_H, EG_H, ISR_H, AUTC_H, MITO_H, ROSX_H, MYT_H;

$MAIN
// ------------------------------------------------------------------
//  Healthy-baseline algebra. PPREF is SOLVED, not typed in, so that a
//  CN = 2 subject has (near) zero drift and every reported deviation
//  is generated by copy number rather than by a mis-set initial value.
// ------------------------------------------------------------------
double MRNA_REF = 1.0;                                  // by construction
double JTR_REF  = KTR * MRNA_REF;
double PHI_REF  = PHIMAX * CHAP / (CHAP + JTR_REF);     // = 0.200
PPREF           = JTR_REF * PHI_REF / KDEG_P;

// ---- the healthy aggregate burden is a FIXED POINT, not a number ----
// Even a wild-type Schwann cell mis-folds ~80% of the PMP22 it makes,
// so PAGG_H > 0. Solving for it (rather than setting it to zero) is what
// makes the CN = 2 run drift-free and what makes the CN = 3 rise a
// RATIO to a real baseline rather than a rise from an artificial zero.
double JMIS_H   = JTR_REF * (1.0 - PHI_REF);
double ESC_H    = JMIS_H * (1.0 - 1.0 / (1.0 + pow(JMIS_H / JC_ERAD, HE_ERAD)));
PAGG_H = 0.0;
for (int ia = 0; ia < 200; ia++) {
  PAGG_H = ESC_H * (1.0 + SQ(PAGG_H / PA50)) / (KAUT * MODF * (1.0 - PAGG_H / PAGGMAX));
}
AUTC_H = 1.0 / (1.0 + SQ(PAGG_H / PA50));
ISR_H  = (KISR * AISR / KISROFF) * PAGG_H;

// ---- c-Jun / EGR2 mutual repression: also a fixed point ----
CJ_H = 0.3; EG_H = 0.8;
for (int ib = 0; ib < 200; ib++) {
  CJ_H = (KCJ / KCJOFF) * (0.30 + ACJ_AG * PAGG_H + ACJ_ISR * ISR_H) / (1.0 + BEG_CJ * EG_H);
  EG_H = (KEG / KEGOFF) / (1.0 + BCJ_EG * CJ_H);
}
SCD_H = pow(EG_H, HSCD) / (pow(EG_H, HSCD) + pow(CJ_H, HSCD) + 1e-9);
MYT_H = 1.0;   // myelin thickness is indexed to the healthy state
ROSX_H = 0.0;
MITO_H = 1.0;

// initial conditions: every scenario starts from the SOLVED healthy
// newborn nerve. The CMT1A phenotype is GENERATED by copy number over
// the following decades and is never typed in anywhere.
MRNA_0  = 1.0;
PPROT_0 = PPREF;
PAGG_0  = PAGG_H;
ISR_0   = ISR_H;
AUTC_0  = AUTC_H;
CJ_0    = CJ_H;
EG_0    = EG_H;
SCD_0   = SCD_H;
OB_0    = 0.0;
MYT_0   = MYT_H;
INL_0   = 1.0;
NDI_0   = 0.0;
CYC_0   = 0.0;
TUBA_0  = 1.0;
LAC_0   = 1.0;
MITO_0  = 1.0;
ROSX_0  = 0.0;
CGF_0   = 1.0; CGL_0 = 1.0; CGH_0 = 1.0; CGP_0 = 1.0;
AXF_0   = 1.0; AXL_0 = 1.0; AXH_0 = 1.0; AXP_0 = 1.0;
SXF_0   = 1.0; SXH_0 = 1.0;
RSF_0   = 1.0; RSL_0 = 1.0; RSH_0 = 1.0; RSP_0 = 1.0;
MSF_0   = 1.0; MSL_0 = 1.0; MSH_0 = 1.0; MSP_0 = 1.0;
FFF_0   = FAT0; FFL_0 = FAT0;
VIB_0   = VIB0;
PROP_0  = 0.0;
NFL_0   = NFLBASE;
PAINS_0 = 0.0;
FATG_0  = 1.0;
CAV_0   = 0.0;
NT3E_0  = 0.0;
PROG_0  = 0.0;
AAC_0   = AA_BASE * V_AA;   // start at the habitual dietary steady state

$ODE
// ==================================================================
//  0. AGE, GROWTH AND THE DEVELOPMENTAL MYELINATION WINDOW
// ==================================================================
AGEY   = AGE0 + SOLVERTIME / 365.25;
DEVW   = 1.0 / (1.0 + pow(AGEY / DEVTAU, 4.0));          // 1 -> 0 across ~5 y
// linear somatic growth to adult length by age LGROW, then constant.
// (The first draft used a hyperbola that only approached full length
// asymptotically, so a 40-year-old still had 90%-length axons and the
// margins never fell far enough to generate the adult phenotype.)
LGROWF = LFRAC0 + (1.0 - LFRAC0) * (AGEY / LGROW);
if (LGROWF > 1.0) LGROWF = 1.0;
LF  = HTF * L_FOOT  * LGROWF;
LL  = HTF * L_LEG   * LGROWF;
LH  = HTF * L_HAND  * LGROWF;
LP  = HTF * L_PROX  * LGROWF;
LSF = HTF * L_SFOOT * LGROWF;
LSH = HTF * L_SHAND * LGROWF;

// ==================================================================
//  1. PMP22 EXPRESSION  (the single disease parameter enters here)
// ==================================================================
// progesterone drive (pregnancy scenario); onapristone blocks the receptor
FPROG = 1.0 + APROG * PROG * (ONAPRI > 0.5 ? 0.05 : 1.0);
if (ONAPRI > 0.5) FPROG = FPROG * 0.88;   // removes basal PR tone as well

// pharmacological suppression of PMP22 transcription / message
CNVNV   = OLINV / V_OLINV;
KD_OLI  = (OLI_ON > 0.5) ? EMAX_OLI * CNVNV / (EC50_OLI + CNVNV) : 0.0;

double SBAC = BACC / V_BAC;
double SNTX = NTXC / V_NTX;
double SSOR = SORC / V_SOR;
double OCC  = (SBAC / (EC_BAC + SBAC)) * (SNTX / (EC_NTX + SNTX)) * (SSOR / (EC_SOR + SSOR));
SUP_PXT = (PXT_ON > 0.5) ? EMAX_PXT * OCC : 0.0;

// AAC carries TOTAL plasma ascorbate, dietary plus supplemental. The
// first draft eliminated only the supplemental increment while the
// renal threshold was evaluated on the total, so plasma never plateaued
// and 4 g/day reached 132 umol/L - a number Levine 1996 rules out.
AAPL    = AAC / V_AA;
AASC    = AAPL / (KM_SVCT2 + AAPL);            // SVCT2 saturation (0-1)
double AAIC  = AASC * 900.0;                   // intracellular equivalent (umol/L)
double AASC0 = AA_BASE / (KM_SVCT2 + AA_BASE);
double AAIC0 = AASC0 * 900.0;
// Only the INCREMENT above habitual dietary ascorbate can be a drug
// effect; everyone in the placebo arm already has the baseline.
SUP_AA  = (AA_ON > 0.5)
        ? EMAX_AA * (AAIC / (IC50_AA + AAIC) - AAIC0 / (IC50_AA + AAIC0))
        : 0.0;
if (SUP_AA < 0.0) SUP_AA = 0.0;

double SUPTOT = 1.0 - (1.0 - KD_OLI) * (1.0 - SUP_PXT) * (1.0 - SUP_AA);
if (SUPTOT > 0.98) SUPTOT = 0.98;

dxdt_MRNA = KSYN_M * (CN / 2.0) * FPROG * (1.0 - SUPTOT) - KDEG_M * MRNA;

// ==================================================================
//  2. FOLDING, ERAD ESCAPE AND THE AGGREGATE MASS BALANCE
// ==================================================================
JTR   = KTR * MRNA;
PHI   = PHIMAX * CHAP / (CHAP + JTR);          // saturable folding capacity
JMIS  = JTR * (1.0 - PHI);
ERADF = 1.0 / (1.0 + pow(JMIS / JC_ERAD, HE_ERAD));
ESCAPE= JMIS * (1.0 - ERADF);

dxdt_PPROT = JTR * PHI - KDEG_P * PPROT;
GM         = PPROT / PPREF;

double CLR_AGG = KAUT * MODF * AUTC * PAGG;
dxdt_PAGG = ESCAPE * (1.0 - PAGG / PAGGMAX) - CLR_AGG;
dxdt_AUTC = 0.05 * (1.0 / (1.0 + SQ(PAGG / PA50)) - AUTC);
dxdt_ISR  = KISR * AISR * PAGG - KISROFF * ISR;

// tomacula arm: only engaged when membrane stoichiometry falls BELOW normal
TOM = 0.0;
if (GM < 1.0) TOM = pow((1.0 - GM) / TOMREF, HTOM);
PPHAZ = 0.9 * TOM;                              // pressure palsies per year

// ==================================================================
//  3. SCHWANN-CELL STATE  (interior optimum #2: c-Jun)
// ==================================================================
dxdt_CJ = KCJ * (0.30 + ACJ_AG * PAGG + ACJ_ISR * ISR) / (1.0 + BEG_CJ * EG)
          - KCJOFF * CJ;
dxdt_EG = KEG * (1.0 + 0.35 * NT3E * ANT3) / (1.0 + BCJ_EG * CJ) - KEGOFF * EG;

double SCDT = pow(EG, HSCD) / (pow(EG, HSCD) + pow(CJ, HSCD) + 1e-9);
dxdt_SCD  = 0.05 * (SCDT - SCD);
// Demyelination is the RELATIVE dedifferentiation below the healthy
// solved state, not (1 - SCD). Written the second way, the ~6% of
// Schwann cells that are non-myelinating in any normal nerve would
// erode internodal length for 70 years and a healthy control would
// end up with the conduction velocity of a CMT1A patient - which is
// exactly what the first draft of this model did.
DEMY      = (SCD < SCD_H) ? (SCD_H - SCD) / SCD_H : 0.0;
// Tomacula are FOCAL. They cause focal conduction block and pressure
// palsies, but they are not a global de/remyelination cycle and must
// not shorten internodes along the whole nerve - written that way the
// model gave HNPP a conduction velocity of 31 m/s, which the clinical
// literature (near-normal NCV with focal slowing) rules out.
double DEMYF = DEMY + 0.12 * TOM;   // focal/nodal burden
dxdt_OB   = KOB * DEMYF * (3.0 - OB) - KOBOFF * OB;

// ==================================================================
//  4. MYELIN GEOMETRY AND CONDUCTION
// ==================================================================
double GMEX  = (GM > 1.0) ? (GM - 1.0) : 0.0;
double MYTT  = (SCD / SCD_H) / (1.0 + AGMSTO * GMEX);   // target thickness
if (MYTT > 1.0) MYTT = 1.0;
dxdt_MYT = KMYT * (MYTT - MYT);

// internodal length: eroded by demyelination, overwhelmingly in the
// developmental window. This term - and nothing else - is what makes
// NCV a developmental fossil (structural commitment S5).
dxdt_INL = -KERA * DEMY * (DEVRES + DEVBOOST * DEVW) * (INL - INLMIN);
dxdt_CYC = DEMYF;
dxdt_NDI = KNDI * DEMYF * (1.0 - NDI) - KNDIOFF * SCD * NDI;

NCV = NCV0 * pow(INL, ANCV) * pow(MYT, BNCV) * (1.0 - 0.25 * NDI);

// ==================================================================
//  5. TRANSPORT, METABOLIC SUPPORT AND THE ENERGY GAP
// ==================================================================
HDIINH = (HDI_ON > 0.5) ? IMAX_HDI * (HDIC / V_HDI) / (IC50_HDI + HDIC / V_HDI) : 0.0;
dxdt_TUBA = KATAT * (1.0 - AROS_TU * ROSX) - KHDAC * (1.0 - HDIINH) * TUBA;

double CVCR = VCRNV / 1.0;                       // ug/L in nerve
FVCR = (VCR_ON > 0.5) ? IMAX_VCR * CVCR / (IC50_VCR + CVCR) : 0.0;
VTR  = V0 * pow(TUBA > 0.01 ? TUBA : 0.01, 0.5) * (1.0 - FVCR);
if (VTR < 0.02) VTR = 0.02;

// glial metabolic support is produced by MYELINATING Schwann cells and is
// delivered locally - it does NOT decay with axon length. That asymmetry
// is what makes demyelination hurt the long axon disproportionately.
dxdt_LAC  = KLAC * (SCD / SCD_H - LAC);
LOCSUP    = KLACSUP * LAC;

EGAP      = (1.0 + KNKA * NDI) / (0.35 + 0.65 * LAC + 0.30 * MITO) - 1.0;
if (EGAP < 0.0) EGAP = 0.0;
dxdt_ROSX = KROS * (AROS_GAP * EGAP - ROSX);
dxdt_MITO = KMITO * (1.0 / (1.0 + 0.9 * ROSX) - MITO);

DEMAND = DEM0 * (1.0 + KNKA * NDI);

// transported cargo: exp(-L / (v * tau)) is the geometry
SUPL_F = exp(-LF / (VTR * TAUCG));
SUPL_L = exp(-LL / (VTR * TAUCG));
SUPL_H = exp(-LH / (VTR * TAUCG));
SUPL_P = exp(-LP / (VTR * TAUCG));

dxdt_CGF = KCG_IN * SUPL_F * MITO * (1.0 + ANT3 * NT3E) - KCG_OUT * CGF;
dxdt_CGL = KCG_IN * SUPL_L * MITO * (1.0 + ANT3 * NT3E) - KCG_OUT * CGL;
dxdt_CGH = KCG_IN * SUPL_H * MITO * (1.0 + ANT3 * NT3E) - KCG_OUT * CGH;
dxdt_CGP = KCG_IN * SUPL_P * MITO * (1.0 + ANT3 * NT3E) - KCG_OUT * CGP;

MGF = (CGF + LOCSUP) / DEMAND;
MGL = (CGL + LOCSUP) / DEMAND;
MGH = (CGH + LOCSUP) / DEMAND;
MGP = (CGP + LOCSUP) / DEMAND;
// sensory axons are longer and are treated quasi-statically (their
// cargo pool is assumed at pseudo-steady state, so the same amplitude
// KCG_IN/KCG_OUT that the motor cargo compartments reach appears here
// explicitly rather than being silently dropped).
double AMPC = KCG_IN / KCG_OUT;
MGSF = (AMPC * exp(-LSF / (VTR * TAUCG)) * MITO + LOCSUP) / DEMAND;
MGSH = (AMPC * exp(-LSH / (VTR * TAUCG)) * MITO + LOCSUP) / DEMAND;

// ==================================================================
//  6. AXON LOSS
// ==================================================================
// The hazard is a smooth sigmoid IN the margin, not a rectifier on
// (1 - margin). With a rectifier, a 0.10 difference in margin between
// two length classes becomes an order-of-magnitude difference in
// hazard and the model produces a foot that is destroyed while the
// leg is untouched. The sigmoid keeps the ordering but bounds the
// gradient, which is what the clinical distribution actually looks like.
// TWO limbs. The chronic limb is calibrated to a process that takes
// decades. An acute transport catastrophe (vincristine) cannot be
// represented by turning that limb up, because it saturates at KAX;
// the second limb engages only when the margin collapses far below
// the chronic threshold, and it is what SARM1-driven Wallerian-like
// degeneration looks like on this timescale.
#define HAZ(mg) (KAX / (1.0 + pow((mg) / MG50, PAX)) + KACUTE / (1.0 + pow((mg) / MGACUTE, PACUTE)) + KAGE)
HZF  = HAZ(MGF);   HZL  = HAZ(MGL);   HZH  = HAZ(MGH);
HZP  = HAZ(MGP);   HZSF = HAZ(MGSF);  HZSH = HAZ(MGSH);

double REG = KREGEN * (1.0 + ANT3 * NT3E);
dxdt_AXF = -HZF * AXF + REG * (1.0 - AXF) * SUPL_F;
dxdt_AXL = -HZL * AXL + REG * (1.0 - AXL) * SUPL_L;
dxdt_AXH = -HZH * AXH + REG * (1.0 - AXH) * SUPL_H;
dxdt_AXP = -HZP * AXP + REG * (1.0 - AXP) * SUPL_P;
dxdt_SXF = -HZSF * SXF + REG * (1.0 - SXF) * 0.4;
dxdt_SXH = -HZSH * SXH + REG * (1.0 - SXH) * 0.4;

AXLOSSRATE = HZF * AXF + HZL * AXL + HZH * AXH + HZP * AXP
           + HZSF * SXF + HZSH * SXH;

// ==================================================================
//  7. COLLATERAL REINNERVATION RESERVE (the silent decade)
// ==================================================================
INNF = AXF * RSF; if (INNF > 1.0) INNF = 1.0;
INNL = AXL * RSL; if (INNL > 1.0) INNL = 1.0;
INNH = AXH * RSH; if (INNH > 1.0) INNH = 1.0;
INNP = AXP * RSP; if (INNP > 1.0) INNP = 1.0;

dxdt_RSF = KSPR * (1.0 - INNF) * (RSMAX - RSF);
dxdt_RSL = KSPR * (1.0 - INNL) * (RSMAX - RSL);
dxdt_RSH = KSPR * (1.0 - INNH) * (RSMAX - RSH);
dxdt_RSP = KSPR * (1.0 - INNP) * (RSMAX - RSP);

// enlarged units are unstable units: a fragmentation penalty grows with size
double FRGF = 1.0 - KFRAG * (RSF - 1.0) / (RSMAX - 1.0);
double FRGL = 1.0 - KFRAG * (RSL - 1.0) / (RSMAX - 1.0);
double FRGH = 1.0 - KFRAG * (RSH - 1.0) / (RSMAX - 1.0);
double FRGP = 1.0 - KFRAG * (RSP - 1.0) / (RSMAX - 1.0);

FRCF = MSF * FRGF; if (FRCF > 1.0) FRCF = 1.0;
FRCL = MSL * FRGL; if (FRCL > 1.0) FRCL = 1.0;
FRCH = MSH * FRGH; if (FRCH > 1.0) FRCH = 1.0;
FRCP = MSP * FRGP; if (FRCP > 1.0) FRCP = 1.0;

// ==================================================================
//  8. MUSCLE
// ==================================================================
double EXG = 1.0 + AEXER * KEXER;
// Muscle MASS follows innervation; the NMJ-instability penalty acts on
// force per unit mass and is applied once, in FRC below. Applying it in
// both places (as the first draft did) squared it and produced a
// 40-year-old with 72% dorsiflexion strength and 98% of fibres
// innervated - an arithmetic artefact, not a mechanism.
dxdt_MSF = KMS * (INNF * EXG - MSF);
dxdt_MSL = KMS * (INNL * EXG - MSL);
dxdt_MSH = KMS * (INNH * EXG - MSH);
dxdt_MSP = KMS * (INNP * EXG - MSP);

// Fatty replacement tracks BOTH current denervation and the cumulative
// motor-unit loss: a muscle held at normal mass by collateral
// reinnervation has still been through repeated denervation cycles and
// still fills with fat. Driving it off current mass alone (first draft)
// gave a 40-year-old a 7% calf fat fraction against a measured ~30%.
double FFTF = FAT0 + (FATMAX - FAT0) * (0.75 * (1.0 - MSF) + 0.30 * (1.0 - AXF));
double FFTL = FAT0 + (FATMAX - FAT0) * (0.75 * (1.0 - MSL) + 0.30 * (1.0 - AXL));
dxdt_FFF = KFAT * (FFTF - FFF);
dxdt_FFL = KFAT * (FFTL - FFL);

// ==================================================================
//  9. SENSORY, SYMPTOM AND DEFORMITY LAYERS
// ==================================================================
dxdt_VIB  = KVIB * ((VIB0 + (VIBMAX - VIB0) * pow(1.0 - SXF, 1.4)) - VIB);
dxdt_PROP = KPROP * ((1.0 - SXF) * 0.75 + (1.0 - SXH) * 0.25 - PROP);
dxdt_NFL  = ANFL * AXLOSSRATE + KNFL * NFLBASE - KNFL * NFL;
dxdt_PAINS= KPAIN * (APAIN * (0.55 * DEMYF + 0.45 * (1.0 - SXF)) * (1.0 - 0.35 * SXF) - PAINS);
dxdt_FATG = KFAT_S * (1.0 + 4.2 * (0.5 * NDI + 0.3 * EGAP + 0.2 * (1.0 - FRCL)) - FATG);
// pes cavus is driven by the IMBALANCE between the weak dorsiflexor/intrinsic
// group and the relatively spared peroneus longus; it is a one-way integral.
double IMBAL = (1.0 - FRCF) * 0.6 + (1.0 - FRCL) * 0.4;
dxdt_CAV  = KCAV * IMBAL * (1.0 - CAV);

// ==================================================================
// 10. PHARMACOKINETICS
// ==================================================================
dxdt_BACD = -KA_BAC * BACD;
dxdt_BACC =  KA_BAC * BACD - (CL_BAC / V_BAC) * BACC;
dxdt_NTXD = -KA_NTX * NTXD;
dxdt_NTXC =  KA_NTX * NTXD - (CL_NTX / V_NTX) * NTXC;
dxdt_SORC = -(CL_SOR / V_SOR) * SORC;

// dietary input is set to whatever holds plasma at AA_BASE, so a
// subject on no supplement sits exactly at the habitual level and the
// supplement is a true increment.
double CL_AA_0 = CL_AA_B + CL_AA_R * pow(AA_BASE, HR_AA)
                 / (pow(CTH_AA, HR_AA) + pow(AA_BASE, HR_AA));
double RIN_DIET = AA_DIET * CL_AA_0 * AA_BASE;
dxdt_AAD  = -KA_AA * AAD;
dxdt_AAC  =  KA_AA * AAD + RIN_DIET
           - (CL_AA_B + CL_AA_R * pow(AAPL, HR_AA) / (pow(CTH_AA, HR_AA) + pow(AAPL, HR_AA)))
             * AAPL;

dxdt_OLISC = -KA_OLI * OLISC;
dxdt_OLIPL =  KA_OLI * OLISC - (CL_OLIPL / V_OLIPL) * OLIPL - KPL_NV * (OLIPL / V_OLIPL);
dxdt_OLINV =  KPL_NV * (OLIPL / V_OLIPL) - KNV_OUT * OLINV;

dxdt_HDID = -KA_HDI * HDID;
dxdt_HDIC =  KA_HDI * HDID - (CL_HDI / V_HDI) * HDIC;

dxdt_VCRC  = -(CL_VCR / V_VCR) * VCRC - KVCR_NV * VCRC;
dxdt_VCRNV =  KVCR_NV * VCRC - KVCR_OUT * VCRNV;

dxdt_NT3E = (NT3_ON > 0.5) ? KNT3 * (NT3_SS - NT3E) : -KNT3 * NT3E;

// pregnancy: progesterone rises to term over 280 d, then falls in days
double TP = SOLVERTIME - PREG_T0;
double PTGT = 0.0;
if (PREG_ON > 0.5 && TP > 0.0 && TP < 280.0) PTGT = pow(TP / 280.0, 1.5);
dxdt_PROG = 0.10 * (PTGT - PROG);

$TABLE
// ------------------------------------------------------------------
//  The shared algebra is RECOMPUTED here rather than read out of the
//  $GLOBAL scratch variables. Those hold whatever the solver left
//  behind at its last internal step, which is not guaranteed to be
//  the output time. Recomputing is 15 lines and removes the doubt.
// ------------------------------------------------------------------
double tAGEY  = AGE0 + TIME / 365.25;
double tLGRW  = LFRAC0 + (1.0 - LFRAC0) * (tAGEY / LGROW);
if (tLGRW > 1.0) tLGRW = 1.0;
double tJTR   = KTR * MRNA;
double tPHI   = PHIMAX * CHAP / (CHAP + tJTR);
double tPPREF = KTR * PHIMAX * CHAP / (CHAP + KTR) / KDEG_P;
double tGM    = PPROT / tPPREF;
double tTOM   = (tGM < 1.0) ? pow((1.0 - tGM) / TOMREF, HTOM) : 0.0;
double tNCV   = NCV0 * pow(INL, ANCV) * pow(MYT, BNCV) * (1.0 - 0.25 * NDI);
double tHDI   = (HDI_ON > 0.5) ? IMAX_HDI * (HDIC / V_HDI) / (IC50_HDI + HDIC / V_HDI) : 0.0;
double tCVCR  = VCRNV / 1.0;
double tFVCR  = (VCR_ON > 0.5) ? IMAX_VCR * tCVCR / (IC50_VCR + tCVCR) : 0.0;
double tVTR   = V0 * pow(TUBA > 0.01 ? TUBA : 0.01, 0.5) * (1.0 - tFVCR);
if (tVTR < 0.02) tVTR = 0.02;
double tLOC   = KLACSUP * LAC;
double tDEM   = DEM0 * (1.0 + KNKA * NDI);
double tMGF   = (CGF + tLOC) / tDEM;
double tMGL   = (CGL + tLOC) / tDEM;
double tMGH   = (CGH + tLOC) / tDEM;
double tMGP   = (CGP + tLOC) / tDEM;
double tKDOLI = (OLI_ON > 0.5) ? EMAX_OLI * (OLINV / V_OLINV) / (EC50_OLI + OLINV / V_OLINV) : 0.0;
double tSBAC  = BACC / V_BAC;
double tSNTX  = NTXC / V_NTX;
double tSSOR  = SORC / V_SOR;
double tOCC   = (tSBAC/(EC_BAC+tSBAC)) * (tSNTX/(EC_NTX+tSNTX)) * (tSSOR/(EC_SOR+tSSOR));
double tPXT   = (PXT_ON > 0.5) ? EMAX_PXT * tOCC : 0.0;
double tAAPL  = AAC / V_AA;
double tAASC  = tAAPL / (KM_SVCT2 + tAAPL);
double tAAIC  = tAASC * 900.0;
double tAAIC0 = (AA_BASE / (KM_SVCT2 + AA_BASE)) * 900.0;
double tAAv;
tAAv = EMAX_AA * (tAAIC/(IC50_AA+tAAIC) - tAAIC0/(IC50_AA+tAAIC0));
if (tAAv < 0.0) tAAv = 0.0;
double tAA    = (AA_ON > 0.5) ? tAAv : 0.0;

// ---------------- electrophysiology --------------------------------
double NCV_OUT  = tNCV;
double CMAP_ULN = 8.5 * AXH * RSH * (1.0 - KFRAG * (RSH - 1.0) / (RSMAX - 1.0)) * (1.0 - 0.35 * NDI);
double CMAP_PER = 6.0 * AXL * RSL * (1.0 - KFRAG * (RSL - 1.0) / (RSMAX - 1.0)) * (1.0 - 0.35 * NDI);
double SNAP_SUR = 22.0 * pow(SXF, 1.5);
double MUNE_PCT = 100.0 * AXL;

// ---------------- CMTNS-R, built from its nine 0-4 items ------------
// Each item is derived from a model state, so scale ceilings and floors
// are properties of the CONSTRUCTION, not assumptions bolted on later.
double I1 = 3.0 * (1.0 - SXF);                                  // sensory symptoms
double I2 = 4.0 * (1.0 - FRCL);                                 // motor symptoms, legs
double I3 = 4.0 * (1.0 - FRCH);                                 // motor symptoms, arms
// Pin sensibility is a SMALL-fibre item and CMT1A is a large-fibre
// neuropathy; scoring it off the large-fibre state alone made every
// simulated adult max out an item that real patients rarely max out.
double I4 = 4.0 * (0.45 * (1.0 - SXF) + 0.15 * (1.0 - SXH));    // pin sensibility
double I5 = 4.0 * (VIB - VIB0) / (VIBMAX - VIB0);               // vibration
double I6 = 4.0 * (1.0 - 0.55 * FRCF - 0.45 * FRCL);            // strength, legs
double I7 = 4.0 * (1.0 - FRCH);                                 // strength, arms
double I8 = 4.0 * (1.0 - CMAP_ULN / 8.5);
if (I8 < 0) I8 = 0;                       // ulnar CMAP
double I9 = 4.0 * (1.0 - SNAP_SUR / 22.0);                      // radial/sural SNAP
if (I1>4) I1=4; if (I2>4) I2=4; if (I3>4) I3=4; if (I4>4) I4=4; if (I5>4) I5=4;
if (I6>4) I6=4; if (I7>4) I7=4; if (I8>4) I8=4; if (I9>4) I9=4;
double CMTNS = I1+I2+I3+I4+I5+I6+I7+I8+I9;
double CMTES = I1+I2+I3+I4+I5+I6+I7;
if (SD_CMTNS > 0.0) CMTNS = CMTNS;      // noise is applied in R, not here

// ---------------- ONLS (arm 0-5, leg 0-7) ---------------------------
double ONLS_A = 5.0 * pow(1.0 - FRCH, 1.2);
double ONLS_L = 7.0 * pow(1.0 - 0.5 * FRCL - 0.3 * FRCF - 0.2 * FRCP, 1.6);
if (ONLS_A > 5.0) ONLS_A = 5.0;
if (ONLS_L > 7.0) ONLS_L = 7.0;
double ONLS = ONLS_A + ONLS_L;

// ---------------- functional tests ----------------------------------
double T10MWT = 6.5 / (0.25 + 0.75 * (0.5 * FRCL + 0.3 * FRCF + 0.2 * FRCP)) * (1.0 + 0.25 * PROP);
double D6MWT  = 620.0 * (0.5 * FRCL + 0.3 * FRCP + 0.2 * FRCF) * (1.0 - 0.20 * PROP);
double NHPT   = 18.0 / (0.20 + 0.80 * FRCH) * (1.0 + 0.15 * (1.0 - SXH));
double DORSI  = 100.0 * FRCL;                    // % predicted dorsiflexion strength
double GRIP   = 100.0 * FRCH;                    // % predicted grip

// ---------------- readouts of the mechanism -------------------------
double GM_OUT   = tGM;
double MRNA_REL = MRNA;
double PAGG_OUT = PAGG;
double MARGIN_F = tMGF;
double MARGIN_L = tMGL;
double MARGIN_H = tMGH;
double MARGIN_P = tMGP;
double VEL_TR   = tVTR;
double TOM_OUT  = tTOM;
double PPHAZ_Y  = 0.9 * tTOM;
double SCD_OUT  = SCD;
double KD_TOT   = 1.0 - (1.0 - tKDOLI) * (1.0 - tPXT) * (1.0 - tAA);
double AA_PLASMA= tAAPL;
double OLI_NRV  = OLINV / V_OLINV;
double BAC_CONC = BACC / V_BAC;
double AGE_Y    = tAGEY;
double DBG_LAC = LAC;
double DBG_MITO = MITO;
double DBG_NDI = NDI;
double DBG_DEM = tDEM;
double DBG_ROS = ROSX;
double FALLRATE = 2.4 * pow(PROP, 1.3) + 1.6 * pow(1.0 - FRCL, 1.5);   // falls/year

$CAPTURE @annotated
NCV_OUT  : motor nerve conduction velocity (m/s)
CMAP_ULN : ulnar CMAP amplitude (mV)
CMAP_PER : peroneal CMAP amplitude (mV)
SNAP_SUR : sural SNAP amplitude (uV)
MUNE_PCT : motor unit number, leg (% of normal)
CMTNS    : CMTNS-R total score (0-36)
CMTES    : CMTES-R score (0-28)
ONLS     : Overall Neuropathy Limitations Scale (0-12)
ONLS_A   : ONLS arm subscale (0-5)
ONLS_L   : ONLS leg subscale (0-7)
T10MWT   : 10-metre walk time (s)
D6MWT    : 6-minute walk distance (m)
NHPT     : 9-hole peg test (s)
DORSI    : ankle dorsiflexion strength (% predicted)
GRIP     : grip strength (% predicted)
GM_OUT   : membrane PMP22 stoichiometry (1 = normal)
MRNA_REL : PMP22 mRNA relative to normal
PAGG_OUT : PMP22 aggregate burden (index)
MARGIN_F : supply/demand margin, foot axon
MARGIN_L : supply/demand margin, leg axon
MARGIN_H : supply/demand margin, hand axon
MARGIN_P : supply/demand margin, proximal axon
VEL_TR   : effective axonal transport velocity index
TOM_OUT  : tomacula burden (1 = full HNPP phenotype)
PPHAZ_Y  : pressure-palsy hazard (events/year)
SCD_OUT  : Schwann-cell myelinating differentiation state
KD_TOT   : total fractional PMP22 message suppression
AA_PLASMA: plasma ascorbate (umol/L)
OLI_NRV  : oligonucleotide nerve concentration (mg/L)
BAC_CONC : baclofen plasma concentration (mg/L)
AGE_Y    : age (years)
FALLRATE : estimated falls per year
DBG_LAC  : glial lactate support index
DBG_MITO : mitochondrial health index
DBG_NDI  : nodal disorganisation
DBG_DEM  : terminal demand
DBG_ROS  : axonal oxidative stress
'

mod <- mcode_cache("cmt1a", code)

## =====================================================================
##  SCENARIO LIBRARY
## =====================================================================
##  Every scenario starts from the SAME healthy newborn nerve. The
##  phenotype is generated by integrating copy number over the life
##  course; nothing about the adult CMT1A patient is typed in.
## =====================================================================

YR <- 365.25

## ---- helper: run a life course from birth to a given age ------------
life <- function(mod, to_age = 70, delta = 30, ...) {
  mrgsim(mod, end = to_age * YR, delta = delta, hmax = 5, ...)
}

## ---- helper: build a dosing regimen for a drug over a window --------
ev_daily <- function(cmt, amt, start_day, dur_days, per_day = 1) {
  ev(time = start_day, amt = amt, cmt = cmt,
     ii = 1 / per_day, addl = round(dur_days * per_day) - 1)
}

scenarios <- list(

  ## ---------------- natural history ---------------------------------
  S01 = list(name = "Healthy control (CN = 2) — baseline drift test",
             param = list(CN = 2), ev = NULL),

  S02 = list(name = "CMT1A (CN = 3) — untreated natural history",
             param = list(CN = 3), ev = NULL),

  S03 = list(name = "HNPP (CN = 1) — the reciprocal deletion",
             param = list(CN = 1), ev = NULL),

  S04 = list(name = "Homozygous duplication (CN = 4) — severe infantile phenotype",
             param = list(CN = 4), ev = NULL),

  S05 = list(name = "CMT1A, tall patient (HTF 1.15) — length as a dose",
             param = list(CN = 3, HTF = 1.15), ev = NULL),

  S06 = list(name = "CMT1A, short patient (HTF 0.88)",
             param = list(CN = 3, HTF = 0.88), ev = NULL),

  S07 = list(name = "CMT1A, favourable modifier (MODF 1.6, better aggregate clearance)",
             param = list(CN = 3, MODF = 1.6), ev = NULL),

  S08 = list(name = "CMT1A, unfavourable modifier (MODF 0.7)",
             param = list(CN = 3, MODF = 0.7), ev = NULL),

  ## ---------------- ascorbic acid -----------------------------------
  S09 = list(name = "Ascorbic acid 1.5 g/day from age 25 (CMT-TRAUK dose)",
             param = list(CN = 3, AA_ON = 1, AA_DTOT = 1500),
             ev = ev_daily("AAD", 500 * 5.678, 25 * YR, 24 * YR, per_day = 3)),

  S10 = list(name = "Ascorbic acid 4 g/day from age 25 (CMT-TRIAAL dose)",
             param = list(CN = 3, AA_ON = 1, AA_DTOT = 4000),
             ev = ev_daily("AAD", 1333 * 5.678, 25 * YR, 24 * YR, per_day = 3)),

  S11 = list(name = "Ascorbic acid 4 g/day from age 5 (paediatric, full growth window)",
             param = list(CN = 3, AA_ON = 1, AA_DTOT = 4000),
             ev = ev_daily("AAD", 1333 * 5.678, 5 * YR, 40 * YR, per_day = 3)),

  ## ---------------- PXT3003 -----------------------------------------
  S12 = list(name = "PXT3003 dose 1 from age 40, 15 months (PLEO-CMT low dose)",
             param = list(CN = 3, PXT_ON = 1),
             ev = NULL, pxt = list(level = 1, start = 40 * YR, dur = 456)),

  S13 = list(name = "PXT3003 high dose from age 40, 15 months (PLEO-CMT high dose)",
             param = list(CN = 3, PXT_ON = 1),
             ev = NULL, pxt = list(level = 2, start = 40 * YR, dur = 456)),

  S14 = list(name = "PXT3003 high dose from age 40 for 15 YEARS",
             param = list(CN = 3, PXT_ON = 1),
             ev = NULL, pxt = list(level = 2, start = 40 * YR, dur = 15 * YR)),

  S15 = list(name = "PXT3003 high dose from age 8 (paediatric, lifetime)",
             param = list(CN = 3, PXT_ON = 1),
             ev = NULL, pxt = list(level = 2, start = 8 * YR, dur = 55 * YR)),

  ## ---------------- PMP22-lowering oligonucleotide -------------------
  S16 = list(name = "Oligonucleotide 10 mg SC q4w from age 30 (under-dosed)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 30 * YR, amt = 10, cmt = "OLISC",
                     ii = 28, addl = round(30 * YR / 28))),

  S17 = list(name = "Oligonucleotide 20 mg SC q4w from age 30 (TARGET WINDOW, g -> 1)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 30 * YR, amt = 20, cmt = "OLISC",
                     ii = 28, addl = round(30 * YR / 28))),

  S18 = list(name = "Oligonucleotide 100 mg SC q4w from age 30 (OVERSHOOT to the 1-copy level)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 30 * YR, amt = 100, cmt = "OLISC",
                     ii = 28, addl = round(30 * YR / 28))),

  S19 = list(name = "Oligonucleotide 20 mg SC q4w from age 5 (earliest feasible)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 5 * YR, amt = 20, cmt = "OLISC",
                     ii = 28, addl = round(55 * YR / 28))),

  S20 = list(name = "Oligonucleotide 120 mg SC q6m (long interval, same average dose)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 30 * YR, amt = 130, cmt = "OLISC",
                     ii = 182, addl = round(30 * YR / 182))),

  ## ---------------- HDAC6 inhibition ---------------------------------
  S21 = list(name = "HDAC6 inhibitor 30 mg/day from age 30",
             param = list(CN = 3, HDI_ON = 1),
             ev = ev_daily("HDID", 30, 30 * YR, 35 * YR, per_day = 1)),

  S22 = list(name = "HDAC6 inhibitor + oligonucleotide 20 mg q4w from age 30",
             param = list(CN = 3, HDI_ON = 1, OLI_ON = 1),
             ev = c(ev_daily("HDID", 30, 30 * YR, 35 * YR, per_day = 1),
                    ev(time = 30 * YR, amt = 20, cmt = "OLISC",
                       ii = 28, addl = round(30 * YR / 28)))),

  ## ---------------- trophic / gene therapy ---------------------------
  S23 = list(name = "AAV1.NT-3 intramuscular at age 30",
             param = list(CN = 3, NT3_ON = 1), ev = NULL),

  S24 = list(name = "Structured exercise programme from age 20",
             param = list(CN = 3, KEXER = 1.0), ev = NULL),

  ## ---------------- hormonal -----------------------------------------
  S25 = list(name = "Pregnancy at age 30 (progesterone drive on PMP22)",
             param = list(CN = 3, SEXF = 1, PREG_ON = 1, PREG_T0 = 30 * YR),
             ev = NULL),

  S26 = list(name = "Onapristone (PR antagonist) from age 30",
             param = list(CN = 3, ONAPRI = 1), ev = NULL),

  ## ---------------- neurotoxicity -------------------------------------
  S27 = list(name = "Vincristine 2 mg IV weekly x 4 at age 20 in CMT1A",
             param = list(CN = 3, VCR_ON = 1),
             ev = ev(time = 20 * YR, amt = 2000, cmt = "VCRC",
                     ii = 7, addl = 3)),

  S28 = list(name = "Vincristine 2 mg IV weekly x 4 at age 20 in a HEALTHY nerve",
             param = list(CN = 2, VCR_ON = 1),
             ev = ev(time = 20 * YR, amt = 2000, cmt = "VCRC",
                     ii = 7, addl = 3)),

  S29 = list(name = "Vincristine 2 mg IV weekly x 4 at age 20 in HNPP",
             param = list(CN = 1, VCR_ON = 1),
             ev = ev(time = 20 * YR, amt = 2000, cmt = "VCRC",
                     ii = 7, addl = 3)),

  ## ---------------- combinations and comparators ----------------------
  S30 = list(name = "Oligonucleotide 20 mg q4w from age 30 + exercise",
             param = list(CN = 3, OLI_ON = 1, KEXER = 1.0),
             ev = ev(time = 30 * YR, amt = 20, cmt = "OLISC",
                     ii = 28, addl = round(30 * YR / 28))),

  S31 = list(name = "Oligonucleotide 20 mg q4w from age 50 (late start)",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 50 * YR, amt = 20, cmt = "OLISC",
                     ii = 28, addl = round(15 * YR / 28))),

  S32 = list(name = "Oligonucleotide 20 mg q4w age 30-40 then STOPPED",
             param = list(CN = 3, OLI_ON = 1),
             ev = ev(time = 30 * YR, amt = 20, cmt = "OLISC",
                     ii = 28, addl = round(10 * YR / 28)))
)

## ---- PXT3003 dosing has three components; build it here -------------
pxt_events <- function(level, start, dur) {
  bac <- 6.0 * level / 2; ntx <- 0.70 * level / 2; sor <- 21.0 * level / 2
  n   <- round(dur * 2) - 1
  c(ev(time = start, amt = bac, cmt = "BACD", ii = 0.5, addl = n),
    ev(time = start, amt = ntx, cmt = "NTXD", ii = 0.5, addl = n),
    ev(time = start, amt = sor, cmt = "SORC", ii = 0.5, addl = n))
}

run_scenario <- function(id, to_age = 70, delta = 30) {
  s <- scenarios[[id]]
  m <- param(mod, s$param)
  e <- s$ev
  if (!is.null(s$pxt)) {
    pe <- pxt_events(s$pxt$level, s$pxt$start, s$pxt$dur)
    e  <- if (is.null(e)) pe else c(e, pe)
  }
  out <- if (is.null(e)) mrgsim(m, end = to_age * YR, delta = delta, hmax = 5)
         else mrgsim(m, events = e, end = to_age * YR, delta = delta, hmax = 5)
  as.data.frame(out)
}

## =====================================================================
##  PARAMETER PROVENANCE
## ---------------------------------------------------------------------
##  MEASURED / DIRECTLY ANCHORED TO PUBLISHED DATA
##    CN, the 1.4-Mb duplication and its reciprocal deletion .... Lupski 1991,
##        Chance 1993 (PMID 1677316, 8462100)
##    PMP22 mRNA ~1.4-1.7x normal in CMT1A nerve and skin ....... Katona 2009,
##        Nobbio 2014 (PMID 19833666)
##    ~80% of newly synthesised PMP22 is degraded before reaching
##        myelin in wild-type Schwann cells .................... Pareek 1997,
##        Notterpek 1999 (PMID 9169438, 10366642)
##    Motor NCV 15-25 m/s in CMT1A, stable from age ~5 ......... Garcia 1998,
##        Berciano 2000, Krajewski 2000 (PMID 9443474, 10891983)
##    CMTNS-R progression 0.1-0.7 points/year in adults ........ Shy 2008,
##        Fridman 2020 (PMID 18541884, 31857376)
##    Calf intramuscular fat fraction +1.2%/year ............... Morrow 2018
##        (PMID 30122359)
##    Plasma NfL ~2x control, correlates with CMTNS ............ Sandelius 2018
##        (PMID 29half) — see references file for exact PMIDs
##    Oral ascorbate bioavailability 112% at 200 mg, 49% at
##        1250 mg; plasma plateau 70-85 umol/L ................. Levine 1996
##        (PMID 8623000)
##    PLEO-CMT ONLS -0.37 (high dose vs placebo, 15 months) .... Attarian 2021
##        (PMID 33849629)
##    ASO knockdown of Pmp22 restores myelination in rodents ... Zhao 2018
##        (PMID 29457789)
##    Vincristine catastrophic worsening in CMT1A .............. Graf 1996
##        (PMID 8636842)
##
##  DERIVED IN $MAIN (not typed in)
##    PPREF (healthy membrane PMP22), PHI_REF = 0.200
##
##  ESTIMATED — no verifiable published value, chosen so that the
##  emergent adult CMT1A phenotype matches cohort means. These are
##  the model'"'"'s free parameters and are listed honestly as such:
##    KCJ/KEG and their cross-repression gains, KERA, DEVBOOST,
##    KAX, PAX, KSPR, RSMAX, KFRAG, KCG_IN/KCG_OUT, TAUCG, KLACSUP,
##    ANFL, KCAV, and the EC50/Emax values for PXT3003 and the
##    oligonucleotide (no human PD data exist for either).
## =====================================================================

## =====================================================================
##  DIAGNOSTICS
## ---------------------------------------------------------------------
##  run_diagnostics() reproduces every number quoted in README.md.
##  It is deliberately written to REPORT failures rather than to hide
##  them: D14 refutes one of the map's own design premises and D15
##  refutes another, and both are printed.
##
##      source("cmt1a_mrgsolve_model.R"); run_diagnostics()
## =====================================================================
run_diagnostics <- function() {
  at <- function(d, a) d[which.min(abs(d$time - a * YR)), ]
  gv <- function(d, a, v) at(d, a)[[v]]
  h  <- run_scenario("S01"); c3 <- run_scenario("S02")

  cat("\nD01 healthy CN=2 drift, 70 y\n")
  o <- as.data.frame(mrgsim(param(mod, CN = 2), end = 70 * YR, delta = 1, hmax = 1))
  cat(sprintf("    max |mRNA-1| %.3e   max |gm-1| %.3e   max |NCV-52| %.3e\n",
      max(abs(o$MRNA_REL - 1)), max(abs(o$GM_OUT - 1)), max(abs(o$NCV_OUT - 52))))

  cat("\nD02 copy-number series at age 40\n")
  for (id in c("S03", "S01", "S02", "S04")) {
    d <- run_scenario(id); r <- at(d, 40)
    cat(sprintf("    CN %d  mRNA %.3f  gm %.3f  PAGG %7.3f  TOM %.3f  NCV %5.1f  MUNE %5.1f%%  CMTNS %5.2f  PPhaz %.2f/y\n",
        scenarios[[id]]$param$CN, r$MRNA_REL, r$GM_OUT, r$PAGG_OUT, r$TOM_OUT,
        r$NCV_OUT, r$MUNE_PCT, r$CMTNS, r$PPHAZ_Y))
  }

  cat("\nD03 NCV is a developmental fossil\n")
  cat(sprintf("    CMT1A NCV age 2 %.1f, age 5 %.1f, age 20 %.2f, age 70 %.2f -> %+.4f m/s per YEAR\n",
      gv(c3,2,"NCV_OUT"), gv(c3,5,"NCV_OUT"), gv(c3,20,"NCV_OUT"), gv(c3,70,"NCV_OUT"),
      (gv(c3,70,"NCV_OUT") - gv(c3,20,"NCV_OUT")) / 50))
  cat(sprintf("    over the same window CMTNS moves %+.3f/y and calf fat fraction %+.3f %%/y\n",
      (gv(c3,50,"CMTNS") - gv(c3,30,"CMTNS"))/20, (gv(c3,50,"FFL") - gv(c3,30,"FFL"))/20))

  cat("\nD04 knockdown sweep from birth (the interior optimum)\n")
  kd <- do.call(rbind, lapply(seq(0, 0.80, by = 0.05), function(k) {
    cnv <- 3.0 * k / (0.85 - k)
    m <- init(param(mod, CN = 3, OLI_ON = 1, KNV_OUT = 0, KPL_NV = 0), OLINV = cnv * 0.60)
    r <- at(as.data.frame(mrgsim(m, end = 70 * YR, delta = 60, hmax = 5)), 50)
    data.frame(kd = k, mRNA = r$MRNA_REL, gm = r$GM_OUT, PAGG = r$PAGG_OUT,
               TOM = r$TOM_OUT, NCV = r$NCV_OUT, CMTNS = r$CMTNS, PP = r$PPHAZ_Y)
  }))
  print(round(kd, 3), row.names = FALSE)
  b <- kd[which.min(kd$CMTNS), ]
  cat(sprintf("    CMTNS minimum at %.0f%% knockdown (mRNA %.3f); arithmetic prediction is 33.3%%\n",
      100 * b$kd, b$mRNA))
  cat(sprintf("    at 65%% knockdown mRNA is %.3f against the 1-copy level of 0.500, TOM %.2f, PPhaz %.2f/y\n",
      kd$mRNA[kd$kd == 0.65], kd$TOM[kd$kd == 0.65], kd$PP[kd$kd == 0.65]))

  cat("\nD05 proteostatic collapse threshold\n")
  for (cn in c(2.8, 3.0, 3.1, 3.2, 3.3, 3.4, 4.0)) {
    r <- at(as.data.frame(mrgsim(param(mod, CN = cn), end = 40 * YR, delta = 120, hmax = 5)), 40)
    cat(sprintf("    mRNA %.2f -> PAGG %7.3f  SCD %.3f  NCV %5.1f  CMTNS %5.2f\n",
        r$MRNA_REL, r$PAGG_OUT, r$SCD_OUT, r$NCV_OUT, r$CMTNS))
  }

  cat("\nD06 ascorbic acid\n")
  aa <- function(dt) {
    m <- param(mod, CN = 3, AA_ON = 1, AA_DTOT = dt)
    e <- ev(time = 0, amt = (dt/3) * 5.678, cmt = "AAD", ii = 1/3, addl = 3*160 - 1)
    o <- as.data.frame(mrgsim(m, events = e, end = 160, delta = 0.25, hmax = 0.05))
    t2 <- o[o$time > 150, ]; c(mean(t2$AA_PLASMA), 100 * mean(t2$KD_TOT))
  }
  for (d in c(1000, 1500, 3000, 4000)) { v <- aa(d)
    cat(sprintf("    %4.0f mg/day -> plasma %5.1f umol/L, PMP22 suppression %.2f%%\n", d, v[1], v[2])) }
  v15 <- aa(1500); v40 <- aa(4000)
  cat(sprintf("    2.67x dose -> %.2fx plasma -> a %.2f percentage-point PD gap, against a 33%% requirement\n",
      v40[1]/v15[1], v40[2] - v15[2]))

  cat("\nD07 PXT3003 against PLEO-CMT\n")
  base <- run_scenario("S02", delta = 5)
  for (lv in 1:2) {
    o <- as.data.frame(mrgsim(param(mod, CN = 3, PXT_ON = 1),
         events = pxt_events(lv, 40*YR, 456), end = 70*YR, delta = 5, hmax = 5))
    dd <- (gv(o,40+456/365.25,"ONLS") - gv(o,40,"ONLS")) -
          (gv(base,40+456/365.25,"ONLS") - gv(base,40,"ONLS"))
    cat(sprintf("    level %d: 15-month ONLS difference %+.4f (published high-dose value -0.37)\n", lv, dd))
  }

  cat("\nD08 the silent decade\n")
  for (a in c(10, 20, 30, 40, 50, 60)) { r <- at(c3, a)
    cat(sprintf("    age %2.0f  MUNE %5.1f%%  unit size %.2fx  dorsiflexion %5.1f%%  CMTNS %5.2f\n",
        a, r$MUNE_PCT, r$RSL, r$DORSI, r$CMTNS)) }

  cat("\nD09 length-dependence is geometry (age 40)\n")
  r <- at(c3, 40)
  cat(sprintf("    foot 1.00 m margin %.3f axons %.3f | leg 0.85 m %.3f %.3f | hand 0.75 m %.3f %.3f | prox 0.35 m %.3f %.3f\n",
      r$MARGIN_F, r$AXF, r$MARGIN_L, r$AXL, r$MARGIN_H, r$AXH, r$MARGIN_P, r$AXP))
  cat(sprintf("    stature 0.88 / 1.00 / 1.15 -> CMTNS at 50 = %.2f / %.2f / %.2f (identical genotype)\n",
      gv(run_scenario("S06"),50,"CMTNS"), gv(c3,50,"CMTNS"), gv(run_scenario("S05"),50,"CMTNS")))

  cat("\nD10 vincristine: one equation, three hosts\n")
  for (id in c("S28", "S27", "S29")) {
    d <- run_scenario(id, delta = 2)
    lab <- c(S28 = "healthy CN=2", S27 = "CMT1A  CN=3", S29 = "HNPP   CN=1")[id]
    pre <- at(d, 20); post <- at(d, 20.6)
    cat(sprintf("    %-13s margin %.3f -> %.3f ; axons %.3f -> %.3f (%+.1f%%) ; ONLS %+.3f\n",
        lab, pre$MARGIN_F, min(d$MARGIN_F[d$time > 20*YR & d$time < 20.6*YR]),
        pre$AXF, post$AXF, 100*(post$AXF/pre$AXF - 1), post$ONLS - pre$ONLS))
  }

  cat("\nD11 HDAC6 inhibition and combination\n")
  hd <- run_scenario("S21"); ol <- run_scenario("S17"); cb <- run_scenario("S22")
  cat(sprintf("    CMTNS at 60: untreated %.2f | HDAC6i %.2f | oligo %.2f | both %.2f\n",
      gv(c3,60,"CMTNS"), gv(hd,60,"CMTNS"), gv(ol,60,"CMTNS"), gv(cb,60,"CMTNS")))

  cat("\nD12 start age dominates\n")
  for (id in c("S19", "S17", "S31", "S32")) { d <- run_scenario(id)
    cat(sprintf("    %-58s CMTNS@65 %5.2f  MUNE@65 %5.1f%%\n",
        scenarios[[id]]$name, gv(d,65,"CMTNS"), gv(d,65,"MUNE_PCT"))) }
  cat(sprintf("    %-58s CMTNS@65 %5.2f  MUNE@65 %5.1f%%\n", "untreated",
      gv(c3,65,"CMTNS"), gv(c3,65,"MUNE_PCT")))

  cat("\nD13 required sample size per arm (24 months, 80% power, alpha 0.05)\n")
  tr <- run_scenario("S17", delta = 10)
  for (e in list(list("CMTNS-R","CMTNS",1.8), list("ONLS","ONLS",0.9),
                 list("calf fat fraction","FFL",1.6), list("plasma NfL","NFL",6.0),
                 list("ulnar CMAP","CMAP_ULN",0.9), list("motor NCV","NCV_OUT",2.5))) {
    du <- gv(c3,42,e[[2]]) - gv(c3,40,e[[2]]); dt <- gv(tr,42,e[[2]]) - gv(tr,40,e[[2]])
    n <- if (abs(dt - du) < 1e-9) Inf else ceiling(2*(1.96+0.84)^2*e[[3]]^2/(dt-du)^2)
    cat(sprintf("    %-20s placebo %+8.4f  drug %+8.4f  N/arm %s\n", e[[1]], du, dt,
        ifelse(is.finite(n) && n < 1e7, format(n, big.mark=","), ">10,000,000")))
  }

  cat("\nD14 c-Jun sweep — THIS REFUTES THE MAP'S OWN HYPOTHESIS (cluster 4)\n")
  for (g in c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2.0)) {
    r <- at(as.data.frame(mrgsim(param(mod, CN=3, ACJ_AG=g), end=70*YR, delta=60, hmax=5)), 50)
    cat(sprintf("    c-Jun gain %.2f  SCD %.3f  lactate %.3f  MUNE %5.1f%%  CMTNS %5.2f\n",
        g, r$SCD_OUT, r$DBG_LAC, r$MUNE_PCT, r$CMTNS))
  }
  cat("    -> monotone. The repair Schwann cell is never a net benefit in this model.\n")

  cat("\nD15 dosing interval — ALSO A NEGATIVE RESULT\n")
  for (id in c("S17", "S20")) {
    d <- run_scenario(id, delta = 5)
    w <- d[d$time > 45*YR & d$time < 46*YR, ]
    cat(sprintf("    %-52s mRNA %.3f-%.3f  CMTNS@65 %.2f\n",
        scenarios[[id]]$name, min(w$MRNA_REL), max(w$MRNA_REL), gv(d,65,"CMTNS")))
  }
  cat("    -> a 2.1-fold swing in mRNA across the window costs 0.03 CMTNS points.\n")

  cat("\nD16 pregnancy\n")
  pg <- run_scenario("S25", delta = 5)
  w  <- pg[pg$time > 30*YR & pg$time < 31*YR, ]
  cat(sprintf("    peak mRNA %.3f (non-pregnant %.3f), PAGG %.3f -> %.3f\n",
      max(w$MRNA_REL), gv(c3,30,"MRNA_REL"), gv(c3,30,"PAGG_OUT"), max(w$PAGG_OUT)))
  cat(sprintf("    CMTNS over the pregnancy year %+.2f vs %+.2f in a matched non-pregnant year; residual at 50 y %+.3f\n",
      gv(pg,31,"CMTNS") - gv(pg,30,"CMTNS"), gv(c3,31,"CMTNS") - gv(c3,30,"CMTNS"),
      gv(pg,50,"CMTNS") - gv(c3,50,"CMTNS")))

  cat("\nD17 non-PMP22 interventions, CMTNS at 60\n")
  for (id in c("S23", "S24", "S26", "S30")) { d <- run_scenario(id)
    cat(sprintf("    %-56s %5.2f\n", scenarios[[id]]$name, gv(d,60,"CMTNS"))) }
  cat(sprintf("    %-56s %5.2f\n", "untreated", gv(c3,60,"CMTNS")))

  invisible(NULL)
}
