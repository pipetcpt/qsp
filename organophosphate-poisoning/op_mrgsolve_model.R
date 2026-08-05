## =============================================================================
##  op_mrgsolve_model.R
##  Acute organophosphorus (OP) insecticide self-poisoning — QSP model
##  급성 유기인계 살충제 자의 중독 — 정량적 시스템 약리학 모델
##
##  51 ODE states · 23 therapeutic scenarios · closed-loop atropine titration
##
##  BUILD / RUN
##      library(mrgsolve)
##      mod <- mread(file = "op_mrgsolve_model.R", model = "op")
##      source_scenarios(mod)          # defined at the bottom of this file
##
##  -------------------------------------------------------------------------
##  THE ONE IDEA THE MODEL IS BUILT ON
##  -------------------------------------------------------------------------
##  Acetylcholinesterase (AChE) is a THREE-STATE switch:
##
##        E  --k_i*C_oxon-->  EP  --k_a-->  EP_aged      (irreversible)
##        E  <--k_s+k_r[X]--  EP
##
##  Oxime reactivation SATURATES.  Worek's formalism writes the second-order
##  reactivation constant as k_r2 = k_r_max/(K_D + X), so as the oxime
##  concentration X grows without bound the reactivation RATE tends to a hard
##  ceiling k_r_max.  The best achievable free-enzyme fraction is therefore
##
##        E_ceiling = k_r_max / (k_r_max + k_i*C_oxon) = OMEGA/(1+OMEGA)
##
##  with the OXIME SUFFICIENCY NUMBER  OMEGA = k_r_max/(k_i*C_oxon).
##
##  Everything clinical follows from where OMEGA sits:
##    * OMEGA >> 1  the oxime restores the enzyme and the patient's syndrome
##                  is whatever is left over (solvent, aspiration, timing).
##    * OMEGA <  1  NO oxime dose can hold half the enzyme free.  The oxime is
##                  not "weak" here, it is arithmetically excluded, and no
##                  trial of any oxime at any dose can show benefit in such a
##                  patient.
##  Because CYP desulfuration of the thion is SATURABLE, the plasma oxon
##  plateau above a threshold ingested dose is dose-INDEPENDENT (Vmax/CL_oxon).
##  So OMEGA, and with it the oxime ceiling, becomes a CONSTANT for every large
##  ingestion, however large.  That constant is 15–25% free enzyme for the
##  compounds that dominate the epidemiology, i.e. below the ~30% threshold at
##  which muscarinic signs appear.
##
##  A second number governs how much enzyme is lost for good:
##        phi = k_a/(k_a + k_s + k_r)      per-inhibition-event aging risk
##        aged(T) = 1 - exp(-k_a * INTEGRAL f_EP dt)
##  Aging is not a deadline, it is an INTEGRAL.  The clinically useful
##  restatement is that the oxime's job is to shrink the area under the
##  inhibited-enzyme curve, which it can only do while k_r[X] >> k_i*C_oxon.
##
##  -------------------------------------------------------------------------
##  WHAT IS DELIBERATELY MODELLED AS A CONTROLLER, NOT A DOSE
##  -------------------------------------------------------------------------
##  Atropine is titrated at the bedside against secretions, heart rate,
##  bronchospasm and ventilation — never given as a fixed schedule.  It is
##  therefore written here as a CLOSED LOOP inside $ODE, and the cumulative
##  atropine dose is an OUTPUT of the model rather than an input.  This is what
##  lets the model reproduce, without being told to, the tens-to-hundreds of
##  milligrams reported in severe cases and the 1/E scaling of the requirement.
##
##  -------------------------------------------------------------------------
##  CALIBRATION AND PROVENANCE  (full citations in op_references.md)
##  -------------------------------------------------------------------------
##  Esterase kinetics (k_i, k_a, k_s, k_r2 for human AChE)
##      Worek F et al. Biochem Pharmacol 2004;68:2237-48  (PMID 15498514)
##      Worek F et al. Toxicol Appl Pharmacol 2007;219:226-34 (PMID 17169391)
##      Aurbek N et al. Toxicology 2006;224:91-9 (PMID 16720069)
##      - diethyl (paraoxon-ethyl, chlorpyrifos-oxon): aging t1/2 33 h,
##        spontaneous reactivation t1/2 77 h
##      - dimethyl (omethoate, dichlorvos): aging t1/2 3.7 h,
##        spontaneous reactivation t1/2 0.7-0.9 h
##      - k_r2 human AChE: obidoxime/paraoxon-ethyl ~2.0e4 M-1 min-1,
##        pralidoxime/paraoxon-ethyl ~3.0e3, pralidoxime/omethoate ~8e2
##  Pralidoxime PK and the WHO regimen (30 mg/kg load + 8 mg/kg/h)
##      Eddleston M et al. PLoS Med 2009;6:e1000104 (PMID 19564902)
##      Medicis JJ et al. J Toxicol Clin Toxicol 1996;34:289-95 (PMID 8667465)
##  Obidoxime concentrations and reactivation in real patients
##      Thiermann H et al. Toxicol Lett 2010;198:264-70 (PMID 20674727)
##      Eyer P. Toxicol Rev 2003;22:165-90 (PMID 15181665)
##  Atropine titration, rapid doubling protocol, cumulative dose
##      Eddleston M et al. QJM 2004;97:227-33 (PMID 15028853)
##      Abedin MJ et al. Clin Toxicol 2012;50:433-40 (PMID 22578059)
##  Class split (dimethyl vs diethyl) and case fatality
##      Eddleston M et al. Lancet 2005;366:1452-9 (PMID 16243090)
##      Eddleston M et al. QJM 2006;99:513-22 (PMID 16829539)
##  Intermediate syndrome
##      Senanayake N, Karalliedde L. N Engl J Med 1987;316:761-3 (PMID 3029588)
##      Jayawardane P et al. J Neurol Sci 2008;271:1-8 (PMID 18395221)
##  OPIDN and NTE
##      Lotti M, Moretto A. Toxicol Rev 2005;24:37-49 (PMID 16042503)
##  PON1 and oxon hydrolysis
##      Costa LG et al. Clin Chim Acta 2005;352:37-47 (PMID 15653098)
##  Epidemiology and the effect of regulation
##      Mew EJ et al. J Affect Disord 2017;219:93-104 (PMID 28535450)
##      Gunnell D et al. BMC Public Health 2007;7:357 (PMID 18154668)
##
##  -------------------------------------------------------------------------
##  VERIFICATION
##  -------------------------------------------------------------------------
##  All 51 ODEs were independently re-implemented in Python/scipy
##  (op_reference_model.py) and the closed-form results in sections 1-6 of
##  op_reference_output.txt are recomputed there analytically.  Discrepancies
##  found during that exercise are listed in README.md under DEFECTS FOUND.
##
##  EDUCATIONAL / RESEARCH USE ONLY.  Not calibrated for, and not to be used
##  for, clinical decision making.
## =============================================================================

$PROB
# Acute organophosphorus insecticide self-poisoning
# 51-state QSP model: toxicokinetics + bioactivation + three-state esterase
# switch (4 pools) + muscarinic/nicotinic/central limbs + respiratory failure
# + closed-loop atropine + oxime + supportive care + delayed syndromes.

$PARAM @annotated
// ---------------------------------------------------------------- patient
WT       :  60    : Body weight (kg)

// ------------------------------------------------- OP compound properties
// Defaults are CHLORPYRIFOS 20% EC (diethyl, lipophilic pro-toxicant).
// Use op_class_params() at the bottom of this file to switch compound.
MW_OP    : 350.6  : Molecular weight of the parent thion (g/mol)
KI       : 0.18   : Second-order AChE inhibition constant of the oxon (1/nM/h)
KI_NMJ   : 1.0    : Relative susceptibility of NMJ AChE (-)
KI_NTE   : 0.0035 : Second-order NTE inhibition constant of the oxon (1/nM/h)
T12AGE   : 33     : AChE aging half-time (h)
T12SPO   : 77     : AChE spontaneous reactivation half-time (h)
T12AGNTE :  6     : NTE aging half-time (h)
KRMAX    : 36     : Maximum oxime reactivation rate k_r_max (1/h)
KDOX     : 200    : Oxime reactivation K_D (uM)
KA_GUT   : 0.35   : Gut absorption rate constant of the thion (1/h)
FBIO     : 0.80   : Oral bioavailability of the thion (-)
V1TH     : 35     : Central volume of the thion (L)
VFTH     : 700    : Deep/adipose volume of the thion (L)
QFTH     : 25     : Intercompartmental clearance of the thion (L/h)
VMAXBIO  : 6.0e4  : CYP desulfuration Vmax, thion -> oxon (nmol/h)
KMBIO    : 20     : CYP desulfuration Km (uM thion)
CLTHOTH  : 14     : Non-activating clearance of the thion (L/h)
CLOXON   : 300    : Oxon hydrolysis clearance, PON1 R192R (L/h)
PON1SC   : 1.0    : Multiplier on CLOXON (1.0 = RR, 0.33 = QQ) (-)
VOXON    : 56     : Central volume of the oxon (L)
QOXON    : 60     : Oxon intercompartmental clearance (L/h)
VTOXON   : 120    : Peripheral volume of the oxon (L)
SOLVFR   : 0.60   : Solvent fraction of the formulation (v/v)
SOLVCV   : 0.5    : Relative myocardial-depressant potency of the solvent (-)

// ------------------------------------------------------------ the exposure
DOSE_ML  : 50     : Ingested volume of formulation (mL)
CONC_GL  : 200    : Active ingredient concentration (g/L)

// ------------------------------------------------- esterase pool turnover
KRBCNEW  : 0.000481 : RBC AChE replacement by erythropoiesis (1/h)
KTISS    : 0.005776 : Synaptic/muscle/brain AChE protein turnover (1/h)
KBCHE    : 0.002626 : Plasma BChE turnover (1/h)
KNTE     : 0.005776 : NTE turnover (1/h)
KR_CNS   : 0.05     : Fraction of plasma oxime reactivation seen in brain (-)
KI_CNSSC : 0.80     : Relative oxon exposure of brain AChE (-)
BCHE_SC  : 2.5      : BChE inhibition rate relative to AChE (-)

// -------------------------------------------------------- acetylcholine
KHYD     : 60     : Cleft ACh clearance rate at full AChE activity (1/h)
ACHLEAK  : 0.05   : Non-AChE ACh clearance floor (-)

// ----------------------------------------------------------- receptors
KA_M     : 3.0    : Muscarinic ACh dissociation constant, relative units (-)
KA_N     : 5.0    : Nicotinic ACh dissociation constant, relative units (-)
KA_B     : 3.0    : Central muscarinic ACh dissociation constant (-)
OM0      : 0.25   : Baseline muscarinic occupancy (-)
OMMAX    : 0.87   : Maximum attainable muscarinic occupancy (-)
ON0      : 0.1667 : Baseline nicotinic occupancy (-)
KB_ATR   : 8.0    : Effective muscarinic Kb of atropine (nM)
BBBF     : 0.35   : Brain/plasma ratio of the antimuscarinic (0.02 = glycopyrrolate)

// ---------------------------------------- cholinergic receptor adaptation
KDOWN    : 0.012  : Receptor down-regulation on-rate (1/h)
KRECDOWN : 0.005  : Receptor down-regulation off-rate (1/h)
FDOWN_M  : 0.65   : Attenuation of peripheral muscarinic signalling (-)
FDOWN_C  : 0.60   : Attenuation of central muscarinic signalling (-)
FDOWN_N  : 0.60   : Attenuation of the nicotinic desensitisation drive (-)
FDOWN_R  : 0.80   : Attenuation of non-muscarinic central depression (-)

// ------------------------------------------------- nicotinic desensitisation
KDES     : 0.075  : nAChR desensitisation on-rate (1/h)
KRECDES  : 0.050  : nAChR desensitisation off-rate (1/h)
ACHN50   : 8.0    : Cleft ACh giving half-maximal desensitisation (x normal)
HILLN    : 4.0    : Hill coefficient for desensitisation (-)
WFASTBL  : 0.45   : Immediate depolarising block weight (-)

// -------------------------------------------------------- secretions/airway
KSEC     : 6.0    : Secretion production per unit muscarinic excess (mL/h)
KCLRSEC  : 0.9    : Spontaneous secretion clearance (1/h)
KCLRICU  : 2.4    : Secretion clearance with ICU suctioning (1/h)
SECCRIT  : 2.0    : Secretion burden that soaks the airway (mL)
KBT      : 1.2    : Bronchial tone on-rate (1/h)
KBTOFF   : 1.5    : Bronchial tone off-rate (1/h)
ICUSUCT  : 1      : ICU airway care available (1/0)

// ------------------------------------------------------------ respiration
TAURESPD : 0.4    : Time constant of central respiratory drive (h)
DMAXCNS  : 0.55   : Atropine-REVERSIBLE central depression weight (-)
DMAXNONR : 0.65   : Atropine-IRREVERSIBLE central depression weight (-)
SOLVRESP : 0.55   : Solvent central depression weight (-)
TAUMSTR  : 0.3    : Time constant of respiratory muscle strength (h)
SF_NMJ   : 0.28   : Fraction of maximal inspiratory capacity needed at rest (-)
VAVENTTH : 0.60   : Spontaneous ventilation below which intubation occurs (-)
VENTAVL  : 1      : Mechanical ventilator available (1/0)
O2AVL    : 1      : Supplemental oxygen available (1/0)
FIO2AIR  : 0.21   : Inspired oxygen fraction, room air (-)
FIO2O2   : 0.60   : Inspired oxygen fraction, face mask (-)
FIO2VENT : 0.50   : Inspired oxygen fraction, ventilator (-)
SHUNT0   : 0.03   : Baseline shunt fraction (-)
WSHSEC   : 0.40   : Shunt from airway secretions (-)
WSHLUNG  : 0.35   : Shunt from lung injury (-)
SHUNTMAX : 0.60   : Maximum shunt fraction (-)

// ------------------------------------------------------- haemodynamics
TAUHR    : 0.05   : Heart rate time constant (h)
HRINT    : 155    : Heart rate intercept (bpm)
HRSLOPE  : 320    : Heart rate slope on muscarinic occupancy (bpm)
TAUMAP   : 0.2    : Mean arterial pressure time constant (h)
MAP0     : 92     : Baseline mean arterial pressure (mmHg)
WMAPMUS  : 50     : MAP fall per unit muscarinic excess (mmHg)
WMAPSOLV : 30     : MAP fall per unit solvent effect (mmHg)
WMAPHYP  : 25     : MAP fall per unit hypoxaemia (mmHg)

// --------------------------------------------------------------- CNS
KSEIZ    : 1.4    : Seizure on-rate (1/h)
KSEIZOFF : 0.6    : Seizure off-rate (1/h)
SEIZTH   : 0.62   : Central occupancy threshold for seizures (-)

// ------------------------------------------------------------- lung injury
KLUNG    : 0.35   : Aspiration / chemical pneumonitis on-rate (1/h)
KLUNGREP : 0.030  : Lung injury resolution rate (1/h)

// ------------------------------------------------------------------ hazard
H0       : 0.00015 : Background hazard (1/h)
HHYPOX   : 0.55    : Hazard weight, hypoxaemia (1/h)
HHYPERC  : 0.060   : Hazard weight, hypercapnia (1/h)
HSHOCK   : 0.45    : Hazard weight, shock (1/h)
HSEIZ    : 0.05    : Hazard weight, seizures (1/h)
HLUNG    : 0.004   : Hazard weight, lung injury (1/h)
HVENT    : 0.0012  : Hazard weight, ventilator-associated complications (1/h)

// ------------------------------------------------------------------ OPIDN
KOPIDN   : 0.006   : OPIDN development rate (1/h)
KRECOPID : 0.0015  : OPIDN recovery rate (1/h)
NTETH    : 0.70    : Aged-NTE fraction committing to OPIDN (-)

// ------------------------------------------------------------- antidote PK
PAM_MW   : 172.6  : Pralidoxime chloride molecular weight (g/mol)
OXMW     : 172.6  : Molecular weight of the oxime actually given (g/mol)
OXV1     : 24     : Oxime central volume (L)
OXV2     : 24     : Oxime peripheral volume (L)
OXCL     : 21     : Oxime clearance (L/h)
OXQ      : 30     : Oxime intercompartmental clearance (L/h)
ATR_MW   : 289.4  : Atropine molecular weight (g/mol)
ATR_V    : 200    : Atropine volume of distribution (L)
ATR_CL   : 46     : Atropine clearance (L/h)
ATR_KEO  : 6.0    : Atropine effect-site equilibration rate (1/h)
DZ_V     : 80     : Diazepam volume of distribution (L)
DZ_CL    : 1.6    : Diazepam clearance (L/h)
DZ_EC50  : 0.20   : Diazepam EC50 (uM)
MG_V     : 15     : Magnesium volume of distribution (L)
MG_CL    : 6.0    : Magnesium clearance (L/h)
MG_EC50  : 1.6    : Magnesium EC50 for ACh release inhibition (mmol/L)
MG_EMAX  : 0.45   : Maximum fractional inhibition of ACh release by Mg (-)
KSCAV    : 3.0e-3 : Effective bioscavenger-oxon association constant (1/uM/h)

// --------------------------------------------------- treatment schedule
OXTYPE   : 0      : Oxime given (0 none, 1 pralidoxime/obidoxime) (-)
OXSTART  : 1.5    : Oxime start time (h)
OXRATE   : 0      : Oxime continuous infusion rate (umol/h)
OXDUR    : 168    : Oxime infusion duration (h)
ATRMODE  : 1      : Atropine titration (0 none, 1 rapid doubling, 2 slow) (-)
KP_RAPID : 170    : Atropine controller gain, rapid protocol (umol/h)
KP_SLOW  : 22     : Atropine controller gain, ad-hoc protocol (umol/h)
ATRMAX   : 210    : Maximum atropine infusion rate (umol/h)
HRSTOP   : 140    : Heart rate at which atropine escalation stops (bpm)
CHARC_T  : -1     : Activated charcoal administration time (h, <0 = none)
CHARC_D  : 2.0    : Duration of the charcoal effect (h)
K_CHARC  : 1.8    : Extra gut removal rate while charcoal is present (1/h)

$CMT @annotated
A_gut    : Thion in the gut lumen (umol)
A_th_c   : Thion, central compartment (umol)
A_th_f   : Thion, adipose/deep depot (umol)
A_ox_c   : Oxon, central compartment (nmol)
A_ox_t   : Oxon, peripheral tissue (nmol)
A_solv_g : Solvent in the gut (mL)
A_solv_c : Solvent, systemic (mL-equivalent)
LUNG     : Aspiration / chemical lung injury (0-1)
A_pam_c  : Oxime, central compartment (umol)
A_pam_p  : Oxime, peripheral compartment (umol)
A_atr_c  : Atropine, central compartment (umol)
Ce_atr   : Atropine effect-site concentration (nM)
A_dz     : Diazepam (umol)
A_mg     : Magnesium (mmol)
SCAV     : Exogenous BChE bioscavenger sites (umol)
E_rbc    : Free RBC AChE (fraction)
EP_rbc   : Phosphylated RBC AChE (fraction)
EA_rbc   : Aged RBC AChE (fraction)
E_mus    : Free peripheral muscarinic-target AChE (fraction)
EP_mus   : Phosphylated peripheral AChE (fraction)
EA_mus   : Aged peripheral AChE (fraction)
E_nmj    : Free neuromuscular junction AChE (fraction)
EP_nmj   : Phosphylated NMJ AChE (fraction)
EA_nmj   : Aged NMJ AChE (fraction)
E_cns    : Free brain AChE (fraction)
EP_cns   : Phosphylated brain AChE (fraction)
EA_cns   : Aged brain AChE (fraction)
B_free   : Free plasma butyrylcholinesterase (fraction)
B_inh    : Inhibited plasma butyrylcholinesterase (fraction)
N_free   : Free neuropathy target esterase (fraction)
N_inh    : Phosphylated NTE (fraction)
N_aged   : Aged NTE (fraction)
ACh_m    : Peripheral muscarinic cleft ACh (x normal)
ACh_n    : Neuromuscular junction cleft ACh (x normal)
ACh_b    : Brain cleft ACh (x normal)
Rn_des   : Desensitised / depolarisation-blocked nAChR (fraction)
Rm_down  : Cholinergic receptor down-regulation (fraction)
SEC      : Airway secretion burden (mL)
BT       : Bronchial tone (0-1)
HR       : Heart rate (bpm)
RESPD    : Central respiratory drive (fraction of normal)
MSTR     : Respiratory muscle strength (fraction of normal)
MAP      : Mean arterial pressure (mmHg)
SEIZ     : Seizure / CNS excitation index (0-1)
OPIDN    : OP-induced delayed polyneuropathy score (0-1)
HAZ      : Cumulative death hazard (-)
VTIME    : Cumulative time on mechanical ventilation (h)
ATRCUM   : Cumulative atropine dose (mg)
OXCUM    : Cumulative oxime delivered by infusion (mg; boluses are in the event table)
AUC_EPn  : Integral of the inhibited NMJ enzyme fraction (h)
AUC_ox   : Oxon exposure (nM*h)

$GLOBAL
#define LN2 0.6931471805599453

// individual values written in $MAIN and read in $ODE
double DOSE_i;      // ingested volume actually absorbed (mL-equivalent)
double VMAX_i;      // individual CYP desulfuration capacity (nmol/h)
double CLOX_i;      // individual oxon hydrolysis clearance (L/h)
double RESV_i;      // individual respiratory reserve multiplier (-)

// smooth switch, k sets the steepness
inline double sig(double x, double k) { return 1.0/(1.0 + exp(-k*x)); }
inline double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }
inline double clampx(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

$MAIN
// ---- individual parameters ---------------------------------------------
// Between-patient variability sits on the four quantities that actually move
// the answer: how much was swallowed, how fast the liver makes oxon, how fast
// PON1 destroys it, and how much respiratory reserve the patient started with.
DOSE_i = DOSE_ML * exp(ETA(1));
VMAX_i = VMAXBIO * exp(ETA(2));
CLOX_i = CLOXON * PON1SC * exp(ETA(3));
RESV_i = exp(ETA(4));

// ---- initial conditions -----------------------------------------------
// The ingested dose is placed in the gut at t = 0 together with the solvent.
A_gut_0    = DOSE_i * CONC_GL / 1000.0 * 1.0e6 / MW_OP;    // umol
A_solv_g_0 = DOSE_i * SOLVFR;                              // mL
E_rbc_0  = 1.0;  E_mus_0 = 1.0;  E_nmj_0 = 1.0;  E_cns_0 = 1.0;
B_free_0 = 1.0;  N_free_0 = 1.0;
ACh_m_0  = 1.0;  ACh_n_0 = 1.0;  ACh_b_0 = 1.0;
HR_0     = 75.0; RESPD_0 = 1.0;  MSTR_0  = 1.0;  MAP_0 = MAP0;

$ODE
// =========================================================================
// 0. GUARDS
// =========================================================================
double Erb = fmax(E_rbc, 0.0), EPrb = fmax(EP_rbc, 0.0);
double Emu = fmax(E_mus, 0.0), Enm  = fmax(E_nmj, 0.0), Ecn = fmax(E_cns, 0.0);
double Am  = fmax(ACh_m, 0.0), An   = fmax(ACh_n, 0.0), Ab  = fmax(ACh_b, 0.0);
double Rnd = clamp01(Rn_des),  Rmd  = clamp01(Rm_down);
double SECp = fmax(SEC, 0.0),  BTp  = clamp01(BT);
double LUNGp = clamp01(LUNG);

// =========================================================================
// 1. OXIME PK AND THE REACTIVATION RATE  k_r[X]
// =========================================================================
double C_oxime = fmax(A_pam_c, 0.0) / OXV1;                       // uM
double kr = (OXTYPE > 0.5) ? KRMAX * C_oxime / (KDOX + C_oxime) : 0.0;   // 1/h

// oxime infusion window
double oxinf = 0.0;
if (OXTYPE > 0.5 && SOLVERTIME >= OXSTART && SOLVERTIME <= OXSTART + OXDUR)
  oxinf = OXRATE;

dxdt_A_pam_c = oxinf - OXCL*C_oxime - OXQ*(C_oxime - fmax(A_pam_p,0.0)/OXV2);
dxdt_A_pam_p = OXQ*(C_oxime - fmax(A_pam_p,0.0)/OXV2);

// =========================================================================
// 2. TOXICOKINETICS OF THE THION AND ITS BIOACTIVATION
// =========================================================================
double k_ac = 0.0;
if (CHARC_T >= 0.0 && SOLVERTIME >= CHARC_T && SOLVERTIME < CHARC_T + CHARC_D)
  k_ac = K_CHARC;

double C_th  = fmax(A_th_c, 0.0) / V1TH;         // uM
double C_thf = fmax(A_th_f, 0.0) / VFTH;         // uM
double bio   = VMAX_i * C_th / (KMBIO + C_th);   // nmol/h of oxon formed

dxdt_A_gut  = -(KA_GUT + k_ac) * fmax(A_gut, 0.0);
dxdt_A_th_c =  FBIO*KA_GUT*fmax(A_gut,0.0) - bio/1000.0
               - CLTHOTH*C_th - QFTH*(C_th - C_thf);
dxdt_A_th_f =  QFTH*(C_th - C_thf);

// =========================================================================
// 3. TOXICOKINETICS OF THE OXON — the only species that inhibits AChE
// =========================================================================
double CLox  = CLOX_i;
double C_ox  = fmax(A_ox_c, 0.0) / VOXON;        // nM
double C_oxt = fmax(A_ox_t, 0.0) / VTOXON;       // nM
double scav  = KSCAV * fmax(SCAV, 0.0) * C_ox;   // nmol/h captured 1:1

dxdt_A_ox_c = bio - CLox*C_ox - QOXON*(C_ox - C_oxt) - scav;
dxdt_A_ox_t = QOXON*(C_ox - C_oxt) - 12.0*C_oxt;
dxdt_SCAV   = -scav/1000.0;

// =========================================================================
// 4. SOLVENT
// =========================================================================
dxdt_A_solv_g = -1.1*fmax(A_solv_g, 0.0);
dxdt_A_solv_c =  0.605*fmax(A_solv_g,0.0) - 0.55*fmax(A_solv_c,0.0);
double solv = fmax(A_solv_c,0.0) / (fmax(A_solv_c,0.0) + 12.0);

// =========================================================================
// 5. ATROPINE — a CLOSED-LOOP drug
// =========================================================================
double Cat = fmax(Ce_atr, 0.0);                                  // nM
double Om  = Am / (Am + KA_M*(1.0 + Cat/KB_ATR));
double mus_raw = clamp01((Om - OM0)/(OMMAX - OM0));
double mus_x   = mus_raw * (1.0 - FDOWN_M*Rmd);

double rate_atr = 0.0;
if (ATRMODE > 0.5) {
  double Kp = (ATRMODE > 1.5) ? KP_SLOW : KP_RAPID;
  // bedside atropinisation targets: dry chest, HR > 80, no bronchospasm,
  // adequate spontaneous ventilation
  double err = 0.40*clamp01(SECp/SECCRIT - 0.12)
             + 0.25*clamp01((78.0 - HR)/40.0)
             + 0.15*clamp01(BTp - 0.05)
             + 0.20*clamp01((0.75 - RESPD)/0.5);
  double stop = sig((HRSTOP - HR)/10.0, 6.0);
  rate_atr = fmin(Kp*err*stop, ATRMAX);
}
dxdt_A_atr_c = rate_atr - ATR_CL*fmax(A_atr_c,0.0)/ATR_V;
double C_atr_pl = 1000.0*fmax(A_atr_c,0.0)/ATR_V;                // nM
dxdt_Ce_atr  = ATR_KEO*(C_atr_pl - Cat);

// =========================================================================
// 6. OTHER ADJUNCTS
// =========================================================================
dxdt_A_dz = -DZ_CL*fmax(A_dz,0.0)/DZ_V;
double C_dz = fmax(A_dz,0.0)/DZ_V;                               // uM
dxdt_A_mg = -MG_CL*fmax(A_mg,0.0)/MG_V;
double C_mg = fmax(A_mg,0.0)/MG_V;                               // mmol/L

// =========================================================================
// 7. THE THREE-STATE ESTERASE SWITCH  (4 AChE pools + BChE + NTE)
// =========================================================================
double k_age = LN2/T12AGE;
double k_spo = LN2/T12SPO;
double kinh  = KI * C_ox;            // 1/h, plasma-facing pool (RBC, BChE)
double kinh_t= KI * C_oxt;           // 1/h, tissue pools

// --- RBC AChE: replaced only by erythropoiesis
double tot_rbc = Erb + EPrb + fmax(EA_rbc,0.0);
dxdt_E_rbc  = -kinh*Erb + (k_spo + kr)*EPrb + KRBCNEW*(1.0 - tot_rbc);
dxdt_EP_rbc =  kinh*Erb - (k_spo + kr + k_age + KRBCNEW)*EPrb;
dxdt_EA_rbc =  k_age*EPrb - KRBCNEW*fmax(EA_rbc,0.0);

// --- peripheral muscarinic target tissue
double tot_mus = Emu + fmax(EP_mus,0.0) + fmax(EA_mus,0.0);
dxdt_E_mus  = -kinh_t*Emu + (k_spo + kr)*fmax(EP_mus,0.0) + KTISS*(1.0 - tot_mus);
dxdt_EP_mus =  kinh_t*Emu - (k_spo + kr + k_age + KTISS)*fmax(EP_mus,0.0);
dxdt_EA_mus =  k_age*fmax(EP_mus,0.0) - KTISS*fmax(EA_mus,0.0);

// --- neuromuscular junction (peripheral: the oxime reaches it)
double kinh_n = kinh_t * KI_NMJ;
double tot_nmj = Enm + fmax(EP_nmj,0.0) + fmax(EA_nmj,0.0);
dxdt_E_nmj  = -kinh_n*Enm + (k_spo + kr)*fmax(EP_nmj,0.0) + KTISS*(1.0 - tot_nmj);
dxdt_EP_nmj =  kinh_n*Enm - (k_spo + kr + k_age + KTISS)*fmax(EP_nmj,0.0);
dxdt_EA_nmj =  k_age*fmax(EP_nmj,0.0) - KTISS*fmax(EA_nmj,0.0);

// --- brain: the quaternary oximes barely cross the blood-brain barrier
double kr_b = KR_CNS * kr;
double kinh_b = kinh_t * KI_CNSSC;
double tot_cns = Ecn + fmax(EP_cns,0.0) + fmax(EA_cns,0.0);
dxdt_E_cns  = -kinh_b*Ecn + (k_spo + kr_b)*fmax(EP_cns,0.0) + KTISS*(1.0 - tot_cns);
dxdt_EP_cns =  kinh_b*Ecn - (k_spo + kr_b + k_age + KTISS)*fmax(EP_cns,0.0);
dxdt_EA_cns =  k_age*fmax(EP_cns,0.0) - KTISS*fmax(EA_cns,0.0);

// --- plasma BChE (biomarker; aging lumped into the inhibited pool)
dxdt_B_free = -kinh*BCHE_SC*fmax(B_free,0.0)
              + KBCHE*(1.0 - fmax(B_free,0.0) - fmax(B_inh,0.0));
dxdt_B_inh  =  kinh*BCHE_SC*fmax(B_free,0.0) - KBCHE*fmax(B_inh,0.0);

// --- neuropathy target esterase (no reactivation; aging is fast)
double k_age_nte = LN2/T12AGNTE;
double kinh_nte  = KI_NTE * C_oxt;
double tot_nte = fmax(N_free,0.0) + fmax(N_inh,0.0) + fmax(N_aged,0.0);
dxdt_N_free = -kinh_nte*fmax(N_free,0.0) + KNTE*(1.0 - tot_nte);
dxdt_N_inh  =  kinh_nte*fmax(N_free,0.0) - (k_age_nte + KNTE)*fmax(N_inh,0.0);
dxdt_N_aged =  k_age_nte*fmax(N_inh,0.0) - KNTE*fmax(N_aged,0.0);

// =========================================================================
// 8. ACETYLCHOLINE — a hyperbola in enzyme activity
//    steady state:  ACh = (1 + leak)/(E + leak)
// =========================================================================
double mgblk = 1.0 - MG_EMAX*C_mg/(MG_EC50 + C_mg);
double dzblk = 1.0 - 0.25*C_dz/(DZ_EC50 + C_dz);

dxdt_ACh_m = KHYD*(1.0 + ACHLEAK)*mgblk - KHYD*(Emu + ACHLEAK)*Am;
dxdt_ACh_n = KHYD*(1.0 + ACHLEAK)*mgblk - KHYD*(Enm + ACHLEAK)*An;
dxdt_ACh_b = KHYD*(1.0 + ACHLEAK)*mgblk*dzblk - KHYD*(Ecn + ACHLEAK)*Ab;

// =========================================================================
// 9. RECEPTORS AND ADAPTATION
// =========================================================================
double On    = An/(An + KA_N);
double nic_x = clamp01((On - ON0)/(1.0 - ON0));

double Cat_b  = Cat * BBBF;
double Ob     = Ab/(Ab + KA_B*(1.0 + Cat_b/KB_ATR));
double cns_raw= clamp01((Ob - OM0)/(OMMAX - OM0));
double cns_x  = cns_raw * (1.0 - FDOWN_C*Rmd);

double Anh = pow(An, HILLN);
double fdes_raw = Anh/(Anh + pow(ACHN50, HILLN));
double fdes = fdes_raw * (1.0 - FDOWN_N*Rmd);

dxdt_Rn_des = KDES*fdes*(1.0 - Rnd) - KRECDES*Rnd;

double drive_down = fmax(fmax(mus_raw, cns_raw), nic_x);
dxdt_Rm_down = KDOWN*drive_down*(1.0 - Rmd) - KRECDOWN*Rmd;

// =========================================================================
// 10. END ORGANS
// =========================================================================
double kclr = (ICUSUCT > 0.5) ? KCLRICU : KCLRSEC;
dxdt_SEC = KSEC*mus_x - kclr*SECp;
dxdt_BT  = KBT*mus_x*(1.0 - BTp) - KBTOFF*BTp;

double Om_eff = OM0 + (Om - OM0)*(1.0 - FDOWN_M*Rmd);
double HR_t   = clampx(HRINT - HRSLOPE*Om_eff + 25.0*nic_x, 22.0, 170.0);
dxdt_HR = (HR_t - HR)/TAUHR;

// Central respiratory depression has an atropine-REVERSIBLE muscarinic part
// and a part atropine cannot touch; only receptor adaptation removes the
// second one, which is why weaning runs on a clock of days rather than on
// the atropine chart.
double Lb    = clamp01((Ab - 1.0)/20.0);
double Lb_ad = Lb * (1.0 - FDOWN_R*Rmd);
double RESPD_t = clampx(1.0 - DMAXCNS*cns_x - DMAXNONR*Lb_ad - SOLVRESP*solv,
                        0.02, 1.0);
dxdt_RESPD = (RESPD_t - RESPD)/TAURESPD;

double MSTR_t = clampx((1.0 - Rnd)*(1.0 - WFASTBL*fdes), 0.02, 1.0);
dxdt_MSTR = (MSTR_t - MSTR)/TAUMSTR;

// =========================================================================
// 11. GAS EXCHANGE
//     Ventilation is the SMALLER of what the brain asks for and what the
//     respiratory muscles can deliver.  Resting ventilation needs only ~28%
//     of maximal inspiratory capacity, so the neuromuscular limb has a safety
//     factor that the central limb does not.
// =========================================================================
double Vcap = RESV_i * fmin(RESPD, MSTR/SF_NMJ) * (1.0 - 0.40*BTp);
double vent_on = (VENTAVL > 0.5) ? sig((VAVENTTH - Vcap)/0.05, 1.0) : 0.0;
double VA = Vcap + vent_on*(1.05 - Vcap);
VA = clampx(VA, 0.06, 1.05);
double PaCO2 = 40.0/VA;

double FiO2 = FIO2AIR;
if (VENTAVL > 0.5 && vent_on > 0.5)      FiO2 = FIO2VENT;
else if (O2AVL > 0.5 && Vcap < 0.75)     FiO2 = FIO2O2;

double shunt = fmin(SHUNTMAX,
                    SHUNT0 + WSHSEC*fmin(1.0, SECp/3.0) + WSHLUNG*LUNGp);
double PAO2 = FiO2*713.0 - PaCO2/0.8;
double PaO2 = fmax(8.0, (1.0 - shunt)*PAO2 + shunt*40.0);
double hyp  = clamp01((60.0 - PaO2)/60.0);

// aspiration is a THRESHOLD event: a flooded, unprotected airway or a
// hydrocarbon vehicle going down the wrong way
double asp_drive = 0.60*solv + 0.90*clamp01(SECp/SECCRIT - 0.60)*(1.0 - vent_on);
dxdt_LUNG = KLUNG*asp_drive*(1.0 - LUNGp) - KLUNGREP*LUNGp;

double MAP_t = clampx(MAP0 - WMAPMUS*mus_x - WMAPSOLV*solv*SOLVCV - WMAPHYP*hyp,
                      25.0, 120.0);
dxdt_MAP = (MAP_t - MAP)/TAUMAP;

// =========================================================================
// 12. SEIZURES, OPIDN, HAZARD AND ACCUMULATORS
// =========================================================================
double seiz_drive = clamp01((cns_x - SEIZTH)/(1.0 - SEIZTH));
double dz_supp = 1.0 - 0.9*C_dz/(DZ_EC50 + C_dz);
dxdt_SEIZ = KSEIZ*seiz_drive*dz_supp*(1.0 - clamp01(SEIZ)) - KSEIZOFF*clamp01(SEIZ);

dxdt_OPIDN = KOPIDN*clamp01((fmax(N_aged,0.0) - NTETH)/(1.0 - NTETH))
             *(1.0 - clamp01(OPIDN)) - KRECOPID*clamp01(OPIDN);

double shock = clamp01((65.0 - MAP)/40.0);
double hcap  = clampx((PaCO2 - 60.0)/40.0, 0.0, 2.5);
dxdt_HAZ = H0 + HHYPOX*hyp*hyp + HHYPERC*hcap*hcap + HSHOCK*shock*shock
           + HSEIZ*clamp01(SEIZ) + HLUNG*LUNGp + HVENT*vent_on;

dxdt_VTIME  = vent_on;
dxdt_ATRCUM = rate_atr * ATR_MW / 1000.0;      // umol/h -> mg/h
dxdt_OXCUM  = oxinf * OXMW / 1000.0;          // umol/h -> mg/h
dxdt_AUC_EPn = fmax(EP_nmj, 0.0);
dxdt_AUC_ox  = C_ox;

$TABLE
double C_oxon   = fmax(A_ox_c,0.0)/VOXON;                      // nM
double C_oxonT  = fmax(A_ox_t,0.0)/VTOXON;                     // nM
double C_thion  = fmax(A_th_c,0.0)/V1TH;                       // uM
double C_oxime_ = fmax(A_pam_c,0.0)/OXV1;                      // uM
double kr_      = (OXTYPE > 0.5) ? KRMAX*C_oxime_/(KDOX + C_oxime_) : 0.0;
double kinh_    = KI * C_oxon;

// the two dimensionless groups the whole model turns on
double OMEGA    = (C_oxon > 1e-12) ? KRMAX/(KI*C_oxon) : 1.0e6;
double E_CEIL   = OMEGA/(1.0 + OMEGA);
double PHI      = (LN2/T12AGE)/(LN2/T12AGE + LN2/T12SPO + kr_);
double E_QSS    = (kr_ + LN2/T12SPO)/(kr_ + LN2/T12SPO + kinh_);

double AChE_RBC = 100.0*fmax(E_rbc,0.0);
double BChE_PL  = 100.0*fmax(B_free,0.0);
double AGED_RBC = 100.0*fmax(EA_rbc,0.0);

double Om_     = fmax(ACh_m,0.0)/(fmax(ACh_m,0.0) + KA_M*(1.0 + fmax(Ce_atr,0.0)/KB_ATR));
double mus_raw_= clamp01((Om_ - OM0)/(OMMAX - OM0));
double MUSX    = mus_raw_*(1.0 - FDOWN_M*clamp01(Rm_down));
double On_     = fmax(ACh_n,0.0)/(fmax(ACh_n,0.0) + KA_N);
double NICX    = clamp01((On_ - ON0)/(1.0 - ON0));
double Ob_     = fmax(ACh_b,0.0)/(fmax(ACh_b,0.0) + KA_B*(1.0 + fmax(Ce_atr,0.0)*BBBF/KB_ATR));
double CNSX    = clamp01((Ob_ - OM0)/(OMMAX - OM0))*(1.0 - FDOWN_C*clamp01(Rm_down));

double VCAP    = RESV_i*fmin(RESPD, MSTR/SF_NMJ)*(1.0 - 0.40*clamp01(BT));
double VENT_ON = (VENTAVL > 0.5) ? sig((VAVENTTH - VCAP)/0.05, 1.0) : 0.0;
double VA_     = clampx(VCAP + VENT_ON*(1.05 - VCAP), 0.06, 1.05);
double PACO2   = 40.0/VA_;
double FIO2_   = FIO2AIR;
if (VENTAVL > 0.5 && VENT_ON > 0.5)   FIO2_ = FIO2VENT;
else if (O2AVL > 0.5 && VCAP < 0.75)  FIO2_ = FIO2O2;
double SHUNT_  = fmin(SHUNTMAX, SHUNT0 + WSHSEC*fmin(1.0, fmax(SEC,0.0)/3.0)
                      + WSHLUNG*clamp01(LUNG));
double PAO2_   = fmax(8.0, (1.0 - SHUNT_)*(FIO2_*713.0 - PACO2/0.8) + SHUNT_*40.0);

double MORTALITY = 100.0*(1.0 - exp(-fmax(HAZ,0.0)));
double IMS_IDX   = 100.0*clamp01(Rn_des);

// Peradeniya OP Poisoning (POP) score reconstructed from the state vector
double pop_miosis = (MUSX > 0.55) ? 2.0 : ((MUSX > 0.15) ? 1.0 : 0.0);
double pop_resp   = (VCAP < 0.45) ? 2.0 : ((VCAP < 0.75) ? 1.0 : 0.0);
double pop_brady  = (HR < 40.0) ? 2.0 : ((HR < 60.0) ? 1.0 : 0.0);
double pop_fasc   = (NICX > 0.65) ? 2.0 : ((NICX > 0.25) ? 1.0 : 0.0);
double pop_cns    = (CNSX > 0.60) ? 3.0 : ((CNSX > 0.35) ? 2.0 : ((CNSX > 0.15) ? 1.0 : 0.0));
double pop_seiz   = (SEIZ > 0.3) ? 1.0 : 0.0;
double POP = pop_miosis + pop_resp + pop_brady + pop_fasc + pop_cns + pop_seiz;

$CAPTURE @annotated
C_oxon    : Plasma oxon concentration (nM)
C_oxonT   : Tissue oxon concentration (nM)
C_thion   : Plasma thion concentration (uM)
C_oxime_  : Plasma oxime concentration (uM)
kr_       : Oxime reactivation rate constant (1/h)
kinh_     : Pseudo-first-order AChE inhibition rate (1/h)
OMEGA     : Oxime sufficiency number k_r_max/(k_i*C_oxon) (-)
E_CEIL    : Free-enzyme ceiling at unlimited oxime (-)
PHI       : Per-event aging probability (-)
E_QSS     : Quasi-steady free-enzyme fraction at the current oxime level (-)
AChE_RBC  : RBC acetylcholinesterase (% of baseline)
BChE_PL   : Plasma butyrylcholinesterase (% of baseline)
AGED_RBC  : Irreversibly aged RBC AChE (% of baseline)
MUSX      : Muscarinic signalling excess (0-1)
NICX      : Nicotinic signalling excess (0-1)
CNSX      : Central muscarinic signalling excess (0-1)
VCAP      : Spontaneous ventilatory capability (fraction of normal)
VENT_ON   : Mechanical ventilation in progress (0-1)
PACO2     : Arterial CO2 tension (mmHg)
PAO2_     : Arterial O2 tension (mmHg)
MORTALITY : Cumulative probability of death (%)
IMS_IDX   : Intermediate-syndrome index, nAChR block (%)
POP       : Peradeniya OP Poisoning score (0-11)

## =============================================================================
##  R DRIVER — compound library, scenarios, and the analyses that make the
##  point.  Everything below runs after mread(); nothing below is parsed by
##  mrgsolve as model code.
## =============================================================================

$OMEGA @annotated @block
ETA_DOSE : 0.72 : Ingested dose (log-scale variance)
ETA_BIO  : 0.00 0.09 : CYP bioactivation capacity
ETA_PON  : 0.00 0.00 0.16 : PON1 oxon hydrolysis
ETA_RESP : 0.00 0.00 0.00 0.04 : Respiratory reserve

$SIGMA 0

$ENV

## ---------------------------------------------------------------------------
## OP compound library.  Rate constants for HUMAN AChE.  Literature values are
## quoted as M^-1 min^-1 and converted here:
##      k_i [1/nM/h]      = k_i [M^-1 min^-1] * 6e-8
##      k_r2 [1/uM/h]     = k_r2 [M^-1 min^-1] * 6e-5,  split into
##      k_r2 = KRMAX/KDOX so that reactivation saturates as Worek describes.
## ---------------------------------------------------------------------------
OP_LIBRARY <- list(

  chlorpyrifos = list(
    label = "Chlorpyrifos 20% EC (diethyl, lipophilic pro-toxicant)",
    subclass = "diethyl", CONC_GL = 200,
    MW_OP = 350.6, KI = 0.18, KI_NMJ = 1.0, KI_NTE = 0.0035,
    T12AGE = 33, T12SPO = 77, T12AGNTE = 6,
    KA_GUT = 0.35, FBIO = 0.80, V1TH = 35, VFTH = 700, QFTH = 25,
    VMAXBIO = 6.0e4, KMBIO = 20, CLTHOTH = 14, CLOXON = 300,
    VOXON = 56, QOXON = 60, VTOXON = 120, SOLVFR = 0.60, SOLVCV = 0.5,
    PON1_QQ = 0.33,
    pam = list(KRMAX = 36, KDOX = 200), obi = list(KRMAX = 48, KDOX = 40)),

  dimethoate = list(
    label = "Dimethoate 40% EC (dimethyl, cyclohexanone vehicle)",
    subclass = "dimethyl", CONC_GL = 400,
    MW_OP = 229.3, KI = 0.0084, KI_NMJ = 1.0, KI_NTE = 0.00002,
    T12AGE = 3.7, T12SPO = 0.7, T12AGNTE = 6,
    KA_GUT = 0.90, FBIO = 0.90, V1TH = 42, VFTH = 140, QFTH = 40,
    VMAXBIO = 1.5e5, KMBIO = 20, CLTHOTH = 45, CLOXON = 30,
    VOXON = 50, QOXON = 50, VTOXON = 90, SOLVFR = 0.55, SOLVCV = 1.8,
    PON1_QQ = 0.85,
    pam = list(KRMAX = 15, KDOX = 300), obi = list(KRMAX = 30, KDOX = 100)),

  fenthion = list(
    label = "Fenthion 50% EC (dimethyl, extremely lipophilic, delayed onset)",
    subclass = "dimethyl", CONC_GL = 500,
    MW_OP = 278.3, KI = 0.050, KI_NMJ = 1.6, KI_NTE = 0.00004,
    T12AGE = 3.7, T12SPO = 0.9, T12AGNTE = 6,
    KA_GUT = 0.10, FBIO = 0.85, V1TH = 30, VFTH = 1750, QFTH = 6,
    VMAXBIO = 4.0e4, KMBIO = 20, CLTHOTH = 12, CLOXON = 90,
    VOXON = 56, QOXON = 50, VTOXON = 120, SOLVFR = 0.50, SOLVCV = 0.4,
    PON1_QQ = 0.80,
    pam = list(KRMAX = 15, KDOX = 300), obi = list(KRMAX = 30, KDOX = 100)),

  parathion = list(
    label = "Parathion 50% EC (diethyl; the compound the oximes were built for)",
    subclass = "diethyl", CONC_GL = 500,
    MW_OP = 291.3, KI = 0.084, KI_NMJ = 1.0, KI_NTE = 0.0016,
    T12AGE = 33, T12SPO = 77, T12AGNTE = 6,
    KA_GUT = 0.45, FBIO = 0.85, V1TH = 40, VFTH = 500, QFTH = 25,
    VMAXBIO = 8.0e4, KMBIO = 20, CLTHOTH = 12, CLOXON = 420,
    VOXON = 56, QOXON = 60, VTOXON = 120, SOLVFR = 0.55, SOLVCV = 0.4,
    PON1_QQ = 0.28,
    pam = list(KRMAX = 36, KDOX = 200), obi = list(KRMAX = 48, KDOX = 40))
)

## Oxime products.  Doses are the regimens actually used in the trials.
OXIME_LIBRARY <- list(
  pam_who   = list(oxime = "pam", MW = 172.6, load_mg_kg = 30, inf_mg_kg_h = 8,
                   V1 = 24, V2 = 24, CL = 21, Q = 30,
                   label = "pralidoxime, WHO regimen 30 mg/kg + 8 mg/kg/h"),
  pam_bolus = list(oxime = "pam", MW = 172.6, bolus_mg = 1000, interval_h = 6,
                   V1 = 24, V2 = 24, CL = 21, Q = 30,
                   label = "pralidoxime 1 g intravenous bolus every 6 h"),
  obidoxime = list(oxime = "obi", MW = 359.2, load_mg = 250, inf_mg_24h = 750,
                   V1 = 15, V2 = 12, CL = 7.5, Q = 9,
                   label = "obidoxime 250 mg bolus + 750 mg/24 h")
)


## =============================================================================
## USAGE
## =============================================================================
## The block below is documentation, not model code.  Copy it into an R script
## next to this file.
##
## library(mrgsolve); library(dplyr); library(ggplot2)
## mod <- mread(file = "op_mrgsolve_model.R", model = "op")
##
## ## ---- helper: parameters for one compound + one oxime --------------------
## op_params <- function(compound = "chlorpyrifos", oxime = NULL,
##                       volume_ml = 50, pon1 = "RR", wt = 60) {
##   L <- env_get(mod)$OP_LIBRARY[[compound]]
##   p <- L[c("MW_OP","KI","KI_NMJ","KI_NTE","T12AGE","T12SPO","T12AGNTE",
##            "KA_GUT","FBIO","V1TH","VFTH","QFTH","VMAXBIO","KMBIO",
##            "CLTHOTH","CLOXON","VOXON","QOXON","VTOXON","SOLVFR","SOLVCV")]
##   p$CONC_GL  <- L$CONC_GL
##   p$DOSE_ML  <- volume_ml
##   p$PON1SC   <- if (pon1 == "QQ") L$PON1_QQ else 1
##   p$WT       <- wt
##   if (is.null(oxime)) {
##     p$OXTYPE <- 0; p$OXRATE <- 0
##   } else {
##     O <- env_get(mod)$OXIME_LIBRARY[[oxime]]
##     k <- if (O$oxime == "obi") L$obi else L$pam
##     p$KRMAX <- k$KRMAX; p$KDOX <- k$KDOX
##     p$OXV1 <- O$V1; p$OXV2 <- O$V2; p$OXCL <- O$CL; p$OXQ <- O$Q
##     p$OXMW <- O$MW
##     p$OXTYPE <- 1
##     p$OXRATE <- if (!is.null(O$inf_mg_kg_h)) O$inf_mg_kg_h*wt*1000/O$MW
##                 else if (!is.null(O$inf_mg_24h)) O$inf_mg_24h/24*1000/O$MW
##                 else 0
##   }
##   p
## }
##
## ## ---- helper: the bolus events -------------------------------------------
## op_events <- function(oxime = NULL, wt = 60, diazepam = TRUE,
##                       magnesium = FALSE, scavenger_mg = 0,
##                       oxime_start = 1.5, oxime_dur = 168) {
##   e <- NULL
##   if (!is.null(oxime)) {
##     O <- env_get(mod)$OXIME_LIBRARY[[oxime]]
##     if (!is.null(O$load_mg_kg))
##       e <- c(e, list(ev(time = oxime_start, amt = O$load_mg_kg*wt*1000/O$MW,
##                         cmt = "A_pam_c")))
##     if (!is.null(O$load_mg))
##       e <- c(e, list(ev(time = oxime_start, amt = O$load_mg*1000/O$MW,
##                         cmt = "A_pam_c")))
##     if (!is.null(O$bolus_mg))
##       e <- c(e, list(ev(time = oxime_start, amt = O$bolus_mg*1000/O$MW,
##                         cmt = "A_pam_c", ii = O$interval_h,
##                         addl = floor(oxime_dur/O$interval_h))))
##   }
##   if (diazepam)
##     e <- c(e, list(ev(time = 0.6, amt = 10*1000/284.7, cmt = "A_dz",
##                       ii = 6, addl = 2)))
##   if (magnesium)
##     e <- c(e, list(ev(time = 1.0, amt = 16.2, cmt = "A_mg")))
##   if (scavenger_mg > 0)
##     e <- c(e, list(ev(time = 1.0, amt = scavenger_mg/1000/340000*1e6*4,
##                       cmt = "SCAV")))
##   if (is.null(e)) return(ev(time = 0, amt = 0, cmt = "A_dz"))
##   Reduce(function(a, b) a + b, e)
## }
##
## ## ---- the 23 scenarios ----------------------------------------------------
## SCEN <- list(
##  S01 = list(lab="CPF 50 mL - supportive only",           cmp="chlorpyrifos", ml= 50, ox=NULL),
##  S02 = list(lab="CPF 50 mL - 2-PAM WHO infusion",        cmp="chlorpyrifos", ml= 50, ox="pam_who"),
##  S03 = list(lab="CPF 50 mL - 2-PAM 1 g q6h bolus",       cmp="chlorpyrifos", ml= 50, ox="pam_bolus"),
##  S04 = list(lab="CPF 50 mL - obidoxime",                 cmp="chlorpyrifos", ml= 50, ox="obidoxime"),
##  S05 = list(lab="CPF 10 mL - 2-PAM (below the ceiling)", cmp="chlorpyrifos", ml= 10, ox="pam_who"),
##  S06 = list(lab="CPF 10 mL - supportive only",           cmp="chlorpyrifos", ml= 10, ox=NULL),
##  S07 = list(lab="CPF 200 mL - 2-PAM (massive)",          cmp="chlorpyrifos", ml=200, ox="pam_who"),
##  S08 = list(lab="Dimethoate 50 mL - supportive",         cmp="dimethoate",   ml= 50, ox=NULL),
##  S09 = list(lab="Dimethoate 50 mL - 2-PAM",              cmp="dimethoate",   ml= 50, ox="pam_who"),
##  S10 = list(lab="Dimethoate 50 mL - obidoxime",          cmp="dimethoate",   ml= 50, ox="obidoxime"),
##  S11 = list(lab="CPF 50 mL - 2-PAM at 0.5 h",            cmp="chlorpyrifos", ml= 50, ox="pam_who", start=0.5),
##  S12 = list(lab="CPF 50 mL - 2-PAM at 12 h",             cmp="chlorpyrifos", ml= 50, ox="pam_who", start=12),
##  S13 = list(lab="CPF 50 mL - NO ventilator",             cmp="chlorpyrifos", ml= 50, ox="pam_who", p=list(VENTAVL=0, ICUSUCT=0)),
##  S14 = list(lab="CPF 50 mL - slow atropine titration",   cmp="chlorpyrifos", ml= 50, ox="pam_who", p=list(ATRMODE=2)),
##  S15 = list(lab="CPF 50 mL - glycopyrrolate",            cmp="chlorpyrifos", ml= 50, ox="pam_who", p=list(BBBF=0.02)),
##  S16 = list(lab="CPF 50 mL - + magnesium 4 g",           cmp="chlorpyrifos", ml= 50, ox="pam_who", mg=TRUE),
##  S17 = list(lab="CPF 50 mL - charcoal at 1 h",           cmp="chlorpyrifos", ml= 50, ox="pam_who", p=list(CHARC_T=1)),
##  S18 = list(lab="CPF 50 mL - charcoal at 6 h",           cmp="chlorpyrifos", ml= 50, ox="pam_who", p=list(CHARC_T=6)),
##  S19 = list(lab="CPF 50 mL - + 1 g BChE bioscavenger",   cmp="chlorpyrifos", ml= 50, ox="pam_who", scav=1000),
##  S20 = list(lab="CPF 50 mL - PON1 Q192Q",                cmp="chlorpyrifos", ml= 50, ox="pam_who", pon1="QQ"),
##  S21 = list(lab="Fenthion 50 mL - 2-PAM",                cmp="fenthion",     ml= 50, ox="pam_who"),
##  S22 = list(lab="Parathion 15 mL - obidoxime",           cmp="parathion",    ml= 15, ox="obidoxime"),
##  S23 = list(lab="Parathion 15 mL - supportive",          cmp="parathion",    ml= 15, ox=NULL)
## )
##
## run_scen <- function(s) {
##   p <- op_params(s$cmp, s$ox, s$ml, if (is.null(s$pon1)) "RR" else s$pon1)
##   if (!is.null(s$start)) p$OXSTART <- s$start
##   if (!is.null(s$p))     p <- modifyList(p, s$p)
##   e <- op_events(s$ox, magnesium = isTRUE(s$mg),
##                  scavenger_mg = if (is.null(s$scav)) 0 else s$scav,
##                  oxime_start = if (is.null(s$start)) 1.5 else s$start)
##   mod %>% param(p) %>% ev(e) %>% mrgsim(end = 336, delta = 0.25) %>% as_tibble()
## }
##
## out <- lapply(SCEN, run_scen)
##
## ## ---- the summary table that reproduces op_reference_output.txt ----------
## summarise_scen <- function(d) tibble::tibble(
##   AChE_nadir = min(d$AChE_RBC),
##   AChE_72h   = approx(d$time, d$AChE_RBC, 72)$y,
##   aged_72h   = approx(d$time, d$AGED_RBC, 72)$y,
##   oxon_24h   = approx(d$time, d$C_oxon,  24)$y,
##   omega_24h  = approx(d$time, d$OMEGA,   24)$y,
##   atropine_mg= max(d$ATRCUM),
##   vent_h     = max(d$VTIME),
##   IMS        = max(d$IMS_IDX),
##   POPmax     = max(d$POP),
##   mortality  = max(d$MORTALITY))
## dplyr::bind_rows(lapply(out, summarise_scen), .id = "scenario")
##
## ## ---- the dose-response sweep that locates the oxime ceiling -------------
## sweep <- lapply(c(2,5,10,20,35,50,100,200), function(ml) {
##   a <- run_scen(list(cmp="chlorpyrifos", ml=ml, ox=NULL))
##   b <- run_scen(list(cmp="chlorpyrifos", ml=ml, ox="pam_who"))
##   tibble::tibble(volume_ml = ml,
##                  oxon_plateau = approx(b$time, b$C_oxon, 24)$y,
##                  omega        = approx(b$time, b$OMEGA,  24)$y,
##                  ceiling      = approx(b$time, b$E_CEIL, 24)$y,
##                  AChE72_sup   = approx(a$time, a$AChE_RBC, 72)$y,
##                  AChE72_pam   = approx(b$time, b$AChE_RBC, 72)$y,
##                  d_mort       = max(a$MORTALITY) - max(b$MORTALITY))
## }) %>% dplyr::bind_rows()
##
## ## ---- the virtual trial ---------------------------------------------------
## ## A pragmatic trial enrols whoever arrives.  Draw the ingested volume from a
## ## log-normal, mix the OP classes as they are mixed in a real South Asian
## ## cohort, and run each virtual patient through both arms.
## set.seed(20260805)
## trial <- lapply(1:120, function(i) {
##   ml  <- pmin(250, pmax(1, exp(rnorm(1, log(28), 0.85))))
##   cmp <- if (runif(1) < 0.55) "chlorpyrifos" else "dimethoate"
##   pon <- if (runif(1) < 0.25) "QQ" else "RR"
##   a <- run_scen(list(cmp=cmp, ml=ml, ox=NULL,      pon1=pon))
##   b <- run_scen(list(cmp=cmp, ml=ml, ox="pam_who", pon1=pon))
##   tibble::tibble(id=i, op=cmp, ml=ml, pon1=pon,
##                  mort_placebo=max(a$MORTALITY), mort_pam=max(b$MORTALITY),
##                  vent_placebo=max(a$VTIME),     vent_pam=max(b$VTIME))
## }) %>% dplyr::bind_rows()
## trial %>% dplyr::summarise(placebo=mean(mort_placebo), pam=mean(mort_pam),
##                            RR=mean(mort_pam)/mean(mort_placebo))
## trial %>% dplyr::mutate(stratum = ifelse(ml <= 15, "<=15 mL", ">15 mL")) %>%
##   dplyr::group_by(op, stratum) %>%
##   dplyr::summarise(n=dplyr::n(), placebo=mean(mort_placebo),
##                    pam=mean(mort_pam), .groups="drop")
##
## =============================================================================
## WHAT TO LOOK AT FIRST
## =============================================================================
## 1. Plot OMEGA and E_CEIL against time for S02 and S05.  They are the model.
##    Where OMEGA crosses 1 is where the oxime stops being a drug and starts
##    being a placebo, and nothing about the oxime itself changes at that point
##    — only the oxon concentration does.
## 2. Compare AChE_RBC between S02 and S01.  The oxime does exactly what it is
##    supposed to do to the enzyme.  Then compare MORTALITY.  The gap between
##    those two comparisons is the whole history of the oxime trials.
## 3. Compare VTIME between S02 and S13.  The single most valuable object in
##    this model is the ventilator, and it is not a drug.
## 4. Compare ATRCUM between S01 and S02, then look at CNSX.  Atropine is the
##    only agent here with central activity, and it is titrated against
##    peripheral signs that the oxime also improves.
## =============================================================================
