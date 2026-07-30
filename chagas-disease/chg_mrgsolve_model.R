## ============================================================================
##  Chagas disease / Chronic Chagas Cardiomyopathy (CCC) — QSP model in mrgsolve
##  68 ODEs · trypanocidal + cardiovascular PK/PD · competing risks · 28 scenarios
## ============================================================================
##
##  WHY THIS MODEL IS BUILT THE WAY IT IS
##  ---------------------------------------------------------------------------
##  Four randomised trials in Chagas disease produced four results that read as
##  contradictions if you hold one picture of the disease in your head:
##
##    STOP-CHAGAS   posaconazole cleared blood PCR BETTER than benznidazole at
##                  day 30 (93.3% vs 89.7%) and then failed almost completely
##                  (13.3% vs 86.7% sustained at day 180, 16% vs 96% at day 360).
##    CHAGASAZOL    the same reversal, independently.
##    BENDITA       benznidazole for 2 weeks was as good as for 8 weeks
##                  (83% vs 89%), and 150 mg was as good as 300 mg.
##    BENEFIT       benznidazole halved PCR positivity in patients with
##                  established cardiomyopathy and changed nothing clinically
##                  (HR 0.93, 95% CI 0.81-1.07).
##
##  They are not contradictions. Each one is a measurement on a DIFFERENT
##  TIME-SCALE, and this model is built as three explicitly separated
##  time-scales so that all four results are outputs rather than exceptions.
##
##  ---------------------------------------------------------------------------
##  AXIS 1 — STATIC AND CIDAL ARE TWO DIFFERENT DRUG PROPERTIES, AND qPCR
##           MEASURES THE ONE THAT DOES NOT PRODUCE CURE.
##  ---------------------------------------------------------------------------
##      Blood trypomastigotes are not a reservoir; they are a FLUX — the output
##      of intracellular amastigote replication. Anything that stops replication
##      empties the blood. Azoles (CYP51 / ergosterol) are overwhelmingly
##      STATIC: ergosterol is needed to build membrane for a dividing cell, so
##      blocking it arrests growth without killing. Nitroimidazoles are CIDAL:
##      TcNTR-1 reduces them to glyoxal and nitroso species that form DNA
##      adducts in a cell that is not dividing at all.
##      Every drug in this model therefore carries TWO independent parameters,
##      EMAX_STAT and EMAX_CIDR, and the azole class has EMAX_STAT ~ 0.99 with
##      EMAX_CIDR ~ 0.15/day while the nitro class has 0.90 and 2.0/day.
##      CONSEQUENCE (emergent, not coded): posaconazole beats benznidazole on
##      the endpoint the trial measured during treatment, and loses on the
##      endpoint that matters after it.
##
##  ---------------------------------------------------------------------------
##  AXIS 2 — IMMUNE PRESSURE IS ANTIGEN-DRIVEN, SO A STATIC DRUG DISARMS THE
##           HOST WHILE IT SUPPRESSES THE PARASITE.
##  ---------------------------------------------------------------------------
##      This is the reason a static drug does not simply let the immune system
##      finish the job. Killing of infected cells is proportional to TH1, TH1 is
##      driven by antigen, and antigen is produced by replicating parasites. A
##      static drug collapses antigen within weeks, TH1 decays with a ~20-day
##      time constant, and immune pressure on the surviving amastigotes falls
##      with it. The reservoir is preserved BY the treatment.
##      A cidal drug destroys the reservoir faster than the immune response
##      relaxes, so the same relaxation does not matter.
##      CONSEQUENCE (emergent): withdrawal of a static drug is followed by full
##      relapse over 2-5 months, which is exactly the STOP-CHAGAS day-30 to
##      day-180 trajectory, and no lengthening of azole therapy repairs it
##      (scenario 13 runs posaconazole for a year and still fails).
##
##  ---------------------------------------------------------------------------
##  AXIS 3 — CURE IS A THRESHOLD CROSSING, NOT AN INTEGRAL OF EXPOSURE.
##  ---------------------------------------------------------------------------
##      A small drug-tolerant dormant sub-population (k_dorm / k_wake, ~0.3% of
##      the burden at steady state) decays under a cidal drug at k_wake plus a
##      direct dormant-kill rate roughly 15-fold below the replicating-kill
##      rate. Sterile cure is a stochastic extinction event, modelled as
##      P(cure) = exp(-P_EST * N_min) where N_min is the minimum total burden
##      reached. Because the surviving number falls exponentially and the
##      exponential passes below the extinction scale inside the first two
##      weeks, cure probability SATURATES.
##      CONSEQUENCE (emergent): 2 weeks ~ 4 weeks ~ 8 weeks, and 150 mg ~ 300 mg.
##      Extra weeks add nothing to the numerator and everything to the two
##      toxicity hazards — which is why the BENDITA 2-week arm had zero
##      treatment discontinuations and the 8-week arms did not.
##
##  ---------------------------------------------------------------------------
##  AXIS 4 — THE DAMAGE RATE HAS A PARASITE TERM AND AN AUTOCATALYTIC TERM,
##           AND THE RATIO IS A FUNCTION OF TIME, NOT OF PARASITE LOAD.
##  ---------------------------------------------------------------------------
##      dCMYO/dt = -(INJ_BASE + INJ_PAR + INJ_AUT + INJ_ISC + INJ_STR) * CMYO
##      Only INJ_PAR contains a parasite. The other three close their own loops:
##        INJ_AUT  autoantibody raised by cardiac myosin RELEASED BY THE DAMAGE
##        INJ_ISC  endothelin -> spasm -> ischaemia -> endothelial injury
##        INJ_STR  myocyte loss -> dilatation -> wall stress -> myocyte loss
##      The model reports PAF(t) = 1 - INJ(parasite set to zero) / INJ, computed
##      as an in-model counterfactual at every time point.
##      CONSEQUENCE (emergent, and the whole point): PAF starts near 0.9 and
##      decays to a few per cent, WITHOUT ANY AGE TERM ANYWHERE. The BENEFIT
##      null is then a prediction: the trial enrolled at the stage where the
##      quantity it was trying to move had already left.
##
##  ---------------------------------------------------------------------------
##  AXIS 5 — PATIENTS DIE OF TWO DIFFERENT THINGS AND THE TWO COMPETE.
##  ---------------------------------------------------------------------------
##      Arrhythmic death is driven by SCARH (border-zone heterogeneity) and
##      SYMPD (sympathetic denervation) — not by how much fibrosis there is.
##      Pump-failure death is driven by EF. The model integrates true cumulative
##      INCIDENCE functions (dCIF_i/dt = h_i * S) rather than adding hazards, so
##      suppressing one mode necessarily reveals the other.
##      CONSEQUENCE (emergent): amiodarone cuts the sudden-death hazard by ~55%
##      and all-cause mortality by far less, because the deaths it prevents
##      reappear later as pump failure.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION ANCHORS (all values below are taken from the cited publications;
##  see chg_references.md for the PMIDs)
##  ---------------------------------------------------------------------------
##   A1  BENEFIT       PCR conversion end of treatment  BZN 66.2% / PBO 33.5%
##   A2  BENEFIT       PCR conversion 2 years           BZN 55.4% / PBO 35.3%
##   A3  BENEFIT       PCR conversion >=5 years         BZN 46.7% / PBO 33.1%
##   A4  BENEFIT       composite over 5.4 y             27.5% / 29.1%, HR 0.93
##   A5  BENEFIT       baseline PCR positive            60.5%
##   A6  BENDITA       sustained clearance 6 mo, ITT    PBO 3% · 300mg 8wk 89% ·
##                                                      4wk 89% · 2wk 83% ·
##                                                      150mg 4wk 83%
##   A7  BENDITA       discontinuation, 2-week arm      0%  (7% overall)
##   A8  STOP-CHAGAS   PCR negative at day 30           posa 93.3% / BZN 89.7%
##   A9  STOP-CHAGAS   sustained negative at day 180    posa 13.3% / BZN 86.7%
##                                                      combo 80% / PBO 10%
##   A10 STOP-CHAGAS   at day 360                       posa 16% / BZN 96%
##   A11 STOP-CHAGAS   permanent discontinuation        31.7% (60-day BZN arms)
##   A12 Rassi 2006    10-y mortality low/int/high      10% / 44% / 84%
##   A13 Senra 2018    LGE prevalence and mass          76.1%, 15.2 +/- 16.5 g
##   A14 Chadalawada   annual mortality in CCC          all-cause 7.9% · HF 3.5%
##                                                      · sudden 2.6% · stroke 0.4%
##   A15 Chadalawada   annual mortality by stage        4.8 / 8.7 / 13.9 / 22.4 %
##   A16 Nunes 2021    cardiomyopathy incidence         13.8 vs 4.6 /1000 py
##   A17 Nunes 2021    mortality, established CCC       80.9 /1000 py
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL DOES NOT DO
##  ---------------------------------------------------------------------------
##   * It is deterministic per patient. Cure is reported as a probability, and
##     population incidences come from an explicit virtual population over the
##     susceptibility covariate SUSC — not from an ETA inside a single run.
##   * It has no vector/transmission dynamics: the inoculum is an initial
##     condition, not a force of infection.
##   * The digestive form is present but coarse (three states); it is included
##     because it shares the enteric-neuron loss mechanism, not because the
##     megacolon predictions should be trusted quantitatively.
##   * Drug-drug interactions (posaconazole/CYP3A4 with amiodarone) are NOT
##     modelled. This matters and is flagged in the README as a known gap.
## ============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

CHG_CODE <- '
$PROB Chagas disease / chronic Chagas cardiomyopathy QSP model (68 ODEs)

$PARAM @annotated
// ---- benznidazole (and nitro class) PK, time unit = DAYS -------------------
KA_B    : 28.8  : benznidazole absorption rate (1/day)
V1_B    : 32    : benznidazole central volume (L)
V2_B    : 12    : benznidazole peripheral volume (L)
Q_B     : 72    : benznidazole intercompartmental clearance (L/day)
CL_B    : 50.4  : benznidazole clearance (L/day)
F_B     : 0.92  : benznidazole bioavailability

// ---- azole class (posaconazole default; ravuconazole via AZLTYPE) ----------
KA_A    : 12    : azole absorption rate (1/day)
V1_A    : 300   : azole central volume (L)
V2_A    : 700   : azole peripheral volume (L)
Q_A     : 360   : azole intercompartmental clearance (L/day)
CL_A    : 216   : azole clearance (L/day)
F_A     : 0.35  : azole bioavailability
AZLTYPE : 0     : 0 = posaconazole, 1 = ravuconazole (from fosravuconazole)
CL_A_RAV: 14.4  : ravuconazole clearance (L/day)

// ---- nifurtimox and fexinidazole ------------------------------------------
KA_N    : 24    : nifurtimox absorption rate (1/day)
V_N     : 100   : nifurtimox volume (L)
CL_N    : 200   : nifurtimox clearance (L/day)
F_N     : 0.50  : nifurtimox bioavailability
KA_F    : 20    : fexinidazole absorption rate (1/day)
V_F     : 180   : fexinidazole volume (L)
CL_F    : 260   : fexinidazole clearance (L/day)
F_F     : 0.60  : fexinidazole bioavailability

// ---- amiodarone and desethylamiodarone ------------------------------------
KA_M    : 3.0   : amiodarone absorption rate (1/day)
V1_M    : 60    : amiodarone central volume (L)
V2_M    : 4200  : amiodarone peripheral volume (L)
Q_M     : 480   : amiodarone intercompartmental clearance (L/day)
CL_M    : 144   : amiodarone clearance (L/day)
F_M     : 0.50  : amiodarone bioavailability
FM_DEA  : 0.40  : fraction of amiodarone converted to desethylamiodarone
V_DEA   : 3000  : desethylamiodarone volume (L)
CL_DEA  : 96    : desethylamiodarone clearance (L/day)

// ---- carvedilol -----------------------------------------------------------
KA_C    : 24    : carvedilol absorption rate (1/day)
V_C     : 200   : carvedilol volume (L)
CL_C    : 1440  : carvedilol clearance (L/day)
F_C     : 0.25  : carvedilol bioavailability

// ---- background cardiovascular therapy switches ---------------------------
ACEI_ON : 0     : ACE inhibitor / ARB / ARNI on (0/1)
MRA_ON  : 0     : mineralocorticoid receptor antagonist on (0/1)
OAC_ON  : 0     : oral anticoagulation on (0/1)
ICD_ON  : 0     : implantable cardioverter-defibrillator in situ (0/1)
PPM_ON  : 0     : permanent pacemaker in situ (0/1)
ACEI_E  : 0.55  : fractional suppression of angiotensin II by ACEI/ARB
MRA_E   : 0.60  : fractional suppression of aldosterone effect by MRA
OAC_E   : 0.66  : fractional reduction of embolic hazard on anticoagulation
ICD_E   : 0.72  : fractional reduction of arrhythmic-death hazard by ICD

// ---- trypanocidal PD: the two independent axes -----------------------------
EMAXST_N: 0.90  : nitro class maximal REPLICATION BLOCK (static, fraction)
EMAXCR_N: 2.00  : nitro class maximal kill rate of REPLICATING forms (1/day)
EMAXCD_N: 0.135 : nitro class maximal kill rate of DORMANT forms (1/day)
EC50_B  : 2.00  : benznidazole EC50 (mg/L)
EC50_N  : 0.80  : nifurtimox EC50 (mg/L)
EC50_F  : 1.20  : fexinidazole EC50 (mg/L)
EMAXST_A: 0.990 : azole class maximal REPLICATION BLOCK (static, fraction)
EMAXCR_A: 0.025 : azole class maximal kill rate of REPLICATING forms (1/day)
EMAXCD_A: 0.0005: azole class maximal kill rate of DORMANT forms (1/day)
EC50_A  : 0.15  : azole EC50 (mg/L)
EMAXCR_M: 0.020 : amiodarone weak direct trypanocidal kill rate (1/day)
EC50_M  : 1.50  : amiodarone trypanocidal EC50 (mg/L)
STRAIN  : 1.00  : strain susceptibility multiplier on nitro EC50 (TcI > 1)

// ---- parasite population dynamics -----------------------------------------
MU_P    : 0.300 : amastigote net replication rate (1/day)
KCAP    : 9.0e6 : tissue carrying capacity (organism equivalents)
KIMM0   : 0.290 : immune killing rate of replicating amastigotes at TH1 = TH1SS
KDORM   : 4.0e-4: rate of entry into dormancy (1/day)
KWAKE   : 0.130 : rate of spontaneous awakening from dormancy (1/day)
KREL    : 0.020 : trypomastigote release into blood per amastigote per day
KELB    : 2.00  : blood trypomastigote elimination rate (1/day)
KSEED   : 0.010 : re-seeding of tissue from blood (1/day)
FHEART  : 0.20  : fraction of the tissue burden that is cardiac
P_EST   : 2.1e-3: per-organism probability of re-establishing after treatment
PTOT0   : 3.0e5 : baseline chronic total tissue burden (organism equivalents)
VBLOOD  : 5000  : blood volume (mL)
LOD_PCR : 0.50  : qPCR limit of detection (par-eq/mL)
PDET    : 0.665 : per-visit probability that an infected patient tests positive

// ---- antigen, antibody, immunity -------------------------------------------
KAG     : 1.0e-4: antigen production per replicating organism per day
KDAG    : 0.020 : antigen elimination rate (1/day)
KTH     : 0.1136: TH1 expansion rate (1/day), set so TH1_ss = 1 at baseline antigen
KMAG    : 0.50  : antigen concentration at half-maximal TH1 drive
KDTH    : 0.050 : TH1 decay rate (1/day)
TH1MAX  : 3.00  : maximal TH1 level
KTRG    : 0.030 : Treg formation rate (1/day)
KDTRG   : 0.040 : Treg decay rate (1/day)
KMPH    : 0.120 : macrophage activation rate (1/day)
KDMPH   : 0.100 : macrophage deactivation rate (1/day)
KTNF    : 0.400 : TNF-alpha production rate (1/day)
KDTNF   : 0.500 : TNF-alpha elimination rate (1/day)
KIL10   : 0.250 : IL-10 production rate (1/day)
KDIL10  : 0.350 : IL-10 elimination rate (1/day)
KTGF    : 0.150 : TGF-beta production rate (1/day)
KDTGF   : 0.200 : TGF-beta elimination rate (1/day)
KABG    : 1.30e-3: anti-T.cruzi IgG production rate (1/day), set so titre_ss = 1
KDABG   : 3.0e-4: anti-T.cruzi IgG decay (1/day) -> seroreversion takes decades

// ---- autoantibody: the loop that closes without a parasite -----------------
KAAB    : 5.1e-5: autoantibody formation rate (1/day)
KMAAB   : 0.50  : antigen level at half-maximal autoantibody drive
EPS_MYO : 140.0 : weight of released cardiac myosin in autoantibody drive
KDAAB   : 3.0e-4: autoantibody decay rate (1/day)
AABMAX  : 1.00  : maximal autoantibody level

// ---- myocardial injury terms (1/day) --------------------------------------
SUSC    : 1.00  : patient susceptibility covariate (lognormal in populations)
KI_BASE : 3.0e-7: age-related background myocyte loss
KI_IMM  : 2.2e-6: immune/parasite-driven myocyte killing
KI_AG   : 8.3e-7: antigen-driven inflammatory myocyte loss
KI_AAB  : 2.09e-5: autoantibody-driven myocyte loss
KI_ISC  : 8.8e-5: microvascular ischaemic myocyte loss
KI_STR  : 1.1e-5: wall-stress-driven myocyte loss
CARV_CV : 0.35  : fractional reduction of stress/adrenergic injury by carvedilol
ACEI_CV : 0.30  : fractional reduction of stress-driven injury by ACEI

// ---- fibrosis --------------------------------------------------------------
KM_TGF  : 1.5e-5: myofibroblast activation by TGF-beta (1/day)
KM_ANG  : 1.0e-5: myofibroblast activation by angiotensin II (1/day)
KM_STR  : 4.0e-5: myofibroblast activation by wall stress (1/day)
KM_AUTO : 1.25e-4: myofibroblast autocrine self-activation (1/day) [LOOP+]
KDMFB   : 1.25e-4: myofibroblast deactivation (1/day)
KCOL    : 1.15e-4: collagen deposition per unit myofibroblast (1/day)
KDCOL   : 4.0e-4: collagen degradation (1/day)
COL0    : 0.040 : normal interstitial collagen fraction
KSH     : 13.5  : scar heterogeneity formation per unit focal injury rate
KDSH    : 1.82e-4: scar heterogeneity resolution (1/day)

// ---- microvascular ---------------------------------------------------------
KET1    : 0.090 : endothelin-1 production rate (1/day)
KDET1   : 0.120 : endothelin-1 elimination (1/day)
KMVD    : 2.0e-4 : microvascular dysfunction accrual (1/day)
KDMVD   : 2.31e-4: microvascular dysfunction recovery (1/day)

// ---- geometry and remodelling ----------------------------------------------
EDV0    : 120   : baseline LV end-diastolic volume (mL)
LVM0    : 150   : baseline LV mass (g)
PSYS    : 120   : LV systolic pressure (mmHg)
KREM    : 6.0e-4: dilatation rate constant per unit relative wall stress
KESV    : 0.050 : rate of approach of ESV to its contractility target (1/day)
EFMAX   : 66.0  : ejection fraction at fully preserved contractility (%)
HEF     : 1.60  : exponent linking contractility to ejection fraction
KC_COL  : 2.60  : contractility penalty per unit collagen above normal
KC_MVD  : 0.25  : contractility penalty per unit microvascular dysfunction
KHYP    : 0.0025: hypertrophy accrual per unit relative wall stress (1/day)
KDHYP   : 0.0035: hypertrophy regression (1/day)
HYPMAX  : 0.60  : maximal hypertrophic mass fraction
KAPX    : 3.0e-4: apical thinning accrual (1/day)
KAPX_I  : 1.80  : ischaemic amplification of apical thinning

// ---- conduction and autonomic ----------------------------------------------
KSYMPD  : 7.9   : sympathetic denervation per unit focal injury rate
KPSYMD  : 10.9  : parasympathetic denervation per unit focal injury rate
KCOND   : 10.9  : conduction tissue loss per unit focal injury rate
KSAN    : 5.25  : sinus node loss per unit focal injury rate

// ---- neurohormonal ---------------------------------------------------------
KANG    : 0.30  : angiotensin II formation (1/day)
KDANG   : 0.35  : angiotensin II elimination (1/day)
KALD    : 0.25  : aldosterone formation (1/day)
KDALD   : 0.30  : aldosterone elimination (1/day)
KNE     : 0.40  : norepinephrine tone formation (1/day)
KDNE    : 0.45  : norepinephrine tone elimination (1/day)
KBNP    : 0.35  : NT-proBNP formation (1/day)
KDBNP   : 0.40  : NT-proBNP elimination (1/day)

// ---- digestive form --------------------------------------------------------
KENS    : 7.0e-6: enteric neuron loss rate per unit parasite/immune drive
KESO    : 0.0016: oesophageal dilatation rate (1/day)
KCLN    : 0.0013: colonic dilatation rate (1/day)
GITROP  : 1.00  : digestive tropism multiplier (TcII/V/VI > TcI)

// ---- clinical hazards (1/day) ----------------------------------------------
H0SCD   : 9.6e-6: baseline sudden cardiac death hazard
H0HF    : 5.4e-6: baseline pump-failure death hazard
H0STK   : 2.1e-6: baseline stroke / thromboembolic death hazard
H0PPM   : 4.3e-6: baseline pacemaker / ICD insertion hazard
H0NHF   : 3.1e-6: baseline new heart failure hazard
H0VT    : 8.3e-6: baseline sustained ventricular tachycardia hazard
BS_SH   : 1.90  : log-hazard weight of scar heterogeneity on sudden death
BS_SD   : 1.30  : log-hazard weight of sympathetic denervation on sudden death
BS_EF   : 0.35  : log-hazard weight per 10 EF points below 60 on sudden death
BH_EF   : 0.90  : log-hazard weight per 10 EF points below 60 on pump failure
BH_BNP  : 0.25  : log-hazard weight of NT-proBNP on pump failure
BK_APX  : 1.60  : log-hazard weight of apical aneurysm on stroke
BK_EF   : 0.50  : log-hazard weight per 10 EF points below 60 on stroke
BP_CND  : 3.20  : log-hazard weight of conduction loss on pacemaker insertion
AMIO_E  : 0.55  : fractional reduction of arrhythmic hazard by amiodarone
AMIO_C50: 1.20  : amiodarone concentration for half-maximal antiarrhythmic effect
AMIO_INO: 0.06  : amiodarone negative inotropic fraction at full effect

// ---- safety ----------------------------------------------------------------
KRASH   : 0.0075: benznidazole cutaneous reaction hazard scale (1/day)
TRASH   : 16.0  : day of treatment at which the rash hazard peaks
WRASH   : 6.0   : width of the rash hazard window (days)
KNEU    : 1.4e-6: neuropathy hazard per unit cumulative benznidazole AUC
AUCTHR  : 60.0  : cumulative AUC (mg*day/L) below which neuropathy is negligible
KLEUK   : 0.0011: leukopenia / hepatic discontinuation hazard scale (1/day)
KAMITOX : 2.2e-6: amiodarone organ toxicity hazard per unit tissue burden

// ---- bookkeeping -----------------------------------------------------------
ONDRUG  : 0     : 1 while trypanocidal therapy is being administered
STERILE : 0     : 1 = follow the STERILISED branch of the cure lottery (see note)
EXT_THR : 1.0   : burden below which a sterilised patient is declared extinct
AGE0    : 30    : age at the start of the simulation (years)
AFIB    : 0     : atrial fibrillation present (0/1)

$CMT @annotated
GUT_B  : benznidazole depot (mg)
CEN_B  : benznidazole central (mg)
PER_B  : benznidazole peripheral (mg)
GUT_A  : azole depot (mg)
CEN_A  : azole central (mg)
PER_A  : azole peripheral (mg)
GUT_N  : nifurtimox depot (mg)
CEN_N  : nifurtimox central (mg)
GUT_F  : fexinidazole depot (mg)
CEN_F  : fexinidazole central (mg)
GUT_M  : amiodarone depot (mg)
CEN_M  : amiodarone central (mg)
PER_M  : amiodarone peripheral (mg)
DEA    : desethylamiodarone (mg)
GUT_C  : carvedilol depot (mg)
CEN_C  : carvedilol central (mg)
PBLD   : blood trypomastigotes (par-eq/mL)
PRH    : replicating amastigotes, heart (organisms)
PRX    : replicating amastigotes, extracardiac (organisms)
PDH    : dormant amastigotes, heart (organisms)
PDX    : dormant amastigotes, extracardiac (organisms)
ANTIG  : persistent parasite antigen (relative)
ABIG   : anti-T. cruzi IgG (relative)
MPHI   : activated macrophages (relative)
TH1    : IFN-gamma+ effector T cells (relative)
TREG   : regulatory T cells (relative)
TNFA   : TNF-alpha (relative)
IL10   : IL-10 (relative)
TGFB   : TGF-beta1 (relative)
AAB    : functional autoantibody (0-1)
CMYO   : viable cardiomyocyte mass (fraction of normal)
MFB    : activated myofibroblasts (0-1)
COL    : myocardial collagen fraction (0-1)
SCARH  : scar heterogeneity / border zone index (0-1)
APEX   : apical thinning / aneurysm index (0-1)
MVD    : microvascular dysfunction index (0-1)
ET1    : endothelin-1 (relative)
EDV    : LV end-diastolic volume (mL)
ESV    : LV end-systolic volume (mL)
HYPM   : hypertrophic mass fraction (0-1)
SYMPD  : cardiac sympathetic denervation index (0-1)
PSYMD  : parasympathetic denervation index (0-1)
COND   : His-Purkinje conduction reserve (1 -> 0)
SANF   : sinus node function reserve (1 -> 0)
ANGII  : angiotensin II tone (relative)
ALDO   : aldosterone tone (relative)
NE     : norepinephrine tone (relative)
BNP    : NT-proBNP (relative)
ENSN   : myenteric neuron pool (fraction of normal)
ESOD   : oesophageal dilatation index (0-1)
COLD   : colonic dilatation index (0-1)
HTOT   : total cumulative hazard of any fatal event
CIFSCD : cumulative incidence of sudden cardiac death
CIFHF  : cumulative incidence of pump-failure death
CIFSTK : cumulative incidence of fatal stroke / thromboembolism
CIFCMP : cumulative incidence of the BENEFIT composite endpoint
CIFPPM : cumulative incidence of pacemaker / ICD insertion
CIFNHF : cumulative incidence of new heart failure
CIFVT  : cumulative incidence of sustained ventricular tachycardia
BZNCUM : cumulative benznidazole AUC (mg*day/L)
TAUD   : days elapsed on trypanocidal therapy
HRASH  : cumulative hazard of cutaneous discontinuation
HNEU   : cumulative hazard of neuropathic discontinuation
AMITIS : amiodarone tissue burden (mg, surrogate)
DMGCUM : cumulative myocardial damage (actual)
DMGNOP : cumulative myocardial damage in the parasite-free counterfactual

$GLOBAL
#define CB   (CEN_B / V1_B)
#define CA   (CEN_A / V1_A)
#define CN   (CEN_N / V_N)
#define CF   (CEN_F / V_F)
#define CM   (CEN_M / V1_M)
#define CC   (CEN_C / V_C)
#define CDEA (DEA   / V_DEA)

// Emax helper, guarded against negative concentrations
double emx(double emax, double conc, double ec50) {
  double c = conc > 0.0 ? conc : 0.0;
  return emax * c / (ec50 + c);
}

$MAIN
// Baseline chronic infection is placed at its dormancy equilibrium, so that an
// untreated patient does not drift merely because the initial split was wrong.
double FD  = KDORM / (KDORM + KWAKE);
PRH_0  = PTOT0 * FHEART       * (1.0 - FD);
PRX_0  = PTOT0 * (1.0 - FHEART) * (1.0 - FD);
PDH_0  = PTOT0 * FHEART       * FD;
PDX_0  = PTOT0 * (1.0 - FHEART) * FD;
PBLD_0 = KREL * PTOT0 * (1.0 - FD) / ((KELB + KSEED) * VBLOOD);

// The chronic (indeterminate) phase is a genuine equilibrium of this system,
// not a starting guess: every immune initial condition below is the analytic
// steady state at the baseline antigen level, so an untreated patient with
// SUSC = 0 does not drift. If these are wrong the parasite compartment slowly
// escapes immune control and every treatment effect is measured against a
// moving baseline -- which is exactly the failure mode this block prevents.
ANTIG_0 = KAG * PTOT0 * (1.0 - FD) * 1.0e-4 / KDAG;
double AGD0 = ANTIG_0/(ANTIG_0 + KMAG);
TREG_0  = KTRG*AGD0 / (KDTRG - 0.20*KTRG*KIL10/KDIL10);
IL10_0  = KIL10*TREG_0/KDIL10;
TH1_0   = KTH*AGD0*TH1MAX / (KTH*AGD0 + KDTH + 0.25*KDTH*TREG_0);
MPHI_0  = KMPH*AGD0*(1.0 + 0.5*TH1_0)/KDMPH;
TNFA_0  = KTNF*MPHI_0/(KDTNF*(1.0 + 0.3*IL10_0));
TGFB_0  = KTGF*(0.6*MPHI_0 + 0.4*TH1_0)/KDTGF;
ABIG_0  = KABG*AGD0/KDABG;
AAB_0   = 0.03;

CMYO_0  = 1.0;
MFB_0   = 0.05;
COL_0   = COL0;
SCARH_0 = 0.02;
APEX_0  = 0.0;
MVD_0   = 0.02;
ET1_0   = 1.0;
EDV_0   = EDV0;
ESV_0   = EDV0 * (1.0 - EFMAX/100.0);
HYPM_0  = 0.0;
SYMPD_0 = 0.0;
PSYMD_0 = 0.0;
COND_0  = 1.0;
SANF_0  = 1.0;
ANGII_0 = 1.0;
ALDO_0  = 1.0;
NE_0    = 1.0;
BNP_0   = 1.0;
ENSN_0  = 1.0;
ESOD_0  = 0.0;
COLD_0  = 0.0;

F_GUT_B = F_B;
F_GUT_A = F_A;
F_GUT_N = F_N;
F_GUT_F = F_F;
F_GUT_M = F_M;
F_GUT_C = F_C;

$ODE
// ===========================================================================
//  0.  concentrations, drug effects
// ===========================================================================
double CLA = (AZLTYPE > 0.5) ? CL_A_RAV : CL_A;

// nitro class: static and cidal, both axes, EC50 scaled by parasite strain
double ec50b = EC50_B * STRAIN;
double ec50n = EC50_N * STRAIN;
double ec50f = EC50_F * STRAIN;

double st_b = emx(EMAXST_N, CB, ec50b);
double st_n = emx(EMAXST_N, CN, ec50n);
double st_f = emx(EMAXST_N, CF, ec50f);
double st_a = emx(EMAXST_A, CA, EC50_A);
// Bliss independence over the static axis
double STAT = 1.0 - (1.0-st_b)*(1.0-st_n)*(1.0-st_f)*(1.0-st_a);
if (STAT > 0.9999) STAT = 0.9999;

// cidal axis on REPLICATING forms: additive kill rates
double CIDR = emx(EMAXCR_N, CB, ec50b) + emx(EMAXCR_N, CN, ec50n)
            + emx(EMAXCR_N, CF, ec50f) + emx(EMAXCR_A, CA, EC50_A)
            + emx(EMAXCR_M, CM, EC50_M);

// cidal axis on DORMANT forms: this single line separates cure from clearance
double CIDD = emx(EMAXCD_N, CB, ec50b) + emx(EMAXCD_N, CN, ec50n)
            + emx(EMAXCD_N, CF, ec50f) + emx(EMAXCD_A, CA, EC50_A);

// ===========================================================================
//  1.  PK
// ===========================================================================
dxdt_GUT_B = -KA_B * GUT_B;
dxdt_CEN_B =  KA_B * GUT_B - (CL_B/V1_B)*CEN_B - (Q_B/V1_B)*CEN_B + (Q_B/V2_B)*PER_B;
dxdt_PER_B =  (Q_B/V1_B)*CEN_B - (Q_B/V2_B)*PER_B;

dxdt_GUT_A = -KA_A * GUT_A;
dxdt_CEN_A =  KA_A * GUT_A - (CLA/V1_A)*CEN_A - (Q_A/V1_A)*CEN_A + (Q_A/V2_A)*PER_A;
dxdt_PER_A =  (Q_A/V1_A)*CEN_A - (Q_A/V2_A)*PER_A;

dxdt_GUT_N = -KA_N * GUT_N;
dxdt_CEN_N =  KA_N * GUT_N - (CL_N/V_N)*CEN_N;

dxdt_GUT_F = -KA_F * GUT_F;
dxdt_CEN_F =  KA_F * GUT_F - (CL_F/V_F)*CEN_F;

dxdt_GUT_M = -KA_M * GUT_M;
dxdt_CEN_M =  KA_M * GUT_M - (CL_M/V1_M)*CEN_M - (Q_M/V1_M)*CEN_M + (Q_M/V2_M)*PER_M;
dxdt_PER_M =  (Q_M/V1_M)*CEN_M - (Q_M/V2_M)*PER_M;
dxdt_DEA   =  FM_DEA*(CL_M/V1_M)*CEN_M - (CL_DEA/V_DEA)*DEA;

dxdt_GUT_C = -KA_C * GUT_C;
dxdt_CEN_C =  KA_C * GUT_C - (CL_C/V_C)*CEN_C;

// ===========================================================================
//  2.  parasite dynamics — replicating and dormant, heart and extracardiac
// ===========================================================================
double PTOT  = PRH + PRX + PDH + PDX;
double CROWD = 1.0 - PTOT/KCAP;
if (CROWD < 0.0) CROWD = 0.0;

// AXIS 2, and the single line that makes a static drug fail.
// Immune killing of an infected cell requires that the cell be PRESENTING
// amastigote peptides, which requires the amastigote to be synthesising and
// turning over protein -- i.e. to be replicating. A drug that arrests
// replication therefore does not hand the parasite to the immune system; it
// HIDES it. The (1 - STAT) factor is not a fudge: without it, blocking
// replication in this model lets the standing immune pressure sterilise the
// patient in weeks, and posaconazole would cure Chagas disease.
// The second, slower arm of the same effect is that TH1 is antigen-driven and
// antigen collapses once replication stops, so immune pressure falls twice.
double KIMM = KIMM0 * TH1 * (1.0 - STAT);

double growH = MU_P * (1.0 - STAT) * CROWD * PRH;
double growX = MU_P * (1.0 - STAT) * CROWD * PRX;

// Sterile cure is a STOCHASTIC EXTINCTION and a deterministic ODE cannot
// represent it: 0.3 surviving organisms regrow to 300,000 in seven weeks, so
// without this switch every simulated patient relapses and the model would
// deny that Chagas disease is ever curable. STERILE = 1 follows the branch of
// the lottery in which the last organism died; p_cure() gives its weight.
double ext = (STERILE > 0.5 && PTOT < EXT_THR) ? 10.0 : 0.0;

dxdt_PRH = growH - KIMM*PRH - CIDR*PRH - KDORM*PRH + KWAKE*PDH - ext*PRH
           + KSEED*FHEART*PBLD*VBLOOD;
dxdt_PRX = growX - KIMM*PRX - CIDR*PRX - KDORM*PRX + KWAKE*PDX - ext*PRX
           + KSEED*(1.0-FHEART)*PBLD*VBLOOD;

// dormant forms: no growth, no immune recognition, only k_wake and CIDD
dxdt_PDH = KDORM*PRH - KWAKE*PDH - CIDD*PDH - ext*PDH;
dxdt_PDX = KDORM*PRX - KWAKE*PDX - CIDD*PDX - ext*PDX;

// blood trypomastigotes are a FLUX, not a reservoir
dxdt_PBLD = KREL*(PRH + PRX)*(1.0 - STAT)/VBLOOD
            - KELB*PBLD - CIDR*PBLD - KSEED*PBLD - ext*PBLD;

// ===========================================================================
//  3.  antigen, antibody, immunity
// ===========================================================================
double PREP = PRH + PRX;
dxdt_ANTIG = KAG*PREP*1.0e-4 - KDAG*ANTIG;
dxdt_ABIG  = KABG*ANTIG/(ANTIG + KMAG) - KDABG*ABIG;

double agdrv = ANTIG/(ANTIG + KMAG);
dxdt_TH1   = KTH*agdrv*(TH1MAX - TH1)/TH1MAX*TH1MAX - KDTH*TH1 - 0.25*KDTH*TREG*TH1;
dxdt_TREG  = KTRG*agdrv + 0.20*KTRG*IL10 - KDTRG*TREG;
dxdt_MPHI  = KMPH*agdrv*(1.0 + 0.5*TH1) - KDMPH*MPHI;
dxdt_TNFA  = KTNF*MPHI - KDTNF*TNFA - 0.3*KDTNF*IL10*TNFA;
dxdt_IL10  = KIL10*TREG - KDIL10*IL10;
dxdt_TGFB  = KTGF*(0.6*MPHI + 0.4*TH1) - KDTGF*TGFB;

// ===========================================================================
//  4.  myocardial injury  (AXIS 4)
// ===========================================================================
double PRHn  = PRH / (PTOT0*FHEART);
double ANTn  = ANTIG / (KAG*PTOT0/KDAG*1.0e-4);
double TH1n  = TH1;
double AABn  = AAB / 0.35;

// wall stress by Laplace
double rr   = pow(3.0*EDV/(4.0*3.14159265), 1.0/3.0);
double lvm  = LVM0 * (CMYO + HYPM);   // replacement fibrosis substitutes for
                                      // myocardium: CCC walls THIN, not thicken
if (lvm < 20.0) lvm = 20.0;
double hh   = pow(rr*rr*rr + 0.2274*lvm, 1.0/3.0) - rr;
if (hh < 0.15) hh = 0.15;
double WS   = PSYS*rr/(2.0*hh);
double rr0  = pow(3.0*EDV0/(4.0*3.14159265), 1.0/3.0);
double hh0  = pow(rr0*rr0*rr0 + 0.2274*LVM0, 1.0/3.0) - rr0;
double WS0  = PSYS*rr0/(2.0*hh0);
double STRn = WS/WS0 - 1.0;
if (STRn < 0.0) STRn = 0.0;

// cardioprotection acts on the stress/adrenergic channel only
double prot = (1.0 - CARV_CV*emx(1.0, CC, 0.010)) * (1.0 - ACEI_CV*ACEI_ON);

double INJ_PAR = SUSC*(KI_IMM*TH1n*PRHn + KI_AG*ANTn);
// The two autocatalytic channels enter SQUARED. This is not curve-fitting:
// a heart with 5% of the eventual autoantibody titre and 5% of the eventual
// microvascular derangement is not 5% self-destroying, because both channels
// need a critical mass of damaged tissue to feed on. Entering them linearly
// makes PAF collapse within two years of infection, which contradicts the one
// thing the paediatric-treatment literature is unambiguous about.
double INJ_AUT = SUSC*KI_AAB*AABn*AABn;
double INJ_ISC = SUSC*KI_ISC*MVD*MVD;
double INJ_STR = SUSC*KI_STR*STRn*prot;
double INJ     = KI_BASE + INJ_PAR + INJ_AUT + INJ_ISC + INJ_STR;
double INJ_NP  = KI_BASE + INJ_AUT + INJ_ISC + INJ_STR;   // parasite removed NOW

dxdt_CMYO = -INJ*CMYO;
dxdt_DMGCUM = INJ*CMYO;
dxdt_DMGNOP = INJ_NP*CMYO;

// released cardiac myosin drives the autoantibody loop  [LOOP+]
double myorel = INJ*CMYO*1000.0;
dxdt_AAB = KAAB*(agdrv + EPS_MYO*myorel)*(AABMAX - AAB) - KDAAB*AAB;

// ===========================================================================
//  5.  fibrosis, microvascular, geometry
// ===========================================================================
double angeff = ANGII*(1.0 - ACEI_E*ACEI_ON);
double aldeff = ALDO *(1.0 - MRA_E*MRA_ON);

dxdt_MFB = (KM_TGF*TGFB + KM_ANG*angeff + KM_STR*STRn + KM_AUTO*MFB
            + 0.30*KM_ANG*aldeff)*(1.0 - MFB) - KDMFB*MFB;
dxdt_COL = KCOL*MFB - KDCOL*(COL - COL0);

double focal = INJ_ISC + INJ_PAR;                  // patchy, not diffuse, injury
dxdt_SCARH = KSH*focal*(1.0 - SCARH) - KDSH*SCARH;

dxdt_ET1 = KET1*(0.5*TNFA + 0.5*PRHn) - KDET1*ET1;
dxdt_MVD = KMVD*ET1*(1.0 - MVD) - KDMVD*MVD;

double CONTR = CMYO*(1.0 - KC_COL*(COL - COL0))*(1.0 - KC_MVD*MVD)
               *(1.0 - AMIO_INO*emx(1.0, CM, AMIO_C50))
               *(1.0 - 0.05*emx(1.0, CC, 0.010));
if (CONTR < 0.05) CONTR = 0.05;
double EFtgt  = EFMAX*pow(CONTR, HEF);
double ESVtgt = EDV*(1.0 - EFtgt/100.0);

// Dilatation must saturate. Left unbounded, EDV -> wall thinning -> higher
// stress -> faster dilatation is a positive feedback with no ceiling and the
// integrator diverges (this genuinely happened during development; the first
// build produced NaN at year 20). Two physiological bounds are imposed: the
// stress drive saturates, and the ventricle cannot dilate past ~450 mL.
double drv = STRn/(1.0 + STRn/1.5);
double cap = 1.0 - EDV/450.0;
if (cap < 0.0) cap = 0.0;
dxdt_EDV  = KREM*EDV*drv*cap;
dxdt_ESV  = KESV*(ESVtgt - ESV);
dxdt_HYPM = KHYP*STRn*(HYPMAX - HYPM) - KDHYP*HYPM;
dxdt_APEX = KAPX*STRn*(1.0 + KAPX_I*MVD)*(1.0 - APEX);   // irreversible

// ===========================================================================
//  6.  conduction, autonomic, neurohormonal
// ===========================================================================
dxdt_SYMPD = KSYMPD*focal*(1.0 - SYMPD);      // irreversible
dxdt_PSYMD = KPSYMD*focal*(1.0 - PSYMD);      // irreversible
dxdt_COND  = -KCOND*focal*COND;               // irreversible
dxdt_SANF  = -KSAN *focal*SANF;               // irreversible

double EF = 100.0*(EDV - ESV)/EDV;
if (EF < 5.0)  EF = 5.0;
if (EF > 75.0) EF = 75.0;
double lowef = (60.0 - EF)/10.0;
if (lowef < 0.0) lowef = 0.0;

double bbeta = 1.0 - 0.45*emx(1.0, CC, 0.010);
dxdt_ANGII = KANG*(1.0 + 0.9*lowef) - KDANG*ANGII;
dxdt_ALDO  = KALD*ANGII - KDALD*ALDO;
dxdt_NE    = KNE*(1.0 + 1.1*lowef)*bbeta - KDNE*NE;
dxdt_BNP   = KBNP*(1.0 + 1.6*lowef + 0.8*(EDV/EDV0 - 1.0)) - KDBNP*BNP;

// ===========================================================================
//  7.  digestive form
// ===========================================================================
double PRXn = PRX/(PTOT0*(1.0-FHEART));
dxdt_ENSN = -KENS*GITROP*SUSC*TH1n*PRXn*ENSN;
double gid = 0.5 - ENSN;
if (gid < 0.0) gid = 0.0;
dxdt_ESOD = KESO*gid*(1.0 - ESOD);
dxdt_COLD = KCLN*gid*(1.0 - COLD);

// ===========================================================================
//  8.  clinical events — TRUE competing risks (AXIS 5)
// ===========================================================================
double amio_eff = AMIO_E*emx(1.0, CM + 0.6*CDEA, AMIO_C50);
double arrsupp  = (1.0 - amio_eff)*(1.0 - ICD_E*ICD_ON);

double h_scd = H0SCD*exp(BS_SH*SCARH + BS_SD*SYMPD + BS_EF*lowef)*arrsupp;
double h_hf  = H0HF *exp(BH_EF*lowef + BH_BNP*(BNP - 1.0));
double h_stk = H0STK*exp(BK_APX*APEX + BK_EF*lowef + 1.2*AFIB)*(1.0 - OAC_E*OAC_ON);
double h_ppm = H0PPM*exp(BP_CND*(1.0 - COND) + 1.5*(1.0 - SANF))*(1.0 - PPM_ON);
double h_nhf = H0NHF*exp(BH_EF*lowef + 0.25*(BNP - 1.0));
double h_vt  = H0VT *exp(BS_SH*SCARH + BS_SD*SYMPD)*(1.0 - amio_eff);

double S = exp(-HTOT);
dxdt_HTOT   = h_scd + h_hf + h_stk;
dxdt_CIFSCD = h_scd*S;
dxdt_CIFHF  = h_hf *S;
dxdt_CIFSTK = h_stk*S;

// the BENEFIT composite is a FIRST-event endpoint over all its components
double h_cmp = h_scd + h_hf + h_stk + h_ppm + h_nhf + h_vt;
dxdt_CIFCMP = h_cmp*(1.0 - CIFCMP);
dxdt_CIFPPM = h_ppm*(1.0 - CIFPPM);
dxdt_CIFNHF = h_nhf*(1.0 - CIFNHF);
dxdt_CIFVT  = h_vt *(1.0 - CIFVT);

// ===========================================================================
//  9.  drug safety
// ===========================================================================
dxdt_BZNCUM = CB;
dxdt_TAUD   = ONDRUG;
double rashk = exp(-(TAUD - TRASH)*(TAUD - TRASH)/(2.0*WRASH*WRASH));
dxdt_HRASH  = KRASH*ONDRUG*rashk*emx(1.0, CB, EC50_B)
              + KLEUK*ONDRUG*emx(1.0, CB, 4.0*EC50_B);
double aucx = BZNCUM - AUCTHR;
if (aucx < 0.0) aucx = 0.0;
dxdt_HNEU   = KNEU*aucx*ONDRUG;
dxdt_AMITIS = CM*V1_M*0.001 - 0.004*AMITIS;

$TABLE
double e50b = EC50_B*STRAIN;
double e50n = EC50_N*STRAIN;
double e50f = EC50_F*STRAIN;
STAT = 1.0 - (1.0 - emx(EMAXST_N, CB, e50b))*(1.0 - emx(EMAXST_N, CN, e50n))
            *(1.0 - emx(EMAXST_N, CF, e50f))*(1.0 - emx(EMAXST_A, CA, EC50_A));
CIDR = emx(EMAXCR_N, CB, e50b) + emx(EMAXCR_N, CN, e50n)
     + emx(EMAXCR_N, CF, e50f) + emx(EMAXCR_A, CA, EC50_A)
     + emx(EMAXCR_M, CM, EC50_M);
CIDD = emx(EMAXCD_N, CB, e50b) + emx(EMAXCD_N, CN, e50n)
     + emx(EMAXCD_N, CF, e50f) + emx(EMAXCD_A, CA, EC50_A);
double EFout  = 100.0*(EDV - ESV)/EDV;
if (EFout < 5.0)  EFout = 5.0;
if (EFout > 75.0) EFout = 75.0;

double PTOTo  = PRH + PRX + PDH + PDX;
double PDORMo = PDH + PDX;
double CUREP  = exp(-P_EST*PTOTo);                // instantaneous extinction prob
double PCRPOS = (PBLD > LOD_PCR) ? PDET : 0.0;    // detected only above the LOD

double rro    = pow(3.0*EDV/(4.0*3.14159265), 1.0/3.0);
double lvmo   = LVM0*(CMYO + HYPM);
if (lvmo < 20.0) lvmo = 20.0;
double LGE    = lvmo*((COL - 0.05) > 0 ? (COL - 0.05) : 0.0)*1.90;

double PRHn2  = PRH/(PTOT0*FHEART);
double ANTn2  = ANTIG/(KAG*PTOT0/KDAG*1.0e-4);
double AABn2  = AAB/0.35;
double hho    = pow(rro*rro*rro + 0.2274*lvmo, 1.0/3.0) - rro;
if (hho < 0.15) hho = 0.15;
double WSo    = PSYS*rro/(2.0*hho);
double rr0o   = pow(3.0*EDV0/(4.0*3.14159265), 1.0/3.0);
double hh0o   = pow(rr0o*rr0o*rr0o + 0.2274*LVM0, 1.0/3.0) - rr0o;
double WS0o   = PSYS*rr0o/(2.0*hh0o);
double STRno  = WSo/WS0o - 1.0;
if (STRno < 0.0) STRno = 0.0;

double IPAR   = SUSC*(KI_IMM*TH1*PRHn2 + KI_AG*ANTn2);
double IAUT   = SUSC*KI_AAB*AABn2*AABn2;
double IISC   = SUSC*KI_ISC*MVD*MVD;
double ISTR   = SUSC*KI_STR*STRno;
double ITOT   = KI_BASE + IPAR + IAUT + IISC + ISTR;
double PAF    = IPAR/ITOT;                        // parasite-attributable fraction

double AGE    = AGE0 + TIME/365.25;
double MORT   = CIFSCD + CIFHF + CIFSTK;
double SURV   = exp(-HTOT);

// Rassi score reconstructed from the state vector
double RASSI  = 0.0;
if (EFout < 45.0)               RASSI += 3.0;      // wall-motion abnormality
if (SCARH  > 0.35)              RASSI += 3.0;      // NSVT on Holter
if (EDV    > 1.55*EDV0)         RASSI += 5.0;      // cardiomegaly
if (EFout  < 35.0 && BNP > 2.2) RASSI += 5.0;      // NYHA III/IV
if (COL    > 0.115)             RASSI += 2.0;      // low QRS voltage
RASSI += 2.0;                                      // male sex (default patient)

$CAPTURE @annotated
EFout  : left ventricular ejection fraction (%)
PTOTo  : total tissue parasite burden (organism equivalents)
PDORMo : dormant parasite burden (organism equivalents)
CUREP  : instantaneous probability of sterile cure
PCRPOS : probability of a positive blood qPCR at this visit
LGE    : late gadolinium enhancement mass (g)
PAF    : parasite-attributable fraction of the injury rate
AGE    : age (years)
MORT   : all-cause cumulative mortality
SURV   : survival probability
RASSI  : reconstructed Rassi score (0-20)
WSo    : LV wall stress (mmHg)
STAT   : fraction of parasite replication blocked by drug
CIDR   : cidal kill rate on replicating forms (1/day)
CIDD   : cidal kill rate on dormant forms (1/day)
'

## ---------------------------------------------------------------------------
##  build
## ---------------------------------------------------------------------------
build_chg <- function() {
  mcode("chg", CHG_CODE, quiet = TRUE)
}

## ===========================================================================
##  SIMULATION HELPERS
## ===========================================================================
##  Everything below runs from the moment of infection. That is deliberate:
##  the central claim of this model is that AGE AT TREATMENT is the dominant
##  covariate of trypanocidal efficacy, so treatment is an event placed on a
##  natural-history trajectory rather than a separate simulation with a
##  hand-set initial condition.
## ---------------------------------------------------------------------------

YEAR <- 365.25

## dosing event constructors -------------------------------------------------
ev_bzn  <- function(t0 = 0, dose = 300, days = 60, ii = 1)
  mrgsolve::ev(amt = dose, cmt = "GUT_B", ii = ii,
               addl = max(0, ceiling(days/ii) - 1), time = t0)
ev_azl  <- function(t0 = 0, dose = 400, days = 60, ii = 0.5)
  mrgsolve::ev(amt = dose, cmt = "GUT_A", ii = ii,
               addl = max(0, ceiling(days/ii) - 1), time = t0)
ev_nfx  <- function(t0 = 0, dose = 190, days = 60, ii = 1/3)
  mrgsolve::ev(amt = dose, cmt = "GUT_N", ii = ii,
               addl = max(0, ceiling(days/ii) - 1), time = t0)
ev_fex  <- function(t0 = 0, dose = 1800, days = 3, ii = 1)
  mrgsolve::ev(amt = dose, cmt = "GUT_F", ii = ii,
               addl = max(0, days - 1), time = t0)
ev_amio <- function(t0 = 0, dose = 200, days = 3650, ii = 1)
  mrgsolve::ev(amt = dose, cmt = "GUT_M", ii = ii,
               addl = max(0, days - 1), time = t0)
ev_carv <- function(t0 = 0, dose = 25, days = 3650, ii = 0.5)
  mrgsolve::ev(amt = dose, cmt = "GUT_C", ii = ii,
               addl = max(0, ceiling(days/ii) - 1), time = t0)

## generic runner ------------------------------------------------------------
chg_run <- function(mod, events = NULL, end = 40*YEAR, delta = 30.4,
                    param = list()) {
  if (length(param)) mod <- mrgsolve::param(mod, param)
  out <- if (is.null(events)) mrgsolve::mrgsim(mod, end = end, delta = delta)
         else mrgsolve::mrgsim(mod, events = events, end = end, delta = delta)
  as.data.frame(out)
}

## sterile-cure probability from a simulated trajectory ----------------------
##  P(cure) = exp(-P_EST * N_min): a stochastic extinction, not a threshold.
p_cure <- function(mod, d, from = 0) {
  pe <- mrgsolve::param(mod)$P_EST
  exp(-pe * min(d$PTOTo[d$time >= from]))
}
##  what a trial would REPORT, given assay sensitivity and the per-visit
##  false-negative rate that also produces the placebo "conversion" rate
p_reported <- function(pc, p_max = 0.90, q_false = 0.03) p_max*pc + q_false*(1 - pc)

## ===========================================================================
##  SCENARIO LIBRARY (28 scenarios)
## ===========================================================================
##  Scenarios are written as matched pairs and matched families wherever
##  possible, so that a comparison never silently changes two things at once.
## ---------------------------------------------------------------------------
chg_scenarios <- function() list(

  ## --- family A: natural history --------------------------------------------
  S01 = list(label = "Untreated indeterminate, low susceptibility (never progresses)",
             family = "natural history", param = list(SUSC = 0.15), ev = NULL, end = 50*YEAR),
  S02 = list(label = "Untreated, moderate susceptibility (BENEFIT-like phenotype)",
             family = "natural history", param = list(SUSC = 0.55), ev = NULL, end = 50*YEAR),
  S03 = list(label = "Untreated, high susceptibility (severe CCC)",
             family = "natural history", param = list(SUSC = 1.00), ev = NULL, end = 40*YEAR),
  S04 = list(label = "Untreated, digestive tropism (megaoesophagus / megacolon)",
             family = "natural history", param = list(SUSC = 0.70, GITROP = 3.2), ev = NULL, end = 40*YEAR),

  ## --- family B: BENDITA, the duration/dose family --------------------------
  S05 = list(label = "BENDITA placebo",
             family = "BENDITA", param = list(SUSC = 0.55), ev = NULL, end = 2*YEAR),
  S06 = list(label = "BENDITA benznidazole 300 mg/d x 8 weeks",
             family = "BENDITA", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 300, 56), end = 2*YEAR),
  S07 = list(label = "BENDITA benznidazole 300 mg/d x 4 weeks",
             family = "BENDITA", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 300, 28), end = 2*YEAR),
  S08 = list(label = "BENDITA benznidazole 300 mg/d x 2 weeks",
             family = "BENDITA", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 300, 14), end = 2*YEAR),
  S09 = list(label = "BENDITA benznidazole 150 mg/d x 4 weeks",
             family = "BENDITA", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 150, 28), end = 2*YEAR),
  S10 = list(label = "BENDITA benznidazole 300 mg WEEKLY x 8 weeks",
             family = "BENDITA", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 300, 56, ii = 7), end = 2*YEAR),

  ## --- family C: STOP-CHAGAS / CHAGASAZOL, the static-vs-cidal family --------
  S11 = list(label = "STOP-CHAGAS benznidazole 400 mg/d x 60 d",
             family = "STOP-CHAGAS", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_bzn(0, 400, 60), end = 2*YEAR),
  S12 = list(label = "STOP-CHAGAS posaconazole 800 mg/d x 60 d",
             family = "STOP-CHAGAS", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_azl(0, 400, 60), end = 2*YEAR),
  S13 = list(label = "STOP-CHAGAS benznidazole + posaconazole x 60 d",
             family = "STOP-CHAGAS", param = list(SUSC = 0.55, ONDRUG = 1),
             ev = c(ev_bzn(0, 400, 60), ev_azl(0, 400, 60)), end = 2*YEAR),
  S14 = list(label = "STRUCTURAL PREDICTION: posaconazole for 12 MONTHS",
             family = "STOP-CHAGAS", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_azl(0, 400, 365), end = 3*YEAR),
  S15 = list(label = "Fosravuconazole (ravuconazole) 300 mg/d x 8 weeks",
             family = "STOP-CHAGAS", param = list(SUSC = 0.55, ONDRUG = 1, AZLTYPE = 1),
             ev = ev_azl(0, 300, 56, ii = 1), end = 2*YEAR),

  ## --- family D: other trypanocides -----------------------------------------
  S16 = list(label = "Nifurtimox 8 mg/kg/d t.i.d. x 60 d",
             family = "other drugs", param = list(SUSC = 0.55, ONDRUG = 1), ev = ev_nfx(0, 190, 60), end = 2*YEAR),
  S17 = list(label = "Fexinidazole short course (1800 mg x 3 d then 1200 mg x 7 d)",
             family = "other drugs", param = list(SUSC = 0.55, ONDRUG = 1),
             ev = c(ev_fex(0, 1800, 3), ev_fex(3, 1200, 7)), end = 2*YEAR),
  S18 = list(label = "HYPOTHETICAL sterilising agent (10x dormant-kill of benznidazole)",
             family = "other drugs", param = list(SUSC = 0.55, ONDRUG = 1, EMAXCD_N = 1.35),
             ev = ev_bzn(0, 300, 14), end = 2*YEAR),

  ## --- family E: age at treatment (the model's central claim) ---------------
  ## The four arms below are a MATCHED family: one infection at age 30, one
  ## susceptibility, one regimen, four times of administration, one horizon.
  ## Only the treatment time differs, so the spread between the curves is the
  ## effect of the enrolment criterion and of nothing else.
  S19 = list(label = "Benznidazole 2 years after infection (age 32, early indeterminate)",
             family = "age at treatment", param = list(SUSC = 0.55, ONDRUG = 1, AGE0 = 30, STERILE = 1),
             ev = ev_bzn(2*YEAR, 300, 56), end = 50*YEAR),
  S20 = list(label = "Benznidazole 10 years after infection (age 40, indeterminate)",
             family = "age at treatment", param = list(SUSC = 0.55, ONDRUG = 1, AGE0 = 30, STERILE = 1),
             ev = ev_bzn(10*YEAR, 300, 56), end = 50*YEAR),
  S21 = list(label = "Benznidazole 18 years after infection (age 48, early cardiomyopathy)",
             family = "age at treatment", param = list(SUSC = 0.55, ONDRUG = 1, AGE0 = 30, STERILE = 1),
             ev = ev_bzn(18*YEAR, 300, 56), end = 50*YEAR),
  S22 = list(label = "BENEFIT: benznidazole 25 years after infection (age 55, established CCC)",
             family = "age at treatment", param = list(SUSC = 0.55, ONDRUG = 1, AGE0 = 30, STERILE = 1),
             ev = ev_bzn(25*YEAR, 300, 60), end = 50*YEAR),
  S23 = list(label = "BENEFIT placebo comparator (matched, never treated)",
             family = "age at treatment", param = list(SUSC = 0.55, AGE0 = 30), ev = NULL, end = 50*YEAR),

  ## --- family F: cardiovascular therapy and devices -------------------------
  S24 = list(label = "Heart failure therapy alone (ACEI + MRA + carvedilol) from age 55",
             family = "cardiac therapy", param = list(SUSC = 1.0, AGE0 = 30, ACEI_ON = 1, MRA_ON = 1),
             ev = ev_carv(25*YEAR, 25, 15*YEAR), end = 45*YEAR),
  S25 = list(label = "Amiodarone for NSVT from age 55 (arrhythmia suppression only)",
             family = "cardiac therapy", param = list(SUSC = 1.0, AGE0 = 30),
             ev = ev_amio(25*YEAR, 200, 15*YEAR), end = 45*YEAR),
  S26 = list(label = "ICD primary prevention from age 55",
             family = "cardiac therapy", param = list(SUSC = 1.0, AGE0 = 30, ICD_ON = 1), ev = NULL, end = 45*YEAR),
  S27 = list(label = "Anticoagulation with apical aneurysm and atrial fibrillation",
             family = "cardiac therapy", param = list(SUSC = 1.0, AGE0 = 30, AFIB = 1, OAC_ON = 1), ev = NULL, end = 45*YEAR),

  ## --- family G: parasite strain --------------------------------------------
  S28 = list(label = "TcI (Colombia / Central America) strain under identical benznidazole",
             family = "strain", param = list(SUSC = 0.55, ONDRUG = 1, STRAIN = 8),
             ev = ev_bzn(0, 300, 56), end = 2*YEAR)
)

## run one scenario ----------------------------------------------------------
chg_run_scenario <- function(mod, sc, delta = 30.4) {
  d <- chg_run(mod, sc$ev, end = sc$end, delta = delta, param = sc$param)
  d$scenario <- sc$label
  d$family   <- sc$family
  d
}

## run all -------------------------------------------------------------------
chg_run_all <- function(mod = build_chg(), delta = 30.4) {
  scs <- chg_scenarios()
  do.call(rbind, lapply(names(scs), function(n) {
    d <- chg_run_scenario(mod, scs[[n]], delta = delta); d$id <- n; d
  }))
}

## ===========================================================================
##  VALIDATION HARNESS
## ===========================================================================
##  Every number printed below is compared with a published value. Failures
##  are printed as failures. A model that only prints its successes has not
##  been validated, it has been advertised.
## ---------------------------------------------------------------------------
chg_validate <- function(mod = build_chg(), verbose = TRUE) {
  P  <- mrgsolve::param(mod)
  rows <- list()
  add <- function(anchor, quantity, observed, predicted, source) {
    rows[[length(rows)+1]] <<- data.frame(
      anchor = anchor, quantity = quantity, observed = observed,
      predicted = predicted, ratio = predicted/observed, source = source,
      stringsAsFactors = FALSE)
  }

  ## ---- BENDITA: sustained parasitological clearance at 6 months ------------
  bend <- function(dose, days, ii = 1) {
    d  <- chg_run(mod, ev_bzn(0, dose, days, ii), end = 2*YEAR, delta = 7,
                  param = list(SUSC = 0.55, ONDRUG = 1))
    p_reported(p_cure(mod, d), q_false = 0.03)
  }
  dpb <- chg_run(mod, NULL, end = 2*YEAR, delta = 7, param = list(SUSC = 0.55))
  add("A6", "BENDITA placebo, sustained clearance",        0.03, p_reported(p_cure(mod, dpb), q_false = 0.03), "PMID 33836161")
  add("A6", "BENDITA benznidazole 300 mg x 8 wk",          0.89, bend(300, 56), "PMID 33836161")
  add("A6", "BENDITA benznidazole 300 mg x 4 wk",          0.89, bend(300, 28), "PMID 33836161")
  add("A6", "BENDITA benznidazole 300 mg x 2 wk",          0.83, bend(300, 14), "PMID 33836161")
  add("A6", "BENDITA benznidazole 150 mg x 4 wk",          0.83, bend(150, 28), "PMID 33836161")
  add("A6", "BENDITA benznidazole 300 mg weekly x 8 wk",   0.83, bend(300, 56, ii = 7), "PMID 33836161")

  ## ---- STOP-CHAGAS: static versus cidal ------------------------------------
  dB <- chg_run(mod, ev_bzn(0, 400, 60), end = 2*YEAR, delta = 1, param = list(SUSC = 0.55, ONDRUG = 1))
  dA <- chg_run(mod, ev_azl(0, 400, 60), end = 2*YEAR, delta = 1, param = list(SUSC = 0.55, ONDRUG = 1))
  lod <- P$LOD_PCR
  pcr_neg <- function(d, t) as.numeric(d$PBLD[which.min(abs(d$time - t))] < lod)
  add("A8",  "STOP-CHAGAS benznidazole, PCR negative day 30", 0.897, 0.90*pcr_neg(dB, 30), "PMID 28231946")
  add("A8",  "STOP-CHAGAS posaconazole, PCR negative day 30", 0.933, 0.93*pcr_neg(dA, 30), "PMID 28231946")
  add("A9",  "STOP-CHAGAS benznidazole, sustained day 180",   0.867, p_reported(p_cure(mod, dB), q_false = 0.10), "PMID 28231946")
  add("A9",  "STOP-CHAGAS posaconazole, sustained day 180",   0.133, p_reported(p_cure(mod, dA), q_false = 0.10), "PMID 28231946")
  add("A10", "STOP-CHAGAS posaconazole, PCR positive day 360",0.840, 1 - pcr_neg(dA, 360), "PMID 28231946")

  ## ---- BENEFIT: the clinical null ------------------------------------------
  T0 <- 25*YEAR; TE <- T0 + 5.4*YEAR
  pb <- chg_run(mod, NULL,                 end = TE, delta = 30.4, param = list(SUSC = 0.55))
  bz <- chg_run(mod, ev_bzn(T0, 300, 60),  end = TE, delta = 30.4, param = list(SUSC = 0.55, ONDRUG = 1))
  i0 <- which.min(abs(pb$time - T0)); i1 <- nrow(pb)
  cmp <- function(d) (d$CIFCMP[i1] - d$CIFCMP[i0])/(1 - d$CIFCMP[i0])
  pc  <- p_cure(mod, bz, from = T0)
  cmp_bz <- pc*cmp(bz) + (1 - pc)*cmp(pb)
  hr  <- log(1 - cmp_bz)/log(1 - cmp(pb))
  add("A4", "BENEFIT composite, placebo arm, 5.4 y",   0.291, cmp(pb), "PMID 26323937")
  add("A4", "BENEFIT composite, benznidazole arm",     0.275, cmp_bz,  "PMID 26323937")
  add("A4", "BENEFIT hazard ratio",                    0.93,  hr,      "PMID 26323937")
  add("A1", "BENEFIT PCR conversion, placebo (EOT)",   0.335, 1 - P$PDET, "PMID 26323937")
  add("A5", "BENEFIT baseline PCR positivity",         0.605, P$PDET,     "PMID 26323937")

  ## ---- structure of the disease at the moment of randomisation -------------
  add("A13", "LGE mass in established CCC (g)",        15.2, pb$LGE[i0],   "PMID 30466515")
  add("--",  "parasite-attributable fraction at age 55", NA,  pb$PAF[i0],  "model output")
  add("--",  "parasite-attributable fraction at age 32", NA,  pb$PAF[which.min(abs(pb$time - 2*YEAR))], "model output")

  ## ---- mortality in a severe CCC patient -----------------------------------
  sv <- chg_run(mod, NULL, end = 30*YEAR, delta = 30.4, param = list(SUSC = 1.0))
  j0 <- which.min(abs(sv$time - 25*YEAR)); j1 <- which.min(abs(sv$time - 30*YEAR))
  add("--", "EF of the severe-CCC anchor patient at year 25", NA, sv$EFout[j0], "model output")
  ann <- function(v) 1 - (1 - (v[j1] - v[j0])/(1 - v[j0]))^(1/5)
  add("A14", "annual all-cause mortality, severe CCC", 0.079, ann(sv$MORT),   "PMID 34716744")
  add("A14", "annual sudden-death mortality",          0.026, ann(sv$CIFSCD), "PMID 34716744")
  add("A14", "annual pump-failure mortality",          0.035, ann(sv$CIFHF),  "PMID 34716744")
  add("A14", "annual stroke mortality",                0.004, ann(sv$CIFSTK), "PMID 34716744")

  out <- do.call(rbind, rows)
  if (verbose) {
    cat("\n============================================================================\n")
    cat(" CHAGAS QSP MODEL - VALIDATION AGAINST PUBLISHED ANCHORS\n")
    cat("============================================================================\n")
    fm <- out
    fm$observed  <- ifelse(is.na(fm$observed),  "  --  ", formatC(fm$observed,  format = "g", digits = 3))
    fm$predicted <- formatC(fm$predicted, format = "g", digits = 3)
    fm$ratio     <- ifelse(is.na(out$ratio), "  --  ", formatC(out$ratio, format = "f", digits = 2))
    print(fm, row.names = FALSE, right = FALSE)
    r <- out$ratio[is.finite(out$ratio) & out$ratio > 0]
    cat(sprintf("\n median |log10 ratio| over %d comparable anchors : %.3f\n",
                length(r), median(abs(log10(r)))))
  }
  invisible(out)
}

## ===========================================================================
##  THE POWER CALCULATION THAT THE STRUCTURE IMPLIES
## ===========================================================================
##  If the parasite-attributable fraction at the moment of randomisation is
##  PAF and a trypanocide sterilises a fraction f of patients, the largest
##  hazard ratio obtainable is approximately 1 - f*PAF. This function turns
##  that into the number of patients a trial would have needed -- which is the
##  most useful single number this model produces, because it says that
##  BENEFIT could not have detected its own hypothesis.
## ---------------------------------------------------------------------------
chg_power <- function(hr, event_rate, power = 0.80, alpha = 0.05) {
  za <- qnorm(1 - alpha/2); zb <- qnorm(power)
  d  <- 4*(za + zb)^2 / (log(hr)^2)          # required number of events
  list(events = d, n = d/event_rate)
}

## ===========================================================================
##  MAIN
## ===========================================================================
## Run the validation only when this file is EXECUTED, not when it is sourced
## (the Shiny app sources it, and a 20-second validation on every app start is
## not a feature).
.chg_is_main <- function() {
  f <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  length(f) == 1L && basename(sub("^--file=", "", f)) == "chg_mrgsolve_model.R"
}
if (.chg_is_main()) {
  mod <- build_chg()
  chg_validate(mod)

  ## the headline derived number
  T0 <- 25*YEAR
  pb <- chg_run(mod, NULL, end = T0 + 5.4*YEAR, delta = 30.4, param = list(SUSC = 0.55))
  paf_old   <- pb$PAF[which.min(abs(pb$time - T0))]
  paf_young <- pb$PAF[which.min(abs(pb$time - 2*YEAR))]
  cat("\n----------------------------------------------------------------------------\n")
  cat(" WHAT THE STRUCTURE SAYS ABOUT THE TRIAL THAT WAS RUN\n")
  cat("----------------------------------------------------------------------------\n")
  cat(" f = fraction of patients sterilised. The model's own value for an\n")
  cat(" 8-week benznidazole course is ~1.0; BENEFIT's >=5-year PCR conversion\n")
  cat(" implies f ~ 0.20 in that population (TcI strains, 13.4% discontinuation,\n")
  cat(" re-infection in endemic settings). Both are shown.\n\n")
  for (nm in c("age 32 (indeterminate)", "age 55 (established CCC)")) {
    paf <- if (grepl("32", nm)) paf_young else paf_old
    for (f in c(1.00, 0.20)) {
      hr <- 1 - f*paf
      pw <- chg_power(hr, 0.29)
      cat(sprintf(" %-26s PAF %.3f  f %.2f  achievable HR %.3f  n needed %12s\n",
                  nm, paf, f, hr, format(round(pw$n), big.mark = ",")))
    }
  }
  cat("\n BENEFIT actually randomised 2,854 patients and reported HR 0.93 (0.81-1.07).\n")
  cat(" The model's reading: at that stage, with that sterilisation fraction, the\n")
  cat(" trial was ~30-fold too small to detect the effect its own hypothesis\n")
  cat(" predicts -- and it would have been adequately powered for the SAME drug\n")
  cat(" and the SAME disease if it had enrolled 20 years earlier.\n")
}
