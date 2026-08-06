## =============================================================================
##  INVASIVE PULMONARY ASPERGILLOSIS (IPA) — QSP / mrgsolve MODEL
## =============================================================================
##  53 ODEs · 5 antifungals · 2 co-medications · 23 treatment scenarios
##
##  THE ONE IDEA
##  ------------
##      dB/dt = k_grow*(1 - I_azole)*B  -  (k_host*N_eff + k_cidal)*B
##
##  A triazole enters ONLY the first term. It multiplies growth; it never adds
##  to removal. Therefore, if N_eff = 0 (profound neutropenia) and no polyene
##  is given, dB/dt >= k_grow*(1 - Imax)*B > 0 for every azole exposure that is
##  physically achievable, because Imax < 1. An azole cannot clear invasive
##  aspergillosis in a host with no neutrophils — it can only buy time until
##  the marrow recovers. Everything else in this file is a consequence:
##
##    * neutrophil recovery outranks drug choice as a mortality lever
##    * there is a critical time-to-treatment: a bifurcation, not a gradient.
##      I expected angioinvasion-limited drug delivery to be the cause and
##      tested it by deleting that limitation; the cliff did not move. It is
##      set instead by how large the burden already is when a growth-rate
##      inhibitor takes over. See README, "a hypothesis of my own that did
##      not survive"
##    * galactomannan and fungal burden disagree in week 1, because GM is
##      released by the GROWTH flux and by the KILLING flux, not by the stock
##    * voriconazole's therapeutic window is narrow because of nonlinear PK,
##      not because of a steep PD curve
##
##  VERIFICATION
##  ------------
##  Every equation below is implemented a second time, independently, in
##  `ipa_reference_model.py` (Python/scipy). Numbers quoted in README.md come
##  from that file. The dual implementation found ten defects during
##  development; they are listed in README.md and each is annotated in the
##  Python source at the line where it lived.
##
##  UNITS  time h · concentration mg/L · fungal burden CFU-equivalents (CFUe)
##         ANC 10^9/L · cytokines pg/mL · GM optical density index
##
##  Requires: mrgsolve, dplyr, ggplot2, tidyr
## =============================================================================

library(mrgsolve)
suppressMessages({library(dplyr); library(tidyr)})

code <- '
$PROB IPA — invasive pulmonary aspergillosis QSP model (53 ODEs)

$PLUGIN base

$PARAM @annotated
// ---------------------------- voriconazole ---------------------------------
// 2-cpt; elimination = saturable CYP2C19 path + linear CYP3A4/FMO3 path.
// Calibrated to 4 mg/kg IV q12h -> Cmin ~1.1, AUC24 ~48 mg*h/L in a normal
// metaboliser (Purkins 2002; Pascual 2008; Dolton 2012).
VRC_ka    : 1.10  : voriconazole absorption rate (1/h)
VRC_F     : 0.96  : voriconazole oral bioavailability
VRC_Vc    : 50.0  : voriconazole central volume (L)
VRC_Vp    : 150.0 : voriconazole peripheral volume (L)
VRC_Q     : 25.0  : voriconazole intercompartmental clearance (L/h)
VRC_Vmax  : 47.8  : CYP2C19 Vmax at reference genotype (mg/h)
VRC_Km    : 3.00  : CYP2C19 Km (mg/L)
VRC_CLlin : 2.50  : CYP3A4/FMO3 linear clearance (L/h)
VRC_fu    : 0.42  : voriconazole unbound fraction
VRC_Kpelf : 7.10  : ELF-to-plasma partition
VRC_keqelf: 0.60  : ELF equilibration rate (1/h)
VRC_Kpbr  : 0.50  : brain-to-plasma partition
VRC_keqbr : 0.10  : brain equilibration rate (1/h)
VRC_iEmax : 0.62  : max CRP-driven CYP2C19 phenoconversion
VRC_iEC50 : 105.0 : CRP giving half-maximal phenoconversion (mg/L)
CYP2C19   : 1.00  : CYP2C19 activity multiplier (UM 1.85 / NM 1.0 / PM 0.30)

// ---------------------------- isavuconazole --------------------------------
ISA_ka    : 0.90  : isavuconazole absorption rate (1/h)
ISA_F     : 0.98  : isavuconazole bioavailability
ISA_Vc    : 90.0  : isavuconazole central volume (L)
ISA_Vp    : 360.0 : isavuconazole peripheral volume (L)
ISA_Q     : 8.00  : isavuconazole intercompartmental clearance (L/h)
ISA_CL    : 1.75  : isavuconazole clearance (L/h)
ISA_Kpelf : 1.20  : isavuconazole ELF partition
ISA_keqelf: 0.50  : isavuconazole ELF equilibration (1/h)
ISA_Kpbr  : 0.11  : isavuconazole brain partition

// ---------------------------- posaconazole ---------------------------------
POS_ka    : 0.30  : posaconazole absorption rate (1/h)
POS_F     : 0.54  : posaconazole DR-tablet bioavailability
POS_Vc    : 190.0 : posaconazole central volume (L)
POS_Vp    : 290.0 : posaconazole peripheral volume (L)
POS_Q     : 6.00  : posaconazole intercompartmental clearance (L/h)
POS_CL    : 4.50  : posaconazole clearance (L/h)
POS_Kpelf : 1.60  : posaconazole ELF partition
POS_keqelf: 0.40  : posaconazole ELF equilibration (1/h)
POS_Kpbr  : 0.03  : posaconazole brain partition

// ------------------------ liposomal amphotericin B -------------------------
AMB_Vc    : 2.00  : L-AmB central volume (L)
AMB_Vp    : 10.0  : L-AmB peripheral volume (L)
AMB_Q     : 0.90  : L-AmB intercompartmental clearance (L/h)
AMB_CL    : 0.38  : L-AmB clearance (L/h)
AMB_Kpelf : 0.42  : L-AmB ELF partition
AMB_keqelf: 0.06  : L-AmB lung equilibration (1/h)
AMB_Kpbr  : 0.045 : L-AmB brain partition
AMB_keqbr : 0.02  : L-AmB brain equilibration (1/h)

// ---------------------------- echinocandin ---------------------------------
ECH_Vc    : 30.0  : echinocandin central volume (L)
ECH_Vp    : 20.0  : echinocandin peripheral volume (L)
ECH_Q     : 1.00  : echinocandin intercompartmental clearance (L/h)
ECH_CL    : 0.90  : echinocandin clearance (L/h)
ECH_Kpelf : 0.50  : echinocandin ELF partition
ECH_keqelf: 0.10  : echinocandin ELF equilibration (1/h)
ECH_Kpbr  : 0.02  : echinocandin brain partition

// ---------------------------- co-medication --------------------------------
TAC_Vc    : 1000.0 : tacrolimus apparent volume (L)
TAC_CL    : 30.0   : tacrolimus apparent clearance CL/F (L/h)
TAC_KI_VRC: 1.30   : voriconazole Ki for CYP3A4 (mg/L)
TAC_KI_ISA: 10.0   : isavuconazole Ki for CYP3A4 (mg/L)
TAC_KI_POS: 0.45   : posaconazole Ki for CYP3A4 (mg/L)
GCS_CL    : 1.40   : corticosteroid EFFECT-compartment k_out surrogate (L/h)
GCS_Vc    : 45.0   : corticosteroid effect compartment volume (L)
GCSF_CL   : 0.55   : filgrastim clearance (L/h)
GCSF_Vc   : 4.00   : filgrastim volume (L)

// ------------------------------ fungus -------------------------------------
KGROW     : 0.095  : hyphal growth rate constant (1/h), 0.99 log10/day
BMAX      : 3.0e10 : carrying capacity of one lobe (CFUe)
KGERM     : 0.085  : conidium to germling (1/h)
KMAT      : 0.220  : germling to hypha (1/h)
KMUT      : 2.0e-8 : cyp51A resistance per replication
KGROWBR   : 0.100  : cerebral hyphal growth rate (1/h)
BMAXBR    : 1.0e9  : carrying capacity of cerebral parenchyma (CFUe)
INOCULUM  : 1.0e3  : conidia deposited in the alveoli (CFUe)
RESIST_FRAC: 0.0   : fraction of the inoculum already azole-resistant

// --------------------------- host defence ----------------------------------
KHOST_MAX : 0.320  : maximal recruited-neutrophil hyphal killing (1/h)
KHOST_BASE: 0.220  : constitutive alveolar surveillance at normal ANC (1/h)
KN50      : 0.55   : lesion neutrophils giving half-maximal killing
KMAC_CON  : 0.150  : macrophage conidial killing (1/h per unit MAC)
GCS_IC50  : 0.900  : steroid concentration halving the phagocyte burst (mg/L)
BIO_SHIELD: 0.55   : max fraction of killing blocked by the hyphal matrix
KBIO      : 0.60   : matrix giving half-maximal shielding
KBIO_ON   : 2.2e-11: matrix formation per CFUe per h
KBIO_OFF  : 0.020  : matrix turnover (1/h)
KPROL     : 0.052  : marrow production rate (1/h)
KTR       : 0.052  : marrow to blood transit (1/h)
KOUT_ANC  : 0.056  : neutrophil egress from blood (1/h)
ANC_SS    : 4.00   : reference blood ANC (10^9/L)
ANC0      : 0.05   : ANC at inoculation (10^9/L)
T_ANC_REC : 240.0  : time at which marrow recovery starts (h)
KMARG     : 0.011  : margination into the lesion (1/h)
KMARG_EC50: 1.0e5  : CFUe giving half-maximal neutrophil recruitment
GCSF_EMAX : 3.60   : maximal G-CSF effect on marrow production
GCSF_EC50 : 6.00   : filgrastim concentration for half-maximal effect (ng/mL)
KNLES_IN  : 0.020  : lesion neutrophil influx scaling
KNLES_OUT : 0.030  : lesion neutrophil loss (1/h)
MAC_SS    : 1.00   : reference alveolar macrophage pool
KMAC_TURN : 0.010  : macrophage turnover (1/h)

// ----------------------------- drug PD -------------------------------------
// Azoles inhibit GROWTH (Imax < 1 by construction); polyene KILLS.
// The three azole EC50 coefficients are calibrated so that the three LABELLED
// regimens are equipotent against a wild-type isolate — i.e. the model ENCODES
// the non-inferiority the registration trials measured, and PREDICTS how that
// equipotence comes apart as MIC, perfusion, genotype and site move.
IMAX_AZOLE: 0.955  : maximal azole growth inhibition (structurally < 1)
VRC_EC50M : 3.78   : voriconazole EC50 per unit MIC at the site (mg/L)
ISA_EC50M : 0.525  : isavuconazole EC50 per unit MIC at the site (mg/L)
POS_EC50M : 2.33   : posaconazole EC50 per unit MIC at the site (mg/L)
HILL_AZ   : 2.00   : azole Hill coefficient
AMB_EMAX  : 0.185  : maximal polyene killing rate (1/h)
AMB_EC50M : 3.40   : polyene EC50 per unit MIC at the site (mg/L)
HILL_AMB  : 2.40   : polyene Hill coefficient
ECH_IMAX  : 0.520  : maximal echinocandin growth inhibition
ECH_EC50  : 0.30   : echinocandin EC50 at the site (mg/L)
ECH_CPAR  : 12.0   : concentration above which the paradoxical effect appears
HILL_ECH  : 1.60   : echinocandin Hill coefficient
MIC_VRC   : 0.50   : voriconazole MIC of the susceptible population (mg/L)
MIC_ISA   : 1.00   : isavuconazole MIC susceptible (mg/L)
MIC_POS   : 0.12   : posaconazole MIC susceptible (mg/L)
MIC_AMB   : 1.00   : amphotericin B MIC susceptible (mg/L)
MIC_VRC_R : 8.00   : voriconazole MIC of the resistant population (mg/L)
MIC_ISA_R : 8.00   : isavuconazole MIC resistant (mg/L)
MIC_POS_R : 0.75   : posaconazole MIC resistant (mg/L)
MIC_AMB_R : 1.00   : amphotericin B MIC resistant (unchanged)

// ------------------------ lesion / angioinvasion ---------------------------
KVLES     : 2.0e-10: lesion volume formed per CFUe per h (mL)
KRESOLVE_B: 0.020  : passive lesion resolution (1/h)
KRESOLVE  : 0.010  : extra neutrophil-dependent resolution (1/h)
KANG      : 1.2e-12: angioinvasion index formed per CFUe per h
KANG_OFF  : 0.0090 : angioinvasion repair (1/h)
KPERF_LOSS: 0.020  : maximal perfusion loss rate (1/h), saturating in ANG
KANG50    : 0.35   : angioinvasion index at half-maximal perfusion loss
KPERF_REC : 0.012  : perfusion recovery (1/h)
PERF_GAMMA: 1.00   : exponent in the delivery fraction
DELIV_FLOOR: 0.25  : diffusive delivery to the viable rim, perfusion-independent
KNEC      : 0.0025 : infarction of the perfused rim (1/h)
KNEC_OFF  : 0.0018 : infarct resorption (1/h)
KSEED     : 2.6e-8 : cerebral seeding per unit angioinvasion per CFUe per h

// ---------------------------- cytokines ------------------------------------
KIFNG_ON  : 0.90   : IFN-gamma production (pg/mL/h)
KIFNG_OFF : 0.10   : IFN-gamma elimination (1/h)
KIL6_MAX  : 40.0   : maximal IL-6 production (pg/mL/h)
KIL6_EC50 : 2.0e9  : burden giving half-maximal IL-6 production (CFUe)
KIL6_OFF  : 0.16   : IL-6 elimination (1/h)
KCRP_ON   : 0.014  : CRP production per pg/mL IL-6 per h
KCRP_OFF  : 0.0145 : CRP elimination (1/h), t1/2 ~19 h
KIL10_ON  : 0.25   : IL-10 production (pg/mL/h)
KIL10_OFF : 0.13   : IL-10 elimination (1/h)

// ---------------------------- biomarkers -----------------------------------
// GM is released by ACTIVELY ELONGATING hyphae and by lysis — not by the
// standing stock. This is what makes GM and burden dissociate under a
// fungistatic drug.
KGM_GROWTH: 3.5e-10: GM released per unit growth flux
KGM_LYSIS : 3.9e-10: GM released per unit killing flux
KGM_ABS   : 0.045  : lesion-to-serum GM transfer (1/h), scaled by perfusion
KGM_DEG   : 0.030  : GM degradation in the lesion (1/h)
KGM_EL    : 0.030  : serum GM elimination (1/h), t1/2 ~23 h
GM_VD     : 1.00   : serum GM distribution scaling
KBDG_ON   : 2.0e-8 : BDG release per unit flux
KBDG_OFF  : 0.028  : BDG elimination (1/h)
KPCR_ON   : 1.0e-6 : DNA release per unit flux
KPCR_OFF  : 0.075  : circulating DNA elimination (1/h)

// ------------------------- organ and toxicity ------------------------------
ALT_BASE  : 25.0   : baseline ALT (U/L)
KALT_OFF  : 0.011  : ALT normalisation (1/h)
ALT_VRC_S : 1.55   : ALT rise per mg/L voriconazole above threshold per h
ALT_VRC_T : 4.00   : voriconazole hepatotoxicity threshold (mg/L)
SCR_BASE  : 0.90   : baseline creatinine (mg/dL)
KSCR_OFF  : 0.020  : creatinine normalisation (1/h)
SCR_AMB_S : 3.0e-4 : creatinine rise per mg/L L-AmB above threshold per h
SCR_AMB_T : 18.0   : L-AmB tubular injury threshold (mg/L)
KMG_ON    : 1.4e-4 : K/Mg wasting per mg/L L-AmB per h
KMG_OFF   : 0.020  : K/Mg repletion (1/h)
QTC_VRC   : 3.40   : QTc per mg/L voriconazole (ms)
QTC_POS   : 5.00   : QTc per mg/L posaconazole (ms)
KQTC_OFF  : 0.30   : QTc offset rate (1/h)
NEURO_ON  : 0.062  : neurotoxicity accrual per mg/L above threshold per h
NEURO_OFF : 0.10   : neurotoxicity offset (1/h)
NEURO_THR : 4.00   : voriconazole neurotoxicity threshold (mg/L)
PFR_BASE  : 400.0  : baseline PaO2/FiO2 (mmHg)
KPFR_REC  : 0.010  : oxygenation recovery (1/h)
KPFR_LOSS : 0.030  : oxygenation loss per mL of consolidated lung per h
TEMP_BASE : 36.8   : baseline core temperature (degC)
KTEMP     : 0.10   : temperature relaxation (1/h)
TEMP_GAIN : 0.030  : pyrogenic gain

// ---------------------------- mortality ------------------------------------
HAZ0      : 2.4e-5 : background hazard of the underlying disease (1/h)
HAZ_BURDEN: 2.2e-4 : hazard per log10 CFUe above 4 (1/h)
HAZ_PERF  : 1.0e-3 : hazard per unit (1 - perfusion) (1/h)
HAZ_HYPOX : 1.1e-3 : hazard per unit relative PF-ratio loss (1/h)
HAZ_BRAIN : 1.2e-4 : hazard per log10 cerebral CFUe above 2 (1/h)
HAZ_NEUTRO: 2.5e-4 : hazard while ANC < 0.5 (1/h)
HAZ_ALT   : 1.5e-6 : hazard per U/L ALT above 120 (1/h)
HAZ_SCR   : 2.6e-4 : hazard per mg/dL creatinine above 1.5 (1/h)

$CMT @annotated
VRC_gut  : voriconazole gut depot (mg)
VRC_c    : voriconazole central (mg)
VRC_p    : voriconazole peripheral (mg)
VRC_elf  : voriconazole ELF concentration (mg/L)
VRC_brain: voriconazole brain concentration (mg/L)
ISA_gut  : isavuconazole gut depot (mg)
ISA_c    : isavuconazole central (mg)
ISA_p    : isavuconazole peripheral (mg)
ISA_elf  : isavuconazole ELF concentration (mg/L)
POS_gut  : posaconazole gut depot (mg)
POS_c    : posaconazole central (mg)
POS_p    : posaconazole peripheral (mg)
POS_elf  : posaconazole ELF concentration (mg/L)
AMB_c    : liposomal amphotericin B central (mg)
AMB_p    : liposomal amphotericin B peripheral (mg)
AMB_elf  : liposomal amphotericin B ELF concentration (mg/L)
AMB_brain: liposomal amphotericin B brain concentration (mg/L)
ECH_c    : echinocandin central (mg)
ECH_p    : echinocandin peripheral (mg)
ECH_elf  : echinocandin ELF concentration (mg/L)
TAC_c    : tacrolimus central (mg)
GCS_c    : corticosteroid effect compartment (mg)
GCSF_c   : filgrastim central (mg)
CON      : resting and swollen conidia (CFUe)
GERM     : germlings (CFUe)
HYW      : hyphal biomass azole-susceptible (CFUe)
HYR      : hyphal biomass azole-resistant (CFUe)
BRB      : cerebral fungal burden (CFUe)
BIO      : hyphal matrix / biofilm (relative)
VLES     : lesion volume (mL)
NEC      : infarcted lesion volume (mL)
PERF     : lesion perfusion index (0-1)
ANG      : angioinvasion index (relative)
MAR      : marrow granulocyte reserve (relative)
ANC      : blood absolute neutrophil count (10^9/L)
NLES     : neutrophils recruited into the lesion (relative)
MAC      : functional alveolar macrophages (relative)
IFNG     : interferon gamma (pg/mL)
IL6      : interleukin 6 (pg/mL)
CRP      : C-reactive protein (mg/L)
IL10     : interleukin 10 (pg/mL)
GMles    : galactomannan in the lesion (index x mL)
GMser    : serum galactomannan index
BDG      : serum beta-D-glucan (pg/mL)
PCRs     : serum Aspergillus DNA (copies/mL)
ALT      : alanine aminotransferase (U/L)
SCR      : serum creatinine (mg/dL)
KMG      : renal K/Mg wasting index (relative)
QTC      : QTc prolongation above baseline (ms)
NEUROE   : visual / neuropsychiatric effect (relative)
PFR      : PaO2 / FiO2 (mmHg)
TEMP     : core temperature (degC)
HAZ      : cumulative mortality hazard (dimensionless)

$GLOBAL
#define Cvrc  (VRC_c / VRC_Vc)
#define Cisa  (ISA_c / ISA_Vc)
#define Cpos  (POS_c / POS_Vc)
#define Camb  (AMB_c / AMB_Vc)
#define Cech  (ECH_c / ECH_Vc)
#define Cgcs  (GCS_c / GCS_Vc)

// Hill function guarded against overflow: burdens reach 1e10 and a Hill
// exponent of 4 would overflow a double if evaluated naively.
double hillf(double c, double ec50, double n) {
  if (c <= 0.0) return 0.0;
  if (c > 1.0e6 * ec50) return 1.0;
  double cn = pow(c, n);
  return cn / (pow(ec50, n) + cn);
}

$MAIN
CON_0   = INOCULUM;
PERF_0  = 1.0;
MAC_0   = MAC_SS;
ANC_0   = ANC0;
MAR_0   = ANC0 * KOUT_ANC / KTR;
ALT_0   = ALT_BASE;
SCR_0   = SCR_BASE;
PFR_0   = PFR_BASE;
TEMP_0  = TEMP_BASE;

$ODE
// ---------------------------------------------------------------- perfusion
double PERFx = PERF; if (PERFx < 1.0e-4) PERFx = 1.0e-4; if (PERFx > 1.0) PERFx = 1.0;
// Drug delivery keeps a perfusion-INDEPENDENT floor: antifungal reaches the
// viable rim of an infarcted lesion by diffusion; only the core is cut off.
double fdel  = DELIV_FLOOR + (1.0 - DELIV_FLOOR) * pow(PERFx, PERF_GAMMA);
double Svrc = VRC_elf * fdel;
double Sisa = ISA_elf * fdel;
double Spos = POS_elf * fdel;
double Samb = AMB_elf * fdel;
double Sech = ECH_elf * fdel;

// ------------------------------------------------- voriconazole disposition
double infl = 1.0 - VRC_iEmax * hillf(CRP, VRC_iEC50, 1.0);
double vmax = VRC_Vmax * CYP2C19 * infl;
double clv  = vmax * Cvrc / (VRC_Km + Cvrc) + VRC_CLlin * Cvrc;
double qv   = VRC_Q * (Cvrc - VRC_p / VRC_Vp);
dxdt_VRC_gut   = -VRC_ka * VRC_gut;
dxdt_VRC_c     = VRC_F * VRC_ka * VRC_gut - clv - qv;
dxdt_VRC_p     = qv;
dxdt_VRC_elf   = VRC_keqelf * (VRC_Kpelf * Cvrc - VRC_elf);
dxdt_VRC_brain = VRC_keqbr  * (VRC_Kpbr  * Cvrc - VRC_brain);

double qi = ISA_Q * (Cisa - ISA_p / ISA_Vp);
dxdt_ISA_gut = -ISA_ka * ISA_gut;
dxdt_ISA_c   = ISA_F * ISA_ka * ISA_gut - ISA_CL * Cisa - qi;
dxdt_ISA_p   = qi;
dxdt_ISA_elf = ISA_keqelf * (ISA_Kpelf * Cisa - ISA_elf);

double qp = POS_Q * (Cpos - POS_p / POS_Vp);
dxdt_POS_gut = -POS_ka * POS_gut;
dxdt_POS_c   = POS_F * POS_ka * POS_gut - POS_CL * Cpos - qp;
dxdt_POS_p   = qp;
dxdt_POS_elf = POS_keqelf * (POS_Kpelf * Cpos - POS_elf);

double qa = AMB_Q * (Camb - AMB_p / AMB_Vp);
dxdt_AMB_c     = -AMB_CL * Camb - qa;
dxdt_AMB_p     = qa;
dxdt_AMB_elf   = AMB_keqelf * (AMB_Kpelf * Camb - AMB_elf);
dxdt_AMB_brain = AMB_keqbr  * (AMB_Kpbr  * Camb - AMB_brain);

double qe = ECH_Q * (Cech - ECH_p / ECH_Vp);
dxdt_ECH_c   = -ECH_CL * Cech - qe;
dxdt_ECH_p   = qe;
dxdt_ECH_elf = ECH_keqelf * (ECH_Kpelf * Cech - ECH_elf);

// Competitive CYP3A4 inhibition: each azole contributes C/Ki independently.
double inh = 1.0 / (1.0 + Cvrc / TAC_KI_VRC + Cisa / TAC_KI_ISA + Cpos / TAC_KI_POS);
dxdt_TAC_c  = -TAC_CL * inh * (TAC_c / TAC_Vc);
dxdt_GCS_c  = -GCS_CL * Cgcs;
double Cgcsf = GCSF_c / GCSF_Vc * 1000.0;             // ng/mL
dxdt_GCSF_c = -GCSF_CL * Cgcsf / 1000.0;

// ------------------------------------------------------ host defence terms
double fster   = 1.0 / (1.0 + Cgcs / GCS_IC50);
double shield  = 1.0 - BIO_SHIELD * hillf(BIO, KBIO, 1.0);
double Neff    = hillf(NLES, KN50, 1.0) * fster * shield;
double surv    = (ANC / ANC_SS > 1.5 ? 1.5 : ANC / ANC_SS) * fster * shield;
double killh   = KHOST_MAX * Neff + KHOST_BASE * surv;

// ------------------------------------------------------------- drug effect
double Iw = IMAX_AZOLE * hillf(Svrc, VRC_EC50M * MIC_VRC, HILL_AZ)
          + IMAX_AZOLE * hillf(Sisa, ISA_EC50M * MIC_ISA, HILL_AZ)
          + IMAX_AZOLE * hillf(Spos, POS_EC50M * MIC_POS, HILL_AZ);
double Ir = IMAX_AZOLE * hillf(Svrc, VRC_EC50M * MIC_VRC_R, HILL_AZ)
          + IMAX_AZOLE * hillf(Sisa, ISA_EC50M * MIC_ISA_R, HILL_AZ)
          + IMAX_AZOLE * hillf(Spos, POS_EC50M * MIC_POS_R, HILL_AZ);
if (Iw > IMAX_AZOLE) Iw = IMAX_AZOLE;
if (Ir > IMAX_AZOLE) Ir = IMAX_AZOLE;

// echinocandin: apical growth inhibition, attenuated above ECH_CPAR (Eagle effect)
double par  = 1.0 / (1.0 + pow(Sech / ECH_CPAR, 3.0));
double Iech = ECH_IMAX * hillf(Sech, ECH_EC50, HILL_ECH) * par;
Iw = 1.0 - (1.0 - Iw) * (1.0 - Iech); if (Iw > 0.995) Iw = 0.995;
Ir = 1.0 - (1.0 - Ir) * (1.0 - Iech); if (Ir > 0.995) Ir = 0.995;

double killaw = AMB_EMAX * hillf(Samb, AMB_EC50M * MIC_AMB,   HILL_AMB);
double killar = AMB_EMAX * hillf(Samb, AMB_EC50M * MIC_AMB_R, HILL_AMB);

// ---------------------------------------------------------------- the fungus
double HY  = HYW + HYR;
double lg  = 1.0 - HY / BMAX; if (lg < 0.0) lg = 0.0;
double kgerm = KGERM * (1.0 - 0.55 * hillf(MAC * fster, 0.7, 1.0));
double ckill = KMAC_CON * MAC * fster;
dxdt_CON  = -kgerm * CON - ckill * CON;
dxdt_GERM = kgerm * CON - KMAT * GERM - 0.6 * killh * GERM - ckill * 0.35 * GERM;

// Extinction gate: a continuous state can carry 1e-30 "CFUe" and regrow from
// it. Nothing below one organism can divide.
double gw = hillf(HYW, 1.0, 4.0);
double gr = hillf(HYR, 1.0, 4.0);
double grow_w = KGROW * (1.0 - Iw) * lg * gw;
double grow_r = KGROW * (1.0 - Ir) * lg * gr;
double mut    = KMUT * grow_w * HYW;
double lyse_w = (killh + killaw) * HYW;
double lyse_r = (killh + killar) * HYR;
dxdt_HYW = KMAT * GERM * (1.0 - RESIST_FRAC) + grow_w * HYW - lyse_w - mut;
dxdt_HYR = KMAT * GERM * RESIST_FRAC         + grow_r * HYR - lyse_r + mut;

double seed = KSEED * ANG * HY;
double killb = killh * 0.25 + AMB_EMAX * hillf(AMB_brain, AMB_EC50M * MIC_AMB, HILL_AMB);
double Ib = IMAX_AZOLE * hillf(VRC_brain, VRC_EC50M * MIC_VRC, HILL_AZ)
          + IMAX_AZOLE * hillf(Cisa * ISA_Kpbr, ISA_EC50M * MIC_ISA, HILL_AZ);
if (Ib > IMAX_AZOLE) Ib = IMAX_AZOLE;
double lgbr = 1.0 - BRB / BMAXBR; if (lgbr < 0.0) lgbr = 0.0;
dxdt_BRB = seed * lgbr + KGROWBR * (1.0 - Ib) * lgbr * hillf(BRB, 1.0, 4.0) * BRB
           - killb * BRB;
dxdt_BIO = KBIO_ON * HY - KBIO_OFF * BIO;

// -------------------------------------------------- lesion and angioinvasion
dxdt_VLES = KVLES * HY - (KRESOLVE_B + KRESOLVE * hillf(NLES, 0.5, 1.0)) * VLES;
dxdt_ANG  = KANG * HY - KANG_OFF * ANG;
double heal = hillf(NLES, 0.5, 1.0) * (1.0 - hillf(HY, 1.0e7, 1.0));
dxdt_PERF = -KPERF_LOSS * hillf(ANG, KANG50, 1.0) * PERFx
            + KPERF_REC * (1.0 - PERFx) * (0.25 + heal);
double rim = VLES - NEC; if (rim < 0.0) rim = 0.0;
dxdt_NEC  = KNEC * (1.0 - PERFx) * rim - KNEC_OFF * NEC;

// ----------------------------------------------------------- haematopoiesis
double marrow_on = (SOLVERTIME >= T_ANC_REC) ? 1.0 : 0.0;
double egcsf = GCSF_EMAX * hillf(Cgcsf, GCSF_EC50, 1.0);
dxdt_MAR  = KPROL * marrow_on * (1.0 + egcsf) * ANC_SS - KTR * MAR;
double marg = KMARG * hillf(HY, KMARG_EC50, 1.0) * ANC;
dxdt_ANC  = KTR * MAR - KOUT_ANC * ANC - marg;
dxdt_NLES = KNLES_IN * marg * 100.0 * fster - KNLES_OUT * NLES;
dxdt_MAC  = KMAC_TURN * (MAC_SS - MAC);

// ---------------------------------------------------------------- cytokines
dxdt_IFNG = KIFNG_ON * hillf(HY, 1.0e6, 1.0) * fster - KIFNG_OFF * IFNG;
dxdt_IL6  = KIL6_MAX * hillf(HY, KIL6_EC50, 1.0) * (0.3 + 0.7 * fster) - KIL6_OFF * IL6;
dxdt_CRP  = KCRP_ON * IL6 - KCRP_OFF * CRP;
dxdt_IL10 = KIL10_ON * hillf(IL6, 250.0, 1.0) - KIL10_OFF * IL10;

// --------------------------------------------------------------- biomarkers
double gflux = grow_w * HYW + grow_r * HYR;
double lflux = lyse_w + lyse_r;
// GM entering the blood is a TRANSVASCULAR flux, so it takes the bare
// perfusion term and not the diffusive delivery floor.
dxdt_GMles = KGM_GROWTH * gflux + KGM_LYSIS * lflux - KGM_ABS * PERFx * GMles - KGM_DEG * GMles;
dxdt_GMser = KGM_ABS * PERFx * GMles / GM_VD - KGM_EL * GMser;
dxdt_BDG   = KBDG_ON * (gflux + lflux) - KBDG_OFF * BDG;
dxdt_PCRs  = KPCR_ON * (0.25 * gflux + lflux) - KPCR_OFF * PCRs;

// ----------------------------------------------------------------- toxicity
double vex = Cvrc - ALT_VRC_T; if (vex < 0.0) vex = 0.0;
double iex = Cisa - 8.0;       if (iex < 0.0) iex = 0.0;
dxdt_ALT = ALT_VRC_S * vex + 0.9 * iex - KALT_OFF * (ALT - ALT_BASE);
double aex = Camb - SCR_AMB_T; if (aex < 0.0) aex = 0.0;
dxdt_SCR = SCR_AMB_S * aex - KSCR_OFF * (SCR - SCR_BASE);
dxdt_KMG = KMG_ON * Camb - KMG_OFF * KMG;
dxdt_QTC = QTC_VRC * Cvrc + QTC_POS * Cpos - KQTC_OFF * QTC;
double nex = Cvrc - NEURO_THR; if (nex < 0.0) nex = 0.0;
dxdt_NEUROE = NEURO_ON * nex - NEURO_OFF * NEUROE;
dxdt_PFR = KPFR_REC * (PFR_BASE - PFR) - KPFR_LOSS * (VLES + 0.8 * NEC) * PFR / PFR_BASE;
dxdt_TEMP = KTEMP * (TEMP_BASE - TEMP) + TEMP_GAIN * hillf(IL6, 90.0, 1.0) * fster * 10.0;

// ----------------------------------------------------------------- mortality
double lb  = log10(HY  > 1.0 ? HY  : 1.0);
double lbb = log10(BRB > 1.0 ? BRB : 1.0);
double h = HAZ0
         + HAZ_BURDEN * ((lb  - 4.0) > 0.0 ? (lb  - 4.0) : 0.0)
         + HAZ_PERF   * (1.0 - PERFx)
         + HAZ_HYPOX  * ((1.0 - PFR / PFR_BASE) > 0.0 ? (1.0 - PFR / PFR_BASE) : 0.0)
         + HAZ_BRAIN  * ((lbb - 2.0) > 0.0 ? (lbb - 2.0) : 0.0)
         + (ANC < 0.5 ? HAZ_NEUTRO : 0.0)
         + HAZ_ALT    * ((ALT - 120.0) > 0.0 ? (ALT - 120.0) : 0.0)
         + HAZ_SCR    * ((SCR - 1.5)   > 0.0 ? (SCR - 1.5)   : 0.0);
dxdt_HAZ = h;

$TABLE
double CP_VRC  = Cvrc;
double CP_ISA  = Cisa;
double CP_POS  = Cpos;
double CP_AMB  = Camb;
double CP_ECH  = Cech;
double CP_TAC  = TAC_c / TAC_Vc * 1000.0;                 // ng/mL
double BURDEN  = HYW + HYR;
double LOGB    = log10((BURDEN > 1.0 ? BURDEN : 1.0));
double LOGBR   = log10((BRB    > 1.0 ? BRB    : 1.0));
double RFRAC   = (BURDEN > 0.0) ? HYR / BURDEN : 0.0;
double SURV    = exp(-HAZ);
double MORT    = 100.0 * (1.0 - SURV);
double GMPOS   = (GMser >= 0.5) ? 1.0 : 0.0;
// Reporting-only, exactly as computed in the Python reference outside derivs:
double CU_VRC  = Cvrc * VRC_fu;               // unbound voriconazole, drives fAUC/MIC
double CBR_ISA = Cisa * ISA_Kpbr;             // algebraic brain concentrations for
double CBR_POS = Cpos * POS_Kpbr;             // the drugs given no brain compartment
double CBR_ECH = Cech * ECH_Kpbr;

$CAPTURE CP_VRC CP_ISA CP_POS CP_AMB CP_ECH CP_TAC CU_VRC CBR_ISA CBR_POS CBR_ECH
         BURDEN LOGB LOGBR RFRAC SURV MORT GMPOS
'

mod <- mcode_cache("ipa_qsp", code)

## =============================================================================
##  DOSING HELPERS
## =============================================================================
WT <- 70

cmt_of <- function(name) which(mrgsolve::cmt(mod) == name)

.ev <- function(time, cmt, amt, rate = 0, ii = 0, addl = 0)
  data.frame(ID = 1, time = time, cmt = cmt_of(cmt), amt = amt,
             rate = rate, ii = ii, addl = addl, evid = 1)

## Voriconazole: 6 mg/kg IV q12h x2 loading, then 4 mg/kg IV q12h, 2 h infusions
rx_vrc <- function(start_d = 3, days = 42, load = 6, maint = 4) {
  t0 <- start_d * 24
  rbind(
    .ev(t0,      "VRC_c", load  * WT, rate = load  * WT / 2, ii = 12, addl = 1),
    .ev(t0 + 24, "VRC_c", maint * WT, rate = maint * WT / 2, ii = 12, addl = days * 2 - 1))
}
## Isavuconazole: 200 mg q8h x6 then 200 mg q24h (1 h infusions)
rx_isa <- function(start_d = 3, days = 42) {
  t0 <- start_d * 24
  rbind(.ev(t0,      "ISA_c", 200, rate = 200, ii = 8,  addl = 5),
        .ev(t0 + 48, "ISA_c", 200, rate = 200, ii = 24, addl = days - 1))
}
## Posaconazole DR tablet: 300 mg bid day 1 then 300 mg qd
rx_pos <- function(start_d = 0, days = 90) {
  t0 <- start_d * 24
  rbind(.ev(t0,      "POS_gut", 300, ii = 12, addl = 1),
        .ev(t0 + 24, "POS_gut", 300, ii = 24, addl = days - 1))
}
## Liposomal amphotericin B: 2 h infusion q24h
rx_amb <- function(start_d = 3, days = 42, mgkg = 3)
  .ev(start_d * 24, "AMB_c", mgkg * WT, rate = mgkg * WT / 2, ii = 24, addl = days - 1)
## Anidulafungin 200 mg load then 100 mg q24h (1.5 h infusions)
rx_ech <- function(start_d = 3, days = 42) {
  t0 <- start_d * 24
  rbind(.ev(t0,      "ECH_c", 200, rate = 200 / 1.5),
        .ev(t0 + 24, "ECH_c", 100, rate = 100 / 1.5, ii = 24, addl = days - 1))
}
rx_steroid <- function(mg_day = 60, days = 60, start_d = 0)
  .ev(start_d * 24, "GCS_c", mg_day, rate = mg_day, ii = 24, addl = days - 1)
rx_gcsf <- function(start_d = 3, days = 10, mcg = 300)
  .ev(start_d * 24, "GCSF_c", mcg / 1000, rate = mcg / 1000, ii = 24, addl = days - 1)
rx_tac <- function(days = 90, mg_bid = 3)
  .ev(0, "TAC_c", mg_bid, ii = 12, addl = days * 2 - 1)
rx_none <- function() .ev(0, "VRC_c", 0)

## =============================================================================
##  CYP2C19 DIPLOTYPES  (Vmax multipliers; Moriyama 2017 CPIC; Weiss 2009)
## =============================================================================
CYP2C19 <- c("UM (*17/*17)" = 1.85, "RM (*1/*17)" = 1.38, "NM (*1/*1)" = 1.00,
             "IM (*1/*2)"   = 0.62, "PM (*2/*2)"  = 0.30)

sim <- function(dosing = rx_none(), param = list(), end = 84 * 24, delta = 2) {
  m <- mod
  if (length(param)) m <- mrgsolve::param(m, param)
  mrgsim_d(m, data = dosing, end = end, delta = delta, atol = 1e-10, rtol = 1e-8,
           hmax = 0.5) %>% as.data.frame()
}

## =============================================================================
##  23 TREATMENT SCENARIOS
## =============================================================================
##  The scenario set is built to answer five questions the trials cannot:
##    (1) how much of the outcome is the drug and how much is the marrow
##    (2) where the time-to-treatment cliff is
##    (3) what CYP2C19 genotype costs and what TDM buys back
##    (4) what an azole-resistant isolate does to each arm
##    (5) whether galactomannan tracks the burden it is used as a surrogate for
## =============================================================================
SCENARIOS <- list(
  S01 = list(lab = "Untreated, persistent neutropenia",
             dose = rx_none(), par = list(T_ANC_REC = 1e9)),
  S02 = list(lab = "Untreated, ANC recovery day 10",
             dose = rx_none(), par = list(T_ANC_REC = 240)),
  S03 = list(lab = "Voriconazole, persistent neutropenia",
             dose = rx_vrc(), par = list(T_ANC_REC = 1e9)),
  S04 = list(lab = "Voriconazole NM, ANC recovery day 10",
             dose = rx_vrc(), par = list(T_ANC_REC = 240)),
  S05 = list(lab = "Voriconazole UM (*17/*17)",
             dose = rx_vrc(), par = list(T_ANC_REC = 240, CYP2C19 = 1.85)),
  S06 = list(lab = "Voriconazole PM (*2/*2)",
             dose = rx_vrc(), par = list(T_ANC_REC = 240, CYP2C19 = 0.30)),
  S07 = list(lab = "Voriconazole UM, TDM-escalated 7.6 mg/kg",
             dose = rx_vrc(maint = 7.6), par = list(T_ANC_REC = 240, CYP2C19 = 1.85)),
  S08 = list(lab = "Voriconazole PM, TDM-reduced 1.8 mg/kg",
             dose = rx_vrc(maint = 1.8), par = list(T_ANC_REC = 240, CYP2C19 = 0.30)),
  S09 = list(lab = "Isavuconazole 200 mg q24h",
             dose = rx_isa(), par = list(T_ANC_REC = 240)),
  S10 = list(lab = "Liposomal amphotericin B 3 mg/kg",
             dose = rx_amb(mgkg = 3), par = list(T_ANC_REC = 240)),
  S11 = list(lab = "Liposomal amphotericin B 10 mg/kg",
             dose = rx_amb(mgkg = 10), par = list(T_ANC_REC = 240)),
  S12 = list(lab = "Anidulafungin monotherapy",
             dose = rx_ech(), par = list(T_ANC_REC = 240)),
  S13 = list(lab = "Voriconazole + anidulafungin",
             dose = rbind(rx_vrc(), rx_ech()), par = list(T_ANC_REC = 240)),
  S14 = list(lab = "TR34/L98H (VRC MIC 8), voriconazole",
             dose = rx_vrc(), par = list(T_ANC_REC = 240, RESIST_FRAC = 1)),
  S15 = list(lab = "TR34/L98H, switched to L-AmB on day 10",
             dose = rbind(rx_vrc(days = 7), rx_amb(start_d = 10, days = 35)),
             par = list(T_ANC_REC = 240, RESIST_FRAC = 1)),
  S16 = list(lab = "Steroid host (no neutropenia), voriconazole",
             dose = rbind(rx_vrc(), rx_steroid(60, 60)),
             par = list(T_ANC_REC = 0, ANC0 = 4)),
  S17 = list(lab = "Steroid host, steroid tapered on day 14",
             dose = rbind(rx_vrc(), rx_steroid(60, 14)),
             par = list(T_ANC_REC = 0, ANC0 = 4)),
  S18 = list(lab = "Voriconazole + G-CSF from day 3",
             dose = rbind(rx_vrc(), rx_gcsf(3, 12)), par = list(T_ANC_REC = 240)),
  S19 = list(lab = "Late start, voriconazole from day 10",
             dose = rx_vrc(start_d = 10), par = list(T_ANC_REC = 240)),
  S20 = list(lab = "CNS disease, voriconazole",
             dose = rx_vrc(), par = list(T_ANC_REC = 240, KSEED = 2.6e-7)),
  S21 = list(lab = "CNS disease, L-AmB 5 mg/kg",
             dose = rx_amb(mgkg = 5), par = list(T_ANC_REC = 240, KSEED = 2.6e-7)),
  S22 = list(lab = "Transplant: tacrolimus + voriconazole",
             dose = rbind(rx_vrc(), rx_tac(), rx_steroid(20, 84)),
             par = list(T_ANC_REC = 0, ANC0 = 3)),
  S23 = list(lab = "Posaconazole prophylaxis, breakthrough resistant isolate",
             dose = rx_pos(0, 90), par = list(T_ANC_REC = 240, RESIST_FRAC = 1))
)

run_scenarios <- function() {
  bind_rows(lapply(names(SCENARIOS), function(k) {
    s <- SCENARIOS[[k]]
    sim(s$dose, s$par) %>% mutate(scenario = k, label = s$lab)
  }))
}

summarise_scenarios <- function(out) {
  out %>%
    group_by(scenario, label) %>%
    summarise(
      logB_peak    = max(LOGB),
      logB_d14     = LOGB[which.min(abs(time - 14 * 24))],
      logB_d42     = LOGB[which.min(abs(time - 42 * 24))],
      logB_d84     = LOGB[which.min(abs(time - 84 * 24))],
      GM_peak      = max(GMser),
      GM_d14       = GMser[which.min(abs(time - 14 * 24))],
      PERF_min     = min(PERF),
      lesion_mL    = max(VLES),
      ALT_max      = max(ALT),
      SCr_max      = max(SCR),
      mort_d42_pct = MORT[which.min(abs(time - 42 * 24))],
      mort_d84_pct = MORT[which.min(abs(time - 84 * 24))],
      .groups = "drop")
}

## =============================================================================
##  THE FIVE ANALYSES THAT MOTIVATE THE MODEL
## =============================================================================

## (1) The structural claim, in closed form. No simulation needed.
azole_growth_floor <- function() {
  p <- as.list(mrgsolve::param(mod))
  fl <- p$KGROW * (1 - p$IMAX_AZOLE)
  data.frame(
    growth_floor_per_h        = fl,
    growth_floor_log10_per_d  = fl * 24 / log(10),
    doubling_time_h           = log(2) / fl,
    min_Neff_for_clearance    = fl / p$KHOST_MAX)
}

## (2) Time-to-treatment bifurcation.
delay_sweep <- function(days = c(1:8, 10, 12, 14)) {
  bind_rows(lapply(days, function(d) {
    o <- sim(rx_vrc(start_d = d), list(T_ANC_REC = 336))
    data.frame(start_day = d,
               logB_d28  = o$LOGB[which.min(abs(o$time - 28 * 24))],
               logB_d84  = tail(o$LOGB, 1),
               PERF_min  = min(o$PERF),
               mort_d84  = tail(o$MORT, 1))
  }))
}

## (3) Drug lever vs host lever. The comparison the trials are not powered for.
lever_grid <- function() {
  arms <- list(none = rx_none(), VRC = rx_vrc(), ISA = rx_isa(),
               AMB = rx_amb(mgkg = 3), VRC_ECH = rbind(rx_vrc(), rx_ech()))
  recs <- c(5, 10, 14, 21, 28, 1e9)
  bind_rows(lapply(recs, function(r) {
    row <- lapply(arms, function(a) tail(sim(a, list(T_ANC_REC = ifelse(r > 1e8, 1e9, r * 24)))$MORT, 1))
    data.frame(anc_recovery_day = ifelse(r > 1e8, NA, r), as.data.frame(row))
  }))
}

## (4) The galactomannan dissociation.
gm_dissociation <- function() {
  arms <- list(`no therapy` = rx_none(), voriconazole = rx_vrc(),
               `L-AmB 3 mg/kg` = rx_amb(mgkg = 3), anidulafungin = rx_ech())
  bind_rows(lapply(names(arms), function(a) {
    o <- sim(arms[[a]], list(T_ANC_REC = 336), end = 20 * 24, delta = 0.5)
    pick <- function(d, v) o[[v]][which.min(abs(o$time - d * 24))]
    data.frame(arm = a, GM_d3 = pick(3, "GMser"), GM_d5 = pick(5, "GMser"),
               GM_d7 = pick(7, "GMser"), GM_d10 = pick(10, "GMser"),
               GM_d14 = pick(14, "GMser"),
               logB_d7 = pick(7, "LOGB"), logB_d14 = pick(14, "LOGB"))
  }))
}

## (5) MIC ladder and PK/PD target attainment.
mic_ladder <- function(mics = c(0.125, 0.25, 0.5, 1, 2, 4, 8)) {
  bind_rows(lapply(mics, function(m) {
    o <- sim(rx_vrc(), list(T_ANC_REC = 336, MIC_VRC = m))
    data.frame(MIC = m,
               logB_d42 = o$LOGB[which.min(abs(o$time - 42 * 24))],
               mort_d84 = tail(o$MORT, 1))
  }))
}

## Delivery ablation. DELIV_FLOOR = 1 gives the drug full ELF access however
## thrombosed the lesion is. In the reference implementation this does NOT
## move the time-to-treatment cliff -- see README, section "a hypothesis of
## my own that did not survive".
perfusion_ablation <- function(days = c(3, 5, 7, 10, 14)) {
  bind_rows(lapply(days, function(d) {
    a <- sim(rx_vrc(start_d = d), list(T_ANC_REC = 336))
    b <- sim(rx_vrc(start_d = d), list(T_ANC_REC = 336, DELIV_FLOOR = 1))
    data.frame(start_day = d,
               peak_logB_feedback = max(a$LOGB), peak_logB_deleted = max(b$LOGB),
               mort_feedback = tail(a$MORT, 1), mort_deleted = tail(b$MORT, 1))
  }))
}

## PK qualification against the product labels.
pk_qualification <- function() {
  g <- bind_rows(lapply(names(CYP2C19), function(k) {
    o <- sim(rx_vrc(start_d = 0, days = 14), list(CYP2C19 = CYP2C19[[k]], INOCULUM = 0),
             end = 14 * 24, delta = 0.25)
    w <- o[o$time >= 12 * 24 & o$time <= 13 * 24, ]
    data.frame(regimen = paste("VRC 4 mg/kg q12h,", k),
               Cmin = min(w$CP_VRC), Cmax = max(w$CP_VRC),
               AUC24 = sum(diff(w$time) * (head(w$CP_VRC, -1) + tail(w$CP_VRC, -1)) / 2))
  }))
  others <- list(
    `ISA 200 mg q24h`      = list(rx_isa(0, 14), "CP_ISA"),
    `POS DR 300 mg q24h`   = list(rx_pos(0, 14), "CP_POS"),
    `L-AmB 3 mg/kg q24h`   = list(rx_amb(0, 14, 3), "CP_AMB"),
    `L-AmB 10 mg/kg q24h`  = list(rx_amb(0, 14, 10), "CP_AMB"),
    `Anidulafungin 100 mg` = list(rx_ech(0, 14), "CP_ECH"))
  o2 <- bind_rows(lapply(names(others), function(k) {
    o <- sim(others[[k]][[1]], list(INOCULUM = 0), end = 14 * 24, delta = 0.25)
    w <- o[o$time >= 12 * 24 & o$time <= 13 * 24, ]
    v <- w[[others[[k]][[2]]]]
    data.frame(regimen = k, Cmin = min(v), Cmax = max(v),
               AUC24 = sum(diff(w$time) * (head(v, -1) + tail(v, -1)) / 2))
  }))
  bind_rows(g, o2)
}

## Azole-tacrolimus interaction.
ddi_table <- function() {
  arms <- list(`tacrolimus alone` = rx_tac(),
               `+ voriconazole`   = rbind(rx_tac(), rx_vrc(start_d = 0)),
               `+ isavuconazole`  = rbind(rx_tac(), rx_isa(start_d = 0)),
               `1/3 dose + voriconazole` = rbind(rx_tac(mg_bid = 1), rx_vrc(start_d = 0)))
  base <- NA
  bind_rows(lapply(names(arms), function(a) {
    o <- sim(arms[[a]], list(T_ANC_REC = 0, ANC0 = 3), end = 20 * 24, delta = 0.5)
    tr <- o$CP_TAC[which.min(abs(o$time - 14 * 24))]
    if (is.na(base)) base <<- tr
    data.frame(arm = a, trough_ng_mL = tr, fold = tr / base)
  }))
}

## =============================================================================
##  PARAMETER PROVENANCE AND CALIBRATION NOTES
## =============================================================================
##  PK, qualified against product labels and published population analyses:
##    voriconazole  4 mg/kg IV q12h, CYP2C19 NM  -> Cmin 1.13, AUC24 47.9 mg*h/L
##                  the same mg/kg dose gives AUC24 25.0 in a *17/*17 UM and
##                  133.8 in a *2/*2 PM: a 5.4-fold genotype-driven spread that
##                  no weight-based dosing rule can remove (Purkins 2002;
##                  Pascual 2008; Moriyama 2017 CPIC).
##    isavuconazole 200 mg q24h -> Cmin 3.14, AUC24 92.5 mg*h/L (Desai 2016)
##    posaconazole  DR tab 300 mg q24h -> Cmin 1.13, AUC24 34 mg*h/L (Krishna)
##    L-AmB 3 mg/kg -> Cmax 73.3, AUC24 552 mg*h/L (AmBisome label; Walsh 2001)
##    anidulafungin 100 mg q24h -> Cmax 6.5, AUC24 111 mg*h/L
##
##  PD:
##    KGROW 0.095/h gives 0.99 log10/day in the untreated neutropenic host, so
##    a 1e3 alveolar inoculum reaches carrying capacity in ~7.5 days. An
##    earlier value of 0.140/h got there in 5 days and left no clinically
##    meaningful treatment window at all.
##    IMAX_AZOLE 0.955 is the load-bearing parameter of the whole model. It is
##    < 1 because ergosterol depletion slows hyphal elongation but does not
##    lyse hyphae; every in-vitro time-kill curve for a triazole against
##    A. fumigatus is a growth-rate curve, not a killing curve.
##    AMB_EMAX 0.185/h is calibrated to AmBiLoad 12-week survival.
##
##  Mortality:
##    HAZ_* are calibrated on the VIRTUAL POPULATION mean, not on any single
##    scenario: a 320-patient cohort drawing CYP2C19 diplotype, marrow recovery
##    day, therapy start day, MIC, steroid exposure and inoculum from plausible
##    distributions gives 40.8% predicted 12-week mortality, inside the 30-60%
##    range of published series. Individual favourable scenarios sit well below
##    that (voriconazole, treatment on day 3, marrow back on day 10: 15.5%) and
##    unfavourable ones above it (no marrow recovery: >98% on any azole).
##    Untreated persistent neutropenia is uniformly fatal by construction.
##    HAZ0 (2.4e-5/h, ~4.7% over 12 weeks) is the competing risk of the
##    underlying haematological disease, not of the infection.
##
##  WHAT WAS AND WAS NOT EXECUTED — read this before trusting any number
##    The numbers quoted in README.md were produced by `ipa_reference_model.py`,
##    which was run to completion in the authoring environment. THIS FILE was
##    NOT executed there: no R installation was available. It is a faithful
##    transcription of the same equation sheet, checked line by line against
##    the Python implementation, but it has not been compiled by mrgsolve or
##    run. Treat the first `Rscript ipa_mrgsolve_model.R` as a smoke test and
##    compare its output against `ipa_reference_output.txt`; the two should
##    agree to integration tolerance, and any disagreement is a transcription
##    error in this file.
##
##  KNOWN LIMITATIONS — stated because they bound what the model may be used for
##    * The azole EC50 coefficients ENCODE the non-inferiority of the three
##      licensed triazoles rather than predicting it. Comparisons BETWEEN
##      licensed azoles at label dose in a wild-type isolate are therefore not
##      model predictions; comparisons ACROSS MIC, genotype, perfusion and site
##      are.
##    * Echinocandin monotherapy in profound neutropenia is an extrapolation
##      beyond the trial evidence, which is salvage-only and uncontrolled.
##    * The model is deterministic. An extinction gate stops the fungal states
##      regrowing from fractional organisms, but true stochastic clearance
##      near one CFU is not represented.
##    * Immune reconstitution inflammatory syndrome on neutrophil recovery is
##      not modelled; the model treats returning neutrophils as purely
##      beneficial, which is not always what the chest CT shows.
##    * The brain compartment carries TOTAL drug concentration. Liposomal
##      amphotericin B is therefore credited with a higher cerebral
##      concentration than voriconazole, and the model does NOT reproduce the
##      clinical preference for voriconazole in cerebral aspergillosis. Do not
##      use it to choose a drug for CNS disease.
##    * The time-to-treatment threshold is real in the model but its LOCATION
##      (day 6 under the default host) is parameter-dependent. The transferable
##      finding is that a threshold exists, not the date.
## =============================================================================

if (identical(environment(), globalenv()) && !interactive()) {
  cat("\n== PK qualification ==\n");            print(pk_qualification(), digits = 4)
  cat("\n== Azole growth floor (closed form) ==\n"); print(azole_growth_floor(), digits = 4)
  out <- run_scenarios()
  cat("\n== 23 scenarios ==\n");                print(as.data.frame(summarise_scenarios(out)), digits = 4)
  cat("\n== Time-to-treatment bifurcation ==\n"); print(delay_sweep(), digits = 4)
  cat("\n== Drug lever vs host lever ==\n");    print(lever_grid(), digits = 4)
  cat("\n== Galactomannan dissociation ==\n");  print(gm_dissociation(), digits = 4)
  cat("\n== MIC ladder ==\n");                  print(mic_ladder(), digits = 4)
  cat("\n== Perfusion ablation ==\n");          print(perfusion_ablation(), digits = 4)
  cat("\n== Azole-tacrolimus DDI ==\n");        print(ddi_table(), digits = 4)
}
