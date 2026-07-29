##############################################################################
## flu_mrgsolve_model.R
## Influenza A infection — Quantitative Systems Pharmacology model
## ===========================================================================
##
## THE STRUCTURAL CLAIM
## --------------------
## "Start antivirals within 48 hours" is taught as a rule about the clock.
## This model poses it instead as a rule about a QUANTITY THAT IS BEING SPENT.
##
## For a drug started at t_rx, the reduction it can produce in viral AUC is
## bounded above, exactly, by the RESIDUAL AUC of the untreated trajectory
##
##       R(t_rx) = INT_{t_rx}^{inf} ( log10 V(t) - LOD )^+ dt
##
## An antiviral cannot subtract viral load that has already happened.  R is a
## property of the UNTREATED course alone: it contains no drug, no potency and
## no mechanism, and no molecule can cross it.  Everything else in this model
## is an account of how fast R falls and how much of it a given operator takes.
##
## In the calibrated adult (analysis A1 of the reference implementation):
##       at symptom onset          R = 98% of the total viral AUC
##       24 h later (median enrol) R = 75%
##       48 h later (licence edge) R = 46%
##
## and at the median enrolment time the calibrated baloxavir dose already
## takes 66% of what is left, so that
##
##       potency headroom (Emax -> 1)     +0.09 log10.d
##       timing  headroom (dose at onset) +4.34 log10.d
##
## Fifty-fold more of the available benefit lies on the axis a better molecule
## cannot buy.  That is the quantitative content of "treat early".
##
## OPERATOR CLASSIFICATION
## -----------------------
## Every therapy is classified by WHICH TERM of the replication loop it
## touches, because the classes have different signatures in time and are not
## interchangeable (analysis A5):
##
##   ENTRY (beta)              mucosal IgA, neutralising mAb, NAI (weakly)
##   TRANSCRIPTION (E->I)      baloxavir — PA endonuclease, cap-snatching
##   PRODUCTION / RELEASE (p)  oseltamivir, zanamivir, peramivir, baloxavir,
##                             favipiravir
##   MUTAGENESIS               favipiravir — progeny made non-viable
##   VIRION CLEARANCE (c)      neutralising mAb, convalescent plasma
##   TARGET PROTECTION (T->R)  type I/III interferon, ISG antiviral state
##   INFECTED-CELL DEATH       CD8 CTL, NK  (the strongest operator in the
##                             model, and no licensed antiviral uses it)
##   IMMUNOPATHOLOGY DAMPING   corticosteroids — negative here
##
## RESISTANCE IS COMPETITIVE RELEASE
## ---------------------------------
## The resistant subpopulation is never seeded.  It arises from the mutation
## term acting on the wild type from the first replication cycle, and it grows
## when the drug REMOVES ITS COMPETITOR and hands it the target cells the wild
## type would have taken.  Two consequences fall out (A6) that were not put
## in: selection is non-monotone in potency, and it requires a target-cell
## field to still exist — so it needs the same early dosing that clinical
## benefit needs.
##
## 50 ODEs: 18 virus/epithelium (URT + LRT, wild type + resistant), 9 immune,
## 3 clinical, 17 drug PK, 3 bookkeeping.
##
## VERIFICATION
## ------------
## The build environment for this repository has no R toolchain, so this
## system has ALSO been transcribed into `flu_reference_check.py` (numpy /
## scipy LSODA) and integrated there.  Every number quoted in README.md and in
## the comments below is produced by a function in that file and is reproduced
## verbatim in `flu_reference_output.txt`.  The two transcriptions are meant to
## be the same system with the same parameter block; if they disagree, one of
## them is wrong.
##
##   python3 flu_reference_check.py            # all analyses A0-A13
##   python3 flu_reference_check.py --only A1  # the bound
##   python3 flu_reference_check.py --list
##
## Usage:
##   source("flu_mrgsolve_model.R")
##   FLU_baseline()                  # A0  natural history / calibration
##   FLU_residual_bound()            # A1  the bound, and how much is taken
##   FLU_pd_calibration()            # A2  in-vivo vs in-vitro, Hill slope
##   FLU_run_scenarios()             # A3  the 14 shipped scenarios
##   FLU_timing_vs_potency()         # A4  the two headrooms
##   FLU_operator_decomposition()    # A5  which term does each drug touch
##   FLU_competitive_release()       # A6  resistance as competitive release
##   FLU_host_comparison()           # A7  who selects I38T (and a failure)
##   FLU_symptom_requirement()       # A8  what the symptom endpoint requires
##   FLU_steroid()                   # A9  corticosteroid trade
##   FLU_prophylaxis()               # A10 the other side of the peak
##   FLU_lrt_bacterial()             # A11 lower tract and the sequel
##   FLU_immunocompromised()         # A12 when the window never closes
##   FLU_sensitivity()               # A13 local sensitivity
##   FLU_trial_ledger()              # model vs published, in one table
##############################################################################

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
})

##############################################################################
## 1. THE MODEL
##############################################################################
##
## UNITS.  Mixed-unit in-host convention (Baccam 2006 J Virol 80:7590,
## PMID 16840338): cell compartments are absolute cell counts, virus is a
## titre in TCID50/mL of nasal-wash (URT) or BAL (LRT) equivalent, and the
## infectivity constant beta absorbs the sampling volume.  Time is in DAYS.
##
## Two consequences of that convention are load-bearing and easy to get wrong:
##
##  (a) the absorption term -beta*T*V has units of cells/day and CANNOT be
##      subtracted from a titre.  It is omitted here, as in every published
##      human influenza fit; loss of virions to cell entry is absorbed into c.
##
##  (b) beta and p are NOT separately identifiable from titre data.  Only the
##      product beta*T0*p is fixed by the observed growth rate.  p is then set
##      by the observed titre SCALE.  Neither should be quoted alone.

flu_code <- '
$PROB
Influenza A: a replication loop, the operators that act on each of its terms,
and the residual viral AUC that bounds what any of them can achieve.

$PARAM @annotated
// ------------------- URT viral dynamics -------------------
T0U     : 4.00e8  : URT epithelial target cells (cells)             [LIT]
BETA    : 1.63e-5 : infectivity (mL/TCID50/d)                       [CAL]
KECL    : 4.0     : eclipse exit rate, 1/d (6 h)                    [LIT]
DELTA   : 2.0     : productive-cell death rate, 1/d                 [CAL]
PVIR    : 0.114   : virion production, TCID50/mL per cell per d     [CAL]
CVIR    : 2.4     : free-virion clearance, 1/d                      [CAL]
LREG    : 0.10    : epithelial resupply, 1/d                        [ASM]
EXTC    : 1.0     : extinction floor, infected cells                [ASM]

// ------------------- LRT viral dynamics -------------------
T0L     : 1.00e9  : accessible LRT epithelial targets (cells)       [ASM]
BETAL   : 1.05e-6 : LRT infectivity (mucus, surfactant, SP-D)       [ASM]
PVIRL   : 0.018   : LRT virion production                           [ASM]
DELTAL  : 1.6     : LRT productive-cell death, 1/d                  [ASM]
LREGL   : 0.055   : alveolar repair, 1/d (slower than nasal)        [ASM]
ADESC   : 0.008   : URT -> LRT descent / aspiration flux, 1/d       [ASM]

// ------------------- interferon / refractory state --------
QF      : 3.1e-8  : IFN induction per infected cell, 1/d            [CAL]
DF      : 2.0     : IFN decay, 1/d                                  [LIT]
PHIF    : 2.6     : max T -> R conversion rate, 1/d                 [CAL]
KF      : 0.45    : IFN EC50 for the antiviral state                [ASM]
RHOR    : 0.35    : R -> T reversion, 1/d                           [ASM]

// ------------------- adaptive immunity --------------------
QAG     : 6.0     : antigen appearance per infected fraction        [ASM]
DAG     : 1.2     : antigen decay, 1/d                              [ASM]
RCTL    : 3.1     : CD8 expansion rate, 1/d                         [CAL]
KAG     : 0.09    : antigen EC50 for CD8 expansion                  [ASM]
CTLMAX  : 120.0   : CD8 carrying capacity (fold over naive)         [ASM]
DCTL    : 0.32    : CD8 contraction, 1/d                            [ASM]
KKILL   : 0.003   : CTL killing of infected cells, 1/d per unit     [CAL]
QPC     : 8.0     : plasma-cell recruitment, 1/d                    [CAL]
DPC     : 0.25    : plasma-cell loss, 1/d                           [ASM]
KABP    : 3.0     : IgG secretion per plasma cell, 1/d              [ASM]
DAB     : 0.05    : IgG decay, 1/d (t1/2 ~14 d)                     [LIT]
KIGA    : 0.30    : mucosal IgA secretion, 1/d                      [ASM]
DIGA    : 0.03    : IgA decay, 1/d                                  [ASM]
KNEUT   : 0.9     : antibody augmentation of virion clearance       [ASM]
KABENT  : 0.30    : antibody block of entry, per unit               [ASM]
KIGAENT : 0.55    : IgA block of entry, per unit                    [ASM]
QNEU    : 2.4     : neutrophil recruitment per unit LRT damage rate [ASM]
DNEU    : 0.9     : neutrophil clearance, 1/d                       [ASM]

// ------------------- cytokines, symptoms, fever -----------
Q6F     : 200.0   : IL-6 production per unit URT IFN, pg/mL/d       [CAL]
Q6L     : 400.0   : IL-6 production per unit LRT IFN, pg/mL/d       [ASM]
Q6D     : 9.0e-9  : IL-6 per LRT cell killed per day                [ASM]
D6      : 5.5     : IL-6 elimination, 1/d (t1/2 ~3 h)               [LIT]
IL60    : 1.6     : healthy baseline IL-6, pg/mL                    [LIT]
KON     : 3.2     : symptom onset rate, 1/d                         [CAL]
KOFF    : 2.0     : symptom resolution rate, 1/d                    [CAL]
K6      : 11.0    : IL-6 EC50 for the symptom drive, pg/mL          [CAL]
SMAX    : 21.0    : max composite score (7 symptoms x 0-3)          [LIT]
SALLEV  : 7.0     : alleviation ceiling (all 7 symptoms <= 1)       [LIT]
SONSET  : 2.0     : score at which a patient calls themselves ill   [ASM]
WVIR    : 0.60    : symptom drive fraction tracking titre (see A8)  [CAL]
LODS    : 0.5     : titre below which the viral drive is zero       [ASM]
LREFS   : 6.0     : titre span over which the viral drive saturates [ASM]
KTMP    : 8.5     : fever gain, degC/d                              [CAL]
KTMPD   : 2.6     : defervescence, 1/d                              [ASM]
KT6     : 30.0    : IL-6 EC50 for fever, pg/mL                      [ASM]

// ------------------- bacterial superinfection -------------
BAC0    : 2.0     : baseline log10 CFU nasopharyngeal carriage      [ASM]
BMAXB   : 9.0     : log10 CFU carrying capacity                     [ASM]
KBG     : 1.05    : bacterial growth, log10/d at full permissiveness[ASM]
KBADH   : 6.0     : adhesion gain from exposed sialic acid / damage [ASM]
KBCLR   : 0.81667 : clearance; FIXED by KBG*(1-BAC0/BMAXB)          [derived]
KBIFN   : 0.75    : IFN suppression of antibacterial defence        [LIT]
BACTHR  : 6.0     : log10 CFU stipulated as superinfection          [ASM]

// ------------------- mutation and fitness -----------------
MU      : 2.5e-5  : per-site per-replication mutation rate          [LIT]
NU      : 2.5e-5  : back mutation                                   [LIT]
COST    : 0.18    : replicative fitness cost carried by the mutant  [LIT]

// ------------------- oseltamivir PK -----------------------
KA_OS   : 12.0    : absorption, 1/d                                 [LIT]
KCONV   : 26.0    : CES1 prodrug hydrolysis, 1/d                    [LIT]
CLOP    : 100.0   : prodrug non-converting clearance, L/d           [LIT]
VOP     : 30.0    : prodrug volume, L                               [LIT]
CLOC    : 451.2   : carboxylate clearance, L/d (18.8 L/h)           [LIT]
VOC     : 28.0    : carboxylate central volume, L                   [LIT]
VOC2    : 30.0    : carboxylate peripheral volume, L                [ASM]
QOC     : 60.0    : intercompartmental clearance, L/d               [ASM]
FOS     : 0.80    : fraction of dose converted to carboxylate       [LIT]
RELF_OC : 1.00    : ELF : plasma ratio for the carboxylate          [LIT]
KEQ_OC  : 24.0    : ELF equilibration, 1/d                          [ASM]

// ------------------- baloxavir PK -------------------------
KA_BX   : 19.2    : absorption, 1/d (Tmax ~4 h)                     [LIT]
CLBX    : 149.8   : CL/F, L/d (6.24 L/h)                            [LIT]
VBX     : 380.0   : central volume, L                               [CAL]
VBX2    : 420.0   : peripheral volume, L                            [CAL]
QBX     : 500.0   : intercompartmental clearance, L/d               [ASM]
FUBX    : 0.07    : free fraction (93% protein bound)               [LIT]
RELF_BX : 1.00    : free drug equilibrates with ELF                 [ASM]
KEQ_BX  : 24.0    : ELF equilibration, 1/d                          [ASM]

// ------------------- other drug PK ------------------------
CLPR    : 264.0   : peramivir clearance, L/d (11 L/h)               [LIT]
VPR     : 12.6    : peramivir Vss, L                                [LIT]
VPR2    : 8.0     : peramivir peripheral volume, L                  [ASM]
QPR     : 40.0    : peramivir intercompartmental clearance, L/d     [ASM]
KA_FV   : 24.0    : favipiravir absorption, 1/d                     [LIT]
VFV     : 15.0    : favipiravir volume, L                           [LIT]
VMFV    : 9000.0  : favipiravir Vmax, mg/d (auto-inhibited AO)      [LIT]
KMFV    : 45.0    : favipiravir Km, mg/L                            [LIT]
CLMB    : 0.25    : mAb clearance, L/d (t1/2 ~21 d)                 [LIT]
VMB     : 3.2     : mAb central volume, L                           [LIT]
VMB2    : 3.0     : mAb peripheral volume, L                        [ASM]
QMB     : 0.6     : mAb intercompartmental clearance, L/d           [ASM]
KA_ST   : 24.0    : dexamethasone absorption, 1/d                   [LIT]
CLST    : 200.0   : dexamethasone clearance, L/d                    [LIT]
VST     : 60.0    : dexamethasone volume, L                         [LIT]

// ------------------- pharmacodynamics ---------------------
// These EC50s are IN-VIVO EFFECTIVE values calibrated against CAPSTONE-1,
// NOT in-vitro values.  The gap is reported, not hidden: see FLU_pd_calibration.
EC50_NAI : 3.0    : carboxylate ELF conc for 50% release block, ng/mL [CAL]
EMAX_NAI : 0.995  : max fractional blockade of release               [CAL]
HILL_NAI : 1.0    : Hill slope, NAI                                  [ASM]
FBETA_NAI: 0.30   : NAI entry-block as a fraction of its release block[ASM]
EC50_BX  : 0.052  : free baloxavir acid for 50% block, ng/mL         [CAL]
EMAX_BX  : 0.9999 : max fractional blockade                          [CAL]
HILL_BX  : 2.0    : Hill slope -- REQUIRED by the data, see A2        [CAL]
ETA_BX   : 1.00   : transcription block also blocks E -> I fully      [ASM]
FMATU    : 1.00   : idealised productive-fraction operator (A5 only)  [ASM]
EC50_PR  : 3.0    : peramivir shares the NAI effective EC50           [ASM]
EMAX_PR  : 0.995  : peramivir Emax                                    [ASM]
EC50_FV  : 25.0   : favipiravir, mg/L                                 [ASM]
EMAX_FV  : 0.93   : favipiravir Emax                                  [ASM]
FMUTA_FV : 0.55   : fraction of favipiravir progeny made non-viable   [LIT]
KMAB_C   : 0.55   : mAb-driven virion clearance, 1/d per ug/mL        [ASM]
KMAB_B   : 0.020  : mAb entry blockade, per ug/mL                     [ASM]
ECST     : 12.0   : dexamethasone conc for half-max suppression, ng/mL[ASM]
IMAX_ST  : 0.65   : max suppression of IFN induction and CTL killing  [LIT]

// ------------------- resistance profile (set per scenario) -
RF_NAI  : 1.0     : fold EC50 shift carried by the mutant, NAIs
RF_BX   : 1.0     : fold EC50 shift carried by the mutant, baloxavir
RF_FV   : 1.0     : fold EC50 shift carried by the mutant, favipiravir

// ------------------- host phenotype -----------------------
FIFN    : 1.0     : multiplier on IFN induction
FCTL    : 1.0     : multiplier on the CD8 response
FAB     : 1.0     : multiplier on the humoral response
IGA0    : 0.0     : pre-existing mucosal IgA (vaccination / prior infection)
CTL0    : 1.0     : pre-existing cross-reactive CD8 memory
V0      : 0.01    : inoculum, TCID50/mL equivalent
VM0FRAC : 0.0     : pre-existing resistant fraction of the inoculum

$CMT @annotated
// --- upper respiratory tract ---
TU   : susceptible epithelial cells, URT (cells)
RU   : refractory (IFN-protected) cells, URT (cells)
EW   : eclipse cells, wild type, URT (cells)
IW   : productively infected, wild type, URT (cells)
EM   : eclipse cells, resistant mutant, URT (cells)
IM   : productively infected, mutant, URT (cells)
VW   : wild-type virus, URT (TCID50/mL)
VM   : mutant virus, URT (TCID50/mL)
DU   : dead / denuded epithelium, URT (cells)
// --- lower respiratory tract ---
TL   : susceptible epithelial cells, LRT (cells)
RL   : refractory cells, LRT (cells)
ELW  : eclipse, wild type, LRT (cells)
ILW  : productively infected, wild type, LRT (cells)
ELM  : eclipse, mutant, LRT (cells)
ILM  : productively infected, mutant, LRT (cells)
VLW  : wild-type virus, LRT (TCID50/mL)
VLM  : mutant virus, LRT (TCID50/mL)
DL   : dead / denuded alveolar epithelium, LRT (cells)
// --- immune ---
FU   : type I interferon, URT (arbitrary)
FL   : type I interferon, LRT (arbitrary)
AG   : antigen load in the draining node (arbitrary)
CTL  : effector CD8 T cells (fold over naive)
PC   : antibody-secreting plasma cells (arbitrary)
AB   : serum neutralising IgG (arbitrary)
IGA  : mucosal secretory IgA (arbitrary)
NEU  : neutrophil infiltrate, LRT (arbitrary)
IL6  : systemic pro-inflammatory cytokine, IL-6 proxy (pg/mL)
// --- clinical ---
SYM  : composite symptom score (0-21)
TMP  : temperature elevation above 37.0 (degC)
BAC  : bacterial burden (log10 CFU)
// --- drug PK ---
OSd  : oseltamivir phosphate, gut depot (mg base)
OSp  : oseltamivir prodrug, plasma (mg)
OCc  : oseltamivir carboxylate, central (mg)
OCp  : oseltamivir carboxylate, peripheral (mg)
OCe  : oseltamivir carboxylate, ELF / effect site (ng/mL)
BXd  : baloxavir marboxil, gut depot (mg)
BXc  : baloxavir acid, central (mg)
BXp  : baloxavir acid, peripheral (mg)
BXe  : baloxavir acid, free, ELF / effect site (ng/mL)
PRc  : peramivir, central (mg)
PRp  : peramivir, peripheral (mg)
FVd  : favipiravir, gut depot (mg)
FVc  : favipiravir, central (mg)
MBc  : anti-HA monoclonal antibody, central (mg)
MBp  : anti-HA monoclonal antibody, peripheral (mg)
STd  : dexamethasone, gut depot (mg)
STc  : dexamethasone, central (mg)
// --- bookkeeping ---
AUCU : cumulative URT log-titre AUC (log10.d)
AUCL : cumulative LRT log-titre AUC (log10.d)
CUMK : cumulative epithelial cells killed, URT (cells)

$MAIN
TU_0   = T0U;
TL_0   = T0L;
VW_0   = V0 * (1.0 - VM0FRAC);
VM_0   = V0 * VM0FRAC;
CTL_0  = CTL0;
IGA_0  = IGA0;
IL6_0  = IL60;
BAC_0  = BAC0;

$ODE
// ---------------- drug concentrations ----------------------------------
double C_OC  = OCc / VOC * 1000.0;      // ng/mL
double C_OCE = OCe;                     // ng/mL, ELF
double C_BX  = BXc / VBX * 1000.0;      // ng/mL total
double C_BXE = BXe;                     // ng/mL free, ELF
double C_PR  = PRc / VPR * 1000.0;      // ng/mL
double C_FV  = FVc / VFV;               // mg/L
double C_MB  = MBc / VMB;               // mg/L = ug/mL
double C_ST  = STc / VST * 1000.0;      // ng/mL

// ---------------- pharmacodynamic operators ----------------------------
// RELEASE blockade, neuraminidase inhibitors.  The wild type sees EC50 as
// given; the mutant sees it shifted by RF_NAI.
double eo_W = EMAX_NAI * pow(C_OCE, HILL_NAI) /
              (pow(C_OCE, HILL_NAI) + pow(EC50_NAI, HILL_NAI) + 1e-30);
double eo_M = EMAX_NAI * pow(C_OCE, HILL_NAI) /
              (pow(C_OCE, HILL_NAI) + pow(EC50_NAI * RF_NAI, HILL_NAI) + 1e-30);
double ep_W = EMAX_PR * C_PR / (C_PR + EC50_PR + 1e-30);
double ep_M = EMAX_PR * C_PR / (C_PR + EC50_PR * RF_NAI + 1e-30);
double e_nai_W = eo_W + ep_W; if (e_nai_W > EMAX_NAI) e_nai_W = EMAX_NAI;
double e_nai_M = eo_M + ep_M; if (e_nai_M > EMAX_NAI) e_nai_M = EMAX_NAI;

// TRANSCRIPTION blockade, baloxavir (PA endonuclease / cap-snatching).
// The Hill slope is not cosmetic: with HILL_BX = 1 the residual production
// fraction is bounded below by EC50/(C+EC50) and the model cannot fall faster
// than about -3.3 log10 in 24 h at ANY Emax.  CAPSTONE-1 reports -4.8.
double e_bx_W = EMAX_BX * pow(C_BXE, HILL_BX) /
                (pow(C_BXE, HILL_BX) + pow(EC50_BX, HILL_BX) + 1e-30);
double e_bx_M = EMAX_BX * pow(C_BXE, HILL_BX) /
                (pow(C_BXE, HILL_BX) + pow(EC50_BX * RF_BX, HILL_BX) + 1e-30);

// POLYMERASE blockade + lethal mutagenesis, favipiravir
double e_fv_W = EMAX_FV * C_FV / (C_FV + EC50_FV + 1e-30);
double e_fv_M = EMAX_FV * C_FV / (C_FV + EC50_FV * RF_FV + 1e-30);

// combined multiplicative blockade of virion OUTPUT per infected cell
double prod_W = (1.0 - e_nai_W) * (1.0 - e_bx_W) * (1.0 - e_fv_W);
double prod_M = (1.0 - e_nai_M) * (1.0 - e_bx_M) * (1.0 - e_fv_M);

// blockade of the eclipse -> productive transition (transcription only)
double matu_W = FMATU * (1.0 - ETA_BX * e_bx_W);
double matu_M = FMATU * (1.0 - ETA_BX * e_bx_M);

// fraction of released progeny that remain infectious (favipiravir)
double inf_W = 1.0 - FMUTA_FV * e_fv_W;
double inf_M = 1.0 - FMUTA_FV * e_fv_M;

// ENTRY blockade: NAI (weak), IgG, mucosal IgA, mAb
double ent_blk = FBETA_NAI * e_nai_W + KABENT * AB + KIGAENT * IGA
                 + KMAB_B * C_MB;
if (ent_blk > 0.995) ent_blk = 0.995;
double ent = 1.0 - ent_blk;

// VIRION CLEARANCE augmentation: antibody + mAb
double cV = CVIR * (1.0 + KNEUT * AB) + KMAB_C * C_MB;

// IMMUNOPATHOLOGY DAMPING: corticosteroid suppresses IFN induction and killing
double supp      = IMAX_ST * C_ST / (C_ST + ECST + 1e-30);
double ifn_gain  = FIFN * (1.0 - supp);
double kill_gain = FCTL * (1.0 - supp);

// EXTINCTION FLOOR.  A deterministic ODE lets a strain fall to 1e-20 cells
// and grow back when targets regenerate.  Real infections at that burden are
// extinct.  Each strain is scaled by a smooth establishment factor that is
// ~1 above one infected cell and -> 0 below it.  Without it the resupply term
// manufactures a spurious second wave; with it, clearance is final.
double xiW = (EW + IW + ELW + ILW) / ((EW + IW + ELW + ILW) + EXTC);
double xiM = (EM + IM + ELM + ILM) / ((EM + IM + ELM + ILM) + EXTC);

// ---------------- URT viral dynamics -----------------------------------
double infW = BETA * ent * TU * VW;
double infM = BETA * ent * TU * VM;
double Itot = IW + IM;
double ctl_kill = KKILL * kill_gain * CTL;
double to_ref   = PHIF * FU / (FU + KF) * TU;
double from_ref = RHOR * RU;
double occupied = TU + RU + EW + IW + EM + IM;
double regen    = LREG * (T0U - occupied > 0 ? T0U - occupied : 0.0);

dxdt_TU = -infW - infM - to_ref + from_ref + regen;
dxdt_RU = to_ref - from_ref;

// An eclipse cell ALWAYS leaves the eclipse state at rate KECL.  A
// transcription blocker does not freeze it there -- it decides where the cell
// goes: a fraction matu becomes productive, the rest is abortively infected
// and counted as lost epithelium.  Writing it the other way (blocking the
// exit) creates an eclipse reservoir that keeps seeding productive cells for
// days after the drug has taken hold, and the model then cannot reproduce
// baloxavir's 24-hour titre fall at any potency.
dxdt_EW = infW - KECL * EW - DELTA * 0.25 * EW;
dxdt_IW = KECL * matu_W * EW - (DELTA + ctl_kill) * IW;
dxdt_EM = infM - KECL * EM - DELTA * 0.25 * EM;
dxdt_IM = KECL * matu_M * EM - (DELTA + ctl_kill) * IM;
double abortU = KECL * ((1.0 - matu_W) * EW + (1.0 - matu_M) * EM);

double relW = PVIR * prod_W * IW * xiW;
double relM = PVIR * prod_M * (1.0 - COST) * IM * xiM;

// NOTE ON UNITS: the absorption term -beta*T*V has units of cells/day and is
// NOT subtracted from a titre.  See the header.
dxdt_VW = relW * (1.0 - MU) * inf_W + relM * NU * inf_M - cV * VW - ADESC * VW;
dxdt_VM = relM * (1.0 - NU) * inf_M + relW * MU * inf_W - cV * VM - ADESC * VM;

double killedU = (DELTA + ctl_kill) * Itot + abortU;
dxdt_DU   = killedU - LREG * DU;
dxdt_CUMK = killedU;

// ---------------- LRT viral dynamics -----------------------------------
double infLW = BETAL * ent * TL * VLW;
double infLM = BETAL * ent * TL * VLM;
double ILtot = ILW + ILM;
double to_refL   = PHIF * FL / (FL + KF) * TL;
double from_refL = RHOR * RL;
double occupiedL = TL + RL + ELW + ILW + ELM + ILM;
double regenL    = LREGL * (T0L - occupiedL > 0 ? T0L - occupiedL : 0.0);

dxdt_TL  = -infLW - infLM - to_refL + from_refL + regenL;
dxdt_RL  = to_refL - from_refL;
dxdt_ELW = infLW - KECL * ELW - DELTAL * 0.25 * ELW;
dxdt_ILW = KECL * matu_W * ELW - (DELTAL + ctl_kill) * ILW;
dxdt_ELM = infLM - KECL * ELM - DELTAL * 0.25 * ELM;
dxdt_ILM = KECL * matu_M * ELM - (DELTAL + ctl_kill) * ILM;
double abortL = KECL * ((1.0 - matu_W) * ELW + (1.0 - matu_M) * ELM);

double relLW = PVIRL * prod_W * ILW * xiW;
double relLM = PVIRL * prod_M * (1.0 - COST) * ILM * xiM;
dxdt_VLW = relLW * (1.0 - MU) * inf_W + relLM * NU * inf_M - cV * VLW + ADESC * VW;
dxdt_VLM = relLM * (1.0 - NU) * inf_M + relLW * MU * inf_W - cV * VLM + ADESC * VM;

double killedL = (DELTAL + ctl_kill) * ILtot + abortL;
dxdt_DL = killedL - LREGL * DL;

// ---------------- interferon -------------------------------------------
dxdt_FU = QF * ifn_gain * Itot  - DF * FU;
dxdt_FL = QF * ifn_gain * ILtot - DF * FL;

// ---------------- adaptive ---------------------------------------------
double ag_drive = Itot / T0U + ILtot / T0L;
dxdt_AG  = QAG * ag_drive - DAG * AG;
dxdt_CTL = RCTL * AG / (AG + KAG) * CTL * (1.0 - CTL / CTLMAX)
           - DCTL * (CTL - CTL0);
dxdt_PC  = QPC * FAB * AG - DPC * PC;
dxdt_AB  = KABP * PC - DAB * AB;
dxdt_IGA = KIGA * PC - DIGA * (IGA - IGA0);
dxdt_NEU = QNEU * killedL / T0L - DNEU * NEU;

// ---------------- cytokine, symptoms, fever ----------------------------
dxdt_IL6 = Q6F * FU + Q6L * FL + Q6D * killedL - D6 * (IL6 - IL60);
double drive6 = (IL6 - IL60 > 0.0) ? IL6 - IL60 : 0.0;
double f_cyt  = drive6 / (drive6 + K6);

// WVIR is the fraction of the symptom drive that tracks the INSTANTANEOUS
// titre rather than the cytokine state.  At WVIR = 0 the symptom score is a
// pure function of a cytokine response whose size and timing are settled
// before any post-peak drug can act, and then NO antiviral shortens the
// illness at any potency.  See FLU_symptom_requirement().
double vtot = VW + VM; if (vtot < 1e-6) vtot = 1e-6;
double lv    = log10(vtot);
double f_vir = (lv - LODS) / LREFS;
if (f_vir < 0.0) f_vir = 0.0;
if (f_vir > 1.0) f_vir = 1.0;
double f_sym = (1.0 - WVIR) * f_cyt + WVIR * f_vir;

dxdt_SYM = KON * f_sym * (SMAX - SYM) - KOFF * SYM;
dxdt_TMP = KTMP * drive6 / (drive6 + KT6) - KTMPD * TMP;

// ---------------- bacterial superinfection -----------------------------
// Two influenza-created conditions multiply: exposed adhesion sites on
// denuded epithelium, and interferon suppression of antibacterial defence.
double dam_frac = DL / T0L + DU / T0U; if (dam_frac > 1.0) dam_frac = 1.0;
double permiss  = 1.0 + KBADH * dam_frac;
double neu_fn   = 1.0 / (1.0 + KBIFN * (FU + FL) / KF);
dxdt_BAC = KBG * permiss * (1.0 - BAC / BMAXB) - KBCLR * neu_fn * (BAC / BAC0);

// ---------------- PK ----------------------------------------------------
dxdt_OSd = -KA_OS * OSd;
dxdt_OSp = KA_OS * OSd - (KCONV + CLOP / VOP) * OSp;
double conv = KCONV * OSp * FOS * (284.4 / 312.4);
dxdt_OCc = conv - CLOC / VOC * OCc - QOC * (OCc / VOC - OCp / VOC2);
dxdt_OCp = QOC * (OCc / VOC - OCp / VOC2);
dxdt_OCe = KEQ_OC * (C_OC * RELF_OC - OCe);

dxdt_BXd = -KA_BX * BXd;
dxdt_BXc = KA_BX * BXd - CLBX / VBX * BXc - QBX * (BXc / VBX - BXp / VBX2);
dxdt_BXp = QBX * (BXc / VBX - BXp / VBX2);
dxdt_BXe = KEQ_BX * (C_BX * FUBX * RELF_BX - BXe);

dxdt_PRc = -CLPR / VPR * PRc - QPR * (PRc / VPR - PRp / VPR2);
dxdt_PRp = QPR * (PRc / VPR - PRp / VPR2);

dxdt_FVd = -KA_FV * FVd;
dxdt_FVc = KA_FV * FVd - VMFV * C_FV / (KMFV + C_FV);

dxdt_MBc = -CLMB / VMB * MBc - QMB * (MBc / VMB - MBp / VMB2);
dxdt_MBp = QMB * (MBc / VMB - MBp / VMB2);

dxdt_STd = -KA_ST * STd;
dxdt_STc = KA_ST * STd - CLST / VST * STc;

// ---------------- bookkeeping ------------------------------------------
double vt  = (VW + VM  > 1e-6) ? VW + VM   : 1e-6;
double vlt = (VLW + VLM > 1e-6) ? VLW + VLM : 1e-6;
dxdt_AUCU = log10(vt / 1e-6);
dxdt_AUCL = log10(vlt / 1e-6);

$TABLE
double VTOT   = VW + VM;
double LOGV   = log10(VTOT   > 1e-6 ? VTOT   : 1e-6);
double LOGVW  = log10(VW     > 1e-6 ? VW     : 1e-6);
double LOGVM  = log10(VM     > 1e-6 ? VM     : 1e-6);
double LOGVL  = log10(VLW+VLM> 1e-6 ? VLW+VLM: 1e-6);
// Mutant share is defined only where the total titre is still detectable.
// Below the assay floor the ratio is arithmetically valid and clinically
// empty, and reporting it there manufactures 100% mutant readings out of
// virus nobody could sample.
double MUTFR  = (VTOT > pow(10.0, LODS)) ? VM / (VTOT + 1e-30) : 0.0;
double TFRAC  = TU / T0U;
double EPILOST= CUMK / T0U;
double TEMPC  = 37.0 + TMP;
double DAMFR  = DL / T0L + DU / T0U;
double COCE   = OCe;
double CBXE   = BXe;
double CBXTOT = BXc / VBX * 1000.0;
double SPO2   = 98.0 - 22.0 * pow(DAMFR, 2.0) / (pow(DAMFR, 2.0) + pow(0.35, 2.0));
double ALLEV  = (SYM <= SALLEV) ? 1.0 : 0.0;
double SUPER  = (BAC >= BACTHR) ? 1.0 : 0.0;

$CAPTURE LOGV LOGVW LOGVM LOGVL MUTFR TFRAC EPILOST TEMPC DAMFR
         COCE CBXE CBXTOT SPO2 ALLEV SUPER
'

flu <- mcode_cache("influenza_qsp", flu_code)

##############################################################################
## 2. CONSTANTS, DOSING AND READOUTS
##############################################################################

LOD        <- 0.5      # log10 TCID50/mL assay floor (CAPSTONE-1)
RX_DELAY_H <- 24       # median time from symptom onset to first dose

## Resistance profiles.  These are switched on per scenario; the DEFAULT model
## carries no resistance shift, because I38T emerges in roughly one treated
## adult in ten and the ledger should describe the other nine.
I38T  <- list(RF_BX = 50,  RF_NAI = 1,   COST = 0.18)   # baloxavir, PA/I38T
H275Y <- list(RF_NAI = 300, RF_BX = 1,   COST = 0.28)   # oseltamivir, NA/H275Y

CMT_NUM <- function(name) which(mrgsolve::cmt(flu) == name)

.dose <- function(time, cmt, amt) {
  data.frame(ID = 1, time = time, cmt = CMT_NUM(cmt), amt = amt, evid = 1)
}

rx_oseltamivir <- function(t0, days = 5, mg = 75)
  .dose(t0 + seq(0, by = 0.5, length.out = days * 2), "OSd", mg)
rx_baloxavir   <- function(t0, mg = 40)  .dose(t0, "BXd", mg)
rx_peramivir   <- function(t0, mg = 600) .dose(t0, "PRc", mg)
rx_favipiravir <- function(t0, days = 5)
  rbind(.dose(c(t0, t0 + 0.5), "FVd", 1800),
        .dose(t0 + seq(1, by = 0.5, length.out = days * 2 - 2), "FVd", 800))
rx_mab         <- function(t0, mg = 3600) .dose(t0, "MBc", mg)
rx_dexamethasone <- function(t0, days = 5, mg = 6)
  .dose(t0 + seq(0, by = 1, length.out = days), "STd", mg)
## mrgsim_d needs at least one record; a zero-amount, zero-evid row is the
## idiomatic way to run an untreated arm through the same code path.
rx_none        <- function() data.frame(ID = 1, time = 0, cmt = 1, amt = 0, evid = 0)

sim <- function(doses = NULL, tmax = 30, delta = 1/96, ...) {
  m <- flu
  ov <- list(...)
  if (length(ov)) m <- param(m, ov)
  ev <- if (is.null(doses) || nrow(doses) == 0) rx_none() else doses
  as.data.frame(mrgsim_d(m, data = ev, end = tmax, delta = delta,
                         atol = 1e-10, rtol = 1e-8, maxsteps = 1e6))
}

symptom_onset <- function(out, SONSET = 2.0) {
  i <- which(out$SYM >= SONSET)
  if (!length(i)) return(NA_real_)
  out$time[i[1]]
}

placebo_onset <- function(...) symptom_onset(sim(NULL, tmax = 30, ...))

## CAPSTONE alleviation rule: all seven symptoms mild or absent, SUSTAINED for
## 21.5 h.  Implemented literally rather than as "first time below threshold".
time_to_alleviation <- function(out, t_start, SALLEV = 7.0) {
  dt   <- out$time[2] - out$time[1]
  need <- round((21.5 / 24) / dt)
  ok   <- out$SYM <= SALLEV
  lo   <- which(out$time >= t_start)[1]
  run  <- 0
  for (i in seq(lo, nrow(out))) {
    run <- if (ok[i]) run + 1 else 0
    if (run >= need) return((out$time[i - need + 1] - t_start) * 24)
  }
  NA_real_
}

time_to_shed_stop <- function(out, t_start, thr = LOD) {
  lo <- which(out$time >= t_start)[1]
  for (i in seq(lo, nrow(out)))
    if (out$LOGV[i] <= thr && all(out$LOGV[i:nrow(out)] <= thr))
      return((out$time[i] - t_start) * 24)
  NA_real_
}

titre_drop_at <- function(out, t_start, hours) {
  i0 <- which.min(abs(out$time - t_start))
  i1 <- which.min(abs(out$time - (t_start + hours / 24)))
  out$LOGV[i1] - out$LOGV[i0]
}

viral_auc <- function(out, thr = LOD) {
  y <- pmax(out$LOGV - thr, 0)
  sum((head(y, -1) + tail(y, -1)) / 2 * diff(out$time))
}

## R(t): the residual AUC, the EXACT upper bound on what a drug started at t
## can remove.  This is the single most important readout in the model.
residual_auc <- function(out, thr = LOD) {
  y  <- pmax(out$LOGV - thr, 0)
  dt <- out$time[2] - out$time[1]
  rev(cumsum(rev(y))) * dt
}

fever_duration <- function(out, thr = 0.7)
  sum(out$TMP >= thr) * (out$time[2] - out$time[1]) * 24

summarise_run <- function(out, t_ref) {
  list(
    peak_log = max(out$LOGV),
    t_peak   = out$time[which.max(out$LOGV)],
    auc_log  = viral_auc(out),
    drop24   = titre_drop_at(out, t_ref, 24),
    shed_h   = time_to_shed_stop(out, t_ref),
    ttas_h   = time_to_alleviation(out, t_ref),
    fever_h  = fever_duration(out),
    peak_sym = max(out$SYM),
    peak_il6 = max(out$IL6),
    lrt_peak = max(out$LOGVL),
    mut_peak = max(out$LOGVM[out$time >= t_ref]),
    mut_frac = max(out$MUTFR[out$time >= t_ref]),
    epi_lost = max(out$EPILOST),
    tmin     = min(out$TFRAC),
    bac_peak = max(out$BAC),
    spo2_min = min(out$SPO2)
  )
}

##############################################################################
## 3. SCENARIOS
##############################################################################
## Each entry is (dose table, parameter overrides, reference time).  The
## reference time is the arm's OWN first dose, so that a timing variant is
## scored from when its drug was actually given rather than from a shared
## randomisation clock.

flu_scenarios <- function(t_rx) {
  onset <- t_rx - RX_DELAY_H / 24
  list(
    "1. Placebo"                        = list(NULL, list(), t_rx),
    "2. Oseltamivir 75 mg BID x5d"      = list(rx_oseltamivir(t_rx), list(), t_rx),
    "3. Baloxavir 40 mg single"         = list(rx_baloxavir(t_rx), list(), t_rx),
    "4. Peramivir 600 mg IV single"     = list(rx_peramivir(t_rx), list(), t_rx),
    "5. Favipiravir 1800/800 BID x5d"   = list(rx_favipiravir(t_rx), list(), t_rx),
    "6. Baloxavir + oseltamivir"        = list(rbind(rx_baloxavir(t_rx),
                                                     rx_oseltamivir(t_rx)),
                                               list(), t_rx),
    "7. Anti-HA mAb single IV"          = list(rx_mab(t_rx), list(), t_rx),
    "8. Baloxavir at symptom onset"     = list(rx_baloxavir(onset), list(), onset),
    "9. Baloxavir late (72 h)"          = list(rx_baloxavir(onset + 3), list(),
                                               onset + 3),
    "10. Oseltamivir + dexamethasone"   = list(rbind(rx_oseltamivir(t_rx),
                                                     rx_dexamethasone(t_rx)),
                                               list(), t_rx),
    "11. Baloxavir, immunocompromised"  = list(rx_baloxavir(t_rx),
                                               list(FCTL = 0.12, FAB = 0.15,
                                                    FIFN = 0.55), t_rx),
    "12. Baloxavir, PA/I38T profile"    = list(rx_baloxavir(t_rx), I38T, t_rx),
    "13. Oseltamivir, H275Y profile"    = list(rx_oseltamivir(t_rx), H275Y, t_rx),
    "14. Vaccinated host, no drug"      = list(NULL, list(IGA0 = 0.55,
                                                          CTL0 = 6.0), t_rx)
  )
}

##############################################################################
## 4. ANALYSES
##   Every number in the comments below is the output of the corresponding
##   function in flu_reference_check.py, reproduced in flu_reference_output.txt.
##############################################################################

## --- A0 -------------------------------------------------------------------
## Untreated natural history.  Reference values:
##   peak URT titre 6.62 log10 at 51.5 h    (target 6.0-7.0 at 48-72 h)
##   shedding above LOD 6.3 d               (target 4.0-5.5 d -- LONG, reported)
##   symptom onset 24.2 h, peak score 12.4/21 at 62.0 h
##   peak temperature 38.99 C, fever 67 h, peak IL-6 60 pg/mL
##   min susceptible pool 0.1% of T0        <- the infection IS target-limited
##   peak LRT titre 3.87, peak CD8 99x naive
FLU_baseline <- function() {
  out <- sim(NULL, tmax = 30)
  data.frame(
    quantity = c("peak URT titre (log10)", "time to peak (h)",
                 "shedding >LOD (d)", "symptom onset (h)",
                 "peak symptom score", "time to peak symptoms (h)",
                 "peak temperature (C)", "fever duration (h)",
                 "peak IL-6 (pg/mL)", "min susceptible pool (%)",
                 "peak LRT titre (log10)", "peak CD8 (x naive)"),
    model = c(max(out$LOGV), out$time[which.max(out$LOGV)] * 24,
              sum(out$LOGV > LOD) * (out$time[2] - out$time[1]),
              symptom_onset(out) * 24, max(out$SYM),
              out$time[which.max(out$SYM)] * 24, max(out$TEMPC),
              fever_duration(out), max(out$IL6), min(out$TFRAC) * 100,
              max(out$LOGVL), max(out$CTL)),
    target = c("6.0-7.0", "48-72", "4.0-5.5", "24-48", "10-15", "48-72",
               "38.5-39.5", "48-96", "20-100", "<20", "2-4", "20-100")
  )
}

## --- A1 -------------------------------------------------------------------
## The bound.  At symptom onset R = 98% of the total viral AUC; at +24 h it is
## 75%; at +48 h, 46%.  The 48-hour rule is the interval over which the
## quantity a drug could act on falls by roughly an order of magnitude, and it
## would still be there if the drug were infinitely potent.
FLU_residual_bound <- function(hours = c(0, 12, 24, 36, 48, 72, 96)) {
  base  <- sim(NULL, tmax = 30)
  ons   <- symptom_onset(base)
  R     <- residual_auc(base)
  total <- R[1]
  do.call(rbind, lapply(hours, function(h) {
    trx <- ons + h / 24
    i   <- which.min(abs(base$time - trx))
    o2  <- sim(rx_baloxavir(trx), tmax = 30)
    ach <- total - viral_auc(o2)
    data.frame(h_post_onset = h, T_left_pct = base$TFRAC[i] * 100,
               bound_R = R[i], pct_of_total = R[i] / total * 100,
               achieved = ach, efficiency_pct = ach / R[i] * 100)
  }))
}

## --- A2 -------------------------------------------------------------------
## DISCREPANCY 1.  The CAPSTONE-1 day-2 titre fall cannot be produced by a
## Michaelis (Hill = 1) concentration-response at ANY Emax: with h = 1 the
## residual production fraction is bounded below by EC50/(C+EC50), so the fall
## saturates near -3.3 log10 in 24 h.  The published value is -4.8.  The trial
## datum is therefore evidence about the SHAPE of the concentration-response,
## not only about potency, and a model fitted on Emax alone would have
## absorbed a shape error into a potency parameter.
## Alongside it, the calibrated in-vivo EC50 is 13x LOWER than the in-vitro
## EC50 on a free-drug basis (0.11 nM vs 1.4 nM).
FLU_pd_calibration <- function() {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  os  <- sim(rx_oseltamivir(trx), tmax = 30)
  bx  <- sim(rx_baloxavir(trx),   tmax = 30)
  p   <- as.list(param(flu))
  MW_BXA <- 483.9; ec50_invitro_nM <- 1.4
  ec50_invitro_ng <- ec50_invitro_nM * MW_BXA / 1000
  e_bx <- p$EMAX_BX * bx$CBXE^p$HILL_BX / (bx$CBXE^p$HILL_BX + p$EC50_BX^p$HILL_BX)
  list(
    oc_elf_cmax_ng_mL      = max(os$COCE),
    oc_ec50_invivo_ng_mL   = p$EC50_NAI,
    bxa_total_cmax_ng_mL   = max(bx$CBXTOT),
    bxa_free_cmax_ng_mL    = max(bx$CBXE),
    bxa_free_cmax_nM       = max(bx$CBXE) / MW_BXA * 1000,
    bxa_ec50_invitro_ng_mL = ec50_invitro_ng,
    bxa_ec50_invivo_ng_mL  = p$EC50_BX,
    fold_invitro_over_invivo = ec50_invitro_ng / p$EC50_BX,
    hill_slope             = p$HILL_BX,
    peak_block_pct         = max(e_bx) * 100,
    residual_production    = 1 - max(e_bx)
  )
}

## --- A3 -------------------------------------------------------------------
## The ledger.  Model vs CAPSTONE-1 (Hayden 2018 NEJM 379:913, PMID 30184455):
##   TTAS   published  80.2 / 53.8 / 53.7 h    model  78.7 / 63.5 / 44.0
##   shed   published  96.0 / 72.0 / 24.0 h    model 118.0 / 85.2 / 47.0
##   day-2  published  -1.3 / -2.8 / -4.8 log  model  -1.32 / -2.59 / -4.76
## (placebo / oseltamivir / baloxavir).  The virological column is close; the
## symptom column is the model's weakest, and analysis A8 says why.
FLU_run_scenarios <- function() {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  sc  <- flu_scenarios(trx)
  do.call(rbind, lapply(names(sc), function(nm) {
    s   <- sc[[nm]]
    out <- do.call(sim, c(list(doses = s[[1]], tmax = 30), s[[2]]))
    r   <- summarise_run(out, s[[3]])
    data.frame(scenario = nm, TTAS_h = r$ttas_h, shed_h = r$shed_h,
               drop24 = r$drop24, AUC = r$auc_log, peak = r$peak_log,
               mut_peak = r$mut_peak, fever_h = r$fever_h,
               LRT = r$lrt_peak, epi_lost = r$epi_lost)
  }))
}

## --- A4 -------------------------------------------------------------------
## The two headrooms, at the median enrolment time:
##   bound R(t_rx) = 14.27 log10.d; the calibrated drug takes 9.37 (65.6%)
##   potency headroom (Emax -> 1)      +0.09 log10.d
##   timing  headroom (dose at onset)  +4.34 log10.d
## Fifty-fold more benefit sits on the axis a better molecule cannot buy.
FLU_timing_vs_potency <- function() {
  base  <- sim(NULL, tmax = 30); ons <- symptom_onset(base)
  R     <- residual_auc(base);   total <- R[1]
  timing <- do.call(rbind, lapply(c(0, 6, 12, 24, 36, 48, 72, 96), function(h) {
    trx <- ons + h / 24; i <- which.min(abs(base$time - trx))
    o   <- sim(rx_baloxavir(trx), tmax = 30)
    ach <- total - viral_auc(o); r <- summarise_run(o, trx)
    data.frame(h_post_onset = h, T_at_dose = base$TFRAC[i] * 100,
               bound = R[i], achieved = ach, left_on_table = R[i] - ach,
               TTAS_h = r$ttas_h, shed_h = r$shed_h)
  }))
  trx  <- ons + 1; i_rx <- which.min(abs(base$time - trx))
  potency <- do.call(rbind, lapply(
    c(0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.997, 0.9999, 0.999999), function(em) {
      o <- sim(rx_baloxavir(trx), tmax = 30, EMAX_BX = em)
      ach <- total - viral_auc(o); r <- summarise_run(o, trx)
      data.frame(Emax = em, achieved = ach, pct_of_bound = ach / R[i_rx] * 100,
                 TTAS_h = r$ttas_h, shed_h = r$shed_h)
    }))
  i_ons <- which.min(abs(base$time - ons))
  list(timing = timing, potency = potency,
       potency_headroom = tail(potency$achieved, 1) -
                          potency$achieved[potency$Emax == 0.9999],
       timing_headroom  = R[i_ons] - R[i_rx])
}

## --- A5 -------------------------------------------------------------------
## Each operator switched on alone at 95% from t_rx, no PK, so they can be
## compared on equal terms.  Reference values (TTAS h / shed h / d24 / AUC):
##   production-release  66.2 /  88.5 / -2.63 / 13.89
##   transcription       73.5 / 114.5 / -1.52 / 18.23
##   entry               78.2 / 117.5 / -1.34 / 18.92
##   virion clearance    65.8 /  88.0 / -2.70 / 13.24
##   infected-cell death 32.2 /  29.0 / -4.16 /  9.56   <- strongest by far
##   target protection   78.5 / 117.5 / -1.33 / 18.94
##   (placebo)           78.7 / 118.0 / -1.32 / 19.00
## Infected-cell death is the only operator that removes the SOURCE rather
## than throttling it -- and no licensed influenza antiviral works that way.
FLU_operator_decomposition <- function(E = 0.95) {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  p0  <- as.list(param(flu))
  ops <- list(
    "production/release (p)"  = list(PVIR = p0$PVIR * (1 - E),
                                     PVIRL = p0$PVIRL * (1 - E)),
    "transcription (E->I)"    = list(FMATU = 1 - E),
    "entry (beta)"            = list(BETA = p0$BETA * (1 - E),
                                     BETAL = p0$BETAL * (1 - E)),
    "virion clearance (c)"    = list(CVIR = p0$CVIR / (1 - E)),
    "infected-cell death"     = list(DELTA = p0$DELTA / (1 - E),
                                     DELTAL = p0$DELTAL / (1 - E)),
    "target protection (T->R)"= list(PHIF = p0$PHIF / (1 - E)),
    "(none: placebo)"         = list()
  )
  ## NOTE.  mrgsolve has no time-varying parameter switch in this formulation,
  ## so the operator is applied from t = 0 rather than from t_rx; the Python
  ## reference applies it at t_rx by restarting the integrator.  The ORDERING
  ## of the operators is the same; the magnitudes are not directly comparable
  ## between the two, and the reference implementation is the one quoted.
  do.call(rbind, lapply(names(ops), function(nm) {
    out <- do.call(sim, c(list(doses = NULL, tmax = 30), ops[[nm]]))
    r   <- summarise_run(out, trx)
    data.frame(operator = nm, TTAS_h = r$ttas_h, shed_h = r$shed_h,
               drop24 = r$drop24, AUC = r$auc_log, epi_lost = r$epi_lost)
  }))
}

## --- A6 -------------------------------------------------------------------
## Peak mutant titre over (dose time x potency).  Two results that were not
## put in: selection is non-monotone in potency, and it REQUIRES a target-cell
## field -- dosing late selects almost nothing, because the wild type has
## already eaten the epithelium and there is nothing left to hand over.
## Resistance emergence and clinical benefit share a prerequisite.
## DISCREPANCY 4: this is deterministic, so a lineage at 1e-4 is "present" in
## every simulated patient.  Read the table as selection pressure, not as a
## predicted viral load in any individual.
FLU_competitive_release <- function() {
  ons <- placebo_onset()
  hs  <- c(0, 12, 24, 36, 48, 72)
  ems <- c(0, 0.90, 0.99, 0.999, 0.9999, 0.999999)
  do.call(rbind, lapply(hs, function(h) {
    trx <- ons + h / 24
    row <- vapply(ems, function(em) {
      o <- do.call(sim, c(list(doses = rx_baloxavir(trx), tmax = 30,
                               EMAX_BX = em), I38T))
      max(o$LOGVM[o$time >= trx])
    }, numeric(1))
    setNames(data.frame(h, t(row)), c("h_post_onset", paste0("Emax_", ems)))
  }))
}

## --- A7 -------------------------------------------------------------------
## DISCREPANCY 5, and the model's clearest failure.  The ordering comes out
## BACKWARDS: the adult with pre-existing immunity selects the most mutant and
## the child the least, the reverse of the paediatric-vs-adult signal.  The
## mechanism is visible and is not a coding accident -- release is governed by
## T(t_rx), and anything that SLOWS the wild type leaves more of it, including
## prior immunity.  What the model does get right is the immunocompromised
## host, whose mutant AUC is the largest by a wide margin.
FLU_host_comparison <- function() {
  hosts <- list(
    "adult, prior exposure"   = list(CTL0 = 6.0, IGA0 = 0.30, FAB = 1.00, T0U = 4e8),
    "adult, naive"            = list(CTL0 = 1.0, IGA0 = 0.00, FAB = 1.00, T0U = 4e8),
    "child (naive, larger T)" = list(CTL0 = 1.0, IGA0 = 0.00, FAB = 0.65, T0U = 6e8),
    "immunocompromised"       = list(CTL0 = 1.0, IGA0 = 0.00, FAB = 0.15,
                                     FCTL = 0.12, FIFN = 0.55, T0U = 4e8)
  )
  do.call(rbind, lapply(names(hosts), function(nm) {
    ov  <- hosts[[nm]]
    ons <- do.call(placebo_onset, ov)
    trx <- ons + 0.5
    out <- do.call(sim, c(list(doses = rx_baloxavir(trx), tmax = 35),
                          modifyList(ov, I38T)))
    post <- out$time >= trx
    dt   <- out$time[2] - out$time[1]
    y    <- pmax(out$LOGVM[post] - LOD, 0)
    data.frame(host = nm,
               shed_h    = time_to_shed_stop(out, trx),
               mut_peak  = max(out$LOGVM[post]),
               mut_above_LOD_d = sum(out$LOGVM[post] > LOD) * dt,
               mut_AUC   = sum(y) * dt,
               T_at_dose = out$TFRAC[which.min(abs(out$time - trx))] * 100)
  }))
}

## --- A8 -------------------------------------------------------------------
## What the symptom endpoint requires.  Sweeping WVIR:
##   WVIR   placebo TTAS   oseltamivir benefit   baloxavir benefit
##   0.00        55.5              0.0 h               0.8 h
##   0.60        78.7             15.2 h              34.7 h
##   1.00        86.2             32.7 h              54.2 h
## CAPSTONE-1 measured 26.5 h for both.  At WVIR = 0 -- symptoms driven purely
## by the cytokine state -- the model gives ZERO benefit at any potency,
## because the cytokine peak is settled before the drug arrives.  The published
## symptom endpoint therefore cannot be explained by a cytokine-only symptom
## model; WVIR = 0.60 is the price of matching it, and it is a falsifiable
## claim about influenza symptoms rather than a fitting convenience.
## DISCREPANCY 3: even at the calibrated WVIR the model splits the two drugs,
## whereas CAPSTONE-1 found them indistinguishable despite a 2-log difference
## in day-2 titre.
FLU_symptom_requirement <- function(ws = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) {
  do.call(rbind, lapply(ws, function(w) {
    ons <- placebo_onset(WVIR = w); trx <- ons + RX_DELAY_H / 24
    pl  <- sim(NULL, tmax = 30, WVIR = w)
    os  <- sim(rx_oseltamivir(trx), tmax = 30, WVIR = w)
    bx  <- sim(rx_baloxavir(trx),   tmax = 30, WVIR = w)
    tp  <- time_to_alleviation(pl, trx)
    data.frame(WVIR = w, onset_h = ons * 24, placebo_TTAS = tp,
               oselt_TTAS = time_to_alleviation(os, trx),
               balox_TTAS = time_to_alleviation(bx, trx),
               oselt_benefit = tp - time_to_alleviation(os, trx),
               balox_benefit = tp - time_to_alleviation(bx, trx))
  }))
}

## --- A9 -------------------------------------------------------------------
## The direction is right and the size is small, and both matter.  The steroid
## lowers the symptom score and shortens TTAS while raising viral AUC and
## lengthening shedding -- the shape of the observational mortality signal
## (Lansbury 2019 Cochrane, PMID 30798570) from mechanism alone.  But in THIS
## calibration neither IFN-driven target protection nor CD8 killing is
## rate-limiting for clearance, so suppressing them costs little.  A host in
## which the CD8 response IS rate-limiting -- the elderly, the immunosuppressed,
## exactly where the observational signal comes from -- would show a much
## larger penalty.  Read the small effect as a property of the calibration.
FLU_steroid <- function() {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  arms <- list(
    "oseltamivir alone"        = rx_oseltamivir(trx),
    "oseltamivir + dexameth."  = rbind(rx_oseltamivir(trx), rx_dexamethasone(trx)),
    "dexamethasone alone"      = rx_dexamethasone(trx),
    "placebo"                  = NULL
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    out <- sim(arms[[nm]], tmax = 30); r <- summarise_run(out, trx)
    data.frame(arm = nm, AUC = r$auc_log, shed_h = r$shed_h,
               peak_sym = r$peak_sym, TTAS_h = r$ttas_h, bac_peak = r$bac_peak)
  }))
}

## --- A10 ------------------------------------------------------------------
## Prophylaxis is the same molecule on the other side of the peak.  Dosed at
## or before 12 h post-infection the model never registers an illness at all;
## dosed at 48 h it shortens one.  BLOCKSTONE (Ikematsu 2020 NEJM 383:309,
## PMID 32640124) reported 1.9% vs 13.6% symptomatic influenza with single-dose
## baloxavir prophylaxis -- a reduction in the PROBABILITY of illness, not in
## its duration.  The switch between the two regimes is the sign of
## (peak time - dose time).
FLU_prophylaxis <- function(hours = c(-12, 0, 12, 24, 36, 48, 72)) {
  do.call(rbind, lapply(hours, function(h) {
    trx <- max(0, h / 24)
    out <- sim(rx_baloxavir(trx), tmax = 25); r <- summarise_run(out, 0)
    data.frame(dose_h_pi = h, peak_log = r$peak_log,
               shed_d = sum(out$LOGV > LOD) * (out$time[2] - out$time[1]),
               peak_sym = r$peak_sym, epi_lost = r$epi_lost,
               illness = ifelse(r$peak_sym >= 2.0, "yes", "NO"))
  }))
}

## --- A11 ------------------------------------------------------------------
## Nothing in this table crosses the superinfection threshold, including the
## deliberately worsened arms.  The threshold has no independent calibration --
## it is a stipulation -- so the module can ORDER arms by bacterial risk and
## cannot say who develops pneumonia.  The one arm that changes the endpoint is
## the one that aborts the infection outright, which is the honest version of
## "antivirals prevent secondary pneumonia".
FLU_lrt_bacterial <- function() {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  arms <- list(
    "placebo"                     = list(NULL, list()),
    "oseltamivir 24 h post-onset" = list(rx_oseltamivir(trx), list()),
    "baloxavir 24 h post-onset"   = list(rx_baloxavir(trx), list()),
    "baloxavir at onset"          = list(rx_baloxavir(ons), list()),
    "baloxavir 72 h post-onset"   = list(rx_baloxavir(ons + 3), list()),
    "placebo, 5x LRT seeding"     = list(NULL, list(ADESC = 0.040)),
    "placebo, blunted interferon" = list(NULL, list(FIFN = 0.45))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    a   <- arms[[nm]]
    out <- do.call(sim, c(list(doses = a[[1]], tmax = 25), a[[2]]))
    data.frame(arm = nm, LRT_peak = max(out$LOGVL),
               LRT_damage_pct = max(out$DL) / as.list(param(flu))$T0L * 100,
               bac_peak = max(out$BAC), spo2_min = min(out$SPO2),
               crosses = ifelse(max(out$SUPER) > 0, "YES", "no"))
  }))
}

## --- A12 ------------------------------------------------------------------
## DISCREPANCY 2.  In the immunocompetent host the drug competes with an
## immune response that would have cleared the virus anyway, so its marginal
## value decays with time (benefit 9.37 -> 1.90 -> 0.00 log10.d at 24 / 96 /
## 168 h).  Remove that response and the decay slows (11.89 -> 4.95 -> 0.08).
## Consistent with treating these patients regardless of symptom duration, but
## no randomised trial has tested it: the prediction is unfalsified, not
## validated.
FLU_immunocompromised <- function(hours = c(0, 24, 48, 96, 168)) {
  hosts <- list("immunocompetent"   = list(),
                "immunocompromised" = list(FCTL = 0.12, FAB = 0.15, FIFN = 0.55))
  do.call(rbind, lapply(names(hosts), function(nm) {
    ov   <- hosts[[nm]]
    ons  <- do.call(placebo_onset, ov)
    base <- do.call(sim, c(list(doses = NULL, tmax = 35), ov))
    b    <- viral_auc(base)
    do.call(rbind, lapply(hours, function(h) {
      trx <- ons + h / 24
      out <- do.call(sim, c(list(doses = rx_baloxavir(trx), tmax = 35), ov))
      data.frame(host = nm, dose_h_post_onset = h,
                 shed_d = sum(out$LOGV > LOD) * (out$time[2] - out$time[1]),
                 AUC = viral_auc(out), benefit = b - viral_auc(out))
    }))
  }))
}

## --- A13 ------------------------------------------------------------------
## Local sensitivity, +20% on each parameter, relative sensitivity
## (dY/Y)/(dX/X).  Values above 1 in magnitude mean the output is amplifying
## the parameter's uncertainty.
FLU_sensitivity <- function(keys = c("BETA", "PVIR", "DELTA", "CVIR", "KECL",
                                     "QF", "PHIF", "KKILL", "RCTL", "T0U",
                                     "LREG", "EC50_BX", "COST", "MU", "K6",
                                     "KOFF", "WVIR")) {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  p0  <- as.list(param(flu))
  b   <- summarise_run(do.call(sim, c(list(doses = rx_baloxavir(trx),
                                           tmax = 30), I38T)), trx)
  rel <- function(new, old) if (is.na(new) || is.na(old) || old == 0) NA
                            else (new - old) / old / 0.2
  do.call(rbind, lapply(keys, function(k) {
    ov <- modifyList(I38T, setNames(list(p0[[k]] * 1.2), k))
    r  <- summarise_run(do.call(sim, c(list(doses = rx_baloxavir(trx),
                                            tmax = 30), ov)), trx)
    data.frame(parameter = k, sens_TTAS = rel(r$ttas_h, b$ttas_h),
               sens_AUC = rel(r$auc_log, b$auc_log),
               sens_mutfrac = rel(r$mut_frac, b$mut_frac))
  }))
}

## --- the ledger, in one table ---------------------------------------------
FLU_trial_ledger <- function() {
  ons <- placebo_onset(); trx <- ons + RX_DELAY_H / 24
  arms <- list(placebo = NULL, oseltamivir = rx_oseltamivir(trx),
               baloxavir = rx_baloxavir(trx))
  m <- lapply(arms, function(d) summarise_run(sim(d, tmax = 30), trx))
  data.frame(
    endpoint = c("TTAS (h)", "shedding cessation (h)", "day-2 titre (log10)"),
    published_placebo     = c(80.2, 96.0, -1.3),
    model_placebo         = c(m$placebo$ttas_h, m$placebo$shed_h, m$placebo$drop24),
    published_oseltamivir = c(53.8, 72.0, -2.8),
    model_oseltamivir     = c(m$oseltamivir$ttas_h, m$oseltamivir$shed_h,
                              m$oseltamivir$drop24),
    published_baloxavir   = c(53.7, 24.0, -4.8),
    model_baloxavir       = c(m$baloxavir$ttas_h, m$baloxavir$shed_h,
                              m$baloxavir$drop24)
  )
}

##############################################################################
## 5. WHAT THIS MODEL IS NOT
##############################################################################
## - It is deterministic and well-mixed.  Influenza replicates across a
##   spatially structured mucosa, and both stochastic extinction and spatial
##   segregation of strains matter for exactly the resistance question the
##   model is used to ask (see DISCREPANCY 4 and 5).
## - Symptoms are a one-compartment composite score with a calibrated coupling
##   to viral load.  The coupling (WVIR) is the model's least defensible
##   parameter and its most important, because every symptom endpoint in every
##   influenza trial depends on it.
## - The superinfection threshold is stipulated, not calibrated (A11).
## - The mAb, favipiravir and peramivir PD parameters are assumed rather than
##   fitted; only oseltamivir and baloxavir are calibrated against trial data.
## - Nothing here is validated for clinical use.
##############################################################################
