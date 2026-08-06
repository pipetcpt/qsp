# =============================================================================
# cgd_mrgsolve_model.R
# Chronic Granulomatous Disease (만성 육아종병) — QSP model, mrgsolve
#
# WHAT THIS FILE IS
# -----------------
# A mirror, equation for equation, of cgd_python_reference.py, which is the
# implementation that has actually been executed.  Where the two disagree, the
# Python file is right.  Numeric parameters are NOT hand-copied: they are
# injected into the $PARAM block below by
#
#     python3 sync_r_params.py
#
# from cgd_calibration.json, which cgd_analysis.py writes.  Everything between
# the AUTO-SYNC markers is machine-written.
#
# THE ONE IDEA
# ------------
# CGD is a disease of one number: the electron flux phi through NOX2.  The
# phagosomal chemistry that turns that flux into a dead bacterium lives in
# cgd_kernel.py and is NOT reproduced here — it is summarised by the single
# function K(phi), the log10 kill of one ingested organism in 60 minutes,
# fitted as a Hill curve.  That curve, plus the rule for averaging it over a
# MIXED neutrophil population, is the entire interface.
#
#   s(phi)  = 10^-K(phi)                     survival of one ingested organism
#   s_enc   = f*s(1) + (1-f)*s(r)            f  = fraction of oxidase-competent
#                                            r  = residual activity of the rest
#
# The averaging happens in SURVIVAL space, not in log-kill space.  Averaging
# log kills is the arithmetic error that makes an X-linked carrier and a
# hypomorphic hemizygote — who report the same DHR mean — look like the same
# patient.  They are not; see cgd_reference_output.txt, Section D.
#
# WHAT THE MODEL COMPUTES THAT MATTERS
# ------------------------------------
#   * a bacterial focus with a SHELTERED compartment (BACTi): organisms that
#     were ingested by a phagocyte which failed to kill them, are hidden from
#     every competent neutrophil in the lesion, divide, and are released alive.
#     This term is where the clinical correction threshold comes from — NOT
#     from K(phi), which has no threshold at all.
#   * a NECROSIS compartment that walls the focus off from both neutrophils and
#     antibiotic once it outruns recruitment.  Without it the outcome does not
#     depend on inoculum size and there is no critical inoculum to find.
#   * an inflammatory arm driven by UNCLEARED APOPTOTIC NEUTROPHILS, which
#     produces granuloma and colitis in a patient who has never been infected
#     with anything — and which antimicrobials therefore cannot touch.
#
# USAGE
#   library(mrgsolve); library(dplyr)
#   mod <- mread("cgd_mrgsolve_model", ".", quiet = TRUE)
#   # X-CGD null on triple prophylaxis for a year
#   mod %>% param(f_normal = 0, phi_res = 0) %>%
#     ev(cgd_regimen(tmpsmx = TRUE, itra = TRUE, ifng = TRUE, days = 365)) %>%
#     mrgsim(end = 365, delta = 1) %>% plot()
#
# NOTE ON EXECUTION: no R toolchain was available in the environment where this
# model was built, so this file has not been run.  Its equations mirror a
# Python implementation that has been run extensively; its R syntax has not
# been exercised by an interpreter.
# =============================================================================

$PROB
# Chronic Granulomatous Disease QSP model — 53 ODEs
# Layer 1 (phagosome chemistry) enters only as the Hill parameters K0/Kmax/
# phi50/hill.  Layers 2 and 3 are here in full.

$PARAM @annotated
// ---- BEGIN AUTO-SYNC (written by sync_r_params.py) ------------------------
K0      : 0.22185 : log10 kill at phi=0, from the phagosome kernel
Kmax    : 1.30103 : log10 kill at phi=1, from the phagosome kernel
phi50   : 0.54482 : Hill midpoint of K(phi)
hillK   : 2.73506 : Hill coefficient of K(phi)
// ---- END AUTO-SYNC ---------------------------------------------------------

// ---- patient descriptors ---------------------------------------------------
f_normal : 0.0  : fraction of neutrophils with a fully competent oxidase
phi_res  : 0.0  : residual oxidase activity of the remaining neutrophils
chim_tgt : 0.0  : donor myeloid chimerism target (0-1) after HSCT
vcn_tgt  : 0.0  : gene-corrected myeloid fraction target (0-1)
WT       : 30   : kg
extra_hz : 0.0  : additional per-year hazard (e.g. transplant TRM)
ido_nox  : 0    : 1 = IDO requires NADPH oxidase (Romani 2008; did not replicate)

// ---- granulopoiesis and neutrophil traffic ---------------------------------
tau_prog   : 6.0    : d, marrow mitotic pool transit
tau_res    : 5.0    : d, marrow post-mitotic reserve transit
kout_ANC   : 2.3765 : /d, log(2)/(7 h)
ANC0       : 4000   : cells/uL
GCSF0      : 0.025  : ng/mL
kin_GCSF   : 1.0    : /d
kout_GCSF  : 24.0   : /d
GCSF_EC50  : 0.20   : ng/mL
GCSF_Emax  : 6.0    : -
kcl_GCSF_N : 0.004  : /d per cell/uL, receptor-mediated G-CSF clearance
MONO0      : 400    : cells/uL
kout_MONO  : 1.0    : /d
k_rec      : 40.0   : /d, maximal neutrophil recruitment
Krec       : 50.0   : CFU, half-maximal recruitment
kout_NEUT  : 1.2    : /d
k_rec_mono : 4.0    : /d
kout_MAC   : 0.10   : /d

// ---- bacterial focus --------------------------------------------------------
mu_b       : 4.0    : /d, S. aureus net growth in tissue
Bmax       : 1e10   : CFU
k_up       : 0.03   : /d per cell/uL, phagocytic uptake
k_release  : 0.50   : /d, release when the host phagocyte dies
mu_bi      : 0.10   : /d, replication inside a failed phagosome
Km_up      : 1e7    : CFU, uptake saturation
N_resident : 200    : cells/uL-equivalent, resident phagocytes
k_nec      : 3.0    : /d
Knec       : 1e5    : CFU
kout_nec   : 0.15   : /d

// ---- fungal focus ------------------------------------------------------------
k_germ     : 0.35   : /d, conidium -> hypha
mu_h       : 2.2    : /d, hyphal extension
Hmax       : 1e9    : hyphal units
k_up_c     : 0.25   : /d per cell/uL, macrophage conidial uptake
k_dam_h    : 0.03   : /d per cell/uL, neutrophil hyphal damage
k_GM       : 3e-8   : galactomannan index per hypha per day
kout_GM    : 1.4    : /d

// ---- cytokines ---------------------------------------------------------------
IL1B0 : 2.0 : pg/mL
IL180 : 60  : pg/mL
TNF0  : 5.0 : pg/mL
IL60  : 2.0 : pg/mL
IL170 : 5.0 : pg/mL
IL100 : 4.0 : pg/mL
IFNg0 : 3.0 : pg/mL
CRP0  : 1.0 : mg/L
kout_IL1B : 24 : /d
kout_IL18 : 12 : /d
kout_TNF  : 40 : /d
kout_IL6  : 24 : /d
kout_IL17 : 6  : /d
kout_IL10 : 24 : /d
kout_IFNg : 30 : /d
kout_CRP  : 0.37 : /d
k_CRP_IL6 : 0.185 : mg/L per pg/mL per d
k_IL1_stim : 0.030 : -
k_IL1_apop : 0.60  : -
k_TNF_stim : 0.020 : -
k_IL6_TNF  : 0.35  : -
k_IL6_IL1  : 0.55  : -
k_IL17_IL1 : 0.020 : -
k_IL10_eff : 0.40  : -
k_IFNg_IL18 : 0.010 : -
amp_IL1 : 4.0 : fold amplification of IL-1b per unit lost oxidase

// ---- resolution failure --------------------------------------------------------
k_apop     : 0.55 : /d
k_effero   : 3.0  : /d at full oxidase competence
effero_min : 0.25 : residual efferocytosis at phi = 0
N_tissue0  : 50   : resident tissue phagocytes turning over everywhere
gut_stim   : 1.0  : constant translocated microbial stimulus (everyone has it)
k_NET      : 0.30 : /d
kout_NET   : 1.0  : /d

// ---- granuloma / colitis --------------------------------------------------------
k_gran        : 0.135 : granuloma formation gain (normalisation)
kout_GRAN     : 0.020 : /d
IL1_gran_EC50 : 30    : pg/mL above baseline
k_fib         : 0.010 : /d
kout_FIB      : 0.0015 : /d
k_COL         : 0.286 : colitis gain (normalisation)
kout_COL      : 0.10  : /d

// ---- tryptophan / kynurenine ------------------------------------------------------
TRP0 : 55 : uM
kin_TRP : 55 : uM/d
kout_TRP : 1.0 : /d
k_IDO : 0.35 : /d
IDO_IFNg_EC50 : 20 : pg/mL
kout_KYN : 3.0 : /d

// ---- pharmacokinetics.  amounts mg, volumes L, so concentrations mg/L. -------------
// Interferon gamma-1b alone is ug and L, so its concentration is ug/L = ng/mL.
TMP_F : 0.90 : -
TMP_ka : 12 : /d
TMP_Vc : 30 : L
TMP_Vp : 18 : L
TMP_Q  : 30 : L/d
TMP_CL : 82 : L/d
TMP_fu : 0.55 : -
SMX_F : 0.85 : -
SMX_ka : 10 : /d
SMX_V : 6.3 : L
SMX_CL : 10.5 : L/d
ITZ_F : 0.55 : -
ITZ_ka : 4 : /d
ITZ_Vc : 120 : L
ITZ_Vp : 250 : L
ITZ_Q : 150 : L/d
ITZ_CL : 220 : L/d
OHITZ_fm : 0.75 : -
OHITZ_V : 300 : L
OHITZ_CL : 150 : L/d
VOR_F : 0.95 : -
VOR_ka : 12 : /d
VOR_V : 138 : L
VOR_Vmax : 1400 : mg/d
VOR_Km : 3.0 : mg/L
IFN_F : 0.90 : -
IFN_ka : 6.0 : /d
IFN_V : 30 : L
IFN_CL : 100 : L/d
PRED_F : 0.85 : -
PRED_ka : 20 : /d
PRED_V : 18 : L
PRED_CL : 100 : L/d
PRED_fu : 0.25 : -
ANA_F : 0.95 : -
ANA_ka : 6.0 : /d
ANA_V : 18 : L
ANA_CL : 60 : L/d

// ---- pharmacodynamics ---------------------------------------------------------------
TMP_Emax : 6.0 : /d
TMP_EC50 : 0.50 : mg/L free trimethoprim
ITZ_Emax : 2.4 : /d
ITZ_EC50 : 0.25 : mg/L
VOR_Emax : 3.4 : /d
VOR_EC50 : 0.50 : mg/L
IFN_Emax_nonox : 1.9 : fold boost of the OXIDASE-INDEPENDENT arm
IFN_EC50 : 0.30 : ng/mL
IFN_Emax_rec : 0.7 : -
PRED_Emax_cyt : 0.85 : -
PRED_EC50 : 12 : ng/mL free prednisolone
PRED_Emax_mig : 0.45 : -
ANA_Imax : 0.92 : -
ANA_IC50 : 250 : ng/mL

// ---- engraftment and mortality -------------------------------------------------------
k_engraft : 0.12 : /d
k_vcn_loss : 0.0011 : /d
k_h_bact : 2e-5 : /d per log10 CFU above 4
k_h_fung : 1e-4 : /d per log10 hyphae above 3
h_base : 0.0008 : /yr
kout_ALT : 0.35 : /d
ALT0 : 25 : U/L

$CMT @annotated
PROG    : marrow mitotic pool (cells/uL/d equivalent)
RES     : marrow post-mitotic reserve
ANC     : blood neutrophils (cells/uL)
GCSF    : G-CSF (ng/mL)
MONO    : blood monocytes (cells/uL)
NEUT_T  : neutrophils at the focus
MAC     : lesion macrophages
EPI     : epithelioid / giant cells
BACT    : extracellular bacteria at the focus (CFU)
BACTi   : bacteria sheltering inside failed phagocytes (CFU)
NEC     : necrotic / walled-off tissue
CONID   : Aspergillus conidia
HYPH    : Aspergillus hyphae
GM      : galactomannan index
IL1B    : IL-1 beta (pg/mL)
IL18    : IL-18 (pg/mL)
TNF     : TNF-alpha (pg/mL)
IL6     : IL-6 (pg/mL)
IL17    : IL-17 (pg/mL)
IL10    : IL-10 (pg/mL)
IFNg    : endogenous IFN-gamma (pg/mL)
CRP     : C-reactive protein (mg/L)
APOP    : uncleared apoptotic neutrophils
NET     : NET material
GRAN    : granuloma burden
FIB     : fibrosis / stricture index
COL     : colitis activity index (0-10)
TRP     : tryptophan (uM)
KYN     : kynurenine (uM)
TMP_g   : trimethoprim gut (mg)
TMP_c   : trimethoprim central (mg)
TMP_p   : trimethoprim peripheral (mg)
SMX_g   : sulfamethoxazole gut (mg)
SMX_c   : sulfamethoxazole central (mg)
ITZ_g   : itraconazole gut (mg)
ITZ_c   : itraconazole central (mg)
ITZ_p   : itraconazole peripheral (mg)
OHITZ   : hydroxy-itraconazole (mg)
VOR_g   : voriconazole gut (mg)
VOR_c   : voriconazole central (mg)
IFN_sc  : interferon gamma-1b SC depot (ug)
IFN_c   : interferon gamma-1b central (ug)
PRED_g  : prednisolone gut (mg)
PRED_c  : prednisolone central (mg)
ANA_sc  : anakinra SC depot (mg)
ANA_c   : anakinra central (mg)
CHIM    : donor myeloid chimerism (0-1)
VCN     : gene-corrected myeloid fraction (0-1)
CUMBACT : cumulative bacterial burden-days
CUMFUNG : cumulative fungal burden-days
SURV    : survival probability
STEROID : cumulative prednisolone exposure (mg.d)
ALT     : alanine aminotransferase (U/L)

$GLOBAL
#define SAT(x, k) ((x) / ((x) + (k)))
// Below one organism there is nothing left to divide.  Without this factor an
// integrator round-off of 1e-30 grows to 1e10 CFU over a year at mu_b = 4/d,
// and the model reports an abscess in a patient who was never inoculated.
#define ALIVE(x) ((x) / ((x) + 1.0))
#define POS(x) ((x) > 0.0 ? (x) : 0.0)

$MAIN
double rel0 = 1.0 + GCSF_Emax * GCSF0 / (GCSF_EC50 + GCSF0);
double prod0 = ANC0 * kout_ANC;

PROG_0 = prod0 * tau_prog;
RES_0  = prod0 * tau_res;
ANC_0  = ANC0;
GCSF_0 = GCSF0;
MONO_0 = MONO0;
IL1B_0 = IL1B0;
IL18_0 = IL180;
TNF_0  = TNF0;
IL6_0  = IL60;
IL17_0 = IL170;
IL10_0 = IL100;
IFNg_0 = IFNg0;
CRP_0  = CRP0;
TRP_0  = TRP0;
SURV_0 = 1.0;
ALT_0  = ALT0;

$ODE
// ---------------------------------------------------------------- drug levels
double CTMP  = POS(TMP_c) / TMP_Vc;
double CTMPf = CTMP * TMP_fu;
double CTMPp = POS(TMP_p) / TMP_Vp;
double CSMX  = POS(SMX_c) / SMX_V;
double CITZ  = POS(ITZ_c) / ITZ_Vc;
double CITZp = POS(ITZ_p) / ITZ_Vp;
double COHI  = POS(OHITZ) / OHITZ_V;
double CVOR  = POS(VOR_c) / VOR_V;
double CIFN  = POS(IFN_c) / IFN_V;                       // ng/mL
double CPREDf = POS(PRED_c) / PRED_V * 1000.0 * PRED_fu; // ng/mL free
double CANA  = POS(ANA_c) / ANA_V * 1000.0;              // ng/mL

double E_TMP  = TMP_Emax * CTMPf / (TMP_EC50 + CTMPf);
double azol   = CITZ + 0.6 * COHI;
double E_AZOL = ITZ_Emax * azol / (ITZ_EC50 + azol)
              + VOR_Emax * CVOR / (VOR_EC50 + CVOR);
double E_IFNn = 1.0 + IFN_Emax_nonox * CIFN / (IFN_EC50 + CIFN);
double E_IFNr = 1.0 + IFN_Emax_rec   * CIFN / (IFN_EC50 + CIFN);
double I_PRED = PRED_Emax_cyt * CPREDf / (PRED_EC50 + CPREDf);
double I_PMIG = PRED_Emax_mig * CPREDf / (PRED_EC50 + CPREDf);
double I_ANA  = ANA_Imax * CANA / (ANA_IC50 + CANA);

// ------------------------------------------------- the oxidase state of the marrow
// Three sources of oxidase-competent neutrophils, and as far as the phagosome
// is concerned they are the same variable: what you were born with, what a
// donor gave you, and what a vector put back.
double f_norm = f_normal + POS(CHIM) + POS(VCN);
if (f_norm > 1.0) f_norm = 1.0;

// K(phi) from the phagosome kernel, and the interferon boost.  K0 is by
// construction the whole NON-OXIDATIVE log-kill (at phi = 0 no HOCl exists),
// so boosting that arm by E_IFNn adds (E_IFNn - 1) * K0 log to EVERY cell,
// competent or not.  Additive, and identical in both subpopulations: that is
// precisely the claim the 1991 interferon trial forces, since it found no
// reproducible restoration of superoxide.
double extra = (E_IFNn - 1.0) * K0;
double K_one = K0 + (Kmax - K0) * pow(1.0, hillK)
             / (pow(phi50, hillK) + pow(1.0, hillK));
double pr    = POS(phi_res);
double K_res = K0 + (Kmax - K0) * pow(pr, hillK)
             / (pow(phi50, hillK) + pow(pr, hillK) + 1e-30);
double s1 = pow(10.0, -(K_one + extra));
double sr = pow(10.0, -(K_res + extra));

// THE MIXTURE.  Averaged in SURVIVAL space, never in log-kill space.
double s_enc = f_norm * s1 + (1.0 - f_norm) * sr;
if (s_enc < 1e-12) s_enc = 1e-12;
if (s_enc > 1.0)   s_enc = 1.0;

// phi_eff is the population MEAN, and it is used only by the arms that
// genuinely answer to a mean — efferocytosis and inflammasome restraint —
// as opposed to killing, which answers to the distribution.
double phi_eff = f_norm + (1.0 - f_norm) * pr;

// ----------------------------------------------------------- granulopoiesis
double relx = 1.0 + GCSF_Emax * POS(GCSF) / (GCSF_EC50 + POS(GCSF));
double rel0b = 1.0 + GCSF_Emax * GCSF0 / (GCSF_EC50 + GCSF0);
double burden = POS(BACT) + POS(BACTi) + POS(HYPH);
double recruit = k_rec * E_IFNr * (1.0 - I_PMIG) * SAT(burden, Krec)
               * POS(ANC) / ANC0 * 20.0;

dxdt_PROG = ANC0 * kout_ANC - PROG / tau_prog;
dxdt_RES  = PROG / tau_prog - RES / tau_res * relx / rel0b;
dxdt_ANC  = RES / tau_res * relx / rel0b - kout_ANC * ANC - recruit;
dxdt_GCSF = kin_GCSF * (1.0 + 3.0 * SAT(POS(IL1B) + POS(TNF), 40.0))
          - kout_GCSF * GCSF - kcl_GCSF_N * GCSF * POS(ANC);
dxdt_MONO = MONO0 * kout_MONO * (1.0 + 1.5 * SAT(POS(IL6), 50.0))
          - kout_MONO * MONO;

// --------------------------------------------------------------- the focus
double pen = 1.0 / (1.0 + POS(NEC));
double NT  = (POS(NEUT_T) + N_resident) * pen;
double uptake = k_up * NT * POS(BACT) / (1.0 + POS(BACT) / Km_up);

dxdt_BACT  = mu_b * POS(BACT) * ALIVE(POS(BACT)) * (1.0 - POS(BACT) / Bmax)
           - uptake + k_release * POS(BACTi) - E_TMP * pen * POS(BACT);
dxdt_BACTi = uptake * s_enc + mu_bi * POS(BACTi) * ALIVE(POS(BACTi))
           - k_release * POS(BACTi);
dxdt_NEC   = k_nec * SAT(burden, Knec) - kout_nec * NEC;

// Conidia are ingestible and see s_enc directly; hyphae are far too large to
// swallow and can only be attacked from outside.  That asymmetry, and not any
// fitted constant, is why Aspergillus is the leading cause of death in CGD.
double hyph_comp = (f_norm + (1.0 - f_norm) * (0.10 + 0.90 * pr)) * E_IFNn;
double ingest = k_up_c * POS(MAC) * pen;
dxdt_CONID = -(ingest + k_germ) * POS(CONID);
dxdt_HYPH  = (ingest * s_enc + k_germ) * POS(CONID)
           + mu_h * POS(HYPH) * ALIVE(POS(HYPH)) * (1.0 - POS(HYPH) / Hmax)
           - k_dam_h * NT * POS(HYPH) * hyph_comp - E_AZOL * pen * POS(HYPH);
dxdt_GM    = k_GM * POS(HYPH) - kout_GM * GM;

dxdt_NEUT_T = recruit - kout_NEUT * NEUT_T;
dxdt_MAC = k_rec_mono * SAT(POS(CONID) + POS(HYPH) + POS(BACT), Krec)
         * POS(MONO) / MONO0 * 20.0 - kout_MAC * MAC;

// ---------------------------------------------- efferocytosis: the resolution defect
// Phosphatidylserine externalisation on the dying neutrophil is itself an
// oxidation event, so the CGD neutrophil dies without being recognised.  The
// uncleared corpse is BOTH an IL-1beta stimulus AND the missing source of the
// IL-10 that would have switched the lesion off.
double eff = effero_min + (1.0 - effero_min) * phi_eff;
dxdt_APOP = k_apop * (POS(NEUT_T) + POS(MAC) + N_tissue0) - k_effero * eff * POS(APOP);
dxdt_NET  = k_NET * POS(NEUT_T) * (0.25 + 0.75 * phi_eff) - kout_NET * NET;

// ------------------------------------------------------------------ cytokines
double amp = 1.0 + amp_IL1 * (1.0 - phi_eff);
double stim = log10(1.0 + burden + POS(CONID)) + gut_stim;
dxdt_IL1B = IL1B0 * kout_IL1B
          + amp * (k_IL1_stim * stim * 100.0 + k_IL1_apop * POS(APOP)) * (1.0 - I_PRED)
          - kout_IL1B * IL1B;
dxdt_IL18 = IL180 * kout_IL18
          + amp * k_IL1_stim * stim * 120.0 * (1.0 - I_PRED) - kout_IL18 * IL18;
dxdt_TNF  = TNF0 * kout_TNF + k_TNF_stim * stim * 100.0 * (1.0 - I_PRED)
          - kout_TNF * TNF;
double IL1free = POS(IL1B) * (1.0 - I_ANA);
dxdt_IL6  = IL60 * kout_IL6
          + (k_IL6_TNF * POS(TNF) + k_IL6_IL1 * IL1free) * (1.0 - I_PRED)
          - kout_IL6 * IL6;
dxdt_IL17 = IL170 * kout_IL17 + k_IL17_IL1 * IL1free * 100.0 * (1.0 - I_PRED)
          - kout_IL17 * IL17;
dxdt_IL10 = IL100 * kout_IL10 + k_IL10_eff * k_effero * eff * POS(APOP)
          - kout_IL10 * IL10;
dxdt_IFNg = IFNg0 * kout_IFNg + k_IFNg_IL18 * POS(IL18) * 10.0 - kout_IFNg * IFNg;
dxdt_CRP  = k_CRP_IL6 * POS(IL6) - kout_CRP * CRP;

// ------------------------------------------------------- granuloma and colitis
// These answer to the INCREMENT above the healthy baseline, not to the absolute
// level.  Driving them from the absolute level gives a healthy control a
// granuloma, which an earlier draft duly reported.
double dIL1 = POS(IL1free - IL1B0);
double dIL17x = POS(POS(IL17) - IL170);
dxdt_EPI  = 0.35 * SAT(dIL1 + dIL17x, 60.0) * POS(MAC) - 0.02 * EPI;
dxdt_GRAN = k_gran * SAT(dIL1, IL1_gran_EC50) * (1.0 + 0.01 * POS(EPI))
          - kout_GRAN * GRAN * (1.0 + 2.0 * I_PRED + 2.0 * I_ANA);
dxdt_FIB  = k_fib * POS(GRAN) - kout_FIB * FIB;
dxdt_COL  = k_COL * (SAT(dIL1, 25.0) * 10.0 + 0.02 * POS(GRAN))
          - kout_COL * COL * (1.0 + 3.0 * I_PRED + 3.5 * I_ANA);

// ------------------------------------------------ tryptophan (switchable arm)
double ido = k_IDO * SAT(POS(IFNg), IDO_IFNg_EC50) * (ido_nox > 0.5 ? phi_eff : 1.0);
dxdt_TRP = kin_TRP - kout_TRP * TRP - ido * POS(TRP);
dxdt_KYN = ido * POS(TRP) - kout_KYN * KYN;

// -------------------------------------------------------------------------- PK
dxdt_TMP_g = -TMP_ka * TMP_g;
dxdt_TMP_c = TMP_F * TMP_ka * TMP_g - TMP_CL * CTMP - TMP_Q * (CTMP - CTMPp);
dxdt_TMP_p = TMP_Q * (CTMP - CTMPp);
dxdt_SMX_g = -SMX_ka * SMX_g;
dxdt_SMX_c = SMX_F * SMX_ka * SMX_g - SMX_CL * CSMX;
dxdt_ITZ_g = -ITZ_ka * ITZ_g;
dxdt_ITZ_c = ITZ_F * ITZ_ka * ITZ_g - ITZ_CL * CITZ - ITZ_Q * (CITZ - CITZp);
dxdt_ITZ_p = ITZ_Q * (CITZ - CITZp);
dxdt_OHITZ = OHITZ_fm * ITZ_CL * CITZ - OHITZ_CL * COHI;
dxdt_VOR_g = -VOR_ka * VOR_g;
dxdt_VOR_c = VOR_F * VOR_ka * VOR_g - VOR_Vmax * CVOR / (VOR_Km + CVOR);
dxdt_IFN_sc = -IFN_ka * IFN_sc;
dxdt_IFN_c  = IFN_F * IFN_ka * IFN_sc - IFN_CL * CIFN;
dxdt_PRED_g = -PRED_ka * PRED_g;
dxdt_PRED_c = PRED_F * PRED_ka * PRED_g - PRED_CL * POS(PRED_c) / PRED_V;
dxdt_ANA_sc = -ANA_ka * ANA_sc;
dxdt_ANA_c  = ANA_F * ANA_ka * ANA_sc - ANA_CL * POS(ANA_c) / ANA_V;

// ------------------------------------------------------ corrected myelopoiesis
dxdt_CHIM = k_engraft * (chim_tgt - CHIM);
dxdt_VCN  = k_engraft * (vcn_tgt - VCN) - k_vcn_loss * VCN;

// -------------------------------------------------------------------- outcomes
double lb = log10(1.0 + POS(BACT) + POS(BACTi));
double lh = log10(1.0 + POS(HYPH));
double haz = h_base / 365.0 + k_h_bact * POS(lb - 4.0) + k_h_fung * POS(lh - 3.0)
           + 0.00002 * POS(FIB) + extra_hz / 365.0;
dxdt_CUMBACT = POS(lb - 4.0);
dxdt_CUMFUNG = POS(lh - 3.0);
dxdt_SURV    = -haz * SURV;
dxdt_STEROID = POS(PRED_c);
double tox = 0.30 * (CITZ + COHI) + 0.25 * CVOR + 0.002 * (CTMP + CSMX);
dxdt_ALT = kout_ALT * ALT0 * (1.0 + tox) - kout_ALT * ALT;

$TABLE
double DHRmean = f_normal + POS(CHIM) + POS(VCN)
               + (1.0 - (f_normal + POS(CHIM) + POS(VCN))) * POS(phi_res);
if (DHRmean > 1.0) DHRmean = 1.0;
double logB = log10(1.0 + POS(BACT) + POS(BACTi));
double logH = log10(1.0 + POS(HYPH));
double CTMPo = POS(TMP_c) / TMP_Vc;
double CITZo = POS(ITZ_c) / ITZ_Vc;
double COHIo = POS(OHITZ) / OHITZ_V;
double CVORo = POS(VOR_c) / VOR_V;
double CIFNo = POS(IFN_c) / IFN_V;
double CANAo = POS(ANA_c) / ANA_V * 1000.0;
double CPREDo = POS(PRED_c) / PRED_V * 1000.0;

$CAPTURE @annotated
DHRmean : DHR-123 mean, i.e. the single number the flow lab reports
logB    : log10 total bacterial burden at the focus
logH    : log10 hyphal burden
CTMPo   : trimethoprim (mg/L)
CITZo   : itraconazole (mg/L)
COHIo   : hydroxy-itraconazole (mg/L)
CVORo   : voriconazole (mg/L)
CIFNo   : interferon gamma-1b (ng/mL)
CANAo   : anakinra (ng/mL)
CPREDo  : prednisolone (ng/mL)

# =============================================================================
# HELPERS — source this file's companion below, or paste into your session.
# =============================================================================
# cgd_regimen <- function(days = 365, wt = 30, bsa = 1.1,
#                         tmpsmx = FALSE, itra = FALSE, vori = FALSE,
#                         ifng = FALSE, pred = FALSE, pred_dose = 1.0,
#                         anakinra = FALSE) {
#   ev_list <- list()
#   add <- function(cmt, amt, ii, n, start = 0) {
#     ev_list[[length(ev_list) + 1]] <<-
#       mrgsolve::ev(time = start, cmt = cmt, amt = amt, ii = ii, addl = n - 1)
#   }
#   if (tmpsmx) { add("TMP_g", 5 * wt / 2, 0.5, days * 2)
#                 add("SMX_g", 25 * wt / 2, 0.5, days * 2) }
#   if (itra)   add("ITZ_g", if (wt < 40) 100 else 200, 0.5, days * 2)
#   if (vori)   add("VOR_g", 9 * wt, 0.5, days * 2)
#   if (ifng) {                       # 50 ug/m2 SC on days 0, 2, 4 of each week
#     for (w in 0:floor(days / 7)) for (off in c(0, 2, 4)) {
#       tt <- w * 7 + off
#       if (tt <= days) add("IFN_sc", 50 * bsa, 1e9, 1, start = tt)
#     }
#   }
#   if (pred)     add("PRED_g", pred_dose * wt, 1, days)
#   if (anakinra) add("ANA_sc", 2 * wt, 1, days)
#   Reduce(function(a, b) a + b, ev_list)
# }
#
# GENOTYPES <- data.frame(
#   name    = c("healthy", "X-CGD null", "X-CGD hypomorph 5%",
#               "X-CGD hypomorph 20%", "p47phox AR", "p67phox AR",
#               "X-CGD carrier 50%", "X-CGD carrier 20%", "X-CGD carrier 5%"),
#   f_normal = c(1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.20, 0.05),
#   phi_res  = c(1.00, 0.00, 0.05, 0.20, 0.12, 0.03, 0.00, 0.00, 0.00))
#
# NOTE that rows 4 and 8 have the SAME DHR mean (0.20) and are different
# patients.  That is the point of the model.
