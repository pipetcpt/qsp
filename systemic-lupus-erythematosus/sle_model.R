################################################################################
# Systemic Lupus Erythematosus (SLE) — Quantitative Systems Pharmacology Model
# mrgsolve-compatible ODE model with PK/PD for:
#   - Hydroxychloroquine (HCQ)  — anti-malarial / TLR7/9 blocker
#   - Belimumab                  — anti-BAFF/BLyS mAb (TMDD)
#   - Anifrolumab                — anti-IFNAR1 mAb (type I IFN blocker)
#   - Mycophenolate Mofetil (MMF/MPA) — IMPDH inhibitor
#   - Voclosporin / Tacrolimus   — calcineurin inhibitors (LN)
#   - Corticosteroids            — broad immunosuppression
#
# Reference: Ding et al. CPT PSP 2020; Forde et al. J PK/PD 2021;
#            Navarra et al. Lancet 2011; Morand et al. NEJM 2020
#
# Disease PD state variables:
#   IFN pathway: TLR_act, pDC_IFN, IFNa, IFNscore
#   B cell axis: BAFF_free, Bcell, Plasmablast, LLPC_cells
#   Autoantibodies: AntiDsDNA
#   Complement: C3_serum, C4_serum
#   Immune complexes: IC_burden
#   Kidney: Proteinuria, eGFR
#   Clinical: SLEDAI_score (calculated)
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ==============================================================================
# MODEL CODE
# ==============================================================================

sle_model_code <- '
$PROB
 SLE QSP Model v1.0
 Multi-drug PK/PD integrating IFN-BAFF-B cell-autoAb-complement-kidney axis

$PARAM @annotated
// ---- HCQ PK ----
Ka_HCQ    : 0.020  : /h  HCQ absorption rate constant (t_lag~2h corrected)
F_HCQ     : 0.40   : -   HCQ oral bioavailability (~40%)
Vd_HCQ    : 5200   : L   HCQ apparent volume of distribution (large tissue Vd)
CL_HCQ    : 100    : L/h HCQ total clearance (t1/2 ~40-50 days)
dose_HCQ  : 400    : mg  HCQ daily oral dose (400 mg/day; 200 mg BID)

// ---- Belimumab PK (2-compartment IV, TMDD) ----
CL_BELI   : 0.215  : L/day  Belimumab clearance
Q_BELI    : 0.45   : L/day  Belimumab inter-compartmental clearance
V1_BELI   : 5.29   : L      Belimumab central volume
V2_BELI   : 3.46   : L      Belimumab peripheral volume
// TMDD parameters for BAFF target
Rin_BAFF  : 0.42   : ng/mL/day  BAFF synthesis rate at baseline
kout_BAFF : 0.15   : /day       BAFF natural elimination rate constant
kon_BELI  : 0.60   : 1/(ng/mL*day)  Belimumab-BAFF association rate
koff_BELI : 0.002  : /day       Belimumab-BAFF dissociation rate
kint_BELI : 0.20   : /day       Internalization rate of Belimumab-BAFF complex
dose_BELI : 0      : mg        Belimumab IV dose (0 = off; 10 mg/kg = ~700 mg q4w)

// ---- Anifrolumab PK (1-compartment IV) ----
CL_ANIF   : 0.170  : L/day  Anifrolumab clearance
V_ANIF    : 6.8    : L      Anifrolumab central volume
t50_ANIF  : 0.012  : -      Target occupancy EC50 as fraction of dose-normalized conc
dose_ANIF : 0      : mg     Anifrolumab IV dose (0=off; 300 mg q4w)

// ---- MMF/MPA PK ----
Ka_MMF    : 2.5    : /h   MMF absorption rate (rapid hydrolysis to MPA)
F_MMF     : 0.93   : -    MMF/MPA oral bioavailability
Vd_MPA    : 3.6    : L/kg MPA apparent volume (use 70kg → 252 L)
CL_MPA    : 10.5   : L/h  MPA total clearance (t1/2 ~17 h)
dose_MMF  : 0      : mg   MMF oral dose per administration (0=off; 1500 mg BID=3000 mg/day)

// ---- Voclosporin/Calcineurin Inhibitor PK ----
Ka_VOC    : 1.2    : /h   Voclosporin absorption
F_VOC     : 0.29   : -    Voclosporin bioavailability
Vd_VOC    : 2300   : L    Voclosporin Vd (highly lipophilic)
CL_VOC    : 63     : L/h  Voclosporin clearance
IC50_VOC  : 0.8    : ng/mL Voclosporin IC50 calcineurin inhibition
Imax_VOC  : 0.95   : -     Voclosporin maximum calcineurin inhibition
dose_VOC  : 0      : mg    Voclosporin dose (0=off; 23.7 mg BID)

// ---- Corticosteroid PD (simplified Emax) ----
CS_dose   : 0      : mg/day  Prednisone equivalent daily dose
EC50_CS   : 5.0    : mg/day  CS half-maximal anti-inflammatory effect
Emax_CS   : 0.80   : -       CS maximum suppression of inflammation

// ---- IFN Pathway PD ----
kin_TLR   : 0.10   : /day  TLR7/9 baseline activation signal
kout_TLR  : 0.10   : /day  TLR7/9 signal decay
kin_pDC   : 0.15   : /day  pDC IFN production baseline
kout_pDC  : 0.12   : /day  pDC IFN decay rate
IFNa_ss   : 1.0    : IU/mL  Baseline IFN-alpha (quiescent)
IFNa_SLE  : 8.0    : IU/mL  SLE disease IFN-alpha (relative units)
k_IFNA    : 0.20   : /day   IFN-alpha clearance
k_pDC_IFN : 0.50  : /day   pDC → IFN-α production coupling
EC50_TLR_HCQ : 500 : ng/mL HCQ blood conc EC50 for TLR inhibition (whole blood ~1000 ng/mL target)
Imax_HCQ  : 0.85   : -      HCQ maximum TLR inhibition fraction
k_IFN_score : 1.2  : -      IFN score sensitivity coefficient (IFN→score)

// ---- BAFF-B cell axis ----
kout_Bcell : 0.030 : /day  B cell natural turnover rate
kin_Bcell  : 3.0   : cells/uL/day  B cell production rate at baseline
Bcell_ss   : 100   : cells/uL  Baseline B cell count in SLE (reduced)
k_IFN_BAFF : 0.015 : /day/(IU/mL) IFN-α → BAFF upregulation coupling
k_BAFF_Bcell : 0.10 : /day  BAFF → B cell survival amplification
EC50_BAFF_B : 2.5  : ng/mL  BAFF EC50 for B cell survival
k_Blast    : 0.05  : /day  B cell → Plasmablast differentiation rate
kout_Blast : 0.25  : /day  Plasmablast natural decay
k_LLPC     : 0.02  : /day  Plasmablast → LLPC conversion
kout_LLPC  : 0.008 : /day  LLPC slow natural turnover

// ---- Anti-dsDNA & Immune Complexes ----
k_Ab       : 0.15  : /day  LLPC → anti-dsDNA antibody production rate
kout_Ab    : 0.035 : /day  Anti-dsDNA clearance rate
Anti_ss    : 100   : IU/mL  Normal anti-dsDNA titer
Anti_SLE   : 800   : IU/mL  Typical active SLE titer
k_IC       : 0.002 : /day/(IU/mL)  Anti-dsDNA → IC formation rate
kout_IC    : 0.15  : /day  IC clearance rate (C1q/phagocytosis)

// ---- Complement Consumption ----
C3_normal   : 1.0  : g/L   Normal C3 serum level
C4_normal   : 0.25 : g/L   Normal C4 serum level
kin_C3      : 0.12 : g/L/day  C3 liver synthesis rate
kin_C4      : 0.030 : g/L/day  C4 liver synthesis rate
kout_C3_base : 0.12 : /day   C3 baseline consumption/turnover
kout_C4_base : 0.12 : /day   C4 baseline consumption/turnover
k_IC_C3     : 0.15 : /day/(AU)  IC-driven C3 consumption
k_IC_C4     : 0.10 : /day/(AU)  IC-driven C4 consumption

// ---- Lupus Nephritis: Proteinuria & eGFR ----
UPCR_baseline : 0.05 : mg/mg  Normal UPCR
UPCR_SLE   : 2.0    : mg/mg  Active LN UPCR
k_IC_prot  : 0.80   : /day/(AU)  IC burden → proteinuria driving force
kout_prot  : 0.35   : /day   Proteinuria resolution rate
eGFR_ss    : 75     : mL/min/1.73m2  Baseline eGFR in active LN
k_prot_GFR : 0.008  : mL/min/1.73m2/(mg/mg)  Proteinuria→GFR worsening
k_GFR_recov: 0.003  : /day   GFR recovery rate (slow fibrosis reversal)
GFR_min    : 15     : mL/min/1.73m2  Floor for GFR (ESRD threshold marker)

$CMT @annotated
// HCQ PK
HCQ_gut     : HCQ gut compartment (mg)
HCQ_cent    : HCQ central compartment (mg)
// Belimumab PK+TMDD
BELI_cent   : Belimumab central (ug/mL * V1)
BELI_periph : Belimumab peripheral (ug/mL * V2)
BAFF_free   : Free BAFF (ng/mL; TMDD target)
BELI_complex: Belimumab-BAFF complex (ng/mL equivalents)
// Anifrolumab PK
ANIF_cent   : Anifrolumab central (ug/mL * V)
// MMF/MPA PK
MMF_gut     : MMF gut compartment (mg)
MPA_cent    : MPA central compartment (mg)
// Voclosporin PK
VOC_gut     : Voclosporin gut compartment (mg)
VOC_cent    : Voclosporin central compartment (mg)
// IFN Pathway PD
TLR_act     : TLR7/9 activation signal (relative units, RU)
IFNa_conc   : IFN-alpha concentration (IU/mL relative)
IFNscore    : IFN gene signature score (z-score units)
// BAFF-B cell axis
Bcell_count : Peripheral B cell count (cells/uL)
Plasmablast : Circulating plasmablast count (relative units)
LLPC_count  : Long-lived plasma cell pool (relative units)
// Autoantibodies & IC
AntiDsDNA   : Anti-dsDNA IgG titer (IU/mL)
IC_burden   : Circulating immune complex burden (arbitrary units, AU)
// Complement
C3_serum    : C3 serum concentration (g/L)
C4_serum    : C4 serum concentration (g/L)
// Kidney
Proteinuria : UPCR (mg/mg)
eGFR_ode    : eGFR (mL/min/1.73m2)

$MAIN
// HCQ: convert dose to mg, absorbed fraction
D_HCQ = dose_HCQ * F_HCQ;  // daily dose available

// Initialize states at SLE disease steady state
if (NEWIND <= 1) {
  // PK compartments start at 0 (build up with dosing)
  HCQ_gut_0    = 0;
  HCQ_cent_0   = 0;
  BELI_cent_0  = 0;
  BELI_periph_0 = 0;
  BELI_complex_0 = 0;
  BAFF_free_0  = Rin_BAFF / kout_BAFF;  // BAFF SS ~ 2.8 ng/mL
  ANIF_cent_0  = 0;
  MMF_gut_0    = 0;
  MPA_cent_0   = 0;
  VOC_gut_0    = 0;
  VOC_cent_0   = 0;
  // Disease state: SLE active at baseline
  TLR_act_0    = 1.5;   // elevated TLR activation in SLE
  IFNa_conc_0  = IFNa_SLE; // elevated IFN-α
  IFNscore_0   = 4.0;   // high IFN score (typical active SLE)
  Bcell_count_0 = 80;   // reduced B cells in active SLE (consumed)
  Plasmablast_0 = 15;   // elevated plasmablasts
  LLPC_count_0 = 40;    // elevated LLPC pool
  AntiDsDNA_0  = Anti_SLE;
  IC_burden_0  = 3.5;   // high IC burden at baseline
  C3_serum_0   = 0.6;   // low C3 (consumed)
  C4_serum_0   = 0.10;  // low C4 (consumed)
  Proteinuria_0 = UPCR_SLE;  // active LN at baseline
  eGFR_ode_0   = eGFR_ss;
}

// Derived PK concentrations
double C_HCQ   = HCQ_cent / Vd_HCQ * 1000;  // ng/mL (whole blood equivalent proxy)
double C_BELI  = BELI_cent / V1_BELI;         // ug/mL
double C_ANIF  = ANIF_cent / V_ANIF;          // ug/mL
double C_MPA   = MPA_cent / (Vd_MPA * 70);    // ug/mL
double C_VOC   = VOC_cent / Vd_VOC * 1e6;     // ng/mL

// ---- Drug Effect Calculations ----

// 1. HCQ: TLR inhibition (Emax, whole blood concentration)
double E_HCQ_TLR = Imax_HCQ * C_HCQ / (EC50_TLR_HCQ + C_HCQ);

// 2. Anifrolumab: IFNAR1 occupancy → IFN signal blockade
double ANIF_conc_ngmL = C_ANIF * 1000;  // ug/mL → ng/mL
double ANIF_occupancy = ANIF_conc_ngmL / (t50_ANIF * 150000 + ANIF_conc_ngmL);
// (rough approximation; t50_ANIF scales ~EC50 in PD units)
double E_ANIF = ANIF_occupancy;  // fractional IFNAR blockade 0-1

// 3. MPA: IMPDH inhibition → lymphocyte proliferation inhibition
double E_MPA_IMPDH = C_MPA / (0.25 + C_MPA);  // IC50 MPA on IMPDH ~0.25 ug/mL

// 4. Voclosporin: calcineurin inhibition → NFAT/IL-2/T cell suppression
double E_VOC = Imax_VOC * C_VOC / (IC50_VOC + C_VOC);

// 5. Corticosteroid Emax
double E_CS = Emax_CS * CS_dose / (EC50_CS + CS_dose);

// 6. BAFF inhibition by belimumab (TMDD effectively reduces BAFF_free)
double BAFF_total_eff = BAFF_free;  // TMDD tracks free BAFF directly
double E_BELI_BAFF = 1.0 - BAFF_free / (Rin_BAFF / kout_BAFF + BAFF_free);
// Approximation: fractional reduction of free BAFF relative to no-drug SS

// Combined immunosuppression factor (additive on key cytokine axes)
double ImmunoSuppr = 1.0 - (0.3*E_HCQ_TLR + 0.25*E_ANIF + 0.20*E_MPA_IMPDH
                             + 0.15*E_VOC + 0.10*E_CS);
ImmunoSuppr = (ImmunoSuppr < 0.05) ? 0.05 : ImmunoSuppr;

$ODE
// ===========================================================================
// HCQ PK (once-daily dosing handled via $EVENT)
// ===========================================================================
dxdt_HCQ_gut  = -Ka_HCQ * HCQ_gut;
dxdt_HCQ_cent = Ka_HCQ * HCQ_gut - (CL_HCQ / Vd_HCQ) * HCQ_cent;

// ===========================================================================
// Belimumab PK: 2-compartment + TMDD
// ===========================================================================
double k10_B = CL_BELI / V1_BELI;
double k12_B = Q_BELI  / V1_BELI;
double k21_B = Q_BELI  / V2_BELI;

dxdt_BELI_cent   = -(k10_B + k12_B) * BELI_cent + k21_B * BELI_periph
                    - kon_BELI * (BELI_cent/V1_BELI) * BAFF_free * V1_BELI
                    + koff_BELI * BELI_complex * V1_BELI;
dxdt_BELI_periph = k12_B * BELI_cent - k21_B * BELI_periph;

// TMDD: free BAFF
dxdt_BAFF_free   = Rin_BAFF - kout_BAFF * BAFF_free
                    - kon_BELI * (BELI_cent/V1_BELI) * BAFF_free
                    + koff_BELI * BELI_complex
                    + k_IFN_BAFF * IFNa_conc;  // IFN→BAFF upregulation

dxdt_BELI_complex = kon_BELI * (BELI_cent/V1_BELI) * BAFF_free
                    - (koff_BELI + kint_BELI) * BELI_complex;

// ===========================================================================
// Anifrolumab PK: 1-compartment IV
// ===========================================================================
dxdt_ANIF_cent = -(CL_ANIF / V_ANIF) * ANIF_cent;

// ===========================================================================
// MMF/MPA PK
// ===========================================================================
dxdt_MMF_gut  = -Ka_MMF * MMF_gut;
dxdt_MPA_cent = Ka_MMF * MMF_gut * F_MMF - (CL_MPA / (Vd_MPA * 70)) * MPA_cent;

// ===========================================================================
// Voclosporin PK
// ===========================================================================
dxdt_VOC_gut  = -Ka_VOC * VOC_gut;
dxdt_VOC_cent = Ka_VOC * VOC_gut * F_VOC - (CL_VOC / Vd_VOC) * VOC_cent;

// ===========================================================================
// IFN Pathway PD
// ===========================================================================
// TLR activation: disease-driven signal (NETs/chromatin → TLR7/9)
// HCQ blocks TLR endosomal signaling; CS broadly suppresses
double TLR_drive = 1.5;  // baseline disease drive in active SLE
double TLR_inh   = E_HCQ_TLR + 0.3 * E_CS;
TLR_inh = (TLR_inh > 0.95) ? 0.95 : TLR_inh;

dxdt_TLR_act = kin_TLR * TLR_drive * (1.0 - TLR_inh) - kout_TLR * TLR_act;

// IFN-α: driven by TLR/pDC; blocked by anifrolumab (downstream feedback)
// Anifrolumab blocks IFNAR signaling (not IFN-α itself but its effects)
double IFN_prod = k_pDC_IFN * TLR_act;
double IFN_signal_block = E_ANIF;  // anifrolumab blocks downstream IFNAR signal

dxdt_IFNa_conc = IFN_prod - k_IFNA * IFNa_conc;

// IFN score: reflects ISG expression; directly proportional to IFNa
// Anifrolumab blocks IFNAR → suppresses IFN score (even while IFNa may persist)
double IFN_effective = IFNa_conc * (1.0 - E_ANIF);
dxdt_IFNscore = k_IFN_score * IFN_effective - 0.20 * IFNscore;

// ===========================================================================
// B Cell Axis
// ===========================================================================
// BAFF effect on B cell survival (Emax)
double BAFF_eff = BAFF_free / (EC50_BAFF_B + BAFF_free);

// B cell count: production + BAFF-driven survival - turnover - MPA suppression
double Bcell_prod = kin_Bcell * (1.0 - 0.5 * E_MPA_IMPDH) * (1.0 - 0.3 * E_CS);
double Bcell_death = kout_Bcell * Bcell_count * (1.0 - k_BAFF_Bcell * BAFF_eff);

dxdt_Bcell_count = Bcell_prod - Bcell_death;

// Plasmablast: from B cell differentiation; IFN-α amplifies
double k_blast_eff = k_Blast * (1.0 + 0.5 * IFNa_conc / IFNa_SLE)
                     * (1.0 - 0.4 * E_MPA_IMPDH) * (1.0 - 0.3 * E_CS);
dxdt_Plasmablast = k_blast_eff * Bcell_count - kout_Blast * Plasmablast;

// LLPC: from plasmablasts; calcineurin inhibition impairs IL-2-driven differentiation
double k_LLPC_eff = k_LLPC * (1.0 - 0.5 * E_VOC) * (1.0 - 0.2 * E_CS);
dxdt_LLPC_count  = k_LLPC_eff * Plasmablast - kout_LLPC * LLPC_count;

// ===========================================================================
// Autoantibodies & Immune Complexes
// ===========================================================================
dxdt_AntiDsDNA = k_Ab * LLPC_count - kout_Ab * AntiDsDNA;

dxdt_IC_burden = k_IC * AntiDsDNA - kout_IC * IC_burden;

// ===========================================================================
// Complement Consumption
// ===========================================================================
double C3_consume = (kout_C3_base + k_IC_C3 * IC_burden) * C3_serum;
dxdt_C3_serum = kin_C3 - C3_consume;

double C4_consume = (kout_C4_base + k_IC_C4 * IC_burden) * C4_serum;
dxdt_C4_serum = kin_C4 - C4_consume;

// ===========================================================================
// Lupus Nephritis: Proteinuria & eGFR
// ===========================================================================
// Proteinuria driven by IC burden → glomerular injury
// Voclosporin: stabilizes podocyte cytoskeleton (additional effect)
double prot_drive = k_IC_prot * IC_burden;
double prot_resol = kout_prot * Proteinuria * (1.0 + 1.5 * E_VOC + 0.5 * E_CS);

dxdt_Proteinuria = prot_drive - prot_resol;
if (Proteinuria < 0.04) dxdt_Proteinuria = 0;  // floor at near-normal

// eGFR: declines with sustained proteinuria; partial recovery if proteinuria resolves
double GFR_loss = k_prot_GFR * Proteinuria * eGFR_ode;
double GFR_gain = k_GFR_recov * (90.0 - eGFR_ode);  // recovery toward 90 if treated
dxdt_eGFR_ode = GFR_gain - GFR_loss;
if (eGFR_ode <= GFR_min) dxdt_eGFR_ode = 0;

$TABLE
// Derived PK concentrations for output
double C_HCQ_ng   = HCQ_cent / Vd_HCQ * 1000;    // HCQ ng/mL (whole blood proxy)
double C_BELI_ug  = BELI_cent / V1_BELI;           // Belimumab ug/mL
double C_ANIF_ug  = ANIF_cent / V_ANIF;            // Anifrolumab ug/mL
double C_MPA_ug   = MPA_cent / (3.6 * 70);         // MPA ug/mL
double C_VOC_ng   = VOC_cent / Vd_VOC * 1e6;       // Voclosporin ng/mL

// SLEDAI-2K simplified dynamic score (0-105 range approximation)
// Each SLE domain weighted: anti-dsDNA (4), complement C3/C4 (2+2), proteinuria (4),
// arthritis (2), immunologic (2), ...; simplified to 3 main drivers here
double dsDNA_pts  = (AntiDsDNA > 200) ? 4.0 : ((AntiDsDNA > 100) ? 2.0 : 0.0);
double C3C4_pts   = ((C3_serum < 0.8) ? 2.0 : 0.0) + ((C4_serum < 0.16) ? 2.0 : 0.0);
double renal_pts  = (Proteinuria > 0.5) ? 4.0 : ((Proteinuria > 0.2) ? 2.0 : 0.0);
double SLEDAI_calc = dsDNA_pts + C3C4_pts + renal_pts + 4.0; // +4 baseline CNS/skin
SLEDAI_calc = (SLEDAI_calc > 30) ? 30.0 : SLEDAI_calc;

// IFN score normalized (0-10 scale for display)
double IFN_score_norm = IFNscore * 2.5;

// B cell percent reference range marker
double Bcell_pct = Bcell_count / 250.0 * 100;  // % relative to healthy 250/uL

// Renal response classification
double CR_renal  = (Proteinuria < 0.5 && eGFR_ode > 60) ? 1.0 : 0.0;
double PR_renal  = (Proteinuria < 1.5 && (Proteinuria < UPCR_SLE * 0.5)) ? 1.0 : 0.0;

capture C_HCQ_ng C_BELI_ug C_ANIF_ug C_MPA_ug C_VOC_ng
capture SLEDAI_calc IFN_score_norm Bcell_pct
capture AntiDsDNA C3_serum C4_serum Proteinuria eGFR_ode
capture BAFF_free IC_burden CR_renal PR_renal

$CAPTURE
C_HCQ_ng C_BELI_ug C_ANIF_ug C_MPA_ug C_VOC_ng
SLEDAI_calc IFN_score_norm Bcell_pct
AntiDsDNA C3_serum C4_serum Proteinuria eGFR_ode
BAFF_free IC_burden CR_renal PR_renal
'

# ==============================================================================
# COMPILE MODEL
# ==============================================================================
sle_mod <- mrgsolve::mcode("SLE_QSP", sle_model_code)

cat("Model compiled successfully.\n")
cat("Compartments:", length(init(sle_mod)), "\n")

# ==============================================================================
# DOSING REGIMENS
# ==============================================================================

# Helper: create dosing events (mrgsolve ev object)
make_HCQ_dose <- function(start_day = 0, n_days = 365, dose_mg = 200) {
  # 200 mg BID → 400 mg/day in two doses
  ev(amt = dose_mg, cmt = "HCQ_gut", ii = 12, addl = 2*n_days - 1, time = start_day * 24) %>%
    mrgsolve::as_data_frame()
}

# Scenarios (all times in HOURS; model time unit = hours for PK, days for PD)
# We use time in days (model parameterized in days for PD; PK rates also in /day)
# Use time = days throughout

week_hrs <- function(w) w * 7  # weeks → days

# Scenario definitions
scenarios <- list(
  "No Treatment (SLE active)" = list(
    HCQ = 0, BELI = 0, ANIF = 0, MMF = 0, VOC = 0, CS = 0
  ),
  "HCQ monotherapy" = list(
    HCQ = 400, BELI = 0, ANIF = 0, MMF = 0, VOC = 0, CS = 0
  ),
  "HCQ + MMF (standard of care, mild-moderate)" = list(
    HCQ = 400, BELI = 0, ANIF = 0, MMF = 3000, VOC = 0, CS = 5
  ),
  "HCQ + MMF + Belimumab (biologic add-on)" = list(
    HCQ = 400, BELI = 10, ANIF = 0, MMF = 3000, VOC = 0, CS = 5
  ),
  "HCQ + MMF + Anifrolumab (IFN-high)" = list(
    HCQ = 400, BELI = 0, ANIF = 300, MMF = 3000, VOC = 0, CS = 5
  ),
  "HCQ + MMF + Voclosporin (lupus nephritis)" = list(
    HCQ = 400, BELI = 0, ANIF = 0, MMF = 3000, VOC = 23.7 * 2, CS = 10
  ),
  "Triple therapy LN (MMF + Belimumab + Voclosporin)" = list(
    HCQ = 400, BELI = 10, ANIF = 0, MMF = 3000, VOC = 23.7 * 2, CS = 5
  )
)

# ==============================================================================
# SIMULATION FUNCTION
# ==============================================================================

simulate_scenario <- function(scen_params, scen_name, sim_days = 365) {

  # Build dose events
  ev_list <- list()

  # HCQ: BID dosing (every 12 h)
  if (scen_params$HCQ > 0) {
    hcq_amt <- scen_params$HCQ / 2  # split into BID
    ev_list[["HCQ"]] <- ev(amt = hcq_amt, cmt = "HCQ_gut",
                            ii = 0.5,  # 0.5 days = 12h (day units)
                            addl = sim_days * 2 - 1,
                            time = 0)
  }

  # Belimumab: IV q4w (every 28 days), dose = 10 mg/kg × 70 kg = 700 mg
  # administered as bolus to BELI_cent directly
  if (scen_params$BELI > 0) {
    beli_dose_total <- scen_params$BELI * 70  # mg/kg × 70 kg body weight
    n_beli_doses <- floor(sim_days / 28)
    ev_list[["BELI"]] <- ev(amt = beli_dose_total * 1000 / V1_BELI_val,  # convert
                             # Actually just add to central in ug: amt in ug/mL * L = ug
                             # dose_BELI = mg → mg/V1_BELI = mg/5.29L = ug/mL * 5.29
                             # Simplification: use parameter dose_BELI directly
                             amt = beli_dose_total,
                             cmt = "BELI_cent",
                             ii = 28,
                             addl = n_beli_doses - 1,
                             time = 0)
  }

  # Anifrolumab: 300 mg IV q4w
  if (scen_params$ANIF > 0) {
    n_anif_doses <- floor(sim_days / 28)
    ev_list[["ANIF"]] <- ev(amt = scen_params$ANIF,
                             cmt = "ANIF_cent",
                             ii = 28,
                             addl = n_anif_doses - 1,
                             time = 0)
  }

  # MMF: BID dosing → MPA via MMF_gut
  if (scen_params$MMF > 0) {
    mmf_bid <- scen_params$MMF / 2
    ev_list[["MMF"]] <- ev(amt = mmf_bid, cmt = "MMF_gut",
                            ii = 0.5,
                            addl = sim_days * 2 - 1,
                            time = 0)
  }

  # Voclosporin: BID dosing
  if (scen_params$VOC > 0) {
    voc_bid <- scen_params$VOC / 2
    ev_list[["VOC"]] <- ev(amt = voc_bid, cmt = "VOC_gut",
                            ii = 0.5,
                            addl = sim_days * 2 - 1,
                            time = 0)
  }

  # Combine dose events
  if (length(ev_list) > 0) {
    dose_event <- Reduce(mrgsolve::c, ev_list)
  } else {
    dose_event <- NULL
  }

  # Update model parameters
  mod_upd <- param(sle_mod,
                   CS_dose = scen_params$CS)

  # Simulation time: daily output
  tgrid <- seq(0, sim_days, by = 1)

  # Run
  if (!is.null(dose_event)) {
    out <- mrgsim(mod_upd, events = dose_event, tgrid = tgrid,
                  recover = "scenario", delta = 0.25) %>%
      as_tibble() %>%
      mutate(scenario = scen_name)
  } else {
    out <- mrgsim(mod_upd, tgrid = tgrid, delta = 0.25) %>%
      as_tibble() %>%
      mutate(scenario = scen_name)
  }

  return(out)
}

# ==============================================================================
# RUN ALL SCENARIOS
# ==============================================================================

cat("\nRunning", length(scenarios), "treatment scenarios over 365 days...\n")

results_list <- lapply(names(scenarios), function(sn) {
  cat("  Simulating:", sn, "\n")
  tryCatch(
    simulate_scenario(scenarios[[sn]], sn, sim_days = 365),
    error = function(e) {
      message("  Warning - ", sn, ": ", e$message)
      NULL
    }
  )
})

# Filter out NULLs and combine
results_list <- Filter(Negate(is.null), results_list)

if (length(results_list) > 0) {
  results_all <- bind_rows(results_list)
  cat("Simulation complete. Total rows:", nrow(results_all), "\n")
} else {
  stop("All simulations failed - please check mrgsolve model code.")
}

# ==============================================================================
# PHARMACOKINETIC PLOTS
# ==============================================================================

plot_pk <- function(data) {
  # HCQ PK
  p1 <- data %>%
    filter(scenario %in% grep("HCQ", unique(data$scenario), value = TRUE)) %>%
    ggplot(aes(x = time, y = C_HCQ_ng, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "HCQ Whole Blood Concentration",
         x = "Time (days)", y = "HCQ (ng/mL)", color = "Scenario") +
    geom_hline(yintercept = c(500, 1000, 1200), linetype = "dashed",
               color = c("orange", "green", "red"), alpha = 0.7) +
    annotate("text", x = 350, y = 1050, label = "Therapeutic target\n1000-1200 ng/mL",
             size = 3, color = "darkgreen") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  # Belimumab PK
  p2 <- data %>%
    filter(grepl("Belimumab", scenario)) %>%
    ggplot(aes(x = time, y = C_BELI_ug, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "Belimumab Serum Concentration",
         x = "Time (days)", y = "Belimumab (µg/mL)", color = "Scenario") +
    theme_bw() + theme(legend.position = "bottom")

  # Anifrolumab PK
  p3 <- data %>%
    filter(grepl("Anifrolumab", scenario)) %>%
    ggplot(aes(x = time, y = C_ANIF_ug, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "Anifrolumab Serum Concentration",
         x = "Time (days)", y = "Anifrolumab (µg/mL)", color = "Scenario") +
    theme_bw() + theme(legend.position = "bottom")

  # MPA PK
  p4 <- data %>%
    filter(grepl("MMF", scenario)) %>%
    ggplot(aes(x = time, y = C_MPA_ug, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "MPA (Mycophenolic Acid) Concentration",
         x = "Time (days)", y = "MPA (µg/mL)", color = "Scenario") +
    geom_hline(yintercept = c(1, 3.5), linetype = "dashed",
               color = c("orange", "green")) +
    theme_bw() + theme(legend.position = "bottom")

  (p1 + p2) / (p3 + p4) +
    plot_annotation(title = "SLE Drug PK Profiles",
                    subtitle = "Key drug concentrations over 52 weeks")
}

# ==============================================================================
# PHARMACODYNAMIC PLOTS
# ==============================================================================

plot_ifn_axis <- function(data) {
  p1 <- ggplot(data, aes(x = time, y = IFN_score_norm, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "IFN Gene Signature Score",
         x = "Time (days)", y = "IFN Score (normalized)", color = NULL) +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray50") +
    annotate("text", x = 300, y = 2.2, label = "Low IFN threshold", size = 3) +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p2 <- ggplot(data, aes(x = time, y = BAFF_free, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "Free BAFF/BLyS (TMDD)",
         x = "Time (days)", y = "BAFF (ng/mL)", color = NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p1 / p2
}

plot_bcell_ab <- function(data) {
  p1 <- ggplot(data, aes(x = time, y = Bcell_pct, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "B Cell Count (% relative to normal)",
         x = "Time (days)", y = "B cells (%)", color = NULL) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p2 <- ggplot(data, aes(x = time, y = AntiDsDNA, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "Anti-dsDNA IgG Titer",
         x = "Time (days)", y = "Anti-dsDNA (IU/mL)", color = NULL) +
    geom_hline(yintercept = 200, linetype = "dashed", color = "red") +
    annotate("text", x = 300, y = 220, label = "Positive threshold", size = 3, color = "red") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p1 / p2
}

plot_complement <- function(data) {
  p1 <- ggplot(data, aes(x = time, y = C3_serum, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "C3 Serum Level",
         x = "Time (days)", y = "C3 (g/L)", color = NULL) +
    geom_hline(yintercept = c(0.9, 1.8), linetype = "dashed",
               color = c("red", "gray60")) +
    annotate("text", x = 300, y = 0.75, label = "Low (<0.9)", size = 3, color = "red") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p2 <- ggplot(data, aes(x = time, y = C4_serum, color = scenario)) +
    geom_line(linewidth = 0.8) +
    labs(title = "C4 Serum Level",
         x = "Time (days)", y = "C4 (g/L)", color = NULL) +
    geom_hline(yintercept = 0.16, linetype = "dashed", color = "red") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p1 / p2
}

plot_nephritis <- function(data) {
  p1 <- ggplot(data, aes(x = time, y = Proteinuria, color = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = c(0.5, 1.0), linetype = c("dashed", "dotted"),
               color = c("green4", "orange2")) +
    annotate("text", x = 280, y = 0.55, label = "CR threshold (<0.5)", size = 3, color = "green4") +
    annotate("text", x = 280, y = 1.05, label = "PR threshold (<1.0)", size = 3, color = "orange2") +
    labs(title = "Proteinuria (UPCR)",
         x = "Time (days)", y = "UPCR (mg/mg)", color = NULL) +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p2 <- ggplot(data, aes(x = time, y = eGFR_ode, color = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = c(60, 30), linetype = "dashed",
               color = c("orange", "red")) +
    labs(title = "eGFR Trajectory",
         x = "Time (days)", y = "eGFR (mL/min/1.73m²)", color = NULL) +
    annotate("text", x = 300, y = 65, label = "CKD G2 (<60)", size = 3, color = "orange") +
    ylim(0, 100) +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))

  p1 / p2
}

plot_sledai <- function(data) {
  ggplot(data, aes(x = time, y = SLEDAI_calc, color = scenario)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = c(4, 6), linetype = c("dashed", "dotted"),
               color = c("green4", "orange")) +
    annotate("text", x = 280, y = 4.5, label = "LLDAS (≤4)", size = 3, color = "green4") +
    annotate("text", x = 280, y = 6.5, label = "Mild activity", size = 3, color = "orange") +
    labs(title = "SLEDAI-2K Dynamic Score",
         subtitle = "Composite of anti-dsDNA, complement, renal, constitutional domains",
         x = "Time (days)", y = "SLEDAI-2K", color = "Treatment Scenario") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
}

# ==============================================================================
# DOSE-RESPONSE: BELIMUMAB EFFECT ON ANTI-dsDNA AT WEEK 52
# ==============================================================================

dr_belimumab <- function() {
  doses <- c(0, 1, 3, 5, 10, 15) * 70  # mg/kg × 70kg
  dr_results <- lapply(doses, function(d) {
    mod_dr <- param(sle_mod, CS_dose = 5)
    if (d > 0) {
      ev_dr <- ev(amt = d, cmt = "BELI_cent", ii = 28, addl = 12, time = 0)
      out <- mrgsim(mod_dr, events = ev_dr, end = 365, delta = 1) %>%
        as_tibble() %>% filter(time == 364) %>%
        select(AntiDsDNA, SLEDAI_calc, Proteinuria) %>%
        mutate(dose_mgkg = d / 70)
    } else {
      out <- mrgsim(mod_dr, end = 365, delta = 1) %>%
        as_tibble() %>% filter(time == 364) %>%
        select(AntiDsDNA, SLEDAI_calc, Proteinuria) %>%
        mutate(dose_mgkg = 0)
    }
    out
  })
  bind_rows(dr_results)
}

cat("\nRunning belimumab dose-response analysis...\n")
dr_beli <- tryCatch(dr_belimumab(), error = function(e) {
  message("DR analysis skipped: ", e$message); NULL
})

# ==============================================================================
# GENERATE ALL PLOTS (if simulation succeeded)
# ==============================================================================

if (exists("results_all") && nrow(results_all) > 0) {
  cat("\nGenerating plots...\n")

  p_sledai   <- plot_sledai(results_all)
  p_ifn      <- plot_ifn_axis(results_all)
  p_bcell_ab <- plot_bcell_ab(results_all)
  p_compl    <- plot_complement(results_all)
  p_nephritis <- plot_nephritis(results_all)

  # Assemble dashboard
  dashboard <- (p_sledai) /
    (p_ifn[[1]] + p_bcell_ab[[1]]) /
    (p_compl[[1]] + p_nephritis[[1]]) +
    plot_annotation(
      title = "SLE QSP Model — Treatment Scenario Dashboard",
      subtitle = paste0("7 treatment regimens over 52 weeks | ",
                        "IFN-α axis, BAFF-B cell axis, LN endpoints"),
      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                    plot.subtitle = element_text(size = 11))
    )

  # Print
  print(p_sledai)
  cat("\nKey endpoint summary at Week 52:\n")
  summary_tbl <- results_all %>%
    filter(time == 364) %>%
    select(scenario, SLEDAI_calc, AntiDsDNA, C3_serum, C4_serum,
           Proteinuria, eGFR_ode, IFN_score_norm) %>%
    mutate(across(where(is.double), ~ round(.x, 2)))
  print(as.data.frame(summary_tbl))
}

# ==============================================================================
# POPULATION VARIABILITY (Monte Carlo, N = 200)
# ==============================================================================

pop_sim <- function(n_subj = 200, scenario_name = "HCQ + MMF + Belimumab (biologic add-on)") {

  # IIV parameters (% CV → lognormal)
  set.seed(20260616)
  iiv_cv <- data.frame(
    ID          = 1:n_subj,
    CL_BELI     = rlnorm(n_subj, 0, 0.28),   # 28% CV
    Rin_BAFF    = rlnorm(n_subj, 0, 0.35),
    UPCR_SLE    = rlnorm(n_subj, 0, 0.40),
    Anti_SLE    = rlnorm(n_subj, 0, 0.45),
    eGFR_ss     = rnorm(n_subj, 75, 15) %>% pmax(20)
  )

  scen <- scenarios[[scenario_name]]
  dose_events <- ev(amt = 700, cmt = "BELI_cent", ii = 28, addl = 12, time = 0)
  mmf_ev      <- ev(amt = 1500, cmt = "MMF_gut", ii = 0.5, addl = 729, time = 0)
  hcq_ev      <- ev(amt = 200, cmt = "HCQ_gut", ii = 0.5, addl = 729, time = 0)
  all_ev      <- mrgsolve::c(dose_events, mmf_ev, hcq_ev)

  pop_out <- lapply(1:n_subj, function(i) {
    row <- iiv_cv[i, ]
    mod_i <- param(sle_mod,
                   CL_BELI = 0.215 * row$CL_BELI,
                   Rin_BAFF = 0.42 * row$Rin_BAFF,
                   UPCR_SLE = 2.0 * row$UPCR_SLE,
                   Anti_SLE = 800 * row$Anti_SLE,
                   eGFR_ss  = row$eGFR_ss,
                   CS_dose  = scen$CS)
    mrgsim(mod_i, events = all_ev, end = 365, delta = 4) %>%
      as_tibble() %>%
      mutate(ID = i)
  })
  bind_rows(pop_out)
}

cat("\nRunning population simulation (N=200)...\n")
pop_results <- tryCatch(
  pop_sim(200),
  error = function(e) { message("Pop sim skipped: ", e$message); NULL }
)

if (!is.null(pop_results)) {
  cat("Population simulation complete:", nrow(pop_results), "rows\n")

  # 5th/50th/95th percentile ribbon
  pop_summary <- pop_results %>%
    group_by(time) %>%
    summarise(
      p5_prot  = quantile(Proteinuria, 0.05, na.rm = TRUE),
      p50_prot = quantile(Proteinuria, 0.50, na.rm = TRUE),
      p95_prot = quantile(Proteinuria, 0.95, na.rm = TRUE),
      p5_egfr  = quantile(eGFR_ode, 0.05, na.rm = TRUE),
      p50_egfr = quantile(eGFR_ode, 0.50, na.rm = TRUE),
      p95_egfr = quantile(eGFR_ode, 0.95, na.rm = TRUE),
      cr_rate  = mean(CR_renal, na.rm = TRUE),
      .groups = "drop"
    )

  p_pop_prot <- ggplot(pop_summary, aes(x = time)) +
    geom_ribbon(aes(ymin = p5_prot, ymax = p95_prot), fill = "#2196F3", alpha = 0.25) +
    geom_line(aes(y = p50_prot), color = "#1565C0", linewidth = 1.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "green4") +
    labs(title = "Population Variability: Proteinuria (UPCR)",
         subtitle = "Median ± 5th–95th percentile | N=200 virtual patients",
         x = "Time (days)", y = "UPCR (mg/mg)") +
    theme_bw()

  p_pop_cr <- ggplot(pop_summary, aes(x = time, y = cr_rate * 100)) +
    geom_line(color = "#2E7D32", linewidth = 1.2) +
    geom_ribbon(aes(ymin = pmax(0, cr_rate * 100 - 5),
                    ymax = pmin(100, cr_rate * 100 + 5)),
                fill = "#4CAF50", alpha = 0.2) +
    labs(title = "Complete Renal Response Rate Over Time",
         x = "Time (days)", y = "CR Rate (%)") +
    theme_bw()

  print(p_pop_prot + p_pop_cr)
}

cat("\n=== SLE QSP Model — Simulation Complete ===\n")
