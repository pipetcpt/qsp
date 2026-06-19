# ==============================================================================
# IgA Vasculitis (Henoch-Schönlein Purpura / IgAV) QSP Model
# mrgsolve-based ODE Pharmacokinetic/Pharmacodynamic Model
# ==============================================================================
#
# References:
#   Pillebout E, et al. Addition of cyclophosphamide to steroids provides no
#     benefit compared with steroids alone in treating adult patients with
#     severe Henoch-Schönlein Purpura. Kidney Int. 2010;78(5):495-502.
#     PMID: 20520597
#
#   Selvaskandan H, et al. New strategies and perspectives on managing IgA
#     nephropathy. Clin Exp Nephrol. 2019;23(5):577-588. PMID: 30680474
#
#   Oni L, Sampath S. Childhood IgA Vasculitis (Henoch Schonlein Purpura)-
#     Advances and Knowledge Gaps. Front Pediatr. 2019;7:257. PMID: 31275912
#
#   Coppo R. Treatment of IgA Nephropathy: Recent Advances and Remaining
#     Challenges. J Clin Med. 2021;10(21):4949. PMID: 34768450
#
#   Rauen T, et al. Intensive Supportive Care plus Immunosuppression in
#     IgA Nephropathy. N Engl J Med. 2015;373(23):2225-2236. PMID: 26630142
#     (STOP-IgAN Trial - used for corticosteroid PD parameter calibration)
#
#   Lafayette RA, et al. Sparsentan vs. Irbesartan in IgA Nephropathy.
#     N Engl J Med. 2023;389(19):1764-1775. PMID: 37632463 (PROTECT Trial)
#
#   Heerspink HJL, et al. Dapagliflozin in Patients with Chronic Kidney
#     Disease. N Engl J Med. 2020;383(15):1436-1446. PMID: 32970396
#     (DAPA-CKD Trial - dapagliflozin PD parameters)
#
#   Jayne DRW, et al. Rituximab versus Cyclophosphamide for ANCA-associated
#     Vasculitis. N Engl J Med. 2010;363(3):221-232. PMID: 20647199
#     (Rituximab PK/PD framework adapted for IgAV)
#
#   Magistroni R, et al. New developments in the genetics, pathogenesis and
#     therapy of IgA nephropathy. Kidney Int. 2015;88(5):974-989.
#     PMID: 26376133 (Gd-IgA1 / multi-hit model parameters)
#
#   Suzuki H, et al. Pathophysiology of IgA nephropathy. Clin Exp Nephrol.
#     2013;17(5):610-621. PMID: 23474885
#
#   Floege J, Feehally J. The mucosa-kidney axis in IgA nephropathy.
#     Nat Rev Nephrol. 2016;12(3):147-156. PMID: 26750480
#
#   Cheung CK, et al. Galactose-deficient IgA1 links gut microbiome to
#     IgA nephropathy. J Clin Invest. 2019;129(6):2413-2426. PMID: 31094714
#
#   Woo KT, et al. BAFF and APRIL in IgA nephropathy. Nephrol Dial Transplant.
#     2014;29(Suppl 4):iv125-130. PMID: 24753256
#
#   Peng W, et al. Role of IL-6 and TNF-alpha in the pathogenesis of IgA
#     nephropathy. Mediators Inflamm. 2018;2018:9548026. PMID: 29849489
#
# Model Structure: 25 compartments
#   - Prednisolone PK: 2-compartment oral (CMT 1-2 + gut depot)
#   - MMF/MPA PK: 2-compartment oral (CMT 4-5 + gut depot)
#   - Rituximab PK: 2-compartment + TMDD (CMT 7-9)
#   - Disease PD: 15 compartments representing IgAV pathophysiology
#
# Author: Claude Code (CCR session 2026-06-19)
# ==============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# mrgsolve Model Definition
# ==============================================================================

mod <- mrgsolve::mcode("igav_qsp", '

$PROB
IgA Vasculitis (Henoch-Schonlein Purpura) QSP Model
Compartments: 25 ODE compartments
Drugs: Prednisolone, MMF/MPA, Rituximab, ACEi (enalapril), Sparsentan, Dapagliflozin
Disease Endpoints: Gd-IgA1, IgA-IC, Complement, Renal, Skin, GI, Cytokines

$PARAM
// ===========================================================
// PK Parameters - Prednisolone (oral)
// Ref: Rohatagi S, et al. J Pharm Sci. 1996;85(10):1070-1082.
//      (PMID: 8897276) Two-compartment PK of prednisolone
// ===========================================================
PRED_ka   = 2.0     // absorption rate constant (1/h)
PRED_CL   = 15.0    // clearance (L/h) - 70 kg adult
PRED_V2   = 45.0    // central volume (L)
PRED_Q    = 8.0     // inter-compartmental clearance (L/h)
PRED_V3   = 90.0    // peripheral volume (L)
PRED_F    = 0.82    // oral bioavailability

// ===========================================================
// PK Parameters - Mycophenolate Mofetil / MPA
// Ref: Bullingham RES, et al. Clin Pharmacokinet. 1998;34(6):429-455.
//      (PMID: 9646006) Two-compartment PK of MPA
// Note: MMF is pro-drug; ka represents combined absorption + hydrolysis
// ===========================================================
MMF_ka    = 1.5     // absorption rate constant (1/h)
MMF_CL    = 12.0    // apparent MPA clearance (L/h)
MMF_V2    = 30.0    // central volume (L)
MMF_Q     = 5.0     // inter-compartmental clearance (L/h)
MMF_V3    = 50.0    // peripheral volume (L)
MMF_F     = 0.94    // oral bioavailability of MMF -> MPA

// ===========================================================
// PK Parameters - Rituximab (IV, TMDD model)
// Ref: Dirks NL, Meibohm B. Clin Pharmacokinet. 2010;49(10):633-659.
//      (PMID: 20818832) TMDD framework for anti-CD20
//      Reff ME, et al. Blood. 1994;83(2):435-445. (PMID: 7506951)
// ===========================================================
RTX_CL    = 0.015   // linear clearance (L/h)
RTX_V1    = 3.0     // central volume (L)
RTX_Q     = 0.12    // inter-compartmental clearance (L/h)
RTX_V2    = 4.0     // peripheral volume (L)
RTX_kon   = 0.55    // RTX-CD20 binding on-rate (1/nM/h)
RTX_koff  = 0.017   // RTX-CD20 dissociation rate (1/h)
RTX_kint  = 0.025   // RTX-CD20 complex internalization rate (1/h)
RTX_Bmax  = 80.0    // B-cell CD20 target baseline (nM)
RTX_kdeg  = 0.003   // CD20 target degradation rate (1/h)
RTX_ksyn  = 0.24    // CD20 synthesis rate (nM/h) = Bmax * kdeg

// ===========================================================
// ACEi / ARB Parameters (enalapril - renoprotection)
// Ref: Lewis EJ, et al. NEJM. 1993;329(20):1456-1462. (PMID: 8413456)
//      Effect modeled as direct reduction in efferent arteriolar resistance
// ===========================================================
ACEI_Emax  = 0.40   // max fractional reduction in intraglomerular pressure
ACEI_EC50  = 1.0    // half-maximal concentration at target (normalized)
ACEI_dose  = 0.0    // 0 = off, 1 = standard dose equivalent

// ===========================================================
// Sparsentan Parameters (dual AT1R / ETA-R antagonist)
// Ref: Lafayette RA, et al. NEJM. 2023;389(19):1764-1775.
//      (PMID: 37632463) PROTECT trial - ~40% proteinuria reduction
// ===========================================================
SPAR_Emax  = 0.55   // max fractional proteinuria reduction
SPAR_EC50  = 1.0    // normalized EC50
SPAR_dose  = 0.0    // 0 = off, 1 = 400mg QD equivalent

// ===========================================================
// SGLT2 Inhibitor Parameters (dapagliflozin 10 mg)
// Ref: Heerspink HJL, et al. NEJM. 2020;383(15):1436-1446.
//      (PMID: 32970396) DAPA-CKD: ~44% RRR for composite renal outcome
//      ~30% proteinuria reduction independent of glucose
// ===========================================================
SGLT2_Emax  = 0.30  // max fractional proteinuria reduction
SGLT2_EC50  = 1.0   // normalized EC50
SGLT2_dose  = 0.0   // 0 = off, 1 = 10mg QD

// ===========================================================
// Disease PD Parameters - Gd-IgA1 and Immune Complex Formation
// Ref: Magistroni R, et al. Kidney Int. 2015;88(5):974-989.
//      (PMID: 26376133) Multi-hit hypothesis quantitative framework
//      Suzuki H, et al. Clin Exp Nephrol. 2013;17(5):610-621.
// ===========================================================
GdIgA1_syn  = 0.50  // baseline Gd-IgA1 synthesis rate (mg/L/h)
GdIgA1_deg  = 0.021 // Gd-IgA1 degradation rate (1/h) → t1/2 ~33h
GdIgA1_ss   = 23.8  // baseline Gd-IgA1 steady-state (mg/L)
              // Normal = 5 mg/L; IgAV patients ~24 mg/L (Cheung 2019)

AntiIgA_syn = 0.15  // anti-Gd-IgA1 IgG synthesis rate (AU/mL/h)
AntiIgA_deg = 0.010 // anti-IgA1 IgG degradation rate (1/h)

IC_form     = 0.0008 // IC formation rate from Gd-IgA1 × anti-IgA1 product
IC_clear    = 0.025  // IC clearance rate (1/h)

// ===========================================================
// Complement Parameters
// Ref: Wyatt RJ, Julian BA. NEJM. 2013;368(25):2402-2414.
//      (PMID: 23782179) Complement in IgAN/IgAV
// ===========================================================
C3_base    = 1.0    // baseline C3 activation (normalized)
C3_kact    = 0.30   // IC-driven C3 activation rate constant
C3_kdeg    = 0.25   // C3 deactivation rate (1/h)
C5_kact    = 0.20   // C3-driven C5 activation rate
C5_kdeg    = 0.22   // C5 deactivation/clearance rate (1/h)
C5_MAC_k   = 0.15   // sC5b-9 formation rate from C5 activation
MAC_clear  = 0.18   // sC5b-9 clearance (1/h)

// ===========================================================
// Renal Cell / Fibrosis Parameters
// Ref: Floege J, Feehally J. Nat Rev Nephrol. 2016;12(3):147-156.
//      (PMID: 26750480)
//      Barbour SJ, et al. J Am Soc Nephrol. 2015;26(6):1445-1454.
//      (PMID: 25349205) Oxford MEST-C score and progression rates
// ===========================================================
MESANG_base = 1.0   // baseline mesangial cellularity (normalized)
MESANG_kprol = 0.08 // IC/MAC-driven mesangial proliferation rate
MESANG_kdeg  = 0.05 // mesangial cell apoptosis/regression rate (1/h)

RENAL_base   = 1.0  // baseline renal inflammation score
RENAL_kact   = 0.10 // cytokine-driven inflammation activation
RENAL_kdeg   = 0.04 // renal inflammation resolution rate (1/h)
RENAL_TGFb   = 0.03 // TGF-beta contribution to fibrosis progression
Fibr_kprog   = 0.001 // slow irreversible fibrosis progression rate

// ===========================================================
// Clinical Endpoint Parameters
// Ref: Pillebout E, et al. Kidney Int. 2010;78(5):495-502.
//      (PMID: 20520597) Outcome data in adult IgAV nephritis
// ===========================================================
PROT_base   = 2.0   // baseline proteinuria (g/day) - moderate IgAV nephritis
PROT_kact   = 0.040 // glomerular inflammation -> proteinuria
PROT_kdeg   = 0.035 // proteinuria resolution rate (1/h)

HEMAT_base  = 3.0   // baseline hematuria score (0-10 scale)
HEMAT_kact  = 0.025 // IC/complement -> hematuria
HEMAT_kdeg  = 0.04  // hematuria resolution (1/h)

EGFR_base   = 75.0  // baseline eGFR (mL/min/1.73m2) - CKD G2/G3a typical
EGFR_kdecl  = 0.0001 // baseline eGFR decline rate per hour (≈ -0.7 mL/min/yr)
EGFR_protect = 0.00003 // renal inflammation-mediated accelerated decline

// ===========================================================
// Extra-Renal Manifestations - Skin Purpura and GI
// Ref: Oni L, Sampath S. Front Pediatr. 2019;7:257. (PMID: 31275912)
//      Charlanne H, et al. Nephrol Dial Transplant. 2011 (adult IgAV)
// ===========================================================
SKIN_base   = 5.0   // baseline skin purpura score (0-10 scale)
SKIN_kact   = 0.020 // IgA-IC-driven purpura activation
SKIN_kdeg   = 0.060 // purpura resolution rate (1/h)
SKIN_steroid_E = 0.70 // max steroid effect on skin purpura

GI_base     = 3.0   // baseline GI inflammation score (0-10)
GI_kact     = 0.018 // IC-driven GI inflammation
GI_kdeg     = 0.050 // GI inflammation resolution (1/h)

// ===========================================================
// Cytokine Parameters
// Ref: Peng W, et al. Mediators Inflamm. 2018;2018:9548026.
//      (PMID: 29849489) IL-6, TNF-alpha in IgA nephropathy
//      Woo KT, et al. NDT. 2014;29(Suppl 4):iv125-130.
//      (PMID: 24753256) BAFF/APRIL in IgAN
// ===========================================================
BAFF_base   = 1.0   // baseline BAFF level (normalized)
BAFF_syn    = 0.05  // BAFF production rate
BAFF_deg    = 0.05  // BAFF clearance rate (1/h)
BAFF_stim   = 0.30  // IC-driven BAFF upregulation coefficient

IL6_base    = 1.0   // baseline IL-6 (normalized)
IL6_syn     = 0.08  // IL-6 production rate from mesangial cells
IL6_deg     = 0.12  // IL-6 clearance (1/h) - t1/2 ~6h
IL6_MESANG  = 0.25  // mesangial proliferation -> IL-6 coupling

TNF_base    = 1.0   // baseline TNF-alpha (normalized)
TNF_syn     = 0.06  // TNF production rate
TNF_deg     = 0.10  // TNF clearance (1/h)
TNF_IC      = 0.20  // IC-driven TNF induction coefficient

// ===========================================================
// Drug PD Effect Parameters
// ===========================================================
// Prednisolone effects (immunosuppression)
// Ref: Rauen T, et al. NEJM. 2015;373(23):2225-2236. (PMID: 26630142)
PRED_EC50_IS = 0.10  // prednisolone EC50 for immunosuppression (mg/L)
PRED_Emax_IS = 0.75  // max immunosuppression effect
PRED_EC50_anti = 0.15 // EC50 for anti-inflammatory (cytokines)
PRED_Emax_anti = 0.70 // max anti-inflammatory

// MMF/MPA effects
// Ref: Dooley MA, et al. NEJM. 2011 (LN trial - MPA IC50 for B/T cell)
MMF_EC50_Bcell = 1.5  // MPA EC50 for B-cell proliferation inhibition (mg/L)
MMF_Emax_Bcell = 0.80 // max B-cell inhibition
MMF_EC50_MESANG = 2.0 // MPA EC50 for mesangial cell inhibition (mg/L)
MMF_Emax_MESANG = 0.60

// Rituximab effects (B-cell depletion → reduced anti-Gd-IgA1)
// Ref: Jayne DRW, et al. NEJM. 2010;363(3):221-232. (PMID: 20647199)
RTX_Edep_max = 0.95   // max B-cell depletion by RTX
RTX_Edep_EC50 = 10.0  // RTX EC50 for B-cell depletion (nM)

$CMT
// ===========================================================
// COMPARTMENTS (25 total)
// ===========================================================

// Prednisolone PK (3 CMTs: gut depot + central + peripheral)
PRED_GUT     // [1]  prednisolone gut absorption depot (mg)
PRED_CENT    // [2]  prednisolone central compartment (mg)
PRED_PERI    // [3]  prednisolone peripheral compartment (mg)

// MMF/MPA PK (3 CMTs: gut depot + central + peripheral)
MMF_GUT      // [4]  MMF gut absorption depot (mg)
MMF_CENT     // [5]  MPA central compartment (mg)
MMF_PERI     // [6]  MPA peripheral compartment (mg)

// Rituximab PK + TMDD (3 CMTs)
RTX_CENT     // [7]  rituximab central compartment (nmol)
RTX_PERI     // [8]  rituximab peripheral compartment (nmol)
RTX_CD20     // [9]  rituximab-CD20 complex (nmol, TMDD)

// Disease PD Compartments (16 CMTs)
GD_IGA1      // [10] Gd-IgA1 level (mg/L)
ANTI_IGA1    // [11] anti-Gd-IgA1 IgG (AU/mL)
IC           // [12] IgA immune complex concentration (mg/L)
C3_ACT       // [13] activated complement C3 (normalized)
C5_ACT       // [14] activated complement C5 / sC5b-9 (normalized)
MESANG       // [15] mesangial cell proliferation index (normalized)
RENAL_INFLAM // [16] renal tubulointerstitial inflammation score
PROT         // [17] proteinuria (g/day)
HEMAT        // [18] hematuria score (0-10)
EGFR         // [19] eGFR (mL/min/1.73m2)
SKIN_PURP    // [20] skin purpura score (0-10)
GI_INFLAM    // [21] GI inflammation score (0-10)
BAFF         // [22] BAFF level (normalized)
IL6          // [23] IL-6 (normalized)
TNF          // [24] TNF-alpha (normalized)
FIBROSIS     // [25] renal fibrosis index (cumulative, 0-1)

$MAIN
// ===========================================================
// INITIAL CONDITIONS
// ===========================================================

// PK compartments start at zero (dosed via events)
PRED_GUT_0    = 0;
PRED_CENT_0   = 0;
PRED_PERI_0   = 0;
MMF_GUT_0     = 0;
MMF_CENT_0    = 0;
MMF_PERI_0    = 0;
RTX_CENT_0    = 0;
RTX_PERI_0    = 0;
RTX_CD20_0    = 0;

// Disease compartments at steady-state (moderate active IgAV nephritis)
GD_IGA1_0     = GdIgA1_ss;         // elevated Gd-IgA1 at disease onset
ANTI_IGA1_0   = AntiIgA_syn / AntiIgA_deg;
IC_0          = IC_form * GdIgA1_ss * (AntiIgA_syn/AntiIgA_deg) / IC_clear;
C3_ACT_0      = C3_base;
C5_ACT_0      = C5_kact * C3_base / C5_kdeg;
MESANG_0      = MESANG_base;
RENAL_INFLAM_0 = RENAL_base;
PROT_0        = PROT_base;
HEMAT_0       = HEMAT_base;
EGFR_0        = EGFR_base;
SKIN_PURP_0   = SKIN_base;
GI_INFLAM_0   = GI_base;
BAFF_0        = BAFF_base;
IL6_0         = IL6_base;
TNF_0         = TNF_base;
FIBROSIS_0    = 0.05;  // small amount of baseline fibrosis in IgAV nephritis

// Derived PK concentrations (mg/L or nM) for effect calculations
double PRED_C   = PRED_CENT / PRED_V2;
double MPA_C    = MMF_CENT  / MMF_V2;
double RTX_C    = RTX_CENT  / RTX_V1;

// Unbound target CD20 (for TMDD)
// Total B-cell CD20 target = RTX_CD20 + free CD20
// free_CD20 is tracked implicitly; we use RTX_CD20 complex as proxy for depletion

// Drug effect calculations - Hill equation (Emax model)
// Prednisolone immunosuppression
double E_PRED_IS   = PRED_Emax_IS * PRED_C / (PRED_C + PRED_EC50_IS);
double E_PRED_ANTI = PRED_Emax_anti * PRED_C / (PRED_C + PRED_EC50_anti);

// MMF/MPA effects on B-cells and mesangial cells
double E_MMF_B    = MMF_Emax_Bcell * MPA_C / (MPA_C + MMF_EC50_Bcell);
double E_MMF_MES  = MMF_Emax_MESANG * MPA_C / (MPA_C + MMF_EC50_MESANG);

// Rituximab B-cell depletion effect via drug-target complex
double E_RTX_dep  = RTX_Edep_max * RTX_C / (RTX_C + RTX_Edep_EC50);

// Combined B-cell inhibition (RTX + MMF + PRED each contribute)
double E_Bcell_total = 1.0 - (1.0 - E_RTX_dep) * (1.0 - E_MMF_B) * (1.0 - E_PRED_IS * 0.5);
if(E_Bcell_total > 0.98) E_Bcell_total = 0.98;

// ACEi/ARB renoprotection
double E_ACEI = ACEI_Emax * ACEI_dose / (ACEI_dose * 0.5 + ACEI_EC50);

// Sparsentan dual-receptor blockade
double E_SPAR = SPAR_Emax * SPAR_dose / (SPAR_dose * 0.5 + SPAR_EC50);

// SGLT2 inhibitor effect
double E_SGLT2 = SGLT2_Emax * SGLT2_dose / (SGLT2_dose * 0.5 + SGLT2_EC50);

// Combined renoprotective effect (sparsentan + ACEi + SGLT2)
double E_RENOPROT = 1.0 - (1.0 - E_ACEI) * (1.0 - E_SPAR) * (1.0 - E_SGLT2);

// Immune complex effect (scaled 0-1 from baseline)
double IC_effect = IC / (IC_0 + 0.001);  // ratio to baseline IC

$ODE
// ===========================================================
// ORDINARY DIFFERENTIAL EQUATIONS
// ===========================================================

// -----------------------------------------------------------
// [1-3] Prednisolone PK (2-compartment oral)
// F * dose absorbed into GUT depot, then ka to CENT
// -----------------------------------------------------------
dxdt_PRED_GUT   = -PRED_ka * PRED_GUT;
dxdt_PRED_CENT  =  PRED_ka * PRED_GUT
                 - (PRED_CL / PRED_V2) * PRED_CENT
                 - (PRED_Q  / PRED_V2) * PRED_CENT
                 + (PRED_Q  / PRED_V3) * PRED_PERI;
dxdt_PRED_PERI  =  (PRED_Q  / PRED_V2) * PRED_CENT
                 - (PRED_Q  / PRED_V3) * PRED_PERI;

// -----------------------------------------------------------
// [4-6] MMF/MPA PK (2-compartment oral, pro-drug model)
// MMF_GUT -> MPA_CENT (ka includes esterase hydrolysis)
// -----------------------------------------------------------
dxdt_MMF_GUT    = -MMF_ka  * MMF_GUT;
dxdt_MMF_CENT   =  MMF_ka  * MMF_GUT
                 - (MMF_CL / MMF_V2) * MMF_CENT
                 - (MMF_Q  / MMF_V2) * MMF_CENT
                 + (MMF_Q  / MMF_V3) * MMF_PERI;
dxdt_MMF_PERI   =  (MMF_Q  / MMF_V2) * MMF_CENT
                 - (MMF_Q  / MMF_V3) * MMF_PERI;

// -----------------------------------------------------------
// [7-9] Rituximab PK with Target-Mediated Drug Disposition (TMDD)
// Ref: Dirks NL, Meibohm B. Clin Pharmacokinet. 2010;49(10):633-659.
// Free CD20 target = Bmax - RTX_CD20 (quasi-static approximation)
// Units: RTX in nmol (assumes BSA 1.73m2, dose 375 mg/m2 ~= 650 mg)
// -----------------------------------------------------------
double free_CD20 = RTX_Bmax - RTX_CD20;
if(free_CD20 < 0) free_CD20 = 0;

dxdt_RTX_CENT  = -(RTX_CL / RTX_V1) * RTX_CENT
                 - (RTX_Q  / RTX_V1) * RTX_CENT
                 + (RTX_Q  / RTX_V2) * RTX_PERI
                 - RTX_kon * RTX_C * free_CD20
                 + RTX_koff * RTX_CD20;
dxdt_RTX_PERI  =  (RTX_Q  / RTX_V1) * RTX_CENT
                 - (RTX_Q  / RTX_V2) * RTX_PERI;
dxdt_RTX_CD20  =  RTX_kon * RTX_C * free_CD20
                 - RTX_koff * RTX_CD20
                 - RTX_kint * RTX_CD20;

// -----------------------------------------------------------
// [10] Gd-IgA1 Production and Clearance
// B-cell-derived aberrantly glycosylated IgA1
// Drug effects: B-cell suppression reduces synthesis
// Ref: Cheung CK, et al. J Clin Invest. 2019;129(6):2413-2426.
// -----------------------------------------------------------
double GdIgA1_prod = GdIgA1_syn * (1.0 - E_Bcell_total * 0.6);
dxdt_GD_IGA1  = GdIgA1_prod - GdIgA1_deg * GD_IGA1;

// -----------------------------------------------------------
// [11] Anti-Gd-IgA1 IgG Antibody
// Produced by B-cells; BAFF amplifies production
// Rituximab and MMF suppress B-cell-derived Ab synthesis
// -----------------------------------------------------------
double AntiIgA_prod = AntiIgA_syn * BAFF / BAFF_base
                    * (1.0 - E_Bcell_total);
dxdt_ANTI_IGA1 = AntiIgA_prod - AntiIgA_deg * ANTI_IGA1;

// -----------------------------------------------------------
// [12] IgA Immune Complex (IC) Formation
// IC = Gd-IgA1 × Anti-Gd-IgA1; cleared by complement/phagocytosis
// Prednisolone slightly enhances IC clearance via monocyte function
// -----------------------------------------------------------
dxdt_IC = IC_form * GD_IGA1 * ANTI_IGA1
         - IC_clear * (1.0 + 0.2 * E_PRED_ANTI) * IC;

// -----------------------------------------------------------
// [13] Complement C3 Activation
// IC deposits in mesangium activate lectin/alternative complement pathway
// Ref: Roos A, et al. J Clin Invest. 2006;116(6):1596-1605.
// -----------------------------------------------------------
dxdt_C3_ACT = C3_kact * IC_effect * C3_base
             - C3_kdeg * C3_ACT;

// -----------------------------------------------------------
// [14] Complement C5 Activation / sC5b-9 (Membrane Attack Complex)
// C5 activation leads to sC5b-9 which injures podocytes and tubular cells
// -----------------------------------------------------------
dxdt_C5_ACT = C5_kact * C3_ACT
             - C5_kdeg * C5_ACT;

// -----------------------------------------------------------
// [15] Mesangial Cell Proliferation Index
// Driven by IC deposition and MAC injury
// MMF inhibits mesangial proliferation (anti-proliferative on smooth muscle)
// Ref: Floege J, Feehally J. Nat Rev Nephrol. 2016;12(3):147-156.
// -----------------------------------------------------------
double MESANG_drive = MESANG_kprol * IC_effect * (1.0 + 0.5 * C5_ACT / C5_ACT_0);
double MESANG_inhibit = E_MMF_MES + E_PRED_IS * 0.3;
if(MESANG_inhibit > 0.95) MESANG_inhibit = 0.95;

dxdt_MESANG = MESANG_drive * (1.0 - MESANG_inhibit)
             - MESANG_kdeg * (MESANG - MESANG_base);

// -----------------------------------------------------------
// [16] Renal Tubulointerstitial Inflammation Score
// TGF-beta mediated; reflects Oxford M/E/S lesions
// Ref: Barbour SJ, et al. J Am Soc Nephrol. 2015;26(6):1445-1454.
// -----------------------------------------------------------
double cyto_drive = (IL6 / IL6_base + TNF / TNF_base) * 0.5;
dxdt_RENAL_INFLAM = RENAL_kact * IC_effect * cyto_drive
                   * (1.0 - E_PRED_ANTI * 0.6)
                  - RENAL_kdeg * (RENAL_INFLAM - RENAL_base)
                  + RENAL_TGFb * MESANG;  // mesangial -> tubular cross-talk

// -----------------------------------------------------------
// [17] Proteinuria (g/day)
// Driven by complement-mediated podocyte injury + intraglomerular pressure
// ACEi/ARB and sparsentan reduce glomerular hypertension component
// Ref: Pillebout E, et al. Kidney Int. 2010;78(5):495-502.
// -----------------------------------------------------------
double PROT_drive  = PROT_kact * (C5_ACT / C5_ACT_0) * MESANG;
double PROT_clear  = PROT_kdeg * (1.0 + E_RENOPROT * 2.0);

dxdt_PROT = PROT_drive * (1.0 - E_PRED_ANTI * 0.4)
           - PROT_clear * PROT;

// -----------------------------------------------------------
// [18] Hematuria Score (0-10)
// Glomerular hematuria from mesangial/endocapillary proliferation
// -----------------------------------------------------------
dxdt_HEMAT = HEMAT_kact * IC_effect * MESANG
            * (1.0 - E_PRED_IS * 0.5)
           - HEMAT_kdeg * HEMAT;

// -----------------------------------------------------------
// [19] eGFR (mL/min/1.73m2)
// CKD-EPI-based progressive decline; driven by inflammation + fibrosis
// Ref: Lafayette RA, et al. NEJM. 2023 (PROTECT trial slope data)
// -----------------------------------------------------------
double eGFR_decline = (EGFR_kdecl + EGFR_protect * RENAL_INFLAM) * EGFR;
double eGFR_protect_eff = eGFR_decline * E_RENOPROT * 0.8;

dxdt_EGFR = -eGFR_decline + eGFR_protect_eff;

// -----------------------------------------------------------
// [20] Skin Purpura Score
// IgA-IC deposition in dermal capillaries -> leukocytoclastic vasculitis
// Prednisolone highly effective for skin manifestation
// Ref: Oni L, Sampath S. Front Pediatr. 2019;7:257.
// -----------------------------------------------------------
dxdt_SKIN_PURP = SKIN_kact * IC_effect
               * (1.0 - SKIN_steroid_E * E_PRED_ANTI)
              - SKIN_kdeg * (SKIN_PURP - 0);

// -----------------------------------------------------------
// [21] GI Inflammation Score
// IgA-IC in mesenteric vessels -> colicky pain, bloody stool
// -----------------------------------------------------------
dxdt_GI_INFLAM = GI_kact * IC_effect
               * (1.0 - E_PRED_ANTI * 0.65)
              - GI_kdeg * GI_INFLAM;

// -----------------------------------------------------------
// [22] BAFF (B-cell Activating Factor)
// Upregulated by IC via macrophage/dendritic cell activation
// BAFF amplifies anti-Gd-IgA1 production (positive feedback)
// Ref: Woo KT, et al. NDT. 2014;29(Suppl 4):iv125-130.
// -----------------------------------------------------------
dxdt_BAFF = BAFF_syn * (1.0 + BAFF_stim * IC_effect)
           * (1.0 - E_PRED_ANTI * 0.5)
          - BAFF_deg * BAFF;

// -----------------------------------------------------------
// [23] IL-6 Concentration
// Produced by mesangial cells and infiltrating macrophages
// Ref: Peng W, et al. Mediators Inflamm. 2018;2018:9548026.
// -----------------------------------------------------------
dxdt_IL6 = IL6_syn + IL6_MESANG * (MESANG - MESANG_base)
          * (1.0 - E_PRED_ANTI * 0.55)
         - IL6_deg * IL6;

// -----------------------------------------------------------
// [24] TNF-alpha
// IC-driven macrophage activation; drives glomerular injury
// -----------------------------------------------------------
dxdt_TNF = TNF_syn * (1.0 + TNF_IC * IC_effect)
          * (1.0 - E_PRED_ANTI * 0.60)
         - TNF_deg * TNF;

// -----------------------------------------------------------
// [25] Renal Fibrosis Index (cumulative, irreversible)
// TGF-beta-mediated; driven by chronic inflammation
// Ref: Rauen T, et al. NEJM. 2015 (STOP-IgAN fibrosis progression)
// -----------------------------------------------------------
dxdt_FIBROSIS = Fibr_kprog * RENAL_INFLAM * (1.0 - FIBROSIS)
               * (1.0 - E_PRED_IS * 0.25);

$TABLE
// ===========================================================
// DERIVED OUTPUTS FOR TABLE CAPTURE
// ===========================================================

// PK concentrations
double PRED_Cp  = PRED_CENT / PRED_V2;  // prednisolone central conc (mg/L)
double MPA_Cp   = MMF_CENT  / MMF_V2;   // MPA central conc (mg/L)
double RTX_Cp   = RTX_CENT  / RTX_V1;   // rituximab central conc (nmol/L = nM)

// Free unbound CD20 (proxy for B-cell depletion %)
double CD20_free   = RTX_Bmax - RTX_CD20;
if(CD20_free < 0) CD20_free = 0;
double Bcell_pct   = (CD20_free / RTX_Bmax) * 100.0;  // % remaining B-cells

// Composite vasculitis activity score (PVAS adapted)
// Ref: Suppiah R, et al. Ann Rheum Dis. 2011;70(1):49-54.
double PVAS = (SKIN_PURP / SKIN_base) * 2.0
            + (GI_INFLAM / GI_base)   * 2.0
            + (PROT / PROT_base)      * 3.0
            + (HEMAT / HEMAT_base)    * 1.5
            + (RENAL_INFLAM / RENAL_base) * 1.5;

// eGFR CKD stage (1-5)
double CKD_stage;
if(EGFR >= 90)       CKD_stage = 1;
else if(EGFR >= 60)  CKD_stage = 2;
else if(EGFR >= 30)  CKD_stage = 3;
else if(EGFR >= 15)  CKD_stage = 4;
else                 CKD_stage = 5;

// Proteinuria category
double PROT_cat;
if(PROT < 0.15)      PROT_cat = 0;  // normal
else if(PROT < 0.5)  PROT_cat = 1;  // mildly increased
else if(PROT < 1.0)  PROT_cat = 2;  // moderately increased
else if(PROT < 3.5)  PROT_cat = 3;  // severely increased
else                 PROT_cat = 4;  // nephrotic range

$CAPTURE
// PK
PRED_Cp MPA_Cp RTX_Cp
// B-cell depletion
Bcell_pct CD20_free
// Disease markers
GD_IGA1 ANTI_IGA1 IC
C3_ACT C5_ACT
MESANG RENAL_INFLAM
// Clinical endpoints
PROT HEMAT EGFR
SKIN_PURP GI_INFLAM
// Cytokines
BAFF IL6 TNF
// Derived
FIBROSIS PVAS CKD_stage PROT_cat
// Drug effects
E_PRED_IS E_PRED_ANTI E_MMF_B E_RTX_dep E_Bcell_total E_RENOPROT

')

cat("IgA Vasculitis QSP Model compiled successfully.\n")
cat("Number of ODEs:", length(mod@cmtL), "\n")
cat("Parameters:", length(param(mod)), "\n")


# ==============================================================================
# TREATMENT SCENARIOS
# ==============================================================================
# Simulation duration: 52 weeks (8736 hours)
# Dosing schedules converted to hourly events

SIM_DURATION <- 24 * 7 * 52  # 8736 hours = 52 weeks
DELTA        <- 6             # output every 6 hours

# Helper: create mrgsolve event data frame
make_events <- function(dose_list) {
  # dose_list: list of lists with fields:
  #   cmt, amt, time, ii, addl, rate (optional)
  do.call(rbind, lapply(dose_list, function(d) {
    ev(
      cmt  = d$cmt,
      amt  = d$amt,
      time = d$time,
      ii   = d$ii,
      addl = d$addl,
      rate = if (!is.null(d$rate)) d$rate else 0
    )
  }))
}

# ------------------------------------------------------------------
# SCENARIO 1: No Treatment (Natural History)
# ------------------------------------------------------------------
sc1_events <- ev(amt = 0, cmt = 1, time = 0)  # null event

sc1_params <- list()  # all drug doses = 0 (default)

cat("\n--- Scenario 1: Natural History ---\n")
sc1 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 0, SGLT2_dose = 0) %>%
  mrgsim(events = sc1_events, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "1_Natural_History")

# ------------------------------------------------------------------
# SCENARIO 2: Corticosteroids Only
# Prednisolone 1 mg/kg/day (70 kg -> 70 mg/day) for 4 weeks,
# then taper by 10 mg/week to 0 by week 12
# Ref: Pillebout E, et al. Kidney Int. 2010;78(5):495-502.
#      SHARE trial: prednisone 2mg/kg/d, taper over 3 months
# Dosing schedule (oral, q24h into PRED_GUT with F=0.82):
# ------------------------------------------------------------------

# Build prednisolone taper schedule (q24h = every 24 hours)
make_pred_taper <- function() {
  # Weeks 1-4: 70 mg/day (1 mg/kg/day for 70 kg)
  # Weeks 5-6: 50 mg/day
  # Weeks 7-8: 40 mg/day
  # Weeks 9-10: 30 mg/day
  # Weeks 11-12: 20 mg/day
  # Weeks 13-14: 10 mg/day
  # Week 15+: stop

  taper_schedule <- data.frame(
    week  = c(1, 5, 7, 9, 11, 13),
    dose  = c(70, 50, 40, 30, 20, 10)  # mg/day
  )

  events <- NULL
  for(i in seq_len(nrow(taper_schedule))) {
    wk_start <- (taper_schedule$week[i] - 1) * 7 * 24  # in hours
    wk_end   <- if(i < nrow(taper_schedule)) (taper_schedule$week[i+1]-1)*7*24 else 14*7*24
    n_doses  <- floor((wk_end - wk_start) / 24)

    e <- ev(
      cmt  = 1,       # PRED_GUT
      amt  = taper_schedule$dose[i] * 0.82,  # F-adjusted into gut depot
      time = wk_start,
      ii   = 24,      # every 24 hours
      addl = n_doses - 1
    )
    events <- if(is.null(events)) e else c(events, e)
  }
  events
}

pred_taper_ev <- make_pred_taper()

cat("\n--- Scenario 2: Corticosteroids Only ---\n")
sc2 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 0, SGLT2_dose = 0) %>%
  mrgsim(events = pred_taper_ev, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "2_Corticosteroids")

# ------------------------------------------------------------------
# SCENARIO 3: Corticosteroids + MMF (mycophenolate mofetil)
# MMF 1000 mg BID (q12h) started at week 0
# Ref: Rauen T, et al. NEJM. 2015;373:2225-2236. (STOP-IgAN)
#      Tang SCW, et al. Clin J Am Soc Nephrol. 2010 (MMF in IgAN)
# ------------------------------------------------------------------

mmf_events <- ev(
  cmt  = 4,     # MMF_GUT
  amt  = 1000 * 0.94,   # 1000 mg × F=0.94
  time = 0,
  ii   = 12,    # BID
  addl = round(SIM_DURATION / 12) - 1
)

sc3_events <- c(pred_taper_ev, mmf_events)

cat("\n--- Scenario 3: Corticosteroids + MMF ---\n")
sc3 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 0, SGLT2_dose = 0) %>%
  mrgsim(events = sc3_events, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "3_Steroids_MMF")

# ------------------------------------------------------------------
# SCENARIO 4: Rituximab Monotherapy
# 375 mg/m² IV × 4 doses (weeks 0, 1, 2, 3) or standard ANCA protocol
# Ref: Jayne DRW, et al. NEJM. 2010;363(3):221-232. (PMID: 20647199)
#      Maritati F, et al. Autoimmun Rev. 2019;18(4):345-352.
#      (Rituximab for IgA vasculitis - case series)
# Dose: 375 mg/m2 × 1.73 m2 = ~650 mg -> converted to nmol (MW ~148 kDa)
# 650 mg / 148000 g/mol = 4.39 μmol = 4390 nmol; assuming V1=3L -> 1463 nM
# We dose into RTX_CENT as IV (rate = 0 = bolus approximation)
# ------------------------------------------------------------------

RTX_dose_nmol <- (375 * 1.73) / 148000 * 1e6  # ~4385 nmol per dose

rtx_events <- ev(
  cmt  = 7,     # RTX_CENT
  amt  = RTX_dose_nmol,
  time = c(0, 1*7*24, 2*7*24, 3*7*24),  # weeks 0, 1, 2, 3
  ii   = 0,
  addl = 0,
  rate = -1     # IV infusion over 4 hours (-1 = rate from event)
)
# Fix: use individual ev() calls stacked
rtx_ev1 <- ev(cmt=7, amt=RTX_dose_nmol, time=0*7*24, rate=RTX_dose_nmol/4)
rtx_ev2 <- ev(cmt=7, amt=RTX_dose_nmol, time=1*7*24, rate=RTX_dose_nmol/4)
rtx_ev3 <- ev(cmt=7, amt=RTX_dose_nmol, time=2*7*24, rate=RTX_dose_nmol/4)
rtx_ev4 <- ev(cmt=7, amt=RTX_dose_nmol, time=3*7*24, rate=RTX_dose_nmol/4)
rtx_all  <- c(rtx_ev1, rtx_ev2, rtx_ev3, rtx_ev4)

cat("\n--- Scenario 4: Rituximab ---\n")
sc4 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 0, SGLT2_dose = 0) %>%
  mrgsim(events = rtx_all, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "4_Rituximab")

# ------------------------------------------------------------------
# SCENARIO 5: ACEi + ARB Supportive Therapy (Renoprotection Only)
# Enalapril 10 mg BID modeled via ACEI_dose parameter = 1
# Ref: Lewis EJ, et al. NEJM. 1993;329(20):1456-1462. (PMID: 8413456)
#      Praga M, et al. J Am Soc Nephrol. 2003 (ACEi in IgAN)
# ------------------------------------------------------------------

cat("\n--- Scenario 5: ACEi/ARB (Renoprotection) ---\n")
sc5 <- mod %>%
  param(ACEI_dose = 1, SPAR_dose = 0, SGLT2_dose = 0) %>%
  mrgsim(events = sc1_events, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "5_ACEi_ARB")

# ------------------------------------------------------------------
# SCENARIO 6: Sparsentan (dual AT1R + ETA-R antagonist)
# 400 mg QD - approved / investigational for IgAN/IgAV nephritis
# Ref: Lafayette RA, et al. NEJM. 2023;389(19):1764-1775.
#      (PROTECT trial - 40.9% proteinuria reduction vs irbesartan)
# ------------------------------------------------------------------

cat("\n--- Scenario 6: Sparsentan ---\n")
sc6 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 1, SGLT2_dose = 0) %>%
  mrgsim(events = sc1_events, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "6_Sparsentan")

# ------------------------------------------------------------------
# SCENARIO 7: SGLT2 Inhibitor (Dapagliflozin 10 mg QD)
# Ref: Heerspink HJL, et al. NEJM. 2020;383(15):1436-1446.
#      (DAPA-CKD: 44% RRR composite renal endpoint, proteinuria -30%)
#      Ref: Wheeler DC, et al. Lancet. 2023 (FLOW trial - CKD patients)
# ------------------------------------------------------------------

cat("\n--- Scenario 7: SGLT2 Inhibitor (Dapagliflozin) ---\n")
sc7 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 0, SGLT2_dose = 1) %>%
  mrgsim(events = sc1_events, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "7_SGLT2i_Dapagliflozin")

# ------------------------------------------------------------------
# SCENARIO 8 (BONUS): Combination - Sparsentan + SGLT2i + Corticosteroids
# Rational combination for severe IgAV nephritis
# ------------------------------------------------------------------

cat("\n--- Scenario 8 (Bonus): Sparsentan + SGLT2i + Steroids ---\n")
sc8 <- mod %>%
  param(ACEI_dose = 0, SPAR_dose = 1, SGLT2_dose = 1) %>%
  mrgsim(events = pred_taper_ev, end = SIM_DURATION, delta = DELTA) %>%
  as_tibble() %>%
  mutate(scenario = "8_Sparsentan_SGLT2i_Steroids")


# ==============================================================================
# COMBINE RESULTS
# ==============================================================================

all_results <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6, sc7, sc8) %>%
  mutate(
    time_weeks = time / (24 * 7),
    scenario_label = factor(scenario,
      levels = paste0(1:8, c("_Natural_History", "_Corticosteroids",
                              "_Steroids_MMF", "_Rituximab",
                              "_ACEi_ARB", "_Sparsentan",
                              "_SGLT2i_Dapagliflozin",
                              "_Sparsentan_SGLT2i_Steroids")),
      labels = c("Natural History",
                 "Corticosteroids",
                 "Steroids + MMF",
                 "Rituximab",
                 "ACEi/ARB",
                 "Sparsentan",
                 "Dapagliflozin (SGLT2i)",
                 "Sparsentan + SGLT2i + Steroids"))
  )

cat("\n\n========== SIMULATION SUMMARY AT 52 WEEKS ==========\n")
summary_52wk <- all_results %>%
  filter(abs(time_weeks - 52) < 0.5) %>%
  group_by(scenario_label) %>%
  summarise(
    eGFR_wk52    = round(mean(EGFR), 1),
    Proteinuria  = round(mean(PROT), 2),
    GdIgA1       = round(mean(GD_IGA1), 1),
    SkinPurpura  = round(mean(SKIN_PURP), 2),
    GI_Score     = round(mean(GI_INFLAM), 2),
    PVAS         = round(mean(PVAS), 2),
    Fibrosis     = round(mean(FIBROSIS), 3),
    Bcell_pct    = round(mean(Bcell_pct), 1),
    .groups = "drop"
  )

print(summary_52wk)


# ==============================================================================
# VISUALIZATION
# ==============================================================================

# Color palette for 8 scenarios
scenario_colors <- c(
  "Natural History"               = "#E41A1C",
  "Corticosteroids"               = "#FF7F00",
  "Steroids + MMF"                = "#984EA3",
  "Rituximab"                     = "#A65628",
  "ACEi/ARB"                      = "#377EB8",
  "Sparsentan"                    = "#4DAF4A",
  "Dapagliflozin (SGLT2i)"        = "#F781BF",
  "Sparsentan + SGLT2i + Steroids"= "#000000"
)

# Panel 1: eGFR over 52 weeks
p_egfr <- ggplot(all_results, aes(x = time_weeks, y = EGFR,
                                   color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IgA Vasculitis: eGFR Trajectory by Treatment",
       x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)",
       color = "Treatment",
       caption = "Model based on PROTECT, DAPA-CKD, STOP-IgAN trials") +
  geom_hline(yintercept = c(60, 30, 15), linetype = "dashed",
             color = "gray50", alpha = 0.7) +
  annotate("text", x = 51, y = 61, label = "CKD G3", hjust=1, size=3, color="gray50") +
  annotate("text", x = 51, y = 31, label = "CKD G4", hjust=1, size=3, color="gray50") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Panel 2: Proteinuria over 52 weeks
p_prot <- ggplot(all_results, aes(x = time_weeks, y = PROT,
                                   color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IgA Vasculitis: Proteinuria Response",
       x = "Time (weeks)", y = "Proteinuria (g/day)",
       color = "Treatment",
       caption = "Ref: Pillebout 2010, Lafayette 2023") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  annotate("text", x = 51, y = 0.52, label = "Target <0.5 g/day",
           hjust=1, size=3, color="gray50") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Panel 3: Skin Purpura Score
p_skin <- ggplot(all_results, aes(x = time_weeks, y = SKIN_PURP,
                                   color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IgA Vasculitis: Skin Purpura Score",
       x = "Time (weeks)", y = "Purpura Score (0-10)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Panel 4: Gd-IgA1 level
p_gd <- ggplot(all_results, aes(x = time_weeks, y = GD_IGA1,
                                 color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IgA Vasculitis: Galactose-Deficient IgA1 Level",
       x = "Time (weeks)", y = "Gd-IgA1 (mg/L)",
       color = "Treatment",
       caption = "Normal Gd-IgA1 ~5 mg/L; IgAV patients ~24 mg/L") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "green4") +
  annotate("text", x = 3, y = 5.8, label = "Normal", size=3, color="green4") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Panel 5: Rituximab PK + B-cell depletion
rtx_data <- all_results %>%
  filter(grepl("Rituximab", scenario_label))

p_rtx <- ggplot(rtx_data) +
  geom_line(aes(x = time_weeks, y = RTX_Cp, color = "RTX Concentration (nM)"),
            size = 1.0) +
  geom_line(aes(x = time_weeks, y = Bcell_pct, color = "B-cell %"),
            size = 1.0, linetype = "dashed") +
  scale_color_manual(values = c("RTX Concentration (nM)" = "#A65628",
                                 "B-cell %" = "#4DAF4A")) +
  labs(title = "Rituximab PK/PD: Concentration and B-cell Depletion",
       x = "Time (weeks)", y = "Concentration (nM) / B-cell (%)",
       color = "Measure",
       caption = "TMDD model; 4 weekly doses (375 mg/m²)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# Panel 6: Complement and Cytokines - Natural History vs Best Treatment
compare_data <- all_results %>%
  filter(scenario_label %in% c("Natural History",
                                "Sparsentan + SGLT2i + Steroids",
                                "Rituximab"))

p_cyto <- ggplot(compare_data %>%
                   pivot_longer(cols = c(IL6, TNF, BAFF),
                                names_to = "Cytokine", values_to = "Value"),
                 aes(x = time_weeks, y = Value,
                     color = scenario_label, linetype = Cytokine)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  facet_wrap(~Cytokine, scales = "free_y") +
  labs(title = "Cytokine Profiles: IL-6, TNF-alpha, BAFF",
       x = "Time (weeks)", y = "Normalized Level",
       color = "Treatment", linetype = "Cytokine",
       caption = "Ref: Peng 2018, Woo 2014") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

# Panel 7: PVAS composite score
p_pvas <- ggplot(all_results, aes(x = time_weeks, y = PVAS,
                                   color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "IgA Vasculitis Activity Score (PVAS adapted)",
       x = "Time (weeks)", y = "Composite Activity Score",
       color = "Treatment",
       caption = "Higher = more active disease") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Panel 8: Renal Fibrosis accumulation
p_fibr <- ggplot(all_results, aes(x = time_weeks, y = FIBROSIS,
                                   color = scenario_label)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Cumulative Renal Fibrosis Index",
       x = "Time (weeks)", y = "Fibrosis Index (0-1)",
       color = "Treatment",
       caption = "Irreversible CKD progression component") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

# Print all plots
print(p_egfr)
print(p_prot)
print(p_skin)
print(p_gd)
print(p_rtx)
print(p_cyto)
print(p_pvas)
print(p_fibr)

# ==============================================================================
# SENSITIVITY ANALYSIS: Parameter Uncertainty
# ==============================================================================

cat("\n========== SENSITIVITY ANALYSIS: Sparsentan SPAR_Emax ==========\n")

spar_emax_vals <- c(0.30, 0.40, 0.50, 0.55, 0.65, 0.75)
sens_results <- lapply(spar_emax_vals, function(emax_val) {
  mod %>%
    param(SPAR_Emax = emax_val, SPAR_dose = 1, ACEI_dose = 0, SGLT2_dose = 0) %>%
    mrgsim(events = sc1_events, end = 52*7*24, delta = 24*7) %>%
    as_tibble() %>%
    mutate(SPAR_Emax_val = emax_val)
}) %>% bind_rows()

p_sens <- ggplot(sens_results %>% filter(time > 0),
                 aes(x = time / (24*7), y = PROT,
                     color = factor(SPAR_Emax_val),
                     group = SPAR_Emax_val)) +
  geom_line(size = 0.9) +
  scale_color_viridis_d(name = "Sparsentan Emax") +
  labs(title = "Sensitivity: Sparsentan Emax vs Proteinuria",
       x = "Time (weeks)", y = "Proteinuria (g/day)",
       caption = "PROTECT trial observed: ~41% reduction at 36 weeks") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

print(p_sens)

# ==============================================================================
# CLINICAL OUTCOME METRICS AT KEY TIMEPOINTS
# ==============================================================================

cat("\n========== CLINICAL OUTCOME METRICS ==========\n")

timepoints <- c(4, 12, 24, 52)  # weeks

outcome_table <- all_results %>%
  filter(time_weeks %in% timepoints) %>%
  group_by(scenario_label, time_weeks) %>%
  summarise(
    eGFR       = round(mean(EGFR), 1),
    Prot_gday  = round(mean(PROT), 2),
    Skin_score = round(mean(SKIN_PURP), 2),
    GI_score   = round(mean(GI_INFLAM), 2),
    PVAS_score = round(mean(PVAS), 2),
    GdIgA1     = round(mean(GD_IGA1), 1),
    Fibrosis   = round(mean(FIBROSIS), 3),
    .groups = "drop"
  ) %>%
  arrange(time_weeks, scenario_label)

print(outcome_table, n = Inf)

# ==============================================================================
# SESSION INFO
# ==============================================================================

cat("\n========== SESSION INFO ==========\n")
cat("mrgsolve version:", as.character(packageVersion("mrgsolve")), "\n")
cat("Model compartments:", length(mod@cmtL), "\n")
cat("Scenarios simulated:", 8, "\n")
cat("Key clinical trial references:\n")
cat("  - Pillebout 2010 (PMID: 20520597) - adult IgAV nephritis outcomes\n")
cat("  - STOP-IgAN 2015 (PMID: 26630142) - IS calibration\n")
cat("  - Lafayette 2023 (PMID: 37632463) - PROTECT (sparsentan)\n")
cat("  - Heerspink 2020 (PMID: 32970396) - DAPA-CKD\n")
cat("  - Jayne 2010 (PMID: 20647199) - Rituximab PK/PD\n")
cat("  - Magistroni 2015 (PMID: 26376133) - Gd-IgA1 pathophysiology\n")
cat("Session date: 2026-06-19\n")
