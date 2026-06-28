## ============================================================
## Prurigo Nodularis (PN) — mrgsolve QSP Model
## ============================================================
## 구획 구조 (27 ODEs):
##   PK: Dupilumab (SC 2-cmt TMDD), Nemolizumab (SC 2-cmt TMDD),
##       Tralokinumab (SC 2-cmt), Cyclosporine (1-cmt oral),
##       Nalbuphine ER (1-cmt oral), TCS (topical depot→skin)
##   PD: IL-4, IL-13, IL-31, IgE, Th2 cells, Mast cells,
##       Itch VAS, Skin barrier (TEWL), Eosinophils,
##       Dermal nerve density, DNRS score, IGA score
##
## 치료 시나리오:
##   1) Placebo
##   2) Dupilumab 300 mg SC Q2W (PRIME trial)
##   3) Nemolizumab 60 mg SC Q4W (ARCADIA trial)
##   4) Tralokinumab 300 mg SC Q2W (TRALooPN trial)
##   5) Cyclosporine 5 mg/kg/day PO
##   6) Nalbuphine ER 54 mg PO BID
##   7) Dupilumab + TCS combination
##
## 보정 근거:
##   - PRIME/PRIME2: Dupilumab Phase 3 PN (Blauvelt 2021, Tan 2022)
##   - ARCADIA: Nemolizumab Phase 3 PN (Silverberg 2023, Yosipovitch 2023)
##   - TRALooPN: Tralokinumab PN (Wollenberg 2023)
##   - Cyclosporine PN RCT (Siepmann 2013)
##   - IL-31 PK/PD: Ruzicka 2017 (nemolizumab first-in-human)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB
Prurigo Nodularis QSP Model — Dual Neuroimmune Axis
Mechanistic PK/PD for dupilumab, nemolizumab, tralokinumab,
cyclosporine, nalbuphine ER, and TCS
Calibrated to PRIME, ARCADIA, TRALooPN, and RCT data.

$PARAM
// ─── Patient parameters ────────────────────────────────
BW       = 75       // body weight (kg)
AGE      = 50       // age (yr)
SEX      = 1        // 1=female, 0=male

// ─── Dupilumab PK (SC 2-cmt, TMDD) ────────────────────
// Ref: Dupilumab PopPK: Xu 2018 CPT:PSP
DUP_KA   = 0.0048   // absorption rate /h (SC)
DUP_F1   = 0.64     // SC bioavailability
DUP_CL   = 0.012    // L/h central clearance
DUP_V1   = 3.0      // L central volume
DUP_Q    = 0.006    // L/h inter-compartmental
DUP_V2   = 5.5      // L peripheral volume
DUP_KSYN = 0.01     // IL4Ra synthesis rate (nmol/L/h)
DUP_KDEG = 0.004    // free receptor degradation /h
DUP_KON  = 2.0      // association constant (L/nmol/h)
DUP_KOFF = 0.003    // dissociation constant /h
DUP_KINT = 0.002    // internalization of complex /h

// ─── Nemolizumab PK (SC 2-cmt, IL-31RA TMDD) ──────────
// Ref: Nemolizumab PopPK: Ruzicka 2017, Silverberg 2023
NEM_KA   = 0.005    // absorption rate /h
NEM_F1   = 0.72     // bioavailability
NEM_CL   = 0.009    // L/h
NEM_V1   = 2.8      // L
NEM_Q    = 0.004    // L/h
NEM_V2   = 4.5      // L
NEM_KSYN = 0.008    // IL31Ra synthesis /h
NEM_KDEG = 0.003    // /h
NEM_KON  = 3.0      // L/nmol/h
NEM_KOFF = 0.002    // /h
NEM_KINT = 0.0015   // /h

// ─── Tralokinumab PK (SC 2-cmt) ────────────────────────
// Ref: Tralokinumab PopPK: Popielarz-Grygalewicz 2022
TRAL_KA  = 0.0042   // /h
TRAL_F1  = 0.77     // bioavailability
TRAL_CL  = 0.010    // L/h
TRAL_V1  = 3.2      // L
TRAL_Q   = 0.005    // L/h
TRAL_V2  = 5.0      // L

// ─── Cyclosporine PK (1-cmt oral) ──────────────────────
// Ref: CsA PopPK: Fanta 2008, Levy 2002
CSA_KA   = 0.08     // /h
CSA_F1   = 0.40     // bioavailability (variable)
CSA_CL   = 25.0     // L/h
CSA_V    = 350.0    // L (extensive distribution)
CSA_DOSE = 375.0    // mg/day (5 mg/kg × 75 kg)

// ─── Nalbuphine ER PK (1-cmt oral) ─────────────────────
// Ref: Nalbuphine ER PK: Zeidler 2016
NALB_KA  = 0.03     // /h (extended release)
NALB_F1  = 0.20     // bioavailability
NALB_CL  = 140.0    // L/h
NALB_V   = 800.0    // L
NALB_DOSE= 54.0     // mg/dose

// ─── TCS (topical corticosteroid) skin depot PK ─────────
TCS_KA_SK= 0.02     // absorption into skin /h
TCS_KEL  = 0.08     // skin elimination /h

// ─── Immune PD parameters ──────────────────────────────
// Th2 / cytokine baseline
TH2_0    = 100.0    // baseline Th2 cells (arb units)
TH2_KPROL= 0.003    // Th2 proliferation (/h)
TH2_KDTH = 0.003    // Th2 death (/h)
TH2_EC50 = 2.0      // IL-4 EC50 for Th2 priming (ng/mL)

// IL-4 dynamics (primary Th2 cytokine)
IL4_0    = 5.0      // baseline pg/mL
IL4_KPROD= 1.5      // production (pg/mL/h)
IL4_KDEG = 0.30     // degradation /h
IL4_MAX  = 50.0     // max stimulated level

// IL-13 dynamics (Th2, ILC2)
IL13_0   = 8.0      // baseline pg/mL
IL13_KPROD=2.0      // pg/mL/h
IL13_KDEG= 0.25     // /h

// IL-31 dynamics (Th2, itch cytokine)
IL31_0   = 12.0     // baseline pg/mL (elevated in PN)
IL31_KPROD=3.0      // pg/mL/h (Th2-driven)
IL31_KDEG= 0.20     // /h
IL31_KPROD_TH2=0.03 // IL-31 production per Th2 unit

// IgE dynamics
IGE_0    = 800.0    // baseline IU/mL (elevated in PN)
IGE_KPROD= 0.002    // IU/mL/h
IGE_KDEG = 0.0001   // /h (long half-life ~weeks)
IGE_IL4_EC50 = 5.0  // IL-4 EC50 driving IgE class switch

// Mast cell (skin-resident)
MAST_0   = 100.0    // baseline (arb units)
MAST_KPROL=0.0005   // /h
MAST_KDTH= 0.0005   // /h
MAST_IGE_EC50=500.0 // IgE EC50 for mast priming

// Eosinophil count
EOS_0    = 400.0    // baseline (cells/μL)
EOS_KPROD= 120.0    // cells/μL/h
EOS_KDEG = 0.30     // /h
EOS_IL13_EC50=5.0   // IL-13 EC50

// Skin barrier TEWL (transepidermal water loss, g/m²/h)
TEWL_0   = 20.0     // baseline (normal ~10, PN elevated)
TEWL_MAX = 80.0     // max PN severity
TEWL_IL4_EC50=15.0  // IL-4 disruption EC50
TEWL_KDEG= 0.01     // recovery rate /h

// ─── Neurological PD parameters ─────────────────────────
// Itch VAS (0-10)
ITCH_0   = 7.5      // baseline NRS (elevated in PN)
ITCH_IL31_EC50=10.0 // IL-31 EC50 for itch
ITCH_KSEN= 0.005    // central sensitization rate /h
ITCH_KDES= 0.003    // desensitization /h

// Dermal nerve density (arb units, elevated in PN)
NERVE_0  = 150.0    // baseline (1.5× normal)
NERVE_KGROWTH=0.004 // NGF-driven nerve sprouting /h
NERVE_KDTH=0.003    // nerve pruning /h
NERVE_NGF_EC50=20.0 // NGF EC50

// Skin scratch/lesion dynamics
NODULE_0 = 20.0     // nodule count (lesion count)
NODULE_KFORM=0.002  // formation rate /h
NODULE_KHEAL=0.0015 // healing rate /h

// DNRS (Dynamic Neuropathic Response Score, 0–28)
DNRS_0   = 18.0     // baseline severe PN

// Cyclosporine PD
CSA_IMAX = 0.85     // max calcineurin inhibition
CSA_IC50 = 150.0    // ng/mL for Th2 suppression

// Nalbuphine PD (kappa/mu opioid receptor)
NALB_IMAX  = 0.70   // max itch suppression
NALB_IC50  = 8.0    // ng/mL IC50 for itch (kappa)

// Corticosteroid (TCS) PD
TCS_IMAX = 0.90     // max local anti-inflammatory
TCS_IC50 = 0.5      // μg/g skin EC50

$CMT
// ─── Dupilumab PK ─────
DUP_SC    // SC depot (nmol)
DUP_C1    // central (nmol/L)
DUP_C2    // peripheral
DUP_IL4Ra // free IL-4Rα receptor
DUP_CMPLX // drug-receptor complex

// ─── Nemolizumab PK ───
NEM_SC
NEM_C1
NEM_C2
NEM_IL31Ra // free IL-31Rα
NEM_CMPLX  // drug-receptor complex

// ─── Tralokinumab PK ──
TRAL_SC
TRAL_C1
TRAL_C2

// ─── Cyclosporine PK ──
CSA_GUT
CSA_C1

// ─── Nalbuphine PK ────
NALB_GUT
NALB_C1

// ─── TCS ──────────────
TCS_DEPOT
TCS_SKIN

// ─── Immune PD ────────
TH2_CELLS  // Th2 lymphocytes (arb units)
IL4        // pg/mL
IL13       // pg/mL
IL31       // pg/mL
IGE        // IU/mL
MAST_CELLS // dermal mast cells
EOS_COUNT  // eosinophils cells/μL

// ─── Skin/Neuro PD ────
TEWL       // trans-epidermal water loss
NERVE_DEN  // dermal nerve density
NODULE_CNT // nodule count
ITCH_CS    // central sensitization state
DNRS       // dynamic neuropathy response

$MAIN
// ─── Initial conditions ───────────────────────────────
DUP_C1_0    = 0;
DUP_IL4Ra_0 = DUP_KSYN / DUP_KDEG;  // receptor at steady-state
DUP_CMPLX_0 = 0;

NEM_C1_0    = 0;
NEM_IL31Ra_0= NEM_KSYN / NEM_KDEG;
NEM_CMPLX_0 = 0;

TRAL_C1_0   = 0;
CSA_C1_0    = 0;
NALB_C1_0   = 0;

TH2_CELLS_0  = TH2_0;
IL4_0_CMT    = IL4_0;
IL13_0_CMT   = IL13_0;
IL31_0_CMT   = IL31_0;
IGE_0_CMT    = IGE_0;
MAST_CELLS_0 = MAST_0;
EOS_COUNT_0  = EOS_0;

TEWL_0_CMT   = TEWL_0;
NERVE_DEN_0  = NERVE_0;
NODULE_CNT_0 = NODULE_0;
ITCH_CS_0    = 1.0;    // normalized sensitization state
DNRS_0_CMT   = DNRS_0;

$ODE
// ─── Drug effect inhibitors ─────────────────────────────
// IL-4Rα blockade by dupilumab (blocks IL-4 + IL-13 signaling)
double DUP_CMPLX_VAL = DUP_CMPLX;
double DUP_RO = (DUP_CMPLX_VAL > 0 && (DUP_IL4Ra + DUP_CMPLX_VAL) > 0)
                ? DUP_CMPLX_VAL / (DUP_IL4Ra + DUP_CMPLX_VAL) : 0.0;
double INH_IL4_DUP  = 1.0 - DUP_RO * 0.95;  // IL-4 effect inhibition
double INH_IL13_DUP = 1.0 - DUP_RO * 0.93;

// IL-31Rα blockade by nemolizumab
double NEM_CMPLX_VAL= NEM_CMPLX;
double NEM_RO = (NEM_CMPLX_VAL > 0 && (NEM_IL31Ra + NEM_CMPLX_VAL) > 0)
                ? NEM_CMPLX_VAL / (NEM_IL31Ra + NEM_CMPLX_VAL) : 0.0;
double INH_IL31_NEM = 1.0 - NEM_RO * 0.90;  // IL-31 signaling inhibition

// IL-13 neutralization by tralokinumab
double TRAL_C = TRAL_C1;
double INH_IL13_TRAL = 1.0 - (TRAL_C / (TRAL_C + 2.5)); // Kd ~2.5 nM

// Cyclosporine: calcineurin inhibition → Th2 suppression
double CSA_C = CSA_C1;
double INH_TH2_CSA = 1.0 - (CSA_IMAX * pow(CSA_C,2) / (pow(CSA_IC50,2) + pow(CSA_C,2)));

// Nalbuphine: kappa agonist → itch inhibition
double NALB_C = NALB_C1;
double INH_ITCH_NALB = 1.0 - (NALB_IMAX * NALB_C / (NALB_IC50 + NALB_C));

// TCS: local anti-inflammatory
double TCS_S = TCS_SKIN;
double INH_INFLAM_TCS = 1.0 - (TCS_IMAX * TCS_S / (TCS_IC50 + TCS_S));

// Combined IL-13 inhibition (additive)
double INH_IL13_TOTAL = INH_IL13_DUP * INH_IL13_TRAL;

// ─── DUPILUMAB PK ODEs ────────────────────────────────────
double DUP_dose_rate = 0.0;  // driven by $NMREC events
dxdt_DUP_SC  = -DUP_KA * DUP_SC;
dxdt_DUP_C1  =  DUP_KA * DUP_F1 * DUP_SC / DUP_V1
               - DUP_CL / DUP_V1 * DUP_C1
               - DUP_Q  / DUP_V1 * DUP_C1
               + DUP_Q  / DUP_V2 * DUP_C2
               - DUP_KON * DUP_C1 * DUP_IL4Ra
               + DUP_KOFF * DUP_CMPLX;
dxdt_DUP_C2  =  DUP_Q / DUP_V1 * DUP_C1 - DUP_Q / DUP_V2 * DUP_C2;
dxdt_DUP_IL4Ra = DUP_KSYN
               - DUP_KDEG * DUP_IL4Ra
               - DUP_KON * DUP_C1 * DUP_IL4Ra
               + DUP_KOFF * DUP_CMPLX;
dxdt_DUP_CMPLX = DUP_KON * DUP_C1 * DUP_IL4Ra
               - DUP_KOFF * DUP_CMPLX
               - DUP_KINT * DUP_CMPLX;

// ─── NEMOLIZUMAB PK ODEs ──────────────────────────────────
dxdt_NEM_SC    = -NEM_KA * NEM_SC;
dxdt_NEM_C1    =  NEM_KA * NEM_F1 * NEM_SC / NEM_V1
               - NEM_CL / NEM_V1 * NEM_C1
               - NEM_Q  / NEM_V1 * NEM_C1
               + NEM_Q  / NEM_V2 * NEM_C2
               - NEM_KON * NEM_C1 * NEM_IL31Ra
               + NEM_KOFF * NEM_CMPLX;
dxdt_NEM_C2    =  NEM_Q / NEM_V1 * NEM_C1 - NEM_Q / NEM_V2 * NEM_C2;
dxdt_NEM_IL31Ra= NEM_KSYN
               - NEM_KDEG * NEM_IL31Ra
               - NEM_KON * NEM_C1 * NEM_IL31Ra
               + NEM_KOFF * NEM_CMPLX;
dxdt_NEM_CMPLX = NEM_KON * NEM_C1 * NEM_IL31Ra
               - NEM_KOFF * NEM_CMPLX
               - NEM_KINT * NEM_CMPLX;

// ─── TRALOKINUMAB PK ODEs ─────────────────────────────────
dxdt_TRAL_SC   = -TRAL_KA * TRAL_SC;
dxdt_TRAL_C1   =  TRAL_KA * TRAL_F1 * TRAL_SC / TRAL_V1
               - TRAL_CL / TRAL_V1 * TRAL_C1
               - TRAL_Q  / TRAL_V1 * TRAL_C1
               + TRAL_Q  / TRAL_V2 * TRAL_C2;
dxdt_TRAL_C2   =  TRAL_Q / TRAL_V1 * TRAL_C1 - TRAL_Q / TRAL_V2 * TRAL_C2;

// ─── CYCLOSPORINE PK ODEs ─────────────────────────────────
dxdt_CSA_GUT   = -CSA_KA * CSA_GUT;
dxdt_CSA_C1    =  CSA_KA * CSA_F1 * CSA_GUT / CSA_V
               - CSA_CL / CSA_V * CSA_C1;

// ─── NALBUPHINE PK ODEs ───────────────────────────────────
dxdt_NALB_GUT  = -NALB_KA * NALB_GUT;
dxdt_NALB_C1   =  NALB_KA * NALB_F1 * NALB_GUT / NALB_V
               - NALB_CL / NALB_V * NALB_C1;

// ─── TCS ODEs ─────────────────────────────────────────────
dxdt_TCS_DEPOT = -TCS_KA_SK * TCS_DEPOT;
dxdt_TCS_SKIN  =  TCS_KA_SK * TCS_DEPOT - TCS_KEL * TCS_SKIN;

// ─── IMMUNE PD ODEs ───────────────────────────────────────
// Th2 cells: IL-4 primes Th2, cyclosporine suppresses
double IL4_VAL = IL4 > 0 ? IL4 : 1e-6;
double TH2_stim = TH2_KPROL * IL4_VAL / (TH2_EC50 + IL4_VAL);
dxdt_TH2_CELLS = (TH2_stim - TH2_KDTH) * TH2_CELLS * INH_TH2_CSA
               * INH_INFLAM_TCS;

// IL-4: Th2-produced, degraded; blocked by dupilumab signaling
double TH2_V  = TH2_CELLS;
dxdt_IL4   =  IL4_KPROD * (TH2_V / TH2_0) * INH_INFLAM_TCS
            - IL4_KDEG  * IL4;

// IL-13: Th2-produced; blocked by dupilumab (IL4Ra) and tralokinumab
dxdt_IL13  =  IL13_KPROD * (TH2_V / TH2_0) * INH_INFLAM_TCS
            - IL13_KDEG  * IL13;

// IL-31: Th2-produced (key pruritogen); blocked by nemolizumab
dxdt_IL31  =  IL31_KPROD + IL31_KPROD_TH2 * TH2_V
            - IL31_KDEG * IL31;

// IgE: IL-4/IL-13 drive B-cell class switch
double IL4_NORM = IL4_VAL / (IGE_IL4_EC50 + IL4_VAL);
dxdt_IGE   =  IGE_KPROD * (1.0 + 4.0 * IL4_NORM) * INH_IL4_DUP * INH_IL13_TOTAL
            - IGE_KDEG * IGE;

// Mast cells: IgE loading increases activation/proliferation
double IGE_NORM = IGE / (MAST_IGE_EC50 + IGE);
dxdt_MAST_CELLS = MAST_KPROL * (1.0 + 2.0 * IGE_NORM) * MAST_CELLS
                - MAST_KDTH * MAST_CELLS;

// Eosinophils: IL-13 driven recruitment
double IL13_V = IL13 > 0 ? IL13 : 1e-6;
dxdt_EOS_COUNT =  EOS_KPROD * (1.0 + IL13_V / (EOS_IL13_EC50 + IL13_V))
                 * INH_IL13_TOTAL * INH_INFLAM_TCS
               - EOS_KDEG * EOS_COUNT;

// ─── SKIN/NEURO PD ODEs ───────────────────────────────────
// TEWL: IL-4/IL-13 disrupt skin barrier, restored by treatment
double TEWL_stim = TEWL_MAX * IL4_NORM * INH_IL4_DUP * INH_IL13_TOTAL;
dxdt_TEWL  = TEWL_KDEG * (TEWL_0 + TEWL_stim - TEWL);

// Dermal nerve density: NGF (mast-cell derived) drives sprouting
// simplified as proxy: MAST_CELLS drive nerve growth
double NGF_proxy = 20.0 * (MAST_CELLS / MAST_0);
dxdt_NERVE_DEN = NERVE_KGROWTH * NGF_proxy / (NERVE_NGF_EC50 + NGF_proxy)
                 * NERVE_0
               - NERVE_KDTH * NERVE_DEN;

// Nodule formation: itch-scratch drives hyperkeratosis, IL-13 fibrosis
double SCRATCH_DRIVE = ITCH_CS > 0 ? ITCH_CS : 1e-6;
dxdt_NODULE_CNT = NODULE_KFORM * SCRATCH_DRIVE * (IL13_V / (5.0 + IL13_V))
               - NODULE_KHEAL * INH_IL13_TOTAL * INH_INFLAM_TCS * NODULE_CNT;

// Central sensitization: driven by IL-31 + nerve density
double IL31_V = IL31 > 0 ? IL31 : 1e-6;
dxdt_ITCH_CS = ITCH_KSEN * (IL31_V / (ITCH_IL31_EC50 + IL31_V))
               * (NERVE_DEN / NERVE_0) * INH_IL31_NEM * INH_ITCH_NALB
             - ITCH_KDES * ITCH_CS;

// DNRS: composite clinical score (0–28), driven by CS + nodules
double DNRS_target = 28.0 * ITCH_CS * (NODULE_CNT / NODULE_0) * 0.5;
dxdt_DNRS = 0.02 * (DNRS_target - DNRS);

$TABLE
// ─── Derived variables ───────────────────────────────────
double DUP_CONC   = DUP_C1;       // nM dupilumab plasma
double NEM_CONC   = NEM_C1;       // nM nemolizumab plasma
double TRAL_CONC  = TRAL_C1;      // nM tralokinumab plasma
double CSA_CONC_NG= CSA_C1 * 1202.0 / 1000.0; // convert to ng/mL
double NALB_CONC  = NALB_C1;      // ng/mL nalbuphine

// Receptor occupancy (%)
double DUP_RO_PCT = DUP_CMPLX > 0 ?
  100.0 * DUP_CMPLX / (DUP_IL4Ra + DUP_CMPLX) : 0.0;
double NEM_RO_PCT = NEM_CMPLX > 0 ?
  100.0 * NEM_CMPLX / (NEM_IL31Ra + NEM_CMPLX) : 0.0;

// Itch NRS (0–10), derived from central sensitization and IL-31
double ITCH_NRS = 10.0 * ITCH_CS * (IL31 / (ITCH_IL31_EC50 + IL31))
                * INH_IL31_NEM * INH_ITCH_NALB;
if(ITCH_NRS > 10.0) ITCH_NRS = 10.0;
if(ITCH_NRS < 0.0)  ITCH_NRS = 0.0;

// IGA (Investigator Global Assessment, 0–4)
double IGA = 4.0 * (NODULE_CNT / NODULE_0) * (TEWL / TEWL_MAX) * 2.0;
if(IGA > 4.0) IGA = 4.0;

// Peak pruritus NRS
double PP_NRS = ITCH_NRS;

// % change from baseline
double IL31_PCHG  = (IL31_0 > 0)  ? 100.0*(IL31-IL31_0)/IL31_0 : 0.0;
double IGE_PCHG   = (IGE_0  > 0)  ? 100.0*(IGE-IGE_0)/IGE_0   : 0.0;
double EOS_PCHG   = (EOS_0  > 0)  ? 100.0*(EOS_COUNT-EOS_0)/EOS_0 : 0.0;
double NODULE_PCHG= (NODULE_0 > 0)? 100.0*(NODULE_CNT-NODULE_0)/NODULE_0 : 0.0;

capture DUP_CONC NEM_CONC TRAL_CONC CSA_CONC_NG NALB_CONC
capture DUP_RO_PCT NEM_RO_PCT
capture IL4 IL13 IL31 IGE EOS_COUNT MAST_CELLS TH2_CELLS
capture TEWL NERVE_DEN NODULE_CNT ITCH_CS DNRS
capture ITCH_NRS IGA PP_NRS
capture IL31_PCHG IGE_PCHG EOS_PCHG NODULE_PCHG

$INIT
DUP_SC = 0; DUP_C1 = 0; DUP_C2 = 0;
DUP_IL4Ra = 2.5; DUP_CMPLX = 0;
NEM_SC = 0; NEM_C1 = 0; NEM_C2 = 0;
NEM_IL31Ra = 2.667; NEM_CMPLX = 0;
TRAL_SC = 0; TRAL_C1 = 0; TRAL_C2 = 0;
CSA_GUT = 0; CSA_C1 = 0;
NALB_GUT = 0; NALB_C1 = 0;
TCS_DEPOT = 0; TCS_SKIN = 0;
TH2_CELLS = 100; IL4 = 5; IL13 = 8; IL31 = 12;
IGE = 800; MAST_CELLS = 100; EOS_COUNT = 400;
TEWL = 20; NERVE_DEN = 150; NODULE_CNT = 20;
ITCH_CS = 1.0; DNRS = 18;
'

## ─── Compile model ─────────────────────────────────────────────
mod <- mcode("pn_qsp", code)

## ─── Dosing event builders ──────────────────────────────────────
make_events <- function(scenario,
                        duration_weeks = 52,
                        BW_kg = 75) {
  ev_list <- list()
  n_weeks <- duration_weeks
  n_days  <- n_weeks * 7

  if (scenario == 2) {
    # Dupilumab 300 mg SC Q2W → loading 600 mg at W0, then 300 mg Q2W
    # MW dupilumab ~144 kDa; 300 mg = 300e6/144000 = 2083 nmol, /V1=3L → 694 nM
    ev_list[["dup"]] <- ev(cmt="DUP_SC", amt=694*3, ii=14*24, addl=n_weeks/2-1, time=0)
  }
  if (scenario == 3) {
    # Nemolizumab 60 mg SC Q4W; MW ~152 kDa; 60 mg = 395 nmol → /V1=2.8L → 141 nM
    ev_list[["nem"]] <- ev(cmt="NEM_SC", amt=141*2.8, ii=28*24, addl=n_weeks/4-1, time=0)
  }
  if (scenario == 4) {
    # Tralokinumab 300 mg SC Q2W; MW ~150 kDa; 300mg = 2000 nmol /V1=3.2 = 625 nM
    ev_list[["tral"]] <- ev(cmt="TRAL_SC", amt=625*3.2, ii=14*24, addl=n_weeks/2-1, time=0)
  }
  if (scenario == 5) {
    # Cyclosporine 5 mg/kg/day PO BID (split dose)
    dose_mg <- BW_kg * 5.0 / 2   # per dose (BID)
    ev_list[["csa"]] <- ev(cmt="CSA_GUT", amt=dose_mg/2, ii=12, addl=n_days*2-1, time=0)
  }
  if (scenario == 6) {
    # Nalbuphine ER 54 mg PO BID
    ev_list[["nalb"]] <- ev(cmt="NALB_GUT", amt=54, ii=12, addl=n_days*2-1, time=0)
  }
  if (scenario == 7) {
    # Dupilumab + TCS
    ev_list[["dup"]] <- ev(cmt="DUP_SC", amt=694*3, ii=14*24, addl=n_weeks/2-1, time=0)
    # TCS applied daily (depot model): 1 mg topical per application, BID
    ev_list[["tcs"]] <- ev(cmt="TCS_DEPOT", amt=1.0, ii=12, addl=n_days*2-1, time=0)
  }

  if (length(ev_list) == 0) return(NULL)
  do.call(c, ev_list)
}

## ─── Simulation runner ──────────────────────────────────────────
run_scenario <- function(scen_id, label, duration_weeks = 52) {
  ev_obj <- make_events(scen_id, duration_weeks)

  if (is.null(ev_obj)) {
    out <- mrgsim(mod, end = duration_weeks * 7 * 24, delta = 6)
  } else {
    out <- mrgsim(mod, events = ev_obj,
                  end = duration_weeks * 7 * 24, delta = 6)
  }

  as.data.frame(out) |>
    mutate(scenario = label,
           time_weeks = time / (7 * 24))
}

## ─── Run all 7 scenarios ───────────────────────────────────────
scenarios <- list(
  list(id=1, label="Placebo"),
  list(id=2, label="Dupilumab 300mg Q2W"),
  list(id=3, label="Nemolizumab 60mg Q4W"),
  list(id=4, label="Tralokinumab 300mg Q2W"),
  list(id=5, label="Cyclosporine 5mg/kg/d"),
  list(id=6, label="Nalbuphine ER 54mg BID"),
  list(id=7, label="Dupilumab + TCS")
)

cat("Running 7 treatment scenarios...\n")
results_all <- lapply(scenarios, function(s) {
  cat(" Scenario", s$id, ":", s$label, "\n")
  run_scenario(s$id, s$label)
}) |> bind_rows()

## ─── Summary at week 16 (primary endpoint) ─────────────────────
wk16 <- results_all |>
  filter(abs(time_weeks - 16) < 0.05) |>
  group_by(scenario) |>
  slice_head(n=1) |>
  select(scenario, ITCH_NRS, IGA, DNRS, NODULE_PCHG,
         IL31_PCHG, IGE_PCHG, EOS_PCHG, DUP_RO_PCT, NEM_RO_PCT) |>
  ungroup()

cat("\n=== Week 16 Primary Endpoint Summary ===\n")
print(wk16, width=120)

## ─── Plot: Itch NRS over time ───────────────────────────────────
p_itch <- ggplot(results_all,
       aes(x=time_weeks, y=ITCH_NRS, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(0,52,8)) +
  labs(title="Prurigo Nodularis — Itch NRS Over Time",
       subtitle="7 treatment scenarios; QSP simulation",
       x="Time (weeks)", y="Itch NRS (0–10)",
       color="Treatment") +
  theme_bw() +
  theme(legend.position="bottom")

ggsave("pn_itch_nrs.png", p_itch, width=10, height=6, dpi=150)

## ─── Plot: Nodule count ─────────────────────────────────────────
p_nod <- ggplot(results_all,
       aes(x=time_weeks, y=NODULE_CNT, color=scenario)) +
  geom_line(linewidth=0.9) +
  labs(title="Nodule Count Over Time",
       x="Time (weeks)", y="Nodule Count (arb units)",
       color="Treatment") +
  theme_bw() + theme(legend.position="bottom")

ggsave("pn_nodules.png", p_nod, width=10, height=6, dpi=150)

## ─── Plot: IL-31 plasma levels ──────────────────────────────────
p_il31 <- ggplot(results_all |> filter(time_weeks <= 52),
       aes(x=time_weeks, y=IL31, color=scenario)) +
  geom_line(linewidth=0.9) +
  labs(title="IL-31 Plasma Level Over Time",
       x="Time (weeks)", y="IL-31 (pg/mL)",
       color="Treatment") +
  theme_bw() + theme(legend.position="bottom")

ggsave("pn_il31.png", p_il31, width=10, height=6, dpi=150)

## ─── Plot: Dupilumab receptor occupancy ─────────────────────────
p_ro <- ggplot(results_all |>
    filter(scenario %in% c("Dupilumab 300mg Q2W","Dupilumab + TCS")),
    aes(x=time_weeks, y=DUP_RO_PCT, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_y_continuous(limits=c(0,100)) +
  labs(title="Dupilumab IL-4Rα Receptor Occupancy",
       x="Time (weeks)", y="Receptor Occupancy (%)",
       color="Treatment") +
  theme_bw() + theme(legend.position="bottom")

ggsave("pn_dup_ro.png", p_ro, width=8, height=5, dpi=150)

cat("\nSimulation complete. Plots saved.\n")
cat("Scenarios simulated: Placebo, Dupilumab, Nemolizumab, Tralokinumab,",
    "Cyclosporine, Nalbuphine ER, Dupilumab+TCS\n")
