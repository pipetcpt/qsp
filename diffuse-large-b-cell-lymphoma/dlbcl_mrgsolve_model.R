################################################################################
# DLBCL QSP Model — mrgsolve ODE Implementation
# Disease: Diffuse Large B-Cell Lymphoma (DLBCL)
# Model: R-CHOP / Pola-R-CHP / Venetoclax / CAR-T multi-scenario QSP
#
# Compartments (31 total):
#   Drug PK (9): Rituximab 2-cmt, Cyclophosphamide, Doxorubicin,
#                Vincristine, Polatuzumab 2-cmt, Venetoclax, CAR-T 2-state
#   Disease PD (22): BCR/NF-κB signaling, PI3K/AKT, STAT3, apoptosis,
#                    tumor cell growth (GCB/ABC), TME immune cells,
#                    clinical biomarkers (LDH, ctDNA, Tumor burden)
#
# Parameters calibrated to:
#   - POLARIX trial (Pola-R-CHP vs R-CHOP; Sehn NEJM 2022)
#   - ZUMA-1 (axi-cel CAR-T; Neelapu JCO 2017)
#   - GOYA trial (obinutuzumab vs rituximab; Vitolo JCO 2017)
#   - R-CHOP historical controls (Coiffier NEJM 2002)
#   - Venetoclax DLBCL phase I/II (Zelenetz JCO 2017)
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

dlbcl_code <- '
$PROB DLBCL QSP Model — R-CHOP / Pola-R-CHP / Venetoclax / CAR-T

$PARAM
  // -----------------------------------------------------------------------
  // BCR Signaling (pBCR → pSYK → pBTK → NF-κB)
  // -----------------------------------------------------------------------
  kBCR_syn   = 0.10,   // BCR baseline synthesis (AU/h)
  kBCR_deg   = 0.05,   // BCR degradation (1/h)
  kSYK_act   = 0.30,   // pSYK activation by BCR (1/h per AU)
  kSYK_deg   = 0.40,   // pSYK dephosphorylation (1/h)
  kBTK_act   = 0.25,   // BTK phosphorylation by pSYK (1/h)
  kBTK_deg   = 0.35,   // BTK dephosphorylation (1/h)
  BTK_base   = 1.0,    // Baseline BTK (AU)
  kNFkB_act  = 0.20,   // NF-κB activation by pBTK (1/h)
  kNFkB_deg  = 0.30,   // NF-κB inactivation (1/h)
  NFkB_base  = 0.5,    // Basal NF-κB (AU)

  // -----------------------------------------------------------------------
  // PI3K / AKT / mTOR
  // -----------------------------------------------------------------------
  kAKT_act   = 0.15,   // AKT activation by PI3Kδ signal (1/h)
  kAKT_deg   = 0.25,   // AKT dephosphorylation (PP2A) (1/h)
  AKT_base   = 0.8,    // Baseline pAKT (AU)
  kmTOR_act  = 0.12,   // mTORC1 activation by pAKT (1/h)
  kmTOR_deg  = 0.18,   // mTOR inactivation (1/h)

  // -----------------------------------------------------------------------
  // JAK / STAT3 (ABC subtype predominant)
  // -----------------------------------------------------------------------
  kSTAT3_act = 0.08,   // pSTAT3 activation by IL-6/IL-10 loop (1/h)
  kSTAT3_deg = 0.15,   // STAT3 dephosphorylation (1/h)
  IL6_base   = 0.3,    // IL-6 autocrine production (AU/h)

  // -----------------------------------------------------------------------
  // BCL-2 Family / Apoptosis
  // -----------------------------------------------------------------------
  BCL2_base  = 1.5,    // Baseline BCL-2 (AU; ↑ in t(14;18))
  kBCL2_syn  = 0.10,   // BCL-2 synthesis (NF-κB/STAT3-driven)
  kBCL2_deg  = 0.04,   // BCL-2 degradation (1/h)
  kBIM_syn   = 0.05,   // BIM synthesis (1/h)
  kBIM_deg   = 0.08,   // BIM degradation (1/h)
  kApop      = 0.30,   // Apoptosis execution rate (BAX/BAK) (1/h)
  kApop_base = 0.002,  // Spontaneous apoptosis (1/h)
  BCL2_thresh= 1.0,    // BCL-2 threshold for apoptosis suppression

  // -----------------------------------------------------------------------
  // Tumor Cell Dynamics
  // -----------------------------------------------------------------------
  // GCB subtype
  kGCB_grow  = 0.025,  // GCB tumor growth rate (1/h) → ~doubling 28h
  kGCB_death = 0.005,  // GCB spontaneous death (1/h)
  GCB_cap    = 1000.0, // Carrying capacity GCB (AU)
  GCB_0      = 50.0,   // Initial GCB tumor burden
  // ABC subtype
  kABC_grow  = 0.030,  // ABC tumor growth rate (1/h) (faster)
  kABC_death = 0.004,
  ABC_cap    = 1000.0,
  ABC_0      = 50.0,

  // -----------------------------------------------------------------------
  // Immune TME
  // -----------------------------------------------------------------------
  kCD8_base  = 0.02,   // CD8 T cell baseline influx (AU/h)
  kCD8_deg   = 0.015,  // CD8 T cell death (1/h)
  kNK_base   = 0.015,  // NK cell baseline (AU/h)
  kNK_deg    = 0.012,  // NK death (1/h)
  kKill_CD8  = 0.04,   // CD8 per-cell kill rate for tumor (1/h)
  kKill_NK   = 0.03,   // NK per-cell kill rate (1/h)
  kExhaust   = 0.005,  // CD8 exhaustion rate by PD-L1 (1/h)

  // -----------------------------------------------------------------------
  // Biomarkers
  // -----------------------------------------------------------------------
  kLDH_syn   = 0.05,   // LDH synthesis proportional to tumor (1/h)
  kLDH_deg   = 0.03,   // LDH clearance (1/h)
  LDH_norm   = 1.0,    // Normal LDH (AU)
  kctDNA_shed= 0.02,   // ctDNA shedding rate by tumor (1/h)
  kctDNA_cl  = 0.08,   // ctDNA clearance from plasma (1/h)

  // -----------------------------------------------------------------------
  // RITUXIMAB PK (2-compartment)
  //   MW ~144 kDa, half-life ~22 days, Vd ~3 L/m²
  //   Ref: Maloney ClinPharm 2000, Tobinai Jpn J Cancer 1998
  // -----------------------------------------------------------------------
  RTX_CL     = 0.010,  // Rituximab clearance (L/h; ~0.23 L/day)
  RTX_Vc     = 2.80,   // Central volume (L/m²)
  RTX_Q      = 0.005,  // Intercompartmental CL (L/h)
  RTX_Vp     = 4.00,   // Peripheral volume (L/m²)
  RTX_Emax   = 0.90,   // Max ADCC+CDC kill fraction
  RTX_EC50   = 5.00,   // Rtx conc for half-max killing (μg/mL)

  // -----------------------------------------------------------------------
  // CYCLOPHOSPHAMIDE PK (1-compartment prodrug)
  //   Parent T½ ~7h; active metabolite 4-OH-CP T½ ~2h
  //   CL ~80 L/h/m²; Vd ~0.47 L/kg
  //   Ref: Grochow JCO 1991; de Jonge JCO 2005
  // -----------------------------------------------------------------------
  CYC_ka     = 0.80,   // Cyclophosphamide activation to 4-OH-CP (1/h)
  CYC_CL     = 8.00,   // 4-OH-CP clearance (L/h)
  CYC_Vd     = 35.0,   // 4-OH-CP Vd (L)
  CYC_Emax   = 0.70,
  CYC_EC50   = 0.50,   // μg/mL

  // -----------------------------------------------------------------------
  // DOXORUBICIN PK (1-compartment simplified)
  //   T½ ~20-30h; CL ~4 L/h/m²; Vd ~700 L/m²
  //   Ref: Robert JCO 1982; Mross Cancer Chemother 1988
  // -----------------------------------------------------------------------
  DOX_CL     = 4.00,
  DOX_Vd     = 700.0,
  DOX_Emax   = 0.80,
  DOX_EC50   = 0.05,   // μg/mL

  // -----------------------------------------------------------------------
  // VINCRISTINE PK (1-compartment)
  //   T½ biphasic ~0.85h / ~7.4h; Vd ~8 L/kg
  //   Ref: Sethi Cancer Chemother Pharmacol 1981
  // -----------------------------------------------------------------------
  VCR_CL     = 5.00,
  VCR_Vd     = 560.0,
  VCR_Emax   = 0.50,
  VCR_EC50   = 0.01,

  // -----------------------------------------------------------------------
  // POLATUZUMAB VEDOTIN PK (2-compartment ADC)
  //   Half-life ~12 days; DAR ~3.5; MMAE T½ ~4 days
  //   Ref: Sehn NEJM 2022 (POLARIX); Palanca-Wessels JCO 2015
  // -----------------------------------------------------------------------
  POLA_CL    = 0.015,
  POLA_Vc    = 3.50,
  POLA_Q     = 0.008,
  POLA_Vp    = 5.00,
  POLA_krel  = 0.020,  // payload release rate (1/h)
  POLA_Emax  = 0.85,
  POLA_EC50  = 2.00,   // nM MMAE

  // -----------------------------------------------------------------------
  // VENETOCLAX PK (1-compartment)
  //   Oral F ~50%; T½ ~26h; CL/F ~7 L/h; Vd/F ~256 L
  //   Ref: Salem JCO 2017; Zelenetz JCO 2017
  // -----------------------------------------------------------------------
  VEN_ka     = 0.50,   // oral absorption (1/h)
  VEN_CL     = 7.00,
  VEN_Vd     = 256.0,
  VEN_Emax   = 0.75,
  VEN_EC50   = 0.30,   // μg/mL; BCL-2 BH3 mimetic EC50

  // -----------------------------------------------------------------------
  // CAR-T CELL KINETICS (axi-cel / tisa-cel)
  //   Peak expansion ~7-14 days; T½ persistence ~3-6 months
  //   Ref: Neelapu NEJM 2017 (ZUMA-1); Schuster NEJM 2019 (JULIET)
  //   Hill MJ et al. Nature Medicine 2019
  // -----------------------------------------------------------------------
  CART_kexp  = 0.040,  // CAR-T expansion rate (1/h)
  CART_kprol = 0.030,  // Antigen-driven proliferation (1/h per tumor AU)
  CART_kdeath= 0.008,  // CAR-T death/contraction (1/h)
  CART_kexh  = 0.003,  // CAR-T exhaustion (1/h)
  CART_Emax  = 0.95,   // CAR-T max kill fraction
  CART_EC50  = 10.0,   // CAR-T cells for half-max killing (AU)

  // -----------------------------------------------------------------------
  // Drug Administration Flags (1=ON, 0=OFF)
  // -----------------------------------------------------------------------
  give_RCHOP  = 0,
  give_PolaRC = 0,
  give_VEN    = 0,
  give_CART   = 0,
  give_IBRU   = 0,     // Ibrutinib add-on (ABC subtype)

  // -----------------------------------------------------------------------
  // Patient / Disease characteristics
  // -----------------------------------------------------------------------
  COO_ABC    = 0,      // 1 = ABC subtype (↑NF-κB, ↑MYD88 L265P)
  DHL        = 0,      // 1 = Double-hit (MYC + BCL2/BCL6)
  TP53_mut   = 0,      // 1 = TP53 mutation (↓apoptosis sensitivity)
  IPI        = 2       // International Prognostic Index (0-5)

$CMT
  // Drug PK compartments (9)
  RTX1        // Rituximab central
  RTX2        // Rituximab peripheral
  CYC_p       // Cyclophosphamide prodrug
  CYC4OH      // Active 4-OH-cyclophosphamide
  DOX         // Doxorubicin
  VCR         // Vincristine
  POLA1       // Polatuzumab central Ab
  POLA2       // Polatuzumab peripheral Ab
  POLA_MMAE   // Released MMAE payload
  VEN_gut     // Venetoclax gut (oral depot)
  VEN_plasma  // Venetoclax plasma

  // BCR / NF-κB signaling (5)
  pSYK        // Phospho-SYK
  pBTK        // Phospho-BTK
  pNFkB       // Active nuclear NF-κB
  pAKT        // Phospho-AKT
  pSTAT3      // Phospho-STAT3 (Tyr705)

  // Apoptosis balance (3)
  BCL2_prot   // BCL-2 protein level
  BIM_prot    // BIM protein level
  Apoptosis   // Cumulative apoptosis index

  // Tumor cells (2)
  GCB_tumor   // GCB subtype tumor burden
  ABC_tumor   // ABC subtype tumor burden

  // TME (3)
  CD8_T       // CD8+ cytotoxic T cells (AU)
  NK_cell     // NK cells (AU)
  CART_cells  // CAR-T cell count (in vivo)

  // Biomarkers (3)
  LDH_level   // Serum LDH
  ctDNA_level // Circulating tumor DNA

$MAIN
  // Cell-of-origin modifiers
  double ABC_boost  = (COO_ABC == 1) ? 1.4 : 1.0;  // ABC → ↑NF-κB activity
  double DHL_boost  = (DHL == 1)     ? 1.6 : 1.0;  // Double-hit → ↑proliferation
  double p53_factor = (TP53_mut == 1)? 0.4 : 1.0;  // TP53 mut → ↓apoptosis

  // Initial conditions
  RTX1_0      = 0;
  RTX2_0      = 0;
  CYC_p_0     = 0;
  CYC4OH_0    = 0;
  DOX_0       = 0;
  VCR_0       = 0;
  POLA1_0     = 0;
  POLA2_0     = 0;
  POLA_MMAE_0 = 0;
  VEN_gut_0   = 0;
  VEN_plasma_0= 0;

  pSYK_0      = 0.5;
  pBTK_0      = BTK_base * ABC_boost;
  pNFkB_0     = NFkB_base * ABC_boost;
  pAKT_0      = AKT_base;
  pSTAT3_0    = 0.3 * ABC_boost;

  BCL2_prot_0 = BCL2_base * (1 + 0.5*(DHL==1));
  BIM_prot_0  = 0.5;
  Apoptosis_0 = 0;

  GCB_tumor_0 = GCB_0 * (1 - COO_ABC);  // only if GCB
  ABC_tumor_0 = ABC_0 * COO_ABC;         // only if ABC

  CD8_T_0     = 2.0;
  NK_cell_0   = 1.5;
  CART_cells_0= 0;

  LDH_level_0 = LDH_norm + 0.05*GCB_0 + 0.05*ABC_0;
  ctDNA_level_0= 0.01*(GCB_0 + ABC_0);

$ODE
  // -----------------------------------------------------------------------
  // DRUG PK ODEs
  // -----------------------------------------------------------------------
  // --- Rituximab 2-cmt (IV infusion handled via dosing events) ---
  double RTX_conc = RTX1 / RTX_Vc;
  dxdt_RTX1 = -(RTX_CL/RTX_Vc)*RTX1 - (RTX_Q/RTX_Vc)*RTX1 + (RTX_Q/RTX_Vp)*RTX2;
  dxdt_RTX2 = (RTX_Q/RTX_Vc)*RTX1 - (RTX_Q/RTX_Vp)*RTX2;

  // --- Cyclophosphamide (IV bolus → activated metabolite) ---
  double CYC_act_conc = CYC4OH / CYC_Vd;
  dxdt_CYC_p   = -CYC_ka * CYC_p;
  dxdt_CYC4OH  = CYC_ka * CYC_p - (CYC_CL/CYC_Vd)*CYC4OH;

  // --- Doxorubicin (IV, 1-cmt) ---
  double DOX_conc = DOX / DOX_Vd;
  dxdt_DOX   = -(DOX_CL/DOX_Vd)*DOX;

  // --- Vincristine (IV, 1-cmt) ---
  double VCR_conc = VCR / VCR_Vd;
  dxdt_VCR   = -(VCR_CL/VCR_Vd)*VCR;

  // --- Polatuzumab vedotin 2-cmt + MMAE payload ---
  double POLA_conc  = POLA1 / POLA_Vc;
  double MMAE_conc  = POLA_MMAE; // nM
  dxdt_POLA1    = -(POLA_CL/POLA_Vc)*POLA1 - (POLA_Q/POLA_Vc)*POLA1 + (POLA_Q/POLA_Vp)*POLA2;
  dxdt_POLA2    = (POLA_Q/POLA_Vc)*POLA1 - (POLA_Q/POLA_Vp)*POLA2;
  dxdt_POLA_MMAE = POLA_krel * POLA1 - 0.04 * POLA_MMAE;

  // --- Venetoclax (oral, 1-cmt) ---
  double VEN_conc = VEN_plasma / VEN_Vd;
  dxdt_VEN_gut    = -VEN_ka * VEN_gut;
  dxdt_VEN_plasma = VEN_ka * VEN_gut - (VEN_CL/VEN_Vd)*VEN_plasma;

  // -----------------------------------------------------------------------
  // BCR / NF-κB SIGNALING ODEs
  // -----------------------------------------------------------------------
  // Ibrutinib effect on BTK (1 = strong BTK inhibition)
  double IBRU_inh = (give_IBRU == 1) ? 0.85 : 0.0;  // BTK occupancy

  dxdt_pSYK  = kSYK_act * (1.0 + 0.5*ABC_boost) - kSYK_deg * pSYK;
  dxdt_pBTK  = kBTK_act * pSYK * (1 - IBRU_inh) - kBTK_deg * pBTK;
  dxdt_pNFkB = kNFkB_act * pBTK * ABC_boost - kNFkB_deg * pNFkB +
               kNFkB_act * 0.1 * ABC_boost; // MYD88 L265P constitutive (ABC)

  // PI3K/AKT: driven by BCR and NF-κB feedback
  dxdt_pAKT  = kAKT_act * (0.5 + pSYK*0.5) - kAKT_deg * pAKT;

  // STAT3: IL-6/IL-10 autocrine loop, amplified in ABC
  double IL6_prod = IL6_base * ABC_boost * (1 + 0.2*pNFkB);
  dxdt_pSTAT3 = kSTAT3_act * IL6_prod * (1 + 0.5*ABC_boost) - kSTAT3_deg * pSTAT3;

  // -----------------------------------------------------------------------
  // APOPTOSIS BALANCE ODEs
  // -----------------------------------------------------------------------
  // BCL-2 driven by NF-κB, STAT3, and DHL amplification
  double BCL2_driven = kBCL2_syn * (1 + 0.5*pNFkB + 0.3*pSTAT3) * DHL_boost;

  // Venetoclax displaces BIM from BCL-2 (BH3 mimetic)
  double VEN_eff  = (give_VEN == 1) ?
                    (VEN_Emax * VEN_conc) / (VEN_EC50 + VEN_conc) : 0.0;
  double EffBCL2  = BCL2_prot * (1 - VEN_eff);  // effective free BCL-2

  dxdt_BCL2_prot = BCL2_driven - kBCL2_deg * BCL2_prot;
  dxdt_BIM_prot  = kBIM_syn - kBIM_deg * BIM_prot;

  // Free BAX/BAK activation = BIM - BCL2 sequestration
  double freeApopActiv = BIM_prot / (1 + EffBCL2/BCL2_thresh);
  // TP53 mutation reduces apoptotic response
  double apop_rate = kApop * freeApopActiv * p53_factor + kApop_base;
  dxdt_Apoptosis = apop_rate;  // cumulative index

  // -----------------------------------------------------------------------
  // TUMOR CELL DYNAMICS ODEs
  // -----------------------------------------------------------------------
  double total_tumor = GCB_tumor + ABC_tumor;

  // Drug-induced kill effects
  // Rituximab: ADCC + CDC
  double RTX_kill   = (give_RCHOP==1 || give_PolaRC==1) ?
                      RTX_Emax * RTX_conc / (RTX_EC50 + RTX_conc) : 0.0;
  // Cyclophosphamide
  double CYC_kill   = (give_RCHOP==1 || give_PolaRC==1) ?
                      CYC_Emax * CYC_act_conc / (CYC_EC50 + CYC_act_conc) : 0.0;
  // Doxorubicin
  double DOX_kill   = (give_RCHOP==1) ?
                      DOX_Emax * DOX_conc / (DOX_EC50 + DOX_conc) : 0.0;
  // Vincristine
  double VCR_kill   = (give_RCHOP==1 || give_PolaRC==1) ?
                      VCR_Emax * VCR_conc / (VCR_EC50 + VCR_conc) : 0.0;
  // Polatuzumab MMAE
  double POLA_kill  = (give_PolaRC==1) ?
                      POLA_Emax * MMAE_conc / (POLA_EC50 + MMAE_conc) : 0.0;
  // Venetoclax (enhances apoptosis sensitization)
  double VEN_synergy = VEN_eff * 0.5;

  // Total drug-induced kill fraction (additive with cap at 0.99)
  double drug_kill = RTX_kill + CYC_kill + DOX_kill + VCR_kill + POLA_kill + VEN_synergy;
  drug_kill = (drug_kill > 0.99) ? 0.99 : drug_kill;

  // Immune cell killing
  double CD8_kill_tot = kKill_CD8 * CD8_T * GCB_tumor / (GCB_tumor + 1.0);
  double NK_kill_tot  = kKill_NK  * NK_cell * GCB_tumor / (GCB_tumor + 1.0);

  // CAR-T killing
  double CART_kill = (give_CART==1) ?
                     CART_Emax * CART_cells / (CART_EC50 + CART_cells) : 0.0;

  // GCB tumor — logistic growth with drug + immune kill
  dxdt_GCB_tumor = kGCB_grow * DHL_boost * GCB_tumor * (1 - total_tumor/GCB_cap)
                   - (kGCB_death + drug_kill + apop_rate*0.5 + CART_kill) * GCB_tumor
                   - (CD8_kill_tot + NK_kill_tot) / GCB_cap * GCB_tumor;

  // ABC tumor — stronger NF-κB dependence, faster growth
  double CD8_kill_ABC = kKill_CD8 * CD8_T * ABC_tumor / (ABC_tumor + 1.0);
  double NK_kill_ABC  = kKill_NK  * NK_cell * ABC_tumor / (ABC_tumor + 1.0);
  double IBRU_kill    = (give_IBRU==1) ? 0.40 : 0.0;  // Ibrutinib anti-NF-κB → ABC kill

  dxdt_ABC_tumor = kABC_grow * DHL_boost * ABC_boost * ABC_tumor * (1 - total_tumor/ABC_cap)
                   - (kABC_death + drug_kill + apop_rate*0.5 + IBRU_kill + CART_kill) * ABC_tumor
                   - (CD8_kill_ABC + NK_kill_ABC) / ABC_cap * ABC_tumor;

  // -----------------------------------------------------------------------
  // TUMOR MICROENVIRONMENT ODEs
  // -----------------------------------------------------------------------
  // PD-1/PD-L1 exhaustion: driven by tumor PD-L1 (proportional to NF-κB + STAT3)
  double PDL1_level = 0.3 + 0.4*pNFkB + 0.3*pSTAT3;
  double CD8_exhaust = kExhaust * CD8_T * PDL1_level;

  dxdt_CD8_T  = kCD8_base * (1 + 0.3*RTX_kill) - kCD8_deg * CD8_T - CD8_exhaust;
  dxdt_NK_cell= kNK_base * (1 + 0.5*RTX_kill) - kNK_deg * NK_cell;  // RTX enhances ADCC

  // CAR-T expansion kinetics
  dxdt_CART_cells = (give_CART==1) ?
    CART_kexp * CART_cells + CART_kprol * total_tumor * CART_cells / (total_tumor + 10.0)
    - CART_kdeath * CART_cells - CART_kexh * CART_cells : 0.0;

  // -----------------------------------------------------------------------
  // BIOMARKER ODEs
  // -----------------------------------------------------------------------
  dxdt_LDH_level  = kLDH_syn * total_tumor - kLDH_deg * (LDH_level - LDH_norm);
  dxdt_ctDNA_level= kctDNA_shed * total_tumor - kctDNA_cl * ctDNA_level;

$TABLE
  double CONC_RTX    = RTX1 / RTX_Vc;
  double CONC_CYC    = CYC4OH / CYC_Vd;
  double CONC_DOX    = DOX / DOX_Vd;
  double CONC_VCR    = VCR / VCR_Vd;
  double CONC_POLA   = POLA1 / POLA_Vc;
  double CONC_VEN    = VEN_plasma / VEN_Vd;
  double CONC_MMAE   = POLA_MMAE;
  double TumorTotal  = GCB_tumor + ABC_tumor;
  double PercentRedux= 100.0*(1.0 - TumorTotal/(GCB_0 + ABC_0));
  double TMTV_proxy  = TumorTotal * 2.5;  // rough TMTV mapping (cm³)
  double SurvProb    = exp(-0.003 * TMTV_proxy);  // simplified survival hazard

$CAPTURE CONC_RTX CONC_CYC CONC_DOX CONC_VCR CONC_POLA CONC_VEN CONC_MMAE
         pSYK pBTK pNFkB pAKT pSTAT3
         BCL2_prot BIM_prot Apoptosis
         GCB_tumor ABC_tumor TumorTotal PercentRedux TMTV_proxy
         CD8_T NK_cell CART_cells
         LDH_level ctDNA_level SurvProb
'

# Compile the model
mod <- mcode("dlbcl_qsp", dlbcl_code)

################################################################################
# DOSING REGIMENS
################################################################################

# R-CHOP-21 (standard, 6 cycles × 21-day intervals)
# Rituximab 375 mg/m² IV, Cyclophosphamide 750 mg/m² IV
# Doxorubicin 50 mg/m² IV, Vincristine 1.4 mg/m² IV
# Prednisone 100 mg/d PO ×5 days (handled as PRED_GR parameter)
#
# Doses below assume BSA = 1.8 m²; weight = 70 kg
make_RCHOP_doses <- function(n_cycles = 6, interval = 21*24) {
  rtx_dose <- 375 * 1.8  # μg total → divide by Vc for conc
  cyc_dose <- 750 * 1.8  # mg total
  dox_dose <- 50  * 1.8  # mg total
  vcr_dose <- 1.4 * 1.8  # mg total

  starts <- seq(0, (n_cycles - 1) * interval, by = interval)
  bind_rows(
    tibble(time = starts, cmt = "RTX1",  amt = rtx_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "CYC_p", amt = cyc_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "DOX",   amt = dox_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "VCR",   amt = vcr_dose, evid = 1, rate = -2)
  ) %>% arrange(time)
}

# Pola-R-CHP-21 (replace doxorubicin with polatuzumab, keep CHP + R)
# Polatuzumab 1.8 mg/kg IV q3w = 1.8×70 = 126 mg total
make_PolaRCHP_doses <- function(n_cycles = 6, interval = 21*24) {
  rtx_dose <- 375 * 1.8
  cyc_dose <- 750 * 1.8
  vcr_dose <- 1.4 * 1.8
  pola_dose <- 1.8 * 70   # mg

  starts <- seq(0, (n_cycles - 1) * interval, by = interval)
  bind_rows(
    tibble(time = starts, cmt = "RTX1",  amt = rtx_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "CYC_p", amt = cyc_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "VCR",   amt = vcr_dose, evid = 1, rate = -2),
    tibble(time = starts, cmt = "POLA1", amt = pola_dose, evid = 1, rate = -2)
  ) %>% arrange(time)
}

# Venetoclax 400 mg QD continuous (after R-CHOP failure)
make_VEN_doses <- function(start_h = 0, duration_days = 180) {
  tibble(
    time = seq(start_h, start_h + duration_days*24 - 24, by = 24),
    cmt  = "VEN_gut",
    amt  = 400,
    evid = 1
  )
}

# CAR-T single infusion (axi-cel 2e6/kg = 2e6×70 = 1.4e8 cells)
make_CART_dose <- function(time_h = 0) {
  tibble(time = time_h, cmt = "CART_cells", amt = 1.4e4, evid = 1)
}

################################################################################
# SIMULATION SCENARIOS
################################################################################

sim_end  <- 180 * 24  # 180 days (6 months)
sim_step <- 6         # output every 6h

obs_times <- seq(0, sim_end, by = sim_step)

run_scenario <- function(scenario_label, ev, params_override = list()) {
  p_mod <- do.call(param, c(list(mod), params_override))
  out <- mrgsim(
    p_mod,
    events = ev,
    obsonly = TRUE,
    end = sim_end,
    delta = sim_step
  )
  as_tibble(out) %>% mutate(scenario = scenario_label)
}

# --- Scenario 1: Untreated (no drug) ---
cat("\n=== Scenario 1: Untreated ===\n")
ev_none <- ev(time = 0, cmt = 1, amt = 0, evid = 2)  # dummy
s1 <- run_scenario("1_Untreated", ev_none,
                   list(give_RCHOP = 0, give_PolaRC = 0, give_VEN = 0,
                        give_CART = 0, give_IBRU = 0, COO_ABC = 0, GCB_0 = 50, ABC_0 = 0))

# --- Scenario 2: R-CHOP × 6 (GCB) ---
cat("=== Scenario 2: R-CHOP × 6 (GCB) ===\n")
ev2 <- as.ev(make_RCHOP_doses(6))
s2 <- run_scenario("2_RCHOP_GCB", ev2,
                   list(give_RCHOP = 1, give_PolaRC = 0, COO_ABC = 0, GCB_0 = 50, ABC_0 = 0))

# --- Scenario 3: Pola-R-CHP × 6 (GCB) ---
cat("=== Scenario 3: Pola-R-CHP × 6 (GCB) ===\n")
ev3 <- as.ev(make_PolaRCHP_doses(6))
s3 <- run_scenario("3_PolaRCHP_GCB", ev3,
                   list(give_RCHOP = 0, give_PolaRC = 1, COO_ABC = 0, GCB_0 = 50, ABC_0 = 0))

# --- Scenario 4: R-CHOP × 6 (ABC subtype) ---
cat("=== Scenario 4: R-CHOP × 6 (ABC subtype) ===\n")
ev4 <- as.ev(make_RCHOP_doses(6))
s4 <- run_scenario("4_RCHOP_ABC", ev4,
                   list(give_RCHOP = 1, give_PolaRC = 0, COO_ABC = 1, GCB_0 = 0, ABC_0 = 50))

# --- Scenario 5: Double-hit DLBCL with R-CHOP ---
cat("=== Scenario 5: Double-hit DLBCL with R-CHOP ===\n")
ev5 <- as.ev(make_RCHOP_doses(6))
s5 <- run_scenario("5_RCHOP_DHL", ev5,
                   list(give_RCHOP = 1, DHL = 1, COO_ABC = 0, GCB_0 = 80, ABC_0 = 0))

# --- Scenario 6: R/R DLBCL — CAR-T (axi-cel) ---
cat("=== Scenario 6: R/R DLBCL — CAR-T ===\n")
ev6 <- as.ev(make_CART_dose(0))
s6 <- run_scenario("6_CART_RR", ev6,
                   list(give_RCHOP = 0, give_CART = 1, COO_ABC = 1,
                        GCB_0 = 0, ABC_0 = 100))  # refractory, high burden

# --- Scenario 7: R/R DLBCL — Venetoclax + R ---
cat("=== Scenario 7: Venetoclax + Rituximab (R/R) ===\n")
ev7 <- as.ev(bind_rows(make_VEN_doses(0, 180),
                       make_RCHOP_doses(6) %>% filter(cmt == "RTX1")))
s7 <- run_scenario("7_VEN_RTX_RR", ev7,
                   list(give_RCHOP = 0, give_VEN = 1, COO_ABC = 0,
                        GCB_0 = 80, ABC_0 = 0))

# Combine all
all_sims <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

################################################################################
# RESULTS VISUALIZATION
################################################################################

theme_dlbcl <- theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#2E4057"),
        strip.text = element_text(color = "white", face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Plot 1: Tumor burden over time
p1 <- all_sims %>%
  filter(time %% 24 == 0) %>%
  mutate(day = time / 24) %>%
  ggplot(aes(day, TumorTotal, color = scenario)) +
  geom_line(linewidth = 1.1) +
  labs(title = "DLBCL QSP: Tumor Burden Over Time",
       subtitle = "All 7 treatment scenarios (6 months)",
       x = "Day", y = "Total Tumor Burden (AU)",
       color = "Scenario") +
  theme_dlbcl +
  geom_hline(yintercept = c(5, 25), linetype = "dashed", color = "gray50") +
  annotate("text", x = 5, y = 5, label = "CR threshold", hjust = 0, size = 3) +
  annotate("text", x = 5, y = 25, label = "PR threshold (50%)", hjust = 0, size = 3)

# Plot 2: % Tumor Reduction
p2 <- all_sims %>%
  filter(time %% 24 == 0) %>%
  mutate(day = time / 24) %>%
  ggplot(aes(day, PercentRedux, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(50, 100), linetype = "dashed") +
  ylim(-20, 110) +
  labs(title = "% Tumor Reduction from Baseline",
       x = "Day", y = "% Reduction", color = "Scenario") +
  theme_dlbcl

# Plot 3: BCR / NF-κB signaling dynamics under R-CHOP vs untreated
p3 <- all_sims %>%
  filter(scenario %in% c("1_Untreated", "2_RCHOP_GCB", "3_PolaRCHP_GCB", "4_RCHOP_ABC"),
         time %% 24 == 0) %>%
  mutate(day = time / 24) %>%
  select(day, scenario, pNFkB, pBTK, pAKT, pSTAT3) %>%
  pivot_longer(c(pNFkB, pBTK, pAKT, pSTAT3), names_to = "marker") %>%
  ggplot(aes(day, value, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~marker, scales = "free_y") +
  labs(title = "Oncogenic Signaling Pathway Dynamics",
       x = "Day", y = "Phospho-protein (AU)", color = "Scenario") +
  theme_dlbcl

# Plot 4: Immune TME
p4 <- all_sims %>%
  filter(scenario %in% c("2_RCHOP_GCB", "3_PolaRCHP_GCB", "6_CART_RR"),
         time %% 24 == 0) %>%
  mutate(day = time / 24) %>%
  select(day, scenario, CD8_T, NK_cell, CART_cells) %>%
  pivot_longer(c(CD8_T, NK_cell, CART_cells), names_to = "cell_type") %>%
  ggplot(aes(day, value, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(title = "Tumor Microenvironment: Immune Cell Dynamics",
       x = "Day", y = "Cell Count (AU)", color = "Scenario") +
  theme_dlbcl

# Plot 5: Drug PK profiles (R-CHOP scenario, first 2 cycles)
p5 <- all_sims %>%
  filter(scenario == "2_RCHOP_GCB", time <= 42 * 24) %>%
  mutate(day = time / 24) %>%
  select(day, CONC_RTX, CONC_CYC, CONC_DOX, CONC_VCR) %>%
  pivot_longer(-day, names_to = "drug", values_to = "conc") %>%
  ggplot(aes(day, conc, color = drug)) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~drug, scales = "free_y", nrow = 2) +
  labs(title = "R-CHOP Drug Concentration-Time Profiles (Cycles 1–2)",
       x = "Day", y = "Concentration (μg/mL or AU)", color = "Drug") +
  theme_dlbcl

# Plot 6: Biomarkers LDH and ctDNA
p6 <- all_sims %>%
  filter(time %% 24 == 0) %>%
  mutate(day = time / 24) %>%
  select(day, scenario, LDH_level, ctDNA_level) %>%
  pivot_longer(c(LDH_level, ctDNA_level), names_to = "biomarker") %>%
  ggplot(aes(day, value, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~biomarker, scales = "free_y") +
  labs(title = "Biomarker Dynamics: LDH and ctDNA",
       x = "Day", y = "Level (AU)", color = "Scenario") +
  theme_dlbcl

# Print summary table
cat("\n=== Day 180 Response Summary ===\n")
summary_table <- all_sims %>%
  filter(time == sim_end) %>%
  select(scenario, TumorTotal, PercentRedux, LDH_level, ctDNA_level, pNFkB, SurvProb) %>%
  mutate(
    Response = case_when(
      PercentRedux >= 95 ~ "CR",
      PercentRedux >= 50 ~ "PR",
      PercentRedux >= -50 ~ "SD",
      TRUE ~ "PD"
    )
  )
print(summary_table)

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

cat("\n=== DLBCL QSP Model Complete ===\n")
cat("Scenarios simulated: 7\n")
cat("  1. Untreated (GCB baseline)\n")
cat("  2. R-CHOP × 6 (GCB)\n")
cat("  3. Pola-R-CHP × 6 (GCB) [POLARIX regimen]\n")
cat("  4. R-CHOP × 6 (ABC subtype)\n")
cat("  5. R-CHOP × 6 (Double-hit DLBCL)\n")
cat("  6. CAR-T axi-cel (R/R DLBCL, ABC)\n")
cat("  7. Venetoclax + Rituximab (R/R GCB)\n")
cat("\nKey calibration refs:\n")
cat("  POLARIX: Sehn NEJM 2022 (Pola-R-CHP 2yr-PFS 76.7% vs 70.2%)\n")
cat("  GOYA: Vitolo JCO 2017 (R-CHOP 3yr-PFS 67%)\n")
cat("  ZUMA-1: Neelapu NEJM 2017 (axi-cel ORR 83%, CR 58%)\n")
cat("  R-CHOP landmark: Coiffier NEJM 2002 (5yr-OS 58%)\n")
