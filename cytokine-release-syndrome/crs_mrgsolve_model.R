## =============================================================================
## Cytokine Release Syndrome (CRS) — QSP Model (mrgsolve)
## =============================================================================
## Disease: Cytokine Release Syndrome following CAR-T cell therapy /
##          bispecific antibody therapy
##
## Model structure (22 ODE compartments):
##   - CAR-T cell dynamics (4 compartments)
##   - T cell–derived cytokines (4 compartments)
##   - Macrophage activation cascade (3 compartments)
##   - IL-6 / JAK-STAT signaling (3 compartments)
##   - Endothelial / vascular response (2 compartments)
##   - Organ damage / CRS severity (2 compartments)
##   - Drug PK: Tocilizumab (2-cmt), Siltuximab (1-cmt),
##              Dexamethasone (2-cmt), Ruxolitinib (2-cmt), Anakinra (1-cmt)
##
## Treatment scenarios:
##   1. Untreated / no intervention (natural CRS progression)
##   2. Tocilizumab (Grade 2 CRS, single dose Day 3)
##   3. Tocilizumab + Dexamethasone (Grade 3 CRS, Day 2)
##   4. Ruxolitinib + Dexamethasone (Refractory/Grade 4 CRS)
##   5. Anakinra + Dexamethasone (NLRP3/macrophage-dominant / pediatric)
##
## Clinical calibration references:
##   - Maude SL et al. NEJM 2014 (tisagenlecleucel)
##   - Lee DW et al. Blood 2014 (CRS grading criteria)
##   - Neelapu SS et al. NEJM 2017 (axicabtagene)
##   - Grupp SA et al. NEJM 2013 (ALL CAR-T)
##   - Teachey DT et al. Cancer Discov 2016 (cytokine biomarkers)
##   - Shah BD et al. Lancet 2021 (KymRAH Phase 3)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# mrgsolve Model Definition
# ==============================================================================

crs_model_code <- '
$PROB CRS QSP Model — CAR-T Immunotherapy

$PARAM
  // ---- CAR-T cell dynamics ----
  kCART_act     = 0.6    // CAR-T activation rate (/day; antigen-dependent)
  kCART_prolif  = 1.8    // CAR-T expansion rate (/day; IL-2 driven)
  kCART_death   = 0.25   // CAR-T death/exhaustion rate (/day)
  kCART_exh     = 0.12   // exhaustion rate (/day; TOX-dependent)
  kIL2_prolif   = 20.0   // IL-2 EC50 for T cell proliferation (pg/mL)
  fNK_frac      = 0.0    // NK dysfunction fraction (0=normal, 0.9=primary CRS)
  TUMOR0        = 1.0    // initial tumor antigen level (relative units)
  kTumor_kill   = 0.4    // tumor killing rate by CAR-T (/day)
  kTumor_grow   = 0.05   // residual tumor growth (/day)

  // ---- Cytokine production & clearance ----
  kIFNg_prod    = 0.8    // IFN-γ production rate by CART_ACT (pg/mL/cell-unit/day)
  kIL2_prod     = 0.5    // IL-2 production rate
  kGMCSF_prod   = 0.4    // GM-CSF production rate
  kTNFa_T_prod  = 0.3    // TNF-α (T cell) production rate
  dIFNg         = 1.5    // IFN-γ clearance rate (/day; t½~11h)
  dIL2          = 6.0    // IL-2 clearance rate (/day; t½~2.8h)
  dGMCSF        = 2.4    // GM-CSF clearance rate (/day; t½~7h)
  dTNFa_T       = 3.0    // TNF-α clearance rate (/day)

  // ---- Macrophage activation ----
  kMAC_act      = 0.5    // macrophage activation rate by IFN-γ/GM-CSF
  kMAC_deact    = 0.3    // macrophage deactivation rate (/day)
  EC50_IFNg_MAC = 50.0   // IFN-γ EC50 for macrophage activation (pg/mL)
  EC50_GMCSF    = 20.0   // GM-CSF EC50
  kIL6_prod     = 2.0    // IL-6 production rate by activated macrophage
  kIL1B_prod    = 0.6    // IL-1β production rate (NLRP3)
  kTNFa_M_prod  = 0.8    // TNF-α (macrophage) production rate
  dIL6          = 3.0    // IL-6 clearance (/day; t½~5.5h)
  dIL1B         = 4.0    // IL-1β clearance (/day; t½~4h)
  dTNFa_M       = 3.0    // TNF-α (macrophage) clearance
  IL6_feedback  = 0.3    // IL-6 macrophage autocrine amplification

  // ---- IL-6 downstream signaling ----
  kSTAT3_act    = 0.8    // STAT3 activation rate by IL-6/gp130
  kSTAT3_deact  = 1.2    // STAT3 deactivation (SOCS3) rate (/day)
  EC50_IL6_STAT3= 30.0   // IL-6 EC50 for STAT3 activation
  kCRP_prod     = 0.4    // CRP production rate (IL-6 driven)
  dCRP          = 0.7    // CRP clearance (/day; t½~24h)
  kFerr_prod    = 0.25   // Ferritin production rate by MAC_ACT
  dFerr         = 0.15   // Ferritin clearance (/day; t½~4.6d)

  // ---- Endothelial / vascular ----
  kENDO_act     = 0.4    // endothelial activation rate by IL-6/TNF-α
  kENDO_deact   = 0.6    // endothelial recovery rate (/day)
  kVASC_leak    = 0.5    // vascular permeability driving force
  dVASC         = 0.8    // vascular permeability recovery (/day)
  EC50_IL6_ENDO = 80.0   // IL-6 EC50 for endothelial activation

  // ---- CRS severity / organ damage ----
  kCRS_onset    = 0.6    // CRS severity accumulation rate
  kCRS_resolve  = 0.4    // natural CRS resolution rate (/day)
  kORGAN_dmg    = 0.08   // organ damage accumulation rate
  dORGAN        = 0.05   // organ damage repair (/day)

  // ---- Tocilizumab PK (2-cmt IV) ----
  CL_TOCI       = 0.22   // clearance (L/day; 8 mg/kg IV, ~70 kg adult)
  V1_TOCI       = 3.1    // central volume (L)
  Q_TOCI        = 0.29   // intercompartmental clearance (L/day)
  V2_TOCI       = 2.9    // peripheral volume (L)
  Emax_TOCI     = 0.92   // max IL-6R blockade
  EC50_TOCI     = 1.5    // EC50 for IL-6R occupancy (μg/mL)

  // ---- Siltuximab PK (1-cmt IV) ----
  CL_SILT       = 0.015  // clearance (L/day; t½~21d)
  V1_SILT       = 5.5    // volume of distribution (L)
  Emax_SILT     = 0.90   // max IL-6 neutralization
  EC50_SILT     = 2.0    // EC50 (μg/mL)

  // ---- Dexamethasone PK (2-cmt IV) ----
  CL_DEX        = 12.0   // clearance (L/h → 288 L/day; rapid)
  V1_DEX        = 25.0   // central volume (L)
  Q_DEX         = 8.0    // intercompartmental
  V2_DEX        = 65.0   // peripheral volume
  Emax_DEX_IL6  = 0.75   // max IL-6 suppression by DEX
  EC50_DEX      = 0.3    // EC50 (ng/mL)
  Emax_DEX_IFNg = 0.60   // max IFN-γ suppression

  // ---- Ruxolitinib PK (2-cmt oral) ----
  ka_RUX        = 6.0    // absorption rate (/day)
  CL_RUX        = 18.0   // clearance (L/h → 432 L/day; t½~3h)
  V1_RUX        = 72.0   // central volume (L)
  Q_RUX         = 12.0   // intercompartmental
  V2_RUX        = 144.0  // peripheral volume
  F_RUX         = 0.95   // bioavailability
  Emax_RUX_JAK  = 0.85   // max JAK1/2 inhibition (IC50 = 3.3 nM JAK1)
  EC50_RUX_JAK  = 80.0   // EC50 (ng/mL; ~ IC50 in ng/mL units)

  // ---- Anakinra PK (1-cmt SC) ----
  ka_ANA        = 4.0    // absorption rate (/day; SC t½~5h)
  CL_ANA        = 35.0   // clearance (L/day; t½~4-6h)
  V1_ANA        = 18.0   // volume (L)
  F_ANA         = 0.95   // bioavailability
  Emax_ANA_IL1  = 0.90   // max IL-1β blockade
  EC50_ANA      = 30.0   // EC50 (ng/mL)

  // ---- Treatment flags (set per scenario) ----
  use_TOCI      = 0      // 1 = tocilizumab given
  use_SILT      = 0      // 1 = siltuximab given
  use_DEX       = 0      // 1 = dexamethasone given
  use_RUX       = 0      // 1 = ruxolitinib given
  use_ANA       = 0      // 1 = anakinra given
  TOCI_dose_mg  = 560.0  // 8 mg/kg × 70 kg (mg) → converts to μg/mL in V1
  SILT_dose_mg  = 770.0  // 11 mg/kg × 70 kg
  DEX_dose_mg   = 10.0   // 10 mg IV per dose (mg) → ng/mL in V1
  RUX_dose_mg   = 10.0   // 10 mg BID → 20 mg/day
  ANA_dose_mg   = 100.0  // 100 mg SC QD (2 mg/kg ~50 kg pedi; or 100 mg adult)

$CMT
  // CAR-T dynamics (cells, relative units)
  CART_ACT      // Activated/expanding CAR-T cells
  CART_EXH      // Exhausted CAR-T (sink)
  TUMOR         // Tumor burden (relative antigen units)

  // Cytokines — T cell derived
  IFNg          // IFN-γ (pg/mL)
  IL2           // IL-2 (pg/mL)
  GMCSF         // GM-CSF (pg/mL)
  TNFa_T        // TNF-α (T cell) (pg/mL)

  // Macrophage activation cascade
  MAC_ACT       // Macrophage activation score (0-1 normalized)
  IL6           // IL-6 (pg/mL)
  IL1B          // IL-1β (pg/mL)
  TNFa_M        // TNF-α (macrophage) (pg/mL)

  // IL-6 downstream
  STAT3         // pSTAT3 activation (relative units)
  CRP           // CRP (mg/L)
  FERRITIN      // Ferritin (ng/mL)

  // Endothelial / vascular
  ENDO_ACT      // Endothelial activation score (0-1)
  VASC_PERM     // Vascular permeability (relative units)

  // CRS severity / organ damage
  CRS_SEV       // CRS severity score (0-4 scale, continuous)
  ORGAN_DMG     // Cumulative organ damage (0-1)

  // Drug PK compartments
  TOCI_C1       // Tocilizumab central (μg/mL)
  TOCI_C2       // Tocilizumab peripheral (μg/mL)
  SILT_C1       // Siltuximab central (μg/mL)
  DEX_C1        // Dexamethasone central (ng/mL)
  DEX_C2        // Dexamethasone peripheral (ng/mL)
  RUX_GUT       // Ruxolitinib gut depot
  RUX_C1        // Ruxolitinib central (ng/mL)
  RUX_C2        // Ruxolitinib peripheral (ng/mL)
  ANA_SC        // Anakinra SC depot
  ANA_C1        // Anakinra central (ng/mL)

$MAIN
  // Initial conditions
  CART_ACT_0 = 0.1;     // Small initial CAR-T expansion (just infused)
  CART_EXH_0 = 0.0;
  TUMOR_0    = TUMOR0;
  IFNg_0     = 5.0;     // Baseline IFN-γ (pg/mL)
  IL2_0      = 8.0;
  GMCSF_0    = 2.0;
  TNFa_T_0   = 3.0;
  MAC_ACT_0  = 0.05;    // Basal macrophage activity
  IL6_0      = 10.0;    // Baseline IL-6 (pg/mL)
  IL1B_0     = 3.0;
  TNFa_M_0   = 5.0;
  STAT3_0    = 0.1;
  CRP_0      = 3.0;     // Baseline CRP (mg/L)
  FERRITIN_0 = 150.0;   // Baseline ferritin (ng/mL)
  ENDO_ACT_0 = 0.05;
  VASC_PERM_0 = 0.1;
  CRS_SEV_0  = 0.0;
  ORGAN_DMG_0 = 0.0;

$ODE
  // ---- Drug PK ----
  // Tocilizumab (2-cmt IV)
  double ke_TOCI  = CL_TOCI / V1_TOCI;
  double k12_TOCI = Q_TOCI  / V1_TOCI;
  double k21_TOCI = Q_TOCI  / V2_TOCI;
  dxdt_TOCI_C1 = -(ke_TOCI + k12_TOCI) * TOCI_C1 + k21_TOCI * TOCI_C2;
  dxdt_TOCI_C2 = k12_TOCI * TOCI_C1 - k21_TOCI * TOCI_C2;

  // Siltuximab (1-cmt IV)
  double ke_SILT  = CL_SILT / V1_SILT;
  dxdt_SILT_C1 = -ke_SILT * SILT_C1;

  // Dexamethasone (2-cmt IV)
  double ke_DEX   = CL_DEX  / V1_DEX;
  double k12_DEX  = Q_DEX   / V1_DEX;
  double k21_DEX  = Q_DEX   / V2_DEX;
  dxdt_DEX_C1  = -(ke_DEX + k12_DEX) * DEX_C1 + k21_DEX * DEX_C2;
  dxdt_DEX_C2  = k12_DEX * DEX_C1 - k21_DEX * DEX_C2;

  // Ruxolitinib (2-cmt oral)
  double ke_RUX   = CL_RUX  / V1_RUX;
  double k12_RUX  = Q_RUX   / V1_RUX;
  double k21_RUX  = Q_RUX   / V2_RUX;
  dxdt_RUX_GUT = -ka_RUX * RUX_GUT;
  dxdt_RUX_C1  = ka_RUX * RUX_GUT - (ke_RUX + k12_RUX) * RUX_C1 + k21_RUX * RUX_C2;
  dxdt_RUX_C2  = k12_RUX * RUX_C1 - k21_RUX * RUX_C2;

  // Anakinra (1-cmt SC)
  double ke_ANA   = CL_ANA  / V1_ANA;
  dxdt_ANA_SC  = -ka_ANA * ANA_SC;
  dxdt_ANA_C1  = ka_ANA * ANA_SC - ke_ANA * ANA_C1;

  // ---- Drug PD — inhibition terms ----
  // Tocilizumab: IL-6R blockade
  double Inh_TOCI = (TOCI_C1 > 0) ? Emax_TOCI * TOCI_C1 / (EC50_TOCI + TOCI_C1) : 0.0;
  // Siltuximab: IL-6 neutralization
  double Inh_SILT = (SILT_C1 > 0) ? Emax_SILT * SILT_C1 / (EC50_SILT + SILT_C1) : 0.0;
  // Dexamethasone: cytokine suppression (dual: IL-6 and IFN-γ)
  double Inh_DEX_IL6  = (DEX_C1 > 0) ? Emax_DEX_IL6  * DEX_C1 / (EC50_DEX + DEX_C1) : 0.0;
  double Inh_DEX_IFNg = (DEX_C1 > 0) ? Emax_DEX_IFNg * DEX_C1 / (EC50_DEX + DEX_C1) : 0.0;
  // Ruxolitinib: JAK1/2 inhibition → downstream STAT3 and macrophage
  double Inh_RUX  = (RUX_C1 > 0) ? Emax_RUX_JAK * RUX_C1 / (EC50_RUX_JAK + RUX_C1) : 0.0;
  // Anakinra: IL-1β blockade
  double Inh_ANA  = (ANA_C1 > 0) ? Emax_ANA_IL1 * ANA_C1 / (EC50_ANA + ANA_C1) : 0.0;

  // Combined IL-6 signaling inhibition (TOCI + SILT + DEX + RUX)
  double Inh_IL6_sig = 1.0 - (1.0 - Inh_TOCI) * (1.0 - Inh_SILT) *
                              (1.0 - Inh_DEX_IL6) * (1.0 - Inh_RUX);

  // ---- CAR-T dynamics ----
  // Tumor drives activation; IL-2 drives expansion
  double IL2_prolif  = IL2 / (kIL2_prolif + IL2);
  double tumor_stim  = TUMOR / (0.5 + TUMOR);     // antigen availability
  double CART_prod   = kCART_act   * tumor_stim;
  double CART_expand = kCART_prolif * IL2_prolif * CART_ACT;
  double CART_loss   = (kCART_death + kCART_exh * CART_ACT) * CART_ACT;
  dxdt_CART_ACT = CART_prod + CART_expand - CART_loss;
  dxdt_CART_EXH = kCART_exh * CART_ACT * CART_ACT;

  // Tumor dynamics
  double tumor_kill  = kTumor_kill * CART_ACT * TUMOR;
  double tumor_grow  = kTumor_grow * TUMOR;
  dxdt_TUMOR = tumor_grow - tumor_kill;

  // ---- T cell–derived cytokine ODEs ----
  // IFN-γ: suppressed by DEX and RUX (via STAT1)
  double IFNg_prod = kIFNg_prod * CART_ACT * (1.0 - Inh_DEX_IFNg) * (1.0 - Inh_RUX * 0.5);
  dxdt_IFNg = IFNg_prod + 2.0 - dIFNg * IFNg;

  // IL-2: autocrine CAR-T growth signal
  dxdt_IL2  = kIL2_prod * CART_ACT + 2.0 - dIL2 * IL2;

  // GM-CSF: key macrophage activator
  dxdt_GMCSF = kGMCSF_prod * CART_ACT + 0.5 - dGMCSF * GMCSF;

  // TNF-α (T cell): suppressed by DEX
  dxdt_TNFa_T = kTNFa_T_prod * CART_ACT * (1.0 - Inh_DEX_IL6 * 0.7) + 1.0 - dTNFa_T * TNFa_T;

  // ---- Macrophage activation ----
  // IFN-γ and GM-CSF drive M1 macrophage polarization
  double IFNg_stim = IFNg / (EC50_IFNg_MAC + IFNg);
  double GMCSF_stim = GMCSF / (EC50_GMCSF + GMCSF);
  double MAC_stim  = 0.6 * IFNg_stim + 0.4 * GMCSF_stim;
  double MAC_prod  = kMAC_act * MAC_stim * (1.0 - MAC_ACT);  // saturates at 1
  double MAC_decay = kMAC_deact * MAC_ACT;
  dxdt_MAC_ACT = MAC_prod - MAC_decay;

  // IL-6: produced by activated macrophage; self-amplifying loop
  double IL6_amp   = 1.0 + IL6_feedback * MAC_ACT;
  double IL6_prod  = kIL6_prod * MAC_ACT * IL6_amp * (1.0 - Inh_DEX_IL6 * 0.8);
  dxdt_IL6 = IL6_prod + 3.0 - dIL6 * IL6;

  // IL-1β: NLRP3 inflammasome pathway (IFN-γ primed + TNF co-stimulation)
  double NLRP3_act = (IFNg / (IFNg + 100.0)) * (TNFa_T / (TNFa_T + 20.0));
  double IL1B_prod = kIL1B_prod * MAC_ACT * (1.0 + NLRP3_act) * (1.0 - Inh_ANA * 0.8);
  dxdt_IL1B = IL1B_prod + 0.5 - dIL1B * IL1B;

  // TNF-α (macrophage): suppressed by DEX
  double TNFaM_prod = kTNFa_M_prod * MAC_ACT * (1.0 - Inh_DEX_IL6 * 0.7);
  dxdt_TNFa_M = TNFaM_prod + 1.0 - dTNFa_M * TNFa_M;

  // ---- IL-6 downstream: STAT3 / CRP / Ferritin ----
  // STAT3: activated by IL-6 signaling; blocked by TOCI/SILT/RUX
  double IL6_sig   = (IL6 / (EC50_IL6_STAT3 + IL6)) * (1.0 - Inh_IL6_sig);
  dxdt_STAT3 = kSTAT3_act * IL6_sig * (1.0 - STAT3) - kSTAT3_deact * STAT3;

  // CRP: acute phase protein, IL-6 driven (hepatocyte response)
  dxdt_CRP = kCRP_prod * STAT3 * 50.0 - dCRP * CRP;

  // Ferritin: macrophage activation + IL-6 drive
  double Ferr_stim  = 1.0 + (IL6 / (IL6 + 100.0)) + (TNFa_M / (TNFa_M + 80.0));
  dxdt_FERRITIN = kFerr_prod * MAC_ACT * Ferr_stim * 500.0 - dFerr * FERRITIN;

  // ---- Endothelial activation ----
  double ENDO_stim = (IL6 / (EC50_IL6_ENDO + IL6)) + 0.5 * (TNFa_M / (50.0 + TNFa_M));
  dxdt_ENDO_ACT  = kENDO_act * ENDO_stim * (1.0 - ENDO_ACT) - kENDO_deact * ENDO_ACT;

  // Vascular permeability
  dxdt_VASC_PERM = kVASC_leak * ENDO_ACT - dVASC * VASC_PERM;

  // ---- CRS severity score (0-4 continuous) ----
  // Driven by: fever (IL-6→IL-1β), hypotension (VASC_PERM), IFN-γ
  double CRS_driver = (IFNg / (IFNg + 200.0)) * 2.0 +
                      (IL6  / (IL6  + 100.0)) * 1.5 +
                      VASC_PERM * 1.5;
  dxdt_CRS_SEV = kCRS_onset * CRS_driver * (4.0 - CRS_SEV) / 4.0 -
                 kCRS_resolve * CRS_SEV * (1.0 - Inh_IL6_sig * 0.8);

  // ---- Organ damage ----
  double dmg_driver = (CRS_SEV > 2.5) ? (CRS_SEV - 2.5) * 0.4 : 0.0;
  dxdt_ORGAN_DMG = kORGAN_dmg * dmg_driver * (1.0 - ORGAN_DMG) - dORGAN * ORGAN_DMG;

$TABLE
  // ---- Derived outputs ----
  // CRS grade (integer equivalent)
  int CRS_Grade = 0;
  if (CRS_SEV >= 0.5 && CRS_SEV < 1.5) CRS_Grade = 1;
  else if (CRS_SEV >= 1.5 && CRS_SEV < 2.5) CRS_Grade = 2;
  else if (CRS_SEV >= 2.5 && CRS_SEV < 3.5) CRS_Grade = 3;
  else if (CRS_SEV >= 3.5) CRS_Grade = 4;

  // Total IL-6 (pg/mL) — IL6 compartment
  double IL6_total = IL6;

  // Fever proxy (scale: 37°C baseline + IL-1β/IL-6 effect)
  double Temp_C = 37.0 + 2.5 * (IL1B / (IL1B + 20.0)) + 1.5 * (IL6 / (IL6 + 100.0));

  // Survival probability (exponential decline from organ damage)
  double Survival = exp(-2.0 * ORGAN_DMG);

  // CAR-T peak expansion (cells/μL surrogate)
  double CART_conc = CART_ACT * 100.0;  // scale to cells/μL

  // ICANS proxy (IL-6 CNS + IFN-γ BBB disruption)
  double ICANS_score = (IFNg / (IFNg + 500.0)) * 4.0 + (IL6 / (IL6 + 200.0)) * 3.0;
  double ICANS_Grade = (ICANS_score > 3.0) ? 4.0 :
                       (ICANS_score > 2.0) ? 3.0 :
                       (ICANS_score > 1.0) ? 2.0 :
                       (ICANS_score > 0.3) ? 1.0 : 0.0;

  // Biomarker panel
  double Total_TNFa = TNFa_T + TNFa_M;
  double Fibrinogen = 4.0 - 2.5 * (CRS_SEV / 4.0);  // ↓ with CRS severity (DIC)

  // CRS biomarker index (CARTOX-10 inspired)
  double CRS_Index = (IFNg / 500.0) + (IL6 / 200.0) + (CRP / 100.0) + (FERRITIN / 5000.0);

  double capture_IFNg   = IFNg;
  double capture_IL6    = IL6;
  double capture_IL1B   = IL1B;
  double capture_CRP    = CRP;
  double capture_FERRITIN = FERRITIN;
  double capture_TOCI   = TOCI_C1;
  double capture_DEX    = DEX_C1;
  double capture_RUX    = RUX_C1;
  double capture_ANA    = ANA_C1;
  double capture_SILT   = SILT_C1;

$CAPTURE
  CRS_Grade ICANS_Grade CRS_SEV ORGAN_DMG Survival
  capture_IFNg capture_IL6 capture_IL1B capture_CRP capture_FERRITIN
  CART_conc CART_ACT TUMOR Temp_C Fibrinogen CRS_Index
  VASC_PERM ENDO_ACT MAC_ACT STAT3
  capture_TOCI capture_DEX capture_RUX capture_ANA capture_SILT
  Total_TNFa ICANS_score GMCSF IL2
'

# ==============================================================================
# Build and compile model
# ==============================================================================
crs_mod <- mcode("crs_model", crs_model_code)

# ==============================================================================
# Helper: build dosing events
# ==============================================================================
make_doses <- function(scenario,
                       TOCI_dose = 560,   # mg → μg/mL in V1=3.1 L → 560000/3100 ≈ 180 μg/mL
                       SILT_dose = 770,   # mg
                       DEX_dose  = 10,    # mg → ng/mL in V1=25 L → 10000000/25000 ≈ 400 ng/mL
                       RUX_dose  = 10,    # mg per dose (BID)
                       ANA_dose  = 100) { # mg per dose (QD)

  # Convert doses to concentration equivalents (dose/volume × unit conversions)
  # TOCI: mg → μg/mL in 3.1 L = 3100 mL: 560 mg × 1000 μg/mg / 3100 mL = 180.6 μg/mL
  TOCI_conc <- TOCI_dose * 1000 / 3100      # μg/mL
  SILT_conc <- SILT_dose * 1000 / 5500      # μg/mL
  DEX_conc  <- DEX_dose  * 1e6  / 25000     # ng/mL (10 mg × 1e6 ng/mg / 25000 mL)
  RUX_amt   <- RUX_dose  * 1e6  / 72000     # ng/mL gut depot (10 mg BID)
  ANA_amt   <- ANA_dose  * 1e6  / 18000     # ng/mL SC depot

  ev <- NULL

  if (scenario == 1) {
    # No treatment
    ev <- ev(amt = 0, cmt = "TOCI_C1", time = 9999)  # placeholder

  } else if (scenario == 2) {
    # Tocilizumab alone (Grade 2 CRS, Day 3)
    ev <- ev(amt = TOCI_conc, cmt = "TOCI_C1", time = 3, tinf = 0.04) # 1h infusion
    # Second dose allowed at Day 5 if no improvement
    ev <- ev + ev(amt = TOCI_conc * 0.5, cmt = "TOCI_C1", time = 5, tinf = 0.04)

  } else if (scenario == 3) {
    # Tocilizumab + Dexamethasone (Grade 3 CRS, Day 2)
    ev_toci <- ev(amt = TOCI_conc, cmt = "TOCI_C1", time = 2, tinf = 0.04)
    ev_dex  <- ev(amt = DEX_conc,  cmt = "DEX_C1",  time = 2, tinf = 0.01)
    # Daily DEX for 5 days
    for (d in 3:6) {
      ev_dex <- ev_dex + ev(amt = DEX_conc * 0.8, cmt = "DEX_C1", time = d, tinf = 0.01)
    }
    # Repeat TOCI if needed
    ev_toci2 <- ev(amt = TOCI_conc * 0.5, cmt = "TOCI_C1", time = 4, tinf = 0.04)
    ev <- ev_toci + ev_dex + ev_toci2

  } else if (scenario == 4) {
    # Ruxolitinib + Dexamethasone (Refractory Grade 4 CRS, starting Day 2)
    ev_rux <- NULL
    for (d in seq(2, 20, by = 0.5)) {  # BID dosing
      ev_rux <- if (is.null(ev_rux)) ev(amt = RUX_amt, cmt = "RUX_GUT", time = d) else
                ev_rux + ev(amt = RUX_amt, cmt = "RUX_GUT", time = d)
    }
    ev_dex <- ev(amt = DEX_conc * 2, cmt = "DEX_C1", time = 2, tinf = 0.01)  # pulse
    for (d in 3:10) {
      ev_dex <- ev_dex + ev(amt = DEX_conc, cmt = "DEX_C1", time = d, tinf = 0.01)
    }
    # TOCI for immediate control
    ev_toci <- ev(amt = TOCI_conc, cmt = "TOCI_C1", time = 2, tinf = 0.04)
    ev <- ev_rux + ev_dex + ev_toci

  } else if (scenario == 5) {
    # Anakinra + Dexamethasone (MAS-type / NLRP3 dominant, QD SC)
    ev_ana <- NULL
    for (d in seq(2, 20, by = 1)) {
      ev_ana <- if (is.null(ev_ana)) ev(amt = ANA_amt, cmt = "ANA_SC", time = d) else
                ev_ana + ev(amt = ANA_amt, cmt = "ANA_SC", time = d)
    }
    ev_dex <- ev(amt = DEX_conc, cmt = "DEX_C1", time = 2, tinf = 0.01)
    for (d in 3:8) {
      ev_dex <- ev_dex + ev(amt = DEX_conc * 0.75, cmt = "DEX_C1", time = d, tinf = 0.01)
    }
    ev <- ev_ana + ev_dex
  }

  return(ev)
}

# ==============================================================================
# Run 5 scenarios
# ==============================================================================

scenarios <- list(
  list(id=1, name="Untreated",
       params=list(use_TOCI=0, use_DEX=0, use_RUX=0, use_ANA=0)),
  list(id=2, name="Tocilizumab (Grade 2 CRS)",
       params=list(use_TOCI=1, use_DEX=0, use_RUX=0, use_ANA=0)),
  list(id=3, name="Tocilizumab + DEX (Grade 3 CRS)",
       params=list(use_TOCI=1, use_DEX=1, use_RUX=0, use_ANA=0)),
  list(id=4, name="Ruxolitinib + DEX (Refractory/Grade 4)",
       params=list(use_TOCI=0, use_DEX=1, use_RUX=1, use_ANA=0)),
  list(id=5, name="Anakinra + DEX (MAS-type/Pediatric)",
       params=list(use_TOCI=0, use_DEX=1, use_RUX=0, use_ANA=1))
)

sim_results <- lapply(scenarios, function(sc) {
  dose_events <- make_doses(sc$id)
  mod_sc <- crs_mod %>% param(sc$params)
  out <- mod_sc %>%
    ev(dose_events) %>%
    mrgsim(start=0, end=28, delta=0.25) %>%
    as.data.frame()
  out$Scenario <- sc$name
  out$ScenarioID <- sc$id
  return(out)
})

all_results <- bind_rows(sim_results)
all_results$Scenario <- factor(all_results$Scenario,
  levels = sapply(scenarios, function(x) x$name))

# ==============================================================================
# Visualization
# ==============================================================================

scenario_colors <- c(
  "Untreated"                       = "#d32f2f",
  "Tocilizumab (Grade 2 CRS)"       = "#1976d2",
  "Tocilizumab + DEX (Grade 3 CRS)" = "#388e3c",
  "Ruxolitinib + DEX (Refractory/Grade 4)" = "#f57c00",
  "Anakinra + DEX (MAS-type/Pediatric)"    = "#7b1fa2"
)

# Figure 1: Cytokine time courses
fig1 <- all_results %>%
  select(time, Scenario, capture_IFNg, capture_IL6, capture_IL1B, capture_CRP) %>%
  pivot_longer(cols = c(capture_IFNg, capture_IL6, capture_IL1B, capture_CRP),
               names_to = "Biomarker", values_to = "Concentration") %>%
  mutate(Biomarker = recode(Biomarker,
    capture_IFNg = "IFN-γ (pg/mL)",
    capture_IL6  = "IL-6 (pg/mL)",
    capture_IL1B = "IL-1β (pg/mL)",
    capture_CRP  = "CRP (mg/L)"
  )) %>%
  ggplot(aes(x=time, y=Concentration, color=Scenario)) +
  geom_line(size=1.1) +
  facet_wrap(~Biomarker, scales="free_y", nrow=2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="CRS Cytokine Time Courses (5 Treatment Scenarios)",
       x="Time (days post CAR-T infusion)",
       y="Concentration",
       color="Treatment") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom", legend.direction="vertical")

print(fig1)

# Figure 2: CRS severity and organ damage
fig2 <- all_results %>%
  select(time, Scenario, CRS_SEV, ORGAN_DMG, Survival, capture_FERRITIN) %>%
  pivot_longer(cols=c(CRS_SEV, ORGAN_DMG, Survival, capture_FERRITIN),
               names_to="Endpoint", values_to="Value") %>%
  mutate(Endpoint = recode(Endpoint,
    CRS_SEV        = "CRS Severity (0-4)",
    ORGAN_DMG      = "Organ Damage (0-1)",
    Survival       = "Survival Probability",
    capture_FERRITIN = "Ferritin (ng/mL)"
  )) %>%
  ggplot(aes(x=time, y=Value, color=Scenario)) +
  geom_line(size=1.1) +
  facet_wrap(~Endpoint, scales="free_y", nrow=2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="CRS Clinical Endpoints",
       x="Time (days)", y="Value", color="Treatment") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom", legend.direction="vertical")

print(fig2)

# Figure 3: CAR-T expansion and drug PK
fig3_a <- all_results %>%
  filter(Scenario %in% c("Untreated","Tocilizumab + DEX (Grade 3 CRS)")) %>%
  ggplot(aes(x=time, y=CART_conc, color=Scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=scenario_colors) +
  labs(title="CAR-T Cell Expansion (cells/μL surrogate)",
       x="Time (days)", y="CAR-T cells/μL") +
  theme_bw(base_size=12)

fig3_b <- all_results %>%
  filter(ScenarioID == 2) %>%
  select(time, capture_TOCI) %>%
  ggplot(aes(x=time, y=capture_TOCI)) +
  geom_line(color="#1976d2", size=1.3) +
  labs(title="Tocilizumab PK (Scenario 2)",
       x="Time (days)", y="Tocilizumab (μg/mL)") +
  theme_bw(base_size=12)

print(fig3_a)
print(fig3_b)

# Figure 4: Sensitivity analysis — tumor burden at infusion vs CRS Grade
tumor_levels <- seq(0.2, 2.0, by=0.2)
sa_results <- lapply(tumor_levels, function(tb) {
  mod_sa <- crs_mod %>%
    param(list(TUMOR0=tb, use_TOCI=0, use_DEX=0, use_RUX=0, use_ANA=0))
  out <- mod_sa %>%
    mrgsim(start=0, end=14, delta=0.25) %>%
    as.data.frame()
  data.frame(
    TumorBurden = tb,
    MaxCRS = max(out$CRS_SEV),
    MaxIFNg = max(out$capture_IFNg),
    MaxIL6  = max(out$capture_IL6)
  )
})
sa_df <- bind_rows(sa_results)

fig4 <- ggplot(sa_df, aes(x=TumorBurden)) +
  geom_line(aes(y=MaxCRS * 25, color="Max CRS Grade (×25)"), size=1.3) +
  geom_line(aes(y=MaxIFNg / 20, color="Max IFN-γ/20 (pg/mL÷20)"), size=1.3, linetype=2) +
  geom_line(aes(y=MaxIL6 / 10, color="Max IL-6/10 (pg/mL÷10)"), size=1.3, linetype=3) +
  scale_color_manual(values=c("Max CRS Grade (×25)"="#d32f2f",
                               "Max IFN-γ/20 (pg/mL÷20)"="#ff7043",
                               "Max IL-6/10 (pg/mL÷10)"="#9c27b0")) +
  labs(title="Sensitivity Analysis: Tumor Burden vs CRS Severity (Untreated)",
       x="Initial Tumor Burden (relative units)",
       y="Scaled value",
       color="Output") +
  theme_bw(base_size=12)

print(fig4)

# ==============================================================================
# Summary statistics at Day 7 and Day 14
# ==============================================================================
cat("\n====== CRS Model Summary (Day 7) ======\n")
day7 <- all_results %>%
  filter(abs(time - 7) < 0.3) %>%
  group_by(Scenario) %>%
  summarise(
    CRS_Grade = round(mean(CRS_SEV), 2),
    IFNg_pg_mL = round(mean(capture_IFNg), 0),
    IL6_pg_mL  = round(mean(capture_IL6), 0),
    CRP_mg_L   = round(mean(capture_CRP), 0),
    Ferritin_ng_mL = round(mean(capture_FERRITIN), 0),
    Survival_pct = round(mean(Survival) * 100, 1),
    .groups = "drop"
  )
print(as.data.frame(day7))

cat("\n====== CRS Model Summary (Day 14) ======\n")
day14 <- all_results %>%
  filter(abs(time - 14) < 0.3) %>%
  group_by(Scenario) %>%
  summarise(
    CRS_Grade = round(mean(CRS_SEV), 2),
    IFNg_pg_mL = round(mean(capture_IFNg), 0),
    IL6_pg_mL  = round(mean(capture_IL6), 0),
    Ferritin_ng_mL = round(mean(capture_FERRITIN), 0),
    Survival_pct = round(mean(Survival) * 100, 1),
    .groups = "drop"
  )
print(as.data.frame(day14))

# ==============================================================================
# END: Run the Shiny app instead for interactive exploration:
# shiny::runApp("crs_shiny_app.R")
# ==============================================================================
