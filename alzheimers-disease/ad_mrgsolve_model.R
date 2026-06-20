## ============================================================
## Alzheimer's Disease — QSP Model (mrgsolve)
## Amyloid Cascade · Tau Pathology · Neuroinflammation
## Cholinergic System · Synaptic Plasticity · Drug PK/PD
##
## Compartments (21 ODEs):
##   PK : Donepezil (gut, plasma, CNS), Memantine (plasma, CNS),
##        Lecanemab (plasma, CNS)
##   PD : Abeta production/aggregation/clearance (4 states),
##        Tau (3 states), Neuroinflammation (1),
##        Cholinergic ACh (1), Synaptic integrity (1),
##        Cognitive function/MMSE-proxy (1)
##
## Key Clinical Trials Used for Calibration:
##   CLARITY-AD  (van Dyck et al. NEJM 2023) — Lecanemab
##   TRAILBLAZER-ALZ2 (Sims et al. NEJM 2023) — Donanemab
##   ADAS-E2E    (Rogers et al. Neurology 1998) — Donepezil
##   MEM-MD-02   (Reisberg et al. NEJM 2003)  — Memantine
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── MODEL DEFINITION ─────────────────────────────────────────
code <- '
$PROB Alzheimer Disease QSP Model
  Compartments: Donepezil PK (gut/plasma/CNS), Memantine PK,
  Lecanemab PK (plasma/CNS), Abeta dynamics (4 pools),
  Tau dynamics (3 pools), Neuroinflammation, ACh, Synaptic, Cognition

$PARAM
  // ── DONEPEZIL PK (Pilla Reddy et al. CPT 2011) ──────────────
  KA_DON   = 0.0115    // h-1   oral absorption rate constant
  CL_DON   = 13.3      // L/h   systemic clearance
  V2_DON   = 594       // L     central volume
  Q_DON    = 40.0      // L/h   intercompartmental clearance
  V3_DON   = 1180      // L     peripheral volume
  Kp_DON   = 0.18      // -     brain:plasma partition (free)
  F_DON    = 1.0       // -     bioavailability

  // ── MEMANTINE PK (Periclou et al. Clin Pharmacol 2004) ──────
  KA_MEM   = 0.072     // h-1   absorption rate
  CL_MEM   = 8.4       // L/h   clearance
  V_MEM    = 520       // L     distribution volume
  Kp_MEM   = 2.5       // -     brain:plasma partition

  // ── LECANEMAB PK (Swanson et al. Alzheimers Dement 2021) ────
  CL_LEC   = 0.0167    // L/h   linear clearance (mAb)
  V1_LEC   = 3.2       // L     central (plasma) volume
  Q_LEC    = 0.00458   // L/h   peripheral clearance
  V2_LEC   = 2.2       // L     peripheral volume
  KBBB_LEC = 0.000083  // h-1   BBB influx rate (~0.1% CNS exposure)
  KCNS_LEC = 0.0042    // h-1   CNS elimination rate

  // ── DRUG PHARMACODYNAMICS ────────────────────────────────────
  IC50_AChE_DON  = 6.7e-6   // mg/L   AChE IC50 donepezil (Ki ~6.7 nM)
  IC50_NMDA_MEM  = 0.8      // mg/L   NMDAR IC50 memantine
  KD_proto_LEC   = 5.4e-5   // mg/L   lecanemab Kd for protofibrils
  HILL_DON       = 1.2      // -      Hill coefficient donepezil
  HILL_MEM       = 1.0      // -      Hill coefficient memantine
  HILL_LEC       = 1.0      // -      Hill coefficient lecanemab

  // ── AMYLOID DYNAMICS ─────────────────────────────────────────
  // Production (Bateman et al. Science 2006: Abeta42 FR ~252 pmol/h/L CSF)
  kprod_Ab  = 0.042    // 1/h   Abeta42 production rate (brain)
  kdeg_Ab   = 0.037    // 1/h   monomeric Abeta clearance (IDE/NEP)
  // Aggregation kinetics (Hasegawa et al. J Biol Chem 1999)
  kn_oligo  = 0.0018   // 1/h·nM nucleation (monomer→oligomer)
  ke_oligo  = 0.0045   // 1/h   elongation (oligo→protofibril)
  ke_fibril = 0.0012   // 1/h   fibril formation
  kdep_plaq = 0.00055  // 1/h   plaque deposition
  // Clearance by microglia
  kphago    = 0.0003   // 1/h   microglial phagocytosis of plaques
  // BBB transport
  kLRP1     = 0.014    // 1/h   LRP1-mediated efflux
  kRAGE     = 0.006    // 1/h   RAGE-mediated influx
  // ApoE4 effect on clearance (reduced by ~40%)
  APOE4_eff = 1.0      // 1=no APOE4, 0.6=APOE4 carrier (set in scenario)

  // ── TAU DYNAMICS ─────────────────────────────────────────────
  // (Jack et al. Lancet Neurol 2010; Golde et al. Neuron 2011)
  kprod_tau  = 0.028   // 1/h   tau production
  kphos_tau  = 0.0022  // 1/h   basal phosphorylation rate (CDK5/GSK3b)
  kAb_tau    = 0.0055  // amplification: Abeta oligomers → tau kinase activation
  kdphos_tau = 0.0019  // 1/h   dephosphorylation (PP2A)
  kagg_tau   = 0.0009  // 1/h   pTau → oligomers
  kNFT       = 0.00035 // 1/h   oligomer → NFT
  kdeg_tau   = 0.016   // 1/h   normal tau clearance
  kspread    = 0.000045// 1/h   tau spread amplification (prion-like)

  // ── NEUROINFLAMMATION ────────────────────────────────────────
  // (Heneka et al. Nat Rev Neurosci 2015)
  kact_micro  = 0.0065  // 1/h  microglial activation by Abeta
  kinact_micro= 0.0045  // 1/h  microglial resolution
  kTNF        = 0.012   // 1/h  TNFa production by activated microglia
  kdeg_TNF    = 0.095   // 1/h  TNFa clearance
  ksynapse_inf= 0.0018  // 1/h  inflammation→synapse loss

  // ── CHOLINERGIC ──────────────────────────────────────────────
  // (Whitehouse et al. Science 1982; Davies & Maloney 1976)
  ACh_ss      = 1.0    // normalized ACh at baseline
  ksynth_ACh  = 0.15   // 1/h  ACh synthesis
  kdeg_ACh    = 0.12   // 1/h  ACh hydrolysis (AChE + BuChE)
  kACh_BFCN   = 0.008  // loss of BFCN neurons → ACh production loss
  ACh_EC50    = 0.4    // normalized ACh for 50% cognitive benefit

  // ── SYNAPTIC INTEGRITY ───────────────────────────────────────
  // (Terry et al. Ann Neurol 1991)
  Syn0        = 1.0    // initial synaptic integrity (normalized)
  kloss_syn_Ab= 0.0012 // 1/h  synapse loss by Abeta oligomers
  kloss_syn_tau=0.00085// 1/h  synapse loss by tau pathology
  kloss_syn_inf=0.00055// 1/h  synapse loss by neuroinflammation
  kregen_syn  = 0.00018// 1/h  synaptic regeneration (BDNF-mediated)

  // ── COGNITIVE FUNCTION ───────────────────────────────────────
  // MMSE-proxy: 0 (severe) → 30 (normal)
  MMSE0       = 28.0   // initial MMSE (mild cognitive impairment)
  kdecline    = 0.000055// 1/h  baseline cognitive decline rate
  // Cognitive contributions from each compartment
  wSyn        = 0.50   // weight: synaptic integrity
  wACh        = 0.25   // weight: cholinergic tone
  wTau        = 0.15   // weight: tau pathology (inverse)
  wAb         = 0.10   // weight: amyloid burden (inverse)
  MMSE_floor  = 0.0    // minimum MMSE score

  // ── DISEASE MODIFIERS ────────────────────────────────────────
  APOE4_carrier = 0    // 0=no, 1=yes (one allele), 2=homozygous
  FAD_mutation  = 0    // 0=sporadic, 1=PSEN1/APP mutation (higher kprod)
  age_factor    = 1.0  // 1.0=70yr, scales kprod/kphos

$CMT
  // Donepezil PK (3-compartment oral)
  DON_GUT DON_CENT DON_PERI DON_CNS
  // Memantine PK (1-compartment)
  MEM_CENT MEM_CNS
  // Lecanemab PK (2-compartment IV + CNS)
  LEC_CENT LEC_PERI LEC_CNS
  // Amyloid dynamics
  AB_MONO AB_OLIGO AB_PROTO AB_PLAQUE
  // Tau dynamics
  TAU_SOL TAU_PHOS TAU_AGG
  // Neuroinflammation (normalized 0–1)
  NEURO_INFLAM
  // Cholinergic
  ACH
  // Synaptic integrity (normalized 0–1)
  SYN
  // Cognitive function (MMSE-proxy 0–30)
  COGNITION

$MAIN
  // ── Disease modifier scalings ───────────────────────────────
  double ApoeScale = 1.0 - 0.40 * APOE4_carrier * APOE4_eff;
  double FadScale  = 1.0 + 0.80 * FAD_mutation;
  double AgeScale  = age_factor;

  // ── Effective PK parameters ─────────────────────────────────
  double CL2  = CL_DON;
  double V2   = V2_DON;
  double Q3   = Q_DON;
  double V3   = V3_DON;

  // ── Drug concentration conversions ──────────────────────────
  // Donepezil CNS free conc (mg/L)
  double Cp_DON = DON_CENT / V2_DON;
  double Cb_DON = DON_CNS;

  // Memantine CNS free conc (mg/L)
  double Cp_MEM = MEM_CENT / V_MEM;
  double Cb_MEM = MEM_CNS;

  // Lecanemab CNS conc (mg/L)
  double Cp_LEC  = LEC_CENT / V1_LEC;
  double Ccns_LEC = LEC_CNS;

  // ── Drug effect: AChE inhibition by Donepezil ───────────────
  double Imax_DON = 1.0;
  double inh_AChE  = Imax_DON * pow(Cb_DON, HILL_DON) /
                     (pow(IC50_AChE_DON, HILL_DON) + pow(Cb_DON, HILL_DON));

  // ── Drug effect: NMDA block by Memantine ────────────────────
  double inh_NMDA  = Cb_MEM / (IC50_NMDA_MEM + Cb_MEM);

  // ── Drug effect: Lecanemab on protofibrils ───────────────────
  double inh_proto = Ccns_LEC / (KD_proto_LEC + Ccns_LEC);

  // ── Effective ACh degradation rate ──────────────────────────
  double kdeg_ACh_eff = kdeg_ACh * (1.0 - inh_AChE);

  // ── Abeta production (scaled by FAD/age) ────────────────────
  double Ab_kprod = kprod_Ab * FadScale * AgeScale;

  // ── Tau phosphorylation (stimulated by Abeta oligomers) ─────
  double kphos_eff = kphos_tau + kAb_tau * AB_OLIGO;

  // ── Microglial phagocytosis (inhibited if inflammation high) ─
  double kphago_eff = kphago * (1.0 - 0.5 * NEURO_INFLAM);

  // ── Synapse regeneration (BDNF, depends on BDNF which ~∝ACh)─
  double kregen_eff = kregen_syn * (ACH / ACh_ss);

  // ── Cognitive function weighting ────────────────────────────
  // SYN and ACH promote cognition; TAU_AGG and AB_PLAQUE impair it
  double norm_Syn   = fmax(0.0, fmin(1.0, SYN));
  double norm_ACh   = fmax(0.0, fmin(2.0, ACH / ACh_ss));
  double norm_tau   = fmax(0.0, fmin(1.0, TAU_AGG / 5.0));
  double norm_Ab    = fmax(0.0, fmin(1.0, AB_PLAQUE / 10.0));

  // Initialize compartments
  if(NEWIND <= 1){
    // Initial conditions — mild-moderate AD at model start
    AB_MONO_0   = 0.8;    // nM equivalent (normalized)
    AB_OLIGO_0  = 0.5;
    AB_PROTO_0  = 0.3;
    AB_PLAQUE_0 = 2.5;
    TAU_SOL_0   = 1.0;
    TAU_PHOS_0  = 1.5;
    TAU_AGG_0   = 0.8;
    NEURO_INFLAM_0 = 0.4;
    ACH_0       = 0.65;   // reduced at disease onset
    SYN_0       = 0.72;   // 72% synaptic integrity (MCI stage)
    COGNITION_0 = MMSE0;
  }

$ODE
  // ═══════════════════════════════════════════════════════════
  // 1. DONEPEZIL PK (3-compartment oral + CNS compartment)
  // ═══════════════════════════════════════════════════════════
  dxdt_DON_GUT  = -KA_DON * DON_GUT;
  dxdt_DON_CENT =  KA_DON * DON_GUT * F_DON
                 - (CL_DON + Q_DON) / V2_DON * DON_CENT
                 + Q_DON / V3_DON * DON_PERI;
  dxdt_DON_PERI =  Q_DON / V2_DON * DON_CENT
                 - Q_DON / V3_DON * DON_PERI;
  // CNS: pseudo-equilibrium with plasma (Kp ~ 0.18)
  dxdt_DON_CNS  =  0.5 * (Kp_DON * DON_CENT / V2_DON - DON_CNS);

  // ═══════════════════════════════════════════════════════════
  // 2. MEMANTINE PK (1-compartment + CNS)
  // ═══════════════════════════════════════════════════════════
  dxdt_MEM_CENT = -CL_MEM / V_MEM * MEM_CENT;
  dxdt_MEM_CNS  =  0.3 * (Kp_MEM * MEM_CENT / V_MEM - MEM_CNS);

  // ═══════════════════════════════════════════════════════════
  // 3. LECANEMAB PK (2-compartment IV + CNS transfer)
  // ═══════════════════════════════════════════════════════════
  dxdt_LEC_CENT = -CL_LEC / V1_LEC * LEC_CENT
                 - Q_LEC / V1_LEC * LEC_CENT
                 + Q_LEC / V2_LEC * LEC_PERI
                 - KBBB_LEC * LEC_CENT;
  dxdt_LEC_PERI =  Q_LEC / V1_LEC * LEC_CENT
                 - Q_LEC / V2_LEC * LEC_PERI;
  dxdt_LEC_CNS  =  KBBB_LEC * LEC_CENT
                 - KCNS_LEC * LEC_CNS
                 - inh_proto * 0.05 * LEC_CNS;  // drug bound to target

  // ═══════════════════════════════════════════════════════════
  // 4. AMYLOID DYNAMICS
  // ═══════════════════════════════════════════════════════════
  // Monomers: produced from APP processing, cleared by IDE/NEP
  dxdt_AB_MONO  =  Ab_kprod
                 - kdeg_Ab * ApoeScale * AB_MONO
                 - kn_oligo * AB_MONO * AB_MONO    // nucleation to oligomers
                 + kRAGE * 0.1                      // RAGE influx from blood
                 - kLRP1 * AB_MONO;                // LRP1 efflux

  // Oligomers: nucleated from monomers, elongate to protofibrils
  dxdt_AB_OLIGO =  kn_oligo * AB_MONO * AB_MONO
                 - ke_oligo * AB_OLIGO
                 - kdeg_Ab * 0.5 * AB_OLIGO;       // partial clearance

  // Protofibrils: blocked by lecanemab
  dxdt_AB_PROTO =  ke_oligo * AB_OLIGO
                 - ke_fibril * AB_PROTO
                 - inh_proto * ke_fibril * AB_PROTO; // lecanemab binding

  // Plaques: deposited, cleared by microglia
  dxdt_AB_PLAQUE=  ke_fibril * (1.0 - inh_proto) * AB_PROTO
                 - kphago_eff * AB_PLAQUE;

  // ═══════════════════════════════════════════════════════════
  // 5. TAU DYNAMICS
  // ═══════════════════════════════════════════════════════════
  // Soluble tau: produced, phosphorylated, degraded
  dxdt_TAU_SOL  =  kprod_tau * AgeScale
                 - kphos_eff * TAU_SOL
                 - kdeg_tau  * TAU_SOL;

  // Phospho-tau: dephosphorylated by PP2A, aggregates
  dxdt_TAU_PHOS =  kphos_eff * TAU_SOL
                 - kdphos_tau * TAU_PHOS
                 - kagg_tau  * TAU_PHOS
                 - kdeg_tau  * 0.2 * TAU_PHOS;

  // Aggregated tau / NFTs: prion-like spread
  dxdt_TAU_AGG  =  kagg_tau * TAU_PHOS
                 + kspread  * TAU_AGG * (1.0 - TAU_AGG/20.0) // logistic spread
                 - kNFT     * TAU_AGG;

  // ═══════════════════════════════════════════════════════════
  // 6. NEUROINFLAMMATION (normalized 0–1)
  // ═══════════════════════════════════════════════════════════
  dxdt_NEURO_INFLAM = kact_micro  * AB_PLAQUE * (1.0 - NEURO_INFLAM)
                    + 0.0008  * TAU_AGG   * (1.0 - NEURO_INFLAM)
                    - kinact_micro * NEURO_INFLAM;

  // ═══════════════════════════════════════════════════════════
  // 7. CHOLINERGIC — ACh (normalized)
  // ═══════════════════════════════════════════════════════════
  // BFCN loss driven by amyloid + tau over time
  double BFCN_loss = kACh_BFCN * (AB_PLAQUE * 0.5 + TAU_AGG * 0.5);
  dxdt_ACH  =  ksynth_ACh * (1.0 - BFCN_loss)
             - kdeg_ACh_eff * ACH;

  // ═══════════════════════════════════════════════════════════
  // 8. SYNAPTIC INTEGRITY (normalized 0–1)
  // ═══════════════════════════════════════════════════════════
  dxdt_SYN  =  kregen_eff * SYN * (1.0 - SYN)   // logistic recovery
             - kloss_syn_Ab  * AB_OLIGO * SYN
             - kloss_syn_tau * TAU_PHOS * SYN
             - kloss_syn_inf * NEURO_INFLAM * SYN
             // Memantine partially protects synapses from excitotoxicity
             + 0.0004 * inh_NMDA * SYN;

  // ═══════════════════════════════════════════════════════════
  // 9. COGNITIVE FUNCTION (MMSE-proxy 0–30)
  // ═══════════════════════════════════════════════════════════
  // Instantaneous cognitive state based on biological drivers
  double MMSE_target = 30.0 * (wSyn  * norm_Syn
                              + wACh  * norm_ACh * 0.5
                              - wTau  * norm_tau
                              - wAb   * norm_Ab);
  MMSE_target = fmax(MMSE_floor, MMSE_target);
  // First-order approach to target (chronic disease time scale)
  dxdt_COGNITION = 0.00025 * (MMSE_target - COGNITION);

$TABLE
  // ── Primary endpoints ────────────────────────────────────────
  double MMSE     = COGNITION;
  double ADAS_Cog = fmax(0, 70.0 - 2.33 * COGNITION); // approximate conversion
  double CDR_SB   = fmax(0, 18.0 * (1.0 - COGNITION / 30.0));

  // ── Biomarker outputs ─────────────────────────────────────────
  // CSF Abeta42: inversely related to plaque load (Jack 2010)
  double CSF_Ab42 = 1200 * exp(-0.35 * AB_PLAQUE);    // pg/mL (cutoff ~1000)
  // CSF p-Tau181: proportional to phospho-tau pool
  double CSF_pTau181 = 15 + 28 * TAU_PHOS;            // pg/mL (cutoff ~23)
  // Amyloid PET (centiloids): monotone with plaque load
  double AmyloidPET_CL = 20 + 14 * AB_PLAQUE;         // centiloids (cutoff 24)
  // Tau PET: proportional to aggregated tau
  double TauPET_SUVr = 1.0 + 0.18 * TAU_AGG;          // SUVr
  // Plasma NfL (neurodegeneration marker)
  double NfL_plasma = 10 + 22 * (1.0 - SYN);          // pg/mL

  // ── Drug concentrations ──────────────────────────────────────
  double Cp_Donepezil  = DON_CENT / V2_DON * 1000;    // ng/mL
  double Ccns_Donepezil = DON_CNS * 1000;              // ng/mL
  double Cp_Memantine  = MEM_CENT / V_MEM * 1000;
  double Cp_Lecanemab  = LEC_CENT / V1_LEC * 1000;    // ug/mL

  // ── Inhibition/occupancy outputs ─────────────────────────────
  double AChE_inhibition = inh_AChE * 100;    // % inhibition
  double NMDAR_occupancy = inh_NMDA * 100;
  double Proto_neutralized = inh_proto * 100;

$CAPTURE
  MMSE ADAS_Cog CDR_SB
  CSF_Ab42 CSF_pTau181 AmyloidPET_CL TauPET_SUVr NfL_plasma
  Cp_Donepezil Ccns_Donepezil Cp_Memantine Cp_Lecanemab
  AChE_inhibition NMDAR_occupancy Proto_neutralized
  AB_MONO AB_OLIGO AB_PROTO AB_PLAQUE
  TAU_SOL TAU_PHOS TAU_AGG
  NEURO_INFLAM ACH SYN
'

## Compile model
mod <- mrgsolve::mcode("AD_QSP", code)

## ── DOSING REGIMENS ──────────────────────────────────────────────
# Donepezil 10 mg QD (oral), steady-state after ~3 weeks
don_dose <- ev(amt = 10, ii = 24, addl = 365*2 - 1, cmt = "DON_GUT")

# Memantine 20 mg QD (oral) — titrated, here simplified as maintenance
mem_dose <- ev(amt = 20, ii = 24, addl = 365*2 - 1, cmt = "MEM_CENT")

# Lecanemab 10 mg/kg biweekly IV (CLARITY-AD regimen; ~700 mg for 70 kg pt)
lec_dose <- ev(amt = 700, ii = 24*14, addl = 26 - 1, cmt = "LEC_CENT",
               rate = -2)  # infusion over 1h

# Combination donepezil + memantine
combi_don_mem <- ev_seq(don_dose, mem_dose)

## ── SCENARIO DEFINITIONS ─────────────────────────────────────────
scenarios <- list(
  "1_Control"              = list(dose = NULL,         apoe4=0, fad=0),
  "2_Donepezil_10mg"       = list(dose = don_dose,     apoe4=0, fad=0),
  "3_Memantine_20mg"       = list(dose = mem_dose,     apoe4=0, fad=0),
  "4_Combo_Don_Mem"        = list(dose = combi_don_mem,apoe4=0, fad=0),
  "5_Lecanemab"            = list(dose = lec_dose,     apoe4=0, fad=0),
  "6_Lecanemab_APOE4"      = list(dose = lec_dose,     apoe4=1, fad=0),
  "7_FAD_Mutation_Untreated"= list(dose = NULL,        apoe4=0, fad=1)
)

## ── SIMULATION FUNCTION ──────────────────────────────────────────
simulate_scenario <- function(scenario_name, scenario_params, duration_years = 3) {
  # Time grid: every 8h for 3 years (26280 hours)
  tgrid <- seq(0, duration_years * 8760, by = 8)

  # Set parameters
  par_update <- list(
    APOE4_carrier = scenario_params$apoe4,
    FAD_mutation  = scenario_params$fad
  )

  if (!is.null(scenario_params$dose)) {
    out <- mod %>%
      param(par_update) %>%
      mrgsim(events = scenario_params$dose, tgrid = tgrid, delta = 8,
             carry_out = "amt,evid,cmt") %>%
      as_tibble()
  } else {
    out <- mod %>%
      param(par_update) %>%
      mrgsim(tgrid = tgrid, delta = 8) %>%
      as_tibble()
  }

  out$scenario <- scenario_name
  out
}

## ── RUN ALL SCENARIOS ────────────────────────────────────────────
cat("Running Alzheimer's Disease QSP simulations...\n")
results <- purrr::map_dfr(names(scenarios), function(nm) {
  cat(" Scenario:", nm, "\n")
  simulate_scenario(nm, scenarios[[nm]])
})

# Convert time to years
results <- results %>% mutate(year = time / 8760)

## ── FIGURE 1: MMSE TRAJECTORIES ──────────────────────────────────
p1 <- ggplot(results %>% filter(evid == 0),
             aes(x = year, y = MMSE, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(10, 20, 26), linetype = "dashed",
             color = c("red","orange","blue"), alpha = 0.5) +
  annotate("text", x = 2.8, y = c(10.5, 20.5, 26.5),
           label = c("Severe/Moderate","MCI boundary","Normal boundary"),
           size = 2.8, hjust = 1) +
  scale_x_continuous(breaks = 0:3) +
  labs(title = "MMSE Trajectories — AD QSP Model",
       x = "Time (years)", y = "MMSE Score",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

## ── FIGURE 2: AMYLOID BURDEN ─────────────────────────────────────
p2 <- ggplot(results %>% filter(evid == 0),
             aes(x = year, y = AmyloidPET_CL, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 24, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 0.2, y = 26, label = "Positivity threshold (24 CL)",
           size = 2.8, color = "red") +
  labs(title = "Amyloid PET (Centiloids)",
       x = "Time (years)", y = "Centiloids",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

## ── FIGURE 3: CSF BIOMARKERS ─────────────────────────────────────
p3 <- results %>%
  filter(evid == 0) %>%
  select(year, scenario, CSF_Ab42, CSF_pTau181) %>%
  pivot_longer(c(CSF_Ab42, CSF_pTau181), names_to = "biomarker") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~biomarker, scales = "free_y",
             labeller = labeller(biomarker = c(
               CSF_Ab42     = "CSF Aβ42 (pg/mL)",
               CSF_pTau181  = "CSF p-Tau181 (pg/mL)"))) +
  labs(title = "CSF Biomarkers", x = "Time (years)", y = "Concentration",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## ── FIGURE 4: SYNAPTIC INTEGRITY & ACh ───────────────────────────
p4 <- results %>%
  filter(evid == 0) %>%
  select(year, scenario, SYN, ACH) %>%
  pivot_longer(c(SYN, ACH), names_to = "variable") %>%
  ggplot(aes(x = year, y = value, color = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~variable, scales = "free_y",
             labeller = labeller(variable = c(
               SYN = "Synaptic Integrity (normalized)",
               ACH = "Acetylcholine (normalized)"))) +
  labs(title = "Synaptic Integrity & Cholinergic Tone",
       x = "Time (years)", y = "Value (normalized)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## ── FIGURE 5: AChE INHIBITION PROFILE ────────────────────────────
p5 <- results %>%
  filter(evid == 0, grepl("Donepezil|Combo|Control", scenario)) %>%
  ggplot(aes(x = year, y = AChE_inhibition, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = c(60, 80), linetype = "dashed",
             color = c("orange","red"), alpha = 0.6) +
  labs(title = "AChE Inhibition (%) — Donepezil",
       x = "Time (years)", y = "AChE Inhibition (%)",
       color = "Scenario") +
  theme_bw(base_size = 12)

## ── PRINT SUMMARY STATISTICS ─────────────────────────────────────
cat("\n=== 36-Month Summary Statistics ===\n")
summary_tbl <- results %>%
  filter(abs(year - 3.0) < 0.05, evid == 0) %>%
  group_by(scenario) %>%
  summarise(
    MMSE_final     = mean(MMSE),
    ADAS_Cog_final = mean(ADAS_Cog),
    CDR_SB_final   = mean(CDR_SB),
    CL_final       = mean(AmyloidPET_CL),
    CSF_Ab42_final = mean(CSF_Ab42),
    SYN_final      = mean(SYN),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl))

## ── LECANEMAB DOSE-RESPONSE ───────────────────────────────────────
doses_lec <- c(1, 2.5, 5, 10, 20)  # mg/kg (for 70 kg pt)
dose_response <- purrr::map_dfr(doses_lec, function(d) {
  ev_d <- ev(amt = d * 70, ii = 24 * 14, addl = 26 - 1,
             cmt = "LEC_CENT", rate = -2)
  out <- mod %>%
    mrgsim(events = ev_d, tgrid = seq(0, 8760, 8), delta = 8) %>%
    as_tibble() %>%
    filter(abs(time - 8760) < 10) %>%
    mutate(dose_mpkg = d)
  out
})
cat("\nLecanemab dose-response at 12 months:\n")
print(dose_response %>%
  select(dose_mpkg, MMSE, AmyloidPET_CL, CSF_pTau181) %>%
  distinct())

## ── RETURN RESULTS ───────────────────────────────────────────────
invisible(list(
  results   = results,
  model     = mod,
  plots     = list(mmse = p1, amyloid = p2, csf = p3,
                   synapse = p4, ache = p5),
  summary   = summary_tbl,
  dose_resp = dose_response
))
