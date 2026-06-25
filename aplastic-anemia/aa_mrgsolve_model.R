## ============================================================
##  Aplastic Anemia (AA) — QSP Model (mrgsolve)
##  19 ODE Compartments: Disease PD (13) + Drug PK (6)
##
##  Key references for parameter calibration:
##  - Young NS et al. NEJM 2012 (hATG vs rATG, 6-mo response 68% vs 37%)
##  - Townsley DM et al. NEJM 2017 (hATG+CsA+EPAG: 6-mo CR ~58%)
##  - Scheinberg P et al. Blood 2012 (IST review)
##  - Desmond R et al. J Clin Invest 2013 (EPAG restores multilineage)
##  - Townsley DM et al. NEJM 2016 (Danazol in telomere diseases)
##
##  Clinical thresholds for Severe AA (sAA):
##    ANC  < 0.5 × 10⁹/L
##    PLT  < 20  × 10⁹/L
##    Hgb  < 8   g/dL (or reticulocytes < 20×10⁹/L with both above)
##
##  Complete Response (CR): ANC>1.0, PLT>100, Hgb>10
##  Partial Response (PR):  Transfusion-independent; no longer meets sAA
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model definition ----------------------------------------
code <- '
$PROB
  Aplastic Anemia QSP Model
  Disease PD compartments (13):
    CD8  - Autoreactive CD8+ T cells (cells/uL)
    TREG - Regulatory T cells (cells/uL)
    CYTO - Effective cytokine index (IFN-g/TNF-a, normalized)
    HSC  - HSC pool (AU, 100 = normal)
    PROG - Myeloid progenitor pool (AU, 100 = normal)
    CFU_E- Erythroid progenitors (AU)
    RETIC- Reticulocytes (x10^9/L)
    RBC  - Hemoglobin equivalent (g/dL)
    MKP  - Megakaryocyte progenitors (AU)
    PLT  - Platelets (x10^9/L)
    NEUP - Neutrophil precursors (AU)
    NEU  - Mature neutrophils (x10^9/L)
    TELO - Telomere length (relative, 1.0 = normal)
  Drug PK compartments (6):
    CATG - ATG central compartment (ug/mL)
    PATG - ATG peripheral compartment (ug/mL)
    CGUT - CsA gut compartment (mg/L)
    CCSA - CsA central compartment (ng/mL equivalent)
    EGUT - Eltrombopag gut (mg)
    CEPG - Eltrombopag central (ug/mL)

$PARAM @annotated
  // --- Disease parameters ---
  CD8_0    : 800   : Baseline CD8+ autoreactive T cells (cells/uL)
  TREG_0   : 100   : Baseline Treg cells (cells/uL)
  CYTO_0   : 1.0   : Baseline cytokine index (normalized)
  HSC_0    : 100   : Baseline HSC pool (AU = 100%)
  kCD8p    : 0.03  : CD8 proliferation rate (day^-1)
  kCD8d    : 0.03  : CD8 death rate (day^-1)
  kTRp     : 3.0   : Treg production rate (cells/uL/day)
  kTRd     : 0.03  : Treg death rate (day^-1)
  kCYp     : 0.5   : Cytokine production rate by CD8
  kCYd     : 2.0   : Cytokine degradation rate (day^-1)
  kHSCsr   : 0.08  : HSC self-renewal rate (day^-1)
  kHSCd    : 0.02  : HSC baseline death rate (day^-1)
  kHSCdf   : 0.05  : HSC differentiation rate (day^-1)
  kCYkill  : 0.25  : Cytokine-driven HSC kill rate coefficient
  kPROGdf  : 0.15  : Progenitor differentiation rate (day^-1)
  kCFUEdf  : 0.20  : CFU-E differentiation rate (day^-1)
  kRETmat  : 0.14  : Reticulocyte maturation (day^-1, ~7 days)
  kRBCd    : 0.0083: RBC death rate (day^-1, t1/2~120d)
  kMKPdf   : 0.15  : MKP differentiation rate (day^-1)
  kPLTd    : 0.10  : Platelet death rate (day^-1, t1/2~10d)
  kNEUPdf  : 0.33  : Neutrophil precursor maturation (day^-1)
  kNEUd    : 2.0   : Neutrophil death rate (day^-1, t1/2~12h)
  kTELOs   : 0.0005: Telomere shortening per HSC division
  // --- PK parameters: ATG ---
  kATGel   : 0.116 : ATG elimination rate (day^-1; hATG t1/2~6d)
  kATG12   : 0.5   : ATG distribution rate (day^-1)
  kATG21   : 0.2   : ATG redistribution rate (day^-1)
  kATGkill : 0.008 : ATG T-cell killing rate constant (mL/ug/day)
  // --- PK parameters: CsA ---
  kCSAabs  : 1.5   : CsA absorption rate (day^-1)
  kCSAel   : 0.35  : CsA elimination rate (day^-1)
  kCSA12   : 0.8   : CsA distribution (day^-1)
  kCSA21   : 0.4   : CsA redistribution (day^-1)
  IC50CSA  : 150   : CsA IC50 for IL-2/IFN-g suppression (ng/mL)
  VCSA     : 5.0   : CsA apparent Vd relative (L/kg scaling factor)
  // --- PK parameters: Eltrombopag ---
  kEPGabs  : 0.55  : EPAG absorption rate (day^-1)
  kEPGel   : 0.77  : EPAG elimination rate (day^-1, t1/2~21h)
  EmaxHSC  : 0.8   : EPAG max effect on HSC self-renewal
  EC50HSC  : 2.0   : EPAG EC50 for HSC expansion (ug/mL)
  EmaxPLT  : 0.9   : EPAG max effect on MKP production
  EC50PLT  : 1.5   : EPAG EC50 for platelet production (ug/mL)
  EmaxTR   : 0.3   : EPAG max Treg induction effect
  // --- Dosing flags (0 = off, 1 = on) ---
  hATG_on  : 0     : 1 = use hATG PK (ATGAM)
  rATG_on  : 0     : 1 = use rATG PK (Thymoglobulin; slower CL)
  CSA_on   : 0     : 1 = give CsA
  EPAG_on  : 0     : 1 = give Eltrombopag
  GCSF_on  : 0     : 1 = give G-CSF (neutrophil rescue factor)

$CMT @annotated
  CD8  : Autoreactive CD8 T cells (cells/uL)
  TREG : Regulatory T cells (cells/uL)
  CYTO : Cytokine index (normalized)
  HSC  : HSC pool (AU)
  PROG : Myeloid progenitor pool (AU)
  CFU_E: Erythroid progenitors (AU)
  RETIC: Reticulocytes (x10^9/L)
  RBC  : Hemoglobin (g/dL)
  MKP  : Megakaryocyte progenitors (AU)
  PLT  : Platelets (x10^9/L)
  NEUP : Neutrophil precursors (AU)
  NEU  : Neutrophils (x10^9/L)
  TELO : Telomere length (relative)
  CATG : ATG central (ug/mL)
  PATG : ATG peripheral (ug/mL)
  CGUT : CsA gut (mg/L)
  CCSA : CsA central (ng/mL)
  EGUT : EPAG gut (mg)
  CEPG : EPAG central (ug/mL)

$MAIN
  // Disease initial conditions at sAA diagnosis
  // (severe depletion of marrow; brisk T-cell activity)
  CD8_0  = 2000;    // elevated from normal ~800 (cells/uL)
  TREG_0 = 20;      // depleted from normal ~100 (cells/uL)
  HSC_0  = 5;       // 5% of normal (consistent with severe aplasia)
  double PROG_0 = 5;
  double CFUE_0 = 5;
  double RETIC_0 = 5;    // ~5 x10^9/L (ARC <20 in sAA)
  double RBC_0  = 7.0;   // Hgb ~7 g/dL
  double MKP_0  = 5;
  double PLT_0  = 10;    // PLT ~10 x10^9/L (threshold sAA <20k)
  double NEUP_0 = 3;
  double NEU_0  = 0.35;  // ANC ~0.35 x10^9/L (<0.5 = sAA)
  double TELO_0 = 0.90;  // slightly short telomeres
  double CYTO_0b= 3.0;   // elevated cytokine index at diagnosis

  CD8_0  = 2000;
  TREG_0 = 20;

$INIT
  CD8  = 2000,
  TREG = 20,
  CYTO = 3.0,
  HSC  = 5,
  PROG = 5,
  CFU_E= 5,
  RETIC= 5,
  RBC  = 7.0,
  MKP  = 5,
  PLT  = 10,
  NEUP = 3,
  NEU  = 0.35,
  TELO = 0.90,
  CATG = 0,
  PATG = 0,
  CGUT = 0,
  CCSA = 0,
  EGUT = 0,
  CEPG = 0

$ODE
  // ---- Drug effects ----------------------------------------
  // ATG: rATG has 50% slower elimination (longer exposure)
  double kATGel_eff = hATG_on*kATGel + rATG_on*(kATGel*0.23);
  double ATG_Tkill  = kATGkill * CATG;  // T-cell depletion by ATG (1/day per ug/mL)

  // CsA: inhibition of calcineurin → NFATc1 → IL-2/IFN-γ
  double CsA_inh   = CSA_on * CCSA / (CCSA + IC50CSA);  // 0→1 (0=no drug, 1=full inhibition)

  // EPAG: stimulate HSC self-renewal and MKP via MPL
  double EPAG_HSCeff = EPAG_on * EmaxHSC * CEPG / (CEPG + EC50HSC);
  double EPAG_PLTeff = EPAG_on * EmaxPLT * CEPG / (CEPG + EC50PLT);
  double EPAG_TReff  = EPAG_on * EmaxTR  * CEPG / (CEPG + EC50PLT);

  // G-CSF: boost neutrophil progenitor production
  double GCSF_boost = GCSF_on * 3.0;

  // ---- Cytokine & immune dynamics --------------------------
  // Treg-mediated suppression of CD8 (Hill function)
  double TReg_sup = 1.0 / (1.0 + (TREG / 50.0));  // 50 cells/uL = half-suppression

  // Cytokine production proportional to CD8 activity; suppressed by CsA
  double CYTO_prod = kCYp * (CD8/CD8_0) * TReg_sup * (1 - CsA_inh*0.7);

  // ATG depletes both CD8 and Treg
  dxdt_CD8  = kCD8p*CD8*(1.0 - CD8/(CD8_0*5)) - kCD8d*CD8
              - ATG_Tkill*CD8 - CsA_inh*0.4*kCD8p*CD8;
  dxdt_TREG = kTRp - kTRd*TREG - ATG_Tkill*0.7*TREG + EPAG_TReff*kTRp;
  dxdt_CYTO = CYTO_prod - kCYd*CYTO;

  // ---- HSC dynamics ----------------------------------------
  // IFN-γ/cytokine kills HSC; EPAG promotes self-renewal
  double HSC_kill = kCYkill * CYTO * HSC;
  double HSC_cap  = 100.0;  // normal = 100 AU
  // Telomere shortening reduces self-renewal capacity
  double TELO_eff = TELO;   // relative effect (1 = normal)

  dxdt_HSC  = kHSCsr * HSC * (1.0 - HSC/HSC_cap) * TELO_eff * (1 + EPAG_HSCeff)
             - kHSCd * HSC - kHSCdf * HSC - HSC_kill;

  // ---- Progenitor hierarchy --------------------------------
  // HSC differentiation feeds progenitor pool
  double PROG_cap = 100.0;
  dxdt_PROG = kHSCdf * HSC * 10 - kPROGdf * PROG;

  // Erythroid arm
  dxdt_CFU_E = kPROGdf * PROG * 0.35 - kCFUEdf * CFU_E;
  dxdt_RETIC = kCFUEdf * CFU_E * 2   - kRETmat * RETIC;
  dxdt_RBC   = kRETmat * RETIC * 0.07 - kRBCd  * RBC;
  // (scaling: 0.07 maps AU reticulocytes to g/dL Hgb units)

  // Megakaryocyte arm
  dxdt_MKP  = kPROGdf * PROG * 0.20 * (1 + EPAG_PLTeff) - kMKPdf * MKP;
  dxdt_PLT  = kMKPdf  * MKP * 1.2   - kPLTd * PLT;

  // Neutrophil arm
  dxdt_NEUP = kPROGdf * PROG * 0.45 * (1 + GCSF_boost) - kNEUPdf * NEUP;
  dxdt_NEU  = kNEUPdf * NEUP * 0.04  - kNEUd  * NEU;
  // (scaling: 0.04 maps AU precursors to x10^9/L neutrophils)

  // ---- Telomere dynamics -----------------------------------
  // Slow shortening with HSC cycling; stabilized if TELO_eff is high
  dxdt_TELO = -kTELOs * (kHSCsr * HSC / HSC_cap);  // proportional to HSC cycling

  // ---- ATG PK (2-compartment IV) --------------------------
  dxdt_CATG = -(kATGel_eff + kATG12) * CATG + kATG21 * PATG;
  dxdt_PATG = kATG12 * CATG - kATG21 * PATG;

  // ---- CsA PK (2-compartment oral) ------------------------
  dxdt_CGUT = -kCSAabs * CGUT;
  dxdt_CCSA = kCSAabs * CGUT / VCSA - (kCSAel + kCSA12)*CCSA + kCSA21*(CCSA*kCSA12/kCSA21);
  // Note: simplified; in practice CCSA unit-matched to ng/mL

  // ---- Eltrombopag PK (1-compartment oral) ----------------
  dxdt_EGUT = -kEPGabs * EGUT;
  dxdt_CEPG = kEPGabs * EGUT - kEPGel * CEPG;

$TABLE
  // Clinical outputs
  double ANC   = NEU;      // x10^9/L (= cells/uL / 1000)
  double HGB   = RBC;      // g/dL
  double PLTc  = PLT;      // x10^9/L
  double BM_cell = HSC;    // % of normal (AU / 100 * 100%)
  // Response classification
  double CR  = (ANC > 1.0 && PLTc > 100 && HGB > 10) ? 1 : 0;
  double PR  = (ANC > 0.5 && PLTc > 20  && HGB > 8 && !CR) ? 1 : 0;
  double NR  = (1 - CR - PR);
  // Derived
  double Hct = HGB / 3.0;  // approximate Hct (g/dL × 3 = %)
  double ATG_conc = CATG;
  double CsA_trough = CCSA;
  double EPAG_conc  = CEPG;
  double CD8_count  = CD8;
  double Treg_count = TREG;
  double Cytokine   = CYTO;

  capture ANC CsA_trough EPAG_conc ATG_conc
  capture HGB PLTc BM_cell Hct CR PR NR
  capture CD8_count Treg_count Cytokine TELO

$CAPTURE ANC HGB PLTc BM_cell CR PR NR CD8_count Treg_count Cytokine TELO
         ATG_conc CsA_trough EPAG_conc
'

## ---- Compile ------------------------------------------------
mod <- mcode("AplasticAnemia_QSP", code)

## ============================================================
## TREATMENT SCENARIO DEFINITIONS
## ============================================================

## Helper: build event table for a scenario
make_events <- function(scenario = 1) {

  ev_list <- list()

  if (scenario == 1) {
    ## Scenario 1: Untreated Severe AA (natural history)
    # No drug events; observe spontaneous progression over 1 year
    return(NULL)
  }

  if (scenario == 2) {
    ## Scenario 2: hATG + CsA (standard IST, CASG protocol)
    ## hATG 40 mg/kg/d × 4 days IV; CsA 5 mg/kg/d continuous oral
    ## Calibration target: 6-mo CR ~68% (Young 2012 NEJM)
    # hATG: 4 daily doses of ~2800 mg (70kg patient), as central-compartment bolus
    ev_hATG <- ev(ID=1, time=0:3, cmt="CATG", amt=2800, rate=-2, addl=0)
    # CsA: 350 mg/d oral continuously for 180 days (2 equal daily doses → simplify as QD)
    ev_CSA  <- ev(ID=1, time=0, cmt="CGUT", amt=350, rate=0, addl=179, ii=1)
    return(list(hATG_on=1, rATG_on=0, CSA_on=1, EPAG_on=0, GCSF_on=0,
                events=ev_hATG + ev_CSA))
  }

  if (scenario == 3) {
    ## Scenario 3: rATG + CsA (CASG protocol with Thymoglobulin)
    ## rATG 3.5 mg/kg/d × 5d IV (70 kg → ~245 mg/d)
    ## Calibration target: 6-mo CR ~37% (inferior to hATG; Young 2012 NEJM)
    ev_rATG <- ev(ID=1, time=0:4, cmt="CATG", amt=245, rate=-2)
    ev_CSA  <- ev(ID=1, time=0, cmt="CGUT", amt=350, addl=179, ii=1)
    return(list(hATG_on=0, rATG_on=1, CSA_on=1, EPAG_on=0, GCSF_on=0,
                events=ev_rATG + ev_CSA))
  }

  if (scenario == 4) {
    ## Scenario 4: hATG + CsA + Eltrombopag (Townsley 2017 NEJM protocol)
    ## EPAG 150 mg/d × 6 months then 75 mg/d × 3 months
    ## Calibration target: 6-mo CR ~58%
    ev_hATG <- ev(ID=1, time=0:3, cmt="CATG", amt=2800, rate=-2)
    ev_CSA  <- ev(ID=1, time=0, cmt="CGUT", amt=350, addl=179, ii=1)
    ev_EPAG_high <- ev(ID=1, time=0, cmt="EGUT", amt=150, addl=179, ii=1)
    ev_EPAG_low  <- ev(ID=1, time=180, cmt="EGUT", amt=75, addl=89, ii=1)
    return(list(hATG_on=1, rATG_on=0, CSA_on=1, EPAG_on=1, GCSF_on=0,
                events=ev_hATG + ev_CSA + ev_EPAG_high + ev_EPAG_low))
  }

  if (scenario == 5) {
    ## Scenario 5: Eltrombopag monotherapy (refractory/elderly)
    ## EPAG 75-150 mg/d; Desmond 2013 J Clin Invest: ~44% response
    ev_EPAG <- ev(ID=1, time=0, cmt="EGUT", amt=150, addl=179, ii=1)
    return(list(hATG_on=0, rATG_on=0, CSA_on=0, EPAG_on=1, GCSF_on=0,
                events=ev_EPAG))
  }

  if (scenario == 6) {
    ## Scenario 6: CsA monotherapy (non-severe AA, elderly/frail)
    ## Response rate ~30-40% in nsAA; not recommended for sAA
    ev_CSA <- ev(ID=1, time=0, cmt="CGUT", amt=350, addl=179, ii=1)
    return(list(hATG_on=0, rATG_on=0, CSA_on=1, EPAG_on=0, GCSF_on=0,
                events=ev_CSA))
  }
}

## ============================================================
## RUN ALL SCENARIOS
## ============================================================
sim_duration <- 365  # 12 months
delta_t      <- 1    # daily outputs

run_scenario <- function(sc_id, label) {
  sc <- make_events(sc_id)
  if (is.null(sc)) {
    # Untreated
    mod2 <- param(mod, hATG_on=0, rATG_on=0, CSA_on=0, EPAG_on=0, GCSF_on=0)
    out  <- mrgsim(mod2, end=sim_duration, delta=delta_t)
  } else {
    mod2 <- param(mod,
      hATG_on = sc$hATG_on, rATG_on = sc$rATG_on,
      CSA_on  = sc$CSA_on,  EPAG_on  = sc$EPAG_on,
      GCSF_on = sc$GCSF_on)
    out <- mrgsim(mod2, events=sc$events, end=sim_duration, delta=delta_t)
  }
  as.data.frame(out) %>% mutate(Scenario = label)
}

results <- bind_rows(
  run_scenario(1, "1: Untreated sAA"),
  run_scenario(2, "2: hATG + CsA"),
  run_scenario(3, "3: rATG + CsA"),
  run_scenario(4, "4: hATG + CsA + EPAG"),
  run_scenario(5, "5: EPAG Monotherapy"),
  run_scenario(6, "6: CsA Monotherapy")
)

## ============================================================
## RESPONSE RATES AT 6 MONTHS (calibration check)
## ============================================================
resp_6mo <- results %>%
  filter(time == 180) %>%
  group_by(Scenario) %>%
  summarise(
    ANC_x10_9  = round(mean(ANC), 3),
    Hgb_gdL    = round(mean(HGB), 2),
    PLT_x10_9  = round(mean(PLTc), 1),
    CR_rate    = round(mean(CR) * 100, 1),
    PR_rate    = round(mean(PR) * 100, 1),
    NR_rate    = round(mean(NR) * 100, 1)
  )

cat("\n=== 6-Month Response Summary ===\n")
print(resp_6mo)

cat("\n=== Clinical Calibration Targets ===\n")
cat("Scenario 2 (hATG+CsA): target CR ~68% | Young 2012 NEJM\n")
cat("Scenario 3 (rATG+CsA): target CR ~37% | Young 2012 NEJM\n")
cat("Scenario 4 (hATG+CsA+EPAG): target CR ~58% | Townsley 2017 NEJM\n")
cat("Scenario 5 (EPAG mono): target resp ~44% | Desmond 2013 J Clin Invest\n")

## ============================================================
## VISUALIZATION
## ============================================================
theme_qsp <- theme_bw() +
  theme(strip.background = element_rect(fill="#E0E0F0"),
        legend.position  = "bottom",
        legend.key.width = unit(1.5, "cm"))

colors6 <- c("#999999","#2166AC","#4DAF4A","#FF7F00","#E41A1C","#984EA3")

## Figure 1: Blood count recovery
p1 <- results %>%
  select(time, Scenario, ANC, HGB, PLTc) %>%
  pivot_longer(cols=c(ANC, HGB, PLTc), names_to="Measure", values_to="Value") %>%
  mutate(Measure = recode(Measure,
    ANC="ANC (×10⁹/L)", HGB="Hemoglobin (g/dL)", PLTc="Platelets (×10⁹/L)")) %>%
  ggplot(aes(x=time, y=Value, color=Scenario)) +
  geom_line(size=0.9) +
  facet_wrap(~Measure, scales="free_y", ncol=1) +
  scale_color_manual(values=colors6) +
  geom_vline(xintercept=180, linetype=2, color="grey40") +
  labs(title="Aplastic Anemia: Blood Count Recovery by Treatment",
       subtitle="Dashed line = 6-month response assessment",
       x="Time (days)", y="Value", color="Treatment") +
  theme_qsp
print(p1)

## Figure 2: HSC dynamics and cytokine
p2 <- results %>%
  select(time, Scenario, BM_cell, Cytokine, CD8_count) %>%
  pivot_longer(cols=c(BM_cell, Cytokine, CD8_count),
               names_to="Measure", values_to="Value") %>%
  mutate(Measure = recode(Measure,
    BM_cell="BM HSC Pool (% normal)", Cytokine="Cytokine Index",
    CD8_count="CD8+ T Cells (cells/μL)")) %>%
  ggplot(aes(x=time, y=Value, color=Scenario)) +
  geom_line(size=0.9) +
  facet_wrap(~Measure, scales="free_y", ncol=1) +
  scale_color_manual(values=colors6) +
  geom_vline(xintercept=180, linetype=2, color="grey40") +
  labs(title="Aplastic Anemia: Immune & BM Dynamics",
       x="Time (days)", y="Value", color="Treatment") +
  theme_qsp
print(p2)

## Figure 3: Drug PK profiles (Scenario 4: triple therapy)
pk_data <- results %>%
  filter(Scenario == "4: hATG + CsA + EPAG") %>%
  select(time, ATG_conc, CsA_trough, EPAG_conc) %>%
  pivot_longer(-time, names_to="Drug", values_to="Conc") %>%
  mutate(Drug = recode(Drug,
    ATG_conc="ATG (μg/mL)", CsA_trough="CsA Cmin (ng/mL)", EPAG_conc="EPAG (μg/mL)"))

p3 <- ggplot(pk_data, aes(x=time, y=Conc, color=Drug)) +
  geom_line(size=1.0) +
  facet_wrap(~Drug, scales="free_y", ncol=1) +
  labs(title="Drug PK Profiles – Scenario 4: hATG + CsA + Eltrombopag",
       x="Time (days)", y="Concentration") +
  scale_color_brewer(palette="Set1") +
  theme_qsp
print(p3)

## ============================================================
## SENSITIVITY: EPAG DOSE ESCALATION ON PLATELET RECOVERY
## ============================================================
epag_doses <- c(0, 50, 75, 100, 150)
epag_sens <- bind_rows(lapply(epag_doses, function(dose) {
  ev_EPAG <- ev(time=0, cmt="EGUT", amt=dose, addl=179, ii=1)
  mod_e <- param(mod, hATG_on=1, rATG_on=0, CSA_on=1, EPAG_on=ifelse(dose>0,1,0), GCSF_on=0)
  out   <- mrgsim(mod_e, events=ev_EPAG, end=270, delta=1)
  as.data.frame(out) %>% mutate(EPAG_dose_mg = dose)
}))

p4 <- ggplot(epag_sens, aes(x=time, y=PLTc, color=factor(EPAG_dose_mg))) +
  geom_line(size=0.9) +
  scale_color_viridis_d(name="EPAG dose (mg/d)") +
  geom_hline(yintercept=100, linetype=2, color="red", size=0.8) +
  annotate("text", x=260, y=105, label="CR threshold (PLT=100)", size=3, color="red") +
  labs(title="EPAG Dose-Response: Platelet Recovery (hATG+CsA backbone)",
       x="Time (days)", y="Platelets (×10⁹/L)") +
  theme_qsp
print(p4)

cat("\n=== Model run complete ===\n")
cat("All 6 treatment scenarios simulated over 365 days.\n")
cat("Output: results data frame with 13 captured variables.\n")
