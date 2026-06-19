## =============================================================================
## PSC QSP Shiny App — Primary Sclerosing Cholangitis Interactive Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Drug Pharmacokinetics (PK)
##   3. Bile Acid & FXR Signaling
##   4. Immune & Inflammatory Dynamics
##   5. Fibrosis & Clinical Biomarkers
##   6. Treatment Scenarios & Long-term Outcomes
##   7. Biomarker Dashboard & CCA Risk
##   8. Dose-Response & Sensitivity Analysis
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(patchwork)

## ---------------------------------------------------------------------------
## Inline model code (self-contained Shiny)
## ---------------------------------------------------------------------------
psc_model_code <- '
$PROB PSC QSP Shiny Model

$PARAM @annotated
Ka_UDCA   : 0.8   : UDCA absorption rate (1/h)
F_UDCA    : 0.50  : UDCA bioavailability
CL_UDCA   : 12.0  : UDCA clearance (L/h)
Vd_UDCA   : 25.0  : UDCA Vd (L)
Kbile_UDCA: 0.15  : UDCA biliary transfer (1/h)
Ka_OCA    : 0.6   : OCA absorption rate (1/h)
F_OCA     : 0.60  : OCA bioavailability
CL_OCA    : 8.0   : OCA clearance (L/h)
Vd_OCA    : 20.0  : OCA Vd (L)
Kbile_OCA : 0.12  : OCA biliary transfer (1/h)
IC50_OCA_FXR: 0.10 : OCA FXR IC50 (μmol/L)
Ka_BEZ    : 1.0   : BEZ absorption (1/h)
CL_BEZ    : 15.0  : BEZ clearance (L/h)
Vd_BEZ    : 18.0  : BEZ Vd (L)
k_LPS_prod   : 0.05 : LPS production rate
k_LPS_clear  : 0.10 : LPS clearance (1/h)
k_barrier    : 0.02 : Gut barrier repair (1/h)
k_IBD_damage : 0.005 : IBD barrier damage
GutBarrier0  : 1.0  : Baseline gut barrier
k_FXR_act   : 0.30  : FXR activation rate
k_FXR_decay : 0.05  : FXR decay (1/h)
FXR0        : 0.30  : Baseline FXR activity
k_BA_synth  : 0.10  : BA synthesis rate
k_BA_excrete: 0.08  : BA excretion rate
k_hydro_form : 0.04 : Hydrophobic BA formation
k_hydro_clear: 0.06 : Hydrophobic BA clearance
HydroIndex0 : 0.30  : Baseline hydrophobic index
UDCA_hydro_eff: 0.50 : UDCA hydrophobic BA reduction
k_IL17_base  : 0.002 : Basal IL-17A
k_IL17_LPS   : 0.05  : LPS → IL-17A
k_IL17_BA    : 0.04  : HydroBA → IL-17A
k_IL17_decay : 0.08  : IL-17A decay
k_TNFa_base  : 0.003 : Basal TNF-α
k_TNFa_LPS   : 0.06  : LPS → TNF-α
k_TNFa_decay : 0.10  : TNF-α decay
k_IL6_IL17   : 0.03  : IL-17A → IL-6
k_IL6_TNFa   : 0.04  : TNF-α → IL-6
k_IL6_decay  : 0.12  : IL-6 decay
k_Treg_base  : 0.001 : Basal Treg
k_Treg_decay : 0.06  : Treg decay
k_Treg_inhib : 0.03  : IL-6 inhibits Treg
Cholangio0   : 1.0   : Baseline cholangiocyte health
k_cholangio_damage_IL17: 0.02 : IL-17A damages cholangiocyte
k_cholangio_damage_BA  : 0.03 : HydroBA damages cholangiocyte
k_cholangio_repair     : 0.01 : Cholangiocyte repair
k_senesc_form: 0.015  : Senescence formation
k_senesc_clear: 0.008 : Senescence clearance
k_HSC_act    : 0.01  : HSC activation rate
k_HSC_resol  : 0.005 : HSC resolution
k_Col_synth  : 0.008 : Collagen synthesis
k_Col_degrad : 0.004 : Collagen degradation
LOXL2_0      : 0.20  : Baseline LOXL2
k_LOXL2_form : 0.005 : LOXL2 formation
k_LOXL2_clear: 0.008 : LOXL2 clearance
ALP0         : 200.0 : Baseline ALP (IU/L)
k_ALP_chol   : 5.0   : Cholestasis → ALP
k_ALP_repair : 0.30  : ALP normalization
ALP_normal   : 100.0 : Normal ALP
Bili0        : 1.5   : Baseline bilirubin
k_Bili_chol  : 0.20  : Cholestasis → bilirubin
k_Bili_synth : 0.05  : Normal bilirubin production
k_Bili_excrete: 0.10 : Bilirubin excretion
Fibro0       : 8.0   : Baseline liver stiffness kPa
k_Fibro_form : 0.10  : Fibrosis → stiffness
k_Fibro_resol: 0.005 : Stiffness resolution
PP0          : 5.0   : Baseline portal pressure mmHg
k_PP_fibrosis: 0.50  : Col → portal pressure
k_PP_decay   : 0.02  : Portal pressure resolution
k_CCA_chol   : 0.0001 : CCA risk accumulation
k_CCA_senesc : 0.0002 : Senescence CCA risk
IBD_status   : 1.0   : IBD co-existence flag

$CMT @annotated
UDCA_gut : UDCA gut (mg)
UDCA_plasma : UDCA plasma (mg)
UDCA_bile : UDCA bile (mg)
OCA_gut : OCA gut (mg)
OCA_plasma : OCA plasma (mg)
OCA_bile : OCA bile (mg)
BEZ_plasma : Bezafibrate plasma (mg)
LPS : Portal LPS (AU)
GutBarrier : Gut barrier (0-1)
FXR_act : FXR activation (0-1)
BilePool : Bile pool (AU)
HydroIndex : Hydrophobic BA (AU)
IL17A : IL-17A (AU)
TNFa : TNF-α (AU)
IL6 : IL-6 (AU)
Treg_IL10 : Treg/IL-10 (AU)
Cholangio_health : Cholangiocyte health (0-1)
Senescence : Senescent cells (AU)
HSC_act : Activated HSC (0-1)
Col1a1 : Collagen index (AU)
LOXL2 : LOXL2 (AU)
ALP : ALP (IU/L)
Bilirubin : Bilirubin (mg/dL)
Fibroscan : Liver stiffness (kPa)
PortalPressure : Portal pressure (mmHg)
CCA_risk : CCA risk (0-1)

$MAIN
double OCA_bile_conc  = OCA_bile / 0.5;
double UDCA_bile_conc = UDCA_bile / 2.0;
double FXR_OCA   = OCA_bile_conc  / (IC50_OCA_FXR + OCA_bile_conc);
double FXR_UDCA  = 0.05 * UDCA_bile_conc / (5.0 + UDCA_bile_conc);
double FXR_input = FXR0 + (1.0 - FXR0) * (FXR_OCA + FXR_UDCA);
FXR_input = (FXR_input > 1.0) ? 1.0 : FXR_input;
double Cholestasis = 1.0 - Cholangio_health * BilePool;
Cholestasis = (Cholestasis < 0.0) ? 0.0 : Cholestasis;
Cholestasis = (Cholestasis > 1.0) ? 1.0 : Cholestasis;

$ODE
dxdt_UDCA_gut    = -Ka_UDCA * UDCA_gut;
dxdt_UDCA_plasma = F_UDCA * Ka_UDCA * UDCA_gut - (CL_UDCA/Vd_UDCA)*UDCA_plasma;
dxdt_UDCA_bile   = Kbile_UDCA * UDCA_plasma - 0.05 * UDCA_bile;
dxdt_OCA_gut    = -Ka_OCA * OCA_gut;
dxdt_OCA_plasma = F_OCA * Ka_OCA * OCA_gut - (CL_OCA/Vd_OCA)*OCA_plasma;
dxdt_OCA_bile   = Kbile_OCA * OCA_plasma - 0.04 * OCA_bile;
dxdt_BEZ_plasma = Ka_BEZ * 1.0 - (CL_BEZ/Vd_BEZ)*BEZ_plasma;
double LPS_prod = k_LPS_prod * (2.0 - GutBarrier) * (1.0 + IBD_status * 0.5);
dxdt_LPS = LPS_prod - k_LPS_clear * LPS;
double barrier_damage = k_IBD_damage * IBD_status + 0.05 * HydroIndex * LPS + 0.02 * IL17A;
dxdt_GutBarrier = k_barrier * (GutBarrier0 - GutBarrier) - barrier_damage * GutBarrier;
dxdt_FXR_act = k_FXR_act * (FXR_input - FXR_act) - k_FXR_decay * (FXR_act - FXR0);
double BA_synth = k_BA_synth * (1.0 - 0.7 * FXR_act);
dxdt_BilePool = BA_synth - k_BA_excrete * BilePool;
double hydro_decrease = k_hydro_clear * HydroIndex * (1.0 + UDCA_hydro_eff * UDCA_bile_conc / (2.0 + UDCA_bile_conc));
dxdt_HydroIndex = k_hydro_form * (1.0 - FXR_act) * LPS - hydro_decrease;
double IL17_prod = k_IL17_base + k_IL17_LPS * LPS * (1.0 - Treg_IL10 * 0.5) + k_IL17_BA * HydroIndex;
dxdt_IL17A = IL17_prod - k_IL17_decay * IL17A;
dxdt_TNFa = k_TNFa_base + k_TNFa_LPS * LPS - k_TNFa_decay * TNFa;
dxdt_IL6 = k_IL6_IL17 * IL17A + k_IL6_TNFa * TNFa - k_IL6_decay * IL6;
dxdt_Treg_IL10 = k_Treg_base - k_Treg_decay * Treg_IL10 - k_Treg_inhib * IL6 * Treg_IL10;
double cholangio_damage = k_cholangio_damage_IL17 * IL17A + k_cholangio_damage_BA * HydroIndex + 0.01 * TNFa;
dxdt_Cholangio_health = k_cholangio_repair * (1.0 - Cholangio_health) - cholangio_damage * Cholangio_health;
dxdt_Senescence = k_senesc_form * IL17A * HydroIndex - k_senesc_clear * Senescence;
double HSC_activation = k_HSC_act * (IL17A + Senescence) * (1.0 - HSC_act);
dxdt_HSC_act = HSC_activation - k_HSC_resol * HSC_act * Treg_IL10;
double BEZ_antifib = 1.0 - 0.3 * BEZ_plasma / (5.0 + BEZ_plasma);
dxdt_Col1a1 = k_Col_synth * HSC_act * BEZ_antifib - k_Col_degrad * Col1a1 * (1.0 - LOXL2 * 0.3);
dxdt_LOXL2 = k_LOXL2_form * HSC_act - k_LOXL2_clear * LOXL2;
double OCA_ALP_eff = 0.25 * OCA_bile_conc / (2.0 + OCA_bile_conc);
double UDCA_ALP_eff = 0.15 * UDCA_bile_conc / (5.0 + UDCA_bile_conc);
double BEZ_ALP_eff = 0.20 * BEZ_plasma / (10.0 + BEZ_plasma);
dxdt_ALP = k_ALP_chol * Cholestasis * ALP0 - k_ALP_repair * (ALP - ALP_normal) - (OCA_ALP_eff + UDCA_ALP_eff + BEZ_ALP_eff) * ALP;
dxdt_Bilirubin = k_Bili_synth + k_Bili_chol * Cholestasis - k_Bili_excrete * Bilirubin;
dxdt_Fibroscan = k_Fibro_form * Col1a1 - k_Fibro_resol * Fibroscan;
dxdt_PortalPressure = k_PP_fibrosis * Col1a1 - k_PP_decay * (PortalPressure - PP0);
dxdt_CCA_risk = k_CCA_chol * Cholestasis * ALP / ALP0 + k_CCA_senesc * Senescence;

$INIT
UDCA_gut = 0
UDCA_plasma = 0
UDCA_bile = 0
OCA_gut = 0
OCA_plasma = 0
OCA_bile = 0
BEZ_plasma = 0
LPS = 0.6
GutBarrier = 0.55
FXR_act = 0.25
BilePool = 0.70
HydroIndex = 0.55
IL17A = 0.25
TNFa = 0.20
IL6 = 0.15
Treg_IL10 = 0.05
Cholangio_health = 0.50
Senescence = 0.30
HSC_act = 0.35
Col1a1 = 0.40
LOXL2 = 0.35
ALP = 380.0
Bilirubin = 2.8
Fibroscan = 12.5
PortalPressure = 8.0
CCA_risk = 0.02

$CAPTURE
Cholestasis
FXR_input
OCA_bile_conc
UDCA_bile_conc
'

## Compile model once at app load
psc_mod <- mcode("PSC_Shiny", psc_model_code)

## ---------------------------------------------------------------------------
## Simulation helper
## ---------------------------------------------------------------------------
run_sim <- function(
    udca_dose = 0, oca_dose = 0, bez_dose = 0,
    ibd_status = 1,
    sim_years = 3,
    freq_h = 12
) {
  params_list <- list(IBD_status = ibd_status)
  mod <- param(psc_mod, params_list)

  sim_hours <- sim_years * 365 * 24
  obs_times <- seq(0, sim_hours, by = 24 * 7)

  events_list <- list()

  if (udca_dose > 0) {
    times_u <- seq(0, sim_hours - freq_h, by = freq_h)
    events_list[["UDCA"]] <- ev(
      time = times_u,
      amt  = udca_dose / (24 / freq_h),
      cmt  = "UDCA_gut",
      evid = 1
    )
  }
  if (oca_dose > 0) {
    times_o <- seq(0, sim_hours - 24, by = 24)
    events_list[["OCA"]] <- ev(
      time = times_o,
      amt  = oca_dose,
      cmt  = "OCA_gut",
      evid = 1
    )
  }
  if (bez_dose > 0) {
    times_b <- seq(0, sim_hours - 24, by = 24)
    events_list[["BEZ"]] <- ev(
      time = times_b,
      amt  = bez_dose,
      cmt  = "BEZ_plasma",
      evid = 1
    )
  }

  if (length(events_list) == 0) {
    out <- mrgsim(mod, tgrid = obs_times)
  } else {
    all_ev <- do.call(c, events_list)
    out <- mrgsim(mod, events = all_ev, tgrid = obs_times)
  }

  as.data.frame(out) %>%
    mutate(time_years = time / (365 * 24))
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "PSC QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",               tabName = "tab_pk",       icon = icon("capsules")),
      menuItem("Bile Acid / FXR",       tabName = "tab_ba",       icon = icon("vial")),
      menuItem("Immune Dynamics",       tabName = "tab_immune",   icon = icon("shield-virus")),
      menuItem("Fibrosis & Biomarkers", tabName = "tab_fibrosis", icon = icon("chart-line")),
      menuItem("Treatment Outcomes",    tabName = "tab_outcomes",  icon = icon("heartbeat")),
      menuItem("CCA Risk & Biomarkers", tabName = "tab_cca",      icon = icon("exclamation-triangle")),
      menuItem("Dose-Response",         tabName = "tab_dr",       icon = icon("sliders-h"))
    ),

    hr(),
    h5("Treatment Parameters", style = "color:white; padding-left:15px"),
    sliderInput("udca_dose",  "UDCA (mg/day):",    0, 1500, 900, step = 150),
    sliderInput("oca_dose",   "OCA (mg/day):",     0, 50,   10,  step = 5),
    sliderInput("bez_dose",   "Bezafibrate (mg):", 0, 400,  0,   step = 100),
    hr(),
    radioButtons("ibd_status", "IBD Status:",
                 choices = c("PSC+IBD" = 1, "PSC only" = 0),
                 selected = 1, inline = TRUE),
    sliderInput("sim_years", "Simulation (years):", 1, 10, 5, step = 1),
    actionButton("run_btn", "Run Simulation", icon = icon("play"),
                 class = "btn-success", width = "90%"),
    br()
  ),

  dashboardBody(
    tabItems(

      ## -----------------------------------------------------------------------
      ## TAB 1: Patient Profile
      ## -----------------------------------------------------------------------
      tabItem("tab_patient",
        fluidRow(
          box(title = "Disease Overview — Primary Sclerosing Cholangitis", width = 12,
              status = "primary", solidHeader = TRUE,
              HTML("
              <h4>What is PSC?</h4>
              <p>Primary Sclerosing Cholangitis (PSC) is a chronic, progressive cholestatic liver disease
              characterized by inflammation, fibrosis, and stricturing of both intrahepatic and extrahepatic bile ducts.
              Median transplant-free survival is ~12–21 years from diagnosis.</p>
              <hr>
              <div style='display:flex; gap:20px; flex-wrap:wrap'>
              <div style='flex:1; min-width:200px'>
              <h5>Epidemiology</h5>
              <ul>
                <li>Prevalence: ~10–16/100,000 (Northern Europe/US)</li>
                <li>Male predominance (~65%)</li>
                <li>Peak onset: 30–40 years</li>
                <li>~70% have concurrent IBD (UC > CD)</li>
                <li>Lifetime CCA risk: 10–20%</li>
              </ul>
              </div>
              <div style='flex:1; min-width:200px'>
              <h5>Key Pathomechanisms</h5>
              <ul>
                <li>Gut dysbiosis → LPS → TLR4 → Kupffer cell activation</li>
                <li>Th17/CD8+ T-cell mediated cholangiocyte damage</li>
                <li>Cholangiocyte senescence (SASP) → peribiliary fibrosis</li>
                <li>FXR/TGR5 signaling dysregulation</li>
                <li>Portal fibroblast activation → collagen cross-linking (LOXL2)</li>
              </ul>
              </div>
              <div style='flex:1; min-width:200px'>
              <h5>Drug Targets Modeled</h5>
              <ul>
                <li><b>UDCA</b>: ↓ bile hydrophobicity, ↓ ALP</li>
                <li><b>OCA (obeticholic acid)</b>: FXR agonist → ↓ BA synthesis</li>
                <li><b>Bezafibrate</b>: PPARα → ↓ TGF-β, ↓ ALP</li>
                <li><b>Simtuzumab</b>: Anti-LOXL2 (anti-fibrotic)</li>
                <li><b>Vedolizumab</b>: Anti-α4β7 → ↓ gut-homing T cells</li>
              </ul>
              </div>
              </div>
              ")
          )
        ),
        fluidRow(
          valueBoxOutput("box_ALP",     width = 3),
          valueBoxOutput("box_Bili",    width = 3),
          valueBoxOutput("box_Fibro",   width = 3),
          valueBoxOutput("box_PP",      width = 3)
        ),
        fluidRow(
          box(title = "Baseline Disease State (Simulated)", width = 12,
              status = "info", solidHeader = TRUE,
              DTOutput("tbl_baseline"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 2: Drug PK
      ## -----------------------------------------------------------------------
      tabItem("tab_pk",
        fluidRow(
          box(title = "UDCA Plasma & Biliary Concentrations", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_udca_pk", height = "350px")),
          box(title = "OCA Plasma & Biliary Concentrations", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_oca_pk", height = "350px"))
        ),
        fluidRow(
          box(title = "Bezafibrate Plasma Concentration", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plt_bez_pk", height = "300px")),
          box(title = "PK Parameter Summary", width = 6,
              status = "info", solidHeader = TRUE,
              DTOutput("tbl_pk"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 3: Bile Acid & FXR
      ## -----------------------------------------------------------------------
      tabItem("tab_ba",
        fluidRow(
          box(title = "FXR Activation (0-1 scale)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_fxr", height = "300px")),
          box(title = "Bile Pool & Hydrophobic BA Index", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_ba_pool", height = "300px"))
        ),
        fluidRow(
          box(title = "Gut-Liver Axis — LPS & Gut Barrier", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plt_gut", height = "300px")),
          box(title = "Bile Acid Pathway Explanation", width = 6,
              status = "primary", solidHeader = TRUE,
              HTML("<h5>FXR Signaling in PSC</h5>
              <p>The farnesoid X receptor (FXR) is the master regulator of bile acid homeostasis.
              In PSC, FXR activity is impaired, leading to:</p>
              <ul>
              <li><b>↑ CYP7A1</b> → excess BA synthesis</li>
              <li><b>↓ BSEP / MRP2</b> → bile acid retention (cholestasis)</li>
              <li><b>↑ Hydrophobic BA</b> → cholangiocyte toxicity</li>
              </ul>
              <h5>Drug Interventions:</h5>
              <ul>
              <li><b>OCA</b>: Potent FXR agonist (IC50 ~0.1 μmol/L), restores FXR signaling</li>
              <li><b>UDCA</b>: Dilutes hydrophobic BA pool, partial FXR agonism</li>
              <li><b>norUDCA</b>: Cholehepatic shunting → HCO₃⁻ umbrella protection</li>
              <li><b>Bezafibrate</b>: PPARα-mediated FGF21 → indirect CYP7A1 suppression</li>
              </ul>"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 4: Immune Dynamics
      ## -----------------------------------------------------------------------
      tabItem("tab_immune",
        fluidRow(
          box(title = "IL-17A (Th17) Kinetics", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_il17", height = "300px")),
          box(title = "TNF-α & IL-6", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_tnf_il6", height = "300px"))
        ),
        fluidRow(
          box(title = "Regulatory T Cells (Treg/IL-10)", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plt_treg", height = "300px")),
          box(title = "Cholangiocyte Health & Senescence", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plt_cholangio", height = "300px"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 5: Fibrosis & Biomarkers
      ## -----------------------------------------------------------------------
      tabItem("tab_fibrosis",
        fluidRow(
          box(title = "Hepatic Collagen (Fibrosis Index)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_col", height = "300px")),
          box(title = "LOXL2 (Collagen Cross-linking)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_loxl2", height = "300px"))
        ),
        fluidRow(
          box(title = "Liver Stiffness (FibroScan kPa)", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plt_fibroscan", height = "300px")),
          box(title = "Portal Pressure (mmHg)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_pp", height = "300px"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 6: Treatment Outcomes
      ## -----------------------------------------------------------------------
      tabItem("tab_outcomes",
        fluidRow(
          box(title = "ALP (Alkaline Phosphatase, IU/L)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_alp", height = "300px")),
          box(title = "Total Bilirubin (mg/dL)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_bili", height = "300px"))
        ),
        fluidRow(
          box(title = "ALP Response Rate Over Time", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plt_alp_resp", height = "300px")),
          box(title = "Outcome Table at Years 1, 3, 5", width = 6,
              status = "info", solidHeader = TRUE,
              DTOutput("tbl_outcomes"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 7: CCA Risk & Biomarker Dashboard
      ## -----------------------------------------------------------------------
      tabItem("tab_cca",
        fluidRow(
          box(title = "Cumulative CCA Risk Index (0-1)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_cca", height = "300px")),
          box(title = "Senescent Cell Burden", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_senesc", height = "300px"))
        ),
        fluidRow(
          box(title = "HSC Activation (Myofibroblast)", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plt_hsc", height = "300px")),
          box(title = "Multi-Biomarker Radar at Year 2", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_radar", height = "300px"))
        )
      ),

      ## -----------------------------------------------------------------------
      ## TAB 8: Dose-Response
      ## -----------------------------------------------------------------------
      tabItem("tab_dr",
        fluidRow(
          box(title = "OCA Dose-Response: ALP at Year 1", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_oca_dr", height = "350px")),
          box(title = "UDCA Dose-Response: ALP at Year 1", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plt_udca_dr", height = "350px"))
        ),
        fluidRow(
          box(title = "IBD Impact on PSC Progression (ALP)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plt_ibd", height = "300px")),
          box(title = "Sensitivity Analysis: Parameters affecting ALP", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plt_sensitivity", height = "300px"))
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Running PSC simulation...", {
      run_sim(
        udca_dose  = input$udca_dose,
        oca_dose   = input$oca_dose,
        bez_dose   = input$bez_dose,
        ibd_status = as.numeric(input$ibd_status),
        sim_years  = input$sim_years
      )
    })
  }, ignoreNULL = FALSE)

  # --- Value Boxes (Tab 1) ---
  output$box_ALP <- renderValueBox({
    d <- sim_data()
    alp_now <- round(tail(d$ALP, 1), 0)
    color <- if (alp_now < 150) "green" else if (alp_now < 300) "yellow" else "red"
    valueBox(paste0(alp_now, " IU/L"), "Current ALP", icon = icon("flask"), color = color)
  })
  output$box_Bili <- renderValueBox({
    d <- sim_data()
    bili_now <- round(tail(d$Bilirubin, 1), 2)
    color <- if (bili_now < 1.2) "green" else if (bili_now < 3) "yellow" else "red"
    valueBox(paste0(bili_now, " mg/dL"), "Bilirubin", icon = icon("tint"), color = color)
  })
  output$box_Fibro <- renderValueBox({
    d <- sim_data()
    fibro_now <- round(tail(d$Fibroscan, 1), 1)
    color <- if (fibro_now < 9.5) "green" else if (fibro_now < 15) "yellow" else "red"
    valueBox(paste0(fibro_now, " kPa"), "Liver Stiffness", icon = icon("lungs"), color = color)
  })
  output$box_PP <- renderValueBox({
    d <- sim_data()
    pp_now <- round(tail(d$PortalPressure, 1), 1)
    color <- if (pp_now < 10) "green" else if (pp_now < 14) "yellow" else "red"
    valueBox(paste0(pp_now, " mmHg"), "Portal Pressure", icon = icon("heartbeat"), color = color)
  })

  # Baseline table
  output$tbl_baseline <- renderDT({
    d <- sim_data()
    first <- d[1, ]
    df <- data.frame(
      Variable = c("ALP (IU/L)", "Bilirubin (mg/dL)", "Liver Stiffness (kPa)",
                   "Portal Pressure (mmHg)", "IL-17A (AU)", "TNF-α (AU)",
                   "Cholangiocyte Health", "HSC Activation", "Collagen Index", "CCA Risk"),
      Baseline = c(first$ALP, first$Bilirubin, first$Fibroscan, first$PortalPressure,
                   round(first$IL17A, 3), round(first$TNFa, 3),
                   round(first$Cholangio_health, 3), round(first$HSC_act, 3),
                   round(first$Col1a1, 3), round(first$CCA_risk, 4))
    )
    datatable(df, options = list(pageLength = 10, dom = "t"), rownames = FALSE)
  })

  ## PK Plots
  make_pk_plot <- function(d, var1, var2, label1, label2, title_str) {
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~.data[[var1]], name = label1, line = list(color = "#3498DB", width = 2)) %>%
      add_lines(y = ~.data[[var2]], name = label2, line = list(color = "#E74C3C", width = 2, dash = "dash")) %>%
      layout(title = title_str,
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "Concentration (mg)"),
             legend = list(orientation = "h"))
  }

  output$plt_udca_pk <- renderPlotly({
    d <- sim_data()
    make_pk_plot(d, "UDCA_plasma", "UDCA_bile", "UDCA Plasma", "UDCA Bile", "UDCA PK")
  })
  output$plt_oca_pk <- renderPlotly({
    d <- sim_data()
    make_pk_plot(d, "OCA_plasma", "OCA_bile", "OCA Plasma", "OCA Bile", "OCA PK")
  })
  output$plt_bez_pk <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~BEZ_plasma, type = "scatter", mode = "lines",
            line = list(color = "#E67E22", width = 2)) %>%
      layout(title = "Bezafibrate Plasma", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "Plasma Conc (mg)"))
  })
  output$tbl_pk <- renderDT({
    df <- data.frame(
      Drug = c("UDCA", "OCA", "Bezafibrate"),
      `Dose (mg/day)` = c(input$udca_dose, input$oca_dose, input$bez_dose),
      `F (%)` = c(50, 60, 100),
      `CL (L/h)` = c(12.0, 8.0, 15.0),
      `Vd (L)` = c(25.0, 20.0, 18.0),
      `Target` = c("Bile hydrophobicity", "FXR agonist", "PPARα agonist")
    )
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  ## Bile Acid / FXR plots
  output$plt_fxr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~FXR_act, name = "FXR Activity", line = list(color = "#9B59B6", width = 2)) %>%
      add_lines(y = ~FXR_input, name = "FXR Target", line = list(color = "#F39C12", width = 2, dash = "dot")) %>%
      layout(title = "FXR Activation", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "FXR Activity (0-1)", range = c(0, 1)))
  })
  output$plt_ba_pool <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~BilePool, name = "Bile Pool", line = list(color = "#27AE60", width = 2)) %>%
      add_lines(y = ~HydroIndex, name = "Hydrophobic BA Index", line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "Bile Pool & Hydrophobic BA",
             xaxis = list(title = "Time (years)"), yaxis = list(title = "AU (normalized)"))
  })
  output$plt_gut <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~LPS, name = "Portal LPS", line = list(color = "#E74C3C", width = 2)) %>%
      add_lines(y = ~GutBarrier, name = "Gut Barrier", line = list(color = "#27AE60", width = 2)) %>%
      layout(title = "Gut-Liver Axis", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "AU"))
  })

  ## Immune plots
  output$plt_il17 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~IL17A, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "IL-17A", xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })
  output$plt_tnf_il6 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~TNFa, name = "TNF-α", line = list(color = "#C0392B", width = 2)) %>%
      add_lines(y = ~IL6, name = "IL-6", line = list(color = "#E67E22", width = 2)) %>%
      layout(title = "TNF-α & IL-6", xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })
  output$plt_treg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~Treg_IL10, type = "scatter", mode = "lines",
            line = list(color = "#27AE60", width = 2)) %>%
      layout(title = "Treg / IL-10", xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })
  output$plt_cholangio <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years) %>%
      add_lines(y = ~Cholangio_health, name = "Cholangiocyte Health", line = list(color = "#27AE60", width = 2)) %>%
      add_lines(y = ~Senescence, name = "Senescent Cells", line = list(color = "#E67E22", width = 2)) %>%
      layout(title = "Cholangiocyte Health & Senescence",
             xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })

  ## Fibrosis plots
  output$plt_col <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~Col1a1, type = "scatter", mode = "lines",
            line = list(color = "#8E44AD", width = 2)) %>%
      layout(title = "Collagen Index", xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })
  output$plt_loxl2 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~LOXL2, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "LOXL2 (Cross-linking)", xaxis = list(title = "Time (years)"), yaxis = list(title = "AU"))
  })
  output$plt_fibroscan <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~Fibroscan, type = "scatter", mode = "lines",
            line = list(color = "#2C3E50", width = 2)) %>%
      add_lines(y = rep(9.5, nrow(d)), name = "F2 threshold (9.5 kPa)",
                line = list(color = "#E74C3C", dash = "dash")) %>%
      layout(title = "Liver Stiffness (kPa)", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "kPa"))
  })
  output$plt_pp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~PortalPressure, type = "scatter", mode = "lines",
            line = list(color = "#C0392B", width = 2)) %>%
      add_lines(y = rep(12, nrow(d)), name = "Varices threshold (12 mmHg)",
                line = list(color = "orange", dash = "dash")) %>%
      layout(title = "Portal Pressure", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "mmHg"))
  })

  ## Outcome plots
  output$plt_alp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~ALP, type = "scatter", mode = "lines",
            line = list(color = "#2980B9", width = 2)) %>%
      add_lines(y = rep(225, nrow(d)), name = "ALP < 1.5x ULN target",
                line = list(color = "green", dash = "dash")) %>%
      add_lines(y = rep(150, nrow(d)), name = "ALP normal ULN",
                line = list(color = "gray", dash = "dot")) %>%
      layout(title = "ALP (IU/L)", xaxis = list(title = "Time (years)"), yaxis = list(title = "IU/L"))
  })
  output$plt_bili <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~Bilirubin, type = "scatter", mode = "lines",
            line = list(color = "#F39C12", width = 2)) %>%
      add_lines(y = rep(1.2, nrow(d)), name = "Normal ULN", line = list(color = "green", dash = "dash")) %>%
      layout(title = "Total Bilirubin", xaxis = list(title = "Time (years)"), yaxis = list(title = "mg/dL"))
  })
  output$plt_alp_resp <- renderPlotly({
    d <- sim_data()
    d <- d %>% mutate(ALP_resp = as.numeric(ALP < 225))
    d_monthly <- d %>%
      mutate(month = round(time_years * 12)) %>%
      group_by(month) %>%
      summarise(resp_rate = mean(ALP_resp) * 100, time_y = mean(time_years))
    plot_ly(d_monthly, x = ~time_y, y = ~resp_rate, type = "scatter", mode = "lines+markers",
            line = list(color = "#27AE60", width = 2)) %>%
      layout(title = "ALP Responder Rate (ALP<225 IU/L)",
             xaxis = list(title = "Time (years)"), yaxis = list(title = "%", range = c(0, 105)))
  })
  output$tbl_outcomes <- renderDT({
    d <- sim_data()
    yrs <- c(1, 3, 5)
    tbl <- lapply(yrs, function(y) {
      row <- d %>% filter(abs(time_years - y) == min(abs(time_years - y))) %>% head(1)
      data.frame(
        Year = y,
        ALP  = round(row$ALP, 0),
        Bili = round(row$Bilirubin, 2),
        Stiffness = round(row$Fibroscan, 1),
        PP = round(row$PortalPressure, 1),
        IL17 = round(row$IL17A, 3),
        CCA_risk = round(row$CCA_risk, 4)
      )
    })
    datatable(bind_rows(tbl), rownames = FALSE, options = list(dom = "t"))
  })

  ## CCA / Biomarker plots
  output$plt_cca <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~CCA_risk, type = "scatter", mode = "lines",
            line = list(color = "#C0392B", width = 2)) %>%
      layout(title = "CCA Risk Index", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "Risk (0-1)", range = c(0, max(d$CCA_risk) * 1.2)))
  })
  output$plt_senesc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~Senescence, type = "scatter", mode = "lines",
            line = list(color = "#E67E22", width = 2)) %>%
      layout(title = "Senescent Cholangiocytes", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "AU"))
  })
  output$plt_hsc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_years, y = ~HSC_act, type = "scatter", mode = "lines",
            line = list(color = "#8E44AD", width = 2)) %>%
      layout(title = "HSC Activation (0-1)", xaxis = list(title = "Time (years)"),
             yaxis = list(title = "Fraction", range = c(0, 1)))
  })
  output$plt_radar <- renderPlotly({
    d <- sim_data()
    end_pt <- d %>% filter(abs(time_years - min(2, max(time_years))) == min(abs(time_years - min(2, max(time_years))))) %>% head(1)
    categories <- c("ALP norm", "Bilirubin inv", "Stiffness inv", "Cholangio health",
                    "Treg", "BA control")
    values <- c(
      1 - min(1, end_pt$ALP / 500),
      1 - min(1, end_pt$Bilirubin / 10),
      1 - min(1, end_pt$Fibroscan / 25),
      end_pt$Cholangio_health,
      min(1, end_pt$Treg_IL10 * 5),
      end_pt$FXR_act
    )
    plot_ly(
      type = "scatterpolar",
      r = c(values, values[1]),
      theta = c(categories, categories[1]),
      fill = "toself",
      mode = "lines+markers"
    ) %>%
      layout(title = "Multi-Biomarker Profile at Year 2",
             polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))))
  })

  ## Dose-Response Tab
  output$plt_oca_dr <- renderPlotly({
    withProgress(message = "OCA dose-response...", {
      oca_doses <- c(0, 5, 10, 25, 50)
      dr <- lapply(oca_doses, function(dose) {
        ev_oca <- if (dose > 0) {
          ev(time = seq(0, 365 * 24 - 24, by = 24), amt = dose, cmt = "OCA_gut", evid = 1)
        } else NULL
        out <- if (is.null(ev_oca)) mrgsim(psc_mod, tgrid = seq(0, 365 * 24, by = 24 * 7))
               else mrgsim(psc_mod, events = ev_oca, tgrid = seq(0, 365 * 24, by = 24 * 7))
        as.data.frame(out) %>%
          filter(abs(time - 365 * 24) < 168) %>%
          summarise(dose = dose, ALP_Y1 = mean(ALP))
      })
      dr_df <- bind_rows(dr)
      plot_ly(dr_df, x = ~dose, y = ~ALP_Y1, type = "scatter", mode = "lines+markers",
              line = list(color = "#9B59B6", width = 2),
              marker = list(size = 8, color = "#9B59B6")) %>%
        add_lines(y = rep(225, nrow(dr_df)), name = "Target ALP",
                  line = list(color = "green", dash = "dash")) %>%
        layout(title = "OCA Dose-Response (ALP at Year 1)",
               xaxis = list(title = "OCA Dose (mg/day)"),
               yaxis = list(title = "ALP (IU/L)"))
    })
  })

  output$plt_udca_dr <- renderPlotly({
    withProgress(message = "UDCA dose-response...", {
      udca_doses <- c(0, 300, 600, 900, 1200, 1500)
      dr <- lapply(udca_doses, function(dose) {
        ev_u <- if (dose > 0) {
          ev(time = seq(0, 365 * 24 - 12, by = 12), amt = dose / 2, cmt = "UDCA_gut", evid = 1)
        } else NULL
        out <- if (is.null(ev_u)) mrgsim(psc_mod, tgrid = seq(0, 365 * 24, by = 24 * 7))
               else mrgsim(psc_mod, events = ev_u, tgrid = seq(0, 365 * 24, by = 24 * 7))
        as.data.frame(out) %>%
          filter(abs(time - 365 * 24) < 168) %>%
          summarise(dose = dose, ALP_Y1 = mean(ALP))
      })
      dr_df <- bind_rows(dr)
      plot_ly(dr_df, x = ~dose, y = ~ALP_Y1, type = "scatter", mode = "lines+markers",
              line = list(color = "#3498DB", width = 2),
              marker = list(size = 8, color = "#3498DB")) %>%
        add_lines(y = rep(225, nrow(dr_df)), name = "Target ALP",
                  line = list(color = "green", dash = "dash")) %>%
        layout(title = "UDCA Dose-Response (ALP at Year 1)",
               xaxis = list(title = "UDCA Dose (mg/day)"),
               yaxis = list(title = "ALP (IU/L)"))
    })
  })

  output$plt_ibd <- renderPlotly({
    withProgress(message = "IBD sensitivity...", {
      sims <- lapply(c(0, 1), function(ibd) {
        ev_u <- ev(time = seq(0, 3 * 365 * 24 - 12, by = 12), amt = 450, cmt = "UDCA_gut", evid = 1)
        mod2 <- param(psc_mod, list(IBD_status = ibd))
        out <- mrgsim(mod2, events = ev_u, tgrid = seq(0, 3 * 365 * 24, by = 24 * 7))
        as.data.frame(out) %>% mutate(IBD = ifelse(ibd == 1, "PSC+IBD", "PSC only"),
                                      time_years = time / (365 * 24))
      })
      d_ibd <- bind_rows(sims)
      plot_ly(d_ibd, x = ~time_years, y = ~ALP, color = ~IBD, type = "scatter", mode = "lines") %>%
        layout(title = "IBD Impact on ALP (UDCA 900 mg/day)",
               xaxis = list(title = "Time (years)"), yaxis = list(title = "ALP (IU/L)"))
    })
  })

  output$plt_sensitivity <- renderPlotly({
    withProgress(message = "Sensitivity analysis...", {
      params_test <- list(
        list(name = "k_IL17_LPS +50%", k_IL17_LPS = 0.075),
        list(name = "k_Col_synth +50%", k_Col_synth = 0.012),
        list(name = "FXR0 +50%", FXR0 = 0.45),
        list(name = "k_cholangio_repair +50%", k_cholangio_repair = 0.015),
        list(name = "IBD_status = 0", IBD_status = 0),
        list(name = "Baseline", nothing = 0)
      )
      base_alp <- {
        out <- mrgsim(psc_mod, tgrid = seq(0, 365 * 24, by = 24 * 7))
        mean(as.data.frame(out) %>% filter(abs(time - 365 * 24) < 168) %>% pull(ALP))
      }
      sens_res <- lapply(params_test, function(p) {
        nm <- p$name; p$name <- NULL
        if (length(p) == 0) {
          alp <- base_alp
        } else {
          mod2 <- param(psc_mod, p)
          out <- mrgsim(mod2, tgrid = seq(0, 365 * 24, by = 24 * 7))
          alp <- mean(as.data.frame(out) %>% filter(abs(time - 365 * 24) < 168) %>% pull(ALP))
        }
        data.frame(param = nm, ALP_pct_change = (alp - base_alp) / base_alp * 100)
      })
      df_sens <- bind_rows(sens_res)
      plot_ly(df_sens, x = ~ALP_pct_change, y = ~param,
              type = "bar", orientation = "h",
              marker = list(color = ifelse(df_sens$ALP_pct_change > 0, "#E74C3C", "#27AE60"))) %>%
        layout(title = "ALP Sensitivity at Year 1 (%Δ from baseline)",
               xaxis = list(title = "% Change in ALP vs baseline"),
               yaxis = list(title = ""))
    })
  })
}

## ---------------------------------------------------------------------------
## Launch
## ---------------------------------------------------------------------------
shinyApp(ui, server)
