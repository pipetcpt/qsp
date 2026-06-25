## ============================================================
## Beta-Thalassemia QSP Shiny Dashboard
## ============================================================
## 6 Tabs:
##   1. Patient Profile & Disease Severity
##   2. PK Profiles (Luspatercept, Deferasirox, Hydroxyurea)
##   3. Erythropoiesis Dynamics (Progenitors → Hb)
##   4. Iron Metabolism (LIC, Ferritin, Hepcidin, ERFE)
##   5. Clinical Endpoints & Scenario Comparison
##   6. Biomarker Dashboard & Correlation Explorer
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(shinydashboard)
library(DT)

## ============================================================
## MODEL CODE (inline for Shiny)
## ============================================================
bth_code <- '
$PARAM
k_BFU_prod=1.5, k_BFU_diff=0.28, k_CFU_diff=0.42
k_pro_diff=0.50, k_baso_diff=0.45, k_poly_diff=0.40
k_ortho_diff=0.55, k_retic_mat=2.0, k_rbc_elim=0.0083
ie_frac=0.70, HGB_per_RBC=28e-6, BLOOD_VOL=4.5, RBC_ref=4500.0
EPO_base=15.0, EPO_max=8000.0, HGB_target=14.0
kEPO_elim=0.35, EC50_EPO=8.0, Emax_EPO=3.5
k_ERFE_prod=0.15, k_ERFE_elim=0.35
HEPCIDIN_base=25.0, k_HEPC_prod=0.08, k_HEPC_elim=0.35
IC50_ERFE_HEPC=2.0, Emax_ERFE_HEPC=0.90
k_Fe_absorb=1.2, k_Fe_RBC=0.25, k_Fe_liver=0.04
k_Fe_release=0.01, k_FERR_form=0.005, k_FERR_elim=0.04
k_FeCARD_in=0.002, k_FeCARD_out=0.008
Ka_L=0.27, F_L=0.60, CL_L=0.35, Vc_L=8.0, Vp_L=20.0, Q_L=0.90
EC50_L=0.5, Emax_L=0.65
Ka_DFX=1.50, F_DFX=0.70, CL_DFX=14.0, V_DFX=100.0
k_DFX_Fe=0.12
Ka_HU=6.0, F_HU=0.80, CL_HU=175.0, V_HU=42.0
EC50_HU=10.0, Emax_HU=0.40

$INIT
BFU_E=5.0, CFU_E=3.5, PRO_E=7.0, BASO_E=14.0, POLY_E=14.0
ORTHO_E=14.0, RETIC=50.0, RBC_MAT=4500.0, EPO_CMT=15.0
ERFE_CMT=1.0, HEPC_CMT=25.0
FE_PL=18.0, FE_LIV=1.5, FERR_CMT=100.0, FE_CARD=0.1
LUSPAT_SC=0.0, LUSPAT_C1=0.0, LUSPAT_C2=0.0
DFX_GUT=0.0, DFX_CENT=0.0
HU_GUT=0.0, HU_CENT=0.0

$ODE
double C_LUSPAT = LUSPAT_C1 / Vc_L;
double C_DFX    = DFX_CENT  / V_DFX;
double C_HU     = HU_CENT   / V_HU;
double IE_eff_L = Emax_L * C_LUSPAT / (EC50_L + C_LUSPAT);
double ie_eff   = ie_frac * (1.0 - IE_eff_L);
double HU_HbF   = Emax_HU * C_HU / (EC50_HU + C_HU);
double ie_eff2  = ie_eff * (1.0 - 0.3 * HU_HbF);
double EPO_stim = 1.0 + Emax_EPO * EPO_CMT / (EC50_EPO + EPO_CMT);
double ERFE_sup = 1.0 - Emax_ERFE_HEPC * ERFE_CMT / (IC50_ERFE_HEPC + ERFE_CMT);
double FPN_act  = HEPCIDIN_base / (HEPCIDIN_base + HEPC_CMT);
double Fe_abs   = k_Fe_absorb * (1.0 - FPN_act);
double DFX_rem  = k_DFX_Fe * C_DFX;
double HGB_now  = (RBC_MAT / RBC_ref) * 14.0;
double POOL     = PRO_E + BASO_E + POLY_E + ORTHO_E;
double ERFE_p   = k_ERFE_prod * POOL;
double EPO_p    = EPO_base + (EPO_max - EPO_base) * pow(HGB_target, 4) /
                              (pow(HGB_target, 4) + pow(HGB_now + 0.01, 4));
double NTBI_g   = (FE_LIV > 7.0) ? (FE_LIV - 7.0) * 0.05 : 0.0;

dxdt_BFU_E   = k_BFU_prod * EPO_stim - k_BFU_diff  * BFU_E;
dxdt_CFU_E   = k_BFU_diff * BFU_E    - k_CFU_diff  * CFU_E;
dxdt_PRO_E   = k_CFU_diff * CFU_E    - k_pro_diff  * PRO_E * (1.0 + ie_eff2);
dxdt_BASO_E  = k_pro_diff * PRO_E * (1.0 - ie_eff2*0.3) - k_baso_diff * BASO_E * (1.0 + ie_eff2);
dxdt_POLY_E  = k_baso_diff * BASO_E * (1.0 - ie_eff2*0.3) - k_poly_diff * POLY_E * (1.0 + ie_eff2*0.5);
dxdt_ORTHO_E = k_poly_diff * POLY_E  * (1.0 - ie_eff2*0.15) - k_ortho_diff * ORTHO_E;
dxdt_RETIC   = k_ortho_diff * ORTHO_E - k_retic_mat * RETIC;
dxdt_RBC_MAT = k_retic_mat * RETIC   - k_rbc_elim  * RBC_MAT;
dxdt_EPO_CMT = EPO_p - kEPO_elim * EPO_CMT;
dxdt_ERFE_CMT= ERFE_p - k_ERFE_elim * ERFE_CMT;
dxdt_HEPC_CMT= k_HEPC_prod * HEPCIDIN_base * ERFE_sup - k_HEPC_elim * HEPC_CMT;
dxdt_FE_PL   = Fe_abs + k_Fe_release * FE_LIV - k_Fe_RBC * FE_PL - k_Fe_liver * FE_PL;
dxdt_FE_LIV  = k_Fe_liver * FE_PL - k_Fe_release * FE_LIV - DFX_rem * FE_LIV;
dxdt_FERR_CMT= k_FERR_form * FE_LIV - k_FERR_elim * FERR_CMT;
dxdt_FE_CARD = k_FeCARD_in * NTBI_g - k_FeCARD_out * FE_CARD;
dxdt_LUSPAT_SC = -Ka_L * LUSPAT_SC;
dxdt_LUSPAT_C1 = Ka_L * F_L * LUSPAT_SC - (CL_L + Q_L) / Vc_L * LUSPAT_C1 + Q_L / Vp_L * LUSPAT_C2;
dxdt_LUSPAT_C2 = Q_L / Vc_L * LUSPAT_C1 - Q_L / Vp_L * LUSPAT_C2;
dxdt_DFX_GUT  = -Ka_DFX * DFX_GUT;
dxdt_DFX_CENT =  Ka_DFX * F_DFX * DFX_GUT - CL_DFX / V_DFX * DFX_CENT;
dxdt_HU_GUT   = -Ka_HU * HU_GUT;
dxdt_HU_CENT  =  Ka_HU * F_HU * HU_GUT - CL_HU / V_HU * HU_CENT;

$TABLE
double Hb_gdL   = (RBC_MAT / RBC_ref) * 14.0;
double LIC_mgFe = FE_LIV;
double FERR_ugL = FERR_CMT;
double T2star   = 50.0 / (FE_CARD + 0.01);
double C_Luspa  = LUSPAT_C1 / Vc_L;
double C_DFX_pl = DFX_CENT  / V_DFX;
double C_HU_pl  = HU_CENT   / V_HU;
double IE_now   = ie_frac * (1.0 - Emax_L * C_Luspa / (EC50_L + C_Luspa));

$CAPTURE Hb_gdL LIC_mgFe FERR_ugL T2star C_Luspa C_DFX_pl C_HU_pl IE_now
         EPO_CMT ERFE_CMT HEPC_CMT BFU_E CFU_E PRO_E BASO_E POLY_E ORTHO_E RETIC RBC_MAT FE_PL FE_CARD
'

## Compile model (suppress output)
suppressMessages({
  mod_bth <- mcode("bth_shiny", bth_code)
})

## ============================================================
## Helper: build events and simulate
## ============================================================
simulate_bth <- function(
  ie_frac    = 0.70,
  epo_base   = 50,
  luspat_dose_mgkg = 0,     # 0 = no luspatercept
  dfx_dose_mgkg    = 0,
  hu_dose_mgkg     = 0,
  tx_every         = 0,     # 0 = no transfusion; else days
  wt_kg            = 70,
  duration_days    = 365
) {
  m <- param(mod_bth, ie_frac = ie_frac, EPO_base = epo_base)

  events <- NULL

  if (luspat_dose_mgkg > 0) {
    n_doses <- floor(duration_days / 21)
    ev_L <- ev(amt = luspat_dose_mgkg * wt_kg * 1000, cmt = "LUSPAT_SC",
               addl = max(0, n_doses - 1), ii = 21)
    events <- if (is.null(events)) ev_L else ev_seq(events, ev_L)
  }

  if (dfx_dose_mgkg > 0) {
    ev_D <- ev(amt = dfx_dose_mgkg * wt_kg, cmt = "DFX_GUT",
               addl = duration_days - 1, ii = 1)
    events <- if (is.null(events)) ev_D else ev_seq(events, ev_D)
  }

  if (hu_dose_mgkg > 0) {
    ev_H <- ev(amt = hu_dose_mgkg * wt_kg, cmt = "HU_GUT",
               addl = duration_days - 1, ii = 1)
    events <- if (is.null(events)) ev_H else ev_seq(events, ev_H)
  }

  if (tx_every > 0) {
    n_tx <- floor(duration_days / tx_every)
    ev_T <- ev(amt = 0.30, cmt = "FE_LIV", addl = max(0, n_tx - 1), ii = tx_every)
    events <- if (is.null(events)) ev_T else ev_seq(events, ev_T)
  }

  if (is.null(events)) {
    out <- mrgsim(m, end = duration_days, delta = 1)
  } else {
    out <- mrgsim(m, events = events, end = duration_days, delta = 1)
  }

  as_tibble(out)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(
    title = "Beta-Thalassemia QSP",
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user")),
      menuItem("PK Profiles",          tabName = "tab_pk",       icon = icon("chart-line")),
      menuItem("Erythropoiesis",       tabName = "tab_eryth",    icon = icon("droplet")),
      menuItem("Iron Metabolism",      tabName = "tab_iron",     icon = icon("atom")),
      menuItem("Clinical Endpoints",   tabName = "tab_clinical", icon = icon("stethoscope")),
      menuItem("Biomarker Dashboard",  tabName = "tab_bm",       icon = icon("microscope"))
    ),

    hr(),
    tags$div(style = "padding: 0 15px;",
      tags$h5("Global Settings", style = "color: white; margin-bottom: 5px;"),

      sliderInput("ie_frac", "IE Fraction (Disease Severity):",
                  min = 0.05, max = 0.95, value = 0.75, step = 0.05,
                  post = ""),
      helpText("0.05 = near-normal; 0.80 = severe TDT", style = "color: #ccc; font-size: 11px;"),

      sliderInput("epo_base", "Baseline EPO (IU/L):",
                  min = 15, max = 5000, value = 500, step = 50),

      sliderInput("wt_kg", "Body Weight (kg):",
                  min = 20, max = 100, value = 70, step = 5),

      sliderInput("duration", "Simulation Duration (days):",
                  min = 90, max = 730, value = 365, step = 30),

      hr(),
      tags$h5("Therapies", style = "color: white;"),

      sliderInput("luspat_dose", "Luspatercept (mg/kg SC q21d):",
                  min = 0, max = 1.75, value = 0, step = 0.25),

      sliderInput("dfx_dose", "Deferasirox (mg/kg/day PO):",
                  min = 0, max = 40, value = 0, step = 5),

      sliderInput("hu_dose", "Hydroxyurea (mg/kg/day):",
                  min = 0, max = 35, value = 0, step = 5),

      selectInput("tx_interval", "Transfusion Interval:",
                  choices = c("None (0)" = 0,
                              "Every 14 days" = 14,
                              "Every 21 days" = 21,
                              "Every 28 days" = 28),
                  selected = 0),

      actionButton("run_sim", "Run Simulation",
                   class = "btn-danger btn-block",
                   style = "margin-top: 10px;")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".box-title { font-weight: bold; }
       .small-box .icon-large { font-size: 60px; }
       .shiny-output-error { color: red; }"
    ))),

    tabItems(

      ## ========================================================
      ## TAB 1: Patient Profile
      ## ========================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(
            title = "Disease Severity Profile", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_severity_radar", height = 350),
            helpText("Spider/radar chart of disease burden metrics at simulation start vs end")
          ),
          box(
            title = "Baseline Characteristics Summary", width = 6, status = "warning", solidHeader = TRUE,
            tableOutput("tbl_baseline"),
            hr(),
            tags$h5("Genotype Reference Guide"),
            tableOutput("tbl_genotype")
          )
        ),
        fluidRow(
          box(
            title = "Beta-Thalassemia Disease Overview", width = 12, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(4,
                tags$h5("Pathophysiology"),
                tags$ul(
                  tags$li("β-globin chain deficiency → excess free α-chains"),
                  tags$li("α-chain precipitation → erythroblast apoptosis (60–90%)"),
                  tags$li("Massive ineffective erythropoiesis → severe anemia"),
                  tags$li("ERFE ↑↑↑ → hepcidin suppression → iron overload"),
                  tags$li("LIC ↑↑, cardiac T2* ↓ → organ damage")
                )
              ),
              column(4,
                tags$h5("Classification"),
                tags$ul(
                  tags$li(tags$b("TDT (Transfusion-Dependent Thalassemia):")),
                  tags$li("  β0/β0 or β0/β+ severe: Hb <7 g/dL without Tx"),
                  tags$li(tags$b("NTDT (Non-Transfusion-Dependent):")),
                  tags$li("  β+/β+, β0/β+mild, HbE/β0: Hb 7–10 g/dL"),
                  tags$li(tags$b("Thalassemia Minor (trait):")),
                  tags$li("  β/β0 or β/β+: mild microcytosis, asymptomatic")
                )
              ),
              column(4,
                tags$h5("Key Treatment Goals"),
                tags$ul(
                  tags$li("Pre-transfusion Hb ≥9.5–10.5 g/dL (TDT)"),
                  tags$li("Hb ≥9 g/dL for ≥12 consecutive wk → TI (NTDT)"),
                  tags$li("LIC <7 mg Fe/g dw (chelation target)"),
                  tags$li("Cardiac T2* >20 ms (MRI safety threshold)"),
                  tags$li("Serum ferritin <1000 µg/L (optimal)"),
                  tags$li("Growth & development normalized (pediatric)")
                )
              )
            )
          )
        )
      ),

      ## ========================================================
      ## TAB 2: PK Profiles
      ## ========================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(
            title = "Luspatercept Plasma Concentration", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_pk_luspat", height = 320)
          ),
          box(
            title = "Deferasirox Plasma Concentration", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("p_pk_dfx", height = 320)
          )
        ),
        fluidRow(
          box(
            title = "Hydroxyurea Plasma Concentration", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("p_pk_hu", height = 300)
          ),
          box(
            title = "PK Summary Table", width = 6, status = "info", solidHeader = TRUE,
            tableOutput("tbl_pk_summary"),
            hr(),
            tags$h5("PK Parameters Reference"),
            tableOutput("tbl_pk_params")
          )
        )
      ),

      ## ========================================================
      ## TAB 3: Erythropoiesis Dynamics
      ## ========================================================
      tabItem(tabName = "tab_eryth",
        fluidRow(
          box(
            title = "Erythroid Progenitor Cascade", width = 7, status = "primary", solidHeader = TRUE,
            plotOutput("p_eryth_cascade", height = 380)
          ),
          box(
            title = "Hemoglobin Time-Course", width = 5, status = "danger", solidHeader = TRUE,
            plotOutput("p_hb_timecourse", height = 380)
          )
        ),
        fluidRow(
          box(
            title = "EPO & Reticulocytes", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("p_epo_retic", height = 280)
          ),
          box(
            title = "Ineffective Erythropoiesis (IE) Rate", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_ie_rate", height = 280)
          )
        )
      ),

      ## ========================================================
      ## TAB 4: Iron Metabolism
      ## ========================================================
      tabItem(tabName = "tab_iron",
        fluidRow(
          box(
            title = "Liver Iron Content (LIC)", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("p_lic", height = 320)
          ),
          box(
            title = "Serum Ferritin", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_ferritin", height = 320)
          )
        ),
        fluidRow(
          box(
            title = "Hepcidin & ERFE Dynamics", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("p_hepc_erfe", height = 300)
          ),
          box(
            title = "Cardiac Iron (T2* proxy)", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_cardiac_iron", height = 300)
          )
        )
      ),

      ## ========================================================
      ## TAB 5: Clinical Endpoints & Scenario Comparison
      ## ========================================================
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(
            title = "Multi-Scenario Comparison — Hemoglobin", width = 8, status = "success", solidHeader = TRUE,
            plotOutput("p_compare_hb", height = 360)
          ),
          box(
            title = "Scenario Parameters", width = 4, status = "info", solidHeader = TRUE,
            tableOutput("tbl_scenarios")
          )
        ),
        fluidRow(
          box(
            title = "Multi-Scenario — LIC", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("p_compare_lic", height = 300)
          ),
          box(
            title = "Multi-Scenario — T2*", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_compare_t2star", height = 300)
          )
        ),
        fluidRow(
          box(
            title = "Clinical Endpoint Summary Table", width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("tbl_endpoint_summary")
          )
        )
      ),

      ## ========================================================
      ## TAB 6: Biomarker Dashboard
      ## ========================================================
      tabItem(tabName = "tab_bm",
        fluidRow(
          valueBoxOutput("vbox_hb",     width = 3),
          valueBoxOutput("vbox_lic",    width = 3),
          valueBoxOutput("vbox_ferrit", width = 3),
          valueBoxOutput("vbox_t2star", width = 3)
        ),
        fluidRow(
          box(
            title = "Hb vs LIC Correlation", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("p_hb_lic_corr", height = 320)
          ),
          box(
            title = "ERFE vs Hepcidin Correlation", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("p_erfe_hepc_corr", height = 320)
          )
        ),
        fluidRow(
          box(
            title = "Erythroid Cascade Waterfall (Endpoint)", width = 6, status = "info", solidHeader = TRUE,
            plotOutput("p_eryth_waterfall", height = 280)
          ),
          box(
            title = "Iron Compartment Distribution", width = 6, status = "danger", solidHeader = TRUE,
            plotOutput("p_iron_pie", height = 280)
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## -- Reactive: run single-arm simulation
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Simulating...  ", {
      simulate_bth(
        ie_frac          = input$ie_frac,
        epo_base         = input$epo_base,
        luspat_dose_mgkg = input$luspat_dose,
        dfx_dose_mgkg    = input$dfx_dose,
        hu_dose_mgkg     = input$hu_dose,
        tx_every         = as.numeric(input$tx_interval),
        wt_kg            = input$wt_kg,
        duration_days    = input$duration
      )
    })
  }, ignoreNULL = FALSE)

  ## -- Reactive: multi-scenario for comparison tab
  compare_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running 6 scenarios...", {
      scenarios <- list(
        list(label = "1. Natural History (no Rx)",
             params = list(ie_frac = input$ie_frac, epo_base = input$epo_base),
             luspat = 0, dfx = 0, hu = 0, tx = 0),
        list(label = "2. Transfusions Only",
             params = list(ie_frac = input$ie_frac, epo_base = input$epo_base),
             luspat = 0, dfx = 0, hu = 0, tx = 21),
        list(label = "3. Transfusions + DFX 30 mg/kg",
             params = list(ie_frac = input$ie_frac, epo_base = input$epo_base),
             luspat = 0, dfx = 30, hu = 0, tx = 21),
        list(label = "4. Luspatercept 1.0 mg/kg (NTDT)",
             params = list(ie_frac = max(0.40, input$ie_frac - 0.25), epo_base = max(50, input$epo_base / 3)),
             luspat = 1.0, dfx = 0, hu = 0, tx = 0),
        list(label = "5. Luspat + Tx + DFX (TDT)",
             params = list(ie_frac = input$ie_frac, epo_base = input$epo_base),
             luspat = 1.0, dfx = 30, hu = 0, tx = 21),
        list(label = "6. Gene Therapy (engrafted)",
             params = list(ie_frac = 0.05, epo_base = 15),
             luspat = 0, dfx = 0, hu = 0, tx = 0)
      )

      all_res <- lapply(scenarios, function(s) {
        res <- simulate_bth(
          ie_frac          = s$params$ie_frac,
          epo_base         = s$params$epo_base,
          luspat_dose_mgkg = s$luspat,
          dfx_dose_mgkg    = s$dfx,
          hu_dose_mgkg     = s$hu,
          tx_every         = s$tx,
          wt_kg            = input$wt_kg,
          duration_days    = input$duration
        )
        res$Scenario <- s$label
        res
      })
      bind_rows(all_res)
    })
  }, ignoreNULL = FALSE)

  ## ============================================================
  ## TAB 1 OUTPUTS
  ## ============================================================
  output$tbl_baseline <- renderTable({
    dat <- sim_result()
    hb_ss  <- round(mean(dat$Hb_gdL[dat$time > max(dat$time) * 0.8]),   2)
    lic_ss <- round(mean(dat$LIC_mgFe[dat$time > max(dat$time) * 0.8]), 2)
    ferr   <- round(mean(dat$FERR_ugL[dat$time > max(dat$time) * 0.8]), 0)
    t2s    <- round(mean(dat$T2star[dat$time > max(dat$time) * 0.8]),   1)
    epo_ss <- round(mean(dat$EPO_CMT[dat$time > max(dat$time) * 0.8]),  0)
    erfe_ss<- round(mean(dat$ERFE_CMT[dat$time > max(dat$time) * 0.8]), 2)

    data.frame(
      Parameter = c("Hemoglobin (g/dL)",
                    "LIC (mg Fe/g dw)", "Serum Ferritin (µg/L)",
                    "Cardiac T2* (ms)", "EPO (IU/L)", "ERFE (rel.)"),
      `Model (SS)` = c(hb_ss, lic_ss, ferr, t2s, epo_ss, erfe_ss),
      `Normal Range` = c("12–16", "<2", "<150", ">20", "10–30", "~1"),
      check.names = FALSE
    )
  })

  output$tbl_genotype <- renderTable({
    data.frame(
      Genotype   = c("β0/β0", "β0/β+sev", "β+/β+", "HbE/β0", "β/β0 or β/β+ (trait)"),
      `Disease Type` = c("TDT Severe", "TDT", "NTDT", "NTDT/TDT", "Thalassemia Minor"),
      `IE Fraction` = c("~0.80–0.90", "~0.70–0.80", "~0.50–0.65", "~0.60–0.75", "~0.05–0.15"),
      `Hb w/o Tx` = c("<6 g/dL", "6–7", "7–10", "7–9", ">11"),
      check.names = FALSE
    )
  })

  output$p_severity_radar <- renderPlot({
    dat <- sim_result()
    t_end <- max(dat$time)
    early <- dat %>% filter(time <= 30)
    late  <- dat %>% filter(time > t_end * 0.85)

    bm_early <- data.frame(
      Variable = c("Hb\n(g/dL)", "LIC\n(mg/g)", "EPO×10\n(IU/L)", "ERFE×5", "Ferritin/500\n(µg/L)"),
      Value    = c(
        mean(early$Hb_gdL,   na.rm=TRUE),
        mean(early$LIC_mgFe, na.rm=TRUE),
        mean(early$EPO_CMT,  na.rm=TRUE) / 10,
        mean(early$ERFE_CMT, na.rm=TRUE) * 5,
        mean(early$FERR_ugL, na.rm=TRUE) / 500
      ), Time = "Baseline"
    )
    bm_late <- data.frame(
      Variable = bm_early$Variable,
      Value    = c(
        mean(late$Hb_gdL,   na.rm=TRUE),
        mean(late$LIC_mgFe, na.rm=TRUE),
        mean(late$EPO_CMT,  na.rm=TRUE) / 10,
        mean(late$ERFE_CMT, na.rm=TRUE) * 5,
        mean(late$FERR_ugL, na.rm=TRUE) / 500
      ), Time = "End"
    )
    bm_all <- bind_rows(bm_early, bm_late)

    ggplot(bm_all, aes(x = Variable, y = Value, fill = Time, group = Time)) +
      geom_col(position = "dodge", alpha = 0.75) +
      scale_fill_manual(values = c("Baseline" = "#E53935", "End" = "#1E88E5")) +
      labs(title = "Disease Burden: Baseline vs End-of-Simulation",
           x = "", y = "Normalised Value") +
      theme_bw(base_size = 12) +
      theme(legend.position = "top")
  })

  ## ============================================================
  ## TAB 2 OUTPUTS: PK
  ## ============================================================
  output$p_pk_luspat <- renderPlot({
    dat <- sim_result()
    if (max(dat$C_Luspa, na.rm=TRUE) < 0.001) {
      ggplot() + annotate("text", x=0.5, y=0.5,
        label="No Luspatercept administered.\nSet dose > 0 in sidebar.", size=6) +
        theme_void()
    } else {
      ggplot(dat, aes(x=time, y=C_Luspa)) +
        geom_line(color="#C62828", linewidth=1) +
        geom_hline(yintercept=0.5, linetype="dashed", color="gray40") +
        annotate("text", x=max(dat$time)*0.7, y=0.6,
                 label="EC50 = 0.5 µg/mL", size=3.5) +
        labs(title="Luspatercept Plasma Concentration",
             x="Time (days)", y="Concentration (µg/mL)") +
        theme_bw(base_size=12)
    }
  })

  output$p_pk_dfx <- renderPlot({
    dat <- sim_result()
    if (max(dat$C_DFX_pl, na.rm=TRUE) < 0.001) {
      ggplot() + annotate("text", x=0.5, y=0.5,
        label="No Deferasirox administered.\nSet dose > 0 in sidebar.", size=6) +
        theme_void()
    } else {
      ggplot(dat %>% filter(time <= 14), aes(x=time, y=C_DFX_pl)) +
        geom_line(color="#EF6C00", linewidth=1) +
        labs(title="Deferasirox Plasma (first 14 days shown)",
             x="Time (days)", y="Concentration (µg/mL)") +
        theme_bw(base_size=12)
    }
  })

  output$p_pk_hu <- renderPlot({
    dat <- sim_result()
    if (max(dat$C_HU_pl, na.rm=TRUE) < 0.001) {
      ggplot() + annotate("text", x=0.5, y=0.5,
        label="No Hydroxyurea administered.\nSet dose > 0 in sidebar.", size=6) +
        theme_void()
    } else {
      ggplot(dat %>% filter(time <= 10), aes(x=time, y=C_HU_pl)) +
        geom_line(color="#7B1FA2", linewidth=1) +
        geom_hline(yintercept=10, linetype="dashed", color="gray50") +
        annotate("text", x=7, y=11, label="EC50 HbF induction = 10 µg/mL", size=3.5) +
        labs(title="Hydroxyurea Plasma (first 10 days shown)",
             x="Time (days)", y="Concentration (µg/mL)") +
        theme_bw(base_size=12)
    }
  })

  output$tbl_pk_params <- renderTable({
    data.frame(
      Drug         = c("Luspatercept", "Deferasirox", "Hydroxyurea"),
      Bioavail     = c("60%", "70%", "80%"),
      t_half       = c("~11 days", "~12 hours", "~4 hours"),
      CL           = c("0.35 L/day", "14 L/day", "175 L/day"),
      EC50_target  = c("0.5 µg/mL", "IC50 Fe = 1 nM", "EC50 HbF = 10 µg/mL"),
      check.names  = FALSE
    )
  })

  output$tbl_pk_summary <- renderTable({
    dat <- sim_result()
    data.frame(
      Drug       = c("Luspatercept", "Deferasirox", "Hydroxyurea"),
      Cmax_ugmL  = c(
        round(max(dat$C_Luspa,   na.rm=TRUE), 3),
        round(max(dat$C_DFX_pl, na.rm=TRUE), 2),
        round(max(dat$C_HU_pl,  na.rm=TRUE), 2)
      ),
      C_avg_ugmL = c(
        round(mean(dat$C_Luspa[dat$C_Luspa > 0.001],   na.rm=TRUE), 3),
        round(mean(dat$C_DFX_pl[dat$C_DFX_pl > 0.001], na.rm=TRUE), 2),
        round(mean(dat$C_HU_pl[dat$C_HU_pl > 0.001],   na.rm=TRUE), 2)
      ),
      check.names = FALSE
    )
  })

  ## ============================================================
  ## TAB 3 OUTPUTS: Erythropoiesis
  ## ============================================================
  output$p_hb_timecourse <- renderPlot({
    dat <- sim_result()
    ggplot(dat, aes(x=time, y=Hb_gdL)) +
      geom_line(color="#C62828", linewidth=1.2) +
      geom_hline(yintercept=9.5, linetype="dashed", color="orange") +
      geom_hline(yintercept=7.0, linetype="dotted", color="red") +
      annotate("text", x=max(dat$time)*0.6, y=10.2,
               label="Pre-Tx target ≥9.5 g/dL", size=3.5, color="orange4") +
      annotate("text", x=max(dat$time)*0.6, y=7.7,
               label="Tx threshold <7 g/dL", size=3.5, color="red") +
      labs(title="Hemoglobin Time-Course", x="Time (days)", y="Hb (g/dL)") +
      coord_cartesian(ylim=c(0, 16)) +
      theme_bw(base_size=12)
  })

  output$p_eryth_cascade <- renderPlot({
    dat <- sim_result()
    stages_long <- dat %>%
      select(time, BFU_E, CFU_E, PRO_E, BASO_E, POLY_E, ORTHO_E, RETIC) %>%
      pivot_longer(-time, names_to="Stage", values_to="Count") %>%
      mutate(Stage = factor(Stage,
                            levels=c("BFU_E","CFU_E","PRO_E","BASO_E","POLY_E","ORTHO_E","RETIC"),
                            labels=c("BFU-E","CFU-E","Pro-EB","Baso-EB","Poly-EB","Ortho-EB","Retic")))

    ggplot(stages_long, aes(x=time, y=Count, color=Stage)) +
      geom_line(linewidth=0.9) +
      scale_color_brewer(palette="Spectral") +
      labs(title="Erythroid Progenitor Cascade",
           x="Time (days)", y="Cells (per µL)") +
      theme_bw(base_size=12) +
      theme(legend.position="right")
  })

  output$p_epo_retic <- renderPlot({
    dat <- sim_result()
    p1 <- ggplot(dat, aes(x=time, y=EPO_CMT)) +
      geom_line(color="#1565C0", linewidth=1) +
      labs(y="EPO (IU/L)") +
      theme_bw(base_size=11)
    p2 <- ggplot(dat, aes(x=time, y=RETIC)) +
      geom_line(color="#2E7D32", linewidth=1) +
      labs(y="Reticulocytes (cells/µL)") +
      theme_bw(base_size=11)
    p1 / p2 + plot_annotation(title="EPO & Reticulocytes")
  })

  output$p_ie_rate <- renderPlot({
    dat <- sim_result()
    ggplot(dat, aes(x=time, y=IE_now)) +
      geom_line(color="#B71C1C", linewidth=1.2) +
      geom_hline(yintercept=input$ie_frac, linetype="dashed", color="gray50") +
      annotate("text", x=max(dat$time)*0.6, y=input$ie_frac + 0.03,
               label="Untreated IE fraction", size=3.5) +
      labs(title="Effective IE Fraction (drug-modified)",
           x="Time (days)", y="IE Fraction (0–1)") +
      coord_cartesian(ylim=c(0, 1)) +
      theme_bw(base_size=12)
  })

  ## ============================================================
  ## TAB 4 OUTPUTS: Iron
  ## ============================================================
  output$p_lic <- renderPlot({
    dat <- sim_result()
    ggplot(dat, aes(x=time, y=LIC_mgFe)) +
      geom_line(color="#E65100", linewidth=1.2) +
      geom_hline(yintercept=7,  linetype="dashed",  color="orange") +
      geom_hline(yintercept=15, linetype="dotted",  color="red") +
      annotate("text", x=max(dat$time)*0.5, y=7.8, label="Chelation threshold: 7 mg/g",  size=3.5, color="orange4") +
      annotate("text", x=max(dat$time)*0.5, y=15.8,label="High risk: 15 mg/g", size=3.5, color="red") +
      labs(title="Liver Iron Content (LIC)", x="Time (days)", y="LIC (mg Fe/g dw)") +
      theme_bw(base_size=12)
  })

  output$p_ferritin <- renderPlot({
    dat <- sim_result()
    ggplot(dat, aes(x=time, y=FERR_ugL)) +
      geom_line(color="#F57C00", linewidth=1.2) +
      geom_hline(yintercept=1000, linetype="dashed", color="orange") +
      geom_hline(yintercept=2500, linetype="dotted", color="red") +
      annotate("text", x=max(dat$time)*0.5, y=1100, label="Target: <1000 µg/L", size=3.5, color="orange4") +
      labs(title="Serum Ferritin", x="Time (days)", y="Ferritin (µg/L)") +
      theme_bw(base_size=12)
  })

  output$p_hepc_erfe <- renderPlot({
    dat <- sim_result()
    p1 <- ggplot(dat, aes(x=time, y=ERFE_CMT)) +
      geom_line(color="#880E4F", linewidth=1) +
      labs(title="Erythroferrone (ERFE)", y="ERFE (relative units)") +
      theme_bw(base_size=11)
    p2 <- ggplot(dat, aes(x=time, y=HEPC_CMT)) +
      geom_line(color="#4A148C", linewidth=1) +
      labs(title="Hepcidin", x="Time (days)", y="Hepcidin (nM)") +
      theme_bw(base_size=11)
    p1 / p2
  })

  output$p_cardiac_iron <- renderPlot({
    dat <- sim_result()
    ggplot(dat, aes(x=time, y=T2star)) +
      geom_line(color="#B71C1C", linewidth=1.2) +
      geom_hline(yintercept=20, linetype="dashed", color="red") +
      annotate("text", x=max(dat$time)*0.5, y=19,
               label="T2* < 20 ms = cardiac iron overload risk", size=3.5, color="red") +
      labs(title="Cardiac T2* (MRI proxy, higher = safer)",
           x="Time (days)", y="T2* (ms)") +
      theme_bw(base_size=12)
  })

  ## ============================================================
  ## TAB 5 OUTPUTS: Clinical Scenarios
  ## ============================================================
  output$p_compare_hb <- renderPlot({
    dat <- compare_data()
    ggplot(dat, aes(x=time, y=Hb_gdL, color=Scenario, linetype=Scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=9.5, linetype="dashed", color="gray40") +
      scale_color_brewer(palette="Dark2") +
      labs(title="Hemoglobin — 6-Scenario Comparison",
           x="Time (days)", y="Hb (g/dL)") +
      theme_bw(base_size=12) +
      theme(legend.position="bottom", legend.title=element_blank()) +
      guides(color=guide_legend(nrow=3), linetype=guide_legend(nrow=3))
  })

  output$p_compare_lic <- renderPlot({
    dat <- compare_data()
    ggplot(dat, aes(x=time, y=LIC_mgFe, color=Scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=7, linetype="dashed", color="orange") +
      scale_color_brewer(palette="Dark2") +
      labs(title="LIC Comparison", x="Time (days)", y="LIC (mg Fe/g dw)") +
      theme_bw(base_size=11) +
      theme(legend.position="bottom", legend.title=element_blank())
  })

  output$p_compare_t2star <- renderPlot({
    dat <- compare_data()
    ggplot(dat, aes(x=time, y=T2star, color=Scenario)) +
      geom_line(linewidth=0.9) +
      geom_hline(yintercept=20, linetype="dashed", color="red") +
      scale_color_brewer(palette="Dark2") +
      labs(title="Cardiac T2* Comparison", x="Time (days)", y="T2* (ms)") +
      theme_bw(base_size=11) +
      theme(legend.position="bottom", legend.title=element_blank())
  })

  output$tbl_scenarios <- renderTable({
    data.frame(
      `#` = 1:6,
      Description = c("Natural History",
                       "Transfusions q21d",
                       "Tx + DFX 30 mg/kg",
                       "Luspat 1.0 (NTDT)",
                       "Luspat + Tx + DFX",
                       "Gene Therapy"),
      check.names = FALSE
    )
  })

  output$tbl_endpoint_summary <- renderDT({
    dat <- compare_data()
    t_end <- max(dat$time)
    dat %>%
      filter(time > t_end * 0.8) %>%
      group_by(Scenario) %>%
      summarise(
        `Hb (g/dL)`    = round(mean(Hb_gdL,   na.rm=TRUE), 2),
        `LIC (mg/g)`   = round(mean(LIC_mgFe,  na.rm=TRUE), 1),
        `Ferritin (µg/L)` = round(mean(FERR_ugL,na.rm=TRUE), 0),
        `T2* (ms)`     = round(mean(T2star,     na.rm=TRUE), 1),
        `EPO (IU/L)`   = round(mean(EPO_CMT,    na.rm=TRUE), 0),
        `Hepcidin (nM)`= round(mean(HEPC_CMT,   na.rm=TRUE), 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength=6, dom='t'), rownames=FALSE)
  })

  ## ============================================================
  ## TAB 6 OUTPUTS: Biomarker Dashboard
  ## ============================================================
  get_ss <- reactive({
    dat <- sim_result()
    t_end <- max(dat$time)
    dat %>% filter(time > t_end * 0.8) %>%
      summarise(across(c(Hb_gdL, LIC_mgFe, FERR_ugL, T2star), mean, na.rm=TRUE))
  })

  output$vbox_hb <- renderValueBox({
    ss <- get_ss()
    color <- ifelse(ss$Hb_gdL >= 9.5, "green", ifelse(ss$Hb_gdL >= 7, "yellow", "red"))
    valueBox(round(ss$Hb_gdL, 1), "Hb (g/dL)", icon=icon("droplet"), color=color)
  })
  output$vbox_lic <- renderValueBox({
    ss <- get_ss()
    color <- ifelse(ss$LIC_mgFe <= 7, "green", ifelse(ss$LIC_mgFe <= 15, "yellow", "red"))
    valueBox(round(ss$LIC_mgFe, 1), "LIC (mg Fe/g dw)", icon=icon("liver"), color=color)
  })
  output$vbox_ferrit <- renderValueBox({
    ss <- get_ss()
    color <- ifelse(ss$FERR_ugL <= 1000, "green", ifelse(ss$FERR_ugL <= 2500, "yellow", "red"))
    valueBox(round(ss$FERR_ugL, 0), "Ferritin (µg/L)", icon=icon("flask"), color=color)
  })
  output$vbox_t2star <- renderValueBox({
    ss <- get_ss()
    color <- ifelse(ss$T2star >= 20, "green", ifelse(ss$T2star >= 10, "yellow", "red"))
    valueBox(round(ss$T2star, 1), "Cardiac T2* (ms)", icon=icon("heart"), color=color)
  })

  output$p_hb_lic_corr <- renderPlot({
    dat <- sim_result()
    ggplot(dat[seq(1, nrow(dat), by=5),], aes(x=LIC_mgFe, y=Hb_gdL, color=time)) +
      geom_point(alpha=0.6, size=1.5) +
      geom_smooth(method="loess", se=FALSE, color="red") +
      scale_color_viridis_c(option="plasma") +
      labs(title="Hb vs Liver Iron Content", x="LIC (mg Fe/g dw)", y="Hb (g/dL)",
           color="Day") +
      theme_bw(base_size=12)
  })

  output$p_erfe_hepc_corr <- renderPlot({
    dat <- sim_result()
    ggplot(dat[seq(1, nrow(dat), by=5),], aes(x=ERFE_CMT, y=HEPC_CMT, color=time)) +
      geom_point(alpha=0.6, size=1.5) +
      geom_smooth(method="loess", se=FALSE, color="blue") +
      scale_color_viridis_c(option="magma") +
      labs(title="ERFE vs Hepcidin (Inverse Relationship)",
           x="ERFE (relative units)", y="Hepcidin (nM)", color="Day") +
      theme_bw(base_size=12)
  })

  output$p_eryth_waterfall <- renderPlot({
    dat <- sim_result()
    t_end <- max(dat$time)
    ep <- dat %>% filter(time > t_end * 0.85) %>%
      summarise(across(c(BFU_E, CFU_E, PRO_E, BASO_E, POLY_E, ORTHO_E, RETIC), mean, na.rm=TRUE)) %>%
      pivot_longer(everything(), names_to="Stage", values_to="Count") %>%
      mutate(Stage = factor(Stage,
                            levels=c("BFU_E","CFU_E","PRO_E","BASO_E","POLY_E","ORTHO_E","RETIC"),
                            labels=c("BFU-E","CFU-E","Pro-EB","Baso-EB","Poly-EB","Ortho-EB","Retic")))

    ggplot(ep, aes(x=Stage, y=Count, fill=Stage)) +
      geom_col(show.legend=FALSE) +
      scale_fill_brewer(palette="Spectral") +
      labs(title="Erythroid Cascade (Steady-State)", x="", y="Cells/µL") +
      theme_bw(base_size=12) +
      theme(axis.text.x=element_text(angle=30, hjust=1))
  })

  output$p_iron_pie <- renderPlot({
    dat <- sim_result()
    t_end <- max(dat$time)
    iron_ss <- dat %>% filter(time > t_end * 0.85) %>%
      summarise(Plasma_Fe = mean(FE_PL, na.rm=TRUE),
                Liver_Fe  = mean(LIC_mgFe * 100, na.rm=TRUE),
                Cardiac_Fe= mean(FE_CARD * 20,   na.rm=TRUE))
    iron_df <- data.frame(
      Compartment = c("Plasma", "Liver", "Cardiac"),
      Value       = c(iron_ss$Plasma_Fe, iron_ss$Liver_Fe, iron_ss$Cardiac_Fe)
    )
    ggplot(iron_df, aes(x="", y=Value, fill=Compartment)) +
      geom_col(width=1) +
      coord_polar("y") +
      scale_fill_manual(values=c("Plasma"="#FF9800","Liver"="#E64A19","Cardiac"="#B71C1C")) +
      labs(title="Iron Distribution (Relative)") +
      theme_void(base_size=12) +
      theme(legend.position="right")
  })
}

## ============================================================
## LAUNCH APP
## ============================================================
shinyApp(ui, server)
