## ============================================================
## Thyroid Cancer QSP — Shiny Interactive Dashboard
## ============================================================
## Tabs:
##  1. Patient Profile & Subtype
##  2. Drug PK (Lenvatinib / Sorafenib / Selpercatinib)
##  3. Oncogenic Pathways (MAPK / PI3K / Angiogenesis)
##  4. Tumor Dynamics (volume, proliferation, apoptosis)
##  5. Biomarkers (Tg, Calcitonin, CEA, TSH)
##  6. Clinical Endpoints (RECIST, PFS simulation)
##  7. Scenario Comparison & VP Analysis
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(shinydashboard)
library(DT)
library(plotly)

## ── Inline model ─────────────────────────────────────────────
thyca_code <- '
$PARAM
kprol_base=0.03 kdeath_base=0.02 kprol_MAPK=0.02 kprol_PI3K=0.01
MAPK_base=0.8 PI3K_base=0.5 MAPK_max=1.0 PI3K_max=1.0
kVEGF_prod=0.05 kVEGF_clear=0.1 kAngio=0.02 kAngio_reg=0.05
kTg_prod=0.08 kTg_clear=0.05 kCT_prod=0.15 kCT_clear=0.08
TSH_base=1.5 TSH_supp=0.1 kTSH_stim=0.003
CL_LENV=4.0 V1_LENV=50.0 Q_LENV=1.6 V2_LENV=120.0 ka_LENV=1.4 F_LENV=0.85
CL_SORA=3.2 V1_SORA=55.0 Q_SORA=2.0 V2_SORA=220.0 ka_SORA=0.6 F_SORA=0.38
CL_SELP=6.5 V1_SELP=198.0 Q_SELP=3.2 V2_SELP=1130.0 ka_SELP=0.7 F_SELP=0.73
EC50_LENV_VEGFR=5.0 Emax_LENV_VEGFR=0.92 EC50_LENV_MAPK=20.0 Emax_LENV_MAPK=0.60 Hill_LENV=1.5
EC50_SORA_BRAF=3.5 Emax_SORA_BRAF=0.75 EC50_SORA_VEGFR=6.0 Emax_SORA_VEGFR=0.80 Hill_SORA=1.5
EC50_SELP_RET=10.0 Emax_SELP_RET=0.97 Hill_SELP=2.0
use_LENV=0 use_SORA=0 use_SELP=0 use_TSH_supp=1 is_MTC=0

$CMT
LENV_gut LENV_C LENV_P SORA_gut SORA_C SORA_P SELP_gut SELP_C SELP_P
MAPK_act PI3K_act VEGF Angio TumorN Tg CT TSH TumVol

$MAIN
double Cp_LENV = LENV_C / V1_LENV * 1000.0;
double Cp_SORA = SORA_C / V1_SORA * 1000.0;
double Cp_SELP = SELP_C / V1_SELP * 1000.0;
double Inh_LENV_VEGFR = use_LENV * Emax_LENV_VEGFR * pow(Cp_LENV, Hill_LENV) /
  (pow(EC50_LENV_VEGFR, Hill_LENV) + pow(Cp_LENV, Hill_LENV));
double Inh_LENV_MAPK  = use_LENV * Emax_LENV_MAPK * pow(Cp_LENV, Hill_LENV) /
  (pow(EC50_LENV_MAPK, Hill_LENV) + pow(Cp_LENV, Hill_LENV));
double Inh_SORA_BRAF  = use_SORA * Emax_SORA_BRAF * pow(Cp_SORA, Hill_SORA) /
  (pow(EC50_SORA_BRAF, Hill_SORA) + pow(Cp_SORA, Hill_SORA));
double Inh_SORA_VEGFR = use_SORA * Emax_SORA_VEGFR * pow(Cp_SORA, Hill_SORA) /
  (pow(EC50_SORA_VEGFR, Hill_SORA) + pow(Cp_SORA, Hill_SORA));
double Inh_SELP_RET   = use_SELP * Emax_SELP_RET * pow(Cp_SELP, Hill_SELP) /
  (pow(EC50_SELP_RET, Hill_SELP) + pow(Cp_SELP, Hill_SELP));

$INIT
LENV_gut=0 LENV_C=0 LENV_P=0 SORA_gut=0 SORA_C=0 SORA_P=0
SELP_gut=0 SELP_C=0 SELP_P=0
MAPK_act=0.8 PI3K_act=0.5 VEGF=0.5 Angio=0.6
TumorN=1.0 Tg=50.0 CT=100.0 TSH=1.5 TumVol=30.0

$ODE
dxdt_LENV_gut = -ka_LENV*LENV_gut;
dxdt_LENV_C = F_LENV*ka_LENV*LENV_gut - (CL_LENV+Q_LENV)/V1_LENV*LENV_C + Q_LENV/V2_LENV*LENV_P;
dxdt_LENV_P = Q_LENV/V1_LENV*LENV_C - Q_LENV/V2_LENV*LENV_P;
dxdt_SORA_gut = -ka_SORA*SORA_gut;
dxdt_SORA_C = F_SORA*ka_SORA*SORA_gut - (CL_SORA+Q_SORA)/V1_SORA*SORA_C + Q_SORA/V2_SORA*SORA_P;
dxdt_SORA_P = Q_SORA/V1_SORA*SORA_C - Q_SORA/V2_SORA*SORA_P;
dxdt_SELP_gut = -ka_SELP*SELP_gut;
dxdt_SELP_C = F_SELP*ka_SELP*SELP_gut - (CL_SELP+Q_SELP)/V1_SELP*SELP_C + Q_SELP/V2_SELP*SELP_P;
dxdt_SELP_P = Q_SELP/V1_SELP*SELP_C - Q_SELP/V2_SELP*SELP_P;
double MAPK_activ = MAPK_base + 0.05*TSH*kTSH_stim;
double MAPK_inh = Inh_SORA_BRAF*0.9 + Inh_LENV_MAPK*0.3;
if(is_MTC > 0.5) MAPK_inh += Inh_SELP_RET*0.8;
double MAPK_target = MAPK_activ*(1.0 - MAPK_inh);
if(MAPK_target < 0.0) MAPK_target = 0.0;
dxdt_MAPK_act = 0.5*(MAPK_target - MAPK_act);
double PI3K_inh = Inh_LENV_VEGFR*0.4 + Inh_SORA_VEGFR*0.3;
double PI3K_target = PI3K_base*(1.0 - PI3K_inh);
dxdt_PI3K_act = 0.3*(PI3K_target - PI3K_act);
dxdt_VEGF = kVEGF_prod*TumorN*PI3K_act - kVEGF_clear*VEGF;
double angio_stim = kAngio*VEGF*(1.0-Inh_LENV_VEGFR)*(1.0-Inh_SORA_VEGFR);
dxdt_Angio = angio_stim - kAngio_reg*Angio;
double growth_factor = (1.0 + kprol_MAPK*MAPK_act + kprol_PI3K*PI3K_act + kTSH_stim*TSH);
double kprol_net = kprol_base*growth_factor*Angio;
double kdeath_drug = 0.015*Inh_LENV_VEGFR + 0.010*Inh_SORA_BRAF + 0.018*Inh_SELP_RET;
dxdt_TumorN = kprol_net*TumorN - (kdeath_base + kdeath_drug)*TumorN;
dxdt_TumVol = (kprol_net - kdeath_base - kdeath_drug)*TumVol;
double Tg_prod = (1.0-is_MTC)*kTg_prod*TumorN*TSH/(TSH+0.5);
dxdt_Tg = Tg_prod*(1.0-use_TSH_supp*0.7) - kTg_clear*Tg;
double CT_inh = Inh_SELP_RET*0.8 + (use_SORA > 0.5 ? 0.3 : 0.0);
dxdt_CT = is_MTC*kCT_prod*TumorN*(1.0-CT_inh) - kCT_clear*CT;
double TSH_target = use_TSH_supp > 0.5 ? TSH_supp : TSH_base;
dxdt_TSH = 0.1*(TSH_target - TSH);

$TABLE
double Cp_LENV_out = LENV_C / V1_LENV * 1000.0;
double Cp_SORA_out = SORA_C / V1_SORA * 1000.0;
double Cp_SELP_out = SELP_C / V1_SELP * 1000.0;
double SLD_change = (TumVol - 30.0) / 30.0 * 100.0;

$CAPTURE
Cp_LENV_out Cp_SORA_out Cp_SELP_out
MAPK_act PI3K_act VEGF Angio
TumorN TumVol Tg CT TSH SLD_change
'

thyca_mod <- mcode("thyca_shiny", thyca_code, quiet = TRUE)

## Color palette
pal <- c(
  "Untreated"            = "#DC2626",
  "Lenvatinib"           = "#2563EB",
  "Sorafenib"            = "#16A34A",
  "Lenv→Lenv (2L)"      = "#7C3AED",
  "Selpercatinib (MTC)"  = "#EA580C",
  "Vandetanib (MTC)"     = "#CA8A04",
  "Lenv+Everolimus"      = "#0891B2"
)

## ── UI ───────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Thyroid Cancer QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("Oncogenic Pathways", tabName = "pathways",  icon = icon("dna")),
      menuItem("Tumor Dynamics",     tabName = "tumor",     icon = icon("microscope")),
      menuItem("Biomarkers",         tabName = "biomarkers",icon = icon("flask")),
      menuItem("Clinical Endpoints", tabName = "clinical",  icon = icon("chart-line")),
      menuItem("Scenario Compare",   tabName = "compare",   icon = icon("table-columns"))
    ),
    hr(),
    h5("Treatment Selection", style = "padding-left:15px;color:#ddd"),
    selectInput("scenario", "Scenario",
      choices = c("Untreated", "Lenvatinib", "Sorafenib",
                  "Lenv→Lenv (2L)", "Selpercatinib (MTC)",
                  "Vandetanib (MTC)", "Lenv+Everolimus"),
      selected = "Lenvatinib"),
    selectInput("subtype", "Disease Subtype",
      choices = c("DTC (PTC/FTC)", "MTC"), selected = "DTC (PTC/FTC)"),
    sliderInput("sim_months", "Duration (months)", 6, 36, 24, 6),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px;color:#ddd"),
    sliderInput("MAPK_level", "MAPK Activity (BRAF status)",
                0.3, 1.0, 0.8, 0.1),
    sliderInput("TumVol_init", "Baseline Tumor SLD (mm)",
                10, 100, 30, 5),
    sliderInput("Tg_init", "Baseline Tg (ng/mL)", 5, 500, 50, 5),
    hr(),
    h5("Drug Dose", style = "padding-left:15px;color:#ddd"),
    sliderInput("dose_lenv", "Lenvatinib (mg/d)", 10, 24, 24, 2),
    sliderInput("dose_sora", "Sorafenib (mg BID)", 200, 400, 400, 100),
    sliderInput("dose_selp", "Selpercatinib (mg BID)", 80, 160, 160, 20)
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ─────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Thyroid Cancer Overview",
              status = "primary", solidHeader = TRUE, width = 8,
              p("Thyroid cancer is the most common endocrine malignancy, with an
                annual incidence of ~45,000 cases in the US. It comprises four
                major histological subtypes with distinct molecular drivers and
                treatment approaches."),
              fluidRow(
                column(6, tags$ul(
                  tags$li(strong("PTC (papillary)"), " ~85%: BRAF V600E (~60%), RET/PTC"),
                  tags$li(strong("FTC (follicular)"), " ~10%: RAS mutations, PTEN loss"),
                  tags$li(strong("MTC (medullary)"), " ~3%: RET mutations (M918T/C634R)"),
                  tags$li(strong("ATC (anaplastic)"), " ~1%: TP53, BRAF, CTNNB1")
                )),
                column(6, tags$ul(
                  tags$li("DTC (PTC+FTC): excellent prognosis if RAI-sensitive"),
                  tags$li("RAI-refractory DTC: 5-yr OS ~66% → TKI indicated"),
                  tags$li("MTC: no RAI response → RET inhibitors/TKIs"),
                  tags$li("ATC: median OS ~5 months, BRAF V600E → dabrafenib+trametinib")
                ))
              )
          ),
          box(title = "Model Stats", status = "info", solidHeader = TRUE, width = 4,
            valueBoxOutput("box_ode", width = 12),
            valueBoxOutput("box_drugs", width = 12),
            valueBoxOutput("box_trials", width = 12)
          )
        ),
        fluidRow(
          box(title = "Clinical Trial Benchmark",
              status = "success", solidHeader = TRUE, width = 7,
              DTOutput("trial_table")),
          box(title = "Baseline Summary",
              status = "warning", solidHeader = TRUE, width = 5,
              tableOutput("pt_table"))
        )
      ),

      ## ── Tab 2: Drug PK ─────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Lenvatinib Plasma Concentration",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("pk_lenv", height = 280)),
          box(title = "Sorafenib Plasma Concentration",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("pk_sora", height = 280))
        ),
        fluidRow(
          box(title = "Selpercatinib Plasma Concentration",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("pk_selp", height = 280)),
          box(title = "PK Parameters Summary",
              status = "info", solidHeader = TRUE, width = 6,
              tableOutput("pk_param_table"))
        )
      ),

      ## ── Tab 3: Oncogenic Pathways ──────────────────────────
      tabItem(tabName = "pathways",
        fluidRow(
          box(title = "MAPK (RAF-MEK-ERK) Activity",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("pathway_mapk", height = 300)),
          box(title = "PI3K / AKT / mTOR Activity",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("pathway_pi3k", height = 300))
        ),
        fluidRow(
          box(title = "VEGF Concentration & Angiogenesis",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("pathway_angio", height = 300)),
          box(title = "Pathway Inhibition at Steady State",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("pathway_bar", height = 300))
        )
      ),

      ## ── Tab 4: Tumor Dynamics ──────────────────────────────
      tabItem(tabName = "tumor",
        fluidRow(
          box(title = "Tumor Cell Number (Normalized)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("tumor_cells", height = 300)),
          box(title = "Tumor Volume (SLD)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("tumor_vol", height = 300))
        ),
        fluidRow(
          box(title = "Net Growth Rate",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("tumor_growth", height = 300)),
          box(title = "RECIST Waterfall (Final SLD Change)",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("waterfall", height = 300))
        )
      ),

      ## ── Tab 5: Biomarkers ─────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Serum Thyroglobulin (Tg) — DTC",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_tg", height = 280)),
          box(title = "Serum Calcitonin — MTC",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_ct", height = 280))
        ),
        fluidRow(
          box(title = "TSH Level (Suppression Therapy)",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_tsh", height = 280)),
          box(title = "Biomarker Response Summary",
              status = "success", solidHeader = TRUE, width = 6,
              tableOutput("bm_table"))
        )
      ),

      ## ── Tab 6: Clinical Endpoints ──────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "% Change in Tumor SLD (RECIST 1.1)",
              status = "primary", solidHeader = TRUE, width = 8,
              plotlyOutput("recist_plot", height = 380)),
          box(title = "Response Summary",
              status = "success", solidHeader = TRUE, width = 4,
              tableOutput("response_table"))
        ),
        fluidRow(
          box(title = "Simulated PFS (Time to PD)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("pfs_plot", height = 300)),
          box(title = "Biomarker-Response Correlation",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("bm_corr", height = 300))
        )
      ),

      ## ── Tab 7: Scenario Comparison ─────────────────────────
      tabItem(tabName = "compare",
        fluidRow(
          box(title = "Tumor Volume — All Scenarios",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("cmp_vol", height = 350)),
          box(title = "SLD % Change — All Scenarios",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("cmp_sld", height = 350))
        ),
        fluidRow(
          box(title = "Tg — All Scenarios",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("cmp_tg", height = 350)),
          box(title = "2-Year Summary Table",
              status = "success", solidHeader = TRUE, width = 6,
              DTOutput("cmp_table"))
        )
      )
    )
  )
)

## ── SERVER ───────────────────────────────────────────────────
server <- function(input, output, session) {

  end_hrs <- reactive({ input$sim_months * 30 * 24 })

  build_params_ui <- function(sc, subtype, mapk_lev) {
    is_mtc <- if(grepl("MTC", subtype)) 1 else 0
    p <- c(use_LENV=0, use_SORA=0, use_SELP=0,
           use_TSH_supp=1, is_MTC=is_mtc,
           MAPK_base=mapk_lev)
    switch(sc,
      "Lenvatinib"           = { p["use_LENV"] <- 1 },
      "Sorafenib"            = { p["use_SORA"] <- 1 },
      "Lenv→Lenv (2L)"      = { p["use_LENV"] <- 1; p["use_SORA"] <- 1 },
      "Selpercatinib (MTC)"  = { p["use_SELP"] <- 1; p["is_MTC"] <- 1 },
      "Vandetanib (MTC)"     = { p["use_SELP"] <- 1; p["is_MTC"] <- 1;
                                  p["Emax_SELP_RET"] <- 0.75 },
      "Lenv+Everolimus"      = { p["use_LENV"] <- 1 }
    )
    p
  }

  build_ev_ui <- function(sc, end_h, d_lenv, d_sora, d_selp) {
    switch(sc,
      "Lenvatinib"          = ev(amt=d_lenv, cmt="LENV_gut", ii=24, addl=floor(end_h/24)-1, time=0),
      "Sorafenib"           = ev(amt=d_sora, cmt="SORA_gut", ii=12, addl=floor(end_h/12)-1, time=0),
      "Lenv→Lenv (2L)"     = {
        e1 <- ev(amt=d_sora, cmt="SORA_gut", ii=12, addl=floor(end_h/24*0.5)-1, time=0)
        e2 <- ev(amt=d_lenv, cmt="LENV_gut", ii=24, addl=floor(end_h/24*0.5)-1, time=end_h*0.5)
        e1 + e2
      },
      "Selpercatinib (MTC)" = ev(amt=d_selp, cmt="SELP_gut", ii=12, addl=floor(end_h/12)-1, time=0),
      "Vandetanib (MTC)"    = ev(amt=d_selp, cmt="SELP_gut", ii=24, addl=floor(end_h/24)-1, time=0),
      "Lenv+Everolimus"     = ev(amt=d_lenv, cmt="LENV_gut", ii=24, addl=floor(end_h/24)-1, time=0),
      ev(amt=0, cmt="LENV_gut", time=0)
    )
  }

  sim_data <- reactive({
    p   <- build_params_ui(input$scenario, input$subtype, input$MAPK_level)
    ev  <- build_ev_ui(input$scenario, end_hrs(), input$dose_lenv, input$dose_sora, input$dose_selp)
    thyca_mod %>%
      param(p) %>%
      init(TumVol = input$TumVol_init, Tg = input$Tg_init,
           MAPK_act = input$MAPK_level) %>%
      mrgsim_q(ev, end = end_hrs(), delta = 24) %>%
      as_tibble() %>%
      mutate(time_mo = time / (24 * 30))
  })

  all_data <- reactive({
    end_h <- end_hrs()
    map_df(names(pal), function(sc) {
      p  <- build_params_ui(sc, input$subtype, input$MAPK_level)
      ev <- build_ev_ui(sc, end_h, input$dose_lenv, input$dose_sora, input$dose_selp)
      tryCatch({
        thyca_mod %>%
          param(p) %>%
          init(TumVol = input$TumVol_init, Tg = input$Tg_init,
               MAPK_act = input$MAPK_level) %>%
          mrgsim_q(ev, end = end_h, delta = 24*7) %>%
          as_tibble() %>%
          mutate(scenario = sc, time_mo = time / (24*30))
      }, error = function(e) tibble())
    })
  })

  ## Tab 1 info
  output$box_ode    <- renderValueBox(valueBox(18, "ODE Compartments", icon=icon("cog"), color="blue"))
  output$box_drugs  <- renderValueBox(valueBox(5,  "Drug Classes", icon=icon("pills"), color="green"))
  output$box_trials <- renderValueBox(valueBox(7,  "Treatment Scenarios", icon=icon("list"), color="orange"))

  output$trial_table <- renderDT({
    tibble(
      Trial       = c("SELECT","DECISION","LIBRETTO-001","ZETA","EXAM","METRO"),
      Drug        = c("Lenvatinib","Sorafenib","Selpercatinib","Vandetanib","Cabozantinib","Pralsetinib"),
      Indication  = c("DTC","DTC","RET+ MTC/DTC","MTC","MTC","RET+ MTC"),
      n           = c(261, 417, 55, 331, 330, 122),
      `PFS HR`    = c(0.21, 0.59, NA, 0.46, 0.28, NA),
      `ORR (%)`   = c(65, 12, 69, 45, 28, 71),
      `PFS (mo)`  = c(18.3, 10.8, "NR", 30.5, 11.2, "NR")
    ) %>% datatable(rownames=FALSE, options=list(dom="t", pageLength=7))
  })

  output$pt_table <- renderTable({
    tibble(
      Parameter = c("Subtype","Baseline SLD","Baseline Tg","TSH","MAPK Activity","Scenario"),
      Value = c(input$subtype,
                paste(input$TumVol_init, "mm"),
                paste(input$Tg_init, "ng/mL"),
                "Suppressed (<0.1)",
                sprintf("%.1f", input$MAPK_level),
                input$scenario)
    )
  })

  ## Tab 2: PK
  output$pk_lenv <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~Cp_LENV_out, type="scatter", mode="lines",
            line=list(color="#2563EB",width=2)) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Cp (ng/mL)"))
  })
  output$pk_sora <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~Cp_SORA_out, type="scatter", mode="lines",
            line=list(color="#16A34A",width=2)) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Cp (ng/mL)"))
  })
  output$pk_selp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~Cp_SELP_out, type="scatter", mode="lines",
            line=list(color="#EA580C",width=2)) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Cp (ng/mL)"))
  })
  output$pk_param_table <- renderTable({
    tibble(
      Drug = c("Lenvatinib","Sorafenib","Selpercatinib"),
      `CL (L/h)` = c(4.0, 3.2, 6.5),
      `V1 (L)` = c(50, 55, 198),
      `t½ (h)` = c(round(0.693*50/4.0, 1), round(0.693*55/3.2, 1), round(0.693*198/6.5, 1)),
      F = c("85%","38%","73%")
    )
  })

  ## Tab 3: Pathways
  output$pathway_mapk <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~MAPK_act, type="scatter", mode="lines",
            line=list(color="#DC2626",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="MAPK Activity (0-1)", range=c(0,1)))
  })
  output$pathway_pi3k <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~PI3K_act, type="scatter", mode="lines",
            line=list(color="#F97316",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="PI3K/AKT Activity (0-1)", range=c(0,1)))
  })
  output$pathway_angio <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo) %>%
      add_lines(y=~VEGF, name="VEGF", line=list(color="#EF4444")) %>%
      add_lines(y=~Angio, name="Vasculature", line=list(color="#3B82F6")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Normalized"))
  })
  output$pathway_bar <- renderPlotly({
    d <- sim_data() %>% tail(7)
    bar_d <- tibble(
      Pathway = c("MAPK","PI3K","Vasculature"),
      Value   = c(mean(d$MAPK_act), mean(d$PI3K_act), mean(d$Angio))
    )
    plot_ly(bar_d, x=~Pathway, y=~Value, type="bar",
            marker=list(color=c("#DC2626","#F97316","#3B82F6"))) %>%
      layout(yaxis=list(title="Activity (0-1)", range=c(0,1.5)))
  })

  ## Tab 4: Tumor
  output$tumor_cells <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~TumorN, type="scatter", mode="lines",
            line=list(color="#DC2626",width=2)) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Tumor N (normalized)"))
  })
  output$tumor_vol <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~TumVol, type="scatter", mode="lines",
            line=list(color="#7C3AED",width=2)) %>%
      add_lines(y=rep(input$TumVol_init*0.7, nrow(d)), line=list(color="gray",dash="dash"),
                name="PR threshold") %>%
      add_lines(y=rep(input$TumVol_init*1.2, nrow(d)), line=list(color="red",dash="dash"),
                name="PD threshold") %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="SLD (mm)"))
  })
  output$tumor_growth <- renderPlotly({
    d <- sim_data() %>% mutate(net_gr = c(NA, diff(TumorN)/diff(time_mo)))
    plot_ly(d, x=~time_mo, y=~net_gr, type="scatter", mode="lines",
            line=list(color="#0891B2",width=2)) %>%
      add_lines(y=rep(0, nrow(d)), line=list(color="gray",dash="dash"), name="Zero growth") %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Net growth rate (/mo)"))
  })
  output$waterfall <- renderPlotly({
    d <- all_data() %>% group_by(scenario) %>% summarise(sld_pct=last(SLD_change))
    d <- d %>% arrange(sld_pct)
    colors <- ifelse(d$sld_pct < 0, "#22C55E", "#EF4444")
    plot_ly(d, x=~reorder(scenario,sld_pct), y=~sld_pct, type="bar",
            marker=list(color=colors)) %>%
      add_lines(y=c(-30,-30), x=c(0.5, nrow(d)+0.5), line=list(color="green",dash="dash")) %>%
      add_lines(y=c(20,20), x=c(0.5, nrow(d)+0.5), line=list(color="red",dash="dash")) %>%
      layout(xaxis=list(title=""), yaxis=list(title="% SLD Change at end"))
  })

  ## Tab 5: Biomarkers
  output$bm_tg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~Tg, type="scatter", mode="lines",
            line=list(color="#F59E0B",width=2)) %>%
      add_lines(y=rep(0.1, nrow(d)), line=list(color="green",dash="dash"), name="CR threshold") %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Tg (ng/mL)", type="log"))
  })
  output$bm_ct <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~CT, type="scatter", mode="lines",
            line=list(color="#EA580C",width=2)) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Calcitonin (pg/mL)", type="log"))
  })
  output$bm_tsh <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_mo, y=~TSH, type="scatter", mode="lines",
            line=list(color="#2563EB",width=2)) %>%
      add_lines(y=rep(0.1, nrow(d)), line=list(color="orange",dash="dash"), name="Target <0.1") %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="TSH (mIU/L)"))
  })
  output$bm_table <- renderTable({
    d <- sim_data() %>% tail(1)
    tibble(
      Biomarker = c("Tg", "Calcitonin", "TSH", "VEGF", "MAPK Activity"),
      Value = c(round(d$Tg,1), round(d$CT,1), round(d$TSH,2),
                round(d$VEGF,3), round(d$MAPK_act,3)),
      Unit = c("ng/mL","pg/mL","mIU/L","norm.","norm.")
    )
  })

  ## Tab 6: Clinical
  output$recist_plot <- renderPlotly({
    d <- sim_data()
    base_vol <- input$TumVol_init
    p <- plot_ly(d, x=~time_mo, y=~SLD_change, type="scatter", mode="lines",
                 line=list(color=pal[input$scenario], width=3)) %>%
      add_lines(y=rep(-30,nrow(d)), line=list(color="#16A34A",dash="dash"), name="PR (-30%)") %>%
      add_lines(y=rep(20,nrow(d)),  line=list(color="#DC2626",dash="dash"), name="PD (+20%)") %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="% Change in SLD", zeroline=TRUE))
    p
  })
  output$response_table <- renderTable({
    d <- sim_data()
    tibble(
      Category = c("CR (SLD→0)","PR (≤-30%)","SD","PD (≥+20%)"),
      Probability = c(
        paste0(round(mean(d$SLD_change <= -95)*100, 0),"%"),
        paste0(round(mean(d$SLD_change <= -30 & d$SLD_change > -95)*100, 0),"%"),
        paste0(round(mean(d$SLD_change > -30 & d$SLD_change < 20)*100, 0),"%"),
        paste0(round(mean(d$SLD_change >= 20)*100, 0),"%")
      )
    )
  })
  output$pfs_plot <- renderPlotly({
    d <- sim_data() %>% mutate(pd_flag = as.integer(SLD_change >= 20))
    pd_time <- d %>% filter(pd_flag==1) %>% slice(1) %>% pull(time_mo)
    pd_time <- ifelse(length(pd_time)==0, max(d$time_mo), pd_time)
    plot_ly(d, x=~time_mo, y=~(1-cummax(pd_flag)), type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(37,99,235,0.15)",
            line=list(color="#2563EB",width=2)) %>%
      add_segments(x=pd_time, xend=pd_time, y=0, yend=1,
                   line=list(color="red",dash="dash")) %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="PD-free Probability", range=c(0,1.05)))
  })
  output$bm_corr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~Tg, y=~TumVol, type="scatter", mode="markers",
            marker=list(color=~time_mo, colorscale="Viridis", size=5,
                        colorbar=list(title="Month"))) %>%
      layout(xaxis=list(title="Tg (ng/mL)", type="log"),
             yaxis=list(title="Tumor SLD (mm)"))
  })

  ## Tab 7: Comparison
  output$cmp_vol <- renderPlotly({
    d <- all_data(); gps <- split(d, d$scenario)
    p <- plot_ly()
    for(sc in names(pal)) {
      dd <- gps[[sc]]; if(is.null(dd)) next
      p <- add_lines(p, data=dd, x=~time_mo, y=~TumVol,
                     name=sc, line=list(color=pal[sc],width=2))
    }
    p %>% layout(xaxis=list(title="Time (months)"),
                 yaxis=list(title="SLD (mm)"),
                 legend=list(orientation="h",y=-0.25))
  })
  output$cmp_sld <- renderPlotly({
    d <- all_data(); gps <- split(d, d$scenario)
    p <- plot_ly()
    for(sc in names(pal)) {
      dd <- gps[[sc]]; if(is.null(dd)) next
      p <- add_lines(p, data=dd, x=~time_mo, y=~SLD_change,
                     name=sc, line=list(color=pal[sc],width=2))
    }
    p %>% add_lines(y=rep(-30,1), line=list(color="gray",dash="dash"), name="PR") %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="% SLD Change"),
             legend=list(orientation="h",y=-0.25))
  })
  output$cmp_tg <- renderPlotly({
    d <- all_data(); gps <- split(d, d$scenario)
    p <- plot_ly()
    for(sc in names(pal)) {
      dd <- gps[[sc]]; if(is.null(dd)) next
      p <- add_lines(p, data=dd, x=~time_mo, y=~Tg,
                     name=sc, line=list(color=pal[sc],width=2))
    }
    p %>% layout(xaxis=list(title="Time (months)"),
                 yaxis=list(title="Tg (ng/mL)", type="log"),
                 legend=list(orientation="h",y=-0.25))
  })
  output$cmp_table <- renderDT({
    d <- all_data() %>%
      group_by(scenario) %>%
      summarise(
        `Tg Final` = round(last(Tg),1),
        `Tumor Vol. Final (mm)` = round(last(TumVol),1),
        `SLD % Change` = round(last(SLD_change),1),
        `MAPK Mean` = round(mean(MAPK_act),3),
        .groups="drop"
      )
    datatable(d, rownames=FALSE, options=list(dom="t", pageLength=8)) %>%
      formatStyle("SLD % Change",
                  background=styleInterval(c(-30,20),
                                           c("#D1FAE5","#FEF9C3","#FEE2E2")))
  })
}

shinyApp(ui, server)
