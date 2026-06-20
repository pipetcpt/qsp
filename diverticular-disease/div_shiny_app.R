## =============================================================================
## Diverticular Disease (게실병) QSP Shiny Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & Risk Assessment (환자 프로파일)
##   2. Pharmacokinetics – Drug Concentration (PK)
##   3. PD / Biomarkers – Inflammation & Microbiome (PD 바이오마커)
##   4. Clinical Endpoints – Disease Progression (임상 결과)
##   5. Scenario Comparison (시나리오 비교)
##   6. Microbiome & Barrier Health (장내 미생물 & 점막)
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## Inline mrgsolve Model (minimal, fast)
## ─────────────────────────────────────────────────────────────────────────────

div_model_code <- '
$PARAM @annotated
Fiber_base  : 15.0 : Baseline dietary fiber (g/day)
Fiber_supp  : 0.0  : Fiber supplement (g/day)
NSAID_flag  : 0.0  : NSAID use (0/1)
kPress_base : 0.05 : Pressure build-up rate
kPress_decay: 0.10 : Pressure decay rate
Press_ss    : 45.0 : Pressure steady-state
kMuc_repair : 0.15 : Mucosal repair rate
kMuc_damage : 0.008: Mucosal damage rate
kNSAID_dmg  : 0.05 : NSAID mucosal damage
kProt_grow  : 0.20 : Protective bacteria growth
kProt_die   : 0.08 : Protective bacteria death
kPath_grow  : 0.12 : Pathogenic bacteria growth
kPath_die   : 0.05 : Pathogenic bacteria death
kCarry      : 1.5  : Microbiome capacity
kLPS_prod   : 0.30 : LPS production rate
kLPS_clear  : 0.50 : LPS clearance rate
kNFkB_act   : 0.80 : NF-kB activation
kNFkB_decay : 0.40 : NF-kB decay
kTNF_prod   : 1.20 : TNF production
kTNF_decay  : 0.35 : TNF decay
kIL6_prod   : 0.90 : IL-6 production
kIL6_decay  : 0.30 : IL-6 decay
kIL1b_prod  : 0.70 : IL-1b production
kIL1b_decay : 0.25 : IL-1b decay
kNeut_recr  : 0.60 : Neutrophil recruitment
kNeut_decay : 0.20 : Neutrophil decay
kCRP_prod   : 0.50 : CRP production
kCRP_decay  : 0.15 : CRP decay
kDivert_form: 0.002: Diverticula formation rate
kCollagen_loss: 0.003: Collagen degradation
kCollagen_repair: 0.05: Collagen repair
kChron_accum: 0.015: Chronic infl. accumulation
kChron_resol: 0.008: Chronic infl. resolution
kViscHyp_dev: 0.010: Visceral hypersens. dev.
kViscHyp_res: 0.005: Visceral hypersens. resol.
kRifax_elim : 0.69 : Rifaximin GI elimination
kMesa_elim  : 0.35 : Mesalamine GI elimination
kCipro_elim : 0.20 : Ciprofloxacin elimination
kMetro_elim : 0.12 : Metronidazole elimination
EC50_Rifax  : 8.0  : Rifaximin EC50
Emax_Rifax  : 0.65 : Rifaximin Emax
EC50_Mesa   : 12.0 : Mesalamine EC50
Emax_Mesa   : 0.60 : Mesalamine Emax
EC50_Cipro  : 1.5  : Ciprofloxacin EC50
Emax_Cipro  : 0.80 : Ciprofloxacin Emax
EC50_Metro  : 4.0  : Metronidazole EC50
Emax_Metro  : 0.75 : Metronidazole Emax
Emax_Fiber  : 0.70 : Fiber prebiotic Emax
EC50_Fiber  : 10.0 : Fiber prebiotic EC50
Acute_flag  : 0.0  : Acute diverticulitis trigger
Acute_onset : 30.0 : Acute onset day

$INIT
Fiber    : 15.0
Press    : 45.0
Mucosal  : 0.75
ProtBact : 0.60
PathBact : 0.40
LPS      : 0.25
NFkB     : 8.0
TNF      : 5.0
IL6      : 4.0
IL1b     : 3.0
Neut     : 5.0
CRP      : 3.0
DivertN  : 2.0
ChronInfl: 15.0
Rifax    : 0.0
Mesa     : 0.0
Cipro    : 0.0
Metro    : 0.0
ViscHyp  : 10.0
Collagen : 0.80

$MAIN
double acute_active = (Acute_flag > 0.5 && TIME > Acute_onset) ? 1.0 : 0.0;
double eff_Rifax = (Rifax > 0) ? Emax_Rifax * Rifax / (EC50_Rifax + Rifax) : 0.0;
double eff_Mesa  = (Mesa  > 0) ? Emax_Mesa  * Mesa  / (EC50_Mesa  + Mesa)  : 0.0;
double eff_Cipro = (Cipro > 0) ? Emax_Cipro * Cipro / (EC50_Cipro + Cipro) : 0.0;
double eff_Metro = (Metro > 0) ? Emax_Metro * Metro / (EC50_Metro + Metro) : 0.0;
double eff_antibio = 1.0 - (1.0 - eff_Cipro)*(1.0 - eff_Metro)*(1.0 - eff_Rifax);
double eff_Fiber = Emax_Fiber * Fiber / (EC50_Fiber + Fiber);
double butyrate_support = 0.5 * ProtBact;
double muc_damage = kMuc_damage * LPS * (1.0 + NSAID_flag * kNSAID_dmg / kMuc_damage);
double col_damage = kCollagen_loss * (NFkB / 50.0) * (1.0 - Collagen * 0.3);

$ODE
dxdt_Fiber    = -0.05 * Fiber + 0.05 * (Fiber_base + Fiber_supp);
dxdt_Press    = kPress_base * (Press_ss - Press) * (1.0 - Fiber/40.0)
                - kPress_decay * (Press - 20.0) * (Fiber/20.0)
                + 5.0 * acute_active;
dxdt_Mucosal  = kMuc_repair * (1.0 - Mucosal) * butyrate_support
                - muc_damage * Mucosal
                - Mucosal * acute_active * 0.1;
dxdt_ProtBact = kProt_grow * ProtBact * (1.0 - ProtBact/kCarry) * (1.0 + eff_Fiber)
                - kProt_die * ProtBact - 0.3 * eff_antibio * ProtBact;
dxdt_PathBact = kPath_grow * PathBact * (1.0 - PathBact/kCarry)
                * (1.0 - 0.5 * ProtBact/kCarry) * (1.0 + 2.0 * acute_active)
                - kPath_die * PathBact - eff_antibio * PathBact;
dxdt_LPS      = kLPS_prod * PathBact * (1.0 + acute_active * 3.0) / (Mucosal + 0.1)
                - kLPS_clear * LPS;
dxdt_NFkB     = kNFkB_act * LPS * (1.0 - eff_Mesa)
                + 0.05 * TNF + 0.03 * IL1b - kNFkB_decay * NFkB;
dxdt_TNF      = kTNF_prod * (NFkB/100.0) * (1.0 - eff_Mesa) - kTNF_decay * TNF;
dxdt_IL6      = kIL6_prod * (NFkB/100.0) * (1.0 - eff_Mesa) - kIL6_decay * IL6;
dxdt_IL1b     = kIL1b_prod * (NFkB/100.0) * (1.0 - eff_Mesa) - kIL1b_decay * IL1b;
dxdt_Neut     = kNeut_recr * (TNF + IL1b) / 20.0 - kNeut_decay * Neut;
dxdt_CRP      = kCRP_prod * IL6 * (1.0 - eff_Mesa * 0.4) - kCRP_decay * CRP;
dxdt_DivertN  = kDivert_form * (Press - 25.0) * (1.0 - Collagen) * (Press > 25.0 ? 1.0 : 0.0);
dxdt_ChronInfl= kChron_accum * (NFkB/100.0) * ChronInfl * (1.0 - ChronInfl/100.0)
                - kChron_resol * ChronInfl * (1.0 + eff_Mesa)
                + kChron_accum * 2.0 * acute_active;
dxdt_Rifax    = -kRifax_elim * Rifax;
dxdt_Mesa     = -kMesa_elim  * Mesa;
dxdt_Cipro    = -kCipro_elim * Cipro;
dxdt_Metro    = -kMetro_elim * Metro;
dxdt_ViscHyp  = kViscHyp_dev * ChronInfl/100.0 * (1.0 - ViscHyp/100.0)
                - kViscHyp_res * ViscHyp;
dxdt_Collagen = kCollagen_repair * (1.0 - Collagen) * butyrate_support - col_damage;

$TABLE
capture Pain_score = fmin(10.0, fmax(0.0, 0.5 + 0.08*NFkB + 0.05*ViscHyp));
capture Microbiome_index = ProtBact / (ProtBact + PathBact + 0.001);
capture Complication_risk = fmin(1.0, 0.01*DivertN + 0.005*ChronInfl + 0.02*(1.0-Mucosal));
capture eff_Rifax eff_Mesa eff_Cipro eff_Metro eff_antibio eff_Fiber
capture acute_active butyrate_support
'

# Compile model
div_mod <- mcode("div_shiny_qsp", div_model_code, quiet = TRUE)

# Helper: run simulation
run_sim <- function(mod, fiber_supp, rifax_flag, mesa_flag,
                    cipro_flag, metro_flag, acute_flag, nsaid_flag,
                    end_time = 365, delta = 1, acute_onset = 30) {

  events_list <- list()

  if (rifax_flag) {
    for (m in 0:11) {
      for (d in seq(0, 6, by = 1/3)) {
        t_dose <- m * 30 + d
        if (t_dose <= end_time) {
          events_list[[length(events_list) + 1]] <-
            data.frame(time = t_dose, amt = 400, cmt = "Rifax", evid = 1)
        }
      }
    }
  }
  if (mesa_flag) {
    for (d in 0:end_time) {
      events_list[[length(events_list) + 1]] <-
        data.frame(time = d, amt = 1600, cmt = "Mesa", evid = 1)
    }
  }
  if (cipro_flag) {
    for (i in 0:19) {
      t_dose <- acute_onset + i * 0.5
      if (t_dose <= end_time)
        events_list[[length(events_list) + 1]] <-
          data.frame(time = t_dose, amt = 500, cmt = "Cipro", evid = 1)
    }
  }
  if (metro_flag) {
    for (i in 0:29) {
      t_dose <- acute_onset + i * 0.333
      if (t_dose <= end_time)
        events_list[[length(events_list) + 1]] <-
          data.frame(time = t_dose, amt = 500, cmt = "Metro", evid = 1)
    }
  }

  events_df <- if (length(events_list) > 0)
    bind_rows(events_list) %>% arrange(time)
  else
    data.frame(time = 0, amt = 0, cmt = "Rifax", evid = 0)

  mod %>%
    param(Fiber_supp    = fiber_supp,
          Acute_flag    = ifelse(acute_flag, 1, 0),
          Acute_onset   = acute_onset,
          NSAID_flag    = ifelse(nsaid_flag, 1, 0)) %>%
    mrgsim(end = end_time, delta = delta,
           events = events_df) %>%
    as.data.frame()
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  titlePanel(
    div(
      h2("게실병 (Diverticular Disease) QSP Dashboard",
         style = "color:#2c3e50; font-weight:bold; margin:0"),
      h5("Quantitative Systems Pharmacology Model — Diverticulosis / Diverticulitis",
         style = "color:#7f8c8d; margin:0")
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # Patient Profile
      wellPanel(
        h4("Patient Profile", style = "color:#2980b9; font-weight:bold"),
        sliderInput("fiber_supp", "Fiber Supplement Add-on (g/day)",
                    min = 0, max = 20, value = 0, step = 1),
        sliderInput("age_factor", "Patient Age (years)",
                    min = 40, max = 80, value = 55, step = 5),
        checkboxInput("nsaid_flag", "Chronic NSAID Use", FALSE),
        sliderInput("end_time", "Simulation Duration (days)",
                    min = 90, max = 730, value = 365, step = 30)
      ),

      # Treatment Options
      wellPanel(
        h4("Treatment Options", style = "color:#27ae60; font-weight:bold"),
        checkboxInput("rifax_flag",  "Rifaximin 400mg TID × 7d/month", FALSE),
        checkboxInput("mesa_flag",   "Mesalamine (5-ASA) 1.6g/day",     FALSE),
        checkboxInput("cipro_flag",  "Ciprofloxacin 500mg BID",          FALSE),
        checkboxInput("metro_flag",  "Metronidazole 500mg TID",          FALSE),
        checkboxInput("acute_flag",  "Simulate Acute Diverticulitis",    FALSE),
        conditionalPanel(
          condition = "input.acute_flag == true",
          sliderInput("acute_onset", "Acute Event Day",
                      min = 7, max = 180, value = 30, step = 1)
        )
      ),

      # Comparison Scenario
      wellPanel(
        h4("Comparison Scenario", style = "color:#9b59b6; font-weight:bold"),
        selectInput("comp_scen", "Compare Against:",
                    choices = c(
                      "Natural History (no Tx)"   = "natural",
                      "High Fiber (+15g/day)"      = "fiber",
                      "Rifaximin Cyclic"           = "rifax",
                      "Mesalamine 1.6g/day"        = "mesa",
                      "Combination (Fiber+Rifax+Mesa)" = "combo"
                    ), selected = "natural")
      ),

      actionButton("run_sim", "Run Simulation",
                   class = "btn-primary btn-block",
                   style = "margin-top:10px"),

      hr(),
      tags$small(
        "Model: 20-ODE QSP | Reference: Tursi A et al. (2020) Nat Rev Dis Primers 6:20"
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # ── Tab 1: Patient Profile & Risk ──
        tabPanel(
          "Patient Profile",
          icon = icon("user-md"),
          br(),
          fluidRow(
            column(4, wellPanel(
              h4("Disease Stage Assessment"),
              tableOutput("stage_table")
            )),
            column(4, wellPanel(
              h4("Risk Factor Summary"),
              tableOutput("risk_table")
            )),
            column(4, wellPanel(
              h4("Complication Risk"),
              plotOutput("gauge_plot", height = "220px")
            ))
          ),
          fluidRow(
            column(12, plotlyOutput("profile_plot", height = "350px"))
          )
        ),

        # ── Tab 2: PK ──
        tabPanel(
          "Pharmacokinetics (PK)",
          icon = icon("pills"),
          br(),
          fluidRow(
            column(12,
              h4("Drug Concentration Profiles", style = "color:#2980b9"),
              plotlyOutput("pk_plot", height = "450px")
            )
          ),
          fluidRow(
            column(12,
              h4("PK Summary Table"),
              DTOutput("pk_table")
            )
          )
        ),

        # ── Tab 3: PD Biomarkers ──
        tabPanel(
          "PD / Biomarkers",
          icon = icon("vial"),
          br(),
          fluidRow(
            column(6, plotlyOutput("inflam_plot", height = "350px")),
            column(6, plotlyOutput("cytokine_plot", height = "350px"))
          ),
          fluidRow(
            column(6, plotlyOutput("crp_neut_plot", height = "300px")),
            column(6, plotlyOutput("collagen_plot", height = "300px"))
          )
        ),

        # ── Tab 4: Clinical Endpoints ──
        tabPanel(
          "Clinical Endpoints",
          icon = icon("hospital"),
          br(),
          fluidRow(
            column(6, plotlyOutput("divert_plot", height = "320px")),
            column(6, plotlyOutput("pain_plot",   height = "320px"))
          ),
          fluidRow(
            column(6, plotlyOutput("chron_plot",  height = "320px")),
            column(6, plotlyOutput("risk_plot",   height = "320px"))
          ),
          br(),
          DTOutput("endpoint_table")
        ),

        # ── Tab 5: Scenario Comparison ──
        tabPanel(
          "Scenario Comparison",
          icon = icon("chart-bar"),
          br(),
          fluidRow(
            column(12,
              h4("Selected Patient vs. Comparison Scenario"),
              plotlyOutput("compare_plot1", height = "400px")
            )
          ),
          fluidRow(
            column(6, plotlyOutput("compare_plot2", height = "350px")),
            column(6, plotlyOutput("compare_plot3", height = "350px"))
          ),
          br(),
          DTOutput("compare_table")
        ),

        # ── Tab 6: Microbiome & Barrier ──
        tabPanel(
          "Microbiome & Gut Barrier",
          icon = icon("bacteria"),
          br(),
          fluidRow(
            column(6, plotlyOutput("micro_plot",   height = "350px")),
            column(6, plotlyOutput("barrier_plot", height = "350px"))
          ),
          fluidRow(
            column(6, plotlyOutput("lps_plot",    height = "300px")),
            column(6, plotlyOutput("visceral_plot",height = "300px"))
          ),
          br(),
          h4("Microbiome Health Summary", style = "color:#27ae60"),
          DTOutput("micro_table")
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive: run patient simulation
  patient_sim <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", {
      run_sim(
        mod        = div_mod,
        fiber_supp = input$fiber_supp,
        rifax_flag = input$rifax_flag,
        mesa_flag  = input$mesa_flag,
        cipro_flag = input$cipro_flag,
        metro_flag = input$metro_flag,
        acute_flag = input$acute_flag,
        nsaid_flag = input$nsaid_flag,
        end_time   = input$end_time,
        acute_onset= if (input$acute_flag) input$acute_onset else 9999
      )
    })
  }, ignoreNULL = FALSE)

  # Reactive: comparison scenario
  comp_sim <- eventReactive(input$run_sim, {
    params <- switch(input$comp_scen,
      "natural" = list(fs = 0,  rf = FALSE, ms = FALSE, cf = FALSE, mf = FALSE),
      "fiber"   = list(fs = 15, rf = FALSE, ms = FALSE, cf = FALSE, mf = FALSE),
      "rifax"   = list(fs = 0,  rf = TRUE,  ms = FALSE, cf = FALSE, mf = FALSE),
      "mesa"    = list(fs = 0,  rf = FALSE, ms = TRUE,  cf = FALSE, mf = FALSE),
      "combo"   = list(fs = 15, rf = TRUE,  ms = TRUE,  cf = FALSE, mf = FALSE)
    )
    run_sim(
      mod        = div_mod,
      fiber_supp = params$fs,
      rifax_flag = params$rf,
      mesa_flag  = params$ms,
      cipro_flag = params$cf,
      metro_flag = params$mf,
      acute_flag = FALSE,
      nsaid_flag = FALSE,
      end_time   = input$end_time,
      acute_onset= 9999
    )
  }, ignoreNULL = FALSE)

  # ── Tab 1: Patient Profile ──

  output$stage_table <- renderTable({
    df <- patient_sim()
    last <- tail(df, 1)
    data.frame(
      Parameter  = c("CRP (mg/L)", "Diverticula (N)", "Chronic Infl.", "Pain Score",
                     "Visceral Hypersens.", "Collagen Index"),
      Value      = round(c(last$CRP, last$DivertN, last$ChronInfl,
                           last$Pain_score, last$ViscHyp, last$Collagen), 2),
      Status     = c(
        ifelse(last$CRP < 5, "Normal", ifelse(last$CRP < 50, "Elevated", "High")),
        ifelse(last$DivertN < 5, "Mild", ifelse(last$DivertN < 15, "Moderate", "Severe")),
        ifelse(last$ChronInfl < 20, "Mild", ifelse(last$ChronInfl < 50, "Moderate", "Severe")),
        ifelse(last$Pain_score < 3, "Mild", ifelse(last$Pain_score < 6, "Moderate", "Severe")),
        ifelse(last$ViscHyp < 20, "Normal", ifelse(last$ViscHyp < 50, "Moderate", "High")),
        ifelse(last$Collagen > 0.7, "Good", ifelse(last$Collagen > 0.4, "Fair", "Poor"))
      )
    )
  })

  output$risk_table <- renderTable({
    data.frame(
      Risk_Factor = c("Dietary fiber", "NSAID use", "Acute event", "Age"),
      Status      = c(
        ifelse(input$fiber_supp >= 10, "Favorable", "Risk"),
        ifelse(input$nsaid_flag, "Risk", "None"),
        ifelse(input$acute_flag, "Active", "None"),
        ifelse(input$age_factor >= 65, "High risk", "Moderate")
      )
    )
  })

  output$gauge_plot <- renderPlot({
    df   <- patient_sim()
    last <- tail(df, 1)
    risk <- round(last$Complication_risk * 100, 1)
    col  <- ifelse(risk < 20, "#27ae60", ifelse(risk < 50, "#f39c12", "#e74c3c"))
    par(mar = c(0,0,2,0))
    pie(c(risk, 100 - risk),
        col   = c(col, "#ecf0f1"),
        border= "white",
        labels= NA,
        main  = paste0("Complication Risk\n", risk, "%"),
        cex.main = 1.1)
  })

  output$profile_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~CRP, name = "CRP (mg/L)", line = list(color="#e74c3c")) %>%
      add_lines(y = ~NFkB, name = "NF-κB (AU)", line = list(color="#3498db")) %>%
      add_lines(y = ~Pain_score * 5, name = "Pain Score × 5", line = list(color="#9b59b6")) %>%
      add_lines(y = ~Collagen * 100, name = "Collagen × 100", line = list(color="#27ae60")) %>%
      layout(title  = "Overall Disease Profile Over Time",
             xaxis  = list(title = "Day"),
             yaxis  = list(title = "Value"),
             legend = list(orientation = "h"))
  })

  # ── Tab 2: PK ──

  output$pk_plot <- renderPlotly({
    df <- patient_sim()
    p  <- plot_ly(df, x = ~time)
    if (any(df$Rifax > 0.01))
      p <- add_lines(p, y = ~Rifax, name = "Rifaximin GI (mg/L)",
                     line = list(color = "#16a085"))
    if (any(df$Mesa > 0.01))
      p <- add_lines(p, y = ~Mesa,  name = "Mesalamine GI (mg/L)",
                     line = list(color = "#8e44ad"))
    if (any(df$Cipro > 0.01))
      p <- add_lines(p, y = ~Cipro, name = "Ciprofloxacin Plasma (mg/L)",
                     line = list(color = "#e67e22"))
    if (any(df$Metro > 0.01))
      p <- add_lines(p, y = ~Metro, name = "Metronidazole Plasma (mg/L)",
                     line = list(color = "#c0392b"))
    p %>% layout(title = "Drug Concentration Profiles",
                 xaxis = list(title = "Day"),
                 yaxis = list(title = "Concentration (mg/L)"),
                 legend = list(orientation = "h"))
  })

  output$pk_table <- renderDT({
    df <- patient_sim()
    data.frame(
      Drug        = c("Rifaximin (GI)", "Mesalamine (GI)",
                      "Ciprofloxacin (plasma)", "Metronidazole (plasma)"),
      Max_Conc    = round(c(max(df$Rifax), max(df$Mesa),
                             max(df$Cipro), max(df$Metro)), 3),
      Mean_Conc   = round(c(mean(df$Rifax), mean(df$Mesa),
                             mean(df$Cipro), mean(df$Metro)), 3),
      Active_days = c(sum(df$Rifax > 0.1), sum(df$Mesa > 0.1),
                       sum(df$Cipro > 0.1), sum(df$Metro > 0.1)),
      Eff_Mean    = round(c(mean(df$eff_Rifax), mean(df$eff_Mesa),
                             mean(df$eff_Cipro), mean(df$eff_Metro)), 3)
    )
  }, options = list(pageLength = 5, dom = "t"))

  # ── Tab 3: PD / Biomarkers ──

  output$inflam_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~NFkB, name = "NF-κB (AU)", line = list(color = "#c0392b")) %>%
      add_lines(y = ~TNF,  name = "TNF-α (pg/mL)", line = list(color = "#e74c3c")) %>%
      layout(title = "NF-κB & TNF-α Dynamics",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Value"))
  })

  output$cytokine_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~IL6,  name = "IL-6 (pg/mL)",  line = list(color = "#3498db")) %>%
      add_lines(y = ~IL1b, name = "IL-1β (pg/mL)", line = list(color = "#9b59b6")) %>%
      layout(title = "Cytokine Profile (IL-6, IL-1β)",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Concentration (pg/mL)"))
  })

  output$crp_neut_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~CRP,  name = "CRP (mg/L)",    line = list(color = "#e74c3c")) %>%
      add_lines(y = ~Neut, name = "Neutrophil Score", line = list(color = "#e67e22")) %>%
      add_hline(y = 5, line = list(dash = "dot", color = "gray"),
                annotation = list(text = "CRP normal", x = 10)) %>%
      layout(title = "CRP & Neutrophil Score",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Value"))
  })

  output$collagen_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~Collagen, name = "Collagen Index (0-1)",
                line = list(color = "#27ae60")) %>%
      add_lines(y = ~Mucosal, name = "Mucosal Integrity (0-1)",
                line = list(color = "#1abc9c")) %>%
      layout(title = "Structural Integrity Indices",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Index (0-1)", range = c(0, 1.2)))
  })

  # ── Tab 4: Clinical Endpoints ──

  output$divert_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time, y = ~DivertN, type = "scatter", mode = "lines",
            line = list(color = "#e67e22", width = 2.5),
            name = "Diverticula Count") %>%
      layout(title = "Cumulative Diverticula Formation",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Diverticula (N)"))
  })

  output$pain_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~Pain_score, name = "Pain Score (VAS 0-10)",
                line = list(color = "#e74c3c")) %>%
      add_lines(y = ~ViscHyp / 10, name = "Visceral Hypersens. / 10",
                line = list(color = "#9b59b6", dash = "dash")) %>%
      layout(title = "Pain Score & Visceral Hypersensitivity",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score"))
  })

  output$chron_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time, y = ~ChronInfl, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(52,152,219,0.2)",
            line = list(color = "#2980b9"),
            name = "Chronic Infl.") %>%
      layout(title = "Chronic Inflammation Score (SUDD)",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score (0-100)"))
  })

  output$risk_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time, y = ~Complication_risk * 100, type = "scatter",
            mode = "lines", fill = "tozeroy",
            fillcolor = "rgba(231,76,60,0.15)",
            line = list(color = "#e74c3c"),
            name = "Complication Risk (%)") %>%
      layout(title = "Complication Risk Over Time",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Risk (%)"))
  })

  output$endpoint_table <- renderDT({
    df <- patient_sim()
    selected_days <- c(30, 90, 180, 365)
    selected_days <- selected_days[selected_days <= max(df$time)]
    df %>%
      filter(time %in% selected_days) %>%
      select(Day = time, CRP, Diverticula = DivertN,
             ChronInfl, Pain_score, ViscHyp, Collagen, Mucosal,
             Microbiome_index, Complication_risk) %>%
      mutate(across(where(is.numeric), ~round(.x, 3)))
  }, options = list(pageLength = 5, scrollX = TRUE))

  # ── Tab 5: Scenario Comparison ──

  output$compare_plot1 <- renderPlotly({
    df1 <- patient_sim() %>% mutate(Group = "Patient Scenario")
    df2 <- comp_sim()   %>% mutate(Group = input$comp_scen)
    combined <- bind_rows(df1, df2)
    plot_ly(combined, x = ~time, y = ~CRP, color = ~Group,
            type = "scatter", mode = "lines") %>%
      layout(title = "CRP Comparison: Patient vs. Reference Scenario",
             xaxis = list(title = "Day"),
             yaxis = list(title = "CRP (mg/L)"),
             legend = list(orientation = "h"))
  })

  output$compare_plot2 <- renderPlotly({
    df1 <- patient_sim() %>% mutate(Group = "Patient")
    df2 <- comp_sim()   %>% mutate(Group = "Reference")
    combined <- bind_rows(df1, df2)
    plot_ly(combined, x = ~time, y = ~DivertN, color = ~Group,
            type = "scatter", mode = "lines") %>%
      layout(title = "Diverticula Count Comparison",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Diverticula (N)"),
             legend = list(orientation = "h"))
  })

  output$compare_plot3 <- renderPlotly({
    df1 <- patient_sim() %>% mutate(Group = "Patient")
    df2 <- comp_sim()   %>% mutate(Group = "Reference")
    combined <- bind_rows(df1, df2)
    plot_ly(combined, x = ~time, y = ~ChronInfl, color = ~Group,
            type = "scatter", mode = "lines") %>%
      layout(title = "Chronic Inflammation Score Comparison",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score (0-100)"),
             legend = list(orientation = "h"))
  })

  output$compare_table <- renderDT({
    df1 <- patient_sim() %>%
      filter(time == max(time)) %>%
      mutate(Scenario = "Patient") %>%
      select(Scenario, CRP, DivertN, ChronInfl, Pain_score, Microbiome_index,
             Collagen, Complication_risk)
    df2 <- comp_sim() %>%
      filter(time == max(time)) %>%
      mutate(Scenario = "Reference") %>%
      select(Scenario, CRP, DivertN, ChronInfl, Pain_score, Microbiome_index,
             Collagen, Complication_risk)
    bind_rows(df1, df2) %>%
      mutate(across(where(is.numeric), ~round(.x, 3)))
  }, options = list(pageLength = 5, dom = "t"))

  # ── Tab 6: Microbiome & Barrier ──

  output$micro_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~ProtBact, name = "Protective Bacteria",
                line = list(color = "#27ae60")) %>%
      add_lines(y = ~PathBact, name = "Pathogenic Bacteria",
                line = list(color = "#e74c3c")) %>%
      add_lines(y = ~Microbiome_index, name = "Microbiome Health Index",
                line = list(color = "#3498db", dash = "dash")) %>%
      layout(title = "Gut Microbiome Dynamics",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Normalized Units"),
             legend = list(orientation = "h"))
  })

  output$barrier_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~Mucosal, name = "Mucosal Integrity (0-1)",
                line = list(color = "#8e44ad")) %>%
      add_lines(y = ~Collagen, name = "Collagen Index (0-1)",
                line = list(color = "#e67e22")) %>%
      add_lines(y = ~butyrate_support, name = "Butyrate Support",
                line = list(color = "#27ae60", dash = "dot")) %>%
      layout(title = "Gut Barrier & Structural Integrity",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Index (0-1)", range = c(0, 1.2)))
  })

  output$lps_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~LPS, name = "LPS (ng/mL)",
                line = list(color = "#c0392b")) %>%
      layout(title = "Circulating LPS / Endotoxin",
             xaxis = list(title = "Day"),
             yaxis = list(title = "LPS (ng/mL)"))
  })

  output$visceral_plot <- renderPlotly({
    df <- patient_sim()
    plot_ly(df, x = ~time) %>%
      add_lines(y = ~ViscHyp, name = "Visceral Hypersensitivity Score",
                line = list(color = "#9b59b6")) %>%
      add_lines(y = ~Pain_score * 10, name = "Pain Score × 10",
                line = list(color = "#e74c3c", dash = "dash")) %>%
      layout(title = "Visceral Hypersensitivity & Pain",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score"))
  })

  output$micro_table <- renderDT({
    df <- patient_sim()
    selected_days <- c(30, 90, 180, 365)
    selected_days <- selected_days[selected_days <= max(df$time)]
    df %>%
      filter(time %in% selected_days) %>%
      select(Day = time, ProtBact, PathBact, LPS, Mucosal, Collagen,
             Microbiome_index, butyrate_support) %>%
      mutate(across(where(is.numeric), ~round(.x, 4)))
  }, options = list(pageLength = 5, dom = "t"))
}

## ─────────────────────────────────────────────────────────────────────────────
## Launch App
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
