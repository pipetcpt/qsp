##############################################################################
# Pernicious Anemia (악성 빈혈) — QSP Interactive Shiny Dashboard
#
# Tabs:
#  1. 환자 프로파일   — Patient profile, risk factors, PA stage
#  2. 약물 PK         — Cobalamin PK (plasma B12, HoloTC, depot)
#  3. 혈액 PD         — Hematological PD (Hb, MCV, retic crisis)
#  4. 신경계 PD       — Neurological PD (NEURO score, SCD timeline)
#  5. 바이오마커       — Functional markers (MMA, Hcy, gastrin)
#  6. 시나리오 비교    — Multi-scenario comparison (5 treatment arms)
#
# Dependencies: shiny, mrgsolve, ggplot2, dplyr, tidyr, plotly, DT
##############################################################################

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ============================================================
# MODEL CODE (inline)
# ============================================================
pa_model_code <- '
$PARAM
KA_IM=0.40 KA_ORAL=0.08 F_PASSIVE=0.01 KIF_ABS=2.50 KIF_MAX=1.50
CL_PLASMA=28.0 V_PLASMA=6.0 K12_LIVER=0.30 K21_LIVER=0.015
K12_BM=0.08 K21_BM=0.025 K12_NERVE=0.04 K21_NERVE=0.002
KREN=0.25 ENTERO_F=0.005 KAB_PROD=0.0005 KAB_DECAY=0.0002
AB0=0.80 KIF_DECAY=0.003 KIF_REGEN=0.0001
HGB0=7.5 HGB_MAX=14.5 EC50_HGB=300.0 EMAX_HGB=1.0 KOUT_HGB=0.004
MCV0=115.0 MCV_NL=88.0 K_MCV=0.005 EC50_MCV=200.0
RETIC_BASE=30.0 RETIC_MAX=400.0 RETIC_EC50=250.0 K_RETIC=0.05
NEURO0=4.0 NEURO_MAX=8.0 NEURO_MIN=1.0 KNEURO_PROG=0.0008 KNEURO_RECOV=0.003
EC50_NEURO=400.0 MMA_BASE=3.5 MMA_NORM=0.15 K_MMA=0.05
HCY_BASE=28.0 HCY_NORM=9.0 K_HCY=0.04

$CMT DEPOT ORAL_GI IF_POOL CBA_PORT PLASMA HOLOTC LIVER BONE_MRW NERVE AUTO_AB PARIETAL HGB MCV RETIC NEURO

$INIT
DEPOT=0 ORAL_GI=0 IF_POOL=0.20 CBA_PORT=0 PLASMA=480 HOLOTC=12 LIVER=1000
BONE_MRW=0.8 NERVE=0.5 AUTO_AB=0.80 PARIETAL=0.20 HGB=7.5 MCV=115 RETIC=30 NEURO=4.0

$ODE
double Cp=PLASMA/V_PLASMA;
double IF_eff=IF_POOL*(1.0-AUTO_AB*0.9);
double IF_abs_rate=KIF_ABS*IF_eff*KIF_MAX;
double IF_med=(ORAL_GI>0)?IF_abs_rate*ORAL_GI/(ORAL_GI+0.1):0.0;
double passive=KA_ORAL*F_PASSIVE*ORAL_GI;
dxdt_DEPOT=-KA_IM*DEPOT;
dxdt_ORAL_GI=-KA_ORAL*ORAL_GI;
dxdt_IF_POOL=KIF_REGEN*(PARIETAL-IF_POOL)-KIF_DECAY*AUTO_AB*IF_POOL;
double entero=ENTERO_F*LIVER;
dxdt_CBA_PORT=IF_med+passive+entero-0.5*CBA_PORT;
double lu=K12_LIVER*PLASMA; double lr=K21_LIVER*LIVER;
double bu=K12_BM*PLASMA;    double br=K21_BM*BONE_MRW;
double nu=K12_NERVE*PLASMA;  double nr=K21_NERVE*NERVE;
double ren=(Cp>300.0)?KREN*(Cp-300.0)*V_PLASMA:0.0;
dxdt_PLASMA=KA_IM*DEPOT+0.5*CBA_PORT+lr+br+nr-lu-bu-nu-ren;
dxdt_HOLOTC=0.10*(0.20*Cp*0.738-HOLOTC);
dxdt_LIVER=lu-lr-entero;
dxdt_BONE_MRW=bu-br;
dxdt_NERVE=nu-nr;
dxdt_AUTO_AB=KAB_PROD*(1.0-AUTO_AB)-KAB_DECAY*AUTO_AB;
dxdt_PARIETAL=-KIF_DECAY*AUTO_AB*PARIETAL+KIF_REGEN*(0.05-PARIETAL);
double HGB_target=HGB0+(HGB_MAX-HGB0)*(Cp/(Cp+EC50_HGB))*EMAX_HGB;
dxdt_HGB=KOUT_HGB*(HGB_target-HGB);
double MCV_target=MCV_NL+(MCV0-MCV_NL)*(1.0-Cp/(Cp+EC50_MCV));
dxdt_MCV=K_MCV*(MCV_target-MCV);
double RETIC_target=RETIC_BASE+(RETIC_MAX-RETIC_BASE)*(Cp/(Cp+RETIC_EC50));
dxdt_RETIC=K_RETIC*(RETIC_target-RETIC);
double np=Cp/(Cp+EC50_NEURO);
dxdt_NEURO=KNEURO_PROG*(NEURO_MAX-NEURO)*(1.0-np)-KNEURO_RECOV*(NEURO-NEURO_MIN)*np;

$TABLE
double Cp_pgmL=PLASMA/V_PLASMA;
double MMA=MMA_NORM+(MMA_BASE-MMA_NORM)*exp(-K_MMA*(Cp_pgmL/100.0));
double HCY=HCY_NORM+(HCY_BASE-HCY_NORM)*exp(-K_HCY*(Cp_pgmL/200.0));
capture Cp_pgmL MMA HCY
$CAPTURE Cp_pgmL HGB MCV RETIC NEURO HOLOTC MMA HCY PARIETAL AUTO_AB LIVER BONE_MRW NERVE
'

pa_mod <- mrgsolve::mcode("PA_Shiny", pa_model_code)

# ============================================================
# HELPER: build events from inputs
# ============================================================
build_events <- function(route, dose_mg, freq_h, n_doses, load_dose, load_n, load_freq_h) {
  ev_list <- list()
  # Loading events
  if (load_n > 0) {
    load_times <- seq(0, (load_n-1)*load_freq_h, by=load_freq_h)
    ev_list <- c(ev_list, list(data.frame(
      time=load_times, amt=load_dose*1000, cmt=ifelse(route=="IM", 1, 2), evid=1
    )))
  }
  # Maintenance events
  start_main <- if (load_n > 0) (load_n-1)*load_freq_h + freq_h else 0
  main_times  <- seq(start_main, start_main + (n_doses-1)*freq_h, by=freq_h)
  ev_list <- c(ev_list, list(data.frame(
    time=main_times, amt=dose_mg*1000, cmt=ifelse(route=="IM", 1, 2), evid=1
  )))
  ev_df <- bind_rows(ev_list)
  ev_df <- ev_df[ev_df$time <= 365*24, ]
  ev_df
}

run_sim <- function(ev_df, tsim=365*24, delta=12, params=list()) {
  mod <- do.call(param, c(list(pa_mod), params))
  if (nrow(ev_df)==0) {
    out <- mod %>% mrgsim(end=tsim, delta=delta) %>% as.data.frame()
  } else {
    out <- mod %>% ev(ev_df) %>% mrgsim(end=tsim, delta=delta) %>% as.data.frame()
  }
  out$time_days <- out$time/24
  out
}

# Scenario presets
scenario_presets <- list(
  "S1: Untreated PA"  = list(route="IM", dose=0,    freq=720, n=0,  load=0,    load_n=0, load_freq=24),
  "S2: IM Standard"   = list(route="IM", dose=1000, freq=720, n=24, load=1000, load_n=7, load_freq=24),
  "S3: IM Maintenance"= list(route="IM", dose=1000, freq=720, n=24, load=0,    load_n=0, load_freq=24),
  "S4: Oral HD"       = list(route="Oral", dose=1000, freq=24, n=365, load=0,  load_n=0, load_freq=24),
  "S5: Aggressive IM" = list(route="IM", dose=2000, freq=720, n=24, load=2000, load_n=14,load_freq=24)
)

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Arial, sans-serif; background: #F8F9FA; }
      .well { background: #FFFFFF; border: 1px solid #DEE2E6; border-radius: 8px; }
      .nav-pills > li > a { border-radius: 6px; }
      .nav-pills > li.active > a { background-color: #2C3E50; }
      h4 { color: #2C3E50; font-weight: 700; }
      .badge-info { background-color: #3498DB; }
    "))
  ),

  titlePanel(
    div(
      h3("Pernicious Anemia QSP Dashboard", style="color:#2C3E50; font-weight:700; margin:0"),
      p("악성 빈혈 정량적 시스템 약리학 시뮬레이터",
        style="color:#7F8C8D; margin:0; font-size:14px")
    )
  ),
  hr(),

  sidebarLayout(
    sidebarPanel(
      width=3,
      wellPanel(
        h4("Patient Profile (환자 프로파일)"),
        sliderInput("baseline_hgb", "Baseline Hb (g/dL)", 4, 12, 7.5, 0.5),
        sliderInput("baseline_neuro", "Baseline Neuro Score (0-10)", 0, 8, 4, 0.5),
        sliderInput("baseline_mcv", "Baseline MCV (fL)", 100, 130, 115, 1),
        selectInput("pa_severity", "PA Severity",
                    choices=c("Mild (parietal 40%)"="mild",
                              "Moderate (parietal 20%)"="moderate",
                              "Severe (parietal 5%)"="severe"))
      ),
      wellPanel(
        h4("Treatment (치료 설정)"),
        selectInput("route", "Route", choices=c("IM Injection"="IM", "Oral"="Oral")),
        numericInput("dose_ug", "Maintenance Dose (µg)", 1000, 100, 5000, 100),
        selectInput("freq", "Maintenance Frequency",
                    choices=c("Daily (24h)"=24, "Weekly (168h)"=168,
                              "Biweekly (336h)"=336, "Monthly (720h)"=720)),
        hr(),
        h5("Loading Regimen"),
        numericInput("load_dose", "Loading Dose (µg)", 1000, 0, 5000, 500),
        numericInput("load_n", "# Loading Doses", 7, 0, 30, 1),
        selectInput("load_freq", "Loading Frequency",
                    choices=c("Daily (24h)"=24, "Every 12h"=12)),
        sliderInput("sim_months", "Simulation Duration (months)", 3, 24, 12, 1)
      ),
      wellPanel(
        h4("PK/PD Parameters"),
        sliderInput("f_passive", "Passive Oral F (%)", 0.1, 5, 1.0, 0.1),
        sliderInput("ec50_hgb", "EC50 Hgb (pg/mL)", 50, 800, 300, 25),
        sliderInput("ec50_neuro", "EC50 Neuro (pg/mL)", 100, 1000, 400, 25)
      ),
      actionButton("simulate", "Run Simulation", class="btn btn-primary btn-block",
                   icon=icon("play-circle")),
      br(),
      downloadButton("download_data", "Download Results (CSV)")
    ),

    mainPanel(
      width=9,
      tabsetPanel(
        id="main_tabs",

        # TAB 1: Patient Profile
        tabPanel(
          "1. 환자 프로파일",
          icon=icon("user-md"),
          br(),
          fluidRow(
            column(6,
              wellPanel(
                h4("Pathophysiology Status"),
                tableOutput("profile_table")
              )
            ),
            column(6,
              wellPanel(
                h4("PA Staging & Risk"),
                uiOutput("risk_panel")
              )
            )
          ),
          fluidRow(
            column(12,
              wellPanel(
                h4("Disease Mechanism Overview"),
                plotOutput("mechanism_plot", height="350px")
              )
            )
          )
        ),

        # TAB 2: Drug PK
        tabPanel(
          "2. 약물 PK",
          icon=icon("pills"),
          br(),
          fluidRow(
            column(8, plotlyOutput("pk_plasma_plot", height="350px")),
            column(4,
              wellPanel(
                h5("PK Summary"),
                tableOutput("pk_summary")
              )
            )
          ),
          fluidRow(
            column(6, plotlyOutput("holotc_plot", height="300px")),
            column(6, plotlyOutput("liver_stores_plot", height="300px"))
          )
        ),

        # TAB 3: Hematological PD
        tabPanel(
          "3. 혈액 PD",
          icon=icon("tint"),
          br(),
          fluidRow(
            column(6, plotlyOutput("hgb_plot",   height="320px")),
            column(6, plotlyOutput("mcv_plot",   height="320px"))
          ),
          fluidRow(
            column(6, plotlyOutput("retic_plot", height="320px")),
            column(6,
              wellPanel(
                h5("Hematological Milestones"),
                tableOutput("heme_milestones")
              )
            )
          )
        ),

        # TAB 4: Neurological PD
        tabPanel(
          "4. 신경계 PD",
          icon=icon("brain"),
          br(),
          fluidRow(
            column(8, plotlyOutput("neuro_plot", height="380px")),
            column(4,
              wellPanel(
                h5("SCD Staging"),
                uiOutput("scd_stage"),
                hr(),
                h5("Nerve B12 Pool"),
                plotlyOutput("nerve_b12_plot", height="200px")
              )
            )
          ),
          fluidRow(
            column(12,
              wellPanel(
                h5("Neurological Recovery Assessment (at 6 & 12 months)"),
                tableOutput("neuro_recovery_table")
              )
            )
          )
        ),

        # TAB 5: Biomarkers
        tabPanel(
          "5. 바이오마커",
          icon=icon("flask"),
          br(),
          fluidRow(
            column(6, plotlyOutput("mma_plot",     height="300px")),
            column(6, plotlyOutput("hcy_plot",     height="300px"))
          ),
          fluidRow(
            column(6, plotlyOutput("if_func_plot", height="300px")),
            column(6, plotlyOutput("autoab_plot",  height="300px"))
          ),
          fluidRow(
            column(12,
              wellPanel(
                h5("Biomarker Summary Table"),
                DTOutput("biomarker_table")
              )
            )
          )
        ),

        # TAB 6: Scenario Comparison
        tabPanel(
          "6. 시나리오 비교",
          icon=icon("chart-bar"),
          br(),
          fluidRow(
            column(4,
              wellPanel(
                h5("Select Scenarios to Compare"),
                checkboxGroupInput("compare_scenarios",
                                   label=NULL,
                                   choices=names(scenario_presets),
                                   selected=names(scenario_presets))
              )
            ),
            column(8,
              selectInput("compare_endpoint", "Endpoint to Compare",
                          choices=c("Hemoglobin (Hb)"="HGB",
                                    "MCV (fL)"="MCV",
                                    "Plasma B12 (pg/mL)"="Cp_pgmL",
                                    "Neurological Score"="NEURO",
                                    "MMA (µmol/L)"="MMA",
                                    "Homocysteine (µmol/L)"="HCY",
                                    "HoloTC (pmol/L)"="HOLOTC"),
                          selected="HGB"),
              plotlyOutput("compare_plot", height="360px")
            )
          ),
          fluidRow(
            column(12,
              wellPanel(
                h5("Outcome Summary at 6 and 12 Months"),
                DTOutput("compare_table")
              )
            )
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: patient parameters
  pa_params <- reactive({
    parietal_fn <- switch(input$pa_severity,
      "mild"     = 0.40,
      "moderate" = 0.20,
      "severe"   = 0.05
    )
    list(
      HGB0       = input$baseline_hgb,
      NEURO0     = input$baseline_neuro,
      MCV0       = input$baseline_mcv,
      PARIETAL0  = parietal_fn,
      F_PASSIVE  = input$f_passive / 100,
      EC50_HGB   = input$ec50_hgb,
      EC50_NEURO = input$ec50_neuro
    )
  })

  # ---- Reactive: simulation
  sim_result <- eventReactive(input$simulate, {
    req(input$dose_ug, input$freq, input$sim_months)
    tsim <- as.numeric(input$sim_months) * 30 * 24

    pars <- pa_params()

    ev_df <- build_events(
      route     = input$route,
      dose_mg   = as.numeric(input$dose_ug),
      freq_h    = as.numeric(input$freq),
      n_doses   = ceiling(tsim / as.numeric(input$freq)),
      load_dose = input$load_dose,
      load_n    = input$load_n,
      load_freq_h = as.numeric(input$load_freq)
    )

    init_vals <- list(
      HGB    = pars$HGB0,
      MCV    = pars$MCV0,
      NEURO  = input$baseline_neuro
    )

    mod_custom <- pa_mod %>%
      param(F_PASSIVE=pars$F_PASSIVE, EC50_HGB=pars$EC50_HGB,
            EC50_NEURO=pars$EC50_NEURO, HGB0=pars$HGB0,
            MCV0=pars$MCV0, NEURO0=input$baseline_neuro)

    if (nrow(ev_df)==0) {
      out <- mod_custom %>% mrgsim(end=tsim, delta=12) %>% as.data.frame()
    } else {
      out <- mod_custom %>% ev(ev_df) %>% mrgsim(end=tsim, delta=12) %>% as.data.frame()
    }
    out$time_days <- out$time/24
    out
  })

  # ---- Reactive: all scenarios
  all_scenarios <- eventReactive(input$simulate, {
    tsim <- 365*24
    lapply(names(scenario_presets), function(sname) {
      sp <- scenario_presets[[sname]]
      ev_df <- build_events(sp$route, sp$dose, sp$freq, sp$n,
                            sp$load, sp$load_n, sp$load_freq)
      out <- run_sim(ev_df, tsim=tsim, delta=12)
      out$Scenario <- sname
      out
    }) %>% bind_rows()
  })

  # ============================================================
  # TAB 1: Patient Profile
  # ============================================================
  output$profile_table <- renderTable({
    pars <- pa_params()
    data.frame(
      Parameter = c("Baseline Hb", "Baseline MCV", "Baseline Neuro Score",
                    "Parietal Cell Function", "Passive Oral F", "PA Severity"),
      Value     = c(paste(pars$HGB0, "g/dL"), paste(pars$MCV0, "fL"),
                    paste(input$baseline_neuro, "/10"),
                    paste0(pars$PARIETAL0*100, "%"),
                    paste0(pars$F_PASSIVE*100, "%"),
                    input$pa_severity),
      Normal    = c(">12 g/dL", "80-100 fL", "0/10", "~100%", "1%", "—")
    )
  })

  output$risk_panel <- renderUI({
    sev <- input$pa_severity
    risk_col <- switch(sev, "mild"="#F39C12", "moderate"="#E74C3C", "severe"="#922B21")
    tagList(
      div(style=paste0("background:", risk_col, "; color:white; padding:10px; border-radius:6px; margin-bottom:10px;"),
        h5(paste("Severity:", toupper(sev)), style="margin:0; font-weight:bold")
      ),
      tags$ul(
        tags$li("Autoimmune antibody burden: HIGH (80%)"),
        tags$li("IF production: significantly reduced"),
        tags$li("Risk of gastric carcinoid: elevated"),
        tags$li("Neurological involvement: potential SCD"),
        tags$li("Associated autoimmune diseases: thyroid, T1DM")
      )
    )
  })

  output$mechanism_plot <- renderPlot({
    # Simple schematic bar chart of B12 pathway efficiency in PA vs. Normal
    df <- data.frame(
      Stage     = factor(c("Dietary\nIntake","Gastric\nRelease","IF\nProduction",
                           "IF-B12\nComplex","Ileal\nAbsorption","Plasma\nB12",
                           "Tissue\nDelivery"),
                         levels=c("Dietary\nIntake","Gastric\nRelease","IF\nProduction",
                                  "IF-B12\nComplex","Ileal\nAbsorption","Plasma\nB12","Tissue\nDelivery")),
      Normal_pct= c(100, 95, 100, 95, 90, 100, 90),
      PA_pct    = c(100, 40, 20, 10, 5, 3, 2)
    ) %>% pivot_longer(c(Normal_pct, PA_pct), names_to="Condition", values_to="Efficiency")
    df$Condition <- recode(df$Condition, "Normal_pct"="Normal", "PA_pct"="PA (Untreated)")
    ggplot(df, aes(x=Stage, y=Efficiency, fill=Condition)) +
      geom_col(position="dodge", alpha=0.85) +
      scale_fill_manual(values=c("Normal"="#2ECC71","PA (Untreated)"="#E74C3C")) +
      labs(title="Cobalamin Absorption Efficiency Along the Pathway",
           x=NULL, y="Relative Efficiency (%)", fill=NULL) +
      theme_bw(base_size=11) + theme(legend.position="bottom")
  })

  # ============================================================
  # TAB 2: Drug PK
  # ============================================================
  output$pk_plasma_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=Cp_pgmL)) +
      geom_line(color="#2980B9", linewidth=1) +
      geom_hline(yintercept=200, linetype="dashed", color="red", alpha=0.7) +
      geom_hline(yintercept=900, linetype="dashed", color="green4", alpha=0.7) +
      annotate("text", x=max(df$time_days)*0.8, y=230, label="Deficiency (<200)", color="red", size=3) +
      labs(title="Plasma Cobalamin (B12)", x="Time (days)", y="Plasma B12 (pg/mL)") +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$pk_summary <- renderTable({
    req(sim_result())
    df <- sim_result()
    df_pk <- df %>% filter(time_days > 0)
    data.frame(
      Metric      = c("Peak Plasma B12 (pg/mL)","Trough Plasma B12",
                      "Time to B12 >200 pg/mL (days)","Final HoloTC (pmol/L)"),
      Value       = c(round(max(df_pk$Cp_pgmL), 0),
                      round(min(df_pk$Cp_pgmL[df_pk$time_days > 10]), 0),
                      round(min(df_pk$time_days[df_pk$Cp_pgmL > 200], na.rm=TRUE), 1),
                      round(tail(df$HOLOTC, 1), 1))
    )
  })

  output$holotc_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=HOLOTC)) +
      geom_line(color="#9B59B6", linewidth=1) +
      geom_hline(yintercept=35, linetype="dashed", color="red") +
      labs(title="HoloTranscobalamin II (Active B12)", x="Time (days)", y="HoloTC (pmol/L)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$liver_stores_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=LIVER)) +
      geom_line(color="#27AE60", linewidth=1) +
      geom_hline(yintercept=2000, linetype="dashed", color="gray40") +
      labs(title="Hepatic B12 Stores", x="Time (days)", y="Liver B12 (µg)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  # ============================================================
  # TAB 3: Hematological PD
  # ============================================================
  output$hgb_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=HGB)) +
      geom_line(color="#E74C3C", linewidth=1) +
      geom_hline(yintercept=12, linetype="dashed", color="gray40") +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=12, alpha=0.05, fill="red") +
      coord_cartesian(ylim=c(4, 16)) +
      labs(title="Hemoglobin", x="Time (days)", y="Hb (g/dL)") +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$mcv_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=MCV)) +
      geom_line(color="#E67E22", linewidth=1) +
      geom_hline(yintercept=100, linetype="dashed", color="gray40") +
      geom_hline(yintercept=80, linetype="dotted", color="gray60") +
      coord_cartesian(ylim=c(75, 130)) +
      labs(title="Mean Corpuscular Volume (MCV)", x="Time (days)", y="MCV (fL)") +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$retic_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% filter(time_days <= 30)
    p <- ggplot(df, aes(x=time_days, y=RETIC)) +
      geom_line(color="#3498DB", linewidth=1.2) +
      geom_vline(xintercept=7, linetype="dotted", color="darkred") +
      annotate("text", x=7.5, y=max(df$RETIC, na.rm=TRUE)*0.9,
               label="Crisis peak\n(day 4-10)", color="darkred", size=3) +
      labs(title="Reticulocyte Crisis (First 30 Days)", x="Time (days)",
           y="Reticulocytes (×10⁹/L)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$heme_milestones <- renderTable({
    req(sim_result())
    df <- sim_result()
    t_hgb12 <- min(df$time_days[df$HGB >= 12], na.rm=TRUE)
    t_mcv_nl <- min(df$time_days[df$MCV <= 100], na.rm=TRUE)
    retic_pk <- max(df$RETIC[df$time_days <= 30], na.rm=TRUE)
    data.frame(
      Milestone       = c("Reticulocyte Peak", "Time to Hb ≥12 g/dL", "Time to MCV ≤100 fL"),
      Value           = c(paste0(round(retic_pk,0), " ×10⁹/L"),
                          ifelse(is.finite(t_hgb12), paste0(round(t_hgb12,0), " days"), "Not reached"),
                          ifelse(is.finite(t_mcv_nl), paste0(round(t_mcv_nl,0), " days"), "Not reached"))
    )
  })

  # ============================================================
  # TAB 4: Neurological PD
  # ============================================================
  output$neuro_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=NEURO)) +
      geom_line(color="#1ABC9C", linewidth=1) +
      geom_hline(yintercept=5, linetype="dashed", color="orange", linewidth=0.8) +
      geom_hline(yintercept=2, linetype="dashed", color="green4", linewidth=0.8) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=5, ymax=10, alpha=0.08, fill="red") +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=2, alpha=0.08, fill="green") +
      coord_cartesian(ylim=c(0, 10)) +
      labs(title="Neurological Disability Score (0=none, 10=severe)",
           x="Time (days)", y="Neuro Score (0-10)") +
      theme_bw(base_size=11)
    ggplotly(p) %>% layout(hovermode="x unified")
  })

  output$scd_stage <- renderUI({
    req(sim_result())
    df <- sim_result()
    score_now <- df$NEURO[df$time_days == max(df$time_days)][1]
    stage <- if (score_now < 2) "Mild / Subclinical"
             else if (score_now < 4) "Mild-Moderate"
             else if (score_now < 6) "Moderate SCD"
             else "Severe SCD"
    col <- if (score_now < 2) "#2ECC71" else if (score_now < 4) "#F39C12"
            else if (score_now < 6) "#E74C3C" else "#7B241C"
    div(style=paste0("background:", col, "; color:white; padding:8px; border-radius:5px;"),
      strong(stage), br(), paste0("Score: ", round(score_now, 1), "/10")
    )
  })

  output$nerve_b12_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=NERVE)) +
      geom_line(color="#1ABC9C", linewidth=1) +
      labs(title="Nerve B12 Pool", x="Time (days)", y="Nerve B12 (µg)") +
      theme_bw(base_size=9)
    ggplotly(p)
  })

  output$neuro_recovery_table <- renderTable({
    req(sim_result())
    df <- sim_result()
    get_neuro <- function(t_days) {
      nearest <- df[which.min(abs(df$time_days - t_days)), ]
      round(nearest$NEURO, 2)
    }
    data.frame(
      Timepoint   = c("Baseline", "1 Month", "3 Months", "6 Months", "12 Months"),
      Neuro_Score = sapply(c(0, 30, 90, 180, 365), get_neuro),
      Nerve_B12   = sapply(c(0, 30, 90, 180, 365), function(t) {
        round(df$NERVE[which.min(abs(df$time_days - t))], 3)
      })
    )
  })

  # ============================================================
  # TAB 5: Biomarkers
  # ============================================================
  output$mma_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=MMA)) +
      geom_line(color="#8E44AD", linewidth=1) +
      geom_hline(yintercept=0.4, linetype="dashed", color="red") +
      labs(title="Methylmalonic Acid (MMA)", x="Time (days)", y="MMA (µmol/L)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$hcy_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=HCY)) +
      geom_line(color="#D35400", linewidth=1) +
      geom_hline(yintercept=15, linetype="dashed", color="red") +
      labs(title="Total Homocysteine (tHcy)", x="Time (days)", y="Hcy (µmol/L)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$if_func_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=PARIETAL*100)) +
      geom_line(color="#27AE60", linewidth=1) +
      coord_cartesian(ylim=c(0, 100)) +
      labs(title="Parietal Cell Function (%)", x="Time (days)", y="% Function") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$autoab_plot <- renderPlotly({
    req(sim_result())
    df <- sim_result()
    p <- ggplot(df, aes(x=time_days, y=AUTO_AB*100)) +
      geom_line(color="#C0392B", linewidth=1) +
      coord_cartesian(ylim=c(0, 100)) +
      labs(title="Autoimmune Antibody Burden (%)", x="Time (days)",
           y="Relative Ab burden (%)") +
      theme_bw(base_size=11)
    ggplotly(p)
  })

  output$biomarker_table <- renderDT({
    req(sim_result())
    df <- sim_result() %>%
      filter(time_days %in% c(0, 14, 30, 90, 180, 365)) %>%
      mutate(across(c(Cp_pgmL, HGB, MCV, RETIC, NEURO, HOLOTC, MMA, HCY,
                      PARIETAL, AUTO_AB), ~round(.x, 2))) %>%
      select(time_days, Cp_pgmL, HGB, MCV, RETIC, NEURO, HOLOTC, MMA, HCY) %>%
      rename("Day"="time_days", "B12 (pg/mL)"="Cp_pgmL", "Hb (g/dL)"="HGB",
             "MCV (fL)"="MCV", "Retic (×10⁹/L)"="RETIC", "Neuro"="NEURO",
             "HoloTC (pmol/L)"="HOLOTC", "MMA (µmol/L)"="MMA", "Hcy (µmol/L)"="HCY")
    datatable(df, options=list(dom="t", pageLength=10), rownames=FALSE) %>%
      formatStyle("Hb (g/dL)", backgroundColor=styleInterval(12, c("#FADBD8","#D5F5E3"))) %>%
      formatStyle("B12 (pg/mL)", backgroundColor=styleInterval(200, c("#FADBD8","#D5F5E3"))) %>%
      formatStyle("MMA (µmol/L)", backgroundColor=styleInterval(0.4, c("#D5F5E3","#FADBD8")))
  }, server=FALSE)

  # ============================================================
  # TAB 6: Scenario Comparison
  # ============================================================
  output$compare_plot <- renderPlotly({
    req(all_scenarios())
    df <- all_scenarios() %>%
      filter(Scenario %in% input$compare_scenarios)
    ep <- input$compare_endpoint
    colors_s <- c(
      "S1: Untreated PA"   = "#E74C3C",
      "S2: IM Standard"    = "#2ECC71",
      "S3: IM Maintenance" = "#3498DB",
      "S4: Oral HD"        = "#F39C12",
      "S5: Aggressive IM"  = "#9B59B6"
    )
    ref_lines <- list(
      "HGB"="12", "MCV"="100", "Cp_pgmL"="200",
      "HOLOTC"="35", "MMA"="0.4", "HCY"="15"
    )
    p <- ggplot(df, aes_string(x="time_days", y=ep, color="Scenario")) +
      geom_line(linewidth=0.9) +
      scale_color_manual(values=colors_s) +
      labs(title=paste("Scenario Comparison:", ep),
           x="Time (days)", y=ep, color=NULL) +
      theme_bw(base_size=11) + theme(legend.position="bottom")
    if (!is.null(ref_lines[[ep]])) {
      p <- p + geom_hline(yintercept=as.numeric(ref_lines[[ep]]),
                          linetype="dashed", color="gray40", linewidth=0.7)
    }
    ggplotly(p) %>% layout(hovermode="x unified", legend=list(orientation="h"))
  })

  output$compare_table <- renderDT({
    req(all_scenarios())
    df <- all_scenarios() %>%
      filter(Scenario %in% input$compare_scenarios)
    tbl <- df %>%
      filter(abs(time_days - 180) < 1 | abs(time_days - 365) < 1) %>%
      mutate(Timepoint=ifelse(time_days < 200, "6 months", "12 months")) %>%
      group_by(Scenario, Timepoint) %>%
      summarise(across(c(Cp_pgmL, HGB, MCV, NEURO, MMA, HCY, HOLOTC), mean), .groups="drop") %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(tbl, options=list(pageLength=15, dom="t"), rownames=FALSE)
  }, server=FALSE)

  # ---- Download
  output$download_data <- downloadHandler(
    filename=function() paste0("PA_QSP_results_", Sys.Date(), ".csv"),
    content=function(file) {
      req(sim_result())
      write.csv(sim_result(), file, row.names=FALSE)
    }
  )
}

# ============================================================
# RUN APP
# ============================================================
shinyApp(ui=ui, server=server)
