## ============================================================
## Primary Biliary Cholangitis (PBC) — Interactive Shiny Dashboard
## Author: Claude Code Routine (CCR)  |  Date: 2026-06-17
##
## 6 Tabs:
##   1. Patient Profile & Disease Stage
##   2. PK Profiles & Drug Concentrations
##   3. PD Biomarkers (ALP, GGT, Bilirubin, IgM)
##   4. Disease Biology (AMA, Th1, BEC, Fibrosis)
##   5. Clinical Endpoints & Risk Scores
##   6. Scenario Comparison & Treatment Summary
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────
## MODEL SIMULATION FUNCTION
## ─────────────────────────────────────────────────────────────

simulate_pbc <- function(
  # Patient parameters
  body_weight = 65,
  ALP0 = 280, BILI0 = 1.2, GGT0 = 120,
  AMA0 = 1.0, FIB0 = 1.5, IgM0 = 350,
  ULN_ALP = 120, ULN_BILI = 1.0, ULN_GGT = 55,

  # Treatment flags
  use_UDCA = TRUE, use_OCA = FALSE,
  use_ELF = FALSE, use_SEL = FALSE, use_BEZ = FALSE,

  # OCA dose (5 or 10 mg)
  OCA_dose_mg = 5,

  # Simulation duration (months)
  duration = 24,

  # Scenario label
  scenario_name = "Treatment"
) {
  # Steady-state concentrations (mg/L)
  C_UDCA <- if (use_UDCA) 0.60 * 975 / 40     else 0
  C_OCA  <- if (use_OCA)  0.35 * OCA_dose_mg / 8 else 0
  C_ELF  <- if (use_ELF)  0.72 * 80  / 15     else 0
  C_SEL  <- if (use_SEL)  0.60 * 10  / 12     else 0
  C_BEZ  <- if (use_BEZ)  1.00 * 400 / 5      else 0

  # Fixed model parameters
  kFIB     <- 0.004; kBEC <- 0.02; kBEC_rep <- 0.005
  kAMA     <- 0.005; kTh1 <- 0.03; kTh1_die <- 0.025; kTreg <- 0.015
  kBA_syn  <- 0.5;  kBA_detox <- 0.1
  kFGF19   <- 0.02; FGF19_ss <- 100
  EC50_UDCA_ALP <- 150; Emax_UDCA_ALP <- 0.40
  EC50_OCA_FXR  <- 0.5; Emax_OCA_ALP  <- 0.38
  EC50_OCA_FGF  <- 0.3
  EC50_ELF_ALP  <- 0.8; Emax_ELF_ALP  <- 0.55
  EC50_SEL_ALP  <- 0.3; Emax_SEL_ALP  <- 0.30
  EC50_BEZ_ALP  <- 2.0; Emax_BEZ_ALP  <- 0.67
  EC50_SEL_GGT  <- 0.2; Emax_SEL_GGT  <- 0.50
  Emax_UDCA_FIB <- 0.20; Emax_OCA_FIB <- 0.25; Emax_ELF_FIB <- 0.30
  Emax_ELF_ITCH <- 0.70; EC50_ELF_ITCH <- 0.6

  # PD effects
  E_UDCA_ALP  <- Emax_UDCA_ALP * C_UDCA / (EC50_UDCA_ALP + C_UDCA)
  E_OCA_ALP   <- Emax_OCA_ALP  * C_OCA  / (EC50_OCA_FXR  + C_OCA)
  E_ELF_ALP   <- Emax_ELF_ALP  * C_ELF  / (EC50_ELF_ALP  + C_ELF)
  E_SEL_ALP   <- Emax_SEL_ALP  * C_SEL  / (EC50_SEL_ALP  + C_SEL)
  E_BEZ_ALP   <- Emax_BEZ_ALP  * C_BEZ  / (EC50_BEZ_ALP  + C_BEZ)
  E_ALP_total  <- min(E_UDCA_ALP + E_OCA_ALP + E_ELF_ALP + E_SEL_ALP + E_BEZ_ALP, 0.85)

  E_SEL_GGT  <- Emax_SEL_GGT * C_SEL / (EC50_SEL_GGT + C_SEL)
  E_ELF_GGT  <- 0.35 * C_ELF / (0.9 + C_ELF)
  E_GGT_total <- min(E_SEL_GGT + E_ELF_GGT + 0.2 * E_UDCA_ALP, 0.70)

  E_FIB_total <- min(
    Emax_UDCA_FIB * C_UDCA / (EC50_UDCA_ALP + C_UDCA) +
    Emax_OCA_FIB  * C_OCA  / (EC50_OCA_FXR  + C_OCA)  +
    Emax_ELF_FIB  * C_ELF  / (EC50_ELF_ALP  + C_ELF), 0.55)

  E_ELF_itch <- Emax_ELF_ITCH * C_ELF / (EC50_ELF_ITCH + C_ELF)
  E_SEL_itch <- 0.50 * C_SEL / (0.4 + C_SEL)
  OCA_itch_pen <- 0.30 * C_OCA / (0.5 + C_OCA)  # OCA causes pruritus

  # ODE simulation
  dt <- 0.25  # months
  t  <- seq(0, duration, by = dt)
  n  <- length(t)

  AMA      <- numeric(n); AMA[1]      <- AMA0
  Th1      <- numeric(n); Th1[1]      <- 1.0
  BEC      <- numeric(n); BEC[1]      <- 0.30
  BA_toxic <- numeric(n); BA_toxic[1] <- 1.0
  FGF19    <- numeric(n); FGF19[1]    <- FGF19_ss
  ALP      <- numeric(n); ALP[1]      <- ALP0
  BILI     <- numeric(n); BILI[1]     <- BILI0
  GGT      <- numeric(n); GGT[1]      <- GGT0
  FIB      <- numeric(n); FIB[1]      <- FIB0
  IgM      <- numeric(n); IgM[1]      <- IgM0

  for (i in 2:n) {
    FGF19_new <- FGF19[i-1] + dt * (
      FGF19_ss * kFGF19 * (1 + 2.5 * C_OCA / (EC50_OCA_FGF + C_OCA) + 0.3 * E_UDCA_ALP) -
      kFGF19 * FGF19[i-1]
    )
    FGF19_effect <- FGF19_new / (FGF19_ss + FGF19_new)

    BA_synth <- kBA_syn * (1 - 0.4 * FGF19_effect)
    BA_detox  <- kBA_detox * (1 + 1.5 * (E_ELF_ALP + E_SEL_ALP + E_BEZ_ALP) + 0.5 * E_UDCA_ALP)
    BA_new    <- max(0.1, BA_toxic[i-1] + dt * (BA_synth - BA_detox * BA_toxic[i-1]))

    AMA_new <- max(0.01, AMA[i-1] + dt * (
      0.3 * BEC[i-1] - kAMA * AMA[i-1] - 0.1 * E_UDCA_ALP * AMA[i-1]
    ))

    Th1_new <- max(0.01, Th1[i-1] + dt * (
      kTh1 * AMA_new * (1 + 0.5 * BEC[i-1]) -
      (kTh1_die + kTreg + 0.2 * E_UDCA_ALP) * Th1[i-1]
    ))

    BEC_dmg_rate <- kBEC * (0.5 * Th1_new + 0.3 * AMA_new + 0.2 * BA_new)
    BEC_rep_rate <- kBEC_rep * (1 - BEC[i-1]) + 0.25 * E_UDCA_ALP * BEC[i-1]
    BEC_new <- max(0, min(1, BEC[i-1] + dt * (BEC_dmg_rate * (1 - BEC[i-1]) - BEC_rep_rate)))

    ALP_target <- ALP0 * (1 + 0.5 * BEC_new + 0.3 * FIB[i-1] / 4) * (1 - E_ALP_total)
    ALP_new    <- max(40, ALP[i-1] + dt * (ALP_target - ALP[i-1]) / 3)

    BILI_target <- BILI0 * (1 + 0.8 * FIB[i-1] / 4 + 0.4 * BEC_new) *
                   (1 - 0.25 * (E_UDCA_ALP + E_OCA_ALP + E_ELF_ALP))
    BILI_new    <- max(0.3, BILI[i-1] + dt * (BILI_target - BILI[i-1]) / 4)

    GGT_target <- GGT0 * (1 + 0.6 * BEC_new + 0.4 * BA_new) * (1 - E_GGT_total)
    GGT_new    <- max(20, GGT[i-1] + dt * (GGT_target - GGT[i-1]) / 2)

    FIB_prog    <- kFIB * (0.4 * BEC_new + 0.3 * BA_new + 0.2 * Th1_new + 0.1) * max(0, 4 - FIB[i-1])
    FIB_regress <- 0.001 * E_FIB_total * max(FIB[i-1] - 1, 0)
    FIB_new     <- max(0, min(4, FIB[i-1] + dt * (FIB_prog - FIB_regress)))

    IgM_drive  <- IgM0 * (0.8 + 0.4 * AMA_new)
    IgM_target <- IgM_drive * (1 - 0.15 * E_UDCA_ALP - 0.1 * E_OCA_ALP)
    IgM_new    <- max(50, IgM[i-1] + dt * (IgM_target - IgM[i-1]) / 6)

    FGF19[i] <- max(50, FGF19_new)
    BA_toxic[i] <- BA_new
    AMA[i]   <- AMA_new
    Th1[i]   <- Th1_new
    BEC[i]   <- BEC_new
    ALP[i]   <- ALP_new
    BILI[i]  <- BILI_new
    GGT[i]   <- GGT_new
    FIB[i]   <- FIB_new
    IgM[i]   <- IgM_new
  }

  # Derived quantities
  Pruritus_base <- 4.0 * BA_toxic * (1 - 0.2 * BEC)
  Pruritus_treat <- Pruritus_base * (1 - E_ELF_itch - E_SEL_itch) + 2.0 * OCA_itch_pen * Pruritus_base
  Pruritus <- pmax(0, pmin(10, Pruritus_treat))

  # PK profiles (simplified daily oscillation for display)
  t_fine <- seq(0, min(3, duration), by = 1/24)  # first 3 months, hourly
  pk_data <- data.frame(
    time_h = t_fine * 24 * 30,  # convert to hours (approximate)
    UDCA   = if (use_UDCA) C_UDCA * (1 + 0.4 * sin(2 * pi * t_fine * 24 * 30 / 24)) else 0,
    OCA    = if (use_OCA)  C_OCA  * (1 + 0.5 * sin(2 * pi * t_fine * 24 * 30 / 24)) else 0,
    ELF    = if (use_ELF)  C_ELF  * (1 + 0.3 * sin(2 * pi * t_fine * 24 * 30 / 24)) else 0,
    SEL    = if (use_SEL)  C_SEL  * (1 + 0.4 * sin(2 * pi * t_fine * 24 * 30 / 24)) else 0,
    BEZ    = if (use_BEZ)  C_BEZ  * (1 + 0.35 * sin(2 * pi * t_fine * 24 * 30 / 24)) else 0
  )

  list(
    time     = t,
    ALP      = ALP,     ALP_xULN = ALP / ULN_ALP,
    BILI     = BILI,    BILI_xULN = BILI / ULN_BILI,
    GGT      = GGT,     GGT_xULN = GGT / ULN_GGT,
    IgM      = IgM,     FGF19    = FGF19,
    AMA      = AMA,     Th1      = Th1,
    BEC_dmg  = BEC,     BA_toxic = BA_toxic,
    Fibrosis = FIB,     Pruritus = Pruritus,
    ALP_norm = as.numeric(ALP / ULN_ALP <= 1.0),
    Paris2   = as.numeric((ALP / ULN_ALP <= 1.5) & (BILI / ULN_BILI <= 1.0)),
    GLOBE_surrogate = 0.5 * (ALP / ULN_ALP - 1) + 0.8 * log(BILI + 0.1) + 0.3 * FIB,
    C_UDCA = C_UDCA, C_OCA = C_OCA, C_ELF = C_ELF, C_SEL = C_SEL, C_BEZ = C_BEZ,
    pk_data  = pk_data,
    scenario = scenario_name,
    ULN_ALP  = ULN_ALP, ULN_BILI = ULN_BILI, ULN_GGT = ULN_GGT
  )
}

## ─────────────────────────────────────────────────────────────
## SHINY UI
## ─────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = "PBC QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",       tabName = "tab_patient",   icon = icon("user")),
      menuItem("② PK Profiles",           tabName = "tab_pk",        icon = icon("chart-line")),
      menuItem("③ PD Biomarkers",         tabName = "tab_pd",        icon = icon("flask")),
      menuItem("④ Disease Biology",       tabName = "tab_biology",   icon = icon("dna")),
      menuItem("⑤ Clinical Endpoints",    tabName = "tab_endpoints", icon = icon("hospital")),
      menuItem("⑥ Scenario Comparison",   tabName = "tab_compare",   icon = icon("balance-scale"))
    ),

    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#ccc;"),

    sliderInput("ALP0",   "Baseline ALP (IU/L)",    60,  1000, 280, step = 10),
    sliderInput("BILI0",  "Baseline Bilirubin (mg/dL)", 0.3, 6, 1.2, step = 0.1),
    sliderInput("GGT0",   "Baseline GGT (IU/L)",    30,  500,  120, step = 10),
    sliderInput("FIB0",   "Baseline Fibrosis (F0-F4)", 0, 3.5, 1.5, step = 0.5),
    sliderInput("IgM0",   "Baseline IgM (mg/dL)",   80,  800,  350, step = 10),

    hr(),
    h5("Treatment Selection", style = "padding-left:15px; color:#ccc;"),
    checkboxInput("use_UDCA", "UDCA 975 mg/day", value = TRUE),
    checkboxInput("use_OCA",  "Obeticholic acid (OCA)", value = FALSE),
    radioButtons("OCA_dose", NULL,
                 choices = c("5 mg/day" = 5, "10 mg/day" = 10), inline = TRUE, selected = 5),
    checkboxInput("use_ELF", "Elafibranor 80 mg/day (FDA 2024)", value = FALSE),
    checkboxInput("use_SEL", "Seladelpar 10 mg/day (FDA 2024)",  value = FALSE),
    checkboxInput("use_BEZ", "Bezafibrate 400 mg/day",           value = FALSE),

    hr(),
    sliderInput("duration", "Simulation Duration (months)", 6, 60, 24, step = 6),

    actionButton("run", "▶ Run Simulation", class = "btn-primary btn-block",
                 style = "margin: 10px 15px; width: 250px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-title { font-size: 14px; }
      .info-box-content { font-size: 13px; }
    "))),

    tabItems(

      # ── Tab 1: Patient Profile ──────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("box_ALP",  width = 3),
          valueBoxOutput("box_BILI", width = 3),
          valueBoxOutput("box_GGT",  width = 3),
          valueBoxOutput("box_FIB",  width = 3)
        ),
        fluidRow(
          box(title = "Disease Stage Summary", width = 6, status = "purple", solidHeader = TRUE,
            tableOutput("tbl_stage")
          ),
          box(title = "PBC Disease Overview", width = 6, status = "primary", solidHeader = TRUE,
            tags$ul(
              tags$li(strong("Epidemiology:"), " 1/1000 women >40yr; 9:1 F:M ratio"),
              tags$li(strong("Autoantigen:"), " PDC-E2 (pyruvate dehydrogenase complex)"),
              tags$li(strong("Hallmark:"), " AMA ≥1:40 (≥95% sensitivity)"),
              tags$li(strong("Cholestasis:"), " ALP, GGT, IgM elevation; pruritus"),
              tags$li(strong("Staging:"), " Metavir F0-F4; ductopenia → cirrhosis"),
              tags$li(strong("FDA 2024:"), " Elafibranor (ELATIVE) + Seladelpar (RESPONSE)"),
              tags$li(strong("GLOBE score:"), " Predicts transplant-free survival at 5yr")
            )
          )
        ),
        fluidRow(
          box(title = "Drug Information", width = 12, status = "success", solidHeader = TRUE,
            DT::dataTableOutput("tbl_drugs")
          )
        )
      ),

      # ── Tab 2: PK Profiles ─────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Concentration Profiles (first 72 hours)", width = 12,
              status = "info", solidHeader = TRUE,
            plotlyOutput("plt_pk", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Steady-State Concentrations & PD Parameters", width = 12,
              status = "primary", solidHeader = TRUE,
            DT::dataTableOutput("tbl_pk_summary")
          )
        )
      ),

      # ── Tab 3: PD Biomarkers ──────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "ALP Over Time (×ULN)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plt_ALP", height = "300px")
          ),
          box(title = "Bilirubin Over Time (×ULN)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_BILI", height = "300px")
          )
        ),
        fluidRow(
          box(title = "GGT Over Time (×ULN)", width = 6, status = "info", solidHeader = TRUE,
            plotlyOutput("plt_GGT", height = "300px")
          ),
          box(title = "Serum IgM (mg/dL)", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plt_IgM", height = "300px")
          )
        ),
        fluidRow(
          box(title = "FGF19 (pg/mL)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plt_FGF19", height = "300px"),
            helpText("FGF19 rises with FXR activation (OCA); suppresses CYP7A1 → reduced BA synthesis")
          ),
          box(title = "Pruritus NRS (0-10)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_pruritus", height = "300px"),
            helpText("OCA may worsen pruritus; PPARδ agonists (ELF, SEL) improve pruritus")
          )
        )
      ),

      # ── Tab 4: Disease Biology ────────────────────────────
      tabItem(tabName = "tab_biology",
        fluidRow(
          box(title = "AMA Titer (normalized)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plt_AMA", height = "300px"),
            helpText("Anti-mitochondrial antibody — hallmark of PBC; correlates with BEC damage")
          ),
          box(title = "Th1 Cell Activity (normalized)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_Th1", height = "300px"),
            helpText("CD4+ Th1 cells drive IFN-γ → BEC cytotoxicity")
          )
        ),
        fluidRow(
          box(title = "Biliary Epithelial Cell Damage Index (0-1)", width = 6, status = "info", solidHeader = TRUE,
            plotlyOutput("plt_BEC", height = "300px"),
            helpText("0 = intact, 1 = complete ductopenia")
          ),
          box(title = "Toxic BA Pool (normalized)", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plt_BA", height = "300px"),
            helpText("Hydrophobic bile acids drive cholestasis & NLRP3 activation")
          )
        ),
        fluidRow(
          box(title = "Hepatic Fibrosis Score (Metavir 0-4)", width = 12,
              status = "purple", solidHeader = TRUE,
            plotlyOutput("plt_FIB", height = "300px"),
            helpText("F3-F4: portal hypertension risk. Treatment slows/arrests progression; partial regression possible at F2+")
          )
        )
      ),

      # ── Tab 5: Clinical Endpoints ──────────────────────────
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          valueBoxOutput("vbox_paris2",     width = 4),
          valueBoxOutput("vbox_alp_norm",   width = 4),
          valueBoxOutput("vbox_globe",      width = 4)
        ),
        fluidRow(
          box(title = "Paris II Biochemical Response Rate Over Time", width = 6,
              status = "success", solidHeader = TRUE,
            plotlyOutput("plt_paris2", height = "300px"),
            helpText("Paris II: ALP <1.5×ULN AND bilirubin ≤ULN at 12 months")
          ),
          box(title = "GLOBE Score Surrogate", width = 6,
              status = "primary", solidHeader = TRUE,
            plotlyOutput("plt_globe", height = "300px"),
            helpText("Lower GLOBE = better transplant-free survival")
          )
        ),
        fluidRow(
          box(title = "12-Month Biochemical Summary", width = 12,
              status = "info", solidHeader = TRUE,
            tableOutput("tbl_m12")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Multi-Scenario ALP Comparison", width = 12,
              status = "danger", solidHeader = TRUE,
            plotlyOutput("plt_compare_ALP", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Multi-Scenario Fibrosis Comparison", width = 6,
              status = "purple", solidHeader = TRUE,
            plotlyOutput("plt_compare_FIB", height = "350px")
          ),
          box(title = "Multi-Scenario Pruritus Comparison", width = 6,
              status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_compare_pruritus", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Comparison Summary Table at Month 12", width = 12,
              status = "primary", solidHeader = TRUE,
            DT::dataTableOutput("tbl_compare")
          )
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────
## SHINY SERVER
## ─────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Reactive: run primary simulation on button click
  sim_result <- eventReactive(input$run, {
    simulate_pbc(
      ALP0       = input$ALP0,
      BILI0      = input$BILI0,
      GGT0       = input$GGT0,
      FIB0       = input$FIB0,
      IgM0       = input$IgM0,
      use_UDCA   = input$use_UDCA,
      use_OCA    = input$use_OCA,
      OCA_dose_mg = as.numeric(input$OCA_dose),
      use_ELF    = input$use_ELF,
      use_SEL    = input$use_SEL,
      use_BEZ    = input$use_BEZ,
      duration   = input$duration,
      scenario_name = "Custom Treatment"
    )
  }, ignoreNULL = FALSE)

  # Reactive: run all preset scenarios for comparison
  all_scenarios <- eventReactive(input$run, {
    presets <- list(
      list(n="No Treatment",          udca=F, oca=F, elf=F, sel=F, bez=F),
      list(n="UDCA Monotherapy",      udca=T, oca=F, elf=F, sel=F, bez=F),
      list(n="UDCA + OCA 5mg",        udca=T, oca=T, elf=F, sel=F, bez=F),
      list(n="UDCA + OCA 10mg",       udca=T, oca=T, elf=F, sel=F, bez=F),
      list(n="UDCA + Elafibranor",    udca=T, oca=F, elf=T, sel=F, bez=F),
      list(n="UDCA + Seladelpar",     udca=T, oca=F, elf=F, sel=T, bez=F),
      list(n="UDCA + Bezafibrate",    udca=T, oca=F, elf=F, sel=F, bez=T)
    )
    oca_doses <- c(5, 5, 5, 10, 5, 5, 5)  # OCA dose per scenario index

    results <- lapply(seq_along(presets), function(j) {
      cfg <- presets[[j]]
      s   <- simulate_pbc(
        ALP0=input$ALP0, BILI0=input$BILI0, GGT0=input$GGT0,
        FIB0=input$FIB0, IgM0=input$IgM0,
        use_UDCA=cfg$udca, use_OCA=cfg$oca,
        OCA_dose_mg=oca_doses[j],
        use_ELF=cfg$elf, use_SEL=cfg$sel, use_BEZ=cfg$bez,
        duration=input$duration, scenario_name=cfg$n
      )
      data.frame(
        time     = s$time,
        scenario = s$scenario,
        ALP_xULN = s$ALP_xULN,
        Fibrosis = s$Fibrosis,
        Pruritus = s$Pruritus,
        Paris2   = s$Paris2,
        ALP_norm = s$ALP_norm,
        BILI_xULN = s$BILI_xULN,
        GGT_xULN = s$GGT_xULN
      )
    })
    do.call(rbind, results)
  }, ignoreNULL = FALSE)

  # ── Tab 1: Patient Profile ──────────────────────────────────
  output$box_ALP <- renderValueBox({
    s <- sim_result()
    valueBox(sprintf("%.0f IU/L", s$ALP[1]),
             sprintf("ALP (%.1f×ULN)", s$ALP_xULN[1]),
             icon = icon("chart-bar"), color = "red")
  })
  output$box_BILI <- renderValueBox({
    s <- sim_result()
    valueBox(sprintf("%.2f mg/dL", s$BILI[1]),
             sprintf("Bilirubin (%.1f×ULN)", s$BILI_xULN[1]),
             icon = icon("tint"), color = "yellow")
  })
  output$box_GGT <- renderValueBox({
    s <- sim_result()
    valueBox(sprintf("%.0f IU/L", s$GGT[1]),
             sprintf("GGT (%.1f×ULN)", s$GGT_xULN[1]),
             icon = icon("microscope"), color = "orange")
  })
  output$box_FIB <- renderValueBox({
    s <- sim_result()
    fib_label <- c("F0 (No fibrosis)", "F1 (Mild)", "F2 (Moderate)",
                   "F3 (Severe)", "F4 (Cirrhosis)")[floor(s$Fibrosis[1]) + 1]
    valueBox(sprintf("F%.1f", s$Fibrosis[1]),
             paste("Metavir:", fib_label),
             icon = icon("layer-group"), color = "purple")
  })

  output$tbl_stage <- renderTable({
    s <- sim_result()
    tribble(
      ~Parameter,          ~Value,             ~Status,
      "ALP",               sprintf("%.0f IU/L (%.1f×ULN)", s$ALP[1], s$ALP_xULN[1]),
                           ifelse(s$ALP_xULN[1]>1.67,"Abnormal (>1.67×ULN)",ifelse(s$ALP_xULN[1]>1,"Elevated","Normal")),
      "Bilirubin",         sprintf("%.2f mg/dL", s$BILI[1]),
                           ifelse(s$BILI_xULN[1]>1, "Elevated (poor prognosis)", "Normal"),
      "GGT",               sprintf("%.0f IU/L", s$GGT[1]),
                           ifelse(s$GGT_xULN[1]>1, "Elevated", "Normal"),
      "Fibrosis",          sprintf("F%.1f (Metavir)", s$Fibrosis[1]),
                           ifelse(s$Fibrosis[1]>=3, "Advanced", ifelse(s$Fibrosis[1]>=2, "Moderate", "Mild")),
      "AMA status",        "Positive (assumed)",  "≥1:40",
      "IgM",               sprintf("%.0f mg/dL", s$IgM[1]),
                           ifelse(input$IgM0>230,"Elevated","Normal")
    )
  })

  output$tbl_drugs <- DT::renderDataTable({
    DT::datatable(
      tribble(
        ~Drug,              ~Mechanism,        ~Dose,         ~Approval,  ~`ALP ↓`,     ~Trial,
        "UDCA",             "FXR partial / BSEP↑ / bicarbonate umbrella", "13-15 mg/kg/day", "SOC", "30-40%", "Lindor 1994",
        "Obeticholic acid", "FXR agonist (100×)",  "5-10 mg/day", "FDA 2016", "38%", "POISE",
        "Elafibranor",     "PPARα/δ dual",    "80 mg/day",    "FDA Jun 2024", "51% norm", "ELATIVE",
        "Seladelpar",      "PPARδ selective", "10 mg/day",    "FDA Aug 2024","25% norm", "RESPONSE",
        "Bezafibrate",     "pan-PPAR (α>δ>γ)","400 mg/day",  "Off-label", "67% norm", "BEZURSO",
        "Budesonide",      "GR agonist",      "9 mg/day",    "Off-label (non-cirrhotic)", "Variable", "Rautiainen 2005"
      ),
      options = list(dom = 't', pageLength = 10), rownames = FALSE
    )
  })

  # ── Tab 2: PK ─────────────────────────────────────────────
  output$plt_pk <- renderPlotly({
    s <- sim_result()
    pk <- s$pk_data

    fig <- plot_ly()
    if (input$use_UDCA && any(pk$UDCA > 0))
      fig <- fig %>% add_lines(data=pk, x=~time_h, y=~UDCA, name="UDCA (mg/L)", line=list(color="#1f77b4"))
    if (input$use_OCA && any(pk$OCA > 0))
      fig <- fig %>% add_lines(data=pk, x=~time_h, y=~OCA, name="OCA (mg/L)", line=list(color="#ff7f0e"), yaxis="y2")
    if (input$use_ELF && any(pk$ELF > 0))
      fig <- fig %>% add_lines(data=pk, x=~time_h, y=~ELF, name="Elafibranor (mg/L)", line=list(color="#9467bd"))
    if (input$use_SEL && any(pk$SEL > 0))
      fig <- fig %>% add_lines(data=pk, x=~time_h, y=~SEL, name="Seladelpar (mg/L)", line=list(color="#2ca02c"), yaxis="y2")
    if (input$use_BEZ && any(pk$BEZ > 0))
      fig <- fig %>% add_lines(data=pk, x=~time_h, y=~BEZ, name="Bezafibrate (mg/L)", line=list(color="#d62728"), yaxis="y3")

    fig %>% layout(
      title = "Drug Plasma Concentration Profiles (72h)",
      xaxis = list(title = "Time (hours)"),
      yaxis = list(title = "Concentration (mg/L)", side="left"),
      hovermode = "x unified"
    )
  })

  output$tbl_pk_summary <- DT::renderDataTable({
    s <- sim_result()
    DT::datatable(
      tribble(
        ~Drug,         ~`Css (mg/L)`, ~`T½ (h)`,  ~`F (%)`, ~`EC50 (mg/L)`, ~`Emax (ALP↓)`,
        "UDCA",        round(s$C_UDCA, 1), "24-48 (EHC)", "60%", "150",    "40%",
        "OCA",         round(s$C_OCA, 3),  "24-48",       "35%", "0.5",    "38%",
        "Elafibranor", round(s$C_ELF, 2),  "~20",         "72%", "0.8",    "55%",
        "Seladelpar",  round(s$C_SEL, 3),  "~14",         "60%", "0.3",    "30%",
        "Bezafibrate", round(s$C_BEZ, 1),  "~2.1",        "100%","2.0",    "67%"
      ),
      options = list(dom='t'), rownames=FALSE
    )
  })

  # ── Tab 3: PD Biomarkers ──────────────────────────────────
  make_gg_plotly <- function(x, y, ylab, thresh = NULL, thresh_label = NULL) {
    s <- sim_result()
    df <- data.frame(time = s[[x]], y = s[[y]])
    fig <- plot_ly(df, x=~time, y=~y, type="scatter", mode="lines",
                   line=list(color="#7b2d8b", width=2), name=ylab)
    if (!is.null(thresh)) {
      fig <- fig %>% add_segments(x=0, xend=max(df$time), y=thresh, yend=thresh,
                                  line=list(color="red",dash="dash",width=1.5),
                                  name=thresh_label, showlegend=TRUE)
    }
    fig %>% layout(xaxis=list(title="Time (months)"), yaxis=list(title=ylab))
  }

  output$plt_ALP      <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$ALP_xULN)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines",
            line=list(color="#d62728",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=1.67,yend=1.67,line=list(color="orange",dash="dash"),name="POISE threshold") %>%
      add_segments(x=0,xend=max(df$time),y=1.0,yend=1.0,line=list(color="red",dash="dot"),name="ULN") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="ALP (×ULN)"),hovermode="x")
  })
  output$plt_BILI     <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$BILI_xULN)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#ff7f0e",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=1,yend=1,line=list(color="red",dash="dash"),name="ULN") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Bilirubin (×ULN)"))
  })
  output$plt_GGT      <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$GGT_xULN)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#2196F3",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=1,yend=1,line=list(color="red",dash="dash"),name="ULN") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="GGT (×ULN)"))
  })
  output$plt_IgM      <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$IgM)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#e91e63",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=230,yend=230,line=list(color="red",dash="dash"),name="ULN 230 mg/dL") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="IgM (mg/dL)"))
  })
  output$plt_FGF19    <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$FGF19)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#4CAF50",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=100,yend=100,line=list(color="gray",dash="dash"),name="Baseline") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="FGF19 (pg/mL)"))
  })
  output$plt_pruritus <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$Pruritus)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#FF9800",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Pruritus NRS (0-10)", range=c(0,10)))
  })

  # ── Tab 4: Disease Biology ───────────────────────────────
  output$plt_AMA <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$AMA)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#E91E63",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="AMA (normalized)"))
  })
  output$plt_Th1 <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$Th1)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#FF9800",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Th1 (normalized)"))
  })
  output$plt_BEC <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$BEC_dmg)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#3F51B5",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="BEC Damage (0-1)",range=c(0,1)))
  })
  output$plt_BA <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$BA_toxic)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#009688",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Toxic BA pool (norm)"))
  })
  output$plt_FIB <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$Fibrosis)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#9C27B0",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=3,yend=3,line=list(color="red",dash="dash"),name="F3 (severe)") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Fibrosis (Metavir 0-4)",range=c(0,4)))
  })

  # ── Tab 5: Endpoints ──────────────────────────────────────
  output$vbox_paris2 <- renderValueBox({
    s <- sim_result()
    idx <- which.min(abs(s$time - 12))
    pct <- round(s$Paris2[idx] * 100)
    valueBox(paste0(pct, "%"),
             "Paris II Response at Month 12",
             icon=icon("check-circle"), color=if(pct>0)"green" else "red")
  })
  output$vbox_alp_norm <- renderValueBox({
    s <- sim_result()
    idx <- which.min(abs(s$time - 12))
    val <- round(s$ALP_xULN[idx], 2)
    valueBox(sprintf("%.2f×ULN", val),
             "ALP at Month 12",
             icon=icon("flask"), color=if(val<=1)"green" else if(val<=1.67)"yellow" else "red")
  })
  output$vbox_globe <- renderValueBox({
    s <- sim_result()
    idx <- which.min(abs(s$time - 24))
    val <- round(s$GLOBE_surrogate[idx], 2)
    valueBox(sprintf("%.2f", val),
             "GLOBE Score Surrogate (Month 24)",
             icon=icon("globe"), color=if(val<0)"green" else if(val<1)"yellow" else "red")
  })

  output$plt_paris2 <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$Paris2)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#4CAF50",width=2)) %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="Paris II Response (0/1)",range=c(-0.1,1.1)))
  })
  output$plt_globe <- renderPlotly({
    s <- sim_result()
    df <- data.frame(time=s$time, y=s$GLOBE_surrogate)
    plot_ly(df, x=~time, y=~y, type="scatter", mode="lines", line=list(color="#2196F3",width=2)) %>%
      add_segments(x=0,xend=max(df$time),y=0,yend=0,line=list(color="green",dash="dash"),name="Favorable boundary") %>%
      layout(xaxis=list(title="Time (months)"),yaxis=list(title="GLOBE Surrogate Score"))
  })
  output$tbl_m12 <- renderTable({
    s <- sim_result()
    idx <- which.min(abs(s$time - 12))
    tribble(
      ~Endpoint,              ~Value,                                    ~Interpretation,
      "ALP",                  sprintf("%.0f IU/L (%.1f×ULN)", s$ALP[idx], s$ALP_xULN[idx]),
                              ifelse(s$ALP_xULN[idx]<=1.0,"Normalized",ifelse(s$ALP_xULN[idx]<=1.67,"Partial response","No response")),
      "Bilirubin",            sprintf("%.2f mg/dL (%.1f×ULN)", s$BILI[idx], s$BILI_xULN[idx]),
                              ifelse(s$BILI_xULN[idx]<=1.0,"Normal","Elevated"),
      "GGT",                  sprintf("%.0f IU/L", s$GGT[idx]),
                              ifelse(s$GGT_xULN[idx]<=1,"Normal","Elevated"),
      "Paris II Response",    ifelse(s$Paris2[idx]==1,"YES","NO"),
                              "ALP<1.5×ULN + Bilirubin≤ULN",
      "ALP Normalization",    ifelse(s$ALP_norm[idx]==1,"YES","NO"),
                              "ALP≤ULN (ELATIVE/RESPONSE primary)",
      "Pruritus NRS",         sprintf("%.1f /10", s$Pruritus[idx]),
                              ifelse(s$Pruritus[idx]<4,"Mild",ifelse(s$Pruritus[idx]<7,"Moderate","Severe")),
      "Fibrosis",             sprintf("F%.1f (Metavir)", s$Fibrosis[idx]),
                              ifelse(s$Fibrosis[idx]<2,"Mild","Moderate/Advanced")
    )
  })

  # ── Tab 6: Scenario Comparison ────────────────────────────
  cols7 <- c(
    "No Treatment"       = "#7f7f7f",
    "UDCA Monotherapy"   = "#1f77b4",
    "UDCA + OCA 5mg"     = "#ff7f0e",
    "UDCA + OCA 10mg"    = "#ffbb78",
    "UDCA + Elafibranor" = "#9467bd",
    "UDCA + Seladelpar"  = "#2ca02c",
    "UDCA + Bezafibrate" = "#d62728"
  )

  output$plt_compare_ALP <- renderPlotly({
    df <- all_scenarios()
    fig <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df[df$scenario == sc, ]
      fig <- fig %>% add_lines(data=sub, x=~time, y=~ALP_xULN,
                               name=sc, line=list(width=2))
    }
    fig %>%
      add_segments(x=0, xend=max(df$time), y=1.67, yend=1.67,
                   line=list(color="red",dash="dash",width=1.5), name="POISE threshold", showlegend=TRUE) %>%
      add_segments(x=0, xend=max(df$time), y=1.0, yend=1.0,
                   line=list(color="darkred",dash="dot",width=1.5), name="ULN", showlegend=TRUE) %>%
      layout(title="ALP Response by Treatment (×ULN)",
             xaxis=list(title="Time (months)"),
             yaxis=list(title="ALP (×ULN)"),
             hovermode="x unified")
  })

  output$plt_compare_FIB <- renderPlotly({
    df <- all_scenarios()
    fig <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df[df$scenario == sc, ]
      fig <- fig %>% add_lines(data=sub, x=~time, y=~Fibrosis, name=sc, line=list(width=2))
    }
    fig %>% layout(title="Fibrosis by Treatment",
                   xaxis=list(title="Time (months)"),
                   yaxis=list(title="Metavir Score", range=c(0,4)))
  })

  output$plt_compare_pruritus <- renderPlotly({
    df <- all_scenarios()
    fig <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df[df$scenario == sc, ]
      fig <- fig %>% add_lines(data=sub, x=~time, y=~Pruritus, name=sc, line=list(width=2))
    }
    fig %>% layout(title="Pruritus NRS by Treatment (OCA worsen; PPARδ improve)",
                   xaxis=list(title="Time (months)"),
                   yaxis=list(title="Pruritus NRS (0-10)", range=c(0,10)))
  })

  output$tbl_compare <- DT::renderDataTable({
    df <- all_scenarios()
    m12 <- df %>%
      group_by(scenario) %>%
      filter(abs(time - 12) == min(abs(time - 12))) %>%
      slice(1) %>%
      ungroup() %>%
      select(
        Scenario = scenario,
        `ALP ×ULN (M12)` = ALP_xULN,
        `Bili ×ULN (M12)` = BILI_xULN,
        `GGT ×ULN (M12)` = GGT_xULN,
        `Paris II` = Paris2,
        `ALP Norm` = ALP_norm,
        `Pruritus NRS` = Pruritus
      ) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))

    DT::datatable(m12,
      options = list(dom = 'tpf', pageLength = 10, order = list(list(2, 'asc'))),
      rownames = FALSE
    )
  })
}

## ─────────────────────────────────────────────────────────────
## LAUNCH APP
## ─────────────────────────────────────────────────────────────

shinyApp(ui, server)
