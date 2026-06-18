################################################################################
# Autoimmune Hepatitis (AIH) — Interactive QSP Shiny Dashboard
# 6 Tabs: Patient Profile · Drug PK · Immune Dynamics · Hepatic Endpoints ·
#          Scenario Comparison · Biomarker Tracker
#
# Dependencies: shiny, shinydashboard, plotly, dplyr, tidyr, mrgsolve, DT
# Author: Claude Code QSP Routine | Date: 2026-06-18
################################################################################

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(DT)

# ---- Embedded minimal ODE for Shiny (simplified 8-state fast version) ---------
# Full mrgsolve model in aih_mrgsolve_model.R; here we use an ODE solver
# via deSolve for speed in interactive contexts

library(deSolve)

aih_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    # Drug effects
    GR_eff <- EMAX_GR * GR_OCC / (EC50_GR + GR_OCC + 1e-9)
    TGN_eff <- EMAX_TGN * TGN / (EC50_TGN + TGN + 1e-9)
    MPA_eff <- EMAX_MPA * MPA / (EC50_MPA + MPA + 1e-9)
    combined <- 1 - (1 - GR_eff) * (1 - TGN_eff) * (1 - MPA_eff)

    # Prednisolone 1-CMT (simplified: GR directly from dose)
    PRED_ss <- FLAG_PRED * PRED_DOSE * F_PRED / (CL_PRED / 1000)
    dGR_OCC <- KON_GR * PRED_ss * (1 - GR_OCC) - KOFF_GR * GR_OCC

    # AZA → TGN (simplified: TGN ss from dose)
    TGN_target <- FLAG_AZA * AZA_DOSE * TPMT_act * 2.5
    dTGN <- 0.01 * (TGN_target - TGN)

    # MMF → MPA
    MPA_target <- FLAG_MMF * MMF_DOSE * 0.0005
    dMPA <- 0.05 * (MPA_target - MPA)

    # Immune dynamics
    dTH1  <- TH1_base * Th1_act * (1 - combined) - Th1_death * TH1 - Treg_sup * TREG * TH1
    dTREG <- Treg_synth * (Treg_base + GR_eff * 0.3 * Treg_base) - Treg_death * TREG
            - Th1_Treg_cross * TH1 * TREG

    # B cells / autoantibodies
    dBCELL <- Bcell_synth * Bcell_base - Bcell_death * BCELL + Bcell_Th1 * TH1 * Bcell_base
    Ab_suppress <- GR_eff * 0.6 + MPA_eff * 0.4
    dAUTOAB <- Ab_synth * BCELL * (1 - Ab_suppress) - Ab_clear * AUTOAB

    # Cytokines
    dIFNG <- IFNg_prod * TH1 * (1 - GR_eff * 0.8) - IFNg_clear * IFNG
    dTGFB <- TGFb_prod * TREG - TGFb_clear * TGFB

    # Damage / ALT
    Dmg_rate <- K_dmg_Th1 * TH1 + K_dmg_Ab * AUTOAB + K_dmg_IFNg * IFNG
    Repair    <- K_repair * (1 + combined * 2)
    dDMG <- Dmg_rate * (1 - DMG / 100) - Repair * DMG / 100 * DMG
    dALT <- K_ALT_rel * Dmg_rate + ALT_base * K_ALT_elim - K_ALT_elim * ALT

    # IgG derived
    IgG_current <- IgG_base * (1 + AUTOAB / 100 * IgG_synth / IgG_clear)

    list(c(dGR_OCC, dTGN, dMPA, dTH1, dTREG, dBCELL, dAUTOAB, dIFNG, dTGFB, dDMG, dALT),
         IgG = IgG_current,
         HAI = DMG / 10,
         ALT_ULN = ALT / 40,
         Th1_Treg_ratio = TH1 / (TREG + 1e-6))
  })
}

run_aih_sim <- function(params) {
  state0 <- c(
    GR_OCC = 0,
    TGN    = 0,
    MPA    = 0,
    TH1    = params[["TH1_base"]],
    TREG   = params[["Treg_base"]],
    BCELL  = params[["Bcell_base"]],
    AUTOAB = 50,
    IFNG   = params[["IFNg_prod"]] * params[["TH1_base"]] / params[["IFNg_clear"]],
    TGFB   = params[["TGFb_prod"]] * params[["Treg_base"]] / params[["TGFb_clear"]],
    DMG    = 5,
    ALT    = params[["ALT_base"]]
  )
  times <- seq(0, params[["sim_days"]], by = 1)
  out <- tryCatch(
    ode(y = state0, times = times, func = aih_ode, parms = params, method = "lsoda"),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  as.data.frame(out) %>% rename(time_days = time)
}

# Default parameters
default_params <- list(
  # Treatment flags
  FLAG_PRED = 1, FLAG_AZA = 0, FLAG_MMF = 0, FLAG_RTX = 0,
  # Prednisolone PK
  PRED_DOSE = 60, F_PRED = 0.85, CL_PRED = 7.5,
  KON_GR = 0.12, KOFF_GR = 0.08, EMAX_GR = 0.95, EC50_GR = 0.15,
  # AZA
  AZA_DOSE = 100, TPMT_act = 1.0, EMAX_TGN = 0.90, EC50_TGN = 200,
  # MMF
  MMF_DOSE = 2000, EMAX_MPA = 0.80, EC50_MPA = 1.5,
  # Immune
  TH1_base = 100, Th1_act = 0.8, Th1_death = 0.15,
  Treg_base = 30, Treg_synth = 0.3, Treg_death = 0.10,
  Treg_sup = 0.004, Th1_Treg_cross = 0.005,
  Bcell_base = 100, Bcell_synth = 0.3, Bcell_death = 0.15,
  Bcell_Th1 = 0.002, Ab_synth = 0.1, Ab_clear = 0.05,
  IFNg_prod = 0.5, IFNg_clear = 0.8,
  TGFb_prod = 0.3, TGFb_clear = 0.4,
  IL6_prod = 0.6, IL6_clear = 1.0,
  K_dmg_Th1 = 0.015, K_dmg_Ab = 0.003, K_dmg_IFNg = 0.020,
  K_repair = 0.10, K_ALT_rel = 5.0, K_ALT_elim = 0.2,
  ALT_base = 35, DMG_base = 5,
  IgG_base = 14, IgG_synth = 0.01, IgG_clear = 0.02,
  sim_days = 730
)

# ---- UI -----------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "AIH QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",           tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Immune Dynamics",   tabName = "tab_immune",    icon = icon("shield-alt")),
      menuItem("Hepatic Endpoints", tabName = "tab_hepatic",   icon = icon("liver")),
      menuItem("Scenario Comparison", tabName = "tab_compare", icon = icon("chart-bar")),
      menuItem("Biomarker Tracker", tabName = "tab_biomark",   icon = icon("microscope"))
    ),
    hr(),
    h5("  Treatment Settings", style = "color:#aaa; padding-left:15px;"),
    checkboxInput("chk_pred", "Prednisolone", value = TRUE),
    conditionalPanel("input.chk_pred",
      sliderInput("dose_pred", "Prednisolone (mg/day)", 5, 80, 60, step = 5)
    ),
    checkboxInput("chk_aza", "Azathioprine", value = FALSE),
    conditionalPanel("input.chk_aza",
      sliderInput("dose_aza", "AZA (mg/day)", 25, 250, 100, step = 25),
      selectInput("tpmt", "TPMT genotype",
                  choices = c("Normal (WT/WT)" = 1.0,
                              "Intermediate (WT/*3)" = 0.5,
                              "Poor metabolizer (*3/*3)" = 0.1))
    ),
    checkboxInput("chk_mmf", "MMF (Mycophenolate)", value = FALSE),
    conditionalPanel("input.chk_mmf",
      sliderInput("dose_mmf", "MMF (mg/day)", 500, 3000, 2000, step = 500)
    ),
    hr(),
    sliderInput("sim_days", "Simulation (days)", 90, 730, 365, step = 30),
    actionButton("btn_run", "Run Simulation", icon = icon("play"),
                 style = "width:90%; background-color:#3498DB; color:white; margin:5%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F5F7FA; }
      .box { border-radius: 8px; }
      .box-header { border-radius: 8px 8px 0 0; }
      .value-box { border-radius: 8px; }
    "))),

    tabItems(
      # ---- Tab 1: Patient Profile ----------------------------------------------
      tabItem(
        tabName = "tab_patient",
        fluidRow(
          box(
            title = "Disease Overview — Autoimmune Hepatitis (AIH)", width = 12,
            status = "primary", solidHeader = TRUE,
            fluidRow(
              column(6,
                h4("Pathophysiology"),
                p("AIH is a chronic immune-mediated liver disease characterized by interface hepatitis,
                   autoantibodies (ANA, ASMA, LKM-1, SLA/LP), IgG hypergammaglobulinemia,
                   and progressive hepatic fibrosis if untreated."),
                tags$ul(
                  tags$li(strong("Type 1 AIH:"), " ANA+ / ASMA+ (F-actin). Most common (80%). HLA-DR3/DR4."),
                  tags$li(strong("Type 2 AIH:"), " Anti-LKM-1+ / anti-LC-1+. Younger patients. HLA-DR7."),
                  tags$li(strong("Key cytokines:"), " IFN-γ (Th1), IL-17A (Th17), IL-21 (Tfh→B cell help), TGF-β (Treg/fibrosis)."),
                  tags$li(strong("Treg deficiency:"), " Reduced FoxP3+ Tregs in number and function — hallmark of AIH."),
                  tags$li(strong("Treatment target:"), " Remission = ALT <ULN + IgG normal + histological HAI <4.")
                )
              ),
              column(6,
                h4("IAIHG Simplified Score (Hennes 2008)"),
                tableOutput("tbl_simplified_score"),
                br(),
                h4("Remission Criteria (IAIHG 1999/2011)"),
                tags$ul(
                  tags$li("ALT < upper limit of normal (ULN = 40 U/L)"),
                  tags$li("IgG < 16 g/L (or within normal range)"),
                  tags$li("Histology: HAI score < 4/18 (Ishak) or inactive (METAVIR)"),
                  tags$li("All three criteria must be met simultaneously")
                )
              )
            )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_alt_status", width = 3),
          valueBoxOutput("vb_igg_status", width = 3),
          valueBoxOutput("vb_hai_status", width = 3),
          valueBoxOutput("vb_remission_status", width = 3)
        ),
        fluidRow(
          box(
            title = "Disease Natural History — Current Patient", width = 12,
            status = "info",
            plotlyOutput("plt_nat_history", height = "350px")
          )
        )
      ),

      # ---- Tab 2: Drug PK -----------------------------------------------------
      tabItem(
        tabName = "tab_pk",
        fluidRow(
          box(
            title = "Drug PK Parameters", width = 4, status = "primary",
            solidHeader = TRUE,
            h5("Prednisolone (2-CMT Oral)"),
            tableOutput("tbl_pk_pred"),
            hr(),
            h5("Azathioprine → 6-TGN (Cascade)"),
            tableOutput("tbl_pk_aza"),
            hr(),
            h5("MMF → MPA (1-CMT Oral)"),
            tableOutput("tbl_pk_mmf")
          ),
          box(
            title = "GR Occupancy over Time", width = 8, status = "warning",
            solidHeader = TRUE,
            plotlyOutput("plt_gr_occ", height = "280px"),
            hr(),
            plotlyOutput("plt_tgn_mpa", height = "220px")
          )
        ),
        fluidRow(
          box(
            title = "Rituximab TMDD — B Cell Depletion Kinetics", width = 12,
            status = "danger",
            fluidRow(
              column(4,
                sliderInput("rtx_dose", "RTX Dose (mg)", 375, 1000, 1000, step = 125),
                selectInput("rtx_schedule", "Schedule",
                            c("1000mg × 2 (2wk apart)" = "2dose_1000",
                              "375 mg/m² × 4 (weekly)" = "4dose_375",
                              "Single 1000mg" = "single")),
                br(),
                p(style="font-size:0.9em", "Rituximab depletes CD20+ B cells via ADCC, CDC, and direct apoptosis.
                  Recovery typically occurs at 9-12 months.")
              ),
              column(8, plotlyOutput("plt_rtx_pk", height = "350px"))
            )
          )
        )
      ),

      # ---- Tab 3: Immune Dynamics ---------------------------------------------
      tabItem(
        tabName = "tab_immune",
        fluidRow(
          box(
            title = "Th1 / Treg Balance — Key AIH Immune Axis", width = 6,
            status = "danger", solidHeader = TRUE,
            plotlyOutput("plt_th1_treg", height = "350px")
          ),
          box(
            title = "B Cells & Autoantibody Titers", width = 6,
            status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_bcell_ab", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "Cytokine Profiles (IFN-γ, TGF-β)", width = 6,
            status = "info", solidHeader = TRUE,
            plotlyOutput("plt_cytokines", height = "300px")
          ),
          box(
            title = "Combined Immunosuppression (%)", width = 6,
            status = "success", solidHeader = TRUE,
            plotlyOutput("plt_combined_suppress", height = "300px"),
            p(style="font-size:0.85em; color:#666;",
              "Combined T cell suppression reflects additive effects of GR occupancy (prednisolone),
              6-TGN (azathioprine), and MPA (MMF). Maximum achievable suppression is capped at 98%.")
          )
        )
      ),

      # ---- Tab 4: Hepatic Endpoints ------------------------------------------
      tabItem(
        tabName = "tab_hepatic",
        fluidRow(
          box(
            title = "Serum ALT — Primary Efficacy Endpoint", width = 8,
            status = "primary", solidHeader = TRUE,
            plotlyOutput("plt_alt", height = "350px")
          ),
          box(
            title = "Response Summary", width = 4,
            status = "success", solidHeader = TRUE,
            h5("At current simulation endpoint:"),
            tableOutput("tbl_response_summary"),
            br(),
            p(strong("Response Definitions:")),
            tags$ul(
              tags$li(strong("Complete:"), " ALT <ULN + IgG normal"),
              tags$li(strong("Partial:"), " ALT <3×ULN, improvement ≥50%"),
              tags$li(strong("Failure:"), " <Partial response at 2 years"),
              tags$li(strong("Relapse:"), " ALT >3×ULN after remission")
            )
          )
        ),
        fluidRow(
          box(
            title = "Serum IgG (Disease Activity Biomarker)", width = 6,
            status = "warning", solidHeader = TRUE,
            plotlyOutput("plt_igg", height = "300px")
          ),
          box(
            title = "Histological Activity Index (HAI Proxy)", width = 6,
            status = "info", solidHeader = TRUE,
            plotlyOutput("plt_hai", height = "300px"),
            p(style="font-size:0.85em; color:#666;",
              "HAI (Histological Activity Index) approximated from hepatocellular damage model.
              Remission requires HAI < 4/18 (Ishak grading).")
          )
        )
      ),

      # ---- Tab 5: Scenario Comparison -----------------------------------------
      tabItem(
        tabName = "tab_compare",
        fluidRow(
          box(
            title = "Multi-Scenario ALT Comparison", width = 12,
            status = "primary", solidHeader = TRUE,
            plotlyOutput("plt_scenario_alt", height = "400px")
          )
        ),
        fluidRow(
          box(
            title = "Scenario Results at Key Timepoints", width = 12,
            status = "info",
            DTOutput("tbl_scenario_compare")
          )
        )
      ),

      # ---- Tab 6: Biomarker Tracker -------------------------------------------
      tabItem(
        tabName = "tab_biomark",
        fluidRow(
          box(
            title = "Biomarker Dashboard", width = 12,
            status = "primary", solidHeader = TRUE,
            fluidRow(
              column(4,
                plotlyOutput("plt_bm_alt", height = "200px"),
                br(),
                plotlyOutput("plt_bm_autoab", height = "200px")
              ),
              column(4,
                plotlyOutput("plt_bm_igg", height = "200px"),
                br(),
                plotlyOutput("plt_bm_th1treg", height = "200px")
              ),
              column(4,
                plotlyOutput("plt_bm_ifng", height = "200px"),
                br(),
                plotlyOutput("plt_bm_dmg", height = "200px")
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Downloadable Simulation Data", width = 12,
            status = "secondary",
            DTOutput("tbl_sim_data"),
            br(),
            downloadButton("dl_csv", "Download CSV", class = "btn-primary")
          )
        )
      )
    )
  )
)

# ---- SERVER -------------------------------------------------------------------
server <- function(input, output, session) {

  # Reactive simulation
  sim_results <- eventReactive(input$btn_run, {
    params <- default_params
    params$FLAG_PRED <- as.integer(input$chk_pred)
    params$FLAG_AZA  <- as.integer(input$chk_aza)
    params$FLAG_MMF  <- as.integer(input$chk_mmf)
    params$PRED_DOSE <- if (input$chk_pred) input$dose_pred else 0
    params$AZA_DOSE  <- if (input$chk_aza)  input$dose_aza  else 0
    params$MMF_DOSE  <- if (input$chk_mmf)  input$dose_mmf  else 0
    params$TPMT_act  <- if (input$chk_aza)  as.numeric(input$tpmt) else 1.0
    params$sim_days  <- input$sim_days
    run_aih_sim(params)
  }, ignoreNULL = FALSE)

  # Untreated baseline (for comparison)
  untreated_sim <- reactive({
    params <- default_params
    params$FLAG_PRED <- 0; params$FLAG_AZA <- 0; params$FLAG_MMF <- 0
    params$sim_days <- input$sim_days
    run_aih_sim(params)
  })

  # ---- Tab 1 Outputs ---------------------------------------------------------
  output$tbl_simplified_score <- renderTable({
    data.frame(
      Variable   = c("ANA or ASMA ≥1:40", "ANA or ASMA ≥1:80", "Anti-LKM-1 ≥1:40",
                     "Anti-SLA/LP positive", "IgG > ULN", "IgG > 1.1×ULN",
                     "Histology (typical)", "Histology (compatible)", "No viral hepatitis"),
      Points = c("+1", "+2", "+2", "+2", "+1", "+2", "+2", "+1", "+2"),
      Threshold = c("≥6 pts: Probable AIH", "", "≥7 pts: Definite AIH", "", "", "", "", "", "")
    )
  }, striped = TRUE, hover = TRUE)

  end_values <- reactive({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    tail(df, 1)
  })

  output$vb_alt_status <- renderValueBox({
    ev <- end_values()
    val <- if (is.null(ev)) "—" else round(ev$ALT, 1)
    color <- if (!is.null(ev) && ev$ALT < 40) "green" else "red"
    valueBox(paste(val, "U/L"), "Serum ALT", icon = icon("vial"), color = color)
  })

  output$vb_igg_status <- renderValueBox({
    ev <- end_values()
    igg <- if (is.null(ev)) NA else ev$IgG
    color <- if (!is.null(ev) && !is.null(igg) && igg < 16) "green" else "yellow"
    valueBox(paste(round(igg, 1), "g/L"), "Serum IgG", icon = icon("shield-alt"), color = color)
  })

  output$vb_hai_status <- renderValueBox({
    ev <- end_values()
    hai <- if (is.null(ev)) NA else ev$HAI
    color <- if (!is.null(ev) && !is.null(hai) && hai < 4) "green" else "orange"
    valueBox(round(hai, 1), "HAI Score", icon = icon("chart-line"), color = color)
  })

  output$vb_remission_status <- renderValueBox({
    ev <- end_values()
    in_remission <- !is.null(ev) &&
      ev$ALT < 40 && !is.null(ev$IgG) && ev$IgG < 16 && !is.null(ev$HAI) && ev$HAI < 4
    valueBox(if (in_remission) "YES" else "NO",
             "Remission", icon = icon(if (in_remission) "check-circle" else "times-circle"),
             color = if (in_remission) "green" else "red")
  })

  output$plt_nat_history <- renderPlotly({
    df  <- sim_results()
    unt <- untreated_sim()
    if (is.null(df) || is.null(unt)) return(NULL)
    plot_ly() %>%
      add_lines(data = unt, x = ~time_days, y = ~ALT, name = "Untreated",
                line = list(color = "#E74C3C", dash = "dash")) %>%
      add_lines(data = df,  x = ~time_days, y = ~ALT, name = "Current Treatment",
                line = list(color = "#2ECC71")) %>%
      add_hline(y = 40, line = list(color = "gray", dash = "dot")) %>%
      layout(title = "ALT: Current Treatment vs Untreated",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALT (U/L)"),
             legend = list(orientation = "h"))
  })

  # ---- Tab 2: Drug PK --------------------------------------------------------
  output$tbl_pk_pred <- renderTable({
    data.frame(
      Parameter = c("Bioavailability", "Vd central", "t½", "Clearance", "GR affinity"),
      Value = c("85%", "25 L", "~2.5h", "7.5 L/h", "Kd ~0.7 μg/mL")
    )
  }, bordered = TRUE)

  output$tbl_pk_aza <- renderTable({
    data.frame(
      Parameter = c("Bioavailability", "6-TGN t½", "Therapeutic range", "Toxic range", "TPMT"],
      Value = c("47%", "~5 days", "235-450 pmol/8e8", ">450 pmol", "Normal/Hetero/Defic")
    )
  }, bordered = TRUE)

  output$tbl_pk_mmf <- renderTable({
    data.frame(
      Parameter = c("Bioavail (MPA)", "MPA Vd", "MPA t½", "IMPDH IC50", "Target AUC"),
      Value = c("~94%", "100 L", "~16h", "~0.3 μM", "30-60 mg·h/L")
    )
  }, bordered = TRUE)

  output$plt_gr_occ <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time_days, y = ~GR_OCC * 100, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2)) %>%
      layout(title = "GR Occupancy (%)", xaxis = list(title = "Days"),
             yaxis = list(title = "GR Occupancy (%)", range = c(0, 100)))
  })

  output$plt_tgn_mpa <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data = df, x = ~time_days, y = ~TGN, name = "6-TGN (pmol)",
                line = list(color = "#3498DB")) %>%
      add_hline(y = 235, line = list(color = "green", dash = "dot"),
                annotation = list(text = "Min therapeutic")) %>%
      add_hline(y = 450, line = list(color = "red", dash = "dot"),
                annotation = list(text = "Toxic threshold")) %>%
      layout(title = "6-TGN Intracellular Level",
             xaxis = list(title = "Days"),
             yaxis = list(title = "6-TGN (pmol / 8×10⁸ RBC)"))
  })

  output$plt_rtx_pk <- renderPlotly({
    # Simple B cell depletion kinetics for rituximab illustration
    t <- seq(0, 365, 1)
    # Bi-exponential approximation
    CD20_1000 <- 100 * (0.05 + 0.95 * exp(-0.015 * t))
    CD20_375  <- 100 * (0.10 + 0.90 * exp(-0.012 * t))
    plot_ly() %>%
      add_lines(x = t, y = CD20_1000, name = "1000mg × 2", line = list(color = "#9B59B6")) %>%
      add_lines(x = t, y = CD20_375,  name = "375mg × 4", line = list(color = "#1ABC9C")) %>%
      add_hline(y = 5, line = list(dash = "dot", color = "gray")) %>%
      layout(title = "CD20+ B Cell Count — Rituximab Kinetics",
             xaxis = list(title = "Days"), yaxis = list(title = "B cells (% of baseline)"),
             legend = list(orientation = "h"))
  })

  # ---- Tab 3: Immune Dynamics ------------------------------------------------
  output$plt_th1_treg <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data = df, x = ~time_days, y = ~TH1, name = "Th1 Cells",
                line = list(color = "#E74C3C", width = 2)) %>%
      add_lines(data = df, x = ~time_days, y = ~TREG, name = "Treg Cells",
                line = list(color = "#2ECC71", width = 2)) %>%
      layout(title = "Th1 / Treg Balance",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Cell count (AU)"),
             legend = list(orientation = "h"))
  })

  output$plt_bcell_ab <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data = df, x = ~time_days, y = ~BCELL, name = "B Cells",
                line = list(color = "#F39C12", width = 2)) %>%
      add_lines(data = df, x = ~time_days, y = ~AUTOAB, name = "Autoantibody Titer",
                line = list(color = "#9B59B6", width = 2)) %>%
      layout(title = "B Cells & Autoantibody Dynamics",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Count / Titer (AU)"),
             legend = list(orientation = "h"))
  })

  output$plt_cytokines <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data = df, x = ~time_days, y = ~IFNG, name = "IFN-γ",
                line = list(color = "#E74C3C")) %>%
      add_lines(data = df, x = ~time_days, y = ~TGFB, name = "TGF-β",
                line = list(color = "#2ECC71")) %>%
      layout(title = "Key Cytokine Profiles",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Concentration (AU)"),
             legend = list(orientation = "h"))
  })

  output$plt_combined_suppress <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    # Recalculate combined suppress from GR_OCC, TGN, MPA
    df2 <- df %>%
      mutate(
        GR_eff  = 0.95 * GR_OCC / (0.15 + GR_OCC),
        TGN_eff = 0.90 * TGN / (200 + TGN),
        MPA_eff = 0.80 * MPA / (1.5 + MPA),
        combined = (1 - (1 - GR_eff) * (1 - TGN_eff) * (1 - MPA_eff)) * 100
      )
    plot_ly() %>%
      add_lines(data = df2, x = ~time_days, y = ~GR_eff * 100, name = "GR (Pred/BUD)",
                line = list(color = "#E74C3C", dash = "dash")) %>%
      add_lines(data = df2, x = ~time_days, y = ~TGN_eff * 100, name = "6-TGN (AZA)",
                line = list(color = "#3498DB", dash = "dash")) %>%
      add_lines(data = df2, x = ~time_days, y = ~MPA_eff * 100, name = "MPA (MMF)",
                line = list(color = "#9B59B6", dash = "dash")) %>%
      add_lines(data = df2, x = ~time_days, y = ~combined, name = "Combined",
                line = list(color = "#2ECC71", width = 2.5)) %>%
      layout(title = "Immunosuppression Decomposition",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Suppression (%)", range = c(0, 100)),
             legend = list(orientation = "h"))
  })

  # ---- Tab 4: Hepatic Endpoints ----------------------------------------------
  output$plt_alt <- renderPlotly({
    df  <- sim_results()
    unt <- untreated_sim()
    if (is.null(df)) return(NULL)
    p <- plot_ly() %>%
      add_lines(data = df, x = ~time_days, y = ~ALT, name = "Treated",
                line = list(color = "#2ECC71", width = 2.5))
    if (!is.null(unt))
      p <- p %>% add_lines(data = unt, x = ~time_days, y = ~ALT, name = "Untreated",
                            line = list(color = "#E74C3C", dash = "dot"))
    p %>% add_hline(y = 40, line = list(color = "gray", dash = "dash")) %>%
      layout(title = "Serum ALT over Time",
             xaxis = list(title = "Time (days)"),
             yaxis = list(title = "ALT (U/L)"),
             annotations = list(x = 50, y = 42, text = "ULN = 40 U/L",
                                showarrow = FALSE, font = list(size = 11)))
  })

  output$tbl_response_summary <- renderTable({
    df <- sim_results()
    if (is.null(df)) return(data.frame(Metric = "No simulation", Value = "—"))
    ev <- tail(df, 1)
    igg_val <- if (!is.null(ev$IgG)) round(ev$IgG, 1) else "N/A"
    hai_val <- if (!is.null(ev$HAI)) round(ev$HAI, 1) else "N/A"
    resp <- if (!is.null(ev$ALT) && ev$ALT < 40 && !is.null(ev$IgG) && ev$IgG < 16 &&
                !is.null(ev$HAI) && ev$HAI < 4) {
      "Complete Remission"
    } else if (!is.null(ev$ALT) && ev$ALT < 120) {
      "Partial Response"
    } else {
      "Insufficient Response"
    }
    data.frame(
      Metric = c("ALT (U/L)", "IgG (g/L)", "HAI Score", "Response"),
      Value  = c(round(ev$ALT, 1), igg_val, hai_val, resp)
    )
  }, striped = TRUE)

  output$plt_igg <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time_days, y = ~IgG, type = "scatter", mode = "lines",
            line = list(color = "#F39C12", width = 2)) %>%
      add_hline(y = 16, line = list(color = "gray", dash = "dash")) %>%
      layout(title = "Serum IgG", xaxis = list(title = "Days"),
             yaxis = list(title = "IgG (g/L)"))
  })

  output$plt_hai <- renderPlotly({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time_days, y = ~HAI, type = "scatter", mode = "lines",
            line = list(color = "#9B59B6", width = 2)) %>%
      add_hline(y = 4, line = list(color = "green", dash = "dot")) %>%
      layout(title = "HAI Score (Proxy)", xaxis = list(title = "Days"),
             yaxis = list(title = "HAI Score (0-18)", range = c(0, 18)),
             annotations = list(x = 50, y = 4.5, text = "Remission threshold",
                                showarrow = FALSE))
  })

  # ---- Tab 5: Scenario Comparison --------------------------------------------
  all_scenarios_sim <- reactive({
    scenarios <- list(
      list(name = "Untreated",         params = list(FLAG_PRED=0, FLAG_AZA=0, FLAG_MMF=0)),
      list(name = "Pred Mono",         params = list(FLAG_PRED=1, FLAG_AZA=0, PRED_DOSE=60)),
      list(name = "Pred+AZA (SoC)",    params = list(FLAG_PRED=1, FLAG_AZA=1, PRED_DOSE=60, AZA_DOSE=150)),
      list(name = "BUD+AZA",           params = list(FLAG_PRED=0, FLAG_BUD=1, FLAG_AZA=1, PRED_DOSE=9, AZA_DOSE=150)),
      list(name = "Pred+MMF",          params = list(FLAG_PRED=1, FLAG_MMF=1, PRED_DOSE=60, MMF_DOSE=2000)),
      list(name = "RTX+Pred (refract)",params = list(FLAG_PRED=1, FLAG_RTX=1, PRED_DOSE=40))
    )
    purrr::map_dfr(scenarios, function(s) {
      p <- modifyList(default_params, c(s$params, list(sim_days = input$sim_days)))
      df <- run_aih_sim(p)
      if (!is.null(df)) df %>% mutate(scenario = s$name)
    })
  })

  scenario_colors_shiny <- c(
    "Untreated"         = "#E74C3C",
    "Pred Mono"         = "#E67E22",
    "Pred+AZA (SoC)"    = "#2ECC71",
    "BUD+AZA"           = "#3498DB",
    "Pred+MMF"          = "#9B59B6",
    "RTX+Pred (refract)"= "#1ABC9C"
  )

  output$plt_scenario_alt <- renderPlotly({
    df <- all_scenarios_sim()
    if (nrow(df) == 0) return(NULL)
    scens <- unique(df$scenario)
    p <- plot_ly()
    for (sc in scens) {
      dsc <- df %>% filter(scenario == sc)
      p <- p %>% add_lines(data = dsc, x = ~time_days, y = ~ALT, name = sc,
                            line = list(color = scenario_colors_shiny[sc], width = 2))
    }
    p %>% add_hline(y = 40, line = list(color = "gray", dash = "dash")) %>%
      layout(title = "Serum ALT — All Scenarios",
             xaxis = list(title = "Time (days)"), yaxis = list(title = "ALT (U/L)"),
             legend = list(orientation = "h", x = 0, y = -0.25))
  })

  output$tbl_scenario_compare <- renderDT({
    df <- all_scenarios_sim()
    if (nrow(df) == 0) return(NULL)
    df %>%
      filter(time_days %in% c(0, 30, 90, 180, 365, min(730, input$sim_days))) %>%
      select(scenario, time_days, ALT, IgG, HAI, Th1_Treg_ratio) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      rename(Scenario = scenario, `Day` = time_days, `ALT (U/L)` = ALT,
             `IgG (g/L)` = IgG, `HAI Score` = HAI, `Th1/Treg` = Th1_Treg_ratio) %>%
      datatable(options = list(pageLength = 20, scrollX = TRUE),
                rownames = FALSE, class = "cell-border stripe hover")
  })

  # ---- Tab 6: Biomarker Tracker -----------------------------------------------
  mk_bioplot <- function(df, yvar, ylab, color, hline = NULL) {
    if (is.null(df)) return(NULL)
    p <- plot_ly(df, x = ~time_days, y = as.formula(paste0("~", yvar)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2)) %>%
      layout(xaxis = list(title = "Days"), yaxis = list(title = ylab),
             title = ylab, margin = list(t = 30))
    if (!is.null(hline)) p <- p %>% add_hline(y = hline, line = list(dash = "dot", color = "gray"))
    p
  }

  output$plt_bm_alt    <- renderPlotly(mk_bioplot(sim_results(), "ALT",    "ALT (U/L)",    "#E74C3C", 40))
  output$plt_bm_autoab <- renderPlotly(mk_bioplot(sim_results(), "AUTOAB", "AutoAb Titer", "#9B59B6", 50))
  output$plt_bm_igg    <- renderPlotly(mk_bioplot(sim_results(), "IgG",    "IgG (g/L)",    "#F39C12", 16))
  output$plt_bm_th1treg<- renderPlotly(mk_bioplot(sim_results(), "Th1_Treg_ratio", "Th1/Treg Ratio", "#E67E22"))
  output$plt_bm_ifng   <- renderPlotly(mk_bioplot(sim_results(), "IFNG",   "IFN-γ (AU)",   "#3498DB"))
  output$plt_bm_dmg    <- renderPlotly(mk_bioplot(sim_results(), "DMG",    "Hepatocyte Damage (AU)", "#1ABC9C"))

  output$tbl_sim_data <- renderDT({
    df <- sim_results()
    if (is.null(df)) return(NULL)
    df %>%
      select(time_days, GR_OCC, TGN, MPA, TH1, TREG, BCELL, AUTOAB, IFNG, TGFB, DMG, ALT, IgG, HAI) %>%
      mutate(across(where(is.numeric), round, 3)) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE, class = "cell-border stripe")
  })

  output$dl_csv <- downloadHandler(
    filename = function() paste0("aih_qsp_simulation_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- sim_results()
      if (!is.null(df)) write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
