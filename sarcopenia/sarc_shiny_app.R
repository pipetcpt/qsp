################################################################################
##  Sarcopenia QSP — Shiny Dashboard
##  Eight tabs: Patient profile · Drug PK · Anabolic/Catabolic PD ·
##              Muscle mass & function · Clinical endpoints (EWGSOP2) ·
##              Scenario comparison · Biomarkers · References
##
##  Dependencies:
##    install.packages(c("shiny","shinydashboard","mrgsolve","tidyverse",
##                       "plotly","DT","markdown"))
##
##  Launch:  shiny::runApp("sarc_shiny_app.R")
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(tidyverse)
library(plotly)
library(DT)

#### -- Model loader (graceful fallback if rendering fails) -- ###################
load_model <- function() {
  tryCatch(
    mrgsolve::mread(
      "sarc_mrgsolve_model.R",
      project = system.file(package = "mrgsolve") %>% dirname() %>%
        file.path("..", "..", "qsp", "sarcopenia") %>% normalizePath(mustWork = FALSE)
    ),
    error = function(e) NULL
  )
}
mod <- tryCatch(mrgsolve::mread("sarc_mrgsolve_model.R", project = "."), error = function(e) NULL)

#### -- UI -- ####################################################################
ui <- dashboardPage(
  dashboardHeader(title = "Sarcopenia QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient profile",   tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",           tabName = "pk",      icon = icon("pills")),
      menuItem("3. Anabolic / Catabolic PD", tabName = "pd", icon = icon("balance-scale")),
      menuItem("4. Muscle mass & function",  tabName = "muscle", icon = icon("dumbbell")),
      menuItem("5. Clinical endpoints (EWGSOP2)", tabName = "clinical", icon = icon("stethoscope")),
      menuItem("6. Scenario comparison", tabName = "scenarios", icon = icon("layer-group")),
      menuItem("7. Biomarkers",        tabName = "bio",     icon = icon("flask")),
      menuItem("8. References",        tabName = "ref",     icon = icon("book"))
    )
  ),
  dashboardBody(
    tabItems(
      # 1. Patient profile -------------------------------------------------------
      tabItem(tabName = "patient",
              fluidRow(
                box(width = 4, title = "Demographics",
                    numericInput("age", "Age (yr)", 75, min = 40, max = 100),
                    selectInput("sex", "Sex", c("Male" = 0, "Female" = 1)),
                    numericInput("wt",  "Weight (kg)", 70),
                    numericInput("ht",  "Height (m)",  1.70, step = 0.01),
                    sliderInput("act", "Activity level (0 = bedrest, 1 = active)",
                                0, 1, 0.5, 0.1)),
                box(width = 4, title = "Baseline phenotype",
                    numericInput("alm",  "Baseline ALM (kg)", 18),
                    numericInput("grip", "Grip (kg)", 25),
                    numericInput("gait", "Gait speed (m/s)", 0.95, step = 0.05),
                    numericInput("sppb", "SPPB (0-12)", 9, min = 0, max = 12)),
                box(width = 4, title = "Endocrine / inflammation",
                    numericInput("igf1", "Serum IGF-1 (ng/mL)", 130),
                    numericInput("vitd", "25(OH)D (ng/mL)", 22),
                    numericInput("mstn", "Serum myostatin (ng/mL)", 6),
                    numericInput("il6",  "IL-6 (pg/mL)", 3.5),
                    numericInput("gdf",  "GDF-15 (pg/mL)", 1200))
              ),
              fluidRow(
                box(width = 12, title = "Computed sarcopenia status (EWGSOP2)",
                    verbatimTextOutput("dx_status"))
              )
      ),
      # 2. Drug PK ---------------------------------------------------------------
      tabItem(tabName = "pk",
              fluidRow(
                box(width = 4, title = "Bimagrumab",
                    numericInput("bima_dose", "Dose (mg/kg)", 30),
                    numericInput("bima_ii",   "Interval (days)", 28),
                    numericInput("bima_n",    "Doses", 12)),
                box(width = 4, title = "Apitegromab",
                    numericInput("apit_dose", "Dose (mg/kg)", 20),
                    numericInput("apit_ii",   "Interval (days)", 28),
                    numericInput("apit_n",    "Doses", 12)),
                box(width = 4, title = "Testosterone / Vit D / EAA / RT",
                    numericInput("test_dose", "Testosterone IM (mg/wk)", 0),
                    numericInput("vitd_dose", "Vitamin D₃ (IU/day)", 0),
                    numericInput("leu_dose",  "Leucine (mmol/meal)", 0),
                    numericInput("rt_sess",   "Resistance sessions/wk", 0))
              ),
              fluidRow(
                box(width = 12, title = "Drug concentrations vs time",
                    plotlyOutput("pk_plot", height = 480))
              )
      ),
      # 3. PD: anabolic/catabolic ------------------------------------------------
      tabItem(tabName = "pd",
              fluidRow(
                box(width = 6, plotlyOutput("anab_plot")),
                box(width = 6, plotlyOutput("catab_plot"))
              ),
              fluidRow(
                box(width = 12, title = "IGF-1, myostatin, IL-6, GDF-15 trajectories",
                    plotlyOutput("endo_plot"))
              )
      ),
      # 4. Muscle mass & function ------------------------------------------------
      tabItem(tabName = "muscle",
              fluidRow(
                box(width = 6, plotlyOutput("alm_plot")),
                box(width = 6, plotlyOutput("grip_plot"))
              ),
              fluidRow(
                box(width = 6, plotlyOutput("gait_plot")),
                box(width = 6, plotlyOutput("sppb_plot"))
              )
      ),
      # 5. Clinical endpoints ----------------------------------------------------
      tabItem(tabName = "clinical",
              fluidRow(
                box(width = 6, title = "ALMI vs EWGSOP2 cutoff",
                    plotlyOutput("almi_plot")),
                box(width = 6, title = "Falls cumulative & frailty",
                    plotlyOutput("fall_plot"))
              ),
              fluidRow(
                box(width = 12, title = "Diagnostic flag trajectory",
                    DTOutput("dx_table"))
              )
      ),
      # 6. Scenario comparison ---------------------------------------------------
      tabItem(tabName = "scenarios",
              fluidRow(
                box(width = 12, title = "Six treatment arms — ALM and gait speed at 1 year",
                    plotlyOutput("scen_alm", height = 360),
                    plotlyOutput("scen_gait", height = 360))
              ),
              fluidRow(
                box(width = 12, title = "Summary table",
                    DTOutput("scen_table"))
              )
      ),
      # 7. Biomarkers ------------------------------------------------------------
      tabItem(tabName = "bio",
              fluidRow(
                box(width = 12,
                    title = "Biomarker panel: IGF-1, MSTN, IL-6, GDF-15, 25(OH)D, CAF",
                    plotlyOutput("bio_panel", height = 520))
              )
      ),
      # 8. References ------------------------------------------------------------
      tabItem(tabName = "ref",
              fluidRow(
                box(width = 12, title = "Curated references (PubMed)",
                    includeMarkdown("sarc_references.md"))
              )
      )
    )
  )
)

#### -- Server -- ##############################################################
server <- function(input, output, session) {

  sim_arm <- function(arm_name = "baseline") {
    if (is.null(mod)) return(NULL)
    cov_vec <- c(
      AGE = input$age, SEX = as.numeric(input$sex), WT = input$wt, HT_M = input$ht,
      ACTIVITY = input$act, BASE_ALM = input$alm, BASE_GRIP = input$grip,
      BASE_GAIT = input$gait, BASE_SPPB = input$sppb, BASE_IGF1 = input$igf1,
      BASE_VITD = input$vitd, BASE_MSTN = input$mstn, BASE_IL6 = input$il6,
      BASE_GDF15 = input$gdf, RT_SESSIONS = input$rt_sess
    )
    events <- list()
    if (input$bima_dose > 0)
      events[["bima"]] <- ev(amt = input$bima_dose * input$wt, ii = 24*input$bima_ii,
                             addl = input$bima_n - 1, cmt = "BIMA_CENT")
    if (input$apit_dose > 0)
      events[["apit"]] <- ev(amt = input$apit_dose * input$wt, ii = 24*input$apit_ii,
                             addl = input$apit_n - 1, cmt = "APIT_CENT")
    if (input$test_dose > 0)
      events[["test"]] <- ev(amt = input$test_dose, ii = 24*7, addl = 51,
                             cmt = "TEST_DEPOT")
    if (input$vitd_dose > 0)
      events[["vitd"]] <- ev(amt = input$vitd_dose, ii = 24, addl = 364,
                             cmt = "VITD_GUT")
    if (input$leu_dose > 0)
      events[["leu"]]  <- ev(amt = input$leu_dose, ii = 6, addl = 4*365,
                             cmt = "LEU_GUT")

    e_all <- if (length(events) == 0) ev(time = 0, amt = 0, cmt = "LEU_GUT")
             else Reduce(`+`, events)
    out <- mod %>% param(cov_vec) %>% ev(e_all) %>%
      mrgsim(end = 24 * 365, delta = 24) %>% as_tibble()
    out$arm <- arm_name
    out
  }

  base_sim <- reactive(sim_arm("user"))

  output$dx_status <- renderPrint({
    cutoff_almi <- if (as.numeric(input$sex) == 1) 6.0 else 7.0
    cutoff_grip <- if (as.numeric(input$sex) == 1) 16  else 27
    almi <- input$alm / (input$ht^2)
    cat("ALMI =", round(almi,2), "kg/m² (cutoff:", cutoff_almi, ")\n",
        "Grip =", input$grip, "kg (cutoff:", cutoff_grip, ")\n",
        "Gait =", input$gait, "m/s (severe sarcopenia <0.8)\n",
        "→",
        if (almi < cutoff_almi && input$grip < cutoff_grip) {
          if (input$gait < 0.8) "Severe sarcopenia" else "Confirmed sarcopenia"
        } else "Below diagnostic threshold")
  })

  output$pk_plot <- renderPlotly({
    df <- base_sim()
    if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~CBIMA, type="scatter", mode="lines",
            name="Bimagrumab (mg/L)") %>%
      add_lines(y = ~CAPIT, name="Apitegromab (mg/L)") %>%
      add_lines(y = ~CTEST/10, name="Testosterone (×0.1 ng/mL)") %>%
      add_lines(y = ~CVITD/2, name="25(OH)D (×0.5 ng/mL)") %>%
      layout(xaxis = list(title = "Weeks"), yaxis = list(title = "Concentration"))
  })

  output$anab_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~ANABDRIVE, type="scatter", mode="lines",
            name="Anabolic drive") %>%
      layout(xaxis = list(title="Weeks"), yaxis = list(title="Anabolic drive (a.u.)"))
  })
  output$catab_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~CATDRIVE, type="scatter", mode="lines",
            name="Catabolic drive") %>%
      layout(xaxis = list(title="Weeks"), yaxis = list(title="Catabolic drive (a.u.)"))
  })
  output$endo_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~IGF1, type="scatter", mode="lines",
            name="IGF-1 (ng/mL)") %>%
      add_lines(y = ~MSTN_TOT*10, name="MSTN ×10 (ng/mL)") %>%
      add_lines(y = ~IL6*10,      name="IL-6 ×10 (pg/mL)") %>%
      add_lines(y = ~GDF15/10,    name="GDF-15 /10 (pg/mL)") %>%
      layout(xaxis = list(title="Weeks"))
  })

  output$alm_plot  <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~ALM, mode = "lines", type="scatter") %>%
      layout(title = "Appendicular Lean Mass (kg)",
             xaxis = list(title="Weeks"), yaxis = list(title="ALM"))
  })
  output$grip_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~GRIP, mode = "lines", type="scatter") %>%
      layout(title = "Grip strength (kg)", xaxis = list(title="Weeks"))
  })
  output$gait_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~GAIT, mode = "lines", type="scatter") %>%
      layout(title = "Gait speed (m/s)", xaxis = list(title="Weeks"))
  })
  output$sppb_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~SPPB_S, mode = "lines", type="scatter") %>%
      layout(title = "SPPB score", xaxis = list(title="Weeks"))
  })

  output$almi_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    cutoff <- if (as.numeric(input$sex) == 1) 6.0 else 7.0
    plot_ly(df, x = ~time/24/7, y = ~ALMI, mode = "lines", type="scatter",
            name="ALMI") %>%
      add_lines(y = rep(cutoff, nrow(df)), name="EWGSOP2 cutoff",
                line=list(dash="dash")) %>%
      layout(xaxis = list(title="Weeks"))
  })
  output$fall_plot <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~CUM_FALL, name="Cum. falls", mode = "lines",
            type="scatter") %>%
      add_lines(y = ~CUM_FRAIL*5, name="Frailty index ×5",
                line=list(dash="dot")) %>%
      layout(xaxis = list(title="Weeks"))
  })
  output$dx_table <- renderDT({
    df <- base_sim(); if (is.null(df)) return(NULL)
    df %>% mutate(week = round(time/24/7, 1)) %>%
      select(week, ALM, ALMI, GRIP, GAIT, SPPB_S, SARC_DX, SEVERE) %>%
      filter(week %% 4 == 0)
  })

  # Scenario tab
  scen_data <- reactive({
    if (is.null(mod)) return(NULL)
    base_cov <- c(AGE = input$age, SEX = as.numeric(input$sex), WT = input$wt,
                  HT_M = input$ht, ACTIVITY = input$act, BASE_ALM = input$alm)
    do_sim <- function(label, evt, rt = 0) {
      mod %>% param(c(base_cov, RT_SESSIONS = rt)) %>% ev(evt) %>%
        mrgsim(end = 24*365, delta = 24) %>% as_tibble() %>%
        mutate(arm = label)
    }
    bind_rows(
      do_sim("No treatment", ev(time = 0, amt = 0, cmt = "LEU_GUT")),
      do_sim("Bimagrumab 30 mg/kg q4w",
             ev(amt = 30*input$wt, ii = 24*28, addl = 11, cmt = "BIMA_CENT")),
      do_sim("Apitegromab 20 mg/kg q4w",
             ev(amt = 20*input$wt, ii = 24*28, addl = 11, cmt = "APIT_CENT")),
      do_sim("Testosterone 100 mg/wk",
             ev(amt = 100, ii = 24*7, addl = 51, cmt = "TEST_DEPOT")),
      do_sim("Vit D 2000 IU + RT 3x/wk",
             ev(amt = 2000, ii = 24, addl = 364, cmt = "VITD_GUT"), rt = 3),
      do_sim("Combo: Bima + RT 3x/wk",
             ev(amt = 30*input$wt, ii = 24*28, addl = 11, cmt = "BIMA_CENT"),
             rt = 3)
    )
  })

  output$scen_alm  <- renderPlotly({
    df <- scen_data(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~ALM, color = ~arm, type="scatter",
            mode = "lines") %>% layout(xaxis = list(title="Weeks"))
  })
  output$scen_gait <- renderPlotly({
    df <- scen_data(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~GAIT, color = ~arm, type="scatter",
            mode = "lines") %>% layout(xaxis = list(title="Weeks"))
  })
  output$scen_table <- renderDT({
    df <- scen_data(); if (is.null(df)) return(NULL)
    df %>% group_by(arm) %>%
      summarize(ALM_year = last(ALM), Gait_year = last(GAIT),
                Grip_year = last(GRIP), SPPB_year = last(SPPB_S),
                Falls = last(CUM_FALL))
  })

  output$bio_panel <- renderPlotly({
    df <- base_sim(); if (is.null(df)) return(NULL)
    plot_ly(df, x = ~time/24/7, y = ~IGF1, name = "IGF-1", mode = "lines",
            type="scatter") %>%
      add_lines(y = ~MSTN_TOT*10, name = "MSTN ×10") %>%
      add_lines(y = ~IL6*10,      name = "IL-6 ×10") %>%
      add_lines(y = ~GDF15/10,    name = "GDF-15 /10") %>%
      add_lines(y = ~CVITD,       name = "25(OH)D") %>%
      layout(xaxis = list(title="Weeks"),
             yaxis = list(title="Biomarker (scaled)"))
  })
}

#### -- Run -- ##################################################################
shinyApp(ui, server)
