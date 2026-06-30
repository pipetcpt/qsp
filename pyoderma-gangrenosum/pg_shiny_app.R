# =============================================================================
# Pyoderma Gangrenosum (PG) — QSP Shiny Dashboard
# =============================================================================
# Tabs:
#   1. Patient Profile     — demographics, comorbidities, baseline ulcer
#   2. Drug PK              — concentration–time profiles for 6 agents
#   3. PD — Cytokines       — TNFα, IL-1β, IL-17A, IL-6, IL-23, IL-8
#   4. Cellular Inflammation — Neutrophils, NET, Th17/Treg, MMP9, ROS, CRP, Calprotectin
#   5. Clinical Endpoints   — Ulcer area, PARACELSUS, Pain VAS, DLQI, %healing
#   6. Scenario Comparison  — comparative panel + endpoint table
#   7. Virtual Population   — VPOP distribution of complete-healing time
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(plotly)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Load model from the same directory
source("pg_mrgsolve_model.R")

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- dashboardPage(skin = "purple",
  dashboardHeader(title = "Pyoderma Gangrenosum QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",        tabName = "profile",   icon = icon("user-injured")),
      menuItem("Drug PK",                 tabName = "pk",        icon = icon("syringe")),
      menuItem("PD — Cytokines",          tabName = "cyto",      icon = icon("dna")),
      menuItem("Cellular Inflammation",   tabName = "cells",     icon = icon("microscope")),
      menuItem("Clinical Endpoints",      tabName = "endpoints", icon = icon("notes-medical")),
      menuItem("Scenario Comparison",     tabName = "scenarios", icon = icon("balance-scale")),
      menuItem("Virtual Population",      tabName = "vpop",      icon = icon("users"))
    ),
    hr(),
    selectInput("scenario", "Treatment scenario:",
                choices = c("NoTreatment","Prednisone_SOC","Cyclosporine",
                            "Infliximab","Adalimumab","Anakinra","Ustekinumab",
                            "Combo_PRED_CSA","Combo_IFX_low_CsA"),
                selected = "Infliximab"),
    sliderInput("dsev",        "Disease severity:", min = 1, max = 5, value = 3,  step = 0.1),
    sliderInput("baseUlcer",   "Baseline ulcer area (cm²):", min = 2, max = 40, value = 12),
    checkboxInput("pathergy",  "Recent pathergy event",  TRUE),
    checkboxInput("ibd",       "IBD comorbidity",        FALSE),
    sliderInput("tmax",        "Simulation horizon (days):", min = 28, max = 364, value = 168, step = 14),
    actionButton("run",        "Run simulation", icon = icon("play"))
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .info-box-icon { padding-top: 18px; }
      .small-box .icon-large { font-size: 60px; }
    "))),
    tabItems(

      # ---------------- TAB 1: Patient profile ----------------
      tabItem("profile",
        fluidRow(
          box(width = 6, title = "Disease & Patient Inputs", status = "primary", solidHeader = TRUE,
              textInput("age",  "Age (years):", "45"),
              selectInput("sex","Sex:", c("Female","Male"), "Female"),
              numericInput("weight","Weight (kg):", 70, min = 40, max = 150),
              selectInput("variant","PG variant:",
                          c("Classic ulcerative","Bullous","Pustular","Vegetative","Peristomal","PAPA/PASH")),
              selectInput("location","Anatomic site:",
                          c("Lower leg","Trunk","Peristomal","Genital","Face","Hand")),
              numericInput("priorRx","Number of prior immunosuppressants failed:", 1, 0, 5)),
          box(width = 6, title = "Baseline disease severity", status = "warning", solidHeader = TRUE,
              h4(paste0("Disease severity index: ")), verbatimTextOutput("severityIdx"),
              p("Severity = baseline ulcer cm² × cytokine burden × comorbidity coefficient."),
              tags$ul(
                tags$li("PAPA syndrome → PSTPIP1-driven pyrin inflammasome"),
                tags$li("IBD comorbidity → bacterial translocation, ↑TLR4 drive"),
                tags$li("Pathergy → trauma-induced amplification of NETosis"),
                tags$li("Monoclonal gammopathy → paraprotein-induced neutrophil activation"))),
        ),
        fluidRow(
          box(width = 12, title = "PARACELSUS scoring guide (Jockenhöfer 2019)", solidHeader = TRUE,
              p("Score 0–60.  Components: progression, exclusion of differential dx, erythematous violaceous edge, ",
                "improvement on immunosuppression, characteristic irregular shape, pain VAS, ",
                "size, undermined wound margin, deep ulcer, pathergy."))
        )
      ),

      # ---------------- TAB 2: Drug PK ----------------
      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Drug concentration-time profile", solidHeader = TRUE, status = "primary",
              plotlyOutput("pkPlot", height = "500px"))
        ),
        fluidRow(
          box(width = 12, title = "Drug PK parameters",
              DTOutput("pkTable"))
        )
      ),

      # ---------------- TAB 3: Cytokines ----------------
      tabItem("cyto",
        fluidRow(
          box(width = 12, title = "Cytokine trajectories (pg/mL)", solidHeader = TRUE, status = "danger",
              plotlyOutput("cytokinePlot", height = "550px"))
        ),
        fluidRow(
          box(width = 6, title = "Final cytokine state @ end of simulation",
              tableOutput("cytoTable")),
          box(width = 6, title = "% suppression vs. NoTreatment (Week 12)",
              plotlyOutput("supBar", height = "320px"))
        )
      ),

      # ---------------- TAB 4: Cellular inflammation ----------------
      tabItem("cells",
        fluidRow(
          box(width = 6, title = "Neutrophils, NET, MMP-9", solidHeader = TRUE, status = "warning",
              plotlyOutput("neutPlot", height = "350px")),
          box(width = 6, title = "Th17 / Treg / M1 macrophages", solidHeader = TRUE, status = "warning",
              plotlyOutput("th17Plot", height = "350px"))
        ),
        fluidRow(
          box(width = 6, title = "ROS / NET (oxidative axis)",
              plotlyOutput("rosPlot", height = "300px")),
          box(width = 6, title = "Serum biomarkers — CRP & Calprotectin",
              plotlyOutput("crpPlot", height = "300px"))
        )
      ),

      # ---------------- TAB 5: Endpoints ----------------
      tabItem("endpoints",
        fluidRow(
          valueBoxOutput("vbUlcer", width = 3),
          valueBoxOutput("vbPara",  width = 3),
          valueBoxOutput("vbPain",  width = 3),
          valueBoxOutput("vbHeal",  width = 3)
        ),
        fluidRow(
          box(width = 12, title = "Ulcer area & cumulative healed fraction",
              solidHeader = TRUE, status = "success",
              plotlyOutput("ulcerPlot", height = "400px"))
        ),
        fluidRow(
          box(width = 6, title = "PARACELSUS / Pain VAS / DLQI", plotlyOutput("scorePlot", height = "350px")),
          box(width = 6, title = "% achieving complete healing (binary trajectory)",
              plotlyOutput("compHealPlot", height = "350px"))
        )
      ),

      # ---------------- TAB 6: Scenario comparison ----------------
      tabItem("scenarios",
        fluidRow(
          box(width = 12, title = "Comparative ulcer trajectory across all scenarios",
              solidHeader = TRUE, status = "primary",
              plotlyOutput("scenPlot", height = "500px"))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint table @ Week 24",
              DTOutput("scenTable"))
        )
      ),

      # ---------------- TAB 7: Virtual population ----------------
      tabItem("vpop",
        fluidRow(
          box(width = 6, title = "Virtual population settings",
              numericInput("nvpop", "Number of patients:", 50, 10, 500),
              selectInput("vp_scenario", "Treatment for VPOP:",
                          choices = c("Infliximab","Adalimumab","Cyclosporine","Anakinra","Ustekinumab"),
                          selected = "Infliximab"),
              actionButton("runVPop", "Run VPOP simulation", icon = icon("rocket"))),
          box(width = 6, title = "Notes",
              p("Sources of variability: CL_ADA, CL_IFX, CL_CSA (log-normal CV 25–30%); ",
                "baseline ulcer area (log-normal CV 50%); disease severity (Gaussian); ",
                "pathergy / IBD (Bernoulli)."))
        ),
        fluidRow(
          box(width = 6, title = "Distribution of healed fraction @ Week 24",
              plotlyOutput("vpopDistr", height = "400px")),
          box(width = 6, title = "PARACELSUS distribution over time",
              plotlyOutput("vpopParacel", height = "400px"))
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  simResult <- reactiveVal()

  observeEvent(input$run, {
    m <- pg_qsp_model %>% param(DSEV          = input$dsev,
                                UlcerArea_0   = input$baseUlcer,
                                F_pathergy    = as.numeric(input$pathergy),
                                F_comorb_IBD  = as.numeric(input$ibd))
    df <- simulate_pg(m, input$scenario, tmax = input$tmax)
    simResult(df)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  output$severityIdx <- renderText({
    sprintf("Index = %.1f", input$dsev * input$baseUlcer/10 * (1 + 0.4*input$ibd + 0.5*input$pathergy))
  })

  # ----- TAB 2: PK -----
  output$pkPlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>%
      select(time, CONC_ADA, CONC_IFX, CONC_ANA, CONC_CSA, CONC_UST, CONC_PRED) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc")
    plot_ly(df, x = ~time, y = ~conc, color = ~drug, type = "scatter", mode = "lines") %>%
      layout(yaxis = list(title = "Concentration", type = "log"),
             xaxis = list(title = "Time (days)"))
  })
  output$pkTable <- renderDT({
    data.frame(
      Drug = c("Adalimumab","Infliximab","Anakinra","Cyclosporine","Ustekinumab","Prednisone"),
      Route = c("SC","IV","SC","PO","SC","PO"),
      Dose  = c("80→40 mg q1w","5 mg/kg @0,2,6 wk then q8w","100 mg/d","4 mg/kg/d","90 mg q12w","60 mg→taper"),
      CL_L_d= c(0.31,0.32,17.0,27.0,0.45,90.0),
      V_L   = c(4.7,3.5,20.0,85.0,4.6,35.0),
      t_half_d = c(14,9,0.3,0.8,21,0.2),
      Target = c("TNF-α","TNF-α","IL-1Rα","Calcineurin/NFAT","IL-12/23 p40","GR (NF-κB)")
    )
  }, options = list(pageLength = 6))

  # ----- TAB 3: Cytokines -----
  output$cytokinePlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>%
      select(time, TNFa, IL1b, IL17A, IL6, IL23, IL8) %>%
      pivot_longer(-time, names_to = "cytokine", values_to = "conc")
    plot_ly(df, x = ~time, y = ~conc, color = ~cytokine, type = "scatter", mode = "lines") %>%
      layout(yaxis = list(title = "Concentration (pg/mL)"),
             xaxis = list(title = "Time (days)"))
  })
  output$cytoTable <- renderTable({
    req(simResult())
    df <- simResult()
    tail(df,1) %>% select(TNFa, IL1b, IL17A, IL6, IL23, IL8) %>% round(2)
  })
  output$supBar <- renderPlotly({
    req(simResult())
    df <- simResult() %>% filter(abs(time - 84) < 1) %>% slice(1)
    baseline <- c(TNFa = 18, IL1b = 22, IL17A = 9, IL6 = 14, IL23 = 6, IL8 = 30)
    final <- as.numeric(df[1, c("TNFa","IL1b","IL17A","IL6","IL23","IL8")])
    sup <- pmax(0, 100*(baseline - final)/baseline)
    plot_ly(x = names(baseline), y = sup, type = "bar", marker=list(color="#8B008B")) %>%
      layout(yaxis = list(title = "% suppression vs. untreated reference"))
  })

  # ----- TAB 4: Cellular -----
  output$neutPlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, Neutroph, NET, MMP9_act) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines") %>%
      layout(yaxis = list(title="Index (au)"))
  })
  output$th17Plot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, Th17, Treg, M1) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines") %>%
      layout(yaxis = list(title="Index (au)"))
  })
  output$rosPlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, ROS_les, NET) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines")
  })
  output$crpPlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, CRP, Calprot) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines")
  })

  # ----- TAB 5: Endpoints -----
  output$ulcerPlot <- renderPlotly({
    req(simResult())
    df <- simResult()
    plot_ly() %>%
      add_lines(data = df, x = ~time, y = ~Ulcer, name = "Ulcer area (cm²)") %>%
      add_lines(data = df, x = ~time, y = ~Healed * 100, name = "% healed", yaxis = "y2") %>%
      layout(
        yaxis  = list(title = "Ulcer area (cm²)"),
        yaxis2 = list(title = "% healed", overlaying = "y", side = "right"))
  })
  output$scorePlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, PARACELSUS, Pain_VAS, DLQI) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines")
  })
  output$compHealPlot <- renderPlotly({
    req(simResult())
    df <- simResult() %>% select(time, CompleteHeal, HiSCRpseudo) %>% pivot_longer(-time)
    plot_ly(df, x = ~time, y = ~value, color = ~name, type="scatter", mode="lines") %>%
      layout(yaxis = list(title = "Response status (0/1)"))
  })

  # Value boxes
  output$vbUlcer <- renderValueBox({
    req(simResult())
    v <- tail(simResult()$Ulcer, 1)
    valueBox(sprintf("%.1f cm²", v), "Ulcer area (end)", icon = icon("notes-medical"), color = "red")
  })
  output$vbPara  <- renderValueBox({
    req(simResult())
    valueBox(sprintf("%.1f", tail(simResult()$PARACELSUS, 1)),
             "PARACELSUS", icon = icon("clipboard"), color = "yellow")
  })
  output$vbPain  <- renderValueBox({
    req(simResult())
    valueBox(sprintf("%.1f / 10", tail(simResult()$Pain_VAS, 1)),
             "Pain VAS", icon = icon("hand-holding-medical"), color = "orange")
  })
  output$vbHeal  <- renderValueBox({
    req(simResult())
    valueBox(sprintf("%.0f%%", 100*tail(simResult()$Healed, 1)),
             "% healed", icon = icon("heartbeat"), color = "green")
  })

  # ----- TAB 6: Scenario comparison -----
  scenAll <- reactive({
    isolate({
      m <- pg_qsp_model %>% param(DSEV          = input$dsev,
                                  UlcerArea_0   = input$baseUlcer,
                                  F_pathergy    = as.numeric(input$pathergy),
                                  F_comorb_IBD  = as.numeric(input$ibd))
      bind_rows(lapply(
        c("NoTreatment","Prednisone_SOC","Cyclosporine","Infliximab",
          "Adalimumab","Anakinra","Ustekinumab","Combo_PRED_CSA","Combo_IFX_low_CsA"),
        function(s) simulate_pg(m, s, tmax = input$tmax)))
    })
  })
  output$scenPlot <- renderPlotly({
    df <- scenAll()
    plot_ly(df, x = ~time, y = ~Ulcer, color = ~scenario, type="scatter", mode="lines") %>%
      layout(yaxis = list(title = "Ulcer area (cm²)"),
             xaxis = list(title = "Time (days)"))
  })
  output$scenTable <- renderDT({
    df <- scenAll() %>% filter(abs(time - input$tmax) < 0.6)
    df %>%
      group_by(scenario) %>%
      summarise(Ulcer_cm2 = mean(Ulcer),
                Healed_pct = 100*mean(Healed),
                PARACELSUS = mean(PARACELSUS),
                Pain = mean(Pain_VAS),
                CRP = mean(CRP)) %>%
      arrange(Ulcer_cm2) %>%
      mutate(across(where(is.numeric), ~ round(.,1)))
  })

  # ----- TAB 7: VPop -----
  vpopRes <- reactiveVal()
  observeEvent(input$runVPop, {
    withProgress(message = "Running virtual population...", value = 0.3, {
      vp <- run_vpop_scenario(pg_qsp_model, input$vp_scenario,
                              n = input$nvpop, tmax = input$tmax)
      vpopRes(vp)
    })
  })
  output$vpopDistr <- renderPlotly({
    req(vpopRes())
    df <- vpopRes() %>% filter(abs(time - input$tmax) < 0.6) %>%
      group_by(ID) %>% summarise(healed = max(Healed))
    plot_ly(df, x = ~healed, type = "histogram", nbinsx = 20,
            marker = list(color = "#228B22"))
  })
  output$vpopParacel <- renderPlotly({
    req(vpopRes())
    df <- vpopRes()
    plot_ly(df, x = ~time, y = ~PARACELSUS, split = ~ID,
            type="scatter", mode="lines",
            line = list(color="rgba(139,0,139,0.2)")) %>%
      layout(showlegend = FALSE, yaxis = list(title = "PARACELSUS"))
  })
}

# -----------------------------------------------------------------------------
# Launch
# -----------------------------------------------------------------------------
shinyApp(ui, server)
