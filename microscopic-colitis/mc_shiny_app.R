# =====================================================================
# Microscopic Colitis QSP — Shiny Dashboard
# Tabs: Patient profile · Drug PK · Mucosal immunology · Histology
#       Clinical endpoints · Scenario comparison · Long-term outcomes
# =====================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("mc_mrgsolve_model.R", local = TRUE)

# ---- Helpers --------------------------------------------------------
scen_choices <- names(mc_scenarios)

ui <- dashboardPage(
  dashboardHeader(title = "Microscopic Colitis QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient",        tabName = "pt",     icon = icon("user")),
      menuItem("② Drug PK",        tabName = "pk",     icon = icon("flask")),
      menuItem("③ Mucosal immun.", tabName = "imm",    icon = icon("dna")),
      menuItem("④ Histology",      tabName = "hist",   icon = icon("microscope")),
      menuItem("⑤ Clinical",       tabName = "clin",   icon = icon("notes-medical")),
      menuItem("⑥ Scenario cmp",   tabName = "scen",   icon = icon("balance-scale")),
      menuItem("⑦ Long-term",      tabName = "long",   icon = icon("chart-line"))
    ),
    hr(),
    selectInput("subtype", "Subtype", choices = c("Collagenous (CC)" = "CC",
                                                  "Lymphocytic (LC)" = "LC")),
    selectInput("scenario", "Treatment scenario", choices = scen_choices,
                selected = "02_budesonide_taper"),
    sliderInput("days", "Simulation horizon (days)", 30, 730, 365, step = 30),
    actionButton("run", "Run simulation", icon = icon("play"),
                 class = "btn-primary")
  ),
  dashboardBody(
    tabItems(
      tabItem("pt",
        fluidRow(
          box(width = 6, title = "Baseline patient profile", status = "primary",
              numericInput("age", "Age (yr)",   65, 18, 95),
              selectInput("sex", "Sex", c("Female", "Male"), "Female"),
              numericInput("stools_bl", "Baseline stools/day", 6, 1, 20),
              numericInput("col_bl",    "Baseline collagen band (μm)", 18, 0, 40),
              numericInput("iel_bl",    "Baseline CD8+ IEL per 100 EC", 35, 0, 100),
              checkboxInput("nsaid", "Concomitant NSAID use", FALSE),
              checkboxInput("ppi",   "Concomitant PPI use",   TRUE),
              checkboxInput("ssri",  "Concomitant SSRI use",  FALSE)
          ),
          box(width = 6, title = "Summary at chosen horizon", status = "info",
              tableOutput("summary_tbl"))
        )
      ),
      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Drug concentration vs time",
              plotOutput("plot_pk", height = 420))
        )
      ),
      tabItem("imm",
        fluidRow(
          box(width = 6, title = "Cytokines",  plotOutput("plot_cyto", height = 380)),
          box(width = 6, title = "CD8+ IEL",   plotOutput("plot_iel",  height = 380))
        )
      ),
      tabItem("hist",
        fluidRow(
          box(width = 6, title = "Subepithelial collagen band (μm)",
              plotOutput("plot_col", height = 380)),
          box(width = 6, title = "Barrier integrity (0–1)",
              plotOutput("plot_bar", height = 380))
        )
      ),
      tabItem("clin",
        fluidRow(
          box(width = 6, title = "Stool frequency (per day)",
              plotOutput("plot_stool", height = 380)),
          box(width = 6, title = "Hjortswang remission probability",
              plotOutput("plot_hj",    height = 380))
        ),
        fluidRow(
          box(width = 6, title = "HRQoL (composite)",
              plotOutput("plot_qol", height = 320)),
          box(width = 6, title = "Net colonic water flux (mL/day)",
              plotOutput("plot_wat", height = 320))
        )
      ),
      tabItem("scen",
        fluidRow(
          box(width = 12, title = "Scenario comparison — stool frequency",
              plotOutput("plot_scen", height = 500))
        )
      ),
      tabItem("long",
        fluidRow(
          box(width = 6, title = "BMD trajectory (T-score Δ)",
              plotOutput("plot_bmd", height = 320)),
          box(width = 6, title = "HPA axis recovery",
              plotOutput("plot_hpa", height = 320))
        ),
        fluidRow(
          box(width = 12, title = "Year-end status table",
              tableOutput("long_tbl"))
        )
      )
    )
  )
)

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    pset <- list(
      SUBTYPE_CC = ifelse(input$subtype == "CC", 1, 0),
      IEL_BASE   = input$iel_bl,
      COL_BASE   = ifelse(input$subtype == "CC", input$col_bl, input$col_bl/3),
      STOOL_BASE = input$stools_bl
    )
    evd <- mc_scenarios[[input$scenario]]
    out <- mc_model %>% param(pset) %>%
      mrgsim(events = evd, end = input$days * 24, delta = 24) %>%
      as.data.frame()
    out$day <- out$time / 24
    out
  }, ignoreNULL = FALSE)

  output$summary_tbl <- renderTable({
    df <- sim()
    last <- tail(df, 1)
    data.frame(
      Variable = c("Stools/day", "Hjortswang score", "Barrier",
                   "Collagen band (μm)", "Bile-acid load", "HRQoL"),
      Value    = round(c(last$STOOL_CLIP, last$HJ, last$BAR,
                          last$COL, last$BA, last$QOL), 2)
    )
  })

  output$plot_pk <- renderPlot({
    df <- sim()
    df_long <- df %>%
      select(day, BUD_CONC, IFX_C, ADA_C, VDZ_C, TGN) %>%
      pivot_longer(-day, names_to = "Compound", values_to = "Conc")
    ggplot(df_long, aes(day, Conc, colour = Compound)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~Compound, scales = "free_y") +
      labs(x = "Day", y = "Concentration (units depend on drug)") +
      theme_minimal()
  })

  output$plot_cyto <- renderPlot({
    df <- sim()
    df_long <- df %>%
      select(day, IFNG, TNF, IL6, IL17, TGFB) %>%
      pivot_longer(-day, names_to = "Cytokine", values_to = "Conc")
    ggplot(df_long, aes(day, Conc, colour = Cytokine)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "pg/mL") + theme_minimal()
  })

  output$plot_iel <- renderPlot({
    ggplot(sim(), aes(day, IEL)) +
      geom_line(colour = "#1f618d", linewidth = 1) +
      geom_hline(yintercept = 20, linetype = 2) +
      labs(x = "Day", y = "CD8+ IEL per 100 epithelial cells") +
      theme_minimal()
  })

  output$plot_col <- renderPlot({
    ggplot(sim(), aes(day, COL)) +
      geom_line(colour = "#a04000", linewidth = 1) +
      geom_hline(yintercept = 10, linetype = 2) +
      labs(x = "Day", y = "Subepithelial collagen (μm)") +
      theme_minimal()
  })

  output$plot_bar <- renderPlot({
    ggplot(sim(), aes(day, BAR)) +
      geom_line(colour = "#1e8449", linewidth = 1) +
      ylim(0, 1) +
      labs(x = "Day", y = "Barrier integrity (0–1)") +
      theme_minimal()
  })

  output$plot_stool <- renderPlot({
    ggplot(sim(), aes(day, STOOL_CLIP)) +
      geom_line(colour = "#922b21", linewidth = 1) +
      geom_hline(yintercept = 3, linetype = 2) +
      labs(x = "Day", y = "Stools/day") + theme_minimal()
  })

  output$plot_hj <- renderPlot({
    ggplot(sim(), aes(day, HJ)) +
      geom_line(colour = "#7d3c98", linewidth = 1) +
      labs(x = "Day", y = "Hjortswang active score") +
      theme_minimal()
  })

  output$plot_qol <- renderPlot({
    ggplot(sim(), aes(day, QOL)) +
      geom_line(colour = "#0e6655", linewidth = 1) +
      labs(x = "Day", y = "HRQoL composite") + theme_minimal()
  })

  output$plot_wat <- renderPlot({
    ggplot(sim(), aes(day, WAT)) +
      geom_line(colour = "#117a65", linewidth = 1) +
      labs(x = "Day", y = "Net colonic water flux (mL/day)") +
      theme_minimal()
  })

  output$plot_scen <- renderPlot({
    pset <- list(SUBTYPE_CC = ifelse(input$subtype == "CC", 1, 0))
    runs <- lapply(scen_choices, function(s) {
      out <- mc_model %>% param(pset) %>%
        mrgsim(events = mc_scenarios[[s]],
               end = input$days * 24, delta = 24) %>%
        as.data.frame()
      out$day <- out$time/24
      out$scenario <- s
      out
    })
    df <- do.call(rbind, runs)
    ggplot(df, aes(day, STOOL_CLIP, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "Stools/day", colour = "Scenario") +
      theme_minimal()
  })

  output$plot_bmd <- renderPlot({
    ggplot(sim(), aes(day, BMD)) +
      geom_line(colour = "#7d6608", linewidth = 1) +
      labs(x = "Day", y = "ΔT-score (proxy)") + theme_minimal()
  })

  output$plot_hpa <- renderPlot({
    ggplot(sim(), aes(day, HPA)) +
      geom_line(colour = "#5b2c6f", linewidth = 1) +
      ylim(0, 1.05) +
      labs(x = "Day", y = "HPA function (1 = intact)") +
      theme_minimal()
  })

  output$long_tbl <- renderTable({
    df <- sim()
    last <- tail(df, 1)
    data.frame(
      Endpoint = c("Stools/day", "Hjortswang", "BMD Δ", "HPA",
                   "Collagen (μm)", "Bile-acid load", "HRQoL"),
      Value    = round(c(last$STOOL_CLIP, last$HJ, last$BMD, last$HPA,
                          last$COL, last$BA, last$QOL), 3)
    )
  })
}

shinyApp(ui, server)
