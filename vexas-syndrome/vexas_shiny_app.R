# ============================================================================
# VEXAS Syndrome QSP — Shiny app
# Interactive simulator wrapping the mrgsolve model in `vexas_mrgsolve_model.R`.
# 8 tabs: Patient · UBA1 clone · Cytokines · Hematology · Clinical activity
#         · Drug PK · Scenarios · Biomarker dashboard
# ============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("vexas_mrgsolve_model.R")  # provides mod_vexas, scenarios()

# ----------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "VEXAS QSP Simulator"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient profile",  tabName = "pt",      icon = icon("user")),
      menuItem("UBA1 clone & VAF", tabName = "clone",   icon = icon("dna")),
      menuItem("Cytokine storm",   tabName = "cyto",    icon = icon("fire")),
      menuItem("Hematology",       tabName = "heme",    icon = icon("tint")),
      menuItem("Clinical activity",tabName = "clin",    icon = icon("notes-medical")),
      menuItem("Drug PK",          tabName = "pk",      icon = icon("pills")),
      menuItem("Scenario compare", tabName = "scen",    icon = icon("balance-scale")),
      menuItem("Biomarker panel",  tabName = "biomark", icon = icon("chart-line"))
    ),
    hr(),
    h4("Therapy", style = "padding-left:15px;color:white"),
    selectInput("scenario", "Regimen",
                choices = c("Untreated"        = "1_untreated",
                            "Prednisone taper" = "2_pred_taper",
                            "Tocilizumab + GC" = "3_toci_lowGC",
                            "Anakinra"         = "4_anakinra",
                            "Ruxolitinib + GC" = "5_ruxo_pred",
                            "Azacitidine"      = "6_azacitidine",
                            "Allogeneic HSCT"  = "7_HSCT"),
                selected = "3_toci_lowGC"),
    sliderInput("duration", "Simulation (weeks)", 4, 52, 24, step = 4),
    sliderInput("VAF0",   "Baseline UBA1 VAF (%)", 10, 95, 60, step = 5),
    sliderInput("IL60",   "Baseline IL-6 (pg/mL)", 20, 500, 150, step = 10),
    sliderInput("HB0",    "Baseline Hb (g/dL)",    6,   14,  9,  step = 0.5),
    sliderInput("CRP0",   "Baseline CRP (mg/L)",   5,  300, 110, step = 5),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background-color:#7B1FA2;border-color:#4A148C")
  ),
  dashboardBody(
    tabItems(
      tabItem("pt",
              fluidRow(
                box(width = 12, title = "Patient profile inputs", status = "primary", solidHeader = TRUE,
                    DT::dataTableOutput("patient_tbl"),
                    helpText("VEXAS predominantly affects men >50 with somatic UBA1 mutations in HSCs; ",
                             "macrocytic anemia and CRP/ferritin elevation are diagnostic anchors."))
              )),
      tabItem("clone",
              fluidRow(
                box(width = 6, title = "UBA1 clonal dynamics (VAF)", plotOutput("plot_VAF")),
                box(width = 6, title = "Misfolded protein / ER stress / ROS", plotOutput("plot_proteo"))
              )),
      tabItem("cyto",
              fluidRow(
                box(width = 6, title = "IL-6 / IL-1β / TNF-α", plotOutput("plot_cyto1")),
                box(width = 6, title = "IFN-α · CXCL8 · CCL2", plotOutput("plot_cyto2"))
              )),
      tabItem("heme",
              fluidRow(
                box(width = 6, title = "Hemoglobin & Platelets", plotOutput("plot_heme")),
                box(width = 6, title = "ANC & Bone marrow vacuolation surrogate", plotOutput("plot_anc"))
              )),
      tabItem("clin",
              fluidRow(
                box(width = 6, title = "Fever / Skin-chondritis activity", plotOutput("plot_act")),
                box(width = 6, title = "VTE risk index",                   plotOutput("plot_vte"))
              )),
      tabItem("pk",
              fluidRow(
                box(width = 6, title = "Steroid + small-molecule PK",          plotOutput("plot_pk1")),
                box(width = 6, title = "Biologic PK & receptor occupancy",     plotOutput("plot_pk2"))
              )),
      tabItem("scen",
              fluidRow(
                box(width = 12, title = "Cross-regimen comparison",
                    plotOutput("plot_scenarios", height = 600))
              )),
      tabItem("biomark",
              fluidRow(
                box(width = 6, title = "CRP & Ferritin", plotOutput("plot_acu")),
                box(width = 6, title = "Composite VEXAS activity score",  plotOutput("plot_score")),
                box(width = 12, title = "Endpoint table at week 12",      DT::dataTableOutput("tbl_endpoint"))
              ))
    )
  )
)

# ----------------------------------------------------------------------------
server <- function(input, output, session) {

  patient_df <- reactive({
    data.frame(
      Variable = c("Sex", "Age (yr)", "VAF baseline (%)", "Hb (g/dL)",
                   "Platelets (×10⁹/L)", "ANC (×10⁹/L)", "CRP (mg/L)",
                   "Ferritin (ng/mL)", "Phenotype dominance"),
      Value = c("Male", 67, input$VAF0, input$HB0,
                130, 5.5, input$CRP0,
                1800, "Auricular chondritis + neutrophilic dermatosis")
    )
  })
  output$patient_tbl <- DT::renderDataTable(patient_df(), options = list(dom = "t"))

  run_sim <- eventReactive(input$run, {
    evts <- scenarios()[[input$scenario]]
    mod_local <- mod_vexas %>%
      param(VAF0 = input$VAF0/100) %>%
      init(VAF = input$VAF0/100, IL6 = input$IL60, HB = input$HB0, CRP = input$CRP0)
    if (input$scenario == "7_HSCT") {
      mod_local <- mod_local %>% param(k_HSCT = 0.5)
    }
    if (input$scenario == "6_azacitidine") {
      mod_local <- mod_local %>% param(k_Aza = 0.002)
    }
    if (is.null(evts)) {
      mod_local %>% mrgsim(end = input$duration*168, delta = 6) %>% as.data.frame()
    } else {
      mod_local %>% ev(evts) %>% mrgsim(end = input$duration*168, delta = 6) %>% as.data.frame()
    }
  })

  long <- function(df, cols) {
    df %>% select(time, all_of(cols)) %>%
      mutate(day = time/24) %>% select(-time) %>%
      pivot_longer(-day, names_to = "var", values_to = "y")
  }

  output$plot_VAF <- renderPlot({
    df <- run_sim()
    ggplot(df, aes(time/24, VAF*100)) + geom_line(color = "#6A1B9A", lwd = 1.1) +
      labs(x = "Day", y = "UBA1 VAF (%)", title = "Clonal dynamics") +
      theme_minimal(base_size = 13)
  })
  output$plot_proteo <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("MISF","ERST","ROS")), aes(day, y, color = var)) +
      geom_line(lwd = 1) + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "Stress index (au)", color = "Pathway")
  })
  output$plot_cyto1 <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("IL6","IL1B","TNFa")), aes(day, y, color = var)) +
      geom_line(lwd = 1) + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "pg/mL")
  })
  output$plot_cyto2 <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("IFNa","CXCL8","CCL2")), aes(day, y, color = var)) +
      geom_line(lwd = 1) + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "pg/mL")
  })
  output$plot_heme <- renderPlot({
    df <- run_sim()
    df2 <- long(df, c("HB","PLT"))
    ggplot(df2, aes(day, y, color = var)) + geom_line(lwd = 1.1) +
      facet_wrap(~var, scales = "free_y") + theme_minimal(base_size = 13) +
      labs(x = "Day", y = NULL)
  })
  output$plot_anc <- renderPlot({
    df <- run_sim()
    ggplot(df, aes(time/24, ANC)) + geom_line(color = "#5D4037", lwd = 1.1) +
      labs(x = "Day", y = "ANC (×10⁹/L)") + theme_minimal(base_size = 13)
  })
  output$plot_act <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("FEV","SKIN")), aes(day, y, color = var)) +
      geom_line(lwd = 1.1) + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "Activity index (0-100)")
  })
  output$plot_vte <- renderPlot({
    df <- run_sim()
    ggplot(df, aes(time/24, VTE)) + geom_line(color = "#C62828", lwd = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed") +
      labs(x = "Day", y = "VTE risk index", caption = "Threshold for prophylaxis at 50") +
      theme_minimal(base_size = 13)
  })
  output$plot_pk1 <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("C_PRED","C_RUX","C_AZA")), aes(day, y, color = var)) +
      geom_line() + facet_wrap(~var, scales = "free_y") + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "Concentration")
  })
  output$plot_pk2 <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("C_TOC","C_ANA","C_CAN","occ_TOC","occ_ANA","occ_RUX")),
           aes(day, y, color = var)) +
      geom_line() + facet_wrap(~var, scales = "free_y") + theme_minimal(base_size = 13) +
      labs(x = "Day", y = "")
  })
  output$plot_acu <- renderPlot({
    df <- run_sim()
    ggplot(long(df, c("CRP","FER")), aes(day, y, color = var)) +
      geom_line(lwd = 1.1) + facet_wrap(~var, scales = "free_y") +
      theme_minimal(base_size = 13) + labs(x = "Day", y = NULL)
  })
  output$plot_score <- renderPlot({
    df <- run_sim()
    df$score <- 0.3*df$FEV + 0.3*df$SKIN + 0.1*df$VTE +
                0.15*pmin(df$CRP/2,100) + 0.15*pmin(df$FER/100, 100)
    ggplot(df, aes(time/24, score)) + geom_line(color = "#4A148C", lwd = 1.2) +
      labs(x = "Day", y = "Composite activity score (0-100)") +
      theme_minimal(base_size = 13)
  })
  output$tbl_endpoint <- DT::renderDataTable({
    df <- run_sim()
    end_day <- input$duration * 7
    snap <- df %>% mutate(day = time/24) %>%
      filter(abs(day - end_day) <= 0.5) %>%
      slice(1) %>%
      select(VAF, HB, PLT, ANC, CRP, FER, IL6, IL1B, FEV, SKIN, VTE)
    snap %>%
      mutate(across(everything(), ~round(.x, 2)))
  }, options = list(dom = "t"))

  output$plot_scenarios <- renderPlot({
    snames <- names(scenarios())
    df_all <- lapply(snames, function(s) {
      mod_local <- mod_vexas %>%
        param(VAF0 = input$VAF0/100,
              k_HSCT = ifelse(s == "7_HSCT", 0.5, 0),
              k_Aza  = ifelse(s == "6_azacitidine", 0.002, 0)) %>%
        init(VAF = input$VAF0/100, IL6 = input$IL60, HB = input$HB0, CRP = input$CRP0)
      evts <- scenarios()[[s]]
      out <- if (is.null(evts)) mod_local %>% mrgsim(end = input$duration*168, delta = 24)
             else mod_local %>% ev(evts) %>% mrgsim(end = input$duration*168, delta = 24)
      as.data.frame(out) %>% mutate(scenario = s)
    }) %>% bind_rows()
    df_all %>%
      pivot_longer(c(IL6, CRP, HB, FER), names_to = "var", values_to = "y") %>%
      ggplot(aes(time/24, y, color = scenario)) +
      geom_line(lwd = 0.9) + facet_wrap(~var, scales = "free_y") +
      theme_minimal(base_size = 13) +
      labs(x = "Day", y = NULL, color = "Regimen")
  })
}

shinyApp(ui, server)
