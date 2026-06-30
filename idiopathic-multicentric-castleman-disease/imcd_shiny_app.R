# =============================================================================
# Idiopathic Multicentric Castleman Disease (iMCD) — Shiny Dashboard
# Interactive QSP exploration: drug PK / cytokine PD / clinical endpoints
# 8 tabs : Patient · PK · IL-6 axis · Acute phase · TAFRO/VEGF · Lymph node ·
#          Scenario comparison · Biomarker panel
# =============================================================================
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# Source the model (assumes companion file in same directory)
source_imcd <- function() {
  source("imcd_mrgsolve_model.R", local = TRUE, chdir = TRUE)
}

# Pre-load model
suppressMessages(source_imcd())

# Helper: run user-defined regimen
build_events <- function(input) {
  ev_all <- list()
  if (input$silt_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$silt_dose * input$weight, cmt = "SILT_C",
       time = seq(0, input$horizon, input$silt_int))
  if (input$tocz_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$tocz_dose * input$weight, cmt = "TOCZ_C",
       time = seq(0, input$horizon, input$tocz_int))
  if (input$siro_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$siro_dose, cmt = "SIRO_GUT",
       ii = 1, addl = input$horizon - 1)
  if (input$rtx_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$rtx_dose * 1.5, cmt = "RTX_C",
       time = seq(0, input$rtx_n * 7 - 1, 7))
  if (input$ana_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$ana_dose, cmt = "ANA_SC", ii = 1, addl = input$horizon - 1)
  if (input$pred_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$pred_dose, cmt = "PRED_GUT", ii = 1, addl = input$horizon - 1)
  if (input$rux_on) ev_all[[length(ev_all)+1]] <-
    ev(amt = input$rux_dose, cmt = "RUX_GUT", ii = 0.5, addl = (input$horizon-1)*2)
  if (length(ev_all) == 0) {
    return(ev(amt = 0, cmt = "SILT_C", time = 0))
  }
  Reduce(`+`, ev_all)
}

# ---------------------- UI -----------------------------------------------
ui <- fluidPage(
  titlePanel("Idiopathic Multicentric Castleman Disease (iMCD) — QSP Simulator"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Patient profile"),
      numericInput("weight", "Weight (kg)", 70, min = 30, max = 150),
      selectInput("subtype", "Disease subtype",
                  c("Hyperplastic GC", "Plasmacytic", "Mixed", "TAFRO"),
                  selected = "Mixed"),
      sliderInput("baseline_IL6", "Baseline serum IL-6 (pg/mL)",
                  10, 300, 60, step = 5),
      sliderInput("baseline_CRP", "Baseline CRP (mg/L)",
                  10, 250, 120, step = 5),
      sliderInput("baseline_LN", "Baseline LN composite (cm)",
                  3, 20, 8, step = 1),
      sliderInput("horizon", "Simulation horizon (days)",
                  90, 730, 365, step = 30),
      hr(),
      h4("Treatments"),
      checkboxInput("silt_on", "Siltuximab IV", TRUE),
      conditionalPanel("input.silt_on",
        numericInput("silt_dose", "mg/kg", 11, min = 3, max = 15),
        numericInput("silt_int",  "Interval (d)", 21, min = 14, max = 42)
      ),
      checkboxInput("tocz_on", "Tocilizumab IV", FALSE),
      conditionalPanel("input.tocz_on",
        numericInput("tocz_dose", "mg/kg", 8, min = 4, max = 12),
        numericInput("tocz_int",  "Interval (d)", 14, min = 14, max = 28)
      ),
      checkboxInput("siro_on", "Sirolimus PO", FALSE),
      conditionalPanel("input.siro_on",
        numericInput("siro_dose", "mg QD", 2, min = 1, max = 6)
      ),
      checkboxInput("rtx_on", "Rituximab IV induction", FALSE),
      conditionalPanel("input.rtx_on",
        numericInput("rtx_dose", "mg/m2", 375, min = 200, max = 500),
        numericInput("rtx_n",    "Weekly doses", 4,  min = 1, max = 8)
      ),
      checkboxInput("ana_on", "Anakinra SC", FALSE),
      conditionalPanel("input.ana_on",
        numericInput("ana_dose", "mg QD", 100, min = 50, max = 200)
      ),
      checkboxInput("rux_on", "Ruxolitinib PO BID", FALSE),
      conditionalPanel("input.rux_on",
        numericInput("rux_dose", "mg per dose", 20, min = 5, max = 25)
      ),
      checkboxInput("pred_on", "Prednisone PO", FALSE),
      conditionalPanel("input.pred_on",
        numericInput("pred_dose", "mg QD", 60, min = 5, max = 100)
      ),
      hr(),
      actionButton("run", "Run simulation", class = "btn-primary"),
      br(), br(),
      helpText("CDCN 2017 diagnostic criteria · van Rhee 2018 consensus treatment ·",
               "CONCERT trial (siltuximab) · Fajgenbaum 2019 JCI (sirolimus TAFRO)")
    ),
    mainPanel(width = 9,
      tabsetPanel(id = "tabs",
        tabPanel("① Patient profile",
          h4("Selected patient"),
          verbatimTextOutput("patient_text"),
          h4("CDCN 2017 Diagnostic Criteria"),
          tags$ul(
            tags$li("MAJOR (both required): Castleman lymphadenopathy ≥2 stations & LN histology"),
            tags$li("MINOR (≥2, ≥1 lab): CRP↑, ESR↑, IL-6↑, Hb↓, Plt↑/↓, IgG↑↑, Alb↓ /"),
            tags$li("                         Fever, splenomegaly, hepatomegaly, edema/effusion, eruptive hemangiomas, LIP")
          ),
          h4("Approved & off-label drugs modeled"),
          tableOutput("drug_summary")
        ),
        tabPanel("② Drug PK",
          plotOutput("plot_PK", height = 480)
        ),
        tabPanel("③ IL-6 axis (PD)",
          plotOutput("plot_IL6", height = 380),
          helpText("Siltuximab traps IL-6 → total IL-6 rises in serum (clearance block) but free IL-6 falls"),
          plotOutput("plot_mTOR", height = 280)
        ),
        tabPanel("④ Acute-phase response",
          plotOutput("plot_CRP",  height = 320),
          fluidRow(
            column(6, plotOutput("plot_Hb",  height = 280)),
            column(6, plotOutput("plot_IgG", height = 280))
          )
        ),
        tabPanel("⑤ TAFRO / VEGF / Anasarca",
          plotOutput("plot_VEGF",      height = 280),
          plotOutput("plot_Anasarca",  height = 280),
          plotOutput("plot_Plt",       height = 280)
        ),
        tabPanel("⑥ Lymph node burden",
          plotOutput("plot_LN", height = 360),
          plotOutput("plot_PB", height = 280),
          plotOutput("plot_Bmem", height = 280)
        ),
        tabPanel("⑦ Scenario comparison",
          h4("Pre-defined scenarios (running on default 70 kg patient)"),
          plotOutput("plot_scenarios_CRP", height = 320),
          plotOutput("plot_scenarios_LN", height = 320),
          plotOutput("plot_scenarios_response", height = 320)
        ),
        tabPanel("⑧ Biomarker panel & response",
          h4("Composite CDCN response (lower is better)"),
          plotOutput("plot_CDCN", height = 320),
          h4("Endpoint summary table"),
          DT::DTOutput("biomarker_table"),
          h4("Survival (cumulative)"),
          plotOutput("plot_surv", height = 280)
        )
      )
    )
  )
)

# ---------------------- Server ------------------------------------------
server <- function(input, output, session) {
  sim_data <- eventReactive(input$run, {
    ev <- build_events(input)
    # adjust initial conditions to match user inputs
    out <- mod %>%
      init(IL6_T = input$baseline_IL6, IL6_F = input$baseline_IL6 * 0.5,
           CRP = input$baseline_CRP, LN = input$baseline_LN) %>%
      mrgsim(events = ev, end = input$horizon, delta = 1) %>%
      as_tibble()
    out
  }, ignoreNULL = FALSE)

  output$patient_text <- renderText({
    paste0("Subtype: ", input$subtype, "  ·  Weight: ", input$weight, " kg",
           "\nBaseline IL-6: ", input$baseline_IL6, " pg/mL",
           "\nBaseline CRP: ", input$baseline_CRP, " mg/L",
           "\nBaseline LN composite: ", input$baseline_LN, " cm",
           "\n\nActive treatments:\n",
           paste(c(
             if (input$silt_on) sprintf("• Siltuximab %g mg/kg q%dd", input$silt_dose, input$silt_int),
             if (input$tocz_on) sprintf("• Tocilizumab %g mg/kg q%dd", input$tocz_dose, input$tocz_int),
             if (input$siro_on) sprintf("• Sirolimus %g mg QD",      input$siro_dose),
             if (input$rtx_on)  sprintf("• Rituximab %g mg/m² wk×%d", input$rtx_dose, input$rtx_n),
             if (input$ana_on)  sprintf("• Anakinra %g mg SC QD",    input$ana_dose),
             if (input$rux_on)  sprintf("• Ruxolitinib %g mg PO BID", input$rux_dose),
             if (input$pred_on) sprintf("• Prednisone %g mg PO QD",   input$pred_dose)
           ), collapse = "\n"))
  })

  output$drug_summary <- renderTable({
    data.frame(
      Drug = c("Siltuximab", "Tocilizumab", "Sirolimus", "Rituximab",
               "Anakinra", "Ruxolitinib", "CHOP-doxo/cyc", "Prednisone", "Thalidomide"),
      Target = c("IL-6 (direct)", "IL-6Rα (mIL6R + sIL6R)", "mTORC1 (FKBP12)",
                 "CD20+ B cell", "IL-1Rα", "JAK1/2",
                 "DNA / microtubule (plasmablast)", "GR → NF-κB",
                 "Cereblon (TNF↓ angiogenesis↓)"),
      Status = c("FDA-approved 2014 (HHV-8 neg)", "Off-label / JP-approved",
                 "Off-label TAFRO rescue", "Off-label",
                 "Off-label (TAFRO/inflammasome)", "Off-label investigational",
                 "Salvage cytotoxic", "Adjunct anti-inflammatory",
                 "Salvage refractory TAFRO")
    )
  })

  output$plot_PK <- renderPlot({
    d <- sim_data() %>%
      select(time, Siltuximab = Csilt, Tocilizumab = Ctocz,
             Sirolimus = Csiro, Rituximab = Crtx, Anakinra = Cana,
             Ruxolitinib = Crux) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc")
    ggplot(d, aes(time, Conc + 1e-3, color = Drug)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(title = "Drug concentration-time profiles (log scale)",
           x = "Day", y = "Concentration (mg/L or ng/mL)") +
      theme_bw()
  })

  output$plot_IL6 <- renderPlot({
    d <- sim_data() %>% select(time, Total = IL6_serum, Free = IL6_free) %>%
      pivot_longer(-time, names_to = "Form", values_to = "IL6")
    ggplot(d, aes(time, IL6, color = Form)) +
      geom_line(linewidth = 1) +
      labs(title = "Serum IL-6: total vs free (siltuximab → free IL-6 ↓ but total can rise)",
           x = "Day", y = "IL-6 (pg/mL)") + theme_bw()
  })

  output$plot_mTOR <- renderPlot({
    ggplot(sim_data(), aes(time, mTOR_act)) +
      geom_line(linewidth = 1, color = "#2980b9") +
      labs(title = "mTORC1 activity (sirolimus target — TAFRO driver)",
           x = "Day", y = "mTOR (0–1)") + theme_bw()
  })

  output$plot_CRP <- renderPlot({
    ggplot(sim_data(), aes(time, CRP_lab)) + geom_line(color = "firebrick", linewidth = 1) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "darkgreen") +
      labs(title = "CRP trajectory (treatment goal: < 10 mg/L)",
           x = "Day", y = "CRP (mg/L)") + theme_bw()
  })

  output$plot_Hb <- renderPlot({
    ggplot(sim_data(), aes(time, Hb_lab)) +
      geom_line(color = "#c0392b", linewidth = 1) +
      geom_hline(yintercept = 12, linetype = "dashed") +
      labs(title = "Hemoglobin recovery", x = "Day", y = "Hb (g/dL)") + theme_bw()
  })

  output$plot_IgG <- renderPlot({
    ggplot(sim_data(), aes(time, IgG_lab)) +
      geom_line(color = "#8e44ad", linewidth = 1) +
      labs(title = "Polyclonal IgG (g/dL)", x = "Day", y = "IgG") + theme_bw()
  })

  output$plot_VEGF <- renderPlot({
    ggplot(sim_data(), aes(time, VEGF_lab)) +
      geom_line(color = "#16a085", linewidth = 1) +
      labs(title = "VEGF-A serum", x = "Day", y = "VEGF (pg/mL)") + theme_bw()
  })

  output$plot_Anasarca <- renderPlot({
    ggplot(sim_data(), aes(time, Anasarca_lab)) +
      geom_line(color = "#1f618d", linewidth = 1) +
      labs(title = "Anasarca / third-space fluid score (TAFRO)",
           x = "Day", y = "Score 0–10") + theme_bw()
  })

  output$plot_Plt <- renderPlot({
    ggplot(sim_data(), aes(time, Plt_lab)) +
      geom_line(color = "#7d3c98", linewidth = 1) +
      geom_hline(yintercept = 150, linetype = "dashed", color = "red") +
      labs(title = "Platelet count", x = "Day", y = "Plt (×10⁹/L)") + theme_bw()
  })

  output$plot_LN <- renderPlot({
    ggplot(sim_data(), aes(time, LN_size)) +
      geom_line(color = "#d35400", linewidth = 1.2) +
      labs(title = "Lymph node composite size",
           x = "Day", y = "LN (cm sum)") + theme_bw()
  })

  output$plot_PB <- renderPlot({
    ggplot(sim_data(), aes(time, PB)) +
      geom_line(color = "#e67e22", linewidth = 1) +
      labs(title = "Plasmablast fraction in LN",
           x = "Day", y = "PB fraction") + theme_bw()
  })

  output$plot_Bmem <- renderPlot({
    ggplot(sim_data(), aes(time, Bmem)) +
      geom_line(color = "#2c3e50", linewidth = 1) +
      labs(title = "Memory B-cell pool (rituximab effect)",
           x = "Day", y = "% baseline") + theme_bw()
  })

  output$plot_CDCN <- renderPlot({
    ggplot(sim_data(), aes(time, CDCN_resp)) +
      geom_line(color = "#922b21", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "CDCN composite response (lower = better)",
           x = "Day", y = "Composite score") + theme_bw()
  })

  output$plot_surv <- renderPlot({
    ggplot(sim_data(), aes(time, Survival)) +
      geom_line(color = "#1a5276", linewidth = 1.2) +
      ylim(0, 1) +
      labs(title = "Modeled cumulative survival",
           x = "Day", y = "S(t)") + theme_bw()
  })

  output$biomarker_table <- DT::renderDT({
    d <- sim_data()
    final <- tail(d, 1)
    initial <- d[1, ]
    biomarkers <- data.frame(
      Endpoint = c("IL-6 free (pg/mL)", "CRP (mg/L)", "Hb (g/dL)",
                   "IgG polyclonal (g/dL)", "LN composite (cm)",
                   "VEGF (pg/mL)", "Anasarca (0-10)", "Platelet (×10⁹/L)",
                   "Survival probability", "CDCN composite"),
      Baseline = round(c(initial$IL6_free, initial$CRP_lab, initial$Hb_lab,
                         initial$IgG_lab, initial$LN_size, initial$VEGF_lab,
                         initial$Anasarca_lab, initial$Plt_lab,
                         initial$Survival, initial$CDCN_resp), 2),
      Endpoint_val = round(c(final$IL6_free, final$CRP_lab, final$Hb_lab,
                             final$IgG_lab, final$LN_size, final$VEGF_lab,
                             final$Anasarca_lab, final$Plt_lab,
                             final$Survival, final$CDCN_resp), 2)
    )
    DT::datatable(biomarkers, rownames = FALSE,
                  options = list(dom = "t", pageLength = 12))
  })

  # ----- Scenario comparison (uses static all_sim from source) ------
  output$plot_scenarios_CRP <- renderPlot({
    ggplot(all_sim, aes(time, CRP_lab, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = "CRP across pre-defined treatment scenarios",
           x = "Day", y = "CRP (mg/L)") + theme_bw()
  })

  output$plot_scenarios_LN <- renderPlot({
    ggplot(all_sim, aes(time, LN_size, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = "Lymph node composite — scenario comparison",
           x = "Day", y = "LN (cm)") + theme_bw()
  })

  output$plot_scenarios_response <- renderPlot({
    ggplot(all_sim, aes(time, CDCN_resp, color = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      labs(title = "CDCN composite response — lower is better",
           x = "Day", y = "Composite score") + theme_bw()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
