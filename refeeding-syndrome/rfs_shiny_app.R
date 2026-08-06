## ============================================================================
##  rfs_shiny_app.R
##  REFEEDING SYNDROME  -  interactive QSP dashboard
## ============================================================================
##
##  The app is built around the claim the model exists to make:
##
##      refeeding syndrome is a FLUX MISMATCH, Lambda_P = J_demand/J_supply,
##      in which the CLINICIAN sets the numerator (through the glucose
##      infusion rate) and the PATIENT'S HISTORY sets the denominator -
##      and the blood test everybody looks at samples 0.06 % of the pool
##      that is actually at risk.
##
##  So the left-hand panel is deliberately split in two.  The top half is
##  the HISTORY, which you cannot change and cannot see in the bloods; the
##  bottom half is the PRESCRIPTION, which is the only thing you control.
##  Every tab is an attempt to show what those two halves do to each other.
##
##  Run with:   shiny::runApp("rfs_shiny_app.R")
##  Requires:   rfs_mrgsolve_model.R in the same directory.
##
##  EDUCATIONAL / RESEARCH MODEL.  NOT FOR CLINICAL USE.
## ============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

source("rfs_mrgsolve_model.R", local = TRUE)

## ---- reference intervals, used everywhere ---------------------------------
REF <- list(
  PSER  = c(0.80, 1.45), KSER  = c(3.50, 5.10),
  MGSER = c(0.70, 1.00), CASER = c(1.15, 1.30)
)
SEVERE <- list(PSER = 0.32, KSER = 2.50, MGSER = 0.50)

theme_rfs <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

#' A time-course panel with the reference interval shaded behind it.
band_plot <- function(d, col, title, sub, ylab, ref = NULL, severe = NULL) {
  p <- ggplot(d, aes(time / 24, .data[[col]]))
  if (!is.null(ref))
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf,
                      ymin = ref[1], ymax = ref[2],
                      fill = "#2e8b57", alpha = 0.10)
  if (!is.null(severe))
    p <- p + geom_hline(yintercept = severe, linetype = "dashed",
                        colour = "#b03030", linewidth = 0.5)
  p + geom_line(linewidth = 1.05, colour = "#1d4a76") +
    labs(title = title, subtitle = sub, x = "days of refeeding", y = ylab) +
    theme_rfs()
}

## ============================================================================
##  UI
## ============================================================================

ui <- fluidPage(
  titlePanel("Refeeding Syndrome · QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         tags$em(paste("A flux model. The clinician sets the numerator of",
                       "Lambda_P; the starvation history sets the denominator.",
                       "Educational model - not for clinical use."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      tags$h4("1 · The history"),
      tags$p(style = "font-size:11px;color:#777",
             "You cannot change any of this, and none of it shows up in the",
             "admission blood panel. It sets J_supply."),
      sliderInput("starve_days", "Days of reduced intake", 7, 150, 60, 1),
      sliderInput("starve_frac", "Intake as a fraction of requirement",
                  0.05, 0.90, 0.32, 0.01),
      sliderInput("quality", "Micronutrient quality of that intake",
                  0, 1, 0.25, 0.05),
      sliderInput("alcohol",  "Alcohol use", 0, 1, 0.00, 0.05),
      sliderInput("diuretic", "Diuretic exposure", 0, 1, 0.00, 0.05),
      sliderInput("gi_loss",  "Vomiting / purging / diarrhoea", 0, 1, 0.20, 0.05),

      tags$hr(),
      tags$h4("2 · The prescription"),
      tags$p(style = "font-size:11px;color:#777",
             "This is the only half you control. It sets J_demand."),
      sliderInput("kcal_start", "Starting energy (kcal/kg/d)", 5, 40, 10, 1),
      sliderInput("kcal_goal",  "Target energy (kcal/kg/d)", 10, 45, 30, 1),
      sliderInput("adv_frac",   "Advance per day (fraction of goal)",
                  0, 0.60, 0.20, 0.01),
      sliderInput("cho_frac",   "Carbohydrate fraction of energy",
                  0.10, 1.00, 0.50, 0.05),
      radioButtons("route", "Route",
                   c("Enteral formula" = 0, "Intravenous dextrose" = 1),
                   selected = 0, inline = FALSE),

      tags$hr(),
      tags$h4("3 · Repletion"),
      sliderInput("p_dose",  "IV phosphate (mmol/kg/d)", 0, 2.0, 0.5, 0.05),
      sliderInput("k_dose",  "IV potassium (mmol/kg/d)", 0, 5.0, 2.5, 0.1),
      sliderInput("mg_dose", "IV magnesium (mmol/kg/d)", 0, 1.0, 0.3, 0.05),
      checkboxInput("mg_bolus", "Give the magnesium as a 2 h bolus", FALSE),
      sliderInput("th_iv", "IV thiamine (mg/d)", 0, 1000, 200, 25),
      sliderInput("th_po", "Oral thiamine (mg/d)", 0, 1500, 0, 50),
      sliderInput("th_delay", "Hours by which thiamine is started late",
                  0, 96, 0, 6),
      sliderInput("p_delay", "Days by which phosphate is started late",
                  0, 7, 0, 1),
      sliderInput("gut_fail", "Impaired absorption (vomiting, ileus)",
                  0, 1, 0, 0.05),
      checkboxInput("rescue", "Reactive rescue at the usual thresholds", TRUE),
      sliderInput("sim_days", "Days simulated", 3, 28, 14, 1),

      tags$hr(),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---------------------------------------------------------------
        tabPanel(
          "1 · Patient",
          br(),
          fluidRow(column(12, wellPanel(
            tags$h4("What two months of starvation did, and what the blood test shows"),
            tableOutput("adm_tbl")))),
          fluidRow(column(6, plotOutput("p_pools", height = 320)),
                   column(6, plotOutput("p_serum_adm", height = 320))),
          fluidRow(column(12, wellPanel(
            htmlOutput("adm_note"))))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "2 · Glucose & insulin",
          br(),
          fluidRow(column(6, plotOutput("p_glu", height = 280)),
                   column(6, plotOutput("p_ins", height = 280))),
          fluidRow(column(6, plotOutput("p_gir", height = 280)),
                   column(6, plotOutput("p_uins", height = 280))),
          helpText(paste("The glycolytic flux, not the calorie count and not",
                         "the insulin concentration, is what sets the cellular",
                         "phosphate demand. Intravenous dextrose bypasses the",
                         "incretin signal, so a gram of it raises insulin less",
                         "than a gram of enteral carbohydrate."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "3 · Lambda_P (the flux balance)",
          br(),
          fluidRow(column(12, plotOutput("p_lambda", height = 330))),
          fluidRow(column(6, plotOutput("p_pser", height = 300)),
                   column(6, plotOutput("p_picf", height = 300))),
          helpText(paste("Lambda_P above 1 means the cells are asking for more",
                         "phosphate per day than every supply line can deliver.",
                         "The buffer that absorbs the difference is the entire",
                         "extracellular pool - about 14 mmol, which is a few",
                         "hours of demand."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "4 · Electrolytes",
          br(),
          fluidRow(column(6, plotOutput("p_k", height = 280)),
                   column(6, plotOutput("p_mg", height = 280))),
          fluidRow(column(6, plotOutput("p_ca", height = 280)),
                   column(6, plotOutput("p_caxp", height = 280))),
          helpText(paste("Magnesium sits upstream of both of the others: it is",
                         "required for PTH secretion, and PTH is what makes",
                         "bone release phosphate; and its absence releases the",
                         "intracellular block on ROMK, so the kidney wastes",
                         "potassium."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "5 · Thiamine & the pyruvate gate",
          br(),
          fluidRow(column(6, plotOutput("p_tht", height = 280)),
                   column(6, plotOutput("p_tk", height = 280))),
          fluidRow(column(6, plotOutput("p_lac", height = 280)),
                   column(6, plotOutput("p_we", height = 280))),
          helpText(paste("The thiamine store empties in about three weeks and",
                         "refills in hours, because at supra-physiological",
                         "plasma levels it reaches tissue by passive diffusion.",
                         "Phosphate does the opposite. That asymmetry, not",
                         "tradition, is why the two are prescribed differently."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "6 · Organ energetics",
          br(),
          fluidRow(column(6, plotOutput("p_atp", height = 280)),
                   column(6, plotOutput("p_dpg", height = 280))),
          fluidRow(column(6, plotOutput("p_qtc", height = 280)),
                   column(6, plotOutput("p_oed", height = 280))),
          helpText(paste("A starved heart has lost mass on a three-week time",
                         "constant. Insulin then retains sodium. A volume load",
                         "delivered to a small ventricle with a low",
                         "phosphocreatine content is the cardiac-failure",
                         "pathway in this model."))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "7 · Clinical endpoints",
          br(),
          fluidRow(column(12, tableOutput("endpoints"))),
          fluidRow(column(6, plotOutput("p_mort", height = 300)),
                   column(6, plotOutput("p_kcal", height = 300)))
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "8 · Scenario comparison",
          br(),
          checkboxGroupInput(
            "scen", "Regimens to compare", names(rfs_scenarios),
            selected = names(rfs_scenarios)[c(1, 2, 3, 6, 8)], inline = FALSE),
          actionButton("go_scen", "Run selected", class = "btn-primary"),
          br(), br(),
          tableOutput("scen_tbl"),
          plotOutput("scen_plot", height = 340)
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "9 · Sweeps",
          br(),
          radioButtons("sweep", "Sweep",
                       c("Glucose infusion rate (energy held fixed)" = "gir",
                         "Starvation history (regimen held fixed)"   = "hist",
                         "Phosphate dose"                            = "pdose",
                         "Timing: the two clocks"                    = "time"),
                       inline = TRUE),
          actionButton("go_sweep", "Run sweep", class = "btn-primary"),
          br(), br(),
          plotOutput("sweep_plot", height = 360),
          tableOutput("sweep_tbl")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "10 · Model notes",
          br(),
          htmlOutput("notes")
        )
      )
    )
  )
)

## ============================================================================
##  SERVER
## ============================================================================

server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Starving, then refeeding...", value = 0.3, {
      rfs_run(
        starve_days = input$starve_days,
        sim_days    = input$sim_days,
        STARVE_FRAC = input$starve_frac,
        QUALITY     = input$quality,
        ALCOHOL     = input$alcohol,
        DIURETIC    = input$diuretic,
        GI_LOSS     = input$gi_loss,
        KCAL_START  = input$kcal_start,
        KCAL_GOAL   = input$kcal_goal,
        ADV_FRAC    = input$adv_frac,
        ADV_START   = 2,
        CHO_FRAC    = input$cho_frac,
        FAT_FRAC    = (1 - input$cho_frac) * 0.64,
        PRO_FRAC    = (1 - input$cho_frac) * 0.36,
        ROUTE       = as.numeric(input$route),
        P_DOSE      = input$p_dose,
        K_DOSE      = input$k_dose,
        MG_DOSE     = input$mg_dose,
        MG_BOLUS    = as.numeric(input$mg_bolus),
        TH_IV       = input$th_iv,
        TH_PO       = input$th_po,
        TH_DELAY    = input$th_delay,
        P_DELAY     = input$p_delay,
        GUT_FAIL    = input$gut_fail,
        RESCUE      = as.numeric(input$rescue)
      )
    })
  })

  d  <- reactive(sim()$refeed)
  ad <- reactive(sim()$admission)

  ## ---- tab 1 ---------------------------------------------------------------
  output$adm_tbl <- renderTable({
    a <- ad()
    data.frame(
      Quantity = c("Body mass index", "Fat-free mass (kg)",
                   "Serum phosphate (mmol/L)", "Serum potassium (mmol/L)",
                   "Serum magnesium (mmol/L)", "Ionised calcium (mmol/L)",
                   "Intracellular phosphate (% of normal)",
                   "Intracellular potassium (% of normal)",
                   "Exchangeable magnesium (% of normal)",
                   "Thiamine store (% of normal)"),
      Value = c(sprintf("%.1f", a$BMI), sprintf("%.1f", a$FFM),
                sprintf("%.2f", a$PSER), sprintf("%.2f", a$KSER),
                sprintf("%.2f", a$MGSER), sprintf("%.2f", a$CASER),
                sprintf("%.0f", 100 * a$PICF), sprintf("%.0f", 100 * a$KICF),
                sprintf("%.0f", 100 * a$MGICF),
                sprintf("%.0f", 100 * a$THT / 100)),
      Visible = c("yes", "yes", "MEASURED", "MEASURED", "MEASURED", "MEASURED",
                  "invisible", "invisible", "invisible", "invisible"))
  }, striped = TRUE, width = "100%")

  output$p_pools <- renderPlot({
    a <- ad()
    df <- data.frame(
      pool = factor(c("ECF (measured)", "intracellular", "bone"),
                    levels = c("ECF (measured)", "intracellular", "bone")),
      mmol = c(a$PE, a$PI, a$PB))
    ggplot(df, aes(pool, mmol, fill = pool)) +
      geom_col(width = 0.65) +
      scale_y_log10() +
      scale_fill_manual(values = c("#b03030", "#7fa8cc", "#a9a9a9"),
                        guide = "none") +
      labs(title = "The three phosphate pools, on a log scale",
           subtitle = sprintf(
             "The measured compartment is %.3f %% of the total",
             100 * a$PE / (a$PE + a$PI + a$PB)),
           x = NULL, y = "mmol (log scale)") +
      theme_rfs()
  })

  output$p_serum_adm <- renderPlot({
    a <- ad()
    df <- data.frame(
      analyte = c("phosphate", "potassium", "magnesium", "calcium"),
      value = c(a$PSER, a$KSER, a$MGSER, a$CASER),
      lo = c(REF$PSER[1], REF$KSER[1], REF$MGSER[1], REF$CASER[1]),
      hi = c(REF$PSER[2], REF$KSER[2], REF$MGSER[2], REF$CASER[2]))
    df$frac <- (df$value - df$lo) / (df$hi - df$lo)
    ggplot(df, aes(analyte, frac)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1,
               fill = "#2e8b57", alpha = 0.12) +
      geom_point(size = 5, colour = "#1d4a76") +
      geom_hline(yintercept = c(0, 1), linetype = "dotted") +
      labs(title = "The admission panel, as a position in the reference range",
           subtitle = "0 = lower limit, 1 = upper limit. This is all you get to see.",
           x = NULL, y = "position in reference interval") +
      theme_rfs()
  })

  output$adm_note <- renderUI({
    a <- ad()
    HTML(sprintf(
      "<p>After <b>%d days</b> at <b>%.0f %%</b> of energy requirement this
       patient has lost <b>%.0f %%</b> of their intracellular potassium and
       <b>%.0f %%</b> of their thiamine store, and their serum phosphate is
       <b>%.2f mmol/L</b> &mdash; inside the reference interval.<br><br>
       That is not a quirk of these slider settings. During starvation,
       catabolic release from shrinking cells and maximal renal conservation
       hold the serum concentration steady while the stores empty behind it.
       A normal admission phosphate in a starved patient is the
       <i>expected</i> finding, not a reassuring one, which is why the NICE
       criteria triage on history &mdash; body mass index, weight loss, days
       without intake, alcohol &mdash; rather than on this panel.</p>",
      input$starve_days, 100 * input$starve_frac,
      100 * (1 - a$KICF), 100 * (1 - a$THT / 100), a$PSER))
  })

  ## ---- tab 2 ---------------------------------------------------------------
  output$p_glu  <- renderPlot(band_plot(d(), "GLU", "Plasma glucose",
    "refeeding hyperglycaemia is transient: sensitivity recovers over a week",
    "mmol/L"))
  output$p_ins  <- renderPlot(band_plot(d(), "INS", "Plasma insulin",
    "glucose-driven plus an incretin term that only enteral feed switches on",
    "pmol/L"))
  output$p_gir  <- renderPlot({
    dd <- d(); dd$GIR <- dd$CUMKCAL  # placeholder replaced below
    ggplot(d(), aes(time / 24, X)) + geom_line(linewidth = 1.05, colour = "#b03878") +
      labs(title = "Insulin effect compartment X",
           subtitle = "the signal, which is NOT what drives the phosphate demand",
           x = "days of refeeding", y = "X (1 = well-fed reference)") + theme_rfs()
  })
  output$p_uins <- renderPlot({
    ggplot(d(), aes(time / 24, LAC)) +
      geom_line(linewidth = 1.05, colour = "#c05a2a") +
      labs(title = "Lactate",
           subtitle = "pyruvate that cannot pass a TPP-depleted PDH gate",
           x = "days of refeeding", y = "mmol/L") + theme_rfs()
  })

  ## ---- tab 3 ---------------------------------------------------------------
  output$p_lambda <- renderPlot({
    dd <- d()
    ggplot(dd, aes(time / 24, PSER)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = REF$PSER[1],
               ymax = REF$PSER[2], fill = "#2e8b57", alpha = 0.10) +
      geom_hline(yintercept = 0.65, linetype = "dashed", colour = "#c08a1e") +
      geom_hline(yintercept = SEVERE$PSER, linetype = "dashed",
                 colour = "#b03030") +
      geom_line(linewidth = 1.2, colour = "#1d4a76") +
      annotate("text", x = max(dd$time) / 24 * 0.02, y = 0.67, hjust = 0,
               label = "0.65 = refeeding hypophosphataemia", size = 3.2,
               colour = "#8a6114") +
      annotate("text", x = max(dd$time) / 24 * 0.02, y = 0.34, hjust = 0,
               label = "0.32 = severe", size = 3.2, colour = "#b03030") +
      labs(title = "Serum phosphate through the refeed",
           subtitle = "the whole disease, read off a compartment that holds 14 mmol",
           x = "days of refeeding", y = "mmol/L") + theme_rfs()
  })
  output$p_pser <- renderPlot(band_plot(d(), "TBP",
    "Total body phosphorus", "ECF + intracellular + bone", "mmol"))
  output$p_picf <- renderPlot(band_plot(d(), "PICF",
    "Intracellular phosphate density",
    "fraction of the normal 67 mmol per kg fat-free mass", "fraction"))

  ## ---- tab 4 ---------------------------------------------------------------
  output$p_k  <- renderPlot(band_plot(d(), "KSER", "Serum potassium",
    "insulin drives it into cells; hypomagnesaemia makes the kidney waste it",
    "mmol/L", REF$KSER, SEVERE$KSER))
  output$p_mg <- renderPlot(band_plot(d(), "MGSER", "Serum magnesium",
    "renal handling is a threshold, so a bolus is largely excreted",
    "mmol/L", REF$MGSER, SEVERE$MGSER))
  output$p_ca <- renderPlot(band_plot(d(), "CASER", "Ionised calcium",
    "falls when magnesium is too low for PTH secretion", "mmol/L", REF$CASER))
  output$p_caxp <- renderPlot({
    ggplot(d(), aes(time / 24, CAXP)) +
      geom_hline(yintercept = 4.4, linetype = "dashed", colour = "#b03030") +
      geom_line(linewidth = 1.05, colour = "#7a5aa0") +
      labs(title = "Calcium × phosphate product",
           subtitle = "above 4.4 mmol²/L² the model precipitates it",
           x = "days of refeeding", y = "mmol²/L²") + theme_rfs()
  })

  ## ---- tab 5 ---------------------------------------------------------------
  output$p_tht <- renderPlot(band_plot(d(), "THT", "Tissue thiamine store",
    "26.5 mg when full; a 14 d half-life predicts the 1.1-1.4 mg/d RDA", "umol"))
  output$p_tk  <- renderPlot(band_plot(d(), "TKACT",
    "Transketolase activity", "the TPP-dependent enzyme the assay reports",
    "fraction of normal"))
  output$p_lac <- renderPlot(band_plot(d(), "LAC", "Lactate",
    "type B acidosis: glucose load meeting a blocked pyruvate gate", "mmol/L"))
  output$p_we  <- renderPlot(band_plot(d(), "WERN",
    "Cumulative Wernicke incidence", "hazard scales with the glucose load", "%"))

  ## ---- tab 6 ---------------------------------------------------------------
  output$p_atp <- renderPlot(band_plot(d(), "ATPM", "Myocardial ATP index",
    "limited by BOTH inorganic phosphate and TPP", "fraction of normal"))
  output$p_dpg <- renderPlot(band_plot(d(), "DPG", "Erythrocyte 2,3-DPG",
    "low DPG left-shifts the oxygen curve: tissue hypoxia at a normal SaO2",
    "fraction of normal"))
  output$p_qtc <- renderPlot(band_plot(d(), "QTC", "QTc",
    "the arrhythmic consequence of the potassium, magnesium and calcium",
    "ms") + geom_hline(yintercept = 500, linetype = "dashed", colour = "#b03030"))
  output$p_oed <- renderPlot(band_plot(d(), "OEDEMA", "Refeeding oedema",
    "insulin-driven renal sodium retention loading an atrophic ventricle", "L"))

  ## ---- tab 7 ---------------------------------------------------------------
  output$endpoints <- renderTable({
    s <- rfs_summary(sim())
    data.frame(Endpoint = names(s), Value = sprintf("%.3g", unlist(s)))
  }, striped = TRUE, width = "100%")
  output$p_mort <- renderPlot(band_plot(d(), "MORT", "Cumulative mortality",
    "arrhythmia + cardiac failure + hypercapnic respiratory failure", "%"))
  output$p_kcal <- renderPlot(band_plot(d(), "CUMKCAL",
    "Cumulative energy delivered",
    "the benefit that caloric restriction trades away", "kcal"))

  ## ---- tab 8 ---------------------------------------------------------------
  scen <- eventReactive(input$go_scen, {
    withProgress(message = "Running scenarios...", value = 0.3, {
      do.call(rbind, lapply(input$scen, function(nm) {
        r <- do.call(rfs_run, c(list(starve_days = input$starve_days,
                                     sim_days = input$sim_days),
                                rfs_scenarios[[nm]]))
        cbind(scenario = nm, rfs_summary(r), row.names = NULL)
      }))
    })
  })
  output$scen_tbl  <- renderTable(scen(), digits = 3, striped = TRUE)
  output$scen_plot <- renderPlot({
    s <- scen()
    s %>%
      select(scenario, P_nadir, K_d3, Mg_d3, mortality) %>%
      pivot_longer(-scenario) %>%
      ggplot(aes(reorder(scenario, value), value, fill = name)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x") +
      labs(x = NULL, y = NULL, title = "Regimens compared in one patient") +
      theme_rfs()
  })

  ## ---- tab 9 ---------------------------------------------------------------
  sw <- eventReactive(input$go_sweep, {
    withProgress(message = "Sweeping...", value = 0.3, {
      switch(input$sweep,
             gir   = rfs_sweep_gir(),
             hist  = rfs_sweep_history(),
             pdose = rfs_sweep_pdose(),
             time  = rfs_sweep_timing())
    })
  })
  output$sweep_tbl <- renderTable(sw(), digits = 3, striped = TRUE)
  output$sweep_plot <- renderPlot({
    s <- sw()
    xv <- switch(input$sweep, gir = "CHO_pct", hist = "days",
                 pdose = "P_dose", time = "delay_h")
    ggplot(s, aes(.data[[xv]], P_nadir)) +
      geom_line(linewidth = 1.1, colour = "#1d4a76") +
      geom_point(size = 2.4, colour = "#1d4a76") +
      geom_hline(yintercept = 0.65, linetype = "dashed", colour = "#c08a1e") +
      labs(title = "Phosphate nadir across the sweep",
           x = xv, y = "nadir (mmol/L)") + theme_rfs()
  })

  ## ---- tab 10 --------------------------------------------------------------
  output$notes <- renderUI(HTML('
    <h4>What this model claims</h4>
    <p>Refeeding syndrome is written here as a <b>flux mismatch</b>,
       &Lambda;<sub>P</sub> = J<sub>demand</sub>/J<sub>supply</sub>, rather
       than as a low number on a blood test. The numerator is set by the
       clinician through the glucose infusion rate; the denominator is set by
       the patient&rsquo;s history. Three pieces of arithmetic make the ratio
       dangerous, and none of them is a fitted parameter:</p>
    <ol>
      <li>extracellular phosphate is 14 mmol, which is <b>0.06 %</b> of total
          body phosphorus, so a 60 mmol/d demand is four times the entire
          measured compartment per day;</li>
      <li>urinary phosphate excretion is a <b>threshold</b>, so the kidney&rsquo;s
          only contribution is to stop excreting &mdash; worth about 14 mmol/d,
          spent within a day, and then gone;</li>
      <li>bone efflux is slow and PTH-driven, and <b>hypomagnesaemia switches
          PTH off</b>, so the last endogenous supply line is disabled by a
          second electrolyte the same syndrome depletes.</li>
    </ol>

    <h4>What was calibrated on refeeding syndrome</h4>
    <p><b>Three numbers.</b> The size of the flux-driven rise in the cellular
       phosphate set-point, the rate at which cells approach it, and one
       global mortality-hazard scale. Everything else comes from normal
       physiology or from non-refeeding pharmacology. As a check that the
       normal-physiology layer is not being quietly fitted: the whole-body
       thiamine store and its 14-day half-life were taken from tracer studies,
       and the model then <i>computes</i> a dietary requirement of 1.42 mg/d
       against a published RDA of 1.1&ndash;1.4 mg/d.</p>

    <h4>Findings that ran against expectation, reported as such</h4>
    <ul>
      <li><b>Oral thiamine is not the weak link it is usually said to be.</b>
          With an intact gut, oral and intravenous converge within a day,
          because saturated absorption still delivers about thirty times the
          requirement. The intravenous route earns its place when absorption
          is impaired &mdash; which is the same population, so the
          recommendation survives with its justification replaced.</li>
      <li><b>The phosphate nadir is not monotone in starvation duration.</b> It
          deepens out to about 60&ndash;75 days and then becomes shallower,
          because the most wasted patients cannot mount the glycolytic flux
          that generates the demand. Their mortality still rises. A shallow
          nadir in a very wasted patient is not reassurance.</li>
      <li><b>The U-shaped phosphate dose&ndash;response is real but clinically
          irrelevant.</b> Benefit saturates by 0.5 mmol/kg/d and the
          calcium&ndash;phosphate product only bites at twenty times that.
          The model has no soft-tissue calcification and fixed renal function,
          so read the flat top as &ldquo;this model cannot see the harm&rdquo;,
          not as &ldquo;there is none&rdquo;.</li>
      <li><b>Restricting calories without repleting electrolytes buys little</b>
          and prolongs the thiamine-deficient window. In this model the active
          ingredient in the NICE/ASPEN disagreement is the repletion, not the
          calorie cap.</li>
    </ul>

    <h4>Verification</h4>
    <p>Every equation was independently re-implemented in Python/scipy. That
       exercise found five real defects &mdash; a set-point referenced to the
       wrong baseline so the disease mechanism was inert; the same set-point
       then driven by insulin, which has the wrong sign; renal potassium
       excretion with no depletion adaptation; an unbounded PaCO<sub>2</sub>
       inside an exponential hazard; and an insulin secretion curve that was
       three-quarters glucose-independent basal. All five are fixed. The
       healthy state is a numerically exact steady state
       (max |dy/dt| = 6&times;10<sup>-15</sup>).</p>

    <h4>Limitations</h4>
    <p>No sepsis, organ failure, delirium, hepatic steatosis or vitamin
       deficiency other than thiamine. Fixed renal function. Bone is one
       well-mixed compartment, so the skeletal buffer is probably overstated
       in the most chronic cases. The worst scenario over-shoots: real nadirs
       on phosphate-free dextrose are 0.10&ndash;0.30 mmol/L and the model
       reaches 0.07. Mortality carries a single fitted scale and its split
       across three mechanisms is structural, not separately validated.</p>

    <p style="color:#b03030"><b>Educational and research model. Not validated
       for clinical decision-making, prescribing, or regulatory use.</b></p>'))
}

shinyApp(ui, server)
