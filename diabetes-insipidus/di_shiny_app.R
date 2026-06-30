## =============================================================================
##  Shiny dashboard for the Diabetes Insipidus (AVP-D / AVP-R) QSP model
##  8 tabs: Patient profile · Disease severity · Desmopressin PK · AQP2/V2R PD ·
##          Electrolytes & osmolality · Urine & thirst endpoints · Scenarios ·
##          Biomarkers & safety
##  Run from the directory root:
##    shiny::runApp("diabetes-insipidus/di_shiny_app.R")
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
})

source_local <- function() {
  here <- "di_mrgsolve_model.R"
  if (!file.exists(here)) here <- file.path("diabetes-insipidus", here)
  source(here, local = TRUE)
}
source_local()

# Convenient palette
pal <- c("#1c7ed6", "#c92a2a", "#2f9e44", "#e67700", "#6741d9", "#0c8599",
         "#fd7e14", "#c2255c", "#495057")

scenarios_full <- c(
  "untreated_CDI"               = "Untreated complete CDI (AVP-D)",
  "DDAVP_SC_2ug_BID"            = "DDAVP 2 µg SC q12h",
  "DDAVP_IN_10ug_BID"           = "DDAVP 10 µg intranasal q12h",
  "DDAVP_PO_200ug_TID"          = "DDAVP 200 µg PO q8h",
  "DDAVP_SL_120ug_TID"          = "DDAVP 120 µg sublingual lyophilisate q8h",
  "NDI_lithium_HCTZ"            = "Lithium NDI + HCTZ 25 mg/d",
  "NDI_lithium_amiloride"       = "Lithium NDI + amiloride 10 mg/d",
  "NDI_indomethacin"            = "NDI + indomethacin 50 mg q8h",
  "tolvaptan_SIADH_comparator"  = "Tolvaptan 15 mg/d (V2 antagonist comparator)",
  "primary_polydipsia"          = "Primary polydipsia (+6 L/d psych drive)",
  "gestational_DDAVP"           = "Gestational DI + DDAVP IN",
  "pediatric_DDAVP_SC"          = "Pediatric (20 kg) DDAVP 0.3 µg SC q12h"
)

# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Diabetes Insipidus (AVP-D / AVP-R) — QSP Simulator"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("Patient profile"),
      numericInput("WT", "Body weight (kg)", 70, 5, 150, 1),
      sliderInput("AGE", "Age (yr)", 0, 90, 40),
      radioButtons("SEX", "Sex", choices = c("Male" = 1, "Female" = 0), inline = TRUE),
      numericInput("eGFR", "eGFR (mL/min/1.73 m²)", 95, 10, 130, 1),
      tags$hr(),

      tags$h4("Disease phenotype"),
      sliderInput("AVP_D_sev", "AVP-D severity (central DI)",
                  0, 1, 1.0, step = 0.05),
      sliderInput("AVP_R_sev", "AVP-R severity (nephrogenic DI)",
                  0, 1, 0.0, step = 0.05),
      sliderInput("PSYCH_drive", "Primary polydipsia drive (L/d)",
                  0, 10, 0, step = 0.5),
      radioButtons("GEST", "Pregnancy (vasopressinase)",
                   choices = c("No" = 0, "Yes" = 1), inline = TRUE),
      tags$hr(),

      tags$h4("Therapy scenario"),
      selectInput("scen", "Pre-built scenario", choices = scenarios_full,
                  selected = "DDAVP_SC_2ug_BID"),
      numericInput("dur", "Simulation duration (h)", 96, 24, 720, 12),
      tags$hr(),
      actionButton("run", "Run simulation", class = "btn-primary",
                   icon = icon("play")),
      tags$br(), tags$br(),
      tags$small("Calibrated against Robertson 1976, Vavra 1968, Christ-Crain 2019 NEJM, Bichet 2019 NRDP and Bedford 2008 JASN.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient profile",
                 tags$h3("Demographic & physiological summary"),
                 tableOutput("profileTbl"),
                 tags$h4("Disease classification"),
                 tableOutput("dxTbl"),
                 tags$p("Mapping: severity 1.0 ≈ complete deficiency/resistance; 0.5 ≈ partial. ",
                        "Primary polydipsia adds drinking drive without AVP defect.")),

        tabPanel("2. Disease severity",
                 tags$h3("AVP-osmolality response curve"),
                 plotOutput("avpOsmPlot", height = "320px"),
                 tags$h3("Predicted basal labs"),
                 tableOutput("baselineTbl")),

        tabPanel("3. Desmopressin PK",
                 tags$h3("Plasma desmopressin (pg/mL)"),
                 plotOutput("ddavpPK", height = "360px"),
                 tags$p("F: SC 100% · IN ~4% · PO ~0.16% · SL ~0.25% (Vavra 1968, ",
                        "Lottermoser 1997, Steiner 2007).")),

        tabPanel("4. AQP2 / V2R PD",
                 tags$h3("V2 occupancy & AQP2 apical fraction"),
                 plotOutput("aqp2plot", height = "320px"),
                 tags$h3("Urinary osmolality"),
                 plotOutput("uosmPlot", height = "260px")),

        tabPanel("5. Electrolytes & osmolality",
                 tags$h3("Plasma Na+ and osmolality"),
                 plotOutput("naPlot", height = "300px"),
                 tags$h3("Total body water"),
                 plotOutput("tbwPlot", height = "260px")),

        tabPanel("6. Urine & thirst endpoints",
                 tags$h3("Daily urine output (L/d, instantaneous)"),
                 plotOutput("uflowPlot", height = "300px"),
                 tags$h3("Cumulative urine vs intake"),
                 plotOutput("cumPlot", height = "260px"),
                 tags$h3("Thirst (mOsm/kg above threshold)"),
                 plotOutput("thirstPlot", height = "260px")),

        tabPanel("7. Scenario comparison",
                 tags$h3("Compare untreated vs treated"),
                 plotOutput("scenCompare", height = "420px"),
                 tags$p("Solid: current scenario. Dashed: untreated AVP-D reference.")),

        tabPanel("8. Biomarkers & safety",
                 tags$h3("Copeptin trajectory"),
                 plotOutput("copPlot", height = "260px"),
                 tags$h3("Hyponatremia hazard accumulation"),
                 plotOutput("hazPlot", height = "260px"),
                 tags$p("Hazard accumulates when plasma Na+ < 130 mmol/L — surrogate for over-DDAVP."))
      )
    )
  )
)

# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  simulate <- eventReactive(input$run, {
    mod <- param(di_mod,
                 WT       = input$WT,
                 AGE      = input$AGE,
                 SEX      = as.numeric(input$SEX),
                 eGFR     = input$eGFR,
                 AVP_D_sev = input$AVP_D_sev,
                 AVP_R_sev = input$AVP_R_sev,
                 PSYCH_drive = input$PSYCH_drive,
                 GEST = as.numeric(input$GEST))
    ev_scen <- run_scenario(input$scen, duration_h = input$dur)
    list(custom = ev_scen, dur = input$dur)
  }, ignoreNULL = FALSE)

  # ------ Tab 1: profile -----------------------------------------------------
  output$profileTbl <- renderTable({
    tibble::tibble(
      Variable = c("Weight (kg)", "Age", "Sex", "eGFR"),
      Value    = c(input$WT, input$AGE,
                   ifelse(input$SEX == 1, "Male", "Female"), input$eGFR)
    )
  })
  output$dxTbl <- renderTable({
    typ <- if (input$AVP_D_sev > 0.5 && input$AVP_R_sev < 0.3) "AVP-D (central DI)"
           else if (input$AVP_R_sev > 0.5) "AVP-R (nephrogenic DI)"
           else if (input$PSYCH_drive > 3) "Primary polydipsia"
           else "Mixed / partial"
    tibble::tibble(
      Phenotype     = typ,
      `AVP-D sev`   = input$AVP_D_sev,
      `AVP-R sev`   = input$AVP_R_sev,
      `Polydipsia drive (L/d)` = input$PSYCH_drive,
      Gestational   = ifelse(input$GEST == 1, "Yes", "No")
    )
  })

  # ------ Tab 2: AVP-osmolality curve ----------------------------------------
  output$avpOsmPlot <- renderPlot({
    Posm <- seq(265, 310, 1)
    avp <- pmax(0, (Posm - 280) * 0.4) * (1 - input$AVP_D_sev)
    df <- data.frame(Posm = Posm, AVP = avp)
    ggplot(df, aes(Posm, AVP)) +
      geom_line(color = pal[1], linewidth = 1.2) +
      geom_vline(xintercept = 280, lty = 2, color = pal[2]) +
      geom_vline(xintercept = 290, lty = 2, color = pal[3]) +
      annotate("text", 281, max(avp) * 0.9, label = "AVP threshold", hjust = 0) +
      annotate("text", 291, max(avp) * 0.75, label = "Thirst threshold", hjust = 0) +
      labs(x = "Plasma osmolality (mOsm/kg)",
           y = "Plasma AVP (pmol/L)",
           title = "Endogenous AVP-osmolality response (Robertson slope)") +
      theme_minimal(base_size = 13)
  })

  output$baselineTbl <- renderTable({
    sim <- simulate()
    s   <- as.data.frame(sim$custom$sim)
    last <- tail(s, 1)
    tibble::tibble(
      `Plasma Na+ (mmol/L)`   = round(last$Plasma_Na, 1),
      `Plasma Osm (mOsm/kg)`  = round(last$Plasma_Osm, 0),
      `Urine osm (mOsm/kg)`   = round(last$Uosm_mosm, 0),
      `Urine output (L/d)`    = round(last$Urine_Lday, 2),
      `Copeptin (pmol/L)`     = round(last$Copeptin_pmol, 2),
      `V2 occupancy`          = round(last$V2_occupancy, 3),
      `Apical AQP2`           = round(last$AQP2_apical, 3)
    )
  })

  # ------ Tab 3: PK ----------------------------------------------------------
  output$ddavpPK <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, DDAVP_conc_pgmL)) +
      geom_line(color = pal[1], linewidth = 1) +
      labs(x = "Time (h)", y = "Plasma desmopressin (pg/mL)",
           title = paste("PK trajectory —", scenarios_full[input$scen])) +
      theme_minimal(base_size = 13)
  })

  # ------ Tab 4: AQP2 --------------------------------------------------------
  output$aqp2plot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    df <- tidyr::pivot_longer(s, c("V2_occupancy", "AQP2_apical"),
                              names_to = "var", values_to = "value")
    ggplot(df, aes(time, value, color = var)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(pal[1], pal[3])) +
      labs(x = "Time (h)", y = "Fraction (0-1)",
           color = NULL, title = "V2R occupancy & AQP2 apical fraction") +
      theme_minimal(base_size = 13)
  })

  output$uosmPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, Uosm_mosm)) +
      geom_line(color = pal[6], linewidth = 1) +
      geom_hline(yintercept = 300, lty = 2, color = pal[2]) +
      labs(x = "Time (h)", y = "Urinary osmolality (mOsm/kg)",
           title = "Concentrating capacity") +
      theme_minimal(base_size = 13)
  })

  # ------ Tab 5: Electrolytes ------------------------------------------------
  output$naPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, Plasma_Na)) +
      geom_line(color = pal[2], linewidth = 1) +
      geom_hline(yintercept = 135, lty = 2, color = "grey50") +
      geom_hline(yintercept = 145, lty = 2, color = "grey50") +
      labs(x = "Time (h)", y = "Plasma Na+ (mmol/L)",
           title = "Sodium trajectory") +
      theme_minimal(base_size = 13)
  })
  output$tbwPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, TBW)) +
      geom_line(color = pal[5], linewidth = 1) +
      labs(x = "Time (h)", y = "Total body water (L)") +
      theme_minimal(base_size = 13)
  })

  # ------ Tab 6: Urine -------------------------------------------------------
  output$uflowPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, Urine_Lday)) +
      geom_line(color = pal[4], linewidth = 1) +
      geom_hline(yintercept = 3, lty = 2, color = "grey50") +
      labs(x = "Time (h)", y = "Urine flow (L/d, instantaneous)",
           title = "Polyuria control") +
      theme_minimal(base_size = 13)
  })
  output$cumPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, CUM_URINE)) +
      geom_line(color = pal[7], linewidth = 1) +
      labs(x = "Time (h)", y = "Cumulative urine (L)") +
      theme_minimal(base_size = 13)
  })
  output$thirstPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, Thirst_score)) +
      geom_line(color = pal[8], linewidth = 1) +
      labs(x = "Time (h)", y = "Thirst (mOsm/kg above threshold)") +
      theme_minimal(base_size = 13)
  })

  # ------ Tab 7: Scenario comparison -----------------------------------------
  output$scenCompare <- renderPlot({
    s    <- as.data.frame(simulate()$custom$sim)
    ref  <- as.data.frame(run_scenario("untreated_CDI",
                                       duration_h = input$dur)$sim)
    df <- rbind(
      transform(s,   group = "Current scenario"),
      transform(ref, group = "Untreated CDI")
    )
    ggplot(df, aes(time, Urine_Lday, color = group, linetype = group)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(pal[1], pal[2])) +
      labs(x = "Time (h)", y = "Urine flow (L/d)",
           color = NULL, linetype = NULL,
           title = "Urine flow: current scenario vs untreated reference") +
      theme_minimal(base_size = 13)
  })

  # ------ Tab 8: Biomarkers / safety -----------------------------------------
  output$copPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, Copeptin_pmol)) +
      geom_line(color = pal[5], linewidth = 1) +
      geom_hline(yintercept = 4.9, lty = 2, color = pal[2]) +
      labs(x = "Time (h)", y = "Plasma copeptin (pmol/L)",
           title = "Stimulated copeptin — Christ-Crain cut-off 4.9") +
      theme_minimal(base_size = 13)
  })
  output$hazPlot <- renderPlot({
    s <- as.data.frame(simulate()$custom$sim)
    ggplot(s, aes(time, CUM_HAZ)) +
      geom_line(color = pal[2], linewidth = 1) +
      labs(x = "Time (h)", y = "Cumulative hyponatremia hazard",
           title = "Iatrogenic hyponatremia risk integral") +
      theme_minimal(base_size = 13)
  })
}

shinyApp(ui, server)
