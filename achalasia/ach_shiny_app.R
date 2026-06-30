# =============================================================================
# Achalasia QSP — Shiny dashboard
# Tabs:
#   1. Patient profile (covariates + subtype)
#   2. Drug PK (ISDN, Nifedipine, Sildenafil, Tadalafil, Botox)
#   3. LES/IRP physiology (manometric endpoints)
#   4. Clinical endpoints (Eckardt components, TBE, QoL)
#   5. Scenario comparison (10 arms: drugs/procedures/combo)
#   6. Procedure outcomes (PD vs LHM vs POEM long-term)
#   7. Safety / Adverse events
#   8. Calibration & trial anchors
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(DT)
})

source("ach_mrgsolve_model.R", local = TRUE)

scn_choices <- c(
  "No Tx"                          = "S01_NoTx",
  "Isosorbide dinitrate 10 mg SL TID" = "S02_ISDN",
  "Nifedipine 20 mg SL TID"        = "S03_Nifedipine",
  "Sildenafil 50 mg TID"           = "S04_Sildenafil",
  "Tadalafil 20 mg QD"             = "S05_Tadalafil",
  "Botox 100 U intra-LES q6 mo"    = "S06_Botox",
  "Pneumatic dilation 35 mm"       = "S07_PneumaticDilation",
  "Heller myotomy + Dor"           = "S08_HellerMyotomy",
  "POEM"                           = "S09_POEM",
  "ISDN + Sildenafil combo"        = "S10_ISDN_plus_Sildenafil"
)

ui <- fluidPage(
  titlePanel("Achalasia (식도이완불능증) — QSP Dashboard"),
  sidebarLayout(
    sidebarPanel(width = 3,
      h4("Patient profile"),
      sliderInput("AGE", "Age (y)", 18, 90, 50),
      selectInput("SEX", "Sex", c("Female" = 0, "Male" = 1)),
      sliderInput("BWT", "Body weight (kg)", 40, 130, 70),
      selectInput("SUBTYPE", "Chicago subtype",
                  c("Type I (classic)" = 1, "Type II (panesophageal)" = 2,
                    "Type III (spastic)" = 3)),
      sliderInput("DURATION", "Disease duration (y)", 0, 30, 3),
      selectInput("CYP3A4", "CYP3A4 phenotype",
                  c("PM 0.5×" = 0.5, "EM 1×" = 1, "UM 2×" = 2),
                  selected = 1),
      hr(),
      h4("Scenario"),
      selectInput("scenario", "Treatment", scn_choices,
                  selected = "S09_POEM"),
      sliderInput("horizon", "Simulation horizon (months)", 1, 60, 24),
      hr(),
      h4("Multi-scenario comparison"),
      checkboxGroupInput("multi", "Pick arms (Scenario tab)",
                         choices = scn_choices,
                         selected = c("S01_NoTx", "S03_Nifedipine",
                                      "S06_Botox", "S07_PneumaticDilation",
                                      "S09_POEM")),
      actionButton("run", "Run / Re-run", class = "btn-primary")
    ),
    mainPanel(width = 9,
      tabsetPanel(
        tabPanel("Patient",
                 h3("Predicted untreated trajectory"),
                 plotOutput("p_traj_untreat", height = 420),
                 verbatimTextOutput("baseline_info")),
        tabPanel("Drug PK",
                 h3("Plasma concentration-time"),
                 plotOutput("p_pk", height = 420),
                 helpText("ISDN+ISMN, Nifedipine, Sildenafil, Tadalafil tracks; ",
                          "Botox shown as local LES units on log scale.")),
        tabPanel("LES / IRP",
                 h3("LES pressure and IRP over time"),
                 plotOutput("p_les", height = 420),
                 h4("Subtype-stratified comparison"),
                 plotOutput("p_les_subtype", height = 320)),
        tabPanel("Clinical endpoints",
                 h3("Eckardt symptom score & components"),
                 plotOutput("p_eckardt", height = 380),
                 h4("TBE column height and QoL"),
                 plotOutput("p_tbe_qol", height = 320)),
        tabPanel("Scenario comparison",
                 h3("Eckardt score across selected arms"),
                 plotOutput("p_multi_eck", height = 380),
                 h4("LES pressure"),
                 plotOutput("p_multi_les", height = 320)),
        tabPanel("Procedure outcomes",
                 h3("PD vs LHM vs POEM 2-year remission"),
                 plotOutput("p_proc", height = 360),
                 DTOutput("tbl_proc")),
        tabPanel("Safety",
                 h3("Drug headache & postural hypotension"),
                 plotOutput("p_ae", height = 360),
                 helpText("Headache from nitrate / PDE5 vasodilation; hypotension ",
                          "from nifedipine + nitrate. Botox AEs (heartburn, chest ",
                          "pain) are modeled separately as procedural flags.")),
        tabPanel("Calibration",
                 h3("Trial calibration table"),
                 DTOutput("tbl_calib"))
      )
    )
  )
)

server <- function(input, output, session) {
  patient_params <- reactive({
    list(AGE = as.numeric(input$AGE), SEX = as.numeric(input$SEX),
         BWT = as.numeric(input$BWT), SUBTYPE = as.numeric(input$SUBTYPE),
         DURATION = as.numeric(input$DURATION), CYP3A4 = as.numeric(input$CYP3A4))
  })

  run_scn <- function(scn_name, horizon_d = 720, subtype = 2,
                      extra_params = list()) {
    scn <- scenarios[[which(sapply(scenarios, function(x) x$name) == scn_name)]]
    m   <- mod_ach %>% param(c(extra_params, list(SUBTYPE = subtype)))
    if (!is.null(scn$proc)) {
      m <- apply_procedure(m, scn$proc, subtype)
      out <- mrgsim(m, end = horizon_d, delta = 1)
    } else {
      out <- mrgsim(m, events = scn$regimen, end = horizon_d, delta = 1)
    }
    df <- as.data.frame(out)
    df$scenario <- scn_name
    df$subtype  <- subtype
    df
  }

  active <- eventReactive(input$run, {
    pp <- patient_params()
    horizon_d <- input$horizon * 30
    list(main  = run_scn(input$scenario, horizon_d, pp$SUBTYPE, pp),
         multi = do.call(rbind, lapply(input$multi, run_scn,
                                       horizon_d = horizon_d,
                                       subtype = pp$SUBTYPE,
                                       extra_params = pp)),
         subt  = do.call(rbind, lapply(1:3, function(s)
                          run_scn(input$scenario, horizon_d, s, pp))),
         pp = pp, horizon_d = horizon_d)
  }, ignoreNULL = FALSE)

  output$p_traj_untreat <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(LES_PRESS, IRP_C, ESO_DIL, ECKARDT),
                   names_to = "var", values_to = "value") %>%
      ggplot(aes(time / 30, value, colour = var)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Months from baseline", y = "Value",
           title = paste("Trajectory —", input$scenario)) +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  output$baseline_info <- renderPrint({
    a <- active(); pp <- a$pp
    cat("Predicted baseline (untreated achalasia, Chicago subtype",
        pp$SUBTYPE, ", duration", pp$DURATION, "y)\n")
    cat("  LES pressure   :", round(40 + 2 * pp$DURATION, 1), "mmHg (target <15 after Tx)\n")
    cat("  IRP            :", round(25 + 0.5 * pp$DURATION, 1), "mmHg (target <15)\n")
    cat("  Esophageal Ø   :", round(3.0 + 0.2 * pp$DURATION, 2), "cm\n")
    cat("  Eckardt        :", "~6-9 (untreated), success <=3\n")
  })

  output$p_pk <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(Cp_ISDN, Cp_NIF, Cp_SIL, Cp_TAD, BOTOX_LES),
                   names_to = "drug", values_to = "Cp") %>%
      ggplot(aes(time / 24, pmax(Cp, 1e-3), colour = drug)) +
      geom_line(linewidth = 0.7) +
      scale_y_log10() +
      labs(x = "Days", y = "Concentration (ng/mL) or Botox U (local)",
           title = "Drug exposure tracks") +
      theme_minimal(base_size = 12)
  })

  output$p_les <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(LES_PRESS, IRP_C), names_to = "var", values_to = "value") %>%
      ggplot(aes(time / 30, value, colour = var)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 15, lty = 2, colour = "grey50") +
      labs(x = "Months", y = "Pressure (mmHg)",
           title = "LES resting pressure & IRP") +
      theme_minimal(base_size = 12)
  })

  output$p_les_subtype <- renderPlot({
    a <- active(); df <- a$subt
    df$Subtype <- paste0("Type ", df$subtype)
    ggplot(df, aes(time / 30, LES_PRESS, colour = Subtype)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Months", y = "LES pressure (mmHg)",
           title = "Subtype response to selected therapy") +
      theme_minimal(base_size = 12)
  })

  output$p_eckardt <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(DYS_S, REG_S, CP_S, WT_S, ECKARDT),
                   names_to = "var", values_to = "value") %>%
      ggplot(aes(time / 30, value, colour = var)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 3, lty = 2, colour = "grey50") +
      labs(x = "Months", y = "Score",
           title = "Eckardt component scores & total (Eckardt ≤3 = success)") +
      theme_minimal(base_size = 12)
  })

  output$p_tbe_qol <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(TBE5, QOL), names_to = "var", values_to = "value") %>%
      ggplot(aes(time / 30, value, colour = var)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Months", y = "") +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
  })

  output$p_multi_eck <- renderPlot({
    a <- active(); df <- a$multi
    ggplot(df, aes(time / 30, ECKARDT, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 3, lty = 2) +
      labs(x = "Months", y = "Eckardt score",
           title = "Eckardt score by treatment arm") +
      theme_minimal(base_size = 12)
  })

  output$p_multi_les <- renderPlot({
    a <- active(); df <- a$multi
    ggplot(df, aes(time / 30, LES_PRESS, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(15, 35), lty = 3) +
      labs(x = "Months", y = "LES pressure (mmHg)") +
      theme_minimal(base_size = 12)
  })

  output$p_proc <- renderPlot({
    proc_df <- data.frame(
      treatment = factor(c("PD", "LHM", "POEM"), levels = c("PD", "LHM", "POEM")),
      remission_2y = c(86, 90, 83),                       # Boeckxstaens 2011, Werner 2019
      remission_5y = c(60, 88, 79),                       # Ponds 2019 / Werner 2019
      gerd_rate    = c(25, 15, 45)
    )
    proc_long <- pivot_longer(proc_df, -treatment,
                              names_to = "endpoint", values_to = "pct")
    ggplot(proc_long, aes(treatment, pct, fill = endpoint)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = paste0(pct, "%")),
                position = position_dodge(width = 0.8), vjust = -0.3) +
      labs(y = "Patients (%)", x = NULL,
           title = "Procedural outcomes — published cohort anchors") +
      theme_minimal(base_size = 12)
  })

  output$tbl_proc <- renderDT({
    data.frame(
      Procedure  = c("Pneumatic dilation (35 mm)", "Heller myotomy + Dor",
                     "POEM"),
      `2y success` = c("86%", "90%", "83%"),
      `5y success` = c("60%", "88%", "79%"),
      `Post-Tx GERD` = c("25%", "15%", "45%"),
      Anchor = c("Boeckxstaens 2011 NEJM", "Werner 2019 NEJM",
                 "Ponds 2019 JAMA / Werner 2019")
    )
  }, options = list(dom = "t"))

  output$p_ae <- renderPlot({
    a <- active(); df <- a$main
    df %>%
      pivot_longer(c(AE_HA, AE_HYPO), names_to = "AE", values_to = "score") %>%
      ggplot(aes(time / 24, score, colour = AE)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "Severity (0-3)",
           title = "Drug-related adverse events") +
      theme_minimal(base_size = 12)
  })

  output$tbl_calib <- renderDT({
    data.frame(
      Anchor = c("Eckardt score validation (Eckardt 1992 GE)",
                 "Chicago Classification v4.0 (Yadlapati 2021 NGM)",
                 "Boeckxstaens 2011 NEJM (PD vs LHM)",
                 "Werner 2019 NEJM (POEM vs LHM)",
                 "Ponds 2019 JAMA (POEM vs PD)",
                 "Pasricha 1995/1996 (Botox pilot/RCT)",
                 "Triadafilopoulos 1991 DDS (Nifedipine SL)",
                 "Bortolotti 2000 GE (Sildenafil 50 mg)",
                 "Pandolfino 2008 GE (subtype response)",
                 "ACG Vaezi 2020 / ISDE Zaninotto 2018"),
      Use = c("0-12 symptom score, success ≤3",
              "IRP cutoff 15 mmHg; subtype I/II/III",
              "PD 86% vs LHM 90% 2y",
              "POEM 83% vs LHM 81% 2y (non-inferior)",
              "POEM 92% vs PD 54% 2y",
              "Botox 70% response 1 mo, ~32% at 6 mo",
              "Nifedipine 10–20 mg SL ↓LES 30–40%",
              "Sildenafil 50 mg ↓LES 35% × 2 h",
              "Type II best response across all Tx",
              "Guideline algorithm and Tx ordering")
    )
  }, options = list(dom = "t", pageLength = 12))
}

shinyApp(ui, server)
