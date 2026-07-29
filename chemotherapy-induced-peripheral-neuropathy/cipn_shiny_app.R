## =============================================================================
## cipn_shiny_app.R
## Interactive dashboard for the CIPN QSP model
## 항암화학요법 유발 말초신경병증 QSP 대시보드
##
## 12 tabs:
##   1  Patient & regimen      — build the exposure and the host
##   2  Drug PK                — plasma profiles, Cmax vs AUC vs Tc>threshold
##   3  Biophase & mechanism    — DRG platinum, adducts, transport, mitochondria
##   4  Axon & coasting         — the SARM1 program and post-treatment worsening
##   5  Clinical endpoints      — CIPN20, CTCAE grade, TNSc, BPI
##   6  Biomarkers              — plasma NfL, IENFD, and the NfL lead time
##   7  Acute vs chronic        — the two clocks of oxaliplatin
##   8  Scenario comparison     — all 18 protocols side by side
##   9  Population incidence    — exact grade>=2/3 incidence by arm
##  10  Route & schedule        — Cmax vs AUC (bortezomib SC/IV), cycle spacing
##  11  Therapeutic index       — the dose-intensity vs tumour-control optimum
##  12  Calibration & validation— what was fitted vs what was predicted
##
## Run with:  shiny::runApp("cipn_shiny_app.R")
## Requires:  cipn_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("cipn_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey92"))

GRADE_BANDS <- function() {
  pp <- param(mod)
  list(
    geom_hline(yintercept = c(pp$G1, pp$G2, pp$G3), linetype = 3,
               colour = c("grey50", "orange", "red"))
  )
}

## =============================================================================
## UI
## =============================================================================

ui <- fluidPage(
  titlePanel("CIPN QSP 대시보드 — Chemotherapy-Induced Peripheral Neuropathy"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. 항암제 요법 (Regimen)"),
      selectInput("drug", "Neurotoxic drug",
                  c("Oxaliplatin (FOLFOX/CAPOX)" = "OXA",
                    "Paclitaxel"                 = "PAC",
                    "Bortezomib"                 = "BTZ")),
      conditionalPanel(
        "input.drug == 'OXA'",
        sliderInput("oxa_dose", "Dose per cycle (mg/m2)", 40, 130, 85, 5),
        sliderInput("oxa_n", "Number of cycles", 1, 16, 12, 1),
        sliderInput("oxa_ii", "Cycle interval (days)", 7, 35, 14, 7)
      ),
      conditionalPanel(
        "input.drug == 'PAC'",
        sliderInput("pac_dose", "Dose per cycle (mg/m2)", 60, 200, 80, 5),
        sliderInput("pac_n", "Number of cycles", 1, 16, 12, 1),
        sliderInput("pac_ii", "Cycle interval (days)", 7, 28, 7, 7),
        sliderInput("pac_hr", "Infusion duration (h)", 1, 24, 1, 1)
      ),
      conditionalPanel(
        "input.drug == 'BTZ'",
        sliderInput("btz_dose", "Dose (mg/m2)", 0.7, 1.6, 1.3, 0.1),
        sliderInput("btz_cyc", "Cycles", 1, 12, 8, 1),
        radioButtons("btz_route", "Route", c("IV bolus" = "IV",
                                             "Subcutaneous" = "SC")),
        checkboxInput("btz_weekly", "Weekly schedule (d1/8/15/22)", FALSE)
      ),

      hr(),
      h4("2. 예방 전략 (Prevention)"),
      sliderInput("cryo", "Cryotherapy / compression: distal delivery blocked",
                  0, 0.8, 0, 0.05),
      helpText("Acts on the DRG/distal biophase only — systemic exposure, and",
               "therefore tumour kill, is unchanged."),

      hr(),
      h4("3. 증상 치료 (Symptomatic)"),
      checkboxInput("dul", "Duloxetine 60 mg/day", FALSE),
      checkboxInput("pgb", "Pregabalin 300 mg/day", FALSE),
      sliderInput("tx_start", "Start day of symptomatic treatment", 100, 350,
                  250, 10),
      sliderInput("tx_dur", "Duration (days)", 14, 180, 35, 7),

      hr(),
      h4("4. 환자 (Host)"),
      sliderInput("S", "Susceptibility multiplier S", 0.2, 4, 1, 0.05),
      sliderInput("axon0", "Baseline axon density (% of norm)", 55, 100, 100, 1),
      sliderInput("reserve", "Bioenergetic reserve", 0.5, 1.1, 1, 0.02),
      sliderInput("regen", "Regenerative capacity", 0.3, 1.2, 1, 0.05),
      sliderInput("oct2", "DRG platinum uptake (OCT2)", 0.4, 1.8, 1, 0.05),
      sliderInput("repair", "NER repair capacity (ERCC1)", 0.4, 1.8, 1, 0.05),
      selectInput("preset", "Host preset",
                  c("Reference" = "ref", "Diabetic" = "dm", "Age 75" = "old",
                    "Prior taxane exposure" = "prior",
                    "CMT1A carrier" = "cmt")),
      hr(),
      sliderInput("end", "Simulation horizon (days)", 200, 1500, 500, 50),
      helpText("Horizon must exceed treatment end + ~120 d to capture coasting.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "1 · 환자·요법",
          h4("Exposure summary"),
          tableOutput("summary_tbl"),
          h4("Where this patient sits"),
          plotOutput("overview_plot", height = "420px"),
          helpText("Dotted lines: CTCAE grade 1 / 2 / 3 thresholds on the",
                   "severity index. The index is fixed a priori at 10 / 25 / 45%",
                   "sensory-fibre-loss equivalent — it is not fitted.")
        ),

        tabPanel(
          "2 · 약물 PK",
          h4("Plasma concentration"),
          plotOutput("pk_plot", height = "320px"),
          h4("Cmax vs AUC vs time-above-threshold"),
          tableOutput("pk_tbl"),
          helpText("Bortezomib SC and IV have EQUAL AUC and an ~11-fold",
                   "difference in Cmax — the model's central natural experiment.")
        ),

        tabPanel(
          "3 · 바이오페이즈·기전",
          plotOutput("mech_plot", height = "560px"),
          helpText("DRG platinum accumulates over cycles and washes out with",
                   "t1/2 = 8 d; adducts suppress somal protein synthesis;",
                   "tubulin occupancy and proteasome inhibition degrade axonal",
                   "transport. ENERGY = MITO x transport x reserve is the",
                   "single quantity the axon-death program reads.")
        ),

        tabPanel(
          "4 · 축삭·Coasting",
          plotOutput("coast_plot", height = "480px"),
          h4("Coasting quantified"),
          tableOutput("coast_tbl"),
          helpText("The SARM1 program switches off with t1/2 = 23 d and blocks",
                   "regeneration while it runs, so axon density keeps falling",
                   "for weeks after the last dose. Nothing in the model forces",
                   "this — it is a consequence of those two rate constants.")
        ),

        tabPanel(
          "5 · 임상 엔드포인트",
          plotOutput("endpoint_plot", height = "480px"),
          tableOutput("endpoint_tbl")
        ),

        tabPanel(
          "6 · 바이오마커",
          plotOutput("biomarker_plot", height = "420px"),
          h4("NfL lead time over the clinical grade"),
          tableOutput("nfl_tbl"),
          helpText("Plasma NfL is released by degenerating axons and so rises",
                   "while the CTCAE grade is still 0-1. The lead time is the",
                   "window in which stopping oxaliplatin costs little tumour",
                   "kill — see the therapeutic-index tab.")
        ),

        tabPanel(
          "7 · 급성 vs 만성",
          plotOutput("acute_plot", height = "460px"),
          tableOutput("acute_tbl"),
          helpText("Two clocks, one drug: the acute channelopathy resolves in",
                   "days and grows with cycle number, while the chronic",
                   "axonopathy is at its worst once the acute syndrome is gone.")
        ),

        tabPanel(
          "8 · 시나리오 비교",
          h4("All 18 protocols"),
          DTOutput("scen_tbl"),
          plotOutput("scen_plot", height = "480px")
        ),

        tabPanel(
          "9 · 집단 발생률",
          h4("Exact grade >=2 / >=3 incidence by arm"),
          helpText("Peak severity is monotone in S, so P(grade >= g) is obtained",
                   "by bisecting for the boundary susceptibility and evaluating",
                   "the lognormal survival function — no Monte-Carlo noise."),
          actionButton("run_pop", "Compute incidence (≈30 sims per arm)"),
          DTOutput("pop_tbl"),
          plotOutput("pop_plot", height = "360px")
        ),

        tabPanel(
          "10 · 투여경로·스케줄",
          h4("Bortezomib: same AUC, 1/11 the Cmax"),
          plotOutput("route_plot", height = "400px"),
          h4("Cycle interval at matched cumulative dose"),
          plotOutput("sched_plot", height = "340px"),
          helpText("Convexity, not exposure, is what the subcutaneous route",
                   "buys: occupancy AUC is nearly identical but the injury flux",
                   "is a Hill function of occupancy, so the brief IV spike",
                   "dominates the damage integral.")
        ),

        tabPanel(
          "11 · 치료지수",
          h4("Optimum cumulative oxaliplatin dose"),
          sliderInput("w_cipn", "Weight on CIPN burden (utility units per pp)",
                      0.02, 0.5, 0.15, 0.01),
          actionButton("run_ti", "Compute therapeutic index"),
          plotOutput("ti_plot", height = "420px"),
          DTOutput("ti_tbl"),
          helpText("Neurotoxicity integrates linearly with cumulative dose;",
                   "the anti-tumour exposure-response saturates by ~500 mg/m2.",
                   "An optimum therefore exists, and it moves with RISK, not",
                   "with neurotoxicity — the CIPN cost curve is the same in",
                   "both strata.")
        ),

        tabPanel(
          "12 · 보정·검증",
          h4("Fitted parameters (6) vs trial observations (6)"),
          tableOutput("calib_tbl"),
          h4("Out-of-sample predictions"),
          tableOutput("valid_tbl"),
          helpText("Only the rows marked FIT were used to set parameters.",
                   "Every other row is a prediction of the calibrated model.")
        )
      )
    )
  )
)

## =============================================================================
## SERVER
## =============================================================================

server <- function(input, output, session) {

  observeEvent(input$preset, {
    v <- switch(input$preset,
      ref   = list(100, 1.00, 1.00),
      dm    = list(82,  0.90, 0.70),
      old   = list(96,  0.94, 0.70),
      prior = list(90,  0.97, 0.90),
      cmt   = list(70,  0.88, 0.55))
    updateSliderInput(session, "axon0",   value = v[[1]])
    updateSliderInput(session, "reserve", value = v[[2]])
    updateSliderInput(session, "regen",   value = v[[3]])
  })

  ## ---- events ---------------------------------------------------------------
  events <- reactive({
    e <- switch(input$drug,
      OXA = ev_oxa(input$oxa_n, dose = input$oxa_dose, ii = input$oxa_ii),
      PAC = ev_pac(input$pac_n, dose = input$pac_dose, ii = input$pac_ii,
                   hours = input$pac_hr),
      BTZ = ev_btz(input$btz_cyc, input$btz_route, dose = input$btz_dose,
                   weekly = input$btz_weekly))
    if (input$dul) e <- c(e, ev_duloxetine(input$tx_start,
                                           input$tx_start + input$tx_dur))
    if (input$pgb) e <- c(e, ev_pregabalin(input$tx_start,
                                           input$tx_start + input$tx_dur))
    e
  })

  pars <- reactive(list(
    S = input$S, AXON0 = input$axon0, RESERVE = input$reserve,
    REGEN = input$regen, OCT2 = input$oct2, REPAIR = input$repair,
    CRYO = input$cryo))

  sim <- reactive({
    mrgsim(param(mod, pars()), events = events(), end = input$end,
           delta = 0.25) %>% as_tibble()
  })

  last_dose <- reactive({
    switch(input$drug,
      OXA = (input$oxa_n - 1) * input$oxa_ii,
      PAC = (input$pac_n - 1) * input$pac_ii,
      BTZ = { d <- if (input$btz_weekly) 21 else 10
              cy <- if (input$btz_weekly) 35 else 21
              (input$btz_cyc - 1) * cy + d })
  })

  ## ---- 1 overview -----------------------------------------------------------
  output$summary_tbl <- renderTable({
    s <- sim()
    tibble(
      Quantity = c("Cumulative platinum (mg/m2)", "Cumulative taxane (mg/m2)",
                   "Last dose (day)", "Peak severity index",
                   "Day of peak severity", "Coasting (days after last dose)",
                   "Peak CTCAE grade", "Nadir IENFD (fibres/mm)",
                   "CIPN20 at day 365", "Predicted 3-year DFS (%)"),
      Value = c(sprintf("%.0f", max(s$CUMPT)), sprintf("%.0f", max(s$CUMTAX)),
                sprintf("%.0f", last_dose()), sprintf("%.3f", max(s$CS)),
                sprintf("%.0f", s$time[which.max(s$CS)]),
                sprintf("%+.0f", s$time[which.max(s$CS)] - last_dose()),
                sprintf("%.0f", max(s$GRADE)), sprintf("%.2f", min(s$IENFD)),
                sprintf("%.1f", s$CIPN20[which.min(abs(s$time - 365))]),
                sprintf("%.1f", 100 * last(s$DFS3))))
  })

  output$overview_plot <- renderPlot({
    s <- sim()
    s %>% select(time, CS, CIPN20, IENFD, AXON) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = "steelblue") +
      geom_vline(xintercept = last_dose(), linetype = 2, colour = "red") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "Red dashed line = last neurotoxic dose") + THEME
  })

  ## ---- 2 PK ----------------------------------------------------------------
  output$pk_plot <- renderPlot({
    s <- sim()
    v <- switch(input$drug, OXA = "C_OXA", PAC = "C_PAC", BTZ = "C_BTZ")
    tmax <- min(input$end, last_dose() + 30)
    s %>% filter(time <= tmax) %>%
      select(time, all_of(v)) %>% pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      labs(x = "Day", y = "Plasma concentration",
           title = paste("Plasma", v)) + THEME
  })

  output$pk_tbl <- renderTable({
    s <- sim()
    v <- switch(input$drug, OXA = s$C_OXA, PAC = s$C_PAC, BTZ = s$C_BTZ)
    dt <- 0.25
    tibble(Metric = c("Cmax", "AUC (conc x day)", "Peak DRG occupancy",
                      "Paclitaxel time above 0.05 umol/L (h)"),
           Value = c(sprintf("%.3f", max(v)), sprintf("%.3f", sum(v) * dt),
                     sprintf("%.3f", max(pmax(s$TOCCo, s$PI_DRGo))),
                     sprintf("%.1f", 24 * max(s$TCTHR))))
  })

  ## ---- 3 mechanism ---------------------------------------------------------
  output$mech_plot <- renderPlot({
    sim() %>%
      select(time, CE_PT, CE_TAX, CE_BTZ, ADDUCT, PSYN, MITO, ATRANS, ROS,
             ENERGYo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.7, colour = "darkgreen") +
      geom_vline(xintercept = last_dose(), linetype = 2, colour = "red") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL) + THEME
  })

  ## ---- 4 coasting ----------------------------------------------------------
  output$coast_plot <- renderPlot({
    sim() %>% select(time, SARM, AXON, NEURON, CS) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8) +
      geom_vline(xintercept = last_dose(), linetype = 2, colour = "red") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "Axon density keeps falling after the last dose") + THEME
  })

  output$coast_tbl <- renderTable({
    s <- sim(); ld <- last_dose()
    i_ld <- which.min(abs(s$time - ld))
    tibble(
      Quantity = c("SARM1 activity at last dose", "SARM1 at peak severity",
                   "AXON at last dose (%)", "AXON at nadir (%)",
                   "Further loss after the last dose (pp)",
                   "Day of nadir", "Coasting (days)"),
      Value = c(sprintf("%.3f", s$SARM[i_ld]),
                sprintf("%.3f", s$SARM[which.max(s$CS)]),
                sprintf("%.1f", s$AXON[i_ld]), sprintf("%.1f", min(s$AXON)),
                sprintf("%.1f", s$AXON[i_ld] - min(s$AXON)),
                sprintf("%.0f", s$time[which.min(s$AXON)]),
                sprintf("%+.0f", s$time[which.min(s$AXON)] - ld)))
  })

  ## ---- 5 endpoints ---------------------------------------------------------
  output$endpoint_plot <- renderPlot({
    s <- sim()
    p1 <- s %>% select(time, CIPN20, TNSc, BPI, GRADE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = "firebrick") +
      geom_vline(xintercept = last_dose(), linetype = 2) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
    p1
  })

  output$endpoint_tbl <- renderTable({
    s <- sim()
    tibble(Endpoint = c("Peak CIPN20", "Peak TNSc", "Peak BPI", "Peak grade",
                        "CIPN20 at 12 months", "Grade at 12 months",
                        "Grade at 24 months"),
           Value = c(sprintf("%.1f", max(s$CIPN20)), sprintf("%.1f", max(s$TNSc)),
                     sprintf("%.2f", max(s$BPI)), sprintf("%.0f", max(s$GRADE)),
                     sprintf("%.1f", s$CIPN20[which.min(abs(s$time - 365))]),
                     sprintf("%.0f", s$GRADE[which.min(abs(s$time - 365))]),
                     sprintf("%.0f", s$GRADE[which.min(abs(s$time - 730))])))
  })

  ## ---- 6 biomarkers --------------------------------------------------------
  output$biomarker_plot <- renderPlot({
    sim() %>% select(time, NFL, IENFD, MAC, IL1B) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = "purple") +
      geom_vline(xintercept = last_dose(), linetype = 2) +
      facet_wrap(~name, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
  })

  output$nfl_tbl <- renderTable({
    s <- sim(); b <- first(s$NFL)
    t_nfl <- suppressWarnings(min(s$time[s$NFL >= 2 * b]))
    t_g1  <- suppressWarnings(min(s$time[s$GRADE >= 1]))
    t_g2  <- suppressWarnings(min(s$time[s$GRADE >= 2]))
    f <- function(x) if (is.finite(x)) sprintf("%.0f", x) else "not reached"
    tibble(Quantity = c("Baseline NfL (pg/mL)", "Peak NfL (pg/mL)",
                        "Day NfL doubles", "Day grade 1 reached",
                        "Day grade 2 reached", "NfL lead time over grade 2 (d)"),
           Value = c(sprintf("%.1f", b), sprintf("%.1f", max(s$NFL)),
                     f(t_nfl), f(t_g1), f(t_g2),
                     if (is.finite(t_g2) && is.finite(t_nfl))
                       sprintf("%.0f", t_g2 - t_nfl) else "n/a"))
  })

  ## ---- 7 acute vs chronic --------------------------------------------------
  output$acute_plot <- renderPlot({
    s <- sim()
    s %>% select(time, COLDA, EXCITC, CENTS, AXON) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.7) +
      geom_vline(xintercept = last_dose(), linetype = 2, colour = "red") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "COLDA spikes and resolves each cycle; AXON only falls") +
      THEME
  })

  output$acute_tbl <- renderTable({
    s <- sim(); ld <- last_dose()
    after <- s %>% filter(time > ld + 40)
    tibble(Quantity = c("Peak acute cold allodynia",
                        "Acute state 40 d after last dose",
                        "Axon density 40 d after last dose (%)"),
           Value = c(sprintf("%.3f", max(s$COLDA)),
                     sprintf("%.4f", if (nrow(after)) first(after$COLDA) else NA),
                     sprintf("%.1f", if (nrow(after)) first(after$AXON) else NA)))
  })

  ## ---- 8 scenarios ---------------------------------------------------------
  scen <- reactive(summarise_scenarios())

  output$scen_tbl <- renderDT({
    datatable(scen() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 20, scrollX = TRUE))
  })

  output$scen_plot <- renderPlot({
    scen() %>%
      ggplot(aes(cum_Pt + cum_taxane, peak_CS, label = scenario)) +
      geom_point(size = 3, colour = "steelblue") +
      geom_hline(yintercept = c(param(mod)$G2, param(mod)$G3), linetype = 3,
                 colour = c("orange", "red")) +
      labs(x = "Cumulative neurotoxic dose (mg/m2)", y = "Peak severity index",
           title = "Peak severity vs cumulative dose across all 18 protocols") +
      THEME
  })

  ## ---- 9 population --------------------------------------------------------
  pop <- eventReactive(input$run_pop, {
    arms <- list(
      list("FOLFOX 6 months", ev_oxa(12), 500),
      list("FOLFOX 3 months", ev_oxa(6), 400),
      list("CAPOX 3 months",  ev_oxa(4, dose = 130, ii = 21), 400),
      list("CAPOX 6 months",  ev_oxa(8, dose = 130, ii = 21), 500),
      list("Paclitaxel weekly x12", ev_pac(12), 400),
      list("Paclitaxel q3w x4", ev_pac(4, dose = 175, ii = 21, hours = 3), 400),
      list("Bortezomib IV", ev_btz(8, "IV"), 500),
      list("Bortezomib SC", ev_btz(8, "SC"), 500))
    withProgress(message = "Computing exact incidence", value = 0, {
      bind_rows(lapply(seq_along(arms), function(i) {
        a <- arms[[i]]; incProgress(1 / length(arms), detail = a[[1]])
        tibble(arm = a[[1]],
               grade2 = 100 * incidence(a[[2]], a[[3]], 2),
               grade3 = 100 * incidence(a[[2]], a[[3]], 3))
      }))
    })
  })

  output$pop_tbl <- renderDT(datatable(
    pop() %>% mutate(across(where(is.numeric), ~round(.x, 1)))))

  output$pop_plot <- renderPlot({
    pop() %>% pivot_longer(-arm) %>%
      ggplot(aes(reorder(arm, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      labs(x = NULL, y = "Incidence (%)") + THEME
  })

  ## ---- 10 route & schedule -------------------------------------------------
  output$route_plot <- renderPlot({
    bind_rows(run_scenario(scenarios$S10), run_scenario(scenarios$S11)) %>%
      filter(time <= 14) %>%
      select(time, scenario, C_BTZ, PI_DRGo) %>%
      pivot_longer(c(C_BTZ, PI_DRGo)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "Day", y = NULL) + THEME
  })

  output$sched_plot <- renderPlot({
    bind_rows(
      run_scenario(scenarios$S01) %>% mutate(arm = "q2w x12 (6 months)"),
      run_scenario(scenarios$S05) %>% mutate(arm = "q3w x12 (9 months)"),
      run_scenario(scenarios$S06) %>% mutate(arm = "stop-and-go")) %>%
      ggplot(aes(time, CS, colour = arm)) + geom_line(linewidth = 0.8) +
      GRADE_BANDS() +
      labs(x = "Day", y = "Severity index",
           title = "Identical cumulative dose (1020 mg/m2), different spacing") +
      THEME
  })

  ## ---- 11 therapeutic index ------------------------------------------------
  ti <- eventReactive(input$run_ti, {
    withProgress(message = "Computing therapeutic index", value = 0.5,
                 therapeutic_index(w = input$w_cipn))
  })

  output$ti_plot <- renderPlot({
    ti() %>% select(cum_Pt, grade2_pct, utility_high, utility_low) %>%
      pivot_longer(-cum_Pt) %>%
      ggplot(aes(cum_Pt, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      labs(x = "Cumulative oxaliplatin (mg/m2)", y = "percentage points",
           title = "Neurotoxicity integrates; anti-tumour benefit saturates") +
      THEME
  })

  output$ti_tbl <- renderDT(datatable(
    ti() %>% mutate(across(where(is.numeric), ~round(.x, 2)))))

  ## ---- 12 calibration ------------------------------------------------------
  output$calib_tbl <- renderTable({
    pp <- param(mod)
    tibble(
      Parameter = c("PHI_RELIEF", "SIGMA_S", "KDAM_PT", "KDAM_TAX",
                    "KDAM_BTZ", "BTZ_JH"),
      Value = sprintf("%.4f", c(pp$PHI_RELIEF, pp$SIGMA_S, pp$KDAM_PT,
                                pp$KDAM_TAX, pp$KDAM_BTZ, pp$BTZ_JH)),
      `Fitted to` = c("IDEA FOLFOX 3 months grade>=2 = 16.6%",
                      "MOSAIC FOLFOX 6 months grade>=3 = 12.4%",
                      "IDEA FOLFOX 6 months grade>=2 = 47.7%",
                      "ECOG 1199 weekly paclitaxel grade>=2 = 27%",
                      "MMY-3021 bortezomib IV grade>=2 = 41%",
                      "MMY-3021 bortezomib SC grade>=2 = 24%"))
  })

  output$valid_tbl <- renderTable({
    tibble(
      Prediction = c("FOLFOX 3 months grade>=3",
                     "CAPOX 3 months / 6 months grade>=2",
                     "Paclitaxel q3w 175 x4 grade>=2",
                     "Bortezomib IV weekly schedule grade>=2",
                     "MOSAIC grade 3 at 12 / 48 months",
                     "IDEA 3-year DFS, low and high risk",
                     "Duloxetine change in BPI over 5 weeks",
                     "Coasting: peak severity after last dose"),
      Observed = c("2.7%", "~15% / ~45%", "20%", "~24%", "1.1% / 0.7%",
                   "83.1/83.3% and 62.7/64.4%", "-0.73 net", "4-12 weeks"),
      Source = c("IDEA 2018", "IDEA 2018", "ECOG 1199", "Bringhen 2010",
                 "MOSAIC 2009", "IDEA 2018", "Smith 2013 JAMA",
                 "Pachman 2015 / Briani 2014"))
  })
}

shinyApp(ui, server)
