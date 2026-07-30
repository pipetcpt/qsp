## =============================================================================
##  zes_shiny_app.R
##  Zollinger-Ellison Syndrome / gastrinoma — interactive QSP dashboard
##
##  Eleven tabs.  The design principle is that the app should make the model's
##  ONE claim checkable by hand rather than merely displayed:
##
##      BAO = BAOCAP x PCM x PUMPA x ACTP
##            (capacity) (factor 1) (factor 2) (factor 3)
##
##  So the "Three factors" tab plots the three factors on one panel with the
##  product beside them, and the "Headroom" tab shows the three counterfactual
##  integrators (normal parietal mass / normal drive / no pump inhibition)
##  running alongside the real trajectory, so the user can see WHICH factor is
##  responsible for the residual acid rather than being told.  The "Matched
##  contrasts" tab runs pairs of regimens that differ in exactly one thing --
##  the same 60 mg a day split three ways, the same dose with and without a
##  gastrinoma, the same drug in four CYP2C19 phenotypes.
##
##  Run:   shiny::runApp("zes_shiny_app.R")
##  Needs: shiny, mrgsolve, ggplot2, dplyr, tidyr  (DT optional)
##
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY --
##  not validated for clinical or regulatory use.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("zes_mrgsolve_model.R")   # builds `mod`, the phenotypes and the helpers

HAS_DT <- requireNamespace("DT", quietly = TRUE)

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "#eef2f7", colour = NA),
                legend.position  = "bottom"))

PAL <- c("#2b6cb0", "#c05621", "#2f855a", "#805ad5", "#b83280",
         "#4a5568", "#d69e2e", "#319795")

`%||%` <- function(a, b) if (is.null(a)) b else a
day <- function(x) x / 24

## night shading for the 24-hour panels
night_band <- function(days) {
  n <- ceiling(days)
  data.frame(xmin = 22 / 24 + seq_len(n) - 1 - 1,
             xmax = 30 / 24 + seq_len(n) - 1 - 1)
}

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Zollinger-Ellison syndrome / gastrinoma — QSP model explorer"),
  tags$p(style = "color:#666;margin-top:-8px",
         "59 ODEs. Acid output is a PRODUCT of parietal-cell mass x active-pump ",
         "fraction x secretagogue drive. Hypergastrinaemia raises the first and ",
         "third; every acid-suppressing drug acts only on the second."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      radioButtons("pheno", "Phenotype", selected = "ZES_SPOR",
        choices = c("Healthy stomach" = "HEALTHY",
                    "Sporadic ZES" = "ZES_SPOR",
                    "MEN1-ZES (with pHPT)" = "ZES_MEN1",
                    "Metastatic gastrinoma" = "ZES_MET")),
      sliderInput("tum",  "Primary gastrinoma (cm3)", 0, 20, 3, step = 0.5),
      sliderInput("met",  "Hepatic metastases (cm3)", 0, 60, 0, step = 1),
      sliderInput("grade", "Grade multiplier on growth (Ki-67 proxy)",
                  0.3, 3.0, 1.0, step = 0.1),
      sliderInput("sstr", "SSTR2 density (fold of reference)",
                  0.1, 2.0, 1.0, step = 0.05),
      sliderInput("gland", "Parathyroid functional mass (MEN1 pHPT)",
                  1.0, 3.0, 1.0, step = 0.1),
      selectInput("cyp", "CYP2C19 phenotype",
                  c("Ultrarapid (1.9)" = 1.9, "Normal (1.0)" = 1.0,
                    "Intermediate (0.5)" = 0.5, "Poor (0.28)" = 0.28),
                  selected = 1.0),
      sliderInput("renf", "Gastrin clearance factor (renal function)",
                  0.3, 1.2, 1.0, step = 0.05),
      checkboxInput("nsaid", "NSAID / aspirin exposure", FALSE),
      checkboxInput("hp", "H. pylori co-infection", FALSE),
      checkboxInput("meals", "Meals (08:00 / 13:00 / 19:00)", TRUE),

      hr(),
      h4("Acid suppression"),
      radioButtons("acid", "Class", selected = "ppi", inline = FALSE,
        choices = c("None" = "none",
                    "PPI (omeprazole-equivalent)" = "ppi",
                    "P-CAB (vonoprazan)" = "von",
                    "H2 antagonist (famotidine)" = "h2")),
      sliderInput("dose", "Dose per administration (mg)", 5, 160, 30, step = 5),
      selectInput("ii", "Dosing interval (h)",
                  c("24 (once daily)" = 24, "12 (twice daily)" = 12,
                    "8 (three times daily)" = 8, "6 (four times daily)" = 6),
                  selected = 12),

      hr(),
      h4("Tumour-directed therapy"),
      checkboxGroupInput("anti", NULL,
        choices = c("Octreotide LAR 30 mg q28d" = "oct",
                    "Everolimus 10 mg od" = "eve",
                    "Sunitinib 37.5 mg od" = "sun",
                    "CAPTEM q28d" = "captem",
                    "177Lu-DOTATATE 7.4 GBq x4" = "prrt",
                    "Netazepide 25 mg od" = "netaz")),
      checkboxInput("aalys", "PRRT amino-acid renal protection", TRUE),

      hr(),
      sliderInput("days", "Simulation length (days)", 3, 1825, 14, step = 1),
      sliderInput("runin", "Disease run-in before t=0 (weeks)",
                  4, 60, 40, step = 4),
      actionButton("go", "Simulate", class = "btn-primary"),
      tags$p(style = "color:#888;font-size:11px;margin-top:8px",
             "The run-in establishes the chronic phenotype: the trophic loops ",
             "settle for the tumour burden you set, with tumour growth frozen ",
             "so the burden you asked for is the burden at t = 0.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient profile",     plotOutput("p_profile", height = 620),
                 tableOutput("t_profile")),
        tabPanel("2 · Three factors",       plotOutput("p_factors", height = 660)),
        tabPanel("3 · Intragastric pH",     plotOutput("p_ph", height = 620),
                 tableOutput("t_ph")),
        tabPanel("4 · Drug PK",             plotOutput("p_pk", height = 620)),
        tabPanel("5 · Pump lifecycle",      plotOutput("p_pump", height = 620)),
        tabPanel("6 · Gastrin axis",        plotOutput("p_gastrin", height = 620)),
        tabPanel("7 · Clinical endpoints",  plotOutput("p_clin", height = 660)),
        tabPanel("8 · Headroom",            plotOutput("p_head", height = 520),
                 tableOutput("t_head")),
        tabPanel("9 · Dose titration",      plotOutput("p_titr", height = 460),
                 if (HAS_DT) DT::dataTableOutput("t_titr") else tableOutput("t_titr")),
        tabPanel("10 · Matched contrasts",  plotOutput("p_match", height = 700),
                 tableOutput("t_match")),
        tabPanel("11 · Tumour & harms",     plotOutput("p_tum", height = 660),
                 tableOutput("t_tum")),
        tabPanel("12 · Diagnostics",        tableOutput("t_secr"),
                 plotOutput("p_secr", height = 420))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  base_pheno <- reactive({
    p <- switch(input$pheno, HEALTHY = HEALTHY, ZES_SPOR = ZES_SPOR,
                ZES_MEN1 = ZES_MEN1, ZES_MET = ZES_MET)
    .pp(p, list(TUM0 = input$tum, MET0 = input$met, GRADEF = input$grade,
                SSTR2D = input$sstr, GLANDM = input$gland,
                MENFLG = if (input$pheno == "ZES_MEN1") 1 else 0))
  })

  extras <- reactive({
    list(CYPF = as.numeric(input$cyp), RENF = input$renf,
         NSAID = as.numeric(input$nsaid), HPYL = as.numeric(input$hp),
         MEALS = as.numeric(input$meals), AALYS = as.numeric(input$aalys))
  })

  events <- reactive({
    d <- input$days; ii <- as.numeric(input$ii); ev_list <- list()
    if (input$acid == "ppi") ev_list[[length(ev_list) + 1]] <-
      ev_ppi(input$dose, ii, d)
    if (input$acid == "von") ev_list[[length(ev_list) + 1]] <-
      ev_von(input$dose, ii, d)
    if (input$acid == "h2")  ev_list[[length(ev_list) + 1]] <-
      ev_h2(input$dose, ii, d)
    a <- input$anti %||% character(0)
    if ("oct" %in% a)    ev_list[[length(ev_list) + 1]] <-
      ev_oct(30, max(1, ceiling(d / 28)))
    if ("eve" %in% a)    ev_list[[length(ev_list) + 1]] <- ev_eve(10, d)
    if ("sun" %in% a)    ev_list[[length(ev_list) + 1]] <- ev_sun(37.5, d)
    if ("captem" %in% a) ev_list[[length(ev_list) + 1]] <-
      ev_captem(max(1, floor(d / 28)))
    if ("prrt" %in% a)   ev_list[[length(ev_list) + 1]] <- ev_prrt(7.4, 4)
    if ("netaz" %in% a)  ev_list[[length(ev_list) + 1]] <- ev_netaz(25, d)
    if (!length(ev_list)) NULL else do.call(c, ev_list)
  })

  init_state <- eventReactive(input$go, {
    zes_init(base_pheno(), weeks = input$runin, extra = extras())
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$go, {
    d  <- input$days
    dl <- if (d <= 20) 0.25 else if (d <= 120) 1 else 6
    zes_sim(base_pheno(), events(), d, init_state(), dl, extras())
  }, ignoreNULL = FALSE)

  ## an untreated reference run of the same patient, for every contrast panel
  sim_ref <- eventReactive(input$go, {
    d  <- input$days
    dl <- if (d <= 20) 0.25 else if (d <= 120) 1 else 6
    zes_sim(base_pheno(), NULL, d, init_state(), dl, extras())
  }, ignoreNULL = FALSE)

  ## ---- 1. patient profile -------------------------------------------------
  output$p_profile <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `Acid output (mEq/h)` = BAO, `Intragastric pH` = PH,
      `Fasting gastrin (pg/mL)` = FSG,
      `Parietal cell mass (fold)` = PCM, `ECL cell mass (fold)` = ECL,
      `Tumour volume (cm3)` = TUMTOT) %>%
      pivot_longer(-t)
    ggplot(d, aes(t, value)) +
      geom_line(colour = PAL[1], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "days", y = NULL,
           title = "Patient profile — the six numbers a clinic would follow")
  })

  output$t_profile <- renderTable({
    x <- sim(); r <- sim_ref()
    data.frame(
      Quantity = c("Basal (fasting) acid output, mEq/h",
                   "24-h mean acid output, mEq/h",
                   "Peak acid output in last day, mEq/h",
                   "Time with pH > 4, % of last 24 h",
                   "Nocturnal hours with pH < 4",
                   "Fasting serum gastrin, pg/mL",
                   "Pump pool inhibited, fraction",
                   "Duodenal ulcer index (0-1)",
                   "Diarrhoea index (0-1)"),
      Treated = c(round(bao_fasting(x), 1), round(bao_ss(x), 1),
                  round(bao_trough(x), 1), round(ph_hold(x), 0),
                  round(nab_hours(x), 1), round(fsg_fasting(x), 0),
                  round(mean(tail(x$PINH, 96)), 3),
                  round(tail(x$ULCD, 1), 3), round(tail(x$DIARR, 1), 3)),
      `Same patient, untreated` =
        c(round(bao_fasting(r), 1), round(bao_ss(r), 1),
          round(bao_trough(r), 1), round(ph_hold(r), 0),
          round(nab_hours(r), 1), round(fsg_fasting(r), 0),
          round(mean(tail(r$PINH, 96)), 3),
          round(tail(r$ULCD, 1), 3), round(tail(r$DIARR, 1), 3)),
      check.names = FALSE)
  }, digits = 3)

  ## ---- 2. the three factors ----------------------------------------------
  output$p_factors <- renderPlot({
    x <- sim()
    fac <- x %>% transmute(t = day(time),
      `FACTOR 1 · parietal cell mass (fold)` = PCM,
      `FACTOR 2 · active pump pool (fraction)` = PUMPA,
      `FACTOR 3 · per-pump activation ACTP` = ACTP,
      `PRODUCT · acid output (mEq/h)` = BAO) %>%
      pivot_longer(-t)
    fac$name <- factor(fac$name, levels = c(
      "FACTOR 1 · parietal cell mass (fold)",
      "FACTOR 2 · active pump pool (fraction)",
      "FACTOR 3 · per-pump activation ACTP",
      "PRODUCT · acid output (mEq/h)"))
    ggplot(fac, aes(t, value)) +
      geom_line(aes(colour = name), linewidth = 0.75, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(3, 1, 2, 6)]) +
      labs(x = "days", y = NULL,
           title = "The product and its three factors",
           subtitle = paste("Only factor 2 is reachable by any licensed",
                            "acid-suppressing drug. Factors 1 and 3 move only",
                            "if the gastrin signal moves."))
  })

  ## ---- 3. intragastric pH -------------------------------------------------
  output$p_ph <- renderPlot({
    x <- sim(); r <- sim_ref()
    d <- bind_rows(mutate(x, arm = "treated"), mutate(r, arm = "untreated")) %>%
      transmute(t = day(time), arm, PH, BAO) %>% pivot_longer(c(PH, BAO))
    d$name <- recode(d$name, PH = "Intragastric pH",
                     BAO = "Acid output (mEq/h)")
    tail_days <- max(1, min(3, input$days))
    d <- d[d$t > max(d$t) - tail_days, ]
    hl <- data.frame(name = c("Intragastric pH", "Acid output (mEq/h)"),
                     y = c(4, 10))
    ggplot(d, aes(t, value, colour = arm)) +
      geom_line(linewidth = 0.7) +
      geom_hline(data = hl, aes(yintercept = y), linetype = 2,
                 colour = "#888") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "days", y = NULL, colour = NULL,
           title = paste("Last", tail_days, "day(s): pH and acid output"),
           subtitle = paste("Dashed lines: the pH 4 holding-time threshold and",
                            "the 10 mEq/h treatment target."))
  })

  output$t_ph <- renderTable({
    x <- sim()
    data.frame(
      `pH > 4 (% of last 24 h)` = round(ph_hold(x), 0),
      `pH > 3 (%)` = round(ph_hold(x, thr = 3), 0),
      `Mean pH (last 24 h)` = round(mean(tail(x$PH, 96)), 2),
      `Nocturnal hours pH < 4` = round(nab_hours(x), 1),
      `Luminal [H+] (mmol/L, fasting)` = round(mean(tail(x$HCONC, 24)), 0),
      check.names = FALSE)
  })

  ## ---- 4. drug PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    x <- sim()
    cols <- c(CPPI = "PPI plasma (mg/L)", CVON = "Vonoprazan (mg/L)",
              CH2 = "Famotidine (mg/L)", COCT = "Octreotide (ng/mL)",
              CEVE = "Everolimus (mg/L)", CSUN = "Sunitinib (mg/L)",
              PPICAN = "Activated canalicular species (a.u.)",
              PRRTT = "Tumour-bound 177Lu (GBq)")
    keep <- names(cols)[vapply(names(cols),
                               function(k) max(x[[k]], na.rm = TRUE) > 1e-9,
                               logical(1))]
    if (!length(keep)) return(ggplot() + labs(title = "No drug on board"))
    d <- x[, c("time", keep)] %>% mutate(t = day(time)) %>% select(-time) %>%
      pivot_longer(-t)
    d$name <- cols[d$name]
    ggplot(d, aes(t, value)) + geom_line(colour = PAL[4], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "days", y = NULL, title = "Drug exposures")
  })

  ## ---- 5. pump lifecycle --------------------------------------------------
  output$p_pump <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `inactive (tubulovesicular)` = PUMPI, `ACTIVE (canalicular)` = PUMPA,
      `covalently blocked` = PUMPB, `reversibly blocked` = PUMPR) %>%
      pivot_longer(-t)
    tail_days <- max(2, min(10, input$days))
    d <- d[d$t > max(d$t) - tail_days, ]
    ggplot(d, aes(t, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.9) +
      scale_fill_manual(values = PAL[c(6, 3, 2, 5)]) +
      labs(x = "days", y = "fraction of the healthy total pump pool",
           fill = NULL, title = "H+/K+-ATPase pool, by state",
           subtitle = paste("The pool total is itself a variable: gastrin",
                            "up-regulates synthesis, so a hypergastrinaemic",
                            "stomach refills faster and reaches its trough",
                            "sooner."))
  })

  ## ---- 6. gastrin axis ----------------------------------------------------
  output$p_gastrin <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `G17 (pg/mL)` = G17, `G34 (pg/mL)` = G34,
      `Bioactive gastrin GBIO` = GBIO,
      `Antral G-cell output (fold)` = GANT,
      `Acid brake on the G cell (BRK)` = BRK,
      `Tumour secretory phenotype SPHEN` = SPHEN) %>% pivot_longer(-t)
    ggplot(d, aes(t, value)) + geom_line(colour = PAL[5], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "days", y = NULL,
           title = "Gastrin axis",
           subtitle = paste("BRK gates the ANTRAL G cell only. Nothing gates",
                            "the tumour -- that single missing edge is the",
                            "disease, and it is also why a PPI raises gastrin",
                            "in someone with no gastrinoma at all."))
  })

  ## ---- 7. clinical endpoints ---------------------------------------------
  output$p_clin <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `Duodenal ulcer index` = ULCD, `Oesophageal injury index` = ESOPH,
      `Diarrhoea index` = DIARR, `Mucosal defence` = MUCUS,
      `Cumulative bleed/perforation risk` = BLDRSK,
      `Excess acid load (mEq/h over defence)` = ALOAD) %>% pivot_longer(-t)
    ggplot(d, aes(t, value)) + geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "days", y = NULL, title = "Clinical endpoints")
  })

  ## ---- 8. headroom --------------------------------------------------------
  output$p_head <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `attributable to the enlarged parietal mass` = HFPCM,
      `attributable to the raised secretagogue drive` = HFDRV,
      `removed by pump inhibition` = HFPMP) %>% pivot_longer(-t)
    ggplot(d, aes(t, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL[c(3, 5, 1)]) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "days", y = "fraction of cumulative acid output", colour = NULL,
           title = "Measured headroom, not assumed headroom",
           subtitle = paste("Each curve compares the real cumulative acid",
                            "output with a counterfactual integrator running",
                            "alongside it with exactly one factor set to",
                            "normal."))
  })

  output$t_head <- renderTable({
    x <- sim(); e <- x[nrow(x), ]
    data.frame(
      Factor = c("1 · parietal cell mass", "3 · secretagogue drive",
                 "2 · pump inhibition"),
      `Counterfactual` = c("mass forced to normal", "drive forced to normal",
                           "blocked pumps returned to service"),
      `Share of acid output` = paste0(round(100 * c(e$HFPCM, e$HFDRV, e$HFPMP)),
                                      "%"),
      Reachable = c("only by lowering gastrin", "only by lowering gastrin",
                    "by PPI / P-CAB / H2RA"),
      check.names = FALSE)
  })

  ## ---- 9. dose titration --------------------------------------------------
  titr <- eventReactive(input$go, {
    titrate_ppi(base_pheno(), init_state(), interval = as.numeric(input$ii),
                extra = extras())
  }, ignoreNULL = FALSE)

  output$p_titr <- renderPlot({
    d <- titr()
    ggplot(d, aes(total_daily_mg)) +
      geom_line(aes(y = BAO_fasting, colour = "fasting (target)"),
                linewidth = 0.8) +
      geom_point(aes(y = BAO_fasting, colour = "fasting (target)")) +
      geom_line(aes(y = BAO_24h_mean, colour = "24-h mean"), linewidth = 0.8) +
      geom_line(aes(y = BAO_peak, colour = "worst hour"), linewidth = 0.8) +
      geom_hline(yintercept = 10, linetype = 2, colour = "#888") +
      scale_colour_manual(values = PAL[c(1, 3, 2)]) +
      labs(x = "total daily dose (omeprazole-equivalent mg)",
           y = "acid output (mEq/h)", colour = NULL,
           title = paste0("Titration at a ", input$ii,
                          "-hour interval — the maintenance dose as an OUTPUT"),
           subtitle = paste("Dashed line: the < 10 mEq/h treatment target.",
                            "Change the interval in the sidebar and the whole",
                            "curve moves, at unchanged total milligrams."))
  })

  output$t_titr <- if (HAS_DT) DT::renderDataTable({
    DT::datatable(titr(), rownames = FALSE, options = list(dom = "t"))
  }) else renderTable({ titr() })

  ## ---- 10. matched contrasts ---------------------------------------------
  match_runs <- eventReactive(input$go, {
    st <- init_state(); ph <- base_pheno(); ex <- extras(); d <- 10
    runs <- list(
      "60 mg once daily"        = ev_ppi(60, 24, d),
      "30 mg twice daily"       = ev_ppi(30, 12, d),
      "20 mg three times daily" = ev_ppi(20,  8, d),
      "15 mg four times daily"  = ev_ppi(15,  6, d))
    out <- lapply(names(runs), function(nm) {
      x <- zes_sim(ph, runs[[nm]], d, st, 0.25, ex); x$arm <- nm; x
    })
    do.call(rbind, out)
  }, ignoreNULL = FALSE)

  output$p_match <- renderPlot({
    a <- match_runs()
    a <- a[a$time > max(a$time) - 48, ]
    d <- a %>% transmute(t = day(time), arm,
                         `Acid output (mEq/h)` = BAO,
                         `Intragastric pH` = PH,
                         `Pump pool inhibited` = PINH) %>%
      pivot_longer(c(-t, -arm))
    ggplot(d, aes(t, value, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[1:4]) +
      labs(x = "days", y = NULL, colour = NULL,
           title = "60 mg a day, split four ways — identical total dose",
           subtitle = paste("If the daily dose were what mattered these four",
                            "curves would coincide. What separates them is",
                            "the rate at which the pump pool refills."))
  })

  output$t_match <- renderTable({
    a <- match_runs()
    do.call(rbind, lapply(split(a, a$arm), function(x) data.frame(
      Arm = x$arm[1], `Total mg/day` = 60,
      `Fasting BAO` = round(bao_fasting(x), 1),
      `24-h mean BAO` = round(bao_ss(x), 1),
      `Worst hour` = round(bao_trough(x), 1),
      `pH > 4 (%)` = round(ph_hold(x), 0),
      `Nocturnal h pH < 4` = round(nab_hours(x), 1),
      `At target` = bao_fasting(x) < 10,
      check.names = FALSE, row.names = NULL)))
  })

  ## ---- 11. tumour and long-term harms ------------------------------------
  output$p_tum <- renderPlot({
    x <- sim()
    d <- x %>% transmute(t = day(time),
      `Tumour volume (cm3)` = TUMTOT, `RECIST change (%)` = RECST,
      `Resistant clone fraction` = RESCL,
      `Cobalamin store (fold)` = B12, `Serum magnesium (mmol/L)` = MGS,
      `Bone mineral density (fold)` = BMD,
      `Gastric NET burden (cm3)` = GNET,
      `Cumulative renal dose (Gy)` = DOSEK) %>% pivot_longer(-t)
    ggplot(d, aes(t, value)) + geom_line(colour = PAL[7], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "days", y = NULL,
           title = "Tumour trajectory and the price of decades of achlorhydria",
           subtitle = paste("Set the simulation length to 1825 days to see the",
                            "harms develop; the 23 Gy renal constraint applies",
                            "to the PRRT panel."))
  })

  output$t_tum <- renderTable({
    x <- sim(); e <- x[nrow(x), ]
    data.frame(
      `Tumour volume (cm3)` = round(e$TUMTOT, 2),
      `RECIST change (%)` = round(e$RECST, 1),
      `Best RECIST (%)` = round(min(x$RECST), 1),
      `Tumour dose (Gy)` = round(max(x$DOSET), 1),
      `Renal dose (Gy)` = round(max(x$DOSEK), 1),
      `B12 (fold)` = round(e$B12, 2), `Mg (mmol/L)` = round(e$MGS, 2),
      `BMD (fold)` = round(e$BMD, 3), `Gastric NET (cm3)` = round(e$GNET, 2),
      check.names = FALSE)
  })

  ## ---- 12. diagnostics: the secretin sign inversion ----------------------
  secr <- eventReactive(input$go, {
    st <- init_state(); ph <- base_pheno()
    o <- zes_sim(ph, ev_ppi(40, 24, 28), 28, st, 1, extras())
    ppi_state <- .state_of(o, zero_integrators = FALSE)
    rbind(cbind(Subject = "this patient, off PPI",
                secretin_test(ph, st)),
          cbind(Subject = "this patient, after 4 weeks of PPI",
                secretin_test(ph, ppi_state)),
          cbind(Subject = "healthy reference",
                secretin_test(HEALTHY, zes_init(HEALTHY, weeks = 12))))
  }, ignoreNULL = FALSE)

  output$t_secr <- renderTable({
    d <- secr()
    names(d) <- c("Subject", "Basal gastrin (pg/mL)", "Peak (pg/mL)",
                  "Delta (pg/mL)", "Positive (Delta > 120)")
    d
  })

  output$p_secr <- renderPlot({
    st <- init_state(); ph <- base_pheno()
    a <- zes_sim(ph, ev_secretin(140), 0.5, st, 1 / 60,
                 .pp(extras(), list(MEALS = 0)))
    b <- zes_sim(HEALTHY, ev_secretin(140), 0.5,
                 zes_init(HEALTHY, weeks = 12), 1 / 60, list(MEALS = 0))
    d <- rbind(data.frame(t = a$time * 60, FSG = a$FSG, who = "this patient"),
               data.frame(t = b$time * 60, FSG = b$FSG, who = "healthy"))
    ggplot(d, aes(t, FSG, colour = who)) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL[c(2, 1)]) +
      labs(x = "minutes after 2 U/kg secretin IV", y = "serum gastrin (pg/mL)",
           colour = NULL, title = "The secretin sign inversion",
           subtitle = paste("Secretin INHIBITS the normal antral G cell and",
                            "STIMULATES the gastrinoma. Same hormone, opposite",
                            "sign, two compartments."))
  })
}

shinyApp(ui, server)
