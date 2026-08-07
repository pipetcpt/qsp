## =====================================================================
##  lep_shiny_app.R
##  Leprosy (Hansen's disease) QSP model — interactive dashboard
##  한센병 QSP 모델 — 인터랙티브 대시보드
##
##  10 tabs, built around the one thing the model exists to show:
##  the bacterial index and the morphological index are reading two
##  DIFFERENT clocks, and every therapeutic argument in leprosy turns on
##  which of the two a given decision is actually using.
##
##    Tab 1  Patient           spectrum, burden, immunological set-point
##    Tab 2  Drug PK           the four PK profiles on one time axis
##    Tab 3  Two clocks        BI vs MI vs viable burden — the core figure
##    Tab 4  Antigen           liberation rate, cumulative area, conservation
##    Tab 5  Reactions         ENL and reversal reaction, with the drivers
##    Tab 6  Nerve             reversible pool, permanent sink, the window
##    Tab 7  Scenarios         compare any subset of the 17 regimens
##    Tab 8  Relapse           residual viable burden and the escape index
##    Tab 9  Safety            methaemoglobin, haemoglobin, pigmentation, HPA
##    Tab 10 Model             equations, parameters, anchors, caveats
##
##  Run:  shiny::runApp("lep_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, ggplot2, tidyr, DT
##  DISCLAIMER: educational and research model.  Not for clinical use.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("lep_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

PAL <- c("#1f6fb2", "#b22222", "#2e8b57", "#d2691e",
         "#6a51a3", "#8c8c8c", "#c9a227", "#0f7b7b")

## ---------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Leprosy (Hansen's disease) — QSP model / 한센병 정량적 시스템 약리학 모델"),
  tags$p(style = "color:#555;",
         "The bacterial index reads the clearance clock (host-set, t½ ≈ 137 d).",
         "The morphological index reads the kill clock (drug-set, hours to weeks).",
         "Educational model — not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("pt", "Ridley-Jopling pole",
                  choices = c("Lepromatous (LL)"            = "LL",
                              "Borderline lepromatous (BL)"  = "BL",
                              "Borderline tuberculoid (BT)"  = "BT",
                              "Tuberculoid (TT)"             = "TT",
                              "Polar LL below containment"   = "POLAR"),
                  selected = "LL"),
      sliderInput("bi0", "Bacterial index at diagnosis", 0, 6.5, 5.9, 0.1),
      sliderInput("spec", "Immunological set-point (SPEC)",
                  0.005, 1.0, 0.05, 0.005),
      helpText(HTML("Containment threshold <b>SPEC* = (MUMAX − KNAT)/KHOST
                     = 0.0225</b>. Below it no residual clone can be held.")),

      h4("Regimen"),
      checkboxInput("rif", "Rifampicin 600 mg monthly", TRUE),
      checkboxInput("dds", "Dapsone 100 mg daily", TRUE),
      checkboxInput("clo", "Clofazimine 50 mg daily + 300 mg monthly", TRUE),
      checkboxInput("cload", "…with a 1-month 300 mg/d loading phase", FALSE),
      checkboxInput("rom", "Add ofloxacin + minocycline daily", FALSE),
      sliderInput("months", "Duration of MDT (months)", 0, 24, 12, 1),

      h4("Reaction management"),
      checkboxInput("pdn", "Prednisolone taper", FALSE),
      sliderInput("pdn_start", "Start (day)", 0, 400, 20, 10),
      sliderInput("pdn_weeks", "Length (weeks)", 8, 78, 20, 2),
      sliderInput("pdn_floor", "Maintenance floor (mg)", 5, 20, 5, 5),
      checkboxInput("tha", "Thalidomide 300 mg for 12 weeks", FALSE),

      h4("Host factors"),
      checkboxInput("g6pd", "G6PD deficiency", FALSE),
      sliderInput("fres", "Pre-existing resistant fraction (log10)",
                  -9, -3, -9, 0.5),
      selectInput("resto", "Resistance phenotype",
                  choices = c("rifampicin (rpoB)" = "rif",
                              "dapsone (folP1)"   = "dds",
                              "quinolone (gyrA)"  = "rom")),
      sliderInput("kmaxr_scale", "Bactericidal rate multiplier (experiment)",
                  0.02, 2, 1, 0.02),
      numericInput("end", "Follow-up (days)", 1825, 365, 5475, 365)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient",
                 plotOutput("p_patient", height = 460),
                 tableOutput("t_patient")),
        tabPanel("2 · Drug PK",
                 plotOutput("p_pk", height = 520),
                 helpText("Clofazimine is plotted as its tissue depot, the slow
                           compartment that carries both the antimycobacterial
                           and the anti-TNF effect.")),
        tabPanel("3 · Two clocks",
                 plotOutput("p_clocks", height = 520),
                 tableOutput("t_clocks"),
                 helpText(HTML("<b>The point of the model.</b> The upper panel is
                    what the programme measures; the lower panel is what the drug
                    is doing. A regimen change moves the lower panel by orders of
                    magnitude and the upper panel by less than the reading error
                    of a slit-skin smear."))),
        tabPanel("4 · Antigen",
                 plotOutput("p_ag", height = 520),
                 helpText(HTML("Antigen becomes available on DIGESTION, not on
                    killing, and live bacilli suppress digestion. Killing the
                    population therefore un-blocks the whole standing load at
                    once. The cumulative curve is set by the burden at diagnosis;
                    the regimen moves only its distribution in time."))),
        tabPanel("5 · Reactions",
                 plotOutput("p_react", height = 560),
                 tableOutput("t_react")),
        tabPanel("6 · Nerve",
                 plotOutput("p_nerve", height = 460),
                 plotOutput("p_window", height = 320),
                 helpText("The lower panel re-runs the current patient with the
                           steroid course started at a range of delays.")),
        tabPanel("7 · Scenarios",
                 checkboxGroupInput("scen", "Scenarios", inline = TRUE,
                                    choices = names(SCEN),
                                    selected = c("S01","S02","S04","S06","S09")),
                 plotOutput("p_scen", height = 520),
                 DT::dataTableOutput("t_scen")),
        tabPanel("8 · Relapse",
                 plotOutput("p_relapse", height = 460),
                 tableOutput("t_relapse"),
                 helpText(HTML("The continuum ODEs cannot represent stochastic
                    extinction of a few thousand organisms, so relapse is
                    reported as a Poisson escape index computed from the
                    residual VIABLE burden at the end of treatment, not read
                    off the deterministic trajectory."))),
        tabPanel("9 · Safety",
                 plotOutput("p_safety", height = 560)),
        tabPanel("10 · Model",
                 h4("Structure"),
                 verbatimTextOutput("txt_struct"),
                 h4("Parameters"),
                 DT::dataTableOutput("t_par"),
                 h4("Caveats"),
                 verbatimTextOutput("txt_caveat"))
      )
    )
  )
)

## ---------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$pt, {
    p <- PATIENT[[input$pt]]
    updateSliderInput(session, "bi0",
                      value = round(log10(p$B0 / (1 - 0.19) ) + 3, 1))
    updateSliderInput(session, "spec", value = p$SPEC)
  })

  pars <- reactive({
    b0 <- 10^(input$bi0 - 3) * 0.19        # BI is total; ~19% of it is live
    p <- list(B0 = max(b0, 1e-9), SPEC = input$spec,
              KMAXR = 12.0 * input$kmaxr_scale,
              FRES = 10^input$fres,
              RRIF = ifelse(input$resto == "rif", 0, 1),
              RDAP = ifelse(input$resto == "dds", 0, 1),
              RROM = ifelse(input$resto == "rom", 0, 1),
              RCLO = 1)
    if (input$g6pd) { p$GRED <- 0.35; p$GHEM <- 6.0 }
    p
  })

  events <- reactive({
    m <- input$months
    d <- NULL
    if (m > 0) {
      if (input$rif) d <- rbind(d, dose_seq("RIFG", 600, 0, 30, m))
      if (input$dds) d <- rbind(d, dose_seq("DAPG", 100, 0, 1, 30 * m))
      if (input$clo) d <- rbind(d, dose_seq("CLOG", 300, 0, 30, m),
                                   dose_seq("CLOG",  50, 0,  1, 30 * m))
      if (input$clo && input$cload) d <- rbind(d, dose_seq("CLOG", 250, 0, 1, 30))
      if (input$rom) d <- rbind(d, dose_seq("ROMG", 400, 0, 1, 30 * m))
    }
    if (input$pdn) d <- rbind(d, pred_taper(input$pdn_start,
                                            weeks = input$pdn_weeks,
                                            floor = input$pdn_floor))
    if (input$tha) d <- rbind(d, thalidomide(30, 84, 300))
    d
  })

  sim <- reactive({
    m <- lep_mod %>% param(pars())
    ee <- make_ev(events())
    out <- if (is.null(ee)) m %>% mrgsim(end = input$end, delta = 1)
           else              m %>% mrgsim(events = ee, end = input$end, delta = 1)
    as.data.frame(out)
  })

  long <- function(d, cols, labs) {
    d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time, names_to = "v", values_to = "y") %>%
      mutate(v = factor(v, levels = cols, labels = labs))
  }

  ## ---- 1. patient ----------------------------------------------------
  output$p_patient <- renderPlot({
    d <- sim()
    ggplot(long(d, c("BI", "MI", "LES", "CMI"),
                c("Bacterial index (log10/g)", "Morphological index (%)",
                  "Skin lesion activity (0-100)", "Cell-mediated index (-)")),
           aes(time, y)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "Disease course for the selected patient") + THEME
  })
  output$t_patient <- renderTable({
    d <- sim(); at <- function(t, c) d[[c]][which.min(abs(d$time - t))]
    data.frame(
      Quantity = c("BI at diagnosis", "BI at 1 year", "BI fall in year 1 (log10)",
                   "MI at diagnosis (%)", "Viable burden at diagnosis (log10/g)",
                   "Immunological set-point", "Containment threshold SPEC*"),
      Value = c(round(at(0, "BI"), 2), round(at(365, "BI"), 2),
                round(at(0, "BI") - at(365, "BI"), 2),
                round(at(0, "MI"), 1), round(at(0, "BILIV"), 2),
                input$spec, 0.0225))
  })

  ## ---- 2. PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time <= min(input$end, 400))
    ggplot(long(d, c("CRIF", "CDAP", "CLOE", "CPDN"),
                c("Rifampicin (mg/L)", "Dapsone (mg/L)",
                  "Clofazimine tissue depot (x standard)",
                  "Prednisolone (mg/L)")),
           aes(time, y, colour = v)) +
      geom_line(linewidth = 0.8, show.legend = FALSE) +
      facet_wrap(~v, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "Day", y = NULL,
           title = "Pharmacokinetics: three drugs on three different time scales",
           subtitle = "rifampicin t½ 4 h  ·  dapsone t½ 28 h  ·  clofazimine t½ 70 DAYS") +
      THEME
  })

  ## ---- 3. the two clocks ---------------------------------------------
  output$p_clocks <- renderPlot({
    d <- sim()
    a <- d %>% transmute(time, `Bacterial index (measured)` = BI,
                         `Viable burden (log10/g, unmeasured)` = BILIV) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(a, aes(time, y, colour = v)) +
      geom_line(linewidth = 1.0) +
      scale_colour_manual(values = c(PAL[2], PAL[1]), name = NULL) +
      coord_cartesian(ylim = c(-6, 7)) +
      labs(x = "Day", y = "log10 organisms per gram",
           title = "Clock 2 (what is measured) against clock 1 (what the drug does)",
           subtitle = "The gap between the two lines is the model's whole argument") +
      THEME
  })
  output$t_clocks <- renderTable({
    d <- sim(); at <- function(t, c) d[[c]][which.min(abs(d$time - t))]
    data.frame(Day = c(0, 21, 90, 365, 730),
               BI = round(sapply(c(0, 21, 90, 365, 730), at, "BI"), 2),
               MI_pct = signif(sapply(c(0, 21, 90, 365, 730), at, "MI"), 3),
               log10_viable = round(sapply(c(0, 21, 90, 365, 730), at, "BILIV"), 2),
               ENL = round(sapply(c(0, 21, 90, 365, 730), at, "ENL"), 1))
  })

  ## ---- 4. antigen ----------------------------------------------------
  output$p_ag <- renderPlot({
    d <- sim()
    ggplot(long(d, c("AGRATE", "AGC", "AG", "AB"),
                c("Antigen liberation RATE (AU/d)",
                  "Cumulative antigen (AU) — the conserved quantity",
                  "Free antigen pool (AU)", "Anti-PGL-1 (AU)")),
           aes(time, y)) +
      geom_line(colour = PAL[4], linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "Antigen: the rate is drug-modifiable, the area is not") + THEME
  })

  ## ---- 5. reactions --------------------------------------------------
  output$p_react <- renderPlot({
    d <- sim()
    ggplot(long(d, c("ENL", "T1R", "TNF", "IC", "NEU", "ENLC"),
                c("ENL score (0-100)", "Reversal reaction activity (0-100)",
                  "TNF-alpha (x baseline)", "Immune complexes (AU)",
                  "Neutrophil infiltrate (AU)", "Cumulative ENL (score-days)")),
           aes(time, y)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Day", y = NULL,
           title = "Type 2 reaction (ENL) and type 1 reaction (reversal)",
           subtitle = "ENL needs antigen RELEASE x antibody, so it belongs to the lepromatous pole; reversal needs restored T-cell function x residual antigen, so it belongs to the borderline zone") +
      THEME
  })
  output$t_react <- renderTable({
    d <- sim()
    data.frame(Quantity = c("Peak ENL score", "Day of peak ENL",
                            "Cumulative ENL (score-days)",
                            "Peak reversal-reaction activity",
                            "Day of peak reversal reaction"),
               Value = c(round(max(d$ENL), 1), d$time[which.max(d$ENL)],
                         round(max(d$ENLC), 1), round(max(d$T1R), 1),
                         d$time[which.max(d$T1R)]))
  })

  ## ---- 6. nerve ------------------------------------------------------
  output$p_nerve <- renderPlot({
    d <- sim()
    a <- d %>% transmute(time,
                         `Reversible pool (drug target)` = NFIR,
                         `Permanent axonal loss (sink)`  = NFIP,
                         `Total impairment`              = NFI) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(a, aes(time, y, colour = v)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(PAL[3], PAL[2], PAL[6]), name = NULL) +
      labs(x = "Day", y = "Impairment (points, 0-100)",
           title = "Nerve function: a reversible pool draining into an irreversible sink",
           subtitle = "k_fix = 1/180 d⁻¹. Prednisolone acts on the pool, never on the sink.") +
      THEME
  })
  output$p_window <- renderPlot({
    delays <- c(0, 30, 60, 90, 120, 180, 270, 365)
    base <- events()
    base <- base[!(base$cmt == "PDNG"), , drop = FALSE]
    res <- sapply(delays, function(dl) {
      ev <- rbind(base, pred_taper(20 + dl, weeks = max(input$pdn_weeks, 52),
                                   floor = 10))
      o <- lep_mod %>% param(pars()) %>%
        mrgsim(events = make_ev(ev), end = 1200, delta = 4)
      tail(as.data.frame(o)$NFIP, 1)
    })
    ggplot(data.frame(delay = delays, nfip = res), aes(delay, nfip)) +
      geom_line(colour = PAL[2], linewidth = 1) + geom_point(size = 2) +
      labs(x = "Delay in starting prednisolone (days after reaction onset)",
           y = "Permanent deficit at day 1200 (points)",
           title = "The window",
           subtitle = "Half of the salvageable deficit is lost in the first ~117 days") +
      THEME
  })

  ## ---- 7. scenarios --------------------------------------------------
  scen_data <- reactive({
    req(length(input$scen) > 0)
    bind_rows(lapply(input$scen, run_scenario))
  })
  output$p_scen <- renderPlot({
    d <- scen_data()
    a <- d %>% select(time, scenario, BI, BILIV, ENL, NFI) %>%
      pivot_longer(c(BI, BILIV, ENL, NFI), names_to = "v", values_to = "y") %>%
      mutate(v = recode(v, BI = "Bacterial index",
                        BILIV = "Viable burden (log10/g)",
                        ENL = "ENL score", NFI = "Nerve impairment"))
    ggplot(a, aes(time, y, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~v, scales = "free_y") +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "Day", y = NULL, colour = NULL,
           title = "Scenario comparison") + THEME
  })
  output$t_scen <- DT::renderDataTable({
    d <- scen_data()
    DT::datatable(bind_rows(lapply(split(d, d$scenario), summarise_scenario)),
                  options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 8. relapse ----------------------------------------------------
  output$p_relapse <- renderPlot({
    ms <- c(3, 6, 9, 12, 18, 24)
    res <- sapply(ms, function(m) {
      ev <- mdt_mb(m)
      o <- lep_mod %>% param(pars()) %>%
        mrgsim(events = make_ev(ev), end = 30 * m + 1, delta = 5)
      tail(as.data.frame(o)$PREL, 1)
    })
    ggplot(data.frame(months = ms, risk = res), aes(months, risk)) +
      geom_line(colour = PAL[1], linewidth = 1) + geom_point(size = 2) +
      labs(x = "Duration of MDT (months)", y = "Relapse escape index (%)",
           title = "Residual viable burden translated into a relapse index",
           subtitle = "Poisson single-escape from the viable burden at the end of treatment") +
      THEME
  })
  output$t_relapse <- renderTable({
    d <- sim()
    i <- which.min(abs(d$time - 30 * input$months))
    data.frame(Quantity = c("Duration of MDT (months)",
                            "Viable burden at end of treatment (log10/g)",
                            "Viable organisms remaining",
                            "Relapse escape index (%)"),
               Value = c(input$months, round(d$BILIV[i], 2),
                         signif(d$NVIAB[i], 3), signif(d$PREL[i], 3)))
  })

  ## ---- 9. safety -----------------------------------------------------
  output$p_safety <- renderPlot({
    d <- sim()
    ggplot(long(d, c("METHB", "HB", "PIG", "HPA", "ALT", "CLOE"),
                c("Methaemoglobin (%)", "Haemoglobin (g/dL)",
                  "Clofazimine pigmentation (0-100)", "HPA axis integrity (-)",
                  "ALT (U/L)", "Clofazimine tissue depot (x standard)")),
           aes(time, y)) +
      geom_line(colour = PAL[5], linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      labs(x = "Day", y = NULL, title = "Treatment-emergent toxicity") + THEME
  })

  ## ---- 10. model -----------------------------------------------------
  output$txt_struct <- renderText(paste(
    "38 ODE compartments, 17 scenarios.",
    "",
    "  BI = log10(live + dead per gram)        <- clearance clock, t1/2 137 d",
    "  MI = 100 x live / (live + dead)         <- kill clock",
    "  antigen liberation = YB*death_flux + YD*k_clr*B_dead",
    "  k_clr = KCLR0 * (1 + CLRBOOST*(1 - B_live/(B_live+KSUB)))",
    "  ENL   = Hill(TNF, neutrophils), driven by immune complexes = Ag x Ab",
    "  T1R   = CMI x f(intraneural antigen), CMI released from anergy by",
    "          the fall in LIVE burden",
    "  nerve : dNFIR/dt = damage - recovery - k_fix*NFIR ; dNFIP/dt = k_fix*NFIR",
    "  containment threshold SPEC* = (MUMAX - KNAT)/KHOST = 0.0225",
    "",
    "All 38 ODEs are re-implemented in lep_verify_python.py and checked",
    "against 54 published anchors (54/54 pass).",
    sep = "\n"))
  output$t_par <- DT::renderDataTable({
    p <- as.list(param(lep_mod))
    DT::datatable(data.frame(parameter = names(p),
                             value = unlist(p), row.names = NULL),
                  options = list(pageLength = 25), rownames = FALSE)
  })
  output$txt_caveat <- renderText(paste(
    "1. The clofazimine loading result (49% reduction in cumulative ENL",
    "   burden from one month of 300 mg/d) is the model's most exposed",
    "   claim. No trial has run it.",
    "2. Relapse is a Poisson index computed from the residual viable",
    "   burden, not a probability read off the deterministic trajectory;",
    "   the continuum ODEs cannot extinguish a population.",
    "3. The dead-bacillus degradation rate and its de-repression by killing",
    "   are inferred from BI decline data, not measured directly.",
    "4. Nerve impairment is a 0-100 index calibrated to the shape of the",
    "   WHO disability grades, not to individual nerve conduction studies.",
    "5. Educational and research model. Not for clinical use.",
    sep = "\n"))
}

shinyApp(ui, server)
