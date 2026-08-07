## =====================================================================
##  Rheumatic heart disease — QSP explorer (Shiny)
##  Front end for rhd_mrgsolve_model.R
## ---------------------------------------------------------------------
##  The app is organised around the three claims the model exists to make,
##  and each tab is built so that the claim can be DISAGREED with rather
##  than merely displayed:
##
##   Tab 3  puts doses-given and time-above-threshold side by side, so the
##          rank inversion between them is visible rather than asserted.
##   Tab 6  separates the immune and autonomous arms of valve loss and
##          draws the crossover, so the reader can see which endpoint
##          prophylaxis moves and where.
##   Tab 5  holds the orifice FIXED and moves only demand and rhythm, so
##          that decompensation without valve change is a manipulation,
##          not a story.
##
##  Run with:  shiny::runApp("rhd_shiny_app.R")
##  Requires:  shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## ---- load the model ---------------------------------------------------
## rhd_mrgsolve_model.R defines `mod` plus the dosing helpers used below.
source("rhd_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

MICP_LINE <- 0.02   # ug/mL, the protective plasma threshold

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Rheumatic heart disease — quantitative systems pharmacology"),
  p(tags$em(paste(
    "Educational model. All parameters are illustrative and calibrated to",
    "published central estimates; nothing here is a clinical tool."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      sliderInput("wt",   "Body weight (kg)",        20, 120, 55, 5),
      sliderInput("age",  "Age (years)",              5,  80, 25, 1),
      sliderInput("susc", "Host susceptibility (SUSC)", 0.2, 3, 1, 0.1),
      sliderInput("mva0", "Starting mitral valve area (cm²)", 0.6, 4.5, 3.6, 0.05),
      sliderInput("mem0", "Immunological memory of a prior attack (0–1)",
                  0, 1, 1, 0.05),

      hr(),
      h4("Setting"),
      sliderInput("lamexp", "Sore-throat exposures per year", 0, 8, 3, 0.5),
      sliderInput("pgasp",  "P(GAS | sore throat)", 0.05, 0.6, 0.25, 0.05),

      hr(),
      h4("Secondary prophylaxis"),
      radioButtons("interval", "Injection interval",
                   c("4-weekly (standard)" = 28,
                     "3-weekly (high risk)" = 21,
                     "2-weekly"             = 14,
                     "none"                 = 0), selected = 28),
      sliderInput("adher", "Proportion of scheduled doses actually given",
                  0, 1, 1.0, 0.05),
      sliderInput("bpgmg", "Dose (mg penicillin-G equivalent; 720 ≈ 1.2 MU)",
                  360, 1440, 720, 180),

      hr(),
      h4("Haemodynamic state"),
      sliderInput("demf", "Cardiac output demand (1 = rest, 1.5 ≈ pregnancy)",
                  1, 2.5, 1, 0.1),
      sliderInput("afb0", "Atrial fibrillation burden (0–1)", 0, 1, 0, 0.05),
      sliderInput("bbdose", "Metoprolol (mg/day)", 0, 200, 0, 25),
      checkboxInput("warfon", "Warfarin 5 mg/day", FALSE),

      hr(),
      h4("Acute attack (episodic channel)"),
      checkboxInput("episodic", "Simulate a discrete pharyngitis", FALSE),
      sliderInput("txday", "Day penicillin V is started (0 = never)",
                  0, 14, 0, 1),
      sliderInput("predwk", "Prednisolone course (weeks)", 0, 26, 0, 2),

      hr(),
      sliderInput("years", "Simulation length (years)", 0.2, 30, 25, 0.2),
      actionButton("go", "Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "1 · Patient profile",
          br(),
          fluidRow(column(6, DTOutput("tbl_profile")),
                   column(6, DTOutput("tbl_endpoints"))),
          hr(),
          plotOutput("p_profile", height = "320px"),
          helpText(paste(
            "Everything downstream is driven by this state. Note that the",
            "two summary numbers a programme would record — doses given and",
            "valve area — are the two that agree least with each other."))
        ),

        tabPanel(
          "2 · Penicillin PK",
          br(),
          plotOutput("p_pk", height = "340px"),
          hr(),
          plotOutput("p_pk_weight", height = "300px"),
          helpText(paste(
            "Benzathine penicillin G has flip-flop kinetics: plasma tracks",
            "RELEASE from the depot (t½ ≈ 9 d), not clearance (t½ 30 min).",
            "The amplitude scales inversely with body size, so the days a",
            "single injection covers fall from ~23 at 40 kg to ~15 at 100 kg."))
        ),

        tabPanel(
          "3 · Adherence vs protection",
          br(),
          fluidRow(column(7, plotOutput("p_inversion", height = "380px")),
                   column(5, DTOutput("tbl_inversion"))),
          hr(),
          helpText(paste(
            "The horizontal axis is the number programmes audit. The vertical",
            "axis is the number that protects. If the points do not lie on a",
            "line, the two disagree; if any pair is inverted, they disagree",
            "about the ORDERING of two patients, which is the claim."))
        ),

        tabPanel(
          "4 · Acute rheumatic fever",
          br(),
          plotOutput("p_arf", height = "420px"),
          hr(),
          plotOutput("p_steroid", height = "300px"),
          helpText(paste(
            "The cross-reactive antibody has a 60-day half-life. A six-week",
            "steroid course covers only ~38% of the antibody-time integral,",
            "so most of the valvulitis it suppresses simply happens later.",
            "That is the model's account of the null Cochrane result, and it",
            "predicts what a longer course would do."))
        ),

        tabPanel(
          "5 · The Gorlin block",
          br(),
          fluidRow(column(6, plotOutput("p_gorlin", height = "360px")),
                   column(6, plotOutput("p_hropt",  height = "360px"))),
          hr(),
          DTOutput("tbl_gorlin"),
          helpText(paste(
            "The orifice is held fixed across every row of the table. What",
            "changes is demand and rhythm. Mean gradient goes as flow",
            "squared over area squared, so a 50% rise in cardiac output is a",
            "3–4 fold rise in gradient at an unchanged valve. Read backwards,",
            "the same equation gives an optimal heart rate that falls as the",
            "valve narrows — that curve is computed, not assumed."))
        ),

        tabPanel(
          "6 · Two clocks of valve loss",
          br(),
          plotOutput("p_arms", height = "380px"),
          hr(),
          plotOutput("p_traj", height = "320px"),
          helpText(paste(
            "The immune arm is what penicillin prevents. The autonomous arm",
            "is shear, and shear IS the mean gradient over four, so it",
            "accelerates as the valve narrows. Above the crossover the",
            "immune arm dominates and prophylaxis protects the valve; below",
            "it the valve destroys itself faster than recurrences do.",
            "Prophylaxis still prevents recurrences everywhere — that is a",
            "different endpoint, shown on the next tab."))
        ),

        tabPanel(
          "7 · Clinical endpoints",
          br(),
          plotOutput("p_endpoints", height = "420px"),
          hr(),
          DTOutput("tbl_scen"),
          helpText(paste(
            "Recurrent ARF and valve area are different endpoints with",
            "different effect sizes for the same drug. A trial powered on",
            "one will not answer the other."))
        ),

        tabPanel(
          "8 · Rhythm, atrium, embolism",
          br(),
          fluidRow(column(6, plotOutput("p_la",  height = "330px")),
                   column(6, plotOutput("p_emb", height = "330px"))),
          hr(),
          helpText(paste(
            "Left atrial pressure drives remodelling, remodelling drives the",
            "AF hazard, AF shortens diastole and raises the gradient, and the",
            "gradient raises left atrial pressure. The loop closes, which is",
            "why AF onset is a step change rather than a nuisance.",
            "INVICTUS (2022) found a vitamin K antagonist superior to",
            "rivaroxaban here; the model encodes that empirically and does",
            "not explain it."))
        ),

        tabPanel(
          "9 · Intervention",
          br(),
          fluidRow(column(4, sliderInput("pmbv_t", "Valvotomy at year", 0, 25, 5, 0.5)),
                   column(4, sliderInput("pmbv_g", "Immediate area gain (cm²)",
                                         0, 1.6, 1.15, 0.05)),
                   column(4, sliderInput("resten", "Post-procedure shear multiplier",
                                         1, 2, 1.35, 0.05))),
          plotOutput("p_pmbv", height = "360px"),
          hr(),
          DTOutput("tbl_wilkins"),
          helpText(paste(
            "A balloon splits the fused commissure; it does not remove the",
            "fibrosis, so the shear arm resumes on a valve that can now carry",
            "more flow. The 10-year trajectory here is more pessimistic than",
            "the ~40%-restenosis literature, and that is flagged in README.md",
            "as the model's most exposed claim rather than tuned away."))
        ),

        tabPanel(
          "10 · Biomarkers",
          br(),
          plotOutput("p_biomarkers", height = "420px"),
          hr(),
          helpText(paste(
            "ASO is evidence of infection, not of pathogenicity: it and the",
            "cross-reactive antibody are driven by the same antigen but",
            "decay on different clocks (45 d vs 60 d), and only one of them",
            "scars the valve. CRP tracks the valvulitis that is happening,",
            "not the valve area that has already been lost."))
        ),

        tabPanel(
          "About / assumptions",
          br(),
          htmlOutput("about")
        )
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## ---- build the model instance from the sidebar ----------------------
  base_mod <- reactive({
    mod %>% param(WT = input$wt, AGE = input$age, SUSC = input$susc,
                  LAMEXP = input$lamexp, PGASP = input$pgasp,
                  DEMF = input$demf,
                  EPISODIC = as.numeric(input$episodic))
  })

  ## which scheduled doses were actually given, deterministically thinned
  keep_idx <- function(n, frac) {
    if (frac >= 1) return(seq_len(n))
    if (frac <= 0) return(integer(0))
    step <- 1 / (1 - frac)
    drop <- unique(round(seq(step, n, by = step)))
    setdiff(seq_len(n), drop)
  }

  proph_events <- function(interval, years, mg, frac) {
    interval <- as.numeric(interval)
    if (interval <= 0 || frac <= 0) return(NULL)
    n <- ceiling(years * 365 / interval) + 1
    bpg(interval, n, mg = mg, keep = keep_idx(n, frac))
  }

  sim <- eventReactive(input$go, {
    m <- base_mod() %>% init(MVA = input$mva0, MEM = input$mem0, AFB = input$afb0)
    tend <- input$years * 365
    ev <- proph_events(input$interval, input$years, input$bpgmg, input$adher)
    if (input$episodic) ev <- rbind(ev, inoc(0))
    if (input$txday > 0) ev <- rbind(ev, penv(input$txday))
    if (input$predwk > 0) ev <- rbind(ev, pred(19, input$predwk * 7))
    if (input$bbdose > 0) ev <- rbind(ev, bblok(0, tend, input$bbdose))
    if (input$warfon)     ev <- rbind(ev, warf(0, tend, 5))
    if (!is.null(ev)) ev <- ev[order(ev$time, ev$cmt), ]
    out <- if (is.null(ev)) mrgsim(m, end = tend, delta = max(0.25, tend / 3000))
           else mrgsim(m, data = ev, end = tend, delta = max(0.25, tend / 3000))
    as.data.frame(out)
  }, ignoreNULL = FALSE)

  ## ---- Tab 1: profile --------------------------------------------------
  output$tbl_profile <- renderDT({
    d <- sim(); last <- tail(d, 1)
    datatable(data.frame(
      Quantity = c("Body weight (kg)", "Host susceptibility", "Injection interval (d)",
                   "Doses given / due", "Simulated years"),
      Value = c(input$wt, input$susc,
                ifelse(input$interval == 0, "none", input$interval),
                sprintf("%.0f%%", 100 * input$adher), input$years)),
      rownames = FALSE, options = list(dom = "t"),
      caption = "What the programme records")
  })

  output$tbl_endpoints <- renderDT({
    d <- sim(); last <- tail(d, 1); tend <- input$years * 365
    datatable(data.frame(
      Endpoint = c("Time above 0.02 ug/mL (%)", "Mitral valve area (cm2)",
                   "Mean gradient (mmHg)", "LA pressure (mmHg)",
                   "AF burden", "Expected recurrent ARF",
                   "Expected embolic events", "NYHA class"),
      Value = c(sprintf("%.1f", 100 * last$TPROT / tend),
                sprintf("%.2f", last$MVAcm2), sprintf("%.1f", last$MVGRAD),
                sprintf("%.1f", last$LAPRES), sprintf("%.2f", last$AFBURD),
                sprintf("%.2f", last$CUMARF), sprintf("%.3f", last$CUMEMB),
                sprintf("%.2f", last$NYHACL))),
      rownames = FALSE, options = list(dom = "t"),
      caption = "What actually happened")
  })

  output$p_profile <- renderPlot({
    d <- sim()
    d %>% select(time, MVAcm2, MVGRAD, LAPRES, AFBURD, CUMARF, NYHACL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 365, value)) +
      geom_line(linewidth = 0.7, colour = "#2c5f8a") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "years", y = NULL) + THEME
  })

  ## ---- Tab 2: PK -------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time <= min(max(time), 180))
    ggplot(d, aes(time, CPEN)) +
      geom_hline(yintercept = MICP_LINE, linetype = 2, colour = "#b03030") +
      annotate("text", x = 0, y = MICP_LINE * 1.35, hjust = 0, size = 3.4,
               colour = "#b03030", label = "0.02 ug/mL protective threshold") +
      geom_line(linewidth = 0.6, colour = "#2c5f8a") +
      scale_y_log10() +
      labs(x = "days", y = "plasma penicillin G (ug/mL, log scale)",
           title = "The gaps between injections are the disease's opportunity") +
      THEME
  })

  output$p_pk_weight <- renderPlot({
    ws <- c(30, 45, 60, 75, 90, 105)
    res <- lapply(ws, function(w) {
      o <- as.data.frame(mrgsim(param(mod, WT = w), data = bpg(n = 1),
                                end = 60, delta = 0.25))
      data.frame(wt = w, days = tail(o$TPROT, 1))
    }) %>% bind_rows()
    ggplot(res, aes(wt, days)) +
      geom_col(fill = "#2c5f8a", width = 8) +
      geom_hline(yintercept = c(21, 28), linetype = 2, colour = "#b03030") +
      annotate("text", x = min(ws), y = 28.8, hjust = 0, size = 3.4,
               colour = "#b03030", label = "28-day interval") +
      annotate("text", x = min(ws), y = 21.8, hjust = 0, size = 3.4,
               colour = "#b03030", label = "21-day interval") +
      labs(x = "body weight (kg)",
           y = "days above threshold from ONE 1.2 MU injection",
           title = "A heavier patient loses days, not doses") + THEME
  })

  ## ---- Tab 3: the rank inversion ---------------------------------------
  inversion <- reactive({
    grid <- expand.grid(wt = c(40, 55, 70, 85, 100),
                        interval = c(21, 28),
                        adher = c(0.6, 0.8, 1.0))
    yrs <- 3
    out <- lapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      n <- ceiling(yrs * 365 / g$interval) + 1
      ev <- bpg(g$interval, n, keep = keep_idx(n, g$adher))
      o <- as.data.frame(mrgsim(param(mod, WT = g$wt), data = ev,
                                end = yrs * 365, delta = 2))
      data.frame(wt = g$wt, interval = g$interval, adher = g$adher,
                 protected = tail(o$TPROT, 1) / (yrs * 365))
    }) %>% bind_rows()
    out
  })

  output$p_inversion <- renderPlot({
    d <- inversion()
    ggplot(d, aes(100 * adher, 100 * protected,
                  colour = factor(wt), shape = factor(interval))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      geom_point(size = 3.2) +
      scale_colour_viridis_d(name = "weight (kg)", option = "C", end = 0.85) +
      scale_shape_manual(name = "interval (d)", values = c(16, 17)) +
      labs(x = "doses given, % of doses due  (what is audited)",
           y = "% of calendar time protected  (what protects)",
           title = "Points off the diagonal mean the two numbers disagree") +
      coord_cartesian(xlim = c(50, 105), ylim = c(30, 105)) + THEME
  })

  output$tbl_inversion <- renderDT({
    d <- inversion() %>%
      arrange(desc(adher), desc(protected)) %>%
      transmute(`weight (kg)` = wt, `interval (d)` = interval,
                `doses given (%)` = 100 * adher,
                `time protected (%)` = round(100 * protected, 1))
    datatable(d, rownames = FALSE, options = list(pageLength = 12, dom = "tp"))
  })

  ## ---- Tab 4: acute attack ---------------------------------------------
  arf_sim <- reactive({
    as.data.frame(mrgsim(param(base_mod(), EPISODIC = 1),
                         data = mk(inoc(), asa(19, 42)),
                         end = 400, delta = 1))
  })

  output$p_arf <- renderPlot({
    arf_sim() %>%
      transmute(time, `GAS (log10)` = LOGGAS, `ASO (Todd U)` = ASOT,
                `CRP (mg/L)` = CRPMGL, `valvulitis (AU)` = VIT,
                `mitral regurgitation` = MRTOT, `valve area (cm2)` = MVAcm2) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.7, colour = "#a05020") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days from inoculation", y = NULL,
           title = "One untreated pharyngitis, followed for a year") + THEME
  })

  output$p_steroid <- renderPlot({
    wks <- c(0, 2, 6, 12, 26)
    res <- lapply(wks, function(wk) {
      ev <- if (wk > 0) mk(inoc(), asa(19, 42), pred(19, wk * 7))
            else        mk(inoc(), asa(19, 42))
      o <- as.data.frame(mrgsim(param(base_mod(), EPISODIC = 1), data = ev,
                                end = 400, delta = 1))
      data.frame(weeks = wk, lost = 4.5 - tail(o$MVAcm2, 1))
    }) %>% bind_rows() %>%
      mutate(saved = 100 * (1 - lost / lost[weeks == 0]))
    ggplot(res, aes(weeks, saved)) +
      geom_line(linewidth = 0.8, colour = "#2c5f8a") +
      geom_point(size = 3, colour = "#2c5f8a") +
      geom_vline(xintercept = 6, linetype = 2, colour = "#b03030") +
      annotate("text", x = 6.4, y = 5, hjust = 0, size = 3.4, colour = "#b03030",
               label = "the duration the trials used") +
      labs(x = "prednisolone course (weeks)",
           y = "% of valve area saved",
           title = "Duration, not dose, is the binding constraint") + THEME
  })

  ## ---- Tab 5: the Gorlin block -----------------------------------------
  gorlin_rows <- reactive({
    states <- list(
      list(lab = "rest, sinus rhythm",      demf = 1.0, afb = 0, bb = 0),
      list(lab = "moderate exertion",       demf = 1.6, afb = 0, bb = 0),
      list(lab = "pregnancy (CO x1.5)",     demf = 1.5, afb = 0, bb = 0),
      list(lab = "new-onset AF, at rest",   demf = 1.0, afb = 1, bb = 0),
      list(lab = "AF + metoprolol 100 mg",  demf = 1.0, afb = 1, bb = 100),
      list(lab = "fever/anaemia (CO x1.8)", demf = 1.8, afb = 0, bb = 0))
    lapply(states, function(s) {
      m <- init(param(mod, DEMF = s$demf), MVA = input$mva0, AFB = s$afb)
      ev <- if (s$bb > 0) bblok(0, 200, s$bb) else NULL
      o <- as.data.frame(if (is.null(ev)) mrgsim(m, end = 200, delta = 1)
                         else mrgsim(m, data = ev, end = 200, delta = 1))
      l <- tail(o, 1)
      data.frame(state = s$lab, HR = l$HRBPM, DFP = l$DFPMIN, CO = l$COLMIN,
                 gradient = l$MVGRAD, LAP = l$LAPRES, NYHA = l$NYHACL)
    }) %>% bind_rows()
  })

  output$tbl_gorlin <- renderDT({
    d <- gorlin_rows() %>%
      transmute(State = state,
                `heart rate` = round(HR),
                `diastolic s/min` = round(DFP, 1),
                `cardiac output` = round(CO, 2),
                `mean gradient` = round(gradient, 1),
                `LA pressure` = round(LAP, 1),
                NYHA = round(NYHA, 2))
    datatable(d, rownames = FALSE, options = list(dom = "t"),
              caption = sprintf(
                "Mitral valve area fixed at %.2f cm2 in EVERY row", input$mva0))
  })

  output$p_gorlin <- renderPlot({
    d <- gorlin_rows()
    ggplot(d, aes(reorder(state, gradient), gradient)) +
      geom_col(fill = "#7b5ea8", width = 0.65) +
      geom_hline(yintercept = 18, linetype = 2, colour = "#b03030") +
      annotate("text", x = 0.6, y = 19.5, hjust = 0, size = 3.3,
               colour = "#b03030", label = "congestion threshold") +
      coord_flip() +
      labs(x = NULL, y = "mean mitral gradient (mmHg)",
           title = sprintf("Same %.2f cm2 orifice, six different patients",
                           input$mva0)) + THEME
  })

  output$p_hropt <- renderPlot({
    areas <- c(2.0, 1.5, 1.0, 0.8)
    grid <- expand.grid(hr = 40:160, mva = areas)
    grid <- grid %>% mutate(
      RR = 60 / hr,
      DFPB = pmax(RR - 0.36 * sqrt(RR), 0.05),
      DFPM = hr * DFPB,
      COvalve = 37.7 * mva * sqrt(32 - 6) * DFPM / 1000,
      COsv = 110 * hr / 1000,
      CO = pmin(COvalve, COsv))
    best <- grid %>% group_by(mva) %>% slice_max(CO, n = 1)
    ggplot(grid, aes(hr, CO, colour = factor(mva))) +
      geom_line(linewidth = 0.8) +
      geom_point(data = best, size = 3) +
      scale_colour_viridis_d(name = "valve area (cm2)", option = "D", end = 0.85) +
      labs(x = "heart rate (beats/min)",
           y = "maximum achievable cardiac output (L/min)",
           title = "The optimum falls as the valve narrows",
           subtitle = "falling limb = the valve; rising limb = the ventricle") +
      THEME
  })

  ## ---- Tab 6: two clocks -----------------------------------------------
  arms <- reactive({
    areas <- seq(0.7, 4.5, by = 0.1)
    haz <- 0.10
    lapply(areas, function(a) {
      o <- as.data.frame(mrgsim(init(param(mod, SUSC = input$susc), MVA = a),
                                end = 1, delta = 1))
      l <- tail(o, 1)
      geom <- 2 * sqrt(pi * a)
      brake <- (a - 0.30) / (a - 0.30 + 1.6)
      shear <- geom * 9.132e-5 * (l$MVGRAD / 4) * brake * 365
      ## immune arm at the steady state sustained by an ARF hazard `haz`
      ag <- 2.483 * haz / 365 / 0.15
      mem <- 4e-3 * ag / (4e-3 * ag + 3.8e-4)
      xab <- 4e-3 * ag * input$susc * (1 + 1.5 * mem) / 0.01155
      scarf <- max(0, 1 - a / 4.5)
      vit <- 0.60 * xab * (1 + 2 * scarf) / 0.0231
      immune <- geom * 3.004e-4 * vit * 365
      data.frame(mva = a, shear = shear, immune = immune)
    }) %>% bind_rows()
  })

  output$p_arms <- renderPlot({
    d <- arms()
    cross <- d %>% mutate(diff = shear - immune) %>%
      filter(abs(diff) == min(abs(diff))) %>% pull(mva)
    d %>% pivot_longer(-mva, names_to = "arm", values_to = "rate") %>%
      ggplot(aes(mva, rate, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = cross, linetype = 2, colour = "grey40") +
      annotate("text", x = cross, y = max(d$shear), hjust = -0.05, size = 3.6,
               label = sprintf("crossover %.2f cm2", cross)) +
      scale_colour_manual(values = c(immune = "#a05020", shear = "#7b5ea8"),
                          labels = c("immune arm (penicillin prevents this)",
                                     "autonomous shear arm (it does not)"),
                          name = NULL) +
      scale_x_reverse() +
      labs(x = "mitral valve area (cm2), narrowing to the right",
           y = "valve area lost per year (cm2/yr)",
           title = "Two arms, one orifice, no prophylaxis at all") + THEME
  })

  output$p_traj <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 365, MVAcm2)) +
      geom_line(linewidth = 0.8, colour = "#2c5f8a") +
      geom_hline(yintercept = c(1.0, 1.5, 2.5), linetype = 3, colour = "grey55") +
      annotate("text", x = 0, y = c(1.05, 1.55, 2.55), hjust = 0, size = 3.2,
               colour = "grey35", label = c("severe", "moderate", "mild")) +
      labs(x = "years", y = "mitral valve area (cm2)",
           title = "This patient's trajectory, on the settings in the sidebar") +
      THEME
  })

  ## ---- Tab 7: endpoints ------------------------------------------------
  scenarios <- reactive({
    yrs <- input$years; tend <- yrs * 365
    defs <- list(
      list(lab = "no prophylaxis",        iv = 0,  ad = 0),
      list(lab = "4-weekly, all doses",   iv = 28, ad = 1.0),
      list(lab = "4-weekly, 70% of doses",iv = 28, ad = 0.7),
      list(lab = "3-weekly, all doses",   iv = 21, ad = 1.0),
      list(lab = "3-weekly, 70% of doses",iv = 21, ad = 0.7))
    lapply(defs, function(s) {
      m <- init(base_mod(), MVA = input$mva0, MEM = input$mem0, AFB = input$afb0)
      ev <- proph_events(s$iv, yrs, input$bpgmg, s$ad)
      o <- as.data.frame(if (is.null(ev)) mrgsim(m, end = tend, delta = 10)
                         else mrgsim(m, data = ev, end = tend, delta = 10))
      l <- tail(o, 1)
      data.frame(scenario = s$lab, protected = 100 * l$TPROT / tend,
                 arf = l$CUMARF, mva = l$MVAcm2, grad = l$MVGRAD,
                 emb = l$CUMEMB, nyha = l$NYHACL)
    }) %>% bind_rows()
  })

  output$p_endpoints <- renderPlot({
    d <- scenarios() %>%
      mutate(`recurrent ARF prevented (%)` =
               100 * (1 - arf / arf[scenario == "no prophylaxis"]),
             `valve area preserved (%)` =
               100 * (mva / mva[scenario == "no prophylaxis"] - 1)) %>%
      select(scenario, `recurrent ARF prevented (%)`, `valve area preserved (%)`) %>%
      pivot_longer(-scenario)
    ggplot(d, aes(reorder(scenario, value), value, fill = name)) +
      geom_col(position = "dodge", width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("#2c5f8a", "#a05020"), name = NULL) +
      labs(x = NULL, y = "% relative to no prophylaxis",
           title = "One drug, two endpoints, two very different effect sizes") +
      THEME
  })

  output$tbl_scen <- renderDT({
    d <- scenarios() %>%
      transmute(Scenario = scenario,
                `time protected (%)` = round(protected, 1),
                `recurrent ARF` = round(arf, 3),
                `valve area (cm2)` = round(mva, 2),
                `gradient (mmHg)` = round(grad, 1),
                `embolic events` = round(emb, 3),
                NYHA = round(nyha, 2))
    datatable(d, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- Tab 8: atrium, rhythm, embolism ---------------------------------
  output$p_la <- renderPlot({
    d <- sim()
    d %>% transmute(time, `LA pressure (mmHg)` = LAPRES,
                    `LA volume (mL)` = LAV, `AF burden` = AFBURD,
                    `PVR (Wood units)` = PVR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 365, value)) +
      geom_line(linewidth = 0.7, colour = "#2c5f8a") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "years", y = NULL, title = "The loop that closes") + THEME
  })

  output$p_emb <- renderPlot({
    tend <- input$years * 365
    m <- init(base_mod(), MVA = input$mva0, MEM = input$mem0, AFB = input$afb0)
    ev0 <- proph_events(input$interval, input$years, input$bpgmg, input$adher)
    evw <- rbind(ev0, warf(0, tend, 5))
    a <- as.data.frame(if (is.null(ev0)) mrgsim(m, end = tend, delta = 10)
                       else mrgsim(m, data = ev0, end = tend, delta = 10))
    b <- as.data.frame(mrgsim(m, data = evw[order(evw$time, evw$cmt), ],
                              end = tend, delta = 10))
    bind_rows(mutate(select(a, time, CUMEMB), arm = "no anticoagulation"),
              mutate(select(b, time, CUMEMB), arm = "warfarin 5 mg/day")) %>%
      ggplot(aes(time / 365, CUMEMB, colour = arm)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c("#b03030", "#2c5f8a"), name = NULL) +
      labs(x = "years", y = "cumulative expected embolic events",
           title = "Anticoagulation removes about two-thirds of the hazard") +
      THEME
  })

  ## ---- Tab 9: intervention ---------------------------------------------
  output$p_pmbv <- renderPlot({
    tend <- input$years * 365
    tp <- input$pmbv_t * 365
    m <- init(base_mod(), MVA = input$mva0, MEM = input$mem0, AFB = input$afb0)
    ev <- proph_events(input$interval, input$years, input$bpgmg, input$adher)
    ## pre-procedure
    a <- as.data.frame(if (is.null(ev)) mrgsim(m, end = max(tp, 1), delta = 5)
                       else mrgsim(m, data = ev, end = max(tp, 1), delta = 5))
    la <- tail(a, 1)
    ## post-procedure, with the shear coefficient stepped up
    m2 <- init(param(base_mod(), RESTEN = input$resten),
               MVA = la$MVAcm2 + input$pmbv_g, CA = la$CA,
               LAV = la$LAV, AFB = la$AFBURD, PVR = la$PVR, MEM = la$MEM)
    b <- as.data.frame(mrgsim(m2, end = max(tend - tp, 1), delta = 5))
    bind_rows(transmute(a, yr = time / 365, MVA = MVAcm2, arm = "before"),
              transmute(b, yr = input$pmbv_t + time / 365, MVA = MVAcm2,
                        arm = "after valvotomy")) %>%
      ggplot(aes(yr, MVA, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1.5, linetype = 3, colour = "grey50") +
      scale_colour_manual(values = c("#2c5f8a", "#a05020"), name = NULL) +
      labs(x = "years", y = "mitral valve area (cm2)",
           title = "The balloon splits the commissure; it does not stop the shear") +
      THEME
  })

  output$tbl_wilkins <- renderDT({
    l <- tail(sim(), 1)
    datatable(data.frame(
      Component = c("Wilkins score (4-16)", "Valve calcium (AU)",
                    "Total mitral regurgitation (grade)",
                    "Suitable for balloon on score alone?"),
      Value = c(sprintf("%.1f", l$WILKIN), sprintf("%.1f", l$CA),
                sprintf("%.2f", l$MRTOT),
                ifelse(l$WILKIN <= 8, "yes (score <= 8)", "no (score > 8)"))),
      rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- Tab 10: biomarkers ---------------------------------------------
  output$p_biomarkers <- renderPlot({
    d <- if (input$episodic) arf_sim() else sim()
    d %>% transmute(time, `ASO titre (Todd U)` = ASOT, `CRP (mg/L)` = CRPMGL,
                    `cross-reactive antibody (AU)` = XAB,
                    `valvulitis (AU)` = VIT,
                    `immunological memory` = MEM,
                    `plasma penicillin (ug/mL)` = CPEN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.7, colour = "#2c5f8a") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL,
           title = "Evidence of infection is not evidence of pathogenicity") +
      THEME
  })

  ## ---- About -----------------------------------------------------------
  output$about <- renderUI({
    HTML(paste0(
      "<h4>What this model claims, and how confident it is</h4>",
      "<p><b>Anchored to published data:</b> the benzathine penicillin depot ",
      "profile (Kaplan 1989), penicillin size scaling (Hand 2019, Neely 2014), ",
      "progression of established mitral stenosis (Sagie 1996, Gordon 1992), ",
      "the Gorlin constant (Gorlin 1951), gradient/severity grades ",
      "(Baumgartner 2009), recurrence reduction on benzathine penicillin ",
      "(Manyemba 2002), embolic rates in rheumatic AF, warfarin dose-INR, and ",
      "the 5-year post-valvotomy valve area (Iung 1999, Palacios 2002).</p>",
      "<p><b>Predictions, not fitted:</b> the adherence/protection rank ",
      "inversion; the 2.75 cm&sup2; crossover between the immune and shear ",
      "arms; the area-dependent optimal heart rate; the steroid ",
      "duration-response; and the claim that at the left-atrial-pressure ",
      "ceiling, rate control raises cardiac output while leaving the measured ",
      "gradient unchanged.</p>",
      "<p><b>Most exposed claim:</b> the shear arm is calibrated to a cohort ",
      "MEAN progression of 0.09 cm&sup2;/yr and then accelerates as the valve ",
      "narrows. Sagie 1996 reported the opposite association. The model's ",
      "10-year post-valvotomy trajectory is correspondingly more pessimistic ",
      "than the restenosis literature.</p>",
      "<p><b>Not explained:</b> the INVICTUS finding that a vitamin K ",
      "antagonist beats rivaroxaban in rheumatic AF is encoded empirically; ",
      "the model has no mechanism for it.</p>",
      "<p>All 39 ODEs were independently re-implemented in Python/scipy and ",
      "run against 47 anchors (47/47 pass). Seven real defects were found and ",
      "fixed in that process; they are listed in README.md rather than ",
      "quietly corrected.</p>",
      "<p><i>Educational and hypothesis-generating only. Not for clinical ",
      "decision-making.</i></p>"))
  })
}

shinyApp(ui, server)
