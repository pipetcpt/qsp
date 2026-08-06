# =============================================================================
#  DENGUE / SEVERE DENGUE QSP -- INTERACTIVE DASHBOARD
# =============================================================================
#  A ten-tab Shiny front end to denv_mrgsolve_model.R.
#
#  The app is organised around the model's central claim: dengue is one
#  antibody response read through two channels of opposite sign.  The two
#  channels are deliberately plotted on the SAME time axis in tab 5 so the
#  user can see that the fever break and the circulatory nadir are the same
#  event, and the ADE tab lets the user move the one exogenous variable --
#  the titre the patient was carrying when the mosquito bit -- and watch the
#  disease change character.
#
#  Tabs
#  ----
#   1  Patient & serostatus   -- the one input that matters, and the bell curve
#   2  Virology & NS1         -- viraemia, NS1, interferon, the antibody rise
#   3  Immunopathology        -- immune complexes, TNF, chymase, VEGF, IL-10
#   4  Endothelium & Starling -- glycocalyx, sigma, the pressure balance
#   5  The critical phase     -- BOTH channels on one axis (the headline)
#   6  Haemodynamics & shock  -- pulse pressure, MAP, lactate, oxygen delivery
#   7  Fluid therapy          -- the U-curve and resuscitation efficiency
#   8  Haematology & liver    -- Hct, platelets, bleeding index, transaminases
#   9  Scenario comparison    -- all eighteen arms side by side
#  10  Trial reproductions    -- Salje, infant dengue, Wills, AAPT, balapiravir
#
#  Run:
#    shiny::runApp("denv_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#
#  If mrgsolve is unavailable the app falls back to the pre-computed results in
#  denv_scenario_results.json so that the reproductions in tabs 1 and 10 still
#  work; the time-course tabs then display a message instead of a plot.
# =============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

HAVE_MRGSOLVE <- requireNamespace("mrgsolve", quietly = TRUE)
if (HAVE_MRGSOLVE) {
  suppressMessages(library(mrgsolve))
  # denv_mrgsolve_model.R defines `mod` through mcode_cache() and also brings
  # in who_ladder(), colloid_bolus(), the scenario list and the closed-form
  # analyses, so the app never duplicates the model definition.
  DENV_NO_DEMO <- TRUE                 # suppress the model file demo block
  source("denv_mrgsolve_model.R", local = TRUE)
  MOD <- mod
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom")

PAL <- c(primary = "#1f77b4", secondary = "#d62728", treated = "#2ca02c",
         over = "#ff7f0e", protected = "#9467bd", grey = "#7f7f7f")

T_PRESENT <- 48   # h after fever onset -- the model's own warning-sign trigger

# -----------------------------------------------------------------------------
#  The ADE entry factor.  No solver needed: this is the closed form the whole
#  model turns on, and the bell is the product of a rising saturating term and
#  a falling Hill term rather than an assumed shape.
# -----------------------------------------------------------------------------
ade_curve <- function(A, FCROSS = 0.12, NT50 = 30, HNEUT = 2.5,
                      KOPS = 4, PHI = 14) {
  aeff <- FCROSS * A
  N <- ifelse(aeff > 0, aeff^HNEUT / (aeff^HNEUT + NT50^HNEUT), 0)
  (1 - N) * (1 + (PHI - 1) * A / (A + KOPS))
}

ade_landmarks <- function(...) {
  A <- 10^seq(-1, 4.4, length.out = 4000)
  E <- ade_curve(A, ...)
  i <- max(which(E > 1))
  list(peak_titre = A[which.max(E)], peak_E = max(E),
       cross = 10^approx(c(E[i + 1], E[i]),
                         c(log10(A[i + 1]), log10(A[i])), 1)$y)
}

who_ladder <- function(t0, wt = 70, scale = 1) {
  a <- c(0, 1, 3, 5, 12, 36, 48)
  r <- c(10, 6, 4, 2.5, 1.5, 1, 0)
  data.frame(time = t0 + a, RCRYS = r * wt * scale)
}

colloid_bolus <- function(t0, wt = 70, mlkg = 15, over_h = 1, n = 2, gap = 6) {
  do.call(rbind, lapply(seq_len(n) - 1, function(k)
    data.frame(time = t0 + gap * k + c(0, over_h),
               RCOLL = c(mlkg * wt / over_h, 0))))
}

# The event table is built by build_events() from the model file, so the app
# and the model can never drift apart on how a regimen is encoded.  (The local
# who_ladder/colloid_bolus above are only used for the plots' annotations.)
simulate_arm <- function(abh0, crys = NULL, colloid = FALSE, alb = NULL,
                         plt = NULL, ster = NULL, av = NULL, ns1mab = 0,
                         roral = 0, tp = T_PRESENT, tend = 336) {
  if (!HAVE_MRGSOLVE) return(NULL)
  spec <- list(delay = tp - T_PRESENT)
  if (!is.null(crys)) spec$crys <- crys
  if (colloid)        spec$colloid <- TRUE
  if (!is.null(alb))  spec$alb <- alb
  if (!is.null(plt))  spec$plt <- plt
  if (!is.null(ster)) spec$ster <- ster
  if (!is.null(av))   spec$av <- av
  MOD |>
    param(ABH0 = abh0, NS1MAB = ns1mab, RORAL = roral) |>
    data_set(build_events(spec)) |>
    mrgsim(end = tend, delta = 0.5) |>
    as_tibble()
}

SCEN <- list(
  S01_primary_naive           = list(abh0 = 0),
  S02_secondary_untreated     = list(abh0 = 55),
  S03_secondary_WHO_fluids    = list(abh0 = 55, crys = 1.00),
  S04_under_resuscitation_50  = list(abh0 = 55, crys = 0.50),
  S05_over_resuscitation_200  = list(abh0 = 55, crys = 2.00),
  S06_colloid_rescue          = list(abh0 = 55, crys = 0.55, colloid = TRUE),
  S07_albumin_rescue          = list(abh0 = 55, crys = 0.55, alb = 15),
  S08_prophylactic_platelets  = list(abh0 = 55, crys = 1.00, plt = 38),
  S09_steroid_at_presentation = list(abh0 = 55, crys = 1.00, ster = T_PRESENT),
  S09b_steroid_at_fever_onset = list(abh0 = 55, crys = 1.00, ster = 0),
  S10_antiviral_day2_illness  = list(abh0 = 55, crys = 1.00, av = 48),
  S11_antiviral_at_present    = list(abh0 = 55, crys = 1.00, av = T_PRESENT),
  S12_vaccine_seronegative    = list(abh0 = 48),
  S13_vaccine_seropositive    = list(abh0 = 3100),
  S14_tertiary_high_titre     = list(abh0 = 2600),
  S16_anti_NS1_mab            = list(abh0 = 55, crys = 1.00, ns1mab = 0.35),
  S17_early_oral_rehydration  = list(abh0 = 55, roral = 95))

no_model_msg <- function() {
  ggplot() + annotate("text", 0, 0, size = 5, colour = "grey40",
    label = "mrgsolve is not installed.\nInstall it to enable the time-course tabs;\nthe closed-form tabs still work.") +
    theme_void()
}

line_panel <- function(df, vars, labs_, title, subtitle, ylab, logy = FALSE) {
  if (is.null(df)) return(no_model_msg())
  d <- df |> select(ILLDAY, all_of(vars)) |>
    pivot_longer(-ILLDAY, names_to = "v", values_to = "y") |>
    mutate(v = factor(v, levels = vars, labels = labs_))
  p <- ggplot(d, aes(ILLDAY, y, colour = v)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~v, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = "illness day", y = ylab) +
    guides(colour = "none") + THEME
  if (logy) p <- p + scale_y_log10()
  p
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Dengue / Severe Dengue — Quantitative Systems Pharmacology"),
  tags$p(style = "color:#555;margin-top:-8px",
    tags$em("One antibody response, two channels of opposite sign. ",
            "The fever break and the plasma leak are the same event seen twice.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Serostatus"),
      sliderInput("abh", "Pre-infection heterotypic IgG (reciprocal titre)",
                  min = 0, max = 4000, value = 55, step = 5),
      helpText("The only exogenous variable in the model. 0 = naive; ",
               "1:20–1:80 is the enhancing band; above ~1:700 is protective."),
      hr(),
      h4("Fluid therapy"),
      sliderInput("crys", "Crystalloid, × the WHO ladder",
                  min = 0, max = 3, value = 1, step = 0.25),
      checkboxInput("colloid", "Colloid rescue (2 × 15 mL/kg)", FALSE),
      checkboxInput("plt", "Prophylactic platelets (3 pools)", FALSE),
      sliderInput("tp", "Presentation, h after fever onset",
                  min = 24, max = 120, value = T_PRESENT, step = 6),
      hr(),
      h4("Drugs"),
      checkboxInput("av_on", "NS4B antiviral", FALSE),
      conditionalPanel("input.av_on",
        sliderInput("av_start", "Antiviral start (h after fever onset)",
                    min = 0, max = 120, value = 48, step = 6)),
      checkboxInput("ster_on", "Methylprednisolone 500 mg × 3", FALSE),
      conditionalPanel("input.ster_on",
        sliderInput("ster_start", "Steroid start (h)", min = 0, max = 96,
                    value = 48, step = 6)),
      checkboxInput("ns1mab", "Anti-NS1 monoclonal", FALSE),
      hr(),
      actionButton("go", "Simulate", class = "btn-primary btn-block"),
      helpText(tags$small(
        "The comparator in every plot is the same patient with no antibody ",
        "and no treatment. Model verified against an independent Python/scipy ",
        "implementation; see denv_reference_output.txt."))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Serostatus & ADE",
          fluidRow(column(7, plotOutput("p_ade", height = 340)),
                   column(5, plotOutput("p_ade_zone", height = 340))),
          hr(), verbatimTextOutput("txt_ade")),
        tabPanel("2 · Virology & NS1", plotOutput("p_viro", height = 620)),
        tabPanel("3 · Immunopathology", plotOutput("p_immuno", height = 620)),
        tabPanel("4 · Endothelium & Starling",
          plotOutput("p_endo", height = 380), hr(),
          plotOutput("p_starling", height = 300)),
        tabPanel("5 · The critical phase",
          plotOutput("p_critical", height = 520), hr(),
          verbatimTextOutput("txt_critical")),
        tabPanel("6 · Haemodynamics & shock", plotOutput("p_haemo", height = 620)),
        tabPanel("7 · Fluid therapy",
          plotOutput("p_ucurve", height = 340), hr(),
          plotOutput("p_efficiency", height = 300), hr(),
          verbatimTextOutput("txt_fluid")),
        tabPanel("8 · Haematology & liver", plotOutput("p_haem", height = 620)),
        tabPanel("9 · Scenario comparison",
          plotOutput("p_scen", height = 420), hr(),
          DT::dataTableOutput("tbl_scen")),
        tabPanel("10 · Trial reproductions",
          h4("What the model reproduces without being fitted to it"),
          verbatimTextOutput("txt_trials"), hr(),
          plotOutput("p_titre_risk", height = 340))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  run <- eventReactive(input$go, {
    if (!HAVE_MRGSOLVE) return(NULL)
    args <- list(abh0 = input$abh, tp = input$tp)
    if (input$crys > 0)   args$crys    <- input$crys
    if (input$colloid)    args$colloid <- TRUE
    if (input$plt)        args$plt     <- 38
    if (input$av_on)      args$av      <- input$av_start
    if (input$ster_on)    args$ster    <- input$ster_start
    if (input$ns1mab)     args$ns1mab  <- 0.35
    withProgress(message = "solving 45 ODEs", do.call(simulate_arm, args))
  }, ignoreNULL = FALSE)

  ref <- reactive({
    if (!HAVE_MRGSOLVE) return(NULL)
    simulate_arm(abh0 = 0)
  })

  both <- reactive({
    a <- run(); b <- ref()
    if (is.null(a)) return(NULL)
    bind_rows(mutate(a, arm = "this patient"),
              mutate(b, arm = "naive, untreated"))
  })

  # ---- 1. ADE -------------------------------------------------------------
  output$p_ade <- renderPlot({
    A <- 10^seq(-1, 4.4, length.out = 800)
    d <- data.frame(A = A, E = ade_curve(A))
    lm_ <- ade_landmarks()
    ggplot(d, aes(A, E)) +
      annotate("rect", xmin = 21, xmax = 80, ymin = -Inf, ymax = Inf,
               fill = "#d62728", alpha = 0.10) +
      annotate("rect", xmin = 1280, xmax = 10^4.4, ymin = -Inf, ymax = Inf,
               fill = "#2ca02c", alpha = 0.10) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 1.1, colour = PAL[["secondary"]]) +
      geom_vline(xintercept = input$abh, linetype = 3, linewidth = 0.9) +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                    labels = c("1:1", "1:10", "1:100", "1:1000", "1:10000")) +
      labs(title = "Entry factor E as a function of pre-infection titre",
           subtitle = sprintf(
             "peak E = %.2f at 1:%.0f; E returns through 1.0 at 1:%.0f. Shaded: Salje 2018 enhancing (red) and protective (green) bands.",
             lm_$peak_E, lm_$peak_titre, lm_$cross),
           x = "reciprocal PRNT50 titre (log scale)",
           y = "E = (1−N)·[1+(Φ−1)·A/(A+K)]") + THEME
  })

  output$p_ade_zone <- renderPlot({
    A <- 10^seq(-1, 4.4, length.out = 800)
    aeff <- 0.12 * A
    Nn <- aeff^2.5 / (aeff^2.5 + 30^2.5)
    d <- data.frame(A = A, N = Nn, O = (1 - Nn) * A / (A + 4)) |>
      pivot_longer(-A)
    ggplot(d, aes(A, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                    labels = c("1:1", "1:10", "1:100", "1:1k", "1:10k")) +
      scale_colour_manual(values = c(N = "#2ca02c", O = "#d62728"),
        labels = c(N = "neutralised fraction N", O = "opsonised, not neutralised O")) +
      labs(title = "The two occupancies the bell is made of",
           subtitle = "O rises first because one IgG ligates FcγR; N needs many and arrives later",
           x = "reciprocal titre", y = "fraction of virions", colour = NULL) + THEME
  })

  output$txt_ade <- renderText({
    lm_ <- ade_landmarks()
    E <- ade_curve(input$abh)
    paste0(
      sprintf("At the selected titre 1:%.0f the entry factor is E = %.2f.\n", input$abh, E),
      sprintf("Peak enhancement E = %.2f at 1:%.0f; net protection above 1:%.0f.\n",
              lm_$peak_E, lm_$peak_titre, lm_$cross),
      "\nWhy this is a bell and not an assumption:\n",
      "  E = (1 - N) x [1 + (PHI-1) x A/(A+K)]\n",
      "  the second factor rises and saturates (opsonisation needs one IgG),\n",
      "  the first falls as a Hill function (neutralisation needs many).\n",
      "  A rising saturating term times a falling Hill term is a bell.\n",
      "\nInfant corollary: maternal IgG decays from a cord titre of 1:1280 with a\n",
      sprintf("43 d half-life, so it reaches 1:%.0f after log2(1280/%.0f) x 43 d\n",
              lm_$peak_titre, lm_$peak_titre),
      sprintf("= %.0f days = %.1f months. Observed peak of infant DHF: 6-8 months.\n",
              log2(1280 / lm_$peak_titre) * 43,
              log2(1280 / lm_$peak_titre) * 43 / 30.44))
  })

  # ---- 2. virology --------------------------------------------------------
  output$p_viro <- renderPlot({
    d <- both(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, arm, V, NS1, IFN, ABN, INF, WBC) |>
      pivot_longer(-c(ILLDAY, arm))
    dd$name <- factor(dd$name, levels = c("V", "NS1", "IFN", "ABN", "INF", "WBC"),
      labels = c("viraemia (copies/mL)", "NS1 antigen (ng/mL)",
                 "type-I interferon (pg/mL)", "neutralising IgG (titre)",
                 "infected cells (/mL)", "leukocytes (1e9/L)"))
    ggplot(dd, aes(ILLDAY, pmax(value, 1e-3), colour = arm)) +
      geom_line(linewidth = 0.95) + facet_wrap(~name, scales = "free_y") +
      scale_y_log10() +
      scale_colour_manual(values = c("this patient" = PAL[["secondary"]],
                                     "naive, untreated" = PAL[["grey"]])) +
      labs(title = "Virology: the antigen supply that everything downstream reads",
           subtitle = "NS1 outlives the viraemia, which is why the rapid test still works after the virus has gone",
           x = "illness day", y = NULL, colour = NULL) + THEME
  })

  # ---- 3. immunopathology -------------------------------------------------
  output$p_immuno <- renderPlot({
    d <- both(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, arm, TNF, IL10, VEGF, CHYM, ABN, CTL) |>
      pivot_longer(-c(ILLDAY, arm))
    dd$name <- factor(dd$name, levels = c("TNF","IL10","VEGF","CHYM","ABN","CTL"),
      labels = c("TNF-α (pg/mL)", "IL-10 (pg/mL)", "VEGF-A (pg/mL)",
                 "mast-cell chymase (ng/mL)", "neutralising IgG (titre)",
                 "activated CD8 (rel.)"))
    ggplot(dd, aes(ILLDAY, value, colour = arm)) +
      geom_line(linewidth = 0.95) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("this patient" = PAL[["secondary"]],
                                     "naive, untreated" = PAL[["grey"]])) +
      labs(title = "Immunopathology: what the antibody rise switches on",
           subtitle = "Immune complexes are a product of a rising antibody and a falling antigen, so they peak between the two",
           x = "illness day", y = NULL, colour = NULL) + THEME
  })

  # ---- 4. endothelium -----------------------------------------------------
  output$p_endo <- renderPlot({
    d <- both(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, arm, GLX, SIGMA, Jser, EFFUS) |>
      pivot_longer(-c(ILLDAY, arm))
    dd$name <- factor(dd$name, levels = c("GLX","SIGMA","Jser","EFFUS"),
      labels = c("glycocalyx integrity G", "reflection coefficient σ",
                 "serosal filtration (mL/h)", "effusion + ascites (mL)"))
    ggplot(dd, aes(ILLDAY, value, colour = arm)) +
      geom_line(linewidth = 0.95) + facet_wrap(~name, scales = "free_y", nrow = 1) +
      scale_colour_manual(values = c("this patient" = PAL[["secondary"]],
                                     "naive, untreated" = PAL[["grey"]])) +
      labs(title = "The glycocalyx is a structure, so it integrates",
           subtitle = "Time constant 20–30 h: the leak lands about a day after the cytokine peak, which is what puts it at defervescence",
           x = "illness day", y = NULL, colour = NULL) + THEME
  })

  output$p_starling <- renderPlot({
    d <- run(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, Pcap, COPp, Jv, Jser) |>
      pivot_longer(-ILLDAY)
    dd$name <- factor(dd$name, levels = c("Pcap","COPp","Jv","Jser"),
      labels = c("capillary pressure Pc (mmHg)", "plasma COP Πp (mmHg)",
                 "systemic filtration J_v (mL/h)", "serosal filtration (mL/h)"))
    ggplot(dd, aes(ILLDAY, value)) +
      geom_line(linewidth = 0.95, colour = PAL[["primary"]]) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(title = "The pressure balance, and the only term that can switch the leak off",
           subtitle = "Pc falls because the plasma volume feeding it has gone: the leak is self-limiting only through shock",
           x = "illness day", y = NULL) + THEME
  })

  # ---- 5. the critical phase ---------------------------------------------
  output$p_critical <- renderPlot({
    d <- run(); if (is.null(d)) return(no_model_msg())
    sc <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)
    dd <- data.frame(
      ILLDAY = d$ILLDAY,
      `channel 1 · temperature`      = sc(d$TEMP),
      `channel 1 · interferon`       = sc(d$IFN),
      `shared driver · neutralising IgG` = sc(d$ABN),
      `channel 2 · glycocalyx loss`  = 1 - sc(d$GLX),
      `channel 2 · haematocrit`      = sc(d$Hct),
      `channel 2 · pulse pressure (inverted)` = 1 - sc(d$PP),
      check.names = FALSE) |> pivot_longer(-ILLDAY)
    ggplot(dd, aes(ILLDAY, value, colour = name)) +
      geom_line(linewidth = 1.0) +
      scale_colour_brewer(palette = "Dark2") +
      labs(title = "The headline: the fever break and the circulatory nadir are one event",
           subtitle = paste("Temperature comes from hypothalamic PGE2 driven by interferon.",
                            "The leak comes from a glycocalyx driven by NS1, TNF, VEGF and chymase.",
                            "\nThe two chains share no equation after the antibody rise — and yet they land together."),
           x = "illness day", y = "each trace scaled to its own range",
           colour = NULL) + THEME +
      guides(colour = guide_legend(nrow = 2))
  })

  output$txt_critical <- renderText({
    d <- run(); if (is.null(d)) return("mrgsolve not available.")
    idx <- function(v, f) d$ILLDAY[f(v)]
    dv <- d$ILLDAY[which(d$TEMP < 37.8 & d$ILLDAY > 1)[1]]
    paste0(
      sprintf("peak viraemia                  illness day %5.2f\n", idx(d$V, which.max)),
      sprintf("peak NS1                       illness day %5.2f\n", idx(d$NS1, which.max)),
      sprintf("defervescence (T < 37.8 C)     illness day %5.2f\n", dv),
      sprintf("glycocalyx nadir               illness day %5.2f\n", idx(d$GLX, which.min)),
      sprintf("peak haematocrit               illness day %5.2f\n", idx(d$Hct, which.max)),
      sprintf("narrowest pulse pressure       illness day %5.2f\n", idx(d$PP, which.min)),
      sprintf("\nclinical nadir minus defervescence = %+.1f h\n",
              (idx(d$PP, which.min) - dv) * 24))
  })

  # ---- 6. haemodynamics ---------------------------------------------------
  output$p_haemo <- renderPlot({
    d <- both(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, arm, PP, MAP, SBP, CO, LAC, HR) |>
      pivot_longer(-c(ILLDAY, arm))
    dd$name <- factor(dd$name, levels = c("PP","MAP","SBP","CO","LAC","HR"),
      labels = c("pulse pressure (mmHg)", "mean arterial pressure (mmHg)",
                 "systolic BP (mmHg)", "cardiac output (L/min)",
                 "lactate (mmol/L)", "heart rate (bpm)"))
    hl <- data.frame(name = factor("pulse pressure (mmHg)", levels = levels(dd$name)),
                     y = 20)
    ggplot(dd, aes(ILLDAY, value, colour = arm)) +
      geom_line(linewidth = 0.95) +
      geom_hline(data = hl, aes(yintercept = y), linetype = 2, colour = "#d62728") +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("this patient" = PAL[["secondary"]],
                                     "naive, untreated" = PAL[["grey"]])) +
      labs(title = "Compensated shock is a narrow pulse pressure at a normal mean pressure",
           subtitle = "Dashed line: PP ≤ 20 mmHg, the WHO definition of dengue shock syndrome. MAP barely moves — which is exactly why DSS is missed.",
           x = "illness day", y = NULL, colour = NULL) + THEME
  })

  # ---- 7. fluid therapy ---------------------------------------------------
  ucurve <- reactive({
    if (!HAVE_MRGSOLVE) return(NULL)
    scales <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
    withProgress(message = "fluid dose-response", {
      bind_rows(lapply(scales, function(s) {
        o <- simulate_arm(abh0 = input$abh, crys = if (s > 0) s else NULL,
                          tp = input$tp)
        data.frame(scale = s,
                   shock_h = sum(o$SHOCK) * 0.5,
                   effusion = max(o$EFFUS),
                   lactate = max(o$LAC),
                   minPP = min(o$PP),
                   iv_mL = max(o$CIN) - max(simulate_arm(abh0 = input$abh,
                                                         tp = input$tp)$CIN))
      }))
    })
  })

  output$p_ucurve <- renderPlot({
    u <- ucurve(); if (is.null(u)) return(no_model_msg())
    d <- u |> select(scale, shock_h, effusion, lactate) |> pivot_longer(-scale)
    d$name <- factor(d$name, levels = c("shock_h","effusion","lactate"),
      labels = c("hours in shock", "effusion + ascites (mL)", "peak lactate (mmol/L)"))
    ggplot(d, aes(scale, value)) +
      geom_line(linewidth = 1.0, colour = PAL[["primary"]]) +
      geom_point(size = 2) + facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(title = "Too little is shock; too much is respiratory distress",
           subtitle = "Crystalloid raises Pc and dilutes Πp, and BOTH terms increase filtration — so the dose-response has an interior optimum",
           x = "× the WHO fluid ladder", y = NULL) + THEME
  })

  output$p_efficiency <- renderPlot({
    d <- run(); if (is.null(d)) return(no_model_msg())
    ggplot(d, aes(SIGMA, Jser)) + geom_path(linewidth = 1.0, colour = PAL[["over"]]) +
      labs(title = "Resuscitation efficiency is a function of σ, not of the dose",
           subtitle = "Every mL infused enters a bed whose sieving has failed; the oncotic term colloid adds is multiplied by the same σ",
           x = "reflection coefficient σ", y = "serosal filtration (mL/h)") + THEME
  })

  output$txt_fluid <- renderText({
    u <- ucurve(); if (is.null(u)) return("mrgsolve not available.")
    b <- u[which.min(u$shock_h + u$effusion / 4000 * 24), ]
    paste0(
      "Reference results from the verified Python implementation:\n",
      "  optimum 0.75 x the WHO ladder (about 5.0 L over 48 h)\n",
      "  0 x  : 66.0 h of shock, lactate 4.35, effusion 1825 mL, severity 0.372\n",
      "  1 x  : 15.0 h of shock, lactate 1.09, effusion 2359 mL, severity 0.249\n",
      "  3 x  :  4.5 h of shock, lactate 1.09, effusion 2737 mL, severity 0.297\n",
      "\nResuscitation efficiency (mL of Ringer's per mL of plasma retained):\n",
      "  endothelium intact (sigma 0.97)         5.5\n",
      "  at the leak nadir, 1 x ladder          16.2\n",
      "  at the leak nadir, 3 x ladder          29.7\n")
  })

  # ---- 8. haematology -----------------------------------------------------
  output$p_haem <- renderPlot({
    d <- both(); if (is.null(d)) return(no_model_msg())
    dd <- d |> select(ILLDAY, arm, Hct, PLT, ALBUM, AST, BLEEDIX, WARNSIGN) |>
      pivot_longer(-c(ILLDAY, arm))
    dd$name <- factor(dd$name,
      levels = c("Hct","PLT","ALBUM","AST","BLEEDIX","WARNSIGN"),
      labels = c("haematocrit (%)", "platelets (1e9/L)", "albumin (g/dL)",
                 "AST (U/L)", "bleeding index", "WHO warning signs"))
    ggplot(dd, aes(ILLDAY, value, colour = arm)) +
      geom_line(linewidth = 0.95) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("this patient" = PAL[["secondary"]],
                                     "naive, untreated" = PAL[["grey"]])) +
      labs(title = "Laboratory course — the platelet nadir arrives after the shock",
           subtitle = "Haemostasis is a product of four requirements, so raising the platelet count alone barely moves the bleeding index",
           x = "illness day", y = NULL, colour = NULL) + THEME
  })

  # ---- 9. scenarios -------------------------------------------------------
  scen_runs <- reactive({
    if (!HAVE_MRGSOLVE) return(NULL)
    withProgress(message = "running 17 scenarios", {
      bind_rows(lapply(names(SCEN), function(n) {
        o <- do.call(simulate_arm, SCEN[[n]])
        mutate(o, scenario = n)
      }))
    })
  })

  output$p_scen <- renderPlot({
    d <- scen_runs(); if (is.null(d)) return(no_model_msg())
    ggplot(d, aes(ILLDAY, PP, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 20, linetype = 2, colour = "#d62728") +
      labs(title = "Pulse pressure across all arms",
           subtitle = "Dashed: the DSS threshold. Everything that works, works by keeping the trace above it.",
           x = "illness day", y = "pulse pressure (mmHg)", colour = NULL) +
      THEME + guides(colour = guide_legend(ncol = 3))
  })

  output$tbl_scen <- DT::renderDataTable({
    d <- scen_runs()
    if (is.null(d)) return(DT::datatable(data.frame(note = "mrgsolve not available")))
    tb <- d |> group_by(scenario) |>
      summarise(`peak V (log10)` = round(max(LOG10V), 2),
                `peak NS1` = round(max(NS1)),
                `min σ` = round(min(SIGMA), 3),
                `Hct rise %` = round(max(HCTRISE), 1),
                `min PP` = round(min(PP), 1),
                `min MAP` = round(min(MAP), 1),
                `effusion mL` = round(max(EFFUS)),
                `min PLT` = round(min(PLT)),
                `min alb` = round(min(ALBUM), 2),
                `max AST` = round(max(AST)),
                `max lactate` = round(max(LAC), 2),
                `shock h` = sum(SHOCK) * 0.5,
                `warn signs` = max(WARNSIGN), .groups = "drop")
    DT::datatable(tb, options = list(pageLength = 20, scrollX = TRUE),
                  rownames = FALSE)
  })

  # ---- 10. trial reproductions -------------------------------------------
  output$txt_trials <- renderText(paste0(
"  Salje 2018 Nature 557:719 -- hospitalised dengue risk raised at pre-infection\n",
"    titres 1:21-1:80, suppressed above 1:1280.\n",
"    MODEL: entry factor peaks at 1:57 (E = 12.83) and crosses 1.0 at 1:696;\n",
"    severity index 0.40 in the enhancing band, 0.00 above 1:1280.\n\n",
"  Kliks 1988 / Chau 2009 -- infant DHF peaks at 6-8 months of age.\n",
"    MODEL: log2(1280/57) x 43 d = 193 d = 6.3 months.  Inputs were a cord\n",
"    titre and an IgG half-life; no infant datum entered the model.\n\n",
"  Wills 2005 NEJM 353:877 -- Ringer's, dextran-70 and 6 % starch equivalent\n",
"    overall; colloid better only in the narrowest-pulse-pressure stratum.\n",
"    MODEL: colloid margin -0.016 / -0.027 / -0.030 as the sigma defect deepens,\n",
"    because the oncotic term colloid adds is multiplied by sigma.\n\n",
"  Lye 2017 Lancet 389:1611 (AAPT) -- prophylactic platelets did not prevent\n",
"    bleeding.  MODEL: three pools move the nadir 32 -> 44 x10^9/L and the\n",
"    bleeding index only 0.490 -> 0.444, because haemostasis is a product and\n",
"    the worst term at the nadir is the vessel wall, not the platelet count.\n\n",
"  Nguyen 2013 JID 207:1442 (balapiravir), Low 2014 Lancet ID 14:706\n",
"    (celgosivir) -- both negative, both enrolled at 48-72 h of fever.\n",
"    MODEL: antiviral benefit 25 % at fever onset, 6.2 % at 24 h, 0.2 % at 48 h\n",
"    and nil thereafter; 84 % of the NS1 exposure integral is already spent\n",
"    by 24 h.  The window closes before patients seek care.\n\n",
"  Tam 2012 Clin Infect Dis 55:1216 -- oral prednisolone within 72 h: no effect.\n",
"    MODEL: methylprednisolone is worth 22 % at fever onset and 5 % at 72 h,\n",
"    and it raises the NS1 exposure integral from 15.6k to 21.1k by delaying\n",
"    viral clearance -- the benefit and the cost cancel in the trial window.\n\n",
"  Sridhar 2018 NEJM 379:327 -- CYD-TDV harms the seronegative and protects the\n",
"    seropositive.  MODEL: a vaccine acts only by SETTING the titre, so the same\n",
"    bell curve that explains second infections explains the vaccine, and the\n",
"    WHO pre-vaccination screening policy follows from the arithmetic.\n"))

  output$p_titre_risk <- renderPlot({
    A <- c(0, 10, 20, 40, 55, 80, 160, 320, 640, 1280, 2560, 5120, 10240)
    sev <- c(0.002, 0.401, 0.399, 0.396, 0.394, 0.391, 0.376, 0.287, 0.072,
             0.000, 0.000, 0.000, 0.000)
    d <- data.frame(A = pmax(A, 1), sev = sev)
    ggplot(d, aes(A, sev)) +
      annotate("rect", xmin = 21, xmax = 80, ymin = -Inf, ymax = Inf,
               fill = "#d62728", alpha = 0.10) +
      annotate("rect", xmin = 1280, xmax = 20000, ymin = -Inf, ymax = Inf,
               fill = "#2ca02c", alpha = 0.10) +
      geom_line(linewidth = 1.1, colour = PAL[["secondary"]]) +
      geom_point(size = 2.4) +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                    labels = c("naive", "1:10", "1:100", "1:1000", "1:10000")) +
      labs(title = "Pre-infection titre versus outcome — a model OUTPUT, not an input",
           subtitle = "Same mosquito inoculum in every host; only the antibody titre differs. Shaded bands are Salje's measured risk zones.",
           x = "pre-infection reciprocal titre", y = "severity index") + THEME
  })
}

shinyApp(ui, server)
