## =============================================================================
##  ros_shiny_app.R — Rosacea QSP Dashboard
##
##  An interactive front-end for ros_mrgsolve_model.R. The whole app is built
##  around one question:
##
##      WHICH of the four effector states does this intervention read,
##      and which endpoint did you happen to measure?
##
##        STATE 1  TONE  (tau ~ 1 h,       no memory)  -> flushing, blanchable
##        STATE 2  VDEN  (tau ~ 3 months,  memory)     -> persistent erythema
##        STATE 3  PAP   (tau ~ 3 weeks,   memory)     -> papules / pustules
##        STATE 4  FIB   (tau ~ years,     HYSTERESIS) -> phyma
##
##  10 tabs:
##    1  Patient & phenotype      — four susceptibility knobs, subtype falls out
##    2  Four states, four clocks — the central plot of the whole app
##    3  Drug PK & target engagement
##    4  Innate amplifier         — KLK5 / LL-37 / TLR2 / MMP-9 and its loops
##    5  Demodex ecology          — density, reservoir, and the relapse clock
##    6  Erythema decomposition   — reversible vs structural, and the rebound
##    7  Clinical endpoints       — CEA / PSA / IGA / ILC / TELSC / PHYGR / DLQI
##    8  Scenario comparison      — the 18 scenarios side by side
##    9  Ocular arm               — MGD, OSDI, lid hygiene, doxycycline
##   10  Virtual population       — subtypes emerging from continuous parameters
##
##  Run with:
##    setwd("rosacea"); shiny::runApp("ros_shiny_app.R")
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("ros_mrgsolve_model.R")

## ---- palette: one colour per STATE, used everywhere ------------------------
PAL <- c(S1 = "#2E6FA7",   # tone / flushing      (blue)
         S2 = "#1B5E8C",   # vessel structure     (deep blue)
         S3 = "#C0392B",   # infiltrate / lesions (red)
         S4 = "#8B6A3F",   # fibrosis / phyma     (tan)
         AMP = "#D9822B",  # innate amplifier     (amber)
         MITE = "#4C8C5A", # Demodex              (green)
         DRUG = "#7A5AA8", # drug exposure        (violet)
         GREY = "#7A7A7A")

theme_ros <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "#555555", size = 11),
          strip.text    = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title  = element_blank())
}

PHENOS <- ros_phenotypes()$phenotype
SCEN   <- ros_scenarios()
SCEN_CHOICES <- setNames(names(SCEN), vapply(SCEN, function(x) x$label, ""))

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well { background: #FBFBFD; }
    h4 { margin-top: 4px; }
    .note { color:#666; font-size: 12px; }
    .big  { font-size: 15px; font-weight: 600; }
  "))),
  titlePanel("Rosacea QSP Dashboard — one amplifier, four states, four clocks"),
  p(class = "note",
    "ETR / PPR / phyma are not settings in this model. Four continuous ",
    "susceptibility parameters plus TIME produce them, and every drug is ",
    "wired to exactly one state — which is what makes the endpoint ",
    "dissociations appear."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("pheno", "Phenotype preset", choices = PHENOS,
                  selected = "PPR-moderate"),
      sliderInput("SPROT", "SPROT — KLK5 / protease set-point",
                  0.5, 5, 3.0, 0.1),
      sliderInput("SNEUR", "SNEUR — neurovascular (TRPV1-CGRP) gain",
                  0.5, 5, 1.8, 0.1),
      sliderInput("SMITE", "SMITE — Demodex carrying capacity",
                  0.5, 25, 14, 0.5),
      sliderInput("SFIBR", "SFIBR — fibrogenic propensity",
                  0.5, 6, 1.5, 0.1),
      sliderInput("ANDROG", "ANDROG — androgen / sebaceous drive",
                  0.8, 2.0, 1.0, 0.05),
      sliderInput("DESENS", "DESENS — alpha-2A desensitisation rate",
                  0.2, 3, 1, 0.1),
      hr(),
      h4("Exposome"),
      sliderInput("TRIGB", "Baseline trigger load", 0, 1, 0.40, 0.05),
      sliderInput("UVLOAD", "Chronic UV load", 0, 1, 0.30, 0.05),
      sliderInput("AVOID", "Trigger avoidance achieved", 0, 1, 0, 0.05),
      sliderInput("STRESSL", "Stress load", 0, 1, 0.20, 0.05),
      checkboxInput("SKINCARE", "Barrier skincare programme", FALSE),
      hr(),
      h4("Treatment"),
      checkboxGroupInput(
        "rx", NULL,
        choices = c("Ivermectin 1% od"          = "ivm",
                    "Metronidazole 0.75% bid"   = "mtz",
                    "Azelaic acid 15% bid"      = "aza",
                    "Minocycline 1.5% foam od"  = "min",
                    "Brimonidine 0.33% gel od"  = "brm",
                    "Oxymetazoline 1% od"       = "oxy",
                    "Doxycycline 40 mg MR"      = "dox40",
                    "Doxycycline 100 mg"        = "dox100",
                    "Isotretinoin 20 mg/day"    = "iso",
                    "PDL x3 (q4w)"              = "laser",
                    "Lid hygiene"               = "lid")),
      sliderInput("weeks", "Treatment duration (weeks)", 2, 52, 16, 1),
      sliderInput("follow", "Follow-up after stopping (weeks)", 0, 52, 12, 1),
      hr(),
      actionButton("go", "Simulate", class = "btn-primary btn-block"),
      p(class = "note",
        "Every run burns in for 10 model years first, so STATE 2 and STATE 4 ",
        "are at their chronic values before treatment starts.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "1 · Patient & phenotype",
          br(),
          fluidRow(
            column(6, h4("Chronic (burn-in) state"), DTOutput("t_pheno")),
            column(6, h4("Where this patient sits"), plotOutput("p_pheno", height = 330))
          ),
          hr(),
          htmlOutput("pheno_txt")
        ),

        tabPanel(
          "2 · Four states, four clocks",
          br(),
          p(class = "big", "The central plot: the same intervention moves these four curves on four different timescales."),
          plotOutput("p_states", height = 520),
          hr(),
          plotOutput("p_states_norm", height = 260)
        ),

        tabPanel(
          "3 · Drug PK & target engagement",
          br(),
          fluidRow(
            column(6, plotOutput("p_pk", height = 320)),
            column(6, plotOutput("p_te", height = 320))
          ),
          hr(),
          h4("Sub-antimicrobial dissociation"),
          p(class = "note",
            "MMP-9 inhibition (IC50 0.15 mg/L) is engaged at 40 mg; the ",
            "antibacterial term (IC50 2 mg/L) is not. Watch the two curves ",
            "separate as the dose changes."),
          plotOutput("p_dox", height = 300)
        ),

        tabPanel(
          "4 · Innate amplifier",
          br(),
          plotOutput("p_amp", height = 420),
          hr(),
          p(class = "note",
            "Loop 1: LL-37 -> TLR2 -> KLK5 -> LL-37.  Loop 2: KLK5 -> IL-1beta ",
            "-> MMP-9 -> KLK5.  Azelaic acid and doxycycline enter at the TOP ",
            "of both loops; metronidazole enters at ROS; ivermectin damps LL-37 ",
            "and TLR2 directly."),
          plotOutput("p_amp_bar", height = 280)
        ),

        tabPanel(
          "5 · Demodex ecology & the relapse clock",
          br(),
          fluidRow(
            column(7, plotOutput("p_demo", height = 340)),
            column(5, plotOutput("p_demo_kill", height = 340))
          ),
          hr(),
          h4("Withdrawal: ivermectin vs metronidazole"),
          p(class = "note",
            "Clearing the reservoir buys TIME TO RELAPSE rather than speed of ",
            "response, because DEMO has a small re-immigration term. Run the ",
            "button below to reproduce the ATTRACT-extension shape."),
          actionButton("go_relapse", "Run relapse experiment"),
          plotOutput("p_relapse", height = 340)
        ),

        tabPanel(
          "6 · Erythema decomposition & rebound",
          br(),
          plotOutput("p_ery", height = 340),
          hr(),
          h4("Brimonidine: two adaptation states, one rebound"),
          p(class = "note",
            "A2AR falls as chronic occupancy internalises the receptor; VDILC ",
            "builds underneath as an unopposed vasodilator drive. Neither is a ",
            "coded side effect — the overshoot after withdrawal is what they do."),
          actionButton("go_rebound", "Run rebound experiment"),
          plotOutput("p_rebound", height = 380)
        ),

        tabPanel(
          "7 · Clinical endpoints",
          br(),
          plotOutput("p_end", height = 480),
          hr(),
          fluidRow(
            column(6, h4("Week-by-week endpoint table"), DTOutput("t_end")),
            column(6, h4("Percent change from baseline"), plotOutput("p_pct", height = 320))
          )
        ),

        tabPanel(
          "8 · Scenario comparison",
          br(),
          selectInput("scen", "Scenarios", choices = SCEN_CHOICES,
                      multiple = TRUE,
                      selected = c("S3", "S4", "S5", "S7", "S14")),
          actionButton("go_scen", "Run selected scenarios"),
          br(), br(),
          plotOutput("p_scen", height = 460),
          hr(),
          DTOutput("t_scen")
        ),

        tabPanel(
          "9 · Ocular arm",
          br(),
          plotOutput("p_ocul", height = 380),
          hr(),
          p(class = "note",
            "The ocular arm shares the amplifier but has its own slow ",
            "meibomian state, which is why ocular severity tracks the skin ",
            "score so poorly. Doxycycline acts through MMP-9; lid hygiene and ",
            "lid IPL act on the gland itself."),
          DTOutput("t_ocul")
        ),

        tabPanel(
          "10 · Virtual population",
          br(),
          sliderInput("npop", "Population size", 40, 400, 120, 20),
          actionButton("go_pop", "Sweep susceptibility space"),
          br(), br(),
          plotOutput("p_pop", height = 420),
          hr(),
          fluidRow(
            column(5, h4("Derived subtypes"), DTOutput("t_pop_tab")),
            column(7, h4("Population"), DTOutput("t_pop"))
          )
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## keep sliders in sync with the phenotype preset -------------------------
  observeEvent(input$pheno, {
    p <- ros_pheno(input$pheno)
    updateSliderInput(session, "SPROT",  value = p$SPROT)
    updateSliderInput(session, "SNEUR",  value = p$SNEUR)
    updateSliderInput(session, "SMITE",  value = p$SMITE)
    updateSliderInput(session, "SFIBR",  value = p$SFIBR)
    updateSliderInput(session, "ANDROG", value = p$ANDROG)
    updateSliderInput(session, "TRIGB",  value = p$TRIGB)
    updateSliderInput(session, "UVLOAD", value = p$UVLOAD)
  })

  user_pars <- reactive({
    list(SPROT = input$SPROT, SNEUR = input$SNEUR, SMITE = input$SMITE,
         SFIBR = input$SFIBR, ANDROG = input$ANDROG, DESENS = input$DESENS,
         TRIGB = input$TRIGB, UVLOAD = input$UVLOAD, AVOID = input$AVOID,
         STRESSL = input$STRESSL,
         SKINCARE = as.numeric(isTRUE(input$SKINCARE)),
         LIDHYG = as.numeric("lid" %in% input$rx))
  })

  user_rx <- reactive({
    W <- input$weeks
    rx <- list()
    if ("ivm"    %in% input$rx) rx <- c(rx, list(dose_ivermectin(weeks = W)))
    if ("mtz"    %in% input$rx) rx <- c(rx, list(dose_metronidazole(weeks = W)))
    if ("aza"    %in% input$rx) rx <- c(rx, list(dose_azelaic(weeks = W)))
    if ("min"    %in% input$rx) rx <- c(rx, list(dose_minocycline_foam(weeks = W)))
    if ("brm"    %in% input$rx) rx <- c(rx, list(dose_brimonidine(weeks = W)))
    if ("oxy"    %in% input$rx) rx <- c(rx, list(dose_oxymetazoline(weeks = W)))
    if ("dox40"  %in% input$rx) rx <- c(rx, list(dose_doxycycline(40, weeks = W)))
    if ("dox100" %in% input$rx) rx <- c(rx, list(dose_doxycycline(100, weeks = W)))
    if ("iso"    %in% input$rx) rx <- c(rx, list(dose_isotretinoin(20, weeks = W)))
    if ("laser"  %in% input$rx) rx <- c(rx, list(dose_laser(3, 28)))
    if (!length(rx)) return(NULL)
    do.call(combine_rx, rx)
  })

  ## the main simulation ----------------------------------------------------
  sim <- eventReactive(input$go, {
    days <- (input$weeks + input$follow) * 7
    ## fine grid if a fast (STATE 1) drug is on board, otherwise daily
    fast <- any(c("brm", "oxy") %in% input$rx)
    ros_run(phenotype = NULL, regimen = user_rx(), days = days,
            delta = if (fast) 0.1 else 1, pars = user_pars())
  }, ignoreNULL = FALSE)

  base_sim <- eventReactive(input$go, {
    days <- (input$weeks + input$follow) * 7
    fast <- any(c("brm", "oxy") %in% input$rx)
    df <- ros_run(phenotype = NULL, regimen = NULL, days = days,
                  delta = if (fast) 0.1 else 1, pars = user_pars())
    df$arm <- "untreated"
    df
  }, ignoreNULL = FALSE)

  both <- reactive({
    a <- sim();  a$arm <- "treated"
    rbind(a, base_sim())
  })

  ## ---- tab 1 -------------------------------------------------------------
  output$t_pheno <- renderDT({
    ss <- ros_steady(ros, user_pars(), days = 3650)
    v <- data.frame(
      quantity = c("Demodex density (mites/cm2)", "CEA (0-4)", "PSA (0-4)",
                   "IGA (0-4)", "Lesion count", "Telangiectasia (0-3)",
                   "Phyma grade (0-3)", "Flush frequency (/day)",
                   "Stinging (0-10)", "DLQI (0-30)", "OSDI (0-100)",
                   "KLK5 (x healthy)", "LL-37 (x healthy)",
                   "IL-17A (x healthy)", "Vessel density (x healthy)",
                   "Tone (0-1)"),
      value = round(c(ss$DEMO, ss$CEA, ss$PSA, ss$IGA, ss$ILC, ss$TELSC,
                      ss$PHYGR, ss$FLFREQ, ss$STING, ss$DLQI, ss$OSDI,
                      ss$KLK, ss$LL37, ss$IL17, ss$VDEN, ss$TONE), 2))
    datatable(v, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$p_pheno <- renderPlot({
    ss <- ros_steady(ros, user_pars(), days = 3650)
    d <- data.frame(
      state = factor(c("STATE 1 tone", "STATE 2 vessels", "STATE 3 lesions",
                       "STATE 4 phyma"),
                     levels = c("STATE 1 tone", "STATE 2 vessels",
                                "STATE 3 lesions", "STATE 4 phyma")),
      value = c((ss$TONE - 0.30) / 0.35, (ss$VDEN - 1) / 1.0,
                ss$ILC / 40, (0.6 * ss$FIB + 0.4 * ss$GLND) / 0.5))
    d$value <- pmin(1, pmax(0, d$value))
    ggplot(d, aes(state, value, fill = state)) +
      geom_col(width = .62) +
      scale_fill_manual(values = unname(PAL[c("S1", "S2", "S3", "S4")])) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(title = "Chronic occupancy of the four effector states",
           subtitle = "the subtype label is a description of this bar chart, not an input",
           x = NULL, y = "fraction of maximal") +
      theme_ros() + theme(legend.position = "none")
  })

  output$pheno_txt <- renderUI({
    ss <- ros_steady(ros, user_pars(), days = 3650)
    sub <- if (ss$ILC >= 10 && ss$CEA >= 2) "mixed ETR + PPR" else
           if (ss$ILC >= 10) "papulopustular (PPR)" else
           if (ss$CEA >= 2) "erythematotelangiectatic (ETR)" else "subclinical"
    HTML(paste0(
      "<p class='big'>Derived subtype: <b>", sub, "</b></p>",
      "<p class='note'>Derived, not set: the model has no subtype switch. ",
      "Raising SMITE moves this patient toward PPR, raising SNEUR toward ETR, ",
      "and SFIBR plus years toward phyma.</p>"))
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$p_states <- renderPlot({
    d <- both() %>%
      select(time, arm, TONE, VDEN, ILC, PHYGR) %>%
      pivot_longer(c(TONE, VDEN, ILC, PHYGR))
    d$name <- recode(d$name,
                     TONE = "STATE 1 — tone (tau ~ 1 h)",
                     VDEN = "STATE 2 — vessel density (tau ~ 3 months)",
                     ILC  = "STATE 3 — lesion count (tau ~ 3 weeks)",
                     PHYGR = "STATE 4 — phyma grade (tau ~ years)")
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = .9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["S3"]])) +
      labs(title = "Four states, four clocks",
           subtitle = "same simulation, four different time constants — this is why endpoints dissociate",
           x = "days", y = NULL) +
      theme_ros()
  })

  output$p_states_norm <- renderPlot({
    d <- sim()
    n <- data.frame(time = d$time,
                    `STATE 1` = (d$TONE - 0.30) / 0.35,
                    `STATE 2` = (d$VDEN - 1) / 1.0,
                    `STATE 3` = d$ILC / max(1, max(d$ILC)),
                    `STATE 4` = (0.6 * d$FIB + 0.4 * d$GLND) / 0.5,
                    check.names = FALSE) %>%
      pivot_longer(-time)
    ggplot(n, aes(time, pmin(1, pmax(0, value)), colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = unname(PAL[c("S1", "S2", "S3", "S4")])) +
      labs(title = "The same four states on one normalised axis",
           x = "days", y = "fraction of maximal") +
      theme_ros()
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>%
      select(time, IVMFO, MTZSK, AZASK, BRMEF, OXYEF, MINSK, CDOXO, CISOO) %>%
      pivot_longer(-time) %>% filter(value > 1e-8)
    if (!nrow(d)) return(
      ggplot() + labs(title = "No drug on board") + theme_ros())
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = .8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Drug exposure", subtitle = "topical in applications, systemic in mg/L",
           x = "days", y = NULL) +
      theme_ros() + theme(legend.position = "none")
  })

  output$p_te <- renderPlot({
    d <- sim() %>%
      select(time, `alpha-2A occupancy` = OCC_A2,
             `realised alpha-2A constriction` = TE_A2,
             `MMP-9 inhibition` = TE_MMP,
             `antibacterial engagement` = TE_ABX) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = .9) +
      labs(title = "Target engagement",
           subtitle = "occupancy is not effect: realised constriction falls as A2AR internalises",
           x = "days", y = "fraction") +
      theme_ros()
  })

  output$p_dox <- renderPlot({
    d <- ros_doxy_dissociation(weeks = min(input$weeks, 16))
    s <- d %>% group_by(dose_mg) %>%
      summarise(ILC = last(ILC), DEMO = last(DEMO),
                MMP9_inh = last(TE_MMP), ABX = last(TE_ABX), .groups = "drop") %>%
      pivot_longer(-dose_mg)
    ggplot(s, aes(factor(dose_mg), value, fill = name)) +
      geom_col(position = "dodge") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Doxycycline: anti-protease, not anti-parasitic",
           subtitle = "lesion count falls with dose; Demodex density does not move",
           x = "daily dose (mg)", y = NULL) +
      theme_ros() + theme(legend.position = "none")
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$p_amp <- renderPlot({
    d <- both() %>%
      select(time, arm, KLK, LL37, TLR2, IL1B, MMP9, ROS, MC, IL17) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = .85) +
      geom_hline(yintercept = 1, linetype = 3, colour = PAL[["GREY"]]) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["AMP"]])) +
      labs(title = "The innate amplifier and its two positive loops",
           subtitle = "1.0 = healthy reference",
           x = "days", y = "x healthy") +
      theme_ros()
  })

  output$p_amp_bar <- renderPlot({
    a <- sim(); b <- base_sim()
    nm <- c("KLK", "LL37", "TLR2", "IL1B", "MMP9", "ROS", "MC", "NEU", "IL17")
    d <- data.frame(node = nm,
                    treated = as.numeric(a[nrow(a), nm]),
                    untreated = as.numeric(b[nrow(b), nm])) %>%
      pivot_longer(-node)
    ggplot(d, aes(node, value, fill = name)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 1, linetype = 3) +
      scale_fill_manual(values = c(untreated = PAL[["GREY"]],
                                   treated = PAL[["AMP"]])) +
      labs(title = "End-of-treatment amplifier state", x = NULL, y = "x healthy") +
      theme_ros()
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$p_demo <- renderPlot({
    d <- both()
    ggplot(d, aes(time, DEMO, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.8, linetype = 3) +
      annotate("text", x = 0, y = 0.95, hjust = 0, size = 3.2,
               colour = PAL[["GREY"]], label = "healthy density ~0.7-0.8/cm2") +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["MITE"]])) +
      labs(title = "Demodex density", x = "days", y = "mites / cm2") +
      theme_ros()
  })

  output$p_demo_kill <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = TE_KILL, colour = "kill rate (1/day)"), linewidth = .9) +
      geom_line(aes(y = BOL, colour = "Bacillus antigen load"), linewidth = .9) +
      scale_colour_manual(values = c(`kill rate (1/day)` = PAL[["MITE"]],
                                     `Bacillus antigen load` = PAL[["S3"]])) +
      labs(title = "Mite kill and the antigen burst that follows",
           subtitle = "killing mites transiently RAISES antigen exposure (ABURST)",
           x = "days", y = NULL) +
      theme_ros()
  })

  relapse <- eventReactive(input$go_relapse, {
    ros_relapse(weeks = input$weeks, follow = max(60, input$follow * 7))
  })

  output$p_relapse <- renderPlot({
    d <- relapse()
    ggplot(d, aes(time, ILC, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = input$weeks * 7, linetype = 2) +
      annotate("text", x = input$weeks * 7, y = Inf, vjust = 1.4, hjust = -0.05,
               size = 3.3, label = "treatment stops") +
      scale_colour_manual(values = c(ivermectin = PAL[["MITE"]],
                                     metronidazole = PAL[["DRUG"]])) +
      labs(title = "On-treatment difference is small; the relapse difference is not",
           subtitle = "the mite reservoir sets the relapse clock (IMMIG)",
           x = "days", y = "inflammatory lesion count") +
      theme_ros()
  })

  ## ---- tab 6 -------------------------------------------------------------
  output$p_ery <- renderPlot({
    d <- both() %>%
      select(time, arm, `reversible (tone)` = ERYS1,
             `structural (vessels)` = ERYS2, CEA) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = .9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["S2"]])) +
      labs(title = "Erythema is two things measured as one number",
           subtitle = "an alpha-agonist can only move the left panel; a laser only the middle",
           x = "days", y = NULL) +
      theme_ros()
  })

  rebound <- eventReactive(input$go_rebound, {
    ros_rebound(weeks = min(input$weeks, 12),
                follow = max(28, input$follow * 7))
  })

  output$p_rebound <- renderPlot({
    d <- rebound()
    ggplot(d, aes(time, CEA, colour = arm)) +
      geom_line(linewidth = .8) +
      geom_vline(xintercept = min(input$weeks, 12) * 7, linetype = 2) +
      labs(title = "Brimonidine: the day-1 win and the week-9 overshoot",
           subtitle = "dashed line = withdrawal; higher DESENS gives a larger rebound",
           x = "days", y = "CEA (0-4)") +
      theme_ros()
  })

  ## ---- tab 7 -------------------------------------------------------------
  output$p_end <- renderPlot({
    d <- both() %>%
      select(time, arm, CEA, PSA, IGA, ILC, TELSC, PHYGR, FLFREQ, DLQI) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = .9) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["S3"]])) +
      labs(title = "Clinical endpoints", x = "days", y = NULL) +
      theme_ros()
  })

  output$t_end <- renderDT({
    d <- sim()
    w <- d[d$time %% 7 == 0, ]
    w <- data.frame(week = w$time / 7,
                    CEA = round(w$CEA, 2), PSA = round(w$PSA, 2),
                    IGA = round(w$IGA, 2), ILC = round(w$ILC, 1),
                    TELSC = round(w$TELSC, 2), PHYGR = round(w$PHYGR, 2),
                    flush = round(w$FLFREQ, 2), DLQI = round(w$DLQI, 1))
    datatable(w, rownames = FALSE, options = list(pageLength = 12))
  })

  output$p_pct <- renderPlot({
    d <- both()
    p <- do.call(rbind, lapply(split(d, d$arm), function(s) {
      data.frame(arm = s$arm[1], time = s$time,
                 ILC = 100 * (s$ILC - s$ILC[1]) / max(s$ILC[1], 1e-6),
                 CEA = 100 * (s$CEA - s$CEA[1]) / max(s$CEA[1], 1e-6))
    })) %>% pivot_longer(c(ILC, CEA))
    ggplot(p, aes(time, value, colour = arm, linetype = name)) +
      geom_line(linewidth = .9) +
      geom_hline(yintercept = 0, linetype = 3) +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = PAL[["S3"]])) +
      labs(title = "Percent change from baseline",
           subtitle = "the way phase-3 papers report lesion counts",
           x = "days", y = "% change") +
      theme_ros()
  })

  ## ---- tab 8 -------------------------------------------------------------
  scen <- eventReactive(input$go_scen, {
    req(input$scen)
    ros_run_all(keys = input$scen, weeks = input$weeks, delta = 1)
  })

  output$p_scen <- renderPlot({
    d <- scen() %>%
      select(time, key, CEA, ILC, DEMO, TELSC, FLFREQ, DLQI) %>%
      pivot_longer(c(-time, -key))
    ggplot(d, aes(time, value, colour = key)) +
      geom_line(linewidth = .9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Scenario comparison", x = "days", y = NULL) +
      theme_ros()
  })

  output$t_scen <- renderDT({
    s <- ros_summary(scen())
    num <- vapply(s, is.numeric, TRUE)
    s[num] <- lapply(s[num], round, 2)
    datatable(s, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })

  ## ---- tab 9 -------------------------------------------------------------
  output$p_ocul <- renderPlot({
    d <- both() %>%
      select(time, arm, `meibomian dysfunction` = MGDX,
             `ocular inflammation` = OCUL, OSDI) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = .9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(untreated = PAL[["GREY"]],
                                     treated = "#1E8C82")) +
      labs(title = "Ocular arm", x = "days", y = NULL) +
      theme_ros()
  })

  output$t_ocul <- renderDT({
    d <- sim(); b <- base_sim()
    v <- data.frame(
      quantity = c("MGD (0-1)", "Ocular inflammation (0-1)", "OSDI (0-100)",
                   "CEA (0-4)", "ILC"),
      untreated = round(as.numeric(b[nrow(b), c("MGDX", "OCUL", "OSDI",
                                               "CEA", "ILC")]), 2),
      treated = round(as.numeric(d[nrow(d), c("MGDX", "OCUL", "OSDI",
                                              "CEA", "ILC")]), 2))
    datatable(v, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 10 ------------------------------------------------------------
  pop <- eventReactive(input$go_pop, {
    withProgress(message = "sweeping susceptibility space", {
      ros_vpop(n = input$npop)
    })
  })

  output$p_pop <- renderPlot({
    d <- pop()
    ggplot(d, aes(CEA, ILC, colour = subtype, size = SMITE)) +
      geom_point(alpha = .75) +
      scale_size_continuous(range = c(1.2, 5)) +
      labs(title = "Subtypes are an emergent partition of a continuous space",
           subtitle = "colour is derived from the simulation output, never assigned as an input",
           x = "CEA (0-4)", y = "inflammatory lesion count") +
      theme_ros()
  })

  output$t_pop_tab <- renderDT({
    d <- as.data.frame(table(pop()$subtype))
    names(d) <- c("derived subtype", "n")
    datatable(d, rownames = FALSE, options = list(dom = "t"))
  })

  output$t_pop <- renderDT({
    datatable(pop(), rownames = FALSE,
              options = list(pageLength = 8, scrollX = TRUE)) %>%
      formatRound(columns = c("SPROT", "SNEUR", "SMITE", "SFIBR",
                              "TRIGB", "UVLOAD"), digits = 2)
  })
}

shinyApp(ui, server)
