## =============================================================================
##  cio_shiny_app.R
##  Cisplatin-Induced Ototoxicity (CIO) — interactive QSP dashboard
##  시스플라틴 유발 이독성 QSP 대시보드
##
##  10 tabs:
##    1  Patient & regimen      — who is being treated and with what
##    2  Pharmacokinetics       — the cascade and the lag that creates the window
##    3  The 6-hour rule        — what is still ahead of you at time t
##    4  Audiogram              — threshold shift across the eight bands
##    5  Tonotopy               — reserve vs uptake, and the critical load PT*
##    6  Cochlear cell state    — OHC / IHC / ganglion / glutathione / EP
##    7  Grading scales         — Brock, SIOP Boston, ASHA, CTCAE as OUTPUTS
##    8  Otoprotection          — delay sweep, route comparison, efficacy cost
##    9  Kidney & tumour        — the amplifier loop and the efficacy trade-off
##   10  Scenario comparison    — the 20-scenario library side by side
##
##  Run with:  shiny::runApp("cio_shiny_app.R")
##  Requires:  shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##  Loads the model from cio_mrgsolve_model.R in the same directory.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## -- load the model and the dosing helpers ------------------------------------
source("cio_mrgsolve_model.R", local = TRUE)

FREQ_HZ    <- c(250, 500, 1000, 2000, 4000, 6000, 8000, 12500)
FREQ_LAB   <- c("0.25k", "0.5k", "1k", "2k", "4k", "6k", "8k", "12.5k")
TS_COLS    <- paste0("TS", 1:8)
XPOS       <- log10(FREQ_HZ/165.4 + 0.88)/2.1

## Audiograms are conventionally drawn with threshold increasing downwards.
audiogram_scale <- scale_y_reverse(limits = c(110, -5), breaks = seq(0, 110, 10))

theme_cio <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Cisplatin-Induced Ototoxicity — QSP Dashboard (시스플라틴 유발 이독성)"),
  tags$p(style = "color:#555;margin-top:-8px",
         "73-ODE model - 8 tonotopic bands at their Greenwood positions - ",
         "grading scales are computed as OUTPUTS, not inputs."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("age",  "Age (years)", 0.5, 40, 8, step = 0.5),
      sliderInput("bsa",  "BSA (m²)", 0.4, 2.2, 1.73, step = 0.01),
      sliderInput("gfr0", "Baseline GFR (mL/min)", 20, 130, 100, step = 5),

      h4("요법 (Regimen)"),
      selectInput("drug", "Platinum", c("Cisplatin" = "cis", "Carboplatin" = "carb")),
      sliderInput("dose",   "Dose per administration (mg/m²)", 10, 1600, 100, step = 10),
      sliderInput("cycles", "Number of cycles", 1, 12, 6, step = 1),
      sliderInput("cyclen", "Cycle length (days)", 7, 28, 21, step = 7),
      sliderInput("ndays",  "Consecutive days per cycle", 1, 5, 1, step = 1),
      sliderInput("follow", "Off-treatment follow-up (months)", 0, 48, 0, step = 3),

      h4("이독성 보호 (Otoprotection)"),
      radioButtons("prot", NULL,
                   c("None" = "none",
                     "Systemic sodium thiosulfate" = "sts",
                     "Intratympanic thiosulfate" = "it",
                     "IV N-acetylcysteine" = "nac"), selected = "none"),
      conditionalPanel("input.prot == 'sts'",
        sliderInput("sts_dose",  "Thiosulfate dose (g/m²)", 4, 30, 20, step = 1),
        sliderInput("sts_delay", "Delay after infusion (h)", 0, 24, 6, step = 0.5)),
      conditionalPanel("input.prot == 'it'",
        sliderInput("it_umol",  "Intratympanic dose (µmol)", 20, 800, 161, step = 10),
        sliderInput("it_delay", "Delay after infusion (h)", 0, 24, 6, step = 0.5)),
      conditionalPanel("input.prot == 'nac'",
        sliderInput("nac_mg", "NAC dose (mg)", 500, 15000, 2250, step = 250)),

      h4("동반 노출 (Co-exposures)"),
      checkboxInput("amgly", "Concurrent aminoglycoside", FALSE),
      checkboxInput("noise", "Concurrent noise ≥85 dBA", FALSE),
      checkboxInput("furo",  "Loop diuretic (furosemide)", FALSE),

      h4("구조 검정 (Structural switches)"),
      checkboxInput("nogsh",   "Remove the reserve gradient (BGSH = 0)", FALSE),
      checkboxInput("noupt",   "Remove the uptake gradient (BUPT = 0)", FALSE),
      checkboxInput("noneph",  "Clamp GFR (disable the amplifier loop)", FALSE),
      sliderInput("tlab", "Labile platinum half-life (days)", 5, 365, 60, step = 5),
      sliderInput("tret", "Bound platinum half-life (days)", 30, 7300, 730, step = 30),
      helpText("The bound pool carries nearly all the retained mass but none of ",
               "the damage: moving its half-life should not move the audiogram.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. 환자·요법 (Profile)",
          h4("Regimen summary"), tableOutput("t_profile"),
          h4("Cumulative platinum exposure"), plotOutput("p_cumdose", height = 260),
          h4("End-of-treatment summary"), tableOutput("t_summary")),

        tabPanel("2. 약동학 (Pharmacokinetics)",
          h4("Cisplatin cascade over one cycle: plasma → stria → perilymph"),
          plotOutput("p_pk", height = 330),
          h4("Same data on a log scale — the lag is four orders of magnitude wide"),
          plotOutput("p_pklog", height = 300),
          tableOutput("t_pk")),

        tabPanel("3. 6시간 규칙 (The 6-hour rule)",
          h4("Fraction of each exposure integral still AHEAD of you at time t"),
          plotOutput("p_ahead", height = 350),
          tableOutput("t_ahead"),
          tags$p(style = "color:#444",
                 "The therapeutic window of delayed thiosulfate is not a ",
                 "pharmacological trick. It is the difference between a 22-minute ",
                 "plasma half-life and a perilymph peak that arrives at 7 hours.")),

        tabPanel("4. 청력도 (Audiogram)",
          h4("Threshold shift across the eight bands"),
          plotOutput("p_audio", height = 400),
          h4("Evolution over the treatment course"),
          plotOutput("p_audio_time", height = 320)),

        tabPanel("5. 음조지형 (Tonotopy)",
          h4("Reserve, uptake and the critical labile load PT*"),
          plotOutput("p_grad", height = 330),
          h4("Retained platinum against its own per-band threshold"),
          plotOutput("p_ptstar", height = 330),
          tableOutput("t_grad")),

        tabPanel("6. 와우 세포 상태 (Cochlear state)",
          fluidRow(
            column(6, h4("Outer hair cell survival"), plotOutput("p_ohc", height = 290)),
            column(6, h4("Glutathione reserve"), plotOutput("p_gsh", height = 290))),
          fluidRow(
            column(6, h4("Inner hair cells and ganglion"), plotOutput("p_ihc", height = 290)),
            column(6, h4("Endocochlear potential and inflammation"),
                   plotOutput("p_ep", height = 290)))),

        tabPanel("7. 등급 체계 (Grading scales)",
          h4("Grades are a staircase quantising a continuous variable"),
          plotOutput("p_grades", height = 330),
          h4("Grade against cumulative dose"),
          plotOutput("p_stair", height = 320),
          DTOutput("t_grades")),

        tabPanel("8. 이독성 보호 (Otoprotection)",
          h4("Thiosulfate delay sweep: otoprotection against efficacy cost"),
          plotOutput("p_delay", height = 350),
          h4("Route comparison — systemic is flat, intratympanic is base-weighted"),
          plotOutput("p_route", height = 320),
          DTOutput("t_delay")),

        tabPanel("9. 신장·종양 (Kidney & tumour)",
          fluidRow(
            column(6, h4("GFR and tubular injury"), plotOutput("p_gfr", height = 300)),
            column(6, h4("Tumour adducts and cumulative log-kill"),
                   plotOutput("p_tum", height = 300))),
          h4("The amplifier loop: falling GFR raises platinum exposure"),
          tableOutput("t_amp")),

        tabPanel("10. 시나리오 비교 (Scenarios)",
          h4("The 20-scenario library"),
          DTOutput("t_scen"),
          h4("Speech-frequency average across scenarios"),
          plotOutput("p_scen", height = 400))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## -- assemble the parameter set from the sidebar ---------------------------
  pset <- reactive({
    p <- list(AGE = input$age, BSA = input$bsa, GFR0 = input$gfr0,
              AMGLY = as.numeric(input$amgly), NOISE = as.numeric(input$noise),
              FURO = as.numeric(input$furo),
              TLAB_D = input$tlab, TRET_D = input$tret)
    if (input$nogsh)  p$BGSH  <- 0
    if (input$noupt)  p$BUPT  <- 0
    if (input$noneph) { p$KLOSS <- 0; p$KPERM <- 0 }
    if (input$drug == "carb") {
      p$FADC <- 0.25; p$KINBL <- 0.085*0.060; p$KKID <- 0.55*0.25
    }
    p
  })

  mod_r <- reactive(param(mod, pset()))

  ev_r <- reactive({
    cio_ev(cycles = input$cycles, dose_mgm2 = input$dose,
           cycle_h = input$cyclen*24, ndays = input$ndays,
           sts_cycles = if (input$prot == "sts") seq_len(input$cycles) - 1L else NULL,
           sts_gm2  = if (input$prot == "sts") input$sts_dose else 20,
           sts_delay = if (input$prot == "sts") input$sts_delay
                       else if (input$prot == "it") input$it_delay else 6,
           it_cycles = if (input$prot == "it") seq_len(input$cycles) - 1L else NULL,
           it_umol   = input$it_umol,
           nac_cycles = if (input$prot == "nac") seq_len(input$cycles) - 1L else NULL,
           nac_mg = input$nac_mg, bsa = input$bsa,
           carbo = input$drug == "carb")
  })

  sim <- reactive({
    as.data.frame(cio_run(mod_r(), ev_r(), cycles = input$cycles,
                          cycle_h = input$cyclen*24,
                          follow_h = input$follow*30.4*24, delta = 6))
  })

  last <- reactive({ s <- sim(); s[nrow(s), ] })

  ## -- one-cycle fine-grained run, for the PK and 6-hour-rule tabs -----------
  sim1 <- reactive({
    as.data.frame(mrgsim(mod_r(),
                         events = cio_ev(1, input$dose, bsa = input$bsa,
                                         carbo = input$drug == "carb"),
                         end = 72, delta = 0.05, atol = 1e-10, rtol = 1e-8))
  })

  ## ---- tab 1 --------------------------------------------------------------
  output$t_profile <- renderTable({
    data.frame(
      Field = c("Age (y)", "BSA (m²)", "Baseline GFR", "Platinum", "Dose",
                "Cycles", "Cycle length", "Cumulative dose", "Otoprotection",
                "Co-exposures"),
      Value = c(sprintf("%.1f", input$age), sprintf("%.2f", input$bsa),
                sprintf("%d mL/min", input$gfr0),
                if (input$drug == "cis") "Cisplatin" else "Carboplatin",
                sprintf("%d mg/m² × %d day(s)", input$dose, input$ndays),
                as.character(input$cycles), sprintf("%d days", input$cyclen),
                sprintf("%d mg/m²", input$dose*input$cycles*input$ndays),
                switch(input$prot, none = "none",
                       sts = sprintf("STS %g g/m² at %.1f h", input$sts_dose,
                                     input$sts_delay),
                       it  = sprintf("Intratympanic %g µmol", input$it_umol),
                       nac = sprintf("NAC %g mg", input$nac_mg)),
                paste(c(if (input$amgly) "aminoglycoside",
                        if (input$noise) "noise",
                        if (input$furo) "furosemide"), collapse = ", ")))
  })

  output$p_cumdose <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/24, PTAUC)) + geom_line(linewidth = 1, colour = "#6C63C7") +
      labs(x = "Days", y = "Cumulative labile cochlear platinum (AUC)",
           title = "The damaging variable is an integral, not a concentration") +
      theme_cio
  })

  output$t_summary <- renderTable({
    e <- last()
    data.frame(Metric = c("8 kHz shift (dB)", "4 kHz shift (dB)",
                          "Speech-frequency PTA (dB)", "Brock grade",
                          "SIOP Boston grade", "CTCAE grade",
                          "Speech-in-noise penalty (dB SNR)",
                          "Endocochlear potential (mV)", "GFR (mL/min)",
                          "Tumour log-kill"),
               Value = c(sprintf("%.1f", e$TS7), sprintf("%.1f", e$TS5),
                         sprintf("%.2f", e$PTA), as.character(e$BROCK),
                         as.character(e$SIOP), as.character(e$CTCAE),
                         sprintf("%.1f", e$SNRLOSS), sprintf("%.1f", e$EP),
                         sprintf("%.1f", e$GFROUT), sprintf("%.2f", e$TUMLK)))
  })

  ## ---- tab 2 --------------------------------------------------------------
  pkdat <- reactive({
    s <- sim1()
    data.frame(time = rep(s$time, 3),
               conc = c(s$CISC, s$STRPT, s$PERI),
               pool = rep(c("free plasma platinum", "stria vascularis",
                            "perilymph"), each = nrow(s)))
  })

  output$p_pk <- renderPlot({
    ggplot(pkdat(), aes(time, conc, colour = pool)) +
      geom_line(linewidth = 1) + xlim(0, 48) +
      scale_colour_manual(values = c("free plasma platinum" = "#4A90D9",
                                     "stria vascularis" = "#8E44AD",
                                     "perilymph" = "#C0392B")) +
      labs(x = "Hours after the start of infusion", y = "Concentration (µM)",
           colour = NULL,
           title = "Free plasma platinum is gone before perilymph platinum peaks") +
      theme_cio
  })

  output$p_pklog <- renderPlot({
    d <- pkdat(); d <- d[d$conc > 1e-8, ]
    ggplot(d, aes(time, conc, colour = pool)) + geom_line(linewidth = 1) +
      scale_y_log10() + xlim(0, 48) +
      geom_vline(xintercept = 6, linetype = "dashed", colour = "#B8860B") +
      annotate("text", x = 6.4, y = 1e-4, hjust = 0, colour = "#B8860B",
               label = "systemic thiosulfate given here") +
      labs(x = "Hours", y = "Concentration (µM, log scale)", colour = NULL) +
      theme_cio
  })

  output$t_pk <- renderTable({
    s <- sim1()
    data.frame(Pool = c("free plasma platinum", "stria vascularis", "perilymph"),
               Cmax = sprintf("%.4g", c(max(s$CISC), max(s$STRPT), max(s$PERI))),
               Tmax_h = sprintf("%.2f", c(s$time[which.max(s$CISC)],
                                          s$time[which.max(s$STRPT)],
                                          s$time[which.max(s$PERI)])))
  })

  ## ---- tab 3 --------------------------------------------------------------
  ahead <- reactive({
    s <- sim1()
    f <- function(v) {
      cum <- rev(cumsum(rev(v)))          # trailing integral, uniform grid
      cum/cum[1]*100
    }
    data.frame(time = s$time,
               plasma = f(s$CISC),
               cochlea = f(s$PERI),
               adduct = f(s$TUMPT))
  })

  output$p_ahead <- renderPlot({
    d <- ahead() %>% pivot_longer(-time)
    d$name <- factor(d$name, c("plasma", "adduct", "cochlea"),
                     c("free plasma platinum", "tumour adduct formation",
                       "cochlear platinum uptake"))
    ggplot(d[d$time <= 36, ], aes(time, value, colour = name)) +
      geom_line(linewidth = 1) + scale_y_log10() +
      geom_vline(xintercept = 6, linetype = "dashed", colour = "#B8860B") +
      scale_colour_manual(values = c("#4A90D9", "#4CAF50", "#C0392B")) +
      labs(x = "Hours after the start of infusion",
           y = "% of the integral still ahead (log scale)", colour = NULL,
           title = "At 6 h the cochlea has barely started; the tumour has finished") +
      theme_cio
  })

  output$t_ahead <- renderTable({
    d <- ahead()
    idx <- sapply(c(0, 1, 2, 4, 6, 8, 12, 24), function(t) which.min(abs(d$time - t)))
    data.frame(`t (h)` = sprintf("%.0f", d$time[idx]),
               `plasma ahead %` = sprintf("%.4f", d$plasma[idx]),
               `tumour ahead %` = sprintf("%.2f", d$adduct[idx]),
               `cochlea ahead %` = sprintf("%.1f", d$cochlea[idx]),
               check.names = FALSE)
  })

  ## ---- tab 4 --------------------------------------------------------------
  ts_last <- reactive({
    e <- last()
    data.frame(freq = factor(FREQ_LAB, FREQ_LAB),
               khz = FREQ_HZ/1000,
               shift = as.numeric(e[TS_COLS]))
  })

  output$p_audio <- renderPlot({
    ggplot(ts_last(), aes(khz, shift)) +
      geom_line(linewidth = 1.1, colour = "#C0392B") +
      geom_point(size = 3, colour = "#C0392B") +
      geom_hline(yintercept = 20, linetype = "dotted", colour = "#B8860B") +
      geom_hline(yintercept = 40, linetype = "dashed", colour = "#8E44AD") +
      annotate("text", x = 0.26, y = 22, hjust = 0, size = 3.4,
               colour = "#B8860B", label = "SIOP Boston criterion (20 dB)") +
      annotate("text", x = 0.26, y = 42, hjust = 0, size = 3.4,
               colour = "#8E44AD", label = "Brock criterion (40 dB)") +
      scale_x_log10(breaks = FREQ_HZ/1000, labels = FREQ_LAB) +
      audiogram_scale +
      labs(x = "Frequency (kHz)", y = "Threshold shift (dB HL)",
           title = "End-of-treatment audiogram") +
      theme_cio
  })

  output$p_audio_time <- renderPlot({
    s <- sim()
    d <- s[, c("time", TS_COLS)] %>% pivot_longer(-time)
    d$freq <- factor(FREQ_LAB[as.integer(sub("TS", "", d$name))], FREQ_LAB)
    ggplot(d, aes(time/24, value, colour = freq)) + geom_line(linewidth = 0.9) +
      scale_colour_viridis_d(option = "plasma", end = 0.9) +
      labs(x = "Days", y = "Threshold shift (dB HL)", colour = "Frequency",
           title = "The loss marches apex-ward as cumulative dose rises") +
      theme_cio
  })

  ## ---- tab 5 --------------------------------------------------------------
  grads <- reactive({
    p <- pset()
    bgsh <- if (!is.null(p$BGSH)) p$BGSH else 1.35
    bupt <- if (!is.null(p$BUPT)) p$BUPT else 0.95
    agef <- (0.62 + 0.38*input$age/(4 + input$age))/(0.62 + 0.38*25/(4 + 25))
    gmax <- exp(-bgsh*XPOS)*agef
    upt  <- exp(bupt*(XPOS - 0.5))
    data.frame(freq = factor(FREQ_LAB, FREQ_LAB), khz = FREQ_HZ/1000,
               reserve = gmax, uptake = upt,
               PTstar = 0.045*gmax/(0.02720*1.0),
               vulnerability = (gmax/upt)/(gmax[7]/upt[7]))
  })

  output$p_grad <- renderPlot({
    d <- grads() %>%
      select(khz, reserve, uptake, vulnerability) %>% pivot_longer(-khz)
    ggplot(d, aes(khz, value, colour = name)) + geom_line(linewidth = 1) +
      geom_point(size = 2.4) + scale_x_log10(breaks = FREQ_HZ/1000, labels = FREQ_LAB) +
      labs(x = "Frequency (kHz)", y = "Relative value", colour = NULL,
           title = "Vulnerability is the reserve gradient divided by the uptake gradient") +
      theme_cio
  })

  output$p_ptstar <- renderPlot({
    e <- last()
    pt <- as.numeric(e[paste0("PT", 1:8)])
    d <- data.frame(khz = rep(FREQ_HZ/1000, 2),
                    value = c(pt, grads()$PTstar),
                    what = rep(c("retained labile platinum",
                                 "critical load PT*"), each = 8))
    ggplot(d, aes(khz, value, colour = what)) + geom_line(linewidth = 1) +
      geom_point(size = 2.4) +
      scale_x_log10(breaks = FREQ_HZ/1000, labels = FREQ_LAB) +
      scale_colour_manual(values = c("critical load PT*" = "#B8860B",
                                     "retained labile platinum" = "#C0392B")) +
      labs(x = "Frequency (kHz)", y = "Labile platinum (normalised)", colour = NULL,
           title = "Bands whose platinum exceeds their own threshold lose their reserve") +
      theme_cio
  })

  output$t_grad <- renderTable({
    e <- last(); g <- grads()
    data.frame(Frequency = FREQ_LAB,
               `Greenwood x` = sprintf("%.3f", XPOS),
               Reserve = sprintf("%.3f", g$reserve),
               Uptake = sprintf("%.3f", g$uptake),
               `PT*` = sprintf("%.3f", g$PTstar),
               `retained PT` = sprintf("%.3f", as.numeric(e[paste0("PT", 1:8)])),
               `over threshold` = ifelse(as.numeric(e[paste0("PT", 1:8)]) > g$PTstar,
                                         "yes", "no"),
               check.names = FALSE)
  })

  ## ---- tab 6 --------------------------------------------------------------
  band_plot <- function(prefix, ylab, ttl) {
    s <- sim()
    d <- s[, c("time", paste0(prefix, 1:8))] %>% pivot_longer(-time)
    d$freq <- factor(FREQ_LAB[as.integer(sub(prefix, "", d$name))], FREQ_LAB)
    ggplot(d, aes(time/24, value, colour = freq)) + geom_line(linewidth = 0.9) +
      scale_colour_viridis_d(option = "plasma", end = 0.9) +
      labs(x = "Days", y = ylab, colour = "Frequency", title = ttl) + theme_cio
  }

  output$p_ohc <- renderPlot(band_plot("OHC", "Fraction surviving",
                                       "Outer hair cells (no regeneration: a ratchet)"))
  output$p_gsh <- renderPlot(band_plot("GSH", "Glutathione (normalised)",
                                       "Reserve collapses base-first"))
  output$p_ihc <- renderPlot({
    s <- sim()
    d <- data.frame(time = rep(s$time, 2),
                    value = c(rowMeans(s[, paste0("IHC", 2:5)]),
                              rowMeans(s[, paste0("SGN", 2:5)])),
                    what = rep(c("inner hair cells", "spiral ganglion"),
                               each = nrow(s)))
    ggplot(d, aes(time/24, value, colour = what)) + geom_line(linewidth = 1) +
      labs(x = "Days", y = "Fraction surviving (speech frequencies)", colour = NULL,
           title = "The substrate of the speech-in-noise penalty") + theme_cio
  })
  output$p_ep <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/24)) +
      geom_line(aes(y = EP, colour = "endocochlear potential (mV)"), linewidth = 1) +
      geom_line(aes(y = INFL*85, colour = "inflammation index (scaled)"),
                linewidth = 1) +
      scale_y_continuous(sec.axis = sec_axis(~./85, name = "Inflammation index")) +
      labs(x = "Days", y = "EP (mV)", colour = NULL) + theme_cio
  })

  ## ---- tab 7 --------------------------------------------------------------
  output$p_grades <- renderPlot({
    s <- sim()
    d <- s[, c("time", "BROCK", "SIOP", "CTCAE")] %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) +
      geom_step(linewidth = 1) +
      labs(x = "Days", y = "Grade", colour = NULL,
           title = "Grades quantise a continuous variable into a staircase") +
      theme_cio
  })

  output$p_stair <- renderPlot({
    res <- lapply(1:10, function(n) {
      o <- as.data.frame(cio_run(mod_r(),
              cio_ev(n, input$dose, cycle_h = input$cyclen*24, bsa = input$bsa),
              cycles = n, cycle_h = input$cyclen*24))
      e <- o[nrow(o), ]
      data.frame(cum = n*input$dose, PTA = e$PTA, TS8k = e$TS7,
                 Brock = e$BROCK, SIOP = e$SIOP)
    })
    d <- do.call(rbind, res)
    ggplot(d, aes(cum)) +
      geom_line(aes(y = TS8k, colour = "8 kHz shift (dB)"), linewidth = 1) +
      geom_line(aes(y = PTA, colour = "speech PTA (dB)"), linewidth = 1) +
      geom_step(aes(y = Brock*15, colour = "Brock grade × 15"), linewidth = 1) +
      labs(x = "Cumulative dose (mg/m²)", y = "dB", colour = NULL,
           title = "Cumulative dose is the exposure variable") + theme_cio
  })

  output$t_grades <- renderDT({
    e <- last()
    datatable(data.frame(
      Scale = c("Brock", "SIOP Boston", "ASHA", "CTCAE v5.0"),
      Criterion = c("40 dB at 8/4/2/1 kHz", "20 dB, apex-ward march",
                    "≥20 dB at one, or ≥10 dB at two adjacent",
                    "shift at 2–8 kHz and speech PTA"),
      Result = c(e$BROCK, e$SIOP, ifelse(e$ASHA == 1, "positive", "negative"),
                 e$CTCAE)), options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 8 --------------------------------------------------------------
  delay_sweep <- reactive({
    delays <- c(0, 1, 2, 4, 6, 8, 12, 24)
    base <- as.data.frame(cio_run(mod_r(), cio_ev(input$cycles, input$dose,
                                                  bsa = input$bsa),
                                  cycles = input$cycles))
    b <- base[nrow(base), ]
    do.call(rbind, lapply(delays, function(d) {
      o <- as.data.frame(cio_run(mod_r(),
             cio_ev(input$cycles, input$dose, sts_cycles = seq_len(input$cycles) - 1L,
                    sts_gm2 = input$sts_dose, sts_delay = d, bsa = input$bsa),
             cycles = input$cycles))
      e <- o[nrow(o), ]
      data.frame(delay = d, PTA = e$PTA, TS8k = e$TS7,
                 otoprotection = (b$PTA - e$PTA)/b$PTA*100,
                 logkill_lost = (b$TUMLK - e$TUMLK)/b$TUMLK*100)
    }))
  })

  output$p_delay <- renderPlot({
    d <- delay_sweep()
    ggplot(d, aes(delay)) +
      geom_line(aes(y = otoprotection, colour = "otoprotection (%)"), linewidth = 1.1) +
      geom_line(aes(y = logkill_lost, colour = "tumour log-kill lost (%)"),
                linewidth = 1.1) +
      geom_point(aes(y = otoprotection, colour = "otoprotection (%)"), size = 2.4) +
      geom_point(aes(y = logkill_lost, colour = "tumour log-kill lost (%)"), size = 2.4) +
      geom_vline(xintercept = 6, linetype = "dashed", colour = "#B8860B") +
      scale_colour_manual(values = c("otoprotection (%)" = "#17A2A2",
                                     "tumour log-kill lost (%)" = "#C0392B")) +
      labs(x = "Delay between cisplatin and thiosulfate (h)", y = "%", colour = NULL,
           title = "The two curves separate because the cochlea and tumour read different clocks") +
      theme_cio
  })

  output$p_route <- renderPlot({
    n  <- as.data.frame(cio_run(mod_r(), cio_ev(input$cycles, input$dose,
                                                bsa = input$bsa),
                                cycles = input$cycles))
    sy <- as.data.frame(cio_run(mod_r(),
            cio_ev(input$cycles, input$dose, sts_cycles = seq_len(input$cycles) - 1L,
                   bsa = input$bsa), cycles = input$cycles))
    it <- as.data.frame(cio_run(mod_r(),
            cio_ev(input$cycles, input$dose, it_cycles = seq_len(input$cycles) - 1L,
                   it_umol = input$it_umol, bsa = input$bsa), cycles = input$cycles))
    g <- function(x) as.numeric(x[nrow(x), TS_COLS])
    b <- g(n)
    d <- data.frame(khz = rep(FREQ_HZ/1000, 2),
                    prot = c((b - g(sy))/b*100, (b - g(it))/b*100),
                    route = rep(c("systemic thiosulfate", "intratympanic"), each = 8))
    ggplot(d, aes(khz, prot, colour = route)) + geom_line(linewidth = 1.1) +
      geom_point(size = 2.4) +
      scale_x_log10(breaks = FREQ_HZ/1000, labels = FREQ_LAB) +
      scale_colour_manual(values = c("systemic thiosulfate" = "#4A90D9",
                                     "intratympanic" = "#17A2A2")) +
      labs(x = "Frequency (kHz)", y = "Protection (% of the unprotected shift)",
           colour = NULL,
           title = "The round window delivers a gradient that matches the damage gradient") +
      theme_cio
  })

  output$t_delay <- renderDT({
    datatable(delay_sweep() %>%
                mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 9 --------------------------------------------------------------
  output$p_gfr <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/24)) +
      geom_line(aes(y = GFROUT, colour = "GFR (mL/min)"), linewidth = 1) +
      geom_line(aes(y = TUBI*100, colour = "tubular injury × 100"), linewidth = 1) +
      labs(x = "Days", y = NULL, colour = NULL) + theme_cio
  })

  output$p_tum <- renderPlot({
    s <- sim()
    ggplot(s, aes(time/24)) +
      geom_line(aes(y = TUMAD, colour = "Pt-DNA adducts"), linewidth = 1) +
      geom_line(aes(y = TUMLK, colour = "cumulative log-kill"), linewidth = 1) +
      labs(x = "Days", y = NULL, colour = NULL) + theme_cio
  })

  output$t_amp <- renderTable({
    run <- function(p) {
      m <- param(mod_r(), p)
      o <- as.data.frame(cio_run(m, ev_r(), cycles = input$cycles,
                                 cycle_h = input$cyclen*24))
      o[nrow(o), ]
    }
    a <- run(list()); b <- run(list(KLOSS = 0, KPERM = 0))
    data.frame(Condition = c("renal feedback intact", "GFR clamped at baseline"),
               `GFR at end` = sprintf("%.1f", c(a$GFROUT, b$GFROUT)),
               `PTA (dB)` = sprintf("%.2f", c(a$PTA, b$PTA)),
               `8 kHz (dB)` = sprintf("%.2f", c(a$TS7, b$TS7)),
               check.names = FALSE)
  })

  ## ---- tab 10 -------------------------------------------------------------
  scen_tab <- reactive(cio_table(mod))

  output$t_scen <- renderDT({
    datatable(scen_tab() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(pageLength = 22, dom = "t"), rownames = FALSE)
  })

  output$p_scen <- renderPlot({
    d <- scen_tab()
    d$scenario <- factor(d$scenario, rev(d$scenario))
    ggplot(d, aes(PTA, scenario)) +
      geom_col(fill = "#6C63C7") +
      labs(x = "Speech-frequency pure-tone average (dB)", y = NULL,
           title = "Scenario library, end of treatment") + theme_cio
  })
}

shinyApp(ui, server)
