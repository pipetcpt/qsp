## =====================================================================
##  ANDROGENETIC ALOPECIA QSP — Shiny dashboard
##  10 tabs.  Requires aga_mrgsolve_model.R in the same directory.
##
##  The app is organised around the model's claim: AGA is a disease of
##  follicle SIZE, executed only at telogen exit, so every clinical number
##  is downstream of two things — how much inhibitory drive there is, and
##  how many cycles have happened.  Tab 3 is the one that matters.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("aga_mrgsolve_model.R")     # defines `mod` and the ev_* helpers

PAL <- c(placebo = "#7f7f7f", arm1 = "#1f77b4", arm2 = "#d62728",
         arm3 = "#2ca02c", arm4 = "#9467bd")
theme_set(theme_bw(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95")))

## ------------------------------------------------------------------ UI
ui <- fluidPage(
  titlePanel("Androgenetic alopecia — QSP model (50 ODE states)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("gs", "Genetic susceptibility GS", 0.5, 1.3, 0.90, 0.05),
      sliderInput("target", "Presenting hair count (per 5.07 cm²)",
                  500, 1250, 876, 25),
      helpText("The presenting patient is grown, not assumed: the untreated ",
               "model runs from an intact 18-year-old scalp until it reaches ",
               "this count. GS sets how long that takes."),
      radioButtons("sex", "Phenotype",
                   c("Male pattern" = "M", "Female pattern (post-menopausal)" = "F")),
      conditionalPanel("input.sex == 'F'",
        sliderInput("arind", "Androgen-INDEPENDENT drive (ARIND)",
                    0.0, 0.9, 0.62, 0.02)),
      hr(),
      h4("Arm 1"),
      selectInput("d1", NULL,
        c("none", "finasteride 1 mg", "finasteride 0.2 mg", "finasteride 5 mg",
          "dutasteride 0.5 mg", "topical finasteride 0.25%",
          "topical minoxidil 5% BID", "oral minoxidil 1 mg",
          "oral minoxidil 5 mg", "spironolactone 100 mg",
          "finasteride 1 mg + minoxidil 5%", "dutasteride + minoxidil 5%"),
        selected = "finasteride 1 mg"),
      h4("Arm 2"),
      selectInput("d2", NULL,
        c("none", "finasteride 1 mg", "finasteride 0.2 mg", "finasteride 5 mg",
          "dutasteride 0.5 mg", "topical finasteride 0.25%",
          "topical minoxidil 5% BID", "oral minoxidil 1 mg",
          "oral minoxidil 5 mg", "spironolactone 100 mg",
          "finasteride 1 mg + minoxidil 5%", "dutasteride + minoxidil 5%"),
        selected = "topical minoxidil 5% BID"),
      hr(),
      sliderInput("sult", "Follicular SULT1A1 activity", 0.1, 2.5, 2.0, 0.05),
      helpText("Minoxidil is a prodrug. Below ~0.5 the patient is a ",
               "non-responder and the trial mean is a mixture of two populations."),
      sliderInput("adh", "Adherence (doses taken per 7)", 1, 7, 7, 1),
      sliderInput("stopd", "Stop treatment at (months, 0 = never)",
                  0, 60, 0, 6),
      sliderInput("years", "Follow-up (years)", 1, 5, 3, 1),
      checkboxInput("seti", "Block CRTH2 (setipiprant)", FALSE),
      actionButton("go", "Run", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("1 · Patient",        plotOutput("p_pat", height = 520),
                                       tableOutput("t_pat")),
        tabPanel("2 · Drug PK",        plotOutput("p_pk",  height = 560)),
        tabPanel("3 · Target & DHT",   plotOutput("p_dht", height = 560),
                                       tableOutput("t_dht")),
        tabPanel("4 · Hair cycle",     plotOutput("p_cyc", height = 560)),
        tabPanel("5 · Miniaturisation",plotOutput("p_min", height = 560)),
        tabPanel("6 · Clinical endpoints", plotOutput("p_end", height = 560),
                                       tableOutput("t_end")),
        tabPanel("7 · Shedding",       plotOutput("p_shed", height = 520),
                                       htmlOutput("x_shed")),
        tabPanel("8 · Scenario table", tableOutput("t_scen")),
        tabPanel("9 · Systemic effects", plotOutput("p_sys", height = 560),
                                       tableOutput("t_sys")),
        tabPanel("10 · Dose–response", plotOutput("p_dr", height = 560),
                                       tableOutput("t_dr"))
      )
    )
  )
)

## -------------------------------------------------------------- server
server <- function(input, output, session) {

  patient <- eventReactive(input$go, {
    m <- mod %>% param(GS = input$gs, SULT = input$sult,
                       SETI = as.numeric(input$seti))
    if (input$sex == "F")
      m <- m %>% param(ARIND = input$arind, LOC5AR = 0.30, AROM = 1.6,
                       W_LOC = 0.30, W_SER = 0.30, W_ALT = 0.40)
    b <- burn_in(m, gs = input$gs, target = input$target)
    list(mod = m, y0 = b, age = b$age, trace = b$out)
  }, ignoreNULL = FALSE)

  make_ev <- function(lab, days) {
    if (lab == "none") return(NULL)
    switch(lab,
      "finasteride 1 mg"    = ev_fin(1, days),
      "finasteride 0.2 mg"  = ev_fin(0.2, days),
      "finasteride 5 mg"    = ev_fin(5, days),
      "dutasteride 0.5 mg"  = ev_dut(0.5, days),
      "topical finasteride 0.25%" = ev_topfin(2.5, days),
      "topical minoxidil 5% BID"  = ev_mxt(50, days),
      "oral minoxidil 1 mg" = ev_mxo(1, days),
      "oral minoxidil 5 mg" = ev_mxo(5, days),
      "spironolactone 100 mg" = ev_spi(100, days),
      "finasteride 1 mg + minoxidil 5%" = c(ev_fin(1, days), ev_mxt(50, days)),
      "dutasteride + minoxidil 5%"      = c(ev_dut(0.5, days), ev_mxt(50, days)))
  }

  sims <- eventReactive(input$go, {
    P <- patient(); end <- input$years * 365
    days <- if (input$stopd > 0) min(round(input$stopd*30.4), end) else end
    thin <- function(e) {           # adherence filter: keep k of every 7 days
      if (is.null(e) || input$adh >= 7) return(e)
      e <- as.data.frame(e)
      e[(floor(e$time) %% 7) < input$adh, , drop = FALSE]
    }
    run <- function(lab) {
      e <- make_ev(lab, days)
      if (!is.null(e) && grepl("finasteride|dutasteride|spirono", lab)) e <- thin(e)
      as_tibble(run_arm(P$mod, P$y0, e, end = end)) %>% mutate(arm = lab)
    }
    bind_rows(run("none") %>% mutate(arm = "placebo"),
              run(input$d1), run(input$d2)) %>%
      mutate(arm = factor(arm, levels = unique(arm)), year = time/365)
  }, ignoreNULL = FALSE)

  base <- reactive(sims() %>% group_by(arm) %>% slice(1) %>% ungroup())
  delta <- reactive({
    b <- base() %>% select(arm, TAHC0 = TAHC, HMI0 = HMI)
    sims() %>% left_join(b, by = "arm") %>%
      mutate(dTAHC = TAHC - TAHC0, dHMI = HMI - HMI0)
  })

  ## --- 1 patient -----------------------------------------------------
  output$p_pat <- renderPlot({
    tr <- patient()$trace %>% mutate(age = 18 + time/365)
    tr %>% select(age, TAHC, HMI, VELL, LOST) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value)) + geom_line(linewidth = 1, colour = "#1f77b4") +
      geom_vline(xintercept = patient()$age, linetype = 2, colour = "#d62728") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (years)", y = NULL,
           title = "Natural history from an intact 18-year-old scalp",
           subtitle = sprintf("GS = %.2f reaches the presenting count of %d hairs at age %.1f years",
                              input$gs, input$target, patient()$age))
  })
  output$t_pat <- renderTable({
    y <- patient()$y0; s <- setNames(as.list(y$state), mod@cmtn)
    data.frame(
      quantity = c("presenting age (y)", "target-area hair count",
                   "terminal", "intermediate", "small intermediate",
                   "vellus — recent (RECOVERABLE)", "vellus — established",
                   "follicles lost (irreversible)", "hair mass index",
                   "telogen %"),
      value = c(sprintf("%.1f", y$age), sprintf("%.0f", y$tahc),
        sprintf("%.0f", s$A1+s$C1+s$T1), sprintf("%.0f", s$A2+s$C2+s$T2),
        sprintf("%.0f", s$A3+s$C3+s$T3), sprintf("%.0f", s$A4A+s$C4A+s$T4A),
        sprintf("%.0f", s$A4B+s$C4B+s$T4B), sprintf("%.0f", s$LOST),
        sprintf("%.1f", tail(patient()$trace$HMI[patient()$trace$TAHC >= y$tahc],1)),
        sprintf("%.1f", tail(patient()$trace$TELPC[patient()$trace$TAHC >= y$tahc],1))))
  })

  ## --- 2 PK ----------------------------------------------------------
  output$p_pk <- renderPlot({
    sims() %>% filter(time <= 40) %>%
      select(time, arm, FINC, FINT, DUTC, MXSF, MXDC, SPIC) %>%
      pivot_longer(-c(time, arm)) %>% filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = "concentration (ng/mL) or follicular units",
           title = "Drug pharmacokinetics — first 40 days",
           subtitle = "Finasteride plasma swings ~9-fold within a day. Nothing downstream does.")
  })

  ## --- 3 target and DHT ----------------------------------------------
  output$p_dht <- renderPlot({
    sims() %>% select(time, arm, E2S, E1S, DHTC, DHTS, ARN, DKK1) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(name, c("E2S","E1S","DHTC","DHTS","ARN","DKK1"))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = "fraction of baseline",
           title = "The target integrates: free 5AR, then DHT, then AR signal",
           subtitle = "E2S/E1S = free scalp type II / type I 5-alpha-reductase")
  })
  output$t_dht <- renderTable({
    d <- sims() %>% filter(abs(time - 365) < 3) %>% group_by(arm) %>% slice(1)
    data.frame(arm = as.character(d$arm),
      `scalp DHT %` = sprintf("%+.1f", (d$DHTC-1)*100),
      `serum DHT %` = sprintf("%+.1f", (d$DHTS-1)*100),
      `free scalp 5AR-II` = sprintf("%.4f", d$E2S),
      `free scalp 5AR-I`  = sprintf("%.3f", d$E1S),
      `testosterone %`    = sprintf("%+.1f", (d$TST-1)*100),
      check.names = FALSE)
  })

  ## --- 4 hair cycle ---------------------------------------------------
  output$p_cyc <- renderPlot({
    sims() %>%
      mutate(anagen = A1+A2+A3+A4A+A4B, catagen = C1+C2+C3+C4A+C4B,
             telogen = T1+T2+T3+T4A+T4B) %>%
      select(time, arm, anagen, catagen, telogen, TELPC) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = "follicles per 5.07 cm² (TELPC in %)",
           title = "Cycle-phase populations",
           subtitle = "Telogen % is an anagen-duration assay: T_T/(T_A+T_C+T_T)")
  })

  ## --- 5 miniaturisation ----------------------------------------------
  output$p_min <- renderPlot({
    sims() %>%
      transmute(time, arm,
                terminal = N1, intermediate = N2, `small intermediate` = N3,
                `vellus (recent)` = N4A, `vellus (established)` = N4B,
                lost = LOST) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(name, c("terminal","intermediate","small intermediate",
                                   "vellus (recent)","vellus (established)","lost"))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = "follicles per 5.07 cm²",
           title = "The miniaturisation ratchet, class by class",
           subtitle = "Only the RECENT vellus pool is recoverable; it drains one way into 'established' and 'lost'")
  })

  ## --- 6 clinical endpoints -------------------------------------------
  output$p_end <- renderPlot({
    delta() %>% select(year, arm, dTAHC, dHMI, MDIA, VELL) %>%
      pivot_longer(-c(year, arm)) %>%
      mutate(name = recode(name, dTAHC = "Δ target-area hair count",
                           dHMI = "Δ hair mass index",
                           MDIA = "mean shaft diameter (µm)",
                           VELL = "vellus follicles")) %>%
      ggplot(aes(year, value, colour = arm)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_line(linewidth = 1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "years", y = NULL,
           title = "Clinical endpoints",
           subtitle = "Mass moves ~1.9× as far as count, because mass goes as d²")
  })
  output$t_end <- renderTable({
    pick <- function(d) delta() %>% filter(abs(time - d) < 3) %>%
      group_by(arm) %>% slice(1) %>% ungroup()
    a <- pick(182); b <- pick(365); cc <- pick(730)
    pb <- b %>% filter(arm == "placebo") %>% pull(TAHC)
    data.frame(arm = as.character(b$arm),
      `Δ count 6 m`  = sprintf("%+.1f", a$dTAHC),
      `Δ count 12 m` = sprintf("%+.1f", b$dTAHC),
      `Δ count 24 m` = if (nrow(cc)) sprintf("%+.1f", cc$dTAHC) else NA,
      `vs placebo 12 m` = sprintf("%+.1f", b$TAHC - pb),
      `Δ mass 12 m %` = sprintf("%+.1f", 100*b$dHMI/b$HMI0),
      check.names = FALSE)
  })

  ## --- 7 shedding ------------------------------------------------------
  output$p_shed <- renderPlot({
    sims() %>% filter(time <= 400) %>%
      ggplot(aes(time, SHEDSC, colour = arm)) + geom_line(linewidth = 1) +
      labs(x = "days", y = "whole-scalp-equivalent shed (hairs/day)",
           title = "The initial shed is the mechanism working",
           subtitle = "Minoxidil shortens telogen; telogen exit requires releasing the club hair first")
  })
  output$x_shed <- renderUI(HTML(paste0(
    "<p><b>Why 100 hairs a day is normal.</b> Daily shed = ",
    "N<sub>total</sub> × f<sub>telogen</sub> / T<sub>telogen</sub> = ",
    "100&nbsp;000 × 0.090 / 100&nbsp;d = <b>90 hairs/day</b>. The textbook number is ",
    "an identity between three parameters, not an observation. In AGA the telogen ",
    "fraction rises only because anagen is shorter in the smaller classes, so ",
    "'increased shedding' is a consequence of miniaturisation, not a second disease.</p>")))

  ## --- 8 scenario table -------------------------------------------------
  output$t_scen <- renderTable({
    P <- patient(); end <- input$years*365
    labs <- c("finasteride 1 mg","finasteride 0.2 mg","finasteride 5 mg",
              "dutasteride 0.5 mg","topical finasteride 0.25%",
              "topical minoxidil 5% BID","oral minoxidil 5 mg",
              "finasteride 1 mg + minoxidil 5%","dutasteride + minoxidil 5%")
    p0 <- as_tibble(run_arm(P$mod, P$y0, NULL, end = end))
    g  <- function(df, d) df$TAHC[which.min(abs(df$time - d))]
    do.call(rbind, lapply(labs, function(l) {
      r <- as_tibble(run_arm(P$mod, P$y0, make_ev(l, end), end = end))
      data.frame(arm = l,
        `Δ 6 m`  = sprintf("%+.1f", g(r,182)  - r$TAHC[1]),
        `Δ 12 m` = sprintf("%+.1f", g(r,365)  - r$TAHC[1]),
        `vs placebo 12 m` = sprintf("%+.1f", g(r,365) - g(p0,365)),
        `vs placebo, end` = sprintf("%+.1f", g(r,end) - g(p0,end)),
        check.names = FALSE)
    }))
  })

  ## --- 9 systemic ------------------------------------------------------
  output$p_sys <- renderPlot({
    sims() %>% select(time, arm, PSA, TST, EST, SEB, MAP, HTR) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = "normalised (MAP in mmHg)",
           title = "Systemic and off-target consequences",
           subtitle = "PSA halves on a 5ARI — the 'double the value' rule is this curve")
  })
  output$t_sys <- renderTable({
    d <- sims() %>% filter(abs(time - 365) < 3) %>% group_by(arm) %>% slice(1)
    data.frame(arm = as.character(d$arm),
      `PSA (× baseline)` = sprintf("%.2f", d$PSA),
      `PSA-adjusted`     = sprintf("%.2f", d$PSAADJ),
      `testosterone`     = sprintf("%.2f", d$TST),
      `estradiol`        = sprintf("%.2f", d$EST),
      `sebum`            = sprintf("%.2f", d$SEB),
      `MAP (mmHg)`       = sprintf("%.1f", d$MAP),
      `hypertrichosis`   = sprintf("%.1f", d$HTR),
      `sexual AE (%)`    = sprintf("%.1f", d$SEXPCT),
      check.names = FALSE)
  })

  ## --- 10 dose-response --------------------------------------------------
  output$p_dr <- renderPlot({
    P <- patient()
    doses <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)
    dr <- do.call(rbind, lapply(doses, function(d) {
      r <- as_tibble(run_arm(P$mod, P$y0, ev_fin(d, 60), end = 50))
      k <- which.min(abs(r$time - 42))
      data.frame(dose = d, scalp = (r$DHTC[k]-1)*100, serum = (r$DHTS[k]-1)*100)
    }))
    obsv <- data.frame(dose = c(0.01,0.05,0.2,1,5),
                       scalp = c(-14.9,-61.6,-56.5,-64.1,-69.4),
                       serum = c(NA,-49.5,-68.6,-71.4,-72.2))
    ggplot(dr, aes(dose)) +
      geom_line(aes(y = scalp, colour = "scalp DHT (model)"), linewidth = 1) +
      geom_line(aes(y = serum, colour = "serum DHT (model)"), linewidth = 1) +
      geom_point(data = obsv, aes(y = scalp, colour = "scalp DHT (Drake 1999)"), size = 3) +
      geom_point(data = obsv, aes(y = serum, colour = "serum DHT (Drake 1999)"), size = 3) +
      scale_x_log10(breaks = doses) +
      labs(x = "finasteride dose (mg/day, log scale)", y = "change at day 42 (%)",
           colour = NULL,
           title = "The dose–response is flat because the target integrates",
           subtitle = "0.2 mg already buys ~90% of what 5 mg buys. Points: Drake 1999, PMID 10495374")
  })
  output$t_dr <- renderTable({
    data.frame(
      dose = c("0.01","0.05","0.2","1","5"),
      `scalp DHT — observed` = c("-14.9","-61.6","-56.5","-64.1","-69.4"),
      `scalp DHT — model`    = c("-18.3","-43.1","-57.5","-63.6","-67.4"),
      `serum DHT — observed` = c("n/a","-49.5","-68.6","-71.4","-72.2"),
      `serum DHT — model`    = c("-20.1","-47.3","-63.1","-70.3","-74.9"),
      check.names = FALSE)
  })
}

shinyApp(ui, server)
