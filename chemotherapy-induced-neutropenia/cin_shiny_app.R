##  Chemotherapy-Induced Neutropenia / Febrile Neutropenia — Shiny dashboard
##  ============================================================================
##  Companion app for cin_mrgsolve_model.R.
##
##  DESIGN RULE, and it is the reason the tabs are arranged the way they are:
##  the ANC is never shown on its own.  Every panel that displays a neutrophil
##  count also displays the quantity that count is standing in for --
##
##      the DURATION below 0.5, not the nadir           (tab 3)
##      the effective defence, count x function         (tab 6)
##      the marrow output the count is over-reading     (tab 5)
##      the dose intensity the count is costing         (tab 7)
##
##  A dashboard that plots ANC(t) and stops is the visual equivalent of quoting
##  a nadir, which is the specific mistake the model exists to argue against.
##
##  Run with:
##      library(shiny); library(mrgsolve); library(ggplot2); library(dplyr)
##      shiny::runApp("cin_shiny_app.R")
##
##  Disclaimer: research / education / hypothesis generation only.
##  ============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

MOD <- mread("cin_mrgsolve_model.R")

## ---------------------------------------------------------------------------
## Agent and regimen library (mirrors cin_mrgsolve_model.R's usage block)
## ---------------------------------------------------------------------------
CIN_AGENTS <- list(
  docetaxel = list(slot="A", mw=807.9, per_bsa=TRUE, tinf=1.0,
                   CL=21.6, V1=7.4, V2=25, V3=42, Q2=45, Q3=5,
                   slope=108.7, slopeP=0.020, slopeE=0.35,
                   fcyc=0.90, hsc=0.55, muc=0.85, cthr=0, taxane=TRUE),
  paclitaxel = list(slot="A", mw=853.9, per_bsa=TRUE, tinf=3.0,
                   CL=14, V1=10, V2=30, V3=180, Q2=40, Q3=8,
                   slope=4.468, slopeP=0.018, slopeE=0.30,
                   fcyc=0.90, hsc=0.35, muc=0.55, cthr=0.05, taxane=TRUE),
  doxorubicin = list(slot="B", mw=543.5, per_bsa=TRUE, tinf=0.25,
                   CL=50, V1=12, V2=700, Q2=38,
                   slope=25.52, slopeP=0.055, slopeE=0.55,
                   fcyc=0.60, hsc=0.75, muc=1.00, cthr=0),
  gemcitabine = list(slot="B", mw=263.2, per_bsa=TRUE, tinf=0.5,
                   CL=92, V1=14, V2=6, Q2=8,
                   slope=0.5853, slopeP=0.360, slopeE=0.60,
                   fcyc=0.95, hsc=0.30, muc=0.35, cthr=0),
  cyclophosphamide = list(slot="C", mw=261.1, per_bsa=FALSE, tinf=0.5,
                   CL=5.6, V1=32, V2=18, Q2=6,
                   slope=0.0463, slopeP=0.070, slopeE=0.60,
                   fcyc=0.40, hsc=1.00, muc=0.55, cthr=0),
  carboplatin = list(slot="C", mw=371.3, per_bsa=FALSE, tinf=0.5,
                   CL=NA, V1=16, V2=8, Q2=3,
                   slope=0.1046, slopeP=0.9375, slopeE=0.70,
                   fcyc=0.35, hsc=0.85, muc=0.40, cthr=0, auc_dosed=TRUE),
  etoposide = list(slot="D", mw=588.6, per_bsa=TRUE, tinf=1.0,
                   CL=1.35, V1=8, V2=10, Q2=2.2,
                   slope=0.1449, slopeP=0.060, slopeE=0.35,
                   fcyc=0.92, hsc=0.40, muc=0.60, cthr=0),
  topotecan = list(slot="A", mw=421.4, per_bsa=TRUE, tinf=0.5,
                   CL=26, V1=17, V2=30, V3=60, Q2=18, Q3=3,
                   slope=555.2, slopeP=0.090, slopeE=0.55,
                   fcyc=0.95, hsc=0.55, muc=0.65, cthr=0)
)

CIN_REGIMENS <- list(
  "Docetaxel 100 mg/m2 q3w (Vogel 2005)" = list(cycle=504, drugs=list(
      list(agent="docetaxel", dose=100, days=0))),
  "Docetaxel 75 mg/m2 q3w" = list(cycle=504, drugs=list(
      list(agent="docetaxel", dose=75, days=0))),
  "AT: doxorubicin 60 + docetaxel 75 (Green 2003)" = list(cycle=504, drugs=list(
      list(agent="docetaxel", dose=75, days=0),
      list(agent="doxorubicin", dose=60, days=0))),
  "CAE: cyclo 1000 + doxo 45 + etop 100 d1-3 (Crawford 1991)" =
    list(cycle=504, drugs=list(
      list(agent="cyclophosphamide", dose=1000, days=0),
      list(agent="doxorubicin", dose=45, days=0),
      list(agent="etoposide", dose=100, days=c(0,24,48)))),
  "TAC: docetaxel 75 + doxo 50 + cyclo 500" = list(cycle=504, drugs=list(
      list(agent="docetaxel", dose=75, days=0),
      list(agent="doxorubicin", dose=50, days=0),
      list(agent="cyclophosphamide", dose=500, days=0))),
  "AC: doxorubicin 60 + cyclophosphamide 600 q3w" = list(cycle=504, drugs=list(
      list(agent="doxorubicin", dose=60, days=0),
      list(agent="cyclophosphamide", dose=600, days=0))),
  "Dose-dense AC q2w" = list(cycle=336, drugs=list(
      list(agent="doxorubicin", dose=60, days=0),
      list(agent="cyclophosphamide", dose=600, days=0))),
  "EP: carboplatin AUC5 + etoposide 100 d1-3 (SCLC)" = list(cycle=504, drugs=list(
      list(agent="carboplatin", dose=5, days=0, auc=TRUE),
      list(agent="etoposide", dose=100, days=c(0,24,48)))),
  "Gemcitabine 1000 d1,8 + carboplatin AUC5" = list(cycle=504, drugs=list(
      list(agent="gemcitabine", dose=1000, days=c(0,168)),
      list(agent="carboplatin", dose=5, days=0, auc=TRUE))),
  "Topotecan 1.5 mg/m2 d1-5 q3w" = list(cycle=504, drugs=list(
      list(agent="topotecan", dose=1.5, days=c(0,24,48,72,96)))),
  "Paclitaxel 175 mg/m2 q3w (3 h)" = list(cycle=504, drugs=list(
      list(agent="paclitaxel", dose=175, days=0))),
  "Paclitaxel 80 mg/m2 weekly" = list(cycle=504, drugs=list(
      list(agent="paclitaxel", dose=80, days=c(0,168,336))))
)
SLOT_CMT <- c(A = 1, B = 4, C = 6, D = 8)

slot_params <- function(a, bsa, gfr) {
  s <- a$slot; scl <- if (isTRUE(a$per_bsa)) bsa else 1
  CL <- if (isTRUE(a$auc_dosed)) (gfr + 25) * 60/1000 else a$CL * scl
  p <- list()
  p[[paste0(s,"_CL")]] <- CL
  p[[paste0(s,"_V1")]] <- a$V1 * scl
  p[[paste0(s,"_V2")]] <- a$V2 * scl
  p[[paste0(s,"_Q2")]] <- a$Q2 * scl
  if (s == "A") { p$A_V3 <- a$V3 * scl; p$A_Q3 <- a$Q3 * scl }
  p[[paste0(s,"_SLOPE")]]  <- a$slope
  p[[paste0(s,"_SLOPEP")]] <- a$slopeP
  p[[paste0(s,"_SLOPEE")]] <- a$slopeE
  p[[paste0(s,"_FCYC")]]   <- a$fcyc
  p[[paste0(s,"_HSC")]]    <- a$hsc
  p[[paste0(s,"_MUC")]]    <- a$muc
  p[[paste0(s,"_CTHR")]]   <- a$cthr
  p
}

build_run <- function(reg_name, n_cycles, bsa, wt, gfr, anc0, plt0, hb0,
                      reserve, sens, mttmult, gammult,
                      gcsf, gcsf_start, gcsf_days, tril, dex, abx, dose_mult,
                      delta = 2) {
  reg <- CIN_REGIMENS[[reg_name]]
  pars <- list(BSA = bsa, WT = wt, GFR = gfr, ANC0 = anc0, PLT0 = plt0,
               HB0 = hb0, RESERVE = reserve, SENS = sens,
               MTTMULT = mttmult, GAMMULT = gammult)
  ev <- NULL
  add <- function(e) ev <<- if (is.null(ev)) e else rbind(ev, e)
  has_taxane <- FALSE
  for (d in reg$drugs) {
    a <- CIN_AGENTS[[d$agent]]
    if (isTRUE(a$taxane)) has_taxane <- TRUE
    pars <- modifyList(pars, slot_params(a, bsa, gfr))
    dose_mg <- if (isTRUE(d$auc)) d$dose * (gfr + 25) else d$dose * bsa
    amt <- dose_mg * 1000 / a$mw * dose_mult
    for (ci in seq_len(n_cycles) - 1) for (dd in d$days)
      add(data.frame(ID=1, time=ci*reg$cycle + dd, cmt=SLOT_CMT[[a$slot]],
                     amt=amt, rate=amt/a$tinf, evid=1))
  }
  if (!has_taxane) pars$SENS_TAXANE <- 1
  if (gcsf != "none") {
    pars$GTYPE <- if (gcsf == "filgrastim") 1 else 2
    amt <- if (gcsf == "filgrastim") 5 * wt * 1000 else 6e6
    nd  <- if (gcsf == "filgrastim") gcsf_days else 1
    for (ci in seq_len(n_cycles) - 1) for (k in seq_len(nd) - 1)
      add(data.frame(ID=1, time=ci*reg$cycle + gcsf_start + 24*k, cmt=10,
                     amt=amt, rate=0, evid=1))
  }
  if (tril) for (d in reg$drugs) for (ci in seq_len(n_cycles) - 1)
    for (dd in unique(d$days))
      add(data.frame(ID=1, time=ci*reg$cycle + dd - 0.5, cmt=14,
                     amt=240*bsa, rate=240*bsa/0.5, evid=1))
  if (dex) for (ci in seq_len(n_cycles) - 1) for (h in c(0,12,24,36,48))
    add(data.frame(ID=1, time=ci*reg$cycle + h, cmt=16, amt=8, rate=0, evid=1))
  if (abx) for (ci in seq_len(n_cycles) - 1) for (k in 0:10)
    add(data.frame(ID=1, time=ci*reg$cycle + 120 + 24*k, cmt=17, amt=500,
                   rate=0, evid=1))
  ev <- ev[ev$time >= 0, ]
  ev <- ev[order(ev$time), ]
  end <- (n_cycles - 1) * reg$cycle + reg$cycle
  MOD %>% param(pars) %>% data_set(ev) %>%
    mrgsim(end = end, delta = delta, hmax = 0.25) %>% as.data.frame()
}

summarise_cycle <- function(df, cycle_h) {
  df %>% mutate(cyc = floor(time / cycle_h) + 1) %>%
    group_by(cyc) %>%
    summarise(`nadir ANC` = round(min(ANC), 3),
              `nadir day` = round(time[which.min(ANC)] %% cycle_h / 24, 1),
              `DSN (d)` = round(max(DSN_D) - min(DSN_D), 2),
              `days ANC<1` = round(sum(ANC < 1) * (time[2]-time[1]) / 24, 2),
              `P(FN) cycle` = round(1 - exp(-(max(CUMHAZ) - min(CUMHAZ))), 3),
              `nadir plt` = round(min(PLT), 0),
              `nadir Hb` = round(min(HB), 1),
              `peak temp` = round(max(TEMP), 2),
              `peak CRP` = round(max(CRP), 1),
              `min barrier` = round(min(BAR), 2),
              `min function` = round(min(FUNC), 2),
              .groups = "drop")
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Chemotherapy-Induced Neutropenia / Febrile Neutropenia — QSP explorer"),
  tags$p(style = "color:#555;",
    HTML("The ANC is never shown alone. Each panel pairs the count with what it ",
         "is standing in for: the <b>duration</b> below 0.5, the <b>effective ",
         "defence</b> (count x function), the <b>marrow output</b> it over-reads ",
         "on a growth factor, and the <b>dose intensity</b> it costs.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("regimen", "Regimen", names(CIN_REGIMENS),
                  selected = names(CIN_REGIMENS)[4]),
      sliderInput("n_cycles", "Cycles", 1, 6, 1, step = 1),
      sliderInput("dose_mult", "Delivered dose (fraction of planned)",
                  0.5, 1.2, 1.0, step = 0.05),
      tags$hr(),
      tags$b("Growth factor"),
      selectInput("gcsf", NULL,
                  c("none", "pegfilgrastim", "filgrastim", "eflapegrastim"),
                  selected = "none"),
      sliderInput("gcsf_start", "Start (h after day-1 dose)", 0, 240, 24,
                  step = 6),
      sliderInput("gcsf_days", "Daily doses (filgrastim only)", 1, 14, 10,
                  step = 1),
      tags$hr(),
      tags$b("Other interventions"),
      checkboxInput("tril", "Trilaciclib 240 mg/m2 before each dose", FALSE),
      checkboxInput("dex", "Dexamethasone 8 mg BID x 3 d", FALSE),
      checkboxInput("abx", "Levofloxacin 500 mg daily d5-15", FALSE),
      tags$hr(),
      tags$b("Patient"),
      sliderInput("bsa", "BSA (m2)", 1.4, 2.2, 1.75, step = 0.05),
      sliderInput("wt", "Weight (kg)", 40, 130, 70, step = 1),
      sliderInput("gfr", "GFR (mL/min)", 25, 130, 90, step = 5),
      sliderInput("anc0", "Baseline ANC", 1.5, 9.0, 5.0, step = 0.1),
      sliderInput("plt0", "Baseline platelets", 100, 450, 250, step = 10),
      sliderInput("hb0", "Baseline Hb (g/dL)", 8, 17, 13.5, step = 0.1),
      sliderInput("reserve", "Marrow reserve (prior chemo/RT)", 0.5, 1.0, 1.0,
                  step = 0.05),
      sliderInput("sens", "Myelotoxic sensitivity", 0.5, 2.0, 1.0, step = 0.05),
      sliderInput("mttmult", "Transit-time multiplier", 0.7, 1.4, 1.0,
                  step = 0.05),
      sliderInput("gammult", "Feedback-strength multiplier", 0.6, 1.5, 1.0,
                  step = 0.05)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1. Patient & summary",
          h4("Cycle-by-cycle summary"),
          tableOutput("cyc_tbl"),
          h4("Course-level endpoints"),
          tableOutput("course_tbl"),
          tags$p(style="color:#666;",
            "The nadir is shown because everyone asks for it. The column that",
            " drives the febrile-neutropenia hazard is DSN.")),
        tabPanel("2. Pharmacokinetics",
          h4("Cytotoxic and growth-factor concentrations"),
          plotOutput("pk_plot", height = 320),
          h4("Pegfilgrastim clears itself"),
          plotOutput("tmdd_plot", height = 300),
          tags$p(style="color:#666;",
            "PEGylation abolishes the renal route, so the only clearance left",
            " is internalisation through receptors on neutrophils. The G-CSF",
            " concentration therefore collapses exactly as the ANC recovers --",
            " the drug is cleared by the cells it created, which is why a",
            " single fixed 6 mg dose self-titrates and needs no weight",
            " adjustment.")),
        tabPanel("3. ANC and the duration below 0.5",
          h4("ANC with the grade-4 threshold and the shaded time below it"),
          plotOutput("anc_plot", height = 340),
          h4("Cumulative days below 0.5 (DSN) and the FN hazard it generates"),
          plotOutput("dsn_plot", height = 280),
          tags$p(style="color:#666;",
            "The nadir is a point; the hazard is an area. Two patients with",
            " the same lowest value can carry very different risk.")),
        tabPanel("4. Marrow compartments",
          h4("The delay line: proliferative pool through to circulating"),
          plotOutput("chain_plot", height = 360),
          h4("Storage and marginated pools, and the stem reserve"),
          plotOutput("pool_plot", height = 300),
          tags$p(style="color:#666;",
            "The kill happens in the proliferative pool on days 0-2. What the",
            " blood count shows on day 8-10 is that event after it has",
            " traversed a fixed-length pipeline.")),
        tabPanel("5. Count vs marrow output",
          h4("Measured ANC against the flux actually leaving the marrow"),
          plotOutput("output_plot", height = 340),
          tags$p(style="color:#666;",
            "On a growth factor the two separate. Part of the rise is",
            " production, part is release of the storage pool, and part is",
            " simply slower clearance of cells that already existed. Only the",
            " first is marrow recovery.")),
        tabPanel("6. Count vs function",
          h4("ANC, neutrophil function, and the effective defence"),
          plotOutput("func_plot", height = 340),
          h4("Bacterial load, CRP and temperature"),
          plotOutput("infect_plot", height = 300),
          tags$p(style="color:#666;",
            "Corticosteroid raises the count and lowers the function.",
            " Levofloxacin changes no haematological number at all and still",
            " moves the hazard, because it acts on the other side of the",
            " race.")),
        tabPanel("7. Other lineages",
          h4("Three lineages, three clocks"),
          plotOutput("lineage_plot", height = 380),
          tableOutput("lineage_tbl"),
          tags$p(style="color:#666;",
            "Same insult, same term in the equations. What separates the nadir",
            " days is the transit time and the cell lifespan: 7.9 h for a",
            " neutrophil, 10 days for a platelet, 120 days for a red cell.",
            " That is why anaemia is the cumulative toxicity and neutropenia",
            " is the acute one.")),
        tabPanel("8. Dose intensity",
          h4("Delivered dose intensity and tumour burden"),
          plotOutput("rdi_plot", height = 320),
          tableOutput("rdi_tbl"),
          tags$p(style="color:#666;",
            "The tumour layer is a deliberately thin log-kill. Read it as a",
            " price list for dose reduction, not as a survival prediction --",
            " the dose-intensity/survival association is observational.")),
        tabPanel("9. Scenario comparison",
          h4("Compare interventions on the current regimen and patient"),
          checkboxGroupInput("cmp", NULL, inline = TRUE,
            choices = c("no support", "pegfilgrastim d2", "pegfilgrastim same-day",
                        "filgrastim d2-11", "filgrastim from d8",
                        "trilaciclib", "levofloxacin", "dose 80%", "dose 60%"),
            selected = c("no support", "pegfilgrastim d2", "dose 80%")),
          plotOutput("cmp_plot", height = 340),
          tableOutput("cmp_tbl"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  cycle_h <- reactive(CIN_REGIMENS[[input$regimen]]$cycle)

  sim <- reactive({
    build_run(input$regimen, input$n_cycles, input$bsa, input$wt, input$gfr,
              input$anc0, input$plt0, input$hb0, input$reserve, input$sens,
              input$mttmult, input$gammult, input$gcsf, input$gcsf_start,
              input$gcsf_days, input$tril, input$dex, input$abx,
              input$dose_mult)
  })

  output$cyc_tbl <- renderTable(summarise_cycle(sim(), cycle_h()))

  output$course_tbl <- renderTable({
    df <- sim()
    data.frame(
      quantity = c("total days ANC < 0.5", "P(at least one FN episode)",
                   "worst ANC", "worst platelets", "final Hb (g/dL)",
                   "marrow reserve at end", "tumour burden at end",
                   "delivered dose intensity"),
      value = c(round(max(df$DSN_D), 2),
                round(1 - exp(-max(df$CUMHAZ)), 3),
                round(min(df$ANC), 3), round(min(df$PLT), 0),
                round(tail(df$HB, 1), 1),
                round(tail(df$RESERVE_N, 1), 3),
                round(tail(df$TUMOUR, 1), 2),
                paste0(round(input$dose_mult * 100), "%")))
  })

  output$pk_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `slot A` = CDOCE, `G-CSF (ng/mL)` = CGCSF) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL) + theme_bw()
  })

  output$tmdd_plot <- renderPlot({
    df <- sim()
    if (all(df$CGCSF == 0)) {
      return(ggplot() + annotate("text", 0, 0,
        label = "Select a growth factor to see the target-mediated loop") +
        theme_void())
    }
    sc <- max(df$CGCSF) / max(df$ANC)
    ggplot(df, aes(time/24)) +
      geom_line(aes(y = CGCSF, colour = "G-CSF (ng/mL)"), linewidth = 0.8) +
      geom_line(aes(y = ANC * sc, colour = "ANC (rescaled)"), linewidth = 0.8) +
      scale_y_continuous(sec.axis = sec_axis(~./sc, name = "ANC (10^9/L)")) +
      labs(x = "days", y = "G-CSF (ng/mL)", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$anc_plot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time/24, ANC)) +
      geom_ribbon(aes(ymin = pmin(ANC, 0.5), ymax = 0.5), fill = "#f4a6a6",
                  alpha = 0.6) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(0.5, 1.0, 1.5), linetype = c(1, 2, 3),
                 colour = c("#cc0000", "#999999", "#999999")) +
      annotate("text", x = 0, y = 0.55, hjust = 0, size = 3,
               label = "grade 4 threshold 0.5", colour = "#cc0000") +
      scale_y_log10() + labs(x = "days", y = "ANC (10^9/L), log scale") +
      theme_bw()
  })

  output$dsn_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `cumulative days ANC<0.5` = DSN_D,
                       `P(FN) cumulative` = P_FN) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL) + theme_bw()
  })

  output$chain_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, PROL, TR1, TR2, TR3, CIRC) %>% pivot_longer(-time)
    d$name <- factor(d$name, levels = c("PROL","TR1","TR2","TR3","CIRC"))
    ggplot(d, aes(time/24, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "days", y = "10^9/L", colour = NULL,
           subtitle = "the deficit is created in PROL and arrives in CIRC one transit time later") +
      theme_bw()
  })

  output$pool_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `storage (TR3)` = TR3, marginated = MARG,
                       circulating = CIRC, `stem reserve x5` = RESERVE_N) %>%
      mutate(`stem reserve x5` = `stem reserve x5` * 5) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "days", y = "10^9/L (reserve scaled x5)", colour = NULL) +
      theme_bw()
  })

  output$output_plot <- renderPlot({
    df <- sim()
    ## marrow OUTPUT proxy: the flux leaving the storage pool
    d <- df %>% mutate(`measured ANC` = ANC,
                       `marrow output (storage pool)` = TR3) %>%
      select(time, `measured ANC`, `marrow output (storage pool)`) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "days", y = "10^9/L", colour = NULL,
           subtitle = "on a growth factor the count can rise while the marrow is still empty") +
      theme_bw() + theme(legend.position = "top")
  })

  output$func_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `ANC (count)` = ANC,
                       `effective defence` = ANC_EFF,
                       `function (0-1)` = FUNC) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value, colour = name)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL, colour = NULL) +
      theme_bw() + theme(legend.position = "none")
  })

  output$infect_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `barrier integrity` = BAR, `bacteria` = BACT,
                       `CRP (mg/L)` = CRP, `temperature (C)` = TEMP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = NULL) + theme_bw()
  })

  output$lineage_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `ANC (10^9/L)` = ANC, `platelets (10^9/L)` = PLT,
                       `Hb (g/dL)` = HB) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL) + theme_bw()
  })

  output$lineage_tbl <- renderTable({
    df <- sim()
    data.frame(lineage = c("neutrophil", "platelet", "erythroid"),
               `transit time (h)` = c(89.3, 200, 150),
               `lifespan` = c("7.9 h", "10 d", "120 d"),
               `nadir` = c(round(min(df$ANC), 3), round(min(df$PLT), 0),
                           round(min(df$HB), 1)),
               `nadir day` = c(round(df$time[which.min(df$ANC)]/24, 1),
                               round(df$time[which.min(df$PLT)]/24, 1),
                               round(df$time[which.min(df$HB)]/24, 1)),
               check.names = FALSE)
  })

  output$rdi_plot <- renderPlot({
    df <- sim()
    d <- df %>% select(time, `tumour burden` = TUMOUR,
                       `day-1 gate open` = GATE_D1) %>% pivot_longer(-time)
    ggplot(d, aes(time/24, value)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "days", y = NULL) + theme_bw()
  })

  output$rdi_tbl <- renderTable({
    df <- sim()
    data.frame(
      quantity = c("delivered dose intensity (set by slider)",
                   "tumour burden at end",
                   "fraction of the course with the day-1 gate open",
                   "total days ANC < 0.5"),
      value = c(paste0(round(input$dose_mult * 100), "%"),
                round(tail(df$TUMOUR, 1), 2),
                round(mean(df$GATE_D1), 3),
                round(max(df$DSN_D), 2)))
  })

  cmp_sim <- reactive({
    opts <- list(
      "no support"             = list(gcsf="none", st=24, nd=1, tril=FALSE, abx=FALSE, dm=1.0),
      "pegfilgrastim d2"       = list(gcsf="pegfilgrastim", st=24, nd=1, tril=FALSE, abx=FALSE, dm=1.0),
      "pegfilgrastim same-day" = list(gcsf="pegfilgrastim", st=0.5, nd=1, tril=FALSE, abx=FALSE, dm=1.0),
      "filgrastim d2-11"       = list(gcsf="filgrastim", st=24, nd=10, tril=FALSE, abx=FALSE, dm=1.0),
      "filgrastim from d8"     = list(gcsf="filgrastim", st=168, nd=7, tril=FALSE, abx=FALSE, dm=1.0),
      "trilaciclib"            = list(gcsf="none", st=24, nd=1, tril=TRUE, abx=FALSE, dm=1.0),
      "levofloxacin"           = list(gcsf="none", st=24, nd=1, tril=FALSE, abx=TRUE, dm=1.0),
      "dose 80%"               = list(gcsf="none", st=24, nd=1, tril=FALSE, abx=FALSE, dm=0.8),
      "dose 60%"               = list(gcsf="none", st=24, nd=1, tril=FALSE, abx=FALSE, dm=0.6))
    sel <- input$cmp
    if (length(sel) == 0) return(NULL)
    do.call(rbind, lapply(sel, function(nm) {
      o <- opts[[nm]]
      d <- build_run(input$regimen, input$n_cycles, input$bsa, input$wt,
                     input$gfr, input$anc0, input$plt0, input$hb0,
                     input$reserve, input$sens, input$mttmult, input$gammult,
                     o$gcsf, o$st, o$nd, o$tril, input$dex, o$abx, o$dm)
      d$scenario <- nm
      d
    }))
  })

  output$cmp_plot <- renderPlot({
    d <- cmp_sim()
    if (is.null(d)) return(NULL)
    ggplot(d, aes(time/24, ANC, colour = scenario)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "#cc0000") +
      scale_y_log10() + labs(x = "days", y = "ANC (10^9/L), log scale",
                             colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$cmp_tbl <- renderTable({
    d <- cmp_sim()
    if (is.null(d)) return(NULL)
    d %>% group_by(scenario) %>%
      summarise(`nadir ANC` = round(min(ANC), 3),
                `nadir day` = round(time[which.min(ANC)]/24, 1),
                `DSN (d)` = round(max(DSN_D), 2),
                `P(FN)` = round(1 - exp(-max(CUMHAZ)), 3),
                `nadir plt` = round(min(PLT), 0),
                `final Hb` = round(last(HB), 1),
                `reserve` = round(last(RESERVE_N), 3),
                `tumour` = round(last(TUMOUR), 2),
                .groups = "drop")
  })
}

shinyApp(ui, server)
