# =====================================================================
# Autoimmune Pulmonary Alveolar Proteinosis (aPAP) — QSP Shiny dashboard
# ---------------------------------------------------------------------
# Ten tabs over the model in pap_mrgsolve_model.R:
#   1  Patient          - generate a patient from seroconversion
#   2  Mass balance     - the two large fluxes and the small difference
#   3  Antibody buffer  - the titration curve and where this patient sits
#   4  Drug delivery    - ELF pharmacology, molar ratio, coverage, reach
#   5  Gas exchange     - PaO2 / A-aDO2 / DLCO computed from shunt
#   6  Clinical         - DSS, SGRQ, 6MWD, oxygen requirement
#   7  Biomarkers       - KL-6, SP-D, CEA, LDH, CT density
#   8  Scenarios        - side-by-side comparison of regimens
#   9  Trial            - IMPALA-2 / PAGE replication in a virtual population
#  10  Covariates       - what predicts severity, and what does not
#
# Run with:  Rscript -e 'shiny::runApp("pap_shiny_app.R", port = 8080)'
# =====================================================================

library(shiny)
Sys.setenv(PAP_NORUN = "1")
source("pap_mrgsolve_model.R")

# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------
PAL <- c(drug = "#a63f7a", ctrl = "#5b6f85", pool = "#b8a12a",
         gas = "#3d8a4f", ab = "#7a4fa5", mac = "#c07a2a", warn = "#b03030")

lineplot <- function(d, cols, labs = cols, ylab = "", main = "", cols_pal = NULL,
                     xlab = "days", ylim = NULL, hline = NULL, log = "") {
  cols <- cols[cols %in% names(d)]
  if (!length(cols)) return(invisible(NULL))
  yy <- unlist(d[cols])
  if (is.null(ylim)) ylim <- range(c(yy, hline), na.rm = TRUE)
  cp <- if (is.null(cols_pal)) grDevices::hcl.colors(length(cols), "Dark 3") else cols_pal
  par(mar = c(4, 4.4, 3, 1))
  plot(d$time, d[[cols[1]]], type = "n", ylim = ylim, xlab = xlab, ylab = ylab,
       main = main, las = 1, log = log)
  grid(col = "#e8e8e8", lty = 1)
  if (!is.null(hline)) abline(h = hline, col = "#999999", lty = 2)
  for (i in seq_along(cols)) lines(d$time, d[[cols[i]]], col = cp[i], lwd = 2.2)
  if (length(cols) > 1)
    legend("topright", legend = labs[seq_along(cols)], col = cp, lwd = 2.2,
           bty = "n", cex = 0.85)
}

barcmp <- function(v, labs, ylab = "", main = "", col = PAL[["pool"]],
                   ref = NULL, reflab = "published") {
  par(mar = c(7, 4.4, 3, 1))
  bp <- barplot(v, names.arg = labs, las = 2, ylab = ylab, main = main,
                col = col, border = NA, cex.names = 0.8)
  if (!is.null(ref)) {
    points(bp, ref, pch = 18, cex = 1.6, col = PAL[["warn"]])
    legend("topleft", legend = reflab, pch = 18, col = PAL[["warn"]], bty = "n")
  }
  abline(h = 0, col = "#666666")
}

# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body {font-family: -apple-system, Helvetica, Arial, sans-serif;}
    .well {background:#fbfbfc; border-color:#e4e4e8;}
    h4 {color:#2c3e50; margin-top:2px;}
    .note {font-size:12px; color:#5a6570; line-height:1.5;}
    .warn {font-size:12px; color:#8a1a1a; line-height:1.5;}
  "))),
  titlePanel("Autoimmune Pulmonary Alveolar Proteinosis — QSP model explorer"),
  p(class = "note",
    HTML("A surfactant mass balance broken by loss of a <b>signal</b>, not loss of a cell. ",
         "Production is normal; the dominant catabolic sink is switched off by a ",
         "stoichiometric antibody buffer. A-aDO2 and DLCO are computed from shunt, the ",
         "alveolar gas equation and the Severinghaus curve - they are physics here, not scores.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("capfloor", "Residual GM-CSF-independent catabolic floor",
                  0.24, 0.42, 0.3095, step = 0.002),
      helpText(class = "note",
               "The covariate that sets severity. Below ~0.30 the patient reaches",
               "respiratory failure; above ~0.33 they never meet the enrolment",
               "criterion and are found on screening."),
      sliderInput("titre", "Serum GMAb titre (ug/mL)", 0, 300, 25, step = 1),
      sliderInput("fneut", "Neutralising fraction of that titre", 0.05, 0.85, 0.32,
                  step = 0.01),
      selectInput("class", "PAP class",
                  c("autoimmune" = "auto", "hereditary (no receptor)" = "hered",
                    "secondary (monocytopenia)" = "sec")),
      sliderInput("years", "Years since seroconversion at enrolment", 1, 25, 6,
                  step = 1),
      hr(),
      h4("Treatment"),
      selectInput("regimen", "Regimen",
                  c("none", "molgramostim 300 ug QD",
                    "molgramostim 300 ug QD every other week",
                    "sargramostim 125 ug BID 7/14 (PAGE)",
                    "sargramostim SC 5 ug/kg/d", "sargramostim SC 20 ug/kg/d",
                    "atorvastatin 80 mg", "molgramostim + statin",
                    "rituximab 1000 mg x2", "plasmapheresis x10",
                    "FcRn inhibitor")),
      checkboxInput("wll", "Whole lung lavage at day 30", FALSE),
      sliderInput("dur", "Follow-up (days)", 60, 1460, 336, step = 28),
      hr(),
      h4("Delivery and reach"),
      sliderInput("fdep", "Alveolar deposition fraction", 0.05, 0.9, 0.40,
                  step = 0.05),
      sliderInput("edgef", "Interface clearance coefficient", 0, 0.15, 0.020,
                  step = 0.005),
      helpText(class = "note",
               "Aerosol follows ventilation, and the burden sits in units that do",
               "not ventilate. Reach, not dose, is what limits inhaled GM-CSF."),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient",
                 fluidRow(column(6, plotOutput("p1a", height = 300)),
                          column(6, plotOutput("p1b", height = 300))),
                 fluidRow(column(6, plotOutput("p1c", height = 280)),
                          column(6, verbatimTextOutput("p1txt")))),
        tabPanel("2 Mass balance",
                 fluidRow(column(7, plotOutput("p2a", height = 320)),
                          column(5, plotOutput("p2b", height = 320))),
                 fluidRow(column(7, plotOutput("p2c", height = 280)),
                          column(5, htmlOutput("p2txt")))),
        tabPanel("3 Antibody buffer",
                 fluidRow(column(6, plotOutput("p3a", height = 330)),
                          column(6, plotOutput("p3b", height = 330))),
                 fluidRow(column(12, htmlOutput("p3txt")))),
        tabPanel("4 Drug delivery",
                 fluidRow(column(6, plotOutput("p4a", height = 300)),
                          column(6, plotOutput("p4b", height = 300))),
                 fluidRow(column(6, plotOutput("p4c", height = 300)),
                          column(6, plotOutput("p4d", height = 300)))),
        tabPanel("5 Gas exchange",
                 fluidRow(column(6, plotOutput("p5a", height = 300)),
                          column(6, plotOutput("p5b", height = 300))),
                 fluidRow(column(6, plotOutput("p5c", height = 300)),
                          column(6, plotOutput("p5d", height = 300)))),
        tabPanel("6 Clinical",
                 fluidRow(column(6, plotOutput("p6a", height = 300)),
                          column(6, plotOutput("p6b", height = 300))),
                 fluidRow(column(6, plotOutput("p6c", height = 300)),
                          column(6, plotOutput("p6d", height = 300)))),
        tabPanel("7 Biomarkers",
                 fluidRow(column(6, plotOutput("p7a", height = 300)),
                          column(6, plotOutput("p7b", height = 300))),
                 fluidRow(column(6, plotOutput("p7c", height = 300)),
                          column(6, plotOutput("p7d", height = 300)))),
        tabPanel("8 Scenarios",
                 fluidRow(column(12, plotOutput("p8a", height = 340))),
                 fluidRow(column(6, plotOutput("p8b", height = 320)),
                          column(6, plotOutput("p8c", height = 320)))),
        tabPanel("9 Trial",
                 fluidRow(column(4, numericInput("npop", "Virtual patients per arm",
                                                 16, 6, 40, 2)),
                          column(8, htmlOutput("p9note"))),
                 fluidRow(column(6, plotOutput("p9a", height = 330)),
                          column(6, plotOutput("p9b", height = 330))),
                 fluidRow(column(12, verbatimTextOutput("p9txt")))),
        tabPanel("10 Covariates",
                 fluidRow(column(6, plotOutput("p10a", height = 320)),
                          column(6, plotOutput("p10b", height = 320))),
                 fluidRow(column(12, htmlOutput("p10txt"))))
      )
    )
  )
)

# ---------------------------------------------------------------------
# server
# ---------------------------------------------------------------------
server <- function(input, output, session) {

  basepars <- reactive({
    p <- list(CAPFLOOR = input$capfloor, TITRE0 = input$titre,
              FNEUT = input$fneut, FDEP = input$fdep, EDGEF = input$edgef)
    if (input$class == "hered") p$RECFRAC <- 0
    if (input$class == "sec") { p$ABON <- 0; p$MONOF <- 0.34 }
    p
  })

  # natural history to the enrolment point
  patient <- eventReactive(input$go, {
    p <- basepars()
    nh <- natural_history(p, years = max(input$years, 2), delta = 5)
    i <- which.min(abs(nh$time - input$years * 365))
    st <- as.list(nh[i, CMTN]); names(st) <- CMTN
    list(nh = nh, state = st, pars = p)
  }, ignoreNULL = FALSE)

  dosing <- reactive({
    d <- input$dur
    switch(input$regimen,
      "none" = NULL,
      "molgramostim 300 ug QD" = ev_inh(300, d),
      "molgramostim 300 ug QD every other week" = ev_inh_cycle(300, d, 1, 7, 7),
      "sargramostim 125 ug BID 7/14 (PAGE)" = ev_inh_cycle(125, d, 2, 7, 7),
      "sargramostim SC 5 ug/kg/d" = ev_sc(5, 70, d),
      "sargramostim SC 20 ug/kg/d" = ev_sc(20, 70, d),
      "atorvastatin 80 mg" = ev_statin(80, d),
      "molgramostim + statin" = ev_bind(ev_inh(300, d), ev_statin(80, d)),
      "rituximab 1000 mg x2" = ev_rtx(),
      "plasmapheresis x10" = NULL,
      "FcRn inhibitor" = NULL)
  })

  extrapars <- reactive({
    p <- list()
    if (input$wll) p$WLLT1 <- 30
    if (input$regimen == "plasmapheresis x10") { p$PLEXT <- 30; p$PLEXN <- 10 }
    if (input$regimen == "FcRn inhibitor") { p$FCRNT <- 7; p$FCRNX <- 3 }
    p
  })

  run <- reactive({
    pt <- patient()
    pr <- c(pt$pars, extrapars())
    dl <- if (input$dur > 400) 2 else 1
    list(tx = sim(pr, end = input$dur, delta = dl, ev = dosing(),
                  init = pt$state),
         ctl = sim(c(pt$pars, if (input$wll) list(WLLT1 = 30) else NULL),
                   end = input$dur, delta = dl, init = pt$state),
         pt = pt)
  })

  # ---- 1 Patient ----------------------------------------------------
  output$p1a <- renderPlot({
    nh <- patient()$nh; nh$yr <- nh$time / 365
    d <- nh; d$time <- d$yr
    lineplot(d, c("POOLR"), "burden / healthy pool", "x normal",
             "Burden accumulates as the integral of a small imbalance",
             cols_pal = PAL[["pool"]], xlab = "years since seroconversion")
    abline(v = input$years, col = PAL[["warn"]], lty = 3, lwd = 2)
  })
  output$p1b <- renderPlot({
    nh <- patient()$nh; d <- nh; d$time <- nh$time / 365
    lineplot(d, c("AADO2", "DLCO"), c("A-aDO2 (mmHg)", "DLCO (%pred)"), "",
             "Gas exchange over the natural history",
             cols_pal = c(PAL[["gas"]], PAL[["ab"]]),
             xlab = "years since seroconversion", hline = 30)
    abline(v = input$years, col = PAL[["warn"]], lty = 3, lwd = 2)
  })
  output$p1c <- renderPlot({
    nh <- patient()$nh; d <- nh; d$time <- nh$time / 365
    lineplot(d, c("TITRE", "SIGREL", "TOTCAP"),
             c("titre (ug/mL)", "signalling (rel)", "catabolic capacity (rel)"),
             "", "Antibody rises first, then capacity falls, then the pool fills",
             cols_pal = c(PAL[["ab"]], PAL[["gas"]], PAL[["mac"]]),
             xlab = "years", log = "")
  })
  output$p1txt <- renderText({
    r <- run(); d <- r$tx; nh <- r$pt$nh
    pres <- which(nh$PAO2A <= 70)
    sprintf(paste0(
      "PATIENT AT ENROLMENT (year %.0f after seroconversion)\n",
      "-------------------------------------------------\n",
      "burden            %8.0f mg  (%.1f x healthy pool)\n",
      "of which sequestered %5.0f mg  (%.0f%%)\n",
      "fully filled units   %7.1f %%   partially filled %5.1f %%\n",
      "PaO2              %8.1f mmHg   A-aDO2 %6.1f mmHg\n",
      "DLCO              %8.1f %%pred  DSS    %6.0f\n",
      "serum GMAb        %8.1f ug/mL  ELF GMAb %5.3f ug/mL\n",
      "free ELF GM-CSF   %8.3f pM     (healthy %.2f pM)\n",
      "signalling        %8.3f of healthy\n",
      "catabolic capacity%8.3f of healthy\n",
      "net imbalance     %+8.1f mg/d\n",
      "infection hazard  %8.3f per year\n",
      "%s"),
      input$years, d$PTOTo[1], d$POOLR[1], d$SEQ[1], 100 * d$SEQ[1] / d$PTOTo[1],
      100 * d$FFUL[1], 100 * d$FPART[1], d$PAO2A[1], d$AADO2[1], d$DLCO[1],
      d$DSS[1], d$TITRE[1], d$TITELF[1], d$LFREE[1],
      as.list(param(mod))$CGMH * 1000 / as.list(param(mod))$MWGM,
      d$SIGREL[1], d$TOTCAP[1], -d$NETBAL[1], d$INFYR[1],
      if (length(pres))
        sprintf("met the trial enrolment criterion (PaO2 <= 70) at year %.1f\n",
                nh$time[pres[1]] / 365)
      else "NEVER meets the enrolment criterion: screening-detected disease\n")
  })

  # ---- 2 Mass balance ----------------------------------------------
  output$p2a <- renderPlot({
    d <- run()$tx
    lineplot(d, c("PROD", "NETMAC", "IIDEG", "MUCO"),
             c("de novo production", "net macrophage", "type-II", "mucociliary"),
             "mg/d", "Two large fluxes; the disease is their small difference",
             cols_pal = c(PAL[["pool"]], PAL[["mac"]], "#4a6f9a", "#7f8c8d"))
  })
  output$p2b <- renderPlot({
    d <- run()$tx
    barcmp(c(d$NETMAC[1], d$IIDEG[1], d$MUCO[1], d$PROD[1]),
           c("macrophage", "type-II", "mucociliary", "production"),
           "mg/d", "Sinks vs source at enrolment",
           col = c(PAL[["mac"]], "#4a6f9a", "#7f8c8d", PAL[["pool"]]))
  })
  output$p2c <- renderPlot({
    r <- run()
    d <- r$tx; c0 <- r$ctl
    par(mar = c(4, 4.4, 3, 1))
    ylim <- range(c(d$PTOTo, c0$PTOTo, d$SEQ))
    plot(d$time, d$PTOTo, type = "n", ylim = ylim, las = 1, xlab = "days",
         ylab = "mg", main = "Alveolar burden: treated vs untreated")
    grid(col = "#e8e8e8", lty = 1)
    lines(c0$time, c0$PTOTo, col = PAL[["ctrl"]], lwd = 2.4)
    lines(d$time, d$PTOTo, col = PAL[["drug"]], lwd = 2.4)
    lines(d$time, d$SEQ, col = PAL[["pool"]], lwd = 1.8, lty = 2)
    legend("topright", c("untreated", "treated", "sequestered (treated)"),
           col = c(PAL[["ctrl"]], PAL[["drug"]], PAL[["pool"]]),
           lwd = 2.2, lty = c(1, 1, 2), bty = "n", cex = 0.85)
  })
  output$p2txt <- renderUI({
    d <- run()$tx
    HTML(sprintf(paste0(
      "<p class='note'><b>Why the pool reaches tens of times normal.</b> ",
      "Type-II re-uptake and re-secretion is a large flux (%.0f mg/d here) with ",
      "zero net effect at steady state. Lump it into clearance and the model ",
      "predicts a pool that settles at 2-3x normal; patients reach 30-100x. ",
      "Only true catabolism and mucociliary egress are sinks.</p>",
      "<p class='note'><b>Mass-balance audit.</b> residual %+.2e mg against ",
      "%.0f mg turned over in this run - the phospholipid balance is closed and ",
      "checked, not assumed.</p>",
      "<p class='note'>Digestion is capacity-limited, not substrate-limited. ",
      "That is what makes the foamy macrophage a <i>stalled</i> macrophage: ",
      "written first-order in lipid, a cell carrying six times the load would ",
      "digest as much as a healthy one and the disease could not be generated.</p>"),
      3 * d$PROD[1], fin(d, "MBAL"), fin(d, "CUMPROD")))
  })

  # ---- 3 Antibody buffer -------------------------------------------
  output$p3a <- renderPlot({
    p <- basepars()
    tt <- c(0.3, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 25, 40, 60, 100, 200, 300)
    sg <- vapply(tt, function(x) {
      d <- sim(c(p, list(TITRE0 = x, ABON = 1)), end = 400, delta = 40)
      fin(d, "SIGREL")
    }, numeric(1))
    par(mar = c(4, 4.4, 3, 1))
    plot(tt, sg, type = "b", log = "x", pch = 19, lwd = 2.2, las = 1,
         col = PAL[["ab"]], ylim = c(0, 1.05),
         xlab = "serum GMAb titre (ug/mL, log scale)",
         ylab = "GM-CSF signalling / healthy",
         main = "The buffer titration curve")
    grid(col = "#e8e8e8", lty = 1)
    abline(v = 5, col = PAL[["warn"]], lty = 2, lwd = 2)
    abline(h = 0.5, col = "#999999", lty = 3)
    points(input$titre, approx(tt, sg, xout = input$titre, rule = 2)$y,
           pch = 21, bg = "white", col = PAL[["warn"]], cex = 2.2, lwd = 2.5)
    text(5, 0.95, "published critical\nthreshold 5 ug/mL", pos = 4, cex = 0.8,
         col = PAL[["warn"]])
  })
  output$p3b <- renderPlot({
    d <- run()$tx
    lineplot(d, c("TITRE", "TITELF"),
             c("serum (ug/mL)", "ELF (ug/mL)"), "ug/mL",
             "Antibody in the two compartments",
             cols_pal = c(PAL[["ab"]], "#c9a8e8"), hline = 5)
  })
  output$p3txt <- renderUI({
    d <- run()$tx
    HTML(sprintf(paste0(
      "<p class='note'><b>The antibody is a stoichiometric buffer, not an ",
      "inhibitor with an IC50.</b> Alveolar GM-CSF sits at tens of pM; ",
      "neutralising sites at hundreds to thousands of pM. Free ligand is solved ",
      "from the 1:1 binding equilibrium - in the numerically stable form, ",
      "because at these ratios the textbook root cancels to zero in double ",
      "precision. This patient: %.0f ug/mL serum, %.3f ug/mL in lining fluid, ",
      "%.1f pM of neutralising sites, free GM-CSF %.4f pM against a healthy ",
      "%.2f pM.</p>",
      "<p class='note'>Two clinical facts come out of the one equation: a steep ",
      "titre <b>threshold</b> near 5 ug/mL, and <b>no titre-severity ",
      "correlation above it</b> - past the threshold free GM-CSF is already ",
      "near zero and more antibody changes nothing.</p>",
      "<p class='warn'>The model only partly reproduces the observed absence of ",
      "a titre-severity correlation (see D05). Reported titre measures binding, ",
      "not neutralisation; varying the neutralising fraction between patients ",
      "weakens the correlation but does not abolish it.</p>"),
      d$TITRE[1], d$TITELF[1],
      d$TITELF[1] * 30 * (2e6 / 150000) * input$fneut / 0.03 / 1000,
      d$LFREE[1], as.list(param(mod))$CGMH * 1000 / as.list(param(mod))$MWGM))
  })

  # ---- 4 Drug delivery ---------------------------------------------
  output$p4a <- renderPlot({
    pt <- patient()
    d <- sim(c(pt$pars, extrapars()), end = 8, delta = 1 / 96, ev = dosing(),
             init = pt$state)
    lineplot(d, "LFREE", "free GM-CSF", "pM (log)",
             "Free GM-CSF in lining fluid, first 8 days",
             cols_pal = PAL[["drug"]], log = "y")
  })
  output$p4b <- renderPlot({
    pt <- patient()
    d <- sim(c(pt$pars, extrapars()), end = 8, delta = 1 / 96, ev = dosing(),
             init = pt$state)
    lineplot(d, "SIGREL", "signalling", "relative to healthy",
             "Time above the signalling threshold is the efficacy variable",
             cols_pal = PAL[["gas"]], hline = c(0.5, 1))
  })
  output$p4c <- renderPlot({
    pt <- patient()
    doses <- c(10, 30, 75, 150, 300, 600, 1200, 3000)
    v <- vapply(doses, function(dd) {
      d <- sim(c(pt$pars, extrapars()), end = 168, ev = ev_inh(dd, 168),
               init = pt$state)
      at(d, 168, "DLCO") - d$DLCO[1]
    }, numeric(1))
    par(mar = c(4, 4.4, 3, 1))
    plot(doses, v, type = "b", pch = 19, log = "x", lwd = 2.2, las = 1,
         col = PAL[["drug"]], xlab = "nominal inhaled dose (ug/d, log)",
         ylab = "dDLCO at 24 wk (%pred)",
         main = "Dose-response: threshold, then plateau")
    grid(col = "#e8e8e8", lty = 1)
    abline(v = 300, col = "#999999", lty = 2)
  })
  output$p4d <- renderPlot({
    d <- run()$tx
    lineplot(d, c("REACH", "WFILL", "TOTCAP"),
             c("drug reach", "burden in filled units", "catabolic capacity"),
             "fraction / relative",
             "Reach is what limits inhaled therapy",
             cols_pal = c(PAL[["drug"]], PAL[["ctrl"]], PAL[["mac"]]))
  })

  # ---- 5 Gas exchange ----------------------------------------------
  output$p5a <- renderPlot({
    r <- run()
    par(mar = c(4, 4.4, 3, 1))
    plot(r$tx$time, r$tx$AADO2, type = "n", las = 1, xlab = "days",
         ylab = "A-aDO2 (mmHg)", main = "A-aDO2: treated vs untreated",
         ylim = range(c(r$tx$AADO2, r$ctl$AADO2)))
    grid(col = "#e8e8e8", lty = 1)
    lines(r$ctl$time, r$ctl$AADO2, col = PAL[["ctrl"]], lwd = 2.4)
    lines(r$tx$time, r$tx$AADO2, col = PAL[["drug"]], lwd = 2.4)
    legend("topright", c("untreated", "treated"),
           col = c(PAL[["ctrl"]], PAL[["drug"]]), lwd = 2.4, bty = "n")
  })
  output$p5b <- renderPlot({
    r <- run()
    par(mar = c(4, 4.4, 3, 1))
    plot(r$tx$time, r$tx$DLCO, type = "n", las = 1, xlab = "days",
         ylab = "DLCO (% predicted)", main = "DLCO: treated vs untreated",
         ylim = range(c(r$tx$DLCO, r$ctl$DLCO)))
    grid(col = "#e8e8e8", lty = 1)
    lines(r$ctl$time, r$ctl$DLCO, col = PAL[["ctrl"]], lwd = 2.4)
    lines(r$tx$time, r$tx$DLCO, col = PAL[["drug"]], lwd = 2.4)
    legend("bottomright", c("untreated", "treated"),
           col = c(PAL[["ctrl"]], PAL[["drug"]]), lwd = 2.4, bty = "n")
  })
  output$p5c <- renderPlot({
    d <- run()$tx
    lineplot(d, c("PAO2A", "PAO2AX"), c("rest", "exercise"), "PaO2 (mmHg)",
             "Exercise amplifies the same shunt",
             cols_pal = c(PAL[["gas"]], PAL[["warn"]]), hline = c(60, 70))
  })
  output$p5d <- renderPlot({
    d <- run()$tx
    lineplot(d, c("SHUNTF", "FFUL", "FPART"),
             c("shunt fraction", "fully filled", "partially filled"), "fraction",
             "Shunt comes from fully filled units only",
             cols_pal = c(PAL[["warn"]], PAL[["ctrl"]], "#a9b8c9"))
  })

  # ---- 6 Clinical ---------------------------------------------------
  output$p6a <- renderPlot({
    d <- run()$tx
    lineplot(d, "DSS", "DSS", "score (1-5)",
             "Disease severity score (Inoue 2008 definition)",
             cols_pal = PAL[["warn"]], ylim = c(0.8, 5.2))
  })
  output$p6b <- renderPlot({
    d <- run()$tx
    lineplot(d, "SGRQ", "SGRQ-T", "points",
             "St George Respiratory Questionnaire (total)",
             cols_pal = PAL[["ab"]])
  })
  output$p6c <- renderPlot({
    d <- run()$tx
    lineplot(d, "SIXMWD", "6MWD", "metres", "Exercise capacity",
             cols_pal = PAL[["gas"]])
  })
  output$p6d <- renderPlot({
    d <- run()$tx
    lineplot(d, c("O2REQ", "MMRC", "INFYR"),
             c("oxygen requirement index", "mMRC", "infection hazard /y"), "",
             "Oxygen need, dyspnoea, infection risk",
             cols_pal = c(PAL[["warn"]], PAL[["ctrl"]], PAL[["mac"]]))
  })

  # ---- 7 Biomarkers -------------------------------------------------
  output$p7a <- renderPlot({
    d <- run()$tx
    lineplot(d, "KL6", "KL-6", "U/mL", "Serum KL-6", cols_pal = PAL[["pool"]],
             hline = 500)
  })
  output$p7b <- renderPlot({
    d <- run()$tx
    lineplot(d, c("SPD", "CEA"), c("SP-D (ng/mL)", "CEA (ng/mL)"), "",
             "Serum SP-D and CEA", cols_pal = c(PAL[["gas"]], PAL[["ab"]]))
  })
  output$p7c <- renderPlot({
    d <- run()$tx
    lineplot(d, "CTHU", "mean lung density", "HU",
             "Quantitative CT density (endpoint in PAGE)",
             cols_pal = PAL[["ctrl"]], hline = -850)
  })
  output$p7d <- renderPlot({
    d <- run()$tx
    lineplot(d, c("NEUT", "EOS", "HB"),
             c("neutrophils (10^9/L)", "eosinophils", "haemoglobin (g/dL)"), "",
             "Systemic effects: myeloid response and erythrocytosis",
             cols_pal = c(PAL[["mac"]], PAL[["drug"]], PAL[["warn"]]))
  })

  # ---- 8 Scenarios --------------------------------------------------
  scen_cmp <- reactive({
    pt <- patient(); st <- pt$state; pr <- pt$pars
    L <- list(
      "untreated"            = list(p = pr, e = NULL),
      "molgra 300 QD"        = list(p = pr, e = ev_inh(300, 336)),
      "molgra 300 alt week"  = list(p = pr, e = ev_inh_cycle(300, 336, 1, 7, 7)),
      "sargra 125 BID 7/14"  = list(p = pr, e = ev_inh_cycle(125, 336, 2, 7, 7)),
      "SC 5 ug/kg/d"         = list(p = pr, e = ev_sc(5, 70, 336)),
      "SC 20 ug/kg/d"        = list(p = pr, e = ev_sc(20, 70, 336)),
      "statin 80 mg"         = list(p = pr, e = ev_statin(80, 336)),
      "WLL alone"            = list(p = c(pr, list(WLLT1 = 30)), e = NULL),
      "WLL + molgra"         = list(p = c(pr, list(WLLT1 = 30)),
                                    e = ev_inh(300, 336)),
      "rituximab"            = list(p = pr, e = ev_rtx()),
      "plasmapheresis"       = list(p = c(pr, list(PLEXT = 30, PLEXN = 10)),
                                    e = NULL),
      "FcRn inhibitor"       = list(p = c(pr, list(FCRNT = 7, FCRNX = 3)),
                                    e = NULL))
    do.call(rbind, lapply(names(L), function(nm) {
      d <- sim(L[[nm]]$p, end = 336, ev = L[[nm]]$e, init = st)
      data.frame(name = nm, dDLCO = at(d, 168, "DLCO") - d$DLCO[1],
                 dAADO2 = at(d, 168, "AADO2") - d$AADO2[1],
                 dpool = 100 * (at(d, 336, "PTOTo") / d$PTOTo[1] - 1))
    }))
  })
  output$p8a <- renderPlot({
    s <- scen_cmp()
    barcmp(s$dDLCO, s$name, "dDLCO at 24 wk (%pred)",
           "Regimens compared in the same patient", col = PAL[["drug"]])
  })
  output$p8b <- renderPlot({
    s <- scen_cmp()
    barcmp(s$dAADO2, s$name, "dA-aDO2 at 24 wk (mmHg)",
           "A-aDO2 change (negative is better)", col = PAL[["gas"]])
  })
  output$p8c <- renderPlot({
    s <- scen_cmp()
    barcmp(s$dpool, s$name, "burden change at 48 wk (%)",
           "Burden change", col = PAL[["pool"]])
  })

  # ---- 9 Trial ------------------------------------------------------
  trial <- eventReactive(input$go, {
    pop <- vpop(n = input$npop, pars = list(FDEP = input$fdep,
                                            EDGEF = input$edgef))
    list(pop = pop,
         drug = trial_arm(pop, ev_inh(300, 336), 336),
         pbo  = trial_arm(pop, NULL, 336))
  }, ignoreNULL = FALSE)

  output$p9note <- renderUI({
    HTML(paste0("<p class='note'>The placebo arm here is not a placebo ",
      "parameter. It is <b>enrolment at the low point of a deterioration</b> ",
      "that then resolves, plus the DLCO measurement-learning effect. PAGE ran ",
      "a 12-week observation period and excluded improvers - and PAGE is the ",
      "trial whose placebo arm did not move.</p>"))
  })
  output$p9a <- renderPlot({
    t <- trial()
    barcmp(c(mean(t$drug$dDLCO24), mean(t$pbo$dDLCO24),
             mean(t$drug$dDLCO48), mean(t$pbo$dDLCO48)),
           c("drug 24wk", "placebo 24wk", "drug 48wk", "placebo 48wk"),
           "dDLCO (%pred)", "IMPALA-2 replication",
           col = c(PAL[["drug"]], PAL[["ctrl"]], PAL[["drug"]], PAL[["ctrl"]]),
           ref = c(9.8, 3.8, 11.6, 4.7), reflab = "IMPALA-2 observed")
  })
  output$p9b <- renderPlot({
    t <- trial()
    par(mar = c(4, 4.4, 3, 1))
    plot(t$drug$dlco0, t$drug$dDLCO24, pch = 19, col = PAL[["drug"]], las = 1,
         xlab = "baseline DLCO (%pred)", ylab = "dDLCO at 24 wk",
         main = "Who responds")
    grid(col = "#e8e8e8", lty = 1)
    points(t$pbo$dlco0, t$pbo$dDLCO24, pch = 1, col = PAL[["ctrl"]])
    legend("topright", c("molgramostim", "placebo"), pch = c(19, 1),
           col = c(PAL[["drug"]], PAL[["ctrl"]]), bty = "n")
  })
  output$p9txt <- renderText({
    t <- trial()
    sprintf(paste0(
      "VIRTUAL TRIAL (n = %d screened in per arm)\n",
      "baseline DLCO %.1f %%pred, A-aDO2 %.1f mmHg   (observed ~45 and ~39)\n",
      "-----------------------------------------------------------------\n",
      "dDLCO   24 wk : drug %+6.2f   placebo %+6.2f   difference %+6.2f  (obs +6.0)\n",
      "dDLCO   48 wk : drug %+6.2f   placebo %+6.2f   difference %+6.2f  (obs +6.9)\n",
      "dA-aDO2 24 wk : drug %+6.2f   placebo %+6.2f   difference %+6.2f  (obs -6.2)\n",
      "dSGRQ   24 wk : drug %+6.2f   placebo %+6.2f   difference %+6.2f  (obs -6.6)\n",
      "-----------------------------------------------------------------\n",
      "The 24-week differences land on two trials and two endpoints. The\n",
      "48-week difference is OVER-PREDICTED, because the model contains a\n",
      "recovery-side positive feedback: clearing the lung opens units, which\n",
      "lets more aerosol reach the rest. Quote the 24-week number.\n"),
      nrow(t$drug), mean(t$drug$dlco0), mean(t$drug$aado20),
      mean(t$drug$dDLCO24), mean(t$pbo$dDLCO24),
      mean(t$drug$dDLCO24) - mean(t$pbo$dDLCO24),
      mean(t$drug$dDLCO48), mean(t$pbo$dDLCO48),
      mean(t$drug$dDLCO48) - mean(t$pbo$dDLCO48),
      mean(t$drug$dAADO24), mean(t$pbo$dAADO24),
      mean(t$drug$dAADO24) - mean(t$pbo$dAADO24),
      mean(t$drug$dSGRQ24), mean(t$pbo$dSGRQ24),
      mean(t$drug$dSGRQ24) - mean(t$pbo$dSGRQ24))
  })

  # ---- 10 Covariates ------------------------------------------------
  output$p10a <- renderPlot({
    p <- basepars()
    cf <- seq(0.27, 0.37, by = 0.01)
    v <- vapply(cf, function(cc) {
      d <- natural_history(c(p, list(CAPFLOOR = cc)), years = 20, delta = 40)
      c(fin(d, "AADO2"), fin(d, "DLCO"))
    }, numeric(2))
    par(mar = c(4, 4.4, 3, 4.4))
    plot(cf, v[1, ], type = "b", pch = 19, lwd = 2.2, col = PAL[["gas"]], las = 1,
         xlab = "residual catabolic floor", ylab = "A-aDO2 at 20 y (mmHg)",
         main = "Severity vs catabolic floor")
    grid(col = "#e8e8e8", lty = 1)
    abline(v = input$capfloor, col = PAL[["warn"]], lty = 3, lwd = 2)
  })
  output$p10b <- renderPlot({
    p <- basepars()
    tt <- c(2, 5, 10, 20, 40, 80, 160, 300)
    v <- vapply(tt, function(x) {
      d <- natural_history(c(p, list(TITRE0 = x)), years = 20, delta = 40)
      fin(d, "AADO2")
    }, numeric(1))
    par(mar = c(4, 4.4, 3, 1))
    plot(tt, v, type = "b", pch = 19, log = "x", lwd = 2.2, col = PAL[["ab"]],
         las = 1, xlab = "serum GMAb titre (ug/mL, log)",
         ylab = "A-aDO2 at 20 y (mmHg)", main = "Severity vs titre: saturates")
    grid(col = "#e8e8e8", lty = 1)
    abline(v = 5, col = PAL[["warn"]], lty = 2)
  })
  output$p10txt <- renderUI({
    HTML(paste0(
      "<p class='note'><b>What predicts severity.</b> Above the critical titre ",
      "free GM-CSF is already near zero, so raising the titre further changes ",
      "little - the curve on the right flattens. The curve on the left does ",
      "not: a change of one part in ten in the residual GM-CSF-independent ",
      "catabolic floor spans asymptomatic disease to respiratory failure. That ",
      "steepness is this model's account of why aPAP severity is so ",
      "heterogeneous and why the 31.8% found on health screening exist at ",
      "all.</p>",
      "<p class='warn'><b>Two honest failures.</b> (1) The model does NOT ",
      "produce bistability: for a given floor there is one equilibrium burden, ",
      "reached from either direction. The hypothesis the map was drawn on is ",
      "refuted. (2) The titre-severity correlation is weakened but not ",
      "abolished, so buffer saturation is not a complete account of Inoue's ",
      "negative finding.</p>"))
  })
}

shinyApp(ui, server)
