##  Infective Endocarditis QSP — Shiny Dashboard
##  ============================================================================
##  Interactive front end for ie_mrgsolve_model.R.
##
##  Run with:
##      shiny::runApp("ie_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp("ie_shiny_app.R", launch.browser = TRUE)'
##
##  The dashboard is organised around the model's argument rather than around
##  its compartment list. The argument is that infective endocarditis is a
##  GEOMETRY problem wearing the clothes of a susceptibility problem, so:
##
##    Tab 1  The ledger        — logs required vs logs delivered, per shell.
##                               If you read one tab, read this one.
##    Tab 2  Patient & lesion  — the presenting state is EARNED by simulating
##                               the natural history, not typed in.
##    Tab 3  Pharmacokinetics  — plasma, and then the same drug at three depths.
##    Tab 4  The exponential   — the depth profile itself, and what changing
##                               lambda does to it. This is the whole model in
##                               one plot.
##    Tab 5  Bacteriology      — shells, phenotypes, and the rpoB subpopulation.
##    Tab 6  Embolism & valve  — why the hazard falls before the lesion shrinks.
##    Tab 7  Toxicity          — the AKI feedback loop, CK, ALT.
##    Tab 8  Regimen compare   — arbitrary head-to-head, including the negatives.
##    Tab 9  Exposure sweeps   — vancomycin AUC, penetration length, duration.
##
##  Everything is computed live from the ODE model; nothing here is a stored
##  result.
##  ============================================================================

library(shiny)
library(mrgsolve)

source("ie_mrgsolve_model.R")
MOD <- ie_model()

PAL <- c(base = "#2c6fbb", drug = "#c0392b", alt  = "#27ae60",
         warn = "#e67e22", grey = "#7f8c8d", purp = "#8e44ad",
         core = "#8b0000", mid  = "#d35400", surf = "#16a085")

PHENOS <- names(IE_PHENO)
RXS    <- c("none", "nafcillin", "cefazolin", "ampicillin", "ceftriaxone",
            "vanco", "dapto6", "dapto8", "dapto10", "gent", "rif", "rif0",
            "dalba", "oral", "bl_continuous")

## ---------------------------------------------------------------------------
## small plotting helpers (base graphics — no extra dependencies)
## ---------------------------------------------------------------------------
gpar <- function(...) par(mar = c(4.2, 4.4, 2.6, 1.2), mgp = c(2.6, 0.7, 0),
                         cex.lab = 1.02, cex.axis = 0.92, bty = "l", ...)

lineplot <- function(x, ys, cols, labs, xlab = "day", ylab = "", main = "",
                     lty = NULL, ylim = NULL, hline = NULL, hlab = NULL,
                     log = "") {
  gpar()
  if (is.null(ylim)) ylim <- range(unlist(ys), na.rm = TRUE)
  if (is.null(lty))  lty  <- rep(1, length(ys))
  plot(NA, xlim = range(x), ylim = ylim, xlab = xlab, ylab = ylab,
       main = main, log = log)
  grid(col = "#e8e8e8", lty = 1)
  if (!is.null(hline))
    for (i in seq_along(hline)) abline(h = hline[i], col = "#999999", lty = 3)
  if (!is.null(hlab))
    for (i in seq_along(hline))
      text(max(x), hline[i], hlab[i], pos = 3, cex = 0.75, col = "#777777")
  for (i in seq_along(ys)) lines(x, ys[[i]], col = cols[i], lwd = 2.2, lty = lty[i])
  legend("topright", legend = labs, col = cols, lwd = 2.2, lty = lty,
         bty = "n", cex = 0.85)
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Infective Endocarditis QSP — geometry, growth-rate tolerance, and the log-kill ledger"),
  tags$p(style = "color:#666;margin-top:-8px",
         paste("The MIC is measured on 10^5 CFU/mL of log-phase planktonic cells in broth.",
               "This lesion is 10^9-10^10 stationary-phase cells behind 2-5 mm of avascular",
               "platelet-fibrin. Everything below follows from that mismatch.")),
  hr(),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient & lesion"),
      selectInput("pheno", "Organism / valve", choices = PHENOS,
                  selected = "MSSA_native"),
      sliderInput("dxday", "Lesion age at diagnosis (days)",
                  min = 3, max = 35, value = 21, step = 1),
      helpText(style = "font-size:11px",
               "The presenting vegetation is not typed in. The model grows it",
               "from a 10^3 CFU seed for this many days, and the size, the",
               "dormant fraction and the pre-existing rpoB mutant count are",
               "all consequences of that."),

      h4("Regimen A"),
      selectizeInput("rxA", NULL, choices = RXS, selected = "nafcillin",
                     multiple = TRUE),
      selectInput("blA", "beta-lactam PK slot",
                  choices = c("(from phenotype)", names(IE_BLPK)),
                  selected = "(from phenotype)"),

      h4("Regimen B (comparator)"),
      selectizeInput("rxB", NULL, choices = RXS, selected = "vanco",
                     multiple = TRUE),
      selectInput("blB", "beta-lactam PK slot",
                  choices = c("(from phenotype)", names(IE_BLPK)),
                  selected = "(from phenotype)"),

      h4("Course"),
      sliderInput("days", "Duration of therapy (days)", 7, 84, 42, step = 7),
      sliderInput("follow", "Follow-up after stopping (days)", 0, 180, 0, step = 30),
      checkboxInput("surg", "Surgical debridement", FALSE),
      conditionalPanel("input.surg",
        sliderInput("surgday", "Surgery on therapy day", 1, 28, 7, step = 1),
        sliderInput("surglog", "log10 removed by debridement", 2, 6, 4, step = 0.5)),

      h4("Levers worth pulling"),
      sliderInput("lamvan", "lambda vancomycin (mm)", 0.15, 3.0, 0.45, step = 0.05),
      sliderInput("lamrif", "lambda rifampicin (mm)", 0.3, 4.0, 2.5, step = 0.1),
      sliderInput("mutrate", "log10 rpoB mutation frequency", -10, -6, -8, step = 0.5),
      checkboxInput("asa", "Aspirin 325 mg from diagnosis", FALSE),
      checkboxInput("noinoc", "Switch OFF the inoculum effect", FALSE),
      sliderInput("crcl", "Baseline CrCL (mL/min)", 20, 160, 100, step = 10)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. The ledger",
          br(),
          p(strong("Cure means driving the entire vegetation below one organism."),
            "From 10^9-10^10 CFU that is 9-10 log10 of sustained net killing,",
            "and it has to be delivered", em("in every shell"), "— an average",
            "over the lesion is not a cure. This table is the model's central",
            "accounting identity; everything else in the dashboard explains one",
            "of its rows."),
          verbatimTextOutput("ledger_txt"),
          plotOutput("p_ledger", height = "380px"),
          br(),
          h4("Where the shortfall is"),
          plotOutput("p_shells", height = "330px")),

        tabPanel("2. Patient & lesion",
          br(),
          p("The presenting state is", strong("earned"), "by simulating the",
            "untreated natural history from a seeding inoculum. Move the",
            "'lesion age at diagnosis' slider and watch a different disease",
            "appear: a 2 mm vegetation has essentially no protected core."),
          verbatimTextOutput("present_txt"),
          plotOutput("p_natural", height = "420px"),
          plotOutput("p_biomarker", height = "330px")),

        tabPanel("3. Pharmacokinetics",
          br(),
          p("Top: plasma. Bottom: the same drug at three depths inside the",
            "lesion. Note that the peak-to-trough excursion that dominates the",
            "plasma curve is almost gone by the mid shell — the vegetation",
            "sees a smoothed, delayed and heavily attenuated version of what",
            "the pharmacy delivered."),
          plotOutput("p_pk", height = "330px"),
          plotOutput("p_pkdepth", height = "380px"),
          verbatimTextOutput("pen_txt")),

        tabPanel("4. The exponential",
          br(),
          p(strong("This is the model in one figure."), "Concentration falls as",
            "C(x) = C_surface x exp(-x/lambda). Doubling the dose shifts the",
            "whole curve up by a factor of two; it does not change the slope.",
            "Beyond a few lambda, extra plasma exposure buys almost nothing —",
            "which is why there is an optimum plasma AUC and not a maximum."),
          plotOutput("p_profile", height = "430px"),
          br(),
          p("The second panel is the same physics applied to nutrient and",
            "oxygen. It is why the core is stationary-phase, why beta-lactam",
            "killing collapses there, and why gentamicin — which needs a",
            "proton-motive force to get inside a cell at all — does nothing",
            "in the place it is most needed."),
          plotOutput("p_avail", height = "330px")),

        tabPanel("5. Bacteriology",
          br(),
          p("Twelve bacterial states: three shells x growing/dormant x",
            "wild-type/rpoB-mutant. The mutants are not placed at t = 0; they",
            "are generated by a per-replication mutation term while the lesion",
            "grows, which is why there are already ~100 of them before the",
            "first dose."),
          plotOutput("p_bugs", height = "400px"),
          plotOutput("p_resist", height = "340px"),
          verbatimTextOutput("rif_txt")),

        tabPanel("6. Embolism & valve",
          br(),
          p("Embolic hazard = f(diameter) x", strong("fragility"),
            ", and fragility tracks the surface", em("bacterial"), "population,",
            "not the matrix. So the hazard collapses within days of effective",
            "therapy while the vegetation is still plainly visible on echo for",
            "weeks. That single structural choice is what makes the model",
            "reproduce the aspirin trial as a null."),
          plotOutput("p_emb", height = "400px"),
          plotOutput("p_valve", height = "360px"),
          verbatimTextOutput("surg_txt")),

        tabPanel("7. Toxicity",
          br(),
          p("AKI -> CrCL down -> vancomycin clearance down -> exposure up ->",
            "AKI. A genuine positive feedback with loop gain above one at high",
            "exposure, which is why the model can be pushed into a runaway and",
            "why the drive term has to be capped."),
          plotOutput("p_tox", height = "400px"),
          plotOutput("p_renal", height = "340px")),

        tabPanel("8. Regimen compare",
          br(),
          p("Regimen A against regimen B on the same patient. The three",
            "negative trials the model is required to reproduce are one click",
            "away below."),
          plotOutput("p_compare", height = "420px"),
          hr(),
          h4("The informative negatives"),
          p("A model that cannot reproduce a negative trial is a fitted curve,",
            "not a mechanism. Each of these is produced by a structural",
            "property put in for an independent reason."),
          actionButton("run_neg", "Run T5 gentamicin / T6 rifampicin / T7 aspirin"),
          verbatimTextOutput("neg_txt")),

        tabPanel("9. Exposure sweeps",
          br(),
          h4("Vancomycin AUC: a saturating benefit against an accelerating harm"),
          actionButton("run_auc", "Run the AUC sweep"),
          verbatimTextOutput("auc_txt"),
          hr(),
          h4("How much rides on the penetration length"),
          p("lambda is the least directly measurable parameter in the model and",
            "the most load-bearing. This is the honest place to look."),
          actionButton("run_lam", "Run the lambda sweep"),
          verbatimTextOutput("lam_txt"),
          hr(),
          h4("Duration and relapse"),
          actionButton("run_dur", "Run the duration sweep"),
          verbatimTextOutput("dur_txt"))
      )
    )
  ),
  hr(),
  tags$p(style = "color:#888;font-size:11px",
         "Research / education / hypothesis generation only. Not for clinical",
         "decision-making, prescribing, or regulatory use.")
)

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## parameters applied AFTER the lesion is grown -> they change the treatment
  post <- reactive({
    p <- list(LAM_VAN = input$lamvan, LAM_RIF = input$lamrif,
              MUTRATE = 10^input$mutrate)
    if (input$asa)    p$ASA_ON <- 1
    if (input$noinoc) p$INOC_BL_ON <- 0
    p
  })
  ## parameters applied BEFORE -> they change the patient
  pre <- reactive({
    p <- list(CRCL0 = input$crcl, MUTRATE = 10^input$mutrate)
    p
  })

  blpk <- function(which) {
    v <- if (which == "A") input$blA else input$blB
    if (v == "(from phenotype)") list() else IE_BLPK[[v]]
  }

  runarm <- function(rx, bl) {
    rx <- setdiff(rx, "none")
    if (!length(rx)) rx <- NULL
    sd <- if (input$surg) input$surgday else NULL
    pp <- c(post(), bl)
    if (!is.null(sd)) pp$SURG_LOG <- input$surglog
    IE_sim(MOD, rx, days = input$days, delta = 2, pheno = input$pheno,
           dx_day = input$dxday, param = pre(), param_post = pp,
           surgery_day = sd, follow_days = input$follow)
  }

  A <- reactive(runarm(input$rxA, blpk("A")))
  B <- reactive(runarm(input$rxB, blpk("B")))
  N <- reactive({
    IE_sim(MOD, NULL, days = input$days, delta = 2, pheno = input$pheno,
           dx_day = input$dxday, param = pre(), param_post = post(),
           follow_days = input$follow)
  })

  labA <- reactive(paste(setdiff(input$rxA, "none"), collapse = " + "))
  labB <- reactive(paste(setdiff(input$rxB, "none"), collapse = " + "))

  ## ---- Tab 1 : the ledger --------------------------------------------------
  output$ledger_txt <- renderPrint({
    mk <- function(o, nm) {
      a <- o[1, ]; z <- o[nrow(o), ]
      data.frame(regimen = nm,
                 required = round(a$LOGN, 2),
                 delivered = round(a$LOGN - z$LOGN, 2),
                 surface = round(a$LOGN_S - z$LOGN_S, 2),
                 mid = round(a$LOGN_M - z$LOGN_M, 2),
                 core = round(a$LOGN_C - z$LOGN_C, 2),
                 end_log10 = round(z$LOGN, 2),
                 sterile = ifelse(z$STERILE > 0.5, "YES", "no"),
                 cured = ifelse(z$CURED > 0.5, "YES", "no"),
                 stringsAsFactors = FALSE)
    }
    out <- rbind(mk(N(), "no therapy"),
                 mk(A(), paste0("A: ", labA())),
                 mk(B(), paste0("B: ", labB())))
    cat("Required = log10 CFU present at diagnosis. 'delivered' must reach it\n")
    cat("IN EVERY SHELL. 'sterile' is the vegetation; 'cured' also requires the\n")
    cat("intracellular reservoir to be gone.\n\n")
    print(out, row.names = FALSE)
  })

  output$p_ledger <- renderPlot({
    a <- A(); b <- B(); n <- N()
    lineplot(a$time / 24,
             list(n$LOGN, a$LOGN, b$LOGN),
             c(PAL["grey"], PAL["base"], PAL["drug"]),
             c("no therapy", paste0("A: ", labA()), paste0("B: ", labB())),
             ylab = "log10 total viable CFU in the vegetation",
             main = "Total burden — the quantity that has to reach zero",
             hline = c(0), hlab = c("1 CFU = cure"),
             ylim = c(-3.4, max(c(n$LOGN, a$LOGN, b$LOGN)) + 0.4))
  })

  output$p_shells <- renderPlot({
    a <- A()
    lineplot(a$time / 24,
             list(a$LOGN_S, a$LOGN_M, a$LOGN_C),
             c(PAL["surf"], PAL["mid"], PAL["core"]),
             c("surface shell", "mid shell", "CORE shell"),
             ylab = "log10 CFU", main = paste0("Regimen A by depth — ", labA()),
             hline = 0, hlab = "1 CFU", ylim = c(-3.4, 11))
  })

  ## ---- Tab 2 : patient -----------------------------------------------------
  output$present_txt <- renderPrint({
    a <- A()[1, ]
    s <- IE_seed(MOD, input$pheno, dx_day = input$dxday, param = pre())
    r <- with(s$init, RSG + RSD + RMG + RMD + RCG + RCD)
    cat("THE PATIENT AT THE MOMENT OF DIAGNOSIS (all of it emergent)\n\n")
    cat(sprintf("  vegetation diameter ............ %.1f mm\n", a$VEGMM))
    cat(sprintf("  total viable burden ............ 10^%.2f CFU\n", a$LOGN))
    cat(sprintf("  dormant fraction ............... %.3f\n", a$FRAC_DORM))
    cat(sprintf("  pre-existing rpoB mutants ...... %.0f CFU\n", r))
    cat(sprintf("  bacteraemia .................... %.2f CFU/mL (P(set +) = %.2f)\n",
                a$CBLD, a$P_BCPOS))
    cat(sprintf("  CRP / NT-proBNP ................ %.0f mg/L / %.0f pg/mL\n",
                a$CRP, a$NTBNP))
    cat(sprintf("  regurgitant fraction ........... %.2f\n", a$REG))
    cat(sprintf("  free vancomycin in the core .... %.3f mg/L (%.2f x MIC)\n",
                a$C_VAN_C, a$CORE_VAN_MIC))
    cat(sprintf("  surgical indications ........... HF %s | uncontrolled %s | embolic %s\n",
                ifelse(a$SI_HF > .5, "YES", "-"), ifelse(a$SI_UNC > .5, "YES", "-"),
                ifelse(a$SI_EMB > .5, "YES", "-")))
  })

  output$p_natural <- renderPlot({
    o <- IE_sim(MOD, NULL, days = 45, delta = 6, pheno = input$pheno,
                dx_day = 1, param = pre(), param_post = post())
    gpar(mfrow = c(1, 2))
    plot(o$time / 24, o$VEGMM, type = "l", lwd = 2.4, col = PAL["drug"],
         xlab = "day of the lesion", ylab = "vegetation diameter (mm)",
         main = "Untreated natural history"); grid(col = "#e8e8e8", lty = 1)
    abline(v = input$dxday, lty = 2, col = PAL["base"])
    abline(h = c(10, 15), lty = 3, col = "#999999")
    text(2, 10, "10 mm", pos = 3, cex = 0.75, col = "#777777")
    text(2, 15, "15 mm", pos = 3, cex = 0.75, col = "#777777")
    text(input$dxday, par("usr")[3] + 1, "diagnosis", pos = 4, cex = 0.8,
         col = PAL["base"])
    plot(o$time / 24, o$LOGN, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day of the lesion", ylab = "log10 CFU / dormant fraction x 10",
         main = "Burden and the dormant fraction", ylim = c(0, 11))
    grid(col = "#e8e8e8", lty = 1)
    lines(o$time / 24, o$FRAC_DORM * 10, lwd = 2.4, col = PAL["purp"], lty = 2)
    abline(v = input$dxday, lty = 2, col = PAL["base"])
    legend("bottomright", c("log10 CFU", "dormant fraction x 10"),
           col = c(PAL["base"], PAL["purp"]), lwd = 2.2, lty = c(1, 2), bty = "n",
           cex = 0.85)
  })

  output$p_biomarker <- renderPlot({
    a <- A(); n <- N()
    gpar(mfrow = c(1, 3))
    for (v in list(c("CRP", "CRP (mg/L)"), c("PCT", "procalcitonin (ng/mL)"),
                   c("NTBNP", "NT-proBNP (pg/mL)"))) {
      plot(a$time / 24, a[[v[1]]], type = "l", lwd = 2.4, col = PAL["base"],
           xlab = "day", ylab = v[2], main = v[2],
           ylim = range(c(a[[v[1]]], n[[v[1]]])))
      grid(col = "#e8e8e8", lty = 1)
      lines(n$time / 24, n[[v[1]]], lwd = 2, col = PAL["grey"], lty = 2)
    }
  })

  ## ---- Tab 3 : PK ----------------------------------------------------------
  output$p_pk <- renderPlot({
    a <- A()
    vs <- c(CVANP = "vancomycin", CDAPP = "daptomycin", CBLP = "beta-lactam",
            CGENP = "gentamicin", CRIFP = "rifampicin", CDALP = "dalbavancin",
            CORLP = "oral")
    keep <- names(vs)[sapply(names(vs), function(v) max(a[[v]]) > 1e-6)]
    if (!length(keep)) { plot.new(); title("no drug in regimen A"); return() }
    cols <- c(PAL["base"], PAL["drug"], PAL["alt"], PAL["warn"], PAL["purp"],
              PAL["grey"], PAL["core"])[seq_along(keep)]
    lineplot(a$time / 24, lapply(keep, function(v) a[[v]]), cols, vs[keep],
             ylab = "total plasma concentration (mg/L)",
             main = paste0("Plasma PK — ", labA()))
  })

  output$p_pkdepth <- renderPlot({
    a <- A()
    cand <- list(c("CVANP", "C_VAN_S", "C_VAN_C", "vancomycin"),
                 c("CDAPP", "C_DAP_C", "C_DAP_C", "daptomycin"),
                 c("CBLP",  "C_BL_C",  "C_BL_C",  "beta-lactam"),
                 c("CRIFP", "C_RIF_C", "C_RIF_C", "rifampicin"),
                 c("CGENP", "C_GEN_C", "C_GEN_C", "gentamicin"))
    act <- Filter(function(k) max(a[[k[1]]]) > 1e-6, cand)
    if (!length(act)) { plot.new(); title("no drug in regimen A"); return() }
    gpar(mfrow = c(1, min(3, length(act))))
    for (k in act[seq_len(min(3, length(act)))]) {
      yy <- cbind(a[[k[1]]], a[[k[3]]])
      plot(a$time / 24, yy[, 1], type = "l", lwd = 2.2, col = PAL["grey"],
           log = "y", xlab = "day", ylab = "mg/L (log scale)",
           main = paste0(k[4], ": plasma vs CORE"),
           ylim = c(max(1e-5, min(yy[yy > 0])), max(yy) * 1.4))
      grid(col = "#e8e8e8", lty = 1)
      lines(a$time / 24, yy[, 2], lwd = 2.4, col = PAL["core"])
      legend("bottomleft", c("total plasma", "free, vegetation CORE"),
             col = c(PAL["grey"], PAL["core"]), lwd = 2.2, bty = "n", cex = 0.8)
    }
  })

  output$pen_txt <- renderPrint({
    invisible(capture.output(x <- IE_penetration(MOD, pheno = input$pheno,
                                                 days = 7, dx_day = input$dxday)))
    cat("Steady-state exposure ratios, geometry held fixed at the presenting\n")
    cat("size. The last three columns contain no dose.\n\n")
    print(x, row.names = FALSE)
    cat("\ncore_over_freeplasma is what the surviving bacteria experience.\n")
    cat("homogenate_over_free is what a ground-up vegetation assay measures --\n")
    cat("a volume-weighted average dominated by the outer shells. The published\n")
    cat("0.2-0.4 vegetation:serum ratios are higher still because they are TOTAL\n")
    cat("drug, including everything bound to the fibrin -- which is the very\n")
    cat("binding that makes lambda short in the first place.\n")
  })

  ## ---- Tab 4 : the exponential --------------------------------------------
  output$p_profile <- renderPlot({
    a <- A(); z <- a[nrow(a), ]; a0 <- a[1, ]
    pp <- as.list(param(MOD))
    lam <- c(vancomycin = input$lamvan, daptomycin = pp$LAM_DAP,
             `beta-lactam` = pp$LAM_BL, gentamicin = pp$LAM_GEN,
             rifampicin = input$lamrif)
    R <- a0$VEGMM / 2
    x <- seq(0, R, length.out = 200)
    gpar()
    plot(NA, xlim = c(0, R), ylim = c(1e-4, 1.2), log = "y",
         xlab = "depth from the vegetation surface (mm)",
         ylab = "concentration relative to the surface",
         main = sprintf("C(x) = C_surface x exp(-x/lambda)   [vegetation %.1f mm at diagnosis]",
                        a0$VEGMM))
    grid(col = "#e8e8e8", lty = 1)
    cols <- c(PAL["base"], PAL["drug"], PAL["alt"], PAL["warn"], PAL["purp"])
    for (i in seq_along(lam)) lines(x, exp(-x / lam[i]), lwd = 2.4, col = cols[i])
    dd <- c(0.13, 0.40, 0.72) * R
    abline(v = dd, lty = 3, col = "#888888")
    text(dd, 1.0, c("surface", "mid", "CORE"), pos = 4, cex = 0.78, col = "#666666")
    legend("bottomleft", sprintf("%s (lambda = %.2f mm)", names(lam), lam),
           col = cols, lwd = 2.4, bty = "n", cex = 0.85)
  })

  output$p_avail <- renderPlot({
    a0 <- A()[1, ]; pp <- as.list(param(MOD))
    R <- a0$VEGMM / 2; x <- seq(0, R, length.out = 200)
    av <- exp(-x / pp$LAM_NUT); ox <- av^1.5
    gpar()
    plot(x, av, type = "l", lwd = 2.6, col = PAL["alt"], ylim = c(0, 1.05),
         xlab = "depth from the surface (mm)", ylab = "relative availability",
         main = "Nutrient and oxygen fall on the same physics, a different length scale")
    grid(col = "#e8e8e8", lty = 1)
    lines(x, ox, lwd = 2.6, col = PAL["warn"])
    lines(x, 1 - av, lwd = 2.2, col = PAL["purp"], lty = 2)
    abline(v = c(0.13, 0.40, 0.72) * R, lty = 3, col = "#888888")
    legend("topright",
           c("nutrient availability", "oxygen (drives aminoglycoside uptake)",
             "-> drive into DORMANCY"),
           col = c(PAL["alt"], PAL["warn"], PAL["purp"]), lwd = 2.4,
           lty = c(1, 1, 2), bty = "n", cex = 0.85)
  })

  ## ---- Tab 5 : bacteriology ------------------------------------------------
  output$p_bugs <- renderPlot({
    a <- A()
    gpar(mfrow = c(1, 2))
    plot(NA, xlim = range(a$time / 24), ylim = c(-3.4, 11), xlab = "day",
         ylab = "log10 CFU", main = paste0("By shell — ", labA()))
    grid(col = "#e8e8e8", lty = 1)
    lines(a$time / 24, a$LOGN_S, lwd = 2.4, col = PAL["surf"])
    lines(a$time / 24, a$LOGN_M, lwd = 2.4, col = PAL["mid"])
    lines(a$time / 24, a$LOGN_C, lwd = 2.4, col = PAL["core"])
    abline(h = 0, lty = 3)
    legend("topright", c("surface", "mid", "core"),
           col = c(PAL["surf"], PAL["mid"], PAL["core"]), lwd = 2.4, bty = "n")
    plot(a$time / 24, a$FRAC_DORM, type = "l", lwd = 2.6, col = PAL["purp"],
         ylim = c(0, 1), xlab = "day", ylab = "fraction",
         main = "Dormant fraction — tolerance, not resistance")
    grid(col = "#e8e8e8", lty = 1)
    legend("bottomright", "dormant fraction of survivors", col = PAL["purp"],
           lwd = 2.4, bty = "n", cex = 0.85)
  })

  output$p_resist <- renderPlot({
    a <- A(); b <- B()
    lineplot(a$time / 24, list(a$LOGR, b$LOGR),
             c(PAL["base"], PAL["drug"]),
             c(paste0("A: ", labA()), paste0("B: ", labB())),
             ylab = "log10 rpoB-mutant CFU",
             main = "The resistant subpopulation that was already there",
             hline = 0, hlab = "1 CFU", ylim = c(-3.4, 11))
  })

  output$rif_txt <- renderPrint({
    s <- IE_seed(MOD, input$pheno, dx_day = input$dxday, param = pre())
    n <- with(s$init, BSG + BSD + BMG + BMD + BCG + BCD)
    r <- with(s$init, RSG + RSD + RMG + RMD + RCG + RCD)
    f <- 10^input$mutrate
    cat("PRE-EXISTING RESISTANCE — arithmetic, not bad luck\n\n")
    cat(sprintf("  wild-type burden ............. %.3e CFU\n", n))
    cat(sprintf("  rpoB mutants (emergent) ...... %.3e CFU\n", r))
    cat(sprintf("  observed mutant fraction ..... %.3e\n", r / (n + r)))
    cat(sprintf("  mutation frequency f_mut ..... %.3e\n", f))
    cat(sprintf("  naive f_mut x N .............. %.3e CFU\n", f * n))
    cat("\nNothing is placed at t = 0: the mutation term is f_mut x mu x B, so the\n")
    cat("mutants are a by-product of growing the lesion. The excess over the\n")
    cat("naive product is Luria-Delbruck accumulation over the growth phase.\n")
    cat("Rifampicin monotherapy therefore fails DETERMINISTICALLY. Select 'rif0'\n")
    cat("alone as regimen A and watch it.\n")
  })

  ## ---- Tab 6 : embolism ----------------------------------------------------
  output$p_emb <- renderPlot({
    a <- A(); b <- B(); n <- N()
    gpar(mfrow = c(1, 2))
    plot(NA, xlim = range(a$time / 24), ylim = c(0, max(n$P_EMB) * 1.05),
         xlab = "day of therapy", ylab = "cumulative P(embolic event)",
         main = "Embolic probability")
    grid(col = "#e8e8e8", lty = 1)
    lines(n$time / 24, n$P_EMB, lwd = 2.2, col = PAL["grey"], lty = 2)
    lines(a$time / 24, a$P_EMB, lwd = 2.4, col = PAL["base"])
    lines(b$time / 24, b$P_EMB, lwd = 2.4, col = PAL["drug"])
    legend("topleft", c("no therapy", paste0("A: ", labA()), paste0("B: ", labB())),
           col = c(PAL["grey"], PAL["base"], PAL["drug"]), lwd = 2.2,
           lty = c(2, 1, 1), bty = "n", cex = 0.85)
    plot(a$time / 24, a$VEGMM, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day of therapy", ylab = "vegetation diameter (mm)",
         main = "…and the lesion is still there",
         ylim = c(0, max(c(a$VEGMM, n$VEGMM))))
    grid(col = "#e8e8e8", lty = 1)
    lines(n$time / 24, n$VEGMM, lwd = 2.2, col = PAL["grey"], lty = 2)
    lines(a$time / 24, 10 * a$LOGN_S / max(a$LOGN_S), lwd = 2, col = PAL["surf"],
          lty = 3)
    legend("topright", c("diameter, treated", "diameter, untreated",
                         "surface burden (scaled)"),
           col = c(PAL["base"], PAL["grey"], PAL["surf"]), lwd = 2.2,
           lty = c(1, 2, 3), bty = "n", cex = 0.8)
  })

  output$p_valve <- renderPlot({
    a <- A(); n <- N()
    gpar(mfrow = c(1, 3))
    plot(a$time / 24, a$REG, type = "l", lwd = 2.4, col = PAL["drug"],
         xlab = "day", ylab = "regurgitant fraction",
         main = "Valve destruction is a RATCHET",
         ylim = c(0, max(c(a$REG, n$REG))))
    grid(col = "#e8e8e8", lty = 1)
    lines(n$time / 24, n$REG, lwd = 2.2, col = PAL["grey"], lty = 2)
    plot(a$time / 24, a$NTBNP, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day", ylab = "NT-proBNP (pg/mL)", main = "Wall stress",
         ylim = c(0, max(c(a$NTBNP, n$NTBNP))))
    grid(col = "#e8e8e8", lty = 1)
    lines(n$time / 24, n$NTBNP, lwd = 2.2, col = PAL["grey"], lty = 2)
    plot(a$time / 24, a$SI_HF + 2 * a$SI_UNC + 4 * a$SI_EMB, type = "s",
         lwd = 2.4, col = PAL["warn"], xlab = "day", ylab = "indication code",
         main = "Surgical indications (1 HF, 2 uncontrolled, 4 embolic)",
         ylim = c(0, 7))
    grid(col = "#e8e8e8", lty = 1)
  })

  output$surg_txt <- renderPrint({
    a <- A()
    pre_h <- a$P_EMB[1]
    tot <- a$P_EMB[nrow(a)]
    cat(sprintf("P(embolism) accrued BEFORE diagnosis ......... %.4f\n", pre_h))
    cat(sprintf("P(embolism) accrued during therapy .......... %.4f\n",
                1 - exp(log(1 - tot) - log(1 - pre_h))))
    cat(sprintf("Vegetation %.1f mm at diagnosis -> %.1f mm at the end\n",
                a$VEGMM[1], a$VEGMM[nrow(a)]))
    cat("\nMost of the hazard a patient carries has already been discharged by the\n")
    cat("time anyone knows they have endocarditis. That is not a nihilistic\n")
    cat("observation; it is the reason the surgical window for EMBOLISM\n")
    cat("PREVENTION is measured in days and not in weeks.\n")
  })

  ## ---- Tab 7 : toxicity ----------------------------------------------------
  output$p_tox <- renderPlot({
    a <- A(); b <- B()
    gpar(mfrow = c(1, 3))
    plot(a$time / 24, a$SCR, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day", ylab = "serum creatinine (mg/dL)", main = "Renal function",
         ylim = c(0.6, max(c(a$SCR, b$SCR)) * 1.1))
    grid(col = "#e8e8e8", lty = 1)
    lines(b$time / 24, b$SCR, lwd = 2.4, col = PAL["drug"])
    abline(h = a$SCR[1] * c(1.5, 2, 3), lty = 3, col = "#999999")
    legend("topleft", c(paste0("A: ", labA()), paste0("B: ", labB())),
           col = c(PAL["base"], PAL["drug"]), lwd = 2.2, bty = "n", cex = 0.8)
    plot(a$time / 24, a$CK, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day", ylab = "CK (U/L)", main = "Creatine kinase (daptomycin)",
         ylim = c(0, max(c(a$CK, b$CK)) * 1.1))
    grid(col = "#e8e8e8", lty = 1); lines(b$time / 24, b$CK, lwd = 2.4, col = PAL["drug"])
    plot(a$time / 24, a$ALT, type = "l", lwd = 2.4, col = PAL["base"],
         xlab = "day", ylab = "ALT (U/L)", main = "ALT (rifampicin, nafcillin)",
         ylim = c(0, max(c(a$ALT, b$ALT)) * 1.1))
    grid(col = "#e8e8e8", lty = 1); lines(b$time / 24, b$ALT, lwd = 2.4, col = PAL["drug"])
  })

  output$p_renal <- renderPlot({
    a <- A(); b <- B()
    lineplot(a$time / 24, list(a$CRCLT, b$CRCLT),
             c(PAL["base"], PAL["drug"]),
             c(paste0("A: ", labA()), paste0("B: ", labB())),
             ylab = "CrCL (mL/min)",
             main = "The feedback loop: CrCL is a STATE, not a covariate")
  })

  ## ---- Tab 8 : compare -----------------------------------------------------
  output$p_compare <- renderPlot({
    a <- A(); b <- B(); n <- N()
    gpar(mfrow = c(2, 2))
    for (v in list(c("LOGN", "log10 total CFU"),
                   c("LOGN_C", "log10 CORE CFU"),
                   c("P_BCPOS", "P(blood culture set positive)"),
                   c("VEGMM", "vegetation diameter (mm)"))) {
      yl <- range(c(a[[v[1]]], b[[v[1]]], n[[v[1]]]))
      plot(a$time / 24, a[[v[1]]], type = "l", lwd = 2.4, col = PAL["base"],
           xlab = "day", ylab = v[2], main = v[2], ylim = yl)
      grid(col = "#e8e8e8", lty = 1)
      lines(b$time / 24, b[[v[1]]], lwd = 2.4, col = PAL["drug"])
      lines(n$time / 24, n[[v[1]]], lwd = 2, col = PAL["grey"], lty = 2)
      if (v[1] %in% c("LOGN", "LOGN_C")) abline(h = 0, lty = 3)
      if (v[1] == "P_BCPOS") abline(h = 0.25, lty = 3)
    }
  })

  neg <- eventReactive(input$run_neg,
    paste(capture.output(IE_negatives(MOD, dx_day = input$dxday)), collapse = "\n"))
  output$neg_txt <- renderText(neg())

  ## ---- Tab 9 : sweeps ------------------------------------------------------
  auc <- eventReactive(input$run_auc,
    paste(capture.output(IE_aucsweep(MOD, pheno = input$pheno,
                                     days = input$days, dx_day = input$dxday)),
          collapse = "\n"))
  output$auc_txt <- renderText(auc())

  lam <- eventReactive(input$run_lam,
    paste(capture.output(IE_lambda(MOD, pheno = input$pheno,
                                   days = input$days, dx_day = input$dxday)),
          collapse = "\n"))
  output$lam_txt <- renderText(lam())

  dur <- eventReactive(input$run_dur, {
    rx <- setdiff(input$rxA, "none"); if (!length(rx)) rx <- "dapto8"
    paste(capture.output(IE_duration(MOD, pheno = input$pheno, regimen = rx,
                                     dx_day = input$dxday)), collapse = "\n")
  })
  output$dur_txt <- renderText(dur())
}

shinyApp(ui, server)
