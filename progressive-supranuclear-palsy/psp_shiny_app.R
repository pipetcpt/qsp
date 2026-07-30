## ============================================================================
##  Progressive Supranuclear Palsy (PSP) QSP model -- Shiny dashboard
##  ---------------------------------------------------------------------------
##  psp_shiny_app.R          12 tabs
##
##  RUN:
##    setwd("<this directory>")
##    shiny::runApp("psp_shiny_app.R")
##
##  The app is a front end for psp_mrgsolve_model.R and adds nothing to the
##  model.  Its purpose is to make the six quantities the model exists to
##  separate visible at the same time:
##
##    1. the two tau pools that differ by 10^3-10^4               (tab 2)
##    2. the front moving through eight nuclei                     (tab 3)
##    3. antibody exposure vs the pool it can reach                (tab 5)
##    4. the two-integration lag between drug and PSPRS            (tab 6)
##    5. the biomarker that misled and the ones that would not     (tab 7)
##    6. the same drug given before vs after symptom onset         (tab 10)
##
##  Base graphics only (no ggplot2 dependency), so the app runs anywhere
##  mrgsolve does.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
})

## ---- load the model --------------------------------------------------------
.here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
if (is.null(.here) || !nzchar(.here)) .here <- "."
.mfile <- file.path(.here, "psp_mrgsolve_model.R")
if (!file.exists(.mfile)) .mfile <- "psp_mrgsolve_model.R"
source(.mfile, local = FALSE)

SCEN   <- psp_scenarios()
SC_IDS <- names(SCEN)
SC_LAB <- vapply(SCEN, function(x) x$label, character(1))
Y      <- 365.25

REGION <- c("R1 STN", "R2 GPi/SNr", "R3 SNc", "R4 midbrain tegmentum",
            "R5 PPN", "R6 dentate", "R7 frontal cortex", "R8 LC/bulbar")
RCOL   <- c("#8c510a", "#bf812d", "#dfc27d", "#c51b7d",
            "#4d9221", "#35978f", "#2166ac", "#762a83")

## ---- plotting helpers ------------------------------------------------------
pal <- c("#2166ac", "#b2182b", "#1b7837", "#762a83", "#d95f02",
         "#7570b3", "#e7298a", "#66a61e")

vline_marks <- function() {
  abline(v = TSYMPT/Y, lty = 3, col = "grey50")
  abline(v = TENROL/Y, lty = 2, col = "grey20")
}
mark_legend <- function() {
  legend("topleft", bty = "n", cex = 0.75, lty = c(3, 2),
         col = c("grey50", "grey20"),
         legend = c("symptom onset", "trial enrolment (PSPRS 38)"))
}

lp <- function(d, cols, ylab, main, log = "", labels = cols, lwd = 2,
               cols_pal = pal, legpos = "topleft") {
  x <- d$time/Y
  ymat <- as.matrix(d[, cols, drop = FALSE])
  if (log == "y") ymat[ymat <= 0] <- NA
  yl <- range(ymat, na.rm = TRUE)
  if (!all(is.finite(yl))) yl <- c(0, 1)
  plot(NA, xlim = range(x), ylim = yl, xlab = "years from first tau aggregate",
       ylab = ylab, main = main, log = log, las = 1)
  grid(col = "grey92")
  for (i in seq_along(cols))
    lines(x, ymat[, i], col = cols_pal[(i - 1) %% length(cols_pal) + 1], lwd = lwd)
  vline_marks()
  if (length(cols) > 1)
    legend(legpos, legend = labels, col = cols_pal[seq_along(cols)],
           lwd = lwd, bty = "n", cex = 0.8)
  invisible(NULL)
}

## ============================================================================
##  UI
## ============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, Helvetica, Arial, sans-serif; }
    .note { font-size: 12px; color: #444; background: #f6f6f6;
            border-left: 3px solid #b2182b; padding: 8px 10px; margin: 6px 0 12px 0; }
    .num  { font-family: ui-monospace, Menlo, monospace; font-size: 12px; }
    h4 { margin-top: 4px; }
  "))),

  titlePanel("Progressive Supranuclear Palsy -- QSP model explorer"),
  div(class = "note",
      strong("The question this model exists to answer. "),
      "Gosuranemab removed 98% of CSF unbound N-terminal tau and changed the PSPRS by ",
      "-0.2 points. NIO752 removed 20% of CSF total tau -- five times less -- and produced ",
      "the only downstream signal in this disease's history. A model with one ",
      "\"tau lowering\" axis cannot hold both facts, so this one separates tau into ",
      "pools that differ by three to four orders of magnitude in size and by everything ",
      "in causal role."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("scen", "Scenario", choices = setNames(SC_IDS, SC_LAB),
                  selected = "S02"),
      selectInput("ref", "Comparator", choices = setNames(SC_IDS, SC_LAB),
                  selected = "S00"),
      sliderInput("endy", "Simulate to (years)", min = 8, max = 20,
                  value = 15, step = 1),
      hr(),
      h4("Patient profile"),
      sliderInput("gmapt", "MAPT H1 haplotype dose", min = 0, max = 2,
                  value = 2, step = 1),
      sliderInput("trim11", "TRIM11 disaggregase capacity", min = 0.3, max = 2,
                  value = 1, step = 0.1),
      radioButtons("orig", "Wave origin (phenotype)",
                   choices = c("STN + midbrain tegmentum (PSP-RS)" = 1,
                               "SNc / nigrostriatal (PSP-P)" = 0),
                   selected = 1),
      hr(),
      h4("Structural parameters"),
      sliderInput("phiacc", HTML("&Phi;<sub>ACC</sub> -- transfer flux exposed to bulk ISF"),
                  min = 0, max = 1, value = 0.45, step = 0.05),
      sliderInput("asodeep", "ASO knockdown reaching the deep-nucleus front",
                  min = 0.1, max = 1, value = 0.4, step = 0.05),
      sliderInput("kpab", "mAb brain:plasma partition (%)",
                  min = 0.05, max = 5, value = 0.15, step = 0.05),
      sliderInput("klag", "Committed-damage transit rate (1/day)",
                  min = 0.002, max = 0.03, value = 0.0055, step = 0.0005),
      hr(),
      checkboxInput("peg", "PEG gastrostomy at enrolment", FALSE),
      helpText(HTML("Dotted line = symptom onset. Dashed line = trial enrolment,
                     defined as the model's own PSPRS = 38 crossing rather than
                     assumed."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1 | Patient & course",
                 div(class = "note",
                     "The pre-symptomatic interval is an OUTPUT: enrolment is placed at the
                      model's own PSPRS = 38 crossing, and symptom onset 3.5 y before it,
                      because that is the median onset-to-enrolment interval in PASSPORT
                      and ARISE. Everything before the dotted line is disease no trial
                      has ever seen."),
                 fluidRow(column(6, plotOutput("p_course", height = 300)),
                          column(6, plotOutput("p_dom", height = 300))),
                 verbatimTextOutput("t_profile")),

        tabPanel("2 | Two tau pools",
                 div(class = "note",
                     HTML("<b>AXIS 1.</b> Intraneuronal monomer is ~2000 nM; the ISF pool an
                      antibody can reach is ~1 nM. F_EXT is their ratio. Clearing 98% of
                      the smaller pool removes a fraction of a per cent of brain tau, and
                      that is the ceiling on any extracellular-only mechanism at ANY dose
                      or affinity.")),
                 fluidRow(column(6, plotOutput("p_pools", height = 320)),
                          column(6, plotOutput("p_fext", height = 320))),
                 verbatimTextOutput("t_pools")),

        tabPanel("3 | The travelling wave",
                 div(class = "note",
                     HTML("<b>AXIS 4.</b> Each nucleus aggregates LOGISTICALLY and saturates in
                      ~1.5-2 y. A sum of saturating sigmoids is linear only if regions are
                      recruited sequentially, so the famous ~11 PSPRS points/year is
                      evidence of a front -- and it means every already-involved nucleus
                      sits where changing the local rate constant does nothing.")),
                 fluidRow(column(6, plotOutput("p_wave", height = 340)),
                          column(6, plotOutput("p_neur", height = 340))),
                 plotOutput("p_order", height = 210)),

        tabPanel("4 | Pharmacokinetics",
                 fluidRow(column(6, plotOutput("p_pkmab", height = 300)),
                          column(6, plotOutput("p_pkaso", height = 300))),
                 fluidRow(column(12, plotOutput("p_pksm", height = 300)))),

        tabPanel("5 | Target engagement",
                 div(class = "note",
                     HTML("<b>AXIS 2 and 3.</b> The -98% is not fitted: 2000 mg q4w at a 0.15%
                      brain partition gives a low-nanomolar ISF antibody concentration,
                      which saturates a ~1 nM antigen and is irrelevant to a 2000 nM one.
                      <b>98% engagement is itself the evidence that the pool was small.</b>
                      Note that TOTAL CSF tau goes UP while unbound goes down -- one
                      mechanism, two readouts, opposite signs.")),
                 fluidRow(column(6, plotOutput("p_eng", height = 320)),
                          column(6, plotOutput("p_updown", height = 320))),
                 verbatimTextOutput("t_eng")),

        tabPanel("6 | Clinical endpoints",
                 div(class = "note",
                     HTML("<b>AXIS 7.</b> Two integrations and one reserve threshold sit between
                      the drug and the score. Vertical saccade velocity is the exception:
                      it reads riMLF survival almost without a threshold, which is why it
                      moves first and why it FLOORS while the disease continues.")),
                 fluidRow(column(6, plotOutput("p_psprs", height = 320)),
                          column(6, plotOutput("p_vsv", height = 320))),
                 fluidRow(column(6, plotOutput("p_domains", height = 300)),
                          column(6, verbatimTextOutput("t_endpoint")))),

        tabPanel("7 | Fluid biomarkers",
                 div(class = "note",
                     HTML("The decoy and the discriminators. <b>CSF unbound N-terminal tau</b>
                      is a mass measurement of an abundant, non-seeding fragment pool --
                      exactly what an N-terminal antibody binds. <b>MTBR-tau243</b>,
                      <b>CSF seeding activity</b> and <b>NfL</b> track quantities the
                      mechanism actually passes through.")),
                 fluidRow(column(6, plotOutput("p_nfl", height = 300)),
                          column(6, plotOutput("p_tauassays", height = 300))),
                 fluidRow(column(6, plotOutput("p_gfap", height = 280)),
                          column(6, verbatimTextOutput("t_bm")))),

        tabPanel("8 | Imaging",
                 fluidRow(column(6, plotOutput("p_mid", height = 300)),
                          column(6, plotOutput("p_mrpi", height = 300))),
                 div(class = "note",
                     "MRPI beyond ~20 is outside the range in which the index was validated;
                      the model reports it anyway rather than clipping it.")),

        tabPanel("9 | Scenario comparison",
                 div(class = "note",
                     HTML("<b>vs_pbo</b> is the trial endpoint (change from baseline).
                      <b>abs_vs_pbo</b> is the absolute score. They disagree in sign for
                      any intervention that DELAYS the disease, because a delayed patient
                      sits lower on the curve at the fixed enrolment time and therefore
                      has more room to change. Look at S09 and S27.")),
                 sliderInput("chz", "Treatment horizon (weeks)", min = 26, max = 156,
                             value = 52, step = 2, width = "60%"),
                 actionButton("go_cmp", "Run all 29 scenarios"),
                 tableOutput("t_compare")),

        tabPanel("10 | When to treat",
                 div(class = "note",
                     HTML("<b>The single sharpest number in the model.</b> The same MAPT
                      knockdown is worth ~nothing if it starts at enrolment and years of
                      delay if it starts at biological onset. The treatment-phase
                      elasticity dln(slope)/dln(M) is ~0.02, so halving the PSPRS slope is
                      unreachable at ANY knockdown once symptoms exist.")),
                 actionButton("go_el", "Compute elasticity (~30 s)"),
                 fluidRow(column(6, plotOutput("p_el1", height = 320)),
                          column(6, plotOutput("p_el2", height = 320))),
                 tableOutput("t_el")),

        tabPanel("11 | Survival",
                 div(class = "note",
                     HTML("Competing risks written properly: dCIF/dt = h&middot;S, not summed
                      hazards. Aspiration pneumonia via bulbar failure is the dominant
                      channel, and PEG moves that hazard without touching the PSPRS.")),
                 fluidRow(column(6, plotOutput("p_surv", height = 320)),
                          column(6, plotOutput("p_cif", height = 320))),
                 verbatimTextOutput("t_surv")),

        tabPanel("12 | Validation",
                 div(class = "note",
                     "Every row is a number published by someone else. FIT rows were used to
                      calibrate; held-out rows were not. Two failures are reported rather
                      than repaired -- see README.md."),
                 actionButton("go_val", "Run validation (~30 s)"),
                 tableOutput("t_val"))
      )
    )
  )
)

## ============================================================================
##  SERVER
## ============================================================================
server <- function(input, output, session) {

  user_par <- reactive({
    list(GMAPT   = input$gmapt,
         TRIM11F = input$trim11,
         ORIG_RS = as.numeric(input$orig),
         PHI_ACC = input$phiacc,
         ASO_DEEP= input$asodeep,
         KP_AB   = input$kpab/100,
         KLAG    = input$klag,
         PEG     = as.numeric(input$peg))
  })

  simA <- reactive({
    psp_run(input$scen, end = round(input$endy*Y), delta = 7,
            extra_par = user_par())
  })
  simB <- reactive({
    psp_run(input$ref, end = round(input$endy*Y), delta = 7,
            extra_par = user_par())
  })

  at <- function(d, t, c) d[[c]][which.min(abs(d$time - t))]

  two <- function(col, ylab, main, log = "") {
    a <- simA(); b <- simB()
    x <- a$time/Y
    ym <- cbind(b[[col]], a[[col]])
    if (log == "y") ym[ym <= 0] <- NA
    yl <- range(ym, na.rm = TRUE); if (!all(is.finite(yl))) yl <- c(0, 1)
    plot(NA, xlim = range(x), ylim = yl, log = log, las = 1,
         xlab = "years from first tau aggregate", ylab = ylab, main = main)
    grid(col = "grey92")
    lines(x, ym[, 1], col = "grey45", lwd = 2, lty = 2)
    lines(x, ym[, 2], col = "#b2182b", lwd = 2.4)
    vline_marks()
    legend("topleft", bty = "n", cex = 0.78, lwd = c(2, 2.4), lty = c(2, 1),
           col = c("grey45", "#b2182b"),
           legend = c(paste0(input$ref, " (comparator)"), paste0(input$scen, " (scenario)")))
  }

  ## ---- 1. course -----------------------------------------------------------
  output$p_course <- renderPlot(two("PSPRS", "PSPRS (0-100)", "Disease course"))
  output$p_dom <- renderPlot({
    a <- simA()
    lp(a, c("ATOT", "NTOT"), "fraction",
       "mean aggregate load and surviving neurons",
       labels = c("mean aggregate load A", "mean surviving neurons N"))
  })
  output$t_profile <- renderPrint({
    a <- simA()
    ms <- { i <- which(a$SURV <= 0.5); if (length(i)) a$time[i[1]]/Y else NA }
    cat("Model-implied timeline (OUTPUTS, not inputs)\n")
    cat("  first tau aggregate           t = 0\n")
    cat(sprintf("  symptom onset                 %.2f y  (pre-symptomatic tauopathy)\n", TSYMPT/Y))
    cat(sprintf("  trial-entry severity, PSPRS38 %.2f y\n", TENROL/Y))
    cat(sprintf("  median survival               %.2f y from onset of tau pathology\n", ms))
    cat(sprintf("                                %.2f y from symptom onset\n", ms - TSYMPT/Y))
    cat("\nAt enrolment in this scenario\n")
    cat(sprintf("  PSPRS %.1f | mPSPRS-10 %.1f | SEADL %.0f%% | saccade %.0f deg/s\n",
                at(a, TENROL, "PSPRS"), at(a, TENROL, "MP10"),
                at(a, TENROL, "SEADL"), at(a, TENROL, "VSV")))
    cat(sprintf("  midbrain %.0f mm2 | MRPI %.1f | CSF NfL %.0f pg/mL\n",
                at(a, TENROL, "MIDA"), at(a, TENROL, "MRPI"), at(a, TENROL, "NFLC")))
  })

  ## ---- 2. pools ------------------------------------------------------------
  output$p_pools <- renderPlot({
    a <- simA()
    lp(a, c("TMON", "ETN", "ETS", "CTN"), "nM (log scale)",
       "the pools that differ by 10^3-10^4", log = "y",
       labels = c("intraneuronal monomer (~2000 nM)",
                  "ISF N-terminal fragments (~1 nM)",
                  "ISF seed-competent (~0.02 nM)",
                  "CSF N-terminal tau"))
  })
  output$p_fext <- renderPlot(two("F_EXT", "ISF tau / intraneuronal tau",
                                  "F_EXT -- the ceiling on extracellular mechanisms"))
  output$t_pools <- renderPrint({
    a <- simA()
    cat(sprintf("At enrolment:  intraneuronal monomer %.0f nM\n", at(a, TENROL, "TMON")))
    cat(sprintf("               ISF N-terminal        %.3f nM\n", at(a, TENROL, "ETN")))
    cat(sprintf("               ISF seed-competent    %.4f nM\n", at(a, TENROL, "ETS")))
    cat(sprintf("               F_EXT                 %.2e\n", at(a, TENROL, "F_EXT")))
    cat(sprintf("\nMaximum reachable fraction of brain tau = F_EXT x PHI_ACC = %.2e\n",
                at(a, TENROL, "F_EXT")*input$phiacc))
    cat("An extracellular-only mechanism cannot exceed this at any dose or affinity.\n")
  })

  ## ---- 3. wave -------------------------------------------------------------
  output$p_wave <- renderPlot({
    a <- simA()
    lp(a, paste0("A", 1:8), "regional aggregate load (0-1)",
       "the front moving through eight nuclei",
       labels = REGION, cols_pal = RCOL, legpos = "bottomright")
  })
  output$p_neur <- renderPlot({
    a <- simA()
    lp(a, paste0("N", 1:8), "surviving neuron fraction",
       "the integrator: regional neuronal loss",
       labels = REGION, cols_pal = RCOL, legpos = "bottomleft")
  })
  output$p_order <- renderPlot({
    a <- simA()
    tt <- vapply(1:8, function(i) {
      j <- which(a[[paste0("A", i)]] >= 0.5)
      if (length(j)) a$time[j[1]]/Y else NA_real_
    }, numeric(1))
    o <- order(tt, na.last = TRUE)
    par(mar = c(4, 12, 3, 2))
    barplot(rev(tt[o]), horiz = TRUE, names.arg = rev(REGION[o]), las = 1,
            col = rev(RCOL[o]), xlab = "years to 50% regional aggregate load",
            main = "recruitment order -- the wave's itinerary", cex.names = 0.85)
    abline(v = TENROL/Y, lty = 2, col = "grey20")
    abline(v = TSYMPT/Y, lty = 3, col = "grey50")
  })

  ## ---- 4. PK ---------------------------------------------------------------
  output$p_pkmab <- renderPlot({
    a <- simA()
    if (max(a$MABC) <= 0) { plot.new(); title("no mAb in this scenario"); return() }
    lp(a, c("MABC", "AB_ISF"), "mg (central) / nM (ISF)",
       "anti-tau mAb: plasma vs ISF", log = "y",
       labels = c("central compartment (mg)", "ISF antibody (nM)"))
  })
  output$p_pkaso <- renderPlot({
    a <- simA()
    if (max(a$ASOT) <= 0) { plot.new(); title("no ASO in this scenario"); return() }
    lp(a, c("ASOC", "ASOT", "KNOCK"), "mg / %",
       "MAPT ASO: CSF, tissue, knockdown",
       labels = c("CSF ASO (mg)", "brain tissue ASO (mg)", "monomer knockdown (%)"))
  })
  output$p_pksm <- renderPlot({
    a <- simA()
    cand <- c(TIDC = "tideglusib", OGAC = "OGA inhibitor", DAVC = "davunetide",
              SALC = "salicylate", FASC = "fasudil", TPNC = "TPN-101",
              LMC = "LM11A-31", EZEC = "ezeprogind", RILC = "riluzole",
              LDPC = "levodopa", ZOLC = "zolpidem", VTIT = "AADvac1 titre")
    keep <- names(cand)[vapply(names(cand), function(k) max(a[[k]]) > 0, logical(1))]
    if (!length(keep)) { plot.new(); title("no small molecule in this scenario"); return() }
    lp(a, keep, "amount (mg) / titre", "small molecules and vaccine",
       labels = unname(cand[keep]))
  })

  ## ---- 5. engagement -------------------------------------------------------
  output$p_eng <- renderPlot(two("CSF_NUB", "nM", "CSF UNBOUND N-terminal tau (the -98%)", log = "y"))
  output$p_updown <- renderPlot({
    a <- simA(); b <- simB(); x <- a$time/Y
    r1 <- a$CSF_NUB/b$CSF_NUB; r2 <- a$CSF_TT/b$CSF_TT; r3 <- a$CSF_MTBR/b$CSF_MTBR
    ym <- cbind(r1, r2, r3); ym[!is.finite(ym) | ym <= 0] <- NA
    plot(NA, xlim = range(x), ylim = range(ym, na.rm = TRUE), log = "y", las = 1,
         xlab = "years from first tau aggregate", ylab = "ratio to comparator",
         main = "three tau assays, three answers")
    grid(col = "grey92"); abline(h = 1, col = "grey40")
    lines(x, r1, col = "#b2182b", lwd = 2.4)
    lines(x, r2, col = "#2166ac", lwd = 2.4)
    lines(x, r3, col = "#1b7837", lwd = 2.4)
    vline_marks()
    legend("left", bty = "n", cex = 0.8, lwd = 2.4,
           col = c("#b2182b", "#2166ac", "#1b7837"),
           legend = c("unbound N-terminal tau (DOWN)",
                      "TOTAL tau (UP -- sink effect)",
                      "MTBR-tau243 (the causal species)"))
  })
  output$t_eng <- renderPrint({
    a <- simA(); b <- simB(); tt <- TENROL + 364
    f <- function(col) 100*(at(a, tt, col)/at(b, tt, col) - 1)
    cat(sprintf("ISF antibody at steady state           %.2f nM\n", at(a, tt, "AB_ISF")))
    cat(sprintf("CSF unbound N-terminal tau            %+.2f %%\n", f("CSF_NUB")))
    cat(sprintf("CSF total tau                         %+.1f %%\n", f("CSF_TT")))
    cat(sprintf("CSF MTBR-tau243 (causal species)      %+.2f %%\n", f("CSF_MTBR")))
    cat(sprintf("CSF seeding activity                  %+.2f %%\n", f("CSF_SEED")))
    cat(sprintf("PSPRS at 52 weeks                     %+.2f points vs comparator\n",
                at(a, tt, "PSPRS") - at(b, tt, "PSPRS")))
  })

  ## ---- 6. endpoints --------------------------------------------------------
  output$p_psprs <- renderPlot(two("PSPRS", "PSPRS", "PSPRS -- the primary endpoint"))
  output$p_vsv <- renderPlot(two("VSV", "deg/s",
                                 "vertical saccade peak velocity (floors early)"))
  output$p_domains <- renderPlot({
    a <- simA(); x <- a$time/Y
    plot(NA, xlim = range(x), ylim = c(0, 100), las = 1,
         xlab = "years from first tau aggregate", ylab = "score",
         main = "PSPRS vs mPSPRS-10 vs SEADL")
    grid(col = "grey92")
    lines(x, a$PSPRS, col = "#b2182b", lwd = 2.4)
    lines(x, a$MP10,  col = "#d95f02", lwd = 2, lty = 2)
    lines(x, a$SEADL, col = "#2166ac", lwd = 2)
    vline_marks()
    legend("left", bty = "n", cex = 0.8, lwd = 2,
           col = c("#b2182b", "#d95f02", "#2166ac"),
           legend = c("PSPRS", "mPSPRS-10", "SEADL (%)"))
  })
  output$t_endpoint <- renderPrint({
    a <- simA(); b <- simB()
    d1 <- at(a, TENROL + 364, "PSPRS") - at(a, TENROL, "PSPRS")
    d0 <- at(b, TENROL + 364, "PSPRS") - at(b, TENROL, "PSPRS")
    cat("52-week PSPRS change\n")
    cat(sprintf("  comparator %s : %+.2f points\n", input$ref, d0))
    cat(sprintf("  scenario   %s : %+.2f points\n", input$scen, d1))
    cat(sprintf("  difference        : %+.2f points\n\n", d1 - d0))
    cat("absolute PSPRS at enrolment + 52 weeks\n")
    cat(sprintf("  comparator %.2f | scenario %.2f | difference %+.2f\n",
                at(b, TENROL + 364, "PSPRS"), at(a, TENROL + 364, "PSPRS"),
                at(a, TENROL + 364, "PSPRS") - at(b, TENROL + 364, "PSPRS")))
    cat("\nA disease-DELAYING drug can look worse on change-from-baseline and\n")
    cat("better on the absolute score. Both numbers are shown for that reason.\n")
  })

  ## ---- 7. fluid biomarkers -------------------------------------------------
  output$p_nfl <- renderPlot({
    a <- simA(); b <- simB(); x <- a$time/Y
    plot(NA, xlim = range(x), ylim = range(c(a$NFLC, b$NFLC)), las = 1,
         xlab = "years from first tau aggregate", ylab = "CSF NfL (pg/mL)",
         main = "CSF neurofilament light")
    grid(col = "grey92")
    lines(x, b$NFLC, col = "grey45", lwd = 2, lty = 2)
    lines(x, a$NFLC, col = "#b2182b", lwd = 2.4)
    vline_marks(); mark_legend()
  })
  output$p_tauassays <- renderPlot({
    a <- simA()
    lp(a, c("CSF_NUB", "CSF_TT", "CSF_P181", "CSF_MTBR", "CSF_SEED"),
       "nM (log scale)", "CSF tau assays", log = "y",
       labels = c("unbound N-terminal", "total tau", "p-tau181",
                  "MTBR-tau243", "seeding activity"), legpos = "bottomright")
  })
  output$p_gfap <- renderPlot({
    a <- simA()
    lp(a, c("NFLP", "GFAPP"), "pg/mL", "plasma NfL and GFAP",
       labels = c("plasma NfL", "plasma GFAP"))
  })
  output$t_bm <- renderPrint({
    a <- simA(); b <- simB()
    r <- function(col, t0, t1) 100*(at(a, t1, col)/at(a, t0, col) - 1)
    cat("Change over the 52 weeks after enrolment (this scenario)\n")
    for (k in c("NFLC", "NFLP", "GFAPP", "CSF_TT", "CSF_MTBR", "CSF_SEED"))
      cat(sprintf("  %-9s %+7.1f %%\n", k, r(k, TENROL, TENROL + 364)))
    cat("\nComparator NfL change over the same window: ")
    cat(sprintf("%+.1f %%\n", 100*(at(b, TENROL + 364, "NFLC")/at(b, TENROL, "NFLC") - 1)))
  })

  ## ---- 8. imaging ----------------------------------------------------------
  output$p_mid <- renderPlot(two("MIDA", "mm2", "midbrain area (hummingbird sign)"))
  output$p_mrpi <- renderPlot({
    a <- simA(); b <- simB(); x <- a$time/Y
    plot(NA, xlim = range(x), ylim = range(c(a$MRPI, b$MRPI)), las = 1,
         xlab = "years from first tau aggregate", ylab = "MRPI",
         main = "MR parkinsonism index")
    grid(col = "grey92"); abline(h = 13.55, col = "#b2182b", lty = 3)
    text(min(x), 13.55, "diagnostic cut-off 13.55", pos = 4, cex = 0.75, col = "#b2182b")
    lines(x, b$MRPI, col = "grey45", lwd = 2, lty = 2)
    lines(x, a$MRPI, col = "#2166ac", lwd = 2.4)
    vline_marks()
  })

  ## ---- 9. comparison -------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    withProgress(message = "running all scenarios", value = 0.1, {
      psp_compare(horizon = round(input$chz*7), verbose = FALSE)
    })
  })
  output$t_compare <- renderTable(cmp(), digits = 3)

  ## ---- 10. when to treat ---------------------------------------------------
  el <- eventReactive(input$go_el, {
    withProgress(message = "computing elasticity", value = 0.1, {
      psp_elasticity(verbose = FALSE)
    })
  })
  output$p_el1 <- renderPlot({
    e <- el()$table
    plot(e$knockdown_pct, e$tx_rel, type = "b", pch = 19, col = "#b2182b",
         lwd = 2, las = 1, ylim = c(0, 1.05),
         xlab = "MAPT knockdown (%)", ylab = "slope / placebo slope",
         main = "started AT ENROLMENT (what a trial measures)")
    grid(col = "grey92"); abline(h = 0.5, lty = 3, col = "grey40")
    text(45, 0.55, "halving the slope -- unreachable", cex = 0.8, col = "grey30")
  })
  output$p_el2 <- renderPlot({
    e <- el()$table
    plot(e$knockdown_pct, e$prev_delay_y, type = "b", pch = 19, col = "#1b7837",
         lwd = 2, las = 1, xlab = "MAPT knockdown (%)",
         ylab = "delay to trial-entry severity (years)",
         main = "started at BIOLOGICAL ONSET")
    grid(col = "grey92")
  })
  output$t_el <- renderTable(el()$table, digits = 3)

  ## ---- 11. survival --------------------------------------------------------
  output$p_surv <- renderPlot({
    a <- simA(); b <- simB(); x <- a$time/Y
    plot(NA, xlim = range(x), ylim = c(0, 1), las = 1,
         xlab = "years from first tau aggregate", ylab = "S(t)",
         main = "overall survival")
    grid(col = "grey92"); abline(h = 0.5, lty = 3, col = "grey40")
    lines(x, b$SURV, col = "grey45", lwd = 2, lty = 2)
    lines(x, a$SURV, col = "#b2182b", lwd = 2.4)
    vline_marks(); mark_legend()
  })
  output$p_cif <- renderPlot({
    a <- simA(); x <- a$time/Y
    plot(NA, xlim = range(x), ylim = c(0, 1), las = 1,
         xlab = "years from first tau aggregate", ylab = "cumulative incidence",
         main = "competing risks: dCIF/dt = h * S")
    grid(col = "grey92")
    lines(x, a$CIFP, col = "#b2182b", lwd = 2.4)
    lines(x, a$CIFO, col = "#2166ac", lwd = 2.4)
    lines(x, a$DYS/4, col = "#d95f02", lwd = 2, lty = 2)
    vline_marks()
    legend("topleft", bty = "n", cex = 0.8, lwd = 2,
           col = c("#b2182b", "#2166ac", "#d95f02"),
           legend = c("aspiration pneumonia death", "other-cause death",
                      "dysphagia score / 4"))
  })
  output$t_surv <- renderPrint({
    a <- simA(); b <- simB()
    ms <- function(d) { i <- which(d$SURV <= 0.5); if (length(i)) d$time[i[1]]/Y else NA }
    cat(sprintf("median survival, comparator %s : %.2f y from tau onset (%.2f y from symptoms)\n",
                input$ref, ms(b), ms(b) - TSYMPT/Y))
    cat(sprintf("median survival, scenario   %s : %.2f y from tau onset (%.2f y from symptoms)\n",
                input$scen, ms(a), ms(a) - TSYMPT/Y))
    cat(sprintf("\ncause-specific share at end of simulation: pneumonia %.0f%%, other %.0f%%\n",
                100*tail(a$CIFP, 1)/(tail(a$CIFP, 1) + tail(a$CIFO, 1)),
                100*tail(a$CIFO, 1)/(tail(a$CIFP, 1) + tail(a$CIFO, 1))))
  })

  ## ---- 12. validation ------------------------------------------------------
  val <- eventReactive(input$go_val, {
    withProgress(message = "validating", value = 0.1, psp_validate(verbose = FALSE))
  })
  output$t_val <- renderTable(val(), digits = 4)
}

shinyApp(ui, server)
