# =====================================================================
# Retinopathy of Prematurity (ROP) — interactive QSP dashboard
#   Companion to rop_mrgsolve_model.R
#
#   Eight tabs, organised around the model's structure rather than
#   around the drug list:
#
#     1  Patient & oxygen policy  — who this infant is, what the
#        ventilator is doing, and what the oximeter is therefore showing
#     2  Retinal oxygen           — the SUPPLY/DEMAND collision:
#        PCHOR against PCRIT, the crossing date, the diffusion depth
#     3  VEGF pools               — front pool vs vitreous pool, and the
#        three gates (growth KG, survival KSURV, neovascular KNV)
#     4  Disease course           — ICROP zone/stage/plus and the ETROP
#        type 1 flag over postmenstrual age
#     5  Anti-VEGF PK/PD          — vitreous and serum drug, free VEGF
#        in both compartments, ocular/systemic asymmetry
#     6  Dose ladder              — PEDIG de-escalation: what each dose
#        buys ocularly and costs systemically
#     7  Treatment comparison     — laser vs bevacizumab vs ranibizumab
#        vs aflibercept in one treatment-requiring eye
#     8  Cohort & trials          — virtual cohort against NeOProM, and
#        the achieved-vs-prescribed separation sweep
#
#   Run:  shiny::runApp("rop_shiny_app.R")
# =====================================================================

library(shiny)
suppressMessages(library(mrgsolve))
suppressMessages(library(dplyr))

# Source the model without triggering its scenario run.
local({
  f <- "rop_mrgsolve_model.R"
  if (!file.exists(f)) f <- file.path("retinopathy-of-prematurity", f)
  if (!file.exists(f)) stop("rop_mrgsolve_model.R not found")
  source(f, local = FALSE)
})

PAL <- list(ox = "#2f6fb5", dem = "#c0504d", front = "#4c9a4c",
            vit = "#c0504d", drug = "#7b52a8", grey = "#8a8a8a",
            gold = "#d9a419", ink = "#2b2b2b")

pnl <- function(...) tags$div(style = "padding:6px 2px", ...)
note <- function(...) tags$p(style = "color:#555;font-size:12px;margin:4px 0", ...)

# ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{font-family:Helvetica,Arial,sans-serif}
    h4{margin-top:14px;color:#2b2b2b}
    .well{background:#fafafa}
    table{font-size:12px}
  "))),
  titlePanel("Retinopathy of Prematurity — QSP model explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Infant"),
      sliderInput("GA", "Gestational age (wk)", 23, 28, 25, 0.5),
      sliderInput("BW", "Birth weight (g)", 400, 1400, 700, 25),
      sliderInput("NUTR", "Nutritional adequacy (IGF-1 drive)", 0.4, 1.4, 1.0, 0.05),
      h4("Oxygen policy"),
      sliderInput("SPO2TGT", "ACHIEVED SpO2 midpoint (%)", 84, 99, 93.5, 0.5),
      note("This is the achieved midpoint, not the prescribed range. ",
           "Reproducing NeOProM needs an achieved separation of ~2 points."),
      sliderInput("PMAWEAN", "PMA weaned to room air (wk)", 32, 44, 36, 0.5),
      checkboxInput("O2SUPP", "Supplemental O2 from a set time (STOP-ROP design)", FALSE),
      conditionalPanel("input.O2SUPP",
        sliderInput("SPO2SUP", "Supplemental SpO2 target (%)", 92, 99.4, 97.5, 0.5),
        sliderInput("O2SUPPPMA", "Start at PMA (wk)", 30, 40, 34.5, 0.5)),
      h4("Exposures"),
      checkboxInput("TX", "Red-cell transfusion(s)", FALSE),
      conditionalPanel("input.TX",
        sliderInput("FTX", "HbF fraction replaced per unit", 0.05, 0.5, 0.20, 0.05)),
      checkboxInput("SEPSIS", "Sepsis episode (day 20)", FALSE),
      checkboxInput("AADHA", "AA:DHA enteral supplement", FALSE),
      h4("Treatment"),
      selectInput("TRT", "Intervention",
                  c("none", "laser", "bevacizumab", "ranibizumab", "aflibercept")),
      conditionalPanel("input.TRT != 'none'",
        sliderInput("TRTPMA", "Treat at PMA (wk)", 31, 42, 35, 0.5)),
      conditionalPanel("input.TRT != 'none' && input.TRT != 'laser'",
        selectInput("DOSE", "Dose (mg)",
                    c("0.625", "0.25", "0.125", "0.063", "0.031",
                      "0.016", "0.008", "0.004", "0.002", "0.2", "0.1", "0.4"),
                    selected = "0.625")),
      conditionalPanel("input.TRT == 'laser'",
        sliderInput("LFRAC", "Fraction of avascular retina ablated", 0.3, 0.98, 0.85, 0.01)),
      h4("Model detail"),
      sliderInput("END", "Follow-up (days)", 100, 300, 200, 10)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & oxygen",
          pnl(h4("Oxygen carriage: the oximeter reading is not the oxygen tension"),
              plotOutput("p_ox", height = "330px"),
              note("Left: SpO2 and the P50 that the haemoglobin switch and any ",
                   "transfusion impose. Right: the resulting PaO2. The two panels ",
                   "can move in opposite directions."),
              h4("Same SpO2, two haemoglobins"),
              tableOutput("t_oximetry"),
              note("dPaO2/dSpO2 rises from ~1.3 mmHg per point at 87% to ~45 mmHg ",
                   "per point at 99%: above 95% the instrument cannot resolve the ",
                   "oxygen tension that decides the therapy in tab 2."))),
        tabPanel("2 · Retinal oxygen",
          pnl(h4("Supply against demand in the avascular retina"),
              plotOutput("p_krogh", height = "340px"),
              note("PCRIT = MO2·h²/2k is the choroidal pO2 required to oxygenate ",
                   "the full retinal thickness. It is set by maturation and no drug ",
                   "moves it. Where PCHOR crosses PCRIT the disease changes sign."),
              h4("Oxygen requirement by postmenstrual age"),
              tableOutput("t_krogh"),
              note("The square-root law: doubling the oxygenated depth needs four ",
                   "times the choroidal pO2."))),
        tabPanel("3 · VEGF pools",
          pnl(h4("Two pools, two gates, opposite oxygen sensitivities"),
              plotOutput("p_vegf", height = "340px"),
              h4("Vascularised and avascular retinal fraction"),
              plotOutput("p_vasc", height = "260px"),
              note("The front pool sets growth (above KG) and endothelial survival ",
                   "(above KSURV). The vitreous pool sets neovascularisation (above ",
                   "KNV). Anti-VEGF suppresses both, because they share a compartment."))),
        tabPanel("4 · Disease course",
          pnl(h4("ICROP stage, zone and plus disease"),
              plotOutput("p_stage", height = "340px"),
              h4("Milestones"),
              tableOutput("t_course"),
              note("ETROP type 1 = zone I any stage with plus, zone I stage 3, or ",
                   "zone II stage 2-3 with plus."))),
        tabPanel("5 · Anti-VEGF PK/PD",
          pnl(h4("Vitreous and plasma exposure after one intravitreal dose"),
              plotOutput("p_pk", height = "340px"),
              h4("Free VEGF in both compartments"),
              plotOutput("p_pd", height = "260px"),
              note("Bevacizumab keeps an Fc, so FcRn recycling gives it a serum ",
                   "half-life of ~21 days in infants; ranibizumab is a Fab and has ",
                   "none. Measured free serum VEGF falls far less than 1:1 binding ",
                   "predicts because drug-bound VEGF is protected from clearance."),
              h4("Ocular / systemic asymmetry"),
              tableOutput("t_asym"))),
        tabPanel("6 · Dose ladder",
          pnl(h4("PEDIG de-escalation: efficacy is logarithmic, exposure is linear"),
              plotOutput("p_ladder", height = "340px"),
              tableOutput("t_ladder"),
              note("Observed 4-week success: 0.031 mg 9/9, 0.016 mg 13/13, ",
                   "0.008 mg 9/9, 0.004 mg 9/10, 0.002 mg 17/23."))),
        tabPanel("7 · Treatment comparison",
          pnl(h4("One treatment-requiring eye, six managements"),
              actionButton("go_trt", "Run comparison"),
              tableOutput("t_trt"),
              plotOutput("p_trt", height = "320px"),
              note("Laser deletes the VEGF source and cannot reactivate in this ",
                   "model; anti-VEGF silences it and always can. The model does NOT ",
                   "reproduce BEAT-ROP's 22% laser recurrence — see F1 in the ",
                   "model header."))),
        tabPanel("8 · Cohort & trials",
          pnl(h4("Virtual cohort against NeOProM"),
              fluidRow(
                column(4, numericInput("N", "Cohort size", 150, 40, 800, 10)),
                column(4, numericInput("SEP", "Achieved separation (SpO2 points)",
                                       2.0, 0.5, 8, 0.5)),
                column(4, br(), actionButton("go_cohort", "Run cohort"))),
              tableOutput("t_cohort"),
              note("NeOProM observed: treatment-requiring ROP 14.9% (higher target) ",
                   "vs 10.9% (lower), RR 0.74; death 17.1% vs 19.9%, RR 1.17; ",
                   "severe NEC 6.9% vs 9.2%."),
              h4("The competing-risk exchange rate (trial arithmetic, no model)"),
              tableOutput("t_exch")))
      )
    )
  )
)

# ---------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- reactive({
    p <- list(GA = input$GA, BW = input$BW, NUTR = input$NUTR,
              SPO2TGT = input$SPO2TGT, PMAWEAN = input$PMAWEAN)
    if (isTRUE(input$O2SUPP)) {
      p$O2SUPP <- 1
      p$O2SUPPT <- (input$O2SUPPPMA - input$GA) * 7
      p$SPO2SUP <- input$SPO2SUP
    }
    if (isTRUE(input$TX)) { p$TXT1 <- 18; p$TXT2 <- 34; p$FTX <- input$FTX }
    if (isTRUE(input$SEPSIS)) p$SEPT <- 20
    if (isTRUE(input$AADHA)) p$KAADHA <- 1
    if (identical(input$TRT, "laser")) {
      p$LASERT <- (input$TRTPMA - input$GA) * 7
      p$LFRAC <- input$LFRAC
    }
    p
  })

  sim <- reactive({
    p <- pars()
    if (input$TRT %in% c("bevacizumab", "ranibizumab", "aflibercept")) {
      do.call(sim_rop, c(p, list(dose_mg = as.numeric(input$DOSE),
                                 drug = input$TRT,
                                 tdose = (input$TRTPMA - input$GA) * 7,
                                 end = input$END, delta = 0.25)))
    } else {
      do.call(sim_rop, c(p, list(end = input$END, delta = 0.25)))
    }
  })

  # ---------------- tab 1 -------------------------------------------
  output$p_ox <- renderPlot({
    d <- sim(); op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 4))
    plot(d$PMAo, d$spo2o, type = "l", lwd = 2.5, col = PAL$ox, ylim = c(82, 100),
         xlab = "PMA (weeks)", ylab = "SpO2 (%)", main = "Oximeter and P50")
    par(new = TRUE)
    plot(d$PMAo, d$P50o, type = "l", lwd = 2, lty = 2, col = PAL$dem,
         axes = FALSE, xlab = "", ylab = "", ylim = c(18, 27))
    axis(4); mtext("P50 (mmHg)", 4, line = 2.6, col = PAL$dem)
    legend("bottomright", c("SpO2", "P50"), col = c(PAL$ox, PAL$dem),
           lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    plot(d$PMAo, d$PAO2o, type = "l", lwd = 2.5, col = PAL$ox,
         xlab = "PMA (weeks)", ylab = "PaO2 (mmHg)",
         main = "Arterial oxygen tension")
    abline(h = 30, lty = 3, col = PAL$grey)
    text(max(d$PMAo), 32, "intrauterine", adj = 1, cex = 0.75, col = PAL$grey)
    par(op)
  })

  output$t_oximetry <- renderTable({
    x <- scenario_oximetry()
    names(x) <- c("SpO2 (%)", "PaO2 if HbF (mmHg)", "PaO2 if HbA (mmHg)",
                  "ratio", "dPaO2 per SpO2 point")
    x
  }, digits = 2)

  # ---------------- tab 2 -------------------------------------------
  output$p_krogh <- renderPlot({
    d <- sim(); op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 1))
    ylim <- range(c(d$PCHORo, d$PCRITo), na.rm = TRUE)
    plot(d$PMAo, d$PCHORo, type = "l", lwd = 2.5, col = PAL$ox, ylim = ylim,
         xlab = "PMA (weeks)", ylab = "mmHg", main = "Supply vs demand")
    lines(d$PMAo, d$PCRITo, lwd = 2.5, col = PAL$dem)
    lines(d$PMAo, d$PINNo, lwd = 2, col = PAL$gold)
    j <- which(diff(sign(d$PCHORo - d$PCRITo)) != 0)
    if (length(j)) {
      abline(v = d$PMAo[j[1]], lty = 2, col = PAL$ink)
      text(d$PMAo[j[1]], ylim[2], sprintf(" switch %.1f wk", d$PMAo[j[1]]),
           adj = 0, cex = 0.8)
    }
    legend("topleft", c("PCHOR (supply)", "PCRIT (demand)", "PINNER (inner pO2)"),
           col = c(PAL$ox, PAL$dem, PAL$gold), lwd = 2, bty = "n", cex = 0.8)
    plot(d$PMAo, d$LPENo, type = "l", lwd = 2.5, col = PAL$ox,
         xlab = "PMA (weeks)", ylab = "microns",
         main = "Oxygen penetration depth vs retinal thickness")
    lines(d$PMAo, 1e4 * (0.0105 + (0.0185 - 0.0105) /
          (1 + exp(-(d$PMAo - 33) / 2.5))), lwd = 2.5, col = PAL$dem, lty = 2)
    legend("topright", c("penetration L", "retinal thickness h"),
           col = c(PAL$ox, PAL$dem), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
    par(op)
  })

  output$t_krogh <- renderTable({
    x <- scenario_diffusion_ceiling()
    x$MO2 <- x$MO2 * 1e4
    names(x) <- c("PMA (wk)", "MO2 (x1e-4)", "h (um)", "PCRIT (mmHg)",
                  "PaO2 required", "SpO2 required if HbF", "SpO2 required if HbA",
                  "L at PaO2 60 (um)", "L at PaO2 110 (um)")
    x
  }, digits = 2)

  # ---------------- tab 3 -------------------------------------------
  output$p_vegf <- renderPlot({
    d <- sim(); op <- par(mar = c(4, 4.2, 2, 1))
    ymax <- max(c(d$VFRo, d$VEGFFo, 1400), na.rm = TRUE)
    plot(d$PMAo, d$VFRo, type = "l", lwd = 2.5, col = PAL$front, ylim = c(0, ymax),
         xlab = "PMA (weeks)", ylab = "free VEGF (pg/mL)",
         main = "Front pool (green) and vitreous pool (red)")
    lines(d$PMAo, d$VEGFFo, lwd = 2.5, col = PAL$vit)
    abline(h = 300, lty = 2, col = PAL$front)
    abline(h = 180, lty = 3, col = PAL$front)
    abline(h = 1170, lty = 2, col = PAL$vit)
    text(min(d$PMAo), 300, " KG growth gate", adj = 0, cex = 0.75, col = PAL$front)
    text(min(d$PMAo), 180, " KSURV survival gate", adj = 0, cex = 0.75, col = PAL$front)
    text(min(d$PMAo), 1170, " KNV neovascular gate", adj = 0, cex = 0.75, col = PAL$vit)
    par(op)
  })

  output$p_vasc <- renderPlot({
    d <- sim(); op <- par(mar = c(4, 4.2, 2, 1))
    plot(d$PMAo, d$VASC, type = "l", lwd = 2.5, col = PAL$front, ylim = c(0, 1),
         xlab = "PMA (weeks)", ylab = "fraction of retina",
         main = "Vascularised, avascular and ablated retina")
    lines(d$PMAo, d$AVASCo, lwd = 2.5, col = PAL$vit)
    lines(d$PMAo, d$ABLA, lwd = 2, col = PAL$drug, lty = 2)
    legend("right", c("vascularised", "avascular", "ablated"),
           col = c(PAL$front, PAL$vit, PAL$drug), lty = c(1, 1, 2),
           lwd = 2, bty = "n", cex = 0.85)
    par(op)
  })

  # ---------------- tab 4 -------------------------------------------
  output$p_stage <- renderPlot({
    d <- sim(); op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 1))
    plot(d$PMAo, d$STAGE, type = "s", lwd = 2.5, col = PAL$vit, ylim = c(0, 5),
         xlab = "PMA (weeks)", ylab = "ICROP stage", main = "Stage and zone")
    lines(d$PMAo, d$ZONE, type = "s", lwd = 2, col = PAL$ox, lty = 2)
    legend("topright", c("stage", "zone"), col = c(PAL$vit, PAL$ox),
           lty = c(1, 2), lwd = 2, bty = "n", cex = 0.85)
    plot(d$PMAo, d$PLUS, type = "l", lwd = 2.5, col = PAL$gold, ylim = c(0, 1),
         xlab = "PMA (weeks)", ylab = "index",
         main = "Plus disease, NV and detachment")
    lines(d$PMAo, d$NV, lwd = 2.5, col = PAL$vit)
    lines(d$PMAo, d$DETACH, lwd = 2, col = PAL$ink, lty = 2)
    abline(h = 0.5, lty = 3, col = PAL$gold)
    if (any(d$TYPE1 > 0.5)) abline(v = d$PMAo[which(d$TYPE1 > 0.5)[1]],
                                   lty = 2, col = "red")
    legend("topright", c("plus", "NV", "detachment"),
           col = c(PAL$gold, PAL$vit, PAL$ink), lty = c(1, 1, 2),
           lwd = 2, bty = "n", cex = 0.85)
    par(op)
  })

  output$t_course <- renderTable({
    d <- sim()
    first <- function(v) if (any(v)) sprintf("%.1f", d$PMAo[which(v)[1]]) else "never"
    data.frame(
      milestone = c("stage 1", "stage 2", "stage 3", "plus disease",
                    "ETROP type 1", "threshold surrogate", "phase 2 switch",
                    "peak vitreous VEGF (pg/mL)", "peak NV index",
                    "final spherical equivalent (D)", "field loss",
                    "cumulative death probability"),
      value = c(first(d$STAGE >= 1), first(d$STAGE >= 2), first(d$STAGE >= 3),
                first(d$PLUSD > 0.5), first(d$TYPE1 > 0.5), first(d$THRESH > 0.5),
                first(d$PINNo <= 0),
                sprintf("%.0f", max(d$VEGFFo)), sprintf("%.3f", max(d$NV)),
                sprintf("%.2f", tail(d$MYOP, 1)), sprintf("%.3f", tail(d$VFLOSS, 1)),
                sprintf("%.3f", tail(d$PDEATH, 1))),
      stringsAsFactors = FALSE)
  })

  # ---------------- tab 5 -------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim(); op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 1))
    if (max(d$CVITo) <= 0) {
      plot.new(); title("select an anti-VEGF agent in the sidebar")
      plot.new(); par(op); return(invisible())
    }
    plot(d$PMAo, pmax(d$CVITo, 1e-4), type = "l", lwd = 2.5, col = PAL$drug,
         log = "y", xlab = "PMA (weeks)", ylab = "free drug (nM)",
         main = "Vitreous")
    abline(h = 0.3, lty = 2, col = PAL$grey)
    text(min(d$PMAo), 0.3, " Kd", adj = 0, cex = 0.75)
    plot(d$PMAo, pmax(d$CSYSNG, 1e-3), type = "l", lwd = 2.5, col = PAL$drug,
         log = "y", xlab = "PMA (weeks)", ylab = "free drug (ng/mL)",
         main = "Plasma")
    par(op)
  })

  output$p_pd <- renderPlot({
    d <- sim(); op <- par(mar = c(4, 4.2, 2, 4))
    plot(d$PMAo, d$VEGFFo, type = "l", lwd = 2.5, col = PAL$vit,
         xlab = "PMA (weeks)", ylab = "free vitreous VEGF (pg/mL)",
         main = "Free VEGF: vitreous (red) and serum (grey)")
    par(new = TRUE)
    plot(d$PMAo, d$VEGFSFo, type = "l", lwd = 2.5, col = PAL$grey,
         axes = FALSE, xlab = "", ylab = "")
    axis(4); mtext("free serum VEGF (pg/mL)", 4, line = 2.6, col = PAL$grey)
    par(op)
  })

  output$t_asym <- renderTable({
    data.frame(
      quantity = c("preterm vitreous volume (mL)", "adult vitreous volume (mL)",
                   "ocular concentration, preterm / adult",
                   "systemic dose per kg, preterm / adult",
                   "bevacizumab 0.625 mg (nmol)", "ranibizumab 0.2 mg (nmol)",
                   "molar ratio of the two standard doses",
                   "binding-site capacity ratio (bev / ranib)",
                   "vitreous binding sites at 0.625 mg (nM)",
                   "vitreous VEGF at 1500 pg/mL (nM)",
                   "molar excess of drug over target"),
      value = c("1.10", "4.00", "1.8x", "35x",
                sprintf("%.3f", mg_to_nmol(0.625, 149)),
                sprintf("%.3f", mg_to_nmol(0.2, 48)),
                sprintf("%.3f", mg_to_nmol(0.625, 149) / mg_to_nmol(0.2, 48)),
                "2.01", "3813", "0.033", "~229,000x"),
      stringsAsFactors = FALSE)
  })

  # ---------------- tab 6 -------------------------------------------
  ladder <- reactive(scenario_dose_ladder())

  output$p_ladder <- renderPlot({
    x <- ladder(); op <- par(mar = c(4.2, 4.4, 2, 4.4))
    plot(x$dose_mg, x$vitreous_nM_d28, log = "xy", type = "b", pch = 19,
         lwd = 2, col = PAL$drug, xlab = "bevacizumab dose (mg)",
         ylab = "vitreous free drug at 4 weeks (nM)",
         main = "Ocular coverage (log) vs systemic exposure (linear)")
    abline(h = c(0.058, 1.1), lty = c(2, 3), col = PAL$grey)
    text(max(x$dose_mg), 0.058, "Kd range", adj = 1, cex = 0.75, col = PAL$grey)
    par(new = TRUE)
    plot(x$dose_mg, x$serumVEGF_pct_d14, log = "x", type = "b", pch = 17,
         lwd = 2, col = PAL$dem, axes = FALSE, xlab = "", ylab = "",
         ylim = c(0, 105))
    axis(4); mtext("free serum VEGF at day 14 (% of baseline)", 4,
                   line = 2.8, col = PAL$dem)
    legend("bottomright", c("vitreous trough", "serum VEGF"),
           col = c(PAL$drug, PAL$dem), pch = c(19, 17), lwd = 2,
           bty = "n", cex = 0.8)
    par(op)
  })

  output$t_ladder <- renderTable({
    x <- ladder()
    obs <- c("-", "11/11", "14/14", "21/24", "9/9", "13/13", "9/9", "9/10", "17/23")
    data.frame(`dose (mg)` = x$dose_mg,
               `sites (pmol)` = round(x$sites_pmol, 1),
               `vitreous nM d1` = round(x$vitreous_nM_d1, 1),
               `vitreous nM d28` = round(x$vitreous_nM_d28, 3),
               `serum ng/mL d14` = round(x$serum_ngmL_d14, 1),
               `serum VEGF % d14` = round(x$serumVEGF_pct_d14, 1),
               `days VEGF < 50%` = round(x$days_VEGF_below_half, 1),
               `PEDIG 4-wk success` = obs,
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  # ---------------- tab 7 -------------------------------------------
  trt <- eventReactive(input$go_trt, {
    scenario_treatment(GA = input$GA, NUTR = min(input$NUTR, 0.8))
  })

  output$t_trt <- renderTable({
    x <- trt()
    if (!"arm" %in% names(x) || nrow(x) == 0) return(x)
    data.frame(arm = x$arm,
               `treat PMA` = round(x$treat_PMA, 1),
               reactivated = ifelse(x$reactivation > 0.5, "yes", "no"),
               `weeks to reactivation` = round(x$react_wk_after, 1),
               `weeks VEGF < KG` = round(x$wk_VEGF_below_KG, 1),
               `residual avascular` = round(x$residual_avascular, 3),
               `spherical equiv (D)` = round(x$myopia_D, 2),
               `field loss` = round(x$field_loss, 3),
               `max detachment` = round(x$max_detach, 3),
               check.names = FALSE, stringsAsFactors = FALSE)
  })

  output$p_trt <- renderPlot({
    x <- trt(); op <- par(mfrow = c(1, 2), mar = c(7, 4.4, 2, 1))
    if (!"residual_avascular" %in% names(x)) { plot.new(); plot.new(); par(op); return() }
    barplot(x$residual_avascular, names.arg = x$arm, las = 2, col = PAL$vit,
            ylab = "residual avascular retina", main = "Persistent avascular retina")
    barplot(-x$myopia_D, names.arg = x$arm, las = 2, col = PAL$drug,
            ylab = "myopia (−D)", main = "Refractive cost")
    par(op)
  })

  # ---------------- tab 8 -------------------------------------------
  coh <- eventReactive(input$go_cohort, {
    hi <- 93.5
    cohort_rop(n = input$N, targets = c(hi - input$SEP, hi))
  })

  output$t_cohort <- renderTable({
    x <- coh()
    x <- x[order(x$SpO2_target), ]
    rr <- function(v) if (length(v) == 2 && v[2] != 0) v[1] / v[2] else NA
    out <- data.frame(
      `SpO2 target` = x$SpO2_target, n = x$n,
      `treated ROP %` = round(x$treated_ROP_pct, 2),
      `plus %` = round(x$plus_pct, 2),
      `stage 3+ %` = round(x$stage3plus_pct, 2),
      `death %` = round(x$death_pct, 2),
      check.names = FALSE)
    rbind(out, data.frame(`SpO2 target` = "risk ratio (low/high)", n = NA,
      `treated ROP %` = round(rr(x$treated_ROP_pct), 3),
      `plus %` = round(rr(x$plus_pct), 3),
      `stage 3+ %` = round(rr(x$stage3plus_pct), 3),
      `death %` = round(rr(x$death_pct), 3), check.names = FALSE))
  })

  output$t_exch <- renderTable({
    d_lo <- 0.199; d_hi <- 0.171; r_lo <- 0.109; r_hi <- 0.149
    ed <- 1000 * (d_lo - d_hi); fr <- 1000 * (r_hi - r_lo)
    data.frame(
      quantity = c("extra deaths per 1000", "fewer ROP treatments per 1000",
                   "extra severe NEC per 1000",
                   "unfavourable structural outcomes avoided (ETROP 9.1%)",
                   "unfavourable acuity outcomes avoided (ETROP 14.5%)",
                   "deaths per structural outcome avoided",
                   "deaths per acuity outcome avoided"),
      value = c(sprintf("%.1f", ed), sprintf("%.1f", fr),
                sprintf("%.1f", 1000 * (0.092 - 0.069)),
                sprintf("%.2f", fr * 0.091), sprintf("%.2f", fr * 0.145),
                sprintf("%.1f", ed / (fr * 0.091)),
                sprintf("%.1f", ed / (fr * 0.145))),
      stringsAsFactors = FALSE)
  })
}

shinyApp(ui, server)
