# =====================================================================
#  Hereditary Haemorrhagic Telangiectasia (HHT) — QSP explorer
#  Shiny front end for hht_mrgsolve_model.R
# ---------------------------------------------------------------------
#  The app is organised around the one identity the model turns on:
#
#        Hb_ss = A_net / (0.0347 * B)
#
#  Every tab is there to make one consequence of it visible and
#  manipulable.  In particular the IRON POLICY control is not a
#  convenience: switching it between "demand-driven" (what clinicians
#  actually do) and "protocol-fixed" (what a trial should do) is what
#  reproduces or removes the PATH-HHT haemoglobin paradox, and the app
#  exists largely so a user can flip that switch and watch it happen.
#
#  Tabs
#    1  Patient          genotype, bleeding, iron policy, baseline state
#    2  Drug PK          concentrations, target occupancy, molar excess
#    3  Bleeding & lesion nasal / GI loss, lesion burden, mural coverage
#    4  Iron & Hb        the mass balance and the hyperbola
#    5  Endpoints        ESS, HHT-QoL, transfusion, IV iron
#    6  Cardiac          Fick decomposition and the alpha sweep
#    7  Pulmonary shunt  SpO2 versus embolic hazard
#    8  Trial replica    PATH-HHT two-arm reproduction
#    9  Hb frontier      Hb_ss surface over (B, A_net) with B_crit
#   10  Scenarios        side-by-side comparison table
#   11  Sensitivity      tornado on the 2-year endpoints
#   12  Model notes      assumptions, calibration targets, known limits
#
#  Run:  Rscript -e 'shiny::runApp("hht_shiny_app.R", port = 8080)'
# =====================================================================

suppressMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
})

# ---------------------------------------------------------------------
# load the model and its helpers WITHOUT running the whole analysis
# ---------------------------------------------------------------------
.here <- function(f) {
  cand <- c(f, file.path(dirname(sys.frame(1)$ofile %||% "."), f))
  for (p in cand) if (file.exists(p)) return(p)
  f
}
`%||%` <- function(a, b) if (is.null(a)) b else a

MODEL_FILE <- if (file.exists("hht_mrgsolve_model.R")) "hht_mrgsolve_model.R" else
  file.path("hereditary-hemorrhagic-telangiectasia", "hht_mrgsolve_model.R")

local({
  src <- readLines(MODEL_FILE, warn = FALSE)
  cut <- grep("^if \\(identical\\(environment\\(\\), globalenv\\(\\)\\)", src)[1]
  if (is.na(cut)) cut <- length(src) + 1
  eval(parse(text = paste(src[seq_len(cut - 1)], collapse = "\n")),
       envir = globalenv())
})

PAL <- c(nose = "#c62828", gi = "#ad1457", iron = "#f9a825",
         hb = "#d84315", card = "#1565c0", lung = "#00838f",
         drug = "#0277bd", ess = "#5d4037", grey = "#78909c",
         green = "#2e7d32", purple = "#6a1b9a")

lineplot <- function(x, ys, cols, labs, ylab, main = "", ylim = NULL,
                     lty = NULL, vlines = NULL, hlines = NULL,
                     xlab = "days") {
  if (is.null(ylim)) {
    rng <- range(unlist(ys), na.rm = TRUE)
    if (!all(is.finite(rng))) rng <- c(0, 1)
    if (diff(rng) == 0) rng <- rng + c(-1, 1) * max(0.5, abs(rng[1]) * 0.1)
    ylim <- rng + c(-0.04, 0.10) * diff(rng)
  }
  if (is.null(lty)) lty <- rep(1, length(ys))
  op <- par(mar = c(4.2, 4.4, 2.6, 1.0), cex.axis = 0.9, cex.lab = 0.95)
  on.exit(par(op), add = TRUE)
  plot(x, ys[[1]], type = "n", ylim = ylim, xlab = xlab, ylab = ylab,
       main = main, las = 1, bty = "l")
  grid(col = "#e8ebee", lty = 1)
  if (!is.null(hlines)) abline(h = hlines, col = "#b0bec5", lty = 3)
  if (!is.null(vlines)) abline(v = vlines, col = "#cfd8dc", lty = 2)
  for (i in seq_along(ys))
    lines(x, ys[[i]], col = cols[i], lwd = 2.2, lty = lty[i])
  legend("topright", legend = labs, col = cols, lwd = 2.2, lty = lty,
         bty = "n", cex = 0.85)
}

kv <- function(...) {
  d <- list(...)
  tags$table(class = "table table-sm",
    lapply(names(d), function(n)
      tags$tr(tags$td(tags$b(n)), tags$td(style = "text-align:right", d[[n]]))))
}

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, 'Helvetica Neue', Arial, sans-serif; }
    .well { background:#fbfcfd; border-color:#e3e8ec; }
    h4 { margin-top: 4px; color:#263238; }
    .note { font-size: 12px; color:#546e7a; line-height:1.5; }
    .keynum { font-size: 13px; }
    table.table-sm td { padding: 2px 8px; font-size: 12.5px;
                        border-top: 1px solid #eceff1; }
  "))),
  titlePanel("Hereditary Haemorrhagic Telangiectasia — QSP explorer"),
  p(class = "note", HTML(paste0(
    "Central identity: <b>Hb<sub>ss</sub> = A<sub>net</sub> / (0.0347 &middot; B)</b>. ",
    "Haemoglobin is a ratio of net iron supply to blood-loss rate, so it cannot ",
    "distinguish an antiangiogenic drug from an iron infusion. Calibrated to ",
    "PATH-HHT (pomalidomide, n=144), Dupuis-Girod 2012 (bevacizumab, n=25), ",
    "Azzopardi 2015 exposure-response and NOSE (topical, n=121)."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("geno", "Genotype",
                  c("ACVRL1 / HHT2" = 2, "ENG / HHT1" = 1,
                    "SMAD4 / JP-HHT" = 3), selected = 2),
      sliderInput("age", "Age (yr)", 18, 85, 59, 1),
      sliderInput("wt", "Weight (kg)", 45, 130, 75, 1),
      sliderInput("qbleed", "Bleeding flow per episode (mL/min)",
                  0.3, 6, 1.84, 0.02),
      sliderInput("bgi", "GI blood loss (mL/day)", 0, 60, 9.9, 0.1),
      sliderInput("lam0", "Episode frequency (per week)", 0, 21, 7, 0.5),
      sliderInput("pv0", "Pulmonary AVM burden (relative)", 0, 3, 1, 0.05),
      sliderInput("gs0", "Hepatic shunt burden (relative)", 0, 25, 1, 0.1),

      h4("Iron policy"),
      radioButtons("ivpol", NULL,
                   c("Demand-driven (clinical practice)" = "demand",
                     "Protocol-fixed rate (trial)" = "fixed"),
                   selected = "demand"),
      conditionalPanel("input.ivpol == 'demand'",
        sliderInput("ivtrig", "Hb trigger for IV iron (g/dL)", 8, 13, 11, 0.1),
        sliderInput("ivmax", "Max IV iron (mg/day)", 0, 60, 26, 1)),
      conditionalPanel("input.ivpol == 'fixed'",
        sliderInput("ivfix", "Fixed IV iron (mg/day)", 0, 60, 12, 0.5)),
      checkboxInput("oralfe", "Oral iron 200 mg/day", FALSE),
      checkboxInput("txallow", "Allow RBC transfusion", TRUE),

      h4("Therapy"),
      checkboxInput("pom", "Pomalidomide 4 mg/day", FALSE),
      checkboxInput("bev", "Bevacizumab 5 mg/kg induction x6 q14d", FALSE),
      conditionalPanel("input.bev == true",
        selectInput("bevmaint", "Maintenance interval",
                    c("none" = 0, "monthly" = 30, "2-monthly" = 61,
                      "3-monthly" = 91), selected = 30)),
      checkboxInput("paz", "Pazopanib", FALSE),
      conditionalPanel("input.paz == true",
        sliderInput("pazdose", "Pazopanib dose (mg/day)", 50, 800, 50, 50)),
      checkboxInput("txa", "Tranexamic acid 1 g tid", FALSE),
      checkboxInput("close", "Young nasal closure", FALSE),
      sliderInput("embol", "Pulmonary AVM occlusion achieved", 0, 1, 0, 0.05),
      sliderInput("giendo", "GI ablation (fractional lesion loss/day)",
                  0, 0.1, 0, 0.005),

      h4("Physiology / score assumptions"),
      sliderInput("alpha", "alpha — anaemia-compensation exponent",
                  0, 1, 1, 0.05),
      sliderInput("wfe", "ESS anaemia+transfusion weight share",
                  0, 0.35, 0.18, 0.01),
      sliderInput("fresl", "Hepatic shunt VEGF-independent share",
                  0.1, 0.95, 0.30, 0.05),
      sliderInput("horizon", "Horizon (days)", 180, 1825, 730, 30),
      checkboxInput("blind", "Blinded study drug (placebo response on)",
                    FALSE),
      hr(),
      p(class = "note", "Educational / research model only. Not for clinical use.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 Patient",
          fluidRow(column(5, h4("Baseline equilibrium"), uiOutput("baseKV")),
                   column(7, h4("Where this patient sits on the hyperbola"),
                          plotOutput("hypPlot", height = 320))),
          hr(), p(class = "note", HTML(paste0(
            "The baseline is a true equilibrium of all 40 ODEs, reached by a ",
            "burn-in before any therapy starts. <b>Hb predicted</b> is the ",
            "closed form; it agrees with the solved ODEs to well under 1% ",
            "whenever erythropoiesis is iron-limited, because the identity is ",
            "a property of mass balance rather than a fit."))),
          h4("Blood-loss budget"), plotOutput("budgetPlot", height = 240)),

        tabPanel("2 Drug PK",
          fluidRow(column(6, h4("Concentrations"),
                          plotOutput("pkPlot", height = 300)),
                   column(6, h4("Bevacizumab target occupancy"),
                          plotOutput("occPlot", height = 300))),
          hr(), p(class = "note", HTML(paste0(
            "The occupancy panel is the point of this tab. Bevacizumab sits in ",
            "a 10<sup>4</sup>-10<sup>5</sup> molar excess over plasma VEGF, so ",
            "free VEGF stays suppressed for many half-lives after the ",
            "concentration has become negligible. The exposure metric for a ",
            "stoichiometric neutralising antibody is concentration relative to ",
            "<i>target molarity</i>, not concentration and not AUC."))),
          h4("Molar excess over target (log scale)"),
          plotOutput("excessPlot", height = 240)),

        tabPanel("3 Bleeding & lesion",
          fluidRow(column(6, h4("Blood loss"), plotOutput("bldPlot", height = 300)),
                   column(6, h4("Lesion burden and mural coverage"),
                          plotOutput("lesPlot", height = 300))),
          hr(),
          fluidRow(column(6, h4("Episode frequency and duration"),
                          plotOutput("epiPlot", height = 280)),
                   column(6, h4("f_GI — the share ESS cannot see"),
                          plotOutput("fgiPlot", height = 280))),
          p(class = "note", HTML(paste0(
            "Pomalidomide is implemented as a <i>vessel-maturation</i> agent, ",
            "not an anti-VEGF one, because the PATH-HHT nasal biopsies showed ",
            "no inhibition of endothelial proliferation but increased mural-cell ",
            "coverage. Coverage therefore acts on all three bleeding dimensions ",
            "at once: rupture frequency, flow and duration.")))),

        tabPanel("4 Iron & Hb",
          fluidRow(column(6, h4("Iron fluxes"), plotOutput("fePlot", height = 300)),
                   column(6, h4("Haemoglobin and stores"),
                          plotOutput("hbPlot", height = 300))),
          hr(),
          fluidRow(column(6, h4("Hepcidin / erythroferrone"),
                          plotOutput("hepPlot", height = 260)),
                   column(6, h4("ODE Hb versus the closed form"),
                          plotOutput("idPlot", height = 260))),
          p(class = "note", HTML(paste0(
            "Iron balance: <b>dFe/dt = A_net &minus; 0.0347&middot;Hb&middot;B</b>. ",
            "Because the iron cost of bleeding is proportional to Hb, the drain ",
            "shrinks as Hb falls and the patient settles at a low haemoglobin ",
            "instead of exsanguinating: anaemia in HHT is an equilibrium, not a ",
            "decompensation.")))),

        tabPanel("5 Endpoints",
          fluidRow(column(6, h4("ESS and HHT-QoL"),
                          plotOutput("essPlot", height = 300)),
                   column(6, h4("Transfusion and IV iron burden"),
                          plotOutput("burdenPlot", height = 300))),
          hr(), h4("ESS item decomposition at the horizon"),
          plotOutput("essItems", height = 260),
          p(class = "note", HTML(paste0(
            "ESS contains anaemia and transfusion items (Hoag 2010), so part of ",
            "any ESS response is reachable by iron alone. Move the ",
            "<i>ESS anaemia+transfusion weight share</i> slider and watch how ",
            "much of the score an iron change can buy. The item weights were ",
            "never published in closed form; this is a transparent surrogate, ",
            "not the instrument.")))),

        tabPanel("6 Cardiac",
          fluidRow(column(6, h4("Cardiac index, decomposed"),
                          plotOutput("ciPlot", height = 300)),
                   column(6, h4("alpha sweep — the identifiability problem"),
                          plotOutput("alphaPlot", height = 300))),
          hr(), h4("NT-proBNP and RV reserve"),
          plotOutput("bnpPlot", height = 240),
          p(class = "note", HTML(paste0(
            "By Fick, C(a&minus;v)O<sub>2</sub> is proportional to Hb, so ",
            "CI scales as Hb<sup>&minus;alpha</sup>. At alpha = 1 correcting ",
            "anaemia lowers cardiac index with no change in shunt anatomy at ",
            "all; at alpha = 0 it does nothing. alpha is not identifiable from ",
            "any published HHT dataset, which makes the shunt-attributable ",
            "share of a cardiac-index response unidentifiable too.")))),

        tabPanel("7 Pulmonary shunt",
          fluidRow(column(6, h4("SpO2 and shunt fraction"),
                          plotOutput("spo2Plot", height = 300)),
                   column(6, h4("SpO2 versus embolic hazard across burden"),
                          plotOutput("screenPlot", height = 300))),
          hr(), h4("Shunt fraction needed to reach a given SpO2"),
          tableOutput("sstarTab"),
          p(class = "note", HTML(paste0(
            "Embolic risk (paradoxical stroke, cerebral abscess) is set by the ",
            "anatomy. SpO<sub>2</sub> is set by the anatomy <i>and</i> ",
            "haemoglobin <i>and</i> cardiac output. The two do not co-vary, ",
            "which is the quantitative reason the guideline screens with ",
            "contrast echocardiography rather than a pulse oximeter.")))),

        tabPanel("8 Trial replica",
          h4("PATH-HHT reproduced: pomalidomide 4 mg vs placebo, 24 weeks"),
          p(class = "note", HTML(paste0(
            "Both arms run the <b>same</b> demand-driven iron policy, so the ",
            "fall in iron use is an output of the model rather than an ",
            "assumption. Switch the sidebar iron policy to <i>protocol-fixed</i> ",
            "to run the counterfactual trial that was never done."))),
          plotOutput("trialPlot", height = 340),
          hr(), h4("Model versus published"), tableOutput("trialTab"),
          h4("Iron and transfusion burden by arm"), tableOutput("trialIron")),

        tabPanel("9 Hb frontier",
          h4("Hb_ss over blood loss and net iron supply"),
          plotOutput("frontPlot", height = 420),
          hr(), h4("B_crit — blood loss a regimen can just compensate"),
          tableOutput("bcritTab"),
          p(class = "note", HTML(paste0(
            "Contours are the closed form. The dashed line is the frontier ",
            "A_net = 0.0347&middot;Hb_set&middot;B: above it there is no ",
            "equilibrium in iron stores at all and the patient becomes ",
            "iron-loaded rather than polycythaemic.")))),

        tabPanel("10 Scenarios",
          h4("Standard arms at the chosen horizon"),
          p(class = "note", paste("All arms start from the same equilibrium",
            "patient defined in the sidebar, then diverge.")),
          tableOutput("scenTab"),
          hr(), h4("Ranked by haemoglobin, and by ESS"),
          fluidRow(column(6, tableOutput("rankHb")),
                   column(6, tableOutput("rankESS"))),
          p(class = "note", paste("The two rankings disagree. Arms that buy",
            "haemoglobin with iron rank well on Hb and poorly on ESS; arms",
            "that remove bleeding rank well on both but need time."))),

        tabPanel("11 Sensitivity",
          h4("Local sensitivity, +/-25%, at the horizon"),
          plotOutput("tornado", height = 460),
          tableOutput("sensTab")),

        tabPanel("12 Model notes",
          h4("What is calibrated, and to what"),
          uiOutput("notesUI"))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  ironpar <- reactive({
    if (identical(input$ivpol, "demand"))
      list(IVFE_ON = 1, IVFERATE = 0, IVFETRIG = input$ivtrig,
           IVFEMAX = input$ivmax)
    else
      list(IVFE_ON = 0, IVFERATE = input$ivfix)
  })

  parlist <- reactive({
    c(list(GENO = as.numeric(input$geno), AGE0 = input$age, WT = input$wt,
           BSA = 0.0235 * input$wt^0.51456 * 170^0.42246,
           QBLEED = input$qbleed, BGI0 = input$bgi, LAMEP0 = input$lam0,
           PV0 = input$pv0, GS0 = input$gs0, CLOSEF = as.numeric(input$close),
           EMBOLF = input$embol, GIENDOF = input$giendo,
           ALPHA = input$alpha, W_FE_SHARE = input$wfe, FRES_L = input$fresl,
           KII = input$qbleed, TXALLOW = as.numeric(input$txallow),
           BLIND = as.numeric(input$blind), ONPLAC = as.numeric(input$blind)),
      ironpar())
  })

  evlist <- reactive({
    H <- input$horizon; e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    if (input$pom)    e <- add(e, ev_pom(H))
    if (input$paz)    e <- add(e, ev_paz(H, input$pazdose))
    if (input$txa)    e <- add(e, ev_txa(H))
    if (input$oralfe) e <- add(e, ev_oral_fe(H, 200))
    if (input$bev) {
      e <- add(e, ev_bev(6, 14, input$wt))
      ii <- as.numeric(input$bevmaint)
      if (ii > 0) {
        m <- ev_bev_maint(ii, H, input$wt)
        if (!is.null(m)) e <- add(e, m)
      }
    }
    e
  })

  sim <- reactive({
    do.call(sim_eq, c(list(end = input$horizon, delta = 2,
                           events = evlist()), parlist()))
  })
  base <- reactive({ d <- sim(); d[1, ] })
  last <- reactive({ d <- sim(); d[nrow(d), ] })

  # ---- tab 1 --------------------------------------------------------
  output$baseKV <- renderUI({
    b <- base(); l <- last()
    kv("Haemoglobin (g/dL)"       = sprintf("%.2f  ->  %.2f", b$HB, l$HB),
       "Hb predicted, closed form"= sprintf("%.2f", b$HB_PRED),
       "Blood loss B (mL/day)"    = sprintf("%.1f  ->  %.1f", b$B_TOT, l$B_TOT),
       "  nasal"                  = sprintf("%.1f  ->  %.1f", b$B_NOSE, l$B_NOSE),
       "  gastrointestinal"       = sprintf("%.1f  ->  %.1f", b$B_GI, l$B_GI),
       "f_GI"                     = sprintf("%.2f", b$F_GI),
       "A_net (mg Fe/day)"        = sprintf("%.2f  ->  %.2f", b$ANET, l$ANET),
       "IV iron (mg per 4 wk)"    = sprintf("%.0f  ->  %.0f",
                                            b$IVFERT * 28, l$IVFERT * 28),
       "B_crit at Hb 11.5"        = sprintf("%.1f mL/day", b$B_CRIT),
       "Ferritin (ng/mL)"         = sprintf("%.0f  ->  %.0f", b$FERRITIN, l$FERRITIN),
       "ESS (0-10)"               = sprintf("%.2f  ->  %.2f", b$ESSC, l$ESSC),
       "HHT-QoL (0-16)"           = sprintf("%.2f  ->  %.2f", b$QOLC, l$QOLC),
       "Cardiac index"            = sprintf("%.2f  ->  %.2f", b$CI, l$CI),
       "SpO2 (%)"                 = sprintf("%.1f  ->  %.1f", b$SPO2, l$SPO2),
       "Shunt fraction"           = sprintf("%.3f", b$SFRACT),
       "RBC units / yr"           = sprintf("%.2f", l$TXPERYR))
  })

  output$hypPlot <- renderPlot({
    b <- base()
    B <- seq(2, 140, length.out = 300)
    op <- par(mar = c(4.2, 4.4, 2, 1)); on.exit(par(op))
    plot(B, b$ANET / (0.0347 * B), type = "l", lwd = 2.4, col = PAL["hb"],
         ylim = c(0, 16), las = 1, bty = "l",
         xlab = "blood loss B (mL/day)", ylab = "equilibrium Hb (g/dL)",
         main = sprintf("Hb_ss = %.1f / (0.0347 B)", b$ANET))
    grid(col = "#e8ebee", lty = 1)
    abline(h = c(11.5, 13), col = "#b0bec5", lty = 3)
    points(b$B_TOT, b$HB, pch = 19, col = PAL["card"], cex = 1.8)
    points(last()$B_TOT, last()$HB, pch = 17, col = PAL["green"], cex = 1.8)
    legend("topright", c("baseline", "at horizon", "Hb 11.5 / 13"),
           pch = c(19, 17, NA), lty = c(NA, NA, 3),
           col = c(PAL["card"], PAL["green"], "#b0bec5"), bty = "n", cex = 0.85)
  })

  output$budgetPlot <- renderPlot({
    d <- sim()
    op <- par(mar = c(4.2, 4.4, 2, 1)); on.exit(par(op))
    ymax <- max(d$B_TOT, na.rm = TRUE) * 1.15 + 1
    plot(d$time, d$B_TOT, type = "n", ylim = c(0, ymax), las = 1, bty = "l",
         xlab = "days", ylab = "blood loss (mL/day)",
         main = "nasal + gastrointestinal")
    grid(col = "#e8ebee", lty = 1)
    polygon(c(d$time, rev(d$time)), c(d$B_GI, rev(rep(0, nrow(d)))),
            col = "#f8bbd0", border = NA)
    polygon(c(d$time, rev(d$time)), c(d$B_TOT, rev(d$B_GI)),
            col = "#ffcdd2", border = NA)
    lines(d$time, d$B_TOT, lwd = 2.2, col = PAL["nose"])
    lines(d$time, d$B_GI, lwd = 2.0, col = PAL["gi"])
    legend("topright", c("total", "GI only"), col = c(PAL["nose"], PAL["gi"]),
           lwd = 2.2, bty = "n", cex = 0.85)
  })

  # ---- tab 2 --------------------------------------------------------
  output$pkPlot <- renderPlot({
    d <- sim()
    ys <- list(); cols <- c(); labs <- c()
    if (max(d$CBEVT) > 0) { ys <- c(ys, list(d$CBEVT)); cols <- c(cols, PAL["drug"]); labs <- c(labs, "bevacizumab mg/L") }
    if (max(d$CPOMT) > 0) { ys <- c(ys, list(d$CPOMT * 1000)); cols <- c(cols, PAL["green"]); labs <- c(labs, "pomalidomide ng/mL") }
    if (max(d$CPAZT) > 0) { ys <- c(ys, list(d$CPAZT)); cols <- c(cols, PAL["purple"]); labs <- c(labs, "pazopanib mg/L") }
    if (max(d$CTXAT) > 0) { ys <- c(ys, list(d$CTXAT)); cols <- c(cols, PAL["ess"]); labs <- c(labs, "TXA mg/L") }
    if (!length(ys)) { plot.new(); title("no drug selected"); return() }
    lineplot(d$time, ys, cols, labs, "concentration")
  })

  output$occPlot <- renderPlot({
    d <- sim()
    if (max(d$CBEVT) <= 0) { plot.new(); title("bevacizumab not selected"); return() }
    lt <- d$VEGFT * 2.618e-5; bt <- d$CBEVT * 6.711; KD <- mod$KD_BEV
    bb <- lt - bt - KD
    lf <- 0.5 * (bb + sqrt(pmax(0, bb * bb + 4 * KD * lt)))
    occ <- 100 * (1 - lf / max(lf[1], 1e-12))
    lineplot(d$time, list(occ, 100 * d$CBEVT / max(d$CBEVT)),
             c(PAL["green"], PAL["drug"]),
             c("VEGF neutralised (%)", "drug, % of peak"),
             "%", "occupancy outlasts concentration", ylim = c(0, 108),
             hlines = c(50, 100))
  })

  output$excessPlot <- renderPlot({
    d <- sim()
    if (max(d$CBEVT) <= 0) { plot.new(); title("bevacizumab not selected"); return() }
    exc <- (d$CBEVT * 6.711) / pmax(d$VEGFT * 2.618e-5, 1e-12)
    op <- par(mar = c(4.2, 4.6, 2, 1)); on.exit(par(op))
    plot(d$time, pmax(exc, 1e-2), type = "l", log = "y", lwd = 2.2,
         col = PAL["drug"], las = 1, bty = "l", xlab = "days",
         ylab = "molar excess (antibody : VEGF)",
         main = "engagement persists while this exceeds 1")
    grid(col = "#e8ebee", lty = 1); abline(h = 1, col = PAL["nose"], lty = 2)
  })

  # ---- tab 3 --------------------------------------------------------
  output$bldPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$B_TOT, d$B_NOSE, d$B_GI),
             c(PAL["hb"], PAL["nose"], PAL["gi"]),
             c("total", "nasal", "GI"), "mL/day", "blood loss")
  })
  output$lesPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$TNOSE, d$TGIB, d$GSHUNT / max(d$GSHUNT[1], 1e-9),
                          d$MURALC / mod$MURSS),
             c(PAL["nose"], PAL["gi"], PAL["green"], PAL["card"]),
             c("nasal lesion", "GI lesion", "hepatic shunt (rel)",
               "mural coverage (rel)"),
             "relative to baseline", "lesion burden and vessel maturation",
             hlines = 1)
  })
  output$epiPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$LAMEP, d$EPDURC),
             c(PAL["nose"], PAL["ess"]),
             c("episodes / week", "duration (min)"), "value")
  })
  output$fgiPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$F_GI), PAL["gi"], "f_GI", "fraction",
             "GI share of total blood loss", ylim = c(0, 1), hlines = 0.5)
  })

  # ---- tab 4 --------------------------------------------------------
  output$fePlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$ANET, d$FELOSS, d$IVFERT),
             c(PAL["green"], PAL["nose"], PAL["iron"]),
             c("A_net supply", "0.0347*Hb*B loss", "IV iron"),
             "mg Fe/day", "iron in versus iron out")
  })
  output$hbPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$HB, d$FERRITIN / 10),
             c(PAL["hb"], PAL["iron"]),
             c("Hb (g/dL)", "ferritin / 10 (ng/mL)"), "value",
             hlines = c(11.5, 13))
  })
  output$hepPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$HEPCC, d$ERFEC, d$RETICC),
             c(PAL["ess"], PAL["green"], PAL["hb"]),
             c("hepcidin", "erythroferrone", "marrow drive"),
             "relative", hlines = 1)
  })
  output$idPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$HB, d$HB_PRED), c(PAL["hb"], PAL["grey"]),
             c("Hb, solved ODEs", "Hb, closed form"), "g/dL",
             "the identity, checked continuously", lty = c(1, 2))
  })

  # ---- tab 5 --------------------------------------------------------
  output$essPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$ESSC, d$QOLC), c(PAL["ess"], PAL["purple"]),
             c("ESS (0-10)", "HHT-QoL (0-16)"), "score",
             hlines = c(4, 7))
  })
  output$burdenPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$CUMTX, d$CUMIVFE / 1000),
             c(PAL["nose"], PAL["iron"]),
             c("cumulative RBC units", "cumulative IV iron (g)"), "cumulative")
  })
  output$essItems <- renderPlot({
    l <- last()
    lam <- l$LAMEP; dur <- l$EPDURC; q <- l$QEFFT; hb <- l$HB
    iF <- lam / (lam + mod$KIF); iD <- dur / (dur + mod$KID)
    iI <- q / (q + mod$KII)
    iA <- 1 / (1 + exp((hb - mod$HBANEM) / mod$SHB))
    iT <- l$TXPERYR / 365.25 / (l$TXPERYR / 365.25 + mod$KITX)
    wn <- 1 - input$wfe
    ws <- mod$W_FREQ + mod$W_DUR + mod$W_INT + mod$W_ATT
    contrib <- c(
      frequency = 10 * wn * mod$W_FREQ * iF / ws,
      duration  = 10 * wn * mod$W_DUR * iD / ws,
      intensity = 10 * wn * mod$W_INT * iI / ws,
      attention = 10 * wn * mod$W_ATT * (0.5 * iF + 0.5 * iD) / ws,
      anaemia   = 10 * input$wfe * 0.6 * iA,
      transfusion = 10 * input$wfe * 0.4 * iT)
    op <- par(mar = c(4.2, 8, 2.4, 1)); on.exit(par(op))
    cols <- c(rep(PAL["nose"], 4), PAL["iron"], PAL["iron"])
    barplot(rev(contrib), horiz = TRUE, las = 1, col = rev(cols),
            border = NA, xlab = "ESS points contributed",
            main = sprintf("ESS %.2f — iron-sensitive items in amber (%.0f%% of weight)",
                           sum(contrib), 100 * input$wfe))
    grid(col = "#e8ebee", lty = 1)
  })

  # ---- tab 6 --------------------------------------------------------
  output$ciPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$CI, d$CI_TIS, d$CI_SH),
             c(PAL["card"], PAL["green"], PAL["nose"]),
             c("cardiac index", "Fick / tissue term", "shunt term"),
             "L/min/m2", "CI = Fick(Hb) + shunt", hlines = 4)
  })
  output$alphaPlot <- renderPlot({
    ps <- parlist(); ev <- evlist(); H <- input$horizon
    al <- c(0, 0.25, 0.5, 0.75, 1)
    res <- sapply(al, function(a) {
      p <- ps; p$ALPHA <- a
      d <- do.call(sim_eq, c(list(end = H, delta = 30, events = ev), p))
      c(d$CI[1], d$CI[nrow(d)])
    })
    op <- par(mar = c(4.2, 4.6, 2.4, 1)); on.exit(par(op))
    plot(al, res[2, ] - res[1, ], type = "b", pch = 19, lwd = 2.2,
         col = PAL["card"], las = 1, bty = "l",
         xlab = "alpha (anaemia-compensation exponent)",
         ylab = "change in cardiac index",
         main = "the same drug, a different answer")
    grid(col = "#e8ebee", lty = 1); abline(h = 0, col = "#b0bec5", lty = 3)
  })
  output$bnpPlot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$NTBNPC, 1000 * d$RVF),
             c(PAL["purple"], PAL["card"]),
             c("NT-proBNP (pg/mL)", "RV reserve x1000"), "value")
  })

  # ---- tab 7 --------------------------------------------------------
  output$spo2Plot <- renderPlot({
    d <- sim()
    lineplot(d$time, list(d$SPO2, 100 * d$SFRACT),
             c(PAL["lung"], PAL["nose"]),
             c("SpO2 (%)", "shunt fraction x100"), "value", hlines = c(90, 95))
  })
  output$screenPlot <- renderPlot({
    ps <- parlist()
    pv <- seq(0, 3, by = 0.1)
    out <- t(sapply(pv, function(p) {
      q <- ps; q$PV0 <- p; q$TXALLOW <- 0
      s <- do.call(eqstate, q)
      c(s$SPO2, s$SFRACT, mod$HAZPAV * s$PAVMB * (1 - ps$EMBOLF))
    }))
    op <- par(mar = c(4.2, 4.6, 2.4, 4.6)); on.exit(par(op))
    plot(pv, out[, 1], type = "l", lwd = 2.4, col = PAL["lung"], las = 1,
         bty = "l", ylim = c(60, 100), xlab = "pulmonary AVM burden",
         ylab = "SpO2 (%)", main = "oximetry is blind to the low-shunt zone")
    grid(col = "#e8ebee", lty = 1)
    abline(h = 95, col = "#b0bec5", lty = 3)
    ok <- out[, 1] >= 95
    if (any(ok)) rect(0, 60, max(pv[ok]), 100, col = "#ffebee", border = NA)
    lines(pv, out[, 1], lwd = 2.4, col = PAL["lung"])
    par(new = TRUE)
    plot(pv, out[, 3], type = "l", lwd = 2.2, lty = 2, col = PAL["nose"],
         axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1); mtext("embolic hazard / yr", side = 4, line = 3, cex = 0.9)
    legend("bottomleft", c("SpO2", "embolic hazard", "SpO2 >= 95% zone"),
           col = c(PAL["lung"], PAL["nose"], "#ffcdd2"),
           lwd = c(2.4, 2.2, 8), lty = c(1, 2, 1), bty = "n", cex = 0.85)
  })
  output$sstarTab <- renderTable({
    sc <- mod$SCAP
    f <- function(ds, tgt) { g <- sc - tgt; g / (ds + g) }
    data.frame(
      `dSat` = c(0.20, 0.25, 0.30, 0.35),
      `s for SpO2 97%` = sapply(c(0.20, 0.25, 0.30, 0.35), f, tgt = 0.97),
      `s for SpO2 95%` = sapply(c(0.20, 0.25, 0.30, 0.35), f, tgt = 0.95),
      `s for SpO2 90%` = sapply(c(0.20, 0.25, 0.30, 0.35), f, tgt = 0.90),
      check.names = FALSE)
  }, digits = 3)

  # ---- tab 8 --------------------------------------------------------
  trial <- reactive({
    D <- 24 * 7; ps <- parlist()
    ps$BLIND <- 1; ps$ONPLAC <- 1
    pl <- do.call(sim_eq, c(list(end = D + 28, delta = 2, events = NULL), ps))
    po <- do.call(sim_eq, c(list(end = D + 28, delta = 2,
                                 events = ev_pom(D)), ps))
    list(pl = pl, po = po, D = D)
  })
  output$trialPlot <- renderPlot({
    t <- trial(); op <- par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1))
    on.exit(par(op))
    for (v in c("ESSC", "HB", "IVFERT")) {
      y1 <- t$pl[[v]]; y2 <- t$po[[v]]
      if (v == "IVFERT") { y1 <- y1 * 28; y2 <- y2 * 28 }
      ttl <- c(ESSC = "ESS", HB = "haemoglobin (g/dL)",
               IVFERT = "IV iron (mg per 4 wk)")[v]
      plot(t$pl$time, y1, type = "l", lwd = 2.4, col = PAL["grey"],
           ylim = range(c(y1, y2)), las = 1, bty = "l", xlab = "days",
           ylab = "", main = ttl)
      grid(col = "#e8ebee", lty = 1)
      lines(t$po$time, y2, lwd = 2.4, col = PAL["green"])
      abline(v = t$D, lty = 2, col = "#cfd8dc")
      legend("topright", c("placebo", "pomalidomide"),
             col = c(PAL["grey"], PAL["green"]), lwd = 2.4, bty = "n",
             cex = 0.9)
    }
  })
  output$trialTab <- renderTable({
    t <- trial(); D <- t$D
    at <- function(d, tt) d[which.min(abs(d$time - tt)), ]
    dl <- function(d, v, tt) at(d, tt)[[v]] - at(d, 0)[[v]]
    data.frame(
      endpoint = c("ESS change, pomalidomide", "ESS change, placebo",
                   "ESS difference", "HHT-QoL difference",
                   "Hb difference, week 24",
                   "Hb difference, 4 weeks post-treatment"),
      published = c("-1.84", "-0.90", "-0.94", "-1.40",
                    "not significant", "+1.09"),
      model = c(dl(t$po, "ESSC", D), dl(t$pl, "ESSC", D),
                dl(t$po, "ESSC", D) - dl(t$pl, "ESSC", D),
                dl(t$po, "QOLC", D) - dl(t$pl, "QOLC", D),
                dl(t$po, "HB", D) - dl(t$pl, "HB", D),
                dl(t$po, "HB", D + 28) - dl(t$pl, "HB", D + 28)))
  }, digits = 3)
  output$trialIron <- renderTable({
    t <- trial(); D <- t$D
    at <- function(d, tt) d[which.min(abs(d$time - tt)), ]
    data.frame(
      quantity = c("IV iron at week 24 (mg per 4 wk)", "A_net (mg Fe/day)",
                   "B_total (mL/day)", "Hb (g/dL)",
                   "implied blood loss from iron use (mL/day)"),
      placebo = c(at(t$pl, D)$IVFERT * 28, at(t$pl, D)$ANET,
                  at(t$pl, D)$B_TOT, at(t$pl, D)$HB,
                  at(t$pl, D)$IVFERT / (0.0347 * at(t$pl, D)$HB)),
      pomalidomide = c(at(t$po, D)$IVFERT * 28, at(t$po, D)$ANET,
                       at(t$po, D)$B_TOT, at(t$po, D)$HB,
                       at(t$po, D)$IVFERT / (0.0347 * at(t$po, D)$HB)))
  }, digits = 2)

  # ---- tab 9 --------------------------------------------------------
  output$frontPlot <- renderPlot({
    B <- seq(2, 140, length.out = 120)
    A <- seq(0.5, 45, length.out = 120)
    Z <- outer(B, A, function(b, a) pmin(mod$HBSET, a / (0.0347 * b)))
    op <- par(mar = c(4.4, 4.6, 2.6, 1)); on.exit(par(op))
    filled.contour(B, A, Z, nlevels = 18,
      color.palette = colorRampPalette(c("#7f0000", "#ef9a9a", "#fff59d",
                                         "#a5d6a7", "#2e7d32")),
      plot.title = title(xlab = "blood loss B (mL/day)",
                         ylab = "A_net (mg Fe/day)",
                         main = "equilibrium haemoglobin (g/dL), capped at the set-point"),
      plot.axes = {
        axis(1); axis(2)
        contour(B, A, Z, levels = c(8, 10, 11.5, 13), add = TRUE,
                col = "#263238", lwd = 1.1, labcex = 0.8)
        lines(B, 0.0347 * mod$HBSET * B, col = "#000000", lty = 2, lwd = 2)
        points(base()$B_TOT, base()$ANET, pch = 19, cex = 2, col = "#0d47a1")
      })
  })
  output$bcritTab <- renderTable({
    ps <- parlist(); ps$IVFE_ON <- 0; ps$TXALLOW <- 0
    do.call(rbind, lapply(c(0, 5, 12, 25, 50), function(iv) {
      q <- ps; q$IVFERATE <- iv
      s <- do.call(eqstate, q)
      data.frame(`IV iron (mg/day)` = iv, `A_net (mg Fe/day)` = s$ANET,
                 `B_crit at Hb 11.5 (mL/day)` = s$B_CRIT,
                 `equilibrium Hb (g/dL)` = s$HB, check.names = FALSE)
    }))
  }, digits = 2)

  # ---- tab 10 -------------------------------------------------------
  output$scenTab <- renderTable({ scen()$tab }, digits = 3)
  scen <- reactive({
    H <- input$horizon; ps <- parlist()
    arms <- list(
      "no therapy"                 = list(ev = NULL, p = list()),
      "oral iron 200 mg/day"       = list(ev = ev_oral_fe(H, 200), p = list()),
      "IV iron 30 mg/day (fixed)"  = list(ev = NULL,
                                          p = list(IVFE_ON = 0, IVFERATE = 30)),
      "tranexamic acid"            = list(ev = ev_txa(H), p = list()),
      "pomalidomide 4 mg"          = list(ev = ev_pom(H), p = list()),
      "pazopanib 50 mg"            = list(ev = ev_paz(H, 50), p = list()),
      "bevacizumab, monthly maint" = list(
        ev = c(ev_bev(6, 14, input$wt), ev_bev_maint(30, H, input$wt)),
        p = list()),
      "bevacizumab, 3-monthly"     = list(
        ev = c(ev_bev(6, 14, input$wt), ev_bev_maint(91, H, input$wt)),
        p = list()),
      "Young nasal closure"        = list(ev = NULL, p = list(CLOSEF = 1)))
    tab <- do.call(rbind, lapply(names(arms), function(nm) {
      a <- arms[[nm]]
      p <- modifyList(ps, a$p)
      d <- do.call(sim_eq, c(list(end = H, delta = 14, events = a$ev), p))
      s0 <- d[1, ]; s1 <- d[nrow(d), ]
      data.frame(scenario = nm, ESS = s1$ESSC, dESS = s1$ESSC - s0$ESSC,
                 Hb = s1$HB, dHb = s1$HB - s0$HB, B = s1$B_TOT,
                 CI = s1$CI, `RBC units` = s1$CUMTX,
                 `IV iron (g)` = s1$CUMIVFE / 1000, QoL = s1$QOLC,
                 check.names = FALSE)
    }))
    list(tab = tab)
  })
  output$rankHb <- renderTable({
    t <- scen()$tab; t <- t[order(-t$Hb), c("scenario", "Hb", "ESS")]
    t }, digits = 3)
  output$rankESS <- renderTable({
    t <- scen()$tab; t <- t[order(t$ESS), c("scenario", "ESS", "Hb")]
    t }, digits = 3)

  # ---- tab 11 -------------------------------------------------------
  sens <- reactive({
    H <- input$horizon; ps <- parlist(); ev <- evlist()
    pars <- c("KOUT_N", "KOUT_L", "AMAXO", "KHEP", "BGI0", "QBLEED",
              "MURSS", "GFRAG", "ALPHA", "KSHUNT", "FRES", "FRES_L",
              "EMAX_PM", "W_FE_SHARE", "IVFETRIG", "KSUP", "IL6G")
    b <- do.call(sim_eq, c(list(end = H, delta = 60, events = ev), ps))
    bl <- b[nrow(b), ]
    do.call(rbind, lapply(pars, function(pn) {
      v <- if (!is.null(ps[[pn]])) ps[[pn]] else mod[[pn]]
      if (is.null(v) || !is.finite(v) || v == 0)
        return(data.frame(parameter = pn, Hb_pct = 0, ESS_abs = 0, CI_abs = 0))
      lo <- modifyList(ps, setNames(list(v * 0.75), pn))
      hi <- modifyList(ps, setNames(list(v * 1.25), pn))
      dl <- do.call(sim_eq, c(list(end = H, delta = 60, events = ev), lo))
      dh <- do.call(sim_eq, c(list(end = H, delta = 60, events = ev), hi))
      l <- dl[nrow(dl), ]; h <- dh[nrow(dh), ]
      data.frame(parameter = pn,
                 Hb_pct = 100 * (h$HB - l$HB) / bl$HB,
                 ESS_abs = h$ESSC - l$ESSC, CI_abs = h$CI - l$CI)
    }))
  })
  output$tornado <- renderPlot({
    s <- sens(); s <- s[order(abs(s$Hb_pct)), ]
    op <- par(mar = c(4.4, 8.5, 2.6, 1)); on.exit(par(op))
    barplot(s$Hb_pct, horiz = TRUE, names.arg = s$parameter, las = 1,
            col = ifelse(s$Hb_pct > 0, PAL["green"], PAL["nose"]),
            border = NA, xlab = "% change in haemoglobin at the horizon",
            main = "+/-25% parameter moves")
    grid(col = "#e8ebee", lty = 1); abline(v = 0, col = "#455a64")
  })
  output$sensTab <- renderTable({
    s <- sens(); s[order(-abs(s$Hb_pct)), ] }, digits = 3)

  # ---- tab 12 -------------------------------------------------------
  output$notesUI <- renderUI({
    tagList(
      h4("Calibration targets"),
      tags$ul(
        tags$li(HTML("<b>PATH-HHT</b> (Al-Samkari 2024 NEJM 391:1015, PMID 39292928), n=144, randomised 2:1. Baseline ESS 5.0&plusmn;1.5, HHT-QoL 6.3&plusmn;3.1, age 58.8&plusmn;12.2; 69% anaemic; 84% IV iron and 19% transfusion in the prior 6 months. ESS &minus;1.84 vs &minus;0.90, difference &minus;0.94 (&minus;1.57 to &minus;0.31). HHT-QoL &minus;2.7 vs &minus;1.2. Median IV iron 0 vs 333 mg per 4 weeks. Haemoglobin not different at 24 weeks; +1.09 g/dL at 4 weeks post-treatment. Nasal biopsy: no inhibition of endothelial proliferation, increased mural-cell coverage.")),
        tags$li(HTML("<b>Dupuis-Girod 2012</b> (JAMA 307:948, PMID 22396517), single arm, n=25. Cardiac index 5.05 &rarr; 4.2 &rarr; 4.1; epistaxis 221 &rarr; 134 &rarr; 43 min/month; last dose day 70.")),
        tags$li(HTML("<b>Azzopardi 2015</b> (MAbs 7:630, PMID 25751241). Maintenance q3/q2/q1 month at 24 months: cardiac index &lt;4 in 41/45/50%, epistaxis &lt;20 min/month in 34/43/60%.")),
        tags$li(HTML("<b>NOSE</b> (Whitehead 2016 JAMA 316:943, PMID 27599329), n=121. Topical bevacizumab, estriol and tranexamic acid all no better than saline; every arm improved ESS significantly.")),
        tags$li(HTML("<b>ESS instrument</b> (Hoag 2010 Laryngoscope 120:838, PMID 20087969); MCID 0.71 (Yin 2016, PMID 26393959)."))),
      h4("Known limits — please read before quoting a number"),
      tags$ol(
        tags$li("The ESS item weights were never published in closed form. The six-item surrogate here is transparent and monotone in the right variables, but absolute ESS values are not on the same scale as trial ESS values."),
        tags$li("B_nose, B_gi and f_GI are inferred through the iron balance, never observed. No HHT trial measures blood loss in mL/day."),
        tags$li(HTML("A hypothesis was <b>rejected</b>: that the 6-month bevacizumab epistaxis nadir was pharmacologically impossible. Stoichiometric molar excess keeps the target engaged for many half-lives after the concentration becomes negligible, so the premise fails. A side effect is that the model over-predicts the duration of benefit after a single induction course.")),
        tags$li("alpha is not identifiable from published HHT data, so the shunt-attributable share of any cardiac-index response is not identifiable either."),
        tags$li("The model does not reproduce Azzopardi's dosing-interval dependence, because full target occupancy at every interval leaves no mechanism for it. The discriminating measurement is trough free VEGF, which has not been reported."),
        tags$li("AVM formation is deterministic here. A real second hit is stochastic, so the model cannot predict the appearance of a new AVM after surgery or pregnancy."),
        tags$li("HHT-PAH (true pulmonary arterial hypertension in ACVRL1 carriers) is out of scope; the pulmonary hypertension in this model is the post-capillary, high-output kind only."),
        tags$li("The anticoagulation paradox is represented but not calibrated, so no scenario here quantifies it.")),
      h4("The experiment this model most wants run"),
      p(HTML(paste0("Reduce bleeding while holding iron and transfusion ",
        "<b>protocol-fixed</b>, and measure haemoglobin. If Hb still does not ",
        "rise, the exchangeability argument at the heart of this model is ",
        "wrong and something else is carrying the effect."))),
      hr(),
      p(class = "note", "Educational and research use only. Not validated for clinical decision-making, prescribing, or regulatory submission.")
    )
  })
}

shinyApp(ui, server)
