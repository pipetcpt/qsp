# =============================================================================
#  ts_shiny_app.R
#  THYROID STORM (갑상선 폭풍) — interactive QSP dashboard
# =============================================================================
#
#  The app is built around ONE question, which is also the model's thesis:
#
#      Two patients have the same total T4 and the same total T3.
#      One is a severe thyrotoxic. The other is in a storm.
#      What is different?
#
#  Tab 1 lets you set the hormone level and the precipitant INDEPENDENTLY, and
#  shows that the storm verdict follows the precipitant, not the hormone.
#  Every later tab is a consequence of that.
#
#  Requires: shiny, mrgsolve, dplyr, and ts_mrgsolve_model.R in the same folder.
#
#      shiny::runApp("ts_shiny_app.R")
#
# =============================================================================

suppressMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
})

source("ts_mrgsolve_model.R", local = TRUE)

# ---------------------------------------------------------------------------
#  cached baselines -- continuation in TRAb is slow, so memoise it
# ---------------------------------------------------------------------------
.BASE <- new.env(parent = emptyenv())
base_for <- function(trab) {
  key <- sprintf("%.3f", trab)
  if (is.null(.BASE[[key]])) .BASE[[key]] <- ts_thyrotoxic(trab)
  .BASE[[key]]
}

TT_LEVELS <- c("euthyroid (T3 ~1.8)"        = 1.0,
               "mild (T3 ~2.6)"             = 2.0,
               "moderate (T3 ~3.9)"         = 3.0,
               "FLORID (T3 ~5.4) — index"   = 4.19,
               "very severe (T3 ~7.3)"      = 6.0,
               "extreme (T3 ~9.2)"          = 9.0)

DRUGS <- c("PTU 600 mg then 250 mg q4h"            = "ptu",
           "Methimazole 40 mg then 25 mg q6h"      = "mmi",
           "Propranolol 80 mg q4h"                 = "pro",
           "Esmolol infusion (t½ 9 min)"           = "esm",
           "SSKI 250 mg iodide q6h"                = "iod",
           "Hydrocortisone 100 mg q8h"             = "hc",
           "Iopanoic acid 1 g q8h"                 = "iop",
           "Cholestyramine 4 g q6h"                = "chol",
           "Acetaminophen 650 mg q6h"              = "apap",
           "Aspirin 650 mg q4h  (CONTRAINDICATED)" = "asa",
           "Plasma exchange at 12 h"               = "tpe")

build_spec <- function(drugs, iodide_time, cool, fluid) {
  evs <- list()
  p <- list()
  if ("ptu"  %in% drugs) evs <- c(evs, list(ev_ptu()))
  if ("mmi"  %in% drugs) evs <- c(evs, list(ev_mmi()))
  if ("pro"  %in% drugs) evs <- c(evs, list(ev_pro()))
  if ("iod"  %in% drugs) evs <- c(evs, list(ev_iod(iodide_time)))
  if ("hc"   %in% drugs) evs <- c(evs, list(ev_hc()))
  if ("iop"  %in% drugs) evs <- c(evs, list(ev_iop()))
  if ("chol" %in% drugs) evs <- c(evs, list(ev_chol()))
  if ("apap" %in% drugs) evs <- c(evs, list(ev_apap()))
  if ("asa"  %in% drugs) evs <- c(evs, list(ev_asa()))
  if ("esm"  %in% drugs) p$ESMRATE <- 1.60
  if ("tpe"  %in% drugs) p$TPE_ON  <- 1
  p$COOL  <- cool
  p$FLUID <- fluid
  if (!length(evs)) evs <- list(ts_none())
  list(ev = do.call(c, evs), p = p)
}

GUIDELINE <- c("ptu", "iod", "pro", "hc")

# ---------------------------------------------------------------------------
#  UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  title = "Thyroid Storm QSP",
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, 'Helvetica Neue', Arial, sans-serif; }
    .thesis { background:#FFEBEE; border-left:6px solid #B00020; padding:10px 14px;
              margin-bottom:12px; }
    .verdict-storm { background:#B00020; color:white; padding:14px; font-size:20px;
                     font-weight:bold; text-align:center; border-radius:6px; }
    .verdict-ok    { background:#1B5E20; color:white; padding:14px; font-size:20px;
                     font-weight:bold; text-align:center; border-radius:6px; }
    .factorbox { background:#F1F8E9; border:1px solid #C5E1A5; padding:8px;
                 border-radius:4px; }
    .warn { color:#AD1457; font-weight:bold; }
  "))),

  h2("갑상선 폭풍 — Thyroid Storm QSP model"),
  div(class = "thesis",
      strong("Storm is a loop-gain disease, not a hormone-level disease."),
      br(),
      "Serum T4 and T3 do not separate storm from uncomplicated thyrotoxicosis. ",
      "Set the hormone level and the precipitant independently below and watch ",
      "which one moves the verdict. Every drug is then classified by ",
      em("which factor of the heat-balance ratio it multiplies"), "."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1 · Patient"),
      selectInput("trab", "Thyrotoxicosis severity (TRAb drive)",
                  choices = TT_LEVELS, selected = 4.19),
      sliderInput("prec", "Precipitant intensity (sepsis / surgery / …)",
                  min = 0, max = 2, value = 1.30, step = 0.05),
      sliderInput("pgefrac", "…of which febrile (PGE₂) fraction",
                  min = 0, max = 1, value = 0.6, step = 0.05),
      sliderInput("cr0", "Cardiac reserve at onset",
                  min = 0.3, max = 1, value = 1.0, step = 0.05),
      sliderInput("ar", "Adrenal reserve", min = 0.3, max = 1,
                  value = 1.0, step = 0.05),
      hr(),
      h4("2 · Treatment"),
      checkboxGroupInput("drugs", NULL, choices = DRUGS, selected = character(0)),
      actionButton("guide", "Load ATA/JTA bundle", class = "btn-sm"),
      actionButton("clear", "Clear", class = "btn-sm"),
      sliderInput("iodt", "Iodide timing relative to the thionamide (h)",
                  min = -4, max = 8, value = 1, step = 1),
      sliderInput("cool", "External cooling (fractional rise in h)",
                  min = 0, max = 3, value = 0, step = 0.1),
      sliderInput("fluid", "IV fluids (volume fraction / h)",
                  min = 0, max = 0.01, value = 0, step = 0.0005),
      hr(),
      sliderInput("hours", "Simulate (h)", min = 24, max = 336,
                  value = 168, step = 24),
      actionButton("run", "RUN", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      uiOutput("verdict"),
      br(),
      tabsetPanel(
        id = "tabs",

        # ---------------- 1. the controlled comparison ----------------
        tabPanel(
          "1 · Same hormone, two fates",
          br(),
          p(strong("This is the whole argument."), " The left column varies the ",
            "hormone level with NO precipitant. The right column varies the ",
            "precipitant at a FIXED hormone level. Only one of them crosses ",
            "into a storm."),
          fluidRow(
            column(6, h5("Hormone sweep — no precipitant ever"),
                   tableOutput("tblHormone")),
            column(6, h5("Precipitant sweep — hormone held constant"),
                   tableOutput("tblPrecip"))
          ),
          plotOutput("plotThreshold", height = "320px")
        ),

        # ---------------- 2. patient profile ----------------
        tabPanel(
          "2 · Patient profile",
          br(),
          fluidRow(
            column(6, h5("Baseline before the precipitant"),
                   tableOutput("tblBase")),
            column(6, h5("Where the hormone actually is"),
                   tableOutput("tblPools"),
                   p(class = "warn",
                     "Note the reservoir line: most of this patient's hormone ",
                     "is already made. That is why no synthesis blocker is an ",
                     "acute drug."))
          )
        ),

        # ---------------- 3. drug PK ----------------
        tabPanel(
          "3 · Drug PK",
          br(),
          p("Plasma concentrations and the target engagement each produces. ",
            "Note how fast PTU turns over (t½ 1.5 h) relative to the 2-month ",
            "reservoir it is trying to empty."),
          plotOutput("plotPK", height = "620px")
        ),

        # ---------------- 4. hormone axes ----------------
        tabPanel(
          "4 · Hormone axes: total vs FREE",
          br(),
          p(strong("The second axis."), " Total T3 is what the laboratory ",
            "reports. Free T3 is what the nucleus sees. A β-blocker moves the ",
            "second without moving the first, because the NEFA that displaces ",
            "T3 from TBG is produced by β-driven lipolysis."),
          plotOutput("plotHormone", height = "620px")
        ),

        # ---------------- 5. the heat balance ----------------
        tabPanel(
          "5 · The heat balance (Λ and its factors)",
          br(),
          div(class = "factorbox",
              withMathJax(helpText(
                "$$\\Lambda=\\frac{Q_{prod}}{Q_{loss}}=
                  \\frac{80\\,\\text{W}\\cdot M_{thy}\\cdot M_{temp}\\cdot
                  M_{sns}\\cdot M_{unc}}
                  {(8+52E)\\cdot Vol^{1.5}\\cdot(1+cool)\\cdot(T_c-33)}$$"))),
          p("Both sides are products of factors, and every drug multiplies one ",
            "or two of them. β-blockade is the only drug that lowers TWO ",
            "(M", tags$sub("sns"), " directly, and M", tags$sub("thy"),
            " through NEFA → displacement → free T3). Aspirin RAISES two."),
          plotOutput("plotLambda", height = "560px"),
          h5("Factor decomposition at 6 h"),
          tableOutput("tblFactors")
        ),

        # ---------------- 6. clinical endpoints ----------------
        tabPanel(
          "6 · Clinical endpoints",
          br(),
          p("Burch–Wartofsky is computed from the simulated signs, not asserted. ",
            "The dashed line is the diagnostic threshold of 45."),
          plotOutput("plotEndpoints", height = "620px"),
          tableOutput("tblEndpoints")
        ),

        # ---------------- 7. scenario comparison ----------------
        tabPanel(
          "7 · Scenario comparison",
          br(),
          p("The canonical monotherapies against your current patient. Read the ",
            "runaway column, then read the ΔtotalT3 column, and notice that ",
            "they do not agree."),
          tableOutput("tblScenarios"),
          plotOutput("plotScenarios", height = "420px")
        ),

        # ---------------- 8. biomarkers ----------------
        tabPanel(
          "8 · Biomarkers",
          br(),
          p("rT3 is the laboratory proof that a D1 block is working: it rises ",
            "on PTU and does not move on methimazole, because rT3 is cleared ",
            "mainly BY D1. NEFA is the (rarely measured) marker of ring A."),
          plotOutput("plotBio", height = "620px")
        ),

        # ---------------- 9. iodide ----------------
        tabPanel(
          "9 · Iodide: Wolff–Chaikoff vs Jod-Basedow",
          br(),
          p("The same drug with two opposite arms. Move the ",
            strong("iodide timing"), " slider and the ",
            strong("gland autoregulation"), " selector below."),
          selectInput("wcfail", "Gland autoregulation",
                      choices = c("normal (Wolff–Chaikoff intact)" = "normal",
                                  "partial failure" = "partial",
                                  "autonomous nodule / iodine-deficient" = "auto"),
                      selected = "normal", width = "420px"),
          plotOutput("plotIodide", height = "560px"),
          tableOutput("tblIodide")
        ),

        # ---------------- 10. safety ----------------
        tabPanel(
          "10 · Safety: the β-blockade trap",
          br(),
          p("β-blockade of a heart that depends on its rate for output can ",
            "precipitate shock. Move the ", strong("cardiac reserve"),
            " slider and compare propranolol (t½ 4 h) with esmolol (t½ 9 min)."),
          plotOutput("plotSafety", height = "520px"),
          tableOutput("tblSafety")
        ),

        # ---------------- 11. provenance ----------------
        tabPanel(
          "11 · What was fitted",
          br(),
          h4("Fitted to NORMAL physiology (20 numbers)"),
          p("total and free T4/T3, rT3, the 80/20 peripheral/thyroidal split of ",
            "T3 production, the T4 and T3 half-lives, the 60–90 day glandular ",
            "store, TSH, dietary iodide, BMR 80 W, Q10 = 2, core temperature ",
            "37.0, heart rate 70, NEFA 0.4, cortisol 400. The euthyroid heat ",
            "balance and iodine mass balance are algebraically consistent, not ",
            "fitted."),
          h4("Fitted to NON-STORM drug data (7 numbers)"),
          p("thionamide TPO potency, PTU D1 potency, propranolol β potency, ",
            "iodide release-inhibition Emax, NIS down-regulation time constant, ",
            "glucocorticoid and iopanoic D1/D2 potencies."),
          h4("Fitted to STORM data"),
          p(strong("Exactly one number"), " — the hazard scale, set so that the ",
            "untreated fulminant arm gives 7-day mortality 85%."),
          h4(class = "warn", "Therefore read as PREDICTIONS, not fits"),
          tags$ul(
            tags$li("no hormone level alone storms"),
            tags$li("the boundary is a threshold in the precipitant"),
            tags$li("thionamides cannot stop a storm (reservoir arithmetic)"),
            tags$li("β-blockade halves FREE T3 while total T3 barely moves"),
            tags$li("the glandular store falls from 66 d to ~18 d in florid Graves"),
            tags$li("the sign of iodide is set by the gland's autoregulation"),
            tags$li("aspirin worsens the storm through two independent factors")),
          h4(class = "warn", "Honest limitations"),
          tags$ul(
            tags$li(strong("Treated-arm mortality is under-predicted."),
                    " Modern series report ~8–25% with full treatment; this ",
                    "model gives well under 1%, because it contains no ",
                    "comorbidity, no multi-organ failure, no thromboembolism ",
                    "and no drug-toxicity deaths."),
            tags$li("The runaway is capped at 43 °C by a numerical ceiling; ",
                    "temperatures at the ceiling mean lethal hyperthermia, ",
                    "not a predicted measurement."),
            tags$li("One scalar precipitant axis only — real precipitants ",
                    "differ qualitatively."),
            tags$li("No plasma/extravascular split, so plasma exchange is ",
                    "represented only by its net whole-body removal."),
            tags$li("Clinical Wolff–Chaikoff escape is not reproduced within ",
                    "21 days at SSKI doses.")),
          hr(),
          p(em("Educational and research use only. Not for clinical decisions."))
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
#  SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$guide, updateCheckboxGroupInput(session, "drugs",
                                                     selected = GUIDELINE))
  observeEvent(input$guide, {
    updateSliderInput(session, "cool", value = 1.2)
    updateSliderInput(session, "fluid", value = 0.0045)
    updateSliderInput(session, "iodt", value = 1)
  })
  observeEvent(input$clear, {
    updateCheckboxGroupInput(session, "drugs", selected = character(0))
    updateSliderInput(session, "cool", value = 0)
    updateSliderInput(session, "fluid", value = 0)
  })

  `%||%` <- function(a, b) if (is.null(a)) b else a

  trab <- reactive(as.numeric(input$trab))

  wcparam <- reactive(switch(input$wcfail %||% "normal",
    normal  = list(),
    partial = list(ImWCrel = 0.40, ImWCorg = 0.45),
    auto    = list(ImWCrel = 0.15, ImWCorg = 0.20)))

  base <- reactive({
    b <- base_for(trab())
    b$CR <- input$cr0
    b
  })

  sim <- eventReactive(input$run, {
    # the thionamide is fixed at t = 0; a negative iodide time means the iodide
    # went in FIRST, which is the classic error
    spec <- build_spec(input$drugs, max(input$iodt, 0), input$cool, input$fluid)
    if (input$iodt < 0 && "iod" %in% input$drugs) {
      spec$ev <- c(spec$ev, ev(amt = 1970 / 17, cmt = "Ipl", time = 0,
                               ii = 6, addl = 27))
      spec$ev <- spec$ev[spec$ev$time != max(input$iodt, 0), ]
    }
    spec$p <- c(spec$p, wcparam(), list(AR = input$ar))
    withProgress(message = "integrating…", {
      ts_run(spec, base = base(), TRAb = trab(), hours = input$hours,
             prec = input$prec, pge = input$pgefrac * input$prec)
    })
  }, ignoreNULL = FALSE)

  # ------------------------- verdict banner -------------------------
  output$verdict <- renderUI({
    d <- sim()
    Tmax <- max(d$Tc); bw <- max(d$BWPS); mo <- 100 * max(d$mort)
    storm <- Tmax >= 40
    div(class = if (storm) "verdict-storm" else "verdict-ok",
        sprintf("%s   |   peak core temp %.2f °C   |   peak BWPS %.0f   |   %s   |   7-day mortality %.2f%%",
                if (storm) "THERMAL RUNAWAY — STORM NOT CONTROLLED"
                else "CONTROLLED — a stable operating point exists",
                Tmax, bw,
                if (bw >= 45) "BWPS ≥ 45 (storm by score)" else "BWPS < 45",
                mo))
  })

  # ------------------------- tab 1 -------------------------
  output$tblHormone <- renderTable({
    withProgress(message = "hormone sweep…",
                 ts_sweep_hormone(trabs = c(1, 2, 3, 4.19, 6, 9),
                                  hours = 72)) %>%
      transmute(TRAb, totT3 = round(totT3, 2), freeT3 = round(freeT3, 1),
                Tc = round(Tc, 2), HR = round(HR), BWPS,
                verdict = ifelse(storm, "STORM", "no storm"))
  }, digits = 2)

  output$tblPrecip <- renderTable({
    withProgress(message = "precipitant sweep…",
                 ts_sweep_precipitant(precs = c(0, .6, .9, 1.2, 1.3, 1.5),
                                      TRAb = trab(), hours = 168)) %>%
      transmute(Prec, Tmax = round(Tmax, 2), HRmax = round(HRmax),
                BWPSmax, mort = round(mort_pct, 2),
                verdict = ifelse(storm, "STORM", "no storm"))
  }, digits = 2)

  output$plotThreshold <- renderPlot({
    s <- withProgress(message = "threshold…",
                      ts_sweep_precipitant(
                        precs = seq(0, 2, by = 0.1), TRAb = trab(), hours = 168))
    par(mar = c(4.5, 4.5, 3, 4.5))
    plot(s$Prec, s$Tmax, type = "b", pch = 19, lwd = 2, col = "#B00020",
         xlab = "precipitant intensity (hormone level held constant)",
         ylab = "peak core temperature (°C)",
         main = "The boundary is a THRESHOLD in the precipitant, not a gradient")
    abline(h = 40, lty = 3)
    par(new = TRUE)
    plot(s$Prec, s$BWPSmax, type = "l", lwd = 2, col = "#1565C0",
         axes = FALSE, xlab = "", ylab = "")
    axis(4); mtext("peak Burch–Wartofsky", side = 4, line = 3)
    abline(h = 45, lty = 2, col = "#1565C0")
    legend("topleft", c("peak core temperature", "peak BWPS",
                        "storm thresholds"),
           col = c("#B00020", "#1565C0", "black"), lwd = 2,
           lty = c(1, 1, 3), bty = "n")
  })

  # ------------------------- tab 2 -------------------------
  output$tblBase <- renderTable({
    b <- base()
    d <- ts_mod %>% update(init = b) %>% param(TRAb = trab()) %>%
      mrgsim(end = 0) %>% as.data.frame()
    data.frame(
      quantity = c("total T4 (nmol/L)", "total T3 (nmol/L)",
                   "free T3 (pmol/L)", "displacement Fd",
                   "TR signal Sig", "TSH (mIU/L)",
                   "core temperature (°C)", "heart rate (bpm)",
                   "NEFA (mmol/L)", "β1 density (× normal)",
                   "cortisol (nmol/L)", "cardiac work index",
                   "Burch–Wartofsky (no precipitant)"),
      value = c(round(d$T4, 1), round(d$T3, 2), round(d$fT3, 2),
                round(d$Fd, 3), round(d$Sig, 3), signif(d$TSH, 3),
                round(d$Tc, 2), round(d$HR), round(d$NEFA, 2),
                round(b$Rb, 2), round(b$Cort), round(d$Widx, 2),
                round(d$BWPS)))
  })

  output$tblPools <- renderTable({
    b <- base()
    d <- ts_mod %>% update(init = b) %>% param(TRAb = trab()) %>%
      mrgsim(end = 0) %>% as.data.frame()
    data.frame(
      pool = c("colloid store S (nmol T4-eq)",
               "RESERVOIR — days of supply already made",
               "plasma T4 pool (nmol)", "plasma T3 pool (nmol)",
               "fraction of T3 bound to protein (%)",
               "fraction of T3 made outside the thyroid (%)"),
      value = c(round(d$S), round(d$storeD, 1),
                round(d$T4 * 10), round(d$T3 * 40),
                round(100 * (1 - 0.003 * d$Fd), 2), 80))
  })

  # ------------------------- tab 3 -------------------------
  output$plotPK <- renderPlot({
    d <- sim()
    par(mfrow = c(3, 2), mar = c(4, 4.2, 2.5, 1))
    pk <- list(c("Cptu", "PTU (mg/L)"), c("Cmmi", "methimazole (mg/L)"),
               c("Cpro", "propranolol (ng/mL)"), c("Ipl", "plasma iodide (µmol/L)"),
               c("Cgc", "glucocorticoid (mg/L)"), c("Casa", "salicylate (mg/dL)"))
    for (v in pk) {
      plot(d$time, d[[v[1]]], type = "l", lwd = 2, col = "#2E7D32",
           xlab = "time (h)", ylab = v[2], main = v[2])
    }
    par(mfrow = c(1, 1))
  })

  # ------------------------- tab 4 -------------------------
  output$plotHormone <- renderPlot({
    d <- sim()
    par(mfrow = c(3, 2), mar = c(4, 4.2, 2.5, 1))
    plot(d$time, d$T4, type = "l", lwd = 2, col = "#1565C0",
         xlab = "time (h)", ylab = "nmol/L", main = "TOTAL T4")
    plot(d$time, d$T3, type = "l", lwd = 2, col = "#1565C0",
         xlab = "time (h)", ylab = "nmol/L", main = "TOTAL T3 (what is measured)")
    plot(d$time, d$fT3, type = "l", lwd = 2, col = "#B00020",
         xlab = "time (h)", ylab = "pmol/L",
         main = "FREE T3 (what the nucleus sees)")
    plot(d$time, d$Fd, type = "l", lwd = 2, col = "#B00020",
         xlab = "time (h)", ylab = "× normal",
         main = "displacement Fd = free-fraction multiplier")
    plot(d$time, d$S, type = "l", lwd = 2, col = "#E65100",
         xlab = "time (h)", ylab = "nmol T4-eq",
         main = "the RESERVOIR (colloid store)")
    plot(d$time, d$Secr, type = "l", lwd = 2, col = "#00695C",
         xlab = "time (h)", ylab = "nmol/h",
         main = "secretion (release from the reservoir)")
    par(mfrow = c(1, 1))
  })

  # ------------------------- tab 5 -------------------------
  output$plotLambda <- renderPlot({
    d <- sim()
    par(mfrow = c(2, 3), mar = c(4, 4.2, 2.5, 1))
    plot(d$time, d$Lambda, type = "l", lwd = 2, col = "#B00020",
         xlab = "time (h)", ylab = "Q_prod / Q_loss",
         main = "Λ (1.0 = balanced)")
    abline(h = 1, lty = 3)
    for (v in list(c("M_thy", "M_thy  (hormone)"),
                   c("M_sns", "M_sns  (β / shivering)"),
                   c("M_unc", "M_unc  (salicylate uncoupling)"),
                   c("Eeff", "E  (heat-loss effector)"),
                   c("hcond", "h  (heat-loss coefficient, W/K)"))) {
      plot(d$time, d[[v[1]]], type = "l", lwd = 2, col = "#1B5E20",
           xlab = "time (h)", ylab = "", main = v[2])
    }
    par(mfrow = c(1, 1))
  })

  output$tblFactors <- renderTable({
    d <- sim(); i <- which.min(abs(d$time - 6))
    data.frame(factor = c("M_thy", "M_temp", "M_sns", "M_unc",
                          "E (effector)", "Vol^1.5", "1 + cool",
                          "h (W/K)", "Λ"),
               `at 6 h` = c(d$M_thy[i], d$M_temp[i], d$M_sns[i], d$M_unc[i],
                            d$Eeff[i], (pmin(d$Vol[i], 1))^1.5,
                            1 + input$cool, d$hcond[i], d$Lambda[i]),
               `moved by` = c("thionamide, iodide, D1/D2 blockers, β-BLOCKADE",
                              "temperature itself (Q10 = 2)",
                              "β-BLOCKADE, PGE₂",
                              "ASPIRIN raises this",
                              "acetaminophen raises, CNS injury destroys",
                              "IV fluids raise this",
                              "external cooling",
                              "product of E, Vol and cooling",
                              "the ratio of the two columns above"),
               check.names = FALSE)
  }, digits = 3)

  # ------------------------- tab 6 -------------------------
  output$plotEndpoints <- renderPlot({
    d <- sim()
    par(mfrow = c(3, 2), mar = c(4, 4.2, 2.5, 1))
    plot(d$time, d$Tc, type = "l", lwd = 2, col = "#B00020",
         xlab = "time (h)", ylab = "°C", main = "core temperature")
    abline(h = 40, lty = 3)
    plot(d$time, d$HR, type = "l", lwd = 2, col = "#BF360C",
         xlab = "time (h)", ylab = "bpm", main = "heart rate")
    plot(d$time, d$BWPS, type = "l", lwd = 2, col = "#F57F17",
         xlab = "time (h)", ylab = "points", main = "Burch–Wartofsky")
    abline(h = 45, lty = 2)
    plot(d$time, d$CNSx, type = "l", lwd = 2, col = "#4A148C",
         xlab = "time (h)", ylab = "0–1", main = "CNS dysfunction (ring C)")
    plot(d$time, d$Bili, type = "l", lwd = 2, col = "#F9A825",
         xlab = "time (h)", ylab = "mg/dL", main = "bilirubin (jaundice at 3)")
    abline(h = 3, lty = 3)
    plot(d$time, 100 * d$mort, type = "l", lwd = 2, col = "black",
         xlab = "time (h)", ylab = "%", main = "cumulative mortality")
    par(mfrow = c(1, 1))
  })

  output$tblEndpoints <- renderTable({
    d <- sim()
    ix <- sapply(c(0, 6, 12, 24, 48, 72, 168),
                 function(h) which.min(abs(d$time - h)))
    ix <- ix[ix <= nrow(d)]
    data.frame(`time (h)` = round(d$time[ix]),
               `core T` = round(d$Tc[ix], 2),
               HR = round(d$HR[ix]),
               BWPS = d$BWPS[ix],
               `total T3` = round(d$T3[ix], 2),
               `free T3` = round(d$fT3[ix], 2),
               bilirubin = round(d$Bili[ix], 2),
               `perfusion` = round(d$perf[ix], 2),
               `mortality %` = round(100 * d$mort[ix], 2),
               check.names = FALSE)
  })

  # ------------------------- tab 7 -------------------------
  output$tblScenarios <- renderTable({
    b <- base()
    arms <- list(
      `untreated`             = list(ev = ts_none(), p = list()),
      `PTU alone`             = list(ev = ev_ptu(), p = list()),
      `methimazole alone`     = list(ev = ev_mmi(), p = list()),
      `iodide (after PTU)`    = list(ev = c(ev_ptu(), ev_iod(1)), p = list()),
      `propranolol alone`     = list(ev = ev_pro(), p = list()),
      `esmolol alone`         = list(ev = ts_none(), p = list(ESMRATE = 1.6)),
      `cooling + fluids`      = list(ev = ts_none(),
                                     p = list(COOL = 1.2, FLUID = 0.0045)),
      `FULL BUNDLE`           = list(
          ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()),
          p = list(COOL = 1.2, FLUID = 0.0045)))
    withProgress(message = "scenarios…", {
      ref <- NULL
      bind_rows(lapply(names(arms), function(nm) {
        d <- ts_run(arms[[nm]], base = b, TRAb = trab(), hours = 168,
                    prec = input$prec, pge = input$pgefrac * input$prec)
        i <- which.min(abs(d$time - 24))
        r <- data.frame(arm = nm, Tmax = round(max(d$Tc), 2),
                        `BWPS 24h` = d$BWPS[i],
                        `totT3 24h` = round(d$T3[i], 2),
                        `freeT3 24h` = round(d$fT3[i], 2),
                        `HR 24h` = round(d$HR[i]),
                        runaway = ifelse(max(d$Tc) >= 40, "YES", "no"),
                        `mortality %` = round(100 * max(d$mort), 2),
                        check.names = FALSE)
        r
      }))
    })
  })

  output$plotScenarios <- renderPlot({
    b <- base()
    arms <- list(untreated = list(ev = ts_none(), p = list()),
                 `PTU alone` = list(ev = ev_ptu(), p = list()),
                 `propranolol alone` = list(ev = ev_pro(), p = list()),
                 `FULL BUNDLE` = list(
                     ev = c(ev_ptu(), ev_iod(1), ev_pro(), ev_hc()),
                     p = list(COOL = 1.2, FLUID = 0.0045)))
    cols <- c("#B00020", "#2E7D32", "#1565C0", "black")
    dl <- lapply(names(arms), function(nm)
      ts_run(arms[[nm]], base = b, TRAb = trab(), hours = 168,
             prec = input$prec, pge = input$pgefrac * input$prec))
    par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 1))
    for (v in list(c("Tc", "core temperature (°C)"),
                   c("T3", "TOTAL T3 (nmol/L)"),
                   c("fT3", "FREE T3 (pmol/L)"))) {
      rng <- range(unlist(lapply(dl, function(d) d[[v[1]]])))
      plot(NA, xlim = c(0, 168), ylim = rng, xlab = "time (h)",
           ylab = v[2], main = v[2])
      for (i in seq_along(dl)) lines(dl[[i]]$time, dl[[i]][[v[1]]],
                                     col = cols[i], lwd = 2)
      if (v[1] == "Tc") abline(h = 40, lty = 3)
    }
    legend("topright", names(arms), col = cols, lwd = 2, bty = "n", cex = 0.9)
    par(mfrow = c(1, 1))
  })

  # ------------------------- tab 8 -------------------------
  output$plotBio <- renderPlot({
    d <- sim()
    par(mfrow = c(3, 2), mar = c(4, 4.2, 2.5, 1))
    plot(d$time, d$rT3, type = "l", lwd = 2, col = "#00838F",
         xlab = "time (h)", ylab = "nmol/L",
         main = "reverse T3 — rises only if D1 is blocked")
    plot(d$time, 100 * d$inhD1, type = "l", lwd = 2, col = "#00838F",
         xlab = "time (h)", ylab = "%", main = "D1 inhibition")
    plot(d$time, d$NEFA, type = "l", lwd = 2, col = "#B00020",
         xlab = "time (h)", ylab = "mmol/L", main = "NEFA — the ring A marker")
    plot(d$time, d$Bsig, type = "l", lwd = 2, col = "#EF6C00",
         xlab = "time (h)", ylab = "× normal", main = "β signal")
    plot(d$time, d$TSH, type = "l", lwd = 2, col = "#4527A0", log = "y",
         xlab = "time (h)", ylab = "mIU/L", main = "TSH (log scale)")
    plot(d$time, d$Cort, type = "l", lwd = 2, col = "#455A64",
         xlab = "time (h)", ylab = "nmol/L", main = "endogenous cortisol")
    par(mfrow = c(1, 1))
  })

  # ------------------------- tab 9 -------------------------
  output$plotIodide <- renderPlot({
    b <- base()
    regs <- list(`normal (WC intact)` = list(),
                 `partial failure` = list(ImWCrel = 0.40, ImWCorg = 0.45),
                 `autonomous nodule` = list(ImWCrel = 0.15, ImWCorg = 0.20))
    cols <- c("#1B5E20", "#EF6C00", "#B00020")
    withProgress(message = "iodide regimes…", {
      dl <- lapply(regs, function(pp)
        ts_run(list(ev = ev_iod(0), p = pp), base = b, TRAb = trab(),
               hours = 24 * 21, prec = 0, pge = 0))
    })
    par(mfrow = c(2, 2), mar = c(4, 4.2, 2.5, 1))
    for (v in list(c("T4", "total T4 (nmol/L)"), c("S", "colloid store (nmol)"),
                   c("Synth", "synthesis (nmol/h)"), c("NIS", "NIS activity"))) {
      rng <- range(unlist(lapply(dl, function(d) d[[v[1]]])))
      plot(NA, xlim = c(0, 21), ylim = rng, xlab = "days",
           ylab = v[2], main = v[2])
      for (i in seq_along(dl)) lines(dl[[i]]$time / 24, dl[[i]][[v[1]]],
                                     col = cols[i], lwd = 2)
    }
    legend("topright", names(regs), col = cols, lwd = 2, bty = "n", cex = 0.85)
    par(mfrow = c(1, 1))
  })

  output$tblIodide <- renderTable({
    withProgress(message = "Wolff–Chaikoff…", ts_wolff_chaikoff(trab()))
  }, digits = 3)

  # ------------------------- tab 10 -------------------------
  output$plotSafety <- renderPlot({
    s <- withProgress(message = "reserve sweep…",
                      ts_sweep_reserve(TRAb = trab()))
    par(mfrow = c(1, 3), mar = c(4, 4.5, 2.5, 1))
    for (v in list(c("CO_nadir", "perfusion nadir"),
                   c("shock_max", "peak shock index"),
                   c("BWPS_24h", "BWPS at 24 h"))) {
      pp <- s[s$arm == "propranolol", ]; ee <- s[s$arm == "esmolol", ]
      plot(pp$CR0, pp[[v[1]]], type = "b", pch = 19, lwd = 2, col = "#B00020",
           ylim = range(c(pp[[v[1]]], ee[[v[1]]])),
           xlab = "cardiac reserve at onset", ylab = v[2], main = v[2])
      lines(ee$CR0, ee[[v[1]]], type = "b", pch = 17, lwd = 2, col = "#1565C0")
    }
    legend("topright", c("propranolol (t½ 4 h)", "esmolol (t½ 9 min)"),
           col = c("#B00020", "#1565C0"), lwd = 2, pch = c(19, 17), bty = "n")
    par(mfrow = c(1, 1))
  })

  output$tblSafety <- renderTable({
    withProgress(message = "reserve sweep…", ts_sweep_reserve(TRAb = trab()))
  }, digits = 3)
}

shinyApp(ui, server)
