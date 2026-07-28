##  Preterm Labor / Spontaneous Preterm Birth — QSP Shiny Dashboard
##  ============================================================================
##  Companion app for `ptl_mrgsolve_model.R`.
##
##  Run with:
##      library(shiny); library(mrgsolve)
##      shiny::runApp("ptl_shiny_app.R")
##
##  The app is organised around ONE argument, and every tab is a different way
##  of looking at it:
##
##      The delivery hazard is a PRODUCT of three limbs —
##          lambda = LAMBDA0 x f_contr x f_cerv x f_memb x f_term
##      — and every tocolytic in clinical use acts on f_contr alone. So the
##      benefit of tocolysis is arithmetically bounded by what the OTHER two
##      limbs are carrying, and its clinical value is not the days themselves
##      but what you put INTO those days: antenatal corticosteroids, magnesium
##      for the fetal brain, and transfer to a level-III unit.
##
##  Tabs:
##      1. Patient & pregnancy  — set up the case and the four ignition arms
##      2. Drug PK              — every drug's concentration-time profile
##      3. Uterine PD           — CAP programme, Ca2+, contraction index
##      4. Cervix & membranes   — the two limbs tocolysis does not touch
##      5. Hazard decomposition — the central argument, as a stacked picture
##      6. Clinical endpoints   — GA at delivery, RDS, IVH, NEC, sepsis, CP
##      7. Scenario comparison  — the 11 prebuilt scenarios side by side
##      8. Biomarkers           — fFN, amniotic IL-6, cervical length, FIRS
##      9. Safety ledger        — Mg toxicity, ductal risk, hypotension, etc.
##     10. Model notes          — assumptions, calibration anchors, disclaimer
##
##  Disclaimer: research / education only. Not a tocolysis protocol, not a
##  dosing calculator, not validated against patient data.
##  ============================================================================

library(shiny)
library(mrgsolve)

# ---------------------------------------------------------------------------
# Load the model and pull the $ENV helpers into this session
# ---------------------------------------------------------------------------
MODEL_FILE <- "ptl_mrgsolve_model.R"

mod <- mrgsolve::mread_cache("ptl", file = MODEL_FILE, quiet = TRUE)
helpers <- as.list(mrgsolve::env_get(mod))
list2env(helpers, environment())

PALETTE <- c("#2166ac", "#d6604d", "#1b7837", "#762a83", "#b8860b",
             "#00838f", "#c2185b", "#546e7a")

# ---------------------------------------------------------------------------
# Small plotting helper: multi-series line plot with a legend
# ---------------------------------------------------------------------------
lineplot <- function(dlist, xvar, yvar, xlab, ylab, main,
                     hline = NULL, hlab = NULL, log_y = FALSE) {
  xs <- unlist(lapply(dlist, function(d) d[[xvar]]))
  ys <- unlist(lapply(dlist, function(d) d[[yvar]]))
  ys <- ys[is.finite(ys)]
  if (!length(ys)) ys <- c(0, 1)
  yl <- range(c(ys, hline), na.rm = TRUE)
  if (diff(yl) == 0) yl <- yl + c(-1, 1) * max(1e-6, abs(yl[1]) * 0.1)
  plot(NA, xlim = range(xs, na.rm = TRUE), ylim = yl,
       xlab = xlab, ylab = ylab, main = main, las = 1,
       log = if (log_y) "y" else "")
  grid(col = "grey90")
  for (i in seq_along(dlist)) {
    lines(dlist[[i]][[xvar]], dlist[[i]][[yvar]],
          col = PALETTE[(i - 1) %% length(PALETTE) + 1], lwd = 2.2)
  }
  if (!is.null(hline)) {
    abline(h = hline, lty = 2, col = "grey40")
    if (!is.null(hlab)) mtext(hlab, side = 4, at = hline, cex = 0.65, las = 1)
  }
  if (length(dlist) > 1)
    legend("topright", names(dlist), col = PALETTE[seq_along(dlist)],
           lwd = 2.2, bty = "n", cex = 0.8)
}

# ===========================================================================
#  UI
# ===========================================================================
ui <- fluidPage(
  titlePanel("Preterm Labor — QSP / PK-PD Explorer"),
  tags$p(style = "color:#555; margin-top:-8px;",
         tags$em(paste("Four ignition arms, three effector limbs, one",
                       "multiplicative hazard. Research/education only —",
                       "not a tocolysis protocol."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Pregnancy"),
      sliderInput("ga0", "Gestational age at t = 0 (wk)",
                  min = 18, max = 36, value = 29, step = 0.5),
      sliderInput("cl0", "Cervical length at t = 0 (mm)",
                  min = 5, max = 45, value = 14, step = 1),
      sliderInput("horizon", "Simulation horizon (days)",
                  min = 3, max = 140, value = 21, step = 1),
      selectInput("nfetus", "Fetuses", c("Singleton" = 1, "Twins" = 2),
                  selected = 1),
      checkboxInput("prior", "Prior spontaneous preterm birth", FALSE),
      checkboxInput("pprom", "Membranes already ruptured (PPROM)", FALSE),
      sliderInput("crcl", "Maternal CrCl (mL/min)",
                  min = 20, max = 180, value = 130, step = 5),

      hr(),
      h4("Ignition arms"),
      sliderInput("clock", "CRH clock advancement (wk)",
                  min = 0, max = 6, value = 0, step = 0.5),
      sliderInput("infect", "Intra-amniotic microbial inoculum",
                  min = 0, max = 1, value = 0, step = 0.05),
      sliderInput("sterile", "Sterile (DAMP) inflammatory drive",
                  min = 0, max = 1, value = 0.35, step = 0.05),
      sliderInput("thromb", "Decidual haemorrhage / thrombin",
                  min = 0, max = 1, value = 0.15, step = 0.05),
      sliderInput("fragile", "Intrinsic cervical insufficiency",
                  min = 0, max = 1, value = 0, step = 0.05),
      sliderInput("ottone", "Endogenous oxytocin tone (fold)",
                  min = 0.5, max = 3, value = 1.6, step = 0.1),

      hr(),
      h4("Interventions"),
      checkboxGroupInput(
        "drugs", NULL,
        choices = c("Vaginal progesterone 200 mg qhs" = "p4v",
                    "17-OHPC 250 mg IM weekly"        = "ohpc",
                    "Atosiban (bolus + 48 h)"         = "ato",
                    "Nifedipine 20 mg q6h x 48 h"     = "nif",
                    "Indomethacin 48 h"               = "ind",
                    "MgSO4 (load + infusion)"         = "mg",
                    "Terbutaline 0.25 mg SC q4h"      = "terb",
                    "Betamethasone 12 mg x 2"         = "bet",
                    "Erythromycin (PPROM)"            = "ery",
                    "Cerclage"                        = "cerc",
                    "Arabin pessary"                  = "pess"),
        selected = c("ato", "bet")),
      sliderInput("mg_rate", "MgSO4 maintenance (g/h)",
                  min = 0.5, max = 3, value = 1, step = 0.5),
      sliderInput("mg_hours", "MgSO4 duration (h)",
                  min = 6, max = 72, value = 24, step = 6),
      sliderInput("toco_h", "Oral tocolytic duration (h)",
                  min = 6, max = 168, value = 48, step = 6),

      hr(),
      actionButton("reset", "Reset to default acute case", class = "btn-sm")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        # ---- 1 --------------------------------------------------------
        tabPanel(
          "1 · Patient",
          br(),
          fluidRow(
            column(6, h4("Case summary"), verbatimTextOutput("caseTxt")),
            column(6, h4("Baseline risk picture"), plotOutput("riskPlot",
                                                              height = 300))
          ),
          hr(),
          h4("Endocrine ignition — the placental clock"),
          plotOutput("clockPlot", height = 320)
        ),

        # ---- 2 --------------------------------------------------------
        tabPanel(
          "2 · Drug PK",
          br(),
          fluidRow(
            column(6, plotOutput("pkToco",  height = 300)),
            column(6, plotOutput("pkMg",    height = 300))
          ),
          fluidRow(
            column(6, plotOutput("pkBet",   height = 300)),
            column(6, plotOutput("pkP4",    height = 300))
          ),
          tags$p(style = "color:#555;",
                 tags$em(paste0(
                   "Note the vaginal-progesterone panel: plasma concentration ",
                   "is unimpressive while the LOCAL (first-uterine-pass) ",
                   "exposure is an order of magnitude higher. Plasma is the ",
                   "wrong exposure metric for this route.")))
        ),

        # ---- 3 --------------------------------------------------------
        tabPanel(
          "3 · Uterine PD",
          br(),
          fluidRow(
            column(6, plotOutput("pdCAP",   height = 300)),
            column(6, plotOutput("pdCa",    height = 300))
          ),
          fluidRow(
            column(6, plotOutput("pdContr", height = 300)),
            column(6, plotOutput("pdPRW",   height = 300))
          )
        ),

        # ---- 4 --------------------------------------------------------
        tabPanel(
          "4 · Cervix & membranes",
          br(),
          tags$p(tags$strong(paste0(
            "These are the two limbs no tocolytic touches. Watch them keep ",
            "running while the contraction index falls."))),
          fluidRow(
            column(6, plotOutput("cxLen",  height = 300)),
            column(6, plotOutput("cxMMP",  height = 300))
          ),
          fluidRow(
            column(6, plotOutput("mbTens", height = 300)),
            column(6, plotOutput("cxColl", height = 300))
          )
        ),

        # ---- 5 --------------------------------------------------------
        tabPanel(
          "5 · Hazard decomposition",
          br(),
          tags$p(tags$strong(
            "lambda = LAMBDA0 x f_contr x f_cerv x f_memb x f_term")),
          tags$p(style = "color:#555;", tags$em(paste0(
            "A tocolytic moves f_contr and nothing else. The table below ",
            "prints the arithmetic floor: the fraction of the untreated ",
            "hazard that survives even if contractions are abolished ",
            "completely."))),
          fluidRow(
            column(7, plotOutput("hazLimbs", height = 340)),
            column(5, plotOutput("hazLambda", height = 340))
          ),
          hr(),
          h4("Limb decomposition across tocolytic arms (at 24 h)"),
          tableOutput("limbTable")
        ),

        # ---- 6 --------------------------------------------------------
        tabPanel(
          "6 · Clinical endpoints",
          br(),
          fluidRow(
            column(6, plotOutput("epPdel", height = 300)),
            column(6, plotOutput("epSurf", height = 300))
          ),
          hr(),
          h4("Expected neonatal outcomes at the model's expected GA at delivery"),
          tableOutput("outcomeTable"),
          tags$p(style = "color:#555;", tags$em(paste0(
            "Outcome probabilities come from gestational-age logistics ",
            "calibrated to Stoll 2015 JAMA and Manuck 2016 AJOG, shifted by ",
            "the simulated lung-maturity index, fetal inflammatory state and ",
            "magnesium exposure.")))
        ),

        # ---- 7 --------------------------------------------------------
        tabPanel(
          "7 · Scenario comparison",
          br(),
          tags$p(paste0(
            "The eleven prebuilt scenarios. Scenario 10 is the one to read ",
            "carefully: tocolysis gains days AND worsens the composite ",
            "outcome, because the days are bought inside an inflamed ",
            "amniotic cavity.")),
          actionButton("runScen", "Run all 11 scenarios", class = "btn-primary"),
          br(), br(),
          fluidRow(
            column(6, plotOutput("scPrev",  height = 320)),
            column(6, plotOutput("scAcute", height = 320))
          ),
          hr(),
          h4("Ledger — days gained vs OUTCOME gained"),
          tableOutput("ledgerTable")
        ),

        # ---- 8 --------------------------------------------------------
        tabPanel(
          "8 · Biomarkers",
          br(),
          fluidRow(
            column(6, plotOutput("bmFFN", height = 300)),
            column(6, plotOutput("bmIL6", height = 300))
          ),
          fluidRow(
            column(6, plotOutput("bmFIRS", height = 300)),
            column(6, plotOutput("bmBact", height = 300))
          )
        ),

        # ---- 9 --------------------------------------------------------
        tabPanel(
          "9 · Safety ledger",
          br(),
          fluidRow(
            column(6, plotOutput("safeMg",   height = 300)),
            column(6, plotOutput("safeDuct", height = 300))
          ),
          fluidRow(
            column(6, plotOutput("safeMap",  height = 300)),
            column(6, plotOutput("safeB2",   height = 300))
          ),
          hr(),
          h4("Tocolytic head-to-head — efficacy and harm in the same table"),
          actionButton("runH2H", "Run head-to-head", class = "btn-primary"),
          br(), br(),
          tableOutput("h2hTable")
        ),

        # ---- 10 -------------------------------------------------------
        tabPanel(
          "10 · Model notes",
          br(),
          htmlOutput("notes")
        )
      )
    )
  )
)

# ===========================================================================
#  SERVER
# ===========================================================================
server <- function(input, output, session) {

  observeEvent(input$reset, {
    updateSliderInput(session, "ga0", value = 29)
    updateSliderInput(session, "cl0", value = 14)
    updateSliderInput(session, "horizon", value = 21)
    updateSliderInput(session, "sterile", value = 0.35)
    updateSliderInput(session, "thromb", value = 0.15)
    updateSliderInput(session, "infect", value = 0)
    updateSliderInput(session, "ottone", value = 1.6)
    updateCheckboxGroupInput(session, "drugs", selected = c("ato", "bet"))
  })

  # -- parameter set from the sidebar ---------------------------------------
  parlist <- reactive({
    list(GA0_WK       = input$ga0,
         CL0          = input$cl0,
         NFETUS       = as.numeric(input$nfetus),
         PRIOR_PTB    = as.numeric(input$prior),
         PPROM_FLAG   = as.numeric(input$pprom),
         CRCL         = input$crcl,
         CLOCK_ADV    = input$clock,
         INFECT_DRIVE = input$infect,
         BACTBASE     = if (input$infect > 0) 2.0 else 0.0,
         STERILE_DRIVE= input$sterile,
         THROMBIN_DRV = input$thromb,
         CERV_FRAGILE = input$fragile,
         OT_TONE      = input$ottone,
         ERY_ON       = as.numeric("ery"  %in% input$drugs),
         CERC_SUPPORT = as.numeric("cerc" %in% input$drugs),
         PESS_SUPPORT = as.numeric("pess" %in% input$drugs))
  })

  # -- dosing events from the sidebar ---------------------------------------
  evlist <- reactive({
    d <- input$drugs
    e <- list()
    if ("p4v"  %in% d) e <- c(e, list(ev_vaginal_P4(0, input$horizon)))
    if ("ohpc" %in% d) e <- c(e, list(ev_17OHPC(0, input$horizon)))
    if ("ato"  %in% d) e <- c(e, list(ev_atosiban(0)))
    if ("nif"  %in% d) e <- c(e, list(ev_nifedipine(0, input$toco_h)))
    if ("ind"  %in% d) e <- c(e, list(ev_indomethacin(0, input$toco_h)))
    if ("mg"   %in% d) e <- c(e, list(ev_magnesium(0, input$mg_hours,
                                                   4, input$mg_rate)))
    if ("terb" %in% d) e <- c(e, list(ev_terbutaline(0, input$toco_h)))
    if ("bet"  %in% d) e <- c(e, list(ev_betamethasone(0)))
    if (!length(e)) return(ev_none())
    Reduce(c, e)
  })

  # -- the two runs we always keep side by side -----------------------------
  sim <- reactive({
    m  <- mrgsolve::param(mod, parlist())
    dl <- max(0.01, input$horizon / 900)
    treated <- as.data.frame(
      mrgsolve::mrgsim_e(m, evlist(), end = input$horizon, delta = dl))
    untreated <- as.data.frame(
      mrgsolve::mrgsim_e(mrgsolve::param(m, ERY_ON = 0, CERC_SUPPORT = 0,
                                         PESS_SUPPORT = 0),
                         ev_none(), end = input$horizon, delta = dl))
    list("selected regimen" = treated, "no intervention" = untreated)
  })

  # =========================== 1. Patient ==================================
  output$caseTxt <- renderText({
    p <- parlist()
    paste0(
      sprintf("Gestational age at t = 0 : %.1f weeks\n", p$GA0_WK),
      sprintf("Cervical length          : %.0f mm%s\n", p$CL0,
              if (p$CL0 < 25) "   (SHORT — below the 25 mm threshold)" else ""),
      sprintf("Fetuses                  : %d\n", as.integer(p$NFETUS)),
      sprintf("Prior sPTB               : %s\n", if (p$PRIOR_PTB > 0) "yes" else "no"),
      sprintf("Membranes                : %s\n",
              if (p$PPROM_FLAG > 0) "RUPTURED (PPROM)" else "intact"),
      sprintf("Maternal CrCl            : %.0f mL/min%s\n", p$CRCL,
              if (p$CRCL < 60)
                "   (! magnesium clearance is reduced proportionally)" else ""),
      "\nIgnition arms\n",
      sprintf("  CRH clock advancement  : %.1f wk\n", p$CLOCK_ADV),
      sprintf("  Microbial inoculum     : %.2f%s\n", p$INFECT_DRIVE,
              if (p$INFECT_DRIVE > 0.3)
                "   (! tocolysis buys days inside an infected cavity)" else ""),
      sprintf("  Sterile inflammation   : %.2f\n", p$STERILE_DRIVE),
      sprintf("  Thrombin / haemorrhage : %.2f\n", p$THROMBIN_DRV),
      sprintf("  Cervical insufficiency : %.2f\n", p$CERV_FRAGILE),
      sprintf("  Oxytocin tone          : %.1f x", p$OT_TONE),
      if (p$OT_TONE > 2)
        "   (! high tone erodes competitive OXTR blockade)" else "")
  })

  output$riskPlot <- renderPlot({
    p <- parlist()
    v <- c("CRH clock"   = p$CLOCK_ADV / 6,
           "Infection"   = p$INFECT_DRIVE,
           "Sterile infl"= p$STERILE_DRIVE,
           "Thrombin"    = p$THROMBIN_DRV,
           "Stretch"     = (p$NFETUS - 1) * 0.55,
           "Cervix"      = max(p$CERV_FRAGILE, (25 - p$CL0) / 25),
           "PPROM"       = p$PPROM_FLAG)
    v[v < 0] <- 0
    op <- par(mar = c(4, 8, 2, 1)); on.exit(par(op))
    barplot(rev(v), horiz = TRUE, las = 1, xlim = c(0, 1),
            col = rev(PALETTE[seq_along(v)]), border = NA,
            xlab = "Relative driver strength (0-1)",
            main = "Ignition-arm profile")
  })

  output$clockPlot <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1)); on.exit(par(op))
    lineplot(s, "GAW", "CRH", "GA (wk)", "pg/mL", "Placental CRH")
    lineplot(s, "GAW", "E3",  "GA (wk)", "fold",  "Estriol (E3)")
    lineplot(s, "GAW", "CORTF", "GA (wk)", "fold", "Fetal cortisol")
  })

  # =========================== 2. Drug PK ==================================
  output$pkToco <- renderPlot({
    d <- sim()[["selected regimen"]]
    op <- par(mar = c(4, 4, 3, 4)); on.exit(par(op))
    plot(d$time, d$CATOo, type = "l", lwd = 2.2, col = PALETTE[1],
         xlab = "Days", ylab = "ng/mL", las = 1,
         main = "Tocolytic plasma concentrations",
         ylim = c(0, max(1, d$CATOo, d$CNIFo, d$CINDo, na.rm = TRUE)))
    grid(col = "grey90")
    lines(d$time, d$CNIFo, lwd = 2.2, col = PALETTE[2])
    lines(d$time, d$CINDo, lwd = 2.2, col = PALETTE[3])
    lines(d$time, d$CTERBo, lwd = 2.2, col = PALETTE[4])
    legend("topright", c("atosiban", "nifedipine", "indomethacin", "terbutaline"),
           col = PALETTE[1:4], lwd = 2.2, bty = "n", cex = 0.8)
  })

  output$pkMg <- renderPlot({
    d <- sim()[["selected regimen"]]
    plot(d$time, d$CMGo, type = "l", lwd = 2.4, col = PALETTE[5],
         xlab = "Days", ylab = "mg/dL", las = 1,
         main = "Serum magnesium (with the toxicity ladder)",
         ylim = c(0, max(8, max(d$CMGo, na.rm = TRUE) * 1.1)))
    grid(col = "grey90")
    rect(par("usr")[1], 4, par("usr")[2], 7,
         col = adjustcolor("#4caf50", 0.12), border = NA)
    abline(h = c(4, 7, 9, 12), lty = c(3, 3, 2, 2),
           col = c("grey50", "grey50", "#e65100", "#b71c1c"))
    text(par("usr")[1], 9,  " DTR loss",            adj = 0, cex = 0.7, col = "#e65100")
    text(par("usr")[1], 12, " respiratory depression", adj = 0, cex = 0.7, col = "#b71c1c")
    text(par("usr")[1], 5.5, " therapeutic 4-7",     adj = 0, cex = 0.7, col = "#2e7d32")
    lines(d$time, d$CMGo, lwd = 2.4, col = PALETTE[5])
  })

  output$pkBet <- renderPlot({
    d <- sim()[["selected regimen"]]
    op <- par(mar = c(4, 4, 3, 4)); on.exit(par(op))
    plot(d$time, d$CBETo, type = "l", lwd = 2.2, col = PALETTE[1],
         xlab = "Days", ylab = "ng/mL", las = 1,
         main = "Betamethasone: maternal vs fetal",
         ylim = c(0, max(1, d$CBETo, na.rm = TRUE)))
    grid(col = "grey90")
    lines(d$time, d$BET_F, lwd = 2.2, col = PALETTE[2])
    legend("topright", c("maternal plasma", "fetal plasma"),
           col = PALETTE[1:2], lwd = 2.2, bty = "n", cex = 0.8)
    mtext("betamethasone escapes placental 11b-HSD2 — cortisol does not",
          side = 3, cex = 0.7, col = "grey40")
  })

  output$pkP4 <- renderPlot({
    d <- sim()[["selected regimen"]]
    plot(d$time, d$CP4LOCo, type = "l", lwd = 2.4, col = PALETTE[3],
         xlab = "Days", ylab = "ng/mL (or ng/mL-eq)", las = 1,
         main = "Progestogen: plasma vs LOCAL uterine exposure",
         ylim = c(0, max(1, d$CP4LOCo, d$COHPCo, na.rm = TRUE)))
    grid(col = "grey90")
    lines(d$time, d$CP4o,   lwd = 2.2, col = PALETTE[1])
    lines(d$time, d$COHPCo, lwd = 2.2, col = PALETTE[2])
    legend("topright",
           c("vaginal P4 — LOCAL (first uterine pass)",
             "vaginal P4 — plasma", "17-OHPC — plasma"),
           col = PALETTE[c(3, 1, 2)], lwd = 2.2, bty = "n", cex = 0.75)
  })

  # =========================== 3. Uterine PD ===============================
  output$pdCAP <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    lineplot(s, "time", "OXTR", "Days", "fold", "OXTR density")
    lineplot(s, "time", "CX43", "Days", "fold", "Connexin-43")
  })
  output$pdCa <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    lineplot(s, "time", "CAI",  "Days", "fold", "Myometrial Ca2+ index")
    lineplot(s, "time", "CAMP", "Days", "fold", "Myometrial cAMP")
  })
  output$pdContr <- renderPlot({
    s <- sim()
    lineplot(s, "time", "CONTR", "Days", "MVU-like index",
             "Contraction index — the ONLY limb tocolysis moves")
  })
  output$pdPRW <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    lineplot(s, "time", "PRW", "Days", "0-1",
             "Functional P4 withdrawal (PR-A/PR-B)")
    lineplot(s, "time", "PG", "Days", "fold", "Prostaglandin pool")
  })

  # ===================== 4. Cervix & membranes =============================
  output$cxLen <- renderPlot({
    s <- sim()
    lineplot(s, "time", "CLEN", "Days", "mm", "Cervical length",
             hline = 25, hlab = "25 mm")
  })
  output$cxMMP <- renderPlot({
    s <- sim()
    lineplot(s, "time", "MMP", "Days", "fold", "Cervical / membrane MMP activity")
  })
  output$mbTens <- renderPlot({
    s <- sim()
    lineplot(s, "time", "MEMB", "Days", "fold",
             "Membrane tensile strength", hline = 0.30, hlab = "rupture zone")
  })
  output$cxColl <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
    lineplot(s, "time", "COLL", "Days", "fold", "Cross-linked collagen")
    lineplot(s, "time", "HA",   "Days", "0-1",  "Hyaluronan / hydration")
  })

  # ===================== 5. Hazard decomposition ===========================
  output$hazLimbs <- renderPlot({
    d <- sim()[["selected regimen"]]
    plot(NA, xlim = range(d$time),
         ylim = c(0, max(d$FCONTR, d$FCERV, d$FMEMB, na.rm = TRUE) * 1.05),
         xlab = "Days", ylab = "hazard factor", las = 1,
         main = "The three limbs of the hazard")
    grid(col = "grey90")
    lines(d$time, d$FCONTR, lwd = 2.6, col = PALETTE[2])
    lines(d$time, d$FCERV,  lwd = 2.6, col = PALETTE[4])
    lines(d$time, d$FMEMB,  lwd = 2.6, col = PALETTE[3])
    legend("topright",
           c("f_contr  (tocolytics act HERE, and only here)",
             "f_cerv   (untouched)", "f_memb  (untouched)"),
           col = PALETTE[c(2, 4, 3)], lwd = 2.6, bty = "n", cex = 0.8)
  })

  output$hazLambda <- renderPlot({
    s <- sim()
    lineplot(s, "time", "LAMBDA", "Days", "1/day",
             "Instantaneous delivery hazard")
  })

  output$limbTable <- renderTable({
    out <- PTL_limb_decomposition(mod, ga0 = input$ga0, cl0 = input$cl0)
    out
  }, digits = 3)

  # ===================== 6. Clinical endpoints =============================
  output$epPdel <- renderPlot({
    s <- sim()
    lineplot(s, "time", "PDELIV", "Days", "probability",
             "Cumulative probability of delivery")
  })
  output$epSurf <- renderPlot({
    s <- sim()
    lineplot(s, "time", "SURF", "Days", "0-1",
             "Fetal lung-maturity index")
  })

  output$outcomeTable <- renderTable({
    s <- sim()
    rows <- lapply(names(s), function(nm) {
      r <- PTL_summarise_run(s[[nm]])
      o <- PTL_neonatal_outcomes(r$E_GA_wk, surf = r$surf_at_birth,
                                 firs = r$FIRS_at_birth, mgbrain = r$MgBrain)
      data.frame(arm = nm, E_GA_wk = r$E_GA_wk,
                 P_48h = r$P_deliver_48h, P_7d = r$P_deliver_7d,
                 surf = r$surf_at_birth, FIRS = r$FIRS_at_birth,
                 RDS = o$RDS, severe_IVH = o$severe_IVH, NEC = o$NEC,
                 sepsis = o$sepsis, death = o$death, CP = o$CP,
                 composite = o$composite, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }, digits = 3)

  # ===================== 7. Scenario comparison ============================
  scen <- eventReactive(input$runScen, {
    PTL_simulate_scenarios(mod)
  })

  output$scPrev <- renderPlot({
    sim_all <- scen()
    keep <- c("2_short_cervix_untreated", "3_short_cervix_vaginalP4",
              "5_cerclage_plus_P4")
    dl <- setNames(lapply(keep, function(k) sim_all[sim_all$scenario == k, ]),
                   c("untreated", "vaginal P4", "cerclage + P4"))
    lineplot(dl, "GAW", "CLEN", "GA (wk)", "mm",
             "Prevention arm — cervical length", hline = 25)
  })

  output$scAcute <- renderPlot({
    sim_all <- scen()
    keep <- c("6_acute_PTL_untreated", "7_acute_PTL_atosiban_ACS",
              "10_intraamniotic_infection")
    dl <- setNames(lapply(keep, function(k) {
      d <- sim_all[sim_all$scenario == k, ]; d[d$time <= 14, ]
    }), c("untreated", "atosiban + ACS", "MIAC + atosiban + ACS"))
    lineplot(dl, "time", "PDELIV", "Days", "probability",
             "Acute arm — cumulative delivery probability")
  })

  output$ledgerTable <- renderTable({
    req(input$runScen)
    led <- PTL_ledger(mod)
    led[, c("scenario", "E_GA_wk", "days_gained", "P_deliver_48h",
            "P_deliver_7d", "surf_at_birth", "FIRS_at_birth",
            "RDS", "severe_IVH", "sepsis", "CP", "composite")]
  }, digits = 3)

  # ===================== 8. Biomarkers =====================================
  output$bmFFN <- renderPlot({
    s <- sim()
    lineplot(s, "time", "FFN", "Days", "ng/mL",
             "Cervicovaginal fetal fibronectin", hline = 50, hlab = "50")
  })
  output$bmIL6 <- renderPlot({
    s <- sim()
    lineplot(s, "time", "IL6", "Days", "ng/mL",
             "Amniotic fluid IL-6", hline = 2.6, hlab = "2.6")
  })
  output$bmFIRS <- renderPlot({
    s <- sim()
    lineplot(s, "time", "FIRS", "Days", "0-1",
             "Fetal inflammatory response index")
  })
  output$bmBact <- renderPlot({
    s <- sim()
    lineplot(s, "time", "BACT", "Days", "log10 CFU/mL",
             "Intra-amniotic bacterial load")
  })

  # ===================== 9. Safety ledger ==================================
  output$safeMg <- renderPlot({
    d <- sim()[["selected regimen"]]
    plot(d$time, d$CMGo, type = "l", lwd = 2.4, col = PALETTE[5],
         xlab = "Days", ylab = "mg/dL", las = 1,
         main = sprintf("Serum Mg (CrCl %.0f mL/min)", input$crcl),
         ylim = c(0, max(8, max(d$CMGo, na.rm = TRUE) * 1.1)))
    grid(col = "grey90")
    abline(h = c(7, 9, 12), lty = 2, col = c("grey50", "#e65100", "#b71c1c"))
  })
  output$safeDuct <- renderPlot({
    d <- sim()[["selected regimen"]]
    op <- par(mar = c(4, 4, 3, 1)); on.exit(par(op))
    plot(d$time, d$DUCTRISK, type = "l", lwd = 2.4, col = PALETTE[2],
         xlab = "Days", ylab = "risk index", las = 1,
         main = sprintf("Indomethacin ductal-constriction risk (GA %.1f wk)",
                        input$ga0),
         ylim = c(0, max(0.05, max(d$DUCTRISK, na.rm = TRUE) * 1.1)))
    grid(col = "grey90")
    mtext("susceptibility rises steeply after 32 weeks", side = 3,
          cex = 0.7, col = "grey40")
  })
  output$safeMap <- renderPlot({
    s <- sim()
    lineplot(s, "time", "MAP_DROP", "Days", "mmHg",
             "Nifedipine-associated MAP reduction")
  })
  output$safeB2 <- renderPlot({
    s <- sim()
    lineplot(s, "time", "B2_LOST", "Days", "fraction",
             "Beta2 responsiveness lost (tachyphylaxis)")
  })

  h2h <- eventReactive(input$runH2H, {
    PTL_tocolytic_head_to_head(mrgsolve::param(mod, CRCL = input$crcl),
                               ga0 = input$ga0)
  })
  output$h2hTable <- renderTable({ h2h() }, digits = 3)

  # ===================== 10. Model notes ===================================
  output$notes <- renderUI({
    HTML(paste0(
      "<h4>What this model claims</h4>",
      "<p>The delivery hazard is written as a <b>product</b> of three limbs: ",
      "<code>lambda = LAMBDA0 &times; f_contr &times; f_cerv &times; f_memb ",
      "&times; f_term</code>. Every tocolytic in clinical use acts on ",
      "<code>f_contr</code> alone — atosiban at OXTR, nifedipine at Cav1.2, ",
      "indomethacin at COX, magnesium at the myofilament. Four targets, one ",
      "limb. So the achievable reduction in hazard is bounded by what the ",
      "cervical and membrane limbs are already carrying, and that bound is ",
      "printed as a number on the <i>Hazard decomposition</i> tab.</p>",

      "<p>The corollary is the clinically important half: the interventions ",
      "that change outcome do not treat the uterus at all. Antenatal ",
      "corticosteroids mature the fetal lung, magnesium protects the fetal ",
      "brain, and transfer moves the baby to a level-III unit. Tocolysis is ",
      "worth exactly the 48 hours those three need — and no more, which is ",
      "why maintenance tocolysis has never shown benefit.</p>",

      "<h4>The structural choice you should argue with first</h4>",
      "<p>The multiplicative hazard is a <b>modelling choice, not a measured ",
      "fact</b>. It was chosen because an additive hazard would let a ",
      "sufficiently strong tocolytic drive lambda to zero — precisely the ",
      "prediction forty years of trials have falsified. Anyone re-using this ",
      "model should treat that choice, rather than any parameter value, as ",
      "its principal assumption.</p>",

      "<h4>Calibration anchors</h4><ul>",
      "<li>McLean 1995 <i>Nat Med</i> — CRH doubling time 3.3 wk ",
      "(the placental clock)</li>",
      "<li>Mesiano 2002 <i>JCEM</i> / Merlino 2007 — functional progesterone ",
      "withdrawal is a PR-A/PR-B ratio shift, not a fall in plasma P4</li>",
      "<li>Iams 1996 <i>NEJM</i> — mid-trimester cervical-length distribution ",
      "and the continuous CL-to-risk gradient</li>",
      "<li>Romero 2014 <i>AJRI</i> — sterile intra-amniotic inflammation is ",
      "about twice as common as culture-proven MIAC</li>",
      "<li>EPPPIC 2021 <i>Lancet</i> — vaginal progesterone RR ~0.78 for ",
      "birth &lt; 34 wk in singletons at risk</li>",
      "<li>Goodwin 1995 — atosiban PK (Vc 18 L, CL 42 L/h)</li>",
      "<li>Flenady 2014 Cochrane — tocolytics delay birth 48 h to 7 d with no ",
      "neonatal outcome benefit</li>",
      "<li>Roberts 2020 Cochrane — antenatal corticosteroids: RDS RR 0.66</li>",
      "<li>Rouse 2008 <i>NEJM</i> (BEAM) — MgSO4 neuroprotection, CP RR ~0.68, ",
      "NNT ~63</li>",
      "<li>Moise 1988 <i>NEJM</i> — indomethacin ductal constriction, steeply ",
      "GA-dependent</li>",
      "</ul>",

      "<h4>Known limitations</h4><ul>",
      "<li>Contractions are modelled as a smooth index, not as discrete ",
      "events; there is no uterine electromyography layer.</li>",
      "<li>Fast states (Ca<sup>2+</sup>, cAMP) are given large rate constants ",
      "rather than being solved on their true second-to-minute timescale — a ",
      "timescale-separation approximation on a day-scale model.</li>",
      "<li>Normalised indices (receptor densities, tissue integrity) are ",
      "relative-scale QSP states, not assayable concentrations.</li>",
      "<li>The neonatal outcome layer is a set of gestational-age logistics ",
      "with additive shifts, not a mechanistic neonatal model.</li>",
      "<li>The MEIS/PROLONG disagreement over 17-OHPC is <i>encoded</i>, not ",
      "resolved. Do not read the 17-OHPC scenario as evidence either way.</li>",
      "</ul>",

      "<hr><p style='color:#b71c1c;'><b>Disclaimer.</b> Research, education ",
      "and hypothesis generation only. This is not a tocolysis protocol, not ",
      "a dosing calculator, and has not been validated against patient data. ",
      "Threatened preterm labour is an obstetric emergency requiring clinical ",
      "judgment.</p>"
    ))
  })
}

shinyApp(ui, server)
