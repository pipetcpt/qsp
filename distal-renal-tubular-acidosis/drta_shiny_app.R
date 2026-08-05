## ===========================================================================
##  dRTA QSP EXPLORER — Shiny dashboard for drta_mrgsolve_model.R
##  10 tabs · run with:  shiny::runApp("drta_shiny_app.R")
## ===========================================================================
##
##  The app is organised around the model's four structural claims, so that a
##  user can try to BREAK each one rather than merely admire a curve:
##
##   Tab 2  "Acid balance & the 3 sinks"  — watch the acid gap partition into
##          ECF / buffer / bone and see that plasma HCO3- is a RATIO
##   Tab 4  "Urine chemistry"             — urine pH is SOLVED, so NH4+, TA and
##          bicarbonaturia all move together
##   Tab 6  "Schedule laboratory"         — matched daily mEq, different
##          schedules: the responder rate moves, the mean barely does
##   Tab 7  "Opposing kinetics"           — the fast-citrate / slow-bicarbonate
##          optimum, computed rather than asserted
##
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

MODEL_FILE <- "drta_mrgsolve_model.R"
mod <- mread_cache("drta", MODEL_FILE)
E <- as.list(mrgsolve::env(mod))          # SUBJ, LESION, SCENARIOS, helpers ...

`%||%` <- function(a, b) if (is.null(a)) b else a

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f5", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## reference bands used across tabs -----------------------------------------
BANDS <- list(
  HCO3_e   = c(22, 27,  "plasma HCO3- normal range"),
  K_pl     = c(3.5, 5.1, "plasma K+ normal range"),
  UpH_s    = c(5.0, 6.5, "urine pH, ad libitum"),
  UCit_s   = c(1.7, 6.5, "urine citrate normal (adult)"),
  SS_s     = c(0, 1,     "below the brushite threshold"),
  BMDz     = c(-1, 1,    "BMD z-score within 1 SD"),
  Hz       = c(-1, 1,    "height z-score within 1 SD")
)

band_layer <- function(v) {
  b <- BANDS[[v]]
  if (is.null(b)) return(NULL)
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = as.numeric(b[1]), ymax = as.numeric(b[2]),
           fill = "#2e7d32", alpha = 0.07)
}

## =========================================================================
##  UI
## =========================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well { background:#f7f9fa; border-color:#dde4e8; }
    h4 { color:#1a3d5c; font-weight:600; margin-top:2px; }
    .claim { background:#fff8e1; border-left:4px solid #f9a825;
             padding:8px 12px; margin:6px 0 12px 0; font-size:13px; }
    .metric { display:inline-block; min-width:132px; margin:3px 10px 3px 0;
              padding:6px 10px; background:#eef2f5; border-radius:6px; }
    .metric b { display:block; font-size:17px; color:#0d47a1; }
    .bad b  { color:#c62828; } .good b { color:#2e7d32; }
  "))),
  titlePanel("Distal Renal Tubular Acidosis — QSP Explorer"),
  p(tags$i(paste("A saturating acid-excretion actuator, three acid sinks on",
                 "three timescales, and alkali pharmacology in which the",
                 "delivery RATE is the active variable."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("subj", "Body size / age", names(E$SUBJ), selected = "child"),
      selectInput("lesion", "Lesion severity", names(E$LESION),
                  selected = "complete"),
      helpText(tags$small(paste("The lesion enters as exactly two numbers:",
                                "LES = retained pump Vmax (rate defect),",
                                "LES_grad = retained maximal blood-to-urine",
                                "pH gradient (gradient defect)."))),
      sliderInput("les", "LES — retained pump Vmax", 0.05, 1.0, 0.18, 0.01),
      sliderInput("lesg", "LES_grad — retained pH gradient", 0.40, 1.0, 0.58, 0.01),
      selectInput("geno", "Genotype preset (extrarenal only)",
                  c("(none)", names(E$GENOTYPE))),

      h4("Diet"),
      sliderInput("diet", "Dietary acid load (x normal NEAP)", 0.5, 2.6, 1.0, 0.05),

      h4("Alkali regimen"),
      selectInput("kind", "Formulation",
                  c("none", "bicIR", "citIR", "mixed", "ADV"), selected = "ADV"),
      selectInput("cation", "Salt cation", c("K", "NaK", "Na"), selected = "K"),
      sliderInput("dose", "Daily dose (mEq/kg/day)", 0, 4, 1.0, 0.05),
      selectInput("sched", "Intake schedule",
                  c("BID 08:00/20:00"      = "8,20",
                    "TID 07/13/19"         = "7,13,19",
                    "QID 07/12/18/23"      = "7,12,18,23",
                    "OD 08:00"             = "8",
                    "BID split 07/19"      = "7,19"), selected = "8,20"),
      sliderInput("fcit", "Fraction of alkali given as CITRATE", 0, 1, 0.35, 0.05),
      sliderInput("krbic", "Bicarbonate release rate (/h)", 0.05, 1.6, 0.215, 0.005),
      sliderInput("krcit", "Citrate release rate (/h)", 0.05, 1.6, 0.90, 0.005),

      h4("Adjuncts"),
      sliderInput("hctz", "Thiazide (mg/day)", 0, 50, 0, 5),
      sliderInput("kcl", "KCl supplement (mmol/day)", 0, 80, 0, 5),
      sliderInput("vitd", "Cholecalciferol (IU/day)", 0, 4000, 0, 200),

      h4("Adherence"),
      sliderInput("adh", "Adherence ceiling", 0.3, 1.0, 0.92, 0.02),

      h4("Simulation"),
      sliderInput("days", "Duration (days)", 30, 2190, 180, 30),
      actionButton("go", "Simulate", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 --------------------------------------------------------
        tabPanel(
          "1 · Patient profile",
          div(class = "claim", HTML(paste(
            "<b>The lesion is two numbers.</b> Everything else in this",
            "dashboard is a consequence. Note the <i>acid-excretion",
            "reserve</i>: in health it is 4-5x the daily acid load; when it",
            "reaches 1x the controller is railed and the patient has no",
            "capacity left for a dietary challenge or an intercurrent",
            "illness."))),
          uiOutput("kpi"),
          h4("Where the day's acid goes"),
          plotOutput("sankeyish", height = "300px"),
          h4("Derived subject parameters"),
          tableOutput("subjtab")
        ),

        ## ---- 2 --------------------------------------------------------
        tabPanel(
          "2 · Acid balance & the 3 sinks",
          div(class = "claim", HTML(paste(
            "<b>Claim 1 — plasma bicarbonate is a ratio, not a flux.</b>",
            "The acid gap is disposed of into ECF bicarbonate (minutes),",
            "non-bicarbonate buffer (hours) and BONE CARBONATE (years).",
            "Only the first is measured. Watch the bone flux stay positive",
            "while plasma HCO3- looks acceptable."))),
          plotOutput("pl_acid", height = "620px"),
          h4("Cumulative bone base withdrawn — the unmeasured integral"),
          plotOutput("pl_cum", height = "250px")
        ),

        ## ---- 3 --------------------------------------------------------
        tabPanel(
          "3 · Alkali PK",
          div(class = "claim", HTML(paste(
            "<b>Claim 3 — efficiency is delivery-rate matching.</b>",
            "The dashed line is the endogenous acid production rate. Base",
            "arriving <i>above</i> that line cannot be used and crosses the",
            "proximal reabsorptive threshold, leaving in the urine."))),
          plotOutput("pl_pk", height = "560px"),
          h4("Last three days, hourly resolution"),
          plotOutput("pl_pk_zoom", height = "300px")
        ),

        ## ---- 4 --------------------------------------------------------
        tabPanel(
          "4 · Urine chemistry",
          div(class = "claim", HTML(paste(
            "<b>Urine pH is solved, not fitted.</b> It is the root of the",
            "luminal proton balance, so NH4+ trapping, titratable acid and",
            "bicarbonaturia are forced to move consistently. A positive urine",
            "anion gap in dRTA is <i>low ammonium</i>, and the model produces",
            "it rather than asserting it."))),
          plotOutput("pl_urine", height = "620px"),
          h4("Net acid excretion, decomposed"),
          plotOutput("pl_nae", height = "270px")
        ),

        ## ---- 5 --------------------------------------------------------
        tabPanel(
          "5 · Stone & nephrocalcinosis risk",
          div(class = "claim", HTML(paste(
            "<b>dRTA is the one acidosis whose urine is alkaline</b>, which is",
            "why it makes calcium-PHOSPHATE stones. Alkali therapy pushes two",
            "arms in opposite directions: it raises urine citrate (protective)",
            "but also raises urine pH (raising HPO4(2-), harmful). In COMPLETE",
            "dRTA the pH arm is already near saturation, so alkali wins",
            "outright; in INCOMPLETE dRTA it has headroom and the benefit is",
            "not automatic."))),
          plotOutput("pl_stone", height = "620px"),
          fluidRow(
            column(6, h4("Supersaturation decomposition"),
                   plotOutput("pl_ssdec", height = "290px")),
            column(6, h4("Urine pH sensitivity of brushite SS"),
                   plotOutput("pl_phsens", height = "290px")))
        ),

        ## ---- 6 --------------------------------------------------------
        tabPanel(
          "6 · Schedule laboratory",
          div(class = "claim", HTML(paste(
            "<b>Same daily mEq/kg, different schedule.</b> This is the B21CS",
            "experiment. Because the alkali waste is worst <i>at</i> the",
            "therapeutic target, and because the bone term is rectified and",
            "therefore convex, a schedule change moves the responder rate and",
            "the bone flux far more than it moves the mean bicarbonate."))),
          fluidRow(
            column(4, sliderInput("sl_dose", "Matched daily dose (mEq/kg/day)",
                                  0.2, 3.0, 1.0, 0.05)),
            column(4, sliderInput("sl_days", "Duration (days)", 60, 730, 180, 30)),
            column(4, br(), actionButton("go_sched", "Compare schedules",
                                         class = "btn-primary"))),
          tableOutput("sched_tab"),
          plotOutput("pl_sched", height = "420px"),
          h4("24-hour bicarbonate profiles at matched daily dose"),
          plotOutput("pl_sched_prof", height = "320px")
        ),

        ## ---- 7 --------------------------------------------------------
        tabPanel(
          "7 · Opposing kinetics",
          div(class = "claim", HTML(paste(
            "<b>Claim 4 — the two endpoints want opposite kinetics.</b>",
            "Systemic alkalinisation wants SLOW delivery (rate matching).",
            "Citraturia wants FAST delivery, because NaDC1 is Tm-limited and a",
            "bolus escapes reabsorption. Sweep the two release rates and see",
            "that the joint optimum is a fast-citrate / slow-bicarbonate",
            "combination — which is what a two-granule formulation is."))),
          fluidRow(
            column(4, sliderInput("op_dose", "Daily dose (mEq/kg/day)",
                                  0.3, 3.0, 1.0, 0.05)),
            column(4, sliderInput("op_n", "Release rates per axis", 3, 7, 5, 1)),
            column(4, br(), actionButton("go_opp", "Run the sweep",
                                         class = "btn-primary"))),
          fluidRow(
            column(6, plotOutput("pl_opp_hco3", height = "340px")),
            column(6, plotOutput("pl_opp_cit", height = "340px"))),
          h4("Joint objective (bicarbonate control x citraturia)"),
          plotOutput("pl_opp_joint", height = "340px"),
          tableOutput("opp_tab")
        ),

        ## ---- 8 --------------------------------------------------------
        tabPanel(
          "8 · Bone, growth & CKD",
          div(class = "claim", HTML(paste(
            "<b>The dissociation that names the disease.</b> In the B22CS",
            "cohort plasma HCO3- was 22.0 mmol/L (near normal) while lumbar",
            "BMD z was -1.1 (clearly abnormal) at the same visit. The model",
            "produces both from one lesion because the two readouts have",
            "different gearing on the same acid gap."))),
          plotOutput("pl_bone", height = "620px"),
          h4("Long-run trajectories against the B22CS 6-year anchors"),
          tableOutput("b22_tab")
        ),

        ## ---- 9 --------------------------------------------------------
        tabPanel(
          "9 · Diagnostic tests",
          div(class = "claim", HTML(paste(
            "<b>Incomplete dRTA is derived, not assumed.</b> Give an acute",
            "NH4Cl load and watch the urine pH nadir. A normal subject reaches",
            "below 5.45; a patient whose actuator is already railed cannot,",
            "even with a completely normal resting bicarbonate."))),
          fluidRow(
            column(4, sliderInput("nh4", "NH4Cl load (mEq/kg)", 0.5, 3.0, 1.87, 0.1)),
            column(4, sliderInput("dx_days", "Baseline before load (days)",
                                  20, 120, 60, 10)),
            column(4, br(), actionButton("go_dx", "Run the loading test",
                                         class = "btn-primary"))),
          plotOutput("pl_dx", height = "440px"),
          tableOutput("dx_tab"),
          h4("Static discriminators at the current steady state"),
          tableOutput("disc_tab")
        ),

        ## ---- 10 -------------------------------------------------------
        tabPanel(
          "10 · Scenario library",
          div(class = "claim", HTML(paste(
            "All 28 predefined scenarios from the mrgsolve model, with the",
            "note that states what each one is FOR and what it must",
            "reproduce."))),
          selectInput("sc", "Scenario", names(E$SCENARIOS), width = "60%"),
          verbatimTextOutput("sc_note"),
          actionButton("go_sc", "Run scenario", class = "btn-primary"),
          br(), br(),
          plotOutput("pl_sc", height = "560px"),
          tableOutput("sc_tab")
        )
      )
    )
  )
)

## =========================================================================
##  SERVER
## =========================================================================
server <- function(input, output, session) {

  ## keep the lesion sliders in step with the preset ----------------------
  observeEvent(input$lesion, {
    L <- E$LESION[[input$lesion]]
    updateSliderInput(session, "les",  value = L$LES)
    updateSliderInput(session, "lesg", value = L$LES_grad)
  })
  observeEvent(input$geno, {
    if (input$geno != "(none)") {
      G <- E$GENOTYPE[[input$geno]]
      updateSliderInput(session, "les",  value = G$LES)
      updateSliderInput(session, "lesg", value = G$LES_grad)
    }
  })
  observeEvent(input$kind, {
    fc <- switch(input$kind, bicIR = 0, citIR = 1, mixed = 0.5, ADV = 0.35, 0)
    updateSliderInput(session, "fcit", value = fc)
    if (input$kind == "ADV") {
      updateSliderInput(session, "krbic", value = 0.215)
      updateSliderInput(session, "krcit", value = 0.90)
      updateSelectInput(session, "sched", selected = "8,20")
    } else {
      updateSliderInput(session, "krbic", value = 1.50)
      updateSliderInput(session, "krcit", value = 1.05)
    }
  })

  times_of <- reactive(as.numeric(strsplit(input$sched, ",")[[1]]))

  pars <- reactive({
    s <- E$SUBJ[[input$subj]]
    p <- c(s, list(LES = input$les, LES_grad = input$lesg, DIET = input$diet,
                   ADH_ref = input$adh, NINTAKE = length(times_of()),
                   f_Kcat = E$f_Kcat_of(input$cation),
                   f_cit_ADV = input$fcit,
                   kr_bicPR = input$krbic, kr_citPR = input$krcit))
    if (input$geno != "(none)")
      p$kHEAR <- E$GENOTYPE[[input$geno]]$kHEAR
    p
  })

  ## ---- generic runner --------------------------------------------------
  run_one <- function(p, kind, dose, times, cation, days, extra = list(),
                      acid = NULL, init = NULL, delta = 0.5, fcit = NULL) {
    p[names(extra)] <- extra
    BW <- p$BW
    ev <- E$make_regimen(kind = kind, mEq_kg_day = dose, times = times,
                         cation = cation, days = days, BW = BW,
                         f_cit = fcit,
                         KCl_mmol_day = input$kcl %||% 0,
                         hctz_mg = input$hctz %||% 0,
                         vitD_IU = input$vitd %||% 0)
    if (!is.null(acid))
      ev <- rbind(ev, data.frame(time = acid[1], cmt = E$CMT["AG_acid"],
                                 amt = acid[2] * BW, evid = 1))
    ev <- ev[order(ev$time), ]
    m <- param(mod, p)
    if (!is.null(init)) m <- init(m, init)
    as.data.frame(mrgsim(m, data = ev, end = days * 24, delta = delta,
                         maxsteps = 5e6, hmax = 0.5))
  }

  sim_main <- eventReactive(input$go, {
    run_one(pars(), input$kind, input$dose, times_of(), input$cation,
            input$days, fcit = input$fcit)
  }, ignoreNULL = FALSE)

  last_day <- reactive({ d <- sim_main(); d[d$time >= max(d$time) - 24, ] })

  ## ---- tab 1: KPIs -----------------------------------------------------
  output$kpi <- renderUI({
    d <- last_day()
    box <- function(lab, val, cls = "") {
      div(class = paste("metric", cls), tags$span(lab), tags$b(val))
    }
    ok <- function(x, lo, hi) if (x >= lo && x <= hi) "good" else "bad"
    h <- mean(d$HCO3_e); k <- mean(d$K_pl); u <- mean(d$urine_pH)
    tagList(
      box("plasma HCO3-", sprintf("%.1f mmol/L", h), ok(h, 22, 27)),
      box("blood pH", sprintf("%.3f", mean(d$pH_blood)), ok(mean(d$pH_blood), 7.35, 7.45)),
      box("urine pH", sprintf("%.2f", u), if (u < 5.5) "good" else "bad"),
      box("plasma K+", sprintf("%.2f mmol/L", k), ok(k, 3.5, 5.1)),
      box("plasma Cl-", sprintf("%.0f mmol/L", mean(d$Cl_pl))),
      box("NAE", sprintf("%.0f mEq/day", mean(d$NAE_day))),
      box("NEAP", sprintf("%.0f mEq/day", mean(d$NEAP_day))),
      box("unmet acid gap", sprintf("%.1f mEq/day", mean(d$acid_gap)),
          if (mean(d$acid_gap) < 1) "good" else "bad"),
      box("actuator VH", sprintf("%.2f", mean(d$VH_eff)),
          if (mean(d$VH_eff) >= 0.999) "bad" else "good"),
      box("excretion reserve", sprintf("%.2f x", mean(d$reserve)),
          if (mean(d$reserve) > 2) "good" else "bad"),
      box("bone base flux", sprintf("%.1f mEq/day", mean(d$bone_day)),
          if (mean(d$bone_day) < 1) "good" else "bad"),
      box("urine Ca", sprintf("%.1f mg/kg/d", mean(d$UCa_mgkg)),
          if (mean(d$UCa_mgkg) < 4) "good" else "bad"),
      box("urine citrate", sprintf("%.2f mmol/d", mean(d$UCit_day))),
      box("Ca/citrate", sprintf("%.2f", mean(d$CaCit)),
          if (mean(d$CaCit) < 0.33) "good" else "bad"),
      box("brushite SS", sprintf("%.2f", mean(d$SS_inst)),
          if (mean(d$SS_inst) < 1) "good" else "bad"),
      box("FE HCO3-", sprintf("%.2f %%", mean(d$FE_HCO3))),
      box("alkali wasted", sprintf("%.0f %%", 100 * tail(d$waste_frac, 1))),
      box("adherence", sprintf("%.2f", tail(d$ADH, 1))),
      box("BMD z", sprintf("%.2f", tail(d$BMDz, 1)), ok(tail(d$BMDz, 1), -1, 3)),
      box("height z", sprintf("%.2f", tail(d$Hz, 1)), ok(tail(d$Hz, 1), -1, 3)),
      box("eGFR", sprintf("%.0f", tail(d$eGFR, 1)))
    )
  })

  output$sankeyish <- renderPlot({
    d <- last_day()
    v <- c(`NEAP (acid in)` = mean(d$NEAP_day),
           `alkali absorbed` = mean(d$alkali_day),
           `renal NAE` = mean(d$NAE_day),
           `bone base` = mean(d$bone_day),
           `urinary HCO3- (wasted)` = mean(d$HCO3u_day),
           `unmet gap` = mean(d$acid_gap))
    df <- data.frame(term = factor(names(v), levels = names(v)), value = as.numeric(v))
    ggplot(df, aes(term, value, fill = term)) +
      geom_col(width = .68, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.4, size = 4) +
      labs(x = NULL, y = "mEq/day") +
      scale_fill_manual(values = c("#ef6c00", "#2e7d32", "#1565c0",
                                   "#6a1b9a", "#c62828", "#f9a825")) +
      THEME + theme(axis.text.x = element_text(angle = 18, hjust = 1))
  })

  output$subjtab <- renderTable({
    d <- last_day(); p <- pars()
    data.frame(
      quantity = c("HCO3- setpoint (age)", "proximal HCO3- threshold",
                   "NEAP", "max distal H+ secretion", "acid-excretion reserve",
                   "minimum achievable urine pH", "apparent HCO3- space"),
      value = c(sprintf("%.2f mmol/L", tail(d$HCO3_target, 1)),
                sprintf("%.2f mmol/L", tail(d$HCO3_target, 1) + 1.05),
                sprintf("%.0f mEq/day", mean(d$NEAP_day)),
                sprintf("%.1f mEq/day", 11.0 * p$BSA / 1.73 * p$LES * 24),
                sprintf("%.2f x", mean(d$reserve)),
                sprintf("%.2f", mean(d$pH_blood) - 2.95 * p$LES_grad),
                sprintf("%.2f L/kg", 0.20 + 0.30)))
  }, striped = TRUE)

  ## ---- long-format plot helper ----------------------------------------
  tsplot <- function(d, vars, labs_ = NULL, xdays = TRUE, ncol = 2) {
    dd <- d[, c("time", vars), drop = FALSE]
    long <- pivot_longer(dd, -time, names_to = "v", values_to = "y")
    if (!is.null(labs_)) long$v <- factor(long$v, levels = vars, labels = labs_)
    else long$v <- factor(long$v, levels = vars)
    long$x <- if (xdays) long$time / 24 else long$time
    p <- ggplot(long, aes(x, y)) +
      geom_line(colour = "#0d47a1", linewidth = .55) +
      facet_wrap(~v, scales = "free_y", ncol = ncol) +
      labs(x = if (xdays) "days" else "hours", y = NULL) + THEME
    p
  }

  ## ---- tab 2 -----------------------------------------------------------
  output$pl_acid <- renderPlot({
    d <- sim_main()
    v <- c("HCO3_e", "pH_blood", "PaCO2", "acid_gap", "bone_day", "VH_eff",
           "BUF", "Cl_pl")
    l <- c("plasma HCO3- (mmol/L)  [SINK 1]", "blood pH", "PaCO2 (mmHg)",
           "unmet acid gap (mEq/day)", "bone base flux (mEq/day)  [SINK 3]",
           "actuator VH (1.0 = RAILED)", "buffer base donated (mEq)  [SINK 2]",
           "plasma Cl- (mmol/L)")
    tsplot(d, v, l)
  })
  output$pl_cum <- renderPlot({
    d <- sim_main()
    ggplot(d, aes(time / 24, CUMBASE / 1000)) +
      geom_area(fill = "#6a1b9a", alpha = .25) +
      geom_line(colour = "#6a1b9a", linewidth = .7) +
      labs(x = "days", y = "cumulative bone base withdrawn (Eq)") + THEME
  })

  ## ---- tab 3 -----------------------------------------------------------
  output$pl_pk <- renderPlot({
    d <- sim_main()
    v <- c("AG_bicPR", "AG_citPR", "AG_bicIR", "CIT_pl", "alkali_day",
           "HCO3u_day", "waste_frac", "GI")
    l <- c("slow granules remaining (mEq)", "fast granules remaining (mmol)",
           "IR alkali in gut (mEq)", "plasma citrate (mmol/L)",
           "alkali absorbed (mEq/day)", "urinary HCO3- (mEq/day)",
           "cumulative wasted fraction", "GI irritation index")
    tsplot(d, v, l)
  })
  output$pl_pk_zoom <- renderPlot({
    d <- sim_main(); d <- d[d$time >= max(d$time) - 72, ]
    dd <- data.frame(t = (d$time - min(d$time)),
                     alkali = d$alkali_day / 24, NEAP = d$NEAP_day / 24)
    long <- pivot_longer(dd, -t)
    ggplot(long, aes(t, value, colour = name, linetype = name)) +
      geom_line(linewidth = .75) +
      scale_colour_manual(values = c(alkali = "#2e7d32", NEAP = "#ef6c00"),
                          name = NULL) +
      scale_linetype_manual(values = c(alkali = 1, NEAP = 2), name = NULL) +
      labs(x = "hours (last 3 days)", y = "mEq/h",
           subtitle = paste("Base above the acid-production line cannot be",
                            "used and is voided as bicarbonaturia")) + THEME
  })

  ## ---- tab 4 -----------------------------------------------------------
  output$pl_urine <- renderPlot({
    d <- sim_main()
    v <- c("urine_pH", "NH4_day", "TA_day", "HCO3u_day",
           "UVol_s", "Pi_urine", "FE_HCO3", "UCit_day")
    l <- c("urine pH (SOLVED)", "urinary NH4+ (mEq/day)",
           "titratable acid (mEq/day)", "urinary HCO3- (mEq/day)",
           "urine volume (L/day)", "urinary phosphate (mmol/day)",
           "FE HCO3- (%)  <3% = distal", "urine citrate (mmol/day)")
    tsplot(d, v, l)
  })
  output$pl_nae <- renderPlot({
    d <- sim_main(); d <- d[seq(1, nrow(d), length.out = min(600, nrow(d))), ]
    dd <- data.frame(t = d$time / 24, NH4 = d$NH4_day, TA = d$TA_day,
                     `minus HCO3` = -d$HCO3u_day, check.names = FALSE)
    long <- pivot_longer(dd, -t)
    ggplot(long, aes(t, value, fill = name)) +
      geom_area(alpha = .8) +
      geom_line(data = d, aes(time / 24, NAE_day), inherit.aes = FALSE,
                colour = "black", linewidth = .6) +
      scale_fill_manual(values = c(NH4 = "#2e7d32", TA = "#f9a825",
                                   `minus HCO3` = "#c62828"), name = NULL) +
      labs(x = "days", y = "mEq/day",
           subtitle = "black line = net acid excretion") + THEME
  })

  ## ---- tab 5 -----------------------------------------------------------
  output$pl_stone <- renderPlot({
    d <- sim_main()
    v <- c("SS_inst", "UCa_mgkg", "UCit_day", "CaCit",
           "NC", "STONE", "eGFR", "FIB")
    l <- c("brushite supersaturation (1 = threshold)",
           "urine Ca (mg/kg/day; >4 = hypercalciuria)",
           "urine citrate (mmol/day)", "urine Ca/citrate (>0.33 lithogenic)",
           "nephrocalcinosis burden", "stone burden",
           "eGFR (mL/min/1.73m2)", "interstitial fibrosis")
    tsplot(d, v, l)
  })
  output$pl_ssdec <- renderPlot({
    d <- last_day()
    v <- c(`free Ca2+ (mmol/L)` = mean(d$UCa_day) / mean(d$UVol_s),
           `HPO4(2-) (mmol/L)` = mean(d$Pi_urine) / mean(d$UVol_s) /
             (1 + 10^(6.8 - mean(d$urine_pH))),
           `citrate (mmol/L)` = mean(d$UCit_day) / mean(d$UVol_s))
    df <- data.frame(term = names(v), value = as.numeric(v))
    ggplot(df, aes(term, value, fill = term)) +
      geom_col(show.legend = FALSE, width = .6) +
      geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.4) +
      scale_fill_manual(values = c("#c2185b", "#5d4037", "#558b2f")) +
      labs(x = NULL, y = "mmol/L") + THEME
  })
  output$pl_phsens <- renderPlot({
    d <- last_day()
    pi_c <- mean(d$Pi_urine) / mean(d$UVol_s)
    ca_c <- mean(d$UCa_day) / mean(d$UVol_s)
    cit_c <- mean(d$UCit_day) / mean(d$UVol_s)
    ph <- seq(5.0, 8.0, 0.02)
    ss <- (ca_c / (1 + 0.28 * cit_c)) * (pi_c / (1 + 10^(6.8 - ph))) / 4.25
    df <- data.frame(ph, ss)
    ggplot(df, aes(ph, ss)) + geom_line(colour = "#5d4037", linewidth = .8) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#c62828") +
      geom_vline(xintercept = mean(d$urine_pH), linetype = 3) +
      annotate("text", x = mean(d$urine_pH), y = Inf, vjust = 1.4, hjust = -0.1,
               label = "current urine pH", size = 3.4) +
      labs(x = "urine pH", y = "brushite SS",
           subtitle = paste("How much headroom does urine pH still have?",
                            "Steep here = alkali carries a real cost")) + THEME
  })

  ## ---- tab 6: schedule laboratory -------------------------------------
  sched_res <- eventReactive(input$go_sched, {
    p <- pars(); dose <- input$sl_dose; days <- input$sl_days
    grid <- list(
      list(lab = "KHCO3 IR, 3x/day", kind = "bicIR", t = c(7, 13, 19), cat = "K", fc = 0),
      list(lab = "NaHCO3 IR, 3x/day", kind = "bicIR", t = c(7, 13, 19), cat = "Na", fc = 0),
      list(lab = "K-citrate IR, 3x/day", kind = "citIR", t = c(7, 13, 19), cat = "K", fc = 1),
      list(lab = "KHCO3 IR, 2x/day", kind = "bicIR", t = c(8, 20), cat = "K", fc = 0),
      list(lab = "KHCO3 IR, 4x/day + night", kind = "bicIR", t = c(7, 12, 18, 23), cat = "K", fc = 0),
      list(lab = "ADV7103, 2x/day", kind = "ADV", t = c(8, 20), cat = "K", fc = 0.35),
      list(lab = "ADV7103, 1x/day", kind = "ADV", t = c(8), cat = "K", fc = 0.35))
    do.call(rbind, lapply(grid, function(g) {
      pp <- p
      pp$NINTAKE <- length(g$t)
      pp$f_Kcat <- E$f_Kcat_of(g$cat)
      pp$f_cit_ADV <- g$fc
      if (g$kind == "ADV") { pp$kr_bicPR <- 0.215; pp$kr_citPR <- 0.90 }
      else { pp$kr_bicPR <- 1.50; pp$kr_citPR <- 1.05 }
      d <- run_one(pp, g$kind, dose, g$t, g$cat, days, fcit = g$fc)
      ld <- d[d$time >= max(d$time) - 24, ]
      data.frame(regimen = g$lab, intakes = length(g$t),
                 HCO3 = mean(ld$HCO3_e), trough = min(ld$HCO3_e),
                 swing = max(ld$HCO3_e) - min(ld$HCO3_e),
                 responder = mean(ld$responder) >= .5,
                 urine_pH = mean(ld$urine_pH),
                 wasted_pct = 100 * tail(ld$waste_frac, 1),
                 UCit = mean(ld$UCit_day), CaCit = mean(ld$CaCit),
                 SS = mean(ld$SS_inst), bone = mean(ld$bone_day),
                 K = mean(ld$K_pl), adherence = tail(ld$ADH, 1),
                 BMDz = tail(ld$BMDz, 1))
    }))
  })
  output$sched_tab <- renderTable(sched_res(), digits = 2, striped = TRUE)
  output$pl_sched <- renderPlot({
    r <- sched_res()
    long <- pivot_longer(r[, c("regimen", "HCO3", "trough", "bone", "UCit",
                               "SS", "wasted_pct")], -regimen)
    long$name <- factor(long$name,
      levels = c("HCO3", "trough", "wasted_pct", "bone", "UCit", "SS"),
      labels = c("mean HCO3-", "trough HCO3-", "alkali wasted (%)",
                 "bone base flux (mEq/d)", "urine citrate (mmol/d)",
                 "brushite SS"))
    ggplot(long, aes(reorder(regimen, value), value, fill = regimen)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x", ncol = 3) +
      labs(x = NULL, y = NULL) + THEME
  })
  output$pl_sched_prof <- renderPlot({
    p <- pars(); dose <- input$sl_dose
    grid <- list(list(lab = "IR 3x/day", kind = "bicIR", t = c(7, 13, 19), fc = 0),
                 list(lab = "IR 4x/day", kind = "bicIR", t = c(7, 12, 18, 23), fc = 0),
                 list(lab = "ADV7103 2x/day", kind = "ADV", t = c(8, 20), fc = .35))
    dd <- do.call(rbind, lapply(grid, function(g) {
      pp <- p; pp$NINTAKE <- length(g$t); pp$f_cit_ADV <- g$fc
      if (g$kind == "ADV") { pp$kr_bicPR <- .215; pp$kr_citPR <- .90 }
      else { pp$kr_bicPR <- 1.5; pp$kr_citPR <- 1.05 }
      d <- run_one(pp, g$kind, dose, g$t, "K", 90, delta = 0.25, fcit = g$fc)
      d <- d[d$time >= max(d$time) - 48, ]
      data.frame(regimen = g$lab, h = d$time - min(d$time), HCO3 = d$HCO3_e)
    }))
    ggplot(dd, aes(h, HCO3, colour = regimen)) + geom_line(linewidth = .75) +
      geom_hline(yintercept = 22, linetype = 2, colour = "#c62828") +
      scale_colour_manual(values = c("#c62828", "#f9a825", "#1565c0"), name = NULL) +
      labs(x = "hours (last 2 days)", y = "plasma HCO3- (mmol/L)",
           subtitle = "dashed = responder threshold") + THEME
  })

  ## ---- tab 7: opposing kinetics ---------------------------------------
  opp_res <- eventReactive(input$go_opp, {
    p <- pars(); n <- input$op_n
    rates <- exp(seq(log(0.10), log(1.50), length.out = n))
    grid <- expand.grid(kr_bic = rates, kr_cit = rates)
    res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      pp <- p
      pp$kr_bicPR <- grid$kr_bic[i]; pp$kr_citPR <- grid$kr_cit[i]
      pp$f_cit_ADV <- 0.35; pp$NINTAKE <- 2
      d <- run_one(pp, "ADV", input$op_dose, c(8, 20), "K", 120, fcit = 0.35)
      ld <- d[d$time >= max(d$time) - 24, ]
      data.frame(kr_bic = grid$kr_bic[i], kr_cit = grid$kr_cit[i],
                 HCO3 = mean(ld$HCO3_e), UCit = mean(ld$UCit_day),
                 SS = mean(ld$SS_inst), bone = mean(ld$bone_day),
                 wasted = 100 * tail(ld$waste_frac, 1))
    }))
    res$joint <- scale(res$HCO3)[, 1] + scale(res$UCit)[, 1]
    res
  })
  heat <- function(r, z, ttl, sub) {
    ggplot(r, aes(kr_bic, kr_cit, fill = .data[[z]])) +
      geom_raster(interpolate = TRUE) +
      scale_x_log10() + scale_y_log10() +
      scale_fill_viridis_c(option = "D", name = NULL) +
      labs(x = "bicarbonate release rate (/h)  <- slower | faster ->",
           y = "citrate release rate (/h)", title = ttl, subtitle = sub) + THEME
  }
  output$pl_opp_hco3 <- renderPlot(
    heat(opp_res(), "HCO3", "Mean plasma HCO3-",
         "best towards the LEFT: slow bicarbonate release"))
  output$pl_opp_cit <- renderPlot(
    heat(opp_res(), "UCit", "Urine citrate",
         "best towards the TOP: fast citrate release"))
  output$pl_opp_joint <- renderPlot({
    r <- opp_res(); b <- r[which.max(r$joint), ]
    heat(r, "joint", "Joint objective (z(HCO3-) + z(urine citrate))",
         sprintf(paste("optimum at bicarbonate %.2f /h, citrate %.2f /h",
                       "— i.e. SLOW bicarbonate with FAST citrate"),
                 b$kr_bic, b$kr_cit)) +
      annotate("point", x = b$kr_bic, y = b$kr_cit, colour = "red",
               size = 4, shape = 4, stroke = 2)
  })
  output$opp_tab <- renderTable({
    r <- opp_res()
    rbind(
      cbind(objective = "best HCO3-", r[which.max(r$HCO3), ]),
      cbind(objective = "best urine citrate", r[which.max(r$UCit), ]),
      cbind(objective = "best joint", r[which.max(r$joint), ]),
      cbind(objective = "lowest brushite SS", r[which.min(r$SS), ]))
  }, digits = 3, striped = TRUE)

  ## ---- tab 8 -----------------------------------------------------------
  output$pl_bone <- renderPlot({
    d <- sim_main()
    v <- c("BMDz", "Hz", "BMIN", "bone_day", "OC", "bALP", "OSM", "IGF1")
    l <- c("lumbar BMD z-score", "height z-score",
           "bone mineral (fraction expected)", "bone base flux (mEq/day)",
           "osteoclast activity", "bone ALP (U/L)",
           "osteomalacia index", "IGF-1 (ng/mL)")
    tsplot(d, v, l)
  })
  output$b22_tab <- renderTable({
    d <- sim_main(); ld <- last_day()
    data.frame(
      endpoint = c("plasma HCO3- (mmol/L)", "height z-score",
                   "lumbar BMD z-score", "eGFR (mL/min/1.73m2)"),
      model_start = c(sprintf("%.1f", d$HCO3_e[1]), sprintf("%.2f", d$Hz[1]),
                      sprintf("%.2f", d$BMDz[1]), sprintf("%.0f", d$eGFR[1])),
      model_end = c(sprintf("%.1f", mean(ld$HCO3_e)),
                    sprintf("%.2f", tail(ld$Hz, 1)),
                    sprintf("%.2f", tail(ld$BMDz, 1)),
                    sprintf("%.0f", tail(ld$eGFR, 1))),
      B22CS_6yr = c("22.0 -> 22.6 (NS)", "-0.6 -> -0.3 (p=0.03)",
                    "-1.1 -> -0.8 (p=0.005)", "105 -> 104 (NS)"))
  }, striped = TRUE)

  ## ---- tab 9 -----------------------------------------------------------
  dx_res <- eventReactive(input$go_dx, {
    p <- pars(); d0 <- input$dx_days
    do.call(rbind, lapply(c("none", "incomplete", "complete", "gradonly"),
                          function(ln) {
      L <- E$LESION[[ln]]
      pp <- p; pp$LES <- L$LES; pp$LES_grad <- L$LES_grad
      d <- run_one(pp, input$kind, if (ln == "none") 0 else input$dose,
                   times_of(), input$cation, d0 + 2,
                   acid = c(d0 * 24 + 8, input$nh4), delta = 0.1,
                   fcit = input$fcit)
      w <- d[d$time >= d0 * 24 + 6 & d$time <= d0 * 24 + 22, ]
      w$lesion <- ln; w$h <- w$time - d0 * 24
      w
    }))
  })
  output$pl_dx <- renderPlot({
    d <- dx_res()
    long <- pivot_longer(d[, c("h", "lesion", "urine_pH", "HCO3_e")],
                         c(urine_pH, HCO3_e))
    long$name <- factor(long$name, c("urine_pH", "HCO3_e"),
                        c("urine pH", "plasma HCO3- (mmol/L)"))
    ggplot(long, aes(h, value, colour = lesion)) +
      geom_line(linewidth = .8) +
      geom_hline(data = data.frame(name = factor("urine pH",
                    levels = levels(long$name)), y = 5.45),
                 aes(yintercept = y), linetype = 2, colour = "#c62828") +
      geom_vline(xintercept = 8, linetype = 3) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_brewer(palette = "Dark2", name = NULL) +
      labs(x = "hours on the test day (load at 08:00)", y = NULL,
           subtitle = "dashed = 5.45, the normal urine-pH nadir") + THEME
  })
  output$dx_tab <- renderTable({
    d <- dx_res()
    do.call(rbind, lapply(split(d, d$lesion), function(x) {
      w <- x[x$h >= 8 & x$h <= 18, ]
      data.frame(lesion = x$lesion[1],
                 pre_HCO3 = x$HCO3_e[which.min(abs(x$h - 7))],
                 pre_urine_pH = x$urine_pH[which.min(abs(x$h - 7))],
                 nadir_urine_pH = min(w$urine_pH),
                 acidifies = min(w$urine_pH) < 5.45)
    }))
  }, digits = 2, striped = TRUE)
  output$disc_tab <- renderTable({
    d <- last_day()
    data.frame(
      discriminator = c("urine pH during acidaemia", "FE HCO3-",
                        "urine anion gap direction", "serum anion gap",
                        "plasma K+"),
      model = c(sprintf("%.2f", mean(d$urine_pH)),
                sprintf("%.2f %%", mean(d$FE_HCO3)),
                if (mean(d$NH4_day) < 30) "positive (low NH4+)" else "negative",
                "normal (hyperchloraemic)",
                sprintf("%.2f mmol/L", mean(d$K_pl))),
      interpretation = c(">5.5 with acidaemia = distal",
                         "<3% distal, >15% proximal",
                         "positive = impaired NH4+ excretion",
                         "normal gap excludes ketoacidosis / lactate",
                         "low in classic dRTA, high in type 4"))
  }, striped = TRUE)

  ## ---- tab 10 ----------------------------------------------------------
  output$sc_note <- renderText({
    sc <- E$SCENARIOS[[input$sc]]
    paste0(input$sc, "\n\nsubject : ", sc$subj,
           "\nlesion  : ", sc$lesion,
           "\nduration: ", sc$days, " days\n\n",
           gsub("\\s+", " ", sc$note))
  })
  sc_res <- eventReactive(input$go_sc, {
    as.data.frame(E$sim_scenario(mod, input$sc))
  })
  output$pl_sc <- renderPlot({
    d <- sc_res()
    v <- c("HCO3_e", "urine_pH", "K_pl", "UCit_day", "UCa_mgkg", "SS_inst",
           "BMDz", "Hz", "eGFR")
    l <- c("plasma HCO3-", "urine pH", "plasma K+", "urine citrate (mmol/d)",
           "urine Ca (mg/kg/d)", "brushite SS", "BMD z", "height z", "eGFR")
    tsplot(d, v, l, ncol = 3)
  })
  output$sc_tab <- renderTable({
    d <- sc_res(); ld <- d[d$time >= max(d$time) - 24, ]
    data.frame(t(sapply(list(
      `HCO3- mean` = mean(ld$HCO3_e), `HCO3- min` = min(ld$HCO3_e),
      `blood pH` = mean(ld$pH_blood), `urine pH` = mean(ld$urine_pH),
      `NAE mEq/d` = mean(ld$NAE_day), `gap mEq/d` = mean(ld$acid_gap),
      `VH` = mean(ld$VH_eff), `reserve` = mean(ld$reserve),
      `K+` = mean(ld$K_pl), `urine Ca mg/kg/d` = mean(ld$UCa_mgkg),
      `urine citrate` = mean(ld$UCit_day), `Ca/citrate` = mean(ld$CaCit),
      `brushite SS` = mean(ld$SS_inst), `bone mEq/d` = mean(ld$bone_day),
      `BMD z` = tail(ld$BMDz, 1), `height z` = tail(ld$Hz, 1),
      `eGFR` = tail(ld$eGFR, 1)), function(x) round(x, 3))))
  }, striped = TRUE, rownames = FALSE)
}

shinyApp(ui, server)
