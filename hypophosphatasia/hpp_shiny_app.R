## ============================================================================
##  Hypophosphatasia (HPP) QSP Model — Shiny Dashboard
## ============================================================================
##  Interactive front-end for hpp_mrgsolve_model.R.
##
##  The app is organised around the four structural claims of the model rather
##  than around organ systems, because that is what makes it useful:
##
##    Tab 2  PK / target delivery ...... plasma exposure vs HA-BOUND enzyme
##    Tab 3  PPi economy ............... the causal variable you cannot measure
##    Tab 4  Biochemical markers ....... the variables you CAN measure, and how
##                                       early they saturate  (T2)
##    Tab 5  Skeletal endpoints ........ RSS, osteoid, growth plate, stature
##    Tab 6  Respiratory & survival .... what actually kills infants
##    Tab 7  Mineral & renal ........... paradoxical hypercalcaemia
##    Tab 8  CNS / vitamin B6 .......... substrate excess, product deficiency
##    Tab 9  Scenario comparison ....... 17 predefined arms
##    Tab 10 Threshold explorer ........ T1: the severity spectrum is ONE curve
##    Tab 11 Dose-ranging .............. T2: which readout saturates first?
##    Tab 12 Model notes ............... assumptions, calibration, limitations
##
##  Run:
##    setwd("hypophosphatasia"); shiny::runApp("hpp_shiny_app.R")
##
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr
##  Disclaimer: research / education only. NOT for clinical use.
## ============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("hpp_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

PHENO <- c("Perinatal (severe, ~1% activity)"        = "perinatal",
           "Infantile (~3% activity)"                = "infantile",
           "Childhood (~12% activity)"               = "childhood",
           "Adult / odonto (~35% activity)"          = "adult",
           "Heterozygous carrier (~50% activity)"    = "carrier")

## helper: long-format for faceted plots -------------------------------------
tolong <- function(d, vars, labels = vars) {
  d %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    mutate(var = factor(var, levels = vars, labels = labels))
}

lineplot <- function(d, vars, labels, title, ylab = "", logy = FALSE) {
  p <- ggplot(tolong(d, vars, labels), aes(time, value)) +
    geom_line(linewidth = 0.8, colour = "#1f6f8b") +
    facet_wrap(~var, scales = "free_y") +
    labs(title = title, x = "time (days)", y = ylab) + THEME
  if (logy) p <- p + scale_y_log10()
  p
}

# ============================================================================
# UI
# ============================================================================
ui <- fluidPage(
  titlePanel("Hypophosphatasia (HPP) — QSP Model Explorer"),
  p(em(paste("ALPL / TNSALP deficiency → perivesicular pyrophosphate",
             "excess → competitive inhibition of hydroxyapatite growth.",
             "Research / education only — not for clinical use."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("pheno", "Phenotype (residual LOCAL TNSALP activity)",
                  choices = PHENO, selected = "infantile"),
      sliderInput("fracenz", "Override residual activity (fraction of normal)",
                  min = 0.005, max = 1.0, value = 0.03, step = 0.005),
      checkboxInput("useoverride", "Use the override slider", FALSE),
      numericInput("wt", "Body weight (kg)", 5, min = 2, max = 120, step = 0.5),
      numericInput("age0", "Age at t = 0 (days)", 14, min = 0, max = 25000),
      numericInput("days", "Simulation horizon (days)", 365, min = 30,
                   max = 5475, step = 30),
      hr(),
      h4("Asfotase alfa (SC)"),
      numericInput("dose", "Dose (mg/kg per injection)", 2, min = 0, max = 10,
                   step = 0.25),
      selectInput("perweek", "Injections per week",
                  choices = c("3 (label, TIW)" = 3, "6 (label)" = 6,
                              "7 (daily)" = 7, "1 (weekly)" = 1),
                  selected = 3),
      numericInput("dstart", "Start day", 14, min = 0, max = 3650),
      numericInput("dstop", "Stop day (withdrawal)", 3650, min = 0, max = 5475),
      sliderInput("adamax", "Anti-drug antibody ceiling", 0, 1, 0.2, 0.05),
      hr(),
      h4("Adjunctive / other therapy"),
      numericInput("pn", "Pyridoxine (mg/day)", 0, min = 0, max = 600, step = 25),
      checkboxInput("tptd", "Teriparatide (adult, off-label)", FALSE),
      checkboxInput("bp", "Bisphosphonate exposure (HARM arm)", FALSE),
      sliderInput("cain", "Dietary Ca / vitamin-D intake (x reference)",
                  0.3, 1.6, 1.0, 0.1),
      hr(),
      helpText(paste("Every scenario starts from the untreated steady state of",
                     "the chosen genotype, obtained by burn-in — no",
                     "“diseased baseline” is hand-set."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---- 1. patient profile ------------------------------------------
        tabPanel(
          "1. Patient profile",
          br(),
          fluidRow(column(6, tableOutput("profile")),
                   column(6, plotOutput("profpos", height = "320px"))),
          h4("Where this patient sits on the mineralization hyperbola"),
          p(paste("The curve is the model's closed-form relation between LOCAL",
                  "TNSALP activity and relative mineral apposition. The point is",
                  "this patient untreated; the triangle is the same patient on",
                  "the selected regimen at the end of the horizon.")),
          verbatimTextOutput("profnote")
        ),

        # ---- 2. PK ---------------------------------------------------------
        tabPanel(
          "2. PK / target delivery",
          br(), plotOutput("pk", height = "560px"),
          p(strong("The point of this tab:"),
            paste("plasma exposure and serum ALP are what a clinic measures,",
                  "but the effect compartment is the hydroxyapatite-BOUND pool",
                  "(EB / occupancy). The D10 tag is what puts the enzyme there,",
                  "and the bound pool's ~3.5-day half-life is why dosing",
                  "interval matters much less than total weekly dose."))
        ),

        # ---- 3. PPi economy ------------------------------------------------
        tabPanel(
          "3. PPi economy & mineralization",
          br(), plotOutput("ppi", height = "560px"),
          p(strong("The point of this tab:"),
            paste("perivesicular PPi is the causal variable and is not",
                  "measurable in patients; plasma PPi is measurable and moves",
                  "much less. Mrel = 1 is healthy mineral apposition capacity;",
                  "values above 1 mean PPi has been suppressed BELOW normal,",
                  "which is where the ectopic-calcification term switches on."))
        ),

        # ---- 4. biochemical markers ---------------------------------------
        tabPanel(
          "4. Biochemical markers (T2)",
          br(), plotOutput("markers", height = "560px"),
          tableOutput("markertab"),
          p(strong("The point of this tab:"),
            paste("PLP, PEA and plasma PPi all report PLASMA enzyme activity.",
                  "They normalise early and then saturate. A normal PLP",
                  "therefore does not establish that the skeletal dose is",
                  "adequate — see the Dose-ranging tab."))
        ),

        # ---- 5. skeletal ---------------------------------------------------
        tabPanel(
          "5. Skeletal endpoints",
          br(), plotOutput("skel", height = "620px"),
          p(strong("Two clocks (T3):"),
            paste("RSS and osteoid volume are reversible. GPD (irreversible",
                  "growth-plate damage) and the craniosynostosis index only",
                  "ever increase; the height-Z trace inherits that floor."))
        ),

        # ---- 6. respiratory / survival ------------------------------------
        tabPanel(
          "6. Respiratory & survival",
          br(), plotOutput("resp", height = "520px"),
          p(paste("RESP is an index; below the ventilation threshold (0.5) the",
                  "model assumes invasive ventilation. The mortality hazard is",
                  "gated to infancy, so childhood- and adult-onset phenotypes",
                  "carry essentially no respiratory mortality here."))
        ),

        # ---- 7. mineral / renal -------------------------------------------
        tabPanel(
          "7. Mineral & renal",
          br(), plotOutput("mineral", height = "600px"),
          p(strong("Paradoxical hypercalcaemia:"),
            paste("calcium that is not deposited into bone stays in the",
                  "extracellular fluid. PTH is suppressed, urinary calcium",
                  "rises, and nephrocalcinosis accrues. Restoring",
                  "mineralization corrects it because the skeleton IS the",
                  "calcium sink; dietary restriction corrects the chemistry",
                  "without fixing the bone (and slightly worsens it)."))
        ),

        # ---- 8. CNS / B6 ---------------------------------------------------
        tabPanel(
          "8. CNS / vitamin B6",
          br(), plotOutput("cns", height = "480px"),
          p(strong("Substrate excess, product deficiency:"),
            paste("PLP cannot cross the blood-brain barrier; it must first be",
                  "dephosphorylated by TNSALP. So plasma PLP is HIGH while the",
                  "brain cofactor pool is LOW, which is why the seizures are",
                  "pyridoxine-responsive. Asfotase alfa is bone-targeted and",
                  "CNS-excluded, so it lowers plasma PLP without supplying the",
                  "brain — a mechanistic argument for continuing B6",
                  "monitoring on ERT. Model-generated, not demonstrated."))
        ),

        # ---- 9. scenario comparison ---------------------------------------
        tabPanel(
          "9. Scenario comparison",
          br(),
          checkboxGroupInput(
            "scn", "Predefined scenarios",
            choices = setNames(names(HPP_scenarios()),
                               vapply(HPP_scenarios(), function(z) z$label, "")),
            selected = c("S2", "S5", "S8", "S10"), inline = FALSE),
          selectInput("scnvar", "Readout",
                      choices = c("RSS" = "RSS", "Mrel" = "MRELN",
                                  "osteoid volume fraction" = "OST",
                                  "perivesicular PPi (uM)" = "PPIL",
                                  "plasma PLP (nM)" = "PLPP",
                                  "respiratory index" = "RESP",
                                  "survival" = "SURV",
                                  "height-Z change" = "HTZ",
                                  "irreversible GP damage" = "GPD",
                                  "muscle function" = "MUS",
                                  "serum ALP (U/L)" = "ALPUL"),
                      selected = "RSS"),
          plotOutput("scnplot", height = "480px"),
          tableOutput("scntab")
        ),

        # ---- 10. threshold explorer ---------------------------------------
        tabPanel(
          "10. Threshold explorer (T1)",
          br(), plotOutput("thresh", height = "440px"),
          tableOutput("threshtab"),
          p(paste("One hyperbola, read at different residual activities,",
                  "reproduces the whole clinical spectrum. The steep region",
                  "lies below ~0.2 of normal activity: that is why a",
                  "heterozygous carrier at 0.5 is nearly normal while an",
                  "infant at 0.03 is not, and why small absolute gains in",
                  "activity produce large gains only in the severe range."))
        ),

        # ---- 11. dose ranging ---------------------------------------------
        tabPanel(
          "11. Dose-ranging (T2)",
          br(),
          numericInput("swdays", "Read-out day", 182, min = 28, max = 1095,
                       step = 14),
          plotOutput("sweep", height = "460px"),
          tableOutput("sweeptab"),
          p(paste("Each readout is normalised to its own maximum achievable",
                  "change across the sweep. Plasma PLP reaches ~99% of its",
                  "maximum at roughly a tenth of the label weekly dose;",
                  "mineralization capacity is still climbing at five times the",
                  "label dose. The radiographic endpoint sits in between."))
        ),

        # ---- 12. notes -----------------------------------------------------
        tabPanel(
          "12. Model notes",
          br(),
          h4("Structure"),
          tags$ul(
            tags$li("32 differential states; 4 fast species (perivesicular PPi, plasma PPi, plasma PLP, plasma PEA) solved in closed form at quasi-steady state."),
            tags$li("Genotype enters as ONE number: residual local TNSALP activity."),
            tags$li("Baselines are obtained by burn-in of the untreated genotype, not hand-set.")
          ),
          h4("What is calibrated, and to what"),
          tags$ul(
            tags$li("Asfotase alfa F = 0.458, tmax 24-48 h, terminal t1/2 = 2.28 d: label values."),
            tags$li("Healthy plasma PPi ~3 uM; Ki(PPi) for hydroxyapatite growth low-micromolar: classical calcification literature."),
            tags$li("Osteoid volume fraction 2% healthy vs 9-11% infantile: emerges, not fitted."),
            tags$li("One-year survival: ~12% untreated perinatal, ~51% untreated infantile, ~89% on label dose (reported ~95% treated vs 27-42% historical control).")
          ),
          h4("What is NOT calibrated"),
          tags$ul(
            tags$li("Bone-binding parameters (KON/KOFF/KDEGB/BMAXKG) are unpublished; treat HA occupancy as a relative scale."),
            tags$li("Craniosynostosis is a rate law, not a mechanism."),
            tags$li("The calcium homeostat runs near its excretory capacity and is the most parameter-sensitive module."),
            tags$li("The prediction that ERT lowers the CNS B6 substrate gradient is model-generated and unproven.")
          ),
          h4("Files"),
          tags$ul(
            tags$li("hpp_qsp_model.dot / .svg / .png — mechanistic map (167 nodes, 16 clusters)"),
            tags$li("hpp_mrgsolve_model.R — this model"),
            tags$li("hpp_reference_model.py — stdlib-only reference implementation (numerical authority)"),
            tags$li("hpp_model_report.txt — full computed output"),
            tags$li("hpp_references.md — literature with PubMed links")
          )
        )
      )
    )
  )
)

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  ## keep weight / age sensible when the phenotype changes
  observeEvent(input$pheno, {
    g <- HPP_genotype(input$pheno)
    updateNumericInput(session, "wt", value = g$WT)
    updateNumericInput(session, "age0", value = g$AGE0)
    updateSliderInput(session, "fracenz", value = g$FRACENZ)
  })

  pars <- reactive({
    g <- HPP_genotype(input$pheno)
    if (isTRUE(input$useoverride)) g$FRACENZ <- input$fracenz
    g$WT <- input$wt
    g$AGE0 <- input$age0
    g
  })

  sim <- reactive({
    g <- pars()
    HPP_run(sev = input$pheno, days = input$days, dose = input$dose,
            per_week = as.numeric(input$perweek), dstart = input$dstart,
            dstop = input$dstop, pn = input$pn, ada_max = input$adamax,
            tptd = if (isTRUE(input$tptd)) 0.8 else 0,
            bp = if (isTRUE(input$bp)) 0.014 else 0,
            cain_mult = input$cain,
            extra = list(FRACENZ = g$FRACENZ, WT = g$WT, AGE0 = g$AGE0))
  })

  simbase <- reactive({
    g <- pars()
    HPP_run(sev = input$pheno, days = input$days, dose = 0,
            cain_mult = input$cain,
            extra = list(FRACENZ = g$FRACENZ, WT = g$WT, AGE0 = g$AGE0))
  })

  # ---- 1. profile ----------------------------------------------------------
  output$profile <- renderTable({
    d <- sim(); b <- simbase(); e <- d[nrow(d), ]; e0 <- b[1, ]; eb <- b[nrow(b), ]
    data.frame(
      Quantity = c("residual LOCAL activity (fraction)",
                   "serum ALP, untreated (U/L)",
                   "perivesicular PPi, untreated (uM)",
                   "plasma PPi, untreated (uM)",
                   "plasma PLP, untreated (nM)",
                   "osteoid volume, untreated (%)",
                   "RSS, untreated (end of horizon)",
                   "--- on the selected regimen ---",
                   "serum ALP (U/L)", "HA occupancy",
                   "LOCAL activity ELOC", "perivesicular PPi (uM)",
                   "plasma PLP (nM)", "Mrel", "osteoid volume (%)",
                   "RSS", "respiratory index", "survival"),
      Value = c(sprintf("%.3f", pars()$FRACENZ),
                sprintf("%.0f", e0$ALPUL), sprintf("%.1f", e0$PPIL),
                sprintf("%.2f", e0$PPIP), sprintf("%.0f", e0$PLPP),
                sprintf("%.1f", 100 * e0$OST), sprintf("%.2f", eb$RSS), "",
                sprintf("%.0f", e$ALPUL), sprintf("%.3f", e$OCC),
                sprintf("%.2f", e$ELOC), sprintf("%.2f", e$PPIL),
                sprintf("%.1f", e$PLPP), sprintf("%.3f", e$MRELN),
                sprintf("%.1f", 100 * e$OST), sprintf("%.2f", e$RSS),
                sprintf("%.2f", e$RESP), sprintf("%.3f", e$SURV)),
      check.names = FALSE)
  })

  output$profpos <- renderPlot({
    p <- as.list(param(HPP_mod))
    E <- 10 ^ seq(log10(0.005), log10(2), length.out = 200)
    mrel <- (1 + p$PPI0 / p$KIPPI) /
      (1 + (p$JPPI + p$KOUT * p$PPI0) / (p$KIPPI * (p$KCATP * E + p$KOUT)))
    d <- sim(); e <- d[nrow(d), ]
    ggplot(data.frame(E, mrel), aes(E, mrel)) +
      geom_line(linewidth = 1, colour = "#1f6f8b") +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
      geom_point(data = data.frame(E = pars()$FRACENZ,
                                   mrel = simbase()[1, ]$MRELN),
                 size = 4, colour = "#b03030") +
      geom_point(data = data.frame(E = e$ELOC, mrel = e$MRELN),
                 size = 4, shape = 17, colour = "#2c7a3f") +
      scale_x_log10() +
      labs(title = "Mineralization capacity vs LOCAL enzyme activity",
           x = "LOCAL TNSALP activity (fraction of normal, log scale)",
           y = "Mrel (1 = healthy)") + THEME
  })

  output$profnote <- renderText({
    d <- sim(); e <- d[nrow(d), ]; b <- simbase()[1, ]
    sprintf(paste0("Untreated: local activity %.3f, Mrel %.3f, perivesicular PPi %.1f uM ",
                   "(%.1fx normal) while plasma PPi is only %.1fx normal.\n",
                   "On regimen: local activity %.2f, Mrel %.3f, perivesicular PPi %.2f uM.\n",
                   "Irreversible growth-plate damage accrued: %.4f. ",
                   "Craniosynostosis index: %.3f."),
            pars()$FRACENZ, b$MRELN, b$PPIL, b$PPIL / 3, b$PPIP / 3,
            e$ELOC, e$MRELN, e$PPIL, e$GPD, e$CRAN)
  })

  # ---- 2. PK ---------------------------------------------------------------
  output$pk <- renderPlot({
    lineplot(sim(),
             c("Cp", "AC", "EB", "OCC", "ALPUL", "ELOC", "ADA", "EPL"),
             c("plasma asfotase alfa (mg/L)", "plasma amount (mg)",
               "HA-BOUND enzyme (mg)", "HA occupancy (0-1)",
               "serum ALP (U/L)", "LOCAL activity ELOC (1 = healthy)",
               "ADA titre", "PLASMA activity EPL (1 = healthy)"),
             "Asfotase alfa: exposure vs target delivery")
  })

  # ---- 3. PPi --------------------------------------------------------------
  output$ppi <- renderPlot({
    lineplot(sim(),
             c("PPIL", "PPIP", "MINH", "MRELN", "MAR", "OPN", "SSATN", "ECT"),
             c("PERIVESICULAR PPi (uM)", "plasma PPi (uM)",
               "inhibition factor MINH", "Mrel (1 = healthy)",
               "mineral apposition rate", "phospho-osteopontin",
               "Ca x Pi supersaturation", "ectopic calcification index"),
             "Pyrophosphate economy and the mineralization front")
  })

  # ---- 4. markers ----------------------------------------------------------
  output$markers <- renderPlot({
    lineplot(sim(), c("PLPP", "PEAP", "PPIP", "ALPUL", "PPIL", "MRELN"),
             c("plasma PLP (nM)", "plasma PEA (uM)", "plasma PPi (uM)",
               "serum ALP (U/L)", "PERIVESICULAR PPi (uM)", "Mrel"),
             "Measurable markers (top) vs causal quantities (bottom)")
  })

  output$markertab <- renderTable({
    d <- sim(); e <- d[nrow(d), ]; b <- simbase()[1, ]
    data.frame(
      Marker = c("plasma PLP (nM)", "plasma PEA (uM)", "plasma PPi (uM)",
                 "serum ALP (U/L)", "PERIVESICULAR PPi (uM)", "Mrel"),
      Untreated = sprintf("%.2f", c(b$PLPP, b$PEAP, b$PPIP, b$ALPUL, b$PPIL, b$MRELN)),
      OnTherapy = sprintf("%.2f", c(e$PLPP, e$PEAP, e$PPIP, e$ALPUL, e$PPIL, e$MRELN)),
      FoldChange = sprintf("%.2f", c(e$PLPP / b$PLPP, e$PEAP / b$PEAP,
                                     e$PPIP / b$PPIP, e$ALPUL / b$ALPUL,
                                     e$PPIL / b$PPIL, e$MRELN / b$MRELN)),
      check.names = FALSE)
  })

  # ---- 5. skeletal ---------------------------------------------------------
  output$skel <- renderPlot({
    lineplot(sim(),
             c("RSS", "OST", "BMIN", "GPREL", "GPD", "HTZ", "DENTREL", "FX",
               "CRAN"),
             c("RSS (0-10)", "osteoid volume fraction",
               "mineralised matrix (1 = normal)",
               "growth-plate integrity (rel.)",
               "IRREVERSIBLE GP damage", "height-Z change",
               "dental attachment (rel.)", "fracture burden",
               "craniosynostosis index"),
             "Skeletal, dental and stature endpoints")
  })

  # ---- 6. respiratory ------------------------------------------------------
  output$resp <- renderPlot({
    lineplot(sim(), c("RIBREL", "RESP", "VENT", "SURV", "MUS", "SIXMWT"),
             c("thoracic mineralization (rel.)", "respiratory index",
               "ventilation assumed (0/1)", "survival probability",
               "muscle function index", "6MWT-like (% predicted)"),
             "Thorax, respiration and survival")
  })

  # ---- 7. mineral ----------------------------------------------------------
  output$mineral <- renderPlot({
    lineplot(sim(), c("CAS", "PIS", "PTHS", "VITD", "JURN", "NEPH", "SSATN",
                      "ECT"),
             c("serum Ca (mmol/L)", "serum Pi (mmol/L)", "PTH (pmol/L)",
               "1,25(OH)2D (pmol/L)", "urinary Ca (mmol/kg/d)",
               "nephrocalcinosis", "Ca x Pi supersaturation",
               "ectopic calcification"),
             "Systemic mineral homeostasis and the kidney")
  })

  # ---- 8. CNS --------------------------------------------------------------
  output$cns <- renderPlot({
    lineplot(sim(), c("PLPP", "PLBR", "SEIZ", "MUS"),
             c("plasma PLP (nM) — HIGH",
               "CNS pyridoxal cofactor (1 = healthy) — LOW",
               "seizure burden index", "muscle function index"),
             "Vitamin B6: substrate excess with product deficiency")
  })

  # ---- 9. scenarios --------------------------------------------------------
  scnruns <- reactive({
    req(length(input$scn) > 0)
    scn <- HPP_scenarios()
    out <- lapply(input$scn, function(k) {
      d <- do.call(HPP_run, scn[[k]]$args)
      d$scenario <- scn[[k]]$label
      d
    })
    bind_rows(out)
  })

  output$scnplot <- renderPlot({
    d <- scnruns()
    ggplot(d, aes(time, .data[[input$scnvar]], colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = paste("Scenario comparison:", input$scnvar),
           x = "time (days)", y = input$scnvar, colour = NULL) +
      THEME + guides(colour = guide_legend(ncol = 2))
  })

  output$scntab <- renderTable({
    d <- scnruns()
    d %>%
      group_by(scenario) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      transmute(scenario, day = time, PPi_local = PPIL, Mrel = MRELN,
                osteoid = OST, RSS, dHeightZ = HTZ, RESP, survival = SURV,
                PLP = PLPP, seizure = SEIZ) %>%
      as.data.frame()
  }, digits = 3)

  # ---- 10. threshold -------------------------------------------------------
  threshdat <- reactive({
    p <- as.list(param(HPP_mod))
    E <- 10 ^ seq(log10(0.005), log10(2), length.out = 120)
    ppi_lin <- (p$JPPI + p$KOUT * p$PPI0) / (p$KCATP * E + p$KOUT)
    mrel_lin <- (1 + p$PPI0 / p$KIPPI) / (1 + ppi_lin / p$KIPPI)
    ## exact quadratic QSS at the untreated plasma activity of that genotype
    epl <- E * (1 - p$FISO) + p$FISO
    c1 <- p$FVOL * p$KOUT; c2 <- p$KCATPP * epl + p$CLPPI
    al <- p$KOUT * c2 / (c1 + c2); be <- p$KOUT * p$JPSYS / (c1 + c2)
    Vm <- p$KCATP * p$KMPPI; S <- p$JPPI + be
    bq <- Vm * E + al * p$KMPPI - S
    ppi_ex <- (-bq + sqrt(bq ^ 2 + 4 * al * S * p$KMPPI)) / (2 * al)
    mrel_ex <- (1 + p$PPI0 / p$KIPPI) / (1 + ppi_ex / p$KIPPI)
    data.frame(E, ppi_lin, ppi_ex, mrel_lin, mrel_ex)
  })

  output$thresh <- renderPlot({
    d <- threshdat()
    dd <- bind_rows(
      data.frame(E = d$E, value = d$mrel_ex, form = "exact (quadratic QSS)"),
      data.frame(E = d$E, value = d$mrel_lin, form = "linearised hyperbola"))
    marks <- data.frame(
      E = c(0.01, 0.03, 0.12, 0.35, 0.50),
      lab = c("perinatal", "infantile", "childhood", "adult", "carrier"))
    ggplot(dd, aes(E, value, linetype = form)) +
      geom_line(linewidth = 0.9, colour = "#1f6f8b") +
      geom_vline(data = marks, aes(xintercept = E), colour = "grey60",
                 linetype = 3) +
      geom_text(data = marks, aes(x = E, y = 1.65, label = lab), angle = 90,
                vjust = -0.3, hjust = 1, size = 3, inherit.aes = FALSE) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      scale_x_log10() +
      labs(title = "T1: one hyperbola generates the whole severity spectrum",
           x = "LOCAL TNSALP activity (fraction of normal, log scale)",
           y = "Mrel (1 = healthy)", linetype = NULL) + THEME
  })

  output$threshtab <- renderTable({
    d <- threshdat()
    idx <- vapply(c(0.01, 0.03, 0.12, 0.35, 0.50, 1.0),
                  function(x) which.min(abs(d$E - x)), 1L)
    data.frame(phenotype = c("perinatal", "infantile", "childhood", "adult",
                             "carrier", "healthy"),
               fracEnz = d$E[idx],
               PPi_local_exact = d$ppi_ex[idx],
               PPi_local_linear = d$ppi_lin[idx],
               Mrel_exact = d$mrel_ex[idx],
               Mrel_linear = d$mrel_lin[idx])
  }, digits = 3)

  # ---- 11. dose ranging ----------------------------------------------------
  sweepdat <- reactive({
    weekly <- c(0, 0.375, 0.75, 1.5, 3, 6, 9, 12, 18, 30)
    g <- pars()
    rows <- lapply(weekly, function(w) {
      d <- HPP_run(sev = input$pheno, days = input$swdays, dose = w / 3,
                   per_week = 3, delta = 7,
                   extra = list(FRACENZ = g$FRACENZ, WT = g$WT, AGE0 = g$AGE0))
      e <- d[nrow(d), ]
      data.frame(weekly = w, ALP = e$ALPUL, occupancy = e$OCC, ELOC = e$ELOC,
                 PLP = e$PLPP, PPi_plasma = e$PPIP, PPi_local = e$PPIL,
                 Mrel = e$MRELN, RSS = e$RSS, osteoid = e$OST)
    })
    bind_rows(rows)
  })

  output$sweep <- renderPlot({
    d <- sweepdat()
    fr <- function(v) (v - v[1]) / (v[length(v)] - v[1])
    dd <- bind_rows(
      data.frame(weekly = d$weekly, f = fr(d$PLP), readout = "plasma PLP"),
      data.frame(weekly = d$weekly, f = fr(d$PPi_plasma),
                 readout = "plasma PPi"),
      data.frame(weekly = d$weekly, f = fr(-d$RSS), readout = "RSS"),
      data.frame(weekly = d$weekly, f = fr(d$Mrel),
                 readout = "perivesicular Mrel"))
    ggplot(dd, aes(weekly, f, colour = readout)) +
      geom_line(linewidth = 1) + geom_point() +
      geom_vline(xintercept = 6, linetype = 2, colour = "grey40") +
      annotate("text", x = 6, y = 0.05, label = " label dose 6 mg/kg/wk",
               hjust = 0, size = 3.4) +
      labs(title = "T2: fraction of each readout's own maximum effect",
           x = "total weekly dose (mg/kg/week, given TIW)",
           y = "fraction of maximum achievable change", colour = NULL) + THEME
  })

  output$sweeptab <- renderTable({ sweepdat() }, digits = 3)
}

shinyApp(ui, server)
