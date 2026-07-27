##  Acute Respiratory Distress Syndrome (ARDS) — QSP Simulator (Shiny)
##  ============================================================================
##  Front-end for ards_mrgsolve_model.R (38-ODE ARDS QSP model).
##
##  Run with:
##      shiny::runApp("ards_shiny_app.R")
##  from inside the acute-respiratory-distress-syndrome/ directory, or
##      shiny::runApp("acute-respiratory-distress-syndrome/ards_shiny_app.R")
##  from the repository root.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##
##  Design intent
##  -------------
##  The clinically important thing about ARDS is that the ventilator setting and
##  the disease are the same problem: the aerated lung the ventilator is
##  breathing into is itself a function of how the ventilator has been breathing
##  into it. The app is therefore laid out so that the ventilator controls sit
##  next to the mechanics readout on every tab, and so that the "Mechanics &
##  VILI" tab shows driving pressure, strain and mechanical power as the three
##  exposures that feed the injury loop.
##
##  Tabs
##    1  Patient & phenotype   — presenting PaO2/FiO2, insult type, latent class
##    2  Ventilator & mechanics— Vt/PEEP/RR/FiO2, prone, ECMO, the baby lung
##    3  Gas exchange          — PaO2/FiO2, shunt, dead space, PaCO2, Berlin class
##    4  Lung water & barrier  — EVLW, permeability, AFC, epithelium/endothelium
##    5  Inflammation          — cytokines, neutrophils, NETs, macrophages
##    6  Drug PK / PD          — steroid effect, NMB depth, tocilizumab, diuretic
##    7  Clinical endpoints    — mortality, ventilator-free days, SOFA, weakness
##    8  Scenario comparison   — run several arms side by side with a summary
##  ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
MODEL_FILE <- "ards_mrgsolve_model.R"
if (!file.exists(MODEL_FILE)) {
  alt <- file.path("acute-respiratory-distress-syndrome", MODEL_FILE)
  if (file.exists(alt)) MODEL_FILE <- alt
}
mod <- mrgsolve::mread_cache("ards", project = dirname(normalizePath(MODEL_FILE)),
                             file = basename(MODEL_FILE))

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

PAL <- c("#2c6fbb", "#c0392b", "#27884a", "#8e44ad", "#d98a00",
         "#16a085", "#7f8c8d", "#c2185b")

# ---------------------------------------------------------------------------
# Intervention presets — everything a clinician actually turns
# ---------------------------------------------------------------------------
VENT_PRESETS <- list(
  "Lung-protective (6 mL/kg, PEEP 12)" = list(VT_KG = 6,  PEEP = 12, RR = 28),
  "Conventional 1990s (12 mL/kg)"      = list(VT_KG = 12, PEEP = 8,  RR = 16),
  "Ultraprotective (4 mL/kg)"          = list(VT_KG = 4,  PEEP = 12, RR = 32),
  "High PEEP (6 mL/kg, PEEP 16)"       = list(VT_KG = 6,  PEEP = 16, RR = 28)
)

build_events <- function(input) {
  ev_list <- list()
  if (isTRUE(input$use_dex)) {
    ev_list <- c(ev_list, list(
      ev(amt = input$dex_hi, cmt = "DEX", ii = 24, addl = input$dex_hi_d - 1, time = 0),
      ev(amt = input$dex_lo, cmt = "DEX", ii = 24, addl = input$dex_lo_d - 1,
         time = 24 * input$dex_hi_d)))
  }
  if (isTRUE(input$use_mp)) {
    ev_list <- c(ev_list, list(
      ev(amt = input$mp_dose, cmt = "MPRED", ii = 24, addl = input$mp_days - 1,
         time = 24 * input$mp_start)))
  }
  if (isTRUE(input$use_nmba)) {
    ev_list <- c(ev_list, list(
      ev(amt = 37.5 * input$nmba_h, rate = 37.5, cmt = "CIS", time = 0)))
  }
  if (isTRUE(input$use_tcz)) {
    ev_list <- c(ev_list, list(ev(amt = input$tcz_dose, cmt = "TCZ", time = 0)))
  }
  if (isTRUE(input$use_fur)) {
    ev_list <- c(ev_list, list(
      ev(amt = input$fur_dose, cmt = "FUR", ii = 6, addl = 4 * input$fur_days - 1,
         time = 0)))
  }
  if (length(ev_list) == 0) return(NULL)
  Reduce(c, ev_list)
}

build_params <- function(input) {
  vp <- VENT_PRESETS[[input$vent_preset]]
  list(
    PF0      = input$pf0,
    PBW      = input$pbw,
    DIRECT   = input$direct,
    HYPER    = as.numeric(input$hyper == "Hyperinflammatory"),
    INSULT0  = input$insult0,
    KCLR_P   = input$kclr,
    RECRUIT  = input$recruit,
    ANTIOX   = input$antiox,
    AGE_RISK = input$age_risk,
    VT_KG    = if (input$vent_manual) input$vt_kg else vp$VT_KG,
    PEEP     = if (input$vent_manual) input$peep  else vp$PEEP,
    RR       = if (input$vent_manual) input$rr    else vp$RR,
    FIO2     = input$fio2,
    SEDATION = input$sedation,
    PRONE_H  = if (isTRUE(input$use_prone)) input$prone_h else 0,
    ECMO     = as.numeric(isTRUE(input$use_ecmo)),
    NO_PPM   = if (isTRUE(input$use_ino)) input$ino_ppm else 0,
    FLUID_IN = input$fluid_in,
    BARI     = as.numeric(isTRUE(input$use_bari)),
    BETA2    = as.numeric(isTRUE(input$use_beta2)),
    MOBILISE = input$mobilise
  )
}

simulate_arm <- function(pars, events, end_h = 672) {
  m <- param(mod, pars)
  out <- if (is.null(events)) mrgsim(m, end = end_h, delta = 2)
         else mrgsim(m, events = events, end = end_h, delta = 2)
  as.data.frame(out) %>% mutate(day = time / 24)
}

endpoint_row <- function(d, label) {
  # mrgsim emits duplicated time points at dose records; keep the last of each
  d <- d[!duplicated(d$time, fromLast = TRUE), ]
  n <- nrow(d)
  data.frame(
    Arm               = label,
    `PF at 24 h`      = round(approx(d$time, d$PaO2_FiO2, 24)$y, 0),
    `PF at day 7`     = round(approx(d$time, d$PaO2_FiO2, 168)$y, 0),
    `Driving P (t0)`  = round(d$DrivingP_cmH2O[1], 1),
    `Mech power (t0)` = round(d$MechPower_Jmin[1], 1),
    `Peak VILI`       = round(max(d$VILI_burden), 2),
    `EVLW day 7`      = round(approx(d$time, d$EVLW_mlkg, 168)$y, 1),
    `Fluid bal day 7` = round(approx(d$time, d$Fluid_balance_L, 168)$y, 1),
    `IL-6 day 7`      = round(approx(d$time, d$IL6_pgml, 168)$y, 0),
    `SOFA day 7`      = round(approx(d$time, d$SOFA_nonpulm, 168)$y, 1),
    `28-d mortality`  = sprintf("%.1f%%", 100 * (1 - d$Survival_prob[n])),
    `VFD (expected)`  = round(d$VFD_expected[n], 1),
    `Weakness`        = round(d$Weakness[n], 2),
    check.names = FALSE
  )
}

lineplot <- function(d, vars, labels, ylab, hlines = NULL) {
  dd <- d %>% select(day, all_of(vars)) %>%
    pivot_longer(-day, names_to = "var", values_to = "value") %>%
    mutate(var = factor(var, levels = vars, labels = labels))
  p <- ggplot(dd, aes(day, value, colour = var)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = PAL, name = NULL) +
    labs(x = "Day", y = ylab) + THEME
  if (!is.null(hlines)) {
    for (h in hlines) p <- p + geom_hline(yintercept = h, linetype = 2,
                                          colour = "grey40", linewidth = 0.4)
  }
  p
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("ARDS — Quantitative Systems Pharmacology Simulator"),
  helpText(HTML(
    "38-compartment QSP model of the acute respiratory distress syndrome. ",
    "Barrier failure sets the lung water; the lung water sets the aerated ",
    "<i>baby lung</i>; the baby lung sets the driving pressure, strain and ",
    "mechanical power that the ventilator delivers; and that mechanical ",
    "exposure feeds back into the inflammation as VILI. ",
    "<b>Educational / research use only — not for clinical decision-making.</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Presenting patient"),
      sliderInput("pf0", "PaO2/FiO2 at enrolment (mmHg)", 60, 300, 130, step = 10),
      sliderInput("pbw", "Predicted body weight (kg)", 40, 90, 60, step = 1),
      sliderInput("direct", "Injury pattern (0 = indirect / endothelial, 1 = direct / epithelial)",
                  0, 1, 0.8, step = 0.05),
      radioButtons("hyper", "Latent inflammatory class",
                   c("Hypoinflammatory", "Hyperinflammatory"), inline = TRUE),
      sliderInput("recruit", "Lung recruitability", 0.1, 0.95, 0.55, step = 0.05),
      sliderInput("insult0", "Insult burden at t0", 0.2, 3, 1, step = 0.1),
      sliderInput("kclr", "Insult clearance rate (1/h)", 0.002, 0.06, 0.02, step = 0.002),
      sliderInput("antiox", "Antioxidant reserve", 0.4, 1.2, 1.0, step = 0.05),
      sliderInput("age_risk", "Age / comorbidity hazard multiplier", 0.5, 2.5, 1.0, step = 0.1),

      hr(), h4("Ventilator"),
      selectInput("vent_preset", "Preset", names(VENT_PRESETS)),
      checkboxInput("vent_manual", "Override preset manually", FALSE),
      conditionalPanel("input.vent_manual == true",
        sliderInput("vt_kg", "Tidal volume (mL/kg PBW)", 3, 14, 6, step = 0.5),
        sliderInput("peep", "PEEP (cmH2O)", 0, 24, 12, step = 1),
        sliderInput("rr", "Respiratory rate (/min)", 8, 40, 28, step = 1)),
      sliderInput("fio2", "FiO2", 0.21, 1.0, 0.7, step = 0.05),
      sliderInput("sedation", "Sedation depth (suppresses P-SILI)", 0, 1, 0.5, step = 0.05),
      checkboxInput("use_prone", "Prone positioning", FALSE),
      conditionalPanel("input.use_prone == true",
        sliderInput("prone_h", "Hours prone per day", 4, 20, 16, step = 1)),
      checkboxInput("use_ecmo", "VV-ECMO (ultraprotective settings)", FALSE),
      sliderInput("fluid_in", "Net fluid intake (mL/h)", 30, 200, 125, step = 5),
      sliderInput("mobilise", "Early mobilisation intensity", 0, 1, 0.3, step = 0.05),

      hr(), h4("Pharmacotherapy"),
      checkboxInput("use_dex", "Dexamethasone IV (DEXA-ARDS schedule)", FALSE),
      conditionalPanel("input.use_dex == true",
        numericInput("dex_hi", "High dose (mg/day)", 20, 4, 40, 2),
        numericInput("dex_hi_d", "High-dose days", 5, 1, 14, 1),
        numericInput("dex_lo", "Step-down dose (mg/day)", 10, 0, 20, 2),
        numericInput("dex_lo_d", "Step-down days", 5, 0, 14, 1)),
      checkboxInput("use_mp", "Methylprednisolone IV", FALSE),
      conditionalPanel("input.use_mp == true",
        numericInput("mp_dose", "Dose (mg/day)", 120, 20, 500, 10),
        numericInput("mp_days", "Days", 14, 1, 28, 1),
        numericInput("mp_start", "Start on day", 0, 0, 21, 1)),
      checkboxInput("use_nmba", "Cisatracurium infusion 37.5 mg/h", FALSE),
      conditionalPanel("input.use_nmba == true",
        sliderInput("nmba_h", "Infusion duration (h)", 12, 96, 48, step = 12)),
      checkboxInput("use_ino", "Inhaled NO / epoprostenol", FALSE),
      conditionalPanel("input.use_ino == true",
        sliderInput("ino_ppm", "Dose (ppm-equivalent)", 1, 40, 20, step = 1)),
      checkboxInput("use_tcz", "Tocilizumab (IL-6R blockade)", FALSE),
      conditionalPanel("input.use_tcz == true",
        numericInput("tcz_dose", "Dose (mg)", 480, 80, 800, 40)),
      checkboxInput("use_bari", "Baricitinib 4 mg/day (JAK1/2)", FALSE),
      checkboxInput("use_fur", "Furosemide 20 mg q6h", FALSE),
      conditionalPanel("input.use_fur == true",
        numericInput("fur_dose", "Dose per administration (mg)", 20, 5, 80, 5),
        numericInput("fur_days", "Days", 7, 1, 28, 1)),
      checkboxInput("use_beta2", "Inhaled beta2-agonist (BALTI-2)", FALSE),

      hr(),
      actionButton("go", "Simulate", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · Patient & phenotype",
          br(),
          fluidRow(column(6, h4("Solved presenting state")),
                   column(6, h4("Berlin class over time"))),
          fluidRow(column(6, tableOutput("baseline_tbl")),
                   column(6, plotOutput("berlin_plot", height = 260))),
          hr(),
          p(HTML("The model does not ask for the barrier integrities — it ",
                 "<b>solves for them</b>. Given the presenting PaO2/FiO2 it inverts the ",
                 "gas-exchange relation for the shunt, splits the non-aerated lung ",
                 "between flooding and atelectasis according to the direct/indirect ",
                 "injury pattern, and then solves the Starling balance backwards for the ",
                 "epithelial and endothelial integrities that hold that much lung water.")),
          plotOutput("insult_plot", height = 300)
        ),

        tabPanel("2 · Ventilator & mechanics",
          br(),
          fluidRow(column(4, uiOutput("box_dp")),
                   column(4, uiOutput("box_mp")),
                   column(4, uiOutput("box_aer"))),
          plotOutput("mech_plot", height = 300),
          hr(), h4("The injury loop: mechanical exposure and VILI burden"),
          plotOutput("vili_plot", height = 300),
          helpText("Dashed lines mark the thresholds the literature treats as ",
                   "inflection points: driving pressure 15 cmH2O, plateau 30 cmH2O, ",
                   "mechanical power 12 J/min, global strain 1.0.")
        ),

        tabPanel("3 · Gas exchange",
          br(),
          plotOutput("pf_plot", height = 300),
          fluidRow(column(6, plotOutput("shunt_plot", height = 280)),
                   column(6, plotOutput("co2_plot", height = 280))),
          helpText("Shunt drives hypoxaemia; dead space drives the CO2 and is an ",
                   "independent mortality signal (Nuckton 2002). On ECMO the reported ",
                   "PaO2/FiO2 includes the circuit — the ventilator-only trace is the ",
                   "one that decides decannulation.")
        ),

        tabPanel("4 · Lung water & barrier",
          br(),
          fluidRow(column(6, plotOutput("evlw_plot", height = 290)),
                   column(6, plotOutput("barrier_plot", height = 290))),
          hr(),
          fluidRow(column(6, plotOutput("afc_plot", height = 280)),
                   column(6, plotOutput("fluid_plot", height = 280)))
        ),

        tabPanel("5 · Inflammation",
          br(),
          plotOutput("cyto_plot", height = 300),
          fluidRow(column(6, plotOutput("neut_plot", height = 280)),
                   column(6, plotOutput("mac_plot", height = 280))),
          helpText("Note that tocilizumab RAISES measured IL-6 while lowering IL-6 ",
                   "signalling — the model reports both, because the assay and the ",
                   "biology diverge exactly where the drug works.")
        ),

        tabPanel("6 · Drug PK / PD",
          br(),
          fluidRow(column(6, plotOutput("gc_plot", height = 280)),
                   column(6, plotOutput("nmb_plot", height = 280))),
          hr(),
          fluidRow(column(6, plotOutput("tcz_plot", height = 280)),
                   column(6, plotOutput("adverse_plot", height = 280)))
        ),

        tabPanel("7 · Clinical endpoints",
          br(),
          fluidRow(column(3, uiOutput("box_mort")),
                   column(3, uiOutput("box_vfd")),
                   column(3, uiOutput("box_sofa")),
                   column(3, uiOutput("box_weak"))),
          plotOutput("outcome_plot", height = 300),
          hr(),
          fluidRow(column(6, plotOutput("sofa_plot", height = 280)),
                   column(6, plotOutput("rv_plot", height = 280))),
          hr(), h4("Endpoint summary"), tableOutput("endpoint_tbl")
        ),

        tabPanel("8 · Scenario comparison",
          br(),
          p("Runs the current patient through a standard set of arms, so the ",
            "comparison is always within the same virtual patient."),
          checkboxGroupInput("cmp_arms", "Arms to compare",
            choices = c("Conventional Vt 12", "LTVV Vt 6", "LTVV + high PEEP",
                        "LTVV + prone 16 h", "LTVV + NMBA 48 h",
                        "LTVV + dexamethasone", "Conservative fluid",
                        "Inhaled NO 20 ppm", "VV-ECMO", "Full bundle"),
            selected = c("Conventional Vt 12", "LTVV Vt 6", "LTVV + prone 16 h",
                         "LTVV + dexamethasone"),
            inline = TRUE),
          actionButton("cmp_go", "Run comparison", class = "btn-success"),
          br(), br(),
          plotOutput("cmp_plot", height = 380),
          hr(),
          DT::dataTableOutput("cmp_tbl")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    simulate_arm(build_params(input), build_events(input))
  })

  # --- Tab 1 -----------------------------------------------------------------
  output$baseline_tbl <- renderTable({
    d <- sim()[1, ]
    data.frame(
      Quantity = c("PaO2/FiO2 (mmHg)", "Shunt fraction", "Aerated lung fraction",
                   "Static compliance (mL/cmH2O)", "Driving pressure (cmH2O)",
                   "EVLW (mL/kg PBW)", "Epithelial integrity", "Endothelial integrity",
                   "Surfactant pool", "IL-6 (pg/mL)", "Non-pulmonary SOFA"),
      Value = c(sprintf("%.0f", d$PaO2_FiO2), sprintf("%.2f", d$Shunt_fraction),
                sprintf("%.2f", d$Aerated_frac), sprintf("%.1f", d$Compliance_Crs),
                sprintf("%.1f", d$DrivingP_cmH2O), sprintf("%.1f", d$EVLW_mlkg),
                sprintf("%.2f", d$Epi_integrity), sprintf("%.2f", d$Endo_integrity),
                sprintf("%.2f", d$Surfactant), sprintf("%.0f", d$IL6_pgml),
                sprintf("%.1f", d$SOFA_nonpulm)))
  }, striped = TRUE, width = "100%")

  output$berlin_plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, Berlin_class)) +
      geom_step(colour = PAL[2], linewidth = 1) +
      scale_y_continuous(breaks = 0:3,
        labels = c("none", "mild", "moderate", "severe"), limits = c(0, 3)) +
      labs(x = "Day", y = NULL, title = "Berlin severity class") + THEME
  })

  output$insult_plot <- renderPlot({
    lineplot(sim(),
      c("Epi_integrity", "Endo_integrity", "Surfactant", "AT2_pool"),
      c("Epithelial integrity", "Endothelial integrity", "Surfactant pool",
        "AT2 progenitor pool"),
      "Fraction of normal") +
      labs(title = "Barrier and repair-capacity trajectories")
  })

  # --- Tab 2 -----------------------------------------------------------------
  kpi <- function(label, value, sub, colour) {
    HTML(sprintf(
      "<div style='border:1px solid #ddd;border-left:6px solid %s;border-radius:6px;
                   padding:10px 14px;margin-bottom:6px;'>
         <div style='font-size:12px;color:#666'>%s</div>
         <div style='font-size:26px;font-weight:600'>%s</div>
         <div style='font-size:11px;color:#888'>%s</div></div>",
      colour, label, value, sub))
  }
  output$box_dp <- renderUI({
    d <- sim(); v <- d$DrivingP_cmH2O[1]
    kpi("Driving pressure at t0", sprintf("%.1f cmH2O", v),
        "target < 15", if (v > 15) PAL[2] else PAL[3])
  })
  output$box_mp <- renderUI({
    d <- sim(); v <- d$MechPower_Jmin[1]
    kpi("Mechanical power at t0", sprintf("%.1f J/min", v),
        "inflection ~12 J/min", if (v > 17) PAL[2] else PAL[3])
  })
  output$box_aer <- renderUI({
    d <- sim(); v <- d$Aerated_frac[1]
    kpi("Aerated lung ('baby lung')", sprintf("%.0f%%", 100 * v),
        "of total parenchyma", if (v < 0.3) PAL[2] else PAL[3])
  })

  output$mech_plot <- renderPlot({
    lineplot(sim(),
      c("Compliance_Crs", "DrivingP_cmH2O", "Plateau_cmH2O", "VT_delivered_ml"),
      c("Compliance (mL/cmH2O)", "Driving pressure (cmH2O)",
        "Plateau pressure (cmH2O)", "Delivered Vt (mL)"),
      "Value", hlines = c(15, 30)) +
      facet_wrap(~var, scales = "free_y") +
      theme(legend.position = "none") +
      labs(title = "Respiratory mechanics")
  })

  output$vili_plot <- renderPlot({
    lineplot(sim(),
      c("Strain_global", "MechPower_Jmin", "VILI_burden", "Effort_PSILI"),
      c("Global strain", "Mechanical power (J/min)", "Cumulative VILI burden",
        "Spontaneous effort (P-SILI)"),
      "Value", hlines = c(1, 12)) +
      facet_wrap(~var, scales = "free_y") +
      theme(legend.position = "none")
  })

  # --- Tab 3 -----------------------------------------------------------------
  output$pf_plot <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, PaO2_FiO2, PaO2_FiO2_vent) %>%
      pivot_longer(-day, names_to = "var", values_to = "value") %>%
      mutate(var = recode(var, PaO2_FiO2 = "Arterial (incl. ECMO circuit)",
                          PaO2_FiO2_vent = "Ventilator-only (native lung)"))
    ggplot(dd, aes(day, value, colour = var)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = c(100, 200, 300), linetype = 2, colour = "grey50") +
      annotate("text", x = 0.4, y = c(100, 200, 300) + 9, hjust = 0, size = 3,
               colour = "grey35", label = c("severe", "moderate", "mild")) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "Day", y = "PaO2/FiO2 (mmHg)", title = "Oxygenation") + THEME
  })
  output$shunt_plot <- renderPlot({
    lineplot(sim(), c("Shunt_fraction", "DeadSpace_VdVt"),
             c("Shunt fraction Qs/Qt", "Dead-space fraction Vd/Vt"),
             "Fraction", hlines = 0.6) +
      labs(title = "Shunt and dead space")
  })
  output$co2_plot <- renderPlot({
    lineplot(sim(), c("PaCO2_mmHg", "pH_arterial"),
             c("PaCO2 (mmHg)", "Arterial pH"), "Value") +
      facet_wrap(~var, scales = "free_y", ncol = 1) +
      theme(legend.position = "none") +
      labs(title = "CO2 and acid-base")
  })

  # --- Tab 4 -----------------------------------------------------------------
  output$evlw_plot <- renderPlot({
    lineplot(sim(), c("EVLW_mlkg", "Permeability"),
             c("EVLW (mL/kg PBW)", "Permeability multiplier"), "Value",
             hlines = 10) +
      facet_wrap(~var, scales = "free_y", ncol = 1) +
      theme(legend.position = "none") +
      labs(title = "Lung water and barrier permeability")
  })
  output$barrier_plot <- renderPlot({
    lineplot(sim(), c("Epi_integrity", "Endo_integrity", "Surfactant"),
             c("Epithelium", "Endothelium", "Surfactant"), "Fraction of normal") +
      labs(title = "Alveolar-capillary barrier")
  })
  output$afc_plot <- renderPlot({
    lineplot(sim(), c("AFC_per_h"), c("Alveolar fluid clearance (1/h)"),
             "Fraction cleared per hour", hlines = 0.03) +
      labs(title = "Alveolar fluid clearance (< 3%/h = impaired)")
  })
  output$fluid_plot <- renderPlot({
    lineplot(sim(), c("Fluid_balance_L"), c("Cumulative fluid balance (L)"),
             "Litres", hlines = 0) +
      labs(title = "Fluid balance")
  })

  # --- Tab 5 -----------------------------------------------------------------
  output$cyto_plot <- renderPlot({
    lineplot(sim(),
      c("IL6_pgml", "IL6_signalling", "IL8_pgml", "IL1b_pgml", "TNFa_pgml", "IL10_pgml"),
      c("IL-6 (measured)", "IL-6 signalling", "IL-8", "IL-1beta", "TNF-alpha", "IL-10"),
      "pg/mL") +
      scale_y_log10() +
      labs(title = "Cytokine trajectories (log scale)")
  })
  output$neut_plot <- renderPlot({
    lineplot(sim(), c("Neutrophils_alv", "NETs_burden", "ROS_burden"),
             c("Alveolar neutrophils", "NET / elastase burden", "Oxidant burden"),
             "Index") + labs(title = "Neutrophil-driven injury")
  })
  output$mac_plot <- renderPlot({
    lineplot(sim(), c("M1_programme", "M2_programme", "Ang2_ngml", "PAI1_ngml"),
             c("M1 macrophage", "M2 macrophage", "Ang-2 (ng/mL)", "PAI-1 (ng/mL)"),
             "Value") +
      facet_wrap(~var, scales = "free_y") +
      theme(legend.position = "none") +
      labs(title = "Resolution programme and endothelial / coagulation markers")
  })

  # --- Tab 6 -----------------------------------------------------------------
  output$gc_plot <- renderPlot({
    lineplot(sim(), c("GC_effect"), c("Glucocorticoid genomic effect"),
             "0-1") + labs(title = "Corticosteroid effect compartment")
  })
  output$nmb_plot <- renderPlot({
    lineplot(sim(), c("NMB_depth", "Effort_PSILI"),
             c("Neuromuscular block depth", "Spontaneous effort"), "0-1") +
      labs(title = "Neuromuscular blockade and P-SILI")
  })
  output$tcz_plot <- renderPlot({
    lineplot(sim(), c("TCZ_conc"), c("Tocilizumab (mg/L)"), "mg/L") +
      labs(title = "Tocilizumab concentration")
  })
  output$adverse_plot <- renderPlot({
    lineplot(sim(), c("Glucose_mgdl", "Weakness"),
             c("Glucose (mg/dL)", "ICU-acquired weakness"), "Value",
             hlines = 180) +
      facet_wrap(~var, scales = "free_y", ncol = 1) +
      theme(legend.position = "none") +
      labs(title = "Iatrogenic burden")
  })

  # --- Tab 7 -----------------------------------------------------------------
  output$box_mort <- renderUI({
    d <- sim(); v <- 1 - d$Survival_prob[nrow(d)]
    kpi("28-day mortality", sprintf("%.1f%%", 100 * v), "model estimate",
        if (v > 0.35) PAL[2] else PAL[3])
  })
  output$box_vfd <- renderUI({
    d <- sim(); v <- d$VFD_expected[nrow(d)]
    kpi("Ventilator-free days", sprintf("%.1f", v), "survival-weighted, day 28",
        if (v < 8) PAL[2] else PAL[3])
  })
  output$box_sofa <- renderUI({
    d <- sim(); v <- max(d$SOFA_nonpulm)
    kpi("Peak non-pulmonary SOFA", sprintf("%.1f", v), "points",
        if (v > 8) PAL[2] else PAL[3])
  })
  output$box_weak <- renderUI({
    d <- sim(); v <- d$Weakness[nrow(d)]
    kpi("ICU-acquired weakness", sprintf("%.2f", v), "index at day 28",
        if (v > 1.0) PAL[2] else PAL[3])
  })
  output$outcome_plot <- renderPlot({
    lineplot(sim(), c("Mortality_prob", "VFD_raw", "VFD_expected"),
             c("Cumulative mortality probability", "Ventilator-free days (raw)",
               "Ventilator-free days (survival-weighted)"), "Value") +
      facet_wrap(~var, scales = "free_y") +
      theme(legend.position = "none") +
      labs(title = "Outcome accrual")
  })
  output$sofa_plot <- renderPlot({
    lineplot(sim(), c("SOFA_nonpulm", "Murray_score"),
             c("Non-pulmonary SOFA", "Murray lung injury score"), "Points") +
      labs(title = "Severity scores")
  })
  output$rv_plot <- renderPlot({
    lineplot(sim(), c("PVR_dyn", "CardiacOutput"),
             c("PVR (dyn.s.cm-5)", "Cardiac output (L/min)"), "Value") +
      facet_wrap(~var, scales = "free_y", ncol = 1) +
      theme(legend.position = "none") +
      labs(title = "Right heart")
  })
  output$endpoint_tbl <- renderTable({
    endpoint_row(sim(), "Current settings")
  }, striped = TRUE, width = "100%")

  # --- Tab 8 -----------------------------------------------------------------
  arm_defs <- function(base_par) {
    ltvv <- modifyList(base_par, list(VT_KG = 6, PEEP = 12, RR = 28))
    list(
      "Conventional Vt 12"   = list(p = modifyList(base_par,
                                     list(VT_KG = 12, PEEP = 8, RR = 16)), e = NULL),
      "LTVV Vt 6"            = list(p = ltvv, e = NULL),
      "LTVV + high PEEP"     = list(p = modifyList(ltvv, list(PEEP = 16)), e = NULL),
      "LTVV + prone 16 h"    = list(p = modifyList(ltvv, list(PRONE_H = 16)), e = NULL),
      "LTVV + NMBA 48 h"     = list(p = modifyList(ltvv, list(SEDATION = 0.9)),
                                    e = ev(amt = 1800, rate = 37.5, cmt = "CIS", time = 0)),
      "LTVV + dexamethasone" = list(p = ltvv,
                                    e = c(ev(amt = 20, cmt = "DEX", ii = 24, addl = 4, time = 0),
                                          ev(amt = 10, cmt = "DEX", ii = 24, addl = 4, time = 120))),
      "Conservative fluid"   = list(p = modifyList(ltvv, list(FLUID_IN = 60)),
                                    e = ev(amt = 20, cmt = "FUR", ii = 6, addl = 27, time = 0)),
      "Inhaled NO 20 ppm"    = list(p = modifyList(ltvv, list(NO_PPM = 20)), e = NULL),
      "VV-ECMO"              = list(p = modifyList(ltvv, list(ECMO = 1)), e = NULL),
      "Full bundle"          = list(p = modifyList(ltvv,
                                     list(PRONE_H = 16, SEDATION = 0.9, FLUID_IN = 60)),
                                    e = c(ev(amt = 20, cmt = "DEX", ii = 24, addl = 4, time = 0),
                                          ev(amt = 10, cmt = "DEX", ii = 24, addl = 4, time = 120),
                                          ev(amt = 1800, rate = 37.5, cmt = "CIS", time = 0),
                                          ev(amt = 20, cmt = "FUR", ii = 6, addl = 27, time = 0)))
    )
  }

  cmp <- eventReactive(input$cmp_go, {
    base_par <- build_params(input)
    defs <- arm_defs(base_par)
    sel <- intersect(input$cmp_arms, names(defs))
    validate(need(length(sel) > 0, "Select at least one arm."))
    res <- lapply(sel, function(nm) {
      d <- simulate_arm(defs[[nm]]$p, defs[[nm]]$e)
      d$arm <- nm
      d
    })
    bind_rows(res)
  })

  output$cmp_plot <- renderPlot({
    d <- cmp()
    dd <- d %>%
      select(day, arm, PaO2_FiO2, DrivingP_cmH2O, EVLW_mlkg, VILI_burden,
             Mortality_prob, VFD_expected) %>%
      pivot_longer(-c(day, arm), names_to = "var", values_to = "value") %>%
      mutate(var = recode(var,
        PaO2_FiO2 = "PaO2/FiO2 (mmHg)", DrivingP_cmH2O = "Driving pressure (cmH2O)",
        EVLW_mlkg = "EVLW (mL/kg)", VILI_burden = "VILI burden",
        Mortality_prob = "Cumulative mortality", VFD_expected = "Ventilator-free days"))
    ggplot(dd, aes(day, value, colour = arm)) +
      geom_line(linewidth = 0.85) +
      facet_wrap(~var, scales = "free_y") +
      scale_colour_manual(values = rep(PAL, 2), name = NULL) +
      labs(x = "Day", y = NULL) + THEME
  })

  output$cmp_tbl <- DT::renderDataTable({
    d <- cmp()
    rows <- lapply(split(d, d$arm), function(x) endpoint_row(x, x$arm[1]))
    DT::datatable(bind_rows(rows), rownames = FALSE,
                  options = list(dom = "t", pageLength = 20, scrollX = TRUE))
  })
}

shinyApp(ui, server)
