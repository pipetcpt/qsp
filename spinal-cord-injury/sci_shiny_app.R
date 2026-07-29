## ============================================================================
## Traumatic Spinal Cord Injury (SCI) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 12 tabs: Patient & injury profile · Drug PK · Secondary cascade ·
##          Neuroinflammation · Tissue & imaging surrogates · Neurological
##          recovery · Therapeutic-window explorer · Surrogate-endpoint
##          dissociation · Secondary complications · Scenario comparison ·
##          Steroid safety · References
##
## The two tabs that carry the model's argument are "Therapeutic window"
## (efficacy collapses with start time while toxicity does not) and
## "Surrogate vs endpoint" (the same tissue benefit is worth a different
## number of ISNCSCI points depending on where the patient sits on the Hill
## curve). Everything else is instrumentation for them.
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("sci_shiny_app.R")   (model file in the same dir)
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load ----------
get_model <- function() {
  if (!exists(".SCI_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.SCI_MOD)) {
    assign(".SCI_MOD", mread_cache("sci_mrgsolve_model", project = "."),
           envir = .GlobalEnv)
  }
  .GlobalEnv$.SCI_MOD
}

LEVELS <- c("C1" = 1, "C2" = 2, "C3" = 3, "C4" = 4, "C5" = 5, "C6" = 6,
            "C7" = 7, "C8" = 8, "T1" = 9, "T2" = 10, "T3" = 11, "T4" = 12,
            "T5" = 13, "T6" = 14, "T7" = 15, "T8" = 16, "T9" = 17, "T10" = 18,
            "T11" = 19, "T12" = 20, "L1" = 21, "L2" = 22, "L3" = 23,
            "L4" = 24, "L5" = 25)

SCENARIO_LIST <- c(
  "1. Natural history (no decompression, no rehab)",
  "2. Standard care (decompression 24 h + rehab)",
  "3. Early decompression < 12 h",
  "4. Late decompression 72 h",
  "5. Methylprednisolone within 3 h (NASCIS II)",
  "6. Methylprednisolone at 9 h (outside the window)",
  "7. Methylprednisolone 48-h infusion (NASCIS III)",
  "8. Riluzole 14 d (RISCIS)",
  "9. Minocycline 7 d",
  "10. Glibenclamide 3 d (SUR1-TRPM4)",
  "11. MAP >= 85 mmHg x 7 d",
  "12. Bundle: decompression 8 h + MAP protocol + riluzole",
  "13. Anti-Nogo-A + intensive rehab from d14",
  "14. Intensive rehab EARLY (d14-180)",
  "15. Intensive rehab LATE (d120-286, same dose)",
  "16. Symptomatic: baclofen + pregabalin from d30",
  "17. Incomplete AIS C, standard care",
  "18. AIS C + epidural stimulation from d120"
)

## ---------- Regimen builders (times in DAYS) ----------
ev_mp <- function(start_h = 3, hours = 23) {
  t0 <- start_h / 24
  ev(amt = 30, cmt = "MP_CENT", time = t0) +
    ev(amt = 5.4 * hours, rate = 5.4 * 24, cmt = "MP_CENT", time = t0)
}
ev_riluzole <- function(start_h = 12, days = 14) {
  t0 <- start_h / 24
  ev(amt = 100, cmt = "RIL_DEPOT", time = t0, ii = 0.5, addl = 1) +
    ev(amt = 50, cmt = "RIL_DEPOT", time = t0 + 1, ii = 0.5,
       addl = max(floor((days - 1) / 0.5) - 1, 0))
}
ev_minocycline <- function(start_h = 12, days = 7) {
  t0 <- start_h / 24
  ev(amt = 800, cmt = "MINO_CENT", time = t0, ii = 0.5, addl = 3) +
    ev(amt = 400, cmt = "MINO_CENT", time = t0 + 2, ii = 0.5,
       addl = max(floor((days - 2) / 0.5) - 1, 0))
}
ev_glibenclamide <- function(start_h = 6, days = 3)
  ev(amt = 1, cmt = "GLY_CENT", time = start_h / 24, ii = 0.25,
     addl = floor(days / 0.25))
ev_antinogo <- function(start_d = 14, n = 6, interval = 7)
  ev(amt = 1, cmt = "NOGO_ITH", time = start_d, ii = interval, addl = n - 1)
ev_baclofen <- function(start_d = 30, end_d = 365)
  ev(amt = 20, cmt = "BAC_CENT", time = start_d, ii = 1 / 3,
     addl = floor((end_d - start_d) * 3))
ev_pregabalin <- function(start_d = 30, end_d = 365)
  ev(amt = 150, cmt = "PGB_CENT", time = start_d, ii = 0.5,
     addl = floor((end_d - start_d) * 2))

NULL_EV <- function() ev(amt = 0, cmt = "MP_CENT", time = 0)

## ---------- Scenario -> (params, events) ----------
scenario_spec <- function(scenario, horizon) {
  std <- list(TDECOMP = 1.0, REHAB = 0.5)
  switch(
    substr(scenario, 1, 3),
    "1. " = list(par = list(TDECOMP = 999, REHAB = 0), ev = NULL_EV()),
    "2. " = list(par = std, ev = NULL_EV()),
    "3. " = list(par = list(TDECOMP = 0.5, REHAB = 0.5), ev = NULL_EV()),
    "4. " = list(par = list(TDECOMP = 3.0, REHAB = 0.5), ev = NULL_EV()),
    "5. " = list(par = std, ev = ev_mp(3)),
    "6. " = list(par = std, ev = ev_mp(9)),
    "7. " = list(par = std, ev = ev_mp(3, hours = 47)),
    "8. " = list(par = std, ev = ev_riluzole()),
    "9. " = list(par = std, ev = ev_minocycline()),
    "10." = list(par = std, ev = ev_glibenclamide()),
    "11." = list(par = c(std, list(VASO = 1)), ev = NULL_EV()),
    "12." = list(par = list(TDECOMP = 8 / 24, VASO = 1, REHAB = 0.5),
                 ev = ev_riluzole(start_h = 8)),
    "13." = list(par = list(TDECOMP = 1, REHAB = 1.0), ev = ev_antinogo()),
    "14." = list(par = list(TDECOMP = 1, REHAB = 1.0, TREHAB_START = 14,
                            TREHAB_END = 180), ev = NULL_EV()),
    "15." = list(par = list(TDECOMP = 1, REHAB = 1.0, TREHAB_START = 120,
                            TREHAB_END = 286), ev = NULL_EV()),
    "16." = list(par = std,
                 ev = ev_baclofen(end_d = horizon) + ev_pregabalin(end_d = horizon)),
    "17." = list(par = c(std, list(AIS_INIT = 3)), ev = NULL_EV()),
    "18." = list(par = c(std, list(AIS_INIT = 3, ESTIM = 1)), ev = NULL_EV()),
    list(par = std, ev = NULL_EV())
  )
}

run_sim <- function(scenario, horizon, subj, delta = 0.25) {
  spec <- scenario_spec(scenario, horizon)
  par <- modifyList(
    list(LEVEL_IDX = subj$level, AIS_INIT = subj$ais, COMP0 = subj$comp,
         SUBLESIONAL_LMN = subj$lmn, MAP_BASE = subj$map_base),
    spec$par
  )
  ## a scenario that does not itself set the AIS grade inherits the sidebar one
  get_model() %>% param(par) %>%
    mrgsim(events = spec$ev, end = horizon, delta = delta, hmax = 0.01) %>%
    as_tibble() %>% mutate(scenario = scenario)
}

THEME <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

long_plot <- function(df, vars, labels = NULL, ylab = "", logx = FALSE) {
  d <- df %>% select(time, scenario, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "var", values_to = "value")
  if (!is.null(labels)) d$var <- factor(d$var, levels = vars, labels = labels)
  p <- ggplot(d, aes(time, value, colour = var)) + geom_line(linewidth = 0.7) +
    labs(x = "days after injury", y = ylab) + THEME
  if (length(unique(df$scenario)) > 1) p <- p + facet_wrap(~scenario)
  if (logx) p <- p + scale_x_log10()
  p
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  dashboardHeader(title = "Spinal Cord Injury QSP", titleWidth = 320),
  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient & injury profile", tabName = "profile", icon = icon("user-injured")),
      menuItem("Drug PK", tabName = "pk", icon = icon("pills")),
      menuItem("Secondary cascade", tabName = "cascade", icon = icon("bolt")),
      menuItem("Neuroinflammation", tabName = "inflam", icon = icon("shield-virus")),
      menuItem("Tissue & imaging", tabName = "tissue", icon = icon("dna")),
      menuItem("Neurological recovery", tabName = "recovery", icon = icon("person-walking")),
      menuItem("Therapeutic window", tabName = "window", icon = icon("clock")),
      menuItem("Surrogate vs endpoint", tabName = "surrogate", icon = icon("chart-line")),
      menuItem("Secondary complications", tabName = "compl", icon = icon("triangle-exclamation")),
      menuItem("Scenario comparison", tabName = "compare", icon = icon("layer-group")),
      menuItem("Steroid safety", tabName = "safety", icon = icon("skull-crossbones")),
      menuItem("References", tabName = "refs", icon = icon("book"))
    ),
    hr(),
    selectInput("level", "Neurological level", choices = names(LEVELS), selected = "C5"),
    selectInput("ais", "Initial AIS grade",
                choices = c("A (complete)" = 1, "B (sensory incomplete)" = 2,
                            "C (motor incomplete, <3)" = 3,
                            "D (motor incomplete, >=3)" = 4), selected = 1),
    sliderInput("comp", "Residual cord compression (0-1)", 0, 1, 0.60, 0.05),
    sliderInput("map_base", "Pre-injury MAP (mmHg)", 70, 110, 90, 1),
    checkboxInput("lmn", "Suprasacral lesion (reflex arc intact)", TRUE),
    sliderInput("horizon", "Simulation horizon (days)", 30, 730, 365, 5),
    selectInput("scenario", "Scenario", choices = SCENARIO_LIST,
                selected = SCENARIO_LIST[2]),
    helpText("Cascade tabs default to a 14-day view; recovery tabs use the",
             "full horizon.")
  ),
  dashboardBody(
    tabItems(
      ## ---- 1. profile ----
      tabItem(
        "profile",
        fluidRow(
          valueBoxOutput("vb_motor0"), valueBoxOutput("vb_motor1"),
          valueBoxOutput("vb_lesion")
        ),
        fluidRow(
          box(title = "Where this patient sits on the drive->motor curve",
              width = 7, status = "primary", solidHeader = TRUE,
              plotOutput("p_profile_hill", height = 330)),
          box(title = "Injury summary", width = 5, status = "primary",
              solidHeader = TRUE, tableOutput("t_profile"))
        ),
        fluidRow(
          box(title = "Haemodynamics: MAP, ISP and spinal cord perfusion pressure",
              width = 12, status = "info", solidHeader = TRUE,
              plotOutput("p_profile_scpp", height = 300),
              helpText("SCPP = MAP - ISP. Neurogenic shock (lesion >= T6) lowers",
                       "MAP exactly while swelling raises ISP; the ischemia term",
                       "is a deficit relative to intact perfusion."))
        )
      ),
      ## ---- 2. PK ----
      tabItem(
        "pk",
        fluidRow(
          box(title = "Acute-phase drug exposure (first 14 days)", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_pk_acute", height = 340))
        ),
        fluidRow(
          box(title = "Chronic / symptomatic agents", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("p_pk_chronic", height = 300)),
          box(title = "Cumulative methylprednisolone exposure", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("p_pk_auc", height = 300),
              helpText("MP_AUC is a function of DOSE alone — the same integral",
                       "whether the drug is given at 1 h or 24 h."))
        )
      ),
      ## ---- 3. cascade ----
      tabItem(
        "cascade",
        fluidRow(
          box(title = "Secondary injury cascade (first 14 days)", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_cascade", height = 360),
              helpText("Glutamate peaks within hours, oxidative burden within a",
                       "day, edema at day 2-3; the cascade flux is essentially",
                       "over by 72 h. This is the window every neuroprotectant",
                       "has to hit."))
        ),
        fluidRow(
          box(title = "The compression-edema-ischemia vicious cycle", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("p_cycle", height = 300)),
          box(title = "Cascade flux (weighted death-signal input)", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("p_flux", height = 300))
        )
      ),
      ## ---- 4. inflammation ----
      tabItem(
        "inflam",
        fluidRow(
          box(title = "Neutrophils, cytokines and the M1/M2 balance", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_inflam", height = 360))
        ),
        fluidRow(
          box(title = "Debris generation and phagocytic clearance", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("p_debris", height = 300),
              helpText("Debris sustains M1 activation; M2-mediated clearance is",
                       "what allows the inflammatory phase to resolve at all.")),
          box(title = "Caspase-3 activity vs the destruction threshold", width = 6,
              status = "info", solidHeader = TRUE,
              plotOutput("p_apop", height = 300))
        )
      ),
      ## ---- 5. tissue ----
      tabItem(
        "tissue",
        fluidRow(
          box(title = "Tissue integrity over the full horizon", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_tissue", height = 350),
              helpText("AXON only ever falls. OLIG recovers, but never past the",
                       "OLIG_CAP ceiling set by how many oligodendrocyte-lineage",
                       "cells survived."))
        ),
        fluidRow(
          box(title = "Imaging surrogate: lesion volume", width = 6,
              status = "warning", solidHeader = TRUE,
              plotOutput("p_lesion", height = 300)),
          box(title = "Glial scar (GFAP, CSPG) and the plasticity window",
              width = 6, status = "warning", solidHeader = TRUE,
              plotOutput("p_scar", height = 300))
        )
      ),
      ## ---- 6. recovery ----
      tabItem(
        "recovery",
        fluidRow(
          box(title = "ISNCSCI total motor score", width = 8, status = "primary",
              solidHeader = TRUE, plotOutput("p_motor", height = 340),
              helpText("The dip at day 2-3 is reversible conduction block, not",
                       "tissue loss — which is why patients improve from their",
                       "admission exam.")),
          box(title = "Recovery decomposition at end of horizon", width = 4,
              status = "primary", solidHeader = TRUE,
              tableOutput("t_recovery"))
        ),
        fluidRow(
          box(title = "Effective descending drive and its three components",
              width = 12, status = "info", solidHeader = TRUE,
              plotOutput("p_conn", height = 300),
              helpText("CONN = AXON x myelin efficiency x conduction factor x",
                       "(1 + plasticity + stimulation gain)."))
        )
      ),
      ## ---- 7. therapeutic window ----
      tabItem(
        "window",
        fluidRow(
          box(width = 12, status = "danger", solidHeader = TRUE,
              title = "Timing dominates dose",
              helpText("Left: benefit of the NASCIS methylprednisolone regimen",
                       "as a function of start time. Right: the complication",
                       "index for the same regimen. Efficacy is an integral",
                       "against a decaying flux; toxicity is not."),
              radioButtons("win_agent", NULL, inline = TRUE,
                           choices = c("Methylprednisolone start time" = "mp",
                                       "Time to surgical decompression" = "decomp"),
                           selected = "mp"),
              selectInput("win_ais", "Baseline AIS grade for the sweep",
                          choices = c("A" = 1, "B" = 2, "C" = 3, "D" = 4),
                          selected = 3),
              actionButton("run_window", "Run sweep", icon = icon("play"),
                           class = "btn-danger"))
        ),
        fluidRow(
          box(title = "Benefit vs delay", width = 6, status = "danger",
              solidHeader = TRUE, plotOutput("p_window", height = 330)),
          box(title = "Sweep table", width = 6, status = "danger",
              solidHeader = TRUE, DTOutput("t_window"))
        )
      ),
      ## ---- 8. surrogate vs endpoint ----
      tabItem(
        "surrogate",
        fluidRow(
          box(width = 12, status = "danger", solidHeader = TRUE,
              title = "The surrogate and the endpoint are not the same variable",
              helpText("MOTOR is a steep saturating function of spared drive.",
                       "The SAME fractional gain in drive is worth a different",
                       "number of ISNCSCI points at every starting point — which",
                       "is why an intervention can cut lesion volume in every",
                       "patient and only be detectable in some of them."))
        ),
        fluidRow(
          box(title = "Drive -> motor score (with this patient marked)",
              width = 7, status = "danger", solidHeader = TRUE,
              plotOutput("p_hill", height = 340)),
          box(title = "Marginal value of a doubling of spared drive", width = 5,
              status = "danger", solidHeader = TRUE, DTOutput("t_hill"))
        ),
        fluidRow(
          box(title = "Grade interaction: identical tissue effect, different endpoint effect",
              width = 12, status = "warning", solidHeader = TRUE,
              sliderInput("tissue_gain", "Fractional reduction in tissue loss",
                          0, 0.5, 0.10, 0.01),
              plotOutput("p_grade", height = 300))
        )
      ),
      ## ---- 9. complications ----
      tabItem(
        "compl",
        fluidRow(
          box(title = "Reflex reorganization and what it produces", width = 12,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_reflex", height = 340),
              helpText("Spasticity, detrusor overactivity and dysreflexia",
                       "susceptibility all rise as sublesional reflex circuits",
                       "reorganize — i.e. they are consequences of the recovery",
                       "process, not of the original injury alone."))
        ),
        fluidRow(
          box(title = "Autonomic dysreflexia challenge", width = 6,
              status = "danger", solidHeader = TRUE,
              sliderInput("ad_day", "Bladder-distension trigger day", 3, 365, 180, 1),
              actionButton("run_ad", "Apply trigger", icon = icon("bolt")),
              plotOutput("p_ad", height = 260)),
          box(title = "Organ systems: respiratory function and sublesional bone",
              width = 6, status = "info", solidHeader = TRUE,
              plotOutput("p_organ", height = 320))
        )
      ),
      ## ---- 10. comparison ----
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Compare scenarios",
              checkboxGroupInput("cmp", NULL, choices = SCENARIO_LIST,
                                 selected = SCENARIO_LIST[c(2, 5, 6, 11, 14, 15)],
                                 inline = FALSE),
              actionButton("run_cmp", "Run comparison", icon = icon("play"),
                           class = "btn-primary"))
        ),
        fluidRow(
          box(title = "ISNCSCI motor score", width = 6, status = "primary",
              solidHeader = TRUE, plotOutput("p_cmp_motor", height = 320)),
          box(title = "Lesion volume (imaging surrogate)", width = 6,
              status = "primary", solidHeader = TRUE,
              plotOutput("p_cmp_lesion", height = 320))
        ),
        fluidRow(
          box(title = "End-of-horizon summary", width = 12, status = "primary",
              solidHeader = TRUE, DTOutput("t_cmp"))
        )
      ),
      ## ---- 11. safety ----
      tabItem(
        "safety",
        fluidRow(
          box(title = "Benefit and complication index side by side", width = 12,
              status = "danger", solidHeader = TRUE,
              plotOutput("p_safety", height = 340),
              helpText("The 24-h and 48-h NASCIS regimens are compared against",
                       "no steroid. The benefit curve saturates; the exposure",
                       "curve does not."))
        ),
        fluidRow(
          box(title = "Steroid exposure ledger", width = 12, status = "danger",
              solidHeader = TRUE, DTOutput("t_safety"))
        )
      ),
      ## ---- 12. references ----
      tabItem(
        "refs",
        box(width = 12, status = "primary", solidHeader = TRUE,
            title = "Key references",
            HTML(paste0(
              "<ul>",
              "<li>Bracken MB et al. NASCIS II. <i>N Engl J Med</i> 1990;322:1405.</li>",
              "<li>Bracken MB et al. NASCIS III. <i>JAMA</i> 1997;277:1597.</li>",
              "<li>Fehlings MG et al. STASCIS. <i>PLoS One</i> 2012;7:e32037.</li>",
              "<li>Badhiwala JH et al. Decompression timing (pooled IPD). ",
              "<i>Lancet Neurol</i> 2021;20:117.</li>",
              "<li>Casha S et al. Minocycline phase 2. <i>Brain</i> 2012;135:1224.</li>",
              "<li>Grossman RG et al. NACTN riluzole. <i>J Neurotrauma</i> 2014;31:239.</li>",
              "<li>Squair JW et al. Spinal cord perfusion pressure. ",
              "<i>Neurology</i> 2017;89:1660.</li>",
              "<li>Angeli CA et al. Epidural stimulation. ",
              "<i>N Engl J Med</i> 2018;379:1244.</li>",
              "<li>Wagner FB et al. Targeted neurotechnology. ",
              "<i>Nature</i> 2018;563:65.</li>",
              "<li>Fawcett JW et al. ICCP clinical-trial guidelines. ",
              "<i>Spinal Cord</i> 2007;45:190.</li>",
              "<li>Kirshblum SC et al. ISNCSCI 2011 revision. ",
              "<i>J Spinal Cord Med</i> 2011;34:535.</li>",
              "<li>Ahuja CS et al. Traumatic spinal cord injury. ",
              "<i>Nat Rev Dis Primers</i> 2017;3:17018.</li>",
              "</ul>",
              "<p>The full annotated bibliography is in ",
              "<code>sci_references.md</code>; the mechanistic map is ",
              "<code>sci_qsp_model.dot</code>; every quantitative claim in the ",
              "directory README is reproduced by ",
              "<code>sci_reference_model.py</code>.</p>",
              "<p><b>Research and education only. Not for clinical use.</b></p>"
            )))
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  subj <- reactive(list(
    level = as.numeric(LEVELS[[input$level]]),
    ais = as.numeric(input$ais),
    comp = input$comp,
    lmn = if (isTRUE(input$lmn)) 1 else 0.4,
    map_base = input$map_base
  ))

  sim <- reactive({
    run_sim(input$scenario, input$horizon, subj())
  })

  acute <- reactive(sim() %>% filter(time <= min(14, input$horizon)))

  ## ---- 1. profile ----
  output$vb_motor0 <- renderValueBox({
    v <- sim()$ISNCSCI_motor[1]
    valueBox(round(v, 1), "ISNCSCI motor at presentation",
             icon = icon("stethoscope"), color = "orange")
  })
  output$vb_motor1 <- renderValueBox({
    d <- sim()
    valueBox(round(tail(d$ISNCSCI_motor, 1), 1),
             sprintf("ISNCSCI motor at day %d (%+.1f)", input$horizon,
                     tail(d$ISNCSCI_motor, 1) - d$ISNCSCI_motor[1]),
             icon = icon("person-walking"), color = "green")
  })
  output$vb_lesion <- renderValueBox({
    valueBox(sprintf("%.2f mL", tail(sim()$Lesion_volume_mL, 1)),
             "Lesion volume (imaging surrogate)", icon = icon("magnifying-glass"),
             color = "red")
  })

  output$t_profile <- renderTable({
    d <- sim()
    data.frame(
      Quantity = c("Neurological level", "Initial AIS grade",
                   "Motor points preserved above lesion",
                   "Spared descending axons (t0 -> end)",
                   "Effective descending drive (end)",
                   "Peak ischemic burden", "Peak cord edema",
                   "Vital capacity (end)", "Sublesional BMD (end)"),
      Value = c(input$level,
                c("A", "B", "C", "D")[as.numeric(input$ais)],
                sprintf("%.0f", d$Motor_above_lesion[1]),
                sprintf("%.3f -> %.3f", d$Axon_sparing[1], tail(d$Axon_sparing, 1)),
                sprintf("%.3f", tail(d$Descending_drive, 1)),
                sprintf("%.2f", max(d$Ischemic_burden)),
                sprintf("%.2f", max(d$EDEMA)),
                sprintf("%.0f %% predicted", tail(d$Vital_capacity, 1)),
                sprintf("%.2f", tail(d$BMD_fraction, 1)))
    )
  })

  hill_curve <- function(th50 = 0.22, hill = 2.2, above = 10) {
    x <- seq(0.001, 0.9, length.out = 400)
    data.frame(conn = x,
               motor = above + (100 - above) * x^hill / (th50^hill + x^hill))
  }

  output$p_profile_hill <- renderPlot({
    d <- sim()
    above <- d$Motor_above_lesion[1]
    cur <- hill_curve(above = above)
    pt <- data.frame(conn = c(d$Descending_drive[1], tail(d$Descending_drive, 1)),
                     motor = c(d$ISNCSCI_motor[1], tail(d$ISNCSCI_motor, 1)),
                     lab = c("presentation", "end of horizon"))
    ggplot(cur, aes(conn, motor)) + geom_line(linewidth = 0.9) +
      geom_point(data = pt, aes(colour = lab), size = 4) +
      geom_vline(xintercept = 0.22, linetype = "dashed", colour = "grey50") +
      annotate("text", x = 0.23, y = 5, hjust = 0, size = 3.2,
               label = "TH50 = 0.22") +
      labs(x = "effective descending drive (CONN)",
           y = "ISNCSCI total motor score") + THEME
  })

  output$p_profile_scpp <- renderPlot({
    long_plot(acute(), c("MAP", "ISP", "SCPP"),
              c("MAP (mmHg)", "ISP (mmHg)", "SCPP (mmHg)"), "mmHg") +
      geom_hline(yintercept = 85, linetype = "dotted")
  })

  ## ---- 2. PK ----
  output$p_pk_acute <- renderPlot({
    long_plot(acute(), c("MP_CENT", "RIL_CENT", "MINO_CENT", "GLY_CENT"),
              c("methylprednisolone", "riluzole", "minocycline", "glibenclamide"),
              "amount (au)")
  })
  output$p_pk_chronic <- renderPlot({
    long_plot(sim(), c("NOGO_ITH", "BAC_CENT", "PGB_CENT"),
              c("anti-Nogo-A (intrathecal)", "baclofen", "pregabalin"),
              "amount (au)")
  })
  output$p_pk_auc <- renderPlot({
    ggplot(acute(), aes(time, MP_AUC)) + geom_line(linewidth = 0.8) +
      labs(x = "days after injury", y = "cumulative MP exposure (au*day)") + THEME
  })

  ## ---- 3. cascade ----
  output$p_cascade <- renderPlot({
    long_plot(acute(), c("Ischemic_burden", "EDEMA", "GLU", "CAI", "ROS"),
              c("ischemic burden", "cord edema", "glutamate", "Ca2+ overload",
                "oxidative burden"), "normalized au")
  })
  output$p_cycle <- renderPlot({
    long_plot(acute(), c("ISP", "Ischemic_burden", "EDEMA"),
              c("ISP (mmHg)", "ischemia", "edema"), "value")
  })
  output$p_flux <- renderPlot({
    ggplot(acute(), aes(time, Cascade_flux)) + geom_line(linewidth = 0.8) +
      labs(x = "days after injury", y = "cascade flux (au/day)") + THEME
  })

  ## ---- 4. inflammation ----
  output$p_inflam <- renderPlot({
    d <- sim() %>% filter(time <= min(60, input$horizon))
    long_plot(d, c("NEUT", "CYTO", "M1", "M2"),
              c("neutrophils", "cytokines (au)", "M1", "M2"), "au")
  })
  output$p_debris <- renderPlot({
    d <- sim() %>% filter(time <= min(60, input$horizon))
    long_plot(d, c("DEBRIS", "M2"), c("debris", "M2 (clearance)"), "au")
  })
  output$p_apop <- renderPlot({
    d <- sim() %>% filter(time <= min(30, input$horizon))
    ggplot(d, aes(time, APOP)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.05, linetype = "dashed", colour = "red") +
      annotate("text", x = max(d$time) * 0.6, y = 0.09, size = 3.2,
               label = "APOP_THR: below this, no tissue is destroyed") +
      labs(x = "days after injury", y = "caspase-3 activity (au)") + THEME
  })

  ## ---- 5. tissue ----
  output$p_tissue <- renderPlot({
    long_plot(sim(), c("AXON", "OLIG", "OLIG_CAP", "NEURON"),
              c("spared axons", "myelin integrity", "remyelination ceiling",
                "perilesional neurons"), "fraction")
  })
  output$p_lesion <- renderPlot({
    ggplot(sim(), aes(time, Lesion_volume_mL)) + geom_line(linewidth = 0.8) +
      labs(x = "days after injury", y = "lesion volume (mL-equivalent)") + THEME
  })
  output$p_scar <- renderPlot({
    long_plot(sim(), c("GFAP", "Scar_CSPG", "Critical_period", "Plasticity_pool"),
              c("GFAP (gliosis)", "CSPG (scar)", "critical period (open)",
                "plasticity pool"), "au / fraction")
  })

  ## ---- 6. recovery ----
  output$p_motor <- renderPlot({
    ggplot(sim(), aes(time, ISNCSCI_motor)) + geom_line(linewidth = 0.9) +
      labs(x = "days after injury", y = "ISNCSCI total motor score (0-100)") +
      THEME
  })
  output$t_recovery <- renderTable({
    d <- sim()
    data.frame(
      Component = c("motor at presentation", "motor at end", "change",
                    "axon sparing (end)", "myelin efficiency (end)",
                    "conduction factor (end)", "plasticity pool (end)",
                    "critical period remaining (end)"),
      Value = c(sprintf("%.1f", d$ISNCSCI_motor[1]),
                sprintf("%.1f", tail(d$ISNCSCI_motor, 1)),
                sprintf("%+.1f", tail(d$ISNCSCI_motor, 1) - d$ISNCSCI_motor[1]),
                sprintf("%.3f", tail(d$Axon_sparing, 1)),
                sprintf("%.3f", tail(d$Myelin_efficiency, 1)),
                sprintf("%.3f", tail(d$Conduction_factor, 1)),
                sprintf("%.3f", tail(d$Plasticity_pool, 1)),
                sprintf("%.4f", tail(d$Critical_period, 1)))
    )
  })
  output$p_conn <- renderPlot({
    long_plot(sim(), c("Descending_drive", "AXON", "Myelin_efficiency",
                       "Conduction_factor", "Plasticity_pool"),
              c("CONN (effective drive)", "spared axons", "myelin efficiency",
                "conduction factor", "plasticity pool"), "fraction")
  })

  ## ---- 7. therapeutic window ----
  window_res <- eventReactive(input$run_window, {
    ais <- as.numeric(input$win_ais)
    s <- subj(); s$ais <- ais
    ref <- get_model() %>%
      param(LEVEL_IDX = s$level, AIS_INIT = ais, COMP0 = s$comp, REHAB = 0.5) %>%
      mrgsim(end = 365, delta = 1, hmax = 0.01) %>% as_tibble()
    m_ref <- tail(ref$ISNCSCI_motor, 1)
    les_ref <- tail(ref$Lesion_volume_mL, 1)

    if (input$win_agent == "mp") {
      hours <- c(1, 3, 6, 9, 12, 24)
      out <- lapply(hours, function(h) {
        o <- get_model() %>%
          param(LEVEL_IDX = s$level, AIS_INIT = ais, COMP0 = s$comp, REHAB = 0.5) %>%
          mrgsim(events = ev_mp(h), end = 365, delta = 1, hmax = 0.01) %>%
          as_tibble()
        data.frame(delay_h = h,
                   motor = tail(o$ISNCSCI_motor, 1),
                   d_motor = tail(o$ISNCSCI_motor, 1) - m_ref,
                   lesion_mL = tail(o$Lesion_volume_mL, 1),
                   d_lesion_pct = 100 * (tail(o$Lesion_volume_mL, 1) - les_ref) / les_ref,
                   complication_index = tail(o$MP_complication_index, 1))
      })
    } else {
      hours <- c(4, 8, 12, 24, 48, 72)
      out <- lapply(hours, function(h) {
        o <- get_model() %>%
          param(LEVEL_IDX = s$level, AIS_INIT = ais, COMP0 = s$comp,
                REHAB = 0.5, TDECOMP = h / 24) %>%
          mrgsim(end = 365, delta = 0.25, hmax = 0.01) %>% as_tibble()
        data.frame(delay_h = h,
                   motor = tail(o$ISNCSCI_motor, 1),
                   d_motor = tail(o$ISNCSCI_motor, 1) - m_ref,
                   lesion_mL = tail(o$Lesion_volume_mL, 1),
                   d_lesion_pct = 100 * (tail(o$Lesion_volume_mL, 1) - les_ref) / les_ref,
                   complication_index = NA_real_)
      })
    }
    bind_rows(out)
  })

  output$p_window <- renderPlot({
    d <- window_res()
    ggplot(d, aes(delay_h, d_motor)) +
      geom_line(linewidth = 0.9) + geom_point(size = 3) +
      { if (all(!is.na(d$complication_index)))
          geom_line(aes(y = complication_index * max(abs(d$d_motor)) /
                          max(d$complication_index)),
                    linetype = "dashed", colour = "red", linewidth = 0.9) } +
      labs(x = "delay from injury to intervention (h)",
           y = "change in ISNCSCI motor at 1 year",
           caption = paste("solid = benefit;",
                           "dashed red (methylprednisolone only) = complication",
                           "index, rescaled - it does not fall with delay")) +
      THEME
  })
  output$t_window <- renderDT({
    datatable(window_res() %>% mutate(across(where(is.numeric), ~round(., 3))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ---- 8. surrogate vs endpoint ----
  output$p_hill <- renderPlot({
    d <- sim()
    above <- d$Motor_above_lesion[1]
    cur <- hill_curve(above = above)
    c0 <- tail(d$Descending_drive, 1)
    seg <- data.frame(x = c0, xend = min(2 * c0, 0.9))
    seg$y <- above + (100 - above) * seg$x^2.2 / (0.22^2.2 + seg$x^2.2)
    seg$yend <- above + (100 - above) * seg$xend^2.2 / (0.22^2.2 + seg$xend^2.2)
    ggplot(cur, aes(conn, motor)) + geom_line(linewidth = 0.9) +
      geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                   colour = "red", linewidth = 1.4,
                   arrow = arrow(length = unit(0.18, "cm"))) +
      geom_point(data = seg, aes(x = x, y = y), size = 4, colour = "red") +
      labs(x = "effective descending drive (CONN)",
           y = "ISNCSCI total motor score",
           caption = sprintf(paste("red arrow: doubling this patient's spared",
                                   "drive (%.3f -> %.3f) is worth %+.1f motor points"),
                             seg$x, seg$xend, seg$yend - seg$y)) + THEME
  })
  output$t_hill <- renderDT({
    above <- sim()$Motor_above_lesion[1]
    f <- function(c) above + (100 - above) * c^2.2 / (0.22^2.2 + c^2.2)
    c0 <- c(0.03, 0.06, 0.11, 0.15, 0.22, 0.30, 0.45)
    datatable(data.frame(CONN = c0, motor = round(f(c0), 1),
                         motor_if_doubled = round(f(2 * c0), 1),
                         gain = round(f(2 * c0) - f(c0), 1),
                         detectable = ifelse(f(2 * c0) - f(c0) >= 5,
                                             "yes (>5 pts)", "no (< noise)")),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
  output$p_grade <- renderPlot({
    g <- input$tissue_gain
    grades <- data.frame(ais = c("A", "B", "C", "D"),
                         axon0 = c(0.06, 0.14, 0.30, 0.55),
                         above = c(10, 10, 10, 10))
    ## secondary loss of ~28% of spared axons, reduced by the tissue gain
    grades$axon_ctl <- grades$axon0 * (1 - 0.28)
    grades$axon_trt <- grades$axon0 * (1 - 0.28 * (1 - g))
    f <- function(a, above) {
      c <- pmin(0.999, a * 0.82 * 1.4)
      above + (100 - above) * c^2.2 / (0.22^2.2 + c^2.2)
    }
    grades$control <- f(grades$axon_ctl, grades$above)
    grades$treated <- f(grades$axon_trt, grades$above)
    d <- grades %>% mutate(delta = treated - control)
    ggplot(d, aes(ais, delta, fill = ais)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "red") +
      annotate("text", x = 1.2, y = 5.6, size = 3.2, colour = "red",
               label = "~ISNCSCI measurement noise") +
      labs(x = "baseline AIS grade", y = "gain in ISNCSCI motor points",
           caption = sprintf(paste("identical %.0f%% reduction in tissue loss in",
                                   "every grade"), 100 * g)) +
      THEME + theme(legend.position = "none")
  })

  ## ---- 9. complications ----
  output$p_reflex <- renderPlot({
    long_plot(sim(), c("REFLEX", "Ashworth", "Pain_NRS", "Bladder_index"),
              c("reflex reorganization (0-1)", "modified Ashworth (0-4)",
                "pain NRS (0-10)", "detrusor overactivity (0-1)"), "value")
  })
  ad_res <- eventReactive(input$run_ad, {
    s <- subj()
    d <- input$ad_day
    get_model() %>%
      param(LEVEL_IDX = s$level, AIS_INIT = s$ais, COMP0 = s$comp, REHAB = 0.5) %>%
      mrgsim(events = ev(amt = 1, cmt = "AD_TRIG", time = d),
             end = d + 0.5, delta = 0.002, hmax = 0.002) %>%
      as_tibble() %>% filter(time >= d - 0.05)
  })
  output$p_ad <- renderPlot({
    d <- ad_res()
    ggplot(d, aes((time - min(time)) * 24 * 60, AD_peak_SBP)) +
      geom_line(linewidth = 0.9, colour = "red") +
      geom_hline(yintercept = 150, linetype = "dashed") +
      labs(x = "minutes after the noxious stimulus", y = "systolic BP (mmHg)",
           caption = sprintf("peak SBP %.0f mmHg at reflex maturity %.2f",
                             max(d$AD_peak_SBP), tail(d$REFLEX, 1))) + THEME
  })
  output$p_organ <- renderPlot({
    long_plot(sim(), c("Vital_capacity", "BMD_fraction", "ATRO"),
              c("vital capacity (% predicted)", "BMD (fraction)",
                "muscle atrophy (fraction)"), "value") +
      facet_wrap(~var, scales = "free_y", ncol = 1)
  })

  ## ---- 10. comparison ----
  cmp_res <- eventReactive(input$run_cmp, {
    req(length(input$cmp) > 0)
    bind_rows(lapply(input$cmp, run_sim, horizon = input$horizon,
                     subj = subj(), delta = 1))
  })
  output$p_cmp_motor <- renderPlot({
    ggplot(cmp_res(), aes(time, ISNCSCI_motor, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      labs(x = "days after injury", y = "ISNCSCI total motor") + THEME +
      theme(legend.text = element_text(size = 7))
  })
  output$p_cmp_lesion <- renderPlot({
    ggplot(cmp_res(), aes(time, Lesion_volume_mL, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      labs(x = "days after injury", y = "lesion volume (mL-equivalent)") + THEME +
      theme(legend.text = element_text(size = 7))
  })
  output$t_cmp <- renderDT({
    d <- cmp_res() %>% group_by(scenario) %>%
      summarise(motor_start = round(first(ISNCSCI_motor), 1),
                motor_end = round(last(ISNCSCI_motor), 1),
                change = round(last(ISNCSCI_motor) - first(ISNCSCI_motor), 1),
                lesion_mL = round(last(Lesion_volume_mL), 3),
                axon = round(last(Axon_sparing), 3),
                peak_ischemia = round(max(Ischemic_burden), 3),
                Ashworth = round(last(Ashworth), 2),
                pain_NRS = round(last(Pain_NRS), 2),
                MP_index = round(last(MP_complication_index), 3),
                .groups = "drop")
    datatable(d, options = list(dom = "t", pageLength = 20, scrollX = TRUE),
              rownames = FALSE)
  })

  ## ---- 11. safety ----
  safety_res <- reactive({
    s <- subj()
    base <- function(par, e) {
      get_model() %>% param(par) %>%
        mrgsim(events = e, end = 365, delta = 1, hmax = 0.01) %>% as_tibble()
    }
    par <- list(LEVEL_IDX = s$level, AIS_INIT = s$ais, COMP0 = s$comp, REHAB = 0.5)
    ref <- base(par, NULL_EV())
    arms <- list("no steroid" = NULL_EV(),
                 "MPSS 24 h from 3 h" = ev_mp(3, 23),
                 "MPSS 24 h from 9 h" = ev_mp(9, 23),
                 "MPSS 48 h from 3 h" = ev_mp(3, 47))
    bind_rows(lapply(names(arms), function(nm) {
      o <- base(par, arms[[nm]])
      data.frame(arm = nm,
                 motor = tail(o$ISNCSCI_motor, 1),
                 d_motor = tail(o$ISNCSCI_motor, 1) - tail(ref$ISNCSCI_motor, 1),
                 MP_AUC = tail(o$MP_AUC, 1),
                 complication_index = tail(o$MP_complication_index, 1))
    }))
  })
  output$p_safety <- renderPlot({
    d <- safety_res() %>%
      pivot_longer(c(d_motor, complication_index), names_to = "metric")
    ggplot(d, aes(arm, value, fill = metric)) +
      geom_col(position = "dodge") +
      labs(x = NULL, y = "motor-point gain / complication index (0-1)") +
      THEME + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })
  output$t_safety <- renderDT({
    datatable(safety_res() %>% mutate(across(where(is.numeric), ~round(., 3))),
              options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
