## ============================================================================
## MOG Antibody-Associated Disease (MOGAD) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 12 tabs: Patient & attack · Drug PK / target engagement · Antibody source ·
##          Barrier & CNS effectors · Tissue (myelin vs axon) ·
##          Clinical endpoints · Relapse hazard · Biomarkers ·
##          Scenario comparison · Taper explorer · NMOSD comparator · Notes
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run from this directory:  shiny::runApp("mogad_shiny_app.R")
##
## The model file `mogad_mrgsolve_model.R` must be in the same directory.
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

## ---------- Lazy model load -------------------------------------------------
get_model <- function() {
  if (!exists(".MOGAD_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.MOGAD_MOD)) {
    assign(".MOGAD_MOD",
           mread_cache("mogad", project = ".", file = "mogad_mrgsolve_model.R"),
           envir = .GlobalEnv)
  }
  .GlobalEnv$.MOGAD_MOD
}

## ---------- Dosing conventions ---------------------------------------------
## MP_C   1 g IV methylprednisolone = amt 1250 (prednisone-equivalent mg)
## IVIG_C 1 g/kg IVIG               = amt 12 (g/L increment)
## PLEX_D one plasma exchange       = amt 1
## RTX_C / TCZ_C / FCRN_C / C5I_C   = amt 1 per administration
## TRIG   one attack trigger        = amt 1
AMT_IVMP <- 1250
AMT_IVIG <- 12

ev_attack <- function(t = 0) ev(amt = 1, cmt = "TRIG", time = t)

ev_ivmp <- function(day0 = 4, n = 5) {
  if (n <= 0) return(NULL)
  ev(amt = AMT_IVMP, cmt = "MP_C", time = day0, ii = 1, addl = n - 1)
}

ev_plex <- function(day0 = 6, n = 5, spacing = 2) {
  if (n <= 0) return(NULL)
  ev(amt = 1, cmt = "PLEX_D", time = day0, ii = spacing, addl = n - 1)
}

ev_ivig_acute <- function(day0 = 5) {
  ev(amt = AMT_IVIG, cmt = "IVIG_C", time = day0, ii = 1, addl = 1)
}

ev_repeat <- function(cmt, day0, ii, end, amt = 1) {
  if (day0 > end) return(NULL)
  ev(amt = amt, cmt = cmt, time = day0, ii = ii,
     addl = max(floor((end - day0) / ii), 0))
}

ev_rtx <- function(day0 = 28, end = 365) {
  induction <- ev(amt = 1, cmt = "RTX_C", time = day0, ii = 14, addl = 1)
  if (end > day0 + 182) {
    c(induction, ev_repeat("RTX_C", day0 + 182, 182, end))
  } else {
    induction
  }
}

## Combine possibly-NULL event objects
ev_bind <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(ev(amt = 0, cmt = "TRIG", time = 0))
  Reduce(function(a, b) c(a, b), parts)
}

## ---------- Named maintenance regimens -------------------------------------
MAINT <- c("None (steroid wean only)",
           "Maintenance IVIG 1 g/kg q4wk",
           "Maintenance IVIG 0.4 g/kg q8wk (under-dosed)",
           "Rituximab (1 g x2, then q6mo)",
           "IL-6R blockade (tocilizumab / satralizumab q4wk)",
           "Mycophenolate / azathioprine (continuous)",
           "FcRn inhibitor (weekly SC)",
           "C5 inhibitor (q2wk)",
           "IVIG + IL-6R blockade (combination)")

maint_events <- function(which, start, end) {
  switch(which,
    "None (steroid wean only)" = NULL,
    "Maintenance IVIG 1 g/kg q4wk" =
      ev_repeat("IVIG_C", start, 28, end, AMT_IVIG),
    "Maintenance IVIG 0.4 g/kg q8wk (under-dosed)" =
      ev_repeat("IVIG_C", start, 56, end, 0.4 * AMT_IVIG),
    "Rituximab (1 g x2, then q6mo)" = ev_rtx(start, end),
    "IL-6R blockade (tocilizumab / satralizumab q4wk)" =
      ev_repeat("TCZ_C", start, 28, end),
    "Mycophenolate / azathioprine (continuous)" = NULL,  # via MMF_RATE param
    "FcRn inhibitor (weekly SC)" = ev_repeat("FCRN_C", start, 7, end),
    "C5 inhibitor (q2wk)" = ev_repeat("C5I_C", start, 14, end),
    "IVIG + IL-6R blockade (combination)" =
      ev_bind(ev_repeat("IVIG_C", start, 28, end, AMT_IVIG),
              ev_repeat("TCZ_C", start, 28, end)),
    NULL)
}

maint_params <- function(which) {
  if (identical(which, "Mycophenolate / azathioprine (continuous)")) {
    list(MMF_RATE = 0.42)
  } else {
    list(MMF_RATE = 0)
  }
}

## ---------- Simulation driver ----------------------------------------------
run_sim <- function(pars, events, horizon) {
  get_model() %>%
    param(pars) %>%
    ev(events) %>%
    mrgsim(end = horizon, delta = 0.5) %>%
    as_tibble()
}

build_case <- function(inp, horizon, maint = NULL, extra_par = list()) {
  maint <- if (is.null(maint)) inp$maint else maint
  pars <- c(list(SEV = inp$sev, SITE_ON = inp$site_on, SITE_SC = inp$site_sc,
                 PRED0 = inp$pred0, TAPER_DAYS = inp$taper,
                 PRED_START = inp$ivmp_day + 5,
                 FRAC_ESC = inp$frac_esc, AQP4_MODE = 0),
            maint_params(maint), extra_par)
  evs <- ev_bind(ev_attack(0),
                 ev_ivmp(inp$ivmp_day, inp$ivmp_n),
                 if (inp$use_plex) ev_plex(inp$ivmp_day + 2, inp$plex_n, 2) else NULL,
                 if (inp$use_ivig_acute) ev_ivig_acute(inp$ivmp_day + 1) else NULL,
                 maint_events(maint, inp$maint_start, horizon),
                 if (inp$trigger2 > 0 && inp$trigger2 < horizon)
                   ev_attack(inp$trigger2) else NULL)
  run_sim(pars, evs, horizon)
}

THEME <- theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

lineplot <- function(df, ycol, title, ylab, logy = FALSE) {
  p <- ggplot(df, aes(time, .data[[ycol]])) +
    geom_line(linewidth = 0.9, colour = "#2b6cb0") +
    labs(title = title, x = "Time (days)", y = ylab) + THEME
  if (logy) p <- p + scale_y_log10()
  p
}

multiplot <- function(df, cols, labels, title, ylab) {
  long <- df %>% select(time, all_of(cols)) %>%
    pivot_longer(-time, names_to = "series", values_to = "value") %>%
    mutate(series = factor(series, levels = cols, labels = labels))
  ggplot(long, aes(time, value, colour = series)) +
    geom_line(linewidth = 0.9) +
    labs(title = title, x = "Time (days)", y = ylab) + THEME
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "MOGAD QSP Model", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("1 · Patient & attack",   tabName = "patient",  icon = icon("user")),
      menuItem("2 · Drug PK / target",   tabName = "pk",       icon = icon("pills")),
      menuItem("3 · Antibody source",    tabName = "source",   icon = icon("dna")),
      menuItem("4 · Barrier & effectors", tabName = "effector", icon = icon("shield-halved")),
      menuItem("5 · Tissue: myelin/axon", tabName = "tissue",  icon = icon("bolt")),
      menuItem("6 · Clinical endpoints", tabName = "endpoint", icon = icon("eye")),
      menuItem("7 · Relapse hazard",     tabName = "hazard",   icon = icon("chart-line")),
      menuItem("8 · Biomarkers",         tabName = "biomark",  icon = icon("vial")),
      menuItem("9 · Scenario comparison", tabName = "compare", icon = icon("layer-group")),
      menuItem("10 · Taper explorer",    tabName = "taper",    icon = icon("stairs")),
      menuItem("11 · NMOSD comparator",  tabName = "nmosd",    icon = icon("code-compare")),
      menuItem("12 · Notes & references", tabName = "notes",   icon = icon("book")),
      hr(),
      sliderInput("horizon", "Simulation horizon (days)", 90, 730, 365, 5),
      h5(HTML("&nbsp;<b>Attack</b>")),
      sliderInput("sev", "Attack drive (SEV)", 0.4, 1.6, 1.0, 0.05),
      sliderInput("site_on", "Optic-nerve weight", 0, 1.2, 1.0, 0.05),
      sliderInput("site_sc", "Spinal-cord weight", 0, 1.2, 0.45, 0.05),
      sliderInput("trigger2", "2nd trigger (day; 0 = none)", 0, 700, 0, 10),
      h5(HTML("&nbsp;<b>Acute treatment</b>")),
      sliderInput("ivmp_day", "IVMP start (treatment delay, d)", 0, 21, 4, 1),
      sliderInput("ivmp_n", "IVMP 1 g doses", 0, 7, 5, 1),
      checkboxInput("use_plex", "Add plasma exchange", FALSE),
      sliderInput("plex_n", "Number of exchanges", 3, 7, 5, 1),
      checkboxInput("use_ivig_acute", "Add acute IVIG 2 g/kg", FALSE),
      h5(HTML("&nbsp;<b>Oral steroid wean</b>")),
      sliderInput("pred0", "Starting prednisone (mg/d)", 0, 80, 60, 5),
      sliderInput("taper", "Taper duration (days)", 7, 240, 90, 7),
      h5(HTML("&nbsp;<b>Maintenance</b>")),
      selectInput("maint", NULL, choices = MAINT, selected = MAINT[2]),
      sliderInput("maint_start", "Maintenance start (day)", 0, 120, 28, 7),
      h5(HTML("&nbsp;<b>Structural assumption</b>")),
      sliderInput("frac_esc", "FRAC_ESC (CD20-negative share)", 0, 0.95, 0.60, 0.05)
    )
  ),
  dashboardBody(
    tags$style(HTML(".content-wrapper {background-color:#f7f9fb;}")),
    tabItems(
      ## ---- 1 patient ------------------------------------------------------
      tabItem("patient",
        fluidRow(
          box(width = 8, title = "Attack course — what the current settings produce",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_overview", height = 380)),
          box(width = 4, title = "Key readouts", status = "primary",
              solidHeader = TRUE, tableOutput("t_summary"))),
        fluidRow(
          box(width = 12, title = "Model in one paragraph", status = "info",
              solidHeader = TRUE,
              HTML("<p>Class-switched <b>IgG1 against MOG</b> is produced by a
                    <b>short-lived, CD20-negative plasmablast</b> pool. It reaches
                    the CNS only through a permeabilised blood-brain barrier, where
                    it opsonises the <b>outermost myelin lamella</b> and destroys
                    myelin by complement-dependent lysis and antibody-dependent
                    cellular phagocytosis. The <b>oligodendrocyte cell body
                    survives</b>, so remyelination proceeds and the recovery
                    ceiling is high; the astrocyte/AQP4 compartment is untouched,
                    so serum GFAP stays low. Relapse hazard is a
                    <b>threshold function of functionally active antibody</b>,
                    which is why the ranking of therapies does not follow the
                    ranking of titre reduction.</p>")))),

      ## ---- 2 PK -----------------------------------------------------------
      tabItem("pk",
        fluidRow(
          box(width = 6, title = "Glucocorticoid receptor occupancy",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_gc", height = 300)),
          box(width = 6, title = "Cumulative steroid dose and bone density",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_ster", height = 300))),
        fluidRow(
          box(width = 6, title = "Target engagement of maintenance agents",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_engage", height = 300)),
          box(width = 6, title = "Total serum IgG (IVIG raises it; RTX and FcRn blockade lower it)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_igg", height = 300)))),

      ## ---- 3 source -------------------------------------------------------
      tabItem("source",
        fluidRow(
          box(width = 6, title = "Memory B cells vs plasmablasts (commitment 1)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_bcell", height = 320)),
          box(width = 6, title = "MOG-IgG titre: measured vs functionally active",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_titer", height = 320))),
        fluidRow(
          box(width = 12, status = "warning", solidHeader = TRUE,
              title = "Why the FRAC_ESC slider matters",
              HTML("<p>FRAC_ESC is the share of plasmablast generation that is
                    <b>CD20-negative</b> and therefore structurally out of reach of
                    rituximab. Set it to 0 and rituximab becomes a highly effective
                    drug in this model; set it to 0.6 (the fitted default) and
                    rituximab reduces the titre while barely moving the relapse
                    hazard. It is a <b>fitted, not measured</b> parameter — the
                    single assumption the source-versus-clearance conclusion rests
                    on.</p>")))),

      ## ---- 4 effectors ----------------------------------------------------
      tabItem("effector",
        fluidRow(
          box(width = 6, title = "Blood-brain barrier permeability (the gate)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_bbb", height = 300)),
          box(width = 6, title = "CNS MOG-IgG, complement and myeloid activation",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_eff", height = 300))),
        fluidRow(
          box(width = 6, title = "IL-6 (measured rises on IL-6R blockade)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_il6", height = 300)),
          box(width = 6, title = "Lesion oedema", status = "primary",
              solidHeader = TRUE, plotOutput("p_edema", height = 300)))),

      ## ---- 5 tissue -------------------------------------------------------
      tabItem("tissue",
        fluidRow(
          box(width = 8, title = "Myelin is lost and recovered; axons are not (commitment 2)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_tissue", height = 360)),
          box(width = 4, title = "Surviving oligodendrocytes set the recovery ceiling",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_ol", height = 360))),
        fluidRow(
          box(width = 12, title = "Peripapillary RNFL thickness (structural loss)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_rnfl", height = 260)))),

      ## ---- 6 endpoints ----------------------------------------------------
      tabItem("endpoint",
        fluidRow(
          box(width = 6, title = "Visual acuity (logMAR — lower is better)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_va", height = 320)),
          box(width = 6, title = "EDSS", status = "primary", solidHeader = TRUE,
              plotOutput("p_edss", height = 320))),
        fluidRow(
          box(width = 12, title = "Endpoint summary", status = "primary",
              solidHeader = TRUE, tableOutput("t_endpoint")))),

      ## ---- 7 hazard -------------------------------------------------------
      tabItem("hazard",
        fluidRow(
          box(width = 6, title = "Instantaneous annualised relapse rate",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_arr", height = 320)),
          box(width = 6, title = "Cumulative probability of relapse",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_prel", height = 320))),
        fluidRow(
          box(width = 12, status = "warning", solidHeader = TRUE,
              title = "Steroid dependency is emergent, not imposed",
              HTML("<p>Nothing in the model says &ldquo;relapse when steroids
                    stop&rdquo;. The hazard depends on antibody, and the antibody
                    depends on a plasmablast pool with a 7-day half-life that
                    glucocorticoids suppress <i>reversibly</i>. When the wean ends,
                    plasmablasts recover within a couple of weeks and the titre
                    follows over roughly a month — so the hazard climbs back on its
                    own. Move the taper slider and watch <i>when</i>, not
                    <i>whether</i>, that happens.</p>")))),

      ## ---- 8 biomarkers ---------------------------------------------------
      tabItem("biomark",
        fluidRow(
          box(width = 6, title = "Serum neurofilament light (axonal injury)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_nfl", height = 300)),
          box(width = 6, title = "Serum GFAP (stays LOW — this is not astrocytopathy)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_gfap", height = 300))),
        fluidRow(
          box(width = 6, title = "CSF white cell count", status = "primary",
              solidHeader = TRUE, plotOutput("p_csf", height = 300)),
          box(width = 6, title = "MOG-IgG titre (live cell-based assay dilution)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_titer2", height = 300)))),

      ## ---- 9 comparison ---------------------------------------------------
      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Maintenance strategies under the current attack settings",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_cmp", height = 420))),
        fluidRow(
          box(width = 12, title = "Titre reduction does not predict relapse protection",
              status = "primary", solidHeader = TRUE,
              DT::dataTableOutput("t_cmp")))),

      ## ---- 10 taper -------------------------------------------------------
      tabItem("taper",
        fluidRow(
          box(width = 12, title = "Matched cumulative dose, different taper lengths",
              status = "primary", solidHeader = TRUE,
              sliderInput("total_mg", "Total prednisone-equivalent dose (mg)",
                          300, 3000, 840, 60),
              plotOutput("p_taper", height = 360))),
        fluidRow(
          box(width = 12, title = "One-year relapse probability by taper length",
              status = "primary", solidHeader = TRUE,
              tableOutput("t_taper")))),

      ## ---- 11 NMOSD -------------------------------------------------------
      tabItem("nmosd",
        fluidRow(
          box(width = 12, title = "Same attack, same steroids, one switch (AQP4_MODE)",
              status = "primary", solidHeader = TRUE,
              plotOutput("p_nmosd", height = 420))),
        fluidRow(
          box(width = 12, status = "info", solidHeader = TRUE,
              title = "The built-in control experiment",
              HTML("<p>Setting <code>AQP4_MODE = 1</code> multiplies
                    oligodendrocyte death and astrocytic GFAP release. Nothing else
                    changes — same antibody kinetics, same barrier, same effectors,
                    same steroids. The MOGAD-versus-NMOSD divergence in recovery and
                    in serum GFAP therefore comes entirely from
                    <b>whether the target cell survives the attack</b>.</p>")))),

      ## ---- 12 notes -------------------------------------------------------
      tabItem("notes",
        fluidRow(
          box(width = 12, title = "Calibration anchors", status = "primary",
              solidHeader = TRUE,
              HTML("<ul>
                <li>Vilaseca 2026 <i>JAMA Neurol</i> (PMID 42440328): relapsing
                    MOGAD ARR 0.64 (0.58-0.70) before preventive therapy; 0.09
                    (0.06-0.14) on IL-6R blockade; 0.22 (0.15-0.32) on IVIG; IL-6R
                    blockade not significantly different from IVIG &ge;1 g/kg q4wk.</li>
                <li>Chen 2022 <i>JAMA Neurol</i> (PMID 35377395): relapse in 17% on
                    IVIG &ge;1 g/kg q4wk vs 50% on lower/less frequent dosing.</li>
                <li>Chen 2020 <i>Neurology</i> (PMID 32554760) and Thakolwiboon 2021
                    (PMID 34634625): on-treatment ARR ~0.59-0.63 rituximab,
                    ~0.67-0.84 mycophenolate, ~0.08 maintenance IVIG.</li>
                <li>Chen 2018 <i>Am J Ophthalmol</i> (PMID 30055153): MOG-ON nadir
                    ~count-fingers, average final VA 20/30, 6% &le;20/200.</li>
                <li>Stiebel-Kalish 2017 <i>PLoS One</i> (PMID 28125740): RNFL better
                    preserved after MOG-IgG than AQP4-IgG optic neuritis.</li>
                <li>Marignier 2025 <i>JNNP</i> (PMID 39939136): serum GFAP low in
                    MOGAD, unlike AQP4-IgG NMOSD.</li>
              </ul>
              <p>Full list with PubMed links:
                 <code>mogad_references.md</code>. Every simulated number quoted in
                 <code>README.md</code> is reproduced by
                 <code>python3 mogad_reference_impl.py</code>.</p>"))),
        fluidRow(
          box(width = 12, title = "Known misses", status = "danger",
              solidHeader = TRUE,
              HTML("<ol>
                <li>IL-6R blockade is under-predicted (model 0.23 vs observed
                    0.09) — either the retrospective incidence-rate ratio is
                    inflated by regression to the mean, or the IL-6 axis is more
                    load-bearing than the model allows.</li>
                <li>Azathioprine and mycophenolate share one arm, so their reported
                    ordering cannot be reproduced.</li>
                <li>FRAC_ESC is fitted, not measured — the source-versus-clearance
                    conclusion is contingent on it.</li>
                <li>Relapses are a hazard, not events: no attack actually fires, so
                    step-wise disability accrual is not simulated.</li>
                <li>No paediatric ADEM physiology; encephalopathy, seizures and
                    cognition appear in the mechanistic map but have no ODE.</li>
              </ol>
              <p><b>Research and education only.</b> Not for clinical decisions.</p>")))
      )
    )
  )
)

## ============================================================================
## Server
## ============================================================================
server <- function(input, output, session) {

  sim <- reactive({
    req(input$horizon)
    build_case(input, input$horizon)
  })

  sim_nmosd <- reactive({
    build_case(input, input$horizon, extra_par = list(AQP4_MODE = 1))
  })

  at_day <- function(df, d) df[which.min(abs(df$time - d)), ]

  ## ---- 1 overview ----------------------------------------------------------
  output$p_overview <- renderPlot({
    df <- sim()
    multiplot(df, c("Myelin_optic_nerve", "Axon_optic_nerve",
                    "Oligodendrocytes", "VA_logMAR"),
              c("Optic-nerve myelin", "Optic-nerve axons",
                "Oligodendrocytes", "VA (logMAR)"),
              "Attack, injury and recovery", "Fraction surviving / logMAR")
  })

  output$t_summary <- renderTable({
    df <- sim(); last <- at_day(df, input$horizon)
    data.frame(
      Readout = c("VA nadir (logMAR)", "VA at end (logMAR)",
                  "VA at end (Snellen 20/x)", "RNFL at end (um)",
                  "Min oligodendrocyte fraction", "Peak serum NfL (pg/mL)",
                  "Peak serum GFAP (pg/mL)", "Peak CSF cells (/uL)",
                  "ARR at end (/year)", "P(relapse) by end",
                  "Cumulative steroid (g)", "BMD at end (fraction)"),
      Value = c(sprintf("%.2f", max(df$VA_logMAR)),
                sprintf("%.2f", last$VA_logMAR),
                sprintf("20/%.0f", last$VA_snellen_denom),
                sprintf("%.0f", last$RNFL_thickness),
                sprintf("%.2f", min(df$Oligodendrocytes)),
                sprintf("%.0f", max(df$Serum_NfL)),
                sprintf("%.0f", max(df$Serum_GFAP)),
                sprintf("%.0f", max(df$CSF_cells)),
                sprintf("%.2f", last$ARR),
                sprintf("%.2f", last$Relapse_prob),
                sprintf("%.2f", last$Steroid_dose_g),
                sprintf("%.3f", last$BMD_fraction)),
      check.names = FALSE)
  })

  ## ---- 2 PK ---------------------------------------------------------------
  output$p_gc <- renderPlot({
    lineplot(sim(), "GC_occupancy",
             "Glucocorticoid receptor occupancy", "Fraction occupied")
  })
  output$p_ster <- renderPlot({
    multiplot(sim(), c("Steroid_dose_g", "BMD_fraction"),
              c("Cumulative dose (g)", "BMD (fraction of baseline)"),
              "Steroid exposure and its cost", "Value")
  })
  output$p_engage <- renderPlot({
    multiplot(sim(), c("CD20_occupancy", "IL6R_blockade",
                       "FcRn_occupancy", "C5_blockade"),
              c("CD20 occupancy", "IL-6R blockade",
                "FcRn occupancy", "C5 blockade"),
              "Target engagement", "Fraction")
  })
  output$p_igg <- renderPlot({
    lineplot(sim(), "Total_IgG", "Total serum IgG", "g/L")
  })

  ## ---- 3 source -----------------------------------------------------------
  output$p_bcell <- renderPlot({
    multiplot(sim(), c("Memory_B_cells", "Plasmablasts"),
              c("CD20+ memory B cells", "CD20- plasmablasts"),
              "Antibody source compartments", "Relative to baseline")
  })
  output$p_titer <- renderPlot({
    df <- sim() %>% mutate(active_scaled = MOG_IgG_active * 100)
    multiplot(df, c("MOG_IgG_titer", "active_scaled"),
              c("Measured titre (1:x)", "Functionally active (same scale)"),
              "Measured versus functionally active antibody",
              "Titre (assay dilution equivalent)")
  })

  ## ---- 4 effectors --------------------------------------------------------
  output$p_bbb <- renderPlot({
    lineplot(sim(), "BBB_permeability", "BBB permeability index", "0-1")
  })
  output$p_eff <- renderPlot({
    multiplot(sim(), c("CNS_MOG_IgG", "Complement_MAC", "Myeloid_activation"),
              c("CNS MOG-IgG", "Complement / MAC", "Microglia / macrophages"),
              "CNS effector arm", "Arbitrary units")
  })
  output$p_il6 <- renderPlot({
    lineplot(sim(), "IL6_measured", "Measured serum IL-6", "pg/mL")
  })
  output$p_edema <- renderPlot({
    lineplot(sim(), "Lesion_edema", "Lesion oedema", "Arbitrary units")
  })

  ## ---- 5 tissue -----------------------------------------------------------
  output$p_tissue <- renderPlot({
    multiplot(sim(), c("Myelin_optic_nerve", "Axon_optic_nerve",
                       "Myelin_cord", "Axon_cord"),
              c("Optic-nerve myelin", "Optic-nerve axons",
                "Cord myelin", "Cord axons"),
              "Myelin recovers, axons do not", "Fraction surviving")
  })
  output$p_ol <- renderPlot({
    lineplot(sim(), "Oligodendrocytes",
             "Surviving oligodendrocytes", "Fraction")
  })
  output$p_rnfl <- renderPlot({
    lineplot(sim(), "RNFL_thickness",
             "Peripapillary RNFL thickness", "um")
  })

  ## ---- 6 endpoints --------------------------------------------------------
  output$p_va <- renderPlot({
    lineplot(sim(), "VA_logMAR", "Visual acuity", "logMAR")
  })
  output$p_edss <- renderPlot({
    lineplot(sim(), "EDSS", "Expanded Disability Status Scale", "EDSS")
  })
  output$t_endpoint <- renderTable({
    df <- sim()
    days <- c(7, 14, 30, 60, 90, 180, min(365, input$horizon))
    days <- sort(unique(days[days <= input$horizon]))
    rows <- lapply(days, function(d) {
      r <- at_day(df, d)
      data.frame(Day = d,
                 `VA logMAR` = sprintf("%.2f", r$VA_logMAR),
                 Snellen = sprintf("20/%.0f", r$VA_snellen_denom),
                 EDSS = sprintf("%.1f", r$EDSS),
                 `RNFL (um)` = sprintf("%.0f", r$RNFL_thickness),
                 `Titre (1:x)` = sprintf("%.0f", r$MOG_IgG_titer),
                 `ARR (/y)` = sprintf("%.2f", r$ARR),
                 `NfL` = sprintf("%.0f", r$Serum_NfL),
                 `GFAP` = sprintf("%.0f", r$Serum_GFAP),
                 check.names = FALSE)
    })
    do.call(rbind, rows)
  })

  ## ---- 7 hazard -----------------------------------------------------------
  output$p_arr <- renderPlot({
    lineplot(sim(), "ARR", "Instantaneous annualised relapse rate", "ARR (/year)")
  })
  output$p_prel <- renderPlot({
    lineplot(sim(), "Relapse_prob", "Cumulative relapse probability",
             "Probability")
  })

  ## ---- 8 biomarkers -------------------------------------------------------
  output$p_nfl <- renderPlot({
    lineplot(sim(), "Serum_NfL", "Serum neurofilament light", "pg/mL")
  })
  output$p_gfap <- renderPlot({
    lineplot(sim(), "Serum_GFAP", "Serum GFAP", "pg/mL")
  })
  output$p_csf <- renderPlot({
    lineplot(sim(), "CSF_cells", "CSF white cell count", "cells/uL")
  })
  output$p_titer2 <- renderPlot({
    lineplot(sim(), "MOG_IgG_titer", "Serum MOG-IgG titre",
             "Assay dilution (1:x)")
  })

  ## ---- 9 comparison -------------------------------------------------------
  cmp <- reactive({
    horizon <- input$horizon
    bind_rows(lapply(MAINT, function(m) {
      build_case(input, horizon, maint = m) %>% mutate(strategy = m)
    }))
  })

  output$p_cmp <- renderPlot({
    ggplot(cmp(), aes(time, ARR, colour = strategy)) +
      geom_line(linewidth = 0.9) +
      labs(title = "Annualised relapse rate by maintenance strategy",
           x = "Time (days)", y = "ARR (/year)") + THEME +
      guides(colour = guide_legend(nrow = 3))
  })

  output$t_cmp <- DT::renderDataTable({
    horizon <- input$horizon
    ref <- NULL
    rows <- lapply(MAINT, function(m) {
      d <- build_case(input, horizon, maint = m)
      r <- d[which.min(abs(d$time - (horizon - 35))), ]
      data.frame(Strategy = m,
                 `Titre (1:x)` = round(r$MOG_IgG_titer),
                 `ARR (/y)` = round(r$ARR, 2),
                 `P(relapse)` = round(d[nrow(d), ]$Relapse_prob, 2),
                 `Total IgG (g/L)` = round(r$Total_IgG, 1),
                 `VA logMAR (end)` = round(d[nrow(d), ]$VA_logMAR, 2),
                 check.names = FALSE)
    })
    tab <- do.call(rbind, rows)
    base <- tab[tab$Strategy == MAINT[1], ]
    tab$`Titre vs none` <- round(tab$`Titre (1:x)` / base$`Titre (1:x)`, 2)
    tab$`ARR vs none` <- round(tab$`ARR (/y)` / base$`ARR (/y)`, 2)
    DT::datatable(tab, rownames = FALSE,
                  options = list(dom = "t", pageLength = 20))
  })

  ## ---- 10 taper -----------------------------------------------------------
  taper_tab <- reactive({
    tds <- c(14, 28, 45, 60, 90, 120, 180)
    horizon <- max(365, input$horizon)
    bind_rows(lapply(tds, function(td) {
      pars <- list(SEV = input$sev, SITE_ON = input$site_on,
                   SITE_SC = input$site_sc, FRAC_ESC = input$frac_esc,
                   PRED0 = 2 * input$total_mg / td, TAPER_DAYS = td,
                   PRED_START = input$ivmp_day + 5, MMF_RATE = 0)
      evs <- ev_bind(ev_attack(0), ev_ivmp(input$ivmp_day, input$ivmp_n))
      run_sim(pars, evs, horizon) %>% mutate(taper_days = factor(td))
    }))
  })

  output$p_taper <- renderPlot({
    ggplot(taper_tab(), aes(time, ARR, colour = taper_days)) +
      geom_line(linewidth = 0.9) +
      labs(title = paste0("Relapse hazard by taper length at a matched total dose of ",
                          input$total_mg, " mg"),
           x = "Time (days)", y = "ARR (/year)") + THEME
  })

  output$t_taper <- renderTable({
    taper_tab() %>%
      group_by(taper_days) %>%
      summarise(`Starting dose (mg/d)` = sprintf("%.1f", 2 * input$total_mg /
                                                   as.numeric(as.character(first(taper_days)))),
                `P(relapse) 6 mo` = sprintf("%.3f",
                    Relapse_prob[which.min(abs(time - 180))]),
                `P(relapse) 1 y` = sprintf("%.3f",
                    Relapse_prob[which.min(abs(time - 365))]),
                `BMD at 1 y` = sprintf("%.3f",
                    BMD_fraction[which.min(abs(time - 365))]),
                .groups = "drop") %>%
      rename(`Taper (days)` = taper_days)
  })

  ## ---- 11 NMOSD comparator ------------------------------------------------
  output$p_nmosd <- renderPlot({
    a <- sim() %>% mutate(arm = "MOGAD (AQP4_MODE = 0)")
    b <- sim_nmosd() %>% mutate(arm = "AQP4-IgG NMOSD (AQP4_MODE = 1)")
    both <- bind_rows(a, b) %>%
      select(time, arm, VA_logMAR, Oligodendrocytes, Serum_GFAP, RNFL_thickness) %>%
      pivot_longer(c(VA_logMAR, Oligodendrocytes, Serum_GFAP, RNFL_thickness),
                   names_to = "readout", values_to = "value") %>%
      mutate(readout = recode(readout,
                              VA_logMAR = "Visual acuity (logMAR)",
                              Oligodendrocytes = "Surviving oligodendrocytes",
                              Serum_GFAP = "Serum GFAP (pg/mL)",
                              RNFL_thickness = "RNFL (um)"))
    ggplot(both, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~readout, scales = "free_y") +
      labs(title = "MOGAD versus AQP4-IgG NMOSD — one parameter switch",
           x = "Time (days)", y = NULL) + THEME
  })
}

shinyApp(ui, server)
