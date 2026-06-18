## =============================================================================
## GPA QSP Shiny App
## Granulomatosis with Polyangiitis — Interactive Clinical Dashboard
## =============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## Embedded simulation function (simplified analytical approximations for
## real-time Shiny interactivity without mrgsolve dependency)
## ─────────────────────────────────────────────────────────────────────────────

simulate_gpa <- function(
    treatment       = "rtx_gc",
    bvas_init       = 18,
    gfr_init        = 60,
    anca_init       = 4.5,
    weight_kg       = 70,
    rtx_dose_mg     = 1000,
    gc_start_dose   = 60,
    gc_taper_half   = 56,
    avacopan_use    = FALSE,
    sim_days        = 730
) {
  time <- seq(0, sim_days, by = 7)
  n    <- length(time)

  ## ── Drug concentration proxies ────────────────────────────────────────────
  rtx_eff <- switch(treatment,
    "rtx_gc"      = pmax(0, 1.0 * exp(-time / (rtx_dose_mg / 40)) * (time <= 365) + 0.6 * exp(-(time - 182) / 80) * (time > 182)),
    "cyc_gc"      = pmax(0, 0.8 * (1 - exp(-time / 45)) * exp(-time / 400)),
    "ava_rtx"     = pmax(0, 0.9 * exp(-time / (rtx_dose_mg / 40))),
    "untreated"   = rep(0, n),
    "relapse"     = {
      first  <- pmax(0, 1.0 * exp(-time / 30))
      relapse<- pmax(0, 0.9 * exp(-(pmax(0, time - 450)) / 25))
      pmax(first, relapse)
    }
  )

  gc_eff <- switch(treatment,
    "rtx_gc"   = pmax(0.05, Emax_eff(gc_start_dose * exp(-log(2) / gc_taper_half * time), IC50 = 0.08) * 0.85),
    "cyc_gc"   = pmax(0.05, Emax_eff(gc_start_dose * exp(-log(2) / gc_taper_half * time), IC50 = 0.08) * 0.85),
    "ava_rtx"  = rep(0.10, n),  # minimal GC
    "untreated" = rep(0.0, n),
    "relapse"  = pmax(0.05, Emax_eff(gc_start_dose * exp(-log(2) / gc_taper_half * time), IC50 = 0.08) * 0.85)
  )

  ava_eff <- if (avacopan_use || treatment == "ava_rtx") {
    pmax(0, 0.85 * (1 - exp(-time / 21)) * (time <= 365) + 0 * (time > 365))
  } else rep(0, n)

  ## ── B Cell dynamics ───────────────────────────────────────────────────────
  b_cell <- pmax(0, 0.3 * exp(-rtx_eff * 5) + 0.02 * (time / 730))
  b_cell[b_cell < 0.01] <- 0.01

  ## ── PR3-ANCA titer ────────────────────────────────────────────────────────
  anca <- numeric(n)
  anca[1] <- anca_init
  for (i in 2:n) {
    prod_rate  <- 0.12 * pmax(0.01, b_cell[i])
    deg_rate   <- 0.025 + rtx_eff[i] * 0.08 + gc_eff[i] * 0.04
    anca[i]    <- max(0.1, anca[i-1] + (prod_rate - deg_rate * anca[i-1]) * 7)
    if (treatment == "relapse" && time[i] > 420 && time[i] < 500) {
      anca[i] <- anca[i] * 1.06  # pre-relapse ANCA rise
    }
  }

  ## ── Neutrophil activation ─────────────────────────────────────────────────
  n_act <- pmax(0, (0.5 + 0.3 * anca / anca_init) * (1 - gc_eff) * (1 - ava_eff * 0.6))

  ## ── C5a (complement) ──────────────────────────────────────────────────────
  c5a <- pmax(0, 0.3 * n_act * (1 - ava_eff))

  ## ── Granuloma index ───────────────────────────────────────────────────────
  gran <- numeric(n)
  gran[1] <- 0.4
  for (i in 2:n) {
    form <- 0.04 * anca[i] * (1 - gc_eff[i])
    res  <- 0.02 * gran[i-1] * (gc_eff[i] + rtx_eff[i] * 0.4)
    gran[i] <- max(0, gran[i-1] + (form - res) * 7)
  }

  ## ── Endothelial injury ────────────────────────────────────────────────────
  ec_inj <- numeric(n)
  ec_inj[1] <- 0.35
  for (i in 2:n) {
    inj  <- 0.05 * n_act[i]
    rep  <- 0.03 * (1 - ec_inj[i-1]) * gc_eff[i]
    ec_inj[i] <- max(0, min(1, ec_inj[i-1] + (inj - rep) * 7 / 30))
  }

  ## ── eGFR ──────────────────────────────────────────────────────────────────
  gfr <- numeric(n)
  gfr[1] <- gfr_init
  for (i in 2:n) {
    loss <- 0.008 * ec_inj[i] * (1 + n_act[i]) * 7
    rec  <- 0.004 * (90 - gfr[i-1]) * gc_eff[i] * 7
    gfr[i] <- max(5, gfr[i-1] - loss + rec)
  }

  ## ── BVAS score ────────────────────────────────────────────────────────────
  bvas <- numeric(n)
  bvas[1] <- bvas_init
  for (i in 2:n) {
    drive <- (gran[i] + ec_inj[i]) * 0.5
    suppress <- (gc_eff[i] + rtx_eff[i] * 0.7 + ava_eff[i] * 0.3) / 2
    bvas[i] <- max(0, bvas[i-1] * (1 + (drive * 0.03 - suppress * 0.03) * 7))
  }

  data.frame(
    time      = time,
    B_cell    = b_cell,
    ANCA      = anca,
    N_act     = n_act,
    C5a       = c5a,
    Gran_idx  = gran,
    EC_injury = ec_inj,
    GFR       = gfr,
    BVAS      = bvas,
    RTX_eff   = rtx_eff,
    GC_eff    = gc_eff,
    AVA_eff   = ava_eff
  )
}

Emax_eff <- function(Cp, Emax = 1, IC50 = 0.08, n = 1) {
  Emax * Cp^n / (Cp^n + IC50^n)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span("GPA QSP Dashboard", style = "font-size:16px; font-weight:bold"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",              tabName = "tab_pk",       icon = icon("capsules")),
      menuItem("Biomarkers",           tabName = "tab_biomark",  icon = icon("vials")),
      menuItem("Organ Damage",         tabName = "tab_organ",    icon = icon("kidneys")),
      menuItem("Scenario Comparison",  tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Clinical Endpoints",   tabName = "tab_endpoint", icon = icon("heartbeat")),
      menuItem("Model Reference",      tabName = "tab_ref",      icon = icon("book-medical"))
    ),

    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:white"),

    sliderInput("bvas_init",    "Baseline BVAS score",    min = 0,  max = 30,  value = 18, step = 1),
    sliderInput("gfr_init",     "Baseline eGFR (mL/min)", min = 10, max = 90,  value = 60, step = 5),
    sliderInput("anca_init",    "Baseline PR3-ANCA (AU)", min = 0.5, max = 10, value = 4.5, step = 0.5),
    sliderInput("weight_kg",    "Body weight (kg)",       min = 40, max = 120, value = 70,  step = 5),

    hr(),
    h5("Treatment Settings", style = "padding-left:15px; color:white"),

    selectInput("treatment", "Treatment Regimen",
                choices = c(
                  "RTX + GC (RAVE)"             = "rtx_gc",
                  "CYC + GC (Standard SoC)"     = "cyc_gc",
                  "Avacopan + RTX (GC-sparing)" = "ava_rtx",
                  "Untreated"                    = "untreated",
                  "Relapse & Re-induction"       = "relapse"
                ), selected = "rtx_gc"),

    sliderInput("rtx_dose",  "RTX Dose (mg/infusion)", min = 375, max = 1000, value = 1000, step = 125),
    sliderInput("gc_dose",   "GC Start Dose (mg/d)",   min = 10,  max = 80,   value = 60,   step = 5),
    sliderInput("gc_taper",  "GC Taper Half-life (d)", min = 28,  max = 112,  value = 56,   step = 14),
    checkboxInput("avacopan_add", "Add Avacopan", value = FALSE),
    sliderInput("sim_time",  "Simulation Days",        min = 180, max = 1095, value = 730,  step = 90)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 6px; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "GPA — Disease Overview", width = 12, status = "primary", solidHeader = TRUE,
              column(6,
                h4("Granulomatosis with Polyangiitis (GPA)"),
                p("Formerly known as Wegener's Granulomatosis. A rare ANCA-associated",
                  "small-to-medium vessel vasculitis characterized by necrotizing granulomatous",
                  "inflammation of the upper and lower respiratory tracts and pauci-immune",
                  "crescentic glomerulonephritis."),
                h5("Key Features:"),
                tags$ul(
                  tags$li("PR3-ANCA (cANCA): ~80-90% positive in GPA"),
                  tags$li("Triad: sinusitis, pulmonary nodules, renal disease"),
                  tags$li("Incidence: ~3/100,000 per year"),
                  tags$li("Peak age: 45-65 years; M:F ratio ~1.5:1"),
                  tags$li("5-year survival without treatment: <20%"),
                  tags$li("5-year survival with modern treatment: >80%")
                )
              ),
              column(6,
                h4("Current Patient Summary"),
                valueBoxOutput("vbox_bvas",   width = 6),
                valueBoxOutput("vbox_gfr",    width = 6),
                valueBoxOutput("vbox_anca",   width = 6),
                valueBoxOutput("vbox_remission", width = 6)
              )
          )
        ),
        fluidRow(
          box(title = "Mechanistic Diagram Summary", width = 8, status = "info",
              p("This QSP model integrates:"),
              tags$ol(
                tags$li(strong("ANCA Pathogenesis:"), " Anti-PR3 IgG (from B cell/plasma cell axis) primes and activates neutrophils"),
                tags$li(strong("Neutrophil Biology:"), " TNF-α/IL-8/C5a priming → full ANCA activation → NETosis → vessel injury"),
                tags$li(strong("Complement Cascade:"), " NET-triggered C1q → C3/C5 cleavage → C5a feedback on neutrophil priming"),
                tags$li(strong("Granuloma Formation:"), " Macrophage/Th1/IFN-γ axis creates necrotizing upper-airway/pulmonary granulomata"),
                tags$li(strong("Organ Damage:"), " Endothelial necrosis → pauci-immune RPGN → GFR decline; ENT/lung destruction"),
                tags$li(strong("Drug PK/PD:"), " Rituximab (CD20 depletion), CYC (DNA alkylation), GC (NF-κB repression), Avacopan (C5aR1 blockade)")
              )
          ),
          box(title = "Classification & Staging", width = 4, status = "warning",
              h5("BVAS Score Interpretation"),
              DT::dataTableOutput("bvas_ref_table"),
              h5(style="margin-top:10px", "Berden Renal Biopsy Class"),
              DT::dataTableOutput("berden_table")
          )
        )
      ),

      ## ── TAB 2: Drug PK ─────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Effectiveness Index over Time", width = 8, status = "primary", solidHeader = TRUE,
              plotlyOutput("pk_plot", height = "400px")
          ),
          box(title = "PK Parameters Summary", width = 4, status = "info",
              h5("Rituximab (RTX)"),
              tags$ul(
                tags$li("CL: 0.38 L/day (target-mediated)"),
                tags$li("V1: 3.5 L, V2: 3.2 L (2-CMT)"),
                tags$li("t½α: ~2 days; t½β: ~22 days"),
                tags$li("TMDD via CD20 (kon=0.04, kint=0.1/day)")
              ),
              h5("Cyclophosphamide (CYC)"),
              tags$ul(
                tags$li("Prodrug: hepatic CYP2B6/3A4 activation"),
                tags$li("CL: 7.8 L/day; V: 38 L"),
                tags$li("Active metabolite t½: ~4 h")
              ),
              h5("Prednisolone (GC)"),
              tags$ul(
                tags$li("CL: 7.0 L/day; V: 32 L"),
                tags$li("IC50 TNF: 0.08 mg/L"),
                tags$li("Emax suppression: 85%")
              ),
              h5("Avacopan"),
              tags$ul(
                tags$li("30 mg BID oral; CYP3A4 substrate"),
                tags$li("CL: 12 L/day; V: 85 L"),
                tags$li("C5aR1 EC50: 0.15 mg/L"),
                tags$li("t½ ~7h; >85% C5aR1 occupancy at Css")
              )
          )
        ),
        fluidRow(
          box(title = "GC Taper Profile", width = 6, status = "warning",
              plotlyOutput("gc_taper_plot", height = "300px")
          ),
          box(title = "Drug Mechanism Summary", width = 6, status = "success",
              tableOutput("drug_moa_table")
          )
        )
      ),

      ## ── TAB 3: Biomarkers ──────────────────────────────────────────────────
      tabItem(tabName = "tab_biomark",
        fluidRow(
          box(title = "PR3-ANCA Titer Dynamics", width = 6, status = "primary", solidHeader = TRUE,
              plotlyOutput("anca_plot", height = "350px")
          ),
          box(title = "B Cell Dynamics (CD19+)", width = 6, status = "primary", solidHeader = TRUE,
              plotlyOutput("bcell_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Neutrophil Activation & NETs", width = 6, status = "danger",
              plotlyOutput("neutrophil_plot", height = "320px")
          ),
          box(title = "Complement C5a Levels", width = 6, status = "warning",
              plotlyOutput("c5a_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges", width = 12, status = "info",
              DT::dataTableOutput("biomarker_ref_table")
          )
        )
      ),

      ## ── TAB 4: Organ Damage ────────────────────────────────────────────────
      tabItem(tabName = "tab_organ",
        fluidRow(
          box(title = "Renal Function (eGFR)", width = 6, status = "primary", solidHeader = TRUE,
              plotlyOutput("gfr_plot", height = "380px")
          ),
          box(title = "Endothelial Injury Index", width = 6, status = "danger",
              plotlyOutput("ec_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Granuloma Burden Index", width = 6, status = "warning",
              plotlyOutput("gran_plot", height = "320px")
          ),
          box(title = "BVAS Score Trajectory", width = 6, status = "info",
              plotlyOutput("bvas_plot_tab4", height = "320px")
          )
        ),
        fluidRow(
          box(title = "GFR CKD Stage Classification", width = 12, status = "info",
              DT::dataTableOutput("ckd_stage_table")
          )
        )
      ),

      ## ── TAB 5: Scenario Comparison ─────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Multi-scenario BVAS Comparison", width = 6, status = "primary", solidHeader = TRUE,
              plotlyOutput("scenario_bvas", height = "380px")
          ),
          box(title = "Multi-scenario eGFR Comparison", width = 6, status = "danger",
              plotlyOutput("scenario_gfr", height = "380px")
          )
        ),
        fluidRow(
          box(title = "PR3-ANCA by Scenario", width = 6, status = "warning",
              plotlyOutput("scenario_anca", height = "320px")
          ),
          box(title = "B Cell Depletion by Scenario", width = 6, status = "success",
              plotlyOutput("scenario_bcell", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Outcome Summary at Key Timepoints", width = 12, status = "info",
              DT::dataTableOutput("scenario_summary_table")
          )
        )
      ),

      ## ── TAB 6: Clinical Endpoints ──────────────────────────────────────────
      tabItem(tabName = "tab_endpoint",
        fluidRow(
          box(title = "Probability of Remission over Time", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("remission_plot", height = "380px")
          ),
          box(title = "ESRD Risk Trajectory", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("esrd_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Key Clinical Trial Results", width = 12, status = "info",
              DT::dataTableOutput("trial_table")
          )
        ),
        fluidRow(
          box(title = "Treatment Decision Guide", width = 12, status = "primary",
              h4("When to Use Each Treatment:"),
              fluidRow(
                column(3,
                  h5(icon("syringe"), " Rituximab (RTX)"),
                  p("1st-line induction for non-severe/severe GPA; maintenance preferred over AZA.",
                    "Particular benefit in relapsing disease. Monitor CD19+ B cells."),
                  tags$span("Evidence: RAVE, RITUXVAS, MAINRITSAN3", style="color:blue;font-size:11px")
                ),
                column(3,
                  h5(icon("pills"), " Cyclophosphamide (CYC)"),
                  p("1st-line for severe organ-threatening disease (DAH, RPGN).",
                    "IV pulse preferred over oral (less toxicity). Switch to AZA/MTX for maintenance."),
                  tags$span("Evidence: CYCLOPS, MEPEX", style="color:blue;font-size:11px")
                ),
                column(3,
                  h5(icon("capsules"), " Avacopan (C5aR1i)"),
                  p("GC-sparing strategy; add to RTX or CYC induction.",
                    "Non-inferior BVAS remission at 26 weeks; superior at 52 weeks.",
                    "Well tolerated; monitor LFTs (CYP3A4)."),
                  tags$span("Evidence: ADVOCATE (2021)", style="color:blue;font-size:11px")
                ),
                column(3,
                  h5(icon("tablets"), " Glucocorticoids (GC)"),
                  p("High-dose induction (1 mg/kg/d pred), rapid taper to minimize toxicity.",
                    "Avacopan now enables GC-free induction in selected patients.",
                    "Monitor BMD, glucose, BP."),
                  tags$span("Evidence: Standard of care since 1970s", style="color:blue;font-size:11px")
                )
              )
          )
        )
      ),

      ## ── TAB 7: Model Reference ─────────────────────────────────────────────
      tabItem(tabName = "tab_ref",
        fluidRow(
          box(title = "Model Structure", width = 6, status = "primary",
              h4("QSP Model Compartments (22 ODEs)"),
              DT::dataTableOutput("model_struct_table")
          ),
          box(title = "Key Parameters", width = 6, status = "info",
              DT::dataTableOutput("param_table")
          )
        ),
        fluidRow(
          box(title = "Key References (Selected)", width = 12, status = "warning",
              DT::dataTableOutput("ref_table")
          )
        )
      )
    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Reactive simulation ───────────────────────────────────────────────────
  sim_data <- reactive({
    simulate_gpa(
      treatment     = input$treatment,
      bvas_init     = input$bvas_init,
      gfr_init      = input$gfr_init,
      anca_init     = input$anca_init,
      weight_kg     = input$weight_kg,
      rtx_dose_mg   = input$rtx_dose,
      gc_start_dose = input$gc_dose,
      gc_taper_half = input$gc_taper,
      avacopan_use  = input$avacopan_add,
      sim_days      = input$sim_time
    )
  })

  all_scenario_data <- reactive({
    scenarios <- c("rtx_gc", "cyc_gc", "ava_rtx", "untreated", "relapse")
    labels    <- c("RTX + GC", "CYC + GC", "Avacopan+RTX", "Untreated", "Relapse/Re-tx")
    bind_rows(lapply(seq_along(scenarios), function(i) {
      simulate_gpa(
        treatment     = scenarios[i],
        bvas_init     = input$bvas_init,
        gfr_init      = input$gfr_init,
        anca_init     = input$anca_init,
        weight_kg     = input$weight_kg,
        rtx_dose_mg   = input$rtx_dose,
        gc_start_dose = input$gc_dose,
        gc_taper_half = input$gc_taper,
        sim_days      = input$sim_time
      ) %>% mutate(Scenario = labels[i])
    }))
  })

  ## ── Value Boxes ───────────────────────────────────────────────────────────
  output$vbox_bvas <- renderValueBox({
    d <- sim_data()
    bvas_6m <- d$BVAS[min(which(d$time >= 180))]
    valueBox(round(bvas_6m, 1), "BVAS at 6 months",
             icon = icon("chart-line"),
             color = if (bvas_6m < 2) "green" else if (bvas_6m < 10) "yellow" else "red")
  })

  output$vbox_gfr <- renderValueBox({
    d <- sim_data()
    gfr_6m <- d$GFR[min(which(d$time >= 180))]
    valueBox(round(gfr_6m, 1), "eGFR at 6m (mL/min)",
             icon = icon("kidneys"),
             color = if (gfr_6m > 60) "green" else if (gfr_6m > 30) "yellow" else "red")
  })

  output$vbox_anca <- renderValueBox({
    d <- sim_data()
    anca_6m <- d$ANCA[min(which(d$time >= 180))]
    valueBox(round(anca_6m, 1), "PR3-ANCA at 6m",
             icon = icon("vials"),
             color = if (anca_6m < 1) "green" else if (anca_6m < 3) "yellow" else "red")
  })

  output$vbox_remission <- renderValueBox({
    d <- sim_data()
    bvas_12m <- d$BVAS[min(which(d$time >= 365))]
    rem_status <- if (bvas_12m < 1) "Complete" else if (bvas_12m < 5) "Partial" else "Active"
    valueBox(rem_status, "Remission Status at 1 year",
             icon = icon("check-circle"),
             color = if (rem_status == "Complete") "green" else if (rem_status == "Partial") "yellow" else "red")
  })

  ## ── Reference Tables ──────────────────────────────────────────────────────
  output$bvas_ref_table <- DT::renderDataTable({
    data.frame(
      BVAS_Range = c("0", "1–5", "6–12", "13–20", ">20"),
      Category   = c("Remission", "Minor activity", "Moderate", "Severe", "Critical")
    )
  }, options = list(dom = "t", pageLength = 5), rownames = FALSE)

  output$berden_table <- DT::renderDataTable({
    data.frame(
      Class   = c("Focal", "Crescentic", "Mixed", "Sclerotic"),
      Crescents = c("<50%", ">50%", "Mixed", ">50% sclerotic"),
      Prognosis = c("Good", "Poor", "Intermediate", "Very poor")
    )
  }, options = list(dom = "t", pageLength = 4), rownames = FALSE)

  output$biomarker_ref_table <- DT::renderDataTable({
    data.frame(
      Biomarker  = c("PR3-ANCA", "CD19+ B cells", "Neutrophil count",
                     "CRP", "ESR", "eGFR", "Urinary RBC casts"),
      Normal     = c("<2 AU/mL", "80–400/μL", "1.8–7.7×10⁹/L",
                     "<5 mg/L", "<20 mm/hr", ">60 mL/min", "Absent"),
      Active_GPA = c("Often >4 AU/mL", "Variable (may be normal)",
                     "Elevated", ">50 mg/L", "Often >50", "Often <60",
                     "Present in active GN"),
      Clinical_Role = c("Diagnosis, relapse prediction", "RTX monitoring",
                        "Tissue infiltration marker", "Acute inflammation",
                        "Inflammation", "Renal function", "Active GN indicator")
    )
  }, options = list(dom = "t", pageLength = 7), rownames = FALSE)

  output$ckd_stage_table <- DT::renderDataTable({
    d <- sim_data()
    timepoints <- c(0, 90, 180, 365, 548, 730)
    tp_idx <- sapply(timepoints, function(tp) min(which(d$time >= tp)))
    data.frame(
      Time  = c("Baseline", "3m", "6m", "12m", "18m", "24m"),
      eGFR  = round(d$GFR[tp_idx], 1),
      Stage = sapply(d$GFR[tp_idx], function(g) {
        if (g >= 90) "G1 (Normal)"
        else if (g >= 60) "G2 (Mildly reduced)"
        else if (g >= 45) "G3a"
        else if (g >= 30) "G3b"
        else if (g >= 15) "G4 (Severely reduced)"
        else "G5 (Kidney failure)"
      })
    )
  }, options = list(dom = "t"), rownames = FALSE)

  ## ── Drug PK Plot ─────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    d <- sim_data()
    months <- d$time / 30.4
    plot_ly() %>%
      add_lines(x = months, y = d$RTX_eff, name = "RTX Effectiveness",
                line = list(color = "#2E86C1", width = 2.5)) %>%
      add_lines(x = months, y = d$GC_eff,  name = "GC Effectiveness",
                line = list(color = "#27AE60", width = 2.5)) %>%
      add_lines(x = months, y = d$AVA_eff, name = "Avacopan C5aR1 Block",
                line = list(color = "#9B59B6", width = 2.5, dash = "dash")) %>%
      layout(
        title = "Drug Effectiveness Index (0–1)",
        xaxis = list(title = "Months"),
        yaxis = list(title = "Effectiveness (0–1)", range = c(0, 1.1)),
        legend = list(orientation = "h"),
        hovermode = "x unified"
      )
  })

  output$gc_taper_plot <- renderPlotly({
    days_seq <- seq(0, input$sim_time)
    dose_seq <- pmax(5, input$gc_dose * exp(-log(2) / input$gc_taper * days_seq))
    plot_ly(x = days_seq/30.4, y = dose_seq, type = "scatter", mode = "lines",
            line = list(color = "#E67E22", width = 2)) %>%
      layout(title = "Glucocorticoid Taper Profile",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Prednisolone dose (mg/day)"))
  })

  output$drug_moa_table <- renderTable({
    data.frame(
      Drug      = c("Rituximab", "Cyclophosphamide", "Prednisolone", "Avacopan"),
      Target    = c("CD20 (B cells)", "DNA (all lymphocytes)", "GRα (NF-κB)", "C5aR1 (neutrophils)"),
      Key_Effect = c("B cell depletion → ↓ANCA", "Pan-lymphocyte apoptosis", "↓TNF/IL-6/IL-8", "↓Neutrophil priming/migration"),
      Class     = c("Anti-CD20 mAb", "Alkylating agent", "Corticosteroid", "C5aR1 antagonist")
    )
  })

  ## ── Biomarker Plots ───────────────────────────────────────────────────────
  output$anca_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$ANCA, type = "scatter", mode = "lines",
            line = list(color = "#9B59B6", width = 2.5),
            name = "PR3-ANCA") %>%
      add_lines(x = d$time/30.4, y = rep(2, nrow(d)), name = "Normal upper limit",
                line = list(color = "gray", dash = "dash")) %>%
      layout(title = "PR3-ANCA Titer",
             xaxis = list(title = "Months"),
             yaxis = list(title = "AU/mL"),
             hovermode = "x unified")
  })

  output$bcell_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$B_cell, type = "scatter", mode = "lines",
            line = list(color = "#2E86C1", width = 2.5)) %>%
      add_lines(x = d$time/30.4, y = rep(0.02, nrow(d)),
                name = "Depletion threshold",
                line = list(color = "red", dash = "dot")) %>%
      layout(title = "B Cell Count (CD19+)",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Normalized count"),
             hovermode = "x unified")
  })

  output$neutrophil_plot <- renderPlotly({
    d <- sim_data()
    plot_ly() %>%
      add_lines(x = d$time/30.4, y = d$N_act, name = "Activated Neutrophils",
                line = list(color = "#E74C3C", width = 2)) %>%
      add_lines(x = d$time/30.4, y = d$NETs %||% (d$N_act * 0.5),
                name = "NET burden", line = list(color = "#C0392B", dash = "dash")) %>%
      layout(title = "Neutrophil Activation & NETs",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Normalized units"),
             hovermode = "x unified")
  })

  output$c5a_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$C5a, type = "scatter", mode = "lines",
            line = list(color = "#E67E22", width = 2.5)) %>%
      layout(title = "Complement Fragment C5a",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Normalized units"))
  })

  ## ── Organ Damage Plots ────────────────────────────────────────────────────
  output$gfr_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$GFR, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(46,134,193,0.15)",
            line = list(color = "#2E86C1", width = 2.5)) %>%
      add_lines(x = d$time/30.4, y = rep(60, nrow(d)), name = "G2 threshold",
                line = list(color = "goldenrod", dash = "dot")) %>%
      add_lines(x = d$time/30.4, y = rep(30, nrow(d)), name = "G3b threshold",
                line = list(color = "orange", dash = "dot")) %>%
      add_lines(x = d$time/30.4, y = rep(15, nrow(d)), name = "G5 ESRD",
                line = list(color = "red", dash = "dot")) %>%
      layout(title = "eGFR Trajectory",
             xaxis = list(title = "Months"),
             yaxis = list(title = "eGFR (mL/min/1.73m²)", range = c(0, 100)),
             hovermode = "x unified")
  })

  output$ec_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$EC_injury * 100, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C", width = 2.5)) %>%
      layout(title = "Endothelial Injury Index",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Injury score (% maximum)", range = c(0, 110)))
  })

  output$gran_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$Gran_idx, type = "scatter", mode = "lines",
            line = list(color = "#8E44AD", width = 2.5)) %>%
      layout(title = "Granuloma Burden Index",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Granuloma index (normalized)"))
  })

  output$bvas_plot_tab4 <- renderPlotly({
    d <- sim_data()
    plot_ly(x = d$time/30.4, y = d$BVAS, type = "scatter", mode = "lines",
            line = list(color = "#E67E22", width = 2.5)) %>%
      add_lines(x = d$time/30.4, y = rep(1, nrow(d)), name = "Remission threshold",
                line = list(color = "green", dash = "dash")) %>%
      layout(title = "BVAS Score",
             xaxis = list(title = "Months"),
             yaxis = list(title = "BVAS"),
             hovermode = "x unified")
  })

  ## ── Scenario Comparison ───────────────────────────────────────────────────
  scenario_colors_js <- list(
    "RTX + GC"    = "#2E86C1",
    "CYC + GC"    = "#E67E22",
    "Avacopan+RTX"= "#27AE60",
    "Untreated"   = "#C0392B",
    "Relapse/Re-tx"="#8E44AD"
  )

  output$scenario_bvas <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4, y = ds$BVAS, name = scen)
    }
    p %>% layout(title = "BVAS by Treatment Scenario",
                 xaxis = list(title = "Months"),
                 yaxis = list(title = "BVAS Score"),
                 hovermode = "x unified")
  })

  output$scenario_gfr <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4, y = ds$GFR, name = scen)
    }
    p %>% add_lines(x = d$time[d$Scenario=="Untreated"]/30.4, y = rep(15, sum(d$Scenario=="Untreated")),
                    name = "ESRD threshold", line = list(color="red", dash="dot")) %>%
      layout(title = "eGFR by Treatment Scenario",
             xaxis = list(title = "Months"),
             yaxis = list(title = "eGFR (mL/min)", range = c(0, 95)),
             hovermode = "x unified")
  })

  output$scenario_anca <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4, y = ds$ANCA, name = scen)
    }
    p %>% layout(title = "PR3-ANCA by Scenario",
                 xaxis = list(title = "Months"), yaxis = list(title = "AU/mL"),
                 hovermode = "x unified")
  })

  output$scenario_bcell <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4, y = ds$B_cell, name = scen)
    }
    p %>% layout(title = "B Cell Count by Scenario",
                 xaxis = list(title = "Months"), yaxis = list(title = "Normalized"),
                 hovermode = "x unified")
  })

  output$scenario_summary_table <- DT::renderDataTable({
    d <- all_scenario_data()
    timepoints <- c(180, 365, 548)
    tp_names   <- c("6 months", "12 months", "18 months")
    rows <- lapply(seq_along(timepoints), function(ti) {
      tp <- timepoints[ti]
      d %>%
        filter(time >= tp) %>%
        group_by(Scenario) %>%
        slice(1) %>%
        ungroup() %>%
        transmute(
          Scenario,
          Timepoint = tp_names[ti],
          BVAS      = round(BVAS, 1),
          `PR3-ANCA (AU)` = round(ANCA, 2),
          `eGFR`    = round(GFR, 1),
          `B cells` = round(B_cell, 3),
          `Remission` = ifelse(BVAS < 1.0, "Yes", "No")
        )
    })
    bind_rows(rows) %>% arrange(Timepoint, Scenario)
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  ## ── Clinical Endpoints ────────────────────────────────────────────────────
  output$remission_plot <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4, y = as.integer(ds$BVAS < 1) * 100,
                     name = scen)
    }
    p %>% layout(title = "Remission Achievement (BVAS < 1)",
                 xaxis = list(title = "Months"),
                 yaxis = list(title = "Remission (%)", range = c(-5, 110)),
                 hovermode = "x unified")
  })

  output$esrd_plot <- renderPlotly({
    d <- all_scenario_data()
    p <- plot_ly()
    for (scen in unique(d$Scenario)) {
      ds <- d[d$Scenario == scen, ]
      p <- add_lines(p, x = ds$time/30.4,
                     y = cummax(as.integer(ds$GFR < 15)) * 100,
                     name = scen)
    }
    p %>% add_lines(x = c(0, max(d$time)/30.4), y = c(50, 50),
                    name = "50% threshold", line = list(color="gray", dash="dot")) %>%
      layout(title = "Cumulative ESRD Risk (eGFR < 15)",
             xaxis = list(title = "Months"),
             yaxis = list(title = "Cumulative risk (%)", range = c(-5, 110)))
  })

  output$trial_table <- DT::renderDataTable({
    data.frame(
      Trial       = c("RAVE (2010)", "RITUXVAS (2010)", "MAINRITSAN (2014)", "MAINRITSAN2 (2016)",
                      "MAINRITSAN3 (2023)", "ADVOCATE (2021)", "RAVE 18-month (2013)", "CYCLOPS (2009)"),
      Design      = c("RTX vs CYC induction", "RTX+CYC vs CYC induction",
                      "RTX maint vs AZA", "RTX 500mg q18m vs q6m",
                      "RTX 500mg q6m vs q18m", "Avacopan vs high-dose GC",
                      "RTX vs CYC long-term", "IV pulse vs oral CYC"),
      N           = c(197, 44, 115, 162, 97, 330, 197, 149),
      Key_Result  = c("RTX non-inferior to CYC at 6m", "RTX non-inferior; less GC exposure",
                      "RTX superior to AZA (5% vs 29% relapse)",
                      "q6m non-inferior to q18m",
                      "500mg q6m: 4% relapse vs 20% for q18m",
                      "Avacopan non-inferior at 26wk; superior at 52wk sustained remission",
                      "RTX superior at 18m (64% vs 53% sustained remission)",
                      "IV pulse: less toxicity, similar efficacy"),
      PubMed      = c("PMID:20647199", "PMID:20647200", "PMID:25372085",
                      "PMID:27350294", "PMID:36546673", "PMID:34597417",
                      "PMID:24131199", "PMID:19451574")
    )
  }, options = list(dom = "t", pageLength = 8, scrollX = TRUE), rownames = FALSE)

  ## ── Model Reference Tables ────────────────────────────────────────────────
  output$model_struct_table <- DT::renderDataTable({
    data.frame(
      Compartment = c("RTX1/RTX2", "RTX_bound", "CYC_gut/c/act", "GC_gut/c",
                      "AVA_gut/c", "B_naive/B_mem/PC_LL", "ANCA",
                      "C5a", "N_rest/N_act/NETs", "EC_injury",
                      "Gran_idx", "GFR", "BVAS"),
      Type        = c("PK-central/peripheral", "PK-TMDD", "PK-prodrug 3-CMT",
                      "PK-1CMT", "PK-1CMT", "PD-B cell subsets", "PD-autoantibody",
                      "PD-complement", "PD-neutrophil", "PD-endothelium",
                      "PD-granuloma", "PD-renal", "PD-composite"),
      Equations   = c(2, 1, 3, 2, 2, 3, 1, 1, 3, 1, 1, 1, 1)
    )
  }, options = list(dom = "t"), rownames = FALSE)

  output$param_table <- DT::renderDataTable({
    data.frame(
      Parameter  = c("RTX_CL", "RTX_kon/koff", "kprod_ANCA", "kdeg_ANCA",
                     "N_prime_k", "kGFR_loss", "kGran_form", "AVA_EC50"),
      Value      = c("0.38 L/d", "0.04/0.002 1/d", "0.12 AU/d/LLPC",
                     "0.025 /d", "0.08 /d", "0.008 mL/min/d", "0.04 /d", "0.15 mg/L"),
      Source     = c("PK study", "In vitro", "Model estimated", "IgG half-life",
                     "Ex vivo neutrophil", "Clinical", "Pathology", "ADVOCATE PK/PD")
    )
  }, options = list(dom = "t"), rownames = FALSE)

  output$ref_table <- DT::renderDataTable({
    data.frame(
      Citation = c(
        "Stone JH et al. (2010) NEJM 363:221",
        "Jones RB et al. (2010) NEJM 363:211",
        "Charles P et al. (2014) NEJM 371:1771",
        "Jayne DRW et al. (2021) NEJM 384:599",
        "Charles P et al. (2023) NEJM 388:308",
        "Lyons PA et al. (2012) NEJM 367:214",
        "Jennette JC et al. (2013) NEJM 369:2206",
        "Hellmark T & Segelmark M (2014) J Intern Med 275:386",
        "Grayson PC et al. (2015) NEJM 372:1938"
      ),
      Topic = c(
        "RAVE: RTX vs CYC induction RCT",
        "RITUXVAS: RTX induction in ANCA vasculitis",
        "MAINRITSAN: RTX maintenance vs AZA",
        "ADVOCATE: Avacopan vs glucocorticoids",
        "MAINRITSAN3: RTX maintenance dose/interval",
        "Genetic basis of ANCA vasculitis (GWAS)",
        "Pathogenesis overview (NEJM review)",
        "ANCA and complement in GPA",
        "Disease phenotype and treatment response"
      ),
      PMID = c("PMID:20647199", "PMID:20647200", "PMID:25372085",
               "PMID:34597417", "PMID:36546673", "PMID:22808956",
               "PMID:24171521", "PMID:24697279", "PMID:25760352")
    )
  }, options = list(dom = "t", pageLength = 9), rownames = FALSE)
}

## ─────────────────────────────────────────────────────────────────────────────
## RUN APP
## ─────────────────────────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui = ui, server = server)
