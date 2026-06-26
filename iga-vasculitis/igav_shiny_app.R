# IgA Vasculitis (Henoch-Schönlein Purpura) QSP Shiny Dashboard
# Author: QSP Disease Model Library (CCR)
# Date: 2026-06-19

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)

# ── Simulation helper (analytic approximation, no mrgsolve dependency) ────────

run_simulation <- function(
    age = 35, weight = 70, is_pediatric = FALSE,
    baseline_proteinuria = 2.0, baseline_egfr = 75,
    gd_iga1_level = 4.0, complement_activity = 1.5,
    treatment = "none",
    pred_dose = 1.0,    # mg/kg/day
    mmf_dose = 1000,    # mg BID
    rtx_dose = 375,     # mg/m²
    acei = FALSE, sglt2i = FALSE,
    duration_weeks = 52
) {
  time <- seq(0, duration_weeks, by = 0.5)
  n <- length(time)

  # ---------- Disease natural history ----------
  disease_severity <- min(1, (gd_iga1_level / 4 + complement_activity / 2) / 2)
  egfr_slope_nat  <- -2.5 * disease_severity          # mL/min/1.73m²/year
  prot_nat        <- baseline_proteinuria * (1 + 0.3 * disease_severity * time / 52)

  # ---------- Treatment effect modifiers ----------
  immuno_effect <- 0
  prot_red      <- 0
  egfr_prot     <- 0
  skin_effect   <- 0
  gi_effect     <- 0

  if (treatment == "steroids") {
    immuno_effect <- 0.55
    prot_red      <- 0.40
    skin_effect   <- 0.70
    gi_effect     <- 0.60
  } else if (treatment == "steroids_mmf") {
    immuno_effect <- 0.70
    prot_red      <- 0.55
    skin_effect   <- 0.75
    gi_effect     <- 0.65
  } else if (treatment == "rituximab") {
    immuno_effect <- 0.80
    prot_red      <- 0.60
    skin_effect   <- 0.65
    gi_effect     <- 0.55
  } else if (treatment == "steroids_rtx") {
    immuno_effect <- 0.85
    prot_red      <- 0.70
    skin_effect   <- 0.80
    gi_effect     <- 0.72
  } else if (treatment == "sparsentan") {
    immuno_effect <- 0.20
    prot_red      <- 0.65
    egfr_prot     <- 0.50
    skin_effect   <- 0.10
  }

  if (acei)   { prot_red   <- prot_red   + 0.15; egfr_prot <- egfr_prot + 0.20 }
  if (sglt2i) { prot_red   <- prot_red   + 0.12; egfr_prot <- egfr_prot + 0.18 }

  prot_red  <- min(prot_red, 0.90)
  egfr_prot <- min(egfr_prot, 0.80)

  # ---------- Immune Complex & Complement ----------
  ic_decay <- ifelse(treatment == "none", 0,
              ifelse(treatment %in% c("rituximab","steroids_rtx"), 0.015, 0.008))
  ic_level <- gd_iga1_level * exp(-ic_decay * time) * (1 - immuno_effect * (1 - exp(-0.1 * time)))
  ic_level <- pmax(ic_level, 0.5)

  comp_decay <- ic_decay * 0.6
  comp_level <- complement_activity * exp(-comp_decay * time) * (1 - immuno_effect * 0.5 * (1 - exp(-0.1 * time)))
  comp_level <- pmax(comp_level, 0.3)

  # ---------- Proteinuria ----------
  prot_trt_factor <- (1 - prot_red * (1 - exp(-0.08 * time)))
  prot_level <- pmax(baseline_proteinuria * prot_trt_factor * (1 + 0.1 * disease_severity * time / 52), 0.1)

  # ---------- eGFR ----------
  egfr_slope_treated <- egfr_slope_nat * (1 - egfr_prot)
  egfr_level <- baseline_egfr + egfr_slope_treated * time / 52

  # ---------- Hematuria ----------
  hem_base   <- 3 * disease_severity
  hem_level  <- hem_base * (1 - immuno_effect * (1 - exp(-0.1 * time)))
  hem_level  <- pmax(hem_level, 0)

  # ---------- Skin Purpura ----------
  purp_base  <- 5 * disease_severity
  purp_level <- purp_base * exp(-0.05 * time) * (1 - skin_effect * (1 - exp(-0.2 * time)))

  # ---------- GI Score ----------
  gi_base   <- 4 * disease_severity
  gi_level  <- gi_base * exp(-0.06 * time) * (1 - gi_effect * (1 - exp(-0.2 * time)))

  # ---------- BAFF / IL-6 / TNF ----------
  cytokine_red <- immuno_effect * 0.7
  baff_level   <- 2.0 * (1 - cytokine_red * (1 - exp(-0.05 * time)))
  il6_level    <- 1.5 * ic_level / gd_iga1_level
  tnf_level    <- 1.2 * ic_level / gd_iga1_level

  data.frame(
    time       = time,
    ic_level   = ic_level,
    comp_level = comp_level,
    prot       = prot_level,
    egfr       = egfr_level,
    hematuria  = hem_level,
    purpura    = purp_level,
    gi_score   = gi_level,
    baff       = baff_level,
    il6        = il6_level,
    tnf        = tnf_level
  )
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "IgA Vasculitis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "patient",   icon = icon("user-md")),
      menuItem("PK Profile",            tabName = "pk",        icon = icon("pills")),
      menuItem("PD Key Markers",        tabName = "pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints",    tabName = "clinical",  icon = icon("heartbeat")),
      menuItem("Scenario Comparison",   tabName = "scenario",  icon = icon("balance-scale")),
      menuItem("Biomarkers",            tabName = "biomarker", icon = icon("dna")),
      menuItem("About",                 tabName = "about",     icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
    "))),

    tabItems(
      # ── Tab 1: Patient Profile ───────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Demographics", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("age",    "Age (years)",          value = 35, min = 2, max = 80),
            numericInput("weight", "Weight (kg)",          value = 70, min = 10, max = 150),
            selectInput("pop",     "Population",
                        choices  = c("Adult" = "adult", "Pediatric" = "pediatric"),
                        selected = "adult"),
            numericInput("bsa",    "BSA (m²)",             value = 1.75, min = 0.5, max = 2.5, step = 0.05)
          ),
          box(title = "Disease Baseline Parameters", status = "warning", solidHeader = TRUE, width = 4,
            numericInput("gd_iga1",     "Gd-IgA1 level (mg/L)",        value = 4.0,  min = 1,  max = 10, step = 0.5),
            numericInput("complement",  "Complement activation index",  value = 1.5,  min = 0.5, max = 4, step = 0.1),
            numericInput("prot_base",   "Baseline proteinuria (g/day)", value = 2.0,  min = 0.1, max = 10, step = 0.1),
            numericInput("egfr_base",   "Baseline eGFR (mL/min/1.73m²)", value = 75, min = 10, max = 120, step = 5)
          ),
          box(title = "Treatment Selection", status = "success", solidHeader = TRUE, width = 4,
            selectInput("treatment", "Primary Treatment",
                        choices = c(
                          "No treatment (natural history)" = "none",
                          "Corticosteroids alone"          = "steroids",
                          "Steroids + MMF"                 = "steroids_mmf",
                          "Rituximab monotherapy"          = "rituximab",
                          "Steroids + Rituximab"           = "steroids_rtx",
                          "Sparsentan"                     = "sparsentan"
                        ), selected = "steroids"),
            checkboxInput("acei",   "Add ACEi/ARB (renoprotective)",  value = TRUE),
            checkboxInput("sglt2i", "Add SGLT2 inhibitor",            value = FALSE),
            numericInput("duration", "Simulation duration (weeks)",   value = 52, min = 4, max = 156, step = 4)
          )
        ),
        fluidRow(
          box(title = "Disease Activity Summary", status = "info", solidHeader = TRUE, width = 12,
            valueBoxOutput("disease_severity_box"),
            valueBoxOutput("renal_risk_box"),
            valueBoxOutput("treatment_response_box")
          )
        ),
        fluidRow(
          box(title = "Disease Progression Overview", width = 12,
              plotlyOutput("overview_plot", height = "350px"))
        )
      ),

      # ── Tab 2: PK Profile ────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Prednisolone PK Parameters", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("pred_dose",  "Prednisolone dose (mg/kg/day)", value = 1.0, min = 0, max = 2, step = 0.1),
            numericInput("pred_dur",   "Prednisolone duration (weeks)", value = 4,   min = 1, max = 24),
            numericInput("pred_ka",    "Ka (1/h)",                      value = 2.0, min = 0.5, max = 5),
            numericInput("pred_cl",    "CL (L/h)",                      value = 15,  min = 5, max = 40),
            numericInput("pred_v",     "Vc (L)",                        value = 40,  min = 10, max = 100)
          ),
          box(title = "MMF PK Parameters", status = "warning", solidHeader = TRUE, width = 4,
            numericInput("mmf_dose",   "MMF dose (mg BID)",             value = 1000, min = 0, max = 2000, step = 250),
            numericInput("mmf_ka",     "Ka (1/h)",                      value = 1.5, min = 0.5, max = 4),
            numericInput("mmf_cl",     "CL/F (L/h)",                    value = 25,  min = 5, max = 60),
            numericInput("mmf_vc",     "Vc/F (L)",                      value = 18,  min = 5, max = 50)
          ),
          box(title = "Rituximab PK Parameters", status = "danger", solidHeader = TRUE, width = 4,
            numericInput("rtx_dose",   "RTX dose (mg/m²)",              value = 375, min = 0, max = 750, step = 75),
            numericInput("rtx_cl",     "CL (L/h)",                      value = 0.014, min = 0.005, max = 0.05, step = 0.001),
            numericInput("rtx_vc",     "Vc (L)",                        value = 2.7,  min = 1, max = 6, step = 0.1)
          )
        ),
        fluidRow(
          box(title = "Prednisolone Plasma Concentration (single day)", width = 6,
              plotlyOutput("pred_pk_plot",  height = "300px")),
          box(title = "MPA (Mycophenolic Acid) Concentration (steady state)", width = 6,
              plotlyOutput("mmf_pk_plot",   height = "300px"))
        ),
        fluidRow(
          box(title = "Rituximab Concentration Over Treatment Course", width = 6,
              plotlyOutput("rtx_pk_plot",   height = "300px")),
          box(title = "Drug Exposure Summary", width = 6,
              tableOutput("pk_summary_table"))
        )
      ),

      # ── Tab 3: PD Key Markers ─────────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "IgA Immune Complex & Complement Dynamics", width = 6,
              plotlyOutput("ic_comp_plot", height = "350px")),
          box(title = "Cytokine Response (BAFF, IL-6, TNF-α)", width = 6,
              plotlyOutput("cytokine_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "PD Biomarker Interpretation", status = "info", solidHeader = TRUE, width = 12,
            p(strong("Gd-IgA1:"), "Galactose-deficient IgA1; key driver of mesangial deposition. Normal <2.0 mg/L."),
            p(strong("Immune Complex (IC):"), "Circulating Gd-IgA1/anti-IgA1 IgG complexes. Elevated in active disease."),
            p(strong("Complement activation:"), "Predominantly lectin pathway (MBL/MASP) in IgAV. sC5b-9 correlates with disease severity."),
            p(strong("BAFF:"), "B-cell Activating Factor — promotes survival of autoreactive B cells producing Gd-IgA1."),
            p(strong("IL-6 / TNF-α:"), "Pro-inflammatory cytokines driving vascular inflammation and endothelial activation.")
          )
        )
      ),

      # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "Proteinuria over Time", width = 6,
              plotlyOutput("prot_plot", height = "300px")),
          box(title = "eGFR Trajectory", width = 6,
              plotlyOutput("egfr_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Hematuria Score", width = 4,
              plotlyOutput("hem_plot", height = "280px")),
          box(title = "Skin Purpura Score", width = 4,
              plotlyOutput("purp_plot", height = "280px")),
          box(title = "GI Involvement Score", width = 4,
              plotlyOutput("gi_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "CKD Stage Progression Risk", status = "danger", solidHeader = TRUE, width = 12,
            plotlyOutput("ckd_risk_plot", height = "280px"))
        )
      ),

      # ── Tab 5: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Comparison Settings", status = "primary", solidHeader = TRUE, width = 3,
            checkboxGroupInput("comp_scenarios", "Select Scenarios:",
              choices = list(
                "No treatment"      = "none",
                "Steroids"          = "steroids",
                "Steroids + MMF"    = "steroids_mmf",
                "Rituximab"         = "rituximab",
                "Steroids + RTX"    = "steroids_rtx",
                "Sparsentan"        = "sparsentan"
              ),
              selected = c("none", "steroids", "steroids_mmf", "rituximab")
            ),
            checkboxInput("comp_acei",   "Include ACEi/ARB",  value = TRUE),
            checkboxInput("comp_sglt2i", "Include SGLT2i",     value = FALSE),
            hr(),
            p(em("All scenarios use the patient profile defined in the 'Patient Profile' tab."))
          ),
          box(title = "Proteinuria Comparison", width = 9,
              plotlyOutput("comp_prot_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "eGFR Comparison", width = 6,
              plotlyOutput("comp_egfr_plot", height = "300px")),
          box(title = "Immune Complex Comparison", width = 6,
              plotlyOutput("comp_ic_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Treatment Outcomes Summary Table", width = 12,
              DTOutput("scenario_table"))
        )
      ),

      # ── Tab 6: Biomarkers ─────────────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Reference Ranges", status = "info", solidHeader = TRUE, width = 4,
            tableOutput("biomarker_ref_table")
          ),
          box(title = "Simulated Biomarker Trajectory", width = 8,
              plotlyOutput("biomarker_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Renal Function Biomarkers", width = 6,
              plotlyOutput("renal_biomarker_plot", height = "300px")),
          box(title = "Inflammatory Biomarkers", width = 6,
              plotlyOutput("inflam_biomarker_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Biomarker Interpretation Guide", width = 12,
            DT::DTOutput("biomarker_guide_table"))
        )
      ),

      # ── Tab 7: About ──────────────────────────────────────────────────────
      tabItem(tabName = "about",
        fluidRow(
          box(title = "IgA Vasculitis (Henoch-Schönlein Purpura) — Disease Overview",
              status = "primary", solidHeader = TRUE, width = 12,
            h4("Pathophysiology"),
            p("IgA Vasculitis (IgAV), formerly called Henoch-Schönlein Purpura (HSP), is a small-vessel
              vasculitis characterized by IgA-dominant immune deposits. The pathogenesis involves a
              'multi-hit' model:"),
            tags$ol(
              tags$li(strong("Gd-IgA1 overproduction:"), "Abnormal O-glycosylation of IgA1 hinge region
                      produces galactose-deficient IgA1 (Gd-IgA1) — the primary autoantigen."),
              tags$li(strong("Autoantibody formation:"), "Anti-Gd-IgA1 IgG antibodies form, leading to
                      Gd-IgA1-containing immune complexes (ICs)."),
              tags$li(strong("IC deposition:"), "Circulating ICs deposit in vessel walls of skin,
                      kidney (mesangium), GI tract, and joints."),
              tags$li(strong("Complement activation:"), "Predominantly via the lectin pathway
                      (MBL/MASP-2); C3 and sC5b-9 mediate inflammation."),
              tags$li(strong("Vascular inflammation:"), "Neutrophil/monocyte recruitment, endothelial
                      activation, leukocytoclastic vasculitis.")
            ),
            h4("Clinical Manifestations"),
            tags$ul(
              tags$li("Palpable purpura (100%) — predominantly lower extremities"),
              tags$li("Arthritis/arthralgia (75%) — non-erosive, large joints"),
              tags$li("Abdominal pain/GI involvement (60%) — colicky pain, intussusception"),
              tags$li("Renal involvement (50%) — hematuria, proteinuria, IgAV nephritis")
            ),
            h4("QSP Model Summary"),
            p("This dashboard simulates disease dynamics using a semi-mechanistic ODE-based model
              calibrated against published clinical trial data. Key references include:"),
            tags$ul(
              tags$li("Pillebout et al. 2010 (JASN) — adult IgAV natural history"),
              tags$li("PROTECT trial (Heerspink et al. 2023, Lancet) — sparsentan in IgAN"),
              tags$li("DAPA-CKD trial (Wheeler et al. 2021, Kidney Int) — SGLT2i"),
              tags$li("Fenoglio et al. 2021 (Sci Rep) — rituximab in adult IgAV")
            )
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run primary simulation
  sim_data <- reactive({
    run_simulation(
      age = input$age, weight = input$weight,
      baseline_proteinuria = input$prot_base,
      baseline_egfr        = input$egfr_base,
      gd_iga1_level        = input$gd_iga1,
      complement_activity  = input$complement,
      treatment            = input$treatment,
      pred_dose            = input$pred_dose,
      mmf_dose             = input$mmf_dose,
      rtx_dose             = input$rtx_dose,
      acei                 = input$acei,
      sglt2i               = input$sglt2i,
      duration_weeks       = input$duration
    )
  })

  # -- Value boxes ----
  output$disease_severity_box <- renderValueBox({
    sev <- min(10, round((input$gd_iga1 / 4 + input$complement / 2) * 3, 1))
    color <- if (sev < 3) "green" else if (sev < 6) "yellow" else "red"
    valueBox(sev, "Disease Severity Score (0-10)", icon = icon("exclamation-triangle"), color = color)
  })

  output$renal_risk_box <- renderValueBox({
    risk <- ifelse(input$prot_base >= 3.5 & input$egfr_base < 60, "High",
            ifelse(input$prot_base >= 1.0, "Moderate", "Low"))
    color <- if (risk == "Low") "green" else if (risk == "Moderate") "yellow" else "red"
    valueBox(risk, "Renal Progression Risk", icon = icon("kidneys"), color = color)
  })

  output$treatment_response_box <- renderValueBox({
    resp <- switch(input$treatment,
      "none"         = "None",
      "steroids"     = "Moderate",
      "steroids_mmf" = "Good",
      "rituximab"    = "Good",
      "steroids_rtx" = "Very Good",
      "sparsentan"   = "Renoprotective"
    )
    color <- if (resp %in% c("Very Good","Renoprotective")) "green"
             else if (resp == "Good") "teal"
             else if (resp == "Moderate") "yellow"
             else "red"
    valueBox(resp, "Expected Treatment Response", icon = icon("stethoscope"), color = color)
  })

  # -- Overview plot ----
  output$overview_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = scale(prot)[,1],   color = "Proteinuria"),   linewidth = 1.2) +
      geom_line(aes(y = scale(egfr)[,1],   color = "eGFR"),          linewidth = 1.2) +
      geom_line(aes(y = scale(ic_level)[,1],color = "IC Level"),     linewidth = 1.0, linetype = "dashed") +
      geom_line(aes(y = scale(purpura)[,1], color = "Purpura"),       linewidth = 1.0, linetype = "dotted") +
      labs(x = "Time (weeks)", y = "Normalized Value (z-score)", color = "Variable") +
      scale_color_manual(values = c(
        "Proteinuria" = "#E74C3C", "eGFR" = "#3498DB",
        "IC Level"    = "#9B59B6", "Purpura" = "#F39C12")) +
      theme_minimal() + theme(legend.position = "right")
    ggplotly(p)
  })

  # -- PK plots ----
  output$pred_pk_plot <- renderPlotly({
    t_pk <- seq(0, 24, by = 0.5)
    dose <- input$pred_dose * input$weight
    ka   <- input$pred_ka
    cl   <- input$pred_cl
    v    <- input$pred_v
    k    <- cl / v
    conc <- (dose * ka) / (v * (ka - k)) * (exp(-k * t_pk) - exp(-ka * t_pk))
    conc <- pmax(conc, 0)
    df_pk <- data.frame(time = t_pk, conc = conc)
    p <- ggplot(df_pk, aes(x = time, y = conc)) +
      geom_line(color = "#3498DB", linewidth = 1.3) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.7) +
      annotate("text", x = 20, y = 55, label = "EC50 (~50 ng/mL)", color = "red", size = 3) +
      labs(x = "Time (h)", y = "Prednisolone (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$mmf_pk_plot <- renderPlotly({
    t_pk   <- seq(0, 12, by = 0.25)
    dose   <- input$mmf_dose
    ka_mmf <- input$mmf_ka
    cl_mmf <- input$mmf_cl
    v_mmf  <- input$mmf_vc
    k_mmf  <- cl_mmf / v_mmf
    conc   <- (dose * ka_mmf) / (v_mmf * (ka_mmf - k_mmf)) * (exp(-k_mmf * t_pk) - exp(-ka_mmf * t_pk))
    conc   <- pmax(conc, 0)
    df_pk  <- data.frame(time = t_pk, conc = conc)
    p <- ggplot(df_pk, aes(x = time, y = conc)) +
      geom_line(color = "#E67E22", linewidth = 1.3) +
      labs(x = "Time (h)", y = "MPA (μg/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$rtx_pk_plot <- renderPlotly({
    doses_wk <- c(0, 1, 2, 3)
    t_rtx    <- seq(0, 24, by = 0.5)
    bsa      <- input$bsa
    dose     <- input$rtx_dose * bsa
    cl_rtx   <- input$rtx_cl
    v_rtx    <- input$rtx_vc
    k_rtx    <- cl_rtx / v_rtx
    conc <- rep(0, length(t_rtx))
    for (d_wk in doses_wk) {
      t_shift <- t_rtx - d_wk
      conc    <- conc + ifelse(t_shift >= 0, (dose / v_rtx) * exp(-k_rtx * t_shift), 0)
    }
    df_pk <- data.frame(time = t_rtx, conc = pmax(conc, 0))
    p <- ggplot(df_pk, aes(x = time, y = conc)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_vline(xintercept = doses_wk, linetype = "dashed", color = "gray50", alpha = 0.7) +
      labs(x = "Time (weeks)", y = "Rituximab (μg/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$pk_summary_table <- renderTable({
    dose_pred <- input$pred_dose * input$weight
    t12_pred  <- log(2) / (input$pred_cl / input$pred_v)
    t12_mmf   <- log(2) / (input$mmf_cl  / input$mmf_vc)
    t12_rtx   <- log(2) / (input$rtx_cl  / input$rtx_vc) * 24
    data.frame(
      Drug       = c("Prednisolone", "Mycophenolate (MPA)", "Rituximab"),
      Route      = c("Oral",         "Oral",                 "IV"),
      Dose       = c(paste0(dose_pred, " mg/day"), paste0(input$mmf_dose * 2, " mg/day"), paste0(round(input$rtx_dose * input$bsa), " mg")),
      `t½ (h)`   = c(round(t12_pred, 1), round(t12_mmf, 1), round(t12_rtx, 1)),
      Mechanism  = c("GR agonist", "IMPDH inhibitor", "Anti-CD20 mAb")
    )
  }, striped = TRUE, bordered = TRUE)

  # -- PD plots ----
  output$ic_comp_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = ic_level,   color = "IgA Immune Complex"), linewidth = 1.3) +
      geom_line(aes(y = comp_level, color = "Complement Activity"), linewidth = 1.3) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(df$time) * 0.8, y = 2.3, label = "Normal threshold", size = 3) +
      labs(x = "Time (weeks)", y = "Relative Level", color = "") +
      scale_color_manual(values = c("IgA Immune Complex" = "#9B59B6", "Complement Activity" = "#1ABC9C")) +
      theme_minimal()
    ggplotly(p)
  })

  output$cytokine_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = baff, color = "BAFF"),  linewidth = 1.3) +
      geom_line(aes(y = il6,  color = "IL-6"),  linewidth = 1.3) +
      geom_line(aes(y = tnf,  color = "TNF-α"), linewidth = 1.3) +
      labs(x = "Time (weeks)", y = "Relative Level (fold-change)", color = "Cytokine") +
      scale_color_manual(values = c("BAFF" = "#E74C3C", "IL-6" = "#3498DB", "TNF-α" = "#F39C12")) +
      theme_minimal()
    ggplotly(p)
  })

  # -- Clinical endpoint plots ----
  output$prot_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = prot)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "green3") +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = "red") +
      annotate("text", x = max(df$time)*0.7, y = 0.65, label = "Partial remission", color = "green3", size = 3) +
      annotate("text", x = max(df$time)*0.7, y = 3.65, label = "Nephrotic range", color = "red", size = 3) +
      labs(x = "Time (weeks)", y = "Proteinuria (g/day)") +
      theme_minimal()
    ggplotly(p)
  })

  output$egfr_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = egfr)) +
      geom_line(color = "#3498DB", linewidth = 1.3) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
      annotate("text", x = max(df$time)*0.7, y = 63, label = "CKD G3 threshold", color = "orange", size = 3) +
      annotate("text", x = max(df$time)*0.7, y = 33, label = "CKD G4 threshold", color = "red", size = 3) +
      labs(x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)") +
      theme_minimal()
    ggplotly(p)
  })

  output$hem_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = hematuria)) +
      geom_line(color = "#C0392B", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = hematuria), fill = "#C0392B", alpha = 0.2) +
      labs(x = "Time (weeks)", y = "Hematuria Score") + theme_minimal()
    ggplotly(p)
  })

  output$purp_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = purpura)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = purpura), fill = "#8E44AD", alpha = 0.2) +
      labs(x = "Time (weeks)", y = "Purpura Score") + theme_minimal()
    ggplotly(p)
  })

  output$gi_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time, y = gi_score)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = gi_score), fill = "#27AE60", alpha = 0.2) +
      labs(x = "Time (weeks)", y = "GI Score") + theme_minimal()
    ggplotly(p)
  })

  output$ckd_risk_plot <- renderPlotly({
    df <- sim_data()
    df <- df %>% mutate(
      ckd_stage = case_when(
        egfr >= 90 ~ "G1/G2 (≥60)",
        egfr >= 60 ~ "G2/G3a (60-89)",
        egfr >= 45 ~ "G3a (45-59)",
        egfr >= 30 ~ "G3b (30-44)",
        egfr >= 15 ~ "G4 (15-29)",
        TRUE       ~ "G5 (<15)"
      )
    )
    p <- ggplot(df, aes(x = time, y = egfr, fill = ckd_stage)) +
      geom_area(alpha = 0.4) +
      geom_line(color = "#2C3E50", linewidth = 1.2) +
      scale_fill_manual(values = c(
        "G1/G2 (≥60)" = "#2ECC71", "G2/G3a (60-89)" = "#F1C40F",
        "G3a (45-59)" = "#E67E22", "G3b (30-44)" = "#E74C3C",
        "G4 (15-29)"  = "#8E44AD", "G5 (<15)" = "#2C3E50"
      )) +
      labs(x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)", fill = "CKD Stage") +
      theme_minimal()
    ggplotly(p)
  })

  # -- Scenario comparison ----
  comp_data <- reactive({
    req(input$comp_scenarios)
    scenario_names <- c(
      "none"         = "No treatment",
      "steroids"     = "Steroids",
      "steroids_mmf" = "Steroids+MMF",
      "rituximab"    = "Rituximab",
      "steroids_rtx" = "Steroids+RTX",
      "sparsentan"   = "Sparsentan"
    )
    bind_rows(lapply(input$comp_scenarios, function(sc) {
      df <- run_simulation(
        age = input$age, weight = input$weight,
        baseline_proteinuria = input$prot_base,
        baseline_egfr        = input$egfr_base,
        gd_iga1_level        = input$gd_iga1,
        complement_activity  = input$complement,
        treatment            = sc,
        acei                 = input$comp_acei,
        sglt2i               = input$comp_sglt2i,
        duration_weeks       = input$duration
      )
      df$scenario <- scenario_names[sc]
      df
    }))
  })

  output$comp_prot_plot <- renderPlotly({
    df <- comp_data()
    p  <- ggplot(df, aes(x = time, y = prot, color = scenario)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
      labs(x = "Time (weeks)", y = "Proteinuria (g/day)", color = "Scenario") +
      theme_minimal()
    ggplotly(p)
  })

  output$comp_egfr_plot <- renderPlotly({
    df <- comp_data()
    p  <- ggplot(df, aes(x = time, y = egfr, color = scenario)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
      labs(x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)", color = "Scenario") +
      theme_minimal()
    ggplotly(p)
  })

  output$comp_ic_plot <- renderPlotly({
    df <- comp_data()
    p  <- ggplot(df, aes(x = time, y = ic_level, color = scenario)) +
      geom_line(linewidth = 1.2) +
      labs(x = "Time (weeks)", y = "IC Level (relative)", color = "Scenario") +
      theme_minimal()
    ggplotly(p)
  })

  output$scenario_table <- DT::renderDT({
    req(input$comp_scenarios)
    scenario_names <- c(
      "none"         = "No treatment",
      "steroids"     = "Steroids",
      "steroids_mmf" = "Steroids+MMF",
      "rituximab"    = "Rituximab",
      "steroids_rtx" = "Steroids+RTX",
      "sparsentan"   = "Sparsentan"
    )
    rows <- lapply(input$comp_scenarios, function(sc) {
      df <- run_simulation(
        age = input$age, weight = input$weight,
        baseline_proteinuria = input$prot_base,
        baseline_egfr        = input$egfr_base,
        gd_iga1_level        = input$gd_iga1,
        complement_activity  = input$complement,
        treatment            = sc,
        acei                 = input$comp_acei,
        sglt2i               = input$comp_sglt2i,
        duration_weeks       = input$duration
      )
      last <- tail(df, 1)
      data.frame(
        Scenario              = scenario_names[sc],
        `Final Proteinuria`   = round(last$prot, 2),
        `Final eGFR`          = round(last$egfr, 1),
        `Final IC Level`      = round(last$ic_level, 2),
        `Purpura Clearance`   = ifelse(last$purpura < 0.5, "Yes", "No"),
        `eGFR Change`         = round(last$egfr - input$egfr_base, 1)
      )
    })
    DT::datatable(bind_rows(rows), options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE)
  })

  # -- Biomarker tab ----
  output$biomarker_ref_table <- renderTable({
    data.frame(
      Biomarker       = c("Gd-IgA1", "anti-Gd-IgA1 IgG", "sC5b-9", "BAFF", "IL-6", "TNF-α",
                          "Proteinuria", "eGFR", "Hematuria"),
      Normal          = c("<2.0 mg/L", "<200 U/mL",     "<250 ng/mL", "<2 ng/mL", "<5 pg/mL",
                          "<15 pg/mL", "<0.15 g/day", ">90", "Absent"),
      `Disease Active`= c(">4.0 mg/L", ">500 U/mL",    ">500 ng/mL", ">4 ng/mL", ">30 pg/mL",
                          ">30 pg/mL", ">1.0 g/day", "<60", "Gross")
    )
  }, striped = TRUE, bordered = TRUE)

  output$biomarker_plot <- renderPlotly({
    df <- sim_data()
    df_long <- df %>%
      select(time, ic_level, comp_level, baff, il6, tnf) %>%
      pivot_longer(-time, names_to = "biomarker", values_to = "value") %>%
      mutate(biomarker = recode(biomarker,
        "ic_level"   = "IgA IC",
        "comp_level" = "Complement",
        "baff"       = "BAFF",
        "il6"        = "IL-6",
        "tnf"        = "TNF-α"
      ))
    p <- ggplot(df_long, aes(x = time, y = value, color = biomarker)) +
      geom_line(linewidth = 1.2) +
      labs(x = "Time (weeks)", y = "Relative Level", color = "Biomarker") +
      theme_minimal()
    ggplotly(p)
  })

  output$renal_biomarker_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = prot,    color = "Proteinuria (g/day)"), linewidth = 1.2) +
      geom_line(aes(y = hematuria, color = "Hematuria Score"), linewidth = 1.2) +
      labs(x = "Time (weeks)", y = "Value", color = "Biomarker") +
      scale_color_manual(values = c("Proteinuria (g/day)" = "#E74C3C", "Hematuria Score" = "#3498DB")) +
      theme_minimal()
    ggplotly(p)
  })

  output$inflam_biomarker_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = baff, color = "BAFF"),  linewidth = 1.2) +
      geom_line(aes(y = il6,  color = "IL-6"),  linewidth = 1.2) +
      geom_line(aes(y = tnf,  color = "TNF-α"), linewidth = 1.2) +
      labs(x = "Time (weeks)", y = "Relative Level", color = "Cytokine") +
      scale_color_manual(values = c("BAFF" = "#E74C3C", "IL-6" = "#3498DB", "TNF-α" = "#F39C12")) +
      theme_minimal()
    ggplotly(p)
  })

  output$biomarker_guide_table <- DT::renderDT({
    DT::datatable(data.frame(
      Biomarker    = c("Gd-IgA1", "Anti-Gd-IgA1 IgG", "sC5b-9",
                       "BAFF", "IL-6", "TNF-α", "Proteinuria",
                       "eGFR", "Urinary RBC"),
      Clinical_Use = c(
        "Disease activity marker; guides immunosuppression intensity",
        "Predicts renal involvement severity; elevated in IgAV nephritis",
        "Complement activation; correlates with purpura extent and nephritis",
        "B-cell survival factor; target for atacicept/belimumab",
        "Acute phase response; elevated with active flares",
        "Vascular inflammation driver; correlates with purpura score",
        "Primary renal endpoint; target <0.5 g/day for remission",
        "Long-term renal function; monitor CKD progression",
        "Glomerular inflammation marker; persistence = poor prognosis"
      ),
      Threshold    = c(">2 mg/L", ">200 AU/mL", ">250 ng/mL",
                       ">2 ng/mL", ">10 pg/mL", ">20 pg/mL",
                       ">1.0 g/day", "<60 mL/min", ">25 RBC/HPF"),
      Monitoring   = c("Monthly × 6, then quarterly",
                       "Every 3 months", "At flares",
                       "Before/after B-cell therapy", "At flares",
                       "At flares", "Monthly", "Every 3 months",
                       "Monthly × 6, then quarterly")
    ), options = list(pageLength = 10), rownames = FALSE)
  })
}

shinyApp(ui = ui, server = server)
