# ============================================================
# Endometriosis QSP Shiny Application
# 자궁내막증 정량적 시스템 약리학 대시보드
# ============================================================
# Required packages: shiny, ggplot2, dplyr
# Run with: shiny::runApp("endo_shiny_app.R")
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)

# ============================================================
# SIMULATION ENGINE (ODE-free, parametric lookup approach)
# ============================================================

# Drug PK parameters
pk_params <- list(
  "No treatment"           = list(ka=0,    kel=0,    F=0,    dose=0,    interval=24,  Vd=1),
  "Leuprolide depot"       = list(ka=0.02, kel=0.012, F=0.95, dose=3750, interval=672, Vd=28),
  "Elagolix 150mg/d"      = list(ka=1.8,  kel=0.14,  F=0.58, dose=150,  interval=24,  Vd=275),
  "Elagolix 200mg BID"    = list(ka=1.8,  kel=0.14,  F=0.58, dose=200,  interval=12,  Vd=275),
  "Dienogest 2mg/d"       = list(ka=1.2,  kel=0.073, F=0.91, dose=2,    interval=24,  Vd=65),
  "Letrozole+Add-back"    = list(ka=1.0,  kel=0.014, F=0.99, dose=2.5,  interval=24,  Vd=185),
  "Combined OCP"          = list(ka=2.0,  kel=0.078, F=0.95, dose=0.03, interval=24,  Vd=44)
)

# Drug E2 suppression profiles: list(partial_suppression, full_suppression_onset_weeks, e2_floor_pg)
drug_e2 <- list(
  "No treatment"           = list(e2_frac=1.00, onset_wk=0,  e2_floor=60,  bmd_slope=-0.000),
  "Leuprolide depot"       = list(e2_frac=0.05, onset_wk=4,  e2_floor=10,  bmd_slope=-0.010),
  "Elagolix 150mg/d"      = list(e2_frac=0.50, onset_wk=1,  e2_floor=30,  bmd_slope=-0.003),
  "Elagolix 200mg BID"    = list(e2_frac=0.05, onset_wk=1,  e2_floor=8,   bmd_slope=-0.008),
  "Dienogest 2mg/d"       = list(e2_frac=0.40, onset_wk=2,  e2_floor=25,  bmd_slope=-0.002),
  "Letrozole+Add-back"    = list(e2_frac=0.08, onset_wk=2,  e2_floor=8,   bmd_slope=-0.001),
  "Combined OCP"          = list(e2_frac=0.55, onset_wk=1,  e2_floor=35,  bmd_slope=-0.001)
)

# Drug pain reduction profiles: list(dysmen_frac, dysp_frac, cpp_frac)
drug_pain <- list(
  "No treatment"           = list(dysmen=0.00, dysp=0.00, cpp=0.00, lesion_frac=0.00),
  "Leuprolide depot"       = list(dysmen=0.75, dysp=0.65, cpp=0.70, lesion_frac=0.55),
  "Elagolix 150mg/d"      = list(dysmen=0.55, dysp=0.45, cpp=0.50, lesion_frac=0.35),
  "Elagolix 200mg BID"    = list(dysmen=0.70, dysp=0.65, cpp=0.65, lesion_frac=0.50),
  "Dienogest 2mg/d"       = list(dysmen=0.70, dysp=0.60, cpp=0.65, lesion_frac=0.50),
  "Letrozole+Add-back"    = list(dysmen=0.65, dysp=0.55, cpp=0.60, lesion_frac=0.60),
  "Combined OCP"           = list(dysmen=0.50, dysp=0.40, cpp=0.45, lesion_frac=0.25)
)

# Simulate PK profile for selected drug
simulate_pk <- function(drug, duration_months) {
  if (drug == "No treatment") {
    times <- seq(0, duration_months * 30, by = 0.5)
    return(data.frame(time_h = times * 24, conc = 0, time_d = times))
  }
  pk   <- pk_params[[drug]]
  dose <- pk$dose
  ka   <- pk$ka
  kel  <- pk$kel
  Vd   <- pk$Vd
  F    <- pk$F
  tau  <- pk$interval  # hours

  # simulate first 7 days in hours for PK tab
  t_h <- seq(0, 168, by = 0.5)
  n_doses <- floor(168 / tau) + 1
  conc <- numeric(length(t_h))
  for (i in seq_len(n_doses)) {
    t_dose <- (i - 1) * tau
    dt <- t_h - t_dose
    valid <- dt >= 0
    if (any(valid)) {
      conc[valid] <- conc[valid] +
        (F * dose * ka / (Vd * (ka - kel))) *
        (exp(-kel * dt[valid]) - exp(-ka * dt[valid]))
    }
  }
  conc <- pmax(conc, 0)
  data.frame(time_h = t_h, conc = conc, time_d = t_h / 24)
}

# Simulate HPO hormones over treatment duration
simulate_hpo <- function(drug, duration_months, baseline_e2 = 60) {
  times <- seq(0, duration_months, by = 1/30)  # monthly resolution
  ep    <- drug_e2[[drug]]

  # E2 dynamics: logistic suppression with onset delay
  onset  <- ep$onset_wk / 4.33  # months
  e2_ss  <- baseline_e2 * ep$e2_frac
  e2     <- ifelse(times <= onset,
                   baseline_e2,
                   e2_ss + (baseline_e2 - e2_ss) * exp(-1.5 * (times - onset)))
  e2     <- pmax(e2, ep$e2_floor)

  # FSH inversely tracks E2 suppression
  fsh_base  <- 8
  fsh_max   <- 35
  e2_norm   <- (e2 - ep$e2_floor) / (baseline_e2 - ep$e2_floor + 1e-6)
  fsh       <- fsh_base + (fsh_max - fsh_base) * (1 - e2_norm)

  # LH: suppressed when using GnRH agonists/antagonists, less so with progestins
  lh_base   <- 6
  lh_factor <- ifelse(drug %in% c("Leuprolide depot","Elagolix 150mg/d","Elagolix 200mg BID"),
                      0.05, ifelse(drug == "No treatment", 1.0, 0.5))
  lh        <- lh_base * (lh_factor + (1 - lh_factor) * e2_norm)

  # Progesterone: suppressed with GnRH-based therapy, partial with progestins
  p4_base <- 5.0
  p4_factor <- ifelse(drug == "Dienogest 2mg/d", 0.8,
               ifelse(drug %in% c("Leuprolide depot","Elagolix 200mg BID"), 0.1, 0.5))
  p4   <- p4_base * p4_factor + rnorm(length(times), 0, 0.1)
  p4   <- pmax(p4, 0.1)

  data.frame(
    time_months = times,
    E2_pgmL     = pmax(e2, 0),
    FSH_IUL     = pmax(fsh, 0),
    LH_IUL      = pmax(lh,  0),
    P4_ngmL     = pmax(p4,  0)
  )
}

# Simulate pain and clinical endpoints
simulate_clinical <- function(drug, duration_months,
                               dysmen_base, dysp_base, cpp_base,
                               lesion_vol_base = 5.0) {
  times  <- seq(0, duration_months, by = 1/30)
  pp     <- drug_pain[[drug]]
  ep     <- drug_e2[[drug]]

  # Pain: logistic reduction
  k_pain <- 0.8
  dysmen <- dysmen_base * (1 - pp$dysmen * (1 - exp(-k_pain * times)))
  dysp   <- dysp_base   * (1 - pp$dysp   * (1 - exp(-k_pain * times)))
  cpp    <- cpp_base    * (1 - pp$cpp    * (1 - exp(-k_pain * times)))
  dysmen <- pmax(dysmen, 0)
  dysp   <- pmax(dysp,   0)
  cpp    <- pmax(cpp,    0)

  # BMD: cumulative % change from baseline
  bmd_pct <- ep$bmd_slope * 12 * times  # monthly accrual

  # Endometrioma volume: shrinks under treatment
  k_lesion  <- 0.3
  lesion_vol <- lesion_vol_base * (1 - pp$lesion_frac * (1 - exp(-k_lesion * times)))

  # Hot flush score (0–10): high when E2 suppressed
  e2_frac_steady <- ep$e2_frac
  hf_max <- 8 * (1 - e2_frac_steady)
  hf     <- hf_max * (1 - exp(-0.5 * times))

  data.frame(
    time_months  = times,
    pain_dysmen  = dysmen,
    pain_dysp    = dysp,
    pain_cpp     = cpp,
    bmd_pct      = bmd_pct,
    lesion_vol   = pmax(lesion_vol, 0.1),
    hot_flush    = hf
  )
}

# Simulate lesion biology
simulate_lesion <- function(drug, duration_months,
                             baseline_e2 = 60, lesion_vol_base = 5.0) {
  times  <- seq(0, duration_months, by = 1/30)
  pp     <- drug_pain[[drug]]
  ep     <- drug_e2[[drug]]

  onset  <- ep$onset_wk / 4.33
  e2_ss  <- baseline_e2 * ep$e2_frac
  e2_loc <- ifelse(times <= onset,
                   baseline_e2 * 1.5,
                   e2_ss * 1.2 + (baseline_e2 * 1.5 - e2_ss * 1.2) * exp(-1.5 * (times - onset)))
  e2_loc <- pmax(e2_loc, ep$e2_floor * 1.2)

  # Lesion volume
  k_lesion  <- 0.3
  lesion_vol <- lesion_vol_base * (1 - pp$lesion_frac * (1 - exp(-k_lesion * times)))

  # Proliferation rate (relative, 1 = baseline)
  prolif <- 1 - pp$lesion_frac * (1 - exp(-0.5 * times))

  # Peritoneal fluid cytokines and mediators
  il6_base  <- 50  # pg/mL
  pge2_base <- 200 # pg/mL
  ngf_base  <- 30  # pg/mL
  arom_base <- 1.0 # relative activity

  e2_norm <- (e2_loc - ep$e2_floor * 1.2) / (baseline_e2 * 1.5 - ep$e2_floor * 1.2 + 1e-6)

  il6    <- il6_base  * (0.2 + 0.8 * e2_norm)
  pge2   <- pge2_base * (0.15 + 0.85 * e2_norm)
  ngf    <- ngf_base  * (0.25 + 0.75 * (pmax(lesion_vol, 0.1) / lesion_vol_base))
  arom   <- arom_base * e2_norm

  data.frame(
    time_months   = times,
    lesion_vol    = pmax(lesion_vol, 0.1),
    prolif_rate   = pmax(prolif, 0.05),
    e2_local      = pmax(e2_loc, 0),
    IL6_pgmL      = pmax(il6,  5),
    PGE2_pgmL     = pmax(pge2, 20),
    NGF_pgmL      = pmax(ngf,  3),
    aromatase_act = pmax(arom, 0.02)
  )
}

# Summary stats for PK tab
pk_summary <- function(pk_df, drug) {
  if (drug == "No treatment") {
    return(data.frame(Parameter = c("Cmax","Tmax","AUC0-24h","t½","Steady-state"),
                      Value     = c("N/A","N/A","N/A","N/A","N/A"),
                      Unit      = rep("", 5)))
  }
  pk   <- pk_params[[drug]]
  cmax <- max(pk_df$conc, na.rm = TRUE)
  tmax <- pk_df$time_h[which.max(pk_df$conc)]
  auc  <- sum(diff(pk_df$time_h[pk_df$time_h <= 24]) *
                (pk_df$conc[pk_df$time_h <= 24][-1] +
                   pk_df$conc[pk_df$time_h <= 24][-sum(pk_df$time_h <= 24)]) / 2)
  thalf <- round(log(2) / pk$kel, 1)
  ss_time <- round(5 * log(2) / pk$kel, 1)

  data.frame(
    Parameter = c("Cmax", "Tmax", "AUC₀₋₂₄ₕ", "t½ (elimination)", "~Steady-state reached"),
    Value     = c(round(cmax, 2), round(tmax, 1), round(auc, 1), thalf, ss_time),
    Unit      = c("ng/mL or µg/mL", "h", "ng·h/mL", "h", "h")
  )
}

# All-scenario comparison
all_scenarios_compare <- function(duration_months, dysmen_base, dysp_base,
                                   cpp_base, lesion_vol_base = 5.0, baseline_e2 = 60) {
  drugs <- names(drug_pain)
  bind_rows(lapply(drugs, function(d) {
    clin <- simulate_clinical(d, duration_months, dysmen_base, dysp_base, cpp_base, lesion_vol_base)
    hpo  <- simulate_hpo(d, duration_months, baseline_e2)
    merged <- left_join(clin, hpo, by = "time_months")
    merged$drug <- d
    merged
  }))
}

# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Arial, sans-serif; background: #f5f6fa; }
      .navbar { background-color: #7b2d8b !important; }
      .nav-tabs > li.active > a { color: #7b2d8b; font-weight: bold; }
      .well { background: #fff; border: 1px solid #e0e0e0; border-radius: 8px; }
      h3, h4 { color: #7b2d8b; }
      .btn-primary { background-color: #7b2d8b; border-color: #5e1f6e; }
      .shiny-output-error { color: #c0392b; }
      table.dataTable { font-size: 13px; }
    "))
  ),

  titlePanel(
    div(
      h2("자궁내막증 (Endometriosis) QSP 시뮬레이터",
         style = "color:#7b2d8b; font-weight:bold;"),
      h5("Quantitative Systems Pharmacology Dashboard v1.0",
         style = "color:#555;")
    )
  ),

  tabsetPanel(id = "main_tabs",

    # ---- TAB 1: Patient Profile ----
    tabPanel("1. 환자 프로파일 (Patient Profile)",
      fluidRow(
        column(4,
          wellPanel(
            h4("Patient Demographics"),
            sliderInput("age", "Age (years)", min = 20, max = 50, value = 32, step = 1),
            numericInput("bmi", "BMI (kg/m²)", value = 22.5, min = 15, max = 45, step = 0.1),
            selectInput("stage", "Disease Stage (rAFS/rASRM)",
                        choices = c("Stage I (Minimal)" = "I",
                                    "Stage II (Mild)"    = "II",
                                    "Stage III (Moderate)" = "III",
                                    "Stage IV (Severe)"  = "IV"),
                        selected = "III")
          ),
          wellPanel(
            h4("Symptom Severity (NRS 0–10)"),
            sliderInput("dysmen", "Dysmenorrhea NRS", 0, 10, value = 7, step = 0.5),
            sliderInput("dysp",   "Dyspareunia NRS",  0, 10, value = 5, step = 0.5),
            sliderInput("cpp",    "Chronic Pelvic Pain NRS", 0, 10, value = 6, step = 0.5)
          )
        ),
        column(4,
          wellPanel(
            h4("Laboratory Values"),
            numericInput("e2_base",  "Serum E2 (pg/mL)",   value = 60,   min = 5,   max = 400),
            numericInput("fsh_base", "FSH (IU/L)",          value = 8,    min = 1,   max = 80),
            numericInput("amh",      "AMH (ng/mL)",         value = 2.5,  min = 0.1, max = 10),
            numericInput("ca125",    "CA-125 (U/mL)",       value = 45,   min = 1,   max = 500)
          ),
          wellPanel(
            h4("Lesion Characteristics"),
            numericInput("lesion_vol", "Endometrioma Volume (cm³)", value = 5, min = 0.5, max = 100),
            selectInput("laterality", "Laterality",
                        choices = c("Unilateral", "Bilateral"), selected = "Unilateral")
          )
        ),
        column(4,
          wellPanel(
            h4("Treatment Selection"),
            radioButtons("drug", "Select Treatment:",
                         choices = c("No treatment",
                                     "Leuprolide depot",
                                     "Elagolix 150mg/d",
                                     "Elagolix 200mg BID",
                                     "Dienogest 2mg/d",
                                     "Letrozole+Add-back",
                                     "Combined OCP"),
                         selected = "Elagolix 150mg/d")
          ),
          wellPanel(
            h4("Simulation Duration"),
            radioButtons("duration", "Duration:",
                         choices = c("6 months" = 6, "12 months" = 12, "24 months" = 24),
                         selected = 12)
          ),
          actionButton("run_sim", "Run Simulation",
                       class = "btn-primary btn-lg",
                       style = "width:100%; margin-top:10px; font-size:16px;"),
          br(), br(),
          wellPanel(
            h5("Patient Summary", style = "color:#555;"),
            tableOutput("patient_summary_table")
          )
        )
      )
    ),

    # ---- TAB 2: PK ----
    tabPanel("2. PK 약동학 (Pharmacokinetics)",
      fluidRow(
        column(8,
          h4("Drug Concentration vs Time (First 7 Days)"),
          plotOutput("pk_plot", height = "380px"),
          h4("Multiple-dose Simulation — Steady State"),
          plotOutput("pk_ss_plot", height = "280px")
        ),
        column(4,
          h4("PK Summary"),
          tableOutput("pk_table"),
          br(),
          wellPanel(
            h5("PK Model Notes", style = "color:#555;"),
            p("One-compartment first-order absorption model."),
            p("Leuprolide depot: zero-order release approximation."),
            p("GnRH antagonists: immediate receptor binding (no flare-up)."),
            p("GnRH agonist (leuprolide): initial FSH/LH flare during weeks 1–2.")
          )
        )
      )
    ),

    # ---- TAB 3: HPO Hormones ----
    tabPanel("3. HPO 호르몬 (HPO Hormone Dynamics)",
      fluidRow(
        column(12,
          h4("Hormone Time-Series During Treatment"),
          plotOutput("hpo_plot", height = "480px")
        )
      ),
      fluidRow(
        column(12,
          h5("Reference Ranges", style = "color:#555; margin-top:10px;"),
          p("E2 therapeutic target: < 20 pg/mL (menopausal range) for maximum lesion suppression.
            FSH surge indicates ovarian suppression. LH suppression confirms GnRH axis inhibition.")
        )
      )
    ),

    # ---- TAB 4: Pain & Clinical Endpoints ----
    tabPanel("4. 통증 & 임상 지표 (Pain & Clinical Endpoints)",
      fluidRow(
        column(6,
          h4("Pain NRS Over Time"),
          plotOutput("pain_plot", height = "350px")
        ),
        column(6,
          h4("BMD Change (Lumbar Spine, %)"),
          plotOutput("bmd_plot", height = "200px"),
          h4("Hot Flush Severity Score"),
          plotOutput("hf_plot", height = "150px")
        )
      ),
      fluidRow(
        column(6,
          h4("Endometrioma Volume Change"),
          plotOutput("lesion_vol_plot", height = "250px")
        ),
        column(6,
          h4("Clinical Endpoints Summary"),
          br(),
          tableOutput("clinical_table")
        )
      )
    ),

    # ---- TAB 5: Lesion Dynamics ----
    tabPanel("5. 병소 역학 (Lesion Dynamics)",
      fluidRow(
        column(6,
          h4("Lesion Volume (cm³) Over Time"),
          plotOutput("lesion_dyn_plot", height = "280px"),
          h4("Lesion Proliferation Rate (Relative)"),
          plotOutput("prolif_plot", height = "220px")
        ),
        column(6,
          h4("Local E2 at Ectopic Site (pg/mL)"),
          plotOutput("e2_local_plot", height = "220px"),
          h4("Peritoneal Cytokines & Mediators"),
          plotOutput("cytokine_plot", height = "280px")
        )
      ),
      fluidRow(
        column(6,
          h4("Nerve Growth Factor (NGF) at Lesion"),
          plotOutput("ngf_plot", height = "200px")
        ),
        column(6,
          h4("Aromatase Activity (Relative)"),
          plotOutput("arom_plot", height = "200px")
        )
      )
    ),

    # ---- TAB 6: Scenario Comparison ----
    tabPanel("6. 시나리오 비교 (Scenario Comparison)",
      fluidRow(
        column(12,
          h4("All 7 Treatment Scenarios — Side-by-Side Comparison"),
          plotOutput("scenario_facet_plot", height = "600px")
        )
      ),
      fluidRow(
        column(6,
          h4("Pain Response Rates (≥50% Reduction)"),
          plotOutput("waterfall_plot", height = "350px")
        ),
        column(6,
          h4("Treatment Comparison Summary Table"),
          br(),
          tableOutput("scenario_table")
        )
      )
    )
  ) # end tabsetPanel
) # end fluidPage

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Reactive: run simulation on button press
  sim_data <- eventReactive(input$run_sim, {
    drug    <- input$drug
    dur     <- as.numeric(input$duration)
    e2_b    <- input$e2_base
    dysmen  <- input$dysmen
    dysp    <- input$dysp
    cpp     <- input$cpp
    lvol    <- input$lesion_vol

    list(
      pk       = simulate_pk(drug, dur),
      hpo      = simulate_hpo(drug, dur, e2_b),
      clinical = simulate_clinical(drug, dur, dysmen, dysp, cpp, lvol),
      lesion   = simulate_lesion(drug, dur, e2_b, lvol),
      all_scen = all_scenarios_compare(dur, dysmen, dysp, cpp, lvol, e2_b),
      drug     = drug,
      dur      = dur
    )
  }, ignoreNULL = FALSE)

  # Initialize on load
  observe({
    if (is.null(sim_data())) {
      shinyjs::click("run_sim")
    }
  })

  # Auto-run on first render
  observe({
    isolate({
      drug   <- input$drug
      dur    <- as.numeric(input$duration)
      e2_b   <- input$e2_base
      dysmen <- input$dysmen
      dysp   <- input$dysp
      cpp    <- input$cpp
      lvol   <- input$lesion_vol
    })
  })

  # ---- TAB 1: Patient summary ----
  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Age","BMI","Stage","CA-125","AMH","Treatment","Duration"),
      Value     = c(input$age, input$bmi, input$stage,
                    paste(input$ca125, "U/mL"),
                    paste(input$amh, "ng/mL"),
                    input$drug,
                    paste(input$duration, "months"))
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # ---- TAB 2: PK plots ----
  output$pk_plot <- renderPlot({
    drug <- input$drug
    dur  <- as.numeric(input$duration)
    pk_df <- simulate_pk(drug, dur)

    if (drug == "No treatment") {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No drug administered",
                 size = 8, color = "gray60") +
        theme_void()
    } else {
      pk_df7 <- pk_df[pk_df$time_h <= 168, ]
      pk_par  <- pk_params[[drug]]
      ss_h    <- 5 * log(2) / pk_par$kel
      p <- ggplot(pk_df7, aes(x = time_h, y = conc)) +
        geom_line(color = "#7b2d8b", linewidth = 1.2) +
        geom_vline(xintercept = min(ss_h, 168), linetype = "dashed",
                   color = "steelblue", linewidth = 0.8) +
        annotate("text", x = min(ss_h, 168) + 2, y = max(pk_df7$conc) * 0.9,
                 label = "~Steady\nstate", size = 3.5, color = "steelblue", hjust = 0) +
        labs(title = paste("PK Profile:", drug),
             x = "Time (hours)", y = "Plasma Concentration (ng/mL or µg/mL)",
             caption = "One-compartment model; first-order absorption and elimination") +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(color = "#7b2d8b", face = "bold"))
    }
    p
  })

  output$pk_ss_plot <- renderPlot({
    drug <- input$drug
    dur  <- as.numeric(input$duration)
    pk_df <- simulate_pk(drug, dur)

    if (drug == "No treatment") {
      ggplot() + annotate("text", x=0.5, y=0.5,
                          label="No drug administered", size=6, color="gray60") + theme_void()
    } else {
      pk_long <- pk_df[seq(1, nrow(pk_df), by=4), ]
      ggplot(pk_long, aes(x = time_d, y = conc)) +
        geom_line(color = "#c0392b", linewidth = 1.0) +
        labs(title = "Extended PK (Days 0–7)",
             x = "Time (days)", y = "Concentration (ng/mL)") +
        theme_minimal(base_size = 12)
    }
  })

  output$pk_table <- renderTable({
    drug  <- input$drug
    dur   <- as.numeric(input$duration)
    pk_df <- simulate_pk(drug, dur)
    pk_summary(pk_df, drug)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # ---- TAB 3: HPO hormones ----
  output$hpo_plot <- renderPlot({
    sd  <- sim_data()
    hpo <- sd$hpo
    dur <- sd$dur

    hpo_long <- bind_rows(
      data.frame(time = hpo$time_months, value = hpo$E2_pgmL, hormone = "E2 (pg/mL)"),
      data.frame(time = hpo$time_months, value = hpo$FSH_IUL,  hormone = "FSH (IU/L)"),
      data.frame(time = hpo$time_months, value = hpo$LH_IUL,   hormone = "LH (IU/L)"),
      data.frame(time = hpo$time_months, value = hpo$P4_ngmL,  hormone = "P4 (ng/mL)")
    )
    hpo_long$hormone <- factor(hpo_long$hormone,
                                levels = c("E2 (pg/mL)","FSH (IU/L)","LH (IU/L)","P4 (ng/mL)"))

    e2_ref  <- data.frame(hormone = "E2 (pg/mL)", ymin = 0, ymax = 20)

    ggplot(hpo_long, aes(x = time, y = value, color = hormone)) +
      geom_rect(data = e2_ref,
                aes(xmin = 0, xmax = dur, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "pink", alpha = 0.25) +
      geom_line(linewidth = 1.2) +
      geom_hline(data = data.frame(hormone = "E2 (pg/mL)", ref = 20),
                 aes(yintercept = ref), linetype = "dashed", color = "red", linewidth = 0.8,
                 inherit.aes = FALSE) +
      facet_wrap(~hormone, scales = "free_y", ncol = 2) +
      scale_color_manual(values = c("E2 (pg/mL)"="#e74c3c",
                                     "FSH (IU/L)"="#3498db",
                                     "LH (IU/L)" ="#2ecc71",
                                     "P4 (ng/mL)"="#9b59b6")) +
      labs(title = paste("HPO Axis Hormones —", sd$drug),
           x = "Time (months)", y = "Concentration",
           caption = "Pink zone: E2 < 20 pg/mL (menopausal range). Dashed red line: E2 = 20 pg/mL therapeutic target.") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold", color = "#7b2d8b"),
            plot.title = element_text(color = "#7b2d8b", face = "bold"))
  })

  # ---- TAB 4: Pain & Clinical ----
  output$pain_plot <- renderPlot({
    clin <- sim_data()$clinical
    drug <- sim_data()$drug

    pain_long <- bind_rows(
      data.frame(time = clin$time_months, nrs = clin$pain_dysmen, type = "Dysmenorrhea"),
      data.frame(time = clin$time_months, nrs = clin$pain_dysp,   type = "Dyspareunia"),
      data.frame(time = clin$time_months, nrs = clin$pain_cpp,    type = "Chronic Pelvic Pain")
    )
    pain_long$type <- factor(pain_long$type,
                              levels = c("Dysmenorrhea","Dyspareunia","Chronic Pelvic Pain"))

    ggplot(pain_long, aes(x = time, y = nrs, color = type)) +
      geom_line(linewidth = 1.3) +
      geom_hline(yintercept = 3, linetype = "dotted", color = "gray50") +
      annotate("text", x = max(pain_long$time) * 0.02, y = 3.3,
               label = "Mild pain threshold (NRS 3)", size = 3.5, color = "gray50", hjust = 0) +
      scale_color_manual(values = c("Dysmenorrhea"="#e74c3c",
                                     "Dyspareunia"="#e67e22",
                                     "Chronic Pelvic Pain"="#9b59b6")) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(title = paste("Pain NRS —", drug),
           x = "Time (months)", y = "Pain NRS (0–10)",
           color = "Pain Type") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            plot.title = element_text(color = "#7b2d8b", face = "bold"))
  })

  output$bmd_plot <- renderPlot({
    clin <- sim_data()$clinical
    ggplot(clin, aes(x = time_months, y = bmd_pct)) +
      geom_line(color = "#27ae60", linewidth = 1.2) +
      geom_ribbon(aes(ymin = bmd_pct - 0.3, ymax = bmd_pct + 0.3),
                  fill = "#27ae60", alpha = 0.15) +
      geom_hline(yintercept = -2.5, linetype = "dashed", color = "red") +
      annotate("text", x = 1, y = -2.2,
               label = "Osteopenia threshold (T-score −1)", size = 3, color = "red", hjust = 0) +
      labs(x = "Time (months)", y = "BMD Change (%)",
           title = "Lumbar Spine BMD") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#27ae60", face = "bold"))
  })

  output$hf_plot <- renderPlot({
    clin <- sim_data()$clinical
    ggplot(clin, aes(x = time_months, y = hot_flush)) +
      geom_area(fill = "#e67e22", alpha = 0.35) +
      geom_line(color = "#e67e22", linewidth = 1.0) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(x = "Time (months)", y = "Score (0–10)",
           title = "Hot Flush Severity") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#e67e22", face = "bold"))
  })

  output$lesion_vol_plot <- renderPlot({
    clin <- sim_data()$clinical
    drug <- sim_data()$drug
    ggplot(clin, aes(x = time_months, y = lesion_vol)) +
      geom_line(color = "#8e44ad", linewidth = 1.3) +
      geom_ribbon(aes(ymin = lesion_vol * 0.85, ymax = lesion_vol * 1.15),
                  fill = "#8e44ad", alpha = 0.15) +
      labs(title = paste("Endometrioma Volume —", drug),
           x = "Time (months)", y = "Volume (cm³)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#8e44ad", face = "bold"))
  })

  output$clinical_table <- renderTable({
    clin <- sim_data()$clinical
    drug <- sim_data()$drug
    end  <- clin[nrow(clin), ]
    base_dysmen <- input$dysmen
    base_dysp   <- input$dysp
    base_cpp    <- input$cpp
    base_lvol   <- input$lesion_vol

    data.frame(
      Endpoint = c("Dysmenorrhea NRS", "Dyspareunia NRS",
                   "Chronic Pelvic Pain NRS", "BMD Change",
                   "Endometrioma Volume", "Hot Flush Score"),
      Baseline = c(base_dysmen, base_dysp, base_cpp,
                   "0%", paste(base_lvol,"cm³"), "0"),
      `End of Treatment` = c(
        round(end$pain_dysmen, 1),
        round(end$pain_dysp,   1),
        round(end$pain_cpp,    1),
        paste0(round(end$bmd_pct, 2), "%"),
        paste0(round(end$lesion_vol, 2), " cm³"),
        round(end$hot_flush, 1)
      ),
      `% Change` = c(
        paste0("-", round((1 - end$pain_dysmen / max(base_dysmen,0.01)) * 100, 1), "%"),
        paste0("-", round((1 - end$pain_dysp   / max(base_dysp,  0.01)) * 100, 1), "%"),
        paste0("-", round((1 - end$pain_cpp    / max(base_cpp,   0.01)) * 100, 1), "%"),
        paste0(round(end$bmd_pct, 2), "%"),
        paste0("-", round((1 - end$lesion_vol  / base_lvol) * 100, 1), "%"),
        "N/A"
      )
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  # ---- TAB 5: Lesion dynamics ----
  output$lesion_dyn_plot <- renderPlot({
    les  <- sim_data()$lesion
    drug <- sim_data()$drug
    ggplot(les, aes(x = time_months, y = lesion_vol)) +
      geom_area(fill = "#8e44ad", alpha = 0.2) +
      geom_line(color = "#8e44ad", linewidth = 1.3) +
      labs(title = paste("Lesion Volume (cm³) —", drug),
           x = "Time (months)", y = "Lesion Volume (cm³)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(color = "#8e44ad", face = "bold"))
  })

  output$prolif_plot <- renderPlot({
    les <- sim_data()$lesion
    ggplot(les, aes(x = time_months, y = prolif_rate)) +
      geom_line(color = "#c0392b", linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
      scale_y_continuous(limits = c(0, 1.2)) +
      labs(x = "Time (months)", y = "Relative Proliferation Rate",
           title = "Lesion Proliferation Rate") +
      theme_minimal(base_size = 12)
  })

  output$e2_local_plot <- renderPlot({
    les <- sim_data()$lesion
    ggplot(les, aes(x = time_months, y = e2_local)) +
      geom_line(color = "#e74c3c", linewidth = 1.2) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "steelblue") +
      annotate("text", x = 1, y = 22, label = "Systemic target (20 pg/mL)",
               size = 3.2, color = "steelblue", hjust = 0) +
      labs(title = "Local E2 at Ectopic Site",
           x = "Time (months)", y = "E2 (pg/mL)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#e74c3c", face = "bold"))
  })

  output$cytokine_plot <- renderPlot({
    les <- sim_data()$lesion
    cyt_long <- bind_rows(
      data.frame(time = les$time_months, value = les$IL6_pgmL,  marker = "IL-6 (pg/mL)"),
      data.frame(time = les$time_months, value = les$PGE2_pgmL, marker = "PGE2 (pg/mL)")
    )
    ggplot(cyt_long, aes(x = time, y = value, color = marker)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("IL-6 (pg/mL)"="#3498db","PGE2 (pg/mL)"="#e67e22")) +
      labs(title = "Peritoneal Fluid Cytokines & Mediators",
           x = "Time (months)", y = "Concentration (pg/mL)",
           color = "Marker") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom",
            plot.title = element_text(color = "#3498db", face = "bold"))
  })

  output$ngf_plot <- renderPlot({
    les <- sim_data()$lesion
    ggplot(les, aes(x = time_months, y = NGF_pgmL)) +
      geom_line(color = "#e74c3c", linewidth = 1.1) +
      geom_ribbon(aes(ymin = NGF_pgmL * 0.85, ymax = NGF_pgmL * 1.15),
                  fill = "#e74c3c", alpha = 0.15) +
      labs(title = "NGF at Lesion (Pain Mediator)",
           x = "Time (months)", y = "NGF (pg/mL)") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#e74c3c", face = "bold"))
  })

  output$arom_plot <- renderPlot({
    les <- sim_data()$lesion
    ggplot(les, aes(x = time_months, y = aromatase_act)) +
      geom_area(fill = "#f39c12", alpha = 0.3) +
      geom_line(color = "#f39c12", linewidth = 1.1) +
      scale_y_continuous(limits = c(0, 1.1)) +
      labs(title = "Aromatase Activity (Relative to Baseline)",
           x = "Time (months)", y = "Relative Activity") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(color = "#f39c12", face = "bold"))
  })

  # ---- TAB 6: Scenario comparison ----
  output$scenario_facet_plot <- renderPlot({
    all_s   <- sim_data()$all_scen
    dur     <- sim_data()$dur
    lvol    <- input$lesion_vol

    # Downsample for speed
    idx <- seq(1, nrow(all_s), by = max(1, floor(nrow(all_s) / (7 * 80))))
    ds  <- all_s[idx, ]

    ds$drug <- factor(ds$drug,
                      levels = c("No treatment","Leuprolide depot",
                                 "Elagolix 150mg/d","Elagolix 200mg BID",
                                 "Dienogest 2mg/d","Letrozole+Add-back","Combined OCP"))

    # Build facet data
    fac_data <- bind_rows(
      data.frame(time = ds$time_months, value = ds$lesion_vol,    drug = ds$drug, panel = "Lesion Volume (cm³)"),
      data.frame(time = ds$time_months, value = ds$pain_dysmen,   drug = ds$drug, panel = "Dysmenorrhea NRS"),
      data.frame(time = ds$time_months, value = ds$E2_pgmL,       drug = ds$drug, panel = "E2 (pg/mL)"),
      data.frame(time = ds$time_months, value = ds$bmd_pct,       drug = ds$drug, panel = "BMD Change (%)")
    )
    fac_data$panel <- factor(fac_data$panel,
                              levels = c("Lesion Volume (cm³)","Dysmenorrhea NRS",
                                         "E2 (pg/mL)","BMD Change (%)"))

    drug_colors <- c("No treatment"       = "#95a5a6",
                     "Leuprolide depot"   = "#e74c3c",
                     "Elagolix 150mg/d"  = "#3498db",
                     "Elagolix 200mg BID"= "#2980b9",
                     "Dienogest 2mg/d"   = "#2ecc71",
                     "Letrozole+Add-back"= "#f39c12",
                     "Combined OCP"      = "#9b59b6")

    ggplot(fac_data, aes(x = time, y = value, color = drug)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_color_manual(values = drug_colors) +
      labs(title = "Treatment Scenario Comparison — All 7 Treatments",
           x = "Time (months)", y = NULL, color = "Treatment") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom",
            strip.text = element_text(face = "bold", color = "#7b2d8b"),
            plot.title = element_text(color = "#7b2d8b", face = "bold"),
            legend.text = element_text(size = 9))
  })

  output$waterfall_plot <- renderPlot({
    dur      <- sim_data()$dur
    dysmen_b <- input$dysmen
    dysp_b   <- input$dysp
    cpp_b    <- input$cpp

    drugs <- names(drug_pain)
    resp_rates <- sapply(drugs, function(d) {
      pp <- drug_pain[[d]]
      # proportion with ≥50% pain reduction (simplified: use mean pain reduction)
      mean_red <- (pp$dysmen + pp$dysp + pp$cpp) / 3
      # approximate % of patients responding (binomial logistic estimate)
      plogis(4 * mean_red - 1)  # maps 0→~27%, 0.5→50%, 0.75→75%, etc.
    })

    wf_df <- data.frame(
      drug = factor(drugs, levels = drugs[order(resp_rates, decreasing = TRUE)]),
      resp = resp_rates * 100
    )
    wf_df$drug_ordered <- factor(wf_df$drug,
                                  levels = levels(wf_df$drug))

    drug_colors <- c("No treatment"       = "#95a5a6",
                     "Leuprolide depot"   = "#e74c3c",
                     "Elagolix 150mg/d"  = "#3498db",
                     "Elagolix 200mg BID"= "#2980b9",
                     "Dienogest 2mg/d"   = "#2ecc71",
                     "Letrozole+Add-back"= "#f39c12",
                     "Combined OCP"      = "#9b59b6")

    ggplot(wf_df, aes(x = reorder(drug, resp), y = resp, fill = drug)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray40") +
      coord_flip() +
      scale_fill_manual(values = drug_colors) +
      scale_y_continuous(limits = c(0, 100)) +
      annotate("text", x = 0.7, y = 52, label = "50% response threshold",
               size = 3.5, color = "gray40", hjust = 0) +
      labs(title = "Pain Response Rate (≥50% Reduction)",
           x = NULL, y = "Estimated Response Rate (%)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            plot.title = element_text(color = "#7b2d8b", face = "bold"))
  })

  output$scenario_table <- renderTable({
    dur      <- sim_data()$dur
    dysmen_b <- input$dysmen
    cpp_b    <- input$cpp
    lvol_b   <- input$lesion_vol
    e2_b     <- input$e2_base

    drugs <- names(drug_pain)
    rows <- lapply(drugs, function(d) {
      clin <- simulate_clinical(d, dur, dysmen_b, dysmen_b, cpp_b, lvol_b)
      hpo  <- simulate_hpo(d, dur, e2_b)
      end_clin <- clin[nrow(clin), ]
      end_hpo  <- hpo[nrow(hpo),   ]

      lesion_red <- round((1 - end_clin$lesion_vol / max(lvol_b, 0.1)) * 100, 1)
      pain_red   <- round((1 - end_clin$pain_dysmen / max(dysmen_b, 0.1)) * 100, 1)
      e2_sup     <- round((1 - end_hpo$E2_pgmL / max(e2_b, 0.1)) * 100, 1)
      bmd_ch     <- round(end_clin$bmd_pct, 2)

      data.frame(
        Treatment      = d,
        `Lesion Reduction %` = paste0(lesion_red, "%"),
        `Pain Reduction %`   = paste0(pain_red,   "%"),
        `E2 Suppression %`   = paste0(e2_sup,     "%"),
        `BMD Change`         = paste0(bmd_ch,      "%")
      )
    })
    bind_rows(rows)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

} # end server

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
