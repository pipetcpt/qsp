## ============================================================================
## Bipolar Disorder QSP – Interactive Shiny Dashboard
## ============================================================================
## Tabs:
##   1. Patient Profile     – demographics, diagnosis, episode history
##   2. Pharmacokinetics    – PK curves, Css, therapeutic window
##   3. PD Biomarkers       – GSK-3β, BDNF, IL-6, cortisol dynamics
##   4. Clinical Endpoints  – YMRS / MADRS trajectories, response rates
##   5. Scenario Comparison – side-by-side 6-scenario analysis
##   6. Safety Monitor      – toxicity flags, QTc, metabolic risk
## ============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## ---- Embedded ODE solver (pure R, no mrgsolve dependency in app) ----------
## Simple Euler integrator for demonstration
run_bd_model <- function(
    dose_Li = 8.1, freq_Li = 8,    # mmol, h
    dose_VPA = 0,  freq_VPA = 12,
    dose_QTP = 0,  freq_QTP = 24,
    dose_LTG = 0,  freq_LTG = 24,
    YMRS_init = 25, MADRS_init = 5,
    duration_d = 56,
    dt = 0.5
) {
  # Parameters (same as mrgsolve model)
  p <- list(
    ka_Li=0.80, Vc_Li=30, CL_Li=1.80,
    ka_VPA=1.50, Vc_VPA=14, CL_VPA=0.55, fu_VPA0=0.10, Km_fu=50,
    ka_QTP=1.10, Vc_QTP=900, CL_QTP=250, F_QTP=0.09, km_QTP=0.25, CL_NQT=80,
    ka_LTG=0.45, Vc_LTG=105, CL_LTG=1.85,
    kDA_syn=0.20, kDA_deg=0.20,
    k5HT_syn=0.15, k5HT_deg=0.15,
    kGSK_syn=0.05, kGSK_deg=0.05,
    IC50_Li_GSK=0.70, IC50_VPA_GSK=60, Emax_GSK=0.90,
    kBDNF_syn=0.03, kBDNF_deg=0.03, Emax_BDNF_Li=0.50, EC50_BDNF_Li=0.40,
    kIL6_syn=0.04, kIL6_deg=0.04, Emax_IL6_Li=0.40, EC50_IL6_Li=0.50,
    kCort_prod=0.08, kCort_deg=0.08,
    omega=2*pi/24, Amp_circ=0.30,
    Emax_YMRS_Li=18, EC50_YMRS_Li=0.65,
    Emax_YMRS_VPA=20, EC50_YMRS_VPA=65,
    Emax_YMRS_QTP=16, EC50_YMRS_QTP=120,
    Emax_MADRS_QTP=15, EC50_MADRS_QTP=80,
    Emax_MADRS_LTG=10, EC50_MADRS_LTG=2.5,
    ALPHA_combo=0.4, kYMRS_nat=0.01, kMADRS_nat=0.005,
    kWt_gain=0.003,
    Imax_D2=0.8
  )

  # Time vector
  times <- seq(0, duration_d * 24, by = dt)
  n <- length(times)

  # State variables
  s <- list(
    Li_gut=0, Li_central=0,
    VPA_gut=0, VPA_central=0,
    QTP_gut=0, QTP_central=0, NQT_central=0,
    LTG_gut=0, LTG_central=0,
    DA=1, HT5=1, GSK3=1, BDNF=1, IL6=1, Cort=1,
    YMRS=YMRS_init, MADRS=MADRS_init,
    Wt=0, GAF=50
  )

  # Storage
  out <- data.frame(matrix(NA, nrow=n, ncol=20))
  names(out) <- c("time","Li_mEqL","VPA_ugmL","VPA_free","QTP_ngmL","NQT_ngmL","LTG_ugmL",
                  "DA","HT5","GSK3","BDNF","IL6","Cortisol",
                  "YMRS","MADRS","GAF","Wt",
                  "YMRS_resp","MADRS_resp","MADRS_remit")

  for (i in seq_along(times)) {
    t <- times[i]

    # ---- Dosing ----
    if (freq_Li  > 0 && dose_Li  > 0 && i > 1) {
      if ((t %% freq_Li)  < dt) s$Li_gut  <- s$Li_gut  + dose_Li
    }
    if (freq_VPA > 0 && dose_VPA > 0 && i > 1) {
      if ((t %% freq_VPA) < dt) s$VPA_gut <- s$VPA_gut + dose_VPA
    }
    if (freq_QTP > 0 && dose_QTP > 0 && i > 1) {
      if ((t %% freq_QTP) < dt) s$QTP_gut <- s$QTP_gut + dose_QTP
    }
    if (freq_LTG > 0 && dose_LTG > 0 && i > 1) {
      if ((t %% freq_LTG) < dt) s$LTG_gut <- s$LTG_gut + dose_LTG
    }

    # ---- Derived concentrations ----
    Li_c   <- s$Li_central / p$Vc_Li
    VPA_c  <- s$VPA_central / p$Vc_VPA
    VPA_f  <- VPA_c * p$fu_VPA0 * (1 + VPA_c / p$Km_fu)
    QTP_c  <- s$QTP_central / (p$Vc_QTP / 1000)
    NQT_c  <- s$NQT_central / (p$Vc_QTP / 1000)
    LTG_c  <- s$LTG_central / p$Vc_LTG

    # ---- Drug effects ----
    iGSK_Li  <- p$Emax_GSK * Li_c  / (p$IC50_Li_GSK  + Li_c)
    iGSK_VPA <- p$Emax_GSK * VPA_f / (p$IC50_VPA_GSK + VPA_f)
    GSK_inh  <- min(0.95, 1 - (1 - iGSK_Li) * (1 - iGSK_VPA))
    BDNF_stim <- p$Emax_BDNF_Li * Li_c / (p$EC50_BDNF_Li + Li_c)
    IL6_inh   <- p$Emax_IL6_Li  * Li_c / (p$EC50_IL6_Li  + Li_c)
    D2_occ   <- QTP_c / (p$EC50_YMRS_QTP + QTP_c)
    DA_inh   <- p$Imax_D2 * D2_occ

    E_YMRS_Li  <- p$Emax_YMRS_Li  * Li_c  / (p$EC50_YMRS_Li  + Li_c)
    E_YMRS_VPA <- p$Emax_YMRS_VPA * VPA_f / (p$EC50_YMRS_VPA + VPA_f)
    E_YMRS_QTP <- p$Emax_YMRS_QTP * QTP_c / (p$EC50_YMRS_QTP + QTP_c)
    E_YMRS     <- E_YMRS_Li + E_YMRS_VPA + E_YMRS_QTP

    E_MADRS_QTP   <- p$Emax_MADRS_QTP * QTP_c / (p$EC50_MADRS_QTP + QTP_c)
    E_MADRS_LTG   <- p$Emax_MADRS_LTG * LTG_c / (p$EC50_MADRS_LTG + LTG_c)
    combo_bonus   <- p$ALPHA_combo * (if (Li_c > 0.3 && QTP_c > 50) 4.0 else 0)
    E_MADRS       <- E_MADRS_QTP + E_MADRS_LTG + combo_bonus

    Circ <- p$Amp_circ * sin(p$omega * t)

    # ---- ODEs (Euler) ----
    d_Li_gut     <- -p$ka_Li  * s$Li_gut
    d_Li_central <-  p$ka_Li  * s$Li_gut - (p$CL_Li/p$Vc_Li) * s$Li_central
    d_VPA_gut    <- -p$ka_VPA * s$VPA_gut
    d_VPA_c      <-  p$ka_VPA * s$VPA_gut - (p$CL_VPA/p$Vc_VPA) * s$VPA_central
    d_QTP_gut    <- -p$ka_QTP * s$QTP_gut
    d_QTP_c      <-  p$ka_QTP * p$F_QTP * s$QTP_gut - (p$CL_QTP/(p$Vc_QTP/1000)) * s$QTP_central
    d_NQT_c      <-  p$km_QTP * (p$CL_QTP/(p$Vc_QTP/1000)) * s$QTP_central - (p$CL_NQT/(p$Vc_QTP/1000)) * s$NQT_central
    d_LTG_gut    <- -p$ka_LTG * s$LTG_gut
    d_LTG_c      <-  p$ka_LTG * s$LTG_gut - (p$CL_LTG/p$Vc_LTG) * s$LTG_central

    d_DA    <- p$kDA_syn - p$kDA_deg * s$DA - DA_inh * s$DA
    d_HT5   <- p$k5HT_syn * (1 + 0.3 * NQT_c/(80 + NQT_c)) - p$k5HT_deg * s$HT5
    d_GSK3  <- p$kGSK_syn * (1 - GSK_inh) - p$kGSK_deg * s$GSK3
    d_BDNF  <- p$kBDNF_syn * (1 + BDNF_stim) - p$kBDNF_deg * s$BDNF
    d_IL6   <- p$kIL6_syn * (1 - IL6_inh) - p$kIL6_deg * s$IL6
    d_Cort  <- p$kCort_prod * (1 + Circ) - p$kCort_deg * s$Cort

    YMRS_d <- max(0, s$YMRS)
    MADRS_d <- max(0, s$MADRS)
    d_YMRS  <- YMRS_d * (s$DA - 1) * 0.05 - p$kYMRS_nat * YMRS_d -
               E_YMRS * (YMRS_d / (YMRS_init + 1e-6)) * p$kYMRS_nat * 50
    d_MADRS <- MADRS_d * ((s$IL6 - 1) * 0.05 + (s$Cort - 1) * 0.03 + (1 - s$BDNF) * 0.04) -
               p$kMADRS_nat * MADRS_d - E_MADRS * (MADRS_d / (MADRS_init + 1e-6)) * p$kMADRS_nat * 50
    d_Wt   <- p$kWt_gain * NQT_c
    d_GAF  <- 0.01 * (70 - s$GAF) - 0.05 * (YMRS_d + MADRS_d) / 30 * s$GAF * 0.01

    # Update state (Euler)
    s$Li_gut     <- max(0, s$Li_gut     + d_Li_gut     * dt)
    s$Li_central <- max(0, s$Li_central + d_Li_central * dt)
    s$VPA_gut    <- max(0, s$VPA_gut    + d_VPA_gut    * dt)
    s$VPA_central<- max(0, s$VPA_central+ d_VPA_c      * dt)
    s$QTP_gut    <- max(0, s$QTP_gut    + d_QTP_gut    * dt)
    s$QTP_central<- max(0, s$QTP_central+ d_QTP_c      * dt)
    s$NQT_central<- max(0, s$NQT_central+ d_NQT_c      * dt)
    s$LTG_gut    <- max(0, s$LTG_gut    + d_LTG_gut    * dt)
    s$LTG_central<- max(0, s$LTG_central+ d_LTG_c      * dt)
    s$DA   <- max(0.01, s$DA   + d_DA   * dt)
    s$HT5  <- max(0.01, s$HT5  + d_HT5  * dt)
    s$GSK3 <- max(0.01, s$GSK3 + d_GSK3 * dt)
    s$BDNF <- max(0.01, s$BDNF + d_BDNF * dt)
    s$IL6  <- max(0.01, s$IL6  + d_IL6  * dt)
    s$Cort <- max(0.01, s$Cort + d_Cort * dt)
    s$YMRS <- max(0, s$YMRS + d_YMRS * dt)
    s$MADRS<- max(0, s$MADRS + d_MADRS * dt)
    s$Wt   <- s$Wt + d_Wt * dt
    s$GAF  <- min(100, max(0, s$GAF + d_GAF * dt))

    out[i, ] <- list(t, Li_c, VPA_c, VPA_f, QTP_c, NQT_c, LTG_c,
                     s$DA, s$HT5, s$GSK3, s$BDNF, s$IL6, s$Cort,
                     s$YMRS, s$MADRS, s$GAF, s$Wt,
                     as.integer(s$YMRS <= YMRS_init * 0.5),
                     as.integer(s$MADRS <= MADRS_init * 0.5),
                     as.integer(s$MADRS <= 12))
  }
  out$day <- out$time / 24
  out
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Bipolar Disorder QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",   icon = icon("user")),
      menuItem("Pharmacokinetics",    tabName = "tab_pk",        icon = icon("flask")),
      menuItem("PD Biomarkers",       tabName = "tab_pd",        icon = icon("dna")),
      menuItem("Clinical Endpoints",  tabName = "tab_endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "tab_scenario",  icon = icon("balance-scale")),
      menuItem("Safety Monitor",      tabName = "tab_safety",    icon = icon("shield-alt"))
    ),

    ## Global controls
    hr(),
    h5("  Treatment Setup", style = "color:#ccc; padding-left:10px"),
    sliderInput("dose_Li",  "Lithium dose (mmol/dose)", 0, 16.2, 8.1, 0.9,
                post = " mmol"),
    selectInput("freq_Li",  "Lithium frequency", c("BID(12h)"=12,"TID(8h)"=8,"QD(24h)"=24),
                selected = 8),
    sliderInput("dose_QTP", "Quetiapine (mg/dose)", 0, 800, 300, 25, post = " mg"),
    sliderInput("dose_VPA", "Valproate (mg/dose)", 0, 1000, 0, 50, post = " mg"),
    sliderInput("dose_LTG", "Lamotrigine (mg/dose)", 0, 200, 0, 25, post = " mg"),
    hr(),
    h5("  Episode Type", style = "color:#ccc; padding-left:10px"),
    selectInput("episode", "Episode",
                c("Acute Mania"="mania",
                  "Bipolar Depression"="bdep",
                  "Mixed/Dysphoric"="mixed")),
    sliderInput("duration_d", "Simulation (days)", 14, 365, 56, 7),
    actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tabItems(

      ## ------------------------------------------------------------------
      ## TAB 1 – Patient Profile
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics", width = 4, status = "primary",
            numericInput("age",    "Age (years)",          35, 18, 80, 1),
            selectInput("sex",     "Sex", c("Male","Female")),
            numericInput("weight", "Weight (kg)",          70, 40, 150, 1),
            numericInput("creatinine", "Serum Creatinine (mg/dL)", 0.9, 0.5, 5.0, 0.1),
            hr(),
            h5("Calculated eGFR:"),
            verbatimTextOutput("eGFR_out")
          ),
          box(title = "Diagnosis & History", width = 4, status = "warning",
            selectInput("bd_type", "BD Subtype",
                        c("BD-I (manic)","BD-II (hypomanic/depressive)","BD-NOS","Schizoaffective BD")),
            selectInput("episode_onset", "Episode Onset",
                        c("First episode","2nd–3rd episode","≥4 recurrences","Rapid cycling")),
            numericInput("episodes_year", "Episodes per year (lifetime avg)", 1.5, 0, 12, 0.5),
            selectInput("psychosis_history", "Psychosis history", c("No","Yes – resolved","Yes – current")),
            checkboxGroupInput("comorbid", "Comorbidities",
                               c("Anxiety disorder","ADHD","Substance use",
                                 "Hypothyroidism","Metabolic syndrome","Migraine"))
          ),
          box(title = "Current Episode Severity", width = 4, status = "danger",
            sliderInput("ymrs_init",  "YMRS (baseline)",  0, 60, 25, 1),
            sliderInput("madrs_init", "MADRS (baseline)", 0, 60, 5,  1),
            hr(),
            h5("Episode Severity Classification:"),
            verbatimTextOutput("severity_out"),
            hr(),
            h5("Recommended Tier (CANMAT 2018):"),
            verbatimTextOutput("canmat_rec")
          )
        ),
        fluidRow(
          box(title = "Pharmacogenomics", width = 6, status = "info",
            selectInput("cacna1c", "CACNA1C (rs1006737) genotype",
                        c("GG (reference)","GA (1 risk allele)","AA (2 risk alleles)")),
            selectInput("bdnf_val66met", "BDNF Val66Met (rs6265)",
                        c("Val/Val (WT)","Val/Met","Met/Met")),
            selectInput("comt", "COMT Val158Met (rs4680)",
                        c("Val/Val (high COMT)","Val/Met","Met/Met (low COMT)")),
            selectInput("cyp3a4", "CYP3A4 phenotype",
                        c("Normal metaboliser","Poor metaboliser","Ultra-rapid metaboliser"))
          ),
          box(title = "PGx Implications Summary", width = 6, status = "info",
            verbatimTextOutput("pgx_summary")
          )
        )
      ),

      ## ------------------------------------------------------------------
      ## TAB 2 – Pharmacokinetics
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Concentration – Time Profiles", width = 12, status = "primary",
            plotOutput("pk_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Steady-State PK Summary", width = 6, status = "info",
            tableOutput("pk_summary_tbl")
          ),
          box(title = "Therapeutic Window", width = 6, status = "warning",
            plotOutput("pk_window_plot", height = "260px")
          )
        )
      ),

      ## ------------------------------------------------------------------
      ## TAB 3 – PD Biomarkers
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "Signal Transduction Biomarkers", width = 12, status = "success",
            plotOutput("pd_biomarker_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Neuroplasticity Indices", width = 6, status = "success",
            plotOutput("pd_neuroplast_plot", height = "260px")
          ),
          box(title = "Neuroinflammation & HPA", width = 6, status = "danger",
            plotOutput("pd_inflam_plot", height = "260px")
          )
        )
      ),

      ## ------------------------------------------------------------------
      ## TAB 4 – Clinical Endpoints
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title = "YMRS (Mania Score) – Time Course", width = 6, status = "danger",
            plotOutput("ymrs_plot", height = "300px")
          ),
          box(title = "MADRS (Depression Score) – Time Course", width = 6, status = "info",
            plotOutput("madrs_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "GAF / Functioning Score", width = 6, status = "success",
            plotOutput("gaf_plot", height = "260px")
          ),
          box(title = "Clinical Outcome KPIs", width = 6, status = "primary",
            valueBoxOutput("vbox_ymrs_resp",   width = 12),
            valueBoxOutput("vbox_madrs_resp",  width = 12),
            valueBoxOutput("vbox_remission",   width = 12),
            valueBoxOutput("vbox_gaf_end",     width = 12)
          )
        )
      ),

      ## ------------------------------------------------------------------
      ## TAB 5 – Scenario Comparison
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Predefined Scenarios (CANMAT-aligned)", width = 12, status = "primary",
            fluidRow(
              column(4,
                h5("Acute Mania Scenarios"),
                checkboxGroupInput("mania_scen", NULL,
                  choices = c("Lithium 900 mg/d"="sc1",
                              "Valproate 1000 mg/d"="sc2",
                              "Quetiapine 600 mg/d"="sc3"),
                  selected = c("sc1","sc2"))
              ),
              column(4,
                h5("BD Depression Scenarios"),
                checkboxGroupInput("dep_scen", NULL,
                  choices = c("Quetiapine 300 mg/d"="sc4",
                              "Li + Quetiapine 300"="sc5",
                              "Lamotrigine titration"="sc6"),
                  selected = c("sc4","sc5"))
              ),
              column(4,
                h5("Maintenance Scenarios"),
                checkboxGroupInput("maint_scen", NULL,
                  choices = c("Lithium maintenance"="sc7",
                              "VPA maintenance"="sc8"),
                  selected = "sc7")
              )
            ),
            actionButton("run_compare", "Compare Scenarios", class = "btn-success btn-block",
                         icon = icon("chart-bar"))
          )
        ),
        fluidRow(
          box(title = "YMRS Comparison", width = 6, status = "danger",
            plotOutput("compare_ymrs", height = "300px")
          ),
          box(title = "MADRS Comparison", width = 6, status = "info",
            plotOutput("compare_madrs", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Scenario Summary Table", width = 12,
            DT::dataTableOutput("compare_tbl")
          )
        )
      ),

      ## ------------------------------------------------------------------
      ## TAB 6 – Safety Monitor
      ## ------------------------------------------------------------------
      tabItem(tabName = "tab_safety",
        fluidRow(
          valueBoxOutput("vbox_li_safe",  width = 4),
          valueBoxOutput("vbox_vpa_safe", width = 4),
          valueBoxOutput("vbox_qtc",      width = 4)
        ),
        fluidRow(
          box(title = "Lithium Safety: Concentration vs Therapeutic Window", width = 6, status = "warning",
            plotOutput("li_safety_plot", height = "280px")
          ),
          box(title = "Metabolic Risk: Weight & Insulin Resistance", width = 6, status = "danger",
            plotOutput("metabolic_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Safety Summary & Monitoring Checklist", width = 12, status = "warning",
            tableOutput("safety_tbl")
          )
        )
      )
    )
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## ---- Reactive simulation result ----
  sim_result <- eventReactive(input$run_sim, {
    # Map episode to initial scores
    YMRS_init  <- if (input$episode == "mania")  input$ymrs_init  else 5
    MADRS_init <- if (input$episode == "bdep")   input$madrs_init else 5
    if (input$episode == "mixed") {
      YMRS_init  <- max(input$ymrs_init, 12)
      MADRS_init <- max(input$madrs_init, 12)
    }
    run_bd_model(
      dose_Li  = input$dose_Li,  freq_Li  = as.numeric(input$freq_Li),
      dose_VPA = input$dose_VPA, freq_VPA = 12,
      dose_QTP = input$dose_QTP, freq_QTP = 24,
      dose_LTG = input$dose_LTG, freq_LTG = 24,
      YMRS_init  = YMRS_init,
      MADRS_init = MADRS_init,
      duration_d = input$duration_d,
      dt = 0.5
    )
  }, ignoreNULL = FALSE)

  ## ---- TAB 1: Patient Profile ----
  output$eGFR_out <- renderText({
    # CKD-EPI simplified
    age <- input$age; wt <- input$weight; cr <- input$creatinine
    sex_factor <- if (input$sex == "Female") 0.85 else 1.0
    egfr <- 186 * (cr^(-1.154)) * (age^(-0.203)) * sex_factor
    sprintf("eGFR = %.0f mL/min/1.73m²\n%s",
            egfr,
            if (egfr >= 60) "Normal / CKD G1-2 – standard Li dosing"
            else if (egfr >= 30) "CKD G3 – reduce Li dose, monitor closely"
            else "CKD G4-5 – Li generally contraindicated")
  })

  output$severity_out <- renderText({
    ym <- input$ymrs_init; md <- input$madrs_init
    mania_sev <- if (ym >= 35) "Severe Mania" else if (ym >= 20) "Moderate Mania"
                 else if (ym >= 12) "Mild Mania" else "Euthymic/Subthreshold"
    dep_sev   <- if (md >= 35) "Severe Depression" else if (md >= 22) "Moderate Depression"
                 else if (md >= 13) "Mild Depression" else "Euthymic"
    sprintf("%s (YMRS=%d)\n%s (MADRS=%d)", mania_sev, ym, dep_sev, md)
  })

  output$canmat_rec <- renderText({
    ep <- input$episode
    if (ep == "mania")
      "1st line: Li / VPA / QTP / risperidone\n2nd line: Li+VPA, Li+QTP, asenapine"
    else if (ep == "bdep")
      "1st line: QTP, Li, LTG (adjunct), lurasidone\n2nd line: Li+LTG, Li+QTP"
    else
      "Mixed: QTP, aripiprazole, asenapine\nAvoid AD monotherapy"
  })

  output$pgx_summary <- renderText({
    lines <- character(0)
    if (grepl("GA|AA", input$cacna1c))
      lines <- c(lines, "CACNA1C risk allele: ↑ L-type Ca²⁺, may benefit from VPA")
    if (grepl("Val/Met|Met/Met", input$bdnf_val66met))
      lines <- c(lines, "BDNF Val66Met: ↓ activity-dep secretion; Li BDNF benefit attenuated")
    if (grepl("Met/Met", input$comt))
      lines <- c(lines, "COMT Met/Met: lower COMT → ↑ PFC DA; avoid high-dose D2 blockade")
    if (grepl("Ultra", input$cyp3a4))
      lines <- c(lines, "CYP3A4 ultra-rapid: QTP doses may need ↑ 30–50%")
    if (grepl("Poor", input$cyp3a4))
      lines <- c(lines, "CYP3A4 poor: QTP exposures ↑; start at 25–50 mg and titrate slowly")
    if (length(lines) == 0) lines <- "No major PGx alerts for current genotype selections."
    paste(lines, collapse = "\n")
  })

  ## ---- TAB 2: Pharmacokinetics ----
  output$pk_plot <- renderPlot({
    df <- sim_result()
    df_long <- df %>%
      select(day, Li_mEqL, VPA_ugmL, QTP_ngmL, LTG_ugmL) %>%
      pivot_longer(-day, names_to = "Drug", values_to = "Conc") %>%
      mutate(Drug = recode(Drug,
        "Li_mEqL"  = "Lithium (mEq/L)",
        "VPA_ugmL" = "Valproate (μg/mL)",
        "QTP_ngmL" = "Quetiapine (ng/mL) /10",
        "LTG_ugmL" = "Lamotrigine (μg/mL)"),
        Conc = ifelse(grepl("Quetiapine", Drug), Conc / 10, Conc))
    ggplot(df_long, aes(day, Conc, color = Drug)) +
      geom_line(size = 0.9) +
      labs(title = "Plasma Concentration – Time Profiles",
           x = "Day", y = "Concentration (units per legend)") +
      theme_bw(base_size = 13) + scale_color_brewer(palette = "Set1")
  })

  output$pk_summary_tbl <- renderTable({
    df <- sim_result()
    last_24h <- df %>% filter(day >= max(day) - 1)
    tibble(
      Drug = c("Lithium","Valproate","Quetiapine","Lamotrigine"),
      "Css mean" = c(
        round(mean(last_24h$Li_mEqL), 2),
        round(mean(last_24h$VPA_ugmL), 1),
        round(mean(last_24h$QTP_ngmL), 0),
        round(mean(last_24h$LTG_ugmL), 2)
      ),
      Unit = c("mEq/L","μg/mL","ng/mL","μg/mL"),
      "Therapeutic range" = c("0.6–1.2","50–100","100–400","1–4")
    )
  })

  output$pk_window_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, Li_mEqL)) +
      geom_line(color = "steelblue", size = 1) +
      geom_ribbon(aes(ymin = 0.6, ymax = 1.2), alpha = 0.2, fill = "green3") +
      geom_hline(yintercept = 1.5, color = "red3", linetype = "dashed") +
      annotate("text", x = max(df$day) * 0.7, y = 1.55,
               label = "Toxic (>1.5)", color = "red3", size = 3.5) +
      labs(title = "Lithium Therapeutic Window",
           x = "Day", y = "Li (mEq/L)") +
      theme_bw(base_size = 13)
  })

  ## ---- TAB 3: PD Biomarkers ----
  output$pd_biomarker_plot <- renderPlot({
    df <- sim_result()
    df_bio <- df %>%
      select(day, GSK3, BDNF, IL6, Cortisol) %>%
      pivot_longer(-day, names_to = "Biomarker", values_to = "Value")
    ggplot(df_bio, aes(day, Value, color = Biomarker)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "gray50") +
      labs(title = "Signal Transduction & Biomarker Dynamics",
           x = "Day", y = "Normalised Index (1 = baseline)") +
      theme_bw(base_size = 13) +
      scale_color_manual(values = c(
        GSK3="firebrick", BDNF="green4", IL6="steelblue", Cortisol="purple"))
  })

  output$pd_neuroplast_plot <- renderPlot({
    df <- sim_result()
    df_np <- df %>% select(day, BDNF, DA, HT5) %>%
      pivot_longer(-day, names_to = "Marker", values_to = "Value")
    ggplot(df_np, aes(day, Value, color = Marker)) +
      geom_line(size = 1) +
      labs(title = "Neuroplasticity: BDNF / DA / 5-HT",
           x = "Day", y = "Normalised Index") +
      theme_bw(base_size = 13) +
      scale_color_brewer(palette = "Dark2")
  })

  output$pd_inflam_plot <- renderPlot({
    df <- sim_result()
    df_inf <- df %>% select(day, IL6, Cortisol, GSK3) %>%
      pivot_longer(-day, names_to = "Marker", values_to = "Value")
    ggplot(df_inf, aes(day, Value, color = Marker)) +
      geom_line(size = 1) +
      labs(title = "Neuroinflammation & HPA Axis",
           x = "Day", y = "Normalised Index") +
      theme_bw(base_size = 13) +
      scale_color_manual(values = c(IL6="tomato", Cortisol="darkorange", GSK3="darkred"))
  })

  ## ---- TAB 4: Clinical Endpoints ----
  output$ymrs_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, YMRS)) +
      geom_line(color = "firebrick", size = 1.2) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = input$ymrs_init * 0.5, linetype = "dashed", color = "blue") +
      annotate("text", x = max(df$day) * 0.1, y = 13.5, label = "Remission", size = 3.5) +
      labs(title = "YMRS – Mania Score", x = "Day", y = "YMRS") +
      theme_bw(base_size = 13) + ylim(0, max(df$YMRS) + 2)
  })

  output$madrs_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, MADRS)) +
      geom_line(color = "steelblue", size = 1.2) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(df$day) * 0.1, y = 13.5, label = "Remission", size = 3.5) +
      labs(title = "MADRS – Depression Score", x = "Day", y = "MADRS") +
      theme_bw(base_size = 13) + ylim(0, max(df$MADRS) + 2)
  })

  output$gaf_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, GAF)) +
      geom_line(color = "green4", size = 1.2) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "gray40") +
      labs(title = "Global Assessment of Functioning", x = "Day", y = "GAF Score") +
      theme_bw(base_size = 13) + ylim(0, 100)
  })

  output$vbox_ymrs_resp <- renderValueBox({
    df <- sim_result()
    resp <- tail(df$YMRS_resp, 1)
    valueBox(if (resp == 1) "Yes ✓" else "No", "YMRS Response (≥50%↓)",
             color = if (resp == 1) "green" else "red", icon = icon("arrow-down"))
  })

  output$vbox_madrs_resp <- renderValueBox({
    df <- sim_result()
    resp <- tail(df$MADRS_resp, 1)
    valueBox(if (resp == 1) "Yes ✓" else "No", "MADRS Response (≥50%↓)",
             color = if (resp == 1) "green" else "orange", icon = icon("arrow-down"))
  })

  output$vbox_remission <- renderValueBox({
    df <- sim_result()
    remit <- tail(df$MADRS_remit, 1)
    valueBox(if (remit == 1) "Yes ✓" else "No", "MADRS Remission (≤12)",
             color = if (remit == 1) "blue" else "yellow", icon = icon("check-circle"))
  })

  output$vbox_gaf_end <- renderValueBox({
    df <- sim_result()
    gaf <- round(tail(df$GAF, 1), 1)
    valueBox(gaf, "End GAF Score",
             color = if (gaf >= 60) "green" else if (gaf >= 40) "yellow" else "red",
             icon = icon("user-check"))
  })

  ## ---- TAB 5: Scenario Comparison ----
  compare_data <- eventReactive(input$run_compare, {
    scen_defs <- list(
      sc1 = list(label="Li 900mg/d",    Li=8.1, freq_Li=8,  VPA=0,   QTP=0,   LTG=0,   YMRS=25, MADRS=5),
      sc2 = list(label="VPA 1000mg/d",  Li=0,   freq_Li=8,  VPA=500, QTP=0,   LTG=0,   YMRS=25, MADRS=5),
      sc3 = list(label="QTP 600mg/d",   Li=0,   freq_Li=8,  VPA=0,   QTP=600, LTG=0,   YMRS=25, MADRS=5),
      sc4 = list(label="QTP 300 BDdep", Li=0,   freq_Li=24, VPA=0,   QTP=300, LTG=0,   YMRS=5,  MADRS=30),
      sc5 = list(label="Li+QTP BDdep",  Li=8.1, freq_Li=8,  VPA=0,   QTP=300, LTG=0,   YMRS=5,  MADRS=30),
      sc6 = list(label="LTG 200 BDdep", Li=0,   freq_Li=24, VPA=0,   QTP=0,   LTG=200, YMRS=5,  MADRS=30),
      sc7 = list(label="Li maint",       Li=8.1, freq_Li=8,  VPA=0,   QTP=0,   LTG=0,   YMRS=10, MADRS=10),
      sc8 = list(label="VPA maint",      Li=0,   freq_Li=8,  VPA=500, QTP=0,   LTG=0,   YMRS=10, MADRS=10)
    )
    sel_scen <- c(input$mania_scen, input$dep_scen, input$maint_scen)
    lapply(sel_scen, function(sc) {
      d <- scen_defs[[sc]]
      res <- run_bd_model(dose_Li=d$Li, freq_Li=d$freq_Li,
                          dose_VPA=d$VPA, dose_QTP=d$QTP, dose_LTG=d$LTG,
                          YMRS_init=d$YMRS, MADRS_init=d$MADRS,
                          duration_d = input$duration_d, dt = 1)
      res$Scenario <- d$label
      res
    })
  })

  output$compare_ymrs <- renderPlot({
    df_list <- compare_data()
    df <- bind_rows(df_list)
    ggplot(df, aes(day, YMRS, color = Scenario)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      labs(title = "YMRS Comparison", x = "Day", y = "YMRS") +
      theme_bw(base_size = 13)
  })

  output$compare_madrs <- renderPlot({
    df_list <- compare_data()
    df <- bind_rows(df_list)
    ggplot(df, aes(day, MADRS, color = Scenario)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      labs(title = "MADRS Comparison", x = "Day", y = "MADRS") +
      theme_bw(base_size = 13)
  })

  output$compare_tbl <- DT::renderDataTable({
    df_list <- compare_data()
    tbl <- lapply(df_list, function(df) {
      last_d <- tail(df, 24)
      tibble(
        Scenario      = unique(df$Scenario),
        `End YMRS`    = round(tail(df$YMRS, 1), 1),
        `End MADRS`   = round(tail(df$MADRS, 1), 1),
        `End GAF`     = round(tail(df$GAF, 1), 1),
        `YMRS resp`   = ifelse(tail(df$YMRS_resp, 1) == 1, "Yes", "No"),
        `MADRS resp`  = ifelse(tail(df$MADRS_resp, 1) == 1, "Yes", "No"),
        `MADRS remit` = ifelse(tail(df$MADRS_remit, 1) == 1, "Yes", "No"),
        `Li Css (mEq/L)` = round(mean(last_d$Li_mEqL), 2)
      )
    })
    bind_rows(tbl)
  }, options = list(pageLength = 10))

  ## ---- TAB 6: Safety Monitor ----
  output$vbox_li_safe <- renderValueBox({
    df <- sim_result()
    li_max <- max(df$Li_mEqL)
    status <- if (li_max > 1.5) "TOXIC" else if (li_max > 1.2) "High – monitor" else "Safe"
    valueBox(sprintf("%.2f mEq/L", li_max), "Peak Lithium",
             color = if (li_max > 1.5) "red" else if (li_max > 1.2) "orange" else "green",
             icon = icon("tint"))
  })

  output$vbox_vpa_safe <- renderValueBox({
    df <- sim_result()
    vpa_max <- max(df$VPA_ugmL)
    valueBox(sprintf("%.0f μg/mL", vpa_max), "Peak VPA",
             color = if (vpa_max > 125) "red" else if (vpa_max > 100) "orange" else "green",
             icon = icon("pills"))
  })

  output$vbox_qtc <- renderValueBox({
    df <- sim_result()
    qtp_max <- max(df$QTP_ngmL)
    # Simplified QTc risk estimate
    qtc_est <- 420 + qtp_max * 0.04  # rough linear approximation
    valueBox(sprintf("~%.0f ms", qtc_est), "Estimated QTc",
             color = if (qtc_est > 500) "red" else if (qtc_est > 470) "orange" else "green",
             icon = icon("heartbeat"))
  })

  output$li_safety_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, Li_mEqL)) +
      geom_line(color = "steelblue", size = 1.1) +
      geom_ribbon(aes(ymin = 0.6, ymax = 1.2), alpha = 0.15, fill = "green3") +
      geom_hline(yintercept = 1.5, color = "red",    linetype = "dashed") +
      geom_hline(yintercept = 0.6, color = "orange", linetype = "dashed") +
      annotate("text", x = max(df$day)*0.6, y = 1.55, label="Toxic", color="red3", size=3.5) +
      labs(title = "Lithium Safety", x = "Day", y = "Li Concentration (mEq/L)") +
      theme_bw(base_size = 13)
  })

  output$metabolic_plot <- renderPlot({
    df <- sim_result()
    ggplot(df, aes(day, Wt)) +
      geom_line(color = "coral", size = 1.1) +
      labs(title = "Body Weight Change (QTP-driven)",
           x = "Day", y = "Weight change (kg)") +
      theme_bw(base_size = 13)
  })

  output$safety_tbl <- renderTable({
    df <- sim_result()
    last <- tail(df, 1)
    tibble(
      Parameter        = c("Lithium (mEq/L)", "Valproate (μg/mL)", "Quetiapine (ng/mL)",
                           "BDNF index", "IL-6 index", "Weight Δ (kg)"),
      "End value"      = round(c(last$Li_mEqL, last$VPA_ugmL, last$QTP_ngmL,
                                 last$BDNF, last$IL6, last$Wt), 2),
      "Normal range"   = c("0.6–1.2", "50–100", "100–400", "~1.0", "~1.0", "<3 kg"),
      "Action"         = c(
        if (last$Li_mEqL > 1.2) "Reduce dose / hydrate" else "OK",
        if (last$VPA_ugmL > 100) "Monitor LFTs / platelets" else "OK",
        if (last$QTP_ngmL > 400) "ECG – QTc" else "OK",
        if (last$BDNF < 0.8) "Consider BDNF-enhancing add-on" else "OK",
        if (last$IL6 > 1.2) "Consider anti-inflammatory monitoring" else "OK",
        if (last$Wt > 5) "Metabolic monitoring (fasting glucose, lipids)" else "OK"
      )
    )
  })
}

## ============================================================================
## Run App
## ============================================================================
shinyApp(ui, server)
