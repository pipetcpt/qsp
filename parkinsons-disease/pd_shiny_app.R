## =============================================================================
## Parkinson's Disease QSP — Shiny Interactive Dashboard
## 6 Tabs: Patient Profile · PK · PD Biomarkers · Motor Endpoints ·
##         Scenario Comparison · Neuroprotection Explorer
## =============================================================================
## Run: shiny::runApp("pd_shiny_app.R")
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## Minimal inline ODE solver (avoids mrgsolve dependency for Shiny demo)
## Euler integration of simplified PD QSP ODEs
## ─────────────────────────────────────────────────────────────────────────────

pd_ode_euler <- function(params, times) {
  # Unpack parameters
  kSN   <- params$kSN_death
  kProt <- params$kSN_protect
  kAO   <- params$kaSyn_nuc
  kAF   <- params$kaSyn_elong
  kAcl  <- params$kaSyn_clear
  aSyn0 <- params$aSyn0
  kTH   <- params$kTH
  kMAOB <- params$kMAOB * (1 - params$MAOBinh)
  kCOMT <- params$kCOMT * (1 - params$COMTinh)
  kDAT  <- params$kDAT_reup
  EC50  <- params$EC50_D2
  Emax  <- params$Emax_D2
  kUP   <- params$kUPDRS_prog
  LDdose  <- params$LD_daily_dose  # total mg/day levodopa
  PRAMdose<- params$PRAM_daily_dose
  RASdose <- params$RAS_daily_dose
  # Bioavailability/brain penetration
  LD_brain_frac <- 0.05   # fraction reaching brain
  PRAM_DA_equiv <- 0.3    # DA-equivalent per mg pramipexole

  dt <- times[2] - times[1]
  n  <- length(times)

  # State variables
  SNpc   <- 1.0; ASyn_M <- aSyn0; ASyn_O <- 0.05; ASyn_F <- 0.001
  ROS    <- 0.1;  NEUROINF <- 0.05
  DA_brain <- 10.0; DA_syn <- 2.0
  UPDRS_III <- 5.0; LID_risk <- 0.0

  results <- matrix(0, nrow = n, ncol = 14)
  colnames(results) <- c(
    "time","SNpc","ASyn_O","ASyn_F","ROS","NEUROINF",
    "DA_syn","DA_brain","UPDRS_III","LID_risk",
    "D2R_stim","GPi_output","MotorDrive","DA_eff"
  )

  for (i in seq_along(times)) {
    # Exogenous DA sources
    DA_exo_LD   <- (LDdose / 1000) * LD_brain_frac
    DA_exo_PRAM <- PRAMdose * PRAM_DA_equiv

    DA_eff <- DA_syn + DA_exo_LD + DA_exo_PRAM
    D2R_stim   <- Emax * DA_eff / (DA_eff + EC50)
    GPi_output <- max(0, 1.0 + 0.4 * (1 - D2R_stim * 0.5) - D2R_stim * 0.6)
    MotorDrive <- 1.0 / (1.0 + 0.8 * GPi_output)

    results[i,] <- c(
      times[i], SNpc, ASyn_O, ASyn_F, ROS, NEUROINF,
      DA_syn, DA_brain, UPDRS_III, LID_risk,
      D2R_stim, GPi_output, MotorDrive, DA_eff
    )

    # ODE step (simplified Euler)
    dASyn_M <- -kAO * ASyn_M^2 - kAF * ASyn_M * ASyn_F + kAcl * (aSyn0 - ASyn_M)
    dASyn_O <- kAO * ASyn_M^2 - kAF * ASyn_M * ASyn_O - kAcl * ASyn_O
    dASyn_F <- kAF * ASyn_M * ASyn_O - 0.001 * ASyn_F
    dROS    <- 0.05 * (1 + ASyn_O / 5) * (1 - SNpc) - 0.5 * ROS
    dNEUROINF <- 0.03 * (ASyn_O + ROS) - 0.1 * NEUROINF

    dSNpc   <- -kSN * SNpc * (1 + ASyn_O/10 + ROS/5 + NEUROINF/3) +
               kProt * SNpc * (1 - SNpc)
    if (RASdose > 0) dSNpc <- dSNpc + 0.00002 * SNpc  # rasagiline neuroprotection

    dDA_brain <- kTH * SNpc + DA_exo_LD - (kMAOB + kCOMT) * DA_brain
    dDA_syn   <- 0.5 * DA_brain - kDAT * DA_syn - (kMAOB + kCOMT) * DA_syn / 2 +
                 DA_exo_PRAM * 0.1

    # UPDRS driven by neuronal loss and GPi overactivation
    UPDRS_target <- 5 + 60 * (1 - SNpc) / 0.6
    UPDRS_benefit <- 20 * D2R_stim * MotorDrive
    dUPDRS <- kUP * (1 - SNpc) * 0.1 - 0.001 * UPDRS_benefit
    dLID   <- 0.0003 * (DA_exo_LD^2) * (1 - SNpc)

    # Update states
    ASyn_M    <- max(0, ASyn_M + dASyn_M * dt)
    ASyn_O    <- max(0, ASyn_O + dASyn_O * dt)
    ASyn_F    <- max(0, ASyn_F + dASyn_F * dt)
    ROS       <- max(0, ROS + dROS * dt)
    NEUROINF  <- max(0, NEUROINF + dNEUROINF * dt)
    SNpc      <- max(0.01, min(1, SNpc + dSNpc * dt))
    DA_brain  <- max(0, DA_brain + dDA_brain * dt)
    DA_syn    <- max(0, DA_syn + dDA_syn * dt)
    UPDRS_III <- max(0, UPDRS_III + dUPDRS * dt)
    LID_risk  <- max(0, LID_risk + dLID * dt)
  }

  as.data.frame(results)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = "Parkinson's Disease QSP Dashboard",
    titleWidth = 380
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",        tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Pharmacokinetics (PK)",  tabName = "tab_pk",       icon = icon("chart-line")),
      menuItem("PD Biomarkers",          tabName = "tab_pd",       icon = icon("brain")),
      menuItem("Motor Endpoints",        tabName = "tab_motor",    icon = icon("walking")),
      menuItem("Scenario Comparison",    tabName = "tab_scen",     icon = icon("vials")),
      menuItem("Neuroprotection",        tabName = "tab_neuro",    icon = icon("shield-alt"))
    ),
    hr(),
    h4("Global Parameters", style = "margin-left:10px; color:#CE93D8"),
    sliderInput("kSN_death",  "Neuron Death Rate",  0.0001, 0.0006, 0.00019, step=0.00001, width="90%"),
    sliderInput("aSyn0",      "Baseline α-Syn (nM)", 1, 20, 5, step=0.5, width="90%"),
    sliderInput("sim_years",  "Simulation Years",   10, 30, 25, step=1, width="90%"),
    hr(),
    h4("Treatment (daily dose, mg)", style = "margin-left:10px; color:#80CBC4"),
    sliderInput("LD_daily_dose",   "Levodopa (mg/day)",     0, 1500, 750, step=50, width="90%"),
    sliderInput("PRAM_daily_dose", "Pramipexole (mg/day)",  0, 4.5,  1.5, step=0.25, width="90%"),
    sliderInput("RAS_daily_dose",  "Rasagiline (mg/day)",   0, 1,    0,   step=0.5, width="90%"),
    sliderInput("MAOBinh",         "MAO-B inhibition (0-1)", 0, 1, 0, step=0.05, width="90%"),
    sliderInput("COMTinh",         "COMT inhibition (0-1)",  0, 1, 0, step=0.05, width="90%"),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 style = "margin:10px; background:#7B1FA2; color:white; width:88%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #FAFAFA; }
      .box { border-radius: 8px; }
    "))),
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient & Disease Parameters", width = 12, status = "primary",
              fluidRow(
                column(4,
                  h4("Demographics"),
                  numericInput("age_dx", "Age at Diagnosis (yr)", 62, 30, 90),
                  selectInput("sex", "Sex", c("Male", "Female")),
                  numericInput("weight_kg", "Weight (kg)", 70, 40, 130),
                  selectInput("pd_type", "PD Subtype",
                              c("Tremor-dominant", "Akinetic-rigid", "Mixed/Other"))
                ),
                column(4,
                  h4("Disease Stage"),
                  selectInput("hoehn_yahr", "Hoehn & Yahr Stage",
                              c("1 – Unilateral only",
                                "2 – Bilateral, no balance",
                                "3 – Mild bilateral, impaired balance",
                                "4 – Severe disability, can still walk",
                                "5 – Wheelchair/bedridden")),
                  numericInput("dx_delay_yr", "Years from Sx Onset to Dx", 2, 0, 10),
                  numericInput("updrs_baseline", "UPDRS-III at Baseline", 28, 0, 100)
                ),
                column(4,
                  h4("Genetic Risk"),
                  checkboxInput("has_LRRK2", "LRRK2 G2019S mutation", FALSE),
                  checkboxInput("has_GBA", "GBA variant", FALSE),
                  checkboxInput("has_SNCA_mult", "SNCA multiplication", FALSE),
                  checkboxInput("has_PINK1", "PINK1/Parkin mutation", FALSE),
                  h4("Comorbidities"),
                  checkboxInput("has_depression", "Depression", FALSE),
                  checkboxInput("has_RBD", "REM Sleep Behavior Disorder", FALSE),
                  checkboxInput("has_anosmia", "Hyposmia/Anosmia", FALSE)
                )
              )
          )
        ),
        fluidRow(
          valueBoxOutput("vb_snpc"),
          valueBoxOutput("vb_updrs"),
          valueBoxOutput("vb_da_syn")
        ),
        fluidRow(
          box(title = "Braak Staging & Pathological Progression", width = 12, status = "warning",
              plotOutput("plot_braak", height = "220px"))
        )
      ),

      ## ── Tab 2: Pharmacokinetics ─────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Levodopa PK Settings", width = 4, status = "info",
              sliderInput("LD_single_dose", "Single L-DOPA Dose (mg)", 50, 500, 250),
              sliderInput("LD_freq_h",      "Dosing Interval (h)", 3, 12, 8),
              selectInput("LD_formulation", "Formulation",
                          c("Immediate Release (Sinemet)",
                            "Extended Release (Rytary)",
                            "Intestinal Gel (Duopa)")),
              checkboxInput("with_carbidopa", "Add Carbidopa (DDC inhibitor)", TRUE),
              checkboxInput("high_protein_meal", "High-protein meal (LAA competition)", FALSE)
          ),
          box(title = "Levodopa Plasma PK Profile", width = 8, status = "info",
              plotlyOutput("plot_LD_pk", height = "350px"))
        ),
        fluidRow(
          box(title = "Pramipexole PK", width = 6, status = "success",
              sliderInput("PRAM_single", "Single Pramipexole Dose (mg)", 0.125, 1.5, 0.5),
              plotlyOutput("plot_PRAM_pk", height = "250px")
          ),
          box(title = "Rasagiline MAO-B Inhibition", width = 6, status = "warning",
              plotlyOutput("plot_RAS_MAOB", height = "250px")
          )
        )
      ),

      ## ── Tab 3: PD Biomarkers ────────────────────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "α-Synuclein Dynamics", width = 6, status = "danger",
              plotlyOutput("plot_asyn", height = "300px")),
          box(title = "Neuroinflammation & Oxidative Stress", width = 6, status = "warning",
              plotlyOutput("plot_neuro_inflam", height = "300px"))
        ),
        fluidRow(
          box(title = "SNpc Neuron Survival", width = 6, status = "primary",
              plotlyOutput("plot_snpc", height = "300px")),
          box(title = "Synaptic Dopamine", width = 6, status = "success",
              plotlyOutput("plot_da_syn", height = "300px"))
        )
      ),

      ## ── Tab 4: Motor Endpoints ──────────────────────────────────────────
      tabItem(tabName = "tab_motor",
        fluidRow(
          box(title = "UPDRS-III Motor Score", width = 6, status = "primary",
              plotlyOutput("plot_updrs", height = "320px")),
          box(title = "LID Risk Index", width = 6, status = "danger",
              plotlyOutput("plot_lid", height = "320px"))
        ),
        fluidRow(
          box(title = "GPi Output & Motor Drive", width = 6, status = "info",
              plotlyOutput("plot_gpi", height = "280px")),
          box(title = "D2R Receptor Stimulation", width = 6, status = "success",
              plotlyOutput("plot_d2r", height = "280px"))
        )
      ),

      ## ── Tab 5: Scenario Comparison ──────────────────────────────────────
      tabItem(tabName = "tab_scen",
        fluidRow(
          box(title = "Select Scenarios to Compare", width = 12, status = "primary",
              checkboxGroupInput("compare_scens", NULL,
                choices = c(
                  "No Treatment"                   = "s1",
                  "Levodopa 250mg TID"             = "s2",
                  "Pramipexole 0.75mg TID"         = "s3",
                  "Rasagiline 1mg QD"              = "s4",
                  "Levodopa + Entacapone"          = "s5",
                  "Triple Therapy (LD+PRAM+RAS)"   = "s6",
                  "Continuous Delivery (LD CR+PRAM)"= "s7"
                ),
                selected = c("s1","s2","s3","s6"),
                inline = TRUE
              )
          )
        ),
        fluidRow(
          box(title = "Neuron Survival by Scenario", width = 6, status = "primary",
              plotlyOutput("plot_scen_snpc", height = "320px")),
          box(title = "UPDRS-III by Scenario", width = 6, status = "info",
              plotlyOutput("plot_scen_updrs", height = "320px"))
        ),
        fluidRow(
          box(title = "LID Risk by Scenario", width = 6, status = "danger",
              plotlyOutput("plot_scen_lid", height = "280px")),
          box(title = "Scenario Outcome Table", width = 6, status = "success",
              DTOutput("table_scen"))
        )
      ),

      ## ── Tab 6: Neuroprotection Explorer ────────────────────────────────
      tabItem(tabName = "tab_neuro",
        fluidRow(
          box(title = "Disease Modification Strategies", width = 4, status = "warning",
              h4("Intervention Parameters"),
              sliderInput("np_aSyn_clear", "α-Syn Clearance Boost (fold)", 1, 5, 1, step=0.5),
              sliderInput("np_LRRK2_inh",  "LRRK2 Inhibition (%)", 0, 100, 0),
              sliderInput("np_GDNF",        "GDNF/Neurotrophic Support (%)", 0, 100, 0),
              sliderInput("np_GLP1",        "GLP-1 (Exenatide) Effect (%)", 0, 100, 0),
              sliderInput("np_GBA_boost",   "GBA Activator Effect (%)", 0, 100, 0),
              checkboxInput("np_DBS",       "DBS (STN/GPi)", FALSE),
              actionButton("run_np", "Calculate Neuroprotection",
                           style="background:#E65100;color:white;margin-top:10px")
          ),
          box(title = "Neuroprotection Outcome", width = 8, status = "warning",
              plotlyOutput("plot_np_snpc", height = "350px")
          )
        ),
        fluidRow(
          box(title = "Braak Stage Progression with Intervention", width = 6, status = "danger",
              plotlyOutput("plot_np_braak", height = "280px")),
          box(title = "Time to H&Y Stage 3 (years from Dx)", width = 6, status = "primary",
              DTOutput("table_np_endpoints"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## Server
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Base simulation (reactive) ─────────────────────────────────────────
  base_sim <- eventReactive(input$run_sim, {
    times <- seq(0, input$sim_years * 365, by = 7)
    params <- list(
      kSN_death    = input$kSN_death,
      kSN_protect  = 0.0001,
      kaSyn_nuc    = 0.0005,
      kaSyn_elong  = 0.002,
      kaSyn_clear  = 0.15,
      aSyn0        = input$aSyn0,
      kTH          = 0.8,
      kMAOB        = 0.6,
      kCOMT        = 0.2,
      kDAT_reup    = 1.5,
      EC50_D2      = 0.15,
      Emax_D2      = 1.0,
      kUPDRS_prog  = 0.002,
      LD_daily_dose  = input$LD_daily_dose,
      PRAM_daily_dose= input$PRAM_daily_dose,
      RAS_daily_dose = input$RAS_daily_dose,
      MAOBinh        = input$MAOBinh,
      COMTinh        = input$COMTinh
    )
    pd_ode_euler(params, times) %>% mutate(time_years = time / 365)
  }, ignoreNULL = FALSE)

  ## ── Scenario simulations ───────────────────────────────────────────────
  make_scen_sim <- function(LD = 0, PRAM = 0, RAS = 0, MAOBi = 0, COMTi = 0,
                             aSyn_boost = 1) {
    times <- seq(0, input$sim_years * 365, by = 14)
    params <- list(
      kSN_death = input$kSN_death, kSN_protect = 0.0001,
      kaSyn_nuc = 0.0005, kaSyn_elong = 0.002,
      kaSyn_clear = 0.15 * aSyn_boost, aSyn0 = input$aSyn0,
      kTH = 0.8, kMAOB = 0.6, kCOMT = 0.2, kDAT_reup = 1.5,
      EC50_D2 = 0.15, Emax_D2 = 1.0, kUPDRS_prog = 0.002,
      LD_daily_dose = LD, PRAM_daily_dose = PRAM, RAS_daily_dose = RAS,
      MAOBinh = MAOBi, COMTinh = COMTi
    )
    pd_ode_euler(params, times) %>% mutate(time_years = time/365)
  }

  scen_sims <- reactive({
    list(
      s1 = make_scen_sim() %>% mutate(scenario = "No Treatment"),
      s2 = make_scen_sim(LD = 750) %>% mutate(scenario = "Levodopa 250mg TID"),
      s3 = make_scen_sim(PRAM = 2.25) %>% mutate(scenario = "Pramipexole 0.75mg TID"),
      s4 = make_scen_sim(RAS = 1, MAOBi = 0.98) %>% mutate(scenario = "Rasagiline 1mg QD"),
      s5 = make_scen_sim(LD = 750, COMTi = 0.65) %>% mutate(scenario = "Levodopa + Entacapone"),
      s6 = make_scen_sim(LD = 750, PRAM = 1.5, RAS = 1, MAOBi = 0.98) %>%
           mutate(scenario = "Triple Therapy"),
      s7 = make_scen_sim(LD = 600, PRAM = 1.5) %>% mutate(scenario = "Continuous Delivery")
    )
  })

  ## ── Value Boxes ────────────────────────────────────────────────────────
  output$vb_snpc <- renderValueBox({
    sim <- base_sim()
    last <- tail(sim, 1)
    valueBox(
      paste0(round(last$SNpc * 100, 1), "%"),
      "Surviving SNpc Neurons",
      icon = icon("brain"),
      color = if (last$SNpc > 0.6) "green" else if (last$SNpc > 0.4) "yellow" else "red"
    )
  })

  output$vb_updrs <- renderValueBox({
    sim <- base_sim()
    last <- tail(sim, 1)
    val <- round(last$UPDRS_III, 1)
    valueBox(val, "UPDRS-III Score", icon = icon("chart-bar"),
             color = if (val < 20) "green" else if (val < 40) "yellow" else "red")
  })

  output$vb_da_syn <- renderValueBox({
    sim <- base_sim()
    last <- tail(sim, 1)
    valueBox(
      round(last$DA_syn, 3),
      "Synaptic DA (nM)",
      icon = icon("atom"),
      color = "purple"
    )
  })

  ## ── Braak staging plot ─────────────────────────────────────────────────
  output$plot_braak <- renderPlot({
    stages <- data.frame(
      stage = factor(1:6),
      region = c("Olfactory bulb\n& ENS", "Medulla\n(dorsal vagus)",
                  "Pons\n(locus coeruleus)", "Midbrain\n(SNpc)",
                  "Limbic cortex\n(temporal)", "Neocortex\n(prefrontal)"),
      year_approx = c(-5, -2, 0, 3, 8, 15),
      color_val = c("#FFCDD2","#EF9A9A","#E57373","#EF5350","#E53935","#B71C1C")
    )
    ggplot(stages, aes(year_approx, stage, fill = color_val)) +
      geom_tile(width = 5, height = 0.8, color = "white", linewidth = 1.2) +
      geom_text(aes(label = paste0("Stage ", stage, "\n", region)),
                size = 3, fontface = "bold") +
      scale_fill_identity() +
      geom_vline(xintercept = 0, color = "#3F51B5", linetype = "dashed", linewidth = 1.2) +
      annotate("text", x = 0.3, y = 6.5, label = "Symptom\nOnset", color = "#3F51B5", size = 3) +
      labs(title = "Braak Staging of α-Synuclein Pathology",
           x = "Years relative to motor symptom onset", y = "Braak Stage") +
      theme_minimal(base_size = 11)
  })

  ## ── PK plots ───────────────────────────────────────────────────────────
  output$plot_LD_pk <- renderPlotly({
    # Simple 1-compartment PK for demo
    Ka <- 1.2; Vd <- 63; CL <- 98  # for 70kg
    if (input$LD_formulation == "Extended Release (Rytary)") Ka <- 0.4
    if (input$high_protein_meal) Ka <- Ka * 0.6

    t_sim <- seq(0, 24, by = 0.1)
    dose  <- input$LD_single_dose
    freq  <- input$LD_freq_h
    dosing_times <- seq(0, 24, by = freq)

    cp <- numeric(length(t_sim))
    for (td in dosing_times) {
      t_after <- pmax(t_sim - td, 0)
      cp <- cp + (dose * Ka) / (Vd * (Ka - CL/Vd)) *
            (exp(-(CL/Vd) * t_after) - exp(-Ka * t_after))
    }

    df <- data.frame(time = t_sim, Cp = pmax(cp, 0))
    thr_on  <- 200; thr_lih <- 800
    p <- ggplot(df, aes(time, Cp)) +
      geom_area(fill = "#B2EBF2", alpha = 0.5) +
      geom_line(color = "#0097A7", linewidth = 1.2) +
      geom_hline(yintercept = thr_on, color = "#43A047", linetype = "dashed") +
      geom_hline(yintercept = thr_lih, color = "#E91E63", linetype = "dashed") +
      annotate("text", x = 1, y = thr_on + 30, label = "Therapeutic threshold",
               color = "#43A047", size = 3) +
      annotate("text", x = 1, y = thr_lih + 30, label = "Dyskinesia risk",
               color = "#E91E63", size = 3) +
      labs(title = "L-DOPA Plasma Concentration (µg/L)",
           x = "Time (h)", y = "Cp (µg/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_PRAM_pk <- renderPlotly({
    Ka <- 0.8; Vd <- 490; CL <- 28
    t_sim <- seq(0, 24, by = 0.1)
    dose  <- input$PRAM_single * 1000  # µg
    freq  <- 8
    cp <- numeric(length(t_sim))
    for (td in seq(0, 24, by = freq)) {
      t_after <- pmax(t_sim - td, 0)
      cp <- cp + (dose * Ka) / (Vd * (Ka - CL/Vd)) *
            (exp(-(CL/Vd) * t_after) - exp(-Ka * t_after))
    }
    df <- data.frame(time = t_sim, Cp = pmax(cp, 0))
    p <- ggplot(df, aes(time, Cp)) +
      geom_line(color = "#43A047", linewidth = 1.2) +
      labs(title = "Pramipexole Plasma (ng/L)", x = "Time (h)", y = "Cp (ng/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_RAS_MAOB <- renderPlotly({
    # Rasagiline irreversible MAO-B inhibition
    t_days <- seq(0, 14, by = 0.2)
    dose_mg <- input$RAS_daily_dose
    Imax <- 1.0; IC50_approx <- 0.98  # 1mg dose → ~98% inhibition
    inhibition_steady <- Imax * dose_mg / (dose_mg + 0.02)  # simplified
    recovery_half_life_days <- 7  # MAO-B enzyme recovery

    # Steady-state scenario: washout after 7 days
    maob_inhib <- ifelse(t_days <= 7,
                         inhibition_steady * (1 - exp(-0.2 * t_days)),
                         inhibition_steady * exp(-(t_days - 7) * log(2) / recovery_half_life_days))

    df <- data.frame(time = t_days, inhibition = maob_inhib * 100)
    p <- ggplot(df, aes(time, inhibition)) +
      geom_line(color = "#9C27B0", linewidth = 1.2) +
      geom_vline(xintercept = 7, color = "red", linetype = "dashed") +
      annotate("text", x = 7.2, y = 50, label = "Drug\nwithdrawn", color = "red", size = 3) +
      ylim(0, 100) +
      labs(title = "MAO-B Inhibition by Rasagiline (%)",
           x = "Time (days)", y = "MAO-B Inhibition (%)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── PD Biomarker plots ─────────────────────────────────────────────────
  output$plot_asyn <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years)) +
      geom_line(aes(y = ASyn_O, color = "Oligomers"), linewidth = 1) +
      geom_line(aes(y = ASyn_F * 10, color = "Fibrils (×10)"), linewidth = 1) +
      scale_color_manual(values = c("Oligomers" = "#E91E63", "Fibrils (×10)" = "#880E4F")) +
      labs(title = "α-Synuclein Aggregation", x = "Time (years)", y = "Conc. (nM)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_neuro_inflam <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years)) +
      geom_line(aes(y = NEUROINF, color = "Neuroinflammation"), linewidth = 1) +
      geom_line(aes(y = ROS, color = "ROS"), linewidth = 1) +
      scale_color_manual(values = c("Neuroinflammation" = "#E65100", "ROS" = "#EF9A9A")) +
      labs(title = "Neuroinflammation & Oxidative Stress",
           x = "Time (years)", y = "Index (a.u.)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_snpc <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, SNpc * 100)) +
      geom_area(fill = "#7986CB", alpha = 0.3) +
      geom_line(color = "#3F51B5", linewidth = 1.2) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "#E91E63") +
      annotate("text", x = 2, y = 38, label = "Sx onset threshold (~60% loss)",
               color = "#E91E63", size = 3) +
      labs(title = "SNpc DA Neuron Survival", x = "Time (years)", y = "Surviving Neurons (%)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_da_syn <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, DA_syn)) +
      geom_line(color = "#9C27B0", linewidth = 1.2) +
      labs(title = "Synaptic Dopamine Concentration", x = "Time (years)", y = "DA (nM)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Motor endpoint plots ───────────────────────────────────────────────
  output$plot_updrs <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, UPDRS_III)) +
      geom_area(fill = "#1565C0", alpha = 0.2) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      labs(title = "UPDRS-III Motor Score", x = "Time (years)", y = "UPDRS-III (0-108)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_lid <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, LID_risk)) +
      geom_area(fill = "#E91E63", alpha = 0.2) +
      geom_line(color = "#E91E63", linewidth = 1.2) +
      labs(title = "Levodopa-Induced Dyskinesia Risk Index",
           x = "Time (years)", y = "Cumulative LID Risk (a.u.)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_gpi <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, GPi_output)) +
      geom_line(color = "#FF8F00", linewidth = 1.2) +
      labs(title = "GPi Output (PD: overactive → suppresses thalamus)",
           x = "Time (years)", y = "GPi Activity (a.u.)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_d2r <- renderPlotly({
    sim <- base_sim()
    p <- ggplot(sim, aes(time_years, D2R_stim * 100)) +
      geom_line(color = "#43A047", linewidth = 1.2) +
      labs(title = "D2 Receptor Stimulation (%)",
           x = "Time (years)", y = "D2R Stim (% max)") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Scenario comparison plots ──────────────────────────────────────────
  scen_data_combined <- reactive({
    sims <- scen_sims()
    sel  <- input$compare_scens
    bind_rows(sims[sel])
  })

  output$plot_scen_snpc <- renderPlotly({
    df <- scen_data_combined()
    p <- ggplot(df, aes(time_years, SNpc * 100, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = "Neuron Survival by Scenario", x = "Time (years)",
           y = "Surviving Neurons (%)", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_scen_updrs <- renderPlotly({
    df <- scen_data_combined()
    p <- ggplot(df, aes(time_years, UPDRS_III, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = "UPDRS-III by Scenario", x = "Time (years)",
           y = "UPDRS-III", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$plot_scen_lid <- renderPlotly({
    df <- scen_data_combined()
    p <- ggplot(df, aes(time_years, LID_risk, color = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = "LID Risk", x = "Time (years)", y = "LID Risk (a.u.)", color = "") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$table_scen <- renderDT({
    sims <- scen_sims()
    sel  <- input$compare_scens
    bind_rows(sims[sel]) %>%
      filter(round(time_years, 0) %in% c(10, 15, 20, 25)) %>%
      group_by(scenario, Year = round(time_years, 0)) %>%
      summarise(
        SNpc_pct   = round(mean(SNpc) * 100, 1),
        UPDRS_III  = round(mean(UPDRS_III), 1),
        LID_risk   = round(mean(LID_risk), 4),
        DA_syn     = round(mean(DA_syn), 3),
        .groups    = "drop"
      ) %>%
      datatable(options = list(pageLength = 12), rownames = FALSE)
  })

  ## ── Neuroprotection tab ────────────────────────────────────────────────
  np_sim <- eventReactive(input$run_np, {
    times <- seq(0, input$sim_years * 365, by = 14)
    base_params <- list(
      kSN_death = input$kSN_death,
      kSN_protect = 0.0001 * (1 + input$np_GDNF / 100 * 2),
      kaSyn_nuc = 0.0005 * (1 - input$np_LRRK2_inh / 200),
      kaSyn_elong = 0.002,
      kaSyn_clear = 0.15 * input$np_aSyn_clear,
      aSyn0 = input$aSyn0,
      kTH = 0.8, kMAOB = 0.6, kCOMT = 0.2, kDAT_reup = 1.5,
      EC50_D2 = 0.15, Emax_D2 = 1.0, kUPDRS_prog = 0.002,
      LD_daily_dose  = input$LD_daily_dose,
      PRAM_daily_dose= input$PRAM_daily_dose,
      RAS_daily_dose = input$RAS_daily_dose,
      MAOBinh = input$MAOBinh,
      COMTinh = input$COMTinh
    )
    intervention_sim <- pd_ode_euler(base_params, times) %>%
      mutate(time_years = time/365, group = "With Neuroprotection")

    baseline_params <- base_params
    baseline_params$kSN_protect  <- 0.0001
    baseline_params$kaSyn_nuc    <- 0.0005
    baseline_params$kaSyn_clear  <- 0.15
    baseline_sim <- pd_ode_euler(baseline_params, times) %>%
      mutate(time_years = time/365, group = "Baseline (Symptomatic Only)")

    bind_rows(intervention_sim, baseline_sim)
  })

  output$plot_np_snpc <- renderPlotly({
    df <- np_sim()
    p <- ggplot(df, aes(time_years, SNpc * 100, color = group)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("With Neuroprotection" = "#43A047",
                                    "Baseline (Symptomatic Only)" = "#B71C1C")) +
      labs(title = "Neuroprotection: SNpc Neuron Survival",
           x = "Time (years)", y = "Surviving Neurons (%)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_np_braak <- renderPlotly({
    df <- np_sim()
    p <- ggplot(df, aes(time_years, ASyn_O, color = group)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("With Neuroprotection" = "#43A047",
                                    "Baseline (Symptomatic Only)" = "#B71C1C")) +
      labs(title = "α-Syn Oligomers (Pathological Burden)",
           x = "Time (years)", y = "α-Syn Oligomers (nM)", color = "") +
      theme_bw()
    ggplotly(p)
  })

  output$table_np_endpoints <- renderDT({
    df <- np_sim()
    df %>%
      filter(round(time_years, 0) %in% c(5, 10, 15, 20, 25)) %>%
      group_by(group, Year = round(time_years, 0)) %>%
      summarise(
        SNpc_pct   = round(mean(SNpc) * 100, 1),
        UPDRS_III  = round(mean(UPDRS_III), 1),
        LID_risk   = round(mean(LID_risk), 4),
        .groups    = "drop"
      ) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })
}

## ─────────────────────────────────────────────────────────────────────────────
## Launch
## ─────────────────────────────────────────────────────────────────────────────

shinyApp(ui, server)
