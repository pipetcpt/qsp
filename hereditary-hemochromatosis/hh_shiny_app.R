## ============================================================
## Hereditary Hemochromatosis (HH) — Shiny QSP Dashboard
## 6 Tabs: Patient Profile | PK | Iron Kinetics | Organ Loading
##         Scenario Comparison | Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(plotly)
library(DT)

## ── Simulation function (lightweight ODE, no mrgsolve dependency) ──
simulate_HH <- function(genotype, init_ferritin, init_LIC, init_HEART_Fe,
                        phlebotomy_on, phlebotomy_freq_days,
                        DFO_on, DFO_dose,
                        DFX_on, DFX_dose,
                        DFP_on, DFP_dose,
                        HEPC_agonist,
                        sim_years, body_wt=70) {

  dt    <- 7      # step size (days)
  Nstep <- round(sim_years * 365 / dt)
  Fe_intake <- 15 # mg/day

  # Hepcidin based on genotype
  hep_base <- switch(as.character(genotype),
    "0" = 25, "1" = 5, "2" = 15, 25)
  if (HEPC_agonist) hep_base <- hep_base * 5

  # Initial states
  LIVER_Fe  <- init_LIC * 100   # mg (LIC mg/g × 100 g liver)
  FERRITIN  <- init_ferritin
  HEART_Fe  <- init_HEART_Fe
  PANCR_Fe  <- init_HEART_Fe * 0.3
  TBI       <- 5
  NTBI      <- 0.1
  HEPC      <- hep_base
  FPN       <- 1.0
  LIV_FIB   <- pmax(0, (init_LIC - 5) * 0.15)
  BCELL     <- pmax(0.3, 1 - PANCR_Fe * 0.05)
  HBA1C     <- 3.5 + 0.85 * (5.5 + 4 * (1-BCELL))
  T2STAR    <- 35 * exp(-0.12 * HEART_Fe)
  EF        <- 65

  out <- data.frame(
    time=numeric(Nstep), FERRITIN=numeric(Nstep), LIC=numeric(Nstep),
    T2STAR=numeric(Nstep), EF=numeric(Nstep), BCELL=numeric(Nstep),
    HBA1C=numeric(Nstep), LIV_FIB=numeric(Nstep), HEPC=numeric(Nstep),
    TSAT=numeric(Nstep), LIVER_Fe=numeric(Nstep), HEART_Fe=numeric(Nstep),
    NTBI=numeric(Nstep), TBI=numeric(Nstep)
  )

  for (i in seq_len(Nstep)) {
    t_now <- (i-1) * dt

    # Hepcidin dynamics (target-based, slow equilibration)
    liver_norm <- LIVER_Fe / 400
    HEPC <- HEPC + dt * 0.1 * (hep_base * liver_norm - HEPC)
    HEPC <- pmax(1, HEPC)

    # FPN activity (hepcidin suppresses FPN)
    FPN <- FPN + dt * (0.3 * (1 - FPN) - 0.3 * HEPC / (10 + HEPC) * FPN)
    FPN <- pmax(0.01, FPN)

    # Absorption
    abs_frac <- 0.10 + 0.25 * (1 - HEPC^2 / (10^2 + HEPC^2)) * FPN
    abs_frac <- pmin(0.40, pmax(0.05, abs_frac))
    Fe_abs   <- Fe_intake * abs_frac * dt

    # TBI & NTBI
    TSAT <- pmin(100, TBI / 10 * 100)
    NTBI_form <- ifelse(TSAT > 75, 0.8 * (TBI - 7.5) * dt, 0)
    TBI  <- TBI  + Fe_abs - 0.3 * TBI * dt - NTBI_form
    NTBI <- NTBI + NTBI_form - 0.19 * NTBI * dt

    # Phlebotomy: remove 250 mg Fe every N days
    phlebotomy_Fe <- 0
    if (phlebotomy_on && (t_now %% phlebotomy_freq_days) < dt) {
      phlebotomy_Fe <- 250
      LIVER_Fe <- pmax(0, LIVER_Fe - phlebotomy_Fe * 0.6)
      TBI      <- pmax(0, TBI - phlebotomy_Fe * 0.001)
    }

    # Chelation iron removal
    DFO_Fe_remove <- ifelse(DFO_on, DFO_dose * body_wt * 0.002 * dt, 0)
    DFX_Fe_remove <- ifelse(DFX_on, DFX_dose * body_wt * 0.0015 * dt, 0)
    DFP_Fe_remove <- ifelse(DFP_on, DFP_dose * body_wt * 0.001 * dt, 0)

    # Liver iron
    LIVER_Fe <- LIVER_Fe + (0.08 * TBI + 0.12 * NTBI) * dt
    LIVER_Fe <- LIVER_Fe - DFO_Fe_remove * 0.6 - DFX_Fe_remove * 0.7
    LIVER_Fe <- pmax(0, LIVER_Fe)

    # Cardiac iron
    HEART_Fe <- HEART_Fe + 0.04 * NTBI * dt
    HEART_Fe <- HEART_Fe - DFP_Fe_remove * 0.6 - DFO_Fe_remove * 0.1
    HEART_Fe <- pmax(0, HEART_Fe)

    # Pancreatic iron
    PANCR_Fe <- PANCR_Fe + 0.03 * NTBI * dt
    PANCR_Fe <- pmax(0, PANCR_Fe)

    # Ferritin
    ferritin_tgt <- LIVER_Fe / 8
    FERRITIN <- FERRITIN + 0.05 * (ferritin_tgt - FERRITIN) * dt

    # Fibrosis (Metavir)
    LIC_now <- LIVER_Fe / 100
    if (LIC_now > 7) {
      LIV_FIB <- LIV_FIB + 0.002 * (LIC_now - 7) * dt
    } else {
      LIV_FIB <- LIV_FIB - 0.001 * LIV_FIB * dt
    }
    LIV_FIB <- pmax(0, pmin(4, LIV_FIB))

    # β-cell
    BCELL <- BCELL - 0.01 * PANCR_Fe * BCELL * dt + 0.001 * (1 - BCELL) * dt
    BCELL <- pmax(0.1, pmin(1, BCELL))
    glu   <- 5.5 + 4 * (1 - BCELL)
    HBA1C <- HBA1C + 0.02 * (3.5 + 0.85 * glu - HBA1C) * dt

    # T2* & EF
    T2STAR <- pmax(3, 35 * exp(-0.12 * HEART_Fe))
    EF_drop <- (65 - 20) * HEART_Fe^2 / (5^2 + HEART_Fe^2)
    EF <- EF + 0.002 * (65 - EF_drop - EF) * dt
    EF <- pmax(15, pmin(75, EF))

    out[i,] <- list(t_now, FERRITIN, LIVER_Fe/100, T2STAR, EF, BCELL,
                    HBA1C, LIV_FIB, HEPC, TSAT, LIVER_Fe, HEART_Fe, NTBI, TBI)
  }
  out
}

## ── UI ────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "HH QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",   tabName="patient",  icon=icon("user")),
      menuItem("② Drug PK",           tabName="pk",       icon=icon("pills")),
      menuItem("③ Iron Kinetics",     tabName="iron",     icon=icon("tint")),
      menuItem("④ Organ Iron Load",   tabName="organ",    icon=icon("heart")),
      menuItem("⑤ Scenario Compare",  tabName="scenario", icon=icon("chart-line")),
      menuItem("⑥ Biomarkers",        tabName="biomark",  icon=icon("vial"))
    ),
    hr(),
    tags$div(style="padding:10px; color:white; font-size:11px;",
      "HH QSP Model v1.0", br(),
      "HFE C282Y/H63D", br(),
      "mrgsolve / ODE system"
    )
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────
      tabItem(tabName="patient",
        fluidRow(
          box(title="Patient Genetics & Phenotype", width=4, solidHeader=TRUE, status="success",
            selectInput("genotype", "HFE Genotype:",
              choices=list("Wild Type (WT)"=0,
                           "C282Y Homozygous"=1,
                           "C282Y/H63D Compound Het"=2),
              selected=1),
            numericInput("age",       "Age (years):",       value=45, min=18, max=80),
            selectInput("sex",        "Sex:",
              choices=list("Male"=1, "Female"=0), selected=1),
            numericInput("body_wt",   "Body Weight (kg):",  value=70, min=40, max=150),
            numericInput("sim_years", "Simulation (years):", value=10, min=1, max=20)
          ),
          box(title="Iron Status at Baseline", width=4, solidHeader=TRUE, status="warning",
            numericInput("init_ferritin", "Ferritin (ng/mL):",    value=800, min=10, max=5000),
            numericInput("init_LIC",      "LIC (mg/g dry wt):",   value=8,   min=0.5, max=50),
            numericInput("init_HEART_Fe", "Cardiac Iron (mg):",   value=8,   min=0, max=50),
            hr(),
            tags$p("Normal ranges:", style="font-weight:bold"),
            tags$ul(
              tags$li("Ferritin: <300 ng/mL (M), <200 (F)"),
              tags$li("LIC: <3 mg/g dry wt (target)"),
              tags$li("Cardiac T2*: >20 ms (normal)")
            )
          ),
          box(title="Disease Summary", width=4, solidHeader=TRUE, status="info",
            tags$h4("Hereditary Hemochromatosis"),
            tags$p("Autosomal recessive iron overload disorder due to HFE mutations."),
            tags$hr(),
            tags$table(class="table table-condensed",
              tags$tr(tags$th("Feature"), tags$th("HH")),
              tags$tr(tags$td("Prevalence"), tags$td("1/200–400 (Northern Europeans)")),
              tags$tr(tags$td("Gene"), tags$td("HFE (chr 6p21)")),
              tags$tr(tags$td("Key mutation"), tags$td("C282Y (~85%)")),
              tags$tr(tags$td("Penetrance"), tags$td("15–50% biochemical")),
              tags$tr(tags$td("Iron overload"), tags$td("20–40 g (vs. 3–5 g normal)")),
              tags$tr(tags$td("Diagnosis"), tags$td("TS% > 45%, ferritin elevated")),
              tags$tr(tags$td("Treatment"), tags$td("Phlebotomy (gold standard)"))
            )
          )
        ),
        fluidRow(
          box(title="Current Biomarker Status", width=12,
            DTOutput("baseline_table")
          )
        )
      ),

      ## ── Tab 2: Drug PK ─────────────────────────────────────
      tabItem(tabName="pk",
        fluidRow(
          box(title="Treatment Selection", width=3, solidHeader=TRUE, status="primary",
            checkboxInput("phlebotomy_on", "Phlebotomy", value=TRUE),
            conditionalPanel("input.phlebotomy_on",
              numericInput("phlebotomy_freq", "Frequency (every N days):", value=14, min=7, max=180)
            ),
            hr(),
            checkboxInput("DFO_on", "Deferoxamine (DFO)", value=FALSE),
            conditionalPanel("input.DFO_on",
              numericInput("DFO_dose", "DFO Dose (mg/kg/day):", value=40, min=10, max=60)
            ),
            hr(),
            checkboxInput("DFX_on", "Deferasirox (DFX)", value=FALSE),
            conditionalPanel("input.DFX_on",
              numericInput("DFX_dose", "DFX Dose (mg/kg/day):", value=20, min=5, max=40)
            ),
            hr(),
            checkboxInput("DFP_on", "Deferiprone (DFP)", value=FALSE),
            conditionalPanel("input.DFP_on",
              numericInput("DFP_dose", "DFP Dose (mg/kg/day TID):", value=75, min=25, max=100)
            ),
            hr(),
            checkboxInput("HEPC_agonist", "Hepcidin Agonist (exp.)", value=FALSE),
            actionButton("run_sim", "Run Simulation", class="btn-success btn-block")
          ),
          box(title="PK Profile — Oral Iron Chelators", width=9, solidHeader=TRUE,
            plotlyOutput("pk_plot", height="450px"),
            tags$p("Simulated plasma concentration profiles (normalized) over 24h single dose.",
                   style="color:gray; font-size:11px")
          )
        ),
        fluidRow(
          box(title="PK Parameters Summary", width=12,
            DTOutput("pk_table")
          )
        )
      ),

      ## ── Tab 3: Iron Kinetics ────────────────────────────────
      tabItem(tabName="iron",
        fluidRow(
          box(title="Iron Pool Dynamics", width=12, solidHeader=TRUE, status="success",
            plotlyOutput("iron_kinetics", height="550px")
          )
        ),
        fluidRow(
          box(title="Ferritin Trajectory", width=6,
            plotlyOutput("ferritin_plot", height="350px")
          ),
          box(title="Transferrin Saturation (TSAT)", width=6,
            plotlyOutput("tsat_plot", height="350px")
          )
        )
      ),

      ## ── Tab 4: Organ Iron Loading ───────────────────────────
      tabItem(tabName="organ",
        fluidRow(
          box(title="Liver Iron Concentration (LIC) & Fibrosis", width=6, solidHeader=TRUE, status="danger",
            plotlyOutput("lic_fib_plot", height="400px")
          ),
          box(title="Cardiac Iron — T2* MRI & Ejection Fraction", width=6, solidHeader=TRUE, status="warning",
            plotlyOutput("cardiac_plot", height="400px")
          )
        ),
        fluidRow(
          box(title="Pancreatic β-cell Function & HbA1c", width=6, solidHeader=TRUE, status="info",
            plotlyOutput("pancr_plot", height="350px")
          ),
          box(title="Hepcidin & FPN Dynamics", width=6, solidHeader=TRUE,
            plotlyOutput("hepc_plot", height="350px")
          )
        )
      ),

      ## ── Tab 5: Scenario Comparison ─────────────────────────
      tabItem(tabName="scenario",
        fluidRow(
          box(title="All Scenarios — Serum Ferritin", width=6,
            plotlyOutput("sc_ferr", height="350px")
          ),
          box(title="All Scenarios — LIC", width=6,
            plotlyOutput("sc_lic", height="350px")
          )
        ),
        fluidRow(
          box(title="All Scenarios — Cardiac T2*", width=6,
            plotlyOutput("sc_t2star", height="350px")
          ),
          box(title="All Scenarios — Liver Fibrosis", width=6,
            plotlyOutput("sc_fib", height="350px")
          )
        ),
        fluidRow(
          box(title="Scenario Parameter Summary", width=12,
            DTOutput("scenario_table")
          )
        )
      ),

      ## ── Tab 6: Biomarkers ──────────────────────────────────
      tabItem(tabName="biomark",
        fluidRow(
          box(title="Biomarker Dashboard (Year-End Values)", width=12, solidHeader=TRUE, status="success",
            DTOutput("biomarker_table")
          )
        ),
        fluidRow(
          box(title="HbA1c Trajectory (Bronze Diabetes)", width=6,
            plotlyOutput("hba1c_plot", height="350px")
          ),
          box(title="Iron Balance — Absorbed vs. Removed", width=6,
            plotlyOutput("fe_balance", height="350px")
          )
        ),
        fluidRow(
          box(title="Risk Stratification", width=12,
            plotlyOutput("risk_radar", height="400px")
          )
        )
      )
    )
  )
)

## ── SERVER ──────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation
  sim_result <- eventReactive(input$run_sim, {
    simulate_HH(
      genotype       = as.integer(input$genotype),
      init_ferritin  = input$init_ferritin,
      init_LIC       = input$init_LIC,
      init_HEART_Fe  = input$init_HEART_Fe,
      phlebotomy_on  = input$phlebotomy_on,
      phlebotomy_freq_days = input$phlebotomy_freq,
      DFO_on         = input$DFO_on,
      DFO_dose       = input$DFO_dose,
      DFX_on         = input$DFX_on,
      DFX_dose       = input$DFX_dose,
      DFP_on         = input$DFP_on,
      DFP_dose       = input$DFP_dose,
      HEPC_agonist   = input$HEPC_agonist,
      sim_years      = input$sim_years,
      body_wt        = input$body_wt
    )
  }, ignoreNULL=FALSE)

  ## Scenario data (all 6 pre-computed)
  scenario_data <- reactive({
    scens <- list(
      list(name="S1: Untreated",        phl=FALSE,14, DFO=FALSE,0, DFX=FALSE,0, DFP=FALSE,0, HEPC=FALSE,
           ferr=1500, LIC=15, heart=12),
      list(name="S2: Phlebotomy",       phl=TRUE, 14, DFO=FALSE,0, DFX=FALSE,0, DFP=FALSE,0, HEPC=FALSE,
           ferr=1500, LIC=15, heart=8),
      list(name="S3: DFO 40 mg/kg",     phl=FALSE,14, DFO=TRUE, 40,DFX=FALSE,0, DFP=FALSE,0, HEPC=FALSE,
           ferr=1500, LIC=15, heart=12),
      list(name="S4: DFX 20 mg/kg",     phl=FALSE,14, DFO=FALSE,0, DFX=TRUE, 20,DFP=FALSE,0, HEPC=FALSE,
           ferr=1500, LIC=15, heart=8),
      list(name="S5: DFP 75 mg/kg",     phl=FALSE,14, DFO=FALSE,0, DFX=FALSE,0, DFP=TRUE, 75,HEPC=FALSE,
           ferr=1200, LIC=12, heart=20),
      list(name="S6: DFP+DFO Combo",    phl=FALSE,14, DFO=TRUE, 30,DFX=FALSE,0, DFP=TRUE, 75,HEPC=FALSE,
           ferr=2000, LIC=20, heart=25)
    )
    do.call(rbind, lapply(seq_along(scens), function(j) {
      sc <- scens[[j]]
      df <- simulate_HH(
        genotype=1, init_ferritin=sc$ferr, init_LIC=sc$LIC, init_HEART_Fe=sc$heart,
        phlebotomy_on=sc$phl, phlebotomy_freq_days=14,
        DFO_on=sc$DFO, DFO_dose=sc[[5]],
        DFX_on=sc$DFX, DFX_dose=sc[[7]],
        DFP_on=sc$DFP, DFP_dose=sc[[9]],
        HEPC_agonist=FALSE, sim_years=5, body_wt=70
      )
      df$scenario <- sc$name
      df
    }))
  })

  ## ── Tab 1 outputs ────────────────────────────────────────────────
  output$baseline_table <- renderDT({
    df <- data.frame(
      Biomarker=c("Serum Ferritin (ng/mL)","LIC (mg/g)","Cardiac T2* (ms)",
                  "TSAT (%)","HbA1c (%)","HFE Genotype"),
      Value=c(input$init_ferritin, input$init_LIC,
              round(35*exp(-0.12*input$init_HEART_Fe),1),
              round(pmin(100, 5/10*100+input$init_LIC*5),1), 5.6,
              c("Wild Type","C282Y Homozygous","C282Y/H63D")[as.integer(input$genotype)+1]),
      Status=c(
        ifelse(input$init_ferritin > 300, "HIGH", "Normal"),
        ifelse(input$init_LIC > 3, "HIGH", "Normal"),
        ifelse(35*exp(-0.12*input$init_HEART_Fe) < 20, "LOW", "Normal"),
        ifelse(input$init_LIC*5+25 > 45, "HIGH", "Normal"),
        "Normal", "Genetic Risk"
      )
    )
    datatable(df, rownames=FALSE, options=list(pageLength=10, dom='t')) %>%
      formatStyle("Status",
        backgroundColor=styleEqual(c("HIGH","LOW","Normal","Genetic Risk"),
                                   c("#FFCDD2","#FFF9C4","#E8F5E9","#FFF3E0")))
  })

  ## ── Tab 2 outputs ────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    t_h <- seq(0, 24, by=0.5)
    pk_df <- bind_rows(
      data.frame(time=t_h,
        conc=30*(1-exp(-0.8*t_h))*exp(-0.3*t_h),  # DFO SC (simplified)
        drug="Deferoxamine (DFO) — SC", ka="0.8/h", t12="3-4h"),
      data.frame(time=t_h,
        conc=20*(1-exp(-0.8*t_h))*exp(-0.06*t_h), # DFX oral
        drug="Deferasirox (DFX) — oral", ka="0.8/h", t12="8-16h"),
      data.frame(time=t_h,
        conc=25*(1-exp(-18*t_h/24))*exp(-2.4*t_h/24), # DFP oral (short t½)
        drug="Deferiprone (DFP) — oral TID", ka="18/d", t12="2-3h")
    )
    plot_ly(pk_df, x=~time, y=~conc, color=~drug,
            type="scatter", mode="lines",
            line=list(width=2)) %>%
      layout(title="PK Profiles (Normalized Concentration)",
             xaxis=list(title="Time (h)"),
             yaxis=list(title="Normalized Concentration (AU)"),
             legend=list(orientation="h", y=-0.2))
  })

  output$pk_table <- renderDT({
    df <- data.frame(
      Drug=c("Deferoxamine (DFO)","Deferasirox (DFX)","Deferiprone (DFP)"),
      Route=c("SC/IV infusion","Oral once-daily","Oral TID"),
      Dose=c("40–50 mg/kg/day","20–40 mg/kg/day","75–100 mg/kg/day"),
      F_pct=c("N/A (parenteral)","~70%","~70%"),
      Tmax=c("~1–2h (SC)","~1.5–4h","~45 min"),
      t_half=c("~3–4h","~8–16h","~2–3h"),
      Vd=c("~0.8 L/kg","~14 L/kg","~1.1 L/kg"),
      Elimination=c("Renal + Fecal","Fecal (~84%)","Renal (glucuronide)"),
      Main_Target=c("Liver/NTBI (LIP)","Liver/NTBI","Cardiac iron"),
      stringsAsFactors=FALSE
    )
    datatable(df, rownames=FALSE, options=list(dom='t', scrollX=TRUE))
  })

  ## ── Tab 3 outputs ────────────────────────────────────────────────
  output$iron_kinetics <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    df$yr <- df$time / 365

    p <- plot_ly() %>%
      add_lines(data=df, x=~yr, y=~LIVER_Fe, name="Liver Fe (mg)", line=list(color="red")) %>%
      add_lines(data=df, x=~yr, y=~HEART_Fe*20, name="Cardiac Fe ×20 (mg)", line=list(color="purple", dash="dash")) %>%
      add_lines(data=df, x=~yr, y=~TBI*50, name="TBI ×50 (mg)", line=list(color="blue")) %>%
      add_lines(data=df, x=~yr, y=~NTBI*200, name="NTBI ×200 (mg)", line=list(color="orange", dash="dot")) %>%
      layout(title="Iron Pool Kinetics Over Time",
             xaxis=list(title="Time (years)"),
             yaxis=list(title="Iron (mg, scaled)"),
             legend=list(orientation="h", y=-0.2))
    p
  })

  output$ferritin_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly(df, x=~time/365, y=~FERRITIN, type="scatter", mode="lines",
            line=list(color="#FF5722", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(300,300),
                line=list(dash="dash", color="red"), name="Alarm (300)") %>%
      add_lines(x=range(df$time/365), y=c(50,50),
                line=list(dash="dash", color="green"), name="Target (50)") %>%
      layout(title="Serum Ferritin", xaxis=list(title="Years"),
             yaxis=list(title="Ferritin (ng/mL)"))
  })

  output$tsat_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly(df, x=~time/365, y=~TSAT, type="scatter", mode="lines",
            line=list(color="#2196F3", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(45,45),
                line=list(dash="dash", color="red"), name="Alarm (45%)") %>%
      layout(title="Transferrin Saturation (TSAT%)",
             xaxis=list(title="Years"), yaxis=list(title="TSAT (%)"))
  })

  ## ── Tab 4 outputs ────────────────────────────────────────────────
  output$lic_fib_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data=df, x=~time/365, y=~LIC, name="LIC (mg/g)",
                line=list(color="darkred", width=2)) %>%
      add_lines(data=df, x=~time/365, y=~LIV_FIB*5, name="Fibrosis×5 (Metavir)",
                line=list(color="brown", dash="dash", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(7,7),
                line=list(dash="dot", color="red"), name="LIC threshold") %>%
      layout(title="LIC & Liver Fibrosis",
             xaxis=list(title="Years"), yaxis=list(title="LIC (mg/g) / Fibrosis×5"))
  })

  output$cardiac_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data=df, x=~time/365, y=~T2STAR, name="T2* (ms)",
                line=list(color="purple", width=2)) %>%
      add_lines(data=df, x=~time/365, y=~EF, name="EF (%)",
                line=list(color="blue", dash="dash", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(20,20),
                line=list(dash="dot", color="red"), name="T2* threshold") %>%
      layout(title="Cardiac T2* & Ejection Fraction",
             xaxis=list(title="Years"), yaxis=list(title="T2* (ms) / EF (%)"))
  })

  output$pancr_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly() %>%
      add_lines(data=df, x=~time/365, y=~BCELL*100, name="β-cell Function (%)",
                line=list(color="orange", width=2)) %>%
      add_lines(data=df, x=~time/365, y=~HBA1C*10, name="HbA1c×10 (%)",
                line=list(color="red", dash="dash", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(65,65),
                line=list(dash="dot", color="red"), name="DM threshold HbA1c 6.5%") %>%
      layout(title="Pancreatic β-cell & HbA1c",
             xaxis=list(title="Years"), yaxis=list(title="β-cell (%) / HbA1c×10"))
  })

  output$hepc_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly(df, x=~time/365, y=~HEPC, type="scatter", mode="lines",
            line=list(color="darkgreen", width=2)) %>%
      add_lines(x=range(df$time/365), y=c(25,25),
                line=list(dash="dash", color="green"), name="Normal hepcidin") %>%
      layout(title="Plasma Hepcidin Dynamics",
             xaxis=list(title="Years"), yaxis=list(title="Hepcidin (ng/mL)"))
  })

  ## ── Tab 5 outputs ────────────────────────────────────────────────
  output$sc_ferr <- renderPlotly({
    df <- scenario_data()
    plot_ly(df, x=~time/365, y=~FERRITIN, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(title="Ferritin — All Scenarios",
             xaxis=list(title="Years"), yaxis=list(title="Ferritin (ng/mL)"),
             legend=list(orientation="h", y=-0.3))
  })

  output$sc_lic <- renderPlotly({
    df <- scenario_data()
    plot_ly(df, x=~time/365, y=~LIC, color=~scenario,
            type="scatter", mode="lines") %>%
      add_lines(x=c(0,5), y=c(3,3), line=list(dash="dash", color="green"),
                name="Target", inherit=FALSE) %>%
      layout(title="LIC — All Scenarios",
             xaxis=list(title="Years"), yaxis=list(title="LIC (mg/g)"),
             legend=list(orientation="h", y=-0.3))
  })

  output$sc_t2star <- renderPlotly({
    df <- scenario_data()
    plot_ly(df, x=~time/365, y=~T2STAR, color=~scenario,
            type="scatter", mode="lines") %>%
      add_lines(x=c(0,5), y=c(20,20), line=list(dash="dash", color="red"),
                name="Threshold", inherit=FALSE) %>%
      layout(title="Cardiac T2* — All Scenarios",
             xaxis=list(title="Years"), yaxis=list(title="T2* (ms)"),
             legend=list(orientation="h", y=-0.3))
  })

  output$sc_fib <- renderPlotly({
    df <- scenario_data()
    plot_ly(df, x=~time/365, y=~LIV_FIB, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(title="Liver Fibrosis — All Scenarios",
             xaxis=list(title="Years"), yaxis=list(title="Metavir Score"),
             legend=list(orientation="h", y=-0.3))
  })

  output$scenario_table <- renderDT({
    df <- data.frame(
      Scenario=c("S1: Untreated","S2: Phlebotomy","S3: DFO","S4: DFX","S5: DFP","S6: DFP+DFO"),
      Treatment=c("None","500 mL q2wk","DFO 40 mg/kg/d SC","DFX 20 mg/kg/d oral",
                  "DFP 75 mg/kg/d TID","DFP 75 + DFO 30 mg/kg/d"),
      Route=c("—","Venesection","SC infusion","Oral once-daily","Oral TID","Oral+SC"),
      `Iron Removed (mg/day)`=c(0, "~18 (avg)", "~70", "~30", "~25","~90"),
      `Cardiac Focus`=c("No","Limited","Moderate","Moderate","Yes","Yes"),
      `Liver Focus`=c("—","Primary","Primary","Primary","Secondary","Primary"),
      `Key AE`=c("Disease progression","Anemia","Retinopathy, ototoxicity",
                 "Nephrotoxicity, GI","Agranulocytosis","Combined AE"),
      stringsAsFactors=FALSE
    )
    datatable(df, rownames=FALSE, options=list(dom='t', scrollX=TRUE))
  })

  ## ── Tab 6 outputs ────────────────────────────────────────────────
  output$biomarker_table <- renderDT({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    last <- tail(df, 1)
    bm <- data.frame(
      Biomarker=c("Serum Ferritin","LIC","Cardiac T2*","EF","β-cell Function",
                  "HbA1c","Liver Fibrosis","Hepcidin","TSAT","Cardiac Fe"),
      Unit=c("ng/mL","mg/g","ms","%","%","%","Metavir","ng/mL","%","mg"),
      Value=round(c(last$FERRITIN, last$LIC, last$T2STAR, last$EF,
                    last$BCELL*100, last$HBA1C, last$LIV_FIB,
                    last$HEPC, last$TSAT, last$HEART_Fe), 2),
      Target=c("<50","<3",">20",">55",">80","<5.7","0","15–25","<45","<5"),
      Status=c(
        ifelse(last$FERRITIN<50,"Target",ifelse(last$FERRITIN<300,"Elevated","High")),
        ifelse(last$LIC<3,"Normal",ifelse(last$LIC<7,"Elevated","High")),
        ifelse(last$T2STAR>20,"Normal","Overloaded"),
        ifelse(last$EF>55,"Normal",ifelse(last$EF>40,"Reduced","Severely Reduced")),
        ifelse(last$BCELL>0.8,"Normal","Impaired"),
        ifelse(last$HBA1C<5.7,"Normal",ifelse(last$HBA1C<6.5,"Pre-DM","DM")),
        ifelse(last$LIV_FIB<1,"None",ifelse(last$LIV_FIB<3,"Fibrosis","Cirrhosis")),
        ifelse(last$HEPC>10,"Normal","Low"),
        ifelse(last$TSAT<45,"Normal","Elevated"),
        ifelse(last$HEART_Fe<5,"Normal","Elevated")
      )
    )
    datatable(bm, rownames=FALSE) %>%
      formatStyle("Status",
        backgroundColor=styleEqual(
          c("Target","Normal","Elevated","High","Overloaded","Impaired",
            "Pre-DM","DM","Fibrosis","Cirrhosis","Low","Severely Reduced","Reduced"),
          c("#E8F5E9","#E8F5E9","#FFF9C4","#FFCDD2","#FFCDD2","#FFF9C4",
            "#FFF9C4","#FFCDD2","#FFF9C4","#FFCDD2","#FFCDD2","#FFCDD2","#FFF9C4")
        )
      )
  })

  output$hba1c_plot <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    plot_ly(df, x=~time/365, y=~HBA1C, type="scatter", mode="lines",
            line=list(color="darkorange", width=2), name="HbA1c") %>%
      add_lines(x=range(df$time/365), y=c(5.7,5.7),
                line=list(dash="dash", color="orange"), name="Pre-DM (5.7%)") %>%
      add_lines(x=range(df$time/365), y=c(6.5,6.5),
                line=list(dash="dash", color="red"), name="DM (6.5%)") %>%
      layout(title="HbA1c — Bronze Diabetes Risk",
             xaxis=list(title="Years"), yaxis=list(title="HbA1c (%)"))
  })

  output$fe_balance <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    # Approximate balance: absorption minus removal
    phl_remove <- ifelse(input$phlebotomy_on, 250 / input$phlebotomy_freq * 7, 0)
    dfo_remove <- ifelse(input$DFO_on, input$DFO_dose * input$body_wt * 0.002 * 7, 0)
    dfx_remove <- ifelse(input$DFX_on, input$DFX_dose * input$body_wt * 0.0015 * 7, 0)
    dfp_remove <- ifelse(input$DFP_on, input$DFP_dose * input$body_wt * 0.001 * 7, 0)
    total_remove <- phl_remove + dfo_remove + dfx_remove + dfp_remove

    fe_abs_avg <- 15 * 0.25  # HH typical ~25% absorption → ~3.75 mg/day × 7 = 26.25/week
    balance_df <- data.frame(
      Component=c("Dietary Absorption","Phlebotomy","DFO","DFX","DFP"),
      Amount=c(fe_abs_avg*7, -phl_remove, -dfo_remove, -dfx_remove, -dfp_remove),
      Direction=c("In","Out","Out","Out","Out")
    )
    balance_df <- balance_df[balance_df$Amount != 0,]
    plot_ly(balance_df, x=~Component, y=~Amount, type="bar",
            color=~Direction, colors=c("Out"="#EF5350","In"="#4CAF50")) %>%
      layout(title="Weekly Iron Balance (mg/week)",
             xaxis=list(title="Component"),
             yaxis=list(title="Iron (mg/week)"))
  })

  output$risk_radar <- renderPlotly({
    df <- sim_result()
    if (is.null(df)) return(NULL)
    last <- tail(df, 1)

    # Normalize risks to 0–10 scale
    risks <- c(
      "Liver Risk"   = pmin(10, last$LIC * 0.5 + last$LIV_FIB * 1.5),
      "Cardiac Risk" = pmin(10, (35 - last$T2STAR) / 3.2),
      "DM Risk"      = pmin(10, (last$HBA1C - 5) * 4),
      "Iron Overload"= pmin(10, last$FERRITIN / 300),
      "Hepcidin Def" = pmin(10, (25 - last$HEPC) / 2),
      "Hypogonadism" = pmin(10, (last$LIVER_Fe / 2000) * 5)
    )

    plot_ly(
      type="scatterpolar", mode="lines+markers",
      r=c(risks, risks[1]), theta=c(names(risks), names(risks)[1]),
      fill="toself", fillcolor="rgba(255,87,34,0.2)",
      line=list(color="#FF5722")
    ) %>%
      layout(
        title="Multi-Organ Risk Radar",
        polar=list(radialaxis=list(range=c(0,10), angle=90)),
        showlegend=FALSE
      )
  })
}

shinyApp(ui, server)
