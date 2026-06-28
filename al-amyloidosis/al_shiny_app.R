## ============================================================
## AL Amyloidosis QSP Shiny Application
## Tabs: Patient Profile | PK Profiles | Hematologic PD |
##       Organ Biomarkers | Clinical Endpoints | Scenario Comparison |
##       Biomarker Staging | Safety Monitor
## ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(dplyr)

## ============================================================
## SIMULATION ENGINE (Pure R, no mrgsolve dependency)
## ============================================================

simulate_al <- function(
  scenario     = "Dara-CyBorD",
  bw           = 75,
  bsa          = 1.7,
  age          = 65,
  mayo_stage   = 3,
  FLC_init     = 150,
  NTproBNP_init = 3500,
  TnT_init     = 0.07,
  eGFR_init    = 35,
  n_days       = 365,
  dara_dose    = 16,
  btz_dose     = 1.3,
  cy_dose      = 300,
  dex_dose     = 20,
  cyp2c19_pm   = FALSE
) {

  dt    <- 1
  times <- seq(0, n_days, by = dt)
  N     <- length(times)

  ## --- Parameters ---
  CL_DARA <- 3.1; V1_DARA <- 56; kint_DARA <- 0.15
  kon_DARA <- 0.005; koff_DARA <- 0.002; Rss <- 120; kdeg_R <- 0.08
  ksyn_R   <- Rss * kdeg_R

  kprolif   <- 0.045; kdeath <- 0.045
  kmax_DARA <- 0.85; EC50_DARA <- 0.15
  kmax_BTZ  <- 0.75; EC50_BTZ  <- 6.5
  kmax_CY   <- 0.55; EC50_CY   <- 3.0
  kmax_DEX  <- 0.65; EC50_DEX  <- 8.0

  kFLC_prod  <- 0.2; kFLC_elim <- 0.001333
  k_dep_card <- 0.003; k_res_card <- 0.0005
  k_dep_ren  <- 0.002; k_res_ren  <- 0.0008
  kBNP_prod  <- 1.5;  kBNP_elim  <- 0.00043
  kTnT_prod  <- 0.002; kTnT_elim <- 0.006
  k_GFR_loss <- 0.0008; k_Prot <- 1.2
  k_NK_kill  <- 0.3; k_NK_rec <- 0.04

  kel_BTZ   <- ifelse(cyp2c19_pm, 0.018 * 0.65, 0.018)
  kon_BTZ   <- 1.5; koff_BTZ <- 0.006
  ka_CY     <- 1.2; CL_CY <- 5.8; V1_CY <- 38; F_CY <- 0.74
  ka_DEX    <- 1.5; CL_DEX <- 15.4; V1_DEX <- 64; F_DEX <- 0.80

  CL_BW_scaled <- CL_DARA * (bw / 70)^0.75

  ## --- Build dose schedule ---
  dose_dara <- rep(0, N)
  dose_btz  <- rep(0, N)
  dose_cy_gut  <- rep(0, N)
  dose_dex_gut <- rep(0, N)

  dara_amt_ugmL <- dara_dose * bw * 1000 / (bw * V1_DARA)

  if (scenario %in% c("Dara Mono", "Dara-CyBorD")) {
    qw_times   <- seq(0, 7*7, by=7)
    q2w_times  <- seq(56, 56 + 15*14, by=14)
    q4w_times  <- seq(56 + 16*14, n_days, by=28)
    dara_times <- unique(c(qw_times, q2w_times, q4w_times))
    dara_times <- dara_times[dara_times <= n_days] + 1
    dara_times <- dara_times[dara_times <= N]
    dose_dara[dara_times] <- dara_amt_ugmL
  }

  btz_amt  <- btz_dose * bsa * 1e6 / (635.84 * 498)
  cy_amt   <- cy_dose  * bsa
  dex_amt  <- dex_dose * 1e6 / (392.46 * V1_DEX)
  cyborg_days <- if (scenario %in% c("CyBorD", "Dara-CyBorD", "VCD")) {
    apply(expand.grid(seq(0, 5*28, by=28), c(1,8,15,22)), 1, sum) + 1
  } else c()

  if (length(cyborg_days) > 0) {
    valid_idx <- cyborg_days[cyborg_days <= N]
    dose_btz[valid_idx]     <- btz_amt
    dose_cy_gut[valid_idx]  <- cy_amt
    dex_dose_use <- ifelse(scenario == "VCD", 40, dex_dose)
    dex_amt_use  <- dex_dose_use * 1e6 / (392.46 * V1_DEX)
    dose_dex_gut[valid_idx] <- dex_amt_use
  }

  ## --- State vectors ---
  DARA_C1 <- DARA_C2 <- RC <- numeric(N)
  RFree    <- numeric(N); RFree[1] <- Rss
  BTZ_C1 <- BTZ_PROT <- numeric(N)
  CY_GUT <- CY_C1 <- numeric(N)
  DEX_GUT <- DEX_C1 <- numeric(N)
  PC      <- numeric(N); PC[1] <- 1.0
  FLC     <- numeric(N); FLC[1] <- FLC_init
  AmylCard <- numeric(N); AmylCard[1] <- 1.0
  AmylRen  <- numeric(N); AmylRen[1]  <- 1.0
  NTpBNP   <- numeric(N); NTpBNP[1]   <- NTproBNP_init
  TnT      <- numeric(N); TnT[1]      <- TnT_init
  GFR      <- numeric(N); GFR[1]      <- eGFR_init
  NK       <- numeric(N); NK[1]       <- 1.0

  ## --- Euler integration ---
  for (i in seq_len(N - 1)) {
    # Daratumumab
    DARA_C1[i+1] <- max(0, DARA_C1[i] + dt * (
      dose_dara[i]
      - (CL_BW_scaled / V1_DARA) * DARA_C1[i]
      - (3.0 / V1_DARA) * (DARA_C1[i] - DARA_C2[i])
      - kon_DARA * DARA_C1[i] * RFree[i]
      + koff_DARA * RC[i]
    ))
    DARA_C2[i+1] <- max(0, DARA_C2[i] + dt * (
      (3.0 / V1_DARA) * (DARA_C1[i] - DARA_C2[i])
      - (3.0 / 40) * DARA_C2[i]
    ))
    RC[i+1]    <- max(0, RC[i] + dt * (
      kon_DARA * DARA_C1[i] * RFree[i] - koff_DARA * RC[i] - kint_DARA * RC[i]
    ))
    RFree[i+1] <- max(0, RFree[i] + dt * (
      ksyn_R - kdeg_R * RFree[i]
      - kon_DARA * DARA_C1[i] * RFree[i]
      + koff_DARA * RC[i]
    ))

    # Bortezomib
    BTZ_bind   <- kon_BTZ * BTZ_C1[i] * pmax(0, 100 - BTZ_PROT[i])
    BTZ_unbind <- koff_BTZ * BTZ_PROT[i]
    BTZ_C1[i+1]   <- max(0, BTZ_C1[i] + dt * (dose_btz[i] - kel_BTZ * BTZ_C1[i] - BTZ_bind + BTZ_unbind))
    BTZ_PROT[i+1] <- max(0, BTZ_PROT[i] + dt * (BTZ_bind - BTZ_unbind))

    # CY
    CY_GUT[i+1] <- max(0, CY_GUT[i] + dt * (dose_cy_gut[i] - ka_CY * CY_GUT[i]))
    CY_C1[i+1]  <- max(0, CY_C1[i]  + dt * (ka_CY * CY_GUT[i] * F_CY / V1_CY - (CL_CY/V1_CY) * CY_C1[i]))

    # DEX
    DEX_GUT[i+1] <- max(0, DEX_GUT[i] + dt * (dose_dex_gut[i] - ka_DEX * DEX_GUT[i]))
    DEX_C1[i+1]  <- max(0, DEX_C1[i]  + dt * (ka_DEX * DEX_GUT[i] * F_DEX / V1_DEX - (CL_DEX/V1_DEX) * DEX_C1[i]))

    # Drug effects
    E_DARA <- kmax_DARA * DARA_C1[i] / (EC50_DARA + DARA_C1[i])
    E_BTZ  <- kmax_BTZ  * BTZ_PROT[i] / (EC50_BTZ  + BTZ_PROT[i])
    E_CY   <- kmax_CY   * CY_C1[i]   / (EC50_CY   + CY_C1[i])
    E_DEX  <- kmax_DEX  * DEX_C1[i]  / (EC50_DEX  + DEX_C1[i])
    NK_eff <- NK[i]
    E_kill <- pmin(0.98, E_DARA * NK_eff + E_BTZ + E_CY + E_DEX)

    # NK cells
    NK[i+1] <- max(0, NK[i] + dt * (-k_NK_kill * E_DARA * NK[i] + k_NK_rec * (1 - NK[i])))

    # Plasma cells
    PC[i+1] <- max(0, PC[i] + dt * (kprolif * PC[i] * (1 - E_kill) - kdeath * PC[i]))

    # FLC
    kFLC_adj <- kFLC_elim * (GFR[i] / eGFR_init)
    FLC[i+1] <- max(0, FLC[i] + dt * (kFLC_prod * PC[i] * FLC_init - kFLC_adj * FLC[i]))

    # Amyloid
    AmylCard[i+1] <- max(0, AmylCard[i] + dt * (k_dep_card * FLC[i] - k_res_card * AmylCard[i]))
    AmylRen[i+1]  <- max(0, AmylRen[i]  + dt * (k_dep_ren  * FLC[i] - k_res_ren  * AmylRen[i]))

    # Biomarkers
    NTpBNP[i+1] <- max(0, NTpBNP[i] + dt * (kBNP_prod * AmylCard[i] * NTproBNP_init - kBNP_elim * NTpBNP[i]))
    TnT[i+1]    <- max(0, TnT[i]    + dt * (kTnT_prod * AmylCard[i] - kTnT_elim * TnT[i]))
    GFR[i+1]    <- max(5,  GFR[i]   + dt * (-k_GFR_loss * AmylRen[i] * GFR[i]))
  }

  data.frame(
    time       = times,
    DARA_Cp    = DARA_C1,
    BTZ_bound  = BTZ_PROT,
    CY_Cp      = CY_C1,
    DEX_Cp     = DEX_C1,
    PC         = PC,
    dFLC       = FLC,
    AmylCard   = AmylCard,
    AmylRen    = AmylRen,
    NTproBNP   = NTpBNP,
    TnT        = TnT,
    eGFR       = GFR,
    NK         = NK,
    ProtUria   = k_Prot * AmylRen,
    CR_flag    = as.integer(FLC < 40),
    CardResp   = as.integer((NTproBNP_init - NTpBNP) / NTproBNP_init > 0.30),
    RenResp    = as.integer((k_Prot - k_Prot * AmylRen) / k_Prot > 0.30)
  )
}

ALL_SCENARIOS <- c(
  "Untreated",
  "Dara Mono",
  "CyBorD",
  "Dara-CyBorD",
  "VCD",
  "Dara-CyBorD CYP2C19 PM"
)

SCENARIO_COLORS <- c(
  "Untreated"             = "#e74c3c",
  "Dara Mono"             = "#3498db",
  "CyBorD"                = "#f39c12",
  "Dara-CyBorD"           = "#2ecc71",
  "VCD"                   = "#9b59b6",
  "Dara-CyBorD CYP2C19 PM" = "#1abc9c"
)

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "AL Amyloidosis QSP Dashboard",
    titleWidth = 310
  ),
  dashboardSidebar(
    width = 310,
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "patient",   icon = icon("user")),
      menuItem("PK Profiles",           tabName = "pk",        icon = icon("chart-line")),
      menuItem("Hematologic PD",        tabName = "heme",      icon = icon("vial")),
      menuItem("Organ Biomarkers",      tabName = "organs",    icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",    tabName = "endpoints", icon = icon("flag-checkered")),
      menuItem("Scenario Comparison",   tabName = "compare",   icon = icon("bar-chart")),
      menuItem("Biomarker Staging",     tabName = "staging",   icon = icon("layer-group")),
      menuItem("Safety Monitor",        tabName = "safety",    icon = icon("shield-alt"))
    ),
    hr(),
    h4("Treatment Scenario", style="color:white; padding-left:15px;"),
    selectInput("scenario", NULL,
                choices = ALL_SCENARIOS, selected = "Dara-CyBorD"),
    h4("Patient Parameters", style="color:white; padding-left:15px;"),
    sliderInput("bw",   "Body Weight (kg)",     50, 120,  75, step=5),
    sliderInput("bsa",  "BSA (m²)",             1.2, 2.5, 1.7, step=0.05),
    sliderInput("age",  "Age (years)",          30, 85,   65, step=5),
    sliderInput("mayo_stage", "Mayo 2012 Stage (1-4)", 1, 4, 3, step=1),
    hr(),
    h4("Baseline Biomarkers", style="color:white; padding-left:15px;"),
    sliderInput("FLC_init",    "Baseline dFLC (mg/L)",      50, 500, 150, step=10),
    sliderInput("NTproBNP_init","Baseline NT-proBNP (pg/mL)",500, 10000, 3500, step=100),
    sliderInput("TnT_init",    "Baseline hs-TnT (ng/mL)",  0.01, 0.30, 0.07, step=0.01),
    sliderInput("eGFR_init",   "Baseline eGFR (mL/min)",   10,   90,   35,   step=5),
    hr(),
    h4("Simulation", style="color:white; padding-left:15px;"),
    sliderInput("n_days", "Follow-up (days)", 90, 730, 365, step=30),
    checkboxInput("cyp2c19_pm", "CYP2C19 Poor Metabolizer (BTZ +35% exposure)", FALSE),
    actionButton("run_sim", "Run Simulation", icon=icon("play"),
                 class="btn-success btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f6f9; }
      .box-header { background: linear-gradient(135deg, #1a1a2e, #16213e); color: white; }
      .nav-tabs-custom>.tab-content { padding: 10px; }
    "))),
    tabItems(

      ## TAB 1: Patient Profile
      tabItem(tabName = "patient",
        fluidRow(
          box(title="Patient Overview", width=12, solidHeader=TRUE, status="primary",
            fluidRow(
              valueBoxOutput("vbox_mayo",    width=3),
              valueBoxOutput("vbox_dflc",    width=3),
              valueBoxOutput("vbox_bnp",     width=3),
              valueBoxOutput("vbox_egfr",    width=3)
            )
          )
        ),
        fluidRow(
          box(title="Mayo 2012 Staging System", width=6, solidHeader=TRUE, status="info",
            DTOutput("mayo_table")
          ),
          box(title="AL Amyloidosis Disease Overview", width=6, solidHeader=TRUE, status="info",
            p(strong("Pathophysiology:"), "Clonal plasma cells produce amyloidogenic
              monoclonal free light chains (FLC) that misfold into β-sheet amyloid fibrils,
              depositing in organs including heart, kidney, liver, and peripheral nerves."),
            p(strong("Epidemiology:"), "Incidence ~10 per million/year; median age ~63 years;
              λ-type light chains more common than κ in AL amyloidosis."),
            p(strong("Prognosis:"), "Highly dependent on cardiac involvement (Mayo Stage).
              Stage IV: median OS ~5-6 months without treatment."),
            p(strong("Treatment Goal:"), "Deep and rapid hematologic response (CR/VGPR)
              to halt FLC production and allow organ recovery."),
            p(strong("Key Biomarkers:"),
              tags$ul(
                tags$li("dFLC (difference FLC): primary PD marker; CR = dFLC <40 mg/L"),
                tags$li("NT-proBNP: cardiac amyloid burden; cardiac response = ≥30% reduction"),
                tags$li("hs-TnT: cardiomyocyte injury; Mayo 2012 risk factor if ≥0.025 ng/mL"),
                tags$li("eGFR + 24h proteinuria: renal response = ≥30% proteinuria reduction")
              )
            )
          )
        ),
        fluidRow(
          box(title="Drug Mechanism Summary", width=12, solidHeader=TRUE, status="warning",
            DTOutput("drug_table")
          )
        )
      ),

      ## TAB 2: PK Profiles
      tabItem(tabName = "pk",
        fluidRow(
          box(title="Daratumumab Plasma Concentration (TMDD PK)", width=6,
              solidHeader=TRUE, status="primary",
            plotlyOutput("plt_dara_pk", height=300)),
          box(title="Bortezomib-Proteasome Complex (Bound BTZ)", width=6,
              solidHeader=TRUE, status="primary",
            plotlyOutput("plt_btz_pk", height=300))
        ),
        fluidRow(
          box(title="Cyclophosphamide Active Metabolite", width=6,
              solidHeader=TRUE, status="info",
            plotlyOutput("plt_cy_pk", height=300)),
          box(title="Dexamethasone Plasma Concentration", width=6,
              solidHeader=TRUE, status="info",
            plotlyOutput("plt_dex_pk", height=300))
        ),
        fluidRow(
          box(title="PK Parameter Summary", width=12, solidHeader=TRUE, status="success",
            DTOutput("pk_table")
          )
        )
      ),

      ## TAB 3: Hematologic PD
      tabItem(tabName = "heme",
        fluidRow(
          box(title="Plasma Cell Pool (Normalized)", width=6,
              solidHeader=TRUE, status="danger",
            plotlyOutput("plt_pc", height=300)),
          box(title="dFLC Over Time (Hematologic Response)", width=6,
              solidHeader=TRUE, status="danger",
            plotlyOutput("plt_dflc", height=300))
        ),
        fluidRow(
          box(title="NK Cell Pool (Daratumumab ADCC)", width=6,
              solidHeader=TRUE, status="warning",
            plotlyOutput("plt_nk", height=300)),
          box(title="Hematologic Response Classification", width=6,
              solidHeader=TRUE, status="warning",
            DTOutput("heme_resp_table")
          )
        )
      ),

      ## TAB 4: Organ Biomarkers
      tabItem(tabName = "organs",
        fluidRow(
          box(title="NT-proBNP (Cardiac Biomarker)", width=6,
              solidHeader=TRUE, status="danger",
            plotlyOutput("plt_bnp", height=300)),
          box(title="hs-Troponin T (Cardiomyocyte Injury)", width=6,
              solidHeader=TRUE, status="danger",
            plotlyOutput("plt_tnt", height=300))
        ),
        fluidRow(
          box(title="eGFR Trajectory (Renal Function)", width=6,
              solidHeader=TRUE, status="info",
            plotlyOutput("plt_egfr", height=300)),
          box(title="24h Proteinuria", width=6,
              solidHeader=TRUE, status="info",
            plotlyOutput("plt_prot", height=300))
        )
      ),

      ## TAB 5: Clinical Endpoints
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title="Cardiac Amyloid Burden (Normalized)", width=6,
              solidHeader=TRUE, status="danger",
            plotlyOutput("plt_amyl_card", height=300)),
          box(title="Renal Amyloid Burden (Normalized)", width=6,
              solidHeader=TRUE, status="info",
            plotlyOutput("plt_amyl_ren", height=300))
        ),
        fluidRow(
          box(title="Organ Response Flags Over Time", width=6,
              solidHeader=TRUE, status="success",
            plotlyOutput("plt_org_resp", height=300)),
          box(title="Key Clinical Trial Benchmarks", width=6,
              solidHeader=TRUE, status="primary",
            DTOutput("trial_table")
          )
        )
      ),

      ## TAB 6: Scenario Comparison
      tabItem(tabName = "compare",
        fluidRow(
          box(title="Select Scenarios to Compare", width=12, solidHeader=TRUE, status="primary",
            checkboxGroupInput("compare_scenarios", NULL,
              choices = ALL_SCENARIOS,
              selected = c("Untreated","CyBorD","Dara-CyBorD"),
              inline = TRUE)
          )
        ),
        fluidRow(
          box(title="dFLC Comparison", width=6, solidHeader=TRUE, status="primary",
            plotlyOutput("cmp_dflc", height=280)),
          box(title="NT-proBNP Comparison", width=6, solidHeader=TRUE, status="primary",
            plotlyOutput("cmp_bnp", height=280))
        ),
        fluidRow(
          box(title="eGFR Comparison", width=6, solidHeader=TRUE, status="info",
            plotlyOutput("cmp_egfr", height=280)),
          box(title="Plasma Cell Pool Comparison", width=6, solidHeader=TRUE, status="info",
            plotlyOutput("cmp_pc", height=280))
        ),
        fluidRow(
          box(title="Day-180 Summary Table", width=12, solidHeader=TRUE, status="success",
            DTOutput("compare_table")
          )
        )
      ),

      ## TAB 7: Biomarker Staging
      tabItem(tabName = "staging",
        fluidRow(
          box(title="Mayo 2012 Stage Evolution Over Time", width=8,
              solidHeader=TRUE, status="warning",
            plotlyOutput("plt_stage", height=350)),
          box(title="Risk Factor Status", width=4,
              solidHeader=TRUE, status="warning",
            plotlyOutput("plt_rf_gauge", height=350))
        ),
        fluidRow(
          box(title="Biomarker Threshold Monitor", width=12,
              solidHeader=TRUE, status="info",
            DTOutput("staging_monitor_table")
          )
        )
      ),

      ## TAB 8: Safety Monitor
      tabItem(tabName = "safety",
        fluidRow(
          box(title="Key Safety Events", width=12, solidHeader=TRUE, status="danger",
            DTOutput("safety_table")
          )
        ),
        fluidRow(
          box(title="Daratumumab Infusion-Related Reaction Risk", width=6,
              solidHeader=TRUE, status="warning",
            plotlyOutput("plt_irr", height=280)),
          box(title="Immunoparesis (NK Depletion)", width=6,
              solidHeader=TRUE, status="warning",
            plotlyOutput("plt_nk_safety", height=280))
        ),
        fluidRow(
          box(title="REMS Monitoring Schedule & DDI Table", width=12,
              solidHeader=TRUE, status="primary",
            DTOutput("ddi_table")
          )
        )
      )

    ) ## end tabItems
  ) ## end dashboardBody
)

## ============================================================
## SERVER
## ============================================================

server <- function(input, output, session) {

  sim_result <- eventReactive(input$run_sim, {
    cyp_pm <- input$cyp2c19_pm || grepl("PM", input$scenario)
    simulate_al(
      scenario      = input$scenario,
      bw            = input$bw,
      bsa           = input$bsa,
      age           = input$age,
      mayo_stage    = input$mayo_stage,
      FLC_init      = input$FLC_init,
      NTproBNP_init = input$NTproBNP_init,
      TnT_init      = input$TnT_init,
      eGFR_init     = input$eGFR_init,
      n_days        = input$n_days,
      cyp2c19_pm    = cyp_pm
    )
  }, ignoreNULL = FALSE)

  compare_results <- reactive({
    req(input$compare_scenarios)
    scens <- input$compare_scenarios
    lapply(scens, function(s) {
      cyp_pm <- grepl("PM", s) || (s == input$scenario && input$cyp2c19_pm)
      df <- simulate_al(
        scenario      = s,
        bw            = input$bw,
        bsa           = input$bsa,
        FLC_init      = input$FLC_init,
        NTproBNP_init = input$NTproBNP_init,
        TnT_init      = input$TnT_init,
        eGFR_init     = input$eGFR_init,
        n_days        = input$n_days,
        cyp2c19_pm    = cyp_pm
      )
      df$Scenario <- s
      df
    }) |> bind_rows()
  })

  ## Value boxes
  output$vbox_mayo <- renderValueBox({
    s <- input$mayo_stage
    col <- c("green","yellow","orange","red")[s]
    valueBox(paste0("Stage ", s), "Mayo 2012 Stage", icon=icon("layer-group"), color=col)
  })
  output$vbox_dflc <- renderValueBox({
    valueBox(paste0(input$FLC_init, " mg/L"), "Baseline dFLC",
             icon=icon("vial"), color="light-blue")
  })
  output$vbox_bnp <- renderValueBox({
    valueBox(paste0(format(input$NTproBNP_init, big.mark=","), " pg/mL"), "Baseline NT-proBNP",
             icon=icon("heartbeat"), color="red")
  })
  output$vbox_egfr <- renderValueBox({
    ckd <- ifelse(input$eGFR_init >= 60, "CKD 1-2", ifelse(input$eGFR_init >= 30, "CKD 3", "CKD 4-5"))
    valueBox(paste0(input$eGFR_init, " mL/min (", ckd, ")"), "Baseline eGFR",
             icon=icon("kidneys"), color="olive")
  })

  ## Mayo staging table
  output$mayo_table <- renderDT({
    df <- data.frame(
      Stage = c("I","II","III","IV"),
      Risk_Factors = c("0","1","2","3 or 4"),
      Criteria = c(
        "hs-TnT <0.025 ng/mL AND NT-proBNP <1800 pg/mL AND dFLC <18 mg/dL AND eGFR ≥50",
        "Any 1 of: TnT ≥0.025 / BNP ≥1800 / dFLC ≥18 / eGFR <50",
        "Any 2 of the above criteria",
        "3 or 4 risk factors met"
      ),
      Median_OS = c(">10 years", "~6 years", "~4 years", "~5-12 months")
    )
    datatable(df, options=list(dom="t", pageLength=4), rownames=FALSE)
  })

  ## Drug mechanism table
  output$drug_table <- renderDT({
    df <- data.frame(
      Drug = c("Daratumumab","Bortezomib","Cyclophosphamide","Dexamethasone","Melphalan"),
      Class = c("Anti-CD38 mAb","Proteasome inhibitor","Alkylating agent","Glucocorticoid","Alkylating agent"),
      Mechanism = c(
        "Binds CD38 on plasma cells → ADCC (NK), CDC (complement), ADCP (macrophage)",
        "Inhibits 20S proteasome → UPR → PC apoptosis",
        "DNA alkylation → clonal PC death (active metabolite 4-OH-CY)",
        "Glucocorticoid receptor → pro-apoptotic gene activation in PCs",
        "DNA cross-linking → PC death; used in SCT conditioning (200 mg/m²)"
      ),
      Key_Dosing = c(
        "16 mg/kg IV: QW×8, Q2W×16, Q4W (ANDROMEDA schedule)",
        "1.3 mg/m² SC on D1,8,15,22 of 28-day cycle",
        "300 mg/m² PO on D1,8,15,22",
        "20-40 mg PO on D1,8,15,22",
        "0.22 mg/kg/day PO D1-4 (maintenance) or 200 mg/m² IV (SCT conditioning)"
      )
    )
    datatable(df, options=list(dom="t", pageLength=5, scrollX=TRUE), rownames=FALSE)
  })

  ## PK plots
  mk_line <- function(df, yvar, ytitle, hline=NULL, hline_label=NULL, color="#3498db") {
    p <- plot_ly(df, x=~time, y=~get(yvar), type="scatter", mode="lines",
                 line=list(color=color, width=2.5)) |>
      layout(xaxis=list(title="Time (days)"), yaxis=list(title=ytitle),
             plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")
    if (!is.null(hline)) {
      p <- p |> add_segments(x=0, xend=max(df$time), y=hline, yend=hline,
                              line=list(color="red", dash="dash", width=1.5),
                              name=hline_label, showlegend=FALSE)
    }
    p
  }

  output$plt_dara_pk <- renderPlotly({
    df <- sim_result()
    mk_line(df, "DARA_Cp", "Daratumumab (µg/mL)", color="#2980b9")
  })
  output$plt_btz_pk <- renderPlotly({
    df <- sim_result()
    mk_line(df, "BTZ_bound", "Bortezomib-20S Complex (nM)", color="#8e44ad")
  })
  output$plt_cy_pk <- renderPlotly({
    df <- sim_result()
    mk_line(df, "CY_Cp", "4-OH-CY (µM)", color="#f39c12")
  })
  output$plt_dex_pk <- renderPlotly({
    df <- sim_result()
    mk_line(df, "DEX_Cp", "Dexamethasone (nM)", color="#27ae60")
  })

  output$pk_table <- renderDT({
    df <- data.frame(
      Drug = c("Daratumumab","Bortezomib","Cyclophosphamide","Dexamethasone","Melphalan"),
      PK_Model = c("2-CMT TMDD","1-CMT + proteasome binding","1-CMT (active metabolite)","1-CMT","1-CMT"),
      t_half = c("18-23d (TMDD-dependent)","~9-15h","~4-6h (4-OH-CY)","~3-4h","~1.5h"),
      Vd = c("56+40 mL/kg","498 L","38 L","64 L","28 L"),
      CL = c("3.1 mL/day/kg","9.2 L/h","5.8 L/h","15.4 L/h","8.5 L/h")
    )
    datatable(df, options=list(dom="t"), rownames=FALSE)
  })

  ## Hematologic PD
  output$plt_pc <- renderPlotly({
    df <- sim_result()
    mk_line(df, "PC", "Plasma Cell Pool (1 = Baseline)", color="#e74c3c")
  })
  output$plt_dflc <- renderPlotly({
    df <- sim_result()
    mk_line(df, "dFLC", "dFLC (mg/L)", hline=40, hline_label="CR threshold (40 mg/L)", color="#c0392b")
  })
  output$plt_nk <- renderPlotly({
    df <- sim_result()
    mk_line(df, "NK", "NK Cell Pool (1 = Baseline)", color="#16a085")
  })
  output$heme_resp_table <- renderDT({
    df <- sim_result()
    last_val <- tail(df, 1)
    cr  <- last_val$dFLC < 40
    resp <- data.frame(
      Category = c("CR","VGPR","PR","NR"),
      Definition = c("dFLC <40 mg/L","dFLC <40 mg/L OR >90% reduction","≥50% dFLC reduction","<50% dFLC reduction"),
      Achieved = c(
        ifelse(cr, "YES","No"),
        ifelse(last_val$dFLC/input$FLC_init < 0.1 || cr, "YES","No"),
        ifelse(last_val$dFLC/input$FLC_init < 0.5, "YES","No"),
        ifelse(last_val$dFLC/input$FLC_init >= 0.5, "YES","No")
      )
    )
    datatable(df, options=list(dom="t"), rownames=FALSE)
  })

  ## Organ biomarkers
  output$plt_bnp <- renderPlotly({
    df <- sim_result()
    mk_line(df, "NTproBNP", "NT-proBNP (pg/mL)",
            hline=input$NTproBNP_init * 0.70, hline_label="30% reduction (organ response)",
            color="#e74c3c")
  })
  output$plt_tnt <- renderPlotly({
    df <- sim_result()
    mk_line(df, "TnT", "hs-Troponin T (ng/mL)",
            hline=0.025, hline_label="Mayo RF threshold (0.025 ng/mL)", color="#c0392b")
  })
  output$plt_egfr <- renderPlotly({
    df <- sim_result()
    mk_line(df, "eGFR", "eGFR (mL/min/1.73m²)",
            hline=15, hline_label="ESKD threshold (15 mL/min)", color="#2980b9")
  })
  output$plt_prot <- renderPlotly({
    df <- sim_result()
    mk_line(df, "ProtUria", "24h Proteinuria (g/day)",
            hline=3.5, hline_label="Nephrotic threshold (3.5 g/day)", color="#27ae60")
  })

  ## Clinical endpoints
  output$plt_amyl_card <- renderPlotly({
    df <- sim_result()
    mk_line(df, "AmylCard", "Cardiac Amyloid Burden (1 = Baseline)", color="#c0392b")
  })
  output$plt_amyl_ren <- renderPlotly({
    df <- sim_result()
    mk_line(df, "AmylRen", "Renal Amyloid Burden (1 = Baseline)", color="#2980b9")
  })
  output$plt_org_resp <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time) |>
      add_lines(y=~CardResp, name="Cardiac Organ Response", line=list(color="#e74c3c", width=2)) |>
      add_lines(y=~RenResp,  name="Renal Organ Response",  line=list(color="#2980b9", width=2)) |>
      add_lines(y=~CR_flag,  name="Hematologic CR",        line=list(color="#27ae60", width=2)) |>
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Response (1=Yes)", range=c(-0.05,1.2)))
  })

  output$trial_table <- renderDT({
    df <- data.frame(
      Trial = c("ANDROMEDA","ANDROMEDA","MDex","VCD/CyBorD","ISA220"),
      Regimen = c("Dara+CyBorD","CyBorD","Melphalan+Dex","CyBorD","Isa+VCd"),
      CR_Rate = c("53.3%","18.1%","~33%","~29%","~45%"),
      CardResp = c("~42%","~22%","~26%","~18%","~40%"),
      Source = c("NEJM 2021","NEJM 2021","Blood 2004","Blood 2012","Lancet 2023")
    )
    datatable(df, options=list(dom="t"), rownames=FALSE)
  })

  ## Scenario comparison
  output$cmp_dflc <- renderPlotly({
    df <- compare_results()
    cols <- SCENARIO_COLORS[unique(df$Scenario)]
    p <- plot_ly()
    for (s in unique(df$Scenario)) {
      sub <- filter(df, Scenario==s)
      p <- add_lines(p, x=sub$time, y=sub$dFLC, name=s,
                     line=list(color=cols[s], width=2.5))
    }
    p |> add_segments(x=0, xend=max(df$time), y=40, yend=40,
                      line=list(color="black", dash="dash"), showlegend=FALSE) |>
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="dFLC (mg/L)"),
             legend=list(orientation="h"))
  })
  output$cmp_bnp <- renderPlotly({
    df <- compare_results()
    cols <- SCENARIO_COLORS[unique(df$Scenario)]
    p <- plot_ly()
    for (s in unique(df$Scenario)) {
      sub <- filter(df, Scenario==s)
      p <- add_lines(p, x=sub$time, y=sub$NTproBNP, name=s,
                     line=list(color=cols[s], width=2.5))
    }
    p |> layout(xaxis=list(title="Time (days)"), yaxis=list(title="NT-proBNP (pg/mL)"),
                legend=list(orientation="h"))
  })
  output$cmp_egfr <- renderPlotly({
    df <- compare_results()
    cols <- SCENARIO_COLORS[unique(df$Scenario)]
    p <- plot_ly()
    for (s in unique(df$Scenario)) {
      sub <- filter(df, Scenario==s)
      p <- add_lines(p, x=sub$time, y=sub$eGFR, name=s,
                     line=list(color=cols[s], width=2.5))
    }
    p |> layout(xaxis=list(title="Time (days)"), yaxis=list(title="eGFR (mL/min)"),
                legend=list(orientation="h"))
  })
  output$cmp_pc <- renderPlotly({
    df <- compare_results()
    cols <- SCENARIO_COLORS[unique(df$Scenario)]
    p <- plot_ly()
    for (s in unique(df$Scenario)) {
      sub <- filter(df, Scenario==s)
      p <- add_lines(p, x=sub$time, y=sub$PC, name=s,
                     line=list(color=cols[s], width=2.5))
    }
    p |> layout(xaxis=list(title="Time (days)"), yaxis=list(title="Plasma Cell Pool"),
                legend=list(orientation="h"))
  })

  output$compare_table <- renderDT({
    df <- compare_results() |>
      filter(time == min(input$n_days, 180)) |>
      group_by(Scenario) |>
      summarise(
        dFLC_mgL    = round(mean(dFLC), 1),
        NT_proBNP   = round(mean(NTproBNP), 0),
        TnT_ng      = round(mean(TnT), 4),
        eGFR        = round(mean(eGFR), 1),
        PC_fraction = round(mean(PC), 3),
        CR_pct      = paste0(round(mean(CR_flag)*100, 1), "%"),
        CardResp_pct= paste0(round(mean(CardResp)*100, 1), "%"),
        .groups="drop"
      )
    datatable(df, options=list(dom="t", scrollX=TRUE), rownames=FALSE)
  })

  ## Staging
  output$plt_stage <- renderPlotly({
    df <- sim_result()
    df <- df |> mutate(
      stage_score =
        as.integer(TnT > 0.025) +
        as.integer(NTproBNP > 1800) +
        as.integer(dFLC > 180) +
        as.integer(eGFR < 50),
      Mayo_Stage = stage_score + 1
    )
    plot_ly(df, x=~time, y=~Mayo_Stage, type="scatter", mode="lines",
            line=list(color="#e67e22", width=2.5)) |>
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="Mayo 2012 Stage", range=c(0.8, 4.2),
                        tickvals=1:4, ticktext=paste("Stage", 1:4)))
  })

  output$plt_rf_gauge <- renderPlotly({
    df <- sim_result()
    last <- tail(df, 1)
    rf_values <- c(
      "hs-TnT ≥0.025"    = as.integer(last$TnT > 0.025),
      "NT-proBNP ≥1800"  = as.integer(last$NTproBNP > 1800),
      "dFLC ≥18 mg/dL"   = as.integer(last$dFLC > 180),
      "eGFR <50"         = as.integer(last$eGFR < 50)
    )
    plot_ly(
      x = names(rf_values),
      y = rf_values,
      type = "bar",
      marker = list(color = ifelse(rf_values == 1, "#e74c3c", "#2ecc71"))
    ) |>
      layout(title="Risk Factor Status (End of Simulation)",
             yaxis=list(title="Present (1=Yes, 0=No)", range=c(0,1.2)),
             xaxis=list(title=""))
  })

  output$staging_monitor_table <- renderDT({
    df <- sim_result()
    timepoints <- c(30, 60, 90, 180, 270, 365)
    timepoints <- timepoints[timepoints <= input$n_days]
    monitor <- lapply(timepoints, function(t) {
      row <- df[which.min(abs(df$time - t)), ]
      data.frame(
        Day        = t,
        dFLC_mgL   = round(row$dFLC, 1),
        NTproBNP   = round(row$NTproBNP, 0),
        TnT_ng     = round(row$TnT, 4),
        eGFR       = round(row$eGFR, 1),
        CR         = ifelse(row$dFLC < 40, "YES", "No"),
        CardResp   = ifelse(row$CardResp == 1, "YES", "No"),
        RenResp    = ifelse(row$RenResp == 1, "YES", "No")
      )
    })
    datatable(do.call(rbind, monitor), options=list(dom="t"), rownames=FALSE)
  })

  ## Safety
  output$safety_table <- renderDT({
    df <- data.frame(
      AE = c(
        "Infusion-Related Reaction (IRR)",
        "Neutropenia (G3-4)",
        "Thrombocytopenia",
        "Peripheral Neuropathy (BTZ)",
        "Herpes Zoster Reactivation",
        "DVT/PE (thromboembolic)",
        "Hyperglycemia (DEX)",
        "Immunoparesis"
      ),
      Drug = c(
        "Daratumumab",
        "Cyclophosphamide / Dara",
        "Bortezomib",
        "Bortezomib (SC preferred to ↓PN)",
        "Bortezomib",
        "Dexamethasone + IMiDs",
        "Dexamethasone",
        "Daratumumab (long-term)"
      ),
      Frequency_ANDROMEDA = c(
        "48% G1-2 / 2% G3",
        "17% G3-4",
        "15% G3-4",
        "8% (SC route preferred)",
        "6% (prophylaxis required)",
        "2%",
        "15% G2-3",
        "16% (prolonged therapy)"
      ),
      Management = c(
        "Pre-medicate: corticosteroid, antihistamine, APAP; slow infusion rate",
        "G-CSF support; dose reduction if G4",
        "BTZ dose hold/reduce",
        "SC route; dose reduction; vitamin B6",
        "Acyclovir/valacyclovir prophylaxis mandatory",
        "Aspirin 81mg; LMWH if high risk",
        "Glucose monitoring; insulin if needed",
        "IVIg if recurrent severe infections"
      )
    )
    datatable(df, options=list(dom="t", scrollX=TRUE), rownames=FALSE)
  })

  output$plt_irr <- renderPlotly({
    df <- sim_result()
    irr_risk <- ifelse(df$DARA_Cp > 0, 0.48 * exp(-df$time / 30), 0)
    irr_risk <- pmin(irr_risk, 0.50)
    plot_ly(df, x=~time, y=irr_risk, type="scatter", mode="lines",
            line=list(color="#e74c3c", width=2)) |>
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="Cumulative IRR Risk", range=c(0,0.55)))
  })
  output$plt_nk_safety <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~NK, type="scatter", mode="lines",
            line=list(color="#16a085", width=2)) |>
      add_segments(x=0, xend=max(df$time), y=0.5, yend=0.5,
                   line=list(color="orange", dash="dash"), showlegend=FALSE) |>
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="NK Cell Pool (1=Normal)", range=c(0,1.1)))
  })

  output$ddi_table <- renderDT({
    df <- data.frame(
      Interaction = c(
        "Daratumumab + blood type serology",
        "Bortezomib + CYP3A4 inhibitors (azoles)",
        "Dexamethasone + warfarin",
        "Cyclophosphamide + allopurinol",
        "Daratumumab + anti-CD38 Ab (interference)"
      ),
      Mechanism = c(
        "Dara binds CD38 on RBCs → indirect antiglobulin test (IAT) interference",
        "CYP3A4 inhibition → ↑BTZ AUC → ↑neurotoxicity",
        "DEX induces CYP2C9 → ↑warfarin metabolism → need dose adjustment",
        "Allopurinol inhibits CY activation → ↓efficacy",
        "Cross-reactivity with anti-CD38 mAbs (isatuximab, daratumumab)"
      ),
      Clinical_Action = c(
        "Notify blood bank; type/screen BEFORE daratumumab; use phenotyping",
        "Avoid strong CYP3A4 inhibitors; dose-reduce BTZ if unavoidable",
        "INR monitoring weekly for 4 weeks after steroid start/stop",
        "Avoid concurrent allopurinol; use rasburicase for TLS prophylaxis",
        "24-month treatment-free interval before switching anti-CD38 agent"
      )
    )
    datatable(df, options=list(dom="t", scrollX=TRUE), rownames=FALSE)
  })

}

## ============================================================
## RUN APP
## ============================================================
shinyApp(ui = ui, server = server)
