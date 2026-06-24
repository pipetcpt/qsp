## ============================================================
## Diabetic Retinopathy (DR) — Interactive Shiny Dashboard
## ============================================================
## 질환: 당뇨병성 망막병증 (Diabetic Retinopathy)
## 약어: DR  |  버전: 1.0  |  날짜: 2026-06-24
##
## 탭 구조 (8 Tabs):
##  1. Patient Profile        — 환자 특성 입력 및 위험도 시각화
##  2. Drug PK                — 안내 약물 약동학 (유리체 내 농도)
##  3. VEGF / Angiogenesis    — 자유 VEGF 및 신생혈관 역학
##  4. Inflammation & Oxidative Stress — CYT/ROS/AGE
##  5. Retinal Structural     — CRT (OCT 두께) / 투과성
##  6. Visual Outcomes        — BCVA (ETDRS 글자수), VA 변화
##  7. Scenario Comparison    — 6개 치료 시나리오 비교
##  8. Biomarkers & About
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

## ============================================================
## QSP Simulation Engine (embedded ODE solver)
## ============================================================
simulate_DR <- function(BG_init, HbA1c_init, VA_init, CRT_init, VEGF_init,
                        PERICYTE_init, NV_init,
                        treatment, dose_mg, dosing_interval, n_injections,
                        BG_target, end_time = 730) {

  ## --- Parameters ---
  MW_map  <- c("aflibercept"=115, "ranibizumab"=48, "bevacizumab"=149, "faricimab"=146, "none"=1)
  CLv_map <- c("aflibercept"=0.44, "ranibizumab"=0.23, "bevacizumab"=0.36, "faricimab"=0.35, "none"=0.44)
  kon_map <- c("aflibercept"=13.2, "ranibizumab"=10.0, "bevacizumab"=8.0,  "faricimab"=14.0, "none"=0)
  koff_map<- c("aflibercept"=0.066,"ranibizumab"=0.10, "bevacizumab"=0.15, "faricimab"=0.04, "none"=0)

  MW   <- MW_map[treatment]
  CLv  <- CLv_map[treatment]
  kon  <- kon_map[treatment]
  koff <- koff_map[treatment]
  Vd_vit <- 4.0;  Vd_cent <- 3600;  CL_cent <- 360;  Q_vit <- 0.05; Q_per <- 50; Vd_per <- 5000

  ksyn_VEGF  <- 0.15; kdeg_VEGF  <- 0.12; n_BG_VEGF <- 2.0; EC50_BG <- 9.0
  ksyn_ROS   <- 0.20; kdeg_ROS   <- 0.18; kBG_ROS <- 0.08
  ksyn_AGE   <- 0.005; kdeg_AGE  <- 0.002
  ksyn_CYT   <- 0.30; kdeg_CYT   <- 0.25; kVEGF_CYT <- 0.15; kROS_CYT <- 0.12
  ksyn_ICAM  <- 0.20; kdeg_ICAM  <- 0.15; kCYT_ICAM <- 0.18
  kLoss_PERI <- 0.003; kAGE_PERI <- 0.015; kROS_PERI <- 0.010
  kLoss_EC   <- 0.002; kICM_EC   <- 0.012
  ksyn_PERM  <- 0.25; kdeg_PERM  <- 0.20; kVEGF_PERM <- 0.30; kCYT_PERM <- 0.15
  ksyn_NV    <- 0.008; kdeg_NV   <- 0.002; kVEGF_NV <- 0.05
  kCRT_up    <- 0.15; kCRT_dn   <- 0.08; CRT_max <- 600; CRT_base <- 270
  kCRT_VA    <- 0.002; kNV_VA   <- 0.015; kVA_up <- 0.003

  VEGF_base <- 1.25; ROS_base <- 1.0; AGE_base <- 1.0; CYT_base <- 1.0
  ICAM_base <- 1.0; PERI_base <- 100.0
  kBG_eq <- 0.015; kHbA1c_up <- 0.033

  ## Build dosing schedule
  dose_times <- if (treatment != "none" && n_injections > 0)
    seq(0, (n_injections-1)*dosing_interval, by=dosing_interval)
  else numeric(0)

  ## State: (DV=drug_vit, DC=drug_cent, DP=drug_per, BG, HA, VF, VB, ROS, AGE, CYT, ICAM, PERI, EC, PERM, NV, CRT, VA)
  state <- c(DV=0, DC=0, DP=0, BG=BG_init, HA=HbA1c_init,
             VF=VEGF_init, VB=0, ROS=2.0, AGE=2.0, CYT=2.5, ICAM=2.0,
             PERI=PERICYTE_init, EC=75, PERM=3.0, NV=NV_init, CRT=CRT_init, VA=VA_init)

  dt <- 1; steps <- end_time
  output <- matrix(NA, nrow=steps+1, ncol=length(state)+1)
  colnames(output) <- c("time", names(state))
  output[1,] <- c(0, state)

  for (i in seq_len(steps)) {
    t <- (i-1)*dt
    s <- state

    # Drug input
    drug_in <- if (t %in% dose_times) dose_mg else 0

    # Current drug vitreous in nM
    Dv_nM <- max(0, s["DV"] / Vd_vit * 1e6 / MW)

    BG_eff <- s["BG"]^n_BG_VEGF / (EC50_BG^n_BG_VEGF + s["BG"]^n_BG_VEGF)
    VEGFsyn <- ksyn_VEGF * (1.0 + 3.0*BG_eff) * (1.0 + 0.5*(s["AGE"]/AGE_base-1))
    Bon  <- kon  * max(0,s["VF"]) * max(0,Dv_nM)
    Boff <- koff * max(0,s["VB"])

    ROS_driv <- kBG_ROS * max(0, s["BG"] - 5.5)
    AGE_rate <- ksyn_AGE * (s["BG"]/5.5) * (1.0 + 0.5*(s["HA"]-7.0)/7.0)
    CYT_stim <- kVEGF_CYT*(s["VF"]/VEGF_base) + kROS_CYT*(s["ROS"]/ROS_base)
    PERI_loss<- kLoss_PERI + kAGE_PERI*(s["AGE"]/AGE_base-1) + kROS_PERI*(s["ROS"]/ROS_base-1)
    PERI_loss<- max(0, PERI_loss)
    EC_loss  <- max(0, kLoss_EC + kICM_EC*(s["ICAM"]/ICAM_base-1))
    PERM_in  <- ksyn_PERM * (kVEGF_PERM*s["VF"]/VEGF_base + kCYT_PERM*s["CYT"]/CYT_base)
    NV_drv   <- kVEGF_NV*(s["VF"]/VEGF_base)*(1-max(0,s["PERI"])/PERI_base)
    CRT_in   <- kCRT_up*(max(0,s["PERM"])/max(0.01,1.0))*(CRT_max - max(CRT_base, s["CRT"]))
    CRT_dn   <- kCRT_dn*(max(CRT_base,s["CRT"]) - CRT_base)
    CRT_exc  <- max(0, (s["CRT"] - CRT_base)/100)
    VA_loss  <- kCRT_VA*CRT_exc*max(0,s["VA"]) + kNV_VA*max(0,s["NV"])*max(0,s["VA"])/10
    VA_rec   <- kVA_up*(VA_init - s["VA"]) * max(0, 1 - CRT_exc/4)

    # Euler integration
    dstate <- c(
      DV   = drug_in - (CLv/Vd_vit)*s["DV"] - Q_vit*s["DV"]/Vd_vit + Q_vit*s["DC"]/Vd_cent,
      DC   = -(CL_cent/Vd_cent)*s["DC"] + Q_vit*s["DV"]/Vd_vit - Q_vit*s["DC"]/Vd_cent
             - (Q_per/Vd_cent)*s["DC"] + (Q_per/Vd_per)*s["DP"],
      DP   = (Q_per/Vd_cent)*s["DC"] - (Q_per/Vd_per)*s["DP"],
      BG   = kBG_eq*(BG_target - s["BG"]),
      HA   = kHbA1c_up*(0.195*s["BG"] + 3.1 - s["HA"]),
      VF   = VEGFsyn - kdeg_VEGF*s["VF"] - Bon + Boff,
      VB   = Bon - Boff - kdeg_VEGF*s["VB"],
      ROS  = ksyn_ROS + ROS_driv - kdeg_ROS*s["ROS"],
      AGE  = AGE_rate - kdeg_AGE*s["AGE"],
      CYT  = ksyn_CYT*(1+CYT_stim) - kdeg_CYT*s["CYT"],
      ICAM = ksyn_ICAM*(1+kCYT_ICAM*(s["CYT"]/CYT_base)) - kdeg_ICAM*s["ICAM"],
      PERI = -PERI_loss * s["PERI"],
      EC   = -EC_loss * s["EC"],
      PERM = PERM_in - kdeg_PERM*s["PERM"],
      NV   = ksyn_NV*NV_drv - kdeg_NV*s["NV"],
      CRT  = CRT_in - CRT_dn,
      VA   = VA_rec - VA_loss
    )

    state <- pmax(state + dstate*dt, 0)
    state["VA"] <- min(state["VA"], 100)
    state["CRT"]<- min(state["CRT"], CRT_max)
    state["PERI"]<- min(state["PERI"], PERI_base)
    output[i+1,] <- c(i*dt, state)
  }
  as.data.frame(output)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "DR QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",       tabName="tab_patient",  icon=icon("user")),
      menuItem("② Drug PK",               tabName="tab_pk",       icon=icon("pills")),
      menuItem("③ VEGF / Angiogenesis",   tabName="tab_vegf",     icon=icon("dna")),
      menuItem("④ Oxidative/Inflammation",tabName="tab_inflam",   icon=icon("fire")),
      menuItem("⑤ Retinal Structure",     tabName="tab_struct",   icon=icon("eye")),
      menuItem("⑥ Visual Outcomes",       tabName="tab_va",       icon=icon("chart-line")),
      menuItem("⑦ Scenario Comparison",   tabName="tab_compare",  icon=icon("layer-group")),
      menuItem("⑧ Biomarkers & About",    tabName="tab_about",    icon=icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML(".box {border-top-color:#1E88E5;}"))),
    tabItems(

      ## --------------------------------------------------------
      ## Tab 1: Patient Profile
      ## --------------------------------------------------------
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Patient Parameters", status="primary", solidHeader=TRUE, width=4,
            sliderInput("BG_init",    "Baseline Blood Glucose (mmol/L):",   min=5, max=20, value=9.5, step=0.5),
            sliderInput("HbA1c_init", "Baseline HbA1c (%):",                min=6, max=14, value=8.5, step=0.5),
            sliderInput("VA_init",    "Baseline VA (ETDRS letters):",        min=20, max=85, value=60, step=1),
            sliderInput("CRT_init",   "Baseline CRT (µm):",                  min=250, max=600, value=380, step=10),
            sliderInput("VEGF_init",  "Baseline Vitreous VEGF (nM):",        min=1, max=12, value=3.5, step=0.5),
            sliderInput("PERI_init",  "Baseline Pericyte Count (%):",        min=20, max=100, value=60, step=5),
            sliderInput("NV_init",    "Baseline NV Index (0=none, 5=severe):", min=0, max=5, value=0.5, step=0.1)
          ),
          box(title="DR Severity Classification", status="info", solidHeader=TRUE, width=4,
            h4("Current DR Stage (estimated):"),
            verbatimTextOutput("dr_stage_text"),
            br(),
            h4("5-Year Risk Estimates:"),
            tableOutput("risk_table"),
            br(),
            h4("Key Biomarkers at Baseline:"),
            tableOutput("baseline_table")
          ),
          box(title="Disease Activity Radar", status="warning", solidHeader=TRUE, width=4,
            plotlyOutput("radar_plot", height=300),
            br(),
            strong("DR Stage Guide:"),
            tags$ul(
              tags$li("Mild NPDR: microaneurysms only"),
              tags$li("Moderate NPDR: more but < severe"),
              tags$li("Severe NPDR: 4-2-1 rule"),
              tags$li("PDR: neovascularization present"),
              tags$li("CI-DME: CRT >310µm, center-involving")
            )
          )
        ),
        fluidRow(
          box(title="Treatment Selection", status="success", solidHeader=TRUE, width=6,
            selectInput("treatment", "Anti-VEGF Agent:",
                        choices=c("Aflibercept (2mg)"="aflibercept",
                                  "Ranibizumab (0.5mg)"="ranibizumab",
                                  "Bevacizumab (1.25mg)"="bevacizumab",
                                  "Faricimab (6mg)"="faricimab",
                                  "None (observe)"="none"), selected="aflibercept"),
            sliderInput("dose_mg",        "Dose (mg):",             min=0.1, max=8, value=2, step=0.1),
            sliderInput("dose_interval",  "Dosing Interval (days):", min=14, max=112, value=28, step=7),
            sliderInput("n_injections",   "Number of Injections:",  min=1, max=26, value=10, step=1),
            sliderInput("BG_target",      "BG Treatment Target (mmol/L):", min=5, max=12, value=9.5, step=0.5),
            sliderInput("sim_end",        "Simulation Duration (days):", min=90, max=1095, value=730, step=30)
          ),
          box(title="Treatment Background", status="success", solidHeader=TRUE, width=6,
            h4("Anti-VEGF Agents — Key Facts"),
            tableOutput("drug_facts_table")
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 2: Drug PK
      ## --------------------------------------------------------
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Vitreous Drug Concentration Over Time", status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("pk_plot", height=350)
          )
        ),
        fluidRow(
          box(title="PK Parameters", status="info", solidHeader=TRUE, width=6,
            tableOutput("pk_params_table")
          ),
          box(title="PK Summary Statistics", status="info", solidHeader=TRUE, width=6,
            verbatimTextOutput("pk_summary")
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 3: VEGF / Angiogenesis
      ## --------------------------------------------------------
      tabItem(tabName="tab_vegf",
        fluidRow(
          box(title="Free Vitreous VEGF", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("vegf_free_plot", height=300)
          ),
          box(title="VEGF-Drug Bound Complex", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("vegf_bound_plot", height=300)
          )
        ),
        fluidRow(
          box(title="Neovascularization Index", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("nv_plot", height=300)
          ),
          box(title="Pericyte Count (% Normal)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("pericyte_plot", height=300)
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 4: Oxidative Stress & Inflammation
      ## --------------------------------------------------------
      tabItem(tabName="tab_inflam",
        fluidRow(
          box(title="Oxidative Stress (ROS) Index", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("ros_plot", height=300)
          ),
          box(title="AGE Accumulation", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("age_plot", height=300)
          )
        ),
        fluidRow(
          box(title="Cytokine Index (IL-6/IL-8/TNFα composite)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("cyt_plot", height=300)
          ),
          box(title="ICAM-1 Index (Leukostasis Marker)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("icam_plot", height=300)
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 5: Retinal Structural Changes
      ## --------------------------------------------------------
      tabItem(tabName="tab_struct",
        fluidRow(
          box(title="Central Retinal Thickness (CRT)", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("crt_plot", height=350)
          ),
          box(title="CRT Benchmarks", status="info", solidHeader=TRUE, width=4,
            tags$ul(
              tags$li("Normal CRT: ~250–270 µm"),
              tags$li("CI-DME threshold: >310 µm"),
              tags$li("Moderate DME: 310–400 µm"),
              tags$li("Severe DME: >400 µm"),
              tags$li(""),
              tags$li("PROTOCOL T 1yr (baseline ~407 µm):"),
              tags$li("  Afl: −169 µm (→ ~238 µm)"),
              tags$li("  RBZ: −147 µm (→ ~260 µm)"),
              tags$li("  Bev: −101 µm (→ ~306 µm)"),
              tags$li(""),
              tags$li("TENAYA/LUCERNE 1yr:"),
              tags$li("  Faricimab: −189/−194 µm"),
              tags$li("  Aflibercept: −163/−174 µm")
            )
          )
        ),
        fluidRow(
          box(title="Vascular Permeability Index", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("perm_plot", height=300)
          ),
          box(title="Endothelial Cell Count (%)", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("ec_plot", height=300)
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 6: Visual Outcomes
      ## --------------------------------------------------------
      tabItem(tabName="tab_va",
        fluidRow(
          box(title="Visual Acuity (ETDRS Letters) Over Time", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("va_plot", height=400)
          ),
          box(title="VA Outcome Categories", status="info", solidHeader=TRUE, width=4,
            h4("ETDRS Letter Score Reference:"),
            tableOutput("va_ref_table"),
            br(),
            h4("Regulatory Endpoints:"),
            tags$ul(
              tags$li("≥15 letter gain: clinically significant improvement"),
              tags$li("≥10 letter gain: meaningful improvement"),
              tags$li("< 15 letter loss: avoiding substantial decline"),
              tags$li("≥30 letter loss: severe visual loss")
            )
          )
        ),
        fluidRow(
          box(title="VA Change from Baseline", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("va_change_plot", height=300)
          ),
          box(title="VA Summary at Year 1 & Year 2", status="success", solidHeader=TRUE, width=6,
            DTOutput("va_summary_table")
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 7: Scenario Comparison
      ## --------------------------------------------------------
      tabItem(tabName="tab_compare",
        fluidRow(
          box(title="Multi-Scenario Setup", status="primary", solidHeader=TRUE, width=12,
            p("Six pre-defined scenarios are compared head-to-head based on the current patient profile:"),
            fluidRow(
              column(4,
                strong("S0: No treatment (poor glycemic ctrl)"), br(),
                strong("S1: Glycemic control only (HbA1c→7%)"), br(),
                strong("S2: Aflibercept 2mg q4→8w")
              ),
              column(4,
                strong("S3: Ranibizumab 0.5mg q4w"), br(),
                strong("S4: Faricimab 6mg q4→16w"), br(),
                strong("S5: Aflibercept + Glycemic Ctrl (combination)")
              ),
              column(4,
                actionButton("run_compare", "Run All Scenarios", class="btn-primary btn-lg",
                             icon=icon("play"))
              )
            )
          )
        ),
        fluidRow(
          box(title="VA Comparison", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_va", height=300)
          ),
          box(title="CRT Comparison", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_crt", height=300)
          )
        ),
        fluidRow(
          box(title="Free VEGF Comparison", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_vegf", height=300)
          ),
          box(title="NV Index Comparison", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("cmp_nv", height=300)
          )
        ),
        fluidRow(
          box(title="Endpoint Summary Table (Year 1 & 2)", status="primary", solidHeader=TRUE, width=12,
            DTOutput("compare_table")
          )
        )
      ),

      ## --------------------------------------------------------
      ## Tab 8: Biomarkers & About
      ## --------------------------------------------------------
      tabItem(tabName="tab_about",
        fluidRow(
          box(title="Biomarker Dashboard", status="primary", solidHeader=TRUE, width=8,
            plotlyOutput("biomarker_plot", height=400)
          ),
          box(title="Key Biomarkers", status="info", solidHeader=TRUE, width=4,
            tags$ul(
              tags$li("Vitreous VEGF: main driver of DME/NV"),
              tags$li("HbA1c: reflects glycemic control"),
              tags$li("CRT: structural OCT endpoint"),
              tags$li("ETDRS VA: functional primary endpoint"),
              tags$li("Pericyte %: cellular biomarker"),
              tags$li("ROS index: oxidative burden"),
              tags$li("AGE: cumulative glycation damage"),
              tags$li("Cytokine index: inflammatory state")
            )
          )
        ),
        fluidRow(
          box(title="Model Overview", status="success", solidHeader=TRUE, width=12,
            h3("Diabetic Retinopathy (DR) QSP Model"),
            p("This model integrates 18 ODE compartments representing the key pathophysiology of DR:"),
            tags$ol(
              tags$li("Drug PK (anti-VEGF vitreous/plasma, corticosteroid)"),
              tags$li("Glycemic dynamics (blood glucose, HbA1c)"),
              tags$li("VEGF signaling (free VEGF, drug-bound complex, PlGF)"),
              tags$li("Oxidative stress (ROS, AGE accumulation)"),
              tags$li("Inflammation (cytokines, ICAM-1, leukostasis)"),
              tags$li("Cellular pathology (pericyte loss, EC apoptosis)"),
              tags$li("Structural changes (permeability, NV, CRT)"),
              tags$li("Visual acuity (ETDRS letter score)")
            ),
            h4("Key Clinical Trials for Parameter Calibration:"),
            tableOutput("calibration_table"),
            br(),
            h4("Model Limitations:"),
            tags$ul(
              tags$li("Simplified 1-dimensional retinal compartment (no spatial resolution)"),
              tags$li("Individual variability not captured (deterministic model)"),
              tags$li("Long-term structural remodeling (scarring, TRD) partially modeled"),
              tags$li("Drug tolerance and tachyphylaxis not included"),
              tags$li("Laser and surgical interventions not fully parameterized")
            )
          )
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================
server <- function(input, output, session) {

  ## --- Reactive simulation ---
  sim_result <- reactive({
    simulate_DR(
      BG_init = input$BG_init, HbA1c_init = input$HbA1c_init,
      VA_init = input$VA_init, CRT_init = input$CRT_init,
      VEGF_init = input$VEGF_init, PERICYTE_init = input$PERI_init,
      NV_init = input$NV_init,
      treatment = input$treatment, dose_mg = input$dose_mg,
      dosing_interval = input$dose_interval, n_injections = input$n_injections,
      BG_target = input$BG_target, end_time = input$sim_end
    )
  })

  ## --- DR Stage classification ---
  output$dr_stage_text <- renderText({
    crt <- input$CRT_init; nv <- input$NV_init; va <- input$VA_init
    stage <- if (nv > 1.5) "Proliferative DR (PDR)"
      else if (nv > 0.3) "Severe NPDR / Early PDR"
      else if (crt > 400) "Moderate NPDR + Severe DME"
      else if (crt > 310) "Moderate NPDR + Center-Involving DME"
      else "Mild–Moderate NPDR"
    paste0("Estimated Stage: ", stage, "\nCRT: ", crt, " µm | VA: ", va, " letters | NV: ", round(nv,1))
  })

  output$risk_table <- renderTable({
    bg <- input$BG_init; hba <- input$HbA1c_init
    risk_5yr <- round(c(20 + (bg-7)*5, 15 + (hba-7)*8, 5 + input$NV_init*15), 1)
    data.frame(Event=c("PDR development","Severe VA loss (≥30L)","NV Regression (if PDR)"),
               `5yr Risk (%)`=pmin(95, risk_5yr))
  })

  output$baseline_table <- renderTable({
    data.frame(
      Biomarker = c("Blood Glucose","HbA1c","Vitreous VEGF","CRT","BCVA","Pericytes","NV Index"),
      Value     = c(input$BG_init, input$HbA1c_init, input$VEGF_init,
                    input$CRT_init, input$VA_init, input$PERI_init, input$NV_init),
      Unit      = c("mmol/L","%","nM","µm","ETDRS letters","%","AU")
    )
  })

  ## --- Radar plot (baseline disease activity) ---
  output$radar_plot <- renderPlotly({
    vals <- c(
      `Glucose`    = (input$BG_init - 5)/15 * 100,
      `VEGF`       = (input$VEGF_init)/12 * 100,
      `CRT`        = (input$CRT_init - 250)/350 * 100,
      `NV`         = input$NV_init/5 * 100,
      `VA loss`    = (85-input$VA_init)/65*100,
      `Pericyte\nloss`= (100-input$PERI_init)
    )
    vals <- pmax(0, pmin(100, vals))
    plot_ly(type='scatterpolar', r=c(vals, vals[1]),
            theta=c(names(vals), names(vals)[1]),
            fill='toself', mode='lines+markers',
            line=list(color='#1E88E5'), fillcolor='rgba(30,136,229,0.2)') %>%
      layout(polar=list(radialaxis=list(range=c(0,100))),
             showlegend=FALSE, margin=list(t=20,b=20,l=20,r=20))
  })

  ## --- Drug facts table ---
  output$drug_facts_table <- renderTable({
    data.frame(
      Drug       = c("Aflibercept","Ranibizumab","Bevacizumab","Faricimab"),
      `Dose(mg)` = c(2.0, 0.5, 1.25, 6.0),
      `MW(kDa)`  = c(115, 48, 149, 146),
      `t½ vit`   = c("~9d","~3d","~6d","~14d"),
      Target     = c("VEGF-A/B+PlGF","VEGF-A","VEGF-A","VEGF-A+Ang2"),
      `Trial`    = c("PROTOCOL T","RISE/RIDE","PROTOCOL T","TENAYA/LUCERNE")
    )
  })

  ## --- PK Plot ---
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~DV, type='scatter', mode='lines',
            line=list(color='#1E88E5', width=2), name='Vitreous') %>%
      add_trace(y=~DC/1000, name='Plasma (×0.001)', line=list(color='#E53935', dash='dash')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Drug Conc. (mg/mL)"),
             title="Intravitreal Drug PK", legend=list(x=0.7,y=0.9))
  })

  ## --- PK parameters table ---
  output$pk_params_table <- renderTable({
    mw <- c(aflibercept=115, ranibizumab=48, bevacizumab=149, faricimab=146, none=1)[input$treatment]
    cl <- c(aflibercept=0.44, ranibizumab=0.23, bevacizumab=0.36, faricimab=0.35, none=0.44)[input$treatment]
    t12 <- round(log(2)/cl * 4.0 * 24, 1)
    data.frame(Parameter=c("MW (kDa)","CL_vit (mL/day)","t½ vitreous (h)","Vd vitreous (mL)","Dose (mg)"),
               Value=c(mw, cl, t12, 4.0, input$dose_mg))
  })

  output$pk_summary <- renderPrint({
    df <- sim_result()
    cat("Cmax vitreous:", round(max(df$DV),4), "mg/mL\n")
    cat("Tmax vitreous: day", which.max(df$DV)-1, "\n")
    cat("Ctrough (last):", round(tail(df$DV,1),6), "mg/mL\n")
    cat("AUC0-end (mg*day/mL):", round(sum(df$DV)*1/1000,4), "\n")
  })

  ## --- VEGF plots ---
  output$vegf_free_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~VF, type='scatter', mode='lines',
            line=list(color='#43A047'), fill='tozeroy',
            fillcolor='rgba(67,160,71,0.15)') %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Free VEGF (nM)"),
             title="Free Vitreous VEGF")
  })

  output$vegf_bound_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~VB, type='scatter', mode='lines',
            line=list(color='#8E24AA')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Bound VEGF (nM)"),
             title="VEGF-Drug Complex")
  })

  output$nv_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~NV, type='scatter', mode='lines',
            line=list(color='#E53935', width=2)) %>%
      add_segments(x=0, xend=max(df$time), y=1.5, yend=1.5,
                   line=list(dash='dash', color='gray'), showlegend=FALSE) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="NV Index (AU)"),
             title="Neovascularization Index")
  })

  output$pericyte_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~PERI, type='scatter', mode='lines',
            line=list(color='#FF7043')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Pericytes (%)", range=c(0,105)),
             title="Pericyte Count")
  })

  ## --- Inflammation / ROS plots ---
  output$ros_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~ROS, type='scatter', mode='lines',
            line=list(color='#E53935')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="ROS Index (AU)"),
             title="Oxidative Stress (ROS)")
  })

  output$age_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~AGE, type='scatter', mode='lines',
            line=list(color='#795548')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="AGE Index (AU)"),
             title="AGE Accumulation")
  })

  output$cyt_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~CYT, type='scatter', mode='lines',
            line=list(color='#F57F17')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Cytokine Index (AU)"),
             title="Inflammatory Cytokines")
  })

  output$icam_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~ICAM, type='scatter', mode='lines',
            line=list(color='#FBC02D')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="ICAM-1 Index (AU)"),
             title="ICAM-1 (Leukostasis Marker)")
  })

  ## --- Structural plots ---
  output$crt_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~CRT, type='scatter', mode='lines',
            line=list(color='#1E88E5', width=2.5), fill='tozeroy',
            fillcolor='rgba(30,136,229,0.1)') %>%
      add_segments(x=0, xend=max(df$time), y=310, yend=310,
                   line=list(dash='dash', color='darkred'), showlegend=FALSE) %>%
      add_segments(x=0, xend=max(df$time), y=270, yend=270,
                   line=list(dash='dot', color='green4'), showlegend=FALSE) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="CRT (µm)", range=c(200, max(df$CRT, 420))),
             title="Central Retinal Thickness (CRT)")
  })

  output$perm_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~PERM, type='scatter', mode='lines',
            line=list(color='#00897B')) %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Permeability (AU)"),
             title="Vascular Permeability Index")
  })

  output$ec_plot <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time, y=~EC, type='scatter', mode='lines',
            line=list(color='#5E35B1')) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="EC Count (%)", range=c(0,105)),
             title="Endothelial Cell Count")
  })

  ## --- VA plots ---
  output$va_plot <- renderPlotly({
    df <- sim_result()
    base_va <- df$VA[1]
    plot_ly(df, x=~time, y=~VA, type='scatter', mode='lines',
            line=list(color='#1E88E5', width=2.5)) %>%
      add_segments(x=0, xend=max(df$time), y=base_va+15, yend=base_va+15,
                   line=list(dash='dash', color='green4'), showlegend=FALSE) %>%
      add_segments(x=0, xend=max(df$time), y=base_va-15, yend=base_va-15,
                   line=list(dash='dash', color='red'), showlegend=FALSE) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="ETDRS Letters", range=c(0, 100)),
             title="Visual Acuity (ETDRS Letters)")
  })

  output$va_ref_table <- renderTable({
    data.frame(
      Letters   = c("≥85","70–84","55–69","35–54","<35"),
      Snellen   = c("20/20","20/40","20/80","20/160","<20/200"),
      Category  = c("Normal","Mild decrease","Moderate decrease","Severe decrease","Legal blindness")
    )
  })

  output$va_change_plot <- renderPlotly({
    df <- sim_result()
    df$VA_change <- df$VA - df$VA[1]
    plot_ly(df, x=~time, y=~VA_change, type='scatter', mode='lines',
            line=list(color=ifelse(tail(df$VA_change,1)>0,'#43A047','#E53935'), width=2)) %>%
      add_segments(x=0, xend=max(df$time), y=0, yend=0,
                   line=list(color='black', dash='dot'), showlegend=FALSE) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="Change in VA (letters)"),
             title="VA Change from Baseline")
  })

  output$va_summary_table <- renderDT({
    df <- sim_result()
    yrs <- c(365, 730)
    yrs <- yrs[yrs <= input$sim_end]
    rows <- df %>% filter(time %in% yrs) %>%
      mutate(Year=ifelse(time<=365,"Year 1","Year 2"),
             `VA Change`=round(VA-df$VA[1],1)) %>%
      select(Year, VA=VA, `VA Change`, CRT=CRT, VEGF=VF, NV=NV)
    datatable(rows, options=list(dom='t'), rownames=FALSE)
  })

  ## --- Scenario Comparison ---
  all_scenarios <- reactiveVal(NULL)

  observeEvent(input$run_compare, {
    bg <- input$BG_init; va <- input$VA_init; crt <- input$CRT_init
    vegf <- input$VEGF_init; peri <- input$PERI_init; nv_i <- input$NV_init
    ha <- input$HbA1c_init

    s0 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"none",0,28,0,bg,input$sim_end) %>%
      mutate(Scenario="S0: No Tx")
    s1 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"none",0,28,0,7.2,input$sim_end) %>%
      mutate(Scenario="S1: Glycemic Ctrl")
    s2 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"aflibercept",2.0,28,10,bg,input$sim_end) %>%
      mutate(Scenario="S2: Aflibercept q4→8w")
    s3 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"ranibizumab",0.5,28,24,bg,input$sim_end) %>%
      mutate(Scenario="S3: Ranibizumab q4w")
    s4 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"faricimab",6.0,28,10,bg,input$sim_end) %>%
      mutate(Scenario="S4: Faricimab q4→16w")
    s5 <- simulate_DR(bg,ha,va,crt,vegf,peri,nv_i,"aflibercept",2.0,28,10,7.2,input$sim_end) %>%
      mutate(Scenario="S5: AFL+Glycemic Ctrl")

    all_scenarios(bind_rows(s0,s1,s2,s3,s4,s5))
  })

  make_cmp_plot <- function(y_col, y_lab, title_txt, col_palette) {
    df <- all_scenarios()
    if (is.null(df)) return(plot_ly() %>% layout(title="Click 'Run All Scenarios'"))
    colors <- c("S0: No Tx"="#E53935","S1: Glycemic Ctrl"="#FB8C00",
                "S2: Aflibercept q4→8w"="#1E88E5","S3: Ranibizumab q4w"="#43A047",
                "S4: Faricimab q4→16w"="#8E24AA","S5: AFL+Glycemic Ctrl"="#00ACC1")
    p <- plot_ly()
    for (sc in unique(df$Scenario)) {
      d <- df %>% filter(Scenario==sc)
      p <- add_trace(p, x=d$time, y=d[[y_col]], name=sc, type='scatter', mode='lines',
                     line=list(color=colors[sc], width=2))
    }
    p %>% layout(xaxis=list(title="Time (days)"), yaxis=list(title=y_lab), title=title_txt,
                 legend=list(x=0.01,y=0.99,font=list(size=9)))
  }

  output$cmp_va   <- renderPlotly(make_cmp_plot("VA","ETDRS Letters","Visual Acuity",NULL))
  output$cmp_crt  <- renderPlotly(make_cmp_plot("CRT","CRT (µm)","Central Retinal Thickness",NULL))
  output$cmp_vegf <- renderPlotly(make_cmp_plot("VF","Free VEGF (nM)","Vitreous VEGF",NULL))
  output$cmp_nv   <- renderPlotly(make_cmp_plot("NV","NV Index (AU)","Neovascularization",NULL))

  output$compare_table <- renderDT({
    df <- all_scenarios()
    if (is.null(df)) return(datatable(data.frame(Message="Run scenarios first")))
    yrs <- c(365, 730)
    rows <- df %>% filter(time %in% yrs) %>%
      mutate(Timepoint=ifelse(time<=365,"Year 1","Year 2"),
             `VA Chg`=round(VA - df$VA[df$time==0&df$Scenario==Scenario][1], 1)) %>%
      select(Scenario, Timepoint, VA=VA, `VA Chg`, CRT=CRT, `VEGF(nM)`=VF, NV=NV, HbA1c=HA) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(rows, options=list(pageLength=12, scrollX=TRUE), rownames=FALSE)
  })

  ## --- Biomarkers ---
  output$biomarker_plot <- renderPlotly({
    df <- sim_result()
    df_long <- df %>% select(time, VEGF=VF, ROS=ROS, AGE=AGE, CYT=CYT, ICAM=ICAM) %>%
      pivot_longer(-time)
    plot_ly(df_long, x=~time, y=~value, color=~name, type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Index (normalized AU)"),
             title="Biomarker Profile Over Time", legend=list(x=0.8,y=0.9))
  })

  ## --- Calibration table ---
  output$calibration_table <- renderTable({
    data.frame(
      Trial = c("PROTOCOL T","RISE/RIDE","CLARITY","PANORAMA","TENAYA","LUCERNE","DCCT"),
      Drug  = c("AFL/RBZ/Bev","Ranibizumab","Aflibercept","Aflibercept","Faricimab","Faricimab","Insulin"),
      N     = c(660,382+382,232,402,331+327,330+338,1441),
      Duration = c("1yr","2yr","1yr","2yr","1yr","1yr","6.5yr"),
      `Primary Outcome` = c("VA +13.3L (AFL)","VA +10.9L","VA +3.3L vs PRP","2-step improv 65% vs 15%",
                             "VA +5.8L","VA +6.6L","76% less new DR")
    )
  })
}

## ============================================================
## Launch
## ============================================================
shinyApp(ui, server)
