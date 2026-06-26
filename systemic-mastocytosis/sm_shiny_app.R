################################################################################
# Systemic Mastocytosis QSP – Shiny Interactive Dashboard
# 8 Tabs: Patient Profile · Drug PK · BM MC Dynamics · Serum Tryptase ·
#         Clinical Endpoints · Scenario Comparison · Bone Disease · Biomarkers
################################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ── Embedded ODE solver (RK4) – no mrgsolve needed in Shiny ──────────────────
sm_ode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    Cm <- CENT_M / Vc_M
    Ca <- CENT_A / Vc_A

    inh_M <- Cm^gamma_M / (IC50_M^gamma_M + Cm^gamma_M)
    inh_A <- Ca^gamma_A / (IC50_A^gamma_A + Ca^gamma_A)
    KIT_inh <- 1 - (1 - inh_M) * (1 - inh_A)
    KIT_sig  <- KIT_stim * (1 - KIT_inh)

    dGUT_M  <- -ka_M * GUT_M
    dCENT_M <- F_M * ka_M * GUT_M - (CL_M + Q_M)/Vc_M * CENT_M

    dGUT_A  <- -ka_A * GUT_A
    dCENT_A <- F_A * ka_A * GUT_A - (CL_A + Q_A)/Vc_A * CENT_A + Q_A/Vp_A * PERI_A
    dPERI_A <- Q_A/Vc_A * CENT_A - Q_A/Vp_A * PERI_A

    MCP_prod  <- k_prod * (1 + KIT_sig)
    MCP_diff  <- k_diff * MCP
    MCP_apop  <- k_death_p * MCP
    dMCP <- MCP_prod - MCP_diff - MCP_apop

    BM_prol  <- k_prol * MC_BM * (1 + KIT_sig) * (1 - MC_BM/Kmax)
    BM_death <- k_death_BM * MC_BM
    BM_egress <- k_egress * MC_BM
    dMC_BM <- MCP_diff + BM_prol - BM_death - BM_egress

    dMC_SK <- k_SK_in * MC_BM - k_SK_death * MC_SK
    dMC_VS <- k_VS_in * MC_BM - k_VS_death * MC_VS

    dTRYP <- k_tryprel * MC_BM - k_tryp_el * TRYP
    dHIST <- k_histrel * (MC_SK + MC_VS) * k_hist_act - k_hist_el * HIST
    dPGD2 <- k_PGD2rel * (MC_SK + MC_VS) - k_PGD2_el * PGD2

    bres  <- (k_bres + IL6_BMD) * MC_VS
    bform <- k_bform
    dBMD  <- bform - bres

    dSYM  <- k_sym_hist * HIST + k_sym_PGD2 * PGD2 + k_sym_base * (SYM0 - SYM)

    dSPLV <- k_splv * MC_VS - k_splv_norm * (SPLV - SPLV0)
    dHEMO <- k_hemo_prod - k_hemo_loss * MC_BM

    list(c(dGUT_M, dCENT_M, dGUT_A, dCENT_A, dPERI_A,
           dMCP, dMC_BM, dMC_SK, dMC_VS,
           dTRYP, dHIST, dPGD2, dBMD, dSYM, dSPLV, dHEMO))
  })
}

run_sm_sim <- function(
    mido_dose = 0, ava_dose = 0, sim_wk = 24,
    MCP0 = 100, MC_BM0 = 500, MC_SK0 = 200, MC_VS0 = 300,
    TRYP0 = 80, SYM0 = 55, BMD0 = 1.10, SPLV0 = 600, HEMO0 = 11
) {
  params <- list(
    ka_M=0.50, F_M=0.85, CL_M=28, Vc_M=180, Q_M=12, Vp_M=400,
    IC50_M=268, gamma_M=1.2,
    ka_A=0.40, F_A=0.73, CL_A=8, Vc_A=250, Q_A=25, Vp_A=1310,
    IC50_A=0.094, gamma_A=1.5,
    k_prod=0.010, k_diff=0.005, k_death_p=0.003, KIT_stim=3.0,
    k_prol=0.008, k_death_BM=0.006, k_egress=0.003, Kmax=5000,
    k_SK_in=0.002, k_SK_death=0.004,
    k_VS_in=0.003, k_VS_death=0.004,
    k_tryprel=0.0050, k_tryp_el=0.040,
    k_histrel=0.0002, k_hist_el=0.50, k_hist_act=0.005,
    k_PGD2rel=0.0010, k_PGD2_el=0.30,
    k_bres=0.000020, k_bform=0.000018, IL6_BMD=0.000005,
    k_sym_hist=0.008, k_sym_PGD2=0.003, k_sym_base=0.002,
    SYM0=SYM0, SPLV0=SPLV0,
    k_splv=0.00005, k_splv_norm=0.00003,
    k_hemo_loss=0.00002, k_hemo_prod=0.00005
  )
  times <- seq(0, sim_wk * 168, by = 2)  # hours
  dt    <- 2

  state <- c(GUT_M=0, CENT_M=0, GUT_A=0, CENT_A=0, PERI_A=0,
             MCP=MCP0, MC_BM=MC_BM0, MC_SK=MC_SK0, MC_VS=MC_VS0,
             TRYP=TRYP0, HIST=8, PGD2=120, BMD=BMD0,
             SYM=SYM0, SPLV=SPLV0, HEMO=HEMO0)

  results <- list()
  dose_M_ng <- mido_dose * 1e6
  dose_A_ng <- ava_dose  * 1e6

  for (i in seq_along(times)) {
    t <- times[i]
    # Dosing events
    if (mido_dose > 0 && t > 0 && (t %% 12) < dt) state["GUT_M"] <- state["GUT_M"] + dose_M_ng
    if (ava_dose  > 0 && t > 0 && (t %% 24) < dt) state["GUT_A"] <- state["GUT_A"] + dose_A_ng

    # RK4 integration
    k1 <- unlist(sm_ode(t,      state,         params))
    k2 <- unlist(sm_ode(t+dt/2, state+dt/2*k1, params))
    k3 <- unlist(sm_ode(t+dt/2, state+dt/2*k2, params))
    k4 <- unlist(sm_ode(t+dt,   state+dt*k3,   params))
    state <- state + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    state[state < 0] <- 0

    inh_M <- (state["CENT_M"]/params$Vc_M)^params$gamma_M /
      (params$IC50_M^params$gamma_M + (state["CENT_M"]/params$Vc_M)^params$gamma_M)
    inh_A <- (state["CENT_A"]/params$Vc_A)^params$gamma_A /
      (params$IC50_A^params$gamma_A + (state["CENT_A"]/params$Vc_A)^params$gamma_A)

    results[[i]] <- data.frame(
      time_h   = t,
      time_wk  = t / 168,
      Cm       = state["CENT_M"] / params$Vc_M,
      Ca       = state["CENT_A"] / params$Vc_A,
      KIT_inh  = (1 - (1-inh_M)*(1-inh_A)) * 100,
      MCP      = state["MCP"],
      MC_BM    = state["MC_BM"],
      MC_SK    = state["MC_SK"],
      MC_VS    = state["MC_VS"],
      TRYP     = state["TRYP"],
      HIST     = state["HIST"],
      PGD2     = state["PGD2"],
      BMD      = state["BMD"],
      SYM      = state["SYM"],
      SPLV     = state["SPLV"],
      HEMO     = state["HEMO"]
    )
  }
  bind_rows(results)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Systemic Mastocytosis QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",    icon = icon("user")),
      menuItem("Drug PK",             tabName = "pk",         icon = icon("pills")),
      menuItem("BM MC Dynamics",      tabName = "mc",         icon = icon("dna")),
      menuItem("Serum Tryptase",      tabName = "tryptase",   icon = icon("vial")),
      menuItem("Clinical Endpoints",  tabName = "endpoints",  icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "comparison", icon = icon("balance-scale")),
      menuItem("Bone Disease",        tabName = "bone",       icon = icon("bone")),
      menuItem("Biomarker Panel",     tabName = "biomarkers", icon = icon("microscope"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f4f4; }
      .box { border-top: 3px solid #6A1B9A; }
    "))),
    tabItems(

      # ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem("profile",
        fluidRow(
          box(title = "Patient Parameters", width = 4, status = "primary",
            selectInput("sm_subtype", "SM Subtype:",
              choices = c("Indolent SM (ISM)" = "ISM",
                          "Smoldering SM (SSM)" = "SSM",
                          "Aggressive SM (ASM)" = "ASM",
                          "SM-AHN" = "SM_AHN",
                          "MC Leukemia (MCL)" = "MCL"),
              selected = "ASM"),
            sliderInput("age",      "Age (years):", 20, 80, 52),
            sliderInput("weight",   "Weight (kg):", 40, 120, 72),
            radioButtons("kit_status", "KIT Mutation:",
              choices = c("D816V Positive" = "d816v", "Wild-type" = "wt"),
              selected = "d816v"),
            numericInput("tryp_base", "Baseline Tryptase (ng/mL):", 85, 20, 500),
            numericInput("bm_mc_pct", "BM MC % (biopsy):", 40, 5, 95),
            checkboxGroupInput("c_findings", "C-Findings:",
              choices = c("Cytopenias" = "cyto",
                          "Hepatomegaly" = "hep",
                          "Splenomegaly" = "spl",
                          "GI malabsorption" = "gi",
                          "Osteolysis" = "os"))
          ),
          box(title = "Disease Summary", width = 8, status = "info",
            fluidRow(
              valueBoxOutput("vb_subtype", width = 6),
              valueBoxOutput("vb_trypt",   width = 6)
            ),
            fluidRow(
              valueBoxOutput("vb_bm",    width = 4),
              valueBoxOutput("vb_crit",  width = 4),
              valueBoxOutput("vb_risk",  width = 4)
            ),
            h4("WHO 2022 Diagnostic Criteria"),
            tableOutput("diag_criteria")
          )
        )
      ),

      # ── Tab 2: Drug PK ─────────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Drug Selection & Dose", width = 4, status = "primary",
            sliderInput("mido_dose", "Midostaurin dose (mg BID):", 0, 200, 100, step = 25),
            sliderInput("ava_dose",  "Avapritinib dose (mg QD):", 0, 300, 200, step = 25),
            sliderInput("sim_wk",    "Simulation duration (weeks):", 4, 52, 24),
            actionButton("run_sim", "Run Simulation", class = "btn-primary btn-lg",
                         icon = icon("play"), width = "100%")
          ),
          box(title = "PK Concentration–Time Profile", width = 8, status = "info",
            plotOutput("pk_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "KIT D816V Inhibition over Time", width = 6,
            plotOutput("kit_inh_plot", height = "300px")),
          box(title = "PK Parameters Summary", width = 6,
            tableOutput("pk_params"))
        )
      ),

      # ── Tab 3: BM MC Dynamics ──────────────────────────────────────────────
      tabItem("mc",
        fluidRow(
          box(title = "BM Mast Cell Burden", width = 8, status = "warning",
            plotOutput("mc_bm_plot", height = "380px")
          ),
          box(title = "MC Compartment Breakdown", width = 4,
            plotOutput("mc_comp_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "MC Progenitor (MCP) Dynamics", width = 6,
            plotOutput("mcp_plot", height = "280px")),
          box(title = "Tissue MC Distribution", width = 6,
            plotOutput("mc_dist_plot", height = "280px"))
        )
      ),

      # ── Tab 4: Serum Tryptase ──────────────────────────────────────────────
      tabItem("tryptase",
        fluidRow(
          box(title = "Serum Tryptase over Time", width = 8, status = "danger",
            plotOutput("trypt_plot", height = "380px")
          ),
          box(title = "Tryptase Response Metrics", width = 4,
            valueBoxOutput("vb_trypt_red", width = 12),
            valueBoxOutput("vb_trypt_norm", width = 12),
            plotOutput("trypt_bar", height = "200px")
          )
        ),
        fluidRow(
          box(title = "Tryptase Threshold Analysis", width = 12,
            p("Serum tryptase ≥20 ng/mL is a WHO minor criterion for SM.",
              style = "color:grey"),
            plotOutput("trypt_thresh", height = "250px"))
        )
      ),

      # ── Tab 5: Clinical Endpoints ──────────────────────────────────────────
      tabItem("endpoints",
        fluidRow(
          box(title = "Symptom Score (MISS 0-100)", width = 6,
            plotOutput("sym_plot", height = "320px")),
          box(title = "Spleen Volume (cm³)", width = 6,
            plotOutput("splv_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Hemoglobin (g/dL)", width = 6,
            plotOutput("hemo_plot", height = "280px")),
          box(title = "Clinical Response Table (Week 24)", width = 6,
            tableOutput("response_table"))
        )
      ),

      # ── Tab 6: Scenario Comparison ─────────────────────────────────────────
      tabItem("comparison",
        fluidRow(
          box(title = "Treatment Scenario Comparison", width = 12, status = "success",
            fluidRow(
              column(4,
                checkboxGroupInput("scen_select", "Select Scenarios:",
                  choices = c(
                    "Untreated"               = "none",
                    "Midostaurin 100mg BID"   = "mido100",
                    "Avapritinib 200mg QD"    = "ava200",
                    "Avapritinib 25mg QD"     = "ava25",
                    "Cladribine 3×Q4W"        = "clad"
                  ),
                  selected = c("none", "mido100", "ava200", "ava25")
                )
              ),
              column(8,
                plotOutput("compare_plot", height = "400px")
              )
            )
          )
        ),
        fluidRow(
          box(title = "Endpoint Comparison Table (Week 24)", width = 12,
            DTOutput("compare_table"))
        )
      ),

      # ── Tab 7: Bone Disease ────────────────────────────────────────────────
      tabItem("bone",
        fluidRow(
          box(title = "Bone Mineral Density (BMD) over Time", width = 7, status = "warning",
            plotOutput("bmd_plot", height = "380px")
          ),
          box(title = "Bone Disease Risk Factors", width = 5,
            sliderInput("bmd_base_in", "Baseline BMD (g/cm²):", 0.7, 1.4, 1.10, step = 0.05),
            sliderInput("mc_vs_scale", "Visceral MC severity (×):", 0.5, 3, 1, step = 0.1),
            actionButton("update_bone", "Update Bone Model", class = "btn-warning"),
            hr(),
            h5("T-score interpretation:"),
            tags$ul(
              tags$li("T-score ≥ −1.0 → Normal"),
              tags$li("−2.5 < T-score < −1.0 → Osteopenia"),
              tags$li("T-score ≤ −2.5 → Osteoporosis")
            ),
            plotOutput("tscore_plot", height = "200px")
          )
        ),
        fluidRow(
          box(title = "RANKL/OPG Balance", width = 6,
            plotOutput("rankl_plot", height = "250px")),
          box(title = "Bisphosphonate / Denosumab Effect", width = 6,
            radioButtons("bone_tx", "Bone Treatment:",
              c("None", "Zoledronate", "Denosumab"), "None"),
            plotOutput("bone_tx_plot", height = "250px"))
        )
      ),

      # ── Tab 8: Biomarker Panel ──────────────────────────────────────────────
      tabItem("biomarkers",
        fluidRow(
          box(title = "Mediator Panel over Time", width = 8, status = "info",
            plotOutput("biomarker_panel", height = "420px")
          ),
          box(title = "Biomarker Correlation Matrix", width = 4,
            plotOutput("corr_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Summary Biomarker Table", width = 12,
            DTOutput("biomarker_table"))
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run simulation on button click (or initial load)
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", {
      run_sm_sim(
        mido_dose = input$mido_dose,
        ava_dose  = input$ava_dose,
        sim_wk    = input$sim_wk,
        TRYP0     = input$tryp_base,
        MC_BM0    = input$bm_mc_pct * 50,  # scale
        SYM0      = ifelse(input$sm_subtype %in% c("ASM","MCL"), 65, 50)
      )
    })
  }, ignoreNULL = FALSE)

  # Default sim (auto-run)
  observe({
    isolate({ if (!is.null(sim_data())) return() })
  })

  # ── Tab 1 outputs ──────────────────────────────────────────────────────────
  output$vb_subtype <- renderValueBox({
    cols <- c(ISM="#2196F3", SSM="#FF9800", ASM="#F44336", SM_AHN="#9C27B0", MCL="#000000")
    valueBox(input$sm_subtype, "Disease Subtype",
             icon = icon("disease"), color = "purple")
  })
  output$vb_trypt <- renderValueBox({
    col <- if (input$tryp_base > 200) "red" else if (input$tryp_base > 20) "orange" else "green"
    valueBox(paste0(input$tryp_base, " ng/mL"), "Baseline Tryptase",
             icon = icon("vial"), color = col)
  })
  output$vb_bm <- renderValueBox({
    valueBox(paste0(input$bm_mc_pct, "%"), "BM MC%",
             icon = icon("dna"), color = "purple")
  })
  output$vb_crit <- renderValueBox({
    n_c <- length(input$c_findings)
    valueBox(n_c, "C-Findings", icon = icon("exclamation-triangle"),
             color = if (n_c > 0) "red" else "green")
  })
  output$vb_risk <- renderValueBox({
    risk <- switch(input$sm_subtype,
      ISM = "Low", SSM = "Intermediate", ASM = "High", SM_AHN = "Very High", MCL = "Critical")
    valueBox(risk, "Risk Category", icon = icon("chart-bar"), color = "navy")
  })
  output$diag_criteria <- renderTable({
    tibble(
      Criterion = c("Major: ≥25% atypical MC in BM aggregates",
                    "Minor 1: <25% MC spindle-shaped morphology",
                    "Minor 2: KIT D816V mutation",
                    "Minor 3: CD25 ± CD2 co-expression",
                    "Minor 4: Tryptase >20 ng/mL"),
      Status = c(
        if (input$bm_mc_pct >= 25) "✅ Positive" else "❌ Negative",
        if (input$bm_mc_pct >= 25) "✅ Positive" else "❓ Unknown",
        if (input$kit_status == "d816v") "✅ Positive" else "❌ Negative",
        "✅ Positive (assumed)",
        if (input$tryp_base > 20) "✅ Positive" else "❌ Negative"
      )
    )
  })

  # ── Tab 2 outputs ──────────────────────────────────────────────────────────
  output$pk_plot <- renderPlot({
    d <- sim_data()
    req(d)
    d_sub <- d %>% filter(time_wk <= min(4, max(time_wk)))
    ggplot(d_sub) +
      geom_line(aes(time_wk, Cm, color = "Midostaurin"), size = 1.1) +
      geom_line(aes(time_wk, Ca * 1000, color = "Avapritinib ×1000"), size = 1.1) +
      scale_color_manual(values = c("Midostaurin"="#1565C0","Avapritinib ×1000"="#00695C")) +
      labs(title = "Drug Concentrations (first 4 weeks)",
           x = "Time (weeks)", y = "Concentration (ng/mL)", color = "Drug") +
      theme_bw(base_size = 13)
  })

  output$kit_inh_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, KIT_inh)) +
      geom_line(color = "#880E4F", size = 1.2) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "red", alpha = 0.7) +
      labs(title = "KIT D816V Inhibition (%)",
           x = "Time (weeks)", y = "KIT Inhibition (%)") +
      theme_bw(base_size = 13)
  })

  output$pk_params <- renderTable({
    tibble(
      Parameter = c("Midostaurin CL (L/h)", "Midostaurin Vc (L)", "Midostaurin t½ (h)",
                    "Avapritinib CL (L/h)", "Avapritinib Vc (L)", "Avapritinib t½ (h)",
                    "IC50_Mido vs KIT D816V (ng/mL)", "IC50_Ava vs KIT D816V (ng/mL)"),
      Value = c("28", "180", "~45", "8", "250", "~32", "268", "0.094")
    )
  })

  # ── Tab 3 outputs ──────────────────────────────────────────────────────────
  output$mc_bm_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, MC_BM)) +
      geom_line(color = "#4A148C", size = 1.3) +
      geom_hline(yintercept = 500, linetype = "dashed", color = "grey50") +
      labs(title = "Bone Marrow Mast Cell Burden over Time",
           x = "Time (weeks)", y = "MC_BM (AU)", subtitle = "Dashed = baseline") +
      theme_bw(base_size = 13)
  })

  output$mc_comp_plot <- renderPlot({
    d <- sim_data(); req(d)
    d_long <- d %>%
      select(time_wk, MCP, MC_BM, MC_SK, MC_VS) %>%
      pivot_longer(-time_wk, names_to = "Compartment", values_to = "Count")
    ggplot(d_long, aes(time_wk, Count, color = Compartment)) +
      geom_line(size = 1) +
      scale_color_manual(values = c(MCP="#9C27B0",MC_BM="#3F51B5",
                                    MC_SK="#E91E63",MC_VS="#FF5722")) +
      labs(title = "MC by Compartment", x = "Time (weeks)", y = "MC (AU)") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$mcp_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, MCP)) +
      geom_line(color = "#7B1FA2", size = 1.2) +
      labs(x = "Time (weeks)", y = "MCP (AU)", title = "Mast Cell Progenitors") +
      theme_bw(base_size = 12)
  })

  output$mc_dist_plot <- renderPlot({
    d_end <- sim_data() %>% filter(time_wk == max(time_wk)) %>% slice(1); req(d_end)
    df <- tibble(
      Compartment = c("Bone Marrow", "Skin", "Visceral"),
      MC_Count    = c(d_end$MC_BM, d_end$MC_SK, d_end$MC_VS)
    )
    ggplot(df, aes(Compartment, MC_Count, fill = Compartment)) +
      geom_col() +
      scale_fill_manual(values = c("#3F51B5","#E91E63","#FF5722")) +
      labs(title = "MC Distribution at Simulation End",
           x = "", y = "MC (AU)") +
      theme_bw(base_size = 12) + theme(legend.position = "none")
  })

  # ── Tab 4 outputs ──────────────────────────────────────────────────────────
  output$trypt_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, TRYP)) +
      geom_line(color = "#C62828", size = 1.3) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 200, linetype = "dotted", color = "red") +
      labs(title = "Serum Tryptase over Time",
           subtitle = "Orange: normal limit (20 ng/mL)  Red: advanced SM threshold (200 ng/mL)",
           x = "Time (weeks)", y = "Tryptase (ng/mL)") +
      theme_bw(base_size = 13)
  })

  output$vb_trypt_red <- renderValueBox({
    d <- sim_data(); req(d)
    t_end <- d %>% filter(time_wk == max(time_wk)) %>% pull(TRYP) %>% mean()
    t0    <- d %>% filter(time_wk == 0) %>% pull(TRYP) %>% mean()
    pct   <- round((1 - t_end/t0)*100, 1)
    valueBox(paste0(pct, "%"), "Tryptase Reduction",
             icon = icon("arrow-down"), color = if (pct > 50) "green" else "orange")
  })

  output$vb_trypt_norm <- renderValueBox({
    d <- sim_data(); req(d)
    t_end <- d %>% filter(time_wk == max(time_wk)) %>% pull(TRYP) %>% mean()
    valueBox(round(t_end, 1), "Final Tryptase (ng/mL)",
             icon = icon("vial"), color = if (t_end < 20) "green" else "red")
  })

  output$trypt_bar <- renderPlot({
    d <- sim_data(); req(d)
    pts <- d %>% filter(time_wk %in% c(4, 12, 24)) %>%
      group_by(time_wk) %>% summarise(TRYP = mean(TRYP))
    ggplot(pts, aes(factor(time_wk), TRYP, fill = factor(time_wk))) +
      geom_col() + scale_fill_brewer(palette = "Reds") +
      labs(x = "Week", y = "Tryptase", title = "Landmark Tryptase") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
  })

  output$trypt_thresh <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, TRYP)) +
      geom_line(color = "#C62828") +
      geom_ribbon(aes(ymin = 0, ymax = pmin(TRYP, 20)), fill = "green", alpha = 0.2) +
      geom_ribbon(aes(ymin = pmax(TRYP, 20), ymax = TRYP), fill = "red", alpha = 0.2) +
      geom_hline(yintercept = 20, color = "orange") +
      labs(x = "Time (weeks)", y = "Tryptase (ng/mL)",
           title = "Tryptase vs. Diagnostic Threshold") +
      theme_bw(base_size = 12)
  })

  # ── Tab 5 outputs ──────────────────────────────────────────────────────────
  output$sym_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, SYM)) +
      geom_line(color = "#E65100", size = 1.3) +
      ylim(0, 100) +
      labs(title = "Symptom Score (MISS 0-100)",
           x = "Time (weeks)", y = "Symptom Score") +
      theme_bw(base_size = 13)
  })

  output$splv_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, SPLV)) +
      geom_line(color = "#1B5E20", size = 1.3) +
      geom_hline(yintercept = 500, linetype = "dashed", color = "grey50") +
      labs(title = "Spleen Volume (cm³)",
           x = "Time (weeks)", y = "Spleen Vol (cm³)") +
      theme_bw(base_size = 13)
  })

  output$hemo_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, HEMO)) +
      geom_line(color = "#B71C1C", size = 1.2) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "orange") +
      labs(title = "Hemoglobin (g/dL)", subtitle = "Orange: anemia threshold",
           x = "Time (weeks)", y = "Hgb (g/dL)") +
      theme_bw(base_size = 13)
  })

  output$response_table <- renderTable({
    d <- sim_data(); req(d)
    t0  <- d %>% filter(time_wk == 0)  %>% summarise(across(c(TRYP,SYM,SPLV,BMD,HEMO), mean))
    t24 <- d %>% filter(time_wk >= min(d$time_wk[d$time_wk >= input$sim_wk - 0.1])) %>%
                 slice(1) %>% summarise(across(c(TRYP,SYM,SPLV,BMD,HEMO), mean))
    tibble(
      Endpoint     = c("Serum Tryptase", "Symptom Score", "Spleen Vol", "BMD", "Hemoglobin"),
      Baseline     = round(as.numeric(t0),  2),
      `End of Sim` = round(as.numeric(t24), 2),
      Change_pct   = round((as.numeric(t24)/as.numeric(t0)-1)*100, 1)
    )
  })

  # ── Tab 6 outputs ──────────────────────────────────────────────────────────
  compare_sims <- reactive({
    scens <- input$scen_select
    params <- list(
      none    = list(mido=0,   ava=0),
      mido100 = list(mido=100, ava=0),
      ava200  = list(mido=0,   ava=200),
      ava25   = list(mido=0,   ava=25),
      clad    = list(mido=0,   ava=0)
    )
    bind_rows(lapply(scens, function(s) {
      p <- params[[s]]
      run_sm_sim(mido_dose=p$mido, ava_dose=p$ava, sim_wk=input$sim_wk,
                 TRYP0=input$tryp_base) %>% mutate(scenario=s)
    }))
  })

  output$compare_plot <- renderPlot({
    d <- compare_sims(); req(d)
    colors <- c(none="#616161",mido100="#1565C0",ava200="#00695C",ava25="#2E7D32",clad="#BF360C")
    ggplot(d, aes(time_wk, TRYP, color = scenario)) +
      geom_line(size = 1.1) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "orange") +
      scale_color_manual(values = colors,
        labels = c(none="Untreated",mido100="Midostaurin",ava200="Ava 200mg",
                   ava25="Ava 25mg",clad="Cladribine")) +
      labs(title = "Serum Tryptase Comparison", x = "Time (weeks)",
           y = "Tryptase (ng/mL)", color = "Scenario") +
      theme_bw(base_size = 13) + theme(legend.position = "bottom")
  })

  output$compare_table <- renderDT({
    d <- compare_sims(); req(d)
    t0 <- d %>% group_by(scenario) %>% filter(time_wk == 0) %>%
          summarise(TRYP0=mean(TRYP), SYM0=mean(SYM), .groups="drop")
    d %>% filter(time_wk >= max(d$time_wk) - 0.5) %>%
      group_by(scenario) %>%
      summarise(Tryptase_end = round(mean(TRYP),1),
                Symptom_end  = round(mean(SYM),1),
                BM_MC_end    = round(mean(MC_BM),0),
                BMD_end      = round(mean(BMD),3),
                .groups      = "drop") %>%
      datatable(options = list(pageLength = 10))
  })

  # ── Tab 7 outputs ──────────────────────────────────────────────────────────
  output$bmd_plot <- renderPlot({
    d <- sim_data(); req(d)
    ggplot(d, aes(time_wk, BMD)) +
      geom_line(color = "#827717", size = 1.3) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "orange",
                 size = 0.8) +
      geom_hline(yintercept = 0.9, linetype = "dotted", color = "red") +
      labs(title = "Bone Mineral Density over Time",
           subtitle = "Orange: osteopenia threshold (~T-score −1) | Red: osteoporosis",
           x = "Time (weeks)", y = "BMD (g/cm²)") +
      theme_bw(base_size = 13)
  })

  output$tscore_plot <- renderPlot({
    d <- sim_data(); req(d)
    d <- d %>% mutate(Tscore = (BMD - 1.05) / 0.11)
    ggplot(d, aes(time_wk, Tscore)) +
      geom_line(color = "#827717") +
      geom_hline(yintercept = -1.0, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = -2.5, color = "red", linetype = "dashed") +
      labs(x = "Week", y = "T-score") + theme_bw(base_size = 11)
  })

  output$rankl_plot <- renderPlot({
    d <- sim_data(); req(d)
    d <- d %>% mutate(RANKL = MC_VS * 0.003, OPG = 50 * (1 - MC_VS/5000))
    d_long <- d %>% select(time_wk, RANKL, OPG) %>%
      pivot_longer(-time_wk, names_to = "Molecule", values_to = "Level")
    ggplot(d_long, aes(time_wk, Level, color = Molecule)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(RANKL="#F44336",OPG="#2196F3")) +
      labs(title = "RANKL vs OPG", x = "Week", y = "Level (AU)") +
      theme_bw(base_size = 12)
  })

  output$bone_tx_plot <- renderPlot({
    d <- sim_data(); req(d)
    reduction <- switch(input$bone_tx, Zoledronate=0.30, Denosumab=0.40, `None`=0)
    d <- d %>% mutate(BMD_tx = BMD * (1 + reduction * (1 - exp(-time_wk/12))))
    ggplot() +
      geom_line(data=d, aes(time_wk, BMD), color="grey50", linetype="dashed") +
      geom_line(data=d, aes(time_wk, BMD_tx), color="#4CAF50", size=1.2) +
      labs(title = paste("BMD with", input$bone_tx),
           x = "Week", y = "BMD (g/cm²)") + theme_bw(base_size = 12)
  })

  # ── Tab 8 outputs ──────────────────────────────────────────────────────────
  output$biomarker_panel <- renderPlot({
    d <- sim_data(); req(d)
    d_long <- d %>%
      select(time_wk, TRYP, HIST, PGD2) %>%
      pivot_longer(-time_wk, names_to = "Biomarker", values_to = "Level")
    ggplot(d_long, aes(time_wk, Level, color = Biomarker)) +
      geom_line(size = 1.1) +
      facet_wrap(~Biomarker, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c(TRYP="#C62828",HIST="#1565C0",PGD2="#00695C")) +
      labs(title = "Mast Cell Mediator Panel",
           x = "Time (weeks)", y = "Concentration") +
      theme_bw(base_size = 12) + theme(legend.position = "none")
  })

  output$corr_plot <- renderPlot({
    d <- sim_data(); req(d)
    d_sub <- d %>% select(MC_BM, TRYP, HIST, PGD2, SYM, BMD)
    if (requireNamespace("corrplot", quietly=TRUE)) {
      corrplot::corrplot(cor(d_sub), method = "color", type = "upper",
                         tl.cex = 0.9, addCoef.col = "black", number.cex = 0.7)
    } else {
      ggplot() + annotate("text", x=0.5, y=0.5,
        label="Install 'corrplot' for correlation matrix") + theme_void()
    }
  })

  output$biomarker_table <- renderDT({
    d <- sim_data(); req(d)
    key_times <- c(0, 4, 12, 24)
    d %>% filter(time_wk %in% key_times | time_wk == max(time_wk)) %>%
      group_by(Week = round(time_wk)) %>%
      summarise(
        Tryptase = round(mean(TRYP), 1),
        Histamine = round(mean(HIST), 2),
        PGD2 = round(mean(PGD2), 1),
        BM_MC = round(mean(MC_BM), 0),
        Symptom = round(mean(SYM), 1),
        BMD = round(mean(BMD), 3),
        KIT_inh_pct = round(mean(KIT_inh), 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 10))
  })
}

shinyApp(ui = ui, server = server)
