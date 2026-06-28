## ============================================================
## Raynaud's Phenomenon — Shiny Dashboard
## 7 Tabs: Patient Profile · PK · Inflammation & Vasoactive
##          Vasomotor Response · Clinical Endpoints
##          Scenario Comparison · Biomarker Explorer
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

# ── Load model ───────────────────────────────────────────────
mod_base <- mread("raynaud_mrgsolve_model.R")

# ── Helper ───────────────────────────────────────────────────
sim_raynaud <- function(params, drug_events, end_h = 2016, delta = 2) {
  mod_base %>%
    param(params) %>%
    mrgsim(ev = drug_events, end = end_h, delta = delta) %>%
    as.data.frame()
}

build_ev <- function(drugs) {
  events <- vector("list", length(drugs))
  for (i in seq_along(drugs)) {
    d <- drugs[[i]]
    if (d$use && d$dose > 0) {
      events[[i]] <- ev(cmt = d$cmt, amt = d$dose * 1000,
                        ii = d$ii, addl = d$addl, time = 0)
    }
  }
  events <- Filter(Negate(is.null), events)
  if (length(events) == 0) return(ev(time = 0, cmt = "GUT_NIF", amt = 0))
  do.call(c, events)
}

# ── UI ───────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Raynaud's Phenomenon QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Vasoactive Mediators", tabName = "tab_vasoact",   icon = icon("dna")),
      menuItem("Vasomotor Response",   tabName = "tab_vasomotor", icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",   tabName = "tab_clinical",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",   icon = icon("balance-scale")),
      menuItem("Biomarker Explorer",   tabName = "tab_biomarker", icon = icon("microscope"))
    )
  ),

  dashboardBody(
    tabItems(

      # ──────────────────────────────────────────────────
      # TAB 1: Patient Profile
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Disease Subtype", width = 4, status = "primary",
            radioButtons("subtype", "Raynaud's Subtype:",
              choices = c("Primary (idiopathic)" = "primary",
                          "Secondary (SSc-associated)" = "secondary",
                          "Secondary (SLE/MCTD)" = "secondary_sle"),
              selected = "primary"),
            hr(),
            sliderInput("alpha2_sens", "α2-AR Sensitivity:", min = 0.8, max = 2.0, value = 1.2, step = 0.05),
            sliderInput("et1_baseline", "Baseline ET-1 (pg/mL):", min = 0.3, max = 3.0, value = 0.8, step = 0.1)
          ),
          box(title = "Trigger Settings", width = 4, status = "warning",
            checkboxInput("cold_ch", "Cold Challenge Simulation", value = FALSE),
            sliderInput("cold_intensity", "Cold Stimulus Intensity:", min = 0, max = 2.0, value = 0.8, step = 0.1),
            sliderInput("stress_lvl", "Emotional Stress Level (0-2):", min = 0, max = 2.0, value = 0, step = 0.1),
            hr(),
            numericInput("sim_weeks", "Simulation Duration (weeks):", value = 12, min = 1, max = 52)
          ),
          box(title = "Comorbidities", width = 4, status = "danger",
            checkboxGroupInput("comorbid", "Select Comorbidities:",
              choices = c("Systemic sclerosis (SSc)" = "ssc",
                          "SLE" = "sle",
                          "MCTD" = "mctd",
                          "Hypothyroidism" = "hypo",
                          "Thoracic outlet syndrome" = "tos",
                          "Vibration exposure" = "vibration")),
            hr(),
            selectInput("severity", "Baseline Severity:",
              choices = c("Mild (1-3 attacks/wk)" = "mild",
                          "Moderate (4-7/wk)" = "moderate",
                          "Severe (>7/wk + DU)" = "severe")),
            valueBoxOutput("risk_box", width = 12)
          )
        ),
        fluidRow(
          box(title = "Patient Summary", width = 12, status = "info",
            verbatimTextOutput("patient_summary")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 2: Drug PK
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Selection & Dosing", width = 4, status = "primary",
            h4("Nifedipine (CCB)"),
            checkboxInput("use_NIF", "Use Nifedipine", value = FALSE),
            numericInput("dose_NIF", "Dose (mg QD):", value = 30, min = 10, max = 90),
            hr(),
            h4("Sildenafil (PDE5i)"),
            checkboxInput("use_SIL", "Use Sildenafil", value = FALSE),
            numericInput("dose_SIL", "Dose (mg BID):", value = 50, min = 25, max = 100),
            hr(),
            h4("Bosentan (ERA)"),
            checkboxInput("use_BOS", "Use Bosentan", value = FALSE),
            numericInput("dose_BOS", "Dose (mg BID):", value = 125, min = 62.5, max = 250),
            hr(),
            h4("Iloprost (IV — 5-day course)"),
            checkboxInput("use_ILO", "Use Iloprost", value = FALSE),
            numericInput("dose_ILO", "Dose (ng, 5-day course):", value = 50, min = 10, max = 200),
            hr(),
            h4("Prazosin (α1-blocker)"),
            checkboxInput("use_PRA", "Use Prazosin", value = FALSE),
            numericInput("dose_PRA", "Dose (mg BID):", value = 1, min = 0.5, max = 5),
            actionButton("run_pk", "Simulate PK", class = "btn-primary btn-block")
          ),
          box(title = "Plasma Concentration-Time Profiles", width = 8, status = "success",
            plotlyOutput("pk_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12,
            DTOutput("pk_table")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 3: Vasoactive Mediators
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_vasoact",
        fluidRow(
          box(title = "Simulation Controls", width = 3, status = "warning",
            actionButton("run_vasoact", "Run Simulation", class = "btn-warning btn-block"),
            hr(),
            checkboxInput("show_ET1", "Show ET-1", value = TRUE),
            checkboxInput("show_ROS", "Show ROS", value = TRUE),
            checkboxInput("show_NE", "Show NE Level", value = TRUE),
            checkboxInput("show_RhoA", "Show RhoA", value = TRUE)
          ),
          box(title = "Vasoactive Mediator Dynamics", width = 9, status = "warning",
            plotlyOutput("vasoact_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "Second Messenger Dynamics", width = 6, status = "info",
            plotlyOutput("second_mess_plot", height = "350px")
          ),
          box(title = "VSMC Contraction Index", width = 6, status = "danger",
            plotlyOutput("contraction_plot", height = "350px")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 4: Vasomotor Response
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_vasomotor",
        fluidRow(
          box(title = "Cold Challenge Parameters", width = 4, status = "primary",
            actionButton("run_cold", "Run Cold Challenge", class = "btn-danger btn-block"),
            hr(),
            sliderInput("cold_dur", "Cold Exposure Duration (min):", min = 5, max = 120, value = 30),
            sliderInput("cold_temp", "Temperature (°C):", min = 0, max = 15, value = 5),
            checkboxInput("pretreat_cold", "Pre-treat with Nifedipine (30mg)", value = FALSE),
            hr(),
            p("Simulates digital blood flow response to cold challenge — the key diagnostic test in Raynaud's.")
          ),
          box(title = "Digital Blood Flow — Cold Challenge", width = 8, status = "primary",
            plotlyOutput("cold_dbf_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Vasospasm Episode Frequency", width = 6, status = "warning",
            plotlyOutput("epi_freq_plot", height = "300px")
          ),
          box(title = "Attack Duration Distribution", width = 6, status = "info",
            plotlyOutput("attack_dur_plot", height = "300px")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 5: Clinical Endpoints
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Run Clinical Simulation", width = 3, status = "success",
            actionButton("run_clin", "Simulate", class = "btn-success btn-block"),
            hr(),
            p("Primary endpoints:"),
            tags$ul(
              tags$li("Raynaud Condition Score (RCS 0-10)"),
              tags$li("Vasospasm episodes/week"),
              tags$li("Digital blood flow (mL/min/100g)"),
              tags$li("VAS pain score"),
              tags$li("Digital ulcer risk (secondary)")
            )
          ),
          box(title = "Raynaud Condition Score Over Time", width = 9, status = "success",
            plotlyOutput("rcs_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Vasospasm Frequency", width = 6,
            plotlyOutput("vaso_freq_plot", height = "300px")
          ),
          box(title = "VAS Pain Score", width = 6,
            plotlyOutput("vas_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Digital Ulcer Risk (Secondary Only)", width = 6,
            plotlyOutput("du_plot", height = "300px")
          ),
          box(title = "Week 12 Endpoint Summary", width = 6,
            DTOutput("endpoint_table")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 6: Scenario Comparison
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Scenarios to Compare", width = 3, status = "primary",
            checkboxGroupInput("compare_scen", "Select Scenarios:",
              choices = c(
                "Untreated Primary"           = "s1",
                "Nifedipine 30mg QD"          = "s2",
                "Sildenafil 50mg BID"         = "s3",
                "Bosentan 125mg BID (2°)"     = "s4",
                "Iloprost IV ×5d (2°)"        = "s5",
                "Prazosin 1mg BID"            = "s6",
                "Nifedipine+Sildenafil"       = "s7",
                "Untreated Secondary (SSc)"   = "s8"
              ),
              selected = c("s1", "s2", "s3")
            ),
            actionButton("run_compare", "Run Comparison", class = "btn-primary btn-block")
          ),
          box(title = "RCS Comparison", width = 9, status = "success",
            plotlyOutput("compare_rcs_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Digital Blood Flow Comparison", width = 6,
            plotlyOutput("compare_dbf_plot", height = "300px")
          ),
          box(title = "ET-1 Comparison", width = 6,
            plotlyOutput("compare_et1_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Week 12 Treatment Comparison", width = 12,
            DTOutput("compare_table")
          )
        )
      ),

      # ──────────────────────────────────────────────────
      # TAB 7: Biomarker Explorer
      # ──────────────────────────────────────────────────
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Biomarker Selection", width = 3, status = "warning",
            selectInput("bm_x", "X-axis Biomarker:",
              choices = c("Time (weeks)" = "time_wk",
                          "Plasma NE (nmol/L)" = "NE_lvl",
                          "cGMP (nM)" = "cGMP_VSMC",
                          "cAMP (nM)" = "cAMP_VSMC",
                          "ET-1 (pg/mL)" = "ET1_plasma",
                          "ROS (AU)" = "ROS_lvl")),
            selectInput("bm_y", "Y-axis Biomarker:",
              choices = c("Digital Blood Flow" = "DBF",
                          "Vasospasm Freq/wk" = "VasoEp_wk",
                          "RCS" = "RCS",
                          "VAS Pain" = "VAS_pain",
                          "VSMC Contraction" = "VSMC_contraction",
                          "RhoA-GTP (AU)" = "RhoA_GTP")),
            actionButton("run_bm", "Plot Biomarker", class = "btn-warning btn-block")
          ),
          box(title = "Biomarker Correlation", width = 9, status = "warning",
            plotlyOutput("bm_scatter", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Pathway Mediator Dashboard", width = 12, status = "info",
            fluidRow(
              valueBoxOutput("bm_NE",    width = 2),
              valueBoxOutput("bm_ET1",   width = 2),
              valueBoxOutput("bm_cGMP",  width = 2),
              valueBoxOutput("bm_cAMP",  width = 2),
              valueBoxOutput("bm_ROS",   width = 2),
              valueBoxOutput("bm_DBF",   width = 2)
            )
          )
        ),
        fluidRow(
          box(title = "Nailfold Capillaroscopy Index (simulated)", width = 6, status = "danger",
            plotlyOutput("capillaro_plot", height = "300px")
          ),
          box(title = "Temperature Recovery Time", width = 6, status = "primary",
            plotlyOutput("temp_recov_plot", height = "300px")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# ── Server ───────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive: base parameters ──────────────────────────────
  base_params <- reactive({
    sec <- ifelse(input$subtype %in% c("secondary", "secondary_sle"), 1, 0)
    list(
      secondary      = sec,
      alpha2_sens    = input$alpha2_sens,
      ET1_base       = input$et1_baseline,
      cold_challenge = ifelse(isTRUE(input$cold_ch), 1, 0),
      k_cold_NE      = input$cold_intensity
    )
  })

  # ── Reactive: drug events ───────────────────────────────────
  drug_events_reactive <- reactive({
    drugs <- list(
      list(use = isTRUE(input$use_NIF), cmt = "GUT_NIF",
           dose = input$dose_NIF, ii = 24, addl = input$sim_weeks * 7 - 1),
      list(use = isTRUE(input$use_SIL), cmt = "GUT_SIL",
           dose = input$dose_SIL, ii = 12, addl = input$sim_weeks * 14 - 1),
      list(use = isTRUE(input$use_BOS), cmt = "GUT_BOS",
           dose = input$dose_BOS, ii = 12, addl = input$sim_weeks * 14 - 1),
      list(use = isTRUE(input$use_ILO), cmt = "PLASMA_ILO",
           dose = input$dose_ILO, ii = 24, addl = 4),
      list(use = isTRUE(input$use_PRA), cmt = "GUT_PRA",
           dose = input$dose_PRA, ii = 12, addl = input$sim_weeks * 14 - 1)
    )
    build_ev(drugs)
  })

  # ── Simulate on button click ────────────────────────────────
  sim_result <- eventReactive(c(input$run_pk, input$run_vasoact,
                                 input$run_clin, input$run_cold), {
    params <- base_params()
    evs    <- drug_events_reactive()
    sim_raynaud(params, evs, end_h = input$sim_weeks * 168, delta = 2) %>%
      mutate(time_wk = time / 168)
  }, ignoreNULL = FALSE)

  # ── TAB 1: Patient summary ────────────────────────────────
  output$patient_summary <- renderText({
    sec <- ifelse(input$subtype == "primary", "Primary (idiopathic)", "Secondary")
    paste0("Subtype: ", sec, "\n",
           "α2-AR sensitivity: ", input$alpha2_sens, "x\n",
           "Baseline ET-1: ", input$et1_baseline, " pg/mL\n",
           "Cold challenge: ", ifelse(input$cold_ch, "ON", "OFF"), "\n",
           "Simulation: ", input$sim_weeks, " weeks\n",
           "Comorbidities: ", paste(input$comorbid, collapse=", "))
  })

  output$risk_box <- renderValueBox({
    sev <- input$severity
    color <- switch(sev, mild="green", moderate="yellow", severe="red", "blue")
    lbl   <- switch(sev, mild="Mild", moderate="Moderate", severe="Severe", "Unknown")
    valueBox(lbl, "Disease Severity", icon=icon("stethoscope"), color=color)
  })

  # ── TAB 2: PK plot ────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    # Gather drug concentrations
    pk_df <- df %>%
      select(time_wk, Cp_NIF_out, Cp_SIL_out, Cp_BOS_out, Cp_ILO_out, Cp_PRA_out) %>%
      pivot_longer(-time_wk, names_to = "drug", values_to = "conc") %>%
      mutate(drug = recode(drug,
        Cp_NIF_out = "Nifedipine", Cp_SIL_out = "Sildenafil",
        Cp_BOS_out = "Bosentan",  Cp_ILO_out = "Iloprost",
        Cp_PRA_out = "Prazosin"))
    p <- ggplot(pk_df, aes(x=time_wk, y=conc, color=drug)) +
      geom_line(size=0.8) +
      labs(x="Time (weeks)", y="Plasma Conc (ng/mL)", color="Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    data.frame(
      Drug       = c("Nifedipine","Sildenafil","Bosentan","Iloprost","Prazosin"),
      Route      = c("PO","PO","PO","IV","PO"),
      `F (%)` = c(85,40,50,100,68),
      `Vd (L)`   = c(120,105,18,25,97),
      `CL (L/h)` = c(60,41,4,15,40),
      `t1/2 (h)` = c(2.4,4.0,5.4,0.4,2.9),
      Mechanism  = c("L-VGCC block","PDE5 inhibit","ETA/B block","IP-R agonist","α1 block")
    )
  }, options=list(pageLength=5))

  # ── TAB 3: Vasoactive mediators ──────────────────────────
  output$vasoact_plot <- renderPlotly({
    df <- sim_result()
    vars <- c()
    if(input$show_ET1) vars <- c(vars, "ET1_plasma")
    if(input$show_ROS) vars <- c(vars, "ROS_lvl")
    if(input$show_NE)  vars <- c(vars, "NE_lvl")
    if(input$show_RhoA) vars <- c(vars, "RhoA_GTP")
    if(length(vars)==0) vars <- "ET1_plasma"
    va_df <- df %>% select(time_wk, all_of(vars)) %>%
      pivot_longer(-time_wk, names_to="biomarker", values_to="value")
    p <- ggplot(va_df, aes(x=time_wk, y=value, color=biomarker)) +
      geom_line(size=0.9) +
      labs(x="Time (weeks)", y="Level (AU / pg/mL)", color="Mediator") +
      theme_bw()
    ggplotly(p)
  })

  output$second_mess_plot <- renderPlotly({
    df <- sim_result()
    sm_df <- df %>%
      select(time_wk, cGMP_VSMC, cAMP_VSMC) %>%
      pivot_longer(-time_wk, names_to="messenger", values_to="nM")
    p <- ggplot(sm_df, aes(x=time_wk, y=nM, color=messenger)) +
      geom_line(size=0.9) +
      labs(x="Time (weeks)", y="Concentration (nM)", color="Messenger") +
      theme_bw()
    ggplotly(p)
  })

  output$contraction_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=VSMC_contraction)) +
      geom_line(color="red", size=1.0) +
      geom_hline(yintercept=0.55, linetype="dashed", color="darkred") +
      labs(x="Time (weeks)", y="VSMC Contraction Index (AU)",
           title="Red dashed = vasospasm threshold") +
      theme_bw()
    ggplotly(p)
  })

  # ── TAB 4: Cold challenge / vasomotor ─────────────────────
  output$cold_dbf_plot <- renderPlotly({
    params_cold <- base_params()
    params_cold$cold_challenge <- 1
    params_cold$k_cold_NE <- input$cold_intensity
    ev_cold <- ev(time=0, cmt="GUT_NIF", amt=0)
    ev_pretreat <- ev(time=0, cmt="GUT_NIF", amt=30*1000, ii=24, addl=1)
    df_cold <- sim_raynaud(params_cold, ev_cold, end_h=24, delta=0.1) %>%
      mutate(Condition="Cold (no Rx)", time_h=time)
    df_results <- df_cold
    if(input$pretreat_cold) {
      df_pre <- sim_raynaud(params_cold, ev_pretreat, end_h=24, delta=0.1) %>%
        mutate(Condition="Cold + Nifedipine", time_h=time)
      df_results <- bind_rows(df_cold, df_pre)
    }
    p <- ggplot(df_results, aes(x=time_h, y=DBF, color=Condition)) +
      geom_line(size=1.0) +
      labs(x="Time (hours)", y="Digital Blood Flow (mL/min/100g)",
           title="Cold Challenge — Digital Blood Flow") +
      theme_bw()
    ggplotly(p)
  })

  output$epi_freq_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=VasoEp_wk)) +
      geom_line(color="#E91E63", size=1.0) +
      labs(x="Time (weeks)", y="Vasospasm Episodes/Week") + theme_bw()
    ggplotly(p)
  })

  output$attack_dur_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=AttackDur_min)) +
      geom_line(color="#FF9800", size=1.0) +
      labs(x="Time (weeks)", y="Attack Duration (min)") + theme_bw()
    ggplotly(p)
  })

  # ── TAB 5: Clinical endpoints ──────────────────────────────
  output$rcs_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=RCS)) +
      geom_line(color="#1565C0", size=1.2) +
      ylim(0, 10) +
      labs(x="Time (weeks)", y="Raynaud Condition Score (0-10)") + theme_bw()
    ggplotly(p)
  })

  output$vaso_freq_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=VasoEp_wk)) +
      geom_line(color="#E91E63", size=1.0) +
      labs(x="Time (weeks)", y="Episodes/week") + theme_bw()
    ggplotly(p)
  })

  output$vas_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=VAS_pain)) +
      geom_line(color="orange", size=1.0) + ylim(0,10) +
      labs(x="Time (weeks)", y="VAS Pain (0-10)") + theme_bw()
    ggplotly(p)
  })

  output$du_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time_wk, y=DU_risk*100)) +
      geom_area(fill="#FFCDD2", alpha=0.6) +
      geom_line(color="red", size=1.0) +
      labs(x="Time (weeks)", y="Digital Ulcer Risk (%)") + theme_bw()
    ggplotly(p)
  })

  output$endpoint_table <- renderDT({
    df <- sim_result()
    last_wk <- df %>% filter(time_wk >= max(time_wk) - 0.5)
    data.frame(
      Endpoint = c("RCS (0-10)", "Vasospasm Ep/wk", "DBF (mL/min/100g)",
                   "VAS Pain (0-10)", "DU Risk (%)", "ET-1 (pg/mL)"),
      Value    = c(round(mean(last_wk$RCS),2),
                   round(mean(last_wk$VasoEp_wk),1),
                   round(mean(last_wk$DBF),2),
                   round(mean(last_wk$VAS_pain),2),
                   round(mean(last_wk$DU_risk)*100,1),
                   round(mean(last_wk$ET1_plasma),2))
    )
  }, options=list(pageLength=6, dom="t"))

  # ── TAB 6: Scenario comparison ─────────────────────────────
  compare_result <- eventReactive(input$run_compare, {
    scenario_defs <- list(
      s1 = list(params=list(secondary=0), evs=ev(time=0,cmt="GUT_NIF",amt=0), label="Untreated Primary"),
      s2 = list(params=list(secondary=0), evs=ev(cmt="GUT_NIF",amt=30000,ii=24,addl=83), label="Nifedipine 30mg QD"),
      s3 = list(params=list(secondary=0), evs=ev(cmt="GUT_SIL",amt=50000,ii=12,addl=167), label="Sildenafil 50mg BID"),
      s4 = list(params=list(secondary=1), evs=ev(cmt="GUT_BOS",amt=125000,ii=12,addl=167), label="Bosentan 125mg BID (2°)"),
      s5 = list(params=list(secondary=1),
                evs=do.call(c, lapply(0:4, function(d) ev(time=d*24,cmt="PLASMA_ILO",amt=50))),
                label="Iloprost IV ×5d (2°)"),
      s6 = list(params=list(secondary=0), evs=ev(cmt="GUT_PRA",amt=1000,ii=12,addl=167), label="Prazosin 1mg BID"),
      s7 = list(params=list(secondary=0),
                evs=c(ev(cmt="GUT_NIF",amt=30000,ii=24,addl=83),
                      ev(cmt="GUT_SIL",amt=25000,ii=12,addl=167)), label="Nifedipine+Sildenafil"),
      s8 = list(params=list(secondary=1), evs=ev(time=0,cmt="GUT_NIF",amt=0), label="Untreated Secondary (SSc)")
    )
    sel <- input$compare_scen
    if(length(sel)==0) sel <- c("s1","s2","s3")
    bind_rows(lapply(sel, function(s) {
      def <- scenario_defs[[s]]
      tryCatch(
        sim_raynaud(def$params, def$evs, end_h=2016, delta=4) %>%
          mutate(scenario=def$label, time_wk=time/168),
        error=function(e) NULL
      )
    }))
  })

  output$compare_rcs_plot <- renderPlotly({
    df <- compare_result()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x=time_wk, y=RCS, color=scenario)) +
      geom_line(size=0.9) + ylim(0,10) +
      labs(x="Time (weeks)", y="RCS (0-10)", color="Scenario") + theme_bw()
    ggplotly(p)
  })

  output$compare_dbf_plot <- renderPlotly({
    df <- compare_result()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x=time_wk, y=DBF, color=scenario)) +
      geom_line(size=0.8) +
      labs(x="Time (weeks)", y="DBF (mL/min/100g)") + theme_bw()
    ggplotly(p)
  })

  output$compare_et1_plot <- renderPlotly({
    df <- compare_result()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x=time_wk, y=ET1_plasma, color=scenario)) +
      geom_line(size=0.8) +
      labs(x="Time (weeks)", y="ET-1 (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    df <- compare_result()
    req(nrow(df) > 0)
    df %>% filter(time_wk >= 11) %>%
      group_by(scenario) %>%
      summarise(
        `RCS (wk12)` = round(mean(RCS),2),
        `VasoEp/wk`  = round(mean(VasoEp_wk),1),
        `DBF`        = round(mean(DBF),2),
        `ET-1 pg/mL` = round(mean(ET1_plasma),2),
        `VAS Pain`   = round(mean(VAS_pain),2)
      )
  }, options=list(pageLength=10, dom="t"))

  # ── TAB 7: Biomarker explorer ──────────────────────────────
  output$bm_scatter <- renderPlotly({
    df <- sim_result()
    req(input$bm_x %in% names(df), input$bm_y %in% names(df))
    p <- ggplot(df, aes_string(x=input$bm_x, y=input$bm_y)) +
      geom_point(alpha=0.4, color="#1565C0", size=0.8) +
      geom_smooth(method="loess", se=FALSE, color="red") +
      labs(x=input$bm_x, y=input$bm_y) + theme_bw()
    ggplotly(p)
  })

  last_val <- reactive({
    df <- sim_result()
    df[nrow(df), ]
  })

  output$bm_NE   <- renderValueBox(valueBox(round(last_val()$NE_lvl,3), "NE (nmol/L)", icon=icon("bolt"), color="orange"))
  output$bm_ET1  <- renderValueBox(valueBox(round(last_val()$ET1_plasma,2), "ET-1 (pg/mL)", icon=icon("wind"), color="red"))
  output$bm_cGMP <- renderValueBox(valueBox(round(last_val()$cGMP_VSMC,2), "cGMP (nM)", icon=icon("circle"), color="green"))
  output$bm_cAMP <- renderValueBox(valueBox(round(last_val()$cAMP_VSMC,2), "cAMP (nM)", icon=icon("circle"), color="teal"))
  output$bm_ROS  <- renderValueBox(valueBox(round(last_val()$ROS_lvl,3), "ROS (AU)", icon=icon("fire"), color="purple"))
  output$bm_DBF  <- renderValueBox(valueBox(round(last_val()$DBF,2), "DBF (mL/min)", icon=icon("tint"), color="blue"))

  output$capillaro_plot <- renderPlotly({
    df <- sim_result()
    # Capillaroscopy score inversely related to DBF (structural damage in secondary)
    sec <- ifelse(input$subtype == "primary", 0, 1)
    df$cap_score <- 10 - df$DBF * (1 + sec * 0.3)
    df$cap_score <- pmax(pmin(df$cap_score, 10), 0)
    p <- ggplot(df, aes(x=time_wk, y=cap_score)) +
      geom_line(color="#E91E63", size=1.0) +
      labs(x="Time (weeks)", y="Capillaroscopy Score (0-10)") + theme_bw()
    ggplotly(p)
  })

  output$temp_recov_plot <- renderPlotly({
    df <- sim_result()
    # Temperature recovery time inversely proportional to DBF
    df$temp_recov <- 20.0 * exp(-df$DBF / 3.0) + 5.0
    p <- ggplot(df, aes(x=time_wk, y=temp_recov)) +
      geom_line(color="#0288D1", size=1.0) +
      labs(x="Time (weeks)", y="Temp Recovery Time (min)") + theme_bw()
    ggplotly(p)
  })

}

shinyApp(ui = ui, server = server)
