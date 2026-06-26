## =============================================================================
## SMA QSP Model — Interactive Shiny Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & Disease Subtype
##   2. Pharmacokinetics (PK)
##   3. Pharmacodynamics — SMN Biology
##   4. Motor Neuron & NMJ Dynamics
##   5. Clinical Endpoints (CMAP, HFMSE, FVC)
##   6. Treatment Scenario Comparison
##   7. Biomarker Dashboard
##   8. Virtual Population & Variability
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("sma_mrgsolve_model.R")

## ─────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = tags$span(
      tags$img(src = "https://upload.wikimedia.org/wikipedia/commons/thumb/1/17/Muscular-dystrophy.svg/40px-Muscular-dystrophy.svg.png",
               height = "30px", style = "margin-right:8px"),
      "SMA QSP Dashboard"
    ),
    titleWidth = 300
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "patient",   icon = icon("user")),
      menuItem("Pharmacokinetics",   tabName = "pk",        icon = icon("chart-line")),
      menuItem("SMN Biology (PD)",   tabName = "smn",       icon = icon("dna")),
      menuItem("Motor Neuron & NMJ", tabName = "mn_nmj",   icon = icon("brain")),
      menuItem("Clinical Endpoints", tabName = "clinical",  icon = icon("stethoscope")),
      menuItem("Scenario Comparison",tabName = "scenario",  icon = icon("balance-scale")),
      menuItem("Biomarker Dashboard",tabName = "biomarker", icon = icon("vial")),
      menuItem("Population Variability", tabName = "vpc",   icon = icon("users"))
    ),

    hr(),
    tags$div(style = "padding:10px; color:#ccc; font-size:11px;",
             tags$b("Global Controls"),
             sliderInput("sim_duration", "Simulation Duration (days)",
                         min = 90, max = 1825, value = 730, step = 30),
             selectInput("disease_type", "SMA Type",
                         choices = c("Type I (2 SMN2 copies)" = "type1",
                                     "Type II (3 SMN2 copies)" = "type2",
                                     "Type III (4 SMN2 copies)" = "type3",
                                     "Presymptomatic (2 copies)" = "presym")),
             selectInput("drug", "Treatment",
                         choices = c("No Treatment"    = "none",
                                     "Nusinersen"      = "nusinersen",
                                     "Risdiplam"       = "risdiplam",
                                     "Zolgensma (AAV9)"= "zolgensma")),
             numericInput("start_day", "Treatment Start Day", value = 0, min = 0, max = 365),
             actionButton("simulate", "Run Simulation",
                          class = "btn-primary btn-block",
                          icon = icon("play"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-radius: 8px; }
      .value-box { min-height: 80px; }
      .nav-tabs-custom { border-radius: 8px; }
    "))),

    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient & Disease Parameters", width = 6, status = "primary",
              sliderInput("smn2_copies", "SMN2 Copy Number", min=1, max=4, value=2, step=1),
              sliderInput("mn_death_rate", "Motor Neuron Death Rate (/day)",
                          min=0.0005, max=0.008, value=0.002, step=0.0005),
              sliderInput("smn_thresh", "SMN Threshold (fraction)",
                          min=0.1, max=0.6, value=0.3, step=0.05),
              numericInput("body_weight", "Body Weight (kg)", value=15, min=3, max=80),
              selectInput("age_group", "Age Group",
                          c("Neonatal (<3 mo)"="neonatal",
                            "Infant (3–12 mo)"="infant",
                            "Toddler (1–3 yr)"="toddler",
                            "Child (3–12 yr)"="child",
                            "Adult (>18 yr)"="adult"))
          ),
          box(title = "SMA Disease Classification", width = 6, status = "info",
              tags$table(class="table table-bordered table-hover", style="font-size:12px;",
                tags$thead(tags$tr(
                  tags$th("Type"), tags$th("Onset"), tags$th("SMN2"), tags$th("Max Ability"),
                  tags$th("Survival (untreated)")
                )),
                tags$tbody(
                  tags$tr(class="danger",
                    tags$td("0"), tags$td("Prenatal/neonatal"), tags$td("1"),
                    tags$td("None"), tags$td("<1 month")),
                  tags$tr(class="danger",
                    tags$td("I"), tags$td("<6 months"), tags$td("1–2"),
                    tags$td("Never sit"), tags$td("<2 years")),
                  tags$tr(class="warning",
                    tags$td("II"), tags$td("6–18 months"), tags$td("3"),
                    tags$td("Sit, never stand"), tags$td("Adult (reduced)")),
                  tags$tr(class="success",
                    tags$td("III"), tags$td(">18 months"), tags$td("3–4"),
                    tags$td("Walk (may lose)"), tags$td("Normal")),
                  tags$tr(class="success",
                    tags$td("IV"), tags$td(">21 years"), tags$td("≥4"),
                    tags$td("Walk throughout"), tags$td("Normal"))
                )
              ),
              hr(),
              h4("Disease Summary Metrics"),
              valueBoxOutput("vb_smn_protein", width=6),
              valueBoxOutput("vb_mn_pool",     width=6)
          )
        ),
        fluidRow(
          box(title = "Disease Mechanism Overview", width = 12, status = "warning",
              tags$img(src = "sma_qsp_model.svg",
                       style = "max-width:100%; border:1px solid #ccc; border-radius:4px;",
                       alt  = "SMA Mechanistic Map")
          )
        )
      ),

      ## ── Tab 2: Pharmacokinetics ─────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug-Specific PK Parameters", width = 4, status = "success",
              conditionalPanel("input.drug == 'nusinersen'",
                h4("Nusinersen PK"),
                numericInput("nus_dose", "IT Dose (mg)", 12, 1, 20),
                selectInput("nus_regime", "Dosing Regimen",
                            c("ENDEAR/CHERISH (4 loading + maintenance)"="standard",
                              "Loading only"="loading",
                              "Maintenance only (q4m)"="maintenance")),
                checkboxInput("nus_show_plasma", "Show Plasma Conc.", FALSE)
              ),
              conditionalPanel("input.drug == 'risdiplam'",
                h4("Risdiplam PK"),
                radioButtons("ris_dose_type", "Dose Type",
                             c("Fixed (5 mg/day)"="fixed", "Weight-based (0.2 mg/kg)"="wt")),
                numericInput("ris_dose_fixed", "Fixed Dose (mg)", 5, 0.5, 10, 0.5)
              ),
              conditionalPanel("input.drug == 'zolgensma'",
                h4("Zolgensma PK"),
                numericInput("zol_dose_vg_kg", "Dose (× 10^14 vg/kg)", 1.1, 0.5, 2.0, 0.1),
                sliderInput("aav9_ab", "Anti-AAV9 Ab Neutralization (0=none, 1=full)",
                            0, 1, 0, 0.1)
              )
          ),
          box(title = "PK Concentration-Time Profile", width = 8, status = "success",
              plotlyOutput("pk_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "PK Summary Table", width = 12, status = "info",
              DTOutput("pk_table"))
        )
      ),

      ## ── Tab 3: SMN Biology ──────────────────────────────────
      tabItem(tabName = "smn",
        fluidRow(
          box(title = "SMN PD Parameters", width = 4, status = "primary",
              h4("Exon-7 Splicing"),
              sliderInput("e7i_base",  "Baseline E7 Inclusion (SMN2)", 0.05, 0.20, 0.10, 0.01),
              sliderInput("emax_nus",  "Nusinersen Emax",  0.3, 0.9, 0.6, 0.05),
              sliderInput("ec50_nus",  "Nusinersen EC50 (ng/g)",  1, 20, 5, 0.5),
              sliderInput("emax_ris",  "Risdiplam Emax",   0.3, 0.8, 0.5, 0.05),
              sliderInput("ec50_ris",  "Risdiplam EC50 (ng/mL)", 20, 200, 80, 5),
              hr(),
              h4("SMN Protein Turnover"),
              sliderInput("k_prot_deg", "Protein Degradation Rate (/day)",
                          0.05, 0.5, 0.231, 0.01)
          ),
          tabBox(title = "SMN Biology Plots", width = 8,
            tabPanel("Exon 7 Inclusion",
              plotlyOutput("e7i_plot", height = "350px")),
            tabPanel("FL-SMN mRNA & Protein",
              plotlyOutput("smn_protein_plot", height = "350px")),
            tabPanel("Dose-Response",
              plotlyOutput("smn_dr_plot", height = "350px"))
          )
        )
      ),

      ## ── Tab 4: Motor Neuron & NMJ ────────────────────────────
      tabItem(tabName = "mn_nmj",
        fluidRow(
          box(title = "MN/NMJ Parameters", width = 4, status = "warning",
              sliderInput("mn_death",   "MN Death Rate (/day)", 0.0005, 0.01, 0.002, 0.0005),
              sliderInput("mn_rescue",  "MN Rescue Rate (/day)", 0.01, 0.1, 0.04, 0.005),
              sliderInput("nmj_mat",    "NMJ Maturation Rate (/day)", 0.001, 0.05, 0.01, 0.001),
              sliderInput("muscle_at",  "Muscle Atrophy Rate (/day)", 0.001, 0.01, 0.004, 0.001)
          ),
          tabBox(title = "Motor Neuron & NMJ Plots", width = 8,
            tabPanel("MN Pool Dynamics",
              plotlyOutput("mn_plot", height = "350px")),
            tabPanel("NMJ Maturation",
              plotlyOutput("nmj_plot", height = "350px")),
            tabPanel("Skeletal Muscle Mass",
              plotlyOutput("muscle_plot", height = "350px"))
          )
        )
      ),

      ## ── Tab 5: Clinical Endpoints ────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          valueBoxOutput("vb_cmap",   width=3),
          valueBoxOutput("vb_hfmse",  width=3),
          valueBoxOutput("vb_chop",   width=3),
          valueBoxOutput("vb_fvc",    width=3)
        ),
        fluidRow(
          box(title = "CMAP Amplitude", width = 6, status = "info",
              plotlyOutput("cmap_plot", height = "300px")),
          box(title = "HFMSE / CHOP-INTEND Score", width = 6, status = "info",
              plotlyOutput("motor_score_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "FVC % Predicted", width = 6, status = "warning",
              plotlyOutput("fvc_plot", height = "300px")),
          box(title = "RULM Score", width = 6, status = "warning",
              plotlyOutput("rulm_plot", height = "300px"))
        )
      ),

      ## ── Tab 6: Scenario Comparison ───────────────────────────
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Multi-Scenario Settings", width = 3, status = "primary",
              checkboxGroupInput("sc_drugs",
                "Treatments to Compare:",
                choices = c(
                  "No Treatment"           = "none",
                  "Nusinersen (early)"     = "nus_early",
                  "Nusinersen (late, 1yr)" = "nus_late",
                  "Risdiplam"              = "risdiplam",
                  "Zolgensma (presym)"     = "zolgensma"
                ),
                selected = c("none","nus_early","risdiplam","zolgensma")
              ),
              selectInput("sc_endpoint",
                "Primary Endpoint:",
                choices = c("SMN Protein"  = "SMN_protein",
                            "MN Pool"      = "MN_fraction",
                            "CMAP"         = "CMAP",
                            "CHOP-INTEND"  = "CHOP_INTEND",
                            "HFMSE"        = "HFMSE",
                            "FVC"          = "FVC",
                            "E7 Inclusion" = "E7_inclusion")),
              selectInput("sc_sma_type",
                "Disease Type for Comparison:",
                choices = c("Type I"="type1","Type II"="type2","Type III"="type3")),
              actionButton("sc_run", "Run All Scenarios", class="btn-warning btn-block",
                           icon=icon("play"))
          ),
          box(title = "Scenario Comparison Plot", width = 9, status = "primary",
              plotlyOutput("sc_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "Summary at Key Timepoints (Day 90, 180, 365, 730)", width = 12,
              DTOutput("sc_table"))
        )
      ),

      ## ── Tab 7: Biomarker Dashboard ───────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Overview", width = 12, status = "info",
              fluidRow(
                column(3,
                  h4("Diagnostic Biomarkers"),
                  tags$ul(
                    tags$li("SMN1/SMN2 MLPA copy number"),
                    tags$li("Newborn screening (DNA)"),
                    tags$li("Carrier testing")
                  )
                ),
                column(3,
                  h4("Pharmacodynamic Biomarkers"),
                  tags$ul(
                    tags$li("FL-SMN protein (blood/CSF)"),
                    tags$li("Exon 7 inclusion ratio"),
                    tags$li("snRNP assembly efficiency")
                  )
                ),
                column(3,
                  h4("Neuronal Injury Biomarkers"),
                  tags$ul(
                    tags$li("Plasma NF-L (neurofilament light)"),
                    tags$li("CSF NF-H (heavy chain)"),
                    tags$li("CMAP amplitude & latency")
                  )
                ),
                column(3,
                  h4("Muscle Biomarkers"),
                  tags$ul(
                    tags$li("Urine creatinine/creatine ratio"),
                    tags$li("DXA lean body mass"),
                    tags$li("MRI muscle fat fraction")
                  )
                )
              )
          )
        ),
        fluidRow(
          box(title = "SMN Protein Level Over Time", width = 6, status = "primary",
              plotlyOutput("bm_smn_plot", height = "300px")),
          box(title = "NF-L (Proxy) Over Time", width = 6, status = "danger",
              plotlyOutput("bm_nfl_plot", height = "300px"))
        ),
        fluidRow(
          box(title = "Biomarker Correlation Matrix", width = 6, status = "info",
              plotlyOutput("bm_corr_plot", height = "300px")),
          box(title = "Biomarker Response to Treatment", width = 6, status = "success",
              DTOutput("bm_table"))
        )
      ),

      ## ── Tab 8: Virtual Population ────────────────────────────
      tabItem(tabName = "vpc",
        fluidRow(
          box(title = "Virtual Population Settings", width = 4, status = "warning",
              numericInput("vpc_n", "Number of Patients", 50, 10, 200, 10),
              sliderInput("iiv_mn_death", "IIV on MN Death Rate (CV%)", 10, 50, 30, 5),
              sliderInput("iiv_ec50", "IIV on EC50 (CV%)", 20, 80, 40, 5),
              selectInput("vpc_drug", "Treatment for VPC",
                          c("Nusinersen"="nusinersen", "Risdiplam"="risdiplam",
                            "No Treatment"="none")),
              actionButton("vpc_run", "Generate Population",
                           class="btn-warning btn-block", icon=icon("users"))
          ),
          box(title = "Population Spread — CMAP Over Time", width = 8, status = "warning",
              plotlyOutput("vpc_plot", height = "400px"))
        ),
        fluidRow(
          box(title = "Population Responder Analysis at Day 365", width = 6, status = "info",
              plotlyOutput("vpc_resp_plot", height = "300px")),
          box(title = "Population Summary Statistics", width = 6, status = "info",
              DTOutput("vpc_table"))
        )
      )
    ) # end tabItems
  ) # end dashboardBody
)

## ─────────────────────────────────────────────────────────────
## Server
## ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive model
  base_mod <- reactive({ sma_model() })

  ## Map disease type to parameters
  type_params <- reactive({
    switch(input$disease_type,
      type1  = list(SMN2_copies=2, k_MN_death=0.004, k_MN_rescue=0.04),
      type2  = list(SMN2_copies=3, k_MN_death=0.002, k_MN_rescue=0.04),
      type3  = list(SMN2_copies=4, k_MN_death=0.001, k_MN_rescue=0.04),
      presym = list(SMN2_copies=2, k_MN_death=0.004, k_MN_rescue=0.08)
    )
  })

  ## Build dosing events
  dosing_events <- reactive({
    drug <- input$drug
    start <- input$start_day
    if (drug == "none") {
      return(mrgsolve::ev(time=99999, amt=0, cmt=1))
    } else if (drug == "nusinersen") {
      return(ev_nusinersen(start_day = start))
    } else if (drug == "risdiplam") {
      dose <- if (input$ris_dose_type == "wt") 0.2 * input$body_weight else input$ris_dose_fixed
      return(ev_risdiplam(start_day = start, end_day = input$sim_duration, dose_mg = dose))
    } else if (drug == "zolgensma") {
      dose_vg <- input$zol_dose_vg_kg * 1e14 * input$body_weight
      return(ev_zolgensma(dose_vg = dose_vg))
    }
  })

  ## Run main simulation (triggered by button)
  sim_result <- eventReactive(input$simulate, {
    mod2 <- mrgsolve::param(base_mod(),
             c(type_params(),
               list(E7I_base    = input$e7i_base,
                    Emax_NUS    = input$emax_nus,
                    EC50_NUS    = input$ec50_nus,
                    Emax_RIS    = input$emax_ris,
                    EC50_RIS    = input$ec50_ris,
                    k_prot_deg  = input$k_prot_deg,
                    k_MN_death  = input$mn_death,
                    k_MN_rescue = input$mn_rescue,
                    k_NMJ_mat   = input$nmj_mat,
                    k_muscle_at = input$muscle_at)))
    ev  <- dosing_events()
    df  <- as.data.frame(mrgsolve::mrgsim(mod2, ev, delta=1, end=input$sim_duration, obsonly=TRUE))
    df
  }, ignoreNULL = FALSE)

  ## Value Boxes
  last_val <- function(col) {
    df <- sim_result()
    if (is.null(df) || nrow(df)==0) return(NA)
    tail(df[[col]], 1)
  }

  output$vb_smn_protein <- renderValueBox({
    v <- round(last_val("SMN_protein"), 3)
    valueBox(v, "SMN Protein (norm)", icon=icon("dna"), color="blue")
  })
  output$vb_mn_pool <- renderValueBox({
    v <- round(last_val("MN_fraction") * 100, 1)
    valueBox(paste0(v, "%"), "Motor Neuron Pool", icon=icon("brain"), color="green")
  })
  output$vb_cmap <- renderValueBox({
    v <- round(last_val("CMAP"), 1)
    valueBox(paste0(v, " mV"), "CMAP Amplitude", icon=icon("bolt"), color="blue")
  })
  output$vb_hfmse <- renderValueBox({
    v <- round(last_val("HFMSE"), 1)
    valueBox(v, "HFMSE Score", icon=icon("walking"), color="green")
  })
  output$vb_chop <- renderValueBox({
    v <- round(last_val("CHOP_INTEND"), 1)
    valueBox(v, "CHOP-INTEND", icon=icon("baby"), color="yellow")
  })
  output$vb_fvc <- renderValueBox({
    v <- round(last_val("FVC"), 1)
    valueBox(paste0(v, "%"), "FVC % Predicted", icon=icon("lungs"), color="orange")
  })

  ## PK plot
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    drug <- input$drug
    if (drug == "nusinersen") {
      p <- ggplot(df, aes(time, C_CSF_lumbar)) +
        geom_line(color="#2196F3", size=1.2) +
        labs(title="Nusinersen — CSF Lumbar Concentration",
             x="Day", y="Conc. (ng/mL)") + theme_bw(13)
    } else if (drug == "risdiplam") {
      p <- ggplot(df, aes(time, C_plasma_risdiplam)) +
        geom_line(color="#4CAF50", size=1.2) +
        labs(title="Risdiplam — Plasma Concentration",
             x="Day", y="Conc. (ng/mL)") + theme_bw(13)
    } else if (drug == "zolgensma") {
      p <- ggplot(df, aes(time, Transgene_mRNA)) +
        geom_line(color="#9C27B0", size=1.2) +
        labs(title="Zolgensma — Transgene FL-SMN mRNA",
             x="Day", y="mRNA (relative units)") + theme_bw(13)
    } else {
      p <- ggplot() + annotate("text", x=0.5, y=0.5, label="No drug selected") + theme_void()
    }
    ggplotly(p)
  })

  ## SMN plots
  output$e7i_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, E7_inclusion)) +
      geom_line(color="#F44336", size=1.2) +
      geom_hline(yintercept=0.10, linetype="dashed", color="grey50") +
      annotate("text", x=10, y=0.12, label="Baseline (10%)", color="grey40", size=3) +
      labs(title="Exon 7 Inclusion Rate (SMN2)",
           x="Day", y="Exon 7 Inclusion Fraction") + theme_bw(13)
    ggplotly(p)
  })

  output$smn_protein_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, SMN_protein)) +
      geom_line(color="#673AB7", size=1.2) +
      geom_hline(yintercept=0.30, linetype="dashed", color="red") +
      annotate("text", x=10, y=0.32, label="SMN threshold", color="red", size=3) +
      labs(title="Full-Length SMN Protein Pool",
           x="Day", y="SMN Protein (normalized)") + theme_bw(13)
    ggplotly(p)
  })

  output$smn_dr_plot <- renderPlotly({
    conc_range <- 10^seq(-1, 3, length.out=100)
    emax_n <- input$emax_nus; ec50_n <- input$ec50_nus
    emax_r <- input$emax_ris; ec50_r <- input$ec50_ris
    e7_base <- input$e7i_base
    df_dr <- data.frame(
      conc = conc_range,
      E7I_NUS = e7_base + (0.90-e7_base) * emax_n*conc_range^1.5 / (ec50_n^1.5 + conc_range^1.5),
      E7I_RIS = e7_base + (0.90-e7_base) * emax_r*conc_range^1.2 / (ec50_r^1.2 + conc_range^1.2)
    ) %>% pivot_longer(-conc, names_to="Drug", values_to="E7I")
    p <- ggplot(df_dr, aes(conc, E7I, color=Drug)) +
      geom_line(size=1.2) + scale_x_log10() +
      scale_color_manual(values=c("E7I_NUS"="#2196F3","E7I_RIS"="#4CAF50")) +
      labs(title="Concentration-Exon 7 Inclusion Response", x="Conc. (ng/mL or ng/g)",
           y="E7 Inclusion Fraction") + theme_bw(13)
    ggplotly(p)
  })

  ## Motor Neuron plots
  output$mn_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, MN_fraction)) +
      geom_line(color="#E91E63", size=1.2) +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="Alpha Motor Neuron Pool", x="Day",
           y="MN Pool (fraction of baseline)") + theme_bw(13)
    ggplotly(p)
  })

  output$nmj_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, NMJ_maturity)) +
      geom_line(color="#FF5722", size=1.2) +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="NMJ Maturation Score", x="Day", y="NMJ Score (0–1)") + theme_bw(13)
    ggplotly(p)
  })

  output$muscle_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, Muscle_mass)) +
      geom_line(color="#FF9800", size=1.2) +
      scale_y_continuous(limits=c(0,1.05)) +
      labs(title="Skeletal Muscle Mass", x="Day",
           y="Muscle Mass (normalized)") + theme_bw(13)
    ggplotly(p)
  })

  ## Clinical endpoint plots
  output$cmap_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, CMAP)) + geom_line(color="#1565C0", size=1.2) +
      labs(title="CMAP Amplitude", x="Day", y="CMAP (mV)") + theme_bw(13)
    ggplotly(p)
  })

  output$motor_score_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df) +
      geom_line(aes(time, HFMSE, color="HFMSE"), size=1.2) +
      geom_line(aes(time, CHOP_INTEND, color="CHOP-INTEND"), size=1.2) +
      scale_color_manual(values=c("HFMSE"="#4CAF50","CHOP-INTEND"="#FF9800")) +
      labs(title="Motor Function Scores", x="Day", y="Score") + theme_bw(13)
    ggplotly(p)
  })

  output$fvc_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, FVC)) + geom_line(color="#00796B", size=1.2) +
      geom_hline(yintercept=40, linetype="dashed", color="red") +
      annotate("text", x=10, y=42, label="Ventilator threshold", color="red", size=3) +
      labs(title="FVC % Predicted", x="Day", y="FVC (%)") + theme_bw(13)
    ggplotly(p)
  })

  output$rulm_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    p <- ggplot(df, aes(time, RULM)) + geom_line(color="#7B1FA2", size=1.2) +
      labs(title="RULM Score", x="Day", y="RULM (0–37)") + theme_bw(13)
    ggplotly(p)
  })

  ## Scenario comparison
  sc_result <- eventReactive(input$sc_run, {
    req(length(input$sc_drugs) > 0)
    mod2 <- base_mod()
    sc_type_params <- switch(input$sc_sma_type,
      type1=list(SMN2_copies=2, k_MN_death=0.004),
      type2=list(SMN2_copies=3, k_MN_death=0.002),
      type3=list(SMN2_copies=4, k_MN_death=0.001))

    sc_list <- list(
      none      = list(name="No Treatment",      events=mrgsolve::ev(time=99999,amt=0,cmt=1)),
      nus_early = list(name="Nusinersen (early)", events=ev_nusinersen(0)),
      nus_late  = list(name="Nusinersen (late)",  events=ev_nusinersen(365)),
      risdiplam = list(name="Risdiplam",          events=ev_risdiplam(0,730,5)),
      zolgensma = list(name="Zolgensma",          events=ev_zolgensma(1.65e15))
    )

    results <- lapply(input$sc_drugs, function(sc_id) {
      sc  <- sc_list[[sc_id]]
      m   <- mrgsolve::param(mod2, sc_type_params)
      df  <- as.data.frame(mrgsolve::mrgsim(m, sc$events, delta=1, end=730, obsonly=TRUE))
      df$Scenario <- sc$name
      df
    })
    bind_rows(results)
  })

  output$sc_plot <- renderPlotly({
    df <- sc_result()
    req(!is.null(df))
    ep <- input$sc_endpoint
    p <- ggplot(df, aes_string("time", ep, color="Scenario")) +
      geom_line(size=1.2) +
      labs(title=paste("Scenario Comparison —", ep), x="Day", y=ep) +
      theme_bw(13) + theme(legend.position="bottom")
    ggplotly(p)
  })

  output$sc_table <- renderDT({
    df <- sc_result()
    req(!is.null(df))
    ep <- input$sc_endpoint
    df %>% filter(time %in% c(90,180,365,730)) %>%
      group_by(Scenario, time) %>% slice(1) %>%
      select(Scenario, time, all_of(ep)) %>%
      pivot_wider(names_from=time, values_from=all_of(ep), names_prefix="Day_") %>%
      datatable(options=list(dom="t", pageLength=10), rownames=FALSE)
  })

  ## Biomarker plots
  output$bm_smn_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    ggplotly(ggplot(df, aes(time, SMN_protein)) + geom_line(color="#7B1FA2", size=1.2) +
      labs(title="SMN Protein", x="Day", y="SMN (normalized)") + theme_bw(13))
  })

  output$bm_nfl_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df))
    df$NF_L_proxy <- 100 * (1 - df$MN_fraction)   # proxy: MN loss → NF-L
    ggplotly(ggplot(df, aes(time, NF_L_proxy)) + geom_line(color="#F44336", size=1.2) +
      labs(title="NF-L (Proxy — Scaled MN Loss)", x="Day", y="NF-L (a.u.)") + theme_bw(13))
  })

  output$bm_corr_plot <- renderPlotly({
    df <- sim_result()
    req(!is.null(df) && nrow(df) > 5)
    bm_df <- df %>% select(SMN_protein, MN_fraction, NMJ_maturity, Muscle_mass, CMAP, HFMSE, FVC)
    cor_mat <- cor(bm_df, use="complete.obs")
    p <- plotly::plot_ly(z=cor_mat, type="heatmap",
                         colorscale="RdBu", zmin=-1, zmax=1,
                         x=colnames(cor_mat), y=rownames(cor_mat)) %>%
      plotly::layout(title="Biomarker Correlation Matrix")
    p
  })

  output$bm_table <- renderDT({
    df <- sim_result()
    req(!is.null(df))
    df %>% filter(time %in% c(0, 90, 180, 365, 730)) %>%
      select(time, SMN_protein, MN_fraction, CMAP, HFMSE, FVC) %>%
      mutate(across(where(is.numeric), ~round(., 3))) %>%
      datatable(rownames=FALSE, options=list(dom="t"))
  })

  ## Virtual Population
  vpc_data <- eventReactive(input$vpc_run, {
    n <- input$vpc_n
    set.seed(2024)
    mod2 <- base_mod()
    idata <- data.frame(
      ID           = 1:n,
      SMN2_copies  = 2,
      k_MN_death   = rlnorm(n, log(0.002), input$iiv_mn_death/100),
      EC50_NUS     = rlnorm(n, log(5.0), input$iiv_ec50/100),
      Emax_NUS     = pmin(0.9, pmax(0.3, rnorm(n, 0.60, 0.08)))
    )
    ev <- if (input$vpc_drug == "nusinersen") ev_nusinersen() else
          if (input$vpc_drug == "risdiplam")  ev_risdiplam(dose_mg=5) else
          mrgsolve::ev(time=99999, amt=0, cmt=1)

    out <- mrgsolve::mrgsim(mod2, idata=idata, events=ev, delta=14, end=730, obsonly=TRUE)
    as.data.frame(out)
  })

  output$vpc_plot <- renderPlotly({
    df <- vpc_data()
    req(!is.null(df))
    stats <- df %>% group_by(time) %>%
      summarise(med=median(CMAP), lo=quantile(CMAP, 0.05), hi=quantile(CMAP, 0.95))
    p <- ggplot(stats, aes(time)) +
      geom_ribbon(aes(ymin=lo, ymax=hi), fill="#2196F3", alpha=0.3) +
      geom_line(aes(y=med), color="#1565C0", size=1.5) +
      labs(title="VPC: CMAP Amplitude (Median ± 90% PI)", x="Day", y="CMAP (mV)") +
      theme_bw(13)
    ggplotly(p)
  })

  output$vpc_resp_plot <- renderPlotly({
    df <- vpc_data()
    req(!is.null(df))
    day365 <- df %>% filter(abs(time-364) < 8) %>% group_by(ID) %>% slice(1)
    p <- ggplot(day365, aes(x=CMAP)) +
      geom_histogram(bins=20, fill="#4CAF50", color="white", alpha=0.8) +
      labs(title="Distribution of CMAP at Day 365 (Population)", x="CMAP (mV)", y="Count") +
      theme_bw(13)
    ggplotly(p)
  })

  output$vpc_table <- renderDT({
    df <- vpc_data()
    req(!is.null(df))
    df %>% filter(time %in% c(0, 182, 365, 730)) %>%
      group_by(time) %>%
      summarise(
        N=n(),
        CMAP_med=round(median(CMAP),2), CMAP_p5=round(quantile(CMAP,0.05),2),
        CMAP_p95=round(quantile(CMAP,0.95),2),
        HFMSE_med=round(median(HFMSE),1),
        FVC_med=round(median(FVC),1)
      ) %>%
      datatable(rownames=FALSE, options=list(dom="t"))
  })

  ## PK table
  output$pk_table <- renderDT({
    df <- sim_result()
    req(!is.null(df))
    df %>% filter(time %in% c(0,14,28,63,180,365,730)) %>%
      select(time, C_CSF_lumbar, C_CNS_nusinersen, C_plasma_risdiplam, Transgene_mRNA) %>%
      mutate(across(where(is.numeric), ~round(.,3))) %>%
      datatable(rownames=FALSE, options=list(dom="t"))
  })
}

## ─────────────────────────────────────────────────────────────
## Run App
## ─────────────────────────────────────────────────────────────
shinyApp(ui, server)
