# ==============================================================================
# Rheumatoid Arthritis QSP Shiny Dashboard
# 류마티스 관절염 QSP 인터랙티브 대시보드
# ==============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(tidyverse)
library(plotly)
library(DT)

# ---- mrgsolve model code (embedded) -----------------------------------------
ra_code <- '
$PARAM
CL_TCZ=0.52, V1_TCZ=3.72, V2_TCZ=2.35, Q_TCZ=0.56,
F_SC=0.80, KA_TCZ=0.108,
KON=0.0033, KOFF=0.019, KSYN_R=0.1, KDEG_R=0.007, RTOT0=14.3,
CL_ADA=0.0114, V1_ADA=2.84, KA_ADA=0.0073,
KA_MTX=0.3, CL_MTX=5.4, V_MTX=28, KPOLY=0.02, KDEPOLY=0.005,
KA_BARI=0.69, CL_BARI=12.1, V_BARI=76, IC50_JAK=5.9,
KOUT_IL6=0.0289, KIN_IL6=0.0722, KOUT_STAT3=0.1,
IL6_EC50=1.5, STAT3MAX=100,
KOUT_CRP=0.0289, KIN_CRP0=2.89, CRP_EC50=0.5,
K_INFLAM_IN=0.01, K_INFLAM_OUT=0.001, INFLAM_BASE=10,
K_STAT3_INFLAM=0.2, K_TNF_INFLAM=0.15,
KFLS_IN=0.001, KFLS_OUT=0.0005,
KVDH_PROG=0.00015, KVDH_INHIB=0.7,
K_CART_DESTR=0.00008

$CMT SC_TCZ CENT_TCZ PERI_TCZ IL6R_FREE DRUG_RL
     SC_ADA CENT_ADA
     GUT_MTX CENT_MTX POLY_MTX
     GUT_BARI CENT_BARI
     IL6 STAT3_P CRP INFLAM FLS_ACT VDH CART

$MAIN
double C_TCZ  = CENT_TCZ / V1_TCZ;
double C_ADA  = CENT_ADA / V1_ADA;
double C_MTX  = CENT_MTX / V_MTX;
double C_BARI = CENT_BARI / V_BARI;
double RTOT   = IL6R_FREE + DRUG_RL;
double OCC_TCZ = DRUG_RL / (RTOT + 1e-6);
double JAK_INH = C_BARI / (C_BARI + IC50_JAK);
double IL6_sig = (IL6/(IL6+IL6_EC50))*(1-OCC_TCZ)*(1-0.8*JAK_INH);
double MTX_EFF = 0.35*POLY_MTX/(POLY_MTX+0.5);
double TNF_INH = C_ADA/(C_ADA+0.035);
double DRUG_AI = 1-(1-OCC_TCZ)*(1-TNF_INH)*(1-MTX_EFF)*(1-0.6*JAK_INH);

$ODE
dxdt_SC_TCZ   = -KA_TCZ*SC_TCZ*F_SC;
dxdt_CENT_TCZ = KA_TCZ*SC_TCZ*F_SC-(CL_TCZ/V1_TCZ)*CENT_TCZ
                -(Q_TCZ/V1_TCZ)*CENT_TCZ+(Q_TCZ/V2_TCZ)*PERI_TCZ
                -V1_TCZ*(KON*C_TCZ*IL6R_FREE-KOFF*DRUG_RL);
dxdt_PERI_TCZ = (Q_TCZ/V1_TCZ)*CENT_TCZ-(Q_TCZ/V2_TCZ)*PERI_TCZ;
dxdt_IL6R_FREE= KSYN_R-KDEG_R*IL6R_FREE-KON*C_TCZ*IL6R_FREE+KOFF*DRUG_RL;
dxdt_DRUG_RL  = KON*C_TCZ*IL6R_FREE-KOFF*DRUG_RL-KDEG_R*DRUG_RL;
dxdt_SC_ADA   = -KA_ADA*SC_ADA;
dxdt_CENT_ADA = KA_ADA*SC_ADA-CL_ADA*C_ADA;
dxdt_GUT_MTX  = -KA_MTX*GUT_MTX;
dxdt_CENT_MTX = KA_MTX*GUT_MTX-CL_MTX*C_MTX;
dxdt_POLY_MTX = KPOLY*C_MTX-KDEPOLY*POLY_MTX;
dxdt_GUT_BARI = -KA_BARI*GUT_BARI;
dxdt_CENT_BARI= KA_BARI*GUT_BARI-(CL_BARI/V_BARI)*CENT_BARI;
dxdt_IL6      = KIN_IL6*(1+0.5*INFLAM/INFLAM_BASE)-KOUT_IL6*IL6;
dxdt_STAT3_P  = STAT3MAX*IL6_sig-KOUT_STAT3*STAT3_P;
dxdt_CRP      = KIN_CRP0*STAT3_P/(STAT3_P+CRP_EC50*STAT3MAX)-KOUT_CRP*CRP;
dxdt_INFLAM   = K_INFLAM_IN*(K_STAT3_INFLAM*STAT3_P/STAT3MAX+K_TNF_INFLAM*(1-TNF_INH))*INFLAM_BASE
                -(K_INFLAM_OUT+K_INFLAM_IN*DRUG_AI)*INFLAM;
dxdt_FLS_ACT  = KFLS_IN*INFLAM-KFLS_OUT*FLS_ACT*(1+2*DRUG_AI);
dxdt_VDH      = KVDH_PROG*INFLAM*(1-KVDH_INHIB*DRUG_AI);
dxdt_CART     = -K_CART_DESTR*INFLAM*(1-DRUG_AI)*(CART/100.0);

$TABLE
double TJC = 14.0*(INFLAM/20.0);
double SJC = 10.0*(INFLAM/20.0);
double GH  = 50.0*(INFLAM/20.0);
double DAS28 = 0.56*sqrt(TJC+0.01)+0.28*sqrt(SJC+0.01)+0.36*log(CRP+1)+0.014*GH;
double HAQ_DI = 1.5*(INFLAM/20.0)+0.02*VDH;
double RO_pct = 100*OCC_TCZ;
double JAK_pct = 100*C_BARI/(C_BARI+IC50_JAK);

$CAPTURE DAS28 CRP INFLAM VDH CART HAQ_DI RO_pct JAK_pct
         C_TCZ C_ADA C_MTX C_BARI STAT3_P FLS_ACT
'

# Build the model once at startup
ra_mod <- mcode("RA_Shiny", ra_code, quiet = TRUE)

# ---- Simulation function -----------------------------------------------------
run_sim <- function(scenario, weight_kg = 70, baseline_inflam = 10, baseline_crp = 15,
                    baseline_das28 = 6.2, sim_weeks = 52) {

  sim_end <- sim_weeks * 168  # hours

  init_vals <- init(ra_mod,
    IL6R_FREE = 14.3,
    IL6       = 2.5,
    STAT3_P   = 50 * (baseline_crp / 15),
    CRP       = baseline_crp,
    INFLAM    = baseline_inflam,
    FLS_ACT   = baseline_inflam * 0.5,
    VDH       = 10,
    CART      = max(100 - (baseline_das28 - 2.6) * 3, 60)
  )

  # Build dosing event based on scenario
  e <- switch(scenario,
    "None"        = ev(time = 0, amt = 0, cmt = 1),
    "MTX"         = ev(time = seq(0, sim_end - 1, 168), amt = 20, cmt = "GUT_MTX"),
    "TCZ_IV"      = ev(time = seq(0, sim_end - 1, 672), amt = weight_kg * 8 / 148 * 1000, cmt = "CENT_TCZ"),
    "TCZ_SC"      = ev(time = seq(0, sim_end - 1, 336), amt = 1095, cmt = "SC_TCZ"),
    "TCZ_MTX"     = ev(time = seq(0, sim_end - 1, 672), amt = weight_kg * 8 / 148 * 1000, cmt = "CENT_TCZ") +
                    ev(time = seq(0, sim_end - 1, 168), amt = 20, cmt = "GUT_MTX"),
    "Baricitinib" = ev(time = seq(0, sim_end - 1, 24), amt = 4, cmt = "GUT_BARI"),
    "Adalimumab"  = ev(time = seq(0, sim_end - 1, 336), amt = 270, cmt = "SC_ADA"),
    ev(time = 0, amt = 0, cmt = 1)
  )

  ra_mod %>%
    init(init_vals) %>%
    param(INFLAM_BASE = baseline_inflam) %>%
    mrgsim(ev = e, end = sim_end, delta = 24, quiet = TRUE) %>%
    as_tibble() %>%
    mutate(week = time / 168, scenario = scenario)
}

# ---- UI ----------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "RA QSP Dashboard | 류마티스 관절염"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "patient",    icon = icon("user")),
      menuItem("Drug PK",            tabName = "pk",         icon = icon("pills")),
      menuItem("IL-6 & STAT3",       tabName = "signaling",  icon = icon("dna")),
      menuItem("Disease Activity",   tabName = "disease",    icon = icon("chart-line")),
      menuItem("Structural Damage",  tabName = "structural", icon = icon("bone")),
      menuItem("Scenario Comparison",tabName = "comparison", icon = icon("table")),
      menuItem("Biomarker Panel",    tabName = "biomarkers", icon = icon("vial"))
    ),
    hr(),
    h5("Simulation Settings", style = "color:white; margin-left:15px"),
    sliderInput("weight",        "Patient Weight (kg)",       40,  120, 70,  step = 5),
    sliderInput("baseline_das28","Baseline DAS28",            3.2,  9.0, 6.2, step = 0.1),
    sliderInput("baseline_crp",  "Baseline CRP (mg/L)",       5,    80,  15,  step = 1),
    sliderInput("sim_weeks",     "Simulation Duration (weeks)",12,  104, 52,  step = 4),
    checkboxGroupInput("scenarios", "Scenarios to Compare",
      choices = list(
        "No Treatment"         = "None",
        "MTX 20mg/wk"          = "MTX",
        "TCZ IV 8mg/kg"        = "TCZ_IV",
        "TCZ SC 162mg q2w"     = "TCZ_SC",
        "TCZ IV + MTX"         = "TCZ_MTX",
        "Baricitinib 4mg QD"   = "Baricitinib",
        "Adalimumab 40mg q2w"  = "Adalimumab"
      ),
      selected = c("None", "MTX", "TCZ_MTX", "Baricitinib")
    ),
    actionButton("run_sim", "Run Simulation", class = "btn-success btn-block")
  ),

  dashboardBody(
    tabItems(

      # ---- Tab 1: Patient Profile -------------------------------------------
      tabItem("patient",
        fluidRow(
          box(title = "Patient Characteristics", status = "primary", solidHeader = TRUE, width = 6,
            valueBoxOutput("vbox_das28",  width = 12),
            valueBoxOutput("vbox_crp",    width = 12),
            valueBoxOutput("vbox_weight", width = 12)
          ),
          box(title = "Disease State at Baseline", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("patient_radar", height = 300),
            p("Baseline inflammatory state estimation based on DAS28 and CRP inputs.")
          )
        ),
        fluidRow(
          box(title = "About Rheumatoid Arthritis", status = "info", solidHeader = TRUE, width = 12,
            p("Rheumatoid Arthritis (RA) is a chronic, systemic autoimmune disease characterised by
              persistent synovial inflammation, progressive joint destruction, and extra-articular manifestations.
              Key pathogenic drivers include TNF-alpha, IL-6, and IL-17A, acting through JAK-STAT signalling
              to activate synovial fibroblasts and osteoclastogenesis."),
            p("This QSP model integrates IL-6 receptor occupancy dynamics (TMDD), STAT3 phosphorylation,
              CRP kinetics, synovial inflammation, and structural damage (Sharp/vdH score) into a unified
              computational framework supporting seven treatment strategies: no treatment, MTX monotherapy,
              tocilizumab IV/SC, TCZ+MTX combination, baricitinib, and adalimumab.")
          )
        )
      ),

      # ---- Tab 2: Drug PK ---------------------------------------------------
      tabItem("pk",
        fluidRow(
          box(title = "Drug Concentration-Time Profiles", status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("pk_plot", height = 450)
          )
        ),
        fluidRow(
          box(title = "TCZ IL-6R Occupancy", status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("ro_plot", height = 300)
          ),
          box(title = "MTX Polyglutamate Accumulation", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("mtx_poly_plot", height = 300)
          )
        )
      ),

      # ---- Tab 3: IL-6 & STAT3 Signalling -----------------------------------
      tabItem("signaling",
        fluidRow(
          box(title = "IL-6 Dynamics", status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("il6_plot", height = 350)
          ),
          box(title = "pSTAT3 Level", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("stat3_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "FLS Activation", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("fls_plot", height = 300)
          ),
          box(title = "Signalling Pathway Summary", status = "info", solidHeader = TRUE, width = 6,
            p(strong("IL-6 -> JAK1/JAK2 -> STAT3 -> CRP / MMP / Th17 differentiation")),
            p("Tocilizumab and sarilumab block membrane and soluble IL-6 receptor,
              preventing gp130 signalling and STAT3 phosphorylation."),
            p("Baricitinib (JAK1/2 inhibitor) and tofacitinib (JAK1/3)
              directly block downstream JAK kinases, effective against multiple cytokines."),
            plotlyOutput("jak_inh_plot", height = 200)
          )
        )
      ),

      # ---- Tab 4: Disease Activity ------------------------------------------
      tabItem("disease",
        fluidRow(
          box(title = "DAS28-CRP over Time", status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("das28_plot", height = 400)
          ),
          box(title = "Response Summary at Key Timepoints", status = "success", solidHeader = TRUE, width = 4,
            DTOutput("response_table")
          )
        ),
        fluidRow(
          box(title = "CRP Trajectory", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("crp_plot", height = 300)
          ),
          box(title = "Synovial Inflammation Index", status = "danger", solidHeader = TRUE, width = 6,
            plotlyOutput("inflam_plot", height = 300)
          )
        )
      ),

      # ---- Tab 5: Structural Damage -----------------------------------------
      tabItem("structural",
        fluidRow(
          box(title = "Sharp/van der Heijde Score Progression", status = "danger", solidHeader = TRUE, width = 8,
            plotlyOutput("sharp_plot", height = 400)
          ),
          box(title = "Radiographic Progression by Scenario", status = "info", solidHeader = TRUE, width = 4,
            p("Sharp/vdH scores at end of simulation:"),
            DTOutput("sharp_table")
          )
        ),
        fluidRow(
          box(title = "Cartilage Integrity (%)", status = "warning", solidHeader = TRUE, width = 6,
            plotlyOutput("cart_plot", height = 300)
          ),
          box(title = "HAQ-DI (Functional Disability)", status = "primary", solidHeader = TRUE, width = 6,
            plotlyOutput("haq_plot", height = 300)
          )
        )
      ),

      # ---- Tab 6: Scenario Comparison ---------------------------------------
      tabItem("comparison",
        fluidRow(
          box(title = "Head-to-Head Comparison at End of Simulation", status = "primary",
              solidHeader = TRUE, width = 12,
            DTOutput("comparison_table")
          )
        ),
        fluidRow(
          box(title = "Remission Rate (DAS28 <2.6)", status = "success", solidHeader = TRUE, width = 6,
            plotlyOutput("remission_bar", height = 300)
          ),
          box(title = "ACR Response Rates", status = "info", solidHeader = TRUE, width = 6,
            plotlyOutput("acr_bar", height = 300)
          )
        )
      ),

      # ---- Tab 7: Biomarker Panel -------------------------------------------
      tabItem("biomarkers",
        fluidRow(
          box(title = "Multi-Biomarker Dashboard", status = "primary", solidHeader = TRUE, width = 12,
            plotlyOutput("biomarker_panel", height = 500)
          )
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges", status = "info", solidHeader = TRUE, width = 6,
            tableOutput("biomarker_ref")
          ),
          box(title = "Drug-Biomarker Mechanistic Summary", status = "warning", solidHeader = TRUE, width = 6,
            p(strong("Key Biomarkers in RA:")),
            tags$ul(
              tags$li("CRP: acute phase reactant driven by IL-6/STAT3 signalling to hepatocytes"),
              tags$li("ESR: elevated by fibrinogen (IL-6 driven) and IgG immune complexes"),
              tags$li("Anti-CCP (ACPA): highly specific (>95%), precedes symptoms by years"),
              tags$li("RF (IgM): ~70% sensitivity, also elevated in other inflammatory conditions"),
              tags$li("Serum IL-6: accumulates after TCZ receptor blockade (paradoxical rise)"),
              tags$li("TNF-alpha: blocked by adalimumab/etanercept/certolizumab"),
              tags$li("pSTAT3: intracellular signalling marker (detectable in synovial biopsy)"),
              tags$li("MMP-3: FLS-derived metalloproteinase, correlates with structural damage progression")
            )
          )
        )
      )
    )
  )
)

# ---- Server ------------------------------------------------------------------
server <- function(input, output, session) {

  # Colour palette consistent across all plots
  scen_colors <- c(
    "None"        = "tomato",
    "MTX"         = "darkorange",
    "TCZ_IV"      = "steelblue",
    "TCZ_SC"      = "skyblue",
    "TCZ_MTX"     = "navy",
    "Baricitinib" = "purple",
    "Adalimumab"  = "darkgreen"
  )

  # ---- Reactive simulation --------------------------------------------------
  sim_data <- eventReactive(input$run_sim, {
    req(input$scenarios)

    withProgress(message = "Running QSP simulation ...", value = 0, {
      results <- map_dfr(input$scenarios, function(sc) {
        incProgress(1 / length(input$scenarios), detail = paste("Scenario:", sc))
        tryCatch(
          run_sim(
            sc,
            weight_kg       = input$weight,
            baseline_inflam = 8 + (input$baseline_das28 - 3.2) * 0.7,
            baseline_crp    = input$baseline_crp,
            baseline_das28  = input$baseline_das28,
            sim_weeks       = input$sim_weeks
          ),
          error = function(e) {
            message("Simulation error for scenario ", sc, ": ", conditionMessage(e))
            tibble()
          }
        )
      })
    })
    results
  }, ignoreNULL = FALSE)

  # ---- Value boxes ----------------------------------------------------------
  output$vbox_das28 <- renderValueBox({
    col <- if (input$baseline_das28 > 5.1) "red" else if (input$baseline_das28 > 3.2) "yellow" else "green"
    valueBox(input$baseline_das28, "Baseline DAS28-CRP",
             icon = icon("chart-bar"), color = col)
  })

  output$vbox_crp <- renderValueBox(
    valueBox(paste0(input$baseline_crp, " mg/L"), "Baseline CRP",
             icon = icon("vial"), color = "orange")
  )

  output$vbox_weight <- renderValueBox(
    valueBox(paste0(input$weight, " kg"), "Patient Weight",
             icon = icon("weight"), color = "blue")
  )

  # ---- Patient radar (polar bar) -------------------------------------------
  output$patient_radar <- renderPlotly({
    bm <- tibble(
      Marker = c("DAS28", "CRP/10", "Inflammation", "HAQ est."),
      Value  = c(
        input$baseline_das28,
        input$baseline_crp / 10,
        8 + (input$baseline_das28 - 3.2) * 0.7,
        1.0
      )
    )
    p <- ggplot(bm, aes(Marker, Value, fill = Marker)) +
      geom_col() +
      coord_polar() +
      labs(title = "Patient Disease Burden Profile") +
      theme_minimal() +
      theme(legend.position = "none")
    ggplotly(p)
  })

  # ---- DAS28 ----------------------------------------------------------------
  output$das28_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    max_wk <- max(d$week, na.rm = TRUE)
    p <- ggplot(d, aes(week, DAS28, color = scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 2.6, linetype = "dashed", color = "green4") +
      geom_hline(yintercept = 3.2, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 5.1, linetype = "dashed", color = "red", alpha = 0.5) +
      annotate("text", x = max_wk * 0.92, y = 2.4, label = "Remission <2.6", color = "green4", size = 2.8) +
      annotate("text", x = max_wk * 0.92, y = 3.0, label = "LDA <3.2",       color = "orange",  size = 2.8) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "DAS28-CRP", color = "Scenario", title = "DAS28-CRP Trajectory") +
      theme_minimal(base_size = 12)
    ggplotly(p, tooltip = c("week", "DAS28", "scenario"))
  })

  # ---- CRP -----------------------------------------------------------------
  output$crp_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, CRP, color = scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "blue", alpha = 0.6) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "CRP (mg/L)", color = "", title = "CRP Dynamics") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- Inflammation --------------------------------------------------------
  output$inflam_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, INFLAM, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "Inflammation Index", color = "", title = "Synovial Inflammation") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- PK ------------------------------------------------------------------
  output$pk_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    d_long <- d %>%
      select(week, scenario, C_TCZ, C_ADA, C_MTX, C_BARI) %>%
      pivot_longer(cols = c(C_TCZ, C_ADA, C_MTX, C_BARI), names_to = "drug", values_to = "conc") %>%
      filter(conc > 1e-6)
    drug_labels <- c(
      C_TCZ  = "Tocilizumab (nmol/L)",
      C_ADA  = "Adalimumab (nmol/L)",
      C_MTX  = "MTX (mg/L)",
      C_BARI = "Baricitinib (ng/mL)"
    )
    p <- ggplot(d_long, aes(week, conc, color = scenario)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~drug, scales = "free_y",
                 labeller = labeller(drug = drug_labels)) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "Concentration", color = "Scenario",
           title = "Drug Concentration Profiles by Compartment") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  # ---- IL-6R Occupancy -----------------------------------------------------
  output$ro_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, RO_pct, color = scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "gray40") +
      ylim(0, 100) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "IL-6R Occupancy (%)", color = "", title = "TCZ IL-6R Receptor Occupancy") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- MTX Polyglutamates --------------------------------------------------
  output$mtx_poly_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, POLY_MTX, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "MTX Polyglutamates (rel. units)", color = "",
           title = "Intracellular MTX Polyglutamate Pool") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- IL-6 ----------------------------------------------------------------
  output$il6_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, IL6, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "IL-6 (nmol/L)", color = "", title = "Free IL-6 Dynamics") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- pSTAT3 --------------------------------------------------------------
  output$stat3_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, STAT3_P, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "pSTAT3 (rel. units)", color = "", title = "pSTAT3 Signalling") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- FLS -----------------------------------------------------------------
  output$fls_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, FLS_ACT, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "FLS Activation Index", color = "",
           title = "Fibroblast-like Synoviocyte Activation") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- JAK inhibition ------------------------------------------------------
  output$jak_inh_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, JAK_pct, color = scenario)) +
      geom_line(linewidth = 0.8) +
      ylim(0, 100) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "JAK Inhibition (%)", color = "") +
      theme_minimal(base_size = 10)
    ggplotly(p)
  })

  # ---- Sharp/vdH -----------------------------------------------------------
  output$sharp_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, VDH, color = scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "Sharp/vdH Score", color = "", title = "Radiographic Progression") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- Cartilage -----------------------------------------------------------
  output$cart_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, CART, color = scenario)) +
      geom_line(linewidth = 1) +
      ylim(0, 100) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "Cartilage Integrity (%)", color = "", title = "Cartilage Preservation") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- HAQ-DI --------------------------------------------------------------
  output$haq_plot <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    p <- ggplot(d, aes(week, HAQ_DI, color = scenario)) +
      geom_line(linewidth = 1) +
      ylim(0, 3) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "HAQ-DI (0-3)", color = "", title = "Functional Disability Index") +
      theme_minimal()
    ggplotly(p)
  })

  # ---- Response table (wks 12 / 24 / sim_end) ------------------------------
  output$response_table <- renderDT({
    d <- sim_data(); req(nrow(d) > 0)
    checkpoints <- c(12, 24, input$sim_weeks)
    tbl <- map_dfr(checkpoints, function(wk) {
      d %>%
        group_by(scenario) %>%
        filter(abs(week - wk) == min(abs(week - wk))) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(Week = wk)
    }) %>%
      mutate(
        Remission = ifelse(DAS28 < 2.6, "Yes", ""),
        LDA       = ifelse(DAS28 < 3.2, "Yes", "")
      ) %>%
      select(Scenario = scenario, Week, DAS28, CRP, Remission, LDA) %>%
      mutate(DAS28 = round(DAS28, 2), CRP = round(CRP, 1))
    datatable(tbl, options = list(pageLength = 20, dom = "t"), rownames = FALSE)
  })

  # ---- Head-to-head comparison table ---------------------------------------
  output$comparison_table <- renderDT({
    d <- sim_data(); req(nrow(d) > 0)
    wk <- input$sim_weeks
    tbl <- d %>%
      group_by(scenario) %>%
      filter(abs(week - wk) == min(abs(week - wk))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(across(c(DAS28, CRP, INFLAM, VDH, CART, HAQ_DI, RO_pct), ~ round(.x, 2))) %>%
      select(
        Scenario      = scenario,
        DAS28,
        CRP,
        Inflammation  = INFLAM,
        Sharp_vdH     = VDH,
        Cartilage_pct = CART,
        HAQ_DI,
        TCZ_RO_pct    = RO_pct
      )
    datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })

  # ---- Remission bar -------------------------------------------------------
  output$remission_bar <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    wk <- input$sim_weeks
    tbl <- d %>%
      group_by(scenario) %>%
      filter(abs(week - wk) == min(abs(week - wk))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(Remission_pct = as.numeric(DAS28 < 2.6) * 100)
    p <- ggplot(tbl, aes(scenario, Remission_pct, fill = scenario)) +
      geom_col() +
      scale_fill_manual(values = scen_colors) +
      ylim(0, 100) +
      labs(x = "", y = "Remission Rate (%)", title = paste0("Remission (DAS28 <2.6) at Wk ", wk)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
    ggplotly(p)
  })

  # ---- ACR response bar ----------------------------------------------------
  output$acr_bar <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    wk <- input$sim_weeks
    tbl <- d %>%
      group_by(scenario) %>%
      filter(abs(week - wk) == min(abs(week - wk))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        ACR20 = as.numeric(DAS28 < 5.1),
        ACR50 = as.numeric(DAS28 < 3.8),
        ACR70 = as.numeric(DAS28 < 3.2)
      ) %>%
      pivot_longer(c(ACR20, ACR50, ACR70), names_to = "Criterion", values_to = "Response")
    p <- ggplot(tbl, aes(scenario, Response * 100, fill = Criterion)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(ACR20 = "tomato", ACR50 = "orange", ACR70 = "darkred")) +
      ylim(0, 100) +
      labs(x = "", y = "Response (%)", title = paste0("ACR Response at Wk ", wk)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggplotly(p)
  })

  # ---- Multi-biomarker panel -----------------------------------------------
  output$biomarker_panel <- renderPlotly({
    d <- sim_data(); req(nrow(d) > 0)
    d_long <- d %>%
      select(week, scenario, CRP, INFLAM, STAT3_P, FLS_ACT, VDH, CART) %>%
      pivot_longer(-c(week, scenario), names_to = "Biomarker", values_to = "Value")
    bm_labels <- c(
      CRP     = "CRP (mg/L)",
      INFLAM  = "Inflammation Index",
      STAT3_P = "pSTAT3 (rel)",
      FLS_ACT = "FLS Activation",
      VDH     = "Sharp/vdH Score",
      CART    = "Cartilage (%)"
    )
    p <- ggplot(d_long, aes(week, Value, color = scenario)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~Biomarker, scales = "free_y",
                 labeller = labeller(Biomarker = bm_labels), ncol = 3) +
      scale_color_manual(values = scen_colors) +
      labs(x = "Week", y = "", color = "Scenario") +
      theme_minimal(base_size = 10) +
      theme(strip.text = element_text(face = "bold"))
    ggplotly(p)
  })

  # ---- Biomarker reference table -------------------------------------------
  output$biomarker_ref <- renderTable({
    tibble(
      Biomarker    = c("CRP", "ESR", "RF", "Anti-CCP", "IL-6 serum", "TNF-alpha serum", "DAS28-CRP", "HAQ-DI"),
      Normal       = c("<5 mg/L", "<20 mm/h", "<14 IU/mL", "<7 U/mL",  "<7 pg/mL",  "<8 pg/mL",  "<2.6", "0"),
      Active_RA    = c(">15 mg/L", ">30 mm/h", ">40 IU/mL", ">60 U/mL", ">20 pg/mL", ">20 pg/mL", "4-8", "0.5-2.0"),
      Drug_Effect  = c(
        "Decreased by anti-IL-6",
        "Decreased by anti-IL-6",
        "Decreased by rituximab",
        "Decreased by rituximab",
        "Rises after TCZ (receptor blocked)",
        "Decreased by ADA/ETN",
        "Primary endpoint",
        "Functional outcome"
      )
    )
  })

  # ---- Sharp table ---------------------------------------------------------
  output$sharp_table <- renderDT({
    d <- sim_data(); req(nrow(d) > 0)
    wk <- input$sim_weeks
    tbl <- d %>%
      group_by(scenario) %>%
      filter(abs(week - wk) == min(abs(week - wk))) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(Delta_VDH = round(VDH - 10, 2)) %>%
      select(Scenario = scenario, VDH_final = VDH, Delta_VDH) %>%
      arrange(Delta_VDH)
    datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
