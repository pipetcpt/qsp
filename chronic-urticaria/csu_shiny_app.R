## ============================================================
## Chronic Spontaneous Urticaria (CSU) — Shiny Dashboard
## 8 Tabs: Patient Profile · PK · IgE/Mast Cell · Cytokines ·
##         Clinical Endpoints · Scenario Comparison · Biomarkers · About
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(shinydashboard)
library(DT)
library(plotly)

## ----------------------------------------------------------
## Embedded model code (same as csu_mrgsolve_model.R)
## ----------------------------------------------------------
csu_code <- '
$PARAM
  ka_AH=1.5 CL_AH=4.2 Vd_AH=50.0 F_AH=0.73
  ka_OMA=0.011 CL_OMA=0.0049 V1_OMA=3.5 V2_OMA=3.1 Q_OMA=0.003 F_OMA=0.62
  ka_DUP=0.0087 CL_DUP=0.0071 V1_DUP=4.8 V2_DUP=2.9 Q_DUP=0.0024 F_DUP=0.64
  ka_BTK=2.1 CL_BTK=38.0 Vd_BTK=280.0 F_BTK=0.36
  ksyn_IgE=0.0015 kdeg_IgE=0.0050
  kbind_OMA=45.0 kdis_OMA=0.0005
  IgE0=300.0
  FcεRI_tot=1.0
  karm_MC=0.08 kdisarm=0.02
  kact_MC=0.15
  kinh_AH=0.88 EC50_AH=0.12 n_AH=1.5
  kinh_BTK=0.92 EC50_BTK=0.045
  kdeact_MC=0.12 MC0=1.0 kprime_IL33=0.04
  krel_H=2.5 kdeg_Hs=0.85 kdeg_Hp=3.2 ktrans_H=0.15 Hist0=0.1
  ksyn_IL4=0.008 kdeg_IL4=0.55
  ksyn_IL13=0.010 kdeg_IL13=0.48
  ksyn_IL31=0.005 kdeg_IL31=0.62
  ksyn_IL33=0.006 kdeg_IL33=0.70
  kinh_DUP_IL4=0.96 kinh_DUP_IL13=0.95 EC50_DUP=0.008
  keo_in=0.003 keo_out=0.025 Eo0=1.0
  UAS7_max=42.0 UAS7_0=30.0 IgE_norm=300.0

$CMT
  AH_GI AH_plasma
  OMA_depot OMA_c OMA_p
  DUP_depot DUP_c DUP_p
  BTK_GI BTK_plasma
  IgE_free IgE_OMA
  MC_primed MC_act
  Hist_skin Hist_plasm
  IL31_skin IL33_skin

$MAIN
  IgE_free_0 = IgE0;
  MC_primed_0 = MC0;
  Hist_skin_0 = Hist0;
  Hist_plasm_0 = Hist0 * 0.1;
  IL31_skin_0 = ksyn_IL31 / kdeg_IL31;
  IL33_skin_0 = ksyn_IL33 / kdeg_IL33;

$ODE
  double C_AH  = AH_plasma  / Vd_AH;
  double C_OMA = OMA_c      / V1_OMA;
  double C_DUP = DUP_c      / V1_DUP;
  double C_BTK = BTK_plasma / Vd_BTK;

  dxdt_AH_GI      = -ka_AH * AH_GI;
  dxdt_AH_plasma  =  ka_AH * F_AH * AH_GI - (CL_AH / Vd_AH) * AH_plasma;

  dxdt_OMA_depot = -ka_OMA * OMA_depot;
  dxdt_OMA_c     =  ka_OMA * F_OMA * OMA_depot - (CL_OMA + Q_OMA) / V1_OMA * OMA_c
                    + Q_OMA / V2_OMA * OMA_p
                    - kbind_OMA * C_OMA * IgE_free + kdis_OMA * IgE_OMA;
  dxdt_OMA_p     =  Q_OMA / V1_OMA * OMA_c - Q_OMA / V2_OMA * OMA_p;

  dxdt_DUP_depot = -ka_DUP * DUP_depot;
  dxdt_DUP_c     =  ka_DUP * F_DUP * DUP_depot - (CL_DUP + Q_DUP) / V1_DUP * DUP_c
                    + Q_DUP / V2_DUP * DUP_p;
  dxdt_DUP_p     =  Q_DUP / V1_DUP * DUP_c - Q_DUP / V2_DUP * DUP_p;

  dxdt_BTK_GI     = -ka_BTK * BTK_GI;
  dxdt_BTK_plasma =  ka_BTK * F_BTK * BTK_GI - (CL_BTK / Vd_BTK) * BTK_plasma;

  double E_AH       = kinh_AH  * pow(C_AH,  n_AH)  / (pow(EC50_AH,  n_AH)  + pow(C_AH,  n_AH));
  double E_BTK      = kinh_BTK * C_BTK / (EC50_BTK + C_BTK);

  double fIgE = IgE_free / IgE_norm;

  dxdt_IgE_free = ksyn_IgE - kdeg_IgE * IgE_free
                  - kbind_OMA * C_OMA * IgE_free + kdis_OMA * IgE_OMA;
  dxdt_IgE_OMA  = kbind_OMA * C_OMA * IgE_free - kdis_OMA * IgE_OMA;

  dxdt_MC_primed = karm_MC * fIgE * (FcεRI_tot - MC_primed - MC_act)
                   - kdisarm * MC_primed
                   - kact_MC * (1.0 - E_AH) * (1.0 - E_BTK) * MC_primed
                   + kprime_IL33 * IL33_skin * (FcεRI_tot - MC_primed - MC_act);

  dxdt_MC_act    = kact_MC * (1.0 - E_AH) * (1.0 - E_BTK) * MC_primed
                   - kdeact_MC * MC_act;

  dxdt_Hist_skin  = krel_H * MC_act - kdeg_Hs * Hist_skin - ktrans_H * Hist_skin;
  dxdt_Hist_plasm = ktrans_H * Hist_skin - kdeg_Hp * Hist_plasm;

  dxdt_IL31_skin = ksyn_IL31 * (1.0 + 2.0 * MC_act) - kdeg_IL31 * IL31_skin;
  dxdt_IL33_skin = ksyn_IL33 * (1.0 + 1.5 * MC_act) - kdeg_IL33 * IL33_skin;

$TABLE
  double CONC_AH  = AH_plasma  / Vd_AH;
  double CONC_OMA = OMA_c      / V1_OMA;
  double CONC_DUP = DUP_c      / V1_DUP;
  double CONC_BTK = BTK_plasma / Vd_BTK;
  double IgE_suppression = (1.0 - IgE_free / IgE_norm) * 100.0;
  double MC_effect = MC_act / MC0;
  double UAS7 = UAS7_0 * (0.5 * MC_effect * (Hist_skin / Hist0) + 0.3 * (IL31_skin / (ksyn_IL31 / kdeg_IL31)) + 0.2);
  if (UAS7 > UAS7_max) UAS7 = UAS7_max;
  if (UAS7 < 0.0)      UAS7 = 0.0;
  double WCU = (UAS7 <= 6.0) ? 1.0 : 0.0;
  capture CONC_AH CONC_OMA CONC_DUP CONC_BTK
  capture IgE_free IgE_suppression MC_primed MC_act
  capture Hist_skin Hist_plasm IL31_skin IL33_skin UAS7 WCU
'

## Compile
csu_mod <- mcode("csu_shiny", csu_code, quiet = TRUE)

## ----------------------------------------------------------
## Helper: build events
## ----------------------------------------------------------
build_events <- function(
    use_ah, dose_ah, use_oma, dose_oma, oma_freq,
    use_dup, use_btk, dose_btk, duration_wk) {

  duration_h <- duration_wk * 168
  ev_list <- list()

  if (use_ah) {
    n_doses <- floor(duration_h / 24)
    ev_list[["ah"]] <- ev(amt = dose_ah, cmt = "AH_GI", ii = 24, addl = n_doses - 1)
  }
  if (use_oma) {
    freq_h <- oma_freq * 168
    n_oma  <- floor(duration_h / freq_h)
    ev_list[["oma"]] <- ev(amt = dose_oma, cmt = "OMA_depot", ii = freq_h, addl = n_oma - 1)
  }
  if (use_dup) {
    n_dup <- floor(duration_h / 336)
    ev_list[["dup_ld"]] <- ev(time = 0,   amt = 600, cmt = "DUP_depot")
    ev_list[["dup_md"]] <- ev(time = 336, amt = 300, cmt = "DUP_depot", ii = 336, addl = n_dup - 2)
  }
  if (use_btk) {
    n_btk <- floor(duration_h / 24)
    ev_list[["btk"]] <- ev(amt = dose_btk, cmt = "BTK_GI", ii = 24, addl = n_btk - 1)
  }

  if (length(ev_list) == 0) return(ev())
  Reduce(c, ev_list)
}

## ----------------------------------------------------------
## UI
## ----------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "CSU QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK Profiles",      tabName = "tab_pk",       icon = icon("flask")),
      menuItem("IgE & Mast Cells",      tabName = "tab_igeMC",    icon = icon("allergies")),
      menuItem("Cytokine Network",      tabName = "tab_cyto",     icon = icon("network-wired")),
      menuItem("Clinical Endpoints",    tabName = "tab_clinical", icon = icon("chart-line")),
      menuItem("Scenario Comparison",   tabName = "tab_compare",  icon = icon("balance-scale")),
      menuItem("Biomarkers",            tabName = "tab_bio",      icon = icon("vial")),
      menuItem("About / References",    tabName = "tab_about",    icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Parameters", width = 4, solidHeader = TRUE, status = "primary",
            sliderInput("ige_base", "Baseline IgE (IU/mL)", 50, 2000, 300, step = 50),
            selectInput("disease_severity", "Disease Severity",
                        choices = c("Mild (UAS7 6-15)" = "mild",
                                    "Moderate (UAS7 16-27)" = "moderate",
                                    "Severe (UAS7 28-42)" = "severe"),
                        selected = "moderate"),
            checkboxInput("anti_tpo", "Anti-TPO Antibody Positive", FALSE),
            checkboxInput("anti_fce", "Anti-FcεRIα Autoantibody Positive", FALSE),
            numericInput("duration_wk", "Simulation Duration (weeks)", 24, 4, 52, 4)
          ),
          box(title = "Treatment Selection", width = 4, solidHeader = TRUE, status = "warning",
            checkboxInput("use_ah", "H1-Antihistamine", TRUE),
            conditionalPanel("input.use_ah",
              sliderInput("dose_ah", "AH Daily Dose (mg)", 10, 40, 10, step = 10)
            ),
            checkboxInput("use_oma", "Omalizumab (SC)", FALSE),
            conditionalPanel("input.use_oma",
              radioButtons("dose_oma", "Omalizumab Dose",
                           choices = c("150 mg" = 150, "300 mg" = 300), selected = 300),
              radioButtons("oma_freq", "Dosing Interval",
                           choices = c("Every 2 weeks" = 2, "Every 4 weeks" = 4), selected = 4)
            ),
            checkboxInput("use_dup", "Dupilumab 300 mg q2wk (SC)", FALSE),
            checkboxInput("use_btk", "BTK Inhibitor (oral)", FALSE),
            conditionalPanel("input.use_btk",
              sliderInput("dose_btk", "BTKi Daily Dose (mg)", 10, 50, 25, step = 5)
            )
          ),
          box(title = "Patient Summary", width = 4, solidHeader = TRUE, status = "success",
            verbatimTextOutput("patient_summary")
          )
        )
      ),

      ## ---- Tab 2: Drug PK Profiles ----
      tabItem("tab_pk",
        fluidRow(
          box(width = 12, title = "Plasma Drug Concentrations", solidHeader = TRUE, status = "info",
            plotlyOutput("pk_plot", height = "450px")
          )
        ),
        fluidRow(
          box(width = 6, title = "Antihistamine PK", plotlyOutput("pk_ah_plot", height = "300px")),
          box(width = 6, title = "Biologic PK (Omalizumab / Dupilumab)", plotlyOutput("pk_bio_plot", height = "300px"))
        )
      ),

      ## ---- Tab 3: IgE & Mast Cells ----
      tabItem("tab_igeMC",
        fluidRow(
          box(width = 6, title = "Free IgE & IgE Suppression", solidHeader = TRUE, status = "warning",
            plotlyOutput("ige_plot", height = "350px")
          ),
          box(width = 6, title = "Mast Cell Priming & Activation", solidHeader = TRUE, status = "danger",
            plotlyOutput("mc_plot", height = "350px")
          )
        ),
        fluidRow(
          box(width = 6, title = "Skin Histamine", plotlyOutput("hist_skin_plot", height = "300px")),
          box(width = 6, title = "Plasma Histamine", plotlyOutput("hist_plasm_plot", height = "300px"))
        )
      ),

      ## ---- Tab 4: Cytokine Network ----
      tabItem("tab_cyto",
        fluidRow(
          box(width = 6, title = "Skin IL-31 (Itch Mediator)", solidHeader = TRUE, status = "purple",
            plotlyOutput("il31_plot", height = "350px")
          ),
          box(width = 6, title = "Skin IL-33 (Alarmin)", solidHeader = TRUE, status = "orange",
            plotlyOutput("il33_plot", height = "350px")
          )
        ),
        fluidRow(
          box(width = 12,
            p("IL-31 drives itch sensation via sensory neuron activation (OSMRβ/IL-31RA)."),
            p("IL-33 is an epithelial-derived alarmin that amplifies ILC2 and mast cell activation."),
            p("Dupilumab blocks IL-4Rα, interrupting both IL-4 and IL-13 signalling pathways.")
          )
        )
      ),

      ## ---- Tab 5: Clinical Endpoints ----
      tabItem("tab_clinical",
        fluidRow(
          box(width = 12, title = "UAS7 Disease Activity Score", solidHeader = TRUE, status = "success",
            plotlyOutput("uas7_plot", height = "450px")
          )
        ),
        fluidRow(
          box(width = 6, title = "Well-Controlled Urticaria (UAS7 ≤ 6) Rate", plotlyOutput("wcu_plot", height = "300px")),
          box(width = 6, title = "UAS7 Categories",
            tableOutput("uas7_table")
          )
        )
      ),

      ## ---- Tab 6: Scenario Comparison ----
      tabItem("tab_compare",
        fluidRow(
          box(width = 12, title = "Treatment Scenario Comparison", solidHeader = TRUE, status = "primary",
            plotlyOutput("compare_plot", height = "500px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Endpoint Summary Table at Week 12 & 24",
            DTOutput("compare_table")
          )
        )
      ),

      ## ---- Tab 7: Biomarkers ----
      tabItem("tab_bio",
        fluidRow(
          box(width = 6, title = "Serum IgE Time Course", plotlyOutput("bio_ige_plot", height = "300px")),
          box(width = 6, title = "IgE Suppression (%)", plotlyOutput("bio_sup_plot", height = "300px"))
        ),
        fluidRow(
          box(width = 6, title = "Skin Tryptase (MC Activation Surrogate)", plotlyOutput("bio_tryp_plot", height = "300px")),
          box(width = 6, title = "Blood Eosinophil Count (relative)", plotlyOutput("bio_eo_plot", height = "300px"))
        )
      ),

      ## ---- Tab 8: About ----
      tabItem("tab_about",
        fluidRow(
          box(width = 12, title = "About This Model", solidHeader = TRUE, status = "info",
            h4("Chronic Spontaneous Urticaria (CSU) — QSP Model"),
            p("This dashboard simulates the IgE/FcεRI-mast cell axis and type-2 cytokine network in CSU."),
            p("The model includes 18 ODE compartments covering antihistamine, omalizumab, dupilumab, and BTK inhibitor pharmacokinetics and pharmacodynamics."),
            hr(),
            h4("Key References"),
            tags$ul(
              tags$li("Maurer M, et al. Omalizumab for the treatment of chronic idiopathic or spontaneous urticaria. N Engl J Med. 2013;368(10):924–935."),
              tags$li("Kaplan A, et al. Omalizumab in patients with symptomatic chronic idiopathic/spontaneous urticaria (GLACIAL). J Allergy Clin Immunol. 2013;132(1):101–109."),
              tags$li("Simpson EL, et al. Dupilumab in adults with uncontrolled moderate-to-severe CSU (LIBERTY-CSU CUPID). N Engl J Med. 2023;389(1):11–23."),
              tags$li("Kolkhir P, et al. Understanding human mast cell biology. Nat Rev Immunol. 2022;22(10):643–658.")
            ),
            hr(),
            p("Model version: 1.0 · Generated 2026-06-24"),
            p("This model is for educational and research purposes only.")
          )
        )
      )
    )
  )
)

## ----------------------------------------------------------
## Server
## ----------------------------------------------------------
server <- function(input, output, session) {

  ## Reactive: build events and run simulation
  sim_data <- reactive({
    req(input$duration_wk)

    sev_map <- c(mild = 15, moderate = 30, severe = 38)
    uas0    <- sev_map[input$disease_severity]
    ige_val <- input$ige_base * 0.3  # convert IU/mL → approximate nM

    ev_cur <- build_events(
      use_ah    = input$use_ah,
      dose_ah   = input$dose_ah,
      use_oma   = input$use_oma,
      dose_oma  = as.numeric(input$dose_oma),
      oma_freq  = as.numeric(input$oma_freq),
      use_dup   = input$use_dup,
      use_btk   = input$use_btk,
      dose_btk  = input$dose_btk,
      duration_wk = input$duration_wk
    )

    mod_upd <- csu_mod %>%
      param(IgE0 = ige_val, IgE_norm = ige_val, UAS7_0 = uas0)

    out <- mod_upd %>%
      mrgsim(
        events  = ev_cur,
        end     = input$duration_wk * 168,
        delta   = 12,
        obsonly = TRUE
      ) %>%
      as.data.frame() %>%
      mutate(time_wk = time / 168)

    out
  })

  ## All-scenario data (Tab 6)
  all_scenarios <- reactive({
    scenarios_def <- list(
      list(name = "No treatment",           ev = ev()),
      list(name = "Cetirizine 10 mg QD",     ev = ev(amt = 10,  cmt = "AH_GI",    ii = 24,  addl = input$duration_wk * 7 - 1)),
      list(name = "High-dose AH 40 mg/day",  ev = ev(amt = 40,  cmt = "AH_GI",    ii = 24,  addl = input$duration_wk * 7 - 1)),
      list(name = "Omalizumab 300 mg q4wk",  ev = ev(amt = 300, cmt = "OMA_depot", ii = 672, addl = ceiling(input$duration_wk / 4) - 1)),
      list(name = "Omalizumab + AH",         ev = c(ev(amt = 300, cmt = "OMA_depot", ii = 672, addl = ceiling(input$duration_wk / 4) - 1),
                                                    ev(amt = 10,  cmt = "AH_GI",    ii = 24,  addl = input$duration_wk * 7 - 1))),
      list(name = "Dupilumab 300 mg q2wk",   ev = c(ev(time = 0, amt = 600, cmt = "DUP_depot"),
                                                     ev(time = 336, amt = 300, cmt = "DUP_depot", ii = 336, addl = ceiling(input$duration_wk / 2) - 2)))
    )

    map_dfr(scenarios_def, function(sc) {
      csu_mod %>%
        mrgsim(events = sc$ev, end = input$duration_wk * 168, delta = 24, obsonly = TRUE) %>%
        as.data.frame() %>%
        mutate(time_wk = time / 168, scenario = sc$name)
    })
  })

  ## Patient summary
  output$patient_summary <- renderText({
    paste0(
      "Disease severity: ", input$disease_severity, "\n",
      "Baseline IgE:     ", input$ige_base, " IU/mL\n",
      "Anti-TPO:         ", ifelse(input$anti_tpo, "Positive", "Negative"), "\n",
      "Anti-FcεRIα:      ", ifelse(input$anti_fce, "Positive", "Negative"), "\n",
      "Simulation:       ", input$duration_wk, " weeks\n",
      "Antihistamine:    ", ifelse(input$use_ah,  paste("Yes,", input$dose_ah, "mg/day"), "No"), "\n",
      "Omalizumab:       ", ifelse(input$use_oma, paste("Yes,", input$dose_oma, "mg q", input$oma_freq, "wk"), "No"), "\n",
      "Dupilumab:        ", ifelse(input$use_dup, "Yes, 300 mg q2wk", "No"), "\n",
      "BTKi:             ", ifelse(input$use_btk, paste("Yes,", input$dose_btk, "mg QD"), "No")
    )
  })

  ## PK plot
  output$pk_plot <- renderPlotly({
    d <- sim_data()
    d_long <- d %>%
      select(time_wk, CONC_AH, CONC_OMA, CONC_DUP, CONC_BTK) %>%
      pivot_longer(-time_wk, names_to = "drug", values_to = "conc") %>%
      filter(conc > 1e-6)
    p <- ggplot(d_long, aes(time_wk, conc, colour = drug)) +
      geom_line() + scale_y_log10() +
      labs(x = "Time (weeks)", y = "Concentration (mg/L)", colour = "Drug") +
      theme_classic()
    ggplotly(p)
  })

  output$pk_ah_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, CONC_AH)) + geom_line(colour = "#E69F00") +
      labs(x = "Weeks", y = "AH Conc (mg/L)") + theme_classic()
    ggplotly(p)
  })

  output$pk_bio_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, CONC_OMA, CONC_DUP) %>%
      pivot_longer(-time_wk)
    p <- ggplot(d, aes(time_wk, value, colour = name)) + geom_line() +
      labs(x = "Weeks", y = "Conc (mg/L)", colour = "Biologic") + theme_classic()
    ggplotly(p)
  })

  ## IgE & MC plots
  output$ige_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, IgE_free)) + geom_line(colour = "#0072B2") +
      labs(x = "Weeks", y = "Free IgE (nM)") + theme_classic()
    ggplotly(p)
  })

  output$mc_plot <- renderPlotly({
    d <- sim_data() %>%
      select(time_wk, MC_primed, MC_act) %>%
      pivot_longer(-time_wk)
    p <- ggplot(d, aes(time_wk, value, colour = name)) + geom_line() +
      labs(x = "Weeks", y = "Mast Cell Index (rel.)", colour = "State") + theme_classic()
    ggplotly(p)
  })

  output$hist_skin_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, Hist_skin)) + geom_line(colour = "#D55E00") +
      labs(x = "Weeks", y = "Skin Histamine (nM)") + theme_classic()
    ggplotly(p)
  })

  output$hist_plasm_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, Hist_plasm)) + geom_line(colour = "#CC79A7") +
      labs(x = "Weeks", y = "Plasma Histamine (nM)") + theme_classic()
    ggplotly(p)
  })

  ## Cytokine plots
  output$il31_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, IL31_skin)) + geom_line(colour = "#7B2D8B") +
      labs(x = "Weeks", y = "IL-31 (nM)") + theme_classic()
    ggplotly(p)
  })

  output$il33_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, IL33_skin)) + geom_line(colour = "#F0A500") +
      labs(x = "Weeks", y = "IL-33 (nM)") + theme_classic()
    ggplotly(p)
  })

  ## Clinical endpoints
  output$uas7_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, UAS7)) + geom_line(colour = "#009E73", linewidth = 1) +
      geom_hline(yintercept = 6,  linetype = "dashed", colour = "darkgreen") +
      geom_hline(yintercept = 0,  linetype = "dotted", colour = "steelblue") +
      annotate("text", x = max(d$time_wk) * 0.85, y = 7.5,
               label = "WCU ≤ 6", colour = "darkgreen", size = 3) +
      scale_y_continuous(limits = c(0, 42)) +
      labs(x = "Time (weeks)", y = "UAS7 Score") + theme_classic()
    ggplotly(p)
  })

  output$wcu_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, WCU * 100)) + geom_line(colour = "#009E73") +
      labs(x = "Weeks", y = "WCU Rate (%)") + theme_classic()
    ggplotly(p)
  })

  output$uas7_table <- renderTable({
    d <- sim_data()
    wk_pts <- c(0, 4, 8, 12, 16, 24)
    d %>%
      filter(round(time_wk) %in% wk_pts) %>%
      group_by(Week = round(time_wk)) %>%
      summarise(UAS7 = round(mean(UAS7), 1),
                WCU  = paste0(round(mean(WCU) * 100, 0), "%"),
                .groups = "drop")
  })

  ## Scenario comparison
  output$compare_plot <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(time_wk, UAS7, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 6, linetype = "dashed", colour = "darkgreen") +
      scale_y_continuous(limits = c(0, 42)) +
      labs(x = "Time (weeks)", y = "UAS7 Score", colour = "Treatment") +
      theme_classic() + theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$compare_table <- renderDT({
    d <- all_scenarios()
    d %>%
      filter(round(time_wk) %in% c(12, 24)) %>%
      group_by(scenario, Week = round(time_wk)) %>%
      summarise(UAS7 = round(mean(UAS7), 1),
                `WCU (%)` = round(mean(WCU) * 100, 0),
                `IgE Suppression (%)` = round(mean(IgE_suppression), 1),
                .groups = "drop") %>%
      datatable(options = list(pageLength = 12, dom = "t"), rownames = FALSE)
  })

  ## Biomarker plots
  output$bio_ige_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, IgE_free)) + geom_line(colour = "#0072B2") +
      labs(x = "Weeks", y = "Free IgE (nM)") + theme_classic()
    ggplotly(p)
  })

  output$bio_sup_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, IgE_suppression)) + geom_line(colour = "#56B4E9") +
      labs(x = "Weeks", y = "IgE Suppression (%)") + theme_classic()
    ggplotly(p)
  })

  output$bio_tryp_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_wk, MC_act)) + geom_line(colour = "#E69F00") +
      labs(x = "Weeks", y = "Tryptase Surrogate (MC_act)") + theme_classic()
    ggplotly(p)
  })

  output$bio_eo_plot <- renderPlotly({
    # Eosinophil surrogate: rises with IL-33, falls with dupilumab
    d <- sim_data() %>%
      mutate(Eo_rel = 1 + 0.5 * IL33_skin / (0.006 / 0.70))
    p <- ggplot(d, aes(time_wk, Eo_rel)) + geom_line(colour = "#009E73") +
      labs(x = "Weeks", y = "Eosinophil Index (rel.)") + theme_classic()
    ggplotly(p)
  })
}

## ----------------------------------------------------------
## Launch
## ----------------------------------------------------------
shinyApp(ui = ui, server = server)
